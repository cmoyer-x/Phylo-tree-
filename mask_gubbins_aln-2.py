#!/usr/bin/env python3
"""
Mask Gubbins Alignment Script
------------------------------
Masks recombinant regions identified by Gubbins from a multiple sequence alignment.

This script reads:
1. A FASTA alignment file (output from Gubbins or harvesttools)
2. A GFF file containing recombinant region predictions from Gubbins

It then masks the recombinant regions by replacing them with 'N' characters.

Usage:
    python3 mask_gubbins_aln.py --aln input.fasta --gff predictions.gff --out masked.fasta

Author: Bacterial Genomics Pipeline
Date: February 2026
"""

import argparse
import sys
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("ERROR: BioPython is not installed.", file=sys.stderr)
    print("Install it with: pip install biopython", file=sys.stderr)
    sys.exit(1)


def parse_gff(gff_file):
    """
    Parse Gubbins GFF file to extract recombinant regions.
    
    Args:
        gff_file (str): Path to Gubbins GFF file
        
    Returns:
        list: List of tuples (start, end) representing recombinant regions (0-based)
    """
    recombinant_regions = []
    
    print(f"Parsing GFF file: {gff_file}")
    
    try:
        with open(gff_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                # Skip header lines
                if line.startswith('#') or line.startswith('##'):
                    continue
                
                # Skip empty lines
                line = line.strip()
                if not line:
                    continue
                
                # Parse GFF fields
                parts = line.split('\t')
                
                if len(parts) < 5:
                    print(f"WARNING: Skipping malformed line {line_num}: {line}", file=sys.stderr)
                    continue
                
                try:
                    # GFF uses 1-based inclusive coordinates
                    # Convert to 0-based for Python indexing
                    start = int(parts[3]) - 1
                    end = int(parts[4])  # Keep end as is for exclusive range
                    
                    # Validate coordinates
                    if start < 0 or end < start:
                        print(f"WARNING: Invalid coordinates at line {line_num}: start={start+1}, end={end}", 
                              file=sys.stderr)
                        continue
                    
                    recombinant_regions.append((start, end))
                    
                except (ValueError, IndexError) as e:
                    print(f"WARNING: Could not parse coordinates at line {line_num}: {e}", file=sys.stderr)
                    continue
    
    except FileNotFoundError:
        print(f"ERROR: GFF file not found: {gff_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to parse GFF file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Sort regions by start position
    recombinant_regions.sort()
    
    print(f"Found {len(recombinant_regions)} recombinant regions")
    
    # Print summary statistics
    if recombinant_regions:
        total_bases = sum(end - start for start, end in recombinant_regions)
        print(f"Total bases in recombinant regions: {total_bases:,}")
        min_length = min(end - start for start, end in recombinant_regions)
        max_length = max(end - start for start, end in recombinant_regions)
        avg_length = total_bases / len(recombinant_regions)
        print(f"Region size: min={min_length}, max={max_length}, avg={avg_length:.1f}")
    
    return recombinant_regions


def mask_alignment(aln_file, regions, out_file, mask_char='N'):
    """
    Mask recombinant regions in a multiple sequence alignment.
    
    Args:
        aln_file (str): Input alignment file in FASTA format
        regions (list): List of (start, end) tuples for regions to mask
        out_file (str): Output file path for masked alignment
        mask_char (str): Character to use for masking (default: 'N')
    """
    print(f"\nReading alignment from: {aln_file}")
    
    try:
        records = list(SeqIO.parse(aln_file, "fasta"))
    except FileNotFoundError:
        print(f"ERROR: Alignment file not found: {aln_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to parse alignment file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if not records:
        print("ERROR: No sequences found in alignment file", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(records)} sequences in alignment")
    
    # Get alignment length
    aln_length = len(records[0].seq)
    print(f"Alignment length: {aln_length:,} bp")
    
    # Verify all sequences have same length
    for i, record in enumerate(records):
        if len(record.seq) != aln_length:
            print(f"ERROR: Sequence {record.id} has different length ({len(record.seq)}) than first sequence ({aln_length})",
                  file=sys.stderr)
            sys.exit(1)
    
    # Mask regions in each sequence
    print(f"\nMasking {len(regions)} regions with '{mask_char}'...")
    
    masked_bases = 0
    
    for record in records:
        # Convert to list for efficient modification
        seq = list(str(record.seq))
        
        for start, end in regions:
            # Ensure we don't exceed sequence length
            mask_start = min(start, len(seq))
            mask_end = min(end, len(seq))
            
            # Mask the region
            for i in range(mask_start, mask_end):
                if seq[i] not in ['-', mask_char]:  # Don't count gaps or already masked
                    masked_bases += 1
                seq[i] = mask_char
        
        # Update sequence
        record.seq = Seq(''.join(seq))
    
    # Write masked alignment
    print(f"\nWriting masked alignment to: {out_file}")
    
    try:
        SeqIO.write(records, out_file, "fasta")
        print(f"Successfully wrote {len(records)} sequences")
    except Exception as e:
        print(f"ERROR: Failed to write output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Print summary
    print("\n" + "="*60)
    print("MASKING SUMMARY")
    print("="*60)
    print(f"Input sequences:        {len(records)}")
    print(f"Alignment length:       {aln_length:,} bp")
    print(f"Recombinant regions:    {len(regions)}")
    print(f"Total bases masked:     {masked_bases:,}")
    print(f"Percentage masked:      {100 * masked_bases / (aln_length * len(records)):.2f}%")
    print(f"Output file:            {out_file}")
    print("="*60)


def main():
    """Main function to parse arguments and run masking."""
    parser = argparse.ArgumentParser(
        description='Mask recombinant regions from Gubbins in a multiple sequence alignment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python3 mask_gubbins_aln.py --aln input.fasta --gff predictions.gff --out masked.fasta
    
    # Use a different masking character
    python3 mask_gubbins_aln.py --aln input.fasta --gff predictions.gff --out masked.fasta --mask X

Notes:
    - Input alignment must be in FASTA format
    - GFF file should be the recombination_predictions.gff from Gubbins
    - Regions are masked with 'N' by default
    - All sequences in alignment must be the same length
        """
    )
    
    parser.add_argument(
        '--aln',
        required=True,
        help='Input alignment FASTA file'
    )
    
    parser.add_argument(
        '--gff',
        required=True,
        help='Gubbins recombination predictions GFF file'
    )
    
    parser.add_argument(
        '--out',
        required=True,
        help='Output masked alignment FASTA file'
    )
    
    parser.add_argument(
        '--mask',
        default='N',
        help='Character to use for masking (default: N)'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.0'
    )
    
    args = parser.parse_args()
    
    # Validate mask character
    if len(args.mask) != 1:
        print("ERROR: Mask character must be a single character", file=sys.stderr)
        sys.exit(1)
    
    # Run masking
    print("="*60)
    print("MASK GUBBINS ALIGNMENT")
    print("="*60)
    
    regions = parse_gff(args.gff)
    mask_alignment(args.aln, regions, args.out, args.mask)
    
    print("\nMasking completed successfully!")


if __name__ == "__main__":
    main()
