# Bacterial Genomics Analysis Pipeline For a Phylogentic Tree  

*IF genomes are too dissimilar the pipeline will kick them out and will not include them*

A comprehensive pipeline for bacterial comparative genomics, including core genome alignment, recombination detection, and phylogenetic reconstruction.

## Overview

This pipeline performs the following steps:

1. **Core Genome Alignment** - Uses Parsnp to align all bacterial genomes against a reference
2. **Format Conversion** - Converts alignment to FASTA format using harvesttools
3. **Recombination Detection** - Identifies recombinant regions using Gubbins
4. **Recombination Masking** - Masks recombinant regions from the alignment
5. **Phylogenetic Reconstruction** - Builds a maximum likelihood tree using IQ-TREE

## Requirements

### Software Dependencies

- **Conda/Miniconda** - For environment management
- **Parsnp** - Core genome alignment (part of Harvest suite)
- **harvesttools** - Format conversion utilities
- **Gubbins** - Recombination detection
- **IQ-TREE** - Phylogenetic tree construction
- **Python 3** - For masking script

### Conda Environment Setup

```bash
# Create environment
conda create -n analyze_mab -c bioconda parsnp harvesttools gubbins iqtree python=3.9

# Activate environment
conda activate analyze_mab
```

## File Structure

Before running the pipeline, ensure your directory structure looks like this:

```
project_directory/
├── ATCC19977.fasta              # Reference genome
├── AllStrains/                   # Directory containing strain genomes
│   ├── strain1.fasta
│   ├── strain2.fasta
│   └── ...
├── mask_gubbins_aln.py          # Masking script (must be present)
└── bacterial_genomics_pipeline.sh
```

## Usage

### Basic Usage

```bash
# Make the script executable
chmod +x bacterial_genomics_pipeline.sh

# Run the pipeline
./bacterial_genomics_pipeline.sh
```

### Configuration

Edit these variables in the script to customize:

```bash
REFERENCE="ATCC19977.fasta"      # Your reference genome
STRAINS_DIR="AllStrains"          # Directory with strain genomes
THREADS=8                         # Threads for Parsnp and Gubbins
IQTREE_THREADS=4                  # Threads for IQ-TREE
CONDA_ENV="analyze_mab"           # Conda environment name
```

## Pipeline Steps Explained

### 1. Parsnp (Core Genome Alignment)

- Aligns all genomes against the reference
- Identifies conserved regions across all strains
- Output: `parsnp_output/parsnp.xmfa`

**Key parameters:**
- `-r`: Reference genome
- `-d`: Directory with query genomes
- `-p`: Number of threads
- `--skip-phylogeny`: Skip initial tree building (we'll use IQ-TREE later)

### 2. Harvesttools (Format Conversion)

- Converts Parsnp's XMFA format to FASTA
- Makes alignment compatible with Gubbins
- Output: `harvesttools/gubbins_input.fasta`

### 3. Gubbins (Recombination Detection)

- Identifies regions affected by recombination
- Uses phylogenetic methods to detect horizontal gene transfer
- Output: `gubbins/gubbins_input.recombination_predictions.gff`

**Why mask recombination?**
Recombination can mislead phylogenetic inference by creating false signals of shared ancestry.

### 4. Masking Script (mask_gubbins_aln.py)

- Masks recombinant regions identified by Gubbins
- Replaces recombinant bases with 'N' or removes them
- Output: `gubbins/masked_fasta_file.fasta`

**Note:** You need to provide this script. If you don't have it, here's a basic version:

```python
#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_gff(gff_file):
    """Parse Gubbins GFF file to extract recombinant regions."""
    recombinant_regions = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])
                recombinant_regions.append((start, end))
    return recombinant_regions

def mask_alignment(aln_file, regions, out_file):
    """Mask recombinant regions in alignment."""
    records = list(SeqIO.parse(aln_file, "fasta"))
    
    for record in records:
        seq = list(str(record.seq))
        for start, end in regions:
            for i in range(start, end):
                if i < len(seq):
                    seq[i] = 'N'
        record.seq = Seq(''.join(seq))
    
    SeqIO.write(records, out_file, "fasta")
    print(f"Masked {len(regions)} recombinant regions")
    print(f"Output written to {out_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Mask recombinant regions from alignment')
    parser.add_argument('--aln', required=True, help='Input alignment FASTA')
    parser.add_argument('--gff', required=True, help='Gubbins GFF file')
    parser.add_argument('--out', required=True, help='Output masked FASTA')
    
    args = parser.parse_args()
    
    regions = parse_gff(args.gff)
    mask_alignment(args.aln, regions, args.out)
```

### 5. IQ-TREE (Phylogenetic Reconstruction)

- Builds maximum likelihood phylogenetic tree
- Uses GTR+G model (General Time Reversible + Gamma distribution)
- Performs model selection and bootstrap analysis
- Output: `iqtree/masked_fasta_file.fasta.treefile`

## Output Files

### Primary Outputs

| File | Description |
|------|-------------|
| `parsnp_output/parsnp.xmfa` | Core genome alignment in XMFA format |
| `harvesttools/gubbins_input.fasta` | Alignment in FASTA format |
| `gubbins/gubbins_input.recombination_predictions.gff` | Recombinant regions |
| `gubbins/masked_fasta_file.fasta` | Alignment with recombination masked |
| `iqtree/masked_fasta_file.fasta.treefile` | Final phylogenetic tree (Newick format) |

### Additional IQ-TREE Outputs

- `.iqtree` - Detailed analysis report
- `.log` - Run log
- `.bionj` - Initial tree
- `.mldist` - ML distances

## Visualization

### Tree Visualization

You can visualize the final tree using:

- **FigTree** - Desktop application
- **iTOL** (Interactive Tree of Life) - Web-based
- **ggtree** (R package) - For publication-quality figures


## Troubleshooting

### Common Issues

**Issue:** Parsnp fails with "No genomes found"
- **Solution:** Check that `.fasta` files are in `AllStrains/` directory
- Ensure files are properly formatted FASTA files

**Issue:** Gubbins runs very slowly
- **Solution:** Reduce thread count or consider using a subset of genomes first
- Gubbins is computationally intensive for large datasets

**Issue:** IQ-TREE reports "Alignment too short"
- **Solution:** Check that core genome alignment is sufficiently long
- May need to adjust Parsnp parameters for more lenient alignment

**Issue:** "conda: command not found"
- **Solution:** Initialize conda in your shell
  ```bash
  # For bash
  source ~/miniconda3/etc/profile.d/conda.sh
  # Or add to ~/.bashrc
  ```

### Log Files

The pipeline creates detailed logs with timestamps. Check for:
- `[ERROR]` messages for failures
- `[WARN]` messages for potential issues
- `[INFO]` messages for progress tracking

## Performance Optimization

### For Large Datasets (>50 genomes)

1. **Increase threads:** Adjust `THREADS` and `IQTREE_THREADS`
2. **Use faster models:** Consider simpler substitution models
3. **Parallel execution:** Run Parsnp on subset batches
4. **Memory:** Ensure adequate RAM (8GB+ recommended)

### Expected Runtime

Approximate times for 20 genomes (~5Mb each):
- Parsnp: 10-30 minutes
- Gubbins: 30-60 minutes  
- IQ-TREE: 10-20 minutes

## Citation

If you use this pipeline, please cite:

- **Parsnp:** Treangen TJ, et al. "The Harvest suite for rapid core-genome alignment and visualization of thousands of intraspecific microbial genomes." Genome Biology (2014).
- **Gubbins:** Croucher NJ, et al. "Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins." Nucleic Acids Research (2015).
- **IQ-TREE:** Nguyen LT, et al. "IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies." Molecular Biology and Evolution (2015).

## License

This pipeline script is provided as-is for academic and research purposes.

## Contact

For issues with:
- **Pipeline script:** Check GitHub issues or contact maintainer
- **Individual tools:** Consult tool-specific documentation

---

**Version:** 1.0  
**Last Updated:** February 2026
