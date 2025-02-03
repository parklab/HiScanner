# HiScanner (HIgh-resolution Single-Cell Allelic copy Number callER)
[![PyPI version](https://badge.fury.io/py/hiscanner.svg)](https://badge.fury.io/py/hiscanner)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

HiScanner is a python package for high-resolution single-cell copy number analysis.


## Table of Contents

- [Prerequisites](#prerequisites)
    - [Environment Setup](#environment-setup)
- [Pipeline Overview](#pipeline-overview)
- [Running the Pipeline](#running-the-pipeline)
    - [Step 1: SNP Calling (Prerequisites)](#step-1-snp-calling-prerequisites)
    - [Steps 2-5: HiScanner Analysis](#steps-2-5-hiscanner-analysis)
- [Output Structure](#output-structure)
- [Required External Files](#required-external-files)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)


## Prerequisites
HiScanner requires [`bcftools`](https://samtools.github.io/bcftools/bcftools.html), which must be included in `PATH`. All other dependencies should be installed automatically with instructions below.




### Installation
```bash
# Create new conda environment with all dependencies

conda create -n hiscanner_test python=3.8
conda activate hiscanner_test

conda install -c conda-forge r-base
conda install -c conda-forge r-mgcv>=1.8
conda install bioconda::snakemake
# conda install -c bioconda samtools>=1.9 bcftools>=1.9 tabix py-bgzip
# conda install -c conda-forge graphviz

# Install HiScanner
pip install .
```


## Pipeline Overview

HiScanner works in a modular fashion with five main steps:

1. SNP Calling (via SCAN2, requires separate environment)
2. Heterozygous SNP Selection & BAF Computation
3. ADO Pattern Analysis
4. Normalization & Segmentation 
5. CNV Calling

## Running the Pipeline

### Step 1: SNP Calling (Prerequisites)

SCAN2 needs to be run separately before using HiScanner. If you have already run SCAN2, ensure you have:
- VCF file with raw variants (`gatk/hc_raw.mmq60.vcf.gz`)
- Phased heterozygous variants (`shapeit/phased_hets.vcf`)

- Additionally, we note that the phased genotype field in `phased_hets.vcf` should be named as `phasedgt`. This is the expected output from the SCAN2 pipeline that we have tested with. If your VCF file has a different field name, please manually rename it to `phasedgt` in the VCF.

The expected location is `scan2_out/` in your project directory.

### Steps 2-5: HiScanner Analysis

1. Initialize project:
```bash
hiscanner init --output ./my_project
cd my_project
```

2. Edit config.yaml with your paths:
```yaml
outdir: "./hiscanner_output"
metadata_path: "./metadata.txt"  # Path to your metadata file
scan2_output: "./scan2_out"      # Path to your SCAN2 results

# External tools and reference files
fasta_folder: "/path/to/reference/split"
mappability_folder_stem: "/path/to/mappability/hg19.CRC.100mer."
bicseq_norm: "/path/to/NBICseq-norm.pl"
bicseq_seg: "/path/to/NBICseq-seg.pl"

# Analysis parameters
binsize: 500000
max_wgd: 1
batch_size: 5
depth_filter: 0
ado_threshold: 0.2
```

3. Prepare metadata file (metadata.txt):
```
bamID    bam    singlecell
bulk1    /path/to/bulk.bam    N
cell1    /path/to/cell1.bam   Y
cell2    /path/to/cell2.bam   Y
```

4. Run the pipeline:
```bash
# Run individual steps
hiscanner run --step phase    # Process SCAN2 results
hiscanner run --step ado      # ADO analysis
hiscanner run --step segment  # Normalization and segmentation
hiscanner run --step cnv      # CNV calling

# Or run all steps after SCAN2
hiscanner run --step all
```

## Output Structure

```
hiscanner_output/
├── phased_hets/       # Processed heterozygous SNPs
├── ado/               # ADO analysis results
├── bins/             # Binned read depth
├── segs/             # Segmentation results
└── final_calls/      # Final CNV calls
```

## Required External Files

1. Reference genome (hg19/v37 recommended)
2. Mappability files
3. SCAN2 output files
4. NBICseq tools (for segmentation)

## Troubleshooting

Common issues:
1. Missing SCAN2 results: Ensure scan2_output directory is correctly specified
2. File permissions: Check access to BAM files and reference data
3. Memory issues: Adjust batch_size in config.yaml

For more detailed information, check the log files in hiscanner_output/logs/


## Support
HiScanner is currently under active development. For support or questions, please open an issue on our [GitHub repository](github.com/parklab/hiscanner).


## Citation

If you use HiScanner in your research, please cite:

@article{zhao2024high,
    title={High-resolution detection of copy number alterations in single cells with HiScanner},
    author={\textbf{Yifan Zhao} and Luquette, Lovelace J and Veit, Alexander D and Wang, Xiaochen and Xi, Ruibin and Viswanadham, Vinayak V and Shao, Diane D and Walsh, Christopher A and Yang, Hong Wei and Johnson, Mark D and Park, Peter J},
    journal={Under revision at Nature Communications},
    year={2024},
}
[Citation information to be added]