# HiScanner (HIgh-resolution Single-Cell Allelic copy Number callER)
[![PyPI version](https://badge.fury.io/py/hiscanner.svg)](https://badge.fury.io/py/hiscanner)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

HiScanner is a python package for high-resolution single-cell copy number analysis.


## Table of Contents

- [Installation](#installation)
- [Pipeline Overview](#pipeline-overview)
- [Running the Pipeline](#running-the-pipeline)
    - [Step 1: SNP Calling (Prerequisites)](#step-1-snp-calling-prerequisites)
    - [Steps 2-5: HiScanner Analysis](#steps-2-5-hiscanner-analysis)
- [Output Structure](#output-structure)
- [Required External Files](#required-external-files)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)



### Installation
```bash
# Create new conda environment with all dependencies

conda create -n hiscanner_test python=3.8
conda activate hiscanner_test
pip install hiscanner --no-cache-dir
```

Ensure R is installed in your system with mgcv package:
```bash
# Option 1: Use system R
R -e "install.packages('mgcv')"

# Option 2: Install via conda
conda install -c conda-forge r-base  
conda install -c bioconda r-mgcv>=1.8
```
Check if R is installed correctly:
```bash
Rscript -e "library(mgcv)"
```


Other dependencies that need to be installed via conda:
```bash
conda install -c bioconda snakemake samtools>=1.9 bcftools>=1.9
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

#### SCAN2 Prerequisites

[SCAN2](https://github.com/parklab/SCAN2) needs to be run separately before using HiScanner. 

> **Note**: If you only need SCAN2 output for HiScanner (and are not interested in SNV calls), you can save time by running SCAN2 with:
> ```bash
> scan2 run --joblimit 5000 --snakemake-args ' --until shapeit/phased_hets.vcf.gz --latency-wait 120'
> ```
> This will stop at the phasing step, which is sufficient for HiScanner's requirements.

If you have already run SCAN2, ensure you have:
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

2. Edit config.yaml with your paths and parameters.

3. Prepare metadata file which must contain the following columns:
```
bamID    bam    singlecell
bulk1    /path/to/bulk.bam    N
cell1    /path/to/cell1.bam   Y
cell2    /path/to/cell2.bam   Y
```

4. Validate the configuration:
```bash
hiscanner --config config.yaml validate
```

5. Run the pipeline:
```bash
# Run individual steps
hiscanner --config config.yaml run --step snp      # Quick check SCAN2 results
hiscanner --config config.yaml run --step phase    # Process SCAN2 results
hiscanner --config config.yaml run --step ado      # ADO analysis to identify best bin size
## edit config.yaml with the suggested bin size
hiscanner --config config.yaml run --step segment  # Normalization and segmentation
hiscanner --config config.yaml run --step cnv      # CNV calling
## review final_calls/ directory for CNV calls. Adjust lambda accordingly

# There is also an option to run all steps after SCAN2
hiscanner --config config.yaml run --step all
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


## Cleaning Up
HiScanner creates several temporary directories during analysis. You can clean these up using the clean command:
```bash
hiscanner clean
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

Zhao, Y., Luquette, L. J., Veit, A. D., Wang, X., Xi, R., Viswanadham, V. V., ... & Park, P. J. (2024). High-resolution detection of copy number alterations in single cells with HiScanner. bioRxiv, 2024-04.