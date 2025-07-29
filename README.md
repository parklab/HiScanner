# HiScanner (HIgh-resolution Single-Cell Allelic copy Number callER)
[![PyPI version](https://badge.fury.io/py/hiscanner.svg)](https://badge.fury.io/py/hiscanner)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

HiScanner is a python package for high-resolution single-cell copy number analysis. It supports two modes of operation:

1. **Standard Pipeline** (RDR + BAF): Full analysis using both read depth ratios (RDR) and B-allele frequencies (BAF)
2. **RDR-only Pipeline**: Simplified analysis using only read depth ratios

We provide a demo dataset and tutorial to help you get started. After installation, see https://github.com/parklab/hiscanner_demo for instructions. The typical run time for the demo without (`--use-cluster` option) is less than 30 minutes.

## Table of Contents

- [Installation](#installation)
- [Required External Files](#required-external-files)
- [Pipeline Options](#pipeline-options)
- [Running the Pipeline](#running-the-pipeline)
- [Output Structure](#output-structure)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

### Installation
```bash
# Create new conda environment with all dependencies
conda create -n hiscanner_test \
    -c conda-forge -c bioconda -c defaults \
    python=3.8 samtools=1.15.1 bcftools=1.13 r-base "r-mgcv>=1.8"
conda activate hiscanner_test
pip install hiscanner --no-cache-dir
```

HiScanner (version 1.3) has been tested with Linux distributions:
- CentOS Linux release 7.9.2009
- Ubuntu 20.04.6 LTS (GNU/Linux 5.4.0-204-generic x86_64)

The typical installation time is less than 5 minutes.


## Required External Files: 

### 1) Mappability Track
<details>
<summary>hg19/GRCh37 (100mer) </summary>

```bash
wget https://www.math.pku.edu.cn/teachers/xirb/downloads/software/BICseq2/Mappability/hg19CRG.100bp.tar.gz --no-check-certificate
tar -xvzf hg19CRG.100bp.tar.gz
```
</details>

<details>
<summary>hg38/GRCh38 (150mer)</summary>

Download from: https://doi.org/10.6084/m9.figshare.28370357.v1

</details>

<details>
<summary>Other reference genomes: instructions on custom track generation</summary>

For other genomes/read length configurations, follow instructions at [CGAP Annotations](https://cgap-annotations.readthedocs.io/en/latest/bic-seq2_mappability.html) to generate mappability tracks.
</details>

<details>
<summary>IMPORTANT: Guidelines when choosing mappability tracks</summary>

- **Shorter mappability track (e.g., 100bp) with longer reads (e.g., 150bp)**: Valid but conservative (some uniquely mappable regions may be missed)
- **Longer mappability track (e.g., 150bp) with shorter reads (e.g., 100bp)**: Not valid, will cause false positives
</details>

### 2) Reference Genome
We require the reference genome fasta to be split into chromosomes, to allow for parallel processing. You can use the following command to split the reference genome:
```bash
samtools faidx /path/to/reference.fasta
mkdir /path/to/reference/split
awk '{print $1}' /path/to/reference.fasta.fai | xargs -I {} sh -c 'samtools faidx /path/to/reference.fasta {} > /path/to/reference/split/{}.fasta'
```


## Pipeline Options

HiScanner offers two pipeline options:

### 1. Standard Pipeline (RDR + BAF)
Uses both read depth ratios and B-allele frequencies for comprehensive CNV analysis.

Steps:
1. SNP Calling (via SCAN2)
2. Phasing & BAF Computation
3. ADO Analysis
4. Segmentation
5. CNV Calling

Requirements:
- SCAN2 output files
- Bulk and single cell BAM files
- Reference genome
- Mappability tracks

### 2. RDR-only Pipeline
Uses only read depth ratios for simplified CNV analysis.

Steps:
1. Segmentation
2. CNV Calling

Requirements:
- Single cell BAM files
- Reference genome
- Mappability tracks

## Running the Pipeline

### Option 1: Standard Pipeline (RDR + BAF)


#### Step 1: SCAN2 (Prerequisites)

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

#### Steps 2-5: HiScanner Pipeline
<details>
<summary>1. Initialize HiScanner project</summary>

```bash
hiscanner init --output ./my_project
cd my_project
```
</details>

<details>
<summary>2. Edit config.yaml</summary>

Edit with your paths and parameters
</details>

<details>
<summary>3. Prepare metadata file</summary>

Must contain the following columns:
```
bamID    bam    singlecell
bulk1    /path/to/bulk.bam    N
cell1    /path/to/cell1.bam   Y 
cell2    /path/to/cell2.bam   Y
```
</details>

<details>
<summary>4. Validate configuration</summary>

```bash
hiscanner validate
```
</details>

<details>
<summary>5. Run the pipeline</summary>

```bash
hiscanner run --step snp      # Check SCAN2 results
hiscanner run --step phase    # Process SCAN2 results
hiscanner run --step ado      # ADO analysis to identify optimal bin size
hiscanner run --step normalize # Normalize read depth ratios
hiscanner run --step segment  # Segmentation
hiscanner run --step cnv      # CNV calling

# Or run all steps at once:
hiscanner run --step all
```
</details>

For ```normalize``` (the most time-consuming step), we provide an option to run with cluster, e.g., 
```bash
hiscanner --config config.yaml run --step normalize --use-cluster
```

### Option 2: RDR-only Pipeline
<details>
<summary>1. Initialize project</summary>

```bash
hiscanner init --output ./my_project
cd my_project
```
</details>

<details>
<summary>2. Edit config.yaml</summary>

- Set `rdr_only: true`
- Configure paths and parameters
</details>

<details>
<summary>3. Prepare metadata file</summary>

Must contain the following columns (bulk samples are not required):
```
bamID    bam    singlecell
cell1    /path/to/cell1.bam   Y 
cell2    /path/to/cell2.bam   Y
```
</details>

<details>
<summary>4. Validate configuration</summary>

```bash
hiscanner validate
```
</details>

<details>
<summary>5. Run the pipeline</summary>

```bash
hiscanner run --step normalize # Normalize read depth ratios
hiscanner run --step segment  # Segmentation
hiscanner run --step cnv      # CNV calling
```
</details>


## Output Structure

```
hiscanner_output/
├── phased_hets/       # Processed heterozygous SNPs (Standard pipeline only)
├── ado/               # ADO analysis results (Standard pipeline only)
├── bins/             # Binned read depth
├── segs/             # Segmentation results
└── final_calls/      # Final CNV calls
```

## Cleaning Up
HiScanner creates several temporary directories during analysis. You can clean these up using the clean command:
```bash
hiscanner clean
```

## Troubleshooting

Common issues:
1. SCAN2 output formatting issues (hg38 users): Newer SCAN2 versions using Eagle phasing (instead of SHAPEIT) may produce incompatible genotype formatting. Solution: Clean the `GT` field in phased_hets.vcf.gz to keep only the phase information (e.g., "0|1", instead of "0|1:5,6:7,8").
2. Missing SCAN2 results: Ensure scan2_output directory is correctly specified. If the vcf is not zip-compressed, you can use `bgzip` to compress it (in scan2 environment).
```bash
bgzip scan2_out/gatk/hc_raw.mmq60.vcf
tabix -p vcf scan2_out/gatk/hc_raw.mmq60.vcf.gz
```
3. File permissions: Check access to BAM files and reference data
4. Memory issues: Adjust batch_size in config.yaml


## Support
HiScanner is currently under active development. For support or questions, please open an issue on our [GitHub repository](github.com/parklab/hiscanner).


## Citation

If you use HiScanner in your research, please cite:

Zhao, Y., Luquette, L. J., Veit, A. D., Wang, X., Xi, R., Viswanadham, V. V., ... & Park, P. J. (2025). High-resolution detection of copy number alterations in single cells with HiScanner. Nature Communications, 16(1), 5477.

## License
HiScanner is freely available for non-commercial use.

