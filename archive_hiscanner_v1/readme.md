### HiScanner Pipeline for O2 (Unpackaged)

#### Environment Setup
```bash
conda create -n hiscanner_env python=3.8
conda activate hiscanner_env
pip install hiscanner seaborn hmmlearn
conda install -c conda-forge snakemake py-bgzip tabix graphviz samtools r-base
```
This covers all dependencies for HiScanner.

#### Step 1: SNP Calling and Phasing (via [SCAN2](https://github.com/parklab/SCAN2), most computationally intensive step)
Follow SCAN2's setup guide to create a conda environment and install the necessary packages.
```bash
conda activate scan2
sbatch run_step1.sh
```

#### Step 2: Select Heterozygous SNPs and Compute BAF

For step 2 and onwards, a metadata file is required. This file should contain the following columns:
- bamID (unique identifier for each sample)
- bamPath (path to the bam file)
- singlecell ("Y" or "N")

Here's an example metadata file:
```bash
bam     bamID   singlecell
/n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/HSB-4442-WGS.bam        HSB-4442-WGS    N
/n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-A3.bam PTA_HSB-4442-WGS-A3     Y
/n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-A4.bam PTA_HSB-4442-WGS-A4     Y
/n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-A5.bam PTA_HSB-4442-WGS-A5     Y
```

```bash
conda activate hiscanner_env # keep the hiscanner_env active for the rest of the steps
bash run_step2.sh # or submit as a slurm job
```


#### Step 3: Analyze ADO Patterns for Optimal Bin Size
```bash
bash run_step3.sh
```

#### Step 4: Normalization and Segmentation (second most computationally intensive step)
First, edit step4_segment/config.yaml to specify the bin size and the path to the metadata and mappability files.

```bash
cd step4_segment
sbatch example_usage.sh
```

#### Step 5: CNV Calling
```bash
sbatch run_step5.sh
```

#### Download external references (only hg19/v37 listed here was tested)
N.B.: Files listed in 2-5 are required for SCAN2. Refer to SCAN2's github repository for most recent updates.

1. Human reference genome
```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.fai.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.dict.gz
gunzip *.gz
```
2. dbSNP common variants
```bash
wget https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
gunzip common_all_20180423.vcf.gz
# Create Tribble index (see SCAN2 README)
```

3. SHAPEIT reference panel
```bash
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz
tar xzvf 1000GP_Phase3.tgz
tar xzvf 1000GP_Phase3_chrX.tgz
mv genetic_map_chrX_* 1000GP_Phase3_chrX* 1000GP_Phase3
```

4. GATK regions file (or create your own)
```bash
wget https://raw.githubusercontent.com/parklab/SCAN2/refs/heads/main/resources/analysis_regions/analysis_regions_file_hs37d5.txt
```

5. Mappability files
```bash
wget https://www.math.pku.edu.cn/teachers/xirb/downloads/software/BICseq2/Mappability/hg19CRG.100bp.tar.gz
tar xzvf hg19CRG.100bp.tar.gz
```