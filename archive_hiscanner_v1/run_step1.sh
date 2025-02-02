#!/bin/bash 
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem 4000
#SBATCH --mail-type=FAIL,END
#SBATCH -n 1

# conda activate scan2
cd /n/scratch/users/y/yiz188/PTA_4442/
scan2 -d scan2_out init
cd scan2_out
scan2 config \
    --verbose \
    --ref /n/data1/hms/dbmi/park/yifan/tools/SCAN2_REF/human_g1k_v37_decoy.fasta \
    --dbsnp /n/data1/hms/dbmi/park/yifan/tools/SCAN2_REF/common_all_20180423.vcf \
    --shapeit-refpanel  /n/data1/hms/dbmi/park/yifan/tools/SCAN2_REF/1000GP_Phase3 \
    --regions-file /n/data1/hms/dbmi/park/yifan/tools/SCAN2_REF/gatk_regions_example.txt \
    --bulk-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/HSB-4442-WGS.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-A3.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-A4.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-A5.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-A6.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-B2.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-B3.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-B4.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-B5.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-B6.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-C1.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-C2.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-C3.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-C4.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-C5.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-C6.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-D1.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-D2.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-D3.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-D4.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-D5.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-D6.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-E1.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-E2.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-E3.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-E4.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-E5.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-E6.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-E7.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-F3.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-F4.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-F5.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-F6.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-F7.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-G1.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-G4.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-G5.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-G6.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-H2.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-H3.bam \
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-H6.bam 

scan2 validate
scan2 run --joblimit 5000 --snakemake-args ' --until shapeit/phased_hets.vcf.gz --latency-wait 120'  \
    --cluster ' sbatch -p park -A park_contrib --mem={resources.mem} -t 120:00:00 -o %logdir/slurm-%A.log'

