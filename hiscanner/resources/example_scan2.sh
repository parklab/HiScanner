#!/bin/bash 
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem 4000
#SBATCH --mail-type=FAIL,END
#SBATCH -n 1

# conda activate scan2
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
    --sc-bam /n/data1/hms/dbmi/park/jbrew/scan2/NF_Kow/PTA_4442/bams/PTA_HSB-4442-WGS-H6.bam 

scan2 validate
scan2 run --joblimit 5000 --snakemake-args ' --until shapeit/phased_hets.vcf.gz --latency-wait 120'  \
    --cluster ' sbatch -p park -A park_contrib --mem={resources.mem} -t 120:00:00 -o %logdir/slurm-%A.log'

