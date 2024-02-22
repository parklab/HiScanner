#!/bin/bash 
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem 16G
#SBATCH --mail-type=FAIL,END

bamdir=/n/data1/hms/dbmi/park/yifan/mening/bams
scansnv \
    --cluster ' sbatch -p park -A park_contrib --mem={resources.mem} -t 120:00:00 -o %logdir/slurm-%A.log' \
    --ref /n/data1/hms/dbmi/park/yifan/data/ref/v37/human_g1k_v37_decoy.fasta \
    --dbsnp /n/data1/hms/dbmi/park/yifan/tools/scansnv/dbsnp_138.b37.vcf \
    --shapeit-panel /n/data1/hms/dbmi/park/yifan/tools/scansnv/1000GP_Phase3 \
    --scripts /home/yiz188/.conda/envs/scansnvtest/lib/scansnv/ \
    --regions-file /n/data1/hms/dbmi/park/yifan/tools/scansnv/gatk_regions_example.txt \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --joblimit 1000 --resume --snakemake-args ' --until shapeit/phased_hsnps.vcf'  \
    --output-dir A_001 \
    --bulk-sample BOS-15Blood_S2_L001 \
    --sc-sample    BOS-15H10_S90_L001              \
    --bam       BOS-15Blood_S2_L001              ${bamdir}/BOS-15Blood_S2_L001.bam                \
    --bam       BOS-15H10_S90_L001               ${bamdir}/BOS-15H10_S90_L001.bam                