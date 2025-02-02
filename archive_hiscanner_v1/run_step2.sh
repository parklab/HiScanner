#!/bin/bash
#SBATCH -n 1
#SBATCH -t 0-12:00                      
#SBATCH -p park,short
#SBATCH -A park 
#SBATCH --mail-type=FAIL                
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2


# conda activate scanner_env
module load bcftools
folder=/n/scratch/users/y/yiz188/PTA_4442
metadata=samples_4442.txt


HISCANNER_PATH=/home/yiz188/my_scratch_space/hiscanner_11052024/scripts/
vcf=${folder}/scan2_out/gatk/hc_raw.mmq60.vcf
if [ ! -f $vcf.gz -o ! -f $vcf.gz.tbi ]; then
    echo "Compressing and indexing $vcf"
    bgzip -c $vcf > ${vcf}.gz
    tabix -p vcf ${vcf}.gz
fi
echo "Running step2_getHetSnp.py"
python $HISCANNER_PATH/step2_getHetSnp.py \
    --gatk_vcf ${vcf}.gz \
    --phase_file ${folder}/scan2_out/shapeit/phased_hets.vcf \
    --out_dir ${folder} \
    --metadata $metadata --batch_size 2 \
    --rerun

