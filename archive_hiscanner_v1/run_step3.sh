#!/bin/bash
#SBATCH -n 1
#SBATCH -t 0-12:00                      
#SBATCH -p park,short
#SBATCH -A park 
#SBATCH --mail-type=FAIL                
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

HISCANNER_PATH=/home/yiz188/my_scratch_space/hiscanner_11052024/scripts/

folder=/n/scratch/users/y/yiz188/PTA_4442
metadata=samples_4442.txt

# # if you want to run eval_ado for selected cells:
# python $HISCANNER_PATH/step3_evalADO.py \
#  --metadata_path ${folder}/${metadata} \
#  --hetsnp_dir ${folder}/phased_hets \
#  --ado_dir ${folder}/ado_twocell_B6A4 \
#  --plot --depth_filter 0 --cell "PTA_HSB-4442-WGS-B6,PTA_HSB-4442-WGS-A4"  --rerun

# ## if you want to run eval_ado for all cells:
# python $HISCANNER_PATH/step3_evalADO.py \
#  --metadata_path ${folder}/${metadata} \
#  --hetsnp_dir ${folder}/phased_hets \
#  --ado_dir ${folder}/ado \
#  --plot --depth_filter 0 --cell "all" --rerun


## if you want to run eval_ado for low-cov data, just turn on --aggregate and specify --k 
python $HISCANNER_PATH/step3_evalADO.py \
 --metadata_path ${folder}/${metadata} \
--hetsnp_dir ${folder}/phased_hets \
 --ado_dir ${folder}/ado \
 --plot --depth_filter 0 --cell "all" --aggregate --k 10 --rerun