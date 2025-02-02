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
# conda activate hiscanner_env
folder=/home/yiz188/my_scratch_space/PTA_4442/
python ${HISCANNER_PATH}/step5_inferCN.py \
    --cell_file ${folder}/samples_4442.txt \
    --batch_size 5 \
    --hetsnp_dir ${folder}/phased_hets \
    --bin_dir ${folder}/bins \
    --seg_dir ${folder}/segs \
    --final_call_dir ${folder}/final_calls \
    --lambda_value 16 --max_wgd 1 \
    --chroms 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y \
    --threads 5 --rerun
    