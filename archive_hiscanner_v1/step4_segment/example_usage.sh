#!/bin/bash
#SBATCH -n 1
#SBATCH -t 0-24:00                      
#SBATCH -p park 
#SBATCH -A park 
#SBATCH --mail-type=FAIL                
#SBATCH --mem=40000

mkdir -p cluster_logs

snakemake --dry-run

snakemake --jobs 100 --rerun-incomplete -p --latency-wait 1200 \
  --cluster-config cluster.yaml \
  --cluster "sbatch --job-name={cluster.job-name} --account={cluster.account} --partition={cluster.partition} --time={cluster.time} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --mem={cluster.mem} --output={cluster.output} --error={cluster.error}"
