import os
import pandas as pd
import subprocess
from pathlib import Path
from typing import List, Dict, Any
from ..logger import logger

def create_config_files(config: Dict[str, Any]) -> None:
    """
    Create configuration files for BICseq2 segmentation.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Configuration dictionary containing paths and parameters
    """
    logger.info("Creating configuration files")
    samples = pd.read_csv(config["metadata_path"], sep="\t")
    outdir = Path(config["outdir"])
    
    # Create necessary directories
    folders = ['cfg', 'bins', 'segs', 'temp', 'segcfg', 'readpos']
    for folder in folders:
        (outdir / folder).mkdir(parents=True, exist_ok=True)
    
    for _, row in samples[["bamID", "bam"]].iterrows():
        bamID, bam = row
        
        # Create directory for bins/{sample}
        (outdir / 'bins' / bamID).mkdir(parents=True, exist_ok=True)
        
        # Create cfg file
        cfg_file = outdir / 'cfg' / f'{bamID}.cfg'
        cfg_data = []
        for chrom in config["chrom_list"]:
            cfg_data.append([
                chrom,
                f"{config['fasta_folder']}/{chrom}.fasta",
                f"{config['mappability_folder_stem']}chr{chrom}.txt",
                f"{outdir}/readpos/{bamID}/{chrom}.readpos.seq",
                f"{outdir}/bins/{bamID}/{chrom}.bin"
            ])
        pd.DataFrame(
            cfg_data, 
            columns=['chrom_name', 'fa_file', 'mappability', 'readPosFile', 'bin_file_normalized']
        ).to_csv(cfg_file, index=None, sep='\t')
        
        # Create seg.cfg file
        segcfg_file = outdir / 'segcfg' / f'{bamID}.seg.cfg'
        segcfg_data = [
            [chrom, f"{outdir}/bins/{bamID}/{chrom}.bin"] 
            for chrom in config["chrom_list"]
        ]
        pd.DataFrame(
            segcfg_data,
            columns=['chromName', 'binFileNorm']
        ).to_csv(segcfg_file, index=None, sep='\t')

def run_segmentation(config: Dict[str, Any]) -> None:
    """
    Run the complete segmentation pipeline.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Configuration dictionary containing paths and parameters
    """
    logger.info("Starting segmentation pipeline")
    outdir = Path(config['outdir'])
    
    try:
        # Validate required tools
        if not os.system('which samtools >/dev/null 2>&1') == 0:
            raise RuntimeError("samtools not found in PATH")
            
        # Create configuration files
        create_config_files(config)
        
        # Process each sample
        samples = pd.read_csv(config["metadata_path"], sep="\t")
        for _, row in samples.iterrows():
            bamID, bam = row['bamID'], row['bam']
            logger.info(f"Processing sample: {bamID}")
            
            # Create readpos files
            try:
                for chrom in config["chrom_list"]:
                    readpos_dir = outdir / 'readpos' / bamID
                    readpos_dir.mkdir(parents=True, exist_ok=True)
                    readpos_file = readpos_dir / f"{chrom}.readpos.seq"
                    
                    if not readpos_file.exists() or config.get('rerun', False):
                        logger.info(f"Creating read position file for {bamID} chromosome {chrom}")
                        cmd = f"samtools view -q 30 -F 1284 {bam} {chrom} | perl -ane 'print $F[3], \"\\n\";' > {readpos_file}"
                        subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"Error creating read position files for {bamID}: {e}")
                raise
            
            # Run BICseq2 normalization
            try:
                cfg_file = outdir / 'cfg' / f'{bamID}.cfg'
                temp_file = outdir / 'temp' / f'{bamID}.temp'
                
                if not (temp_file.parent / f"{temp_file.name}.seq").exists() or config.get('rerun', False):
                    logger.info(f"Running BICseq2 normalization for {bamID}")
                    cmd = f"{config['bicseq_norm']} -b={config['binsize']} --gc_bin -p=0.0001 {cfg_file} {temp_file}"
                    subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"Error in BICseq2 normalization for {bamID}: {e}")
                raise
            
            # Run BICseq2 segmentation
            try:
                segcfg_file = outdir / 'segcfg' / f'{bamID}.seg.cfg'
                seg_dir = outdir / 'segs' / bamID
                seg_dir.mkdir(parents=True, exist_ok=True)
                
                for lambda_val in config["lambda_range"]:
                    out_file = seg_dir / f"lambda{lambda_val}.cnv"
                    if not out_file.exists() or config.get('rerun', False):
                        logger.info(f"Running BICseq2 segmentation for {bamID} with lambda={lambda_val}")
                        cmd = f"{config['bicseq_seg']} --lambda={lambda_val} --bootstrap --detail {segcfg_file} {out_file}"
                        subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"Error in BICseq2 segmentation for {bamID}: {e}")
                raise
            
        logger.info("Segmentation pipeline completed successfully")
        
    except Exception as e:
        logger.error(f"Error in segmentation pipeline: {e}")
        raise

if __name__ == "__main__":
    # Example configuration for testing
    test_config = {
        "outdir": "./output",
        "fasta_folder": "/path/to/reference",
        "mappability_folder_stem": "/path/to/mappability/",
        "metadata_path": "./metadata.txt",
        "bicseq_norm": "/path/to/NBICseq-norm.pl",
        "bicseq_seg": "/path/to/NBICseq-seg.pl",
        "binsize": 500000,
        "chrom_list": [str(i) for i in range(1, 23)] + ['X', 'Y'],
        "lambda_range": [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
    }
    
    # run_segmentation(test_config)  # Uncomment to test