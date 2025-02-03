from typing import Dict, Any
from ..logger import logger
from .workflow import run_segmentation_workflow, WorkflowError

def run_segmentation(config: Dict[str, Any]) -> None:
    """
    Run the segmentation pipeline using Snakemake workflow.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Configuration dictionary containing paths and parameters
    """
    logger.info("Starting segmentation pipeline")
    
    try:
        # Determine if we should use cluster mode
        use_cluster = config.get('use_cluster', False)
        
        # Run the workflow
        run_segmentation_workflow(
            config=config,
            cluster_mode=use_cluster,
            clean=config.get('clean', False)
        )
        
        logger.info("Segmentation pipeline completed successfully")
        
    except WorkflowError as e:
        logger.error(f"Error in segmentation pipeline: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in segmentation pipeline: {e}")
        raise

def create_config_files(config: Dict[str, Any]) -> None:
    """
    Create configuration files for segmentation.
    
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
        
        # Validate required tools
        if not os.system('which samtools >/dev/null 2>&1') == 0:
            raise RuntimeError("samtools not found in PATH")
            
        # Check if R is installed
        if not os.system('which R >/dev/null 2>&1') == 0:
            raise RuntimeError("R not found in PATH. Please install R using: conda install -c conda-forge r-base")
            
        # Check for required R packages
        r_cmd = 'R -q -e "if (!require(mgcv)) quit(status=1)"'
        if not os.system(r_cmd + ' >/dev/null 2>&1') == 0:
            raise RuntimeError("R package 'mgcv' not found. Please install using: conda install -c conda-forge r-mgcv")


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