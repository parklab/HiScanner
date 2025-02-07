from pathlib import Path
import subprocess
import pandas as pd
from typing import Dict, Any
from ..logger import logger
from .workflow import SnakemakeWorkflow, WorkflowError

class NormalizationError(Exception):
    pass

def validate_normalization_inputs(config: Dict[str, Any]) -> None:
    """
    Validate inputs required for normalization step.
    """
    required_params = [
        'bicseq_norm',
        'binsize',
        'outdir',
        'fasta_folder',
        'mappability_folder_stem',
        'metadata_path'
    ]
    
    missing = [param for param in required_params if param not in config]
    if missing:
        raise NormalizationError(
            "Missing required parameters for normalization:\n" +
            "\n".join(f"- {param}" for param in missing)
        )
    
    # Validate paths exist
    paths_to_check = {
        'fasta_folder': "Reference genome folder",
        'metadata_path': "Metadata file"
    }
    
    for param, desc in paths_to_check.items():
        path = Path(config[param])
        if not path.exists():
            raise NormalizationError(f"{desc} not found: {path}")
            
    # Check if reference genome files exist for specified chromosomes
    fasta_dir = Path(config['fasta_folder'])
    chroms = config.get('chrom_list', [str(i) for i in range(1, 23)])
    missing_fastas = []
    for chrom in chroms:
        fasta = fasta_dir / f"{chrom}.fa"
        fasta_alt = fasta_dir / f"{chrom}.fasta"
        if not fasta.exists() and not fasta_alt.exists():
            missing_fastas.append(chrom)
    
    if missing_fastas:
        raise NormalizationError(
            "Missing reference genome files for chromosomes:\n" +
            "\n".join(f"- {chrom}" for chrom in missing_fastas)
        )

    # Check if mappability files exist
    missing_maps = []
    for chrom in chroms:
        map_file = Path(f"{config['mappability_folder_stem']}chr{chrom}.txt")
        if not map_file.exists():
            missing_maps.append(chrom)
    
    if missing_maps:
        raise NormalizationError(
            "Missing mappability files for chromosomes:\n" +
            "\n".join(f"- {chrom}" for chrom in missing_maps)
        )

def run_normalization(config: Dict[str, Any]) -> None:
    """
    Run the normalization step.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Pipeline configuration
    """
    logger.info("Starting normalization pipeline")
    
    try:
        # Validate inputs
        validate_normalization_inputs(config)
        
        # Set up output directory and subdirectories
        outdir = Path(config['outdir'])
        subdirs = ['readpos', 'bins', 'temp', 'cfg', 'segcfg', 'logs']
        for subdir in subdirs:
            (outdir / subdir).mkdir(parents=True, exist_ok=True)

        # Set tmp_dir within outdir for BICseq-norm
        config['tmp_dir'] = str(outdir / 'temp')
        
        # Set default values for optional parameters
        config.setdefault('read_length', 150)
        config.setdefault('fragment_size', 300)

        # Initialize and run Snakemake workflow
        workflow = SnakemakeWorkflow(config)
        workflow.setup_workflow()
        
        # Run normalization with specified cluster mode
        workflow.run_normalization(
            cluster_mode=config.get('use_cluster', False)
        )
        
        # Clean up temporary files if requested
        if config.get('clean', False):
            workflow.clean()
        
        logger.info("Normalization pipeline completed successfully")
        
    except Exception as e:
        raise NormalizationError(f"Error in normalization pipeline: {e}")

def check_normalization_results(config: Dict[str, Any]) -> bool:
    """
    Check if normalization results exist and are valid.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Pipeline configuration
        
    Returns
    -------
    bool
        True if valid normalization results exist
    
    Raises
    ------
    NormalizationError
        If results are incomplete or invalid
    """
    try:
        outdir = Path(config['outdir'])
        bins_dir = outdir / 'bins'
        
        if not bins_dir.exists():
            return False
            
        # Get cell list from metadata
        metadata = pd.read_csv(config['metadata_path'], sep='\t')
        cells = metadata.query('singlecell=="Y"')['bamID'].tolist()
        chroms = config.get('chrom_list', [str(i) for i in range(1, 23)])
        
        # Check if all expected bin files exist and are non-empty
        missing_bins = []
        empty_bins = []
        
        for cell in cells:
            cell_dir = bins_dir / cell
            if not cell_dir.exists():
                missing_bins.extend([f"{cell}/{chrom}" for chrom in chroms])
                continue
                
            for chrom in chroms:
                bin_file = cell_dir / f"{chrom}.bin"
                if not bin_file.exists():
                    missing_bins.append(f"{cell}/{chrom}")
                elif bin_file.stat().st_size == 0:
                    empty_bins.append(f"{cell}/{chrom}")
        
        if missing_bins or empty_bins:
            error_msg = []
            if missing_bins:
                error_msg.append("Missing bin files:")
                error_msg.extend([f"- {b}" for b in missing_bins])
            if empty_bins:
                error_msg.append("\nEmpty bin files:")
                error_msg.extend([f"- {b}" for b in empty_bins])
            raise NormalizationError("\n".join(error_msg))
        
        return True
        
    except Exception as e:
        raise NormalizationError(f"Error checking normalization results: {e}")