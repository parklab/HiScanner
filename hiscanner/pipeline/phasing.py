import os
import subprocess
import pandas as pd
from pathlib import Path
from typing import Dict, Any, List
from concurrent.futures import ThreadPoolExecutor
from ..logger import logger
from .snp_calling import SNPCallingError

class PhasingError(Exception):
    """Custom exception for phasing errors"""
    pass

def validate_inputs(config: Dict[str, Any]) -> Dict[str, Path]:
    """
    Validate input files and paths required for phasing.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Pipeline configuration
        
    Returns
    -------
    Dict[str, Path]
        Dictionary of validated paths
    """
    try:
        # Get paths from config
        output_dir = Path(config['outdir'])
        phased_hets_dir = output_dir / 'phased_hets'
        
        # Create output directory
        phased_hets_dir.mkdir(parents=True, exist_ok=True)
        
        # Validate metadata file
        metadata_path = Path(config['metadata_path'])
        if not metadata_path.exists():
            raise PhasingError(f"Metadata file not found: {metadata_path}")
            
        # Validate SCAN2 outputs through snp_calling module
        from .snp_calling import prepare_scan2_results
        scan2_files = prepare_scan2_results(config)
        
        return {
            'output_dir': output_dir,
            'phased_hets_dir': phased_hets_dir,
            'metadata_path': metadata_path,
            'raw_vcf': scan2_files['raw_vcf'],
            'phased_vcf': scan2_files['phased_vcf']
        }
    except Exception as e:
        raise PhasingError(f"Input validation failed: {e}")

def process_cell_phasing(
    cell: str,
    raw_vcf: Path,
    phased_vcf: Path,
    output_dir: Path,
    rerun: bool = False
) -> None:
    """
    Process phasing for a single cell.
    
    Parameters
    ----------
    cell : str
        Cell identifier
    raw_vcf : Path
        Path to raw VCF from SCAN2
    phased_vcf : Path
        Path to phased VCF from SCAN2
    output_dir : Path
        Output directory for results
    rerun : bool
        Whether to rerun existing results
    """
    cell_output = output_dir / f'{cell}.hetsnp.txt'
    
    if cell_output.exists() and not rerun:
        logger.info(f"Skipping existing results for {cell}")
        return
        
    try:
        # Extract heterozygous SNPs
        cmd = [
            'bcftools', 'query',
            '-s', cell,
            '-f', '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD{0}\\t%AD{1}]\\n',
            str(raw_vcf)
        ]
        
        with open(output_dir / f'{cell}.raw.txt', 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
            
        # Process with phasing information
        raw_data = pd.read_csv(
            output_dir / f'{cell}.raw.txt',
            sep='\t',
            header=None,
            names=['CHROM', 'POS', 'REF', 'ALT', 'AD_REF', 'AD_ALT']
        )
        
        # Read phasing information
        phase_data = pd.read_csv(
            phased_vcf,
            sep='\t',
            comment='#',
            usecols=['CHROM', 'POS', 'GT']
        )
        
        # Merge and process
        merged = pd.merge(
            raw_data,
            phase_data,
            on=['CHROM', 'POS'],
            how='inner'
        )
        
        # Apply phasing
        merged['A'] = merged.apply(
            lambda x: x.AD_ALT if x.GT == '1|0' else x.AD_REF,
            axis=1
        )
        merged['B'] = merged.apply(
            lambda x: x.AD_REF if x.GT == '1|0' else x.AD_ALT,
            axis=1
        )
        
        # Calculate BAF
        merged['TOTAL'] = merged['A'] + merged['B']
        merged = merged[merged.TOTAL > 0]
        merged['BAF'] = merged.apply(
            lambda x: min(x.A, x.B) / x.TOTAL,
            axis=1
        )
        
        # Save results
        result = merged[['CHROM', 'POS', 'A', 'B', 'TOTAL', 'BAF']]
        result.to_csv(cell_output, sep='\t', index=False)
        
        # Clean up intermediate file
        (output_dir / f'{cell}.raw.txt').unlink()
        
    except Exception as e:
        raise PhasingError(f"Error processing cell {cell}: {e}")

def run_phasing(config: Dict[str, Any]) -> None:
    """
    Run the phasing step of the pipeline.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Pipeline configuration
    """
    logger.info("Starting phasing step")
    
    try:
        # Validate inputs
        paths = validate_inputs(config)
        
        # Read metadata
        metadata = pd.read_csv(paths['metadata_path'], sep='\t')
        single_cells = metadata[metadata.singlecell == 'Y']['bamID'].tolist()
        
        if not single_cells:
            raise PhasingError("No single cells found in metadata")
            
        logger.info(f"Processing {len(single_cells)} cells")
        
        # Process cells in parallel
        with ThreadPoolExecutor(max_workers=config.get('threads', 1)) as executor:
            futures = []
            for cell in single_cells:
                futures.append(
                    executor.submit(
                        process_cell_phasing,
                        cell,
                        paths['raw_vcf'],
                        paths['phased_vcf'],
                        paths['phased_hets_dir'],
                        config.get('rerun', False)
                    )
                )
            
            # Wait for all processes and check for errors
            for future in futures:
                future.result()
                
        logger.info("Phasing step completed successfully")
        
    except Exception as e:
        raise PhasingError(f"Error in phasing step: {e}")

if __name__ == "__main__":
    # Example usage
    example_config = {
        'outdir': './output',
        'metadata_path': './metadata.txt',
        'scan2_output': './scan2_out',
        'threads': 4
    }
    
    try:
        run_phasing(example_config)
    except PhasingError as e:
        print(f"Error: {e}")