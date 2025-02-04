import os
import subprocess
import pandas as pd
from pathlib import Path
from typing import Dict, Any, List, Optional
from concurrent.futures import ThreadPoolExecutor
from ..logger import logger
from .snp_calling import SNPCallingError
import numpy as np
import gzip
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
                        config.get('rerun', False),
                        config  # Pass the full config dictionary
                    )
                )
            
            # Wait for all processes and check for errors
            for future in futures:
                future.result()
                
        logger.info("Phasing step completed successfully")
        
    except Exception as e:
        raise PhasingError(f"Error in phasing step: {e}")

def process_cell_phasing(
    cell: str,
    raw_vcf: Path,
    phased_vcf: Path,
    output_dir: Path,
    rerun: bool = False,
    config: Optional[Dict[str, Any]] = None
) -> None:
    """
    Process phasing for a single cell.
    
    Parameters
    ----------
    cell : str
        Cell identifier
    raw_vcf : Path
        Path to raw VCF file
    phased_vcf : Path
        Path to phased VCF file
    output_dir : Path
        Output directory
    rerun : bool, optional
        Whether to rerun existing results
    config : Dict[str, Any], optional
        Configuration dictionary
    """
    cell_output = output_dir / f'{cell}.hetsnp.txt'
    raw_output = output_dir / f'{cell}.raw.txt'
    
    if cell_output.exists() and not rerun:
        logger.info(f"Skipping existing results for {cell}")
        return
        
    try:
        if not raw_output.exists() or rerun:
            # Extract heterozygous SNPs using bcftools
            logger.info(f"Running bcftools query for {cell}")
            cmd = [
                'bcftools', 'query',
                '-s', cell,
                '-f', '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD{0}\\t%AD{1}]\\n',
                str(raw_vcf)
            ]
            
            with open(raw_output, 'w') as f:
                subprocess.run(cmd, stdout=f, check=True)
        else:
            logger.info(f"Using existing raw file for {cell}")
            
        # Process with phasing information
        raw_data = pd.read_csv(
            raw_output,
            sep='\t',
            header=None,
            names=['CHROM', 'POS', 'REF', 'ALT', 'AD_REF', 'AD_ALT'],
            low_memory=False
        )
        
        # Replace '.' with NaN and drop those rows
        raw_data = raw_data.replace('.', np.nan)
        raw_data = raw_data.dropna()
        
        # Convert columns to appropriate types after cleaning
        raw_data = raw_data.astype({
            'CHROM': str,
            'POS': int,
            'REF': str,
            'ALT': str,
            'AD_REF': int,
            'AD_ALT': int
        })
        
        # Count header lines in phased VCF
        skiprows = 0
        with gzip.open(phased_vcf, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    break
                skiprows += 1
        logger.info(f"Skipping {skiprows} rows in phased VCF (header)")
        
        # Read phasing information using gzip
        logger.info(f"Processing phasing information for {cell}")
        phase_data = pd.read_csv(
            phased_vcf,
            compression='gzip',
            sep='\t',
            skiprows=skiprows,
            usecols=['#CHROM', 'POS', 'phasedgt'],
            low_memory=False,
            dtype={
                '#CHROM': str,
                'POS': int,
                'phasedgt': str
            }
        )
        phase_data.columns = ['CHROM', 'POS', 'GT']
        
        logger.info(f"Merging raw and phased data for {cell}")
        # Merge and process
        merged = pd.merge(
            raw_data,
            phase_data,
            on=['CHROM', 'POS'],
            how='inner'
        )
        
        logger.info(f"Applying phasing for {cell}")
        # Apply phasing
        merged['A'] = merged.apply(
            lambda x: x.AD_ALT if x.GT == '1|0' else x.AD_REF,
            axis=1
        )
        merged['B'] = merged.apply(
            lambda x: x.AD_REF if x.GT == '1|0' else x.AD_ALT,
            axis=1
        )
        
        logger.info(f"Calculating BAF for {cell}")
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
        logger.info(f"Results saved for {cell}")
        
        # Clean up intermediate file if needed
        # Make this operation more defensive
        keep_raw_files = config.get('keep_raw_files', False) if config else False
        if not keep_raw_files and raw_output.exists():
            try:
                raw_output.unlink()
            except Exception as e:
                logger.warning(f"Could not remove raw output file for {cell}: {e}")
        
    except Exception as e:
        logger.error(f"Detailed error processing cell {cell}: {str(e)}")
        raise PhasingError(f"Error processing cell {cell}: {str(e)}")
    
    
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