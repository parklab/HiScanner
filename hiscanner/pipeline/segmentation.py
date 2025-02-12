import os
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Any, List, Optional
from concurrent.futures import ThreadPoolExecutor
from ..logger import logger

class SegmentationError(Exception):
    """Custom exception for segmentation errors"""
    pass

def validate_segmentation_inputs(config: Dict[str, Any]) -> None:
    """Validate inputs required for segmentation step."""
    required_params = [
        'outdir',
        'metadata_path',
        'lambda_range'
    ]
    
    # Add tool requirements based on mode
    if config.get('use_multisample_segmentation', False):
        required_params.append('mbicseq_path')
    else:
        required_params.append('bicseq_seg')
    
    missing = [param for param in required_params if param not in config]
    if missing:
        raise SegmentationError(
            "Missing required parameters for segmentation:\n" +
            "\n".join(f"- {param}" for param in missing)
        )
    
    # Check if bins exist
    bins_dir = Path(config['outdir']) / 'bins'
    if not bins_dir.exists():
        raise SegmentationError(
            f"Bins directory not found: {bins_dir}\n"
            "Please run normalization first using:\n"
            "hiscanner --config config.yaml run --step normalize"
        )
    
    # Check if any bins were created
    if not any(bins_dir.iterdir()):
        raise SegmentationError(
            f"No bin files found in {bins_dir}\n"
            "Please ensure normalization completed successfully."
        )

def run_single_sample_segmentation(
    cell: str,
    config: Dict[str, Any],
    lambda_value: int
) -> None:
    try:
        outdir = Path(config['outdir'])
        segcfg_file = outdir / 'segcfg' / f'{cell}.seg.cfg'
        segs_dir = outdir / 'segs' / cell
        segs_dir.mkdir(parents=True, exist_ok=True)
        temp_dir = outdir / 'temp' 
        output_file = segs_dir / f'lambda{lambda_value}.cnv'
        log_file = outdir / 'logs' / f'segmentation_{cell}_lambda{lambda_value}.log'
        
        # Skip if output exists and rerun not requested
        if output_file.exists() and not config.get('rerun', False):
            logger.debug(f"Skipping existing segmentation for {cell} at lambda={lambda_value}")
            return
        
        cmd = [
            config['bicseq_seg'],
            f"--lambda={lambda_value}",
            "--bootstrap",
            "--detail",
            f"--tmp={temp_dir}",
            str(segcfg_file),
            str(output_file)
        ]
        
        logger.debug(f"Running command: {' '.join(cmd)}")
        
        with open(log_file, 'w') as log:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True
            )
            
        logger.info(f"Completed segmentation for {cell} at lambda={lambda_value}")
        
    except subprocess.CalledProcessError as e:
        raise SegmentationError(
            f"Segmentation failed for {cell} at lambda={lambda_value}\n"
            f"Exit code: {e.returncode}\n"
            f"Check log file: {log_file}"
        )
    except Exception as e:
        raise SegmentationError(f"Error processing {cell} at lambda={lambda_value}: {str(e)}")

def prepare_multisample_input(
    cells: List[str],
    chrom: str,
    bins_dir: Path,
    output_dir: Path
) -> Path:
    """Prepare combined input file for multisample segmentation."""
    try:
        # Create output directory
        temp_dir = output_dir / 'temp' / 'multisample'
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        output_file = temp_dir / f'combined_{chrom}.txt'
        
        # Read first cell's bins to get positions
        first_cell = cells[0]
        bin_file = bins_dir / first_cell / f'{chrom}.bin'
        if not bin_file.exists():
            raise SegmentationError(f"Missing bin file: {bin_file}")
            
        base_df = pd.read_csv(bin_file, sep='\t')
        combined_df = pd.DataFrame({
            'start': base_df['start'],
            'end': base_df['end']
        })
        
        # Add data for each cell
        for cell in cells:
            bin_file = bins_dir / cell / f'{chrom}.bin'
            if not bin_file.exists():
                logger.warning(f"Missing bin file for {cell}, chromosome {chrom}")
                continue
                
            cell_df = pd.read_csv(bin_file, sep='\t')
            combined_df[f'A_{cell}'] = cell_df['obs']
            combined_df[f'B_{cell}'] = cell_df['expected']
        
        # Save combined data
        combined_df.to_csv(output_file, sep='\t', index=False)
        return output_file
        
    except Exception as e:
        raise SegmentationError(f"Error preparing multisample input for {chrom}: {e}")

def run_multisample_segmentation(
    config: Dict[str, Any],
    chrom: str,
    input_file: Path,
    lambda_value: int
) -> None:
    """Run multisample segmentation for a chromosome."""
    try:
        outdir = Path(config['outdir'])
        temp_dir = outdir / 'temp' / 'multisample'
        segs_dir = outdir / 'segs' / 'multisample'
        segs_dir.mkdir(parents=True, exist_ok=True)
        
        output_file = temp_dir / f'seg_{chrom}_lambda{lambda_value}.txt'
        log_file = outdir / 'logs' / f'multisample_segmentation_{chrom}_lambda{lambda_value}.log'
        
        cmd = [
            config['mbicseq_path'],
            '-i', str(input_file),
            '-l', str(lambda_value),
            '-o', str(output_file)
        ]
        
        # Add option to remove common CNVs if requested
        if config.get('remove_common_cnv', False):
            cmd.append('-rm')
        
        logger.debug(f"Running command: {' '.join(cmd)}")
        
        with open(log_file, 'w') as log:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True
            )
            
        logger.info(f"Completed multisample segmentation for chromosome {chrom} at lambda={lambda_value}")
        
    except subprocess.CalledProcessError as e:
        raise SegmentationError(
            f"Multisample segmentation failed for chromosome {chrom} at lambda={lambda_value}\n"
            f"Exit code: {e.returncode}\n"
            f"Check log file: {log_file}"
        )
    except Exception as e:
        raise SegmentationError(f"Error in multisample segmentation for {chrom}: {str(e)}")
    
    
def combine_multisample_results(
    config: Dict[str, Any],
    lambda_value: int
) -> None:
    """
    Combine multisample segmentation results across chromosomes and save for each cell.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Configuration dictionary containing output directory and other settings
    lambda_value : int
        Lambda value for the segmentation
        
    Raises
    ------
    SegmentationError
        If no valid segmentation results are found or if there's an error processing results
    """
    try:
        outdir = Path(config['outdir'])
        temp_dir = outdir / 'temp' / 'multisample'
        segs_dir = outdir / 'segs'
        
        # Get chromosome list and cell list
        chrom_list = config.get('chrom_list', [str(i) for i in range(1, 23)])
        metadata = pd.read_csv(config['metadata_path'], sep='\t')
        cells = metadata.query('singlecell=="Y"')['bamID'].tolist()
        
        # Collect segments from all chromosomes
        seg_all = []
        for chrom in chrom_list:
            seg_file = temp_dir / f'seg_{chrom}_lambda{lambda_value}.txt'
            if not seg_file.exists():
                logger.warning(f"Missing segmentation file for chromosome {chrom}")
                continue
                
            seg = pd.read_csv(seg_file, sep='\t')
            seg['chrom'] = chrom
            seg_all.append(seg[['chrom', 'start', 'end']])
        
        if not seg_all:
            raise SegmentationError("No valid segmentation results found")
            
        # Combine segments
        combined_segs = pd.concat(seg_all)
        
        # Save a copy for each cell
        for cell in cells:
            # Create cell directory if it doesn't exist
            cell_dir = segs_dir / cell
            cell_dir.mkdir(parents=True, exist_ok=True)
            
            # Save combined segments for this cell
            output_file = cell_dir / f'lambda{lambda_value}.cnv'
            combined_segs.to_csv(output_file, sep='\t', index=False)
            
            logger.debug(f"Saved combined segments for cell {cell} at lambda={lambda_value}")
        
        logger.info(f"Successfully saved combined segments for {len(cells)} cells at lambda={lambda_value}")
        
    except Exception as e:
        raise SegmentationError(f"Error combining multisample results: {e}")
    
def run_segmentation(config: Dict[str, Any]) -> None:
    """Run segmentation pipeline in either single-sample or multisample mode."""
    logger.info("Starting segmentation")
    
    try:
        # Validate inputs
        validate_segmentation_inputs(config)
        
        # Get cell list from metadata
        metadata = pd.read_csv(config['metadata_path'], sep='\t')
        cells = metadata.query('singlecell=="Y"')['bamID'].tolist()
        
        # Set up output directories
        outdir = Path(config['outdir'])
        (outdir / 'segs').mkdir(parents=True, exist_ok=True)
        (outdir / 'logs').mkdir(parents=True, exist_ok=True)
        
        # Determine mode and run appropriate pipeline
        if config.get('use_multisample_segmentation', False):
            logger.info("Running multisample segmentation")
            
            bins_dir = outdir / 'bins'
            chrom_list = config.get('chrom_list', [str(i) for i in range(1, 23)])
            
            # Process each lambda value
            for lambda_value in config['lambda_range']:
                logger.info(f"Processing lambda={lambda_value}")
                
                # Process each chromosome
                for chrom in chrom_list:
                    # Prepare combined input
                    input_file = prepare_multisample_input(
                        cells=cells,
                        chrom=chrom,
                        bins_dir=bins_dir,
                        output_dir=outdir
                    )
                    
                    # Run segmentation
                    run_multisample_segmentation(
                        config=config,
                        chrom=chrom,
                        input_file=input_file,
                        lambda_value=lambda_value
                    )
                
                # Combine results across chromosomes
                combine_multisample_results(
                    config=config,
                    lambda_value=lambda_value
                )
                
        else:
            logger.info("Running single-sample segmentation")
            
            # Use ThreadPoolExecutor for parallel processing
            with ThreadPoolExecutor(max_workers=config.get('threads', 1)) as executor:
                futures = []
                for cell in cells:
                    for lambda_value in config['lambda_range']:
                        futures.append(
                            executor.submit(
                                run_single_sample_segmentation,
                                cell,
                                config,
                                lambda_value
                            )
                        )
                
                # Wait for all processes and check for errors
                for future in futures:
                    future.result()
        
        logger.info("Segmentation completed successfully")
        
    except Exception as e:
        raise SegmentationError(f"Error in segmentation pipeline: {e}")