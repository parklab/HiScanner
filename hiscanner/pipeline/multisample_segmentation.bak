import os
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Any, List
import subprocess
from ..logger import logger

class MultisampleSegmentationError(Exception):
    """Custom exception for multisample segmentation errors"""
    pass

def prepare_multisample_input(
    bins_dir: Path,
    output_dir: Path,
    cells: List[str],
    chroms: List[str]
) -> Dict[str, Path]:
    """
    Prepare input files for multisample segmentation by combining normalized bin files.
    
    Parameters
    ----------
    bins_dir : Path
        Directory containing normalized bin files
    output_dir : Path
        Directory for output files
    cells : List[str]
        List of cell IDs
    chroms : List[str]
        List of chromosomes
        
    Returns
    -------
    Dict[str, Path]
        Dictionary mapping chromosomes to their combined input files
    """
    try:
        # Create temporary directory for combined files
        temp_dir = output_dir / 'temp' / 'multisample'
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        combined_files = {}
        
        # Process each chromosome separately
        for chrom in chroms:
            logger.info(f"Combining bin data for chromosome {chrom}")
            
            # Initialize dataframe for this chromosome
            combined_data = []
            
            # First cell to get bin positions
            first_cell = cells[0]
            bin_file = bins_dir / first_cell / f'{chrom}.bin'
            if not bin_file.exists():
                logger.warning(f"Missing bin file for chromosome {chrom} in {first_cell}")
                continue
                
            base_df = pd.read_csv(bin_file, sep='\t')
            combined_data = [
                base_df['start'].values,
                base_df['end'].values
            ]
            combined_df = pd.DataFrame(np.array(combined_data).T, columns=['start', 'end'])
            
            # Add data for each cell
            for cell in cells:
                bin_file = bins_dir / cell / f'{chrom}.bin'
                if not bin_file.exists():
                    logger.warning(f"Missing bin file for chromosome {chrom} in {cell}")
                    continue
                    
                cell_df = pd.read_csv(bin_file, sep='\t')
                combined_df[f'A_{cell}'] = cell_df['obs'].values
                combined_df[f'B_{cell}'] = cell_df['expected'].values
            
            # Save combined data
            output_file = temp_dir / f'combined_{chrom}.txt'
            combined_df.to_csv(output_file, sep='\t', index=False)
            combined_files[chrom] = output_file
            
        return combined_files
        
    except Exception as e:
        raise MultisampleSegmentationError(f"Error preparing multisample input: {e}")

def run_mbicseq(
    input_file: Path,
    output_file: Path,
    mbicseq_path: str,
    lambda_value: float = 1.2,
    remove_common_cnv: bool = False
) -> None:
    """
    Run MBICseq on a single chromosome file.
    
    Parameters
    ----------
    input_file : Path
        Input file path
    output_file : Path
        Output file path
    mbicseq_path : str
        Path to MBICseq executable
    lambda_value : float
        Lambda parameter for segmentation
    remove_common_cnv : bool
        Whether to remove common CNVs
    """
    try:
        cmd = [
            mbicseq_path,
            "-i", str(input_file),
            "-l", str(lambda_value),
            "-o", str(output_file)
        ]
        
        if remove_common_cnv:
            cmd.append("--rm")
            
        logger.debug(f"Running MBICseq: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        
        if result.stderr:
            logger.debug(f"MBICseq stderr: {result.stderr}")
            
    except subprocess.CalledProcessError as e:
        raise MultisampleSegmentationError(f"MBICseq failed: {e.stderr}")
    except Exception as e:
        raise MultisampleSegmentationError(f"Error running MBICseq: {e}")

def combine_segmentation_results(
    seg_files: Dict[str, Path],
    output_file: Path
) -> None:
    """
    Combine segmentation results from multiple chromosomes.
    
    Parameters
    ----------
    seg_files : Dict[str, Path]
        Dictionary mapping chromosomes to their segmentation files
    output_file : Path
        Path for combined output file
    """
    try:
        seg_all = []
        for chrom, seg_file in seg_files.items():
            if not seg_file.exists():
                logger.warning(f"Missing segmentation file for chromosome {chrom}")
                continue
                
            seg = pd.read_csv(seg_file, sep='\t')
            seg['chrom'] = chrom
            seg = seg[['chrom', 'start', 'end']]
            seg_all.append(seg)
            
        if not seg_all:
            raise MultisampleSegmentationError("No valid segmentation results found")
            
        combined_segs = pd.concat(seg_all)
        combined_segs.to_csv(output_file, sep='\t', index=False)
        
    except Exception as e:
        raise MultisampleSegmentationError(f"Error combining segmentation results: {e}")

def run_multisample_segmentation(config: Dict[str, Any]) -> None:
    """
    Run complete multisample segmentation pipeline.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Configuration dictionary
    """
    try:
        # Validate configuration
        required_params = ['outdir', 'mbicseq_path', 'cells', 'chrom_list']
        for param in required_params:
            if param not in config:
                raise MultisampleSegmentationError(f"Missing required parameter: {param}")
                
        outdir = Path(config['outdir'])
        bins_dir = outdir / 'bins'
        segs_dir = outdir / 'segs'
        
        if not bins_dir.exists():
            raise MultisampleSegmentationError(
                f"Bins directory not found: {bins_dir}\n"
                "Please ensure normalization step has completed successfully."
            )
            
        # Prepare input files
        logger.info("Preparing multisample input files")
        combined_files = prepare_multisample_input(
            bins_dir=bins_dir,
            output_dir=outdir,
            cells=config['cells'],
            chroms=config['chrom_list']
        )
        
        # Run segmentation for each lambda value
        for lambda_val in config.get('lambda_range', [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]):
            logger.info(f"Running segmentation with lambda={lambda_val}")
            
            # Process each chromosome
            seg_files = {}
            for chrom, input_file in combined_files.items():
                output_file = outdir / 'temp' / 'multisample' / f'seg_{chrom}_lambda{lambda_val}.txt'
                
                run_mbicseq(
                    input_file=input_file,
                    output_file=output_file,
                    mbicseq_path=config['mbicseq_path'],
                    lambda_value=lambda_val,
                    remove_common_cnv=config.get('remove_common_cnv', False)
                )
                
                seg_files[chrom] = output_file
            
            # Combine results
            combined_output = segs_dir / f'lambda{lambda_val}.cnv'
            combine_segmentation_results(seg_files, combined_output)
            
        logger.info("Multisample segmentation completed successfully")
        
    except Exception as e:
        raise MultisampleSegmentationError(f"Error in multisample segmentation: {e}")