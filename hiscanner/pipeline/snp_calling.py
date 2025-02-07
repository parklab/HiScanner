import os
import subprocess
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from ..logger import logger
from ..config import ConfigError

class SNPCallingError(Exception):
    """Custom exception for SNP calling errors"""
    pass

def check_file_readability(path: Path) -> bool:
    """
    Check if a file exists and is readable.
    
    Parameters
    ----------
    path : Path
        Path to file to check
        
    Returns
    -------
    bool
        True if file exists and is readable
    """
    return path.exists() and os.access(path, os.R_OK)

def validate_vcf_files(scan2_dir: Path) -> Tuple[List[str], List[str]]:
    """
    Validate VCF files and their indices.
    
    Parameters
    ----------
    scan2_dir : Path
        SCAN2 output directory
        
    Returns
    -------
    Tuple[List[str], List[str]]
        Lists of missing files and warnings
    """
    missing = []
    warnings = []
    
    # Required VCF files and their indices
    required_files = [
        ('gatk/hc_raw.mmq60.vcf.gz', 'Raw variants VCF'),
        ('gatk/hc_raw.mmq60.vcf.gz.tbi', 'Raw variants VCF index'),
        ('shapeit/phased_hets.vcf.gz', 'Phased hets VCF'),
        ('shapeit/phased_hets.vcf.gz.tbi', 'Phased hets VCF index')
    ]
    
    for file_path, description in required_files:
        full_path = scan2_dir / file_path
        if not check_file_readability(full_path):
            missing.append(f"{description}: {full_path}")
    
    # Additional files that should exist
    recommended_files = [
        ('gatk/hc_raw.mmq60.vcf.gz.md5', 'Raw variants VCF MD5'),
        ('shapeit/phased_hets.vcf.gz.md5', 'Phased hets VCF MD5')
    ]
    
    for file_path, description in recommended_files:
        full_path = scan2_dir / file_path
        if not check_file_readability(full_path):
            warnings.append(f"Missing {description}: {full_path}")
            
    return missing, warnings

def validate_scan2_output(scan2_dir: Path) -> bool:
    """
    Validate SCAN2 output directory and files.
    
    Parameters
    ----------
    scan2_dir : Path
        Path to SCAN2 output directory
        
    Returns
    -------
    bool
        True if validation passes
        
    Raises
    ------
    SNPCallingError
        If validation fails
    """
    if not scan2_dir.exists():
        raise SNPCallingError(
            f"SCAN2 output directory not found: {scan2_dir}\n"
            "Please run SCAN2 first or check the path in your configuration."
        )
    
    # Check required directories
    required_dirs = ['gatk', 'shapeit']
    missing_dirs = [
        d for d in required_dirs 
        if not (scan2_dir / d).is_dir()
    ]
    
    if missing_dirs:
        raise SNPCallingError(
            "Missing required SCAN2 output directories:\n" +
            "\n".join(f"- {d}" for d in missing_dirs)
        )
    
    # Validate VCF files
    missing, warnings = validate_vcf_files(scan2_dir)
    
    if missing:
        raise SNPCallingError(
            "Missing required SCAN2 output files:\n" +
            "\n".join(f"- {f}" for f in missing)
        )
    
    if warnings:
        for warning in warnings:
            logger.warning(warning)
    
    return True

def check_bcftools() -> bool:
    """
    Check if bcftools is available.
    
    Returns
    -------
    bool
        True if bcftools is available
        
    Raises
    ------
    SNPCallingError
        If bcftools is not available
    """
    try:
        result = subprocess.run(
            ['bcftools', '--version'], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
        logger.debug(f"Found bcftools: {result.stdout.split()[1]}")
        return True
    except (subprocess.SubprocessError, FileNotFoundError):
        raise SNPCallingError(
            "bcftools not found in PATH.\n"
            "Please install bcftools and ensure it's available in your PATH."
        )


def index_vcf(vcf_path: Path) -> None:
    """
    Create tabix index for VCF file if needed.
    
    Parameters
    ----------
    vcf_path : Path
        Path to VCF file
        
    Raises
    ------
    SNPCallingError
        If indexing fails
    """
    if not vcf_path.with_suffix('.gz.tbi').exists():
        logger.info(f"Creating index for {vcf_path}")
        try:
            result = subprocess.run(
                ['bcftools', 'index', '--tbi', str(vcf_path)],
                check=True,
                capture_output=True,
                text=True
            )
            logger.debug(f"Index created successfully: {vcf_path}.tbi")
        except subprocess.CalledProcessError as e:
            raise SNPCallingError(
                f"Failed to create index for {vcf_path}: {e.stderr}"
            )

def run_snp_calling(config: Dict[str, Any]) -> Dict[str, Path]:
    """
    Prepare SNP calling results from existing SCAN2 output.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Pipeline configuration
        
    Returns
    -------
    Dict[str, Path]
        Paths to prepared SNP calling files
    
    Raises
    ------
    SNPCallingError
        If required files are missing or invalid
    """
    logger.info("Preparing SNP calling results")
    
    try:
        if 'scan2_output' not in config:
            raise SNPCallingError("scan2_output not specified in configuration")
            
        # Validate and get SCAN2 output paths
        scan2_files = prepare_scan2_results(config)
        
        # Ensure VCF files are indexed
        for vcf_path in scan2_files.values():
            index_vcf(vcf_path)
        
        logger.info("SNP calling results prepared successfully")
        return scan2_files
        
    except Exception as e:
        raise SNPCallingError(f"Error in SNP calling preparation: {str(e)}")

def prepare_scan2_results(config: Dict[str, Any]) -> Dict[str, Path]:
    """
    Prepare and validate SCAN2 results for downstream analysis.
    """
    try:
        if not config or 'scan2_output' not in config:
            raise SNPCallingError("scan2_output not specified in configuration")
            
        scan2_dir = Path(config['scan2_output'])
        logger.debug(f"Using SCAN2 output directory: {scan2_dir}")
        
        # Validate SCAN2 output
        validate_scan2_output(scan2_dir)
        
        # Check bcftools availability
        check_bcftools()
        
        # Return validated paths
        return {
            'raw_vcf': scan2_dir / 'gatk/hc_raw.mmq60.vcf.gz',
            'phased_vcf': scan2_dir / 'shapeit/phased_hets.vcf.gz'
        }
        
    except Exception as e:
        raise SNPCallingError(f"Error preparing SCAN2 results: {str(e)}")
    
    
if __name__ == "__main__":
    # Example usage
    example_config = {
        'scan2_output': './scan2_out',
    }
    
    try:
        results = run_snp_calling(example_config)
        print(f"Prepared files: {results}")
    except SNPCallingError as e:
        print(f"Error: {e}")
        
        
def check_samtools() -> bool:
    """
    Check if samtools is available.
    
    Returns
    -------
    bool
        True if samtools is available
        
    Raises
    ------
    SNPCallingError
        If samtools is not available
    """
    try:
        result = subprocess.run(
            ['samtools', '--version'], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
        logger.debug(f"Found samtools: {result.stdout.split()[1]}")
        return True
    except (subprocess.SubprocessError, FileNotFoundError):
        raise SNPCallingError(
            "samtools not found in PATH.\n"
            "Please install samtools and ensure it's available in your PATH."
        )

def check_r_and_mgcv() -> bool:
    """
    Check if R is available and mgcv package is installed.
    
    Returns
    -------
    bool
        True if R and mgcv are available
        
    Raises
    ------
    SNPCallingError
        If R or mgcv package is not available
    """
    # First check if R is available
    try:
        result = subprocess.run(
            ['R', '--version'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
        logger.debug(f"Found R: {result.stdout.split()[0]}")
    except (subprocess.SubprocessError, FileNotFoundError):
        raise SNPCallingError(
            "R not found in PATH.\n"
            "Please install R and ensure it's available in your PATH."
        )

    # Then check if mgcv package is available
    try:
        result = subprocess.run(
            ['Rscript', '-e', 'library(mgcv)'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
        logger.debug("Found R package: mgcv")
        return True
    except subprocess.CalledProcessError:
        raise SNPCallingError(
            "R package 'mgcv' not found.\n"
            "Please install mgcv package using either:\n"
            "  1. R command: install.packages('mgcv')\n"
            "  2. conda: conda install -c conda-forge r-mgcv"
        )