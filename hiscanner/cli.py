import click
import sys
from pathlib import Path
from typing import Optional
from .config import Config, ConfigError
from .logger import setup_logging, get_logger
from .__version__ import __version__
from .pipeline import (
    run_snp_calling,
    run_phasing,
    run_ado_analysis,
    run_normalization,
    run_segmentation,
    run_cnv_calling
)
from .pipeline.snp_calling import SNPCallingError

# Initialize logger
logger = get_logger(__name__)

# Global config instance
config = Config()

def find_config_file() -> Optional[Path]:
    """Look for config file in standard locations."""
    # Check current directory first
    local_config = Path('config.yaml')
    if local_config.exists():
        return local_config
    return None

def validate_config_file(ctx, param, value: Optional[str]) -> Optional[Path]:
    if value is None:
        # If no config provided via --config, look in standard locations
        config_path = find_config_file()
        if config_path:
            return config_path
        return None
    try:
        path = Path(value)
        if not path.exists():
            raise click.BadParameter(f"Configuration file not found: {value}")
        return path
    except Exception as e:
        raise click.BadParameter(f"Invalid configuration file: {e}")

def ensure_config_loaded() -> None:
    """Ensure configuration is loaded before running commands."""
    if not config.config:
        config_path = find_config_file()
        if config_path:
            try:
                config.update_from_file(config_path)
            except ConfigError as e:
                logger.error(f"Configuration error: {e}")
                sys.exit(1)
        else:
            logger.error("No configuration file found. Please provide a config file using --config option or ensure config.yaml exists in the current directory")
            sys.exit(1)

@click.group()
@click.version_option(version=__version__, message="HiScanner version %(version)s")
@click.option('--config', 
              type=click.Path(exists=True), 
              callback=validate_config_file,
              help='Path to configuration file')
@click.option('--debug', 
              is_flag=True, 
              help='Enable debug output')
def cli(config: Optional[Path], debug: bool):
    """HiScanner: Single-cell Copy Number Variation Analysis Pipeline"""
    # Setup logging first thing
    log_level = 'DEBUG' if debug else 'INFO'
    setup_logging(log_level=log_level)
    
    if config:
        try:
            # Use the global config instance
            globals()['config'].update_from_file(config)
        except ConfigError as e:
            logger.error(f"Configuration error: {e}")
            sys.exit(1)

@cli.command()
@click.option('--step', 
              type=click.Choice(['snp', 'phase', 'ado', 'normalize', 'segment', 'cnv', 'all']), 
              default='all',
              help='Pipeline step to run')
@click.option('--use-cluster',
              is_flag=True,
              default=False,
              help='Run pipeline on a cluster using SLURM')
def run(step: str, use_cluster: bool):
    """Run the HiScanner pipeline."""
    try:
        ensure_config_loaded()
        config.config['use_cluster'] = use_cluster
        
        logger.debug(f"Loaded configuration: {config.config}")
        
        steps = {
            'snp': run_snp_calling,
            'phase': run_phasing,
            'ado': run_ado_analysis,
            'normalize': run_normalization,
            'segment': run_segmentation,
            'cnv': run_cnv_calling
        }
        
        if step == 'all':
            logger.info("Running complete pipeline...")
            for step_name in ['snp', 'phase', 'ado', 'normalize', 'segment', 'cnv']:
                if config.config.get('rdr_only', False) and step_name in ['snp', 'phase', 'ado']:
                    logger.info(f"Skipping {step_name} step (RDR-only mode)")
                    continue
                logger.info(f"Running {step_name} step...")
                steps[step_name](config.config)
        else:
            if config.config.get('rdr_only', False) and step in ['snp', 'phase', 'ado']:
                logger.error(f"Cannot run {step} step in RDR-only mode")
                sys.exit(1)
            logger.info(f"Running {step} step...")
            steps[step](config.config)
            
        logger.info("✓ Pipeline completed successfully")
        
    except Exception as e:
        logger.error(f"Error running pipeline: {e}")
        sys.exit(1)

@cli.command()
@click.argument('config_path', 
                type=click.Path(exists=True),
                required=False)
def validate(config_path: Optional[str]):
    """Validate configuration and required files."""
    try:
        if config_path:
            config.update_from_file(config_path)
        else:
            ensure_config_loaded()
            
        # Run validation
        config._validate_config()
        logger.info("✓ Configuration validated successfully")
        
        # Test output directory creation
        config.create_output_dirs()
        logger.info("✓ Output directories can be created")
        
        # Only validate SCAN2 outputs if not in RDR-only mode
        if not config.config.get('rdr_only', False):
            from .pipeline.snp_calling import validate_scan2_output
            scan2_dir = Path(config.config['scan2_output'])
            validate_scan2_output(scan2_dir)
            logger.info("✓ SCAN2 outputs validated successfully")
        
        # Check external tools - always needed regardless of mode
        from .pipeline.snp_calling import check_bcftools, check_samtools, check_r_and_mgcv
        check_bcftools()
        logger.info("✓ bcftools validated successfully")
        check_samtools()
        logger.info("✓ samtools validated successfully")
        check_r_and_mgcv()
        logger.info("✓ R and mgcv package validated successfully")
        
        logger.info("All validation checks passed successfully!")
        
    except ConfigError as e:
        logger.error(f"Configuration error: {e}")
        sys.exit(1)
    except SNPCallingError as e:
        logger.error(f"Tool validation error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Validation error: {e}")
        sys.exit(1)

@cli.command()
@click.option('--all', 'clean_all', is_flag=True, help='Remove all temporary and output files')
def clean(clean_all: bool):
    """Remove temporary files and directories."""
    try:
        ensure_config_loaded()
        
        # Default directories to clean
        temp_dirs = ['cfg', 'segcfg', 'readpos', 'temp']
        
        # Get output directory from config if available
        outdir = Path('.')
        if config.config:  # If config is loaded
            outdir = Path(config.config['outdir'])
        
        # Add output directories if --all flag is used
        if clean_all:
            temp_dirs.extend(['bins', 'segs', 'phased_hets', 'ado', 'final_calls'])
        
        removed = []
        for dirname in temp_dirs:
            dir_path = outdir / dirname
            if dir_path.exists():
                try:
                    import shutil
                    shutil.rmtree(dir_path)
                    removed.append(dirname)
                except Exception as e:
                    logger.warning(f"Could not remove {dirname}: {e}")
        
        if removed:
            logger.info(f"✓ Removed directories: {', '.join(removed)}")
        else:
            logger.info("No temporary directories found to clean")
            
    except Exception as e:
        logger.error(f"Error cleaning directories: {e}")
        sys.exit(1)

@cli.command()
@click.option('--output', 
              type=click.Path(), 
              help='Output directory for the new project')
def init(output: Optional[str]):
    """Initialize a new HiScanner project."""
    try:
        output_dir = Path(output) if output else Path.cwd()
        
        # Create project structure
        config_template = Path(__file__).parent / "resources/default_config.yaml"
        cluster_template = Path(__file__).parent / "resources/cluster_config.yaml"
        
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / 'config.yaml').write_text(config_template.read_text())
        (output_dir / 'cluster.yaml').write_text(cluster_template.read_text())
        
        logger.info(f"✓ Initialized HiScanner project in {output_dir}")
        
    except Exception as e:
        logger.error(f"Error initializing project: {e}")
        sys.exit(1)

def main():
    """Entry point for the CLI."""
    try:
        cli()
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()