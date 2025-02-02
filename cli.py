import click
import sys
from pathlib import Path
from typing import Optional
from hiscanner.config import Config, load_config
from hiscanner.pipeline import run_pipeline
from hiscanner.logger import setup_logging

# Global config instance
config = Config()

@click.group()
@click.option("--config", type=str, help="Path to config file")
@click.option("--verbose", is_flag=True, help="Enable verbose output")
def cli(config_path, verbose):
    """HiScanner: A pipeline for single-cell CNV analysis"""
    setup_logging(verbose)
    if config_path:
        config.update_from_file(config_path)

@cli.command()
@click.option("--steps", type=str, help="Steps to run (comma-separated)")
def run(steps):
    """Run the HiScanner pipeline"""
    run_pipeline(steps)

@cli.command()
@click.option('--output', type=click.Path(), help='Output directory for the new project')
def init(output: Optional[str]):
    """
    Initialize a new HiScanner project.
    
    Creates a new project directory with required configuration files and
    guides the user through initial setup.
    """
    try:
        output_dir = Path(output) if output else Path.cwd()
        
        # Create project structure
        config_template = Path(__file__).parent.joinpath('resources/default_config.yaml')
        cluster_template = Path(__file__).parent.joinpath('resources/cluster_config.yaml')
        
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / 'config.yaml').write_text(config_template.read_text())
        (output_dir / 'cluster.yaml').write_text(cluster_template.read_text())
        
        # Create metadata template
        metadata_template = """# Metadata file for HiScanner
# Tab-separated file with the following columns:
# bamID: Unique identifier for each sample
# bam: Full path to BAM file
# singlecell: "Y" for single-cell samples, "N" for bulk
bamID	bam	singlecell
bulk1	/path/to/bulk.bam	N
cell1	/path/to/cell1.bam	Y
"""
        (output_dir / 'metadata.txt').write_text(metadata_template)
        
        click.echo(f"✓ Initialized HiScanner project in {output_dir}")
        click.echo("\nNext steps:")
        click.echo("1. Edit config.yaml to set REQUIRED parameters:")
        click.echo("   - scan2_output: Path to SCAN2 results")
        click.echo("   - metadata_path: Path to metadata file")
        click.echo("   - outdir: Output directory")
        click.echo("   - fasta_folder: Reference genome path")
        click.echo("   - mappability_folder_stem: Mappability files path")
        click.echo("   - bicseq_norm and bicseq_seg: Paths to NBICseq tools")
        click.echo("\n2. Update metadata.txt with your sample information")
        click.echo("3. Run 'hiscanner validate' to check your configuration")
        click.echo("4. Run 'hiscanner run --step all' to start analysis")
        
    except Exception as e:
        click.echo(f"Error initializing project: {e}", err=True)
        sys.exit(1)

@cli.command()
@click.option('--config', type=click.Path(exists=True), help='Path to config file to validate')
def validate(config_path: Optional[str]):
    """
    Validate project configuration and required files.
    """
    try:
        if config_path:
            config.update_from_file(config_path)
        elif Path('config.yaml').exists():
            config.update_from_file('config.yaml')
        else:
            click.echo("No configuration found. Please specify --config or run from project directory.", err=True)
            sys.exit(1)
            
        required_params = [
            'scan2_output',
            'metadata_path',
            'outdir',
            'fasta_folder',
            'mappability_folder_stem',
            'bicseq_norm',
            'bicseq_seg'
        ]
        
        missing = []
        for param in required_params:
            if not config.get(param):
                missing.append(param)
                
        if missing:
            click.echo("Missing required parameters in config.yaml:", err=True)
            for param in missing:
                click.echo(f"  - {param}", err=True)
            sys.exit(1)
            
        # Validate paths
        paths_to_check = [
            ('SCAN2 output', Path(config['scan2_output'])),
            ('Metadata file', Path(config['metadata_path'])),
            ('Reference folder', Path(config['fasta_folder'])),
            ('NBICseq-norm', Path(config['bicseq_norm'])),
            ('NBICseq-seg', Path(config['bicseq_seg']))
        ]
        
        for name, path in paths_to_check:
            if not path.exists():
                click.echo(f"Error: {name} not found at {path}", err=True)
                sys.exit(1)
                
        click.echo("✓ Configuration validated successfully")
        
    except Exception as e:
        click.echo(f"Error validating configuration: {e}", err=True)
        sys.exit(1)

def main():
    cli()