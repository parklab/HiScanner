import os
import shutil
import yaml
from pathlib import Path
import subprocess
from typing import Dict, Any, Optional, Union, List
from ..logger import logger
import contextlib
import pkg_resources

class WorkflowError(Exception):
    """Custom exception for workflow errors"""
    pass

@contextlib.contextmanager
def working_directory(path):
    """Context manager for changing directory with safe return"""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)

class SnakemakeWorkflow:
    """Manages Snakemake workflow execution for HiScanner segmentation"""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.outdir = Path(config['outdir']).resolve()
        self.workflow_dir = self.outdir / '.workflow'
        self.snakemake_log_dir = self.workflow_dir / 'logs'

    def _find_resource_file(self, filename: str) -> Path:
        """Find resource file using multiple methods"""
        try:
            # Try pkg_resources first
            try:
                resource_path = pkg_resources.resource_filename('hiscanner', f'resources/{filename}')
                if os.path.exists(resource_path):
                    return Path(resource_path)
            except Exception as e:
                logger.debug(f"pkg_resources lookup failed for {filename}: {e}")

            # Try finding relative to the current file
            current_dir = Path(__file__).resolve().parent
            possible_paths = [
                current_dir / '..' / 'resources' / filename,  # From pipeline dir
                current_dir / '..' / '..' / 'resources' / filename,  # From project root
            ]

            for path in possible_paths:
                if path.exists():
                    return path.resolve()

            raise WorkflowError(f"Could not find resource file: {filename}")
        except Exception as e:
            raise WorkflowError(f"Error locating resource file {filename}: {str(e)}")

    def setup_workflow(self) -> None:
        """Setup Snakemake workflow directory structure and files"""
        try:
            # Create workflow directory structure
            self.workflow_dir.mkdir(parents=True, exist_ok=True)
            self.snakemake_log_dir.mkdir(parents=True, exist_ok=True)
            (self.workflow_dir / 'cluster_logs').mkdir(parents=True, exist_ok=True)
            
            # Copy resource files
            for filename in ['Snakefile', 'cluster.yaml']:
                try:
                    src_path = self._find_resource_file(filename)
                    dst_path = self.workflow_dir / filename
                    
                    logger.debug(f"Copying {src_path} to {dst_path}")
                    shutil.copy2(src_path, dst_path)
                    
                    if not dst_path.exists():
                        raise WorkflowError(f"Failed to copy {filename} to {dst_path}")
                except Exception as e:
                    raise WorkflowError(f"Error copying {filename}: {str(e)}")

            # Generate workflow config
            self._generate_config()
            
            # Create required directories
            for dirname in ['cfg', 'bins', 'segs', 'temp', 'segcfg', 'readpos']:
                (self.outdir / dirname).mkdir(exist_ok=True)
                
            logger.info(f"Workflow setup completed in {self.workflow_dir}")
                
        except Exception as e:
            raise WorkflowError(f"Error setting up workflow: {str(e)}")

    def run(self, cluster_mode: bool = False) -> None:
        """Run Snakemake workflow"""
        if not self.workflow_dir.exists():
            raise WorkflowError("Workflow directory not found. Did you run setup_workflow()?")

        snakefile_path = self.workflow_dir / 'Snakefile'
        if not snakefile_path.exists():
            raise WorkflowError(f"Snakefile not found at {snakefile_path}")

        cmd = [
            'snakemake',
            '--cores', str(self.config.get('threads', 1)),
            '--rerun-incomplete',
            '--latency-wait', '60',
            '--keep-going',
            '--printshellcmds',
            '--verbose',
            f'--directory={self.outdir}',
            f'--snakefile={snakefile_path}',
            f'--configfile={self.workflow_dir}/config.yaml',
            '--stats', str(self.snakemake_log_dir / 'stats.json'),
            '--nolock'
        ]

        if cluster_mode:
            cluster_cmd = (
                "sbatch "
                "--job-name={cluster.job-name} "
                "--account={cluster.account} "
                "--partition={cluster.partition} "
                "--nodes={cluster.nodes} "
                "--ntasks={cluster.ntasks} "
                "--mem={cluster.mem} "
                "--time={cluster.time} "
                "--output={cluster.output} "
                "--error={cluster.error}"
            )
            cmd.extend([
                '--cluster-config', str(self.workflow_dir / 'cluster.yaml'),
                '--cluster', cluster_cmd,
                '--jobs', '100'
            ])

        if self.config.get('dry_run', False):
            cmd.append('--dry-run')

        logger.info(f"Running Snakemake workflow from {self.workflow_dir}")
        logger.debug(f"Command: {' '.join(cmd)}")

        with working_directory(self.workflow_dir):
            try:
                result = subprocess.run(
                    cmd,
                    check=True,
                    capture_output=True,
                    text=True
                )

                if result.stdout:
                    logger.debug(f"Snakemake stdout:\n{result.stdout}")
                if result.stderr:
                    logger.warning(f"Snakemake stderr:\n{result.stderr}")

                if not self.validate_outputs() and not self.config.get('dry_run', False):
                    raise WorkflowError("Workflow completed but expected outputs are missing")

            except subprocess.CalledProcessError as e:
                raise WorkflowError(
                    f"Snakemake workflow failed with exit code {e.returncode}\n"
                    f"stdout: {e.stdout}\n"
                    f"stderr: {e.stderr}"
                )

    def validate_outputs(self) -> bool:
        """
        Validate that expected output files exist
        
        Returns
        -------
        bool
            True if all expected outputs exist, False otherwise
        """
        # Get cell list from metadata
        metadata = Path(self.config['metadata_path'])
        if not metadata.exists():
            raise WorkflowError(f"Metadata file not found: {metadata}")
            
        try:
            import pandas as pd
            cells_df = pd.read_csv(metadata, sep='\t')
            if 'singlecell' in cells_df.columns:
                cells_df = cells_df.query('singlecell=="Y"')
            cells = cells_df['bamID'].tolist()
        except Exception as e:
            raise WorkflowError(f"Error reading metadata file: {e}")

        # Define expected outputs for each cell
        missing_files = []
        
        # Check bins directory
        for cell in cells:
            bin_dir = self.outdir / 'bins' / cell
            for chrom in self.config.get('chrom_list', range(1, 23)):
                bin_file = bin_dir / f'{chrom}.bin'
                if not bin_file.exists():
                    missing_files.append(bin_file)

        # Check segmentation directory
        for cell in cells:
            seg_dir = self.outdir / 'segs' / cell
            for lambda_val in self.config.get('lambda_range', [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]):
                seg_file = seg_dir / f'lambda{lambda_val}.cnv'
                if not seg_file.exists():
                    missing_files.append(seg_file)

        if missing_files:
            logger.debug("Missing output files:")
            for f in missing_files[:10]:  # Show first 10 missing files
                logger.debug(f"  {f}")
            if len(missing_files) > 10:
                logger.debug(f"  ... and {len(missing_files) - 10} more")
            return False
            
        return True

    def _generate_config(self) -> None:
        """Generate Snakemake config.yaml from package config"""
        try:
            # Required configuration fields
            required_fields = [
                'fasta_folder', 
                'mappability_folder_stem',
                'metadata_path', 
                'bicseq_norm', 
                'bicseq_seg'
            ]
            
            # Validate required fields
            for field in required_fields:
                if field not in self.config:
                    raise WorkflowError(f"Missing required configuration field: {field}")
            
            # Create Snakemake configuration
            snakemake_config = {
                'outdir': str(self.outdir),
                'fasta_folder': self.config['fasta_folder'],
                'mappability_folder_stem': self.config['mappability_folder_stem'],
                'metadata_path': self.config['metadata_path'],
                'bicseq_norm': self.config['bicseq_norm'],
                'bicseq_seg': self.config['bicseq_seg'],
                'binsize': self.config.get('binsize', 500000),
                'chrom_list': self.config.get('chrom_list', 
                    [str(i) for i in range(1, 23)] + ['X', 'Y']),
                'lambda_range': self.config.get('lambda_range', 
                    [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048])
            }
            
            # Write configuration to file
            config_path = self.workflow_dir / 'config.yaml'
            with open(config_path, 'w') as f:
                yaml.dump(snakemake_config, f, default_flow_style=False)
                
            logger.debug(f"Generated Snakemake config at {config_path}")
                
        except Exception as e:
            raise WorkflowError(f"Error generating config file: {str(e)}")

    def clean(self, pattern: Optional[str] = None) -> None:
        """
        Clean workflow output files
        
        Parameters
        ----------
        pattern : str, optional
            Specific pattern of files to clean (e.g., "*.temp")
        """
        try:
            if pattern:
                files = list(self.workflow_dir.glob(pattern))
                for f in files:
                    if f.is_file():
                        f.unlink()
            else:
                # Remove everything except config files
                keep = {'config.yaml', 'cluster.yaml', 'Snakefile'}
                for item in self.workflow_dir.iterdir():
                    if item.name not in keep:
                        if item.is_file():
                            item.unlink()
                        elif item.is_dir():
                            shutil.rmtree(item)
                            
            logger.info(f"Cleaned workflow directory: {self.workflow_dir}")
            
        except Exception as e:
            raise WorkflowError(f"Error cleaning workflow directory: {str(e)}")


def run_segmentation_workflow(
    config: Dict[str, Any],
    cluster_mode: bool = False,
    clean: bool = False
) -> None:
    """
    Run segmentation using Snakemake workflow
    
    Parameters
    ----------
    config : Dict[str, Any]
        Configuration dictionary
    cluster_mode : bool
        Whether to run in cluster mode
    clean : bool
        Whether to clean previous workflow files before running
    """
    workflow = SnakemakeWorkflow(config)
    
    if clean:
        workflow.clean()
        
    workflow.setup_workflow()
    workflow.run(cluster_mode=cluster_mode)