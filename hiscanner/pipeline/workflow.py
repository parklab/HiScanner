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
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)

class SnakemakeWorkflow:
    """Manages Snakemake workflow execution for normalization"""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.outdir = Path(config['outdir']).resolve()
        self.workflow_dir = self.outdir / '.workflow'
        self.snakemake_log_dir = self.workflow_dir / 'logs'


    def _find_tool_path(self, tool_name: str) -> Path:
        """Find path to bundled tool executable."""
        try:
            # Determine tool location based on platform
            if tool_name == 'NBICseq-norm.pl':
                tool_path = Path(__file__).parent.parent / 'resources' / 'tools' / 'singlesample_norm' / tool_name
            else:
                raise WorkflowError(f"Unknown tool: {tool_name}")
                
            if not tool_path.exists():
                raise WorkflowError(
                    f"Tool not found: {tool_path}\n"
                    "This might indicate a corrupted or incomplete installation.\n"
                    "Try reinstalling HiScanner: pip install --no-cache-dir --force-reinstall hiscanner"
                )
                    
            # Make executable if needed
            if not tool_path.suffix == '.pl':
                tool_path.chmod(0o755)
                    
            return tool_path
                
        except Exception as e:
            raise WorkflowError(f"Error finding tool {tool_name}: {e}")

    def _generate_config(self) -> None:
        """Generate Snakemake config.yaml from package config"""
        try:
            # Required configuration fields
            required_fields = [
                'outdir',
                'metadata_path',
                'fasta_folder',
                'mappability_folder_stem'
            ]
            
            # Validate required fields with detailed error messages
            missing = []
            for field in required_fields:
                if field not in self.config:
                    missing.append(field)
                elif self.config[field] is None:
                    missing.append(f"{field} (is set to None)")
                elif isinstance(self.config[field], str) and not self.config[field].strip():
                    missing.append(f"{field} (is empty)")
            
            if missing:
                error_msg = [
                    "Missing or invalid configuration fields:",
                    *[f"- {field}" for field in missing],
                    "\nPlease check your config.yaml file and ensure all required fields are set correctly.",
                    "\nRequired fields and their purposes:",
                    "- outdir: Directory where output files will be saved",
                    "- metadata_path: Path to your metadata file with sample information",
                    "- fasta_folder: Directory containing split reference genome files",
                    "- mappability_folder_stem: Path prefix to mappability files"
                ]
                raise WorkflowError("\n".join(error_msg))

            # Get tool paths
            try:
                bicseq_norm = str(self._find_tool_path('NBICseq-norm.pl'))
            except Exception as e:
                raise WorkflowError(f"Error locating bundled tools: {e}")
            
            # Create Snakemake configuration
            snakemake_config = {
                'outdir': str(self.outdir),
                'metadata_path': self.config['metadata_path'],
                'fasta_folder': self.config['fasta_folder'],
                'mappability_folder_stem': self.config['mappability_folder_stem'],
                'bicseq_norm': bicseq_norm,
                'binsize': self.config.get('binsize', 500000),
                'chrom_list': self.config.get('chrom_list', 
                    [str(i) for i in range(1, 23)] + ['X', 'Y'])
            }
            
            # Write configuration to file
            config_path = self.workflow_dir / 'config.yaml'
            with open(config_path, 'w') as f:
                yaml.dump(snakemake_config, f, default_flow_style=False)
                
            logger.debug(f"Generated Snakemake config at {config_path}")
                
        except Exception as e:
            raise WorkflowError(f"Error generating config file: {str(e)}")

    def _find_resource_file(self, filename: str) -> Path:
        """Find resource file using multiple methods"""
        try:
            # First check in current directory
            current_dir_file = Path(filename)
            if current_dir_file.exists():
                logger.debug(f"Found {filename} in current directory")
                return current_dir_file.resolve()

            # Then try pkg_resources
            try:
                resource_path = pkg_resources.resource_filename('hiscanner', f'resources/{filename}')
                if os.path.exists(resource_path):
                    logger.debug(f"Found {filename} in package resources")
                    return Path(resource_path)
            except Exception as e:
                logger.debug(f"pkg_resources lookup failed for {filename}: {e}")

            # Finally try finding relative to the current file
            current_dir = Path(__file__).resolve().parent
            possible_paths = [
                current_dir / '..' / 'resources' / filename,  # From pipeline dir
                current_dir / '..' / '..' / 'resources' / filename,  # From project root
            ]

            for path in possible_paths:
                if path.exists():
                    logger.debug(f"Found {filename} in {path.parent}")
                    return path.resolve()

            raise WorkflowError(
                f"Could not find {filename} in:\n"
                f"- Current directory (./{filename})\n"
                f"- Package resources (hiscanner/resources/{filename})\n"
                f"- Relative paths ({', '.join(str(p) for p in possible_paths)})"
            )
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
            for filename in ['Snakefile']:  # Removed cluster.yaml as it's optional
                try:
                    src_path = self._find_resource_file(filename)
                    dst_path = self.workflow_dir / filename
                    
                    logger.debug(f"Copying {src_path} to {dst_path}")
                    shutil.copy2(src_path, dst_path)
                    
                    if not dst_path.exists():
                        raise WorkflowError(f"Failed to copy {filename} to {dst_path}")
                except Exception as e:
                    raise WorkflowError(f"Error copying {filename}: {str(e)}")

            # Only copy cluster.yaml if it exists in current directory
            cluster_yaml = Path('cluster.yaml')
            if cluster_yaml.exists():
                dst_path = self.workflow_dir / 'cluster.yaml'
                shutil.copy2(cluster_yaml, dst_path)
                logger.debug(f"Copied local cluster.yaml to {dst_path}")

            # Generate workflow config
            self._generate_config()
            
            # Create required directories
            for dirname in ['cfg', 'bins', 'temp', 'segcfg', 'readpos']:
                (self.outdir / dirname).mkdir(exist_ok=True)
                
            logger.info(f"Workflow setup completed in {self.workflow_dir}")
                
        except Exception as e:
            raise WorkflowError(f"Error setting up workflow: {str(e)}")
    
    def run_normalization(self, cluster_mode: bool = False) -> None:
        """Run normalization steps in Snakemake workflow"""
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
            f'--stats={self.snakemake_log_dir}/stats.json',
            '--nolock',
            'all'  # Run all targets (which is just normalization now)
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

        logger.info(f"Running normalization workflow from {self.workflow_dir}")
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

            except subprocess.CalledProcessError as e:
                raise WorkflowError(
                    f"Snakemake workflow failed with exit code {e.returncode}\n"
                    f"stdout: {e.stdout}\n"
                    f"stderr: {e.stderr}"
                )

    def clean(self, pattern: Optional[str] = None) -> None:
        """Clean workflow output files"""
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