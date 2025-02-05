import os
from pathlib import Path
from typing import Dict, Any, Optional, Union, List
import yaml
from .logger import logger

class ConfigError(Exception):
    """Custom exception for configuration errors"""
    pass

class Config:
    """Configuration manager for HiScanner pipeline"""
    
    def __init__(self):
        """Initialize configuration with defaults"""
        self.config: Dict[str, Any] = {}
        self._load_defaults()
        
    def _load_defaults(self) -> None:
        """Load default configuration from resources"""
        try:
            default_config = Path(__file__).parent / "resources/default_config.yaml"
            with open(default_config) as f:
                self.config = yaml.safe_load(f)
                
            # Add tool paths
            self.config['bicseq_norm'] = str(self._find_tool_path('NBICseq-norm.pl'))
            self.config['bicseq_seg'] = str(self._find_tool_path('NBICseq-seg.pl'))
            self.config['mbicseq_path'] = str(self._find_tool_path('mbicseq'))
                
        except Exception as e:
            raise ConfigError(f"Failed to load default configuration: {e}")
    def _find_tool_path(self, tool_name: str) -> Path:
        """Find path to bundled tool executable."""
        try:
            # Check if tool exists in resources/tools directory
            if tool_name.endswith('seg.pl'):
                tool_path = Path(__file__).parent / "resources/tools/singlesample_seg" / tool_name
            elif tool_name.endswith('norm.pl'):
                tool_path = Path(__file__).parent / "resources/tools/singlesample_norm" / tool_name
            elif tool_name == 'mbicseq':
                    tool_path = Path(__file__).parent / "resources/tools/" / tool_name
            
            if not tool_path.exists():
                raise ConfigError(f"Tool not found: {tool_path}")
                
            # Make executable if needed
            if not tool_path.suffix == '.pl':
                tool_path.chmod(0o755)
                
            return tool_path
            
        except Exception as e:
            raise ConfigError(f"Error finding tool {tool_name}: {e}")
    
    def update_from_file(self, config_path: Union[str, Path]) -> None:
        """
        Update configuration from user-provided YAML file
        
        Parameters
        ----------
        config_path : str or Path
            Path to configuration file
        
        Raises
        ------
        ConfigError
            If configuration is invalid or missing required fields
        """
        try:
            with open(config_path) as f:
                user_config = yaml.safe_load(f)
            self.config.update(user_config)
            self._validate_config()
        except Exception as e:
            raise ConfigError(f"Failed to load configuration from {config_path}: {e}")

    def _validate_config(self) -> None:
        """
        Validate configuration parameters and paths
        
        Raises
        ------
        ConfigError
            If configuration is invalid
        """
        # Required parameters that must be set
        required_params = {
            'metadata_path': 'Metadata file',
            'outdir': 'Output directory',
            'fasta_folder': 'Reference genome folder',
            'mappability_folder_stem': 'Mappability files prefix',
            'bicseq_norm': 'NBICseq-norm.pl path',
            'bicseq_seg': 'NBICseq-seg.pl path'
        }

        # Add SCAN2 requirement if not in RDR-only mode
        if not self.config.get('rdr_only', False):
            required_params['scan2_output'] = 'SCAN2 results directory'
                
        # Check for null or missing required parameters
        missing = []
        for param, desc in required_params.items():
            if not self.config.get(param):
                missing.append(f"{param} ({desc})")
        
        if missing:
            raise ConfigError(
                "Missing required parameters in configuration:\n" +
                "\n".join(f"- {param}" for param in missing)
            )

        # Validate paths exist (except output directory which will be created)
        paths_to_check = {
            'scan2_output': ['gatk/hc_raw.mmq60.vcf.gz', 'shapeit/phased_hets.vcf.gz'],
            'metadata_path': [],
            'fasta_folder': [],
            'bicseq_norm': [],
            'bicseq_seg': []
        }
        
        for param, subpaths in paths_to_check.items():
            path = Path(self.config[param])
            
            # Check main path
            if param != 'outdir' and not path.exists():
                raise ConfigError(f"Path does not exist: {path}")
            
            # Check required subpaths for scan2_output
            if subpaths:
                for subpath in subpaths:
                    full_path = path / subpath
                    if not full_path.exists():
                        raise ConfigError(
                            f"Required file not found: {full_path}\n"
                            f"Please ensure SCAN2 has been run successfully."
                        )

        # Validate executables
        for exe in ['bicseq_norm', 'bicseq_seg']:
            exe_path = Path(self.config[exe])
            if not (exe_path.exists() and os.access(exe_path, os.X_OK)):
                raise ConfigError(
                    f"{exe} is not executable: {exe_path}\n"
                    f"Please check file permissions and path."
                )

        # Validate chromosome list
        if not isinstance(self.config.get('chrom_list', []), list):
            raise ConfigError("chrom_list must be a list of chromosome identifiers")

        # Validate numeric parameters
        numeric_params = {
            'binsize': ('positive integer', lambda x: isinstance(x, int) and x > 0),
            'max_wgd': ('positive integer', lambda x: isinstance(x, int) and x > 0),
            'batch_size': ('positive integer', lambda x: isinstance(x, int) and x > 0),
            'depth_filter': ('non-negative integer', lambda x: isinstance(x, int) and x >= 0),
            'ado_threshold': ('float between 0 and 1', lambda x: isinstance(x, float) and 0 <= x <= 1),
            'threads': ('positive integer', lambda x: isinstance(x, int) and x > 0)
        }

        for param, (desc, validator) in numeric_params.items():
            value = self.config.get(param)
            if value is not None and not validator(value):
                raise ConfigError(f"{param} must be a {desc}")

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value
        
        Parameters
        ----------
        key : str
            Configuration key
        default : Any
            Default value if key not found
            
        Returns
        -------
        Any
            Configuration value
        """
        return self.config.get(key, default)

    def __getitem__(self, key: str) -> Any:
        """
        Get configuration value using dictionary syntax
        
        Parameters
        ----------
        key : str
            Configuration key
            
        Returns
        -------
        Any
            Configuration value
            
        Raises
        ------
        KeyError
            If key not found in configuration
        """
        return self.config[key]

    def create_output_dirs(self) -> None:
        """
        Create output directories specified in configuration
        
        Raises
        ------
        ConfigError
            If directory creation fails
        """
        try:
            outdir = Path(self.config['outdir'])
            
            # Main output subdirectories
            subdirs = [
                'phased_hets',  # Heterozygous SNP data
                'ado',          # ADO analysis results
                'bins',         # Binned data
                'segs',         # Segmentation results
                'final_calls',  # Final CNV calls
                'logs'          # Log files
            ]
            
            for subdir in subdirs:
                (outdir / subdir).mkdir(parents=True, exist_ok=True)
                
        except Exception as e:
            raise ConfigError(f"Failed to create output directories: {e}")

# Global configuration instance
config = Config()