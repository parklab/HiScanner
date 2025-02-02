import yaml
from pathlib import Path
from typing import Dict, Any

class Config:
    def __init__(self):
        self.config: Dict[str, Any] = {}
        self._load_defaults()
    
    def _load_defaults(self):
        default_config = Path(__file__).parent / "resources/default_config.yaml"
        with open(default_config) as f:
            self.config = yaml.safe_load(f)
    
    def load_user_config(self, config_path: str):
        with open(config_path) as f:
            user_config = yaml.safe_load(f)
        self.config.update(user_config)
        self._validate_config()
    
    def _validate_config(self):
        # Add validation logic
        pass
