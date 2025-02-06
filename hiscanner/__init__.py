from .__version__ import __version__

# Import all pipeline functions
from .pipeline import (
    run_snp_calling,
    run_phasing,
    run_ado_analysis,
    run_segmentation,
    run_cnv_calling,
    run_normalization
)

__all__ = [
    '__version__',
    'run_snp_calling',
    'run_phasing',
    'run_ado_analysis',
    'run_normalization',
    'run_segmentation',
    'run_cnv_calling'
]