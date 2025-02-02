# hiscanner/__init__.py
from .__version__ import __version__

# Import main functionality
from .pipeline import (
    run_snp_calling,
    run_phasing,
    run_ado_analysis,
    run_segmentation,
    run_cnv_calling
)

from .utils import (
    prep_input_table,
    postprocess,
    annotate_mean,
    annotate_seg,
    draw_track
)

__all__ = [
    '__version__',
    'run_snp_calling',
    'run_phasing',
    'run_ado_analysis',
    'run_segmentation',
    'run_cnv_calling',
    'prep_input_table',
    'postprocess',
    'annotate_mean',
    'annotate_seg',
    'draw_track'
]