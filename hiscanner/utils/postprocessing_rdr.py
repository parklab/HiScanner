import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from ..logger import logger

# Import shared utility functions for general processing
from .postprocessing import (
    get_var,
    get_bic_gauss,
    merge_adjacent_bins
)

# Additional libraries for numeric computations
from scipy import stats
from scipy.signal import find_peaks
from scipy.special import logsumexp

def annotate_mean_rdr(cell_data: pd.DataFrame, 
            seg_all: pd.DataFrame, 
            call_available: bool = False) -> Tuple[pd.DataFrame, int]:
    """
    Annotate cell data with segment means - RDR only version.
    
    Parameters
    ----------
    cell_data : pd.DataFrame
        Cell-level data
    seg_all : pd.DataFrame
        Segmentation data
    call_available : bool
        Whether to include copy number calls
        
    Returns
    -------
    Tuple[pd.DataFrame, int]
        Annotated cell data and number of breakpoints
    """
    cell_data['RDR_MEAN'] = np.nan
    
    df_new = []
    k = 0
    
    for chrom, df in cell_data.groupby('CHROM'):
        df.index = range(df.shape[0])
        seg = seg_all[seg_all.chrom == chrom]
        k = k + seg.shape[0] - 1
        
        for _, row in seg.iterrows():
            mask = (df.START <= row.end) & (df.END >= row.start)
            df.loc[mask, 'RDR_MEAN'] = row['RDR_MEAN']
            
            if call_available:
                df.loc[mask, 'CN'] = row['CN']
            
        df_new.append(df)
    
    df_new = pd.concat(df_new)
    if call_available:
        df_new['gamma'] = seg_all['gamma'].values[0]
        
    assert df_new.shape[0] == cell_data.shape[0]
    return df_new, k

def annotate_seg_rdr(seg_all: pd.DataFrame, cell_data: pd.DataFrame) -> pd.DataFrame:
    """
    Annotate segments with read depth statistics - RDR only version.
    
    Parameters
    ----------
    seg_all : pd.DataFrame
        Segmentation data
    cell_data : pd.DataFrame
        Cell-level data with read depths
        
    Returns
    -------
    pd.DataFrame
        Annotated segmentation data
    """
    cell_data['CHROM'] = cell_data['CHROM'].astype(str)
    seg_all['chrom'] = seg_all['chrom'].astype(str)
    rdr_sum, rdr_mean, n_bin = [], [], []
    
    for _, row in seg_all.iterrows():
        df = cell_data[cell_data['CHROM'] == row.chrom]
        _seg = df[(df.START <= row.end) & (df.END >= row.start)]
        
        n_bin.append(_seg.shape[0])
        rdr_sum.append(_seg.RDR.sum())
        rdr_mean.append(_seg.RDR.mean())
    
    seg_all.loc[:, 'NBIN'] = n_bin
    seg_all.loc[:, 'RDR_SUM'] = rdr_sum
    seg_all.loc[:, 'RDR_MEAN'] = rdr_mean
    
    return seg_all

def compute_bic_rdr(cell_data: pd.DataFrame, k: int) -> float:
    """
    Compute BIC for model selection - RDR only version.
    
    Parameters
    ----------
    cell_data : pd.DataFrame
        Cell data
    k : int
        Number of breakpoints
        
    Returns
    -------
    float
        BIC value
    """
    var_rdr = get_var(cell_data.RDR.values, cell_data.RDR_MEAN.values)
    bic = get_bic_gauss(np.log(var_rdr), k, cell_data.shape[0])
    return bic

def merge_adjacent_bins_all_rdr(seg_all: pd.DataFrame, cell_data: pd.DataFrame) -> pd.DataFrame:
    """
    Merge adjacent bins across all chromosomes - RDR only version.
    """
    updated_segs = []
    for chrom in seg_all.chrom.unique():
        seg = merge_adjacent_bins(seg_all[seg_all.chrom == chrom])
        seg['chrom'] = chrom
        updated_segs.append(seg)
    
    updated_segs = pd.concat(updated_segs)
    updated_segs.columns = ['start', 'end', 'CN', 'chrom']
    seg_all = annotate_seg_rdr(updated_segs, cell_data)
    return seg_all

def postprocess_rdr(seg_path: str, 
                   data_path: str, 
                   restrict_gamma: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Post-process segmentation results - RDR only version.
    
    Parameters
    ----------
    seg_path : str
        Path to segmentation file
    data_path : str
        Path to data file
    restrict_gamma : bool
        Whether to restrict gamma values
        
    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        Original segments, merged segments, and annotated cell data
    """
    seg_all = pd.read_csv(seg_path, sep='\t')
    cell_data = pd.read_csv(data_path, sep='\t')
    
    # Annotate segments
    seg_all = annotate_seg_rdr(seg_all, cell_data)
    seg_all = seg_all[seg_all.NBIN > 2]
    
    # Annotate cell data with segment mean
    cell_data, k = annotate_mean_rdr(cell_data, seg_all)
    cell_data = cell_data[~cell_data.RDR_MEAN.isna()]
    
    # Infer ploidy
    cn, gamma = infer_ploidy(seg_all, restrict_gamma=restrict_gamma)
    seg_all['CN'] = cn
    seg_all['gamma'] = gamma
    
    # Final annotations and merging
    cell_data, k = annotate_mean_rdr(cell_data, seg_all, call_available=True)
    seg_all_merged = merge_adjacent_bins_all_rdr(seg_all, cell_data)
    cell_data, k = annotate_mean_rdr(cell_data, seg_all_merged, call_available=False)
    
    return seg_all, seg_all_merged, cell_data


def infer_ploidy_rdr(cell_data: pd.DataFrame,
                     seg_all: pd.DataFrame,
                     max_wgd: int = 2,
                     restrict_gamma: bool = False) -> Tuple[List[int], float]:
    """
    Infer ploidy from read depth data only.
    
    Parameters
    ----------
    cell_data : pd.DataFrame
        Cell data with read depths
    seg_all : pd.DataFrame
        Segmentation data
    max_wgd : int
        Maximum whole genome duplication level
    restrict_gamma : bool
        Whether to restrict gamma values
        
    Returns
    -------
    Tuple[List[int], float]
        Copy numbers and gamma value
    """
    # Find modal RDR value with more robust binning
    rdr_values = cell_data.RDR[~cell_data.RDR.isna()].values
    if len(rdr_values) == 0:
        raise ValueError("No valid RDR values found")
        
    # Use a more robust method to find the modal RDR
    hist, bin_edges = np.histogram(rdr_values, bins=100)
    modal_idx = np.argmax(hist)
    modal_rdr = (bin_edges[modal_idx] + bin_edges[modal_idx + 1]) / 2
    
    # Add safety check for modal_rdr
    if modal_rdr <= 0:
        logger.warning("Modal RDR is zero or negative, using mean RDR instead")
        modal_rdr = np.mean(rdr_values)
        if modal_rdr <= 0:
            raise ValueError("Unable to determine valid RDR scaling factor")
    
    # Calculate gamma values with safety checks
    gamma_values = []
    for t in range(1, max_wgd + 1):
        gamma = 2 * t / modal_rdr
        if np.isfinite(gamma):  # Only include finite values
            gamma_values.append(gamma)
    
    if not gamma_values:
        raise ValueError("No valid gamma values could be calculated")
    
    if restrict_gamma:
        gamma_values = [g for g in gamma_values if 1.5 <= g <= 3.0]
        if not gamma_values:
            logger.warning("No gamma values within restricted range, using unrestricted values")
            gamma_values = [2 * t / modal_rdr for t in range(1, max_wgd + 1)]
    
    # Select gamma that minimizes variance
    min_var = float('inf')
    best_gamma = None
    best_cn = None
    
    for gamma in gamma_values:
        # Round RDR * gamma to nearest integer for copy numbers
        cn = np.round(seg_all.RDR_MEAN * gamma).astype(int)
        var = np.var(seg_all.RDR_MEAN * gamma - cn)
        
        if var < min_var:
            min_var = var
            best_gamma = gamma
            best_cn = cn.tolist()
    
    if best_cn is None or best_gamma is None:
        raise ValueError("Failed to determine copy numbers and gamma")
        
    return best_cn, best_gamma