import pandas as pd
import numpy as np
from scipy import stats
from sklearn.mixture import GaussianMixture
from collections import Counter
import math
from typing import Tuple, List, Optional

def annotate_seg(seg_all: pd.DataFrame, cell_data: pd.DataFrame) -> pd.DataFrame:
    """
    Annotate segments with read depth and BAF statistics.
    
    Parameters
    ----------
    seg_all : pd.DataFrame
        Segmentation data
    cell_data : pd.DataFrame
        Cell-level data with read depths and BAF values
        
    Returns
    -------
    pd.DataFrame
        Annotated segmentation data
    """
    rdr_sum, rdr_mean, vaf_mean, vaf_estimate, n_bin = [], [], [], [], []
    
    for _, row in seg_all.iterrows():
        df = cell_data[cell_data['#CHROM'] == row.chrom]
        _seg = df[(df.START <= row.end) & (df.END >= row.start)]
        
        n_bin.append(_seg.shape[0])
        rdr_sum.append(_seg.RDR.sum())
        rdr_mean.append(_seg.RDR.mean())
        
        _seg_snp = _seg[_seg.N > 10]
        
        if _seg_snp.shape[0] > 0:
            vaf_estimate.append(estimate_vaf(_seg_snp.A.values, _seg_snp.B.values, shift=2))
            vaf_mean.append(_seg.BAF.mean())
        else:
            vaf_estimate.append(np.nan)
            vaf_mean.append(np.nan)
    
    seg_all.loc[:, 'NBIN'] = n_bin
    seg_all.loc[:, 'RDR_SUM'] = rdr_sum
    seg_all.loc[:, 'RDR_MEAN'] = rdr_mean
    seg_all.loc[:, 'VAF_MEAN'] = vaf_mean
    seg_all.loc[:, 'VAF_ESTIMATE'] = vaf_estimate
    
    return seg_all

def annotate_mean(cell_data: pd.DataFrame, 
                 seg_all: pd.DataFrame, 
                 call_available: bool = False) -> Tuple[pd.DataFrame, int]:
    """
    Annotate cell data with segment means.
    
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
    cell_data['VAF_MEAN'] = np.nan
    cell_data['RDR_MEAN'] = np.nan
    cell_data['VAF_ESTIMATE'] = np.nan
    
    df_new = []
    k = 0
    
    for chrom, df in cell_data.groupby('#CHROM'):
        df.index = range(df.shape[0])
        seg = seg_all[seg_all.chrom == chrom]
        k = k + seg.shape[0] - 1
        
        for _, row in seg.iterrows():
            mask = (df.START <= row.end) & (df.END >= row.start)
            df.loc[mask, 'VAF_MEAN'] = row['VAF_MEAN']
            df.loc[mask, 'RDR_MEAN'] = row['RDR_MEAN']
            df.loc[mask, 'VAF_ESTIMATE'] = row['VAF_ESTIMATE']
            
            if call_available:
                df.loc[mask, 'CN'] = row['CN']
                df.loc[mask, 'pval'] = row['pval']
                df.loc[mask, 'prob'] = row['prob']
                df.loc[mask, 'CN_A'] = row['CN_A']
                df.loc[mask, 'CN_B'] = row['CN_B']
                df.loc[mask, 'CN_total'] = row['CN_total']
                
        df_new.append(df)
    
    df_new = pd.concat(df_new)
    if call_available:
        df_new['gamma'] = seg_all['gamma'].values[0]
        
    assert df_new.shape[0] == cell_data.shape[0]
    return df_new, k

def merge_adjacent_bins_all(seg_all: pd.DataFrame, cell_data: pd.DataFrame) -> pd.DataFrame:
    """
    Merge adjacent bins across all chromosomes.
    
    Parameters
    ----------
    seg_all : pd.DataFrame
        All segmentation data
    cell_data : pd.DataFrame
        Cell-level data
        
    Returns
    -------
    pd.DataFrame
        Merged segments for all chromosomes
    """
    updated_segs = []
    for chrom in seg_all.chrom.unique():
        seg = merge_adjacent_bins(seg_all[seg_all.chrom == chrom])
        seg['chrom'] = chrom
        updated_segs.append(seg)
    
    updated_segs = pd.concat(updated_segs)
    updated_segs.columns = ['start', 'end', 'CN', 'chrom']
    seg_all = annotate_seg(updated_segs, cell_data)
    return seg_all

def merge_adjacent_bins(seg: pd.DataFrame) -> pd.DataFrame:
    """
    Merge adjacent bins with the same copy number state.
    
    Parameters
    ----------
    seg : pd.DataFrame
        Segmentation data
        
    Returns
    -------
    pd.DataFrame
        Merged segments
    """
    df = seg[['start', 'end', 'CN']]
    df = df.sort_values("start")
    
    merged_segments = []
    current_start = df.iloc[0]["start"]
    current_end = df.iloc[0]["end"]
    current_cn = df.iloc[0]["CN"]
    
    for _, row in df.iterrows():
        if row["CN"] == current_cn:
            current_end = row["end"]
        else:
            merged_segments.append({
                "segment_start": current_start,
                "segment_end": current_end,
                "CN": current_cn
            })
            current_start = row["start"]
            current_end = row["end"]
            current_cn = row["CN"]
    
    merged_segments.append({
        "segment_start": current_start,
        "segment_end": current_end,
        "CN": current_cn
    })
    
    return pd.DataFrame(merged_segments)

def postprocess(seg_path: str, data_path: str, MAX_WGD: int = 1) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Post-process segmentation results.
    
    Parameters
    ----------
    seg_path : str
        Path to segmentation file
    data_path : str
        Path to data file
    MAX_WGD : int
        Maximum whole genome duplication level
        
    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        Original segments, merged segments, and annotated cell data
    """
    seg_all = pd.read_csv(seg_path, sep='\t')
    cell_data = pd.read_csv(data_path, sep='\t')
    
    # Annotate segments
    seg_all = annotate_seg(seg_all, cell_data)
    seg_all = seg_all[seg_all.NBIN > 2]
    
    # Annotate cell data with segment mean
    cell_data, k = annotate_mean(cell_data, seg_all)
    cell_data = cell_data[~cell_data.RDR_MEAN.isna()]
    cell_data = cell_data[~cell_data.VAF_ESTIMATE.isna()]
    
    # Identify allelic balanced clusters
    theta = np.arange(1, MAX_WGD + 1)
    AB_best, var_rdr, var_vaf = find_AB_cluster(cell_data, theta)
    gamma_pool = [2 * t * AB_best.shape[0] / AB_best.RDR.sum() for t in theta]
    
    # Infer ploidy
    bic, gamma, allele_cn, pval_statistic, allele_cn_prob = infer_ploidy(
        cell_data, seg_all, gamma_pool, var_rdr, var_vaf, LAMBDA=1
    )
    
    # Update segment information
    seg_all['CN'] = [f'{int(a)}|{int(b)}' for a, b in allele_cn]
    seg_all['pval'] = pval_statistic
    seg_all['prob'] = allele_cn_prob
    seg_all['gamma'] = gamma
    
    seg_all['CN_A'] = [int(cn.split('|')[0]) for cn in seg_all['CN']]
    seg_all['CN_B'] = [int(cn.split('|')[1]) for cn in seg_all['CN']]
    seg_all['CN_total'] = [int(a + b) for a, b in seg_all[['CN_A', 'CN_B']].values]
    
    # Final annotations and merging
    cell_data, k = annotate_mean(cell_data, seg_all, call_available=True)
    seg_all_merged = merge_adjacent_bins_all(seg_all, cell_data)
    cell_data, k = annotate_mean(cell_data, seg_all_merged, call_available=False)
    
    return seg_all, seg_all_merged, cell_data

def pairsum(target: int) -> Tuple[List[Tuple[int, int]], int]:
    """
    Find all pairs of integers that sum to target.
    
    Parameters
    ----------
    target : int
        Target sum
        
    Returns
    -------
    Tuple[List[Tuple[int, int]], int]
        List of pairs and count of pairs
    """
    vals = list(np.arange(0, target+1))
    vals.extend(vals)
    res = sorted([
        (a, b) for a, b in itertools.combinations(vals, 2) 
        if (a + b == target) and (a >= b)
    ])
    res = list(set(res))
    return res, len(res)

# Helper functions
def estimate_vaf(a: np.ndarray, b: np.ndarray, shift: int = 2) -> float:
    """Estimate VAF from allele counts."""
    vaf = a / (b + a + 1e-10)
    mvaf = [min(i, 1-i) for i in vaf]
    ap = np.roll(a, shift)
    vaf_rolled = ap / (b + ap + 1e-10)
    mvaf_rolled = [min(i, 1-i) for i in vaf_rolled]
    
    mvaf_sampled = np.random.choice(mvaf, min(100, len(mvaf)), replace=False)
    mvaf_rolled_sampled = np.random.choice(mvaf_rolled, min(100, len(mvaf)), replace=False)
    
    _, p = stats.ks_2samp(mvaf_sampled, mvaf_rolled_sampled)
    
    if (p > 0.05) and (np.mean(mvaf) < 0.2):
        return 0
    
    return get_peak(mvaf)

def get_peak(values: np.ndarray) -> float:
    """Find the highest peak in a distribution."""
    x, grid = np.histogram(values, bins=np.arange(-0.02, 0.54, 0.02))
    peaks, properties = stats.find_peaks(x, height=0)
    
    if len(peaks) == 0:
        return np.nan
        
    return grid[peaks[properties['peak_heights'].argmax()]]

def get_lnvar(x: np.ndarray, x_bar: np.ndarray) -> float:
    """
    Calculate log variance.
    
    Parameters
    ----------
    x : np.ndarray
        Data values
    x_bar : np.ndarray
        Mean values
        
    Returns
    -------
    float
        Log variance
    """
    n = x.shape[0]
    lnvar = np.log(np.sum((x - x_bar)**2)/n)
    return lnvar

def get_var(x: np.ndarray, x_bar: np.ndarray) -> float:
    """
    Calculate variance.
    
    Parameters
    ----------
    x : np.ndarray
        Data values
    x_bar : np.ndarray
        Mean values
        
    Returns
    -------
    float
        Variance
    """
    n = x.shape[0]
    var = np.sum((x - x_bar)**2)/n
    return var

def get_bic_gauss(lnvar: float, n_bp: int, n: int) -> float:
    """
    Calculate Bayesian Information Criterion for Gaussian model.
    
    Parameters
    ----------
    lnvar : float
        Log variance
    n_bp : int
        Number of breakpoints
    n : int
        Number of data points
        
    Returns
    -------
    float
        BIC value
    """
    return lnvar + n_bp/n * np.log(n)

def compute_bic(cell_data: pd.DataFrame, k: int, vaf_weight: float = 1, rdr_weight: float = 1) -> float:
    """
    Compute BIC for model selection.
    
    Parameters
    ----------
    cell_data : pd.DataFrame
        Cell data
    k : int
        Number of breakpoints
    vaf_weight : float
        Weight for VAF component
    rdr_weight : float
        Weight for RDR component
        
    Returns
    -------
    float
        BIC value
    """
    var_rdr = get_var(cell_data.RDR.values, cell_data.RDR_MEAN.values)
    cell_data_dn = cell_data[~cell_data.VAF_ESTIMATE.isna()]
    var_vaf = get_var(cell_data_dn.BAF.values, cell_data_dn.VAF_ESTIMATE.values)
    bic = get_bic_gauss(
        np.log(var_rdr) * vaf_weight + np.log(var_rdr) * rdr_weight, 
        k, 
        cell_data.shape[0]
    )
    return bic

def infer_ploidy(cell_data: pd.DataFrame, 
                 seg_all3: pd.DataFrame,
                 gamma_pool: np.ndarray,
                 var_rdr: float,
                 var_vaf: float,
                 LAMBDA: float = 1) -> Tuple[float, float, List[Tuple[int, int]], List[float], List[float]]:
    """
    Infer ploidy from cell data.
    
    Parameters
    ----------
    cell_data : pd.DataFrame
        Cell data
    seg_all3 : pd.DataFrame
        Segment data
    gamma_pool : np.ndarray
        Possible gamma values
    var_rdr : float
        RDR variance
    var_vaf : float
        VAF variance
    LAMBDA : float
        Lambda parameter
        
    Returns
    -------
    Tuple[float, float, List[Tuple[int, int]], List[float], List[float]]
        BIC, gamma, allele copy numbers, p-values, probabilities
    """
    epsilon = 1e-20
    assert 'VAF_ESTIMATE' in cell_data.columns
    
    BIC_list = {'BIC': [], 'gamma': [], 'k': [], 'llk': [],
                'allele_cn': [], 'lrt_pvals': [], 'allele_cn_prob': []}
    
    for gamma_p in gamma_pool:
        seg_llk = 0
        allele_cn_record = []
        allele_cn_prob_record = []
        lrt_pvals = []
        
        for _, seg in seg_all3.iterrows():
            total_cn_ceil = np.ceil(gamma_p * seg.RDR_MEAN)
            total_cn_floor = total_cn_ceil-1
            cn_pool = pairsum(total_cn_floor)[0] + pairsum(total_cn_ceil)[0]
            cn_pool_llk = []
            
            for a, b in cn_pool:
                llk = 0
                selected_bins = cell_data[
                    (cell_data['#CHROM']==seg.chrom) & 
                    (cell_data.START<=seg.end) & 
                    (cell_data.END>seg.start)
                ]
                
                for _, BIN in selected_bins.iterrows():
                    if not math.isnan(BIN['VAF_ESTIMATE']):
                        llk += stats.norm.logpdf(
                            BIN.RDR,
                            loc=(a+b)/gamma_p, 
                            scale=np.sqrt(var_rdr)
                        )
                        vaf_llk_corrected = stats.norm.logpdf(
                            BIN.VAF_ESTIMATE, 
                            loc=b/(a+b+1e-10), 
                            scale=np.sqrt(var_vaf)
                        )
                        llk += vaf_llk_corrected
                        
                cn_pool_llk.append(llk)
                
            allele_cn = cn_pool[np.argmax(cn_pool_llk)]
            best_llk = np.max(cn_pool_llk)
            
            max_llk = np.max(cn_pool_llk)
            cn_pool_prob = np.exp(cn_pool_llk - max_llk) / np.exp(stats.logsumexp(cn_pool_llk - max_llk))
            
            allele_cn = cn_pool[np.argmax(cn_pool_llk)]
            allele_cn_prob = cn_pool_prob[np.argmax(cn_pool_llk)]
            allele_cn_prob_record.append(allele_cn_prob)
            
            lrt_stats = 2 * (best_llk - np.array(cn_pool_llk))
            lrt_pvals_seg = stats.chi2.sf(lrt_stats, df=1)
            _, tippett_combined_pval = stats.combine_pvalues(lrt_pvals_seg, method='tippett')
            
            seg_llk += best_llk
            allele_cn_record.append(allele_cn)
            lrt_pvals.append(tippett_combined_pval)
        
        max_total_cn = int(np.ceil(cell_data.RDR.max() * gamma_p))
        k = 0
        for c in range(max_total_cn+1):
            k += np.floor(c/2)+1
            
        bic = LAMBDA * k * np.log(cell_data.shape[0]) - 2*seg_llk
        BIC_list['BIC'].append(bic)
        BIC_list['allele_cn'].append(allele_cn_record)
        BIC_list['allele_cn_prob'].append(allele_cn_prob_record)
        BIC_list['gamma'].append(gamma_p)
        BIC_list['k'].append(k)
        BIC_list['llk'].append(seg_llk)
        BIC_list['lrt_pvals'].append(lrt_pvals)
    
    index = np.argmin(BIC_list['BIC'])
    return (BIC_list['BIC'][index], BIC_list['gamma'][index],
            BIC_list['allele_cn'][index], BIC_list['lrt_pvals'][index],
            BIC_list['allele_cn_prob'][index])

def find_AB_cluster(cell_data: pd.DataFrame, theta: np.ndarray) -> Tuple[pd.DataFrame, float, float]:
    """Find allelic balanced clusters."""
    AB = cell_data[cell_data.VAF_ESTIMATE > 0.45]
    bic_max = np.inf
    AB_best = AB
    
    X = AB.RDR_MEAN.values.reshape(-1, 1)
    for t in theta:
        GMM = GaussianMixture(n_components=t, random_state=0).fit(X)
        y_hat = GMM.predict(X)
        cluster_id, _ = Counter(y_hat).most_common(1)[0]
        bic = GMM.bic(X)
        
        if bic < bic_max:
            bic_max = bic
            AB_best = AB[y_hat == cluster_id]
    
    var_rdr = np.var(AB_best.RDR.values - AB_best.RDR_MEAN.values)
    var_vaf = np.var(AB_best.BAF.values - AB_best.VAF_ESTIMATE.values)
    
    return AB_best, var_rdr, var_vaf