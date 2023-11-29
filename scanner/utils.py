
import os
from joblib import Parallel, delayed
from tqdm import tqdm
import argparse
from scipy.signal import find_peaks
from scipy.stats import ks_2samp as KS 
import heapq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from collections import Counter
import math
import itertools
from scipy.stats import norm, binom_test, beta

import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import json
import joblib
import contextlib

def pairsum(target):
    '''
    input: target sum
    output: the set of all pairs of integerts that add up to the given sum, and the cardinality of the set
    '''
    vals = list(np.arange(0,target+1))
    vals.extend(vals)
    res = sorted([(a, b) for a, b in itertools.combinations(vals, 2) if (a + b == target) and (a>=b)])
    res = list(set(res))
    return res, len(res)

def get_lnvar(x, x_bar):
    n = x.shape[0]
    lnvar = np.log(np.sum((x - x_bar)**2)/n)
    return lnvar

def get_var(x, x_bar):
    n = x.shape[0]
    var = np.sum((x - x_bar)**2)/n
    return var

def get_bic_gauss(lnvar, n_bp, n):
    
    return lnvar + n_bp/n * np.log(n)

def draw_whole_genome(path):
    plt.rcParams['axes.xmargin'] = 0
    # read calls
    print(path)
    df = pd.read_csv(path,sep='\t')
    bins = df[['CHROM', 'START', 'END']].copy(deep=True)
    df.reset_index(inplace=True)
    df.drop(columns=['index'], inplace=True)
    scale_factor = df['gamma'].values[0]

    fig,ax=plt.subplots(2, figsize=(8, 3),sharex=True, dpi=200)

    # plot dotted lines for chromosomes
    chrlines = [0]*len(bins)
    for i in range(1, len(bins)):
        if bins.iloc[i]['CHROM'] != bins.iloc[i-1]['CHROM']:
            chrlines[i] = 1
    # get index
    chrlines = [i for i, x in enumerate(chrlines) if x == 1]
    # plot rectangles for every other chromosome
    for index in range(len(chrlines)):
        if index % 2 == 0:
            if index+1 >= len(chrlines):
                break
            ax[0].axvspan(chrlines[index], chrlines[index+1], facecolor='grey', alpha=0.1)
            ax[1].axvspan(chrlines[index], chrlines[index+1], facecolor='grey', alpha=0.1)

    # plot data and calls
    ax[0].scatter(df.index.values, scale_factor*df['RDR'].values,s=1,color='darkgrey')
    ax[0].plot(df.index.values, df['CN_total'].values,color='black',lw=.5,alpha=1)
    ax[1].scatter(df.index.values, df['pBAF'].values,s=1,color='black',alpha=.1)


    ax[1].set_ylim(-0.05,1.05)
    ax[0].set_ylim(0,10)
    ax[0].set_yticks([0,2,4,6,8,10])
    ax[1].set_ylabel('BAF')
    ax[0].set_ylabel('Copy Number')    
    plt.tight_layout()

    # set x axis margins to 0
    ax[0].set_xlim(0, len(df))
    ax[1].set_xlim(0, len(df))


    # remove x axis ticks
    ax[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax[1].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    return fig,ax

def estimate_vaf(a,b,shift=2):
    '''
    1) KS test. If not signigicant, either homozygous or inconclusive (depending on mean baf)
    2) if KS test significant, find peak. If highest peak is 0, then find second highest peak.
    '''
    vaf = a/(b+a+1e-10)
    mvaf = [min(i,1-i) for i in vaf]
    ap = np.roll(a,shift)
    vaf_rolled = ap/(b+ap+1e-10)
    mvaf_rolled = [min(i,1-i) for i in vaf_rolled]
    mvaf_sampled = np.random.choice(mvaf, min(100, len(mvaf)),replace=False)
    mvaf_rolled_sampled = np.random.choice(mvaf_rolled, min(100, len(mvaf)),replace=False)
    _, p = KS(mvaf_sampled, mvaf_rolled_sampled, alternative='two-sided', mode='auto')
    x, grid = np.histogram(mvaf, bins=np.arange(-0.02,0.54,0.02))
    if (p>0.05) & (np.mean(mvaf)<.2):
        vaf_hat = 0
    else:
        vaf_hat = get_peak(x, grid)
    return vaf_hat

def get_peak(x, grid):
    '''
    Find the highest peak of a signal, excluding zero if possible.
    '''
    peaks, properties = find_peaks(x, height=0)
    if len(peaks) == 0:
        return np.nan
    return get_second_peak(peaks, properties['peak_heights'], grid)

def get_second_peak(peaks, heights, grid):
    '''
    Get the second highest peak when the first one is at zero.
    '''
    if grid[peaks].max() == 0 and len(heights) > 1:
        # Get the indices sorted by peak height
        sorted_indices = np.argsort(heights)
        # Return the second highest peak
        return grid[peaks[sorted_indices[-2]]]
    return grid[peaks[heights.argmax()]]





import scipy.stats as stats
from scipy.stats import combine_pvalues
from scipy.special import logsumexp
def infer_ploidy(cell_data,seg_all3,gamma_pool,var_rdr,var_vaf,LAMBDA=1):
    '''
    requires cell_data to have both VAF and VAF_ESTIMATE (VAF_ESTIMATE is estimated mirrored VAF)
    '''
    epsilon = 1e-20
    assert 'VAF_ESTIMATE' in cell_data.columns
    BIC_list = {'BIC':[],
                'gamma':[],
                'k':[],
                'llk':[],
                'allele_cn':[],
                'lrt_pvals':[],
                'allele_cn_prob':[]}
    for gamma_p in gamma_pool:
        seg_llk = 0
        allele_cn_record = []
        allele_cn_prob_record = []
        lrt_pvals = []
        for _, seg in seg_all3.iterrows():
            # for each segment, find best CN pair and corresponding log likelihood
            total_cn_ceil = np.ceil(gamma_p * seg.RDR_MEAN)
            total_cn_floor = total_cn_ceil-1
            cn_pool, _ = pairsum(total_cn_floor)
            cn_pool2, _ = pairsum(total_cn_ceil)
            cn_pool = cn_pool + cn_pool2
            cn_pool_llk = []
            # iterate over all possible combinations of allele-specific CN (a>=b)
            for index, (a,b) in enumerate(cn_pool):
                llk = 0
                # iterate over all bins in given segment
                selected_bins = cell_data[(cell_data['CHROM']==seg.chrom) & \
                     (cell_data.START<=seg.end) & (cell_data.END>seg.start)]
                for _, BIN in selected_bins.iterrows():
                    if not math.isnan(BIN['VAF_ESTIMATE']):
                        # RDR LLK (normal distribution)
                        llk += norm.logpdf(BIN.RDR,loc=(a+b)/gamma_p, scale=np.sqrt(var_rdr))
                        # VAF LLK 
                        vaf_llk_corrected= norm.logpdf(BIN.VAF_ESTIMATE, loc=b/(a+b+1e-10), scale=np.sqrt(var_vaf))
                        llk+=vaf_llk_corrected
                cn_pool_llk.append(llk)
            allele_cn = cn_pool[np.argmax(cn_pool_llk)]
            best_llk = np.max(cn_pool_llk)

             # normalize the log-likelihoods to get probabilities
            max_llk = np.max(cn_pool_llk)
            cn_pool_prob = np.exp(cn_pool_llk - max_llk) / np.exp(logsumexp(cn_pool_llk - max_llk))

            allele_cn = cn_pool[np.argmax(cn_pool_llk)]
            allele_cn_prob = cn_pool_prob[np.argmax(cn_pool_llk)]
            allele_cn_prob_record.append(allele_cn_prob)

            # get LRT p-values
            ## test statistic: ratio of the maximum ll to the ll under the alt CN. 
            ## this ratio follows a chi2 under the null hypothesis that the alt CN is true CN
            lrt_stats = 2 * (best_llk - np.array(cn_pool_llk))
            lrt_pvals_seg = stats.chi2.sf(lrt_stats, df=1)  # sf = 1 - cdf
            # combine p values 
            _, tippett_combined_pval = combine_pvalues(lrt_pvals_seg, method='tippett')
            # book-keeping
            seg_llk += best_llk
            allele_cn_record.append(allele_cn)
            lrt_pvals.append(tippett_combined_pval)

        max_total_cn = int(np.ceil(cell_data.RDR.max() * gamma_p))
        k = 0
        for c in range(max_total_cn+1):
            k+=np.floor(c/2)+1
        # bookkeeping
        bic = LAMBDA * k * np.log(cell_data.shape[0]) - 2*seg_llk
        BIC_list['BIC'].append(bic)
        BIC_list['allele_cn'].append(allele_cn_record)
        BIC_list['allele_cn_prob'].append(allele_cn_prob_record)
        BIC_list['gamma'].append(gamma_p)
        BIC_list['k'].append(k)
        BIC_list['llk'].append(seg_llk)
        BIC_list['lrt_pvals'].append(lrt_pvals)

    index = np.argmin(BIC_list['BIC'])
    return BIC_list['BIC'][index], BIC_list['gamma'][index],  BIC_list['allele_cn'][index], BIC_list['lrt_pvals'][index],BIC_list['allele_cn_prob'][index]



def annotate_mean(cell_data,seg_all,call_available=False):
    cell_data['VAF_MEAN'] = np.nan
    cell_data['RDR_MEAN'] = np.nan
    cell_data['VAF_ESTIMATE'] = np.nan

    df_new = []
    k = 0
    for chrom, df in cell_data.groupby('CHROM'):
        df.index = range(df.shape[0])
        # get breakpoints
        seg = seg_all[seg_all.chrom==chrom]
        k = k+seg.shape[0]-1
        for index, row in seg.iterrows():
            df.loc[(df.START<=row.end) & (df.END>=row.start),'VAF_MEAN'] = row['VAF_MEAN']
            df.loc[(df.START<=row.end) & (df.END>=row.start),'RDR_MEAN'] = row['RDR_MEAN']
            df.loc[(df.START<=row.end) & (df.END>=row.start),'VAF_ESTIMATE'] = row['VAF_ESTIMATE']
            if call_available:
                df.loc[(df.START<=row.end) & (df.END>=row.start),'CN'] = row['CN']
                df.loc[(df.START<=row.end) & (df.END>=row.start),'pval'] = row['pval']
                df.loc[(df.START<=row.end) & (df.END>=row.start),'prob'] = row['prob']
                df.loc[(df.START<=row.end) & (df.END>=row.start),'CN_A'] = row['CN_A']
                df.loc[(df.START<=row.end) & (df.END>=row.start),'CN_B'] = row['CN_B']
                df.loc[(df.START<=row.end) & (df.END>=row.start),'CN_total'] = row['CN_total']
        df_new.append(df)
    df_new = pd.concat(df_new)
    if call_available:
        df_new['gamma'] = seg_all['gamma'].values[0]   
    assert df_new.shape[0]==cell_data.shape[0]
    return df_new, k

def compute_bic(cell_data,k,vaf_weight=1,rdr_weight=1):
    var_rdr = get_var(cell_data.RDR.values, cell_data.RDR_MEAN.values)
    cell_data_dn = cell_data[~cell_data.VAF_ESTIMATE.isna()]
    var_vaf = get_var(cell_data_dn.BAF.values, cell_data_dn.VAF_ESTIMATE.values)
    bic = get_bic_gauss( np.log(var_rdr) * vaf_weight + np.log(var_rdr) * rdr_weight, k, cell_data.shape[0])
    return bic

def find_AB_cluster(cell_data, theta):
    AB = cell_data[cell_data.VAF_ESTIMATE>.45] 
    bic_max = np.inf 
    t_best = -1
    AB_best = AB
    X = AB.RDR_MEAN.values.reshape(-1,1)
    for t in theta:
        GMM = GaussianMixture(n_components=t, random_state=0).fit(X)
        y_hat = GMM.predict(X)
        cluster_id,_ = Counter(y_hat).most_common(1)[0] 
        bic = GMM.bic(X) 
        if bic<bic_max:
            bic_max = bic
            t_best = t
            AB_best = AB[y_hat==cluster_id]
    print(f'theta: {t_best}')

    var_rdr = get_var(AB_best.RDR.values, AB_best.RDR_MEAN.values)
    var_vaf = get_var(AB_best.BAF.values, AB_best.VAF_ESTIMATE.values)

    return AB_best, var_rdr, var_vaf

def merge_adjacent_bins_all(seg_all,cell_data):
    updated_segs = []
    for chrom in seg_all.chrom.unique():
        seg = merge_adjacent_bins(seg_all[seg_all.chrom==chrom])
        seg['chrom'] = chrom
        updated_segs.append(seg)
    updated_segs = pd.concat(updated_segs)
    updated_segs.columns = ['start','end','CN','chrom']
    seg_all = annotate_seg(updated_segs, cell_data)
    return seg_all

def annotate_seg(seg_all, cell_data):
    rdr_sum,rdr_mean,vaf_mean,vaf_estimate,n_bin = [],[],[],[],[]
    for _, row in seg_all.iterrows():
        df = cell_data[cell_data['CHROM']==row.chrom]
        _seg = df[(df.START<=row.end) & (df.END>=row.start)]

        n_bin.append(_seg.shape[0])
        rdr_sum.append(_seg.RDR.sum())
        rdr_mean.append(_seg.RDR.mean())
        
        _seg_snp = _seg[_seg.N>10]
        
        if _seg_snp.shape[0]>0:
            vaf_estimate.append(estimate_vaf(_seg_snp.A.values,_seg_snp.B.values, shift=2))
            vaf_mean.append(_seg.BAF.mean())
        else: 
            vaf_estimate.append(np.nan)
            vaf_mean.append(np.nan)

    seg_all.loc[:,'NBIN'] = n_bin
    seg_all.loc[:,'RDR_SUM'] = rdr_sum
    seg_all.loc[:,'RDR_MEAN'] = rdr_mean
    seg_all.loc[:,'VAF_MEAN'] = vaf_mean
    seg_all.loc[:,'VAF_ESTIMATE'] = vaf_estimate
    return seg_all

def merge_adjacent_bins(seg):
    df = seg[['start','end','CN']]
    # Sort the dataframe by the 'start' column
    df = df.sort_values('start')
    # Initialize variables to store the merged segments
    merged_segments = []
    current_start = df.iloc[0]['start']
    current_end = df.iloc[0]['end']
    current_cn = df.iloc[0]['CN']
    
    # Iterate through each row in the dataframe
    for index, row in df.iterrows():
        # Check if the current row has the same 'cn' as the previous row
        if row['CN'] == current_cn:
            # Update the current segment's end position
            current_end = row['end']
        else:
            # Add the current segment to the merged_segments list
            merged_segments.append({'segment_start': current_start, 'segment_end': current_end, 'CN': current_cn})
            
            # Update the current segment's start and end positions and 'cn' value
            current_start = row['start']
            current_end = row['end']
            current_cn = row['CN']
    
    # Add the last segment to the merged_segments list
    merged_segments.append({'segment_start': current_start, 'segment_end': current_end, 'CN': current_cn})
    
    # Create a new dataframe from the merged segments list
    merged_df = pd.DataFrame(merged_segments)
    
    return merged_df


def postprocess(seg_path,data_path,MAX_WGD=1):     
    seg_all = pd.read_csv(seg_path,sep='\t')
    cell_data = pd.read_csv(data_path,sep='\t')
    
    rdr_sum,rdr_mean,vaf_mean,n_bin = [],[],[],[]

    # annotate segments
    seg_all = annotate_seg(seg_all, cell_data)
    seg_all = seg_all[seg_all.NBIN>2]

    # annotate cell data with segment mean
    cell_data, k = annotate_mean(cell_data,seg_all)
    cell_data = cell_data[~cell_data.RDR_MEAN.isna()]
    cell_data = cell_data[~cell_data.VAF_ESTIMATE.isna()]
    theta = np.arange(1,MAX_WGD+1) # increase if more WGD expected
    
    # identify allelic balanced clusters
    AB_best, var_rdr, var_vaf = find_AB_cluster(cell_data, theta)
    gamma_pool = [2 * t * AB_best.shape[0] / AB_best.RDR.sum() for t in theta]
    bic = compute_bic(cell_data,k,vaf_weight=1,rdr_weight=1)
    
    assert not math.isnan(var_rdr)
    assert not math.isnan(var_vaf)

    bic, gamma, allele_cn, pval_statistic, allele_cn_prob = infer_ploidy(cell_data,seg_all,gamma_pool,var_rdr,var_vaf,LAMBDA=1)
    seg_all['CN'] = [f'{int(a)}|{int(b)}' for a,b in allele_cn]
    seg_all['pval'] = pval_statistic
    seg_all['prob'] = allele_cn_prob
    seg_all['gamma'] = gamma

    # merge segment by CN and re-annotate
   
    seg_all['CN_A'] = [int(cn.split('|')[0]) for cn in seg_all['CN']]
    seg_all['CN_B'] = [int(cn.split('|')[1]) for cn in seg_all['CN']]
    seg_all['CN_total'] = [int(a+b) for a,b in seg_all[['CN_A','CN_B']].values]

    cell_data, k = annotate_mean(cell_data,seg_all,call_available=True)
    
    seg_all_merged = merge_adjacent_bins_all(seg_all,cell_data)
    cell_data, k = annotate_mean(cell_data,seg_all_merged,call_available=False)

    return seg_all,seg_all_merged, cell_data

def check_file(fn, overwrite):
    if os.path.exists(fn):
        if not overwrite:
            print(f'{fn} exists, skipping')
            return True
        else:
            print(f'{fn} exists, but overwriting')
    return False

@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    '''
    Context manager to patch joblib to report into tqdm progress bar given as argument
    parallelization code adapted from https://stackoverflow.com/questions/24983493/tracking-progress-of-joblib-parallel-execution
    '''
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()
