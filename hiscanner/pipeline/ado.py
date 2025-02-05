import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, Any
from ..logger import logger
from hmmlearn import hmm
import pickle
import scipy
import os

def get_vaf_by_chrom(chrom, vaf, depth_filter=5, aggregate_every_k_snp=False, k=1):
    """
    Get VAF data for a specific chromosome
    
    Parameters
    ----------
    chrom : str
        Chromosome identifier
    vaf : pd.DataFrame
        VAF data
    depth_filter : int, optional
        Minimum depth filter, by default 5
    aggregate_every_k_snp : bool, optional
        Whether to aggregate every k SNPs, by default False
    k : int, optional
        Number of SNPs to aggregate when aggregate_every_k_snp is true, by default 1
    
    Returns
    -------
    pd.DataFrame
        Filtered and potentially aggregated VAF data
    """
    vaf['CHROM'] = vaf['CHROM'].astype(str)
    vaf = vaf[vaf['CHROM'] == chrom]
    
    logger.debug(f'Number of hetSNPs in chromosome {chrom}: {vaf.shape[0]}')
    
    if aggregate_every_k_snp:
        vaf_agg = vaf[['A', 'B']].groupby(vaf.index // k).sum()
        vaf_agg['START'] = vaf['POS'].groupby(vaf.index // k).first()
        vaf_agg['END'] = vaf['POS'].groupby(vaf.index // k).last()
        vaf_agg['TOTAL'] = vaf_agg['A'] + vaf_agg['B']
        vaf_agg['pBAF'] = vaf_agg['B'] / vaf_agg['TOTAL']
        vaf_agg['BAF'] = vaf_agg['pBAF'].apply(lambda x: min(x, 1-x))
        vaf_agg.reset_index(drop=True, inplace=True)
        vaf = vaf_agg
    return vaf[vaf.TOTAL > depth_filter]

def process_cell_ado(cell, config, metadata_all):
    """Process ADO analysis for a single cell"""
    if not os.path.exists(config['hetsnp_dir']):
        raise FileNotFoundError(f"Heterozygous SNP directory not found: {config['hetsnp_dir']}")

    fname = Path(config['hetsnp_dir']) / f'{cell}.hetsnp.txt'
    if not fname.exists():
        raise FileNotFoundError(f"Heterozygous SNP file not found: {fname}")
    try:
        cell_dir = Path(config['ado_dir']) / cell
        cell_dir.mkdir(parents=True, exist_ok=True)
        
        # Load and process data
        fname = Path(config['hetsnp_dir']) / f'{cell}.hetsnp.txt'
        vaf = pd.read_csv(fname, sep='\t', low_memory=False)
        vaf = vaf.rename(columns={'#CHROM': 'CHROM'})
        logger.info(f'Total number of hetSNPs for {cell}: {vaf.shape[0]}')
        
        # Plot BAF distribution if requested
        if config.get('ado_plot_baf_distribution', False):
            if not (cell_dir / 'pBAF_distribution.png').exists() or config.get('rerun', False):
                logger.info(f'Plotting pBAF distribution for {cell}')
                plot_baf_distribution(vaf, cell_dir, config)
        
        # Process data for HMM
        if not (cell_dir / 'data.npy').exists() or config.get('rerun', False):
            data = process_cell_data(vaf, config)
            np.save(cell_dir / 'data.npy', data)
        else:
            logger.info(f'Loading existing data for {cell}')
            data = np.load(cell_dir / 'data.npy')
        
        # Train HMM model
        if not (cell_dir / 'hmm_model.pkl').exists() or config.get('rerun', False):
            train_hmm_model(data, cell_dir, config)
        
        # Process results
        process_cell_results(cell, vaf, cell_dir, config)
        
    except Exception as e:
        logger.error(f"Error processing cell {cell}: {e}")
        raise

def train_hmm_model(data, cell_dir, config):
    """Train HMM model for ADO detection"""
    logger.info('Training HMM model')
    
    data_binary = data > config.get('ado_threshold', 0.2)
    data_binary = np.reshape(data_binary, (-1, 1))
    
    model = hmm.GaussianHMM(n_components=2, covariance_type="diag", verbose=False)
    model.fit(data_binary)
    
    # Save model parameters
    mus = np.ravel(model.means_)
    sigmas = np.ravel(np.sqrt([np.diag(c) for c in model.covars_]))
    P = model.transmat_
    
    logger.info(f'Model parameters - Means: {mus}, Sigmas: {sigmas}, Transition matrix: {P}')
    
    with open(cell_dir / 'hmm_model.pkl', 'wb') as f:
        pickle.dump(model, f)

def process_cell_results(cell, vaf, cell_dir, config):
    """Process and save cell results"""
    model = pickle.load(open(cell_dir / 'hmm_model.pkl', 'rb'))
    
    if not (cell_dir / 'result_gauss.txt').exists() or config.get('rerun', False):
        anno_df_combined = []
        
        for chrom in range(1, 23):
            selected = get_vaf_by_chrom(str(chrom), vaf, 
                        depth_filter=config.get('depth_filter', 5),
                        aggregate_every_k_snp=config.get('aggregate_every_k_snp', False),
                        k=config.get('k', 1))
            
            if selected.shape[0] == 0:
                logger.warning(f'No hetSNPs in chromosome {chrom} for {cell}')
                continue
                
            sequence = selected.BAF.values
            sequence_binary = sequence > config.get('ado_threshold', 0.2)
            sequence_reshaped = sequence_binary.reshape(-1, 1)
            
            _, decoded = model.decode(sequence_reshaped)
            
            # Determine whether states are flipped
            corr = np.corrcoef(decoded.ravel(), sequence_reshaped.ravel())[0, 1]
            convert = {0: 1, 1: 0} if corr > 0 else {0: 0, 1: 1}
            sequence_converted = [convert[d] for d in decoded]
            
            seq_region_indices = scipy.ndimage.find_objects(scipy.ndimage.label(sequence_converted)[0])
            
            anno = []
            if config.get('aggregate_every_k_snp', False):
                START = np.array(selected.START)
                END = np.array(selected.END)
                for region in seq_region_indices:
                    start_idx, end_idx = region[0].start, region[0].stop
                    if start_idx == end_idx - 1:
                        continue
                    start_pos, end_pos = START[start_idx], END[end_idx-1]
                    state = decoded[start_idx]
                    length = end_pos - start_pos
                    baf = sequence[start_idx:end_idx].mean()
                    anno.append([start_pos, end_pos, 
                        decoded[region].shape[0], length, 
                        state, baf])
            else:
                pos = np.array(selected.POS)
                for region in seq_region_indices:
                    start_idx, end_idx = region[0].start, region[0].stop
                    if start_idx == end_idx - 1:
                        continue
                    start_pos, end_pos = pos[start_idx], pos[end_idx-1]
                    state = decoded[start_idx]
                    length = end_pos - start_pos
                    baf = sequence[start_idx:end_idx].mean()
                    anno.append([start_pos, end_pos, 
                            decoded[region].shape[0], length, 
                            state, baf])
            
            anno_df = pd.DataFrame(anno, 
                            columns=['start', 'end', 'n', 'length', 
                            'state', 'mean_baf'])
            anno_df['chrom'] = chrom
            anno_df_combined.append(anno_df)
            
        anno_df_combined = pd.concat(anno_df_combined)
        anno_df_combined.to_csv(cell_dir / 'result_gauss.txt', sep='\t', index=None)
    else:
        anno_df_combined = pd.read_csv(cell_dir / 'result_gauss.txt', sep='\t')
        
    # Calculate minimum bin size
    anno_df_combined = anno_df_combined[anno_df_combined.length > 1]
    anno_df_combined['log10_len'] = np.log10(anno_df_combined['length'].astype(float))
    
    bin_size_min = anno_df_combined.length.quantile(q=.99)
    bin_size_min_round_kb = np.ceil(bin_size_min/1000)
    
    logger.info(f'Suggested minimum bin size (kb): {bin_size_min_round_kb}')
    
    # Save suggested bin size
    with open(cell_dir / 'bin_size_min.txt', 'w') as f:
        f.write(f'{bin_size_min_round_kb} kb')
        
    # Create plot if requested
    if config.get('plot', False):
        plot_path = cell_dir / 'ado_length_distribution_per_chrom.png'
        if not plot_path.exists() or config.get('rerun', False):
            fig, ax = plt.subplots(dpi=500, figsize=(4, 3))
            for chrom in range(1, 23):
                sns.kdeplot(data=anno_df_combined[anno_df_combined.chrom == chrom],
                        x='log10_len', ax=ax, color='black', lw=.5)
            plt.axvline(np.log10(bin_size_min_round_kb*1000), 
                        label='Suggested bin size', color='red', lw=1)
            plt.xlim(0, 9)
            plt.legend()
            plt.xlabel('Allelic dropout length (log10 scale)')
            plt.savefig(plot_path)
            plt.close()


def run_ado_analysis(config: Dict[str, Any]) -> None:
    """
    Run ADO (Allelic Dropout) analysis pipeline.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Configuration dictionary containing analysis parameters
    """
    logger.info("Starting ADO analysis")
    
    if 'outdir' not in config:
        raise ValueError("Output directory not specified in configuration")

    ado_dir = Path(config['outdir']) / 'ado'
    config['ado_dir'] = str(ado_dir)  # Add this to config
    config['hetsnp_dir'] = str(Path(config['outdir']) / 'phased_hets')
    
    try:
        # Get cell list
        metadata_path = config['metadata_path']
        metadata_all = pd.read_csv(metadata_path, sep='\t')
        metadata_all = metadata_all[metadata_all.singlecell == 'Y']

        cells = config.get('cells', 'all')
        if cells == 'all':
            cells = metadata_all.bamID.values
        else:
            cells = cells.split(',')

        n_cells = len(cells)
        logger.info(f'Performing ADO analysis for {n_cells} cells')

        # Process each cell
        ado_dir = Path(config['ado_dir'])
        ado_dir.mkdir(parents=True, exist_ok=True)

        for cell in cells:
            try:
                process_cell_ado(cell, config, metadata_all)
            except Exception as e:
                logger.error(f"Error processing cell {cell}: {e}")
                continue

        # Gather results
        min_bin_size_file = ado_dir / 'min_bin_size_per_cell.txt'
        min_bin_sizes = []
        for cell in cells:
            cell_dir = ado_dir / cell
            try:
                with open(cell_dir / 'bin_size_min.txt', 'r') as f:
                    bin_size = float(f.read().split(' ')[0])
                min_bin_sizes.append([bin_size, cell])
            except Exception as e:
                logger.warning(f"Could not read bin size for {cell}: {e}")
                continue

            min_bin_size_df = pd.DataFrame(min_bin_sizes, columns=['bin_size_min_kb', 'cell'])
            min_bin_size_df.to_csv(min_bin_size_file, sep='\t', index=None)
            suggested_bin_size = min_bin_size_df.dropna().bin_size_min_kb.mean()
        else:
            min_bin_size_df = pd.read_csv(min_bin_size_file, sep='\t')
            suggested_bin_size = min_bin_size_df.bin_size_min_kb.dropna().max()
            suggested_bin_size = np.round(suggested_bin_size,1)
        logger.info('=' * 80)
        logger.info(f'SUGGESTED BIN SIZE FOR ALL CELLS: at least {suggested_bin_size} kb.')
        logger.info(f'Please edit the bin size in the config file to the suggested value.')
        logger.info('=' * 80)


        # Create summary plot
        plt.figure(dpi=200)
        sns.boxplot(data=min_bin_size_df, x='bin_size_min_kb')
        plt.title(f'Suggested bin size: {np.round(suggested_bin_size)} kb | n = {n_cells}', fontsize=10)
        plt.savefig(ado_dir / 'min_bin_size_per_cell.png')
        plt.close()

        logger.info("ADO analysis completed successfully")

    except Exception as e:
        logger.error(f"Error in ADO analysis: {e}")
        raise
def plot_baf_distribution(vaf, cell_dir, config):
    """Plot distribution of pBAF in heterozygous SNPs"""
    logger.info('Plotting pBAF distribution')
    fig, ax = plt.subplots(dpi=200)
    if config.get('aggregate_every_k_snp', False):
        k = config.get('k', 1)
        logger.info(f'Aggregating every {k} SNPs')
        vaf_agg = []
        for chrom_train in range(1, 23):
            selected = get_vaf_by_chrom(
                str(chrom_train), 
                vaf, 
                depth_filter=config.get('depth_filter', 5),
                aggregate_every_k_snp=True,
                k=config.get('k', 1)
            )
            vaf_agg.append(selected)
        vaf = pd.concat(vaf_agg)
    else:
        vaf = vaf[vaf.TOTAL > config.get('depth_filter', 5)]
    
    sns.kdeplot(vaf.pBAF, bw_adjust=0.5, fill=False, color='black', ax=ax, lw=0.5)
    plt.savefig(cell_dir / 'BAF_distribution.png')
    plt.close()

def process_cell_data(vaf, config):
    """Process cell data for HMM analysis"""
    data = []
    vaf['CHROM'] = vaf['CHROM'].astype(str)
    
    for chrom_train in range(1, 23):
        selected = get_vaf_by_chrom(
            str(chrom_train), 
            vaf, 
            depth_filter=config.get('depth_filter', 5),
            aggregate_every_k_snp=config.get('aggregate_every_k_snp', False),
            k=config.get('k', 1)
        )
        sequence = selected.BAF.values
        sequence_reshaped = sequence.reshape(-1, 1)
        data.append(sequence_reshaped)
    
    return np.concatenate(data)