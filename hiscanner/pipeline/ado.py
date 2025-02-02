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


def get_vaf_by_chrom(chrom, vaf, depth_filter=5, aggregate=False, k=1):
    """Get VAF data for a specific chromosome"""
    chrom = str(chrom)
    vaf['#CHROM'] = vaf['#CHROM'].astype(str)
    vaf = vaf[vaf['#CHROM'] == chrom]
    
    logger.debug(f'Number of hetSNPs in chromosome {chrom}: {vaf.shape[0]}')
    
    if aggregate:
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
    try:
        cell_dir = Path(config['ado_dir']) / cell
        cell_dir.mkdir(parents=True, exist_ok=True)
        
        # Load and process data
        fname = Path(config['hetsnp_dir']) / f'{cell}.hetsnp.txt'
        vaf = pd.read_csv(fname, sep='\t')
        logger.info(f'Total number of hetSNPs for {cell}: {vaf.shape[0]}')
        
        # Plot BAF distribution if requested
        if config.get('plot', False):
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
                                      aggregate=config.get('aggregate', False),
                                      k=config.get('k', 1))
            
            if selected.shape[0] == 0:
                logger.warning(f'No hetSNPs in chromosome {chrom} for {cell}')
                continue
def run_ado_analysis(config: Dict[str, Any]) -> None:
    """
    Run ADO (Allelic Dropout) analysis pipeline.
    
    Parameters
    ----------
    config : Dict[str, Any]
        Configuration dictionary containing analysis parameters
    """
    logger.info("Starting ADO analysis")
    
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
        if not min_bin_size_file.exists() or config.get('rerun', False):
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
            suggested_bin_size = min_bin_size_df.bin_size_min_kb.dropna().mean()

        logger.info(f'Suggested bin size for all cells: {suggested_bin_size} kb')

        # Create summary plot
        plt.figure(dpi=200)
        sns.boxplot(data=min_bin_size_df, x='bin_size_min_kb')
        plt.title(f'Suggested bin size: {np.round(suggested_bin_size)} kb | n = {n_cells}', 
                 fontsize=10)
        plt.savefig(ado_dir / 'min_bin_size_per_cell.png')
        plt.close()

        logger.info("ADO analysis completed successfully")

    except Exception as e:
        logger.error(f"Error in ADO analysis: {e}")
        raise