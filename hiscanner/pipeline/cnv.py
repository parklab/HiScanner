from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from ..logger import logger
from ..utils.preprocessing import prep_input_table
from ..utils.postprocessing import (
    annotate_seg,
    annotate_mean,
    find_AB_cluster,
    compute_bic,
    infer_ploidy,
    merge_adjacent_bins_all
)

def process_cell(cell: str,
                config: Dict[str, Any],
                chroms: List[str]) -> None:
    """Process a single cell for CNV calling."""
    try:
        final_call_dir = Path(config['final_call_dir'])
        final_call_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f'Processing cell {cell}')
        
        # Prepare input table
        cell_data = prep_input_table(
            cell,
            config['hetsnp_dir'],
            config['bin_dir'],
            chroms
        )
        cell_data.to_csv(final_call_dir / f'{cell}_input_table.txt', 
                        sep='\t', index=False)

        # Read segments
        logger.info(f'Reading segments for {cell} (lambda={config["lambda_value"]})')
        seg_all = pd.read_csv(
            f'{config["seg_dir"]}/{cell}/lambda{config["lambda_value"]}.cnv',
            sep='\t'
        )

        # Annotate segments
        logger.info('Annotating segments')
        seg_all = annotate_seg(seg_all, cell_data)
        seg_all = seg_all[seg_all.NBIN > 2]

        if seg_all.shape[0] == 0:
            logger.warning(f"No segments left for {cell} after filtering!")
            return

        # Process bins and segments
        cell_data, k = annotate_mean(cell_data, seg_all)
        cell_data = cell_data[~cell_data.RDR_MEAN.isna()]
        cell_data = cell_data[~cell_data.VAF_ESTIMATE.isna()]

        # Find allelic balanced clusters
        theta = np.arange(1, config['max_wgd'] + 1)
        logger.info('Finding allelic balanced clusters')
        AB_best, var_rdr, var_vaf = find_AB_cluster(cell_data, theta)

        if np.isnan(var_rdr) or np.isnan(var_vaf):
            logger.error(f"NaN detected in variance computation for {cell}")
            return

        # Prepare for ploidy inference
        gamma_pool = [2 * t * AB_best.shape[0] / AB_best.RDR.sum() for t in theta]
        bic = compute_bic(cell_data, k, vaf_weight=1, rdr_weight=1)

        # Infer ploidy
        logger.info('Inferring ploidy')
        bic, gamma, allele_cn, pval_statistic, allele_cn_prob = infer_ploidy(
            cell_data, seg_all, gamma_pool, var_rdr, var_vaf, LAMBDA=1
        )
        
        # Update segment information
        seg_all['CN'] = [f'{int(a)}|{int(b)}' for a, b in allele_cn]
        seg_all['pval'] = pval_statistic
        seg_all['prob'] = allele_cn_prob
        seg_all['gamma'] = gamma
        logger.info(f'Inferred ploidy: {gamma}')

        # Process final results
        logger.info('Merging segments')
        seg_all['CN_A'] = [int(cn.split('|')[0]) for cn in seg_all['CN']]
        seg_all['CN_B'] = [int(cn.split('|')[1]) for cn in seg_all['CN']]
        seg_all['CN_total'] = [int(a + b) for a, b in seg_all[['CN_A', 'CN_B']].values]
        
        # Final annotations and merging
        cell_data, k = annotate_mean(cell_data, seg_all, call_available=True)
        seg_all_merged = merge_adjacent_bins_all(seg_all, cell_data)
        cell_data, k = annotate_mean(cell_data, seg_all_merged, call_available=False)

        # Save results
        cell_data.to_csv(final_call_dir / f'{cell}.txt', sep='\t', index=False)
        seg_all_merged.to_csv(final_call_dir / f'{cell}_seg_merged.txt', sep='\t', index=False)
        seg_all.to_csv(final_call_dir / f'{cell}_seg.txt', sep='\t', index=False)

    except Exception as e:
        logger.error(f"Error processing cell {cell}: {e}")
        raise

def draw_track(cell: str, final_call_dir: Path, baf_alpha: float = 0.5) -> None:
    """Draw CNV and BAF tracks for a cell."""
    cell_data_path = final_call_dir / f'{cell}.txt'
    if not cell_data_path.exists():
        logger.error(f"Error: {cell_data_path} does not exist")
        return

    cell_data = pd.read_csv(cell_data_path, sep='\t')
    cell_data = cell_data[~cell_data['#CHROM'].isin(['X', 'Y'])]
    cell_data['#CHROM'] = cell_data['#CHROM'].astype(int)
    cell_data = cell_data.sort_values(by=['#CHROM', 'START'])

    # Create visualization
    fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 3), sharex=True, dpi=200)
    
    # Plot copy number track
    ax1.scatter(cell_data.index.values, 
               cell_data['RDR'].values * cell_data['gamma'].values[0],
               s=1, color='darkgrey')
    ax1.plot(cell_data.index.values, cell_data['CN_total'].values,
             color='black', lw=.5, alpha=1)
    
    # Plot BAF track
    ax2.scatter(cell_data.index.values, cell_data['pBAF'].values,
                s=1, color='darkgrey', alpha=baf_alpha)
    
    # Customize plot appearance
    ax2.set_ylim(-0.05, 1.05)
    ax1.set_ylim(0, 10)
    ax1.set_yticks([0, 2, 4, 6, 8, 10])
    ax2.set_ylabel('BAF')
    ax1.set_ylabel('Copy Number')
    
    plt.tight_layout()
    plt.savefig(final_call_dir / f'{cell}_track.png', dpi=300)
    plt.close()

def run_cnv_calling(config: Dict[str, Any]) -> None:
    """Run CNV calling pipeline."""
    logger.info("Starting CNV calling pipeline")
    
    try:
        # Get cell list
        if 'cell_file' in config:
            cells_df = pd.read_csv(config['cell_file'], sep='\t')
            if 'singlecell' in cells_df.columns:
                cells_df = cells_df.query('singlecell=="Y"')
            cells = cells_df['bamID'].tolist()
        elif 'cells' in config:
            cells = config['cells']
        else:
            raise ValueError("Either 'cells' or 'cell_file' must be provided")

        # Create output directory
        final_call_dir = Path(config['final_call_dir'])
        final_call_dir.mkdir(parents=True, exist_ok=True)

        # Process cells in parallel
        with ThreadPoolExecutor(max_workers=config.get('threads', 1)) as executor:
            futures = []
            for cell in cells:
                if config.get('rerun', False) or not (final_call_dir / f'{cell}.txt').exists():
                    futures.append(
                        executor.submit(process_cell, cell, config, config['chroms'])
                    )

            # Wait for all processes to complete
            for future in futures:
                future.result()

        # Generate visualizations
        logger.info("Generating visualizations")
        for cell in cells:
            draw_track(cell, final_call_dir, config.get('baf_alpha', 0.5))

        logger.info("CNV calling pipeline completed successfully")

    except Exception as e:
        logger.error(f"Error in CNV calling pipeline: {e}")
        raise