from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from ..logger import logger

# Import preprocessing utilities
from ..utils.preprocessing import (
    prep_input_table,
    prep_input_table_rdr_only
)

# Import standard pipeline postprocessing utilities
from ..utils.postprocessing import (
    annotate_seg,
    annotate_mean,
    find_AB_cluster,
    compute_bic,
    infer_ploidy,
    merge_adjacent_bins_all
)

# Import RDR-only pipeline postprocessing utilities
from ..utils.postprocessing_rdr import (
    annotate_seg_rdr,
    annotate_mean_rdr,
    compute_bic_rdr,
    merge_adjacent_bins_all_rdr,
    infer_ploidy_rdr
)


def process_cell(cell: str,
                config: Dict[str, Any],
                chroms: List[str]) -> None:
    try:
        final_call_dir = Path(config['outdir']) / 'final_calls'
        logger.info(f'Processing cell {cell}')
        input_table_path = final_call_dir / f'{cell}_input_table.txt'
        if not input_table_path.exists() or config.get('rerun', False):
            # Prepare input table
            cell_data = prep_input_table(
                cell,
                Path(config['outdir']) / 'phased_hets',
                Path(config['outdir']) / 'bins',
                chroms
            )
            cell_data.to_csv(input_table_path,
                            sep='\t', index=False)
        else:
            cell_data = pd.read_csv(input_table_path, sep='\t')
            
            
        # Read segments
        logger.info(f'Reading segments for {cell} (lambda={config["lambda_value"]})')
        segdir = Path(config['outdir']) / 'segs'
        seg_initial = pd.read_csv(f'{segdir}/{cell}/lambda{config["lambda_value"]}.cnv',
            sep='\t'
        )

        # Annotate segments
        logger.info('Annotating segments')
        seg_initial_anno = annotate_seg(seg_initial, cell_data)
        seg_all = seg_initial_anno[seg_initial_anno.NBIN > 2]

        if seg_all.shape[0] == 0:
            logger.warning(f"No segments left for {cell} after filtering!")
            # save seg_combined
            seg_initial_anno.to_csv(final_call_dir / f'{cell}_seg_initial_anno.txt', sep='\t', index=False)
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


def draw_track(cell: str, final_call_dir: Path, alpha: float = 0.5, rdr_only: bool = False) -> None:
    """
    Draw CNV and optionally BAF tracks for a cell with chromosome borders.
    
    Parameters
    ----------
    cell : str
        Cell identifier
    final_call_dir : Path
        Directory containing final call data
    alpha : float
        Alpha value for plot points
    rdr_only : bool
        If True, only draw RDR track without BAF
    """
    cell_data_path = final_call_dir / f'{cell}.txt'
    if not cell_data_path.exists():
        logger.error(f"Error: {cell_data_path} does not exist")
        return
    
    # Read and prepare data
    cell_data = pd.read_csv(cell_data_path, sep='\t')
    cell_data = cell_data[~cell_data['CHROM'].isin(['X', 'Y'])]
    cell_data['CHROM'] = cell_data['CHROM'].astype(int)
    
    # Sort by chromosome and position
    cell_data = cell_data.sort_values(['CHROM', 'START'])
    
    # Calculate chromosome boundaries
    chrom_boundaries = []
    current_pos = 0
    x_positions = []
    chrom_centers = []
    
    for chrom, group in cell_data.groupby('CHROM'):
        n_bins = len(group)
        x_positions.extend(range(current_pos, current_pos + n_bins))
        chrom_centers.append((current_pos + (current_pos + n_bins - 1)) / 2)
        current_pos += n_bins
        chrom_boundaries.append(current_pos - 0.5)
    
    # Create figure
    if rdr_only:
        fig, ax1 = plt.subplots(figsize=(15, 4))
        ax2 = None
    else:
        fig, (ax1, ax2) = plt.subplots(2, figsize=(15, 6), sharex=True, 
                gridspec_kw={'height_ratios': [2, 1]})
        fig.subplots_adjust(hspace=0.1)
    
    # Plot copy number track
    ax1.scatter(x_positions, cell_data['RDR'].values * cell_data['gamma'].values[0],
                s=1, color='darkgrey', alpha=alpha)
    
    # Use CN for RDR-only mode, CN_total for full mode
    cn_values = cell_data['CN'].values if rdr_only else cell_data['CN_total'].values
    ax1.plot(x_positions, cn_values, color='black', lw=0.5)
    
    # Plot BAF track if not in RDR-only mode
    if not rdr_only and 'pBAF' in cell_data.columns:
        ax2.scatter(x_positions, cell_data['pBAF'].values, 
                   s=1, color='darkgrey', alpha=alpha)
    
    # Add chromosome boundaries and labels
    for boundary in chrom_boundaries[:-1]:  # Don't add line after last chromosome
        ax1.axvline(x=boundary, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
        if not rdr_only:
            ax2.axvline(x=boundary, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
    
    # Add chromosome labels
    chroms = sorted(cell_data['CHROM'].unique())
    if rdr_only:
        ax1.set_xticks(chrom_centers)
        ax1.set_xticklabels(chroms)
    else:
        ax2.set_xticks(chrom_centers)
        ax2.set_xticklabels(chroms)
    
    # Customize plot appearance
    ax1.set_ylim(0, 10)
    ax1.set_yticks([0, 2, 4, 6, 8, 10])
    if not rdr_only:
        ax2.set_ylim(-0.05, 1.05)
    
    # Add labels
    ax1.set_ylabel('Copy Number')
    if rdr_only:
        ax1.set_xlabel('Chromosome')
    else:
        ax2.set_ylabel('BAF')
        ax2.set_xlabel('Chromosome')
    
    # Remove unnecessary spines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    if not rdr_only:
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
    
    plt.savefig(final_call_dir / f'{cell}_track.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    
def run_cnv_calling(config: Dict[str, Any]) -> None:
    """Run CNV calling pipeline."""
    logger.info("Starting CNV calling pipeline")
    
    try:
        # Determine if we're using RDR-only mode
        rdr_only = config.get('rdr_only', False)
        logger.info(f"Using {'RDR-only' if rdr_only else 'standard (RDR + BAF)'} mode")
        
        # Get cell list
        if 'metadata_path' in config:
            cells_df = pd.read_csv(config['metadata_path'], sep='\t')
            if 'singlecell' in cells_df.columns:
                cells_df = cells_df.query('singlecell=="Y"')
            else:
                raise ValueError("singlecell column not found in metadata file")
            cells = cells_df['bamID'].tolist()
        elif 'cells' in config:
            cells = config['cells']
        else:
            raise ValueError("Either 'cells' or 'metadata_path' must be provided")

        # Create output directory
        final_call_dir = Path(config['outdir']) / 'final_calls'
        logger.info(f"Saving final calls to {final_call_dir}")
        final_call_dir.mkdir(parents=True, exist_ok=True)

        # Process cells in parallel
        with ThreadPoolExecutor(max_workers=config.get('threads', 1)) as executor:
            futures = []
            for cell in cells:
                if config.get('rerun', False) or not (final_call_dir / f'{cell}.txt').exists():
                    if rdr_only:
                        from .cnv import process_cell_rdr as process_cell_func
                    else:
                        from .cnv import process_cell as process_cell_func
                    
                    futures.append(
                        executor.submit(
                            process_cell_func, 
                            cell, 
                            config, 
                            config['chrom_list']
                        )
                    )

            # Wait for all processes to complete
            for future in futures:
                future.result()

        # Generate visualizations
        logger.info("Generating visualizations")
        for cell in cells:
            draw_track(
                cell, 
                final_call_dir, 
                alpha=config.get('track_alpha', 0.5),
                rdr_only=rdr_only
            )

        logger.info("HiScanner CNV calling pipeline completed successfully")

    except Exception as e:
        logger.error(f"Error in CNV calling pipeline: {e}")
        raise
    

def process_cell_rdr(cell: str, config: Dict[str, Any], chroms: List[str]) -> None:
    """Process a single cell for CNV calling - RDR only version."""
    try:
        final_call_dir = Path(config['outdir']) / 'final_calls'
        logger.info(f'Processing cell {cell}')
        
        # Prepare input data
        input_table_path = final_call_dir / f'{cell}_input_table.txt'
        if not input_table_path.exists() or config.get('rerun', False):
            cell_data = prep_input_table_rdr_only(
                cell,
                Path(config['outdir']) / 'bins',
                chroms
            )
            cell_data.to_csv(input_table_path, sep='\t', index=False)
        else:
            cell_data = pd.read_csv(input_table_path, sep='\t')

        # Read segments
        logger.info(f'Reading segments for {cell} (lambda={config["lambda_value"]})')
        segdir = Path(config['outdir']) / 'segs'
        seg_all = pd.read_csv(f'{segdir}/{cell}/lambda{config["lambda_value"]}.cnv',
            sep='\t'
        )

        # Annotate segments
        logger.info('Annotating segments')
        seg_all = annotate_seg_rdr(seg_all, cell_data)
        seg_all = seg_all[seg_all.NBIN > 2]

        if seg_all.shape[0] == 0:
            logger.warning(f"No segments left for {cell} after filtering!")
            return

        # Process data
        cell_data, k = annotate_mean_rdr(cell_data, seg_all)
        cell_data = cell_data[~cell_data.RDR_MEAN.isna()]

        # Infer ploidy
        logger.info('Inferring ploidy')
        cn, gamma = infer_ploidy_rdr(
            cell_data,
            seg_all,
            max_wgd=config.get('max_wgd', 2),
            restrict_gamma=config.get('restrict_gamma', False)
        )
        logger.info(f'Inferred ploidy: {gamma}')
        
        # Update segment information
        seg_all['CN'] = cn
        seg_all['gamma'] = gamma
        
        # Process final results
        cell_data, k = annotate_mean_rdr(cell_data, seg_all, call_available=True)
        seg_all_merged = merge_adjacent_bins_all_rdr(seg_all, cell_data)
        cell_data, k = annotate_mean_rdr(cell_data, seg_all_merged, call_available=False)

        # Save results
        cell_data.to_csv(final_call_dir / f'{cell}.txt', sep='\t', index=False)
        seg_all_merged.to_csv(final_call_dir / f'{cell}_seg_merged.txt', sep='\t', index=False)
        seg_all.to_csv(final_call_dir / f'{cell}_seg.txt', sep='\t', index=False)

    except Exception as e:
        logger.error(f"Error processing cell {cell}: {e}")
        raise