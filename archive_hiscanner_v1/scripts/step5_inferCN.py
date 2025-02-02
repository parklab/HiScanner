import os
import pandas as pd
import numpy as np
import argparse
import subprocess
import math
from postprocess_utils import *
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor

def draw(cell, final_call_dir, baf_alpha):
        cell_data_path = f'{final_call_dir}/{cell}.txt'
        if not os.path.exists(cell_data_path):
            print(f"Error: {cell_data_path} does not exist")
            return
        else:
            cell_data = pd.read_csv(cell_data_path, sep='\t')
            cell_data = cell_data[cell_data['#CHROM'] != 'X']
            cell_data = cell_data[cell_data['#CHROM'] != 'Y']
            cell_data['#CHROM'] = cell_data['#CHROM'].astype(int)
            cell_data = cell_data.sort_values(by=['#CHROM', 'START'])
            fig, ax = draw_track(cell_data, baf_alpha=baf_alpha)
            plt.savefig(f'{final_call_dir}/{cell}_track.png', dpi=300)
            return
        
def process_cell(cell, hetsnp_dir, bin_dir, seg_dir, final_call_dir, LAMBDA, MAX_WGD, chroms):
    """
    Process a single cell's data: preparation, annotation, CNV segmentation, and final plotting.
    """
    if not os.path.exists(final_call_dir):
        os.makedirs(final_call_dir, exist_ok=True)
    print(f'Processing cell {cell}')
    try:
        cell_data = prep_input_table(cell, hetsnp_dir, bin_dir, chroms)
        print('Finished preparing input table')
        cell_data.to_csv(f'{final_call_dir}/{cell}_input_table.txt', sep='\t', index=False)

        print(f'Reading segments for {cell} (lambda={LAMBDA})')
        seg_all = pd.read_csv(f'{seg_dir}/{cell}/lambda{LAMBDA}.cnv', sep='\t')

        print('Annotating segments')
        seg_all = annotate_seg(seg_all, cell_data)
        seg_all = seg_all[seg_all.NBIN > 2]

        if seg_all.shape[0] == 0:
            print(f"Warning: No segments left for {cell} after filtering!")
            return

        print('Annotating bins')
        cell_data, k = annotate_mean(cell_data, seg_all)
        cell_data = cell_data[~cell_data.RDR_MEAN.isna()]
        cell_data = cell_data[~cell_data.VAF_ESTIMATE.isna()]

        theta = np.arange(1, MAX_WGD + 1)
        print('Finding allelic balanced clusters')
        AB_best, var_rdr, var_vaf = find_AB_cluster(cell_data, theta)

        if math.isnan(var_rdr) or math.isnan(var_vaf):
            print(f"Error: NaN detected in variance computation for {cell}")
            return

        gamma_pool = [2 * t * AB_best.shape[0] / AB_best.RDR.sum() for t in theta]
        bic = compute_bic(cell_data, k, vaf_weight=1, rdr_weight=1)

        print('Inferring ploidy')
        bic, gamma, allele_cn, pval_statistic, allele_cn_prob = infer_ploidy(
            cell_data, seg_all, gamma_pool, var_rdr, var_vaf, LAMBDA=1)
        seg_all['CN'] = [f'{int(a)}|{int(b)}' for a, b in allele_cn]
        seg_all['pval'] = pval_statistic
        seg_all['prob'] = allele_cn_prob
        seg_all['gamma'] = gamma
        print(f'Inferred ploidy: {gamma}')

        print('Merging segments')
        seg_all['CN_A'] = [int(cn.split('|')[0]) for cn in seg_all['CN']]
        seg_all['CN_B'] = [int(cn.split('|')[1]) for cn in seg_all['CN']]
        seg_all['CN_total'] = [int(a + b) for a, b in seg_all[['CN_A', 'CN_B']].values]
        cell_data, k = annotate_mean(cell_data, seg_all, call_available=True)
        seg_all_merged = merge_adjacent_bins_all(seg_all, cell_data)
        cell_data, k = annotate_mean(cell_data, seg_all_merged, call_available=False)

        cell_data.to_csv(f'{final_call_dir}/{cell}.txt', sep='\t', index=False)
        seg_all_merged.to_csv(f'{final_call_dir}/{cell}_seg_merged.txt', sep='\t', index=False)
        seg_all.to_csv(f'{final_call_dir}/{cell}_seg.txt', sep='\t', index=False)

    except Exception as e:
        print(f"Error processing cell {cell}: {e}")

def process_batch(batch, hetsnp_dir, bin_dir, seg_dir, final_call_dir, LAMBDA, MAX_WGD, chroms, rerun):
    """
    Process a batch of cells.
    """
    for cell in batch:
        # Skip if the final call file already exists and we're not overwriting
        if not rerun and os.path.exists(f'{final_call_dir}/{cell}.txt'):
            print(f'Skipping {cell} because it already exists')
            continue
        process_cell(cell, hetsnp_dir, bin_dir, seg_dir, final_call_dir, LAMBDA, MAX_WGD, chroms)

def main():
    parser = argparse.ArgumentParser(description='Process batches of cells for CNV analysis.')
    parser.add_argument('--batch_size', type=int, default=1, help='Number of cells to process per batch')
    parser.add_argument('--hetsnp_dir', required=True, help='Directory containing phased hetsnp data')
    parser.add_argument('--bin_dir', required=True, help='Directory containing bin data')
    parser.add_argument('--seg_dir', required=True, help='Directory containing segment data')
    parser.add_argument('--final_call_dir', required=True, help='Output directory for final call data')
    parser.add_argument('--lambda_value', type=int, default=16, help='Lambda value for segment data')
    parser.add_argument('--max_wgd', type=int, default=1, help='Maximum Whole Genome Duplication level')
    parser.add_argument('--chroms', nargs='+', type=str, default='22', help='List of chromosomes to process')
    parser.add_argument('--cells', nargs='+', help='List of cell names to process')
    parser.add_argument('--cell_file', help='File containing a list of cells in a column called "bamID"')
    parser.add_argument('--rerun', action='store_true', help='Overwrite existing files')
    parser.add_argument('--draw_only', action='store_true', help='Only draw tracks for existing files')
    parser.add_argument('--baf_alpha', type=float, default=0.5, help='Alpha value for BAF plot transparency')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use for parallel processing')

    args = parser.parse_args()

    
    # Read cell names from file if provided
    if args.cell_file:
        cells_df = pd.read_csv(args.cell_file, sep='\t')  # Assume the file is tab-separated
        if 'singlecell' in cells_df.columns:
            cells_df = cells_df.query('singlecell=="Y"')
        else:
            print("Warning: No 'singlecell' column found in cell file. Assuming all cells are single-cell.")
        cells = cells_df['bamID'].tolist()
    elif args.cells:
        cells = args.cells
    else:
        print("Error: Either --cells or --cell_file must be provided.")
        return

    batch_size = args.batch_size
    hetsnp_dir = args.hetsnp_dir
    bin_dir = args.bin_dir
    seg_dir = args.seg_dir
    final_call_dir = args.final_call_dir
    LAMBDA = args.lambda_value
    MAX_WGD = args.max_wgd
    chroms = args.chroms
    rerun = args.rerun
    draw_only = args.draw_only
    baf_alpha = args.baf_alpha
    num_threads = args.threads

    
    if draw_only:
        # simply loop through the cells and draw the tracks
        for cell in cells:
            print(f'Drawing track for {cell}')
            draw(cell, final_call_dir, baf_alpha)
            plt.clf()
        return
    
    # Split the cells into batches
    batches = [cells[i:i + batch_size] for i in range(0, len(cells), batch_size)]

    # Process batches in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(process_batch, batch, hetsnp_dir, bin_dir, seg_dir, \
            final_call_dir, LAMBDA, MAX_WGD, chroms, rerun) for batch in batches]
        for future in futures:
            try:
                future.result()  # Wait for each batch to complete
            except Exception as e:
                print(f"Error in batch processing: {e}")
    # waiting for all threads to finish
    executor.shutdown(wait=True)
    # draw tracks
    for cell in cells:
        print(f'Drawing track for {cell}')
        draw(cell, final_call_dir, baf_alpha)
        plt.clf()
        
if __name__ == '__main__':
    main()