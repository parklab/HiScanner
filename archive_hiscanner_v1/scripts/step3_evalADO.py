import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from hmmlearn import hmm
import pickle
import argparse
import scipy


def get_vaf_by_chrom(chrom, vaf, depth_filter=5, aggregate=False, k=1):
    chrom = str(chrom)
    vaf['#CHROM'] = vaf['#CHROM'].astype(str)
    vaf = vaf[vaf['#CHROM'] == chrom]
    print('Number of hetSNPs in chromosome', chrom, ':', vaf.shape[0])
    if aggregate:
        vaf_agg = vaf[['A', 'B']].groupby(vaf.index // k).sum()
        vaf_agg['START'] = vaf['POS'].groupby(vaf.index // k).first()
        vaf_agg['END'] = vaf['POS'].groupby(vaf.index // k).last()
        vaf_agg['TOTAL'] = vaf_agg['A'] + vaf_agg['B']
        vaf_agg['pBAF'] = vaf_agg['B'] / vaf_agg['TOTAL']
        vaf_agg['BAF'] = vaf_agg['pBAF'].apply(lambda x: min(x, 1-x))
        vaf_agg.reset_index(drop=True)
        vaf = vaf_agg
    return vaf[vaf.TOTAL > depth_filter]

def main():
    parser = argparse.ArgumentParser(description="Process VAF data and train HMM model")
    parser.add_argument("--metadata_path", required=True, help="Path to metadata file")
    parser.add_argument("--hetsnp_dir", required=True, help="Directory containing hetSNP files")
    parser.add_argument("--ado_dir", required=True, help="Directory for ADO output")
    parser.add_argument("--plot", action="store_true", help="Generate plots")
    parser.add_argument("--depth_filter", type=int, default=5, help="Depth filter value")
    parser.add_argument("--rerun", action="store_true", help="Rerun all steps")
    parser.add_argument("--cell_name", required=True, help="Name of the cell(s) to perform ADO analysis")
    parser.add_argument("--ado_threshold", type=float, default=0.2, help="Threshold for hetSNPs in ADO")
    parser.add_argument("--aggregate", action="store_true", help="Aggregate every k hetSNPs")
    parser.add_argument("--k", type=int, default=1, help="Number of hetSNPs to aggregate")
    args = parser.parse_args()

    if args.aggregate:
        print(f'****Aggregating every {args.k} hetSNPs****')
        if args.k==1:
            print('****No aggregation performed****')
            args.aggregate = False

    os.makedirs(args.ado_dir, exist_ok=True)
    print(f'****Loading metadata from {args.metadata_path}****')
    metadata_all = pd.read_csv(args.metadata_path, sep='\t')
    metadata_all = metadata_all[metadata_all.singlecell == 'Y']

    if args.cell_name == 'all':
        cells = metadata_all.bamID.values
    else:
        cells = args.cell_name.split(',')
    
    n_cells = len(cells)
    print(f'****Performing ADO analysis for {n_cells} cells****')
    print('Cells:', cells)
        
        
    for cell in cells:
        cell_dir = args.ado_dir + f'/{cell}'
        os.makedirs(cell_dir, exist_ok=True)
        if args.plot:
            # check rerun flag
            if not os.path.exists(os.path.join(cell_dir, 'pBAF_distribution.png')) or args.rerun:
                print(f'****Plotting distribution of pBAF in hetSNPs for {cell}****')
                fig, ax = plt.subplots(dpi=200)
                fname = os.path.join(args.hetsnp_dir, f'{cell}.hetsnp.txt')
                vaf = pd.read_csv(fname, sep='\t')
                print('****Loading hetSNP data****')
                print('Total numer of hetSNPs:', vaf.shape[0])
                vaf['#CHROM'] = vaf['#CHROM'].astype(str)
                if args.aggregate:
                    vaf_agg = []
                    for chrom_train in range(1, 23):
                        selected = get_vaf_by_chrom(str(chrom_train), vaf, depth_filter=args.depth_filter, aggregate = args.aggregate, k=args.k)
                        vaf_agg.append(selected)
                    vaf = pd.concat(vaf_agg)
                else:
                    vaf = vaf[vaf.TOTAL > args.depth_filter]
                sns.kdeplot(vaf.pBAF, bw_adjust=0.5, fill=False, color='black', ax=ax, lw = 0.5)
                plt.savefig(os.path.join(cell_dir, 'pBAF_distribution.png'))
                plt.close()
                
        print('****Loading data****')
        if not os.path.exists(os.path.join(cell_dir, 'data.npy')) or args.rerun:
            data = []
            fname = os.path.join(args.hetsnp_dir, f'{cell}.hetsnp.txt')
            vaf = pd.read_csv(fname, sep='\t')
            vaf['#CHROM'] = vaf['#CHROM'].astype(str)
            for chrom_train in range(1, 23):
                selected = get_vaf_by_chrom(str(chrom_train), vaf, depth_filter=args.depth_filter, aggregate = args.aggregate, k=args.k)
                sequence = selected.BAF.values
                sequence_reshaped = sequence.reshape(-1, 1)
                data.append(sequence_reshaped)
            data = np.concatenate(data)
            np.save(os.path.join(cell_dir, 'data.npy'), data)
        else:
            print(f'****Loading data from {cell_dir}/data.npy ****')
            data = np.load(os.path.join(cell_dir, 'data.npy'))
        
        if not os.path.exists(os.path.join(cell_dir, 'hmm_model.pkl')) or args.rerun:
            print('****Training HMM model****')
            data_binary = data>args.ado_threshold
            data_binary = np.reshape(data_binary, (-1, 1))
            model = hmm.GaussianHMM(n_components=2, covariance_type="diag", verbose=False)
            model.fit(data_binary)
            mus = np.ravel(model.means_)
            sigmas = np.ravel(np.sqrt([np.diag(c) for c in model.covars_]))
            P = model.transmat_
            print(mus, sigmas, P)
            with open(os.path.join(cell_dir, 'hmm_model.pkl'), 'wb') as f:   
                pickle.dump(model, f)
            print('****Finished training HMM model****')
            print('Inferred parameters:')
            print('Means:', mus)
            print('Sigmas:', sigmas)
            print('Transition matrix:', P)
            
        else:
            model = pickle.load(open(os.path.join(cell_dir, 'hmm_model.pkl'), 'rb'))
            print('****Loaded HMM model****')
            print('Inferred parameters:')
            print('Means:', model.means_)
            print('Sigmas:', np.sqrt([np.diag(c) for c in model.covars_]))
            print('Transition matrix:', model.transmat_)
        
        del data
        problem_cells = []
        if not os.path.exists(os.path.join(cell_dir, 'result_gauss.txt')) or args.rerun:
            anno_df_combined = []
            fname = os.path.join(args.hetsnp_dir, f'{cell}.hetsnp.txt')
            vaf = pd.read_csv(fname, sep='\t')
            vaf['#CHROM'] = vaf['#CHROM'].astype(str)
            for chrom in range(1, 23):
                selected = get_vaf_by_chrom(str(chrom), vaf, depth_filter=args.depth_filter, aggregate = args.aggregate, k=args.k)
                if selected.shape[0] == 0:
                    print(f'Flagging {cell}: no hetSNPs in chromosome {chrom}')
                    problem_cells.append(cell)
                    continue
                sequence = selected.BAF.values 
                sequence_binary = sequence > args.ado_threshold
                sequence_reshaped = sequence_binary.reshape(-1, 1)
                evaluate, decoded = model.decode(sequence_reshaped)
                # determine whether the states are flipped
                corr = np.corrcoef(decoded.ravel(), sequence_reshaped.ravel())[0, 1]
                if corr > 0:
                    convert = {0: 1, 1: 0}
                else:
                    convert = {0: 0, 1: 1}
                sequence_converted = [convert[d] for d in decoded]
                seq_region_indices = scipy.ndimage.find_objects(scipy.ndimage.label(sequence_converted)[0])
                
                anno = []
                if args.aggregate:
                    START = np.array(selected.START)
                    END = np.array(selected.END)
                    for region in seq_region_indices:
                        n = decoded[region].shape[0]
                        start_idx, end_idx = region[0].start, region[0].stop
                        if start_idx == end_idx-1:
                            continue
                        start_pos, end_pos = START[start_idx], END[end_idx-1]
                        state = decoded[start_idx]
                        length = end_pos - start_pos
                        baf = sequence[start_idx:end_idx].mean()
                        anno.append([start_pos, end_pos, n, length, state, baf])
                else:
                    pos = np.array(selected.POS)
                    for region in seq_region_indices:
                        n = decoded[region].shape[0]
                        start_idx, end_idx = region[0].start, region[0].stop
                        if start_idx == end_idx-1:
                            continue
                        start_pos, end_pos = pos[start_idx], pos[end_idx-1]
                        state = decoded[start_idx]
                        length = end_pos - start_pos
                        baf = sequence[start_idx:end_idx].mean()
                        anno.append([start_pos, end_pos, n, length, state, baf])


                anno_df = pd.DataFrame(anno, columns=['start', 'end', 'n', 'length', 'state', 'mean_baf'])
                anno_df['chrom'] = chrom
                anno_df_combined.append(anno_df)
            anno_df_combined = pd.concat(anno_df_combined)
            anno_df_combined.to_csv(os.path.join(cell_dir, 'result_gauss.txt'), sep='\t', index=None)
        else:
            anno_df_combined = pd.read_csv(os.path.join(cell_dir, 'result_gauss.txt'), sep='\t')
        anno_df_combined = anno_df_combined[anno_df_combined.length > 1]
        
        if cell in problem_cells:
            bin_size_min_round_kb = np.nan
        else:
            anno_df_combined['log10_len'] = np.log10(anno_df_combined['length'].astype(float))
            bin_size_min = anno_df_combined.length.quantile(q=.99)
            bin_size_min_round_kb = np.ceil(bin_size_min/1000)

        print('Suggested minimum bin size (kb): ', bin_size_min_round_kb)
        # save the suggested bin size
        if not os.path.exists(os.path.join(cell_dir, 'bin_size_min.txt')) or args.rerun:    
            with open(os.path.join(cell_dir, 'bin_size_min.txt'), 'w') as f:
                f.write(str(bin_size_min_round_kb) + ' kb')
        else:
            with open(os.path.join(cell_dir, 'bin_size_min.txt'), 'r') as f:
                bin_size_min_round_kb = float(f.read().split(' ')[0])
                print('Loaded suggested bin size:', bin_size_min_round_kb)

        if args.plot:
            if not os.path.exists(os.path.join(cell_dir, 'ado_length_distribution_per_chrom.png')) or args.rerun:
                if np.isnan(bin_size_min_round_kb):
                    print(f'Skip plotting ADO length distribution for {cell} due to reasons above')
                else:
                    print(f'****Plotting ADO length distribution for {cell}****')
                    fig, ax = plt.subplots(dpi=500, figsize=(4, 3))
                    for chrom in range(1, 23):
                        sns.kdeplot(data=anno_df_combined[anno_df_combined.chrom == chrom], x='log10_len',
                                    ax=ax, color='black', lw=.5)
                    plt.axvline(np.log10(bin_size_min_round_kb*1000), label = 'Suggested bin size', color='red', lw=1)
                    plt.xlim(0, 9)
                    plt.legend()
                    plt.xlabel('Allelic dropout length (log10 scale)')
                    plt.savefig(os.path.join(cell_dir, 'ado_length_distribution_per_chrom.png'))
                    plt.close()
                
    # gather results
    if not os.path.exists(os.path.join(args.ado_dir, 'min_bin_size_per_cell.txt')) or args.rerun:
        min_bin_size = []
        for cell in cells:
            cell_dir = args.ado_dir + f'/{cell}'
            with open(os.path.join(cell_dir, 'bin_size_min.txt'), 'r') as f:
                bin_size_min_round_kb = float(f.read().split(' ')[0])
            min_bin_size.append([bin_size_min_round_kb, cell])
        min_bin_size = pd.DataFrame(min_bin_size, columns=['bin_size_min_kb', 'cell'])
        min_bin_size.to_csv(os.path.join(args.ado_dir, 'min_bin_size_per_cell.txt'), sep='\t', index=None)
        suggested_bin_size = min_bin_size.dropna().bin_size_min_kb.mean()
    else:
        min_bin_size = pd.read_csv(os.path.join(args.ado_dir, 'min_bin_size_per_cell.txt'), sep='\t')
        suggested_bin_size = min_bin_size.bin_size_min_kb.dropna().mean()
    print(f'Suggested bin size for all cells: {suggested_bin_size} kb')
    # boxplot of suggested bin size
    fig, ax = plt.subplots(dpi=200)
    ax = sns.boxplot(data=min_bin_size, x='bin_size_min_kb', ax=ax)
    title = f'Suggested bin size for ADO detection: {np.round(suggested_bin_size)} kb | n = {n_cells}'
    ax.set_title(title, fontsize=10)
    plt.savefig(os.path.join(args.ado_dir, 'min_bin_size_per_cell.png'))
    plt.close()
    print('****Finished ADO analysis****')
if __name__ == "__main__":
    main()