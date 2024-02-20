
import os
import json 
import pandas as pd
from .utils import * 

def segment(json_file_path, chroms=range(1,23)):
    with open(json_file_path, 'r') as file:
        args = json.load(file)
    # check if bicseq is in the path
    bicseq = args.get('bicseq_path')
    singlecell = args.get('singlecell').split(',')
    stem = args.get('stem')
    RDR_LAMBDA = args.get('LAMBDA')
    assert os.system('scanner-segment -h > /dev/null 2>&1') == 0

    print('Reshaping data for bicseq segmentation')
    # reshape the data for bicseq segmentation
    for c in singlecell:
        cell_df = pd.read_csv(f'{stem}/{c}.data.txt',sep='\t')
        cell_df = cell_df[['CHROM','START','END','OBS','EXP']]
        cell_df.columns = ['CHROM','START','END',f'A_{c}',f'B_{c}']
        
    for chrom, chrom_df in cell_df.groupby('CHROM'):
        chrom_df = chrom_df.iloc[:,1:]
        chrom_df.index = range(chrom_df.shape[0])
        chrom_df.to_csv(f'{stem}/bicseq_rdr_{chrom}',sep='\t',index=None)


    print('Running bicseq segmentation')
    for c in singlecell:
        for chrom in chroms:
            input_file = f'{stem}/bicseq_rdr_{chrom}'
            cmd = f'scanner-segment -i {input_file} -l {RDR_LAMBDA} -o {input_file}_seg'
            print(cmd)
            os.system(cmd)

    seg_all_rdr = []
    for chrom in chroms:

        seg_rdr = pd.read_csv(f'{stem}/bicseq_rdr_{chrom}_seg', sep='\t')
        seg_rdr['chrom'] = f'{chrom}'
        seg_rdr = seg_rdr[['chrom', 'start', 'end']]
        seg_all_rdr.append(seg_rdr)

    seg_all_rdr = pd.concat(seg_all_rdr)

    seg_all_rdr.to_csv(f'{stem}/rdr_seg.tsv', sep='\t', index=None)

    print(f'Segmentation complete. Segments saved to {stem}/rdr_seg.tsv')
    print('Total number of RDR segments:', seg_all_rdr.shape[0])
    print('Deleting intermediate files...')
    # delete the intermediate files
    os.system(f'rm {stem}/bicseq*')
    return


def infer_copy_number(json_file_path):
    """
    Infers copy number from the given JSON file.

    Args:
        json_file_path (str): The path to the JSON file.

    Returns:
        None
    """
    with open(json_file_path, 'r') as file:
        args = json.load(file)
    stem = args.get('stem')
    MAX_WGD = args.get('MAX_WGD')
    singlecell = args.get('singlecell').split(',')
    for cell in singlecell:
        seg_path = f'{stem}/rdr_seg.tsv'
        data_path = f'{stem}/{cell}.data.txt'
    print('Processing ', stem)
    seg_all, seg_all_merged, cell_data = postprocess(seg_path, data_path, MAX_WGD=1)
    seg_all.to_csv(f'{stem}/{cell}.short.wgd{MAX_WGD}.txt', sep='\t', index=False)
    seg_all_merged.to_csv(f'{stem}/{cell}.short_merged.wgd{MAX_WGD}.txt', sep='\t', index=False)
    cell_data.to_csv(f'{stem}/{cell}.long.wgd{MAX_WGD}.txt', sep='\t', index=False)
    print('saved to :', f'{stem}/{cell}.short.wgd{MAX_WGD}.txt')

    return

