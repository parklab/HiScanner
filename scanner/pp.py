import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from joblib import Parallel, delayed
import json
import contextlib
import subprocess
from .utils import *

eps = 1e-10

def call_bcftools_query(cell: str, vcf: str, out_fn: str):
    '''
    Query the VCF file using bcftools for a specific cell.

    Args:
        cell (str): The name of the cell.
        vcf (str): The path to the VCF file.
        out_fn (str): The path to the output file.

    Returns:
        None
    '''
    print(f'Querying {cell}')
    cmd = f'bcftools query {vcf} -s {cell} -f' + \
            ' "%CHROM\\t%POS\\t%REF\\t%ALT{0}[\\t%AD{0}\\t%AD{1}]\n" ' + \
            f'> {out_fn}'
    print("Running " + cmd)
    os.system(cmd)

def baf_preprocess(cell: str, fn: str, phase_path: str, bin_path: str, out_fn: str, chroms: list = range(1, 23), overwrite: bool = False):
    '''
    Preprocesses BAF (B-Allele Frequency) data.

    Args:
        cell (str): The name of the cell.
        fn (str): The file path of the input data.
        phase_path (str): The file path of the phase data.
        bin_path (str): The directory path of the binned data.
        out_fn (str): The output file name.
        chroms (list, optional): The list of chromosomes to process. Defaults to human autosomes (1-22).
        overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.

    Returns:
        None
    '''
    print(f'Preprocessing {cell}')
    # check if out_fn exists
    if check_file(f'{out_fn}.txt', overwrite):
        return
    
    phase_df = pd.read_csv(phase_path, sep='\t')
    skip_snp_raw = check_file(f'{out_fn}.snp_raw.txt', overwrite)
    if not skip_snp_raw:
        df = pd.read_csv(fn, sep='\t', header=None)
        # concatenate chrom and pos to create a unique identifier
        df['UNIQUE_POS'] = df[0].astype(str) + '_' + df[1].astype(str)
        # filter out non-phased snps
        df = df[df['UNIQUE_POS'].isin(phase_df['UNIQUE_POS'].values)]
        df = df[(df[4] != '.') & (df[5] != '.')]
        df[['A', 'B']] = df[[4, 5]].astype(int)
        df['B_copy'] = df['B'].values.copy()
        # reverse A and B if the snp is phased as 1|0
        reverse_pos = phase_df[phase_df['phasedgt'] == '1|0']['UNIQUE_POS'].values
        df.loc[df['UNIQUE_POS'].isin(reverse_pos), ['A', 'B']] = df.loc[df['UNIQUE_POS'].isin(reverse_pos), ['B_copy', 'A']].values
        df = df[[0, 1, 'A', 'B']].rename(columns={0: 'CHROM', 1: 'POS'})
        df['TOTAL'] = df['A'] + df['B']
        df = df[df['TOTAL'] != 0]
        df['pBAF'] = df['B'] / df['TOTAL']
        df['BAF'] = df[['A', 'B', 'TOTAL']].min(axis=1) / df['TOTAL']
        print(f'Saving unbinned snp data to {out_fn}.snp_raw.txt')
        df.to_csv(f'{out_fn}.snp_raw.txt', sep='\t', index=None)
    else:
        df = pd.read_csv(f'{out_fn}.snp_raw.txt', sep='\t')

    rdr, obs, exp = [], [], []
    grid = {}
    for i in chroms:
        assert os.path.exists(f'{bin_path}/{i}.bin'), f'{bin_path} does not exist'
        b = pd.read_csv(f'{bin_path}/{i}.bin', sep='\t')
        rdr.extend(b['obs'] / b['expected'])
        obs.extend(b['obs'])
        exp.extend(b['expected'])
        grid[i] = b[['start', 'end']].values

    # aggregate into bins, annotate BAF
    df_aux = []
    for chrom, g in grid.items():
        for i in range(len(g)):
            start, end = g[i]
            selected = df[(df['CHROM']==chrom) & (df['POS']>=start) & (df['POS']<=end)]
            df_aux.append([chrom, start, end, selected.A.sum(), selected.B.sum(), selected.shape[0]])
    df_aux = pd.DataFrame(df_aux)
    df_aux.columns = ['CHROM','START','END','A','B','N']
    df_aux.loc[:,'TOTAL'] = df_aux['A'].values+ df_aux['B'].values
    df_aux.loc[:, 'pBAF'] = df_aux['B'] / (eps + df_aux['TOTAL'])
    df_aux.loc[:, 'MAJOR'] = df_aux[['A', 'B']].max(axis=1)
    df_aux.loc[:, 'MINOR'] = df_aux[['A', 'B']].min(axis=1)
    df_aux.loc[:, 'BAF'] = df_aux['MINOR'] / df_aux['TOTAL']
    df_aux.loc[:, 'RDR'], df_aux['OBS'], df_aux['EXP'] = rdr, obs, exp
    print(f'Saving preprocessed data to {out_fn}.txt')
    df_aux.to_csv(f'{out_fn}.txt', sep='\t', index=None)
    
    return 

def preprocess(json_file_path, overwrite=False):
    """
    Main preprocessing function.

    Args:
        json_file_path (str): The path to the JSON file containing the preprocessing arguments.
        overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.
    """
    
    # Read data from the JSON file
    assert os.path.exists(json_file_path), f'JSON file {json_file_path} does not exist'
    with open(json_file_path, 'r') as file:
        args = json.load(file)
    
    # Accessing arguments from the JSON data
    bin_path = args.get('bin_path')
    phase_file = args.get('phase_file')
    germline = args.get('germline')
    gatk_vcf = args.get('gatk_vcf')
    stem = args.get('stem')
    jobs = int(args.get('j', 1))  # Providing a default value if 'j' is not in the JSON
    singlecell = args.get('singlecell', '').split(',')

    # Check paths in JSON file
    assert os.path.exists(bin_path), f'bin_path {bin_path} does not exist'
    assert os.path.exists(phase_file), f'phase_file {phase_file} does not exist'
    assert os.path.exists(gatk_vcf), f'gatk_vcf {gatk_vcf} does not exist'

    # Create the stem directory if it doesn't exist
    os.makedirs(stem, exist_ok=True)

    # Check if bcftools is installed
    assert os.system('bcftools --version') == 0, 'bcftools is not installed or in PATH'

    print('Getting germline snps')

    # check if phased snps have been extracted
    skip_germline = check_file(f'{stem}/phased.pos', overwrite)
    if not skip_germline:
        os.system(f'bcftools view -h {gatk_vcf} > {stem}/header.txt')
        colnames = ['CHROM', 'POS', 'REF', 'ALT']
        cmd = f'bcftools query {gatk_vcf} -s {germline} -f '\
              + '"%CHROM\\t%POS\\t%REF\\t%ALT{0}[\\t%AD{0}\\t%AD{1}]\n"' \
              + f'> {stem}/germline.gatk.txt'
        print("Running " + cmd)
        os.system(cmd)
        print(f'Germline snps saved to {stem}/germline.gatk.txt')
        germline_df = pd.read_csv(f'{stem}/germline.gatk.txt', sep='\t', header=None)
        germline_df = germline_df[(germline_df[4] != '.') & (germline_df[5] != '.')]
        germline_df[[4, 5]] = germline_df[[4, 5]].astype(int)
        germline_df[6] = germline_df[[4, 5]].min(axis=1) / (germline_df[[4, 5]].sum(axis=1) + eps)
        germline_df['UNIQUE_POS'] = germline_df[0].astype(str) + '_' + germline_df[1].astype(str)
        germline_df = germline_df[germline_df[6] > 0.1]
        phase = pd.read_csv(phase_file, comment='#', sep='\t', header=None, usecols=[0, 1, 3, 4, 9])
        phase.columns = ['CHROM', 'POS', 'REF', 'ALT', 'phasedgt']
        phase['UNIQUE_POS'] = phase['CHROM'].astype(str) + '_' + phase['POS'].astype(str)
        phase = phase[phase['UNIQUE_POS'].isin(germline_df['UNIQUE_POS'].values)]
        phase.reset_index(drop=True, inplace=True)
        print(f'Saving phased het germline snps to {stem}/phased.pos')
        phase.to_csv(f'{stem}/phased.pos', sep='\t', index=None)

    f = open(f'{stem}/header.txt','r')
    colnames = ['CHROM', 'POS', 'REF', 'ALT']
    cellnames = f.readlines()[-1][:-1].split('\t')[9:]
    f.close()

    if len(singlecell) > 0:
        cellnames = [c for c in cellnames if c in singlecell]
        
    print('Start querying individual cells')
    with tqdm_joblib(tqdm(desc='Cells processed', total=len(cellnames))) as progress_bar:
        Parallel(n_jobs=jobs)(delayed(call_bcftools_query)(cell, gatk_vcf, f'{stem}/{cell}.gatk.snp.txt') for cell in cellnames)
    print('Start preprocessing BAF')
    with tqdm_joblib(tqdm(desc='Cells processed', total=len(cellnames))) as progress_bar:
        Parallel(n_jobs=jobs)(delayed(baf_preprocess)(cell, f'{stem}/{cell}.gatk.snp.txt', f'{stem}/phased.pos', f'{bin_path}/{cell}/', f'{stem}/{cell}.data', range(1, 23), overwrite) for cell in cellnames if os.path.exists(f'{bin_path}/{cell}/'))

def run_rdr_norm():
    subprocess.run(["perl", "scripts/mbicseq-norm/BICseq2-norm.pl"])
