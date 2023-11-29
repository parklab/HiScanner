import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from joblib import Parallel, delayed
import json
import joblib
import contextlib
from .utils import *

eps=1e-10

def call_bcftools_query(cell:str,
                        vcf:str,
                        outfilename:str):
    '''
    cell: cell name
    vcf: path to vcf file
    outfilename: path to output file
    '''
    print(f'Querying {cell}')
    cmd = 'bcftools query '+ \
            vcf + \
            f' -s {cell}' + \
            " -f '%CHROM\\t%POS\\t%REF\\t%ALT{0}[\\t%AD{0}\\t%AD{1}]\n' > " \
            + outfilename
    os.system(cmd)
    return

def baf_preprocess(cell:str,
                   fn:str, 
                   phase_path:str, 
                   bin_path:str, 
                   outfilename:str, 
                   chroms:list=range(1,23),
                   overwrite:bool=False):
    
    # check if outfilename exists
    if check_file(outfilename + '.txt', overwrite):
        return
    phase_df = pd.read_csv(phase_path,sep='\t')
    skip_snp_raw =  check_file(outfilename + '.snp_raw.txt', overwrite)
    if not skip_snp_raw:  
        df = pd.read_csv(fn,sep='\t',header=None)
        # filter positions
        df['UNIQUE_POS'] = ['{}_{}'.format(c,p) for c,p in df[[0,1]].values]
        df = df[df.UNIQUE_POS.isin(phase_df.UNIQUE_POS.values)]
        # filter genotypes
        df = df[df[4]!='.']
        df = df[df[5]!='.']
        df.loc[:,'A'] = df[4].values.astype(int)
        df.loc[:,'B'] = df[5].values.astype(int)
        df.loc[:,'B_copy'] = df[5].values.astype(int)
        reverse_pos = phase_df[phase_df.phasedgt == '1|0'].UNIQUE_POS.values
        df.loc[df.UNIQUE_POS.isin(reverse_pos),'B'] = df.loc[df.UNIQUE_POS.isin(reverse_pos)]['A'].values
        df.loc[df.UNIQUE_POS.isin(reverse_pos),'A'] = df.loc[df.UNIQUE_POS.isin(reverse_pos)]['B_copy'].values
        df = df[[0,1,'A','B']]
        df.columns = ['CHROM','POS','A','B']
        df.loc[:,'TOTAL'] = df['A'] + df['B']
        # filter positions with no reads
        df = df[df.TOTAL!=0]
        df.loc[:,'pBAF'] = df.B / df.TOTAL
        df.loc[:,'BAF'] = [min(a,b)/t for a,b,t in df[['A','B','TOTAL']].values]
        # save unbinned snp data
        print('Saving unbinned snp data to '+outfilename+'.snp_raw.txt')
        df.to_csv(outfilename+'.snp_raw.txt',sep='\t',index=None)
    else:
        df = pd.read_csv(outfilename+'.snp_raw.txt',sep='\t')

    # get rdr and grids 
    rdr, obs, exp = [],[],[]
    grid = {}
    
    for i in chroms:
        # check if the binned rdr data exists
        assert os.path.exists(f'{bin_path}/{i}.bin'), f'{bin_path} does not exist'
        b = pd.read_csv(f'{bin_path}/{i}.bin',sep='\t')
        rdr.extend(b['obs'].values/b['expected'].values)
        obs.extend(b['obs'].values)
        exp.extend(b['expected'].values)
        grid[i] = b[['start','end']].values
    
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
    df_aux.loc[:,'pBAF'] = df_aux['B'].values / (eps+df_aux['A'].values+ df_aux['B'].values)
    df_aux.loc[:,'BAF'] = [min(1-v,v) for v in df_aux.pBAF.values]
    df_aux.loc[:,'MAJOR'] = [max(a,b) for a,b, in df_aux[['A','B']].values]
    df_aux.loc[:,'MINOR'] = [min(a,b) for a,b, in df_aux[['A','B']].values]
    
    # annotate rdr
    df_aux['RDR'] = rdr
    df_aux['OBS'] = obs
    df_aux['EXP'] = exp
    # save preprocessed data
    print('Saving preprocessed data to '+outfilename+'.txt')
    df_aux.to_csv(outfilename + '.txt', sep='\t',index=None)

    return


def preprocess(json_file_path, overwrite=False):
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
    jobs = args.get('j')  # Providing a default value if 'j' is not in the JSON
    singlecell = args.get('singlecell')
    singlecell = singlecell.split(',')
    # Check paths in JSON file
    assert os.path.exists(bin_path), f'bin_path {bin_path} does not exist'
    assert os.path.exists(phase_file), f'phase_file {phase_file} does not exist'
    
    assert os.path.exists(gatk_vcf), f'gatk_vcf {gatk_vcf} does not exist'

    # Create the stem directory if it doesn't exist
    if not os.path.exists(stem):
        os.mkdir(stem)

    # Check if bcftools is installed
    assert os.system('bcftools --version') == 0, 'bcftools is not installed or in PATH'

    print('Getting germline snps')

    # check if phased snps have been extracted
    skip_germline = check_file(f'{stem}/phased.pos', overwrite)
    if not skip_germline:
        os.system(f'bcftools view -h {gatk_vcf} > {stem}/header.txt')
        f = open(f'{stem}/header.txt','r')
        colnames = ['CHROM','POS','REF','ALT']
        cellnames = f.readlines()[-1][:-1].split('\t')[9:]
        f.close()
        cmd = 'bcftools query '+ \
            gatk_vcf +\
            ' -s '+ germline+\
            ' -f '%CHROM\\t%POS\\t%REF\\t%ALT{0}[\\t%AD{0}\\t%AD{1}]\n' > '+\
            f'{stem}/germline.gatk.txt'
        os.system(cmd)
        print('Germline snps saved to '+f'{stem}/germline.gatk.txt')
        germline_df = pd.read_csv(f'{stem}/germline.gatk.txt', sep='\t', header=None)
        germline_df = germline_df[(germline_df[4]!='.') & (germline_df[5]!='.')]
        germline_df.loc[:,6] = [ (min(int(a),int(b))) /(int(a)+int(b)+eps) 
                                for a,b in germline_df[[4,5]].values ]
        germline_df['UNIQUE_POS'] = ['{}_{}'.format(c,p) for c,p in germline_df[[0,1]].values]
        # filtering out homozygous snps
        germline_df = germline_df[germline_df[6]>.1] 
        phase = pd.read_csv(phase_file,comment='#',sep='\t',header=None)
        phase = phase.loc[:,[0,1,3,4,9]]
        phase.columns = ['CHROM','POS','REF','ALT','phasedgt']
        phase['UNIQUE_POS'] = ['{}_{}'.format(c,p) for c,p in phase[['CHROM','POS']].values]
        phase = phase[phase.UNIQUE_POS.isin(germline_df.UNIQUE_POS.values)]
        phase.index = range(phase.shape[0]) 
        print('Saving phased het germline snps to '+f'{stem}/phased.pos')
        phase.to_csv(f'{stem}/phased.pos',sep='\t',index=None)

    f = open(f'{stem}/header.txt','r')
    colnames = ['CHROM','POS','REF','ALT']
    cellnames = f.readlines()[-1][:-1].split('\t')[9:]
    f.close()

    if len(singlecell)>0:
        cellnames = [c for c in cellnames if c in singlecell]
        
    print('Start querying individual cells')
    with tqdm_joblib(tqdm(desc='Cells processed', total=len(cellnames))) as progress_bar:
        Parallel( n_jobs = jobs )( delayed( call_bcftools_query )(cell,gatk_vcf,f'{stem}/{cell}.gatk.snp.txt') 
                                  for cell in cellnames)
    print('Start preprocessing BAF')
    with tqdm_joblib(tqdm(desc='Cells processed', total=len(cellnames))) as progress_bar:
        Parallel( n_jobs = jobs )( delayed( baf_preprocess )(cell,
                                                    f'{stem}/{cell}.gatk.snp.txt',
                                                    f'{stem}/phased.pos',
                                                    f'{bin_path}/{cell}/',
                                                    f'{stem}/{cell}.data',
                                                    range(1,23),
                                                    overwrite)
                                                        for cell in cellnames if os.path.exists(f'{bin_path}/{cell}/'))



