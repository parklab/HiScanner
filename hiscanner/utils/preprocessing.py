import os
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Union
from ..logger import logger

def prep_input_table(cell: str, 
                    hetsnp_dir: Union[str, Path], 
                    bin_dir: Union[str, Path], 
                    chroms: List[str] = None) -> pd.DataFrame:
    """
    Prepare input table for CNV analysis by combining heterozygous SNP and bin data.
    
    Parameters
    ----------
    cell : str
        Cell identifier
    hetsnp_dir : str or Path
        Directory containing heterozygous SNP data
    bin_dir : str or Path
        Directory containing bin data
    chroms : list of str, optional
        List of chromosomes to process. If None, processes chromosomes 1-22
        
    Returns
    -------
    pd.DataFrame
        Combined data table with CNV analysis input data
    """
    if chroms is None:
        chroms = [str(i) for i in range(1, 23)]

    # Read phased BAF
    fname = Path(hetsnp_dir) / f'{cell}.hetsnp.txt'
    vaf = pd.read_csv(fname, sep='\t', low_memory=False)
    vaf = vaf.rename(columns={'#CHROM':'CHROM'})
    vaf['CHROM'] = vaf['CHROM'].astype(str)

    # Read RDR and grids
    rdr, obs, exp = [], [], []
    grid = {}
    
    for chrom in chroms:
        b = pd.read_csv(f'{bin_dir}/{cell}/{chrom}.bin', sep='\t')
        rdr.extend(b['obs'].values/b['expected'].values)
        obs.extend(b['obs'].values)
        exp.extend(b['expected'].values)
        grid[chrom] = b[['start','end']].values

    # Aggregate into bins and annotate BAF
    cell_table = []
    for chrom, g in grid.items():
        for i in range(len(g)):
            start, end = g[i]
            selected = vaf[(vaf['CHROM']==chrom) & 
                (vaf['POS']>=start) & 
                (vaf['POS']<=end)]
            cell_table.append([
                chrom, start, end, 
                selected.A.sum(), 
                selected.B.sum(), 
                selected.shape[0]
            ])
    
    # Create DataFrame and add derived columns
    cell_table = pd.DataFrame(cell_table)
    cell_table.columns = ['CHROM','START','END','A','B','N']
    cell_table.loc[:,'TOTAL'] = cell_table['A'] + cell_table['B']
    cell_table.loc[:,'pBAF'] = cell_table['B'] / (cell_table['TOTAL'] + 1e-10)
    cell_table.loc[:,'BAF'] = cell_table['pBAF'].apply(lambda x: min(x, 1-x))
    cell_table.loc[:,'MAJOR'] = cell_table[['A','B']].max(axis=1)
    cell_table.loc[:,'MINOR'] = cell_table[['A','B']].min(axis=1)
    
    # Add RDR information
    cell_table['RDR'] = rdr
    cell_table['OBS'] = obs
    cell_table['EXP'] = exp
    cell_table['CHROM'] = cell_table['CHROM'].astype(str)
    
    return cell_table

def prep_input_table_rdr_only(cell: str, 
                            bin_dir: Union[str, Path], 
                            chroms: List[str] = None) -> pd.DataFrame:
    """
    Prepare input table for CNV analysis using bin data only.
    
    Parameters
    ----------
    cell : str
        Cell identifier
    bin_dir : str or Path
        Directory containing bin data
    chroms : list of str, optional
        List of chromosomes to process. If None, processes chromosomes 1-22
        
    Returns
    -------
    pd.DataFrame
        Combined data table with read depth information
    """
    if chroms is None:
        chroms = [str(i) for i in range(1, 23)]
        
    # Initialize data collection
    rdr, obs, exp = [], [], []
    positions = []
    chromosomes = []
    
    # Read and combine bin data
    for chrom in chroms:
        bin_path = Path(bin_dir) / cell / f'{chrom}.bin'
        if not bin_path.exists():
            logger.warning(f"No bin file found for chromosome {chrom}")
            continue
            
        bin_data = pd.read_csv(bin_path, sep='\t')
        rdr.extend(bin_data['obs'].values / bin_data['expected'].values)
        obs.extend(bin_data['obs'].values)
        exp.extend(bin_data['expected'].values)
        positions.extend(zip(bin_data['start'], bin_data['end']))
        chromosomes.extend([chrom] * len(bin_data))
    
    # Create DataFrame with consistent column names
    cell_table = pd.DataFrame({
        'CHROM': chromosomes,
        'START': [p[0] for p in positions],
        'END': [p[1] for p in positions],
        'RDR': rdr,
        'OBS': obs,
        'EXP': exp
    })
    
    cell_table['CHROM'] = cell_table['CHROM'].astype(str)
    return cell_table
