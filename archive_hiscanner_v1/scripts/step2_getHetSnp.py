import os
import pandas as pd
import subprocess
from math import ceil
import sys
import argparse

def process_batch(batch, gatk_vcf, output_dir, phase, rerun=False):
    processes = []
    for cell in batch:
        outfilename = os.path.join(output_dir, f'{cell}')
        if not os.path.exists(outfilename + '.txt') or rerun:
            print(f'****Extracting phased hets for {cell}****')
            cmd = [
                'bcftools', 'query',
                gatk_vcf,
                '-s', cell,
                '-f', '%CHROM\\t%POS\\t%REF\\t%ALT{0}[\\t%AD{0}\\t%AD{1}]\n'
            ]
            with open(outfilename + '.txt', 'w') as outfile:
                process = subprocess.Popen(cmd, stdout=outfile)
                processes.append((process, cell, outfilename))
        else:
            print(f'{outfilename}.txt exists, skipping')
            process_hetsnp(outfilename, phase, rerun)

    for process, cell, outfilename in processes:
        process.wait()
        print(f'Finished extracting phased hets for {cell}')
        process_hetsnp(outfilename, phase, rerun)

def process_hetsnp(outfilename, phase, rerun):
    hetsnp_filename = f'{outfilename}.hetsnp.txt'
    if not os.path.exists(hetsnp_filename) or rerun:
        df = pd.read_csv(outfilename + '.txt', sep='\t', header=None)
        # filter positions
        df['UNIQUE_POS'] = ['{}_{}'.format(c,p) for c,p in df[[0,1]].values]
        df = df[df.UNIQUE_POS.isin(phase.UNIQUE_POS.values)]
        df = df[df[4]!='.']
        df = df[df[5]!='.']
        df.loc[:,'A'] = df[4].values.astype(int)
        df.loc[:,'B'] = df[5].values.astype(int)
        df.loc[:,'B_copy'] = df[5].values.astype(int)
        reverse_pos = phase[phase.phasedgt == '1|0'].UNIQUE_POS.values
        df.loc[df.UNIQUE_POS.isin(reverse_pos),'B'] = df.loc[df.UNIQUE_POS.isin(reverse_pos)]['A'].values
        df.loc[df.UNIQUE_POS.isin(reverse_pos),'A'] = df.loc[df.UNIQUE_POS.isin(reverse_pos)]['B_copy'].values
        df = df[[0,1,'A','B']]
        df.columns = ['#CHROM','POS','A','B']
        df.loc[:,'TOTAL'] = df['A'] + df['B']
        df = df[df.TOTAL!=0]
        df.loc[:,'pBAF'] = df.B / df.TOTAL
        df.loc[:,'BAF'] = [min(a,b)/t for a,b,t in df[['A','B','TOTAL']].values]
        df.to_csv(hetsnp_filename, sep='\t', index=None)
    else:
        print(f'{hetsnp_filename} exists, skipping')

def main(args):
    print(f'****Python version****\n{sys.version}')
    # check bcftools is loaded
    sys_msg = os.system('bcftools -v 2> /dev/null')
    if sys_msg != 0:
        raise Exception('bcftools is not loaded')
    print('****Checking input files****')
    if not os.path.exists(args.gatk_vcf):
        raise Exception('gatk_vcf not found')
    elif not args.gatk_vcf.endswith('.gz'):
        raise Exception('gatk_vcf should be bgzipped')
    else:
        print(f'gatk_vcf: {args.gatk_vcf}')
    if not os.path.exists(args.gatk_vcf + '.tbi'):
        raise Exception('gatk_vcf.tbi not found')
    
    if not os.path.exists(args.phase_file):
        raise Exception('phase_file not found')
    os.chdir(args.out_dir)
    print(f'****Working directory****\n{os.getcwd()}')

    # read metadata
    meta = pd.read_csv(args.metadata, sep='\t')

    if not os.path.exists('header.txt') or args.rerun:
        os.system(f'bcftools view -h {args.gatk_vcf} > header.txt')
    else:
        print('header.txt exists, skipping')
        
    print('****Getting germline snps***')
    f = open(f'header.txt','r')
    colnames = ['#CHROM','POS','REF','ALT']
    cellnames = f.readlines()[-1][:-1].split('\t')[9:]
    f.close()
    print(f'Number of samples: {len(cellnames)}')
    # check if the number of samples is matching with metadata
    if len(cellnames) != len(meta):
        print('Warning: number of samples in metadata does not match with the vcf file.')
    else:
        print('Number of samples in metadata matches with the vcf file.')
        
    germline = meta[meta['singlecell']=='N']['bamID'].tolist()
    # check if there's only one germline sample
    if len(germline) != 1:
        raise Exception('There should be only one bulk sample (either too many or none). Check your metadata.')
    germline = germline[0]
    print(f'Bulk sample name: {germline}')

    os.makedirs('phased_hets', exist_ok=True)
    cmd = 'bcftools query '+ \
        args.gatk_vcf +\
        ' -s '+germline+\
        " -f '%CHROM\\t%POS\\t%REF\\t%ALT{0}[\\t%AD{0}\\t%AD{1}]\n' > "+\
        f'phased_hets/germline.txt'
        
    if not os.path.exists('phased_hets/germline.txt') or args.rerun:
        os.system(cmd)
        print(f'Finished running {cmd}')
    else:
        print('germline.txt exists, skipping')

    print('****Germline snps extracted****')

    if not os.path.exists('phased_hets/phased.pos') or args.rerun: 
        germline_df = pd.read_csv(f'phased_hets/germline.txt', sep='\t', header=None)
        germline_df = germline_df[(germline_df[4]!='.') & (germline_df[5]!='.')].dropna()
        germline_df.loc[:,6] = [ (min(int(a),int(b))) /(int(a)+int(b)+1e-10)
                                for a,b in germline_df[[4,5]].values ]
        germline_df['UNIQUE_POS'] = ['{}_{}'.format(c,p) for c,p in germline_df[[0,1]].values]
        print('filtering out homozygous snps')
        germline_df = germline_df[germline_df[6]>.45] 
        phase = pd.read_csv(args.phase_file,comment='#',sep='\t',header=None)
        phase = phase.loc[:,[0,1,3,4,9]]
        phase.columns = ['#CHROM','POS','REF','ALT','phasedgt']
        phase['UNIQUE_POS'] = ['{}_{}'.format(c,p) for c,p in phase[['#CHROM','POS']].values]
        phase = phase[phase.UNIQUE_POS.isin(germline_df.UNIQUE_POS.values)]
        phase.index = range(phase.shape[0])
        print('getting phased het snps')
        phase.to_csv(f'phased_hets/phased.pos',sep='\t',index=None)
    else:
        print('phased.pos exists, skipping')
        phase = pd.read_csv(f'phased_hets/phased.pos', sep='\t')
    # For each cell, call bcftools to extract the phased hets
    singlecells = meta[meta['singlecell']=='Y']['bamID'].tolist()

    # Process cells in batches
    num_batches = ceil(len(singlecells) / args.batch_size)
    for i in range(num_batches):
        start = i * args.batch_size
        end = min((i + 1) * args.batch_size, len(singlecells))
        batch = singlecells[start:end]
        
        print(f'****Processing batch {i+1}/{num_batches}****')
        process_batch(batch, args.gatk_vcf, 'phased_hets/', phase, args.rerun)
    print('****Finished processing all cells****')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF files for phased heterozygous SNPs")
    parser.add_argument("--gatk_vcf", required=True, help="Path to the GATK VCF file (must be bgzipped)")
    parser.add_argument("--phase_file", required=True, help="Path to the phased VCF file")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument("--metadata", required=True, help="Path to the metadata file")
    parser.add_argument("--rerun", action="store_true", help="Rerun all steps")
    parser.add_argument("--batch_size", type=int, default=5, help="Batch size for processing cells")
    args = parser.parse_args()
    # add example usage
    main(args)
    


