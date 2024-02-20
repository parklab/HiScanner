import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json

def draw_whole_genome(path):
    plt.rcParams['axes.xmargin'] = 0
    # read calls
    print(path)
    df = pd.read_csv(path,sep='\t')
    bins = df[['CHROM', 'START', 'END']].copy(deep=True)
    df.reset_index(inplace=True)
    df.drop(columns=['index'], inplace=True)
    scale_factor = df['gamma'].values[0]

    fig,ax=plt.subplots(2, figsize=(8, 3),sharex=True, dpi=200)

    # plot dotted lines for chromosomes
    chrlines = [0]*len(bins)
    for i in range(1, len(bins)):
        if bins.iloc[i]['CHROM'] != bins.iloc[i-1]['CHROM']:
            chrlines[i] = 1
    # get index
    chrlines = [i for i, x in enumerate(chrlines) if x == 1]
    # plot rectangles for every other chromosome
    for index in range(len(chrlines)):
        if index % 2 == 0:
            if index+1 >= len(chrlines):
                break
            ax[0].axvspan(chrlines[index], chrlines[index+1], facecolor='grey', alpha=0.1)
            ax[1].axvspan(chrlines[index], chrlines[index+1], facecolor='grey', alpha=0.1)

    # plot data and calls
    ax[0].scatter(df.index.values, scale_factor*df['RDR'].values,s=1,color='darkgrey')
    ax[0].plot(df.index.values, df['CN_total'].values,color='black',lw=.5,alpha=1)
    ax[1].scatter(df.index.values, df['pBAF'].values,s=1,color='black',alpha=.1)


    ax[1].set_ylim(-0.05,1.05)
    ax[0].set_ylim(0,10)
    ax[0].set_yticks([0,2,4,6,8,10])
    ax[1].set_ylabel('BAF')
    ax[0].set_ylabel('Copy Number')    
    plt.tight_layout()

    # set x axis margins to 0
    ax[0].set_xlim(0, len(df))
    ax[1].set_xlim(0, len(df))


    # remove x axis ticks
    ax[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax[1].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    return fig,ax

def plot_whole_genome_track(json_file_path, cells=[]):
    with open(json_file_path, 'r') as file:
        args = json.load(file)
    stem = args.get('stem')
    singlecell = args.get('singlecell').split(',')
    if len(cells) > 0:
        singlecell = cells
    WGD = args.get('MAX_WGD')
    for cell in singlecell:
        fig,ax = draw_whole_genome(f'{stem}/{cell}.long.wgd{WGD}.txt')
        ax[0].set_title(f'{cell}')
        fig.savefig(f'{stem}/{cell}.wgd{WGD}.png',
                    dpi=100)