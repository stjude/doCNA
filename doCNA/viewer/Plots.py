from collections import defaultdict, namedtuple
import matplotlib.colors as mcl
import matplotlib.pyplot as plt
import scipy.stats as sts
import pandas as pd
import numpy as np
import argparse as agp

from doCNA import Models
from doCNA import Consts

#colorsCN = {}
colorsCN = defaultdict (lambda: 'purple')
colorsCN['A'] = 'lime'
colorsCN['AA'] = 'blue'
colorsCN['AAB'] = 'cyan'
colorsCN['(AB)n'] = 'black'
colorsCN['AB'] = 'lightgray'
colorsCN['AAAB'] = 'magenta'
colorsCN['AAA'] = 'brown'
colorsCN['AAAA'] = 'darkolivegreen'
colorsCN[np.nan] = 'lightskyblue'
colorsCN['NA'] = 'lightskyblue'
    
    

def meerkat_plot (bed_df, axs, chrom_sizes, model_thr = 5, HE_thr = 3):
    chrs = chrom_sizes.index.values.tolist()
    chrs.sort (key = Consts.CHROM_ORDER.index)
    start = 0
    axs[1].plot ((start, start), (0, 4), 'k:', lw = 0.5)
    axs[0].plot ((start, start), (0, 0.95), 'k:', lw = 0.5)
    mids = []
    #axst = axs[2].twinx() 
   
    for chrom in chrs:
        for _, b in bed_df.loc[bed_df['chrom'] == chrom].iterrows():
            if b['p_model'] > 1/10**model_thr:
                color = colorsCN[b['model']]
            else:
                color = 'yellow'
            
            alpha = 0.1 + 0.9 * (b['score_HE']/HE_thr if b['score_HE'] <= HE_thr else 1)
                
            axs[0].fill_between ((start + b['start'], start + b['end']),
                                 (b['k'], b['k']), color = color, alpha = alpha)
            axs[1].fill_between (x = (start + b['start'], start + b['end']),
                                 y1 = (b['cn'], b['cn']), y2 = (2, 2), color = color, alpha = alpha)
            
            axs[2].fill_between ((start + b['start'], start + b['end']), 
                                     (b['score_HE'], b['score_HE']), color = color, alpha = alpha)
                
                
        end = chrom_sizes[chrom]#bed_df.loc[bed_df['chrom'] == chrom, 'end'].max()
        mids.append (start + end / 2)
        start += end
        axs[2].plot ((start, start), (0, 2.), 'k:', lw = 0.5)
        axs[1].plot ((start, start), (0.0, 4), 'k:', lw = 0.5)
        axs[0].plot ((start, start), (0.0, 0.95), 'k:', lw = 0.5)        
    
    axs[-1].set_xticks (mids)
    axs[-1].set_xticklabels (chrs, rotation = 60)
        
    axs[1].plot ((0, start), (2, 2), 'k--', lw = 1)        
    
    #ranges = bed_df.loc[~(bed_df['k'].isnull()), ['k','m']].agg ([min, max])
    maxk = max (bed_df.loc[~(bed_df['k'].isnull()), 'k'].max(), -bed_df.loc[~(bed_df['k'].isnull()), 'k'].min())
    
    #axs[0].set_ylim ((-0.009, 1.01))
    axs[0].set_ylim ((-0.009,  maxk *1.1))
    axs[0].set_xlim ((-3e7, start + 3e7))
    axs[1].set_ylim (bed_df.cn.agg([min,max]).values*np.array((0.9,1.1)))
        
    axs[0].set_ylabel ('clonality')
    axs[1].set_ylabel ('copy number')
    axs[2].set_ylabel ('HE score')
    #axst.set_ylabel ('score balanced') 
    #axst.set_ylim (0, axst.get_ylim()[1])
    axs[2].set_ylim (0, axs[2].get_ylim()[1])

def reporting_plot (bed_df, axs, chrom_sizes):
    chrs = chrom_sizes.index.values.tolist()
    #chrs.sort (key = lambda x: int(x[3:]))
    chrs.sort(key = Consts.CHROM_ORDER.index)
    
    start = 0
    axs[1].plot ((start, start), (0, 4), 'k:', lw = 0.5)
    axs[0].plot ((start, start), (0, 0.95), 'k:', lw = 0.5)
    mids = []
    
    for chrom in chrs:
        for _, b in bed_df.loc[bed_df['chrom'] == chrom].iterrows():
            a = 1
            color = colorsCN[b['model']]
            k = b['k']
            axs[0].fill_between ((start + b['start'], start + b['end']), (k, k), color = color, alpha = a)
            axs[1].fill_between (x = (start + b['start'], start + b['end']), y1 = (b['cn'], b['cn']), y2 = (2, 2), color = color, alpha = a)
            #if b['model'] != 'cnB':
            
                
                
        end = chrom_sizes[chrom]#bed_df.loc[bed_df['chrom'] == chrom, 'end'].max()
        mids.append (start + end / 2)
        start += end
        
        axs[1].plot ((start, start), (0.0, 4), 'k:', lw = 0.5)
        axs[0].plot ((start, start), (0.0, 0.95), 'k:', lw = 0.5)        
    
    axs[-1].set_xticks (mids)
    axs[-1].set_xticklabels (chrs, rotation = 60)
        
    axs[1].plot ((0, start), (2, 2), 'k--', lw = 1)        
    
    ranges = bed_df.loc[~(bed_df['k'].isnull()), ['k','m']].agg ([min, max])
    maxk = max (bed_df.loc[~(bed_df['k'].isnull()), 'k'].max(), -bed_df.loc[~(bed_df['k'].isnull()), 'k'].min())
    
    #axs[0].set_ylim ((-0.009, 1.01))
    axs[0].set_ylim ((-0.009,  maxk *1.1))
    axs[0].set_xlim ((-3e7, start + 3e7))
    axs[1].set_ylim (bed_df.cn.agg([min,max]).values*np.array((0.9,1.1)))
        
    axs[0].set_ylabel ('clonality')
    axs[1].set_ylabel ('copy number')

def leopard_plot (bed_df, params, ax, highlight = '', color_norm = 'black', color_hit = 'darkred', alpha = 1):
    
    a,b,bt = params
    
    x = np.log10 (bed_df['size'])
    y = np.log10 (np.abs (bed_df['k']))
    ax.plot (x, y, marker = 'o', c = color_norm, lw = 0, alpha = alpha)

    x = np.log10 (bed_df.loc[(bed_df.status != 'norm'), 'size'])
    y = np.log10 (np.abs(bed_df.loc[(bed_df.status != 'norm'), 'k']))
    ax.plot (x, y, marker = 'o', c = color_hit, lw = 0, alpha = alpha)
    

    for chrom in highlight:
        tmp = bed_df.loc[bed_df['chrom'] == chrom].sort_values (by = 'start')
        x = np.log10 (tmp['size'])
        y = np.log10 (np.abs(tmp['k']))
        ax.plot (x, y, marker = 's', c = 'magenta', lw = 1, alpha = 1, fillstyle = 'none')


    xt = np.linspace (-3, 2.5, 10)
    ax.plot (xt, -a*xt - b, c = color_norm)
    ax.plot (xt, -a*xt - bt , c = color_hit, linestyle = ':')

    ax.set_xlabel ('size (MB) / log')
    ax.set_ylabel ('clonality / log')

def plot_cdf (all_values, ax,  all_colors, par = (1,1), n = 100, xscale = 'lin', half = False):
    #l = 0.6*(max(values) - min(values))
    #x = np.linspace ((max(values) + min(values))/2 - l, (max(values) + min(values))/2 + l, n)
    xmax = par[0] + 3*par[1]
    
    if half:
        xmin = par[0]# - 3*par[1]
        x = np.linspace (xmin, xmax, n)
        y = 2*sts.norm.cdf (x, par[0], par[1]) - 1
    else:
        xmin = par[0] - 3*par[1]
        x = np.linspace (xmin, xmax, n)
        y = sts.norm.cdf (x, par[0], par[1])
    
    values = all_values[(all_values >= xmin)&(all_values <= xmax)]
    colors = all_colors[(all_values >= xmin)&(all_values <= xmax)]
    
    if xscale == 'lin':
        ax.scatter (np.sort(values), np.linspace (0.01,0.99, len(values)),
                    c = [colors[i] for i in np.argsort(values)])
        ax.plot (x, y, 'r-')
    elif xscale == 'log':
        ax.scatter (np.log10(np.sort(values)), np.linspace (0.01,0.99, len(values)),
                    c = [colors[i] for i in np.argsort(values)])
        ax.plot (np.log10(x), y, 'r-')
    else:
        raise ('Unknown scale')
        

def earth_worm_plot (data_df, bed_df, params, chrom, axs, markersize = 2, max_k_score = 10):
    chromdata = data_df.loc[data_df.chrom == chrom]

    chromdata.loc[chromdata['symbol'] == 'E'].plot(x = 'position', y = 'vaf', lw = 0, alpha = 0.3, color = 'orange', marker = '.', 
                                                   ms = markersize, ax = axs[0], legend = False)
    chromdata.loc[chromdata['symbol'] == 'U'].plot(x = 'position', y = 'vaf', lw = 0, alpha = 0.3, color = 'darkgray', marker = '.',
                                                   ms = markersize, ax = axs[0], legend = False)
    chromdata.loc[chromdata['symbol'] == 'N'].plot(x = 'position', y = 'vaf', lw = 0, alpha = 0.3, color = 'blue', marker = '.',
                                                   ms = markersize, ax = axs[0], legend = False)


    chromdata.loc[chromdata['symbol'] == 'N'].plot(x = 'position', y = 'cov', lw = 0, alpha = 0.3, color = 'red', marker = '.',
                                                   ms = markersize, ax = axs[1], legend = False)
    chromdata.loc[chromdata['symbol'] == 'U'].plot(x = 'position', y = 'cov', lw = 0, alpha = 0.3, color = 'darkorange', marker = '.',
                                                   ms = markersize, ax = axs[1], legend = False)
    chromdata.loc[chromdata['symbol'] == 'E'].plot(x = 'position', y = 'cov', lw = 0, alpha = 0.3, color = 'darkgray', marker = '.',
                                                   ms = markersize, ax = axs[1], legend = False)
    
    axs[1].plot ((0, chromdata.position.max()), (params['m0'], params['m0']), 'k:')
    axs[1].set_ylim (chromdata['cov'].quantile ([0.005, 0.999]))
    axs[0].set_ylabel ('BAF')
    axs[1].set_ylabel ('cn')
    
    chrombed = bed_df.loc[(bed_df.chrom == chrom)]
    for _, seg in chrombed.loc[(chrombed.cent < 0.5)&(chrombed['size'] > 1)].iterrows():
    
        if (seg['k_score'] <= 0) | (np.isnan(seg['k_score'])):
            a = 0.1
        elif seg['k_score'] >  max_k_score:
            a = 0.9
        else:
            a = 0.1 + 0.8*seg['k_score']/ max_k_score 
    
        if seg['model_score'] < 10:    
            axs[2].plot ((seg.start, seg.end), (seg.k, seg.k), c = colorsCN[seg.model], lw = 10, alpha = a)
            axs[2].plot ((seg.start, seg.end), (seg.k, seg.k), c = colorsCN[seg.model], lw = 1, marker = 'o')
            axs[3].plot ((seg.start, seg.end), (seg.cn, seg.cn), c = colorsCN[seg.model])
        else:
            axs[2].plot ((seg.start, seg.end), (seg.k, seg.k), c = 'magenta', lw = 10, alpha = a)
            axs[2].plot ((seg.start, seg.end), (seg.k, seg.k), c = 'magenta', lw = 1.5, marker = 'o', ls = ':')
            axs[3].plot ((seg.start, seg.end), (seg.cn, seg.cn), c = 'magenta', ls = ':')

    if len(chrombed) > 0:
        axs[3].set_ylim (chrombed['cn'].agg([min, max])*np.array ((0.95,1.05)))
        
    axs[3].plot ((0, chromdata.position.max()), (2, 2), 'k:')
    
    axs[2].set_ylabel ('clonality')
    axs[3].set_ylabel ('cn')

def check_solution_plot_opt (bed, ax, model_thr,
                             highlight = [], xcol = 'cn'):
    
    
    for _, b in bed.loc[bed['model'].notna(),:].iterrows():
        if b['chrom'] == 'chrX':
            ax.scatter (b[xcol],b['ai'], c = colorsCN[b['model']], s = b['size'],
                        edgecolor = 'w', marker = 'X')
        elif b['chrom'] == 'chrY':
            ax.scatter (b[xcol],b['ai'], c = colorsCN[b['model']], s = b['size']*2,
                        edgecolor = 'w', marker = 'v')
        else:
                        ec = 'w' if b['p_model'] > 1/10**model_thr else 'orange'
                        ax.scatter (b[xcol],b['ai'], c = colorsCN[b['model']],
                                    s = b['size'], 
                                    edgecolor = ec,
                                    marker = 'o')

    highlight_filter = [c in highlight for c in bed.chrom.tolist()]
    x = bed.loc[(highlight_filter), xcol].values
    y = bed.loc[(highlight_filter), 'ai'].values
    
    ax.plot (x, y, marker = 's', c = 'darkorange', lw = 0, alpha = 1, fillstyle = 'none')
    
    ax.set_xlabel ('Coverage/copy number')
    ax.set_ylabel ('Allelic imbalance')
    
stat_colors = {}
stat_colors['norm'] = 'blue'
stat_colors['CNVi'] = 'orange'
stat_colors['CNVb'] = 'black'

def verification_plot_CNV (d_ch, ch_bed, ax, par, type = 'CDF', no_bins = 100):
    assert type in ["CDF", "PDF"], "Unknown plot type!"
    for stat, df in ch_bed.groupby (by = 'status'):
        starts = df.start.values
        ends = df.end.values
        pos_filt = ((d_ch.position.values[:, np.newaxis] > starts[np.newaxis,:]) &(d_ch.position.values[:, np.newaxis] < ends[np.newaxis,:])).any (axis = 1) 
        #print (pos_filt.sum())
        tmp = d_ch.loc[(d_ch['vaf'] < 1)&(d_ch['vaf'] > 0)&(pos_filt)]
        stat
        v = np.sort (tmp['vaf'].values)
        if type == "CDF":
            ax.plot(v, np.linspace (0,1,len(v)), '.', markersize = 1, color = stat_colors[stat])
            ax.plot ((),(), lw = 2, label = stat, color = stat_colors[stat])
        else:
            ax.hist (v, bins = np.linspace (0,1, no_bins), lw = 2, 
                     histtype = "step", density = True, color = stat_colors[stat])
            ax.plot ((),(), lw = 2, label = stat, color = stat_colors[stat])
    
    try: 
        x = np.linspace (0,1, 500)
        if type == "CDF":
            ax.plot (x, sts.norm.cdf (x, 0.497, np.sqrt(0.25/par['m0'][0])*par['fb'][0]), '.',
                     markersize = 1, color = 'green')
        else:
            ax.plot (x, sts.norm.pdf (x, 0.497, np.sqrt(0.25/par['m0'][0])*par['fb'][0]), '.',
                     markersize = 1, color = 'green')
        ax.plot ((),(), lw = 2, color = 'green', label = 'diploid reference')
    except:
        pass
    
    ax.legend()
    ax.set_xlabel ('BAF', fontsize = 14)
    ax.set_ylabel (type, fontsize = 14)
