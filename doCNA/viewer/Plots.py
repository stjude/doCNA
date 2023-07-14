from collections import defaultdict, namedtuple
import matplotlib.colors as mcl
import matplotlib.pyplot as plt
import scipy.stats as sts
import pandas as pd
import numpy as np
import argparse as agp


from doCNA import Testing
from doCNA import Consts

#colorsCN = {}
colorsCN = defaultdict (lambda: 'purple')
colorsCN['A'] = 'lime'
colorsCN['AA'] = 'blue'
colorsCN['AAB'] = 'cyan'
colorsCN['(AB)(2+n)'] = 'black'
colorsCN['(AB)(2-n)'] = 'goldenrod'
colorsCN['AB'] = 'darkgray'
colorsCN['AAAB'] = 'magenta'
colorsCN['AAA'] = 'brown'
colorsCN['AAAA'] = 'darkolivegreen'
colorsCN['AA+AAB'] = 'lightcoral'
colorsCN['AAB+AAAB'] = 'peachpuff'
colorsCN['AAB+AABB'] = 'turquoise'
    

def meerkat_plot (bed_df, axs, chrom_sizes, model_thr = 5, HE_thr = 3):
    chrs = chrom_sizes.index.values.tolist()

    chrs.sort (key = Consts.CHROM_ORDER.index)
    start = 0
    axs[1].plot ((start, start), (0, 4), 'k:', lw = 0.5)
    axs[0].plot ((start, start), (0, 0.95), 'k:', lw = 0.5)
    mids = []
    
    max_he = bed_df.loc[np.isfinite(bed_df['score_HE'].values), 'score_HE'].max()

    for chrom in chrs:
        
        for _, b in bed_df.loc[bed_df['chrom'] == chrom].iterrows():
            if b['score_model'] < model_thr:

                if (b['score_HE'] < HE_thr) & (b['model'] != 'AB'):
                    color = 'r'
                else:
                    color = colorsCN[b['model']]

            else:
                color = 'yellow'
            
            alpha = 1#0.1 + 0.9 * (b['score_HE']/HE_thr if b['score_HE'] <= HE_thr else 1)
                
            axs[0].fill_between ((start + b['start'], start + b['end']),
                                 (b['k'], b['k']), color = color, alpha = alpha)
            axs[1].fill_between (x = (start + b['start'], start + b['end']),
                                 y1 = (b['cn'], b['cn']), y2 = (2, 2), color = color, alpha = alpha)
            
            if np.isposinf (b['score_HE']):
                ec = 'r'
                lw = 2
                ys = (max_he, max_he)
                hatch = '*'
                color = None
                alpha = 0.7
                                
            else:
                ec = color
                lw = 0
                ys = (b['score_HE'], b['score_HE'])
                hatch = None
                alpha = 1

            axs[2].fill_between ((start + b['start'], start + b['end']), 
                                 ys, edgecolor = ec, lw = lw,
                                 color = color, hatch  = hatch, alpha = alpha)
                
                
        end = chrom_sizes[chrom]
        mids.append (start + end / 2)
        start += end
        axs[2].plot ((start, start), (0, 2.), 'k:', lw = 0.5)
        axs[1].plot ((start, start), (0.0, 4), 'k:', lw = 0.5)
        axs[0].plot ((start, start), (0.0, 0.95), 'k:', lw = 0.5)        
    
    axs[-1].set_xticks (mids)
    axs[-1].set_xticklabels (chrs, rotation = 60)
        
    axs[1].plot ((0, start), (2, 2), 'k--', lw = 1)        
    
    maxk = max (bed_df.loc[~(bed_df['k'].isnull()), 'k'].max(), -bed_df.loc[~(bed_df['k'].isnull()), 'k'].min())
    
    axs[0].set_ylim ((-0.009,  maxk *1.1))
    axs[0].set_xlim ((-3e7, start + 3e7))
    axs[1].set_ylim (bed_df.cn.agg([min,max]).values*np.array((0.9,1.1)))
        
    axs[0].set_ylabel ('clonality')
    axs[1].set_ylabel ('copy number')
    axs[2].set_ylabel ('HE score')
    axs[2].set_ylim (0, axs[2].get_ylim()[1])

def reporting_plot (bed_df, axs, chrom_sizes):
    chrs = chrom_sizes.index.values.tolist()
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
           
        end = chrom_sizes[chrom]
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
        


def earth_worm_plot (data_df, bed_df, params, chrom, axs, markersize = 2,
                     max_score_HE = 10, model_threshold = 3):

    chromdata = data_df.loc[data_df.chrom == chrom]

    chromdata.loc[chromdata['symbol'] == 'E'].plot(x = 'position', y = 'vaf', lw = 0, alpha = 0.3,
                                                   color = 'orange', marker = '.', 
                                                   ms = markersize, ax = axs[0], legend = False)
    chromdata.loc[chromdata['symbol'] == 'U'].plot(x = 'position', y = 'vaf', lw = 0, alpha = 0.3,
                                                   color = 'darkgray', marker = '.',
                                                   ms = markersize, ax = axs[0], legend = False)
    chromdata.loc[chromdata['symbol'] == 'N'].plot(x = 'position', y = 'vaf', lw = 0, alpha = 0.3,
                                                   color = 'blue', marker = '.',
                                                   ms = markersize, ax = axs[0], legend = False)


    chromdata.loc[chromdata['symbol'] == 'N'].plot(x = 'position', y = 'cov', lw = 0, alpha = 0.3,
                                                   color = 'red', marker = '.',
                                                   ms = markersize, ax = axs[1], legend = False)
    chromdata.loc[chromdata['symbol'] == 'U'].plot(x = 'position', y = 'cov', lw = 0, alpha = 0.3,
                                                   color = 'darkorange', marker = '.',
                                                   ms = markersize, ax = axs[1], legend = False)
    chromdata.loc[chromdata['symbol'] == 'E'].plot(x = 'position', y = 'cov', lw = 0, alpha = 0.3,
                                                   color = 'darkgray', marker = '.',
                                                   ms = markersize, ax = axs[1], legend = False)
    
    axs[1].plot ((0, chromdata.position.max()), (params['m0'], params['m0']), 'k:')
    axs[1].set_ylim (chromdata['cov'].quantile ([0.005, 0.999]))
    axs[0].set_ylabel ('BAF')
    axs[1].set_ylabel ('cn')
    
    chrombed = bed_df.loc[(bed_df.chrom == chrom)]
    for _, seg in chrombed.loc[(chrombed.cent < 0.5)&(chrombed['size'] > 1)].iterrows():
    
        if (seg['score_HE'] <= 0) | (np.isnan(seg['score_HE'])):
            a = 0.1
        elif seg['score_HE'] >  max_score_HE:
            a = 0.9
        else:
            a = 0.1 + 0.8*seg['score_HE']/ max_score_HE 
    
        if seg['score_model'] < model_threshold:    
            axs[2].plot ((seg.start, seg.end), (seg.k, seg.k), c = colorsCN[seg.model], lw = 10, alpha = a)
            axs[2].plot ((seg.start, seg.end), (seg.k, seg.k), c = colorsCN[seg.model], lw = 1, marker = 'o')
            axs[3].plot ((seg.start, seg.end), (seg.cn, seg.cn), c = colorsCN[seg.model])
        else:
            axs[2].plot ((seg.start, seg.end), (seg.k, seg.k), c = colorsCN[seg.model], lw = 10, alpha = a)
            axs[2].plot ((seg.start, seg.end), (seg.k, seg.k), c = 'yellow', lw = 1.5, marker = 'o', ls = ':')
            axs[3].plot ((seg.start, seg.end), (seg.cn, seg.cn), c = 'yellow', ls = ':')

    if len(chrombed) > 0:
        axs[3].set_ylim (chrombed['cn'].agg([min, max])*np.array ((0.95,1.05)))
        
    axs[3].plot ((0, chromdata.position.max()), (2, 2), 'k:')
    
    axs[2].set_ylabel ('clonality')
    axs[3].set_ylabel ('cn')

def check_solution_plot_opt (bed, ax, model_thr,
                             highlight = [], xcol = 'cn'):
    
    
    for _, b in bed.loc[bed['model'].notna(),:].iterrows():

        ec = 'w' if b['score_model'] < model_thr else 'orange'
        if b['chrom'] == 'chrX':
            ax.scatter (b[xcol],b['ai'], c = colorsCN[b['model']], s = b['size']*2,
                        edgecolor = ec, marker = 'X')
        elif b['chrom'] == 'chrY':
            ax.scatter (b[xcol],b['ai'], c = colorsCN[b['model']], s = b['size']*4,
                        edgecolor = ec, marker = 'v')
        else:
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


def verification_plot_CNV (d_ch, ch_bed, ax, par, type = 'CDF', column = 'vaf', no_bins = 100):
    assert type in ["CDF", "PDF"], "Unknown plot type!"
    assert column in ['vaf', 'cov'], "Unknown column!"
    ##Plot reference
    
    vaf_down_lim, vaf_up_lim = sts.norm.ppf ((0.0005,0.9995), par['v0'], np.sqrt(0.25/par['m0']*par['fb']))
     
    if column == 'vaf':
        x = np.linspace (0,1, 500)
        ref_norm = sts.norm.cdf (x, par['v0'], np.sqrt(0.25/par['m0']*par['fb']))
        xlim = (-0.01,1.01)
    elif column == 'cov':
        x = Testing.lambda_ppf (np.linspace (0,1,500), par['m0'], par['l'])
        ref_norm = np.linspace (0,1,500)
        xlim = (d_ch['cov'].quantile ([0.005, 0.999]))
        
    if type == "CDF":
        ax.plot (x, ref_norm,
                 color = 'red', lw = 0.5)
        ax.set_xlim (xlim)
    
    ax.plot ((),(), lw = 2, color = 'red', label = 'diploid reference')
    
    ##Plot AB
    dipl_bed = ch_bed.loc[ch_bed['model'] == 'AB']
    starts = dipl_bed['start'].values
    ends = dipl_bed['end'].values
    pos_filt = ((d_ch.position.values[:, np.newaxis] > starts[np.newaxis,:]) &\
                (d_ch.position.values[:, np.newaxis] < ends[np.newaxis,:])).any (axis = 1) 
    
    tmp = d_ch.loc[(d_ch['symbol'] == Consts.E_SYMBOL)&(pos_filt)]
        
    v,c = np.unique (tmp[column].values, return_counts = True)
    if type == "CDF":
        ax.plot(v, np.cumsum(c)/np.sum(c), '.', markersize = 0.5, color = colorsCN['AB'])
        ax.plot ((),(), lw = 2, label = 'AB', color = colorsCN['AB'])
    else:
        ax.hist (v, bins = np.linspace (0,1, no_bins), lw = 2, 
                 histtype = "step", density = True, color = colorsCN['AB'])
        ax.plot ((),(), lw = 2, label = 'AB', color = colorsCN['AB'])
    
    ##Plot CNVs
    CNV_bed = ch_bed.loc[ch_bed['model'] != 'AB']
    for _, cb in CNV_bed.iterrows():
        #(d_ch['symbol'] == cb['symbol'])&\
        
        tmp = d_ch.loc[(d_ch['vaf'] < vaf_up_lim)&(d_ch['vaf'] > vaf_down_lim)&\
                       (d_ch['position'] >= cb['start'])&\
                       (d_ch['position'] <= cb['end'])]
        v,c = np.unique (tmp[column].values, return_counts = True)
        if type == "CDF":
            ax.plot(v, np.cumsum(c)/np.sum(c), '.', markersize = 1,
                    color = colorsCN[cb['model']])
            ax.plot ((),(), lw = 2, color = colorsCN[cb['model']],
                     label = cb['chrom']+':'+str(cb['start'])+'-'+str(cb['end'])+':'+cb['model'])
        else:
            ax.hist (v, bins = np.linspace (0,1, no_bins), lw = 2, 
                     histtype = "step", density = True,
                     color = colorsCN[cb['model']])
            ax.plot ((),(), lw = 2, color = colorsCN[cb['model']],
                     label = cb['chrom']+':'+str(cb['start'])+'-'+str(cb['end'])+':'+cb['model'])
    
    ax.legend()