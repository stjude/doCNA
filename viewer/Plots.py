from collections import defaultdict, namedtuple
import matplotlib.colors as mcl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse as agp

colorsCN = defaultdict (lambda: 'purple')
colorsCN['AB+A'] = 'lime'
colorsCN['AB+AA'] = 'blue'
colorsCN['AB+AAB'] = 'darkred'
colorsCN['A(AB)B'] = 'black'
colorsCN['AB+AAAB'] = 'magenta'
colorsCN['AA+AAB'] = 'chocolate'
colorsCN['AAB+AABB'] = 'violet'
colorsCN['AB+AAA'] = 'brown'
colorsCN['AB+AAAA'] = 'darkolivegreen'
#copy-paste from Segment.py
def calculate_distance_old (preset, m, ai, m0):
    try:
        d = np.abs (preset.C (m/m0,ai,1) - preset.D (m/m0,ai,1))/np.sqrt (preset.A(m/m0,ai,1)**2 + preset.B(m/m0,ai,1)**2)
    except:
        d = np.inf 
    return d
    
def calculate_distance (preset, m, ai, m0):
    k = preset.k(m,ai,m0)
    if np.isnan(k):
        d = np.inf
    elif (k > 1) | (k < 0):
        ks = np.linspace (0,1,1000)
        ms = preset.m(k,m0)/m0
        d = np.min(np.sqrt((ks-k)**2+(ms-m/m0)**2))
    else:
        d = np.abs (preset.C (m/m0,ai,1) - preset.D (m/m0,ai,1))/np.sqrt (preset.A(m/m0,ai,1)**2 + preset.B(m/m0,ai,1)**2)
    
    return d
    
Preset = namedtuple ('Preset', ['A', 'B', 'C', 'D', 'k', 'm', 'ai'])
model_presets = {'AB+A' : Preset(A = lambda m,dv,m0: -m0/2,
                                 B = lambda m,dv,m0: -1,
                                 C = lambda m,dv,m0: m0,
                                 D = lambda m,dv,m0: m0*(2*dv/(0.5+dv))/2+m,
                                 k = lambda m,dv,m0: 2*dv/(0.5+dv),
                                 m = lambda k,m0: (2-k)*m0/2,
                                 ai = lambda k,m0: k/(2*(2-k))),
                 
                 'AB+AA' : Preset(A = lambda m,dv,m0: 0,
                                  B = lambda m,dv,m0: -1,
                                  C = lambda m,dv,m0: m0,
                                  D = lambda m,dv,m0: m,
                                  k = lambda m,dv,m0: 2*dv,
                                  m = lambda k,m0: np.repeat(m0, len(k)),
                                  ai = lambda k,m0: k/2),
                 
                 'AB+AAB' : Preset(A = lambda m,dv,m0: m0/2,
                                   B = lambda m,dv,m0: -1,
                                   C = lambda m,dv,m0: m0,
                                   D = lambda m,dv,m0: -m0*(2*dv/(0.5-dv))/2+m,
                                   k = lambda m,dv,m0: 2*dv/(0.5-dv),
                                   m = lambda k,m0: (2+k)*m0/2,
                                   ai = lambda k,m0: k/(2*(2+k))),

                 'A(AB)B' : Preset(A = lambda m,dv,m0: 0,
                                   B = lambda m,dv,m0: 1/2,
                                   C = lambda m,dv,m0: dv,
                                   D = lambda m,dv,m0: 0,
                                   k = lambda m,dv,m0: np.abs(m/m0 - 1),
                                   m = lambda k,m0: (1+k)*m0,
                                   ai = lambda k,m0: np.repeat(0, len(k))),#}
#end of copy-paste
#model_presets_4 = {
                 'AB+AAAB' : Preset (A = lambda m,dv,m0 : m0/2,
                                       B = lambda m,dv,m0 : -1,
                                       C = lambda m,dv,m0 : m0,
                                       D = lambda m,dv,m0 : m - 2*m0*dv/(1-2*dv),
                                       k = lambda m,dv,m0 : 2*dv/(1-2*dv),
                                       m = lambda k,m0 : (1+k)*m0,
                                       ai = lambda k,m0 : k/(2+2*k)),
                   
                   #'AAB+AAAB' : Preset (A = lambda m,dv,m0 : m0/2,
                   #                     B = lambda m,dv,m0 : -1,
                   #                     C = lambda m,dv,m0 : 3*m0/2,
                   #                     D = lambda m,dv,m0 : m- m0/((6*dv-1)/(2-4*dv)),
                   #                     k = lambda m,dv,m0 : (6*dv-1)/(1-2*dv),
                   #                     m = lambda k,m0 : (3+k)*m0/2,
                   #                     ai = lambda k,m0 : (1+k)/(6+2*k)),
                   
                   'AA+AAB' : Preset (A = lambda m,dv,m0 : m0/2,
                                      B = lambda m,dv,m0 : -1,
                                      C = lambda m,dv,m0 : m0,
                                      D = lambda m,dv,m0 : m- m0*(1-2*dv)/(2*dv+1),
                                      k = lambda m,dv,m0 : 2*(1-2*dv)/(2*dv+1),
                                      m = lambda k,m0 : (2+k)*m0/2,
                                      ai = lambda k,m0 : (2-k)/(4+2*k)),
                   
                    'AAB+AABB' : Preset (A = lambda m,dv,m0 : m0/2,
                                        B = lambda m,dv,m0 : -1,
                                        C = lambda m,dv,m0 : 3*m0/2,
                                        D = lambda m,dv,m0 : m- m0*(1-6*dv)/(4*dv+1),
                                        k = lambda m,dv,m0 : (1-6*dv)/(2*dv+1),
                                        m = lambda k,m0 : (3+k)*m0/2,
                                        ai = lambda k,m0 : (1-k)/(6+2*k)),
                    
                    'AB+AAA' : Preset (A = lambda m,dv,m0 : m0/2,
                                       B = lambda m,dv,m0 : -1,
                                       C = lambda m,dv,m0 : m0,
                                       D = lambda m,dv,m0 : m - 2*dv*m0/(3-2*dv),
                                       k = lambda m,dv,m0 : 4*dv/(3-2*dv),
                                       m = lambda k,m0 : (2+k)*m0/2,
                                       ai = lambda k,m0 : 3*k/(4+2*k)),
                    
                    'AB+AAAA' : Preset (A = lambda m,dv,m0 : m0/2,
                                        B = lambda m,dv,m0 : -1,
                                        C = lambda m,dv,m0 : m0,
                                        D = lambda m,dv,m0 : m - dv*m0/(1-dv),
                                        k = lambda m,dv,m0 : dv/(1-dv),
                                        m = lambda k,m0 : (1+k)*m0,
                                        ai = lambda k,m0 : k/(1+k))}

def meerkat_plot (bed_df, axs, chrom_sizes, max_k_score = 10, model_thr = 5):
    chrs = chrom_sizes.index.values.tolist()
    chrs.sort (key = lambda x: int(x[3:]))
    
    start = 0
    axs[1].plot ((start, start), (0, 4), 'k:', lw = 0.5)
    axs[0].plot ((start, start), (0, 0.95), 'k:', lw = 0.5)
    mids = []
    
    for chrom in chrs:
        for _, b in bed_df.loc[bed_df['chrom'] == chrom].iterrows():
            try:
                if b['k_score'] <= 0:
                    a = 0.1
                elif b['k_score'] > max_k_score:
                    a = 1
                else:
                    a = 0.1 + 0.9*b['k_score']/max_k_score
                    
                if b['model_score'] < model_thr:
                    color = colorsCN[b['model']]
                else:
                    color = 'yellow'
                k = b['k']
            except:
                color = 'lightskyblue'
                a = 0.1
            
            if np.isnan (b['k'])|np.isnan(a):
                a = 1
                k = 1.1
                color = 'red' #colorsCN[b['model']]
                
            axs[0].fill_between ((start + b['start'], start + b['end']), (k, k), color = color, alpha = a)
            axs[1].fill_between (x = (start + b['start'], start + b['end']), y1 = (b['cn'], b['cn']), y2 = (2, 2), color = color, alpha = a)
            #if b['model'] != 'cnB':
            axs[2].fill_between ((start + b['start'], start + b['end']), (b['k_score'], b['k_score']), color = color, alpha = a)
            
                
                
        end = chrom_sizes[chrom]#bed_df.loc[bed_df['chrom'] == chrom, 'end'].max()
        mids.append (start + end / 2)
        start += end
        axs[2].plot ((start, start), (-2.0, 2.), 'k:', lw = 0.5)
        axs[1].plot ((start, start), (0.0, 4), 'k:', lw = 0.5)
        axs[0].plot ((start, start), (0.0, 0.95), 'k:', lw = 0.5)        
    
    axs[-1].set_xticks (mids)
    axs[-1].set_xticklabels (chrs, rotation = 60)
        
    axs[1].plot ((0, start), (2, 2), 'k--', lw = 1)        
    
    ranges = bed_df.loc[~(bed_df['k'].isnull()), ['k','m']].agg ([min, max])
    axs[0].set_ylim ((-0.009, ranges.loc['max','k']*1.1))
    axs[0].set_xlim ((-3e7, start + 3e7))
    
    axs[1].set_ylim (bed_df.cn.agg([min,max]).values*np.array((0.9,1.1)))
    axs[0].set_ylim ((-0.009, bed_df.k.max()*1.1))  
    
    axs[0].set_ylabel ('clonality')
    axs[1].set_ylabel ('copy number')
    axs[2].set_ylabel ('score') 

def leopard_plot (bed_df, params, ax, highlight = ''):
    
    x = np.log10 (bed_df.loc[(bed_df.model != 'A(AB)B'), 'size'])
    y = np.log10 (bed_df.loc[(bed_df.model != 'A(AB)B'), 'k'])
    ax.plot (x, y, marker = '.', c = 'black', lw = 0, alpha = 0.5)

    x = np.log10 (bed_df.loc[(bed_df.model == 'A(AB)B'), 'size'])
    y = np.log10 (bed_df.loc[(bed_df.model == 'A(AB)B'), 'k'])
    ax.plot (x, y, marker = '.', c = 'lightgray', lw = 0, alpha = 0.5)

    x = np.log10 (bed_df.loc[(bed_df.model != 'cn')&(bed_df.status == 'CNV'), 'size'])
    y = np.log10 (bed_df.loc[(bed_df.model != 'cn')&(bed_df.status == 'CNV'), 'k'])
    ax.plot (x, y, marker = '.', c = 'red', lw = 0, alpha = 1)
    
    x = np.log10 (bed_df.loc[(bed_df.status == 'CNV-b'), 'size'])
    y = np.log10 (bed_df.loc[(bed_df.status == 'CNV-b'), 'k'])
    ax.plot (x, y, marker = '.', c = 'darkorange', lw = 0, alpha = 1)

    highlight_filter = [c in highlight for c in bed_df.chrom.tolist()]
    
    x = np.log10 (bed_df.loc[(highlight_filter), 'size'])
    y = np.log10 (bed_df.loc[(highlight_filter), 'k'])
    ax.plot (x, y, marker = 's', c = 'darkorange', lw = 0, alpha = 1, fillstyle = 'none')


    #ax.plot ((np.log10(5), np.log10(5)), (np.log10(1), np.log10(0.01)), 'k:')

    xt = np.linspace (-3, 2.5, 10)
    ax.plot (xt, -params['a']*xt - params['b'])

    ax.plot (xt, -params['a']*xt - params['bt'] , 'k:')
    #ax.set_ylim (-2.51, 0.1) 
    #ax.set_xlim (-2.01, 2.81) 


    ax.set_xlabel ('size (MB) / log')
    ax.set_ylabel ('clonality / log')

def chicken_feet_plot (bed_df, ax, highlight = '', k_score_column = 'k_score',
                       max_k_score = 10, model_thr = 5,
                       centromere_column = 'cent', centromere_thr = 0.3, size_thr = 1):
    ks = []
    ms = []
    bed_tmp = bed_df.loc[(bed_df[centromere_column] < centromere_thr)&(bed_df['size'] > size_thr)]
    for _,b in bed_tmp.loc[(~bed_tmp['k_score'].isna())].iterrows():
        if b['k_score'] <= 0:
            a = 0.1
        elif b['k_score'] > max_k_score:
            a = 1
        else:
            a = max(0, b[k_score_column]/max_k_score)
            
        s = 20+(b['end'] - b['start'])*3e-6
        
        if b['chrom'] == highlight:
            ax.scatter (b['cn'], b['k'], 
                        s = s, lw = 1.5,
                        marker = 's',
                        c = 'darkorange',
                        facecolor = None)
            
        ax.scatter (b['cn'], b['k'], 
                    s = s,
                    facecolor = mcl.to_rgba(colorsCN[b['model']] , alpha = a),
                    marker = 'o',
                    edgecolor = mcl.to_rgba(colorsCN[b['model']]),
                    lw = 0.5)
        
        
        
        ks.append (b['k'])
        ms.append (b['cn'])
        
    
    kmax = np.max(ks)
    ax.plot ((2, 2), (0,kmax), 'b:')
    k = np.linspace (0, kmax, 2)
    ax.plot ((2+k), k, 'r:')
    ax.plot ((2-k), k, 'g:')
    ax.plot ((2+2*k), k, 'k:')
    ax.plot ((2-2*k), k, 'k:')
    
    ax.set_xlabel ('copy number')
    ax.set_ylabel ('clonality')

def earth_worm_plot (data_df, bed_df, params, chrom, axs, k_score_column = 'k_score',
                     max_k_score = 10, markersize = 2, centromere_column = 'cent',
                     centromere_thr = 0.3, size_thr = 1):
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
    
    axs[1].plot ((0, chromdata.position.max()), (params['m'], params['m']), 'k:')
    axs[1].set_ylim (chromdata['cov'].quantile ([0.005, 0.999]))
    axs[0].set_ylabel ('BAF')
    axs[1].set_ylabel ('cn')
    
    chrombed = bed_df.loc[(bed_df.chrom == chrom)&(bed_df[centromere_column] < centromere_thr)&(bed_df['size'] > size_thr)]
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

    axs[3].set_ylim (chrombed['cn'].agg([min, max])*np.array ((0.95,1.05)))
    axs[3].plot ((0, chromdata.position.max()), (2, 2), 'k:')
    
    axs[2].set_ylabel ('clonality')
    axs[3].set_ylabel ('cn')

def check_solution_plot_opt (bed_df, params, ax, cent_thr = 0.3, size_thr = 1,
                             highlight = []):
    bed = bed_df.loc[(bed_df['cent'] < cent_thr)&(bed_df['size'] > size_thr)]
    
    for _, b in bed.iterrows():
    #    if (b['ai'] > 0.07) & (b['ai'] < 0.01):
    #        color = 'blue'
    #    else:
    #        color = 'darkgray'
        ax.scatter (b['m'],b['ai'], c = colorsCN[b['model']], s = b['size'])

    highlight_filter = [c in highlight for c in bed_df.chrom.tolist()]
    x = bed_df.loc[(highlight_filter), 'm'].values
    y = bed_df.loc[(highlight_filter), 'ai'].values
    ax.plot (x, y, marker = 's', c = 'darkorange', lw = 0, alpha = 1, fillstyle = 'none')
    
    ax.set_xlabel ('Coverage')
    ax.set_ylabel ('Allelic imbalance')
    
    
def verification_plot_CNV (data_chrom, bed_CNV, ax):
    filt = np.repeat (False, len(data_chrom))
    for _, r in bed_CNV.iterrows():
        filt = filt | ((data_chrom.position > r['start'])&(data_chrom.position < r['end']))
        
    tmp = data_chrom.loc[filt&(data_chrom['vaf'] < 1)&(data_chrom['vaf'] > 0)]
    vCNV = np.sort(tmp['vaf'].values)
    
    tmp = data_chrom.loc[~filt&(data_chrom['vaf'] < 1)&(data_chrom['vaf'] > 0)]
    vrest = np.sort(tmp['vaf'].values)

    v,c = np.unique (vCNV, return_counts = True)
    ax.plot (v, np.cumsum(c)/np.sum(c), '.', c = 'red', ms = 0.5, label = 'CNV region')

    v,c = np.unique (vrest, return_counts = True)
    ax.plot (v, np.cumsum(c)/np.sum(c), '.', c = 'blue', ms = 0.5, label = 'control region')

    ax.legend()
    ax.set_xlabel ('VAF')
    ax.set_ylabel ('CDF')