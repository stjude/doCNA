from collections import defaultdict


colorsCN = defaultdict (lambda: 'purple')
colorsCN['cn1'] = 'green'
colorsCN['cnL'] = 'blue'
colorsCN['cn3'] = 'red'
colorsCN['cnB'] = 'black'




def meerkat_plot (bed_df, axs, k_score_column = 'k_score', max_k_score = 10,  model_thr = 3):
    
    chrs = np.unique(bed_df['chrom'].tolist()).tolist()
    chrs.sort (key = lambda x: int(x[3:]))
    
    start = 0
   
    axs[0].plot ((start, start), (0, 0.95), 'k:', lw = 0.5)
    axs[1].plot ((start, start), (0, 4), 'k:', lw = 0.5)
    mids = []
    
    for chrom in chrs:
        for _, b in bed_df.loc[bed_df['chrom'] == chrom].iterrows():
            try:
                if b[k_score_column] <= 0:
                    a = 0.1
                elif b[k_score_column] > max_k_score:
                    a = 1
                else:
                    a = 0.1 + 0.9*b[k_score_column]/max_k_score 
            
                if b['model_score'] < model_thr:
                    color = colorsCN[b['model']]
                else:
                    color = 'magenta'
            
                axs[0].fill_between ((start + b['start'], start + b['end']), (b['k'], b['k']), color = color, alpha = a)
                axs[1].fill_between (x = (start + b['start'], start + b['end']), y1 = (b['cn'], b['cn']), y2 = (2, 2), color = color, alpha = a)
            except:
                pass
                
        end = bed.loc[bed['chrom'] == chrom, 'end'].max()
        mids.append (start + end / 2)
        start += end
        axs[1].plot ((start, start), (0.0, 4), 'k:', lw = 0.5)
        axs[0].plot ((start, start), (0.0, 0.95), 'k:', lw = 0.5)        
    
    axs[1].set_xticks (mids)
    axs[1].set_xticklabels (chrs, rotation = 60)
        
    axs[1].plot ((0, start), (2, 2), 'k--', lw = 1)        
    
    ranges = bed.loc[~(bed['k'].isnull())&((bed['end'] - bed['start']) > 1e6), ['k','m']].agg ([min, max])
    axs[0].set_ylim ((-0.009,ranges.loc['max','k']*1.1))
    axs[0].set_xlim ((-3e7, start + 3e7))
    
    axs[0].set_ylabel ('clonality')
    axs[1].set_ylabel ('coverage')



def chicken_feet_plot (bed_df, ax, k_score_column = 'k_score', max_k_score = 10, min_k_score = 2):
    ks = []
    ms = []
    for _,b in bed_df.loc[(bed_df['k_score'] > min_score)&(~bed_df['k_score'].isna())].iterrows():
              
        if b['k_score'] <= min_k_score:
            a = 0.1
        elif b['k_score'] > max_k_score:
            a = 1
        else:
            a = max(0, (b[k_score_column]-min_k_score)/10)
            
        s = 20+(b['end'] - b['start'])*3e-6
        
        ax.scatter (b['cn'], b['k'], 
                        s = s,
                        facecolor = mcl.to_rgba(colorsCN[b['model']] , alpha = a),
                        marker = 'o',
                        edgecolor = mcl.to_rgba(colorsCN[b['model']]),
                        lw = 0.5)
        
        ks.append (b['k'])
        ms.append (b['cn'])
        
    
    kmax = np.max(ks)
    ax.plot ((2, 2), (0,1.05), 'b:')
    k = np.linspace (0, 1.05, 2)
    ax.plot ((2+k), k, 'r:')
    ax.plot ((2-k), k, 'g:')
    ax.plot ((2+2*k), k, 'k:')
    ax.plot ((2-2*k), k, 'k:')


def leopard_plot ():
   pass #TBD


def eart_worm_plot (data_df, chrom, axs, m0, markersize = 2):
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
    
    axs[1].plot ((0, chromdata.position.max()), (m0, m0), 'k:')
    axs[1].set_ylim (chromdata['cov'].quantile ([0.005, 0.999]))

    

def chameleon_plot (bed_df, chrom, axs, markersize = 2, size_thr = 5, centromere_thr = 0.1):
    chrombed = bed_df.loc[chrom_df == chrom]    
    for _, seg in chrombed.loc[(chrombed.cent < centromere_thr)&(chrombed['size'] > size_thr)].iterrows():
        if (seg['k_score'] <= 0) | (np.isnan(seg['k_score'])):
            a = 0.1
        elif seg['k_score'] > 10:
            a = 0.9
        else:
            a = 0.1 + 0.8*seg['k_score']/10 
    
        if seg['model_score'] < 10:    
            axs[0].plot ((seg.start, seg.end), (seg.k, seg.k), c = colorsCN[seg.model], lw = 10, alpha = a)
            axs[0].plot ((seg.start, seg.end), (seg.k, seg.k), c = colorsCN[seg.model], lw = 1, marker = 'o')
            axs[1].plot ((seg.start, seg.end), (seg.cn, seg.cn), c = colorsCN[seg.model])
        else:
            axs[0].plot ((seg.start, seg.end), (seg.k, seg.k), c = 'magenta', lw = 10, alpha = a)
            axs[0].plot ((seg.start, seg.end), (seg.k, seg.k), c = 'magenta', lw = 1.5, marker = 'o', ls = ':')
            axs[1].plot ((seg.start, seg.end), (seg.cn, seg.cn), c = 'magenta', ls = ':')

    axs[1].set_ylim (chrombed.loc[(chrombed.cent < 0.1)&(chrombed['size'] > 1),'cn'].agg([min, max])*np.array ((0.95,1.05)))
    axs[1].plot ((0, chromdata.position.max()), (2, 2), 'k:')

    axs[0].set_title (name + ': ' + chrom)
    
