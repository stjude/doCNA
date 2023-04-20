from shiny import *
from .Plots import * #need a dot before Plots, like this: .Plots

from doCNA import Models
from doCNA import Consts
from doCNA import Scoring

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as sig
from collections import defaultdict

chromdic = {}
for c in Consts.CHROM_ORDER:
    chromdic[c] = c

def merge_records (all_records, chrom):
    
    filt_indexes = np.where([r.filt for r in all_records])[0]
    #records = [all_records[i] for i in filt_indexes]
    if len(filt_indexes) == 0:
        r = all_records[0]
        records = all_records
        status = 'NA'
        size = np.array ([r.end-r.start for r in records])
        m = np.array([r.m for r in records])
        cn = np.array([r.cn for r in records])
        k = np.abs(np.array([r.k for r in records]))
        #chi2 = (np.array([r.k_score for r in records])*np.log2(10))
        cyto = '-'.join([records[0].cyto.split('-')[0], records[-1].cyto.split('-')[-1]])
        model = 'NA'
        mr = (chrom, r.start, records[-1].end, 
              sum(m*size)/sum(size), sum(cn*size)/sum(size), model,
              sum(k*size)/sum(size), cyto, 'NA')
    else:
        
        records = [all_records[i] for i in filt_indexes]
        i = 0
        try:
            while not records[i].filt:
                i +=1 
        except:
            i = 0

        r = records[i]
        mr = ()
        if len (records) == 1:
            mr = (chrom, r.start, r.end, r.m, r.cn,
                r.model if r.status != 'norm' else 'AB',
                np.abs(r.k) if r.status != 'norm' else 0, r.cyto, r.k_score)
        else:
            status = r.status
            size = np.array ([r.end-r.start for r in records])
            m = np.array([r.m for r in records])
            cn = np.array([r.cn for r in records])
            k = np.abs(np.array([r.k for r in records]))
            chi2 = (np.array([r.k_score for r in records])*np.log2(10))
            cyto = '-'.join([records[0].cyto.split('-')[0], records[-1].cyto.split('-')[-1]])
            if status == 'norm':
                model = 'AB'
            else:
                model = r.model
            
            mr =  (chrom, records[0].start, records[-1].end,
                    sum(m*size)/sum(size), sum(cn*size)/sum(size), model,
                    sum(k*size)/sum(size), cyto, -np.log10(sts.chi2.sf(2*sum(chi2), 2*len(records))))
    
    return mr

  
app_ui = ui.page_fluid(
    ui.h2 ({"style" : "text-align: center;"}, "doCNA results viewer."),
   
    ui.layout_sidebar(ui.panel_sidebar(ui.h4 ("Segments filtering:"),
                                       ui.input_slider ('cent_thr', "Centromere fraction threshold",
                                                        value = 0.3, min = 0, max = 1),
                                       ui.input_slider ('size_thr', "Min segment size (in Mb)",
                                                        value = 5, min = 0, max = 10),
                                       ui.h4 ("Display settings:"),
                                       ui.input_slider ('model_thr', "Model score threshold",
                                                        value = 3, min = 0, max = 10),
                                       ui.input_slider ('HE_max', "Max HE score:",
                                                        value = 2, min = 0, max = 10),
                                       
                                       width = 2),
                      ui.panel_main(
                                    ui.navset_tab (
                                                   ui.nav("Genome-wide view",
                                                          ui.row(ui.column(12,          
                                                                 ui.input_file ('bed_file', "Choose BED file to upload:",
                                                                                multiple = False, accept = '.bed'),
                                                                 ui.output_plot ('genome_plot'),)),
                                                          ui.row(ui.column(12,
                                                                           ui.h5 (''),
                                                                           ui.input_checkbox_group ('chroms_selected',
                                                                                                    "Select chromosomes to highlight",
                                                                                                    chromdic, inline = True)),),   
                                                          ui.row (ui.input_file ('par_file', "Choose PAR file to upload:",
                                                                                          multiple = False, accept = '.par')),
                                                          ui.row(ui.column(6,
                                                                           ui.h5 ('Solution review:'),
                                                                           ui.output_plot ('solution_plot'),
                                                                           ),
                                                                 ui.column(6,
                                                                           ui.row(
                                                                                  ui.h5 ('Scoring review:'),
                                                                                  ui.output_plot ('scoring_plots'),),
                                                                           )),
                                                         ),
                                                   
                                                   ui.nav("LOG", 
                                                          ui.row(ui.column(12,
                                                                 ui.input_file ('log_file', "Choose LOG file to screen:",
                                                                                multiple = False, accept = '.log'),
                                                                 ui.output_ui("dyn_log_ui")))),
 
                                                   ui.nav("CNVs",
                                                         ui.row(ui.column (2,
                                                                           ui.input_radio_buttons ('sort_CNV_by',
                                                                                                   "Sort list by:",
                                                                                                   {'position':'position', 'score':'score'}, inline = True)),
                                                                ui.column (2, 
                                                                           ui.input_radio_buttons ('corrected',
                                                                                                   'Show FDR?:',
                                                                                                   { 'status':'FDR adj', 'status_d':'Score'}, inline = True))),
                                                         ui.row(ui.output_text ("number_CNVs"), ""),
                                                         ui.row(ui.output_table (id = 'CNVs'),)),
                                                   ui.nav("Solution test",
                                                          ui.layout_sidebar(ui.panel_sidebar(ui.h4 ("Optimize settings:"),
                                                                                             ui.input_slider ('min_cn', "Min cov (relative):",
                                                                                                              value = 0.7, min = 0.1, max = 1, step = 0.05),
                                                                                             ui.input_slider ('max_cn', "Max cov (relative):",
                                                                                                              value = 1.1, min = 0.5, max = 1.5, step = 0.05),
                                                                                             ui.row(ui.column(6,ui.input_numeric ("step", 'Step:', 0.5, min = 0.5, max = 5, step = 0.5)),
                                                                                                    ui.column(6,#ui.h6 ("# steps"),
                                                                                                                ui.input_text ("number_points", '# points')),
                                                                                                   ),    
                                                                                             ui.h6("Coverage range:"),
                                                                                             ui.output_text_verbatim ("coverage_range"),
                                                                                             ui.input_checkbox ('all_models', 'Show all models?',
                                                                                                                value = False),
                                                                                             ui.input_action_button ('opt',
                                                                                                                     "Optimize solution"),
                                                                                             width = 2),
                                                                            ui.panel_main(ui.h6("Solution check plot:"),
                                                                                          ui.output_plot('solution_plot_opt', "Solution plot"),
                                                                                          ui.h6("Total weighted distance to solution at coverage:"),
                                                                                          ui.output_plot('opt_plot', "Relative distance plot"),
                                                                                          ui.h6("Pick coverage to plot models:"),
                                                                                          ui.row(ui.input_slider ('m0_cov', "Diploid coverage",
                                                                                                               value = 0.5, min = 0, max = 1, width = '200%')),
                                                                                                 ui.output_text_verbatim ('solutions')))),
                                                                                          
                                                   ui.nav("Chromosome view",
                                                          ui.row(ui.column(12,
                                                                           ui.input_file ('data_file', "Choose data file to upload:",
                                                                                           multiple = False, accept = ('.dat', '.dat.gz')),
                                                                           ui.input_radio_buttons ('chrom_view', "Choose chromosome to inspect",                                                                                    
                                                                                                   chromdic, inline = True),)),
                                                          ui.row(ui.input_radio_buttons('f_to_plot', "Which function plot to compare:", 
                                                                                         {'PDF' : 'Probability density function',
                                                                                          'CDF' : 'Cumulative density function /for the pros/'}, 
                                                                                         selected = 'CDF', inline = True)),                                                          
                                                          ui.row(ui.output_plot ('data_plot'),
                                                                 ui.output_plot ('compare_plot')),
                                                          ui.row(ui.output_table (id = 'chrom_segments'))),
                                                    ui.nav("Report",
                                                           ui.row(ui.column(12,
                                                                            ui.output_plot('report_plot'),
                                                                            ui.row (ui.h6 ("Report diploid regions?"),
                                                                                    ui.input_checkbox ('rep_AB', "Yes", value = False)),
                                                                            ui.output_table('report')))),
                                                    ui.nav("Publication",
                                                           ui.h6 ("AI still in school"))
                                                            
                                                   )
                                   )
          )
        ) 

def server(input, output, session):
    
    bed_full = reactive.Value(pd.DataFrame())
    bed = reactive.Value(pd.DataFrame())
    data = reactive.Value(pd.DataFrame())
    par = reactive.Value({})
    opt_solution = reactive.Value((np.array([]), np.array([]), np.array([])))
    chrom_sizes = reactive.Value(pd.Series())
    m0 = reactive.Value(np.nan)
    m0_opt = reactive.Value(np.nan)
    log_file = reactive.Value ([])
    bed_report = reactive.Value(pd.DataFrame())
    model_presets = reactive.Value (Models.model_presets)
    
          
    @output
    @render.text
    def solutions ():
        m, d, _, _ = opt_solution ()
        ind = sig.argrelmin (d)[0]
        minims = []
        for i in ind:
            minims.append ((m[i], d[i])) 
        minims.sort (key = lambda x: x[1], reverse = False)
        
        return '\n'.join(['m = ' + str(m) + '  d = ' + str(d)  for m,d in minims])
  
   
    @output
    @render.ui
    def dyn_log_ui():
        return ui.TagList (ui.tags.textarea ([l.strip() for l in log_file()], 
                                             cols = "150", rows = "30"))
    
    @output
    @render.text
    def coverage_range ():
        if np.isnan(m0()):
            text = 'n.a.'
        else:
            text = '{:.2f}'.format(int(m0()*input.min_cn())) + ' - ' + '{:.2f}'.format(int(m0()*input.max_cn()))
        return text
    
    @reactive.Effect
    def _():
        try:
            text = str(len(np.arange(m0()*input.min_cn(), m0()*input.max_cn(), input.step())))
            ui.update_text ('number_points', value = text)
        except:
            pass
    
    
    @reactive.Effect
    @reactive.event (input.log_file)
    def _():
        file_input = input.log_file()
        if not file_input:
            return
        with open(file_input[0]['datapath']) as f:
            log = f.readlines()
        log_file.set(log)    
        
    @reactive.Effect
    @reactive.event(input.bed_file)
    def _():
        file_input = input.bed_file()
        
        if not file_input:
            return
        df = pd.read_csv (file_input[0]['datapath'], sep = '\t', header = None, 
                   names = ['chrom', 'start', 'end', 'ai','p_ai', 'm', 'cn',
                            'd_HE', 'score_HE', 'model', 'd_model', 'score_model',
                            'k', 'symbol', 'cyto', 'cent'])
        
        df['size'] = (df['end'] - df['start'])/1e6
        
        print (df.head())
        print (df.columns)
        
        data.set(pd.DataFrame())
        par.set({})
        bed_full.set(df)
        opt_solution.set ((np.array([]), np.array([]), np.array([])))
        log_file.set([])
        bed_report.set(pd.DataFrame())
            
    @reactive.Effect
    @reactive.Calc
    def _():
        tmp = bed_full()
        if len(tmp) > 0:
            chrom_sizes.set(tmp.groupby(by = 'chrom').agg({'end' : 'max'})['end'])
            tmp['filt'] = (tmp['cent'] <= input.cent_thr()) & (tmp['size'] >= input.size_thr())
            
            #bed_full.set(tmp)
            bed.set (tmp.loc[tmp.filt])
    
    #@reactive.Effect
    #@reactive.Calc
    #def _():
    #    bf = bed_full()    
    #    b = bed()
    #    if (len(bf) != 0) & (len(b) != 0):
    #        chrs = chrom_sizes().index.values.tolist()
    #        #chrs.sort (key = lambda x: int(x[3:]))
    #        chrs.sort (key = Consts.CHROM_ORDER.index)
    #        merged_segments = []
    #        for chrom in chrs: #
    #            segments = bf.loc[bf.chrom == chrom] #, segments in bf.groupby (by = 'chrom'):
    #            data = segments.sort_values (by = 'start', ignore_index = True)
    #            
    #            seg_iter = data.itertuples()
    #            
    #            to_merge = [next(seg_iter)]
    #            
    #            try:
    #                while not to_merge[-1].filt:
    #                    to_merge.append (next(seg_iter))
    #            except:
    #                pass
    #            current_record = to_merge[-1]
    #            
    #            last_action = 'merge' 
    #            while True:
    #                try:
    #                    
    #                    next_record = next(seg_iter)
    #                    while not next_record.filt:
    #                        
    #                        to_merge.append (next_record)
    #                        next_record = next(seg_iter)
    #                    
    #                    if (current_record.status == next_record.status)&\
    #                             (current_record.model == next_record.model):
    #                        if (current_record.status == 'norm'):
    #                       
    #                            to_merge.append(next_record)
    #                            last_action = 'merge'
    #                        elif (current_record.model == next_record.model)&\
    #                             np.abs(((current_record.k-next_record.k)/current_record.k) < 0.1):
    #                            to_merge.append(next_record)
    #                            last_action = 'merge'
    #                    
    #                    else: 
    #                        
    #                        merged_segments.append(merge_records(to_merge, chrom))
    #                        to_merge = [next_record]
    #                        current_record = next_record
    #                        last_action = 'no merge'
    #                except StopIteration:
    #                    break
    #                 
    #            if last_action == 'no merge':
    #                merged_segments.append(merge_records(to_merge, chrom))
    #            
    #        bed_report.set(pd.DataFrame.from_records (merged_segments,
    #                                                  columns = ['chrom', 'start', 'end', 'm', 'cn','model', 'k', 'cyto', 'score']))
    
    @reactive.Effect
    @reactive.event(input.par_file)
    def _():
        file_input = input.par_file()
       
        if not file_input:
            return
        pard = {}
        with open(file_input[0]['datapath'],'r') as f:
            for line in f.readlines():
                try:
                    key, value = line.split('\t')
                    pard[key] = float(value)
                except:
                    try:
                        key, values = line.split('\t')
                        pard[key] = [v.strip("',[]") for v in values.split(' ')]
                    except:
                        print ('Line: ' + line + 'not parsed.')
        par.set(pard)
        print ('')
        print(pard)
        print('')
        opt_solution.set ((np.array([]), np.array([]), np.array([])))
        m0.set(pard['m0'])
        m0_opt.set(pard['m0'])
    
            
    @reactive.Effect
    @reactive.event(input.data_file)
    def _():
        file_input = input.data_file()
        
        if not file_input:
            return
        df = pd.read_csv (file_input[0]['datapath'], sep = '\t')    
        data.set(df)
    
    @output
    @render.plot (alt = "Genomewide view")
    def genome_plot ():
        bed_data = bed()
        if  len(bed_data):
            fig, axs = plt.subplots (3, 1, figsize = (16,6), sharex = True)
            meerkat_plot (bed_data, axs, chrom_sizes(),
                          model_thr = input.model_thr(), HE_thr = input.HE_max())
            
            for model in model_presets().keys():
                axs[0].plot ((),(), lw = 10, color = colorsCN[model], label = model)
            axs[0].plot ((),(), lw = 10, color = 'yellow', label = 'complex')
            axs[0].plot ((),(), lw = 10, color = 'red', label = 'fail')
            axs[0].legend (bbox_to_anchor = (0.5, 2), ncol = len(model_presets())+2,
                           loc = 'upper center', title = 'Models of mixed clones: normal (AB) and abnormal karyotypes:')
            return fig
        
    @output
    @render.plot (alt = "Genomewide view")
    def report_plot ():
        bed_data = bed_report ()
        if len (bed_data):
            fig, axs = plt.subplots (2, 1, figsize = (16,4), sharex = True)
            reporting_plot (bed_data, axs, chrom_sizes())
            for model in model_presets().keys():
                axs[0].plot ((),(), lw = 10, color = colorsCN[model], label = model)
            axs[0].plot ((),(), lw = 10, color = 'yellow', label = 'complex')
            axs[0].plot ((),(), lw = 10, color = 'red', label = 'fail')
            axs[0].legend (bbox_to_anchor = (0.5, 2), ncol = len(model_presets())+2,
                           loc = 'upper center', title = 'Models of mixed clones: normal (AB) and abnormal karyotypes:')
            return fig
    
    @output
    @render.plot (alt = "Scoring view")
    def scoring_plots ():
        bed_data = bed()
        par_d = par()
        
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, axs = plt.subplots (2, 1, figsize = (6,6))
                        
            plot_cdf (bed_data['d_HE'].values, axs[0], par = (par_d['m_d'],par_d['s_d']),
                      all_colors = np.array([colorsCN[m] for m in bed_data['model']]), half = True)
            axs[0].set_ylabel ('cdf - HE distance')

            tmp_bed = bed_data.loc[bed_data['model'] != 'AB'].sort_values (by = 'd_model')
            axs[1].scatter (tmp_bed['d_model'].values, np.linspace (0,1, len(tmp_bed)),
                            c = np.array([colorsCN[m] for m in tmp_bed['model']]),
                            s = np.sqrt(tmp_bed['size']))
            x = np.linspace (tmp_bed['d_model'].min(), tmp_bed['d_model'].max(), 100)
            axs[1].plot (x , 1 - np.exp (-par_d['a_d'] * x), 'r-')
            axs[1].set_ylabel ('cdf - Model distance')  
            
            fig.tight_layout()
            
            return fig
    
    @output
    @render.plot (alt = "Solution view")
    def solution_plot ():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, ax = plt.subplots (1, 1, figsize = (6,6))
            check_solution_plot_opt (bed_data, ax, 
                                     highlight = input.chroms_selected(),
                                     model_thr = input.model_thr())
            k = np.linspace (0,1,100)
            m0 = par_d['m0']
            for model in model_presets().keys():
                if model in par()['models']:
                    ls = '-'
                else:
                    ls = ':'
                ax.plot (2*model_presets()[model].m(k, m0)/m0, model_presets()[model].ai(k, m0),  
                         lw = 1, linestyle = ls, color = colorsCN[model], alpha  = 1)
            
            ax.set_xlim (2*0.9*bed_data.m.min()/m0, 2*1.1*bed_data.m.max()/m0)
            ax.set_ylim ((max(-0.02, -0.02*bed_data.ai.max()), bed_data.ai.max()*1.1))
            
            return fig
    
        
    
    @output
    @render.plot (alt = "Solution view")
    def solution_plot_opt ():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
           
            fig, ax = plt.subplots (1, 1, figsize = (6,6))
            check_solution_plot_opt (bed_data, ax, model_thr = input.model_thr(),
                                          highlight = [], xcol = 'm')
            k = np.linspace (0,1,100)
            #m0 = m0_opt()
            for model in model_presets().keys():
                ax.plot (model_presets()[model].m(k, m0_opt()), model_presets()[model].ai(k, m0_opt()), 
                         lw = 1, linestyle = ':', color = colorsCN[model], alpha = 0.6, label = model)
            
            #ax.legend (bbox_to_anchor = (1,0))
            ax.legend (bbox_to_anchor = (1.2,1), loc = 'upper center')
               
            return fig
    
    @reactive.Effect
    @reactive.event(opt_solution)
    def _():
        if len(opt_solution()[0]):
            ui.update_slider ('m0_cov', value = m0_opt(),
                              min = opt_solution()[0][0],
                              max = opt_solution()[0][-1],
                              step = input.step() )
    
    @output
    @render.table
    def chrom_segments ():
        bed_data = bed()
        if (len(bed_data) != 0):
            return bed_data.loc[bed_data.chrom == input.chrom_view()].sort_values(by = 'start')

    @output
    @render.table
    def report():
        report = bed_report()
        if len(report) > 0:
            if not input.rep_AB():
                return report.loc[report.model != 'AB']
            else:
                return report

    @output
    @render.table
    def CNVs ():
        bed_data = bed()
        if (len(bed_data) != 0):
            
            if input.sort_CNV_by() == 'score':
                tmp_bed = bed_data.loc[bed_data['model'] != 'AB'].sort_values(by = 'score_HE', ascending = False)
            else:
                tmp_bed = bed_data.loc[bed_data['model'] != 'AB']
            return tmp_bed

    @output
    @render.text
    def number_CNVs():
        bed_data = bed()
        if len(bed_data) > 0: 
            message = "Number of CNVs found: " + str(len(bed_data.loc[bed_data['model'] != 'AB']))
        else:
            message = ''
        return  message

    @output
    @render.plot
    def data_plot ():
        bed_data = bed()
        par_d = par()
        data_df = data()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0) & (len(data_df) != 0):
            fig, axs = plt.subplots (4, 1, figsize = (12,4), sharex = True)
            earth_worm_plot (data_df, bed_data, par_d, input.chrom_view(), axs)
            return fig
        
    @output
    @render.plot
    def compare_plot():
        bed_data = bed()
        data_df = data()
        if (len(bed_data) != 0) & (len(data_df) != 0):
            CNV_bed = bed_data.loc[bed_data.chrom == input.chrom_view()]
            data_chrom = data_df.loc[data_df.chrom == input.chrom_view()]
            fig, ax = plt.subplots (figsize = (6,8))
            verification_plot_CNV (data_chrom, CNV_bed, ax, par(), input.f_to_plot())
            return fig
        
    ###optimalization below    
    @reactive.Effect
    @reactive.event (input.opt)
    def _():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
           
            ms = np.arange (m0()*input.min_cn(), m0()*input.max_cn(), input.step())
            scorers = []
            #tmp = bed_data.loc[(~bed_data['model'].isna())]
            ai = bed_data['ai'].values
            m_cov = bed_data['m'].values
            
            ai_index = ai < Consts.DIPLOID_AI_THR
            
            with ui.Progress (min = ms[0], max = ms[-1]) as p:
                p.set(message = 'Optimizing solution...', )
                                            
                for m in ms:
                    p.set(m, message = 'Calculating')
                    
                    cn_index = (m_cov > (1-Consts.DIPLOID_dCN_THR)*m/2) & (m_cov < (1+Consts.DIPLOID_dCN_THR)*m/2)
                    #scorers.append (Scoring.Scoring(np.array(ai[], m_cov[]), logger = False))
                    for _, b in bed_data.iterrows ():
                        pass  
                    
                    res = Models.fit_exp_shift (bed_data['ai'].values, 2*bed_data['m'].values/m)
                    #for _, b in bed_data.iterrows():
                        #dt.append (np.sqrt(b['size'])*(np.nanmin([Models.calculate_distance(model, b['m'], b['ai'], m) for model in model_presets().values()])))
                    #    dt.append ((b['size'])*(np.min([Models.calculate_distance_new (b['ai'], 2*b['m']/m, model)['d'] for model in model_presets().values()])))                    
                    #index = np.where(np.isfinite(dt))[0]
                    #sts = np.sum(([st[i] for i in index]))
                    #fts.append (sts/sttotal)
                    #dts.append (np.sum(dt))#/sts)
                                        
                #fraction = np.array (fts)
                #dtsa = np.array(dts)
            
            #opt_solution.set((ms,dtsa, np.array(dist_as), np.array(dist_bs)))
            
            try:
                m0_opt.set (ms[np.where(dts == np.nanmin(dtsa))[0][0]])
                #print ()
                #print (ms[np.where(dts == np.nanmin(dtsa))[0][0]])
                #print ()
            except:
                m0_opt.set(m0())
            
    @output
    @render.plot (alt = 'Total distance to the model')
    def opt_plot ():
        
        if len(opt_solution()[0]):
            opt = opt_solution()
            fig, ax = plt.subplots(1,1, figsize = (6,6))
            ax.plot (opt_solution()[0], opt_solution()[1], 'k-', label = 'Total distance')
            #ax.plot (opt_solution()[0], opt_solution()[1]/opt_solution()[2], 'r-', label = 'Normed total distance')
            axt = ax.twinx()
            axt.plot (opt_solution()[0], opt_solution()[2], 'b:', label = 'Diploid distance')
            
            axtt = ax.twinx()
            axtt.spines.right.set_position(("axes", 1.15))

            axtt.plot (opt_solution()[0], opt_solution()[3], 'r:', label = 'Diploid shift')
            
            ax.set_xlabel ('Covearage')
            ax.set_ylabel ('Relative distance to the model')
            axt.set_ylabel ('Diploid distance')
            axt.yaxis.label.set_color ('b')
            axtt.set_ylabel ('Diploid shift')
            axtt.yaxis.label.set_color ('b')
            
            axt.yaxis.label.set_color('b')
            axtt.yaxis.label.set_color('r')

            axt.tick_params(axis='y', colors = 'b')
            axtt.tick_params(axis='r', colors = 'b')
            
            #for model in colorsCN.keys():
            #    ax.plot ((),(), lw = 2, color = colorsCN[model], label = model)
            ax.legend()
            #ax.legend (bbox_to_anchor = (1.4,1), loc = 'upper center')
            return fig
               
    @reactive.Effect
    @reactive.event (input.m0_cov)
    def _():
        if len(opt_solution()[0]):
            m0_opt.set(input.m0_cov())
    
            
app = App(app_ui, server, debug=True)
