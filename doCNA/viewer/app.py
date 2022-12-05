from shiny import *
from Plots import *

# Import modules for plot rendering
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as sig

chromlist = ['chr' + str (i) for i in range (1,23)]
chromdic = {}
for c in chromlist:
    chromdic[c] = c

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
                                       ui.input_slider ('k_max', "Max clonality score:",
                                                        value = 2, min = 0, max = 10),
                                       
                                       width = 2),
                      ui.panel_main(
                                    ui.navset_tab (
                                                   ui.nav("Genome-wide view",
                                                          #ui.row (ui.column(12, ui.output_plot ('model_legend'), )),
                                                                  #ui.tags.img(src = 'www/models_legend.png', alt = 'models legend', height = '10%', width = '100%'),),
                                                          ui.row(ui.column(12,          
                                                                 #ui.output_plot ('model_legend'),
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
                                                                                  ui.output_plot ('CNV_plot'),),
                                                                           ui.row(
                                                                                  ui.output_plot('scoring_dists')),
                                                                           )),
                                                         ),
                                                   
                                                   ui.nav("LOG", 
                                                          ui.row(ui.column(12,
                                                                 ui.input_file ('log_file', "Choose LOG file to screen:",
                                                                                multiple = False, accept = '.log'),
                                                                 #ui.output_text ('log_text'),
                                                                 
                                                                 #ui.tags.textarea(id = 'log', name = 'log_text', rows = '30', cols = '200',   )
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
                                                         ui.row(ui.output_table (id = 'CNVs'),)),
                                                   ui.nav("Solution test",
                                                          ui.layout_sidebar(ui.panel_sidebar(ui.h4 ("Optimize settings:"),
                                                                                             ui.input_slider ('min_cn', "Min cov (relative):",
                                                                                                              value = 0.4, min = 0.1, max = 1, step = 0.05),
                                                                                             ui.input_slider ('max_cn', "Max cov (relative):",
                                                                                                              value = 1.1, min = 0.5, max = 1.5, step = 0.05),
                                                                                             ui.row(ui.column(6,ui.input_numeric ("step", 'Step:', 0.1, min = 0.01, max = 1, step = 0.1)),
                                                                                                    ui.column(6,#ui.h6 ("# steps"),
                                                                                                                ui.input_text ("number_points", '# points')),
                                                                                                   ),    
                                                                                             ui.h6("Coverage range:"),
                                                                                             ui.output_text_verbatim ("coverage_range"),
                                                                                             #ui.input_checkbox_group ('cn4_models',
                                                                                             #       "Include cn4 models?",
                                                                                             #       {'yes' : 'Please do!'}),
                                                                                             ui.input_slider ('min_ai', "Min allelic imbalance:",
                                                                                                              value = 0.01, min = 0, max = 1),
                                                                                             ui.input_slider ('max_ai', "Max allelic imbalance:",
                                                                                                               value = 0.01, min = 0, max = 1),
                                                                                             ui.input_action_button ('opt',
                                                                                                                    "Optimize solution"),
                                                                                             width = 2),
                                                                            ui.panel_main(ui.h6("Solution check plot:"),
                                                                                          ui.output_plot('solution_plot_opt', "Solution plot"),
                                                                                          ui.h6("Total distance to solution at coverage:"),
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
                                                          ui.row(ui.output_plot ('data_plot'),
                                                                 ui.output_plot ('compare_plot')),
                                                          ui.row(ui.output_table (id = 'chrom_segments'))),
                                                    ui.nav("Report",
                                                           ui.row(ui.column(12,
                                                                            ui.output_plot('report_plot'),
                                                                            ui.output_table('report'))))
                                                   )
                                   )
          )
        ) 

def server(input, output, session):
    bed_full = reactive.Value(pd.DataFrame())
    bed = reactive.Value(pd.DataFrame())
    data = reactive.Value(pd.DataFrame())
    par = reactive.Value({})
    opt_solution = reactive.Value((np.array([]), np.array([])))
    chrom_sizes = reactive.Value(pd.Series())
    m0 = reactive.Value(np.nan)
    m0_opt = reactive.Value(np.nan)
    log_file = reactive.Value ([])
    bed_report = reactive.Value(pd.DataFrame())

    @output
    @render.text
    def solutions ():
        m, d = opt_solution ()
        ind = sig.argrelmin (d)[0]
        minims = []
        for i in ind:
            minims.append ((m[i], d[i])) 
        minims.sort (key = lambda x: x[1], reverse = False)
        
        return '\n'.join(['m = ' + str(m) + '  d = ' + str(d) for m,d in minims])
  
   
    @output
    @render.ui
    def dyn_log_ui():
        return ui.TagList (ui.tags.textarea ([l.strip() for l in log_file()], 
                                             cols = "250", rows = "50"))
    
    @output
    @render.text
    def coverage_range ():
        if np.isnan(m0()):
            text = 'n.a.'
        else:
            text = '{:.2f}'.format(m0()*input.min_cn()) + ' - ' + '{:.2f}'.format(m0()*input.max_cn())
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
        #print (file_input)
        if not file_input:
            return
        df = pd.read_csv (file_input[0]['datapath'], sep = '\t', header = None, 
                   names = ['chrom', 'start', 'end', 'ai', 'm', 'cn','model', 'd', 'model_score',
                            'k', 'k_score','dd', 'cyto', 'cent', 'status_d', 'status'])
        df['size'] = (df['end'] - df['start'])/1e6
        #print (df.head())
        data.set(pd.DataFrame())
        par.set({})
        bed_full.set(df)
        opt_solution.set ((np.array([]), np.array([])))
        log_file.set([])
        bed_report.set(pd.DataFrame())
            
    @reactive.Effect
    @reactive.Calc
    def _():
        tmp = bed_full()
        if len(tmp) > 0:
            chrom_sizes.set(tmp.groupby(by = 'chrom').agg({'end' : 'max'})['end'])
            tmp['filt'] = (tmp['cent'] <= input.cent_thr()) & (tmp['size'] >= input.size_thr())
            bed_full.set(tmp)
            bed.set (tmp.loc[tmp.filt])
    

    @reactive.Effect
    @reactive.Calc
    def _():
        bf = bed_full()    
        b = bed()
        if (len(bf) != 0) & (len(b) != 0):
            records = []
            for chrom, segments in bed_full.groupby (by = 'chrom'):
                data = segments.sort_values (by = 'start')
                current = data.iloc[0,]
                if len (data) > 1:
                    i = 1 
                    while i < len(data):
                        next = data.iloc[i,]
                        if (current['status'] == next)


                        else: #add to records


                    
                else:
                    records = [(chrom, current['start'], current['end'], current['m'], current['cn'],
                               current['model'] if current['status'] != 'norm' else 'AB',
                               current['k'] if current['status'] != 'norm' else 0,
                               current['cyto'])]

            bed_report.set(pd.DataFrame.from_records (records,
                                                      columns = ['chrom', 'start', 'end', 'm', 'cn','model', 'k', 'cyto']))
    


    @reactive.Effect
    @reactive.event(input.par_file)
    def _():
        file_input = input.par_file()
        #print (file_input)
        if not file_input:
            return
        pard = {}
        with open(file_input[0]['datapath'],'r') as f:
            for line in f.readlines():
                try:
                    key, value = line.split('\t')
                    pard[key] = (float(value),)
                except:
                    try:
                        value0, value1 = value.split(' ')
                        pard[key] = (float(value0),float(value1))
                    except:
                        print ('Line: ' + line + 'not parsed.')
        par.set(pard)
        #print (pard)
        opt_solution.set ((np.array([]), np.array([])))
        m0.set(float(pard['m0'][0]))
        m0_opt.set(float(pard['m0'][0]))
        
        
    @reactive.Effect
    @reactive.event(input.data_file)
    def _():
        file_input = input.data_file()
        #print (file_input)
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
                          max_k_score = input.k_max(),
                          model_thr = input.model_thr())
            
            for model in colorsCN.keys():
                axs[0].plot ((),(), lw = 10, color = colorsCN[model], label = model)
            axs[0].plot ((),(), lw = 10, color = 'yellow', label = 'complex')
            axs[0].plot ((),(), lw = 10, color = 'red', label = 'fail')
            axs[0].legend (bbox_to_anchor = (0.5, 2), ncol = len(model_presets)+2,
                           loc = 'upper center')
            return fig
    
    @output
    @render.plot (alt = "Scoring view")
    def CNV_plot ():
        bed_data = bed()
        par_d = par()
        #print (par_d)
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, ax = plt.subplots (1, 1, figsize = (6,3))
            #leopard_plot (bed_data, par_d, ax, highlight = input.chroms_selected())
            leopard_plot (bed_data.loc[bed_data['model'] != 'A(AB)B'], 
                          (par_d['A_i'][0], par_d['C_i'][0], par_d['C_i'][0]-par_d['up_i'][0]),
                          ax, highlight = input.chroms_selected(),
                          color_norm = 'black', color_hit = 'darkred')
            leopard_plot (bed_data.loc[bed_data['model'] == 'A(AB)B'], 
                          (np.nan, np.nan, np.nan),
                          ax, highlight = input.chroms_selected(), 
                          color_norm = 'gray', color_hit = 'darkorange', alpha = 0.3)
            ax.set_xlim ((np.log10(0.95*input.size_thr()), 
                          np.log10(1.05*bed_data['size'].max())))
            k_pos = bed_data.loc[bed_data['k'] > 0, 'k'].values
            ax.set_ylim ((np.log10(k_pos.min()), 0.1))  # type: ignore
            return fig
    
    
    @output
    @render.plot (alt = "Solution view")
    def solution_plot ():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, ax = plt.subplots (1, 1, figsize = (6,6))
            check_solution_plot_opt (bed_data, par_d, ax, 
                                          highlight = input.chroms_selected())
            k = np.linspace (0,1,100)
            m0 = par_d['m0']
            for model in colorsCN.keys():
                ax.plot (model_presets[model].m(k, m0), model_presets[model].ai(k, m0),  
                         lw = 2, linestyle = '-', color = colorsCN[model], alpha  = 0.6)
            
            ax.set_xlim (0.9*bed_data.m.min(), 1.1*bed_data.m.max())
            ax.set_ylim ((max(-0.02, -0.02*bed_data.ai.max()), bed_data.ai.max()*1.1))
            
            return fig
    
    @output
    @render.plot (alt = 'Scoring distribution')
    def scoring_dists ():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, axs = plt.subplots (1, 2, figsize = (6,3), sharey = True)
            tmp = bed_data.loc[bed_data.model != 'A(AB)B']
            if len(tmp) > 0:
                plot_cdf (tmp['dd'].values,
                          ax = axs[0], par = ((par_d['m_i'],),(par_d['s_i'],), (sum(tmp.status == 'norm')/len(tmp),)))
                axs[0].set_title ('Imbalanced')
                axs[0].set_xlabel ('Distance to usual')
                axs[0].set_ylabel ('cdf')
            tmp = bed_data.loc[bed_data.model == 'A(AB)B']
            if len(tmp) > 0:
                if len (par_d['a_b']) == 2:
                    neg_tmp = tmp.loc[tmp.k < 0]
                    fnnorm = sum(neg_tmp.status_d == 'norm')/len (neg_tmp)

                    pos_tmp = tmp.loc[tmp.k >= 0]
                    fpnorm = sum(pos_tmp.status_d == 'norm')/len (pos_tmp)

                    plot_cdf (tmp['k'].values, ax = axs[1], 
                              par = (par_d['m_b'],par_d['s_b'], (par_d['a_b'][0]*fnnorm, par_d['a_b'][1]*fpnorm)),
                              a0 = len(tmp.loc[(tmp.k < 0)& (tmp.status_d != 'norm')])/len(tmp) )

                elif len(par_d['a_b']) == 1:
                    plot_cdf (tmp['k'].values, ax = axs[0], 
                              par = ((par_d['m_b'],),(par_d['s_b'],), (sum(tmp.status != 'norm')/len(tmp),)))
                 
                axs[1].set_title ('Balanced')
                axs[1].set_xlabel ('Sign clonality')
                
            return fig
    
    
    @output
    @render.plot (alt = "Solution view")
    def solution_plot_opt ():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
           
            fig, ax = plt.subplots (1, 1, figsize = (6,6))
            check_solution_plot_opt (bed_data, par_d, ax, 
                                          highlight = [])
            k = np.linspace (0,1,100)
            for model in model_presets.keys():
                ax.plot (model_presets[model].m(k, m0()), model_presets[model].ai(k, m0()), 
                         lw = 2, linestyle = '-', color = colorsCN[model], alpha = 0.6)
                
            
            ax.legend (bbox_to_anchor = (1,0))
               
            return fig
    
    @output
    @render.plot (alt = 'Total distance to the model')
    def opt_plot ():
        #print (opt_solution())
        if len(opt_solution()[0]):
            opt = opt_solution()
            fig, ax = plt.subplots(1,1, figsize = (6,6))
            ax.plot (opt_solution()[0], opt_solution()[1], 'r-')
            ax.set_xlabel ('Covearage')
            ax.set_ylabel ('Total distance to the solution')
            for model in colorsCN.keys():
                ax.plot ((),(), lw = 2, color = colorsCN[model], label = model)
            ax.legend (bbox_to_anchor = (1,1))
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
    def CNVs ():
        bed_data = bed()
        if (len(bed_data) != 0):
            if input.sort_CNV_by() == 'score':
                tmp_bed = bed_data.loc[bed_data[input.corrected()] != 'norm'].sort_values(by = 'k_score', ascending = False)
            else:
                tmp_bed = bed_data.loc[bed_data[input.corrected()] != 'norm']
            return tmp_bed

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
            verification_plot_CNV (data_chrom, CNV_bed, ax, par())
            return fig
        
        
    @reactive.Effect
    @reactive.event (input.opt)
    def _():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            #models = models_dic()
                        
            ms = np.arange (m0()*input.min_cn(), m0()*input.max_cn(), input.step())
            tmp = bed_data.loc[(~bed_data['model'].isna())&((bed_data['ai'] <= input.min_ai())|(bed_data['ai'] >= input.max_ai()))]
            
            with ui.Progress (min = ms[0], max = ms[-1]) as p:
                p.set(message = 'Optimizing solution...', )
                dts = []
                sts = []
                for m in ms:
                    p.set(m, message = 'Calculating')
                    dt = 0
                    st = 0
                    for _, b in tmp.iterrows():
                        dt += min([calculate_distance(model, b['m'], b['ai'], m) for model in model_presets.values()])*b['size']  
                        st += b['size']
                    dts.append(dt)
                    sts.append(st)
                dtsa = np.array(dts)/np.array(sts)
                opt_solution.set((ms,dtsa))
                
            m0_opt.set (ms[np.where(dtsa == dtsa.min())[0][0]])
            
    @reactive.Effect
    @reactive.event (input.max_ai)
    def _():
        if input.min_ai() > input.max_ai():
            ui.update_slider ('min_ai', value = input.max_ai())
            
    @reactive.Effect
    @reactive.event (input.min_ai)
    def _():
        if input.max_ai() < input.min_ai():
            ui.update_slider ('max_ai', value = input.min_ai())
        
    @reactive.Effect
    @reactive.event (input.m0_cov)
    def _():
        if len(opt_solution()[0]):
            m0.set(input.m0_cov())
    
    
    @output
    @render.text
    def m_diploid ():
        par_d = par()
        opt = opt_solution()
        if (len(par_d.keys()) != 0) & (len(opt[0]) > 0):
            min_dist_at = opt[0][np.where(opt[1] == opt[1].min())[0][0]]
            text = ["Experimental feature!!",
                    "Found m0: " + str(par['m0']),
                    "Optimized m0: " + str(min_dist_at),
                    f"The relative difference: {np.abs(min_dist_at - par_d['m0'])/par_d['m0']}"]
            return '\n'.join(text)
        else:
            return
        
        
app = App(app_ui, server, debug=True)
