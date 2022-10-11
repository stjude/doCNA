from shiny import *
from Plots import *

# Import modules for plot rendering
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as sig

__version__ = '0.1.3.1'

chromlist = ['chr' + str (i) for i in range (1,23)]
chromdic = {}
for c in chromlist:
    chromdic[c] = c

app_ui = ui.page_fluid(
    ui.h2 ({"style" : "text-align: center;"}, "doCNA results viewer. v. " + __version__),
    ui.h4 ({"style" : "text-align: center;"}, "for doCNA >= 0.8.3.1"),
    
    ui.layout_sidebar(ui.panel_sidebar(ui.h4 ("Segments filtering:"),
                                       ui.input_slider ('cent_thr', "Centromere fraction threshold",
                                                        value = 0.3, min = 0, max = 1),
                                       ui.input_slider ('size_thr', "Min segment size (in Mb)",
                                                        value = 5, min = 0, max = 10),
                                       ui.h4 ("Display settings:"),
                                       ui.input_slider ('model_thr', "Model score threshold",
                                                        value = 5, min = 0, max = 10),
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
                                                                                multiple = False),
                                                                 ui.output_plot ('genome_plot'),)),
                                                          ui.row(ui.column(12,
                                                                           ui.input_checkbox_group ('chroms_selected',
                                                                                                    "Select chromosomes to highlight",
                                                                                                    chromdic, inline = True)),),   
                                                          ui.row(ui.column(6,
                                                                           ui.input_file ('par_file', "Choose PAR file to upload:",
                                                                                          multiple = False),
                                                                           ui.output_plot ('CNV_plot'),
                                                                           ),
                                                                 ui.column(6,
                                                                           ui.row(ui.output_plot ('solution_plot'),),)),
                                                         ),
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
                                                                                          ui.input_slider ('m0_cov', "Diploid coverage",
                                                                                                               value = 0.5, min = 0, max = 1, width = '200%'),
                                                                                          ui.output_text_verbatim ('solutions')))),
                                                                                          
                                                   ui.nav("Chromosome view",
                                                          ui.row(ui.column(12,
                                                                           ui.input_file ('data_file', "Choose data file to upload:",
                                                                                           multiple = False),
                                                                           ui.input_radio_buttons ('chrom_view', "Choose chromosome to inspect",                                                                                    
                                                                                                   chromdic, inline = True),)),
                                                          ui.row(ui.output_plot ('data_plot'),
                                                                 ui.output_plot ('compare_plot')),
                                                          ui.row(ui.output_table (id = 'chrom_segments'),))
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
    @reactive.event(input.bed_file)
    def _():
        file_input = input.bed_file()
        #print (file_input)
        if not file_input:
            return
        df = pd.read_csv (file_input[0]['datapath'], sep = '\t', header = None, 
                   names = ['chrom', 'start', 'end', 'ai', 'm', 'cn','model', 'd', 'model_score',
                            'k', 'k_score','dd', 'cyto', 'cent', 'status'])
        df['size'] = (df['end'] - df['start'])/1e6
        #print (df.head())
        data.set(pd.DataFrame())
        par.set({})
        bed_full.set(df)
        opt_solution.set ((np.array([]), np.array([])))
        #models_dic().update (model_presets)
            
    @reactive.Effect
    @reactive.Calc
    def _():
        tmp = bed_full()
        if len(tmp) > 0:
            chrom_sizes.set(tmp.groupby(by = 'chrom').agg({'end' : 'max'})['end'])
            bed.set (tmp.loc[(tmp['cent'] <= input.cent_thr()) & (tmp['size'] >= input.size_thr())])
            print (tmp.loc[(tmp['cent'] <= input.cent_thr()) & (tmp['size'] >= input.size_thr())].head()) 
    
        
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
                key, value = line.split('\t')
                pard[key] = float(value)
        par.set(pard)
        opt_solution.set ((np.array([]), np.array([])))
        m0.set(float(pard['m']))
        m0_opt.set(float(pard['m']))
        
        
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
            axs[0].legend (bbox_to_anchor = (0.9, 1.7), ncol = len(model_presets)+2)
            return fig
    
    @output
    @render.plot (alt = "Scoring view")
    def CNV_plot ():
        bed_data = bed()
        par_d = par()
        #print (par_d)
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            fig, ax = plt.subplots (1, 1, figsize = (6,6))
            #leopard_plot (bed_data, par_d, ax, highlight = input.chroms_selected())
            leopard_plot (bed_data.loc[bed_data['model'] != 'A(AB)B'], 
                          (par_d['a_i'], par_d['b_i'], par_d['bt_i']),
                          ax, highlight = input.chroms_selected(),
                          color_norm = 'black', color_hit = 'darkred')
            leopard_plot (bed_data.loc[bed_data['model'] == 'A(AB)B'], 
                          (par_d['a_b'], par_d['b_b'], par_d['bt_b']),
                          ax, highlight = input.chroms_selected(), 
                          color_norm = 'gray', color_hit = 'darkorange')
            ax.set_xlim ((np.log10(0.95*input.size_thr()), 
                          np.log10(1.05*bed_data['size'].max())))
            ax.set_ylim ((np.log10(bed_data['k'].min()), 0.1))  # type: ignore
            return fig
    
    #@output
    #@render.plot (alt = "Models")
    #def model_legend():
    #    fig, ax = plt.subplots (1,1, figsize = (12,0.5))
    #    for model in colorsCN.keys():
    #            ax.plot ((),(), lw = 10, color = colorsCN[model], label = model)
    #    ax.plot ((),(), lw = 10, color = 'yellow', label = 'complex')
    #    ax.plot ((),(), lw = 10, color = 'red', label = 'fail')
    #    ax.legend (loc = "center right", ncol = len(model_presets)+2)
    #    ax.axis('off')
        
    #    return fig
            
    @output
    @render.plot (alt = "Solution view")
    def solution_plot ():
        bed_data = bed()
        par_d = par()
        if (len(bed_data) != 0) & (len(par_d.keys()) != 0):
            #models = models_dic()
            fig, ax = plt.subplots (1, 1, figsize = (6,6))
            check_solution_plot_opt (bed_data, par_d, ax, 
                                          highlight = input.chroms_selected())
            k = np.linspace (0,1,100)
            m0 = par_d['m']
            for model in colorsCN.keys():
                ax.plot (model_presets[model].m(k, m0), model_presets[model].ai(k, m0),  
                         lw = 2, linestyle = '-', color = colorsCN[model], alpha  = 0.6)
            
            ax.set_xlim (0.9*bed_data.m.min(), 1.1*bed_data.m.max())
            ax.set_ylim ((max(-0.02, -0.02*bed_data.ai.max()), bed_data.ai.max()*1.1))
            
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
            bed_CNV = bed_data.loc[(bed_data.chrom == input.chrom_view())&(bed_data['status'] == 'CNV')]
            data_chrom = data_df.loc[data_df.chrom == input.chrom_view()]
            fig, ax = plt.subplots (figsize = (6,8))
            verification_plot_CNV (data_chrom, bed_CNV, ax)
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
                    "Found m0: " + str(par()['m']),
                    "Optimized m0: " + str(min_dist_at),
                    f"The relative difference: {np.abs(min_dist_at - par_d['m'])/par_d['m']}"]
            return '\n'.join(text)
        else:
            return
        
        
app = App(app_ui, server, debug=True)