import numpy as np
import scipy.stats as sts
import scipy.optimize as opt
from collections import namedtuple

from doCNA import Testing
from doCNA import Distribution
from doCNA import Segment
from doCNA import Run
from doCNA.Report import Report

N_SYMBOL = 'N'
E_SYMBOL = 'E'
U_SYMBOL = 'U'
DEFAULT_N_THRESHOLD = 10
DEFAULT_E_THRESHOLD = 3
N_STR_LEN_THR = 100

Run_treshold =  namedtuple('Run_treshold', [N_SYMBOL, E_SYMBOL])


class Chromosome:
    """Class to contain data, and find runs."""
    def __init__ (self, name, data, config, logger, genome_medians, CB):
        self.name = name
        self.data = data
        self.config = config
        self.logger = logger.getChild(f'{self.__class__.__name__}-{self.name}')
        self.genome_medians = genome_medians
        
        self.CB = CB #.loc[CB['gieStain'] != 'acen']
        self.cent = (CB.loc[(CB['gieStain'] == 'acen'),'chromStart'].min(),
                     CB.loc[(CB['gieStain'] == 'acen'),'chromEnd'].max())
         
        self.Eruns = []
        self.Uruns = []
        self.Nruns = []
        self.logger.debug (f"Object chromosome {name} created.")
        
    def markE_onHE(self, he_parameters, threshold = 13.8):
        zv = (self.data['vaf'] - he_parameters['vaf']) / np.sqrt (0.25/(he_parameters['cov']))
        zc = (self.data['cov'] - he_parameters['cov']) / np.sqrt (he_parameters['b']*he_parameters['cov'])
        z = zv**2+zc**2
 
        self.data['symbol'] = N_SYMBOL       
        indexes = self.data.loc[z < threshold, :].index.values.tolist()
        self.data.loc[indexes, 'symbol'] = E_SYMBOL
        self.logger.debug (f"""Chromosome {self.name} marked based on parameters
v = {he_parameters['vaf']}, c = {he_parameters['cov']}.""")
        
    def mark_on_full_model (self, m):
        self.get_fragments (n = int(self.config['Segmenting']['No_SNPs']))
        self.get_vaf_shift_full ()
        
        indexes, merged_string = Run.merge_symbols (self.dv_dist.string, outliers_threshold = 4)
        for r in indexes:
            start = self.windows_positions[r[0]][0]
            end = self.windows_positions [r[1]][1]
            chi2, vaf, fb = Testing.VAF_test (self.data.loc[(self.data['position'] >= start)&(self.data['position'] <= end),], m)
            
            #test chi2
            outlier = (chi2 > float(self.config['VAF']['chi2_high'])) | (chi2 == 0)
                      
            self.logger.info (f'Marking on full model {self.name}:{start}-{end} chi2 = {chi2}')            
            if outlier:
                self.data.loc[(self.data['position'] >= start)&\
                                    (self.data['position'] <= end), 'symbol'] = U_SYMBOL
                self.Uruns.append ((start, end))
    
    def get_fragments (self, n = 1000):
        tmp = self.data
        N = np.max((int(np.floor(len(tmp)/n)),1))
        
        indexes = np.linspace (0,len (tmp), 2*N+2, dtype = int)
        self.indexes = indexes
        self.windows = []
        self.windows_positions = []

        for i in np.arange (2*N):
            tmpi = tmp.iloc[indexes[i]:indexes[i+2], ]
            self.windows.append(tmpi)
            self.windows_positions.append ((tmpi['position'].min(), tmpi['position'].max()))

    def get_vaf_shift_full (self, z_thr = 2.5):
        
        def vaf_cdf (v, dv, a, lerr, f, vaf, b):
            cnai = vaf_cnai (v, dv, f, vaf, b, cov)
            cnHO = vaf_HO (v, lerr)
            return a*cnHO + (1 - a)*cnai 
        
        v0 = self.genome_medians['HE']['vaf']
        b = self.genome_medians['HE']['b']
        cov = self.genome_medians['HE']['cov']
        dvs = []
        v0s = []        
        
        for window in self.windows:
            vaf = window['vaf'].values
            v, c = np.unique(vaf, return_counts = True)

            cnor = np.cumsum(c)/np.sum(c)
            ones0 = c[v >= (cov-1)/cov].sum()
            try:
                f0 = c[v < v0].sum()/(c.sum() - ones0) 
            
                dv0 = v0 - np.median (v[v < v0])
              
                popt, pcov = opt.curve_fit (vaf_cdf, v, cnor, p0 = [dv0, ones0/c.sum(), 2, f0, 0.5, b], 
                                            bounds = ((0,   0,   1, 0, 0.45, 1),
                                                      (0.5, 0.95, 5, 1, 0.55, 10)))
                dvs.append (popt[0])
                v0s.append (popt[-1])
            except (RuntimeError, ValueError):
                dvs.append (0)
                
        dva = np.array (dvs)
        median = np.median(dva[dva > 0])
        dva[dva == 0] = median
        self.dv = dva
        self.v0 = np.array(v0s)
        #p_thr is lower that for sensitive as full is more noisy, but less nosy :D 
        self.dv_dist = Distribution.Distribution (self.dv,
                                                  p_thr = 0.1, thr_z = z_thr)

    def find_runs (self):
        """Method to generate runs. Runs segment themselves."""
        self.find_Nruns ()
                
        self.runs = []
        #print (self.name)
        #print ("N: ", len (self.Nruns))
        #print ("U: ", len (self.Uruns))

        for nr in self.Nruns:
            self.runs.append(Run.Run (self.data.loc[(self.data['position'] >= nr[0])&(self.data['position'] <= nr[1])],
                                      symbol = 'N',
                                      logger = self.logger,
                                      genome_medians = self.genome_medians))
        
        for ur in self.Uruns:
            self.runs.append(Run.Run (self.data.loc[(self.data['position'] >= ur[0])&(self.data['position'] <= ur[1])],
                                      symbol = 'U',
                                      logger = self.logger,
                                      genome_medians = self.genome_medians))
        
        data_view = self.data.loc[self.data['symbol'] == 'E']
        if len(data_view) > 0:
            self.runs.append (Run.Run (self.data.loc[self.data['symbol'] == 'E'],
                                       symbol = 'E',
                                       logger = self.logger,
                                       genome_medians = self.genome_medians))
        
    def generate_segments (self):
        """Method to generate genomic segments of same CNA status, based on Runs."""
        self.segments = []
        unr = []
        for nr in self.Nruns:
            unr.append ((nr[0], nr[1]))
        for ur in self.Uruns:
            unr.append ((ur[0], ur[1]))
        unr.sort (key = lambda x: x[0] )
        #print (self.name, unr)
        startsandends = []
        
        #print (self.name, 'seg')

        for run in self.runs:
            best_solution = run.solutions[0]       
            if run.symbol != 'E':
                startsandends+= best_solution.positions
            else:
                for start, end in best_solution.positions:
                    current_start = start
                    for rstart, rend in unr:
                        if (rstart < end) & (rend > start):
                            startsandends.append((current_start,rstart))
                            current_start = rend
                    if current_start < end:
                        startsandends.append((current_start,end))
        
        #print (self.name, startsandends)

        for start, end in startsandends:
            data_view = self.data.loc[(self.data['position'] >= start) &\
                                      (self.data['position'] <= end)]
            if len(data_view) == 0:
                self.logger.error(f"Wrong segment {start}-{end} in {run.name})")
            else:
                centromere_fraction = min((end - self.cent[0]),(self.cent[1]-start))/(end - start)
                cytobands = self.CB[(self.CB['chromStart'] < end)&(self.CB['chromEnd'] > start)].sort_values (by = 'chromStart')['name'].values
                        
                if len(cytobands) > 1:
                    cytobands_str = cytobands[0] + '-' + cytobands[-1]
                else:
                    cytobands_str = cytobands[0]
                            
                self.segments.append (Segment.Segment (data = data_view, 
                                                       config = self.config, 
                                                       logger = self.logger, 
                                                       genome_medians = self.genome_medians,
                                                       segmentation_score = best_solution.p_norm,
                                                       segmentation_symbol = run.symbol,
                                                       centromere_fraction = 0 if (centromere_fraction < 0) | (centromere_fraction > 1) else centromere_fraction,
                                                       cytobands = cytobands_str))
    
    def find_Nruns (self):
        symbol_list = self.data.loc[(self.data['vaf'] < 1) & (self.data['symbol'] != U_SYMBOL), 'symbol'].tolist()
        if len (symbol_list) >= N_STR_LEN_THR:    
            self.Nruns_indexes, self.Nruns_threshold = analyze_string_N (symbol_list, N = N_SYMBOL, E = E_SYMBOL)
        else:
            self.Nruns_indexes = []
            self.Nruns_threshold = []
        
        for run in self.Nruns_indexes:
            tmp = self.data.loc[self.data['vaf'] < 1,].iloc[run[0]:run[1],:].position.agg((min,max))
            self.data.loc[(self.data.position >= tmp['min'])&(self.data.position <= tmp['max']), 'symbol'] = N_SYMBOL
            self.Nruns.append ((tmp['min'], tmp['max']))
    
    def report (self, report_type = 'bed'):
        return Report(report_type).chromosome_report(self.segments, self.runs)

def vaf_cnai (v, dv, a, vaf,b, cov):
    s = np.sqrt((vaf - dv)*(vaf + dv)/(b*cov))
    return a*sts.norm.cdf (v, vaf - dv, s) + (1-a)*sts.norm.cdf (v, vaf + dv, s)

def vaf_HO (v, lerr):
    err = 10**lerr
    return np.exp ((v-1)*err)

def analyze_string_N (symbol_list, N = 'N', E = 'E'):
    """Finds runs of Ns"""                    
    
    values, counts = Run.rle_encode (''.join(symbol_list))
    threshold = find_runs_thr (values, counts, N = N, E = E)
    runsi = get_N_runs_indexes (values, counts, N = N, E = E, threshold = threshold)
        
    for si,ei in runsi:
        for i in range(si,ei):
            symbol_list[i] = N          
    
    values, counts = Run.rle_encode (symbol_list)
    threshold = find_runs_thr (values, counts, N = N, E = E)
    runsi = get_N_runs_indexes (values, counts, threshold = threshold, N = N, E = E)
        
    return runsi, threshold 
    
def find_runs_thr (values, counts, N = 'N', E = 'E'):
    assert len (values) == len (counts), 'Wrong input!! Go away!'
    hist = np.unique(counts[values == N], return_counts = True)
    x = np.log10 (hist[0])
    y = np.log10 (hist[1])
    
    try:
        ind = np.arange (0, min (4, len(x)))
        popt, _ = opt.curve_fit (lin,  x[ind], y[ind], p0 = [-1,1])
        xt = np.log10 (np.arange(1, hist[0].max()))
        N_thr = 10**(xt[lin(xt, *popt) > np.log10(1)].max())
    except RuntimeError:
        N_thr = DEFAULT_N_THRESHOLD
    
    hist = np.unique(counts[values != N], return_counts = True)
    x = hist[0]
    y = np.log10 (hist[1])
    
    try:
        popt, _ = opt.curve_fit (lin, x[:4], y[:4], p0 = [-1,1])
        if popt[1] < 0:
            xt = np.arange(1, hist[0].max())
            E_thr = xt[lin(xt, *popt) > 0].max()
        else:
            popt, _ = opt.curve_fit (lin, x[:3], y[:3], p0 = [-1,1])
            if popt[1] < 0:
                xt = np.arange(1, hist[0].max())
                E_thr = xt[lin(xt, *popt) > 0].max()
            else:
                E_thr = DEFAULT_E_THRESHOLD
    except RuntimeError:
        E_thr = DEFAULT_E_THRESHOLD
    
    return Run_treshold (N = N_thr, E = E_thr)

def lin (x,a,b):
    return a*x+b

def get_N_runs_indexes (values, counts, threshold, N = N_SYMBOL, E = E_SYMBOL):
    bed = []
    Nindexes = np.where ((values == N)&(counts >= threshold.N))[0]
    for i in Nindexes:
        if any ([(i>=s)&(i <= e) for s,e in bed]):
            continue
        
        s = i
        e = i
        
        if e < len(counts)-2:
            inext = e + 1        
            while counts[inext] < threshold.E:
                if e > len(counts)-4:
                    e = len(counts)-1
                    break
                else:
                    e += 2
                    inext = e + 1
                
        if s > 0:
            iprev = s - 1
            while counts[iprev] < threshold.E:
                if s < 3:
                    s = 0
                    break
                else:
                    s -= 2
                    iprev = s - 1
        bed.append ((s,e))    
    return [[counts[:i].sum(),counts[:j+1].sum()] for i, j in bed] 
