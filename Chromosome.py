import numpy as np
import scipy.stats as sts
import scipy.optimize as opt
from collections import namedtuple

from doCNA import Testing
from doCNA import Distribution
from doCNA import Segment
from doCNA import Run

N_SYMBOL = 'N'
E_SYMBOL = 'E'
U_SYMBOL = 'U'
DEFAULT_N_THRESHOLD = 10
DEFAULT_E_THRESHOLD = 3

Run_treshold =  namedtuple('Run_treshold', [N_SYMBOL, E_SYMBOL])


class Chromosome:
    """Class to contain data, and find runs."""
    def __init__ (self, name, data, config, logger, genome_medians):
        self.name = name
        self.data = data
        self.config = config
        self.logger = logger.getChild(f'{self.__class__.__name__}-{self.name}')
        self.genome_medians = genome_medians
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
        
    def mark_on_full_model (self):
        self.get_fragments (exclude_symbols = [], n = int(self.config['Segment']['No_SNPs']))
        self.get_vaf_shift_full ()
        self.Uruns = []
        
        indexes, merged_string = Run.merge_symbols (self.dv_dist.string, outliers_threshold = 4)

        for r in indexes:
            start = self.positions [r[0]][0]
            end = self.positions [r[1]][1]
            chi2, vaf = Testing.VAF_test (self.data.loc[(self.data['symbol'] != N_SYMBOL)&\
                                                        (self.data['position'] >= start)&\
                                                        (self.data['position'] <= end)])
            
            #test chi2
            outlier = (chi2 > float(self.config['Tests']['VAF']['chi2_high'])) | (chi2 == 0) |\
                      (vaf < self.genome_medians['VAF']['vaf_thr'][0])
                      
            self.logger.info (f'Marking on full model {self.name}:{start}-{end} chi2 = {chi2}')            
            if outlier:
                self.data.loc[(self.data['position'] >= start)&\
                                    (self.data['position'] <= end), 'symbol'] = U_SYMBOL
                self.Uruns.append ((start, end))
    
        
    #or whatever name
    def find_runs (self):
        """Method to generate runs. Runs segment themselves."""
        self.find_Nruns ()
                
        self.runs = []
        for nr in self.Nruns:
            self.runs.append(Run.Run (self.data.loc[(self.data['position'] >= nr[0])&(self.data['position'] <= nr[1])],
                                  symbol = 'N',
                                  logger = self.logger,
                                  genome_medians = self.genome_medians))
        
        for ur in self.Uruns:
            self.runs.append(Run.Run (self.data.loc[(self.data['position'] >= nr[0])&(self.data['position'] <= nr[1])],
                                  symbol = 'U',
                                  logger = self.logger,
                                  genome_medians = self.genome_medians))
            
        #Eruns needs to be found first
        start = self.data['position'].min ()
        end = self.data['position'].max ()
        
        self.runs.append (Run.Run (self.data.loc[self.data['symbol'] == 'E'],
                               symbol = 'E',
                               logger = self.logger,
                               genome_medians = self.genome_medians))
        
    def generate_segments (self):
        """Method to generate genomic segments of same CNA status, based on Runs."""
        pass    
      
    def find_Nruns (self):
        symbol_list = self.data.loc[(self.data['vaf'] < 1) & (self.data['symbol'] != U_SYMBOL), 'symbol'].tolist()
        self.Nruns_indexes, self.Nruns_threshold = analyze_string_N (symbol_list, N = N_SYMBOL, E = E_SYMBOL)
        
        for run in self.Nruns_indexes:
            tmp = self.data.loc[self.data['vaf'] < 1,].iloc[run[0]:run[1]].position.agg((min,max))
            self.data.loc[(self.data.position >= tmp['min'])&(self.data.position <= tmp['max']), 'symbol'] = N_SYMBOL
            self.Nruns.append ((tmp['min'], tmp['max']))
    
    def report (self, type = 'bed'):
        return "blah"
    
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
    #[(counts[:si].sum(), counts[:ei+1].sum()) for s, e in runsi], threshold 
    
def find_runs_thr (values, counts, N = 'N', E = 'E'):
    assert len (values) == len (counts), 'Wrong input!! Go away!'
    hist = np.unique(counts[values == N], return_counts = True)
    # print (hist)
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
    #counts[(counts > threshold.N)& (values == N_SYMBOL)]
    #print (Nindexes)
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