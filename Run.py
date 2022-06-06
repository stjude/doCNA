#the idea here is that after initial segmentation, this stratify run into segments
import scipy.stats as sts
import numpy as np
import pandas as pd
import scipy.optimize as opt

from collections import namedtuple
import warnings as warn
from numpy.random import default_rng

from doCNA import Segment
from doCNA import Distribution
from doCNA import Testing
from doCNA import Chromosome

WINDOWS_THRESHOLD = 9
SNPS_IN_WINDOW = 1000
WINDOWS_TO_TEST_THRESHOLD = 20
UNIFORMITY_THRESHOLD = 1e-5

#Solution keep result of Run segmenting
#Segments are based on /best/ Solution
Solution = namedtuple ('Solution', ['chi2', 'chi2_noO', 'positions', 'p_norm', 'segments', 'merged_segments'])

class Run:
    """Class to segment run of E/N/U"""
    def __init__(self, data, symbol, logger, genome_medians) -> None:
        self.data = data
        self.symbol = symbol
        self.genome_medians = genome_medians
        
        self.chrom = data['chrom'].values[0]
        self.start = data['position'].min()
        self.end = data['position'].max()
        
        self.name = self.chrom+ ':'+str(self.start)+'-'+str(self.end)
        self.logger = logger.getChild(f'{self.__class__.__name__}-{self.name}')
        self.logger.debug ("Object created")
     
        self.analyze ()
        self.logger.info (f"Run analyzed, {len(self.solutions)} solutions")
                
    def analyze (self):
        self.get_windows (n = SNPS_IN_WINDOW)
        self.logger.debug (f'Run divided into {len(self.windows)}')
        if len(self.windows) >= WINDOWS_THRESHOLD:
            self.get_ai ()
            self.get_coverage ()
            self.solve_windows ()
        else:
            self.logger.info (f'Run {self.name} is to short to segment.')
            self.solutions = [Solution (chi2 = np.nan,
                                        chi2_noO = np.nan,
                                        positions = [(self.windows_positions[0][0], self.windows_positions[-1][1])],
                                        p_norm = [(np.nan, np.nan, np.nan)],
                                        segments = '',
                                        merged_segments = '')]
            self.logger.info (f'One solution devised for unsegmented run.')
    
    def get_windows (self, n = 1000):
        tmp = self.data.loc[[s == self.symbol for s in self.data['symbol']], ]
        N = np.max((int(np.floor(len(tmp)/n)),1))
        
        indexes = np.linspace (0,len (tmp), 2*N+2, dtype = int)

        self.windows = []
        self.windows_positions = []

        for i in np.arange (2*N):
            tmpi = tmp.iloc[indexes[i]:indexes[i+2], ]
            self.windows.append(tmpi)
            self.windows_positions.append ((tmpi['position'].min(), tmpi['position'].max()))
        
    def get_ai (self):
        if self.symbol == Chromosome.E_SYMBOL:
            self.get_ai_sensitive()
        else:
            self.get_ai_full()
        
    def get_ai_sensitive (self, zero_thr = 0.01, cov_mult = 1.0, p_thr = 0.5, z_thr = 1.5):
        tmpf = 1
        s0 = np.sqrt (0.25/self.genome_medians['COV']['m'])
        
        vafs = []
        for window in self.windows:
            vafs.append(np.sort(window['vaf'].values))
        
        while tmpf > zero_thr:
            dvl = []
            v0l = []
            s = s0/np.sqrt(cov_mult)
            def make_two_gauss (v, dv, v0):
                a = 0.5
                return a*sts.norm.cdf (v, v0 - dv, s) + (1-a)*sts.norm.cdf (v, v0 + dv, s)
            
            for vaf in vafs:
                
                popt, pcov = opt.curve_fit (make_two_gauss, np.sort(vaf), np.linspace (0,1, len(vaf)),
                                            p0 = [0.05, 0.5],
                                            bounds = ((0, 0.4),(0.5, 0.6)))
                
                dvl.append (popt[0])
                v0l.append (popt[1])
        
            tmpf = sum ([d < zero_thr for d in dvl])/len (dvl)
            cov_mult += 0.01
            
        self.dv = np.array(dvl)
        self.v0 = np.array (v0l)
        self.dv_dist = Distribution.Distribution (self.dv, p_thr = 0.3, thr_z = z_thr)
        self.logger.info ("Vaf shifts calculated. Shrink factor used: {:.2f}.".format (cov_mult))        
            
    def get_ai_full (self, z_thr = 2.5):
        
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
            except RuntimeError:
                dvs.append (0)
            except ValueError:
                dvs.append (0)
        dva = np.array (dvs)
        #median = np.median (dva[dva > 0])
        #dva[dva == 0] = median        
        self.dv = dva
        self.v0 = np.array(v0s)
        #p_thr is lower that for sensitive as full is more noisy, but less nosy :D 
        self.dv_dist = Distribution.Distribution (self.dv,
                                                  p_thr = 0.1, thr_z = z_thr)
    
    def get_coverage (self, z_thr = 1.5):
        ml = []
        ll = []
        
        if self.symbol == Chromosome.E_SYMBOL:
            p_thr = 0.3
        else:
            p_thr = 0.1
            
        for window in self.windows:
            #cov = window['cov'].values
            result = Testing.COV_test (window)
            ml.append (result.m)
            ll.append (result.l)

        self.m = np.array (ml)
        self.m_dist = Distribution.Distribution (self.m, thr_z = z_thr, p_thr = p_thr)
        self.l = np.array (ll)
        self.l_dist = Distribution.Distribution (self.l, thr_z = z_thr, p_thr = p_thr)
    
    def solve_windows (self, chi2_thr = 13.6):
        
        self.solutions = []
        x = np.array([self.dv, self.m, self.l]).T
        
        for m0, s0, labels in zip(*self.get_distributions()):
            self.logger.debug (f'Calculating solution: m = {m0}, s = {s0}') 
            y = ((x[:,:,np.newaxis] - m0[np.newaxis,:,:])/s0[np.newaxis,:,:])**2
            z = y.sum(axis = 1)
            dist_index = np.asarray(z == z.min(axis = 1)[:,np.newaxis]).nonzero()
            segments = [labels[i] for i in dist_index[1]]
             
            for i in np.where(z.min(axis = 1) > chi2_thr)[0]:
                segments [i] = 'O'
            indexes, merged_segments = merge_symbols (''.join(segments))
            chi2 = z[dist_index]
            
            new_indexes = []
            for i in indexes:
                new_indexes.append(self.dv, *i)

            self.logger.debug (f'Segment further divided: {len(new_indexes) > len(indexes)}')
            
            psl = []
            for si, ei in indexes:
                psl.append ((get_norm_p (self.dv[si:ei+1]),
                             get_norm_p (self.m[si:ei+1],
                             get_norm_p (self.l[si:ei+1]))))
            
            noOfilter = [s != 'O' for s in merged_segments]                                   
                   
            df = 2*(np.array([self.dv_dist.fail_normal(), 
                              self.m_dist.fail_normal(),
                              self.l_dist.fail_normal()], dtype = int)+1).sum()
            
            self.solutions.append(Solution (chi2 = chi2.sum()/(3*len(self.dv)-df),
                                            chi2_noO = chi2[noOfilter].sum()/(3*sum(noOfilter)-df),
                                            positions = [(self.windows_positions[si][0], self.windows_positions[ei][1]) for si, ei in indexes],
                                            p_norm = psl,
                                            segments = ''.join(segments),
                                            merged_segments = make_rle_string(''.join(merged_segments))))
        
        self.solutions.sort (key = lambda x: x.chi2_noO)        
        best_runs = ','.join (['(' + str(s)+ ',' + str(e)+ ')' for s,e in self.solutions[0].positions])
        self.logger.info ('Best solution: ' + best_runs + str(df) + ',' + str(sum(noOfilter)))
       
    def get_distributions (self):
                
        dvm, dvs = self.dv_dist.combinations_of_params (size = 1, key = 'single', reverse = False)
        mm, ms = self.m_dist.combinations_of_params (size = 1, key = 'single', reverse = False)
        lm, ls = self.l_dist.combinations_of_params (size = 1, key = 'single', reverse = False)
        zml = [np.array ((dvm, mm, lm))[:, np.newaxis]]
        zsl = [np.array ((dvs, ms, ls))[:, np.newaxis]]
        labels = [('B',)]
        
        ordinal = np.cumsum ([self.dv_dist.fail_normal(), self.m_dist.fail_normal(), self.l_dist.fail_normal()])
        
        if any (ordinal != 0):
            dv_directions = [False]
            m_directions = [False, True] if ordinal[1] > 1 else [False]
            l_directions = [False, True] if ordinal[2] > 1 else [False]
            
            for dv_d in dv_directions:
                for m_d  in m_directions:
                    for l_d in l_directions:
                        dvm, dvs = self.dv_dist.combinations_of_params (size = 2, key = self.dv_dist.key, reverse = dv_d)
                        mm, ms = self.m_dist.combinations_of_params (size = 2, key = self.m_dist.key, reverse = m_d)
                        lm, ls = self.l_dist.combinations_of_params (size = 2, key = self.l_dist.key, reverse = l_d)
                        
                        zml.append (np.array ((dvm, mm, lm)))
                        zsl.append (np.array ((dvs, ms, ls)))
                        
                        labels.append (('C','D'))
        return zml, zsl, labels 

    def tostring (self):
        return self.name + '-' + self.symbol
    
    def report (self, report_type = 'short'):
        report_types = ['short', 'full', 'solution']
        if report_type not in report_types:
            warn.warn ('Unknown report type. Use "short" instead.')
            report_type = 'short'
        
        if report_type == 'short':
            report = ';'.join([self.name, self.symbol, 
                               f'#solutions: {len(self.solutions)}',
                               f'#segments: {len(self.solutions[0].positions)}'])
        elif report_type == 'full':
            report = ';'.join([self.name, self.symbol,
                               f'#solutions: {len(self.solutions)}',
                               f'#segments: {len(self.solutions[0].positions)}'])
        elif report_type == 'solution':
            reports = [self.name]
            for solution in self.solutions:
                sol_str = '    ' + '; '.join([str(solution.chi2), str(solution.chi2_noO), str(solution.positions), solution.merged_segments])
                reports.append (sol_str)
            report = '\n'.join (reports)
        return report
    
    def __repr__(self) -> str:
        return self.tostring()
        
def get_norm_p (values, sinit = 0.05):
    """Tests if normally distributed around zero."""
    def cdf(x,s):
        return sts.norm.cdf (x, 0, s)
    
    try:
        popt, pcov = opt.curve_fit (cdf, np.sort(values), np.linspace (0,1, len(values)), p0 = [sinit])
        pks = sts.kstest (values, cdf, args = popt).pvalue
    except opt.OptimizeWarning:
        pks = np.nan    
    return pks

def vaf_cnai (v, dv, a, vaf,b, cov):
    s = np.sqrt((vaf - dv)*(vaf + dv)/(b*cov))
    return a*sts.norm.cdf (v, vaf - dv, s) + (1-a)*sts.norm.cdf (v, vaf + dv, s)

def vaf_HO (v, lerr):
    err = 10**lerr
    return np.exp ((v-1)*err)

def merge_symbols (in_string, outliers_threshold = 2):
    string = list(in_string)
    symbols, counts = rle_encode (string)

    argsort = np.argsort (counts)[::-1]
    symbolindex = 0

    while symbolindex < len(symbols)-1:  #(min(counts) <= outliers_threshold)&(counter < 10):
        i = argsort[symbolindex]
        symbol = symbols[i]
        s = i
        e = i
        
        sinit = i
        einit = i    
        
        if e < len(counts)-1:
            inext = e + 1        
            while (symbols[inext] == symbol)|(counts[inext] <= outliers_threshold):
                e = inext
                if e == len(counts)-1:
                    break
                inext += 1
            
        if s > 0:
            iprev = s - 1
            while (counts[iprev] <= outliers_threshold)|(symbols[iprev] == symbol):
                s = iprev
                if s == 0:
                    break
                iprev -= 1
        
        if (s == sinit)&(e == einit):
            symbolindex +=1
        else:
            for i in np.arange(counts[:s].sum(),counts[:e+1].sum()):
                string[i] = symbol.copy()
        
            symbols, counts = rle_encode (string)
            argsort = np.argsort (counts)[::-1]
            symbolindex = 0
        
    bed = []
    for i in np.arange(len(counts)):
        s = counts[:i].sum()
        e = counts[:i+1].sum()-1
        bed.append ((s,e))
    return bed, string

def rle_encode (string):
    i = 0
    values = []
    counts = []
    while i < len(string):
        cur_char = string[i]
        count = 1
        while (i < len(string) - 1 and string[i] == string[i+1]):
            count += 1
            i += 1
        values.append (cur_char)
        counts.append (count)
        i += 1
    return np.array(values), np.array(counts)

def make_rle_string(string, sep = ';'):
    values, counts = rle_encode (string)
    rle_string = []
    for v,c in zip(values, counts):
        rle_string.append (str(c)+v)        
    return sep.join(rle_string)
   
def divide_segment (dv, si, ei):
    parameters = Distribution.fit_double_G (dv[si:ei], alpha = 0.05, r = 0.1)
    threshold = get_two_G_threshold (parameters)
    random_length = get_random_lenghts (parameters, ei-si, threshold)
         
    values, counts = rle_encode (['A' if v < threshold else 'B' for v in dv[si:ei]])
    
    new_segments = []
    
    indexes = np.where(counts[values == 'B'] > threshold[1])[0]
    i = 0
    start = si
    while i < len (indexes):
        end = counts[:i].sum()
        new_segments.append((start, end))
        start = end+1
    
    if start < ei:
        new_segments.append ((start, ei))
    
    return new_segments
    
    
    #merge those regions

    #return []

def get_two_G_threshold (params):

    a = params['a']
    m = params['m']
    s = params['s']
    
    x = np.arange (m[0] - 2*s[0], m[1] + 2*s[1], 0.001)
    g0sf = a[0]*sts.norm.sf (x, m[0], s[0])
    g1cdf = a[1]*sts.norm.cdf (x, m[1], s[1])

    thr = x[np.where(g0sf < g1cdf)[1].min[0]]
    return thr

def get_random_lenghts (params, size, thr, tries = 1000):
    maxs = []
    for _ in range(tries):
        rng = default_rng()
        vals = rng.uniform(size = 71)
        rdv = ppf (vals, params)
        values, counts = rle_encode (['A' if v < thr else 'B' for v in rdv])
        maxs.append ((counts[values == 'A'].max(), counts[values == 'B'].max()))
    return np.array(maxs).max(axis = 0)

def ppf (x, p):
    a = p['a']
    m = p['m']
    s = p['s']
    return a[0]*sts.norm.ppf (x, m[0], s[0]) + a[1]*sts.norm.ppf (x, m[1], s[1])
