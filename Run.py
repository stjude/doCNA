#the idea here is that after initial segmentation, this stratify run into segments
import scipy.stats as sts
import numpy as np
import pandas as pd
import scipy.optimize as opt

from collections import namedtuple

from doCNA import Segment
from doCNA import Distribution
from doCNA import Testing
WINDOWS_THRESHOLD = 9
SNPS_IN_WINDOW = 1000


#Solution keep result of Run segmenting
#Segments are based on /best/ Solution
Solution = namedtuple ('Solution', ['chi2', 'chi2_noO', 'positions', 'ps', 'string', 'merged_string'])

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
         self.logger.debug ("Object run created")
                  
    def analyse (self):
        self.get_windows (n = SNPS_IN_WINDOW)
        if len(self.windows) > WINDOWS_THRESHOLD:
            self.get_ai ()
            self.get_coverage ()
            
            #stupid name, but no better idea
            self.find_z0 ()
            #solving windows adds self.solution
            self.solve_windows ()
            
        else:
            self.logger.info(f'Run {self.name} is to short to segment.')
            self.solution = Solution(chi2 = 1, chi2_noO = 1, 
                                     positions = [(self.start, self.end)],
                                     ps = [np.nan, np.nan, np.nan],
                                     string = '', merged_string = '')

    def get_windows (self, n = 1000):
        tmp = self.data.loc[[s == self.symbol for s in self.data['symbol']], ]
        
        N = int(np.floor(len(tmp)/n))
        indexes = np.linspace (0,len (tmp), 2*N+2, dtype = int)

        self.windows = []
        self.windows_positions = []

        for i in np.arange (2*N):
            tmpi = tmp.iloc[indexes[i]:indexes[i+2], ]
            self.windows.append(tmpi)
            self.widows_positions.append ((tmpi['position'].min(), tmpi['position'].max()))
        
    def get_ai (self):
        if self.symbol == 'E':
            self.get_ai_sensitive()     
        else:
            self.get_ai_full()
        
    def get_ai_sensitive (self, zero_thr = 0.01, cov_mult = 1.03, p_thr = 0.5, z_thr = 1.5):
        tmpf = 1
        s0 = np.sqrt (0.25/self.genome_medians['COV']['m'])
                
        while tmpf > zero_thr:
            dvl = []
            v0l = []
            s = s0/np.sqrt(cov_mult)
            def make_two_gauss (v, dv, v0):
                a = 0.5
                return a*sts.norm.cdf (v - dv, s) + (1-a)*sts.norm.cdf (v + dv, s)
            
            for window in self.windows:
                vaf = window['vaf'].values
                popt, pcov = opt.curve_fit (make_two_gauss, np.sort(vaf), np.linspace (0,1, len(vaf)),
                                            p0 = [0.05, 0.5],
                                            bounds = ((0, 0.4),(0.5, 0.6)))
                dvl.append (popt[0])
                v0l.append (popt[1])
        
            tmpf = sum ([d < zero_thr for d in dvl])/len (dvl)
            cov_mult += 0.05
        
        self.dv = np.array (dvl)
        self.v0 = np.array (v0l)

        self.dv_dist = Distribution.Distribution (self.dv, p_thr = 0.5, thr_z = z_thr)
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
        for vaf in self.vafs:
            v, c = np.unique(vaf, return_counts = True)

            cnor = np.cumsum(c)/np.sum(c)
            ones0 = c[v >= (cov-1)/cov].sum()
            f0 = c[v < v0].sum()/(c.sum() - ones0) 

            dv0 = v0 - np.median (v[v < v0])

            try:
                popt, pcov = opt.curve_fit (vaf_cdf, v, cnor, p0 = [dv0, ones0/c.sum(), 2, f0, 0.5, b], 
                                            bounds = ((0,   0,   1, 0, 0.45, 1),
                                                      (0.5, 0.95, 5, 1, 0.55, 10)))
                dvs.append (popt[0])
                v0s.append (popt[-1])
            except:
                #the lenghts must be same
                #adding 0 either makes outlier or not outlir but it is safe value
                dvs.append (0)
                
        self.dv = np.array(dvs)
        self.v0 = np.array(v0s)
        #p_thr is lower that for sensitive as full is more noisy, but less nosy :D 
        self.dv_dist = Distribution.Distribution (self.dv[~np.isnan(self.dv)],
                                                  p_thr = 0.3, thr_z = z_thr)
    
    def get_coverage (self, z_thr = 1.5):
        ml = []
        ll = []
            
        for fragment in self.fragments:
            result = Testing.COV_test (fragment)
            ml.append (result.m)
            ll.append (result.l)

        self.m = np.array (ml)
        self.m_dist = Distribution.Distribution (self.m, thr_z = z_thr, p_thr = 0.3)
        self.l = np.array (ll)
        self.l_dist = Distribution.Distribution (self.l, thr_z = z_thr, p_thr = 0.3)
    
    def solve_windows (self, z_thr = 2.5):
        
        self.solutions = []
        x = np.array([self.dv, self.m, self.l])   
        
        for m0, s0, labels in self.get_distributions():
            z = (((x[:,:,np.newaxis] - m0[np.newaxis,:,:])/s0[np.newaxis,:,:])**2).sum(axis = 1)
            #z.dim = SNPs * solutions 
        
            #axis = 1 as this is a dimention of solutions
            labels = self.points['labels'][(z == z.min(axis = 1)).nonzero()[1]]
            labels [z.min(axis = 1) > z_thr] = 'O'
        
            indexes, merged_string = merge_symbols (''.join(labels.tolist()))
        
            chi2 = 
                
            self.solution.append(Solution (chi2 = chi2.sum(),
                                 chi2_noO = chi2[labels != 'O'].sum(),
                                 positions = [(self.windows_positions[si][0], self.windows_positions[ei][1]) for si, ei in indexes],
                                 ps = ,
                                 string = ''.join(labels.tolist()),
                                 merged_string = merged_string))
        
        #order by chi2_noO
        
        #generate list of segments and return
    
    def get_distributions (self):
        
        
        dvm, dvs = self.dv_dist.combinations_of_params (dim = 1, reverse = False)
        mm, ms = self.dv_dist.combinations_of_params (dim = 1, reverse = False)
        lm, ls = self.dv_dist.combinations_of_params (dim = 1, reverse = False)
        zml = [np.array ([dvm, mm, lm])]
        zsl = [np.array ([dvs, ms, ls])]
        labels = [('B',)]
        
        ordinal = np.cumsum ([self.dv_dist.fail_normal(), self.m_dist.fail_normal(), self.l_dist.fail_normal()])
        
        if any (ordinal == 0):
            dv_directions = [False]
            m_directions = [False, True] if ordinal[1] > 1 else [False]
            l_directions = [False, True] if ordinal[2] > 1 else [False]
            
            for dv_d in dv_directions:
                for m_d  in m_directions:
                    for l_d in l_directions:
                        dvm, dvs = self.dv_dist.combinations_of_params (dim = 2, reverse = dv_d)
                        mm, ms = self.dv_dist.combinations_of_params (dim = 2, reverse = m_d)
                        lm, ls = self.dv_dist.combinations_of_params (dim = 2, reverse = l_d)
                        zml.append (np.array ([dvm, mm, lm]))
                        zsl.append (np.array ([dvs, ms, ls]))
                        labels.append (('C','D'))
        return zml, zsl, labels 

    
    
    def __repr__(self) -> str:
        pass
        
    def generate_segments (self, list_of_other_runs):
        
        pass
    
    
    
    
def get_norm_p (values, sinit = 0.05):
    """Tests if normally distributed around zero."""
    def cdf(x,s):
        return sts.norm.cdf (x, 0, s)
    
    popt, pcov = opt.curve_fit (cdf, np.sort(values), np.linspace (0,1, len(values)), p0 = [sinit])
    pks = sts.kstest (values, cdf, args = popt).pvalue
    return -np.log10 (pks)

def vaf_cnai (v, dv, a, vaf,b, cov):
    s = np.sqrt((vaf - dv)*(vaf + dv)/(b*cov))
    return a*sts.norm.cdf (v, vaf - dv, s) + (1-a)*sts.norm.cdf (v, vaf + dv, s)

def vaf_HO (v, lerr):
    err = 10**lerr
    return np.exp ((v-1)*err)

def merge_symbols (in_string, outliers_threshold = 2):
    string = list(in_string)
    symbols, counts = rle_encode (string)

    argsort = np.argsort (counts)[-1:0:-1]
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
            argsort = np.argsort (counts)[-1:0:-1]
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