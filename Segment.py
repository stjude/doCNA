import scipy.stats as sts
import numpy as np
from collections import namedtuple

import scipy.optimize as opt

from doCNA import Testing
from doCNA import Run


E_SYMBOL = 'E'
MAX_CLON_THRESHOLD_FOR_SENSITIVE = 0.4
MIN_CLON_THRESHOLD_FOR_FULL = 0.2


class Segment:
    """Class to calculate clonality and find the model."""
    def __init__ (self, data, config, logger, genome_medians, segmentation_score, segmentation_symbol) -> None:
        self.data = data
        self.config = config
        self.genome_medians = genome_medians
        self.segmentation_score = segmentation_score
        self.segmentation_symbol = segmentation_symbol
        self.name = data['chrom'].values[0]+ ':'+str(data['position'].min())+'-'+str(data['position'].max())
        #-{self.name}
        self.logger = logger.getChild (f'{self.__class__.__name__}-{self.name}')
        self.logger.debug (f'Segment created: {self.segmentation_symbol} .')
        self.symbols = self.data.loc[self.data['vaf'] < 1, 'symbol'].value_counts()
        self.symbol = self.symbols.index[0]
        self.logger.debug (f'Segment symbol: {self.symbol}')
        self.estimate_parameters ()
        self.logger.debug ('Parameters estimated.')
        self.select_model ()
        self.logger.debug ('Model selected.')
        #self.test_self ()
        
    def tostring(self) -> str:
        return '\n'.join ([self.name, str(self.parameters)])
            
    def __repr__(self) -> str:
        return '\n'.join ([self.name, str(self.parameters)])
    
    def estimate_parameters (self):
        #sometimes, even if classified, it's wrong model
        
        if self.symbol == E_SYMBOL:
            self.parameters = get_sensitive (self.data.loc[self.data['symbol'] == E_SYMBOL,],
                                             self.genome_medians['VAF']['fb'],
                                             self.genome_medians['COV']['m'])
            if self.parameters['ai'] > MAX_CLON_THRESHOLD_FOR_SENSITIVE:
                self.logger.info (f"Estimated ai: {self.parameters['ai']} above threshold for sensitive model: {MAX_CLON_THRESHOLD_FOR_SENSITIVE}")
                self.parameters = get_full (self.data)
                                            #self.genome_medians['VAF']['fb'],
                                            #self.genome_medians['HE']['b'])
            self.logger.info (f"Estimated, by sensitive method, ai: {self.parameters['ai']}")
        else:
            self.parameters = get_full (self.data)
                                        #self.genome_medians['VAF']['fb'],
                                        #self.genome_medians['HE']['b'])
            #if self.parameters['ai'] < MIN_CLON_THRESHOLD_FOR_FULL:
            #    self.logger.info (f"Estimated ai {self.parameters['ai']} below threshold for full model: {MIN_CLON_THRESHOLD_FOR_FULL}")
            #    self.parameters = get_sensitive (self.data.loc[self.data['vaf'] < 1 - 1/self.genome_medians['COV']['m']],
            #                                     self.genome_medians['VAF']['fb'],
            #                                     self.genome_medians['COV']['m'])
            self.logger.info (f"Estimated, by full method, ai: {self.parameters['ai']}.")
    
    def report (self, report_type = 'bed'):
        namestr = self.name.replace (':', '\t').replace ('-', '\t')
        
        if report_type == 'bed':
            gmm = self.genome_medians['clonality']['m']
            gms = self.genome_medians['clonality']['s']
            n = self.parameters['n']/Run.SNPS_IN_WINDOW
            score = np.abs(self.parameters['k'] - gmm)*np.sqrt(n/2)/gms
            report = '\t'.join([str(p) for p in [self.parameters['m'], self.parameters['model'], self.parameters['k'], score]])
        else:
            report = ''
        return namestr + '\t' + report
        
    def select_model (self):        
        if self.parameters['success']:
            m = self.parameters['m']
            v = self.parameters['ai']
            m0 = self.genome_medians['COV']['m']
        
            self.distances = np.array ([calculate_distance (preset, m,v,m0) for preset in model_presets.values()])
        
            #ordering of presets will prioritize models
            # cn2 should be before cn4, for example
            picked = np.where(self.distances == self.distances.min())[0][0]
            
            self.parameters['d'] = self.distances.min()
            self.parameters['model'] = list(model_presets.keys())[picked]
            k = model_presets[self.parameters['model']].k(m,v,m0) 
            self.parameters['k'] = k if k < 1 else np.nan
        else:
            self.parameters['model'] = 'NA'
            self.parameters['k'] = np.nan
            self.logger            

def calculate_distance (preset, m, ai, m0):
    return np.abs (preset.C (m/m0,ai,1) - preset.D (m/m0,ai,1))/np.sqrt (preset.A(m/m0,ai,1)**2 + preset.B(m/m0,ai,1)**2)
    
Preset = namedtuple ('Preset', ['A', 'B', 'C', 'D', 'k'])
model_presets = {'cn1' : Preset(A = lambda m,dv,m0: -m0/2,
                                B = lambda m,dv,m0: -1,
                                C = lambda m,dv,m0: m0,
                                D = lambda m,dv,m0: m0*(2*dv/(0.5+dv))/2+m,
                                k = lambda m,dv,m0: 2*dv/(0.5+dv)),
                 
                 'cnL' : Preset(A = lambda m,dv,m0: 0,
                                B = lambda m,dv,m0: -1,
                                C = lambda m,dv,m0: m0,
                                D = lambda m,dv,m0: m,
                                k = lambda m,dv,m0: 2*dv),
                 
                 'cn3' : Preset(A = lambda m,dv,m0: m0/2,
                                B = lambda m,dv,m0: -1,
                                C = lambda m,dv,m0: m0,
                                D = lambda m,dv,m0: -m0*(2*dv/(0.5-dv))/2+m,
                                k = lambda m,dv,m0: 2*dv/(0.5-dv)),

                 'cnB' : Preset(A = lambda m,dv,m0: 0,
                                B = lambda m,dv,m0: 1,
                                C = lambda m,dv,m0: 2*dv,
                                D = lambda m,dv,m0: 0,
                                k = lambda m,dv,m0: np.abs(m/m0 - 1))}
    
def get_sensitive (data, fb, mG, z_thr = 1.5):
    
    vafs = data['vaf'].values
    covs = data['cov'].values
    
    #this only works for E
    def ai (v, dv, a):
        v0 = 0.5
        return a*sts.norm.cdf (v, v0 - dv, smin) + (1-a)*sts.norm.cdf (v, v0 + dv, smin)
    
    m, dm, l, dl = Testing.COV_test (data)
    smin = np.sqrt (0.25/m)*fb
           
    v,c = np.unique (vafs, return_counts = True)
    try:
        popt, pcov = opt.curve_fit (ai, v, np.cumsum(c)/np.sum(c), p0 = [0.05, 0.5],
                                    bounds = ((0.0, 0.0),
                                              (0.5, 1.0)))
        dv, a = popt
        parameters = {'m': m, 'l': l, 'ai' : dv, 'a': a, 'success' : True, 'n' : len (data)}
        
    except RuntimeError:
        parameters = {'m': m, 'l': l, 'ai' : np.nan, 'success' : False, 'n' : 0}
    
    return parameters 

def get_full (data, b = 1.01):
    
    vafs = data['vaf'].values
    covs = data['cov'].values
    
    m, dm, l, dl = Testing.COV_test (data)
    #cov = mG
    
    def vaf_cdf (v, dv, a, lerr, f, vaf, b):
        return vaf_cdf_c (v, dv, a, lerr, f, vaf, b, m)
    
    v0 = 0.5
    #why on earth there is 0.09?!    
    #s0 = np.sqrt (0.09/m)
    v, c = np.unique(vafs[~np.isnan(vafs)], return_counts = True)

    try:
        cnor = np.cumsum(c)/np.sum(c)
        ones0 = c[v >= (m-1)/m].sum()
        f0 = c[v < v0].sum()/(c.sum() - ones0) 
        dv0 = v0 - np.median (v[v < v0])

        p0 = [dv0, ones0/c.sum(), 2, 0.5, 0.5, b]
        #print ('Initial parameters: ', p0)
        popt, pcov = opt.curve_fit (vaf_cdf, v, cnor, p0 = p0, 
                                    bounds = ((0,   0,   1, 0, 0.45, 1),
                                              (0.5, 0.95, 5, 1, 0.55, 10)))
        dv, a, lerr, f, vaf, b = popt
        parameters = {'m': m, 'l': l, 'ai' : dv, 'v0': v0, 'a': a, 'b' : b, 'success' : True, 'n' : len (data), 'status' : 'valid'} 
    except RuntimeError:
        parameters = {'m': m, 'l': l, 'ai' : np.nan, 'success' : False, 'n' : 0, 'status' : 'Fit failed'}
        #print ('Runtime: Initial parameters: ', p0)
    except ValueError:
        parameters = {'m': m, 'l': l, 'ai' : np.nan, 'success' : False, 'n' : 0, 'status' : 'Parameters failed'}
        #print ('Value: Initial parameters: ', p0)
        
    return parameters

def vaf_cnai (v, dv, a, vaf,b, cov):
    s = np.sqrt((vaf - dv)*(vaf + dv)/(b*cov))
    return a*sts.norm.cdf (v, vaf - dv, s) + (1-a)*sts.norm.cdf (v, vaf + dv, s)

def vaf_HO (v, lerr):
    err = 10**lerr
    return np.exp ((v-1)*err)

def vaf_cdf_c (v, dv, a, lerr, f, vaf, b, cov):
    #cn2 = vaf_cn2 (v, vaf, cov)
    cnai = vaf_cnai (v, dv, f, vaf, b, cov)
    cnHO = vaf_HO (v, lerr)
    
    return a*cnHO + (1 - a)*cnai 
