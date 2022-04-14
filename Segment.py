import scipy.stats as sts
import numpy as np
from collections import namedtuple

import scipy.optimize as opt

from doCNA import Testing

E_SYMBOL = 'E'

class Segment:
    def __init__ (self, data, config, logger, genome_medians):
        self.data = data
        self.config = config
        self.genome_medians = genome_medians
        self.name = data['chrom'][0]+ ':'+str(data['position'].min())+'-'+str(data['position'].max())
        self.logger = logger.getChild (f'{self.__class__.__name__}-{self.name}')
        self.symbols = self.data.loc[self.data['vaf'] < 1, 'symbol'].value_counts()
        self.symbol = self.symbols.index[0]
        self.estimate_parameters ()
        self.logger.info ('Segment created.')
        
    def __repr__(self) -> str:
        return '\n'.join ([self.name, str(self.parameters)])
    
    def estimate_parameters (self):
        if self.symbol == E_SYMBOL:
            self.parameters = get_sensitive (self.data.loc[self.data['symbol'] == E_SYMBOL, 'vaf'].values,
                                             self.data['cov'].values,
                                             self.genome_medians['VAF']['fb'],
                                             self.genome_medians['COV']['m'])
        else:
            self.parameters = get_full (self.data['vaf'].values,
                                        self.data['cov'].values,
                                        self.genome_medians['VAF']['fb'],
                                        self.genome_medians['HE']['b'])
        
    def test_self (self):
        pass
    
    def report (self, type = 'bed'):
        pass
    
    def select_model (self):
        #add 'model' and 'clonality' and 'distance' keys to self.parameters
        
        m = self.parameters.m
        v = self.parameters.ai
        m0 = self.genome_medians['COV']['m']
        
        self.distances = np.array ([calculate_distance (preset, m,v,m0) for preset in model_presets.values()])
        
        picked = np.where(self.distance == self.distance.min())[0][0]
                
        self.parameters['model'] = model_presets.keys()[picked]
        self.parameters['k'] = model_presets[self.parameters['model']].k(m,v,m0)
                

def calculate_distance (preset, m, ai, m0):
    return np.abs (preset.C (m,ai,m0) - preset.D (m,ai,m0))/np.sqrt (preset.A(m,ai,m0)**2 + preset.B(m,ai,m0)**2)
    

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
                                k = lambda m,dv,m0: 2*dv/(0.5-dv))}
    
def get_sensitive (vafs, covs, fb, mG, z_thr = 1.5):
    #this only works for E
    def ai (v, dv, v0, a):
        return a*sts.norm.cdf (v, v0 - dv, smin) + (1-a)*sts.norm.cdf (v, v0 + dv, smin)
    
    m, l = Testing.COV_test (covs)
    smin = np.sqrt (0.25/m.value)*fb
           
    v,c = np.unique (vafs, return_counts = True)
    popt, pcov = opt.curve_fit (ai, v, np.cumsum(c)/np.sum(c), p0 = [0.05, 0.5, 0.5],
                                bounds = ((0.0, 0.4, 0.0),
                                          (0.5, 0.6, 1.0)))
    dv, v0, a = popt
    return {'m': m, 'l': l, 'ai' : dv, 'v0': v0, 'a': a } 

def get_full (vaf, cov, mG, b):
    
    m, l = Testing.COV_test (cov)
    cov = mG
    
    def vaf_cdf (v, dv, a, lerr, f, vaf, b):
        return vaf_cdf_c (v, dv, a, lerr, f, vaf, b, cov)
    
    v0 = 0.5
    s0 = np.sqrt (0.09/m.value)
    v, c = np.unique(vaf[~np.isnan(vaf)], return_counts = True)

    cnor = np.cumsum(c)/np.sum(c)
    ones0 = c[v >= (cov-1)/cov].sum()
    f0 = c[v < v0].sum()/(c.sum() - ones0) 

    dv0 = v0 - np.median (v[v < v0])

    popt, pcov = opt.curve_fit (vaf_cdf, v, cnor, p0 = [dv0, ones0/c.sum(), 2, f0, 0.5, b], 
                                bounds = ((0,   0,   1, 0, 0.45, 1),
                                          (0.5, 0.95, 5, 1, 0.55, 10)))
            
    dv, a, lerr, f, vaf, b = popt
    
    return {'m': m, 'l': l, 'ai' : dv, 'v0': v0, 'a': a, 'b' : b} 
    

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