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
        self.name = data['chrom'].values[0]+ ':'+str(data['position'].min())+'-'+str(data['position'].max())
        #-{self.name}
        self.logger = logger.getChild (f'{self.__class__.__name__}')
        self.symbols = self.data.loc[self.data['vaf'] < 1, 'symbol'].value_counts()
        self.symbol = self.symbols.index[0]
        self.estimate_parameters ()
        self.select_model ()
        self.test_self ()
        self.logger.debug ('Segment created.')
        
    def __repr__(self) -> str:
        return '\n'.join ([self.name, str(self.parameters)])
    
    def estimate_parameters (self):
        if self.symbol == E_SYMBOL:
            self.parameters = get_sensitive (self.data.loc[self.data['symbol'] == E_SYMBOL,],
                                             self.genome_medians['VAF']['fb'],
                                             self.genome_medians['COV']['m'])
                
        else:
            self.parameters = get_full (self.data,
                                        self.genome_medians['VAF']['fb'],
                                        self.genome_medians['HE']['b'])
        
    def test_self (self):
        if self.parameters['success']:
            m0 = self.parameters['m']
            v = self.data['vaf'].values
            m = self.data['cov'].values
        
            self.self_test_p = [get_norm_p (v - np.median(v)),
                                  get_norm_p (m - m0)]
        else:
            self.self_test_p = [np.nan, np.nan]
    
    def report (self, type = 'bed'):
        namestr = self.name.replace (':', '\t').replace ('-', '\t')
        #print (namestr + [self.parameters['m'], self.parameters['ai'],self.parameters['model'], self.parameters['k']])
        string = '\t'.join([namestr] + [str(self.parameters['m']), str(self.parameters['ai']),
                                      self.parameters['model'], str(self.parameters['k'])])
        
        return string
    
    def select_model (self):
        #add 'model' and 'clonality' and 'distance' keys to self.parameters
        if self.parameters['success']:
            m = self.parameters['m']
            v = self.parameters['ai']
            m0 = self.genome_medians['COV']['m']
        
            self.distances = np.array ([calculate_distance (preset, m,v,m0) for preset in model_presets.values()])
        
            picked = np.where(self.distances == self.distances.min())[0][0]
                
            self.parameters['model'] = list(model_presets.keys())[picked]
            self.parameters['k'] = model_presets[self.parameters['model']].k(m,v,m0)
        else:
            self.parameters['model'] = 'NA'
            self.parameters['k'] = np.nan
            
def get_norm_p (values, sinit = 0.05):
    """Tests if normally distributed around zero."""
    def cdf(x,s):
        return sts.norm.cdf (x, 0, s)
    
    popt, pcov = opt.curve_fit (cdf, np.sort(values), np.linspace (0,1, len(values)), p0 = [sinit])
    pks = sts.kstest (values, cdf, args = popt).pvalue
    return -np.log10 (pks)

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
                                k = lambda m,dv,m0: 2*dv/(0.5-dv))}
    
def get_sensitive (data, fb, mG, z_thr = 1.5):
    
    vafs = data['vaf'].values
    covs = data['cov'].values
    
    #this only works for E
    def ai (v, dv, v0, a):
        return a*sts.norm.cdf (v, v0 - dv, smin) + (1-a)*sts.norm.cdf (v, v0 + dv, smin)
    
    m, dm, l, dl = Testing.COV_test (data)
    smin = np.sqrt (0.25/m)*fb
           
    v,c = np.unique (vafs, return_counts = True)
    try:
        popt, pcov = opt.curve_fit (ai, v, np.cumsum(c)/np.sum(c), p0 = [0.05, 0.5, 0.5],
                                    bounds = ((0.0, 0.4, 0.0),
                                              (0.5, 0.6, 1.0)))
        dv, v0, a = popt
        parameters = {'m': m, 'l': l, 'ai' : dv, 'v0': v0, 'a': a, 'success' : True}
        
    except:
        parameters = {'m': m, 'l': l, 'ai' : np.nan, 'success' : False}
    
    return parameters 

def get_full (data, mG, b):
    
    vafs = data['vaf'].values
    covs = data['cov'].values
    
    m, dm, l, dl = Testing.COV_test (data)
    cov = mG
    
    def vaf_cdf (v, dv, a, lerr, f, vaf, b):
        return vaf_cdf_c (v, dv, a, lerr, f, vaf, b, cov)
    
    v0 = 0.5
    s0 = np.sqrt (0.09/m)
    v, c = np.unique(vafs[~np.isnan(vafs)], return_counts = True)

    cnor = np.cumsum(c)/np.sum(c)
    ones0 = c[v >= (cov-1)/cov].sum()
    f0 = c[v < v0].sum()/(c.sum() - ones0) 

    dv0 = v0 - np.median (v[v < v0])

    try:
        popt, pcov = opt.curve_fit (vaf_cdf, v, cnor, p0 = [dv0, ones0/c.sum(), 2, f0, 0.5, b], 
                                    bounds = ((0,   0,   1, 0, 0.45, 1),
                                              (0.5, 0.95, 5, 1, 0.55, 10)))
        dv, a, lerr, f, vaf, b = popt
        parameters = {'m': m, 'l': l, 'ai' : dv, 'v0': v0, 'a': a, 'b' : b, 'success' : True} 
    except:
        parameters = {'m': m, 'l': l, 'ai' : np.nan, 'success' : False}
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