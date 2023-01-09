import scipy.stats as sts
import numpy as np

import scipy.optimize as opt

from doCNA import Testing
from doCNA import Consts
from doCNA import Models
from doCNA.Report import Report

model_presets = {}
model_presets.update (Models.model_presets_2)
model_presets.update (Models.model_presets_4)

class Segment:
    """Class to calculate clonality and find the model."""
    def __init__ (self, data, config, logger, genome_medians, segmentation_score,
                  segmentation_symbol, centromere_fraction, cytobands) -> None:
        self.data = data
        self.config = config
        self.genome_medians = genome_medians
        self.segmentation_score = segmentation_score
        self.segmentation_symbol = segmentation_symbol
        self.centromere_fraction = centromere_fraction
        self.cytobands = cytobands
        
        self.chrom = data['chrom'].values[0]
        self.start = data['position'].min()
        self.end = data['position'].max()
        self.name = data['chrom'].values[0]+ ':'+str(data['position'].min())+'-'+str(data['position'].max())
        self.logger = logger.getChild (f'{self.__class__.__name__}-{self.name}')
        self.logger.debug ('Segment created.')
        self.symbols = self.data.loc[self.data['vaf'] < 1, 'symbol'].value_counts()
        self.symbol = self.symbols.index[0]
        self.logger.debug (f'Segment symbol: {self.symbol}')
        self.estimate_parameters ()
        self.select_model ()
        self.logger.debug ('Segment analyzed.')

 
        
    def tostring(self) -> str:
        return '\n'.join ([self.name, str(self.parameters)])
            
    def __repr__(self) -> str:
        return '\n'.join ([self.name, str(self.parameters)])
    
    def estimate_parameters (self):
        #sometimes, even if classified, it's wrong model
        method = 'unspecified'
        if self.symbol == Consts.E_SYMBOL:
            self.parameters = get_sensitive (self.data.loc[self.data['symbol'] == Consts.E_SYMBOL,],
                                             self.genome_medians['fb'],
                                             self.genome_medians['m'])
            method = 'sensitive'
            if self.parameters['ai'] > Consts.MAX_AI_THRESHOLD_FOR_SENSITIVE:
                self.parameters = get_full (self.data.loc[self.data['symbol'] != 'A',])
                method = 'full'
            
        else:
            self.parameters = get_full (self.data.loc[self.data['symbol'] != 'A',])
            method = 'full'
            
        if self.parameters['success']:
            self.logger.info (f"AI estimated by {method} method, ai = {self.parameters['ai']}")
        else:
            self.logger.info (f"AI not estimated.")
            self.logger.debug (f"Parameters: {self.parameters}")
                
    def report (self, report_type = 'bed'):
        return Report(report_type).segment_report(self)
                                                  
    def select_model (self):        
        if self.parameters['success']:
            m = self.parameters['m']
            v = self.parameters['ai']
            m0 = self.genome_medians['m']
        
            self.distances = np.array ([Models.calculate_distance (preset, m,v,m0) for preset in Models.model_presets.values()])
            picked = np.where(self.distances == self.distances.min())[0][0]
                        
            self.parameters['d'] = self.distances.min()
            self.parameters['model'] = list(model_presets.keys())[picked]
            k = model_presets[self.parameters['model']].k(m,v,m0) 
            self.parameters['k'] = k if k < Consts.K_MAX else np.nan
            self.logger.info (f"Segment identified as {self.parameters['model']}, d = {self.parameters['d']}")
        else:
            self.parameters['d'] = np.nan
            self.parameters['model'] = 'NA'
            self.parameters['k'] = np.nan
            self.logger.info ('No model for this segment.')
                        


def get_sensitive (data, fb, mG, z_thr = 1.5):
    
    vafs = data['vaf'].values
    
    def ai (v, dv, a):
        v0 = 0.5
        return a*sts.norm.cdf (v, v0 - dv, smin) + (1-a)*sts.norm.cdf (v, v0 + dv, smin)
    
    m, dm, l, dl = Testing.COV_test (data)
    smin = np.sqrt (0.25/m)*fb
           
    v,c = np.unique (vafs, return_counts = True)
    try:
        popt, pcov = opt.curve_fit (ai, v, np.cumsum(c)/np.sum(c), p0 = [0.02, 0.5],
                                    bounds = ((0.0, 0.1),
                                              (0.3, 0.9)), check_finite = False)
        dv, a = popt
        parameters = {'m': m, 'l': l, 'ai' : dv, 'a': a, 'success' : True, 'n' : len (data)/Consts.SNPS_IN_WINDOW,
                      'status' : 'valid', 'fraction_1' : np.nan}
        
    except (RuntimeError, ValueError):
        parameters = {'m': m, 'l': l, 'ai' : np.nan, 'success' : False, 'n' : np.nan,
                      'status' : 'valid', 'fraction_1' : np.nan}
        
    
    return parameters 

def get_full (data, b = 1.01):
    
    vafs = data['vaf'].values    
    m, dm, l, dl = Testing.COV_test (data)

    def vaf_cdf (v, dv, a, lerr, f, vaf, b):
        return vaf_cdf_c (v, dv, a, lerr, f, vaf, b, m)
    
    v0 = 0.5
    v, c = np.unique(vafs[~np.isnan(vafs)], return_counts = True)
    try:
        cnor = np.cumsum(c)/np.sum(c)
        ones0 = c[v >= (m-1)/m].sum()        
        f0 = c[v < v0].sum()/(c.sum() - ones0) 
        dv0 = v0 - np.median (v[v < v0])
        p0 = [dv0, ones0/c.sum(), 2, 0.5, 0.5, b]        
        popt, pcov = opt.curve_fit (vaf_cdf, v, cnor, p0 = p0, 
                                    bounds = ((0,   0,   1, 0, 0.45, 1),
                                              (0.5, 0.95, 5, 1, 0.55, 10)))
        dv, a, lerr, f, vaf, b = popt      
        parameters = {'m': m, 'l': l, 'ai' : dv, 'v0': v0, 'a': a, 'b' : b, 'success' : True, 
                      'fraction_1' : ones0/c.sum(), 'n' : len (data)/Consts.SNPS_IN_WINDOW,
                      'status' : 'valid'}
    except RuntimeError:
        parameters = {'m': m, 'l': l, 'ai' : np.nan, 'success' : False, 'n' : np.nan,
                       'fraction_1' : ones0/c.sum(), 'status' : 'Fit failed'}
    except ValueError:
        parameters = {'m': m, 'l': l, 'ai' : np.nan, 'success' : False, 'n' : np.nan,  
                      'fraction_1' : ones0/c.sum(), 'status' : 'Parameters failed'}    
    if ones0/c.sum() > 0.9:
        parameters = {'m': m, 'l': l, 'ai' : 0.5, 'success' : True, 'n' : len (data)/Consts.SNPS_IN_WINDOW,
                      'fraction_1' : ones0/c.sum(), 'status' : 'Parameters guessed'}    
        
    return parameters

def vaf_cnai (v, dv, a, vaf,b, cov):
    s = np.sqrt((vaf - dv)*(vaf + dv)/(b*cov))
    return a*sts.norm.cdf (v, vaf - dv, s) + (1-a)*sts.norm.cdf (v, vaf + dv, s)

def vaf_HO (v, lerr):
    err = 10**lerr
    return np.exp ((v-1)*err)

def vaf_cdf_c (v, dv, a, lerr, f, vaf, b, cov):
    cnai = vaf_cnai (v, dv, f, vaf, b, cov)
    cnHO = vaf_HO (v, lerr)
    return a*cnHO + (1 - a)*cnai
