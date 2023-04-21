import numpy as np
import scipy.stats as sts
import scipy.optimize as opt
import sklearn.linear_model as slm

from doCNA import Consts
from doCNA import Models


class Scoring:
    def __init__(self, initial_data = None, logger = False) -> None:
        
        if initial_data is None:
            self.ai_param = {'m' : np.inf, 's' : 1}
            self.cn_param = {'m' : np.inf, 's' : 1}
            self.dipl_dist = {'m' : np.inf, 's' : 1, 'thr' : 0, 'alpha': np.nan}
        else:
            self.ai_param = fit_QQgauss(initial_data[: ,0])
            self.cn_param = fit_QQgauss(initial_data[: ,1], fit_intercept = False)
            dds =  initial_data - np.array([self.ai_param['m'], self.cn_param['m']])[np.newaxis, :]
            ds = dds/np.array([self.ai_param['s'],self.cn_param['s']])[np.newaxis, :]
            self.dipl_dist = fit_QQgauss (np.sqrt((ds**2).sum(axis = 1)))
            self.dipl_dist['alpha'] = 1/(5*len(initial_data))
            self.dipl_dist['thr'] = sts.norm.ppf (1-self.dipl_dist['alpha'], 
                                                  self.dipl_dist['m'],
                                                  self.dipl_dist['s'])
        
            
        
        if logger:
            self.logger = logger.getChild (f'{self.__class__.__name__}')
            self.logger.info (f"Distribution of distance to diploid (0,2): m = {self.dipl_dist['m']}, s = {self.dipl_dist['s']}")
            self.logger.info (f"Distribution of diploid allelic imbalance: m = {self.ai_param['m']}, s = {self.ai_param['s']}")
            self.logger.info (f"Distribution of diploid copy number: m = {self.cn_param['m']}, s = {self.cn_param['s']}")
            
    def get_d_thr (self):
        return self.dipl_dist['thr']
    
    def score_dipl (self, ai, m, m0, models):
        cn = 2*m/m0
        s_ai = self.ai_param['s']
        s_cn = self.cn_param['s']
        d = np.sqrt ((ai/s_ai)**2 + ((cn-2)/s_cn)**2)
        p_d = sts.norm.sf (d, self.dipl_dist['m'], self.dipl_dist['s'])
        isHE = d < self.dipl_dist['thr']
        
        if isHE:
            #it is diploid
            model_param = {'model' : 'AB', 'd_model' : d,
                           'score_HE' : -np.log10(p_d), 'd_HE' : d,
                           'p_model' : 0,#p_d,
                           'k': cn/2-1}
        else:
            model_param = Models.pick_model (ai, s_ai, cn, s_cn, models)        
            model_param['p_model'] = 0, #sts.norm.sf(model_param['d_model'], self.dipl_dist['m'], self.dipl_dist['s'])
            model_param['score_HE'] = -np.log10(p_d)
            model_param['d_HE'] = d
        return model_param
    
    def analyze_segment (self, segment, models):
        """Convenience wrapper for Segment"""
        ai = segment.parameters['ai']
        m = segment.parameters['m']
        m0 = segment.genome_medians['m0']      
        try:
            segment.parameters.update (self.score_dipl(ai,m,m0,models))
        except IndexError:
            segment.parameters.update ({'model' : 'UN', 'd_model' : np.nan,
                                        'k': np.nan, 'p_model' : np.nan, 'score_HE' : np.nan,
                                        'd_HE' : np.nan})

    
def fit_QQgauss (values, fit_intercept = True):
    x = sts.norm.ppf (np.linspace (0,1,len(values)+2)[1:-1])
    huber = slm.HuberRegressor (fit_intercept = fit_intercept)
    huber.fit (x[:, np.newaxis], np.sort(values))
    return {'m' : huber.intercept_, 's' : huber.coef_[0]}
    
def fit_smallest_gauss (values):
    current_values = np.sort(values)
    thr = current_values.max()+1
    previous_thr = thr + 1    
    
    
    while (previous_thr > thr):
        current_values = current_values[current_values < thr]
        
        popt, pcov = opt.curve_fit (sts.norm.cdf, current_values, np.linspace (0,1,len(current_values)),
                                    p0 = [np.mean (current_values), np.std (current_values)])
        previous_thr = thr
        thr = sts.norm.ppf (1-1/(5*len(current_values)), *popt)
        
    return {'m': popt[0], 's' : popt[1], 'thr' : thr, 'alpha' : 1/(5*len(current_values))}