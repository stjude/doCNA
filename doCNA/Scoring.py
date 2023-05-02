import numpy as np
import scipy.stats as sts
import scipy.optimize as opt
import sklearn.linear_model as slm

from doCNA import Consts
from doCNA import Models


class Scoring:
    def __init__(self, initial_data = None, logger = False) -> None:
        
        if initial_data is None:
            self.ai_param = {'m' : 0, 's' : 1}
            self.cn_param = {'m' : 0, 's' : 1}
            self.dipl_dist = {'m' : 0, 's' : 1, 'thr' : 0, 'alpha': np.nan}
            self.median_size = 1
        else:
            self.median_size = np.median(initial_data[:, 2])            
            self.ai_param = fit_QQgauss(initial_data[: ,0])
            self.cn_param = fit_QQgauss(initial_data[: ,1], fit_intercept = False)
            dds =  initial_data - np.array([self.ai_param['m'], self.cn_param['m']])[np.newaxis, :]
            ds = dds/np.array([self.ai_param['s'],self.cn_param['s']])[np.newaxis, :]
            self.dipl_dist = fit_QQgauss (np.sqrt((ds**2).sum(axis = 1)))
            self.dipl_dist['alpha'] = Consts.SCORE_ALPHA
            self.dipl_dist['thr'] = sts.norm.ppf (1-self.dipl_dist['alpha'], 
                                                  self.dipl_dist['m'],
                                                  self.dipl_dist['s'])
        
        if logger:
            self.logger = logger.getChild (f'{self.__class__.__name__}')
            self.logger.info (f"Median segment size: {self.median_size}")
            self.logger.info (f"Distribution of distance to diploid (0,2): m = {self.dipl_dist['m']}, s = {self.dipl_dist['s']}")
            self.logger.info (f"Distribution of diploid allelic imbalance: m = {self.ai_param['m']}, s = {self.ai_param['s']}")
            self.logger.info (f"Distribution of diploid copy number: m = {self.cn_param['m']}, s = {self.cn_param['s']}")
            
    def get_d_thr (self):
        return self.dipl_dist['thr']
    
    def score_dipl (self, segment): #ai, m, m0, models):
        ai = segment.parameters['ai'] 
        m = segment.parameters['m']
        m0 = segment.genome_medians['m0']
        cn = 2*m/m0
        scale = np.sqrt (self.median_size / segment['n'])
        m_ai = self.ai_param['m']
        s_ai = self.ai_param['s']*scale
        m_cn = self.cn_param['m']
        s_cn = self.cn_param['s']*scale
        d = np.sqrt (((ai-m_ai)/s_ai)**2 + (((cn-2)-m_cn)/s_cn)**2)
        p_d = sts.norm.sf (d, self.dipl_dist['m'], self.dipl_dist['s'])
        
        segment.parameters['d_HE'] = d# np.sqrt ((ai/s_ai)**2 + ((cn-2)/s_cn)**2)
        segment.parameters['p_HE'] = p_d
        segment.parameters['score_HE'] = -np.log10(p_d)
        


    def analyze_segment (self, segment, models):
        m = segment.parameters['m']
        m0 = segment.genome_medians['m0']
        cn = 2*m/m0
        
        if segment.parameters['score_HE'] > self.dipl_dist['thr']:
            ai = segment.parameters['ai']
              
            try:
                segment.parameters.update (Models.pick_model(ai,1, cn,1,models))
            except IndexError:
                segment.parameters.update ({'model' : 'UN', 'd_model' : np.nan,
                                            'k': np.nan, 'p_model' : np.nan,})
        else:
            segment.parameters.update({'model' : 'AB', 'd_model' : np.nan, 'k': cn/2-1})
                
    def set_thr (self, thr):
        self.dipl_dist['thr'] = thr

    
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