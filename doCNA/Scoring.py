import numpy as np
import scipy.stats as sts
import scipy.optimize as opt
import sklearn.linear_model as slm

from doCNA import Consts
from doCNA import Models


class Scoring:
    def __init__(self, initial_data, logger, diploid_ai_thr = 0.1) -> None:
        
        self.logger = logger.getChild (f'{self.__class__.__name__}-{self.name}')
        self.ai_param = {}
        self.logger.info ('')
        self.cn_param = {}
        self.logger.info ('')
        self.dipl_dist = {}
        self.logger.info ('')
        
    def score_n_ktype (self, segments) -> None:
        pass
    
    
    
def fit_QQgauss (values, fit_intercept = True):
    x = sts.norm.ppf (np.linspace (0,1,len(values)+2)[1:-1])
    huber = slm.HuberRegressor (fit_intercept = True)
    huber.fit (x[:, np.newaxis], np.sort(values))
    return {'m' : huber.intercept_, 's' : huber.coef_[0]}
    
def fit_smalles_gauss (values):
    current_values = np.sort(values)
    thr = current_values.max()+1
    previous_thr = thr + 1    
    while (previous_thr > thr) | ():
        current_values = current_values[current_values < thr]
        popt, pcov = opt.curve_fit (sts.norm.cdf, current_values, np.linspace (0,1,len(current_values)))
        previous_thr = thr
        thr = sts.norm.ppf (1-1/(5*len(current_values)), *popt)
    
    return {'m': popt[0], 's' : popt[1], 'thr' : thr, 'alpha' : 1/(5*len(current_values))}