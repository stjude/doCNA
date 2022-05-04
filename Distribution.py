import numpy as np
import scipy.optimize as opt
import scipy.stats as sts
from doCNA.Testing import get_outliers_thrdist


LENGTH_THRESHOLD = 10

class Distribution:
    """
    Class design to describe the data with either single or double Gauss distributions with managing outliers.

    Atributes:
        parameters (dic): Parameters of estimated distribution(s)
        all_parameters (dic): Parameters estimated for both situations, if done
        key (string): Resolved to {single/double} distribution
    """
    
    def __init__ (self, values, p_thr = 0.3, thr_z = 1.0):
        """
        Class constructor. Estimates distribution(s)' parameters based on provided values.

        """
        assert len (values) > LENGTH_THRESHOLD, print ("Not enough data points to consider distribution.")
        single_G_par = fit_single_G (np.sort(values), alpha = 0.01, r = 0.5)
        self.all_parameters = {}
        self.all_parameters['single'] = single_G_par
                
        z = np.abs(values-single_G_par['m'][0])/single_G_par['s'][0]
        string = list ('O' * len(values))
        for i in np.where ((z < thr_z))[0]:
            string[i] = 'B'
        
        self.all_parameters['single']['string'] = string
                
        if single_G_par['p'] < p_thr:
            double_G_par = fit_double_G (np.sort(values), alpha = 0.01, r = 0.5)
            self.key = 'double'
            self.parameters = double_G_par
            self.all_parameters ['double'] = double_G_par
            z0 = np.abs(values-double_G_par['m'][0])/double_G_par['s'][0]
            z1 = np.abs(values-double_G_par['m'][1])/double_G_par['s'][1]
            #generate string
            string = list ('O' * len(values))
            #marking C population
            for i in np.where ((z0 < z1)&(z0 < thr_z))[0]:
                string[i] = 'C'
            #marking D population
            for i in np.where ((z1 < z0)&(z1 < thr_z))[0]:
                string[i] = 'D'
            self.parameters['string'] = ''.join(string)
            self.all_parameters['double']['string'] = string
        else:
            self.key = 'single'
            self.parameters = single_G_par
            
            self.string = ''.join(string)
        
        
    def fail_normal (self):
        return self.key == 'double'

    #To remove after test
    def combinations (self):
        if self.key == 'single':
            comb = [(0,0),]
        else:
            comb = [(0,1),(1,0)]
        return comb
        
    def combinations_of_params (self, dim = 1, key = 'single', reverse = False):
        """Method to generate parameters in desired shape"""
        
        if key == 'double':
            assert key == self.key, "No double distribution!"
            m = self.all_parameters['double']['m']
            s = self.all_parameters['double']['s']
        elif key == 'single':
            if dim == 1:
                m = self.all_parameters['single']['m'][0]
                s = self.all_parameters['single']['s'][0]
            elif dim == 2:
                #ugly
                m = np.array([self.all_parameters['single']['m'][0],self.all_parameters['single']['m'][0]])
                s = np.array([self.all_parameters['single']['s'][0],self.all_parameters['single']['s'][0]])
            else:
                raise (f'Dim can be only 1 or 2. Dim = {dim} does not make much sense.')
                
        if reverse:
            return m[:-1:-1], s[:-1:-1]
        else:
            return m, s
    
    def to_string (self):
        return (self.parameters)
        
    def __str__ (self):
        return (self.parameters)
    
    def __repr__ (self):
        return (self.parameters)
    
    
def fit_single_G (values, alpha, r):
    #outliers
    thr = get_outliers_thrdist (np.sort(values))
    #fit
    a = np.sort(values[(values >= thr[0])&(values <= thr[1])])
    popt, pcov = opt.curve_fit (sts.norm.cdf, a, np.linspace(0,1, len(a)), 
                                p0 = [np.mean(a), np.std (a)])
    ksp = sts.kstest (a, sts.norm.cdf, args = popt)

    #report
    return {'p' : ksp.pvalue, 'm': np.array([popt[0]]),
            's': np.array([popt[1]]), 'thr' : thr, 'a' : np.ones(1)}

def fit_double_G (values_all, alpha, r):
    #p0 = (0.5, np.median(values[:int(len(values)/2)]), np.percentile(values,40)-np.percentile(values,10),
    #      np.median(values[int(len(values)/2):]), np.percentile(values,90)-np.percentile(values,60))
    thr0 = get_outliers_thrdist (np.sort(values_all))
    
    values = values_all[(values_all >= thr0[0]) & (values_all <= thr0[1])]
    
    p0 = (0.5, np.percentile (values, 25), np.percentile(values,40)-np.percentile(values,10),
          np.percentile (values, 75), np.percentile(values,90)-np.percentile(values,60))
        
    
    popti, pcovi = opt.curve_fit (gaus2, np.sort(values), np.linspace (0,1,len(values)), p0 = p0,
                                  bounds = ((0, np.percentile (values, 5), 0.1*np.percentile (values, 40)-0.1*np.percentile (values, 10),
                                                np.percentile (values, 55), 0.1*np.percentile (values, 90)-0.1*np.percentile (values, 60)),
                                            (1, np.percentile (values, 45), 2*np.percentile (values, 40)-2*np.percentile (values, 10),
                                                np.percentile (values, 95), 2*np.percentile (values, 90)-2*np.percentile (values, 60)))) 

    if popti[3] > popti[1]:
        out_max = sts.norm.ppf (1-alpha, popti[3], popti[4])
        out_min = sts.norm.ppf (alpha, popti[1], popti[2])
    else:
        out_max = sts.norm.ppf (1-alpha, popti[1], popti[2])
        out_min = sts.norm.ppf (alpha, popti[3], popti[4])
        tmp = popti
        popti = (tmp[0], tmp[3], tmp[4], tmp[1], tmp[2])


    a = values[(values <= out_max)&(values >= out_min)]
    p0 = (0.5, np.percentile (a, 25), np.percentile(a,40)-np.percentile(a,10),
          np.percentile (a, 75), np.percentile(a,90)-np.percentile(a,60))
    
    popt, pcov = opt.curve_fit (gaus2, np.sort(a), np.linspace (0,1,len(a)), p0 = p0,
                                  bounds = ((0, np.percentile (a, 5), 0.1*np.percentile (a, 40)-0.1*np.percentile (a, 10),
                                                np.percentile (a, 45), 0.1*np.percentile (a, 90)-0.1*np.percentile (a, 60)),
                                            (1, np.percentile (a, 55), 2*np.percentile (a, 40)-2*np.percentile (a, 10),
                                                np.percentile (a, 95), 2.*np.percentile (a, 90)-2*np.percentile (a, 60)))) 
    
    a0,m0,s0,m1,s1 = popt
    da0,dm0,ds0,dm1,ds1 = np.sqrt(np.diag(pcov))
    ksp = sts.kstest (a, gaus2, args = popt)

    return {'p' : ksp.pvalue, 'm': np.array([m0,m1]), 
            's': np.array([s0,s1]), 'a' : np.array([a0, 1-a0]),
            'thr' : (out_min, out_max)}

def gaus2 (v, a, m0, s0, m1, s1):
    return a*sts.norm.cdf (v, m0,s0) + (1-a)*sts.norm.cdf(v, m1, s1)

