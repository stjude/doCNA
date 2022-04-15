import scipy.optimize as opt
import numpy as np
import multiprocessing as mpl
import pandas as pd
from math import ceil,floor
import scipy.stats as sts
from collections import namedtuple


COV_INITIAL_SHAPE = 0.14
COV_SHAPE_RANGE = (-2,1)

HE_COV_PERC_BOUNDS = (0.05, 99.0)
HE_VAF_BOUNDS = (0.4,0.6)
HE_FCOV_BOUNDS = (0.01, 0.8)
HE_FN_BOUNDS = (0.2,1.0)
HE_A_BOUNDS = (0.1,0.9)
HE_B_BOUNDS = (1,10)
HE_LERR_BOUNDS = (2,10)

VAF_VAF_BOUNDS = (0.45,0.55)
VAF_N_THR = 100

FB_F_MAX = 1.4
FB_EPS = 1e-4

COV_results = namedtuple ('COV_results', ['m', 'dm', 'l', 'dl'])
VAF_results = namedtuple ('VAF_results', ['chi2', 'vaf', 'fb'])
HE_results = namedtuple ('HE_results', ['chi2', 'vaf', 'cov', 'b'])

E_SYMBOL = 'E'

class Testing:
    def __init__  (self, test_name, chromosomes, logger):
        assert test_name in ['COV', 'HE', 'VAF'], f"Unknown test: {test_name}!"
        
        i = np.where([t == test_name for t in ['COV', 'HE', 'VAF']])[0][0]
        self.test = [COV_test, HE_test, VAF_test][i]
        self.chromosomes = chromosomes
        self.logger = logger.getChild (f'{self.__class__.__name__}-{self.test.__name__}')
        
    def run_test (self, no_processes = 1, *args):
                
        if no_processes > 1:
            with mpl.Pool (processes = no_processes) as pool:
                results = pool.starmap (self.test, [(self.chromosomes[chrom].data, args) for chrom in self.chromosomes.keys()])
        else:
            results = []
            for chrom in self.chromosomes.keys():
                results.append((self.test (self.chromosomes[chrom].data, args)))
        self.results = pd.DataFrame.from_records (results, columns = results[0]._fields, 
                                                  index = self.chromosomes.keys())
        
    def analyze (self, parameters = {'alpha' : 0.01, 'r' : 0.5}, q = (5,99)):
        #parameters are not done well here
        self.normal_range = {}
        columns = self.results.columns
        status = {}
        
        for column in columns:
            try:
                alpha = parameters[column + '_alpha']
                r = parameters[column + '_r']
            except:
                alpha, r = (0.01, 0.5)
                
            res = self.results[column].values
            self.results[column + '_status'] = 'outlier'
            try:
                range = get_outliers_thrdist (self.results[column].values, alpha, r)
            except:
                self.logger.warning ('Range estimation of {} failed. Using percentiles.'.format (column))
                range = np.percentile (res, q = q)
            
            self.logger.info ('Estimated normal ragne of {} is from {} to {}'.format (column, *range))
            in_or_out = (self.results[column] > range[0]) & (self.results[column] < range[1])
            status[column] = ['inlier' if iou else 'outlier' for iou in in_or_out]
            
            self.normal_range[column] = range
        
        self.status = pd.DataFrame.from_dict (status)
        self.status['chrom'] = self.results.index
        self.status.set_index ('chrom', inplace = True)
        inlier_filter = (self.status == 'inlier').all(axis = 1) 
        self.medians = self.results.loc[inlier_filter, ].median (axis = 0, numeric_only = True)

    def get_parameters (self, chromosome):
        #get chromosome parameters 
        try:
            test_results = self.results['chromosome']
        #or return medians
        except:
            test_results = self.medians
        return test_results
    
    def get_status (self, chromosome):
        try:
            status = self.status[chromosome]
        except:
            status = 'NA'
        return status
    
    def get_genome_medians (self):
        return self.medians
    
    def get_inliers (self):
        return self.status[(self.status == 'inlier').all(axis = 1)].index.values
    
    def get_outliers (self):
        return self.status[(self.status == 'outlier').all(axis = 1)].index.values
    
def COV_test (data, *args, **kwargs):
    
    initial_shape = COV_INITIAL_SHAPE if 'initial_shape' not in kwargs else kwargs['initial_shape']
    shape_range = COV_SHAPE_RANGE if 'shape_range' not in kwargs else kwargs['shape_range']
        
    
    percentiles = np.arange (0.01, 0.99, 0.01)
    covs = data['cov'].values   
    y = np.percentile (covs, percentiles*100)
    cov_range = (floor(y[0]*0.5) , ceil (y[-1]*1.05))
    initial_m = np.median(covs)
             
    #here may be good place for try except
    popt, pcov = opt.curve_fit (lambda_cdf, percentiles, y, p0 = [initial_m, initial_shape],
                                bounds = [(cov_range[0],shape_range[0]),
                                          (cov_range[1],shape_range[1])])
    m, l  = popt
    dm, dl  = np.sqrt(np.diag(pcov)) 
    
    return COV_results (m = m, dm = dm, l = l, dl = dl)
    
def Q (p,l):
    if l == 0:
        return np.log(p/(1-p))
    else:
        return (p**l-(1-p)**l)/l    

def lambda_cdf (p, m, lam):
    return Q (p, lam)*np.sqrt(m)+m

def HE_test (data, *args, **kwargs):
    
    cov_perc_bounds = HE_COV_PERC_BOUNDS if 'cov_perc_bounds' not in kwargs else kwargs['cov_perc_bounds']
    vaf_bounds = HE_VAF_BOUNDS if 'vaf_bounds' not in kwargs else kwargs['vaf_bounds']
    fcov_bounds = HE_FCOV_BOUNDS if 'fcov_bounds' not in kwargs else kwargs['fcov_bounds']
    fN_bounds = HE_FN_BOUNDS if 'fN_bounds' not in kwargs else kwargs['fN_bounds']
    a_bounds = HE_A_BOUNDS if 'a_bounds' not in kwargs else kwargs['a_bounds']
    b_bounds = HE_B_BOUNDS if 'b_bounds' not in kwargs else kwargs['b_bounds']
    lerr_bounds = HE_LERR_BOUNDS if 'lerr_bounds' not in kwargs else kwargs['lerr_bounds']
        
    
    def chi2 (params, counts, N):
        vaf, fcov, fN, a, b, l = params
        fe = 10**(-l)
        cs = np.arange (0, cov_max +1)
        cov = cov_min + fcov*(cov_max-cov_min)
        ns = 2*fN*N*cn2_cov_pdf (cs, cov, b)
        chi2 = 0
        
        for c, cnt, in counts:
            i = np.arange(0,c+1)
            nhe = ns[c]*cn2_vaf_pdf (i/c,vaf,c)
            nho = ns[c]*HO_vaf_pdf (i, c, fe ,b)
        
            na = a*nhe + (1-a)*nho
            chi2 += sum((cnt - na)**2/np.sqrt(na*na+1))/c 
        return chi2/len(counts)

    cov_min, cov_max = np.percentile (data['cov'].values, q = cov_perc_bounds)
   
    counts = []
    for c in np.arange (int(cov_min), int(cov_max+1)):
        d = data.loc[data['cov'] == c]
        h = np.histogram(d['alt_count'].values, bins = np.arange(0,c+2)-0.5)
        counts.append ((c, h[0]))
        del (h)
    
    N = len(data)
    fcov = (data['cov'].median() - cov_min)/(cov_max - cov_min)
    
    res = opt.minimize (chi2, x0 = (0.5, fcov, 0.5,0.8, 1.3, 6), args = (counts, N),
                    bounds = (vaf_bounds, fcov_bounds, fN_bounds, a_bounds, b_bounds, lerr_bounds))
    
    vaf, fcov, fN, a, b, l = res.x
    
    cov = cov_min + fcov*(cov_max-cov_min)
    
    return HE_results(chi2 = res.fun, vaf = vaf, cov = cov, b = b)

def cn2_vaf_pdf (x,v,c):
    p = sts.norm.pdf (x, v, np.sqrt((v*(1-v))/c))
    return p/sum(p)

def cn2_cov_pdf (n,c,b = 1):
    return sts.norm.pdf (n, c, np.sqrt(b*c))

def HO_vaf_pdf (i, n, fe = 10**-6,b=1):
    return sts.binom.pmf(i, n, 1-b*fe)

def VAF_test (data, m, **kwargs):
    
    #assert m in kwargs, 'Median coverage needed for VAF test'
    #m = kwargs['m']
    vaf_bounds = VAF_VAF_BOUNDS if 'vaf_bounds' not in kwargs else kwargs['vaf_bounds'] 
    n_thr = VAF_N_THR if 'n_thr' not in kwargs else kwargs['n_thr']
        
    def chi2 (v, counts):
        chi2 = 0
        cc = 0
        for c, h, n in counts:
            if n >= n_thr:
                cc += c + 1 
                na = n*sts.binom.pmf (np.arange(c+1),c,v)
                chi2 += sum((h-na)**2/np.sqrt(na*na+1))
        
        if cc == 0: 
            cc = 1        
        return chi2/cc
    
    
    counts = []
    for c,d in data.loc[data['symbol'] == E_SYMBOL].groupby (by = 'cov'):
        h = np.histogram(d['alt_count'].values, bins = np.arange(0,c+2)-0.5)
        n = len(d)
        counts.append ((c, h[0], n))
        del (h)
            
    res = opt.minimize (chi2, x0 = (0.5) , args = (counts), method = 'L-BFGS-B',
                        bounds = (vaf_bounds,))        
    
    fb = find_fb (np.sort(data.loc[data['symbol'] == E_SYMBOL, 'vaf'].values), m[0], f_max = FB_F_MAX, eps = FB_EPS)    
    
    return VAF_results (chi2 = res.fun, vaf = res.x[0], fb = fb)

def find_fb (vafs, m, f_max = FB_F_MAX, eps = FB_EPS):
    
    def two_gauss (v, dv, a):
        v0 = 0.5
        return a*sts.norm.cdf (v, v0 + dv, smin)+(1-a)*sts.norm.cdf (v,v0 - dv,smin)
    
    s = np.sqrt (0.25/m)
    f = 1.0
    df = 0.02
    smin = s*f
    cdf = np.linspace (0,1, len(vafs))
    popt, _ = opt.curve_fit (two_gauss, vafs, cdf, p0 = [0.01, 0.5],
                            bounds = ((0, 0.0),
                                      (1, 1)))
    a_prev = popt[1]
    #a_prev = 0.5
    
    f = f + df
    popt, _ = opt.curve_fit (two_gauss, vafs, cdf, p0 = [0.01, 0.5],
                            bounds = ((0, 0.0),
                                      (1, 1)))
    a_next = popt[1]

    while (a_next > eps)&(f < f_max):
        #new f
        f += df
        smin = s*f
        popt, _ = opt.curve_fit (two_gauss, vafs, cdf, p0 = popt,
                            bounds = ((0, 0.0),
                                      (1, 1)))
        a_prev = a_next
        a_next = popt[1]
    
    a_zero = a_next
    a_nzero = a_prev
    
    df = -df/2
    
    while (a_nzero > 2*eps)&(a_nzero - a_zero > eps)&(np.abs(df) > eps)&(f < f_max):
        f = f + df
        smin = s*f
        popt, _ = opt.curve_fit (two_gauss, vafs, cdf, p0 = [0.01, 0.5],
                            bounds = ((0.0, 0.0),
                                      (1.0, 1.0)))
        if popt[1] > eps:
            a_nzero = popt[1]
            df = np.abs(df/2)
        else:
            df = -np.abs(df/2)
        
    if f < f_max:
         ff = f
    else:
        print (f)
        ff = np.nan
    
    return ff

def get_outliers_thrdist (values, alpha = 0.01, r = 0.5):
    n = len (values)
    current = values
    core = values
    mc = np.mean (current)
    sc = np.std (current)
    for i in range (int(np.floor(len(values)*r))):
        m = np.mean (current)
        s = np.std (current)
        R = np.abs ((current -m)/s)
        j = np.where (R == R.max())[0][0]
        
        ni = n - i
        p = 1 - alpha/(ni-1)
        t = sts.t.ppf (p, ni-1)
        l = (ni*t)/np.sqrt((ni - 1 +t*t)*(ni+1))
        
        current = np.concatenate ([current[:j], current[j+1:]])
        if R.max() > l:
            core = current.copy()
            mc = np.mean(current)
            sc = np.std (current)
            
    #return (core.min(), core.max()),(sts.norm.ppf (alpha, mc,sc), sts.norm.ppf (1-alpha, mc,sc)), (mc,sc)
    return (sts.norm.ppf (alpha, mc,sc), sts.norm.ppf (1-alpha, mc,sc))