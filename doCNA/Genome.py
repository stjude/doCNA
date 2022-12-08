from dataclasses import dataclass
import pandas as pd
import numpy as np
import multiprocessing as mpl
import scipy.stats as sts
import scipy.optimize as opt
from sklearn.linear_model import HuberRegressor

from doCNA import Testing
from doCNA import Chromosome
from doCNA import Distribution
from doCNA import Consts
from doCNA import Report
from doCNA import Run


class Genome:
    """Class to run genome wide tests of HE and create chromosomes."""
    def __init__(self, sample_name, logger, config, CB_file, no_processes = 1):
        self.sample_name = sample_name
        self.no_processes = no_processes
        self.config = config
        self.CB = pd.read_csv (CB_file, sep = '\t')
        self.logger = logger.getChild (f'{self.__class__.__name__}')
        self.genome_medians = {}
        self.logger.debug ("Object genome created.")
    
    def retrive_counts_create_chromosomes (self, data_file, columns, SG_file = None):
        """Reads the data in and filters through SuperGood list, if not None"""
        alldata = pd.read_csv (data_file, sep = '\t', usecols = columns)[columns]
        alldata.columns = ['chrom','position','ref_count', 'alt_count', 'Type']
        alldata['chrom'] = np.where(alldata['chrom'].str.contains("chr") == False, "chr"+alldata['chrom'], alldata['chrom'])
        
        self.logger.debug ('Data retrived.') 
        if SG_file is not None:
            SGlist = pd.read_csv (SG_file, compression = 'gzip',
                                header = None, sep = '\t', names = ['chrom','position','Ref', 'Alt', 'Status'],
                                dtype = {'chrom' : str, 'position' : np.int32, 'Ref' : str, 'Alt' : str, 'Status' : str})
            self.data = alldata.loc [alldata['Type'] == 'SNP'].merge (SGlist[['chrom','position']], how = 'inner')
            del (SGlist)
            del (alldata)
            self.logger.info (f"Input file filtered through SG list: {len(self.data)} SNPs to analyze.")
            
        else:
            self.data = alldata
            self.logger.info (f"{len(self.data)} SNPs to analyze. No filtering applied.")

        self.data['cov'] = self.data['ref_count'] + self.data['alt_count']
        self.data['vaf'] = self.data['alt_count']/self.data['cov']
        
        self.chromosomes = {}
        self.logger.debug ("Creating chromosomes...")
        
        for chrom, data in self.data.loc[~self.data['vaf'].isna()].groupby (by = 'chrom'):
            if chrom not in Consts.SEX_CHROMS:
                self.chromosomes[chrom] = Chromosome.Chromosome (chrom, data.copy(), 
                                                                 self.config, self.logger,
                                                                 self.genome_medians, 
                                                                 self.CB.loc[self.CB['chrom'] == chrom])
                self.logger.debug (f"Chromosome {chrom} has {len(data)} markers.")
         

            
    def segment_genome (self, m0 = 0, fb_alpha = Consts.FB_ALPHA):
        
        self.logger.debug ('Starting testing ...')
        
        self.HE = Testing.Testing ('HE', 
                                   self.chromosomes,
                                   self.logger)
        
        self.HE.run_test(no_processes = self.no_processes)

        self.logger.debug ("Genomewide heterozygosity: " + f"\n" + str(self.HE.report_results()))
                
        outliers = self.HE.results.loc[self.HE.results['chi2'] > float(self.config['HE']['max_chi2'])].index.tolist()
        outliers_fraction = len(outliers)/len(self.HE.results)
        if len(outliers):
            self.logger.info (f'HE outliers due to chi2: {outliers}.')
        if outliers_fraction == 1:
            self.logger.critical ('Sample failed HE model. All chromosomes chi2 above threshold.')
            exit (1)
        elif outliers_fraction > 0.5:
            self.logger.warning ('Half or more of the chromosomes above threshold. May be inaccurate.')
               
        self.HE.analyze (parameters = self.config['HE'], outliers = outliers, skip_par = ['cov'])
        self.logger.info ("Genome heterozygosity reference: "+ f"\n" + str(self.HE.get_genome_medians()))
        
        self.genome_medians['HE'] = self.HE.get_genome_medians()
        
        self.logger.debug ('Running N/E marking.')
        
        for chrom in self.chromosomes.keys():
            status = self.HE.get_status (chrom)
            params = self.HE.get_parameters(chrom).copy()
            for  par in status.index.values:
                if ~status[par]:
                    params[par] = self.HE.get_genome_medians()[par]  
            
            self.chromosomes[chrom].markE_onHE (params, float(self.config['HE']['z_thr']))
                    
        self.logger.debug ('Testing N/E marking.')    
        
        self.VAF = Testing.Testing ('VAF', self.chromosomes, self.logger)
        self.VAF.run_test (self.HE.medians['cov'], no_processes = self.no_processes)
        
        self.logger.debug ("Genomewide VAF: " + f"\n" + str(self.VAF.report_results()))
        
        
        VAFS_chi2 = self.VAF.results['chi2']
        thr_chi2 = float(self.config['VAF']['max_chi2'])

        outliers = self.VAF.results.loc[(VAFS_chi2 > thr_chi2)|(VAFS_chi2 == 0)].index.tolist()
        outliers_fraction = len(outliers)/len(self.HE.results)
        if len(outliers):
            self.logger.info (f'VAF outliers due to chi2: {outliers}.')
        if outliers_fraction == 1:
            self.logger.critical ('All chromosoms failed VAF model. All chromosomes chi2 above threshold.')
            exit (1)
        elif outliers_fraction > 0.5:
            self.logger.warning ('Half or more of the chromosomes above threshold. May be inaccurate.')
        self.VAF.analyze (parameters = self.config['VAF'], outliers = outliers)
        
        self.logger.info ("Genome VAF reference: "+ f"\n" + str(self.VAF.get_genome_medians()))
        self.genome_medians['VAF'] = self.VAF.get_genome_medians()
        
        for chrom in self.chromosomes.keys():
            status = self.VAF.get_status (chrom)
            self.logger.debug (f'Chromosome {chrom} inlier: {status}')
            if ~status.T.all(axis = 0):
                params = self.VAF.get_parameters (chrom)
                self.logger.debug (f'Chromosome {chrom} marked on full model.')
                self.chromosomes[chrom].mark_on_full_model (self.HE.medians['cov']) 
                
        self.COV = Testing.Testing ('COV', self.chromosomes, self.logger)
        self.COV.run_test(no_processes = self.no_processes,
                          exclude_symbol = [Consts.N_SYMBOL, Consts.U_SYMBOL])
        self.logger.debug ('Genomewide COV: ' + f"\n" + str(self.COV.report_results()))
        
        self.COV.analyze (parameters = self.config['COV'], outliers = self.VAF.get_outliers())
        
        self.logger.info ("Genome COV reference: " + f"\n" + str(self.COV.get_genome_medians()))
        self.genome_medians['COV'] = self.COV.get_genome_medians()

        inliers = self.VAF.get_inliers()
        inliers_fb = self.VAF.results.loc[inliers,'fb'].values
        
        if len (np.unique(inliers_fb)) < 4:
            self.genome_medians['fb'] = np.percentile (inliers_fb,1-fb_alpha)
            self.logger.warning(f'Widening parameter estimation based on normal approximation not possible. {1-fb_alpha} percentile used.')
        else:
            try:
                self.genome_medians['fb'] = Testing.get_outliers_thrdist (np.unique(inliers_fb), fb_alpha, r = 0.5)[1]
            except:
                self.logger.exception ('Estimation of fb failed.')
                exit (1)
  
        if m0 > 0:
            self.logger.info (f"Using user supplied m0 = {m0}, instead of estimated m0 = {self.genome_medians['COV']['m']}")
            self.genome_medians['m'] = m0
        else:
            self.genome_medians['m'] = self.genome_medians['COV']['m'] 

        if self.COV.medians['m'] < float(self.config['COV']['min_cov']):
            self.logger.critical (f"Coverage is below threshold {self.COV.medians['m']} < {self.config['COV']['min_cov']}")
            exit (1)

        self.logger.debug ("Starting segmentation.")
        if self.no_processes > 1:            
            with mpl.Pool (processes = self.no_processes) as pool:
                segmented_chroms = pool.map (f, self.chromosomes.values())
                for sc in segmented_chroms:
                    self.chromosomes [sc.name] = sc
        else:
            for chrom in self.chromosomes.keys():
                self.chromosomes[chrom].find_runs()
                self.chromosomes[chrom].generate_segments ()
        self.logger.debug ("Segmentation finished.")
        
        self.logger.info (f"Median coverage for the sample: m = {str(self.genome_medians['COV']['m'])}")
        
        self.logger.info ("Scoring segments.")
        
        self.all_segments = []
        for chrom in self.chromosomes.keys():
            for seg in self.chromosomes[chrom].segments:
                self.all_segments.append (seg)
                
        
        self.score_model_distance ()
        self.score_clonality (size_thr = Consts.SIZE_THR, model_thr = Consts.MODEL_THR,
                              dalpha = Consts.DSCORE_ALPHA, kalpha = Consts.KSCORE_ALPHA,
                              k_thr = Consts.K_THR)

        
    def score_model_distance (self):
    
        zs_ns = [(seg.parameters['d'], seg.parameters['n']) for seg in self.all_segments]
        
        z_n_a = np.array(zs_ns)
        z_n = z_n_a[~np.isnan(z_n_a[:,1]) ,:]
        try:
            popt, _ = opt.curve_fit (exp, np.sort (z_n[:,0]), np.linspace (0,1,len(z_n[:,0])),
                                     p0 = (10), sigma = 1/np.sqrt(z_n[:,1])[np.argsort(z_n[:,0])])
            self.logger.info ('Distance from model /d/ distribution: FI(d) = exp(-{:.5f} d)'.format (popt[0]))
            
        except ValueError:
            popt = [np.nan]
            self.logger.warning ("Scoring of models failed. None of the scoring may sense.")
            self.logger.warning ("Consider rerunning with manually set m0.")

        for seg in self.all_segments:
            seg.parameters['model_score'] = -np.log10 (np.exp (-popt[0]*seg.parameters['d']))
        self.genome_medians['model_d'] = {'a' : popt[0]}
        
    def score_clonality (self, size_thr = 5e6, model_thr = 3, dalpha = 0.01, kalpha = 0.01, k_thr = 0.11):
        balanced = [seg.parameters['model'] == 'A(AB)B' for seg in self.all_segments]
        big = [(seg.end - seg.start)/1e6 > size_thr for seg in self.all_segments]
        notHO = [seg.parameters['k'] < k_thr for seg in self.all_segments]
        fit_model = [seg.parameters['model_score'] < model_thr for seg in self.all_segments]
        
        all_data = np.array([(seg.parameters['k'], (seg.end - seg.start)/1e6) for seg in self.all_segments])
                
        balanced_index = np.where ([ba&bi&fi for ba,bi,fi in zip(balanced, big, fit_model)])[0]
        imbalanced_index = np.where ([(~ba)&bi&fi&nh for ba,bi,fi,nh in zip(balanced, big, fit_model, notHO)])[0]
        ed = {'A' : np.nan, 'B' : np.nan, 'C' : np.nan, 'down' : np.nan, 
                  'up' : np.nan, 'm' : np.nan, 's' : np.nan, 'score_FDR' : np.inf}
            
        try:            
            self.genome_medians['clonality_imbalanced'] = fit_huber (all_data[imbalanced_index,:],
                                                                     dalpha)
            
            if self.genome_medians['clonality_imbalanced']['A'] < 0:
                self.logger.warning ("Scoring of imbalanced segments seems to fail. Check before you yell!")
            
        except:
            self.logger.warning ("Scoring of imbalanced segments failed. None of the scoring makes sense.")    
            self.genome_medians['clonality_imbalanced'] = ed
            
        A = self.genome_medians['clonality_imbalanced']['A']
        B = self.genome_medians['clonality_imbalanced']['B']
        C = self.genome_medians['clonality_imbalanced']['C']
        m = self.genome_medians['clonality_imbalanced']['m']
        s = self.genome_medians['clonality_imbalanced']['s']
        up = self.genome_medians['clonality_imbalanced']['up']
        down = self.genome_medians['clonality_imbalanced']['down']
        
        imbalanced_all_index = np.where ([(~ba)&bi&fi for ba,bi,fi in zip(balanced, big, fit_model)])[0] 
        data = all_data[imbalanced_all_index,:]
        ks = np.log10 (data[:,0])
        ss = np.log10 (data[:,1])
        d = (A*ss+B*ks+C)/np.sqrt (A**2+B**2)
        
        self.genome_medians["clonality_imbalanced"]["score_FDR"] = FDR(np.sort(sts.norm.sf (d, m, s)), dalpha)

        self.logger.info ('Score for imbalanced segments:')
        self.logger.info (f'Core usuallness: log(k) = {-A} log(s) + {-C}')
        self.logger.info (f'Normal estimation of distance to usual: m  = {m}, s = {s}.')
        self.logger.info (f'Estimated normal range of distance to usual: from {down} to {up}.')
        self.logger.info (f'FDR corrected score threshold: {self.genome_medians["clonality_imbalanced"]["score_FDR"]}.')
        
        k = all_data[balanced_index,0]
        try: 
            if len(k) < Consts.MIN_LEN_K_BALANCED:
                self.logger.warning (f'Only {len(k)} balance regions.')
                raise ValueError 
            
            gauss = Distribution.fit_single_G(k, alpha = kalpha, r = 0.2)
            if gauss['p'] > 0.3:
                params = gauss['2']
                params['thr'] = gauss['thr']
                params['p'] = gauss['p']
                params['a'] = 1 
            else:
                bounds = [[0,-0.2, 0.01, 0.0, 0.01],
                          [1, 0.0, 0.2, 0.2, 0.2]]
            
                gauss = Distribution.fit_double_G (k, alpha = kalpha, r = 0.2, 
                                                   initial_bounds = bounds,
                                                   initial_p0 = (0.5, -0.1, 0.02, 0.1, 0.02))
            
               
                params = gauss['2']
                params['thr'] = gauss['thr']
                params['p'] = gauss['p']
                

            self.genome_medians['clonality_balanced'] = params
            print (self.genome_medians['clonality_balanced'])
            self.genome_medians['clonality_balanced']['score_FDR'] = FDR (score_double_gauss (k[:,np.newaxis],
                                                                                          params['m'][np.newaxis,:],
                                                                                          params['s'][np.newaxis,:]),
                                                                                          kalpha )

            self.logger.info ('Score for balanced segments:')
            if params['m'][0] != params['m'][1]:
                self.logger.info (f'Estimation fits double normal as: p = {params["p"]}.')
                self.logger.info (f'Double normal estimation of k: m  = {params["m"]}, s = {params["s"]}.')
            else:
                self.logger.info (f'Estimation fits normal as: p = {params["p"]}.')
                self.genome_medians['clonality_balanced']['m'] = params['m'][:1]
                self.genome_medians['clonality_balanced']['s'] = params['s'][:1]
                self.logger.info (f'Single normal estimation of k: m  = {params["m"]}, s = {params["s"]}.')
                
            self.logger.info ('Note on quality')
            self.logger.info ('Distance of balanced distributions to k = 0:')
            self.logger.info (f'Absolute: m = {params["m"]}')
            self.logger.info (f'Relative: z = {params["m"]/params["s"]}')
            self.logger.info (f'FDR corrected score threshold: {self.genome_medians["clonality_balanced"]["score_FDR"]}.')
        except:
            self.logger.exception ('Scoring of balanced regions failed and it makes no sense.')
            params = {'p' : np.nan,
                      'm' : np.array([np.nan, np.nan]),
                      's' : np.array([np.nan, np.nan]),
                      'a' : np.array([np.nan, np.nan]),
                      'thr' : np.array([np.nan, np.nan]),
                      'score_FDR' : np.inf}

        self.genome_medians['clonality_balanced'] = params

        for seg in self.all_segments:
            x = np.log10((seg.end - seg.start)/10**6)
            y = np.log10(seg.parameters['k'])
            if seg.parameters['model'] != 'A(AB)B':
                seg.parameters['k_d'] = (A*x+B*y+C)/np.sqrt (A**2+B**2)
                seg.parameters['clonality_score'] = -np.log10(sts.norm.sf(seg.parameters['k_d'], m, s))
                seg.parameters['call'] = 'CNVi' if seg.parameters['k_d'] > up else 'norm'
                seg.parameters['call_FDR'] = 'CNVi' if seg.parameters['clonality_score'] > self.genome_medians["clonality_imbalanced"]["score_FDR"] else 'norm'
            else:
                k = seg.parameters['k']
                z = (k - params['m'])/params['s']
                p = np.min((sts.norm.cdf (z[0]), sts.norm.sf (z[-1])))
                seg.parameters['k_d'] = np.nan
                seg.parameters['clonality_score'] = -np.log10(p)
                seg.parameters['call'] = 'norm' if (k < params['thr'][1]) & (k > params['thr'][0]) else 'CNVb'
                seg.parameters['call_FDR'] = 'CNVb' if seg.parameters['clonality_score'] > self.genome_medians['clonality_balanced']['score_FDR'] else 'norm'
                   
            
    def report (self, report_type = 'bed'):
        return Report.Report(report_type).genome_report(self)

#k, m, s must have correct dimentions
def score_double_gauss (k, m, s):
    z = (k - m)/s
    ps =np.concatenate((sts.norm.cdf (z[:,0])[:, np.newaxis], sts.norm.sf (z[:,1])[:,np.newaxis]), axis = 1)
    p = np.min(ps , axis = -1)
    return  np.sort(p)

def fit_huber (data, alpha):
    k = np.log10 (data[:,0])
    s = np.log10 (data[:,1])
    huber = HuberRegressor(alpha = 0.0, epsilon = 1.35)
    huber.fit(s[:, np.newaxis], k)

    A = -huber.coef_[0]
    B = 1
    C = -huber.intercept_
    d = (A*s+B*k+C)/np.sqrt (A**2+B**2)
    
    down, up = Testing.get_outliers_thrdist (d, alpha = alpha)
    m, std = sts.norm.fit (d[(d > down)&(d < up)])
    
    score_FDR = FDR (np.sort(sts.norm.sf (d, m, std)), alpha)

    return {'A' : A, 'B' : B, 'C' : C, 'down' : down, 'up' : up, 'm' : m,
            's' : std, 'score_FDR' : score_FDR}

def FDR (p, alpha):
    k = np.arange (1, len(p)+1)
    index = np.where(p <= alpha*k/len(p))[0]
    try:
        return -np.log10(p[np.max(index)])
    except:
        return np.inf

def f (c):
    c.find_runs()
    c.generate_segments ()
    return c    

def exp (x,a):
    return 1 - np.exp(-a*x)
    
    
