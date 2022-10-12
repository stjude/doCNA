from dataclasses import dataclass
import pandas as pd
import numpy as np
import multiprocessing as mpl
import scipy.stats as sts
import scipy.optimize as opt
from sklearn.linear_model import HuberRegressor

from doCNA import Testing
from doCNA import Chromosome
from doCNA import Run
from doCNA import Consts
from doCNA import Report


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
            self.logger.info ("{len(self.data)} SNPs to analyze. No filtering applied.")

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
                              alpha = Consts.DSCORE_ALPHA, k_thr = Consts.K_THR)                
        
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
        
    def score_clonality (self, size_thr = 5e6, model_thr = 3, alpha = 0.01, k_thr = 0.11):
        balanced = [seg.parameters['model'] == 'A(AB)B' for seg in self.all_segments]
        big = [(seg.end - seg.start)/1e6 > size_thr for seg in self.all_segments]
        notHO = [seg.parameters['k'] < k_thr for seg in self.all_segments]
        fit_model = [seg.parameters['model_score'] < model_thr for seg in self.all_segments]
        
        all_data = np.array([(seg.parameters['k'], (seg.end - seg.start)/1e6) for seg in self.all_segments])
                
        balanced_index = np.where ([ba&bi&fi&nh for ba,bi,fi,nh in zip(balanced, big, fit_model, notHO)])[0]
        imbalanced_index = np.where ([(~ba)&bi&fi&nh for ba,bi,fi,nh in zip(balanced, big, fit_model, notHO)])[0]
        ed = {'A' : np.nan, 'B' : np.nan, 'C' : np.nan, 'down' : np.nan, 
                  'up' : np.nan, 'm' : np.nan, 's' : np.nan}
        try:
            self.genome_medians['clonality_balanced'] = fit_huber (all_data[balanced_index,:],
                                                               alpha)
            if self.genome_medians['clonality_balanced']['A'] < 0:
                self.logger.warning ("Scoring of balanced segments seems to fail. Check before you yell!")
        except:
            self.logger.warning ("Scoring of balanced segments failed. None of the scoring may sense.")    
            self.genome_medians['clonality_balanced'] = ed
            
        try:            
            self.genome_medians['clonality_imbalanced'] = fit_huber (all_data[imbalanced_index,:],
                                                                 alpha)
            if self.genome_medians['clonality_imbalanced']['A'] < 0:
                self.logger.warning ("Scoring of imbalanced segments seems to fail. Check before you yell!")
        
        except:
            self.logger.warning ("Scoring of imbalanced segments failed. None of the scoring may sense.")    
            self.genome_medians['clonality_imbalanced'] = ed
            
            
        A = (self.genome_medians['clonality_balanced']['A'], self.genome_medians['clonality_imbalanced']['A'])
        B = (self.genome_medians['clonality_balanced']['B'], self.genome_medians['clonality_imbalanced']['B'])
        C = (self.genome_medians['clonality_balanced']['C'], self.genome_medians['clonality_imbalanced']['C'])
        m = (self.genome_medians['clonality_balanced']['m'], self.genome_medians['clonality_imbalanced']['m'])
        s = (self.genome_medians['clonality_balanced']['s'], self.genome_medians['clonality_imbalanced']['s'])
        up = (self.genome_medians['clonality_balanced']['up'], self.genome_medians['clonality_imbalanced']['up'])
        down = (self.genome_medians['clonality_balanced']['down'], self.genome_medians['clonality_imbalanced']['down'])
        
        
        
        
        i = 0
        self.logger.info ('Score for balanced segments:')
        self.logger.info (f'Core usuallness: log(k) = {-A[i]} log(s) + {-C[i]}')
        self.logger.info (f'Normal estimation of distance to usual: m  = {m[i]}, s = {s[i]}.')
        self.logger.info (f'Estimated normal range of distance to usual: from {down[i]} to {up[i]}.')
        
        i = 1
        self.logger.info ('Score for imbalanced segments:')
        self.logger.info (f'Core usuallness: log(k) = {-A[i]} log(s) + {-C[i]}')
        self.logger.info (f'Normal estimation of distance to usual: m  = {m[i]}, s = {s[i]}.')
        self.logger.info (f'Estimated normal range of distance to usual: from {down[i]} to {up[i]}.')
        
        for seg in self.all_segments:
            x = np.log10((seg.end - seg.start)/10**6)
            y = np.log10(seg.parameters['k'])
            i = 0 if seg.parameters['model'] == 'A(AB)B' else 1
            seg.parameters['k_d'] = (A[i]*x+B[i]*y+C[i])/np.sqrt (A[i]**2+B[i]**2)
            seg.parameters['clonality_score'] = -np.log10(sts.norm.sf(seg.parameters['k_d'], m[i], s[i]))
            seg.parameters['call'] = 'CNV' if seg.parameters['k_d'] > up[i] else 'norm'
            
            
            
            
                   
    def get_clonality_cnB_params (self, percentiles = (1,80)):

        ks = []
        ns = []
        a = self.genome_medians['model_d']['a']
        for chrom in self.chromosomes.keys():
            for seg in self.chromosomes[chrom].segments:
                if seg.parameters['model'] == 'A(AB)B':
                    score = -np.log10 (np.exp (-a*seg.parameters['d']))
                    size = seg.end - seg.start
                    if (size > Consts.SIZE_THR)&(score < Consts.MODEL_THR):
                        ks.append (seg.parameters['k'])
                        ns.append (seg.parameters['n'])
        try:
            z = np.array(ks)*np.sqrt(np.array(ns))
            
            pp = np.percentile (z[~np.isnan(z)], percentiles)
            zz = z[(z >= pp[0])&(z <= pp[1])]
            res, _ = opt.curve_fit (sts.norm.cdf, np.sort(zz), np.linspace (0,1,len(zz)), p0 = [np.mean(zz), np.std(zz)])
            self.logger.info (f'Clonality distribution (for cnB model): FI(k) = G({res[0]}, {res[1]}))')
            down, up = Testing.get_outliers_thrdist (zz, alpha = 0.005)
            self.logger.info (f'Estimated normal range of distance to usual: from {down} to {up}.')
        except (IndexError,TypeError, RuntimeError, ValueError):
            res = (np.nan, np.nan)
            self.logger.warning ('Automatic estimation of clonality distribution for cnB model failed.')
            down = np.nan
            up = np.nan
            
        return {'m' : res[0], 's' : res[1], 'down_thr' : down, 'up_thr' : up}

    def get_k_params (self):
        
        a = self.genome_medians['model_d']['a']
        ks = []
        ss = []
        for chrom in self.chromosomes.keys():
            for seg in self.chromosomes[chrom].segments:
                size = (seg.end - seg.start)/10**6
                score = -np.log10 (np.exp (-a*seg.parameters['d']))

                big_filter = (seg.parameters['k'] < 0.11)|(size < 1)
                num_filter = (~np.isnan(seg.parameters['k']))& (~np.isinf(seg.parameters['k']))  
                score_filter = (~np.isinf(score))&(~np.isnan(score))&(score < Consts.MODEL_THR)
                if (seg.parameters['k']>0)&score_filter&(seg.parameters['model'] != 'A(AB)B')&num_filter&(seg.centromere_fraction < 0.1)&big_filter:
                    ks.append (seg.parameters['k'])
                    ss.append (size)
        k = np.log10 (np.array(ks))
        s = np.log10 (np.array(ss))
        try:
            huber = HuberRegressor(alpha = 0.0, epsilon = 1.35)
            huber.fit(s[:, np.newaxis], k)

            A = -huber.coef_[0]
            B = 1
            C = -huber.intercept_
            d = (A*s+B*k+C)/np.sqrt (A**2+B**2)

            down, up = Testing.get_outliers_thrdist (d, alpha = 0.005)
            m, std = sts.norm.fit (d[(d > down)&(d < up)])
            self.logger.info (f'Core usuallness: log(k) = {-A} log(s) + {-C}')
            self.logger.info (f'Normal estimation of distance to usual: m  = {m}, s = {std}.')
            self.logger.info (f'Estimated normal range of distance to usual: from {down} to {up}.')
        except ValueError:
            A = np.nan
            B = np.nan
            C = np.nan
            down = np.nan
            up = np.nan
            m = np.nan
            std = np.nan
        return {'A' : A, 'B' : B, 'C' : C, 'down' : down, 'up' : up, 'm' : m, 'std' : std}
        
    def report (self, report_type = 'bed'):
        return Report.Report(report_type).genome_report(self)
    
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
    return {'A' : A, 'B' : B, 'C' : C, 'down' : down, 'up' : up, 'm' : m, 's' : std}


def f (c):
    c.find_runs()
    c.generate_segments ()
    return c    

def exp (x,a):
    return 1 - np.exp(-a*x)
    
    
