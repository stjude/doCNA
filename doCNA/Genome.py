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
from doCNA.Report import Report

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
        #self.COV = Testing.Testing ('COV', self.chromosomes, self.logger)
        #self.COV.run_test(no_processes = self.no_processes)
        #self.COV.analyze (parameters = self.config['COV'])
                 
        #if self.COV.medians['m'] < float(self.config['COV']['min_cov']):
        #    self.logger.critical (f"Coverage is below threshold {self.COV.medians['m']} < {self.config['COV']['min_cov']}")
        #    exit (0)
        
        #self.logger.debug ("Genomewide coverage: " + f"\n" + str(self.COV.report_results()))
        #self.logger.info ("Genome coverage medians: "+ f"\n" + str(self.COV.get_genome_medians()))
        
        #self.genome_medians['COV'] = self.COV.get_genome_medians()     
        #if m0 > 0:
        #    self.logger.info (f"Using user supplied m0 = {m0}, instead of estimated m0 = {self.genome_medians['COV']['m']}")
        #    self.genome_medians['m'] = m0
        #else:
        #    self.genome_medians['m'] = self.genome_medians['COV']['m']
        
        #self.HE = Testing.Testing ('HE', 
        #                           self.chromosomes,
        #                           self.logger)
        
        #self.HE.run_test(no_processes = self.no_processes)

        #self.logger.debug ("Genomewide heterozygosity: " + f"\n" + str(self.HE.report_results()))
                
        #outliers = self.HE.results.loc[self.HE.results['chi2'] > float(self.config['HE']['max_chi2'])].index.tolist()
        #outliers_fraction = len(outliers)/len(self.HE.results)
        
        #if outliers_fraction == 1:
        #    self.logger.critical ('Sample failed HE model. All chromosomes chi2 above threshold.')
        #    exit (1)
        #elif outliers_fraction > 0.5:
        #    self.logger.warning ('Half or more of the chromosomes above threshold. May be inaccurate.')
        
                     
        #self.HE.analyze (parameters = self.config['HE'], outliers = outliers)
       
        
        #self.logger.info ("Genome heterozygosity medians: "+ f"\n" + str(self.HE.get_genome_medians()))
        
        #self.genome_medians['HE'] = self.HE.get_genome_medians()        
                
        #self.logger.debug ('First round of N/E marking.')
        
        #for chrom in self.chromosomes.keys():
        #    status = self.HE.get_status (chrom)
        #    if status:
        #        self.chromosomes[chrom].markE_onHE (self.HE.get_parameters(chrom),
        #                                            float(self.config['HE']['z_thr']))
        #    else:
        #        self.chromosomes[chrom].markE_onHE (self.HE.get_genome_medians(),
        #                                            float(self.config['HE']['z_thr']))
        
        #self.logger.debug ('Testing first round of N/E marking.')    
        
        #self.VAF = Testing.Testing ('VAF', self.chromosomes, self.logger)
        #self.VAF.run_test (self.COV.medians['m'], no_processes = self.no_processes)
        #self.logger.debug ("Genomewide VAF: " + f"\n" + str(self.VAF.report_results()))
        #self.VAF.analyze (parameters = self.config['VAF'])
        
        #self.logger.info ("Genome VAF medians: "+ f"\n" + str(self.VAF.get_genome_medians()))
        
        #self.genome_medians['VAF'] = self.VAF.get_genome_medians()
        
        #inliers_fb = self.VAF.results.loc[self.VAF.get_inliers(), 'fb'].values
         
        #if m0 > 0:
        #    self.logger.info (f"Using user supplied m0 = {m0}, instead of estimated m0 = {self.genome_medians['COV']['m']}")
        #    self.genome_medians['m'] = m0
        #else:
        #    self.genome_medians['m'] = self.COV.results.loc[[c in self.VAF.get_inliers() for c in self.COV.results.index.values],'m'].median()

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
            #print (status)
            #print (params)
            for  par in status.index.values:
                if ~status[par]:
                    params[par] = self.HE.get_genome_medians()[par]  
            #if status:
            self.chromosomes[chrom].markE_onHE (params, #self.HE.get_parameters(chrom),
                                                    float(self.config['HE']['z_thr']))
            #else:
            #    self.chromosomes[chrom].markE_onHE (self.HE.get_genome_medians(),
            #                                        float(self.config['HE']['z_thr']))
        
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
        
        self.COV.analyze (parameters = self.config['COV'], outliers = self.VAF.get_outliers())#outliers)
        
        self.logger.info ("Genome COV reference: " + f"\n" + str(self.COV.get_genome_medians()))
        self.genome_medians['COV'] = self.COV.get_genome_medians()

        inliers = self.VAF.get_inliers()
        #print ('inliers: ', inliers)
        inliers_fb = self.VAF.results.loc[inliers,'fb'].values
        #print ('fb: ', inliers_fb)
        
        if len (np.unique(inliers_fb)) < 4:
            self.genome_medians['fb'] = np.percentile (inliers_fb,1-fb_alpha)
            self.logger.warning(f'Widening parameter estimation based on normal approximation not possible. {1-fb_alpha} percentile used.')
        else:
            try:
                self.genome_medians['fb'] = Testing.get_outliers_thrdist (np.unique(inliers_fb), fb_alpha, r = 0.5)[1]
            except:
                self.logger.critical ('Estimation of fb failed.')
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
        
        self.logger.debug ("Scoring segments.")
                
        self.genome_medians['model_d'] = self.get_distance_params ()
        
        self.genome_medians['ai'] = self.get_ai_params ()

        self.genome_medians['k'] = self.get_k_params ()
        
        self.genome_medians['clonality_cnB'] = self.get_clonality_cnB_params ()
        

    def get_clonality_cnB_params (self, percentiles = (1,80)):

        ks = []
        ns = []
        a = self.genome_medians['model_d']['a']
        for chrom in self.chromosomes.keys():
            for seg in self.chromosomes[chrom].segments:
                if seg.parameters['model'] == 'A(AB)B':
                    score = -np.log10 (np.exp (-a*seg.parameters['d']))
                    size = seg.end - seg.start
                    #Consts.SIZE_THR
                    if (size > Consts.SIZE_THR)&(score < Consts.MODEL_THR):
                        ks.append (seg.parameters['k'])
                        ns.append (seg.parameters['n'])
        try:
            z = np.array(ks)*np.sqrt(np.array(ns))
            
            pp = np.percentile (z[~np.isnan(z)], percentiles)
            zz = z[(z >= pp[0])&(z <= pp[1])]
            print (zz)
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
    
    
                
    def get_distance_params (self, percentiles = (10,80)):

        zs = []
        ns = []
        for chrom in self.chromosomes.keys():
            for seg in self.chromosomes[chrom].segments:
                if (seg.symbol == Consts.E_SYMBOL)&(~np.isnan(seg.parameters['d'])):
                    zs.append (seg.parameters['d'])
                    ns.append (seg.parameters['n'])        
        z = np.array(zs)
        s = 1/np.sqrt(np.array(ns))
        try:
            popt, _ = opt.curve_fit (exp, z, np.linspace (0,1,len(z)), p0 = (10), sigma = s)
            self.logger.info ('Distance from model /d/ distribution: FI(d) = exp(-{:.5f} d)'.format (popt[0]))
        except ValueError:
            popt = [np.nan]

        return {'a' : popt[0]}

    def get_k_params (self):
        
        a = self.genome_medians['model_d']['a']
        ks = []
        ss = []
        for chrom in self.chromosomes.keys():
            for seg in self.chromosomes[chrom].segments:
                size = (seg.end - seg.start)/10**6
                #if np.exp (-a*seg.parameters['d']) == 0:
                    #print (seg)
                    #print (-np.log10 (np.exp (-a*seg.parameters['d'])))
                score = -np.log10 (np.exp (-a*seg.parameters['d']))

                big_filter = (seg.parameters['k'] < 0.11)|(size < 1)
                num_filter = (~np.isnan(seg.parameters['k']))& (~np.isinf(seg.parameters['k']))  
                score_filter = (~np.isinf(score))&(~np.isnan(score))&(score < Consts.MODEL_THR)
                if (seg.parameters['k']>0)&score_filter&(seg.parameters['model'] != 'A(AB)B')&num_filter&(seg.centromere_fraction < 0.1)&big_filter:
                    ks.append (seg.parameters['k'])
                    ss.append (size)
        #print (ks)
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
        return {'A' : A, 'B' : B, 'C' : C, 'down_thr' : down, 'up_thr' : up, 'm' : m, 'std' : std}
        
    def get_ai_params (self, percentile = 50):
        zs = []
        ns = []
        for chrom in self.chromosomes.keys():
            for seg in self.chromosomes[chrom].segments:
                size = seg.end - seg.start
                if (size > Consts.SIZE_THR)&(seg.symbol == Consts.E_SYMBOL):
                    zs.append (seg.parameters['ai'])
                    ns.append (seg.parameters['n'])
        z = np.array(zs)
        s = 1/np.sqrt(np.array(ns))
        zsf = z*s
        zs = np.sort(zsf[~np.isnan(zsf)]) 
        try:
            popt, _ = opt.curve_fit (exp, zs, np.linspace (0,1,len(zs)), p0 = (10))
            self.logger.info ('AI distribution (for non-cnB models): FI(ai) = exp(-{:.5f} ai)'.format (popt[0]))
        except ValueError:
            popt = [np.nan]
            
        return {'a' : popt[0]}

    def report (self, report_type = 'bed'):
        return Report(report_type).genome_report(self)
    
def f (c):
    c.find_runs()
    c.generate_segments ()
    return c    

def exp (x,a):
    return 1 - np.exp(-a*x)
    
    