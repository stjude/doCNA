from dataclasses import dataclass
import pandas as pd
import numpy as np
import multiprocessing as mpl
import scipy.stats as sts
import scipy.optimize as opt
import sys
import sklearn.linear_model as slm
from sklearn.linear_model import HuberRegressor

from doCNA import Testing
from doCNA import Chromosome
from doCNA import Consts
from doCNA import Report
from doCNA import Scoring


class Genome:
    """Class to run genome wide tests of HE and create chromosomes."""
    def __init__(self, sample_name, logger, config, CB_file, models, no_processes = 1):
        self.sample_name = sample_name
        self.no_processes = no_processes
        self.config = config
        self.models = models
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
        self.sex_chromosomes = {}
        
        if len(self.data.loc[~self.data['vaf'].isna()]) == 0:
            sys.exit (f'No data read in. Please cross check column names: {columns} with input file.')
            
        self.logger.debug ("Creating chromosomes...")
        
        for chrom, data in self.data.loc[~self.data['vaf'].isna()].groupby (by = 'chrom'):
            if chrom not in Consts.SEX_CHROMS:
                self.chromosomes[chrom] = Chromosome.Chromosome (chrom, data.copy(), 
                                                             self.config, self.logger,
                                                             self.genome_medians, 
                                                             self.CB.loc[self.CB['chrom'] == chrom])
            else:
                self.sex_chromosomes[chrom] = Chromosome.Chromosome (chrom, data.copy(), 
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
        
        ####To Be Tested
        #self.HEn = Testing.Testing ('HEn',
        #                            self.chromosomes,
        #                            self.logger)
        #self.HEn.run_test(no_processes = self.no_processes)

        #self.logger.debug ("Genomewide heterozygosity: " + f"\n" + str(self.HEn.report_results()))
        ####
        
        self.logger.debug ('Running N/E marking.')
        
        
        for chrom in self.chromosomes.keys():
            status = self.HE.get_status (chrom)
            params = self.HE.get_parameters(chrom).copy()
            for  par in status.index.values:
                if ~status[par]:
                    params[par] = self.HE.get_genome_medians()[par]  
            
            self.chromosomes[chrom].markE_onHE (params, float(self.config['HE']['z_thr']))
            
        #for formatting
        for  par in status.index.values:
            params[par] = self.HE.get_genome_medians()[par]  
        
        for chrom in self.sex_chromosomes.keys():
            self.sex_chromosomes[chrom].markE_onHE (params, float(self.config['HE']['z_thr']))
        
        self.logger.debug ('Testing N/E marking.')    
        
        self.VAF = Testing.Testing ('VAF', 
                                    self.chromosomes,
                                    self.logger)
        
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
            #self.logger.debug (f'Chromosome {chrom} inlier: {status}')
            if ~status.T.all(axis = 0):
                params = self.VAF.get_parameters (chrom)
                self.chromosomes[chrom].mark_on_full_model (self.HE.medians['cov']) 
                self.logger.debug (f'Chromosome {chrom} marked on full model.')
        
        for chrom in self.sex_chromosomes.keys():
            self.logger.info (f'Analyzing {chrom}: ...')
            chrom_vaf_results = Testing.VAF_test (self.sex_chromosomes[chrom].data,
                                                     self.HE.medians['cov'])
            self.logger.info (f'VAF test results: {chrom_vaf_results}')
            self.VAF.results.append(pd.DataFrame.from_records ([chrom_vaf_results],
                                                                columns = chrom_vaf_results._fields, 
                                                                index = [chrom]))
        
            self.chromosomes[chrom] = self.sex_chromosomes[chrom]
            self.logger.info(f'Chromosome {chrom} added.')        
        
        self.COV = Testing.Testing ('COV', 
                                    self.chromosomes,
                                    self.logger)
        self.COV.run_test(no_processes = self.no_processes,
                          exclude_symbol = [Consts.N_SYMBOL, Consts.U_SYMBOL])
        self.logger.debug ('Genomewide COV: ' + f"\n" + str(self.COV.report_results()))
        
        self.COV.analyze (parameters = self.config['COV'], outliers = self.VAF.get_outliers()+Consts.SEX_CHROMS)
        
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
  
        self.genome_medians['v0'] = self.genome_medians['VAF']['vaf']
        if m0 > 0:
            self.logger.info (f"Using user supplied m0 = {m0}, instead of estimated m0 = {self.genome_medians['COV']['m']}")
            self.genome_medians['m0'] = m0
        else:
            self.genome_medians['m0'] = self.genome_medians['COV']['m'] 

        self.genome_medians['l'] = self.genome_medians['COV']['l']
        
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
                
        self.segment_filter = [((s.end - s.start)/1e6 > Consts.SIZE_THR) &\
                                (s.centromere_fraction < Consts.CENTROMERE_THR) &\
                                (s.parameters['ai'] < Consts.DIPLOID_AI_THR) &\
                                (np.abs(s.parameters['m']/self.genome_medians['m0']-1) < Consts.DIPLOID_dCN_THR/2)    for s in self.all_segments] 

        
        self.segments = [self.all_segments[i] for i in np.where(self.segment_filter)[0]]
        
        data_for_scoring = np.array([(s.parameters['ai'], 2*s.parameters['m']/self.genome_medians['m0']-2, s.parameters['n']) for s in self.segments])
        self.scorer = Scoring.Scoring (fb = self.genome_medians['fb'], m0 = self.genome_medians['m0'], window_size = Consts.SNPS_IN_WINDOW, 
                                       initial_data = data_for_scoring, logger = self.logger)
        
        ps = np.zeros (len(self.all_segments))
        for i, seg in enumerate (self.all_segments):
            self.scorer.score_dipl(seg)
            ps[i] = seg.parameters['p_HE']
        
        ##FDRing threshold
        thr = FDR(np.sort(ps[np.isfinite(ps)]), alpha = Consts.DIPLOID_ALPHA, score = True)
        self.genome_medians['thr_HE'] = thr
        self.scorer.set_thr (thr)
        for seg in self.all_segments:
            self.scorer.analyze_segment(seg, self.models)
        
        self.score_model_distance ()
            

    def score_model_distance (self):
    
        size_filter = np.array([(seg.end - seg.start)/1e6 > Consts.SIZE_THR  for seg in self.all_segments])
        cent_filter = np.array([seg.centromere_fraction < Consts.CENTROMERE_THR  for seg in self.all_segments])
        model_filter = np.array([seg.parameters['model'] not in ['AB', '(AB)(2-n)', '(AB)(2+n)'] for seg in self.all_segments])
        finite_filter = np.array([np.isfinite(seg.parameters['d_model']) for seg in self.all_segments])
        
        filter = size_filter & cent_filter & model_filter & finite_filter
        
        if sum(filter) >= 3:
            indexes = np.where(filter)[0]
            segments = [self.all_segments[i] for i in indexes]
            zs_ns = [seg.parameters['d_model'] for seg in segments]

    
            z_n = np.array(zs_ns)
            x = sts.expon.ppf (np.linspace (0,1,len(z_n)+2)[1:-1])
            huber = slm.HuberRegressor (fit_intercept = False)
            huber.fit (x[:, np.newaxis], np.sort(z_n))
            a = -1./huber.coef_[0]
            self.logger.info ('Distance from model /d/ distribution: FI(d) = exp({:.5f} d)'.format (a))
            ps = []
            for seg in self.all_segments:
                if seg.parameters['model'] != 'AB':
                    p = np.exp (a*seg.parameters['d_model'])
                    seg.parameters['score_model'] = -np.log10 (p)
                    ps.append(p)
                else:
                    seg.parameters['score_model'] = 0
                
            self.genome_medians['d_model'] = {'a' : a}
            ps = np.array(ps)
            self.genome_medians['thr_model'] = FDR (np.sort(ps[np.isfinite(ps)]), Consts.MODEL_APLHA, score = True)
        else:
            self.logger.info ('Not enough non diploid regions to perform meaningful scoring')
            self.genome_medians['d_model'] = {'a' : np.nan}
            for seg in self.all_segments:
                seg.parameters['score_model'] = 0
            self.genome_medians['thr_model'] = np.nan
            
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
    huber = slm.HuberRegressor(alpha = 0.0, epsilon = 1.35)
    huber.fit(s[:, np.newaxis], k)

    A = -huber.coef_[0]
    B = 1
    C = -huber.intercept_
    d = (A*s+B*k+C)/np.sqrt (A**2+B**2)
    
    
    down, up = Testing.get_outliers_thrdist (d, alpha = alpha)
    inlier_ds = d[(d > down)&(d < up)]
    m, std = sts.norm.fit (inlier_ds)
    std = std/Bolch_correction (len(inlier_ds))

    score_FDR = FDR (np.sort(sts.norm.sf (d, m, std)), alpha, score = True)

    return {'A' : A, 'B' : B, 'C' : C, 'down' : down, 'up' : up, 'm' : m,
            's' : std, 'score_FDR' : score_FDR}

def Bolch_correction (n):
    return 1 - 1/(4*n) - 7/(32*n**2) - 19/(128*n**3)


def FDR (p, alpha, score = False):
    k = np.arange (1, len(p)+1)
    index = np.where(p <= alpha*k/len(p))[0]
    try:
        if score:
            return -np.log10(p[np.max(index)])
        else:
            return p[np.max(index)]
    except:
        return np.inf

def f (c):
    c.find_runs()
    c.generate_segments ()
    return c    

def exp (x,a):
    return 1 - np.exp(-a*x)
    
    
