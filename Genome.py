import pandas as pd
import numpy as np
import multiprocessing as mpl
import scipy.stats as sts
import scipy.optimize as opt
#from Chromosome import E_SYMBOL

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
        for chrom, data in self.data.groupby (by = 'chrom'):
            if chrom not in Consts.SEX_CHROMS:
                self.chromosomes[chrom] = Chromosome.Chromosome (chrom, data.copy(), 
                                                                 self.config, self.logger,
                                                                 self.genome_medians, 
                                                                 self.CB.loc[self.CB['chrom'] == chrom])
            
    def segment_genome (self, fb_alpha = Consts.FB_ALPHA):
        
        self.logger.debug ('Starting testing ...')
        self.COV = Testing.Testing ('COV', self.chromosomes, self.logger)
        self.COV.run_test(no_processes = self.no_processes)
        self.COV.analyze (parameters = self.config['COV'])
                
        if self.COV.medians['m'] < float(self.config['COV']['min_cov']):
            self.logger.critical (f"Coverage is below threshold {self.COV.medians['m']} < {self.config['COV']['min_cov']}")
            exit (0)
        
        self.logger.debug ("Genomewide coverage: " + f"\n" + str(self.COV.report_results()))
        self.logger.info ("Genome coverage medians: "+ f"\n" + str(self.COV.get_genome_medians()))
                 
        self.genome_medians['COV'] = self.COV.get_genome_medians()
        
        self.HE = Testing.Testing ('HE', 
                                   self.chromosomes,
                                   self.logger)
        
        self.HE.run_test(no_processes = self.no_processes)
        self.HE.analyze (parameters = self.config['HE'])
       
        self.logger.debug ("Genomewide heterozygosity: " + f"\n" + str(self.HE.report_results()))
        self.logger.info ("Genome heterozygosity medians: "+ f"\n" + str(self.HE.get_genome_medians()))
        
        self.genome_medians['HE'] = self.HE.get_genome_medians()        
              
        self.logger.debug ('First round of N/E marking.')
        for chrom in self.chromosomes.keys():
            self.chromosomes[chrom].markE_onHE (self.HE.get_parameters(chrom),
                                                float(self.config['HE']['z_thr']))
        
        self.logger.debug ('Testing first round of N/E marking.')    
        
        self.VAF = Testing.Testing ('VAF', self.chromosomes, self.logger)
        self.VAF.run_test (self.COV.medians['m'], no_processes = self.no_processes)
        self.VAF.analyze (parameters = self.config['VAF'])
        
        self.logger.debug ("Genomewide VAF: " + f"\n" + str(self.VAF.report_results()))
        self.logger.info ("Genome VAF medians: "+ f"\n" + str(self.VAF.get_genome_medians()))
        
        self.genome_medians['VAF'] = self.VAF.get_genome_medians()
        
        self.genome_medians['fb'] = Testing.get_outliers_thrdist (self.VAF.results.loc[self.VAF.get_inliers(), 'fb'],
                                                                  alpha = fb_alpha, r = 0.5)[1]
        self.logger.info (f"Widening parameters estimated at: {self.genome_medians['fb']}")

       
        for chrom in self.chromosomes.keys():
            status = self.VAF.get_status (chrom)
            self.logger.debug (f'Chromosome {chrom} inlier: {status}')
            if ~status:
                params = self.VAF.get_parameters (chrom)
                #if params['chi2'] >= float(self.config['VAF']['chi2_high']):
                self.logger.debug (f'Chromosome {chrom} marked on full model.')
                self.chromosomes[chrom].mark_on_full_model (self.COV.medians['m'])

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
        
        self.genome_medians['clonality_cnB'] = self.get_clonality_cnB_params ()
        

    def get_clonality_cnB_params (self, percentiles = (1,80)):

        ks = []
        for chrom in self.chromosomes.keys():
            for seg in self.chromosomes[chrom].segments:
                if seg.parameters['model'] == 'cnB':
                    a = seg.genome_medians['model_d']['a']
                    score = -np.log10 (np.exp (-a*seg.parameters['d']))
                    if (score < 3)&(~np.isnan(seg.parameters['k'])):
                        ks.append (seg.parameters['k'])
        z = np.array(ks)
        pp = np.percentile (z, percentiles)
        zz = z[(z >= pp[0])&(z <= pp[1])]
        res, _ = opt.curve_fit (sts.norm.cdf, np.sort(zz), np.linspace (0,1,len(zz)), p0 = [0.01, 0.01])
        self.logger.info (f'Clonality distribution (for cnB model): FI(k) = G({res[0]}, {res[1]}))')
        return {'m' : res[0], 's' : res[1]}
    
    
                
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
        
        popt, _ = opt.curve_fit (exp, z, np.linspace (0,1,len(z)), p0 = (10), sigma = s)
        self.logger.info ('Distance from model /d/ distribution: FI(d) = exp(-{:.5f} d)'.format (popt[0]))
        return {'a' : popt[0]}
        
    def get_ai_params (self, percentile = 50):
        zs = []
        ns = []
        for chrom in self.chromosomes.keys():
            for seg in self.chromosomes[chrom].segments:
                if seg.symbol == Consts.E_SYMBOL:
                    zs.append (seg.parameters['ai'])
                    ns.append (seg.parameters['n'])
        z = np.array(zs)
        s = np.sqrt(np.array(ns))

        popt, _ = opt.curve_fit (exp, z, np.linspace (0,1,len(z)), p0 = (10))
        self.logger.info ('AI distribution (for non-cnB models): FI(ai) = exp(-{:.5f} ai)'.format (popt[0]))
        return {'a' : popt[0]}

    def report (self, report_type = 'bed'):
        return Report(report_type).genome_report(self.chromosomes)
    
def f (c):
    c.find_runs()
    c.generate_segments ()
    return c    

def exp (x,a):
    return 1 - np.exp(-a*x)
    
    
