import pandas as pd
import numpy as np
import multiprocessing as mpl
import scipy.stats as sts
#from Chromosome import E_SYMBOL

from doCNA import Testing
from doCNA import Chromosome
from doCNA import Run
from doCNA.Report import Report

SEX_CHROMS = ['chrX', 'chrY']

class Genome:
    """Class to run genome wide tests of HE and create chromosomes."""
    def __init__(self, sample_name, logger, config, CB_file = None, no_processes = 1):
        self.sample_name = sample_name
        self.no_processes = no_processes
        self.config = config
        self.CB_file = CB_file
        self.logger = logger.getChild (f'{self.__class__.__name__}')
        self.genome_medians = {}
        self.logger.debug ("Object genome created.")
    
    def retrive_counts_create_chromosomes (self, data_file, columns, SG_file = None):
        """Reads the data in and filters through SuperGood list, if not None"""
        alldata = pd.read_csv (data_file, sep = '\t', usecols = columns)[columns]
        alldata.columns = ['chrom','position','ref_count', 'alt_count', 'Type']
        alldata['chrom'] = np.where(alldata['chrom'].str.contains("chr") == False, "chr"+alldata['chrom'], alldata['chrom'])

        if SG_file is not None:
            SGlist = pd.read_csv (SG_file, compression = 'gzip',
                                header = None, sep = '\t', names = ['chrom','position','Ref', 'Alt', 'Status'],
                                dtype = {'chrom' : str, 'position' : np.int32, 'Ref' : str, 'Alt' : str, 'Status' : str})
            self.data = alldata.loc [alldata['Type'] == 'SNP'].merge (SGlist[['chrom','position']], how = 'inner')
            del (SGlist)
            del (alldata)
            self.logger.info ("Input file filtered through SG list.")
        else:
            self.data = alldata
            self.logger.info ("Entire input file pass for analysis.")

        self.data['cov'] = self.data['ref_count'] + self.data['alt_count']
        self.data['vaf'] = self.data['alt_count']/self.data['cov']
        
        self.chromosomes = {}
        for chrom, data in self.data.groupby (by = 'chrom'):
            if chrom not in SEX_CHROMS:
                self.chromosomes[chrom] = Chromosome.Chromosome (chrom, data.copy(), 
                                                                 self.config, self.logger,
                                                                 self.genome_medians)
            
    #run COV test, run HE test, creates chromosomes with proper parameters
    def segment_genome (self):
        
        self.logger.info ('Starting COV test genome.')
        self.COV = Testing.Testing ('COV', self.chromosomes, self.logger)
        self.COV.run_test(no_processes = self.no_processes)
        self.COV.analyze (parameters = self.config['COV'])
                
        if self.COV.medians['m'] < float(self.config['COV']['min_cov']):
            self.logger.critical (f"Coverage is below threshold {self.COV.medians['m']} < {self.config['COV']['min_cov']}")
            exit (0)
        
        self.logger.info ("Genomewide coverage: " + f"\n" + str(self.COV.results))
        
        self.genome_medians['COV'] = self.COV.get_genome_medians()
        self.logger.info ('Starting HE test genome.')                        
        self.HE = Testing.Testing ('HE', 
                                   #{chrom : self.chromosomes[chrom] for chrom in self.COV.get_inliers()},
                                   self.chromosomes,
                                   self.logger)
        self.HE.run_test(no_processes = self.no_processes)
        #self.logger.info (str(self.HE.results))
        self.HE.analyze (parameters = self.config['HE'])
       
        self.logger.info ("Genomewide heterozygosity:" + "\n" + str(self.HE.results))
        self.genome_medians['HE'] = self.HE.get_genome_medians()        
              
        self.logger.info ('First round of H/E marking.')
        for chrom in self.chromosomes.keys():
            self.chromosomes[chrom].markE_onHE (self.HE.get_parameters(chrom), float(self.config['HE']['z_thr']))
        self.logger.info ('Testing first round of H/E marking.')    
        self.VAF = Testing.Testing ('VAF', self.chromosomes, self.logger)
        self.VAF.run_test (self.COV.medians['m'], no_processes = self.no_processes)
        #self.logger.info("Genomewide VAF:" + " \n" + str(self.VAF.results))
        self.VAF.analyze (parameters = self.config['VAF'])
        self.logger.info("Genomewide VAF:" + " \n" + str(self.VAF.results))
        self.genome_medians['VAF'] = self.VAF.get_genome_medians()
        
        #those that fail need to be marked on full model
        for chrom in self.chromosomes.keys():
            status = self.VAF.get_status (chrom)
            self.logger.info (f'Chromosome {chrom} status: {status}')
            if status:
                params = self.VAF.get_parameters (chrom)
                if params['chi2'] >= float(self.config['VAF']['chi2_high']):
                    self.logger.info (f'Chromosome {chrom} marked on full model.')
                    self.chromosomes[chrom].mark_on_full_model (self.COV.medians['m'])

        #segment_chromosomes
        self.logger.info ("Starting segmentation.")
        if self.no_processes > 1:            
            with mpl.Pool (processes = self.no_processes) as pool:
                segmented_chroms = pool.map (f, self.chromosomes.values())
                for sc in segmented_chroms:
                    self.chromosomes [sc.name] = sc
        else:
            for chrom in self.chromosomes.keys():
                self.chromosomes[chrom].find_runs()
                self.chromosomes[chrom].generate_segments ()
        self.logger.info ("Segmentation finished.")
        
        self.genome_medians['clonality'] = self.get_clonality_thr ()
        
        self.logger.info("Genome medians:" + " \n" + str(self.genome_medians))
        
        
    def get_clonality_thr (self, alpha = 0.05, percentiles = (10,80)):

        zs = []
        for chrom in self.chromosomes.keys():
            for seg in self.chromosomes[chrom].segments:
                if seg.symbol == Chromosome.E_SYMBOL:
                    zs.append (seg.parameters['k']*np.sqrt(seg.parameters['n']/Run.SNPS_IN_WINDOW))
        
        z = np.array(zs)
        pp = np.percentile (z, percentiles)
        res = sts.truncnorm.fit (z[(z >= pp[0])&(z <= pp[1])])
        self.logger.info ('Clonality threshold: min = {:.5f}, max = {:.5f}, m = {:.5f}, s = {:.5f}'.format (*res)) 
        
        return {'m' : res[2], 's' : res[3]}
        
            
    def report (self, report_type = 'bed'):
        return Report(report_type).genome_report(self.chromosomes)
    
def f (c):
    c.find_runs()
    c.generate_segments ()
    return c    

    
    
