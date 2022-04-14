import pandas as pd
import numpy as np
import multiprocessing as mpl

from doCNA import Testing
from doCNA import Chromosome

SEX_CHROMS = ['chrX', 'chrY']

class Genome:
    def __init__(self, sample_name, logger, config, CB_file = None, no_processes = 1):
        self.sample_name = sample_name
        self.no_processes = no_processes
        self.config = config
        self.CB_file = CB_file
        self.logger = logger.getChild (f'{self.__class__.__name__}')
        self.genome_medians = {}
        self.logger.debug ("Object genome created.")
    #reads-in and filters through SuperGood
    def retrive_counts_create_chromosomes (self, data_file, columns, SG_file = None):
        alldata = pd.read_csv (data_file, sep = '\t', usecols = columns)[columns]
        alldata.columns = ['chrom','position','ref_count', 'alt_count', 'Type']           

        if SG_file is not None:
            SGlist = pd.read_csv (SG_file, compression = 'gzip',
                                header = None, sep = '\t', names = ['chrom','position','Ref', 'Alt', 'Status'],
                                dtype = {'chrom' : str, 'position' : np.int32, 'Ref' : str, 'Alt' : str, 'Status' : str})
            data = alldata.loc [alldata['Type'] == 'SNP'].merge (SGlist[['chrom','position']], how = 'inner')
            del (SGlist)
            del (alldata)
            self.logger.info ("Input file filtered through SG list.")
        else:
            data = alldata
            self.logger.info ("Entire input file pass for analysis.")

        self.data = data
        self.data['cov'] = self.data['ref_count'] + self.data['alt_count']
        self.data['vaf'] = self.data['alt_count']/self.data['cov']
        
        self.chromosomes = {}
        for chrom, data in data.groupby (by = 'chrom'):
            if chrom not in SEX_CHROMS:
                #c_start = self.cyto.loc[(self.cyto['#chrom'] == chrom)&((self.cyto['gieStain'] == 'acen')|(self.cyto['gieStain'] == 'gvar')), 'chromStart'].min()
                #c_end = self.cyto.loc[(self.cyto['#chrom'] == chrom)&((self.cyto['gieStain'] == 'acen')|(self.cyto['gieStain'] == 'gvar')), 'chromEnd'].max()                
                #centromere = (c_start, c_end)
                self.chromosomes[chrom] = Chromosome.Chromosome (chrom, data.copy(), 
                                                                 self.config, self.logger,
                                                                 self.genome_medians)
            
    #run COV test, run HE test, creates chromosomes with proper parameters
    def segment_genome (self):
        print ('Starting COV')
        self.COV = Testing.Testing ('COV', self.chromosomes, self.logger)
        self.COV.run_test(self.no_processes)
        self.COV.analyze (parameters = self.config['COV'])
        print ('COV test done')
        
        if self.COV.medians['m'] < float(self.config['COV']['min_cov']):
            self.logger.error ('Coverage is below threshold {} < {}.'.format (self.COV.medians['m'],
                                                                              self.config['COV']['min_cov']))
            exit (0)
        
        self.genome_medians['COV'] = self.COV.get_genome_medians()
        #print (self.genome_medians)
                
        self.HE = Testing.Testing ('HE', 
                                   {chrom : self.chromosomes[chrom] for chrom in self.COV.get_inliers()},
                                   self.logger)
        self.HE.run_test(self.no_processes)
        self.HE.analyze (parameters = self.config['HE'])
        self.genome_medians['HE'] = self.HE.get_genome_medians()        
        
        print ('HE test done')
        
        for chrom in self.chromosomes.keys():
            self.chromosomes[chrom].markE_onHE (self.HE.get_parameters(chrom), float(self.config['HE']['z_thr']))
            
        print ('E/N marking done')
        self.VAF = Testing.Testing ('VAF', self.chromosomes, self.logger)
        self.VAF.run_test (self.no_processes, self.COV.medians['m'])
        self.VAF.analyze (parameters = self.config['VAF'])
        self.genome_medians['VAF'] = self.VAF.get_genome_medians()
        
        self.logger.info ('VAF test done')
        #print (self.genome_medians)
        
        #those that fail need to be marked on full model
        for chrom in self.chromosomes.keys():
            status = self.VAF.get_status (chrom)
            if status == 'outlier':
                params = self.VAF.get_parameters (chrom)
                if params['chi2'] >= float(self.config['VAF']['chi2_high']):
                    self.chromosomes[chrom].mark_on_full_model ()

        #segment_chromosomes
        def f (c):
            c.find_segments()
            c.generate_segments ()
            return c
        
        if self.no_processes > 1:            
            with mpl.Pool (processes = self.no_processes) as pool:
                segmented_chroms = pool.map (f, self.chromosomes.values())
                for sc in segmented_chroms:
                    self.chromosomes [sc.name] = sc
        else:
            for chrom in self.chromosomes.keys():
                self.chromosomes[chrom].find_segments()
                self.chromosomes[chrom].generate_segments ()
            
    def generate_report (self, report_type = 'bed'):
        pass
    

    
    