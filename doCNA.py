import logging
import argparse
import configparser as cfp
import multiprocessing as mpl

import numpy as np
import pandas as pd
import scipy.stats as sts
import scipy.optimize as opt
from collections import namedtuple
from math import floor, ceil



_description = "Scan chromosomes in search for non-HE segments. Assigns copy numbers if can."
__version__ = '0.6.8'


PAR = {'hg19': ((60001,	2699520), (154931044, 155260560)),
       'hg38': ((10001,	2781479), (155701383, 156030895))}

#here it should be easier to add loggers
class WGS:
    def __init__ (self, wgs_file_name,  sample_name, assembly = 'hg19', ini_file = 'doCNA.ini', no_processes = 1):
        #how to best check for file existance
        self.sample_name = sample_name
        self.assembly = assembly
        self.no_processes = no_processes
        self.wgs_file = open (wgs_file_name, 'r')
        self.config = cfp.ConfigParser().read(ini_file)
        self.SG_file = open (assembly + self.config['Input']['SuperGoog_core_name'], 'r')
        self.CB_file = open (assembly + self.config['Input']['CytoBand_core_name'], 'r')
        #loggers
    
    def analyze (self):
        self.genome = Genome (self.sample_name, self.no_processes)
        input_columns = [self.config['Input']['InputColumns']['chrom'],
                        self.config['Input']['InputColumns']['position'],
                        self.config['Input']['InputColumns']['ref_count'],
                        self.config['Input']['InputColumns']['alt_count'],
                        self.config['Input']['InputColumns']['Type']]
        self.genome.retrive_counts_create_chromosomes (data_file = self.wgs_file,
                                                       CB_file = self.CB_file, columns = input_columns)
        self.genome.segment_genome ()
        self.genome.segment_chromosomes ()

    #this probably needs to be split into more functions
    def report (self):
        #return self.genome.bed 
        pass


    def __del__ (self):
        close (self.wgs_file)
        close (self.SG_file)
        close (self.CB_file)

class Genome:
    def __init__(self, sample_name, no_processes = 1):
        self.sample_name = sample_name
        self.no_processes = no_processes
    
    #reads-in and filters through SuperGood
    def retrive_counts_create_chromosomes (self, data_file, columns, CB_file = None, filter_input = True):
        alldata = pd.read_csv (data_file, sep = '\t', usecols = columns)[columns]
        alldata.columns = ['chrom','position','ref_count', 'alt_count', 'Type']           

        if filter_input:
            SGlist = pd.read_csv (CB_file, compression = 'gzip',
                                header = None, sep = '\t', names = ['chrom','position','Ref', 'Alt', 'Status'],
                                dtype = {'chrom' : str, 'position' : np.int32, 'Ref' : str, 'Alt' : str, 'Status' : str})
            self.data = alldata.loc [alldata['Type'] == 'SNP'].merge (SGlist[['chrom','position']], how = 'inner').copy()
            del (SGlist)
            del (alldata)
        else:
            self.data = alldata

        self.data['cov'] = self.data['ref_count'] + self.data['alt_count']
        self.data['vaf'] = self.data['alt_count']/self.data['cov']

    #run COV test, run HE test, creates chromosomes with proper parameters
    def segment_genome (self, config):
        
        self.COV = Testing ('COV', self.chromosomes, config['COV'])
        self.COV.test(self.no_processes)
        self.COV.analyze (self.config['COV'])

        self.HE = Testing ('HE', self.chromosomes[normal_chromosomes])
        self.HE.test(self.no_processes)
        self.HE.analyze (self.config['HE'])
        
    def segment_chromosomes (self):
        #set parameters for each chromosome
        #that triggers self marking 
        
        self.VAF = Testing ('VAF', self.chromosomes, config['VAF'])
        self.VAF.test (self.no_processes)
        self.VAF.analyze (self.config['VAF'])
        #and now chromosomes are ready to segment itself
                                
    def generate_report (self, report_type = 'bed'):
        pass
        
        
    
class Testing:
    def __init__  (self, test_name, chromosomes, config):
        assert test_name not in ['COV', 'HE', 'VAF'] "Unknown test!"
        
        self.parameters = config[test_name]
        i = np.where(['COV', 'HE', 'VAF'] == test_name)[0][0]
        self.test = [COV_test, HE_test, VAF_test][i]
        self.test_columns = [['cov'], ['cov', 'vaf'], ['vaf']][i]
        self.data = chromosomes
        
    def test (self, no_processes = 1):
        
        if self.no_processes > 1:
            with mpl.Pool (processes = processes) as pool:
                indata = []
                for chrom in self.data.keys():
                    indata.append ((chrom, self.data[chrom].data[self.test_coulmns].values))
                results = pool.map(self.test, indata)
        else:
            results = []
            for chrom in self.data.keys():
                results.append((chrom, self.test (self.data[chrom].data[self.test_coulmns].values))
        self.results = pd.DataFrame.from_records (results)
    
    def analyze (self, parameters):
        
        pass

    def get_parameters (self, chromosome):
        pass
    
    def get_status (self, chromosome):
        pass
    
    def get_genome_medians (self):
        pass
    
    
    
class Chromosome:
    pass

    def get_values (self, columns = [], symbol = ''):
        pass

class Segment:
    pass

class Result:
    pass

class Distribution:
    pass

def main ():
    parser = argparse.ArgumentParser (description = _description)
    parser.add_argument ('-i', '--input_sample', required = False, 
                         type = str, default = '',
                         help = 'Input sample name. Default: from file name.')
    parser.add_argument ('-p', '--path', required = False,
                         default = '', help = 'Path to the file, default: current dir.')
    parser.add_argument ('-n', '--no_processes', required = False, default = 1, type = int,
                         help = 'Number of processes. Default: 1')
    parser.add_argument ('-f', '--input_file', required = True,
                         type = str, default = '',
                         help = 'Input file name.')
    #type not needed, as columns are spec in ini
    #parser.add_argument ('-t', '--type', required = False, choices = ["G", "D"],
    #                     type = str, default = 'G',
                         help = 'Germline or Diagnosis set of columns. For high_20.out')
    parser.add_argument ('-a', '--assembly', help = 'Assembly', default = 'hg19', required = False)

    args = parser.parse_args()

if __name__ == '__main__':
    main()
    
        
        