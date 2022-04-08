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

class Genome:
    def __init__(self, sample_name, assembly = 'hg19', no_processes = 1):
        self.logger = logging.getLogger('__main__').getChild(__name__)
        pass
    
    #reads-in and filters through SuperGood
    def retrive_counts_create_chromosomes (self, file_name, columns, filter_input = True):
        
        #pass self to each chromosome, so there is Genome reference 
        Chromosome (....., self.logger.getChild(chrom),)
        pass
    
    
    #run COV test, run HE test, creates chromosomes with proper parameters
    def segment_genome (self):
        
        self.COV = Testing ('COV', self.chromosomes, config['COV'], self.no_processes)
        self.COV.test()
        
        self.HE = Testing ('HE', self.chromosomes[normal_chromosomes], config['HE'], self.no_processes)
        self.HE.test()
        
    def segment_chromosomes (self):
        #set parameters for each chromosome
        #that triggers self marking 
        
        self.VAF = Testing ('VAF', self.chromosomes, config['VAF'], self.no_processes)
        #and now chromosomes are ready to segment itself
        
        if self.no_processes > 1:
            with mpl.Pool (processes = processes) as pool:
                
        else:
            for chrom in self.chromosomes.keys():
                self.chromosomes[chrom].segment()
                
    def generate_report (self, report_type = 'bed'):
        pass
        
        
    
class Testing:
    def __init__  (self, test_name, data, config, no_processes = 1):
        
        
        pass
    
    def test (self):
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
    
    logging.basicConfig(format='%(asctime)s -%(levelname)s:%(message)s', level=logging.DEBUG)
    
    try:
        
        pass
    except:
        print (last message)
        
        