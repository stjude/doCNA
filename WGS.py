import logging
import pickle as pkl
from doCNA import Genome

class WGS:
    """Class to handle WGS read counts file and create the genome."""
    def __init__ (self, wgs_file_name,  sample_name, parameters, assembly = 'hg19',  
                  no_processes = 1, verbosity = 'INFO'):
        
        self.sample_name = sample_name
        self.assembly = assembly
        self.no_processes = no_processes
        self.wgs_file = open (wgs_file_name, 'r')
        self.config = parameters
        self.SG_file = open (assembly + self.config['Input']['SuperGood_core_name'], 'rb')
        self.CB_file = open (assembly + self.config['Input']['cytoband_core_name'], 'r')
        self.logger = self.create_logger (verbosity)
        self.logger.debug ("WGS object created.")
        
    def create_logger (self, verbosity):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(self.sample_name + '.log', mode = 'w')
        fh.setLevel(logging.DEBUG)
        #%(funcName)s:
        fh_formatter = logging.Formatter('%(asctime)s %(name)s: %(levelname)s: %(message)s')
        fh.setFormatter(fh_formatter)
        logger.addHandler(fh)
        
        sh = logging.StreamHandler ()
        sh.setLevel (verbosity)
        sh.setFormatter (fh_formatter)
        logger.addHandler (sh)
        
        return logger
            
    def analyze (self):
        self.logger.info ('Creating genome.')
        self.genome = Genome.Genome (self.sample_name, self.logger, self.config, self.CB_file, self.no_processes)
        input_columns = [self.config['InputColumns']['chrom'],
                         self.config['InputColumns']['position'],
                         self.config['InputColumns']['ref_count'],
                         self.config['InputColumns']['alt_count'],
                         self.config['InputColumns']['Type']]
        self.logger.debug ('Genome object created.')
        self.genome.retrive_counts_create_chromosomes (data_file = self.wgs_file, #SG_file = self.SG_file,
                                                       columns = input_columns)
        self.logger.debug ('Counts data retrived.')
        self.logger.info ('Segmenting genome.')
        self.genome.segment_genome ()
        self.logger.info ('Ready to report.')

    #this probably needs to be split into more functions
    def report (self, type = 'bed'):
        return self.genome.report(type)

    def pickle_genome (self):
        #test is closing file handlers in logger allows pickling, so sort of destructor
        pass
    
    def __del__ (self):
        self.wgs_file.close ()
        if self.SG_file is not None:
            self.SG_file.close ()
        self.CB_file.close ()

