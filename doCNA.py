import logging
import argparse
import configparser as cfp
import multiprocessing as mpl

import warnings

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

COV_INITIAL_SHAPE = 0.14
COV_SHAPE_RANGE = (-2,1)

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
        self.wgs_file.close ()
        self.SG_file.close ()
        self.CB_fileclose ()

class Genome:
    def __init__(self, sample_name, config, no_processes = 1):
        self.sample_name = sample_name
        self.no_processes = no_processes
        self.config = config
        
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

        self.HE = Testing ('HE', self.chromosomes[self.COV.normal_chromosomes])
        self.HE.test(self.no_processes)
        self.HE.analyze (self.config['HE'])
        
    def segment_chromosomes (self):
        #set parameters for each chromosome
        #that triggers self marking 
        
        self.VAF = Testing ('VAF', self.chromosomes, self.config['VAF'])
        self.VAF.test (self.no_processes)
        self.VAF.analyze (self.config['VAF'])
        #and now chromosomes are ready to segment itself
                                
    def generate_report (self, report_type = 'bed'):
        pass
        
class Testing:
    def __init__  (self, test_name, chromosomes, config):
        assert test_name not in ['COV', 'HE', 'VAF'], "Unknown test!"
        
        self.parameters = config[test_name]
        i = np.where(['COV', 'HE', 'VAF'] == test_name)[0][0]
        self.test = [COV_test, HE_test, VAF_test][i]
        self.drop_symbol = ['', '', 'N'][i]
        self.test_columns = [['cov'], ['cov', 'vaf'], ['vaf']][i]
        self.data = chromosomes
        
    def test (self, no_processes = 1):
        
        if self.no_processes > 1:
            with mpl.Pool (processes = no_processes) as pool:
                indata = []
                for chrom in self.data.keys():
                    row_filter = self.data[chrom].data['symbol'] != self.drop_symbol
                    indata.append ((chrom, self.data[chrom].data.loc[row_filter ,self.test_coulmns]))
                results = pool.map(self.test, indata)
        else:
            results = []
            for chrom in self.data.keys():
                results.append((chrom, self.test (self.data[chrom].data[self.test_coulmns].values)))
        self.results = pd.DataFrame.from_records (results)
    
    def analyze (self, parameters):
        
        pass

    def get_parameters (self, chromosome):
        pass
    
    def get_status (self, chromosome):
        pass
    
    def get_genome_medians (self):
        pass
    

def COV_test (covs, initial_shape = COV_INITIAL_SHAPE, 
                  shape_range = COV_SHAPE_RANGE):
    percentiles = np.arange (0.01, 0.99, 0.01)
       
    y = np.percentile (covs, percentiles*100)
    cov_range = (floor(y[0]*0.5) , ceil (y[-1]*1.05))
    initial_m = np.median(covs)
           
    popt, pcov = opt.curve_fit (lambda_cdf, percentiles, y, p0 = [initial_m, initial_shape],
                                bounds = [(cov_range[0],shape_range[0]),
                                          (cov_range[1],shape_range[1])])
    m, l  = popt
    dm, dl  = np.sqrt(np.diag(pcov)) 
    
    return m,dm,l,dl
    
def Q (p,l):
    if l == 0:
        return np.log(p/(1-p))
    else:
        return (p**l-(1-p)**l)/l    

def lambda_cdf (p, m, lam):
    return Q (p, lam)*np.sqrt(m)+m




def HE_test (data, cov_perc_bounds = (0.05, 99.0), vaf_bounds = (0.4,0.6),
             fcov_bounds = (0.01, 0.8), fN_bounds = (0.2,1.0), a_bounds = (0.1,0.9), 
             b_bounds = (1,10), lerr_bounds = (2,10)):
    
    def chi2 (params, counts, N):
        vaf, fcov, fN, a, b, l = params
        fe = 10**(-l)
        cs = np.arange (0, cov_max +1)
        cov = cov_min + fcov*(cov_max-cov_min)
        ns = 2*fN*N*cn2_cov_pdf (cs, cov, b)
        chi2 = 0
        
        for c, cnt, in counts:
            i = np.arange(0,c+1)
            nhe = ns[c]*cn2_vaf_pdf (i/c,vaf,c)
            nho = ns[c]*HO_vaf_pdf (i, c, fe ,b)
        
            na = a*nhe + (1-a)*nho
            chi2 += sum((cnt - na)**2/np.sqrt(na*na+1))/c 
        return chi2/len(counts)

    cov_min, cov_max = np.percentile (data['cov'].values, q = cov_perc_bounds)
   
    counts = []
    for c in np.arange (int(cov_min), int(cov_max+1)):
        d = data.loc[data['cov'] == c]
        h = np.histogram(d['alt_count'].values, bins = np.arange(0,c+2)-0.5)
        counts.append ((c, h[0]))
        del (h)
    
    N = len(data)
    fcov = (data['cov'].median() - cov_min)/(cov_max - cov_min)
    
    res = opt.minimize (chi2, x0 = (0.5, fcov, 0.5,0.8, 1.3, 6), args = (counts, N),
                    bounds = (vaf_bounds, fcov_bounds, fN_bounds, a_bounds, b_bounds, lerr_bounds))
    
    vaf, fcov, fN, a, b, l = res.x
    
    cov = cov_min + fcov*(cov_max-cov_min)
    return res.fun, vaf, cov, 2*fN*N, a, b, l, res.success

def cn2_vaf_pdf (x,v,c):
    p = sts.norm.pdf (x, v, np.sqrt((v*(1-v))/c))
    return p/sum(p)

def cn2_cov_pdf (n,c,b = 1):
    return sts.norm.pdf (n, c, np.sqrt(b*c))

#in general, there needs to be the low end added
def HO_vaf_pdf (i, n, fe = 10**-6,b=1):
    return sts.binom.pmf(i, n, 1-b*fe)

def VAF_test (data, vaf_bounds = (0.45,0.55), n_thr = 100):
    def chi2 (v, counts):
        chi2 = 0
        cc = 0
        for c, h, n in counts:
            if n >= n_thr:
                cc += c + 1 
                na = n*sts.binom.pmf (np.arange(c+1),c,v)
                chi2 += sum((h-na)**2/np.sqrt(na*na+1))
        #temporary solution
        if cc == 0: 
            cc = 1        
        return chi2/cc
        
    counts = []
    for c,d in data.groupby (by = 'cov'):
        h = np.histogram(d['alt_count'].values, bins = np.arange(0,c+2)-0.5)
        n = len(d)
        counts.append ((c, h[0], n))
        del (h)
        
    res = opt.minimize (chi2, x0 = (0.5) , args = (counts), method = 'L-BFGS-B',
                        bounds = (vaf_bounds,))        
    
    return res.fun, res.x[0]
    
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
    #                     help = 'Germline or Diagnosis set of columns. For high_20.out')
    parser.add_argument ('-a', '--assembly', help = 'Assembly', default = 'hg19', required = False)

def COV_test (covs, initial_shape = 0.14, 
                  shape_range = (-2,1)):
    percentiles = np.arange (0.01, 0.99, 0.01)
       
    y = np.percentile (covs, percentiles*100)
    cov_range = (floor(y[0]*0.5) , ceil (y[-1]*1.05))
    initial_m = np.median(covs)
           
    popt, pcov = opt.curve_fit (lambda_cdf, percentiles, y, p0 = [initial_m, initial_shape],
                                bounds = [(cov_range[0],shape_range[0]),
                                          (cov_range[1],shape_range[1])])
                    
    m, l  = popt
    dm, dl  = np.sqrt(np.diag(pcov)) 
        
    return m,dm, l,dl
    
def Q (p,l):
    if l == 0:
        return np.log(p/(1-p))
    else:
        return (p**l-(1-p)**l)/l    

def lambda_cdf (p, m, lam):
    return Q (p, lam)*np.sqrt(m)+m

    args = parser.parse_args()

def HE_test (data, cov_perc_bounds = (0.05, 99.0), vaf_bounds = (0.4,0.6),
             fcov_bounds = (0.01, 0.8), fN_bounds = (0.2,1.0), a_bounds = (0.1,0.9), 
             b_bounds = (1,10), lerr_bounds = (2,10)):
    
    def chi2 (params, counts, N):
        vaf, fcov, fN, a, b, l = params
        fe = 10**(-l)
        cs = np.arange (0, cov_max +1)
        cov = cov_min + fcov*(cov_max-cov_min)
        ns = 2*fN*N*cn2_cov_pdf (cs, cov, b)
        chi2 = 0
        
        for c, cnt, in counts:
            i = np.arange(0,c+1)
            nhe = ns[c]*cn2_vaf_pdf (i/c,vaf,c)
            nho = ns[c]*HO_vaf_pdf (i, c, fe ,b)
        
            na = a*nhe + (1-a)*nho
            chi2 += sum((cnt - na)**2/np.sqrt(na*na+1))/c 
        return chi2/len(counts)

    cov_min, cov_max = np.percentile (data['cov'].values, q = cov_perc_bounds)
   
    counts = []
    for c in np.arange (int(cov_min), int(cov_max+1)):
        d = data.loc[data['cov'] == c]
        h = np.histogram(d['alt_count'].values, bins = np.arange(0,c+2)-0.5)
        counts.append ((c, h[0]))
        del (h)
    
    N = len(data)
    fcov = (data['cov'].median() - cov_min)/(cov_max - cov_min)
    
    res = opt.minimize (chi2, x0 = (0.5, fcov, 0.5,0.8, 1.3, 6), args = (counts, N),
                    bounds = (vaf_bounds, fcov_bounds, fN_bounds, a_bounds, b_bounds, lerr_bounds))
    
    vaf, fcov, fN, a, b, l = res.x
    
    cov = cov_min + fcov*(cov_max-cov_min)
    return res.fun, vaf, cov, 2*fN*N, a, b, l, res.success

def VAF_test (data, vaf_bounds = (0.45,0.55), n_thr = 100):
    def chi2 (v, counts):
        chi2 = 0
        cc = 0
        for c, h, n in counts:
            if n >= n_thr:
                cc += c + 1 
                na = n*sts.binom.pmf (np.arange(c+1),c,v)
                chi2 += sum((h-na)**2/np.sqrt(na*na+1))
        #temporary solution
        if cc == 0: 
            cc = 1        
        return chi2/cc
        
    counts = []
    for c,d in data.groupby (by = 'cov'):
        h = np.histogram(d['alt_count'].values, bins = np.arange(0,c+2)-0.5)
        n = len(d)
        counts.append ((c, h[0], n))
        del (h)
        
    res = opt.minimize (chi2, x0 = (0.5) , args = (counts), method = 'L-BFGS-B',
                        bounds = (vaf_bounds,))        
    
    return res.fun, res.x[0]

def main():
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
    parser.add_argument ('-t', '--type', required = False, choices = ["G", "D"],
                         type = str, default = 'G',
                         help = 'Germline or Diagnosis set of columns. For high_20.out')
    parser.add_argument ('-a', '--assembly', help = 'Assembly', default = 'hg19', required = False)

    args = parser.parse_args()



if __name__ == '__main__':
    main()
    
        
        