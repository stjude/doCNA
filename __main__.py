import argparse
import configparser
import pandas as pd
from doCNA.Run import Solution

from doCNA import WGS

_description = "Scan chromosomes in search for non-HE segments. Assigns copy numbers if can."
__version__ = '0.8.0'

def main():
    parser = argparse.ArgumentParser (description = _description)
    parser.add_argument ('-s', '--sample_name', required = False, 
                         type = str, default = '',
                         help = 'Input sample name. Default: from file name.')
    
    parser.add_argument ('-n', '--no_processes', required = False, default = 1, type = int,
                         help = 'Number of processes. Default: 1')
    
    parser.add_argument ('-i', '--input_file', required = True,
                         type = str, default = '',
                         help = 'Input file name.')
    
    parser.add_argument ('-c', '--config', required = False, default = 'config.ini',
                         help = 'INI file with parameters')
    
    parser.add_argument ('-a', '--assembly', help = 'Assembly', default = 'hg19',
                         required = False)
    
    parser.add_argument ('-l', '--level', default = 'INFO', 
                         choices = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'NOTSET'],
                         help = 'Level of verbosity for std.err logger.')
    parser.add_argument ('-r', '--report_solutions', help = 'Generate report with all solutions.',
                         action = 'store_true')
    
    parser.add_argument ('-m0', '--coverage_diploid', required = False, type = float,
                         help = 'Coverage of diploid.', default = 0)
    
    parser.add_argument ('-v', '--version', help = 'Print version', action = 'version',
                         version = 'doCNA v. {version}'.format(version = __version__))
  
    
    args = parser.parse_args()
    ini = configparser.ConfigParser ()
    ini.read (args.config)
    sample = WGS.WGS (args.input_file,  sample_name = args.sample_name, parameters = ini,
                      assembly = args.assembly, no_processes = args.no_processes, 
                      verbosity = args.level)
    
    sample.analyze (m0 = args.coverage_diploid)
    
    with open (args.sample_name + '.bed', 'w') as bed:
        bed.writelines (sample.report(report_type = 'bed'))
    
    with open (args.sample_name + '.par', 'w') as params:
        params.writelines (sample.report (report_type = 'params'))

    if args.report_solutions:
        with open (args.sample_name + '.solutions', 'w') as full:
            full.writelines (sample.report(report_type = 'solution'))
    
    keys = sample.genome.chromosomes.keys()
    data = pd.concat ([sample.genome.chromosomes[k].data for k in keys])
    
    data.to_csv (args.sample_name + '.dat', index = None, sep = '\t')
    
    print ('All done')
    
if __name__ == '__main__':
    main()
