import argparse
import configparser
import pickle as pkl

from doCNA import WGS

_description = "Scan chromosomes in search for non-HE segments. Assigns copy numbers if can."
__version__ = '0.6.8'

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

    args = parser.parse_args()
    ini = configparser.ConfigParser ()
    ini.read (args.config)
    sample = WGS.WGS (args.input_file,  sample_name = args.sample_name, parameters = ini,
                      assembly = args.assembly, no_processes = args.no_processes)
    
    sample.analyze ()
    
    print ('All done')
    
if __name__ == '__main__':
    main()