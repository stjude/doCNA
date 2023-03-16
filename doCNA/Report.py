from cmath import isnan
import numpy as np
import warnings as warn
import scipy.stats as sts

from doCNA import Run


class Report:
    """ Class that holds reports used by other objects in the program """
    def __init__(self, report_t):
        self._report_type = report_t

    def genome_report(self, genome):
        """ Generates a report for Genome objects """
        if self._report_type == 'bed':
            keys = list(genome.chromosomes.keys())
            keys.sort(key = lambda x: int(x[3:]))
            report = '\n'.join([genome.chromosomes[key].report(report_type=self._report_type) for key in keys])
        elif self._report_type == 'params':
            shift_i = genome.genome_medians['clonality_imbalanced']['up']
            report_list = ['m0\t'+ str(genome.genome_medians['m0']),
                           'fb\t'+ str(genome.genome_medians ['fb']),
                           'a_model\t' + str(genome.genome_medians['model_d']['a']),
                           'Imbalanced:']
            for key in genome.genome_medians['clonality_imbalanced'].keys():
                report_list.append (key + '_i\t' + str(genome.genome_medians['clonality_imbalanced'][key])) 
            report_list.append ('Balanced:')
            for key in genome.genome_medians['clonality_balanced'].keys():
                value = genome.genome_medians['clonality_balanced'][key]
                if hasattr(value, '__iter__'):
                    value_str = ' '.join([str(v) for v in value])
                else:
                    value_str = str(value)
                report_list.append (key + '_b\t' + value_str)
            report = '\n'.join(report_list)
        else:
            report = ""     
        return report

    def chromosome_report(self, segments, runs):
        """ Generates a report for Chromosome objects """
        if self._report_type == 'bed':
            data = '\n'.join([s.report(report_type='bed') for s in segments])
        elif self._report_type == 'solution':
            data = '\n'.join([s.report(report_type='solution') for s in runs])
        return data

    
    def segment_report (self, segment):

        """ Generates a report for Segment objects """
        if self._report_type == 'bed':    
            if np.isnan (segment.parameters['k']):
                if segment.parameters['fraction_1'] > 0.95:
                    k = 1
                else:
                    k = np.nan
            else:
                k = segment.parameters['k']

            report = '\t'.join([str(p) for p in [segment.chrom, segment.start, segment.end,
                                                 segment.parameters['ai'], segment.parameters['m'],
                                                 2*segment.parameters['m']/segment.genome_medians['m0'],
                                                 segment.parameters['model'], segment.parameters['d'], 
                                                 segment.parameters['model_score'],
                                                 k, segment.parameters['clonality_score'],
                                                 segment.parameters['k_d'], 
                                                 segment.cytobands,
                                                 segment.centromere_fraction, 
                                                 segment.parameters['call'], segment.parameters['call_FDR']]])
        else:
            report = ''
        return report

    def run_report(self, run):
        """ Generates a report for Run objects """
        fields = ['chi2', 'chi2_noO', 'positions', 'p_norm', 'merged_segments']
        lines = ['Run: ' + run.name]
        for solution in run.solutions:
            lines.append (str("\t")+'Solution')
            soldic = solution._asdict()
            for f in fields:
                lines.append (str("\t\t") + f + ': '+ str(soldic[f]))
        return '\n'.join(lines)
        
