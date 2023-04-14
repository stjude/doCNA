from cmath import isnan
import numpy as np
import warnings as warn
import scipy.stats as sts

from doCNA import Run
from doCNA import Consts

class Report:
    """ Class that holds reports used by other objects in the program """
    def __init__(self, report_t):
        self._report_type = report_t

    def genome_report (self, genome):
        """ Generates a report for Genome objects """
        if self._report_type == 'bed':
            keys = list(genome.chromosomes.keys())
            keys.sort (key = Consts.CHROM_ORDER.index)
            report = '\n'.join([genome.chromosomes[key].report(report_type=self._report_type) for key in keys])
        elif self._report_type == 'params':
            report_list = ['m0\t'+ str(genome.genome_medians['m']),
                           'm_ai\t'+str(self.scorer.ai_param['m']),
                           's_ai\t'+str(self.scorer.ai_param['s']),
                           'm_cn\t'+str(self.scorer.cn_param['m']),
                           's_ai\t'+str(self.scorer.ai_param['s']),
                           'm_d\t' +str(self.scorer.dipl_dist['m']),
                           's_d\t' +str(self.scorer.dipl_dist['s']),
                           'models\t'+str(self.models)]
                        
            report = '\n'.join(report_list)
        else:
            report = ""
        return report
    
    ##Remove in release
    def genome_report_old(self, genome):
        """ Generates a report for Genome objects """
        if self._report_type == 'bed':
            keys = list(genome.chromosomes.keys())
            #keys.sort(key = lambda x: int(x[3:]))
            keys.sort (key = Consts.CHROM_ORDER.index)
            report = '\n'.join([genome.chromosomes[key].report(report_type=self._report_type) for key in keys])
        elif self._report_type == 'params':
            #shift_b = genome.genome_medians['clonality_balanced']['up']
            shift_i = genome.genome_medians['clonality_imbalanced']['up']
            #report = '\n'.join(['m' + '\t' + str(genome.genome_medians['m']),
            #                    'a_model' + '\t' + str(genome.genome_medians['model_d']['a']),  
            #                    'a_b' + '\t' + str(genome.genome_medians['clonality_balanced']['A']),
            #                    'b_b' + '\t' + str(genome.genome_medians['clonality_balanced']['C']),
            #                    'bt_b' + '\t' + str(genome.genome_medians['clonality_balanced']['C']-shift_b),
            #                    'm_a' + '\t' + str(genome.genome_medians['clonality_balanced']['m']),
            #                    's_a' + '\t' + str(genome.genome_medians['clonality_balanced']['s']),
            #                    'a_i' + '\t' + str(genome.genome_medians['clonality_imbalanced']['A']),
            #                    'b_i' + '\t' + str(genome.genome_medians['clonality_imbalanced']['C']),
            #                    'bt_i' + '\t' + str(genome.genome_medians['clonality_imbalanced']['C']-shift_i),
            #                    'm_a' + '\t' + str(genome.genome_medians['clonality_imbalanced']['m']),
            #                    's_a' + '\t' + str(genome.genome_medians['clonality_imbalanced']['s'])])

            report_list = ['m0\t'+ str(genome.genome_medians['m']),
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
    ##end of remove in release

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
            #if np.isnan (segment.parameters['k']):
            #    if segment.parameters['fraction_1'] > 0.95:
            #        k = 1
            #    else:
            #        k = np.nan
            #else:
            k = segment.parameters['k']

            report = '\t'.join([str(p) for p in [segment.chrom, segment.start, segment.end,
                                                 segment.parameters['ai'], segment.parameters['p_ai'],
                                                 segment.parameters['m'],
                                                 2*segment.parameters['m']/segment.genome_medians['m0'],
                                                 segment.parameters['model'], segment.parameters['d'], 
                                                 segment.parameters['model_score'],
                                                 k, segment.parameters['clonality_score'],
                                                 segment.parameters['k_d'], 
                                                 segment.segmentation_symbol, 
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
        
