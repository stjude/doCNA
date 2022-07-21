import numpy as np
import warnings as warn
import scipy.stats as sts

from doCNA import Run

class Report:
    """ Class that holds reports used by other objects in the program """
    def __init__(self, report_t):
        self._report_type = report_t

    def genome_report(self, chromosomes):
        """ Generates a report for Genome objects """
        keys = list(chromosomes.keys())
        keys.sort(key = lambda x: int(x[3:]))
        return '\n'.join([chromosomes[key].report(report_type=self._report_type) for key in keys])

    def chromosome_report(self, segments, runs):
        """ Generates a report for Chromosome objects """
        if self._report_type == 'bed':
            data = '\n'.join([s.report(report_type='bed') for s in segments])
        elif self._report_type == 'solution':
            data = '\n'.join([s.report(report_type='solution') for s in runs])
        return data

    def segment_report(self, segment):
        """ Generates a report for Segment objects """
        namestr = segment.name.replace(':', '\t').replace ('-', '\t')
        if self._report_type == 'bed':
            
            a = segment.genome_medians['model_d']['a']
            score = -np.log10(np.exp (-a*segment.parameters['d']))

            a = segment.genome_medians['ai']['a']
            ai_score = -np.log10 (np.exp (-a*segment.parameters['ai']*np.sqrt(segment.parameters['n'])))
            
            if segment.parameters['model'] == 'cnB':
                m = segment.genome_medians['clonality_cnB']['m']
                s = segment.genome_medians['clonality_cnB']['s']
                k_score = -np.log10(sts.norm.sf(segment.parameters['k'], m, s))
                score = ai_score
            else:
                k_score = ai_score
                
                
            report = '\t'.join([str(p) for p in [segment.parameters['m'], 
                                                 2*segment.parameters['m']/segment.genome_medians['COV']['m'],
                                                 segment.parameters['model'], score, 
                                                 segment.parameters['k'], k_score, segment.cytobands, 
                                                 segment.centromere_fraction, segment.parameters['d'], 
                                                 segment.parameters['ai'], ai_score]],
                                                 segment.parameters['ddv'])
        else:
            report = ''
        return '\t'.join([namestr, report])

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
        
