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
            shift_b = np.sqrt(genome.genome_medians['clonality_balanced']['A']**2+1)*genome.genome_medians['clonality_balanced']['up']
            shift_i = np.sqrt(genome.genome_medians['clonality_imbalanced']['A']**2+1)*genome.genome_medians['clonality_imbalanced']['up']
            report = '\n'.join(['m' + '\t' + str(genome.genome_medians['m']),
                                'a_model' + '\t' + str(genome.genome_medians['model_d']['a']),  
                                'a_b' + '\t' + str(genome.genome_medians['clonality_balanced']['A']),
                                'b_b' + '\t' + str(genome.genome_medians['clonality_balanced']['C']),
                                'bt_b' + '\t' + str(genome.genome_medians['clonality_balanced']['C']-shift_b),
                                'a_i' + '\t' + str(genome.genome_medians['clonality_imbalanced']['A']),
                                'b_i' + '\t' + str(genome.genome_medians['clonality_imbalanced']['C']),
                                'bt_i' + '\t' + str(genome.genome_medians['clonality_imbalanced']['C']-shift_b)])
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
                                                 2*segment.parameters['m']/segment.genome_medians['m'],
                                                 segment.parameters['model'], segment.parameters['d'], 
                                                 segment.parameters['model_score'],
                                                 k, segment.parameters['clonality_score'],
                                                 segment.parameters['k_d'], 
                                                 segment.cytobands,
                                                 segment.centromere_fraction, segment.parameters['call']]])
        else:
            report = ''
        return report



    def segment_report_old (self, segment):

        """ Generates a report for Segment objects """
        namestr = segment.name.replace(':', '\t').replace ('-', '\t')
        if self._report_type == 'bed':
            
            a = segment.genome_medians['ai']['a']
            try:
                ai_score = -np.log10 (np.exp (-a*segment.parameters['ai']/np.sqrt(segment.parameters['n'])))
            except RuntimeWarning:
                ai_score = np.inf                
            
            if segment.parameters['model'] == 'cnB':
                m = segment.genome_medians['clonality_cnB']['m']
                s = segment.genome_medians['clonality_cnB']['s']
                try:
                    k_score = -np.log10(sts.norm.sf(segment.parameters['k']/np.sqrt(segment.parameters['n']), m, s))
                except RuntimeWarning:
                    k_score = np.inf
                model_score = ai_score
                
            else:
                k_score = ai_score
                a = segment.genome_medians['model_d']['a']
                model_score = -np.log10(np.exp (-a*segment.parameters['d']))
                
            report = '\t'.join([str(p) for p in [segment.parameters['m'], 
                                                 2*segment.parameters['m']/segment.genome_medians['COV']['m'],
                                                 segment.parameters['model'], model_score, 
                                                 segment.parameters['k'], k_score, segment.cytobands, 
                                                 segment.centromere_fraction, segment.parameters['d'], 
                                                 segment.parameters['ai'], ai_score]])
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
        
