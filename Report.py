import numpy as np
import warnings as warn

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
        elif self._report_type == 'run':
            data = '\n'.join([s.report(report_type='short') for s in runs])
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
            ai_score = -np.log10 (np.exp (-a*segment.parameters['ai']))
            
            if segment.parameters['model'] != 'cnB':
                k_score = ai_score
            else:
                m = segment.genome_medians['clonality_cnB']['m']
                s = segment.genome_medians['clonality_cnB']['s']
                k_score = np.abs(segment.parameters['k'] - m)/s
                
            report = '\t'.join([str(p) for p in [segment.parameters['m'], 
                                                 2*segment.parameters['m']/segment.genome_medians['COV']['m'],
                                                 segment.parameters['model'], score, 
                                                 segment.parameters['k'], k_score, segment.cytobands, 
                                                 segment.centromere_fraction, segment.parameters['d'], 
                                                 segment.parameters['ai'], ai_score]])
        else:
            report = ''
        return '\t'.join([namestr, report])

    def run_report(self, name, symbol, solutions):
        """ Generates a report for Run objects """
        report_types = ['short', 'full', 'solution']
        if self._report_type not in report_types:
            warn.warn ('Unknown report type. Use "short" instead.')
            self._report_type = 'short'

        if self._report_type == 'short':
            report = ';'.join([name, symbol,
                               f'#solutions: {len(solutions)}',
                               f'#segments: {len(solutions[0].positions)}'])
        elif self._report_type == 'full':
            report = ';'.join([name, symbol,
                               f'#solutions: {len(solutions)}',
                               f'#segments: {len(solutions[0].positions)}'])
        elif self._report_type == 'solution':
            reports = [name]
            for solution in solutions:
                sol_str = '    ' + '; '.join([str(solution.chi2), str(solution.chi2_noO), str(solution.positions), solution.merged_segments])
                reports.append(sol_str)
            report = '\n'.join(reports)
        return report
