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
            gmm = segment.genome_medians['clonality']['m']
            gms = segment.genome_medians['clonality']['s']
            n = segment.parameters['n']/Run.SNPS_IN_WINDOW
            score = np.abs(segment.parameters['k'] - gmm)*np.sqrt(n/2)/gms
            report = '\t'.join([str(p) for p in [segment.parameters['m'], segment.parameters['model'],
                                                 segment.parameters['k'], score, segment.cytoband, 
                                                 segment.fraction, segment.parameters['d']]])
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
