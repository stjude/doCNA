import numpy as np
import scipy.stats as sts
import scipy.optimize as opt
from collections import namedtuple

from doCNA import Testing
from doCNA import Distribution
from doCNA import Segment
from doCNA import Run

N_SYMBOL = 'N'
E_SYMBOL = 'E'
U_SYMBOL = 'U'
DEFAULT_N_THRESHOLD = 10
DEFAULT_E_THRESHOLD = 3

Run_treshold =  namedtuple('Run_treshold', [N_SYMBOL, E_SYMBOL])


class Chromosome:
    """Class to contain data, and find runs."""
    def __init__ (self, name, data, config, logger, genome_medians):
        self.name = name
        self.data = data
        self.config = config
        self.logger = logger.getChild(f'{self.__class__.__name__}-{self.name}')
        self.genome_medians = genome_medians
        self.Eruns = []
        self.Uruns = []
        self.Nruns = []
        self.logger.debug (f"Object chromosome {name} created.")
        
    def markE_onHE(self, he_parameters, threshold = 13.8):
        zv = (self.data['vaf'] - he_parameters['vaf']) / np.sqrt (0.25/(he_parameters['cov']))
        zc = (self.data['cov'] - he_parameters['cov']) / np.sqrt (he_parameters['b']*he_parameters['cov'])
        z = zv**2+zc**2
 
        self.data['symbol'] = N_SYMBOL       
        indexes = self.data.loc[z < threshold, :].index.values.tolist()
        self.data.loc[indexes, 'symbol'] = E_SYMBOL
        self.logger.debug (f"""Chromosome {self.name} marked based on 
                           parameters v = {he_parameters['vaf']}, c = {he_parameters['cov']}.""")
        
    def mark_on_full_model (self):
        self.get_fragments (exclude_symbols = [], n = int(self.config['Segment']['No_SNPs']))
        self.get_vaf_shift_full ()
        self.Uruns = []
        
        indexes, merged_string = Run.merge_symbols (self.dv_dist.string, outliers_threshold = 4)

        for r in indexes:
            start = self.positions [r[0]][0]
            end = self.positions [r[1]][1]
            chi2, vaf = Testing.VAF_test (self.data.loc[(self.data['symbol'] != N_SYMBOL)&\
                                                        (self.data['position'] >= start)&\
                                                        (self.data['position'] <= end)])
            
            #test chi2
            outlier = (chi2 > float(self.config['Tests']['VAF']['chi2_high'])) | (chi2 == 0) |\
                      (vaf < self.genome_medians['VAF']['vaf_thr'][0])
                      
            self.logger.info (f'Marking on full model {self.name}:{start}-{end} chi2 = {chi2}')            
            if outlier:
                self.data.loc[(self.data['position'] >= start)&\
                                    (self.data['position'] <= end), 'symbol'] = U_SYMBOL
                self.Uruns.append ((start, end))
    
    #or whatever name
    def segment (self):
        
        self.find_Nruns ()
        #now we have nicely marked chromosome
        
        pass
    
    #def get_vaf_shift (self, zero_thr = 0.01, cov_mult = 1.03, p_thr = 0.5, z_thr = 1.5):
    #    tmpf = 1
    #    s0 = np.sqrt (0.25/self.genome_medians['COV']['m'])
    #            
    #    while tmpf > zero_thr:
    #        dvl = []
    #        v0l = []
    #        s = s0/np.sqrt(cov_mult)
    #        def make_two_gauss (v, dv, v0):
    #            return two_gauss (v, dv, v0, s, a = 0.5)
            
            #print (self.fragments[0])
            #print (len(self.fragments))
    #        for fragment in self.fragments:
    #            #print (fragment)
    #            vaf = fragment['vaf'].values
    #            popt, pcov = opt.curve_fit (make_two_gauss, np.sort(vaf), np.linspace (0,1, len(vaf)),
    #                                        p0 = [0.05, 0.5],
    #                                        bounds = ((0, 0.4),(0.5, 0.6)))
    #            dvl.append (popt[0])
    #           v0l.append (popt[1])
    #    
    #        tmpf = sum ([d < zero_thr for d in dvl])/len (dvl)
    #        cov_mult += 0.05
        
    #    self.dv = np.array (dvl)
    #    self.v0 = np.array (v0l)

    #    self.dv_dist = Distribution.Distribution (self.dv, p_thr = 0.5, thr_z = z_thr)

    #    self.logger.info ("Vaf shifts calculated. Shrink factor used: {:.2f}.".format (cov_mult))
    
    #Moved to Run.py    
    #def get_coverage (self, z_thr = 1.5):
    #    ml = []
    #    ll = []
            
    #    for fragment in self.fragments:
    #        result = Testing.COV_test (fragment)
    #        ml.append (result.m)
    #        ll.append (result.l)

    #    self.m = np.array (ml)
    #    self.m_dist = Distribution.Distribution (self.m, thr_z = z_thr, p_thr = 0.3)
    #    self.l = np.array (ll)
    #    self.l_dist = Distribution.Distribution (self.l, thr_z = z_thr, p_thr = 0.3)
                    
    def find_segments (self, z_thr = 2.5):
        #if there are any Uruns already processed
        self.find_Nruns ()
        self.solutions = []
        #process near diploid part of the chromosome
        self.get_fragments (exclude_symbols = [N_SYMBOL, U_SYMBOL])
        if len (self.fragments) >= 10:
            self.get_vaf_shift ()
            self.get_coverage ()
            if any ([self.dv_dist.fail_normal(), self.m_dist.fail_normal(), self.l_dist.fail_normal()]):
                for m,s,a in zip(*self.get_distributions ()):
                    self.segment_on_double (m,s,a, z_thr)

            z_m = np.array ([self.dv_dist.params['single']['m'][0], 
                             self.m_dist.params['single']['m'][0], 
                             self.l_dist.params['single']['m'][0]])
            z_s = np.array ([self.dv_dist.params['single']['s'][0],
                             self.m_dist.params['single']['s'][0],
                             self.l_dist.params['single']['s'][0]])
            self.segment_on_single (z_m, z_s, z_thr)
                    
            if len (self.solutions) > 1:
                #self.solutions.sort (key = lambda x: np.abs(1-x['chi2_noO']))
                self.solutions.sort (key = lambda x: x['chi2_noO'])

            #that is quite ugly
            UNruns = self.Nruns + self.Uruns
            UNruns.sort (key = lambda x: x[0])
            currentNU = 0
            for si, ei in self.solutions[0]['segment']:
                start = self.positions[si][0]
                end = self.positions[ei][0]
                intersect = [(start < s) & (end > e) for s,e in UNruns]
                for i in intersect:
                    self.Eruns.append ((start, UNruns[i][0]))
                    start = UNruns[i][1]
                if start < end:
                    self.Eruns.append ((start, end))
        else:
            self.logger.warning (f"Not enough of near diploid to segment on {self.name}")
            
    def generate_segments (self):
        #once there are N and U runs known as well as near diploid is segmented,
        #we can start generating segments
        allruns = self.Nruns + self.Uruns + self.Eruns
        self.segments = []
        allruns.sort(key = lambda x: x[0])
        
        for seg in allruns:
            start, end = seg
            data = self.data.loc[(self.data['position'] > start) & (self.data['position'] < end)]
            self.segments.append (Segment.Segment(data, self.config, self.logger, self.genome_medians))
      
    #Moved to Run.py     
    def get_distributions (self):
        res_m = []
        res_s = []
        res_a = []
        ordinal = np.cumsum ([self.dv_dist.fail_normal(), self.m_dist.fail_normal(), self.l_dist.fail_normal()])
        
        dv_i  = self.dv_dist.combinations ()[0]
        if ordinal[1] > 1:
            m_is = self.m_dist.combinations ()
        else:
            m_is = self.m_dist.combinations ()[:1]
        if ordinal[2] > 1:
            l_is = self.l_dist.combinations ()
        else:
            l_is = self.l_dist.combinations ()[:1]

        for m_i in m_is:
            for l_i in l_is:
                z_m = (np.array([self.dv_dist.parameters['m'][dv_i[0]], self.m_dist.parameters['m'][m_i[0]], self.l_dist.parameters['m'][l_i[0]]]),
                       np.array([self.dv_dist.parameters['m'][dv_i[1]], self.m_dist.parameters['m'][m_i[1]], self.l_dist.parameters['m'][l_i[1]]]))
                        
                z_s = (np.array([self.dv_dist.parameters['s'][dv_i[0]], self.m_dist.parameters['s'][m_i[0]], self.l_dist.parameters['s'][l_i[0]]]),
                       np.array([self.dv_dist.parameters['s'][dv_i[1]], self.m_dist.parameters['s'][m_i[1]], self.l_dist.parameters['s'][l_i[1]]]))
                
                a = (np.array([self.dv_dist.parameters['a'][dv_i[0]], self.m_dist.parameters['a'][m_i[0]], self.l_dist.parameters['a'][l_i[0]]]),
                     np.array([self.dv_dist.parameters['a'][dv_i[1]], self.m_dist.parameters['a'][m_i[1]], self.l_dist.parameters['a'][l_i[1]]]))
                
                res_m.append (z_m)
                res_s.append (z_s)
                res_a.append (a)
        return (res_m, res_s, res_a)

    #Moved to Run.py
    #def get_fragments (self, exclude_symbols = [], n = 1000):
    #    tmp = self.data.loc[[s not in exclude_symbols for s in self.data['symbol']], ]
    #    #print (tmp)
    #    N = int(np.floor(len(tmp)/n))
    #    #print (self.name, N)
    #    indexes = np.linspace (0,len (tmp), 2*N+2, dtype = int)

        #self.vafs = []
        #self.covs = []
    #    self.fragments = []
    #    self.positions = []

    #    for i in np.arange (2*N):
    #        tmpi = tmp.iloc[indexes[i]:indexes[i+2], ]
    #        self.fragments.append(tmpi)
            #self.vafs.append (np.sort(tmpi['vaf'].values))
            #self.covs.append (np.sort(tmpi['cov'].values))
    #        self.positions.append ((tmpi['position'].min(), tmpi['position'].max()))

    #Moved to Run.py
    #def get_vaf_shift_full (self, z_thr = 2.5):
    #    
    #    def vaf_cdf (v, dv, a, lerr, f, vaf, b):
    #        return vaf_cdf_c (v, dv, a, lerr, f, vaf, b, cov)
                
    #    v0 = self.genome_medians['HE']['vaf']
    #    b = self.genome_medians['HE']['b']
    #    cov = self.genome_medians['HE']['cov']
    #    dvs = []
    #    v0s = []
    #    for vaf in self.vafs:
    #        v, c = np.unique(vaf, return_counts = True)

    #        cnor = np.cumsum(c)/np.sum(c)
    #        ones0 = c[v >= (cov-1)/cov].sum()
    #        f0 = c[v < v0].sum()/(c.sum() - ones0) 

    #        dv0 = v0 - np.median (v[v < v0])

    #        try:
    #            popt, pcov = opt.curve_fit (vaf_cdf, v, cnor, p0 = [dv0, ones0/c.sum(), 2, f0, 0.5, b], 
    #                                        bounds = ((0,   0,   1, 0, 0.45, 1),
    #                                                  (0.5, 0.95, 5, 1, 0.55, 10)))
    #            dvs.append (popt[0])
    #            v0s.append (popt[-1])
    #        except:
    #            print (self.name, dv0, ones0/c.sum(), 2, f0, 0.5, b)
    #            #the lenghts must be same
    #            #adding 0 either makes outlier or not outlir but it is safe value
    #            dvs.append (0)
                
    #    self.dv = np.array(dvs)
    #    self.v0 = np.array(v0s)
    #    self.dv_dist = Distribution.Distribution (self.dv[~np.isnan(self.dv)],
    #                                              p_thr = 0.3, thr_z = z_thr)
    
    def find_Nruns (self):
        symbol_list = self.data.loc[(self.data['vaf'] < 1) & (self.data['symbol'] != U_SYMBOL), 'symbol'].tolist()
        self.Nruns_indexes, self.Nruns_threshold = analyze_string_N (symbol_list, N = N_SYMBOL, E = E_SYMBOL)
        #
        for run in self.Nruns_indexes:
            tmp = self.data.loc[self.data['vaf'] < 1,].iloc[run[0]:run[1]].position.agg((min,max))
            self.data.loc[(self.data.position >= tmp['min'])&(self.data.position <= tmp['max']), 'symbol'] = N_SYMBOL
            self.Nruns.append ((tmp['min'], tmp['max']))
    
    def segment_on_double (self, z_m, z_s, z_a, thr_z = 1.5, cent_intersect = 0.5, N_symbol = 'N'):
        #calculate z
        #the order of parameters is (dv, m, l)
        x = np.array ([self.dv, self.m, self.l])
        z0 = (((x-z_m[0][:,np.newaxis])/z_s[0][:,np.newaxis])**2).sum(axis = 0)
        z1 = (((x-z_m[1][:,np.newaxis])/z_s[1][:,np.newaxis])**2).sum(axis = 0)
        #generate string
        string = list ('O' * len(self.dv))
        #marking C population
        for i in np.where ((z0 < z1)&(z0 < thr_z**2))[0]:
            string[i] = 'C'
        #marking D population
        for i in np.where ((z1 < z0)&(z1 < thr_z**2))[0]:
            string[i] = 'D'
        
        symbol, count = Run.rle_encode (string)
        indexes, merged_string = Run.merge_symbols (string)
        chi2 = 0.0
        chi2_noO = 0.0
        NnoO = 0
        segments = []
        for r in indexes:
            sym = sts.mode(list(string[r[0]:r[1]+1]))[0][0]
            if sym == 'C':
                chi2 += z0[r[0]:r[1]+1].sum()
                chi2_noO += z0[r[0]:r[1]+1].sum()
                NnoO += r[1] - r[0] + 1
            elif sym == 'D':
                chi2 += z1[r[0]:r[1]+1].sum()
                chi2_noO += z1[r[0]:r[1]+1].sum()
                NnoO += r[1] - r[0] + 1
            elif sym == 'O':
                chi2 += min(z0[r[0]:r[1]+1].sum(), z1[r[0]:r[1]+1].sum())
        
        if all (z0 > thr_z) & all (z1 > thr_z): 
            chi2_noO = np.inf
            NnoO = 1
            df = 0
        else:
            df = 2*(np.array([self.dv_dist.fail_normal(), self.m_dist.fail_normal(), self.l_dist.fail_normal()], dtype = int)+1).sum()
                
        self.solutions.append({'chi2' : chi2/(3*len(self.dv)-df),
                               'chi2_noO' : chi2_noO/(3*NnoO-df),
                               'string' : ''.join(string),
                               'merged_string' : ''.join(merged_string),
                               'C(dv, m, l)' : z_m[0],
                               'D(dv, m, l)' : z_m[1],
                               'runs' : (symbol, count),
                               'segment' : indexes})
        
    def segment_on_single (self, z_m, z_s, thr_z = 1.5, cent_intersect = 0.5, N_symbol = 'N'):
        #generate string based on outliers
        
        x = np.array ([self.dv, self.m, self.l])
        z0 = (((x-z_m[:,np.newaxis])/z_s[:,np.newaxis])**2).sum(axis = 0)

        string = list ('O' * len (self.dv))
        for i in np.where (z0 < thr_z**2)[0]:
            string[i] = 'B'

        symbol, count = Run.rle_encode (string)
        indexes, merged_string = Run.merge_symbols (string)
        chi2 = 0.0
        chi2_noO = 0.0
        NnoO = 0
    
        for r in indexes:
            sym = sts.mode(list(string[r[0]:r[1]+1]))[0][0]
            if sym == 'B':
                chi2 += z0[r[0]:r[1]+1].sum()
                chi2_noO += z0[r[0]:r[1]+1].sum()
                NnoO += r[1] - r[0] + 1
            else:
                chi2 += z0[r[0]:r[1]+1].sum()
                
        self.solutions.append({'chi2' : chi2/(3*len(self.dv) - 6),
                               'chi2_noO' : chi2_noO/(3*NnoO - 6),
                               'string' : ''.join(string),
                               'merged_string' : ''.join(merged_string),
                               'B(dv, m, l)' : z_m,
                               'runs' : (symbol, count),
                               'segment' : indexes})
            
    def report (self, type = 'bed'):
        return '\n'.join([seg.report() for seg in self.segments])    

#No longer needed    
#def two_gauss (v, dv, v0, s, a = 0.5):
#    return a*sts.norm.cdf (v, v0 - dv, s)+(1-a)*sts.norm.cdf (v,v0 + dv,s)


#func to analyze N runs
def analyze_string_N (symbol_list, N = 'N', E = 'E'):
                        
    runs = []
    values, counts = Run.rle_encode (symbol_list)
    threshold = find_runs_thr (values, counts, N = N, E = E)
        
    runs = get_N_runs_indexes (values, counts, N = N, E = E, threshold = threshold)
    
    for s,e in runs:
        for i in range(s,e+1):
            symbol_list[i] = N
    
    values, counts = Run.rle_encode (symbol_list)
    threshold = find_runs_thr (values, counts, N = N, E = E)
    runs = get_N_runs_indexes (values, counts, threshold = threshold, N = N, E = E)
        
    return runs, threshold 
    
def find_runs_thr (values, counts, N = 'N', E = 'E'):
    assert len (values) == len (counts), 'Wrong input!! Go away!'
    hist = np.unique(counts[values == N], return_counts = True)
    x = np.log10 (hist[0])
    y = np.log10 (hist[1])
    
    try:
        ind = [0,1,2,3]
        popt, _ = opt.curve_fit (lin, x[ind], y[ind], p0 = [-1,1])
        xt = np.log10 (np.arange(1, hist[0].max()))
        N_thr = 10**(xt[lin(xt, *popt) > np.log10(1)].max())
    except:
        N_thr = DEFAULT_N_THRESHOLD
    
    hist = np.unique(counts[values != N], return_counts = True)
    x = hist[0]
    y = np.log10 (hist[1])
    
    try:
        popt, _ = opt.curve_fit (lin, x[:4], y[:4], p0 = [-1,1])
        if popt[1] < 0:
            xt = np.arange(1, hist[0].max())
            E_thr = xt[lin(xt, *popt) > 0].max()
        else:
            popt, _ = opt.curve_fit (lin, x[:3], y[:3], p0 = [-1,1])
            if popt[1] < 0:
                xt = np.arange(1, hist[0].max())
                E_thr = xt[lin(xt, *popt) > 0].max()
            else:
                E_thr = DEFAULT_E_THRESHOLD
    except:
        E_thr = DEFAULT_E_THRESHOLD
    
    return Run_treshold (N = N_thr, E = E_thr)

def lin (x,a,b):
    return a*x+b

def get_N_runs_indexes (values, counts, threshold, N = N_SYMBOL, E = E_SYMBOL):
    bed = []
    Nindexes = counts[(counts > threshold[0])& (values == N_SYMBOL)]
    
    for i in Nindexes:
        if any ([(i>=s)&(i <= e) for s,e in bed]):
            continue
        
        s = i
        e = i
        
        if e < len(counts)-2:
            inext = e + 1        
            while counts[inext] < threshold[E_SYMBOL]:
                if e > len(counts)-4:
                    e = len(counts)-1
                    break
                else:
                    e += 2
                    inext = e + 1
                
        if s > 0:
            iprev = s - 1
            while counts[iprev] < threshold[E_SYMBOL]:
                if s < 3:
                    s = 0
                    break
                else:
                    s -= 2
                    iprev = s - 1
        bed.append ((s,e))    
        
    return [[counts[:i].sum(),counts[:j+1].sum()] for i, j in bed] 
        
#Moved to Run.py
#def vaf_cdf_c (v, dv, a, lerr, f, vaf, b, cov):
    #cn2 = vaf_cn2 (v, vaf, cov)
#    cnai = vaf_cnai (v, dv, f, vaf, b, cov)
#    cnHO = vaf_HO (v, lerr)
    
#    return a*cnHO + (1 - a)*cnai 

#def vaf_cnai (v, dv, a, vaf,b, cov):
#    s = np.sqrt((vaf - dv)*(vaf + dv)/(b*cov))
#    return a*sts.norm.cdf (v, vaf - dv, s) + (1-a)*sts.norm.cdf (v, vaf + dv, s)

#def vaf_HO (v, lerr):
#    err = 10**lerr
#    return np.exp ((v-1)*err)

#Moved to Run.py
#def merge_symbols (in_string, outliers_threshold = 2):
#    string = list(in_string)
#    symbols, counts = rle_encode (string)#

#    argsort = np.argsort (counts)[-1:0:-1]
#    symbolindex = 0

#    while symbolindex < len(symbols)-1:  #(min(counts) <= outliers_threshold)&(counter < 10):
#        i = argsort[symbolindex]
#        symbol = symbols[i]
#        s = i
#        e = i
        
#        sinit = i
#        einit = i    
        
#        if e < len(counts)-1:
#            inext = e + 1        
#            while (symbols[inext] == symbol)|(counts[inext] <= outliers_threshold):
#                e = inext
#                if e == len(counts)-1:
#                    break
#                inext += 1
            
#        if s > 0:
#            iprev = s - 1
#            while (counts[iprev] <= outliers_threshold)|(symbols[iprev] == symbol):
#                s = iprev
#                if s == 0:
#                    break
#                iprev -= 1
        
#        if (s == sinit)&(e == einit):
#            symbolindex +=1
#        else:
#            for i in np.arange(counts[:s].sum(),counts[:e+1].sum()):
#                string[i] = symbol.copy()
#        
#            symbols, counts = rle_encode (string)
#            argsort = np.argsort (counts)[-1:0:-1]
#            symbolindex = 0
        
#    bed = []
#    for i in np.arange(len(counts)):
#        s = counts[:i].sum()
#        e = counts[:i+1].sum()-1
#        bed.append ((s,e))
#    return bed, string

