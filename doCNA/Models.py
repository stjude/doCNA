##Dictonary with models is define now in one place

from collections import namedtuple
import numpy as np

#Fields of the tuple are functions describing model of karyotype
# A,B,C -  
# D - 
# k - func to calculate clonality
# m - func to describe coverage
# ai - func to describe allelic imbalance
Preset = namedtuple ('Preset', ['A', 'B', 'C', 'D', 'k','m', 'ai'])

#Presets of models with 2 +/- 1 copies
model_presets_2 = {'A'   : Preset(A = lambda m,dv,m0: -m0/2,
                                     B = lambda m,dv,m0: -1,
                                     C = lambda m,dv,m0: m0,
                                     D = lambda m,dv,m0: m0*(2*dv/(0.5+dv))/2+m,
                                     k = lambda m,dv,m0: 2*dv/(0.5+dv),
                                     m = lambda k,m0: (2-k)*m0/2,
                                     ai = lambda k,m0: k/(2*(2-k))),
                 
                   'AA'  : Preset(A = lambda m,dv,m0: 0,
                                     B = lambda m,dv,m0: -1,
                                     C = lambda m,dv,m0: m0,
                                     D = lambda m,dv,m0: m,
                                     k = lambda m,dv,m0: 2*dv,
                                     m = lambda k,m0: np.repeat(m0, len(k)),
                                     ai = lambda k,m0: k/2),
                 
                   'AAB' : Preset(A = lambda m,dv,m0: m0/2,
                                     B = lambda m,dv,m0: -1,
                                     C = lambda m,dv,m0: m0,
                                     D = lambda m,dv,m0: -m0*(2*dv/(0.5-dv))/2+m,
                                     k = lambda m,dv,m0: 2*dv/(0.5-dv),
                                     m = lambda k,m0: (2+k)*m0/2,
                                     ai = lambda k,m0: k/(2*(2+k))),

                   '(AB)n' : Preset(A = lambda m,dv,m0: 0,
                                     B = lambda m,dv,m0: 1/2,
                                     C = lambda m,dv,m0: dv,
                                     D = lambda m,dv,m0: 0,
                                     k = lambda m,dv,m0: (m/m0 - 1),
                                     m = lambda k,m0: (1+k)*m0,
                                     ai = lambda k,m0: np.repeat(0, len(k)))} 
    
#models of more copies, not a strict classification 
model_presets_4 = {'AAAB' : Preset (A = lambda m,dv,m0 : m0/2,
                                       B = lambda m,dv,m0 : -1,
                                       C = lambda m,dv,m0 : m0,
                                       D = lambda m,dv,m0 : m - 2*m0*dv/(1-2*dv),
                                       k = lambda m,dv,m0 : 2*dv/(1-2*dv),
                                       m = lambda k,m0 : (1+k)*m0,
                                       ai = lambda k,m0 : k/(2+2*k)),
                   
                    'AAA' : Preset (A = lambda m,dv,m0 : m0/2,
                                       B = lambda m,dv,m0 : -1,
                                       C = lambda m,dv,m0 : m0,
                                       D = lambda m,dv,m0 : m - 2*dv*m0/(3-2*dv),
                                       k = lambda m,dv,m0 : 4*dv/(3-2*dv),
                                       m = lambda k,m0 : (2+k)*m0/2,
                                       ai = lambda k,m0 : 3*k/(4+2*k)),
                    
                    'AAAA' : Preset (A = lambda m,dv,m0 : m0/2,
                                        B = lambda m,dv,m0 : -1,
                                        C = lambda m,dv,m0 : m0,
                                        D = lambda m,dv,m0 : m - dv*m0/(1-dv),
                                        k = lambda m,dv,m0 : dv/(1-dv),
                                        m = lambda k,m0 : (1+k)*m0,
                                        ai = lambda k,m0 : k/(1+k))}

def calculate_distance (preset, m, ai, m0):
    """Function to calculate the distance of the segment to the model"""
    
    k = preset.k(m,ai,m0)
    if np.isnan(k):
        d = np.inf
    elif ((k > 1) | (k < 0)):
        ks = np.linspace (0,1,1000)
        ms = preset.m(ks,m0)/m0
        d = np.min(np.sqrt((ks-k)**2+(ms-m/m0)**2))
    else:
        d = np.abs (preset.C (m/m0,k,1) - preset.D (m/m0,k,1))/np.sqrt (preset.A(m/m0,k,1)**2 + preset.B(m/m0,k,1)**2)
    
    return d
