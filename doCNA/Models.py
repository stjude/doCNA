##Dictonary with models is define now in one place

from collections import namedtuple
import numpy as np
import scipy.optimize as opt

#Fields of the tuple are functions describing model of karyotype
# A,B,C -  
# D - 
# k - func to calculate clonality
# m - func to describe coverage
# ai - func to describe allelic imbalance
Preset = namedtuple ('Preset', ['A', 'B', 'C', 'D', 'k','m', 'ai'])

#Presets of models with 2 +/- 1 copies
model_presets_2 = {'(AB)n' : Preset(A = lambda m,dv,m0: 0,
                                    B = lambda m,dv,m0: 1/2,
                                    C = lambda m,dv,m0: dv,
                                    D = lambda m,dv,m0: 0,
                                    k = lambda m,dv,m0: np.abs(m/m0 - 1) if (m/m0 > -0.1/2) & (m/m0 < 4.1/2) else np.nan,
                                    m = lambda k,m0: (1+k)*m0,
                                    ai = lambda k,m0:  np.zeros_like(k)),#np.repeat(0, len(k)) if hasattr(k, "shape") else 0.0),
                   
                   'A'   : Preset(A = lambda m,dv,m0: -m0/2,
                                  B = lambda m,dv,m0: -1,
                                  C = lambda m,dv,m0: m0,
                                  D = lambda m,dv,m0: m0*(2*dv/(0.5+dv))/2+m,
                                  k = lambda m,dv,m0: 2*dv/(0.5+dv) if (m/m0 > 0.9/2) & (m/m0 < 2.1/2) else np.nan,
                                  m = lambda k,m0: (2-k)*m0/2,
                                  ai = lambda k,m0: k/(2*(2-k))),
                 
                   'AA'  : Preset(A = lambda m,dv,m0: 0,
                                  B = lambda m,dv,m0: -1,
                                  C = lambda m,dv,m0: m0,
                                  D = lambda m,dv,m0: m,
                                  k = lambda m,dv,m0: 2*dv if (m/m0 > 1.9/2) & (m/m0 < 2.1/2) else np.nan,
                                  m = lambda k,m0: (np.zeros_like(k)+1)*m0, #np.repeat(m0, len(k)) if hasattr(k, "shape") else m0,
                                  ai = lambda k,m0: k/2),
                 
                   'AAB' : Preset(A = lambda m,dv,m0: m0/2,
                                  B = lambda m,dv,m0: -1,
                                  C = lambda m,dv,m0: m0,
                                  D = lambda m,dv,m0: -m0*(2*dv/(0.5-dv))/2+m,
                                  k = lambda m,dv,m0: 2*dv/(0.5-dv) if (m/m0 > 1.9/2) & (m/m0 < 3.1/2) else np.nan,
                                  m = lambda k,m0: (2+k)*m0/2,
                                  ai = lambda k,m0: k/(2*(2+k)))} 
    
#models of more copies, not a strict classification 
model_presets_4 = {'AAAB' : Preset (A = lambda m,dv,m0 : m0/2,
                                    B = lambda m,dv,m0 : -1,
                                    C = lambda m,dv,m0 : m0,
                                    D = lambda m,dv,m0 : m - 2*m0*dv/(1-2*dv),
                                    k = lambda m,dv,m0 : 2*dv/(1-2*dv) if (m/m0 > 1.9/2) & (m/m0 < 4.1/2) else np.nan,
                                    m = lambda k,m0 : (1+k)*m0,
                                    ai = lambda k,m0 : k/(2+2*k)),
                   
                    'AAA' : Preset (A = lambda m,dv,m0 : m0/2,
                                    B = lambda m,dv,m0 : -1,
                                    C = lambda m,dv,m0 : m0,
                                    D = lambda m,dv,m0 : m - 2*dv*m0/(3-2*dv),
                                    k = lambda m,dv,m0 : 4*dv/(3-2*dv) if (m/m0 > 1.9/2) & (m/m0 < 3.1/2) else np.nan,
                                    m = lambda k,m0 : (2+k)*m0/2,
                                    ai = lambda k,m0 : 3*k/(4+2*k)),
                    
                    'AAAA' : Preset (A = lambda m,dv,m0 : m0/2,
                                     B = lambda m,dv,m0 : -1,
                                     C = lambda m,dv,m0 : m0,
                                     D = lambda m,dv,m0 : m - dv*m0/(1-dv),
                                     k = lambda m,dv,m0 : dv/(1-dv) if (m/m0 > 1.9/2) & (m/m0 < 4.1/2) else np.nan,
                                     m = lambda k,m0 : (1+k)*m0,
                                     ai = lambda k,m0 : k/(1+k))}


model_presets_extra = {'AAB+AAAB' : Preset (A = lambda m,dv,m0 : m0/2,
                                            B = lambda m,dv,m0 : -1,
                                            C = lambda m,dv,m0 : 3*m0/2,
                                            D = lambda m,dv,m0 : m- m0/((6*dv-1)/(2-4*dv)),
                                            k = lambda m,dv,m0 : (6*dv-1)/(1-2*dv) if (m/m0 > 2.9/2) & (m/m0 < 4.1/2) else np.nan,
                                            m = lambda k,m0 : (3+k)*m0/2,
                                            ai = lambda k,m0 : (1+k)/(6+2*k)),
                   
                        'AA+AAB' : Preset (A = lambda m,dv,m0 : m0/2,
                                           B = lambda m,dv,m0 : -1,
                                           C = lambda m,dv,m0 : m0,
                                           D = lambda m,dv,m0 : m- m0*(1-2*dv)/(2*dv+1),
                                           k = lambda m,dv,m0 : 2*(1-2*dv)/(2*dv+1) if (m/m0 > 1.9/2) & (m/m0 < 3.1/2) else np.nan,
                                           m = lambda k,m0 : (2+k)*m0/2,
                                           ai = lambda k,m0 : (2-k)/(4+2*k)),
                   
                        'AAB+AABB' : Preset (A = lambda m,dv,m0 : m0/2,
                                             B = lambda m,dv,m0 : -1,
                                             C = lambda m,dv,m0 : 3*m0/2,
                                             D = lambda m,dv,m0 : m- m0*(1-6*dv)/(4*dv+1),
                                             k = lambda m,dv,m0 : (1-6*dv)/(2*dv+1) if (m/m0 > 2.9/2) & (m/m0 < 4.1/2) else np.nan,
                                             m = lambda k,m0 : (3+k)*m0/2,
                                             ai = lambda k,m0 : (1-k)/(6+2*k))}

def calculate_distance (preset, m, ai, m0):
    
    try:
        k = (preset.k(m,ai,m0))
    except ZeroDivisionError:
        k = np.nan

    if np.isnan(k) | (k > 1.1) | (k < -0.1):
        d = np.nan
    elif (k >= 1) | (k <= 0):
        ks = np.linspace (0,1,1000)
        ms = preset.m(k,m0)/m0
        d = np.min(np.sqrt((ks-k)**2+(ms-m/m0)**2))
    else:
        d = np.abs (preset.C (m/m0,ai,1) - preset.D (m/m0,ai,1))/np.sqrt (preset.A(m/m0,ai,1)**2 + preset.B(m/m0,ai,1)**2)
    
    return d

def calculate_distance_new (ai, ci, model):
    
    res = opt.minimize_scalar (dist, bounds = (0,1), args = ((ai, ci, model)), 
                               method = 'bounded', tol = 1e-2)
    
    return {'d' : res.fun,
            'k' : res.x}

def dist (k, ai, ci, model):
    da = ai - model.ai(k,2)
    dc = ci - model.m(k,2)
    return np.sqrt(da**2+dc**2)

def dist_diploid (ai, ci):
    return np.sqrt (ai**2+(ci-2)**2)