##Dictonary with models is define now in one place

from collections import namedtuple
import numpy as np
import scipy.optimize as opt

Preset = namedtuple ('Preset', ['k','m', 'ai'])

#Presets of models with 2 +/- 1 copies

model_presets = {'(AB)(2+n)' : Preset(k = lambda m,dv,m0: np.abs(m/m0 - 1),
                                      m = lambda k,m0: (1+k)*m0,
                                      ai = lambda k,m0:  np.zeros_like(k)),

                   
                   '(AB)(2-n)' : Preset(k = lambda m,dv,m0: np.abs(m/m0 - 1),
                                        m = lambda k,m0: (1-k)*m0,
                                        ai = lambda k,m0:  np.zeros_like(k)),
                   
                   'A'   : Preset(k = lambda m,dv,m0: 2*dv/(0.5+dv),
                                  m = lambda k,m0: (2-k)*m0/2,
                                  ai = lambda k,m0: k/(2*(2-k))),
                 
                   'AA'  : Preset(k = lambda m,dv,m0: 2*dv,
                                  m = lambda k,m0: (np.zeros_like(k)+1)*m0,
                                  ai = lambda k,m0: k/2),
                 
                   'AAB' : Preset(k = lambda m,dv,m0: 2*dv/(0.5-dv),
                                  m = lambda k,m0: (2+k)*m0/2,
                                  ai = lambda k,m0: k/(2*(2+k))),
    
                   'AAAB' : Preset (k = lambda m,dv,m0 : 2*dv/(1-2*dv),
                                    m = lambda k,m0 : (1+k)*m0,
                                    ai = lambda k,m0 : k/(2+2*k)),
                   
                   'AAA' : Preset (k = lambda m,dv,m0 : 4*dv/(3-2*dv),
                                   m = lambda k,m0 : (2+k)*m0/2,
                                   ai = lambda k,m0 : 3*k/(4+2*k)),
                    
                   'AAAA' : Preset (k = lambda m,dv,m0 : dv/(1-dv),
                                    m = lambda k,m0 : (1+k)*m0,
                                    ai = lambda k,m0 : k/(1+k)),
                   
                   'AAB+AAAB' : Preset (k = lambda m,dv,m0 : (6*dv-1)/(1-2*dv),
                                        m = lambda k,m0 : (3+k)*m0/2,
                                        ai = lambda k,m0 : (1+k)/(6+2*k)),
                   
                   'AA+AAB' : Preset (k = lambda m,dv,m0 : 2*(1-2*dv)/(2*dv+1),
                                      m = lambda k,m0 : (2+k)*m0/2,
                                      ai = lambda k,m0 : (2-k)/(4+2*k)),
                   
                   'AAB+AABB' : Preset (k = lambda m,dv,m0 : (1-6*dv)/(2*dv+1),
                                        m = lambda k,m0 : (3+k)*m0/2,
                                        ai = lambda k,m0 : (1-k)/(6+2*k))}

def pick_model (ai, s_ai, cn, s_cn, models):
    assert not np.isnan (ai), "ai is nan"
    assert not np.isnan (cn), "cn in nan"
    
    dsks = [calculate_distance_minim(ai, s_ai, cn, s_cn, model_presets[model]) for model in models]
    ds = np.array ([dk['d'] for dk in dsks])
    ks = np.array ([dk['k'] for dk in dsks])
    
    model_index = np.where(ds == ds.min())[0][0]
       
    return {'model' : models[model_index], 'd_model' : ds[model_index], 'k': ks[model_index]}
    
def calculate_distance_minim (ai, s_ai, cn, s_cn, model):
    
    res = opt.minimize_scalar (dist, bounds = (0,1), args = ((ai, s_ai, cn, s_cn, model)), 
                               method = 'bounded')
    
    return {'d' : res.fun,
            'k' : res.x}

def dist (k, ai, s_ai, cn, s_cn, model):
    da = (ai - model.ai(k,2))/s_ai
    dc = (cn - model.m(k,2))/s_cn
    return np.sqrt(da**2+dc**2)