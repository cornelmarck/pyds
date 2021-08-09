from pyds.utils import *
from pyds.run import *
import pyomo.environ as pyo
from pyomo import dae
from numpy import random
import numpy as np

randomstream = random.default_rng()

def stage0_rule(m, stage=None):
    m.t = dae.ContinuousSet(bounds=(0,1)) #Normalised time
    m.T = pyo.Param(initialize=280, mutable=True)
    m.t_f = pyo.Param(initialize=280, mutable=True)

def stage1_rule(m, stage):
    t = stage[0].t
    t_f = stage[0].t_f
    T = stage[0].T
    m.I = pyo.Set(initialize=['A', 'B', 'C'])
    #Parameters
    
    m.E1 = pyo.Param(mutable=True, initialize=2.5E3)
    m.E2 = pyo.Param(mutable=True, initialize=5E3) 
    m.k1 = pyo.Param(mutable=True, initialize=6.04E-2)
    m.k2 = pyo.Param(mutable=True, initialize=9.93E3) 

    #Variables
    m.c = pyo.Var(m.I, t, bounds=(0, 1E8), initialize=0)
    m.dcdt = dae.DerivativeVar(m.c, wrt=t)

    #Constraints
    R = 1
    def r_mb1(_m, t): 
        return _m.dcdt['A',t]/t_f == -2*_m.k1*pyo.exp(-_m.E1/(R*T))*_m.c['A',t]**2
    def r_mb2(_m, t):
        return _m.dcdt['B',t]/t_f == _m.k1*pyo.exp(-_m.E1/(R*T))*_m.c['A',t]**2 \
            -1*_m.k2*pyo.exp(-_m.E2/(R*T))*_m.c['B',t]
    def r_mb3(_m, t):
        return  _m.dcdt['C',t]/t_f  == _m.k2*pyo.exp(-_m.E2/(R*T))*_m.c['B',t]
    m.mb1 = pyo.Constraint(t, rule=r_mb1)
    m.mb2 = pyo.Constraint(t, rule=r_mb2)
    m.mb3 = pyo.Constraint(t, rule=r_mb3)

    m.g1 = pyo.Var(initialize=0)
    m.g2 = pyo.Var(initialize=0)
    def r_cqa1(_m):
        return _m.g1 == (0.8 -(_m.c['B',1] + 1E-2)/(_m.c['A', 1]+_m.c['B',1]+_m.c['C',1]+ 1E-2))
    def r_cqa2(_m):
        return _m.g2 == -1*(100*_m.c['B',1] - 20*_m.c['A',0] - 128*(t_f + 30))
    m.c_g1 = pyo.Constraint(rule=r_cqa1)
    m.c_g2 = pyo.Constraint(rule=r_cqa2)

    m.bigM_constant = pyo.Param([1,2], initialize={1: 1E2, 2: 1E5})
    add_BigMConstraint(m, 'cqa1', m.bigM_constant[1], expr=m.g1<=0)  
    add_BigMConstraint(m, 'cqa2', m.bigM_constant[2], expr=m.g2<=0)
    m.c['A',0].fix(2E3) #mol/L
    m.c['B',0].fix(0)
    m.c['C',0].fix(0)

def apply_collocation(m):
    discretizer = pyo.TransformationFactory('dae.collocation')
    discretizer.apply_to(m, nfe=8, ncp=3)

m = TwoStageManager([stage0_rule, stage1_rule], 
    {0: ['t_f', 'T'], 1: ['E1', 'E2', 'k1', 'k2']},
    model_transformation=apply_collocation)

n_samples = 1
ran = random.default_rng()
E1 = ran.normal(2.5E3, 2.5E1, size=n_samples)
E2 = ran.normal(5E3, 5E1, size=n_samples)
k1 = ran.normal(6.409E-2, 6.409E-4, size=n_samples)
k2 = ran.normal(9.938E3, 9.938E1, size=n_samples)

d = np.array([[280, 290]])
p = np.array([[E1[i], E2[i], k1[i], k2[i]] for i in range(n_samples)])
m.g(d, p)

print()




