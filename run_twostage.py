from pyds.utils import *
from pyds.solver import *
import pyomo.environ as pyo
from pyomo import dae
from numpy import random
import numpy as np

randomstream = random.default_rng()

def stage0_rule(m, stage=None):
    m.t = dae.ContinuousSet(bounds=(0,1)) #Normalised time
    m.T = pyo.Param(initialize=0, mutable=True)
    m.t_f = pyo.Param(initialize=0, mutable=True)

def stage1_rule(m, stage):
    t = stage[0].t
    t_f = stage[0].t_f
    T = stage[0].T
    m.I = pyo.Set(initialize=['A', 'B', 'C'])
    #Parameters
    
    m.E1 = pyo.Param(mutable=True, initialize=0)
    m.E2 = pyo.Param(mutable=True, initialize=0) 
    m.k1 = pyo.Param(mutable=True, initialize=0)
    m.k2 = pyo.Param(mutable=True, initialize=0) 

    #Variables
    m.c = pyo.Var(m.I, t, bounds=(0, 1E8), initialize=0)
    m.dcdt = dae.DerivativeVar(m.c, wrt=t)

    #Constraints
    R = 1
    def r_mb1(_m, t): 
        return _m.dcdt['A',t]/t_f == -2*_m.k1*pyo.exp(-_m.E1/(R*T))*_m.c['A',t]**2
    def r_mb2(_m, t):
        return _m.dcdt['B',t]/t_f == _m.k1*pyo.exp(-_m.E1/(R*T))*_m.c['A',t]**2 \
            - _m.k2*pyo.exp(-_m.E2/(R*T))*_m.c['B',t]
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
        return _m.g2 == -(100*_m.c['B',1] - 20*_m.c['A',0] - 128*(t_f + 30))
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

m = TwoStageManager(stage_rules=[stage0_rule, stage1_rule], 
    input={'d': [(0, 't_f'), (0, 'T')], 'p': [(1, 'E1'), (1, 'E2'), (1, 'k1'), (1, 'k2')]},
    solver_name='ipopt',
    model_transformation=apply_collocation,
    solver_options={})

n = 10
d = np.random.rand(n, 2)
p = np.ones((n, 4))
m._build_model(4)

m.write_output_to_disk
print()




