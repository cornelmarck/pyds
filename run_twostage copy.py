from time import time
from pyds.buildutils import *
import pyomo.environ as pyo
from pyomo import dae
from numpy import random

from pyds.buildutils import StochParam, add_BigMConstraint
randomstream = random.default_rng()

def stage0_rule(m, stage=None):
    m.t = dae.ContinuousSet(bounds=(0,1)) #Normalised time
    m.ca_in = pyo.Var(m.t, initialize=0, bounds=(0, 0.001))

def stage1_rule(m, stage):
    t = stage[0].t
    m.I = pyo.Set(initialize=['A', 'B', 'C'])

    #Parameters
    m.T = pyo.Param(initialize=0, mutable=True)
    m.t_f = pyo.Param(initialize=0, mutable=True)
    m.E1 = StochParam(randomstream.normal, 2.5E3, 2.5E2)
    m.E2 = StochParam(randomstream.normal, 5E3, 5E2)
    m.k1 = StochParam(randomstream.normal, 6.409E-2, 6.409E-3)
    m.k2 = StochParam(randomstream.normal, 9.938E3, 9.938E2)

    #Variables
    m.c = pyo.Var(m.I, t, bounds=(0, 1E8), initialize=0)
    m.dcdt = dae.DerivativeVar(m.c, wrt=t)
    m.u = pyo.Var(t, initialize=0)
    m.dudt = dae.DerivativeVar(m.u, initialize=0, wrt=t)
     #Defined in natural time, not normalised
    m.u[0].fix(0)
    def r_ca(_m, t):
        return stage[0].ca_in[t]/_m.t_f == _m.dudt[t]
    m.c_ca = pyo.Constraint(t, rule=r_ca)

    #Constraints
    R = 1
    def r_mb1(_m, t): 
        return _m.dcdt['A',t]/_m.t_f == -2*_m.k1*pyo.exp(-_m.E1/(R*_m.T))*_m.c['A',t]**2 + stage[0].ca_in[t]
    def r_mb2(_m, t):
        return _m.dcdt['B',t]/_m.t_f == _m.k1*pyo.exp(-_m.E1/(R*_m.T))*_m.c['A',t]**2 \
            - _m.k2*pyo.exp(-_m.E2/(R*_m.T))*_m.c['B',t]
    def r_mb3(_m, t):
        return  _m.dcdt['C',t]/_m.t_f  == _m.k2*pyo.exp(-_m.E2/(R*_m.T))*_m.c['B',t]
    m.mb1 = pyo.Constraint(t, rule=r_mb1)
    m.mb2 = pyo.Constraint(t, rule=r_mb2)
    m.mb3 = pyo.Constraint(t, rule=r_mb3)

    m.g1 = pyo.Var(initialize=0)
    m.g2 = pyo.Var(initialize=0)
    def r_cqa1(_m):
        return _m.g1 == 1000*(0.8 -(_m.c['B',1] + 1E-2)/(_m.c['A', 1]+_m.c['B',1]+_m.c['C',1]+ 1E-2))
    def r_cqa2(_m):
        return _m.g2 == -(100*_m.c['B',1] - 20*_m.c['A',0] - 128*(_m.t_f + 30))/1E3
    m.c_g1 = pyo.Constraint(rule=r_cqa1)
    m.c_g2 = pyo.Constraint(rule=r_cqa2)

    m.bigM_constant = pyo.Param([1,2], initialize={1: 1E3, 2: 1E5})
    add_BigMConstraint(m, 'cqa1', m.bigM_constant[1], expr=m.g1<=0)  
    add_BigMConstraint(m, 'cqa2', m.bigM_constant[2], expr=m.g2<=0)
    
    def r_obj(model):
        return 1-model._indicator_var
    m.obj = pyo.Objective(rule=r_obj, sense=pyo.maximize)
    
    m.c['A',0].fix(2E3) #mol/L
    m.c['B',0].fix(0)
    m.c['C',0].fix(0)

t0 = time()
ef = create_EF([stage0_rule, stage1_rule], [1000])
print(time()-t0)


