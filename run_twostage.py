import pyomo.environ as pyo
from pyomo import dae
from numpy import random

from pyds.model_builder import *

def create_model():
    m = pyo.ConcreteModel()

    #Sets
    m.t = dae.ContinuousSet(bounds=(0,1)) #Normalised time
    m.I = pyo.Set(initialize=['A', 'B', 'C'])

    #Parameters
    m.d = pyo.Param(['T', 't_f'], initialize={'t_f': 280, 'T': 285}, mutable=True)
    m.p = pyo.Param(['E1', 'E2', 'k1', 'k2'], initialize={'E1': 2.5E3, 'E2':5E3, 
                                                        'k1': 6.409E-2, 'k2': 9.938E3}, \
                                                        mutable=True)

    #Variables
    m.c = pyo.Var(m.I, m.t, bounds=(0, 1E8), initialize=0)
    m.dcdt = dae.DerivativeVar(m.c, wrt=m.t)

    m.u = pyo.Var(m.t, initialize=0)
    m.dudt = dae.DerivativeVar(m.u, initialize=0, wrt=m.t)
    m.ca_in = pyo.Var(m.t, bounds=(0, 2.5)) #Defined in natural time, not normalised
    m.u[0].fix(0)
    def r_ca(_m, t):
        return _m.ca_in[t]/_m.d['t_f'] == _m.dudt[t]
    m.c_ca = pyo.Constraint(m.t, rule=r_ca)

    #Constraints
    R = 1
    def r_mb1(_m, t): 
        return _m.dcdt['A',t]/_m.d['t_f'] == -2*_m.p['k1']*pyo.exp(-_m.p['E1']/(R*_m.d['T']))*_m.c['A',t]**2
    def r_mb2(_m, t):
        return _m.dcdt['B',t]/_m.d['t_f'] == _m.p['k1']*pyo.exp(-_m.p['E1']/(R*_m.d['T']))*_m.c['A',t]**2 \
            - _m.p['k2']*pyo.exp(-_m.p['E2']/(R*_m.d['T']))*_m.c['B',t]
    def r_mb3(_m, t):
        return  _m.dcdt['C',t]/_m.d['t_f']  == _m.p['k2']*pyo.exp(-_m.p['E2']/(R*_m.d['T']))*_m.c['B',t]
    
    m.mb1 = pyo.Constraint(m.t, rule=r_mb1)
    m.mb2 = pyo.Constraint(m.t, rule=r_mb2)
    m.mb3 = pyo.Constraint(m.t, rule=r_mb3)
    
    m.g1 = pyo.Var(initialize=0)
    m.g2 = pyo.Var(initialize=0)

    def r_cqa1(_m):
        return _m.g1 == 1000*(0.8 -(_m.c['B',1] + 1E-2)/(_m.c['A', 1]+_m.c['B',1]+_m.c['C',1]+ 1E-2))
    def r_cqa2(_m):
        _m = _m.model()
        return _m.g2 == -(100*_m.c['B',1] - 20*_m.c['A',0] - 128*(_m.d['t_f'] + 30))/1E3
    m.c_g1 = pyo.Constraint(rule=r_cqa1)
    m.c_g2 = pyo.Constraint(rule=r_cqa2)

    m.cqa1 = pyo.Constraint(expr=m.g1 <= 0)
    m.cqa2 = pyo.Constraint(expr=m.g2 <= 0)

    m.c['A',0].fix(2E3) #mol/L
    m.c['B',0].fix(0)
    m.c['C',0].fix(0)

    m.var_input = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    m.var_input[m.ca_in] = {0: 0}

    #Configure the scenario options
    m._pyds_opts = ScenarioOpts(
                                [m.d, m.t, m.I], 
                                {m.p: random.normal(0,1)},
                                {m.cqa1: 1E3, m.cqa2: 1E2} , 
                                m.var_input,)
    
    return m

m = create_model()
m2 = build_scenarios(m, 1000)
print()


    

