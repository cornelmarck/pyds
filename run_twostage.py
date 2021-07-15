import mpisppy
import pyomo.environ as pyo
from pyomo import dae
from numpy import random
import mpisppy.scenario_tree as scenario_tree
from mpisppy.utils.sputils import create_EF

from pyds.builder import make_bigM_form, set_stochastic
from pyds.solver import resample_all_stoch_params

randomstream = random.default_rng()

def scenario_creator(scenario_name):
    m = pyo.ConcreteModel()
    #Sets
    m.t = dae.ContinuousSet(bounds=(0,1)) #Normalised time
    m.I = pyo.Set(initialize=['A', 'B', 'C'])

    #Parameters
    m.d = pyo.Param(['T', 't_f'], initialize={'t_f': 580, 'T': 285}, mutable=True)
    m.E1 = pyo.Param(initialize=0, mutable=True)
    m.k1 = pyo.Param(initialize=0, mutable=True)
    m.k2 = pyo.Param(initialize=0, mutable=True)
    m.E2 = pyo.Param(initialize=0, mutable=True)
    set_stochastic(m.E1, randomstream.normal, 2.5E3, 2.5E2)
    set_stochastic(m.E2, randomstream.normal, 5E3, 5E2)
    set_stochastic(m.k1, randomstream.normal, 6.409E-2, 6.409E-3)
    set_stochastic(m.k2, randomstream.normal, 9.938E3, 9.938E2)

    #Variables
    m.c = pyo.Var(m.I, m.t, bounds=(0, 1E8), initialize=0)
    m.dcdt = dae.DerivativeVar(m.c, wrt=m.t)
    m.u = pyo.Var(m.t, initialize=0)
    m.dudt = dae.DerivativeVar(m.u, initialize=0, wrt=m.t)
    m.ca_in = pyo.Var(m.t, bounds=(0, 0.001)) #Defined in natural time, not normalised
    m.u[0].fix(0)
    def r_ca(_m, t):
        return _m.ca_in[t]/_m.d['t_f'] == _m.dudt[t]
    m.c_ca = pyo.Constraint(m.t, rule=r_ca)

    #Constraints
    R = 1
    def r_mb1(_m, t): 
        return _m.dcdt['A',t]/_m.d['t_f'] == -2*_m.k1*pyo.exp(-_m.E1/(R*_m.d['T']))*_m.c['A',t]**2 + m.ca_in[t]
    def r_mb2(_m, t):
        return _m.dcdt['B',t]/_m.d['t_f'] == _m.k1*pyo.exp(-_m.E1/(R*_m.d['T']))*_m.c['A',t]**2 \
            - _m.k2*pyo.exp(-_m.E2/(R*_m.d['T']))*_m.c['B',t]
    def r_mb3(_m, t):
        return  _m.dcdt['C',t]/_m.d['t_f']  == _m.k2*pyo.exp(-_m.E2/(R*_m.d['T']))*_m.c['B',t]
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
    m.indicator_var = pyo.Var(bounds=(0,1), initialize=0) #Binary indicator var
    m.bigM_constant = pyo.Param([1,2], initialize={1: 1, 2: 1E3})
    m.cqa1 = make_bigM_form(pyo.Constraint(expr=m.g1<=0), m.indicator_var, m.bigM_constant[1])    
    m.cqa2 = make_bigM_form(pyo.Constraint(expr=m.g2<=0), m.indicator_var, m.bigM_constant[2])  

    m.c['A',0].fix(2E3) #mol/L
    m.c['B',0].fix(0)
    m.c['C',0].fix(0)

    def r_obj(model):
        return 1-model.indicator_var
    m.obj = pyo.Objective(rule=r_obj, sense=pyo.maximize)

    m._mpisppy_node_list = [
        scenario_tree.ScenarioNode(
            name='ROOT',
            cond_prob=1.0, #must be float
            stage=1,
            scen_name_list=None,
            cost_expression=None,
            nonant_list= [m.ca_in],
            scen_model=m
        )
    ]
    return m

scen_count = 1000
scenario_names = ['Scenario' + str(i) for i in range(scen_count)]

ef = create_EF(scenario_names, scenario_creator)
discretizer = pyo.TransformationFactory('dae.collocation')
discretizer.apply_to(ef, nfe=8, ncp=3)
resample_all_stoch_params(ef)




s = pyo.SolverFactory('baron')
s.solve(ef)

print()

