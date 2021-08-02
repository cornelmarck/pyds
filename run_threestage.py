import mpisppy
import pyomo.environ as pyo
from pyomo import dae
from numpy import random
from pyds.buildutils import StochParam, BigMConstraint
import mpisppy.scenario_tree as scenario_tree
import mpisppy.utils.sputils as sputils

randomstream = random.default_rng()

def model_creator(scenario_name, T, t_f):
    m = pyo.ConcreteModel()
    #Sets
    m.t = dae.ContinuousSet(bounds=(0,1)) #Normalised time
    m.I = pyo.Set(initialize=['A', 'B', 'C'])

    #Parameters
    m.T = pyo.Param(initialize=T, mutable=True)
    m.t_f = pyo.Param(initialize=t_f, mutable=True)
    m.E1 = StochParam(randomstream.normal, 2.5E3, 2.5E2)
    m.E2 = StochParam(randomstream.normal, 5E3, 5E2)
    m.k1 = StochParam(randomstream.normal, 6.409E-2, 6.409E-3)
    m.k2 = StochParam(randomstream.normal, 9.938E3, 9.938E2)

    #Variables
    m.c = pyo.Var(m.I, m.t, bounds=(0, 1E8), initialize=0)
    m.dcdt = dae.DerivativeVar(m.c, wrt=m.t)
    m.u = pyo.Var(m.t, initialize=0)
    m.dudt = dae.DerivativeVar(m.u, initialize=0, wrt=m.t)
    m.ca_in = pyo.Var(m.t, initialize=0, bounds=(0, 0.001)) #Defined in natural time, not normalised
    m.u[0].fix(0)
    def r_ca(_m, t):
        return _m.ca_in[t]/_m.t_f == _m.dudt[t]
    m.c_ca = pyo.Constraint(m.t, rule=r_ca)

    #Constraints
    R = 1
    def r_mb1(_m, t): 
        return _m.dcdt['A',t]/_m.t_f == -2*_m.k1*pyo.exp(-_m.E1/(R*_m.T))*_m.c['A',t]**2 + m.ca_in[t]
    def r_mb2(_m, t):
        return _m.dcdt['B',t]/_m.t_f == _m.k1*pyo.exp(-_m.E1/(R*_m.T))*_m.c['A',t]**2 \
            - _m.k2*pyo.exp(-_m.E2/(R*_m.T))*_m.c['B',t]
    def r_mb3(_m, t):
        return  _m.dcdt['C',t]/_m.t_f  == _m.k2*pyo.exp(-_m.E2/(R*_m.T))*_m.c['B',t]
    m.mb1 = pyo.Constraint(m.t, rule=r_mb1)
    m.mb2 = pyo.Constraint(m.t, rule=r_mb2)
    m.mb3 = pyo.Constraint(m.t, rule=r_mb3)

    m.g1 = pyo.Var(initialize=0)
    m.g2 = pyo.Var(initialize=0)
    def r_cqa1(_m):
        return _m.g1 == 1000*(0.8 -(_m.c['B',1] + 1E-2)/(_m.c['A', 1]+_m.c['B',1]+_m.c['C',1]+ 1E-2))
    def r_cqa2(_m):
        _m = _m.model()
        return _m.g2 == -(100*_m.c['B',1] - 20*_m.c['A',0] - 128*(_m.t_f + 30))/1E3
    m.c_g1 = pyo.Constraint(rule=r_cqa1)
    m.c_g2 = pyo.Constraint(rule=r_cqa2)

    m.pyds_indicator_var = pyo.Var(bounds=(0,1), initialize=0) #Binary indicator var
    m.bigM_constant = pyo.Param([1,2], initialize={1: 1E3, 2: 1E5})
    m.cqa1 = BigMConstraint(m.pyds_indicator_var, m.bigM_constant[1], expr=m.g1<=0)  
    m.cqa2 = BigMConstraint(m.pyds_indicator_var, m.bigM_constant[2], expr=m.g2<=0)
    
    m.c['A',0].fix(2E3) #mol/L
    m.c['B',0].fix(0)
    m.c['C',0].fix(0)

    def r_obj(model):
        return 1-model.pyds_indicator_var
    m.obj = pyo.Objective(rule=r_obj, sense=pyo.maximize)
    return m

def make_nodes_for_scenario(m, branching_factors, scenario_number):
    """Initialize the pysp nodes for a single scenario.

    Args:
        m (ConcreteModel): [description]
        branching_factors (list of int): Branching factor for each stage.
        scenario_number (int): Scenario number.

    Returns:
        [type]: [description]
    """
    node_name = "ROOT_"+str((scenario_number-1) // branching_factors[0])
    node_list = [scenario_tree.ScenarioNode(
            name='ROOT',
            cond_prob=1.0, #must be float
            stage=1,
            scen_name_list=None,
            cost_expression=None,
            nonant_list= [m.T, m.t_f],
            scen_model=m
            ),
        
        scenario_tree.ScenarioNode(
            name=node_name,
            cond_prob=1.0/branching_factors[1],
            stage=2,
            scen_name_list=None,
            cost_expression=None,
            nonant_list= [m.ca_in],
            scen_model=m
        )
    ]
    m._mpisppy_node_list = node_list

def scenario_creator(scenario_name, T, t_f, branching_factors=None):
    m = model_creator(scenario_name, T, t_f) 
    scenario_number = sputils.extract_num(scenario_name)
    make_nodes_for_scenario(m, branching_factors, scenario_number)

BFs = [3, 3]
scen_count = BFs[0]*BFs[1]
scenario_names = ['Scenario' + str(i+1) for i in range(scen_count)]

ef = sputils.create_EF(scenario_names, scenario_creator,
    scenario_creator_kwargs={'branching_factors': BFs, 'T': 285, 't_f': 290})
discretizer = pyo.TransformationFactory('dae.collocation')
discretizer.apply_to(ef, nfe=8, ncp=3)


#s = pyo.SolverFactory('baron')

print()
    

