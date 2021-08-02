from mpisppy.utils.sputils import create_EF
from pyomo.core.base import (ConcreteModel, Block, RangeSet, Constraint, ConstraintList, Param, Var)
from pyomo.core.base.objective import Objective, maximize
from itertools import product
from math import prod

def add_BigMConstraint(model, name, bigM_constant, *args, **kwargs):
    """Obtain the big-M form of an (in)equality constraint.

    Args:
        constraint (pyomo.core.base.Constraint): The (in)equality constraint to be deactivated. 
            Constraint list input is not implemented yet. 
        indicator_var (pyomo.core.base.Var): Binary variable (also relaxed) which 0 if all 
            constraints are satisfied.
        bigM_constant (pyomo.core.base.Param): Big-M constant which is scaled to the magnitude 
            of the inequality constraint.
    Returns:
        :pyomo.core.base.Constraint | pyomo.core.base.ConstraintList | None : 
            The constraint in big-m form.
    """
    constraint = Constraint(*args, **kwargs)
    constraint.construct()
    if not hasattr(model, '_indicator_var'):
        model._indicator_var = Var(initialize=0, bounds=(0,1))
    indicator_var = model._indicator_var

    if (constraint.lower is None) & (constraint.upper is not None):
        cstr = Constraint(
            expr= constraint.body <= constraint.upper + bigM_constant*indicator_var)
    elif (constraint.lower is not None) & (constraint.upper is None):
        cstr = Constraint(
            expr= constraint.body >= constraint.lower - bigM_constant*indicator_var)

    #Equality constraint
    elif (constraint.lower is not None) & (constraint.upper is not None):
        lst = ConstraintList()
        lst.construct()
        lst.add(expr= constraint.body <= constraint.upper + bigM_constant*indicator_var)
        lst.add(expr= constraint.body >= constraint.lower - bigM_constant*indicator_var)
        cstr = lst
    else:
        cstr =  None
    
    model.add_component(name, cstr)
    
def create_objective(m):
    def obj_rule(m):
        n_scenarios = prod(m.BFs)
        all_idx = product(*[range(1, i+1) for i in m.BFs])
        return 1/(n_scenarios)*sum(get_scenario(m, idx)._indicator_var for idx in all_idx)

    m.obj = Objective(rule=obj_rule, sense=maximize)

def get_scenario(m, idx):
    """Get the scenario relocating with idx which is 1-based.

    Args:
        m (ConcreteModel): Pyomo model
        idx (list or tuple of int): scenario index
    """
    if not isinstance(idx, (list, tuple)):
        idx = [idx]

    def recursive(obj, stage, idx):
        if not stage < len(idx):
            return obj
        obj = obj.__getattribute__('Stage'+str(stage+1))[idx[stage]]
        return recursive(obj, stage+1, idx)
    return recursive(m, 0, idx)    
        
def create_flattened_model(stage_rules):
    m= ConcreteModel()
    for r in stage_rules:
        r(m, None)
    return m

def create_EF(stage_rules, BFs):
    ef = ConcreteModel()
    ef.n_stages = len(BFs)+1
    ef.BFs = BFs
    stages_dict = {}

    def recursive_block_rule(m, stage, n_stages, BFs, stages_dict):
        stage_rules[stage](m, stages_dict)
        if not stage < n_stages-1: #Check whether last stage has been reached
            return
        stages_dict[stage] = m
        m.add_component('Stage_idx', RangeSet(BFs[stage]))
        m.add_component('Stage'+str(stage+1), Block(m.Stage_idx, 
            rule= lambda m: recursive_block_rule(m, stage+1, n_stages, BFs, stages_dict)))

    recursive_block_rule(ef, 0, ef.n_stages, ef.BFs, stages_dict)
    create_objective(ef)
    return ef

def create_multistage_model_no_branches(scenario_creator):
    """Create multistage model without branches. Is used as a template model for model construction.

    Args:
        scenario_creator ([type]): [description]
        scenario_creator_kwargs ([type]): [description]
    """
    tmp = scenario_creator()
    if not hasattr(tmp, '_pyds_stages'):
        raise Exception('No _pyds_stages dict is found')
    stages_dict = tmp._pyds_stages

    model = ConcreteModel()

    return    

def initialize_empty_substages(root, n_stages):
    pass
    



    


















# OLD FUNCTIONS. DEUS INTERFACE WILL SET THE STOCHASTIC PARAMETERS EXTERNALLY
# class StochModel():
#     """Create the multi-stage stochastic model and associated solver functions.
#     This is a modified mpispy ExtensiveForm object in order to not require an MPI installation.

#     The scenario_creator function must return a pyomo model with the following objects defined at the root.
#         1. _mpisppy_node_list (list): The ScenarioNodes containting e.g. nonanticipative variables
#         2. _mpisppy_probability (float): The probability of a given scenario (1/n_scenarios) for two-stage model.
#         3. StochParam and inequality constraint objects. 
    
#         Do not create hierarchical models (i.e. Pyomo Blocks) in the scenario_creator.

#     Args:
#         ExtensiveForm ([type]): [description]
#     """
#     def __init__(self, options, all_scenario_names, scenario_creator, scenario_creator_kwargs):
#         self.options = options
#         self.all_scenario_names = all_scenario_names
#         self.scenario_creator = scenario_creator

#         self.t = 1
        
#         if "branching_factors" in self.options:
#             self.branching_factors = self.options["branching_factors"]
#         else:
#             self.branching_factors = [len(self.all_scenario_names)]

#         return
# def StochParam(sampling_func, *sampling_func_args, **param_kwargs):
#     param_kwargs.setdefault('mutable', True)
#     param_kwargs.setdefault('initialize', 0)
#     param = Param(**param_kwargs)
#     param.pyds_sampler = Sampler(param, sampling_func, *sampling_func_args)
#     return param

# class Sampler:
#     def __init__(self, param, func, *args, **kwargs):
#         """Store a resampling function and its arguments.
#         TODO: Input error checking

#         Args:
#             param (pyo.Param): Parent of sampler object
#             func (function)): Random samples generator function.
#                 Must return a a single float or dict {key, float}
#             args, kwargs: Additional arguments passed to func, e.g. a stream of random numbers
#         """
#         self.func = func
#         self.param = param
#         self.args = args
#         self.kwargs = kwargs
        
#     def resample(self):
#         """Resample and update the values of the stochastic parameter"""   
#         sample_value = self.__call__()
#         if isinstance(sample_value, dict):
#             for k,v in sample_value.items():
#                 self.param[k].value = v
#         else:
#             self.param.value = sample_value

#     def __call__(self):
#         return self.func(*self.args, **self.kwargs)

# def is_stochastic(obj):
#     return hasattr(obj, "pyds_stoch_sampler")

# def is_ineq_constraint(obj):
#     return hasattr(obj, "pyds_ineq_constraint")
