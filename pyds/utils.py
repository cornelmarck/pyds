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
    
def _create_objective(m):
    def obj_rule(m):
        n_scenarios = prod(m.BFs)
        all_idx = product(*[range(i) for i in m.BFs])
        return 1/(n_scenarios)*sum(1-get_scenario(m, idx)._indicator_var for idx in all_idx)

    m.obj = Objective(rule=obj_rule, sense=maximize)

def get_final_scenarios(m):
    final_stage = m.n_stages - 1
    return [v for v in scenarios_at_stage(m, final_stage)]

def scenarios_at_stage(m, stage):
    """Get an interator which gives index, scenario 

    Args:
        m (ConcreteModel): The multi-stage model
        stage (int): target stage
    """
    if stage==0:
        return [m]

    BFs = m.BFs[0:stage+1]
    all_idx = product(*[range(i) for i in BFs])
    return [get_scenario(m, i) for i in all_idx]

def get_scenario(m, idx):
    """Get a single scenario by idx
    Args:
        m (ConcreteModel): Pyomo model
        idx (list or tuple of int): scenario index
    """
    if isinstance(idx, int):
        idx = (idx,)
    def recursive(obj, stage, idx):
        obj = obj.__getattribute__('Substage')[idx[stage-1]]
        if stage == len(idx):
            return obj
        else:
            return recursive(obj, stage+1, idx)
    return recursive(m, 1, idx)    
        
def create_EF(stage_rules, BFs):
    ef = ConcreteModel()
    ef.n_stages = len(BFs)+1
    ef.BFs = BFs
    stages_dict = {}

    if ef.n_stages == 1:
        ValueError('Multistage model must contain more than one stage')

    def recursive_block_rule(m, stage, n_stages, BFs, stages_dict):
        stage_rules[stage](m, stages_dict)
        if stage == n_stages-1: #Check whether last stage has been reached
            return
        stages_dict[stage] = m
        m.add_component('Substage_idx', RangeSet(0, BFs[stage]-1))
        m.add_component('Substage', Block(m.Substage_idx, 
            rule= lambda m: recursive_block_rule(m, stage+1, n_stages, BFs, stages_dict)))

    recursive_block_rule(ef, 0, ef.n_stages, ef.BFs, stages_dict)
    _create_objective(ef)
    return ef

def create_flattened_model(stage_rules):
    m = ConcreteModel()

    dummy_parent_models = [m]*(len(stage_rules)-1)
    for i in stage_rules:
        i(m, dummy_parent_models)
    return m          

def raw_ouput_writer(m):
    """Return a dict which contains the values of all variables, objectives and parameters. 

    Args:
        m (Pyomo model): [description]

    Returns:
        [dict]: [description]
    """
    output = {}
    for i in m.component_data_objects(Var, active=True):
        output[i.name] = i.value

    for i in m.component_data_objects(Objective, active=True):
        output[i.name] = i.expr()

    for i in m.component_objects(Param):
        if i.is_indexed():
            for k, v in i.items():
                name = '{}[{}]'.format(i.name, k)
                output[name] = v
        else:
            output[i.name] = i.value    
     
    return output


# def output_writer(m):
#     """Return a dict which contains the values of all variables, objectives and parameters. 

#     Args:
#         m ([type]): [description]

#     Returns:
#         [type]: [description]
#     """
#     output = {}
#     for i in m.component_objects(Var, active=True):
#         if i.is_indexed():
#             var = {}
#             for k in i.keys():
#                 var[k] = i[k].value
#             output[i.name] = var
#         else:
#             output[i.name] = i.value

#     for i in m.component_objects(Objective, active=True):
#         output[i.name] = i.expr()

#     for i in m.component_objects(Param):
#         if i.is_indexed():
#             for k, v in i.items():
#                 name = '{}[{}]'.format(i.name, k)
#                 output[name] = v
#         else:
#             output[i.name] = i.value  

#     return output




    


















