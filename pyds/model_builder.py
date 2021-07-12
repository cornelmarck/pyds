from pyomo.core.base.component import Component
from pyomo.core import (Param, Var, Constraint, ConstraintList, \
    Suffix, ConcreteModel, RangeSet, Set, Block, Objective, maximize)
from numpy.random import default_rng
#from pyomo.dae import ContinuousSet
from weakref import ref as weakref_ref
from inspect import signature

class ModelBuilder:
    def __init__(self, input_model):
        self._input_model = input_model
        self._model = None
        self._opts = None
        self._duration = {} #Time measurements

    def build(self, n_scenarios):
        self.build_scenarios(n_scenarios)
        self.apply_bigm()
        return self._model

    def build_scenarios(self, n_scenarios):
        m_tmp = self._input_model.clone()
        if not isinstance(m_tmp._pyds_opts, ScenarioOpts):
            raise ValueError("PyDS scenario options are not found.")
        self._opts = m_tmp._pyds_opts
        
        #Initialize stochastic parameters
        for k,v in self._opts.stoch_param.items():
            k._stoch_func = v

        #Build template scenario
        self._model = ConcreteModel()
        self._model._template_scenario = m_tmp
        for v in m_tmp._pyds_opts.first_stage:
            self._model._template_scenario.del_component(v)
            self._model.add_component(v.name, v)

        #Duplicate the template scenario N times
        def copy_template(block):
            template = block.model()._template_scenario.clone()
            block.transfer_attributes_from(template)
            return
        self._model.n_scenarios = n_scenarios
        self._model.scenario_idx = RangeSet(n_scenarios)
        self._model.scenario = Block(self._model.scenario_idx, rule=copy_template)
        self._model._template_scenario.deactivate()

    def apply_bigm(self):
        bigM = self._model.bigM = Block()
        bigM.cqa_idx = Set([i.local_name for i in self._opts.cqa.keys()])
        bigM.indicator_var = Var(self._model.scenario_idx, initialize=0, bounds=(0,1))
        bigM.constant = Param(initialize=self._opts.cqa, mutable=True)
        
        def r_convert_to_bigM(block, scenario_idx):
            def convert(constr, cqa_idx, scen_idx ):
                b = constr.model().bigM
                if (constr._lower is None) & (constr._upper is not None):
                    return Constraint(expr= \
                        constr._body <= constr._upper + b.constant[cqa_idx]*b.indicator_var[scen_idx])
                elif (constr._lower is not None) & (constr._upper is None):
                    return Constraint(expr= \
                        constr._body >= constr._lower - b.constant[cqa_idx]*b.indicator_var[scen_idx])
                elif (constr._lower is not None) & (constr._upper is not None):
                    lst = ConstraintList()
                    lst.construct()
                    lst.add(expr= constr._body <= constr._upper + b.constant[cqa_idx]*b.indicator_var[scen_idx])
                    lst.add(expr= constr._body >= constr._lower - b.constant[cqa_idx]*b.indicator_var[scen_idx])
                    return lst
            m = block.model()
            for cqa_idx in m.bigM.cqa_idx:
                constraint = m.scenario[scenario_idx].component(cqa_idx)
                if constraint is None:
                    continue
                block.add_component(cqa_idx, convert(constraint, cqa_idx, scenario_idx))
                constraint.deactivate()

            return

        bigM.constraints = Block(self._model.scenario_idx, rule=r_convert_to_bigM)
        self._construct_SF_objective()

    def _construct_SF_objective(self):
        def r_obj(_m):
            return (1/_m.n_scenarios)*sum(1-_m.bigM.indicator_var[i] for i in _m.scenario_idx)
        self._model.obj = Objective(rule=r_obj, sense=maximize)
        return

class StochFunc:
    def __init__(self, func, deterministic, *args, **kwargs):
        """Generator function for given stochastic parameter.
        Args:
            func ([type]): Function taking No or kwargs as argument
            loc (Float or dict )))
        """
        self.func = func
        self.deterministic = deterministic
        self._args = args
        self._kwargs = kwargs

        if isinstance(deterministic, dict):
            self.is_indexed = True
        else:
            self.is_indexed = False

    def __call__(self):
        return self.func(*self._args, **self._kwargs)


class ScenarioOpts:
    def __init__(self, first_stage=None, stoch_param=None, cqa_constraint=None, varying_input=None):
        """Specify the stocastic scenario options and check user-provided input. 
        Args:
            first_stage (list of Component, optional): Components that are shared among scenarios. Defaults to None.
            stoch_param (dict of (Param, function), optional): Mapping of stochastic parameters to corresponding sampling functions. \
                Sampling function cannot take any arguments. Defaults to None.
            cqa
            varying_input (Suffix, optional): Time-varying control trajectories for simulation. Defaults to None
        """
        self._first_stage = first_stage
        self._stoch_param = stoch_param
        self._varying_input = varying_input
        self.cqa = cqa_constraint #TODO: Make getters - setters for attribute

    @property
    def first_stage(self):
        return self._first_stage
    @first_stage.setter
    def first_stage(self, value):
        self._first_stage = []
        #Convert to list if single object is passed
        if not isinstance(value, list):
            value = list(value)
        elif value is None:
            return
        
        for e in value:
            if isinstance(e, Component):
                self._first_stage.append(weakref_ref(e))
            else:
                raise TypeError("First stage component must be a valid Pyomo component")

    @property
    def stoch_param(self):
        return self._stoch_param
    @stoch_param.setter
    def stoch_param(self, value):
        self._stoch_param = {}
        if value is None:
            return
        elif not isinstance(value, dict):
             raise TypeError("Stochastic parameters must be passed in a Python dict.") 
        for (key, f) in value.items():
            if not (isinstance(key, Param) and isinstance (f, function)):
                raise ValueError("Invalid stochastic parameter combination passed")
            elif len(signature(f)) != 0:
                raise ValueError("Sampling function should not take any arguments. \
                        Reformulate e.g. using lambda functions.")
            elif not key._mutable:
                raise ValueError("The stochastic parameter must be mutable.")
            self._stoch_param[key] = f

    @property
    def varying_input(self):
        return self._varying_input
    @varying_input.setter
    def varying_input(self, value):
        if not isinstance(value, Suffix):
            raise TypeError("Time-varying input must be a Pyomo Suffix")
        self._varying_input = value
