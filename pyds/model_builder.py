from pyomo.core.base.component import Component
from pyomo.core import (Param, Suffix, ConcreteModel, RangeSet, Block)
#from pyomo.dae import ContinuousSet

from weakref import ref as weakref_ref
from inspect import signature

class ScenarioOpts:
    def __init__(self, first_stage=None, stoch_param=None, cqa=None, varying_input=None):
        """Specify the stocastic scenario options and check user-provided input. 

        Args:
            first_stage (list of Component, optional): Components that are shared among scenarios. Defaults to None.
            stoch_param (dict of (Param, function), optional): Mapping of stochastic parameters to corresponding sampling functions. \
                Sampling function cannot take any arguments. Defaults to None.
            varying_input (Suffix, optional): Time-varying control trajectories for simulation. Defaults to None
        """
        self._first_stage = first_stage
        self._stoch_param = stoch_param
        self._varying_input = varying_input
        self._cqa = cqa

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

class ModelBuilder:
    def __init__(self, input_model):
        self._input_model = input_model
        self._template_constructed = False
        self._scenarios_constructed = False
        self._bigm_applied = False

        self._model = None
        self.duration = {}
    
    def build_scenarios(self, n_scenarios):
        self._n_scenarios = n_scenarios
        m_tmp = self._input_model.clone()
        if not isinstance(m_tmp._pyds_opts, ScenarioOpts):
            raise ValueError("PyDS scenario options are not found.")
        self._model = ConcreteModel()

        #Build template scenario
        self._model._template_scenario = m_tmp
        for v in m_tmp._pyds_opts.first_stage:
            self._model._template_scenario.del_component(v)
            self._model.add_component(v.name, v)
        self._template_constructed = True

        #Duplicate the template scenario N times
        def copy_template(block):
            template = block.model()._template_scenario.clone()
            block.transfer_attributes_from(template)
            return
        self._model.scenario_idx = RangeSet(n_scenarios)
        self._model.scenario = Block(self._model.scenario_idx, rule=copy_template)
        self._model._template_scenario.deactivate()

    def apply_bigm(self):
        bigM = self._model.bigM = Block()
        bigM.indicator_var = pyo.Var(self._model.I, initialize=0, bounds=(0,1))
        bigM.constant = Param(initialize=bigM_value, mutable=True)
        
        def r_convert_to_bigM(block, I):
            def convert(constr, _i):
                b = constr.model().bigM

                if (constr._lower is None) & (constr._upper is not None):
                    return pyo.Constraint(expr= \
                        constr._body <= constr._upper + b.constant*b.indicator_var[_i])
            
                elif (constr._lower is not None) & (constr._upper is None):
                    return pyo.Constraint(expr= \
                        constr._body >= constr._lower - b.constant*b.indicator_var[_i])

                elif (constr._lower is not None) & (constr._upper is not None):
                    lst = pyo.ConstraintList()
                    lst.construct()
                    lst.add(expr= constr._body <= constr._upper + b.constant*b.indicator_var[_i])
                    lst.add(expr= constr._body >= constr._lower - b.constant*b.indicator_var[_i])
                    return lst
                return

            for name in constraint_names:
                constraint = block.model().scenarios[I].component(name)
                if constraint is None:
                    continue
                block.add_component(name, convert(constraint, I))
                constraint.deactivate()

            return

        bigM.constraints = pyo.Block(self._model.I, rule=r_convert_to_bigM)
        self._construct_SF_objective()

    def _construct_SF_objective(self):
        def r_obj(_m):
            return (1/_m.n_scenarios)*sum(1-_m.bigM.indicator_var[i] for i in _m.I)
        self._model.obj = pyo.Objective(rule=r_obj, sense=pyo.maximize)
        return


def build_scenarios(input_model, n_scenarios):
    m_tmp = input_model.clone()
    if not isinstance(input_model._pyds_opts, ScenarioOpts):
        raise ValueError("Scenario options are not found.")
    opts = m_tmp._pyds_opts
    
    #Create template scenario
    m = ConcreteModel()
    m._template_scenario = m_tmp
    for v in opts.first_stage:
        m._template_scenario.del_component(v)
        m.add_component(v.name, v)

    #Duplicate the template scenario N times
    def copy_template(block):
        template = block.model()._template_scenario.clone()
        block.__setstate__(template.__getstate__())
        return

    m.scenario_idx = RangeSet(n_scenarios)
    m.scenario = Block(m.scenario_idx, rule=copy_template)
    m._template_scenario.deactivate()

    return m


"""TODO: 
        @property
    def time(self):
        return self._time
        
    @time.setter
    def time(self, value):
        if isinstance(value, ContinuousSet):
            self._time = weakref_ref(value)
        elif value is None:
            self._time = None
        else:
            raise TypeError("Invalid continous set passed.")

    @property
    def input_param(self):
        return self._input_param
    
    @input_param.setter
    def input_param(self, value):
        self._input_param = []
        #Convert to list if single object is passed
        if not isinstance(value, list):
            value = list(value)
        elif value is None:
            return
        
        for e in value:
            if not isinstance(e, Parameter):
                raise TypeError("Input parameter must be valid type")
            elif not e._ismutable:
                raise ValueError("Input parameter must be mutable")  
            self._first_stage.append(weakref_ref(e))
"""            

#if len(list(m.component_data_objects(Block, active=True, descend_into=False))):
#    raise Exception("The builder cannot use hierarchical models")