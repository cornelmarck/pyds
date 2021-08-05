from pyomo.core.base import param
from pyds.buildutils import (create_EF, scenarios_at_stage, get_scenario, get_final_scenarios)
from pyomo.opt import SolverFactory
import numpy


class Manager():
    def __init__(self, stage_rules, input_bindings, solver_name, model_transformation=None, solver_options=None):
        self.stage_rules = stage_rules
        self.n_stages = len(stage_rules)
        self.input_bindings = input_bindings
        self.model_transformation = model_transformation
        self.model = None

        self.solver = SolverFactory(solver_name)
        self.solver_options = solver_options
        self.solver_options.setdefault('solver', 'conopt')
        self.solver_options.setdefault('tee', False)
        self.solver_options.setdefault('warmstart', True)


    def g(self, d, p):
        pass

    def update_input(self, input_values):
        """Set the input parameters at a stage input_values.
        Args:
            input_values (dict of tuple, list): keys are tuple (stage_id, parameter_name), values are 2d array of 
        """
        for input_key, values in input_values.items():
            if not input_key in self.input_bindings:
                raise ValueError('Undefined input binding: ' + str(input_key))
                
            for param_idx, param_key in enumerate(self.input_bindings[input_key]):
                stage = param_key[0]
                local_param_name = param_key[1]
                for scen_idx, scen in enumerate(scenarios_at_stage(self.model, stage)): #Scenario
                    scen.component(local_param_name).set_value(values[scen_idx, param_idx])                             

    def solve(self):
        self._reset_indicator_var()
        self.relaxed_result = self.solver.solve(self.model, **self.solver_options)
        self._fix_indicator_var()
        return self.solver.solve(self.model, **self.solver_options)

    def _get_indicator_var_values(self):
        final_scenarios = get_final_scenarios(self.model)
        vars = numpy.empty(len(final_scenarios))
        
        for id, scen in enumerate(final_scenarios):
             if hasattr(scen, '_indicator_var'):
                 vars[id] = scen._indicator_var.value
        return vars


    def _fix_indicator_var(self):
        for scen in get_final_scenarios(self.model):
            if hasattr(scen, '_indicator_var'):
                var = scen._indicator_var
                if var.value==0:
                    var.fix(0)
                else:
                    var.fix(1)
            else:
                continue
    
    def _reset_indicator_var(self):
        for scen in get_final_scenarios(self.model):
            if hasattr(scen, '_indicator_var'):
                var = scen._indicator_var
                var.unfix()
                var.set_value(0)
            else:
                continue

class TwoStageManager(Manager):
    def __init__(self, stage_rules, input, solver_name, model_transformation=None, solver_options=None):
        super().__init__(stage_rules, input, solver_name, model_transformation, solver_options)

    def g(self, d, p):
        n_p = numpy.shape(p)[0]
        if self.model is None or self.model.BFs[0] != n_p:
            self._build_model(n_p)

        g_list = []
        #Iterate the design poi]\nts d
        for i, d_point in enumerate(d):
            g_mat = numpy.empty((n_p, 1))
            self.update_input({'d': numpy.array([d_point]), 'p': p})
            self.solve()
            g_mat[:, 0] = self._get_indicator_var_values()
            g_list.append(g_mat)
        return g_list


    def _build_model(self, n_p):
        self.model = create_EF(self.stage_rules, [n_p])
        self.model_transformation(self.model) 

class ThreeStageManager(Manager):
    def __init__(self, stage_rules, input, solver_name, model_transformation=None, solver_options=None):
        super().__init__(stage_rules, input, solver_name, model_transformation, solver_options)


