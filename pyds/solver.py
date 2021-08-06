from pyomo.core.base import param
from pyds.utils import (create_EF, scenarios_at_stage, get_final_scenarios, raw_ouput_writer)
from pyomo.opt import SolverFactory
import numpy
import pickle
import os

class GamsSolver():
    def __init__(self, parent, io_options=None):
        self.parent = parent
        self.model = parent.model
        self.solver_obj = SolverFactory('gams')

        self.options = io_options
        self.options.setdefault('solver', 'conopt')
        self.options.setdefault('tee', False)
        self.options.setdefault('warmstart', True)
        

    def solve(self):
        self.output = {}

        self._reset_indicator_var()
        self.relaxed_result = self.solver.solve(self.model, **self.solver_options)
        self.output['relaxation'] = {'solver_result': self.relaxed_result, 
                'model_data': raw_ouput_writer(self.model),
                'indicator_vars': self._get_indicator_var_values()}
        
        #Solve model
        self._fix_indicator_var()
        self.solution_result = self.solver.solve(self.model, **self.solver_options)
        self.output['solution'] = {'solver_result': self.solution_result, 
                'model_data': raw_ouput_writer(self.model),
                'indicator_vars': self._get_indicator_var_values()}

        self.output['input_values'] = self.parent.input_values

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
    
    def _get_indicator_var_values(self):
        final_scenarios = get_final_scenarios(self.model)
        vars = numpy.empty(len(final_scenarios))
        
        for id, scen in enumerate(final_scenarios):
             if hasattr(scen, '_indicator_var'):
                 vars[id] = scen._indicator_var.value
        return vars


    