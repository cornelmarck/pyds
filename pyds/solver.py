from pyomo.core.base import param
from pyds.utils import (create_EF, scenarios_at_stage,  get_final_scenarios, raw_ouput_writer)
from pyomo.opt import (SolverFactory, SolverStatus, TerminationCondition)
import numpy
import pickle
import os

class GamsSolver():
    def __init__(self, parent, solve_trajectories=True,
                save_solution_states=False, **io_options):
        self.parent = parent
        #self.solver_obj = SolverFactory('gams')

        self.save_solution_states = save_solution_states
        self.solve_trajectories = solve_trajectories
        self.warn_infeasible = True
        self.no_infeasible = 0

        self.options = io_options
        self.options.setdefault('solver', 'conopt')
        self.options.setdefault('tee', False)
        self.options.setdefault('warmstart', True)
        
    def solve(self):
        self.output = {}
        self.relaxed_result = None
        self.relaxed_result = None
        self.model = self.parent.model

        #http://www.pyomo.org/blog/2015/1/8/accessing-solver
        self._solve_relaxation()
        if self.relaxed_result.termination_condition == TerminationCondition.infeasible:
            self._set_indicator_var_infeasible()
            if self.warn_infeasible: 
                print('Warning: Infeasible relaxation solution')
            self.no_infeasible += 1 
            return

        if self.solve_trajectories:
            self._solve_trajectories()

    def _simulate_intial_state():
        pass

    def _solve_relaxation(self):
        self._reset_indicator_var()
        self.relaxed_result = self.solver_obj.solve(self.model, **self.options)
        if self.save_solution_states:
            self._generate_relaxation_output()

    def _solve_trajectories(self):
        self._fix_indicator_var()
        self.trajectories_result = self.solver_obj.solve(self.model, **self.options)
        if self.save_solution_states:
            self._generate_trajectories_output()

    #Set all indicator variables to 1 if the solution is infeasible
    def _set_infeasible_indicator_var(self): 
        for scen in get_final_scenarios(self.model):
            if hasattr(scen, '_indicator_var'):
                var = scen._indicator_var
                var.fix(1.0)
            else:
                continue

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

    def _generate_relaxation_output(self):
        self.output['relaxation'] = {'solver_result': self.relaxed_result, 
                'model_data': raw_ouput_writer(self.model),
                'indicator_vars': self._get_indicator_var_values()}

    def _generate_trajectories_output(self):
        self.output['trajectories'] = {'solver_result': self.solution_result, 
                'model_data': raw_ouput_writer(self.model),
                'indicator_vars': self._get_indicator_var_values()}
        self.output['input_values'] = self.parent.input_values
