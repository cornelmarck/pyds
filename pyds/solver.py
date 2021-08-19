from pyomo.core.base import param
from pyomo.opt import (SolverFactory, TerminationCondition)
import numpy
import pickle
import os

import pyds.utils as utils

class Solver():
    def __init__(self, parent, config):
        self.parent = parent
        self.tee = config['tee']
        self.io_options = config['io options']
        self.io_options.setdefault('warmstart', True) #Required for solution!
        self.solve_trajectories = config['solve trajectories']
        self.warn_infeasible = config['warn infeasible']
        self.save_output = config['save output']
        self.save_solution_states = config['save solution state']
        self.no_infeasible = 0
        
        self.solver_obj = SolverFactory(config['name'])
        
    def solve(self):
        self.relaxed_result = None
        self.relaxed_result = None
        self.model = self.parent.model

        self.output = {
            'relaxation': None,
            'solution': None,
        }

        #http://www.pyomo.org/blog/2015/1/8/accessing-solver
        self._solve_relaxation()
        self.is_infeasible = (self.relaxed_result['Solver'].termination_condition == TerminationCondition.infeasible)
        if self.is_infeasible:
            self._set_infeasible_indicator_var()
            self.no_infeasible += 1 
            if self.warn_infeasible:    
                print('Warning: Infeasible relaxation. Total number of infeasible solves: {}'.format(self.no_infeasible))
            self.no_infeasible += 1 
            return
        if self.save_output:
            self.output['relaxation'] = self._collect_output().copy()

        if self.solve_trajectories and self.save_output:
            self._solve_trajectories()
            self.output['solution'] = self._collect_output().copy()

    def _solve_relaxation(self):
        self._reset_indicator_var()
        self.relaxed_result = self.solver_obj.solve(self.model, tee=self.tee, io_options=self.io_options)
        if self.save_solution_states:
            self._generate_relaxation_output()

    def _solve_trajectories(self):
        self._fix_indicator_var()
        self.trajectories_result = self.solver_obj.solve(self.model, self.tee, io_options=self.io_options)
        if self.save_solution_states:
            self._generate_trajectories_output()

    #Set all indicator variables to 1 if the solution is infeasible
    def _set_infeasible_indicator_var(self): 
        for scen in utils.get_final_scenarios(self.model):
            if hasattr(scen, '_indicator_var'):
                var = scen._indicator_var
                var.fix(1.0)
            else:
                continue

    def _fix_indicator_var(self):
        for scen in utils.get_final_scenarios(self.model):
            if hasattr(scen, '_indicator_var'):
                var = scen._indicator_var
                if var.value==0:
                    var.fix(0)
                else:
                    var.fix(1)
            else:
                continue
    
    def _reset_indicator_var(self):
        for scen in utils.get_final_scenarios(self.model):
            if hasattr(scen, '_indicator_var'):
                var = scen._indicator_var
                var.unfix()
                var.set_value(0)
            else:
                continue
    
    def _get_indicator_var_values(self):
        final_scenarios = utils.get_final_scenarios(self.model)
        vars = numpy.empty(len(final_scenarios))
        
        for id, scen in enumerate(final_scenarios):
             if hasattr(scen, '_indicator_var'):
                 vars[id] = scen._indicator_var.value
        return vars

    def _collect_output(self):
        container = {
            'data': [],
            'objective': None,
            'infeasible': self.is_infeasible
        }
        container['objective'] = utils.parse_value(self.model, 'obj')
        all_idx = utils.get_all_idx(self.model.BFs)
        for idx in all_idx:
            for stage, names in self.parent.output_map.items():
                scen = utils.get_scenario(idx[:stage+1])
                v = {}
                for n in names:
                    v[n] = utils.parse_value(scen, n)
            container['data'].append(v.copy())
        return container


#Deprecated. Saving the full solution state is too slow. 
    def _generate_relaxation_output(self):
        self.output['relaxation'] = {'solver_result': self.relaxed_result, 
                'model_data': utils.raw_ouput_writer(self.model),
                'indicator_vars': self._get_indicator_var_values()}

    def _generate_trajectories_output(self):
        self.output['trajectories'] = {'solver_result': self.solution_result, 
                'model_data': utils.raw_ouput_writer(self.model),
                'indicator_vars': self._get_indicator_var_values()}
        self.output['input_values'] = self.parent.input_values
