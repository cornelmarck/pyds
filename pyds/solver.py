from pyomo.core.base import param
from pyomo.core.base.transformation import TransformationFactory
from pyomo.opt import (SolverFactory,  TerminationCondition)
import numpy
import pickle
import os

import pyds.utils as utils

class Solver():
    def __init__(self, parent, config):
        self.parent = parent
        self.tee = config['tee']
        config.setdefault('use relaxation', True)
        self.use_relaxation = config['use relaxation']
        self.io_options = config['io options']
        self.io_options.setdefault('warmstart', True) #Required for solution!
        self.solve_trajectories = config['solve trajectories']
        self.warn_infeasible = config['warn infeasible']
        self.save_output = config['save output']
        self.save_solution_states = config['save solution state']
        self.no_infeasible = 0
        
        self.solver_obj = SolverFactory(config['name'])
        
    def solve(self):
        self.result = None
        self.model = self.parent.model

        self.output = {
            'relaxation': None,
            'trajectories': None,
        }
        if self.use_relaxation:
            TransformationFactory('core.relax_integer_vars').apply_to(self.model)
            #http://www.pyomo.org/blog/2015/1/8/accessing-solver
            self._solve_relaxation()
            if (self.result['Solver'].termination_condition == TerminationCondition.infeasible):
                self._set_infeasible_indicator_var()
                self.no_infeasible += 1 
                if self.warn_infeasible:    
                    print('Warning: Infeasible relaxation. Total number of infeasible solves: {}'.format(self.no_infeasible))
                return
        
        if self.solve_trajectories and self.save_output:
            self._solve_trajectories()
            self.output['trajectories'] = self._collect_output().copy()

    def _solve_relaxation(self):
        self._reset_indicator_var()
        self.result = self.solver_obj.solve(self.model, tee=self.tee, io_options=self.io_options)   
        if self.save_output:
            self.output['relaxation'] = self._collect_output().copy()
        self._fix_indicator_var()

    def _solve_trajectories(self):
        self.result = self.solver_obj.solve(self.model, tee=self.tee, io_options=self.io_options)
        if self.save_output:
            self.output['trajectories'] = self._collect_output().copy()

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
            'user time': self.result['Solver'][0]['User time'],
             #'result': self.result, #Disabled due to pickling weakref error
            'infeasible': self.result['Solver'].termination_condition == TerminationCondition.infeasible
        }
        container['objective'] = utils.parse_value(self.model, 'obj')
        all_idx = utils.get_all_idx(self.model.BFs)
        for i, idx in enumerate(all_idx):
            v = {}
            for stage, names in self.parent.output_map.items():
                scen = utils.get_scenario(self.model, idx[0:stage+1])
                for n in names:
                    v[n] = utils.parse_value(scen, n)
            container['data'].append(v.copy())
        return container

