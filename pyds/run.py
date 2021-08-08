from pyds.utils import (create_EF, scenarios_at_stage, get_scenario, get_final_scenarios, raw_ouput_writer)
from pyds.solver import GamsSolver
from pyds.output_manager import OutputManager
from pyomo.opt import SolverFactory
from pyds.simulator import Simulator
import numpy
import os


class Manager():
    def __init__(self, stage_rules, input_bindings, model_transformation=None, \
            solver_io_options=None, simulate_initial_state=True, output_folder=None):
        self.stage_rules = stage_rules
        self.n_stages = len(stage_rules)
        self.input_bindings = input_bindings
        self.model_transformation = model_transformation
        self.simulate_initial_state = simulate_initial_state
        self.model = None
        
        if output_folder is None:
            self.output_folder = os.getcwd()
        self.output_manager = OutputManager(self.output_folder)

        self.solver = GamsSolver(self, solver_io_options)

        self.simulator=None
        if self.simulate_initial_state:
            flat_inputs={}
            for b, v in self.input_bindings.items():
                flat_inputs[b] = [n for s,n in v]
            self.simulator = Simulator(self.stage_rules, flat_inputs, self.model_transformation, simulator_package='casadi')

    def g(self, d, p):
        self.output = {}
        pass

    def scenario(self, idx):
        return get_scenario(self.model, idx)

    def update_input(self, input_values):
        """Set the input parameters at a stage input_values.
        Args:
            input_values (dict of tuple, list): keys are tuple (stage_id, parameter_name), values are 2d array of 
        """
        self.input_values = input_values

        for input_key, values in input_values.items():
            if not input_key in self.input_bindings:
                raise ValueError('Undefined input binding: ' + str(input_key))
                
            for param_idx, param_key in enumerate(self.input_bindings[input_key]):
                stage = param_key[0]
                local_param_name = param_key[1]
                for scen_idx, scen in enumerate(scenarios_at_stage(self.model, stage)): #Scenario
                    scen.component(local_param_name).set_value(values[scen_idx, param_idx])                             

    def write_output_to_disk(self):
        self.output_manager.write_output_to_disk()
    

class TwoStageManager(Manager):
    def __init__(self, stage_rules, input_bindings, model_transformation=None, \
            solver_io_options=None, simulate_initial_state=True, output_folder=None):
        super().__init__(stage_rules, input_bindings, model_transformation, \
            solver_io_options, simulate_initial_state, output_folder)

    def g(self, d, p):
        n_p = numpy.shape(p)[0]
        if self.model is None or self.model.BFs[0] != n_p:
            self._build_model(n_p)

        g_list = []
        #Iterate the design points d
        for i, d_point in enumerate(d):
            g_mat = numpy.empty((n_p, 1))
            self.update_input({'d': numpy.array([d_point]), 'p': p})
            if self.simulate_initial_state:
                self.simulator.simulate_all_scenarios(get_final_scenarios(self.model), self.input_values)  
            self.solver.solve()

            g_mat[:, 0] = -1*self.solver._get_indicator_var_values() #DEUS uses g>=0 inequality constraints, contrary to convention
            g_list.append(g_mat)
            self.output_manager.write_directly(self.solver.output)
        return g_list

    def _build_model(self, n_p):
        self.model = create_EF(self.stage_rules, [n_p])
        self.model_transformation(self.model) 
        

class ThreeStageManager(Manager):
    def __init__(self, stage_rules, input, solver_name, model_transformation=None, solver_options=None):
        super().__init__(stage_rules, input, solver_name, model_transformation, solver_options)
