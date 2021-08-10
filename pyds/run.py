from pyds.utils import (create_EF, load_input, scenarios_at_stage, get_scenario, get_final_scenarios, raw_ouput_writer)
from pyds.solver import Solver
from pyds.output_manager import OutputManager
from pyomo.opt import SolverFactory
from pyds.simulator import Simulator
import numpy
import os


class Manager():
    def __init__(self, config):
        problem = config['problem']

        self.stage_rules = problem['stage rules']
        self.model_transformation = problem['model transformation']
        self.n_stages = len(self.stage_rules)
        self.input_map = problem['input map']
        
        if problem['output folder'] is None:
            self.output_folder = os.getcwd()
        else:
            self.output_folder = problem['output folder']
        self.output_manager = OutputManager(self.output_folder)
        self.solver = Solver(self, config['solver'])
        self.simulator = Simulator(self, config['simulator'])
        self.model = None

    def g(self, d, p):
        self.output = {}
        pass

    def scenario(self, idx):
        return get_scenario(self.model, idx)                           

    def write_output_to_disk(self):
        self.output_manager.write_output_to_disk()
    

class TwoStageManager(Manager):
    def __init__(self, config):
        super().__init__(config)

    def g(self, d, p):
        #Rebuild the model and simulator if the BFs have changed
        n_p = numpy.shape(p)[0]
        if self.model is None or self.model.BFs[0] != n_p:
            self._build_model(n_p)

        g_list = []
        #Iterate the design points d
        for i, d_point in enumerate(d):
            g_mat = numpy.empty((n_p, 1))
            input = {0: numpy.array([d_point]), 1: p}
            load_input(self.model, self.input_map, input)
            self.simulator.simulate_all_scenarios(self.model, input)  
            self.solver.solve()
            g_mat[:, 0] = -self.solver._get_indicator_var_values() #DEUS uses g>=0 inequality constraints, contrary to convention
            g_list.append(g_mat)
            #self.output_manager.write_directly(self.solver.output)
        return g_list

    def _build_model(self, n_p):
        self.model = create_EF(self.stage_rules, [n_p])
        self.model_transformation(self.model) 
        

class ThreeStageManager(Manager):
    def __init__(self, stage_rules, input, solver_name, model_transformation=None, solver_options=None):
        super().__init__(stage_rules, input, solver_name, model_transformation, solver_options)
