from pyds.solver import Solver
from pyds.output_manager import OutputManager
import pyds.utils as utils
from pyomo.opt import SolverFactory
from pyds.simulator import Simulator
import numpy as np
import os


class Manager():
    def __init__(self, config):
        problem = config['problem']

        self.stage_rules = problem['stage rules']
        self.model_transformation = problem['model transformation']
        self.n_stages = len(self.stage_rules)
        self.input_map = problem['input map']
        self.output_map = problem['output map']
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

    def _build_model(self, BFs):
        self.model = utils.create_EF(self.stage_rules, BFs)
        self.model_transformation(self.model) 

    def scenario(self, idx):
        return utils.get_scenario(self.model, idx)                           

    def write_output_to_disk(self):
        self.output_manager.write_output_to_disk()
    

class TwoStageManager(Manager):
    def __init__(self, config):
        super().__init__(config)

    def g(self, d, p):
        #Rebuild the model and simulator if the BFs have changed
        n_p = np.shape(p)[0]
        if self.model is None or self.model.BFs[0] != n_p:
            self._build_model([n_p])

        g_list = []
        #Iterate the design points d
        for i, d_point in enumerate(d):
            g_mat = np.empty((n_p, 1))
            input = {0: np.array([d_point]), 1: p}
            utils.load_input(self.model, self.input_map, input)
            self.output_manager.add_input(input)
            self.simulator.simulate_all_scenarios(self.model, input)  
            self.solver.solve()
            self.output_manager.add_solver_solution(self.solver.output)
            g_mat[:, 0] = -self.solver._get_indicator_var_values() #DEUS uses g>=0 inequality constraints, contrary to convention
            g_list.append(g_mat)
            self.output_manager.write_data_to_disk()
            self.output_manager.clear_buffer()
        return g_list


        
class ThreeStageManager(Manager):
    def __init__(self, config):
        super().__init__(config)
    def g(self, d, p):
        #Rebuild the model and simulator if the BFs have changed
        n_d, d_dim = np.shape(d)
        n_p, p_dim = np.shape(p)
        BFs = np.array([n_d, n_p])

        if self.model is None or np.any(self.model.BFs != BFs):
            self._build_model(BFs)

        g_list = []
        #Iterate the design points d
        input = {1: d, 2: p}
        utils.load_input(self.model, self.input_map, input)
        self.simulator.simulate_all_scenarios(self.model, input)  
        self.solver.solve()
        ind_vars = self.solver._get_indicator_var_values()
        for i, d_point in enumerate(d):
            g_mat = np.empty((n_p, 1))
            g_mat[:, 0] = -self.solver._get_indicator_var_values() #DEUS uses g>=0 inequality constraints, contrary to convention
            g_list.append(g_mat)
            #self.output_manager.write_directly(self.solver.output)
        return g_list

