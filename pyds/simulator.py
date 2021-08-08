from pyds.utils import create_flattened_model

from pyomo.common.collections.component_map import ComponentMap
from pyomo.core.expr.template_expr import IndexTemplate
from pyomo.dae import Simulator as PyomoSimulator
from pyomo.dae.diffvar import DAE_Error

import numpy as np

class Simulator:
    def __init__(self, stage_rules, input_names, model_transformation, simulator_package='scipy', **simulator_kwargs):
        self.stage_rules = stage_rules
        self.model_transformation = model_transformation
        self.input_names = input_names
        self.simulator_package = simulator_package
        self.simulator_kwargs = simulator_kwargs

        self.model = None
        self._build_model()

        self.simulator_obj = PyomoSimulator(self.model, self.simulator_package)

    def simulate_all_scenarios(self, scenario_list, input_values):
        for idx, scen in enumerate(scenario_list):
            #Generate inputs for single scenario
            new_values = {}
            for input_key, values in input_values.items():
                new_values[input_key] = values[idx, :]
            self.simulate(new_values)
            self.export_trajectories_to(scen)


    def simulate(self, input_values):
        self._update_input(input_values)
        self.simulator_obj.simulate(**self.simulator_kwargs)

    def export_trajectories_to(self, model):
        #Adapted from pyomo dae.simulator.initialize_model()
        tvals = list(self.simulator_obj._contset)

        # Build list of state and algebraic variables
        # that can be initialized
        initvars = self.simulator_obj._diffvars + self.simulator_obj._simalgvars

        for idx, v in enumerate(initvars):
            for idx2, i in enumerate(v._args):
                    if type(i) is IndexTemplate:
                        break
            valinit = np.interp(tvals, self.simulator_obj._tsim,
                                self.simulator_obj._simsolution[:, idx])
            for i, t in enumerate(tvals):
                vidx = tuple(v._args[0:idx2]) + (t,) + \
                       tuple(v._args[idx2 + 1:])
                model.component(v._base.local_name)[vidx].set_value(valinit[i])

    def get_sim_result(self):
        return (self.simulator_obj._tsim, self.simulator_obj._simsolution)
   
    def _update_input(self, input_values):
        for name, value in input_values.items():
            if not name in self.input_names:
                raise ValueError('Input {} is not defined in the simulator dict'.format(name))
            for idx, param_name in enumerate(self.input_names[name]):
                self.model.component(param_name).set_value(value[idx])
    
    def _build_model(self):
        self.model = create_flattened_model(self.stage_rules)
        if self.model_transformation is not None:
            self.model_transformation(self.model)

        