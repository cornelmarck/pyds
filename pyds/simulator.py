from pyds.utils import create_flattened_model, get_all_idx, get_scenario, load_input

from pyomo.common.collections.component_map import ComponentMap
from pyomo.core.expr.template_expr import IndexTemplate
from pyomo.dae import Simulator as PyomoSimulator
from pyomo.dae.diffvar import DAE_Error

import numpy as np

class Simulator:
    def __init__(self, parent, config):
        self.parent = parent

        self.stage_rules = parent.stage_rules
        self.model_transformation = parent.model_transformation
        self.enabled = config['enabled']
        self.package = config['package']
        self.suffix_name = config['suffix name']
        self.kwargs = config['kwargs']
        
        self.model = None
        self.simulator_obj = None
        if self.enabled:
            self._build_model()
            self.simulator_obj = PyomoSimulator(self.model, self.package)

    def simulate_all_scenarios(self, output_model, input_values):
        #Warning: This method only initializes time-varying variables, not input parameters
        if not self.enabled:
            return
        model = output_model
        BFs = model.BFs
        n_stages = model.n_stages   
        for s_idx in get_all_idx(BFs):
            input = {}
            input[0] = input_values[0][None, 0,:]
            for s in range(1,n_stages-1):
                input[s] = input_values[s][None, s_idx[s-1], :]
            self._simulate(input)
            self._export_trajectories_to_model(model, s_idx)

    def _simulate(self, input_values):
        load_input(self.model, self.parent.input_map, input_values)
        if self.suffix_name is not None:
            suffix = self.model.component(self.suffix_name)
            self.kwargs['varying_inputs'] = suffix
        self.simulator_obj.simulate(**self.kwargs)

    def _export_trajectories_to_model(self, model, scenario_idx=None):        
        #Adapted from pyomo dae.simulator.initialize_model()
        def stages_to_update(idx):
            a = np.array(idx)
            should_update_mask = [np.all(a[i+1:]==0) for i in range(len(a))] #empty slice of last i evaluates to True
            stages = [i for i in range(len(a)) if should_update_mask[i]]
            return [tuple(idx[:i+1]) for i in stages]
            
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
                
                local_name = v._base.local_name
                if scenario_idx is None:
                    model.component(v._base.local_name)[vidx].set_value(valinit[i])
                else:
                    for idx in stages_to_update(scenario_idx):
                        var = get_scenario(model, idx).component(v._base.local_name)
                        if var is not None:
                            var[vidx].set_value(valinit[i])
                        else:
                            continue

    def get_sim_result(self):
        return (self.simulator_obj._tsim, self.simulator_obj._simsolution)
   
    # def _update_input(self, input_values):
    #     for stage, value in input_values.items():
    #         if not stage in self.parent.input_map.keys():
    #             raise ValueError('Input is not defined for this problem')
    #         for idx, param_name in enumerate(self.parent.input_map[stage]):
    #             self.model.component(param_name).set_value(value[idx]) #TODO: Rewrite to 2d array convention
    
    def _build_model(self):
        self.model = create_flattened_model(self.stage_rules)
        if self.model_transformation is not None:
            self.model_transformation(self.model)

        