from pyomo.core.base import Param
from mpisppy.opt.ef import ExtensiveForm

class StochModel(ExtensiveForm):
    """Create the multi-stage stochastic model and associated solver functions.

    The scenario_creator function must return a pyomo model with the following objects defined at the root.
        1. _mpisppy_node_list (list): The ScenarioNodes containting e.g. nonanticipative variables
        2. _mpisppy_probability (float): The probability of a given scenario (1/n_scenarios) for two-stage model.
        3. StochParam and inequality constraint objects. 
    
        Do not create hierarchical models (i.e. use pyomo Blocks) in the scenario_creator.

    Args:
        ExtensiveForm ([type]): [description]
    """
    def __init__(self, options, all_scenario_names, scenario_creator, scenario_creator_kwargs):
        super().__init__(options, all_scenario_names, scenario_creator, scenario_creator_kwargs)

        return

    def solve(self, resample=True, solver_options=None, tee=False):
        self.resample_all_stoch_params()
        self._reset_indicator_var()
        self.relaxed_result = self.solve_extensive_form(solver_options, tee)
        self._fix_indicator_var()
        return self.solve_extensive_form(solver_options, tee)


    def resample_all_stoch_params(self):
        all_params = self.ef.component_objects(Param, descend_into=True)
        for p in all_params:
            if hasattr(p, 'pyds_sampler'):
                p.pyds_sampler.resample()
            else:
                continue

    def _fix_indicator_var(self):
        for scen in self.scenarios():
            if hasattr(scen, 'pyds_indicator_var'):
                var = scen.pyds_indicator_var
                if var.value==0:
                    var.fix(0)
                else:
                    var.fix(1)
            else:
                continue
    
    def _reset_indicator_var(self):
        for scen in self.scenarios():
            if hasattr(scen, 'pyds_indicator_var'):
                var = scen.pyds_indicator_var
                var.unfix()
                var.setvalue(0)
            else:
                continue

    def set_param_value(self, sample_value):
        """Set the values of the stochastic parameters manually e.g. external sampling.

        Args:
            data (dict): parameter name, values
        """
        for k,v in sample_value.items():
            for scen in self.scenarios():
                scen.param[k].value = v

class TwoStageModel(StochModel):
    def __init__(self, options, all_scenario_names, scenario_creator, scenario_creator_kwargs):
        super().__init__(options, all_scenario_names, scenario_creator, scenario_creator_kwargs)

    def simulate_dae():
        return

    




