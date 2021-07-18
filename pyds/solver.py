from pyomo.core.base import Param
from mpisppy.opt.ef import ExtensiveForm

class StochModel(ExtensiveForm):
    """Create the multi-stage stochastic model and associated solver functions.

    Args:
        ExtensiveForm ([type]): [description]
    """
    def __init__(self, options, all_scenario_names, scenario_creator, scenario_creator_kwargs):
        super().__init__(options, all_scenario_names, scenario_creator, scenario_creator_kwargs)

        return

    def solve(self, resample=True, solver_options=None, tee=False):
        self.resample_all_stoch_params()
        self.reset_indicator_var()
        self.relaxed_result = self.solve_extensive_form(solver_options, tee)
        self.fix_indicator_var()
        return self.solve_extensive_form(solver_options, tee)


    def resample_all_stoch_params(self):
        all_params = self.ef.component_objects(Param, descend_into=True)
        for p in all_params:
            if hasattr(p, 'pyds_sampler'):
                p.pyds_sampler.resample()
            else:
                continue

    def fix_indicator_var(self):
        for scen in self.scenarios():
            if hasattr(scen, 'pyds_indicator_var'):
                var = scen.pyds_indicator_var
                if var.value==0:
                    var.fix(0)
                else:
                    var.fix(1)
            else:
                continue
    
    def reset_indicator_var(self):
        for scen in self.scenarios():
            if hasattr(scen, 'pyds_indicator_var'):
                var = scen.pyds_indicator_var
                var.unfix()
                var.setvalue(0)
            else:
                continue
            




