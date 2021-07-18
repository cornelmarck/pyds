from pyomo.core.base import Param
from mpisppy.opt.ef import ExtensiveForm

class StochModel(ExtensiveForm):
    """Create the multi-stage stochastic model and associated solver functions.

    Args:
        ExtensiveForm ([type]): [description]
    """
    def __init__(self, options, all_scenario_names, scenario_creator, scenario_creator_kwargs):
        options.setdefault('solver':None)
        super().__init__(options, all_scenario_names, scenario_creator, scenario_creator_kwargs)

        return

    def solve_extensive_form(self, solver_options, tee):
        return super().solve_extensive_form(solver_options=solver_options, tee=tee)

    def resample_all_stoch_params(self):
        all_params = self.ef.component_objects(Param, descend_into=True)
        for p in all_params:
            if hasattr(p, '_pyds_sampler'):
                p.pyds_sampler.resample()
            else:
                continue

    def fix_binary_var(self):
        return


