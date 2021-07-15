from pyomo.core.base import Param

def resample_all_stoch_params(model):
    all_params = model.component_objects(Param, descend_into=True)
    for p in all_params:
        if hasattr(p, '_pyds_sampler'):
            p._pyds_sampler.resample_param()
        else:
            continue

def fix_binary_var()