from pyomo.core.base import (Constraint, ConstraintList)

def make_bigM_form(constraint, indicator_var, bigM_constant):
    """Obtain the big-M form of an (in)equality constraint.

    Args:
        constraint (pyomo.core.base.Constraint): The (in)equality constraint to be deactivated. 
            Constraint list input is not implemented yet. 
        indicator_var (pyomo.core.base.Var): Binary variable (also relaxed) which 0 if all 
            constraints are satisfied.
        bigM_constant (pyomo.core.base.Param): Big-M constant which is scaled to the magnitude 
            of the inequality constraint.
    Returns:
        :pyomo.core.base.Constraint | pyomo.core.base.ConstraintList | None : 
            The constraint in big-m form.
    """
    #Inequality constraints
    if not constraint.is_constructed():
        constraint.construct()

    if (constraint.lower is None) & (constraint.upper is not None):
        return Constraint(
            expr= constraint.body <= constraint.upper + bigM_constant*indicator_var)
    elif (constraint.lower is not None) & (constraint.upper is None):
        return Constraint(
            expr= constraint.body >= constraint.lower - bigM_constant*indicator_var)

    #Equality constraint
    elif (constraint.lower is not None) & (constraint.upper is not None):
        lst = ConstraintList()
        lst.construct()
        lst.add(expr= constraint.body <= constraint.upper + bigM_constant*indicator_var)
        lst.add(expr= constraint.body >= constraint.lower - bigM_constant*indicator_var)
        return lst
    else:
        return None

def set_stochastic(param, func, *args, **kwargs):
    if not param.mutable:
        raise ValueError('Parameter is not mutable')
    param._pyds_sampler = Sampler(param, func, *args, **kwargs)

class Sampler:
    def __init__(self, param, func, *args, **kwargs):
        """Store a resampling function and its arguments.
        TODO: Input error checking

        Args:
            param (pyo.Param): [description]
            func (function)): Random samples generator function.
                Must return a a single float or dict {key, float}
            args, kwargs: Additional arguments passed to func, e.g. a stream of random numbers
        """
        self.func = func
        self.param = param
        self.args = args
        self.kwargs = kwargs
        
    def resample_param(self):
        """Resample and update the values of the stochastic parameter"""   
        sample_value = self.__call__()
        if isinstance(sample_value, dict):
            for k,v in sample_value.items():
                self.param[k].value = v
        else:
            self.param.value = sample_value

    def __call__(self):
        return self.func(*self.args, **self.kwargs)

def is_stochastic(obj):
    return hasattr(obj, "_pyds_stoch_sampler")
