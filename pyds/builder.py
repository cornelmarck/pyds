from pyomo.core.base import (Constraint, ConstraintList)

def get_bigM_form(constraint, indicator_var, bigM_constant):
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