import pyomo.environ as pyo
from pyomo import dae
from numpy import random
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import pickle

from deus import DEUS
from pyds.buildutils import *
from pyds.solver import *


def stage0_rule(m, stage=None):
    m.t = dae.ContinuousSet(bounds=(0,1)) #Normalised time
    m.ca_in = pyo.Var(m.t, initialize=0, bounds=(0, 0.001))
    m.T = pyo.Param(initialize=0, mutable=True)
    m.t_f = pyo.Param(initialize=0, mutable=True)

def stage1_rule(m, stage):
    t = stage[0].t
    t_f = stage[0].t_f
    T = stage[0].T
    m.I = pyo.Set(initialize=['A', 'B', 'C'])
    #Parameters
    
    m.E1 = pyo.Param(mutable=True, initialize=0)
    m.E2 = pyo.Param(mutable=True, initialize=0) 
    m.k1 = pyo.Param(mutable=True, initialize=0)
    m.k2 = pyo.Param(mutable=True, initialize=0) 

    #Variables
    m.c = pyo.Var(m.I, t, bounds=(0, 1E8), initialize=0)
    m.dcdt = dae.DerivativeVar(m.c, wrt=t)
    m.u = pyo.Var(t, initialize=0)
    m.dudt = dae.DerivativeVar(m.u, initialize=0, wrt=t)
     #Defined in natural time, not normalised
    m.u[0].fix(0)
    def r_ca(_m, t):
        return stage[0].ca_in[t]/t_f == _m.dudt[t]
    m.c_ca = pyo.Constraint(t, rule=r_ca)

    #Constraints
    R = 1
    def r_mb1(_m, t): 
        return _m.dcdt['A',t]/t_f == -2*_m.k1*pyo.exp(-_m.E1/(R*T))*_m.c['A',t]**2 + stage[0].ca_in[t]
    def r_mb2(_m, t):
        return _m.dcdt['B',t]/t_f == _m.k1*pyo.exp(-_m.E1/(R*T))*_m.c['A',t]**2 \
            - _m.k2*pyo.exp(-_m.E2/(R*T))*_m.c['B',t]
    def r_mb3(_m, t):
        return  _m.dcdt['C',t]/t_f  == _m.k2*pyo.exp(-_m.E2/(R*T))*_m.c['B',t]
    m.mb1 = pyo.Constraint(t, rule=r_mb1)
    m.mb2 = pyo.Constraint(t, rule=r_mb2)
    m.mb3 = pyo.Constraint(t, rule=r_mb3)

    m.g1 = pyo.Var(initialize=0)
    m.g2 = pyo.Var(initialize=0)
    def r_cqa1(_m):
        return _m.g1 == 1000*(0.8 -(_m.c['B',1] + 1E-2)/(_m.c['A', 1]+_m.c['B',1]+_m.c['C',1]+ 1E-2))
    def r_cqa2(_m):
        return _m.g2 == -(100*_m.c['B',1] - 20*_m.c['A',0] - 128*(t_f + 30))/1E3
    m.c_g1 = pyo.Constraint(rule=r_cqa1)
    m.c_g2 = pyo.Constraint(rule=r_cqa2)

    m.bigM_constant = pyo.Param([1,2], initialize={1: 1E3, 2: 1E5})
    add_BigMConstraint(m, 'cqa1', m.bigM_constant[1], expr=m.g1<=0)  
    add_BigMConstraint(m, 'cqa2', m.bigM_constant[2], expr=m.g2<=0)
    m.c['A',0].fix(2E3) #mol/L
    m.c['B',0].fix(0)
    m.c['C',0].fix(0)

def apply_collocation(m):
    discretizer = pyo.TransformationFactory('dae.collocation')
    discretizer.apply_to(m, nfe=8, ncp=3)

m = TwoStageManager(stage_rules=[stage0_rule, stage1_rule], 
    input={'d': [(0, 't_f'), (0, 'T')], 'p': [(1, 'E1'), (1, 'E2'), (1, 'k1'), (1, 'k2')]},
    solver_name='gams',
    model_transformation=apply_collocation,
    solver_options={'solver': 'conopt'})

ran = random.default_rng(12341234)
n_samples = 10

E1 = ran.normal(2.5E3, 2.5E2, size=n_samples)
E2 = ran.normal(5E3, 5E2, size=n_samples)
k1 = ran.normal(6.409E-2, 6.409E-3, size=n_samples)
k2 = ran.normal(9.938E3, 9.938E2, size=n_samples)

p_best = [2.5E3, 5E3, 6.409E-2, 9.938E3]
p_samples = [{'c': [E1[i], E2[i], k1[i], k2[i]], \
    'w': 1.0/n_samples} for i in range(n_samples)]

the_activity_form = {
    "activity_type": "ds",

    "activity_settings": {
        "case_name": "DS_Test1",
        "case_path": os.getcwd(),
        "resume": False,
        "save_period": 10
    },

    "problem": {
        "user_script_filename": "ds_test2_user_script",
        "constraints_func_name": "g_func_cat2",
        "parameters_best_estimate": p_best,
        "parameters_samples": p_samples,
        "target_reliability": 0.95,
        "design_variables": [
            {"t_f": [250, 350]},
            {"T": [250, 300]}
        ]
    },

    "solver": {
        "name": "ds-ns",
        "settings": {
            "score_evaluation": {
                "method": "serial",
                "constraints_func_ptr": m.g,
                # "constraints_func_ptr": None,
                "store_constraints": False
            },
            # "score_evaluation": {
            #     "method": "mppool",
            #     "pool_size": -1,
            #     "constraints_func_ptr": None,
            #     "store_constraints": False
            # },
            "efp_evaluation": {
                "method": "serial",
                "constraints_func_ptr": m.g,
                # "constraints_func_ptr": None,
                "store_constraints": False,
                "acceleration": False
            },
            # "efp_evaluation": {
            #     "method": "mppool",
            #     "pool_size": -1,
            #     "constraints_func_ptr": None,
            #     "store_constraints": False,
            #     "acceleration": False
            # },
            "points_schedule": [
                (.0, 30, 30),
                (.01, 100, 40),
                (.1, 150, 40),
                (.5, 200, 40),
                (.7, 250, 40),
                (.8, 300, 40)
            ],
            "stop_criteria": [
                {"inside_fraction": 1.0}
            ]
        },
        "algorithms": {
            "sampling": {
                "algorithm": "mc_sampling-ns_global",
                "settings": {
                     "nlive": 10,  # This is overridden by points_schedule
                     "nreplacements": 5,  # This is overriden by points_schedule
                     "prng_seed": 1989,
                     "f0": 0.1,
                     "alpha": 0.3,
                     "stop_criteria": [
                         {"max_iterations": 100000}
                     ],
                     "debug_level": 0,
                     "monitor_performance": False
                 },
                "algorithms": {
                    "replacement": {
                        "sampling": {
                            "algorithm": "suob-ellipsoid"
                        }
                    }
                }
            }
        }
    }
}

the_deus = DEUS(the_activity_form)
t0 = time.time()
the_deus.solve()
print(time.time()-t0)

cs_path = the_activity_form["activity_settings"]["case_path"]
cs_name = the_activity_form["activity_settings"]["case_name"]

with open(cs_path + '/' + cs_name + '/' + 'output.pkl', 'rb') \
        as file:
    output = pickle.load(file)

samples = output["solution"]["probabilistic_phase"]["samples"]
inside_samples_coords = np.empty((0, 2))
outside_samples_coords = np.empty((0, 2))
for i, phi in enumerate(samples["phi"]):
    if phi >= 0.95:
        inside_samples_coords = np.append(inside_samples_coords,
                                          [samples["coordinates"][i]], axis=0)
    else:
        outside_samples_coords = np.append(outside_samples_coords,
                                           [samples["coordinates"][i]], axis=0)

fig1 = plt.figure()
x = inside_samples_coords[:, 0]
y = inside_samples_coords[:, 1]
plt.scatter(x, y, s=10, c='r', alpha=1.0, label='inside')

x = outside_samples_coords[:, 0]
y = outside_samples_coords[:, 1]
plt.scatter(x, y, s=10, c='b', alpha=0.5, label='outside')
plt.legend()
plt.show()

