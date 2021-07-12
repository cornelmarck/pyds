import pyomo
from pyomo.core.base import (Var, Param) 

class ModelSolver:
    def __init__(self, model):
        self._model = model
        self.io = OutputWriter()
        return
    
    @property
    def model(self):
        return self._model

    def solve_model(self, input_param=None, solver_opts=None):
        self.output = {}
        self.set_input_param(input_param)
        self.resample_stoch_param()
        self._reset_bigM_indicator_var()

        relaxed_result = self._solve(tee=True, **solver_opts)
        self._model.bigM.indicator_var.pprint()            
        self._fix_bigM_indicator_var()
        #self._model.bigM.indicator_var.pprint() 
        
        result2 = self._solve(**solver_opts)
        self.get_output()
        self.io.write(self.output)
        return

    def get_output(self):
        self.output['obj'] = self._model.obj()
        for v in self._model.component_objects(pyo.Var):
            self.output[v.name] = v.get_values()

        simple_par = (p for p in self._model.component_objects() if isinstance(p, pyomo.core.base.param.SimpleParam)
        for p in simple_par:
            self.output[p.name] = p.value
        
        indexed_par = (p for p in self._model.component_objects() if isinstance(p, pyomo.core.base.param.IndexedParam)
        for p in indexed_par:
            for param_data in p.values():
              self.output[param_data.name] = param_data.value

    def resample_stoch_param(self):
        for _i in m.scenarios:
            self._resample_scenario(_i))
        return

    def set_input_param(self, value):
        m = self._model
        if self._model.input_param_name is None:
            return
        else:
            par = m.component(self._model.input_param_name)
        
        if not par.is_indexed():
            par.set_value(value)
        else:
            for key, val in par.items():
                par[key].set_value(value[key])
        return

    def initialize_manual_solver(self, solver):
        self._solver = solver
        return

    def _solve(self, *args, **kwargs):
        return self._solver.solve(self._model, *args, **kwargs)
        
    def initialize_neos_solver(self, email_address):
        os.environ['NEOS_EMAIL'] = email_address
        self._solver = pyo.SolverManagerFactory('neos')
        return

    def set_bigM_constant(self, bigM_constant):
        m = self._model
        m.bigM.value.set_value=bigM_constant
        return

    def _reset_bigM_indicator_var(self):
        m = self._model
        m.bigM.indicator_var[:].set_value(0)
        m.bigM.indicator_var[:].unfix()
        m.bigM.indicator_var[:].setub(1)
        return

    def _fix_bigM_indicator_var(self):
        m = self._model
        m.bigM.indicator_var[:].setub(1)
        for _i in m.scenarios:
            var = m.bigM.indicator_var[_i]
            if var.value == 0:
                var.fix(0)
            else: 
                var.fix(1)
        return

    def _resample_scenario(self, scenario_idx):
        s = self._model.scenario[scenario_idx]
        for p in s.component_objects(Param):

        if not param.is_indexed():
            param.set_value(m.stoch_sampling_func())
        else:
            values = m.
            for key, par in param.items():
                par.set_value(m.stoch_sampling_func[key]())
        return


class OutputWriter:
    def __init__(self, folder_name=None):
        self._initialize_dir(folder_name)
        self.file_count = 0

    def _initialize_dir(self, folder_name):
        if folder_name is not None:
            self._folder_name = folder_name
        else:
            now = datetime.datetime.now()
            self._folder_name = now.strftime('%Y%m%d %H%M%S')
        
        if not os.path.exists(self._folder_name):
            os.mkdir(self._folder_name)
        self._path = os.path.join(os.getcwd(), self._folder_name)
        return

    def write(self, output):
        filename = str(self.file_count)
        file_path = os.path.join(self._path, '{}.{}'.format(filename, 'p'))
        with open(file_path, 'wb') as output_file:
            try:
                dill.dump(output, output_file)
            except os.error:
                print('IO error! file not saved.')
        self.file_count+=1