import pickle
from os import remove
from os.path import exists, join
import atexit

class OutputManager:
    def __init__(self, folder):
        self.output_filename = 'pyds_output'
        self.output_path_fnx = join(folder, self.output_filename + '.pkl')
        if exists(self.output_path_fnx):
            remove(self.output_path_fnx)
        
        self.clear_buffer()
        
    def clear_buffer(self):
        self.data = {
            'input': None,
            'solver': {},
            'simulator': {}
        }
        

    def add_solver_solution(self, container):
        self.data['solver'].update(container.copy())

    def add_simulator_solution(self, container):
        self.data['simulator'].update(container.copy())

    def add_input(self, data):
        self.data['input'] = data.copy()

    def write_data_to_disk(self):
        with open(self.output_path_fnx, 'ab+') as f:
            pickle.dump(self.data, f)


