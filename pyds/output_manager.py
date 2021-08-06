import pickle
from os import remove
from os.path import exists
import atexit

class OutputManager:
    def __init__(self, folder):
        self.output_filename = 'output'
        self.output_path_fnx = folder + self.output_filename + '.pkl'
        if exists(self.output_path_fnx):
            remove(self.output_path_fnx)
        
        self.output = {
            'runs': {}
        }

    def add_run_to_output(self, data):
        self.output['runs'].extend(data)

    def write_output_to_file(self):
        with open(self.output_path_fnx, 'wb') as f:
            pickle.dump(self.output)