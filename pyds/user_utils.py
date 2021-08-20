import pickle

def read_output_file(filepath):
    objs = []
    with open(filepath, 'rb+') as f:
        while True:
            try:
                objs.append(pickle.load(f))
            except EOFError:
                break
    return objs

def read_last_file(filepath):
    prev_obj = None
    with open(filepath, 'rb+') as f:
        while True:
            try:
                prev_obj = pickle.load(f)
            except EOFError:
                break
    return prev_obj