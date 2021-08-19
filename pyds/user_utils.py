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
