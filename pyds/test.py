import numpy as np

def stages_to_update(idx):
        a = np.array(idx)
        should_update_mask = [np.all(a[i+1:]==0) for i in range(len(a))] #empty slice of last i evaluates to True
        stages = [i for i in range(len(a)) if should_update_mask[i]]
        return [tuple(idx[:i+1]) for i in stages]
        

a = [1,2,0,1,0,0]
print(stages_to_update(a))
