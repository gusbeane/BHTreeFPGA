import numpy as np
from tqdm import tqdm

np.random.seed(9996)

Npart = 10000

# Maximum value for 32-bit unsigned integer
MAX_INT32 = np.iinfo(np.uint32).max

# Generate random 3D positions for particle field
pos = np.random.randint(0, MAX_INT32, size=(Npart, 3), dtype=np.uint32)
pos_float = ((pos.astype(np.float64)) / MAX_INT32).astype(np.float32)
M = np.full(Npart, 1. / Npart)
G = 1.

# Softening length equal to interparticle spacing
eps = np.full(Npart, np.sqrt(3.) / Npart)

# Generate random 3D positions to evaluate acc at
pos_eval = np.random.randint(0, MAX_INT32, size=(Npart, 3), dtype=np.uint32)
pos_eval_float = ((pos_eval.astype(np.float64)) / MAX_INT32).astype(np.float32)

# now we compute the accelerations
acc = np.zeros_like(pos_eval, dtype=np.float32)
for i in tqdm(range(Npart)):
    rsq = np.sum((pos_float - pos_eval_float[i])**2, axis=1)
    rcubed = rsq**(3./2.)

    dacc = np.zeros_like(pos_float)
    np.divide(pos_float - pos_eval_float[i], rcubed[:, None], out=dacc, where=rcubed[:, None]!=0)
    acc[i] = G * np.sum(dacc * M[:, None], axis=0)

# Save all values to a single compressed binary file
np.savez_compressed('nbody_golden.npz',
                   pos=pos,
                   pos_float=pos_float,
                   pos_eval=pos_eval,
                   pos_eval_float=pos_eval_float,
                   acc_eval=acc,
                   M=M,
                   G=G,
                   eps=eps)