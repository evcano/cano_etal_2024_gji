import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
from auxiliary_functions import *


coor_dir1 = '../specfem2d_complete_mesh/DATA'
model_dir1 = './complete_mesh/perturbed'

coor_dir2 = '../specfem2d_truncated_mesh/DATA'
outdir = './truncated_mesh'

nproc = 8

# DO NOT EDIT BELOW THIS LINE
# ===========================

# complete model
xcoor1 = np.array([])
zcoor1 = np.array([])
model1 = np.array([])

for i in range(0, nproc):
    xcoor_proc = read_bin(os.path.join(coor_dir1, f"proc{i:06}_x.bin"))
    zcoor_proc = read_bin(os.path.join(coor_dir1, f"proc{i:06}_z.bin"))
    model_proc = read_bin(os.path.join(model_dir1, f'proc{i:06}_vs.bin'))

    xcoor1 = np.append(xcoor1, xcoor_proc)
    zcoor1 = np.append(zcoor1, zcoor_proc)
    model1 = np.append(model1, model_proc)

# truncated model
xcoor2 = np.array([])
zcoor2 = np.array([])
model2 = np.array([])

proc_ngll = []
for i in range(0, nproc):
    xcoor_proc = read_bin(os.path.join(coor_dir2, f"proc{i:06}_x.bin"))
    zcoor_proc = read_bin(os.path.join(coor_dir2, f"proc{i:06}_z.bin"))

    xcoor2 = np.append(xcoor2, xcoor_proc)
    zcoor2 = np.append(zcoor2, zcoor_proc)

    proc_ngll.append(zcoor_proc.size)

model2 = np.zeros(zcoor2.shape)

for i in range(0, model2.size):
    idx = np.argwhere((xcoor1 == xcoor2[i]) & (zcoor1 == zcoor2[i]))
    a = np.array([])
    for j in idx:
        a = np.append(a, model1[j])
    if np.mean(a) == a[0]:
        model2[i] = model1[j]
    else:
        raise Exception

idx1 = 0
idx2 = 0

for i in range(0, nproc):
    idx2 += proc_ngll[i]
    model_proc = model2[idx1:idx2]

    mfile = 'proc{:06}_vs.bin'.format(i)
    write_bin(model_proc, os.path.join(outdir, mfile))

    idx1 = idx2
