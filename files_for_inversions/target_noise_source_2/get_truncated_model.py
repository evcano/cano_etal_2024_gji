import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
from auxiliary_functions import *


coor_dir1 = '../specfem2d_complete_mesh/OUTPUT_FILES'
model_dir1 = './complete_mesh'

coor_dir2 = '../specfem2d_truncated_mesh/OUTPUT_FILES'
outdir = './truncated_mesh'

model_par = 'mask_noise'
nproc = 8

# DO NOT EDIT BELOW THIS LINE
# ===========================

# complete model
xcoor1 = np.array([])
zcoor1 = np.array([])
model1 = np.array([])

for i in range(0, nproc):
    mesh_glob = np.loadtxt(os.path.join(coor_dir1, f'mesh_glob{i:06}'))

    xcoor_proc = mesh_glob[:, 0]
    zcoor_proc = mesh_glob[:, 1]
    model_proc = read_bin(os.path.join(model_dir1, f'proc{i:06}_mask_noise.bin'))

    xcoor1 = np.append(xcoor1, xcoor_proc)
    zcoor1 = np.append(zcoor1, zcoor_proc)
    model1 = np.append(model1, model_proc)

# truncated model
xcoor2 = np.array([])
zcoor2 = np.array([])
model2 = np.array([])

proc_ngll = []
for i in range(0, nproc):
    mesh_glob = np.loadtxt(os.path.join(coor_dir2, f'mesh_glob{i:06}'))

    xcoor_proc = mesh_glob[:, 0]
    zcoor_proc = mesh_glob[:, 1]
    proc_ngll.append(zcoor_proc.size)

    xcoor2 = np.append(xcoor2, xcoor_proc)
    zcoor2 = np.append(zcoor2, zcoor_proc)

model2 = np.zeros(xcoor2.shape)

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
    xcoor_proc = xcoor2[idx1:idx2]
    zcoor_proc = zcoor2[idx1:idx2]
    model_proc = model2[idx1:idx2]
    idx1 = idx2
    mfile = 'proc{:06}_{}.bin'.format(i, model_par)
    write_bin(model_proc, os.path.join(outdir, mfile))
