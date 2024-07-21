#!/usr/bin/env python3
import numpy as np
import os
import pickle
import sys
import yaml
from mpi4py import MPI
from scipy.interpolate import RegularGridInterpolator
from auxiliary_functions import *


def mute_model_on_pml_elements(xcoor, zcoor, m, par):
    # this only works for rectangular meshes
    pml_xmin = par['MESH']['xmin'] + par['PML_THICKNESS']
    pml_xmax = par['MESH']['xmax'] - par['PML_THICKNESS']
    pml_zmin = par['MESH']['zmin'] + par['PML_THICKNESS']
    pml_zmax = par['MESH']['zmax'] - par['PML_THICKNESS']

    xunique = np.unique(xcoor)
    zunique = np.unique(zcoor)

    # top
    idx_pml_top = np.argwhere(zcoor >= pml_zmax).flatten()
    m[idx_pml_top] = 0.0

    # bottom
    idx_pml_bot = np.argwhere(zcoor <= pml_zmin).flatten()
    m[idx_pml_bot] = 0.0

    # right
    idx_pml_right = np.argwhere(xcoor >= pml_xmax).flatten()
    m[idx_pml_right] = 0.0

    # left
    idx_pml_left = np.argwhere(xcoor <= pml_xmin).flatten()
    m[idx_pml_left] = 0.0
    return m

# MPI
comm = MPI.COMM_WORLD
proc = comm.Get_rank()
nproc = comm.Get_size()

# read parameters
parfile = sys.argv[1]

try:
    with open(parfile, 'r') as _file:
        par = yaml.safe_load(_file)
except IOError:
    print('IOError: parfile not found.')

# read gll coordinates
mesh_file = 'mesh_glob{:06}'.format(proc)
mesh_glob = np.loadtxt(os.path.join(par['MESH_DIR'], mesh_file), dtype='float32')
xcoor = mesh_glob[:, 0]
zcoor = mesh_glob[:, 1]

m_file = 'proc{:06}_{}.bin'.format(proc, 'mask_noise')
m = read_bin(os.path.join(par['MODEL_DIR'], m_file))
print(xcoor.size, m.size)

if par['PML_THICKNESS']:
    m = mute_model_on_pml_elements(xcoor, zcoor, m, par)

write_bin(m, os.path.join(par['OUTPUT_DIR'], f'{m_file}'))
print('proc {} done'.format(proc))
