#!/usr/bin/env python3
import numpy as np
import os
import pickle
import sys
import yaml
from mpi4py import MPI
from scipy.interpolate import RegularGridInterpolator
from auxiliary_functions import *


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
xcoor_file = 'proc{:06}_x.bin'.format(proc)
xcoor = read_bin(os.path.join(par['DATA_DIR'], xcoor_file), dtype='float32')

zcoor_file = 'proc{:06}_z.bin'.format(proc)
zcoor = read_bin(os.path.join(par['DATA_DIR'], zcoor_file), dtype='float32')

dm = np.zeros(xcoor.size)

for P in par['PERTURBATIONS']:
    # get points to perturb
    if P['shape'] == 'square':
        gll_points, idx = points_inside_square(xcoor,
                                               zcoor,
                                               [P['xmin'], P['xmax']],
                                               [P['zmin'], P['zmax']])
    elif P['shape'] == 'circle':
        gll_points, idx = points_inside_circle(xcoor,
                                               zcoor,
                                               [P['x0'], P['z0']],
                                               P['r'])

    if not gll_points.any():
        continue

    # define taper
    if P['taper']:
        if P['type'] == 'dont_perturb':
            reverse_taper = True
        else:
            reverse_taper = False

        if P['shape'] == 'square':
            taper = taper_square(gll_points,
                                 [P['xmin'], P['xmax'], P['zmin'], P['zmax']],
                                 P['taper_side'],
                                 P['taper_width'],
                                 reverse_taper)
        elif P['shape'] == 'circle':
            taper = taper_circle(gll_points,
                                 [P['x0'], P['z0']],
                                 P['r'],
                                 P['taper_r'],
                                 reverse_taper)
    else:
        taper = np.ones(gll_points.shape[0])

    # define perturbations
    if P['type'] == 'gaussian':
        with open(P['gfield'], 'rb') as _file:
            interpolator = pickle.load(_file)

        dm[idx] = interpolator(gll_points) * taper

    elif P['type'] == 'uniform':
        dm[idx] = P['magnitude'] * taper

    elif P['type'] == 'dont_perturb':
        taper[taper == 1.0] = 0.0
        dm[idx] *= taper

# read and perturb  model
for mpar in par['MODEL_PARAMETERS']:
    m_file = 'proc{:06}_{}.bin'.format(proc, mpar)
    m = read_bin(os.path.join(par['DATA_DIR'], m_file), dtype='float32')

    m2 = m * (1.0 + dm)
    write_bin(m2, os.path.join(par['OUTPUT_DIR'], f'{m_file}'))

print('proc {} done'.format(proc))
