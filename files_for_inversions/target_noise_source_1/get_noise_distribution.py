#!/usr/bin/env python3
import numpy as np
import os
import pickle
import sys
import yaml
from mpi4py import MPI
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator
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
mesh_glob = np.loadtxt(os.path.join(par['MESH_DIR'], f'mesh_glob{proc:06}'))
xcoor = mesh_glob[:, 0]
zcoor = mesh_glob[:, 1]

mask_noise = np.zeros(xcoor.size, dtype='float32')

for D in par['DISTRIBUTIONS']:
    # get points to perturb
    if D['shape'] == 'square':
        gll_points, idx = points_inside_square(xcoor,
                                               zcoor,
                                               [D['xmin'], D['xmax']],
                                               [D['zmin'], D['zmax']])
    elif D['shape'] == 'circle':
        gll_points, idx = points_inside_circle(xcoor,
                                               zcoor,
                                               [D['x0'], D['z0']],
                                               D['r'])
    elif D['shape'] == 'ellipse':
        gll_points, idx = points_inside_ellipse(xcoor,
                                                zcoor,
                                                [D['x0'], D['z0']],
                                                [D['rx'], D['rz']],
                                                D['theta'])

    if not gll_points.any():
        continue

    # define taper
    if D['taper']:
        if D['type'] == 'dont_perturb':
            reverse_taper = True
        else:
            reverse_taper = False

        if D['shape'] == 'square':
            taper = taper_square(gll_points,
                                 [D['xmin'], D['xmax'], D['zmin'], D['zmax']],
                                 D['taper_side'],
                                 D['taper_width'],
                                 reverse_taper)
        elif D['shape'] == 'circle':
            taper = taper_circle(gll_points,
                                 [D['x0'], D['z0']],
                                 D['r'],
                                 D['taper_r'],
                                 reverse_taper)
        elif D['shape'] == 'ellipse':
            taper = taper_ellipse(gll_points,
                                  [D['x0'], D['z0']],
                                  [D['rx'], D['rz']],
                                  D['theta'],
                                  D['taper_phi'],
                                  reverse_taper)
    else:
        taper = np.ones(gll_points.shape[0])

    # define perturbations
    if D['type'] == 'gaussian':
        with open(D['gfield'], 'rb') as _file:
            interpolator = pickle.load(_file)

        mask_noise[idx] += interpolator(gll_points) * taper

    elif D['type'] == 'uniform':
        mask_noise[idx] += D['magnitude'] * taper

    elif D['type'] == 'dont_perturb':
        taper[taper == 1.0] = 0.0
        mask_noise[idx] *= taper

print(mask_noise.min(), mask_noise.max())
write_bin(mask_noise, os.path.join(par['OUTPUT_DIR'],
                                   f'proc{proc:06}_mask_noise.bin'))

print('proc {} done'.format(proc))
