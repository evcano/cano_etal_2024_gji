#!/usr/bin/env python3
import numpy as np
import os
import pickle
import sys
import yaml
from matplotlib.path import Path
from mpi4py import MPI
from scipy.interpolate import RegularGridInterpolator


def write_bin(x, filename):
    n = np.array([4 * len(x)], dtype='int32')
    x = np.array(x, dtype='float32')
    with open(filename, 'wb') as file:
        n.tofile(file)
        x.tofile(file)
        n.tofile(file)


def read_bin(filename, dtype='float32'):
    nbytes = os.path.getsize(filename)
    with open(filename, 'rb') as file:
        # read size of record
        file.seek(0)
        n = np.fromfile(file, dtype='int32', count=1)[0]

        if n == nbytes-8:
            file.seek(4)
            data = np.fromfile(file, dtype=dtype)
            return data[:-1]
        else:
            file.seek(0)
            data = np.fromfile(file, dtype=dtype)
            return data


def mask_model(xcoor, zcoor):
    polygon = [(-4.82e5,3.03e5),
               (-4.82e5,1.76e5),
               (-3.667e5,3.1e4),
               (-1.992e5,-1.62e5),
               (6.97e4,-2.44e5),
               (1.773e5,-3.05e5),
               (3.51e5,-3.18e5),
               (4.163e5,-2.97e5),
               (4.163e5,-4.5e4),
               (-4.09e4, 3.03e5),
               (-4.82e5, 3.03e5),
              ]
    poly_path = Path(polygon)
    coors = np.column_stack((xcoor, zcoor))
    mask_idx = np.invert(poly_path.contains_points(coors))
    return mask_idx


# MPI
comm = MPI.COMM_WORLD
proc = comm.Get_rank()
nproc = comm.Get_size()

DATA_DIR = "./inversions_noise_model_1/fdfwani/model_init"
OUTPUT_DIR = "./gradient_mask_truncated"

# read gll coordinates
xcoor_file = 'proc{:06}_x.bin'.format(proc)
zcoor_file = 'proc{:06}_z.bin'.format(proc)

xcoor = read_bin(os.path.join(DATA_DIR, xcoor_file), dtype='float32')
zcoor = read_bin(os.path.join(DATA_DIR, zcoor_file), dtype='float32')

points_outside_array = mask_model(xcoor, zcoor)
mask = np.ones(xcoor.shape)
mask[points_outside_array] = 0.0

m_file = f"proc{proc:06}_gradient_mask.bin"
write_bin(mask, os.path.join(OUTPUT_DIR, f'{m_file}'))
print('proc {} done'.format(proc))
