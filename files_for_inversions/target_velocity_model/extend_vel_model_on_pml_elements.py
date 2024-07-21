#!/usr/bin/env python3
import numpy as np
import os
import pickle
import sys
import yaml
from mpi4py import MPI
from scipy.interpolate import RegularGridInterpolator
from auxiliary_functions import *


def extend_model_on_pml_elements(xcoor, zcoor, m, par):
    # this only works for rectangular meshes
    pml_xmin = par['MESH']['xmin'] + par['PML_THICKNESS']
    pml_xmax = par['MESH']['xmax'] - par['PML_THICKNESS']
    pml_zmin = par['MESH']['zmin'] + par['PML_THICKNESS']
    pml_zmax = par['MESH']['zmax'] - par['PML_THICKNESS']

    xunique = np.unique(xcoor)
    zunique = np.unique(zcoor)

    # top
    idx_pml_top = np.argwhere(zcoor >= pml_zmax).flatten()
    if idx_pml_top.size > 0:
        for x in xunique:
            idx = np.argwhere(xcoor[idx_pml_top] == x)
            if idx.size > 0:
                idx2 = np.argmin(zcoor[idx_pml_top[idx]])
                m[idx_pml_top[idx]] = m[idx_pml_top[idx[idx2]]]

    # bottom
    idx_pml_bot = np.argwhere(zcoor <= pml_zmin).flatten()
    if idx_pml_bot.size > 0:
        for x in xunique:
            idx = np.argwhere(xcoor[idx_pml_bot] == x)
            if idx.size > 0:
                idx2 = np.argmax(zcoor[idx_pml_bot[idx]])
                m[idx_pml_bot[idx]] = m[idx_pml_bot[idx[idx2]]]

    # right
    idx_pml_right = np.argwhere(xcoor >= pml_xmax).flatten()
    if idx_pml_right.size > 0:
        for z in zunique:
            idx = np.argwhere(zcoor[idx_pml_right] == z)
            if idx.size > 0:
                idx2 = np.argmin(xcoor[idx_pml_right[idx]])
                m[idx_pml_right[idx]] = m[idx_pml_right[idx[idx2]]]

    # left
    idx_pml_left = np.argwhere(xcoor <= pml_xmin).flatten()
    if idx_pml_left.size > 0:
        for z in zunique:
            idx = np.argwhere(zcoor[idx_pml_left] == z)
            if idx.size > 0:
                idx2 = np.argmax(xcoor[idx_pml_left[idx]])
                m[idx_pml_left[idx]] = m[idx_pml_left[idx[idx2]]]

    return m


def extend_model_on_pml_elements_circ(xcoor, zcoor, m, par):
    # this only works for circular meshes
    r = np.sqrt(np.square(xcoor-par["X0"]) + np.square(zcoor-par["Z0"]))
    pml_start = par["R"] - par['PML_THICKNESS']

    # calculate azimuth of all gll points
    az = np.arctan2((zcoor-par["Z0"]) , (xcoor-par["X0"])) + np.pi
    az = np.rad2deg(az)

    # loop over azimuth bins
    angles = np.linspace(0, 361, 720)
    for i in range(0, angles.size-1):
        idx_az = np.argwhere((az > angles[i]) & (az < angles[i+1])).flatten()
        if idx_az.size > 0:
            idx_pml = np.argwhere(r[idx_az] > pml_start).flatten()
            if idx_pml.size > 0:
                idx_lim = np.argmin(r[idx_az[idx_pml]]).flatten()
                m[idx_az[idx_pml]] = m[idx_az[idx_pml[idx_lim]]]
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
xcoor_file = 'proc{:06}_x.bin'.format(proc)
xcoor = read_bin(os.path.join(par['DATA_DIR'], xcoor_file), dtype='float32')

zcoor_file = 'proc{:06}_z.bin'.format(proc)
zcoor = read_bin(os.path.join(par['DATA_DIR'], zcoor_file), dtype='float32')

for mpar in par['MODEL_PARAMETERS']:
    m_file = 'proc{:06}_{}.bin'.format(proc, mpar)
    m = read_bin(os.path.join(par['DATA_DIR'], m_file), dtype='float32')

    if par['PML_THICKNESS']:
        if par["CIRCLE"]:
            m = extend_model_on_pml_elements_circ(xcoor, zcoor, m, par)
        else:
            m = extend_model_on_pml_elements(xcoor, zcoor, m, par)

    write_bin(m, os.path.join(par['OUTPUT_DIR'], f'{m_file}'))

print('proc {} done'.format(proc))
