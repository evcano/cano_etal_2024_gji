#!/usr/bin/env python3
# FOR SPECFEM2D !!!
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import yaml
from glob import glob
from scipy.interpolate import griddata, interp1d


def _read(filename):
    """
    Reads Fortran style binary data into numpy array
    """
    nbytes = os.path.getsize(filename)
    with open(filename, 'rb') as file:
        # read size of record
        file.seek(0)
        n = np.fromfile(file, dtype='int32', count=1)[0]

        if n == nbytes-8:
            file.seek(4)
            data = np.fromfile(file, dtype='float32')
            return data[:-1]
        else:
            file.seek(0)
            data = np.fromfile(file, dtype='float32')
            return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mesh_dir', type=str)
    parser.add_argument('model_dir', type=str)
    parser.add_argument('x0', type=float)
    parser.add_argument('z0', type=float)
    args = parser.parse_args()

    ngll_r = 60 * 5
    ngll_boundary = 180 * 5
    nproc = len(glob(os.path.join(args.model_dir, f"*mask_noise.bin")))

    xcoor = np.array([])
    zcoor = np.array([])
    model = np.array([])

    for i in range(0, nproc):
        mesh_glob = np.loadtxt(os.path.join(args.mesh_dir,
                                            f"mesh_glob{i:06}"))
        xcoor_proc = mesh_glob[:, 0]
        zcoor_proc = mesh_glob[:, 1]

        model_proc = _read(os.path.join(args.model_dir,
                                        f"proc{i:06}_mask_noise.bin"))

        xcoor = np.append(xcoor, xcoor_proc)
        zcoor = np.append(zcoor, zcoor_proc)
        model = np.append(model, model_proc)

    ngll = model.size

    # compute radial and angular coordinates wrt the integral reference point
    all_gll_r = np.sqrt((xcoor - args.x0)**2 + (zcoor - args.z0)**2)
    all_gll_theta = (np.arctan2(zcoor - args.z0, xcoor - args.z0) + 2*np.pi) % (2*np.pi)
    rmax = np.max(all_gll_r)

    # obtain noise ring by integrating an existing noise mask
    # generate a polar grid centered at the integral reference point
    dr = rmax / ngll_r
    dtheta = (2*np.pi) / ngll_boundary

    ax_r = np.arange(0.0, rmax + dr, dr)
    ax_theta = np.arange(0.0, 2*np.pi + dtheta, dtheta)

    R, T = np.meshgrid(ax_r, ax_theta)
    X = R * np.cos(T) + args.x0
    Z = R * np.sin(T) + args.z0 

    # interpolate noise mask into the polar grid
    points = np.zeros((ngll, 2))
    points[:, 0] = np.squeeze(xcoor)
    points[:, 1] = np.squeeze(zcoor)

    mask_noise_intp = griddata(points, model, (X, Z),
                               method='linear', fill_value=np.nan)

    # azimuthal integration to obtain the noise-source ring
    r = np.zeros(ax_r.size)
    r[1:] = 0.5 * (ax_r[0:-1] + ax_r[1:])
    da = r * dr * dtheta

    # azimuthal noise intensity from 0 to 360
    noise_az = np.zeros(ax_theta.size)
    for i, theta in enumerate(ax_theta):
        idx = np.where(T == theta)
        noise_az[i] = np.nansum(mask_noise_intp[idx] * da)

    # the ax_theta is in radians measured counterclock wise from positive xaxis
    # convert to degrees measured clockwise from positive yaxis
    az = np.rad2deg(ax_theta)
    az = (450.0 - az) % 360.0

    np.savetxt(os.path.join(args.model_dir, 'azimuthal_noise.ascii'),
               np.array([az, noise_az]).T)
