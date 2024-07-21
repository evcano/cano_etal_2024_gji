#!/usr/bin/env python3
import numpy as np
import os
from scipy.interpolate import interp1d
from scipy.signal.windows import barthann


def points_inside_circle(xcoor, zcoor, center, radius):
    x0 = center[0]
    z0 = center[1]

    rcoor = np.sqrt((xcoor - x0) ** 2 + (zcoor - z0) ** 2)
    idx = np.where(rcoor <= radius)[0]

    if not idx.any():
        return np.array([]), []
    else:
        points = np.zeros((len(idx), 2))

        for i, j in enumerate(idx):
            points[i, 0] = xcoor[j]
            points[i, 1] = zcoor[j]

    return points, idx


def points_inside_ellipse(xcoor, zcoor, center, ax_radius, theta):
    x0 = center[0]
    z0 = center[1]

    rad_x = ax_radius[0]
    rad_z = ax_radius[1]

    theta = np.deg2rad(theta)

    phi = (((xcoor-x0)*np.cos(theta)+(zcoor-z0)*np.sin(theta))/rad_x) ** 2 + \
          (((xcoor-x0)*np.sin(theta)-(zcoor-z0)*np.cos(theta))/rad_z) ** 2

    idx = np.where(phi < 1.0)[0]

    if not idx.any():
        return np.array([]), []
    else:
        points = np.zeros((len(idx), 2))

        for i, j in enumerate(idx):
            points[i, 0] = xcoor[j]
            points[i, 1] = zcoor[j]

    return points, idx

def points_inside_square(xcoor, zcoor, xlim, zlim):
    xmin = xlim[0]
    xmax = xlim[1]
    zmin = zlim[0]
    zmax = zlim[1]

    points = []
    idx = []

    for i in range(0, xcoor.size):
        x = xcoor[i]
        z = zcoor[i]

        if x >= xmin and x <= xmax and z >= zmin and z <= zmax:
            points.append((x, z))
            idx.append(i)

    return np.array(points), idx


def taper_circle(points, circle_center, circle_radius, taper_start_radius,
                 reverse_taper=False):

    x0 = circle_center[0]
    z0 = circle_center[1]

    taper = np.ones(points.shape[0])

    b = barthann(200)
    b = b[100:]

    if reverse_taper:
        b = b[::-1]

    rcoor = np.sqrt((points[:, 0] - x0) ** 2 + (points[:, 1] - z0) ** 2)
    idx = np.argwhere(rcoor >= taper_start_radius)

    if idx.any():
        a = np.linspace(taper_start_radius, circle_radius, 100)
        inter = interp1d(a, b, bounds_error=False, fill_value=0.0)
        taper[idx] = inter(rcoor[idx])

    return taper


def taper_ellipse(points, center, ax_radius, theta, taper_phi,
                  reverse_taper=False):

    xcoor = points[:, 0]
    zcoor = points[:, 1]

    x0 = center[0]
    z0 = center[1]

    rad_x = ax_radius[0]
    rad_z = ax_radius[1]

    theta = np.deg2rad(theta)

    taper = np.ones(points.shape[0])

    b = barthann(200)
    b = b[100:]

    if reverse_taper:
        b = b[::-1]

    phi = (((xcoor-x0)*np.cos(theta)+(zcoor-z0)*np.sin(theta))/rad_x) ** 2 + \
          (((xcoor-x0)*np.sin(theta)-(zcoor-z0)*np.cos(theta))/rad_z) ** 2

    idx = np.where(phi >= taper_phi)[0]

    if idx.any():
        a = np.linspace(taper_phi, 1.0, 100)
        inter = interp1d(a, b, bounds_error=False, fill_value=0.0)
        taper[idx] = inter(phi[idx])

    return taper


def taper_square(points, square_limits, taper_side, taper_width,
                 reverse_taper):

    xmin = square_limits[0]
    xmax = square_limits[1]
    zmin = square_limits[2]
    zmax = square_limits[3]

    taper = np.ones(points.shape[0])

    b1 = barthann(200)
    b1 = b1[100:]

    if reverse_taper:
        b1 = b1[::-1]

    b2 = b1[::-1]

    if taper_side == 'left' or taper_side == 'bottom':
        if taper_side == 'left':
            tmp1 = xmin
            tmp2 = 0
        elif taper_side == 'bottom':
            tmp1 = zmin
            tmp2 = 1

        taper_start = tmp1 + taper_width
        idx = np.argwhere(points[:, tmp2] <= taper_start)

        if idx.any():
            a = np.linspace(tmp1, taper_start, 100)
            ifunc = interp1d(a, b2, bounds_error=False, fill_value=0.0)
            taper[idx] = ifunc(points[idx, tmp2])

    if taper_side == 'right' or taper_side == 'top':
        if taper_side == 'right':
            tmp1 = xmax
            tmp2 = 0
        elif taper_side == 'top':
            tmp1 = zmax
            tmp2 = 1

        taper_start = tmp1 - taper_width
        idx = np.argwhere(points[:, tmp2] >= taper_start)

        if idx.any():
            a = np.linspace(taper_start, tmp1, 100)
            ifunc = interp1d(a, b1, bounds_error=False, fill_value=0.0)
            taper[idx] = ifunc(points[idx, tmp2])

    return taper


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
