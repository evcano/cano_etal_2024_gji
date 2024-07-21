#!/usr/bin/env python3

import cmasher as cmr
import cartopy.crs as ccrs
import h5py
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import sys
from noisi.util.plot import plot_grid


def plot_source(src_file):
    # load h5py file
    src_file = h5py.File(src_file, 'r')

    # get noise distribution and coordinates
    src = np.array(src_file['model']).T[0]
    grd = np.array(src_file['coordinates'])

    plot_grid(grd[0], grd[1], src, sequential=True, cmap=cmr.cosmic, size=100)
    plt.show()
    plt.close()

if __name__ == '__main__':
    src_file = sys.argv[1]
    plot_source(src_file)
