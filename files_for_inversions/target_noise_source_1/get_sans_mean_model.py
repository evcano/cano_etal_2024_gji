#!/usr/bin/env python3
import cmasher as cmr
import h5py
import matplotlib.pyplot as plt
import numpy as np
import sys

from glob import glob
from noisi.util.plot import plot_grid


def plot_source():
    # load h5py file
    files = glob('./sans_maps_global_2022_01_01_days_341/models/*')
    src = np.zeros(29020)

    c = 0
    for f in files:
        src_file = h5py.File(f, 'r')

        # get noise distribution and coordinates
        tmp_grd = np.array(src_file['coordinates'])
        tmp_src = np.array(src_file['model']).T[0]
        tmp_src /= np.max(tmp_src)

        if tmp_src.size == src.size:
            grd = tmp_grd
            src += tmp_src
            c += 1

        src_file.close()

    src /= c

    plot_grid(grd[0], grd[1], src, sequential=True, cmap=cmr.cosmic, size=100)
    plt.show()
    plt.close()

    np.save('mean_model_grd', grd)
    np.save('mean_model', src)

if __name__ == '__main__':
    plot_source()
