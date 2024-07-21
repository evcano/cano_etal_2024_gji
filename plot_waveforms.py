#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob
from scipy.signal import correlation_lags
from obspy.core import Trace


def read_ascii_waveform(file):
    tmp = np.loadtxt(file, dtype="float")
    t = tmp[:, 0]
    a = tmp[:, 1]
    dt = t[-1] - t[-2]
    nt = t.size
    ds = 1./dt
    h = {"npts":nt, "delta":dt, "sampling_rate":ds, "channel":"HXZ"}
    tr = Trace(data=a, header=h)
    return tr


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_path', type=str)
    parser.add_argument('--fil', type=float, default=None, nargs=2)
    parser.add_argument('--nor', action='store_true')
    parser.add_argument('--lim', type=float, default=None, nargs=2)
    parser.add_argument('--lag', action='store_true')
    args = parser.parse_args()

    files = glob(args.input_path)
    files.sort()

    offset = 0.

    for file_ in files:
        tr = read_ascii_waveform(file_)

        if args.fil:
            tr.filter("bandpass", freqmin=args.fil[0], freqmax=args.fil[1],
                      corners=2, zerophase=True)

        if args.nor:
            tr.normalize()
        tr.data += offset

        if args.lag:
            hnpts = (tr.stats.npts + 1.) / 2.
            times = correlation_lags(hnpts, hnpts, mode="full") * tr.stats.delta
        else:
            times = tr.times()

        plt.plot(times, tr.data, color='black', alpha=0.8)
        offset = np.max(tr.data) + 0.5

    if args.lim:
        plt.xlim(args.lim[0], args.lim[1])
    plt.xlabel('time (s)')
    plt.show()
