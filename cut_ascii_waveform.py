#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob


def read_ascii_waveform(file):
    waveform = np.loadtxt(file, dtype="float")
    time = waveform[:, 0]
    amp = waveform[:, 1]
    return time, amp


if __name__ == "__main__":
    datapath = './inversions_noise_model_2/seisflows_fwd_obs/output/solver'
    outdir = './inversions_noise_model_2/fdfwani/data_obs'
    NEW_NT = 3001

    # dont edit below this line
    NEW_BRANCH_NPTS = int((NEW_NT-1)/2)
    events = glob(os.path.join(datapath, '*'))

    for event in events:
        if not os.path.isdir(event):
            continue

        event_outdir = os.path.join(outdir, os.path.basename(event))
        if not os.path.exists(event_outdir):
            os.mkdir(event_outdir)

        files = glob(os.path.join(event, '*'))

        for f in files:
            time, amp = read_ascii_waveform(f)

            OT = time[0]
            DT = abs(time[-1] - time[-2])
            NT = len(time)
            ZLS = int((NT-1)/2)

            amp2 = np.zeros(NEW_NT)
            time2 = np.zeros(NEW_NT)

            amp2 = amp[ZLS-NEW_BRANCH_NPTS:ZLS+NEW_BRANCH_NPTS+1]
            time2 = np.arange(0, NEW_NT) * DT
            time2 = time2 + OT

            dataout = np.vstack((time2, amp2)).T
            outfile = os.path.join(event_outdir, os.path.basename(f))
            np.savetxt(outfile, dataout, ["%17.7f", "%17.7f"])

        print('done ', event)
