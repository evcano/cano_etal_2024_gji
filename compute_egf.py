#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob
from scipy.fft import next_fast_len
from scipy.signal import tukey
from scipy.ndimage import gaussian_filter1d

def read_ascii_waveform(file):
    waveform = np.loadtxt(file, dtype="float")
    time = waveform[:, 0]
    amp = waveform[:, 1]
    return time, amp


if __name__ == "__main__":
    datapath = './inversions_noise_model_2/seisflows_fwd_obs/output/solver'
    outdir = './inversions_noise_model_2/anat/data_obs'
    
    new_nt = 1501
    water_level = 0.05
    white = True

    # load the noise source power spectral density
    noise_fft = np.loadtxt('./specfem2d_complete_mesh/noisespec.txt')
    noise_fft = noise_fft[:,1].real
    nfft = noise_fft.size

    # apply waterlevel
    wl = np.max(noise_fft) * water_level
    noise_fft[np.argwhere(noise_fft < wl)] = wl

    # load synthetics stf
    time_new, syn_stf = read_ascii_waveform(
        "./specfem2d_egf/OUTPUT_FILES/plot_source_time_function.txt")
    syn_stf *= -1.0  # scaling factor as in the parfile

    syn_stf -= np.mean(syn_stf)
    syn_stf *= tukey(syn_stf.size, 0.1)
    syn_stf_fft = np.fft.fft(syn_stf)

    events = glob(os.path.join(datapath, '*'))
    events.sort()

    for event in events:
        if not os.path.isdir(event):
            continue

        event_outdir = os.path.join(outdir, os.path.basename(event))
        if not os.path.exists(event_outdir):
            os.mkdir(event_outdir)

        files = glob(os.path.join(event, '*'))

        for f in files:
            time, amp = read_ascii_waveform(f)

            DT = abs(time[-1] - time[-2])
            NT = len(time)
            ZLS = int((NT-1)/2)

            if white:
                # get spectrum
                amp -= np.mean(amp)
                amp *= tukey(NT, 0.1)
                fft = np.fft.fft(amp, nfft)
    
                # spectral whitening
                fft_white = np.divide(fft, noise_fft)
    
                # back to time domain
                amp2 = np.real(np.fft.ifft(fft_white))
            else:
                amp2 = amp

            # average both branches
            neg = amp2[0:ZLS+1]
            neg = neg[::-1]
            pos = amp2[ZLS:]
            sym = (pos + neg) / 2.0

            # compute egf as negative time derivative 
            egf = np.diff(sym, n=1) * -1.0

            # cut egf
            if new_nt:
                egf = egf[:new_nt]

            # convolve egf with synthetics stf
            egf -= np.mean(egf)
            egf *= tukey(egf.size, 0.1)
            egf_fft = np.fft.fft(egf)
            egf = np.real(np.fft.ifft(egf_fft*syn_stf_fft))

            # save egf
            dataout = np.vstack((time_new, egf)).T
            outfile = os.path.join(event_outdir, os.path.basename(f))
            np.savetxt(outfile, dataout, ["%17.7f", "%17.7f"])

        print('done ', event)
