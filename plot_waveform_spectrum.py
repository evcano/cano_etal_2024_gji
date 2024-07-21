#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.filter import bandpass

def read_ascii_waveform(file):
    waveform = np.loadtxt(file, dtype="float")
    time = waveform[:,0]
    amp = waveform[:,1]
    return time, amp

fontsize = 8

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot one waveform.")
    parser.add_argument('file',help="Waveform file.")
    parser.add_argument('--f',type=float,nargs=2,default=None,help="Bandpass frequencies.")
    parser.add_argument('--fig',type=str,default=None)
    args = parser.parse_args()

    time, amp = read_ascii_waveform(args.file)
    DT = abs(time[1] - time[0])
    DF = 1.0 /  DT
    NS = len(time)

    if args.f:
       amp = bandpass(amp, freqmin=args.f[0], freqmax=args.f[1], df=DF,
                      corners=2, zerophase=True)

    fft = np.fft.rfft(amp)
    spec = np.abs(np.real(fft))
    fqaxis = np.fft.rfftfreq(n=NS, d=DT)

    # FIGURE
    fig, ax = plt.subplots()
    ax.plot(time, amp/abs(amp.max()), color="black", linewidth=0.7)

    ax.set_xlabel('time [s]',fontsize=fontsize)
    ax.set_ylabel('normalized amplitude',fontsize=fontsize)
    ax.tick_params(axis="x", labelsize=fontsize)
    ax.tick_params(axis="y", labelsize=fontsize)
    ax.set_xlim(-1000,1000)

    fig.tight_layout(pad=1.0)
    fig.set_size_inches(3.7, 1.2)

    if args.fig:
        plt.savefig(f"{args.fig}/noise_stf.png", dpi=300, bbox_inches="tight")
    else:
        plt.show()
    plt.close()

    # plot spectrum
    fig, ax = plt.subplots()
    ax.plot(fqaxis, spec/abs(spec.max()), c='k', linewidth=0.7)

    # show period axis
    secax = ax.secondary_xaxis("top", functions=(lambda x:1/x, lambda x:1/x))
    secax.set_xlabel("period [s]",fontsize=fontsize)
    xticks2 = np.append(np.arange(12,20,2), np.arange(20,35,5))
    xticks2 = np.append(xticks2, np.array([40, 50]))
    secax.set_xticks(xticks2)
    secax.tick_params(axis="x",labelsize=fontsize)

    ax.set_xlabel('frequency [Hz]',fontsize=fontsize)
    ax.set_ylabel('normalized amplitude',fontsize=fontsize)
    ax.tick_params(axis="x", labelsize=fontsize)
    ax.tick_params(axis="y", labelsize=fontsize)
    ax.set_xlim(0.02, 0.09)
    ax.yaxis.get_offset_text().set_fontsize(fontsize)

    fig.tight_layout(pad=1.0)
    fig.set_size_inches(3.7, 1.2)

    if args.fig:
        plt.savefig(f"{args.fig}/noise_spec.png", dpi=300, bbox_inches="tight")
    else:
        plt.show()
    plt.close()
