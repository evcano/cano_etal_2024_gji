#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob
from obspy.core import Trace
from scipy.signal import correlation_lags, tukey
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

def compute_window(tr, is_corr):
    wdur = 100.0

    dt = tr.stats.delta
    nt = tr.stats.npts
    whlen = int((wdur/2)/dt)

    if is_corr:
        zlag = int((nt-1)/2)
        neg = tr.data.copy()
        neg[zlag:] *= 0.0
        wcen = np.argmax(np.square(neg))
        wneg = [int(wcen-whlen), int(wcen+whlen)]

        pos = tr.data.copy()
        pos[:zlag] *= 0.0
        wcen = np.argmax(np.square(pos))
        wpos = [int(wcen-whlen), int(wcen+whlen)]
        return [wneg, wpos]
    else:
        wcen = np.argmax(np.square(tr.data))
        win = [int(wcen-whlen), int(wcen+whlen)]
        return win

def check_window(tr_obs, tr_syn, win, is_corr):

    def compute_cc_delay(obs, syn, dt):
        # correlation as in seisflows
        corr = np.convolve(obs, np.flipud(syn))
        delay = (np.argmax(corr) - obs.size + 1) * dt
        #correlation coefficient
        corr /= (np.linalg.norm(syn,ord=2) * np.linalg.norm(obs,ord=2))
        cc = np.max(corr)
        return cc, delay

    def compute_snr(u, win, winnoise):
        signal = u[win[0]:win[1]]
        noise = u[winnoise[0]:winnoise[1]]
        snr = np.max(np.abs(signal)) / np.sqrt(np.mean(np.square(noise)))
        if snr < 0:
            snr = -10.0 * np.log10(-snr)
        elif snr > 0:
            snr = 10.0 * np.log10(snr)
        return snr

    dthr = 10.0
    ccthr = 0.7
    snrthr = 0.0

    tperc = 0.1
    dt = tr_obs.stats.delta
    nt = tr_obs.stats.npts
    zlag = int((nt-1)/2)

    if is_corr:
        # negative branch
        wneg = win[0]
        wdur = wneg[1] - wneg[0]
        noff = int(wdur/2)
        wneg_noise = [wneg[0]-wdur-noff, wneg[0]-noff]

        obs_neg = tr_obs.data[wneg[0]:wneg[1]].copy()
        syn_neg = tr_syn.data[wneg[0]:wneg[1]].copy()
        obs_neg *= tukey(obs_neg.size, tperc)
        syn_neg *= tukey(syn_neg.size, tperc)

        cc_neg, delay_neg = compute_cc_delay(obs=obs_neg, syn=syn_neg, dt=dt)
        snr_neg = compute_snr(tr_obs.data, wneg, wneg_noise)

        if (abs(delay_neg) > dthr or cc_neg < ccthr or snr_neg < snrthr
            or wneg[1]>zlag):
            delay_neg, cc_neg, snr_neg = None, None, None
    
        # positive branch
        wpos = win[1]
        wdur = wpos[1] - wpos[0]
        wpos_noise = [wpos[1]+noff, wpos[1]+noff+wdur]

        obs_pos = tr_obs.data[wpos[0]:wpos[1]].copy()
        syn_pos = tr_syn.data[wpos[0]:wpos[1]].copy()
        obs_pos *= tukey(obs_pos.size, tperc)
        syn_pos *= tukey(syn_pos.size, tperc)

        cc_pos, delay_pos = compute_cc_delay(obs=obs_pos, syn=syn_pos, dt=dt)
        snr_pos = compute_snr(tr_obs.data, wpos, wpos_noise)

        if (abs(delay_pos) > dthr or cc_pos < ccthr or snr_pos < snrthr
            or wpos[0]<zlag):
            delay_pos, cc_pos, snr_pos = None, None, None

        cc = [cc_neg, cc_pos]
        delay = [delay_neg, delay_pos]
        snr = [snr_neg, snr_pos]

        return cc, delay, snr
    else:
        wdur = win[1] - win[0]
        noff = int(wdur/2)
        winnoise = [win[1]+noff, win[1]+noff+wdur]

        obs = tr_obs.data[win[0]:win[1]].copy()
        syn = tr_syn.data[win[0]:win[1]].copy()
        obs *= tukey(obs.size, tperc)
        syn *= tukey(syn.size, tperc)

        cc, delay = compute_cc_delay(obs=obs, syn=syn, dt=dt)
        snr = compute_snr(tr_obs.data, win, winnoise)

        return cc, delay, snr
    
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

#fontsize = 8
fontsize = 10
figx = 4.5
figy = 1.5

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('path1', type=str)
    parser.add_argument('path2', type=str)
    parser.add_argument('--fig', type=str, default=None)
    parser.add_argument('--filt', type=float, default=None, nargs=2)
    parser.add_argument('--tlim', type=float, default=None, nargs=2)
    parser.add_argument('--path3', type=str, default=None)
    parser.add_argument('--corr', action='store_true')
    parser.add_argument('--name', action='store_true')
    parser.add_argument('--norm', action='store_true')
    parser.add_argument('--misf', action='store_true')
    parser.add_argument('--wind', action='store_true')
    parser.add_argument('--title', type=str, default=None)
    args = parser.parse_args()

    files = glob(args.path1)
    files.sort()
    
    fig, ax = plt.subplots()
    fig.set_size_inches(figx,figy)

    offset = 0.

    for file_ in files:
        # read waveforms
        basename = os.path.basename(file_)
        tr1 = read_ascii_waveform(file_)
        tr2 = None
        tr3 = None

        file2_ = os.path.join(args.path2, basename)
        if os.path.exists(file2_):
            tr2 = read_ascii_waveform(file2_)
        else:
            print(f"file {file2_} does not exist\n")

        if args.path3:
            file3_ = os.path.join(args.path3, basename)
            if os.path.exists(file3_):
                tr3 = read_ascii_waveform(file3_)
            else:
                print(f"file {file3_} does not exist\n")

        # apply filter
        if args.filt:
            tr1.filter("bandpass", freqmin=args.filt[0], freqmax=args.filt[1],
                      corners=2, zerophase=True)
            if tr2:
                tr2.filter("bandpass", freqmin=args.filt[0], freqmax=args.filt[1],
                          corners=2, zerophase=True)
            if tr3:
                tr3.filter("bandpass", freqmin=args.filt[0], freqmax=args.filt[1],
                          corners=2, zerophase=True)
        # normalize
        if args.norm:
            fac = 1 / (np.max(tr1.data))
            tr1.data *= fac
            if tr2:
#                fac = (np.sum(tr1.data*tr2.data)*tr1.stats.delta) / (np.sum(tr2.data**2.0)*tr2.stats.delta)
                fac = 1 / (np.max(tr2.data))
                tr2.data *= fac
            if tr3:
#                fac = (np.sum(tr1.data*tr3.data)*tr1.stats.delta) / (np.sum(tr3.data**2.0)*tr3.stats.delta)
                fac = 1 / (np.max(tr3.data))
                tr3.data *= fac

        # define time axis
        if args.corr:
            hnpts = (tr1.stats.npts + 1.) / 2.
            times1 = correlation_lags(hnpts, hnpts, mode="full") * tr1.stats.delta
            if tr2:
                hnpts = (tr2.stats.npts + 1.) / 2.
                times2 = correlation_lags(hnpts, hnpts, mode="full") * tr2.stats.delta
            if tr3:
                hnpts = (tr3.stats.npts + 1.) / 2.
                times3 = correlation_lags(hnpts, hnpts, mode="full") * tr3.stats.delta
        else:
            times1 = tr1.times()
            if tr2:
                times2 = tr2.times()
            if tr3:
                times3 = tr3.times()

        # add offset to waveforms
        tr1.data += offset
        if tr2:
            tr2.data += offset
        if tr3:
            tr3.data += offset

        # plot waveforms
        ax.plot(times1, tr1.data, color='black', alpha=0.9, linewidth=0.7)
        if tr2:
            ax.plot(times2, tr2.data, color='blue', alpha=0.9, linewidth=0.7)
        if tr3:
            ax.plot(times3, tr3.data, color='red', alpha=0.9, linewidth=0.7)

        # plot waveform names
        if args.name:
            sta1 = os.path.basename(os.path.dirname(args.path1))
            sta2 = basename.removesuffix(".BXY.semd")

            if args.tlim:
                text_x = args.tlim[0] + 10 * tr1.stats.delta
            else:
                text_x = times1.min() + 10 * tr1.stats.delta

            text_y = max(np.max(tr1.data),np.max(tr2.data))

            plt.text(text_x, text_y, f"{sta2[5:]} - {sta1[5:]}",
                    fontsize=fontsize)
            if args.wind:
                text_y = min(np.min(tr1.data),np.min(tr2.data))
                minp = 1.0/args.filt[1]
                maxp = 1.0/args.filt[0]
#                plt.text(text_x, -1.5, f"Period band:\n{minp:.0f} - {maxp:.0f} s", fontsize=fontsize)

        # compute and plot misfit information
        if args.misf:
            win = compute_window(tr2, args.corr)
            cc, delay, snr = check_window(tr1, tr2, win, args.corr)
            if tr3:
                cc2, delay2, snr2 = check_window(tr1, tr3, win, args.corr)

            if args.corr:
                wneg = win[0]
                wpos = win[1]
                text_y = min(np.min(tr1.data),np.min(tr2.data))

                # negative branch
                if snr[0]:
                    text_x = times1[int((wneg[0]+wneg[1])*0.5)] + 10.0
                    plt.text(text_x, -1.5, f"$\Delta$T: {delay[0]:.2f} s",
                             color="blue",fontsize=fontsize)
                    if tr3:
                        plt.text(text_x, -1.9, f"$\Delta$T: {delay2[0]:.2f} s",
                                 color="red",fontsize=fontsize)
                    if args.wind:
                        text_x = times1[int((wneg[0]+wneg[1])*0.5)] + 10.0
#                        plt.text(text_x, -1.9, f"CC: {cc[0]:.2f}", color="blue",fontsize=fontsize)
                        ax.axvspan(xmin=times1[wneg[0]], xmax=times1[wneg[1]], color="grey",
                                   edgecolor=None, alpha=0.4)

                # positive branch
                if snr[1]:
                    text_x = times1[int((wpos[0]+wpos[1])*0.5)] + 10.0
                    plt.text(text_x, -1.5, f"$\Delta$T: {delay[1]:.2f} s",
                             color="blue",fontsize=fontsize)
                    if tr3:
                        plt.text(text_x, -1.9, f"$\Delta$T: {delay2[1]:.2f} s",
                                 color="red",fontsize=fontsize)
                    if args.wind:
                        text_x = times1[int((wpos[0]+wpos[1])*0.5)] + 10.0
#                        plt.text(text_x, -1.9, f"CC: {cc[1]:.2f}",color="blue",fontsize=fontsize)
                        ax.axvspan(xmin=times1[wpos[0]], xmax=times1[wpos[1]], color="grey",
                                   edgecolor=None, alpha=0.4)
            else:
                if snr:
                    text_y = np.min(tr1.data)
                    text_x = times1[int((win[0]+win[1])*0.5)] + 10.0
                    plt.text(text_x, -1.5, f"$\Delta$T: {delay:.2f} s",
                             color="blue",fontsize=fontsize)

                    if tr3:
                        plt.text(text_x, -1.9, f"$\Delta$T: {delay2:.2f} s",
                                 color="red",fontsize=fontsize)

                    if args.wind:
#                        plt.text(text_x, -1.9, f"CC: {cc:.2f}", color="blue",fontsize=fontsize)
                        ax.axvspan(xmin=times1[win[0]], xmax=times1[win[1]], color="grey",
                                   edgecolor=None, alpha=0.4)

        # update offset
        offset = np.max(tr1.data)
        if tr2:
            offset = max(offset, np.max(tr2.data))
        if tr3:
            offset = max(offset, np.max(tr3.data))

    # figure settings
    if args.tlim:
        ax.set_xlim(args.tlim[0], args.tlim[1])

    if args.corr:
        ax.set_xlabel('lag [s]',fontsize=fontsize)
    else:
        ax.set_xlabel('time [s]',fontsize=fontsize)

    if args.title:
        ax.set_title(args.title, fontsize=fontsize)

    ax.set_yticks([])
    ax.tick_params(axis="x",labelsize=fontsize)
    ax.set_ylim(-2,1.5)

    if args.fig:
        plt.savefig(args.fig, dpi=300, bbox_inches="tight")
    else:
        plt.show()

    plt.close()
