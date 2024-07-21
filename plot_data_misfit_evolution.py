#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob


parser = argparse.ArgumentParser()
parser.add_argument("indir", type=str)
parser.add_argument("itrange", type=int, nargs=2)
parser.add_argument("--fig", type=str, default=None)
args = parser.parse_args()

iteration = np.arange(args.itrange[0], args.itrange[1]+1)
misfit = np.empty(len(iteration))
nwindows = np.empty(len(iteration))
misfit_decrease = np.empty(len(iteration))

# by seisflows design, the file <residuals_it_0> contains
# the residuals of model <it - 1>, e.g, <residuals_1_0>
# contains the residuals of the initial model

for j, it in enumerate(iteration):
    fname = f'residuals_*_{it+1}_0.txt'

    try:
        files = glob(os.path.join(args.indir, fname))
    except Exception:
        continue

    residuals = np.array([])

    # read all residuals
    for file_ in files:
        tmp = np.loadtxt(file_, dtype=float)
        residuals = np.append(residuals, tmp)

    # number of windows
    if len(residuals) == 0:
        nwindows[j] = np.nan
    else:
        nwindows[j] = len(residuals)

    # misfit
    misfit[j] = 0.5 * (1./nwindows[j]) * np.sum(np.square(residuals))

    # misfit decrease
    if j == 0:
        mdecrease = 1.
    else:
        mdecrease = misfit[j-1] - misfit[j]

    misfit_decrease[j] = mdecrease

    print('')
    print('Iteration', it, 'misfit: ', misfit[j])
    print('Iteration', it, 'misfit decrease: ', mdecrease)
    print('Iteration', it, 'number of windows: ', nwindows[j])
    print('Iteration', it, 'number of events: ', len(files))
    print('Iteration', it, 'residuals mean: ', np.mean(abs(residuals)))
    print('Iteration', it, 'residuals STD: ', np.std(residuals))
    print('Iteration', it, 'residuals min/max: ', residuals.min(), residuals.max())

# save misfit evolution
np.savetxt(os.path.join(args.indir, f"data_misfit_evolution.txt"),
           np.array([iteration, misfit, misfit_decrease]).T)

# misfit evolution
fig, ax = plt.subplots()
ax2 = ax.twinx()
ax.plot(iteration, misfit, '-ob')
ax2.plot(iteration, nwindows, '-or')

ax.set_xlabel('iteration')
ax.set_ylabel('data misfit')
ax2.set_ylabel('windows')
plt.grid()

if args.fig:
    plt.savefig(os.path.join(args.fig, f"data_misfit.png"))
else:
    plt.show()

plt.close()
