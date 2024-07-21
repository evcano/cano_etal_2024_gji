#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob
from scipy.stats import skew

#fontsize = 8
fontsize = 11
figx = 2.6
figy = 3.0

ylim1 = 0

#ylim2 = 1250
ylim2 = 2800

#ytext = 950
ytext = 2100

parser = argparse.ArgumentParser()
parser.add_argument("indir1", type=str)
parser.add_argument("it1", type=int)
parser.add_argument("--indir2", type=str, default=None)
parser.add_argument("--it2", type=int, default=None)
parser.add_argument("--fig", type=str, default=None)
args = parser.parse_args()


# by seisflows design, the file <residuals_it_0> contains
# the residuals of model <it - 1>, e.g, <residuals_1_0>
# contains the residuals of the initial model

# initial residuals
residuals1 = np.array([])

files1 = glob(os.path.join(args.indir1,
                           f"residuals_*_{args.it1+1}_0.txt"))

for file_ in files1:
    tmp = np.loadtxt(file_, dtype=float)
    residuals1 = np.append(residuals1, tmp)

mean1 = np.mean(residuals1)
std1 = np.std(residuals1)

print('total', np.size(residuals1))
print('unique', np.size(np.unique(residuals1)))

print("skew 1: ", skew(residuals1))

# final residuals
if args.indir2:
    residuals2 = np.array([])

    files2 = glob(os.path.join(args.indir2,
                               f"residuals_*_{args.it2+1}_0.txt"))

    for file_ in files2:
        tmp = np.loadtxt(file_, dtype=float)
        residuals2 = np.append(residuals2, tmp)

    mean2 = np.mean(residuals2)
    std2 = np.std(residuals2)

    print("skew 2: ", skew(residuals2))

# plot histogram
tmp1 = max(abs(residuals1.min()), abs(residuals1.max()))
if args.indir2:
    tmp2 = max(abs(residuals2.min()), abs(residuals2.max()))
    maxval = max(tmp1, tmp2)
else:
    maxval = tmp1

bins = np.arange(-np.ceil(maxval), np.ceil(maxval)+0.5, 0.5)

fig, ax = plt.subplots()

ax.hist(residuals1, bins, histtype="bar",
         color="lightskyblue", edgecolor=None, alpha=1.0, linewidth=1.5)
if args.indir2:
    ax.hist(residuals2, bins, histtype="bar",
             color="lightsalmon", edgecolor="red", alpha=1.0, linewidth=1.0)
ax.hist(residuals1, bins, histtype="step",
         color="lightskyblue", edgecolor="blue", alpha=1.0, linewidth=1.5)

ax.text(-5.0, ytext, 
        f"N: {residuals1.size}\n $\mu$: {mean1:.2f}\n $\sigma$: {std1:.2f}",
         color="blue",fontsize=fontsize*.9)

if args.indir2:
    ax.text(1.3, ytext,
            f"N: {residuals2.size}\n $\mu$: {mean2:.2f}\n $\sigma$: {std2:.2f}",
             color="red",fontsize=fontsize*.9)

ax.set_xlabel("$\Delta$T [s]",fontsize=fontsize)
ax.set_ylabel("count",fontsize=fontsize)

ax.set_xticks(np.arange(-6,7,2))

ax.tick_params(axis="x",labelsize=fontsize)
ax.tick_params(axis="y",labelsize=fontsize)

ax.set_xlim(-6, 6)
ax.set_ylim(ylim1,ylim2)

plt.grid(linestyle="--")

fig.set_size_inches(figx,figy)

if args.fig:
    plt.savefig(args.fig, dpi=300, bbox_inches="tight")
else:
    plt.show()
plt.close()
