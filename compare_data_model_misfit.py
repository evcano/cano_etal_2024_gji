#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os


fontsize = 8
figx = 3.5
figy = 2.2

parser = argparse.ArgumentParser()
parser.add_argument("--fig", type=str, default=None)
args = parser.parse_args()

path1 = "inversions_noise_model_2/anat"
path2 = "inversions_noise_model_2/fwani"
path3 = "inversions_noise_model_2/fdfwani"
all_paths = [path1, path2, path3]

legends = ["Inversion I", "Inversion II", "Inversion III"]
colors = ["-Dk", "-Db", "-Dr"]

stop_criterion = 0.00001

# DONT EDIT BELOW THIS LINE
fig1, ax = plt.subplots()
fig2, ax2 = plt.subplots()

for i in range(len(all_paths)):
    it_offset = 0
    legend = [legends[i], None]

    for j, dir_ in enumerate(["output"]):
        # plot data misfit
        try:
            tmp = np.loadtxt(os.path.join(all_paths[i], dir_,
                             "residuals", "data_misfit_evolution.txt"))
        except:
            continue

        iteration = tmp[:,0] + it_offset
        misfit = tmp[:,1]
        nfac = misfit[0]
        misfit /= nfac  # normalize the misfit

        # calculate misfit decrease
        decrease = np.diff(misfit[::-1])[::-1]
        decrease = np.append(1.0, decrease)

        # determine iteration where convergence was reached
        if type(stop_criterion) == int:
            idx = stop_criterion
        else:
            idx = np.argwhere(decrease <= stop_criterion).flatten()

            if idx.size > 0:
                idx = idx[0] + 1
            else:
                idx = iteration.size

        iteration = iteration[:idx]
        misfit = misfit[:idx] * nfac # scale the misfit back

        ax.plot(iteration, misfit, f'{colors[i]}',
                label=legend[j], linewidth=0.9, markersize=1.5)

        # plot model misfit
        try:
            tmp = np.loadtxt(os.path.join(all_paths[i], dir_,
                             "model_misfit_evolution.txt"))
        except:
            continue

        iteration = tmp[:,0]
        misfit = tmp[:,1]
        nfac = misfit[0]
        misfit /= nfac  # normalize the misfit

        iteration = iteration[:idx] + it_offset
        misfit = misfit[:idx] * nfac # scale the misfit back

        ax2.plot(iteration, misfit, f'{colors[i]}',
                label=legend[j], linewidth=0.9, markersize=1.5)

        it_offset = iteration[-1] + 1

xmax = 35

#ax.set_title('data misfit', fontsize=fontsize)
ax.set_xlabel('iteration',fontsize=fontsize)
ax.set_ylabel('misfit [s$^2$]',fontsize=fontsize)
ax.set_xticks(np.arange(0, xmax, 4))
ax.tick_params(axis="x", labelsize=fontsize)
ax.tick_params(axis="y", labelsize=fontsize)
#ax.legend(prop={"size": fontsize}, loc="upper center", ncol=3,
#          frameon=True, framealpha=1.0)
ax.grid()

#ax2.set_title('model misfit', fontsize=fontsize)
ax2.set_xlabel('iteration',fontsize=fontsize)
ax2.set_ylabel('misfit [m/s]',fontsize=fontsize)
ax2.set_xticks(np.arange(0, xmax, 4))
ax2.tick_params(axis="x", labelsize=fontsize)
ax2.tick_params(axis="y", labelsize=fontsize)
#ax2.legend(prop={"size": fontsize}, loc="upper right")
ax2.grid()

fig1.set_size_inches(figx, figy)
fig2.set_size_inches(figx, figy)

if args.fig:
    fig1.savefig(os.path.join(args.fig, "data_misfit.png"), dpi=300, bbox_inches="tight")

    fig2.savefig(os.path.join(args.fig, "model_misfit.png"), dpi=300, bbox_inches="tight")
else:
    plt.show()

plt.close()
