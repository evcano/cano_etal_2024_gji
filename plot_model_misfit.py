#!/usr/bin/env python3
import argparse
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.ticker as pltticker
import numpy as np
import os
from cmcrameri import cm
from glob import glob
from matplotlib.ticker import EngFormatter
from scipy.interpolate import LinearNDInterpolator


def cart2geo(x, z):
    deg2m = 111. * 1000.

    # array centroid
    lon0 = -118.57100841379312
    lat0 = 35.72925262068965

    lon = (x / deg2m) + lon0
    lat = (z / deg2m) + lat0
    return lon,lat


def read_mask(nproc, meshdir, maskdir, lim1, lim2):
    thr = 1.0

    # read mask
    xmask = np.array([])
    zmask = np.array([])
    mask = np.array([])

    for i in range(0, nproc):
        xcoor_proc = _read(os.path.join(meshdir, f'proc{i:06}_x.bin'))
        xmask = np.append(xmask, xcoor_proc)

        zcoor_proc = _read(os.path.join(meshdir, f'proc{i:06}_z.bin'))
        zmask = np.append(zmask, zcoor_proc)

        model_proc = _read(os.path.join(maskdir, f'proc{i:06}_gradient_mask.bin'))
        mask = np.append(mask, model_proc)

    # cut mask to desired limits
    if lim1 and lim2:
        idx1 = np.argwhere((xmask > -lim1) & (xmask < lim1)).flatten()
        idx2 = np.argwhere((zmask > -lim2) & (zmask < lim2)).flatten()
        idx3 = np.intersect1d(idx1, idx2).flatten()

        xmask = xmask[idx3]
        zmask = zmask[idx3]
        mask = mask[idx3]

    # sort mask and gll points
    isort = np.lexsort((xmask,zmask))
    xmask = xmask[isort]
    zmask = zmask[isort]
    mask = mask[isort]

    # apply threshold to mask
    mask[np.where(mask < thr)] = 0.0

    # get gll points on target region
    gll_target = np.argwhere(mask >= thr).flatten()

    # interpolate mask to a regular grid
    interp = LinearNDInterpolator(list(zip(xmask,zmask)), mask)

    nx = 800
    xregular = np.linspace(xmask.min(), xmask.max(), nx)
    zregular = np.linspace(zmask.min(), zmask.max(), nx)
    X, Z = np.meshgrid(xregular, zregular)

    MASKREG = interp(X.flatten(),Z.flatten())
    MASKREG = MASKREG.reshape((nx,nx), order="C")
    return X, Z, MASKREG, gll_target


def _read(filename, dtype='float32'):
    nbytes = os.path.getsize(filename)
    with open(filename, 'rb') as file:
        # read size of record
        file.seek(0)
        n = np.fromfile(file, dtype='int32', count=1)[0]

        if n == nbytes-8:
            file.seek(4)
            data = np.fromfile(file, dtype=dtype)
            return data[:-1]
        else:
            file.seek(0)
            data = np.fromfile(file, dtype=dtype)
            return data

fontsize = 8.0
figx = 4.0
figy = 2.5

majspa = 2.0
minspa = 0.2

# DO NOT EDIT BELOW THIS LINE
# ===========================
parser = argparse.ArgumentParser()
parser.add_argument('m1_mesh', type=str)
parser.add_argument('m1_model', type=str)
parser.add_argument('m2_mesh', type=str)
parser.add_argument('m2_model', type=str)
parser.add_argument('itrange', type=int, nargs=2)
parser.add_argument('lim1', type=float)
parser.add_argument('lim2', type=float)
parser.add_argument('--fig', type=str, default=None)
parser.add_argument('--sta', type=str, default=None)
parser.add_argument('--vmax', type=float, default=None)
parser.add_argument('--vmin', type=float, default=None)
args = parser.parse_args()

xcoor = np.array([])
zcoor = np.array([])
model = np.array([])

nproc = len(glob(os.path.join(args.m1_model, f"*vs.bin")))
iteration = np.arange(args.itrange[0], args.itrange[1]+1)
misfit = np.empty(len(iteration))

# read target model
for i in range(0, nproc):
    xcoor_proc = _read(os.path.join(args.m1_mesh,
                                    'proc{:06}_x.bin'.format(i)))
    zcoor_proc = _read(os.path.join(args.m1_mesh,
                                    'proc{:06}_z.bin'.format(i)))
    model_proc = _read(os.path.join(args.m1_model,
                                    'proc{:06}_vs.bin'.format(i)))

    xcoor = np.append(xcoor, xcoor_proc)
    zcoor = np.append(zcoor, zcoor_proc)
    model = np.append(model, model_proc)

# cut to desired limits
idx1 = np.argwhere((xcoor > -args.lim1) & (xcoor < args.lim1)).flatten()
idx2 = np.argwhere((zcoor > -args.lim2) & (zcoor < args.lim2)).flatten()
idx3 = np.intersect1d(idx1, idx2).flatten()

xcoor = xcoor[idx3]
zcoor = zcoor[idx3]
model = model[idx3]

# sort model and gll points
isort = np.lexsort((xcoor,zcoor))
xcoor = xcoor[isort]
zcoor = zcoor[isort]
model = model[isort]

# convert coor of model
lon, lat = cart2geo(xcoor,zcoor)

# read mask and get gll points where the model difference is computed
maskdir = "./array_mask/gradient_mask_complete/smoothed"
XMASK, ZMASK, MASKREG, gll_target = read_mask(nproc, args.m1_mesh, maskdir, args.lim1, args.lim2)

# convert coor of mask
XMASK, ZMASK = cart2geo(XMASK, ZMASK)

for j, it in enumerate(iteration):
    # read inverted model
    xcoor2 = np.array([])
    zcoor2 = np.array([])
    model2 = np.array([])

    for i in range(0, nproc):
        xcoor_proc = _read(os.path.join(args.m2_mesh,
                                        'proc{:06}_x.bin'.format(i)))
        zcoor_proc = _read(os.path.join(args.m2_mesh,
                                        'proc{:06}_z.bin'.format(i)))
        if it == 0:
            model_name = "MODEL_INIT"
        else:
            model_name = f"MODEL_{it:02d}"

        try:
            model_proc = _read(os.path.join(args.m2_model, model_name,
                                            'proc{:06}_vs.bin'.format(i)))
        except Exception:
            continue

        xcoor2 = np.append(xcoor2, xcoor_proc)
        zcoor2 = np.append(zcoor2, zcoor_proc)
        model2 = np.append(model2, model_proc)

    # cut to desired limits
    idx1 = np.argwhere((xcoor2 > -args.lim1) & (xcoor2 < args.lim1)).flatten()
    idx2 = np.argwhere((zcoor2 > -args.lim2) & (zcoor2 < args.lim2)).flatten()
    idx3 = np.intersect1d(idx1, idx2).flatten()

    xcoor2 = xcoor2[idx3]
    zcoor2 = zcoor2[idx3]
    model2 = model2[idx3]

    # sort model and gll points
    isort = np.lexsort((xcoor2,zcoor2))
    xcoor2 = xcoor2[isort]
    zcoor2 = zcoor2[isort]
    model2 = model2[isort]

    # compute model misfit
    misfit[j] = np.sqrt(np.mean(np.square(model2[gll_target]-model[gll_target])))
    print('Iteration ', it, 'model misfit: ', misfit[j])

    # plot relative perturbations
    mdiff = (model2 - model) / model
    mdiff *= 100.0

    if it == 0:
        mdiff0 = mdiff.copy()

    # FIGURE
    mercator = ccrs.PlateCarree()
    extent = [lon.min(), lon.max(), lat.min(), lat.max()]

    fig, ax = plt.subplots(1, 1, figsize=(figx, figy),
      subplot_kw={"projection": mercator})

    ax.set_extent(extent, mercator)
    ax.gridlines(alpha=0.0)

    # map
    cmap = cm.vik.resampled(31)

    im = ax.tripcolor(lon, lat, mdiff, cmap=cmap,
                      linewidth=0.0, edgecolor='face', shading="gouraud",
                      vmin=args.vmin, vmax=args.vmax, transform=mercator)

    # mask
    ax.pcolormesh(XMASK, ZMASK, MASKREG, alpha=0.1, cmap="gist_gray",
        transform=mercator)

    # stations
    if args.sta:
        sta = np.loadtxt(args.sta, usecols=(2,3), dtype=float)
        sta[:,0], sta[:,1] = cart2geo(sta[:,0],sta[:,1])
        ax.scatter(sta[:,0], sta[:,1], c='magenta', s=1.0, transform=mercator)

    # cities
    cities = ['Los Angeles', 'Sacramento', 'Fresno']
    ccoor = [[-118.24,34.05], [-121.49,38.58], [-119.78,36.73]]

    if not args.lim1 and not args.lim2:
        cities2 = ['Mexico City', 'Vancouver', 'Denver', 'Houston', 'Winnipeg',
             'Edmonton','Chihuahua']

        ccoor2 = [[-99.13, 19.43], [-123.12,49.28], [-104.99,39.73],[-95.36,29.76],
             [-97.13,49.89],[-113.49,53.54],[-106.05,28.64]]

        cities.extend(cities2)
        ccoor.extend(ccoor2)

    for cn, city in enumerate(cities):
        txtof = 0.1
        
        plt.text(ccoor[cn][0], ccoor[cn][1]+txtof, city,
            horizontalalignment='center', fontsize=fontsize*0.63, color="k",
            transform=mercator
            )

        plt.plot(ccoor[cn][0], ccoor[cn][1], marker=".", color="k",
                 markerfacecolor="k", markersize=1.0,
                 transform=mercator
                 )

    # political borders and coastline
    borders = cfeature.NaturalEarthFeature(
        category="cultural",
        name="admin_0_boundary_lines_land",
        scale="10m",
        facecolor="none",
        )

    ax.add_feature(borders, edgecolor="k", lw=0.3, zorder=4)
    ax.coastlines(resolution="10m",lw=0.3, color="k")

    # axes ticks
    ax.set_xticks(np.arange(lon.min(),lon.max()), crs=mercator)
    ax.xaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°", sep=""))
    ax.xaxis.set_major_locator(pltticker.MultipleLocator(base=majspa))
    ax.xaxis.set_minor_locator(pltticker.MultipleLocator(base=minspa))
    ax.tick_params(axis="x", labelsize=fontsize)

    ax.set_yticks(np.arange(lat.min(),lat.max()), crs=mercator)
    ax.yaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°", sep=""))
    ax.yaxis.set_major_locator(pltticker.MultipleLocator(base=majspa))
    ax.yaxis.set_minor_locator(pltticker.MultipleLocator(base=minspa))
    ax.tick_params(axis="y", labelsize=fontsize)

    # colorbar
    cbar = plt.colorbar(mappable=im, pad=0.02, fraction=0.03)
    cbar.ax.tick_params(labelsize=fontsize)

    # title and labels
    cbar.set_label("$\Delta$v/v [%]", fontsize=fontsize)

    ax.set_xlabel("longitude", fontsize=fontsize)
    ax.set_ylabel("latitude", fontsize=fontsize)


    # plot/save figure
    if args.fig:
        plt.savefig(os.path.join(args.fig,
                                 f"model_error_{it:02d}.png"), dpi=300,
                                 bbox_inches="tight")
    else:
        plt.show()
    plt.close()

#    # histogram of model errors
#    mean0 = np.mean(mdiff0[gll_target])
#    std0 = np.std(mdiff0[gll_target])
#    mean1 = np.mean(mdiff[gll_target])
#    std1 = np.std(mdiff[gll_target])
#
#    fig, ax = plt.subplots()
#    fig.set_size_inches(2.6,3.0)
#
#    maxval = np.max(np.abs(mdiff0[gll_target]))
#    bins = np.arange(-np.ceil(maxval),np.ceil(maxval), 1.0)
#
#    ax.hist(mdiff0[gll_target],bins,histtype="bar",
#            color="lightskyblue",edgecolor=None,alpha=1.0,linewidth=1.5)
#    ax.hist(mdiff[gll_target],bins,histtype="bar",
#            color="lightsalmon",edgecolor="red",alpha=1.0,linewidth=1.0)
#    ax.hist(mdiff0[gll_target],bins,histtype="step",
#            color="lightskyblue",edgecolor="blue",alpha=1.0,linewidth=1.5)
#
#    ytext = 3500
#    ax.text(-10,ytext,
#            f"$\mu$: {mean0:.2f}\n$\sigma$: {std0:.2f}",
#            color="blue",fontsize=fontsize*0.95)
#
#    ax.text(5,ytext,
#            f"$\mu$: {mean1:.2f}\n$\sigma$: {std1:.2f}",
#            color="red",fontsize=fontsize*0.95)
#
#    ax.set_xlabel("$\Delta$v/v [%]", fontsize=fontsize)
#    ax.set_ylabel("count", fontsize=fontsize)
#
#    ax.tick_params(axis="x",labelsize=fontsize)
#    ax.tick_params(axis="y",labelsize=fontsize)
#
#    ax.set_xlim(-15,15)
#    ax.set_ylim(0.0, 4000)
#
#    plt.grid(linestyle="--")
#    plt.savefig(os.path.join(args.fig,
#                f"model_error_hist_{it:02d}.png"),dpi=300,bbox_inches="tight")
#    plt.close()
#
## save/plot misifit evolution
#np.savetxt(os.path.join(args.m2_model, f"model_misfit_evolution.txt"),
#           np.array([iteration, misfit]).T)
#
#plt.plot(iteration, misfit, 'k-o')
#plt.xlabel('iteration')
#plt.ylabel('model misfit [m/s]')
#plt.grid()
#
#if args.fig:
#    plt.savefig(os.path.join(args.fig,
#                             f"model_misfit.png"))
#else:
#    plt.show()
#plt.close()
