#!/usr/bin/env python3
import argparse
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmasher as cmr
import matplotlib.colors as colors
import matplotlib.ticker as pltticker
import matplotlib.pyplot as plt
import numpy as np
import os
from cartopy.mpl.ticker import LongitudeLocator, LatitudeLocator
from cmcrameri import cm
from glob import glob
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator


def cart2geo(x, z):
    # array centroid
    lon0 = -118.57100841379312
    lat0 = 35.72925262068965

    deg2m = 111. * 1000.
    lon = (x / deg2m) + lon0
    lat = (z / deg2m) + lat0
    return lon,lat


def hessian_color_map(colormap):
    col = ["white", "yellow", "darkorange", "red", "black"]
    pts = [0.0, 0.25, 0.5, 0.75, 1.0]
    cmap = colors.LinearSegmentedColormap.from_list("model",
                                                    list(zip(pts, col)),
                                                    N=81)
    return cmap


def read_mask(nproc, meshdir, maskdir, lim1, lim2):
    thr = 1.

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

    # apply threshold to the mask
    mask[np.where(mask < thr)] = 0.0

    # interpolate mask to a regular grid
    interp = LinearNDInterpolator(list(zip(xmask, zmask)), mask)

    nx = 3000
    xregular = np.linspace(xmask.min(), xmask.max(), nx)
    zregular = np.linspace(zmask.min(), zmask.max(), nx)
    X, Z = np.meshgrid(xregular, zregular)

    MASKREG = interp(X.flatten(),Z.flatten())
    MASKREG = MASKREG.reshape((nx,nx), order="C")
    return X, Z, MASKREG


def read_model(nproc, suffix, mesh_dir, model_dir):
    xcoor = np.array([])
    zcoor = np.array([])
    model = np.array([])

    for i in range(0, nproc):
        if "mask_noise" in args.suffix:
            mesh_glob = np.loadtxt(os.path.join(args.mesh_dir,
                                                f"mesh_glob{i:06}"))
            xcoor_proc = mesh_glob[:, 0]
            zcoor_proc = mesh_glob[:, 1]
        else:
            xcoor_proc = _read(os.path.join(args.mesh_dir,
                                            f'proc{i:06}_x.bin'))
            zcoor_proc = _read(os.path.join(args.mesh_dir,
                                            f'proc{i:06}_z.bin'))

        model_proc = _read(os.path.join(args.model_dir,
                                        f"proc{i:06}_{args.suffix}.bin"))

        xcoor = np.append(xcoor, xcoor_proc)
        zcoor = np.append(zcoor, zcoor_proc)
        model = np.append(model, model_proc)
    return xcoor, zcoor, model


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

# for all figures
mask_flag = False
trunc_mesh_border = False
map_size = "small"
fontsize = 8.0
bcol = "k"

if map_size == "small":
    # for small figures
    figx = 4.0
    figy = 2.5
    majspa = 2.0
    minspa = 0.2
    marker_size = 1.0
    bordscale="10m"
elif map_size == "big":
    # for big figures
    figx = 4.0
    figy = 3.5
    majspa = 10.0
    minspa = 10.0
    marker_size = 0.3
    bordscale="50m"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('mesh_dir', type=str)
    parser.add_argument('model_dir', type=str)
    parser.add_argument('suffix', type=str)
    parser.add_argument('--fig', type=str, default=None)
    parser.add_argument('--lim1', type=float, default=None)
    parser.add_argument('--lim2', type=float, default=None)
    parser.add_argument('--sta', type=str, default=None)
    parser.add_argument('--vmax', type=float, default=None)
    parser.add_argument('--vmin', type=float, default=None)
    parser.add_argument('--title',type=str, default=None)

    args = parser.parse_args()
    nproc = len(glob(os.path.join(args.model_dir, f"*{args.suffix}.bin")))

    # read model
    xcoor, zcoor, model = read_model(nproc, args.suffix, args.mesh_dir,
                                     args.model_dir)

    # normalize model
    if args.suffix == "Hessian2_kernel":
        model /= np.max(model)
        model += np.max(model)*0.1
        model = 1.0 / model
    elif args.suffix == "mask_noise":
        model /= np.max(model)

    # cut model
    if args.lim1:
        idx1 = np.argwhere((xcoor > -args.lim1) & (xcoor < args.lim1)).flatten()
        idx2 = np.argwhere((zcoor > -args.lim2) & (zcoor < args.lim2)).flatten()
        idx3 = np.intersect1d(idx1, idx2).flatten()
    else:
        idx3 = np.arange(0, model.size)

    xcoor = xcoor[idx3]
    zcoor = zcoor[idx3]
    model = model[idx3]

    isort = np.lexsort((xcoor, zcoor))
    xcoor = xcoor[isort]
    zcoor = zcoor[isort]
    model = model[isort]

    # convert coordinates to cartesian
    xcoor, zcoor = cart2geo(xcoor, zcoor)

    # model information
    print('model ngll: ', model.size)
    print('model norm l2: ', np.linalg.norm(model, ord=2))
    print('model mean: ', np.sum(model))
    print('model norm inf: ', np.max(abs(model)))
    print('model min: ', np.min(model))
    print('model max: ', np.max(model))


    # FIGURE
    # ===================================================
    mercator = ccrs.PlateCarree()
    extent = [xcoor.min(), xcoor.max(), zcoor.min(), zcoor.max()]

    fig, ax = plt.subplots(1, 1, figsize=(figx, figy),
        subplot_kw={"projection": mercator}
        )

    ax.set_extent(extent, mercator)
    ax.gridlines(alpha=0.0)

    # map value limits
    vmax = np.max(np.abs(model))
    vmin = -vmax

    if args.vmax:
        vmax = args.vmax
    if args.vmin:
        vmin = args.vmin
    if args.suffix == "mask_noise":
        vmin = 0.
        vmax = 1.

    # colormap
    ncol = 30
    if args.suffix == "Hessian2_kernel":
        cmap = hessian_color_map("hessian")
    elif 'kernel' in args.suffix:
        cmap = cm.vik
    elif args.suffix == "mask_noise":
        cmap = cm.devon.resampled(ncol)
    else:
        cmap = cm.roma.resampled(ncol)

    # map
    im = ax.tripcolor(xcoor, zcoor, model, cmap=cmap,
                      linewidth=0.0, edgecolor="face", shading="gouraud",
                      vmin=vmin, vmax=vmax, transform=mercator)

    # mask
    if mask_flag:
        maskmesh = "./specfem2d_complete_mesh/DATA"
        maskdir = "./array_mask/gradient_mask_complete/smoothed"
        XMASK, ZMASK, MASKREG = read_mask(nproc, maskmesh, maskdir, args.lim1, args.lim2)
        XMASK, ZMASK = cart2geo(XMASK, ZMASK)

#        MASKREG = np.ma.masked_where(MASKREG > 1.0 - 1.E-3, MASKREG)

        ax.pcolormesh(XMASK, ZMASK, MASKREG, alpha=0.1, cmap="gist_gray", transform=mercator)

    # stations
    if args.sta:
        sta = np.loadtxt(args.sta, usecols=(2,3), dtype=float)
        sta[:,0], sta[:,1] = cart2geo(sta[:,0],sta[:,1])

        ax.scatter(sta[:,0], sta[:,1], c='magenta', s=marker_size,
            transform=mercator)

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
        if map_size == "big":
            txtof = 0.5
        else:
            txtof = 0.1
        
        plt.text(ccoor[cn][0], ccoor[cn][1]+txtof, city,
            horizontalalignment='center', fontsize=fontsize*0.63, color=bcol,
            transform=mercator
            )

        plt.plot(ccoor[cn][0], ccoor[cn][1], marker=".", color=bcol,
                 markerfacecolor=bcol, markersize=marker_size,
                 transform=mercator
                 )

    # political borders and coastline
    borders = cfeature.NaturalEarthFeature(
        category="cultural",
        name="admin_0_boundary_lines_land",
        scale=bordscale,
        facecolor="none",
        )

    ax.add_feature(borders, edgecolor=bcol, lw=0.3, zorder=4)
    ax.coastlines(resolution=bordscale,lw=0.3, color=bcol)

    # truncated noise source location
    # the truncated mesh is a square with sides of 1700 km 
    if trunc_mesh_border:
        c1_x, c1_z = cart2geo(-850.E+3, -850.E+3)
        c2_x, c2_z = cart2geo(-850.E+3, 850.E+3)
        c3_x, c3_z = cart2geo(850.E+3, 850.E+3)
        c4_x, c4_z = cart2geo(850.E+3, -850.E+3)
        plt.plot([c1_x, c2_x], [c1_z, c2_z],color="red",transform=mercator,lw=0.3)
        plt.plot([c2_x, c3_x], [c2_z, c3_z],color="red",transform=mercator,lw=0.3)
        plt.plot([c3_x, c4_x], [c3_z, c4_z],color="red",transform=mercator,lw=0.3)
        plt.plot([c4_x, c1_x], [c4_z, c1_z],color="red",transform=mercator,lw=0.3)

    # axes ticks
    ax.set_xticks(np.arange(xcoor.min(),xcoor.max()), crs=mercator)
    ax.xaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°", sep=""))
    ax.xaxis.set_major_locator(pltticker.MultipleLocator(base=majspa))
    ax.xaxis.set_minor_locator(pltticker.MultipleLocator(base=minspa))
    ax.tick_params(axis="x", labelsize=fontsize)

    ax.set_yticks(np.arange(zcoor.min(),zcoor.max()), crs=mercator)
    ax.yaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°", sep=""))
    ax.yaxis.set_major_locator(pltticker.MultipleLocator(base=majspa))
    ax.yaxis.set_minor_locator(pltticker.MultipleLocator(base=minspa))
    ax.tick_params(axis="y", labelsize=fontsize)

    # colorbar
    cbar = plt.colorbar(mappable=im, pad=0.02, fraction=0.03)
    cbar.ax.tick_params(labelsize=fontsize)

    # title and labels
    if args.title:
        ax.set_title(args.title, fontsize=fontsize)

    if args.suffix == "vs":
        cbar.set_label("shear wave speed [m/s]", fontsize=fontsize)
    elif args.suffix == "mask_noise":
        cbar.set_label("normalized PSD", fontsize=fontsize)
    else:
        cbar.set_label("mask",fontsize=fontsize)

    ax.set_xlabel("longitude", fontsize=fontsize)
    ax.set_ylabel("latitude", fontsize=fontsize)

    # save figure
    if args.fig:
        plt.savefig(args.fig, dpi=300, bbox_inches="tight")
    else:
        plt.show()

    plt.close()
