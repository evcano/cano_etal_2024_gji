import cmasher as cmr
import matplotlib.pyplot as plt
import numpy as np
import pickle
from geo import is_land
from noisi.util.plot import plot_grid
from scipy.interpolate import LinearNDInterpolator,NearestNDInterpolator


grd = np.load('mean_model_grd.npy')
src = np.load('mean_model.npy')

ampfact = 10E+10

# specfem mesh origin (centroid of the array) in geo coordinates
lon0 = -118.57100841379312
lat0 = 35.72925262068965

# region to add points in the continent
lonmin = -157.783
lonmax = -77.783
latmin = 0.249
latmax = 76.249

# DONT EDIT BELOW THIS LINE
grdlon = grd[0, :]
grdlat = grd[1, :]

# add grid points on the continents and set to zero
alllon = np.array([])
alllat = np.array([])
for x in np.linspace(lonmin,lonmax,500):
    for y in np.linspace(latmin,latmax,500):
        alllon = np.append(alllon, x)
        alllat = np.append(alllat, y)
landidx = is_land(alllon, alllat)
landidx = np.argwhere(landidx == 1.0)

for i in landidx:
    grdlon = np.append(grdlon, alllon[i])
    grdlat = np.append(grdlat, alllat[i])
    src = np.append(src, 0.0)

src *= ampfact

plot_grid(grdlon, grdlat, src, cmap=cmr.cosmic, size=100., sequential=True)

# convert to mesh coordinates
grdlon = (grdlon - lon0) * 111. * 1000.
grdlat = (grdlat - lat0) * 111. * 1000.

interpolator = LinearNDInterpolator(list(zip(grdlon, grdlat)), src)

with open(f'sans_interpolator.pkl', 'wb') as _file:
    pickle.dump(interpolator, _file)
