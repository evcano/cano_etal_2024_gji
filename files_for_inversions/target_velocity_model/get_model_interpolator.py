import h5py
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.interpolate import (LinearNDInterpolator,
                               NearestNDInterpolator,
                               CloughTocher2DInterpolator)

# specfem mesh origin (centroid of the array) in geo coordinates
lon0 = -118.57100841379312
lat0 = 35.72925262068965

# array limits (used to compute the mean velocity inside the array)
mmlonmin = -122.814057
mmlonmax = -114.935505
mmlatmin = 33.075909
mmlatmax = 38.312408

# DONT EDIT BELOW THIS LINE
model_file = './model_slice_24000.0.npy'
model = np.load(model_file)

# add missing model point
model = np.vstack((model, np.array([-150, 3.0, 4380])))


# convert coordinates from spherical to mesh
lon = (model[:, 0] - lon0) * 111. * 1000.
lat = (model[:, 1] - lat0) * 111. * 1000.

print(lon.min(), lon.max())
print(lat.min(), lat.max())
print('model min/max: ', model[:, 2].min(), model[:, 2].max())

interpolator = LinearNDInterpolator(list(zip(lon,lat)), model[:, 2])

with open(f'vel_interpolator.pkl', 'wb') as _file:
    pickle.dump(interpolator,_file)

# mean model
idx1 = np.argwhere((model[:, 0] >= mmlonmin) & (model[:, 0] <= mmlonmax))
idx2 = np.argwhere((model[:, 1] >= mmlatmin) & (model[:, 1] <= mmlatmax))
idx = np.intersect1d(idx1, idx2)
print('mean velocity inside the array: ', np.mean(model[idx, 2]))
print('min velocity inside the array: ', np.min(model[idx, 2]))
print('max velocity inside the array: ', np.max(model[idx, 2]))
