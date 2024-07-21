import h5py
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.interpolate import LinearNDInterpolator


def cart2sph(x, y, z):
    """
    Given cartesian coordinates, returns their spherical equivalent.
    :param x: x.
    :param y: y.
    :param z: z.
    :return: colatitude, longitude, and radius
    """

    x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    # Handle division by zero at the core
    with np.errstate(invalid='ignore'):
        c = np.divide(z, r)
        c = np.nan_to_num(c)

    c = np.arccos(c)
    l = np.arctan2(y, x)
    return np.rad2deg(c), np.rad2deg(l), r


slice_depth = 24000.0
depth_tolerance = 500.0
earthradius = 6.371e+06

# read hdf5 file
model_file = './CSEM_for_Eduardo_10s.h5'
h5file = h5py.File(model_file, 'r')
coor = np.array(h5file['MODEL/coordinates'])  # (NEL,NGLL,3)
model = np.array(h5file['MODEL/data'])  # (NEL,9,NGLL)
h5file.close()

# reshape arrays
ngll_glob = coor.shape[0] * coor.shape[1]
coor = np.reshape(coor, (ngll_glob, 3), order='C')
model = np.swapaxes(model,1,2) # (NEL,NGLL,9)
model = np.reshape(model, (ngll_glob, 9), order='C')

vsv = model[:,7]

# convert coordinates from cartesian to spherical
colat, lon, rad = cart2sph(coor[:,0],coor[:,1],coor[:,2])
lat = 90.0 - colat
depth = earthradius - rad

print('lon min/max ',np.min(lon),np.max(lon))
print('lat min/max ',np.min(lat),np.max(lat))
print('depth min/max ',np.min(depth),np.max(depth))

# get model slice
idx = np.argwhere(np.abs(depth - slice_depth) < depth_tolerance)
model_slice = np.zeros((idx.size, 3))
model_slice[:, 0] = lon[idx].flatten()
model_slice[:, 1] = lat[idx].flatten()
model_slice[:, 2] = vsv[idx].flatten()

np.save(f'model_slice_{slice_depth}', model_slice)
