#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import pickle
import sys
import yaml
import gstools as gs
from gstools import Exponential, SRF
from gstools.random import MasterRNG
from scipy.interpolate import RegularGridInterpolator


parfile = sys.argv[1]

try:
    with open(parfile, 'r') as _file:
        par = yaml.safe_load(_file)
except IOError:
    print('IOError: parfile not found.')

x = np.linspace(par['XMIN'], par['XMAX'], par['NX'], endpoint=True)
z = np.linspace(par['ZMIN'], par['ZMAX'], par['NZ'], endpoint=True)

model = Exponential(dim=2, var=par['STD']**2, len_scale=par['LEN_SCALE'])

if not par['SEED']:
    seed = MasterRNG(None)
    seed = seed()
else:
    seed = par['SEED']

srf = SRF(model, mean=par['MEAN'], seed=seed)
perturbations = srf.structured([x, z])

if par['MAXPER']:
    perturbations[perturbations >= 0.0] = par['MAXPER']
    perturbations[perturbations < 0.0] = -par['MAXPER']

if par['ONESIGN'] and par['ONESIGN'] == 'positive':
    perturbations[perturbations < 0.0] = 0.0
if par['ONESIGN'] and par['ONESIGN'] == 'negative':
    perturbations[perturbations > 0.0] = 0.0

print('random field seed: {}'.format(seed))
print('random field minimum value: {}'.format(perturbations.min()))
print('random field maximum value: {}'.format(perturbations.max()))

print('generating perturbations interpolator ...')
interpolator = RegularGridInterpolator((x, z), perturbations, method='linear')

# save iterpolator object
output_file = os.path.join(par['OUTPUT_DIR'], par['NAME'])

with open(f'{output_file}.pkl', 'wb') as _file:
    pickle.dump(interpolator, _file)

# plot field
cmap = cm.get_cmap('seismic')
plt.imshow(perturbations, cmap=cmap)
plt.colorbar()
plt.show()
