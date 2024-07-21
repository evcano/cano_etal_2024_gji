import numpy as np
import os

stations_file = './STATIONS'
outdir = './files_for_inversions/source_files_egf'

allsta = np.loadtxt(stations_file, usecols=(0), dtype=str)
allnet = np.loadtxt(stations_file, usecols=(1), dtype=str)
xcoor = np.loadtxt(stations_file, usecols=(2), dtype=str)
zcoor = np.loadtxt(stations_file, usecols=(3), dtype=str)

factor = -10.E+10
f0 =  0.1
time_function_type = 3

for i in range(0, len(allsta)):
    net = allnet[i]
    sta = allsta[i]
    x = xcoor[i]
    z = zcoor[i]

    force_file_name = f'SOURCE_{net}.{sta}'

    with open(os.path.join(outdir, force_file_name), 'w') as _file:
        _file.write(f'source_surf = .false.\n')
        _file.write(f'xs = {x}\n')
        _file.write(f'zs = {z}\n')
        _file.write(f'source_type = 1\n')
        _file.write(f'time_function_type = {time_function_type}\n')
        _file.write(f'name_of_source_file = YYYYYYYYY\n')
        _file.write(f'burst_band_width = 0.\n')
        _file.write(f'f0 = {f0}\n')
        _file.write(f'tshift = 0.0\n')
        _file.write(f'anglesource = 0.\n')
        _file.write(f'Mxx = 0.d0\n')
        _file.write(f'Mzz = 0.d0\n')
        _file.write(f'Mxz = 0.d0\n')
        _file.write(f'factor = {factor}\n')
        _file.write(f'vx = 0.0\n')
        _file.write(f'vz = 0.0\n')
