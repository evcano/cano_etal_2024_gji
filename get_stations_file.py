import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

stafile = './california-stations.txt'
half_size = 30.0  # half length of the mesh in degrees

def main():
    netcode = np.loadtxt(stafile, usecols=(0),
                         comments='#', delimiter='|', dtype='str')

    stacode = np.loadtxt(stafile, usecols=(1),
                         comments='#', delimiter='|', dtype='str')

    stations = np.loadtxt(stafile, usecols=(2,3),
                          comments='#', delimiter='|') 

    # mesh limits (the origin is the centroid of the array)
    lon0 = np.mean(stations[:, 1])
    lat0 = np.mean(stations[:, 0])

    latmin = lat0 - half_size 
    latmax = lat0 + half_size 

    lonmin = lon0 - half_size
    lonmax = lon0 + half_size

    print('array aperture lon: ', np.max(stations[:,1])-np.min(stations[:,1]))
    print('array aperture lat: ', np.max(stations[:,0])-np.min(stations[:,0]))
    print('array lon min/max: ', np.min(stations[:,1]), np.max(stations[:,1]))
    print('array lat min/max: ', np.min(stations[:,0]), np.max(stations[:,0]))

    print("")
    print('mesh origin (geo)')
    print('lon/lat: ', lon0, lat0)

    print("")
    print('mesh limits (geo)')
    print('min lon/lat: ', lonmin, latmin)
    print('max lon/lat: ', lonmax, latmax)
    print('lonlength/latlength: ', lonmax-lonmin, latmax-latmin)

    # convert station coordinates from geo to cartesian
    # the mesh dimensions do not affect the coordinates of the stations since
    # the origin is the centroid of the array
    deg2m = 111. * 1000.
    with open("STATIONS", "w") as _file:
        for i, stacoor in enumerate(stations):
            x = (stacoor[1] - lon0) * deg2m
            z = (stacoor[0] - lat0) * deg2m
            d = 0.0

            sta = stacode[i]
            net = netcode[i]

            _file.write(f"{sta} {net} {x:.6f} {z:.6f} {d} {d}\n")

    print("")
    print('mesh origin (cart)')
    print('lon/lat: ', (lon0-lon0)*deg2m, (lat0-lat0)*deg2m)

    print("")
    print('mesh limits (cart)')
    print('min lon/lat: ', (lonmin-lon0)*deg2m, (latmin-lat0)*deg2m)
    print('max lon/lat: ', (lonmax-lon0)*deg2m, (latmax-lat0)*deg2m)
    print('xlength/ylength: ', (lonmax-lonmin)*deg2m, (latmax-latmin)*deg2m)

if __name__ == '__main__':
    main()
