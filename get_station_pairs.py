import itertools
import matplotlib.pyplot as plt
import numpy as np
import os
from cmcrameri import cm


def az_two_points(x1, z1, x2, z2):
    dz = z2 - z1
    dx = x2 - x1

    # angle from 0 to 360 measured counterclock wise
    # from the positive x-axis where x1,z1 is the origin
    theta = (np.arctan2(dz, dx) + 2*np.pi) % (2*np.pi)
    theta = np.degrees(theta)

    # angle from 0 to 360 measured clock wise
    # from positive y-axis where x1,z1 is the origin
    az = (450.0 - theta) % 360.0
    return az


stafile = "./STATIONS"

# DONT EDIT BELOW THIS LINE
all_sta = np.loadtxt(stafile, usecols=(0), dtype=str)
all_net = np.loadtxt(stafile, usecols=(1), dtype=str)
all_xcoor = np.loadtxt(stafile, usecols=(2), dtype=float)
all_zcoor = np.loadtxt(stafile, usecols=(3), dtype=float)
nsta = len(all_sta)

sta_dic = {}
for i in range(0, nsta):
    staname = f"{all_net[i]}.{all_sta[i]}"
    sta_dic[staname] = {"x": all_xcoor[i], "z": all_zcoor[i]}

stapairs = itertools.combinations(sta_dic.keys(), 2)

with open("STATION_PAIRS", "w") as _file:
    for pair in stapairs:
        sta1 = pair[0]
        sta2 = pair[1]

        x1 = sta_dic[sta1]["x"]
        z1 = sta_dic[sta1]["z"]

        x2 = sta_dic[sta2]["x"]
        z2 = sta_dic[sta2]["z"]

        dis = np.sqrt((x2-x1)**2 + (z2-z1)**2)
        az = az_two_points(x1, z1, x2, z2)

        print(dis/1000.0)

        _file.write(f"{sta1}_{sta2}  {dis}  {az}\n")
