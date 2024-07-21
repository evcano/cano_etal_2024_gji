import itertools
import matplotlib.pyplot as plt
import numpy as np
import os
from cmcrameri import cm


def cart2polar(x, z):
    r = np.sqrt(x**2 + z**2)
    az = (np.arctan2(z, x) + 2*np.pi) % (2*np.pi)
    az = np.degrees(az)
    az = (450.0 - az) % 360.0
    return r, az

#fontsize = 8
fontsize = 10
figx = 2.6
figy = 3.0

stafile = "./STATIONS"
noisefile = "./files_for_inversions/target_noise_source/complete_mesh/azimuthal_noise.ascii"
sta2plot = ["CAR02","CAR07","CAV09","CAU03","CAV07","CAY12"]

# DONT EDIT BELOW THIS LINE
# read stations
stanames = np.loadtxt(stafile, usecols=(0), dtype=str)
stax = np.loadtxt(stafile, usecols=(2))
staz = np.loadtxt(stafile, usecols=(3))
star, staaz = cart2polar(stax, staz)
staaz = np.deg2rad(staaz)

# read noise intensity
noise_model = np.loadtxt(noisefile)

az = np.deg2rad(noise_model[:,0])
energy = noise_model[:,1] / np.max(noise_model[:,1])

# temp
idx = np.argsort(noise_model[:,0])
tmp1 = noise_model[idx,0]
tmp2 = noise_model[idx,1] / np.max(noise_model[:,1])
print(np.mean(tmp2))
plt.plot(tmp1,tmp2)
plt.show()
plt.close()

energy = np.atleast_2d(energy)
az = az.reshape(1,-1)
radius = np.linspace(star.max()*1.2, star.max()*1.3, energy.shape[0]+1).reshape(-1,1)

# figure
fig, ax = plt.subplots(subplot_kw={"polar": True})
fig.set_size_inches(figx, figy)

# plot stations
for i in range(0, stanames.size):
    if stanames[i] in sta2plot:
        if stanames[i] == "CAV07":
            offset = 200000
        else:
            offset = 0

        ax.text(staaz[i], star[i]*1.1+offset, stanames[i][2:], fontsize=fontsize,
            horizontalalignment="left")
        ax.scatter(staaz[i], star[i], c="magenta")

# plot noise intensity
energy = np.vstack((energy,energy))
cmap = cm.devon.resampled(30)
im = ax.pcolormesh(az, radius, energy, cmap=cmap, shading="gouraud")

# figure settings
ax.set_theta_direction(-1)
ax.set_theta_zero_location("N")
ax.set_rticks([])

ax.tick_params(axis="x", labelsize=fontsize)
ax.tick_params(axis="y", labelsize=fontsize)

cbar = plt.colorbar(im, fraction=0.03, orientation="horizontal")
cbar.set_label("normalized PSD", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

plt.savefig("egf_bias.png", dpi=300, bbox_inches="tight")
plt.show()
plt.close()
