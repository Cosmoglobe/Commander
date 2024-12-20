"""
Script to take the previously generated sky spectra (sky.npy) and plot a map with them.
"""

import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from my_utils import residuals
from scipy.optimize import minimize
from utils.config import gen_nyquistl

reference_maps = "/mn/stornext/d16/cmbco/ola/firas/healpix_maps/"

sky = np.load("../../output/data/sky.npy")
data = h5py.File(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/df_v14.h5",
    "r",
)

idx = np.where(np.array(data["df_data/xcal_pos"][()]) == 2)
pix_gal = np.array(data["df_data/pix_gal"][()])[idx].astype(int)

# frequency mapping
fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

scan_mode = 0  # SS
channel = 3  # LL
frec = 4 * (channel % 2) + scan_mode

f_icm = np.arange(len(sky[0])) * (fnyq["icm"][frec] / 320) + 1
c = 3e8 * 1e2  # cm/s
f_ghz = f_icm * c * 1e-9

NSIDE = 32
npix = hp.nside2npix(NSIDE)

for freq in range(len(f_ghz)):
    hpxmap = np.zeros(npix)
    data_density = np.zeros(npix)

    for i in range(len(pix_gal)):
        hpxmap[pix_gal[i]] += np.abs(sky[i][freq])
        data_density[pix_gal[i]] += 1

    m = np.zeros(npix)
    mask = data_density == 0
    m[~mask] = hpxmap[~mask] / data_density[~mask]

    # mask the map
    m[mask] = hp.UNSEEN

    # hp.mollview(
    #     m,
    #     coord="G",
    #     title=f"{int(f_ghz[freq]):04d} GHz",
    #     unit="MJy/sr",
    #     min=0,
    #     max=500,
    # )
    # # hp.graticule(coord="G")
    # plt.savefig(f"../../output/maps/sky_map_{int(f_ghz[freq]):04d}.png")
    # plt.close()

print(f"pix_gal: {pix_gal.shape}")

# fit bb to each pixel
sum = np.zeros((npix, len(sky[0])))
sum_density = np.zeros(npix)
for i in range(len(pix_gal)):
    sum[pix_gal[i]] += np.abs(sky[i])
    sum_density[pix_gal[i]] += 1

temps = np.zeros(npix)

sum[sum_density != 0] = (
    sum[sum_density != 0] / sum_density[sum_density != 0][:, np.newaxis]
)
t0 = 2.726

for i in range(len(sum)):
    if sum_density[i] != 0:
        fit = minimize(residuals, t0, args=(f_ghz[1:], sum[i][1:]))
        temps[i] = fit.x[0]
    else:
        temps[i] = hp.UNSEEN

print(temps)

hp.mollview(
    temps,
    coord="G",
    title="Temperature map",
    unit="K",
    min=2.7,
    max=3.0,
)
plt.savefig("../../output/maps/temperature_map.png")
plt.close()
