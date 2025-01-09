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

T_CMB = 2.72548  # Fixsen 2009

reference_maps = "/mn/stornext/d16/cmbco/ola/firas/healpix_maps/"

sky = np.load("../../output/data/sky.npy")
data = h5py.File(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v4.1.h5",
    "r",
)

pix_gal = np.array(data["df_data/pix_gal"]).astype(int)
# pix_terr = np.array(data["df_data/pix_terr"]).astype(int)

# to not use ICAL higher than 3
a_ical = np.array(data["df_data/a_ical"][()])
b_ical = np.array(data["df_data/b_ical"][()])
ical = (a_ical + b_ical) / 2
ical_lower_3 = ical < 3
pix_gal = pix_gal[ical_lower_3]
# pix_terr = pix_terr[ical_lower_3]
sky = sky[ical_lower_3]

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

# for freq in range(len(f_ghz)):
hpxmap = np.zeros((npix, len(f_ghz)))
data_density = np.zeros(npix)

for i in range(len(pix_gal)):
    hpxmap[pix_gal[i]] += np.abs(sky[i])
    data_density[pix_gal[i]] += 1
# for i in range(len(pix_terr)):
#     hpxmap[pix_terr[i]] += np.abs(sky[i])
#     data_density[pix_terr[i]] += 1

m = np.zeros((npix, len(f_ghz)))
mask = data_density == 0
m[~mask] = hpxmap[~mask] / data_density[~mask, np.newaxis]

# mask the map
m[mask] = hp.UNSEEN

for freq in range(len(f_ghz)):
    hp.mollview(
        m[:, freq],
        # coord="G",
        title=f"{int(f_ghz[freq]):04d} GHz",
        unit="MJy/sr",
        norm="hist",
        # min=0,
        # max=500,
    )
    # hp.graticule(coord="G")
    plt.savefig(f"../../output/maps/sky_map/{int(f_ghz[freq]):04d}.png")
    plt.close()

temps = np.zeros(npix)

t0 = 2.726

for i in range(len(m)):
    if data_density[i] != 0:
        fit = minimize(residuals, t0, args=(f_ghz[1:], m[i, 1:]))
        temps[i] = fit.x[0]
    else:
        temps[i] = hp.UNSEEN

print(temps)

hp.mollview(
    temps,
    # coord="G",
    title="Temperature map",
    unit="K",
    norm="hist",
    # min=2.7,
    # max=2.8,
)
plt.savefig("../../output/maps/temperature_map.png")
plt.close()
