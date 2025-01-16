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

sky = np.load("../../output/data/sky.npz")["sky"]
pix_gal = np.load("../../output/data/sky.npz")["pix_gal"]
# data = h5py.File(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v4.2.h5",
#     "r",
# )

# pix_gal = np.array(data["df_data/pix_gal"]).astype(int)
# # pix_terr = np.array(data["df_data/pix_terr"]).astype(int)
# mtm_length = np.array(data["df_data/mtm_length"][()])
# mtm_speed = np.array(data["df_data/mtm_speed"][()])

print(f"sky data size: {len(sky)}")
print(f"pix_gal data size: {len(pix_gal)}")

# short_filter = mtm_length == 0
# slow_filter = mtm_speed == 0
# pix_gal = pix_gal[short_filter & slow_filter]
# pix_terr = pix_terr[short_filter & slow_filter]
# sky = sky[short_filter & slow_filter]


# to not use ICAL higher than 3
# a_ical = np.array(data["df_data/a_ical"][()])
# b_ical = np.array(data["df_data/b_ical"][()])
# ical = (a_ical + b_ical) / 2
# ical = ical[short_filter & slow_filter]
# ical_lower_3 = ical < 3
# # pix_gal = pix_gal[ical_lower_3]
# pix_terr = pix_terr[ical_lower_3]
# sky = sky[ical_lower_3]

# remove data that i flagged
# stat_word_1 = np.array(data["df_data/stat_word_1"][()]).astype(int)
# stat_word_12 = np.array(data["df_data/stat_word_12"][()]).astype(int)
# stat_word_9 = np.array(data["df_data/stat_word_9"][()]).astype(int)
# lvdt_stat_b = np.array(data["df_data/lvdt_stat_b"][()]).astype(int)
# stat_word_13 = np.array(data["df_data/stat_word_13"][()]).astype(int)
# stat_word_16 = np.array(data["df_data/stat_word_16"][()]).astype(int)

# stat_word_1 = stat_word_1[short_filter & slow_filter]
# stat_word_12 = stat_word_12[short_filter & slow_filter]
# stat_word_9 = stat_word_9[short_filter & slow_filter]
# lvdt_stat_b = lvdt_stat_b[short_filter & slow_filter]
# stat_word_13 = stat_word_13[short_filter & slow_filter]
# stat_word_16 = stat_word_16[short_filter & slow_filter]

# filter1 = stat_word_1 != 46
# filter2 = stat_word_12 != 19121
# filter3 = stat_word_9 != 45110
# filter4 = stat_word_12 != 18536
# filter5 = stat_word_12 != 54906
# filter6 = lvdt_stat_b != 83
# filter7 = stat_word_12 != 63675
# filter8 = stat_word_13 != 19585
# filter9 = stat_word_16 != 14372

# # pix_gal = pix_gal[
# #     filter1
# #     & filter2
# #     & filter3
# #     & filter4
# #     & filter5
# #     & filter6
# #     & filter7
# #     & filter8
# #     & filter9
# # ]
# pix_terr = pix_terr[
#     filter1
#     & filter2
#     & filter3
#     & filter4
#     & filter5
#     & filter6
#     & filter7
#     & filter8
#     & filter9
# ]
# sky = sky[
#     filter1
#     & filter2
#     & filter3
#     & filter4
#     & filter5
#     & filter6
#     & filter7
#     & filter8
#     & filter9
# ]

print("plotting sky")

# plot the sky
plt.imshow(
    np.abs(sky).T,
    aspect="auto",
    extent=[0, len(sky), 0, len(sky[0])],
    vmax=500,
    vmin=0,
)
plt.savefig("../../output/plots/sky_over_time.png")
plt.clf()

# frequency mapping
fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

scan_mode = 0  # SS
channel = 3  # LL
frec = 4 * (channel % 2) + scan_mode

f_icm = np.arange(len(sky[0])) * (fnyq["icm"][frec] / 320)
c = 3e8 * 1e2  # cm/s
f_ghz = (
    f_icm * c * 1e-9 + 55
)  # this might not be right but it is what matches the initial frequencies of the firas movie


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
        # norm="hist",
        min=300,
        max=400,
    )
    # hp.graticule(coord="G")
    plt.savefig(f"../../output/maps/sky_map/{int(f_ghz[freq]):04d}.png")
    plt.close()

print("Calculating temperature map")

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
    coord="G",
    title="Temperature map",
    unit="K",
    # norm="hist",
    min=2.7,
    max=2.8,
)
plt.savefig("../../output/maps/temperature_map.png")
plt.close()
