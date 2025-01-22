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
modes = {"ss": 0, "lf": 3}

data = np.load("../../output/data/sky.npz")

sky = {}
pix_gal = {}
for mode in modes.keys():
    sky[mode] = data[mode]
    pix_gal[mode] = data[f"pix_gal_{mode}"]

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

# frequency mapping
fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

channel = 3  # LL

c = 3e8 * 1e2  # cm/s
frec = {}
f_icm = {}
f_ghz = {}
for key, value in modes.items():
    frec[key] = 4 * (channel % 2) + value

    f_icm[key] = np.arange(len(sky[key][0])) * (fnyq["icm"][frec[key]] / 320)
    f_ghz[key] = (
        f_icm[key] * c * 1e-9 + 55
    )  # this might not be right but it is what matches the initial frequencies of the firas movie

print("plotting sky")

for mode in modes.keys():
    # plot the sky
    plt.imshow(
        np.abs(sky[mode]).T,
        aspect="auto",
        extent=[0, len(sky[mode]), 0, len(sky[mode][0])],
        vmax=400,
        vmin=0,
    )
    plt.savefig(f"../../output/plots/sky_over_time_{mode}.png")
    plt.clf()


NSIDE = 32
npix = hp.nside2npix(NSIDE)

# for freq in range(len(f_ghz)):
hpxmap = {}
data_density = {}

for mode in modes.keys():
    hpxmap[mode] = np.zeros((npix, len(f_ghz[mode])))
    data_density[mode] = np.zeros(npix)
    for i in range(len(pix_gal[mode])):
        hpxmap[mode][pix_gal[mode][i]] += np.abs(sky[mode][i])
        data_density[mode][pix_gal[mode][i]] += 1
# for i in range(len(pix_terr)):
#     hpxmap[pix_terr[i]] += np.abs(sky[i])
#     data_density[pix_terr[i]] += 1

m = {}
for mode in modes.keys():
    m[mode] = np.zeros((npix, len(f_ghz[mode])))
    mask = data_density[mode] == 0
    m[mode][~mask] = hpxmap[mode][~mask] / data_density[mode][~mask, np.newaxis]

    # mask the map
    m[mode][mask] = hp.UNSEEN

for mode in modes.keys():
    for freq in range(len(f_ghz[mode])):
        hp.mollview(
            m[mode][:, freq],
            # coord="G",
            title=f"{int(f_ghz[mode][freq]):04d} GHz",
            unit="MJy/sr",
            # norm="log",
            min=0,
            max=400,
        )
        # hp.graticule(coord="G")
        plt.savefig(
            f"../../output/maps/sky_map/{mode}/{int(f_ghz[mode][freq]):04d}.png"
        )
        plt.close()

# make maps taking into account both of the modes matching the frequencies
joint_map = hpxmap["lf"]
joint_density = data_density["lf"]
for i in range(len(f_ghz["ss"])):
    joint_map[:, i + 3 * i] += hpxmap["ss"][:, i]
    joint_density += data_density["ss"]

m_joint = np.zeros((npix, len(f_ghz["lf"])))
joint_mask = joint_density == 0
m_joint[~joint_mask] = joint_map[~joint_mask] / joint_density[~joint_mask, np.newaxis]
m_joint[joint_mask] = hp.UNSEEN

for freq in range(len(f_ghz["lf"])):
    hp.mollview(
        m_joint[:, freq],
        # coord="G",
        title=f"{int(f_ghz['lf'][freq]):04d} GHz",
        unit="MJy/sr",
        # norm="hist",
        min=100,
        max=400,
    )
    # hp.graticule(coord="G")
    plt.savefig(f"../../output/maps/sky_map/joint/{int(f_ghz['lf'][freq]):04d}.png")
    plt.close()

print("Calculating temperature map")


t0 = 2.726

fit = {}
temps = {}
for mode in modes:
    temps[mode] = np.zeros(npix)
    for i in range(len(m[mode])):
        fit[mode] = minimize(residuals, t0, args=(f_ghz[mode][1:], m[mode][i, 1:]))
        temps[mode] = fit[mode].x[0]

    print(f"{mode}: {temps}")

    hp.mollview(
        temps[mode],
        coord="G",
        title="Temperature map",
        unit="K",
        # norm="hist",
        min=2.7,
        max=2.8,
    )
    plt.savefig(f"../../output/maps/temperature_map_{mode}.png")
    plt.close()

for i in range(len(m_joint)):
    if joint_density[i] != 0:
        fit = minimize(residuals, t0, args=(f_ghz["lf"][1:], m_joint[i, 1:]))
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
