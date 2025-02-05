"""
Script to take the previously generated sky spectra (sky.npy) and plot a map with them.
"""

import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from my_utils import planck, residuals
from scipy.optimize import minimize
from utils.config import gen_nyquistl

T_CMB = 2.72548  # Fixsen 2009
modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}

data = np.load("../../output/data/sky.npz")

# print(data.files)

sky = {}
pix_gal = {}
scan = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            sky[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]
            pix_gal[mode] = data[f"pix_gal_{mode}"]
            scan[f"{channel}_{mode}"] = data[f"scan_{channel}_{mode}"]

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
# fnyq = gen_nyquistl(
#     "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
# )

# c = 3e8 * 1e2  # cm/s
# frec = {}
# f_icm = {}
# f_ghz = {}
# for channel, channel_value in channels.items():
#     for mode, mode_value in modes.items():
#         if mode == "lf" and (channel == "lh" or channel == "rh"):
#             continue
#         else:
#             frec[f"{channel}_{mode}"] = 4 * (channel_value % 2) + channel_value

#             f_icm[f"{channel}_{mode}"] = np.arange(len(sky[f"{channel}_{mode}"][0])) * (
#                 fnyq["icm"][frec[f"{channel}_{mode}"]] / 320
#             )
#             f_ghz[f"{channel}_{mode}"] = (
#                 f_icm[f"{channel}_{mode}"] * c * 1e-9 + 55
#             )  # this might not be right but it is what matches the initial frequencies of the firas movie

#             print(f"{channel}_{mode}: {f_ghz[f"{channel}_{mode}"]}")

nu0 = {"ss": 68.020812, "lf": 23.807283}
dnu = {"ss": 13.604162, "lf": 3.4010405}
nf = {"lh_ss": 210, "ll_lf": 182, "ll_ss": 43, "rh_ss": 210, "rl_lf": 182, "rl_ss": 43}

f_ghz = {}
for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            f_ghz[f"{channel}_{mode}"] = np.linspace(
                nu0[mode],
                nu0[mode] + dnu[mode] * (nf[f"{channel}_{mode}"] - 1),
                nf[f"{channel}_{mode}"],
            )
            # print(f"{channel}_{mode}: {f_ghz[f"{channel}_{mode}"]}")

print("plotting sky")

for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            # plot the sky
            print(f"plotting {channel}_{mode}")
            plt.imshow(
                np.abs(sky[f"{channel}_{mode}"]).T,
                aspect="auto",
                extent=[
                    0,
                    len(sky[f"{channel}_{mode}"]),
                    0,
                    len(sky[f"{channel}_{mode}"][0]),
                ],
                vmax=400,
                vmin=0,
            )
            plt.title(f"{channel}_{mode}")
            plt.savefig(f"../../output/plots/sky_over_time/{f"{channel}_{mode}"}.png")
            plt.clf()


NSIDE = 32
npix = hp.nside2npix(NSIDE)

# for freq in range(len(f_ghz)):
hpxmap = {}
data_density = {}

for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            hpxmap[f"{channel}_{mode}"] = np.zeros(
                (npix, len(f_ghz[f"{channel}_{mode}"]))
            )
            data_density[f"{channel}_{mode}"] = np.zeros(npix)
            for i in range(len(pix_gal[mode])):
                hpxmap[f"{channel}_{mode}"][pix_gal[mode][i]] += np.abs(
                    sky[f"{channel}_{mode}"][i]
                )
                data_density[f"{channel}_{mode}"][pix_gal[mode][i]] += 1
# for i in range(len(pix_terr)):
#     hpxmap[pix_terr[i]] += np.abs(sky[i])
#     data_density[pix_terr[i]] += 1

# monopole


m = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            m[f"{channel}_{mode}"] = np.zeros((npix, len(f_ghz[f"{channel}_{mode}"])))
            monopole = planck(f_ghz[f"{channel}_{mode}"], np.array(T_CMB))
            mask = data_density[f"{channel}_{mode}"] == 0
            # print(
            #     f"shapes of map and density: {hpxmap[f'{channel}_{mode}'].shape}, {data_density[f'{channel}_{mode}'].shape}, monopole {monopole.shape}"
            # )
            m[f"{channel}_{mode}"][~mask] = (
                hpxmap[f"{channel}_{mode}"][~mask]
                / data_density[f"{channel}_{mode}"][~mask, np.newaxis]
            ) - monopole

            # mask the map
            m[f"{channel}_{mode}"][mask] = np.nan  # hp.UNSEEN

for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            for freq in range(len(f_ghz[f"{channel}_{mode}"])):
                hp.mollview(
                    m[f"{channel}_{mode}"][:, freq],
                    # coord="G",
                    title=f"{int(f_ghz[f"{channel}_{mode}"][freq]):04d} GHz as seen by {channel.upper()}{mode.upper()}",
                    unit="MJy/sr",
                    # norm="hist",
                    min=0,
                    max=50,
                )
                # hp.graticule(coord="G")
                plt.savefig(
                    f"../../output/maps/{f"{channel}_{mode}"}/{int(f_ghz[f"{channel}_{mode}"][freq]):04d}.png"
                )
                plt.close()

print("plotting joint map")

# make maps taking into account both of the modes matching the frequencies
# low frequencies
joint_map = hpxmap["ll_lf"] + hpxmap["rl_lf"]
joint_density = data_density["ll_lf"] + data_density["rl_lf"]
for i in range(len(f_ghz["ll_ss"])):
    joint_map[:, i + 3 * i] += (
        hpxmap["ll_ss"][:, i]
        + hpxmap["rl_ss"][:, i]
        + hpxmap["lh_ss"][:, i]
        + hpxmap["rh_ss"][:, i]
    )
    joint_density += (
        data_density["ll_ss"]
        + data_density["rl_ss"]
        + data_density["lh_ss"]
        + data_density["rh_ss"]
    )

m_joint = np.zeros((npix, len(f_ghz["ll_lf"])))
monopole = planck(f_ghz["ll_lf"], np.array(T_CMB))
joint_mask = joint_density == 0
m_joint[~joint_mask] = (
    joint_map[~joint_mask] / joint_density[~joint_mask, np.newaxis] - monopole
)
m_joint[joint_mask] = np.nan  # hp.UNSEEN

for freq in range(len(f_ghz["ll_lf"])):
    hp.mollview(
        m_joint[:, freq],
        # coord="G",
        title=f"{int(f_ghz['ll_lf'][freq]):04d} GHz",
        unit="MJy/sr",
        # norm="hist",
        min=0,
        max=200,
    )
    # hp.graticule(coord="G")
    plt.savefig(f"../../output/maps/joint/{int(f_ghz['ll_lf'][freq]):04d}.png")
    plt.close()

# high frequencies
joint_map = hpxmap["lh_ss"] + hpxmap["rh_ss"]
joint_density = data_density["lh_ss"] + data_density["rh_ss"]
# print(f"joint density shape: {joint_density.shape}")
m_joint = np.zeros((npix, (len(f_ghz["lh_ss"]) - len(f_ghz["ll_ss"]))))
monopole = planck(f_ghz["lh_ss"], np.array(T_CMB))[len(f_ghz["ll_ss"]) :]
joint_mask = joint_density == 0
# print(
#     f"shapes of joint map and joint density: {joint_map.shape}, {joint_density.shape} and monopole {monopole.shape}"
# )
m_joint[~joint_mask] = (
    joint_map[~joint_mask] / joint_density[~joint_mask, np.newaxis]
)[:, len(f_ghz["ll_ss"]) :] - monopole
m_joint[joint_mask] = np.nan  # hp.UNSEEN

# print(f"shape of m_joint after mask: {m_joint.shape}")

for freq in range(len(f_ghz["ll_ss"]), len(f_ghz["lh_ss"])):
    # print(f"plotting {freq}: {int(f_ghz['lh_ss'][freq])}")
    hp.mollview(
        m_joint[:, (freq - len(f_ghz["ll_ss"]))],
        # coord="G",
        title=f"{int(f_ghz['lh_ss'][freq]):04d} GHz",
        unit="MJy/sr",
        # norm="hist",
        min=0,
        max=50,
    )
    # hp.graticule(coord="G")
    plt.savefig(f"../../output/maps/joint/{int(f_ghz['lh_ss'][freq]):04d}.png")
    plt.close()

print("plotting up/down scan map")

for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            hpxmap_up = np.zeros((npix, len(f_ghz[f"{channel}_{mode}"])))
            hpxmap_down = np.zeros((npix, len(f_ghz[f"{channel}_{mode}"])))
            data_density_up = np.zeros(npix)
            data_density_down = np.zeros(npix)
            for i in range(len(pix_gal[mode])):
                if scan[f"{channel}_{mode}"][i] == 1:
                    hpxmap_up[pix_gal[mode][i]] += np.abs(sky[f"{channel}_{mode}"][i])
                    data_density_up[pix_gal[mode][i]] += 1
                else:
                    hpxmap_down[pix_gal[mode][i]] += np.abs(sky[f"{channel}_{mode}"][i])
                    data_density_down[pix_gal[mode][i]] += 1

            m_up = np.zeros((npix, len(f_ghz[f"{channel}_{mode}"])))
            m_down = np.zeros((npix, len(f_ghz[f"{channel}_{mode}"])))
            monopole = planck(f_ghz[f"{channel}_{mode}"], np.array(T_CMB))
            mask_up = data_density_up == 0
            mask_down = data_density_down == 0
            m_up[~mask_up] = (
                hpxmap_up[~mask_up] / data_density_up[~mask_up, np.newaxis]
            ) - monopole
            m_down[~mask_down] = (
                hpxmap_down[~mask_down] / data_density_down[~mask_down, np.newaxis]
            ) - monopole
            m_up[mask_up] = np.nan  # hp.UNSEEN
            m_down[mask_down] = np.nan

            for freq in range(len(f_ghz[f"{channel}_{mode}"])):
                hp.mollview(
                    m_up[:, freq],
                    title=f"{int(f_ghz[f"{channel}_{mode}"][freq]):04d} GHz as seen by {channel.upper()}{mode.upper()} in up scan",
                    unit="MJy/sr",
                    min=0,
                    max=50,
                )
                plt.savefig(
                    f"../../output/maps/up_down_scan/{f"{channel}_{mode}"}/{int(f_ghz[f"{channel}_{mode}"][freq]):04d}_up.png"
                )
                plt.close()

                hp.mollview(
                    m_down[:, freq],
                    title=f"{int(f_ghz[f"{channel}_{mode}"][freq]):04d} GHz as seen by {channel.upper()}{mode.upper()} in down scan",
                    unit="MJy/sr",
                    min=0,
                    max=50,
                )
                # hp.graticule(coord="G")
                plt.savefig(
                    f"../../output/maps/up_down_scan/{f"{channel}_{mode}"}/{int(f_ghz[f"{channel}_{mode}"][freq]):04d}_down.png"
                )
                plt.close()

# print("Calculating temperature map")


# t0 = 2.726

# fit = {}
# temps = {}
# for channel in ["ll", "rl"]:
#     for mode in modes:
#         temps[f"{channel}_{mode}"] = np.zeros(npix)
#         for i in range(len(m[f"{channel}_{mode}"])):
#             fit[f"{channel}_{mode}"] = minimize(
#                 residuals,
#                 t0,
#                 args=(
#                     f_ghz[f"{channel}_{mode}"],
#                     m[f"{channel}_{mode}"][i],
#                 ),
#             )
#             temps[f"{channel}_{mode}"][i] = fit[f"{channel}_{mode}"].x[0]

#         print(f"{f"{channel}_{mode}"}: {temps[f"{channel}_{mode}"]}")

#         hp.mollview(
#             temps[f"{channel}_{mode}"],
#             coord="G",
#             title="Temperature map",
#             unit="K",
#             # norm="hist",
#             min=2.7,
#             max=2.8,
#         )
#         plt.savefig(f"../../output/maps/temperature_map_{f"{channel}_{mode}"}.png")
#         plt.close()

# temps = np.zeros(npix)
# for i in range(len(m_joint)):
#     if joint_density[i] != 0:
#         fit = minimize(residuals, t0, args=(f_ghz["ll_lf"], m_joint[i]))
#         temps[i] = fit.x[0]
#     else:
#         temps[i] = np.nan  # hp.UNSEEN

# print(temps)

# hp.mollview(
#     temps,
#     coord="G",
#     title="Temperature map",
#     unit="K",
#     # norm="hist",
#     min=2.7,
#     max=2.8,
# )
# plt.savefig("../../output/maps/temperature_map.png")
# plt.close()
