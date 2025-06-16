"""
Script to take the previously generated sky spectra (sky.npy) and plot a map with them.
"""

import os
import sys
from pathlib import Path

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from globals import (CHANNELS, COORDINATES, FITS, JOINT, MODES, NSIDE, OFFSET,
                     PNG, PROCESSED_DATA_PATH, SCANUPDOWN, T_CMB)
from healpy.rotator import Rotator

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from my_utils import generate_frequencies, planck

user = os.environ["USER"]
save_path = f"/mn/stornext/d16/www_cmb/{user}/firas/"
fits_path = "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/fits_files/"

modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}

data = np.load(PROCESSED_DATA_PATH)
# print(data.files)

sky = {}
scan = {}
if COORDINATES == "G":
    pix_gal = {}
    folder = "galactic"

elif COORDINATES == "E":
    pix_ecl = {}
    folder = "ecliptic"


for channel in channels.keys():
    if channel in CHANNELS:
        for mode in modes.keys():
            if mode in MODES:
                if not (mode == "lf" and (channel == "lh" or channel == "rh")):
                    sky[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]
                    scan[mode] = data[f"scan_{mode}"]
                    if COORDINATES == "G":
                        if NSIDE == 32:
                            pix_gal[mode] = data[f"pix_gal_{mode}"]
                        else:
                            gal_lat = data[f"gal_lat_{mode}"]
                            gal_lon = data[f"gal_lon_{mode}"]
                            pix_gal[mode] = hp.ang2pix(NSIDE, gal_lon, gal_lat, lonlat=True).astype(int)
                    elif COORDINATES == "E":
                        if NSIDE == 32:
                            pix_ecl[mode] = data[f"pix_ecl_{mode}"]
                        else:
                            # rotate gal lat and lon to ecliptic
                            r = Rotator(coord=['G','E'])
                            gal_lat = data[f"gal_lat_{mode}"]
                            gal_lon = data[f"gal_lon_{mode}"]
                            ecl_lat, ecl_lon = r(gal_lat, gal_lon)
                            pix_ecl[mode] = hp.ang2pix(NSIDE, ecl_lon, ecl_lat, lonlat=True)

f_ghz = {}
for channel in channels.keys():
    if channel in CHANNELS:
        for mode in modes.keys():
            if mode in MODES:
                if not(mode == "lf" and (channel == "lh" or channel == "rh")):
                    f_ghz[f"{channel}_{mode}"] = generate_frequencies(channel, mode)

print("plotting sky")

for channel in channels.keys():
    if channel in CHANNELS:
        for mode in modes.keys():
            if mode in MODES:
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
                    plt.savefig(f"{save_path}plots/sky_over_time/{f"{channel}_{mode}"}.png")
                    plt.clf()

npix = hp.nside2npix(NSIDE)

# for freq in range(len(f_ghz)):
hpxmap = {}
data_density = {}

curr_path = f"{fits_path}maps/hit_maps/{folder}/"
# check if curr_path exists and create it if it doesn't
if not os.path.exists(curr_path):
    os.makedirs(curr_path)

for channel in channels.keys():
    if channel in CHANNELS:
        for mode in modes.keys():
            if mode in MODES:
                if not (mode == "lf" and (channel == "lh" or channel == "rh")):
                    hpxmap[f"{channel}_{mode}"] = np.zeros(
                        (npix, len(f_ghz[f"{channel}_{mode}"]))
                    )
                    data_density[f"{channel}_{mode}"] = np.zeros(npix)

                    if COORDINATES == "G":
                        for i in range(len(pix_gal[mode])):
                            hpxmap[f"{channel}_{mode}"][pix_gal[mode][i]] += np.abs(
                                sky[f"{channel}_{mode}"][i]
                            )
                            data_density[f"{channel}_{mode}"][pix_gal[mode][i]] += 1
                    elif COORDINATES == "E":
                        for i in range(len(pix_ecl[mode])):
                            hpxmap[f"{channel}_{mode}"][pix_ecl[mode][i]] += np.abs(
                                sky[f"{channel}_{mode}"][i]
                            )
                            data_density[f"{channel}_{mode}"][pix_ecl[mode][i]] += 1

                    # plot hit map
                    if PNG:
                        hp.mollview(
                            data_density[f"{channel}_{mode}"],
                            title=f"Hit map for {channel.upper()}{mode.upper()}",
                            unit="hits",
                            norm="hist",
                        )
                        Path.mkdir(
                            Path(f"{save_path}maps/hit_maps/{folder}"), parents=True, exist_ok=True
                        )
                        plt.savefig(
                            f"{save_path}maps/hit_maps/{folder}/{f"{channel}_{mode}_offset{OFFSET}_nside{NSIDE}"}.png"
                        )
                        plt.close()
                    if FITS:
                        # hp.write_map(
                        #     f"{curr_path}{f"{channel}_{mode}_offset{OFFSET}_nside{NSIDE}"}.fits",
                        #     data_density[f"{channel}_{mode}"],
                        #     overwrite=True,
                        # )
                        fits.writeto(f"{curr_path}{f"{channel}_{mode}_offset{OFFSET}_nside{NSIDE}"}.fits", data_density[f"{channel}_{mode}"], overwrite=True)
# for i in range(len(pix_terr)):
#     hpxmap[pix_terr[i]] += np.abs(sky[i])
#     data_density[pix_terr[i]] += 1

# monopole


m = {}
for channel in channels.keys():
    if channel in CHANNELS:
        for mode in modes.keys():
            if mode in MODES:
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
    if channel in CHANNELS:
        for mode in modes.keys():
            if mode in MODES:
                Path(f"{save_path}maps/frequency_maps/{f"{channel}_{mode}"}/{folder}").mkdir(
                    parents=True, exist_ok=True
                )
                if not (mode == "lf" and (channel == "lh" or channel == "rh")):
                    for freq in range(len(f_ghz[f"{channel}_{mode}"])):
                        if PNG:
                            max_freq = 200
                            if channel == "ll" and mode == "ss":
                                max_freq = 50
                            hp.mollview(
                                m[f"{channel}_{mode}"][:, freq],
                                # coord="G",
                                title=f"{int(f_ghz[f"{channel}_{mode}"][freq]):04d} GHz as seen by {channel.upper()}{mode.upper()}",
                                unit="MJy/sr",
                                # norm="hist",
                                min=0,
                                max=max_freq,
                            )
                            # hp.graticule(coord="G")
                            plt.savefig(
                                f"{save_path}maps/frequency_maps/{f"{channel}_{mode}"}/{folder}/{int(f_ghz[f"{channel}_{mode}"][freq]):04d}_offset{OFFSET}_nside{NSIDE}.png"
                            )
                            plt.close()
                        if FITS:
                            curr_path = f"{fits_path}maps/frequency_maps/{f"{channel}_{mode}"}/{folder}/"
                            if not os.path.exists(curr_path):
                                os.makedirs(curr_path)
                            # hp.write_map(
                            #     f"{curr_path}{int(f_ghz[f'{channel}_{mode}'][freq]):04d}_offset{OFFSET}_nside{NSIDE}.fits",
                            #     m[f"{channel}_{mode}"][:, freq],
                            #     overwrite=True,
                            # )
                            fits.writeto(f"{curr_path}{int(f_ghz[f'{channel}_{mode}'][freq]):04d}_offset{OFFSET}_nside{NSIDE}.fits", m[f"{channel}_{mode}"][:, freq], overwrite=True)


if JOINT:
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

    Path(f"{save_path}maps/joint/{folder}").mkdir(parents=True, exist_ok=True)
    for freq in range(len(f_ghz["ll_lf"])):
        if PNG:
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
            plt.savefig(
                f"{save_path}maps/joint/{folder}/{int(f_ghz['ll_lf'][freq]):04d}_offset{OFFSET}_nside{NSIDE}.png"
            )
            plt.close()
        if FITS:
            curr_path = f"{fits_path}maps/joint/{folder}/"
            if not os.path.exists(curr_path):
                os.makedirs(curr_path)
            # hp.write_map(
            #     f"{curr_path}{int(f_ghz['ll_lf'][freq]):04d}_offset{OFFSET}_nside{NSIDE}.fits",
            #     m_joint[:, freq],
            #     overwrite=True,
            # )
            fits.writeto(f"{curr_path}{int(f_ghz['ll_lf'][freq]):04d}_offset{OFFSET}_nside{NSIDE}.fits", m_joint[:, freq], overwrite=True)

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
    curr_path = f"{fits_path}maps/joint/{folder}/"
    if not os.path.exists(curr_path):
        os.makedirs(curr_path)
    for freq in range(len(f_ghz["ll_ss"]), len(f_ghz["lh_ss"])):
        if PNG:
            # print(f"plotting {freq}: {int(f_ghz['lh_ss'][freq])}")
            hp.mollview(
                m_joint[:, (freq - len(f_ghz["ll_ss"]))],
                # coord="G",
                title=f"{int(f_ghz['lh_ss'][freq]):04d} GHz",
                unit="MJy/sr",
                # norm="hist",
                min=0,
                max=200,
            )
            # hp.graticule(coord="G")
            plt.savefig(
                f"{save_path}maps/joint/{folder}/{int(f_ghz['lh_ss'][freq]):04d}_offset{OFFSET}_nside{NSIDE}.png"
            )
            plt.close()
        if FITS:
            # hp.write_map(
            #     f"{curr_path}{int(f_ghz['lh_ss'][freq]):04d}_offset{OFFSET}_nside{NSIDE}.fits",
            #     m_joint[:, (freq - len(f_ghz["ll_ss"]))],
            #     overwrite=True,
            # )
            fits.writeto(f"{curr_path}{int(f_ghz['lh_ss'][freq]):04d}_offset{OFFSET}_nside{NSIDE}.fits", m_joint[:, (freq - len(f_ghz["ll_ss"]))], overwrite=True)

if SCANUPDOWN:

    print("plotting up/down scan map")

    for channel in channels.keys():
        for mode in modes.keys():
            if not (mode == "lf" and (channel == "lh" or channel == "rh")):
                Path(f"{save_path}maps/up_down_scan/{channel}_{mode}/{folder}").mkdir(
                    parents=True, exist_ok=True
                )
                hpxmap_up = np.zeros((npix, len(f_ghz[f"{channel}_{mode}"])))
                hpxmap_down = np.zeros((npix, len(f_ghz[f"{channel}_{mode}"])))
                data_density_up = np.zeros(npix)
                data_density_down = np.zeros(npix)

                if COORDINATES == "G":
                    for i in range(len(pix_gal[mode])):
                        if scan[mode][i] == 1:
                            hpxmap_up[pix_gal[mode][i]] += np.abs(
                                sky[f"{channel}_{mode}"][i]
                            )
                            data_density_up[pix_gal[mode][i]] += 1
                        elif scan[mode][i] == 0:
                            hpxmap_down[pix_gal[mode][i]] += np.abs(
                                sky[f"{channel}_{mode}"][i]
                            )
                            data_density_down[pix_gal[mode][i]] += 1
                elif COORDINATES == "E":
                    for i in range(len(pix_ecl[mode])):
                        if scan[mode][i] == 1:
                            hpxmap_up[pix_ecl[mode][i]] += np.abs(
                                sky[f"{channel}_{mode}"][i]
                            )
                            data_density_up[pix_ecl[mode][i]] += 1
                        elif scan[mode][i] == 0:
                            hpxmap_down[pix_ecl[mode][i]] += np.abs(
                                sky[f"{channel}_{mode}"][i]
                            )
                            data_density_down[pix_ecl[mode][i]] += 1

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
                    if PNG:
                        hp.mollview(
                            m_up[:, freq],
                            title=f"{int(f_ghz[f"{channel}_{mode}"][freq]):04d} GHz as seen by {channel.upper()}{mode.upper()} in up scan",
                            unit="MJy/sr",
                            min=0,
                            max=200,
                        )
                        plt.savefig(
                            f"{save_path}maps/up_down_scan/{f"{channel}_{mode}"}/{folder}/{int(f_ghz[f"{channel}_{mode}"][freq]):04d}_up_offset{OFFSET}_nside{NSIDE}.png"
                        )
                        plt.close()

                        hp.mollview(
                            m_down[:, freq],
                            title=f"{int(f_ghz[f"{channel}_{mode}"][freq]):04d} GHz as seen by {channel.upper()}{mode.upper()} in down scan",
                            unit="MJy/sr",
                            min=0,
                            max=200,
                        )
                        # hp.graticule(coord="G")
                        plt.savefig(
                            f"{save_path}maps/up_down_scan/{f"{channel}_{mode}"}/{folder}/{int(f_ghz[f"{channel}_{mode}"][freq]):04d}_down_offset{OFFSET}_nside{NSIDE}.png"
                        )
                        plt.close()
                    if FITS:
                        curr_path = f"{fits_path}maps/up_down_scan/{f"{channel}_{mode}"}/{folder}/"
                        if not os.path.exists(curr_path):
                            os.makedirs(curr_path)
                        # hp.write_map(
                        #     f"{curr_path}{int(f_ghz[f'{channel}_{mode}'][freq]):04d}_up_offset{OFFSET}_nside{NSIDE}.fits",
                        #     m_up[:, freq],
                        #     overwrite=True,
                        # )
                        # hp.write_map(
                        #     f"{curr_path}{int(f_ghz[f'{channel}_{mode}'][freq]):04d}_down_offset{OFFSET}_nside{NSIDE}.fits",
                        #     m_down[:, freq],
                        #     overwrite=True,
                        # )
                        fits.writeto(f"{curr_path}{int(f_ghz[f'{channel}_{mode}'][freq]):04d}_up_offset{OFFSET}_nside{NSIDE}.fits", m_up[:, freq], overwrite=True)
                        fits.writeto(f"{curr_path}{int(f_ghz[f'{channel}_{mode}'][freq]):04d}_down_offset{OFFSET}_nside{NSIDE}.fits", m_down[:, freq], overwrite=True)

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
