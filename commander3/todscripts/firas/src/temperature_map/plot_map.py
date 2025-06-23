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
from healpy.rotator import Rotator

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import globals as g
from my_utils import generate_frequencies, planck

fits_path = "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/fits_files/"

data = np.load(g.PROCESSED_DATA_PATH)

if g.COORDINATES == "G":
    folder = "galactic"

elif g.COORDINATES == "E":
    folder = "ecliptic"

curr_path = f"{fits_path}maps/hit_maps/{folder}/"
# check if curr_path exists and create it if it doesn't
if not os.path.exists(curr_path):
    os.makedirs(curr_path)

hpxmap = {}
data_density = {}
f_ghz = {}
for channel in g.CHANNELS:
    for mode in g.MODES: 
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            sky = data[f"{channel}_{mode}"]
            scan = data[f"scan_{mode}"]
            if g.COORDINATES == "G":
                if g.NSIDE == 32:
                    pix_gal = data[f"pix_gal_{mode}"]
                else:
                    gal_lat = data[f"gal_lat_{mode}"]
                    gal_lon = data[f"gal_lon_{mode}"]
                    pix_gal = hp.ang2pix(g.NSIDE, gal_lon, gal_lat, lonlat=True).astype(int)
            elif g.COORDINATES == "E":
                if g.NSIDE == 32:
                    pix_ecl = data[f"pix_ecl_{mode}"]
                else:
                    # rotate gal lat and lon to ecliptic
                    r = Rotator(coord=['G','E'])
                    gal_lat = data[f"gal_lat_{mode}"]
                    gal_lon = data[f"gal_lon_{mode}"]
                    ecl_lat, ecl_lon = r(gal_lat, gal_lon)
                    pix_ecl = hp.ang2pix(g.NSIDE, ecl_lon, ecl_lat, lonlat=True)

            f_ghz[f"{channel}_{mode}"] = generate_frequencies(channel, mode)

            print(f"Plotting sky for {channel.upper()}{mode.upper()}")
            plt.imshow(
                np.abs(sky).T,
                aspect="auto",
                extent=[
                    0,
                    len(sky),
                    0,
                    len(sky[0]),
                ],
                vmax=400,
                vmin=0,
            )
            plt.title(f"{channel}_{mode}")
            plt.savefig(f"{g.SAVE_PATH}plots/sky_over_time/{f"{channel}_{mode}"}.png")
            plt.clf()

            hpxmap[f"{channel}_{mode}"] = np.zeros(
                (g.NPIX, len(f_ghz[f"{channel}_{mode}"]))
            )
            data_density[f"{channel}_{mode}"] = np.zeros(g.NPIX)

            if g.COORDINATES == "G":
                for i in range(len(pix_gal)):
                    hpxmap[f"{channel}_{mode}"][pix_gal[i]] += np.abs(
                        sky[i]
                    )
                    data_density[f"{channel}_{mode}"][pix_gal[i]] += 1
            elif g.COORDINATES == "E":
                for i in range(len(pix_ecl)):
                    hpxmap[f"{channel}_{mode}"][pix_ecl[i]] += np.abs(
                        sky[i]
                    )
                    data_density[f"{channel}_{mode}"][pix_ecl[i]] += 1

            # plot hit map
            if g.PNG:
                hp.mollview(
                    data_density[f"{channel}_{mode}"],
                    title=f"{channel.upper()}{mode.upper()}",
                    unit="hits",
                    norm="hist",
                )
                Path.mkdir(
                    Path(f"{g.SAVE_PATH}maps/hit_maps/{folder}"), parents=True, exist_ok=True
                )
                plt.savefig(
                    f"{g.SAVE_PATH}maps/hit_maps/{folder}/{f"{channel}_{mode}_nside{g.NSIDE}"}.png"
                )
                plt.close()
            if g.FITS:
                fits.writeto(f"{curr_path}{f"{channel}_{mode}_nside{g.NSIDE}"}.fits", data_density[f"{channel}_{mode}"], overwrite=True)

            m = np.zeros((g.NPIX, len(f_ghz[f"{channel}_{mode}"])))
            monopole = planck(f_ghz[f"{channel}_{mode}"], np.array(g.T_CMB))
            mask = data_density[f"{channel}_{mode}"] == 0
            m[~mask] = (
                hpxmap[f"{channel}_{mode}"][~mask]
                / data_density[f"{channel}_{mode}"][~mask, np.newaxis]
            ) - monopole

            # mask the map
            m[mask] = np.nan  # hp.UNSEEN

            for freq in range(len(f_ghz[f"{channel}_{mode}"])):
                if g.PNG:
                    max_amp = 200
                    if channel == "ll" and mode == "ss":
                        max_amp = 25
                    hp.mollview(
                        m[:, freq],
                        title=f"{int(f_ghz[f"{channel}_{mode}"][freq]):04d} GHz as seen by {channel.upper()}{mode.upper()}",
                        unit="MJy/sr",
                        min=0,
                        max=max_amp,
                    )
                    plt.savefig(
                        f"{g.SAVE_PATH}maps/frequency_maps/{f"{channel}_{mode}"}/{folder}/{int(f_ghz[f"{channel}_{mode}"][freq]):04d}_nside{g.NSIDE}.png"
                    )
                    plt.close()
                if g.FITS:
                    fits.writeto(f"{curr_path}{int(f_ghz[f'{channel}_{mode}'][freq]):04d}_nside{g.NSIDE}.fits", m[:, freq], overwrite=True)

            if g.SCANUPDOWN:
                print("Plotting up/down scan map")
                hpxmap_up = np.zeros((g.NPIX, len(f_ghz[f"{channel}_{mode}"])))
                hpxmap_down = np.zeros((g.NPIX, len(f_ghz[f"{channel}_{mode}"])))
                data_density_up = np.zeros(g.NPIX)
                data_density_down = np.zeros(g.NPIX)

                if g.COORDINATES == "G":
                    for i in range(len(pix_gal)):
                        if scan[i] == 1:
                            hpxmap_up[pix_gal[i]] += np.abs(
                                sky[f"{channel}_{mode}"][i]
                            )
                            data_density_up[pix_gal[i]] += 1
                        elif scan[i] == 0:
                            hpxmap_down[pix_gal[i]] += np.abs(
                                sky[f"{channel}_{mode}"][i]
                            )
                            data_density_down[pix_gal[i]] += 1
                elif g.COORDINATES == "E":
                    for i in range(len(pix_ecl)):
                        if scan[i] == 1:
                            hpxmap_up[pix_ecl[i]] += np.abs(
                                sky[f"{channel}_{mode}"][i]
                            )
                            data_density_up[pix_ecl[i]] += 1
                        elif scan[i] == 0:
                            hpxmap_down[pix_ecl[i]] += np.abs(
                                sky[f"{channel}_{mode}"][i]
                            )
                            data_density_down[pix_ecl[i]] += 1

                m_up = np.zeros((g.NPIX, len(f_ghz[f"{channel}_{mode}"])))
                m_down = np.zeros((g.NPIX, len(f_ghz[f"{channel}_{mode}"])))
                monopole = planck(f_ghz[f"{channel}_{mode}"], np.array(g.T_CMB))
                mask_up = data_density_up == 0
                mask_down = data_density_down == 0
                m_up[~mask_up] = (
                    hpxmap_up[~mask_up] / data_density_up[~mask_up, np.newaxis]
                ) - monopole
                m_down[~mask_down] = (
                    hpxmap_down[~mask_down] / data_density_down[~mask_down, np.newaxis]
                ) - monopole
                m_up[mask_up] = np.nan
                m_down[mask_down] = np.nan

                for freq in range(len(f_ghz[f"{channel}_{mode}"])):
                    if g.PNG:
                        hp.mollview(
                            m_up[:, freq],
                            title=f"{int(f_ghz[f"{channel}_{mode}"][freq]):04d} GHz as seen by {channel.upper()}{mode.upper()} in up scan",
                            unit="MJy/sr",
                            min=0,
                            max=200,
                        )
                        plt.savefig(
                            f"{g.SAVE_PATH}maps/up_down_scan/{f"{channel}_{mode}"}/{folder}/{int(f_ghz[f"{channel}_{mode}"][freq]):04d}_up_nside{g.NSIDE}.png"
                        )
                        plt.close()

                        hp.mollview(
                            m_down[:, freq],
                            title=f"{int(f_ghz[f"{channel}_{mode}"][freq]):04d} GHz as seen by {channel.upper()}{mode.upper()} in down scan",
                            unit="MJy/sr",
                            min=0,
                            max=200,
                        )
                        plt.savefig(
                            f"{g.SAVE_PATH}maps/up_down_scan/{f"{channel}_{mode}"}/{folder}/{int(f_ghz[f"{channel}_{mode}"][freq]):04d}_down_nside{g.NSIDE}.png"
                        )
                        plt.close()
                    if g.FITS:
                        fits.writeto(f"{curr_path}{int(f_ghz[f'{channel}_{mode}'][freq]):04d}_up_nside{g.NSIDE}.fits", m_up[:, freq], overwrite=True)
                        fits.writeto(f"{curr_path}{int(f_ghz[f'{channel}_{mode}'][freq]):04d}_down_nside{g.NSIDE}.fits", m_down[:, freq], overwrite=True)


if g.JOINT:
    print("Plotting joint map")

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

    m_joint = np.zeros((g.NPIX, len(f_ghz["ll_lf"])))
    monopole = planck(f_ghz["ll_lf"], np.array(g.T_CMB))
    joint_mask = joint_density == 0
    m_joint[~joint_mask] = (
        joint_map[~joint_mask] / joint_density[~joint_mask, np.newaxis] - monopole
    )
    m_joint[joint_mask] = np.nan  # hp.UNSEEN

    for freq in range(len(f_ghz["ll_lf"])):
        if g.PNG:
            hp.mollview(
                m_joint[:, freq],
                title=f"{int(f_ghz['ll_lf'][freq]):04d} GHz",
                unit="MJy/sr",
                min=0,
                max=200,
            )
            plt.savefig(
                f"{g.SAVE_PATH}maps/joint/{folder}/{int(f_ghz['ll_lf'][freq]):04d}_nside{g.NSIDE}.png"
            )
            plt.close()
        if g.FITS:
            fits.writeto(f"{curr_path}{int(f_ghz['ll_lf'][freq]):04d}_nside{g.NSIDE}.fits", m_joint[:, freq], overwrite=True)

    # high frequencies
    joint_map = hpxmap["lh_ss"] + hpxmap["rh_ss"]
    joint_density = data_density["lh_ss"] + data_density["rh_ss"]
    m_joint = np.zeros((g.NPIX, (len(f_ghz["lh_ss"]) - len(f_ghz["ll_ss"]))))
    monopole = planck(f_ghz["lh_ss"], np.array(g.T_CMB))[len(f_ghz["ll_ss"]) :]
    joint_mask = joint_density == 0

    m_joint[~joint_mask] = (
        joint_map[~joint_mask] / joint_density[~joint_mask, np.newaxis]
    )[:, len(f_ghz["ll_ss"]) :] - monopole
    m_joint[joint_mask] = np.nan  # hp.UNSEEN

    curr_path = f"{fits_path}maps/joint/{folder}/"
    if not os.path.exists(curr_path):
        os.makedirs(curr_path)
    for freq in range(len(f_ghz["ll_ss"]), len(f_ghz["lh_ss"])):
        if g.PNG:
            hp.mollview(
                m_joint[:, (freq - len(f_ghz["ll_ss"]))],
                title=f"{int(f_ghz['lh_ss'][freq]):04d} GHz",
                unit="MJy/sr",
                min=0,
                max=200,
            )
            plt.savefig(
                f"{g.SAVE_PATH}maps/joint/{folder}/{int(f_ghz['lh_ss'][freq]):04d}_nside{g.NSIDE}.png"
            )
            plt.close()
        if g.FITS:
            fits.writeto(f"{curr_path}{int(f_ghz['lh_ss'][freq]):04d}_nside{g.NSIDE}.fits", m_joint[:, (freq - len(f_ghz["ll_ss"]))], overwrite=True)