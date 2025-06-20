import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import globals as g
from utils.config import gen_nyquistl

# m = hp.fitsfunc.read_map(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/maps/fits/0179.fits",
# )

# hp.mollview(m, coord="G", title="179", unit="MJy/sr", min=0, max=500)
# plt.show()

modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}

reference_maps = "/mn/stornext/d16/cmbco/ola/firas/healpix_maps/"

# comparing otf with ical emiss ---------------
# fits_data_ss = fits.open(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/reference/FIRAS_CALIBRATION_MODEL_LLSS.FITS"
# )

# otf_ss = fits_data_ss[1].data["RTRANSFE"][0] + 1j * fits_data_ss[1].data["ITRANSFE"][0]
# otf_ss = otf_ss[np.abs(otf_ss) > 0]

# ical_emiss_ss = fits_data_ss[1].data["RICAL"][0] + 1j * fits_data_ss[1].data["IICAL"][0]
# ical_emiss_ss = ical_emiss_ss[np.abs(ical_emiss_ss) > 0]

# plt.plot(ical_emiss_ss / otf_ss)
# plt.show()

data = np.load(
    g.PROCESSED_DATA_PATH,
    allow_pickle=True,
)

sky = {}
spec = {}
pix_gal = {}
variables = {}
variable_names = [
    # "gmt",
    "a_ical",
    "b_ical",
    "a_dihedral",
    "b_dihedral",
    "a_refhorn",
    "b_refhorn",
    "a_skyhorn",
    "b_skyhorn",
    "a_bol_assem_rh",
    "b_bol_assem_rh",
    "a_bol_assem_rl",
    "b_bol_assem_rl",
    "a_bol_assem_lh",
    "b_bol_assem_lh",
    "a_bol_assem_ll",
    "b_bol_assem_ll",
    # "mtm_length",
    # "mtm_speed",
    # "pix_gal",
    # # "adds_per_group_ll",
    # "bol_cmd_bias_ll",
    # "bol_cmd_bias_lh",
    # "bol_cmd_bias_rl",
    # "bol_cmd_bias_rh",
    # "dwell_stat_a",
    # "dwell_stat_b",
    # "engstat_spares_1",
    # "engstat_spares_2",
    # "engstat_spares_3",
    # "engstat_spares_4",
    # "engstat_spares_5",
    # "engstat_spares_6",
    # "engstat_spares_7",
    # "engstat_spares_8",
    # "engstat_spares_9",
    # "engstat_spares_10",
    # "engstat_spares2_1",
    # "engstat_spares2_2",
    # "engstat_spares2_3",
    # "engstat_spares2_4",
    # "engstat_spares2_5",
    # "ext_cal_temp_a",
    # "ext_cal_temp_b",
    # "stat_word_1",
    # "stat_word_12",
    # "stat_word_13",
    # "stat_word_16",
    # "stat_word_4",
    # "stat_word_5",
    # "stat_word_8",
    # "stat_word_9",
    # "lvdt_stat_a",
    # "lvdt_stat_b",
    # "power_a_status_a",
    # "power_a_status_b",
    # "power_b_status_a",
    # "power_b_status_b",
    # "ref_hrn_temp_a",
    # "ref_hrn_temp_b",
    # "micro_stat_bus_1",
    # "micro_stat_bus_2",
    # "micro_stat_bus_3",
    # "micro_stat_bus_4",
    # "grt_addr_a",
    # "grt_addr_b",
    # "hot_spot_cmd_a",
    # "hot_spot_cmd_b",
    # "int_ref_temp_a",
    # "int_ref_temp_b",
    # "sky_hrn_temp_a",
    # "sky_hrn_temp_b",
]

for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            sky[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]
            print(f"sky[{channel}_{mode}]: {sky[f"{channel}_{mode}"].shape}")
            spec[f"{channel}_{mode}"] = np.abs(sky[f"{channel}_{mode}"])
            print(f"spec[{channel}_{mode}]: {spec[f"{channel}_{mode}"].shape}")
            pix_gal[mode] = data[f"pix_gal_{mode}"]
            # variables = data["variables"]
            for variable_name in variable_names:
                variables[f"{variable_name}_{mode}"] = data[f"{variable_name}_{mode}"]

# data = h5py.File(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v4.2.h5",
#     "r",
# )
fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

mask = hp.read_map("/mn/stornext/d16/cmbco/ola/masks/HI_mask_4e20_n1024.fits")
mask_alm = hp.sphtfunc.map2alm(mask, pol=False)
mask = hp.alm2map(mask_alm, g.NSIDE, pol=False)
mask = np.where(mask < 0.5, 0, 1)

# remove data inside the mask
for mode in modes.keys():
    for name in variable_names:
        variables[f"{name}_{mode}"] = variables[f"{name}_{mode}"][
            mask[pix_gal[mode]] == 1
        ]
        print(f"size of {name}_{mode}: {variables[f'{name}_{mode}'].shape}")
    for channel in channels.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            spec[f"{channel}_{mode}"] = spec[f"{channel}_{mode}"][
                mask[pix_gal[mode]] == 1
            ]
            print(f"size of spec: {spec[f"{channel}_{mode}"].shape}")

# for freq in range(1, len(spec[0])):
#     plt.plot(ical, spec[:, (freq - 1)], ".")
#     plt.xlabel("ical")
#     plt.ylabel("spec")
#     plt.title(f"Frequency {int(f_ghz[(freq)]):04d} GHz")
#     plt.savefig(
#         f"/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/plots/ical_vs_spec/{int(f_ghz[(freq)]):04d}.png"
#     )

#     if int(f_ghz[(freq)]) == 204:
#         plt.show()
#     plt.clf()

print("plotting")

n_cols = 6
n_rows = len(variable_names) // n_cols + 1
for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            print(f"Mode: {f"{channel}_{mode}"}")
            print(f"Shape of spec: {spec[f"{channel}_{mode}"].shape}")

            fig, ax = plt.subplots(n_rows, n_cols, figsize=(20, 20), sharex=True)

            for i, name in enumerate(variable_names):
                ax[i // n_cols, i % n_cols].plot(variables[f"{name}_{mode}"])
                ax[i // n_cols, i % n_cols].set_title(name, fontsize=10)

            ax[-1, -1].imshow(
                np.abs(spec[f"{channel}_{mode}"]).T,
                aspect="auto",
                extent=[
                    0,
                    len(spec[f"{channel}_{mode}"]),
                    0,
                    len(spec[f"{channel}_{mode}"][0]),
                ],
                vmax=500,
                vmin=0,
            )
            ax[-1, -1].set_title("spec", fontsize=10)

            plt.tight_layout()
            plt.show()
