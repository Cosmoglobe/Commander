import h5py
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from utils.config import gen_nyquistl

# m = hp.fitsfunc.read_map(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/maps/fits/0179.fits",
# )

# hp.mollview(m, coord="G", title="179", unit="MJy/sr", min=0, max=500)
# plt.show()

modes = {"ss": 0, "lf": 3}

reference_maps = "/mn/stornext/d16/cmbco/ola/firas/healpix_maps/"

# comparing otf with ical emiss ---------------
fits_data_ss = fits.open(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/reference/FIRAS_CALIBRATION_MODEL_LLSS.FITS"
)

otf_ss = fits_data_ss[1].data["RTRANSFE"][0] + 1j * fits_data_ss[1].data["ITRANSFE"][0]
otf_ss = otf_ss[np.abs(otf_ss) > 0]

ical_emiss_ss = fits_data_ss[1].data["RICAL"][0] + 1j * fits_data_ss[1].data["IICAL"][0]
ical_emiss_ss = ical_emiss_ss[np.abs(ical_emiss_ss) > 0]

# plt.plot(ical_emiss_ss / otf_ss)
# plt.show()

data = np.load(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/data/sky.npz",
    allow_pickle=True,
)

variable_names = [
    # "gmt",
    # "a_ical",
    # "b_ical",
    # "a_dihedral",
    # "b_dihedral",
    # "a_refhorn",
    # "b_refhorn",
    # "a_skyhorn",
    # "b_skyhorn",
    # "mtm_length",
    # "mtm_speed",
    # "pix_gal",
    # "adds_per_group_ll",
    "bol_cmd_bias_ll",
    "dwell_stat_a",
    "dwell_stat_b",
    "engstat_spares_1",
    "engstat_spares_2",
    "engstat_spares_3",
    "engstat_spares_4",
    "engstat_spares_5",
    "engstat_spares_6",
    "engstat_spares_7",
    "engstat_spares_8",
    "engstat_spares_9",
    "engstat_spares_10",
    "engstat_spares2_1",
    "engstat_spares2_2",
    "engstat_spares2_3",
    "engstat_spares2_4",
    "engstat_spares2_5",
    "ext_cal_temp_a",
    "ext_cal_temp_b",
    "stat_word_1",
    "stat_word_12",
    "stat_word_13",
    "stat_word_16",
    "stat_word_4",
    "stat_word_5",
    "stat_word_8",
    "stat_word_9",
    "lvdt_stat_a",
    "lvdt_stat_b",
    "power_a_status_a",
    "power_a_status_b",
    "power_b_status_a",
    "power_b_status_b",
    "ref_hrn_temp_a",
    "ref_hrn_temp_b",
    "micro_stat_bus_1",
    "micro_stat_bus_2",
    "micro_stat_bus_3",
    "micro_stat_bus_4",
    "grt_addr_a",
    "grt_addr_b",
    "hot_spot_cmd_a",
    "hot_spot_cmd_b",
    "int_ref_temp_a",
    "int_ref_temp_b",
    # "sky_hrn_temp_a",
    # "sky_hrn_temp_b",
]

sky = {}
spec = {}
pix_gal = {}
for mode in modes.keys():
    sky[mode] = data[mode]
    print(f"sky[mode]: {sky[mode].shape}")
    spec[mode] = np.abs(sky[mode][:, 1:])
    print(f"spec[mode]: {spec[mode].shape}")
    pix_gal[mode] = data[f"pix_gal_{mode}"]
# variables = data["variables"]
variables = {}
for mode in modes.keys():
    for variable_name in variable_names:
        variables[f"{variable_name}_{mode}"] = data[f"{variable_name}_{mode}"]
        print(f"{variable_name}_{mode}: {variables[f'{variable_name}_{mode}'].shape}")

# data = h5py.File(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v4.2.h5",
#     "r",
# )
fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

# scan_mode = 0  # SS
channel = 3  # LL
frec = {}
for mode, item in modes.items():
    frec[mode] = 4 * (channel % 2) + item

# f_icm = np.arange(len(sky[0])) * (fnyq["icm"][frec] / 320)
# c = 3e8 * 1e2  # cm/s
# f_ghz = f_icm * c * 1e-9 + 55

mask = fits.open("BP_CMB_I_analysis_mask_n1024_v2.fits")
mask = mask[1].data.astype(int)

# remove data inside the mask
for mode in modes:
    spec[mode] = spec[mode][mask[pix_gal[mode]] == 1]
    print(f"size of spec: {spec[mode].shape}")
    for name in variable_names:
        variables[f"{name}_{mode}"] = variables[f"{name}_{mode}"][
            mask[pix_gal[mode]] == 1
        ]
        print(f"size of {name}_{mode}: {variables[f'{name}_{mode}'].shape}")

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
for mode in modes:
    print(f"Mode: {mode}")
    print(f"Shape of spec: {spec[mode].shape}")

    fig, ax = plt.subplots(n_rows, n_cols, figsize=(20, 20), sharex=True)

    for i, name in enumerate(variable_names):
        ax[i // n_cols, i % n_cols].plot(variables[f"{name}_{mode}"])
        ax[i // n_cols, i % n_cols].set_title(name)

    ax[-1, -1].imshow(
        np.abs(spec[mode]).T,
        aspect="auto",
        extent=[0, len(spec[mode]), 0, len(spec[mode][0])],
        vmax=500,
        vmin=0,
    )
    ax[-1, -1].set_title("spec")

    plt.show()
