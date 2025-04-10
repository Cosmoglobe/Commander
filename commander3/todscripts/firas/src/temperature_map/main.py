"""
This script takes the interferograms into spectra.
"""

import os
import sys
from datetime import datetime

import h5py
import healpy as hp
import numpy as np
from astropy.io import fits
from globals import OFFSET, PROCESSED_DATA_PATH

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from my_utils import clean_ifg, filter_junk, ifg_to_spec, planck, tune_pointing
from utils.config import gen_nyquistl

T_CMB = 2.72548  # Fixsen 2009
NSIDE = 32
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}
# channels = ["ll"]
modes = {"ss": 0, "lf": 3}  # can change when i have the new cal models

sky_data = h5py.File(
    "../../data/sky_v4.3_debug.h5",
    "r",
)

# print(f"sky_data keys: {sky_data["df_data"].keys()}")
# cal_data = h5py.File(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/cal_v1.h5",
#     "r",
# )

print("getting variables")
variable_names = [
    "gmt",
    "a_ical",
    "b_ical",
    "a_dihedral",
    "b_dihedral",
    "a_refhorn",
    "b_refhorn",
    "a_skyhorn",
    "b_skyhorn",
    "mtm_length",
    "mtm_speed",
    "pix_gal",
    "pix_terr",
    "pix_ecl",
    "a_bol_assem_rh",
    "b_bol_assem_rh",
    "a_bol_assem_rl",
    "b_bol_assem_rl",
    "a_bol_assem_lh",
    "b_bol_assem_lh",
    "a_bol_assem_ll",
    "b_bol_assem_ll",
    "stat_word_1",
    "stat_word_5",
    "stat_word_9",
    "stat_word_12",
    "stat_word_13",
    "stat_word_16",
    "lvdt_stat_a",
    "lvdt_stat_b",
    "adds_per_group",
    "sweeps",
    "gal_lat",
    "gal_lon",
]
channel_dependent = [
    "ifg",
    "bol_cmd_bias",
    "bol_volt",
    "gain",
]

for channel in channels:
    for variable_name in channel_dependent:
        variable_names = variable_names + [f"{variable_name}_{channel}"]

variables = {}
for variable_name in variable_names:
    variables[variable_name] = np.array(sky_data["df_data/" + variable_name][()])

# get the galactic longitude and latitude into a vector in order to interpolate

print(
    f"min gal_lon: {np.min(variables['gal_lon'])}, max gal_lon: {np.max(variables['gal_lon'])}"
)
print(
    f"min gal_lat: {np.min(variables['gal_lat'])}, max gal_lat: {np.max(variables['gal_lat'])}"
)

gal_vec = tune_pointing(
    gal_lon=variables["gal_lon"],
    gal_lat=variables["gal_lat"],
    gmt=variables["gmt"],
    mtm_length=variables["mtm_length"],
    mtm_speed=variables["mtm_speed"],
    offset=OFFSET,
)

# filter out nans that come from gal_vec into all variables
valid_indices = ~np.isnan(gal_vec).any(axis=1)
for variable in variables.keys():
    variables[variable] = variables[variable][valid_indices]
gal_vec = gal_vec[valid_indices]

variables["gal_lon"], variables["gal_lat"] = hp.vec2ang(gal_vec, lonlat=True)
variables["pix_gal"] = hp.ang2pix(
    NSIDE, variables["gal_lon"], variables["gal_lat"], lonlat=True
).astype(int)

# print sizes of the ifgs
print("sizes of the ifgs")
for channel in channels.keys():
    print(f"{channel}: {len(variables[f'ifg_{channel}'])}")

# filter out bad data (selected "by eye")
filter_bad = filter_junk(
    variables["stat_word_1"],
    variables["stat_word_5"],
    variables["stat_word_9"],
    variables["stat_word_12"],
    variables["stat_word_13"],
    variables["stat_word_16"],
    variables["lvdt_stat_a"],
    variables["lvdt_stat_b"],
)

for variable in variables.keys():
    variables[variable] = variables[variable][filter_bad]

variables["pix_gal"] = variables["pix_gal"].astype(int)
variables["pix_terr"] = variables["pix_terr"].astype(int)
variables["pix_ecl"] = variables["pix_ecl"].astype(int)
variables["gmt"] = variables["gmt"].astype(str)
variables["gmt"] = np.array(
    [datetime.strptime(gmt, "%Y-%m-%d %H:%M:%S") for gmt in variables["gmt"]]
)

# sky_data.close()


# constraining to the recs with ical temp at certain bins
variables["ical"] = (
    # variables["a_ical"] * 0.1 + variables["b_ical"] * 0.9
    # )  # using their weights
    (variables["a_ical"] + variables["b_ical"])
    / 2
)
variables["dihedral"] = (variables["a_dihedral"] + variables["b_dihedral"]) / 2
variables["refhorn"] = (variables["a_refhorn"] + variables["b_refhorn"]) / 2
variables["skyhorn"] = (variables["a_skyhorn"] + variables["b_skyhorn"]) / 2
variables["bolometer_rh"] = (
    variables["a_bol_assem_rh"] + variables["b_bol_assem_rh"]
) / 2
variables["bolometer_rl"] = (
    variables["a_bol_assem_rl"] + variables["b_bol_assem_rl"]
) / 2
variables["bolometer_lh"] = (
    variables["a_bol_assem_lh"] + variables["b_bol_assem_lh"]
) / 2
variables["bolometer_ll"] = (
    variables["a_bol_assem_ll"] + variables["b_bol_assem_ll"]
) / 2

# start = "89-326-1130"
# start_dt = datetime.strptime(start, "%y-%j-%H%M")

# ends of each mission period
# mission_periods = np.array(
#     [
#         "89-327-2359",
#         "89-343-0151",
#         "90-019-0204",
#         "90-080-0114",
#         "90-128-2359",
#         "90-139-1534",
#         "90-193-1849",
#         "90-207-1103",
#         "90-208-1119",
#         "90-220-0459",
#         "90-264-0936",
#     ]
# )

# ical_temps = [
#     [2.789],
#     [2.758, 2.763, 2.789],
#     [2.759, 2.771],
#     [2.758, 2.771],
#     [2.758, 2.771],
#     [2.758, 2.770],
#     [2.7455, 2.755, 2.768],
#     [2.746, 2.757, 2.769],
#     [2.757, 2.769],
#     [2.758, 2.770],
#     [2.758, 2.771],
# ]

# missions_periods_dt = [
#     datetime.strptime(period, "%y-%j-%H%M") for period in mission_periods
# ]

# period_filter = np.empty((len(mission_periods), len(variables["gmt"])), dtype=bool)
# period_filter[0] = (variables["gmt"] < missions_periods_dt[0]) & (
#     variables["gmt"] > start_dt
# )

# for i in range(len(missions_periods_dt) - 1):
#     period_filter[i + 1] = (variables["gmt"] < missions_periods_dt[i + 1]) & (
#         variables["gmt"] > missions_periods_dt[i]
#     )

# get filter for ical temps to be +- 2 mK for each mission period depending on the ical temps
# ical_periods = np.empty((len(mission_periods), len(variables["ical"])), dtype=bool)

# for i in range(len(ical_temps)):
#     for j in range(len(ical_temps[i])):
#         if j == 0:
#             ical_periods[i] = (ical_temps[i][j] - 0.002 < variables["ical"]) & (
#                 variables["ical"] < ical_temps[i][j] + 0.002
#             )
#         else:
#             ical_periods[i] = ical_periods[i] | (
#                 ical_temps[i][j] - 0.002 < variables["ical"]
#             ) & (variables["ical"] < ical_temps[i][j] + 0.002)

# full ical temp and mission filter has to be with and between the two
# ical_filter = (
#     (ical_periods[0] & period_filter[0])
#     | (ical_periods[1] & period_filter[1])
#     | (ical_periods[2] & period_filter[2])
#     | (ical_periods[3] & period_filter[3])
#     | (ical_periods[4] & period_filter[4])
#     | (ical_periods[5] & period_filter[5])
#     | (ical_periods[6] & period_filter[6])
#     | (ical_periods[7] & period_filter[7])
#     | (ical_periods[8] & period_filter[8])
#     | (ical_periods[9] & period_filter[9])
#     | (ical_periods[10] & period_filter[10])
# )

# for variable in variables.keys():
#     variables[variable] = variables[variable][ical_filter]

# filtering the data based on the mode used
short_filter = variables["mtm_length"] == 0
long_filter = variables["mtm_length"] == 1
slow_filter = variables["mtm_speed"] == 0
fast_filter = variables["mtm_speed"] == 1

filters = {}
filters["ss"] = short_filter & slow_filter
filters["lf"] = long_filter & fast_filter

variablesm = {}
for variable in variables.keys():
    for mode in modes.keys():
        variablesm[f"{variable}_{mode}"] = variables[variable][filters[mode]]

for channel in channels.keys():
    for mode in modes.keys():
        print(f"{channel}{mode}: {len(variablesm[f"ifg_{channel}_{mode}"])}")

fits_data = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            fits_data[f"{channel}_{mode}"] = fits.open(
                f"/mn/stornext/d16/cmbco/ola/firas/pub_calibration_model/FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
            )

fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

frec = {}
for channel, channel_value in channels.items():
    for mode, mode_value in modes.items():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            frec[f"{channel}_{mode}"] = 4 * (channel_value % 2) + mode_value

# optical transfer function
otf = {}
for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            otf[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RTRANSFE"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["ITRANSFE"][0]
            )
            # plt.show()
            otf[f"{channel}_{mode}"] = otf[f"{channel}_{mode}"][
                np.abs(otf[f"{channel}_{mode}"]) > 0
            ]

apod = {}
etf = {}
S0 = {}
tau = {}
Jo = {}
Jg = {}
Tbol = {}
T0 = {}
R0 = {}
rho = {}
G1 = {}
beta = {}
C3 = {}
C1 = {}

# bolometer function
for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            apod[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "APODIZAT"
            ][0]
            etf[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RELEX_GA"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["IELEX_GA"][0]
            )
            etf[f"{channel}_{mode}"] = etf[f"{channel}_{mode}"][
                np.abs(etf[f"{channel}_{mode}"]) > 0
            ]
            S0[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "DC_RESPO"
            ][0]
            tau[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "TIME_CON"
            ][0]
            Jo[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM8"
            ][0]
            Jg[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM9"
            ][0]
            Tbol[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLOM_B2"
            ][0]
            T0[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM2"
            ][0]
            R0[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM_"
            ][0]
            rho[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM5"
            ][0]
            G1[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM3"
            ][0]
            beta[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM4"
            ][0]
            C3[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM7"
            ][0]
            C1[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM6"
            ][0]

# print("cleaning interferograms")

# for channel, channel_value in channels.items():
#     for mode in modes.keys():
#         if mode == "lf" and (channel == "lh" or channel == "rh"):
#             continue
#         else:
#             variablesm[f"ifg_{channel}_{mode}"] = clean_ifg(
#                 ifg=variablesm[f"ifg_{channel}_{mode}"],
#                 mtm_length=0 if mode[0] == "s" else 1,
#                 mtm_speed=0 if mode[1] == "s" else 1,
#                 # channel=channel_value,
#                 # adds_per_group=variablesm[f"adds_per_group_{channel}_{mode}"],
#                 gain=variablesm[f"gain_{mode}"],
#                 sweeps=variablesm[f"sweeps_{mode}"],
#                 apod=apod[f"{channel}_{mode}"],
#             )

print("converting interferograms to spectra")

# instrumental gain function normalization (mjy v / w sr)
# norm = {"ss": 7.479336e00, "lf": 2.991734e00}

spec = {}
for channel, channel_value in channels.items():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            print(f"ifg to spec of {channel}_{mode}")
            afreq, spec[f"{channel}_{mode}"] = ifg_to_spec(
                ifg=variablesm[f"ifg_{channel}_{mode}"],
                mtm_speed=0 if mode[1] == "s" else 1,
                channel=channel_value,
                adds_per_group=variablesm[f"adds_per_group_{mode}"],
                # adds_per_group=variablesm[f"adds_per_group_{channel}_{mode}"],
                bol_cmd_bias=variablesm[f"bol_cmd_bias_{channel}_{mode}"] / 25.5,  # needs this factor to put it into volts (from pipeline)
                bol_volt=variablesm[f"bol_volt_{channel}_{mode}"],
                fnyq_icm=fnyq["icm"][frec[f"{channel}_{mode}"]],
                fnyq_hz=fnyq["hz"][frec[f"{channel}_{mode}"]],
                otf=otf[f"{channel}_{mode}"],
                Jo=Jo[f"{channel}_{mode}"],
                Jg=Jg[f"{channel}_{mode}"],
                Tbol=Tbol[f"{channel}_{mode}"],
                rho=rho[f"{channel}_{mode}"],
                R0=R0[f"{channel}_{mode}"],
                T0=T0[f"{channel}_{mode}"],
                beta=beta[f"{channel}_{mode}"],
                G1=G1[f"{channel}_{mode}"],
                # tau=tau[f"{channel}_{mode}"],
                C3=C3[f"{channel}_{mode}"],
                C1=C1[f"{channel}_{mode}"],
                # etf=etf[f"{channel}_{mode}"],
                # S0=S0[f"{channel}_{mode}"],
                # norm=norm[mode],
                mtm_length=0 if mode[0] == "s" else 1,
                # gain=variablesm[f"gain_{mode}"],
                gain=variablesm[f"gain_{channel}_{mode}"],
                sweeps=variablesm[f"sweeps_{mode}"],
                # sweeps=variablesm[f"sweeps_{channel}_{mode}"],
                apod=apod[f"{channel}_{mode}"],
            )
            # print(f"shape of spec: {spec[f"{channel}_{mode}"].shape}")

print("making the diff")

# frequency mapping
nu0 = {"ss": 68.020812, "lf": 23.807283}
dnu = {"ss": 13.604162, "lf": 3.4010405}
nf = {"lh_ss": 210, "ll_lf": 182, "ll_ss": 43, "rh_ss": 210, "rl_lf": 182, "rl_ss": 43}

f_ghz = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            f_ghz[f"{channel}_{mode}"] = np.linspace(
                nu0[mode],
                nu0[mode] + dnu[mode] * (nf[f"{channel}_{mode}"] - 1),
                nf[f"{channel}_{mode}"],
            )

# ical spectrum
bb_ical = {}
ical_emiss = {}
bb_dihedral = {}
dihedral_emiss = {}
bb_refhorn = {}
refhorn_emiss = {}
bb_skyhorn = {}
skyhorn_emiss = {}
bb_bolometer_rh = {}
bb_bolometer_rl = {}
bb_bolometer_lh = {}
bb_bolometer_ll = {}
bolometer_emiss = {}
for channel in channels.keys():
    for mode in modes.keys():
        if mode == "lf" and (channel == "lh" or channel == "rh"):
            continue
        else:
            bb_ical[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"],
                variablesm[f"ical_{mode}"],
            )
            ical_emiss[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RICAL"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["IICAL"][0]
            )
            ical_emiss[f"{channel}_{mode}"] = ical_emiss[f"{channel}_{mode}"][
                np.abs(ical_emiss[f"{channel}_{mode}"]) > 0
            ]

            # dihedral spectrum
            bb_dihedral[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"],
                variablesm[f"dihedral_{mode}"],
            )
            dihedral_emiss[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RDIHEDRA"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["IDIHEDRA"][0]
            )
            dihedral_emiss[f"{channel}_{mode}"] = dihedral_emiss[f"{channel}_{mode}"][
                np.abs(dihedral_emiss[f"{channel}_{mode}"]) > 0
            ]

            bb_refhorn[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"],
                variablesm[f"refhorn_{mode}"],
            )
            refhorn_emiss[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RREFHORN"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["IREFHORN"][0]
            )
            refhorn_emiss[f"{channel}_{mode}"] = refhorn_emiss[f"{channel}_{mode}"][
                np.abs(refhorn_emiss[f"{channel}_{mode}"]) > 0
            ]

            # skyhorn spectrum
            bb_skyhorn[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"],
                variablesm[f"skyhorn_{mode}"],
            )
            skyhorn_emiss[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RSKYHORN"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["ISKYHORN"][0]
            )
            skyhorn_emiss[f"{channel}_{mode}"] = skyhorn_emiss[f"{channel}_{mode}"][
                np.abs(skyhorn_emiss[f"{channel}_{mode}"]) > 0
            ]

            # bolometer spectrum
            bb_bolometer_rh[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"],
                variablesm[f"bolometer_rh_{mode}"],
            )
            bb_bolometer_rl[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"],
                variablesm[f"bolometer_rl_{mode}"],
            )
            bb_bolometer_lh[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"],
                variablesm[f"bolometer_lh_{mode}"],
            )
            bb_bolometer_ll[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"],
                variablesm[f"bolometer_ll_{mode}"],
            )
            bolometer_emiss[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RBOLOMET"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["IBOLOMET"][0]
            )
            bolometer_emiss[f"{channel}_{mode}"] = bolometer_emiss[f"{channel}_{mode}"][
                np.abs(bolometer_emiss[f"{channel}_{mode}"]) > 0
            ]


sky = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            if mode[1] == "s":
                cutoff = 5
            else:
                cutoff = 7
            print(f"channel: {channel}, mode: {mode}")
            print(
                f"shape of bb_ical: {bb_ical[f'{channel}_{mode}'].shape}, ical_emiss: {ical_emiss[f'{channel}_{mode}'].shape}, bb_bolometer_rh: {bb_bolometer_rh[f'{channel}_{mode}'].shape}, bolometer_emiss: {bolometer_emiss[f'{channel}_{mode}'].shape}"
            )
            sky[f"{channel}_{mode}"] = (
                spec[f"{channel}_{mode}"][
                    :, cutoff : (len(otf[f"{channel}_{mode}"]) + cutoff)
                ]
                - (
                    (bb_ical[f"{channel}_{mode}"] * ical_emiss[f"{channel}_{mode}"])
                    + (
                        bb_dihedral[f"{channel}_{mode}"]
                        * dihedral_emiss[f"{channel}_{mode}"]
                    )
                    + (
                        bb_refhorn[f"{channel}_{mode}"]
                        * refhorn_emiss[f"{channel}_{mode}"]
                    )
                    + (
                        bb_skyhorn[f"{channel}_{mode}"]
                        * skyhorn_emiss[f"{channel}_{mode}"]
                    )
                    + (
                        bb_bolometer_rh[f"{channel}_{mode}"]
                        * bolometer_emiss[f"{channel}_{mode}"]
                    )
                    + (
                        bb_bolometer_rl[f"{channel}_{mode}"]
                        * bolometer_emiss[f"{channel}_{mode}"]
                    )
                    + (
                        bb_bolometer_lh[f"{channel}_{mode}"]
                        * bolometer_emiss[f"{channel}_{mode}"]
                    )
                    + (
                        bb_bolometer_ll[f"{channel}_{mode}"]
                        * bolometer_emiss[f"{channel}_{mode}"]
                    )
                )
                / otf[f"{channel}_{mode}"][np.newaxis, :]
            )

print("saving sky")

# other varibles that i just want to save in the same place
extra_variables = [
    "stat_word_1",
    "stat_word_12",
    "stat_word_4",
    "stat_word_8",
    "power_a_status_a",
    "power_a_status_b",
    "power_b_status_a",
    "power_b_status_b",
    "ref_hrn_temp_a",
    "ref_hrn_temp_b",
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
    "micro_stat_bus_1",
    "micro_stat_bus_2",
    "micro_stat_bus_3",
    "micro_stat_bus_4",
    "ext_cal_temp_a",
    "ext_cal_temp_b",
    "grt_addr_a",
    "grt_addr_b",
    "hot_spot_cmd_a",
    "hot_spot_cmd_b",
    "int_ref_temp_a",
    "int_ref_temp_b",
    # "sky_hrn_temp_a",
    # "sky_hrn_temp_b",
    "scan",
]
for variable in extra_variables:
    for mode in modes.keys():
        variablesm[f"{variable}_{mode}"] = np.array(sky_data["df_data/" + variable])[
            valid_indices
        ]

        variablesm[f"{variable}_{mode}"] = variablesm[f"{variable}_{mode}"][filter_bad]
        # variablesm[f"{variable}_{mode}"] = np.array(sky_data["df_data/" + variable])[
        #     ical_filter
        # ]
        # variablesm[f"{variable}_{mode}"] = variablesm[f"{variable}_{mode}"][ical_filter]
        variablesm[f"{variable}_{mode}"] = variablesm[f"{variable}_{mode}"][
            filters[mode]
        ]

# save the sky
np.savez(
    PROCESSED_DATA_PATH,
    **variablesm,
    **sky,
)
