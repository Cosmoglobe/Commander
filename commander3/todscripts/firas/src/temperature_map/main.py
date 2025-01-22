"""
This script takes the interferograms into spectra.
"""

from datetime import datetime

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from my_utils import clean_ifg, ifg_to_spec, planck
from utils.config import gen_nyquistl

T_CMB = 2.72548  # Fixsen 2009
# channels = ["rh", "rl", "lh", "ll"]
channels = ["ll"]
modes = {"ss": 0, "lf": 3}  # can change when i have the new cal models

sky_data = h5py.File(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v4.3.h5",
    "r",
)

print(f"sky_data keys: {sky_data["df_data"].keys()}")
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
    "a_bol_assem_rh",
    "b_bol_assem_rh",
    "a_bol_assem_rl",
    "b_bol_assem_rl",
    "a_bol_assem_lh",
    "b_bol_assem_lh",
    "a_bol_assem_ll",
    "b_bol_assem_ll",
]
channel_dependent = [
    "ifg",
    "adds_per_group",
    "bol_cmd_bias",
    "bol_volt",
    "gain",
    "sweeps",
]

for channel in channels:
    for variable_name in channel_dependent:
        variable_names = variable_names + [f"{variable_name}_{channel}"]

variables = {}
for variable_name in variable_names:
    variables[variable_name] = np.array(sky_data["df_data/" + variable_name][()])

print(f"variables keys: {variables.keys()}")

variables["pix_gal"] = variables["pix_gal"].astype(int)
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

start = "89-326-1130"
start_dt = datetime.strptime(start, "%y-%j-%H%M")

# ends of each mission period
mission_periods = np.array(
    [
        "89-327-2359",
        "89-343-0151",
        "90-019-0204",
        "90-080-0114",
        "90-128-2359",
        "90-139-1534",
        "90-193-1849",
        "90-207-1103",
        "90-208-1119",
        "90-220-0459",
        "90-264-0936",
    ]
)

ical_temps = [
    [2.789],
    [2.758, 2.763, 2.789],
    [2.759, 2.771],
    [2.758, 2.771],
    [2.758, 2.771],
    [2.758, 2.770],
    [2.7455, 2.755, 2.768],
    [2.746, 2.757, 2.769],
    [2.757, 2.769],
    [2.758, 2.770],
    [2.758, 2.771],
]

missions_periods_dt = [
    datetime.strptime(period, "%y-%j-%H%M") for period in mission_periods
]

period_filter = np.empty((len(mission_periods), len(variables["gmt"])), dtype=bool)
period_filter[0] = (variables["gmt"] < missions_periods_dt[0]) & (
    variables["gmt"] > start_dt
)

for i in range(len(missions_periods_dt) - 1):
    period_filter[i + 1] = (variables["gmt"] < missions_periods_dt[i + 1]) & (
        variables["gmt"] > missions_periods_dt[i]
    )

# get filter for ical temps to be +- 2 mK for each mission period depending on the ical temps
ical_periods = np.empty((len(mission_periods), len(variables["ical"])), dtype=bool)

for i in range(len(ical_temps)):
    for j in range(len(ical_temps[i])):
        if j == 0:
            ical_periods[i] = (ical_temps[i][j] - 0.002 < variables["ical"]) & (
                variables["ical"] < ical_temps[i][j] + 0.002
            )
        else:
            ical_periods[i] = ical_periods[i] | (
                ical_temps[i][j] - 0.002 < variables["ical"]
            ) & (variables["ical"] < ical_temps[i][j] + 0.002)

# full ical temp and mission filter has to be with and between the two
ical_filter = (
    (ical_periods[0] & period_filter[0])
    | (ical_periods[1] & period_filter[1])
    | (ical_periods[2] & period_filter[2])
    | (ical_periods[3] & period_filter[3])
    | (ical_periods[4] & period_filter[4])
    | (ical_periods[5] & period_filter[5])
    | (ical_periods[6] & period_filter[6])
    | (ical_periods[7] & period_filter[7])
    | (ical_periods[8] & period_filter[8])
    | (ical_periods[9] & period_filter[9])
    | (ical_periods[10] & period_filter[10])
)

for variable in variables.keys():
    variables[variable] = variables[variable][ical_filter]

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


for mode in modes.keys():
    print(f"{mode}: {len(variablesm['ifg_ll_' + mode])}")

fits_data = {}
for mode in modes.keys():
    fits_data[mode] = fits.open(
        f"/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/reference/FIRAS_CALIBRATION_MODEL_LL{mode.upper()}.FITS"
    )

fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

channel = 3  # LL
frec = {}
for key, value in modes.items():
    frec[key] = 4 * (channel % 2) + value

# optical transfer function
otf = {}
for mode in modes.keys():
    otf[mode] = (
        fits_data[mode][1].data["RTRANSFE"][0]
        + 1j * fits_data[mode][1].data["ITRANSFE"][0]
    )
    plt.plot(np.abs(otf[mode]))
    # plt.show()
    otf[mode] = otf[mode][np.abs(otf[mode]) > 0]

tau = {}
Jo = {}
Jg = {}
Tbol = {}
T0 = {}
R0 = {}
rho = {}
G1 = {}
beta = {}

# bolometer function
for mode in modes.keys():
    tau[mode] = fits_data[mode][1].data["TIME_CON"][0]
    Jo[mode] = fits_data[mode][1].data["BOLPARM8"][0]
    Jg[mode] = fits_data[mode][1].data["BOLPARM9"][0]
    Tbol[mode] = fits_data[mode][1].data["BOLOM_B2"][0]
    T0[mode] = fits_data[mode][1].data["BOLPARM2"][0]
    R0[mode] = fits_data[mode][1].data["BOLPARM_"][0]
    rho[mode] = fits_data[mode][1].data["BOLPARM5"][0]
    G1[mode] = fits_data[mode][1].data["BOLPARM3"][0]
    beta[mode] = fits_data[mode][1].data["BOLPARM4"][0]

# plot ifgs over time
plt.imshow(
    variablesm["ifg_ll_ss"].T,
    aspect="auto",
    extent=[0, len(variablesm["ifg_ll_ss"]), 0, len(variablesm["ifg_ll_ss"][0])],
)
plt.savefig("../../output/plots/ifg_over_time.png")
plt.clf()

# plot random ifg
plt.plot(variablesm["ifg_ll_ss"][np.random.randint(0, len(variablesm["ifg_ll_ss"]))])
plt.savefig("../../output/plots/random_ifg.png")
plt.clf()

print("cleaning interferograms")

for mode in modes:
    variablesm[f"ifg_ll_{mode}"] = clean_ifg(
        ifg=variablesm[f"ifg_ll_{mode}"],
        mtm_length=0 if mode[0] == "s" else 1,
        mtm_speed=0 if mode[1] == "s" else 1,
        channel=3,
        adds_per_group=variablesm[f"adds_per_group_ll_{mode}"],
        gain=variablesm[f"gain_ll_{mode}"],
        sweeps=variablesm[f"sweeps_ll_{mode}"],
    )

# plot cleaned ifgs
plt.imshow(
    variablesm["ifg_ll_ss"].T,
    aspect="auto",
    extent=[0, len(variablesm["ifg_ll_ss"]), 0, len(variablesm["ifg_ll_ss"][0])],
)
plt.savefig("../../output/plots/cleaned_ifg_over_time.png")
plt.clf()

print("converting interferograms to spectra")

spec = {}
for mode in modes.keys():
    spec[mode] = ifg_to_spec(
        ifg=variablesm[f"ifg_ll_{mode}"],
        mtm_speed=0 if mode[1] == "s" else 1,
        channel=channel,
        adds_per_group=variablesm[f"adds_per_group_ll_{mode}"],
        bol_cmd_bias=variablesm[f"bol_cmd_bias_ll_{mode}"],
        bol_volt=variablesm[f"bol_volt_ll_{mode}"],
        fnyq_icm=fnyq["icm"][frec[mode]],
        fnyq_hz=fnyq["hz"][frec[mode]],
        otf=otf[mode],
        Jo=Jo[mode],
        Jg=Jg[mode],
        Tbol=Tbol[mode],
        rho=rho[mode],
        R0=R0[mode],
        T0=T0[mode],
        beta=beta[mode],
        G1=G1[mode],
        tau=tau[mode],
    )

print("making the diff")

# frequency mapping
c = 3e8 * 1e2  # cm/s
f_icm = {}
f_ghz = {}
for mode in modes.keys():
    f_icm[mode] = np.arange(210) * (fnyq["icm"][frec[mode]] / 320) + 2
    f_ghz[mode] = f_icm[mode] * c * 1e-9

# ical spectrum
bb_ical = {}
ical_emiss = {}
for mode in modes.keys():
    bb_ical[mode] = planck(f_ghz[mode][: len(spec["ss"])], variablesm[f"ical_{mode}"])
    ical_emiss[mode] = (
        fits_data[mode][1].data["RICAL"][0] + 1j * fits_data[mode][1].data["IICAL"][0]
    )
    ical_emiss[mode] = ical_emiss[mode][: len(spec[mode])]

# dihedral spectrum
bb_dihedral = {}
dihedral_emiss = {}
for mode in modes.keys():
    bb_dihedral[mode] = planck(
        f_ghz[mode][: len(spec["ss"])], variablesm[f"dihedral_{mode}"]
    )
    dihedral_emiss[mode] = (
        fits_data[mode][1].data["RDIHEDRA"][0]
        + 1j * fits_data[mode][1].data["IDIHEDRA"][0]
    )
    dihedral_emiss[mode] = dihedral_emiss[mode][: len(spec[mode])]

# refhorn spectrum
bb_refhorn = {}
refhorn_emiss = {}
for mode in modes.keys():
    bb_refhorn[mode] = planck(
        f_ghz[mode][: len(spec["ss"])], variablesm[f"refhorn_{mode}"]
    )
    refhorn_emiss[mode] = (
        fits_data[mode][1].data["RREFHORN"][0]
        + 1j * fits_data[mode][1].data["IREFHORN"][0]
    )
    refhorn_emiss[mode] = refhorn_emiss[mode][: len(spec[mode])]

# skyhorn spectrum
bb_skyhorn = {}
skyhorn_emiss = {}
for mode in modes.keys():
    bb_skyhorn[mode] = planck(
        f_ghz[mode][: len(spec["ss"])], variablesm[f"skyhorn_{mode}"]
    )
    skyhorn_emiss[mode] = (
        fits_data[mode][1].data["RSKYHORN"][0]
        + 1j * fits_data[mode][1].data["ISKYHORN"][0]
    )
    skyhorn_emiss[mode] = skyhorn_emiss[mode][: len(spec[mode])]

# bolometer spectrum
bb_bolometer_rh = {}
bb_bolometer_rl = {}
bb_bolometer_lh = {}
bb_bolometer_ll = {}
bolometer_emiss = {}
for mode in modes.keys():
    bb_bolometer_rh[mode] = planck(
        f_ghz[mode][: len(spec["ss"])], variablesm[f"bolometer_rh_{mode}"]
    )
    bb_bolometer_rl[mode] = planck(
        f_ghz[mode][: len(spec["ss"])], variablesm[f"bolometer_rl_{mode}"]
    )
    bb_bolometer_lh[mode] = planck(
        f_ghz[mode][: len(spec["ss"])], variablesm[f"bolometer_lh_{mode}"]
    )
    bb_bolometer_ll[mode] = planck(
        f_ghz[mode][: len(spec["ss"])], variablesm[f"bolometer_ll_{mode}"]
    )
    bolometer_emiss[mode] = (
        fits_data[mode][1].data["RBOLOMET"][0]
        + 1j * fits_data[mode][1].data["IBOLOMET"][0]
    )

sky = {}
for mode in modes.keys():
    sky[mode] = (
        spec[mode]
        - (
            (bb_ical[mode] * ical_emiss[mode])[:, : len(spec[mode][0])]
            + (bb_dihedral[mode] * dihedral_emiss[mode])[:, : len(spec[mode][0])]
            + (bb_refhorn[mode] * refhorn_emiss[mode])[:, : len(spec[mode][0])]
            + (bb_skyhorn[mode] * skyhorn_emiss[mode])[:, : len(spec[mode][0])]
            + (bb_bolometer_rh[mode] * bolometer_emiss[mode])[:, : len(spec[mode][0])]
            + (bb_bolometer_rl[mode] * bolometer_emiss[mode])[:, : len(spec[mode][0])]
            + (bb_bolometer_lh[mode] * bolometer_emiss[mode])[:, : len(spec[mode][0])]
            + (bb_bolometer_ll[mode] * bolometer_emiss[mode])[:, : len(spec[mode][0])]
        )
        / otf[mode][np.newaxis, :]
    )

print("saving sky")

# other varibles that i just want to save in the same place
extra_variables = [
    "adds_per_group_ll",
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
]
for variable in extra_variables:
    for mode in modes.keys():
        variablesm[f"{variable}_{mode}"] = np.array(sky_data["df_data/" + variable])[
            ical_filter
        ]
        variablesm[f"{variable}_{mode}"] = variablesm[f"{variable}_{mode}"][
            filters[mode]
        ]

# save the sky
np.savez(
    "../../output/data/sky.npz",
    **variablesm,
    **sky,
)
