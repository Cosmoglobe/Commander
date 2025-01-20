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
modes = ["ss", "sf", "ls", "lf", "fs", "fl"]
# modes = ["ss"]

sky_data = h5py.File(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v4.2.h5",
    "r",
)
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
    "adds_per_group",
    "bol_cmd_bias",
    "bol_volt",
    "gain",
    "sweeps",
    "ifg",
]

for channel in channels:
    for variable_name in channel_dependent:
        variable_names = variable_names + [f"{variable_name}_{channel}"]

variables = {}
for variable_name in variable_names:
    variables[variable_name] = np.array(sky_data["df_data/" + variable_name][()])

variables["pix_gal"] = variables["pix_gal"].astype(int)
variables["gmt"] = variables["gmt"].astype(str)
variables["gmt"] = np.array(
    [datetime.strptime(gmt, "%Y-%m-%d %H:%M:%S") for gmt in variables["gmt"]]
)

sky_data.close()

# using only ss data for now to test
# make filter with true/false where mtm_length and mtm_speed are 0
short_filter = variables["mtm_length"] == 0
slow_filter = variables["mtm_speed"] == 0

for variable in variables.keys():
    variables[variable] = variables[variable][short_filter & slow_filter]

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

print(f"data loaded: {len(variables['ifg_ll'])}")

fits_data_ss = fits.open(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/reference/FIRAS_CALIBRATION_MODEL_LLSS.FITS"
)
fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

scan_mode = 0  # SS
channel = 3  # LL
frec = 4 * (channel % 2) + scan_mode

# optical transfer function
otf_ss = fits_data_ss[1].data["RTRANSFE"][0] + 1j * fits_data_ss[1].data["ITRANSFE"][0]
otf_ss = otf_ss[np.abs(otf_ss) > 0]

# bolometer function
tau_ss = fits_data_ss[1].data["TIME_CON"][0]
Jo_ss = fits_data_ss[1].data["BOLPARM8"][0]
Jg_ss = fits_data_ss[1].data["BOLPARM9"][0]
Tbol_ss = fits_data_ss[1].data["BOLOM_B2"][0]
T0_ss = fits_data_ss[1].data["BOLPARM2"][0]
R0_ss = fits_data_ss[1].data["BOLPARM_"][0]
rho_ss = fits_data_ss[1].data["BOLPARM5"][0]
G1_ss = fits_data_ss[1].data["BOLPARM3"][0]
beta_ss = fits_data_ss[1].data["BOLPARM4"][0]

# plot ifgs over time
plt.imshow(
    variables["ifg_ll"].T,
    aspect="auto",
    extent=[0, len(variables["ifg_ll"]), 0, len(variables["ifg_ll"][0])],
)
plt.savefig("../../output/plots/ifg_over_time.png")
plt.clf()

# plot random ifg
plt.plot(variables["ifg_ll"][np.random.randint(0, len(variables["ifg_ll"]))])
plt.savefig("../../output/plots/random_ifg.png")
plt.clf()

print("cleaning interferograms")

for variable in variables.keys():
    variables[variable] = variables[variable][ical_filter]

variables["ifg_ll"] = clean_ifg(
    ifg=variables["ifg_ll"],
    mtm_length=variables["mtm_length"],
    mtm_speed=variables["mtm_speed"],
    channel=3,
    adds_per_group=variables["adds_per_group_ll"],
    gain=variables["gain_ll"],
    sweeps=variables["sweeps_ll"],
)

# plot cleaned ifgs
plt.imshow(
    variables["ifg_ll"].T,
    aspect="auto",
    extent=[0, len(variables["ifg_ll"]), 0, len(variables["ifg_ll"][0])],
)
plt.savefig("../../output/plots/cleaned_ifg_over_time.png")
plt.clf()

# check if nans in ifgs
if np.isnan(variables["ifg_ll"]).any():
    print("Nans in interferograms")
    print(variables["ifg_ll"][np.isnan(variables["ifg_ll"])])

print("converting interferograms to spectra")

spec_ss = ifg_to_spec(
    variables["ifg_ll"],
    variables["mtm_speed"],
    channel,
    variables["adds_per_group_ll"],
    variables["bol_cmd_bias_ll"],
    variables["bol_volt_ll"],
    fnyq["icm"][frec],
    fnyq["hz"][frec],
    otf_ss,
    Jo_ss,
    Jg_ss,
    Tbol_ss,
    rho_ss,
    R0_ss,
    T0_ss,
    beta_ss,
    G1_ss,
    tau_ss,
)

print("plotting")

plt.imshow(
    np.abs(spec_ss).T, aspect="auto", extent=[0, len(spec_ss), 0, len(spec_ss[0])]
)
plt.savefig("../../output/tests/spectra_diff_over_time.png")
plt.clf()

print("making the diff")

# frequency mapping
f_icm = np.arange(210) * (fnyq["icm"][frec] / 320) + 2
c = 3e8 * 1e2  # cm/s
f_ghz = f_icm * c * 1e-9

# ical spectrum
bb_ical = planck(f_ghz[: len(spec_ss)], variables["ical"])
ical_emiss_ss = fits_data_ss[1].data["RICAL"][0] + 1j * fits_data_ss[1].data["IICAL"][0]
ical_emiss_ss = ical_emiss_ss[: len(spec_ss)]

# dihedral spectrum
bb_dihedral = planck(f_ghz[: len(spec_ss)], variables["dihedral"])
dihedral_emiss_ss = (
    fits_data_ss[1].data["RDIHEDRA"][0] + 1j * fits_data_ss[1].data["IDIHEDRA"][0]
)
dihedral_emiss_ss = dihedral_emiss_ss[: len(spec_ss)]

# refhorn spectrum
bb_refhorn = planck(f_ghz[: len(spec_ss)], variables["refhorn"])
refhorn_emiss_ss = (
    fits_data_ss[1].data["RREFHORN"][0] + 1j * fits_data_ss[1].data["IREFHORN"][0]
)
refhorn_emiss_ss = refhorn_emiss_ss[: len(spec_ss)]

# skyhorn spectrum
bb_skyhorn = planck(f_ghz[: len(spec_ss)], variables["skyhorn"])
skyhorn_emiss_ss = (
    fits_data_ss[1].data["RSKYHORN"][0] + 1j * fits_data_ss[1].data["ISKYHORN"][0]
)
skyhorn_emiss_ss = skyhorn_emiss_ss[: len(spec_ss)]

# bolometer spectrum
bb_bolometer_rh = planck(f_ghz[: len(spec_ss)], variables["bolometer_rh"])
bb_bolometer_rl = planck(f_ghz[: len(spec_ss)], variables["bolometer_rl"])
bb_bolometer_lh = planck(f_ghz[: len(spec_ss)], variables["bolometer_lh"])
bb_bolometer_ll = planck(f_ghz[: len(spec_ss)], variables["bolometer_ll"])

bolometer_emiss_ss = (
    fits_data_ss[1].data["RBOLOMET"][0] + 1j * fits_data_ss[1].data["IBOLOMET"][0]
)

sky_ss = (
    spec_ss
    - (
        (bb_ical * ical_emiss_ss)[:, : len(spec_ss[0])]
        + (bb_dihedral * dihedral_emiss_ss)[:, : len(spec_ss[0])]
        + (bb_refhorn * refhorn_emiss_ss)[:, : len(spec_ss[0])]
        + (bb_skyhorn * skyhorn_emiss_ss)[:, : len(spec_ss[0])]
        + (bb_bolometer_rh * bolometer_emiss_ss)[:, : len(spec_ss[0])]
        + (bb_bolometer_rl * bolometer_emiss_ss)[:, : len(spec_ss[0])]
        + (bb_bolometer_lh * bolometer_emiss_ss)[:, : len(spec_ss[0])]
        + (bb_bolometer_ll * bolometer_emiss_ss)[:, : len(spec_ss[0])]
    )
    / otf_ss[np.newaxis, :]
)

print("saving sky")

# save the sky
np.savez(
    "../../output/data/sky.npz",
    **variables,
    sky=sky_ss,
    allow_pickle=True,
)
