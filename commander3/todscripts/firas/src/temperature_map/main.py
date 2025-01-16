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
    # "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v1.h5",
    "r",
)
# cal_data = h5py.File(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/cal_v1.h5",
#     "r",
# )

print("parsing gmts")
# getting the xcal data into their own variables
gmts = np.array(sky_data["df_data/gmt"][()]).astype(str)
gmt = np.array([datetime.strptime(gmt, "%Y-%m-%d %H:%M:%S") for gmt in gmts])

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

# dihedral_temps = [ # not using for now because i can't make sense of how they did it
#     [2.14, 2.5, 3.1, 3.7, 4.3, 4.9, 5.5],
#     [2.02, 2.5, 3.1, 3.7, 4.3, 4.9, 5.5],
#     [2.14, 2.5, 3.1, 3.7, 4.3, 4.9, 5.5],
#     [2.14, 2.5, 3.1, 3.7, 4.3, 4.9, 5.5],
#     [1.98, 2.5, 3.1, 3.7, 4.3, 4.9, 5.5],
#     [2.0, 2.5, 3.1, 3.7, 4.3, 4.9, 5.5],
#     [2.03, 2.5, 3.1, 3.7, 4.3, 4.9, 5.5],
#     [2.01, 2.5, 3.1, 3.7, 4.3, 4.9, 5.5],
#     [2.01, 2.5, 3.1, 3.7, 4.3, 4.9, 5.5],
#     [2.0, 2.5, 3.1, 3.7, 4.3, 4.9, 5.5],
#     [2.0, 2.5, 3.1, 3.7, 4.3, 4.9, 5.5],
# ]

missions_periods_dt = [
    datetime.strptime(period, "%y-%j-%H%M") for period in mission_periods
]

print("other variables")
# ical = np.array(sky_data["df_data/ical"][()])
a_ical = np.array(sky_data["df_data/a_ical"][()])
b_ical = np.array(sky_data["df_data/b_ical"][()])
# dihedral = np.array(sky_data["df_data/dihedral"][()])
a_dihedral = np.array(sky_data["df_data/a_dihedral"][()])
b_dihedral = np.array(sky_data["df_data/b_dihedral"][()])
mtm_length = np.array(sky_data["df_data/mtm_length"][()])
mtm_speed = np.array(sky_data["df_data/mtm_speed"][()])
pix_gal = np.array(sky_data["df_data/pix_gal"]).astype(int)

ifg = {}
bol_cmd_bias = {}
bol_volt = {}
adds_per_group = {}
gain = {}
sweeps = {}
for channel in channels:
    ifg[channel] = np.array(sky_data["df_data/ifg_" + channel])
    bol_cmd_bias[channel] = np.array(sky_data["df_data/bol_cmd_bias_" + channel])
    bol_volt[channel] = np.array(sky_data["df_data/bol_volt_" + channel])
    adds_per_group[channel] = np.array(sky_data["df_data/adds_per_group_" + channel])
    gain[channel] = np.array(sky_data["df_data/gain_" + channel])
    sweeps[channel] = np.array(sky_data["df_data/sweeps_" + channel])

sky_data.close()

# using only ss data for now to test
# make filter with true/false where mtm_length and mtm_speed are 0
short_filter = mtm_length == 0
slow_filter = mtm_speed == 0

a_ical = a_ical[short_filter & slow_filter]
b_ical = b_ical[short_filter & slow_filter]
a_dihedral = a_dihedral[short_filter & slow_filter]
b_dihedral = b_dihedral[short_filter & slow_filter]
mtm_length = mtm_length[short_filter & slow_filter]
mtm_speed = mtm_speed[short_filter & slow_filter]
ifg["ll"] = ifg["ll"][short_filter & slow_filter]
bol_cmd_bias["ll"] = bol_cmd_bias["ll"][short_filter & slow_filter]
bol_volt["ll"] = bol_volt["ll"][short_filter & slow_filter]
adds_per_group["ll"] = adds_per_group["ll"][short_filter & slow_filter]
gain["ll"] = gain["ll"][short_filter & slow_filter]
sweeps["ll"] = sweeps["ll"][short_filter & slow_filter]
gmt = gmt[short_filter & slow_filter]
pix_gal = pix_gal[short_filter & slow_filter]

# constraining to the recs with ical temp at certain bins
ical = (
    a_ical + b_ical
) / 2  # doing the mean between the two sides for now, will be fitted later
dihedral = (a_dihedral + b_dihedral) / 2  # same for dihedral

period_filter = np.empty((len(mission_periods), len(gmt)), dtype=bool)
period_filter[0] = (gmt < missions_periods_dt[0]) & (gmt > start_dt)

for i in range(len(missions_periods_dt) - 1):
    period_filter[i + 1] = (gmt < missions_periods_dt[i + 1]) & (
        gmt > missions_periods_dt[i]
    )

# get filter for ical temps to be +- 2 mK for each mission period depending on the ical temps
ical_periods = np.empty((len(mission_periods), len(ical)), dtype=bool)

for i in range(len(ical_temps)):
    for j in range(len(ical_temps[i])):
        if j == 0:
            ical_periods[i] = (ical_temps[i][j] - 0.002 < ical) & (
                ical < ical_temps[i][j] + 0.002
            )
        else:
            ical_periods[i] = ical_periods[i] | (ical_temps[i][j] - 0.002 < ical) & (
                ical < ical_temps[i][j] + 0.002
            )

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

# dihedral filter
# dihedral_periods = np.empty((len(mission_periods), len(dihedral)), dtype=bool)
# for i in range(len(dihedral_temps)):
#     for j in range(len(dihedral_temps[i])):
#         if j == 0:
#             dihedral_periods[i] = dihedral_temps[i][j] == np.round(
#                 dihedral, len(str(dihedral_temps[i][j])[2:])
#             )
#         else:
#             dihedral_periods[i] = dihedral_periods[i] | (
#                 np.round(dihedral, len(str(dihedral_temps[i][j])[2:]))
#                 == dihedral_temps[i][j]
#             )

# dihedral_filter = (
#     (dihedral_periods[0] & period_filter[0])
#     | (dihedral_periods[1] & period_filter[1])
#     | (dihedral_periods[2] & period_filter[2])
#     | (dihedral_periods[3] & period_filter[3])
#     | (dihedral_periods[4] & period_filter[4])
#     | (dihedral_periods[5] & period_filter[5])
#     | (dihedral_periods[6] & period_filter[6])
#     | (dihedral_periods[7] & period_filter[7])
#     | (dihedral_periods[8] & period_filter[8])
#     | (dihedral_periods[9] & period_filter[9])
#     | (dihedral_periods[10] & period_filter[10])
# )

# print(f"dihedral before filter: {len(dihedral)}")

# dihedral = dihedral[dihedral_filter]

# print(f"dihedral after filter: {len(dihedral)}")

print(f"data loaded: {len(ifg['ll'])}")

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
plt.imshow(ifg["ll"].T, aspect="auto", extent=[0, len(ifg["ll"]), 0, len(ifg["ll"][0])])
plt.savefig("../../output/plots/ifg_over_time.png")
plt.clf()

# plot random ifg
plt.plot(ifg["ll"][np.random.randint(0, len(ifg["ll"]))])
plt.savefig("../../output/plots/random_ifg.png")
plt.clf()

print("cleaning interferograms")

ifg["ll"] = ifg["ll"][ical_filter]
mtm_length = mtm_length[ical_filter]
mtm_speed = mtm_speed[ical_filter]
bol_cmd_bias["ll"] = bol_cmd_bias["ll"][ical_filter]
bol_volt["ll"] = bol_volt["ll"][ical_filter]
adds_per_group["ll"] = adds_per_group["ll"][ical_filter]
gain["ll"] = gain["ll"][ical_filter]
sweeps["ll"] = sweeps["ll"][ical_filter]
dihedral = dihedral[ical_filter]
ical = ical[ical_filter]
pix_gal = pix_gal[ical_filter]

ifg["ll"] = clean_ifg(
    ifg=ifg["ll"],
    mtm_length=mtm_length,
    mtm_speed=mtm_speed,
    channel=3,
    adds_per_group=adds_per_group["ll"],
    gain=gain["ll"],
    sweeps=sweeps["ll"],
)

# plot cleaned ifgs
plt.imshow(ifg["ll"].T, aspect="auto", extent=[0, len(ifg["ll"]), 0, len(ifg["ll"][0])])
plt.savefig("../../output/plots/cleaned_ifg_over_time.png")
plt.clf()

# check if nans in ifgs
if np.isnan(ifg["ll"]).any():
    print("Nans in interferograms")
    print(ifg["ll"][np.isnan(ifg["ll"])])

print("converting interferograms to spectra")

spec_ss = ifg_to_spec(
    ifg["ll"],
    mtm_speed,
    channel,
    adds_per_group["ll"],
    bol_cmd_bias["ll"],
    bol_volt["ll"],
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
bb_ical = planck(f_ghz[: len(spec_ss)], ical)
# switching the sign for when the ical has higher temp than T_CMB
# bb_ical[ical > T_CMB] = -1 * bb_ical[ical > T_CMB]
ical_emiss_ss = fits_data_ss[1].data["RICAL"][0] + 1j * fits_data_ss[1].data["IICAL"][0]
ical_emiss_ss = ical_emiss_ss[: len(spec_ss)]

# dihedral spectrum
bb_dihedral = planck(f_ghz[: len(spec_ss)], dihedral)
dihedral_emiss_ss = (
    fits_data_ss[1].data["RDIHEDRA"][0] + 1j * fits_data_ss[1].data["IDIHEDRA"][0]
)
dihedral_emiss_ss = dihedral_emiss_ss[: len(spec_ss)]

sky_ss = (
    spec_ss
    - (
        (bb_ical * ical_emiss_ss)[:, : len(spec_ss[0])]
        + (bb_dihedral * dihedral_emiss_ss)[:, : len(spec_ss[0])]
    )
    / otf_ss[np.newaxis, :]
)

print("saving sky")

# save the sky
np.savez("../../output/data/sky.npz", sky=sky_ss, pix_gal=pix_gal)
