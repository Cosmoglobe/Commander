"""
This script takes the interferograms into spectra.
"""

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
# modes = ["ss", "sf", "ls", "lf", "fs", "fl"]
modes = ["ss"]

sky_data = h5py.File(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v4.1.h5",
    # "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v1.h5",
    "r",
)
# cal_data = h5py.File(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/cal_v1.h5",
#     "r",
# )

# print("parsing gmts")
# # getting the xcal data into their own variables
# gmts = np.array(sky_data["df_data/gmt"][()]).astype(str)
# gmt = pd.to_datetime(gmts).strftime("%Y-%m-%d %H:%M:%S").to_numpy()

print("other variables")
# ical = np.array(sky_data["df_data/ical"][()])
a_ical = np.array(sky_data["df_data/a_ical"][()])
b_ical = np.array(sky_data["df_data/b_ical"][()])
# dihedral = np.array(sky_data["df_data/dihedral"][()])
a_dihedral = np.array(sky_data["df_data/a_dihedral"][()])
b_dihedral = np.array(sky_data["df_data/b_dihedral"][()])
mtm_length = np.array(sky_data["df_data/mtm_length"][()])
mtm_speed = np.array(sky_data["df_data/mtm_speed"][()])

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

print(f"data loaded: {len(ifg['ll'])}")

fits_data = fits.open(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/reference/FIRAS_CALIBRATION_MODEL_LLSS.FITS"
)
fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

scan_mode = 0  # SS
channel = 3  # LL
frec = 4 * (channel % 2) + scan_mode

# optical transfer function
otf = fits_data[1].data["RTRANSFE"][0] + 1j * fits_data[1].data["ITRANSFE"][0]
otf = otf[np.abs(otf) > 0]

# bolometer function
tau = fits_data[1].data["TIME_CON"][0]
Jo = fits_data[1].data["BOLPARM8"][0]
Jg = fits_data[1].data["BOLPARM9"][0]
Tbol = fits_data[1].data["BOLOM_B2"][0]
T0 = fits_data[1].data["BOLPARM2"][0]
R0 = fits_data[1].data["BOLPARM_"][0]
rho = fits_data[1].data["BOLPARM5"][0]
G1 = fits_data[1].data["BOLPARM3"][0]
beta = fits_data[1].data["BOLPARM4"][0]

# plot ifgs over time
plt.imshow(ifg["ll"].T, aspect="auto", extent=[0, len(ifg["ll"]), 0, len(ifg["ll"][0])])
plt.savefig("../../output/plots/ifg_over_time.png")
plt.clf()

# plot random ifg
plt.plot(ifg["ll"][np.random.randint(0, len(ifg["ll"]))])
plt.savefig("../../output/plots/random_ifg.png")
plt.clf()

print("cleaning interferograms")

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

spec = ifg_to_spec(
    ifg["ll"],
    mtm_speed,
    channel,
    adds_per_group["ll"],
    bol_cmd_bias["ll"],
    bol_volt["ll"],
    fnyq["icm"][frec],
    fnyq["hz"][frec],
    otf,
    Jo,
    Jg,
    Tbol,
    rho,
    R0,
    T0,
    beta,
    G1,
    tau,
)

print("plotting")

plt.imshow(np.abs(spec).T, aspect="auto", extent=[0, len(spec), 0, len(spec[0])])
plt.savefig("../../output/tests/spectra_diff_over_time.png")
plt.clf()

print("making the diff")

# frequency mapping
f_icm = np.arange(210) * (fnyq["icm"][frec] / 320) + 2
c = 3e8 * 1e2  # cm/s
f_ghz = f_icm * c * 1e-9

# ical spectrum
ical = (
    a_ical + b_ical
) / 2  # doing the mean between the two sides for now, will be fitted later
bb_ical = planck(f_ghz[: len(spec)], ical)
# switching the sign for when the ical has higher temp than T_CMB
# bb_ical[ical > T_CMB] = -1 * bb_ical[ical > T_CMB]
ical_emiss = fits_data[1].data["RICAL"][0] + 1j * fits_data[1].data["IICAL"][0]
ical_emiss = ical_emiss[: len(spec)]

# dihedral spectrum
dihedral = (a_dihedral + b_dihedral) / 2  # same for dihedral
bb_dihedral = planck(f_ghz[: len(spec)], dihedral)
dihedral_emiss = (
    fits_data[1].data["RDIHEDRA"][0] + 1j * fits_data[1].data["IDIHEDRA"][0]
)
dihedral_emiss = dihedral_emiss[: len(spec)]

sky = (
    spec
    - (
        (bb_ical * ical_emiss)[:, : len(spec[0])]
        + (bb_dihedral * dihedral_emiss)[:, : len(spec[0])]
    )
    / otf[np.newaxis, :]
)

print("saving sky")

# save the sky
np.save("../../output/data/sky.npy", sky)

print("plotting sky")

# plot the sky
plt.imshow(
    np.abs(sky).T, aspect="auto", extent=[0, len(sky), 0, len(sky[0])], vmax=500, vmin=0
)
plt.savefig("../../output/plots/sky_over_time.png")
plt.clf()
