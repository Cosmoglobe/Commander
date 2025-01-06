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

sky_data = h5py.File(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v1.h5",
    "r",
)
# cal_data = h5py.File(
#     "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/cal_v1.h5",
#     "r",
# )

print("parsing gmts")
# getting the xcal data into their own variables
gmts = np.array(sky_data["df_data/gmt"][()]).astype(str)
gmt = pd.to_datetime(gmts).strftime("%Y-%m-%d %H:%M:%S").to_numpy()

print("other variables")
ical = np.array(sky_data["df_data/ical"][()])
dihedral = np.array(sky_data["df_data/dihedral"][()])
ifg_ll = np.array(sky_data["df_data/ifg_ll"][()])
ifg_lh = np.array(sky_data["df_data/ifg_lh"][()])
ifg_rl = np.array(sky_data["df_data/ifg_rl"][()])
ifg_rh = np.array(sky_data["df_data/ifg_rh"][()])
mtm_length = np.array(sky_data["df_data/mtm_length"][()])
mtm_speed = np.array(sky_data["df_data/mtm_speed"][()])
bol_cmd_bias_ll = np.array(sky_data["df_data/bol_cmd_bias_ll"][()])
bol_volt_ll = np.array(sky_data["df_data/bol_volt_ll"][()])
times = np.array(sky_data["df_data/time"][()])
adds_per_group_ll = np.array(sky_data["df_data/adds_per_group_ll"][()])
gain_ll = np.array(sky_data["df_data/gain_ll"][()])
sweeps_ll = np.array(sky_data["df_data/sweeps_ll"][()])

sky_data.close()

print("data loaded")

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

print("cleaning interferograms")

ifg = clean_ifg(
    ifg=ifg_ll,
    mtm_length=mtm_length,
    mtm_speed=mtm_speed,
    channel=3,
    adds_per_group=adds_per_group_ll,
    gain=gain_ll,
    sweeps=sweeps_ll,
)

# check if nans in ifgs
if np.isnan(ifg).any():
    print("Nans in interferograms")
    print(ifg[np.isnan(ifg)])

print("converting interferograms to spectra")

spec = ifg_to_spec(
    ifg,
    mtm_speed,
    channel,
    adds_per_group_ll,
    bol_cmd_bias_ll,
    bol_volt_ll,
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
bb_ical = planck(f_ghz[: len(spec)], ical)
ical_emiss = fits_data[1].data["RICAL"][0] + 1j * fits_data[1].data["IICAL"][0]
ical_emiss = ical_emiss[: len(spec)]

# dihedral spectrum
bb_dihedral = planck(f_ghz[: len(spec)], dihedral)
dihedral_emiss = (
    fits_data[1].data["RDIHEDRA"][0] + 1j * fits_data[1].data["IDIHEDRA"][0]
)
dihedral_emiss = dihedral_emiss[: len(spec)]

sky = (
    spec
    - (bb_ical * ical_emiss)[:, : len(spec[0])] / otf[np.newaxis, :]
    - (bb_dihedral * dihedral_emiss)[:, : len(spec[0])] / otf[np.newaxis, :]
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

for i in range(0, len(sky), 1000):
    plt.plot(f_ghz[1 : len(spec[0])], np.abs(sky[i, 1:]))
    plt.title(f"Sky_{i:03d}")
    plt.savefig(f"../../output/sky/sky_{i:03d}.png")
    plt.clf()
