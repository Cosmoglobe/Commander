"""
This script takes the interferograms into spectra, fits for the monopole temperature
and makes the temperature map.
"""

from datetime import datetime

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from my_utils import clean_ifg, ifg_to_spec, planck, residuals
from scipy.optimize import minimize
from utils.config import gen_nyquistl

data = h5py.File(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/df_v14.h5",
    "r",
)

# indices to run for
# idx_xcal_in = np.where(np.array(data["df_data/xcal_pos"][()]) == 1)
idx = np.where(np.array(data["df_data/xcal_pos"][()]) == 2)

print("parsing gmts")
# getting the xcal data into their own variables
gmts = np.array(data["df_data/gmt"][()])[idx].astype(str)
gmt = pd.to_datetime(gmts).strftime("%Y-%m-%d %H:%M:%S").to_numpy()

print("other variables")
ical_xcal_in = np.array(data["df_data/ical"][()])[idx]
xcal = np.array(data["df_data/xcal"][()])[idx]
ifg_ll = np.array(data["df_data/ifg_ll"][()])[idx]
ifg_lh = np.array(data["df_data/ifg_lh"][()])[idx]
ifg_rl = np.array(data["df_data/ifg_rl"][()])[idx]
ifg_rh = np.array(data["df_data/ifg_rh"][()])[idx]
mtm_length = np.array(data["df_data/mtm_length"][()])[idx]
mtm_speed = np.array(data["df_data/mtm_speed"][()])[idx]
bol_cmd_bias_ll = np.array(data["df_data/bol_cmd_bias_ll"][()])[idx]
bol_volt_ll = np.array(data["df_data/bol_volt_ll"][()])[idx]
times = np.array(data["df_data/time"][()])[idx]
adds_per_group_ll = np.array(data["df_data/adds_per_group_ll"][()])[idx]
gain_ll = np.array(data["df_data/gain_ll"][()])[idx]
sweeps_ll = np.array(data["df_data/sweeps_ll"][()])[idx]

data.close()

print("data loaded")

# find the record with the smallest difference between the two temperatures
# least_diff_idx = np.argmin(np.abs(ical_xcal_in - xcal))

fits_data = fits.open("./../../reference/FIRAS_CALIBRATION_MODEL_LLSS.FITS")
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

# spec = clean_ifg(
#     ifg=ifg_ll[least_diff_idx],
#     mtm_length=mtm_length[least_diff_idx],
#     mtm_speed=mtm_speed[least_diff_idx],
#     channel=3,
#     adds_per_group=adds_per_group[least_diff_idx],
#     bol_cmd_bias=bol_cmd_bias_ll[least_diff_idx],
#     bol_volt=bol_volt_ll[least_diff_idx],
#     fits_data=fits_data,
#     fnyq_icm=fnyq["icm"][frec],
#     fnyq_hz=fnyq["hz"][frec],
#     otf=otf,
# )

# clean_ifg = np.vectorize(clean_ifg, excluded=['channel', 'fits_data', 'fnyq', 'frec', 'otf'])

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
t_ical = ical_xcal_in
bb_ical = planck(f_ghz[: len(spec)], t_ical)

ical_emiss = fits_data[1].data["RICAL"][0] + 1j * fits_data[1].data["IICAL"][0]
ical_emiss = ical_emiss[: len(spec)]

sky = spec - (bb_ical * ical_emiss)[:, : len(spec[0])] / otf[np.newaxis, :]

print("saving sky")

# save the sky
np.save("../../output/data/sky.npy", sky)

print("plotting sky")

# plot the sky
plt.imshow(
    np.abs(sky).T, aspect="auto", extent=[0, len(sky), 0, len(sky[0])], vmax=500, vmin=0
)
plt.savefig("../../output/tests/sky_over_time.png")
plt.clf()

for i in range(0, len(sky), 1000):
    plt.plot(f_ghz[1 : len(spec[0])], np.abs(sky[i, 1:]))
    plt.title(f"Sky_{i:03d}")
    plt.savefig(f"../../output/sky/sky_{i:03d}.png")
    plt.clf()

quit()

t0 = 2.7
fit = minimize(residuals, t0, args=(f_ghz[1 : len(spec[0])], np.abs(sky[:, 1:])))

print(f"fit: {fit}")

fit_temp = fit.x[0]

print(f"fit_temp shape: {fit_temp.shape}")

bb_xcal = planck(f_ghz[1 : len(spec)], xcal)

print(f"Temperature fit: {fit_temp}, actual temperature: {xcal}")
plt.plot(f_ghz[1 : len(spec[0])], np.abs(sky[1:]), label="Sky")
plt.plot(f_ghz[1 : len(spec[0])], np.abs(bb_ical * ical_emiss / otf)[1:], label="ICAL")
plt.plot(f_ghz[1 : len(spec[0])], planck(f_ghz[1 : len(spec)], fit_temp), label="Fit")
plt.plot(f_ghz[1 : len(spec[0])], bb_xcal, label="XCAL")
plt.plot(f_ghz[1 : len(spec[0])], spec[1:], label="Difference (spectrum)")
plt.plot(
    f_ghz[1 : len(spec[0])], np.abs(bb_xcal - bb_ical[1:]), label="Difference (BBs)"
)
plt.legend()
plt.savefig("../../output/tests/fit.png")
