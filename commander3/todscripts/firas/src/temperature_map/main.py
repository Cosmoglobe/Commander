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

# indices of the data where xcal is in
idx_xcal_in = np.where(np.array(data["df_data/xcal_pos"][()]) == 1)

print("parsing gmts")
# getting the xcal data into their own variables
gmts = np.array(data["df_data/gmt"][()])[idx_xcal_in].astype(str)
gmt = pd.to_datetime(gmts).strftime("%Y-%m-%d %H:%M:%S").to_numpy()

print("other variables")
ical_xcal_in = np.array(data["df_data/ical"][()])[idx_xcal_in]
xcal = np.array(data["df_data/xcal"][()])[idx_xcal_in]
ifg_ll = np.array(data["df_data/ifg_ll"][()])[idx_xcal_in]
ifg_lh = np.array(data["df_data/ifg_lh"][()])[idx_xcal_in]
ifg_rl = np.array(data["df_data/ifg_rl"][()])[idx_xcal_in]
ifg_rh = np.array(data["df_data/ifg_rh"][()])[idx_xcal_in]
mtm_length = np.array(data["df_data/mtm_length"][()])[idx_xcal_in]
mtm_speed = np.array(data["df_data/mtm_speed"][()])[idx_xcal_in]
bol_cmd_bias_ll = np.array(data["df_data/bol_cmd_bias_ll"][()])[idx_xcal_in]
bol_volt_ll = np.array(data["df_data/bol_volt_ll"][()])[idx_xcal_in]
times = np.array(data["df_data/time"][()])[idx_xcal_in]
adds_per_group_ll = np.array(data["df_data/adds_per_group_ll"][()])[idx_xcal_in]
gain_ll = np.array(data["df_data/gain_ll"][()])[idx_xcal_in]
sweeps_ll = np.array(data["df_data/sweeps_ll"][()])[idx_xcal_in]

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

# check if the spectra have some nans or weird values
if np.isnan(spec).any():
    print("Nans in spectra")
    print(spec[np.isnan(spec)])

print("plotting")

# make a vertical heatmap of each of the 60740 spectra
# each vertical line should be one spectra, and the x axis should go from 0 to len(spec)
# the y axis should go from 0 to len(spec[0])
# the color should be the absolute value of the spectra
plt.imshow(np.abs(spec).T, aspect="auto", extent=[0, len(spec), 0, len(spec[0])])
# plt.imshow(np.abs(spec))
plt.savefig("../../output/tests/spectra_over_time.png")
plt.clf()

for i in range(0, len(spec), 100):
    plt.plot(np.abs(spec[i]))
    plt.savefig(f"../../output/tests/spectrum_{i}.png")
    plt.clf()

quit()

# frequency mapping
f_icm = np.arange(210) * (fnyq["icm"][frec] / 320) + 2
c = 3e8 * 1e2  # cm/s
f_ghz = f_icm * c * 1e-9

# ical spectrum
t_ical = ical_xcal_in
bb_ical = planck(f_ghz[: len(spec)], t_ical)

ical_emiss = fits_data[1].data["RICAL"][0] + 1j * fits_data[1].data["IICAL"][0]
ical_emiss = ical_emiss[: len(spec)]

print(f"spec: {spec.shape}")

sky = spec - (bb_ical * ical_emiss)[:, : len(spec[0])] / otf[np.newaxis, :]

t0 = 2.7
fit = minimize(residuals, t0, args=(f_ghz[1 : len(spec)], np.abs(sky[1:])))

fit_temp = fit.x[0]

bb_xcal = planck(f_ghz[1 : len(spec)], xcal)

print(f"Temperature fit: {fit_temp}, actual temperature: {xcal}")
# plt.plot(f_ghz[1:len(spec)], np.abs(sky[1:]), label="Sky")
# plt.plot(f_ghz[1:len(spec)], np.abs(bb_ical * ical_emiss / otf)[1:], label="ICAL")
# plt.plot(f_ghz[1:len(spec)], planck(f_ghz[1:len(spec)], fit_temp), label="Fit")
# plt.plot(f_ghz[1:len(spec)], bb_xcal, label="XCAL")
plt.plot(f_ghz[1 : len(spec)], spec[1:], label="Difference (spectrum)")
plt.plot(f_ghz[1 : len(spec)], np.abs(bb_xcal - bb_ical[1:]), label="Difference (BBs)")
plt.legend()
plt.savefig("../../output/tests/fit.png")
