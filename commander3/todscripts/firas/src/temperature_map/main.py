"""
This script takes the interferograms into spectra, fits for the monopole temperature
and makes the temperature map.
"""

import h5py
import numpy as np

from astropy.io import fits

import matplotlib.pyplot as plt

from my_utils import clean_ifg, planck, residuals

from utils.config import gen_nyquistl

from scipy.optimize import minimize

data = h5py.File(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/df_v12.h5",
    "r",
)

# indices of the data where xcal is in
idx_xcal_in = np.where(data["df_data/xcal_pos"][()] == 1)

# getting the xcal data into their own variables
ical_xcal_in = data["df_data/ical"][idx_xcal_in]
xcal = data["df_data/xcal"][idx_xcal_in]
ifg_ll = data["df_data/ifg_ll"][idx_xcal_in]
ifg_lh = data["df_data/ifg_lh"][idx_xcal_in]
ifg_rl = data["df_data/ifg_rl"][idx_xcal_in]
ifg_rh = data["df_data/ifg_rh"][idx_xcal_in]
mtm_length = data["df_data/mtm_length"][idx_xcal_in]
mtm_speed = data["df_data/mtm_speed"][idx_xcal_in]
bol_cmd_bias_ll = data["df_data/bol_cmd_bias_ll"][idx_xcal_in]
bol_volt_ll = data["df_data/bol_volt_ll"][idx_xcal_in]
times = data["df_data/time"][idx_xcal_in]

data.close()

# find the record with the smallest difference between the two temperatures
least_diff_idx = np.argmin(np.abs(ical_xcal_in - xcal))

fits_data = fits.open('./../../reference/FIRAS_CALIBRATION_MODEL_LLSS.FITS')
fnyq = gen_nyquistl('../../reference/fex_samprate.txt', '../../reference/fex_nyquist.txt', 'int')

scan_mode = 0  # SS
channel = 3 # LL
frec = 4 * (channel % 2) + scan_mode

spec = clean_ifg(
    ifg_ll[least_diff_idx],
    mtm_length=mtm_length[least_diff_idx],
    mtm_speed=mtm_speed[least_diff_idx],
    channel=3,
    adds_per_group=1,
    bol_cmd_bias=bol_cmd_bias_ll[least_diff_idx],
    bol_volt=bol_volt_ll[least_diff_idx],
    fits_data=fits_data,
    fnyq=fnyq,
    frec=frec,
)

# optical transfer function
otf = fits_data[1].data["RTRANSFE"][0] + 1j * fits_data[1].data["ITRANSFE"][0]
otf = otf[np.abs(otf) > 0]
spec = spec[:len(otf)] / otf

# frequency mapping
f_icm = np.arange(210)*(fnyq['icm'][frec]/320) + 2
c = 3e8 * 1e2 # cm/s
f_ghz = f_icm * c * 1e-9

# ical spectrum
t_ical = ical_xcal_in[least_diff_idx]
bb_ical = planck(f_ghz[:len(spec)], t_ical)

ical_emiss = fits_data[1].data["RICAL"][0] + 1j * fits_data[1].data["IICAL"][0]
ical_emiss = ical_emiss[:len(spec)]

sky = spec - bb_ical * ical_emiss / otf

t0 = 2.7
fit = minimize(residuals, t0, args=(f_ghz[1:len(spec)], np.abs(sky[1:])))

fit_temp = fit.x[0]

bb_xcal = planck(f_ghz[1:len(spec)], xcal[least_diff_idx])

print(f"Temperature fit: {fit_temp}, actual temperature: {xcal[least_diff_idx]}")
# plt.plot(f_ghz[1:len(spec)], np.abs(sky[1:]), label="Sky")
# plt.plot(f_ghz[1:len(spec)], np.abs(bb_ical * ical_emiss / otf)[1:], label="ICAL")
# plt.plot(f_ghz[1:len(spec)], planck(f_ghz[1:len(spec)], fit_temp), label="Fit")
# plt.plot(f_ghz[1:len(spec)], bb_xcal, label="XCAL")
plt.plot(f_ghz[1:len(spec)], spec[1:], label="Difference (spectrum)")
plt.plot(f_ghz[1:len(spec)], np.abs(bb_xcal - bb_ical[1:]), label="Difference (BBs)")
plt.legend()
plt.savefig("../../output/tests/fit.png")
