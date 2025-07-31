import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

import globals as g
import my_utils as mu

# TODO: generalise for all channels and modes

data = np.load(g.PROCESSED_DATA_PATH_CAL)

ifgs_ll = data["ifg_ll_ss"]

n = np.random.randint(0, ifgs_ll.shape[0])

temp_xcal = data["xcal_ss"]
temp_ical = data["ical_ss"]
temp_dihedral = data["dihedral_ss"]
temp_refhorn = data["refhorn_ss"]
temp_skyhorn = data["skyhorn_ss"]
temp_collimator = data["collimator_ss"]
temp_bolometer_ll = data["bolometer_ll_ss"]
temp_bolometer_lh = data["bolometer_lh_ss"]
temp_bolometer_rl = data["bolometer_rl_ss"]
temp_bolometer_rh = data["bolometer_rh_ss"]

frequencies_ll = mu.generate_frequencies("ll", "ss")

bb_xcal = mu.planck(frequencies_ll, np.array(temp_xcal))
bb_ical = mu.planck(frequencies_ll, np.array(temp_ical))
bb_dihedral = mu.planck(frequencies_ll, np.array(temp_dihedral))
bb_refhorn = mu.planck(frequencies_ll, np.array(temp_refhorn))
bb_skyhorn = mu.planck(frequencies_ll, np.array(temp_skyhorn))
bb_collimator = mu.planck(frequencies_ll, np.array(temp_collimator))
bb_bolometer_ll = mu.planck(frequencies_ll, np.array(temp_bolometer_ll))
bb_bolometer_lh = mu.planck(frequencies_ll, np.array(temp_bolometer_lh))
bb_bolometer_rl = mu.planck(frequencies_ll, np.array(temp_bolometer_rl))
bb_bolometer_rh = mu.planck(frequencies_ll, np.array(temp_bolometer_rh))

adds_per_group = data["adds_per_group_ss"]
sweeps = data["sweeps_ss"]

fits_data = fits.open(f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_LLSS.FITS")
apod = fits_data[1].data["APODIZAT"][0]

emiss_xcal = fits_data[1].data["RTRANSFE"][0][:43] + 1j * fits_data[1].data["ITRANSFE"][0][:43]
emiss_ical = fits_data[1].data["RICAL"][0][:43] + 1j * fits_data[1].data["IICAL"][0][:43]
emiss_dihedral = fits_data[1].data["RDIHEDRA"][0][:43] + 1j * fits_data[1].data["IDIHEDRA"][0][:43]
emiss_refhorn = fits_data[1].data["RREFHORN"][0][:43] + 1j * fits_data[1].data["IREFHORN"][0][:43]
emiss_skyhorn = fits_data[1].data["RSKYHORN"][0][:43] + 1j * fits_data[1].data["ISKYHORN"][0][:43]
emiss_bolometer = fits_data[1].data["RBOLOMET"][0][:43] + 1j * fits_data[1].data["IBOLOMET"][0][:43]
emiss_collimator = fits_data[1].data["RSTRUCTU"][0][:43] + 1j * fits_data[1].data["ISTRUCTU"][0][:43]

spec = bb_xcal + (bb_ical * emiss_ical + bb_dihedral * emiss_dihedral + bb_refhorn * emiss_refhorn + bb_skyhorn * emiss_skyhorn + bb_bolometer_ll * emiss_bolometer + bb_bolometer_lh * emiss_bolometer + bb_bolometer_rl * emiss_bolometer + bb_bolometer_rh * emiss_bolometer + bb_collimator * emiss_collimator) / emiss_xcal

plt.plot(frequencies_ll, np.abs(spec[n]))
plt.xlabel("Frequency (GHz)")
plt.ylabel("Intensity (MJy/sr)")
# plt.show()

