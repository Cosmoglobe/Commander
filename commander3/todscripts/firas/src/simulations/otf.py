import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import globals as g
import utils.my_utils as utils

# TODO: generalize for all channels and modes
channel = "ll"
mode = "ss"

fits_data = fits.open(
    f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
)

emiss_xcal = (
    fits_data[1].data["RTRANSFE"][0][:43] + 1j * fits_data[1].data["ITRANSFE"][0][:43]
)
emiss_ical = (
    fits_data[1].data["RICAL"][0][:43] + 1j * fits_data[1].data["IICAL"][0][:43]
)
emiss_dihedral = (
    fits_data[1].data["RDIHEDRA"][0][:43] + 1j * fits_data[1].data["IDIHEDRA"][0][:43]
)
emiss_refhorn = (
    fits_data[1].data["RREFHORN"][0][:43] + 1j * fits_data[1].data["IREFHORN"][0][:43]
)
emiss_skyhorn = (
    fits_data[1].data["RSKYHORN"][0][:43] + 1j * fits_data[1].data["ISKYHORN"][0][:43]
)
emiss_bolometer = (
    fits_data[1].data["RBOLOMET"][0][:43] + 1j * fits_data[1].data["IBOLOMET"][0][:43]
)
emiss_collimator = (
    fits_data[1].data["RSTRUCTU"][0][:43] + 1j * fits_data[1].data["ISTRUCTU"][0][:43]
)

plt.plot(emiss_xcal, label="Xcal")
plt.plot(emiss_ical, label="Ical")
plt.plot(emiss_dihedral, label="Dihedral")
plt.plot(emiss_refhorn, label="Refhorn")
plt.plot(emiss_skyhorn, label="Skyhorn")
plt.plot(emiss_bolometer, label="Bolometer")
plt.plot(emiss_collimator, label="Collimator")
plt.legend()
plt.title(f"Emissivities for {channel.upper()}{mode.upper()}")
plt.xlabel("Frequency Index")
plt.ylabel("Emissivity")
plt.show()

frequencies_ll = utils.generate_frequencies(channel, mode)
print(frequencies_ll.shape)
frequencies = utils.generate_frequencies(channel, mode, 257)

plt.plot(frequencies_ll, np.abs(emiss_xcal), label="Xcal")
plt.plot(frequencies_ll, np.abs(emiss_ical), label="Ical")
plt.plot(frequencies_ll, np.abs(emiss_dihedral), label="Dihedral")
plt.plot(frequencies_ll, np.abs(emiss_refhorn), label="Refhorn")
plt.plot(frequencies_ll, np.abs(emiss_skyhorn), label="Skyhorn")
plt.plot(frequencies_ll, np.abs(emiss_bolometer), label="Bolometer")
plt.plot(frequencies_ll, np.abs(emiss_collimator), label="Collimator")
plt.xlim(min(frequencies), max(frequencies))
plt.title(f"Emissivities for {channel.upper()}{mode.upper()}")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Emissivity")
plt.legend()
plt.show()

fig, ax = plt.subplots(2, 1)
ax[0].plot(frequencies_ll, emiss_xcal.real, label="XCAL")
ax[0].plot(frequencies_ll, emiss_ical.real, label="ICAL")
ax[0].plot(frequencies_ll, emiss_dihedral.real, label="Dihedral")
ax[0].plot(frequencies_ll, emiss_refhorn.real, label="Refhorn")
ax[0].plot(frequencies_ll, emiss_skyhorn.real, label="Skyhorn")
ax[0].plot(frequencies_ll, emiss_bolometer.real, label="Bolometer")
ax[0].plot(frequencies_ll, emiss_collimator.real, label="Collimator")
ax[0].set_title(f"Real Emissivities for {channel.upper()}{mode.upper()}")
ax[0].set_xlabel("Frequency (Hz)")
ax[0].set_ylabel("Real Emissivity")
ax[0].legend()
ax[0].set_xlim(min(frequencies), max(frequencies))
ax[1].plot(frequencies_ll, emiss_xcal.imag, label="XCAL")
ax[1].plot(frequencies_ll, emiss_ical.imag, label="ICAL")
ax[1].plot(frequencies_ll, emiss_dihedral.imag, label="Dihedral")
ax[1].plot(frequencies_ll, emiss_refhorn.imag, label="Refhorn")
ax[1].plot(frequencies_ll, emiss_skyhorn.imag, label="Skyhorn")
ax[1].plot(frequencies_ll, emiss_bolometer.imag, label="Bolometer")
ax[1].plot(frequencies_ll, emiss_collimator.imag, label="Collimator")
ax[1].set_title(f"Imaginary Emissivities for {channel.upper()}{mode.upper()}")
ax[1].set_xlabel("Frequency (Hz)")
ax[1].set_ylabel("Imaginary Emissivity")
ax[1].legend()
ax[1].set_xlim(min(frequencies), max(frequencies))
plt.tight_layout()
plt.show()
