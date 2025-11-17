"""
Script to simulate the sky as seen by FIRAS, based on given XCAL, ICAL, dihedral, reference and sky horns, bolometer (4 channels) temperatures.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import globals as g
from calibration import bolometer
from pipeline import ifg_spec
from utils import my_utils as utils
from utils.config import gen_nyquistl

modes = {"ss": 0, "lf": 3}
# modes = {"ss": 0}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}
# channels = {"ll": 3}

# temps = {
#     "xcal": np.array([2.70413828]),
#     "ical": np.array([2.71052694]),
#     "dihedral": np.array([5.0066607]),
#     "refhorn": np.array([2.6955471]),
#     "skyhorn": np.array([2.69618714]),
#     "collimator": np.array([2.69618714]),
#     "bolometer_ll": np.array([1.55804986]),
#     "bolometer_lh": np.array([1.55815428]),
#     "bolometer_rl": np.array([1.54744506]),
#     "bolometer_rh": np.array([1.54760718]),
# }
temps = [
    2.70413828,
    2.71052694,
    5.0066607,
    2.6955471,
    2.69618714,
    2.69618714,
    1.55804986,
]


def generate_ifg(
    channel,
    mode,
    temps,
    real_xcal_spec,
    apod,
    adds_per_group=np.array([3]),
    sweeps=np.array([16]),
    bol_cmd_bias=np.array([35]),
    bol_volt=np.array([2.048020362854004]),
    gain=np.array([300]),
):
    frequency = utils.generate_frequencies(channel, mode, 257)

    bb_xcal = utils.planck(frequency, np.array(temps[0]))
    bb_ical = utils.planck(frequency, np.array(temps[1]))
    bb_dihedral = utils.planck(frequency, np.array(temps[2]))
    bb_refhorn = utils.planck(frequency, np.array(temps[3]))
    bb_skyhorn = utils.planck(frequency, np.array(temps[4]))
    bb_collimator = utils.planck(frequency, np.array(temps[5]))
    bb_bolometer = utils.planck(frequency, np.array(temps[6]))
    # bb_bolometer_lh = utils.planck(frequency, np.array(temps[7]))
    # bb_bolometer_rl = utils.planck(frequency, np.array(temps[8]))
    # bb_bolometer_rh = utils.planck(frequency, np.array(temps[9]))
    fits_data = fits.open(
        f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
    )

    (
        emiss_xcal,
        emiss_ical,
        emiss_dihedral,
        emiss_refhorn,
        emiss_skyhorn,
        emiss_bolometer,
        emiss_collimator,
    ) = np.zeros((7, 257), dtype=np.complex128)
    print(f"Processing {channel.upper()}{mode.upper()}...")

    mtm_speed = 0 if mode[1] == "s" else 1
    if mtm_speed == 0:
        cutoff = 5
    else:
        cutoff = 7
    length = len(fits_data[1].data["RTRANSFE"][0])
    emiss_xcal[cutoff : cutoff + length] = (
        fits_data[1].data["RTRANSFE"][0] + 1j * fits_data[1].data["ITRANSFE"][0]
    )
    emiss_ical[cutoff : cutoff + length] = (
        fits_data[1].data["RICAL"][0] + 1j * fits_data[1].data["IICAL"][0]
    )
    emiss_dihedral[cutoff : cutoff + length] = (
        fits_data[1].data["RDIHEDRA"][0] + 1j * fits_data[1].data["IDIHEDRA"][0]
    )
    emiss_refhorn[cutoff : cutoff + length] = (
        fits_data[1].data["RREFHORN"][0] + 1j * fits_data[1].data["IREFHORN"][0]
    )
    emiss_skyhorn[cutoff : cutoff + length] = (
        fits_data[1].data["RSKYHORN"][0] + 1j * fits_data[1].data["ISKYHORN"][0]
    )
    emiss_collimator[cutoff : cutoff + length] = (
        fits_data[1].data["RSTRUCTU"][0] + 1j * fits_data[1].data["ISTRUCTU"][0]
    )
    emiss_bolometer[cutoff : cutoff + length] = (
        fits_data[1].data["RBOLOMET"][0] + 1j * fits_data[1].data["IBOLOMET"][0]
    )

    total_spectra = np.nan_to_num(
        bb_xcal
        + (
            bb_ical * emiss_ical
            + bb_dihedral * emiss_dihedral
            + bb_refhorn * emiss_refhorn
            + bb_skyhorn * emiss_skyhorn
            + bb_collimator * emiss_collimator
            + bb_bolometer * emiss_bolometer
            # + bb_bolometer_lh * emiss_bolometer
            # + bb_bolometer_rl * emiss_bolometer
            # + bb_bolometer_rh * emiss_bolometer
        )
        / emiss_xcal,
        nan=0,
    )

    # noised_spec = total_spectra + (real_xcal_spec - bb_xcal)
    # plt.plot(noised_spec[0], label="Noised Spectrum")
    # plt.plot(total_spectra[0], label="Total Spectrum")
    # plt.legend()
    # plt.title(f"{channel.upper()}{mode.upper()} Spectra Comparison")
    # plt.xlabel("Frequency (GHz)")
    # plt.ylabel("Brightness Temperature (K)")
    # plt.grid()
    # plt.show()

    fnyq = gen_nyquistl(
        "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
    )
    frec = 4 * (channels[channel] % 2) + modes[mode]

    ifg = ifg_spec.spec_to_ifg(
        spec=total_spectra,
        # spec=noised_spec,
        channel=channel,
        mode=mode,
        adds_per_group=adds_per_group,
        bol_cmd_bias=bol_cmd_bias / 25.5,  # convert to volts
        bol_volt=bol_volt,
        Tbol=temps[6],
        gain=gain,
        sweeps=sweeps,
        apod=apod,
        otf=emiss_xcal,
        fnyq_icm=fnyq["icm"][frec],
    )

    # plt.plot(ifg[0], label=f"{channel.upper()}{mode.upper()} IFG")
    # plt.title(f"{channel.upper()}{mode.upper()} IFG")
    # plt.xlabel("Sample")
    # plt.ylabel("Amplitude")
    # plt.legend()
    # plt.grid()
    # plt.show()
    print(f"IFG for {channel.upper()}{mode.upper()} generated.")

    return ifg, total_spectra, bb_xcal
    # return ifg, noised_spec, bb_xcal


if __name__ == "__main__":
    for channel in channels.keys():
        for mode in modes.keys():
            if not (mode == "lf" and (channel == "lh" or channel == "rh")):
                ifg, total_spectra, xcal_spec = generate_ifg(channel, mode, temps)

                # save this ifg to a file
                np.save("/simulations/output/ifgsim.npy", ifg)

                plt.plot((ifg[0]), label=f"{channel.upper()}{mode.upper()} IFG")
                plt.title(f"{channel.upper()}{mode.upper()} IFG")
                plt.xlabel("Sample")
                plt.ylabel("Amplitude")
                plt.legend()
                plt.grid()
                plt.show()

                plt.plot(
                    np.abs(total_spectra[0]),
                    label=f"{channel.upper()}{mode.upper()} Total Spectra",
                )
                plt.title(f"{channel.upper()}{mode.upper()} Total Spectra")
                plt.xlabel("Frequency (GHz)")
                plt.ylabel("Brightness Temperature (K)")
                plt.legend()
                plt.grid()
                plt.show()
