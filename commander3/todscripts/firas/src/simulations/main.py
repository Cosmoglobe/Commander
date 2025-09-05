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

temps = {
    "xcal": np.array([2.70413828]),
    "ical": np.array([2.71052694]),
    "dihedral": np.array([5.0066607]),
    "refhorn": np.array([2.6955471]),
    "skyhorn": np.array([2.69618714]),
    "collimator": np.array([2.69618714]),
    "bolometer_ll": np.array([1.55804986]),
    "bolometer_lh": np.array([1.55815428]),
    "bolometer_rl": np.array([1.54744506]),
    "bolometer_rh": np.array([1.54760718]),
}


def generate_ifg(
    channel,
    mode,
    temps,
    adds_per_group=np.array([3]),
    sweeps=np.array([16]),
    bol_cmd_bias=np.array([35]),
    bol_volt=np.array([2.048020362854004]),
    gain=np.array([300]),
):
    frequency = utils.generate_frequencies(channel, mode, 257)

    bb_xcal = utils.planck(frequency, np.array(temps["xcal"]))
    bb_ical = utils.planck(frequency, np.array(temps["ical"]))
    bb_dihedral = utils.planck(frequency, np.array(temps["dihedral"]))
    bb_refhorn = utils.planck(frequency, np.array(temps["refhorn"]))
    bb_skyhorn = utils.planck(frequency, np.array(temps["skyhorn"]))
    bb_collimator = utils.planck(frequency, np.array(temps["collimator"]))
    bb_bolometer_ll = utils.planck(frequency, np.array(temps["bolometer_ll"]))
    bb_bolometer_lh = utils.planck(frequency, np.array(temps["bolometer_lh"]))
    bb_bolometer_rl = utils.planck(frequency, np.array(temps["bolometer_rl"]))
    bb_bolometer_rh = utils.planck(frequency, np.array(temps["bolometer_rh"]))

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
        (
            bb_xcal * emiss_xcal
            + bb_ical * emiss_ical
            + bb_dihedral * emiss_dihedral
            + bb_refhorn * emiss_refhorn
            + bb_skyhorn * emiss_skyhorn
            + bb_collimator * emiss_collimator
            + bb_bolometer_ll * emiss_bolometer
            + bb_bolometer_lh * emiss_bolometer
            + bb_bolometer_rl * emiss_bolometer
            + bb_bolometer_rh * emiss_bolometer
        )
        / emiss_xcal,
        nan=0,
    )

    fnyq = gen_nyquistl(
        "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
    )
    frec = 4 * (channels[channel] % 2) + modes[mode]

    # apod = fits_data[1].data["APODIZAT"][0]
    apod = np.ones(512, dtype=np.float64)  # No apodization for now
    R0, T0, G1, beta, rho, C1, C3, Jo, Jg = bolometer.get_bolometer_parameters(
        channel, mode
    )

    ifg = ifg_spec.spec_to_ifg(
        spec=total_spectra,
        channel=channel,
        mode=mode,
        adds_per_group=adds_per_group,
        bol_cmd_bias=bol_cmd_bias / 25.5,  # convert to volts
        bol_volt=bol_volt,
        Tbol=temps[f"bolometer_{channel}"],
        gain=gain,
        sweeps=sweeps,
        apod=apod,
        otf=emiss_xcal,
        fnyq_icm=fnyq["icm"][frec],
    )
    print(f"IFG for {channel.upper()}{mode.upper()} generated.")

    return ifg, total_spectra


if __name__ == "__main__":
    for channel in channels.keys():
        for mode in modes.keys():
            if not (mode == "lf" and (channel == "lh" or channel == "rh")):
                ifg, total_spectra = generate_ifg(channel, mode, temps)

                # save this ifg to a file
                np.save("./ifgsim.npy", ifg)

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
