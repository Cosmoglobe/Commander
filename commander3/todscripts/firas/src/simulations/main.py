"""
Script to simulate the sky as seen by FIRAS, based on given XCAL, ICAL, dihedral, reference and sky horns, bolometer (4 channels) temperatures.
"""
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
from utils.config import gen_nyquistl

modes = {"ss": 0, "lf": 3}
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}

temps = {
    "xcal": g.T_CMB,
    "ical": 2.76,
    "dihedral": 2.00,
    "refhorn": 2.72,
    "skyhorn": 2.72,
    "bolometer_ll": 1.53,
    "bolometer_lh": 1.53,
    "bolometer_rl": 1.53,
    "bolometer_rh": 1.53,
}

def generate_ifg(channel, mode, adds_per_group=1, bol_cmd_bias=1, bol_volt=1):
    frequency = mu.generate_frequencies(channel, mode)

    bb_xcal = mu.planck(frequency, np.array(temps["xcal"]))
    bb_ical = mu.planck(frequency, np.array(temps["ical"]))
    bb_dihedral = mu.planck(frequency, np.array(temps["dihedral"]))
    bb_refhorn = mu.planck(frequency, np.array(temps["refhorn"]))
    bb_skyhorn = mu.planck(frequency, np.array(temps["skyhorn"]))
    bb_bolometer_ll = mu.planck(frequency, np.array(temps["bolometer_ll"]))
    bb_bolometer_lh = mu.planck(frequency, np.array(temps["bolometer_lh"]))
    bb_bolometer_rl = mu.planck(frequency, np.array(temps["bolometer_rl"]))
    bb_bolometer_rh = mu.planck(frequency, np.array(temps["bolometer_rh"]))

    fits_data = fits.open(
        f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
    )

    emiss_xcal = fits_data[1].data["RTRANSFE"][0]+ 1j * fits_data[1].data["ITRANSFE"][0]
    emiss_ical = fits_data[1].data["RICAL"][0] + 1j * fits_data[1].data["IICAL"][0]
    emiss_dihedral = fits_data[1].data["RDIHEDRA"][0] + 1j * fits_data[1].data["IDIHEDRA"][0]
    emiss_refhorn = fits_data[1].data["RREFHORN"][0] + 1j * fits_data[1].data["IREFHORN"][0]
    emiss_skyhorn = fits_data[1].data["RSKYHORN"][0] + 1j * fits_data[1].data["ISKYHORN"][0]
    emiss_bolometer = fits_data[1].data["RBOLOMET"][0] + 1j * fits_data[1].data["IBOLOMET"][0]

    total_spectra = (bb_xcal * emiss_xcal+ bb_ical * emiss_ical+ bb_dihedral * emiss_dihedral + bb_refhorn * emiss_refhorn + bb_skyhorn * emiss_skyhorn + bb_bolometer_ll * emiss_bolometer + bb_bolometer_lh * emiss_bolometer + bb_bolometer_rl * emiss_bolometer + bb_bolometer_rh * emiss_bolometer)/emiss_xcal
    plt.plot(frequency, emiss_xcal/emiss_xcal, label="xcal")
    plt.plot(frequency, emiss_ical/emiss_xcal, label="ical")
    plt.plot(frequency, emiss_dihedral/emiss_xcal, label="dihedral")
    plt.plot(frequency, emiss_refhorn/emiss_xcal, label="refhorn")
    plt.plot(frequency, emiss_skyhorn/emiss_xcal, label="skyhorn")
    plt.plot(frequency, emiss_bolometer/emiss_xcal, label="bolometer")
    plt.legend()
    plt.show()
    plt.plot(
        frequency,
        np.abs(total_spectra),
        label=f"{channel.upper()}{mode.upper()}",
    )
    plt.show()

    fnyq = gen_nyquistl(
        "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
    )
    frec = 4 * (channels[channel] % 2) + modes[mode]

    apod = fits_data[1].data["APODIZAT"][0]
    R0 = fits_data[1].data["BOLPARAM_"][0]
    T0 = fits_data[1].data["BOLPARAM2"][0]
    G1 = fits_data[1].data["BOLPARAM3"][0]
    beta = fits_data[1].data["BOLPARAM4"][0]
    rho = fits_data[1].data["BOLPARAM5"][0]
    C1 = fits_data[1].data["BOLPARAM6"][0]
    C3 = fits_data[1].data["BOLPARAM7"][0]
    Jo = fits_data[1].data["JO"][0]
    Jg = fits_data[1].data["JG"][0]

    ifg = mu.spec_to_ifg(spec=total_spectra, mtm_speed=0 if mode[1] == "s" else 1, channel=channel, adds_per_group=adds_per_group, bol_cmd_bias=bol_cmd_bias, bol_volt=bol_volt, Tbol=temps[f"bolometer_{channel}"], gain=gain, sweeps=sweeps, otf=emiss_xcal, apod=apod,fnyq=fnyq["icm"][frec[f"{channel}_{mode}"]], R0=R0, T0=T0, G1=G1, beta=beta, rho=rho, C1=C1, C3=C3, Jo=Jo, Jg=Jg)