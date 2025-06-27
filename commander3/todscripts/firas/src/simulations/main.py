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

# modes = {"ss": 0, "lf": 3}
modes = {"ss": 0}   
# channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}
channels = {"ll": 3}

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

def generate_ifg(channel, mode, temps, adds_per_group=np.array([3]), sweeps=np.array([16]), bol_cmd_bias=np.array([35]), bol_volt=np.array([2.048020362854004]), gain=np.array([300])):
    frequency = mu.generate_frequencies(channel, mode, 257)

    bb_xcal = mu.planck(frequency, np.array(temps["xcal"]))
    bb_ical = mu.planck(frequency, np.array(temps["ical"]))
    bb_dihedral = mu.planck(frequency, np.array(temps["dihedral"]))
    bb_refhorn = mu.planck(frequency, np.array(temps["refhorn"]))
    bb_skyhorn = mu.planck(frequency, np.array(temps["skyhorn"]))
    bb_collimator = mu.planck(frequency, np.array(temps["collimator"]))
    bb_bolometer_ll = mu.planck(frequency, np.array(temps["bolometer_ll"]))
    bb_bolometer_lh = mu.planck(frequency, np.array(temps["bolometer_lh"]))
    bb_bolometer_rl = mu.planck(frequency, np.array(temps["bolometer_rl"]))
    bb_bolometer_rh = mu.planck(frequency, np.array(temps["bolometer_rh"]))

    fits_data = fits.open(
        f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
    )

    emiss_xcal, emiss_ical, emiss_dihedral, emiss_refhorn, emiss_skyhorn, emiss_bolometer, emiss_collimator = np.zeros((7, 257), dtype=np.complex128)
    print(f"Processing {channel.upper()}{mode.upper()}...")
    
    mtm_speed = 0 if mode[1] == "s" else 1
    if mtm_speed == 0:
        cutoff = 5
    else:
        cutoff = 7
    length = len(fits_data[1].data["RTRANSFE"][0])
    emiss_xcal[cutoff:cutoff+length] = fits_data[1].data["RTRANSFE"][0]+ 1j * fits_data[1].data["ITRANSFE"][0]
    emiss_ical[cutoff:cutoff+length] = fits_data[1].data["RICAL"][0] + 1j * fits_data[1].data["IICAL"][0]
    emiss_dihedral[cutoff:cutoff+length] = fits_data[1].data["RDIHEDRA"][0] + 1j * fits_data[1].data["IDIHEDRA"][0]
    emiss_refhorn[cutoff:cutoff+length] = fits_data[1].data["RREFHORN"][0] + 1j * fits_data[1].data["IREFHORN"][0]
    emiss_skyhorn[cutoff:cutoff+length] = fits_data[1].data["RSKYHORN"][0] + 1j * fits_data[1].data["ISKYHORN"][0]
    emiss_collimator[cutoff:cutoff+length] = fits_data[1].data["RSTRUCTU"][0] + 1j * fits_data[1].data["ISTRUCTU"][0]
    emiss_bolometer[cutoff:cutoff+length] = fits_data[1].data["RBOLOMET"][0] + 1j * fits_data[1].data["IBOLOMET"][0]

    total_spectra = np.nan_to_num((bb_xcal * emiss_xcal + bb_ical * emiss_ical + bb_dihedral * emiss_dihedral + bb_refhorn * emiss_refhorn + bb_skyhorn * emiss_skyhorn + bb_collimator * emiss_collimator + bb_bolometer_ll * emiss_bolometer + bb_bolometer_lh * emiss_bolometer + bb_bolometer_rl * emiss_bolometer + bb_bolometer_rh * emiss_bolometer)/emiss_xcal, nan=0)
    #total_spectra = np.nan_to_num((bb_xcal + bb_ical * emiss_ical / emiss_xcal), nan=0)
    # total_spectra = np.nan_to_num((bb_xcal + bb_ical * emiss_ical / emiss_xcal), nan=0)
    
    # plt.plot((emiss_ical / emiss_xcal).real, label="Emissivity ICAL Real")
    # plt.plot((emiss_ical / emiss_xcal).imag, label="Emissivity ICAL Imaginary")
    # plt.plot((emiss_xcal / emiss_xcal).real, label="Emissivity XCAL Real")
    # plt.plot((emiss_xcal / emiss_xcal).imag, label="Emissivity XCAL Imaginary")
    # plt.legend()

    # plt.figure(figsize=(10, 6))
    # plt.plot((-emiss_ical / emiss_xcal).real, label="Emissivity ICAL Real")
    # plt.plot((-emiss_ical / emiss_xcal).imag, label="Emissivity ICAL Imaginary")
    # plt.plot((emiss_xcal / emiss_xcal).real, label="Emissivity XCAL Real")
    # plt.plot((emiss_xcal / emiss_xcal).imag, label="Emissivity XCAL Imaginary")
    # plt.legend()
    # plt.show()

    print(f"Total spectra shape: {total_spectra.shape}")

    # n = np.random.randint(0, len(total_spectra) - 1)
    n = 0

    print(f"total_spectra peak: {np.max(np.abs(total_spectra[n]))}")
    print(f"adds_per_group: {adds_per_group[n]}")
    print(f"sweeps: {sweeps[n]}")
    print(f"bol_cmd_bias: {bol_cmd_bias[n]}")
    print(f"bol_volt: {bol_volt[n]}")
    print(f"gain: {gain[n]}")
    print(f"temps: {temps}")
    
    # plt.plot(frequency, np.abs(total_spectra[n]), label="Total Spectra")
    # plt.plot(frequency, bb_xcal[n] - bb_ical[n], label="XCAL - ICAL")
    # plt.plot(frequency, bb_xcal[n], label=f"XCAL = {temps['xcal'][n]} K")
    # plt.plot(frequency, bb_ical[n], label=f"ICAL = {temps['ical'][n]} K")
    # plt.legend()

    # plt.show()

    fnyq = gen_nyquistl(
        "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
    )
    frec = 4 * (channels[channel] % 2) + modes[mode]

    apod = fits_data[1].data["APODIZAT"][0]
    R0 = fits_data[1].data["BOLPARM_"][0]
    T0 = fits_data[1].data["BOLPARM2"][0]
    G1 = fits_data[1].data["BOLPARM3"][0]
    beta = fits_data[1].data["BOLPARM4"][0]
    rho = fits_data[1].data["BOLPARM5"][0]
    C1 = fits_data[1].data["BOLPARM6"][0]
    C3 = fits_data[1].data["BOLPARM7"][0]
    Jo = fits_data[1].data["BOLPARM8"][0]
    Jg = fits_data[1].data["BOLPARM9"][0]

    ifg = mu.spec_to_ifg(spec=total_spectra, mtm_speed=mtm_speed, channel=channels[channel], adds_per_group=adds_per_group, bol_cmd_bias=bol_cmd_bias/25.5, # convert to volts
                         bol_volt=bol_volt, Tbol=temps[f"bolometer_{channel}"], gain=gain, sweeps=sweeps, apod=apod,otf=emiss_xcal,fnyq_icm=fnyq["icm"][frec], R0=R0, T0=T0, G1=G1, beta=beta, rho=rho, C1=C1, C3=C3, Jo=Jo, Jg=Jg)
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

                plt.plot(np.abs(total_spectra[0]), label=f"{channel.upper()}{mode.upper()} Total Spectra")
                plt.title(f"{channel.upper()}{mode.upper()} Total Spectra")
                plt.xlabel("Frequency (GHz)")
                plt.ylabel("Brightness Temperature (K)")
                plt.legend()
                plt.grid()
                plt.show()