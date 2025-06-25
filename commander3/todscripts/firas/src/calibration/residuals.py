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
from simulations.main import generate_ifg
from utils.config import gen_nyquistl

# get temperatures
data = np.load(g.PROCESSED_DATA_PATH_CAL)

print("Loading calibration data...")

# channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}
channels = {"ll":3}
# modes = {"ss": 0, "lf": 3} 
modes = {"ss":0}

for mode in modes:
    xcal = data[f"xcal_{mode}"][:]
    ical = data[f"ical_{mode}"][:]
    dihedral = data[f"dihedral_{mode}"][:]
    refhorn = data[f"refhorn_{mode}"][:]
    skyhorn = data[f"skyhorn_{mode}"][:]
    bolometer_ll = data[f"bolometer_ll_{mode}"][:]
    bolometer_lh = data[f"bolometer_lh_{mode}"][:]
    bolometer_rl = data[f"bolometer_rl_{mode}"][:]
    bolometer_rh = data[f"bolometer_rh_{mode}"][:]

    temps = {
        "xcal": xcal,
        "ical": ical,
        "dihedral": dihedral,
        "refhorn": refhorn,
        "skyhorn": skyhorn,
        "bolometer_ll": bolometer_ll,
        "bolometer_lh": bolometer_lh,
        "bolometer_rl": bolometer_rl,
        "bolometer_rh": bolometer_rh,
    }

    adds_per_group = data[f"adds_per_group_{mode}"][:]
    sweeps = data[f"sweeps_{mode}"][:]
    
    for channel in channels:
        if not(mode == "lf" and channel[1] == "h"):
            bol_cmd_bias = data[f"bol_cmd_bias_{channel}_{mode}"][:]
            bol_volt = data[f"bol_volt_{channel}_{mode}"][:]
            gain = data[f"gain_{channel}_{mode}"][:]

            simulated_ifgs, simulated_spectra = generate_ifg(
                channel=channel,
                mode=mode,
                temps=temps,
                adds_per_group=adds_per_group,
                sweeps=sweeps,
                bol_cmd_bias=bol_cmd_bias,
                bol_volt=bol_volt,
                gain=gain
            )

            print(f"Simulated IFGs for {channel.upper()} {mode.upper()}")
            print(f"shape of simulated IFGs: {simulated_ifgs.shape}")

            original_ifgs = data[f"ifg_{channel}_{mode}"][:] - np.median(data[f"ifg_{channel}_{mode}"][:], axis=1, keepdims=True)

            # stuff to process the original IFGs into spectra
            fnyq = gen_nyquistl("../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int")
            frec = 4 * (channels[channel] % 2) + modes[mode]
            fits_data = fits.open(
                f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
            )
            otf = (
                fits_data[1].data["RTRANSFE"][0]
                + 1j * fits_data[1].data["ITRANSFE"][0]
            )
            otf = otf[np.abs(otf) > 0]
            # bolometer parameters
            R0 = fits_data[1].data[
                "BOLPARM_"
            ][0]
            T0 = fits_data[1].data[
                "BOLPARM2"
            ][0]
            G1 = fits_data[1].data[
                "BOLPARM3"
            ][0]
            beta = fits_data[1].data[
                "BOLPARM4"
            ][0]
            rho = fits_data[1].data[
                "BOLPARM5"
            ][0]
            C1 = fits_data[1].data[
                "BOLPARM6"
            ][0]
            C3 = fits_data[1].data[
                "BOLPARM7"
            ][0]
            Jo = fits_data[1].data[
                "BOLPARM8"
            ][0]
            Jg = fits_data[1].data[
                "BOLPARM9"
            ][0]

            _, processed_spectra = mu.ifg_to_spec(original_ifgs, mtm_speed=0 if mode[1] == "s" else 1, channel=channels[channel], adds_per_group=adds_per_group, sweeps=sweeps, bol_cmd_bias=bol_cmd_bias, bol_volt=bol_volt, gain=gain, fnyq_icm=fnyq["icm"][frec], otf=otf, Jo=Jo, Jg=Jg, T0=T0, R0=R0, G1=G1, C1=C1, C3=C3, beta=beta, rho=rho, Tbol=temps[f"bolometer_{channel}"], apod=fits_data[1].data["APODIZAT"][0])

            n = 4757

            plt.plot(
                np.abs(simulated_spectra[n, :]),
                label="Simulated Spectra",
                color="blue",
            )
            plt.plot(
                np.abs(processed_spectra[n, :]),
                label="Original Spectra",
                color="orange",
            )
            plt.title(f"{channel.upper()} {mode.upper()} Spectra {n+1}")
            plt.xlabel("Frequency (GHz)")
            plt.ylabel("Brightness Temperature (K)")
            plt.legend()
            plt.show()

            print(f"original ifg peak: {np.max(np.abs(original_ifgs[n]))}")
            print(f"simulated ifg peak: {np.nanmax(np.abs(simulated_ifgs[n]))}")
            print(f"ratio (original/simulated): {np.max(np.abs(original_ifgs[n])) / np.nanmax(np.abs(simulated_ifgs[n]))}")

            # for ifg in range(simulated_ifgs.shape[0]):
            for ifg in range(n, n+1, 1):
            # plot both and residuals
                plt.subplots(3, 1, figsize=(10, 15), sharex=True, sharey=True)
                plt.suptitle(f"{channel.upper()} {mode.upper()} IFG {ifg+1}")

                plt.subplot(3, 1, 1)
                plt.plot(original_ifgs[ifg, :])
                plt.title("Original IFG")
            
                plt.subplot(3, 1, 2)
                plt.plot(simulated_ifgs[ifg, :])
                plt.title("Simulated IFG")

                plt.subplot(3, 1, 3)
                plt.plot(original_ifgs[ifg, :] - simulated_ifgs[ifg, :])
                plt.title("Residuals (original - simulated)")

                plt.show()