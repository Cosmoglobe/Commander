"""
Order:
0 - XCAL = OTF
1 - ICAL
2 - Dihedral
3 - Refhorn
4 - Skyhorn
5 - Collimator
6 - Bolometer_RH
7 - Bolometer_RL
8 - Bolometer_LH
9 - Bolometer_LL
"""

import random

import h5py
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import calibration.fit_otf as fit_otf
import globals as g
import utils.my_utils as utils
from utils.config import gen_nyquistl

fnyq = gen_nyquistl(
    "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
)

for channel in g.CHANNELS.keys():
    data = np.load(f"{g.PREPROCESSED_DATA_PATH}cal_{channel}.npz", "r")

    mtm_length = data["mtm_length"][:]
    mtm_speed = data["mtm_speed"][:]

    for mode in g.MODES.keys():
        if mode == "lf" and (channel[1] == "h"):
            continue
        # open fitted_emissivities.npy and plot the results
        fitted_emissivities = np.load(
            f"./calibration/output/fitted_emissivities_{channel}_{mode}.npy"
        )
        frequencies = utils.generate_frequencies(channel, mode, 257)

        plt.figure(figsize=(10, 6))
        plt.plot(frequencies, fitted_emissivities.real)
        plt.plot(frequencies, fitted_emissivities.imag, linestyle="dashed")
        plt.xlabel("Frequency (GHz)")
        plt.ylabel("Emissivity")
        plt.title("Fitted Emissivities vs Frequency")
        plt.grid()
        plt.savefig("./calibration/output/fitted_emissivities.png")
        plt.close()

        # compare to the published ones
        frequencies_pub = utils.generate_frequencies(channel, mode)
        fits_data = fits.open(
            f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
        )
        otf = fits_data[1].data["RTRANSFE"][0] + 1j * fits_data[1].data["ITRANSFE"][0]
        otf = otf[np.abs(otf) > 0]
        ical = fits_data[1].data["RICAL"][0] + 1j * fits_data[1].data["IICAL"][0]
        ical = ical[np.abs(ical) > 0]
        dihedral = (
            fits_data[1].data["RDIHEDRA"][0] + 1j * fits_data[1].data["IDIHEDRA"][0]
        )
        dihedral = dihedral[np.abs(dihedral) > 0]
        refhorn = (
            fits_data[1].data["RREFHORN"][0] + 1j * fits_data[1].data["IREFHORN"][0]
        )
        refhorn = refhorn[np.abs(refhorn) > 0]
        skyhorn = (
            fits_data[1].data["RSKYHORN"][0] + 1j * fits_data[1].data["ISKYHORN"][0]
        )
        skyhorn = skyhorn[np.abs(skyhorn) > 0]
        collimator = (
            fits_data[1].data["RSTRUCTU"][0] + 1j * fits_data[1].data["ISTRUCTU"][0]
        )
        collimator = collimator[np.abs(collimator) > 0]
        bolometer = (
            fits_data[1].data["RBOLOMET"][0] + 1j * fits_data[1].data["IBOLOMET"][0]
        )
        bolometer = bolometer[np.abs(bolometer) > 0]

        emissivities = np.zeros((len(frequencies_pub), 7), dtype=complex)
        emissivities[:, 0] = otf
        emissivities[:, 1] = ical
        emissivities[:, 2] = dihedral
        emissivities[:, 3] = refhorn
        emissivities[:, 4] = skyhorn
        emissivities[:, 5] = collimator
        emissivities[:, 6] = bolometer

        fig, ax = plt.subplots(5, 2, figsize=(15, 20))
        ax = ax.flatten()
        labels = [
            "OTF",
            "ICAL",
            "Dihedral",
            "Refhorn",
            "Skyhorn",
            "Collimator",
            "Bolometer",
        ]
        for i in range(fitted_emissivities.shape[1]):
            ax[i].plot(
                frequencies, fitted_emissivities[:, i].real, label="Fitted", color="red"
            )
            ax[i].plot(
                frequencies,
                fitted_emissivities[:, i].imag,
                color="red",
                linestyle="dashed",
            )
            ax[i].plot(
                frequencies_pub,
                emissivities[:, i].real,
                label="Published",
                color="black",
            )
            ax[i].plot(
                frequencies_pub,
                emissivities[:, i].imag,
                color="black",
                linestyle="dashed",
            )
            ax[i].set_xlabel("Frequency (GHz)")
            ax[i].set_ylabel("Emissivity")
            ax[i].set_title(f"{labels[i]} Emissivity Comparison")
            ax[i].legend()
            ax[i].grid()
            ax[i].set_ylim(-1, 1)
        plt.tight_layout()
        plt.savefig("./calibration/output/fitted_emissivities_comparison.png")
        plt.close()

        if mode == "lf" and (channel[1] == "h"):
            continue
        if mode[0] == "s":
            length_filter = mtm_length == 0
        else:
            length_filter = mtm_length == 1
        if mode[1] == "s":
            speed_filter = mtm_speed == 0
        else:
            speed_filter = mtm_speed == 1

        mode_filter = length_filter & speed_filter

        # calculate calibration residuals with both sets of emissivities and compare
        # print(f"data keys: {list(data.keys())}")
        ifgs = ifg = data[f"ifg"][mode_filter]
        gain = data[f"gain"][mode_filter]
        sweeps = data["sweeps"][mode_filter]
        bol_cmd_bias = data[f"bol_cmd_bias"][mode_filter]
        bol_volt = data[f"bol_volt"][mode_filter]

        xcal = data["xcal_cone"][mode_filter]  # TODO: update this
        ical = data["ical"][mode_filter]
        dihedral = data["dihedral"][mode_filter]
        refhorn = data["refhorn"][mode_filter]
        skyhorn = data["skyhorn"][mode_filter]
        collimator = data["collimator"][mode_filter]
        # TODO: fit for all bolometers
        bolometer = data["bolometer"][mode_filter]
        temps = np.array(
            [
                xcal,
                ical,
                dihedral,
                refhorn,
                skyhorn,
                collimator,
                bolometer,
            ]
        )
        # print(f"temps shape: {np.array(temps).shape}")
        adds_per_group = data["adds_per_group"][mode_filter]

        frec = 4 * (g.CHANNELS[channel] % 2) + g.MODES[mode]
        fnyq_icm = fnyq["icm"][frec]

        n = random.randint(0, ifgs.shape[0] - 1)

        D = np.zeros(257, dtype=complex)
        R = np.zeros(257, dtype=complex)
        S = np.zeros(257, dtype=complex)
        for freqi, freq in enumerate(frequencies):
            D[freqi] = fit_otf.D(
                fitted_emissivities[freqi, :],
                freqi,
                ifgs[n],
                channel,
                mode,
                gain[n],
                sweeps[n],
                bol_cmd_bias[n],
                bol_volt[n],
                temps[:, n],
                adds_per_group[n],
                fnyq_icm,
            )
            R[freqi] = fit_otf.R(
                fitted_emissivities[freqi, :], freqi, temps[:, n], frequencies
            )
            S[freqi] = fit_otf.S(freqi, frequencies, temps[:, n])
            # print(f"D, R, S: {D}, {R}, {S}")
        residuals_fitted = D - R - S

        residuals_pub = np.zeros(frequencies_pub.shape[0], dtype=complex)
        for freqi, freq in enumerate(frequencies_pub):
            D_pub = fit_otf.D(
                emissivities[freqi, :],
                freqi,
                ifgs[n],
                channel,
                mode,
                gain[n],
                sweeps[n],
                bol_cmd_bias[n],
                bol_volt[n],
                temps[:, n],
                adds_per_group[n],
                fnyq_icm,
            )
            R_pub = fit_otf.R(emissivities[freqi, :], freqi, temps[:, n], frequencies)
            S_pub = fit_otf.S(freqi, frequencies, temps[:, n])
            residuals_pub[freqi] = D_pub - R_pub - S_pub

        plt.figure(figsize=(10, 6))
        # plt.plot(ifgs[n, :], label="Original IFG", color="black")
        plt.plot(
            frequencies,
            residuals_fitted,
            label="Fitted Emissivities Residuals",
            color="red",
        )
        plt.plot(
            frequencies_pub,
            residuals_pub,
            label="Published Emissivities Residuals",
            color="blue",
        )
        plt.plot(
            frequencies, utils.planck(frequencies, g.T_CMB), label="CMB", color="green"
        )
        plt.xlabel("Frequency (GHz)")
        plt.ylabel("Signal")
        plt.title(f"IFG {n+1} Comparison")
        plt.legend()
        plt.grid()
        plt.savefig("./calibration/output/ifg_residuals_comparison.png")
        plt.close()

        # plot the different components of the residual
        plt.plot(frequencies, residuals_fitted.real, label="Residual", color="red")
        plt.plot(frequencies, residuals_fitted.imag, color="red", linestyle="dashed")
        plt.plot(frequencies, D.real, label="D", color="blue")
        plt.plot(frequencies, D.imag, color="blue", linestyle="dashed")
        plt.plot(frequencies, R.real, label="R", color="green")
        plt.plot(frequencies, R.imag, color="green", linestyle="dashed")
        plt.plot(frequencies, S.real, label="S", color="orange")
        plt.plot(frequencies, S.imag, color="orange", linestyle="dashed")
        plt.plot(frequencies, D.real - R.real, label="D-R", color="purple")
        plt.plot(frequencies, (D - R).imag, color="purple", linestyle="dashed")
        plt.xlabel("Frequency (GHz)")
        plt.ylabel("Signal")
        plt.title(f"IFG {n+1} Residual Components")
        plt.legend()
        plt.grid()
        plt.savefig("./calibration/output/residual_components.png")
        plt.close()
