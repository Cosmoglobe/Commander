"""
Script to fit an optical transfer function to the calibration data. For now, it follows the model:
S_sky/XCAL = 1/OTF * 1/ETF * 1/Bol * Y - R,
where Y is the Fourier transformed interferogram and R is defined by
R = 1/OTF * sum over all emitters i (excluding XCAL) of E_i * P(T_i).
"""

import time

import h5py
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize

import globals as g
from calibration import bolometer, electronics
from utils import my_utils as utils
from utils.config import gen_nyquistl


def D(
    Ei,
    nui,
    ifg,
    channel,
    mode,
    gain,
    sweeps,
    bol_cmd_bias,
    bol_volt,
    temps,
    adds_per_group,
    fnyq_icm,
    apod=False,
):
    if apod:
        fits_data = fits.open(
            f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
        )
        apod_func = fits_data[1].data["APODIZAT"][0]

    # subtract dither
    if ifg.ndim == 1:
        ifg = ifg - np.median(ifg)

        ifg = ifg / gain / sweeps

        if apod:
            ifg = ifg * apod_func

        ifg = np.roll(ifg, -g.PEAK_POSITIONS[f"{channel}_{mode}"])
        Y = np.fft.rfft(ifg)
        Y = Y[nui]
    else:
        ifg = ifg - np.median(ifg, axis=1)[:, None]

        # "normalize" ifg by gain and sweeps
        ifg = ifg / gain[:, np.newaxis] / sweeps[:, np.newaxis]

        if apod:
            ifg = ifg * apod_func

        ifg = np.roll(ifg, -g.PEAK_POSITIONS[f"{channel}_{mode}"], axis=1)
        Y = np.fft.rfft(ifg, axis=1)
        Y = Y[:, nui]

    B = bolometer.get_bolometer_response_function(
        channel,
        mode,
        bol_cmd_bias / 25.5,
        bol_volt,
        Tbol=temps[9],  # TODO: generalize for other channels
    )
    B = B[:, nui]

    # electronics transfer function
    samprate = 681.43  # from fex_samprate
    Z = electronics.etfunction(channel, adds_per_group, samprate, nui)
    # plt.plot(Z.real)
    # plt.plot(Z.imag)
    # plt.show()

    H = Ei[0]

    spec_norm = fnyq_icm * g.FAC_ETENDU * g.FAC_ADC_SCALE

    # print(f"Y: {Y}")
    # print(f"H: {H}")
    # print(f"Z: {Z[Z == 0]}")
    # print(f"Z: {Z.shape}")
    # print(f"B: {B}")
    return Y / H / Z / B / spec_norm * g.FAC_ERG_TO_MJY


def R(Ei, nui, temps, frequencies):
    """
    R is the function that weights each of the emitters.

    Parameters
    ----------
    H : array_like
        Optical transfer function a.k.a. emissivity of the XCAL.
    Ei : array_like
        2D array of the emissivities over frequency for each of the other nine emitters.
        1: ICAL, 2: dihedral, 3: refhorn, 4: skyhorn, 5: collimator, 6: bolometer_rh, 7: bolometer_rl, 8: bolometer_lh, 9: bolometer_ll
    """
    sum = np.zeros_like(temps[0], dtype=complex)
    for i in range(1, Ei.shape[0]):
        sum += Ei[i] * utils.planck(temps[i], frequencies[nui])

    H = Ei[0]

    return sum / H


def S(nui, frequencies, temps):
    return utils.planck(frequencies[nui], temps[0])


def full_function(
    Ei,
    nui,
    ifg,
    channel,
    mode,
    gain,
    sweeps,
    bol_cmd_bias,
    bol_volt,
    temps,
    adds_per_group,
    fnyq_icm,
    frequencies,
):
    result = (
        D(
            Ei,
            nui,
            ifg,
            channel,
            mode,
            gain,
            sweeps,
            bol_cmd_bias,
            bol_volt,
            temps,
            adds_per_group,
            fnyq_icm,
        )
        - R(Ei, nui, temps, frequencies)
        - S(nui, frequencies, temps)
    )
    return np.sum(result**2)


# def real_to_complex(z):
#     """
#     Real vector of length 2n -> complex of length n.
#     Taken from https://stackoverflow.com/questions/51211055/can-scipy-optimize-minimize-functions-of-complex-variables-at-all-and-how.
#     """
#     return z[: len(z) // 2] + 1j * z[len(z) // 2 :]


# def complex_to_real(z):
#     """
#     Complex vector of length n -> real of length 2n
#     Taken from https://stackoverflow.com/questions/51211055/can-scipy-optimize-minimize-functions-of-complex-variables-at-all-and-how.
#     """
#     return np.concatenate((np.real(z), np.imag(z)))


if __name__ == "__main__":
    data = h5py.File(g.PREPROCESSED_DATA_PATH_CAL, "r")["df_data"]

    # TODO: generalize
    channel = "ll"
    mode = "ss"

    ifg = data[f"ifg_{channel}"][:]

    # TODO: fit for temp weight coefficients too
    xcal = (data["a_xcal"][:] + data["b_xcal"][:]) / 2
    ical = (data["a_ical"][:] + data["b_ical"][:]) / 2
    dihedral = (data["a_dihedral"][:] + data["b_dihedral"][:]) / 2
    refhorn = (data["a_refhorn"][:] + data["b_refhorn"][:]) / 2
    skyhorn = (data["a_skyhorn"][:] + data["b_skyhorn"][:]) / 2
    collimator = (data["a_collimator"][:] + data["b_collimator"][:]) / 2
    bolometer_ll = (data["a_bol_assem_ll"][:] + data["b_bol_assem_ll"][:]) / 2
    bolometer_lh = (data["a_bol_assem_lh"][:] + data["b_bol_assem_lh"][:]) / 2
    bolometer_rl = (data["a_bol_assem_rl"][:] + data["b_bol_assem_rl"][:]) / 2
    bolometer_rh = (data["a_bol_assem_rh"][:] + data["b_bol_assem_rh"][:]) / 2

    adds_per_group = data["adds_per_group"][:]
    sweeps = data["sweeps"][:]

    bol_cmd_bias = data[f"bol_cmd_bias_{channel}"][:]
    bol_volt = data[f"bol_volt_{channel}"][:]
    gain = data[f"gain_{channel}"][:]

    # we don't want the records where the bol_cmd_bias is lower than 0
    valid = bol_cmd_bias > 0
    ifg = ifg[valid]
    xcal = xcal[valid]
    ical = ical[valid]
    dihedral = dihedral[valid]
    refhorn = refhorn[valid]
    skyhorn = skyhorn[valid]
    collimator = collimator[valid]
    bolometer_rh = bolometer_rh[valid]
    bolometer_rl = bolometer_rl[valid]
    bolometer_lh = bolometer_lh[valid]
    bolometer_ll = bolometer_ll[valid]
    adds_per_group = adds_per_group[valid]
    sweeps = sweeps[valid]
    bol_cmd_bias = bol_cmd_bias[valid]
    bol_volt = bol_volt[valid]
    gain = gain[valid]

    temps = [
        xcal,
        ical,
        dihedral,
        refhorn,
        skyhorn,
        collimator,
        bolometer_rh,
        bolometer_rl,
        bolometer_lh,
        bolometer_ll,
    ]

    frequencies = utils.generate_frequencies(channel, mode, 257)

    fnyq = gen_nyquistl(
        "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
    )
    frec = 4 * (g.CHANNELS[channel] % 2) + g.MODES[mode]
    fnyq_icm = fnyq["icm"][frec]

    solution = np.zeros((g.SPEC_SIZE, 10), dtype=complex)  # 10 emissivities to fit
    start0 = time.time()
    for nui in range(0, g.SPEC_SIZE):
        print(f"Fitting frequency {nui+1}/{g.SPEC_SIZE}...")

        start = time.time()
        Ei = np.ones(
            10, dtype=complex
        )  # First guess for the emissivities of all emitters
        solution[nui] = minimize(
            full_function,
            Ei,
            args=(
                nui,
                ifg,
                channel,
                mode,
                gain,
                sweeps,
                bol_cmd_bias,
                bol_volt,
                temps,
                adds_per_group,
                fnyq_icm,
                frequencies,
            ),
        ).x
        # solution[nui] = minimize(
        #     lambda z, nui: full_function(real_to_complex(z), nui),
        #     x0=complex_to_real(Ei),
        #     args=(nui,),
        # ).x
        # solution[nui] = newton(full_function, Ei, args=(nui, ifg))
        # solution[nui] = dual_annealing(full_function, bounds=[(0, 2)] * 10, args=(nui,))
        print(f"Solution is {solution[nui]}")
        print(f"Time taken: {time.time() - start:.2f} seconds")

    print(f"Total time taken: {time.time() - start0:.2f} seconds")
    # save solution
    np.save("fitted_emissivities.npy", solution)
