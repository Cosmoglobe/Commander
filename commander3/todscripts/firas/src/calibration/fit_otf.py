"""
Script to fit an optical transfer function to the calibration data. For now, it follows the model:
S_sky/XCAL = 1/OTF * 1/ETF * 1/Bol * Y - R,
where Y is the Fourier transformed interferogram and R is defined by
R = 1/OTF * sum over all emitters i (excluding XCAL) of E_i * P(T_i).
"""

import time
from functools import partial
from math import isnan
from multiprocessing import Pool, cpu_count

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
        Tbol=temps[6],  # TODO: generalize for other channels
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
        sum += Ei[i] * utils.planck(frequencies[nui], temps[i])

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


def fit_single_frequency(
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
    """
    Fit emissivities for a single frequency index.
    This function is designed to be called in parallel.
    """
    Ei = np.ones(
        temps.shape[0], dtype=complex
    )  # First guess for the emissivities of all emitters

    # check for nans
    if np.isnan(Ei).any():
        raise ValueError("Initial emissivities contain NaNs")
    if np.isnan(ifg).any():
        raise ValueError("Input interferogram contains NaNs")
    if np.isnan(nui):
        raise ValueError("Frequency index is NaN")
    if np.isnan(gain).any():
        raise ValueError("Gain array contains NaNs")
    if np.isnan(sweeps).any():
        raise ValueError("Sweeps array contains NaNs")
    if np.isnan(bol_cmd_bias).any():
        raise ValueError("Bolometer command bias array contains NaNs")
    if np.isnan(bol_volt).any():
        raise ValueError("Bolometer voltage array contains NaNs")
    if np.isnan(temps).any():
        raise ValueError("Temperatures array contains NaNs")
    if np.isnan(adds_per_group).any():
        raise ValueError("Adds per group array contains NaNs")
    if np.isnan(fnyq_icm).any():
        raise ValueError("Nyquist frequency array contains NaNs")
    if np.isnan(frequencies).any():
        raise ValueError("Frequencies array contains NaNs")

    result = minimize(
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
    )
    # print output message
    print(f"Message for frequency {nui}")
    print(f"Success: {result.success}")
    print(f"Result message: {result.message}")

    return nui, result.x


if __name__ == "__main__":
    # TODO: generalize
    channel = "ll"
    mode = "ss"

    data = np.load(f"{g.PREPROCESSED_DATA_PATH}cal_{channel}.npz", "r")

    ifg = data[f"ifg"][:]

    print(f"Data loaded from {g.PREPROCESSED_DATA_PATH}cal_{channel}.npz")
    print("Preparing data for fitting...")
    # print(f"data keys: {list(data.keys())}")

    # TODO: fit for temp weight coefficients too
    # xcal = (data["a_xcal"][:] + data["b_xcal"][:]) / 2
    # ical = (data["a_ical"][:] + data["b_ical"][:]) / 2
    # dihedral = (data["a_dihedral"][:] + data["b_dihedral"][:]) / 2
    # refhorn = (data["a_refhorn"][:] + data["b_refhorn"][:]) / 2
    # skyhorn = (data["a_skyhorn"][:] + data["b_skyhorn"][:]) / 2
    # collimator = (data["a_collimator"][:] + data["b_collimator"][:]) / 2
    # bolometer_ll = (data["a_bol_assem_ll"][:] + data["b_bol_assem_ll"][:]) / 2
    # bolometer_lh = (data["a_bol_assem_lh"][:] + data["b_bol_assem_lh"][:]) / 2
    # bolometer_rl = (data["a_bol_assem_rl"][:] + data["b_bol_assem_rl"][:]) / 2
    # bolometer_rh = (data["a_bol_assem_rh"][:] + data["b_bol_assem_rh"][:]) / 2
    xcal = data["xcal_cone"][:]  # TODO: update this
    ical = data["ical"][:]
    dihedral = data["dihedral"][:]
    refhorn = data["refhorn"][:]
    skyhorn = data["skyhorn"][:]
    collimator = data["collimator"][:]
    # TODO: fit for all bolometers
    bol = data["bolometer"][:]

    adds_per_group = data["adds_per_group"][:]
    sweeps = data["sweeps"][:]

    bol_cmd_bias = data[f"bol_cmd_bias"][:]
    bol_volt = data[f"bol_volt"][:]
    gain = data[f"gain"][:]

    # check where nans in temperature come from
    nan_temps_indices = np.unique(
        np.where(
            np.isnan(xcal)
            | np.isnan(ical)
            | np.isnan(dihedral)
            | np.isnan(refhorn)
            | np.isnan(skyhorn)
            | np.isnan(collimator)
            | np.isnan(bol)
        )[0]
    )
    if len(nan_temps_indices) > 0:
        print(f"Found NaNs in temperature arrays at indices: {nan_temps_indices}")

    # does it correspond to bol_cmd_bias <0?
    nan_temps_bol_cmd_bias = bol_cmd_bias[nan_temps_indices]
    if np.all(nan_temps_bol_cmd_bias < 0):
        print(
            "All NaNs in temperature arrays correspond to bol_cmd_bias < 0, filtering them out."
        )

    # is it nans in all of the temps or just one object?
    for idx in nan_temps_indices:
        print(
            f"Index {idx}: xcal={xcal[idx]}, ical={ical[idx]}, dihedral={dihedral[idx]}, refhorn={refhorn[idx]}, skyhorn={skyhorn[idx]}, collimator={collimator[idx]}, bolometer={bol[idx]}"
        )

    temps = np.array([xcal, ical, dihedral, refhorn, skyhorn, collimator, bol])

    frequencies = utils.generate_frequencies(channel, mode, 257)

    fnyq = gen_nyquistl(
        "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
    )
    frec = 4 * (g.CHANNELS[channel] % 2) + g.MODES[mode]
    fnyq_icm = fnyq["icm"][frec]

    solution = np.zeros(
        (g.SPEC_SIZE, temps.shape[0]), dtype=complex
    )  # 10 emissivities to fit

    # Determine number of processes to use
    n_processes = cpu_count()
    print(f"Using {n_processes} CPU cores for parallel processing")

    # Create partial function with fixed arguments
    fit_func = partial(
        fit_single_frequency,
        ifg=ifg,
        channel=channel,
        mode=mode,
        gain=gain,
        sweeps=sweeps,
        bol_cmd_bias=bol_cmd_bias,
        bol_volt=bol_volt,
        temps=temps,
        adds_per_group=adds_per_group,
        fnyq_icm=fnyq_icm,
        frequencies=frequencies,
    )

    start0 = time.time()

    # Parallel processing
    with Pool(processes=n_processes) as pool:
        results = pool.map(fit_func, range(g.SPEC_SIZE))

    # Collect results
    for nui, emissivities in results:
        solution[nui] = emissivities
        print(f"Frequency {nui+1}/{g.SPEC_SIZE} - Solution: {emissivities}")

    print(f"Total time taken: {time.time() - start0:.2f} seconds")
    # save solution
    np.save("fitted_emissivities.npy", solution)
