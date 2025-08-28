"""
Script to fit an optical transfer function to the calibration data. For now, it follows the model:
S_sky/XCAL = 1/OTF * 1/ETF * 1/Bol * Y - R,
where Y is the Fourier transformed interferogram and R is defined by
R = 1/OTF * sum over all emitters i (excluding XCAL) of E_i * P(T_i).
"""

import h5py
import numpy as np

import globals as g
from calibration import bolometer, electronics
from utils import my_utils as utils
from scipy.optimize import minimize
import time


def D(Ei, nui):
    ifg = data[f"ifg_{channel}"][:]
    Y = np.fft.rfft(ifg, axis=1)
    Y = Y[:, nui]

    bol_cmd_bias = data[f"bol_cmd_bias_{channel}"][:]
    bol_volt = data[f"bol_volt_{channel}"][:]

    R0, T0, G1, beta, rho, C1, C3, Jo, Jg = bolometer.get_bolometer_parameters(
        channel, mode
    )

    # bolometer response function
    S0 = bolometer.calculate_dc_response(
        bol_cmd_bias=bol_cmd_bias / 25.5,
        bol_volt=bol_volt,
        Jo=Jo,
        Jg=Jg,
        Tbol=temps[9],  # TODO: generalize for other channels
        rho=rho,
        R0=R0,
        T0=T0,
        beta=beta,
        G1=G1,
    )

    tau = bolometer.calculate_time_constant(
        C3=C3,
        Tbol=temps[9],  # TODO: generalize for other channels
        C1=C1,
        G1=G1,
        beta=beta,
        bol_volt=bol_volt,
        Jo=Jo,
        Jg=Jg,
        bol_cmd_bias=bol_cmd_bias / 25.5,
        rho=rho,
        T0=T0,
    )

    omega = utils.get_afreq(0 if mode[1] == "s" else 1, channel, 257)

    B = S0[:, np.newaxis] / (1 + 1j * omega[np.newaxis, :] * tau[:, np.newaxis])
    B = B[:, nui]

    # electronics transfer function
    samprate = 681.43  # from fex_samprate
    Z = electronics.etfunction(channel, adds_per_group, samprate)[:, nui]

    H = Ei[0]

    return Y / H / Z / B


def R(Ei, nui):
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
    sum = np.zeros_like(adds_per_group, dtype=complex)
    for i in range(1, Ei.shape[0]):
        sum += Ei[i] * utils.planck(temps[i], frequencies[nui])

    H = Ei[0]

    return sum / H


def S(nui):
    return utils.planck(frequencies[nui], temps[0])


def full_function(Ei, nui):
    result = D(Ei, nui) - R(Ei, nui) - S(nui)
    return np.sum(result**2)


if __name__ == "__main__":
    data = h5py.File(g.PREPROCESSED_DATA_PATH_CAL, "r")["df_data"]

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

    adds_per_group = data["adds_per_group"][:]
    sweeps = data["sweeps"][:]

    # TODO: generalize
    channel = "ll"
    mode = "ss"

    frequencies = utils.generate_frequencies(channel, mode, 257)

    solution = np.zeros((g.SPEC_SIZE, 10), dtype=complex)  # 10 emissivities to fit
    start0 = time.time()
    for nui in range(g.SPEC_SIZE):
        print(f"Fitting frequency {nui+1}/{g.SPEC_SIZE}...")

        start = time.time()
        Ei = np.zeros(
            10, dtype=complex
        )  # First guess for the emissivities of all emitters
        solution[nui] = minimize(full_function, Ei, args=(nui,)).x
        print(f"Time taken: {time.time() - start:.2f} seconds")

    print(f"Total time taken: {time.time() - start0:.2f} seconds")
    # save solution
    np.save("fitted_emissivities.npy", solution)
