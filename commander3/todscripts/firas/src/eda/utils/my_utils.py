import matplotlib.pyplot as plt
import numpy as np
from numba import jit, njit, prange
from utils.frd import apodl, elex_transfcnl
from utils.fut import apod_recnuml, get_recnum


# @njit(parallel=True)
def calculate_dc_response(bol_cmd_bias, bol_volt, Jo, Jg, Tbol, rho, R0, T0, beta, G1):
    """
    Taken largely from calc_responsivity in fsl.
    """

    rscale = 1.0e-7

    cmd_bias = bol_cmd_bias.astype(
        "double"
    )  # / 25.5  # only use the factor for when it's not in volts? i.e. don't use for the data from lambda

    # print("cmd_bias:", cmd_bias)

    V = (bol_volt - Jo) / Jg

    # print("V:", V)

    RL = 4.0e7
    R = RL * V / (cmd_bias - V)

    # print("R:", R, "RL:", RL)

    X = V * rho
    Y = R / R0 / X

    # print("X:", X, "Y:", Y)

    SQ = np.log(Y * Tbol * np.sinh(X / Tbol))
    Tbol = T0 / SQ**2

    # print("SQ:", SQ, "Tbol:", Tbol)

    SQ = np.log(Y * Tbol * np.sinh(X / Tbol))

    # print("SQ:", SQ)

    Tbol = T0 / (SQ * SQ)

    # print("Tbol:", Tbol)

    G = G1 * Tbol**beta

    # print("G:", G)

    H = Tbol / X * np.tanh(X / Tbol)

    # print("H:", H)

    DT = 1.0 / H - 1.0 - 0.5 * np.sqrt(T0 / Tbol)

    # print("DT:", DT)

    Z = (G * Tbol * R + DT * V**2) / (G * Tbol * R / H - DT * V**2)

    # print("Z:", Z)

    S0 = rscale * R * (Z - H) / (V * (Z * R / RL + 1.0) * (H + 1.0))

    # print("S0:", S0)

    return S0


def calculate_time_constant(
    C3, Tbol, C1, G1, beta, bol_volt, Jo, Jg, bol_cmd_bias, rho, T0
):
    """
    Taken largely from calc_responsivity in fsl.
    """
    C = C3 * Tbol**3 + C1 * Tbol

    G = G1 * Tbol**beta

    V = (bol_volt - Jo) / Jg
    RL = 4.0e7
    R = RL * V / (bol_cmd_bias - V)

    X = V * rho
    H = Tbol / X * np.tanh(X / Tbol)
    DT = 1.0 / H - 1.0 - 0.5 * np.sqrt(T0 / Tbol)

    Z = (G * Tbol * R + DT * V**2) / (G * Tbol * R / H - DT * V**2)

    tau = C / G * (Z + 1.0) * (R * H + RL) / ((Z * R + RL) * (H + 1.0))

    return tau


# @njit(parallel=True)
def planck(freq, temp):
    """
    Planck function returning in units of MJy/sr.
    Input frequency in GHz and temperature in K.
    """
    h = 6.62607015e-34 * 1e9  # J GHz-1
    c = 299792458e-9  # m GHz
    k = 1.380649e-23  # J K-1

    if temp.shape != ():
        freq = freq[np.newaxis, :]
        temp = temp[:, np.newaxis]

    b = 2 * h * freq**3 / c**2 / (np.exp(h * freq / (k * temp)) - 1) * 1e20  # MJy sr-1

    return b


def mjy_to_krj(frequency, intensity):
    """
    Convert intensity in MJy/sr to uK_RJ.
    """
    k = 1.380649e-16  # erg/K
    c = 2.99792458e10  # cm/s

    intensity_erg = intensity * 1e-23 * 1e-9  # MJy/sr to erg/s/cm2/Hz
    frequency_hz = frequency * 1e9  # GHz to Hz

    temperature = intensity_erg / (2 * k * frequency_hz**2 / c**2)
    temperature_uk = temperature * 1e6
    return temperature_uk
