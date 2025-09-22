import numpy as np
from astropy.io import fits

import globals as g
from utils import my_utils as utils


def get_bolometer_parameters(channel, mode):
    """
    Returns the bolometer parameters for the given channel and mode.
    """
    fits_data = fits.open(
        f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
    )

    R0 = fits_data[1].data["BOLPARM_"][0]
    T0 = fits_data[1].data["BOLPARM2"][0]
    G1 = fits_data[1].data["BOLPARM3"][0]
    beta = fits_data[1].data["BOLPARM4"][0]
    rho = fits_data[1].data["BOLPARM5"][0]
    C1 = fits_data[1].data["BOLPARM6"][0]
    C3 = fits_data[1].data["BOLPARM7"][0]
    Jo = fits_data[1].data["BOLPARM8"][0]
    Jg = fits_data[1].data["BOLPARM9"][0]

    return R0, T0, G1, beta, rho, C1, C3, Jo, Jg


def calculate_dc_response(bol_cmd_bias, bol_volt, Tbol, R0, T0, G1, beta, rho, Jo, Jg):
    """
    Taken largely from calc_responsivity in fsl.
    """
    rscale = 1.0e-7

    cmd_bias = bol_cmd_bias.astype("double")

    V = (bol_volt - Jo) / Jg

    RL = 4.0e7
    R = RL * V / (cmd_bias - V)

    X = V * rho
    Y = R / R0 / X

    mult = Y * Tbol * np.sinh(X / Tbol)
    SQ = np.log(mult)
    Tbol = T0 / SQ**2

    SQ = np.log(Y * Tbol * np.sinh(X / Tbol))

    Tbol = T0 / (SQ * SQ)

    G = G1 * Tbol**beta

    H = Tbol / X * np.tanh(X / Tbol)

    DT = 1.0 / H - 1.0 - 0.5 * np.sqrt(T0 / Tbol)

    Z = (G * Tbol * R + DT * V**2) / (G * Tbol * R / H - DT * V**2)

    S0 = rscale * R * (Z - H) / (V * (Z * R / RL + 1.0) * (H + 1.0))

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


def get_bolometer_response_function(channel, mode, bol_cmd_bias, bol_volt, Tbol):
    R0, T0, G1, beta, rho, C1, C3, Jo, Jg = get_bolometer_parameters(channel, mode)

    # bolometer response function
    S0 = calculate_dc_response(
        bol_cmd_bias=bol_cmd_bias,
        bol_volt=bol_volt,
        Tbol=Tbol,
        R0=R0,
        T0=T0,
        G1=G1,
        beta=beta,
        rho=rho,
        Jo=Jo,
        Jg=Jg,
    )

    tau = calculate_time_constant(
        C3=C3,
        Tbol=Tbol,
        C1=C1,
        G1=G1,
        beta=beta,
        bol_volt=bol_volt,
        Jo=Jo,
        Jg=Jg,
        bol_cmd_bias=bol_cmd_bias,
        rho=rho,
        T0=T0,
    )

    omega = utils.get_afreq(0 if mode[1] == "s" else 1, channel, 257)

    if S0.ndim == 1:
        S0 = S0[:, np.newaxis]
        tau = tau[:, np.newaxis]
    B = S0 / (1 + 1j * omega[np.newaxis, :] * tau)
    return B
