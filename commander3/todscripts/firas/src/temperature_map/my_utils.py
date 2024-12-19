from datetime import datetime, timedelta

import numpy as np
from numba import njit, prange
from utils.frd import apodl, elex_transfcnl
from utils.fut import apod_recnuml, get_recnum


# @njit(parallel=True)
def calculate_dc_response(bol_cmd_bias, bol_volt, Jo, Jg, Tbol, rho, R0, T0, beta, G1):
    rscale = 1.0e-7

    cmd_bias = np.double(bol_cmd_bias) / 25.5

    V = (bol_volt - Jo) / Jg

    RL = 4.0e7
    R = RL * V / (cmd_bias - V)

    X = V * rho
    Y = R / R0 / X

    SQ = np.log(Y * Tbol * np.sinh(X / Tbol))
    Tbol = T0 / SQ**2

    SQ = np.log(Y * Tbol * np.sinh(X / Tbol))

    Tbol = T0 / (SQ * SQ)

    G = G1 * Tbol**beta

    H = Tbol / X * np.tanh(X / Tbol)

    DT = 1.0 / H - 1.0 - 0.5 * np.sqrt(T0 / Tbol)

    Z = (G * Tbol * R + DT * V**2) / (G * Tbol * R / H - DT * V**2)

    S0 = rscale * R * (Z - H) / (V * (Z * R / RL + 1.0) * (H + 1.0))

    return S0


@njit(parallel=True)
def my_median(arr):
    sort = np.zeros_like(arr)
    for i in prange(len(arr)):
        tmp = np.sort(arr[i])
        sort[i] = tmp

    n = arr.shape[1]
    # n is always even so we can just use this formula

    median = (sort[:, n // 2] + sort[:, n // 2 - 1]) / 2
    return median


# @njit(parallel=True)
def clean_ifg(
    ifg,
    mtm_length,
    mtm_speed,
    channel,
    adds_per_group,
    gain,
    sweeps,
):
    # subtract dither
    ifg = ifg - my_median(ifg)[:, np.newaxis]

    ifg = ifg / gain[:, np.newaxis] / sweeps[:, np.newaxis]

    # apodize
    sm = 2 * mtm_length + mtm_speed

    arecno = apod_recnuml(channel, sm, adds_per_group).astype(np.int32)
    apodl_all = apodl()
    apod = apodl_all[arecno, :]

    ifg = ifg * apod

    # roll
    peak_pos = 360
    ifg = np.roll(ifg, -peak_pos)

    return ifg


# @njit(parallel=True)
def ifg_to_spec(
    ifg,
    mtm_speed,
    channel,
    adds_per_group,
    bol_cmd_bias,
    bol_volt,
    fnyq_icm,
    fnyq_hz,
    otf,
    Jo,
    Jg,
    Tbol,
    rho,
    R0,
    T0,
    beta,
    G1,
    tau,
):
    # fft
    spec = np.fft.rfft(ifg)

    # etf
    etfl_all = elex_transfcnl(samprate=681.43, nfreq=len(spec[0]))

    erecno = get_recnum(mtm_speed, channel, adds_per_group).astype(np.int32)
    etf = etfl_all[erecno, :]

    fac_etendu = 1.5  # nathan's pipeline
    fac_adc_scale = 204.75  # nathan's pipeline
    spec_norm = fnyq_icm * fac_etendu * fac_adc_scale

    spec = spec / (etf * spec_norm)

    # fcc_spec_length = 321
    spec_len = len(ifg[0]) // 2 + 1
    dw = 2.0 * np.pi * fnyq_hz / spec_len
    # afreq = np.arange(fcc_spec_length) * dw # had to change because i'm not padding for now
    afreq = np.arange(spec_len) * dw

    S0 = calculate_dc_response(
        bol_cmd_bias, bol_volt, Jo, Jg, Tbol, rho, R0, T0, beta, G1
    )
    B = 1.0 + 1j * tau * afreq

    spec = B[np.newaxis, :] * spec / S0[:, np.newaxis]

    # optical transfer function
    spec = spec[:, : len(otf)] / otf

    fac_icm_ghz = 29.9792458
    fac_erg_to_mjy = 1.0e8 / fac_icm_ghz

    spec = spec * fac_erg_to_mjy

    return spec


def planck(freq, temp):
    """
    Planck function returning in units of MJy/sr.
    Input frequency in GHz and temperature in K.
    """
    h = 6.62607015e-34 * 1e9  # J GHz-1
    c = 299792458e-9  # m GHz
    k = 1.380649e-23  # J K-1

    print(f"freq: {freq.shape}, temp: {temp.shape}")
    if temp.shape != ():
        freq = freq[np.newaxis, :]
        temp = temp[:, np.newaxis]

    # b = np.zeros((len(temp), len(freq)))
    b = 2 * h * freq**3 / c**2 / (np.exp(h * freq / (k * temp)) - 1) * 1e20  # MJy sr-1

    return b


def residuals(temperature, frequency_data, sky_data):  # doing least squares for now
    """
    Function to calculate the residuals between the sky data and the Planck function for fitting a sky/XCAL temperature.
    """
    return np.sum((sky_data - planck(frequency_data, temperature[0])) ** 2)
