import matplotlib.pyplot as plt
import numpy as np
from numba import jit, njit, prange
from utils.frd import apodl, elex_transfcnl
from utils.fut import apod_recnuml, get_recnum


# @njit(parallel=True)
def calculate_dc_response(bol_cmd_bias, bol_volt, Jo, Jg, Tbol, rho, R0, T0, beta, G1):
    rscale = 1.0e-7

    cmd_bias = (
        bol_cmd_bias.astype("double") / 25.5
    )  # only use the factor for when it's not in volts? i.e. don't use for the data from lambda

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


# @njit(parallel=True)
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
    median_ifg = np.expand_dims(my_median(ifg), axis=-1)

    # subtract dither
    ifg = ifg - median_ifg

    # Ensure gain and sweeps are reshaped for broadcasting
    gain = np.expand_dims(gain, axis=-1)
    sweeps = np.expand_dims(sweeps, axis=-1)

    ifg = ifg / gain / sweeps

    print("ifg shape after division:", ifg.shape)

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


# @njit(parallel=True) - can't use it because of rfft
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

    # print("spec after rfft:", spec)

    # etf
    etfl_all = elex_transfcnl(samprate=681.43, nfreq=len(spec[0]))

    erecno = get_recnum(mtm_speed, channel, adds_per_group).astype(np.int32)

    etf = etfl_all[erecno, :]

    # print("etf:", etf)

    fac_etendu = 1.5  # nathan's pipeline
    fac_adc_scale = 204.75  # nathan's pipeline
    spec_norm = fnyq_icm * fac_etendu * fac_adc_scale

    spec = spec / etf
    # print("spec after etf:", spec)
    spec = spec / spec_norm

    # plot etf spec over time
    plt.imshow(np.abs(spec).T, aspect="auto", extent=[0, len(spec), 0, len(spec[0])])
    plt.savefig("../../output/plots/etf_spec_over_time.png")
    plt.clf()

    # fcc_spec_length = 321
    spec_len = len(ifg[0]) // 2 + 1
    dw = 2.0 * np.pi * fnyq_hz / spec_len
    # afreq = np.arange(fcc_spec_length) * dw # had to change because i'm not padding for now
    afreq = np.arange(spec_len) * dw

    # print("afreq:", afreq)

    S0 = calculate_dc_response(
        bol_cmd_bias=bol_cmd_bias,
        bol_volt=bol_volt,
        Jo=Jo,
        Jg=Jg,
        Tbol=Tbol,
        rho=rho,
        R0=R0,
        T0=T0,
        beta=beta,
        G1=G1,
    )

    # print("S0:", S0)

    B = 1.0 + 1j * tau * afreq

    # print("B:", B)

    spec = B[np.newaxis, :] * spec / S0[:, np.newaxis]

    # print("spec after dc response:", spec)

    # plot spec after bolometer transfer function
    plt.imshow(np.abs(spec).T, aspect="auto", extent=[0, len(spec), 0, len(spec[0])])
    plt.savefig("../../output/plots/spec_after_bolometer_transfer_function.png")
    plt.clf()

    # optical transfer function
    spec = spec[:, : len(otf)] / otf

    # plot spec after optical transfer function
    plt.imshow(np.abs(spec).T, aspect="auto", extent=[0, len(spec), 0, len(spec[0])])
    plt.savefig("../../output/plots/spec_after_optical_transfer_function.png")
    plt.clf()

    fac_icm_ghz = 29.9792458
    fac_erg_to_mjy = 1.0e8 / fac_icm_ghz

    spec = spec * fac_erg_to_mjy

    return spec


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


# @njit(parallel=True)
def residuals(temperature, frequency_data, sky_data):  # doing least squares for now
    """
    Function to calculate the residuals between the sky data and the Planck function for fitting a sky/XCAL temperature.
    Works with numpy arrays.
    """
    res = np.sum((sky_data - planck(frequency_data, temperature)) ** 2, axis=1)
    return res


def filter_crap(
    stat_word_5, stat_word_9, stat_word_13, stat_word_16, lvdt_stat_a, lvdt_stat_b
):
    filter1 = (stat_word_5 != 16641) & (
        stat_word_5 != 17921
    )  # & (stat_word_5 != 19457)
    print("filter1:", np.count_nonzero(filter1))
    filter2 = stat_word_9 != 15414
    print("filter2:", np.count_nonzero(filter2))
    filter3 = (
        (stat_word_13 != 17345)
        & (stat_word_13 != 25201)
        & (stat_word_13 != 26945)
        & (stat_word_13 != 27073)
        & (stat_word_13 != 27649)
        & (stat_word_13 != 27697)
        & (stat_word_13 != 60465)
        & (stat_word_13 != 60593)
    )
    print("filter3:", np.count_nonzero(filter3))
    filter4 = (
        # (stat_word_16 != 14372)
        # (stat_word_16 != 35584)
        (stat_word_16 != 36032)
        & (stat_word_16 != 52992)
        & (stat_word_16 != 53056)
    )
    print("filter4:", np.count_nonzero(filter4))
    filter5 = (
        (lvdt_stat_a != 4)
        & (lvdt_stat_a != 17)
        & (lvdt_stat_a != 26)
        & (lvdt_stat_a != 31)
        & (lvdt_stat_a != 49)
    )
    print("filter5:", np.count_nonzero(filter5))
    filter6 = (
        (lvdt_stat_b != -127)
        & (lvdt_stat_b != -124)
        & (lvdt_stat_b != -79)
        & (lvdt_stat_b != 82)
        & (lvdt_stat_b != 83)
        & (lvdt_stat_b != 96)
        & (lvdt_stat_b != 111)
    )
    print("filter6:", np.count_nonzero(filter6))

    return filter1 & filter2 & filter3 & filter4 & filter5 & filter6
