import healpy as hp
import numpy as np
from numba import prange
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
from utils.frd import elex_transfcnl
from utils.fut import get_recnum


def calculate_dc_response(bol_cmd_bias, bol_volt, Jo, Jg, Tbol, rho, R0, T0, beta, G1):
    """
    Taken largely from calc_responsivity in fsl.
    """

    rscale = 1.0e-7

    cmd_bias = bol_cmd_bias.astype(
        "double"
    )

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


def calculate_time_constant(
    C3, Tbol, C1, G1, beta, bol_volt, Jo, Jg, bol_cmd_bias, rho, T0
):
    """
    Taken largely from calc_responsivity in fsl.
    """
    # print(
    #     f"shapes of arguments:\nC3: {C3.shape}\nTbol: {Tbol.shape}\nC1: {C1.shape}\nG1: {G1.shape}\nbeta: {beta.shape}\nbol_volt: {bol_volt.shape}\nJo: {Jo.shape}\nJg: {Jg.shape}\nbol_cmd_bias: {bol_cmd_bias.shape}\nrho: {rho.shape}\nT0: {T0.shape}"
    # )
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
def clean_ifg(
    ifg,
    mtm_length,
    mtm_speed,
    # channel,
    # adds_per_group,
    gain,
    sweeps,
    apod,
):
    # print("ifg shape before subtraction:", ifg.shape)

    median_ifg = np.expand_dims(np.median(ifg, axis=1), axis=-1)

    # subtract dither
    ifg = ifg - median_ifg

    # print("ifg shape after subtraction:", ifg.shape)

    # Ensure gain and sweeps are reshaped for broadcasting
    gain = np.expand_dims(gain, axis=-1)
    sweeps = np.expand_dims(sweeps, axis=-1)

    ifg = ifg / gain / sweeps

    # print("ifg shape after division:", ifg.shape)

    # apodize
    sm = 2 * mtm_length + mtm_speed

    ifg = ifg * apod

    # print("ifg shape after apodization:", ifg.shape)

    # roll
    peak_pos = 360
    ifg = np.roll(ifg, -peak_pos)

    # print("ifg shape after roll:", ifg.shape)

    return ifg


def unclean_ifg(
    ifg,
    gain,
    sweeps,
    apod,
):
    gain = np.expand_dims(gain, axis=-1)
    sweeps = np.expand_dims(sweeps, axis=-1)

    ifg = ifg * gain * sweeps

    peak_pos = 360
    ifg = np.roll(ifg, peak_pos)

    ifg = ifg / apod
    ifg[:, apod < 0.1] = np.nan

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
    C3,
    C1,
    mtm_length,
    gain,
    sweeps,
    apod,
):
    '''
    Inputs are the same as spec_to_ifg, except for the ifg argument, which is the interferogram to be converted.

    ifg is the interferogram to be converted. Its frequency values are set specifically by the mode of operation.
    mtm_speed: speed of the mirror transport mechanism, 0 or 1
    channel: channel number, left low, left high, right low, right high
    adds_per_group: number of interferograms added per group.
    bol_cmd_bias: bias command to the bolometer. Found in fdq_eng/en_stat/bol_cmd_bias, (589069, 4) int8 array
    bol_volt: Found in fdq_eng/en_analog/group1/bol_volt, (589069, 4) float32 array
    gain: Found in fdq_sdf_??/sci_head/gain, (589069,) int16 array, values from -1 to 7 inclusive.
    sweeps: Found in fdq_sdf_??/sci_head/sc_head11, (589069,) int16 array, values from -1 to 16 inclusive.

    These parameters all come directly from the calibration model files. We should try to fit these (or in apod, tune).
    apod: the apodization function, length 512 array, specially tuned by FIRAS team.
    otf: optical transfer function, from the published calibration model.
    Jo, Jg, Tbol, rho, R0, T0, beta, G1, C3, C1: parameters for the bolometers, from the published calibration model.
    As an example:
    Jo = data[1].data['BOLPARM8'][0]
    Jg = data[1].data['BOLPARM9'][0]
    Tbol = data[1].data['BOLOM_B2'][0]
    rho = data[1].data['BOLPARM5'][0]
    R0 = data[1].data['BOLPARM_'][0]
    T0 = data[1].data['BOLPARM2'][0]
    beta = data[1].data['BOLPARM4'][0]
    G1 = data[1].data['BOLPARM5'][0]
    C3 = data[1].data['BOLPARM7'][0]
    C1 = data[1].data['BOLPARM6'][0]
    Typical values for LLSS:
    Jo=1.662232
    Jg=0.92
    Tbol=1.5309418
    rho=2.3263578
    R0=11.938669
    T0=371.645
    beta=1.172
    G1=2.3263578
    C3=2.428e-10
    C1=5.8665056e-10 

    mtm_length: length of the mirror transport mechanism, 0 or 1
    gain: gain of the bolometer, from the calibration model
    sweeps: number of sweeps, from the calibration model
    apod: apodization function, from the calibration model
    '''

    ifg = clean_ifg(ifg, mtm_length, mtm_speed, gain, sweeps, apod)

    spec = np.fft.rfft(ifg)

    # print("spec after rfft:", spec)

    # etf from the pipeline
    etfl_all = elex_transfcnl(samprate=681.43, nfreq=len(spec[0]))
    erecno = get_recnum(mtm_speed, channel, adds_per_group).astype(np.int32)
    etf = etfl_all[erecno, :]


    fac_etendu = 1.5  # nathan's pipeline
    fac_adc_scale = 204.75  # nathan's pipeline
    spec_norm = fnyq_icm * fac_etendu * fac_adc_scale

    spec = spec / etf
    spec = spec / spec_norm

    spec_len = len(ifg[0]) // 2 + 1
    dw = 2.0 * np.pi * fnyq_hz / spec_len
    afreq = np.arange(spec_len) * dw

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

    B = 1.0 + 1j * tau[:, np.newaxis] * afreq[np.newaxis, :]

    spec = spec / S0[:, np.newaxis] * B

    if mtm_speed == 0:
        cutoff = 5
    else:
        cutoff = 7

    spec[:, cutoff : (len(otf) + cutoff)] = spec[:, cutoff : (len(otf) + cutoff)] / otf
    spec[:, :cutoff] = 0
    spec[:, (len(otf) + cutoff) :] = 0

    fac_icm_ghz = 29.9792458
    fac_erg_to_mjy = 1.0e8 / fac_icm_ghz

    spec = spec * fac_erg_to_mjy

    return afreq, spec


def spec_to_ifg(
    spec,
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
    C3,
    C1,
    gain,
    sweeps,
    apod,
):
    '''
    Converts spectrum to interferogram using the pipeline's etf and otf. Expects the spectrum to be in units of MJy/sr.

    The arguments are the same as ifg_to_spec, except for the spec argument, which is the spectrum to be converted.

    spec is spectrum to be converted in units of MJy/sr. Its frequency values are set specifically by the mode of operation.

    mtm_speed: speed of the mirror transport mechanism, 0 or 1
    channel: channel number, left low, left high, right low, right high
    adds_per_group: number of interferograms added per group.
    bol_cmd_bias: bias command to the bolometer. Found in fdq_eng/en_stat/bol_cmd_bias, (589069, 4) int8 array
    bol_volt: Found in fdq_eng/en_analog/group1/bol_volt, (589069, 4) float32 array
    gain: Found in fdq_sdf_??/sci_head/gain, (589069,) int16 array, values from -1 to 7 inclusive.
    sweeps: Found in fdq_sdf_??/sci_head/sc_head11, (589069,) int16 array, values from -1 to 16 inclusive.


    These parameters all come directly from the calibration model files. We should try to fit these (or in apod, tune).
    apod: the apodization function, length 512 array, specially tuned by FIRAS team.
    otf: optical transfer function, from the published calibration model.
    Jo, Jg, Tbol, rho, R0, T0, beta, G1, C3, C1: parameters for the bolometers, from the published calibration model.
    As an example:
    Jo = data[1].data['BOLPARM8'][0]
    Jg = data[1].data['BOLPARM9'][0]
    Tbol = data[1].data['BOLOM_B2'][0]
    rho = data[1].data['BOLPARM5'][0]
    R0 = data[1].data['BOLPARM_'][0]
    T0 = data[1].data['BOLPARM2'][0]
    beta = data[1].data['BOLPARM4'][0]
    G1 = data[1].data['BOLPARM5'][0]
    C3 = data[1].data['BOLPARM7'][0]
    C1 = data[1].data['BOLPARM6'][0]
    Typical values for LLSS:
    Jo=1.662232
    Jg=0.92
    Tbol=1.5309418
    rho=2.3263578
    R0=11.938669
    T0=371.645
    beta=1.172
    G1=2.3263578
    C3=2.428e-10
    C1=5.8665056e-10
    '''
    if mtm_speed == 0:
        cutoff = 5
    else:
        cutoff = 7

    #spec_r = np.zeros((len(spec), 257))
    #spec_r[:, cutoff : (len(otf) + cutoff)] = spec
    spec_r = spec

    # Defining constants, getting data model

    fac_icm_ghz = 29.9792458
    fac_erg_to_mjy = 1.0e8 / fac_icm_ghz

    fac_etendu = 1.5  # nathan's pipeline
    fac_adc_scale = 204.75  # nathan's pipeline
    spec_norm = fnyq_icm * fac_etendu * fac_adc_scale

    spec_len = len(spec_r[0])
    dw = 2.0 * np.pi * fnyq_hz / spec_len
    afreq = np.arange(spec_len) * dw

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
    etfl_all = elex_transfcnl(samprate=681.43, nfreq=len(spec_r[0]))
    erecno = get_recnum(mtm_speed, channel, adds_per_group).astype(np.int32)
    etf = etfl_all[erecno, :]

    B = 1.0 + 1j * tau[:, np.newaxis] * afreq[np.newaxis, :]

    spec_r = spec_r / fac_erg_to_mjy

    spec_r[:, cutoff : len(otf) + cutoff] = spec_r[:, cutoff : len(otf) + cutoff] * otf

    spec_r = S0[:, np.newaxis] * spec_r / B

    spec_r = spec_r * spec_norm

    spec_r = spec_r * etf

    ifg = np.fft.irfft(spec_r)

    ifg = unclean_ifg(ifg, gain, sweeps, apod)

    return ifg


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


def filter_junk(
    stat_word_5, stat_word_9, stat_word_13, stat_word_16, lvdt_stat_a, lvdt_stat_b
):
    filter1 = (stat_word_5 != 16641) & (stat_word_5 != 17921) & (stat_word_5 != 17217)
    print("filter1:", np.count_nonzero(filter1))
    filter2 = stat_word_9 != 15414
    print("filter2:", np.count_nonzero(filter2))
    filter3 = (
        (stat_word_13 != 17345)
        & (stat_word_13 != 17393)
        & (stat_word_13 != 19585)  # makes the hole?
        & (stat_word_13 != 25201)
        & (stat_word_13 != 25585)
        & (stat_word_13 != 26945)
        & (stat_word_13 != 27073)
        & (stat_word_13 != 27649)
        & (stat_word_13 != 27697)
        # & (stat_word_13 != 27777)
        & (stat_word_13 != 60465)
        & (stat_word_13 != 60593)
    )
    print("filter3:", np.count_nonzero(filter3))
    filter4 = (
        (stat_word_16 != 0)
        & (stat_word_16 != 14372)  # makes the hole?
        & (stat_word_16 != 35584)  # makes the hole?
        & (stat_word_16 != 36032)
        & (stat_word_16 != 52992)
        & (stat_word_16 != 53056)
    )
    print("filter4:", np.count_nonzero(filter4))
    filter5 = (
        (lvdt_stat_a != 1)
        & (lvdt_stat_a != 2)
        & (lvdt_stat_a != 4)
        & (lvdt_stat_a != 6)
        & (lvdt_stat_a != 13)
        & (lvdt_stat_a != 17)
        & (lvdt_stat_a != 21)
        & (lvdt_stat_a != 24)
        & (lvdt_stat_a != 26)
        & (lvdt_stat_a != 31)
        & (lvdt_stat_a != 32)
        & (lvdt_stat_a != 35)
        & (lvdt_stat_a != 47)
        & (lvdt_stat_a != 49)
    )
    print("filter5:", np.count_nonzero(filter5))
    filter6 = (
        (lvdt_stat_b != -127)
        & (lvdt_stat_b != -124)
        & (lvdt_stat_b != -108)
        & (lvdt_stat_b != -83)
        & (lvdt_stat_b != -79)
        & (lvdt_stat_b != -78)
        & (lvdt_stat_b != -71)
        & (lvdt_stat_b != -70)
        & (lvdt_stat_b != 75)
        & (lvdt_stat_b != 79)
        & (lvdt_stat_b != 82)
        & (lvdt_stat_b != 83)
        & (lvdt_stat_b != 85)
        & (lvdt_stat_b != 86)
        & (lvdt_stat_b != 87)
        & (lvdt_stat_b != 91)
        & (lvdt_stat_b != 93)
        & (lvdt_stat_b != 95)
        & (lvdt_stat_b != 96)
        & (lvdt_stat_b != 99)
        & (lvdt_stat_b != 100)
        & (lvdt_stat_b != 103)
        & (lvdt_stat_b != 107)
        & (lvdt_stat_b != 111)
        & (lvdt_stat_b != 121)
    )
    print("filter6:", np.count_nonzero(filter6))

    return filter1 & filter2 & filter3 & filter4 & filter5 & filter6


def tune_pointing(gal_lon, gal_lat, gmt, mtm_length, mtm_speed, offset=0):
    """
    Interpolate the pointing to different offsets depending on the operating mode.
    Doing this on the main script for now because it is mostly an experiment and doesn't need a lot of precision.
    """
    if offset == 0:
        return gal_lon, gal_lat

    times = {"00": 55.36, "10": 44.92, "01": 39.36, "11": 31.76}

    gal_vec = hp.pixelfunc.ang2vec(gal_lon, gal_lat, lonlat=True)
    new_gal_vec = np.empty_like(gal_vec)
    # gmt = np.array(gmt).timestamp().values

    gmt = np.array(gmt, dtype="datetime64[s]")
    gmt = (gmt - np.datetime64("1970-01-01T00:00:00Z")) / np.timedelta64(1, "s")

    for i in range(len(gmt)):
        target_time = (
            gmt[i] + offset * times[f"{int(mtm_length[i])}{int(mtm_speed[i])}"]
        )
        # print(f"target_time: {target_time}")

        # using two minutes around the target time
        lower_bound = target_time - 120
        upper_bound = target_time + 120

        start_id = np.searchsorted(gmt, lower_bound, side="left")
        end_id = np.searchsorted(gmt, upper_bound, side="right")
        idx = np.arange(start_id, end_id)

        # print(f"idx: {idx}")

        if len(idx) > 1 and gmt[start_id] <= target_time <= gmt[end_id - 1]:
            new_gal_vec[i] = interpolate.interp1d(
                gmt[idx], gal_vec[idx], axis=0, kind="linear"
            )(target_time)
        else:
            # print(f"no interpolation for {i}")
            new_gal_vec[i] = np.nan

    return new_gal_vec
