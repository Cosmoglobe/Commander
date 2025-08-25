import numpy as np

import globals as g
from calibration import bolometer
from utils import frd, fut
from utils import my_utils as utils


def clean_ifg(
    ifg,
    channel,
    mode,
    gain,
    sweeps,
    apod,
):
    median_ifg = np.expand_dims(np.median(ifg, axis=1), axis=-1)

    # subtract dither
    ifg = ifg - median_ifg

    # Ensure gain and sweeps are reshaped for broadcasting
    gain = np.expand_dims(gain, axis=-1)
    sweeps = np.expand_dims(sweeps, axis=-1)

    ifg = ifg / gain / sweeps
    ifg = ifg * apod

    # roll
    peak_pos = g.PEAK_POSITIONS[f"{channel}_{mode}"]
    ifg = np.roll(ifg, -peak_pos, axis=1)

    return ifg


def unclean_ifg(
    ifg,
    channel,
    mode,
    gain,
    sweeps,
    apod,
):
    gain = np.expand_dims(gain, axis=-1)
    sweeps = np.expand_dims(sweeps, axis=-1)

    ifg = ifg * gain * sweeps

    peak_pos = g.PEAK_POSITIONS[f"{channel}_{mode}"]

    # leaving this here for now because
    if mode == "ss" and channel[1] == "l":
        peak_pos = peak_pos + 1
    elif mode == "lf":
        peak_pos = peak_pos - 1
    ifg = np.roll(ifg, peak_pos, axis=1)

    ifg = ifg / apod
    ifg[:, apod < 0.3] = np.nan

    return ifg


def ifg_to_spec(
    ifg,
    channel,
    mode,
    adds_per_group,
    bol_cmd_bias,
    bol_volt,
    fnyq_icm,
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
    apod,
    gain=1,
    sweeps=1,
):
    """
    Inputs are the same as spec_to_ifg, except for the ifg argument, which is the interferogram to be converted.

    ifg is the interferogram to be converted. Its frequency values are set specifically by the mode of operation.
    channel: channel abbreviation: "lh", "ll", "rh" or "rl".
    mode: mode of operation: "ss" or "lf".
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
    """

    ifg = clean_ifg(ifg, channel, mode, gain, sweeps, apod)

    spec = np.fft.rfft(ifg)
    # freqs = np.fft.rfftfreq(constants.ifg_size, 1 / (fnyq_hz * 2)) TODO: find out why this doesn't work

    mtm_speed = 0 if mode[1] == "s" else 1

    # etf from the pipeline
    etfl_all = frd.elex_transfcnl(samprate=681.43, nfreq=len(spec[0]))
    erecno = fut.get_recnum(mtm_speed, utils.channels[channel], adds_per_group).astype(
        np.int32
    )
    etf = etfl_all[erecno, :]

    spec_norm = fnyq_icm * g.FAC_ETENDU * g.FAC_ADC_SCALE

    spec = spec / etf
    spec = spec / spec_norm

    afreq = utils.get_afreq(mtm_speed, utils.channels[channel], 257)

    S0 = bolometer.calculate_dc_response(
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

    tau = bolometer.calculate_time_constant(
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

    spec = spec * g.FAC_ERG_TO_MJY

    return afreq, spec


def spec_to_ifg(
    spec,
    channel,
    mode,
    adds_per_group,
    bol_cmd_bias,
    bol_volt,
    Tbol,
    gain,
    sweeps,
    otf,
    apod,
    fnyq_icm,
    R0,
    T0,
    G1,
    beta,
    rho,
    C1,
    C3,
    Jo,
    Jg,
):
    """
    Converts spectrum to interferogram using the pipeline's etf and otf. Expects the spectrum to be in units of MJy/sr.

    The arguments are the same as ifg_to_spec, except for the spec argument, which is the spectrum to be converted.

    spec is spectrum to be converted in units of MJy/sr. Its frequency values are set specifically by the mode of operation.

    Specific to the spectrum/ifg:
        mtm_speed: speed of the mirror transport mechanism, 0 or 1
        channel: channel number, left low, left high, right low, right high
        adds_per_group: number of interferograms added per group.
        bol_cmd_bias: bias command to the bolometer. Found in fdq_eng/en_stat/bol_cmd_bias, (589069, 4) int8 array
        bol_volt: Found in fdq_eng/en_analog/group1/bol_volt, (589069, 4) float32 array
        gain: Found in fdq_sdf_??/sci_head/gain, (589069,) int16 array, values from -1 to 7 inclusive.
        sweeps: Found in fdq_sdf_??/sci_head/sc_head11, (589069,) int16 array, values from -1 to 16 inclusive.


    These parameters all come directly from the calibration model files. We should try to fit these (or in apod, tune).
    otf: optical transfer function, from the published calibration model.
    apod: the apodization function, length 512 array, specially tuned by FIRAS team.
    R0, T0, G1, beta, rho, C1, C3, Jo, Jg,: parameters for the bolometers, from the published calibration model.
    As an example:
    R0 = fits_data[1].data["BOLPARAM_"][0]
    T0 = fits_data[1].data["BOLPARAM2"][0]
    G1 = fits_data[1].data["BOLPARAM3"][0]
    beta = fits_data[1].data["BOLPARAM4"][0]
    rho = fits_data[1].data["BOLPARAM5"][0]
    C1 = fits_data[1].data["BOLPARAM6"][0]
    C3 = fits_data[1].data["BOLPARAM7"][0]
    Jo = fits_data[1].data["JO"][0]
    Jg = fits_data[1].data["JG"][0]
    Typical values for LLSS:
    R0=11.938669
    T0=371.645
    G1=2.3263578
    beta=1.172
    rho=2.3263578
    C1=5.8665056e-10
    C3=2.428e-10
    Jo=1.662232
    Jg=0.92
    Tbol=1.5309418
    """

    mtm_speed = 0 if mode[1] == "s" else 1

    if mtm_speed == 0:
        cutoff = 5
    else:
        cutoff = 7

    spec_r = np.zeros((len(spec), 257))
    if spec.shape[1] == 257:
        spec_r = spec
    else:
        spec_r[:, cutoff : (len(spec[0]) + cutoff)] = spec

    spec_norm = fnyq_icm * g.FAC_ETENDU * g.FAC_ADC_SCALE

    S0 = bolometer.calculate_dc_response(
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

    tau = bolometer.calculate_time_constant(
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
    etfl_all = frd.elex_transfcnl(samprate=681.43, nfreq=257)
    erecno = fut.get_recnum(mtm_speed, utils.channels[channel], adds_per_group).astype(
        np.int32
    )
    etf = etfl_all[erecno, :]

    afreq = utils.get_afreq(mtm_speed, utils.channels[channel], 257)
    B = 1.0 + 1j * tau[:, np.newaxis] * afreq[np.newaxis, :]

    spec_r = spec_r / g.FAC_ERG_TO_MJY
    spec_r = spec_r * otf
    spec_r = S0[:, np.newaxis] * spec_r / B
    spec_r = spec_r * spec_norm
    spec_r = spec_r * etf

    ifg = np.fft.irfft(
        spec_r,
    )
    ifg = unclean_ifg(
        ifg=ifg, channel=channel, mode=mode, gain=gain, sweeps=sweeps, apod=apod
    )

    return ifg
