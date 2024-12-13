import numpy as np
import matplotlib.pyplot as plt

from utils.fut import apod_recnuml, get_recnum
from utils.frd import apodl, elex_transfcnl

from astropy.io import fits


def calculate_dc_response(bol_cmd_bias, bol_volt, fits_data):
    rscale = 1.0e-7

    cmd_bias = np.double(bol_cmd_bias) / 25.5

    Jo = fits_data[1].data["BOLPARM8"][0]
    Jg = fits_data[1].data["BOLPARM9"][0]
    V = (bol_volt - Jo) / Jg

    RL = 4.0e7
    R = RL * V / (cmd_bias - V)

    Tbol = fits_data[1].data["BOLOM_B2"][0]
    T0 = fits_data[1].data["BOLPARM2"][0]

    R0 = fits_data[1].data["BOLPARM_"][0]
    rho = fits_data[1].data["BOLPARM5"][0]
    X = V * rho
    Y = R / R0 / X

    SQ = np.log(Y * Tbol * np.sinh(X / Tbol))
    Tbol = T0 / SQ**2

    SQ = np.log(Y * Tbol * np.sinh(X / Tbol))

    Tbol = T0 / (SQ * SQ)

    G1 = fits_data[1].data["BOLPARM3"][0]
    beta = fits_data[1].data["BOLPARM4"][0]
    G = G1 * Tbol**beta

    H = Tbol / X * np.tanh(X / Tbol)

    DT = 1.0 / H - 1.0 - 0.5 * np.sqrt(T0 / Tbol)

    Z = (G * Tbol * R + DT * V**2) / (G * Tbol * R / H - DT * V**2)

    S0 = rscale * R * (Z - H) / (V * (Z * R / RL + 1.0) * (H + 1.0))

    return S0


def clean_ifg(ifg, mtm_length, mtm_speed, channel, adds_per_group, bol_cmd_bias, bol_volt, fits_data, fnyq, frec):
    fake_it = 0
    upmode = 4

    # subtract dither
    ifg = ifg - np.median(ifg)

    # apodize
    sm = 2 * mtm_length + mtm_speed

    arecno = int(apod_recnuml(channel, sm, fake_it, upmode, adds_per_group, 0))
    apodl_all = apodl()
    apod = apodl_all[arecno, :]

    ifg = ifg * apod

    # roll
    peak_pos = 360
    ifg = np.roll(ifg, -peak_pos)

    # fft
    spec = np.fft.rfft(ifg)

    # etf
    etfl_all = elex_transfcnl(samprate=681.43, nfreq=len(spec))
    erecno = int(get_recnum(fake_it, mtm_speed, channel, upmode, adds_per_group))
    etf = etfl_all[erecno, :]

    
    fnyq_icm = fnyq["icm"][frec]
    fac_etendu = 1.5  # nathan's pipeline
    fac_adc_scale = 204.75  # nathan's pipeline
    spec_norm = fnyq_icm * fac_etendu * fac_adc_scale

    spec = spec / (etf * spec_norm)

    # bolometer function
    tau = fits_data[1].data['TIME_CON'][0]

    # fcc_spec_length = 321
    spec_len = len(ifg) // 2 + 1
    fnyq_hz = fnyq['hz'][frec]
    dw = 2.0 * np.pi * fnyq_hz / spec_len
    # afreq = np.arange(fcc_spec_length) * dw # had to change because i'm not padding for now
    afreq = np.arange(spec_len) * dw

    S0 = calculate_dc_response(bol_cmd_bias, bol_volt, fits_data)
    B = 1.0 + 1j * tau * afreq

    spec = B * spec / S0

    # optical transfer function
    otf = fits_data[1].data["RTRANSFE"][0] + 1j * fits_data[1].data["ITRANSFE"][0]
    otf = otf[np.abs(otf) > 0]

    spec = spec[:len(otf)] / otf

    return spec