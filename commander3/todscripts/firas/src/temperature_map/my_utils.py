import numpy as np
import matplotlib.pyplot as plt

from utils.fut import apod_recnuml, get_recnum
from utils.frd import apodl, elex_transfcnl
from utils.config import gen_nyquistl


def clean_ifg(ifg, mtm_length, mtm_speed, channel, adds_per_group):
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
    ifg = np.fft.rfft(ifg)

    # etf
    etfl_all = elex_transfcnl(samprate=681.43, nfreq=len(ifg))
    erecno = int(get_recnum(fake_it, mtm_speed, channel, upmode, adds_per_group))
    etf = etfl_all[erecno, :]

    scan_mode = 0  # SS
    frec = 4 * (channel % 2) + scan_mode
    fnyq_icm = gen_nyquistl(
        "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
    )["icm"][frec]
    fac_etendu = 1.5  # nathan's pipeline
    fac_adc_scale = 204.75  # nathan's pipeline
    spec_norm = fnyq_icm * fac_etendu * fac_adc_scale

    ifg = ifg / (etf * spec_norm)

    plt.plot(np.abs(ifg))
    plt.savefig("../../output/tests/test.png")
