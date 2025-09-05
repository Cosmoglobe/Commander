import globals as g
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from utils import frd, fut


def etfunction(channel, adds_per_group, samprate, nui=None):
    """
    This function is adapted from Nathan Miller's python FIRAS pipeline in order to generate the corresponding electronics transfer function, according to:

    Parameters
    ----------
    channel : int or str
        0: RH, 1: RL, 2: LH, 3: LL
    adds_per_group : int
        The number of points added for each point of the interferogram.
    """
    if isinstance(channel, str):
        channel = g.CHANNELS[channel]

    dfhz = samprate / adds_per_group / g.IFG_SIZE
    # dfhz = samprate / adds_per_group / 640 # how to get exactly the same as theirs

    if nui is None:
        ampmax = -1
        ztrans = np.zeros((adds_per_group.shape[0], g.SPEC_SIZE), dtype=complex)
        for k in range(g.SPEC_SIZE):
            freqhz = k * dfhz

            zxfer = compute_etf_per_freq(freqhz, channel, samprate, adds_per_group)

            ampl = np.sqrt(np.abs(zxfer) ** 2)
            ampmax = np.maximum(ampl, ampmax)

            # if ampl < ampmax / 1000.0:
            ztrans[:, k][ampl < ampmax / 1000.0] = 100000.0
            # else:
            ztrans[:, k][ampl >= ampmax / 1000.0] = -zxfer[ampl >= ampmax / 1000.0]
    else:
        freqhz = nui * dfhz

        zxfer = compute_etf_per_freq(freqhz, channel, samprate, adds_per_group)

        ampmax = -1
        ampl = np.sqrt(np.abs(zxfer) ** 2)
        ampmax = max(ampl, ampmax)

        if ampl < ampmax / 1000.0:
            ztrans = 100000.0
        else:
            ztrans = -zxfer

    return ztrans


def compute_etf_per_freq(freqhz, channel, samprate, adds_per_group):
    zdigfil = _digfltr(freqhz, channel, samprate)
    zsmooth = _compress(freqhz, adds_per_group, samprate)
    zdigital = zdigfil * zsmooth

    fixed_gain = 31.0 * 1.3823  # Preamp fixed gain
    bessel_gain = 1.2255 * 1.9099  # Bessel DC gain
    dcgain = bessel_gain * fixed_gain

    bes3db = 100.0  # Bessel corner freq
    tau = [0.00647, 0.03619, 0.00722, 0.04022]  # RH, RL, LH, LL treble boost

    zbesl = _bessel(freqhz, bes3db)
    ztboost = _tboost(freqhz, tau[channel])
    zdcblock = _dcblock(freqhz)
    zanalog = zbesl * ztboost * zdcblock

    return dcgain * zanalog * zdigital


def _digfltr(freqhz, ichan, samplrate):
    """
    Adapted from Nathan's function, removing the micromode parameter since for our data all micromode = 0.
    """
    zi = -2j * np.pi
    z = np.exp(zi * freqhz / samplrate)
    zdigfil = 1.0
    if ichan == 1 or ichan == 3:
        # zdigfil = (1.0 + z*z)**2/(8.0 - 12.5*z*z + 5.0*z**4)/8.0
        zdigfil = (
            (1.0 + 2.0 * z * (1.0 + z + z * z) + z**4) / (8.0 - 8.0 * z * z + z**4)
        ) / 8.0
    else:
        zdigfil = (1.0 + z) ** 2 / (8.0 - 11.0 * z + 5.0 * z * z) / 2.0

    return zdigfil


def _compress(fhz, ncompress, samplrate):
    pif = np.pi * fhz / samplrate
    zsmooth = np.ones_like(fhz)
    if np.all(pif < 0.001):
        return zsmooth
    ampl = np.sin(pif * ncompress) / np.sin(pif)
    zphase = pif * (1 - ncompress) * 1j
    zsmooth = ampl / ncompress * np.exp(zphase)

    return zsmooth


def _bessel(fhz, bes3db):
    zfb = 1j * fhz / bes3db
    zfb2 = zfb * zfb
    zfb3 = zfb * zfb2
    zfb4 = zfb * zfb3
    zfb5 = zfb * zfb4

    zbesl = 1.0 / (
        1.0
        + 2.4275 * zfb
        + 2.6189 * zfb2
        + 1.5894 * zfb3
        + 0.5511 * zfb4
        + 0.0892 * zfb5
    )

    return zbesl


def _tboost(fhz, tau):
    return 1.0 + 2 * np.pi * 1j * fhz * tau


def _dcblock(fhz):
    zs = 3.2j * 2 * np.pi * fhz
    return (zs / (1.0 + zs)) ** 5 * zs / (2.0 + zs)


if __name__ == "__main__":
    channel = 0  # RH
    adds_per_group = 3  # random
    samprate = 681.43  # from fex_samprate

    etf = etfunction(channel, np.array([adds_per_group]), samprate)[0]
    print(etf)
    plt.plot(etf.real)
    plt.plot(etf.imag)
    # plt.show()

    # compare with the published etf
    fits_data = fits.open(f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_RHSS.FITS")[1].data
    fits_etf = fits_data["RELEX_GA"][0] + 1j * fits_data["IELEX_GA"][0]
    print(fits_etf)

    plt.plot(
        np.arange(5, 5 + fits_etf.real.size),
        fits_etf.real,
        color="black",
    )
    plt.plot(
        np.arange(5, 5 + fits_etf.imag.size),
        fits_etf.imag,
        color="black",
        linestyle="--",
    )

    # compare with etf from pipeline
    etfs = frd.elex_transfcnl(samprate=681.43, nfreq=257)
    erecno = fut.get_recnum(0, channel, adds_per_group).astype(np.int32)
    plt.plot(etfs[erecno, :].real, color="red")
    plt.plot(etfs[erecno, :].imag, color="red", linestyle="--")

    plt.show()
