import numpy as np
import globals as g
import matplotlib.pyplot as plt
from astropy.io import fits


def etfunction(channel, adds_per_group, samprate):
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

    fixed_gain = 31.0 * 1.3823  # Preamp fixed gain
    bessel_gain = 1.2255 * 1.9099  # Bessel DC gain
    dcgain = bessel_gain * fixed_gain

    bes3db = 100.0  # Bessel corner freq
    tau = [0.00647, 0.03619, 0.00722, 0.04022]  # RH, RL, LH, LL treble boost

    dfhz = samprate / adds_per_group / g.IFG_SIZE
    zxfer = np.zeros(g.SPEC_SIZE, dtype=complex)
    for k in range(g.SPEC_SIZE):
        freqhz = k * dfhz

        zdigfil = _digfltr(freqhz, channel, samprate)
        zsmooth = _compress(freqhz, adds_per_group, samprate)
        zdigital = zdigfil * zsmooth

        zbesl = _bessel(freqhz, bes3db)
        ztboost = _tboost(freqhz, tau[channel])
        zdcblock = _dcblock(freqhz)
        zanalog = zbesl * ztboost * zdcblock

        zxfer[k] = dcgain * zanalog * zdigital

    ztrans = -zxfer
    return ztrans


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
    zsmooth = 1.0
    if pif < 0.001:
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

    etf = etfunction(channel, adds_per_group, samprate)
    print(etf)
    plt.plot(etf.real)
    plt.plot(etf.imag)
    # plt.show()

    # compare with the published etf
    fits_data = fits.open(f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_RHSS.FITS")[1].data
    fits_etf = fits_data["RELEX_GA"][0] + 1j * fits_data["IELEX_GA"][0]
    print(fits_etf)

    plt.plot(
        # np.arange(7, 7 + fits_etf.real.size),
        fits_etf.real,
        color="black",
    )
    plt.plot(
        # np.arange(7, 7 + fits_etf.imag.size),
        fits_etf.imag,
        color="black",
        linestyle="--",
    )
    plt.show()
