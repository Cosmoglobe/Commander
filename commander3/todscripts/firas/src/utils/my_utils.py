import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

import constants

channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}


def ghz_to_icm(ghz):
    """
    Converts GHz to inverse cm.
    """
    return ghz * 1e9 / constants.c


def generate_frequencies(channel, mode, nfreq=None):
    """
    Generates an array with the frequencies in GHz for the given channel and mode.

    Parameters
    ----------
    channel : str
        The channel to generate frequencies for. Can be "lh", "ll", "rh", or "rl".
    mode : str or int
        The mode to generate frequencies for. Can be "ss" or "lf" for str or 0 or 3 for int.

    Returns
    -------
    f_ghz : np.ndarray
        An array with the frequencies in GHz.
    """

    # check if channel is str or int
    if isinstance(channel, int):
        channel_str = list(constants.channels.keys())[
            list(constants.channels.values()).index(channel)
        ]
    elif isinstance(channel, str):
        channel_str = channel
    else:
        raise ValueError("Channel must be either int or str")

    nu0 = {"ss": 68.020812, "lf": 23.807283}
    dnu = {"ss": 13.604162, "lf": 3.4010405}
    nf = {
        "lh_ss": 210,
        "ll_lf": 182,
        "ll_ss": 43,
        "rh_ss": 210,
        "rl_lf": 182,
        "rl_ss": 43,
    }

    if not (mode == "lf" and (channel_str == "lh" or channel_str == "rh")):
        if nfreq == None:
            nfreq = nf[f"{channel_str}_{mode}"]
            f_ghz = np.linspace(
                nu0[mode],
                nu0[mode] + dnu[mode] * (nfreq - 1),
                nfreq,
            )
        else:
            f_ghz = np.linspace(0, dnu[mode] * (nfreq - 1), nfreq)
    else:
        raise ValueError("Invalid channel and mode combination")

    return f_ghz


def get_afreq(mtm_speed, channel, nfreq=None):
    """
    Returns the afreq for the given channel and mode.

    Parameters
    ----------
    mtm_speed : int
        The speed of the mirror transport mechanism. Can be 0 or 1.
    channel : int or str
        The channel to get the afreq for. Can be 0, 1, 2, or 3 for int or "lh", "ll", "rh", "rl" for str.
    nfreq : int, optional
        The number of frequencies to generate. If None, defaults to the size of the channel.
    """

    # check if channel is str or int
    if isinstance(channel, int):
        channel_str = list(constants.channels.keys())[
            list(constants.channels.values()).index(channel)
        ]
    elif isinstance(channel, str):
        channel_str = channel
    else:
        raise ValueError("Channel must be either int or str")

    speed = constants.speed[mtm_speed]

    mode_str = "ss" if mtm_speed == 0 else "lf"
    f_ghz = generate_frequencies(channel_str, mode_str, nfreq)
    f_icm = ghz_to_icm(f_ghz)

    afreq = speed * f_icm

    return afreq


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


def save_mollview(m, f_ghz, save_path, max_amp=200, norm="linear"):
    for freqi, frequency in enumerate(f_ghz):
        hp.mollview(
            m[:, freqi],
            title=f"{int(frequency):04d} GHz",
            min=1,
            max=max_amp,
            unit="MJy/sr",
            norm=norm,
        )
        plt.savefig(f"{save_path}{int(frequency):04d}.png")
        plt.close()
