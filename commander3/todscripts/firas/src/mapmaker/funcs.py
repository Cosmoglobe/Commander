import numpy as np


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


def dust(nu, tau, nu0, beta, T):
    """
    Returns the dust SED in units of MJy/sr.
    Input frequencies in GHz and temperature in K.

    Parameters
    ----------
    nu : float
        Frequency in GHz.
    nu0 : float
        Reference frequency in GHz.
    beta : float
        Dust spectral index.
    T : float
        Dust temperature in K.

    Returns
    -------
    float
        Dust SED in MJy/sr.
    """
    h = 6.62607015e-34  # J Hz−1
    k_B = 1.380649e-23  # J K-1

    nu = np.array(nu)
    nu0 = np.array(nu0)

    # gamma = h / (k_B * T)
    # sed = np.array(
    #     (nu / nu0) ** (beta + 1)
    #     * (np.exp(gamma * nu0 * 1e6) - 1)
    #     / (np.exp(gamma * nu * 1e6) - 1)
    # )
    sed = tau * planck(nu, np.array(T)) * (nu / nu0) ** beta
    print(f"sed shape: {sed.shape}")

    sed[sed == np.inf] = (
        0  # forcing to ignore the 0 frequency which is calculated to be inf
    )
    return sed
