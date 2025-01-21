import numpy as np

def dSdT(nu):
    h = 6.62607004e-34  # m^2 kg / s
    c = 299792458.0     # m / s
    T = 2.72548         # K
    k = 1.380649e-23    # J / K
    
    x = h * nu / (k * T)
    dIdT = 2 * x ** 4 * k ** 3 * T ** 2 / (h * c) ** 2 / (4 * np.sinh(x / 2) ** 2) * 1e26
    return dIdT


def jy2kcmb(flux, nu, area):
    # flux: flux in Jy 
    # nu  : frequency 
    # area: beam area (i.e. 2 * pi * sigma ** 2)
    return flux / area / dSdT(nu)


