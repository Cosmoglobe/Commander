"""
Downgrades the Planck map to a lower resolution.
"""

import astropy.units as u
import globals_mapmaker as g
import healpy as hp
import numpy as np
from astropy.io import fits

dust_map = "/mn/stornext/d16/cmbco/ola/planck_products/commander/COM_CompMap_ThermalDust-commander_n2048_R2.00.fits"

# print(repr(fits.open(dust_map)[1].header))

dust_map = fits.open(dust_map)[1].data["I_ML_FULL"]

hp.mollview(dust_map, title="Dust map", unit="$\\mu K_{RJ}$", nest=True)
# plt.savefig("tests/dust_map.png")
# plt.close()

# smooth map to 7 degrees
dust_map_smoothed = hp.smoothing(dust_map, fwhm=7 * np.pi / 180, nest=True)
hp.mollview(
    dust_map_smoothed, title="Smoothed dust map", unit="$\\mu K_{RJ}$", nest=True
)
# plt.savefig("tests/dust_map_smoothed.png")
# plt.close()

# downgrade to nside = 32
dust_map_downgraded = hp.ud_grade(
    dust_map_smoothed, nside_out=g.NSIDE, order_in="NESTED", order_out="RING"
)
# hp.mollview(dust_map_downgraded, title="Downgraded dust map", unit="$\\mu K_{RJ}$")
# plt.savefig("tests/dust_map_downgraded.png")
# plt.close()

nu0_dust = 545 # Planck

# let's change the units from uK_RJ to MJy/sr
dust_map_downgraded_mjy = (dust_map_downgraded * u.uK).to(
    u.MJy / u.sr,
    equivalencies=u.brightness_temperature(nu0_dust * u.GHz),
)
# dust_map_downgraded_mjy = (dust_map_downgraded * u.uK / (const.c**2) * 2 * const.k_B * nu0_dust**2).to(u.MJy/u.sr, equivalencies=u.brightness_temperature(nu0_dust))
# print(I_nu.to(u.MJy / u.sr, equivalencies=u.brightness_temperature(nu0_dust)))

# save downgraded map in mjy/sr in order to be able to open it again in python
fits.writeto(
    "test_output/dust_map_downgraded.fits", dust_map_downgraded_mjy.value, overwrite=True
)
