"""
Script to compare the dust simulation with the original FIRAS maps.
"""

import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits

path = "/mn/stornext/d16/cmbco/ola/firas/healpix_maps/"
filename = "FIRAS_map_0544GHz_lowf.fits"

# plot original firas map
firas_map = hp.read_map(path + filename)
hp.mollview(firas_map, title="FIRAS map", unit="MJy/sr", min=0, max=200)
hp.graticule()
plt.savefig("tests/compare/original_map.png")

# plot simulated dust map
dust_map = fits.getdata("/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/src/mapmaker/tests/dust_map_downgraded.fits")
hp.mollview(dust_map, title="Dust map downgraded", unit="MJy/sr", min=0, max=200)
hp.graticule()
plt.savefig("tests/compare/dust_map.png")

# plot simulated 545 map
simulated_map = hp.read_map("tests/dust_maps/0546.fits")
hp.mollview(simulated_map, title="Simulated 545 GHz map", unit="MJy/sr", min=0, max=200)
hp.graticule()
plt.savefig("tests/compare/simulated_map.png")

# plot difference map
difference_map = firas_map - simulated_map
hp.mollview(difference_map, title="Original - Simulated", unit="MJy/sr", min=-200, max=200)
hp.graticule()
plt.savefig("tests/compare/difference_map.png")

# plot difference between input simulated map and output of invert mapmaker
m_invert = hp.read_map("tests/m_invert/0546.fits")
hp.mollview(m_invert, title="Invert map", unit="MJy/sr", min=0, max=200)
hp.graticule()
plt.savefig("tests/compare/invert_map.png")
difference_map = simulated_map - m_invert
hp.mollview(difference_map, title="Simulated - Invert", unit="MJy/sr", min=-200, max=200)
hp.graticule()
plt.savefig("tests/compare/difference_map_invert.png")

# plot difference between input simulated map and output of cg solver
cg_map = hp.read_map("tests/m_cg_per_tod/0546.fits")
hp.mollview(cg_map, title="CG map", unit="MJy/sr", min=0, max=200)
hp.graticule()
plt.savefig("tests/compare/cg_map.png")
difference_map = simulated_map - cg_map
hp.mollview(difference_map, title="Simulated - CG", unit="MJy/sr", min=-200, max=200)
hp.graticule()
plt.savefig("tests/compare/difference_map_cg.png")
