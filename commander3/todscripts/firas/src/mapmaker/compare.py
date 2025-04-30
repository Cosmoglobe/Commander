"""
Script to compare the dust simulation with the original FIRAS maps.
"""

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
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
hp.mollview(difference_map, title="Simulated - Invert", unit="MJy/sr", min=-1, max=1)
hp.graticule()
plt.savefig("tests/compare/difference_map_invert.png")

# check ratio between simulated map and invert map
ratio_map = simulated_map / m_invert
print("Ratio between simulated map and invert map: ", ratio_map)
# plot ratio map
hp.mollview(ratio_map, title="Ratio map", unit="MJy/sr", norm='hist')
hp.graticule()
plt.savefig("tests/compare/ratio_map.png")

# plot difference between input simulated map and output of cg solver
cg_map = hp.read_map("tests/m_cg_per_tod/0546.fits")
hp.mollview(cg_map, title="CG map", unit="MJy/sr", min=0, max=200)
hp.graticule()
plt.savefig("tests/compare/cg_map.png")
difference_map = simulated_map - cg_map
hp.mollview(difference_map, title="Simulated - CG", unit="MJy/sr", min=-1, max=1)
hp.graticule()
plt.savefig("tests/compare/difference_map_cg.png")
# check ratio between simulated map and cg map
ratio_map = simulated_map / cg_map
print("Ratio between simulated map and cg map: ", ratio_map)
# plot ratio map
hp.mollview(ratio_map, title="Ratio map", unit="MJy/sr", min=1, max=3)
hp.graticule()
plt.savefig("tests/compare/ratio_map_cg.png")

# check differences in the maps in the middle and check for weird numbers
# check for nans
print(f"Number of nans in simulated map: {np.isnan(simulated_map).sum()}")
print(f"Number of nans in firas map: {np.isnan(firas_map).sum()}")
print(f"Number of nans in invert map: {np.isnan(m_invert).sum()}")
print(f"Number of nans in cg map: {np.isnan(cg_map).sum()}")
print(f"Number of nans in difference map: {np.isnan(difference_map).sum()}")
print(f"Number of nans in ratio map: {np.isnan(ratio_map).sum()}")

# check numbers in the middle of the map
print(f"Simulated map middle: {simulated_map[3072//2]}")
print(f"Firas map middle: {firas_map[3072//2]}")
print(f"Difference map middle: {difference_map[3072//2]}")
print(f"Ratio map middle: {ratio_map[3072//2]}")
print(f"CG map middle: {cg_map[3072//2]}")
print(f"m_invert map middle: {m_invert[3072//2]}")