import healpy as hp
import numpy as np
from astropy.io import fits

# planck maps paths 
planck_maps_path = "/mn/stornext/d16/cmbco/ola/planck_products/dr3/"
maps = ["LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits", "LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits", "LFI_SkyMap_070-BPassCorrected-tpl_1024_R3.00_full.fits", "HFI_SkyMap_100_2048_R3.01_full.fits", "HFI_SkyMap_143_2048_R3.01_full.fits", "HFI_SkyMap_217_2048_R3.01_full.fits", "HFI_SkyMap_353-psb_2048_R3.01_full.fits", "HFI_SkyMap_545_2048_R3.01_full.fits", "HFI_SkyMap_857_2048_R3.01_full.fits"]
bands = ["030", "044", "070", "100", "143", "217", "353", "545", "857"]

output_path = "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/output/fits_files/"

# smooth planck maps to 7 degrees
for i, map in enumerate(maps):
    map_path = planck_maps_path + map
    map_data = fits.open(map_path)[1].data["I_STOKES"]
    map_smoothed = hp.smoothing(map_data, fwhm=7 * np.pi / 180, nest=True)
    map_smoothed_udgrade = hp.ud_grade(map_smoothed, nside_out=32, order_in="NESTED", order_out="RING")
    map_smoothed_path = output_path + "planck_" + bands[i] + "_smoothed.fits"
    fits.writeto(map_smoothed_path, map_smoothed_udgrade, overwrite=True)
