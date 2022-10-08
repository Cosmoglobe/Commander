import healpy as hp
import cosmoglobe
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt


w = 5

fwhm = 2*np.pi/180

d_DR4_044 = hp.read_map('/mn/stornext/d16/cmbco/ola/npipe/freqmaps/npipe6v20_044_map_K.fits', field=(0,1,2))*1e6
d_BP_044 = hp.read_map('/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2/BP_044_IQU_n0512_v2.fits', field=(0,1,2))

d_DR4_044 = hp.ud_grade(hp.smoothing(d_DR4_044, fwhm=fwhm), 512)
d_BP_044 = hp.smoothing(d_DR4_044, fwhm=fwhm)

d_Q1_orig = hp.read_map('/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Q1_v5.fits',
    field=(0,1,2))*1e3

d_Q1_CG = hp.read_map('/mn/stornext/d16/www_cmb/dwatts/v0/BP_040-WMAP_Q1_IQU_n0512_v0.fits',
    field=(0,1,2))*1e3


d_Q1_orig = hp.smoothing(d_Q1_orig, fwhm=fwhm)
d_Q1_CG = hp.smoothing(d_Q1_CG, fwhm=fwhm)

cosmoglobe.plot(d_DR4_044, sig=1, min=-100, max=100)
cosmoglobe.plot(d_BP_044, sig=1, min=-100, max=100)
cosmoglobe.plot(d_Q1_orig, sig=1, min=-100, max=100)
cosmoglobe.plot(d_Q1_CG, sig=1, min=-100, max=100)


cosmoglobe.plot(d_DR4_044 - d_Q1_orig*0.78, min=-10, max=10,
    sig=1)
cosmoglobe.plot(d_BP_044 - d_Q1_orig*0.78, min=-10, max=10,
    sig=1)
cosmoglobe.plot(d_DR4_044 - d_Q1_CG*0.78, min=-10, max=10, sig=1)
cosmoglobe.plot(d_BP_044 - d_Q1_CG*0.78, min=-10, max=10, sig=1)
plt.show()
