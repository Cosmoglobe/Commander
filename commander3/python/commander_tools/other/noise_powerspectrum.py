import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

# Maps to generate power spectrum from
map1_file = "chains_LB_365d_40GHz_1of_noise/tod_040_res_c0001_k000002.fits"
map2_file = "chains_LB_365d_40GHz_sample_1of_noise_parameters_v3/tod_040_ncorr_c0001_k000006.fits" 
map3_file = "chains_LB_365d_40GHz_sample_1of_noise_parameters_v3/tod_040_res_c0001_k000006.fits"
map4_file = "chains_LB_365d_40GHz_sample_1of_noise_parameters_v3/res_040_c0001_k000006.fits"

lmax = 350
beam_40GHz = 69.3 #beam in maps, in arcmin

# Import maps, both T, Q, U
map1 = hp.fitsfunc.read_map(map1_file, field=(0,1,2))
map2 = hp.fitsfunc.read_map(map2_file, field=(0,1,2))
map3 = hp.fitsfunc.read_map(map3_file, field=(0,1,2))
map4 = hp.fitsfunc.read_map(map4_file, field=(0,1,2))

# vector with l's and gaussian beam
ell = np.arange(lmax+1)
beam = hp.sphtfunc.gauss_beam((beam_40GHz/60)*np.pi/180, lmax=lmax, pol=True)
beam = beam.T

# genereting Cl', dividing with l*(l+1)/(2pi) and deconvlvind beamdeconvolving
spectrum1_cl = hp.sphtfunc.anafast(map1, map1, lmax=lmax)
spectrum1_dl = ell*(ell+1)*spectrum1_cl/2/np.pi  #  divide with % of sky if masked
spectrum1_dl[2]= spectrum1_dl[2]/(beam[2]**2)
spectrum1_dl[1]= spectrum1_dl[1]/(beam[1]**2)

spectrum2_cl = hp.sphtfunc.anafast(map2, map2, lmax=lmax)
spectrum2_dl = ell*(ell+1)*spectrum2_cl/2/np.pi
spectrum2_dl[2]= spectrum2_dl[2]/(beam[2]**2)
spectrum2_dl[1]= spectrum2_dl[1]/(beam[1]**2)

spectrum3_cl = hp.sphtfunc.anafast(map3, map3, lmax=lmax)
spectrum3_dl = ell*(ell+1)*spectrum3_cl/2/np.pi
spectrum3_dl[2]= spectrum3_dl[2]/(beam[2]**2)
spectrum3_dl[1]= spectrum3_dl[1]/(beam[1]**2)

spectrum4_cl = hp.sphtfunc.anafast(map4, map4, lmax=lmax)
spectrum4_dl = ell*(ell+1)*spectrum4_cl/2/np.pi
spectrum4_dl[2]= spectrum4_dl[2]/(beam[2]**2)
spectrum4_dl[1]= spectrum4_dl[1]/(beam[1]**2)

# Plotting Cl
plt.loglog(ell[2:],spectrum1_cl[1][2:], label="tod_res, no fitting")
plt.loglog(ell[2:],spectrum2_cl[1][2:], label="ncorr")
plt.loglog(ell[2:],spectrum3_cl[1][2:], label="tod_res")
plt.loglog(ell[2:],spectrum4_cl[1][2:], label="res")

plt.legend(loc="lower right")
plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_\ell$')

plt.savefig("LB_powerspectrum_noise.png")
plt.show()

# run with %run filename.py
