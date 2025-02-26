from pixell import enmap, enplot, reproject, utils, curvedsky
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import cosmoglobe
import astropy.units as u

OLA = "/mn/stornext/d16/cmbco/ola/actpol"

HFI = "npipe_100_2deg.fits"
HFI = "npipe_100_0.5deg.fits"

# try:
#    HFI = hp.read_map(HFI, field=(0,1,2))
# except FileNotFoundError:
#    data = hp.read_map('/mn/stornext/d16/cmbco/ola/npipe/freqmaps/npipe6v20_100_map_K.fits',
#        field=(0,1,2))
#    data = hp.smoothing(data, fwhm=0.5*np.pi/180)*1e3
#    hp.write_map(HFI, data)
#    HFI = data

nside = 512

model = cosmoglobe.sky_model(nside=nside)
simulated_emission = model(90 * u.GHz, fwhm=60 * u.arcmin)


dec_from, dec_to = np.deg2rad([-20, 10])
ra_from, ra_to = np.deg2rad([-90, 90])
box = [[dec_from, ra_from], [dec_to, ra_to]]

# hires region

# dec_from, dec_to = np.deg2rad([4, 19])
# ra_from, ra_to = np.deg2rad([-175, 124])
# box = [[dec_from,ra_from],[dec_to,ra_to]]

ipix = np.arange(12 * nside**2)
theta, phi = hp.pix2ang(nside, ipix)

lon = np.pi / 2 - theta
lat = phi
lat[lat > np.pi] = np.pi - lat[lat > np.pi]

print(lon.min(), lon.max(), lat.min(), lat.max())

# imap = enmap.read_map(f'{OLA}/act_dr5.01_s08s18_AA_f090_night_map.fits', box=box)
# imap = enmap.read_map('act_planck_dr5.01_s08s18_AA_f090_night_map.fits', box=box)
imap = enmap.read_map(f"{OLA}/act_dr5.01_s08s18_AA_f090_night_map.fits")
# imap = enmap.read_map('act_planck_dr5.01_s08s18_AA_f090_night_map.fits')


alm = curvedsky.map2alm(imap, lmax=3 * nside - 1)
rot_eq2gal = hp.Rotator(coord="CG")
rot_gal2eq = hp.Rotator(coord="GC")
TQU = hp.alm2map(alm.astype(np.complex128), nside)


inds = abs(TQU) < 1e-1

# inds = (lon > box[0][0]) & (lon < box[1][0])
# inds = np.array([inds, inds, inds])

# TQU[inds] = hp.UNSEEN

TQU_GC = rot_eq2gal.rotate_map_alms(TQU)

# cosmoglobe.plot(TQU_GC, min=-10, max=10, sig=1, fwhm=0.1*u.deg)
cosmoglobe.plot(TQU_GC, sig=0, fwhm=0.1 * u.deg)
plt.show()

# emission_GC = rot_gal2eq.rotate_map_alms(simulated_emission)
#
# cosmoglobe.plot(TQU, min=-10, max=10, sig=1, llabel='ACT', cbar=False)
# cosmoglobe.plot(emission_GC,
#    sig=1, min=-10, max=10, llabel='Total', cbar=False)
#
#
# emission_GC[~inds] = TQU[~inds]
#
# cosmoglobe.plot(emission_GC,
#    sig=1, min=-10, max=10, xsize=8000, cbar=False, llabel='Combined')
# plt.savefig('act_90_v_CG.png', bbox_inches='tight', dpi=300)
# plt.show()
