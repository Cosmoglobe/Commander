import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

from astropy import units as u

from glob import glob

fnames = glob("data/wmap_sidelobe*.fits")
fnames.sort()
for i in range(len(fnames)):
    d = hp.read_map(fnames[i])
    slAB = hp.reorder(d, n2r=True)
    hp.gnomview(slAB, rot=(90, -19.5, 0), reso=5, min=-0.5, max=0.5, cmap="coolwarm")
plt.close("all")
plt.show()
d = hp.read_map("data/wmap_sidelobe_map_K1_9yr_v5.fits")
# d = hp.read_map('data/wmap_sidelobe_map_W4_9yr_v5.fits')
slAB = hp.reorder(d, n2r=True)

# The normalization is sum(abs(slAB)) = 2*Npix, which would be true if they hadn't cut out the main beam.
# We need to get the same normalization for the beam.

slA = 0 * slAB
slA[slAB > 0] = slAB[slAB > 0]

# This normalization should be sum(slA) = Npix, which would again be true without the hole.
# My supposition is that if we scaled by Npix_new/Npix_old, that would be the correct normalization.
# Meanwhile, I think it's fair to assume that the main beam would correspond to the difference, i.e.
# sum(beam) = Npix_new - sum(slA)

npix_new = hp.nside2npix(4096)
npix_old = hp.nside2npix(128)
# slA *= npix_new/npix_old

beam = hp.read_map("freq0_hornA.fits")
# beam = hp.read_map('freq9_hornA.fits')

beam = beam * (npix_new - sum(slA)) / sum(beam)

# beam = beam/npix_new
# slA  = slA/npix_new

# from program_pars

vec = np.array([0.00980808738477, 0.93934799994236, -0.34282522723123])
lon, lat = hp.vec2ang(vec, lonlat=True)
rot1 = hp.Rotator(rot=[(0, -90 + lat)])
rot2 = hp.Rotator(rot=[(lon, 0)])
# [89.40177486] [-20.04909596]


rot0 = hp.Rotator(rot=[(0, 0, 0)])
# rot1 = hp.Rotator(rot=[(0, 90+19.5)])
rot1 = hp.Rotator(rot=[(0, 90 + 24.5)])
# rot2 = hp.Rotator(rot=[(-90, 0)])
rot2 = hp.Rotator(rot=[(-90.5, 0)])


beam = hp.ud_grade(beam, 1024)

m_rotated = rot0.rotate_map_pixel(beam)
m_rotated = rot1.rotate_map_pixel(m_rotated)
m_rotated = rot2.rotate_map_pixel(m_rotated)


# bl = hp.anafast(beam)
# sl = hp.anafast(slA)
# plt.semilogy(bl, label='Beam')
# plt.semilogy(sl, label='Sidelobe')
# plt.xscale('log')
# plt.legend(loc='best')
# plt.xlabel(r'$\ell$')
##plt.show()
#
# plt.figure()
# theta = np.linspace(0.01*np.pi/180, 2*np.pi,10000)
# plt.loglog(theta*180/np.pi, abs(hp.bl2beam(bl, theta)))
# plt.plot(theta*180/np.pi, abs(hp.bl2beam(sl, theta)))
#
# plt.ylim([3e-20, 3e-9])
# plt.xlim([0.1, 200])


# I know this isn't totally right, but if I just do a single rotation in
# mollview, it seems to look okay...
# b2 = np.log10(beam)
# b2[~np.isfinite(b2)] = hp.UNSEEN
# hp.mollview(b2, rot=(0,19.5+90,0))
# hp.mollview(slA, rot=(90,0,0))
# hp.mollview(np.log10(m_rotated), rot=(90,0,0))

slA = hp.ud_grade(slA, hp.npix2nside(len(m_rotated)))
hp.mollview(slA + m_rotated, min=-0.5, max=0.5, cmap="coolwarm")
hp.gnomview(
    slA + m_rotated, rot=(90, -19.5, 0), reso=10, min=-0.5, max=0.5, cmap="coolwarm"
)
hp.graticule()
hp.gnomview(m_rotated, rot=(90, -19.5, 0), reso=10, min=-0.5, max=0.5, cmap="coolwarm")
hp.graticule()

plt.show()
