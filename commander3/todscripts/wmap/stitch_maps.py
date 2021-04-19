import numpy as np
import matplotlib.pyplot as plt
import healpy as hp


d = hp.read_map('data/wmap_sidelobe_map_K1_9yr_v5.fits')
d = hp.read_map('data/wmap_sidelobe_map_W4_9yr_v5.fits')
slAB = hp.reorder(d, n2r=True)

# The normalization is sum(abs(slAB)) = 2*Npix, which would be true if they hadn't cut out the main beam.
# We need to get the same normalization for the beam.

slA = 0*slAB
slA[slAB > 0] = slAB[slAB > 0]

# This normalization should be sum(slA) = Npix, which would again be true without the hole.
# My supposition is that if we scaled by Npix_new/Npix_old, that would be the correct normalization. 
# Meanwhile, I think it's fair to assume that the main beam would correspond to the difference, i.e.
# sum(beam) = Npix_new - sum(slA)

npix_new = hp.nside2npix(1024)
npix_old = hp.nside2npix(128)
slA *= npix_new/npix_old

beam = hp.read_map('freq0_hornA.fits')
beam = hp.read_map('freq9_hornA.fits')

beam = beam*(npix_new - sum(slA))/sum(beam)

beam = beam/npix_new
slA  = slA/npix_new


bl = hp.anafast(beam)
sl = hp.anafast(slA)
plt.semilogy(bl, label='Beam')
plt.semilogy(sl, label='Sidelobe')
plt.xscale('log')
plt.legend(loc='best')
plt.xlabel(r'$\ell$')
#plt.show()

plt.figure()
theta = np.linspace(0.01*np.pi/180, 2*np.pi,10000)
plt.loglog(theta*180/np.pi, abs(hp.bl2beam(bl, theta)))
plt.plot(theta*180/np.pi, abs(hp.bl2beam(sl, theta)))

plt.ylim([3e-20, 3e-9])
plt.xlim([0.1, 200])



# I know this isn't totally right, but if I just do a single rotation in
# mollview, it seems to look okay...
hp.mollview(beam, rot=(0,0,90+22.19))
hp.mollview(slA, rot=(0,0,-90-22.19))

beam = hp.ud_grade(beam, 128)

hp.mollview(beam + slA, rot=(0,0,-90-22.19))
plt.show()
