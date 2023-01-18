import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from fits_to_h5 import quat_to_sky_coords

fname = "/mn/stornext/d16/cmbco/ola/wmap/tods/uncalibrated/wmap_tod_20060252355_20060262355_uncalibrated_v5.fits"

data = fits.open(fname)
quat = data[1].data["QUATERN"]
# Return pointing and pol angles for each DA
gal_A, gal_B, pol_A, pol_B = quat_to_sky_coords(quat, lonlat=True, center=True)

gal_A = gal_A[0]
gal_B = gal_B[0]


gal1 = np.hstack((gal_A, gal_B))
np.savetxt("pyth_quot.txt", gal1)
print(gal1.shape)


gal2 = np.loadtxt("wmap_routines/pro/idl_quat_cent.txt")
gal2[:, 0] = (gal2[:, 0] + 360) % 360
gal2[:, 2] = (gal2[:, 2] + 360) % 360

t = np.arange(len(gal1)) * 1.536

# plt.plot(t, gal1[:len(t),0] - gal2[:len(t),0], '.', ms=1, label='A')
plt.plot(t, gal1[: len(t), 2] - gal2[: len(t), 2], ".", ms=1, label="B")
plt.xlabel("t (ms)")
plt.ylabel(r"$\Delta(\mathrm{lon})$ (deg)")
plt.ylim([-6e-3, 6e-3])
plt.savefig("dlon.png", bbox_inches="tight")

plt.figure()
# plt.plot(gal1[:len(t),1] - gal2[:len(t),1], '.', ms=1, label='A')
plt.plot(gal1[: len(t), 3] - gal2[: len(t), 3], ".", ms=1, label="B")
plt.xlabel("t (ms)")
plt.ylabel(r"$\Delta(\mathrm{lat})$ (deg)")
plt.ylim([-6e-3, 6e-3])
plt.savefig("dlat.png", bbox_inches="tight")

plt.figure()


plt.plot(t, gal1[: len(t), 0], ".", ms=1, label="Python")
plt.plot(t, gal2[: len(t), 0], ".", ms=1, label="IDL")
plt.xlabel("t (ms)")
plt.ylabel(r"$\mathrm{lon}$ (deg)")
plt.legend(loc="best")
plt.savefig("lon.png", bbox_inches="tight")


plt.figure()
plt.plot(t, gal1[: len(t), 1], ".", ms=1, label="Python")
plt.plot(t, gal2[: len(t), 1], ".", ms=1, label="IDL")
plt.xlabel("t (ms)")
plt.ylabel(r"$\mathrm{lon}$ (deg)")
plt.legend(loc="best")
plt.savefig("lat.png", bbox_inches="tight")


plt.show()
