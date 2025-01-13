import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from my_utils import planck, residuals
from scipy.optimize import minimize
from utils.config import gen_nyquistl

T_CMB = 2.72548  # Fixsen 2009

sky = np.load("../../output/data/sky.npy")
sky = np.abs(sky)
print(f"sky shape: {sky.shape}")
data = h5py.File(
    "/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/data/sky_v4.1.h5",
    "r",
)
mask = fits.open("BP_CMB_I_analysis_mask_n1024_v2.fits")
mask = mask[1].data.astype(int)

pix_gal = np.array(data["df_data/pix_gal"]).astype(int)
mtm_length = np.array(data["df_data/mtm_length"][()])
mtm_speed = np.array(data["df_data/mtm_speed"][()])

short_filter = mtm_length == 0
slow_filter = mtm_speed == 0
pix_gal = pix_gal[short_filter & slow_filter]

# to not use ICAL higher than 3 temps
a_ical = np.array(data["df_data/a_ical"][()])
b_ical = np.array(data["df_data/b_ical"][()])
ical = (a_ical + b_ical) / 2
ical = ical[short_filter & slow_filter]
ical_lower_3 = ical < 3
pix_gal = pix_gal[ical_lower_3]
sky = sky[ical_lower_3]

# frequency mapping
fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

scan_mode = 0  # SS
channel = 3  # LL
frec = 4 * (channel % 2) + scan_mode

f_icm = np.arange(len(sky[0])) * (fnyq["icm"][frec] / 320)
c = 3e8 * 1e2  # cm/s
f_ghz = (
    f_icm * c * 1e-9 + 55
)  # this might not be right but it is what matches the initial frequencies of the firas movie

NSIDE = 32
npix = hp.nside2npix(NSIDE)

bb_curve = np.zeros(len(f_ghz))

print("Calculating BB curve")

for freq in range(len(f_ghz)):
    hpxmap = np.zeros(npix)
    data_density = np.zeros(npix)

    for i in range(len(pix_gal)):
        hpxmap[pix_gal[i]] += np.abs(sky[i][freq])
        data_density[pix_gal[i]] += 1

    m = np.zeros(npix)
    mask2 = data_density == 0
    m[~mask2] = hpxmap[~mask2] / data_density[~mask2]

    # mask the map with the points of no data
    m[mask2] = hp.UNSEEN

    # mask the galaxy
    m[mask] = hp.UNSEEN

    # average all points at this frequency to get each frequency point in the bb curve
    bb_curve[freq] = np.mean(m[m != hp.UNSEEN])

print(f"Fitting BB curve: {bb_curve}")

t0 = np.array(2.726)
fit = minimize(residuals, t0, args=(f_ghz[1:], bb_curve[1:]))

print(f"Plotting BB curve: {fit.x[0]}")

plt.plot(f_ghz[1:], bb_curve[1:], label="Data")
plt.plot(f_ghz[1:], planck(f_ghz, fit.x[0])[1:], label="Fit")
plt.plot(f_ghz[1:], planck(f_ghz, t0)[1:], label="Original")
plt.xlabel("Frequency [GHz]")
plt.ylabel("Brightness [MJy/sr]")
plt.title("BB curve")
plt.legend()
plt.savefig("../../output/plots/bb_curve.png")
