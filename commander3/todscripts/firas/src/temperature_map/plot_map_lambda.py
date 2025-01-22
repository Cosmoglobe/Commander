import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

modes = {"ss": 0, "lf": 3}
data = np.load("../../output/data/lambda_sky.npz", allow_pickle=True)
sky = {}
for mode in modes.keys():
    sky[mode] = data[mode][:, 1:]
    # print(sky[mode])

for mode in modes.keys():
    plt.imshow(
        np.abs(sky[mode]).T,
        aspect="auto",
        extent=[0, len(sky[mode]), 0, len(sky[mode][0])],
        vmin=0,
        vmax=500,
    )
    plt.savefig(f"../../output/plots/sky_over_time_lambda_{mode}.png")
    plt.clf()


f0 = {}
df = {}
f0["ss"] = 68.020812  # GHz
df["ss"] = 13.604162  # GHz (both from the fits file)
f0["lf"] = 23.807283  # GHz
df["lf"] = 3.4010405  # GHz

f = {}
for mode in modes.keys():
    f[mode] = f0[mode] + np.arange(0, len(sky[mode][0])) * df[mode]

NSIDE = 32
npix = hp.nside2npix(NSIDE)

# get the position on the sky
data_path = "/mn/stornext/d16/cmbco/ola/firas/coadded_interferograms/"
pix_gal = {}
for mode in modes.keys():
    data = fits.open(
        data_path + f"FIRAS_COADDED_SKY_INTERFEROGRAMS_LL{mode.upper()}.FITS"
    )[1].data

    gal_lon = data["GAL_LON"]
    gal_lat = data["GAL_LAT"]
    pix_gal[mode] = hp.ang2pix(NSIDE, gal_lon, gal_lat, lonlat=True).astype(int)

hpxmap = {}
data_density = {}
for mode in modes.keys():
    hpxmap[mode] = np.zeros((npix, len(f[mode])))
    data_density[mode] = np.zeros(npix)
    for i in range(len(pix_gal[mode])):
        hpxmap[mode][pix_gal[mode][i]] += np.abs(sky[mode][i])
        data_density[mode][pix_gal[mode][i]] += 1

m = {}
for mode in modes.keys():
    m[mode] = np.zeros((npix, len(sky[mode][0])))
    mask = data_density[mode] == 0
    m[mode][~mask] = hpxmap[mode][~mask] / data_density[mode][~mask, np.newaxis]
    m[mode][mask] = hp.UNSEEN

for mode in modes.keys():
    for freq in range(len(f[mode])):
        hp.mollview(
            m[mode][:, freq],
            title=f"{f[mode][freq]:.2f} GHz",
            unit="MJy/sr",
            # norm="hist",
            min=0,
            max=400,
        )
        plt.savefig(
            f"../../output/maps/sky_map/lambda_{mode}/{int(f[mode][freq]):04d}.png"
        )
        plt.close()
