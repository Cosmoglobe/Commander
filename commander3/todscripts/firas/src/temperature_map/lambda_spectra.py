import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

data_path = "/mn/stornext/d16/cmbco/ola/firas/spectra_uncombined/"
modes = {"ss": 0, "lf": 3}

f0 = {}
df = {}
f0["ss"] = 68.020812  # GHz
df["ss"] = 13.604162  # GHz
f0["lf"] = 23.807283  # GHz
df["lf"] = 3.4010405  # GHz

data_points = {}
data_points["ss"] = 43
data_points["lf"] = 182

f = {}
f["ss"] = f0["ss"] + np.arange(0, data_points["ss"]) * df["ss"]
f["lf"] = f0["lf"] + np.arange(0, data_points["lf"]) * df["lf"]

NSIDE = 32
npix = hp.nside2npix(NSIDE)

for mode in modes:
    data = fits.open(data_path + f"FIRAS_CALIBRATED_SKY_SPECTRA_LL{mode.upper()}.FITS")

    # print(repr(data[0].header))

    sky = (data[1].data.field("REAL_SPE") + 1j * data[1].data.field("IMAG_SPE"))[
        :, : data_points[mode]
    ]

    gal_lon = data[1].data.field("GAL_LON")
    gal_lat = data[1].data.field("GAL_LAT")
    data.close()

    plt.imshow(
        np.abs(sky).T,
        aspect="auto",
        extent=[0, len(sky), 0, len(sky[0])],
        vmin=0,
        vmax=500,
    )
    plt.title(f"Sky over time for mode {mode}")
    plt.savefig(f"../../output/plots/sky_over_time_lambda_spectra_{mode}.png")

    pix_gal = hp.ang2pix(NSIDE, gal_lon, gal_lat, lonlat=True).astype(int)

    hpxmap = np.zeros((npix, data_points[mode]))
    data_density = np.zeros(npix)
    for i in range(len(pix_gal)):
        hpxmap[pix_gal[i]] += np.abs(sky[i])
        data_density[pix_gal[i]] += 1

    m = np.zeros((npix, data_points[mode]))
    mask = data_density == 0
    m[~mask] = hpxmap[~mask] / data_density[~mask][:, np.newaxis]
    m[mask] = hp.UNSEEN

    for freq in range(data_points[mode]):
        hp.mollview(
            m[:, freq],
            title=f"Sky map at frequency {f[mode][freq]:.2f} GHz for mode {mode}",
            min=0,
            max=400,
        )
        plt.savefig(
            f"../../output/maps/sky_map/lambda_spectra_{mode}/{int(f[mode][freq]):04d}.png"
        )
        plt.close()
