import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from my_utils import planck

T_CMB = 2.72548  # Fixsen 2009

data_path = "/mn/stornext/d16/cmbco/ola/firas/spectra_uncombined/"
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}
modes = {"ss": 0, "lf": 3}

nu0 = {"ss": 68.020812, "lf": 23.807283}
dnu = {"ss": 13.604162, "lf": 3.4010405}
nf = {"lh_ss": 210, "ll_lf": 182, "ll_ss": 43, "rh_ss": 210, "rl_lf": 182, "rl_ss": 43}

f_ghz = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            f_ghz[f"{channel}_{mode}"] = np.linspace(
                nu0[mode],
                nu0[mode] + dnu[mode] * (nf[f"{channel}_{mode}"] - 1),
                nf[f"{channel}_{mode}"],
            )

NSIDE = 32
npix = hp.nside2npix(NSIDE)

for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            data = fits.open(
                data_path
                + f"FIRAS_CALIBRATED_SKY_SPECTRA_{channel.upper()}{mode.upper()}.FITS"
            )

            # print(repr(data[0].header))

            sky = (
                data[1].data.field("REAL_SPE") + 1j * data[1].data.field("IMAG_SPE")
            )[:, : len(f_ghz[f"{channel}_{mode}"])]

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
            plt.title(f"Sky over time for channel {channel} and mode {mode}")
            plt.savefig(
                f"../../output/plots/sky_over_time_lambda_spectra_{channel}_{mode}.png"
            )

            pix_gal = hp.ang2pix(NSIDE, gal_lon, gal_lat, lonlat=True).astype(int)

            hpxmap = np.zeros((npix, len(f_ghz[f"{channel}_{mode}"])))
            data_density = np.zeros(npix)
            for i in range(len(pix_gal)):
                hpxmap[pix_gal[i]] += np.abs(sky[i])
                data_density[pix_gal[i]] += 1

            m = np.zeros((npix, len(f_ghz[f"{channel}_{mode}"])))
            mask = data_density == 0
            monopole = planck(f_ghz[f"{channel}_{mode}"], np.array(T_CMB))
            m[~mask] = hpxmap[~mask] / data_density[~mask][:, np.newaxis] - monopole
            m[mask] = hp.UNSEEN

            for freq in range(len(f_ghz[f"{channel}_{mode}"])):
                hp.mollview(
                    m[:, freq],
                    title=f"Sky map at frequency {f_ghz[f"{channel}_{mode}"][freq]:.2f} GHz for mode {mode}",
                    # min=0,
                    # max=200,
                    norm="hist",
                )
                plt.savefig(
                    f"../../output/maps/lambda_spectra/{channel}{mode}/{int(f_ghz[f"{channel}_{mode}"][freq]):04d}.png"
                )
                plt.close()
