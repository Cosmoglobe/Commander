import astropy.units as u
import cosmoglobe as cg
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

import globals as g
import utils.my_utils as utils

mask_gal = hp.read_map("/mn/stornext/d16/cmbco/ola/masks/HI_mask_4e20_n1024.fits")
mask_alm = hp.sphtfunc.map2alm(mask_gal, pol=False)
mask_gal = hp.alm2map(mask_alm, g.NSIDE, pol=False)
mask_gal = np.where(mask_gal < 0.5, 1, 0)

for channel in g.CHANNELS_PLOT:
    for mode in g.MODES_PLOT:
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            data = np.load(
                f"{g.PROCESSED_DATA_PATH}sky_{channel}_{mode}.npz", allow_pickle=True
            )

            sky = data["sky"]
            ical = data["ical"]

            plt.hist(ical, bins=100, range=(2.7, 2.8))
            plt.savefig("visualizations/split_maps/ical_histogram_" + channel + "_" + mode + ".png")

            # split ICAL temps into 2
            split = 2.76 # K
            sky_low_ical = sky[ical <= split]
            sky_high_ical = sky[ical > split]

            print(f"Data volume in low ical map ({channel}, {mode}): {len(sky_low_ical)/len(sky)*100:.2f}%")
            print(f"Data volume in high ical map ({channel}, {mode}): {len(sky_high_ical)/len(sky)*100:.2f}%")

            gal_lat = data["gal_lat"]
            gal_lon = data["gal_lon"]
            pix_gal = hp.ang2pix(g.NSIDE, gal_lon, gal_lat, lonlat=True).astype(int)

            f_ghz = utils.generate_frequencies(channel, mode)

            # plot maps for high and low ICAL temps separate
            # first the low ICAL
            pix_gal_low = pix_gal[ical <= split]

            cum_maps = np.zeros((g.NPIX, len(f_ghz)))
            data_density = np.zeros(g.NPIX)

            for i, pix in enumerate(pix_gal_low):
                data_density[pix] += 1
                cum_maps[pix] += np.abs(sky_low_ical[i])

            mask = data_density == 0
            monopole = utils.planck(f_ghz, g.T_CMB)

            maps_low = np.zeros((g.NPIX, len(f_ghz)))
            maps_low[~mask] = cum_maps[~mask] / data_density[~mask, None] - monopole
            maps_low[mask] = np.nan

            for nu_i, nu in enumerate(f_ghz):
                cg.plot(maps_low[:, nu_i], title=f"Low ICAL @ {int(nu)} GHz - {channel} {mode}",
                        unit="MJy/sr", min=0, max=200, comp="freqmap", freq=nu * u.GHz,
                        cmap="planck")
                plt.savefig("visualizations/split_maps/low_ical/" + channel + "_" + mode + f"/{int(nu):03d}.png")
                plt.close()

            # now the high ICAL
            pix_gal_high = pix_gal[ical > split]
            cum_maps = np.zeros((g.NPIX, len(f_ghz)))
            data_density = np.zeros(g.NPIX)
            for i, pix in enumerate(pix_gal_high):
                data_density[pix] += 1
                cum_maps[pix] += np.abs(sky_high_ical[i])

            mask = data_density == 0
            
            maps_high = np.zeros((g.NPIX, len(f_ghz)))
            maps_high[~mask] = cum_maps[~mask] / data_density[~mask, None] - monopole
            maps_high[mask] = np.nan
            for nu_i, nu in enumerate(f_ghz):
                cg.plot(maps_high[:, nu_i], title=f"High ICAL @ {int(nu)} GHz - {channel} {mode}",
                        unit="MJy/sr", min=0, max=200, comp="freqmap", freq=nu * u.GHz,
                        cmap="planck")
                plt.savefig("visualizations/split_maps/high_ical/" + channel + "_" + mode + f"/{int(nu):03d}.png")
                plt.close()


            common_mask = (maps_low == np.nan) | (maps_high == np.nan)

            # difference map
            maps_diff = maps_high - maps_low
            maps_diff[common_mask] = np.nan
            for nu_i, nu in enumerate(f_ghz):
                cg.plot(maps_diff[:, nu_i],
                        title=f"High - Low ICAL @ {int(nu)} GHz - {channel} {mode}",
                        unit="MJy/sr", min=-10, max=10, comp="residual", freq=nu * u.GHz,
                        cmap="RdBu_r")
                plt.savefig("visualizations/split_maps/diff/" + channel + "_" + mode + f"/{int(nu):03d}.png")

            # mask the galaxy and plot teh spectrum
            # plot mask_gal
            hp.mollview(mask_gal, title="Galaxy Mask", unit="", min=0, max=1, cmap="gray")
            plt.savefig("debug.png")
            plt.close()

            diff_masked = maps_diff.copy()
            diff_masked[mask_gal == 1] = np.nan
            hp.mollview(diff_masked[:, 3],
                           unit="MJy/sr", min=-10, max=10, cmap="RdBu_r")
            plt.savefig("debug.png")
            plt.close()

            # average for each frequency
            spectrum = np.nanmean(diff_masked, axis=0)
            plt.figure()
            plt.plot(f_ghz, spectrum)
            plt.title(f"High - Low ICAL Spectrum - {channel} {mode}")
            plt.xlabel("Frequency (GHz)")
            plt.ylabel("Difference (MJy/sr)")
            plt.savefig("visualizations/split_maps/spectrum_" + channel + "_" + mode + ".png")
            plt.close()
