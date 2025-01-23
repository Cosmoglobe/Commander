import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from my_utils import planck, residuals
from scipy.optimize import minimize
from utils.config import gen_nyquistl

T_CMB = 2.72548  # Fixsen 2009
modes = {"ss": 0, "lf": 3}

mask = fits.open("BP_CMB_I_analysis_mask_n1024_v2.fits")
mask = mask[1].data.astype(int)

data = np.load("../../output/data/sky.npz")
# print(data.files)

sky = {}
pix_gal = {}
for mode in modes:
    sky[mode] = np.abs(data[f"{mode}"])
    pix_gal[mode] = data[f"pix_gal_{mode}"]

# stat_word_1 = np.array(data["df_data/stat_word_1"][()]).astype(int)
# stat_word_12 = np.array(data["df_data/stat_word_12"][()]).astype(int)
# stat_word_9 = np.array(data["df_data/stat_word_9"][()]).astype(int)
# lvdt_stat_b = np.array(data["df_data/lvdt_stat_b"][()]).astype(int)
# stat_word_13 = np.array(data["df_data/stat_word_13"][()]).astype(int)
# stat_word_16 = np.array(data["df_data/stat_word_16"][()]).astype(int)


# to not use ICAL higher than 3 temps
# a_ical = np.array(data["df_data/a_ical"][()])
# b_ical = np.array(data["df_data/b_ical"][()])
# ical = (a_ical + b_ical) / 2
# ical_lower_3 = ical < 3
# pix_gal = pix_gal[ical_lower_3]
# sky = sky[ical_lower_3]

# remove data that i flagged
# stat_word_1 = stat_word_1[short_filter & slow_filter]
# stat_word_12 = stat_word_12[short_filter & slow_filter]
# stat_word_9 = stat_word_9[short_filter & slow_filter]
# lvdt_stat_b = lvdt_stat_b[short_filter & slow_filter]
# stat_word_13 = stat_word_13[short_filter & slow_filter]
# stat_word_16 = stat_word_16[short_filter & slow_filter]

# filter1 = stat_word_1 != 46
# filter2 = stat_word_12 != 19121
# filter3 = stat_word_9 != 45110
# filter4 = stat_word_12 != 18536
# filter5 = stat_word_12 != 54906
# filter6 = lvdt_stat_b != 83
# filter7 = stat_word_12 != 63675
# filter8 = stat_word_13 != 19585
# filter9 = stat_word_16 != 14372

# pix_gal = pix_gal[
#     filter1
#     & filter2
#     & filter3
#     & filter4
#     & filter5
#     & filter6
#     & filter7
#     & filter8
#     & filter9
# ]
# sky = sky[
#     filter1
#     & filter2
#     & filter3
#     & filter4
#     & filter5
#     & filter6
#     & filter7
#     & filter8
#     & filter9
# ]

# frequency mapping
fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

# scan_mode = 0  # SS
channel = 3  # LL
frec = {}
f_ghz = {}
for mode, item in modes.items():
    frec[mode] = 4 * (channel % 2) + item

    f_icm = np.arange(len(sky[mode][0])) * (fnyq["icm"][frec[mode]] / 320)
    c = 3e8 * 1e2  # cm/s
    f_ghz[mode] = (
        f_icm * c * 1e-9 + 55
    )  # this might not be right but it is what matches the initial frequencies of the firas movie

NSIDE = 32
npix = hp.nside2npix(NSIDE)

# print(f"sky shape: {sky.shape} and pix_gal shape: {pix_gal.shape}")

print("Calculating BB curve")

for mode in modes.keys():
    bb_curve = np.zeros(len(f_ghz[mode]))
    for freq in range(len(f_ghz[mode])):
        hpxmap = np.zeros(npix)
        data_density = np.zeros(npix)

        for i in range(len(pix_gal)):
            hpxmap[pix_gal[mode][i]] += sky[mode][i][freq]
            data_density[pix_gal[mode][i]] += 1

        m = np.zeros(npix)
        mask2 = data_density == 0
        m[~mask2] = hpxmap[~mask2] / data_density[~mask2]

        # mask the map with the points of no data
        m[mask2] = hp.UNSEEN
        # mask the galaxy
        m[mask] = hp.UNSEEN

        # average all points at this frequency to get each frequency point in the bb curve
        bb_curve[freq] = np.mean(m[m != hp.UNSEEN])

    # print(f"Fitting BB curve: {bb_curve}")

    t0 = np.array(2.726)
    fit = minimize(residuals, t0, args=(f_ghz[mode][1:], bb_curve[1:]))

    print(f"Plotting BB curve: {fit.x[0]}")

    plt.plot(f_ghz[mode][1:], bb_curve[1:], label="Data")
    plt.plot(f_ghz[mode][1:], planck(f_ghz[mode], fit.x[0])[1:], label="Fit")
    plt.plot(f_ghz[mode][1:], planck(f_ghz[mode], t0)[1:], label="Original")
    plt.xlabel("Frequency [GHz]")
    plt.ylabel("Brightness [MJy/sr]")
    plt.title("BB curve")
    plt.legend()
    plt.savefig(f"../../output/plots/bb_curve_{mode}.png")
    plt.clf()
