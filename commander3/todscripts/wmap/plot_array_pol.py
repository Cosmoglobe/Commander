import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import healpy as hp

import cmasher as cmr

mpl.rcParams["figure.figsize"][0] *= 2
mpl.rcParams["figure.figsize"][1] *= 2

nside = 512
band = "K1"
data = hp.ud_grade(
    hp.read_map(
        f"data/wmap_iqusmap_r9_9yr_{band}_v5.fits", field=(0, 1, 2), verbose=False
    ),
    nside,
)

deltas = np.load("all_samples_delta_pol.npy")
bla = np.load("all_samples_pol.npy")
monos = []
dipos = []
for i in range(len(deltas)):
    hp.mollview(
        np.split(bla[i], 4)[0],
        min=-3.5,
        max=3.5,
        cmap=cmr.pride,
        title="",
        sub=(2, 2, 1),
    )
    hp.mollview(
        np.split(bla[i], 4)[1],
        min=-0.35,
        max=0.35,
        cmap=cmr.pride,
        title="",
        sub=(2, 2, 2),
    )
    hp.mollview(
        np.split(bla[i], 4)[2],
        min=-0.35,
        max=0.35,
        cmap=cmr.pride,
        title="",
        sub=(2, 2, 3),
    )
    hp.mollview(
        np.split(bla[i], 4)[3],
        min=-0.35,
        max=0.35,
        cmap=cmr.pride,
        title="",
        sub=(2, 2, 4),
    )
    ax = plt.gca()
    ax.set_title(str(np.round(deltas[i], 6)).zfill(6), loc="right")
    plt.savefig(f"plots/sol_{str(i).zfill(3)}")
    plt.close("all")

"""
Bash command to make video
scp -r uio:Commander/commander3/todscripts/wmap/plots/sol_???.png .

ffmpeg -r 10 -f image2 -s 1920x1080 -start_number 1 -i sol_%03d.png -vframes 1000 -vcodec libx264  -crf 25  -pix_fmt yuv420p test.mp4
"""
