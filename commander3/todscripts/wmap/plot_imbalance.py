import matplotlib.pyplot as plt

import numpy as np
import healpy as hp
from glob import glob

from astropy.io import fits

fnames = glob("data/wmap_loss_imbalance_template*.fits")
fnames.sort()

for f in fnames:
    print(f)
    fil = fits.open(f)
    Q, U = np.split(fil[1].data["TEMPERATURE"], 2)
    Q = hp.reorder(Q, n2r=True)
    U = hp.reorder(U, n2r=True)
    lim = 2 * Q.std()
    plt.figure()
    hp.mollview(
        Q,
        min=-lim,
        max=lim,
        title=r"$Q$",
        nest=False,
        cbar=False,
        cmap="coolwarm",
        sub=121,
    )
    hp.mollview(
        U,
        min=-lim,
        max=lim,
        title=r"$U$",
        nest=False,
        cbar=False,
        cmap="coolwarm",
        sub=122,
    )
    plt.suptitle(f)
    m = np.array([0 * Q, Q, U])
    hp.write_map(f[:-5] + "_iqu.fits", m)
    # plt.savefig(f[:-5] + '.png')

plt.show()
