import h5py
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt


toddir = '/home/eirik/data/litebird_sims/'
todfile = 'litebird_60_000001.h5'
outtod = h5py.File(f"{toddir}{todfile}")
nside = 512
outmap = np.zeros((48, 12 * nside ** 2))
nhits = np.zeros((48, 12 * nside ** 2))

scans = [scan for scan in list(outtod) if scan != 'common']

for scan in scans:
    print(scan)
    for i, detector in enumerate(list(outtod[scan])):
        currscan_phi = outtod[scan][detector]['phi']
        currscan_theta = outtod[scan][detector]['theta']
        currscan_pixs = hp.ang2pix(nside, currscan_theta, currscan_phi)
        currscan_tod = outtod[scan][detector]['tod']
        nhits[i, currscan_pixs] += 1
        outmap[i, currscan_pixs] += currscan_tod
outmap[np.where(nhits != 0)] /= nhits[np.where(nhits != 0)]
outmap[np.where(nhits == 0)] = hp.UNSEEN

for i in range(48):
    hp.mollview(outmap[i])
    plt.savefig(f'/home/eirik/temp/sim_{i:03}.png', dpi=800)
    plt.clf()
