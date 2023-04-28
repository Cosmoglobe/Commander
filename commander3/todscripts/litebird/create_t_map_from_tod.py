import h5py
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt


#toddir = '/home/eirik/data/litebird_sims/'
toddir = '/mn/stornext/u3/eirikgje/data/litebird_tods/'
#todfile = 'litebird_60_000001.h5'
todfiles = [f'litebird_60_{i:06}.h5' for i in range(1, 144)]
outtods = [f"{toddir}{todfile}" for todfile in todfiles]
#outtod = h5py.File(f"{toddir}{todfile}")
nside = 512
outmap = np.zeros((48, 12 * nside ** 2))
nhits = np.zeros((48, 12 * nside ** 2))

#scans = [(scan for scan in list(h5py.File(outtod)) if scan != 'common') for outtod in outtods]

for outtod in outtods:
    print(outtod)
    outtod_file = h5py.File(outtod)
    scans = list(outtod_file)
    for scan in scans:
        if scan == 'common': continue
        for i, detector in enumerate(list(outtod_file[scan])):
            if detector == 'common': continue
            currscan_phi = outtod_file[scan][detector]['phi']
            currscan_theta = outtod_file[scan][detector]['theta']
            currscan_pixs = hp.ang2pix(nside, currscan_theta, currscan_phi)
#            currscan_tod = outtod_file[scan][detector]['tod']
            outmap[i, currscan_pixs] += outtod_file[scan][detector]['tod']
            nhits[i, currscan_pixs] += 1
#            outmap[i, currscan_pixs] += currscan_tod
outmap[np.where(nhits != 0)] /= nhits[np.where(nhits != 0)]
outmap[np.where(nhits == 0)] = hp.UNSEEN

for i in range(48):
    hp.write_map(f'/mn/stornext/u3/eirikgje/staging/sim_{i:03}.fits', outmap[i], overwrite=True)
#    hp.mollview(outmap[i])
##    plt.savefig(f'/home/eirik/temp/sim_{i:03}.png', dpi=800)
#    plt.savefig(f'/mn/stornext/u3/eirikgje/staging/sim_{i:03}.png', dpi=800)
#    plt.close()
