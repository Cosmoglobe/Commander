import h5py
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt


#toddir = '/home/eirik/data/litebird_sims/'
toddir = "/Users/maksymb/Desktop/litebird_db/maksymb_lb/commander3/todscripts/litebird/test_litebird_sims"
todfile = "LB_L1-060_000001.h5"
outtod = h5py.File(f"{toddir}/{todfile}")
nside = 512
outmap = np.zeros((48, 12 * nside ** 2))
nhits = np.zeros((48, 12 * nside ** 2))

scans = [scan for scan in list(outtod) if scan != 'common']

print(list(outtod[scans[0]]))

det_names = [det for det in list(outtod[scans[0]]) if det != 'common'] 

print(scans)
for scan in scans:
    print(scan)
    for i, detector in enumerate(det_names):
        #currscan_phi = outtod[scan][detector]['phi']
        #currscan_theta = outtod[scan][detector]['theta']
        #currscan_pixs = hp.ang2pix(nside, currscan_theta, currscan_phi)
        currscan_pixs = outtod[scan][detector]['pix']
        print(currscan_pixs)
        currscan_tod  = outtod[scan][detector]['tod']
        print(currscan_tod)
        nhits[i, currscan_pixs] += 1
        outmap[i, currscan_pixs] += currscan_tod
print("Test message")
outmap[np.where(nhits != 0)] /= nhits[np.where(nhits != 0)]
outmap[np.where(nhits == 0)] = hp.UNSEEN

for i in range(48):
    hp.mollview(outmap[i])
    plt.savefig(
            f"/Users/maksymb/Desktop/litebird_db/maksymb_lb/commander3/todscripts/litebird/test_litebird_sims/sim_{i:03}.png", dpi=800
            )
            #'/home/eirik/temp/sim_{i:03}.png', dpi=800)
    plt.clf()
