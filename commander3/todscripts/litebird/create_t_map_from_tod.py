import h5py
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt


#toddir = '/home/eirik/data/litebird_sims/'
toddir = "/mn/stornext/u3/ragnaur/data/tut/Commander3_LB_TOD/TODS_test/L1-060"
todfile = ["LB_L1-060_000021.h5", "LB_L1-060_000022.h5", "LB_L1-060_000023.h5", "LB_L1-060_000024.h5", "LB_L1-060_000025.h5", "LB_L1-060_000026.h5", "LB_L1-060_000027.h5", "LB_L1-060_000028.h5", "LB_L1-060_000029.h5", "LB_L1-060_000030.h5", "LB_L1-060_000031.h5", "LB_L1-060_000032.h5", "LB_L1-060_000033.h5", "LB_L1-060_000034.h5", "LB_L1-060_000035.h5", "LB_L1-060_000036.h5", "LB_L1-060_000037.h5", "LB_L1-060_000038.h5", "LB_L1-060_000039.h5", "LB_L1-060_000040.h5"]

nside = 512
outmap = np.zeros((48, 12 * nside ** 2))
nhits = np.zeros((48, 12 * nside ** 2))

for file in todfile:
    outtod = h5py.File(f"{toddir}/{file}")
    scans = [scan for scan in list(outtod) if scan != 'common']
    #print(list(outtod[scans[0]]))

    det_names = [det for det in list(outtod[scans[0]]) if det != 'common'] 
    #print(scans)

    for scan in scans:
        print(scan)
        for i, detector in enumerate(det_names):
            #print(detector)
            #currscan_psi = outtod[scan][detector]['psi']
            #currscan_phi = currscan_psi[0,0,
            #currscan_theta = outtod[scan][detector]['theta']
            #currscan_pixs = hp.ang2pix(nside, currscan_theta, currscan_phi)
            currscan_pixs = outtod[scan][detector]['pix']
            #print(currscan_pixs)
            currscan_tod  = outtod[scan][detector]['tod']
            nhits[i, currscan_pixs] += 1
            outmap[i, currscan_pixs] += currscan_tod

print(np.average(outmap))
print("Test message")
map_out = np.sum(outmap,0)
print(np.average(map_out))
nhits_out = np.sum(nhits,0)

#map_out[np.where(nhits_out != 0)] /= nhits_out[np.where(nhits_out != 0)]
map_out[np.where(nhits_out == 0)] = hp.UNSEEN

print(np.average(map_out))

hp.mollview(map_out, min=0, max=0.1)
plt.savefig(
            f"/mn/stornext/u3/ragnaur/data/tut/Commander3_LB_TOD/TODS_test/sim_20_days_2.png", dpi=800
            )
plt.show()
plt.clf()

"""
for i in range(48):
    hp.mollview(outmap[i])
    plt.savefig(
            f"/mn/stornext/u3/ragnaur/data/tut/Commander3_LB_TOD/TODS_test/scan_sim_{todfile}_{i:03}.png", dpi=800
            )
    plt.clf()
"""
