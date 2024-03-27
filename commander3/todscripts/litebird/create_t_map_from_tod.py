import h5py
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
#from astropy import units
#from astropy.coordinates import SkyCoord
#from astropy.coordinates import BarycentricMeanEcliptic
from cosmoglobe.tod_tools import TODLoader

#toddir = '/mn/stornext/d22/cmbco/litebird/e2e_ns512/sim0000/detectors_LFT_L1-060_T+B/tods'
toddir = "/mn/stornext/u3/ragnaur/data/tut/Commander3_LB_TOD/TODS_test/L1-060/"
todfile = ["LB_L1-060_000051", "LB_L1-060_000050", "LB_L1-060_000053", "LB_L1-060_000054", "LB_L1-060_000055", "LB_L1-060_000056", "LB_L1-060_000057"]#, "LB_L1-060_000058", "LB_L1-060_000059"]#, "LB_L1-060_000010.h5", "LB_L1-060_000011.h5", "LB_L1-060_000012.h5", "LB_L1-060_000013.h5", "LB_L1-060_000014.h5", "LB_L1-060_000015.h5", "LB_L1-060_000016.h5", "LB_L1-060_000017.h5", "LB_L1-060_000018.h5", "LB_L1-060_000019.h5", "LB_L1-060_000020.h5"]
#todfile= ["LB_L1-060_000051.h5", "LB_L1-060_000052.h5", "LB_L1-060_000003.h5", "LB_L1-060_000004.h5", "LB_L1-060_000005.h5", "LB_L1-060_000006.h5"]
#todfile = ['LB_LFT_60_obs_rank0000.hdf5']#, 'LB_LFT_60_obs_rank0001.hdf5', 'LB_LFT_60_obs_rank0002.hdf5', 'LB_LFT_60_obs_rank0003.hdf5']
nside = 512
outmap = np.zeros((48, 12 * nside ** 2))
nhits = np.zeros((48, 12 * nside ** 2))

"""
for file in todfile:
    print(file)
    outtod = h5py.File(f"{toddir}/{file}")
    pointings = outtod['pointings']
    currscan_tod  = outtod['tod_fg']
    colat     = pointings[:,:,0]
    lat       = colat - np.pi/2
    lon       = pointings[:,:,1]
    #lbs.coordinates._rotate_coordinates_e2g_for_all(
    #c_ecl     = SkyCoord(lon=lon, lat=lat, unit='rad', frame='barycentricmeanecliptic')
    #c_gal     = c_ecl.transform_to('galactic')
    #theta     = c_gal.l.deg
    #phi       = c_gal.b.deg
    pixels    = hp.ang2pix(nside, theta, phi, lonlat=True)
    currscan_tod  = outtod['tod_fg']
    nhits[:, pixels] += 1
    outmap[:, pixels] += currscan_tod
    
"""
for file in todfile:
    outtod = h5py.File(f"{toddir}{file}.h5")
    scans = [scan for scan in list(outtod) if scan != 'common']
    #print(list(outtod[scans[0]]))
    det_names = [det for det in list(outtod[scans[0]]) if det != 'common'] 
    #print(scans)
    comm_tod = TODLoader(toddir, "")
    comm_tod.init_file(file, "")
    for scan in scans:
        print(scan)
        for i, detector in enumerate(det_names):
            #print(detector)
            #currscan_psi = outtod[scan][detector]['psi']
            #currscan_phi = currscan_psi[0,0,
            #currscan_theta = outtod[scan][detector]['theta']
            #currscan_pixs = hp.ang2pix(nside, currscan_theta, currscan_phi
            currscan_pixs = comm_tod.load_field(f"{scan}/{detector}/pix").astype("float")[()]
            #No comnpression: 
            #currscan_pixs = outtod[scan][detector]['pix']
            #print(currscan_pixs)
            currscan_tod  = outtod[scan][detector]['tod']
            nhits[i, currscan_pixs.astype(int)] += 1
            outmap[i, currscan_pixs.astype(int)] += currscan_tod

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
            f"/mn/stornext/u3/ragnaur/data/tut/Commander3_LB_TOD/TODS_test/sim_9days_uncompressing.png", dpi=800
            )
plt.show()
plt.clf()
"""
hp.mollview(nhits_out, min=0, max=200)
plt.savefig(
            f"/mn/stornext/u3/ragnaur/data/tut/Commander3_LB_TOD/TODS_test/sim_20_days_nhits.png", dpi=800
            )
plt.show()
plt.clf()
"""
"""
for i in range(48):
    hp.mollview(outmap[i])
    plt.savefig(
            f"/mn/stornext/u3/ragnaur/data/tut/Commander3_LB_TOD/TODS_test/scan_sim_{todfile}_{i:03}.png", dpi=800
            )
    plt.clf()
"""
