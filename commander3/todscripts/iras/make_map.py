det_100A = [1, 2, 3, 4, 5, 6, 7]
det_60A  = [8,9,10,11,12,13,14,15]
det_25A  = [16,18,19,21,22]
det_12A  = [23,24,25,26,27,28,29,30]
det_60B  = [31,32,33,34,35,37,38]
det_25B  = [39,40,41,42,43,44,45,46]
det_12B  = [47,48,49,50,51,52,53,54]
det_100B = [55,56,57,58,59,60,61,62]



# solid angle based on NGC6543 scans, table IV.A.1
sA_100A = [14.5, 12.7, 13.0, 11.53, 12, 12.4, 12.6]
sA_100B = [7.1, 14.0, 13.2, 11.2, 11.7, 13.3, 13.5, 10.6]

sA_60A = [7.2, 6.7, 6.6, 2.8, 4.3, 6.6, 6.1, 6.2]
sA_60B = [2.1, 6.4, 5.9, 6.5, 6.3, 6.6, 3.9]


sA_25A = [3.5, 3.6, 2.8, 2.8, 3.1]
sA_25B = [1.4, 3.1, 3.1, 3.4, 3.2, 3.2, 3.2, 2.4]


sA_12A = [2.9, 3, 3.2, 1.2, 2, 3.1, 2.5, 2.8]
sA_12B = [0.77, 3.1, 2.9, 3, 2.7, 2.5, 2.8, 2]


import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

r = hp.Rotator(coord=['C','G'])

from numba import njit

from glob import glob
from tqdm import tqdm

nside = 256
b_map = np.zeros(12*nside**2, dtype='float32')
N_hits = np.zeros(12*nside**2, dtype=int)

@njit()
def bin_map(b_map, N_hits, pix, tod):
    for i in range(len(pix)):
        b_map[pix[i]] += tod[i]
        N_hits[pix[i]] += 1
    return b_map, N_hits

#for det in tqdm(det_12A + det_12B):
#for det in tqdm(det_25A + det_25B):
#for det in tqdm(det_60A + det_60B):
dt = 1/4 # 12, 25
dt = 1/8 # 60
dt = 1/16 # 100
for sA, det in tqdm(zip(sA_100A+sA_100B, det_100A + det_100B)):
#for sA, det in tqdm(zip(sA_60A+sA_60B, det_60A + det_60B)):
#for sA, det in tqdm(zip(sA_25A+sA_25B, det_25A + det_25B)):
#for sA, det in tqdm(zip(sA_12A+sA_12B, det_12A + det_12B)):
    b_map = np.zeros(12*nside**2, dtype='float32')
    N_hits = np.zeros(12*nside**2, dtype=int)
    fnames = glob(f'new_merged_data/iras_chunk????_data_{det:02}.npy')
    for f in fnames:
        data = np.load(f)
        pix = np.zeros(len(data[0]), dtype=int)
        #inds = (data[1] != -999.) & (data[3] != -1e20)
        inds = (~np.isnan(data[1])) & (~np.isnan(data[3]))
        data = data[:,inds]
        lon, lat = r(data[1], data[2], lonlat=True)
        pix = hp.ang2pix(nside, lon, lat, lonlat=True)
        b_map, N_hits = bin_map(b_map, N_hits, pix, data[3]/(sA*1e-7)/dt)

    m = b_map/N_hits
    
    #hp.mollview(m, min=0, max=2e-6, title=f'{det:02}')
    hp.mollview(m, min=0, max=1e-5, title=f'{det:02}')
    #hp.mollview(m, min=0, max=1e-4, title=f'{det:02}')
    plt.savefig(f'{det:02}.png', bbox_inches='tight', dpi=150)
    N_hits = N_hits.astype(float)
    N_hits[N_hits == 0] = hp.UNSEEN
    hp.mollview(N_hits, norm='hist',title=f'{det:02}')
    plt.savefig(f'{det:02}_nhits.png', bbox_inches='tight', dpi=150)
    plt.close('all')
    #hp.gnomview(m, rot=(170, -15), norm='hist', reso=5)
    plt.show()
