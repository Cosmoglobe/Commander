det_100A = [1, 2, 3, 4, 5, 6, 7]
det_60A  = [8,9,10,11,12,13,14,15]
det_25A  = [16,18,19,21,22]
det_12A  = [23,24,25,26,27,28,29,30]
det_60B  = [31,32,33,34,35,37,38]
det_25B  = [39,40,41,42,43,44,45,46]
det_12B  = [47,48,49,50,51,52,53,54]
det_100B = [55,56,57,58,59,60,61,62]


import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

from numba import njit

from glob import glob
from tqdm import tqdm

nside = 2048
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
for det in tqdm(det_60A + det_60B):
#for det in tqdm(det_100A + det_100B):
    fnames = glob(f'merged_data/iras_chunk????_data_{det:02}.npy')
    for f in fnames:
        data = np.load(f)
        pix = np.zeros(len(data[0]), dtype=int)
        inds = (data[1] != -999.) & (data[3] != -1e20)
        data = data[:,inds]
        pix = hp.ang2pix(nside, data[1], data[2], lonlat=True)
        b_map, N_hits = bin_map(b_map, N_hits, pix, data[3])

m = b_map/N_hits

hp.mollview(m, norm='hist', coord=['C','G'])
hp.gnomview(m, rot=(170, -15), norm='hist', reso=5)
plt.show()
