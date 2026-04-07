import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import healpy as hp
from numpy.lib.stride_tricks import sliding_window_view

DIR = '/mn/stornext/d5/data/duncanwa/IRAS'

def run_iras_readout(utc1, plate, sop, obs, npts, det, overlap, chunk):
    plate = plate[chunk]
    sop = sop[chunk]
    obs = obs[chunk]
    npts = npts[chunk]
    utc1 = utc1[chunk]
    L1_data = '/mn/stornext/d16/cmbco/ola/IRAS/IRAS_Level1/data'
    t_tot = []
    lon_tot = []
    lat_tot = []
    flux_tot = []
    for i in range(len(utc1)):
        fname = f'{L1_data}/plate{plate[i]:04}/sop{sop[i]:03}/obs{obs[i]:03}/det{det:02}.tbl'
        lons = []
        lats = []
        fluxs = []
        with open(fname,'r') as f:
            for line in f:
                if 'jband' in line:
                    if '1' in line:
                        dt = 1/4
                    elif '2' in line:
                        dt = 1/8
                    elif ('3' in line) or ('4' in line):
                        dt = 1/16
                    else:
                        print('Where is this even happening?')
                elif '\\' in line:
                    continue
                elif '|' in line:
                    continue
                else:
                    data = line.split()
                    if (float(data[0]) == -999.) | (float(data[1]) == -999.) | (float(data[2]) == -1e20):
                        lons.append(np.nan)
                        lats.append(np.nan)
                        fluxs.append(np.nan)
                    else:
                        lons.append(float(data[0]))
                        lats.append(float(data[1]))
                        fluxs.append(float(data[2]))

        t = np.arange(npts[i])*dt + utc1[i]

        lon_tot  = lon_tot + lons
        lat_tot  = lat_tot + lats
        flux_tot = flux_tot + fluxs
        t_tot    = t_tot + t.tolist()
        ok = True
    
    flux_tot = np.array(flux_tot)#.astype(np.float32)
    lon_tot = np.array(lon_tot)#.astype(np.float32)
    lat_tot = np.array(lat_tot)#.astype(np.float32)
    t_tot = np.array(t_tot)#.astype(np.float4)

    np.save(f'{DIR}/new_merged_data/iras_chunk{chunk+1:04}_data_{det:02}.npy', np.array([t_tot, lon_tot, lat_tot, flux_tot]))


dt = 0.2500128573368466
dt = 0.25

# https://irsa.ipac.caltech.edu/IRASdocs/issa.exp.sup/Ap/B.html
# 16, 16, 8, and 4 samples per second for 12, 25, 60, and 100 microns.


'''
for i in range(len(d1) - 10):
    if np.all(d2[:10,2] == d1[i:i+10,2]):
        print(i)
'''



num_chunks=64

num_chunks = 256

num_chunks = 512

num_chunks = 1024

for det in range(1, 62):
    data = np.loadtxt(f'{DIR}/det{det:02}_files.txt')
   
    try:
        utc1, plate, sop, obs, npts = data.T
    except ValueError:
        continue
    
    inds = np.argsort(utc1)
    utc1   = utc1[inds]
    plate  = plate[inds].astype(int)
    sop    = sop[inds].astype(int)
    obs    = obs[inds].astype(int)
    npts   = npts[inds].astype(int)
    
    
    
    overlap = 50
    
    import multiprocessing
    NUM_PROCS = multiprocessing.cpu_count()
    
    
    with multiprocessing.Pool(processes=NUM_PROCS) as pool:
    
        utc1i = np.array_split(utc1,num_chunks)
        platei = np.array_split(plate,num_chunks)
        sopi = np.array_split(sop, num_chunks)
        obsi = np.array_split(obs, num_chunks)
        nptsi = np.array_split(npts, num_chunks)
        multiple_results = [
                pool.apply_async(
                    run_iras_readout,
                    args=(utc1i, platei, sopi, obsi, nptsi, det, overlap, i),
                    )
            for i in range(len(utc1i))]
        for i in tqdm(range(len(multiple_results))):
            if multiple_results[i]:
                multiple_results[i].get()
