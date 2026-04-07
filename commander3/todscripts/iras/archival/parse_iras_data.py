import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm  import tqdm
import multiprocessing

# See https://irsa.ipac.caltech.edu/IRASdocs/level_1_format.html


L1_data = '/mn/stornext/d16/cmbco/ola/IRAS/IRAS_Level1/data'

from glob import glob

'''
#for i in range(29, 601):
for i in range(29,51):
    # Not sure how many there should be here
    for j in range(10):
        fnames = glob(f'{L1_data}/*/sop{i:03}/obs{j:03}')
        print(i, j,  len(fnames))

'''

fnames = glob(f'{L1_data}/*/sop03?/*/det01.tbl')
fnames.sort()
print(len(fnames))

sopobs = np.zeros(len(fnames))

'''
for i, f in enumerate(fnames):

    ra, dec, flux = np.loadtxt(f, skiprows=14).T
    #plt.plot(flux)
    header = np.loadtxt(f, dtype=str, max_rows=11)
    t1 = header[-2,2].astype(float)
    t2 = header[-1,2].astype(float)
    npnts = header[4, 2].astype(int)
    sop = header[0,2].astype(int)
    obs = header[1,2].astype(int)
    sopobs[i] = t1
    #print(t1)
    #print(t2)
    #print(npnts)
    #print(header)
    #t = np.linspace(t1, t2, npnts)
    #plt.plot(t, flux, 'k')
#plt.show()


new_order = np.argsort(sopobs)

for i in range(5):
    header = np.loadtxt(fnames[new_order[i]], dtype=str, max_rows=11)
    t1 = header[-2,2].astype(float)
    t2 = header[-1,2].astype(float)
    npnts = header[4, 2].astype(int)
    sop = header[0,2].astype(int)
    obs = header[1,2].astype(int)
    print(t1, sop, obs, fnames[new_order[i]])
'''



# I think I need to do a one-time operation of writing down the utc1 and utc2 of
# each file in the nested directory.
# Also need to do it for each detector independently
det = 1

def write_file(det):
    if (det == 1):
        import time
        t0 = time.time()
    with open(f'/mn/stornext/d5/data/duncanwa/IRAS/det{det:02}_files.txt', 'w+') as the_file:
        the_file.write(f'# utcs1         \t plate \t sop\t obs\t npts\n')
        for plate in range(2000):
            for sop in range(601):
                for obs in range(100):
                    fname = f'{L1_data}/plate{plate:04}/sop{sop:03}/obs{obs:03}/det{det:02}.tbl'
                    if os.path.exists(fname):
                        with open(fname, 'r') as f:
                            for line in f:
                                if 'jdet' in line:
                                    jdet = line.split('=')[1].strip()
                                elif 'npts' in line:
                                    npts = int(line.split('=')[1])
                                elif 'utcs1' in line:
                                    utcs1 = float(line.split('=')[1])
                                elif 'utcs2' in line:
                                    break
                        the_file.write(f'{utcs1:8.6f}\t\t{plate:04}\t{sop:03}\t{obs:03}\t{npts}\n')
            if (det == 1):
                dt = (time.time() - t0)/(plate+1)
                time_to_go = dt*(2000 - plate - 1)//60
                print(f'plate {plate} out of 2000, {time_to_go} minutes to go {(time.time() - t0)//60} spent')
    return



#for det in range(1, 63):
#    write_file(det)

N_PROC = multiprocessing.cpu_count()
N_PROC = 62
with multiprocessing.Pool(processes=N_PROC) as pool:
    proc_chunks = [
        pool.apply_async(
            write_file, 
            (det,),
        )
        for det in np.arange(1, 63)
    ]
    for p in proc_chunks:
        p.get()
        


'''
sop   = 559
obs   = 3
jdet  = 1
jband = 4
npts  = 68
idir  = 1
iflag = 0
dist  = -299.579865
angle = 137.041061
utcs1 = 89355220.873004
utcs2 = 89355237.623900
'''
