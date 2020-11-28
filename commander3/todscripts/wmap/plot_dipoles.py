'''
#
#		Nominal radiometer data
#		The gain is in du/mK, the offset is in du, the noise is in mK.
#		These 4x10 arrays are in the order, following CJ's notation:
#		[K1-13, K1-14, K1-23, K1-24,
#		 Ka1-13, Ka1-14, Ka1-23, Ka1-24,
#		 Q1-13, Q1-14, Q1-23, Q1-24,
#		 Q2-13, Q2-14, Q2-23, Q2-24,
#		 V1-13, V1-14, V1-23, V1-24,
#		 V2-13, V2-14, V2-23, V2-24,
#		 W1-13, W1-14, W1-23, W1-24,
#		 W2-13, W2-14, W2-23, W2-24,
#		 W3-13, W3-14, W3-23, W3-24,
#		 W4-13, W4-14, W4-23, W4-24]
#
#               Gains and baselines given signs and values from the first
#		attempt at cal_map for Pass 1.  These are median values
#		from the hourly calibration files.
#
GAIN = [ -0.9700,  0.9938,  1.1745, -1.1200, &
          0.8668, -0.8753, -1.0914,  1.0033, &
          1.0530, -0.9834,  0.4914, -0.5365, &
         -0.9882,  1.0173, -0.8135,  0.7896, &
          0.4896, -0.5380, -0.5840,  0.5840, &
         -0.4948,  0.4872,  0.4096, -0.3802, &
          0.3888, -0.4139,  0.3290, -0.3003, &
         -0.3587,  0.3701,  0.3655, -0.3666, &
         -0.3255,  0.3517, -0.3291,  0.3225, &
          0.2841, -0.2918,  0.3796, -0.3591 ]
'''

import numpy as np
import matplotlib.pyplot as plt

import healpy as hp
from astropy.io import fits
from glob import glob

from get_gain_model import get_gain


file_num = 30

prefix = '/mn/stornext/d16/cmbco/ola/wmap/'
files = glob(prefix + 'tods/uncalibrated/*.fits')
files.sort()
data = fits.open(files[file_num])

version=24

allbands = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']

x_im = np.array([-0.00067, 0.00536, 0.00353, 0.00154,
                 -0.00013, 0.00414, 0.00756, 0.00986,
                  0.00053, 0.00250, 0.00352, 0.00245,
                  0.01134, 0.00173, 0.01017, 0.01142,
                 -0.00122, 0.00463, 0.02311, 0.02054])
x_bar = (x_im[::2] + x_im[1::2])/2
dx_bar = (x_im[::2] - x_im[1::2])/2

xbar = {b:x for b, x in zip(allbands, x_bar)}
dx_bar = {b:dx for b, dx in zip(allbands, dx_bar)}

band = 'K1'

#band = 'V1'

xbar = xbar[band]
dxbar = dx_bar[band]



labels = [f'{band}13', f'{band}14', f'{band}23', f'{band}24']
Ks = []
for l in labels:
    Ks.append(data[2].data[l].flatten())
Ks = np.array(Ks)
#gains = np.array([-0.9700, 0.9938, 1.1745, -1.1200])
gains = np.array([get_gain(data, b)[1][0] for b in labels])
baselines = np.array([32136.98, 31764.96, 31718.19, 32239.29])
#cal = [(Ks[i] - baselines[i])/gains[i] for i in range(4)]
cal = np.array([(Ks[i] - np.median(Ks[i]))/gains[i] for i in range(4)])
print(Ks.shape)

time = data[2].data['time']
dt = np.diff(time)[0]*len(time)/len(cal[0])
time = np.arange(0,dt*(len(cal[0])), dt)
# time in JD
time = time*(24*60)
# time in minutes

d1 = 0.5*(cal[0] + cal[1])
d2 = 0.5*(cal[2] + cal[3])

d = 0.5*(d1 + d2)
p = 0.5*(d1 - d2)

n1 = 0.5*(cal[0] - cal[1])
n2 = 0.5*(cal[2] - cal[3])


n_d = 0.5*(n1 + n2)
n_p = 0.5*(n1 - n2)

data = fits.open('/mn/stornext/u3/hke/xsan/commander3/v2/data_BP8/cmb_init_md_tempdiponly.fits')
i_sol = hp.ud_grade(data[1].data['TEMPERATURE'].flatten(), 512)
q_sol = hp.ud_grade(data[1].data['Q_POLARISATION'].flatten(), 512)
u_sol = hp.ud_grade(data[1].data['U_POLARISATION'].flatten(), 512)
sol = np.array([i_sol, q_sol, u_sol, 0*u_sol])

import h5py
from glob import glob
import huffman
from cg_solver import make_dipole
fnames = glob(f'/mn/stornext/d16/cmbco/bp/wmap/data/wmap_{band}_*v{version}.h5')
fnames.sort()
fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
axs = axes.flatten()
for n in range(4):
    fname = np.random.choice(fnames) 
    print(fname)
    f= h5py.File(fname, 'r')
    DAs = [[], [], [], []]
    sigmas = []
    pixA = []
    pixB = []
    psiA = []
    psiB = []
    
    for i in range(len(list(f.keys()))-1):
        obsid = str(list(f.keys())[i])
        labels = [f'{band}13', f'{band}14',f'{band}23',f'{band}24']
        
        huffTree = f[obsid+'/common/hufftree']
        huffSymb = f[obsid+'/common/huffsymb']
        h = huffman.Huffman(tree=huffTree, symb=huffSymb)
        
        
        
        
        npsi = 2048
        psiBins = np.linspace(0, 2*np.pi, npsi)
        gains = np.zeros(len(labels))
        for num, label in enumerate(labels):
            TODs = np.array(f[obsid + '/' + label + '/tod'])
            scalars = f[obsid + '/' + label + '/scalars']
            gains[num] = scalars[0]
            TODs = TODs - np.median(TODs)
            DAs[num] = DAs[num] + TODs.tolist()
            sigmas.append(TODs.std())
            if label == f'{band}13':
                pixA.append(h.Decoder(np.array(f[obsid + '/' + label + \
                    '/pixA'])).astype('int').tolist())
                pixB.append(h.Decoder(np.array(f[obsid + '/' + label + \
                    '/pixB'])).astype('int').tolist())
                psiA.append(psiBins[h.Decoder(np.array(f[obsid + '/' + label + \
                    '/psiA'])).astype('int')])
                psiB.append(psiBins[h.Decoder(np.array(f[obsid + '/' + label + \
                    '/psiB'])).astype('int')])
    
    pixA = np.array(pixA).flatten()
    pixB = np.array(pixB).flatten()
    psiA = np.array(psiA).flatten()
    psiB = np.array(psiB).flatten()


    d_sol = np.zeros(len(pixA))
    d_solA = np.zeros(len(pixA))
    d_solB = np.zeros(len(pixA))
    for t in range(len(pixA)):
        d_solA[t] = (1+xbar)*(i_sol[pixA[t]])
        d_solB[t] = -(1-xbar)*(i_sol[pixB[t]])

    axs[n].plot(time, d_solA, label='Horn A')
    axs[n].plot(time, d_solB, label='Horn B')
    axs[n].set_title(fname.split('/')[-1])

    bla = np.array([time, d_solA, d_solB])
    np.savetxt(fname.split('/')[-1][:-3:]+'.txt', bla)



plt.show()
