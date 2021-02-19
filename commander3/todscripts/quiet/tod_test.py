import numpy as np
import matplotlib.pyplot as plt
import h5py

from glob import glob
import healpy as hp

nside = 256

m   = np.zeros(hp.nside2npix(nside))
b_T = np.zeros(hp.nside2npix(nside))
b_Q = np.zeros(hp.nside2npix(nside))
b_U = np.zeros(hp.nside2npix(nside))

A_TT = np.zeros(hp.nside2npix(nside))
A_TQ = np.zeros(hp.nside2npix(nside))
A_TU = np.zeros(hp.nside2npix(nside))
A_QQ = np.zeros(hp.nside2npix(nside))
A_QU = np.zeros(hp.nside2npix(nside))
A_UU = np.zeros(hp.nside2npix(nside))

fnames = []
fnames += glob('/mn/stornext/u3/hke/quiet_data/level3/Q/ces/patch_2a/*.hdf')
print(len(fnames))
fnames += glob('/mn/stornext/u3/hke/quiet_data/level3/Q/ces/patch_4a/*.hdf')
print(len(fnames))
fnames += glob('/mn/stornext/u3/hke/quiet_data/level3/Q/ces/patch_6a/*.hdf')
print(len(fnames))
fnames += glob('/mn/stornext/u3/hke/quiet_data/level3/Q/ces/patch_7b/*.hdf')
print(len(fnames))
fnames += glob('/mn/stornext/u3/hke/quiet_data/level3/Q/ces/patch_gb/*.hdf')
print(len(fnames))
fnames += glob('/mn/stornext/u3/hke/quiet_data/level3/Q/ces/patch_gc/*.hdf')
print(len(fnames))



from tqdm import tqdm
for f in tqdm(fnames):
  d = h5py.File(f, 'r')
  
  
  t = d['time'][...]
  tg = d['time_gain'][...]
  tod = d['tod'][...] # 76 x N_TOD
  sigma0 = d['sigma0']
  # There are 19 horns (modules)
  # Each module has 4 radiometers (diodes), 1 and 4 have most of the signal
  # Horn 1 diode 1
  # Horn 1 diode 2
  # Horn 1 diode 3
  # Horn 1 diode 4
  # Horn 2 diode 1
  # Horn 2 diode 2
  # Horn 2 diode 3
  # Horn 2 diode 4
  
  point = d['point']
  phi, theta, psi = point[:,:,0], point[:,:,1], point[:,:,2]
  
  pixels = hp.ang2pix(nside, theta, phi)

  for i in range(19):
    for t in range(len(pixels)):
      m[pixels[t]]    += 1
      for j in range(4):
        b_T[pixels[i,t]]  += tod[i+j,t]/sigma0[i+j]**2
        b_Q[pixels[i,t]]  += tod[i+j,t]*np.cos(2*psi[i,t])/sigma0[i+j]**2
        b_U[pixels[i,t]]  += tod[i+j,t]*np.sin(2*psi[i,t])/sigma0[i+j]**2
    
        A_TT[pixels[i,t]] += 1/sigma0[i+j]**2
        A_TQ[pixels[i,t]] += np.cos(2*psi[i,t])/sigma0[i+j]**2
        A_TU[pixels[i,t]] += np.sin(2*psi[i,t])/sigma0[i+j]**2
        A_QQ[pixels[i,t]] += np.cos(2*psi[i,t])**2/sigma0[i+j]**2
        A_UU[pixels[i,t]] += np.sin(2*psi[i,t])**2/sigma0[i+j]**2
        A_QU[pixels[i,t]] += np.cos(2*psi[i,t])*np.sin(2*psi[i,t])/sigma0[i+j]**2

hp.mollview(m)
hp.mollview(b_Q)
hp.mollview(b_U)

determ = A_QQ*A_UU - A_QU**2

That = b_T/A_TT
Qhat = (b_Q*A_UU - b_U*A_QU)/determ
Uhat = (b_U*A_QQ - b_Q*A_QU)/determ

hp.mollview(That, min=-5e-4, max=5e-4)
hp.mollview(Qhat, min=-5e-4, max=5e-4)
hp.mollview(Uhat, min=-5e-4, max=5e-4)
plt.show()

Mhat = np.array([That,Qhat, Uhat])
hp.write_map('quiet_test.fits', Mhat, overwrite=True)


# tod = q*cos(2*psi) + u*sin(2*psi) + n_corr + n_w
# per module


