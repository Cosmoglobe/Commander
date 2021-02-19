import numpy as np
import matplotlib.pyplot as plt
import h5py

from glob import glob

fnames = glob('/mn/stornext/u3/hke/quiet_data/level3/Q/ces/patch_gc/*.hdf')

for f in fnames:
  d = h5py.File(f, 'r')
  
  
  t = d['time'][...]
  tg = d['time_gain'][...]
  tod = d['tod'][...] # 76 x N_TOD
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
  
  gain = d['gain']
  
  
  #plt.plot(t, tod[0])
  #plt.figure()
  #plt.plot(t, phi[0])
  #plt.figure()
  #plt.plot(t, theta[0])
  #plt.figure()
  #plt.plot(t, psi[0])
  
  
  import healpy as hp
  
  nside = 256
  
  phi = phi[0]
  theta = theta[0]
  psi = psi[0]
  tod = tod[0]
  gain = gain[0]
  
  
  pixels = hp.ang2pix(nside, theta, phi)
  
  # Ax = b
  
  
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
  print(tod)
  print(psi)
  
  for t in range(len(pixels)):
    m[pixels[t]]    += 1
    b_T[pixels[t]]  += tod[t]
    b_Q[pixels[t]]  += tod[t]*np.cos(2*psi[t])
    b_U[pixels[t]]  += tod[t]*np.sin(2*psi[t])
  
    A_TT[pixels[t]] += 1
    A_TQ[pixels[t]] += np.cos(2*psi[t])
    A_TU[pixels[t]] += np.sin(2*psi[t])
    A_QQ[pixels[t]] += np.cos(2*psi[t])**2
    A_UU[pixels[t]] += np.sin(2*psi[t])**2
    A_QU[pixels[t]] += np.cos(2*psi[t])*np.sin(2*psi[t])

hp.gnomview(m,  reso=10)
hp.gnomview(b_Q,  reso=10)
hp.gnomview(b_U,  reso=10)

determ = A_QQ*A_UU - A_QU**2

Qhat = (b_Q*A_UU - b_U*A_QU)/determ
Uhat = (b_U*A_QQ - b_Q*A_QU)/determ
hp.gnomview(Qhat,  reso=10, min=-5e-4, max=5e-4, cmap='RdBu_r')
hp.gnomview(Uhat,  reso=10, min=-5e-4, max=5e-4, cmap='RdBu_r')


# tod = q*cos(2*psi) + u*sin(2*psi) + n_corr + n_w
# per module


plt.show()
