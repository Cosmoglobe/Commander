import matplotlib.pyplot as plt
import numpy as np

import h5py

#data = h5py.File('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_WMAP_all/chain_c0001.h5', 'r')

bands=['023-WMAP_K']
burn = 3
thin = 1

#data = h5py.File('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_WMAP_beamtest/chain_c0001.h5', 'r')
#data = h5py.File('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_WMAP_full_spec_test/chain_c0001.h5', 'r')
#data = h5py.File('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_WMAP_bp_Qband/chain_c0001.h5', 'r')
data = h5py.File('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_sl_test/chain_c0001.h5', 'r')
#
bands=[#'023-WMAP_K', 
       #'030-WMAP_Ka',
       '040-WMAP_Q1']
       #'040-WMAP_Q2']
       #'060-WMAP_V1',
       #'060-WMAP_V2']
       #'090-WMAP_W1',
       #'090-WMAP_W2',
       #'090-WMAP_W3',
       #'090-WMAP_W4']
#burn = 25

x_imw9 = {}
#x_imw9['023-WMAP_K'] = [-0.00067, 0.00536]
#x_imw9['030-WMAP_Ka'] = [0.00353, 0.00154]
x_imw9['040-WMAP_Q1'] = [-0.00013, 0.00414]
#x_imw9['040-WMAP_Q2'] = [0.00756, 0.00986]
#x_imw9['060-WMAP_V1'] = [0.00053, 0.00250]
#x_imw9['060-WMAP_V2'] = [0.00352, 0.00245]
#x_imw9['090-WMAP_W1'] = [0.01134, 0.00173]
#x_imw9['090-WMAP_W2'] = [0.01017, 0.01142]
#x_imw9['090-WMAP_W3'] = [-0.00122, 0.00463]
#x_imw9['090-WMAP_W4'] = [0.02311, 0.02054]
x_imw9u = {}
#x_imw9u['023-WMAP_K'] = [0.00017, 0.00014]
#x_imw9u['030-WMAP_Ka'] = [0.00014, 0.00008]
x_imw9u['040-WMAP_Q1'] = [0.00046, 0.00025]
#x_imw9u['040-WMAP_Q2'] = [0.00052, 0.00115]
#x_imw9u['060-WMAP_V1'] = [0.00020, 0.00057]
#x_imw9u['060-WMAP_V2'] = [0.00033, 0.00098]
#x_imw9u['090-WMAP_W1'] = [0.00199, 0.00036]
#x_imw9u['090-WMAP_W2'] = [0.00216, 0.00121]
#x_imw9u['090-WMAP_W3'] = [0.00062, 0.00041]
#x_imw9u['090-WMAP_W4'] = [0.00380, 0.00202]
for band in bands:
  gain0s = [[],[],[],[],[]]
  x_ims = [[],[]]
  samps = []
  for i in range(burn, len(data.keys())-1, thin):
    gain0 = data[str(i).zfill(6) + '/tod/'+band+'/gain0'][:]
    x_im = data[str(i).zfill(6) + '/tod/'+band+'/x_im'][:]
    for j in range(5):
      gain0s[j].append(gain0[j])
    for j in range(2):
      x_ims[j].append(x_im[j])
    samps.append(i)

  fig, axes = plt.subplots(nrows=7, sharex=True)
  for j in range(5):
    axes[j].plot(samps, gain0s[j])
  for j in range(2):
    axes[j+5].plot(samps, x_ims[j])
    ylim = axes[j+5].get_ylim()
    axes[j+5].axhline(x_imw9[band][j] - x_imw9u[band][j], color='k',
        linestyle='--')
    axes[j+5].axhline(x_imw9[band][j] + x_imw9u[band][j], color='k',
        linestyle='--')
    axes[j+5].set_ylim(ylim)
  plt.suptitle(band)

  axes[0].set_ylabel(r'$g_0$')
  axes[1].set_ylabel(r'$g_1$')
  axes[2].set_ylabel(r'$g_2$')
  axes[3].set_ylabel(r'$g_3$')
  axes[4].set_ylabel(r'$g_4$')
  axes[5].set_ylabel(r'$x^\mathrm{im}_1$')
  axes[6].set_ylabel(r'$x^\mathrm{im}_2$')
  axes[6].set_xlabel('Gibbs sample')

plt.show()

