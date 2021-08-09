import numpy as np
import matplotlib.pyplot as plt

from glob import glob

fnames = glob('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_WMAP_all/tod_K113_pid000021_samp*.dat')

fnames.sort()



delta_t = 1.536/12


fnames.sort()
##data = np.loadtxt('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_WMAP_all/tod_K113_pid000021_samp000048.dat')
#data = np.loadtxt('/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_WMAP_all/tod_K113_pid000021_samp000050.dat')
# Sample   uncal_TOD (mK)  n_corr (mK) cal_TOD (mK)  sky (mK)   s_orb (mK), mask, baseline, sl, bp, gain, sigma0
f_ind = 0
for fname in fnames:
  data = np.loadtxt(fname)

  t, d, ncorr, cal, s_totA, s_orbA, s_totB, s_orbB, mask, baseline, sl, bp, gain, sigma0 = data.T
  
  ind = (t < 5000)
  t *= delta_t
  
  plt.figure(figsize=(16, 8))
  plt.plot(t/3600/24, d - gain*(s_totA - s_totB) - baseline)
  plt.plot(t/3600/24, ncorr)
  plt.xlim([6,7])
  plt.xlabel('Time [days]')
  plt.savefig(f'plots/{f_ind}.png')
  f_ind += 1
  continue
  
  #fig, axes = plt.subplots(figsize=(16,8), sharex=True, nrows=2)
  fig, axes = plt.subplots(figsize=(16,3*16/10), sharex=True, nrows=2)
  axes[0].plot(t[ind], d[ind], 'k.', ms=2, label='TOD')
  axes[0].plot(t[ind], gain[ind]*(s_totA[ind] - s_totB[ind]) + baseline[ind], label=r'$gs+b$')
  axes[0].plot(t[ind], ncorr[ind] + baseline[ind], label=r'$n_\mathrm{corr}+b$')
  axes[0].set_ylabel('Data [du]')
  
  axes[0].set_ylim(baseline[0] - 15, baseline[0] + 15)
  axes[0].legend(loc='best')
  
  #plt.plot(t[ind], d[ind] - gain[ind]*(s_totA[ind] - s_totB[ind]) - baseline[ind] - ncorr[ind], 'k.')
  axes[1].plot(t[ind], d[ind] - gain[ind]*(s_totA[ind] - s_totB[ind]) - baseline[ind] - ncorr[ind], 
               'k.', ms=2)
  axes[1].set_ylim(-4*sigma0[0], 4*sigma0[0])
  axes[1].set_xlabel(r'Time [s]')
  axes[1].set_ylabel('Residual [du]')
  plt.savefig('tod_res.png', bbox_inches='tight', transparent=True, dpi=300)
  
  ind = (mask == 1) & ( t < 500000)
  plt.figure(figsize=(16,5))
  #plt.plot(t[ind], d[ind] - gain[ind]*(s_totA[ind] - s_totB[ind]) - baseline[ind] - ncorr[ind], 'k.')
  plt.plot(t[ind], (d[ind] - gain[ind]*(s_totA[ind] - s_totB[ind]) - baseline[ind] - ncorr[ind])/sigma0[ind], 'k.', ms=3, alpha=0.05)
  plt.ylim(-4, 4)
  plt.xlabel(r'Time [s]')
  plt.ylabel('Residual/$\sigma_0$ [unitless]')
  
  plt.savefig('tod_res_long.png', bbox_inches='tight', transparent=True)
  
  #ind = (mask == 1)
  #delta = d[ind] -gain[ind]*(s_totA[ind] - s_totB[ind]) - baseline[ind] - ncorr[ind]
  
  #plt.figure()
  #plt.hist(delta/sigma0[0], bins=np.linspace(-6, 6, 101), density=True)
  #
  #x = np.linspace(-6, 6,10001)
  #plt.plot(x,np.exp(-x**2/2)/np.sqrt(2*np.pi))
  #plt.xlabel(r'$\delta/\sigma_0$')
  #plt.yscale('log')

plt.close('all')
