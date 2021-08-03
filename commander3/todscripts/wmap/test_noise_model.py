import numpy as np
import matplotlib.pyplot as plt

from glob import glob
from astropy.io import fits
from astropy.table import Table
from scipy.fft import fft, fftfreq

def get_N(t, fsamp, AC, a, b, c, d):
    N = a + b*np.log10(abs(t)) + c*np.log10(abs(t))**2 + d*np.log10(abs(t))**3
    N[t < 1/fsamp] = AC
    N[N < 0] = 0
    return N

data_dir = '/mn/stornext/d16/cmbco/ola/wmap/tods/optimal_filters/'
fnames = glob(data_dir+'wmap_opt_time_domain_filters_yr?_v5.fits')
fnames.sort()

Nobs = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])
Nobs_array = np.zeros(2*len(Nobs))
Nobs_array[::2] = Nobs
Nobs_array[1::2] = Nobs


fig, axes = plt.subplots(nrows=5, ncols=4, sharex=True, sharey=True)
axs = axes.flatten()
for i in range(len(fnames)):
    hdu_list = fits.open(fnames[i])
    for j in range(1, len(hdu_list)):
        data = hdu_list[j].data
        cols = hdu_list[j].columns

        N = len(data[cols.names[0]])
        T = 1/1.536/Nobs_array[j]

        xf = fftfreq(N,T)[:N//2]

        yf = fft(data[cols.names[1]])
        axs[2*(j-1)].loglog(xf[1:],1/np.abs(yf[1:N//2]), color=plt.cm.viridis(i/len(fnames)))
        axs[2*(j-1)].set_title(cols.names[1])

        yf = fft(data[cols.names[2]])
        axs[2*(j-1)+1].loglog(xf[1:],1/np.abs(yf[1:N//2]), color=plt.cm.viridis(i/len(fnames)))
        axs[2*(j-1)+1].set_title(cols.names[2])

        plt.figure('test'+str(i))
        plt.loglog(1/np.abs(yf[1:N//2]), color='k', lw=1)

plt.show()
