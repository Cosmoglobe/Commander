import numpy as np
import matplotlib.pyplot as plt

from glob import glob
from astropy.io import fits
from astropy.table import Table
from scipy.fft import fft, fftfreq
from scipy import interpolate

def get_N(t, fsamp, AC, a, b, c, d):
    N = a + b*np.log10(abs(t)) + c*np.log10(abs(t))**2 + d*np.log10(abs(t))**3
    N[t < 1/fsamp] = AC
    N[N < 0] = 0
    return N

data_dir = '/mn/stornext/d16/cmbco/ola/wmap/tods/optimal_filters/'
data_dir = 'data/tmp/new_'
fnames = glob(data_dir+'wmap_opt_time_domain_filters_yr?_v5.fits')
fnames.sort()

Nobs = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])
Nobs_array = np.zeros(2*len(Nobs))
Nobs_array[::2] = Nobs
Nobs_array[1::2] = Nobs


fig, axes = plt.subplots(nrows=5, ncols=4, sharex=True, sharey=False,
    figsize=(12, 16))
axs = axes.flatten()
fig1, axes1 = plt.subplots(nrows=5, ncols=4, sharex=True, sharey=False,
    figsize=(12, 16))
axs1 = axes1.flatten()
ind = 0

fknees = np.array([0.4, 0.51, 0.71, 0.32, 1.09, 0.35, 5,76, 8.62, 0.09, 1.41,
  0.88, 8.35, 7.88, 0.66, 9.02, 7.47, 0.93, 0.28, 46.5, 26.0])*1e-3
for i in range(len(fnames)):
    hdu_list = fits.open(fnames[i])
    for j in range(1, len(hdu_list)):
        data = hdu_list[j].data
        cols = hdu_list[j].columns

        x_in = data[cols.names[0]]



        axs1[2*(j-1)].plot(data[cols.names[0]], data[cols.names[1]],
            color=plt.cm.viridis(i/len(fnames)))
        axs1[2*(j-1)].set_xlim([-5,5])
        d = data[cols.names[1]]*1.01
        axs1[2*(j-1)].set_ylim([d.min(), -d.min()])
        axs1[2*(j-1)+1].plot(data[cols.names[0]], data[cols.names[2]],
            color=plt.cm.viridis(i/len(fnames)))
        axs1[2*(j-1)+1].set_xlim([-5,5])
        d = data[cols.names[2]]*1.01
        axs1[2*(j-1)+1].set_ylim([d.min(), -d.min()])

        N = len(data[cols.names[0]])
        T = 1/1.536/Nobs_array[j]

        xf = fftfreq(N,T)[:N//2]

        yf = fft(data[cols.names[1]])
        axs[2*(j-1)].loglog(xf[1:],1/np.abs(yf[1:N//2]), color=plt.cm.viridis(i/len(fnames)))
        if i == 0:
            y = (xf[1:]/fknees[ind])**-1
            axs[2*(j-1)].loglog(xf[1:], y/y[0]/np.abs(yf[1]) + 1,
                color='r', zorder=3)
            ind += 1
            axs[2*(j-1)].set_title(cols.names[1])
            axs1[2*(j-1)].set_title(cols.names[1])

        yf = fft(data[cols.names[2]])
        axs[2*(j-1)+1].loglog(xf[1:],1/np.abs(yf[1:N//2]), color=plt.cm.viridis(i/len(fnames)))
        if i == 0:
            y = (xf[1:]/fknees[ind])**-1
            axs[2*(j-1)+1].loglog(xf[1:],
                y/y[0]/np.abs(yf[1]) + 1,
                color='r', zorder=3)
            ind += 1
            axs[2*(j-1)+1].set_title(cols.names[2])
            axs1[2*(j-1)+1].set_title(cols.names[2])

        #plt.figure('test'+str(i))
        #plt.loglog(1/np.abs(yf[1:N//2]), color='k', lw=1)

fig.tight_layout()
fig1.tight_layout()
fig.savefig('test1.png', bbox_inches='tight', dpi=300)
fig1.savefig('test2.png', bbox_inches='tight', dpi=300)
plt.show()
