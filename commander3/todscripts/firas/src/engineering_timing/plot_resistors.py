import h5py
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.time import Time

from scipy.ndimage import gaussian_filter1d

from astropy.visualization import time_support
time_support()



def filt(data, sigma):
    if sigma is None:
        return data
    else:
        return gaussian_filter1d(data, sigma)




sigma = None

data = h5py.File('fdq_eng_new.h5')
data.keys()
ah_cal = data['en_analog/grt/a_hi_cal_resistors'][()]
bh_cal = data['en_analog/grt/b_hi_cal_resistors'][()]
al_cal = data['en_analog/grt/a_lo_cal_resistors'][()]
bl_cal = data['en_analog/grt/b_lo_cal_resistors'][()]

ah_bol = data['en_analog/grt/a_hi_bol_assem'][()]
bh_bol = data['en_analog/grt/b_hi_bol_assem'][()]
al_bol = data['en_analog/grt/a_lo_bol_assem'][()]
bl_bol = data['en_analog/grt/b_lo_bol_assem'][()]
bin_time = data['ct_head/time']*(100*u.ns)

time = Time(bin_time.to('day'), format='mjd')



inds = np.arange(len(al_cal[:,0]))


#inds = np.arange(200_000, 400_000)
inds = np.arange(350_000, 351_000)
#inds = np.arange(350_400, 350_650)
#inds = np.arange(350_525, 350_580)

#inds = np.arange(345_000, 355_000)
print(time[inds].to_value('yday'))
print(bin_time[inds][0].to('s'))
print(bin_time[inds][-1].to('s'))

#time = time.to_value('yday')
fig, axes = plt.subplots(sharex=True, nrows=2, ncols=2, figsize=(12, 8))
axs = axes.flatten()
for i in range(4):
    axs[i].plot(time[inds], filt(al_cal[inds,i], sigma),  color='C0')
    axs[i].plot(time[inds], filt(bl_cal[inds,i], sigma),  color='C1')


#axs[0].set_ylim([1.0225e4, 1.0225e4 + 50])
#axs[1].set_ylim([1.429e4, 1.429e4 + 50])
#axs[2].set_ylim([1.5625e4, 1.5625e4 + 50])
#axs[3].set_ylim([1.596e4, 1.596e4 + 50])
plt.suptitle('Low Current')
plt.tight_layout()
plt.savefig('lo_cal.png', bbox_inches='tight')

plt.close('all')

fig, axes = plt.subplots(sharex=True, nrows=2, ncols=2, figsize=(12, 8))
axs = axes.flatten()
for i in range(4):
    axs[i].plot(time[inds], filt(ah_cal[inds,i], sigma),  color='C0')
    axs[i].plot(time[inds], filt(bh_cal[inds,i], sigma),  color='C1')


#axs[0].set_ylim([4900, 4900 +50])
#axs[1].set_ylim([1.0250e4, 1.0250e4 + 50])
#axs[2].set_ylim([1.43e4, 1.43e4 + 50])
#axs[3].set_ylim([1.5325e4, 1.5325e4 + 50])
plt.suptitle('High Current')
plt.tight_layout()
plt.savefig('hi_cal.png', bbox_inches='tight')
plt.close('all')



fig, axes = plt.subplots(sharex=True, nrows=2, ncols=2, figsize=(12, 8))
axs = axes.flatten()
for i in range(4):
    axs[i].plot(time[inds], filt(al_bol[inds,i], sigma), ms=1, color='C0')
    axs[i].plot(time[inds], filt(bl_bol[inds,i], sigma), ms=1, color='C1')
    #axs[i].plot(inds, bl_bol[inds,i]/al_bol[inds,i], '.')


plt.suptitle('Low Current, bolometer')
plt.tight_layout()
plt.savefig('lo_cal_bol.png', bbox_inches='tight')
plt.close('all')
