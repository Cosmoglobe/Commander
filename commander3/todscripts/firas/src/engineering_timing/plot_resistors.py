import h5py
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter1d as filt

sigma = 5

data = h5py.File('fdq_eng_new.h5')
data.keys()
ah_cal = data['en_analog/grt/a_hi_cal_resistors'][()]
bh_cal = data['en_analog/grt/b_hi_cal_resistors'][()]
al_cal = data['en_analog/grt/a_lo_cal_resistors'][()]
bl_cal = data['en_analog/grt/b_lo_cal_resistors'][()]


fig, axes = plt.subplots(sharex=True, nrows=2, ncols=2)
axs = axes.flatten()
for i in range(4):
    axs[i].plot(filt(al_cal[:,i], sigma))
    axs[i].plot(filt(bl_cal[:,i], sigma))


axs[0].set_ylim([1.0225e4, 1.0225e4 + 50])
axs[1].set_ylim([1.429e4, 1.429e4 + 50])
axs[2].set_ylim([1.5625e4, 1.5625e4 + 50])
axs[3].set_ylim([1.596e4, 1.596e4 + 50])
plt.suptitle('Low Current')
plt.tight_layout()
plt.savefig('lo_cal.png', bbox_inches='tight')

plt.close('all')

fig, axes = plt.subplots(sharex=True, nrows=2, ncols=2)
axs = axes.flatten()
for i in range(4):
    axs[i].plot(ah_cal[:,i])
    axs[i].plot(bh_cal[:,i])


axs[0].set_ylim([4900, 4900 +50])
axs[1].set_ylim([1.0250e4, 1.0250e4 + 50])
axs[2].set_ylim([1.43e4, 1.43e4 + 50])
axs[3].set_ylim([1.5325e4, 1.5325e4 + 50])
plt.suptitle('High Current')
plt.tight_layout()
plt.savefig('hi_cal.png', bbox_inches='tight')
plt.close('all')
