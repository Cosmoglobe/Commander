import h5py
import numpy as np
import matplotlib.pyplot as plt

'''
Trying to find the circles when plotting t_hi versus t_lo, and hte source of
this effect.
'''


data = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5')

t_lo = data['en_analog/grt/b_lo_xcal_cone'][:,0]
t_hi = data['en_analog/grt/b_hi_xcal_cone'][:,0]

t_lo = data['en_analog/grt/b_lo_ical'][:,0]
t_hi = data['en_analog/grt/b_hi_ical'][:,0]

time = np.arange(len(t_lo))

inds = (time > 5.2e5) & ( time < 5.6e5)

inds = time > 0

time = time[inds]

t_lo = t_lo[inds]
t_hi = t_hi[inds]

# circles


# Moving averages


window = 25  # size of the window
t_lo_w = np.lib.stride_tricks.sliding_window_view(t_lo, window)
t_lo_sd = np.std(t_lo_w, axis=-1)/window**0.5
t_lo_mu = np.mean(t_lo_w, axis=-1)/window**0.5

t_hi_w = np.lib.stride_tricks.sliding_window_view(t_hi, window)
t_hi_sd = np.std(t_hi_w, axis=-1)
t_hi_mu = np.mean(t_hi_w, axis=-1)

time_w = np.lib.stride_tricks.sliding_window_view(time, window)
time_mu = np.mean(time_w, axis=-1)


plt.figure('time')

plt.plot(time, t_lo)
plt.plot(time, t_hi)

plt.figure('sd')

plt.plot(time_mu, t_lo_sd)
plt.plot(time_mu, t_hi_sd)

plt.figure('signal-to-noise, low')
plt.plot(time_mu, t_lo_mu/t_lo_sd, '.', alpha=0.5)
plt.yscale('log')
plt.figure('signal-to-noise, hi')
plt.plot(time_mu, t_hi_mu/t_hi_sd, '.', alpha=0.5)
plt.yscale('log')




from scipy.interpolate import interp1d

f1 = interp1d(time_mu, t_hi_mu/t_hi_sd, fill_value='extrapolate')
f2 = interp1d(time_mu, t_lo_mu/t_lo_sd, fill_value='extrapolate')

inds = (f1(time) > 50) & (f2(time) > 50)


plt.figure('time')
plt.plot(time[inds], t_lo[inds])
plt.plot(time[inds], t_hi[inds])


plt.figure('circles')
plt.plot(t_lo, t_hi - t_lo)
t_lo[~inds] = np.nan
t_hi[~inds] = np.nan
plt.plot(t_lo, t_hi - t_lo)


plt.show()
