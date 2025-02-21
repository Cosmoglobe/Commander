'''

 According to the calibrator design paper (Mather 1999) and the design document
 (Mather 1993), the GRT temperatures were read out within a single major
 telemetry frame (one every 32 seconds), which were all cycled through.

 This means there is some sort of timing offset that needs to be accounted for.
 The timing in the data, from COBETRIEVE (ct) is the average of the associated
 science frames, so unless I'm missing something fundamental, the original
 timing of the engineering frames doesn't exist anymore. Probably some sort of
 packet loss.

 Currently I'm trying ot figure this out by looking at very specific regions
 where the temperature shifts quite rapidly, and the lag is very obvious.

'''

def plot_sdf(sci_mode, t_min, t_max, vmin=None, vmax=None):

    time = sci_mode['ct_head/time'][()]
    ifgs = sci_mode['ifg_data/ifg'][()].astype(float)
    xcal = sci_mode['dq_data/xcal_pos'][()]
    inds = (time > t_min) & (time < t_max)

    time = time[inds]
    ifgs = ifgs[inds]
    xcal = xcal[inds]

    ifgs[xcal != 1] = np.nan

    med = np.median(ifgs, axis=1)
    ifgs = ifgs.T
    ifgs -=  med


    plt.imshow(ifgs, extent=[time[0], time[-1], 0, 511], vmin=vmin, vmax=vmax,
            aspect='auto', interpolation='none')
    plt.show()

    return

# cubic interpolation clearly fits to noise fluctuations
kind = 'cubic'
kind = 'linear'


import h5py
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.units as u

from time import time

from astropy.visualization import time_support
time_support()
# https://docs.astropy.org/en/latest/api/astropy.visualization.time_support.html

# FIRAS mission periods, MJD in ns
t_00 = (Time('1989:326:11:30').mjd*u.day).to(u.ns) # First light
t_01 = (Time('1989:328:00:00').mjd*u.day).to(u.ns) # First ICAL nulling
t_02 = (Time('1989:343:01:52').mjd*u.day).to(u.ns) # MTM uses position mode through SAA
t_03 = (Time('1990:019:02:05').mjd*u.day).to(u.ns) # Horns commanded from 2.70 to 2.75 K
t_04 = (Time('1990:080:01:15').mjd*u.day).to(u.ns) # MTM uses power off through SAA
t_05 = (Time('1990:129:00:00').mjd*u.day).to(u.ns) # Eclipse season starts
t_06 = (Time('1990:139:15:35').mjd*u.day).to(u.ns) # Horns commanded to 6 K
t_07 = (Time('1990:193:18:50').mjd*u.day).to(u.ns) # Horns commanded to 4 K
t_08 = (Time('1990:207:11:04').mjd*u.day).to(u.ns) # Sky horn calibration, XCAL out
t_09 = (Time('1990:208:11:20').mjd*u.day).to(u.ns) # Horns commanded to final temperature
t_10 = (Time('1990:220:05:00').mjd*u.day).to(u.ns) # XCAL placed under temperature control
t_11 = (Time('1990:264:09:36').mjd*u.day).to(u.ns) # Final period over


lens_grt =  [
                  1,1,1,1,1,
                  4,1,4,1,1,
                  1,1,1,1,1,
                  4,1,4,1,1,
                  1,1,1,1,1,
                  4,1,4,1,1,
                  1,1,1,1,1,
                  4,1,4,1,1]

# From Appendix H NFS_HKP record (iNgest FIRAS Stripper Housekeeping)
# There are 4 of:
# ex_cal, skyhorn, refhorn, ical, dihedral, bol_assembloy (4),
# mirror, cal_resist (4), xcal, collimator

names_en_analog_grt = [
                   # grts
                   # 'a_lo_grt', 'a_hi_grt', 'b_lo_grt', 'b_hi_grt',
                   'a_lo_xcal_tip', 'a_lo_skyhorn', 'a_lo_refhorn', 'a_lo_ical', 'a_lo_dihedral',
                   'a_lo_bol_assem', 'a_lo_mirror', 'a_lo_cal_resistors', 'a_lo_xcal_cone', 'a_lo_collimator',
                   'a_hi_xcal_tip', 'a_hi_skyhorn', 'a_hi_refhorn', 'a_hi_ical', 'a_hi_dihedral',
                   'a_hi_bol_assem', 'a_hi_mirror', 'a_hi_cal_resistors', 'a_hi_xcal_cone', 'a_hi_collimator',
                   'b_lo_xcal_tip', 'b_lo_skyhorn', 'b_lo_refhorn', 'b_lo_ical', 'b_lo_dihedral',
                   'b_lo_bol_assem', 'b_lo_mirror', 'b_lo_cal_resistors', 'b_lo_xcal_cone', 'b_lo_collimator',
                   'b_hi_xcal_tip', 'b_hi_skyhorn', 'b_hi_refhorn', 'b_hi_ical', 'b_hi_dihedral',
                   'b_hi_bol_assem', 'b_hi_mirror', 'b_hi_cal_resistors', 'b_hi_xcal_cone', 'b_hi_collimator']





data = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng.h5')
sdf = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5')



grts = {}
offsets = {}

# "Binary time"
# ADT Time, in units of 100ns since 1858-11-17
bin_time = data['fdq_eng']['ct_head']['time']*(100*u.ns)

dn = 8
'''
# GMT
gmt_time = data['fdq_eng']['ct_head']['gmt']
times = []
from tqdm import tqdm
for i in tqdm(range(len(gmt_time))):
    t_str = gmt_time[i].decode('UTF-8')
    yr = f'19{t_str[0:2]}'
    day = t_str[2:5]
    hr = t_str[5:7]
    mi = t_str[7:9]
    sec = t_str[9:11]
    dec = t_str[11:]

    t = f'{yr}:{day}:{hr}:{mi}:{sec}.{dec}'
    times.append(t)
time = Time(times, format='yday')
'''

ind0 = 0
for i in range(len(names_en_analog_grt)):
    grts[names_en_analog_grt[i]] = data['fdq_eng']['en_analog']['grt'][:,ind0:ind0+lens_grt[i]].flatten()
    offsets[names_en_analog_grt[i]] = ind0
    ind0 += lens_grt[i]
time = bin_time



dt = 32*u.s / 64

dt = 1*u.s

inds = (np.arange(len(time)) > 527_240) & (np.arange(len(time)) < 527_340)

t0 = (4.1577e18 + 90e12)*u.ns
t1 = (4.1577e18 + 97e12)*u.ns
inds = (time > t0) & (time < t1)

inds = (np.arange(len(time)) > 527_000) & (np.arange(len(time)) < 550_000)
inds = (np.arange(len(time)) > 400_000) & (np.arange(len(time)) < 600_000)

t0 = (4.1582e18+ 2e13)*u.ns
t1 = (4.1582e18+ 4e13)*u.ns
t0 = 4.158223338391165e+18*u.ns
t1 = 4.158230899109607e+18*u.ns
inds = (time > t0) & (time < t1)
'''
t0=4.1568030989473644e+18*u.ns
t1=4.15681471777588e+18*u.ns
inds = (time > t0) & (time < t1)
'''
# In the second period
t0 = (4.1349e18+2.3e13)*u.ns
t1 = (4.1349e18+2.6e13)*u.ns
inds = (time > t_02) & (time < t_03)
# In the third period
t0 = (4.1368e18+6.5e13)*u.ns
t1 = (4.1368e18+7.5e13)*u.ns
inds = (time > t_03) & (time < t_04)
# In the fourth period
t0 = (4.142e18)*u.ns
t1 = (4.142e18+4e13)*u.ns

# Fifth period
t0 = (4.1473e18+5.5e13)*u.ns
t1 = (4.1473e18+6e13)*u.ns


# Sixth period
t0 = (4.1495e18 + 8e13)*u.ns
t1 = (4.1495e18 + 10e13)*u.ns

# Seventh period
t0 = (4.15338e18 + 5e12)*u.ns
t1 = (4.15338e18 + 6e12)*u.ns

t0 = t_07
t1 = t_08
# Eighth period
t0 = (4.1553e18 + 3e13)*u.ns
t1 = (4.1553e18 + 4e13)*u.ns


t0 = t_08
t1 = t_09

t0 = (4.1557e18 + 5.5e13)*u.ns
t1 = (4.1557e18 + 6.5e13)*u.ns


t0 = t_09
t1 = t_10
# Tenth period
t0 = (4.15655e18)*u.ns
t1 = (4.15655e18 + 1e13)*u.ns

t0 = t_10
t1 = t_11
# 11th period
t0 = (4.1594e18 + 3.6e13)*u.ns
t1 = (4.1594e18 + 4.6e13)*u.ns


t0 = t_10
t1 = t_11
t0 = (4.158e18)*u.ns
t1 = (4.159e18)*u.ns
t0 = 4.1582058e+18*u.ns
t1 = 4.1582062e+18*u.ns
t0 = 4.159640058589277e+18*u.ns
t1 = 4.159660079727482e+18*u.ns
inds = (time > t0) & (time < t1)

#inds = (np.arange(len(time)) > 550_000) & (np.arange(len(time)) < 600_000)

#inds = np.ones(len(time), dtype=bool)


grt_names = ['xcal_tip', 'skyhorn', 'refhorn', 'ical', 'dihedral',
        'mirror', 'xcal_cone', 'collimator']
for side in ['a', 'b']:
    for i, grt in enumerate(grt_names):
        not_ok = grts[f'{side}_lo_{grt}'] == -9999
        grts[f'{side}_lo_{grt}'][not_ok] = np.nan
        not_ok = grts[f'{side}_hi_{grt}'] == -9999
        grts[f'{side}_hi_{grt}'][not_ok] = np.nan

fig, axes = plt.subplots(8, 4, sharex=True, sharey=False, figsize=(12, 12))
axs = axes.flatten()
for i, grt in enumerate(grt_names):
    axs[4*i].plot(time[inds], grts[f'a_lo_{grt}'][inds], '.', label='Low current reading')
    axs[4*i+1].plot(time[inds], grts[f'a_hi_{grt}'][inds], '.', label='High current reading')

    axs[4*i+2].plot(time[inds], grts[f'b_lo_{grt}'][inds], '.', label='Low current reading')
    axs[4*i+3].plot(time[inds], grts[f'b_hi_{grt}'][inds], '.', label='High current reading')
    axs[4*i].set_ylabel(grt)
axs[0].set_title('a, low')
axs[1].set_title('a, high')
axs[2].set_title('b, low')
axs[3].set_title('b, high')
plt.savefig('temperature_readings.png')

plt.figure()
ll = sdf['fdq_sdf_ll']
plot_sdf(ll, t0.value/100, t1.value/100, vmin=-100, vmax=100)


plt.show()
plt.close()
asdf

from scipy.interpolate import interp1d

for _, j in enumerate(np.arange(-32, 32, dn)):
    t_lo = (time[inds])
    t_hi = (time[inds] + j*dt)
    fig, axes = plt.subplots(4, 4, sharex=False, sharey=False, figsize=(12, 10))
    axs = axes.flatten()
    for i, grt in enumerate(grt_names):
        T_lo = grts[f'a_lo_{grt}'][inds]
        T_hi = grts[f'a_hi_{grt}'][inds]
        T_lo[T_lo < 0] = np.nan
        T_hi[T_hi < 0] = np.nan

        f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
        f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)
        
        times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()),
                100*len(t_lo))
        inds_ = (f1(times) > 0)
        axs[2*i].plot(f1(times)[inds_], f2(times)[inds_] - f1(times)[inds_], '.', ms=1)

        T_lo = grts[f'b_lo_{grt}'][inds]
        T_hi = grts[f'b_hi_{grt}'][inds]
        T_lo[T_lo < 0] = np.nan
        T_hi[T_hi < 0] = np.nan

        f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
        f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)
        
        times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()),
                100*len(t_lo))
        inds_ = (f1(times) > 0)
        axs[2*i+1].plot(f1(times)[inds_], f2(times)[inds_] - f1(times)[inds_], '.', ms=1)
        axs[2*i].set_title(f"{grt}, a")
        axs[2*i+1].set_title(f"{grt}, b")
    fig.supylabel(r'$T_\mathrm{high} - T_\mathrm{low}$')
    fig.supxlabel(r'$T_\mathrm{low}$')
    plt.suptitle(f'{j} second offset')
    plt.tight_layout()
    plt.savefig(f'offsets_diff_{_:03}.png')
    plt.close()

for _, j in enumerate(np.arange(-32, 32, dn)):
    t_lo = (time[inds])
    t_hi = (time[inds] + j*dt)
    fig, axes = plt.subplots(4, 4, sharex=False, sharey=False, figsize=(12, 10))
    axs = axes.flatten()
    for i, grt in enumerate(grt_names):
        T_lo = grts[f'a_lo_{grt}'][inds]
        T_hi = grts[f'a_hi_{grt}'][inds]
        T_lo[T_lo < 0] = np.nan
        T_hi[T_hi < 0] = np.nan

        f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
        f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)
        
        times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()),
                100*len(t_lo))
        inds_ = (f1(times) > 0)
        axs[2*i].plot(f1(times)[inds_], f2(times)[inds_], '.', ms=1)

        T_lo = grts[f'b_lo_{grt}'][inds]
        T_hi = grts[f'b_hi_{grt}'][inds]
        T_lo[T_lo < 0] = np.nan
        T_hi[T_hi < 0] = np.nan

        f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
        f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)
        
        times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()),
                100*len(t_lo))
        inds_ = (f1(times) > 0)
        axs[2*i+1].plot(f1(times)[inds_], f2(times)[inds_], '.', ms=1)
        axs[2*i].set_title(f"{grt}, a")
        axs[2*i+1].set_title(f"{grt}, b")
    fig.supylabel(r'$T_\mathrm{high}$')
    fig.supxlabel(r'$T_\mathrm{low}$')
    plt.suptitle(f'{j} second offset')
    plt.tight_layout()
    plt.savefig(f'offsets_split_{_:03}.png')
    plt.close()
