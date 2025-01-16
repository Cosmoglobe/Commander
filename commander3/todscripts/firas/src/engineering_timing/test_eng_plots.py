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


lens_grt =  [
                  1,1,1,1,1,
                  4,1,4,1,1,
                  1,1,1,1,1,
                  4,1,4,1,1,
                  1,1,1,1,1,
                  4,1,4,1,1,
                  1,1,1,1,1,
                  4,1,4,1,1]

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
                   'b_hi_bol_assem', 'b_hi_mirror', 'b_hi_cal_resistors',
                   'b_hi_xcal_cone', 'b_hi_collimator']


data = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng.h5')



grts = {}
offsets = {}

# "Binary time"
bin_time = data['fdq_eng']['ct_head']['time']

# "Binary time"
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
ind0 = 0
for i in range(len(names_en_analog_grt)):
    if 'xcal_cone' in names_en_analog_grt[i]:
        grts[names_en_analog_grt[i]] = data['fdq_eng']['en_analog']['grt'][:,ind0:ind0+lens_grt[i]].flatten()
        offsets[names_en_analog_grt[i]] = ind0
    ind0 += lens_grt[i]



dt = 32*u.s / 64

dt = 1*u.s

inds = (np.arange(len(time)) > 527_240) & (np.arange(len(time)) < 527_340)
#inds = (np.arange(len(time)) > 527_000) & (np.arange(len(time)) < 550_000)

plt.figure()
plt.plot(time[inds] - dt*offsets['b_lo_xcal_cone'], grts['b_lo_xcal_cone'][inds], '.')
plt.plot(time[inds] - dt*offsets['b_hi_xcal_cone'], grts['b_hi_xcal_cone'][inds], '.')


t_lo = (time[inds] - dt*offsets['b_lo_xcal_cone']).mjd
t_hi = (time[inds] - dt*offsets['b_hi_xcal_cone']).mjd


T_lo = grts['b_lo_xcal_cone'][inds]
T_hi = grts['b_hi_xcal_cone'][inds]



from scipy.interpolate import interp1d

plt.figure()
print(t_lo.shape, T_lo.shape)
f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
print(t_hi.shape, T_hi.shape)
f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)

times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()),
        100*len(t_lo))
print(times, f1(times))
plt.plot(t_lo, T_lo, 'C0.')
plt.plot(t_hi, T_hi, 'C1.')
plt.plot(times, f1(times), 'C3')
plt.plot(times, f2(times), 'C4')

plt.figure()
plt.plot(T_lo, T_hi - T_lo)
plt.plot(f1(times), f2(times) - f1(times))
plt.figure()
plt.plot(f1(times), f2(times) - f1(times))


plt.figure()

t_min = 0.8
t_max = 1.2

#for i in range(0, t_max, delta_i):
for i in np.linspace(t_min, t_max, 10):
    dt = i*u.s
    t_lo = (time[inds] - dt*offsets['b_lo_xcal_cone']).mjd
    t_hi = (time[inds] - dt*offsets['b_hi_xcal_cone']).mjd
    
    
    T_lo = grts['b_lo_xcal_cone'][inds]
    T_hi = grts['b_hi_xcal_cone'][inds]
    
    
    f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
    f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)
    
    times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()), 10000)
    plt.plot(f1(times), f2(times) - f1(times),
            color=plt.cm.coolwarm((i-t_min)/(t_max-t_min)))
    
plt.show()
