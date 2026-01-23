# Try just printing out the times! We know that there are 64 possible readings between the two chunks of header data. It should be something like half a second between each of them.

'''

 According to the calibrator design paper (Mather 1999) and the design document
 (Mather 1993), the GRT temperatures were read out within a single major
 telemetry frame (one every 32 seconds), which were all cycled through.

 This means there is some sort of timing offset that needs to be accounted for.
 The timing in the data, from COBETRIEVE (ct) is the average of the associated
 science frames, so unless I'm missing something fundamental, the original
 timing of the engineering frames doesn't exist anymore. Probably some sort of
 packet loss.

 Currently I'm trying to figure this out by looking at very specific regions
 where the temperature shifts quite rapidly, and the lag is very obvious.



So there is a 16 second lag between the high and low current readings, with a bit of a sign flip.

I am also looking at the a and b side comparisons, when possible.

As far as I can tell, this is the order things will be read out.

'a_lo_xcal_tip'
'a_lo_skyhorn'
'a_lo_refhorn'
'a_lo_ical'
'a_lo_dihedral'
'a_lo_bol_assem'
'a_lo_mirror'
'a_lo_cal_resistors'
'a_lo_xcal_cone'
'a_lo_collimator'
'a_hi_xcal_tip'
'a_hi_skyhorn'
'a_hi_refhorn'
'a_hi_ical'
'a_hi_dihedral'
'a_hi_bol_assem'
'a_hi_mirror'
'a_hi_cal_resistors'
'a_hi_xcal_cone'
'a_hi_collimator'
'b_lo_xcal_tip'
'b_lo_skyhorn'
'b_lo_refhorn'
'b_lo_ical'
'b_lo_dihedral'
'b_lo_bol_assem'
'b_lo_mirror'
'b_lo_cal_resistors'
'b_lo_xcal_cone'
'b_lo_collimator'
'b_hi_xcal_tip'
'b_hi_skyhorn'                     : 64 seconds before(?) a_hi skyhorn
'b_hi_refhorn'                     : 12 seconds after a_hi_refhorn
'b_hi_ical'                        : 56 seconds after a_hi_ical
'b_hi_dihedral',                   : 24 seconds before(?) a_hi_dihedral
'b_hi_bol_assem'
'b_hi_mirror'                      : 16 seconds after a_hi_mirror
'b_hi_cal_resistors'
'b_hi_xcal_cone'                   : 4 seconds after a_hi_xcal_cone
'b_hi_collimator'


'''

from scipy.interpolate import interp1d
import h5py
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.units as u

from time import time

from astropy.visualization import time_support
time_support()

from tqdm import tqdm

plt.rcParams["xtick.minor.visible"] =  True

def find_epoch(time):
    '''
    Takes time in seconds since MJD=0, returns index of the epoch.
    '''
    t_inds = (time >= time_periods)
    return np.arange(len(time_periods))[t_inds][-1]

def nt(times, epoch):
    newtime = times - time_periods[epoch]
    return newtime.to('s')

DATA_DIR = '/mn/stornext/d16/cmbco/ola/firas/initial_data'
DATA_DIR = '/home/dwatts/Commander/commander3/todscripts/firas/src/engineering_timing'

def plot_sdf(sci_mode1, sci_mode2, t_min, t_max, vmin=None, vmax=None, eng_time=None,
        eng_data1=None, eng_data2=None, eng_data3=None, eng_xcal=None):

    epoch = find_epoch(t_min)

    time1 = sci_mode1['ct_head/time'][()]*(100*u.ns) - 136*u.s + 42*u.s
    time1 = sci_mode1['ct_head/time'][()]*(100*u.ns) #+ 42*u.s
    time1 = sci_mode1['ct_head/time'][()]*(100*u.ns) - 36*u.s
    time1 = sci_mode1['ct_head/time'][()]*(100*u.ns)
    time1 = time1.to('s')
    time2 = sci_mode2['ct_head/time'][()]*(100*u.ns) #- 136*u.s + 42*u.s
    time2 = sci_mode2['ct_head/time'][()]*(100*u.ns) - 36*u.s
    time2 = sci_mode2['ct_head/time'][()]*(100*u.ns)
    time2 = time2.to('s')


    print('time1 min', np.diff(time1).min())
    print('time2 min', np.diff(time2).min())
    data_ready = sci_mode1['sci_head/data_ready'][()] # (N_ifgs x 8)
                                                     # Flagged if anything != 1
                                                     # See line 206 of
                                                     # fut_get_qualflags.for
    data_qual  = sci_mode1['sci_head/data_qual'][()]  # (N_ifgs x 60)
                                                     # Telemetry quality?
                                                     # See line 172 of
                                                     # ftb_list_sci_time.for
                                                     # Also, line 714, flagged
                                                     # if anything != 0
    # These are telemetry issues, so it's unlikely to be a huge problem.
    mtm_length = sci_mode1['sci_head']['mtm_length'][()]
    mtm_speed  = sci_mode1['sci_head']['mtm_speed'][()]
    scan_mode = mtm_length*2 + mtm_speed
    inds1 = (time1 > t_min) & (time1 < t_max)
    inds2 = (time2 > t_min) & (time2 < t_max)


    time1 = nt(time1[inds1],epoch)
    time2 = nt(time2[inds2],epoch)
    xcal1 = sci_mode1['dq_data/xcal_pos'][()]
    xcal1 = xcal1[inds1]
    xcal2 = sci_mode2['dq_data/xcal_pos'][()]
    xcal2 = xcal2[inds2]
    scan_mode = scan_mode[inds1]




    #fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(12, 8))
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(12, 8))
    axs = axes.flatten()

    if (eng_time is not None):
        axs[0].plot(nt(eng_time,epoch), eng_data3, '.')
        #axs[0].plot(eng_time, eng_data1, '.')
        #axs[1].plot(eng_time, eng_data2)
        axs[1].plot(time1, xcal1, '.')
        axs[1].plot(time2, xcal2, '.')
        axs[1].plot(nt(eng_time,epoch), eng_xcal[:,0], '.')
        axs[1].plot(nt(eng_time,epoch), eng_xcal[:,1], '.')
        #axs[3].plot(time1, scan_mode, '.')

        #axs[0].set_title('Existing temp data')
        axs[1].set_title('Xcal position')
        #axs[3].set_title('Scan mode')
    #plt.show()

    #axs[3].imshow(ifgs, extent=[time[0], time[-1], 0, 511], vmin=vmin, vmax=vmax,
    #        aspect='auto', interpolation='none')

    ifgs = sci_mode1['ifg_data/ifg'][inds1].astype(float)
    med = np.median(ifgs, axis=1)
    ifgs1 = ifgs.T
    ifgs1 -=  med
   
    #ifgs1[:,xcal1!=1] = np.nan


    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(12, 8))
    axs = axes.flatten()

    d = np.arange(512)
    xx,yy = np.meshgrid(time1.value, d)
    im = axs[1].pcolormesh(xx, yy, ifgs1, vmin=vmin, vmax=vmax, cmap='RdBu_r')
    fig.colorbar(im, ax=axs[1])

    ifgs = sci_mode2['ifg_data/ifg'][inds2].astype(float)
    med = np.median(ifgs, axis=1)
    ifgs2 = ifgs.T
    ifgs2 -=  med

    #ifgs2[:,xcal2!=1] = np.nan

    xx,yy = np.meshgrid(time2.value, d)
    im = axs[0].pcolormesh(xx, yy, ifgs2, vmin=vmin, vmax=vmax, cmap='RdBu_r')
    fig.colorbar(im, ax=axs[0])
    plt.subplots_adjust(wspace=0, hspace=0)


    #for i in range(352, 358):
    '''
    i = 359
    plt.figure()
    plt.close('all')
    for i in range(355, 365):
        plt.plot(time1, ifgs1[i], label=f'Pos. {i}')
    plt.xlabel('Time')
    plt.ylabel('Signal')
    plt.legend(loc='best')
    plt.show()
    plt.title('IFGs1 times')
    plt.savefig('ifg_at_359.png', bbox_inches='tight')
    plt.close()
    #plt.figure()
    #for i in range(350, 360):
    #    plt.plot(time2, -ifgs2[i], label=f'Pos. {i}')
    #plt.legend(loc='best')
    #plt.title('IFGs2 times')
    #plt.show()
    '''

    return time1, ifgs1, time2, -ifgs2

# cubic interpolation clearly fits to noise fluctuations
kind = 'cubic'
kind = 'linear'


# https://docs.astropy.org/en/latest/api/astropy.visualization.time_support.html

# FIRAS mission periods, MJD in ns
t_launch = (Time('1989:322:14:35').mjd*u.day).to(u.s)
t_apco = (Time('1989:325:11:18').mjd*u.day).to(u.s) # Aperture cover ejection
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

time_periods = np.array([t_launch.value, t_apco.value, t_00.value, t_01.value, t_02.value, t_03.value, t_04.value, t_05.value, t_06.value, t_07.value, t_08.value, t_09.value, t_10.value])*u.ns
time_ranges = np.array([t_launch.value, t_apco.value, t_00.value, t_01.value, t_02.value, t_03.value, t_04.value, t_05.value, t_06.value, t_07.value, t_08.value, t_09.value, t_10.value])*u.ns
period_labels = ['Launch',
                 'APCO ejection',
                 'First light',
                 'First ICAL nulling',
                 'SAA, MTM on',
                 '2.70 to 2.75 K',
                 'SAA, MTM off',
                 'Eclipse season',
                 '6K Horns',
                 '4K Horns',
                 'Sky horn cal',
                 'Final horn temps',
                 'XCAL placed']


t_gc1_start = (Time('1990:081:04:09').mjd*u.day).to(u.ns)
t_gc1_end   = (Time('1990:084:02:23').mjd*u.day).to(u.ns)

t_gc2_start = (Time('1990:254:05:19').mjd*u.day).to(u.ns)
t_gc2_end   = (Time('1990:259:08:34').mjd*u.day).to(u.ns)


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





data = h5py.File(f'{DATA_DIR}/fdq_eng.h5')
eng = h5py.File(f'{DATA_DIR}/fdq_eng_new.h5')
sdf = h5py.File(f'{DATA_DIR}/fdq_sdf_new.h5')



grts = {}
offsets = {}

# "Binary time"
# ADT Time, in units of 100ns since 1858-11-17
bin_time = data['fdq_eng']['ct_head']['time']*(100*u.ns)

reord = np.argsort(bin_time)

# Getting bad flags
stat_word_5 = eng['en_stat/stat_word_5'][()]
stat_word_9 = eng['en_stat/stat_word_9'][()]
stat_word_13 = eng['en_stat/stat_word_13'][()]
stat_word_16 = eng['en_stat/stat_word_16'][()]
lvdt_stat_a, lvdt_stat_b = eng['en_stat/lvdt_stat'][()].T

xcal_eng = eng['en_xcal/pos'][()]

filters = (stat_word_9 == 16185)

filters = filters[reord]



dn = 8
dn = 0.5
dn = 8
dn = 16
dn = 8
dn = 128
'''
# GMT
gmt_time = data['fdq_eng']['ct_head']['gmt']
times = []
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
    grts[names_en_analog_grt[i]] = grts[names_en_analog_grt[i]][reord]
    offsets[names_en_analog_grt[i]] = ind0
    ind0 += lens_grt[i]
time = bin_time[reord]








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

t0 = 4.159641058589277e+18*u.ns
t1 = 4.159645079727482e+18*u.ns
t1 = t0 + 1*u.hr

t0 = t_01
t1 = t0 + 204*u.min

t0 = t_gc2_start - 1*u.hour
t1 = t_gc2_end   + 1*u.hour

# t0 = 2.558e7*u.s + t_apco

# Big slope between 2.558 and 2.559

# t0 = 2.5613e7*u.s + t_apco
# t1 = 2.5616e7*u.s + t_apco

#t1 = t0 + 2*u.hr
#
#t0 = 2.5584e7*u.s + t_apco
#t1 = 2.5587e7*u.s + t_apco

t0 = t_09
t1 = t_10

#t0 = 2.197425e7*u.s + t_apco
#t1 = 2.197475e7*u.s + t_apco
#
t0 = 2.197e7*u.s + t_apco
t1 = 2.198e7*u.s + t_apco

t0 = 546_500*u.s + t_09
t1 = 548_000*u.s + t_09

t0 = 546_750*u.s + t_09
t1 = 547_750*u.s + t_09

t0 = 600_000*u.s + t_09
t1 = 650_000*u.s + t_09

t0 = t_10 + 686_000*u.s
t1 = t_10 + 686_800*u.s

# Good ICAL variation region
#t0 = (1.91e7 + 2000)*u.s
#t1 = (1.91e7 + 5500)*u.s

# First xcal off
t0 = 1.085e6*u.s + t_02
t1 = 4.137e9*u.s

# Second
t0 = (4.139e9 + 500_000)*u.s
t1 = (4.139e9 + 575_000)*u.s

# Third
t0 = (4.142e9 + 20_000)*u.s
t1 = (4.142e9 + 30_000)*u.s

# Fourth
t0 = 5.235e6*u.s + t_03
t1 = (4.144e9 + 700_000)*u.s

# Fifth
t0 = (4.1473e9 + 50_000)*u.s
t1 = (4.1473e9 + 100_000)*u.s
t1 = 2.68e6*u.s + t_04

# Sixth
t0 = (4.1497e9 + 85_000)*u.s
t1 = (4.1497e9 + 120_000)*u.s


# Seventh
t0 = (4.153e9 + 460_000)*u.s
t1 = (4.153e9 + 480_000)*u.s

# Eight
# Mostly calibration, but a lot of small gaps. This is the entire calibration campaign
t0 = (4.156e9 + 300_000)*u.s
t1 = (4.156e9 + 900_000)*u.s
#t0 = (4.156e9 + 400_000)*u.s
#t1 = (4.156e9 + 500_000)*u.s

# Calibration campaign
# t0 = t_09 + 8e5*u.s
# t1 = t_09 + 9e5*u.s

# Calibration campaign, ramping up the ICAL
# t0 = t_09 + 8.4e5*u.s
# t1 = t_09 + 8.6e5*u.s

# Calibration campaign, ICAL step
# t0 = t_09 + 8.41e5*u.s
# t1 = t_09 + 8.44e5*u.s
# 
# t1 = t_09 + 8.54e5*u.s


t0 = t_09
t1 = t_10







'''
# Nine
# Half xcal, half sky
t0 = (4.1576e9 + 10_000)*u.s
t1 = (4.1576e9 + 30_000)*u.s

# Half xcal, half sky
t0 = 805_500*u.s + t_10
t1 = 815_000*u.s + t_10

# Ten
# Sky data, plus a transition
t0 = (4.15805e9 +  2_000)*u.s
t1 = (4.15805e9 + 10_000)*u.s


# 11 (ish)
# Sky data, plus a transition
t0 = 2.455e6*u.s + t_10
t1 = 2.465e6*u.s + t_10

# 12 (ish)
# Sky data, plus a transition
t0 = 2.103e6*u.s + t_10
t1 = 2.110e6*u.s + t_10

# Sky data, plus a transition
t0 = 2.1045e6*u.s + t_10
t1 = 2.108e6*u.s + t_10

# What is this? It's a transition from xcal being in the sky horn to being out of it.
t0 = 805_500*u.s + t_10
t1 = 815_000*u.s + t_10
# Zooming in on the "Bad data" I don't think the flagging entirely covers this.
t0 = 808_500*u.s + t_10
t1 = 811_000*u.s + t_10

# Stable period
t0 = 800_500*u.s + t_10
t1 = 809_800*u.s + t_10
'''

'''
# Putting the xcal in. Do we get the fuzzy region here?
t0 = (4.1586e9 + 50_000)*u.s
t1 = (4.1586e9 + 70_000)*u.s
t0 = 1.853e6*u.s + t_10
t1 = 1.856e6*u.s + t_10



# I think I've found a time when ICAL changes but XCAL stays the same
t0 = (4.159e9 + 89_000)*u.s
t1 = (4.159e9 + 96_000)*u.s

# Ah, now this is a very high ICAL temp in a calibration campaign

t0 = (4.156e9 + 600_000)*u.s
t1 = (4.156e9 + 740_000)*u.s


#t0 = (4.156e9 + 705_000)*u.s
#t0 = (4.156e9 + 608_000)*u.s
#t1 = (4.156e9 + 612_000)*u.s


#t1 = (4.156e9 + 650_000)*u.s
#
#t0 = (4.156e9 + 639_500)*u.s
#
#t0 = (4.156e9 + 642_000)*u.s
#t1 = (4.156e9 + 644_500)*u.s
#
#t0 = (4.156e9 + 642_250)*u.s
#t1 = (4.156e9 + 643_250)*u.s
'''


# This seems to have the ifgs being aligned well with the temperature data.
# I am also not sure that the plateau information is necessary.

'''
t0 = (4.14e9 - 0.3e6)*u.s
t1 = (4.14e9 + 0.9e6)*u.s

t0 = (4.14e9 -75000)*u.s
t1 = (4.14e9 + 250_000)*u.s


t0 = (4.1596e9 +  4_000)*u.s
t0 = (4.1596e9 + 12_000)*u.s
t1 = (4.1596e9 + 14_000)*u.s

# Interesting, the skyhorn has a spike that seems to be picked up in the ifgs...
# The amplitud doesn't quite match though.
t0 = (4.1596e9 + 6_000)*u.s
t1 = (4.1596e9 + 8_000)*u.s


# Zooms in more closely on the transition
t0 = (4.1596e9 + 5_000)*u.s
t1 = (4.1596e9 + 5_500)*u.s





# Trying to look at dihedral mirror now;
t0 = (4.139e9 + 390_000)*u.s

# First half good for dihedral, second half good for others
t0 = (4.139e9 + 425_000)*u.s
t1 = (4.139e9 + 430_000)*u.s



'''



'''
# Mirror and collimator now

t0 = (4.156e9 + 700_000)*u.s
t1 = (4.156e9 + 702_500)*u.s

# Pretty good for the dihedral mirror
t0 = (4.139e9 + 425_000)*u.s
t1 = (4.139e9 + 426_000)*u.s


t0 = (4.139e9 + 427_500)*u.s
t0 = (4.139e9 + 429_250)*u.s
t1 = (4.139e9 + 430_250)*u.s



# Trying to look at skyhorn and refhorn
t0 = (4.1497e9 + 30_000)*u.s
t1 = (4.1497e9 + 90_000)*u.s

#t0 = (4.1497e9 + 33_000)*u.s
#t1 = (4.1497e9 + 34_500)*u.s

t0 = (4.1497e9 + 45_000)*u.s
t1 = (4.1497e9 + 50_000)*u.s
#
#t0 = (4.1497e9 + 39_000)*u.s
#t1 = (4.1497e9 + 40_500)*u.s

#t0 = (4.1497e9 + 88_500)*u.s
#t1 = (4.1497e9 + 90_000)*u.s

# Trying to look at the mirror and collimator.
t0 = (4.138e9 + 895_500)*u.s
t1 = (4.138e9 + 904_000)*u.s

t_start = 4.138e9*u.s + 3e5*u.s + 2.5e5*u.s + 1.25e5*u.s
t0 = t_start + 2.5e4*u.s
t1 = t_start + 5e4*u.s

t_start = 4.138e9*u.s
# There are three spikes in mirror and collimator here that I want to understand, and hopefully use them to match them up.
t0 = t_start + 5e5*u.s
t1 = t_start + 1e6*u.s

# First spike
t_start = 4.138e9*u.s + 6e5*u.s + 7.5e4*u.s + 2.5e4*u.s
t0 = t_start + 8.75e3*u.s
t1 = t_start + 1.0e4*u.s

# Second spike
# Also a big jump, but not quite the same shape
# t_start += 7.5e5*u.s + 1e5*u.s
# t0 = t_start + 4e4*u.s
# t1 = t_start + 6e4*u.s

# Third spike
# Not all look especially similar here
# t_start = 4.138e9*u.s
# t_start += 7.5e5*u.s
#t_start += 2.5e5*u.s + 2.5e4*u.s
#t0 = t_start + 2e3*u.s
#t1 = t_start + 5e3*u.s

print((np.diff(time).to('s').value < 30).sum())

inds = (np.gradient(time).to('s').value < 30)
print(inds)
print(t_00)
t0 = time[inds][50].to('s') - 5_000*u.s
t1 = (t0 + 5_000*u.s).to('s')

print(t0)

# Constant ICAL, varying XCAL; there's a sort of "gaussian bump"
# This is VERY USEFUL for getting the XCAL timings.
t0 = t_09 + 8.36e5*u.s
t1 = t_09 + 8.38e5*u.s
# Should be just the rise
t1 = t_09 + 8.37e5*u.s



# This avoids the spike, but still shows the ICAL transition
t0 = (4.1596e9 + 4_500)*u.s
t0 = (4.1596e9 + 4_750)*u.s
t1 = (4.1596e9 + 5_500)*u.s
'''

'''


# Heating up at the end...

t_start = 4.160625e9*u.s
t0 = t_start + 520*u.s
t1 = t0 + 150*u.s



t0 = 4149995850.804*u.s
t1 = 4150004001.0540004*u.s

t0 = 4_149_978_331.3040004*u.s
t1 = 4_150_020_737.3040004*u.s

t0 = 4.15e9*u.s - 1.5e3*u.s
t1 = 4.15e9*u.s - 0.5e3*u.s


# Beginning of mission
t0 = 4.1344e9*u.s + 1.2e4*u.s
t1 = 4.1344e9*u.s + 1.4e4*u.s

# End of mission
t0 = 4.1606e9*u.s
t1 = t0 + 7e3*u.s


# Interesting shapes in xcal, perhaps due to looking at the galactic plane

t0 = 4.1515e9*u.s + 1e3*u.s
t1 = t0 + 5e3*u.s


t0 = 4.149e9*u.s
t1 = t0 + 1e5*u.s


t0 = 4.15e9*u.s - 5e3*u.s
t1 = 4.15e9*u.s + 5e3*u.s

'''

t0 = t_09
t1 = t_10

epoch = find_epoch(t0)
time = time.to('s')
inds = (time > t0) & (time < t1)





#inds = (np.arange(len(time)) > 550_000) & (np.arange(len(time)) < 600_000)

#inds = np.ones(len(time), dtype=bool)

grt_times = {}

grt_names = ['xcal_tip', 'skyhorn', 'refhorn', 'ical', 'dihedral',
        'mirror', 'xcal_cone', 'collimator']
offset_ind = [0, 1, 2, 3, 4, 9, 14, 15]
ifg_time_offset = 36*u.s
for side in ['a', 'b']:
    for i, grt in enumerate(grt_names):
        not_ok = grts[f'{side}_lo_{grt}'] == -9999
        grts[f'{side}_lo_{grt}'][not_ok] = np.nan
        grts[f'{side}_hi_{grt}'][not_ok] = np.nan
        not_ok = grts[f'{side}_hi_{grt}'] == -9999
        grts[f'{side}_lo_{grt}'][not_ok] = np.nan
        grts[f'{side}_hi_{grt}'][not_ok] = np.nan

        grts[f'{side}_lo_{grt}'][filters] = np.nan
        grts[f'{side}_hi_{grt}'][filters] = np.nan


        # This ordering is an artifact of the time-domain multiplexing, likely wires being plugged in a non-ideal order.
        if (grt == 'collimator') or (grt == 'xcal_cone'):
            grt_times[f'{side}_lo_{grt}'] = time + 16*u.s + ifg_time_offset
            grt_times[f'{side}_hi_{grt}'] = time  + ifg_time_offset
        else:
            grt_times[f'{side}_lo_{grt}'] = time  + ifg_time_offset
            grt_times[f'{side}_hi_{grt}'] = time + 16*u.s + ifg_time_offset

        grt_times[f'{side}_lo_{grt}'] += offset_ind[i]*u.s
        grt_times[f'{side}_hi_{grt}'] += offset_ind[i]*u.s

asdf

'''
for side in ['a', 'b']:
    for curr in ['lo', 'hi']:
        grt_times[f'{side}_{curr}_xcal_cone'] -= 8*u.s
        if side == 'a':
            grt_times[f'{side}_{curr}_xcal_tip'] -= 40*u.s


        # skyhorn appears to change temperature at a different rate....
        grt_times[f'{side}_{curr}_skyhorn'] -= 32*u.s
'''


for key in grt_times:
    print(key, grt_times[key][0] - time[0])


# plotting all temperatures



#inds = np.ones_like(inds)
'''
fig, axes = plt.subplots(sharex=True, nrows=2)
xcal_tip = np.copy(grts['a_hi_xcal_tip'])
ical = np.copy(grts['a_hi_ical'])
grts['a_hi_xcal_tip'][xcal_eng[:,0] != 1] = np.nan
grts['a_hi_ical'][xcal_eng[:,0] != 1] = np.nan
axes[0].plot(time[inds].to('s'), grts['a_hi_xcal_tip'][inds])
axes[0].set_title('XCAL')
axes[1].plot(time[inds].to('s'), grts['a_hi_ical'][inds])
axes[1].set_title('ICAL')

xcal_tip[xcal_eng[:,0] == 1] = np.nan
ical[xcal_eng[:,0] == 1] = np.nan
axes[0].plot(time[inds].to('s'), xcal_tip[inds])
axes[1].plot(time[inds].to('s'), ical[inds])
plt.savefig('t_xcal_ical.png')
#plt.show()
plt.close()
#asdf
'''

not_ok = np.isfinite(grts['a_lo_xcal_tip'])

'''
fig, axes = plt.subplots(sharex=True, nrows=2, ncols=1)
axs = axes.flatten()
axs[0].plot(nt(time[inds],epoch), grts[f'a_hi_xcal_cone'][inds] - grts[f'a_hi_ical'][inds],'.')
axs[1].plot(nt(time[inds],epoch), grts[f'a_lo_xcal_cone'][inds] - grts[f'a_lo_ical'][inds],'.')
axs[0].margins(0)
axs[1].set_xlabel(f'Time since epoch {epoch+1} [s]')
plt.close()
'''

# This is essentially a way to get the peaks to match so that it is easier to see when the changes occur.
#2.2312765 2.3015332 1.5666196 1.5565823
grt_biases = {}
grt_biases['a_hi_xcal_tip'] = 0.007
grt_biases['b_hi_xcal_tip'] = 0
grt_biases['a_hi_skyhorn'] = 0.002
grt_biases['b_hi_skyhorn'] = 0.00
grt_biases['a_hi_refhorn'] = 0.0055
grt_biases['b_hi_refhorn'] = 0.0045 + 4e-3 + 2e-4
grt_biases['a_hi_ical'] = 0.021
grt_biases['b_hi_ical'] = 0.005
grt_biases['a_hi_dihedral'] = 0.26+0.05-0.0025
grt_biases['b_hi_dihedral'] = 0.54 - 0.2 -0.02
grt_biases['a_hi_mirror'] = 2.2328577
grt_biases['b_hi_mirror'] = 2.3026211
grt_biases['a_hi_xcal_cone'] = 0.005
grt_biases['b_hi_xcal_cone'] = 0.0025
grt_biases['a_hi_collimator'] = 2.289626
grt_biases['b_hi_collimator'] = 0


grt_biases['a_lo_xcal_tip'] = 0
grt_biases['b_lo_xcal_tip'] = 0
grt_biases['a_lo_skyhorn'] = 0
grt_biases['b_lo_skyhorn'] = 0
grt_biases['a_lo_refhorn'] = 5e-4
grt_biases['b_lo_refhorn'] = 2e-3
grt_biases['a_lo_ical'] = 0
grt_biases['b_lo_ical'] = 0
grt_biases['a_lo_dihedral'] = 0.3+0.001
grt_biases['b_lo_dihedral'] = 0.3
grt_biases['a_lo_mirror'] = 1.5799628
grt_biases['b_lo_mirror'] = 1.5707934
grt_biases['a_lo_xcal_cone'] = 0
grt_biases['b_lo_xcal_cone'] = 0
grt_biases['a_lo_collimator'] = 1.6264564
grt_biases['b_lo_collimator'] = 0



grt_names = ['xcal_tip', 'skyhorn', 'refhorn', 'ical', 'dihedral',
        'mirror', 'xcal_cone', 'collimator']


#inds = np.ones_like(inds)

#fig, axes = plt.subplots(sharex=True, nrows=4, ncols=4)
#axs = axes.flatten()
fig = plt.figure(figsize=(16, 12))
ax_first = []
ax_second = []
for i in range(8):
    if i == 0:
        ax_i = fig.add_subplot(4, 4, i+1)
    else:
        ax_i = fig.add_subplot(4, 4, i+1, sharex=ax_first[0])
    ax_i8 = fig.add_subplot(4, 4, i+1+8, sharex=ax_i)#, sharey=ax_i)
    ax_first.append(ax_i)
    ax_second.append(ax_i8)

axs = ax_first + ax_second



for i in range(8):
    #grts[f'a_hi_{grt_names[i]}'][xcal_eng[:,0] != 1] = np.nan
    #grts[f'b_hi_{grt_names[i]}'][xcal_eng[:,0] != 1] = np.nan
    if i == 0:
        axs[i].plot(grt_times[f'a_hi_{grt_names[i]}'][inds], grts[f'a_hi_{grt_names[i]}'][inds], label='High')
    else:
        axs[i].plot(grt_times[f'a_hi_{grt_names[i]}'][inds], grts[f'a_hi_{grt_names[i]}'][inds])
    axs[i+8].plot(grt_times[f'b_hi_{grt_names[i]}'][inds], grts[f'b_hi_{grt_names[i]}'][inds])
    if i == 0:
        axs[i].plot(grt_times[f'a_lo_{grt_names[i]}'][inds], grts[f'a_lo_{grt_names[i]}'][inds], label='Low')
        axs[i].legend()
    else:
        axs[i].plot(grt_times[f'a_lo_{grt_names[i]}'][inds], grts[f'a_lo_{grt_names[i]}'][inds])
    axs[i+8].plot(grt_times[f'b_lo_{grt_names[i]}'][inds], grts[f'b_lo_{grt_names[i]}'][inds])

    axs[i].set_title(grt_names[i] + ' A')
    axs[i+8].set_title(grt_names[i] + ' B')

plt.tight_layout()

plt.savefig('temperature_stuff.png')


fig = plt.figure(figsize=(16, 12))
ax_first = []
ax_second = []
for i in range(8):
    if i == 0:
        ax_i = fig.add_subplot(4, 4, i+1)
    else:
        ax_i = fig.add_subplot(4, 4, i+1, sharex=ax_first[0])
    ax_i8 = fig.add_subplot(4, 4, i+1+8, sharex=ax_i)#, sharey=ax_i)
    ax_first.append(ax_i)
    ax_second.append(ax_i8)

axs = ax_first + ax_second



for i in range(8):
    #grts[f'a_hi_{grt_names[i]}'][xcal_eng[:,0] != 1] = np.nan
    #grts[f'b_hi_{grt_names[i]}'][xcal_eng[:,0] != 1] = np.nan
    if i == 0:
        axs[i].plot(grt_times[f'a_hi_{grt_names[i]}'][inds], grts[f'a_hi_{grt_names[i]}'][inds], label='A')
    else:
        axs[i].plot(grt_times[f'a_hi_{grt_names[i]}'][inds], grts[f'a_hi_{grt_names[i]}'][inds])
    axs[i+8].plot(grt_times[f'a_lo_{grt_names[i]}'][inds], grts[f'a_lo_{grt_names[i]}'][inds])

    if i == 0:
        axs[i].plot(grt_times[f'b_hi_{grt_names[i]}'][inds], grts[f'b_hi_{grt_names[i]}'][inds], label='B')
        axs[i].legend()
    else:
        axs[i].plot(grt_times[f'b_hi_{grt_names[i]}'][inds], grts[f'b_hi_{grt_names[i]}'][inds])
    axs[i+8].plot(grt_times[f'b_lo_{grt_names[i]}'][inds], grts[f'b_lo_{grt_names[i]}'][inds])

    axs[i].set_title(grt_names[i] + ' hi')
    axs[i+8].set_title(grt_names[i] + ' lo')

plt.tight_layout()

plt.savefig('temperature_stuff_ab.png')


fig, axes = plt.subplots(sharex=True, nrows=2, ncols=3)
axs = axes.flatten()
axs[0].plot(nt(time[inds],epoch), stat_word_5[inds], '.')
axs[0].set_title('statword5')
axs[1].plot(nt(time[inds],epoch), stat_word_9[inds], '.')
axs[1].set_title('statword9')
axs[2].plot(nt(time[inds],epoch), stat_word_13[inds], '.')
axs[2].set_title('statword13')
axs[3].plot(nt(time[inds],epoch), stat_word_16[inds], '.')
axs[3].set_title('statword16')
axs[4].plot(nt(time[inds],epoch), lvdt_stat_a[inds], '.')
axs[4].set_title('lvdta')
axs[5].plot(nt(time[inds],epoch), lvdt_stat_b[inds], '.')
axs[5].set_title('lvdtb')
axs[0].margins(0)
plt.savefig('stat_words.png', bbox_inches='tight')
plt.close()


fig, axes = plt.subplots(sharex=True, nrows=2, ncols=1)
axs = axes.flatten()
axs[0].plot(nt(time[inds],epoch), lvdt_stat_b[inds], '.')
axs[1].plot(nt(time[inds],epoch), stat_word_9[inds], '.')
axs[0].margins(0)
plt.suptitle('Important flags?')
#plt.savefig('eng_data.png', bbox_inches='tight')
plt.close()

ll = sdf['fdq_sdf_ll']
rh = sdf['fdq_sdf_rh']
lh = sdf['fdq_sdf_lh']
rl = sdf['fdq_sdf_rl']
t0 = t0.to('s')
t1 = t1.to('s')
time_ifg, ifg_values3, time_ifg2, ifg_values4 = plot_sdf(lh, rl, t0, t1, vmin=-100, vmax=100,
        eng_time=time[inds], 
        eng_data1=stat_word_9[inds],
        eng_data2=stat_word_13[inds],
        eng_data3=not_ok[inds],
        #eng_data3=lvdt_stat_b[inds],
        eng_xcal=xcal_eng[inds])
ifg_values3 = ifg_values3[356]
ifg_values4 = ifg_values4[356]

plt.savefig('ifg_with_eng_lhrl.png')
plt.close('all')

time_ifg, ifg_values, time_ifg2, ifg_values2 = plot_sdf(ll, rh, t0, t1, vmin=-100, vmax=100,
        eng_time=time[inds], 
        eng_data1=stat_word_9[inds],
        eng_data2=stat_word_13[inds],
        eng_data3=not_ok[inds],
        #eng_data3=lvdt_stat_b[inds],
        eng_xcal=xcal_eng[inds])

plt.savefig('ifg_with_eng_llrh.png')
plt.close('all')

ifg_values = ifg_values[359]
ifg_values2 = ifg_values2[356]

#plt.show()
#plt.close()


stat_word_5 = eng['en_stat/stat_word_5'][()]
stat_word_9 = eng['en_stat/stat_word_9'][()]
stat_word_13 = eng['en_stat/stat_word_13'][()]
stat_word_16 = eng['en_stat/stat_word_16'][()]
lvdt_stat_a, lvdt_stat_b = eng['en_stat/lvdt_stat'][()].T

'''
plt.figure()
plt.title('ll')
plt.plot(ifg_values)
plt.figure()
plt.title('rh')
plt.plot(ifg_values2)
plt.figure()
plt.title('lh')
plt.plot(ifg_values3)
plt.figure()
plt.title('rl')
plt.plot(ifg_values4)
plt.show()
'''
for ji, j in tqdm(enumerate(np.arange(-128, 129, dn))):
    for i, grt in enumerate(['collimator', 'mirror']):
        if grt == 'collimator':
            fact = 0.5
        else:
            fact = 1
        try:
            #plt.plot(nt(grt_times[f'a_hi_{grt}'][inds] + j*u.s, epoch), grts[f'a_hi_{grt}'][inds] - grt_biases[f'a_hi_{grt}'], '.', label=f'A side, high, {grt}')
            plt.plot(nt(grt_times[f'a_hi_{grt}'][inds] + j*u.s, epoch), grts[f'a_hi_{grt}'][inds] -grts[f'a_hi_{grt}'][inds].max() , '.', label=f'A side, high, {grt}')
        except:
            pass
        try:
            #plt.plot(nt(grt_times[f'b_hi_{grt}'][inds], epoch), grts[f'b_hi_{grt}'][inds] - grt_biases[f'b_hi_{grt}'], '.', label=f'B side, high, {grt}')
            plt.plot(nt(grt_times[f'b_hi_{grt}'][inds], epoch), grts[f'b_hi_{grt}'][inds] - grts[f'b_hi_{grt}'][inds].max(), '.', label=f'B side, high, {grt}')
        except:
            pass
        try:
            #plt.plot(nt(grt_times[f'a_lo_{grt}'][inds] + j*u.s, epoch), (grts[f'a_lo_{grt}'][inds] - grt_biases[f'a_lo_{grt}'])*fact, '.', label=f'A side, low, {grt}')
            plt.plot(nt(grt_times[f'a_lo_{grt}'][inds] + j*u.s, epoch), (grts[f'a_lo_{grt}'][inds] - grts[f'a_lo_{grt}'][inds].max())*fact, '.', label=f'A side, low, {grt}')
        except:
            pass
        try:
            #plt.plot(nt(grt_times[f'b_lo_{grt}'][inds], epoch), grts[f'b_lo_{grt}'][inds] - grt_biases[f'b_lo_{grt}'], '.', label=f'B side, low, {grt}')
            plt.plot(nt(grt_times[f'b_lo_{grt}'][inds], epoch), grts[f'b_lo_{grt}'][inds] - grts[f'b_lo_{grt}'][inds].max(), '.', label=f'B side, low, {grt}')
        except:
            pass
    plt.legend(loc='best')
    plt.title(f'Extra {j} seconds')
    plt.xlabel(f'Time since epoch {epoch+1} [s]')
    plt.xlim([nt(t0, epoch).value, nt(t1, epoch).value])
    plt.ylim([-1, 0.5])
    plt.savefig(f'temperature_comp_mirror_plus_coll_{ji:02}.png', bbox_inches='tight', dpi=150)
    plt.close()
plt.close('all')

for ji, j in tqdm(enumerate(np.arange(-128, 129, dn))):
#for ji, j in tqdm(enumerate(np.arange(-64, 65, dn))):
    plt.figure(figsize=(12, 8))
    for i, grt in enumerate(grt_names):
        print(grt_names[i])
        try:
            #plt.plot(nt(grt_times[f'a_hi_{grt}'][inds] + j*u.s, epoch), grts[f'a_hi_{grt}'][inds] - grt_biases[f'a_hi_{grt}'], '.', label='A side, high')
            plt.plot(nt(grt_times[f'a_hi_{grt}'][inds] + j*u.s, epoch), grts[f'a_hi_{grt}'][inds] - grts[f'a_hi_{grt}'][inds].min() + i, '.', label='A side, high',
                     color='C0')
        except:
            pass
        try:
            #plt.plot(nt(grt_times[f'b_hi_{grt}'][inds], epoch), grts[f'b_hi_{grt}'][inds] - grt_biases[f'b_hi_{grt}'], '.', label='B side, high')
            plt.plot(nt(grt_times[f'b_hi_{grt}'][inds], epoch), grts[f'b_hi_{grt}'][inds] - grts[f'b_hi_{grt}'][inds].min() + i, '.', label='B side, high',
                     color='C1')
        except:
            pass
        try:
            #plt.plot(nt(grt_times[f'a_lo_{grt}'][inds] + j*u.s, epoch), grts[f'a_lo_{grt}'][inds] - grt_biases[f'a_lo_{grt}'], '.', label='A side, low')
            plt.plot(nt(grt_times[f'a_lo_{grt}'][inds] + j*u.s, epoch), grts[f'a_lo_{grt}'][inds] -grts[f'a_lo_{grt}'][inds].min() + i, '.', label='A side, low',
                     color='C2')
        except:
            pass
        try:
            #plt.plot(nt(grt_times[f'b_lo_{grt}'][inds], epoch), grts[f'b_lo_{grt}'][inds] - grt_biases[f'b_lo_{grt}'], '.', label='B side, low')
            plt.plot(nt(grt_times[f'b_lo_{grt}'][inds], epoch), grts[f'b_lo_{grt}'][inds] - grts[f'b_lo_{grt}'][inds].min() + i, '.', label='B side, low',
                     color='C3')
        except:
            pass
        if j == 0:
            print(grt)
            print(grts[f'a_hi_{grt}'][inds][0], grts[f'b_hi_{grt}'][inds][0], grts[f'a_lo_{grt}'][inds][0], grts[f'b_lo_{grt}'][inds][0])
        '''
        if 'xcal' in grt:
            plt.plot(time_ifg, ifg_values/1400/4 + 2.75, label='LL at 359 (scaled)')
            #plt.plot(time_ifg2, (ifg_values2-200)/400/3.8 + 2.75, label='RH at 355 (scaled)')
            #plt.plot(time_ifg2, (ifg_values4)/1000 + 2.75, label='RL at 355 (scaled)')
        if grt == 'ical':
            ifg_v = -ifg_values
            ifg_v2 = -ifg_values2
            ifg_v4 = -ifg_values4
            ifg_v /= ifg_v.max()
            ifg_v2 /= ifg_v2.max()
            ifg_v4 /= ifg_v4.max()

            ifg_v  -= ifg_v.min()
            ifg_v2 -= ifg_v2.min()
            ifg_v4 -= ifg_v4.min()
            # Maybe scale here?
            ifg_v *= 1.15
            ifg_v2 *= 1.15
            ifg_v4 *= 1.15

            ifg_v  += 3.48
            ifg_v2 += 3.45
            ifg_v4 += 3.48

            plt.plot(time_ifg,  ifg_v, label='LL at 359 (scaled)')
            plt.plot(time_ifg2, ifg_v2, label='RH at 355 (scaled)')
            #plt.plot(time_ifg2, ifg_v4, label='RL at 355 (scaled)')
        '''
        #plt.legend(loc='best')
        plt.title(f'Extra {j} seconds')
        plt.xlabel(f'Time since epoch {epoch+1} [s]')
        plt.xlim([nt(t0, epoch).value, nt(t1, epoch).value])
        #plt.ylim([-1, 0.5])
        #plt.savefig(f'temperature_comp_{grt}_{ji:02}.png', bbox_inches='tight', dpi=150)
        #plt.close()
    plt.yticks(np.arange(len(grt_names)), grt_names)
    plt.savefig(f'temperature_comp_all_{ji:02}.png', bbox_inches='tight', dpi=150)
    plt.close()
plt.close('all')


for i, grt in enumerate(grt_names):
    if 'xcal' not in grt:
        continue
    plt.plot(nt(grt_times[f'a_lo_{grt}'][inds],epoch), grts[f'a_lo_{grt}'][inds],  'o-', label=f'{grt} Low current reading a')
    plt.plot(nt(grt_times[f'a_hi_{grt}'][inds] + j*u.s,epoch), grts[f'a_hi_{grt}'][inds] - grt_biases[f'a_hi_{grt}'], 'o-', label=f'{grt} High current reading a')

    plt.plot(nt(grt_times[f'b_lo_{grt}'][inds],epoch), grts[f'b_lo_{grt}'][inds],  'o-', label=f'{grt} Low current reading b')
    plt.plot(nt(grt_times[f'b_hi_{grt}'][inds] + j*u.s,epoch), grts[f'b_hi_{grt}'][inds]- grt_biases[f'b_hi_{grt}'], 'o-', label=f'{grt} High current reading b')
plt.legend(loc='best')
#plt.show()
for i, grt in enumerate(grt_names):
    if 'ical' not in grt:
        continue
    plt.plot(nt(grt_times[f'a_lo_{grt}'][inds],epoch), grts[f'a_lo_{grt}'][inds],  'o-', label=f'{grt} Low current reading a')
    plt.plot(nt(grt_times[f'a_hi_{grt}'][inds] + j*u.s,epoch), grts[f'a_hi_{grt}'][inds] - grt_biases[f'a_hi_{grt}'],  'o-', label=f'{grt} High current reading a')

    plt.plot(nt(grt_times[f'b_lo_{grt}'][inds],epoch), grts[f'b_lo_{grt}'][inds],  'o-', label=f'{grt} Low current reading b')
    plt.plot(nt(grt_times[f'b_hi_{grt}'][inds] + j*u.s,epoch), grts[f'b_hi_{grt}'][inds]- grt_biases[f'b_hi_{grt}'],  'o-', label=f'{grt} High current reading b')
plt.legend(loc='best')
#plt.show()

dn = 8
for _, j in enumerate(np.arange(-8, 8, dn)):
    fig, axes = plt.subplots(4, 4, sharex=False, sharey=False, figsize=(12, 10))
    axs = axes.flatten()
    for i, grt in enumerate(grt_names):
        T_lo = grts[f'a_lo_{grt}'][inds]
        T_hi = grts[f'a_hi_{grt}'][inds]
        T_lo[T_lo < 0] = np.nan
        T_hi[T_hi < 0] = np.nan
        t_lo = grt_times[f'a_lo_{grt}'][inds]
        t_hi = grt_times[f'a_hi_{grt}'][inds] + j*dt

        f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
        f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)
        
        times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()),
                len(t_lo))
        inds_ = (f1(times) > 0)
        axs[2*i].plot(f1(times)[inds_], f2(times)[inds_] - f1(times)[inds_], '.', ms=1)

        T_lo = grts[f'b_lo_{grt}'][inds]
        T_hi = grts[f'b_hi_{grt}'][inds]
        T_lo[T_lo < 0] = np.nan
        T_hi[T_hi < 0] = np.nan
        t_lo = grt_times[f'b_lo_{grt}'][inds]
        t_hi = grt_times[f'b_hi_{grt}'][inds] + j*dt

        f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
        f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)
        
        times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()),
                len(t_lo))
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

for _ in range(len(time_ranges)-1):
    t0 = time_ranges[_]
    t1 = time_ranges[_+1]
    inds = (time > t0) & (time < t1)
    fig, axes = plt.subplots(4, 4, sharex=False, sharey=False, figsize=(12, 10))
    axs = axes.flatten()
    for i, grt in enumerate(grt_names):
        T_lo = grts[f'a_lo_{grt}'][inds]
        T_hi = grts[f'a_hi_{grt}'][inds]
        T_lo[T_lo < 0] = np.nan
        T_hi[T_hi < 0] = np.nan
        t_lo = grt_times[f'a_lo_{grt}'][inds]
        t_hi = grt_times[f'a_hi_{grt}'][inds] + j*dt

        f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
        f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)
        
        times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()),
                len(t_lo))
        inds_ = (f1(times) > 0)
        axs[2*i].plot(f1(times)[inds_], f2(times)[inds_] - f1(times)[inds_], '.', ms=1)

        T_lo = grts[f'b_lo_{grt}'][inds]
        T_hi = grts[f'b_hi_{grt}'][inds]
        T_lo[T_lo < 0] = np.nan
        T_hi[T_hi < 0] = np.nan
        t_lo = grt_times[f'b_lo_{grt}'][inds]
        t_hi = grt_times[f'b_hi_{grt}'][inds] + j*dt

        f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
        f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)
        
        times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()),
                len(t_lo))
        inds_ = (f1(times) > 0)
        axs[2*i+1].plot(f1(times)[inds_], f2(times)[inds_] - f1(times)[inds_], '.', ms=1)
        axs[2*i].set_title(f"{grt}, a")
        axs[2*i+1].set_title(f"{grt}, b")
    fig.supylabel(r'$T_\mathrm{high} - T_\mathrm{low}$')
    fig.supxlabel(r'$T_\mathrm{low}$')
    plt.suptitle(f'Time period {_}')
    plt.tight_layout()
    plt.savefig(f'offsets_diff_timerange_{_:02}.png')
    plt.close()


'''
dn = 8
for _, j in enumerate(np.arange(-8, 8, dn)):
    fig, axes = plt.subplots(4, 4, sharex=False, sharey=False, figsize=(12, 10))
    axs = axes.flatten()
    for i, grt in enumerate(grt_names):
        T_lo = grts[f'a_lo_{grt}'][inds]
        T_hi = grts[f'a_hi_{grt}'][inds]
        T_lo[T_lo < 0] = np.nan
        T_hi[T_hi < 0] = np.nan
        t_lo = grt_times[f'a_lo_{grt}'][inds]
        t_hi = grt_times[f'a_hi_{grt}'][inds] + j*dt

        f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
        f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)
        
        times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()),
                100*len(t_lo))
        inds_ = (f1(times) > 0)
        axs[2*i].plot(times[inds_], f1(times)[inds_], '.', ms=1)
        axs[2*i].plot(times[inds_], f2(times)[inds_], '.', ms=1)

        T_lo = grts[f'b_lo_{grt}'][inds]
        T_hi = grts[f'b_hi_{grt}'][inds]
        T_lo[T_lo < 0] = np.nan
        T_hi[T_hi < 0] = np.nan
        t_lo = grt_times[f'b_lo_{grt}'][inds]
        t_hi = grt_times[f'b_hi_{grt}'][inds] + j*dt

        f1 = interp1d(t_lo, T_lo, fill_value='extrapolate', kind=kind)
        f2 = interp1d(t_hi, T_hi, fill_value='extrapolate', kind=kind)
        
        times = np.linspace(min(t_lo.min(), t_hi.min()), max(t_hi.max(), t_lo.max()),
                100*len(t_lo))
        inds_ = (f1(times) > 0)
        axs[2*i+1].plot(times[inds_], f1(times)[inds_], '.', ms=1)
        axs[2*i+1].plot(times[inds_], f2(times)[inds_], '.', ms=1)
        axs[2*i].set_title(f"{grt}, a")
        axs[2*i+1].set_title(f"{grt}, b")
    fig.supylabel(r'$T$')
    fig.supxlabel(r'time')
    plt.suptitle(f'{j} second offset')
    plt.tight_layout()
    plt.savefig(f'offsets_temps_{_:03}.png')
    plt.close()

'''
