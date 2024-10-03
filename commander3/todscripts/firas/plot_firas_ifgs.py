import h5py
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd


sci = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5')
eng = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5')


'''
a_hi_bol_assem           Dataset {589069, 4}
a_hi_cal_resistors       Dataset {589069, 4}
a_hi_collimator          Dataset {589069, 1}
a_hi_dihedral            Dataset {589069, 1}
a_hi_ical                Dataset {589069, 1}
a_hi_mirror              Dataset {589069, 1}
a_hi_refhorn             Dataset {589069, 1}
a_hi_skyhorn             Dataset {589069, 1}
a_hi_xcal_cone           Dataset {589069, 1}
a_hi_xcal_tip            Dataset {589069, 1}
a_lo_bol_assem           Dataset {589069, 4}
a_lo_cal_resistors       Dataset {589069, 4}
a_lo_collimator          Dataset {589069, 1}
a_lo_dihedral            Dataset {589069, 1}
a_lo_ical                Dataset {589069, 1}
a_lo_mirror              Dataset {589069, 1}
a_lo_refhorn             Dataset {589069, 1}
a_lo_skyhorn             Dataset {589069, 1}
a_lo_xcal_cone           Dataset {589069, 1}
a_lo_xcal_tip            Dataset {589069, 1}
b_hi_bol_assem           Dataset {589069, 4}
b_hi_cal_resistors       Dataset {589069, 4}
b_hi_collimator          Dataset {589069, 1}
b_hi_dihedral            Dataset {589069, 1}
b_hi_ical                Dataset {589069, 1}
b_hi_mirror              Dataset {589069, 1}
b_hi_refhorn             Dataset {589069, 1}
b_hi_skyhorn             Dataset {589069, 1}
b_hi_xcal_cone           Dataset {589069, 1}
b_hi_xcal_tip            Dataset {589069, 1}
b_lo_bol_assem           Dataset {589069, 4}
b_lo_cal_resistors       Dataset {589069, 4}
b_lo_collimator          Dataset {589069, 1}
b_lo_dihedral            Dataset {589069, 1}
b_lo_ical                Dataset {589069, 1}
b_lo_mirror              Dataset {589069, 1}
b_lo_refhorn             Dataset {589069, 1}
b_lo_skyhorn             Dataset {589069, 1}
b_lo_xcal_cone           Dataset {589069, 1}
b_lo_xcal_tip            Dataset {589069, 1}
'''


'''
grts = eng['en_analog/grt']

time = eng['ct_head/time'][()]
gmt = eng['ct_head/gmt'][()].astype(str)
t_eng = pd.to_datetime(gmt, format='%y%j%H%M%S%f')

gmt_lh = pd.to_datetime(sci['fdq_sdf_lh/ct_head/gmt'][()].astype(str),format='%y%j%H%M%S%f')
gmt_ll = pd.to_datetime(sci['fdq_sdf_ll/ct_head/gmt'][()].astype(str),format='%y%j%H%M%S%f')
gmt_rh = pd.to_datetime(sci['fdq_sdf_rh/ct_head/gmt'][()].astype(str),format='%y%j%H%M%S%f')
gmt_rl = pd.to_datetime(sci['fdq_sdf_rl/ct_head/gmt'][()].astype(str),format='%y%j%H%M%S%f')


The science and engineering data are collected asynchronously, and the midpoint
time is used to link them together.

science record times are put into engineering records and vice versa...
'''

'''
plt.plot(t_eng)
plt.plot(gmt_lh)
plt.plot(gmt_ll)
plt.plot(gmt_rh)
plt.plot(gmt_rl)
'''


# In DQ_DATA for science data;
'''
scalar/adt ENG_TIME
! Time of associate
! engineeering record
'''


#plt.figure()
# Empty
#plt.plot(sci['fdq_sdf_lh/dq_data/eng_rec'][()], '.')

eng_times = eng['ct_head/time'][()]
sci_times = eng['en_head/sci_time/bin_time'][()]

label = ['Bad', 'XCAL', 'Sky', 'Transiting']

# Calibration?
i0 = 5_000

#i0 = 0

speeds =  ['S', 'F']
lengths = ['S', 'L']
for i in range(i0, i0+1):
    fig, axes = plt.subplots(sharex=True, nrows=4)
    for k in range(4):
        axes[0].plot([], [], f'C{k}', label=f'{label[k]}')
    axes[0].legend(loc='best')
    xcal_tip = eng['en_analog/grt/a_hi_xcal_tip'][i][0]
    ical     = eng['en_analog/grt/a_hi_ical'][i][0]
    # If you plot this, you can see when they start using the external
    # calibrator more
    xcal_pos = eng['en_xcal/pos'][i].mean()
    for j, mode in enumerate(['ll', 'lh', 'rh', 'rl']):
        axes[j].set_ylabel(mode.upper())
        inds = np.where((sci[f'fdq_sdf_{mode}/dq_data/eng_time'][()] == eng_times[i]))[0]
        print(inds)
        if len(inds) > 0:
            t = pd.to_datetime(sci[f'fdq_sdf_{mode}/ct_head/gmt'][inds[0]].astype(str),format='%y%j%H%M%S%f')
            speed = sci[f'fdq_sdf_{mode}/sci_head/mtm_speed'][inds[0]]
            length = sci[f'fdq_sdf_{mode}/sci_head/mtm_length'][inds[0]]

            gain = eng[f'chan/sci_gain'][i, j]
            speed = eng['chan/xmit_mtm_speed'][i,j]
            length = eng['chan/xmit_mtm_len'][i,j]
            xcal_pos = sci[f'fdq_sdf_{mode}/dq_data/xcal_pos'][inds[0]]
            axes[j].plot(sci[f'fdq_sdf_{mode}/ifg_data/ifg'][inds[0]], 
                    color=f'C{xcal_pos}')
    if xcal_pos == 1:
        fig.suptitle(f'{lengths[length]}{speeds[speed]}, XCAL = {xcal_tip:.3f}, ical = {ical:.3f}')
    else:
        fig.suptitle(f'{lengths[length]}{speeds[speed]}, XCAL out, ical = {ical:.3f}')



plt.show()
