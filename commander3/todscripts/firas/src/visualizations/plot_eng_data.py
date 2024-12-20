import h5py
import matplotlib.pyplot as plt
import numpy as np

# Plots the values that are used for determining plateaus in the FEC routine

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

grts = ['collimator', 'dihedral', 'ical', 'mirror', 'refhorn', 'skyhorn',
        'xcal_cone', 'xcal_tip', 'bol_assem']

data = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5')

xcal_pos = data['en_xcal/pos'][()].mean(axis=1)

bol_cmd_bias = data['en_stat/bol_cmd_bias'][()]
cmd_inds = np.where(np.any(np.gradient(bol_cmd_bias, axis=0) != 0, axis=1))[0]
colors = np.zeros(len(cmd_inds), dtype=int)

cmds = bol_cmd_bias

lengths =   np.array([2, 16, 2, 16, 2, 16, 2, 16])
positions = np.array([12, 0, 12, 0, 12, 0, 12, 0])
words = [1, 4, 8, 9, 12, 13, 16]

c_ind = 0
for word, pos, l in zip(words, positions, lengths):
    word = data[f'en_stat/stat_word_{word}'][()]
    inds = np.where(np.gradient(word & ((2**l - 1) << pos)) != 0)[0]
    cmd_inds = np.concatenate((cmd_inds, inds))
    cmds = np.hstack((cmds, word[:,np.newaxis]))
    c_ind += 1
    colors = np.concatenate((colors, np.ones(len(inds), dtype=int)*c_ind))

labels0 = [f'stat_word_{word}' for word in words]


labels1 = ['int_ref_temp_a', 'int_ref_temp_b', 'ref_hrn_temp_a',
        'ref_hrn_temp_b', 'sky_hrn_temp_a', 'sky_hrn_temp_b', 'ext_cal_temp_a',
        'ext_cal_temp_b']
lengths = np.ones(len(labels1)).astype(int)*12
positions = np.zeros(len(labels1)).astype(int)
for lab, pos, l in zip(labels1, positions, lengths):
    word = data[f'en_stat/{lab}'][()]
    #cmds = np.concatenate((cmds, word))
    cmds = np.hstack((cmds, word[:,np.newaxis]))
    inds = np.where(np.gradient(word & ((2**l - 1) << pos)) != 0)[0]
    cmd_inds = np.concatenate((cmd_inds, inds))
    c_ind += 1
    colors = np.concatenate((colors, np.ones(len(inds), dtype=int)*c_ind))


labels2 = ['dwell_stat', 'grt_addr', 'hot_spot_cmd', 'lvdt_stat',
        'micro_stat_bus', 'power_a_status', 'power_b_status']

for lab in labels2:
    word = data[f'en_stat/{lab}'][()]
    if len(word.shape) == 0:
        cmds = np.hstack((cmds, word[:,np.newaxis]))
    else:
        cmds = np.hstack((cmds, word))

labels2 = ['dwell_stat', 'dwell_stat2', 'grt_addr', 'grt_addr2', 'hot_spot_cmd',
        'hot_spot_cmd2', 'lvdt_stat', 'lvdt_stat2',
        'micro_stat_bus', 'micro_stat_bus2', 'micro_stat_bus3',
        'micro_stat_bus4', 'power_a_status', 'power_a_status2',
        'power_b_status', 'power_b_status2']
#labs = ['bol cmd bias'] + [f'word_{i}' for i in words]


        

labs = [f'bol cmd bias {i}' for i in range(4)] + labels0 + labels1 + labels2


#fig, axes = plt.subplots(sharex=True, sharey='row', nrows=6, ncols=4)
fig, axes = plt.subplots(sharex=True, nrows=7, ncols=5)

axs = axes.flatten()

for i in range(len(cmds[0])):
    axs[i].set_title(labs[i])
    axs[i].plot(cmds[:,i], '.', ms=1)


cmds = bol_cmd_bias

labels3 = ['bol_volt', 'cna_temp', 'dbx_tmp', 'hot_spot', 'ipdu_amp',
        'ipdu_volt', 'ipdu_temp', 'mtm_cal_mtr', 'mtm_pos', 'pamp_chan',
        'pamp_op', 'stat_mon_temp', 'temp_ctrl']

for lab in labels3:
    word = data[f'en_analog/group1/{lab}'][()]
    if len(word.shape) == 0:
        cmds = np.hstack((cmds, word[:,np.newaxis]))
    else:
        cmds = np.hstack((cmds, word))

labels3 = [f'bol_volt_{i}' for i in range(1,5)] + \
          [f'cna_temp_{i}' for i in range(1,5)] + \
          [f'dbx_tmp_{i}' for i in range(1,3)] + \
          [f'hot_spot_{i}' for i in range(1, 3)] +  \
          [f'ipdu_amp_{i:02}' for i in range(1, 13)] + \
          [f'ipdu_volt_{i:02}' for i in range(1, 21)] + \
          [f'ipdu_temp_{i:02}' for i in range(1, 3)] +  \
          [f'mtm_cal_mtr_{i:02}' for i in range(1,3)] + \
          [f'mtm_pos_{i:02}' for i in range(1,3)] +  \
          ['pamp_chan', 'pamp_op', 'stat_mon_temp_1', 'stat_mon_temp_2'] +  \
          [f'temp_ctrl_{i}' for i in range(1,9)]

labs = labels3

fig, axes = plt.subplots(sharex=True, nrows=8, ncols=8)

axs = axes.flatten()

for i in range(len(labs)):
    axs[i].set_title(labs[i])
    axs[i].plot(cmds[:,i+4], '.', ms=1, alpha=0.5)

plt.show()
