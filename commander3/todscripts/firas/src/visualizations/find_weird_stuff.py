import h5py
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# Essentially, this is trying to find housekeeping data that would indicate that
# data is unreliable.


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
print(xcal_pos)

bol_cmd_bias = data['en_stat/bol_cmd_bias'][()]
cmd_inds = np.where(np.any(np.gradient(bol_cmd_bias, axis=0) != 0, axis=1))[0]
colors = np.zeros(len(cmd_inds), dtype=int)

lengths =   np.array([2, 16, 2, 16, 2, 16, 2, 16])
positions = np.array([12, 0, 12, 0, 12, 0, 12, 0])
words = [1, 4, 8, 9, 12, 13, 16]

c_ind = 0
for word, pos, l in zip(words, positions, lengths):
    word = data[f'en_stat/stat_word_{word}'][()]
    inds = np.where(np.gradient(word & ((2**l - 1) << pos)) != 0)[0]
    cmd_inds = np.concatenate((cmd_inds, inds))
    c_ind += 1
    colors = np.concatenate((colors, np.ones(len(inds), dtype=int)*c_ind))

labels0 = [f'stat_word_{word}' for word in words]
labels0 = []


labels = ['int_ref_temp_a', 'int_ref_temp_b', 'ref_hrn_temp_a',
        'ref_hrn_temp_b', 'sky_hrn_temp_a', 'sky_hrn_temp_b', 'ext_cal_temp_a',
        'ext_cal_temp_b']
lengths = np.ones(len(labels)).astype(int)*12
positions = np.zeros(len(labels)).astype(int)
for lab, pos, l in zip(labels, positions, lengths):
    word = data[f'en_stat/{lab}'][()]
    inds = np.where(np.gradient(word & ((2**l - 1) << pos)) != 0)[0]
    cmd_inds = np.concatenate((cmd_inds, inds))
    c_ind += 1
    colors = np.concatenate((colors, np.ones(len(inds), dtype=int)*c_ind))



#labs = ['bol cmd bias'] + [f'word_{i}' for i in words]
labs = ['bol cmd bias'] + labels0 + labels

#fig, axes = plt.subplots(sharex=True, sharey='row', nrows=6, ncols=4)

time = data['ct_head/time'][()]
t_min = 4.157e16 + 5e12
t_max = 4.157e16 + 7.5e12

t_min = 4.144447749785283e+16
t_max = 4.144592924341427e+16

# Dihedral oscillation zoom
t_min = 4.1374414277253224e+16
t_max = 4.1375382536451064e+16


#    t_min = 4.15761328125e+16
#    t_max = 4.157630670362903e+16

#    t_min = 4.159533884432928e+16
#    t_max = 4.159862112663672e+16
#    
#    
#    t_min = 4.136728293892138e+16
#    t_max = 4.137155677051925e+16
#    
#    t_min = 4.136981851668978e+16
#    t_max = 4.1370466771129864e+16
#    
#t_min = time.min()
#t_max = time.max()
inds = (time > t_min) & (time < t_max)

ax_ind = 0
grad_inds = np.array([], dtype=int)
'''
fig, axes = plt.subplots(sharex=True, nrows=6, ncols=4)
axs = axes.flatten()
grad_tol = 1e-3
for side in ['a', 'b']:
    print(side)
    for grt in grts:
        if (grt == 'bol_assem') or (grt == 'cal_resistors'):
            for i in range(4):
                #for j in range(1,len(cmd_inds)):
                #    axs[ax_ind].axvline(x=time[cmd_inds[j]], color='k')
                for curr in ['lo', 'hi']:
                    t = data[f'en_analog/grt/{side}_{curr}_{grt}'][:,i]
                    t[t < 0] = np.nan
                    axs[ax_ind].semilogy(time[inds],
                            abs(np.gradient(t[inds])/t[inds]))
                    grad_inds = np.concatenate((grad_inds, 
                        np.where(abs(np.gradient(t)/t) > grad_tol)[0]))
                axs[ax_ind].set_title(f'{side}_{grt}_{i+1}')
                ax_ind += 1
        else:
            #for j in range(1, len(cmd_inds)):
            #    if (time[cmd_inds[j]] > t_min) & (time[cmd_inds[j]] < t_max):
            #        axs[ax_ind].axvline(x=time[cmd_inds[j]],
            #                color=f'C{colors[j]}')
            for curr in ['lo', 'hi']:
                t = data[f'en_analog/grt/{side}_{curr}_{grt}'][:,0]
                t[t < 0] = np.nan
                #t[np.bitwise_and(stat_word_9, 832) != 0] = np.nan
                #t[xcal_pos != 1.] = np.nan
                axs[ax_ind].semilogy(time[inds], abs(np.gradient(t[inds])/t[inds]))
                grad_inds = np.concatenate((grad_inds, 
                    np.where(abs(np.gradient(t)/t) > grad_tol)[0]))
                #axs[ax_ind].plot(t[inds])
            axs[ax_ind].set_title(f'{side}_{grt}')
            ax_ind += 1


grad_inds = np.unique(grad_inds)
'''


fig, axes = plt.subplots(sharex=True, nrows=6, ncols=4)

axs = axes.flatten()

ax_ind = 0
for i in range(1, len(labs)):
    axs[ax_ind].plot([], [], linestyle='-', color=f'C{i}', label=labs[i])
#axs[ax_ind].legend()


stat_word_9 = data[f'en_stat/stat_word_9'][()]
lvd_stat = np.zeros(len(stat_word_9))

#for i in tqdm(range(window, len(stat_word_9))):
#    lvd_stat[i-window:] = data[f'en_stat/lvdt_stat'][i-window:,0].mean()
N = 20
lvd_stat = np.convolve(data[f'en_stat/lvdt_stat'][:,0], np.ones(N)/N,
        mode='same')

for side in ['a', 'b']:
    print(side)
    for grt in grts:
        if (grt == 'bol_assem') or (grt == 'cal_resistors'):
            for i in range(4):
                #for j in range(1,len(cmd_inds)):
                #    axs[ax_ind].axvline(x=time[cmd_inds[j]], color='k')
                for curr in ['lo', 'hi']:
                    t = data[f'en_analog/grt/{side}_{curr}_{grt}'][:,i]
                    t[t < 0] = np.nan
                    #t[lvd_stat > 1/N] = np.nan
                    t[grad_inds] = np.nan
                    #t[xcal_pos != 1.] = np.nan
                    #t[np.bitwise_and(stat_word_9, 832) != 0] = np.nan
                    axs[ax_ind].plot(time[inds], t[inds])
                    #axs[ax_ind].plot(t[inds])
                axs[ax_ind].set_title(f'{side}_{grt}_{i+1}')
                ax_ind += 1
        else:
            #for j in range(1, len(cmd_inds)):
            #    if (time[cmd_inds[j]] > t_min) & (time[cmd_inds[j]] < t_max):
            #        axs[ax_ind].axvline(x=time[cmd_inds[j]],
            #                color=f'C{colors[j]}')
            for curr in ['lo', 'hi']:
                t = data[f'en_analog/grt/{side}_{curr}_{grt}'][:,0]
                t[t < 0] = np.nan
                #t[lvd_stat > 1/N] = np.nan
                t[grad_inds] = np.nan
                #t[np.bitwise_and(stat_word_9, 832) != 0] = np.nan
                #t[xcal_pos != 1.] = np.nan
                axs[ax_ind].plot(time[inds], t[inds])
                

            #if grt == 'dihedral':
            #    axs[ax_ind].axhline(5.5, color='k')
            axs[ax_ind].set_title(f'{side}_{grt}')
            ax_ind += 1

#un = np.unique(data[f'en_stat/stat_word_9'][()])
#print(un)
#b = [bin(u) for u in un]
#print(b)
axs[12].cla()
axs[12].plot(time[inds], lvd_stat[inds])
axs[12].plot(time[inds], data[f'en_stat/lvdt_stat'][inds,0], '.')
axs[12].set_title('lvdt stat a')
axs[19].cla()
axs[19].plot(time[inds], data[f'en_stat/lvdt_stat'][inds,1], '.')
axs[19].set_title('lvdt stat b')
'''
# Very little correlation obvious here.
axs[12].cla()
axs[12].plot(time[inds], data[f'en_tail/lmac_analog_temp'][inds])
axs[12].set_title('lmac analog temp')
axs[19].cla()
axs[19].plot(time[inds], data[f'en_tail/lmac_digital_temp'][inds])
axs[19].set_title('lmac digital temp')
'''

plt.xlim([t_min, t_max])



# The ideal thing would be to see what the IFGs look like around the times when
# the variation is large. Are the oscillations actually in the data??
plt.xlim([t_min, t_max])

sci_data = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5')

#fig, axes = plt.subplots(sharex=True, sharey=True, nrows=2, ncols=2)
#axs = axes.flatten()
#for i, mode in enumerate(['ll', 'lh', 'rh', 'rl']):

fig, axes = plt.subplots(sharex=True, sharey=True, nrows=1, ncols=1)
axs = [axes]
from scipy.interpolate import interp1d

f = interp1d(time, data['en_stat/stat_word_9'][()], kind='nearest', fill_value='extrapolate')

for i, mode in enumerate(['ll']):

    print(i)

    eng_time = sci_data[f'fdq_sdf_{mode}/dq_data/eng_time'][()]
    ifg      = sci_data[f'fdq_sdf_{mode}/ifg_data/ifg'][()]
    xcal_pos = sci_data[f'fdq_sdf_{mode}/dq_data/xcal_pos'][()]


        

    #print(stat_word_9)
    stat_word_9 = f(eng_time).astype(int)
    print(stat_word_9)
    print(len(stat_word_9))
    print((np.bitwise_and(stat_word_9, 832) != 0).sum())
    print((np.bitwise_and(stat_word_9, 32768) != 0).sum())
    #0b1101000000

    inds = np.where((eng_time > t_min) & (eng_time < t_max))
    
    
    ifg = ifg.astype(float)
    


    # This seems to correspond to funny bumps in the GRT temperature. A few
    # samples after this also shows variation, but it's not obviously bad in the
    # same way.

    ifg[np.bitwise_and(stat_word_9, 832) != 0] = np.nan
    ifg[np.bitwise_and(stat_word_9, 32768) == 0] = np.nan

    # xcal_pos == 1 is "Xcal in"
    # ifg[xcal_pos != 1] = np.nan
    
    ifg = ifg[inds]
    
    
    t = eng_time[inds]
    ifg_ind = np.arange(512)
    
    m = np.median(ifg, axis=1)
    axs[i].pcolormesh(t, ifg_ind, ifg.T - m, vmin=-25, vmax=25, cmap='RdBu_r')
    axs[i].set_xlim([t_min, t_max])

#plt.savefig('example_nocut.png', dpi=300)
plt.show()
