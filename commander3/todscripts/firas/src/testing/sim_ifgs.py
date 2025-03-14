import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('..')
import my_utils as mu
import h5py
from astropy.io import fits

from utils.config import gen_nyquistl

def convert_gain(gain_array):
    conv = {0: 1, 1: 3, 2: 10, 3: 30, 4: 100, 5: 300, 6: 1000, 7: 3000}
    gain_array = gain_array.astype(float)
    for i in range(8):
        gain_array[gain_array.astype(int) == i] = conv[i]
    gain_array[(gain_array < 0)] = np.nan

    return gain_array

sdf = h5py.File('../engineering_timing/fdq_sdf_new.h5')
eng = h5py.File('../engineering_timing/fdq_eng_new.h5')


channels = {'ll':3, 'lh':2, 'rl':1, 'rh':0}
scan_modes = {'ss':0, 'sf': 1, 'ls': 2, 'lf':3, 'fs': 4, 'fl':5}


ch = 'll'
sm = 'ss'

channel = channels[ch] 
scan_mode = scan_modes[sm]


sci_times = sdf[f'fdq_sdf_{ch}/ct_head/time'][()]
ifgs = sdf[f'fdq_sdf_{ch}/ifg_data/ifg'][()]
mtm_speeds = sdf[f'fdq_sdf_{ch}/sci_head/mtm_speed'][()]
mtm_lengths = sdf[f'fdq_sdf_{ch}/sci_head/mtm_length'][()]
gains = convert_gain(sdf[f'fdq_sdf_{ch}/sci_head/gain'][()])
adds_per_groups = sdf[f'fdq_sdf_{ch}/sci_head/sc_head9'][()]
xcal_poss = sdf[f'fdq_sdf_{ch}/dq_data/xcal_pos'][()]
sweepss = sdf[f'fdq_sdf_{ch}/sci_head/sc_head11'][()]



eng_times = sdf[f'fdq_sdf_{ch}/dq_data/eng_time'][()]


pub_model = fits.open(f'FIRAS_CALIBRATION_MODEL_{ch.upper()}{sm.upper()}.FITS')

apod = pub_model[1].data['APODIZAT'][0]
apod = np.ones_like(apod)
etf = pub_model[1].data['RELEX_GA'][0] + 1j*pub_model[1].data['IELEX_GA'][0]
otf = pub_model[1].data['RTRANSFE'][0] + 1j*pub_model[1].data['ITRANSFE'][0]

S0 = pub_model[1].data["DC_RESPO"][0]
tau = pub_model[1].data["TIME_CON"][0]
Jo = pub_model[1].data["BOLPARM8"][0]
Jg = pub_model[1].data["BOLPARM9"][0]
Tbol = pub_model[1].data["BOLOM_B2"][0]
T0 = pub_model[1].data["BOLPARM2"][0]
R0 = pub_model[1].data["BOLPARM_"][0]
rho = pub_model[1].data["BOLPARM5"][0]
G1 = pub_model[1].data["BOLPARM3"][0]
beta = pub_model[1].data["BOLPARM4"][0]
C3 = pub_model[1].data["BOLPARM7"][0]
C1 = pub_model[1].data["BOLPARM6"][0]

ical_emiss = pub_model[1].data['RICAL'][0] + 1j*pub_model[1].data['IICAL'][0]


refhorn_emiss = pub_model[1].data['RREFHORN'][0] + 1j*pub_model[1].data['IREFHORN'][0]
skyhorn_emiss = pub_model[1].data['RSKYHORN'][0] + 1j*pub_model[1].data['ISKYHORN'][0]
dihedral_emiss = pub_model[1].data['RDIHEDRA'][0] + 1j*pub_model[1].data['IDIHEDRA'][0]
bolom_emiss = pub_model[1].data['RBOLOMET'][0] + 1j*pub_model[1].data['IBOLOMET'][0]

ind = np.arange(100_000,100_050)
ind = np.arange(54_825,54_875)
ind = np.arange(529_025,529_100)
#ind = np.arange(529_000,529_100)

#ind = np.arange(54_800, 55_000)

ind = np.arange(529_000,529_500)
#ind = np.arange(1_000,530_000)

#ind = np.arange(114_400,116_400)
#ind = np.arange(170_001,190_000)
#ind = np.arange(500_000,530_000)
#ind = np.arange(186_250,187_500)
#ind = np.arange(len(mtm_speeds))

channel = 3
mtm_speed = mtm_speeds[ind]
mtm_length = mtm_lengths[ind]
sweeps = sweepss[ind]
gain = gains[ind]
xcal_pos = xcal_poss[ind]
adds_per_group = adds_per_groups[ind]
eng_time = eng_times[ind]
sci_time = sci_times[ind]
ifg = ifgs[ind].astype('float32')








#print('xcal position 1 means in sky horn, 2 means out of horn, 3 is moving, 0 is error.')
#print(np.unique(xcal_pos))
#asdf
#bla = (mtm_speed == 0) & (mtm_length == 0) & ((sweeps == 16) | (sweeps == 4)) & (adds_per_group > 0) & (adds_per_group < 13)
bla = ((sweeps == 16) | (sweeps == 4)) & (adds_per_group > 0) & (adds_per_group < 13) & (xcal_pos == 1)
bla = bla & (mtm_length == 0) & (mtm_length == 0)
mtm_speed = mtm_speed[bla]
mtm_length = mtm_length[bla]
sweeps = sweeps[bla]
gain = gain[bla]
xcal_pos = xcal_pos[bla]
adds_per_group = adds_per_group[bla]
eng_time = eng_time[bla]
sci_time = sci_time[bla]
ifg = ifg[bla]

eng_time_array = eng['ct_head/time'][()]
bol_cmd_biass = eng['en_stat/bol_cmd_bias'][()].astype(int)
bol_volts = eng['en_analog/group1/bol_volt'][()]

bol_cmd_biass[bol_cmd_biass < 0] += 256
bol_cmd_biass = np.double(bol_cmd_biass) / 25.5


icals = eng['en_analog/grt/b_lo_ical'][()]
xcals = eng['en_analog/grt/b_lo_xcal_cone'][()]*0.999
skyhorn = eng['en_analog/grt/b_lo_skyhorn'][()]
refhorn = eng['en_analog/grt/b_lo_refhorn'][()]
dihedral = eng['en_analog/grt/b_lo_dihedral'][()]
# lh, ll, rh, rl is the order here.
bol_assem_temps = eng['en_analog/grt/b_lo_bol_assem'][:,1]

stat_word_5 = eng['en_stat/stat_word_5'][()]
stat_word_9 = eng['en_stat/stat_word_9'][()]
stat_word_13 = eng['en_stat/stat_word_13'][()]
stat_word_16 = eng['en_stat/stat_word_16'][()]
lvdt_stats = eng['en_stat/lvdt_stat'][()]
lvdt_stat_a, lvdt_stat_b = lvdt_stats.T

filter_bad = mu.filter_junk(stat_word_5, stat_word_9, stat_word_13, stat_word_16, lvdt_stat_a, lvdt_stat_b)
filter_bad = np.logical_and(filter_bad, (xcals[:,0] < 3))
filter_bad = np.logical_and(filter_bad, (icals[:,0] < 3))

icals_a_lo = eng['en_analog/grt/a_lo_ical'][()]
xcals_cone_a_lo = eng['en_analog/grt/a_lo_xcal_cone'][()]
xcals_tip_a_lo = eng['en_analog/grt/a_lo_xcal_tip'][()]
skyhorn_a_lo = eng['en_analog/grt/a_lo_skyhorn'][()]
refhorn_a_lo = eng['en_analog/grt/a_lo_refhorn'][()]
dihedral_a_lo = eng['en_analog/grt/a_lo_dihedral'][()]
collimator_a_lo = eng['en_analog/grt/a_lo_collimator'][()]
mirror_a_lo = eng['en_analog/grt/a_lo_mirror'][()]

icals_a_hi = eng['en_analog/grt/a_hi_ical'][()]
xcals_cone_a_hi = eng['en_analog/grt/a_hi_xcal_cone'][()]
xcals_tip_a_hi = eng['en_analog/grt/a_hi_xcal_tip'][()]
skyhorn_a_hi = eng['en_analog/grt/a_lo_skyhorn'][()]
refhorn_a_hi = eng['en_analog/grt/a_lo_refhorn'][()]
dihedral_a_hi = eng['en_analog/grt/a_lo_dihedral'][()]
collimator_a_hi = eng['en_analog/grt/a_lo_collimator'][()]
mirror_a_hi = eng['en_analog/grt/a_lo_mirror'][()]

icals_b_lo = eng['en_analog/grt/b_lo_ical'][()]
xcals_cone_b_lo = eng['en_analog/grt/b_lo_xcal_cone'][()]
xcals_tip_b_lo = eng['en_analog/grt/b_lo_xcal_tip'][()]
skyhorn_b_lo = eng['en_analog/grt/b_lo_skyhorn'][()]
refhorn_b_lo = eng['en_analog/grt/b_lo_refhorn'][()]
dihedral_b_lo = eng['en_analog/grt/b_lo_dihedral'][()]
collimator_b_lo = eng['en_analog/grt/b_lo_collimator'][()]
mirror_b_lo = eng['en_analog/grt/b_lo_mirror'][()]

icals_b_hi = eng['en_analog/grt/b_hi_ical'][()]
xcals_cone_b_hi = eng['en_analog/grt/b_hi_xcal_cone'][()]
xcals_tip_b_hi = eng['en_analog/grt/b_hi_xcal_tip'][()]
skyhorn_b_hi = eng['en_analog/grt/b_lo_skyhorn'][()]
refhorn_b_hi = eng['en_analog/grt/b_lo_refhorn'][()]
dihedral_b_hi = eng['en_analog/grt/b_lo_dihedral'][()]
collimator_b_hi = eng['en_analog/grt/b_lo_collimator'][()]
mirror_b_hi = eng['en_analog/grt/b_lo_mirror'][()]

en_xcal = eng['en_xcal/pos'][()]




print('median subtracting')
ifg -= np.median(ifg, axis=1)[:,None]
print('Median subtracted')


#fig, axes = plt.subplots(sharex=True, nrows=2, ncols=1)
#axs = axes.flatten()
#axs[0].plot(sci_time, adds_per_group, '.')
#axs[0].set_title('Adds Per Group')
#axs[1].plot(sci_time, sweeps, '.')
#axs[1].set_title('Sweeps')

inds = (eng_time_array > sci_time.min()) & (eng_time_array < sci_time.max())



#   plt.figure()
#   plt.plot(eng_time_array[inds], filter_bad[inds])
#   plt.xlim(eng_time_array[inds].min(), eng_time_array[inds].max())
eng_inds = []
eng_time_values = np.arange(len(eng_time_array))
print(ifg[:,0].shape)
from tqdm import tqdm
for i in tqdm(range(len(eng_time))):
    eng_ind_matches = eng_time_values[eng_time_array==eng_time[i]]
    if len(eng_ind_matches) == 1:
        eng_inds.append(eng_ind_matches[0])
        if not filter_bad[eng_ind_matches[0]]: 
            ifg[i] = np.nan
    else:
        eng_inds.append(eng_inds[-1])

#st = sci_time[np.isfinite(ifg[:,0])]
#ifg = ifg[np.isfinite(ifg[:,0])]
#XX,YY = np.meshgrid(st, np.arange(512))

plt.figure()
plt.pcolormesh(ifg.T, vmin=-50, vmax=50, cmap='RdBu_r')
plt.title('Data')


#plt.show()
#eng_inds = eng_time_array == eng_time


bol_cmd_bias = bol_cmd_biass[eng_inds,channel]
bol_volt = bol_volts[eng_inds,channel]

icals = icals[eng_inds]
xcals = xcals[eng_inds]
skyhorn = skyhorn[eng_inds]
refhorn = refhorn[eng_inds]

icals_a_lo = icals_a_lo[eng_inds]
icals_a_hi = icals_a_hi[eng_inds]
icals_b_lo = icals_b_lo[eng_inds]
icals_b_hi = icals_b_hi[eng_inds]

xcals_cone_a_lo = xcals_cone_a_lo[eng_inds]
xcals_cone_a_hi = xcals_cone_a_hi[eng_inds]
xcals_cone_b_lo = xcals_cone_b_lo[eng_inds]
xcals_cone_b_hi = xcals_cone_b_hi[eng_inds]

xcals_tip_a_lo = xcals_tip_a_lo[eng_inds]
xcals_tip_a_hi = xcals_tip_a_hi[eng_inds]
xcals_tip_b_lo = xcals_tip_b_lo[eng_inds]
xcals_tip_b_hi = xcals_tip_b_hi[eng_inds]

skyhorn_a_lo = skyhorn_a_lo[eng_inds]
skyhorn_a_hi = skyhorn_a_hi[eng_inds]
skyhorn_b_lo = skyhorn_b_lo[eng_inds]
skyhorn_b_hi = skyhorn_b_hi[eng_inds]

refhorn_a_lo = refhorn_a_lo[eng_inds]
refhorn_a_hi = refhorn_a_hi[eng_inds]
refhorn_b_lo = refhorn_b_lo[eng_inds]
refhorn_b_hi = refhorn_b_hi[eng_inds]

dihedral_a_lo = dihedral_a_lo[eng_inds]
dihedral_a_hi = dihedral_a_hi[eng_inds]
dihedral_b_lo = dihedral_b_lo[eng_inds]
dihedral_b_hi = dihedral_b_hi[eng_inds]

collimator_a_lo = collimator_a_lo[eng_inds]
collimator_a_hi = collimator_a_hi[eng_inds]
collimator_b_lo = collimator_b_lo[eng_inds]
collimator_b_hi = collimator_b_hi[eng_inds]

mirror_a_lo = mirror_a_lo[eng_inds]
mirror_a_hi = mirror_a_hi[eng_inds]
mirror_b_lo = mirror_b_lo[eng_inds]
mirror_b_hi = mirror_b_hi[eng_inds]


fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

frec = 4 * (channel % 2) + scan_mode

fnyq_icm = fnyq['icm'][frec]
fnyq_hz = fnyq['hz'][frec]





fig, axes = plt.subplots(sharex=True, nrows=2, ncols=4)
axs = axes.flatten()
axs[0].plot(icals_a_lo)
axs[1].plot(xcals_cone_a_lo)
axs[2].plot(xcals_tip_a_lo)
axs[3].plot(skyhorn_a_lo)
axs[4].plot(refhorn_a_lo)
axs[5].plot(dihedral_a_lo)
axs[6].plot(collimator_a_lo)
axs[6].plot(mirror_a_lo)




plt.figure()
plt.title('Data ifg')
plt.pcolormesh(ifg)
plt.colorbar()
#plt.show()



# mtm_speed, mtm_length should all be the same
mtm_speed = mtm_speed[0]
mtm_length = mtm_length[0]

afreq_out, spec_out = mu.ifg_to_spec(ifg,
                          mtm_speed,
                          channel,
                          adds_per_group,
                          bol_cmd_bias,
                          bol_volt,
                          fnyq_icm,
                          fnyq_hz,
                          otf,
                          Jo,
                          Jg,
                          Tbol, rho,
                          R0, T0, beta, G1, C3, C1,
                          mtm_length,gain,sweeps,
                          apod)

NU_ZERO =            68.020812
DELTA_NU=            13.604162

spec_out[~np.isfinite(spec_out)] = 0

fig, axes = plt.subplots(sharex=True, nrows=2)
freq = np.arange(len(spec_out[0]))*DELTA_NU
for i in range(len(spec_out)):
    axes[0].plot(freq, spec_out[i].real, 'k', alpha=0.1)
    axes[1].plot(freq, spec_out[i].imag, 'k', alpha=0.1)
axes[1].set_xlabel('Frequency (Ghz')
plt.suptitle('Calibrated spectra (real/imag)')
axes[0].set_xlim([0, 700])


plt.figure()
plt.plot(xcals, label='XCAL')
plt.plot(icals, label='ICAL')
plt.legend(loc='best')
plt.xlabel('Time [sample]')


from astropy.modeling.models import BlackBody
from astropy import units as u

fig, axes = plt.subplots(sharex=True, nrows=3, sharey=False)
nu = np.arange(len(ical_emiss))*DELTA_NU + NU_ZERO
nu *= u.GHz
spec_th = np.zeros((len(xcals)//10, len(spec_out[0])), dtype=complex)
print(spec_th.shape)
spec_xcal = np.zeros_like(spec_th)
spec_ical = np.zeros_like(spec_th)
spec_skyh = np.zeros_like(spec_th)
spec_refh = np.zeros_like(spec_th)
spec_dih = np.zeros_like(spec_th)
for i in tqdm(range(len(spec_th))):
    bb_xcal = BlackBody(temperature=xcals[i]*u.K)
    bb_ical = BlackBody(temperature=icals[i]*u.K)
    bb_refhorn = BlackBody(temperature=refhorn[i]*u.K)
    bb_skyhorn = BlackBody(temperature=skyhorn[i]*u.K)
    bb_dihedr = BlackBody(temperature=dihedral[i]*u.K)
    bb_bolassem = BlackBody(temperature=bol_assem_temps[i]*u.K)

    if i == 0:
        axes[0].plot(nu, (bb_xcal(nu)).to('MJy/sr'), label='XCAL')
        axes[0].plot(nu, (bb_ical(nu)*ical_emiss).to('MJy/sr')/otf, label='ICAL')
        axes[0].plot(nu, (bb_skyhorn(nu)*skyhorn_emiss).to('MJy/sr')/otf, label='SKYHORN')
        axes[0].plot(nu, (bb_refhorn(nu)*refhorn_emiss).to('MJy/sr')/otf, label='REFHORN')
        axes[0].plot(nu, (bb_dihedr(nu)*dihedral_emiss).to('MJy/sr')/otf, label='MIRROR')
    # Equation (14) of the Explanatory Supplement.
    R = (bb_ical(nu)*ical_emiss \
            + bb_refhorn(nu)*refhorn_emiss \
            + bb_skyhorn(nu)*skyhorn_emiss \
            + bb_dihedr(nu)*dihedral_emiss \
            #+ bb_bolassem(nu)*bolom_emiss \
            ).to('MJy/sr')/otf
    R[~np.isfinite(R)] = 0
    # Equation (13) of the Explanatory Supplement, without phase correction.
    D = spec_out[i][5:len(otf)+5]
    S_sky = D - R.value
    theory = R +  bb_xcal(nu).to('MJy/sr')
    axes[1].plot(nu, D, color='k', alpha=0.1)
    axes[1].plot(nu, theory, alpha=0.1, color='r', zorder=5)
    axes[1].set_ylim([-25, 25])
    axes[1].set_title('Calibrated spectra')

    axes[2].plot(nu, theory.real, color='k', alpha=0.1)
    #axes[2].plot(nu, theory.imag, alpha=0.1, color='r', zorder=5)
    #axes[2].set_ylim([-500, 500])
    axes[2].set_title('Theory spectra')
    axes[2].set_ylim([-25, 25])

    spec_th[i][5:len(otf) + 5] = theory.to('MJy/sr')
    spec_xcal[i][5:len(otf) + 5] = (bb_xcal(nu)).to('MJy/sr')
    spec_ical[i][5:len(otf) + 5] = (bb_ical(nu)*ical_emiss).to('MJy/sr')/otf
    spec_skyh[i][5:len(otf) + 5] = (bb_skyhorn(nu)*skyhorn_emiss).to('MJy/sr')/otf
    spec_refh[i][5:len(otf) + 5] = (bb_refhorn(nu)*refhorn_emiss).to('MJy/sr')/otf
    spec_dih[i][5:len(otf) + 5] = (bb_dihedr(nu)*dihedral_emiss).to('MJy/sr')/otf


#axes[1].plot(nu, theory*0 + np.nan, color='r', label='Model')
axes[0].legend()
#axes[1].legend()
axes[-1].set_xlabel('Frequency, GHz')
axes[0].set_xlim([0, 700])

spec_th[~np.isfinite(spec_th)] = 0
spec_ical[~np.isfinite(spec_ical)] = 0
spec_ical[~np.isfinite(spec_ical)] = 0



ifg_th = mu.spec_to_ifg(spec_th, mtm_speed, channel, adds_per_group[:len(spec_th)],
                        bol_cmd_bias[:len(spec_th)], bol_volt[:len(spec_th)], fnyq_icm, fnyq_hz, otf, Jo, Jg,
                        Tbol, rho, R0, T0, beta, G1, C3, C1, gain[:len(spec_th)],sweeps[:len(spec_th)], apod)
ifg_th[~np.isfinite(ifg_th)] = 0


ifg_xcal = mu.spec_to_ifg(spec_xcal, mtm_speed, channel, adds_per_group[:len(spec_xcal)],
                        bol_cmd_bias[:len(spec_th)], bol_volt[:len(spec_th)], fnyq_icm, fnyq_hz, otf, Jo, Jg,
                        Tbol, rho, R0, T0, beta, G1, C3, C1, gain[:len(spec_th)],sweeps[:len(spec_th)], apod)
ifg_xcal[~np.isfinite(ifg_xcal)] = 0

ifg_ical = mu.spec_to_ifg(spec_ical, mtm_speed, channel, adds_per_group[:len(spec_ical)],
                        bol_cmd_bias[:len(spec_th)], bol_volt[:len(spec_th)], fnyq_icm, fnyq_hz, otf, Jo, Jg,
                        Tbol, rho, R0, T0, beta, G1, C3, C1, gain[:len(spec_th)],sweeps[:len(spec_th)], apod)
ifg_ical[~np.isfinite(ifg_ical)] = 0

'''
plt.figure()
plt.plot(ifg_xcal[0], label='XCAL IFG')
plt.plot(ifg_ical[0], label='ICAL IFG')
plt.plot(ifg_th[0], label='Total IFG')
plt.xlabel('Samples')
plt.legend(loc='best')



plt.figure()
plt.plot(spec_xcal[0], label='XCAL')
plt.plot(spec_ical[0], label='ICAL')
plt.plot(spec_th[0], label='Theory')
plt.xlabel('Frequency [GHz]')
plt.legend(loc='best')
#plt.show()

#plt.close('all')


plt.figure()
plt.pcolormesh(ifg_th.T, vmin=-50, vmax=50, cmap='RdBu_r')
plt.title('Theory')

plt.figure()
plt.pcolormesh(ifg_xcal.T, vmin=-400, vmax=400, cmap='RdBu_r')
plt.title('XCAL')

plt.figure()
plt.pcolormesh(ifg_ical.T, vmin=-400, vmax=400, cmap='RdBu_r')
plt.title('ICAL')


'''
plt.figure(figsize=(12, 8))
ifg_med = np.nanmedian(ifg, axis=0)
ifg_th  = np.nanmedian(ifg_th - ifg_xcal - ifg_ical, axis=0)
plt.plot(ifg_med)
plt.plot(ifg_th)
plt.ylim([-70, 50])
plt.title(r'Default, with bolassem')
plt.savefig('theory_versus_med.png')

'''

plt.figure()
plt.title('Model, real')
plt.pcolormesh(spec_th.real, vmin=-10, vmax=10)
plt.figure()
plt.title('Model, imag')
plt.pcolormesh(spec_th.imag, vmin=-10, vmax=10)

plt.figure()
plt.title('Data, real')
plt.pcolormesh(spec_out.real, vmin=-10, vmax=10)
plt.figure()
plt.pcolormesh(spec_out.imag, vmin=-10, vmax=10)
plt.title('Data, imag')


plt.figure()
plt.title('Data - Model, real')
plt.pcolormesh(spec_out.real - spec_th.real, vmin=-10, vmax=10)
plt.figure()
plt.title('Data - Model, imag')
plt.pcolormesh(spec_out.imag - spec_th.imag, vmin=-10, vmax=10)
'''


#plt.show()




# ifg_out = mu.spec_to_ifg()
