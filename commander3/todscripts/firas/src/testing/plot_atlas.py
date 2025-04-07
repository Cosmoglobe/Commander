import sys

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
sys.path.append('..')
import h5py
import my_utils as mu
from astropy import units as u
from astropy.io import fits
from astropy.modeling.models import BlackBody

sys.path.append('..')
import h5py
import my_utils as mu
from astropy.io import fits
from utils.config import gen_nyquistl

sys.path.append('../utils')
import dataio2 as dataio

dt_float = '<f4'
dt_fex_dtrf = np.dtype([('trans', dt_float, (128,))])

fn = os.path.join('../../reference', 'FEX_DTRF.DAT')
dtrf = np.fromfile(fn, dtype=dt_fex_dtrf)
dataio.fix_floats(dtrf)

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

if ch == 'rl':
    vlim = 25
else:
    vlim = 50

dtrf_rec = 4*(channels[ch] % 2) + scan_modes[sm]
dtrf = dtrf[dtrf_rec]['trans']

if sm == 'lf':
    SM = 'FA'
else:
    SM = sm.upper()

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


pub_model = fits.open(f'FIRAS_CALIBRATION_MODEL_{ch.upper()}{SM}.FITS')
NU_ZERO =            pub_model[0].header['NU_ZERO']
DELTA_NU=            pub_model[0].header['DELTA_NU']
NUM_FREQ=            pub_model[0].header['NUM_FREQ']
offset = np.round(NU_ZERO/DELTA_NU).astype(int)

nu = np.arange(NUM_FREQ)*DELTA_NU + NU_ZERO
nu *= u.GHz

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


bolom_emiss = pub_model[1].data['RBOLOMET'][0] + 1j*pub_model[1].data['IBOLOMET'][0]
refhorn_emiss = pub_model[1].data['RREFHORN'][0] + 1j*pub_model[1].data['IREFHORN'][0]
skyhorn_emiss = pub_model[1].data['RSKYHORN'][0] + 1j*pub_model[1].data['ISKYHORN'][0]
dihedral_emiss = pub_model[1].data['RDIHEDRA'][0] + 1j*pub_model[1].data['IDIHEDRA'][0]
struct_emiss = pub_model[1].data['RSTRUCTU'][0] + 1j*pub_model[1].data['ISTRUCTU'][0]


otf = otf[:NUM_FREQ]
etf = etf[:NUM_FREQ]
bolom_emiss = bolom_emiss[:NUM_FREQ]
refhorn_emiss = refhorn_emiss[:NUM_FREQ]
skyhorn_emiss = skyhorn_emiss[:NUM_FREQ]
dihedral_emiss = dihedral_emiss[:NUM_FREQ]
ical_emiss = ical_emiss[:NUM_FREQ]
struct_emiss = struct_emiss[:NUM_FREQ]


ind = (sci_times >= 41574611466410000) & (sci_times <=  41574863373910000)


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
xcals = eng['en_analog/grt/b_lo_xcal_cone'][()]
skyhorn = eng['en_analog/grt/b_lo_skyhorn'][()]
refhorn = eng['en_analog/grt/b_lo_refhorn'][()]
dihedral = eng['en_analog/grt/b_lo_dihedral'][()]
# lh, ll, rh, rl is the order here.
bol_assem_temps = eng['en_analog/grt/b_lo_bol_assem'][:,1]
struct = eng['en_analog/grt/b_lo_mirror'][()]

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




inds = np.arange(512)
meds = np.median(ifg, axis=1)
ifg -= np.median(ifg, axis=1)[:,None]


inds = (eng_time_array > sci_time.min()) & (eng_time_array < sci_time.max())
eng_inds = []
eng_time_values = np.arange(len(eng_time_array))

from tqdm import tqdm

for i in tqdm(range(len(eng_time))):
    eng_ind_matches = eng_time_values[eng_time_array==eng_time[i]]
    if len(eng_ind_matches) == 1:
        eng_inds.append(eng_ind_matches[0])
        if not filter_bad[eng_ind_matches[0]]: 
            ifg[i] = np.nan
    else:
        eng_inds.append(eng_inds[-1])

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
plt.savefig('all_temps.png', bbox_inches='tight')




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


spec_out[~np.isfinite(spec_out)] = 0

plt.close('all')


spec_mbb = np.zeros((5, len(spec_out[0])), dtype=complex)
spec_mbb[0][offset:offset+NUM_FREQ] = BlackBody(temperature=20*u.K)(nu).value*(nu/(353*u.GHz)).to('').value**1.6
spec_mbb[1][offset:offset+NUM_FREQ] = BlackBody(temperature=10*u.K)(nu).value
spec_mbb[2][offset:offset+NUM_FREQ] = BlackBody(temperature=2.73*u.K)(nu).value
spec_mbb[3][offset:offset+NUM_FREQ] = (nu/30*u.GHz).value**-3

x = np.zeros(len(nu))
x[nu.value==299.291566] = 1
spec_mbb[4][offset:offset+NUM_FREQ] = x


labels = ['MBB, 20 K, 1.6', 'BB, 10 K', 'BB 2.73 K', 'Synch, -3', r'$\delta(\nu-300\,\mathrm{GHz})$']

for i in range(len(spec_mbb)):
    spec_mbb[i] /= max(1e-12, max(abs(spec_mbb[i])))

ifg_mbb = mu.spec_to_ifg(spec_mbb, mtm_speed, channel, adds_per_group[:5],
                        bol_cmd_bias[:5], bol_volt[:5], fnyq_icm, fnyq_hz, otf, Jo, Jg,
                        Tbol, rho, R0, T0, beta, G1, C3, C1, gain[:5],sweeps[:5], apod)
#ifg_mbb[~np.isfinite(ifg_mbb)] = 0
fig, axes = plt.subplots(nrows=5, ncols=2, sharex='col', figsize=(12,12))
for i in range(5):
    axes[i,0].plot(nu, spec_mbb[i][offset:offset+NUM_FREQ])
    axes[i,1].plot(ifg_mbb[i])
    axes[i,0].set_ylabel(labels[i])
plt.savefig('ifg_spec_atlas.png', bbox_inches='tight')
plt.close()

spec_mbb = np.zeros((5, len(spec_out[0])), dtype=complex)
spec_mbb[0][offset:offset+NUM_FREQ] = (nu/(300*u.GHz)).value**-1
spec_mbb[1][offset:offset+NUM_FREQ] = (nu/(300*u.GHz)).value**0
spec_mbb[2][offset:offset+NUM_FREQ] = (nu/(300*u.GHz)).value**1


ifg_mbb = mu.spec_to_ifg(spec_mbb, mtm_speed, channel, adds_per_group[:5],
                        bol_cmd_bias[:5], bol_volt[:5], fnyq_icm, fnyq_hz, otf, Jo, Jg,
                        Tbol, rho, R0, T0, beta, G1, C3, C1, gain[:5],sweeps[:5], apod)
fig, axes = plt.subplots(ncols=2)
for i in range(3):
    axes[0].plot(nu, spec_mbb[i][offset:offset+NUM_FREQ])
    axes[1].plot(ifg_mbb[i]/ifg_mbb[i].max())
plt.savefig('ifg_powlaw_atlas.png', bbox_inches='tight')

fig, axes = plt.subplots(sharex=True, nrows=2, sharey=True, figsize=(12, 6))
spec_th = np.zeros((5, len(spec_out[0])), dtype=complex)
spec_xcal = np.zeros_like(spec_th)
spec_ical = np.zeros_like(spec_th)
spec_skyh = np.zeros_like(spec_th)
spec_refh = np.zeros_like(spec_th)
spec_dih = np.zeros_like(spec_th)
spec_bol = np.zeros_like(spec_th)
spec_struct = np.zeros_like(spec_th)
print('xcals[i], icals[i], refhorn[i], skyhorn[i], bol_assem_temps[i], struct[i]')
for i in tqdm(range(len(spec_th))):
    if i == 0:
        print(xcals[i], icals[i], refhorn[i], skyhorn[i], bol_assem_temps[i], struct[i])
    bb_xcal = BlackBody(temperature=xcals[i]*u.K)
    bb_ical = BlackBody(temperature=icals[i]*u.K)
    bb_refhorn = BlackBody(temperature=refhorn[i]*u.K)
    bb_skyhorn = BlackBody(temperature=skyhorn[i]*u.K)
    bb_dihedr = BlackBody(temperature=dihedral[i]*u.K)
    bb_bolassem = BlackBody(temperature=bol_assem_temps[i]*u.K)

    # "Structure" is "mirror and colliamator"?
    bb_struct = BlackBody(temperature=struct[i]*u.K)

    if i == 0:
        #axes[0].plot(nu, (bb_xcal(nu)).to('MJy/sr'), label='XCAL')
        #axes[0].plot(nu, (bb_ical(nu)*ical_emiss).to('MJy/sr')/otf, label='ICAL')
        axes[0].plot(nu, (bb_skyhorn(nu)*skyhorn_emiss).to('MJy/sr')/otf, label='SKYHORN')
        axes[0].plot(nu, (bb_refhorn(nu)*refhorn_emiss).to('MJy/sr')/otf, label='REFHORN')
        axes[0].plot(nu, (bb_dihedr(nu)*dihedral_emiss).to('MJy/sr')/otf, label='MIRROR')
        axes[0].plot(nu, (bb_bolassem(nu)*bolom_emiss).to('MJy/sr')/otf, label='BOL')
        axes[0].plot(nu, (bb_struct(nu)*struct_emiss).to('MJy/sr')/otf, label='STRUCT')
    # Equation (14) of the Explanatory Supplement.
    R = (bb_ical(nu)*ical_emiss \
            + bb_refhorn(nu)*refhorn_emiss \
            + bb_skyhorn(nu)*skyhorn_emiss \
            + bb_dihedr(nu)*dihedral_emiss \
            + bb_bolassem(nu)*bolom_emiss \
            + bb_struct(nu)*struct_emiss \
            ).to('MJy/sr')/otf
    R[~np.isfinite(R)] = 0
    # Equation (13) of the Explanatory Supplement, without phase correction.
    D = spec_out[i][5:len(otf)+5]
    S_sky = D - R.value
    theory = R +  bb_xcal(nu).to('MJy/sr')
    axes[1].plot(nu, D, color='k', alpha=0.1)
    axes[1].plot(nu, theory, alpha=0.1, color='r', zorder=5)
    axes[1].set_ylim([-vlim/2, vlim/2])
    axes[1].set_title('Calibrated spectra')

    #axes[2].plot(nu, theory.real, color='k', alpha=0.1)
    #axes[2].set_title('Theory spectra')
    #axes[2].set_ylim([-25, 25])


    spec_th[i][offset:offset+NUM_FREQ] = theory.to('MJy/sr')

    spec_xcal[i][offset:offset+NUM_FREQ] = (bb_xcal(nu)).to('MJy/sr')
    spec_ical[i][offset:offset+NUM_FREQ] = (bb_ical(nu)*ical_emiss).to('MJy/sr')/otf
    spec_skyh[i][offset:offset+NUM_FREQ] = (bb_skyhorn(nu)*skyhorn_emiss).to('MJy/sr')/otf
    spec_refh[i][offset:offset+NUM_FREQ] = (bb_refhorn(nu)*refhorn_emiss).to('MJy/sr')/otf
    spec_dih[i][offset:offset+NUM_FREQ] = (bb_dihedr(nu)*dihedral_emiss).to('MJy/sr')/otf
    spec_bol[i][offset:offset+NUM_FREQ] = (bb_bolassem(nu)*bolom_emiss).to('MJy/sr')/otf
    spec_struct[i][offset:offset+NUM_FREQ] = (bb_struct(nu)*struct_emiss).to('MJy/sr')/otf


#axes[1].plot(nu, theory*0 + np.nan, color='r', label='Model')
axes[0].legend()
#axes[1].legend()
axes[-1].set_xlabel('Frequency, GHz')
axes[0].set_xlim([0, 700])

plt.savefig('calibrated_versus_model.png', bbox_inches='tight')

spec_th[~np.isfinite(spec_th)] = 0
spec_ical[~np.isfinite(spec_ical)] = 0
spec_ical[~np.isfinite(spec_ical)] = 0

#plt.show()


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

plt.figure()
plt.plot(ifg_xcal[0], label='XCAL IFG')
plt.savefig('xcal_ifg.png')
