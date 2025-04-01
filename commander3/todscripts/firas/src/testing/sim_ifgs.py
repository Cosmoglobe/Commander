import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append('..')
import h5py
import my_utils as mu
from astropy.io import fits
from utils.config import gen_nyquistl


def convert_gain(gain_array):
    conv = {0: 1, 1: 3, 2: 10, 3: 30, 4: 100, 5: 300, 6: 1000, 7: 3000}
    gain_array = gain_array.astype(float)
    for i in range(8):
        gain_array[gain_array.astype(int) == i] = conv[i]
    gain_array[(gain_array < 0)] = np.nan

    return gain_array

sdf = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5')
eng = h5py.File('/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5')


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


pub_model = fits.open(f'/mn/stornext/d16/cmbco/ola/firas/pub_calibration_model/FIRAS_CALIBRATION_MODEL_{ch.upper()}{sm.upper()}.FITS')

apod = pub_model[1].data['APODIZAT'][0]
# apod = np.ones_like(apod)
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

ind = np.arange(100_000,100_050)
ind = np.arange(54_825,54_875)
ind = np.arange(529_025,529_100)
#ind = np.arange(529_000,529_100)

ind = np.arange(54_800, 55_000)


channel = 3
mtm_speed = mtm_speeds[ind]
mtm_length = mtm_lengths[ind]
sweeps = sweepss[ind]
gain = gains[ind]
xcal_pos = xcal_poss[ind]
adds_per_group = adds_per_groups[ind]
eng_time = eng_times[ind]
ifg = ifgs[ind].astype('float32')





#print('xcal position 1 means in sky horn, 2 means out of horn, 3 is moving, 0 is error.')
#print(np.unique(xcal_pos))

bla = (mtm_speed == 0) & (mtm_length == 0)
mtm_speed = mtm_speed[bla]
mtm_length = mtm_length[bla]
sweeps = sweeps[bla]
gain = gain[bla]
xcal_pos = xcal_pos[bla]
adds_per_group = adds_per_group[bla]
eng_time = eng_time[bla]
ifg = ifg[bla]


eng_time_array = eng['ct_head/time'][()]
bol_cmd_biass = eng['en_stat/bol_cmd_bias'][()]
bol_volts = eng['en_analog/group1/bol_volt'][()]


icals = eng['en_analog/grt/a_lo_ical'][()]
xcals = eng['en_analog/grt/a_lo_xcal_cone'][()]
skyhorn = eng['en_analog/grt/a_lo_skyhorn'][()]
refhorn = eng['en_analog/grt/a_lo_refhorn'][()]
dihedral = eng['en_analog/grt/a_lo_dihedral'][()]


en_xcal = eng['en_xcal/pos'][()]


#inds = (en_xcal[:,0] == 1) & (en_xcal[:,1] == 1)
#
fig, axes = plt.subplots(sharex=True, nrows=3)



#secax = axes[0].secondary_xaxis('top', functions=(deg2rad, rad2deg))
axes[0].plot(eng_time_array,icals)
axes[0].plot(eng_time_array,xcals)

axes[1].plot(eng_time_array, en_xcal[:,0])
axes[1].plot(eng_time_array, en_xcal[:,1])

axes[2].plot(sci_times,mtm_speeds)
axes[2].plot(sci_times,mtm_lengths)
plt.show()


eng_inds = []
eng_time_values = np.arange(len(eng_time_array))
for i in range(len(eng_time)):
    eng_ind_matches = eng_time_values[eng_time_array==eng_time[i]]
    if len(eng_ind_matches) == 1:
        eng_inds.append(eng_ind_matches[0])
    else:
        eng_inds.append(eng_inds[-1])
#eng_inds = eng_time_array == eng_time


bol_cmd_bias = bol_cmd_biass[eng_inds,channel]
bol_volt = bol_volts[eng_inds,channel]

icals = icals[eng_inds]
xcals = xcals[eng_inds]
skyhorn = skyhorn[eng_inds]
refhorn = refhorn[eng_inds]


fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

frec = 4 * (channel % 2) + scan_mode

fnyq_icm = fnyq['icm'][frec]
fnyq_hz = fnyq['hz'][frec]



for i in range(len(ifg)):
    ifg[i] -= np.median(ifg[i])


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


from astropy import units as u
from astropy.modeling.models import BlackBody

fig, axes = plt.subplots(sharex=True, nrows=3, sharey=False)
nu = np.arange(len(ical_emiss))*DELTA_NU + NU_ZERO
nu *= u.GHz
spec_th = np.zeros_like(spec_out)
spec_xcal = np.zeros_like(spec_out)
spec_ical = np.zeros_like(spec_out)
spec_skyh = np.zeros_like(spec_out)
spec_refh = np.zeros_like(spec_out)
spec_dih = np.zeros_like(spec_out)
for i in range(len(xcals)):
    bb_xcal = BlackBody(temperature=xcals[i]*u.K)
    bb_ical = BlackBody(temperature=icals[i]*u.K)
    bb_refhorn = BlackBody(temperature=refhorn[i]*u.K)
    bb_skyhorn = BlackBody(temperature=skyhorn[i]*u.K)
    bb_dihedr = BlackBody(temperature=dihedral[i]*u.K)

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
            ).to('MJy/sr')/otf
    R[~np.isfinite(R)] = 0
    # Equation (13) of the Explanatory Supplement, without phase correction.
    D = spec_out[i][5:len(otf)+5]
    S_sky = D - R.value
    theory = R +  bb_xcal(nu).to('MJy/sr')
    axes[1].plot(nu, D, color='k', alpha=0.1)
    #axes[1].plot(nu, theory, alpha=0.1, color='r', zorder=5)
    axes[1].set_ylim([-500, 500])
    axes[1].set_title('Calibrated spectra')

    axes[2].plot(nu, theory.real, color='k', alpha=0.1)
    #axes[2].plot(nu, theory.imag, alpha=0.1, color='r', zorder=5)
    #axes[2].set_ylim([-500, 500])
    axes[2].set_title('Theory spectra')

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



ifg_th = mu.spec_to_ifg(spec_th, mtm_speed, channel, adds_per_group,
                        bol_cmd_bias, bol_volt, fnyq_icm, fnyq_hz, otf, Jo, Jg,
                        Tbol, rho, R0, T0, beta, G1, C3, C1, gain,sweeps, apod)
ifg_th[~np.isfinite(ifg_th)] = 0

ifg_xcal = mu.spec_to_ifg(spec_xcal, mtm_speed, channel, adds_per_group,
                        bol_cmd_bias, bol_volt, fnyq_icm, fnyq_hz, otf, Jo, Jg,
                        Tbol, rho, R0, T0, beta, G1, C3, C1, gain,sweeps, apod)
ifg_xcal[~np.isfinite(ifg_xcal)] = 0

ifg_ical = mu.spec_to_ifg(spec_ical, mtm_speed, channel, adds_per_group,
                        bol_cmd_bias, bol_volt, fnyq_icm, fnyq_hz, otf, Jo, Jg,
                        Tbol, rho, R0, T0, beta, G1, C3, C1, gain,sweeps, apod)

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
plt.show()

plt.figure()
plt.pcolormesh(ifg, vmin=-50, vmax=50, cmap='RdBu_r')
plt.title('Data')

plt.figure()
plt.pcolormesh(ifg_th, vmin=-50, vmax=50, cmap='RdBu_r')
plt.title('Theory')

plt.figure()
plt.title('Model, real')
plt.pcolormesh(spec_th.real, vmin=-10, vmax=10)
plt.figure()
plt.title('Model, imag')
plt.pcolormesh(spec_th.imag, vmin=-10, vmax=10)

plt.figure()
plt.title('Data, real')
plt.pcolormesh(spec_out.real, vmin=-1e2, vmax=1e2)
plt.figure()
plt.pcolormesh(spec_out.imag, vmin=-1e2, vmax=1e2)
plt.title('Data, imag')


plt.figure()
plt.title('Data - Model, real')
plt.pcolormesh(spec_out.real - spec_th.real, vmin=-1e2, vmax=1e2)
plt.figure()
plt.title('Data - Model, imag')
plt.pcolormesh(spec_out.imag - spec_th.imag, vmin=-1e2, vmax=1e2)



#plt.show()




# ifg_out = mu.spec_to_ifg()
