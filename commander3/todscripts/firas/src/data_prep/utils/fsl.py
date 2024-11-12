'''FIRAS Spectra Long

This program drives the FIRAS calibration routines. The program reads coadded IFGs
and write calibrated spectra to disk. The programs reads in the FIRAS calibration
model solution from a reference archive

For a given channel and scan mode, FSL reads in coadded IFGs and transforms them into
voltage spectra. It then applies the calibration model to both the spectra and the
spectrum variances. FSL calibrates both sky and calibration data
'''
import os
import datetime

import numpy as np
from astropy.io import fits
import scipy.interpolate
import h5py

# import firas_pipeline.frd as frd
# import firas_pipeline.fut as fut
# from firas_pipeline.params import *
# import firas_pipeline.misc.util as util
# from firas_pipeline.dataio.config import read_config
# import firas_pipeline.dataio2 as dataio 
# import firas_pipeline.data_types as dt
                
import matplotlib.pyplot as plt

fcc_varflag = -9999.0
trlim = 1.0e-10
sky_temp = 2.726
sky_temp_sig = 2.16659e-5
rscale = 1.0e-7

def produce_spectra_common(coadd_rec, chan, scan_mode, apodl_all, etfl_all):
    '''This function contains the common parts of ffl.produce_spectra and fsl.produce_spectra.
    This differences in the two relate to how the data is stored after producing the spectra.
    FFL has been removed so this could be merged back into produce_spectra.
    '''
   
    #TODO: Why not just read the chan and scan_mode from the coadd record?

    #Initialize the spectrum record
    fftlen = coadd_rec['coad_data']['fft_length']
    speclen = fftlen // 2 + 1
    
    upmode = coadd_rec['coad_spec_data']['sci_mode']
    fakeit = coadd_rec['coad_spec_data']['fakeit']
    speed = coadd_rec['coad_spec_data']['mtm_speed']
    peakpos = coadd_rec['coad_spec_data']['peak_pos']
    
    #Extract the IFG from the coadd record -- for scan_mode = 0-3 extract the full 512 point IFG,
    #for scan_mode = 5 or 6, the IFG has already been decimated or truncated by FIL, so extract
    #the first 128 points of the coadd
    daifg = np.zeros(512)
    if scan_mode <= 3:
        daifg[:] = coadd_rec['coad_data']['ifg']
    else:
        daifg[:128] = coadd_rec['coad_data']['ifg'][:128]
    
    #Read the appropriate electronics transfer function for a new microproc mode
    #Note that for getting the ETF record number for scan modes FS and FL, adds/gp of 2, not 8 must be used
    ngroup = coadd_rec['coad_spec_data']['adds_per_group']
    ngpetf = ngroup
    if scan_mode >= 4: ngpetf = 2
    
    #Generate the voltage spectrum

    #Apodize and rotate the IFG in preparation for the FFT
    aifg = daifg
    arecno = fut.apod_recnuml(chan, scan_mode, fakeit, upmode, ngroup, 0)
    apodl = apodl_all[arecno, :]

    rifg = fut.apod_rotl(aifg, apodl, fftlen, peakpos)

    #Do the FFT
    vspec = np.fft.rfft(rifg)

    #Normalize the spectrum and covert it from counts to volts. Fix the short FFT for the low
    #frequency fast spectra. Shift the phase of the high frequency short fast spectra.
    if (chan == 0 or chan == 2) and scan_mode == 1:
        phase_shift = np.pi / fftlen * 1j
    else:
        phase_shift = 0.0

    erecno = fut.get_recnum(fakeit, speed, chan, upmode, ngpetf)

    etfl = etfl_all[erecno, :]

    #Calculate the normalization for the FFTed spectra.
    #The normalization factor is the product of the Nyquist frequency, the instrument
    #throughput, and the ADC voltage scale factor. The signs of the spectra for the right
    #side channels is flipped so the right side transfer functions will be positive
    fnyq_icm = coadd_rec['coad_spec_data']['nyquist_icm']
    if chan < 2:
        spec_norm = - fnyq_icm * fac_etendu * fac_adc_scale
    else:
        spec_norm = fnyq_icm * fac_etendu * fac_adc_scale
    
    jvals = np.arange(1, speclen)
    
    vspec[1:speclen] *= np.exp(jvals*phase_shift) / (spec_norm * etfl[1:speclen])

    return vspec

def produce_spectra(coadd_rec, chan, scan_mode, apodl_all, etfl_all):
    '''This function apodizes and rotates the coadded IFG, FFTs the IFG into
    a spectrum, and converts the spectrum from counts to volts/cm**2/sr/icm
    '''

    #Generate the voltage spectrum. Using the common function
    vspec = produce_spectra_common(coadd_rec, chan, scan_mode, apodl_all, etfl_all)
    
    #Fill in the voltage spectrum record
    vspec_rec = np.array(0, dtype=dt.dt_fsl_sky)

    #Copying over the records that are common between the input and output data
    vspec_rec['ct_head'] = coadd_rec['ct_head']
    vspec_rec['coad_spec_head'] = coadd_rec['coad_spec_head']
    vspec_rec['coad_spec_data'] = coadd_rec['coad_spec_data']
    vspec_rec['en_stat'] = coadd_rec['en_stat']
    vspec_rec['en_analog'] = coadd_rec['en_analog']
    vspec_rec['en_sigma'] = coadd_rec['en_sigma']
    vspec_rec['en_tempdiff'] = coadd_rec['en_tempdiff']
    vspec_rec['attitude'] = coadd_rec['attitude']

    npts = len(vspec)
    vspec_rec['spec_data']['spec'][:npts] = vspec
    #vspec_rec['spec_data']['real_var'] = vvar[0, :]
    #vspec_rec['spec_data']['imag_var'] = vvar[1, :]
    #vspec_rec['spec_data']['real_imag_var'] = vvar[2, :]
    vspec_rec['spec_data']['real_var'] = coadd_rec['coad_data']['real_var']
    vspec_rec['spec_data']['imag_var'] = coadd_rec['coad_data']['imag_var']
    vspec_rec['spec_data']['real_imag_var'] = coadd_rec['coad_data']['real_imag_var']

    vspec_rec['spec_data']['fft_length'] = coadd_rec['coad_data']['fft_length']

    return vspec_rec

def calibrate_spectra2(vspec_recs, model, config, chan: int, scan_mode: int, single: bool = False, input_type: str = 'sky',
                       tsig0=False, fcc_flv=False, diff: bool = False, dvector=None):
    '''This routine produces calibrated spectra and calibrated spectrum variances
    with the units of ergs/sec/cm**2/sr/icm. The spectra and variances are then
    converted to units of MJr/sr for FSL output records. The sky spectra are in
    the barycentric frame of reference. This function runs on all spectra instead
    of single records
    '''

    nspec = len(vspec_recs)

    vvars = np.zeros([nspec, 3, 361], dtype=np.double)
    vvars[:, 0, :] = vspec_recs['spec_data']['real_var']
    vvars[:, 1, :] = vspec_recs['spec_data']['imag_var']
    vvars[:, 2, :] = vspec_recs['spec_data']['real_imag_var']

    #grtcoawt = config['grtcoawt']
    #grttrans = config['grttrans']

    par = model['bolparm']
    emiss = model['emissivity']

    cspec_recs = np.copy(vspec_recs)
    cspec_recs['spec_data'] = 0

    idxs = np.where(np.abs(emiss[0, :]) > trlim)[0]
    
    nifgs = cspec_recs['coad_spec_head']['num_ifgs']
    temps = cspec_recs['coad_spec_data']['temp']
    tsigs = cspec_recs['coad_spec_data']['temp_sigma']
    
    #Calibrate the variances and the spectrum

    #Calculate the temporal drift corrections
    t0 = np.array(cspec_recs['ct_head']['time'], dtype=np.int64)
    primary_vibs, ticalds = temporal_drift2(t0, par)
    
    #Adjust the Xcal temperature if we are applying the calibration model
    #to the calibration data
    #if input_type == 'cal':
    #    temps[:, 0] += par[20]

    temps[:, 1] -= ticalds

    #Calculate the actual detector temperature and the detector responsitivty and
    #time constants
    bol_volts = cspec_recs['coad_spec_data']['bol_volt']
    bol_cmd_bias = cspec_recs['coad_spec_data']['bol_cmd_bias']
    bol_cmd_bias[bol_cmd_bias < 0] += 256
    cmd_bias = np.double(bol_cmd_bias) / 25.5
    tdets = temps[:, chan+6]
    tbols, qrad, S0, tau = calc_responsivity(bol_volts, cmd_bias, tdets, model['bolparm'])
    temps[:, 6:10] = tbols[:, None]

    fftlen = cspec_recs[0]['spec_data']['fft_length']
    consts = compute_constants(chan, scan_mode, config, fft_len=fftlen)
    #Compute the bolometer delay
    #B = 1.0 + np.outer(tau, consts['afreq'])*1j
    B = 1.0 + consts['afreq'][None, :]*tau[:, None]*1j
    B[:, 0] = 0.0

    if single:
        #Set the calibrated variances to flag values for single IFG spectra
        cvars = fcc_varflag * np.ones([nspec, 3, 361])
        cvars_out = fcc_varflag * np.ones_like(cvars)
    else:
        #If FIL averaged variances are used, weight the FIL variance vectors by the glitch rate
        #adjusted number of IFGs for the coadd and replace the individual coadd voltage variance
        #vectors with the weighted FIL variances
        cvars_out = np.zeros([nspec, 3, 361])
        i0 = nifgs >= 3
        i1 = np.logical_and(nifgs == 2, ~cspec_recs['coad_spec_data']['sec_template']['subtracted'])
        i0 = np.logical_or(i0, i1)
        cvars_out[i0, :, :] = pcalib_variances2(vvars[i0, :, :], emiss, S0[i0], B[i0, :], idxs)
        cvars_out[~i0, :, :] = fcc_varflag

        adjnifgs = vspec_recs['coad_spec_head']['adj_num_ifgs']
        if fcc_flv:
            vvars = fil_var / adjnifgs

        #Apply the calibration model to the variances for coadded IFG spectra
        cvars = pcalib_variances2(vvars, emiss, S0, B, idxs)

    #if input_type == 'cal':
    #    xcalin = True
    #else:
    xcalin = False

    #Apply the calibration model to the spectra. Add adjnifgs to call
    if diff:
        cspecs, pspecs, power, phase_corr = apply_model2(model, idxs, vspec_recs, primary_vibs, temps, tsigs, cvars, config,
                                                         consts, S0, tau, diff=diff, input_type=input_type, dvector=dvector,
                                                         xcalin=xcalin, fcc_flv=fcc_flv)
    else:
        cspecs, power, phase_corr = apply_model2(model, idxs, vspec_recs, primary_vibs, temps, tsigs, cvars, config,
                                                 consts, S0, tau, diff=diff, input_type=input_type, dvector=dvector,
                                                 xcalin=xcalin, fcc_flv=fcc_flv)
   
    if input_type == 'sky':
        #Correct the sky spectra for the projected spacecraft barycentric velocity along the FIRAS
        #line of sight, checking for valid projected barycentric velocities
        vbarys = vspec_recs['attitude']['projected_barycentric_velocity'] / 100.0

        #We do a loop here because theoretically must do a separate interpolation
        #for each record
        for i0, vbary in enumerate(vbarys):
            vblim = 40.0
            if np.abs(vbary) > vblim:
                raise ValueError("vbary too high")
            elif vbary != 0:
                doppler_shift(vbary, cspecs[i0, :], idxs, consts['freq'])

    #Update the spectrum record
    #Convert the spectrum and variances from ergs/sec/cm**2/sr/icm to MJy/sr
    cspecs *= fac_erg_to_mjy
    erg_to_mjy_sq = fac_erg_to_mjy**2
    cvars *= erg_to_mjy_sq
    cvars_out *= erg_to_mjy_sq

    #Set the calibrated flag and the flag for one of the three types of variances used for
    #the autophase corrector computation
    if fcc_flv or (dvector is not None):
        cspec_recs['spec_data']['calibrated'] = 1
    else:
        i0 = nifgs >= 3
        cspec_recs[i0]['spec_data']['calibrated'] = 1
        cspec_recs[~i0]['spec_data']['calibrated'] = 0

    if fcc_flv:
        cspec_recs['spec_data']['fil_vars'] = 1
    elif (dvector is not None):
        cspec_recs['spec_data']['fsl_vars'] = 1
    else:
        cspec_recs['spec_data']['coadd_vars'] = 1

    #Write the calibrated spectrum and variances into the spectrum record. Use the
    #start and stop frequency indices from the FSL_model include file to store only
    #valid spectral values over the range where the model solution exists
    cspec_recs['spec_data']['spec'][:, idxs] = cspecs[:, idxs]

    cvars2 = np.zeros_like(cvars)
    if dvector is not None:
        #TODO: See what I would need to do if I read in dvector
        cvars2[:, 0, idxs] = dvector[idxs]
    else:
        cvars2[:, :, :] = cvars_out[:, :]

    i0 = nifgs == 1
    i1 = np.logical_and(nifgs == 2, cspec_recs['coad_spec_data']['sec_template']['subtracted'])
    i0 = np.logical_or(i0, i1)
    cspec_recs['spec_data']['real_var'][i0, :] = fcc_varflag
    cspec_recs['spec_data']['imag_var'][i0, :] = fcc_varflag
    cspec_recs['spec_data']['real_imag_var'][i0, :] = fcc_varflag
    cspec_recs['spec_data']['real_var'][~i0, :] = cvars2[~i0, 0, :]
    cspec_recs['spec_data']['imag_var'][~i0, :] = cvars2[~i0, 1, :]
    cspec_recs['spec_data']['real_imag_var'][~i0, :] = cvars2[~i0, 2, :]

    #Store field values in output spectrum record
    cspec_recs['spec_data']['fft_length'] = vspec_recs['spec_data']['fft_length']
    cspec_recs['spec_data']['lofreq_bin'] = np.min(idxs)
    cspec_recs['spec_data']['hifreq_bin'] = np.max(idxs)

    #Fill in the model solution fields of the spectrum record. Update the ical
    #and detector temperature fields in the spectrum record
    cspec_recs['spec_data']['responsivity'] = S0
    cspec_recs['spec_data']['time_constant'] = tau
    cspec_recs['spec_data']['phase_corr'] = phase_corr
    cspec_recs['spec_data']['qrad'] = qrad
    cspec_recs['spec_data']['ir_power'] = power
    cspec_recs['coad_spec_data']['temp'][:, 1] = temps[:, 1]
    cspec_recs['coad_spec_data']['temp'][:, 6:10] = temps[:, 6:10]

    return cspec_recs, cvars

def calibrate_spectra(vspec_rec, model, config, chan, scan_mode, single=False, input_type='sky', 
                      tsig0=False, fcc_flv=False, diff=False, dvector=None):
    '''These routines produce calibrated spectra and calibrated spectrum variances
    with the units of ergs/sec/cm**2/sr/icm. The spectra and variances are then
    converted to units of MJr/sr for FSL output records. The sky spectra are in
    the barycentric frame of reference.
    '''

    vvar = np.zeros([3, 361], dtype=np.double)
    vvar[0, :] = vspec_rec['spec_data']['real_var']
    vvar[1, :] = vspec_rec['spec_data']['imag_var']
    vvar[2, :] = vspec_rec['spec_data']['real_imag_var']
    
    grtcoawt = config['grtcoawt']
    grttrans = config['grttrans']

    par = model['bolparm']
    emiss = model['emissivity']

    cspec_rec = np.copy(vspec_rec)
    cspec_rec['spec_data'] = 0

    idxs = np.where(np.abs(emiss[0, :]) > trlim)[0]

    nifgs = cspec_rec['coad_spec_head']['num_ifgs']
    temp = cspec_rec['coad_spec_data']['temp']
    tsig = cspec_rec['coad_spec_data']['temp_sigma']

    #Calibrate the variances and the spectrum

    #Calculate the temporal drift corrections
    primary_vib, ticald = temporal_drift(int(cspec_rec['ct_head']['time']), par)
    #temp[0] += par[20]
    temp[1] -= ticald
    
    #Calculate the actual detector temperature and the detector responsitivty and
    #time constance
    bol_volt = vspec_rec['coad_spec_data']['bol_volt']
    bol_cmd_bias = cspec_rec['coad_spec_data']['bol_cmd_bias']
    if bol_cmd_bias < 0: bol_cmd_bias += 256
    cmd_bias = np.double(bol_cmd_bias) / 25.5
    tdet = temp[chan+6]
    tbol, qrad, S0, tau = calc_responsivity(bol_volt, cmd_bias, tdet, model['bolparm'])

    temp[6:10] = tbol

    fftlen = vspec_rec['spec_data']['fft_length']
    consts = compute_constants(chan, scan_mode, config, fft_len=fftlen)
    #Compute the bolometer delay
    #B = np.zeros_like(afreq)
    #idx == np.abs(emiss[0, :] > trlim)
    #B[idx] = 1.0 + afreq[idx]*tau*1j
    B = 1.0 + consts['afreq']*tau*1j
    B[0] = 0.0

    if single:
        #Set the calibrated variances to flag values for single IFG spectra
        cvar = fcc_varflag * np.ones([3, 361])
        cvar_out = fcc_varflag * np.ones_like(cvar)
    else:
        #If FIL averaged variances are used, weight the FIL variance vectors by the glitch rate
        #adjusted number of IFGs for the coadd and replace the individual coadd voltage variance
        #vectors with the weighted FIL variances
        if nifgs >= 3 or (nifgs == 2 and not cspec_rec['coad_spec_data']['sec_template']['subtracted']):
            cvar_out = pcalib_variances(vvar, emiss, S0, B, idxs)
        else:
            cvar_out = fcc_varflag * np.ones([3, 361])

        adjnifgs = vspec_rec['coad_spec_head']['adj_num_ifgs']
        if fcc_flv:
            vvar = fil_var / adjnifgs

        #Apply the calibration model to the variances for coadded IFG spectra
        cvar = pcalib_variances(vvar, emiss, S0, B, idxs)

    #if input_type == 'cal':
    #    xcalin = True
    #else:
    xcalin = False

    #Apply the calibration model to the spectra. Add adjnifgs to call
    cspec, power, phase_corr = apply_model(model, idxs, vspec_rec, primary_vib, temp, tsig, cvar, config, consts, S0,
                                           tau, diff=diff, input_type=input_type, dvector=dvector, xcalin=xcalin,
                                           fcc_flv=fcc_flv)
    
    if input_type == 'sky':
        #Correct the sky spectra for the projected spacecraft baryventric velocity along the FIRAS
        #line of sight, checking for valid projected barycentric velocities
        vbary = vspec_rec['attitude']['projected_barycentric_velocity'] / 100.0
        
        vblim = 40.0
        if np.abs(vbary) > vblim:
            raise ValueError("vbary too high")
        elif vbary != 0:
            doppler_shift(vbary, cspec, idxs, consts['freq'])

    #Update the spectrum record

    #Convert the spectrum and variances from ergs/sec/cm**2/sr/icm to MJy/sr
    cspec *= fac_erg_to_mjy
    erg_to_mjy_sq = fac_erg_to_mjy**2
    cvar *= erg_to_mjy_sq
    cvar_out *= erg_to_mjy_sq

    #Set the calibrated flag and the flag for one of the three types of variances used for
    #the autophase corrector computation
    if fcc_flv or (dvector is not None) or nifgs >= 3:
        cspec_rec['spec_data']['calibrated'] = 1
    else:
        cspec_rec['spec_data']['calibrated'] = 0

    if fcc_flv:
        cspec_rec['spec_data']['fil_vars'] = 1
    elif (dvector is not None):
        cspec_rec['spec_data']['fsl_vars'] = 1
    else:
        cspec_rec['spec_data']['coadd_vars'] = 1
    
    #Write the calibrated spectrum and variances into the spectrum record. Use the
    #start and stop frequency indices from the FSL_model include file to store only
    #valid spectral values over the range where the model solution exists

    #spec_rec['spec_data']['spec'][start_idx:stop_idx+1] = display.cspec[start_idx:stop_idx+1]
    #cspec_rec['spec_data']['spec'][start_idx:stop_idx+1] = cspec[start_idx:stop_idx+1]
    cspec_rec['spec_data']['spec'][idxs] = cspec[idxs]
    
    cvar2 = np.zeros_like(cvar)
    if dvector is not None:
        cvar2[0, idxs] = dvector[idxs]
    else:
        cvar2[:, :] = cvar_out[:, :]

    if nifgs == 1 or (nifgs == 2 and cspec_rec['coad_spec_data']['sec_template']['subtracted']):
        cspec_rec['spec_data']['real_var'][idxs] = fcc_varflag
        cspec_rec['spec_data']['imag_var'][idxs] = fcc_varflag
        cspec_rec['spec_data']['real_imag_var'][idxs] = fcc_varflag
    else:
        cspec_rec['spec_data']['real_var'] = cvar2[0, :]
        cspec_rec['spec_data']['imag_var'] = cvar2[1, :]
        cspec_rec['spec_data']['real_imag_var'] = cvar2[2, :]

    #Store field values in output spectrum record
    cspec_rec['spec_data']['fft_length'] = vspec_rec['spec_data']['fft_length']
    cspec_rec['spec_data']['lofreq_bin'] = np.min(idxs)
    cspec_rec['spec_data']['hifreq_bin'] = np.max(idxs)

    #Fill in the model solution fields of the spectrum record. Update the ical
    #and detector temperature fields in the spectrum record
    cspec_rec['spec_data']['responsivity'] = S0
    cspec_rec['spec_data']['time_constant'] = tau
    cspec_rec['spec_data']['phase_corr'] = phase_corr
    cspec_rec['spec_data']['qrad'] = qrad
    cspec_rec['spec_data']['ir_power'] = power
    cspec_rec['coad_spec_data']['temp'][1] = temp[1]
    cspec_rec['coad_spec_data']['temp'][6:10] = temp[6:10]
    
    #TODO: Generate model tag and label
    #cspec_rec['spec_data']['model_ttag'] = model_ttag
    #cspec_rec['spec_data']['model_label'] = model_label
 
    return cspec_rec, cvar

def doppler_shift(vbary, cspec, idxs, freq):
    '''This function Doppler shifts a calibrated spectrum from the spacecraft frame of reference
    into the barycentric frame. The frequency shift is performed by a cubic spline interpolation,
    applied to the real and imaginary parts of the calibrated spectrum separately. The frequency
    scaling factor is (1+v/c), while the photometric scaling factor is (1-v/c)**3. Only first-order
    effects are considered since |v/c| << 0.001.
    '''

    vlight = 2.99792458e5

    #Extract the real and imaginary parts of the calibrated spectrum. The emissivity arrays in
    #the FSL_model include file only have 257 frequency points. The non-zero values start and
    #stop indices have already been determined. The spectra arrays are dimensioned for 361 points
    npts = len(idxs)
    offset = np.min(idxs)

    freq0 = freq[idxs]
    spec = cspec[idxs]

    #Calculate the spline interpolation coefficients
    f = scipy.interpolate.interp1d(freq0, spec, kind='cubic')

    #Calculate the frequency points in the barycentric frame. The Doppler-shifted frequency
    #array consists of the original frequencies scaled by (1+v/c)
    beta_plus = 1.0 + vbary/vlight
    freq1 = freq0 * beta_plus
    #freqsp[freqsp >= freq[fcc_spec_length]] = freq[fcc_spec_length]
    
    #NOTE: Can't do spline outside initial region, but Fortran code doesn't seem to have that problem
    #This code makes sure all frequencies are within the range, but should I add some extrapolation?
    freq1[freq1 >= freq0[-1]] = freq0[-1]
    freq1[freq1 <= freq0[0]] = freq0[0]
    
    #Interpolate the calibrated spectrum onto the barycentric frequency array
    dopspec = f(freq1)

    #Rescale the Doppler-shifted spectrum by the photometric scaling factor of (1-v/c)**3.
    #Write the scaled, shifted spectrum back into the calibrated spectrum array.
    beta_minus = (1.0 - vbary/vlight)**3
    cspec[idxs] = dopspec*beta_minus

def apply_model2(model, idxs, vspec_recs, primary_vib, temps, tsigs, cvars, config, consts, S0,
                 tau, diff=False, input_type='sky', dvector=None, xcalin=False, phase_corr=None,
                 fcc_flv=False, derivs=None):

    emiss = model['emissivity']
    param = model['bolparm']

    B = 1.0 + 1j*np.outer(tau, consts['afreq'])

    nifgs = vspec_recs['coad_spec_head']['num_ifgs']
    adjnumifgs = vspec_recs['coad_spec_head']['adj_num_ifgs']
    chan = vspec_recs['coad_spec_data']['chan_id']
    mtm_speed = vspec_recs['coad_spec_data']['mtm_speed']
    mtm_length = vspec_recs['coad_spec_data']['mtm_length']

    vspecs = vspec_recs['spec_data']['spec']
    
    ##iv should go (RHSS, RHSF, RHLF, RLSS, RLFS/FL, RLLF, LHSS, LHSF, LHLF, LLSS, LLFS/FL, LLLF)
    #iv = 3*chan + mtm_speed + mtm_length #only 3 input, speed/length=0, speed=0/length=1, speed=1/length=1

    if input_type == 'sky':
        temps[:, 0] = sky_temp
        tsigs[:, 0] = sky_temp_sig

    #Calibrate the differental spectrum
    specs = np.zeros_like(vspecs)
    cspecs = np.zeros_like(vspecs)
    variances = np.zeros_like(vspecs)
    
    if derivs is not None:
        nrecs, nfreq = vspecs.shape
        npar = len(param)
        ddetparam = np.zeros([nrecs, nfreq, npar], dtype=complex)
    
    #Convert the spectrum from volts/cm**2/sr/icm to ergs/sec/cm**2/sr/icm.
    #Correct the spectrum for the 3rd and 2nd harmonic terms and for the secondary
    #and primary (time-dependent) vibration terms.
    sec_vib = param[14]

    corrections = harm_vib_correction2(vspec_recs, B, config, consts, tau, idxs, derivs=derivs)

    harm_2 = corrections[0]
    harm_3 = corrections[1]
    prim_vib_corr = corrections[2]
    sec_vib_corr = corrections[3]

    denom = np.outer(S0, emiss[0, idxs])

    specs[:, idxs] = (B[:, idxs] * vspecs[:, idxs] - param[9] * harm_3 - param[10] * harm_2 - 
                      sec_vib * sec_vib_corr - primary_vib[:, None] * prim_vib_corr) / denom

    if derivs is not None:
        
        #Harmonic derivative stuff
        Zb = B[:, idxs]
        dA = 1.0/(2.0*Zb) + 1.0/(2.0*(2.0-Zb))-1.0 #dA/A
        dZb = 1.0 - 1.0/Zb #dZb/Zb
        dZb_A = dA + dZb #d(A*Zb)/(A*Zb)

        #Harmonic corrections
        #ddetparam[:, idxs, 10] = -harm_2 / denom
        #ddetparam[:, idxs, 9] = -harm_3 / denom
        ddetparam[:, idxs, 10] = harm_2 / denom
        ddetparam[:, idxs, 9] = harm_3 / denom

        #model = (param[9] * harm_3 + param[10] * harm_2)/denom
        model = param[9] * ddetparam[:, idxs, 9] + param[10] * ddetparam[:, idxs, 10]

        #Keep is separate from XT right now because it needs to be multiplied by
        #the phase correction
        XT1 = param[9] * ddetparam[:, idxs, 9]*2*(dA + dZb*Zb/(2+Zb)) + \
              param[10] * ddetparam[:, idxs, 10]*(dA + dZb*Zb/(1.0+Zb))

        #Vibration corrections
        #dvib is d(primary_vib)/d(vib_param)
        #In FISH spectra is modified as I = I - par(m)*X(m), where X is
        #the derivative stuff
        dvib = derivs['dvib']
        #ddetparam[:, idxs, 14] = -sec_vib_corr / denom
        #ddetparam[:, idxs, 15] = -dvib[:, None, 0] * prim_vib_corr / denom
        #ddetparam[:, idxs, 16] = -dvib[:, None, 1] * prim_vib_corr / denom
        #ddetparam[:, idxs, 17] = -dvib[:, None, 2] * prim_vib_corr / denom
        #ddetparam[:, idxs, 18] = -dvib[:, None, 3] * prim_vib_corr / denom
        #ddetparam[:, idxs, 19] = -dvib[:, None, 4] * prim_vib_corr / denom 
        ddetparam[:, idxs, 14] = sec_vib_corr / denom
        ddetparam[:, idxs, 15] = dvib[:, None, 0] * prim_vib_corr / denom
        ddetparam[:, idxs, 16] = dvib[:, None, 1] * prim_vib_corr / denom
        ddetparam[:, idxs, 17] = dvib[:, None, 2] * prim_vib_corr / denom
        ddetparam[:, idxs, 18] = dvib[:, None, 3] * prim_vib_corr / denom
        ddetparam[:, idxs, 19] = dvib[:, None, 4] * prim_vib_corr / denom

        #model += (sec_vib * sec_vib_corr + primary_vib[:, None] * prim_vib_corr) / denom
        for i in range(14, 20):
            model += param[i] * ddetparam[:, idxs, i]

        #Vibration derivative stuff
        Bvp = derivs['Bvp']
        Bvs = derivs['Bvs']
        XT1 += param[14]*ddetparam[:, idxs, 14]*(1+dA - 1/Bvs)
        for k in range(15, 20):
            XT1 += param[k]*ddetparam[:, idxs, k]*(1+dA - 1/Bvp)

        #Ical correction (Not added yet because we are not varying these parameters)
        ddetparam[:, idxs, 11] = 0
        ddetparam[:, idxs, 12] = 0
        ddetparam[:, idxs, 13] = 0

    pspecs, power = reference_spectra2(emiss, temps, tsigs, idxs, chan, consts, xcalin=xcalin,
                                       derivs=derivs)

    #Determine the variance estimate. Weight Dvector by glitch rate adjusted number
    #of IFGs
    if dvector is not None:
        variances = dvector / adjnifgs
    else:
        variances = cvars[:, 0, :]
    
    cspecs[:, idxs] = specs[:, idxs]
  
    #Only calculate a phase correction if it is not input
    if phase_corr is None:
        #Autophase correct the spectrum
        if dvector is not None or fcc_flv:
            phase_corr = autophase_correct2(specs, pspecs, variances, idxs)
        else:
            phase_corr = np.zeros_like(specs[:, 0])
            i0 = nifgs >= 3
            phase_corr[i0] = autophase_correct2(specs[i0, :], pspecs[i0, :], variances[i0, :], idxs)

    pc = np.exp(np.outer(phase_corr, idxs)*1j)
    cspecs[:, idxs] *= pc
    
    if derivs is not None:
        #This is because the bolometer temperature depends on the detector parameters
        Pkd = np.zeros_like(vspecs, dtype=float)
        it = chan + 6
        for i in range(nrecs):
            Pkd[i, idxs] = util.planck_deriv(consts['freq'][idxs], temps[i, it[i]])
       
        #TODO: Figure out where this comes from, multiplied tau in derivative
        element = derivs['element']
        XT = np.sum(element*dA[:, None, :], axis=1)
        XS = pspecs[:, idxs]

        XT += pc * XT1

        ddetparam[:, idxs, 14:20] *= pc[:, :, None]

        model += XS/(Zb*pc)

        #X[:, :, 0] to X[:, :, 8]
        gS = derivs['gS']
        gT = derivs['gT']
        gtau = derivs['gtau']
        fact = emiss[6, idxs] / emiss[0, idxs]
        for j in range(9):
            ddetparam[:, idxs, j] = gS[:, None, j]*XS
            ddetparam[:, idxs, j] += gtau[:, None, j]*XT
            ddetparam[:, idxs, j] -= gT[:, None, j]*fact[None, :]*Pkd[:, idxs]
            ddetparam[:, idxs, j] -= phase_corr[:, None]*vspecs[:, idxs]*Zb*dZb_A*gtau[:, None, j]

        dpercoadd = derivs['percoadd']
        #TODO: Check minus sign
        cspec2 = (cspecs[:, idxs] - pspecs[:, idxs])
        dpercoadd[:, 4, :] = idxs[None, :]*1j*cspecs[:, idxs]

        derivs['detparam'] = ddetparam[:, idxs, :]
        derivs['dpercoadd'] = dpercoadd
    
    #Calibrate the spectrum
    if not diff:
        cspecs[:, idxs] -= pspecs[:, idxs]
        return cspecs, power, phase_corr
    else:
        return cspecs, pspecs, power, phase_corr

#TODO: Combine some of these inputs into one dictionary(?)
def apply_model(model, idxs, vspec_rec, primary_vib, temp, tsig, cvar, config, consts, S0, tau,
                diff=False, input_type='sky', dvector=None, xcalin=False, phase_corr=None,
                fcc_flv=False):
    '''This function applies the FIRAS calibration model to a voltage spectrum with the units of
    volts/cm**2/sr/icm, producting a calibrated spectrum with the units of ergs/sec/cm**2/sr/icm.
    The routine calibrates the differential spectrum and produces a reference spectrum. It then
    computes a linear autophase corrector for the spectrum. Finally, it calibrates the spectrum.
    '''
   
    emiss = model['emissivity']
    param = model['bolparm']
    
    B = 1.0 + consts['afreq']*tau*1j

    nifgs = vspec_rec['coad_spec_head']['num_ifgs']
    adjnumifgs = vspec_rec['coad_spec_head']['adj_num_ifgs']
    chan = vspec_rec['coad_spec_data']['chan_id']
    mtm_speed = vspec_rec['coad_spec_data']['mtm_speed']
    mtm_length = vspec_rec['coad_spec_data']['mtm_length']

    vspec = vspec_rec['spec_data']['spec']
    
    ##iv should go (RHSS, RHSF, RHLF, RLSS, RLFS/FL, RLLF, LHSS, LHSF, LHLF, LLSS, LLFS/FL, LLLF)
    #iv = 3*chan + mtm_speed + mtm_length #only 3 input, speed/length=0, speed=0/length=1, speed=1/length=1
    
    if input_type == 'sky':
        #External calibrator is not it so location where it is normally is at the sky temperature
        temp[0] = sky_temp
        tsig[0] = sky_temp_sig

    #Calibrate the differential spectrum
    spec = np.zeros_like(vspec)
    cspec = np.zeros_like(vspec)
    variances = np.zeros_like(vspec)

    #Convert the spectrum from volts/cm**2/sr/icm to ergs/sec/cm**2/sr/icm.
    #Correct the spectrum for the 3rd and 2nd harmonic terms and for the secondary
    #and primary (time-dependent) vibration terms.
    sec_vib = param[14]
            
    #spec[idxs] = (B[idxs] * vspec[idxs] - param[9] * (2.0 + B[idxs])/9.0 * (vspec[if3a] + vspec[if3b] + vspec[if3c])
    #              - param[10] * (1.0 + B[idxs])/4.0 * (vspec[if2a] + vspec[if2b])
    #              - sec_vib * fs * Bvs * vspec[ifs]
    #              - primary_vib * fp * Bvp * vspec[ifp]) / (S0 * emiss[0, idxs])

    #This is done in a function so it can be called from the calibration fitting program
    corrections = harm_vib_correction(vspec_rec, B, config, consts, tau, idxs)

    harm_2 = corrections[0]
    harm_3 = corrections[1]
    prim_vib_corr = corrections[2]
    sec_vib_corr = corrections[3]

    spec[idxs] = (B[idxs] * vspec[idxs] - param[9] * harm_3 - param[10] * harm_2 - 
                  sec_vib * sec_vib_corr - primary_vib * prim_vib_corr) / (S0 * emiss[0, idxs])

    #Compute the reference spectrum (sum of emissivities * Planck function) and the IR power incident
    #on the bolometer. Determine the variance estimate for the autophase correction.
    pspec, power = reference_spectra(emiss, temp, tsig, idxs, chan, consts, xcalin=xcalin)

    #Determine the variance estimate. Weight Dvector by glitch rate adjusted number
    #of IFGs
    if dvector is not None:
        variances = dvector / adjnifgs
    else:
        variances = cvar[0, :]
    
    cspec[idxs] = spec[idxs]
  
    #Only calculate a phase correction if it is not input
    if ((dvector is not None) or nifgs >= 3 or fcc_flv) and phase_corr is None:
        #Autophase correct the spectrum
        phase_corr = autophase_correct(spec, pspec, variances, idxs)
    elif phase_corr is None:
        phase_corr = 0

    cspec[idxs] *= np.exp(idxs*phase_corr*1j)

    #Calibrate the spectrum
    if not diff:
        cspec[idxs] -= pspec[idxs]
        return cspec, power, phase_corr
    else:
        #return cspec, pspec, power, phase_corr
        return cspec, power, phase_corr

def reference_spectra(emiss, temp, tsig, idxs, chan, consts, xcalin=False):
    '''Calculate the reference spectra that is subtracted from the measured sky spectra during calibration
    '''
     
    nfreq = np.shape(emiss)[1]

    pspec = np.zeros(nfreq, dtype=complex)
    
    element = np.zeros([7, nfreq], dtype=complex)
    it = chan + 6
    for k in range(6):
        element[k, idxs] = (emiss[k, idxs] / emiss[0, idxs]) * util.planck_dist(consts['freq'][idxs], temp[k], tsig[k])
    element[6, idxs] = (emiss[6, idxs] / emiss[0, idxs]) * util.planck_dist(consts['freq'][idxs], temp[it], tsig[it])
            
    #Compute the reference spectrum and the total IR power
    power = np.sum(np.abs(element[:, idxs]))
    power *= rscale
  
    #TODO: Is this correct? If xcal is in we want to add its contributions to the reference spectra
    if xcalin:
        pspec[idxs] += np.sum(element[:, idxs], axis=0)
    else:
        pspec[idxs] += np.sum(element[1:, idxs], axis=0)

    return pspec, power

def reference_spectra2(emiss, temp, tsig, idxs, chan, consts, xcalin=False,
                       derivs=None):
    '''Calculate the reference spectra that is subtracted from the measured sky spectra during calibration
    '''

    nfreq = np.shape(emiss)[1]
    nrecs = len(chan)

    pspecs = np.zeros([nrecs, nfreq], dtype=complex)

    #element = np.zeros([nrecs, 7, nfreq], dtype=complex)
    planck_vals = np.zeros([nrecs, 7, nfreq])
    it = chan + 6

    bolo_temp = np.zeros(nrecs)
        
    tmp = emiss[0, idxs]
    fact = emiss[:, idxs] / tmp[None, :]
     
    for k in range(6):
        planck_vals[:, k, idxs] = util.planck_dist2(consts['freq'][idxs], temp[:, k], tsig[:, k])
    idx = (np.arange(nrecs), it)
    planck_vals[:, 6, idxs] = util.planck_dist2(consts['freq'][idxs], temp[idx], tsig[idx])

    if derivs is not None:
        planck_derivs = np.zeros([nrecs, 7, nfreq], dtype=complex)
        #Planck function derivatives are only evaluated at input temperature since they
        #will only be used for calibration where tsig=0
        for m in range(nrecs):
            for k in range(6):
                planck_derivs[m, k, idxs] = util.planck_deriv(consts['freq'][idxs], temp[m, k])
            planck_derivs[m, 6, idxs] = util.planck_deriv(consts['freq'][idxs], temp[m, it[m]])

    #Compute the reference spectrum and the total IR power
    #power = np.sum(np.abs(element[:, :, idxs]), axis=(1, 2))
    power = np.sum(np.abs(fact[None, :, :]*planck_vals[:, :, idxs]), axis=(1, 2))
    power *= rscale

    #TODO: Is this correct? If xcal is in we want to add its contributions to the reference spectra
    #if xcalin:
    #    pspecs[:, idxs] += np.sum(element[:, :, idxs], axis=1)
    #else:
    #    pspecs[:, idxs] += np.sum(element[:, 1:, idxs], axis=1)

    if xcalin:
        pspecs[:, idxs] += np.sum(fact[None, :, :] * planck_vals[:, :, idxs], axis=1)
    else:
        pspecs[:, idxs] += np.sum(fact[None, 1:, :] * planck_vals[:, 1:, idxs], axis=1)

    if derivs is not None:
        #NOTE: I am multiplying by emiss[0, idxs], since this is a gain factor that
        #multiplies both spectra and error. The derivative is trivial to calculate
        #after I multiply by the Xcal emissivity
        planck_bol = planck_vals[:, 6, :]

        #This can be used for both the derivative with respect to the real and imaginary
        #parts of the emissivity
        derivs['emiss'] = planck_vals[:, :6, idxs] - planck_bol[:, None, idxs]

        dpercoadd = np.zeros([nrecs, 5, nfreq], dtype=complex)
        dpercoadd[:, :4, idxs] = fact[None, :4, :] * planck_derivs[:, :4, idxs]
        derivs['percoadd'] = dpercoadd[:, :, idxs]

        derivs['element'] = fact[None, :, :] * planck_vals[:, :, idxs]
 
    return pspecs, power

def harm_vib_correction(vspec_rec, B, config, consts, tau, idxs):
    '''Simple function to calculate the values that multiple the harmonic and vibration
    correction parameters
    '''
    
    vspec = vspec_rec['spec_data']['spec']

    #Large non-zero value can mess with vibration correction. This value is not even
    #sent to the Fortran FISH, so I am just zeroing it out.
    vspec[0] = 0

    chan = vspec_rec['coad_spec_data']['chan_id']
    mtm_speed = vspec_rec['coad_spec_data']['mtm_speed']
    mtm_length = vspec_rec['coad_spec_data']['mtm_length']

    #iv should go (RHSS, RHSF, RHLF, RLSS, RLFS/FL, RLLF, LHSS, LHSF, LHLF, LLSS, LLFS/FL, LLLF)
    if mtm_length == 2:
        iv = 3*chan + 1
    else:
        iv = 3*chan + mtm_speed + mtm_length #only 3 input, speed/length=0, speed=X/length=2, speed=1/length=1 
    
    #Calculate the frequencies required for the 2nd and 3rd harmonic corrections
    #I think these indices would match up with the ones in FISH. Shifting by 1 seems
    #to be <1% change in the calculated spectra
    #idx_tmp = idxs
    idx_tmp = idxs - 1

    if2a = idx_tmp // 2
    if2b = (idx_tmp + 1) // 2

    if3a = idx_tmp // 3
    if3b = (idx_tmp + 1) // 3
    if3c = (idx_tmp + 2) // 3

    harm_corr_3 = (2.0 + B[idxs])/9.0 * (vspec[if3a] + vspec[if3b] + vspec[if3c])
    harm_corr_2 = (1.0 + B[idxs])/4.0 * (vspec[if2a] + vspec[if2b])
    
    #Calculate the frequencies required for computing the vibration correction detector delays 
    fp = idxs - config['vibcorr']['primary_offset'][iv]
    fs = idxs - config['vibcorr']['secondary_offset'][iv]
    
    ifp = np.array(np.abs(fp), dtype=int)
    ifs = np.array(np.abs(fs), dtype=int)

    Bvp = 1.0 + np.abs(fp)*consts['dw']*tau*1j
    Bvs = 1.0 + np.abs(fs)*consts['dw']*tau*1j
    
    prim_vib_corr = fp * Bvp * vspec[ifp]
    sec_vib_corr = fs * Bvs * vspec[ifs]

    return harm_corr_2, harm_corr_3, prim_vib_corr, sec_vib_corr

def harm_vib_correction2(vspec_recs, B, config, consts, tau, idxs, derivs=None):

    nrecs = len(vspec_recs)

    vspecs = vspec_recs['spec_data']['spec']
    
    #Large non-zero value can mess with vibration correction. This value is not even
    #sent to the Fortran FISH, so I am just zeroing it out.
    vspecs[:, 0] = 0
    
    chan = vspec_recs['coad_spec_data']['chan_id']
    mtm_speed = vspec_recs['coad_spec_data']['mtm_speed']
    mtm_length = vspec_recs['coad_spec_data']['mtm_length']
    
    #iv should go (RHSS, RHSF, RHLF, RLSS, RLFS/FL, RLLF, LHSS, LHSF, LHLF, LLSS, LLFS/FL, LLLF)
    iv = 3*chan + mtm_speed + mtm_length #only 3 input, speed/length=0, speed=0/length=1, speed=1/length=1

    #Calculate the frequencies required for the 2nd and 3rd harmonic corrections
    #I think these indices would match up with the ones in FISH. Shifting by 1 seems
    #to be <1% change in the calculated spectra
    idx_tmp = idxs

    if2a = idx_tmp // 2
    if2b = (idx_tmp + 1) // 2

    if3a = idx_tmp // 3
    if3b = (idx_tmp + 1) // 3
    if3c = (idx_tmp + 2) // 3
   
    #sf_newd.for has an if statement for when f=1 or 2, but that should never happen so we
    #are ignoring the if statement and only writing down the equation for f>=3
    harm_corr_3 = (2.0 + B[:, idxs])/9.0 * (vspecs[:, if3a] + vspecs[:, if3b] + vspecs[:, if3c])
    harm_corr_2 = (1.0 + B[:, idxs])/4.0 * (vspecs[:, if2a] + vspecs[:, if2b])

    #Calculate the frequencies required for computing the vibration correction detector delays
    fp = idxs[None, :] - config['vibcorr']['primary_offset'][iv][:, None]
    fs = idxs[None, :] - config['vibcorr']['secondary_offset'][iv][:, None]

    ifp = np.array(np.abs(fp), dtype=int)
    ifs = np.array(np.abs(fs), dtype=int)

    Bvp = 1.0 + np.abs(fp)*consts['dw']*tau[:, None]*1j
    Bvs = 1.0 + np.abs(fs)*consts['dw']*tau[:, None]*1j

    #prim_vib_corr = fp * Bvp * vspecs[:, idxs][ifp]
    #sec_vib_corr = fs * Bvs * vspecs[:, idxs][ifs]
    prim_vib_corr = fp * Bvp
    sec_vib_corr = fs * Bvs
    for i in range(nrecs):
        idx1 = ifp[i, :]
        idx2 = ifs[i, :]
        #prim_vib_corr[i, :] *= vspecs[i, idx1]
        #sec_vib_corr[i, :] *= vspecs[i, idx2]

        #NOTE: This seems to be the rough equation used in sf_newd.for. Doesn't
        #seem to make a significant difference
        prim_vib_corr[i, :] *= vspecs[i, idx1+1] * (np.abs(fp[i, :]) - ifp[i, :]) + \
                               vspecs[i, idx1] * (ifp[i, :] + 1 - np.abs(fp[i, :]))
        sec_vib_corr[i, :] *= vspecs[i, idx2+1] * (np.abs(fs[i, :]) - ifs[i, :]) + \
                              vspecs[i, idx2] * (ifs[i, :] + 1 - np.abs(fs[i, :]))

    if derivs is not None:
        derivs['Bvp'] = Bvp
        derivs['Bvs'] = Bvs

    return harm_corr_2, harm_corr_3, prim_vib_corr, sec_vib_corr

def autophase_correct(spec, pspec, variances, idxs):
    '''This function computes a linear autophase corrector for the calibrated spectrum
    '''

    alpha = 1.0e6

    #Compute the autophase corrector
    f = idxs

    x = np.sum(f*f*np.real(spec[idxs])**2 / variances[idxs])
    y = np.sum(f*np.imag(spec[idxs]-pspec[idxs]) * np.real(spec[idxs]) / variances[idxs])

    phase_corr = - y / (x + alpha)

    return phase_corr

def autophase_correct2(specs, pspecs, variances, idxs):
    '''This function computes a linear autophase corrector for the calibrated spectrum
    '''
    
    alpha = 1.0e6

    #Compute the autophase corrector
    f = idxs

    x = np.sum(f[None, :]*f[None, :]*np.real(specs[:, idxs])**2 / variances[:, idxs], axis=1)
    y = np.sum(f[None, :]*np.imag(specs[:, idxs]- pspecs[:, idxs]) * np.real(specs[:, idxs]) /
               variances[:, idxs], axis=1)

    phase_corr = - y / (x + alpha)

    return phase_corr

def pcalib_variances(vvar, emiss, S0, B, idxs):
    '''This function partially applies the FIRAS calibration model solution to the voltage
    variances, producing calibrated variances with the units of (ergs/sec/cm**2/sr/icm)**2.
    The input parameters are the voltage variances, the detector responsivity and delay,
    and the instrument transfer function. The real-real, imaginary-imaginary, and
    real-imaginary variances are all calibrated.
    '''
    
    cvar = np.zeros_like(vvar)

    #Calculate the calibration normalization and phase shifts
    xfer = B[idxs] / emiss[0, idxs] / S0
    xfer2 = np.real(xfer * xfer.conjugate())
    pshift = xfer / np.abs(xfer)
    rphase = np.real(pshift)
    rphase2 = rphase**2
    iphase = - np.imag(pshift)
    iphase2 = iphase**2

    #Calculate the variances (real-real, imag-imag, and real-imag)
    cvar[0, idxs] = (rphase2*vvar[0, idxs] - 2.0*rphase*iphase*vvar[2, idxs] + iphase2*vvar[1, idxs]) * xfer2
    cvar[1, idxs] = (rphase2*vvar[1, idxs] + 2.0*rphase*iphase*vvar[2, idxs] + iphase2*vvar[0, idxs]) * xfer2
    cvar[2, idxs] = ((iphase2-rphase2)*vvar[2, idxs] - rphase*iphase*(vvar[0, idxs] - vvar[1, idxs])) * xfer2

    return cvar

def pcalib_variances2(vvars, emiss, S0, B, idxs):

    cvars = np.zeros_like(vvars)

    xfer = B[:, idxs] / emiss[None, 0, idxs] / S0[:, None]
    xfer2 = np.real(xfer * xfer.conjugate())
    pshift = xfer / np.abs(xfer)
    rphase = np.real(pshift)
    rphase2 = rphase**2
    iphase = - np.imag(pshift)
    iphase2 = iphase**2
    
    #Calculate the variances (real-real, imag-imag, and real-imag)
    cvars[:, 0, idxs] = (rphase2*vvars[:, 0, idxs] - 2.0*rphase*iphase*vvars[:, 2, idxs] +
                         iphase2*vvars[:, 1, idxs]) * xfer2
    cvars[:, 1, idxs] = (rphase2*vvars[:, 1, idxs] + 2.0*rphase*iphase*vvars[:, 2, idxs] +
                         iphase2*vvars[:, 0, idxs]) * xfer2
    cvars[:, 2, idxs] = ((iphase2-rphase2)*vvars[:, 2, idxs] -
                         rphase*iphase*(vvars[:, 0, idxs] - vvars[:, 1, idxs])) * xfer2
 
    return cvars

def calc_responsivity(bol_volt, cmd_bias, Tdet, bolparam, derivs=None):
    '''This function calcualtes the responsivity and time constant for the FIRAS bolometers from
    the linear part of Steve Meyer's detector model. The input parameters are the detector
    voltage, commanded bias, and measured temperature for a given voltage spectrum and the
    derived bolometer parameters from the FIRAS calibration solution. The output parameters
    are the actual bolometer temperature, the IR power incident on the bolometer (in Watts), 
    the detector responsivity (in volts/erg), and the detector time constant for the spectrum

    Notes
    -----
    derivs is used to calculate information needed to calculate the derivative of the chi^2
    function
    '''

    delta_t = 1.0e-6
    #rscale = 1.0e-7

    #Extract the bolometer parameters from the array bolparm
    R0 = bolparam[0] #detector resistance at infinite temperature
    T0 = bolparam[1]#*100.0 #characteristic temperature for detector resistance function
    G1 = bolparam[2]#*1.0e-8 #coefficient for detector thermal conductance
    beta = bolparam[3] #index of temperature dependence of detector
    rho = bolparam[4] #electric field dependence of detector resistance
    #rho = np.abs(bolparam[4]) #electric field dependence of detector resistance
    C3 = bolparam[5]#*1.0e-10 #coefficent of cubic heat capacity term
    C1 = bolparam[6]#*1.0e-10 #coefficient of linear heat capacity term
    Jo = bolparam[7] #JFET offset
    Jg = bolparam[8] #JFET gain

    #This parameter was not fit and its value is not stored in the calibration files
    #RL = bolparam[9] #detector load resistance = 4.0e7
    RL = 4.0e7
  
    bol_volt = np.atleast_1d(bol_volt)
    cmd_bias = np.atleast_1d(cmd_bias)
    Tdet = np.atleast_1d(Tdet)

    ncoadds = len(bol_volt)

    if R0 < 0:
        raise ValueError("R0 must be positive")

    #if rho < 0:
    #    raise ValueError("rho must be positive")

    #Find the bolometer state
    V = (bol_volt - Jo) / Jg
    R = RL*V / (cmd_bias - V)
    X = V*rho
    Y = R/R0/X
    Tbol = Tdet
    
    #Figure out some derivatives
    niter = 8
    if derivs is not None:
        def qq(j, k):
            #TODO: CHECK
            #qq_val = np.max([1-np.abs(j-k), 0])/bolparm[j]
            qq_val = (j==k)/bolparam[j]
            return qq_val
       
        #Setup all the different arrays we need
        gV = np.zeros([ncoadds, 9])
        gR = np.zeros([ncoadds, 9])
        gQ = np.zeros([ncoadds, 9])
        dSQ = np.zeros([ncoadds, 9])
        dTT = np.zeros([ncoadds, niter+1, 9]) #so dTT[:, -1, :] is always 0
        for k in range(9):
            gV[:, k] = -(Jo*qq(k, 7)/Jg + V*qq(k, 8))
            gR[:, k] = R*(1/V + 1/(cmd_bias - V))*gV[:, k]
            gQ[:, k] = rho*qq(k, 4)*V + rho*gV[:, k]

    #Iterate to find the bolometer detector temperature and derivatives
    #for k in range(8):
    for k in range(niter):
        SQ = np.log(Y*Tbol*np.sinh(X/Tbol))
        #Tbol = T0 / np.log(Y*Tbol*np.sinh(X/Tbol))**2

        if derivs is not None:
            for j in range(9):
                dSQ[:, j] = gR[:, j]/R + dTT[:, k-1, j]/Tbol - gQ[:, j]/X - qq(j, 0) + \
                            (1.0/np.tanh(X/Tbol))*(gQ[:, j]/Tbol - dTT[:, k-1, j]*X/Tbol**2)
        
        Tbol = T0 / SQ**2

        if derivs is not None:
            for j in range(9):
                dTT[:, k, j] = Tbol*(qq(j, 1) - 2*dSQ[:, j]/SQ)

    Tconv = Tbol
    SQ = np.log(Y*Tbol*np.sinh(X/Tbol))
    Tbol = T0 / (SQ*SQ)
    #Tbol = T0 / np.log(Y*Tbol*np.sinh(X/Tbol))**2

    if np.any(np.abs(Tbol-Tconv) >= delta_t):
        dt0 = np.abs(Tbol-Tconv)
        idx0 = np.where(dt0 >= delta_t)
        raise ValueError("Couldn't iterate to detector temperature", idx0[0], dt0[idx0], delta_t)

    #Find the IR power incident on the detector
    BP = V * (cmd_bias - V) / RL
    P = (G1/(beta+1.0)) * (Tbol**(beta+1.0) - Tdet**(beta+1.0))
    Qrad = P - BP

    #C = C1*TBol**alpha
    #G = G1*Tbol + G3*Tbol**3
    
    #Find the detector responsivity and time constants
    H = Tbol / X * np.tanh(X/Tbol)
    G = G1*Tbol**beta
    C = C3*Tbol**3 + C1*Tbol

    DT = 1.0/H - 1.0 - 0.5*np.sqrt(T0/Tbol) #dln(R)/dln(T) - dln(R)/dln(V)
    Z = (G*Tbol*R + DT*V**2) / (G*Tbol*R/H - DT*V**2)
    S0 = rscale * R * (Z-H) / (V * (Z*R/RL + 1.0) * (H + 1.0))
    tau = C/G * (Z+1.0) * (R*H + RL) / ((Z*R + RL) * (H + 1.0))
    #avg_tau = np.average(tau, weights=adj_nifgs)

    '''
    #Calculations from sf16a.for
    if derivs is not None: 
        def qq(j, k):
            #qq_val = np.max([1-np.abs(j-k), 0])/bolparm[j]
            qq_val = (j==k)/bolparam[j]
            return qq_val

        #NOTE: gV = (dV/dx_i)/V
        gTs = np.zeros([ncoadds, 9])
        gSs = np.zeros([ncoadds, 9])
        gtaus = np.zeros([ncoadds, 9])
        for j in range(9):
            gV = -qq(j, 8)-qq(j, 7)*Jo/V/Jg
            gR = gV*(1+R/RL)
            gX = gV + qq(j, 4)
            gT = (gR-qq(j, 0)-gX+gX/H-0.5*qq(j, 1)*SQ)/(1/H-1-0.5*SQ)
            gH = (1-1.0/np.cosh(X/Tbol)**2/H)*(gT-gX)
            gDT = (-gH/H + 0.25*SQ*(gT-qq(j, 1)))/DT
            gG = qq(j, 2) + (gT+np.log(Tbol)*qq(j, 3))*beta
            gZ = (G*Tbol*R*(gG+gT+gR)+V*V*DT*(2*gV+gDT))/(G*Tbol*R+V*V*DT) - \
                 (G*Tbol*R/H*(gG+gT+gR-gH)-V*V*DT*(2*gV+gDT)) / (G*Tbol*R/H-V*V*DT)
            gS = gR+(Z*gZ-H*gH)/(Z-H)-gV-(gZ+gR)*(Z*R/RL)/(Z*R/RL+1)-H*gG/(H+1)
            gtau = gT+((qq(j, 5)+2*gT)*C3*Tbol*Tbol+qq(j, 6)*C1) / \
		   (C3*Tbol*Tbol+C1)-gG+Z*gZ/(Z+1)+(gR+gH)*R*H/(R*H+RL)- \
		   (gZ+gR)*Z*R/(Z*R+RL)-H*gH/(H+1)
            gTs[:, j] = gT
            gSs[:, j] = gS
            gtaus[:, j] = gtau
            #Qrad is not used in calibration so we don't need to calculate its derivative

        derivs['gT'] = gTs
        derivs['gS'] = gSs
        derivs['gtau'] = gtaus
    '''

    if derivs is not None:
        #Calculations from sf_newd.for
        gT = np.zeros([ncoadds, 9])
        gS = np.zeros([ncoadds, 9])
        gtau = np.zeros([ncoadds, 9])
        gtau2 = np.zeros([ncoadds, 9])
        for j in range(9):
            #These are NOT logarithmic derivs
            gT[:, j] = dTT[:, 7, j]
            gH = (gT[:, j]/Tbol-gQ[:, j]/X)*(H-1.0/np.cosh(X/Tbol)**2)
            gDT = -gH/H**2 - 0.25*np.sqrt(T0/Tbol)*(qq(j, 1)-gT[:, j]/Tbol)
            gG = G*(qq(j, 2) + beta*qq(j, 3)*np.log(Tbol) + beta*gT[:, j]/Tbol)
            gC = C3*qq(j, 5)*Tbol**3 + C1*qq(j, 6) + (C1 + 3*C3*Tbol**2)*gT[:, j]
            gZ = Z*(((gG*Tbol*R + G*R*gT[:, j] + G*Tbol*gR[:, j] + gDT*V**2 + \
                 2*DT*V*gV[:, j])/(G*Tbol*R+DT*V**2)) - (gG*Tbol*R/H + G*gT[:, j]*R/H + \
                 G*Tbol*gR[:, j]/H - 2*DT*V*gV[:, j] - gDT*V**2 - G*Tbol*R*gH/(H**2)) \
                 /(G*Tbol*R/H - DT*V**2))
	    #These are logarithmic. 
            #First equations for each derivative that is commented out is version in sf_newd.for
            #gS[:, j] = (gR[:, j]/R + (gZ-gH)/(Z-H) - gV[:, j]/V \
            #           -(gZ*R/RL+Z*gR[:, j]/RL)/(1+Z*R/RL))
            gS[:, j] = (gR[:, j]/R + (gZ-gH)/(Z-H) - gV[:, j]/V \
                       -(gZ*R/RL+Z*gR[:, j]/RL)/(1+Z*R/RL)) - gH/(H+1.)
            #gtau[:, j] = (gT[:, j]/Tbol + (C3*qq(5, j)*Tbol**2 + \
            #             2*C3*Tbol*gT[:, j] + C1*qq(6, j))/(C1*C3*Tbol**2) - \
            #             gG/G + gZ/(Z+1.) + (gR[:, j]*H+R*gH)/(R*H + RL) - \
            #             (gZ*R + Z*gR[:, j])/(Z*R+RL) - gH/(H+1.))
            gtau[:, j] = (gT[:, j]/Tbol + (C3*qq(5, j)*Tbol**2 + \
                         2*C3*Tbol*gT[:, j] + C1*qq(6, j))/(C1+C3*Tbol**2) - \
                         gG/G + gZ/(Z+1.) + (gR[:, j]*H+R*gH)/(R*H + RL) - \
                         (gZ*R + Z*gR[:, j])/(Z*R+RL) - gH/(H+1.))
            #gtau[:, j] = (gC/C - gG/G + gZ/(Z+1.) + (gR[:, j]*H+R*gH)/(R*H + RL) - \
            #             (gZ*R + Z*gR[:, j])/(Z*R+RL) - gH/(H+1.))

        #These derivatives should be in the same format as the sf16a derivatives (logarithmic)
        derivs['gT'] = gT
        derivs['gS'] = gS
        derivs['gtau'] = gtau

    return Tbol, Qrad, S0, tau

def temporal_drift(spec_time, param, timesinceapco=False):
    '''This function computes the time in years since the aperture cover ejection for a 
    spectrum. It then calculates the ICAL temporal drift correction and the primary (time
    dependent) vibration correction for the spectrum
    '''

    #ntimes = len(spec_time)

    if timesinceapco:
        time = spec_time
    else:
        time = util.time_since_apco(spec_time)

    #Calculate the ICAL temporal drift correction
    ticald = param[11] * np.exp(-param[12]*time) - param[13]

    #Calculate the primary (time dependent) vibration correction
    primary_vib = 0.0
    for j in range(4, -1, -1):
        primary_vib = primary_vib * time + param[15 + j]

    return primary_vib, ticald

def temporal_drift2(spec_times, param, timesinceapco=False, derivs=None):

    if type(spec_times) is int:
        spec_times = np.array([spec_times], dtype=int)

    if timesinceapco:
        times = spec_times
    else:
        times = np.zeros_like(spec_times, dtype=float)
        for i, spec_time in enumerate(spec_times):
            times[i] = util.time_since_apco(spec_time)
    
    #Calculate the ICAL temporal drift correction
    ticald = param[11] * np.exp(-param[12]*times) - param[13]

    #Calculate the primary (time dependent) vibration correction
    primary_vib = np.zeros_like(times)
    for j in range(4, -1, -1):
        primary_vib = primary_vib * times + param[15 + j]

    if len(primary_vib) == 1:
        primary_vib = primary_vib[0]
        ticald = ticald[0]
    
    if derivs is not None:
        nspec = len(spec_times)
        dvib = np.zeros([nspec, 5])
        dticald = np.zeros([nspec, 3])

        dticald[:, 0] = np.exp(-param[12]*times)
        dticald[:, 1] = -param[11]*times*np.exp(-param[12]*times)
        dticald[:, 2] = -1

        dvib[:, 0] = 1
        dvib[:, 1] = times
        dvib[:, 2] = times**2
        dvib[:, 3] = times**3
        dvib[:, 4] = times**4

        derivs['dticald'] = dticald
        derivs['dvib'] = dvib
    
    return primary_vib, ticald

def compute_constants(chan: int, scan_mode: int, config, fft_len: int = None, nyquist='default'):

    fcc_spec_length = fac_spec_length[scan_mode]

    #Get the nyquist frequencies and frequency intervals in icm and Hz.
    frec = 4*(chan % 2) + scan_mode
    fnyq_icm = config['nyquistl']['icm'][frec]
    spec_len = fcc_spec_length - 1.0

    #NOTE: Applying fac_nyq_correct is not done in Fortran code, but I think
    #it is correct based on published data
    if nyquist is None or nyquist == 'None':
        print("NO NYQUIST CORRECTION FACTOR")
        nyq_fact = 1
    elif nyquist == 'fish':
        print("USING 1.00275 AS NYQUIST CORRECTION FACTOR")
        nyq_fact = 1.00275
    elif nyquist == 'default':
        print("USING", fac_nyq_correct, "AS NYQUIST CORRECTION FACTOR")
        nyq_fact = fac_nyq_correct

    df = fnyq_icm / spec_len * nyq_fact
    fnyq_hz = config['nyquistl']['hz'][frec]
    dw = 2.0 * np.pi * fnyq_hz / spec_len

    #Fill in the frequency arrays
    freq = np.arange(fcc_spec_length)*df
    afreq = np.arange(fcc_spec_length)*dw

    #Calculate the phase shift required to make the high frequency short fast
    #spectra compatable with the high frequency long fast spectra
    if (chan == 0 or chan == 2) and scan_mode == 1:
        phase_shift = np.pi / fft_len*1j
    else:
        phase_shift = 0.0

    consts = {}
    consts['df'] = df
    consts['dw'] = dw
    consts['freq'] = freq
    consts['afreq'] = afreq
    consts['phase_shift'] = phase_shift
    #consts['spec_norm'] = spec_norm

    return consts

def load_model(fn_model, config=None, chan=0, scan_mode=0, norm_emiss=True):

    fn_split = os.path.splitext(fn_model)

    if fn_split[1].lower() == '.fits':
        hdulist = fits.open(fn_model)
        #model = np.array(0, dtype=dt.dt_fex_mod)
        model = {}
        model['bolparm'] = np.zeros(32)
        model['bolparm'][0] = hdulist[1].data['bolparm_'] #R0
        model['bolparm'][1] = hdulist[1].data['bolparm2'] #T0
        model['bolparm'][2] = hdulist[1].data['bolparm3'] #G1
        model['bolparm'][3] = hdulist[1].data['bolparm4'] #beta
        model['bolparm'][4] = hdulist[1].data['bolparm5'] #rho
        model['bolparm'][5] = hdulist[1].data['bolparm7'] #C3
        model['bolparm'][6] = hdulist[1].data['bolparm6'] #C1
        model['bolparm'][7] = hdulist[1].data['bolparm8'] #Jo
        model['bolparm'][8] = hdulist[1].data['bolparm9'] #Jg
        model['bolparm'][9] = hdulist[1].data['optparm2']
        model['bolparm'][10] = hdulist[1].data['optparm_']
        model['bolparm'][11] = hdulist[1].data['drift_am']
        model['bolparm'][12] = 1.0 / hdulist[1].data['drift_tc']
        model['bolparm'][13] = hdulist[1].data['drift_of']
        model['bolparm'][14] = hdulist[1].data['optparm8'] #sec vib off
        model['bolparm'][15] = hdulist[1].data['optparm3'] #prim vib off
        model['bolparm'][16] = hdulist[1].data['optparm4'] #prim vib lin
        model['bolparm'][17] = hdulist[1].data['optparm5'] #prim vib quad
        model['bolparm'][18] = hdulist[1].data['optparm6'] #prim vib cub
        model['bolparm'][19] = hdulist[1].data['optparm7'] #prim vib quart
        model['bolparm'][20] = hdulist[1].data['xcal_cor'] #xcal correction

        #transfer, ical, skyhorn, refhorn, dihedral, struct, bolometer
        npts = 257

        #Based on FISH output, emissivities start at index = 3 (4 in Fortran for FISH)
        nemiss = np.size(hdulist[1].data['rical'])
        
        #Determine starting index of first frequency in published data files
        i0 = np.round(hdulist[0].header['NU_ZERO'] / hdulist[0].header['DELTA_NU']).astype(int)

        model['emissivity'] = np.zeros([7, npts], dtype=complex)
        model['emissivity'][0, i0:i0+nemiss] = hdulist[1].data['rtransfe'] + hdulist[1].data['itransfe']*1j
        model['emissivity'][1, i0:i0+nemiss] = hdulist[1].data['rical'] + hdulist[1].data['iical']*1j
        model['emissivity'][2, i0:i0+nemiss] = hdulist[1].data['rskyhorn'] + hdulist[1].data['iskyhorn']*1j
        model['emissivity'][3, i0:i0+nemiss] = hdulist[1].data['rrefhorn'] + hdulist[1].data['irefhorn']*1j
        model['emissivity'][4, i0:i0+nemiss] = hdulist[1].data['rdihedra'] + hdulist[1].data['idihedra']*1j
        model['emissivity'][5, i0:i0+nemiss] = hdulist[1].data['rstructu'] + hdulist[1].data['istructu']*1j
        model['emissivity'][6, i0:i0+nemiss] = hdulist[1].data['rbolomet'] + hdulist[1].data['ibolomet']*1j

        #model['emissivity'].T.tofile('emissivity_llss.dat')
        #zzz

        hdulist.close()
    elif fn_split[1].lower() == '.hdf5' or fn_split[1] == '.h5':
        h5file = h5py.File(fn_model, 'r')

        model = {}
        model['bolparm'] = h5file['bolparm'][:]
        model['emissivity'] = h5file['emissivity'][:]

    elif fn_split[1].lower() == '.txt':
        pass

    elif fn_split[1].lower() == '.txt_pass4':
        data = np.loadtxt(fn_model, skiprows = 3)
        npts = np.size(data)
        data.shape = (npts,)
        renormalize = np.array([1.0e+00, 1.0e+02, 1.0e-08, 1.0e+00,
                                1.0e+00, 1.0e-10, 1.0e-10, 1.0e+00,
                                1.0e+00, 1.0e-03, 1.0e-03, 1.0e-03,
                                1.0e+00, 1.0e-03, 1.0e-04, 1.0e-04,
                                1.0e-04, 1.0e-04, 1.0e-04, 1.0e-04,
                                1.0e-03, 1.0e+00, 1.0e+00, 1.0e+00,
                                1.0e+00, 1.0e+00, 1.0e+00, 1.0e+00,
                                1.0e+00, 1.0e+00, 1.0e+00, 1.0e+00])

        model = {}
        model['bolparm'] = data[:32] * renormalize
        
        if norm_emiss:
            frec = 4*(chan % 2) + scan_mode
            etendu = 1.5
            fnyq = config['nyquistl']['icm'][frec]
            if chan <= 1:
                norm_fact = - fnyq * etendu
            else:
                norm_fact = fnyq * etendu
        else:
            norm_fact = 1

        data_emiss = data[32:]

        data_emiss.dtype = complex
        
        #data_emiss.shape = (6, -1)
        data_emiss.shape = (-1, 6)
        data_emiss = data_emiss.T
        
        #transfer, ical, skyhorn, refhorn, dihedral, struct, bolometer
        npts = 257

        #Based on FISH output, emissivities start at index = 3 (4 in Fortran for FISH)
        nemiss = data_emiss.shape[1]
        model['emissivity'] = np.zeros([7, npts], dtype=complex)
        model['emissivity'][:6, 1:nemiss+1] = data_emiss / norm_fact
        model['emissivity'][6, :] = - np.sum(model['emissivity'][:6, :], axis=0)

    else:
        raise ValueError("Unknown file type for input model")
        
    model['transfer'] = model['emissivity'][0, :]
    model['ical'] = model['emissivity'][1, :]
    model['skyhorn'] = model['emissivity'][2, :]
    model['refhorn'] = model['emissivity'][3, :]
    model['dihedral'] = model['emissivity'][4, :]
    model['struct'] = model['emissivity'][5, :]
    model['bolometer'] = model['emissivity'][6, :]

    return model

def read_coadd(fn, chan, scan_mode, config, input_type='sky', quality=None):

    file_ext = os.path.splitext(fn)[1].lower()
    if file_ext == '.h5':
        h5file = h5py.File(fn, 'r')
        key = 'fil_' + input_type + '_' + chan + scan_mode
        coadd_recs = h5file[key][:]
        coadd_recs = coadd_recs.view(dt.dt_fil_sky)
    elif file_ext == '.fits':
        #Copy of code in FISH

        #Read in the published FITS data and stick it in the appropriate data
        #structure
        hdulist = fits.open(fn)
        data_fits = hdulist[1].data
        ndata = len(data_fits['NUM_IFGS'])
        data_all = np.zeros(ndata, dtype=dt.dt_fil_sky)
        
        #Determine starting index of first frequency in published data files
        i0 = np.round(hdulist[0].header['NU_ZERO'] / hdulist[0].header['DELTA_NU']).astype(int)
        
        #Not sure these pointing values mean anything since this is more calibration data
        #but they are stored in the FITS file so I am adding them to the record
        data_all['attitude']['pixel_no'] = data_fits['PIXEL']
        data_all['attitude']['pixel_definition'] = 'q'
        data_all['attitude']['skymap_index'] = 6

        lon, lat = util.pointing_conversion(data_fits['GAL_LON'], data_fits['GAL_LAT'])        
        data_all['attitude']['galactic_longitude'] = lon
        data_all['attitude']['galactic_latitude'] = lat
        lon, lat = util.pointing_conversion(data_fits['ECL_LON'], data_fits['ECL_LAT'])        
        data_all['attitude']['ecliptic_longitude'] = lon
        data_all['attitude']['ecliptic_latitude'] = lat
        lon, lat = util.pointing_conversion(data_fits['RA'], data_fits['DEC'])        
        data_all['attitude']['ra'] = lon
        data_all['attitude']['dec'] = lat

        data_all['coad_data']['ifg'] = data_fits['COADDED_']

        nvar = len(data_fits['REAL_VAR'][0, :])
    
        fact = fac_erg_to_mjy

        data_all['coad_data']['real_var'][:, i0:nvar+i0] = data_fits['REAL_VAR'] / fact**2
        data_all['coad_data']['imag_var'][:, i0:nvar+i0] = data_fits['IMAG_VAR'] / fact**2
        data_all['coad_data']['real_imag_var'][:, i0:nvar+i0] = data_fits['REAL_IMA'] / fact**2
        data_all['coad_spec_head']['num_ifgs'] = data_fits['NUM_IFGS']
        data_all['coad_spec_data']['glitch_rate'] = data_fits['GLITCH_R']
        data_all['coad_spec_head']['adj_num_ifgs'] = data_fits['ADJ_NUM_']
        
        #need to calculate time from the time in the FITS file
        time_fits = data_fits['TIME']
        time_ifg_fits2 = data_fits['IFG_TIME']

        data_all['ct_head']['time'] = util.fitstime_to_bintime(time_fits)

        data_all['coad_spec_data']['orphans'] = data_fits['ORPHANS']

        # = data_fits['NUM_TEMP'] #number of templates subtracted from data
        data_all['coad_spec_data']['bol_volt'] = data_fits['BOLOM_VO']
        data_all['coad_spec_data']['bol_cmd_bias'] = data_fits['BOLOM_BI'] * 25.5

        data_all['coad_spec_data']['temp'][:, 0] = data_fits['XCAL_TEM']
        data_all['coad_spec_data']['temp'][:, 1] = data_fits['ICAL_TEM']
        data_all['coad_spec_data']['temp'][:, 2] = data_fits['SKYHORN_']
        data_all['coad_spec_data']['temp'][:, 3] = data_fits['REFHORN_']
        data_all['coad_spec_data']['temp'][:, 4] = data_fits['DIHEDRAL']
        data_all['coad_spec_data']['temp'][:, 5] = data_fits['MIRROR_T']
        data_all['coad_spec_data']['temp'][:, 6:] = data_fits['BATH_TEM'][:, np.newaxis]

        channel = util.determine_chan_num(chan)
        sm = util.determine_scan_mode(scan_mode.lower())

        if scan_mode.lower() == 'fl' or scan_mode.lower() == 'fs':
            data_all['coad_spec_data']['mtm_speed'] = 1
            data_all['coad_spec_data']['mtm_length'] = 0
        else: 
            data_all['coad_spec_data']['mtm_speed'] = sm % 2
            data_all['coad_spec_data']['mtm_length'] = sm // 2

        data_all['coad_data']['fft_length'] = hdulist[0].header['FFT_LENG']
        data_all['coad_spec_data']['peak_pos'] = hdulist[0].header['PEAK_POS'] - 1

        #This seems to be correct for LLSS
        data_all['coad_spec_data']['sci_mode'] = 4
        data_all['coad_spec_data']['fakeit'] = 0
        if scan_mode == 'ss':
            data_all['coad_spec_data']['adds_per_group'] = 3
        #elif scan_mode == 'lf' or scan_mode == 'sf' or scan_mode == 'fs' or scan_mode == 'fl':
        else:
            if chan == 'll' or chan == 'rl':
                data_all['coad_spec_data']['adds_per_group'] = 8
            elif chan == 'rh' or chan == 'lh':
                data_all['coad_spec_data']['adds_per_group'] = 2

        index = 4*(channel % 2) + sm
        nyquist = config['nyquistl']
        data_all['coad_spec_data']['nyquist_icm'] = nyquist['icm'][index]
        data_all['coad_spec_data']['nyquist_hertz'] = nyquist['hz'][index]
        data_all['coad_spec_data']['chan_id'] = channel

        #Remove coadds that seem to have zero uncertainty
        #idx_good = np.std(data_all['coad_data']['real_var'], axis=1) != 0
        #data_all = data_all[idx_good]

        coadd_recs = data_all

    return coadd_recs

def downsample_emiss(emissivity, chan, out_idx=False):
    '''The idea of this function is the same as the fil.short function, except in the frequency space. We want
    to downsample the emissivity so that a high resolution emissivity can be applied to lower resolution data
    (i.e. emissivity with resolution for LLLF can be applied to LLSS, LLFL, LLFS).
    '''

    #nyquistl = config['nyquistl']

    if chan.lower() in ['rh', 'lh']:
        #For LH or RH we don't need to do anything since all the data at the
        #high frequency is at the same resolution
        return emissivity
    elif chan.lower() in ['ll', 'rl']:
        idx0 = 19 #([-1, 0, 1, 2] average and NU_ZERO_LRES / DELTA_NU_HRES = 20
        emiss_out = np.zeros_like(emissivity)
        idx_arr = np.zeros_like(emiss_out[0, :], dtype=int)
        for i in range(43):
            idx1 = idx0 + 4 
            emiss_out[:6, i+5] = np.mean(emissivity[:6, idx0:idx1], axis=1)
            idx_arr[idx0:idx1] = i+5
            idx0 = idx1
        emiss_out[6, :] = -np.sum(emiss_out[:6, :], axis=0)
    else:
        raise ValueError("Bad value for chan:", chan)

    if out_idx:
        return emiss_out, idx_arr
    else:
        return emiss_out

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='FIRAS Spectra Long (replaces FIRAS Calibrate FIRAS)')

    parser.add_argument('-fil_path', dest='fil_path', action='store', help='Path to FIL files')
    parser.add_argument('-ref_path', dest='ref_path', action='store', help='Path to reference files')
    parser.add_argument('-out_path', dest='out_path', action='store', help='Path to save output data')
    parser.add_argument('-jstart', dest='jstart', action='store',
                        help='Input coadded IFGs can be specified by timerange')
    parser.add_argument('-jstop', dest='jstop', action='store', help='Input coadded IFGs can be specified by timerange')
    parser.add_argument('-input', dest='input', default='sky', action='store',
                        help='FSL will process either sky data or calibration data during a single invocation')
    parser.add_argument('-pixel', dest='pixel', nargs='*', type=int, action='store',
                        help='FSL will process the particular sky pixels specified')
    parser.add_argument('-quality', dest='quality', nargs='*', 
                        help='The data quality thresholds for the input coadded ifgs')
    parser.add_argument('-channel', dest='channel', action='store',
                        help='FSL process data for only one channel per invocation')
    parser.add_argument('-scan_mode', dest='scan_mode', action='store', 
                        help='FSL process data for only one scan mode per invocation')
    parser.add_argument('-single_ifg', dest='single_ifg', action='store_true',
                        help='If set, FSL process the spectra as single IFGs, whether or not they actually are single')
    parser.add_argument('-nocalibrate', dest='calibrate', action='store_false',
                        help='''If set, FSL will not calibrate voltage spectra. Will speed up code if producing
                        voltage spectra only''')
    parser.add_argument('-model_ext', dest='model_ext', action='store',
                        help='''The calibration model soltuion that FSL will apply. Input the model extension. This
                        is required unless -nocalibrate is set''')
    parser.add_argument('-dvector', dest='dvector', action='store_true',
                        help='''If this is set, FSL reads the variance estimate from the reference dataset FEX_VAR.
                        Otherwise, FSL uses the individual coadd varainces for the varaince estimates. Cannot be
                        used with -fil_var''')
    parser.add_argument('-differential', dest='differential', action='store_true',
                        help='''If this is set, FSL does not add the reference spectra to the differential spectra to
                        obtain the calibration spectra. Cannot be used in conjunction with -nocalibrate''')
    parser.add_argument('-fil_var', dest='fil_var', action='store_true',
                        help='''If set, FSL reads the variance estimate from the reference dataset FEX_FLV. Cannot be
                        used with -dvector''')
    parser.add_argument('-tsig_zero', dest='tsig_zero', action='store_true', help='''If set, FSL sets the combined
                        temperature sigmas to zero for all coadds before the computation of the reference spectra''')
    parser.add_argument('-write', dest='write', action='store', default='cal', help='''FSL can write voltage
                        and/or differential or calibrated spectra. Sky data are written to skymaps and calibration
                        spectra are written to TODs''')
    parser.add_argument('-make_plots', dest='make_plots', action='store_true', help='''Make and save plots of the
                        calibrated spectra''')

    args = parser.parse_args()

    if args.dvector and args.fil_var:
        raise ValueError('Only one of dvector and fil_var can be set')

    #If no input, we want to process all pixels
    if args.pixel is None:
        pixels = 'all'
    else:
        pixels = args.pixel

    #Process the spectra for the specified channel and scan mode
    if args.channel is None:
        raise ValueError('A channel must be input')
    else:
        chan = util.determine_chan_num(args.channel)

    if args.scan_mode is None:
        raise ValueError('A scan mode must be input')
    else:
        scan_mode = util.determine_scan_mode(args.scan_mode)

    if args.quality is None:
        if args.input == 'sky':
            quality = [3, 3]
        elif args.input == 'cal':
            quality = [3, 6]
    elif len(args.quality) == 2:
        quality = args.quality
    else:
        raise ValueError('quality must be a list of 2 values if input')
    
    if args.dvector:
        raise ValueError("dvector not implemented yet")

    if args.jstart is not None or args.jstop is not None:
        raise ValueError("jstart and jstop are not implemented")
    
    config = read_config(args.ref_path, grtcoawt=True, grttrans=True, samprate=True, nyquistl=True, vibcorr=True)

    if args.calibrate:
        model = load_model(args.model_ext, config=config, chan=chan, scan_mode=scan_mode)
    
    apodl_all = frd.apodl()
    etfl_all = frd.elex_transfcnl(config) #TODO: seems to take 3 seconds to run. Why???
   
    print("Reading coadds")

    #Read in a set of data for a particular scan mode. Calibration data are read sequentially. Sky
    #data are read by pixel number. Each call for a pixel may return several records for the pixel.
    coadd_recs = read_coadd(args.fil_path, args.channel, args.scan_mode, config, input_type=args.input, quality=quality)

    nrecords = len(coadd_recs)

    vspec_recs = np.zeros_like(coadd_recs, dtype=dt.dt_fsl_sky)
    cspec_recs = np.zeros_like(coadd_recs, dtype=dt.dt_fsl_sky)
    
    print("Calibrating Spectra:", nrecords)

    time0 = datetime.datetime.now()
    #TODO: Implement vectorized functions, but it only takes a few seconds to run and not necessarily needed
    #record_list = [467, 468, 469]
    #for j in record_list:
    for j in range(nrecords):
        if j % 1000 == 0:
            print(j, nrecords)

        #Produce the voltage spectra
        vspec_rec  = produce_spectra(coadd_recs[j], chan, scan_mode, apodl_all, etfl_all) 
        vspec_recs[j] = vspec_rec

        '''
        #Calibrate the spectra
        if args.calibrate:
            cspec_rec, cvar = calibrate_spectra(vspec_rec, model, config, chan, scan_mode, input_type=args.input,
                                                tsig0=args.tsig_zero, fcc_flv=args.fil_var, diff=args.differential,
                                                dvector=None)

            cspec_recs[j] = cspec_rec
        '''

        #Plot the spectra
        if args.make_plots:
            #TODO: Do stuff here
            pass
    
    #vspec_recs['spec_data']['spec'][:, 1:] = vtmp['vspec'][:, :] * (a0.real / b0.real)[:, None]
    
    if args.calibrate:
        cspec_recs, cvars = calibrate_spectra2(vspec_recs, model, config, chan, scan_mode, input_type=args.input,
                                               tsig0=args.tsig_zero, fcc_flv=args.fil_var, diff=args.differential,
                                               dvector=None)
   
    time1 = datetime.datetime.now()

    print("Time to calibrate:", time1-time0)

    #TODO: Change HDF5 file keys to channel+scanmode/name
    #Write the voltage spectra to an archive
    if args.write == 'volt' or args.write == 'both':
        if args.type == 'sky':
            fn_out = 'fsl_vsk_'
        elif args.type == 'cal':
            fn_out = 'fsl_vcl_'
        fn_out += args.channel + args.scan_mode
        fn_out = os.path.join(args.out_path, fn_out)

        vspec_recs.tofile(fn_out)

    #Write the calibrated or differential spectra to an archive
    if args.calibrate and (args.write == 'cal' or args.write == 'both'):
        if args.differential:
            if args.input == 'sky':
                fn_out = 'fsl_dsk_'
            elif args.input == 'cal':
                fn_out = 'fsl_dcl_'
        else:
            if args.input == 'sky':
                fn_out = 'fsl_sky_'
            elif args.input == 'cal':
                fn_out = 'fsl_cal_'
        fn_out += args.channel + args.scan_mode

        fn_out_h5 = args.out_path
        h5file = h5py.File(fn_out_h5, 'a')
        if fn_out in h5file.keys():
            print("Deleting:", fn_out)
            del h5file[fn_out]
        h5file[fn_out] = cspec_recs.view(dt.dt_fsl_sky2)
        h5file.close()

    #Write out the number of spectra processed
    if args.input == 'sky':
        print('Number of sky spectra processed:', nrecords)
        npix = len(np.unique(cspec_recs['attitude']['pixel_no']))
        print('Number of pixels containing spectra:', npix)
    elif args.input == 'cal':
        print('Number of calibration spectra processed:', nrecords)

