from glob import glob
import os

map_files = glob('/mn/stornext/d16/cmbco/ola/firas/healpix_maps/FIRAS_map_????GHz_*f.fits')
rms_files = glob('/mn/stornext/d16/cmbco/ola/firas/healpix_maps/FIRAS_rms_????GHz_*f.fits')

map_files.sort()
rms_files.sort()

print(len(map_files))
print(len(rms_files))


# Common

for i in range(len(map_files)):
    freq_str, high_low = map_files[i].split('_')[-2:]
    high_low = high_low.split('.fits')[0]
    freq = int(freq_str[:4])
    print(freq)
    highlow = 'H' if ('high' in high_low) else 'L'
    print(highlow, high_low)

    if freq < 480:
        bp_file = f'FIRAS_bp_{freq_str}_{high_low}_v2_trunc.dat'
    else:
        bp_file = f'FIRAS_bp_{freq_str}_{high_low}_v2.dat'

    if os.path.exists(f'FIRAS_{highlow}{freq:04}_map.defaults'):
        continue
    with open(f'FIRAS_{highlow}{freq:04}_map.defaults', 'w') as f:
        f.writelines([
                f'BAND_LABEL&&&                  = FIRAS_{highlow}{freq:04}\n',
                f'BAND_INSTRUMENT_LABEL&&&       = FIRAS_{highlow}{freq:04}\n',
                f'BAND_MAPFILE&&&                = FIRAS_map_{freq_str}_{high_low}.fits\n',
                f'BAND_NOISEFILE&&&              = FIRAS_rms_{freq_str}_{high_low}.fits\n',
                f'BAND_BANDPASSFILE&&&           = {bp_file}\n',
                f'BAND_NOMINAL_FREQ&&&           = {freq}\n',
                 'BAND_OBS_PERIOD&&&             = 1\n',
                 'BAND_POLARIZATION&&&           = .false.\n',
                 'BAND_NSIDE&&&                  = 16\n',
                 'BAND_LMAX&&&                   = 192\n',
                 'BAND_UNIT&&&                   = MJy/sr\n',
                 'BAND_NOISE_FORMAT&&&           = rms\n',
                 'BAND_NOISE_UNIFORMIZE_FSKY&&&  = 0.0 \n',
                 'BAND_MASKFILE&&&               = FIRAS_mask.fits\n',
                 'BAND_MASKFILE_CALIB&&&         = fullsky\n',
                 'BAND_BEAMTYPE&&&               = FIRAS \n',
                 'BAND_BEAM_B_L_FILE&&&          = beam_7deg.fits\n',
                 'BAND_BEAM_B_PTSRC_FILE&&&      = FIRAS_beam_spline.dat\n',
                 'BAND_PIXEL_WINDOW&&&           = pixel_window_n0016.fits\n',
                 'BAND_SAMP_NOISE_AMP&&&         = .false.\n',
                 'BAND_BANDPASS_TYPE&&&          = HFI_submm\n',
                 'BAND_BANDPASS_MODEL&&&         = additive_shift\n',
                 'BAND_SAMP_BANDPASS&&&          = .false.\n',
                 'BAND_SAMP_GAIN&&&              = .false.\n',
                 'BAND_GAIN_PRIOR_MEAN&&&         = 1.\n',
                 'BAND_GAIN_PRIOR_RMS&&&          = 0\n',
                 'BAND_GAIN_CALIB_COMP&&&        = all\n',
                 'BAND_GAIN_LMIN&&&              = -1\n',
                 'BAND_GAIN_LMAX&&&              = -1\n',
                 'BAND_GAIN_APOD_MASK&&&         = fullsky\n',
                 'BAND_GAIN_APOD_FWHM&&&         = 120.\n',
                 'BAND_DEFAULT_GAIN&&&           =   1.\n',
                 'BAND_DEFAULT_BP_DELTA&&&       =   0.\n',
                 'BAND_DEFAULT_NOISEAMP&&&       =   1.\n',
                 'BAND_COMPONENT_SENSITIVITY&&&  = broadband\n',
                 '\n',
                 'BAND_TOD_TYPE&&&               = none\n',
                 'BAND_REG_NOISEFILE&&&          = none\n',
                 'BAND_NOISE_RMS&&&_SMOOTH01     = none\n',
                 'BAND_NOISE_RMS&&&_SMOOTH02     = none\n',
                 'BAND_NOISE_RMS&&&_SMOOTH03     = none'])

