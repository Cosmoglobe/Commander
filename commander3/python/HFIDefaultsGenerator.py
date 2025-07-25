import os

### Function definitions
def checkFileExists(files, folder=''):
    """Check if a single file or a list of files exist in the specified folder."""
    if isinstance(files, str):
        files = [files]
    for file in files:
        if not os.path.isfile(f"{folder}{file}"):
            print(f"File {file} is not in {folder}")

def write2v(f, textv, val, dim, i, j):
    """Write formatted parameter values to a file."""
    if dim == 0:
        f.write(f"{textv.ljust(35)}= {val}\n")
    elif dim == 1:
        f.write(f"{textv.ljust(35)}= {val[i]}\n")
    elif dim == 2:
        f.write(f"{textv.ljust(35)}= {val[i][j]}\n")
    else:
        print("Error: Dim should be 1, 2, or 3")


# Parameters that are the same for all the HFI bands (across frequency and horns) should have 
# param_name_v=0 to tell the file generator that there is a single value for all bands, if 
# there are different values for different frequencies change to param_name_v=1 and if
# it is unique for each horn param_name_v=2
# This is used for 'dim' (dimension) in write2v


# Formatting parameters
data_folder = '/mn/stornext/d23/cmbco/cg/hfi/data/'
eqsig = "&&&"
debug = True  # Enable file existence checks

# Save folder for default files
homedir='/mn/stornext/u3/raelynsu/code/'
param_folder = homedir+"Commander/commander3/parameter_files/defaults/bands/HFI/"

# Frequency-related parameters
## BAND_NOMINAL_FREQ
nomfreqs = [100, 143, 217, 353, 545, 857]
nomfreqs_v = 1
##BAND_NSIDE
nsides = [2048, 2048, 4096, 4096, 4096, 4096]
nsides_v = 1
##BAND_LMAX
lmaxs = [6000, 6000, 10000, 10000, 10000, 10000]
lmaxs_v = 1


# File names for pixel window functions
band_pixel_window = [f"pixel_window_n{n}.fits" for n in nsides]
band_pixel_window_v = 1
    
if debug:
    checkFileExists(band_pixel_window,folder=data_folder) 
    # checkFileExists(band_beam_b_l_file,folder=data_folder) 


# Data paths for individual horns
common_folder = "/mn/stornext/d23/cmbco/globe/fullres/planck/"
maps_folder, rms_folder, hit_folder = "map/map_", "rms/rms_", "mask/hitmask_"
common_string3 = "_fullres_full_T.fits"

band_label_v=2
band_mapfile_v=2
band_noisefile_v=2
band_bandpassfile_v=2
band_maskfile_v=2

# Horn configurations
band_arr = [
    ['-1', '-2', '-3', '-4'],                           # 100 GHz
    ['-1', '-2', '-3', '-4', '-5', '-6', '-7'],         # 143 GHz
    ['-1', '-2', '-3', '-4', '-5', '-6', '-7', '-8'],   # 217 GHz
    ['-1', '-2', '-3', '-4', '-5', '-6', '-7', '-8'],   # 353 GHz
    ['-1', '-2', '-4'],                                 # 545 GHz
    ['-1', '-2', '-3', '-4']                            # 857 GHz
]


# Initialize band-specific lists
band_label, band_mapfile, band_noisefile, band_bandpassfile, band_maskfile,band_beam_b_l_file = [], [], [], [], [], []

##bandpass name changes for different horns
v1='_SHORN.dat'
v2='.dat'
bp_extra=[[v1,v1,v1,v1],
          [v1,v1,v1,v1,v2,v2,v2],
          [v2,v2,v2,v2,v1,v1,v1,v1],
          [v2,v2,v1,v1,v1,v1,v2,v2],
          [v2,v2,v2],
          [v2,v2,v2,v2]]

for i in range(len(nomfreqs)):
    b_lab = [str(nomfreqs[i]) + h for h in band_arr[i]]
    band_label.append(b_lab)

    common_string2 = f"_npipe_0{nsides[i]}"
    map_file = [f"{common_folder}{maps_folder}planck_{b}{common_string2}{common_string3}" for b in b_lab]
    rms_file = [f"{common_folder}{rms_folder}planck_{b}{common_string2}{common_string3}" for b in b_lab]
    hit_file = [f"{common_folder}{hit_folder}planck_{b}{common_string2}{common_string3}" for b in b_lab]
    band_mapfile.append(map_file)
    band_noisefile.append(rms_file)
    band_maskfile.append(hit_file)
    
    bpfile = [f"bp_RIMO_v2.0_{b}{bp}" for b, bp in zip(b_lab, bp_extra[i])]
    band_bandpassfile.append(bpfile)

    band_beam_b_l= [f"Bl_npipe6v20_{f}x{f}_extended.fits" for f in b_lab]
    band_beam_b_l_file.append(band_beam_b_l)

    if debug:
        checkFileExists(map_file)
        checkFileExists(rms_file)
        checkFileExists(hit_file)
        checkFileExists(bpfile,folder=data_folder)
        checkFileExists(band_beam_b_l_file,folder=data_folder)



band_instrument_label=band_label
band_instrument_label_v=2
band_beam_b_l_file_v = 2 

###constant for all HFI currently
band_obs_period='1' 
band_obs_period_v=0
band_polarization='.false.'
band_polarization_v=0

band_unit='uK_cmb'
band_unit_v=0
band_noise_format='rms'
band_noise_format_v=0
band_reg_noisefile='none'
band_reg_noisefile_v=0
band_noise_uniformize_fsky='0.01'
band_noise_uniformize_fsky_v=0
band_maskfile_calib='fullsky'
band_maskfile_calib_v=0
band_beamtype='b_l' # {b_l, febecop}
band_beamtype_v=0
band_beam_b_ptsrc_file='none'
band_beam_b_ptsrc_file_v=0
band_sample_noise_amp='.false.'
band_sample_noise_amp_v=0

band_bandpass_type='HFI_cmb'
band_bandpass_type_v=0
band_bandpass_model='additive_shift' # {powlaw_tilt, additive_shift}
band_bandpass_model_v=0

band_samp_bandpass='.false'
band_samp_bandpass_v=0
band_samp_gain='.false.'
band_samp_gain_v=0
band_gain_prior_mean='1.'
band_gain_prior_mean_v=0
band_gain_prior_rms='0.'
band_gain_prior_rms_v=0
band_gain_calib_comp='all' #all or cmb
band_gain_calib_comp_v=0
band_gain_lmin='25'
band_gain_lmin_v=0
band_gain_lmax='100'
band_gain_lmax_v=0
band_gain_apod_mask='fullsky'
band_gain_apod_mask_v=0
band_gain_apod_fwhm='120.'
band_gain_apod_fwhm_v=0
band_default_gain='1.'
band_default_gain_v=0
band_default_bp_delta='0.'
band_default_bp_delta_v=0
band_default_noiseamp='1.'
band_default_noiseamp_v=0
band_component_sensitivity='broadband'
band_component_sensitivity_v=0
band_tod_type='none'
band_tod_type_v=0
band_noise_rms_sm1='none'
band_noise_rms_sm1_v=0
band_noise_rms_sm2='none'
band_noise_rms_sm2_v=0
band_noise_rms_sm3='none'
band_noise_rms_sm3_v=0

# band=
# band_v=
# ''


for i in range(len(nomfreqs)):
    for j in range(len(band_arr[i])):
        f = open(f"{param_folder}HFI_{nomfreqs[i]}{band_arr[i][j]}_HFI.defaults",'w') ##will overwrite exiting files
        f.write(f"# {nomfreqs[i]}{band_arr[i][j]} GHz map parameters \n")
        write2v(f,'BAND_LABEL'+eqsig,band_label,band_label_v,i,j)
        write2v(f,'BAND_INSTRUMENT_LABEL'+eqsig,band_instrument_label,band_instrument_label_v,i,j)
        write2v(f,'BAND_OBS_PERIOD'+eqsig,band_obs_period,band_obs_period_v,i,j)
        write2v(f,'BAND_POLARIZATION'+eqsig,band_polarization,band_polarization_v,i,j)
        write2v(f,'BAND_NSIDE'+eqsig,nsides,nsides_v,i,j)
        write2v(f,'BAND_LMAX'+eqsig,lmaxs,lmaxs_v,i,j)
        write2v(f,'BAND_UNIT'+eqsig,band_unit,band_unit_v,i,j)
        write2v(f,'BAND_NOISE_FORMAT'+eqsig,band_noise_format,band_noise_format_v,i,j)
        write2v(f,'BAND_MAPFILE'+eqsig,band_mapfile,band_mapfile_v,i,j)
        write2v(f,'BAND_NOISEFILE'+eqsig,band_noisefile,band_noisefile_v,i,j)
        write2v(f,'BAND_REG_NOISEFILE'+eqsig,band_reg_noisefile,band_reg_noisefile_v,i,j)
        write2v(f,'BAND_NOISE_UNIFORMIZE_FSKY'+eqsig,band_noise_uniformize_fsky,band_noise_uniformize_fsky_v,i,j)
        write2v(f,'BAND_MASKFILE'+eqsig,band_maskfile,band_maskfile_v,i,j)
        write2v(f,'BAND_MASKFILE_CALIB'+eqsig,band_maskfile_calib,band_maskfile_calib_v,i,j)
        write2v(f,'BAND_BEAMTYPE'+eqsig,band_beamtype,band_beamtype_v,i,j)
        write2v(f,'BAND_BEAM_B_L_FILE'+eqsig,band_beam_b_l_file,band_beam_b_l_file_v,i,j)
        write2v(f,'BAND_BEAM_B_PTSRC_FILE'+eqsig,band_beam_b_ptsrc_file,band_beam_b_ptsrc_file_v,i,j)
        write2v(f,'BAND_PIXEL_WINDOW'+eqsig,band_pixel_window,band_pixel_window_v,i,j)
        write2v(f,'BAND_SAMP_NOISE_AMP'+eqsig,band_sample_noise_amp,band_sample_noise_amp_v,i,j)
        write2v(f,'BAND_BANDPASS_TYPE'+eqsig,band_bandpass_type,band_bandpass_type_v,i,j)
        write2v(f,'BAND_BANDPASS_MODEL'+eqsig,band_bandpass_model,band_bandpass_model_v,i,j)
        write2v(f,'BAND_NOMINAL_FREQ'+eqsig,nomfreqs,nomfreqs_v,i,j)
        write2v(f,'BAND_SAMP_BANDPASS'+eqsig,band_samp_bandpass,band_samp_bandpass_v,i,j)
        write2v(f,'BAND_BANDPASSFILE'+eqsig,band_bandpassfile,band_bandpassfile_v,i,j)
        write2v(f,'BAND_SAMP_GAIN'+eqsig,band_samp_gain,band_samp_gain_v,i,j)
        write2v(f,'BAND_GAIN_PRIOR_MEAN'+eqsig,band_gain_prior_mean,band_gain_prior_mean_v,i,j)
        write2v(f,'BAND_GAIN_PRIOR_RMS'+eqsig,band_gain_prior_mean,band_gain_prior_mean_v,i,j)
        write2v(f,'BAND_GAIN_CALIB_COMP'+eqsig,band_gain_calib_comp,band_gain_calib_comp_v,i,j)
        write2v(f,'BAND_GAIN_LMIN'+eqsig,band_gain_lmin,band_gain_lmin_v,i,j)
        write2v(f,'BAND_GAIN_LMAX'+eqsig,band_gain_lmax,band_gain_lmax_v,i,j)
        write2v(f,'BAND_GAIN_APOD_MASK'+eqsig,band_gain_apod_mask,band_gain_apod_mask_v,i,j)
        write2v(f,'BAND_GAIN_APOD_FWHM'+eqsig,band_gain_apod_fwhm,band_gain_apod_fwhm_v,i,j)
        write2v(f,'BAND_DEFAULT_GAIN'+eqsig,band_default_gain,band_default_gain_v,i,j)
        write2v(f,'BAND_DEFAULT_BP_DELTA'+eqsig,band_default_bp_delta,band_default_bp_delta_v,i,j)
        write2v(f,'BAND_DEFAULT_NOISEAMP'+eqsig,band_default_noiseamp,band_default_noiseamp_v,i,j)
        write2v(f,'BAND_COMPONENT_SENSITIVITY'+eqsig,band_component_sensitivity,band_component_sensitivity_v,i,j)
        write2v(f,'BAND_TOD_TYPE'+eqsig,band_tod_type,band_tod_type_v,i,j)
        write2v(f,'BAND_NOISE_RMS'+eqsig+'_SMOOTH01',band_noise_rms_sm1,band_noise_rms_sm1_v,i,j)
        write2v(f,'BAND_NOISE_RMS'+eqsig+'_SMOOTH02',band_noise_rms_sm2,band_noise_rms_sm2_v,i,j)
        write2v(f,'BAND_NOISE_RMS'+eqsig+'_SMOOTH03',band_noise_rms_sm3,band_noise_rms_sm3_v,i,j)
        f.close()

    