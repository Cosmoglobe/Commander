import h5py
import numpy as np
import litebird_sim as lbs
import sys
from pathlib import Path
# TODO: In the future version you should use `cosmoglobe` and so eliminating
# the need to do this workaround in the first place 
# 
# Getting full path to Mathew's library as an object
commander_tools_path = Path(__file__).absolute().parents[2].joinpath(
        'python','commander_tools').resolve()
# Appending the path to `PYTHONPATH`, so no need to 
# modify it externally (in your `.bashrc` etc.)
sys.path.append(str(commander_tools_path))
# Importing necessary modules from Mathew's library 
#from tod_tools import commander_instrument as comm_inst

from tod_tools.litebird_imo import LitebirdImo


def write_chunk_to_file(f, chunk, num_specifier=None):
    lines = []
    for param, param_val in chunk.items():
        if num_specifier is None:
            lines.append(f'{param} = {param_val}\n')
        else:
            if '@' in param:
                t_param = param
                t_param.replace("@", f'{num_specifier:03}')
                lines.append(f'{t_param} = {param_val}\n')
            else:
                lines.append(f'{param}{num_specifier:03} = {param_val}\n')
    lines.append('\n')
    f.writelines(lines)


def copy_component_chunk(f, infile):
    lines = infile.readlines()
    lines_to_copy = []
    num_comps = 0
    for line in lines:
        if line.startswith('COMP'):
            lines_to_copy.append(line)
        if line.startswith('COMP_TYPE'):
            num_comps += 1
    lines_to_copy.insert(0, '\n')
    for i in range(num_comps, 0, -1):
        lines_to_copy.insert(0, f'INCLUDE_COMP{i:02} = .true.\n')
    lines_to_copy.insert(0, f'NUM_SIGNAL_COMPONENTS = {num_comps}\n')
    lines_to_copy.insert(0, 'CMB_DIPOLE_PRIOR = none\n')
    lines_to_copy.insert(0, 'INIT_INSTRUMENT_FROM_HDF = default\n')
    lines_to_copy.insert(0, 'INSTRUMENT_PARAM_FILE = instrument_params_LB.dat\n')
    f.writelines(lines_to_copy)


def create_band_chunk(f, instrument, band_name, imo, band_num):
    beamfile = f'beam_LB_imo{imo.version}_{band_name}.fits'
    center_freq = imo.get_detector_frequency(band_name)
    numscans = 143
    detlist = imo.get_channel_dets(band_name)
    instrument_files = {
        'LFT': 'LFT_instrument.h5',
        'MFT': 'MFT_instrument.h5',
        'HFT': 'HFT_instrument.h5'
    }

    params = {
        'BAND_LABEL': band_name,
        'BAND_OBS_PERIOD': '1',
        'BAND_POLARIZATION': '.true.',
        'BAND_NSIDE' : '512',
        'BAND_LMAX': '1000',
        'BAND_UNIT': 'uK_cmb',
        'BAND_NOISE_FORMAT': 'rms',
        'BAND_MAPFILE': 'sim_LBUF_imov1_cmb_d1s1_noise_0000_LFT_60_n512_51arcmin.fits',
        'BAND_NOISEFILE': 'LBUF_imov1_rmsnoiseLFT_60_n512.fits',
        'BAND_REG_NOISEFILE': 'none',
        'BAND_NOISE_RMS@_SMOOTH01': 'none',
        'BAND_NOISE_UNIFORMIZE_FSKY': '0.00',
        'BAND_MASKFILE': 'fullsky',
        'BAND_MASKFILE_CALIB': 'fullsky',
        'BAND_BEAMTYPE': 'b_l',
        'BAND_BEAM_B_L_FILE': beamfile,
        'BAND_BEAM_B_PTSRC_FILE': 'none',
        'BAND_PIXEL_WINDOW': 'pixel_window_n0512.fits',
        'BAND_SAMP_NOISE_AMP': '.false.',
        'BAND_BANDPASS_TYPE': 'delta',
        'BAND_NOMINAL_FREQ': f'{center_freq:.1f}',
        'BAND_SAMP_BANDPASS': '.false.',
        'BAND_BANDPASSFILE': instrument_files[instrument],
        'BAND_SAMP_GAIN': '.false.',
        'BAND_GAIN_CALIB_COMP': 'cmb',
        'BAND_GAIN_LMIN': '25',
        'BAND_GAIN_LMAX': '100',
        'BAND_GAIN_APOD_MASK': 'fullsky',
        'BAND_GAIN_APOD_FWHM': '120.',
        'BAND_DEFAULT_GAIN': '1.',
        'BAND_DEFAULT_BP_DELTA': '0.',
        'BAND_DEFAULT_NOISEAMP': '1.',
        'BAND_COMPONENT_SENSITIVIY': 'broadband',
        'BAND_TOD_MAIN_PROCMASK': 'mask_fullsky_n0512.fits',
        'BAND_TOD_SMALL_PROCMASK': 'mask_fullsky_n0512.fits',
        'BAND_TOD_BP_INIT_PROP': f'bp_init_{band_name}.dat',
        'BAND_TOD_RIMO': instrument_files[instrument],
        'BAND_TOD_FILELIST': f'filelist_LB_sims_coadded_{band_name}.txt',
        'BAND_TOD_START_SCANID': '1',
        'BAND_TOD_END_SCANID': f'{numscans}',
        'BAND_TOD_TOT_NUMSCAN': f'{numscans}',
        'BAND_TOD_FLAG': '0',
        'BAND_TOD_ORBITAL_ONLY_ABSCAL': '.false.',
        'BAND_TOD_DETECTOR_LIST': ','.join(detlist),
        'BAND_TOD_INIT_FROM_HDF': 'default'
    }

    write_chunk_to_file(f, params, band_num)

def create_algorithm_param_chunk(f):
    output_dir = '/mn/stornext/u3/eirikgje/data/litebird_sim/chains/'
    params = {
        'NUMCHAIN': '1',
        'NUM_GIBBS_ITER': '100',
        'NUM_ITER_WITH_ML_SEARCH': '0',
        'BASE_SEED': '2938109',
        'CHAIN_STATUS': 'append',
        'INIT_CHAIN': 'none',
        'INIT_SAMPLE_NUMBER': '1',
        'NUM_GIBBS_STEPS_PER_TOD_SAMPLE': '1',
        'SAMPLE_ONLY_POLARIZATION': '.true.',
        'SAMPLE_SIGNAL_AMPLITUDES': '.true.',
        'SAMPLE_SPECTRAL_INDICES': '.true.',
        'ENABLE_TOD_ANALYSIS': '.false.',
        'TOD_OUTPUT_4D_MAP_EVERY_NTH_ITER': '0',
        'TOD_INCLUDE_ZODI': '.false.',
        'FFTW3_MAGIC_NUMBERS': 'data_LB/fft3_magic_numbers_230810.txt',
        'NSKIP_FILELIST': 0,
        'RESAMPLE_CMB': '.false.',
        'FIRST_SAMPLE_FOR_CMB_RESAMP': '1',
        'LAST_SAMPLE_FOR_CMB_RESAMP': '15',
        'NUM_SUBSAMP_PER_MAIN_SAMPLE': '10',
        'CG_CONVERGENCE_CRITERION' : 'residual',
        'CG_LMAX_PRECOND': '-1',
        'CG_MAXITER': '500',
        'CG_MINITER': '5',
        'CG_TOLERANCE': '1.d-8',
        'CG_CONV_CHECK_FREQUENCY': '1',
        'CG_PRECOND_TYPE': 'pseudoinv',
        'CG_INIT_AMPS_ON_ZERO': '.false.',
        'SET_ALL_NOISE_MAPS_TO_MEAN': '.false.',
        'NUM_INDEX_CYCLES_PER_ITERATION': '1',
        'IGNORE_GAIN_AND_BANDPASS_CORR': '.false',
        'OUTPUT_DIRECTORY': f'{output_dir}',
        'THINNING_FACTOR': '20',
        'NSIDE_CHISQ': 16,
        'POLARIZATION_CHISQ': '.true.',
        'OUTPUT_MIXING_MATRIX': '.false.',
        'OUTPUT_RESIDUAL_MAPS': '.true.',
        'OUTPUT_CHISQ_MAP': '.true.',
        'OUTPUT_EVERY_NTH_CG_ITERATION': '20',
        'OUTPUT_CG_PRECOND_EIGENVALS': '.false.',
        'OUTPUT_INPUT_MODEL': '.false.',
        'OUTPUT_DEBUG_SEDS': '.false.',
        'OUTPUT_SIGNALS_PER_BAND': '.false.',
    }
    write_chunk_to_file(f, params)

def create_dataset_chunk(f, total_num_channels):
    data_dir = '/mn/stornext/u3/hke/xsan/commander3/BP9/data'
    params = {
        'DATA_DIRECTORY': f'{data_dir}',
        'NUMBAND': f'{total_num_channels}',
        'SOURCE_MASKFILE': 'none',
        'PROCESSING_MASKFILE': 'none',
        'PROCESSING_MASKFILE2': 'none',
        'PROC_SMOOTH_SCALE': '30.',
        'NUM_SMOOTHING_SCALES': '0',
        'TOD_NUM_BP_PROPOSALS_PER_ITER': '0'
    }
    for i in range(1, total_num_channels+1):
        params[f'INCLUDE_BAND{i:03}'] = '.true.'
    for i in range(1, total_num_channels+1):
        params[f'BAND_TOD_TYPE{i:03}'] = 'LB'
    write_chunk_to_file(f, params)


imo_version = 'v1.3'
imo_db_interface = lbs.Imo()
with open('paramfile_LB_commander.txt', 'w') as f:
    create_algorithm_param_chunk(f)
    total_num_channels = 0
    for instrument in ('LFT', 'MFT', 'HFT'):
        imo = LitebirdImo(imo_db_interface, imo_version, instrument) 
        channels = imo.get_channel_names()
        total_num_channels += len(channels)
    create_dataset_chunk(f, total_num_channels)
    k = 0
    for instrument in ('LFT', 'MFT', 'HFT'):
        imo = LitebirdImo(imo_db_interface, imo_version, instrument)
        channels = imo.get_channel_names()
        for i, channel in enumerate(channels):
            k += 1
            create_band_chunk(f, instrument, channel, imo, k)
    with open('param_LB_com2_d1s1_Unni_newbeams.txt', 'w') as infile:
        copy_component_chunk(f, infile)


#for instrument in ('LFT', 'MFT', 'HFT'):
#    imo = LitebirdImo(imo_db_interface, imo_version, instrument) 
#    f = h5py.File(fnames[instrument], 'w')
#    channels = imo.get_channel_names()
##    center_freqs = []
#    for channel in channels:
#        center_freq = imo.get_detector_frequency(channel)
#        bandwidth = imo.get_detector_bandwidth(channel)
#        fwhm = imo.get_detector_fwhm(channel)
##        center_freqs.append(imo.get_detector_frequency(channel))
#        grp = f.create_group(channel)
#        yval = grp.create_dataset("bandpass", (3,), dtype='f')
#        xval = grp.create_dataset("bandpassx", (3,), dtype='f')
#        #Tophat
#        yval[:] = 1.
#        xval[:] = np.array([center_freq - bandwidth/2,
#                            center_freq,
#                            center_freq + bandwidth / 2])
#        for detector in imo.get_channel_dets(channel):
#            grp = f.create_group(detector)
#            yval = grp.create_dataset("bandpass", (3,), dtype='f')
#            xval = grp.create_dataset("bandpassx", (3,), dtype='f')
#            #Tophat
#            yval[:] = 1.
#            xval[:] = np.array([center_freq - bandwidth/2,
#                                center_freq,
#                                center_freq + bandwidth / 2])
#            b = grp.create_group("beam")
#            bT = b.create_dataset("T", (1,), dtype="f")
#            bT[:] = 0.
#            bB = b.create_dataset("B", (1,), dtype="f")
#            bB[:] = 0
#            bE = b.create_dataset("E", (1,), dtype="f")
#            bE[:] = 0
#            blmax = grp.create_dataset("beamlmax", (1,), dtype='f')
#            blmax[0] = 2400
#            bmmax = grp.create_dataset('beammax', (1,), dtype='f')
#            bmmax[0] = 100
#            el = grp.create_dataset("elip", (1,), dtype='f')
#            el[0] = 1.
#            fw = grp.create_dataset("fwhm", (1,), dtype='f')
#            fw[0] = fwhm
#            psi = grp.create_dataset("psi_ell", (1,), dtype='f')
#            psi[0] = 1.
#            sl = grp.create_group("sl")
#            slT = sl.create_dataset("T", (1,), dtype="f")
#            slT[:] = 0
#            slB = sl.create_dataset("B", (1,), dtype="f")
#            slB[:] = 0
#            slE = sl.create_dataset("E", (1,), dtype="f")
#            slE[:] = 0
#            sll= grp.create_dataset("sllmax", (1,), dtype='f')
#            sll[0]=512
#            slm= grp.create_dataset("slmmax", (1,), dtype='f')
#            slm[0]=100
#            eff= grp.create_dataset("mbeam_eff", (1,), dtype='f')
#            eff[0]=1.
#
#    f.close()
