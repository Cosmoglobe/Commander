import numpy as np

#General purpose parameters
fac_channel_ids = ['RH', 'RL', 'LH', 'LL']
fac_scan_mode_ids = ['SS', 'SF', 'LS', 'LF', 'FL']
fac_scan_mode_idsl = ['SS', 'SF', 'LS', 'LF', 'FS', 'FL']
fac_att_soln = ['none', 'simulated', 'predicted', 'coarse', 'fine_without_dirbe',
                'fine_with_dirbe', 'definitive']

fac_jstart_default = '86001000000000'
fac_jstop_default = '99365235959990'

#Skymap parameters
fac_skymap_level = 5
fac_skymap_no_levels = 6
fac_time_offset = 14
fac_num_pixels = 6144

#FIRAS instrument parameters

#IFG peak position by MTM speed, MTM length, and channel
fac_ifg_peak = np.array([356, 354, 358, 354, 359, 89, 361, 89, 356, 354, 358, 354, 359, 89, 361, 89],
                        dtype=int)
fac_ifg_peak.shape = (4, 2, 2)
fac_ifg_peak = fac_ifg_peak.transpose()

#IFG peak position by scan mode (SS=0, SF=1, LS=2, LF=3) and channel
fac_ifg_peak_sm = np.array([356, 358, 354, 354, 359, 361, 89, 89, 356, 358, 354, 354, 359, 361, 89, 89],
                           dtype=int)
fac_ifg_peak_sm.shape = (4, 4)
fac_ifg_peak_sm = fac_ifg_peak_sm.transpose()

#MTM scan times by MTM speed, MTM length, and channel
fac_scan_times = np.array([[2.2, 8.8], [1.5, 6.0]]).transpose()
fac_scan_times_sm = np.array([2.2, 1.5, 8.8, 6.0])

#Etendu, 10 volt adc scale, and sampling rate
fac_etendu = 1.5
fac_adc_scale = 204.75
fac_epoch = 2000.0

#Load resistance, icm to GHz conversion factor, spectrum unit conversion factors, counts to volts
#conversion factor
fac_load_resist = 4.0e7
fac_icm_ghz = 29.9792458
fac_watt_to_mjy = 1.0e15 / fac_icm_ghz
fac_erg_to_mjy = 1.0e8 / fac_icm_ghz
fac_erg_to_watt = 1.0e-7
fac_count_to_volt = 25.5

#Fast Fourier ad spectrum lengths; Nyquist frequency correction
fac_fft_length = [640, 640, 640, 640, 160, 160]
fac_spec_length = [321, 321, 321, 321, 81, 81]
fac_nyq_correct = 1.00159

#Data Quality summary flags
fac_no_dq = 0
fac_good_dq = 1
fac_some_yellow = 2
fac_many_yellow = 3
fac_no_y_some_r = 4
fac_some_y_some_r = 5
fac_many_y_some_r = 6
fac_many_red = 7
fac_no_att_soln = 32
fac_no_hskp = 64
fac_tlm_error = 127

#Attitude parameters
fac_att_conv = 1.0e-4 * 180.0 / np.pi
fac_att_conv_rad = 1.0e-4

#XCal position parameters
fac_xcaltrans = 3
fac_xcalout = 2
fac_xcalin = 1
fac_xcalposerr = 0

#Calibration Parameters
fac_no_find_peak = -2
fac_both = 3
fac_slow = 1
fac_fast = 2
fac_long = 2
fac_short = 1
fac_max_struct = 6
fac_max_periods = 14
fac_left = 14
fac_right = 12
fac_num_model = 6
fac_num_emiss = 5
fac_phase_temp = 2

#Record sizes and the like
fac_max_num = 100 #maximum number of allowable records to coadd
#fac_max_num = 1000

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

fac_sky = 1
fac_cal = 2
fac_raw = 3

fac_min_gain_set = 0
fac_max_gain_set = 7

fac_min_adds = 1
fac_max_adds = 12

fac_max_deglitch = 1000

fac_apco_date = 9626073 #aperture cover ejection date in VAX ADT time
fac_apco_gmt = '893251118' #aperture cover ejection date in GMT
fac_vax_year_len = 73426.0 #length of a year in VAX ADT seconds

fac_gains = np.array([1, 3, 10, 30, 100, 300, 1000, 3000])

fac_channel_ids = ['rh', 'rl', 'lh', 'll']
fac_scan_mode_ids = ['ss', 'sf', 'ls', 'lf', 'fl']
fac_scan_mode_idsl = ['ss', 'sf', 'ls', 'lf', 'fs', 'fl']

