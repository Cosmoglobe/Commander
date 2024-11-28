import os

import numpy as np

import firas_pipeline.misc.util as util

def read_config(path, grtrawwt=False, grtcoawt=False, grttrans=False, mincoadd=False, samprate=False,
                cmdgain=False, nyquistl=False, cth=False, gltchcor=False, vibcorr=False,
                samprate_type='mission'):
 
    config = {}
    if grtrawwt:
        #fn = os.path.join(path, "FIRAS_FEX_GRTRAWWT.TXT")
        fn = os.path.join(path, "fex_grtrawwt.txt")
        config['grtrawwt'] = read_grtwt(fn)
        
    if grtcoawt:
        #fn = os.path.join(path, "FIRAS_FEX_GRTCOAWT.TXT")
        fn = os.path.join(path, "fex_grtcoawt.txt")
        config['grtcoawt'] = read_grtwt(fn)

    if grttrans:
        #fn = os.path.join(path, "FIRAS_FEX_GRTTRANS.TXT")
        fn = os.path.join(path, "fex_grttrans.txt")
        config['grttrans'] = read_grttrans(fn)
    
    if mincoadd:
        #fn = os.path.join(path, "FIRAS_FEX_MINCOADD.TXT")
        fn = os.path.join(path, "fex_mincoadd.txt")
        config['mincoadd'] = read_mincoadd(fn)

    if cmdgain:
        fn = os.path.join(path, "fex_cmdgain.txt")
        config['cmdgain'] = read_cmdgain(fn)

    if samprate:
        fn = os.path.join(path, "fex_samprate.txt")
        config['samprate'] = read_samprate(fn, samprate_type)

    if nyquistl:
        fn_samp = os.path.join(path, "fex_samprate.txt")
        fn_nyquist = os.path.join(path, "fex_nyquist.txt")
        config['nyquistl'] = gen_nyquistl(fn_samp, fn_nyquist, samprate_type)

    if cth:
        fn_cth = os.path.join(path, "fex_cth.txt")
        config['cth'] = read_cth(fn_cth)

    if gltchcor:
        fn_gltchcor = os.path.join(path, "fex_gltchcor.txt")
        config['gltchcor'] = read_gltchcor(fn_gltchcor)

    if vibcorr:
        fn_vibcorr = os.path.join(path, "fex_vibcorrl.txt")
        config['vibcorr'] = read_vibcorr(fn_vibcorr)
 
    return config

def read_grtwt(fn):
    data = np.loadtxt(fn, comments='!')
    
    output = {}
    output['grt_a_weight'] = data[:16]
    output['grt_b_weight'] = data[16:]

    return output

def read_mincoadd(fn):

    data = np.loadtxt(fn, comments='!')

    return {'min_ifg_coadd' : data.astype(int)}

def read_grttrans(fn):
    data = np.loadtxt(fn, comments='!', delimiter=',')
    
    output = {}
    output['grt_a_trans_temp'] = data[:16, 0]
    output['grt_b_trans_temp'] = data[16:, 0]
    output['grt_a_trans_hwid'] = data[:16, 1]
    output['grt_b_trans_hwid'] = data[16:, 1]

    return output

def read_cmdgain(fn):
    data = np.loadtxt(fn, comments='C')

    output = np.transpose(data)

    return output

def read_samprate(fn, samprate_type):
    data = np.loadtxt(fn, comments='!')

    if samprate_type == 'int':
        return data[0]
    else:
        return data[1]

def gen_nyquistl(fn_samp, fn_nyquist, samprate_type):
    '''This function reads the MTM sampling rates from the samprate text file and the optical
    Nyquist frequency correction from the nyquist fuke and then compute the Nyquist frequencies
    in icm and Hz for all channels and scan modes for either I&T data or flight data
    '''

    mtm_speed = [0, 1, 0, 1, 0, 1, 0, 1, 1, 1]
    fakeit = [1.0, 2.0, 3.0, 8.0, 12.0]
    multiplier = [6.0, 4.0]
    ngroup = [3.0, 2.0, 3.0, 2.0, 3.0, 2.0, 12.0, 8.0, 8.0, 8.0]
    
    fringes = 20.00e-4
    optical_path = 4*np.cos(np.pi/6)

    fac_fft_length = [640, 640, 640, 640, 160, 160]
    
    data_samp = np.loadtxt(fn_samp, comments='!')

    if samprate_type == 'int':
        sampling_rate = data_samp[0]
    else:
        sampling_rate = data_samp[1]
    
    freq_shift = np.loadtxt(fn_nyquist, comments='!')
    
    icm = np.zeros(10)
    hz = np.zeros(15)

    #Calculate the Nyquist frequencies
    
    #Calculate the optical Nyquist frequencies
    norm = freq_shift / fringes / optical_path / 2.0
    for j in range(10):
        icm[j] = norm * multiplier[mtm_speed[j]] / ngroup[j]

    #Calculate the scan speeds
    speed = sampling_rate/icm[0]/multiplier

    #Calculate the electronic Nyquist frequencies
    for j in range(10):
        hz[j] = icm[j] * speed[mtm_speed[j]]

    #Calculate the fakeit Nyquist frequencies
    for j in range(10, 15):
        hz[j] = fac_fft_length[0] / fakeit[j-10]

    output = {}
    output['icm'] = icm
    output['hz'] = hz

    return output

def read_cth(fn):

    def get_nextline(f, type=float, nlines=1):

        array = []
        for i in range(nlines):
            line = f.readline()
            line = line.partition('!')[0].rstrip().split(' ')

            if type is not None:
                data = np.array(list(map(type, line)))
            else:
                data = line

            array.append(data)

        data_output = np.array(array)

        if nlines == 1:
            data_output.shape = (data_output.size,)

        return data_output
    
    output = {}

    f = open(fn, 'r')

    #Values for instrument state consistency check
    output['bolometer_voltage_tolerances'] = get_nextline(f)
    output['grt_tolerances'] = get_nextline(f)
    output['spatial_gradients'] =  get_nextline(f, nlines=5)
    data = get_nextline(f, nlines=2, type=int)
    output['temporal_gradients'] = np.reshape(data, (20, ))

    #Values for neighbor selection
    output['gal_lat_cutoff'] = get_nextline(f)

    #Values for secondary template formation
    output['upper_bounds'] = get_nextline(f)
    output['lower_bounds'] = get_nextline(f)
    output['peak_points'] = get_nextline(f, type=int)
    output['prim_temp_amp'] = get_nextline(f, nlines=4)
    output['prim_temp_snr'] = get_nextline(f, nlines=4)
    output['sec_temp_amp'] = get_nextline(f, nlines=4)
    output['sec_temp_snr'] = get_nextline(f, nlines=4)

    #Values for shape consistency check
    output['min_ifg_noise'] = get_nextline(f)
    output['max_ifg_noise'] = get_nextline(f)
    output['mask_ifg'] = get_nextline(f, nlines=2, type=int)
    output['max_point_deviation'] = get_nextline(f)
    output['max_bad_points'] = get_nextline(f, type=int)
    output['min_ifg_shape_coadd'] = get_nextline(f, type=int)

    #Values for binning cutoff for statistics
    output['xcal_temp'] = get_nextline(f)[0]
    output['ical_temp'] = get_nextline(f)[0]
    output['skyhorn_refhorn_avg'] = get_nextline(f)[0]
    output['dihedral_temp'] = get_nextline(f)[0]
    output['bolometer_temp'] = get_nextline(f)
    output['glitch_rate'] = get_nextline(f)
    output['abs_galactic_lat'] = get_nextline(f, type=int)[0]
    output['gmt_time_saa'] = get_nextline(f, type=None)[0]
    output['time_saa'] = util.gmt_to_binary(output['gmt_time_saa'])

    return output

def read_gltchcor(fn):
    
    data = np.loadtxt(fn, comments='!')

    output = {}
    slope = np.zeros(12)
    intercept = np.zeros(12)

    slope[0:3] = data[0, :]
    slope[3:6] = data[1, :]
    slope[6:9] = data[2, :]
    slope[9:12] = data[3, :]
    
    intercept[0:3] = data[4, :]
    intercept[3:6] = data[5, :]
    intercept[6:9] = data[6, :]
    intercept[9:12] = data[7, :]

    output['slope'] = slope
    output['intercept'] = intercept

    return output

def read_vibcorr(fn):

    #Primary offsets
    #Secondary offsets

    #RHSS RHSF RHLF RLSS RLFS/FL RLLF
    #LHSS LHSF LHLF LLSS LLFS/FL LLLF
    
    data = np.loadtxt(fn, comments='!')

    output = {}
    output['primary_offset'] = np.zeros(12)
    output['secondary_offset'] = np.zeros(12)

    output['primary_offset'][:6] = data[0, :]
    output['primary_offset'][6:] = data[1, :]
    output['secondary_offset'][:6] = data[2, :]
    output['secondary_offset'][6:] = data[3, :]
    
    return output

if __name__ == '__main__':
    print("NOT ACTUALLY TESTING ANYTHING")

