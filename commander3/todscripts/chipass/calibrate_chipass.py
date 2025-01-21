import numpy as np 
import matplotlib.pyplot as plt
from astropy.io import fits
# import tools 
import h5py
import multiprocessing
import glob
import os.path
from astropy.time import Time 
import astropy.coordinates as coord 
import astropy.units as u 
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation  
from astropy.coordinates import get_body_barycentric, get_body, get_moon
from astropy.coordinates import AltAz

def median_filter(tod, n_med=25):
    n = tod.shape[0]
    out = np.zeros_like(tod)
    n_half = n_med // 2
    for i in range(n):
        low = max(i-n_half, 0)
        hi = min(i+n_half, n)
        out[i] = np.nanmedian(tod[low:hi], axis=0)
    return out

def calibrate_scan(filename):
    data_raw = fits.open(filename)
    # print(data_raw[1].columns)
    t = data_raw[1].data['time'].reshape(-1,13)[:, 0]  # time, in seconds, since start of day (UTC)
    nt = len(t)  # number of samples in time 
    date = np.array(data_raw[1].data['date-obs']).reshape(-1,13)[:, 0]
    scan = np.array(data_raw[1].data['object']).reshape(-1,13)[0, 0]
    
    t0 = Time(date, format='isot', scale='utc').mjd
    mjd = t0 + t / 3600.0 / 24

    scan_name  = scan.strip(' ') + date[0] + '_%9.3f' % t[0]
    print('Calibrating scan:', scan_name)

    # Check if scan already in file
    if os.path.isfile('chipass_full.h5'):
        f1 = h5py.File('chipass_full.h5', 'r')
        skip = scan_name + '/mjd' in f1
        f1.close()
    else:
        skip = False

    if skip:
        outdata = {
            'skip': skip,
            'scan': scan_name,
        }
        return outdata

    # Pointing info 
    # (at some point we should also add the polarization angle, 
    # although that will requre some care to get right)
    az = data_raw[1].data['AZIMUTH'].reshape(-1,13)
    el = data_raw[1].data['ELEVATIO'].reshape(-1,13)

    # Lat/Lon/h (WGS84) of Parkes telescope 32.9984064, 148.2635101, 414.80 (https://www.narrabri.atnf.csiro.au/observing/users_guide/html/chunked/apg.html)
    loc = coord.EarthLocation(lon=148.2635101 * u.deg, lat=-32.9984064 * u.deg, height=414.80 * u.m)
    time = Time(mjd, format='mjd')
    AltAz = SkyCoord(alt=el*u.deg, az=az*u.deg, obstime=time[:, None], frame='altaz', location=loc)
    point = np.array([AltAz.galactic.l.degree, AltAz.galactic.b.degree, np.zeros_like(AltAz.galactic.b.degree)]).transpose((1, 2, 0))
    tod = data_raw[1].data['Data'][:,0,0,:,:]
    tod = np.reshape(tod, (-1, 13, 2, 1024))
    tsys = data_raw[1].data['Tsys']
    tsys = np.reshape(tsys, (-1, 13, 2))
    flag = data_raw[1].data['flagged'][:,0,0,:,:]
    flag = np.reshape(flag, (-1, 13, 2, 1024))
    mask = np.array(1.0 - flag).astype(float)
    ratio = tod[:, :, :, :] / tsys[:, :, :, None]
    bst = median_filter(ratio)
    bt = median_filter(tsys)
    s_prime = tod / bst - bt[:, :, :, None]  # Calibrated data in Jy
    
    # bandpass integration (using median for now)
    s_cont = np.zeros((nt, 13, 2, 4))
    mask[:, :, :, 85:125] = 0  # remove 
    mask[(mask == 0)] = np.nan
    for i in range(4):
        s_cont[:, :, :, i] = np.nanmedian(mask[:, :, :, i*256:(i+1)*256] * s_prime[:, :, :, i*256:(i+1)*256], axis=3)

    flag = np.isnan(s_cont)
    point = point.astype(np.float32)
    tod = tod.astype(np.float32)
    outdata = {
        'skip': skip,
        'scan': scan_name,
        'mjd': mjd,
        'point_gal':point, 
        'tod':s_cont,
        'flag':flag,
    }
    return outdata
    
def save_to_hdf(outdata):
    if outdata['skip']:
        print('Skipping already calculated scan: ', outdata['scan'])
        return 

    f1 = h5py.File('chipass_test.h5', 'a')
    scan = outdata['scan']
    del outdata['scan']
    del outdata['skip']
    prefactor = str(scan.strip(' '))
    
    skip = prefactor + '/mjd' in f1
    if skip:
        f1.close()
        print('Skipping already calculated scan: ', prefactor)
        return
    
    for key, value in outdata.items():
        f1.create_dataset(prefactor + '/' + key, data=value)
    f1.close()

if __name__ == '__main__':
    foldername = '.'  # '/mn/stornext/d16/cmbco/ola/chipass/l1_data'  
    file_list = glob.glob(foldername + "/**/*.sdfits", recursive=True)

    nproc = 2
    pool = multiprocessing.Pool(nproc)
    for filename in file_list:
        pool.apply_async(calibrate_scan, (filename, ), callback=save_to_hdf)
    pool.close()
    pool.join()
    
    # pool.map(calibrate_scan, file_list)
