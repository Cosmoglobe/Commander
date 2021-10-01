import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import h5py
from astropy.io import fits

dirbe_dir = '/mn/stornext/d16/cmbco/ola/dirbe/'

def dirbe_pix_to_lat_lon(cio_pixel_number):
    # Load data from 'DIRBE_SKYMAP_INFO.FITS
    # take in pixel numbers from a file in 'cio/' and return 
    # ecliptic longitude and latitude
    skymap_info = fits.open("/mn/stornext/d16/cmbco/ola/dirbe/DIRBE_SKYMAP_INFO.FITS")
    qspixels = skymap_info[1].data["QSPIXEL"]
    index = np.where(qspixels == cio_pixel_number)
    elat = skymap_info[1].data["ELAT-CSC"]
    elon = skymap_info[1].data["ELON-CSC"]
    return np.array((elat[index][0], elon[index][0]))

"""
ex_file = fits.open("/mn/stornext/d16/cmbco/ola/dirbe/DIRBE_SKYMAP_INFO.FITS")
ex_qspixels = ex_file[1].data["QSPIXEL"]
coords = np.zeros((len(ex_qspixels),2))

print(dirbe_pix_to_lat_lon(0).shape)

for i in range(1000):
    coords[i] = dirbe_pix_to_lat_lon(ex_qspixels[i])

print("passed for")
plt.figure()
plt.plot(coords[:,0], coords[:,1])
plt.savefig("test.png")
plt.show()
plt.close()
print("tried plotting")
"""

# Function reformatting a single-day fits file to ten hdf5 files, one for each frequency band
def FITS_to_hdf5(fits_filename, n_day):
    # Open fits file
    fits_file = fits.open(fits_filename)
    
    # f5 file names
    f5_filenames = np.array(('Phot01.hdf5', 'Phot02.hdf5', 'Phot03.hdf5', 'Phot04.hdf5', 'Phot05.hdf5', 'Phot06.hdf5', 'Phot07.hdf5', 'Phot08.hdf5', 'Phot09.hdf5', 'Phot10.hdf5'))
    

    # Loop over the file for each band
    for i in range(len(f5_filenames)):
        h5_file = h5py.File('dirbe_hdf5_files/' + f5_filenames[i], 'r+')
        # Create time group
        day_grp = h5_file.create_group(n_day)
        # Create subgroup for each detector
        if (i in np.array((0, 1, 2))):
            # Bands Phot01, Phot02 and Phot03 each have three detectors A, B and C
            detectorA_subgrp = day_grp.create_group('A')
            detectorB_subgrp = day_grp.create_group('B')
            detectorC_subgrp = day_grp.create_group('C')
            
            # Storing the CIO data
            cio_A = detectorA_subgrp.create_dataset('default', data=fits_file[1].data.field(i+6))
            cio_B = detectorB_subgrp.create_dataset('default', data=fits_file[1].data.field(i+7))
            cio_C = detectorC_subgrp.create_dataset('default', data=fits_file[1].data.field(i+8))

        else:
            # Bands Phot04 through Phot10 each have one detector
            detectorA_subgrp = day_grp.create_group('A')
            cio = detectorA_subgrp.create_dataset('default', data=fits_file[1].data.field(i))

# Function parsing through the CIO directory reformatting the data from all the  single-day
# fits files into ten hdf5 files, one for each frequency band
def parse_directory():
    # The measurement days on format yyddd for years 1989 and 1990
    endings89 = np.arange(89345, 89366)
    endings90 = np.arange(90001, 90265)
    endings = np.concatenate((endings89, endings90))
    
    # Loop through fits files for each day
    for i in range(len(endings)):
        FITS_to_hdf5('/mn/stornext/d16/cmbco/ola/dirbe/cio/DIRBE_CIO_P3B_' + str(endings[i]) + '.FITS', '%06d' % (i+1))

#parse_directory()
