import numpy as np
import astropy.time as at
import h5py
from astropy.io import fits
from tqdm import tqdm
from datetime import datetime

# Function reformatting a single-day fits file to ten hdf5 files,
# one for each frequency band
def FITS_to_hdf5(fits_filename, n_day):
    # Open fits file
    fits_data = fits.open(fits_filename)[1].data

    # Time ordered indices
    time_order = np.argsort(fits_data['time'])

    # f5 file names
    f5_filenames = np.array(('Phot01_' + str(nside) + '.hdf5', \
    'Phot02_' + str(nside) + '.hdf5', 'Phot03_' + str(nside) + \
    '.hdf5', 'Phot04_' + str(nside) + '.hdf5', 'Phot05_' + \
    str(nside) + '.hdf5', 'Phot06_' + str(nside) + '.hdf5', \
    'Phot07_' + str(nside) + '.hdf5', 'Phot08_' + str(nside) + \
    '.hdf5', 'Phot09_' + str(nside) + '.hdf5', 'Phot10_' + \
    str(nside) + '.hdf5'))

    # Band names
    bands = np.array(('Phot1', 'Phot2', 'Phot3', 'Phot04', \
    'Phot05', 'Phot06', 'Phot07', 'Phot08', 'Phot09', 'Phot10'))

    # Time array
    time = (at.Time(datetime(day=1,month=1,year=1981)) + \
    at.TimeDelta(fits_data['Time'][time_order], format='sec')).mjd

    # Flags
    rad_zone = fits_data['RadZone'][time_order]   # Radiation zone flags

    oa_flags = fits_data['OA_flags'][time_order]  # Orbit and attitude
                                                  # flags
    xs_noise = fits_data['XSNoise'][time_order]   # Excess noise flags
                                                  # (per detector)

    # Moon and Jupiter flags
    M2LOS = fits_data['Moon2LOS'][time_order]
    J2LOS = fits_data['Jup2LOS'][time_order]

    # HEALpix pixel numbers
    pixel_number           = fits_data['Pixel_no'][time_order]
    pixel_sub_position     = fits_data['PSubPos'][time_order]
    pixel_sub_sub_position = fits_data['PSbSbPos'][time_order]
    hp_pix = csc2healpix(pixel_number, pixel_sub_position, \
             pixel_sub_sub_position, nside)

    # Flagging Moon2LOS and Jupiter2LOS if relevant
    M2LOS_flags = np.zeros(len(M2LOS), dtype=int)
    if np.any(M2LOS <= 10):
        M2LOS_flags[M2LOS <= 10] = 1

    J2LOS_flags = np.zeros(len(J2LOS), dtype=int)
    if np.any(J2LOS <= 1.5):
        J2LOS_flags[J2LOS <= 1.5] = 1

    general_flags = rad_zone + 2**3*oa_flags + 2**11*M2LOS_flags +\
                    2**12*J2LOS_flags

    # Polarisation angle
    pol_theta = polarisation_theta(fits_data)

    # Loop over the file for each band
    for i in range(len(f5_filenames)):
        # Opening the hdf5 file
        h5_file = h5py.File('dirbe_hdf5_files/' + f5_filenames[i], 'r+')
        # Create time group
        day_grp = h5_file.create_group(n_day)

        # Create subgroup for each detector
        if (i in np.array((0, 1, 2))):
            # Bands Phot01, Phot02 and Phot03 each have three \
            # detectors A, B and C
            detectorA_subgrp = day_grp.create_group('A')
            detectorB_subgrp = day_grp.create_group('B')
            detectorC_subgrp = day_grp.create_group('C')

            # Storing the tod, flags, HEALpix pixels, time and
            # polarisation angles for each detector
            tod_A          = detectorA_subgrp.create_dataset('tod', \
                             data=fits_data[bands[i] + 'A'][time_order])
            flags_A        = detectorA_subgrp.create_dataset('flag', \
                    data = general_flags + 2**10*(np.bitwise_and(\
                    xs_noise, 2**i)/2**i).astype(int))
            healpix_A      = detectorA_subgrp.create_dataset('pix', \
                             data=hp_pix)
            time_A         = detectorA_subgrp.create_dataset('time', \
                             data=time)
            polarisation_A = detectorA_subgrp.create_dataset('polang', \
                             data=pol_theta)

            tod_B          = detectorB_subgrp.create_dataset('tod', \
                             data=fits_data[bands[i] + 'B'][time_order])
            flags_B        = detectorB_subgrp.create_dataset('flag', \
                    data = general_flags + 2**10*(np.bitwise_and(\
                    xs_noise, 2**i)/2**i).astype(int))
            healpix_B      = detectorB_subgrp.create_dataset('pix', \
                             data=hp_pix)
            time_B         = detectorB_subgrp.create_dataset('time', \
                             data=time)
            polarisation_B = detectorB_subgrp.create_dataset('polang', \
                             data=pol_theta)

            tod_C          = detectorC_subgrp.create_dataset('tod', \
                             data=fits_data[bands[i] + 'C'][time_order])
            flags_C        = detectorC_subgrp.create_dataset('flag', \
                    data = general_flags + 2**10*(np.bitwise_and(\
                    xs_noise, 2**i)/2**i).astype(int))
            healpix_C      = detectorC_subgrp.create_dataset('pix', \
                             data=hp_pix)
            time_C         = detectorC_subgrp.create_dataset('time', \
                             data=time)
            polarisation_C = detectorC_subgrp.create_dataset('polang', \
                             data=pol_theta)

        else:
            # Bands Phot04 through Phot10 each have one detector
            detectorA_subgrp = day_grp.create_group('A')

            tod              = detectorA_subgrp.create_dataset('tod', \
                               data=fits_data[bands[i]][time_order])
            flags            = detectorA_subgrp.create_dataset('flag', \
                        data = general_flags + 2**10*(np.bitwise_and(\
                        xs_noise, 2**i)/2**i).astype(int))
            healpix          = detectorA_subgrp.create_dataset('pix', \
                               data=hp_pix)
            time             = detectorA_subgrp.create_dataset('time', \
                               data=time)
            polarisation     = detectorA_subgrp.create_dataset('polang'\
                               , data=pol_theta)

# Function parsing through the CIO directory reformatting the data from
# all the single-day fits files into ten hdf5 files, one for each
# frequency band
def parse_directory():
    # The measurement days on format yyddd for years 1989 and 1990
    endings89 = np.arange(89345, 89366)
    endings90 = np.arange(90001, 90265)
    endings   = np.concatenate((endings89, endings90))

    # Loop through fits files for each day
    for i in tqdm(range(len(endings))):
        FITS_to_hdf5('/mn/stornext/d16/cmbco/ola/dirbe/cio/\
        DIRBE_CIO_P3B_' +str(endings[i]) + '.FITS', '%06d' % (i+1))
