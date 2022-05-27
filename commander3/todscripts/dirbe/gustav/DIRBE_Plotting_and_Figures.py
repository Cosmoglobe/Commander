import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
#import astropy.time as at
#import h5py
import cosmoglobe
#import pixfunc as pf
#from astropy.io import fits
from astropy import units as u
#from tqdm import tqdm
#from numba import jit
#from timeit import default_timer as timer
#from datetime import datetime

# Function plotting a single binned map using the cosmoglobe library
def binmap_pretty_plotter(band, detector, rlabel, maximum):
    binmap = np.load('npy_files/binmap' + band + detector + '.npy')
    hitsmap = np.load('npy_files/hitsmap' + band + detector + '.npy')

    # Removing non-observed pixels to avoid division by zero issues
    binmap[hitsmap == 0] = 0
    hitsmap[hitsmap == 0] = 1

    mapp = binmap/hitsmap
    mapp[mapp<=0] = 0.0001  # For log-scale purposes

    wl_label = wavelength(band)

    cosmoglobe.plot(mapp, cmap='afmhot', coord=['E', 'G'], unit=r'\mathrm{MJy/sr}', \
    rlabel=rlabel, llabel=wl_label + r'\mathrm{\mu m}', norm='log', cbar=True, width=5.7, extend='neither', max=maximum)

# Function plotting all binned maps using the binmap_pretty_plotter function and saving in seperate figures
def binmap_pretty_plot_all(save=False):
    bands = np.char.mod('%02d', np.arange(1,11))
    pol_bands = np.array(('01', '02', '03'))
    detectors = np.array(('A', 'B', 'C'))

    # Colorbar maximums
    maximums = np.array((160, 110, 100, 150, 50, 40, 90, 20, 20, 100, 430, 990, 4030, 6700, 10590, 5830))

    for i,band in enumerate(bands):
        if band in pol_bands:
            for j,detector in enumerate(detectors):
                if detector == 'A':
                    rlabel='Photometry'
                else:
                    rlabel='Polarimetry\ ' + detector
                binmap_pretty_plotter(band, detector, rlabel, maximums[3*i + j])
                if save == True:
                    plt.savefig('python_plots/single_binmaps/pretty_map_' + wavelength(band) + detector + '_' + str(nside) + '_gen1.pdf', bbox_inches='tight')
        else:
            binmap_pretty_plotter(band, 'A', 'Photometry', maximums[i + 6])
            if save == True:
                plt.savefig('python_plots/single_binmaps/pretty_map_' + wavelength(band) + 'A_' + str(nside) + '_gen1.pdf', bbox_inches='tight')
        plt.close()

# Function creating and saving relative difference maps for all DIRBE bands
def diff_map_all():
    detectors = np.array(('01A', '02A', '03A', '04A', '05A', '06A', '07A', '08A', '09A', '10A'))
    limit = np.array((0.08, 0.11, 0.08, 0.08, 0.05, 0.03, 0.03, 0.03, 0.16, 0.14))

    rotator = hp.Rotator(coord=['G', 'E'])
    for i, detector in enumerate(detectors):
        dirbe_original_map_256 = hp.fitsfunc.read_map('/mn/stornext/d16/cmbco/ola/dirbe/DIRBE_' + str(i+1) + '_256.fits')
        dirbe_original_map_256 = rotator.rotate_map_alms(dirbe_original_map_256)
        dirbe_original_map_128 = hp.pixelfunc.ud_grade(dirbe_original_map_256, 128)

        my_binmap = hp.fitsfunc.read_map('fits_maps/binmap' + detector + '_gen1.fits')

        diff_map = dirbe_original_map_128 - my_binmap

        # Removing non-observed pixels to avoid division by zero issues
        diff_map[diff_map==0] = hp.UNSEEN

        # Plotting relative difference
        cosmoglobe.plot(diff_map/dirbe_original_map_128, fwhm=1*u.deg, coord=['E', 'G'], llabel=wavelength(detector[:-1]) + r'\mathrm{\mu m}', rlabel=r'\epsilon', max=limit[i], min=-limit[i], width=5.7)#, norm='hist')
        plt.savefig('python_plots/single_diff_maps/difference_map_' + wavelength(detector[:-1]) + '_' + str(nside) + '_gen1.pdf', bbox_inches='tight')
        plt.close()
    #plt.show()

# Function creating an average polarization signal map for a single wavelength
def pretty_Y_plot(band, detector):
    I_a = hp.fitsfunc.read_map('fits_maps/binmap' + band + 'A_gen1.fits')
    I_pol = hp.fitsfunc.read_map('fits_maps/binmap' + band + detector + '_gen1.fits')

    # Removing non-observed pixels to avoid division by zero issues
    I_pol[I_a == 0] = 0
    I_a[I_a == 0] = 1

    # Handling the two plarization channels
    if detector == 'B':
        Y = 1 - I_pol/I_a
    else:
        Y = I_pol/I_a - 1

    # Maps are smoothed to 1 degree FWHM
    Y = hp.smoothing(Y, fwhm=np.pi/180)

    if band == '01':
        minimum = 0.1; maximum = 0.4
    elif band == '02':
        minimum = 0.2; maximum = 0.75
    elif band == '03':
        minimum = 0.55; maximun = 0.9

    cosmoglobe.plot(Y, coord=['E', 'G'], rlabel='Y_' + detector, llabel=wavelength(band) + r'\mathrm{\mu m}', width=5.7)

# Function plotting all Y maps and saving as seperate figures
def pretty_Y_plot_all(save=False):
    bands = np.array(('01', '02', '03'))
    detectors = np.array(('B', 'C'))

    for band in bands:
        for detector in detectors:
            pretty_Y_plot(band, detector)
            if save==True:
                plt.savefig('python_plots/pol_maps/Y_map_' + wavelength(band) + detector + '_' + str(nside) + '_gen1.pdf')
    plt.show()

# Function creating a TT plot with linear regression best-fit line for pixel temperatures
def TT_plot(band):
    my_map = hp.fitsfunc.read_map('fits_maps/binmap' + band + 'A.fits')

    # Rotating and down grading resolution for CADE map
    rotator = hp.Rotator(coord=['G', 'E'])
    dirbe_map_256 = hp.fitsfunc.read_map('/mn/stornext/d16/cmbco/ola/dirbe/DIRBE_' + str(int(band)) + '_256.fits')
    dirbe_map_256 = rotator.rotate_map_alms(dirbe_map_256)
    dirbe_map_128 = hp.pixelfunc.ud_grade(dirbe_map_256, 128)

    # Linear regression
    N, cov = np.polyfit(my_map, dirbe_map_128, 1, cov=True)
    sigma = np.diag(cov)[0]**0.5
    line = np.arange(0, int(np.max(my_map)) + 100)

    plt.plot(my_map, dirbe_map_128, 'k.', alpha=0.1)
    plt.plot(line, N[0]*line, 'r--', label=r'y = (%.3f $\pm$ %.2e)x' % (N[0], sigma))
    plt.title(r'%s $\mu$m' % wavelength(band))
    plt.xlabel('CADE')
    plt.ylabel('Current')
    plt.legend()

# Function creating the TT plots for all bands using the TT_plot function, and saving in a 2by5 grid figure
def TT_plot_all():
    bands = np.char.mod('%02d', np.arange(1,11))
    plt.figure(figsize=(7,10))
    plt.title('TT Correlation Plots')
    for i, band in enumerate(bands):
        print(i+1)
        plt.subplot(5, 2, i+1)
        TT_plot(band)
    plt.tight_layout()
    plt.savefig('python_plots/TT_plots_2by5_gen1.png')
    plt.show()
