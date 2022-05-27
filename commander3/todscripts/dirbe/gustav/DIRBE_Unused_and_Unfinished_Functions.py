import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import astropy.time as at
import h5py
import cosmoglobe
import pixfunc as pf
from astropy.io import fits
from astropy import units as u
from tqdm import tqdm
from numba import jit
from timeit import default_timer as timer
from datetime import datetime

nside = 128

# Function converting an array of QuadCube pixels into HEALpix pixels
def qs_pix2healpix(qs_pixel_numbers):
    # FITS file relating QuadCube pixel number with ecliptic longitude and latitude
    skymap_info = fits.open("/mn/stornext/d16/cmbco/ola/dirbe/auxdata/DIRBE_SKYMAP_INFO.FITS")

    # Ecliptic longitude and latitude corresponding to the QuadCube pixels
    elon = skymap_info[1].data['ELON-CSC'][qs_pixel_numbers]
    elat = skymap_info[1].data['ELAT-CSC'][qs_pixel_numbers]

    # Returning the corresponding HEALpix pixels
    return hp.ang2pix(128, elon, elat, lonlat=True)

# Function creating polarisation maps
def polarisation_mapper():
    pol_bands = np.array(('01', '02', '03'))

    bQ      = np.zeros(12*nside**2)
    bU      = np.zeros(12*nside**2)
    cQ      = np.zeros(12*nside**2)
    cU      = np.zeros(12*nside**2)

    bQhat   = np.zeros((len(pol_bands), 12*nside**2))
    bUhat   = np.zeros((len(pol_bands), 12*nside**2))
    cQhat   = np.zeros((len(pol_bands), 12*nside**2))
    cUhat   = np.zeros((len(pol_bands), 12*nside**2))
    """
    bQhat   = np.zeros(12*nside**2)
    bUhat   = np.zeros(12*nside**2)
    cQhat   = np.zeros(12*nside**2)
    cUhat   = np.zeros(12*nside**2)
    """

    for i, band in enumerate(pol_bands):

        if band != '01':
            break

        h5file = h5py.File('dirbe_hdf5_files/Phot' + band + '_128.hdf5', 'r')
        chunks = np.char.mod('%06d', np.arange(1, len(h5file)+1))
        #chunks = np.char.mod('%06d', np.arange(202, 209))

        for chunk in tqdm(chunks):
            B     = h5file['/' + chunk + '/B/tod'][()]
            C     = h5file['/' + chunk + '/C/tod'][()]
            pix   = h5file['/' + chunk + '/A/pix'][()]
            theta = h5file['/' + chunk + '/A/polang'][()]

            unipix = np.unique(pix)

            bQhat, bUhat, cQhat, cUhat = pol_calc(B, C, pix, theta, unipix, bQ, bU, cQ, cU, bQhat, bUhat, cQhat, cUhat, i)

    result = np.array((bQhat, bUhat, cQhat, cUhat))
    np.save('npy_files/pol_map_full_band1_gen1.npy', result)

#@jit#(nopython=True)
def pol_calc(B, C, pix, theta, unipix, bQ, bU, cQ, cU, bQhat, bUhat, cQhat, cUhat, i):
    #bQ = np.zeros(12*128**2)
    #bU = np.zeros(12*128**2)
    #cQ = np.zeros(12*128**2)
    #cU = np.zeros(12*128**2)

    #unipix = np.unique(pix)

    for p in tqdm(unipix):

        inds = (p == pix) & (B > -1e4)
        if ~np.any(inds):
            continue

        M      = np.zeros((2,2))
        M[0,0] = np.sum(np.cos(2*theta[inds])**2)
        M[1,1] = np.sum(np.sin(2*theta[inds])**2)
        M[0,1] = np.sum(np.sin(2*theta[inds])*np.cos(theta[inds]))
        M[1,0] = M[0,1]

        bQ[p]      = np.sum(B[inds]*np.cos(2*theta[inds]))
        bU[p]      = np.sum(B[inds]*np.sin(2*theta[inds]))
        b          = np.array((bQ[p], bU[p]))
        Mb         = np.linalg.solve(M, b)
        bQhat[i,p] = Mb[0]
        bUhat[i,p] = Mb[1]

        cQ[p]      = np.sum(C[inds]*np.cos(2*theta[inds]))
        cU[p]      = np.sum(C[inds]*np.sin(2*theta[inds]))
        c          = np.array((cQ[p], cU[p]))
        Mc         = np.linalg.solve(M, c) #lstsq
        cQhat[i,p] = Mc[0]
        cUhat[i,p] = Mc[1]

    return bQhat, bUhat, cQhat, cUhat

def binmap_plotter(band, detector, plot_count=None):
    tot_int = np.load('npy_files/binmap' + band + detector + '.npy')
    hits = np.load('npy_files/hitsmap' + band + detector + '.npy')

    hits[hits == 0] = 1

    if (plot_count == None):
        plt.figure(1, figsize=(15,9))
        hp.mollview(tot_int/hits, fig=1, title=wavelength(band) + 'um', coord=['E', 'G'], norm='hist')
        plt.savefig('python_plots/single_binmaps/binmap_' + wavelength(band) + detector + '_' + str(nside) + '_gen1.png', bbox_inches='tight', norm='hist')
        plt.show()
    else:
        hp.mollview(tot_int/hits, title=wavelength(band) + 'um', coord=["E", "G"], sub=(4,4,plot_count), norm='hist')

def pol_map_plot_all():
    bands = np.array(('01', '02', '03'))
    maps = np.load('npy_files/pol_maps.npy')
    bQhat, bUhat, cQhat, cUhat = maps
    #plt.figure(1, figsize=(15,9))
    for i, band in enumerate(bands):

        if band == '02':
            break

        hitsmapb = np.load('npy_files/hitsmap' + band + 'B.npy')
        hitsmapc = np.load('npy_files/hitsmap' + band + 'C.npy')
        #bQhat[bQhat < 0] = 0
        #bUhat[bUhat < 0] = 0
        #cQhat[cQhat < 0] = 0
        #cUhat[cUhat < 0] = 0
        hp.mollview(bQhat[i]/hitsmapb, fig=1, title=wavelength(band) + 'um' + ' bQ', norm='hist')#, sub=(4, 3, i+1), max=0.12, min=-0.12)
        hp.mollview(bUhat[i]/hitsmapb, fig=2, title=wavelength(band) + 'um' + ' bU', norm='hist')#, sub=(4, 3, i+4), max=0.12, min=-0.12)
        hp.mollview(cQhat[i]/hitsmapc, fig=3, title=wavelength(band) + 'um' + ' cQ', norm='hist')#, sub=(4, 3, i+7), max=0.12, min=-0.12)
        hp.mollview(cUhat[i]/hitsmapc, fig=4, title=wavelength(band) + 'um' + ' cU', norm='hist')#, sub=(4, 3, i+10), max=0.12, min=-0.12)
    #plt.savefig('python_plots/pol_maps_full_band1_gen1.png', bbox_inches='tight')
    plt.show()

def binmap_all_plotter(save=False):
    pol_bands  = np.array(('01', '02', '03'))
    bands      = np.char.mod('%02d', np.arange(1,11))
    detectors  = np.array(('A', 'B', 'C'))
    plot_count = 1
    plt.figure(figsize=(15,9))
    for band in bands:
        if (band in pol_bands):
            for detector in detectors:
                binmap_plotter(band, detector, plot_count)
                plot_count += 1
        else:
                binmap_plotter(band, 'A', plot_count)
                plot_count += 1
    if save==True:
        plt.savefig('python_plots/4by4_binmaps_' + str(nside) + '_gen1.png', bbox_inches='tight')

def difference_maps():
    detectors = np.array(('01A', '04A', '10A'))
    dirbe_bandnr = np.array(('1', '4', '10'))

    rotator = hp.Rotator(coord=['G', 'E'])
    plt.figure(figsize=(15,9))
    for i, detector in enumerate(detectors):
        my_map     = hp.fitsfunc.read_map('fits_maps/binmap' + detector + '_gen1.fits')

        dirbe_map_256 = hp.fitsfunc.read_map('/mn/stornext/d16/cmbco/ola/dirbe/DIRBE_' + dirbe_bandnr[i] + '_256.fits')
        dirbe_map_256 = rotator.rotate_map_alms(dirbe_map_256)
        dirbe_map_128 = hp.pixelfunc.ud_grade(dirbe_map_256, 128)

        diff_map = dirbe_map_128 - my_map

        hp.mollview(dirbe_map_128, sub=(3,3,i+1), norm='hist', title='dirbe ' + wavelength(detector[:-1]) + 'um map', coord=['E', 'G'])
        hp.mollview(my_map, sub=(3,3,i+4), norm='hist', title='my ' + wavelength(detector[:-1]) + 'um map', coord=['E', 'G'])
        hp.mollview(diff_map, sub=(3,3,i+7), norm='hist', title='difference map (dirbe - mine)', coord=['E', 'G'])
    plt.savefig('python_plots/3by3_1410_diff_map_gen1.png', bbox_inches='tight')
    plt.show()

def Y_plot(band):
    I_a = hp.fitsfunc.read_map('fits_maps/binmap' + band + 'A_gen1.fits')
    I_b = hp.fitsfunc.read_map('fits_maps/binmap' + band + 'B_gen1.fits')
    I_c = hp.fitsfunc.read_map('fits_maps/binmap' + band + 'C_gen1.fits')

    I_b[I_a == 0] = 0
    I_c[I_a == 0] = 0
    I_a[I_a == 0] = 1

    Y_b = 1 - I_b/I_a
    Y_c = I_c/I_a - 1

    Y_b = hp.smoothing(Y_b, fwhm=np.pi/180)
    Y_c = hp.smoothing(Y_c, fwhm=np.pi/180)

    if band=='01':
        minimum = 0.1; maximum = 0.4
    elif band=='02':
        minimum = 0.2; maximum = 0.75
    elif band=='03':
        minimum = 0.55; maximum = 0.9

    plt.figure(1, figsize=(15,9))
    #plt.title('Y polarization plots')
    hp.mollview(Y_b, sub=(2,1,1), title='Y_b = 1 - I_b/I_a, ' + wavelength(band) + ' um', norm='hist', fig=1)#, min=minimum, max=maximum)
    hp.mollview(Y_c, sub=(2,1,2), title='Y_c = I_c/I_a - 1, ' + wavelength(band) + ' um', norm='hist', fig=1)#, min=minimum, max=maximum)#min=-maximum, max=-minimum)
    #plt.savefig('python_plots/Y_plots_'+ band + '_smoothed_symmetric_gen1.png', bbox_inches='tight')
    plt.show()

def plot_hitsmap():
    hitsmap = np.load('npy_files/hitsmap01A.npy')
    cosmoglobe.plot(hitsmap)#, coord=['E', 'G'], rlabel='\mathrm{DIRBE Hitsmap}')
    plt.show()

if __name__=="__main__":
    #FITS_to_hdf5('/mn/stornext/d16/cmbco/ola/dirbe/cio/DIRBE_CIO_P3B_89345.FITS', 1)
    #parse_directory()
    #file = h5py.File('/mn/stornext/d16/cmbco/bp/gustavbe/master/dirbe_hdf5_files/Phot10_' + str(nside) + '.hdf5', 'r')

    #full_data_binmap(file, '10', 'A')

    #binmap_plotter('10', 'A')
    #binmap_plotter('08', 'A')
    #binmap_plotter('09', 'A')
    #plt.show()

    #TT_plot('10', 'A')
    #plt.show()

    """
    start = timer()
    binmap_all_bands()
    end = timer()
    print('Binmapping all bands and detectors took %.2f minutes' % ((end - start)/60))
    """

    #binmap_all_plotter(save=False)
    #plt.show()

    #polarisation_mapper()

    #pol_map_plot_all()

    #write_maps()

    #difference_maps()

    diff_map_all()

    #TT_plot('03')
    #TT_plot_all()

    #Y_plot('01')
    #pretty_Y_plot('02', 'C')
    #pretty_Y_plot_all(save=True)

    #pol_map_plot_all()

    #binmap_pretty_plotter('01', 'A', 164)
    #plt.show()
    #binmap_pretty_plot_all(save=True)

    #plot_hitsmap()

    #print(np.load('npy_files/binmap01B.npy') - np.load('npy_files/binmap02A.npy'))
    #print(h5py.File('dirbe_hdf5_files/Phot01.hdf5', 'r')['/000001/B/tod'][()] - h5py.File('dirbe_hdf5_files/Phot02.hdf5', 'r')['/000001/A/tod'][()])
