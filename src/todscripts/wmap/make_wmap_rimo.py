'''
Make a WMAP_instrument.h5 file in /mn/stornext/d16/cmbco/bp/dwatts/wmap/data

compare to the LFI RIMO in
/mn/stornext/d16/cmbco/bp/dwatts/wmap/data/LFI_instrument_v4.h5

Top level includes the averaged channel name, plus the different horns,
[18--28](M/S).

030/040/070 include bandpass and bandpassx, which are B(nu) and nu,
respectively.

Each horn contains:
bandpass                 Dataset {251}
bandpassx                Dataset {251}
beam                     Group
beamlmax                 Dataset {1}
beammmax                 Dataset {1}
centFreq                 Dataset {1}
elip                     Dataset {1}
fwhm                     Dataset {1}
mbeam_eff                Dataset {1}
psi_ell                  Dataset {1}
sl                       Group
sllmax                   Dataset {1}
slmmax                   Dataset {1}


Mathew points out that beam and sl are both alm representations of the beam and sidelobes.

mbeam_eff is main beam efficiency, assume it is one.
'''
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import h5py

import healpy as hp

from scipy.integrate import trapz
from scipy.optimize import curve_fit, minimize

import reproject


def gauss(x, sigma2, A):
    if sigma2 < 0:
        return np.inf
    else:
        return A*np.exp(-x**2/(2*sigma2))


fname_out = '/mn/stornext/d16/cmbco/bp/dwatts/wmap/data/WMAP_instrument_v1.h5'
labels = ['K', 'Ka', 'Q', 'V', 'W']



with h5py.File(fname_out, 'a') as f:
    # Bandpasses
    fnames = glob('data/wmap_bandpass_*_v5.cbp')
    fnames.sort()
    
    
    totals = [[],[],[],[],[]]
    for i in range(len(fnames)):
        band = fnames[i].split('_')[2]
        nu, B1, B2 = np.loadtxt(fnames[i]).T
        f.create_dataset(band + '3/bandpassx', data=nu)
        f.create_dataset(band + '4/bandpassx', data=nu)
        f.create_dataset(band + '3/bandpass', data=B1)
        f.create_dataset(band + '4/bandpass', data=B2)
        centFreq1 = trapz(nu*B1, nu)/trapz(nu, B1)
        centFreq2 = trapz(nu*B2, nu)/trapz(nu, B2)
        f.create_dataset(band + '3/centFreq', data=[centFreq1])
        f.create_dataset(band + '4/centFreq', data=[centFreq2])
    
        for j in range(len(labels)):
            if labels[j] == band[:-2]:
                totals[j].append([nu, B1])
                totals[j].append([nu, B2])
    
    for i in range(len(labels)):
        bp = np.array(totals[i])
        nu, B = np.mean(bp, axis=0)
        f.create_dataset(labels[i] + '/bandpassx', data=nu)
        f.create_dataset(labels[i] + '/bandpass', data=B)

    # FWHMs
    ## From radial beam profiles
    DAs = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
    fnames = glob('data/wmap_symm_beam_profile_*_9yr_v5.txt')
    fnames.sort()
    for fname in fnames:
        theta, B = np.loadtxt(fname).T
        popt, pcov = curve_fit(gauss, theta, B, p0=[1, 10])
        sigma = popt[0]**0.5
        fwhm = 2*np.sqrt(2*np.log(2))*sigma
    
        DA = fname.split('_')[4]
        f.create_dataset(DA + '13/fwhm', data=[fwhm])
        f.create_dataset(DA + '14/fwhm', data=[fwhm])
        f.create_dataset(DA + '23/fwhm', data=[fwhm])
        f.create_dataset(DA + '24/fwhm', data=[fwhm])
        """
        Beam ellipticity: WMAP doesn't report the ellipticity directly. They describe
        the beam-symmetrization map-making process in http://arxiv.org/abs/1212.5225.
        The beams for a single pointing are fairly elliptical and non-Gaussian.
        
        For a first pass, I will set these to zero.
        """
        f.create_dataset(DA + '13/elip', data=[0])
        f.create_dataset(DA + '14/elip', data=[0])
        f.create_dataset(DA + '23/elip', data=[0])
        f.create_dataset(DA + '24/elip', data=[0])





fnames = glob('data/wmap_hybrid_beam_maps_*_9yr_v5.fits')
fnames.sort()


hdus = [fits.open(fname)[0] for fname in fnames]
wcs = WCS(hdus[0].header) # hdus are the same


# should be at least 1024, since pixel size is 0.04 ~ 58.6/nside
res_min = np.log(58.6/0.04)/np.log(2)
nside_min = 2**np.ceil(res_min)
nside_beam = 1024

nside_beam = 2048


lmax = 2400
mmax = 100


# I am nearly certain that the projection that they use, 
# X=2*sin(theta/2)*cos(phi), is zenithal equal area, making the coordinate
# system centered at the north pole.
target_header = fits.Header.fromstring("""
NAXIS   =                    2
NAXIS1  =                  600
NAXIS2  =                  600
CTYPE1  = 'RA---ZEA'
CRPIX1  =                0.5
CRVAL1  =                -11.98
CDELT1  =               0.04
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--ZEA'
CRPIX2  =                0.5
CRVAL2  =                -11.98
CDELT2  =                0.04
CUNIT2  = 'deg     '
COORDSYS= 'icrs    '
""", sep='\n')

for fname in fnames:
    data = fits.open(fname)
    
    m_A, footprint_A = reproject.reproject_to_healpix((data[0].data[0], target_header), 'C', 
                                                   nside=nside_beam,
                                                   order=3)
    m_B, footprint_B = reproject.reproject_to_healpix((data[0].data[2], target_header), 'C', 
                                                   nside=nside_beam,
                                                   order=3)
    m_A[footprint_A==0] = 0
    m_B[footprint_B==0] = 0

    r = hp.rotator.Rotator(rot=hp.pix2ang(nside_beam, np.where(m_A==m_A.max())[0][0],
        lonlat=True))
    m_A_rotated = r.rotate_map_pixel(m_A)
    r2 = hp.rotator.Rotator(rot=(0,-90,0))
    m_A_rotated = r2.rotate_map_pixel(m_A_rotated)

    alm_A = hp.map2alm(m_A_rotated, lmax=lmax)

    b_lm_A = np.zeros((lmax+1)**2)
    for i in range(len(alm_A)):
        li, mi = hp.sphtfunc.Alm.getlm(lmax,i=i)
        if (li <= lmax) & (mi <= mmax):
            ind_real = li**2 + li + mi
            ind_imag = li**2 + li - mi
            b_lm_A[ind_real] = alm_A[i].real
            if mi != 0:
                b_lm_A[ind_imag] = alm_A[i].imag
               





    r = hp.rotator.Rotator(rot=hp.pix2ang(nside_beam, np.where(m_B==m_B.max())[0][0],
        lonlat=True))
    m_B_rotated = r.rotate_map_pixel(m_B)
    r2 = hp.rotator.Rotator(rot=(0,-90,0))
    m_B_rotated = r2.rotate_map_pixel(m_B_rotated)

    alm_B = hp.map2alm(m_B_rotated)
    b_lm_B = np.zeros((lmax+1)**2)
    for i in range(len(alm_B)):
        li, mi = hp.sphtfunc.Alm.getlm(lmax,i=i)
        if (li <= lmax) & (mi <= mmax):
            ind_real = li**2 + li + mi
            ind_imag = li**2 + li - mi
            b_lm_B[ind_real] = alm_B[i].real
            if mi != 0:
                b_lm_B[ind_imag] = alm_B[i].imag

    DA = fname.split('_')[4]

    with h5py.File(fname_out, 'a') as f:
        f.create_dataset(DA + '13/beam/T', data=b_lm_A)
        f.create_dataset(DA + '23/beam/T', data=b_lm_A)
        f.create_dataset(DA + '14/beam/T', data=b_lm_B)
        f.create_dataset(DA + '24/beam/T', data=b_lm_B)
        f.create_dataset(DA + '13/beamlmax', data=[lmax])
        f.create_dataset(DA + '23/beamlmax', data=[lmax])
        f.create_dataset(DA + '14/beamlmax', data=[lmax])
        f.create_dataset(DA + '24/beamlmax', data=[lmax])
        f.create_dataset(DA + '13/beammmax', data=[mmax])
        f.create_dataset(DA + '23/beammmax', data=[mmax])
        f.create_dataset(DA + '14/beammmax', data=[mmax])
        f.create_dataset(DA + '24/beammmax', data=[mmax])
        
        f.create_dataset(DA + '13/sl/T', data=b_lm_A*0)
        f.create_dataset(DA + '23/sl/T', data=b_lm_A*0)
        f.create_dataset(DA + '14/sl/T', data=b_lm_B*0)
        f.create_dataset(DA + '24/sl/T', data=b_lm_B*0)
        f.create_dataset(DA + '13/sllmax', data=[lmax])
        f.create_dataset(DA + '23/sllmax', data=[lmax])
        f.create_dataset(DA + '14/sllmax', data=[lmax])
        f.create_dataset(DA + '24/sllmax', data=[lmax])
        f.create_dataset(DA + '13/slmmax', data=[mmax])
        f.create_dataset(DA + '23/slmmax', data=[mmax])
        f.create_dataset(DA + '14/slmmax', data=[mmax])
        f.create_dataset(DA + '24/slmmax', data=[mmax])


        f.create_dataset(DA + '13/mbeam_eff', data=[1])
        f.create_dataset(DA + '23/mbeam_eff', data=[1])
        f.create_dataset(DA + '14/mbeam_eff', data=[1])
        f.create_dataset(DA + '24/mbeam_eff', data=[1])
        f.create_dataset(DA + '13/psi_ell', data=[0])
        f.create_dataset(DA + '23/psi_ell', data=[0])
        f.create_dataset(DA + '14/psi_ell', data=[0])
        f.create_dataset(DA + '24/psi_ell', data=[0])


# a_lm indices
#data = h5py.File('/mn/stornext/d16/cmbco/bp/dwatts/wmap/data/LFI_instrument_v4.h5', 'r')
#lmax = data['18M/beamlmax'][...][0]
#T_lm = data['18M/beam/T'][...]
#T_lmcpx = np.zeros(hp.sphtfunc.Alm.getidx(lmax,lmax,lmax)+1, dtype='complex')
#inds = {}
#for li in range(lmax+1):
#    for mi in range(-li, li+1):
#        inds[li**2+li+mi] = (li,mi)
#for i in range(len(inds)):
#    li, mi = inds[i]
#    if mi >= 0:
#        idx = hp.sphtfunc.Alm.getidx(lmax, li, mi)
#        T_lmcpx[idx] += T_lm[i]
#    else:
#        idx = hp.sphtfunc.Alm.getidx(lmax, li, -mi)
#        T_lmcpx[idx] += 1j*T_lm[i]

plt.show()
"""
# far sidelobes
fnames = glob('data/wmap_sidelobe_map_*_9yr_v5.fits')
nside = 2**7
inds = hp.nest2ring(nside, np.arange(12*nside**2))
for fname in fnames:
    m = fits.open(fname)
    # NB; these are in nest format
    T = m[1].data['TEMPERATURE'][inds]
"""
