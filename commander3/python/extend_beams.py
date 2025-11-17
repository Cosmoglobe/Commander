import healpy as hp
import numpy as np



def extend_beam(bl,fwhm,lmax_bl,lmax_trans,lmax_gauss):
    '''
    bl - the original beam function
    fwhm - the fwhm of the extended beam Gaussian in arcmin
    lmax_bl - the lmax from the original beam 
    lmax_trans - the lmax of the transition from the original beam to the Gaussian beam
    lmax_gauss - the lmax of the extended beam 
    '''
    #generate the gaussian beam and empty array for the extended beam
    gl=hp.gauss_beam(np.deg2rad(fwhm/60),lmax=lmax_gauss)
    bl_ext=np.zeros(lmax_gauss+1)
    #rescale so that the start of the Gaussian beam is the same value as the end of the original beam
    gl=(gl/(gl[lmax_bl]))*bl[lmax_bl]
    #set the lowest ell equal to the original beam, and the highest ell to a pure gaussian beam
    bl_ext[0:lmax_bl]=bl[0:lmax_bl]
    bl_ext[lmax_trans:]=gl[lmax_trans:]
    #calculate and apply apodization function between orginal and fully Gaussian beams
    ell=np.arange(lmax_bl,lmax_trans)
    apod=0.5*(1+np.cos(np.pi*(lmax_trans-ell)/(lmax_trans-lmax_bl)))
    bl_ext[lmax_bl:lmax_trans]=gl[lmax_bl:lmax_trans]*apod+bl[lmax_bl:lmax_trans]*(1-apod)

    bl_ext/=bl_ext[0]
    return bl_ext

freq=np.array([100,143,217,353,545,857])
freq_arr=[[1,2,3,4],
          [1,2,3,4,5,6,7],
          [1,2,3,4,5,6,7,8],
          [1,2,3,4,5,6,7,8],
          [1,2,4],
          [1,2,3,4]]

fwhm=np.array([9.659, 7.22,  4.9,   4.916, 4.675, 4.216]) ##from the FWHMGAUS in the HFI_RIMO_R3.00.fits
lmax_trans=np.array([2200,3000,5200,5500,4300,4300]) ##chosen by eye such that the original beam is still smooth
lmax_bl=lmax_trans-100 
lmax_gauss=10000

folder='/mn/stornext/u3/raelynsu/d23hfi/data/Bl_npipe6v20_'
#353-7x353-7.fits'

for i in range(len(freq)):
    for j in range(len(freq_arr[i])):
        strv=str(freq[i])+'-'+str(freq_arr[i][j])
        bl=hp.read_cl(folder+strv+'x'+strv+'.fits')
        bl_ext=extend_beam(bl,fwhm[i],lmax_bl[i],lmax_trans[i],lmax_gauss)
        hp.write_cl(folder+strv+'x'+strv+'_extended.fits',bl_ext,overwrite=True)


