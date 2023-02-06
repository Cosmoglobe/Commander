
## Global Parameters

In the table below, the **yy** suffix indicates a smoothing scale counter. These smoothing scales may be used to estimate spectral indices at common angular resolutions for each parameter, with different resolutions for each parameter.

| Parameter name                | Description                                          | Allowed Values    | Suggested value |
| ----------------------------- | ---------------------------------------------------- | ----------------- | :-------------: |
| DATA_DIRECTORY                | Path to directory where data sets are available, relative to working directory | Any valid path | data          |
| NUMBAND                       | Maximum number of data sets defined in parameter file; additional data sets may be present, but these are ignored    | >= 1                  | NA            |
| NUM_SMOOTHING_SCALES          | Number of allowed common smoothing scales to be used for local spectral index sampling  | >= 0                  | 0             |
| PROCESSING_MASKFILE           | "Soft" processing mask in HEALPix format; for channels with noise type "rms", the actual RMS of any pixel with mask value smaller than 0.5 is increased by a factor of 20. Useful for suppressing ringing from the Galactic plane or bright sources | valid filename present in DATA_DIRECTORY. Ignored for non-diagonal noise types.                  | none          |
| SOURCE_MASKFILE               | ASCII file with locations to be masked in all channels. Each line must contain (lon, lat, radius), where (lon,lat) is the Galactic coordinates of the source in degrees, and radius is a multiplier relative to the beam FWHM for each channel   | valid filename present in DATA_DIRECTORY  | none          |
| SMOOTHING_SCALE_FWHMyy       | Common angular resolution in arcmin of smoothing scale yy used to estimate some spectral index parameter; must be larger than the largest FWHM of any frequency channel included in the analysis of the component to which this smoothing scale is applied | >0 | 120. |
| SMOOTHING_SCALE_FWHM_POSTPROCyy | Additional smoothing kernel width in arcmin, applied after performing a pixel-by-pixel spectral index fit | > 0. | 60. |
| SMOOTHING_SCALE_LMAXyy       | Maximum multipole used for degrading maps to common angular resolution; should be chosen to match SMOOTHING_SCALE_FWHM | >= 0 | 256 |
| SMOOTHING_SCALE_NSIDEyy      | HEALPix resolution parameter at which pixel-by-pixel spectral index fit is performed | > 0 | 128 |
| SMOOTHING_SCALE_PIXWINyy     | Pixel window file for target resolution | valid filename in DATA_DIRECTORY | pixel_window_n0128.fits |

## Band Definition

For each band, the following parameters are available. Parameters with
a gray background are mandatory, while parameters with a red
background are only needed if ENABLE_TOD_ANALYSIS is enabled. The
**xxx** suffix denotes an incremental counter for each band, e.g., 003
for the third data set defined in the parameter file.

| Parameter name                | Description                               | Allowed Values                | Example |
| ----------------------------- | ----------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-------: |
| BAND_BANDPASS_MODELxxx        | Parametric model type for bandpass corrections | powlaw_tilt={$\tau = \tau_0 (\nu/\nu_c)^{\beta}$}<br> additive_shift={$\tau=\tau_0(\nu+\Delta\nu)$} | additive_shift |
| BAND_BANDPASS_TYPExxx         | Bandpass definition type. Each experiment defines its bandpass ($\tau$) and unit conversion integration formulae slightly differently, adoping different units and/or thresholding for the bandpass   | delta={$\delta$ function}<br>DIRBE={same as HFI_submm, but no lower threshold on bandpass amplitude}<br>LFI={see [Planck V (2013)](https://www.aanda.org/articles/aa/full_html/2014/11/aa21527-13/aa21527-13.html)}<br>HFI_cmb={see [Planck IX (2013)](https://www.aanda.org/articles/aa/abs/2014/11/aa21531-13/aa21531-13.html)}<br>HFI_submm={see [Planck IX (2013)](https://www.aanda.org/articles/aa/abs/2014/11/aa21531-13/aa21531-13.html)}<br>PSM_LFI={alias for HFI_cmb, but no threshold on bandpass amplitude}<br>WMAP={see [Bennett et al. (2013)](https://arxiv.org/abs/1212.5225)} | LFI |
| BAND_BANDPASSFILExxx          | Data file containing bandpass file; may be either ASCII table or an HDF instrument file (with suffix 'h5' or 'hd5') | valid filename in DATA_DIRECTORY <br> none (if 'delta' bandpass type) | bp_030_v1.dat |
| BAND_BEAM_B_L_FILExxx         | Beam data file; for now, only FITS file containing Legendre transform, $b_{\ell}$, is supported | valid filename in DATA_DIRECTORY | beam_5arcmin.fits |
| BAND_BEAM_B_PTSRC_FILEXXX     | HDF data file with FEBeCoP beams per point source pixel | valid filename in DATA_DIRECTORY<br>none | febecop_030_v1.h5 |
| BAND_BEAMTYPExxx              | Beam format for main convolutions; only azimuthally symmetric beams supported for now  | b_l = {$B_{\ell m,\ell' m'} = b_{\ell} \delta_{\ell\ell'}\delta_{mm'}$} | b_l |
| BAND_COMPONENT_SENSITIVITYxxx | Specification of components to which current band is sensitive | broadband={sensitive to all components}<br>*component label*={sensitive only to specified component, e.g., line emission} | broadband |
| BAND_DEFAULT_BP_DELTAxxx      | Default bandpass correction value | any float | 0. |
| BAND_DEFAULT_GAINxxx          | Default map-level calibration factor      | > 0 | 1. |
| BAND_DEFAULT_NOISEAMPxxx      | Default noise RMS amplitude correction value | > 0 | 1. |
| BAND_GAIN_APOD_FWHMxxx        | Apodization FWHM smoothing kernel width for cross-correlation calibration mode in arcmin | >= 0 | 120. |
| BAND_GAIN_APOD_MASKxxx        | Apodization mask used for cross-correlation calibration mode      | *valid filename in DATA_DIRECTORY* = {apply apodization mask}<br>fullsky={do not apply apodization mask}  | fullsky |
| BAND_GAIN_CALIB_COMPxxx       | Calibration source                                          | all={Calibrate with full signal model}<br>*valid component label*={Calibrate with respect to specified component} | cmb |
| BAND_GAIN_LMAXxxx             | Maximum multipole for map-level calibration | <0 = {Pixel space regression}<br> >= 0 = {Cross-correlate power spectra using specified multipole range}   | -1 |
| BAND_GAIN_LMINxxx             | Minimum multipole for map-level calibration  | <0 = {Pixel space regression}<br> >= 0 = {Cross-correlate power spectra using specified multipole range} | -1 |
| BAND_LABELxxx                 | User-defined data set label; used to identify the data set both in input and output files  | Any string shorter than 16 characters | 030, WMAP_K, haslam etc. |
| BAND_LMAXxxx                  | Maximum multipole moment; should be set sufficiently large that the band is strongly noise dominated at $\ell_{\mathrm{max}}$ to avoid ringing   | $0\le \ell_{\mathrm{max}} \le 4N_{\mathrm{side}}$ | $3N_{\mathrm{side}}$|
| BAND_MAPFILExxx               | Data sky map in HEALPix format                     | valid filename in DATA_DIRECTORY | bp_030_map_v1.fits |
| BAND_MASKFILExxx              | Main mask file; pixels with value smaller than 0.5 are given zero statistical weight      | valid filename in DATA_DIRECTORY<br>fullsky | mask_common_v1.fits |
| BAND_MASKFILE_CALIBxxx        | Mask file used for absolute map-level calibration; is internally multiplied with BAND_MASKFILExxx    | valid filename in DATA_DIRECTORY<br>fullsky | fullsky |
| BAND_NOISE_FORMATxxx          | Noise matrix model     | rms={$N_{pp'}^{ss'}=(\sigma^s_p)^2 \delta_{pp'}\delta_{ss'}$} where $s$ denotes pixel and $s$ Stokes parameter  <br>QU_cov = {dense $N_{pp'}^{ss'}$} with $s$={Q,U} (eg., WMAP low-res matrix) |  rms |
| BAND_NOISE_UNIFORMIZE_FSKYxxx | If > 0, white noise is added to the specified sky fraction of pixels with the lowest noise RMS, increasing their noise level to those at $f_{\mathrm{sky}}$. Useful to improve CG convergence for data sets with particularly strong RMS map contrasts, such as Planck | $0 \le f_{\mathrm{sky}} \le 1$ | 0 |
| BAND_NOISEFILExxx             | Pixel-by-pixel noise specification; must match format chosen in BAND_NOISE_FORMAT. If 'diagonal', then this must be a HEALPix sky map. If 'QU_cov', it must be a WMAP-style low-res covariance matrix. (In this case, $N^{-1}$ and $N^{-1/2}$ are pre-computed and stored in "BAND_NOISEFILE_precomp.unf"; later runs will use these to save initialization time.) | valid filename in DATA_DIRECTORY | bp_030_rms_v1.fits |
| BAND_REG_NOISEFILExxx         |        | valid filename in DATA_DIRECTORY <br> none| bp_030_rmsreg_v1.fits |
| BAND_NOISE_RMSxxx_SMOOTHyy    | Noise RMS map for smoothing scale yy; used for smoothed spectral index sampling | valid filename in DATA_DIRECTORY <br> none| bp_030_rms_60arcmin_v1.fits |
| BAND_NSIDExxx                 | HEALPix resolution parameter | >= 1 | 512 |
| BAND_NOMINAL_FREQxxx          | Nominal frequency of current band in GHz                           | >= 0 | 70.1 |
| BAND_OBS_PERIODxxx            | Not active yet; intended to eventually support analysis of time-variable sources | NA | NA|
| BAND_PIXEL_WINDOWxxx          | Pixel window for current band (typically a HEALPix pixel window file)   | valid filename in DATA_DIRECTORY | pixel_window_n0512.fits |
| BAND_POLARIZATIONxxx          | Polarization switch   | .true.<br>.false. | .true. |
| BAND_SAMP_BANDPASSxxx         | Enable sampling of bandpass correction parameters              | .true.<br>.false. | .false. |
| BAND_SAMP_GAINxxx             | Enable sampling of map-level absolute calibration  | .true.<br>.false. | .true. |
| BAND_SAMP_NOISE_AMPxxx        | Not implemented yet                       | .true.<br>.false. | .false. |
| BAND_UNITxxx                  | Physical unit of input data and rms maps  | uK_cmb, mK_cmb, K_cmb<br>MJy/sr<br>uK_RJ<br>K km/s | uK_cmb |
| BAND_TOD_BP_INIT_PROPxxx      | Configuration file for detector-level bandpass corrections | valid filename in DATA_DIRECTORY | bp_init_030_v1.dat |
| BAND_TOD_DETECTOR_LISTxxx     | Detector labels to be included in current data sets         | Comma-separated list of valid detector labels | "27M,27S,28M,28S" |
| BAND_TOD_FILELISTxxx          | TOD scan definition file, defining periods of data with assumed internally stationary properties  | valid filename in DATA_DIRECTORY | filelist_30_v14.txt |
| BAND_TOD_INIT_FROM_HDFxxx     | Initialization mode; if a previous chain file is specified, it must follow the same convention as INIT_CHAIN         | none={use default values}<br>default={use default chain specified in INIT_CHAIN}<br>*previous chain file*={use specified sample} | "chain_init_v2.h5:5" |
| BAND_TOD_MAIN_PROCMASKxxx     | Main TOD processing mask; used for correlated noise, gain and $\chi^2$ estimation | valid filename in DATA_DIRECTORY | mask_proc_v1.fits |
| BAND_TOD_SMALL_PROCMASKxxx    | Secondary TOD processing mask; used for bandpass estimation | valid filename in DATA_DIRECTORY | mask_brightsources_v1.fits |
| BAND_TOD_START_SCANIDxxx      | First scan ID to be included in analysis; may be adjusted during testing and debugging         | > 0 | 3 |
| BAND_TOD_END_SCANIDxxx        | Last scan ID to be included in analysis;  may be adjusted during testing and debugging         | >= BAND_TOD_START_SCANID   | 45860 |
| BAND_TOD_TOT_NUMSCANxxx       | Total number of scan IDs. This defines the array sizes in output chain files. To support chain initialization using previous runs, this should ideally not change between runs | >=  BAND_TOD_END_SCANID | 45860 |
| BAND_TOD_FLAGxxx              | Bit-encoded flag for which data to discard. Any set bits that much the flag are discarded. For example, if bits 1, 3, 4, 5, 7 and 9 are set, then for LFI, the corresponding flag is 10111010100000000000000, which is recorded as 6111232 in decimal.                        | >= 0          | 6111232
| BAND_TOD_RIMOxxx              | Commander-style reduced instrument model file in HDF format   | valid filename in DATA_DIRECTORY | LFI_instrument_v4.h5 |
| INCLUDE_BANDxxx               | Enable channel in analysis                                        | .true.<br>.false. |  .true.             |
