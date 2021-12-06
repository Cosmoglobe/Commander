# File Formats

In addition to the main Commander parameter file, the various Commander inputs may grouped into four main categories of files:
- *High-level data files*: Files used during component separation. Examples include frequency and noise maps; beams and pixel windows; and bandpasses.
- *Low-level data files*: Files used during TOD processing. Examples include uncalibrated compressed TOD files; TOD filelists; and instrument description files.
- *Model-specific files*: Files used to describe individual components. Examples include SED templates; angular power spectrum files; and point source catalogs.
- *MCMC tuning files*: Files used to tune individual Commander sampling steps. Examples include bandpass or spectral index MCMC proposal files.

All necessary files are described in the present section.

## High-level data files

### Sky maps

Sky maps contain a pixelized representation of the sky as seen by a given instrument, including spatial smoothing (ie., beam convolution), frequency averaging (ie., bandpass integration), and pixelization (ie., map making). In Commander, such sky maps are defined in terms of the HEALPix pixelization, and stored as binary FITS files. For full specification of these, see Section 2.1 in the [HEALPix file format specification](https://healpix.sourceforge.io/data/examples/healpix_fits_specs.pdf).

Final Commander products, as produced by [c3pp](https://www.github.com/cosmoglobe/c3pp), conforms to the same header conventions as official [Planck](https://pla.esac.esa.int) sky maps. Basic information such as FITS column content, units, reference frequency, etc. is specified in the FITS header.

Example:
```
$ fitsinfo BP_030_IQU_n0512_v1.fits
2 HDUS:
HDU 1: Image
    SIMPLE  =                    T / conforms to FITS standard
    BITPIX  =                    8 / array data type
    NAXIS   =                    0 / number of array dimensions
    EXTEND  =                    T

HDU 2: Binary table
    XTENSION= 'BINTABLE'           / binary table extension
    BITPIX  =                    8 / array data type
    NAXIS   =                    2 / number of array dimensions
    NAXIS1  =                36864 / length of dimension 1
    NAXIS2  =                 3072 / length of dimension 2
    PCOUNT  =                    0 / number of group parameters
    GCOUNT  =                    1 / number of groups
    TFIELDS =                    9 / number of table fields
    TTYPE1  = 'I_MEAN  '
    TFORM1  = '1024E   '
    TUNIT1  = 'uK      '
    TTYPE2  = 'Q_MEAN  '
    TFORM2  = '1024E   '
    TUNIT2  = 'uK      '
    TTYPE3  = 'U_MEAN  '
    TFORM3  = '1024E   '
    TUNIT3  = 'uK      '
    TTYPE4  = 'I_RMS   '
    TFORM4  = '1024E   '
    TUNIT4  = 'uK      '
    TTYPE5  = 'Q_RMS   '
    TFORM5  = '1024E   '
    TUNIT5  = 'uK      '
    TTYPE6  = 'U_RMS   '
    TFORM6  = '1024E   '
    TUNIT6  = 'uK      '
    TTYPE7  = 'I_STDDEV'
    TFORM7  = '1024E   '
    TUNIT7  = 'uK      '
    TTYPE8  = 'Q_STDDEV'
    TFORM8  = '1024E   '
    TUNIT8  = 'uK      '
    TTYPE9  = 'U_STDDEV'
    TFORM9  = '1024E   '
    TUNIT9  = 'uK      '
    PIXTYPE = 'HEALPIX '           / HEALPIX pixelisation.
    ORDERING= 'RING    '           / Pixel ordering scheme, either RING or NESTED
    COORDSYS= 'GALACTIC'           / Ecliptic, Galactic or Celestial (equatorial)
    EXTNAME = 'xtension'           / name of this binary table extension
    NSIDE   =                  512 / Resolution parameter of HEALPIX
    FIRSTPIX=                    0 / First pixel # (0 based)
    LASTPIX =              3145727 / Last pixel # (0 based)
    INDXSCHM= 'IMPLICIT'           / Indexing: IMPLICIT or EXPLICIT
    OBJECT  = 'FULLSKY '           / Sky coverage, either FULLSKY or PARTIAL
    DATE    = 'Written Sun Nov  1 17:59:25 2020' / Time and date of creation.
    POLAR   =                    T
    BAD_DATA=          -1.6375E+30 / HEALPIX UNSEEN value.
    METHOD  = 'COMMANDER'          / COMMANDER sampling framework
    AST-COMP= '030     '
    FREQ    = '30.0 GHz'
    PROCVER = 'BP8     '           / Release version
    FILENAME= 'BP_030_IQU_full_n0512_BP8.fits'
    BNDCTR  =                   30 / Formal Band Center
    RESTFREQ=               28.456 / Effective Central Frequency
    BNDWID  =    9.898999999999999 / Effective Bandwidth
```

### Noise specification

The statistical properties of the instrumental noise present in the sky maps must be characterized in terms of its covariance structure; its mean is always assumed to be zero. Commander currently supports two different types of noise specifications:

| Noise type | Description | Detailed specification |
| ---------- | ----------- | -------------------- |
| *rms*      | Noise model is specified in terms of an RMS map per channel, $N_{pp'}^{ss'}=(\sigma^s_p)^2 \delta_{pp'}\delta_{ss'}$, where $\sigma^s_p$ is given as a standard HEALPix map file. The unit must match the channel definition. | [HEALPix file format specification](https://healpix.sourceforge.io/data/examples/healpix_fits_specs.pdf) |
| *QU_cov*   |  Noise model is specified in a dense $2N_{\mathrm{pix}}\times 2N_{\mathrm{pix}}$ noise covariance matrix including Stokes $Q$ and $U$ parameters. Note: To save initialization time, this covariance matrix will be pre-processed (masked and inverted) the first time the code is called, and stored in `{data_directory}/{filename}_precomp'. If this file already exists, then that information will be used instead. Warning: Remember to delete the precomputed file when changing processing mask! | [Low-resolution WMAP noise covariance matrix](https://lambda.gsfc.nasa.gov/product/map/dr5/ninv_info.cfm) |

### Beams and pixel windows

Beam convolution represents the effect of spatial smoothing due to the finite angular resolution of an instrument. Pixel window smoothing is the analogous effect arising from binning discrete observations into pixelized maps. The following types of beam and pixel window specifications are used in Commander:

| Beam type | Description |  &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Detailed specification &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;|
| ---------- | ----------- | -------------------- |
| *b_l*      | Beam model specified in terms an azimuthally symmetric transfer function, $b_{\ell}$. Used for full-sky convolution of diffuse components. Uses the HEALPix format of ASCII FITS table | [HEALPix file format specification](https://healpix.sourceforge.io/data/examples/healpix_fits_specs.pdf) |
| *febecop*   | Effective beam per pixel as evaluated by FeBeCOP; used for point source convolution. HDF file containing one group per source present in COMP_CATALOGxx. Each group is denoted by the pixel number of the corresponding source center. Pixel numbers must match the HEALPix resolution (Nside) to which the file is applied | HDF file structure:<br>`[pixel number]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# pixel number of point source location`<br>&nbsp; &nbsp; &nbsp;`indices[n] (int)`&nbsp; &nbsp; &nbsp; `# FeBeCOP pixel centers of sub-beam`<br>&nbsp; &nbsp; &nbsp;`   values[n]  (dp)`&nbsp; &nbsp; &nbsp; `# FeBeCOP effective beam values`<br> where $n$ indicates the number of sub-pixels in the effective beam. One such group must be present for each source in the requested catalog |
| *pixwin*      | Legendre transform of pixel smoothing window; typically given by HEALPix pixel windows (e.g., `pixel_window_n0128.fits` available in the HEALPix distribution), but may be computed separately for each experiment, taking into account the actual pointing. Uses the HEALPix format of ASCII FITS table | [HEALPix file format specification](https://healpix.sourceforge.io/data/examples/healpix_fits_specs.pdf) |

### Bandpass

To account for the effect of frequency averaging of a real-world detector, the user must provide a specification of the relative detector sensitivity as a function of frequency, $\tau(\nu)$. The BAND_BANDPASS_TYPE parameter determines the units of $\tau$, as well as potential low-amplitude thresholds applied during bandpass integration.

| Filetype | Description |   Detailed specification |
| ---------- | ----------- | -------------------- |
| *ASCII*      | Bandpass specified in terms of an ASCII table with $(\nu, \tau(\nu))$ in each row. Frequency is given in GHz. | File structure:<br>`nu tau(nu)`<br>Comments are marked by \#; blank lines are not allowed |
| *HDF*        | Bandpass specified in instrument HDF file; used if suffix is '.h5' or '.hdf' | see "Instrument file" below |

## Low-level data files

### Instrument file

An HDF instrument file is used to collect different types of
instrument parameters in a convenient format. Since this serves a
similar role as the Planck Reduced Instrument MOdel (RIMO), this file
is referred to as the Commander RIMO. A full specification is provided
in the table below. However, not all elements must be present for all
channels, but only those that actually will be used in a given
analysis.

| Quantity      | Description |   HDF path |
| ----------    | ----------- | -------------------- |
| Bandpass specification      | Detector sensitivity as a function of frequency | `[channel label]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# Channel label, e.g., '030'`<br>&nbsp; &nbsp; &nbsp;`centFreq (dp)`&nbsp; &nbsp; &nbsp; `# Band center frequency in GHz`<br>&nbsp; &nbsp; &nbsp;`bandpassx[n] (dp)`&nbsp; &nbsp; &nbsp; `# Frequency in GHz`<br>&nbsp; &nbsp; &nbsp;`   bandpass[n]  (dp)`&nbsp; &nbsp; &nbsp; `# Bandpass amplitude`<br> where $n$ indicates the number of frequency samples |
| 4pi beam    | Spherical harmonics decomposition of full $4\pi$ beam. Real harmonic coefficients are listed using the Libsharp order. Used primarily for orbital dipole convolution | `[channel label]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# Channel label, e.g., '030'`<br>&nbsp; &nbsp; &nbsp;`beamlmax (int)`&nbsp; &nbsp; &nbsp; `# Maximum multipole`<br>&nbsp; &nbsp; &nbsp;`   beammmax  (int)`&nbsp; &nbsp; &nbsp; `# Maximum azimuthal order`<br>&nbsp; &nbsp; &nbsp;`   beam  (group)`&nbsp; &nbsp; &nbsp; `# Group name`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   T(n)  (dp)`&nbsp; &nbsp; &nbsp; `# Intensity beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   E(n)  (dp)`&nbsp; &nbsp; &nbsp; `# E-mode beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   B(n)  (dp)`&nbsp; &nbsp; &nbsp; `# B-mode beam`<br>where $n$ is the number of spherical harmonics determined by `beamlmax` and `beammmax`|
| Far sidelobe beam    | Spherical harmonics decomposition of far sidelobes. Real harmonic coefficients are listed using the Libsharp order. Note that this usually has a lower resolution than the $4\pi$ beam | `[channel label]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# Channel label, e.g., '030'`<br>&nbsp; &nbsp; &nbsp;`sllmax (int)`&nbsp; &nbsp; &nbsp; `# Maximum multipole`<br>&nbsp; &nbsp; &nbsp;`   slmmax  (int)`&nbsp; &nbsp; &nbsp; `# Maximum azimuthal order`<br>&nbsp; &nbsp; &nbsp;`   sl  (group)`&nbsp; &nbsp; &nbsp; `# Group name`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   T(n)  (dp)`&nbsp; &nbsp; &nbsp; `# Intensity beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   E(n)  (dp)`&nbsp; &nbsp; &nbsp; `# E-mode beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   B(n)  (dp)`&nbsp; &nbsp; &nbsp; `# B-mode beam`<br>where $n$ is the number of spherical harmonics determined by `sllmax` and `slmmax`|
| Beam summary statistics    | Various beam summary statistics | `[channel label]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# Channel label, e.g., '030'`<br>&nbsp; &nbsp; &nbsp;`fwhm (dp)`&nbsp; &nbsp; &nbsp; `# FWHM in arcmin`<br>&nbsp; &nbsp; &nbsp;`   elip  (dp)`&nbsp; &nbsp; &nbsp; `# Beam ellipticity`<br>&nbsp; &nbsp; &nbsp;`   mbeam_eff  (dp)`&nbsp; &nbsp; &nbsp; `# Main beam efficiency` <br>&nbsp; &nbsp; &nbsp;`   psi_ell  (dp)`&nbsp; &nbsp; &nbsp; `# Ellipticity rotation angle in degrees`|

### TOD datafile

The TOD data files contain the actual time-ordered information for the detectors. These files 
will contain the majority of the data volume of the entire project. The files are in HDF5 format
and use multiple levels of indexing. The top level references the unique chunk index, and the second
level refers to the detector names. Chunks can have arbitrary length, but should be chosen such 
that the following conditions are true: 1) the noise is stationary over the entire chunk, and 2)
the data volume of a single chunk is in the 1-100 MB range. This seems to be the size that gives
the best read performance on our systems.

| Quantity | Description | Location |
| -------- | ------------| ---------|
| **Common Parameters** | | |
| Detector Labels | A comma separated string of all the detector names that appear in this file | /common/det |
| Sampling Frequency | The sampling frequency of these detectors in Hz | /common/fsamp |
| Main beam angle | The main beam offset polarization angles | /common/mbang |
| Number of psi bins used in compression | This should be populated automatically by the file generation code. For performance reasons we suggest something in the form 2^n | /common/npsi |
| nside of the pointing | The nside at which the pointing is compressed. | /common/nside |
| pid list in this file | The unique pointing IDs (chunks) that are present in this file | /common/pid |
| polarization angles | The polarization angles of the detectors. | /common/polang |
| version information | A unique version identifier so you can tell your files apart | /common/version |
| **Per-chunk parameters**| | |
| Huffman symbol list | This should be generated automatically | /[chunk_num]/common/huffsymb |
| Huffman tree array | This should be generated automatically | /[chunk_num]/common/hufftree |
| Load balancing parameters | An array of shape (2,1) that contains some estimate of this chunk's proximity to other chunks. Chunks with similar parameters will be given to the same core. | [chunk_num]/common/load |
| Number of samples | The length of the tod, flag and pointing arrays in this chunk. | /[chunk_num]/common/ntod |
| Satellite position | The x,y,z positions of the satellite as a function of time in solr system coordinates. Ground based experiments should use the position of the earth. | /[chunk_num]/common/satpos |
| Time | The time at the start of this chunk. Space is given for 3 different units if desired. | /[chunk_num]/common/time |
| Satellite velocity | The x,y,z velocity of the satellite relative to the sun. Used to calculate the orbital dipole. | /[chunk_num]/common/vsun | 
| **Per-detector parameters** | | |
| Flags | The huffman compressed flag field | /[chunk_num]/[detector_name]/flag |
| Pointing | The healpy and huffman compressed pointing information | /[chunk_num]/[detector_name]/pix |
| Psi | The huffman compressed polarization angle | /[chunk_num]/[detector_name]/psi |
| Other scalar information | A length 4 vector that contains default values for (gain, sigma0, fknee, alpha). | /[chunk_num]/[detector_name]/scalars |
| Data | The actual measurements of the sky | /[chunk_num]/[detector_name]/tod |

### TOD filelist
The filelist is a text file which contains a list of all the data chunks that are to be analyzed. 

| Column | Quantity | Descrition |
| ------ | -------- | ---------- |
| Column 1, row 1 | File length | The number of entries in the file. For the full dataset, it should be the number of lines in the file-1. If you only want some of the data, you can lower this number to read only the first n lines. If it is longer than the number of entries, the code will crash. |
| 1 | Chunk ID | A unique list of all the chunk IDs that appear in the hdf file. These do not have to be in order. |
| 2 | Path | A full path to the file that contains that chunk |
| 3 | Time estimate | A time estimate about how long this chunk would take to process. This field is output from Commander during a run, re-writing this file. To build the file in the first place, just use 1 here |
| 4 + 5 | Load balancing parameters | The same load balancing numbers that appear in the HDF files. Chunks with similar numbers here get grouped togeather. Can just used 0, 0 to ignore this optimization. |

## Model-specific files

### Map-level instrument file

INSTRUMENT_PARAM_FILE defines the default absolute calibration and bandpass correction for each frequency channel. The format is
```
[channel]     [gain]      [bandpass params]  
```
Comments are allowed, but not blank lines.

Example:
```
$ cat data/instrument_params_BP7.3.dat
#       Band          Gain         delta_bp(0)%p
            030         1      0.30
            044         1      0.00
            070         1      0.00
            857         1      0.00
     0.4-Haslam         1      0.00
    030-WMAP_Ka         1      0.00
    040-WMAP_Q1         1      0.00
    040-WMAP_Q2         1      0.00
    060-WMAP_V1         1      0.00
    060-WMAP_V2         1      0.00
            353         1      0.00
  033-WMAP_Ka_P         1      0.00
   041-WMAP_Q_P         1      0.00
   061-WMAP_V_P         1      0.00
```

### Angular power spectrum

For diffuse sky components with `COMP_CL_TYPExx = binned` or `single_l`, an input angular power spectrum must be defined through the `COMP_CL_DEFAULT_FILExx` parameter. This file must have the same format as outputted by [CAMB](http://camb.info), i.e., each line must contain
```
ell    D_l_TT  D_l_TE  D_l_EE  D_l_BB      
```
Additional columns are allowed, but not used. Unlike CAMB, however, entries for the monopole and dipole must exist, but may be set to zero.

Example:
```
     #    L    TT             TE             EE             BB             PP
     0   1d3            0              0              0              0
     1   1d3            0              0              0              0
     2   0.101673E+04   0.261753E+01   0.308827E-01   0.181847E-05   0.501352E-07
     3   0.963727E+03   0.293806E+01   0.396903E-01   0.363743E-05   0.609943E-07
     4   0.912608E+03   0.275866E+01   0.344962E-01   0.606345E-05   0.702592E-07
     5   0.874477E+03   0.235185E+01   0.230941E-01   0.909717E-05   0.782921E-07
    ...
```

### Line emission definition

Line emission components must be defined in `COMP_LINE_TEMPLATExx`. Each line must contain the following entries:
```
[channel]   [Line ratio prior mean]   [Line ratio prior RMS]   [line2RJ]    [Poltype]
```
where `[line2RJ]` is the conversion factor between line intensity in K km/s and brightness temperature, uK_RJ, and `[poltype]` specifies the polarization type for the current component (for now, only "temperature=1" is supported).

Example:
```
# Channel   LR mean     LR RMS    line2RJ   Poltype
   100         1.         0.      11.06       1
   217         0.6        0.      14.01       1
   353         0.3        0.      12.24       1

```



### Monopole definition

The `md` (monopole+dipole) diffuse component type must be defined through `COMP_MD_DEFINITION_FILExx`. Each line must contain the following entries:
```
[channel] [monopole default] [X dipole default] [Y dipole default] [Z dipole default] [monopole prior mean] [X dipole prior mean] [Y dipole prior mean] [Z dipole prior mean] [monopole prior RMS] [dipole prior RMS]
```
All values must be given in the same units as the `[channel]` sky map. Note that the dipole is specified in Cartesian coordinates. Also, the dipole prior RMS applies to all three components.

Example:
```
   030        0.00      0.00    0.00   0.00    0.00    0.00    0.00   0.00   1e-5    0.00
   044        0.00      0.00    0.00   0.00    0.00    0.00    0.00   0.00   1e-5    0.00
   070        0.00      0.00    0.00   0.00    0.00    0.00    0.00   0.00   1e-5    0.00
0.4-Haslam   12.9e6     0.89e6  3.2e6  0.70e6  12.9e6  0.89e6  3.2e6  0.70e6 1e-5    1e-5
 030-WMAP_Ka   0.00687   0.00   0.00   0.00  -0.0012   0.00    0.00   0.00   0.001   0.00
```

### Point source catalog

The `ptsrc` (point source) component type must be defined through `COMP_MD_DEFINITION_FILExx`. Each line must contain the following entries:
```
[longitude (deg)] [latitude (deg)] [flux density prior mean] [flux density prior RMS] [alpha prior mean] [beta prior mean] [alpha prior RMS] [beta prior RMS] [chisq] [catalog ID]
```
All flux densities should have the same unit as the defining component (typically mJy). The `[chisq]` entry is not used for input. However, the same file format is used to output point source samples, and this entry then indicates the reduced normalized chisquared for the current source. For power-law SED types, only `alpha` is used.

Example:
```
 # 
 # SED model type      = radio
 # Reference frequency =      30.00 GHz
 # 
 # Glon(deg) Glat(deg)     I(mJy)        I_RMS(mJy)  alpha_I beta_I a_RMS_I  b_RMS_I chisq    ID
 184.543  -5.782           344233.000    68846.600  -0.355  0.000   0.000   0.000   0.000   PCCS2
 287.512  -0.654           200787.859    40157.572  -0.434  0.000   0.000   0.000   0.000   PCCS2
 291.538  -0.565           187444.109    37488.822  -0.640  0.000   0.000   0.000   0.000   PCCS2
```

### Point source template library

To save precomputation time, all pixelized beam profiles are stored in `COMP_PTSRC_TEMPLATExx` if `COMP_OUTPUT_PTSRC_TEMPLATExx = .true.`; in this case, the file cannot already exist. On the other hand, if the latter parameter if set to `.false.`, then `COMP_PTSRC_TEMPLATExx` file must exist, as generated from a previous run. Normally, it should never be necessary to inspect `COMP_PTSRC_TEMPLATExx`, as it is automatically generated when requested by the user.

For completeness, however, the template library is organized in a standard HDF format, with each point source template defined by a mini-map (similar to [FeBeCOP](https://wiki.cosmos.esa.int/planckpla/index.php/Effective_Beams#Comparison_of_the_images_of_compact_sources_observed_by_Planck_with_FEBeCoP_products) beam mini-maps) centered on the specified point source location. Each mini-map has two arrays, namely a list of HEALPix pixel indices (HDF path `[channel]/[center pixel]/indices`) and the corresponding point source template value (HDF path `[channel]/[center pixel]/values`).

### SED template

Some diffuse components require an input spectral energy density template, $S(\nu)$. This template must be defined in `COMP_SED_TEMPLATExx`. Each line in this file must contain the following entries:
```
[nu (GHz)] [SED (MJy/sr)]
```
Note that the normalization of this template (and therefore the specific unit) is irrelevant, as the template is anyway re-normalized to unity at the user-specified reference frequency for the given component. Only the shape is important.

Example:
```
 #        nu            SED
     0.050230790   4.6097699e-28
     0.050695569   4.8651838e-28
     0.051164650   5.1346197e-28
     0.051638070   5.4188411e-28
     0.052115871   5.7186464e-28
     0.052598094   6.0348780e-28
     0.053084778   6.3684237e-28
     0.053575965   6.7202225e-28
        ...
```

## MCMC proposal files
   		 
### Bandpass MCMC file

When sampling bandpass parameters for TOD data, the Metropolis covariance matrix and initial point must be defined in `BAND_TOD_BP_INIT_PROPxxx`. The format of this file is as follows:
```
INIT MEAN             [parameter number] [frequency band offset value] # (for mean value)
or
INIT [detector] [parameter number] [detector differential value] # (for individual detector; sum over detectors should be 0)
or
PROP   MEAN            MEAN           [parameter number]   [Covariance]
or
PROP   [detector 1]    [detector 2]   [parameter number]   [Covariance]
```
The `[parameter number]` refers to the complexity of the bandpass model; for now, only single-parameter models are supported, and this should therefore always be 1. For now, this file must be constructed by hand (using intuition and experience to define the proposal covariance matrix), but support for automatic generation from a pre-existing chain file will eventually be implemented in [c3pp](https://www.github.com/cosmoglobe/c3pp). Covariance matrix entries that are not defined are set to zero. Comments (marked with `#`) and blank lines are allowed.

Example:
```
# Initial point for bandpass parameters for each detector
# Number after label refers to bandpass parameter number,
# most often equal to 1
INIT   MEAN 1    -0.10
INIT   27M  1   -0.0178409
INIT   27S  1    0.2038461
INIT   28M  1   -0.1709084
INIT   28S  1   -0.0150967

# Covariance matrix for each parameter; will be decomposed to sqrt
PROP   MEAN MEAN 1   0.001
PROP   27M  27M  1   0.00001
PROP   27S  27S  1   0.00001
PROP   28M  28M  1   0.00001
PROP   28S  28S  1   0.00001
```



### Spectral index proposal file

When sampling spectral parameters over regions for diffuse components, the Metropolis proposal covariance matrix may be specified through `COMP_BETA_ALMSAMP_INITxx`. If this is set to `none`, then an appropriate proposal matrix will be automatically generated during the first iterations by running a longer Monte Carlo chain, and outputted in `OUTPUT_DIRECTORY`.

The format of this file is
```
[nstep T] [nstep Q] [nstep U]
# Comment
[Covariance matrix T]
# Comment
[Covariance matrix Q]
# Comment
[Covariance matrix U]
```
If a common spectral index is sampled for T+Q+U or Q+U, only the first of these need to be defined.

Example:
```
           0          12           0
 Proposal matrix L for signal           1
   0.00000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.10000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.10000   0.00000   0.00000
   0.00000   0.00000   0.00000   0.10000   0.00000
   0.00000   0.00000   0.00000   0.00000   0.10000
 Proposal matrix L for signal           2
   0.00000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.02051   0.00000   0.00303
   0.00000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00303   0.00000   0.02106
 Proposal matrix L for signal           3
   0.00000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000   0.00000   0.00000
```