title: File Formats

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

### Noise specification

The statistical properties of the instrumental noise present in the sky maps must be characterized in terms of its covariance structure; its mean is always assumed to be zero. Commander currently supports two different types of noise specifications:

| Noise type | Description | Detailed specification |
| ---------- | ----------- | -------------------- |
| *rms*      | Noise model is specified in terms of an RMS map per channel, $N_{pp'}^{ss'}=(\sigma^s_p)^2 \delta_{pp'}\delta_{ss'}$, where $\sigma^s_p$ is given in the HEALPix map file. | [HEALPix file format specification](https://healpix.sourceforge.io/data/examples/healpix_fits_specs.pdf) |
| *QU_cov*   |  Noise model is specified in a dense $2N_{\mathrm{pix}}\times 2N_{\mathrm{pix}}$ noise covariance matrix including Stokes $Q$ and $U$ parameters. | [Low-resolution WMAP noise covariance matrix](https://lambda.gsfc.nasa.gov/product/map/dr5/ninv_info.cfm) |

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
| $4pi$ beam    | Spherical harmonics decomposition of full $4\pi$ beam. Real harmonic coefficients are listed using the Libsharp order. Used primarily for orbital dipole convolution | `[channel label]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# Channel label, e.g., '030'`<br>&nbsp; &nbsp; &nbsp;`beamlmax (int)`&nbsp; &nbsp; &nbsp; `# Maximum multipole`<br>&nbsp; &nbsp; &nbsp;`   beammmax  (int)`&nbsp; &nbsp; &nbsp; `# Maximum azimuthal order`<br>&nbsp; &nbsp; &nbsp;`   beam  (group)`&nbsp; &nbsp; &nbsp; `# Group name`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   T(n)  (dp)`&nbsp; &nbsp; &nbsp; `# Intensity beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   E(n)  (dp)`&nbsp; &nbsp; &nbsp; `# E-mode beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   B(n)  (dp)`&nbsp; &nbsp; &nbsp; `# B-mode beam`<br>where $n$ is the number of spherical harmonics determined by `beamlmax` and `beammmax`|
| Far sidelobe beam    | Spherical harmonics decomposition of far sidelobes. Real harmonic coefficients are listed using the Libsharp order. Note that this usually has a lower resolution than the $4\pi$ beam | `[channel label]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# Channel label, e.g., '030'`<br>&nbsp; &nbsp; &nbsp;`sllmax (int)`&nbsp; &nbsp; &nbsp; `# Maximum multipole`<br>&nbsp; &nbsp; &nbsp;`   slmmax  (int)`&nbsp; &nbsp; &nbsp; `# Maximum azimuthal order`<br>&nbsp; &nbsp; &nbsp;`   sl  (group)`&nbsp; &nbsp; &nbsp; `# Group name`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   T(n)  (dp)`&nbsp; &nbsp; &nbsp; `# Intensity beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   E(n)  (dp)`&nbsp; &nbsp; &nbsp; `# E-mode beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   B(n)  (dp)`&nbsp; &nbsp; &nbsp; `# B-mode beam`<br>where $n$ is the number of spherical harmonics determined by `sllmax` and `slmmax`|
| Beam summary statistics    | Various beam summary statistics | `[channel label]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# Channel label, e.g., '030'`<br>&nbsp; &nbsp; &nbsp;`fwhm (dp)`&nbsp; &nbsp; &nbsp; `# FWHM in arcmin`<br>&nbsp; &nbsp; &nbsp;`   elip  (dp)`&nbsp; &nbsp; &nbsp; `# Beam ellipticity`<br>&nbsp; &nbsp; &nbsp;`   mbeam_eff  (dp)`&nbsp; &nbsp; &nbsp; `# Main beam efficiency` <br>&nbsp; &nbsp; &nbsp;`   psi_ell  (dp)`&nbsp; &nbsp; &nbsp; `# Ellipticity rotation angle in degrees`|

### TOD datafile

### TOD filelist


## Model-specific files

### Map-level instrument file

### Angular power spectrum

### Line emission definition

### Monopole definition

### Point source catalog

### Point source library

### SED template

### Template definition

## MCMC proposal files

### Bandpass MCMC file

### Spectral index MCMC file
