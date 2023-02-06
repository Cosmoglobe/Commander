# Model-specific files

## Map-level instrument file

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

## Angular power spectrum

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

## Line emission definition

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

## Monopole definition

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

## Point source catalog

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

## Point source template library

To save precomputation time, all pixelized beam profiles are stored in `COMP_PTSRC_TEMPLATExx` if `COMP_OUTPUT_PTSRC_TEMPLATExx = .true.`; in this case, the file cannot already exist. On the other hand, if the latter parameter if set to `.false.`, then `COMP_PTSRC_TEMPLATExx` file must exist, as generated from a previous run. Normally, it should never be necessary to inspect `COMP_PTSRC_TEMPLATExx`, as it is automatically generated when requested by the user.

For completeness, however, the template library is organized in a standard HDF format, with each point source template defined by a mini-map (similar to [FeBeCOP](https://wiki.cosmos.esa.int/planckpla/index.php/Effective_Beams#Comparison_of_the_images_of_compact_sources_observed_by_Planck_with_FEBeCoP_products) beam mini-maps) centered on the specified point source location. Each mini-map has two arrays, namely a list of HEALPix pixel indices (HDF path `[channel]/[center pixel]/indices`) and the corresponding point source template value (HDF path `[channel]/[center pixel]/values`).

## SED template

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
