from astropy.io import fits
"""
top level -- 30, 44, 70 <= frequencies for each of these there is a bandpass
030                      Group
044/                     Group
    bandpass                 Dataset {319}
    bandpassx                Dataset {319}
070                      Group
18M                      Group
18S                      Group
19M                      Group
19S                      Group
20M                      Group
20S                      Group
21M                      Group
21S                      Group
22M                      Group
22S                      Group
23M                      Group
23S                      Group
24M                      Group
24S                      Group
25M                      Group
25S                      Group
26M                      Group
26S                      Group
27M                      Group
27S                      Group
28M/                     Group
    bandpass                 Dataset {310} <= response 
    bandpassx                Dataset {310} <= independent variable (frequency), array of all the freqs bandpass was measured at
    beam/                    Group <= a_lm representation of the beams
        B                        Dataset {9006001}
        E                        Dataset {9006001}
        T                        Dataset {9006001}
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
28S                      Group
common/                  Group
    version

Q
W
<det1>
<det2>
...
common/
       version
"""
filepath = "/mn/stornext/u3/hke/quiet_data/auxilliary"
data = fits.open(f'{filepath}/quiet_qband_temp_beam_140910.fits')
print(data[0].header)
print(data[1].header)
