# MCMC proposal files
   		 
## Bandpass MCMC file

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

## Spectral index proposal file

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
