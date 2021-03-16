title: Model

# Model Parameters

Commander supports analysis of observations with different angular
resolution and pixelization, measured with different instruments,
adopting a model that also allows different angular resolution per
physical component. To support such flexibility, each data set (or
"frequency band") must be uniquely and completely specified in the
parameter file. In addition to one such definition per data set, the
user needs to specify a number of global parameters that apply to all
data sets.

The Commander model supports a wide range of physical components, each
with a complete individual specification in terms of spectral energy
density, spatial and spectral priors, angular resolution etc. As in
the case of data sets, this is controlled with a combination of global
and component specific parameters.

At present, three different classes of components are supported:
- *diffuse*: Spatially continuous components; modelled in terms of spherical harmonics coefficients, $a_{\ell}$, with an optional power spectrum prior. Examples include CMB, synchrotron, thermal dust mission etc.
- *ptsrc*: Unresolved components; modelled in terms of catalog of discrete source positions, each with a flux density and spectral parameters. Examples include radio, FIR or SZ sources.
- *template*: Template components; modelled in terms of a fixed spatial template with a single free overall amplitude and a given spectral energy density. Examples include Zodical Light Emission, relativistic CMB quadrupole correction etc.

Currently available components are:

| Name                                   | Type         | Class    | SED parameters         | Parametric SED        | Spatial prior |
| ----                                   | ----------   | -----    | --------------         | :-------------------: | :-----------: |
| CMB                                    | cmb          | diffuse  | none                   |                       |               |
| Far-infrared point sources             | fir          | ptsrc    | T, beta                |                       |               |
| Fixed template                         | template     | template | none                   |                       |               |
| Free-free                              | freefree     | diffuse  | T_e                    |                       |               |
| Modified black-body                    | MBB          | diffuse  | T, beta                |                       |               |
| Monopole+dipole                        | md           | diffuse  | None                   |                       |               |
| Line emission                          | line         | diffuse  | line ratio per channel |                       |               |
| Power-law                              | power_law    | diffuse  | beta                   |                       |               |
| Radio point sources                    | radio        | ptsrc    | Alpha, beta            |                       |               |
| Relativistic kinematic CMB quadrupole  | cmb_relquad  | template | None                   |                       |               |
| Spinning dust                          | spindust     | diffuse  | nu_p                   |                      |               |
| Tilted spinning dust                   | spindust2    | diffuse  | nu_p                   |                       |               |
| Thermal Sunyaev-Zeldovich sources      | sz           | ptsrc    | none                   |                       |               |


## Global Parameters

| Parameter name           | Description                                                                                                                                                                                           | Accepted values       | Example value |
| ------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------- | ------------- |
| INSTRUMENT_PARAM_FILE    | Initialization file for map-level instrumental parameters (absolute calibration and bandpass corrections) | valid filename in DATA_DIRECTORY | instrument_params_bp_v1.dat |
| INIT_INSTRUMENT_FROM_HDF | Initialization mode; if a previous chain file is specified, it must follow the same convention as INIT_CHAIN         | none={use values from INSTRUMENT_PARAM_FILE}<br>default={use default chain specified in INIT_CHAIN}<br>*previous chain file*={use specified sample} | "chain_init_v2.h5:5" |
| NUM_SIGNAL_COMPONENTS    | Maximum number of signal components; additional components may be present, but these are ignored | >= 1 | 12 |
| T_CMB                          | Mean CMB temperature in Kelvin         | > 0      |  2.7255  |
| NUM_CG_SAMPLING_GROUPS   | Number of (user defined) Conjugate Gradient sampling groups. | 0, no (user defined) amplitude sampling <br> >0, sample amplitudes of components specified in the sampling groups | 1 |
| CG_SAMPLING_GROUPxx        | Conjugate Gradient group, ID = xx.<br> Components may be sampled in separate groups to aid CG convergence, or together in the same group. Components will be estimated, conditionally on components not present in the group, through Gibbs sampling.<br> Any active component without a group will be fixed and conditioned upon | Any component label | 'cmb', 'dust' |
| LOCALSAMP_BURN_IN        | Tune local sampling parameters first N Gibbs iterations | >=0      | 3     |

## Common Parameters

The **XX** suffix in the parameter names below is a placeholder for an
incremental integer that uniquely identifies each component.

The following parameters apply to all component types, and must be present for all classes:

| Parameter name                | Description                                               | Accepted values    |  Example value |
| ----------------------------- | --------------------------------------------------------- | ------------------------ | :---------: |
| COMP_CLASSxx                  | Component class                | Valid component class | diffuse<br>ptsrc<br>template | diffuse |
| COMP_INIT_FROM_HDFxx          | Initialization mode; if a previous chain file is specified, it must follow the same convention as INIT_CHAIN         | none={use default values}<br>default={use default chain specified in INIT_CHAIN}<br>*previous chain file*={use specified sample} | "chain_init_v2.h5:5" |
| COMP_LABELxx                  | User-specified component label | String shorter than 16 characters | cmb |
| COMP_TYPExx                   | Component type                 | Valid component type; see above table        | cmb |
| INCLUDE_COMPxx           | Enable component in current analysis | .true.<br>.false. | .true.              |
| LOCALSAMP_BURN_IN        | Number of iterations that the local sampler will tune the sampling parameters (proposal length and number of proposals) | >= 0 | 3 |
## Diffuse Components

All diffuse components (except the monopole+dipole component, 'md') are described through the following general parameters in addition to the component-specific parameters listed for each type below. For the 'md' component, a complete specification is provided separately:

| Parameter name                | Description                                               | Accepted values    |  Example value |
| ----------------------------- | --------------------------------------------------------- | ------------------------ | :---------: |
| COMP_APPLY_POSITIVITY_PRIORxx | Not yet supported                                         | .true.<br>.false.  | .true. |
| COMP_CG_SCALExx               | Multiplicative amplitude scale used in CG search; useful to ensure that all components have comparable magnitude during convergence criterion assessment | > 0  | 1. |
| COMP_CG_SAMP_GROUP_MAXITERxx  | Maximum iterations/number of iterations for the Conjugate Gradient search after marginal sampling component spectral parameters.<br> Number of iterations to compute if "fixed_iter" is used as convergence criterium. <br> Not appliccable for components without spectral indices, e.g. CMB  | >= 0 | 150 |
| COMP_CL_BETA_PRIOR_MEANxx     | Not yet supported                                         | NA | NA |
| COMP_CL_BETA_PRIOR_RMSxx      | Not yet supported                                         | NA | NA |
| COMP_CL_BIN_FILExx            | Used only when COMP_CL_TYPEXX is set tp binned            |                                                                                        |
| COMP_CL_DEFAULT_AMP_Txx       | Initial $TT$ power spectrum amplitude, $D_0$, in the same (squared) physical units at the amplitude map. Set to a large value to impose a weak spatial prior on amplitudes. Only applies to 'power_law', 'exp' and 'gauss' types        | > 0 | 1e6 |                                                                                        |
| COMP_CL_DEFAULT_AMP_Exx       | Same as COMP_CL_DEFAULT_AMP_Txx, but for $EE$ spectrum | > 0 | 1e4 |
| COMP_CL_DEFAULT_AMP_Bxx       | Same as COMP_CL_DEFAULT_AMP_Txx, but for $BB$ spectrum | > 0 | 1e4 |
| COMP_CL_DEFAULT_BETA_Txx      | Initial $TT$ power spectrum shape parameter, $\beta$ | any float | 0. |
| COMP_CL_DEFAULT_BETA_Exx      | Initial $EE$ power spectrum shape parameter, $\beta$ | any float | 0. |
| COMP_CL_DEFAULT_BETA_Bxx      | Initial $BB$ power spectrum shape parameter, $\beta$ | any float | 0. |
| COMP_CL_DEFAULT_FILExx        | Initial power spectrum; used only for COMP_CL_TYPE='binned' | valid filename in DATA_DIRECTORY | 'planck_bf_lcdm_v1.dat' |
| COMP_CL_L_PIVOTxx             | Pivot multipole, $\ell_0$, for angular power spectrum models | >= 0 | 100 |
| COMP_CL_POLTYPExx             | Not active yet  | NA | NA |
| COMP_CL_TYPExx                | Angular power spectrum type; typically use 'binned' for CMB analysis and 'power_law' for foreground analysis. Optionally, 'gauss' is useful for foregrounds with particularly bright compact features, to avoid ringing around compact sources. Note that these power spectrum parameters are *priors* on the amplitudes, not post-analysis smoothing operators; as such, the derived amplitudes may have excess power relative to the specified spectrum.     | none = {No spatial prior, $S^{-1}=0$}<br>binned = {bin $D_{\ell}$ according to specified bin file}<br>power_law = {$D_{\ell}=D_0(\ell/\ell_0)^{\beta}$}<br>exp = {$D_{\ell}=D_0 e^{-\beta(\ell/\ell_0)}$}<br>gauss = {$D_{\ell} = D_0 e^{-\frac{1}{2}\ell(\ell+1)\sigma^2}$ with $\sigma = \beta\cdot\pi/180/\sqrt{8\log 2}$}                               | binned |
| COMP_INDMASKxx                | Mask used used for spectral index estimation; pixels with value smaller than 0.5 are ignored during likelihood evaluations, but derived spectral indices may still apply to those pixels if they are estimated non-locally. Not needed for components with no free spectral parameters  | fullsky<br>valid filename in DATA_DIRECTORY | fullsky |
| COMP_INIT_FROM_HDFxx          | Initialization mode; if a previous chain file is specified, it must follow the same convention as INIT_CHAIN         | none={use values from INSTRUMENT_PARAM_FILE}<br>default={use default chain specified in INIT_CHAIN}<br>*previous chain file*={use specified sample} | "chain_init_v2.h5:5" |
| COMP_INPUT_AMP_MAPxx          | Input amplitude map; should have same units and reference frequency as main component. Beware: The output FWHM is used to deconvolve this map prior to $a_{\ell m}$ initialization! | none <br> valid filename in DATA_DIRECTORY| none |
| COMP_APPLY_JEFFREYS_PRIORxx   | Not yet active; only relevant for components with free spectral parameters                            | .true.<br>.false. | .false. |                                                                     |
| COMP_L_APODxx                 | Enables cosine apodization between $\ell_{\mathrm{apod}}$ and $\ell_{\mathrm{max}}$; useful to suppress ringing around bright point sources  | >= 0 | 1500 |
| COMP_LMAX_AMPxx               | Maximum multipole moment for amplitude map | >= 0 | 1500 |
| COMP_LMAX_INDxx               | Maximum multipole moment for spectral index map <br> (Is to be phased out by COMP_{parameter}_{pol}_LMAXxx to allow for more diverse sampling) | -1 = {Use pixel-by-pixel index map}<br>>=0 = {$\beta = \sum_{\ell=0}^{\ell_{\mathrm{max}}^{\mathrm{ind}}}\sum_{m=-\ell}^{\ell} a_{\ell m} Y_{\ell m}$} | 30                                     |                                                                                        |
| COMP_MASKxx                   | Component mask; pixels with value less than 0.5 are set to zero   | fullsky<br>valid filename in DATA_DIRECTORY | fullsky |
| COMP_NSIDExx                  | HEALPix resolution parameter for output and spectral index maps  | > 0 | 512 |
| COMP_NU_REF_Txx                 | SED reference frequency in GHz; applies only to intensity component, $T$      | > 0 | 30. |
| COMP_NU_REF_Pxx                 | SED reference frequency in GHz; applies only to polarization parameters, $Q$ and $U$      | > 0 | 30. |
| COMP_OUTPUT_EB_MAPXX          | Output $E$ and $B$ maps in addition to $Q$ and $U$                  | .true.<br>.false. | .false. |
| COMP_OUTPUT_FWHMxx            | Gaussian smoothing scale in arcmin applied to output FITS maps; only for visualization | >= 0 | 60. |
| COMP_POLARIZATIONxx           | Polarization switch                                       | .true.<br>.false.  | .true. |
| COMP_PRIOR_AMP_MAPxx          | Prior mean map; if present, the angular power spectrum prior refers to fluctuations relative to this mean map | none <br> valid filename in DATA_DIRECTORY | none |
| COMP_UNITxx                   | Physical unit of amplitude map | uK_cmb, mK_cmb, K_cmb<br>MJy/sr<br>uK_RJ<br>K km/s | uK_cmb |

### Free-free Component

| Parameter name                | Description                                               | Accepted values        | Example value |
| ----------------------------- | --------------------------------------------------------- | ---------------------- | ------------- |
| COMP_INPUT_T_E_MAPxx         | Initial map for the electorn temperature $T_e$             | 'none' = {Default value used} <br> valid filename in DATA_DIRECTORY | 'none' |
| COMP_T_E_POLTYPExx           | Spectral index Stokes parameter coherence switch | 1 = {fit common parameter for T+Q+U}<br>2 = {fit separate parameters for T and Q+U}<br>3 = {fit separate parameters for T, Q and U} | 2 |
| COMP_T_E_SMOOTHING_SCALExx   | Smoothing scale ID for $T_e$; specifies which SMOOTHING_SCALE parameters to use, as defined in Data Set parameters | <= NUM_SMOOTHING_SCALES | 1 |
| COMP_DEFAULT_T_Exx           | Initial value of $T_e$; must lie within uniform prior values  | float | 7000 |
| COMP_PRIOR_GAUSS_T_E_MEANxx  | Mean of Gaussian prior on $T_e$                         | float                     | 7000   |
| COMP_PRIOR_GAUSS_T_E_RMSxx   | Standard deviation (RMS) of Gaussian prior on $T_e$     | <0 = {Do not apply Gaussian prior}<br>0 = {Fix $T_e$ on input}<br> >0 = {Sample with specifed prior RMS} | 500 |
| COMP_T_E_NU_MINxx            | Lowest frequency for estimating $T_e$ in GHz             | > 0  | 1. |  
| COMP_T_E_NU_MAXxx            | Highest frequency for estimating $T_e$ in GHz            | > 0  | 1000. |
| COMP_PRIOR_UNI_T_E_LOWxx     | Absolute lower bound on $T_e$; hard uniform prior      | float                      | 4000   |
| COMP_PRIOR_UNI_T_E_HIGHxx    | Absolute upper bound on $T_e$; hard uniform prior       | > COMP_PRIOR_UNI_T_E_LOW | 12000   |
| COMP_T_E_{pol}_LMAXxx        | Maximum multipole moment of the $T_e$ spectral index map for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| <0 : {Use local pix-by-pix sampling} <br> >=0 : {Use $a_{lm}$ sampling}| -1 <br> 0 |
| COMP_T_E_{pol}_PIXREGxx      | Type of pixel regions for the local sampling of the $T_e$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} | 'pixreg' = pixel region input map <br> 'fullsky' = one fullsky pixel region <br> 'single_pix' = pixel regions as per smoothing scale $N_{\mathrm{side}}$ | 'pixreg' |
| COMP_T_E_{pol}_LNLTYPExx     | Type of log-likelihood evaluation for sampling of the $T_e$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: using ridge og marginal sampling requires separate CG sampling of the component amplitude afterwards. | 'chisq' <br> 'ridge' <br> 'marginal' | 'chisq' |
| COMP_T_E_PIXREG_MAPxx        | Map of defined pixel regions for $T_e$ for the local sampler to use. Pixel region index values starts at 1 | 'fullsky' = {1 pixel region on the whole sky} <br> 'none' = {whole sky set to prior} <br> valid filename in DATA_DIRECTORY| 'fullsky'|
| COMP_T_E_PIXREG_INITVALUE_MAPxx | Map of $T_e$ where each pixelgerion has a uniform value (no smoothing). Used to initialize pixelregions with exact values if the normal spectral index map is smoothed | 'none' = {use original input map to decide pixelregion values} <br> valid filename in DATA_DIRECTORY| 'none'|
| COMP_T_E_{pol}_NUM_PIXREGxx  | Number of pixel regions to sample, all additional regions are set to prior values (pixel region 0). Only needed for the 'pixreg' pixel region type option for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the proposal length map | > 0 | 10 |
| COMP_T_E_MASKxx              | Mask for use of the local sampler to sample $T_e$                       | 'fullsky' <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_T_E_PROPLENxx           | $T_e$ sampling. Map of the standard deviation of each of the Metropolis-Hastings proposals, per pixel region, for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_T_E_{pol}_PROPLEN_INITxx | (Map of) The initial value of the proposal length of the local sampler for all pixel regions ($T_e$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0.0 <br> <= 0.0 {disable} | 0.05d0 <br> -1.d0 |
| COMP_T_E_{pol}_SAMPLE_PROPLENxx | Flag for whether or not to tune the proposal length of the local sampler for the sampling of the $T_e$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_T_E_NPROPxx             | $T_e$ sampling. Map of the number of Metropolis-Hastings proposals per pixel region for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_T_E_{pol}_NPROP_INITxx  | (Map of) The initial value of the number of proposals of the local sampler for all pixel regions ($T_e$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0 <br> <= 0 {disable} | 100 <br> -1 |
| COMP_T_E_{pol}_SAMPLE_NPROPxx | Flag for whether or not to tune number of proposals of the local sampler for the sampling of the $T_e$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_T_E_UNI_NPROP_LOWxx     | Absolute lower limit of number of MH proposals per iteration. Covers all Stokes parameters | >= 0 | 1    |
| COMP_T_E_UNI_NPROP_HIGHxx    | Absolute upper limit of number of MH proposals per iteration. Covers all Stokes parameters | >= 0 | 1000 |


### Line emission component

| Parameter name                | Description                                               | Accepted values        | Example value |
| ----------------------------- | --------------------------------------------------------- | ---------------------- | ------------- |
| COMP_BAND_REFxx               | Component reference band ID, where the ID must match an active data set label | string shorter than 16 characters | dame |
| COMP_LINE_TEMPLATEXX          | Line component definition file; contains information regarding each data set in which the component is observed | valid filename in DATA_DIRECTORY | co10_line_template_v1.dat |                                                 |                                                                                        |

### Modified Blackbody Component

| Parameter name                | Description                                               | Accepted values        | Example value |
| ----------------------------- | --------------------------------------------------------- | ---------------------- | ------------- |
| COMP_INPUT_BETA_MAPxx         | Initial map for the MBB $\beta$ spectral parameter        | 'none' = {Default value used} <br> valid filename in DATA_DIRECTORY | 'none' |
| COMP_BETA_NU_MINxx            | Lowest frequency for estimating $\beta$ in GHz             | > 0  | 1. |  
| COMP_BETA_NU_MAXxx            | Highest frequency for estimating $\beta$ in GHz            | > 0  | 1000. |
| COMP_BETA_POLTYPExx           | Spectral index Stokes parameter coherence switch | 1 = {fit common parameter for T+Q+U}<br>2 = {fit separate parameters for T and Q+U}<br>3 = {fit separate parameters for T, Q and U} | 2 |
| COMP_BETA_SMOOTHING_SCALExx   | Smoothing scale ID for $\beta$; specifies which SMOOTHING_SCALE parameters to use, as defined in Data Set parameters | <= NUM_SMOOTHING_SCALES | 1 || COMP_DEFAULT_BETAxx           | Initial value of $\beta$; must lie within uniform prior values  | float | 1.5 |
| COMP_DEFAULT_Txx           | Initial value of $T$; must lie within uniform prior values  | float | 21. |
| COMP_PRIOR_GAUSS_BETA_MEANxx  | Mean of Gaussian prior on $\beta$                         | float                     | 1.5   |
| COMP_PRIOR_GAUSS_BETA_RMSxx   | Standard deviation (RMS) of Gaussian prior on $\beta$     | <0 = {Do not apply Gaussian prior}<br>0 = {Fix $\beta$ on input}<br> >0 = {Sample with specifed prior RMS} | 0.2 |
| COMP_BETA_{pol}_LMAXxx        | Maximum multipole moment of the $\beta$ spectral index map for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| <0 : {Use local pix-by-pix sampling} <br> >=0 : {Use $a_{lm}$ sampling}| -1 <br> 0 |
| COMP_BETA_{pol}_PIXREGxx      | Type of pixel regions for the local sampling of the $\beta$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} | 'pixreg' = pixel region input map <br> 'fullsky' = one fullsky pixel region <br> 'single_pix' = pixel regions as per smoothing scale $N_{\mathrm{side}}$ | 'pixreg' |
| COMP_BETA_{pol}_LNLTYPExx     | Type of log-likelihood evaluation for sampling of the $\beta$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: using ridge og marginal sampling requires separate CG sampling of the component amplitude afterwards. | 'chisq' <br> 'ridge' <br> 'marginal' | 'chisq' |
| COMP_BETA_PIXREG_MAPxx        | Map of defined pixel regions for $\beta$ for the local sampler to use. Pixel region index values starts at 1 | 'fullsky' = {1 pixel region on the whole sky} <br> 'none' = {whole sky set to prior} <br> valid filename in DATA_DIRECTORY| 'fullsky'|
| COMP_BETA_PIXREG_INITVALUE_MAPxx | Map of $\beta$ where each pixelgerion has a uniform value (no smoothing). Used to initialize pixelregions with exact values if the normal spectral index map is smoothed | 'none' = {use original input map to decide pixelregion values} <br> valid filename in DATA_DIRECTORY| 'none'|
| COMP_BETA_{pol}_NUM_PIXREGxx  | Number of pixel regions to sample, all additional regions are set to prior values (pixel region 0). Only needed for the 'pixreg' pixel region type option for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the proposal length map | > 0 | 10 |
| COMP_BETA_MASKxx              | Mask for use of the local sampler to sample $\beta$                       | 'fullsky' <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_BETA_PROPLENxx           | $\beta$ sampling. Map of the standard deviation of each of the Metropolis-Hastings proposals, per pixel region, for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_BETA_{pol}_PROPLEN_INITxx | (Map of) The initial value of the proposal length of the local sampler for all pixel regions ($\beta$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0.0 <br> <= 0.0 {disable} | 0.05d0 <br> -1.d0 |
| COMP_BETA_{pol}_SAMPLE_PROPLENxx | Flag for whether or not to tune the proposal length of the local sampler for the sampling of the $\beta$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_BETA_NPROPxx             | $\beta$ sampling. Map of the number of Metropolis-Hastings proposals per pixel region for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_BETA_{pol}_NPROP_INITxx  | (Map of) The initial value of the number of proposals of the local sampler for all pixel regions ($\beta$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0 <br> <= 0 {disable} | 100 <br> -1 |
| COMP_BETA_{pol}_SAMPLE_NPROPxx | Flag for whether or not to tune number of proposals of the local sampler for the sampling of the $\beta$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_BETA_UNI_NPROP_LOWxx     | Absolute lower limit of number of MH proposals per iteration for $\beta$. Covers all Stokes parameters | >= 0 | 1    |
| COMP_BETA_UNI_NPROP_HIGHxx    | Absolute upper limit of number of MH proposals per iteration for $\beta$. Covers all Stokes parameters | >= 0 | 1000 |
| COMP_INPUT_T_MAPxx         | Initial map for the MBB temperature $T$ spectral parameter        | 'none' = {Default value used} <br> valid filename in DATA_DIRECTORY | 'none' |
| COMP_PRIOR_GAUSS_T_MEANxx  | Mean of Gaussian prior on $T$                         | float                     | 21.   |
| COMP_PRIOR_GAUSS_T_RMSxx   | Standard deviation (RMS) of Gaussian prior on $T$     | <0 = {Do not apply Gaussian prior}<br>0 = {Fix $T$ on input}<br> >0 = {Sample with specifed prior RMS} | 3. |
| COMP_PRIOR_UNI_BETA_LOWxx     | Absolute lower bound on $\beta$; hard uniform prior      | float                      | 1.2   |
| COMP_PRIOR_UNI_BETA_HIGHxx    | Absolute upper bound on $\beta$; hard uniform prior       | > COMP_PRIOR_UNI_BETA_LOW | 3.0   |
| COMP_PRIOR_UNI_T_LOWxx     | Absolute lower bound on $T$; hard uniform prior      | float                      | 10.   |
| COMP_PRIOR_UNI_T_HIGHxx    | Absolute upper bound on $T$; hard uniform prior       | > COMP_PRIOR_UNI_T_LOW | 30.   |
| COMP_T_NU_MINxx            | Lowest frequency for estimating $T$ in GHz             | > 0  | 1. |  
| COMP_T_NU_MAXxx            | Highest frequency for estimating $T$ in GHz            | > 0  | 1000. |
| COMP_T_POLTYPExx           | Spectral index Stokes parameter coherence switch | 1 = {fit common parameter for T+Q+U}<br>2 = {fit separate parameters for T and Q+U}<br>3 = {fit separate parameters for T, Q and U} | 2 |
| COMP_T_SMOOTHING_SCALExx   | Smoothing scale ID for $T$; specifies which SMOOTHING_SCALE parameters to use, as defined in Data Set parameters | <= NUM_SMOOTHING_SCALES | 1 |
| COMP_T_{pol}_LMAXxx        | Maximum multipole moment of the $T$ spectral index map for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| <0 : {Use local pix-by-pix sampling} <br> >=0 : {Use $a_{lm}$ sampling}| -1 <br> 0 |
| COMP_T_{pol}_PIXREGxx      | Type of pixel regions for the local sampling of the $T$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence inde\
x 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} | 'pixreg' = pixel region input map <br> 'fullsky' = one fullsky pixel region <br> 'single_pix' = pixel regions as per smoothing scale $N_{\mathrm{side}}$ | 'pixreg' |              
| COMP_T_{pol}_LNLTYPExx     | Type of log-likelihood evaluation for sampling of the $T$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence in\
dex 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: using ridge og marginal sampling requires separate CG sampling of the component amplitude afterwards. | 'chisq' <br> 'ridge' <br> 'marginal' | 'chisq' |                
| COMP_T_PIXREG_MAPxx        | Map of defined pixel regions for $T$ for the local sampler to use. Pixel region index values starts at 1 | 'fullsky' = {1 pixel region on the whole sky} <br> 'none' = {whole sky set to prior} <br> valid filename in DATA_DIRECTORY| 'fullsky'|                                                                                                     
| COMP_T_PIXREG_INITVALUE_MAPxx | Map of $T$ where each pixelgerion has a uniform value (no smoothing). Used to initialize pixelregions with exact values if the normal spectral index map is smoothed | 'none' = {use original input map to decide pixelregion values} <br> valid filename in DATA_DIRECTORY| 'none'|                                                               
| COMP_T_{pol}_NUM_PIXREGxx  | Number of pixel regions to sample, all additional regions are set to prior values (pixel region 0). Only needed for the 'pixreg' pixel region type option for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the proposal length map | > 0 | 10 |                                        
| COMP_T_MASKxx              | Mask for use of the local sampler to sample $T$                       | 'fullsky' <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_T_PROPLENxx           | $T$ sampling. Map of the standard deviation of each of the Metropolis-Hastings proposals, per pixel region, for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|                                                                                                                  
| COMP_T_{pol}_PROPLEN_INITxx | (Map of) The initial value of the proposal length of the local sampler for all pixel regions ($T$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0.0 <br> <= 0.0 {disable} | 0.05d0 <br> -1.d0 |                         
| COMP_T_{pol}_SAMPLE_PROPLENxx | Flag for whether or not to tune the proposal length of the local sampler for the sampling of the $T$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_T_NPROPxx             | $T$ sampling. Map of the number of Metropolis-Hastings proposals per pixel region for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|                
| COMP_T_{pol}_NPROP_INITxx  | (Map of) The initial value of the number of proposals of the local sampler for all pixel regions ($T$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0 <br> <= 0 {disable} | 100 <br> -1 |                                
| COMP_T_{pol}_SAMPLE_NPROPxx | Flag for whether or not to tune number of proposals of the local sampler for the sampling of the $T$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_T_UNI_NPROP_LOWxx     | Absolute lower limit of number of MH proposals per iteration for $T$. Covers all Stokes parameters | >= 0 | 1    |
| COMP_T_UNI_NPROP_HIGHxx    | Absolute upper limit of number of MH proposals per iteration for $T$. Covers all Stokes parameters | >= 0 | 1000 |                                                                                                        


### Monopole and dipole component

| Parameter name                | Description                                               | Accepted values        | Example value |
| ----------------------------- | --------------------------------------------------------- | ---------------------- | ------------- |
| COMP_MD_DEFINITION_FILExx     | Monopole and dipole definition file                       | valid filename in DATA_DIRECTORY | md_bp_v1.dat |


### Power-Law Component

| Parameter name                | Description                                               | Accepted values        | Example value |
| ----------------------------- | --------------------------------------------------------- | ---------------------- | ------------- |
| COMP_INPUT_BETA_MAPxx         | Initial map for the power-law $\beta$ spectral parameter        | 'none' = {Default value used} <br> valid filename in DATA_DIRECTORY | 'none' |
| COMP_BETA_POLTYPExx           | Spectral index Stokes parameter coherence switch | 1 = {fit common parameter for T+Q+U}<br>2 = {fit separate parameters for T and Q+U}<br>3 = {fit separate parameters for T, Q and U} | 2 |
| COMP_BETA_SMOOTHING_SCALExx   | Smoothing scale ID for $\beta$; specifies which SMOOTHING_SCALE parameters to use, as defined in Data Set parameters | <= NUM_SMOOTHING_SCALES | 1 || COMP_DEFAULT_BETAxx           | Initial value of $\beta$; must lie within uniform prior values  | float | -2. |
| COMP_PRIOR_GAUSS_BETA_MEANxx  | Mean of Gaussian prior on $\beta$                         | float                     | -2.   |
| COMP_PRIOR_GAUSS_BETA_RMSxx   | Standard deviation (RMS) of Gaussian prior on $\beta$     | <0 = {Do not apply Gaussian prior}<br>0 = {Fix $\beta$ on input}<br> >0 = {Sample with specifed prior RMS} | 0.2 |
| COMP_BETA_NU_MINxx            | Lowest frequency for estimating $\beta$ in GHz             | > 0  | 1. |  
| COMP_BETA_NU_MAXxx            | Highest frequency for estimating $\beta$ in GHz            | > 0  | 100. |
| COMP_PRIOR_UNI_BETA_LOWxx     | Absolute lower bound on $\beta$; hard uniform prior      | float                      | -3.   |
| COMP_PRIOR_UNI_BETA_HIGHxx    | Absolute upper bound on $\beta$; hard uniform prior       | > COMP_PRIOR_UNI_BETA_LOW | -1.   |
| COMP_BETA_{pol}_LMAXxx        | Maximum multipole moment of the $\beta$ spectral index map for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| <0 : {Use local pix-by-pix sampling} <br> >=0 : {Use $a_{lm}$ sampling}| -1 <br> 0 |
| COMP_BETA_{pol}_PIXREGxx      | Type of pixel regions for the local sampling of the $\beta$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} | 'pixreg' = pixel region input map <br> 'fullsky' = one fullsky pixel region <br> 'single_pix' = pixel regions as per smoothing scale $N_{\mathrm{side}}$ | 'pixreg' |
| COMP_BETA_{pol}_LNLTYPExx     | Type of log-likelihood evaluation for sampling of the $\beta$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: using ridge og marginal sampling requires separate CG sampling of the component amplitude afterwards. | 'chisq' <br> 'ridge' <br> 'marginal' | 'chisq' |
| COMP_BETA_PIXREG_MAPxx        | Map of defined pixel regions for $\beta$ for the local sampler to use. Pixel region index values starts at 1 | 'fullsky' = {1 pixel region on the whole sky} <br> 'none' = {whole sky set to prior} <br> valid filename in DATA_DIRECTORY| 'fullsky'|
| COMP_BETA_PIXREG_INITVALUE_MAPxx | Map of $\beta$ where each pixelgerion has a uniform value (no smoothing). Used to initialize pixelregions with exact values if the normal spectral index map is smoothed | 'none' = {use original input map to decide pixelregion values} <br> valid filename in DATA_DIRECTORY| 'none'|
| COMP_BETA_{pol}_NUM_PIXREGxx  | Number of pixel regions to sample, all additional regions are set to prior values (pixel region 0). Only needed for the 'pixreg' pixel region type option for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the proposal length map | > 0 | 10 |
| COMP_BETA_MASKxx              | Mask for use of the local sampler to sample $\beta$                       | 'fullsky' <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_BETA_PROPLENxx           | $\beta$ sampling. Map of the standard deviation of each of the Metropolis-Hastings proposals, per pixel region, for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_BETA_{pol}_PROPLEN_INITxx | (Map of) The initial value of the proposal length of the local sampler for all pixel regions ($\beta$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0.0 <br> <= 0.0 {disable} | 0.05d0 <br> -1.d0 |
| COMP_BETA_{pol}_SAMPLE_PROPLENxx | Flag for whether or not to tune the proposal length of the local sampler for the sampling of the $\beta$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_BETA_NPROPxx             | $\beta$ sampling. Map of the number of Metropolis-Hastings proposals per pixel region for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_BETA_{pol}_NPROP_INITxx  | (Map of) The initial value of the number of proposals of the local sampler for all pixel regions ($\beta$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0 <br> <= 0 {disable} | 100 <br> -1 |
| COMP_BETA_{pol}_SAMPLE_NPROPxx | Flag for whether or not to tune number of proposals of the local sampler for the sampling of the $\beta$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_BETA_UNI_NPROP_LOWxx     | Absolute lower limit of number of MH proposals per iteration for $\beta$. Covers all Stokes parameters | >= 0 | 1    |
| COMP_BETA_UNI_NPROP_HIGHxx    | Absolute upper limit of number of MH proposals per iteration for $\beta$. Covers all Stokes parameters | >= 0 | 1000 |




### Spinning Dust Component

| Parameter name                | Description                                               | Accepted values        | Example value |
| ----------------------------- | --------------------------------------------------------- | ---------------------- | ------------- |
| COMP_INPUT_NU_P_MAPxx         | Initial map for the Spinning dust peak frequency $\nu_\mathrm{p}$ spectral parameter        | 'none' = {Default value used} <br> valid filename in DATA_DIRECTORY | 'none' |
| COMP_NU_P_NU_MINxx            | Lowest frequency for estimating $\nu_\mathrm{p}$ in GHz             | > 0  | 1. |  
| COMP_NU_P_NU_MAXxx            | Highest frequency for estimating $\nu_\mathrm{p}$ in GHz            | > 0  | 100. |
| COMP_NU_P_POLTYPExx           | Spectral index Stokes parameter coherence switch | 1 = {fit common parameter for T+Q+U}<br>2 = {fit separate parameters for T and Q+U}<br>3 = {fit separate parameters for T, Q and U} | 2 |
| COMP_NU_P_SMOOTHING_SCALExx   | Smoothing scale ID for $\nu_\mathrm{p}$; specifies which SMOOTHING_SCALE parameters to use, as defined in Data Set parameters | <= NUM_SMOOTHING_SCALES | 1 || COMP_DEFAULT_NU_Pxx           | Initial value of $\nu_\mathrm{p}$; must lie within uniform prior values  | float | -2. |
COMP_PRIOR_GAUSS_NU_P_MEANxx  | Mean of Gaussian prior on $\nu_\mathrm{p}$                         | float                     | 21.   |
| COMP_PRIOR_GAUSS_NU_P_RMSxx   | Standard deviation (RMS) of Gaussian prior on $\nu_\mathrm{p}$     | <0 = {Do not apply Gaussian prior}<br>0 = {Fix $\nu_\mathrm{p}$ on input}<br> >0 = {Sample with specifed prior RMS} | 0.2 |
| COMP_PRIOR_UNI_NU_P_LOWxx     | Absolute lower bound on $\nu_\mathrm{p}$; hard uniform prior      | float                      | 10.   |
| COMP_PRIOR_UNI_NU_P_HIGHxx    | Absolute upper bound on $\nu_\mathrm{p}$; hard uniform prior       | > COMP_PRIOR_UNI_NU_P_LOW | 30.   |
| COMP_SED_TEMPLATExx           | SED template file; must be defined in flux density units (MJy/sr) | valid filename in DATA_DIRECTORY | spdust2_cnm.dat |
| COMP_NU_P_{pol}_LMAXxx        | Maximum multipole moment of the $\nu_\mathrm{p}$ spectral index map for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| <0 : {Use local pix-by-pix sampling} <br> >=0 : {Use $a_{lm}$ sampling}| -1 <br> 0 |
| COMP_NU_P_{pol}_PIXREGxx      | Type of pixel regions for the local sampling of the $\nu_\mathrm{p}$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} | 'pixreg' = pixel region input map <br> 'fullsky' = one fullsky pixel region <br> 'single_pix' = pixel regions as per smoothing scale $N_{\mathrm{side}}$ | 'pixreg' |
| COMP_NU_P_{pol}_LNLTYPExx     | Type of log-likelihood evaluation for sampling of the $\nu_\mathrm{p}$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: using ridge og marginal sampling requires separate CG sampling of the component amplitude afterwards. | 'chisq' <br> 'ridge' <br> 'marginal' | 'chisq' |          
| COMP_NU_P_PIXREG_MAPxx        | Map of defined pixel regions for $\nu_\mathrm{p}$ for the local sampler to use. Pixel region index values starts at 1 | 'fullsky' = {1 pixel region on the whole sky} <br> 'none' = {whole sky set to prior} <br> valid filename in DATA_DIRECTORY| 'fullsky'|
| COMP_NU_P_PIXREG_INITVALUE_MAPxx | Map of $\nu_\mathrm{p}$ where each pixelgerion has a uniform value (no smoothing). Used to initialize pixelregions with exact values if the normal spectral index map is smoothed | 'none' = {use original input map to decide pixelregion values} <br> valid filename in DATA_DIRECTORY| 'none'|
| COMP_NU_P_{pol}_NUM_PIXREGxx  | Number of pixel regions to sample, all additional regions are set to prior values (pixel region 0). Only needed for the 'pixreg' pixel region type option for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the proposal length map | > 0 | 10 |
| COMP_NU_P_MASKxx              | Mask for use of the local sampler to sample $\nu_\mathrm{p}$                       | 'fullsky' <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_NU_P_PROPLENxx           | $\nu_\mathrm{p}$ sampling. Map of the standard deviation of each of the Metropolis-Hastings proposals, per pixel region, for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_NU_P_{pol}_PROPLEN_INITxx | (Map of) The initial value of the proposal length of the local sampler for all pixel regions ($\nu_\mathrm{p}$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0.0 <br> <= 0.0 {disable} | 0.05d0 <br> -1.d0 |                   
| COMP_NU_P_{pol}_SAMPLE_PROPLENxx | Flag for whether or not to tune the proposal length of the local sampler for the sampling of the $\nu_\mathrm{p}$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_NU_P_NPROPxx             | $\nu_\mathrm{p}$ sampling. Map of the number of Metropolis-Hastings proposals per pixel region for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|     
| COMP_NU_P_{pol}_NPROP_INITxx  | (Map of) The initial value of the number of proposals of the local sampler for all pixel regions ($\nu_\mathrm{p}$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0 <br> <= 0 {disable} | 100 <br> -1 |                          
| COMP_NU_P_{pol}_SAMPLE_NPROPxx | Flag for whether or not to tune number of proposals of the local sampler for the sampling of the $\nu_\mathrm{p}$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_NU_P_UNI_NPROP_LOWxx     | Absolute lower limit of number of MH proposals per iteration for $\nu_\mathrm{p}$. Covers all Stokes parameters | >= 0 | 1    |
| COMP_NU_P_UNI_NPROP_HIGHxx    | Absolute upper limit of number of MH proposals per iteration for $\nu_\mathrm{p}$. Covers all Stokes parameters | >= 0 | 1000 |


### Tilted Spinning Dust Component

| Parameter name                | Description                                               | Accepted values        | Example value |
| ----------------------------- | --------------------------------------------------------- | ---------------------- | ------------- |
| COMP_INPUT_ALPHA_MAPxx        | Initial map for the Spinning dust $\alpha$ (tilt) spectral parameter        | 'none' = {Default value used} <br> valid filename in DATA_DIRECTORY | 'none' |
| COMP_INPUT_NU_P_MAPxx         | Initial map for the Spinning dust peak frequency $\nu_\mathrm{p}$ spectral parameter        | 'none' = {Default value used} <br> valid filename in DATA_DIRECTORY | 'none' |
| COMP_ALPHA_NU_MINxx            | Lowest frequency for estimating $\alpha$ in GHz             | > 0  | 1. |  
| COMP_ALPHA_NU_MAXxx            | Highest frequency for estimating $\alpha$ in GHz            | > 0  | 100. |
| COMP_ALPHA_POLTYPExx           | Spectral index Stokes parameter coherence switch | 1 = {fit common parameter for T+Q+U}<br>2 = {fit separate parameters for T and Q+U}<br>3 = {fit separate parameters for T, Q and U} | 2 |
| COMP_ALPHA_SMOOTHING_SCALExx   | Smoothing scale ID for $\alpha$; specifies which SMOOTHING_SCALE parameters to use, as defined in Data Set parameters | <= NUM_SMOOTHING_SCALES | 1 |
| COMP_DEFAULT_ALPHAxx           | Initial value of $\alpha$; must lie within uniform prior values  | float | 0 |
| COMP_DEFAULT_NU_Pxx           | Initial value of $\nu_\mathrm{p}$; must lie within uniform prior values  | float | -2. |
| COMP_NU_P_NU_MINxx            | Lowest frequency for estimating $\nu_\mathrm{p}$ in GHz             | > 0  | 1. |  
| COMP_NU_P_NU_MAXxx            | Highest frequency for estimating $\nu_\mathrm{p}$ in GHz            | > 0  | 100. |
| COMP_NU_P_POLTYPExx           | Spectral index Stokes parameter coherence switch | 1 = {fit common parameter for T+Q+U}<br>2 = {fit separate parameters for T and Q+U}<br>3 = {fit separate parameters for T, Q and U} | 2 |
| COMP_NU_P_SMOOTHING_SCALExx   | Smoothing scale ID for $\nu_\mathrm{p}$; specifies which SMOOTHING_SCALE parameters to use, as defined in Data Set parameters | <= NUM_SMOOTHING_SCALES | 1 ||COMP_PRIOR_GAUSS_ALPHA_MEANxx  | Mean of Gaussian prior on $\alpha$                         | float                     | 0   |
| COMP_PRIOR_GAUSS_ALPHA_RMSxx   | Standard deviation (RMS) of Gaussian prior on $\alpha$     | <0 = {Do not apply Gaussian prior}<br>0 = {Fix $\alpha$ on input}<br> >0 = {Sample with specifed prior RMS} | 0.5 |
| COMP_PRIOR_GAUSS_NU_P_MEANxx  | Mean of Gaussian prior on $\nu_\mathrm{p}$                         | float                     | 21.   |
| COMP_PRIOR_GAUSS_NU_P_RMSxx   | Standard deviation (RMS) of Gaussian prior on $\nu_\mathrm{p}$     | <0 = {Do not apply Gaussian prior}<br>0 = {Fix $\nu_\mathrm{p}$ on input}<br> >0 = {Sample with specifed prior RMS} | 0.2 |
| COMP_PRIOR_UNI_ALPHA_LOWxx     | Absolute lower bound on $\alpha$; hard uniform prior      | float                      | -1.   |
| COMP_PRIOR_UNI_ALPHA_HIGHxx    | Absolute upper bound on $\alpha$; hard uniform prior       | > COMP_PRIOR_UNI_ALPHA_LOW | 1.   |
| COMP_PRIOR_UNI_NU_P_LOWxx     | Absolute lower bound on $\nu_\mathrm{p}$; hard uniform prior      | float                      | 10.   |
| COMP_PRIOR_UNI_NU_P_HIGHxx    | Absolute upper bound on $\nu_\mathrm{p}$; hard uniform prior       | > COMP_PRIOR_UNI_NU_P_LOW | 30.   |
| COMP_SED_TEMPLATExx           | SED template file; must be defined in flux density units (MJy/sr) | valid filename in DATA_DIRECTORY | spdust2_cnm.dat |
| COMP_NU_P_{pol}_LMAXxx        | Maximum multipole moment of the $\nu_\mathrm{p}$ spectral index map for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| <0 : {Use local pix-by-pix sampling} <br> >=0 : {Use $a_{lm}$ sampling}| -1 <br> 0 |
| COMP_NU_P_{pol}_PIXREGxx      | Type of pixel regions for the local sampling of the $\nu_\mathrm{p}$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} | 'pixreg' = pixel region input map <br> 'fullsky' = one fullsky pixel region <br> 'single_pix' = pixel regions as per smoothing scale $N_{\mathrm{side}}$ | 'pixreg' |
| COMP_NU_P_{pol}_LNLTYPExx     | Type of log-likelihood evaluation for sampling of the $\nu_\mathrm{p}$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: using ridge og marginal sampling requires separate CG sampling of the component amplitude afterwards. | 'chisq' <br> 'ridge' <br> 'marginal' | 'chisq' |          
| COMP_NU_P_PIXREG_MAPxx        | Map of defined pixel regions for $\nu_\mathrm{p}$ for the local sampler to use. Pixel region index values starts at 1 | 'fullsky' = {1 pixel region on the whole sky} <br> 'none' = {whole sky set to prior} <br> valid filename in DATA_DIRECTORY| 'fullsky'|
| COMP_NU_P_PIXREG_INITVALUE_MAPxx | Map of $\nu_\mathrm{p}$ where each pixelgerion has a uniform value (no smoothing). Used to initialize pixelregions with exact values if the normal spectral index map is smoothed | 'none' = {use original input map to decide pixelregion values} <br> valid filename in DATA_DIRECTORY| 'none'|
| COMP_NU_P_{pol}_NUM_PIXREGxx  | Number of pixel regions to sample, all additional regions are set to prior values (pixel region 0). Only needed for the 'pixreg' pixel region type option for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the proposal length map | > 0 | 10 |
| COMP_NU_P_MASKxx              | Mask for use of the local sampler to sample $\nu_\mathrm{p}$                       | 'fullsky' <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_NU_P_PROPLENxx           | $\nu_\mathrm{p}$ sampling. Map of the standard deviation of each of the Metropolis-Hastings proposals, per pixel region, for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_NU_P_{pol}_PROPLEN_INITxx | (Map of) The initial value of the proposal length of the local sampler for all pixel regions ($\nu_\mathrm{p}$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0.0 <br> <= 0.0 {disable} | 0.05d0 <br> -1.d0 |                   
| COMP_NU_P_{pol}_SAMPLE_PROPLENxx | Flag for whether or not to tune the proposal length of the local sampler for the sampling of the $\nu_\mathrm{p}$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_NU_P_NPROPxx             | $\nu_\mathrm{p}$ sampling. Map of the number of Metropolis-Hastings proposals per pixel region for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|     
| COMP_NU_P_{pol}_NPROP_INITxx  | (Map of) The initial value of the number of proposals of the local sampler for all pixel regions ($\nu_\mathrm{p}$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0 <br> <= 0 {disable} | 100 <br> -1 |                          
| COMP_NU_P_{pol}_SAMPLE_NPROPxx | Flag for whether or not to tune number of proposals of the local sampler for the sampling of the $\nu_\mathrm{p}$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_NU_P_UNI_NPROP_LOWxx     | Absolute lower limit of number of MH proposals per iteration for $\nu_\mathrm{p}$. Covers all Stokes parameters | >= 0 | 1    |
| COMP_NU_P_UNI_NPROP_HIGHxx    | Absolute upper limit of number of MH proposals per iteration for $\nu_\mathrm{p}$. Covers all Stokes parameters | >= 0 | 1000 |
| COMP_ALPHA_{pol}_LMAXxx        | Maximum multipole moment of the $\alpha$ spectral index map for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| <0 : {Use local pix-by-pix sampling} <br> >=0 : {Use $a_{lm}$ sampling}| -1 <br> 0 |
| COMP_ALPHA_{pol}_PIXREGxx      | Type of pixel regions for the local sampling of the $\alpha$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} | 'pixreg' = pixel region input map <br> 'fullsky' = one fullsky pixel region <br> 'single_pix' = pixel regions as per smoothing scale $N_{\mathrm{side}}$ | 'pixreg' |               
| COMP_ALPHA_{pol}_LNLTYPExx     | Type of log-likelihood evaluation for sampling of the $\alpha$ spectral index map for the given Stokes coherence switch index: <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: using ridge og marginal sampling requires separate CG sampling of the component amplitude afterwards. | 'chisq' <br> 'ridge' <br> 'marginal' | 'chisq' |                 
| COMP_ALPHA_PIXREG_MAPxx        | Map of defined pixel regions for $\alpha$ for the local sampler to use. Pixel region index values starts at 1 | 'fullsky' = {1 pixel region on the whole sky} <br> 'none' = {whole sky set to prior} <br> valid filename in DATA_DIRECTORY| 'fullsky'|
| COMP_ALPHA_PIXREG_INITVALUE_MAPxx | Map of $\alpha$ where each pixelgerion has a uniform value (no smoothing). Used to initialize pixelregions with exact values if the normal spectral index map is smoothed | 'none' = {use original input map to decide pixelregion values} <br> valid filename in DATA_DIRECTORY| 'none'|
| COMP_ALPHA_{pol}_NUM_PIXREGxx  | Number of pixel regions to sample, all additional regions are set to prior values (pixel region 0). Only needed for the 'pixreg' pixel region type option for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the proposal length map | > 0 | 10 |                                              
| COMP_ALPHA_MASKxx              | Mask for use of the local sampler to sample $\alpha$                       | 'fullsky' <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_ALPHA_PROPLENxx           | $\alpha$ sampling. Map of the standard deviation of each of the Metropolis-Hastings proposals, per pixel region, for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|
| COMP_ALPHA_{pol}_PROPLEN_INITxx | (Map of) The initial value of the proposal length of the local sampler for all pixel regions ($\alpha$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0.0 <br> <= 0.0 {disable} | 0.05d0 <br> -1.d0 |                          
| COMP_ALPHA_{pol}_SAMPLE_PROPLENxx | Flag for whether or not to tune the proposal length of the local sampler for the sampling of the $\alpha$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_ALPHA_NPROPxx             | $\alpha$ sampling. Map of the number of Metropolis-Hastings proposals per pixel region for the local sampler. Covers all Stokes parameters | 'fullsky' = 1.0 <br> valid filename in DATA_DIRECTORY | 'fullsky'|            
| COMP_ALPHA_{pol}_NPROP_INITxx  | (Map of) The initial value of the number of proposals of the local sampler for all pixel regions ($\alpha$ spectral index map), for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$} <br> Note: this overwrites the values read from the input map| > 0 <br> <= 0 {disable} | 100 <br> -1 |                                 
| COMP_ALPHA_{pol}_SAMPLE_NPROPxx | Flag for whether or not to tune number of proposals of the local sampler for the sampling of the $\alpha$ spectral index map, for the given Stokes coherence switch index <br> {pol} = INT {Stokes coherence index 1, $IQU$ or $I$} <br> {pol} = POL {Stokes coherence index 2, $QU$ or $Q$} <br> {pol} = POL3 {Stokes coherence index 3, $U$}| .true. <br> .false. | .false. |
| COMP_ALPHA_UNI_NPROP_LOWxx     | Absolute lower limit of number of MH proposals per iteration for $\alpha$. Covers all Stokes parameters | >= 0 | 1    |
| COMP_ALPHA_UNI_NPROP_HIGHxx    | Absolute upper limit of number of MH proposals per iteration for $\alpha$. Covers all Stokes parameters | >= 0 | 1000 |


## Point Source Components

All compact object components are described through the same set of general parameters, as listed below.

| Parameter name                | Description                                               | Accepted values    |  Example value |
| ----------------------------- | --------------------------------------------------------- | ------------------------ | :---------: |
| COMP_ALPHA_NU_MINxx            | Lowest frequency for estimating $\alpha$ in GHz             | > 0  | 1. |  
| COMP_ALPHA_NU_MAXxx            | Highest frequency for estimating $\alpha$ in GHz            | > 0  | 100. |
| COMP_ALPHA_POLTYPExx           | Spectral index Stokes parameter coherence switch | 1 = {fit common parameter for T+Q+U}<br>2 = {fit separate parameters for T and Q+U}<br>3 = {fit separate parameters for T, Q and U} | 2 |
| COMP_AMP_RMS_SCALE_FACTORxx   | Scale factor for amplitude RMS values; useful for debugging and testing | >= 0 | 1. |
| COMP_APPLY_JEFFREYS_PRIORxx   | Not yet active                             | .true.<br>.false. | .false. |
| COMP_APPLY_POSITIVITY_PRIOR   | Enforce positive flux density amplitudes, $a\ge0$  | .true.<br>.false.  | .true. |
| COMP_BETA_NU_MINxx            | Lowest frequency for estimating $\beta$ in GHz             | > 0  | 1. |  
| COMP_BETA_NU_MAXxx            | Highest frequency for estimating $\beta$ in GHz            | > 0  | 1000. |
| COMP_BETA_POLTYPExx           | Spectral index Stokes parameter coherence switch | 1 = {fit common parameter for T+Q+U}<br>2 = {fit separate parameters for T and Q+U}<br>3 = {fit separate parameters for T, Q and U} | 2 |
| COMP_BURN_IN_ON_FIRST_SAMPLE  | Run MH sampler 10 times longer first iteration than normal; useful for speeding up burn-in | .true.<br>.false. | .true. |
| COMP_CATALOGxx                | Source catalog in ASCII format                            | valid filename in DATA_DIRECTORY | COM_AT20G_GB6_NVSS_PCCS2_v1.dat |
| COMP_INIT_FROM_HDFxx          | Initialization mode; if a previous chain file is specified, it must follow the same convention as INIT_CHAIN         | none={use values from INSTRUMENT_PARAM_FILE}<br>default={use default chain specified in INIT_CHAIN}<br>*previous chain file*={use specified sample} | "chain_init_v2.h5:5" |
| COMP_MIN_DIST_BETWEEN_SRC06   | Minimum allowed distance between two sources in arcmin | >= 0 | 0. |
| COMP_NSIDExx                  | HEALPix resolution parameter for output maps  | > 0 | 512 |
| COMP_NU_REFxx                 | SED reference frequency in GHz      | > 0 | 30. |
| COMP_OUTPUT_EB_MAPXX          | Output $E$ and $B$ maps in addition to $Q$ and $U$                  | .true.<br>.false. | .false. |
| COMP_OUTPUT_FWHMxx            | Gaussian smoothing scale in arcmin applied to output FITS maps; only for visualization | >= 0 | 60. |
| COMP_OUTPUT_PTSRC_TEMPLATEXX  | Request computation of source template; set to .false. if catalog already exists | .true.<br>.false. | .false. |
| COMP_POLARIZATIONxx           | Polarization switch                                       | .true.<br>.false.  | .true. |
| COMP_POLTYPExx                | Spectral index Stokes parameter coherence switch | 1 = {fit common parameter for T+Q+U}<br>2 = {fit separate parameters for T and Q+U}<br>3 = {fit separate parameters for T, Q and U} | 1 |
| COMP_PTSRC_TEMPLATEXX         | HDF file with pre-computed beam templates for each source in the catalog, evaluated for each active frequency channel | valid filename in DATA_DIRECOTRY | COM_AT20G_GB6_NVSS_PCCS2_v1.h5 |
| COMP_UNITxx                   | Physical unit of amplitude map | uK_cmb, mK_cmb, K_cmb<br>MJy/sr<br>uK_RJ<br>K km/s | uK_cmb |


## Template Components

All template components are described through the same set of general parameters, as listed below.

| Parameter name                | Description                                               | Accepted values    |  Example value |
| ----------------------------- | --------------------------------------------------------- | ------------------------ | :---------: |
| COMP_DEFAULT_AMPLITUDExx      | Default template amplitude | float | 1.d0 |
| COMP_PRIOR_GAUSS_MEANxx       | Gaussian template amplitude prior mean  | float | 1.0 |
| COMP_PRIOR_GAUSS_RMSxx        | Gaussian template amplitude prior standard deviation (RMS)   | <0 = {No Gaussian prior}<br>0 = {Fix amplitude on input}<br> >0 = {use specified value as prior RMS} | 0.1 |
| COMP_TEMPLATE_DEFINITION_FILExx | ASCII file specifying whether template is active for each band | valid filename in DATA_DIRECTORY | cmb_relquad_def.txt |
