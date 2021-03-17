
# Infrastructure parameters

| Parameter name                 | Description                            | Allowed Values                                                            |Suggested value  |
| ------------------------------ | --------------------------------------------------------------------| ------------------                           |:-----------------:|
| BASE_SEED                      | Random number generator seed                                        | Any integer                                  | 163425          |
| CG_CONV_CHECK_FREQUENCY        | Check convergence of Conjugate Gradient (CG) solver every $n$'th iteration | >= 0                                         | 1               |
| CG_CONVERGENCE_CRITERION       | CG stopping criterion type for joint component separation map-making solver             | residual = {$\|M^{-1}*(Ax-b)\| < \epsilon$}<br>chisquare = {$\Delta\chi^2 < \epsilon$ since last check}<br>fixed_iter = {perform CG_MAXITER steps} | residual for full-sky compsep;<br> fixed_iter for masked CMB constrained realizations |
| CG_INIT_AMPS_ON_ZERO           | Initialize CG search on zero                                       | .true.={initialize with $a=0$}<br>.false. = {initialize with $a = a^{\mathrm{prev}}$}     | .false.                   |
| CG_LMAX_PRECOND                | $\ell_{\mathrm{max}}$ for low-$\ell$ CG preconditioner; only used when CG_PRECOND_TYPE = 'diagonal' | >= 0                      | 30                        |
| CG_MAXITER                     | Maximum allowed number of CG iterations        | >= 1                      | 3000                  |
| CG_MINITER                     | Minimum allowed number of CG iterations        | >= 1                      | 5                     |
| CG_PRECOND_TYPE                | CG preconditioner type                         |pseudoinv = {$M = U^{+}T^{+}(U^{+})^t$}; see Seljebotn et al. 2019<br>diagonal = {$M = A^{-1}_{lm,l'm'}\delta_{ll'}\delta_{mm'}$}; see Eriksen et al. 2004                       | pseudoinv for full-sky compsep;<br>diagonal for masked CMB constrained realizations     |
| CG_TOLERANCE                   | Fractional CG convergence criterion, $\epsilon$    | > 0                      | 1.d-9                     |
| INIT_CHAIN               | Chain file used for default initialization with full path; format is "filename.h5:$n$", where $n$ denotes sample number within chain file      | HDF chain file and sample<br>none                  |  "data/old_chain.h5:1"             |
| ENABLE_TOD_ANALYSIS      | Enable/disable TOD parameter sampling. If true, BAND_TOD_TYPE must be set for all data sets  | .true.<br>.false. | .true.       |
| FFTW3_MAGIC_NUMBERS      | ASCII file with favorable FFT lengths; input TODs are truncated to match the closest lower number. File format: First line defines number of entries in file; subsequent lines give one integer per line (other entries in each line are ignored) | filename |data/fft3_magic_numbers_230810.txt'                   |               
| IGNORE_GAIN_AND_BANDPASS_CORR  | Enable or disable frequency map gain and bandpass corrections; if true, these are set to 1 and 0, respectively  | .true.<br>.false.     | .false.                   |
| MJYSR_CONVENTION | Convention used for bandpass integration involving MJy/sr units | IRAS = {I_nu*nu = const}<br>PSM = {I_nu = const}   |   IRAS          |
| NSIDE_CHISQ                   | HEALPix Nside value for $\chi^2$ maps            | >= 1                  | 16           |
| NUM_GIBBS_ITER           | Target number of full Gibbs iterations; typically the execution is terminated manually before this is actually reached        | >= 1                  | 3000          |
| NUMCHAIN                       | Number of independent Gibbs chains to be produced in parallel       | >= 1                  | 1             |
| OPERATION                      | Posterior mapping mode                 | sample = {Draw samples from full posterior}<br>optimize = {Conditional steepest descent posterior maximization}  | sample|
| OUTPUT_CG_PRECOND_EIGENVALS    | Output eigenvalues of preconditioner to OUTPUT_DIRECTORY/precond_eigenvals.dat; applies only to diagonal preconditioner type       | .true.<br>.false. | .false.       |
| OUTPUT_CHISQ_MAP               | Output $\chi^2$ map per sample to OUTPUT_DIRECTORY            | .true.<br>.false. | .true.        |
| OUTPUT_DEBUG_SEDS              | Output default SED for all components to 'sed.dat'; terminates the program            | .true.<br>.false. | .false.       |
| OUTPUT_DIRECTORY               | Main output directory            | directory name                  | chains    |
| OUTPUT_EVERY_NTH_CG_ITERATION  | Output components every $n$'th CG iteration; useful for assessing CG convergence during testing. Set to 0 for no output.            | >= 0                  | 0             |
| OUTPUT_INPUT_MODEL             | Output sample immediately after initialization with label 999999; terminates the program            | .true.<br>.false. | .false.       |
| OUTPUT_MIXING_MATRIX           | Output mixing matrix for each sample in one HEALPix map per component and data set; use with care, as it will require significant disk space            | .true.<br>.false. | .false.       |
| OUTPUT_RESIDUAL_MAPS           | Output residual maps (= data - model) for each data set and sample            | .true.<br>.false. | .true.        |
| OUTPUT_SIGNALS_PER_BAND        | Output all component maps for each data set and sample            | .true.<br>.false. | .false.       |
| POLARIZATION_CHISQ             | Outputted $\chi^2$ map is either summed over or separate for each Stokes parameter | .true. = {$\chi^2_T,\chi^2_Q,\chi^2_U$}<br>.false.={$\chi^2_T+\chi^2_Q+\chi^2_U$} | .true.        |
| RESAMPLE_CMB                   | Convenience parameter that disables all other parameters than those directly related to the CMB. Useful to generate a larger ensemble of CMB realizations given an existing sample set | .true.<br>.false. | .false. |
| SAMPLE_ONLY_POLARIZATION       | Convenience parameter to disable intensity sampling            | .true.<br>.false. | .false.       |
| SAMPLE_SIGNAL_AMPLITUDES       | Convenience parameter to disable signal amplitude (CG) sampling           | .true.<br>.false. | .true.        |
| SAMPLE_SPECTRAL_INDICES        | Convenience parameter to disable sampling of non-linear parameters (spectral indices, frequency gains etc.)           | .true.<br>.false. | .true.        |
| SET_ALL_NOISE_MAPS_TO_MEAN     | Convenience parameter to set all noise RMS maps to their mean values; useful for debugging and testing the CG algorithm                                       | .true.<br>.false.     | .false.                   |
| TOD_NUM_BP_PROPOSALS_PER_ITER | Number of bandpass proposals per TOD sample  | 0 = {disable bandpass sampling}<br>1={single-proposal Metropolis sampling}<br> >1 = {faster burn-in, but incorrect uncertainties}                  | 1             |
| THINNING_FACTOR                | Output only every $n$'th sample            | >= 0                  | 0             |
| VERBOSITY                      | Terminal output level                  | 0 (minimal) to 3 (debug)                     |    3            |
