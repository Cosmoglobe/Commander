<a name="top"></a>
<p align="center">
    <img src="_media/commander3-logo.png" height="200">
</p>

# 

# Theoretical Background

[Commander3](https://github.com/Cosmoglobe/Commander) is a Bayesian Markov Chain Monte Carlo sampler for CMB, microwave and sub-mm observations. It fits a user-specified parametric model, $s(\theta)$, to some set of observations, $d_\nu$, by mapping the posterior distribution $P(\theta|d)$ with standard sampling techniques, in particular Gibbs and Metropolis sampling. A concrete example of this is the [BeyondPlanck](https://beyondplanck.science) data model, which reads
$$
d_{j,t} = g_{j,t}\mathsf P_{tp,j}\left[ \mathsf B^{\mathrm{symm}}_{pp',j}\sum_{c} \mathsf M_{cj}(\beta_{p'}, \Delta_\mathrm{bp}^{j})a^c_{p'} + \mathsf B^{\mathrm{asymm}}_{pp',j}\left(s^{\mathrm{orb}}_{j,t} + s^{\mathrm{fsl}}_{j,t}\right)\right] + n^{\mathrm{corr}}_{j,t} + n^{\mathrm{w}}_{j,t}
$$
where $j$ represents a radiometer label, $t$ indicates a single time sample, $p$ denotes a single pixel on the sky, and $c$ represents one single astrophysical signal component, while

- $d_{j,t}$ denotes the measured data value in units of V, this is the calibrated timestream as output from the instrument;
- $g_{j,t}$ denotes the instrumental gain in units of V K$_{\mathrm{cmb}}^{-1}$;
- $\mathsf P_{tp,j}$ is the $N_{\mathrm{TOD}}\times 3N_{\mathrm{pix}}$ pointing matrix stored as compressed pointing and polarization angle timestream per detector;
- $\mathsf B_{pp',j}$ denotes the beam convolution term, where the asymmetric part is only calculated for the orbital dipole and sidelobe terms;
- $\mathsf M_{cj}(\beta_{p}, \Delta_\mathrm{bp})$ denotes element $(c,j)$ of an $N_{\mathrm{comp}}\times N_{\mathrm{comp}}$ mixing matrix, describing the amplitude of component $c$ as seen by radiometer $j$ relative to some reference frequency $j_0$ when assuming some set of bandpass correction parameters $\Delta_\mathrm{bp}$;
- $a^c_{p}$ is the amplitude of component $c$ in pixel $p$, measured at the same reference frequency as the mixing matrix $\mathsf M$;
- $s^{\mathrm{orb}}_{j,t}$ is the orbital CMB dipole signal, including relativistic quadrupole corrections;
- $s^{\mathrm{fsl}}_{j,t}$ denotes the contribution from far sidelobes;
- $n^{\mathrm{corr}}_{j,t}$ denotes correlated instrumental noise;
- $n^{\mathrm{w}}_{j,t}$ is uncorrelated (white) instrumental noise, which is not sampled and is simply left to average down in the maps.

However, the formalism and code are designed to be flexible and extendible, and a wide range of other models may be considered after suitable code modifications.

Following Bayes Theorem, we can write the posterior distribution as
$$
P(\omega\mid \boldsymbol d) = \frac{P(\boldsymbol d\mid \omega)P(\omega)}{P(\boldsymbol d)} \propto \mathcal{L}(\omega)P(\omega),
$$
where $P(\boldsymbol d\mid \omega)\equiv\mathcal{L}(\omega)$ is called the likelihood; $P(\omega)$ is called the prior; and $P(\boldsymbol d)$ is a normalization factor.

Instead of directly explore $P(\boldsymbol d\mid \omega)\equiv\mathcal{L}(\omega)$, [Commander3](https://github.com/Cosmoglobe/Commander) heavily relies on Gibbs sampling theory which states that samples from a joint distribution may be produced by iteratively draw samples from all corresponding conditional distributions. In case of the above data model, this translates into the following Gibbs chain: 
$$
\begin{aligned}
\boldsymbol g &\leftarrow P(\boldsymbol g,\mid \boldsymbol d,\xi_n ,\Delta_\mathrm{bp} ,\boldsymbol a,\beta,C_{\ell})\\
\boldsymbol n_{\mathrm{corr}} &\leftarrow P(\boldsymbol n_{\mathrm{corr}}\mid\boldsymbol d ,\boldsymbol g ,\xi_n, \Delta_\mathrm{bp},\boldsymbol a,\beta ,C_{\ell})\\
\xi_n &\leftarrow P(\xi_n,\mid \boldsymbol d ,\boldsymbol g ,\boldsymbol n_{\mathrm{corr}} ,\Delta_\mathrm{bp} ,\boldsymbol a ,\beta ,C_{\ell})\\
\Delta_\mathrm{bp} &\leftarrow P(\Delta_\mathrm{bp},\mid \boldsymbol d ,\boldsymbol g ,\boldsymbol n_{\mathrm{corr}}, \xi_n, \boldsymbol a ,\beta ,C_{\ell})\\
\beta &\leftarrow P(\beta,\mid \boldsymbol d ,\boldsymbol g ,\boldsymbol n_{\mathrm{corr}} ,\xi_n, \Delta_\mathrm{bp} ,C_{\ell})\\
\boldsymbol a &\leftarrow P(\boldsymbol a,\mid \boldsymbol d ,\boldsymbol g ,\boldsymbol n_{\mathrm{corr}},\xi_n, \Delta_\mathrm{bp}, \beta ,C_{\ell})\\
C_{\ell} &\leftarrow P(C_{\ell},\mid \boldsymbol d ,\boldsymbol g ,\boldsymbol n_{\mathrm{corr}},\xi_n,
\Delta_\mathrm{bp}, \boldsymbol a ,\beta).
\end{aligned}
$$

# Capabilities

The latest version - `Commander3` - brings together critical features such as:

- Modern Linear Solver;
- Map Making;
- Sky and instrumental modelling;
- CMB Component Separation;
- In-memory compression of time-ordered data; 
- FFT optimization; 
- OpenMPI parallelization and load-balancing; 
- I/O efficiency;
- Fast mixing matrix computation through look-up tables.

`Commander3` is written using modern `Fortran` standards such as modules, sub modules, and object oriented derived types. The code is optimized to run on High Performance Computing (HPC) facilities, but it can also be run on your local machine.

The previous incarnation of **Commander**, - `Commander2` - is now an internal part of `Commander3`, while the first version of the code, - `Commander1` - is used mainly for debugging and/or legacy purposes.

# Installation

For installation instructions please refer to the next section.
