<a name="top"></a>
<p align="center">
    <img src="https://github.com/hke/Commander/blob/master/logo/Commander-logo-large-1024x335.png" height="150">
</p>

**Commander** is an **O**ptimal **M**onte-carlo **M**arkov ch**A**i**N** **D**riven **E**stimato**R** which implements fast and efficient end-to-end CMB posterior exploration through Gibbs sampling.

---

| [Main features](#main-features) |  [Installation](#installation) | [Usage](#usage) | [License](#license) | [Releases](#releases)| [Projects](#projects) | [Funding](#funding) | [Citation](#citation) |

---

## Main Features

The latest version - `Commander3` - brings together critical features such as:

- Modern Linear Solver
- Map Making
- Parallelism
- Sky and instrumental modelling
- CMB Component Separation

`Commander3` is written using modern `Fortran` standards such as modules, sub modules, and object oriented derived types. The code is highly tuned and optimized to run on High Performance Computing (HPC) facilities, but it can also be run on your local machine.

The previous incarnation of **Commander**, - `Commander2` - is now an internal part of `Commander3`, while the first version of the code, - `Commander1` - is used mainly for debugging and/or legacy purposes. However, `Commander1` has not been officially released; thus, it doesn't support [CMake](https://cmake.org/) installation, as described in [official documentation](https://docs.beyondplanck.science/#/parameters/intro).

---

## Installation

For the complete installation guide please refer to the [official documentation](https://docs.beyondplanck.science/#/parameters/intro), where you can find how to compile and run `Commander` on different platforms, including HPCs such as NERSC, UNINETT Sigma2, OWLs etc. Below you can find the short summary of how to compile it from source.

### Prerequisites

To successfully run Commander, you need the following libraries:

- [MPI]() - required regardless of installation type;
- [OpenMP]() - required regardless of installation type;
- [BLAS]() - required regardless of installation type;
- [LAPACK](http://www.netlib.org/lapack/) - required regardless of installation type;
- [HDF5](https://www.hdfgroup.org/) - required only if compiled via `Makefile`;
- [FFTW](http://www.fftw.org/) - required only if compiled via `Makefile`;
- [Sharp2](https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/tree/master) - required only if compiled via `Makefile`;
- [Healpix](https://healpix.sourceforge.io/) - required only if compiled via `Makefile`;
- [CFitsio](https://heasarc.gsfc.nasa.gov/fitsio/) - required only if compiled via `Makefile`;

In addition you may want to install/update the following packages:

- Automake version 1.16 or higher - required regardless of installation type;
- Autoconf version 2.69 or higher - required regardless of installation type;
- Libtool version 2.4.6 or higher - required regardless of installation type;

### Compile using CMake

[CMake](https://cmake.org/) is a tool which allows you to compile your code on various platform, via generation of build files (e.g. on Linux are `Makefiles`). It is configured to scan your system and identify present/missing libraries to download and install the missing ones. So, please install CMake before proceeding by this installation type.

Once CMake is installed, the `Commander3` installation procedure is quite simple and consists of the following steps:
```
$ git clone https://github.com/Cosmoglobe/Commander.git
$ cd Commander
$ mkdir build
$ cd build
```
then to configure `Commander3` compillation with, e.g. Intel Fortran compilers, use:
```
$ cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort ..
```
wait while configuration is finished and then run:
```
$ cmake --build . -j n
```
where `n` is the amount of processors you wish to use to speed up the installation.

Because `Commander3` is usually run on HPCs, where users do not have the `sudo`/`root` previleges, the default installation path is configured to be inside `/path/to/cloned/repo/Commander/build/install/`, where `Commander3` binary can be found inside
`bin` folder, under the name `commander3`.

### Compile from source

After you cloned the latest version of the repo, you need to do the following:

1. Determine the locations of your MPI compilers (mpif90, mpif77, mpicc, etc), and ensure that they function correctly.
2. Determine the locations of the CFITSIO and LAPACK libraries, and how to link to these libraries.
3. Look in the `config/` directory and see if a configuration already exists which is similar to your machine.  Copy the `config` (or the `config.example`) to a new file in the same directory.  Call this new file `config.<machine>` where `<machine>` is the name of the system you will be building the software on.
4. Edit the `config` file and specify all the options for your system.
5. `cd` into the top level commander directory and set the `COMMANDER` environment variable to the string you used above as the `<machine>`.  For example, if your config is named `config.mylaptop`, then you would set the environment variable in `.bashrc` as:
    - For Bourne-like shells
    ```
    $ export COMMANDER=mylaptop
    ```
    - For csh-like shells 
    ```
    % setenv COMMANDER mylaptop
    ```
6. To view some additional help about the available make targets:
```
$ make help
```
7. To build and install the software:
```
$ make
$ make install
```

---

## Usage

In short, to run `Commander3`, you need to:
```
$ export OMP_NUM_THREADS=1
```
then create the chain directory as specified in parameter file:
```
$ mkdir chains_dir
```
copy parameter file into it:
```
$ cp param_file.txt chains_dir/ 
```
and, finally, run `Commander3` via following command:
```
$ mpirun -np num_proc ~/Commander/src/commander/commander param_file.txt 2>&1 | tee chains_dir/slurm.txt
```
Here, `num_proc` is the number of processors to use, `slurm.txt` is the file to store output logs.

As stated previously, `Commander1` has not been officially released and is used primarily for debugging. If you wish to run it, however, you can compile it with `Makefile` using:
```
$ cd commander1
$ make
```
and then run the following command:
```
$ mpirun -n num_proc ~/Commander/commander1/src/commander/commander param_file.txt 2>&1 | tee chains_dir/slurm.txt
```

---

## Releases

All Cosmoglobe papers are associated with a specific data and code release. All releases can be found [here](https://github.com/Cosmoglobe/Commander/releases) and are also listed below: 

#### [Cosmoglobe DR2](https://github.com/Cosmoglobe/Commander/releases/tag/DR2.1) (Mar 17, 2026)
>DR2.1 is the code used to run the final DIRBE run in early January 2026. It is the version used to produce the results in the [Cosmoglobe DR2 papers](https://www.cosmoglobe.uio.no/products/cosmoglobe-dr2.html).

#### [DR1 Update](https://github.com/Cosmoglobe/Commander/releases/tag/DR1v2) (Jul 9, 2025)
>Updated version of DR1 code that includes spurious mapmaking for WMAP, as well as minor bug fixes and conforming to Fortran standard.

#### [Cosmoglobe DR2 (preliminary)](https://github.com/Cosmoglobe/Commander/releases/tag/DR2) (Jul 12, 2024)
>Production code used for the DIRBE TOD-processing run for the initial set of Cosmoglobe DR2 papers.

#### [Cosmoglobe DR1 ](https://github.com/Cosmoglobe/Commander/releases/tag/v1.0) (Mar 21, 2023)
>Code associated with the publication [Cosmoglobe DR1 results. I. Improved Wilkinson Microwave Anisotropy Probe maps through Bayesian end-to-end analysis and Cosmoglobe Data Release 1.](https://arxiv.org/abs/2303.08095) and [Cosmoglobe Data Release 1](https://www.cosmoglobe.uio.no/products/cosmoglobe-dr1.html).

#### [BeyondPlanck](https://github.com/Cosmoglobe/Commander/releases/tag/BP10) (Nov 19, 2021)
>Code associated with the publication BeyondPlanck.

---

## Projects

Commander framework is part of the following projects:

<p align="center">
    <img src="./logo/Planck_logo.png" height="100"> 
    <img src="./logo/beyondplanck_logo.png" height="100"> 
    <img src="./logo/LiteBIRD-logo-posi-RGB.png" height="100"> 
    <img src="./logo/Cosmoglobe-logo-vertical-large.png" height="100"> 
</p>

---

## Funding

This work has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreements No 819478 (ERC; Cosmoglobe), 772253 (ERC; bits2cosmology), 101165647 (ERC, Origins), and 101141621 (ERC, Commander). This work is also funded by the Research Council of Norway under grant agreement numbers 344934 (YRT; CosmoglobeHD), 360024 (YRT; OpenHFI), 351037 (FRIPRO; LiteBIRD-Norway) as well as grant agreement no. 776282 (COMPET-4; BeyondPlanck).

<p align="center">
    <img src="./logo/LOGO_ERC-FLAG_EU_.jpg" height="200">
    <img src="./logo/horizon2020_logo.jpg" height="200">
</p>

---

## License

[GNU GPLv3](https://github.com/Cosmoglobe/Commander/blob/master/COPYING)

---

## Citation

If used for published results, please [cite these papers](https://github.com/Cosmoglobe/Commander/blob/master/docs/commander.bib):

- [Jewell et al. 2004, ApJ, 609, 1](https://ui.adsabs.harvard.edu/abs/2004ApJ...609....1J)
- [Wandelt et al. 2004, Phys. Rev. D, 70, 083511](https://ui.adsabs.harvard.edu/abs/2004PhRvD..70h3511W)
- [Eriksen et al. 2004, ApJS, 155, 227 (Commander)](https://ui.adsabs.harvard.edu/abs/2004ApJS..155..227E)
- [Eriksen et al. 2008, ApJ, 676, 10  (Joint FG + CMB)](https://ui.adsabs.harvard.edu/abs/2008ApJ...676...10E)

