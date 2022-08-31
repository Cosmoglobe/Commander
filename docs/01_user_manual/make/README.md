# Manual compilation

> **Note**: Before proceeding, we highly recommend **removing**/**commenting** out 
  everything related to **Anaconda**/**Miniconda** in your `.bashrc` (or other shell).
  The presence of Anaconda in the `PATH` often leads to conflicts for many libraries, 
  and results in compilation failure(s). Therefore, disable it during installation.
  Once Commander binary is produced you can safely put it back.

[**TODO**]: Zlib, SZIP, OpenBLAS+FFTW3 or MKL, HDF5, CFITSIO, HEALPix, CAMB, Commander 

- [LAPACK](http://www.netlib.org/lapack/) - required regardless of installation type;
- [HDF5](https://www.hdfgroup.org/) - required only if compiled via `Makefile`;
- [FFTW](http://www.fftw.org/) - required only if compiled via `Makefile`;
- [Sharp2](https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/tree/master) - required only if compiled via `Makefile`;
- [Healpix](https://healpix.sourceforge.io/) - required only if compiled via `Makefile`;
- [CFitsio](https://heasarc.gsfc.nasa.gov/fitsio/) - required only if compiled via `Makefile`;


##### Zlib 

[**TODO**]: Write this

##### SZip

[**TODO**]: Write this

##### HDF5 

[**TODO**]: Write this

##### CFitsIO 

[**TODO**]: Write this

##### HEALPix

[**TODO**]: Write this

##### LAPACK and FFT

<details>
<summary>
<b>OpenBLAS + FFTW</b>
</summary>

Head over to the official GitHub repository and grap latest release, e.g.:
```
$ wget https://github.com/xianyi/OpenBLAS/archive/v0.3.12.tar.gz
$ tar -xzvf v0.3.12.tar.gz
$ cd OpenBLAS-0.3.12
```
And run Compile and install commands as follows:
```
$ make USE_OPENMP=1 -j 48
$ make PREFIX=$HOME/local/openblas install
```
Lastly, update your `.bashrc` with:
```
export LD_LIBRARY_PATH="$HOME/local/openblas/lib:$LD_LIBRARY_PATH"
```
and restart your terminal session.
</details>

[TODO]: Describe here how to compile HDF5 and otehr libraries to use with Make/CMake 
from source (Compiler and MPi are in the prerequisites)

##### Commander

[**TODO**]: Rewrite/update this

<details>
<summary>
<b>Commander with Intel</b>
</summary>
</details>

<details>
<summary>
<b>Commander with GNU</b>
</summary>
If you have enough time and desire, you can of course compile Commander3 from scratch. For this you will need to do the following:

1. Determine the locations of your MPI compilers (mpif90, mpif77, mpicc, etc), and ensure that they function correctly;
2. Download and compile the following libraries from source: [HDF5](https://www.hdfgroup.org/) (`version 1.10.0` or higher), [FFTW](http://www.fftw.org/) (`version 3.3.8` or higher), [CFitsio](https://heasarc.gsfc.nasa.gov/fitsio/) (`version 3.47` or higher) and [HEALPix](https://healpix.sourceforge.io/) (`version 3.70` or higher);
3. Look in the `config/` directory and see if a configuration already exists which is similar to your machine.  Copy the `config` (or the `config.example`) to a new file in the same directory.  Call this new file `config.<machine>` where `<machine>` is the name of the system you will be building the software on.
4. Edit the `config` file and specify all the options for your system. Here you will need to specify the locations of newly installed libraries as well as linking rules;
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
</details>
