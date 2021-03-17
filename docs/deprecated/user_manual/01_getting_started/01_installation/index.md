title: Installation

## Installation

There are two ways to install `Commander3` - via CMake or compiling from source. The first one is easy to follow and is the main way to install the code. However, it is also possible to compile `Commander3` from source using an "old" `Makefile`, although this approach is not recommended for most users as it is more time consuming and error prone.

Regardless of the installation type you've chosen, there is a set of preliminary requirements to install `Commander3`:

- GET RID OF `ANACONDA/MINICONDA` FROM YOUR PATH, LD_LIBRARY_PATH etc. Note: the presence of Anaconda in the path proved to lead to a buggy behavior for many libraries, which results in the failure of compillation. Thus, just comment it out for the course of installation. Once `Commander3` binary is produced, you can safely get your anaconda back;
- Intel compilers: ifort, icc, icpc. Or set of GNU compilers: gcc, g++, gfortran. Note: we are still testing for other compilers such as PGI and Clang;
- MPI: mpiifort, mpiicc, mpiicpc (in case of Intel) or mpifort, mpicc, mpic++ (in case of GNU);
- OpenMP
- BLAS and LAPACK
- Git
- Automake, `version 1.14 or higher`. Note: required by libsharp2 and other libraries;
- Autoconf, `version 2.69 or higher`. Note: required by libsharp2 and other libraries;
- Libtool, `version 2.4.6 or higher`. Note: required by libsharp2 and others;
- ZLIB

In addition, you may/will also need:
- CMake, `version 3.17 or higher`. Note: needed in case you choose to proceed with recommended installation way;
- Python, `version 2.7 or higher`. Note: not necessary for `Commander3` installation, but is useful if you want to run python scripts.

Please ensure they are installed/exist on the system, before proceeding.

### CMake version

It is now possible to completely automate the installation process of `Commander3` thanks to [CMake](https://cmake.org/) - an open-source, cross-platform family of tools designed to build, test and package software. Our script is configured to scan your system and to identify both present and missing libraries. Once it is done, missing libraries will be downloaded, compiled and, together with the ones already installed on the system, will be linked to `Commander3`.

The full list of libraries to be installed is the following:

- [HDF5](https://www.hdfgroup.org/), `version 1.10.5`, `minimum required 1.10.0`;
- [FFTW](http://www.fftw.org/), `version 3.3.8`, `minimum required 3.3.8`;
- [Sharp2](https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/tree/master). Note: comes as an internal part of `HEALPix v3.60`, but is left for debugging purposes;
- [HEALPix](https://healpix.sourceforge.io/), `version 3.60`, `minimum required 3.50`;
- [CFitsio](https://heasarc.gsfc.nasa.gov/fitsio/), `version 3.47`, `minimum required 3.47`;

CMake will look for these libraries and, if failed to identify suitable version, will download and compile it from source before compiling `Commander3`.

#### How does it work?

From a user's perspective, CMake works similarly to `configure` scripts, meaning that the installation can be split into two stages - the configuration and compillation. Usually, when compiling any library or code, you would, generally, do the following:
```
$ [variable 1] [variable 2] [...] [variable n] ./configure [command 1] [command 2] [...] [command n]
$ make [make command 1] [make command 2] [...] [make command n]
$ make install
```
where `[variable n]` can be `FC=ifort`, `[command n]` can be `--prefix=/path/where/to/install/software`, `[make command 1]` can be `-j 8` (which means that 8 processors will be used to speed up the compillation). The first line will configure, the second and the third will compile and install the project.

In terms of CMake, you would do the following:
```
$ mkdir [build directory] && cd [build directory]
$ cmake [CMake variable 1] [CMake variable 2] [...] [CMake variable n] [Path to the main CMakeLists.txt]
$ cmake --build . --target install [CMake command 1] [...] [CMake command n]
```
Because in-source builds generally considered to be a bad practice, you need first to create `build` (the name can be different, but, as a rule of thumb, the word "build" should be present in the full name) directory and then run CMake configuration from that directory.

The second line in the above command is the main configuration step, where `[CMake variable 1]` can be `-DCMAKE_Fortran_COMPILER=ifort` (if you want to use Intel Fortran compiler), while `[CMake variable 2]` can be `-DCMAKE_INSTALL_PREFIX=/path/where/to/install/software`, `[CMake variable n]` can be `-L` (which will print to the terminal the list of all defined variables of the current project at the bottom of your configuration), `[CMake command 1]` can be `-j 8`, `[CMake command n]` can be `-v` (similar to `VERBOSE=1` from `make` which is very useful for debugging as it shows you all linking/compiling flags). If you created a `build` directory inside the root folder of your project (which is a common practice), you will use `..` and if you created your `build` directory somewhere else, you have to specify the complete path to the main `CMakeLists.txt` (located in the `root` of the project).

The third line is the combination of `make -j 8` and `make install`.

Below you may find the list of currently defined variables. Please note that, although we are trying as best as we can to keep this documentation up to date, due to the constant development this list can be removed/modified/updated. Therefore, use `-L` option mentioned above during configuration to see the full list of available variables.

##### Installation path variables

With CMake, the user can have a fine grained control over the installation path, i.e. the general prefix or even binary and library folders:

- `CMAKE_INSTALL_PREFIX` - standard CMake variable to specify the installation path. The default value depends on your system (e.g. on Linux is `/usr/local`). If you are on a Cluster you need to specify this variable explicitly, as you, most probably, will not have access to the root user. Example usage: `-DCMAKE_INSTALL_PREFIX=$HOME/local`;
- `CMAKE_DOWNLOAD_DIRECTORY` - custom defined variable to specify download path for Commander dependecies. The default value is set to `/build/downloads`. If your `build` directory is named anything other than build, you need to modify this variable. Example usage: `-DCMAKE_DOWNLOAD_DIRECTORY=$HOME/Downloads`;
- `CMAKE_LIBRARY_OUTPUT_DIRECTORY` - standard CMake variable to specify the library installation path. Unless explicitly stated, will get the default value depending on the system (e.g. on Linux it is `lib`) and will be located inside the folder specified by `CMAKE_INSTALL_PREFIX`. Example usage: `-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=$HOME/local/lib`;
- `CMAKE_RUNTIME_OUTPUT_DIRECTORY` - similar to the above, but for produced binary/executable. the default value is system dependent (e.g. on Linux it is `bin`) and the folder located inside the one specified by `CMAKE_INSTALL_PREFIX`. Example usage: `-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=~/local/bin`.

##### Variables, which controls the installation process

In CMake, there are 4 main regimes to compile your code, each one of which has different set of compiler flags depending on the compiler you are using. This is controlled by the CMake standard varibale `CMAKE_BUILD_TYPE`, which can have 4 different values:

- `Debug` for building your library or executable without optimization and with debug symbols;
- `Release` for building your library or executable with an aggressive optimization and without debug symbols;
- `RelWithDebInfo` for building your library or executable with less aggressive optimizations and with debug symbols;
- `MinSizeRel` for building your library or executable with optimizations that do not increase object code size.

For `Intel` compiler, we have specified the following Fortran flags depending on the regime:

- `Release`: `-Ofast -ipo -xHost -parallel -qopenmp -qopt-matmul -assume byterecl -heap-arrays 16384 -fpe0`;
- `Debug`: `-O0 -g -traceback -parallel -qopenmp -C -assume byterecl -heap-arrays 16384 -fpe0`;
- `RelWithDebInfo`: `-O2 -g -DNDEBUG -parallel -qopenmp -C -assume byterecl -heap-arrays 16384 -fpe0`;
- `MinSizeRel`: `-Os -DNDEBUG -parallel -qopenmp -C -assume byterecl -heap-arrays 16384 -fpe0`.

The following Fortran flags are for `GNU` compilers:

- `Release`: `-O3 -march=native -flto -fopenmp -C`;
- `Debug`: `-O0 -g3 -Wall -Wextra -Wconversion -C -pedantic -fbacktrace -fcheck=bounds -ffpe-trap=zero,overflow,underflow -ffunction-sections -pipe`;
- `RelWithDebInfo`: `-O2 -g -DNDEBUG -fopenmp -C`;
- `MinSizeRel`: `-Os -DNDEBUG -fopenmp -C`;

The default value of `CMAKE_BUILD_TYPE` is set to `Release`. To overwrite the default value do the following during configuration: `-DCMAKE_BUILD_TYPE=Debug`.

If you feel the need to add more compiler flags you can populate custom defined variable:

- `COMMANDER3_Fortran_COMPILER_FLAGS`. Usage: `-DCOMMANDER3_Fortran_COMPILER_FLAGS="flag1;flag2;flag3"`. Please note that different versions of `GNU` require slightly different sets of flags. Thus, this variable's default value for `GNU` is either `-ffree-line-length-none -fallow-argument-mismatch` (for `v10.1` and higher) or `-ffree-line-length-none -Wno-argument-mismatch`. If you specify the custom values for this variable you may need to add those flags as well.

In addition, you can also specify additional flags used at link stage via `COMMANDER3_LINKER_FLAGS` variable. Usage is similar to the one above.

Although it is not recommended, you can completely overwrite compiler options above with the following variables:

- `COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE` for `Release` configuration. Usage: `-DCOMMANDER3_Fortran_COMPILER_FLAGS_RELEASE="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG` for `Debug` configuration. Usage: `-DCOMMANDER3_Fortran_COMPILER_FLAGS_DEBUG="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO` for `RelWithDebInfo` configuration. Usage: `-DCOMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL` for `MinSizeRel` configuration. Usage: `-DCOMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL="flag1;flag2;flag3"`;

In addition, it is also possible to completely overwrite linker flags via:

- `COMMANDER3_Fortran_LINKER_FLAGS_RELEASE` for `Release` configuration. Usage: `-DCOMMANDER3_Fortran_LINKER_FLAGS_RELEASE="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_LINKER_FLAGS_DEBUG` for `Debug` configuration. Usage: `-DCOMMANDER3_Fortran_LINKER_FLAGS_DEBUG="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO` for `RelWithDebInfo` configuration. Usage: `-DCOMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL` for `MinSizeRel` configuration. Usage: `-DCOMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL="flag1;flag2;flag3"`;

##### Compiler variables

Following variables enforce CMake to use specific set of compilers:

- `CMAKE_Fortran_COMPILER` - standard CMake variable to specify Fortran compiler. Example usage: `-DCMAKE_Fortran_COMPILER=gfortran`;
- `CMAKE_C_COMPILER` - standard CMake variable to specify C compiler. Example usage: `-DCMAKE_C_COMPILER=gcc`;
- `CMAKE_CXX_COMPILER` - standard CMake variable to specify CXX compielr. Example usage: `-DCMAKE_CXX_COMPILER=g++`;
- `MPI_Fortran_COMPILER` - standard CMake variable to specify MPI Fortran compiler. Example usage: `-DMPI_Fortran_COMPILER=mpifort`;
- `MPI_C_COMPILER` - standard CMake variable to specify MPI C compiler. Example usage: `-DMPI_C_COMPILER=mpicc `;
- `MPI_CXX_COMPILER` - standard CMake variable to specify MPI CXX compiler. Example usage: `-DMPI_CXX_COMPILER=mpicxx `;

##### Additional variables

It is possible to ganerate out-of-source documentation using Doxygen - the standard tool for generating documentation from C, C++, Fortran, Python etc. source codes. The following variables controls the flow:

- `DOXYGEN_BUILD_DOCS` - boolean variable which determines whether to build pdf/html documentation out of `Commander3` code or not. The default value is `FALSE`. If set to `TRUE` it will look for the Doxygen on the system and, if not present, will install it similarly to the libraries described above. Example usage: `-DDOXYGEN_BUILD_DOCS=TRUE`.

#### Commander3 installation with CMake

Keeping everything above in mind we can, finally, proceed to the `Commander3` installation. Roughly speaking, you neeed to execute the following set of steps in order:
```
$ git clone https://github.com/hke/Commander.git
$ cd Commander
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=$HOME/local \
		-DCMAKE_C_COMPILER=icc \
		-DCMAKE_CXX_COMPILER=icpc \
		-DCMAKE_Fortran_COMPILER=ifort \
		-DMPI_C_COMPILER=mpiicc \
		-DMPI_CXX_COMPILER=mpiicpc \
		-DMPI_Fortran_COMPILER=mpiifort \
		..
$ cmake --build . -j 8
```
In this example, we are using Intel Parallel Studio and we want to install `commander3` binary into `~/local/bin` using 8 processors. You can, of course, modify this command by adding more variables to have better control over your installation, as described above.

Once you are done, the installation process will start. You just need to wait a little bit for it to finish and, lastly, update your `.bashrc` adding a new variable:
```
export HEALPIX=$HOME/local/healpix
```

Note: we have recently switched to `Healpix v3.70` which allows to make this process fully automated. However, we may need to switch back to `Healpix v3.50`, which doesn't support autoconfiguration script. So, below you will find guidelines on how to configure the library.

##### Healpix 3.50 and below configuration for Commander 3

At some point during installation, it will prompt user for Healpix configuration which needs to be manually provided. Unfortunately, no
autoconfiguration script exists for Healpix version 3.50 and below. So, here is an example of how to install `Healpix v3.50` on Linux
(empty answers indicate the default values, while 3 dots substituted for the output of configuration script):
```
*************************
 Welcome to Healpix 3.50
*************************

This script takes care of the configuration and compilation
of the C, C++, F90, IDL and Python packages of the Healpix distribution.


Do you want to:
(0): exit
(1): configure Healpix IDL package
(2): configure Healpix C package, and edit Makefile
(3): configure Healpix F90 package, and edit Makefile
(4): configure Healpix C++ package, and edit Makefile
(5): configure Healpix Python (healpy) package, and edit Makefile
(8): see what configuration files have been created so far
(9): edit your shell configuration file to have easier access to Healpix codes
(-1): reset
     (will *REMOVE* the Makefile and configuration files, and exit)
(0): exit

Enter your choice (configuration of packages can be done in any order): 3
```
after that, choose intel compiler:
```
you seem to be running Linux
enter name of your F90 compiler (): ifort
```
then you can use default configuration until the C compiler:
```
enter suffix for directories ():
...
Should I attempt to create these directories (Y|n)?
 enter compilation flags for ifort compiler (-I$(F90_INCDIR) -cm -w -sox -qopt-report=0):
enter optimisation flags for ifort compiler (-O3):
...
enter name of your C compiler (cc): icc
...
```
after that you need to specify the location of the freshly installed cfitsio:
```
enter compilation/optimisation flags for C compiler (-O3 -std=c99 -DINTEL_COMPILER):
...
enter command for library archiving (ar -rsv):
enter full name of cfitsio library (libcfitsio.a):
enter location of cfitsio library (/usr/lib64): /path/to/cloned/repo/Commander/build/install/lib
```
Everything else should be default choices:
```
  The generator of non-gaussian CMB maps (ng_sims) can optionally
produce plots of the maps Prob. Dens. Function using PGPLOT.
Do you want to enable this option ?
(this assumes that PGPLOT is already installed on your computer) (y|N)

 The Spherical Harmonics Transform (C and F90) routines used by
synfast/anafast/smoothing/plmgen
and some routines used by ud_grade and alteralm respectively
have a parallel implementation (based on OpenMP).
Do you want to use :
 0) the standard serial implementation ?
 1) the parallel implementation
Enter choice                                      (1):

 Do you want a Position Independent Compilation  (option  "-fPIC")
(recommended if the Healpix-F90 library is to be linked to external codes)  (Y|n):

 Experimental feature:
A static library is produced by default. Do you rather want a shared/dynamic library ? (y|N)

...

The following line should be inserted into your home shell profile (/path/to/your/profile/here/.profile):


 Where the file .healpix/3_50_Linux/config contains:
# configuration for Healpix 3.50

if [ -r ${HPX_CONF_DIR}/idl.sh ] ; then . ${HPX_CONF_DIR}/idl.sh ; fi
if [ -r ${HPX_CONF_DIR}/gdl.sh ] ; then . ${HPX_CONF_DIR}/gdl.sh ; fi
if [ -r ${HPX_CONF_DIR}/fl.sh ]  ; then . ${HPX_CONF_DIR}/fl.sh  ; fi
if [ -r ${HPX_CONF_DIR}/f90.sh ] ; then . ${HPX_CONF_DIR}/f90.sh ; fi
if [ -r ${HPX_CONF_DIR}/cpp.sh ] ; then . ${HPX_CONF_DIR}/cpp.sh ; fi
if [ -r ${HPX_CONF_DIR}/c.sh ] ;   then . ${HPX_CONF_DIR}/c.sh ;   fi

Do you want this modification to be done (y|N)?
```
At this point the configuration is done and you will see the end result:
```
Do you want this modification to be done (y|N)?

Writing pkgconfig file: .../Commander/build/downloads/healpix/src/healpix/lib/healpix.pc
# HEALPix/F90 pkg-config file
# compiled with ifort

prefix=.../Commander/build/downloads/healpix/src/healpix
suffix=
exec_prefix=${prefix}/bin${suffix}
libdir=${prefix}/lib${suffix}
includedir=${prefix}/include${suffix}

Name: HEALPix
Description: F90 library for HEALPix (Hierarchical Equal-Area iso-Latitude) pixelisation of the sphere
Version: 3_50
URL: https://healpix.sourceforge.io
Requires: cfitsio >= 3.20
Libs: -L${libdir} -lhealpix -lhpxgif
Cflags: -I${includedir} -qopenmp -fPIC

           --------------
```
The last thing you need to do is to save your configuration and exit:
```
Do you want to:
(0): exit
(1): configure Healpix IDL package
(2): configure Healpix C package, and edit Makefile
(3): configure Healpix F90 package, and edit Makefile
(4): configure Healpix C++ package, and edit Makefile
(5): configure Healpix Python (healpy) package, and edit Makefile
(8): see what configuration files have been created so far
(9): edit your shell configuration file to have easier access to Healpix codes
(-1): reset
     (will *REMOVE* the Makefile and configuration files, and exit)
(0): exit

Enter your choice (configuration of packages can be done in any order): 0
```
Once it is installed you need to update the `Healpix` environment variable inside your `.bashrc`:
```
export HEALPIX=/path/to/cloned/repo/Commander/build/install/healpix_build
```
