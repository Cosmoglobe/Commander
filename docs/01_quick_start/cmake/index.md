# CMake Script

There are two ways to install `Commander3` - via CMake or compiling from source. The first one is easier to follow and is the recommended method. However, it is also possible to compile `Commander3` manually using a traditional `Makefile`. This gives greater control of the process, but it is also more time consuming and error prone.

Regardless of installation method, there are several preliminary requirements that must be met:

- GET RID OF `ANACONDA/MINICONDA` FROM YOUR PATH, LD_LIBRARY_PATH etc. Note: The presence of Anaconda in the path often lead to conflicts for many libraries, and results in compilation failure. Thus, we strongly recommend that you disable this during installation. Once the `commander3` binary is produced, you can safely re-enable your Anaconda installation;
- Intel compiler, for instance ifort, icc, icpc, or {`gcc, g++, gfortran`}. PGI and Clang are still being tested.
- MPI: {`mpiifort, mpiicc, mpiicpc`} for Intel or {`mpifort, mpicc, mpic++`} for GNU;
- OpenMP;
- BLAS and LAPACK;
- Git;
- Automake, `version 1.14 or higher` (required by libsharp2 and other libraries);
- Autoconf, `version 2.69 or higher` (required by libsharp2 and other libraries);
- Libtool, `version 2.4.6 or higher` (required by libsharp2 and other libraries);
- ZLIB

In addition:
- CMake, `version 3.17 or higher` is needed for the CMake installation method;
- Python, `version 2.7 or higher` is not strictly necessary for `Commander3`, but is very useful for running post-processing scripts, such as [`c3pp`](https://www.github.com/Cosmoglobe/c3pp/).


Please ensure that these requirements are met before proceeding. You can find some information on how to install these libraries on Linux below.


#### Installing GNU Gfortran on Linux

In this example, we will show how to install GNU Gfortran from source. For this, you need to have a set of working compiler.

First you head over the GNU GCC webpage to get source code:
https://gcc.gnu.org/

Commander3 was sucessfully compiled with GNU Fortran `v.9.3`, so we are going to use this version in our example. You can choose whether to use one of the mirror sites or Git, we will be using first:
https://gcc.gnu.org/mirrors.html

Choose one of the mirrors, which is closer to your own location and navigate into `releases` folder inside which you are interested in `gcc-9.3.0`. Head over there and choose one of the archives to download. Assuming you want to download your code inside `$HOME/local/src`, it can be done from terminal like this (using `wget`, but you can use `curl`):
```
$ cd $HOME/local/src
$ wget ttps://mirror.koddos.net/gcc/releases/gcc-9.3.0/gcc-9.3.0.tar.gz
```
Now, unpack the archive:
```
$ tar -xzvf gcc-9.3.0.tar.gz
```
and move into the created directory:
```
$ cd gcc-9.3.0
```
GNU GCC relies on some additional libraries, which needs to be installed BEFORE the configuration. It also provides script on how to do it. So, inside GCC `root` (your current) folder run this command:
```
$ ./contrib/download_prerequisites
```
Wait until this process finishes. Once it is done, you need to create a `build` directory and run your GCC configuration from it, i.e.:
```
$ mkdir build && cd build
$ ../configure --prefix=$HOME/local/gcc/9.3.0 --enable-languages=all --disable-multilib --enable-threads --enable-checking=release --with-system-zlib --enable-plugin --enable-initfini-array --disable-libgcj --enable-gnu-indirect-function
```
adn finally start compilation as usual:
```
$ make -j 48
$ make install
```
It will take some time to install it. Once it is done, you need to update your `$PATH` and several other variables inside your `.bashrc` to point for compiler installation path, e.g.:
```
export PATH="$HOME/local/gcc/9.3.0/bin:/some/other/paths/:$PATH"
export LD_LIBRARY_PATH="$HOME/local/gcc/9.3.0/lib64:$LD_LIBRARY_PATH"
```
And restart your terminal. You can test your newly installed GNU GCC compiler toolchain via:
```
$ which gfortran
$ gfortran --version
```

#### Installing OpenMPI on Linux

Assuming you have installed GNU GCC as described above, you now need to install OpenMPI. Head over to official website and grab the latest version:
```
$ wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz
$ tar -xzvf openmpi-4.0.5.tar.gz
$ cd openmpi-4.0.5
```
Run configure command, e.g.:
```
$ ./configure --prefix=$HOME/local/openmpi/4.0.5
```
and intall the software:
```
$ make -j 48
$ make install
```
Update your `.bashrc` pointing to your OpenMPI installation, e.g.:
```
export PATH="$HOME/local/openmpi/4.0.5/bin:/some/other/path:$PATH"
export LD_LIBRARY_PATH="$HOME/local/openmpi/4.0.5/lib:$LD_LIBRARY_PATH"
```
and restart your terminal.

#### Installing OpenBLAS on Linux

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

# CMake version

The default `Commander3` installation process is fully automated through [CMake](https://cmake.org/), which is an Open Source, cross-platform family of tools designed to build, test and package software. Our script is configured to scan your system and to identify both present and missing libraries. Once that step is done, missing libraries will be downloaded, compiled and linked with `Commander3`, together with those libraries that may already be installed on the system.

The full list of required libraries is the following:

- [HDF5](https://www.hdfgroup.org/), `version 1.10.5`, `minimum required 1.10.0`;
- [FFTW](http://www.fftw.org/), `version 3.3.8`, `minimum required 3.3.8`;
- [Sharp2](https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/tree/master). Note: This is now provided as an internal part of `HEALPix v3.60`, but is also included here for debugging purposes;
- [HEALPix](https://healpix.sourceforge.io/), `version 3.60`, `minimum required 3.50`;
- [CFitsio](https://heasarc.gsfc.nasa.gov/fitsio/), `version 3.47`, `minimum required 3.47`;

#### How does it work?

From a user's perspective, CMake works similarly to `configure` scripts, in the sense that the installation process is divided into two stages, namely configuration and compilation. Typically, when compiling a library from sources, one performs the following steps:
```
$ [variable 1] [variable 2] [...] [variable n] ./configure [command 1] [command 2] [...] [command n]
$ make [make command 1] [make command 2] [...] [make command n]
$ make install
```
where `[variable n]` might be `FC=ifort`, `[command n]` might be `--prefix=/path/where/to/install/software`, `[make command 1]` can be `-j 8` (in which case 8 processors will be used to speed up the compilation process). The first line will configure the project, while the second and third lines will compile and install.

The similar process in CMake is the following:
```
$ mkdir [build directory] && cd [build directory]
$ cmake [CMake variable 1] [CMake variable 2] [...] [CMake variable n] [Path to the main CMakeLists.txt]
$ cmake --build . --target install [CMake command 1] [...] [CMake command n]
```
Because compiling sources within the source tree is generally considered to be bad practice, you are strongly recommended to create a dedicated `build` directory (the name can be anything, but it is conventional to include "build" in the full name), and then run CMake configuration from that directory.

The second line is the main configuration step, where
- `[CMake variable 1]` can be `-DCMAKE_Fortran_COMPILER=ifort` (if you want to use Intel Fortran compiler)
- `[CMake variable 2]` can be `-DCMAKE_INSTALL_PREFIX=/path/where/to/install/software`,
- `[CMake variable n]` can be `-L` (which will print to the terminal the list of all defined variables of the current project at the bottom of your configuration),
- `[CMake command 1]` can be `-j 8`,
- `[CMake command n]` can be `-v` (similar to `VERBOSE=1` from `make` which is very useful for debugging as it shows you all linking/compiling flags).

If you, as recommended, created a `build` directory inside the root folder of your project, you will specify the main path as `..`, while if you created your `build` directory somewhere else, you need to specify the complete path to the main `CMakeLists.txt` (located in the `root` of the project).

The third line is the combination of `make -j 8` and `make install`.

A list of all currently defined variables is provided below. Please note that, although we are trying as best as we can to keep this documentation up to date, due to the constant development this list can be removed/modified/updated. Therefore, use `-L` option mentioned above during configuration to see the full list of available variables.


##### Installation path variables

With CMake, the user can have a fine grained control over the installation path, i.e. the general prefix or even binary and library folders:

<table style="width:100%">
  <tr>
    <th>Variable</th>
    <th>Type</th>
    <th>Description</th>
	<th>Default value</th>
	<th>Example Usage</th>
  </tr>
  <tr>
    <td><pre><code>CMAKE_INSTALL_PREFIX</code></pre></td>
    <td><pre><code>PATH</code></pre></td>
    <td text-align="justify">Specifies the installation path. The default value is system dependent. Note: If you work on a cluster without root priviliges, you need to specify this variable explicitly</td>
	<td><pre><code>/usr/local</code></pre></td>
	<td><pre><code>-DCMAKE_INSTALL_PREFIX=~/local</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>CMAKE_DOWNLOAD_DIRECTORY</code></pre></td>
    <td><pre><code>PATH</code></pre></td>
    <td text-align="justify">Custom defined variable that specifies the download path for Commander dependecies. </td>
	<td><pre><code>/path/to/commander/repo/build/downloads</code></pre></td>
	<td><pre><code>-DCMAKE_DOWNLOAD_DIRECTORY=~/Downloads</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>CMAKE_LIBRARY_OUTPUT_DIRECTORY</code></pre></td>
    <td><pre><code>PATH</code></pre></td>
    <td text-align="justify">Standard CMake variable that specifies the installation path for libraries. Default value is system dependent, and usual location is inside directory specified by <pre><code>CMAKE_INSTALL_PREFIX</code></pre>.</td>
	<td><pre><code>lib</code></pre>, <pre><code>lib64</code></pre></td>
	<td><pre><code>-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=~/local/lib</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>CMAKE_RUNTIME_OUTPUT_DIRECTORY</code></pre></td>
    <td><pre><code>PATH</code></pre></td>
    <td text-align="justify">Standard CMake variable that specifies the installation path for binaries/executables. Default value is system dependent, and usual location is inside directory specified by <pre><code>CMAKE_INSTALL_PREFIX</code></pre>.</td>
	<td><pre><code>bin</code></pre></td>
	<td><pre><code>-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=~/local/bin</code></pre></td>
  </tr>
</table>

##### Subproject variables

We have defined an additional set of variables to control the compilation.

<table style="width:100%">
  <tr>
    <th>Variable</th>
    <th>Type</th>
    <th>Description</th>
	<th>Default value</th>
	<th>Example Usage</th>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_LIBS</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">Enables search for LIBS on the system. If turned <code>FALSE</code>, the search will be skipped and all libraries will be compiled from source.</td>
	<td><pre><code>ON</pre></td>
	<td><pre><code>-DUSE_SYSTEM_LIBS=OFF</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_BLAS</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">Enables search for BLAS/LAPACK on the system. If turned <code>FALSE</code>, the search will be skipped and OpenBLAS library will be compiled from source.</td>
	<td><pre><code>ON</pre></td>
	<td><pre><code>-DUSE_SYSTEM_BLAS=OFF</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_ZLIB</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">Enables search for ZLIB on the system. If turned <code>FALSE</code>, the search will be skipped and ZLib library will be compiled from source.</td>
	<td><pre><code>ON</pre></td>
	<td><pre><code>-DUSE_SYSTEM_ZLIB=OFF</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_LIBAEC</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">Enables search for LIBAEC on the system. If turned <code>FALSE</code>, the search will be skipped and LibAEC library will be compiled from source.</td>
	<td><pre><code>ON</pre></td>
	<td><pre><code>-DUSE_SYSTEM_LIBAEC=OFF</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_HDF5</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">Enables search for HDF5 on the system. If turned <code>FALSE</code>, the search will be skipped and HDF5 library will be compiled from source.</td>
	<td><pre><code>OFF</pre></td>
	<td><pre><code>-DUSE_SYSTEM_HDF5=ON</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_FFTW</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">Enables search for FFTW on the system. If turned <code>FALSE</code>, the search will be skipped and FFTW library will be compiled from source.</td>
	<td><pre><code>OFF</pre></td>
	<td><pre><code>-DUSE_SYSTEM_FFTW=ON</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_HEALPIX</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">Enables search for HEALPix on the system. If turned <code>FALSE</code>, the search will be skipped and HEALPix library will be compiled from source.</td>
	<td><pre><code>OFF</pre></td>
	<td><pre><code>-DUSE_SYSTEM_HEALPIX=ON</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_CFITSIO</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">Enables search for CFITSIO on the system. If turned <code>FALSE</code>, the search will be skipped and CFITSIO library will be compiled from source.</td>
	<td><pre><code>OFF</pre></td>
	<td><pre><code>-DUSE_SYSTEM_CFITSIO=ON</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>CFITSIO_USE_CURL</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">If <code>TRUE</code>, installs CFITSIO with cURL support.</td>
	<td><pre><code>OFF</pre></td>
	<td><pre><code>-DCFITSIO_USE_CURL=ON</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_CURL</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">If evaluates to <code>TRUE</code> and <code>CFITSIO_USE_CURL</code> turned <code>ON</code>, installs CFITSIO with cURL support. cURL will be installed with LibSSH2 and MbedTLS support.</td>
	<td><pre><code>ON</pre></td>
	<td><pre><code>-DUSE_SYSTEM_CURL=OFF</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_MBEDTLS</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">Enables search for MbedTLS on the system. If turned <code>FALSE</code>, the search will be skipped and MbedTLS library will be compiled from source.</td>
	<td><pre><code>ON</pre></td>
	<td><pre><code>-DUSE_SYSTEM_MBEDTLS=OFF</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_LIBSSH2</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">Enables search for LIBSSH2 on the system. If turned <code>FALSE</code>, the search will be skipped and LibSSH2 library will be compiled from source.</td>
	<td><pre><code>ON</pre></td>
	<td><pre><code>-DUSE_SYSTEM_LIBSSH2=OFF</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>DOXYGEN_BUILD_DOCS</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">Enables usage of [Doxygen](https://www.doxygen.no/) to compile out-of-source documentation for <code>Commander3</code>.</td>
	<td><pre><code>OFF</pre></td>
	<td><pre><code>-DDOXYGEN_BUILD_DOCS=ON</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_DOXYGEN</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">If <code>USE_SYSTEM_DOXYGEN</code> evaluates to <code>TRUE</code>, enables search for [Doxygen](https://www.doxygen.no/) on the system. If turned <code>FALSE</code>, the search will be skipped and Doxygen will be compiled from source.</td>
	<td><pre><code>ON</pre></td>
	<td><pre><code>-DUSE_SYSTEM_DOXYGEN=OFF</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_FLEX</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">If <code>USE_SYSTEM_DOXYGEN</code> evaluates to <code>TRUE</code>, enables search for Flex on the system. If turned <code>FALSE</code>, the search will be skipped and Flex will be compiled from source.</td>
	<td><pre><code>ON</pre></td>
	<td><pre><code>-DUSE_SYSTEM_FLEX=OFF</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>USE_SYSTEM_BISON</code></pre></td>
    <td><pre><code>BOOL</code></pre></td>
    <td text-align="justify">If <code>USE_SYSTEM_DOXYGEN</code> evaluates to <code>TRUE</code>, enables search for Bison on the system. If turned <code>FALSE</code>, the search will be skipped and Bison will be compiled from source.</td>
	<td><pre><code>ON</pre></td>
	<td><pre><code>-DUSE_SYSTEM_BISON=OFF</code></pre></td>
  </tr>
</table>

##### Compiler variables

The following variables enforce CMake to use a specific set of compilers:

<table style="width:100%">
  <tr>
    <th>Variable</th>
    <th>Type</th>
    <th>Description</th>
	<th>Default value</th>
	<th>Example Usage</th>
  </tr>
  <tr>
    <td><pre><code>CMAKE_Fortran_COMPILER</code></pre></td>
    <td><pre><code>FILEPATH</code></pre></td>
    <td text-align="justify">Standard CMake variable to specify Fortran compiler.</td>
	<td><pre><code>None</pre></td>
	<td><pre><code>-DCMAKE_Fortran_COMPILER=gfortran</code></pre><pre><code>-DCMAKE_Fortran_COMPILER=/full/path/to/compiler/gfortran</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>CMAKE_C_COMPILER</code></pre></td>
    <td><pre><code>FILEPATH</code></pre></td>
    <td text-align="justify">Standard CMake variable to specify C compiler.</td>
	<td><pre><code>None</pre></td>
	<td><pre><code>-DCMAKE_C_COMPILER=gcc</code></pre><pre><code>-DCMAKE_Fortran_COMPILER=/full/path/to/compiler/gcc</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>CMAKE_CXX_COMPILER</code></pre></td>
    <td><pre><code>FILEPATH</code></pre></td>
    <td text-align="justify">Standard CMake variable to specify C++ compiler.</td>
	<td><pre><code>None</pre></td>
	<td><pre><code>-DCMAKE_CXX_COMPILER=g++</code></pre><pre><code>-DCMAKE_Fortran_COMPILER=/full/path/to/compiler/g++</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>MPI_Fortran_COMPILER</code></pre></td>
    <td><pre><code>FILEPATH</code></pre></td>
    <td text-align="justify">Standard CMake variable to specify MPI Fortran compiler.</td>
	<td><pre><code>None</pre></td>
	<td><pre><code>-DMPI_Fortran_COMPILER=mpifort</code></pre><pre><code>-DCMAKE_Fortran_COMPILER=/full/path/to/compiler/mpifort</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>MPI_C_COMPILER</code></pre></td>
    <td><pre><code>FILEPATH</code></pre></td>
    <td text-align="justify">Standard CMake variable to specify MPI C compiler.</td>
	<td><pre><code>None</pre></td>
	<td><pre><code>-DMPI_C_COMPILER=mpicc</code></pre><pre><code>-DMPI_C_COMPILER=/full/path/to/compiler/mpicc</code></pre></td>
  </tr>
  <tr>
    <td><pre><code>MPI_CXX_COMPILER</code></pre></td>
    <td><pre><code>FILEPATH</code></pre></td>
    <td text-align="justify">Standard CMake variable to specify MPI C++ compiler.</td>
	<td><pre><code>None</pre></td>
	<td><pre><code>-DMPI_CXX_COMPILER=mpicxx</code></pre><pre><code>-DMPI_CXX_COMPILER=/full/path/to/compiler/mpicxx</code></pre></td>
  </tr>
</table>

##### Installation process variables

In CMake, there are four main regimes that can be used to compile your code, each one with different compiler flags that depend on your compiler. This is controlled by the CMake standard variable `CMAKE_BUILD_TYPE`, which can take the following values:

- `Debug` for building your library or executable without optimization and with debug symbols;
- `Release` for building your library or executable with aggressive optimization and without debug symbols;
- `RelWithDebInfo` for building your library or executable with less aggressive optimizations and with debug symbols;
- `MinSizeRel` for building your library or executable with optimizations that do not increase object code size.

For the `Intel` Fortran compiler, we have specified the following regime-dependent flags:

- `Release`: `-Ofast -ipo -xHost -parallel -qopenmp -qopt-matmul -assume byterecl -heap-arrays 16384 -fpe0`;
- `Debug`: `-O0 -g -traceback -parallel -qopenmp -C -assume byterecl -heap-arrays 16384 -fpe0`;
- `RelWithDebInfo`: `-O2 -g -DNDEBUG -parallel -qopenmp -C -assume byterecl -heap-arrays 16384 -fpe0`;
- `MinSizeRel`: `-Os -DNDEBUG -parallel -qopenmp -C -assume byterecl -heap-arrays 16384 -fpe0`.

For the `GNU` compilers, the flags are:

- `Release`: `-O3 -march=native -flto -fopenmp -C`;
- `Debug`: `-O0 -g3 -Wall -Wextra -Wconversion -C -pedantic -fbacktrace -fcheck=bounds -ffpe-trap=zero,overflow,underflow -ffunction-sections -pipe`;
- `RelWithDebInfo`: `-O2 -g -DNDEBUG -fopenmp -C`;
- `MinSizeRel`: `-Os -DNDEBUG -fopenmp -C`;

The default value of `CMAKE_BUILD_TYPE` is `Release`. To overwrite this default, set this variable to another value during configuration, for instance `-DCMAKE_BUILD_TYPE=Debug`.

If you need to add more compiler flags, you can populate the following custom-defined variable:

- `COMMANDER3_Fortran_COMPILER_FLAGS`. Usage: `-DCOMMANDER3_Fortran_COMPILER_FLAGS="flag1;flag2;flag3"`. Please note that different versions of `GNU` require slightly different sets of flags. Thus, this variable's default value for `GNU` is either `-ffree-line-length-none -fallow-argument-mismatch` (for `v10.1` and higher) or `-ffree-line-length-none -Wno-argument-mismatch`. If you specify the custom values for this variable you may need to add those flags as well.

In addition, you can also specify additional flags used at link stage via `COMMANDER3_LINKER_FLAGS` variable. Usage is similar to the one above.

Although it is not recommended, you can completely overwrite the above compiler options with the following variables:

- `COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE` for `Release` configuration. Usage: `-DCOMMANDER3_Fortran_COMPILER_FLAGS_RELEASE="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG` for `Debug` configuration. Usage: `-DCOMMANDER3_Fortran_COMPILER_FLAGS_DEBUG="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO` for `RelWithDebInfo` configuration. Usage: `-DCOMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL` for `MinSizeRel` configuration. Usage: `-DCOMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL="flag1;flag2;flag3"`;

Finally, it is also possible to completely overwrite linker flags via:

- `COMMANDER3_Fortran_LINKER_FLAGS_RELEASE` for `Release` configuration. Usage: `-DCOMMANDER3_Fortran_LINKER_FLAGS_RELEASE="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_LINKER_FLAGS_DEBUG` for `Debug` configuration. Usage: `-DCOMMANDER3_Fortran_LINKER_FLAGS_DEBUG="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO` for `RelWithDebInfo` configuration. Usage: `-DCOMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO="flag1;flag2;flag3"`;
- `COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL` for `MinSizeRel` configuration. Usage: `-DCOMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL="flag1;flag2;flag3"`;

#### Commander3 installation with CMake

Keeping everything above in mind we can, finally, proceed to the `Commander3` installation. A typical insstallation procedure looks like this:
```
$ git clone https://github.com/Cosmoglobe/Commander.git
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
$ cmake --build . --target install -j 8
```
In this example, we are using Intel Parallel Studio and we want to install the `commander3` binary into `~/local/bin`, while compiling with 8 cores. You can, of course, modify this command by adding more variables to have better control over your installation, as described above.

The installation process will take some time, depending on the network bandwidth, number of missing dependencies, and available computing power. Once the process has successfully completed, you need to add a HEALPIX evironment variable to your `.bashrc` file (or similar):
```
export HEALPIX=$HOME/local/healpix
```

# Installation procedures for specific HPC environments

Some specific systems are used by many Commander users, and we provide detailed instructions for some of these here. Note that all these systems have restricted access.

## owl.uio.no -- CMB&CO cluster at the University of Oslo

System information:
- 4 x 72-core nodes with 1.5 TB RAM
- 2 x 64-core nodes with 1.5 TB RAM
- 8 x 24-core nodes with 768 GB RAM
- 5 x 64-core nodes with 256 GB RAM
- For current load, see [owl.uio.no](http://owl.uio.no)

Procedure:
- BEFORE PROCEEDING PLEASE DISABLE `ANACONDA/MINICONDA`!
- Note that there are two `cmake` versions installed on the Owl cluster, namely versions 2.8 and 3.17. To compile commander you need the latter, which has an alias of `cmake3`. Please check that the command `which cmake3` returns `/usr/bin/cmake3`, and that `cmake3 --version` returns `cmake3 version 3.17.2`.
- Load the following modules:
```
$ module load Intel_parallel_studio/2018/3.051
$ module load gnu
$ module load git/2.9
$ module load hdf5/Intel/1.10.1
```
- Run the CMake configuration process:
```
$ git clone https://github.com/Cosmoglobe/Commander.git
$ cd Commander
$ mkdir build && cd build
$ cmake3 -DCMAKE_INSTALL_PREFIX=$HOME/local -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DMPI_C_COMPILER=mpiicc -DMPI_CXX_COMPILER=mpiicpc -DMPI_Fortran_COMPILER=mpiifort ..
```
- Wait for the configuration step to finish and then start the compilation process
```
$ cmake3 --build . -j n
```
where `n` is the number of processors to use.
- Update your `.bashrc` file by adding the following variable:
```
export HEALPIX=$HOME/local/healpix
```

## fram.uio.no and saga.uio.no -- National Norwegian HPC centers

Fram system information:
- 2 x 64-core nodes with 6 TB RAM
- 8 x 32-core nodes with 512 GB RAM
- 998 x 32-core nodes with 64 GB RAM

Saga system information:
- 8 x 64-core nodes with 3 TB RAM
- 28 x 40-core nodes with 384 GB RAM
- 200 x 40-core nodes with 192 GB RAM

Procedure:
- BEFORE PROCEEDING PLEASE DISABLE ANACONDA/MINICONDA!
- [Fram and Saga](https://documentation.sigma2.no/) use [`EasyBuild`](https://easybuild.readthedocs.io/en/latest/index.html) to manage installation packages, and you must therefore first load this module:
```
$ module load EasyBuild/4.2.2
```
- Then either update your `.bashrc` with relevant `EasyBuild` information:
```
export MODULEPATH=$HOME/local/easybuild/modules/all:$MODULEPATH
export EASYBUILD_PREFIX=$HOME/local/easybuild
export EASYBUILD_INSTALLPATH=$HOME/local/easybuild
export EASYBUILD_MODULES_TOOL=Lmod
```
or edit the EasyBuild configuration file as described in the [documentation](https://easybuild.readthedocs.io/en/latest/Configuration.html).
- Install `Autotools`:
```
$ eb -S Autotools
$ eb Autotools-20150215.eb -Dr --job-cores=8
$ module load Autotools/20150215
```
Note: we are on the login node; thus, we cannot use too many processors for compilation.
- Install `CMake`. (Unfortunately, `EasyBuild` currently only offers version `16.3`, which is not sufficient to compile `Commander3`; this is likely to be fixed in a future release.)
```
$ wget https://github.com/Kitware/CMake/releases/download/v3.17.3/cmake-3.17.3.tar.gz
$ tar -xzvf cmake-3.17.3.tar.gz
$ cd cmake-3.17.3
$ ./bootstrap --prefix=$HOME/local
$ make -j 8
$ make install
```
- Run the installation process:
```
$ module load intel/2019b HDF5/1.10.5-iimpi-2019b
$ git clone https://github.com/Cosmoglobe/Commander.git
$ cd Commander
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=$HOME/local -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DMPI_C_COMPILER=mpiicc -DMPI_CXX_COMPILER=mpiicpc -DMPI_Fortran_COMPILER=mpiifort ..
$ cmake --build . -j 8
```
- Lastly, add `HEALPix` into your `.bashrc`:
```
export HEALPIX=$HOME/local/healpix
```

## Cori, NERSC -- National US HPC center

Cori system information:
- 2 x 64-core nodes with 6 TB RAM
- 8 x 32-core nodes with 512 GB RAM
- 998 x 32-core nodes with 64 GB RAM


BEFORE PROCEEDING PLEASE ENSURE ABSENCE OF `ANACONDA/MINICONDA` IN YOUR PATH!

Check for available modules via:
```
$ module avail
```
and also for the list of loaded modules:
```
$ module list
```
By the time of writing of this guide, `CMake 3.17` wasn't available on NERSC. So you need to install it from source via:
```
$ wget https://github.com/Kitware/CMake/releases/download/v3.17.3/cmake-3.17.3.tar.gz
$ tar -xzvf cmake-3.17.3.tar.gz
$ cd cmake-3.17.3
$ ./bootstrap --prefix=$HOME/local
$ make -j 8
$ make install
```
Note: we are on the login node; thus, we cannot use too many processors for compilation.

Once you are done, you can load Intel and Intel MPI via, e.g.:
```
$ module load intel/19.0.3.199 impi/2020
```
Finally, fetch and install `Commander` as:
```
$ git clone https://github.com/Cosmoglobe/Commander.git
$ cd Commander
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=$HOME/local -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DMPI_C_COMPILER=mpiicc -DMPI_CXX_COMPILER=mpiicpc -DMPI_Fortran_COMPILER=mpiifort ..
```
wait for configuration step to finish and then type
```
$ cmake --build . -j 8
```
Lastly, add `HEALPix` into your `.bashrc`:
```
export HEALPIX=$HOME/local/healpix
```
