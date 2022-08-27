# Install Prerequisites

As any other code, Commander relies on a bunch of different libraries and tools in order to 
compile and run. We are trying to keep the list of requirements to a minimum. So far it 
it includes: 

- System Libraries, such as Linux Math Library 
- Compilers for C, C++ and Fortran with OpenMP support 
- MPI implementation
- Git (in order to clone Commander code from GitHub)
- GNU Autotools: Make, Automake, Autoconf, Autoreconf, Libtool, m4
- Python2

In addition, you man need to install CMake if you wish to automate the build.

Below we are giving detailed description for each of the libraries, the list of vendors and how 
to install them on some Linux systems via package managers or from source. 

## Compilers 

We provide support for the following compiler toolchains:

<table style="width:100%">
  <tr>
    <th>Name</th>
    <th>Version</th>
    <th>Description</th>
    <th>Official Website</th>
  </tr>
  <tr>
    <td>Intel Parallel Studio</td>
    <td></td>
    <td>Intel Compilers and Intel MPI from Intel Parallel Studio for C, C++ and Fortran languages</td>
    <td>https://www.intel.com/</td>
  </tr>
  <tr>
    <td>Intel OneAPI</td>
    <td></td>
    <td>Intel OneAPI compilers and MPI for Fortran, C and C++ languages</td>
    <td>https://www.intel.com/</td>
  </tr>
  <tr>
    <td>GNU Compiler Toolchain</td>
    <td>10.x.x and above</td>
    <td>Free and OpenSource compiler toolchain for C, C++, Fortran and other langueages</td>
    <td>https://gcc.gnu.org/</td>
  </tr>
</table>

<b>Note</b>: Even though it is possible to install Commander with the help of GNU Compiler 
toolchain, we recommend using Intel ones since the majority of our group uses it for debugging
and production purposes; hence, they are better tested. If you are on a big cluster, Intel 
compilers are usually installed by your system administrator and can be loaded via Environment
modules. 

##### Installing GNU Compiler Toolchain

<details>
<summary>
From Source 
</summary>
<p align="justify">
Assuming you have working system C/C++ compilers, you can install GNU GCC from source.

For this, do the following: 

1. Head over to the [GNU GCC webpage](https://gcc.gnu.org/), choose one of the 
   [mirrors](https://gcc.gnu.org/mirrors.html) (which is closer to your location), 
   navigate to the `releases` folder and choose the version you are interested in 
   to download. For Linux, we are using `*.tar.gz` archives; thus, copy the link, 
   choose the directory where you want to store source files, download and unpack 
   the archive there.
 
   For example, for version `11.2.0`, we can do it in the terminal as follows:
   ```
   $ cd $HOME/.local/src \ 
   wget https://bigsearcher.com/mirrors/gcc/releases/gcc-11.2.0/gcc-11.2.0.tar.gz && \
   tar -xzvf gcc-11.2.0.tar.gz
   ```
   You can also verify the identity of the file you downloaded using the SHA512 checksum, i.e.
   ```
   $ sha512sum gcc-11.2.0.tar.gz 
   ```
   If the output matches the one presented in the 
   [sha512.sum](https://bigsearcher.com/mirrors/gcc/releases/gcc-11.2.0/sha512.sum) file 
   for this archive, then everything is good. If not, remove the file and repeat the 
   procedure again. 
2. Install the prerequisites as described in the 
   [GNU GCC official installation instructions](https://gcc.gnu.org/install/):
   ```
   $ cd gcc-11.2.0 && ./contrib/download_prerequisites
   ```
3. Create a build directory and run configuration command from there:
   ```
   $ mkdir build && cd build && ../configure --prefix=$HOME/.local/gcc/11.2.0 --enable-languages=all --disable-multilib --enable-threads --enable-checking=release --with-system-zlib
   ```
   where we stated we want to install GCC into specific directory by specifying the 
   `--prefix` flag. You can choose any directory you want.
4. Build the code by simply running:
   ```
   $ make -j N
   ```
   where N is the number of processors to use.

   **Note**: The building process may easily take several hours depending on your system, the 
   amount CPUs you use etc. So, be patient.
5. Install compilers using:
   ```
   $ make install 
   ```
   The resulting binaries will be located inside `$HOME/.local/gcc/11.2.0/bin`.
6. Update your shell environment via, e.g. adding the correct paths into your `.bashrc`:
   ```
   export PATH=$HOME/local/gcc/11.2.0/bin:$PATH
   export MANPATH=$HOME/local/gcc/11.2.0/share/man:$MANPATH
   export INFOPATH=$HOME/local/gcc/11.2.0/share/info:$INFOPATH
   export LD_LIBRARY_PATH=$HOME/local/gcc/11.2.0/lib64:$LD_LIBRARY_PATH 
   ```
   **Note**: do not forget to source it for the changes to take an effect in the current shell:
   ```
   $ source $HOME/.bashrc
   ```
</p>
</details>
</br>


## MPI 

<table style="width:100%">
  <tr>
    <th>Name</th>
    <th>Version</th>
    <th>Description</th>
    <th>Official Website</th>
  </tr>
  <tr>
    <td>Intel Parallel Studio</td>
    <td></td>
    <td>Intel Compilers  and Intel MPI from Intel Parallel Studio for C, C++ and Fortran languages</td>
    <td>https://www.intel.com/</td>
  </tr>
  <tr>
    <td>Intel OneAPI</td>
    <td></td>
    <td>Intel OneAPI compilers and MPI for Fortran, C and C++ languages</td>
    <td>https://www.intel.com/</td>
  </tr>
  <tr>
    <td>OpenMPI</td>
    <td>4.0.x and above</td>
    <td>Free and OpenSource MPI implementation for C, C++ and Fortran Compilers</td>
    <td>https://www.open-mpi.org/</td>
  </tr>
</table>

##### Installing OpenMPI

<details>
<summary>
From Source 
</summary>
<p align="justify">
Assuming you have working version of GCC <code>10.x.x</code> or above, you can install 
OpenMPI to use it for Commander installation and run. 

We will install OpenMPI version of `4.1.4` with the same GNU GCC compilers `11.2.0` we installed
when were describing installation of GNU GCC compilers from source. 

1. Identify directory where you want to store your source files and download the OpenMPI from 
   the [official website](https://www.open-mpi.org/software/ompi/v4.1/). E.g.:
   ```
   $ wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.4.tar.gz
   ```
2. Unpack the archive and move to the root of the source code:
   ```
   $ tar -xzvf openmpi-4.1.4.tar.gz && cd openmpi-4.1.4
   ```
3. Configure the installation:
   ```
   $ ./configure -prefix=$HOME/.local/gcc/11.2.0 --enable-orterun-prefix-by-default CC=gcc CXX=g++ F77=g77 FC=gfortran
   ```
4. Compile the code: 
   ```
   $ make -j N 
   ```
5. Install the code:
   ```
   $ make install
   ```
</p>
</details>
</br>

## CMake

<details>
<summary>
From Source with Make 
</summary>
<p align="justify">
[TODO]
</p>
</details>
</br>

<details>
<summary>
From Source with CMake 
</summary>
<p align="justify">
[TODO]
</p>
</details>
</br>

<details>
<summary>
With Python and <code>pip</code>
</summary>
<p align="justify">
[TODO]
</p>
</details>
</br>


