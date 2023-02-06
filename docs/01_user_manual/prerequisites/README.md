# Installing Prerequisites

> **Note**: Before proceeding, we highly recommend **removing**/**commenting** out 
  everything related to **Anaconda**/**Miniconda** in your `.bashrc` (or other shell).
  The presence of Anaconda in the `PATH` often leads to conflicts for many libraries, 
  and results in compilation failure(s). Therefore, disable it during installation.
  Once Commander binary is produced you can safely put it back.

As any other code, Commander relies on a bunch of different libraries and tools in order to 
compile and run. We are trying to keep the list of requirements to a minimum. So far it 
includes: 

- System Libraries, such as Linux Math Library 
- Compilers for C, C++ and Fortran with OpenMP support 
- MPI implementation
- Git (in order to clone Commander code from GitHub)
- GNU Autotools: Make, Automake, Autoconf, Autoreconf, Libtool, m4
- Python2

In addition, you may need to install [CMake](https://cmake.org/) if you wish to automate 
the build. [**Recommended**]

Below we are giving detailed description for each of the libraries, the list of vendors and how 
to install them on some Linux systems via package managers or from source. 

## Autotools 

We need GNU Autotools of the following versions:

<table style="width:100%">
  <tr>
    <th>Name</th>
    <th>Minimal Version</th>
    <th>Description</th>
    <th>Official Website</th>
  </tr>
  <tr>
    <td>Make</td>
    <td align="center"><code>3.0</code></td>
    <td></td>
    <td>
    <a href="https://www.gnu.org/software/make/">
    https://www.gnu.org/software/make/
    </a></td>
  </tr>
  <tr>
    <td>Autoconf</td>
    <td align="center"><code>2.69</code></td>
    <td></td>
    <td>
    <a href="https://www.gnu.org/software/autoconf/">
    https://www.gnu.org/software/autoconf/
    </a></td>
  </tr>
  <tr>
    <td>Automake</td>
    <td align="center"><code>1.14</code></td>
    <td></td>
    <td>
    <a href="https://www.gnu.org/software/automake/">
    https://www.gnu.org/software/automake/
    </a></td>
  </tr>
  <tr>
    <td>Libtool</td>
    <td align="center"><code>2.4.6</code></td>
    <td></td>
    <td>
    <a href="https://www.gnu.org/software/libtool/">
    https://www.gnu.org/software/libtool/
    <a></td>
  </tr>
</table>

##### Installing Autotools

<details>
<summary>
<b>From Source</b> 
</summary>
<p align="justify">
We are going to install one of the recent versions of GNU Autools.

1. Get GNU Make version `4.3` archive, unpack, configure and install it via:
   ```
   $ wget https://ftp.gnu.org/gnu/make/make-4.3.tar.gz
   $ tar -xzvf make-4.3.tar.gz && cd make-4.3
   $ ./configure --prefix=$HOME/.local/gcc/11.2.0
   $ sh build.sh
   $ ./make install
   ```
2. Doing the same for GNU Libtool of version `2.4.6`:
   ```
   $ wget -O libtool-2.4.6.tar.gz https://ftpmirror.gnu.org/libtool/libtool-2.4.6.tar.gz
   $ tar -xzvf libtool-2.4.6.tar.gz && cd libtool-2.4.6
   $ ./configure --prefix=$HOME/.local/gcc/11.2.0
   $ make && make install
   ```
3. And for Autoconf of version `2.71`:
   ```
   $ wget https://ftp.gnu.org/gnu/autoconf/autoconf-2.71.tar.gz
   $ tar -xzvf autoconf-2.71.tar.gz && cd autoconf-2.71
   $ ./configure --prefix=$HOME/.local/gcc/11.2.0
   $ make && make install 
   ```
4. Lastly, for Automake of version `1.16.4`:
   ```
   $ wget -O automake-1.16.4.tar.gz https://ftp.gnu.org/gnu/automake/automake-1.16.4.tar.gz
   $ tar -xzvf automake-1.16.4.tar.gz && cd automake-1.16.4
   $ ./configure --prefix=$HOME/.local/gcc/11.2.0
   $ make && make install
   ```
5. Update your shell environment via, e.g. adding the correct paths into your `.bashrc`:
   ```
   export PATH=$HOME/.local/gcc/11.2.0/bin:$PATH
   export MANPATH=$HOME/.local/gcc/11.2.0/share/man:$MANPATH
   export INFOPATH=$HOME/.local/gcc/11.2.0/share/info:$INFOPATH
   export LD_LIBRARY_PATH=$HOME/.local/gcc/11.2.0/lib64:$LD_LIBRARY_PATH 
   ```
   **Note**: do not forget to source it for the changes to take an effect in the current shell:
   ```
   $ source $HOME/.bashrc
   ```
</p>
</details>
</br>

## Compilers 

We provide support for the following compiler toolchains:

<table style="width:100%">
  <tr>
    <th>Name</th>
    <th>Minimal Version</th>
    <th>Description</th>
    <th>Official Website</th>
  </tr>
  <tr>
    <td>Intel Parallel Studio</td>
    <td align="center"></td>
    <td>Intel Compilers and Intel MPI from Intel Parallel Studio for C, C++ and Fortran languages</td>
    <td>
    <a href="https://www.intel.com/">
    https://www.intel.com/
    </a></td>
  </tr>
  <tr>
    <td>Intel OneAPI</td>
    <td align="center"></td>
    <td>Intel OneAPI compilers and MPI for Fortran, C and C++ languages</td>
    <td>
    <a href="https://www.intel.com/">
    https://www.intel.com/
    </a></td>
  </tr>
  <tr>
    <td>GNU Compiler Toolchain</td>
    <td align="center"><code>10.x.x</code></td>
    <td>Free and OpenSource compiler toolchain for C, C++, Fortran and other langueages</td>
    <td>
    <a href="https://gcc.gnu.org/">
    https://gcc.gnu.org
    </a></td>
  </tr>
</table>

<blockquote>
<p align="justify">
<strong>Note</strong>: Even though it is possible to install Commander with the help of GNU Compiler 
toolchain, <strong>we recommend</strong> using <strong>Intel</strong> ones since the majority of our group uses it for debugging
and production purposes; hence, they are better tested. If you are on a big cluster, Intel 
compilers are usually installed by your system administrator and can be loaded via Environment
modules.
</p>
</blockquote>

##### Installing GNU Compiler Toolchain

<details>
<summary>
<b>From Source</b> 
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
   export PATH=$HOME/.local/gcc/11.2.0/bin:$PATH
   export MANPATH=$HOME/.local/gcc/11.2.0/share/man:$MANPATH
   export INFOPATH=$HOME/.local/gcc/11.2.0/share/info:$INFOPATH
   export LD_LIBRARY_PATH=$HOME/.local/gcc/11.2.0/lib64:$LD_LIBRARY_PATH 
   ```
   **Note**: do not forget to source it for the changes to take an effect in the current shell:
   ```
   $ source $HOME/.bashrc
   ```
</p>
</details>
</br>


## MPI 

The following MPI implementations proved to be working:

<table style="width:100%">
  <tr>
    <th>Name</th>
    <th>Minimal Version</th>
    <th>Description</th>
    <th>Official Website</th>
  </tr>
  <tr>
    <td>Intel Parallel Studio</td>
    <td></td>
    <td>Intel Compilers  and Intel MPI from Intel Parallel Studio for C, C++ and Fortran languages</td>
    <td>
    <a href="https://www.intel.com/">
    https://www.intel.com/
    </a></td>
  </tr>
  <tr>
    <td>Intel OneAPI</td>
    <td></td>
    <td>Intel OneAPI compilers and MPI for Fortran, C and C++ languages</td>
    <td>
    <a href="https://www.intel.com/">
    https://www.intel.com/
    </a></td>
  </tr>
  <tr>
    <td>OpenMPI</td>
    <td align="center"><code>4.0.x</code></td>
    <td>Free and OpenSource MPI implementation for C, C++ and Fortran Compilers</td>
    <td>
    <a href="https://www.open-mpi.org/">
    https://www.open-mpi.org
    </a></td>
  </tr>
</table>

##### Installing OpenMPI

<details>
<summary>
<b>From Source</b> 
</summary>
<p align="justify">
Assuming you have working version of GCC <code>10.x.x</code> or above, you can install 
OpenMPI to use it for Commander installation and run.

We will install OpenMPI version of `4.1.4` with the same GNU GCC compilers `11.2.0` we installed
when were describing installation of GNU GCC compilers from source. We will also include 
support for the [OpenUCX](https://github.com/openucx/ucx) library for OpenMPI.

1. Identify directory where you want to store your source files and download the 
   OpenUCX source from:
   ```
   $ wget https://github.com/openucx/ucx/releases/download/v1.12.1/ucx-1.12.1.tar.gz 
   $ tar -xzvf ucx-1.12.1.tar.gz && cd ucx-1.12.1 
   ```
2. Configure it with 
   ```
   $ ./configure CC=gcc CXX=g++ --prefix=$HOME/.local/gcc/11.2.0 --enable-mt
   ```
3. And compile/install it using `make`:
   ```
   $ make -j N && make install
   ```
   where `N` is the number of processors to use.
4. Download OpenMPI sources from the 
   [official website](https://www.open-mpi.org/software/ompi/v4.1/), unpack and move 
   to the root of the source code:
   ```
   $ wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.4.tar.gz
   $ tar -xzvf openmpi-4.1.4.tar.gz && cd openmpi-4.1.4
   ```
5. Configure the installation:
   ```
   $ ./configure --prefix=$HOME/.local/gcc/11.2.0 --enable-orterun-prefix-by-default CC=gcc CXX=g++ F77=g77 FC=gfortran --with-ucx=$HOME/.local/gcc/11.2.0
   ```
6. Compile/install the code: 
   ```
   $ make -j N && make install 
   ```
</p>
</details>
</br>

## CMake

The minimal CMake version we require as of now is `3.21.0`. It is, however, encouraged 
to install the most recent one. 

##### Installing CMake

<details>
<summary>
<b>From Source with Make</b> 
</summary>
<p align="justify">
We will install CMake version <code>3.21.3</code> with the same GNU GCC compilers 
<code>11.2.0</code> we installed when were describing installation of GNU GCC compilers 
from source. 

Overall, the procedure is follows:
1. Choose directory where to store source files:
   ```
   $ cd $HOME/.local/src 
   ```
2. Download minimum required CMake version and unpack it:
   ```
   $ wget https://github.com/Kitware/CMake/releases/download/v3.21.3/cmake-3.21.3.tar.gz 
   $ tar -xzvf cmake-3.21.3.tar.gz && cd cmake-3.21.3
   ```
3. Configure the CMake to be compiled with Release version: 
   ```
   $ ./bootstrap --prefix=$HOME/.local/gcc/11.2.0 -- -DCMAKE_BUILD_TYPE:STRING=Release 
   ```
4. Compile it using <code>make</code>: 
   ```
   $ make -j 8
   ```
5. And install it:
   ```
   $ make install
   ```
   **Note**: You may need to update your `.bashrc` (or other profile) variables to point to 
   the recent CMake installation.
</p>
</details>
</br>

<details>
<summary>
[TODO]: From Source with CMake 
</summary>
<p align="justify">
[TODO]
</p>
</details>
</br>

<details>
<summary>
<b>With Python and <code>pip</code></b>
</summary>
<p align="justify">
In your <code>python3</code> virtual environment type:
<pre><code>$ pip install cmake</code></pre>
</p>
</details>
</br>


