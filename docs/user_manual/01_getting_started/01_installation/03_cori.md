title: Installing on Cori

## Installing on Cori

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
$ git clone https://github.com/hke/Commander.git
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
