title: Installing on OWLs

## Installing on OWLs

BEFORE PROCEEDING PLEASE ENSURE ABSENCE OF `ANACONDA/MINICONDA` IN YOUR PATH!

To succesfully compile Commander 3 on owls you need to load several different modules:
```
$ module load Intel_parallel_studio/2018/3.051
```
This module loads the set of Intel compilers with support for MPI and OpenMP.
```
$ module load gnu
```
This module loads the set of gnu automake, libtool, autoconf, autoreconf etc. necessary to successfully compile LibSharp and other libraries.
```
$ module load git/2.9
```
This module loads Git version 2.9.
```
$ module load hdf5/Intel/1.10.1
```
This module loads hdf5 for intel compilers. In principle, this step can be ommitted as hdf5 can be installed with `Commander3` as well.

Currently, there are 2 `cmake` versions installed - 2.8 and 3.17. To compile commander you need the latter, which has an alias of `cmake3`.

All in all, here is a step-by-step example of building `Commander3` on owls from scratch:
Check the `cmake` version
```
$ which cmake3
```
the output should be `/usr/bin/cmake3`.
```
$ cmake3 --version
```
the output should be `cmake3 version 3.17.2`. If it is correct, you may safely proceed as follows:
```
$ git clone https://github.com/hke/Commander.git
$ cd Commander
$ mkdir build && cd build
$ module load gnu Intel_parallel_studio/2018/3.051 hdf5/Intel/1.10.1
$ cmake3 -DCMAKE_INSTALL_PREFIX=$HOME/local -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DMPI_C_COMPILER=mpiicc -DMPI_CXX_COMPILER=mpiicpc -DMPI_Fortran_COMPILER=mpiifort ..
```
wait for the configuration step to finish and then type
```
$ cmake3 --build . -j n
```
where `n` is the number of processors to use.

Once you compile `Healpix` as described in the previous section, you need to update your `.bashrc` adding a new variable:
```
export HEALPIX=$HOME/local/healpix
```
