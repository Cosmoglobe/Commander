title: Installing on Fram and Saga

## Installing on Fram and Saga

BEFORE PROCEEDING PLEASE ENSURE ABSENCE OF ANACONDA/MINICONDA IN YOUR PATH!

According to official [documentation](https://documentation.sigma2.no/), Fram and Saga are using [`EasyBuild`](https://easybuild.readthedocs.io/en/latest/index.html) to manage all the installation packages. So, first you need to load `EasyBuild` module:
```
$ module load EasyBuild/4.2.2
```
and then either update your `.bashrc` with `EasyBuild` relevant information:
```
export MODULEPATH=$HOME/local/easybuild/modules/all:$MODULEPATH
export EASYBUILD_PREFIX=$HOME/local/easybuild
export EASYBUILD_INSTALLPATH=$HOME/local/easybuild
export EASYBUILD_MODULES_TOOL=Lmod
```
or write the configuration file as described in the [documentation](https://easybuild.readthedocs.io/en/latest/Configuration.html).

Next, we are going to install `Autotools` by first searching via:
```
$ eb -S Autotools
```
and chosing appropriate of suggested list of options:
```
$ eb Autotools-20150215.eb -Dr --job-cores=8
```
Note: we are on the login node; thus, we cannot use too many processors for compilation.

Once it is installed, we can load it as module:
```
$ module load Autotools/20150215
```
Next we need to install `CMake`. Unfortunately, `EasyBuild` offers up to `16.3` version which is still not enough in order to compile `Commander3`. Therefore, we will install it from source:
```
$ wget https://github.com/Kitware/CMake/releases/download/v3.17.3/cmake-3.17.3.tar.gz
$ tar -xzvf cmake-3.17.3.tar.gz
$ cd cmake-3.17.3
$ ./bootstrap --prefix=$HOME/local
$ make -j 8
$ make install
```
With this we are almost ready to compile `Commander3`. What is left is compiler. In this example we will compile it with `Intel`:
```
$ module load intel/2019b HDF5/1.10.5-iimpi-2019b
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
