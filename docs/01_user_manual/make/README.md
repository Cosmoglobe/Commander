# Manual compilation

[TODO]: Describe here how to compile HDF5 and otehr libraries to use with Make/CMake 
from source (Compiler and MPi are in the prerequisites)

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
