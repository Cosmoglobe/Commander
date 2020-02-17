# Commander

An **O**ptimal **M**onte-carlo **M**arkov ch**A**i**N** **D**riven **E**stimato**R** or simply **Commander** is fast and efficient Gibbs sampling code for joint CMB component separation.

The latest version - **Commander 3** - brings together critical features such as:

- Modern Linear Solver
- Map Making
- Parallelism
- Sky and instrumental modelling

## Prerequisites

To successfully run Commander, you need the following libraries:

- [LibSharp](https://github.com/Libsharp/libsharp)
- [Healpix](https://healpix.sourceforge.io/)
- [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)
- [LAPACK](http://www.netlib.org/lapack/)

## Installation

After you cloned the lastest version of the repo, you need to do the following:

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

## Usage

To run **Commander 1** use:
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
and, finally, run the commander via following command:
```
$ mpirun -n num_proc ~/Commander/commander1/src/commander/commander param_file.txt 2>&1 | tee chains_dir/slurm.txt
```
Here, `num_proc` is the number of processors to use, `slurm.txt` is the file to store output logs.

To run **Commander 2 & 3** the above steps are similar, but you need to use appropriate parameter file. The last command is also modified:
```
$ mpirun -np num_proc ~/Commander/src/commander/commander param_file.txt 2>&1 | tee chains_dir/slurm.txt
```

## License


## Contact


## Citation

- Jewell et al. 2004, ApJ, 609, 1                            
- Wandelt et al. 2004, Phys. Rev. D, 70, 083511
- Eriksen et al. 2004, ApJS, 155, 227 (Commander)
- Eriksen et al. 2008, ApJ, 676, 10  (Joint FG + CMB) 

