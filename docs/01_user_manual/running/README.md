# Running Commander

To run [Commander3](https://github.com/Cosmoglobe/Commander) you need to create the `chains_directory`, which is specified in `parameter_file.txt` via
```
$ mkdir chains_directory
```
You will also need to export `HEALPix` as:
```
$ export HEALPIX=/path/to/healpix/root/directory
```
Lastly, as any MPI application, [Commander3](https://github.com/Cosmoglobe/Commander) can be run by the set of simple commands:
```
$ export OMP_NUM_THREADS=1
$ mpirun -np num_proc /path/to/commander/binary/folder/commander3 param_file.txt 2>&1 | tee chains_dir/slurm.txt
```
Here, `num_proc` is the number of processors to use, `slurm.txt` is the file to store output logs. Note: in case you have used `Makefile` and compiled commander and its dependencies from scratch, the binary is called simple `commander` and it is located inside `/path/to/cloned/repo/Commander/commander3/src`.

You can, of course, wrap this commands into bash script, instead of polluting your `.bashrc` file. Here is an example on how to do it:
```
#!/usr/bin/bash

HEALPIX_ROOT=/path/to/healpix/root
COMM3_BIN=/path/to/commander3/binary
PARAM_FILE=param_BP_v8.00_full.txt
CHAINS_DIR=chains
OUT_FILE=slurm.txt
NUM_PROC=64

export HEALPIX=$HEALPIX_ROOT
export OMP_NUM_THREADS=1
mpirun -np $NUM_PROC $COMM3_BIN $PARAM_FILE 2>&1 | tee $CHAINS_DIR/$OUT_FILE
```
And run it as simply as:
```
bash <script_name>.sh
```

