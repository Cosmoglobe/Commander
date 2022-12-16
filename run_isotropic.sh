#!/bin/bash
list1=owl{17..24}.uio.no
list2=owl{25..28}.uio.no
list3=owl{29..30}.uio.no
list4=owl{31..35}.uio.no
list5=owl{36..37}.uio.no
lists=(list1 list2 list3 list4 list5)
ns=(24 72 64 64 128)
tLen=${#lists[@]}
builds=("owl1724" "owl2528" "owl2930" "owl3135" "owl3637")
build=0
n=36
for (( i=0; i<${tLen}; i++ ));
do
  c=${!lists[$i]}
  v=$(eval "echo $c")
  if [[ "$v" == *"$HOSTNAME"* ]]
  then
    build=$((i+1))
    build=${builds[$i]}
    n=${ns[$i]}
  fi
done
echo $build

export LD_LIBRARY_PATH=/mn/stornext/d16/cmbco/bp/johanres/commander_camb/Commander/build/install/lib:/mn/stornext/d16/cmbco/bp/johanres/commander_camb/Commander/build/install/healpix/lib:$LD_LIBRARY_PATH

killall -9 commander3
#COMMANDER_PARAMS_DEFAULT=$HOME"/Commander/commander3/parameter_files/defaults/"
#pfile=param_WMAP_all_create_radio.txt
#mpirun -env I_MPI_FABRICS shm -n $n /mn/stornext/u3/duncanwa/Commander/build_$build/install/bin/commander3 $pfile --OUTPUT_DIRECTORY=$dir 2>&1| tee $dir/slurm.txt

pfile=/mn/stornext/d16/cmbco/bp/johanres/commander_camb/Commander/commander3/parameter_files/param_camb_sim_isotropic_nomask.txt
dir=chains_camb_nside1024_isotropic_nomask

#pfile=param_WMAP_bp_comp_only.txt
#dir=chains_WMAP_amps_220126
mkdir -p $dir
mpiexec -env I_MPI_FABRICS shm -n $n /mn/stornext/d16/cmbco/bp/johanres/commander_camb/Commander/commander3/src/commander $pfile --OUTPUT_DIRECTORY=$dir 2>&1| tee $dir/slurm.txt
