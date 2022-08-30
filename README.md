<a name="top"></a>
<p align="center">
    <img src="https://github.com/hke/Commander/blob/master/logo/Commander-logo-large-1024x335.png" height="150">
</p>

**Commander** is an **O**ptimal **M**onte-carlo **M**arkov ch**A**i**N** **D**riven **E**stimato**R** which implements fast and efficient end-to-end CMB posterior exploration through Gibbs sampling.

---

| [Main features](#main-features) | [Quickstart](#quickstart) | [Projects](#projects) | [License](#license) | [Funding](#funding) | [Citation](#citation) |

---

## Main Features

The latest version - `Commander3` - brings together critical features such as:

- Modern Linear Solver
- Map Making
- Parallelism
- Sky and instrumental modelling
- CMB Component Separation

`Commander3` is written using modern `Fortran` standards such as modules, sub modules, and object oriented derived types. The code is highly tuned and optimized to run on High Performance Computing (HPC) facilities, but it can also be run on your local machine.

The previous incarnation of **Commander**, - `Commander2` - is now an internal part of `Commander3`, while the first version of the code, - `Commander1` - is used mainly for debugging and/or legacy purposes. However, `Commander1` has not been officially released; thus, it doesn't support [CMake](https://cmake.org/) installation, as described in [official documentation](https://docs.beyondplanck.science/#/parameters/intro).

---

## Quickstart

Assuming you have installed all 
[prerequisites](https://cosmoglobe.github.io/Commander/#/01_user_manual/prerequisites/README),
you can run one of the following set of commands to install Commander:
<details>
<summary>
<b>Using Intel Compilers</b>
</summary>
<pre><code>
&#36; git clone https://github.com/Cosmoglobe/Commander.git && cd Commander 
&#36; mkdir build && cd build 
&#36; cmake -DCMAKE_INSTALL_PREFIX=&#36;HOME/.local/commander -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DMPI_C_COMPILER=mpiicc -DMPI_CXX_COMPILER=mpiicpc -DMPI_Fortran_COMPILER=mpiifort ..
&#36; cmake --build . --target install -j N  
</code></pre>
where <code>N</code> is the number of processors to use.
</details>
</br>

<details>
<summary>
<b>Using GNU Compilers</b>
</summary>
<pre><code>
&#36; git clone https://github.com/Cosmoglobe/Commander.git && cd Commander 
&#36; mkdir build && cd build 
&#36; cmake -DCMAKE_INSTALL_PREFIX=&#36;HOME/.local/commander -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpic++ -DMPI_Fortran_COMPILER=mpifort ..
&#36; cmake --build . --target install -j N  
</code></pre>
where <code>N</code> is the number of processors to use.
</details>

<p align="justify">
This step will configure Commander installation inside <code>build</code> directory using 
<a href="https://cmake.org/">CMake</a>. Once the configuration is done, the compilation will 
start and it may take some time to finish depending on your system. The resulting 
binary, <code>commander3</code>, will be stored inside 
<code>&#36;HOME/.local/commander/bin</code> and it can be run like any other MPI 
application using <code>mpirun</code> command.  

In case something went wrong or for the complete installation guide 
please refer to the 
<a href="https://cosmoglobe.github.io/Commander/#/">official documentation</a>, 
where we discuss in detail "what is going on?" and "how does it work?", as 
well as providing information on how to compile and run <code>Commander</code> 
on different platforms, including HPCs such as NERSC, UNINETT Sigma2, OWLs etc. 

We have also put up the 
<a href="https://cosmoglobe.github.io/Commander/#/04_faq/README">FAQ</a> 
section for Troubleshooting.
</p>

---

## Projects

Commander framework is part of the following projects:

<p align="center">
    <img src="./logo/Planck_logo.png" height="100"> 
    <img src="./logo/beyondplanck_logo.png" height="100"> 
    <img src="./logo/LiteBIRD-logo-posi-RGB.png" height="100"> 
    <img src="./logo/Cosmoglobe-logo-vertical-large.png" height="100"> 
</p>

---

## Funding

This work has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreements No 776282 (COMPET-4; BeyondPlanck), 772253 (ERC; bits2cosmology) and 819478 (ERC; Cosmoglobe).

<p align="center">
    <img src="./logo/LOGO_ERC-FLAG_EU_.jpg" height="200">
    <img src="./logo/horizon2020_logo.jpg" height="200">
</p>

---

## License

[GNU GPLv3](https://github.com/Cosmoglobe/Commander/blob/master/COPYING)

---

## Citation

If used for published results, please [cite these papers](https://github.com/Cosmoglobe/Commander/blob/master/docs/commander.bib):

- [Jewell et al. 2004, ApJ, 609, 1](https://ui.adsabs.harvard.edu/abs/2004ApJ...609....1J)
- [Wandelt et al. 2004, Phys. Rev. D, 70, 083511](https://ui.adsabs.harvard.edu/abs/2004PhRvD..70h3511W)
- [Eriksen et al. 2004, ApJS, 155, 227 (Commander)](https://ui.adsabs.harvard.edu/abs/2004ApJS..155..227E)
- [Eriksen et al. 2008, ApJ, 676, 10  (Joint FG + CMB)](https://ui.adsabs.harvard.edu/abs/2008ApJ...676...10E)

