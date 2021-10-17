---
project: Commander
project_github: https://github.com/Cosmoglobe/Commander
project_download: https://github.com/Cosmoglobe/Commander/releases
summary: <p style="text-align: center">
            <img alt="" src="media/commander-logo.png" style="width: 75%; height: 75%">
            <div style="text-align: center; font-size:1.5vw">An Optimal Monte-carlo Markov chAiN Driven EstimatoR which implements fast and efficient end-to-end CMB posterior exploration through Gibbs sampling.</div>
         </p>
author: Cosmoglobe
author_description: Umbrella Project
github: https://github.com/Cosmoglobe
media_dir: ./media
favicon: ./media/commander-favicon.png
src_dir: ../commander3/src
output_dir: ../build/documentation
page_dir: ./user_manual
search: true
graph: true
exclude: comm_mpi_mod.f90
exclude: comm_utils.f90
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

## License

Coming soon...

---

## Projects

Coming soon...

---

## Citation

If used for published results, please cite these papers:

- Jewell et al. 2004, ApJ, 609, 1                            
- Wandelt et al. 2004, Phys. Rev. D, 70, 083511
- Eriksen et al. 2004, ApJS, 155, 227 (Commander)
- Eriksen et al. 2008, ApJ, 676, 10  (Joint FG + CMB)
