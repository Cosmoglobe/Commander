# Getting Started

Generally, to get started with Commander3, you will need three things:

* **Compiled version of Commander**
* **A list of input files to work with**
* **An up-to-date version of parameter file**

To get a **compiled version** of Commander you have the following options:

- [**Automated**]: Use CMake to automatically fetch, download and install both Commander and 
  most of its dependencies. This is a **recommended method** and it requires installation of 
  a couple of prerequisites before proceeding, which are described in the next section. 
- [**Manual**]: This method requires compiling all Commander dependencies manually and then, 
  using the predifined configuration file, compile Commander with `Makefile`. Obviously, this 
  approach is more tedious, time-consuming and error-prone, but it serves as a good practice 
  tool to get to know Linux etc.
- [**Precompiled**]: Perhaps, the easiest, but the least tested, approach is to use the 
  precompiled version of Commander inside [Docker image](/01_user_manual/docker/README.md), 
  based on Ubuntu 20.04 server. 

Regarding the **input files**, we provide [downloader
utility](/01_user_manual/downloader/index.md) that will allow you to
easily get access to all required input files for a successfull
Commander execution.

Lastly, the default and reference **parameter files** are shipped together with the Commander 
itself and can be found in `<commander_root>/commander3/parameter_files` folder. The tricky 
part is to identify missing components (if any) and to include them manually. 

The rest of this documentation provides detailed discussions on all these topics above.
