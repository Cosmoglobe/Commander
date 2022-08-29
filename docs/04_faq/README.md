# Frequently Asked Questions

## About Commander

<details>
<summary>
What is Commander?
</summary>
<p align="justify">
Commander is a <a href="https://en.wikipedia.org/wiki/Bayesian_inference">Bayesian</a> 
<a href="https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo">Markov Chain Monte Carlo</a> 
sampler for CMB and microwave observations that maps out a posterior distribution that is 
defined by a user-defined parametric model and data set through 
<a href="https://en.wikipedia.org/wiki/Gibbs_sampling">Gibbs sampling</a>.
<p>
</details>
</br>

<details>
<summary>
Why is Commander written in Fortran, and not in, say, Python?
</summary>
<p align="justify">
The reasons Commander is written in Fortran are partly historical and partly computational. 
The development of <a href="https://arxiv.org/abs/astro-ph/0407028">Commander1</a> started in 
2004, building directly on the <a href="http://healpix.jpl.nasa.gov">HEALPix</a> pixelization. 
HEALPix was at the time primarily written in Fortran, and Fortran therefore also become the 
language of choice for Commander. Of course, since then Python has become widespread in the 
cosmological community, and an efficient <a href="https://github.com/healpy/">HEALPy</a> 
(Python implementation of HEALPix) has been implemented. In principle, it should therefore 
now be possible to write a Python port of Commander. Still, it is important to remember that 
Commander is a computationally intensive code, with strong requirements in terms of overall 
efficiency and memory management, and it is likely to be quite challenging to develop a 
Python-based code that is competitive in terms of computational speed and resources. If 
somebody is interested in undertaking such a task, however, we would be happy to support the 
initiative!
</p>
</details>
</br>

<details>
<summary>
What is component separation and why should I care?
</summary>
<p align="justify">
Experiments like Planck and WMAP measure the total incoming radiation in a limited number of 
frequency bands, centered on, for instance, 23, 30, 100, or 857 GHz. This radiation contains 
contributions from many different physical sources, including CMB, synchrotron, free-free, 
spinning and thermal dust, CO emission etc. Component separation refers to the process of 
reconstructing these different astrophysical signals based on the observed frequency maps. 
This process is essential for modern CMB analysis, both in order to obtain a clean CMB map, 
which is critically important for cosmology, and because detailed knowledge about the 
foregrounds is needed to calibrate the instrument properly. In modern CMB analysis pipelines,
high-level component separation and low-level data processing have in effect merged into one 
combined operation.  
</p>
</details>
</br>

<details>
<summary>
What is the sky model?
</summary>
<p align="justify">
A sky model is a mathematical model of the incoming radiation. To formulate a good sky model, 
some knowledge regarding the various physical emission mechanisms is required. For instance, 
we know that the CMB frequency spectrum very closely follows that of a blackbody, while the 
synchrotron emission spectrum is known to almost follow a power-law at frequencies above a 
few GHz. In the Commander framework, we typically write the sky model explicitly in terms of 
amplitudes that describe the surface brightness of the radiation, and spectral parameters that 
describe the frequency behaviour of the radiation. For example, the BeyondPlanck sky model 
reads

![BP sky model](BP_skymodel.png ':size=350px')

(see BeyondPlanck I 2020 for details).
</p>
</details>
</br>

<details>
<summary>
What is the difference between Commander1, Commander2 and Commander3?
</summary>
<p align="justify">
Commander1 was the original implementation described by Eriksen et al. (2004,2008), which 
was the first CMB Gibbs sampler to support joint CMB and foreground estimation. This code 
was extensively used for the Planck 2013 and 2015 data releases. However, an important 
shortcoming with Commander1 was a strict requirement that all frequency channels had to 
have the same angular resolution. In practice, that implied that all channels had to be 
smoothed to the lowest resolution of any included frequency channel, for instance 1 degree FWHM. 
Commander2 solved this problem by explicitly deconvolving frequency-dependent beams during the 
component separation stage, and was as such a true multi-resolution component separation code. 
Commander3 is essentially the same code as Commander2, but with additional support for 
time-ordered data processing. As such, Commander is no longer just a component separation 
code, but rather an end-to-end analysis platform for CMB observations.
</p>
</details>
</br>

<details>
<summary>
Are Commander 1 and 2 deprecated?
</summary>
<p align="justify">
Commander1 is still useful for low-resolution analysis, because it is much faster than the 
later codes, and has more matured spectral index sampling algorithms. Commander2 does not 
exist as a separate code anymore, as Commander3 is now the same as Commander2, only with more 
features. In practice, the name "Commander2" is still used to indicate "multi-resolution 
component separation", though, while "Commander3" is typically used to indicate "global 
time-ordered CMB analysis".
</p>
</details>
</br>

<details>
<summary>
What are the projects which use Commander?
</summary>
<p align="justify">
Commander was developed within the Planck collaboration, and has clearly seen the greatest 
mileage within that collaboration. However, with the current public software release we 
believe that the code is now sufficiently mature to be used productively by other teams as well.
</p>
</details>
</br>

<details>
<summary>
Are there any Commander alternatives?
</summary>
<p align="justify">
There are certainly many codes that can do subsets of the Commander operations. For instance, 
DaCapo and MADAM are examples of a TOD calibration and mapmaking codes used by Planck LFI, 
while NILC, SEVEM and SMICA are examples of component separation codes used by Planck. Plik 
and CamSpec are examples of CMB likelihood codes. However, Commander is so far unique in 
merging all these operations into one integrated pipeline.
</p>
</details>
</br>

<details>
<summary>
Can I use Commander in one of my projects? If yes, what should I do?
</summary>
<p align="justify">
Absolutely! You should then first read the BeyondPlanck papers (beyondplanck.science) to 
understand whether it is suitable for you. Then you should download the source codes 
(http://github.com/cosmoglobe/Commander), and try to compile it. Next, it is highly 
recommended that you obtain a working test case, for instance the BeyondPlanck data set, 
and first run that. Then you familiarize yourself with the Commander parameter file, and 
try to add new data sets or components, and play around with that. Once you are getting 
experienced, you can start digging into the source code, and try to add your own modules. 
Be warned, though: Commander is both powerful and complex. It is easy to make mistakes, 
and end up with non-sensical results. In practice, it is a good idea to collaborate with 
experienced users, and Cosmoglobe (http://cosmoglobe.uio.no) is a good platform for such 
collaborations.
</p>
</details>
</br>

<details>
<summary>
How do I install Commander on my local machine or on a cluster?
</summary>
<p align="justify">
Check out the automatic cmake installation procedure in the Quick Start guide. The procedure 
is the same on both a laptop and a cluster. Note that only Linux environments are supported, 
however, and you will need Fortran and MPI compilers.
</p>
</details>
</br>

<details>
<summary>
How do I run Commander?
</summary>
<p align="justify">
We have a separate section on how to run Commander. But, the basic steps are as follows:
<ol>
  <li>Copy an existing parameter file into a new working directory</li>
  <li>Create an output directory, for instance called <code>chains</code></li>
  <li>Run it is as any MPI application:
      <pre><code>$ mpirun -n N <path_to_commander>/commander3 <parameter_file> </code></pre>
      where <code>N</code> is the number of processors to use.
  </li>
</ol>
</p>
</details>
</br>

<details>
<summary>
What is the difference between BeyondPlanck pipeline and Commander?
</summary>
<p align="justify">
The BeyondPlanck pipeline consists of three main components, namely: 
<ol>
  <li>Pre-processing</li>
  <li>Commander</li>
  <li>Post-processing</li>
</ol>
</p>
</details>
</br>

<details>
<summary>
If I want to contribute to the project, what should I do?
</summary>
<p align="justify">
If you have written some new Commander modules or extensions, make a git pull request. 
This will then be reviewed by the Commander managers, and it may or may not be accepted. 
If it is, you will be acknowledged as a co-author. If you want to contribute with new data 
sets or sky models, the best way is to post a message on the Cosmoglobe forum, and initiate 
a project there.
</p>
</details>
</br>

<details>
<summary>
If I have found a bug, what should I do?
</summary>
<p align="justify">
Commander is a research platform, and will as such always remain as a work-in-progress. 
Do not expect smooth problem-free "black-box" operation. The first thing to check out, 
however, is the parameter file: Make sure that the inputs make sense. (It is <em>very</em> easy 
to make mistakes, and Commander may or may not warn you about bad choices!) Second, check 
the documentation and forum. Has the same issue already been reported? If not, send an email 
to the Cosmoglobe and/or Commander email lists.
</p>
</details>
</br>

<details>
<summary>
Is there any place I can post my questions in?
</summary>
<p align="justify">
Yes! Please feel free to use BeyondPlanck/Cosmoglobe forum (forums.beyondplanck.science) actively!
</p>
</details>
</br>

<details>
<summary>
Why is it called Commander and not, say, Ruler or Slave?
</summary>
<p align="justify">
When the first Commander code was developed, the most widely used power spectrum estimator was
a pseudo-Cl code called MASTER (Hivon et al. 2002). Given that a Gibbs sampler necessarily will 
be much more computationally expensive than this approach, we anticipated that most users would 
want to initialize the Gibbs chain on MASTER outputs, and only then start Commander. The 
combined algorithm would therefore be called "MASTER and Commander" -- which, entirely by 
chance, also happens to be the title of a brilliant book series by Patrick O'Brian! (Of course, 
it is also the title of a movie starring Russel Crowe, but that's more forgettable.) Formally 
speaking, Commander is an acronym for "Commander is an Optimal Monte-carlo Markov chAiN Driven 
EstimatoR", but, clearly, more work could have been put into that.  
</p>
</details>
</br>

## Troubleshooting

### Compilation 

#### HEALPix

<details>
<summary>
<b>Error</b>: 
<code>
/bin/rm: cannot remove ‘/tmp/Healpix_autolist.txt’: Operation not permitted
</code>
</summary>
<p align="justify">
The full error code may read:
<pre><code>
/bin/rm: cannot remove ‘/tmp/Healpix_autolist.txt’: Operation not permitted
touch: cannot touch ‘/tmp/Healpix_autolist.txt’: Permission denied
./hpxconfig_functions.sh: line 237: /tmp/Healpix_autolist.txt: Permission denied
./hpxconfig_functions.sh: line 237: /tmp/Healpix_autolist.txt: Permission denied
./hpxconfig_functions.sh: line 237: /tmp/Healpix_autolist.txt: Permission denied
./hpxconfig_functions.sh: line 237: /tmp/Healpix_autolist.txt: Permission denied
./configure: line 173: /tmp/Healpix_autolist.txt: Permission denied
</code></pre>
<b>Solution</b>:

Remove the <code>/tmp/Healpix_autolist.txt</code> and compile again.

<b>Reason</b>: 

Permission issues. HEALPix creates severeal different files during its installation. One of 
them is the above <code>Healpix_autolist.txt</code>. Is say, some user installs HEALPix from 
scratch, then this file will belong to that user, meaning that cmake won't be able to 
rewrite it during automatic installation process.
</p>
</details>
</br>

### Running Commander

<details>
<summary>
<b>Error</b>: 
<code>
forrtl: severe (24): end-of-file during read, unit 10, file 
some/path/here/chains_pw/regdef_c01_p01.unf
</code>
</summary>
<p align="justify">
The full error code may read:
<pre><code>
forrtl: severe (24): end-of-file during read, unit 10, file 
some/path/here/chains_pw/regdef_c01_p01.unf
Image              PC                Routine            Line        Source             
libifcoremt.so.5   00002B78AB782622  for__io_return        Unknown  Unknown
libifcoremt.so.5   00002B78AB7BDB59  for_read_seq          Unknown  Unknown
commander          000000000047F6C5  comm_fg_component        1993  comm_fg_component_mod.f90
commander          0000000000493DA4  comm_fg_component         789  comm_fg_component_mod.f90
commander          00000000004F8437  comm_fg_mod_mp_in          73  comm_fg_mod.f90
commander          0000000000408D7B  MAIN__                    208  commander.f90
commander          0000000000407D2E  Unknown               Unknown  Unknown
libc-2.17.so       00002B78ADF683D5  __libc_start_main     Unknown  Unknown
commander          0000000000407C35  Unknown               Unknown  Unknown
</code></pre>

<b>Solution</b>:

Remove the whole chain directory (all the files in it) and start anew.

</p>
</details>
</br>



<details>
<summary>
<b>Error</b>: 
<code>
forrtl: severe (168): Program Exception - illegal instruction
</code>
</summary>
<p align="justify">
The full error code may read:
<pre><code>
forrtl: severe (168): Program Exception - illegal instruction
Image              PC                Routine            Line        Source             
commander3         00000000010C861B  Unknown               Unknown  Unknown
libpthread-2.17.s  00007F43B9FFB630  Unknown               Unknown  Unknown
commander3         00000000008F5E11  Unknown               Unknown  Unknown
commander3         00000000007F22F2  sharp_mp_sharp_ma         133  sharp.f90
commander3         00000000007818E7  comm_map_mod_mp_c         288  comm_map_mod.f90
commander3         000000000080A6F8  comm_data_mod_mp_         123  comm_data_mod.f90
commander3         00000000004812D3  MAIN__                    154  commander.f90
commander3         0000000000480862  Unknown               Unknown  Unknown
libc-2.17.so       00007F43B7142555  __libc_start_main     Unknown  Unknown
commander3         0000000000480769  Unknown               Unknown  Unknown
</code></pre>

<b>Solution</b>:

[TODO]: Figure this out? Need to recompile Commaner for the same CPU type 

<b>Reason</b>:

Compiled Commander for one CPu ad tried to run it on the other (e.g. Intel and AMD or 
different Intel versions)

</p>
</details>
</br>
