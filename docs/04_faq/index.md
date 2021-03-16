# Frequently Asked Questions

- *What is Commander?*

Commander is a [Bayesian](https://en.wikipedia.org/wiki/Bayesian_inference) [Markov Chain Monte Carlo](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) sampler for CMB and microwave observations that maps out a posterior distribution that is defined by a user-defined parametric model and data set through [Gibbs sampling](https://en.wikipedia.org/wiki/Gibbs_sampling).

- *Why is Commander written in Fortran, and not in, say, Python?*

The reasons Commander is written in Fortran are partly historical and partly computational. The development of [Commander1](https://arxiv.org/abs/astro-ph/0407028) started in 2004, building directly on the [HEALPix](http://healpix.jpl.nasa.gov) pixelization. HEALPix was at the time primarily written in Fortran, and Fortran therefore also become the language of choice for Commander. Of course, since then Python has become widespread in the cosmological community, and an efficient [HEALPy](https://github.com/healpy/) Python implementation of HEALPix has been implemented. In principle, it should therefore now be possible to write a Python port of Commander. Still, it is important to remember that Commander is a computationally intensive code, with strong requirements in terms of overall efficiency and memory management, and it is likely to be quite challenging to develop a Python-based code that is competitive in terms of computational speed and resources. If somebody is interested in undertaking such a task, however, we would be happy to support the initiative!

- *What is component separation and why should I care?*

Experiments like Planck and WMAP measure the total incoming radiation in a limited number of frequency bands, centered on, for instance, 23, 30, 100, or 857 GHz. This radiation contains contributions from many different physical sources, including CMB, synchrotron, free-free, spinning and thermal dust, CO emission etc. Component separation refers to the process of reconstructing these different astrophysical signals based on the observed frequency maps. This process is essential for modern CMB analysis, both in order to obtain a clean CMB map, which is critically important for cosmology, and because detailed knowledge about the foregrounds is needed to calibrate the instrument properly. In modern CMB analysis pipelines, high-level component separation and low-level data processing have in effect merged into one combined operation.  

- *What is the sky model?*

A sky model is a mathematical model of the incoming radiation. To formulate a good sky model, some knowledge regarding the various physical emission mechanisms is required. For instance, we know that the CMB frequency spectrum very closely follows that of a blackbody, while the synchrotron emission spectrum is known to almost follow a power-law at frequencies above a few GHz. In the Commander framework, we typically write the sky model explicitly in terms of amplitudes that describe the surface brightness of the radiation, and spectral parameters that describe the frequency behaviour of the radiation. For example, the BeyondPlanck sky model reads

![BP sky model](BP_skymodel.png ':size=350px')

(see BeyondPlanck I 2020 for details).

- *What is the difference between Commander1, Commander2 and Commander3?*

Commander1 was the original implementation described by Eriksen et al. (2004,2008), which was the first CMB Gibbs sampler to support joint CMB and foreground estimation. This code was extensively used for the Planck 2013 and 2015 data releases. However, an important shortcoming with Commander1 was a strict requirement that all frequency channels had to have the same angular resolution. In practice, that implied that all channels had to be smoothed to the lowest resolution of any included frequency channel, for instance 1 degree FWHM. Commander2 solved this problem by explicitly deconvolving frequency-dependent beams during the component separation stage, and was as such a true multi-resolution component separation code. Commander3 is essentially the same code as Commander2, but with additional support for time-ordered data processing. As such, Commander is no longer just a component separation code, but rather an end-to-end analysis platform for CMB observations.

- *Are Commander 1 and 2 deprecated?*

Commander1 is still useful for low-resolution analysis, because it is much faster than the later codes, and has more matured spectral index sampling algorithms. Commander2 does not exist as a separate code anymore, as Commander3 is now the same as Commander2, only with more features. In practice, the name "Commander2" is still used to indicate "multi-resolution component separation", though, while "Commander3" is typically used to indicate "global time-ordered CMB analysis".

- *What are the projects which use Commander?*

Commander was developed within the Planck collaboration, and has clearly seen the greatest mileage within that collaboration. However, with the current public software release we believe that the code is now sufficiently mature to be used productively by other teams as well.

- *Are there any Commander alternatives?*

There are certainly many codes that can do subsets of the Commander operations. For instance, DaCapo and MADAM are examples of a TOD calibration and mapmaking codes used by Planck LFI, while NILC, SEVEM and SMICA are examples of component separation codes used by Planck. Plik and CamSpec are examples of CMB likelihood codes. However, Commander is so far unique in merging all these operations into one integrated pipeline.

- *Can I use Commander in one of my projects? If yes, what should I do?*

Absolutely! You should then first read the BeyondPlanck papers (beyondplanck.science) to understand whether it is suitable for you. Then you should download the source codes (http://github.com/cosmoglobe/Commander), and try to compile it. Next, it is highly recommended that you obtain a working test case, for instance the BeyondPlanck data set, and first run that. Then you familiarize yourself with the Commander parameter file, and try to add new data sets or components, and play around with that. Once you are getting experienced, you can start digging into the source code, and try to add your own modules. Be warned, though: Commander is both powerful and complex. It is easy to make mistakes, and end up with non-sensical results. In practice, it is a good idea to collaborate with experienced users, and Cosmoglobe (http://cosmoglobe.uio.no) is a good platform for such collaborations.

- *How do I install Commander on my local machine or on a cluster?*

Check out the automatic cmake installation procedure in the Quick Start guide. The procedure is the same on both a laptop and a cluster. Note that only Linux environments are supported, however, and you will need Fortran and MPI compilers.

- *How do I run Commander?*

Very simple! Copy an existing parameter file into a new working directory. Create an output directory, for instance called "chains". Then type "mpirun -n N path_to_commander/commander param.txt". The trick to producing meaningful results, however, is to actually write the parameter file...  

- *What is the difference between BeyondPlanck pipeline and Commander?*

The BeyondPlanck pipeline consists of three main components, namely 1) pre-processing; 2) Commander; and 3) post-processing.

- *If I want to contribute to the project, what should I do?*

If you have written some new Commander modules or extensions, make a git pull request. This will then be reviewed by the Commander managers, and it may or may not be accepted. If it is, you will be acknowledged as a co-author. If you want to contribute with new data sets or sky models, the best way is to post a message on the Cosmoglobe forum, and initiate a project there.

- *If I have found a bug, what should I do?*

Commander is a research platform, and will as such always remain as a work-in-progress. Do not expect smooth problem-free "black-box" operation. The first thing to check out, however, is the parameter file: Make sure that the inputs make sense. (It is *very* easy to make mistakes, and Commander may or may not warn you about bad choices!) Second, check the documentation and forum. Has the same issue already been reported? If not, send an email to the Cosmoglobe and/or Commander email lists.

- *Is there any place I can post my questions in?*

Yes! Please feel free to use BeyondPlanck/Cosmoglobe forum (forums.beyondplanck.science) actively!

- *Why is it called Commander and not, say, Ruler or Slave?*

When the first Commander code was developed, the most widely used power spectrum estimator was a pseudo-Cl code called MASTER (Hivon et al. 2002). Given that a Gibbs sampler necessarily will be much more computationally expensive than this approach, we anticipated that most users would want to initialize the Gibbs chain on MASTER outputs, and only then start Commander. The combined algorithm would therefore be called "MASTER and Commander" -- which, entirely by chance, also happens to be the title of a brilliant book series by Patrick O'Brian! (Of course, it is also the title of a movie starring Russel Crowe, but that's more forgettable.) Formally speaking, Commander is an acronym for "Commander is an Optimal Monte-carlo Markov chAiN Driven EstimatoR", but, clearly, more work could have been put into that.  
