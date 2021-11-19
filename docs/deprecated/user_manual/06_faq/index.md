title: FAQ

# Frequently Asked Questions

- *What is Commander?*

Commander is a [Bayesian](https://en.wikipedia.org/wiki/Bayesian_inference) [Markov Chain Monte Carlo](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) sampler for CMB and microwave observations that maps out a posterior distribution that is defined by a user-defined parametric model and data set through [Gibbs sampling](https://en.wikipedia.org/wiki/Gibbs_sampling).

- *Why is Commander written in Fortran, and not in, say, Python?*

The reasons Commander is written in Fortran are partly historical and partly computational. The development of [Commander1](https://arxiv.org/abs/astro-ph/0407028) started in 2004, building directly on the [HEALPix](http://healpix.jpl.nasa.gov) pixelization. HEALPix was at the time primarily written in Fortran, and Fortran therefore also become the language of choice for Commander. Of course, since then Python has become widespread in the cosmological community, and an efficient [HEALPy](https://github.com/healpy/) Python implementation of HEALPix has been implemented. In principle, it should therefore now be possible to write a Python port of Commander. Still, it is important to remember that Commander is a computationally intensive code, with strong requirements in terms of overall efficiency and memory management, and it is likely to be quite challenging to develop a Python-based code that is competitive in terms of computational speed and resources. If somebody is interested in undertaking such a task, however, we would be happy to support the initiative!

- *What is Component Separation and why should I care?*

- *What is the sky model?*

- *What is the difference between Commander 1, Commander 2 and Commander 3?*

- *Are Commander 1 and 2 depricated?*

- *What are the projects which use Commander?*

- *Are there any Commander alternatives?*

- *Can I use Commander in one of my projects? If yes, what should I do?*

- *How do I install Commander on my local machine?*

- *How do I run Commander on my local machine?*

- *How do I install Commander on a cluster?*

- *How do I run Commander on a cluster?*

- *What is the difference between BeyondPlanck pipeline and Commander?*

- *How do I install BeyondPlanck pipeline?*

- *How do I run BeyondPlanck pipeline?*

- *If I want to contribute to the project, what should I do?*

- *If I have found a bug, what should I do?*

- *Is there any place I can post my questions in?*

- *Why is it called Commander and not, say, Ruler or Slave?*
