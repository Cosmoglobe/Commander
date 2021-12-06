# Docker Image

We provide a precompiled Ubuntu Docker image that has an installed
version of Commander.

This allows you to run Commander through a one line command only.

## Prerequisites

1. You will have to have [Docker](https://www.docker.com/) installed.
   Docker provides [detailed installation
   instructions](https://docs.docker.com/get-docker/) on how to install
   it to your system.

2. You will have to have the required Commander input files already
   downloaded to a folder of your choice.


## Instructions

Download the required input files:

We assume that you have
[configured](/01_quick_start/downloader/index?id=configuration) the `bp`
tool to download input files at `/mnt/bp`, but substitute this with a
folder of your choice.

`bp download l2data,auxcmd3`



The download tool then, will download everything inside the
`/mnt/bp/input` folder.

Inside that same folder (`/mnt/bp/input`) copy a commander configuration
file.  We suggest that you use a sample config file from:

`https://gitlab.com/BeyondPlanck/r13y-helper/-/blob/master/sample_config/sample_cm3_config.txt`

The config file has two parameters that will need to be edited to the
location of your download files (`/mnt/bp/input`, in our particular
example).  These paramaters are:

* `INIT_CHAIN01` and
* `FFTW3_MAGIC_NUMBERS`

You can also change the location of the output files (`OUTPUT_DIRECTORY`
parameter).

Do your edits and then save that file as `/mnt/bp/input/sample_config.txt`

Now create an output folder too.  For consistency, create a folder
similar to the `input` folder, as `/mnt/bp/output`.

Finally you can then run commander with:

```
docker run -it \
  -v /mnt/bp/input:/mnt/bp/input \
  -v /mnt/bp/output:/output \
  registry.gitlab.com/beyondplanck/r13y-helper/cm3 \
  commander3 /mnt/bp/input/some_config.txt
```

The first time you run this docker command, the Commender Docker image
will be downloaded, which is quite big (aprox 990MB).  Concequent
executions will use the cached Docke image, and they will be able to
startup almost instantaneously.
