Setup
The goal of this setup is creating an environment where we can run the jupiter notebook and the code in it. To achieve this you first of all need a working Python 3 environment.
Additionally there is a brief introduction to jupyter notebooks, for those wanting to run the guide and are new to notebooks. Many of you probably have these packages already installed, and can thus skip this setup.

Configuring the environment
Setting up our environment means installing Python and all the packages we will be needing for this project. In this guide we are going to use conda as an environment manager and pip as a package manager. There is, however, a wide variety of options out there, but as long as you are able to run the sanity check you should be good.

Installing conda
Conda is a package, dependency and environment manager for several languages, but in this project we will take advantage of the environment management capabilities. We will be using a version called Miniconda, which is installed by downloading and running a bash-script. Note that both the URL and the name of the script depend on your operating system.

macOS and Linux
macOS: https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
Linux: https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ wget <url> 
$ sh Miniconda3-latest-<OS>-x86_64.sh
Running the script will trigger a bunch of prompts, one of which is
Do you wish the installer to initialize Miniconda3 by running conda init? [yes|no]
where we recommend you to answer yes. Once the installer finishes the installation of Miniconda, Python and pip are ready for use.

Windows
The miniconda installer for windows can be downloaded from https://conda.io/en/latest/miniconda.html
Follow the installation guide until it is completed. Open the "Anaconda Prompt" from your start menu. If you are using linux, you should use the terminal to execute the commands, but if you are using windows you should use the "Anaconda prompt".

Creating the environment
We can create an environment with our newly installed conda installation using the command conda create. We do, however, have to source the .bashrc (or .bash_profile for Mac users) file modified in the previous step. This is not necessary for windows:
$ source .bashrc 
Next, we run the command that creates a new environment. Here we name it "ml" and give it the default python version 3.6.
$ conda create --name fg python=3.6
The new environment has to be activated.
$ conda activate fg
If everything went as intended the command line prompt should now be prefixed with the name of the environment.
(fg) $ .

Installing packages
The following packages are needed to run the code in the notebook:
$ pip install matplotlib
$ pip install numpy
$ pip install scipy
$ pip install tqdm
$ pip install corner
 
We also recommend installing Jupyter to be able to run the guide as a notebook:
$ pip install jupyter

Environment sanity check
We can check that everything works as it should by importing the packages in Python:
$ python -c "import scipi" 
$ python -c "import corner»
If you are able to run these commands without anything failing horribly (warnings are OK!) you are all set up.
Note that whether or not this setup runs smoothly depends heavily on what already exists on your OS. T


Jupyter notebooks
Jupyter notebook is a web application for running and sharing code and documentation in a user friendly and readable format. If you installed jupyter as defined here you are all set up to start using notebooks. To start, run jupyter from the terminal (note: You should be in the folder where the notebook is located)
$ jupyter notebook
If you open http://localhost:8888 you should see the file-structure of the folder where you run the command, and if any .ipynb-files (such as the guide) exists, simply click them to get started.
You could also run jupyter notebook from Anaconda Navigator. 
WINDOWS ISSUES
You might experience some problems starting the kernel when you open the notebook. First, be sure that you are running the environment with administrator rights (right click on cmd, or conda-terminal, and click "run as administrator"). If you get this error, follow the answer given. If not, you can try to uninstall the following packages:
* ipykernel
* ipython
* jupyter_client
* jupyter_core
* traitlets
* ipython_genutils
Clean conda's cache by running conda clean -tipsy. Then, install the packages again.
