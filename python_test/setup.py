# TODO: Use numpy.distutils.core.Extension and numpy.distutils.core.setup
# because it includes addtional informatin about f2py
# Standard packaging tool for python
import setuptools
from setuptools import find_packages, setup, Extension
#from numpy.distutils.core import setup, Extension
# OOP way to handle paths in python. Included in Standard Library from python 3.6
from pathlib import Path


#----------------------------------------------------------------------
root_dir = Path(__file__).parent.absolute()
build_dir = root_dir / "build"
download_dir = build_dir / "downloads"
install_dir = build_dir / "install"
#----------------------------------------------------------------------
"""
Description: Override the base extension class for obstructing setuptools from 
building the project for us. All subprojects will be compiled with either python
from source using custom scripts or CMake.
"""
class CMakeExtension(Extension):
    def __init__(self, name, sources=[]):
        super().__init__(name = name, sources = sources)


#----------------------------------------------------------------------
# Parsing README.md file
readme_file = Path(__file__).parent.absolute()/"README.md"
def get_project_description(readme_file):
    with open(readme_file, "r", encoding="utf-8") as readme_file:
        project_description = readme_file.read()
#----------------------------------------------------------------------
# Creating configure options for setuptools to install the project
setup_config_args = {}
# Populating setuptools options
setup_config_args['name'] = "commander3-cmb" 
setup_config_args['description'] = "Optimal Monte-carlo Markov chAiN Driven EstimatoR" 
setup_config_args['long_description'] = get_project_description(readme_file)
setup_config_args['long_description_content_type'] = "text/markdown"
setup_config_args['author'] = "Maksym Brilenkov" 
setup_config_args['author_email'] = "maksym.brilenkov@astro.uio.no"
setup_config_args['license'] = "GPLv3" 
setup_config_args['url'] = "https://github.com/Cosmoglobe/Commander" 
# TODO: make a unified version file, e.g. versions.cmake and parse it here as well
# This is just a placeholder now
setup_config_args['version'] = "1.0.0" 
setup_config_args['python_requires'] = ">=3.6.0"
setup_config_args['install_requires'] = [
        "cmake >= 3.17",
        "requests",
        # NOT gnupg or pretty_bad_protocols
        "python-gnupg",
        "numpy",
        "packaging",
        "psutil"
        ] 
setup_config_args['classifiers'] = [
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Fortran",
        "Programming Language :: C",
        "Programming Language :: C++"
        ] 
#----------------------------------------------------------------------
# Initiating package setup
setup(**setup_config_args)
