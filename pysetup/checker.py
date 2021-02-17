#----------------------------------------------------------------------
import os
# OOP way to handle paths in python. Included in Standard Library from python 3.6
from pathlib import Path
# To count processors basically
import psutil
# to unpack archavies
import shutil
# to run shell commands
import subprocess
# using regular expressions
import re
# to work with versions
from packaging import version
# for working with lists
from collections import defaultdict
#----------------------------------------------------------------------
#from setuptools import Command
#from setuptools.command.build_ext import build_ext
#from setuptools.command.install import install
#from distutils.command.build import build
#----------------------------------------------------------------------
# Custom modules
from .downloader import Downloader
#----------------------------------------------------------------------

# TODO: 
# [ ] Make check_mpi() work, there is probalem with Intel, it has different versions 

"""
Description: This class contains methods to check for compilers and other necessary prerequisistes

"""
class Checker(object):

    """
    Description: Method to compile sources with the given set of commands
    """
    @staticmethod
    def compile_source(details, env, download_dir, install_dir, num_procs):
        commands = details['command'].format(install_dir=install_dir, num_procs=num_procs).split("&&")
        
        project_archive = Path(Path(details["source"]).stem + Path(details["source"]).suffix) 
        project_dir = download_dir / Path(str(project_archive).replace(".tar.gz",""))
        for command in commands:
            # Spolitting string into a list and removing empty strings
            command = list(filter(None, command.split(" ")))
            #result = subprocess.run(["make", "--version"], 
            #        shell=False, check=True, env=env, cwd=project_dir)#,
            print(command)
            if ((command[0] == "mkdir") and (not Path.is_dir(project_dir/command[1]))):
                Path.mkdir(project_dir/command[1])
            elif(command[0] == "cd"):
                project_dir = project_dir/command[1]
            if (command[0] != "mkdir") and (command[0] != "cd"):
                result = subprocess.run(command, 
                        shell=False, check=True, env=env, cwd=project_dir)#,
                #        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            print(f"--------------------------------------------")

    """
    Description: Method to check for available compilers and their version(s).
    Commander3 uses Fortran 2003 and 2008 features and so we need appropriate 
    compilers to support this. So far these compilers proved to be working
    correctly:
        Intel: icpc, icc, ifort, version 18+;
        GNU: gcc, g++, gfortran, version 9.x+;
    """
    @staticmethod
    def check_compilers(version_string, env):
        # TODO: perform a check on Intel and Intel MPI, and if yes, then no need to download 
        # and install GCC. Comething like:
        # If intel == true:
        compilers = {"Intel": ["icc", "icpc", "ifort"], "GNU": ["gcc", "g++", "gfortran"]} 
        min_versions = {"Intel": "18.0.0", "GNU": "9.3.0"}
        compilers_found = defaultdict(list) 
        # Checking whether we have compilers on the system
        for vendor, executables in compilers.items():
            # Creating regex object to compare versions
            regex = re.compile(version_string)
            for executable in executables:
                # Getting project's location
                exe_location = shutil.which(executable)
                # If present -- check version, else -- say so
                if exe_location is not None:
                    # Getting project details
                    result = subprocess.run([f"{executable}","--version"], 
                            shell=False, check=True, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                    exe_version = regex.search(result.stdout).group()
                    print(f"-- Found -- {executable}, v{exe_version} -- in:\n{exe_location}")
                    if version.parse(exe_version) >= version.parse(min_versions[vendor]):
                        print("-- Current version is suitable for building the project.")
                        compilers_found[vendor].append(True) 
                    else:
                        print("-- Current version is not enough for building the project.")
                        compilers_found[vendor].append(False)
                else:
                    print(f"-- Couldn't find an executable -- {executable}")
                    compilers_found[vendor].append(False)
        #
        compilers_to_use = [] 
        for vendor in compilers.keys():
            if all(compilers_found[vendor]):
                print(f"-- Found suitable compilers from {vendor}")
                compilers_to_use = compilers[vendor]
                print(compilers_to_use)
                compile_compilers = False
                break
            else:
                compile_compilers = True
                compilers_to_use = ["gcc", "g++", "gfortran"]

        # Only if any other compilers were not found, we check for presence 
        # of system shipped compilers -- we need to have some c/c++ compilers
        # to compile compilers
        if compile_compilers:
            system_compilers = ['cc', 'c++']
            for compiler in system_compilers:
                compiler_location = shutil.which(compiler)
                if compiler_location is not None:
                    result = subprocess.run([f"{compiler}","--version"], 
                            shell=False, check=True, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                    print(f"-- Found -- {compiler} -- in:\n{compiler_location}")
                else:
                    raise Exception("Couldn't identify {compiler} compiler!\nPlease, install at least C & C++ compilers!")

        return compile_compilers, compilers_to_use

    """
    Description: Method checks for existence of MPI on the system from different vendors.
    """
    @staticmethod
    def check_mpi(version_string, env):

        # TODO: figure out problem with intel -- basically 
        # it has different version string, so we need to use anotehr version string,
        # OR implement MPI compilation whithin cmake and thus no need to do anything here.
        compile_mpi = False 
        #mpi_to_use = ["mpicc", "mpic++", "mpifort"]
        mpi_to_use = ["mpiicc", "mpiicpc", "mpiifort"]
        return compile_mpi, mpi_to_use

    """
    Description: Method checks all available projects. I starts with the compilers and then goes
    for other prerequisites. If NO COMPILERS WERE FOUND than it will give you an error! There should
    be something on the system beforehand.
    """
    @staticmethod
    def check_prerequisites(prerequisites, env, download_dir, install_dir, num_procs, signatures):
        # Creating regex pattern to match versions
        version_string = '(?:(\d+\.[.\d]*\d+))'
        # Checking if compilers are available on the system
        compile_compilers, compilers_to_use = Checker.check_compilers(version_string, env)
        # Checking if mpi wrappers are available on the system
        compile_mpi, mpi_to_use = Checker.check_mpi(version_string, env)
        # Identify which projects needs to be skipped and which should be build from source:
        # True -- use existing project
        # False -- build from source (e.g. make with Makefile)
        # None -- build from source (e.g. make with build.sh script)
        build_from_src = {}
        prerequisites_to_install = prerequisites
        # if we found compilers, we skip compilation step for GNU GCC
        if not compile_compilers: 
            del prerequisites_to_install['compilers']
        if not compile_mpi:
            del prerequisites_to_install['mpi']
        # 
        #for proj, details in prerequisites.items():
        for proj, details in prerequisites_to_install.items():
            print(f"-- ============================================")
            print(f"-- Start working with project - {proj}...")
            executables = details['executables']
            # Checking if the key is present, we populate its value
            # and if it is not, we just say that any version is enough
            if 'min_version' in details.keys():
                min_version = details['min_version']
            else:
                min_version = ""
            # Creating regex object to compare versions
            regex = re.compile(version_string)
            # Looping through each executable and compare versions
            for executable in executables:
                # Now, performing several checks for the compilers
                # Getting project's location
                exe_location = shutil.which(executable)
                # If present -- check version, else -- say so
                if exe_location is not None:
                    # Getting project details
                    result = subprocess.run([f"{executable}","--version"], 
                            shell=False, check=True, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                    # This just returns None if it is not found
                    print(regex.search(result.stdout))
                    exe_version = regex.search(result.stdout).group()
                    print(f"-- Found -- {executable}, v{exe_version} -- in:\n{exe_location}")

                    if version.parse(exe_version) >= version.parse(min_version):
                        print("-- Current version is suitable for building the project.")
                        build_from_src[proj] = True 
                    else:
                        print("-- Current version is not enough for building the project.")
                        build_from_src[proj] = False 
                else:
                    print(f"-- Couldn't find an executable -- {executable}")
                    build_from_src[proj] = False

            # Downloading and compiling stuff from source
            # TODO: tune GCC installation for Mac?
            if build_from_src[proj] == False:
                Downloader.get_source(details, env, download_dir, signatures)
                Checker.compile_source(details, env, download_dir, install_dir, num_procs)
                print(f"-- Installed project is:")
                subprocess.run([f"{proj}","--version"], shell=False, check=True, env=env)
                print(f"--------------------------------------------")

        return compilers_to_use, mpi_to_use

