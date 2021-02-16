import os
# OOP way to handle paths in python. Included in Standard Library from python 3.6
from pathlib import Path
# parsing command line arguments
import argparse
# to read json files
import json
# to unpack archavies
import shutil
# to run shell commands
import subprocess
# using regular expressions
import re
# to work with versions
from packaging import version
#import multiprocessing
import psutil
# to get data from internet
import requests
# verification of tar.gz.sig files
import gnupg
# find hash values of a file
import hashlib
# for working with lists
from collections import defaultdict
# TODO: Use numpy.distutils.core.Extension and numpy.distutils.core.setup
# because it includes addtional informatin about f2py
# Standard packaging tool for python
from setuptools import find_packages, setup, Extension
from setuptools.command.build_ext import build_ext
# Custom modules

#----------------------------------------------------------------------
"""
Description: This method calculates MD5 hash sum of an archive and compares it to
the given by the authors.
"""
# Method varifies provided MD5 hash and raises an exception if the differ
def verify_md5(archive_md5_hash, file_to_check, source_file):
    print(f"--------------------------------------------")
    md5_hash = hashlib.md5()
    with open(file_to_check, "rb") as current_file:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: current_file.read(4096), b""):
            md5_hash.update(byte_block)
        print(f"Calculated md5 hash is {md5_hash.hexdigest()}")
        print(f"Provided md5 hash is {archive_md5_hash}")
    if (archive_md5_hash == md5_hash.hexdigest()):
        print(f"Hashes are the same for {source_file}. Proceeding...")
    else:
        raise Exception(f"Hashes differ! Something went wrong while downloading {source_file}")

"""
Description: This method calculates SHA1 hash sum of an archive and compares it to
the given by the authors.
"""
def verify_sha1(archive_sha1_hash, file_to_check, source_file): 
    print(f"--------------------------------------------")
    sha1_hash = hashlib.sha1()
    with open(file_to_check, "rb") as current_file:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: current_file.read(4096), b""):
            sha1_hash.update(byte_block)
        print(f"Calculated sha1 hash is {sha1_hash.hexdigest()}")
        print(f"Provided sha1 hash is {archive_sha1_hash}")
    if (archive_sha1_hash == sha1_hash.hexdigest()):
        print(f"Hashes are the same for {source_file}. Proceeding...")
    else:
        raise Exception(f"Hashes differ! Something went wrong while downloading {source_file}")

"""
Description: This method calculates SHA256 hash sum of an archive and compares it to
the given by the authors.
"""
def verify_sha256(archive_sha256_hash, file_to_check, source_file): 
    print(f"--------------------------------------------")
    sha256_hash = hashlib.sha256()
    with open(file_to_check, "rb") as current_file:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: current_file.read(4096), b""):
            sha256_hash.update(byte_block)
        print(f"Calculated sha256 hash is {sha256_hash.hexdigest()}")
        print(f"Provided sha256 hash is {archive_sha256_hash}")
    if (archive_sha256_hash == sha256_hash.hexdigest()):
        print(f"Hashes are the same for {source_file}. Proceeding...")
    else:
        raise Exception(f"Hashes differ! Something went wrong while downloading {source_file}")

"""
Description: This method calculates SHA512 hash sum of an archive and compares it to
the given by the authors. Raises an error if hashes do not match.
"""
def verify_sha512(archive_sha512_hash, file_to_check, source_file): 
    print(f"--------------------------------------------")
    sha512_hash = hashlib.sha512()
    with open(file_to_check, "rb") as current_file:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: current_file.read(4096), b""):
            sha512_hash.update(byte_block)
        print(f"Calculated sha512 hash is {sha512_hash.hexdigest()}")
        print(f"Provided sha512 hash is {archive_sha512_hash}")
    if (archive_sha512_hash == sha512_hash.hexdigest()):
        print(f"Hashes are the same for {source_file}. Proceeding...")
    else:
        raise Exception(f"Hashes differ! Something went wrong while downloading {source_file}")

"""
Description: This method verifies archive signature files (.tar.sig). It raises an error
if signature cannot be verified.
"""
def verify_signature(details, download_dir, file_to_check, signatures, env):
    # Converting Path into string values
    signatures = str(signatures)
    file_to_check = str(file_to_check)
    # Performing download of an archive signature file
    archive_signature_url = details["signature"]
    archive_name = Path(Path(details["signature"]).stem + Path(details["signature"]).suffix) 
    archive_signature_file = download_dir / archive_name 
    perform_download(archive_signature_url, archive_name, archive_signature_file)
    print("Verifying archive signature...")
    gpg = gnupg.GPG(keyring = signatures)
    with open(archive_signature_file, "rb") as signature_file:
        verified = gpg.verify_file(signature_file, file_to_check)
    if not verified: 
        raise ValueError("Signature could not be verified!")
    elif(verified.key_status is not None):
        raise ValueError("There is a problem detected with signature -- key status is not None!")
    else:
        print(f"--------------------------------------------")
        print("Signature has been verified. Details:")
        print(f"Fingerprint: {verified.fingerprint}")
        print(f"Signer: {verified.username}")
        print(f"Sign date: {verified.creation_date}")
        print(f"Signature status: {verified.status}")
        print(f"Trust level: {verified.trust_level}")
        print(f"Trust text: {verified.trust_text}")
        print(f"Key status: {verified.key_status}")
#----------------------------------------------------------------------
# TODO: most probably this needs to be merged into setup.py commands.
# Look into:
# https://github.com/diegoferigo/cmake-build-extension/blob/master/src/cmake_build_extension/build_extension.py
# for more details and examples.
"""
Description: This method parses command line arguments to use them in CMake command.
"""
def get_cmd_args(root_dir):
    # Class instantiation
    parser = argparse.ArgumentParser()
    # Populating cmd arguments 
    default_build_dir = str(root_dir.joinpath("build"))
    parser.add_argument('--build-dir', type=str, action='store', default=default_build_dir, 
            help='Directory where to build (sub)project(s).')
    default_download_dir = str(root_dir.joinpath("build","downloads"))
    parser.add_argument('--download-dir', type=str, action='store', default=default_download_dir, 
            help='Directory where to download (sub)project(s).')
    default_install_dir = str(root_dir.joinpath("build", "install"))
    parser.add_argument('--install-dir', type=str, action='store', default=default_install_dir, 
            help='Directory where to install (sub)project(s).')
    parser.add_argument('--num-procs', type=int, action='store', default=psutil.cpu_count(logical=False), 
            help='Number of processes to use.')
    # Compilers
    default_c_compiler = "gcc"
    parser.add_argument('--c-comp', type=str, action='store', default=default_c_compiler, 
            help=f'C compiler to use. Default: {default_c_compiler}.')
    default_cxx_compiler = "g++"
    parser.add_argument('--cxx-comp', type=str, action='store', default=default_cxx_compiler, 
            help=f'C++ compiler to use. Default: {default_cxx_compiler}.')
    default_fortran_compiler = "gfortran"
    parser.add_argument('--fortran-comp', type=str, action='store', default=default_fortran_compiler, 
            help=f'Fortran compiler to use. Default: {default_fortran_compiler}.')
    #
    args = parser.parse_args()

    return args

"""
Description: This method downlaods source archive if it has not been downloaded yet.
"""
def perform_download(source_url, source_file, full_path_to_source):
    #source_file = str(source_file)
    # Downloading file only if it is not present in specified location
    #if not os.path.isfile(full_path_to_source):
    if not Path.is_file(full_path_to_source):
        #-------------------------------------------------
        # Creating class instance
        # if we set stream=True in requests.get(...) then 
        # headers['Transfer-Encoding'] = 'chunked' is set 
        # in the HTTP headers. Thus specifying the Chunked 
        # transfer encoding. In chunked transfer encoding, 
        # the data stream is divided into a series of 
        # non-overlapping "chunks". The chunks are sent out 
        # independently of one another by the server.
        #-------------------------------------------------
        response = requests.get(source_url, stream = True)
        print(f"Getting -- {source_file} -- from:\n{source_url}")
        if response.status_code == 200: 
            with open(full_path_to_source, 'wb') as current_file:
                current_file.write(response.raw.read())
        else:
            raise Exception(f"Something went wrong while downloading {source_file}")
    else:
        print(f"File -- {source_file} -- already exists! Skipping this step...")

"""
Description: This method veerifies downloaded file via various methods. ALL AVAILABLE METHODS
WILL BE CHECKED. IF ANY OF PROVIDED METHODS DOESN"T WORK, THE ERROR WILL BE RAISED AND SCRIPT
STOPS.
"""
def check_source(details, download_dir, source_file, file_to_check, signatures, env):
    possible_checks = ['signature', 'MD5', 'SHA1', 'SHA256', 'SHA512']
    # This variable is only for pretty message. the important one is file_to_check.
    #source_file = details['source'][1] + details['source'][2]
    #source_file  = str(Path(Path(details["source"]).stem + Path(details["source"]).suffix)) 
    file_to_check = str(file_to_check)
    for check in possible_checks:
        value = details.get(check)
        # If key exists, the value is not "None", "False", "''", or "0"
        if (check in details) and value: #and value.strip()): #details.keys():
            # We will check all hashes and signatures available for us (enchanced secuirty? :))
            if (check == 'signature'):
                verify_signature(details, download_dir, file_to_check, signatures, env)
            if (check == 'MD5'):
                verify_md5(details['MD5'], file_to_check, source_file)
            if (check == 'SHA1'):
                verify_sha1(details['SHA1'], file_to_check, source_file)
            if (check == 'SHA256'):
                verify_sha256(details['SHA256'], file_to_check, source_file)
            if (check == 'SHA512'):
                verify_sha512(details['SHA512'], file_to_check, source_file)
    print(f"--------------------------------------------")

"""
Description: Method to unpack archives.
"""
def unpack_archive(project_file, project_dir, download_dir):
    #project_dir = Path(project_dir)
    #print(project_dir)
    #print(type(project_dir))
    print(type(project_file))
    # Unpacks archive unless the file is already present on disk
    if not Path.is_dir(project_dir):
    #if not os.path.isdir(project_dir):
        print(f"Unpacking file...")
        shutil.unpack_archive(str(project_file), str(download_dir))#, project_dir)
    else:
        print(f"Archive's already unpacked! Skipping this step...")

"""
Description: Method to download source archives, check them with provided options
and unpacks them to a specified location.
"""
def get_source(details, env, download_dir, signatures):
    # Performing download of an archive
    archive_url               = details["source"]
    project_archive           = Path(Path(details["source"]).stem + Path(details["source"]).suffix) 
    full_path_to_project_file = download_dir / Path(project_archive) 
    perform_download(archive_url, project_archive, full_path_to_project_file)
    # Perfoming check for an archive corruption using signature file, md5 or other hash.
    check_source(details['checker'], download_dir, project_archive, full_path_to_project_file, signatures, env)
    # Once archive is verified with its digital signature or hashes, we proceed with installation
    # TODO: finish this part -- unpack archives
    project_dir = download_dir / Path(str(project_archive).replace(".tar.gz",""))
    unpack_archive(full_path_to_project_file, project_dir, download_dir)

"""
Description: Method to compile sources with the given set of commands
"""
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
    GNU: gcc, g++, gfortran, version 9.x;
"""
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
                print(f"Found -- {executable}, v{exe_version} -- in:\n{exe_location}")
                if version.parse(exe_version) >= version.parse(min_versions[vendor]):
                    print("Current version is suitable for building the project.")
                    compilers_found[vendor].append(True) 
                else:
                    print("Current version is not enough for building the project.")
                    compilers_found[vendor].append(False)
            else:
                print(f"Couldn't find an executable -- {executable}")
                compilers_found[vendor].append(False)
    #
    compilers_to_use = [] 
    for vendor in compilers.keys():
        if all(compilers_found[vendor]):
            print(f"Found suitable compilers from {vendor}")
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
                print(f"Found -- {compiler} -- in:\n{compiler_location}")
            else:
                raise Exception("Couldn't identify {compiler} compiler!\nPlease, install at least C & C++ compilers!")

    return compile_compilers, compilers_to_use

"""
Description: Method checks for existence of MPI on the system from different vendors.
"""
def check_mpi(version_string, env):

    # TODO: figure out problem with intel -- basically 
    # it has different version string, so we need to use anotehr version string,
    # OR implement MPI compilation whithin cmake and thus no need to do anything here.
    compile_mpi = False 
    mpi_to_use = ["mpicc", "mpic++", "mpifort"]
    return compile_mpi, mpi_to_use
"""
Description: Method checks all available projects. I starts with the compilers and then goes
for other prerequisites. If NO COMPILERS WERE FOUND than it will give you an error! There should
be something on the system beforehand.
"""
def check_prerequisites(prerequisites, env, download_dir, num_procs, signatures):
    # Creating regex pattern to match versions
    version_string = '(?:(\d+\.[.\d]*\d+))'
    # Checking if compilers are available on the system
    compile_compilers, compilers_to_use = check_compilers(version_string, env)
    # Checking if mpi wrappers are available on the system
    compile_mpi, mpi_to_use = check_mpi(version_string, env)
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
        print(f"============================================")
        print(f"Start working with project - {proj}...")
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
                print(f"Found -- {executable}, v{exe_version} -- in:\n{exe_location}")

                if version.parse(exe_version) >= version.parse(min_version):
                    print("Current version is suitable for building the project.")
                    build_from_src[proj] = True 
                else:
                    print("Current version is not enough for building the project.")
                    build_from_src[proj] = False 
            else:
                print(f"Couldn't find an executable -- {executable}")
                build_from_src[proj] = False

        # Downloading and compiling stuff from source
        # TODO: tune GCC installation for Mac?
        if build_from_src[proj] == False:
            get_source(details, env, download_dir, signatures)
            compile_source(details, env, download_dir, install_dir, num_procs)
            print(f"Installed project is:")
            subprocess.run([f"{proj}","--version"], shell=False, check=True, env=env)
            print(f"--------------------------------------------")

    return compilers_to_use, mpi_to_use

#----------------------------------------------------------------------
# Creating Path object
root_dir = Path(__file__).resolve().parent
# Getting arguments
args             = get_cmd_args(root_dir)
build_dir        = Path(args.build_dir).resolve()
install_dir      = Path(args.install_dir).resolve()
download_dir     = Path(args.download_dir).resolve()
num_procs        = args.num_procs
c_compiler       = args.c_comp
cxx_compiler     = args.cxx_comp
fortran_compiler = args.fortran_comp
# Creating directories if they do not exist
if not Path.is_dir(build_dir):
    Path.mkdir(build_dir)
if not Path.is_dir(install_dir):
    Path.mkdir(install_dir)
if not Path.is_dir(download_dir):
    Path.mkdir(download_dir)

#----------------------------------------------------------------------
print(f"============================================")
print(f"Start analyzing. please, wait a moment")
#----------------------------------------------------------------------
# Getting prerequisites
prerequisites = str(root_dir.joinpath("prerequisites.json"))
# read file
with open(prerequisites, 'r') as myfile:
    data = myfile.read()
# parse file
prerequisites = json.loads(data)
# Updating the environment variables safely:
# Make a copy of the environment and pass is to the childen.
comm3_env = os.environ.copy()
comm3_env["PATH"] = f"{install_dir}/bin:{os.environ.get('PATH')}"
comm3_env["LD_LIBRARY_PATH"] = f"{install_dir}/lib:{os.environ.get('LD_LIBRARY_PATH')}"

#----------------------------------------------------------------------
# Getting GNU signatures for verification
# Downloading GNU GPG keys
gnu_signatures_source_url = "https://ftp.gnu.org/gnu/gnu-keyring.gpg"
gnu_signatures_file = "gnu-keyring.gpg"
full_path_to_gnu_signatures_file = download_dir / gnu_signatures_file 
perform_download(gnu_signatures_source_url, gnu_signatures_file, full_path_to_gnu_signatures_file)
#----------------------------------------------------------------------
# We get compilers with which we will compile commander
compilers_to_use, mpi_to_use = check_prerequisites(prerequisites, comm3_env, download_dir, num_procs, full_path_to_gnu_signatures_file)
print(f"--------------------------------------------")
print(compilers_to_use)
print(mpi_to_use)
print(f"--------------------------------------------")
#----------------------------------------------------------------------
"""
CMake part is below
"""
project_dir = "/Users/maksymb/Desktop/work/Projects/Commander/CommanderSuperbuild/build"
cmake_configure_commands = ["cmake", 
                  f"-DCMAKE_INSTALL_PREFIX={install_dir}",
                  f"-DCMAKE_C_COMPILER={compilers_to_use[0]}",
                  f"-DCMAKE_CXX_COMPILER={compilers_to_use[1]}",
                  f"-DCMAKE_Fortran_COMPILER={compilers_to_use[2]}",
                  ".."
                  ]
cmake_install_commands = ["cmake", "--build", ".", "--target", "install", "-j", f"{num_procs}"]
result = subprocess.run([f"cmake",".."], 
        shell=False, check=True, env=comm3_env, cwd=project_dir)#,
#        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
result = subprocess.run([f"cmake","--build", ".", "-j", f"{num_procs}"], 
        shell=False, check=True, env=comm3_env, cwd=project_dir)#,
#        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
print(f"--------------------------------------------")
#----------------------------------------------------------------------
"""
Description: Override the base extension class for obstructing setuptools from
building the project for us. All subprojects will be compiled with either python
from source using custom scripts or CMake.
"""
class CMakeExtension(Extension):
    def __init__(self, name, sources=[]):
        super().__init__(name = name, sources = sources)


"""
Description: Builds using CMake instead of the python setuptools
"""
class CMakeBuild(build_ext):
    """
    Perform build_cmake before doing the 'normal' stuff
    """
    def run(self):
        for extension in self.extensions:
            if extension.name == 'example_extension':
                self.build_cmake(extension)

        super().run()

    def build_cmake(self, extension: Extension):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build Commander3"
            )


