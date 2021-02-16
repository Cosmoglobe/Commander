import os
# 
import sys
# to get data from internet
import requests
# to unpack archavies
import shutil
# to run shell commands
import subprocess
# to read json files
import json
# find hash values of a file
import hashlib
# verification of tar.gz.sig files
import gnupg
#from pretty_bad_protocol import gnupg
#from pprint import pprint as print
# for working with lists
from collections import defaultdict
# getting compiler (and other) versions
import re
#
#import cmaketools
# OOP way to handle paths in python. Included in Standard Library from python 3.6
from pathlib import Path


#----------------------------------------------------------------------
# Main variables
# TODO: make them as command line variables
num_proc = 6
# Projects to install -- add them to json file or something like that.
# Creating download dir if needed
root_dir = Path(__file__).resolve().parent.absolute()
build_dir = root_dir.joinpath("build")
download_dir = build_dir.joinpath("downloads")
install_dir = build_dir.joinpath("install")
# Creating directories if they do not exist
if not Path.is_dir(build_dir):
    Path.mkdir(build_dir)
if not Path.is_dir(install_dir):
    Path.mkdir(install_dir)
if not Path.is_dir(download_dir):
    Path.mkdir(download_dir)
#if not os.path.isdir(build_dir):
#    os.mkdir(build_dir)
#if not os.path.isdir(install_dir):
#    os.mkdir(install_dir)
#if not os.path.isdir(download_dir):
#    os.mkdir(download_dir)

#----------------------------------------------------------------------
# TODO: Add these to JSON file?
# TODO: The CC, CXX and FC should be properly set before strating this script 
# TODO: Figure out issues with GNU GCC compilers.
projects = {
            "make": {
                "source": ["https://ftp.gnu.org/gnu/make/", "make-4.3", ".tar.gz"],
                "signature": ["https://ftp.gnu.org/gnu/make/","make-4.3.tar.gz.sig"],
                "MD5":    "",
                "SHA1":   "",
                "SHA256": "",
                "commands": [["./configure", f"--prefix={install_dir}"], # "&&", 
                             #["make", "-j", f"{num_proc}"],# "&&",
                             ["sh", "build.sh", "-j", f"{num_proc}"], # "&&",
                             ["./make", "install"]]
                },
            # Note: there is some problems with compiling automake with autoconf 2.70
            # So I switched to 2.69 instead.
            "autoconf": {
                #"source": ["https://ftp.gnu.org/gnu/autoconf/","autoconf-2.70", ".tar.gz"],
                "source": ["https://ftp.gnu.org/gnu/autoconf/","autoconf-2.69", ".tar.gz"],
                #"signature": ["https://ftp.gnu.org/gnu/autoconf/","autoconf-2.70.tar.gz.sig"],
                "signature": ["https://ftp.gnu.org/gnu/autoconf/","autoconf-2.69.tar.gz.sig"],
                "MD5":    "",
                "SHA1":   "",
                "SHA256": "",
                "commands": [#[f"export PATH={install_dir}/bin:$PATH"],
                             #[f"export LD_LIBRARY_PATH={install_dir}/lib:$LD_LIBRARY_PATH"],
                             ["./configure", f"--prefix={install_dir}"],# "&&", 
                             [f"{install_dir}/bin/make", "-j", f"{num_proc}"],# "&&",
                             [f"{install_dir}/bin/make", "install"]]
                },
            "automake": {
                "source": ["https://ftp.gnu.org/gnu/automake/","automake-1.16.3", ".tar.gz"],
                "signature": ["https://ftp.gnu.org/gnu/automake/","automake-1.16.3.tar.gz.sig"],
                "MD5":    "",
                "SHA1":   "",
                "SHA256": "",
                "commands": [#[f"export PATH={install_dir}/bin:$PATH"],
                             #[f"export LD_LIBRARY_PATH={install_dir}/lib:$LD_LIBRARY_PATH"],
                             ["./configure", f"--prefix={install_dir}"],# "&&", 
                             [f"{install_dir}/bin/make", "-j", f"{num_proc}"],# "&&",
                             [f"{install_dir}/bin/make", "install"]]
                },
            "libtool": {
                "source": ["https://ftp.gnu.org/gnu/libtool/","libtool-2.4.6", ".tar.gz"],
                "signature": ["https://ftp.gnu.org/gnu/libtool/","libtool-2.4.6.tar.gz.sig"],
                "MD5":    "",
                "SHA1":   "",
                "SHA256": "",
                "commands": [#[f"export PATH={install_dir}/bin:$PATH"],
                             #[f"export LD_LIBRARY_PATH={install_dir}/lib:$LD_LIBRARY_PATH"],
                             ["./configure", f"--prefix={install_dir}"],# "&&", 
                             [f"{install_dir}/bin/make", "-j", f"{num_proc}"],# "&&",
                             [f"{install_dir}/bin/make", "install"]]
                },
            "gcc": {
                #"source": ["https://mirror.koddos.net/gcc/releases/gcc-9.3.0/","gcc-9.3.0",".tar.gz"],
                "source": ["http://robotlab.itk.ppke.hu/gcc/releases/gcc-9.3.0/","gcc-9.3.0",".tar.gz"],
                #"signature": ["https://mirror.koddos.net/gcc/releases/gcc-9.3.0/","gcc-9.3.0.tar.gz.sig"],
                "signature": ["http://robotlab.itk.ppke.hu/gcc/releases/gcc-9.3.0/","gcc-9.3.0.tar.gz.sig"],
                "MD5":    "",
                "SHA1":   "",
                "SHA256": "",
                "SHA512": "",
                "commands": [#[f"export PATH={install_dir}/bin:$PATH"],
                             #[f"export LD_LIBRARY_PATH={install_dir}/lib:$LD_LIBRARY_PATH"],
                             ["./contrib/download_prerequisites"],
                             ["../configure", f"--prefix={install_dir}", "--enable-languages=default",
                              "--disable-multilib", "--enable-threads"],#, "--enable-checking=release"],
                             [f"{install_dir}/bin/make", "-j", f"{num_proc}"],# "&&",
                             [f"{install_dir}/bin/make", "install"]]
                },
            "openmpi": {
                "source": ["https://download.open-mpi.org/release/open-mpi/v4.1/", "openmpi-4.1.0",".tar.gz"],
                "signature": "",
                "MD5":    "45d272a0541857a40d1808e86833bc15",
                "SHA1":   "760d33ab160370e7cf6262590d3b66d6e34291a8",
                "SHA256": "228467c3dd15339d9b26cf26a291af3ee7c770699c5e8a1b3ad786f9ae78140a",
                "commands": [#[f"export PATH={install_dir}/bin:$PATH"],
                             #[f"export LD_LIBRARY_PATH={install_dir}/lib:$LD_LIBRARY_PATH"],
                             [f"FC={install_dir}/bin/gfortran",
                              f"CC={install_dir}/bin/gcc",
                              f"CXX={install_dir}/bin/g++",
                               "./configure", f"--prefix={install_dir}"],# "&&", 
                             [f"{install_dir}/bin/make", "-j", f"{num_proc}"],# "&&",
                             [f"{install_dir}/bin/make", "install"]]
                }
        }
#----------------------------------------------------------------------
# Placeholder for download function -- it is used extensively throughout the code
"""
Description: This method downlaods source archive if it has not been downloaded yet.
"""
def perform_download(source_url, source_file, full_path_to_source):
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
    # Performing download of an archive signature file
    archive_signature_url = details["signature"][0] + details["signature"][1]
    #archive_signature_file = download_dir + "/" + details["signature"][1]
    archive_signature_file = download_dir / details["signature"][1]
    perform_download(archive_signature_url, details["signature"][1], archive_signature_file)
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

"""
Description: This method veerifies downloaded file via various methods. ALL AVAILABLE METHODS
WILL BE CHECKED. IF ANY OF PROVIDED METHODS DOESN"T WORK, THE ERROR WILL BE RAISED AND SCRIPT
STOPS.
"""
def check_source(details, download_dir, file_to_check, signatures, env):
    possible_checks = ['signature', 'MD5', 'SHA1', 'SHA256', 'SHA512']
    # This variable is only for pretty message. the important one is file_to_check.
    source_file = details['source'][1] + details['source'][2]
    for check in possible_checks:
        value = details.get(check)
        # If  key exists, the value is not "None", "False", "''", or "0"
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

def unpack_archive(project_file, project_dir, download_dir):
    if not os.path.isdir(project_dir):
        print(f"Unpacking file...")
        shutil.unpack_archive(project_file, download_dir)#, project_dir)
    else:
        print(f"Archive's already unpacked! Skipping this step...")

# Downloading source archive
# 1. Download archive source;
# 2. Download signature file;
# 3. Check signature/check MD5/CheckSHA etc.
def get_sources(project, details, download_dir, signatures, env):
    print(f"============================================")
    print(f"Start working with project - {project}...")
    #source = details["source"]

    # Performing download of an archive
    archive_url = details["source"][0] + details["source"][1] + details["source"][2]
    project_archive = details["source"][1] + details["source"][2]
    #full_path_to_project_file = download_dir +"/"+ project_archive 
    full_path_to_project_file = download_dir / project_archive 
    perform_download(archive_url, project_archive, full_path_to_project_file)

    # Perfoming check for an archive corruption using signature file, md5 or other hash.
    #verify_signature(details, download_dir, full_path_to_project_file, signatures)
    check_source(details, download_dir, full_path_to_project_file, signatures, env)
    # Once archive is verified with its digital signature or hashes, we proceed with installation
    #project_dir = download_dir +"/"+ details["source"][1] 
    project_dir = download_dir / details["source"][1] 
    unpack_archive(full_path_to_project_file, project_dir, download_dir)

def compile_sources(project, details, download_dir, env):
    commands = details["commands"]
    source = details["source"]
    # GNU GCC has special commands, so we just account for 
    # it here, otherwise hucky_index is useless.
    hucky_index = 0
    for command in commands:
        separator = " "
        print(f"============================================")
        print(f"Working with project - {project}...")
        print(f"Executing command: {separator.join(command)}")
        #project_dir = os.getcwd() + "/downloads/" + projects[project][1] 
        if project == "gcc":
            if hucky_index == 0:
                #project_dir = download_dir +"/"+ source[1]
                project_dir = download_dir / source[1]
            else:
                #build_dir_gcc = download_dir +"/"+ source[1] + "/build_gcc"
                build_dir_gcc = download_dir / source[1] / build_gcc
                if not os.path.isdir(build_dir_gcc):
                    os.mkdir(build_dir_gcc)
                project_dir = build_dir_gcc 
            hucky_index = hucky_index + 1
        else:
            #project_dir = download_dir +"/"+ source[1] 
            project_dir = download_dir / source[1] 
        # DO NOT USE shell=True <= SECURITY RISK!!! 
        # to capture output on python 3.6 <= we use stdout=subprocess.PIPE 
        # and to capture it as text use universal_newlines=True
        result = subprocess.run(command, shell=False, cwd=project_dir, check=True, env=env) #text=True, check=True, cwd=project, stdout=sproc.PIPE)
        print(f"--------------------------------------------")
        print(f"Installed project is:")
        subprocess.run([f"{project}","--version"], shell=False, cwd=project_dir, check=True, env=env)
        #print(result.stdout)

"""
Description:
"""
def check_gnu_compilers(env):
    pass
"""
Description: Method to check for available compilers and their version(s)
"""
def check_compilers(env):
    # TODO: find version using regular expression
    # We need to have minimum Intel Fortran 18.0, becaus eit has full Fortran 2008 support:
    # https://en.wikipedia.org/wiki/Intel_Fortran_Compiler
    # GCC compilers
    #compilers = ['gcc', 'g++', 'gfortran']
    #compiler_toolchains = {'Intel': ['icc', 'icpc', 'ifort'], 'GNU': ['gcc','g++', 'gfortran']}
    compiler_toolchains = {'GNU': ['gcc','g++', 'gfortran']}
    #compiler_versions = {'GNU': '([0-9]|[1-9][0-9])\.[0-9]\.[0-9]'}
    compiler_versions = {'GNU': '(?:(\d+\.[.\d]*\d+))'}
    # Any version of 9.x and up should be fine.
    compiler_required_versions = {'GNU': [9, 3]}
    #versions = "([0-9]|[1-9][0-9])\.[0-9]\.[0-9]"   
    # variable to store all details about the found compilers. If compilers aren't found it will
    # have value None instead.
    #compiler_toolchain_details = {toolchain: [] for toolchain in compiler_toolchains.keys()}
    compiler_toolchain_details = defaultdict(list) 
    #print(compiler_toolchain_details)
    for toolchain, compilers in compiler_toolchains.items():
        version = re.compile(compiler_versions[toolchain])
        for i, compiler in enumerate(compilers):
            compiler_location = shutil.which(compiler)
            print(f"--------------------------------------------")
            #print(compilers, i)
            if compiler_location is not None:
                print(f"Found -- {compiler} -- in:\n{compiler_location}")
                #result = subprocess.run([f"{compiler}","--version"], 
                compiler_toolchain_details[toolchain].append( subprocess.run([f"{compiler}","--version"], 
                        shell=False, check=True, env=env,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True))
                #print(f"Compiler details are:\n{result.stdout}")
                print(f"Compiler details are:\n{compiler_toolchain_details[toolchain][i].stdout}")
                #print(version.search(compiler_toolchain_details[toolchain][i].stdout).group())
                compiler_version = version.search(compiler_toolchain_details[toolchain][i].stdout).group()
                compiler_version = compiler_version.split('.')
                major = int(compiler_version[0])
                minor = int(compiler_version[1])
                # Comparing minimal and current versions
                if(major > compiler_required_versions[toolchain][0]):
                    print(f"Required version is {compiler_required_versions[toolchain][0]}.")
                    print('Installed version is newer than required.')
                elif(major == compiler_required_versions[toolchain][0]):
                    print('Version is equal to minimal required.')
                else:
                    print('You need to install newer compilers')
            else:
                print(f"Couldn't find {compiler}")
    print(compiler_version)
    print(f"--------------------------------------------")

#----------------------------------------------------------------------
print(f"============================================")
print(f"Start analyzing. please, wait a moment")
# Updating the environment variables safely:
# Make a copy of the environment and pass is to the childen.
#print(f"Old PATH: {os.environ.get('PATH')}")
# TODO: add compiler checks and then set environment variables:
# FC, CC, CXX to respective compilers. Commander3 compiles only 
# with Intel or GNU GCC, but there are tons of others out there 
# as well
comm3_env = os.environ.copy()
comm3_env["PATH"] = f"{install_dir}/bin:{os.environ.get('PATH')}"
comm3_env["LD_LIBRARY_PATH"] = f"{install_dir}/lib:{os.environ.get('LD_LIBRARY_PATH')}"
check_compilers(comm3_env)
#exit()
#----------------------------------------------------------------------
# Getting GNU signatures for verification
# Downloading GNU GPG keys
gnu_signatures_source_url = "https://ftp.gnu.org/gnu/gnu-keyring.gpg"
gnu_signatures_file = "gnu-keyring.gpg"
#full_path_to_gnu_signatures_file = download_dir + "/" + gnu_signatures_file 
full_path_to_gnu_signatures_file = download_dir / gnu_signatures_file 
perform_download(gnu_signatures_source_url, gnu_signatures_file, full_path_to_gnu_signatures_file)
# Performing download and compilation for each project
for project, details in projects.items():
    #print(f"--------------------------------------------")
    #subprocess.run(["autoconf","--version"], shell=False, cwd=".", check=True, env=comm3_env)
    #print(f"--------------------------------------------")
    # TODO: Added for debugging, remove after
    #if project != "openmpi": continue #exit()
    #if project == "gcc": exit()
    #if project != "gcc": continue #exit()
    # Getting source archives
    get_sources(project, details, download_dir, full_path_to_gnu_signatures_file, comm3_env)
    # Compiling and installing projects
    compile_sources(project, details, download_dir, comm3_env)

print(f"============================================")
print(f"Finished installing Prerequisites!")
print(f"============================================")

