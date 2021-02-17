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

#----------------------------------------------------------------------
"""
Description: This module contains miscellanious tools to download. 
check and unpack the project.
"""
class Downloader(object):

    """
    Description: This method calculates MD5 hash sum of an archive and compares it to
    the given by the authors.
    """
    @staticmethod
    def verify_md5(archive_md5_hash, file_to_check, source_file):
        print(f"--------------------------------------------")
        md5_hash = hashlib.md5()
        with open(file_to_check, "rb") as current_file:
            # Read and update hash in chunks of 4K
            for byte_block in iter(lambda: current_file.read(4096), b""):
                md5_hash.update(byte_block)
            print(f"-- Calculated md5 hash is {md5_hash.hexdigest()}")
            print(f"-- Provided md5 hash is {archive_md5_hash}")
        if (archive_md5_hash == md5_hash.hexdigest()):
            print(f"-- Hashes are the same for {source_file}. Proceeding...")
        else:
            raise Exception(f"Hashes differ! Something went wrong while downloading {source_file}")

    """
    Description: This method calculates SHA1 hash sum of an archive and compares it to
    the given by the authors.
    """
    @staticmethod
    def verify_sha1(archive_sha1_hash, file_to_check, source_file): 
        print(f"--------------------------------------------")
        sha1_hash = hashlib.sha1()
        with open(file_to_check, "rb") as current_file:
            # Read and update hash in chunks of 4K
            for byte_block in iter(lambda: current_file.read(4096), b""):
                sha1_hash.update(byte_block)
            print(f"-- Calculated sha1 hash is {sha1_hash.hexdigest()}")
            print(f"-- Provided sha1 hash is {archive_sha1_hash}")
        if (archive_sha1_hash == sha1_hash.hexdigest()):
            print(f"-- Hashes are the same for {source_file}. Proceeding...")
        else:
            raise Exception(f"Hashes differ! Something went wrong while downloading {source_file}")

    """
    Description: This method calculates SHA256 hash sum of an archive and compares it to
    the given by the authors.
    """
    @staticmethod
    def verify_sha256(archive_sha256_hash, file_to_check, source_file): 
        print(f"--------------------------------------------")
        sha256_hash = hashlib.sha256()
        with open(file_to_check, "rb") as current_file:
            # Read and update hash in chunks of 4K
            for byte_block in iter(lambda: current_file.read(4096), b""):
                sha256_hash.update(byte_block)
            print(f"-- Calculated sha256 hash is {sha256_hash.hexdigest()}")
            print(f"-- Provided sha256 hash is {archive_sha256_hash}")
        if (archive_sha256_hash == sha256_hash.hexdigest()):
            print(f"-- Hashes are the same for {source_file}. Proceeding...")
        else:
            raise Exception(f"Hashes differ! Something went wrong while downloading {source_file}")

    """
    Description: This method calculates SHA512 hash sum of an archive and compares it to
    the given by the authors. Raises an error if hashes do not match.
    """
    @staticmethod
    def verify_sha512(archive_sha512_hash, file_to_check, source_file): 
        print(f"--------------------------------------------")
        sha512_hash = hashlib.sha512()
        with open(file_to_check, "rb") as current_file:
            # Read and update hash in chunks of 4K
            for byte_block in iter(lambda: current_file.read(4096), b""):
                sha512_hash.update(byte_block)
            print(f"-- Calculated sha512 hash is {sha512_hash.hexdigest()}")
            print(f"-- Provided sha512 hash is {archive_sha512_hash}")
        if (archive_sha512_hash == sha512_hash.hexdigest()):
            print(f"-- Hashes are the same for {source_file}. Proceeding...")
        else:
            raise Exception(f"Hashes differ! Something went wrong while downloading {source_file}")

    """
    Description: This method verifies archive signature files (.tar.sig). It raises an error
    if signature cannot be verified.
    """
    @staticmethod
    def verify_signature(details, download_dir, file_to_check, signatures, env):
        # Converting Path into string values
        signatures = str(signatures)
        file_to_check = str(file_to_check)
        # Performing download of an archive signature file
        archive_signature_url = details["signature"]
        archive_name = Path(Path(details["signature"]).stem + Path(details["signature"]).suffix) 
        archive_signature_file = download_dir / archive_name 
        Downloader.perform_download(archive_signature_url, archive_name, archive_signature_file)
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
            print(f"-- Signature has been verified. Details:")
            print(f"-- Fingerprint: {verified.fingerprint}")
            print(f"-- Signer: {verified.username}")
            print(f"-- Sign date: {verified.creation_date}")
            print(f"-- Signature status: {verified.status}")
            print(f"-- Trust level: {verified.trust_level}")
            print(f"-- Trust text: {verified.trust_text}")
            print(f"-- Key status: {verified.key_status}")

    """
    Description: This method downlaods source archive if it has not been downloaded yet.
    """
    @staticmethod
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
            print(f"-- Getting -- {source_file} -- from:\n{source_url}")
            if response.status_code == 200: 
                with open(full_path_to_source, 'wb') as current_file:
                    current_file.write(response.raw.read())
            else:
                raise Exception(f"-- Something went wrong while downloading {source_file}")
        else:
            print(f"-- File -- {source_file} -- already exists! Skipping this step...")

    """
    Description: This method veerifies downloaded file via various methods. ALL AVAILABLE METHODS
    WILL BE CHECKED. IF ANY OF PROVIDED METHODS DOESN"T WORK, THE ERROR WILL BE RAISED AND SCRIPT
    STOPS.
    """
    @staticmethod
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
                    Downloader.verify_signature(details, download_dir, file_to_check, signatures, env)
                if (check == 'MD5'):
                    Downloader.verify_md5(details['MD5'], file_to_check, source_file)
                if (check == 'SHA1'):
                    Downloader.verify_sha1(details['SHA1'], file_to_check, source_file)
                if (check == 'SHA256'):
                    Downloader.verify_sha256(details['SHA256'], file_to_check, source_file)
                if (check == 'SHA512'):
                    Downloader.verify_sha512(details['SHA512'], file_to_check, source_file)
        print(f"--------------------------------------------")

    """
    Description: Method to unpack archives.
    """
    @staticmethod
    def unpack_archive(project_file, project_dir, download_dir):
        #project_dir = Path(project_dir)
        #print(project_dir)
        #print(type(project_dir))
        #print(type(project_file))
        # Unpacks archive unless the file is already present on disk
        if not Path.is_dir(project_dir):
        #if not os.path.isdir(project_dir):
            print(f"-- Unpacking file...")
            shutil.unpack_archive(str(project_file), str(download_dir))#, project_dir)
        else:
            print(f"-- Archive's already unpacked! Skipping this step...")

    """
    Description: Method to download source archives, check them with provided options
    and unpacks them to a specified location.
    """
    @staticmethod
    def get_source(details, env, download_dir, signatures):
        # Performing download of an archive
        archive_url               = details["source"]
        project_archive           = Path(Path(details["source"]).stem + Path(details["source"]).suffix) 
        full_path_to_project_file = download_dir / Path(project_archive) 
        Downloader.perform_download(archive_url, project_archive, full_path_to_project_file)
        # Perfoming check for an archive corruption using signature file, md5 or other hash.
        Downloader.check_source(details['checker'], download_dir, project_archive, full_path_to_project_file, signatures, env)
        # Once archive is verified with its digital signature or hashes, we proceed with installation
        # TODO: finish this part -- unpack archives
        project_dir = download_dir / Path(str(project_archive).replace(".tar.gz",""))
        Downloader.unpack_archive(full_path_to_project_file, project_dir, download_dir)
