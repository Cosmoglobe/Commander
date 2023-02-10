# Downloader Utility

We provide a command line utility that downloads all required Commander
files.

## Installation

Download the utility from our releases page at:

https://gitlab.com/BeyondPlanck/r13y-helper/-/releases

We provide two statically compiled binaries, one for Linux/Mac (`bp`)
and one for Windows (`bp.exe`), which we will be refering to as `bp`
from now on. Download what is appropriate for your OS.

## Configuration

The first time `bp` executes it will create in the current folder a
`.bp` subfolder that will contain a `config.toml` configuration file,
with contents:

```
[downloads]
# How many concurrent downloads you would like to run at the same time.
# The default of 3 seems to provide a sane starting point.
max_concurrent_downloads=3

# The location where all BeyondPlanck input files will be downloaded and used from.
# Since this potentially will require a lot of filespace you should put it in a
# partition that has adequate free size.  The `download` command will inform you of
# file space needed.
download_path="/data/disk2"

[app]

# You can increase the verbosity of the output, if you would like to
# debugg certain issues.  Note, that this will create a very busy output
# stream
verbose=false
```

The configuration file is heavily commented and provides some sane
defaults.

> **Note**:  Due to the large size of the files that the download
  tool downloads make sure you set a `download_path` that resides in a
  disk partition with enough disk space.

## Usage

Copy your `bp` executable (or `bp.exe`),  to your current folder or put
it in a folder that is in your path.

The `bp` tool allows you to download multiple filesets:

`bp download`

This command will show you all the available filesets to download and
their download size, without actually downloading anything.

You can select to download `all` files or a specifc list of filesets
only. For a successfull Commander execution the `auxcmd3` and `l2data`
are **required**.

Run:

`bp download auxcmd3,l2data`

The tool will start downloading the specified files.  Note, that
depending on your network connection and your available bandwidth, the
download process can take quite some time.

If your connection drops, you can restart the operation and the tool
should be able to continue downloading from where it left off,
**without** having to re-download everything again.
