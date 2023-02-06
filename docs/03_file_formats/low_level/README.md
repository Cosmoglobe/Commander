# Low-level Data Files

## Instrument file

An HDF instrument file is used to collect different types of
instrument parameters in a convenient format. Since this serves a
similar role as the Planck Reduced Instrument MOdel (RIMO), this file
is referred to as the Commander RIMO. A full specification is provided
in the table below. However, not all elements must be present for all
channels, but only those that actually will be used in a given
analysis.

| Quantity      | Description |   HDF path |
| ----------    | ----------- | -------------------- |
| Bandpass specification      | Detector sensitivity as a function of frequency | `[channel label]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# Channel label, e.g., '030'`<br>&nbsp; &nbsp; &nbsp;`centFreq (dp)`&nbsp; &nbsp; &nbsp; `# Band center frequency in GHz`<br>&nbsp; &nbsp; &nbsp;`bandpassx[n] (dp)`&nbsp; &nbsp; &nbsp; `# Frequency in GHz`<br>&nbsp; &nbsp; &nbsp;`   bandpass[n]  (dp)`&nbsp; &nbsp; &nbsp; `# Bandpass amplitude`<br> where $n$ indicates the number of frequency samples |
| 4pi beam    | Spherical harmonics decomposition of full $4\pi$ beam. Real harmonic coefficients are listed using the Libsharp order. Used primarily for orbital dipole convolution | `[channel label]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# Channel label, e.g., '030'`<br>&nbsp; &nbsp; &nbsp;`beamlmax (int)`&nbsp; &nbsp; &nbsp; `# Maximum multipole`<br>&nbsp; &nbsp; &nbsp;`   beammmax  (int)`&nbsp; &nbsp; &nbsp; `# Maximum azimuthal order`<br>&nbsp; &nbsp; &nbsp;`   beam  (group)`&nbsp; &nbsp; &nbsp; `# Group name`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   T(n)  (dp)`&nbsp; &nbsp; &nbsp; `# Intensity beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   E(n)  (dp)`&nbsp; &nbsp; &nbsp; `# E-mode beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   B(n)  (dp)`&nbsp; &nbsp; &nbsp; `# B-mode beam`<br>where $n$ is the number of spherical harmonics determined by `beamlmax` and `beammmax`|
| Far sidelobe beam    | Spherical harmonics decomposition of far sidelobes. Real harmonic coefficients are listed using the Libsharp order. Note that this usually has a lower resolution than the $4\pi$ beam | `[channel label]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# Channel label, e.g., '030'`<br>&nbsp; &nbsp; &nbsp;`sllmax (int)`&nbsp; &nbsp; &nbsp; `# Maximum multipole`<br>&nbsp; &nbsp; &nbsp;`   slmmax  (int)`&nbsp; &nbsp; &nbsp; `# Maximum azimuthal order`<br>&nbsp; &nbsp; &nbsp;`   sl  (group)`&nbsp; &nbsp; &nbsp; `# Group name`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   T(n)  (dp)`&nbsp; &nbsp; &nbsp; `# Intensity beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   E(n)  (dp)`&nbsp; &nbsp; &nbsp; `# E-mode beam`<br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`   B(n)  (dp)`&nbsp; &nbsp; &nbsp; `# B-mode beam`<br>where $n$ is the number of spherical harmonics determined by `sllmax` and `slmmax`|
| Beam summary statistics    | Various beam summary statistics | `[channel label]` &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; `# Channel label, e.g., '030'`<br>&nbsp; &nbsp; &nbsp;`fwhm (dp)`&nbsp; &nbsp; &nbsp; `# FWHM in arcmin`<br>&nbsp; &nbsp; &nbsp;`   elip  (dp)`&nbsp; &nbsp; &nbsp; `# Beam ellipticity`<br>&nbsp; &nbsp; &nbsp;`   mbeam_eff  (dp)`&nbsp; &nbsp; &nbsp; `# Main beam efficiency` <br>&nbsp; &nbsp; &nbsp;`   psi_ell  (dp)`&nbsp; &nbsp; &nbsp; `# Ellipticity rotation angle in degrees`|

## TOD datafile

The TOD data files contain the actual time-ordered information for the detectors. These files 
will contain the majority of the data volume of the entire project. The files are in HDF5 format
and use multiple levels of indexing. The top level references the unique chunk index, and the second
level refers to the detector names. Chunks can have arbitrary length, but should be chosen such 
that the following conditions are true: 1) the noise is stationary over the entire chunk, and 2)
the data volume of a single chunk is in the 1-100 MB range. This seems to be the size that gives
the best read performance on our systems.

| Quantity | Description | Location |
| -------- | ------------| ---------|
| **Common Parameters** | | |
| Detector Labels | A comma separated string of all the detector names that appear in this file | /common/det |
| Sampling Frequency | The sampling frequency of these detectors in Hz | /common/fsamp |
| Main beam angle | The main beam offset polarization angles | /common/mbang |
| Number of psi bins used in compression | This should be populated automatically by the file generation code. For performance reasons we suggest something in the form 2^n | /common/npsi |
| nside of the pointing | The nside at which the pointing is compressed. | /common/nside |
| pid list in this file | The unique pointing IDs (chunks) that are present in this file | /common/pid |
| polarization angles | The polarization angles of the detectors. | /common/polang |
| version information | A unique version identifier so you can tell your files apart | /common/version |
| **Per-chunk parameters**| | |
| Huffman symbol list | This should be generated automatically | /[chunk_num]/common/huffsymb |
| Huffman tree array | This should be generated automatically | /[chunk_num]/common/hufftree |
| Load balancing parameters | An array of shape (2,1) that contains some estimate of this chunk's proximity to other chunks. Chunks with similar parameters will be given to the same core. | [chunk_num]/common/load |
| Number of samples | The length of the tod, flag and pointing arrays in this chunk. | /[chunk_num]/common/ntod |
| Satellite position | The x,y,z positions of the satellite as a function of time in solr system coordinates. Ground based experiments should use the position of the earth. | /[chunk_num]/common/satpos |
| Time | The time at the start of this chunk. Space is given for 3 different units if desired. | /[chunk_num]/common/time |
| Satellite velocity | The x,y,z velocity of the satellite relative to the sun. Used to calculate the orbital dipole. | /[chunk_num]/common/vsun | 
| **Per-detector parameters** | | |
| Flags | The huffman compressed flag field | /[chunk_num]/[detector_name]/flag |
| Pointing | The healpy and huffman compressed pointing information | /[chunk_num]/[detector_name]/pix |
| Psi | The huffman compressed polarization angle | /[chunk_num]/[detector_name]/psi |
| Other scalar information | A length 4 vector that contains default values for (gain, sigma0, fknee, alpha). | /[chunk_num]/[detector_name]/scalars |
| Data | The actual measurements of the sky | /[chunk_num]/[detector_name]/tod |

## TOD filelist
The filelist is a text file which contains a list of all the data chunks that are to be analyzed. 

| Column | Quantity | Descrition |
| ------ | -------- | ---------- |
| Column 1, row 1 | File length | The number of entries in the file. For the full dataset, it should be the number of lines in the file-1. If you only want some of the data, you can lower this number to read only the first n lines. If it is longer than the number of entries, the code will crash. |
| 1 | Chunk ID | A unique list of all the chunk IDs that appear in the hdf file. These do not have to be in order. |
| 2 | Path | A full path to the file that contains that chunk |
| 3 | Time estimate | A time estimate about how long this chunk would take to process. This field is output from Commander during a run, re-writing this file. To build the file in the first place, just use 1 here |
| 4 + 5 | Load balancing parameters | The same load balancing numbers that appear in the HDF files. Chunks with similar numbers here get grouped togeather. Can just used 0, 0 to ignore this optimization. |
