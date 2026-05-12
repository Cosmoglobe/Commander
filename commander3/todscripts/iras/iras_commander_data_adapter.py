from cosmoglobe.tod_tools import CommanderDataAdapter
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np 
import healpy
import functools
import glob
import os, signal

from iras_native_tod_reader import IRASTODReader, BANDS, BAND_DETS
import iras_native_tod_reader


from astropy.time import Time, TimeDelta
from astropy.coordinates import HeliocentricMeanEcliptic, get_body
from astropy.coordinates import SkyCoord, FK4

"""
This script is essentially an adaption of the script
    - Commander/commander3/todscripts/AKARI/akari_commander_data_adapter.py 
written by Eirik Gjerløw for AKARI. 
"""


"""
FROM "Cosmoglobe/cosmoglobe/tod_tools/commander_data_adapter.py:

    This class represents the interface between an experiment and the
    CommanderHDFWriter, which writes the hdf files for a given experiment needed
    for Commander to run.

    We define each experiment as consisting of several 'bands' which then again
    are subdivided into 'detectors', both of which are represented by strings.

    Each HDF file output represents a 'segment': Typically 500MB-1GB of data that
    consists of several 'chunks'. These are derived from the LFI concepts of
    operational days and scans, but are intended more generally. A chunk could be
    roughly an hour of data, typically.

    Chunks and segments are kept track of using integer indices, and
    it is up to the adapter to make sure these indices are consistent with
    themselves. Both segments and chunks are 1-indexed, NB! So if get_num_segments
    is implemented as returning 479, for example, the segments that will be looped
    over are from 1 to 479, inclusive.
"""


"""
Currently missing/placeholder only:
- satvel
- earthpos, 
    - start and end
- satpos
    - Using earthpos atm
- flags
    - everything
- mbangs, polangs, npsi
    - not sure if relevant 
    
"""

BANDS               = iras_native_tod_reader.BANDS
BAND_DETS           = iras_native_tod_reader.BAND_DETS

class IRASCommanderDataAdapter(CommanderDataAdapter):
    """
    An IRAS-specific implementation of the CommanderDataAdapter 

    All non-underlined functions are defined as in that class.  
    """
    def __init__(
            self,
            data_dir: str, 
            nside: int,
            version: int,
            bands: list[str] = BANDS,
            band_dets: dict[str, list[int]] = BAND_DETS,
            num_files_per_segment: int = 100, 
    ):
        self.data_dir       = data_dir
        self.nside          = nside
        self.version        = version
        self.bands          = bands
        self.band_dets      = band_dets
        self.num_detectors  = {key: len(dets) for key, dets in band_dets.items()} 
        self.num_files_per_segment = num_files_per_segment

        assert type(self.version) is int

        self.flux2MJy_sr = self._get_flux2MJy_sr_conversion_factors()

        self.fsamp = {
            '012': 16.0,  # [Hz],  12 micron
            '025': 16.0,  # [Hz],  25 micron
            '060': 8.0,   # [Hz],  60 micron
            '100': 4.0,   # [Hz], 100 micron
        }

        self.t0 = Time("1981-01-01T00:00:00", scale="utc") # IRAS's t=0
        self.GC = SkyCoord(l=0*u.degree, b=0*u.degree, frame='galactic') # For coord trans.


        # All (alive) detectors have the same number of files
        # Previously used det 30 to compute the number.    
        # self.num_chunks     = len(glob.glob(f"{self.data_dir}/*/*/det30_continuous.npz"))
        # Use manually fixed number for simplicity/safety
        self.num_chunks     = int(5787) 

        self.num_segments   = int(np.ceil(self.num_chunks / self.num_files_per_segment))

        self.todreaders = {}
        self.filelist = {}
        self.chunksizes = {}

        for band in self.bands:
            # print(f"Preparing files for band {band}", end="...\n", flush=True)
            self.todreaders[band] = IRASTODReader(
                iras_npz_dir    = self.data_dir, 
                band            = band, 
                band_dets       = self.band_dets,
                load_filelist   = True,
                )
            self.filelist[band] = self.todreaders[band].filelist
            self.chunksizes[band] = self.todreaders[band].chunksizes
        # print("Done!")
        # print()
        # for i in range(1, 1 + self.get_num_segments(band)):
        #     c = self.get_chunk_indices(band, i)
        #     print(f"{len(c)=}: {c[:5], c[-5:]}")
        # exit()


    def _process_iras_flags(self):
        pass

    def _process_iras_tod_reader_data(self):
        if True:
            raise NotImplementedError("Placeholder method atm.") 

        """Internal function to process the data coming from the IRAS TOD reader and getting it
        into the format needed for the HDFwriter"""

        self.todreaders = None
        self._process_iras_flags = self._process_iras_flags


        todreader_data = {}
        for band in self.bands:
            todreader_data[band] = {}
            files = self.chunk_file_map[band][self.chunk_idx]

            # This does not yet include any code for when the start and end happens
            # within a file rather than at the beginning of the file.
            start_idx = self.todreaders[band].get_file_index_range(files[0])[0]
            end_idx = self.todreaders[band].get_file_index_range(files[-1])[1]
            curr_data = self.todreaders[band].get_data(start_idx, end_idx)
            mode = curr_data['packet_id']
            for det in self.band_dets[band]:
                todreader_data[band][det] = {}
                todreader_data[band][det]['tod'] = curr_data[f'{det}/tod']
                todreader_data[band][det]['flag'] = self._process_iras_flags(
                    curr_data[f'{det}/pixel_flag'],
                    curr_data['frame_flag'],
                    curr_data['status_flag'], mode)
            ra = curr_data['ra']
            dec = curr_data['dec']
            c = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
            todreader_data[band]['pix'] = healpy.ang2pix(self.nside,
                                                         c.galactic.l.value,
                                                         c.galactic.b.value,
                                                         lonlat=True)
            if 'start_time' not in todreader_data:
                todreader_data['start_time'] = curr_data[f'aftime'][0]
                todreader_data['end_time'] = curr_data[f'aftime'][-1]
                todreader_data['start_satpos'] = curr_data['start_satpos']
                todreader_data['end_satpos'] = curr_data['end_satpos']
                todreader_data['start_earthpos'] = curr_data['start_earthpos']
                todreader_data['end_earthpos'] = curr_data['end_earthpos']

        return todreader_data

    def _get_flux2MJy_sr_conversion_factors(self):
        # Create dictionary containing each detector's multiplicative factor 
        # to convert tod from W/m^2 to MJy/sr.
        # Returns dict with
        #   key: detector name [str]
        #   val: conversion factor [float] 

        # Load dict. with det. areas: det_area_dict[det.numb(int)] = area (steradian). 
        # Product of last two columns of Table II.C.3 in IRAS explanatory supplement
        det_file = "/mn/stornext/d5/data/vetleav/IRAS/detector_area.npy"
        det_area_dict = np.load(det_file, allow_pickle=True).item()

        conv_factor_dict = {}
        for band, det_list in self.band_dets.items():
            band_wl = float(band) * u.micron
            band_freq = band_wl.to(u.Hz, equivalencies=u.spectral())
            for det_name in det_list:
                det_num = int(det_name.split("-")[-1])
                det_area = det_area_dict[det_num]
                conv_factor_dict[det_name] = (1 * u.W / u.m**2 / band_freq / det_area).to(u.MJy / u.sr).value
                
        return conv_factor_dict


    def get_num_segments(self, band: str)  ->  int: 
        # For a given band, returns the
        # number of segments (i.e. number of HDF files) for that band.
        # All IRAS detectors have 5787 npz files, i.e. chunks.
        #   - 5787 (for deplated Level 1 data)
        return self.num_segments

    def get_version(self) -> int: 
        #Returns the current version of the HDF generation.

        # Version 1:
        #   - Deplated level1 data from IPAC, in units of W/m^2.
        #   - Minimum working values with proper positions, vels, 
        #       or flags not properly implemented yet. 
        return self.version

    def get_experiment_name(self) -> str: 
        return "IRAS"

    def get_npsi(self) -> int: 
        # Returns the number of psi bins to use for the psi compression.

        # Taken from "/home/vetle/Documents/PhD_astro/Cosmoglobe/cosmoglobe/tod_tools/iras/comm_iras_data_adapter.py"
        # Line 31: self.npsi = 8  # For IRAS this shouldn't matter
        # WILL/SHOULD CHECK LATER 
        return 64
    
    def get_nside(self, band: str) -> int: 
        # Given a band name, returns an integer specifying
        # the nside of the pixelization used for the pointing.
        return self.nside

    def get_polangs(self, band: str) -> list[float]: 
        # Given a band, returns a list of length ndet 
        # which gives the polarization angles for each detector of
        # that band.
        return [0] * self.num_detectors[band]
    
    def get_mbangs(self, band: str) -> list[float]: 
        # Given a band name, returns a list of length ndet which gives 
        # the main beam angle of that each detector of that band.
        return [0] * self.num_detectors[band]

    def get_fsamp(self, band:str) -> float: 
        # Given a band name, returns the sampling
        # rate of that band.
        return self.fsamp[band]
        
    def get_detector_names(self, band:str) -> list[str]: 
        # Given a band name, returns a
        # list of strings specifying the detector names of that band - these will
        # be used as fields in the output HDF file.
        return self.band_dets[band]

    def get_chunk_indices(self, band:str, segment:int) -> list[int]: 
        # Given a band
        # name and a segment number belonging to that band, gives the indices of
        # the chunks (previously:scans) that belong to the segment in question.
        # Called by "commander_hdf_writer.py"
        start_idx   = (segment-1) * self.num_files_per_segment
        end_idx     = segment * self.num_files_per_segment
        if end_idx > self.num_chunks:
            end_idx = self.num_chunks
        return [i for i in range(start_idx+1, end_idx+1)]

    def get_bands(self) -> list[str]: 
        # Returns a list of the band names.
        return self.bands

    def set_chunk_index(self, chunk_idx:int):
        # Since it can be bad,
        # performance-wise, to have to specify the chunk for every low-level data
        # fetch operation, this function allows the reader to signal to the
        # adapter which chunk it is going to process next. After this, low-level
        # data fetch calls will *not* be specifying the chunk number, but will
        # implicitly assume that the adapter knows which chunk number is the
        # current one.

        # Called by "commander_hdf_writer.py"
        self.chunk_idx = chunk_idx
        self.all_chunk_data = None #Reset chunk

    def get_chunk(func):
        if True:
            raise NotImplementedError("Kept for reference atm.") 
        
        """Decorator function to make sure that a chunk is always preloaded"""
        print("NOOOOOOOO")
        @functools.wraps(func)
        def inner(self, *args):
            if self.all_chunk_data is None:
                self.all_chunk_data = self._process_iras_tod_reader_data()
            return func(self, *args)

        return inner

    def get_chunk_data(self, band:str, detector:str): 
        # -> tod(np.array[float/int]), pix(np.array[int]), psi(np.array(float)), flag(np.array[int]):
        # Given a band and detector (with the chunk number assumed set by
        # set_chunk_index), returns the tod, pixel, psi, and flag array
        # corresponding to the detector and chunk in question. The flag array is
        # an array of integers representing bits.

        fname = self.filelist[band][detector][self.chunk_idx - 1]

        ### KEYS: t, ra, dec, flux, n_merged, dedup_seconds, sop, obs, det
        data = np.load(fname)

        # Some bands (60 and 100) have varying ntods betw. detectors for certain cunks
        # N_tod is precomputed, and ensures that all detectors have the same length for all chunks
        # For bands 012 and 025 N_tod=-1 for all chunks.
        N_tod = self.chunksizes[band][self.chunk_idx]

        t   = data["t"][:N_tod]
        self.chunk_start_time   = t[0]  * u.s + self.t0
        self.chunk_end_time     = t[-1] * u.s + self.t0

        tod = data["flux"][:N_tod] * self.flux2MJy_sr[detector] # Convert data from W/m^2 to MJy/sr
        ra  = data["ra"][:N_tod]
        dec = data["dec"][:N_tod]

        # TEMP
        flag_tot = np.zeros_like(tod, dtype=int)
        flag_tot[~np.isfinite(ra)]  += 2**0
        flag_tot[~np.isfinite(tod)] += 2**1

        tod[flag_tot != 0]  = 0
        ra[flag_tot  != 0]   = 0
        dec[flag_tot != 0]  = 0

        psi = np.zeros_like(tod, dtype=int)
        good_data = flag_tot == 0
        sc = SkyCoord(ra=ra[good_data], dec=dec[good_data], unit='deg',
        # sc = SkyCoord(ra=ra, dec=dec, unit='deg',
                equinox='B1950.0', obstime='J1983.5', frame=FK4)

        coords = sc.transform_to(self.GC)
        ra[good_data] = coords.l.value
        dec[good_data] = coords.b.value
        pix = healpy.ang2pix(self.nside, ra, dec, lonlat=True)
        return tod, pix, psi, flag_tot


    def get_chunk_start_time(self): # -> astropy.time.Time: 
        # Returns the time of the
        # beginning of the chunk.
        return self.chunk_start_time

    def get_chunk_end_time(self): # -> astropy.time.Time: 
        # Returns the time of the end of
        # the chunk.
        return self.chunk_end_time

    def get_chunk_start_satpos(self) -> list[float]:
        # Returns the position of the
        # satellite at the beginning of the chunk in galactic coordinates
        # (cartesian, in meters) as a length 3 list.

        # return self.get_chunk_start_earthpos()
        return [0,0,0] # TEMP

    def get_chunk_end_satpos(self) -> list[float]:
        # Returns the position of the
        # satellite at the end of the chunk in galactic coordinates (cartesian,
        # in meters) as a length 3 list. 

        # return self.get_chunk_end_earth()
        return [0,0,0] # TEMP

    def get_chunk_start_earthpos(self) -> list[float]:
        # Returns the position of the
        # Earth at the beginning of the chunk in galactic coordinates
        # (cartesian, in meters) as a length 3 list.
        # earth_pos = get_body("earth", Time(self.chunk_start_time.mjd,
        #                                    format="mjd")).transform_to(
        #                                        HeliocentricMeanEcliptic)
        # earth_pos = earth_pos.cartesian.xyz.to(u.AU).transpose()
        # return earth_pos
        return [0,0,0] # TEMP

    def get_chunk_end_earthpos(self) -> list[float]:
        # Returns the position of the
        # Earth at the end of the chunk in galactic coordinates (cartesian,
        # in meters) as a length 3 list. 
        # earth_pos = get_body("earth", Time(self.chunk_end_time.mjd,
        #                                    format="mjd")).transform_to(
        #                                        HeliocentricMeanEcliptic)
        # earth_pos = earth_pos.cartesian.xyz.to(u.AU).transpose()
        # return earth_pos
        return [0,0,0] # TEMP

    def get_chunk_satvel(self) -> list[float]: 
        # Returns the velocity of the satellite
        # at the beginning of the chunk in m/s as a length 3 list.
        
        return [0, 0, 0] # TEMP

    def get_should_compress_tods(self) -> bool: 
        # Returns True if TODs should be compressed; False otherwise. 
        # If the TOD is given as an int, this should typically be True, and False if it is a float.

        # The currently used IRAS TODs are "deplated" level 1 data, containing floats
        return False 

    def get_gain(self, band:str, detector:str) -> float: 
        # Returns the gain of the specified band and detector.

        # Not sure if relevant 
        # Using same value as AKARI script for now 
        return 1

    def get_alpha(self, band:str, detector:str) -> float: 
        # Returns the alpha (correlated noise spectral index) of the specified band and detector.

        # Not sure if relevant 
        # Using same value as AKARI script for now 
        return -1 

    def get_fknee(self, band:str, detector:str) -> float: 
        # Returns the fknee of the specified band and detector.

        # Not sure if relevant 
        # Using same value as AKARI script for now 
        return  1 / (10 * 60)

    def get_sigma0(self, band:str, detector:str) -> float: 
        # Returns the sigma0 (white noise level) of the specified band and detector.

        # Not sure if relevant 
        # Using same value as AKARI script for now 
        return 1 


if __name__ == '__main__':
    data_dir = "/mn/stornext/d5/data/vetleav/IRAS/tod_data/IPAC_level1"
    bdets = BAND_DETS
    bdets = {}
    for band, dets in BAND_DETS.items():
        bdets[band] = dets

    ICDA = IRASCommanderDataAdapter(
        data_dir        = data_dir,
        nside           = 256,
        version         = 0,
        bands           = BANDS,
        band_dets       = bdets,
    )
