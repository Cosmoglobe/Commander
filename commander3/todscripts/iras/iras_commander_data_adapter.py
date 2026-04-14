from cosmoglobe.tod_tools import CommanderDataAdapter
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np 
import healpy
from astropy.time import Time, TimeDelta



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
!!!!!!! TBD !!!!!!!
List of methods to be implemented:
-get_alpha
-get_bands
-get_chunk_data
-get_chunk_end_earth
-get_chunk_end_satpos
-get_chunk_end_time
-get_chunk_indices
-get_chunk_satvel
-get_chunk_start_earthpos
-get_chunk_start_satpos
-get_chunk_start_time
-get_detector_names
-get_experiment_name
-get_fknee
-get_fsamp
-get_gain
-get_mbangs
-get_npsi
-get_nside
-get_num_segments
-get_polangs
-get_should_compress_tods
-get_sigma0
-get_version
-set_chunk_index
"""

BAND1_DETS = [23, 24, 25, 26, 27, 28, 29, 30, 47, 48, 49, 50, 51, 52, 53, 54]
BAND2_DETS = [16, 17, 18, 19, 20, 21, 22, 39, 40, 41, 42, 43, 44, 45, 46]
BAND3_DETS = [8, 9, 10, 11, 12, 13, 14, 15, 31, 32, 33, 34, 35, 36, 37, 38]
BAND4_DETS = [1, 2, 3, 4, 5, 6, 7, 55, 56, 57, 58, 59, 60, 61, 62]

ALL_DETS = sorted(BAND1_DETS + BAND2_DETS + BAND3_DETS + BAND4_DETS)

DEAD_DETECTORS = [17, 20, 36]
DEGRADED_PERFORMANCE_DETECTORS = [25, 26, 28, 42] # See expl.suppl. Chapter IV.A.2

BANDS = [
    '012',
    '025',
    '060',
    '100'
]

NUM_DETECTORS = {
    '012': 0, # 12 micron
    '025': 0, # 12 micron
    '060': 0, # 12 micron
    '100': 0, # 12 micron
}


class IRASCommanderDataAdapter(CommanderDataAdapter):
    """
    An IRAS-specific implementation of the CommanderDataAdapter 

    All non-underlined functions are defined as in that class.  
    """
    def __init__(self):
        pass 


    def get_num_segments(self, band: str)  ->  int: 
        # For a given band, returns the
        # number of segments (i.e. number of HDF files) for that band.
        pass

    def get_version(self) -> int: 
        #Returns the current version of the HDF generation.
        pass

    def get_experiment_name(self) -> str: 
        # Returns the name of the experiment.
        pass 

    def get_npsi(self) -> int: 
        # Returns the number of psi bins to use for the psi compression.
        pass 
    
    def get_nside(self, band: str) -> int: 
        # Given a band name, returns an integer specifying
        # the nside of the pixelization used for the pointing.
        pass

    def get_polangs(self, band: str) -> list[float]: 
        # Given a band, returns a list of
        # length ndet which gives the polarization angles for each detector of
        # that band.
        pass

    def get_mbangs(self, band: str) -> list[float]: 
        # Given a band name, returns a list
        # of length ndet which gives the main beam angle of that each detector of
        # that band.
        pass

    def get_fsamp(self, band:str) -> float: 
        # Given a band name, returns the sampling
        # rate of that band.
        pass

    def get_detector_names(self, band:str) -> list[str]: 
        # Given a band name, returns a
        # list of strings specifying the detector names of that band - these will
        # be used as fields in the output HDF file.
        pass

    def get_chunk_indices(self, band:str, segment:int) -> list[int]: 
        # Given a band
        # name and a segment number belonging to that band, gives the indices of
        # the chunks (previously:scans) that belong to the segment in question.
        pass

    def get_bands(self) -> list[str]: 
        # Returns a list of the band names.
        pass

    def set_chunk_index(self, chunk_idx:int): 
        # Since it can be bad,
        # performance-wise, to have to specify the chunk for every low-level data
        # fetch operation, this function allows the reader to signal to the
        # adapter which chunk it is going to process next. After this, low-level
        # data fetch calls will *not* be specifying the chunk number, but will
        # implicitly assume that the adapter knows which chunk number is the
        # current one.
        pass

    def get_chunk_data(self, band:str, detector:str): 
        # -> tod(np.array[float/int]), pix(np.array[int]), psi(np.array(float)), flag(np.array[int]):
        # Given a band and detector (with the chunk number assumed set by
        # set_chunk_index), returns the tod, pixel, psi, and flag array
        # corresponding to the detector and chunk in question. The flag array is
        # an array of integers representing bits.
        pass

    def get_chunk_start_time(self): # -> astropy.time.Time: 
        # Returns the time of the
        # beginning of the chunk.
        pass

    def get_chunk_end_time(self): # -> astropy.time.Time: 
        # Returns the time of the end of
        # the chunk.
        pass

    def get_chunk_start_satpos(self) -> list[float]:
        # Returns the position of the
        # satellite at the beginning of the chunk in galactic coordinates
        # (cartesian, in meters) as a length 3 list.
        pass

    def get_chunk_end_satpos(self) -> list[float]:
        # Returns the position of the
        # satellite at the end of the chunk in galactic coordinates (cartesian,
        # in meters) as a length 3 list. 
        pass

    def get_chunk_start_earthpos(self) -> list[float]:
        # Returns the position of the
        # Earth at the beginning of the chunk in galactic coordinates
        # (cartesian, in meters) as a length 3 list.
        pass

    def get_chunk_end_earth(self) -> list[float]:
        # Returns the position of the
        # Earth at the end of the chunk in galactic coordinates (cartesian,
        # in meters) as a length 3 list. 
        pass

    def get_chunk_satvel(self) -> list[float]: 
        # Returns the velocity of the satellite
        # at the beginning of the chunk in m/s as a length 3 list.
        pass

    def get_should_compress_tods(self) -> bool: 
        # Returns True if TODs should be
        # compressed; False otherwise. If the TOD is given as an int, this should
        # typically be True, and False if it is a float.
        pass

    def get_gain(self, band:str, detector:str) -> float: 
        # Returns the gain of the
        # specified band and detector.
        pass

    def get_alpha(self, band:str, detector:str) -> float: 
        # Returns the alpha
        # (correlated noise spectral index) of the specified band and detector.
        pass

    def get_fknee(self, band:str, detector:str) -> float: 
        # Returns the fknee of the
        # specified band and detector.
        pass

    def get_sigma0(self, band:str, detector:str) -> float: 
        # Returns the sigma0 (white
        # noise level) of the specified band and detector.
        pass 

ICDA = IRASCommanderDataAdapter()

