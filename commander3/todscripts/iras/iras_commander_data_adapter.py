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
-get_chunk_data
-get_chunk_end_earth
-get_chunk_end_satpos
-get_chunk_end_time
-get_chunk_indices
-get_chunk_satvel
-get_chunk_start_earthpos
-get_chunk_start_satpos
-get_chunk_start_time
-get_num_segments
-set_chunk_index
"""

BANDS = [
    '012',
    '025',
    '060',
    '100'
]

BAND_DETS = {
    '012': [23, 24, 25, 26, 27, 28, 29, 30, 47, 48, 49, 50, 51, 52, 53, 54],    #  12 micron / Band 1
    '025': [16, 17, 18, 19, 20, 21, 22, 39, 40, 41, 42, 43, 44, 45, 46],        #  25 micron / Band 2
    '060': [8, 9, 10, 11, 12, 13, 14, 15, 31, 32, 33, 34, 35, 36, 37, 38],      #  60 micron / Band 3
    '100': [1, 2, 3, 4, 5, 6, 7, 55, 56, 57, 58, 59, 60, 61, 62],               # 100 micron / Band 4
}
ALL_DETS = sorted(sum(BAND_DETS.values(), []))
# import copy
# FULL_BAND_DETS = copy.deepcopy(BAND_DETS)

DEAD_DETECTOR_LIST = [17, 20, 36]
DEGRADED_PERFORMANCE_DETECTORS = [25, 26, 28, 42] # See expl.suppl. Chapter IV.A.2

DEAD_DETECTORS = {}

# for band, band_dets in BAND_DETS.items():
for band in BANDS:
    dead_dets_band = []
    for dead_det in DEAD_DETECTOR_LIST:
        # if dead_det in band_dets:
        if dead_det in BAND_DETS[band]:
            dead_dets_band.append(dead_det)
            BAND_DETS[band].remove(dead_det)
    DEAD_DETECTORS[band] = dead_dets_band

NUM_DETECTORS = {key: len(dets) for key, dets in BAND_DETS.items()}




# for k, v in BAND_DETS.items():
#     print(k, f"ndet={len(v)}", v)
# print()
# print(NUM_DETECTORS)
# exit()
class IRASCommanderDataAdapter(CommanderDataAdapter):
    """
    An IRAS-specific implementation of the CommanderDataAdapter 

    All non-underlined functions are defined as in that class.  
    """
    def __init__(
            self,
            # data_dir: str | os.PathLike, 
            data_dir: str, 
            nside: int,
            bands: list[str] = BANDS,
            num_detectors: dict[str, int] = NUM_DETECTORS,
            band_dets: dict[str, list[int]] = BAND_DETS,
    ):
        
        self.nside          = nside
        self.bands          = bands
        self.num_detectors  = num_detectors 
        self.band_dets      = band_dets
        

        self.fsamp = {
            '012': 16.0,  # [Hz],  12 micron
            '025': 16.0,  # [Hz],  25 micron
            '060': 8.0,   # [Hz],  60 micron
            '100': 4.0,   # [Hz], 100 micron
        }

        self.detectors = {}

        for band in self.bands:
            self.detectors[band] = []
            for det_num in self.band_dets[band]:
                self.detectors[band].append(f'IRAS_{band}-{det_num:02d}')
            

    def get_num_segments(self, band: str)  ->  int: 
        # For a given band, returns the
        # number of segments (i.e. number of HDF files) for that band.
        pass

    def get_version(self) -> int: 
        #Returns the current version of the HDF generation.

        # Version 1:
        #   - Deplated level1 data from IPAC, in units of W/m^2.
        return 1 

    def get_experiment_name(self) -> str: 
        return "IRAS"

    def get_npsi(self) -> int: 
        # Returns the number of psi bins to use for the psi compression.

        # Taken from "/home/vetle/Documents/PhD_astro/Cosmoglobe/cosmoglobe/tod_tools/iras/comm_iras_data_adapter.py"
        # Line 31: self.npsi = 8  # For IRAS this shouldn't matter
        # WILL/SHOULD CHECK LATER 
        return 8
    
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
        return self.detectors[band]

    def get_chunk_indices(self, band:str, segment:int) -> list[int]: 
        # Given a band
        # name and a segment number belonging to that band, gives the indices of
        # the chunks (previously:scans) that belong to the segment in question.
        pass

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

ICDA = IRASCommanderDataAdapter(
    data_dir        = None,
    nside           = 256,
    # bands           = ,
    # num_detectors   = ,
    # band_dets       = ,
)
