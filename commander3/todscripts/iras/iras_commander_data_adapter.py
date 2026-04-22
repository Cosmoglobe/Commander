from cosmoglobe.tod_tools import CommanderDataAdapter
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np 
import healpy
import functools
import glob

from iras_native_tod_reader import IRASTODReader

from astropy.time import Time, TimeDelta
from astropy.coordinates import HeliocentricMeanEcliptic, get_body
from astropy.coordinates import SkyCoord, FK4
GC = SkyCoord(l=0*u.degree, b=0*u.degree, frame='galactic')

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
!!!!!!! TBD !!!!!!!


-get_chunk_data             Yes
-get_chunk_end_earth        Yes (temp)
-get_chunk_end_satpos       Using earth
-get_chunk_end_time         Yes
-get_chunk_indices          Yes
-get_chunk_satvel           Using zeros
-get_chunk_start_earthpos   Yes (temp)
-get_chunk_start_satpos     Using earth
-get_chunk_start_time       Yes
-set_chunk_index            Yes

TESTING: 
Currently testing the script. Below is a list of what to look after at the moment:
- Chunk indices
    - currently using chunks w/ indices (0, N-1)
    - Should probably be (1, N)
    - Change 'get_chunk_indices' if this is the case. 
    - Change file list indexing accordingly. 
- flags
    - currently copy-paste
- UPDATE
- Main issue:
    - earth positions
    - fix those, redo som stuff and crack on.

    FIX DICTIONARY FUCKUPS!!!
    
"""

BANDS = [
    '012',
    '025',
    '060',
    '100'
]

# BAND_DETS = {
BAND_DET_NUMBERS = {
    '012': [23, 24, 25, 26, 27, 28, 29, 30, 47, 48, 49, 50, 51, 52, 53, 54],    #  12 micron / Band 1
    '025': [16, 17, 18, 19, 20, 21, 22, 39, 40, 41, 42, 43, 44, 45, 46],        #  25 micron / Band 2
    '060': [8, 9, 10, 11, 12, 13, 14, 15, 31, 32, 33, 34, 35, 36, 37, 38],      #  60 micron / Band 3
    '100': [1, 2, 3, 4, 5, 6, 7, 55, 56, 57, 58, 59, 60, 61, 62],               # 100 micron / Band 4
}

# ALL_DETS = sorted(sum(BAND_DET_NUMBERS.values(), []))
# import copy
# FULL_BAND_DETS = copy.deepcopy(BAND_DETS)

DEAD_DETECTOR_LIST = [17, 20, 36]
DEGRADED_PERFORMANCE_DETECTORS = [25, 26, 28, 42] # See expl.suppl. Chapter IV.A.2

DEAD_DETECTORS = {}

BAND_DETS = {}
for band in BANDS:
    DEAD_DETECTORS[band] = []
    BAND_DETS[band] = []
    for det_num in BAND_DET_NUMBERS[band]:
        if det_num in DEAD_DETECTOR_LIST:
            DEAD_DETECTORS[band].append(f'IRAS_{band}-{det_num:02d}')
        else:
            BAND_DETS[band].append(f'IRAS_{band}-{det_num:02d}')


NUM_DETECTORS = {key: len(dets) for key, dets in BAND_DETS.items()}

# print(BAND_DETS)
# exit()

### FITS 2 output hdf. TBD
def fits2output_formatter(): pass # PLACEHOLDER
# # This is a callback function to be used with the AKARITODReader. It formats the output data from
# # the TOD reader in a format that is useful for the HDF file generation.
# def fits2output_formatter(file, start_index, end_index, band,
#                           start_det_inds=START_DET_INDS,
#                           num_detectors=NUM_DETECTORS, band_dets=BAND_DETS):
#     """Formats the data in an AKARI fits file in the way we need for the Commander HDFs.

#     Arguments:
#         file (fitsio.FITS): an open fitsio.FITS file instance, pointing to a TOD fits file.
#         start_index, end_index (int): The indices (relative to the start of the file) to fetch.
#         band (str): The band in question (is needed because each fits file contains data for two
#                     bands)

#     Returns:
#         dictionary containing the data we want from each fits file.
#     """
#     out_data = {}
#     # There is an implicit assumption here that we only want contiguous detectors. This seems to be
#     # an okay assumption for now, but I just note it here.
#     start_det_ind = start_det_inds[band]
#     end_det_ind = start_det_ind + num_detectors[band]
#     tot_flux = file[1]['FLUX'].read()[start_index:end_index, start_det_ind:end_det_ind]
#     tot_pixflag = file[1]['PIX_FLAG'].read()[start_index:end_index,
#                                              start_det_ind:end_det_ind, :].astype(bool)
#     for local_det_idx, detname in enumerate(band_dets[band]):
#         out_data[f'{detname}/tod'] = tot_flux[:, local_det_idx]
#         out_data[f'{detname}/pixel_flag'] = tot_pixflag[:, local_det_idx, :]
#     out_data['status_flag'] = file[1]['STATUS'].read()[start_index:end_index].astype(bool)
#     out_data['frame_flag'] = file[1]['FLAG'].read()[start_index:end_index].astype(bool)
#     out_data['aftime'] = file[1]['AFTIME'].read()[start_index:end_index]
#     out_data['ra'] = file[5]['RA'].read()[start_index:end_index]
#     out_data['dec'] = file[5]['DEC'].read()[start_index:end_index]
#     out_data['packet_id'] = file[1]['PACKETID'].read()[start_index:end_index].astype(int)
#     out_data['start_satpos'] = 0
#     out_data['end_satpos'] = 0
#     out_data['start_earthpos'] = 0
#     out_data['end_earthpos'] = 0

#     return out_data



class IRASCommanderDataAdapter(CommanderDataAdapter):
    """
    An IRAS-specific implementation of the CommanderDataAdapter 

    All non-underlined functions are defined as in that class.  
    """
    def __init__(
            self,
            data_dir: str, 
            nside: int,
            num_files_per_segment: int = 100, 
            bands: list[str] = BANDS,
            band_dets: dict[str, list[int]] = BAND_DETS,
    ):
        
        self.data_dir       = data_dir
        self.nside          = nside
        self.bands          = bands
        self.band_dets      = band_dets
        self.num_detectors  = {key: len(dets) for key, dets in band_dets.items()} 
        self.num_files_per_segment = num_files_per_segment
        
        self.fsamp = {
            '012': 16.0,  # [Hz],  12 micron
            '025': 16.0,  # [Hz],  25 micron
            '060': 8.0,   # [Hz],  60 micron
            '100': 4.0,   # [Hz], 100 micron
        }

        self.t0 = Time("1981-01-01T00:00:00", scale="utc") # IRAS's t=0

        # All (alive) detectors have the same number of files
        # Using det 30 to compute the number.    
        self.num_chunks     = len(glob.glob(f"{self.data_dir}/*/*/det30_continuous.npz"))
        self.num_segments   = int(np.ceil(self.num_chunks / self.num_files_per_segment))

        # self.chunk_idx = 50
        # chunks = self.get_chunk_indices(bands[0], segment=1)
        # print(f"{chunks=}")
        # print(f"{len(chunks)=}")
        # print(f"{type(chunks)=}")
        # print(f"{chunks=}")
        # exit()

        self.todreaders = {}
        self.filelist = {}
        for band in self.bands:
            self.todreaders[band] = IRASTODReader(
                data_dir, 
                band, 
                band_dets = self.band_dets
                )
            if self.todreaders[band].testing:
                raise RuntimeError(
                    f"Error. Testing must be False when calling IRASTODReader from data adapter class." 
                )
            self.filelist[band] = self.todreaders[band].filelist
            

    '''def _calculate_chunk_files(self):
        """
        UPDATE: MAY BE OBSOLETE
        Internal function that calculates a mapping between chunk index and fits files. 
        This is temporary, as we might want to do something more sophisticated about 
        how we define the chunks.
        """
        chunk_file_map = {}
        for band in self.bands:
            chunk_file_map[band] = {}
            curr_chunk_idx = 0
            for file in self.filelist[band]:
                chunk_file_map[band][curr_chunk_idx] =  []
                chunk_file_map[band][curr_chunk_idx].append(file)
                curr_chunk_idx += 1
        return chunk_file_map'''

    def _process_iras_flags(self):
        pass

    def _process_iras_tod_reader_data(self):
        """Internal function to process the data coming from the IRAS TOD reader and getting it
        into the format needed for the HDFwriter"""

        self.todreaders = None
        self._process_akari_flags = self._process_iras_flags


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
                todreader_data[band][det]['flag'] = self._process_akari_flags(
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

    def get_num_segments(self, band: str)  ->  int: 
        # For a given band, returns the
        # number of segments (i.e. number of HDF files) for that band.

        ## IRAS:
        # Number of npz files for a single detector:
        #   - 5787 (for deplated Level 1 data)
        return self.num_segments

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

        start_idx   = (segment-1) * self.num_files_per_segment
        end_idx     = segment * self.num_files_per_segment
        # if end_idx > self.num_chunks[band]:
            # end_idx = self.num_chunks[band]
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
        self.chunk_idx = chunk_idx
        self.all_chunk_data = None #Reset chunk


    '''def get_chunk(func):
        """Decorator function to make sure that a chunk is always preloaded"""
        print("NOOOOOOOO")
        @functools.wraps(func)
        def inner(self, *args):
            if self.all_chunk_data is None:
                self.all_chunk_data = self._process_iras_tod_reader_data()
            return func(self, *args)

        return inner'''

    def get_chunk_data(self, band:str, detector:str): 
        # -> tod(np.array[float/int]), pix(np.array[int]), psi(np.array(float)), flag(np.array[int]):
        # Given a band and detector (with the chunk number assumed set by
        # set_chunk_index), returns the tod, pixel, psi, and flag array
        # corresponding to the detector and chunk in question. The flag array is
        # an array of integers representing bits.

        fname = self.filelist[band][detector][self.chunk_idx - 1]
        # file = np.load(f"{self.data_dir}/sop534/obs013/det30_continuous.npz")

        ### KEYS: t, ra, dec, flux, n_merged, dedup_seconds, sop, obs, det
        data = np.load(fname)

        t = data["t"]
        self.chunk_start_time   = t[0]  * u.s + self.t0
        self.chunk_end_time     = t[-1] * u.s + self.t0

        tod = data["flux"] * self.todreaders[band].flux_to_MJy_sr[detector] # MJy/sr
        ra  = data["ra"]
        dec = data["dec"]

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

        coords = sc.transform_to(GC)
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

        # TEMP
        # return self.get_chunk_start_earthpos()
        return [0,0,0]


    def get_chunk_end_satpos(self) -> list[float]:
        # Returns the position of the
        # satellite at the end of the chunk in galactic coordinates (cartesian,
        # in meters) as a length 3 list. 

        # TEMP
        # return self.get_chunk_end_earth()
        return [0,0,0]


    def get_chunk_start_earthpos(self) -> list[float]:
        # Returns the position of the
        # Earth at the beginning of the chunk in galactic coordinates
        # (cartesian, in meters) as a length 3 list.
        # earth_pos = get_body("earth", Time(self.chunk_start_time.mjd,
        #                                    format="mjd")).transform_to(
        #                                        HeliocentricMeanEcliptic)
        # earth_pos = earth_pos.cartesian.xyz.to(u.AU).transpose()
        # return earth_pos
        return [0,0,0]


    def get_chunk_end_earthpos(self) -> list[float]:
        # Returns the position of the
        # Earth at the end of the chunk in galactic coordinates (cartesian,
        # in meters) as a length 3 list. 
        # earth_pos = get_body("earth", Time(self.chunk_end_time.mjd,
        #                                    format="mjd")).transform_to(
        #                                        HeliocentricMeanEcliptic)
        # earth_pos = earth_pos.cartesian.xyz.to(u.AU).transpose()
        # return earth_pos
        return [0,0,0]

    def get_chunk_satvel(self) -> list[float]: 
        # Returns the velocity of the satellite
        # at the beginning of the chunk in m/s as a length 3 list.
        return [0, 0, 0]

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
        bdets[band] = dets[::3]

    ICDA = IRASCommanderDataAdapter(
        data_dir        = data_dir,
        nside           = 256,
        bands           = [BANDS[1]],
        band_dets       = bdets,
    )
