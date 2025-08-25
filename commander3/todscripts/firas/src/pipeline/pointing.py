import numpy as np
import healpy as hp
import scipy


def tune_pointing(gal_lon, gal_lat, gmt, mtm_length, mtm_speed, offset=0):
    """
    Interpolate the pointing to different offsets depending on the operating mode.
    Doing this on the main script for now because it is mostly an experiment and doesn't need a lot of precision.
    """
    times = {"00": 55.36, "10": 44.92, "01": 39.36, "11": 31.76}

    gal_vec = hp.pixelfunc.ang2vec(gal_lon, gal_lat, lonlat=True)
    if offset == 0:
        return gal_vec
    new_gal_vec = np.empty_like(gal_vec)
    # gmt = np.array(gmt).timestamp().values

    gmt = np.array(gmt, dtype="datetime64[s]")
    gmt = (gmt - np.datetime64("1970-01-01T00:00:00Z")) / np.timedelta64(1, "s")

    for i in range(len(gmt)):
        target_time = (
            gmt[i] + offset * times[f"{int(mtm_length[i])}{int(mtm_speed[i])}"]
        )

        # using two minutes around the target time
        lower_bound = target_time - 120
        upper_bound = target_time + 120

        start_id = np.searchsorted(gmt, lower_bound, side="left")
        end_id = np.searchsorted(gmt, upper_bound, side="right")
        idx = np.arange(start_id, end_id)

        if len(idx) > 1 and gmt[start_id] <= target_time <= gmt[end_id - 1]:
            new_gal_vec[i] = scipy.interpolate.interp1d(
                gmt[idx], gal_vec[idx], axis=0, kind="linear"
            )(target_time)
        else:
            new_gal_vec[i] = np.nan

    return new_gal_vec
