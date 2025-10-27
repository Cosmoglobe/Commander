import math
from datetime import datetime, timedelta

import numpy as np

channels = ["lh", "ll", "rh", "rl"]


def parse_date_string(gmt_nb, include_s=True):
    # gmt_nb = str(gmt_nb)
    # print(gmt_nb)

    # Split the string into components
    year = int(gmt_nb[:2])
    day_of_year = int(gmt_nb[2:5])
    hour = int(gmt_nb[5:7])
    minute = int(gmt_nb[7:9])
    # second = round(float(f"{gmt_nb[9:11]}.{gmt_nb[11:]}"))
    if include_s:
        second = int(gmt_nb[9:11])
        millisecond = int(gmt_nb[11:])

    # if second == 60:
    #     second = 0
    #     minute += 1
    # if minute == 60:
    #     minute = 0
    #     hour += 1
    # if hour == 24:
    #     hour = 0
    #     day_of_year += 1

    # Compute the date from the year and day of the year
    base_date = datetime(year=year + 1900, month=1, day=1)
    date = base_date + timedelta(days=day_of_year - 1)

    # Construct the complete datetime with time components
    if not include_s:
        second = 0
        millisecond = 0
    final_date = datetime(
        year=date.year,
        month=date.month,
        day=date.day,
        hour=hour,
        minute=minute,
        second=second,
        microsecond=millisecond * 1000,
    )

    return final_date


def clean_variable(row, variable):
    """
    Function to clean up variables so that instead of 4 of each variable (one for each channel), we have only one variable after checking that the four channels are all the same (except for the ones that are NaNs).
    """
    values = np.array([row[f"{variable}_{channel}"] for channel in channels])

    # check if all values are NaN
    if np.all(np.isnan(values)):
        return np.nan

    # check if all valid values are the same
    if np.array_equal(values, values, equal_nan=True):
        return values[0]
    else:
        return np.nan


def get_temperature_hl(row, element, side, name_search=None):
    """
    Checks fex_grttrans.txt for which high/low current temperature values to use and returns the weighted temperatures.
    """
    if (
        "xcal" in element
    ):  # for now following what they did with the XCAL, where they weight only the cone GRT. TODO: change this?
        element_search = "xcal s"
        element = "xcal_cone"
    else:
        element_search = element

    # identify which component is passed on and match to line in fex_grttrans.txt
    with open("../reference/fex_grttrans.txt") as f:
        lines = f.readlines()

    if name_search is not None:
        element_search = name_search

    if element == "refhorn" and side == "a":
        temp = row[f"{side}_hi_{element}"]
    if (element == "refhorn" or element == "skyhorn") and side == "b":
        temp = row[f"{side}_lo_{element}"]
    else:
        for line in lines:
            if element_search in line.lower():
                side_search = line.split(" - ")[1]
                if side.upper() in side_search:
                    trantemp = float(line.split(",")[0])
                    tranhwid = float(line.split(",")[1].split("!")[0])

        tranlo = trantemp - tranhwid
        tranhi = trantemp + tranhwid

        if (
            row[f"{side}_lo_{element}"] < tranlo
            and row[f"{side}_hi_{element}"] < tranlo
        ):
            temp = row[f"{side}_lo_{element}"]
        elif (
            row[f"{side}_hi_{element}"] > tranhi
            and row[f"{side}_lo_{element}"] > tranhi
        ):
            temp = row[f"{side}_hi_{element}"]
        else:  # TODO: decide what to do if the temperature is within the uncertainty range
            temp = np.mean([row[f"a_lo_{element}"], row[f"a_hi_{element}"]])

    return temp


# Cache for transition temperatures to avoid repeated file I/O
_TRANSITION_TEMPS = None


def _load_transition_temperatures():
    """
    Load transition temperatures from fex_grttrans.txt once and cache them.
    Returns a dictionary with keys like ('ical', 'a') -> (trantemp, tranhwid)
    """
    global _TRANSITION_TEMPS
    if _TRANSITION_TEMPS is not None:
        return _TRANSITION_TEMPS

    _TRANSITION_TEMPS = {}

    with open("../reference/fex_grttrans.txt") as f:
        lines = f.readlines()

    # Map element names to search strings
    element_map = {
        "ical": "ical grt",
        "dihedral": "dihedral grt",
        "refhorn": "reference horn grt",
        "skyhorn": "skyhorn grt",
        "collimator": "collimator",
        "xcal_cone": "xcal s",
        "bol_assem": "bolometer",
    }

    for line in lines:
        line_lower = line.lower()
        for element_key, search_str in element_map.items():
            if search_str in line_lower:
                # Determine side
                if "a side" in line_lower:
                    side = "a"
                elif "b side" in line_lower:
                    side = "b"
                else:
                    continue

                # Parse temperature and half-width
                parts = line.split(",")
                if len(parts) >= 2:
                    try:
                        trantemp = float(parts[0].strip())
                        tranhwid = float(parts[1].split("!")[0].strip())
                        _TRANSITION_TEMPS[(element_key, side)] = (trantemp, tranhwid)
                    except (ValueError, IndexError):
                        continue

    return _TRANSITION_TEMPS


def get_temperature_hl_vectorized(lo_temps, hi_temps, element, side):
    """
    Vectorized version of get_temperature_hl for numpy arrays.

    Parameters:
    -----------
    lo_temps : np.ndarray
        Array of low-current temperatures
    hi_temps : np.ndarray
        Array of high-current temperatures
    element : str
        Element name ('ical', 'dihedral', 'refhorn', 'skyhorn', 'collimator', 'xcal_cone')
    side : str
        Side identifier ('a' or 'b')

    Returns:
    --------
    np.ndarray
        Selected temperatures based on transition logic
    """
    # Handle special cases first
    if element == "refhorn" and side == "a":
        return hi_temps.copy()
    if (element == "refhorn" or element == "skyhorn") and side == "b":
        return lo_temps.copy()

    # Load transition temperatures
    trans_temps = _load_transition_temperatures()

    # Map 'xcal' to 'xcal_cone' if needed
    element_key = "xcal_cone" if "xcal" in element else element

    # Get transition temperature and half-width
    if (element_key, side) not in trans_temps:
        # If no transition data, return mean
        return (lo_temps + hi_temps) / 2.0

    trantemp, tranhwid = trans_temps[(element_key, side)]
    tranlo = trantemp - tranhwid
    tranhi = trantemp + tranhwid

    # Vectorized selection logic using np.select
    conditions = [
        (lo_temps < tranlo) & (hi_temps < tranlo),  # Both below transition -> use lo
        (hi_temps > tranhi) & (lo_temps > tranhi),  # Both above transition -> use hi
    ]

    choices = [
        lo_temps,
        hi_temps,
    ]

    # Default: within transition range -> use mean
    default = (lo_temps + hi_temps) / 2.0

    return np.select(conditions, choices, default=default)


def convert_gain(row, channel):
    if (
        row[f"gain_{channel}"] == 0
        or row[f"gain_{channel}"] == 1
        or row[f"gain_{channel}"] == 2
        or row[f"gain_{channel}"] == 3
        or row[f"gain_{channel}"] == 4
        or row[f"gain_{channel}"] == 5
        or row[f"gain_{channel}"] == 6
        or row[f"gain_{channel}"] == 7
    ):
        conv = {0: 1, 1: 3, 2: 10, 3: 30, 4: 100, 5: 300, 6: 1000, 7: 3000}

        return conv[row[f"gain_{channel}"]]
    else:
        return np.nan


def scan_up_down(lat):
    """
    Function to determine for each point in the dataframe if the instrument is scanning up or down based on the terrestrial latitude.
    """
    diff = np.diff(lat)
    up = diff > 0
    down = ~up

    scan = np.zeros_like(up, dtype=int)
    scan[up] = 1
    scan[down] = -1
    scan = np.append(scan, scan[-1])  # last value is the same as the second to last

    return scan


def binary_to_gmt(binary):
    """
    Converts input ADT time to 14-element string. Adapted from Nathan's pipeline.
    """

    # ADT is a 64-bit quadword containing the number of 100-nanosecond ticks
    # since November 17, 1858

    epoch = np.datetime64("1858-11-17T00:00:00")
    # Convert binary to microseconds and add to epoch
    microseconds = (0.1 * binary).astype("timedelta64[us]")
    return epoch + microseconds
