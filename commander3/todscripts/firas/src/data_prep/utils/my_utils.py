import math
from datetime import datetime, timedelta

import numpy as np

channels = ["lh", "ll", "rh", "rl"]


def parse_date_string(gmt_nb):
    # gmt_nb = str(gmt_nb)
    # print(gmt_nb)

    # Split the string into components
    year = int(gmt_nb[:2])
    day_of_year = int(gmt_nb[2:5])
    hour = int(gmt_nb[5:7])
    minute = int(gmt_nb[7:9])
    # second = round(float(f"{gmt_nb[9:11]}.{gmt_nb[11:]}"))
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
    with open("../../reference/fex_grttrans.txt") as f:
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


def convert_gain(row):
    if (
        row[f"gain"] == 0
        or row[f"gain"] == 1
        or row[f"gain"] == 2
        or row[f"gain"] == 3
        or row[f"gain"] == 4
        or row["gain"] == 5
        or row["gain"] == 6
        or row["gain"] == 7
    ):
        conv = {0: 1, 1: 3, 2: 10, 3: 30, 4: 100, 5: 300, 6: 1000, 7: 3000}

        return conv[row[f"gain"]]
    else:
        return np.nan


def scan_up_down(lat):
    """
    Function to determine for each point in the dataframe if the instrument is scanning up or down based on the terrestrial latitude.
    """
    up = np.empty_like(lat)
    for i in range(len(lat) - 1):
        if lat[i + 1] > lat[i]:
            up[i] = 1
        elif lat[i + 1] < lat[i]:
            up[i] = 0
        else:
            up[i] = up[i - 1]
    up[-1] = up[-2]
    return up
