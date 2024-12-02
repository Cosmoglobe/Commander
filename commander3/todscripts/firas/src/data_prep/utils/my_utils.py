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
    # filter out channels with non-NaN values
    valid_values = [
        row[f"{variable}_{channel}"]
        for channel in channels
        if not math.isnan(row[f"{variable}_{channel}"])
    ]

    # if there are no valid channels, return None
    if not valid_values:
        return None

    # check if all valid values are the same
    if all(value == valid_values[0] for value in valid_values):
        return valid_values[0]
    else:
        return None


def get_temperature_hl(row, element):
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

    for line in lines:
        if element_search in line.lower():
            side = line.split(" - ")[1]
            if "A" in side:
                trantemp_a = float(line.split(",")[0])
                tranhwid_a = float(line.split(",")[1].split("!")[0])
            elif "B" in side:
                trantemp_b = float(line.split(",")[0])
                tranhwid_b = float(line.split(",")[1].split("!")[0])

    tranlo_a = trantemp_a - tranhwid_a
    tranhi_a = trantemp_a + tranhwid_a
    tranlo_b = trantemp_b - tranhwid_b
    tranhi_b = trantemp_b + tranhwid_b

    if row[f"a_lo_{element}"] < tranlo_a and row[f"a_hi_{element}"] < tranlo_a:
        temp_a = row[f"a_lo_{element}"]
    elif row[f"a_hi_{element}"] > tranhi_a and row[f"a_lo_{element}"] > tranhi_a:
        temp_a = row[f"a_hi_{element}"]
    else:  # TODO: decide what to do if the temperature is within the uncertainty range
        temp_a = np.mean([row[f"a_lo_{element}"], row[f"a_hi_{element}"]])

    if (
        "xcal" != element_search and "collimator" not in element
    ):  # these on the b side are broken
        if row[f"b_lo_{element}"] < tranlo_b and row[f"b_hi_{element}"] < tranlo_b:
            temp_b = row[f"b_lo_{element}"]
        elif row[f"b_hi_{element}"] > tranhi_b and row[f"b_lo_{element}"] > tranhi_b:
            temp_b = row[f"b_hi_{element}"]
        else:  # TODO: decide what to do if the temperature is within the uncertainty range
            temp_b = np.mean([row[f"b_lo_{element}"], row[f"b_hi_{element}"]])

    # TODO: decide how to weight sides, fit for weights?
    if "xcal" != element_search or "collimator" in element:
        return temp_a
    elif "ical" in element:
        return temp_a * 0.1 + temp_b * 0.9
    else:
        return np.mean([temp_a, temp_b])
