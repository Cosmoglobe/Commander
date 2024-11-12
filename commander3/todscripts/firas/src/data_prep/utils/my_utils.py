import math
from datetime import datetime, timedelta

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
