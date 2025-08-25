"""
Data processing script has to be re-written again, not matching eng_time/time now,
and instead interpolating over gmt to get the engineering data corresponding to the science data.

Run directly: python data_prep_v3.py
"""

import time

import h5py
import healpy as hp
import numpy as np
import pandas as pd
import tables as tb
from scipy import interpolate
from utils.my_utils import (
    convert_gain,
    get_temperature_hl,
    parse_date_string,
    scan_up_down,
)

# check how much time the script takes to run
start_time = time.time()

# OPENING ORIGINAL DATA FILES
fdq_sdf = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5")
fdq_eng = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5")

channels = {"lh": 2, "ll": 3, "rh": 0, "rl": 1}

# PARSING THE H5 DATA INTO A PANDAS DATAFRAME TO MANIPULATE DATA EASIER
print("decoding the gmt data")
gmt = {}
for channel in channels:
    gmt[channel] = np.array(fdq_sdf[f"fdq_sdf_{channel}/ct_head/gmt"][()]).astype(
        str
    )  # .astype(int)

gmt_parsed = {"lh": [], "ll": [], "rh": [], "rl": []}

for channel in channels:
    for gmt_nb in gmt[channel]:
        gmt_parsed[channel].append(parse_date_string(gmt_nb))

# adding an offset to the gmt_parsed data
for channel in channels:
    gmt_parsed[channel] = pd.to_datetime(gmt_parsed[channel]) - pd.Timedelta(seconds=60)

ifg = {}
xcal_pos = {}
mtm_length = {}
mtm_speed = {}
fake = {}
upmode = {}
adds_per_group = {}
t = {}
sweeps = {}
gain = {}
sun_angle = {}
earth_limb = {}
moon_angle = {}

NSIDE = 32
pix_terr = {}  # earth coordinates
pix_ecl = {}  # ecliptic coordinates
pix_gal = {}  # galactic coordinates
pix_cel = {}  # celestial coordinates, probably J1950
# Longitudes and latitudes are stored in radians*1e4
fact = 180.0 / np.pi / 1e4
scan = {}

for channel in channels:
    ifg[channel] = fdq_sdf[f"fdq_sdf_{channel}/ifg_data/ifg"][()]
    xcal_pos[channel] = fdq_sdf[f"fdq_sdf_{channel}/dq_data/xcal_pos"][()]
    mtm_length[channel] = fdq_sdf[f"fdq_sdf_{channel}/sci_head/mtm_length"][()]
    mtm_speed[channel] = fdq_sdf[f"fdq_sdf_{channel}/sci_head/mtm_speed"][()]
    fake[channel] = fdq_sdf[f"fdq_sdf_{channel}/dq_data/fake"][()]
    upmode[channel] = fdq_sdf[f"fdq_sdf_{channel}/sci_head/sc_head1a"][()]
    adds_per_group[channel] = fdq_sdf[f"fdq_sdf_{channel}/sci_head/sc_head9"][()]
    # check for nans in adds_per_group
    if np.isnan(adds_per_group[channel]).any():
        print(f"There are nans in adds_per_group (checkpoint 1, channel {channel})")
        # how many nans are there?
        print(np.sum(np.isnan(adds_per_group[channel])))
    else:
        print(f"No nans in adds_per_group (checkpoint 1, channel {channel})")
    t[channel] = fdq_sdf[f"fdq_sdf_{channel}/ct_head/time"][()]
    sweeps[channel] = fdq_sdf[f"fdq_sdf_{channel}/sci_head/sc_head11"][()]
    gain[channel] = fdq_sdf[f"fdq_sdf_{channel}/sci_head/gain"][()]
    sun_angle[channel] = fdq_sdf[f"fdq_sdf_{channel}/attitude/sun_angle"][()]
    earth_limb[channel] = fdq_sdf[f"fdq_sdf_{channel}/attitude/earth_limb"][()]
    moon_angle[channel] = fdq_sdf[f"fdq_sdf_{channel}/attitude/moon_angle"][()]

    lon = fdq_sdf[f"fdq_sdf_{channel}/attitude/terr_longitude"][()] * fact
    lat = fdq_sdf[f"fdq_sdf_{channel}/attitude/terr_latitude"][()] * fact
    scan[channel] = scan_up_down(lat).astype(int)
    pix_terr[channel] = hp.ang2pix(NSIDE, lon, lat, lonlat=True).astype(float)

    lon = fdq_sdf[f"fdq_sdf_{channel}/attitude/ecliptic_longitude"][()] * fact
    lat = fdq_sdf[f"fdq_sdf_{channel}/attitude/ecliptic_latitude"][()] * fact
    pix_ecl[channel] = hp.ang2pix(NSIDE, lon, lat, lonlat=True).astype(float)

    lon = fdq_sdf[f"fdq_sdf_{channel}/attitude/galactic_longitude"][()] * fact
    lat = fdq_sdf[f"fdq_sdf_{channel}/attitude/galactic_latitude"][()] * fact
    pix_gal[channel] = hp.ang2pix(NSIDE, lon, lat, lonlat=True).astype(float)

    lon = fdq_sdf[f"fdq_sdf_{channel}/attitude/ra"][()] * fact
    lat = fdq_sdf[f"fdq_sdf_{channel}/attitude/dec"][()] * fact
    pix_cel[channel] = hp.ang2pix(NSIDE, lon, lat, lonlat=True).astype(float)

print("getting each channel into its own df")

df = {}

for channel in channels:
    print(channel)
    df[channel] = pd.DataFrame(
        {
            "gmt": gmt_parsed[channel],
            "ifg": list(ifg[channel]),
            "xcal_pos": list(xcal_pos[channel]),
            "mtm_length": list(mtm_length[channel]),
            "mtm_speed": list(mtm_speed[channel]),
            "fake": list(fake[channel]),
            "upmode": list(upmode[channel]),
            "adds_per_group": list(adds_per_group[channel]),
            "pix_gal": list(pix_gal[channel]),
            "pix_ecl": list(pix_ecl[channel]),
            "pix_cel": list(pix_cel[channel]),
            "pix_terr": list(pix_terr[channel]),
            "scan": list(scan[channel]),
            "time": list(t[channel]),
            "sweeps": list(sweeps[channel]),
            "gain": list(gain[channel]),
            "sun_angle": list(sun_angle[channel]),
            "earth_limb": list(earth_limb[channel]),
            "moon_angle": list(moon_angle[channel]),
        }
    ).sort_values("gmt")
    print(df[channel].columns)
    # check for nans in adds_per_group
    if np.isnan(df[channel]["adds_per_group"]).any():
        print(f"There are nans in adds_per_group (checkpoint 2, channel {channel})")
        # how many nans are there?
        print(np.sum(np.isnan(df[channel]["adds_per_group"])))
    else:
        print(f"No nans in adds_per_group (checkpoint 2, channel {channel})")


# USING THIS FOR MERGE_ASOF
tolerance = pd.Timedelta(seconds=16)
print("getting all possible gmts so that we can do an outer join using merge_asof")
unified_timestamps = (
    pd.DataFrame(
        pd.concat([df["lh"]["gmt"], df["ll"]["gmt"], df["rh"]["gmt"], df["rl"]["gmt"]])
    )
    .drop_duplicates()
    .sort_values("gmt")
    .reset_index(drop=True)
)

print(f"Timestamps before removing close ones: {unified_timestamps.tail()}")

# Taking initial comparison values from first row
timestamp = unified_timestamps.iloc[0]["gmt"]
# Including first row in result
filters = [True]

# Skipping first row in comparisons
for index, row in unified_timestamps.iloc[1:].iterrows():
    if timestamp - tolerance <= row["gmt"] <= timestamp + tolerance:
        filters.append(False)
    else:
        filters.append(True)
        # Updating values to compare based on latest accepted row
        timestamp = row["gmt"]

non_duplicate_timestamps = (
    unified_timestamps.loc[filters].sort_values("gmt").reset_index(drop=True)
)

print(f"Timestamps after removing close ones: {non_duplicate_timestamps.tail()}")
print(df["lh"].columns)
# outer-join the dataframes on gmt
merged_df = pd.merge_asof(
    non_duplicate_timestamps,
    df["lh"],
    on="gmt",
    direction="nearest",
    tolerance=tolerance,
)

merge_by = [
    "pix_gal",
    "pix_ecl",
    "pix_cel",
    "pix_terr",
    "scan",
    "adds_per_group",
    "gain",
    "sweeps",
    "sun_angle",
    "earth_limb",
    "moon_angle",
    "mtm_length",
    "mtm_speed",
    "fake",
    "upmode",
    "xcal_pos",
]

for channel in channels:
    for merge_key in merge_by:
        df[channel][merge_key] = df[channel][merge_key].astype(float)

merged_df = pd.merge_asof(
    merged_df,
    df["ll"],
    on="gmt",
    direction="nearest",
    suffixes=("_lh", "_ll"),
    tolerance=tolerance,
    by=merge_by,
)
merged_df = pd.merge_asof(
    merged_df,
    df["rh"],
    on="gmt",
    direction="nearest",
    suffixes=("_ll", "_rh"),
    tolerance=tolerance,
    by=merge_by,
)
merged_df = pd.merge_asof(
    merged_df,
    df["rl"],
    on="gmt",
    direction="nearest",
    suffixes=("_rh", "_rl"),
    tolerance=tolerance,
    by=merge_by,
)

# sort by gmt
merged_df = merged_df.sort_values("gmt").reset_index(drop=True)

# drop rows where all ifgs are nans
merged_df = merged_df[
    (merged_df["ifg_lh"].apply(lambda x: np.isnan(x).all() == False))
    | (merged_df["ifg_ll"].apply(lambda x: np.isnan(x).all() == False))
    | (merged_df["ifg_rh"].apply(lambda x: np.isnan(x).all() == False))
    | (merged_df["ifg_rl"].apply(lambda x: np.isnan(x).all() == False))
]

# CLEANING XCAL_POS AND CONSTRAINING FOR ONLY 1 AND 2
# making sure xcal_pos is the same within each record
before = len(merged_df)
print(f"Number of rows removed due to xcal_pos: {before - len(merged_df)}")
merged_df = merged_df[(merged_df["xcal_pos"] == 1) | (merged_df["xcal_pos"] == 2)]

# same thing but for mtm length and speed
before = len(merged_df)
merged_df = merged_df[
    (merged_df["mtm_length"] == 0) | (merged_df["mtm_length"] == 1)
]  # only 0 and 1 for mtm_length
merged_df = merged_df[
    (merged_df["mtm_speed"] == 0) | (merged_df["mtm_speed"] == 1)
]  # only 0 and 1 for mtm_speed

# same thing for fake it flag
before = len(merged_df)
print(f"Number of rows removed due to fake: {before - len(merged_df)}")
# fake-it mode on is 1, so i am guessing that 0 is when it is off
merged_df = merged_df[(merged_df["fake"] == 0)]

# the following variables i don't think necessarily need to be the same in all channels for the record to be valid
# adds per group
merged_df = merged_df[
    (merged_df[f"adds_per_group"] == 1)
    | (merged_df[f"adds_per_group"] == 2)
    | (merged_df[f"adds_per_group"] == 3)
    | (merged_df[f"adds_per_group"] == 8)
    | (merged_df[f"adds_per_group"] == 12)
]
# check for nans in adds_per_group
if np.isnan(merged_df[f"adds_per_group"]).any():
    print(f"There are nans in adds_per_group (checkpoint 3, channel {channel})")
    # how many nans are there?
    print(np.sum(np.isnan(merged_df[f"adds_per_group"])))
else:
    print(f"No nans in adds_per_group (checkpoint 3)")
# sweeps
merged_df = merged_df[
    # (merged_df[f"sweeps_{channel}"] == 1) |
    (merged_df[f"sweeps"] == 4)
    | (merged_df[f"sweeps"] == 16)
]

# gain
merged_df[f"gain"] = merged_df.apply(convert_gain, axis=1)
merged_df = merged_df[
    (merged_df[f"gain"] == 1)
    | (merged_df[f"gain"] == 3)
    | (merged_df[f"gain"] == 10)
    | (merged_df[f"gain"] == 30)
    | (merged_df[f"gain"] == 100)
    | (merged_df[f"gain"] == 300)
    | (merged_df[f"gain"] == 1000)
    | (merged_df[f"gain"] == 3000)
]

# binary time - average
merged_df["time"] = np.mean(
    merged_df[["time_lh", "time_ll", "time_rh", "time_rl"]], axis=1
)
merged_df = merged_df.drop(columns=["time_lh", "time_ll", "time_rh", "time_rl"])

merged_df = merged_df.reset_index(drop=True)

# ENGINEERING DATA - all except xcal which will be added only to the calibration data later

gmt_eng = np.array(fdq_eng["ct_head/gmt"][()]).astype(str)
gmt_eng_parsed = []
for gmt in gmt_eng:
    gmt_eng_parsed.append(parse_date_string(gmt))

currents = ["hi", "lo"]
sides = ["a", "b"]

ical = {}
dihedral = {}

for side in sides:
    for current in currents:
        ical[f"{side}_{current}"] = fdq_eng[f"en_analog/grt/{side}_{current}_ical"][()]
        dihedral[f"{side}_{current}"] = fdq_eng[
            f"en_analog/grt/{side}_{current}_dihedral"
        ][()]


bol_volt = {}
bol_cmd_bias = {}

for channel, chan in channels.items():
    # bolometer voltages
    bol_volt[channel] = fdq_eng["en_analog/group1/bol_volt"][:, chan]
    # commanded bolometer bias?
    bol_cmd_bias[channel] = fdq_eng["en_stat/bol_cmd_bias"][:, chan]

stat_word_1 = fdq_eng["en_stat/stat_word_1"][()]
stat_word_12 = fdq_eng["en_stat/stat_word_12"][()]
stat_word_13 = fdq_eng["en_stat/stat_word_13"][()]
stat_word_16 = fdq_eng["en_stat/stat_word_16"][()]
stat_word_4 = fdq_eng["en_stat/stat_word_4"][()]
stat_word_5 = fdq_eng["en_stat/stat_word_5"][()]
stat_word_8 = fdq_eng["en_stat/stat_word_8"][()]
stat_word_9 = fdq_eng["en_stat/stat_word_9"][()]
int_ref_temp_a = fdq_eng["en_stat/int_ref_temp_a"][()]
int_ref_temp_b = fdq_eng["en_stat/int_ref_temp_b"][()]
ref_hrn_temp_a = fdq_eng["en_stat/ref_hrn_temp_a"][()]
ref_hrn_temp_b = fdq_eng["en_stat/ref_hrn_temp_b"][()]
ext_cal_temp_a = fdq_eng["en_stat/ext_cal_temp_a"][()]
ext_cal_temp_b = fdq_eng["en_stat/ext_cal_temp_b"][()]
lvdt_stat_a = fdq_eng["en_stat/lvdt_stat"][:, 0]
lvdt_stat_b = fdq_eng["en_stat/lvdt_stat"][:, 1]
power_a_status_a = fdq_eng["en_stat/power_a_status"][:, 0]
power_a_status_b = fdq_eng["en_stat/power_a_status"][:, 1]
power_b_status_a = fdq_eng["en_stat/power_b_status"][:, 0]
power_b_status_b = fdq_eng["en_stat/power_b_status"][:, 1]
hot_spot_cmd_a = fdq_eng["en_stat/hot_spot_cmd"][:, 0]
hot_spot_cmd_b = fdq_eng["en_stat/hot_spot_cmd"][:, 1]

# make engineeering data df
df_eng = pd.DataFrame(
    {
        "gmt": gmt_eng_parsed,
        "a_hi_ical": list(ical[f"a_hi"]),
        "a_lo_ical": list(ical[f"a_lo"]),
        "b_hi_ical": list(ical[f"b_hi"]),
        "b_lo_ical": list(ical[f"b_lo"]),
        "a_hi_dihedral": list(dihedral["a_hi"]),
        "a_lo_dihedral": list(dihedral["a_lo"]),
        "b_hi_dihedral": list(dihedral["b_hi"]),
        "b_lo_dihedral": list(dihedral["b_lo"]),
        "bol_volt_rh": list(bol_volt["rh"]),
        "bol_volt_rl": list(bol_volt["rl"]),
        "bol_volt_lh": list(bol_volt["lh"]),
        "bol_volt_ll": list(bol_volt["ll"]),
        "bol_cmd_bias_rh": list(bol_cmd_bias["rh"]),
        "bol_cmd_bias_rl": list(bol_cmd_bias["rl"]),
        "bol_cmd_bias_lh": list(bol_cmd_bias["lh"]),
        "bol_cmd_bias_ll": list(bol_cmd_bias["ll"]),
        "stat_word_1": list(stat_word_1),
        "stat_word_12": list(stat_word_12),
        "stat_word_13": list(stat_word_13),
        "stat_word_16": list(stat_word_16),
        "stat_word_4": list(stat_word_4),
        "stat_word_5": list(stat_word_5),
        "stat_word_8": list(stat_word_8),
        "stat_word_9": list(stat_word_9),
        "int_ref_temp_a": list(int_ref_temp_a),
        "int_ref_temp_b": list(int_ref_temp_b),
        "ref_hrn_temp_a": list(ref_hrn_temp_a),
        "ref_hrn_temp_b": list(ref_hrn_temp_b),
        "ext_cal_temp_a": list(ext_cal_temp_a),
        "ext_cal_temp_b": list(ext_cal_temp_b),
        "lvdt_stat_a": list(lvdt_stat_a),
        "lvdt_stat_b": list(lvdt_stat_b),
        "power_a_status_a": list(power_a_status_a),
        "power_a_status_b": list(power_a_status_b),
        "power_b_status_a": list(power_b_status_a),
        "power_b_status_b": list(power_b_status_b),
        "hot_spot_cmd_a": list(hot_spot_cmd_a),
        "hot_spot_cmd_b": list(hot_spot_cmd_b),
    }
)

print(f"Columns of df_eng: {df_eng.columns}")

df_eng["a_ical"] = df_eng.apply(get_temperature_hl, axis=1, args=("ical", "a"))
df_eng["b_ical"] = df_eng.apply(get_temperature_hl, axis=1, args=("ical", "b"))
df_eng = df_eng.drop(
    columns=[
        "a_hi_ical",
        "a_lo_ical",
        "b_hi_ical",
        "b_lo_ical",
    ]
)
df_eng["a_dihedral"] = df_eng.apply(get_temperature_hl, axis=1, args=("dihedral", "a"))
df_eng["b_dihedral"] = df_eng.apply(get_temperature_hl, axis=1, args=("dihedral", "b"))
df_eng = df_eng.drop(
    columns=[
        "a_hi_dihedral",
        "a_lo_dihedral",
        "b_hi_dihedral",
        "b_lo_dihedral",
    ]
)

# precompute timestamps for merged_df and df_eng
science_times = merged_df["gmt"].apply(lambda x: x.timestamp()).values
engineering_times = df_eng["gmt"].apply(lambda x: x.timestamp()).values

# initialise dict for all of the arrays that need to be interpolated
eng_variables = [
    "a_ical",
    "b_ical",
    "a_dihedral",
    "b_dihedral",
    "bol_volt_rh",
    "bol_volt_rl",
    "bol_volt_lh",
    "bol_volt_ll",
    "bol_cmd_bias_rh",
    "bol_cmd_bias_rl",
    "bol_cmd_bias_lh",
    "bol_cmd_bias_ll",
    "stat_word_1",
    "stat_word_12",
    "stat_word_13",
    "stat_word_16",
    "stat_word_4",
    "stat_word_5",
    "stat_word_8",
    "stat_word_9",
    "int_ref_temp_a",
    "int_ref_temp_b",
    "ref_hrn_temp_a",
    "ref_hrn_temp_b",
    "ext_cal_temp_a",
    "ext_cal_temp_b",
    "lvdt_stat_a",
    "lvdt_stat_b",
    "power_a_status_a",
    "power_a_status_b",
    "power_b_status_a",
    "power_b_status_b",
    "hot_spot_cmd_a",
    "hot_spot_cmd_b",
]

# initialise dict for the new interpolated arrays
new_variables = {}
for variable in eng_variables:
    new_variables[variable] = np.empty(len(science_times))

# iterate over science times and find engineering indices within 1 minute
for target_time, i in zip(science_times, range(len(science_times))):
    # find engineering times within 1 minute of the current science time - TODO: change the constraint?
    lower_bound = target_time - 60
    upper_bound = target_time + 60

    # get indices of engineering times within 1-minute bounds using searchsorted
    start_idx = np.searchsorted(engineering_times, lower_bound, side="left")
    end_idx = np.searchsorted(engineering_times, upper_bound, side="right")
    indices = range(start_idx, end_idx)

    # check if indices are available
    if len(indices) > 1:  # interpolate only if there are multiple points in the range
        for variable in eng_variables:
            f = interpolate.interp1d(
                engineering_times[indices],
                df_eng[variable][indices],
                fill_value="extrapolate",
            )
            new_variables[variable][i] = f(target_time)
    else:
        for variable in eng_variables:
            new_variables[variable][i] = np.nan

for variable in eng_variables:
    merged_df[variable] = new_variables[variable]

# REJECTS
# dihedral temperature must be lower than or equal to 5.5 K
merged_df = merged_df[
    (merged_df["a_dihedral"] <= 5.5) | (merged_df["b_dihedral"] <= 5.5)
]

# drop rows without ical data
merged_df = merged_df[(merged_df["a_ical"].notna()) & (merged_df["b_ical"].notna())]

# SPLIT CALIBRATION AND SKY DATA
calibration_df = merged_df[merged_df["xcal_pos"] == 1]
calibration_df = calibration_df.drop(columns=["xcal_pos"])
calibration_df = calibration_df.reset_index(drop=True)
sky_df = merged_df[merged_df["xcal_pos"] == 2]
sky_df = sky_df.drop(columns=["xcal_pos"])
sky_df = sky_df.reset_index(drop=True)

# clean up for upmode
sky_df = sky_df[(sky_df["upmode"] == 4)]
sky_df = sky_df.drop(columns=["upmode"])

# REJECTS
# sun angle must be greater than 91.2 degrees - removing if any of the channels have a sun angle less than 91.2
sky_df = sky_df[sky_df[f"sun_angle"] * fact > 91.2]
sky_df = sky_df.drop(columns=[f"sun_angle"])

# eart limb must be greater than 87 degrees - removing if any of the channels have an earth limb less than 87
sky_df = sky_df[sky_df[f"earth_limb"] * fact > 87]
sky_df = sky_df.drop(columns=[f"earth_limb"])

# moon angle must be greater than 22 degrees - removing if any of the channels have a moon angle less than 22
sky_df = sky_df[sky_df[f"moon_angle"] * fact > 22]
sky_df = sky_df.drop(columns=[f"moon_angle"])

# calibration dataset doesnt need any of these variables either
calibration_df = calibration_df.drop(
    columns=[
        f"sun_angle",
        f"earth_limb",
        f"moon_angle",
    ]
)


# reset index
sky_df = sky_df.reset_index(drop=True)

print(f"Sky dataset after cleaning science data: {sky_df.tail()}")
print(f"Calibration dataset after cleaning science data: {calibration_df.tail()}")

# ENGINEERING DATA - adding the missing xcal info
xcal_cone = {}
xcal_tip = {}

for side in sides:
    for current in currents:
        xcal_cone[f"{side}_{current}"] = fdq_eng[
            f"en_analog/grt/{side}_{current}_xcal_cone"
        ][()]
        # important to note that all xcal tip b side values are nonsense - for now the tip is actually ignored in get_temperature_hl()
        xcal_tip[f"{side}_{current}"] = fdq_eng[
            f"en_analog/grt/{side}_{current}_xcal_tip"
        ][()]

df_xcal = pd.DataFrame(
    {
        "gmt": gmt_eng_parsed,
        "a_hi_xcal_cone": list(xcal_cone["a_hi"]),
        "a_hi_xcal_tip": list(xcal_tip["a_hi"]),
        "a_lo_xcal_cone": list(xcal_cone["a_lo"]),
        "a_lo_xcal_tip": list(xcal_tip["a_lo"]),
        "b_hi_xcal_cone": list(xcal_cone["b_hi"]),
        "b_hi_xcal_tip": list(xcal_tip["b_hi"]),
        "b_lo_xcal_cone": list(xcal_cone["b_lo"]),
        "b_lo_xcal_tip": list(xcal_tip["b_lo"]),
    }
)

df_xcal["a_xcal"] = df_xcal.apply(get_temperature_hl, axis=1, args=("xcal", "a"))
df_xcal["b_xcal"] = df_xcal.apply(get_temperature_hl, axis=1, args=("xcal", "b"))
df_xcal = df_xcal.drop(
    columns=[
        "a_hi_xcal_cone",
        "a_hi_xcal_tip",
        "a_lo_xcal_cone",
        "a_lo_xcal_tip",
        "b_hi_xcal_cone",
        "b_hi_xcal_tip",
        "b_lo_xcal_cone",
        "b_lo_xcal_tip",
    ]
)

# precompute timestamps for calibration_df and df_xcal
science_times = calibration_df["gmt"].apply(lambda x: x.timestamp()).values
engineering_times = df_xcal["gmt"].apply(lambda x: x.timestamp()).values

a_xcal_new = np.empty(len(science_times))
b_xcal_new = np.empty(len(science_times))

# iterate over science times and find engineering indices within 1 minute
for target_time, i in zip(science_times, range(len(science_times))):
    # find engineering times within 1 minute of the current science time - TODO: change the constraint?
    lower_bound = target_time - 60
    upper_bound = target_time + 60

    # get indices of engineering times within 1-minute bounds using searchsorted
    start_idx = np.searchsorted(engineering_times, lower_bound, side="left")
    end_idx = np.searchsorted(engineering_times, upper_bound, side="right")
    indices = range(start_idx, end_idx)

    # check if indices are available
    if len(indices) > 1:
        # interpolate only if there are multiple points in the range
        f = interpolate.interp1d(
            engineering_times[indices],
            df_xcal["a_xcal"][indices],
            fill_value="extrapolate",
        )
        a_xcal_new[i] = f(target_time)
        f = interpolate.interp1d(
            engineering_times[indices],
            df_xcal["b_xcal"][indices],
            fill_value="extrapolate",
        )
        b_xcal_new[i] = f(target_time)
    else:
        a_xcal_new[i] = np.nan
        b_xcal_new[i] = np.nan

calibration_df["a_xcal"] = a_xcal_new
calibration_df["b_xcal"] = b_xcal_new

# converting gmt to string so we can save
gmt_str_sky = sky_df["gmt"].dt.strftime("%Y-%m-%d %H:%M:%S").to_numpy(dtype="S")
gmt_str_cal = calibration_df["gmt"].dt.strftime("%Y-%m-%d %H:%M:%S").to_numpy(dtype="S")

sky_df = sky_df.reset_index(drop=True)
calibration_df = calibration_df.reset_index(drop=True)

print(f"Dataframe after merging engineering data: {sky_df.tail()}")
print(f"Dataframe after merging engineering data: {calibration_df.tail()}")

print("Column names in sky_df:", sky_df.columns)
print("Column names in calibration_df:", calibration_df.columns)

zero_list = [0] * 512
# filling out ifg nans with zeros
for channel in channels:
    sky_df[f"ifg_{channel}"] = sky_df[f"ifg_{channel}"].apply(
        lambda x: zero_list if (isinstance(x, float) and np.isnan(x)) else x
    )
    calibration_df[f"ifg_{channel}"] = calibration_df[f"ifg_{channel}"].apply(
        lambda x: zero_list if (isinstance(x, float) and np.isnan(x)) else x
    )

sky_variables = eng_variables + [
    "mtm_length",
    "mtm_speed",
    "pix_gal",
    "pix_ecl",
    "pix_cel",
    "pix_terr",
    "scan",
    "adds_per_group",
    "time",
    "gain",
    "sweeps",
]

# saving to a h5 file
with tb.open_file("./../../data/sky_v2_m60sec.h5", mode="w") as h5file:
    group = h5file.create_group("/", "df_data", "Sky Data")

    h5file.create_array(group, "gmt", gmt_str_sky)
    h5file.create_array(group, "ifg_lh", np.stack(sky_df["ifg_lh"].values))
    h5file.create_array(group, "ifg_ll", np.stack(sky_df["ifg_ll"].values))
    h5file.create_array(group, "ifg_rh", np.stack(sky_df["ifg_rh"].values))
    h5file.create_array(group, "ifg_rl", np.stack(sky_df["ifg_rl"].values))

    # for loop to save all the other variables
    for variable in sky_variables:
        h5file.create_array(group, variable, sky_df[variable].values)

cal_variables = sky_variables + [
    "a_xcal",
    "b_xcal",
    "upmode",
]

# saving to a h5 file
with tb.open_file("./../../data/cal_v2_m60sec.h5", mode="w") as h5file:
    group = h5file.create_group("/", "df_data", "Calibration Data")

    h5file.create_array(group, "gmt", gmt_str_sky)
    h5file.create_array(group, "ifg_lh", np.stack(calibration_df["ifg_lh"].values))
    h5file.create_array(group, "ifg_ll", np.stack(calibration_df["ifg_ll"].values))
    h5file.create_array(group, "ifg_rh", np.stack(calibration_df["ifg_rh"].values))
    h5file.create_array(group, "ifg_rl", np.stack(calibration_df["ifg_rl"].values))

    # for loop to save all the other variables
    for variable in cal_variables:
        h5file.create_array(group, variable, calibration_df[variable].values)

end_time = time.time()
print(f"Time taken: {(end_time - start_time)/60} minutes")
