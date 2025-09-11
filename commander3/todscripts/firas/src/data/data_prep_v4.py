"""
After re-reading the Explanatory Supplement I think the engineering data is not completely raw and is already interpolated, i.e. the gmt values it has are not the "true" ones. So I am going back to a different way of pre-processing the data to compare and see which makes more sense.
"""

import time

import h5py
import healpy as hp
import numpy as np
import pandas as pd
import tables as tb

import globals as g
from utils.my_utils import (
    clean_variable,
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

eng_time = {}
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
solution = {}

NSIDE = 32
pix_terr = {}  # earth coordinates
pix_ecl = {}  # ecliptic coordinates
pix_gal = {}  # galactic coordinates
gal_lat = {}
gal_lon = {}
ecl_lat = {}
ecl_lon = {}
pix_cel = {}  # celestial coordinates, probably J1950
# Longitudes and latitudes are stored in radians*1e4
fact = 180.0 / np.pi / 1e4
scan = {}

for channel in channels:
    eng_time[channel] = fdq_sdf[f"fdq_sdf_{channel}/dq_data/eng_time"][()]
    print(f"eng_time: {eng_time[channel].shape}")
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
    sweeps[channel] = fdq_sdf[f"fdq_sdf_{channel}/sci_head/sc_head11"][()]
    gain[channel] = fdq_sdf[f"fdq_sdf_{channel}/sci_head/gain"][()]
    sun_angle[channel] = fdq_sdf[f"fdq_sdf_{channel}/attitude/sun_angle"][()]
    earth_limb[channel] = fdq_sdf[f"fdq_sdf_{channel}/attitude/earth_limb"][()]
    moon_angle[channel] = fdq_sdf[f"fdq_sdf_{channel}/attitude/moon_angle"][()]
    solution[channel] = fdq_sdf[f"fdq_sdf_{channel}/attitude/solution"][()]

    lon = fdq_sdf[f"fdq_sdf_{channel}/attitude/terr_longitude"][()] * fact
    lat = fdq_sdf[f"fdq_sdf_{channel}/attitude/terr_latitude"][()] * fact
    pix_terr[channel] = hp.ang2pix(NSIDE, lon, lat, lonlat=True).astype(float)

    lon = fdq_sdf[f"fdq_sdf_{channel}/attitude/ecliptic_longitude"][()] * fact
    lat = fdq_sdf[f"fdq_sdf_{channel}/attitude/ecliptic_latitude"][()] * fact
    pix_ecl[channel] = hp.ang2pix(NSIDE, lon, lat, lonlat=True).astype(float)
    ecl_lat[channel] = lat
    ecl_lon[channel] = lon
    print(
        f"ecliptic lat for channel {channel} max: {np.max(ecl_lat[channel])} and min: {np.min(ecl_lat[channel])}"
    )
    print(
        f"ecliptic lon for channel {channel} max: {np.max(ecl_lon[channel])} and min: {np.min(ecl_lon[channel])}"
    )
    print(f"getting up/down for channel {channel}")
    scan[channel] = scan_up_down(lat).astype(int)

    lon = fdq_sdf[f"fdq_sdf_{channel}/attitude/galactic_longitude"][()] * fact
    lat = fdq_sdf[f"fdq_sdf_{channel}/attitude/galactic_latitude"][()] * fact
    pix_gal[channel] = hp.ang2pix(NSIDE, lon, lat, lonlat=True).astype(float)
    gal_lat[channel] = lat
    gal_lon[channel] = lon

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
            "eng_time": list(eng_time[channel]),
            "ifg": list(ifg[channel]),
            "xcal_pos": list(xcal_pos[channel]),
            "mtm_length": list(mtm_length[channel]),
            "mtm_speed": list(mtm_speed[channel]),
            "fake": list(fake[channel]),
            "upmode": list(upmode[channel]),
            "adds_per_group": list(adds_per_group[channel]),
            "pix_gal": list(pix_gal[channel]),
            "gal_lat": list(gal_lat[channel]),
            "gal_lon": list(gal_lon[channel]),
            "ecl_lat": list(ecl_lat[channel]),
            "ecl_lon": list(ecl_lon[channel]),
            "pix_ecl": list(pix_ecl[channel]),
            "pix_cel": list(pix_cel[channel]),
            "pix_terr": list(pix_terr[channel]),
            "scan": list(scan[channel]),
            "sweeps": list(sweeps[channel]),
            "gain": list(gain[channel]),
            "sun_angle": list(sun_angle[channel]),
            "earth_limb": list(earth_limb[channel]),
            "moon_angle": list(moon_angle[channel]),
            "solution": list(solution[channel]),
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

# Check for duplicates in the merge keys
print(f"Duplicates in df['lh']['eng_time']: {df['lh']['eng_time'].duplicated().sum()}")
print(f"Duplicates in df['ll']['eng_time']: {df['ll']['eng_time'].duplicated().sum()}")
print(f"Duplicates in df['rh']['eng_time']: {df['rh']['eng_time'].duplicated().sum()}")

# USING THIS FOR MERGE_ASOF - need to keep the merge_asof on gmt because multiple science records correspond to one engineering record
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

merged_df = pd.merge_asof(
    merged_df,
    df["ll"],
    on="gmt",
    direction="nearest",
    suffixes=("_lh", "_ll"),
    tolerance=tolerance,
    by=["pix_gal", "pix_ecl", "pix_cel", "pix_terr"],
)
merged_df = pd.merge_asof(
    merged_df,
    df["rh"],
    on="gmt",
    direction="nearest",
    suffixes=("_ll", "_rh"),
    tolerance=tolerance,
    by=["pix_gal", "pix_ecl", "pix_cel", "pix_terr"],
)
merged_df = pd.merge_asof(
    merged_df,
    df["rl"],
    on="gmt",
    direction="nearest",
    suffixes=("_rh", "_rl"),
    tolerance=tolerance,
    by=["pix_gal", "pix_ecl", "pix_cel", "pix_terr"],
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

# set ids
merged_df["index"] = np.arange(len(merged_df))

# SCIENCE REJECTS

# making sure engineering record is the same for each channel
merged_df["eng_time"] = merged_df.apply(clean_variable, axis=1, args=("eng_time",))
merged_df = merged_df[merged_df["eng_time"] != np.nan]
merged_df = merged_df.drop(
    columns=["eng_time_lh", "eng_time_ll", "eng_time_rh", "eng_time_rl"]
)
# CLEANING XCAL_POS AND CONSTRAINING FOR ONLY 1 AND 2
# making sure xcal_pos is the same within each record
before = len(merged_df)
merged_df["xcal_pos"] = merged_df.apply(clean_variable, axis=1, args=("xcal_pos",))
merged_df = merged_df.drop(
    columns=["xcal_pos_lh", "xcal_pos_ll", "xcal_pos_rh", "xcal_pos_rl"]
)
merged_df = merged_df[merged_df["xcal_pos"] != np.nan]
print(f"Number of rows removed due to xcal_pos: {before - len(merged_df)}")
merged_df = merged_df[(merged_df["xcal_pos"] == 1) | (merged_df["xcal_pos"] == 2)]

# same thing but for mtm length and speed
before = len(merged_df)
merged_df["mtm_length"] = merged_df.apply(clean_variable, axis=1, args=("mtm_length",))
print(f"Number of rows removed due to mtm_length: {before - len(merged_df)}")
before = len(merged_df)
merged_df["mtm_speed"] = merged_df.apply(clean_variable, axis=1, args=("mtm_speed",))
print(f"Number of rows removed due to mtm_speed: {before - len(merged_df)}")
merged_df = merged_df.drop(
    columns=["mtm_length_lh", "mtm_length_ll", "mtm_length_rh", "mtm_length_rl"]
)
merged_df = merged_df.drop(
    columns=["mtm_speed_lh", "mtm_speed_ll", "mtm_speed_rh", "mtm_speed_rl"]
)
merged_df = merged_df[
    (merged_df["mtm_length"] == 0) | (merged_df["mtm_length"] == 1)
]  # only 0 and 1 for mtm_length
merged_df = merged_df[
    (merged_df["mtm_speed"] == 0) | (merged_df["mtm_speed"] == 1)
]  # only 0 and 1 for mtm_speed

# same thing for fake it flag
before = len(merged_df)
merged_df["fake"] = merged_df.apply(clean_variable, axis=1, args=("fake",))
print(f"Number of rows removed due to fake: {before - len(merged_df)}")
# fake-it mode on is 1, so i am guessing that 0 is when it is off
merged_df = merged_df[(merged_df["fake"] == 0)]
merged_df = merged_df.drop(columns=["fake", "fake_lh", "fake_ll", "fake_rh", "fake_rl"])

# upmode making all channels the same
merged_df["upmode"] = merged_df.apply(clean_variable, axis=1, args=("upmode",))
merged_df = merged_df.drop(columns=["upmode_lh", "upmode_ll", "upmode_rh", "upmode_rl"])

# adds per group
merged_df["adds_per_group"] = merged_df.apply(
    clean_variable, axis=1, args=("adds_per_group",)
)
merged_df = merged_df.drop(
    columns=[
        "adds_per_group_lh",
        "adds_per_group_ll",
        "adds_per_group_rh",
        "adds_per_group_rl",
    ]
)

print(f"Amount of lines before merging adds_per_group: {len(merged_df)}")
before = len(merged_df)
merged_df = merged_df[
    (merged_df["adds_per_group"] == 1)
    | (merged_df["adds_per_group"] == 2)
    | (merged_df["adds_per_group"] == 3)
    | (merged_df["adds_per_group"] == 8)
    | (merged_df["adds_per_group"] == 12)
]
print(f"Number of rows removed due to adds_per_group: {before - len(merged_df)}")
print(f"Amount of lines after merging adds_per_group: {len(merged_df)}")

# sweeps
merged_df["sweeps"] = merged_df.apply(clean_variable, axis=1, args=("sweeps",))
merged_df = merged_df.drop(columns=["sweeps_lh", "sweeps_ll", "sweeps_rh", "sweeps_rl"])
merged_df = merged_df[
    # (merged_df["sweeps"] == 1) |
    (merged_df["sweeps"] == 4)
    | (merged_df["sweeps"] == 16)
]

# attitude solution must be available (i.e. non-zero)
merged_df["solution"] = merged_df.apply(clean_variable, axis=1, args=("solution",))
merged_df = merged_df[merged_df["solution"] != 0]
merged_df = merged_df.drop(
    columns=["solution_lh", "solution_ll", "solution_rh", "solution_rl", "solution"]
)

# gain
# merged_df["gain"] = merged_df.apply(clean_variable, axis=1, args=("gain",))
# merged_df = merged_df.drop(columns=["gain_lh", "gain_ll", "gain_rh", "gain_rl"])
# merged_df["gain"] = merged_df.apply(convert_gain, axis=1)
# merged_df = merged_df[
#     (merged_df["gain"] == 1)
#     | (merged_df["gain"] == 3)
#     | (merged_df["gain"] == 10)
#     | (merged_df["gain"] == 30)
#     | (merged_df["gain"] == 100)
#     | (merged_df["gain"] == 300)
#     | (merged_df["gain"] == 1000)
#     | (merged_df["gain"] == 3000)
# ]


# the following variables i don't think necessarily need to be the same in all channels for the record to be valid
for channel in channels:
    # adds per group
    # merged_df = merged_df[
    #     (merged_df[f"adds_per_group_{channel}"] == 1)
    #     | (merged_df[f"adds_per_group_{channel}"] == 2)
    #     | (merged_df[f"adds_per_group_{channel}"] == 3)
    #     | (merged_df[f"adds_per_group_{channel}"] == 8)
    #     | (merged_df[f"adds_per_group_{channel}"] == 12)
    # ]
    # # check for nans in adds_per_group
    # if np.isnan(merged_df[f"adds_per_group_{channel}"]).any():
    #     print(f"There are nans in adds_per_group (checkpoint 3, channel {channel})")
    #     # how many nans are there?
    #     print(np.sum(np.isnan(merged_df[f"adds_per_group_{channel}"])))
    # else:
    #     print(f"No nans in adds_per_group (checkpoint 3, channel {channel})")
    # sweeps
    # merged_df = merged_df[
    #     # (merged_df[f"sweeps_{channel}"] == 1) |
    #     (merged_df[f"sweeps_{channel}"] == 4)
    #     | (merged_df[f"sweeps_{channel}"] == 16)
    # ]

    # gain
    merged_df[f"gain_{channel}"] = merged_df.apply(
        convert_gain, axis=1, args=(channel,)
    )
    merged_df = merged_df[
        (merged_df[f"gain_{channel}"] == 1)
        | (merged_df[f"gain_{channel}"] == 3)
        | (merged_df[f"gain_{channel}"] == 10)
        | (merged_df[f"gain_{channel}"] == 30)
        | (merged_df[f"gain_{channel}"] == 100)
        | (merged_df[f"gain_{channel}"] == 300)
        | (merged_df[f"gain_{channel}"] == 1000)
        | (merged_df[f"gain_{channel}"] == 3000)
    ]

# scan should be the same for all channels
merged_df["scan"] = merged_df.apply(clean_variable, axis=1, args=("scan",))
merged_df = merged_df.drop(columns=["scan_lh", "scan_ll", "scan_rh", "scan_rl"])
merged_df = merged_df[(merged_df["scan"] == -1) | (merged_df["scan"] == 1)]

print(f"Amount of lines before merging scan: {len(merged_df)}")

# galactic latitude and longitude
merged_df["gal_lat"] = merged_df.apply(clean_variable, axis=1, args=("gal_lat",))
merged_df["gal_lon"] = merged_df.apply(clean_variable, axis=1, args=("gal_lon",))
merged_df = merged_df[(merged_df["gal_lat"] != np.nan)]
merged_df = merged_df[(merged_df["gal_lon"] != np.nan)]
merged_df = merged_df.drop(
    columns=[
        "gal_lat_lh",
        "gal_lat_ll",
        "gal_lat_rh",
        "gal_lat_rl",
        "gal_lon_lh",
        "gal_lon_ll",
        "gal_lon_rh",
        "gal_lon_rl",
    ]
)

print(f"Amount of lines before merging ecliptic lat/lon: {len(merged_df)}")

lines_before = len(merged_df)
merged_df["ecl_lat"] = merged_df.apply(clean_variable, axis=1, args=("ecl_lat",))
merged_df["ecl_lon"] = merged_df.apply(clean_variable, axis=1, args=("ecl_lon",))
merged_df = merged_df[(merged_df["ecl_lat"] != np.nan)]
merged_df = merged_df[(merged_df["ecl_lon"] != np.nan)]
merged_df = merged_df.drop(
    columns=[
        "ecl_lat_lh",
        "ecl_lat_ll",
        "ecl_lat_rh",
        "ecl_lat_rl",
        "ecl_lon_lh",
        "ecl_lon_ll",
        "ecl_lon_rh",
        "ecl_lon_rl",
    ]
)

merged_df = merged_df.reset_index(drop=True)
lines_after = len(merged_df)
print(
    f"Merging ecliptic lat/lon reduced the number of lines by {lines_before - lines_after}. Amount now: {lines_after}"
)

# ENGINEERING DATA - all except xcal which will be added only to the calibration data later

currents = ["hi", "lo"]
sides = ["a", "b"]

ical = {}
dihedral = {}
refhorn = {}
skyhorn = {}
collimator = {}
bol_assem = {}

print(f"grts: {fdq_eng['en_analog/grt'].keys()}")

for side in sides:
    # there are no high current readings for the bolometers
    bol_assem[f"{side}_rh"] = fdq_eng[f"en_analog/grt/{side}_lo_bol_assem"][:, 0]
    bol_assem[f"{side}_rl"] = fdq_eng[f"en_analog/grt/{side}_lo_bol_assem"][:, 1]
    bol_assem[f"{side}_lh"] = fdq_eng[f"en_analog/grt/{side}_lo_bol_assem"][:, 2]
    bol_assem[f"{side}_ll"] = fdq_eng[f"en_analog/grt/{side}_lo_bol_assem"][:, 3]
    for current in currents:
        ical[f"{side}_{current}"] = fdq_eng[f"en_analog/grt/{side}_{current}_ical"][()]
        dihedral[f"{side}_{current}"] = fdq_eng[
            f"en_analog/grt/{side}_{current}_dihedral"
        ][()]
        refhorn[f"{side}_{current}"] = fdq_eng[
            f"en_analog/grt/{side}_{current}_refhorn"
        ][()]
        skyhorn[f"{side}_{current}"] = fdq_eng[
            f"en_analog/grt/{side}_{current}_skyhorn"
        ][()]
        collimator[f"{side}_{current}"] = fdq_eng[
            f"en_analog/grt/{side}_{current}_mirror"
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
sky_hrn_temp_a = fdq_eng["en_stat/sky_hrn_temp_a"][()]
sky_hrn_temp_b = fdq_eng["en_stat/sky_hrn_temp_b"][()]
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
t_eng = fdq_eng["ct_head/time"][()]
dwell_stat_a = fdq_eng["en_stat/dwell_stat"][:, 0]
dwell_stat_b = fdq_eng["en_stat/dwell_stat"][:, 1]
engstat_spares_1 = fdq_eng["en_stat/engstat_spares"][:, 0]
engstat_spares_2 = fdq_eng["en_stat/engstat_spares"][:, 1]
engstat_spares_3 = fdq_eng["en_stat/engstat_spares"][:, 2]
engstat_spares_4 = fdq_eng["en_stat/engstat_spares"][:, 3]
engstat_spares_5 = fdq_eng["en_stat/engstat_spares"][:, 4]
engstat_spares_6 = fdq_eng["en_stat/engstat_spares"][:, 5]
engstat_spares_7 = fdq_eng["en_stat/engstat_spares"][:, 6]
engstat_spares_8 = fdq_eng["en_stat/engstat_spares"][:, 7]
engstat_spares_9 = fdq_eng["en_stat/engstat_spares"][:, 8]
engstat_spares_10 = fdq_eng["en_stat/engstat_spares"][:, 9]
engstat_spares2_1 = fdq_eng["en_stat/engstat_spares2"][:, 0]
engstat_spares2_2 = fdq_eng["en_stat/engstat_spares2"][:, 1]
engstat_spares2_3 = fdq_eng["en_stat/engstat_spares2"][:, 2]
engstat_spares2_4 = fdq_eng["en_stat/engstat_spares2"][:, 3]
engstat_spares2_5 = fdq_eng["en_stat/engstat_spares2"][:, 4]
micro_stat_bus_1 = fdq_eng["en_stat/micro_stat_bus"][:, 0]
micro_stat_bus_2 = fdq_eng["en_stat/micro_stat_bus"][:, 1]
micro_stat_bus_3 = fdq_eng["en_stat/micro_stat_bus"][:, 2]
micro_stat_bus_4 = fdq_eng["en_stat/micro_stat_bus"][:, 3]
grt_addr_a = fdq_eng["en_stat/grt_addr"][:, 0]
grt_addr_b = fdq_eng["en_stat/grt_addr"][:, 1]


# make engineeering data df
df_eng = pd.DataFrame(
    {
        "a_hi_ical": list(ical[f"a_hi"]),
        "a_lo_ical": list(ical[f"a_lo"]),
        "b_hi_ical": list(ical[f"b_hi"]),
        "b_lo_ical": list(ical[f"b_lo"]),
        "a_hi_dihedral": list(dihedral["a_hi"]),
        "a_lo_dihedral": list(dihedral["a_lo"]),
        "b_hi_dihedral": list(dihedral["b_hi"]),
        "b_lo_dihedral": list(dihedral["b_lo"]),
        "a_hi_refhorn": list(refhorn["a_hi"]),
        "a_lo_refhorn": list(refhorn["a_lo"]),
        "b_hi_refhorn": list(refhorn["b_hi"]),
        "b_lo_refhorn": list(refhorn["b_lo"]),
        "a_hi_skyhorn": list(skyhorn["a_hi"]),
        "a_lo_skyhorn": list(skyhorn["a_lo"]),
        "b_hi_skyhorn": list(skyhorn["b_hi"]),
        "b_lo_skyhorn": list(skyhorn["b_lo"]),
        "a_hi_collimator": list(collimator["a_hi"]),
        "a_lo_collimator": list(collimator["a_lo"]),
        "b_hi_collimator": list(collimator["b_hi"]),
        "b_lo_collimator": list(collimator["b_lo"]),
        "a_bol_assem_rh": list(bol_assem["a_rh"]),
        "b_bol_assem_rh": list(bol_assem["b_rh"]),
        "a_bol_assem_rl": list(bol_assem["a_rl"]),
        "b_bol_assem_rl": list(bol_assem["b_rl"]),
        "a_bol_assem_lh": list(bol_assem["a_lh"]),
        "b_bol_assem_lh": list(bol_assem["b_lh"]),
        "a_bol_assem_ll": list(bol_assem["a_ll"]),
        "b_bol_assem_ll": list(bol_assem["b_ll"]),
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
        "time": list(t_eng),
        "dwell_stat_a": list(dwell_stat_a),
        "dwell_stat_b": list(dwell_stat_b),
        "engstat_spares_1": list(engstat_spares_1),
        "engstat_spares_2": list(engstat_spares_2),
        "engstat_spares_3": list(engstat_spares_3),
        "engstat_spares_4": list(engstat_spares_4),
        "engstat_spares_5": list(engstat_spares_5),
        "engstat_spares_6": list(engstat_spares_6),
        "engstat_spares_7": list(engstat_spares_7),
        "engstat_spares_8": list(engstat_spares_8),
        "engstat_spares_9": list(engstat_spares_9),
        "engstat_spares_10": list(engstat_spares_10),
        "engstat_spares2_1": list(engstat_spares2_1),
        "engstat_spares2_2": list(engstat_spares2_2),
        "engstat_spares2_3": list(engstat_spares2_3),
        "engstat_spares2_4": list(engstat_spares2_4),
        "engstat_spares2_5": list(engstat_spares2_5),
        "micro_stat_bus_1": list(micro_stat_bus_1),
        "micro_stat_bus_2": list(micro_stat_bus_2),
        "micro_stat_bus_3": list(micro_stat_bus_3),
        "micro_stat_bus_4": list(micro_stat_bus_4),
        "grt_addr_a": list(grt_addr_a),
        "grt_addr_b": list(grt_addr_b),
        "sky_hrn_temp_a": list(sky_hrn_temp_a),
        "sky_hrn_temp_b": list(sky_hrn_temp_b),
    }
)

print(f"Columns of df_eng: {df_eng.columns}")

df_eng["a_ical"] = df_eng.apply(get_temperature_hl, axis=1, args=("ical", "a")).astype(
    np.float64
)
df_eng["b_ical"] = df_eng.apply(get_temperature_hl, axis=1, args=("ical", "b")).astype(
    np.float64
)
df_eng = df_eng.drop(
    columns=[
        "a_hi_ical",
        "a_lo_ical",
        "b_hi_ical",
        "b_lo_ical",
    ]
)
df_eng["a_dihedral"] = df_eng.apply(
    get_temperature_hl, axis=1, args=("dihedral", "a")
).astype(np.float64)
df_eng["b_dihedral"] = df_eng.apply(
    get_temperature_hl, axis=1, args=("dihedral", "b")
).astype(np.float64)
df_eng = df_eng.drop(
    columns=[
        "a_hi_dihedral",
        "a_lo_dihedral",
        "b_hi_dihedral",
        "b_lo_dihedral",
    ]
)

df_eng["a_refhorn"] = df_eng.apply(
    get_temperature_hl, axis=1, args=("refhorn", "a", "reference horn")
).astype(np.float64)
df_eng["b_refhorn"] = df_eng.apply(
    get_temperature_hl, axis=1, args=("refhorn", "b", "reference horn")
).astype(np.float64)
df_eng = df_eng.drop(
    columns=[
        "a_hi_refhorn",
        "a_lo_refhorn",
        "b_hi_refhorn",
        "b_lo_refhorn",
    ]
)

df_eng["a_skyhorn"] = df_eng.apply(
    get_temperature_hl, axis=1, args=("skyhorn", "a")
).astype(np.float64)
df_eng["b_skyhorn"] = df_eng.apply(
    get_temperature_hl, axis=1, args=("skyhorn", "b")
).astype(np.float64)
df_eng = df_eng.drop(
    columns=[
        "a_hi_skyhorn",
        "a_lo_skyhorn",
        "b_hi_skyhorn",
        "b_lo_skyhorn",
    ]
)

df_eng["a_collimator"] = df_eng.apply(
    get_temperature_hl, axis=1, args=("collimator", "a")
).astype(np.float64)
df_eng["b_collimator"] = df_eng.apply(
    get_temperature_hl, axis=1, args=("collimator", "b")
).astype(np.float64)
df_eng = df_eng.drop(
    columns=[
        "a_hi_collimator",
        "a_lo_collimator",
        "b_hi_collimator",
        "b_lo_collimator",
    ]
)

# merge the two dataframes
merged_df = merged_df.merge(
    df_eng,
    left_on="eng_time",
    right_on="time",
    how="inner",
)

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
for channel in channels:
    # sun angle must be greater than 91.2 degrees - removing if any of the channels have a sun angle less than 91.2
    sky_df = sky_df[sky_df[f"sun_angle_{channel}"] * fact > 91.2]
    sky_df = sky_df.drop(columns=[f"sun_angle_{channel}"])

    # eart limb must be greater than 87 degrees - removing if any of the channels have an earth limb less than 87
    sky_df = sky_df[sky_df[f"earth_limb_{channel}"] * fact > 87]
    sky_df = sky_df.drop(columns=[f"earth_limb_{channel}"])

    # moon angle must be greater than 22 degrees - removing if any of the channels have a moon angle less than 22
    sky_df = sky_df[sky_df[f"moon_angle_{channel}"] * fact > 22]
    sky_df = sky_df.drop(columns=[f"moon_angle_{channel}"])

    # calibration dataset doesnt need any of these variables either
    calibration_df = calibration_df.drop(
        columns=[
            f"sun_angle_{channel}",
            f"earth_limb_{channel}",
            f"moon_angle_{channel}",
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
        "time": list(t_eng),
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

df_xcal["a_xcal"] = df_xcal.apply(
    get_temperature_hl, axis=1, args=("xcal", "a")
).astype(np.float64)
df_xcal["b_xcal"] = df_xcal.apply(
    get_temperature_hl, axis=1, args=("xcal", "b")
).astype(np.float64)
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

calibration_df = calibration_df.merge(
    df_xcal,
    left_on="eng_time",
    right_on="time",
    how="inner",
)

# converting gmt to string so we can save
gmt_str_sky = sky_df["gmt"].dt.strftime("%Y-%m-%d %H:%M:%S").to_numpy(dtype="S")
gmt_str_cal = calibration_df["gmt"].dt.strftime("%Y-%m-%d %H:%M:%S").to_numpy(dtype="S")

sky_df = sky_df.reset_index(drop=True)
calibration_df = calibration_df.reset_index(drop=True)

# set up ids for each dataframe
sky_df["id"] = sky_df.index + 1
calibration_df["id"] = calibration_df.index + 1

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

sky_variables = [
    "id",
    "a_ical",
    "b_ical",
    "a_dihedral",
    "b_dihedral",
    "a_refhorn",
    "b_refhorn",
    "a_skyhorn",
    "b_skyhorn",
    "a_collimator",
    "b_collimator",
    "a_bol_assem_rh",
    "b_bol_assem_rh",
    "a_bol_assem_rl",
    "b_bol_assem_rl",
    "a_bol_assem_lh",
    "b_bol_assem_lh",
    "a_bol_assem_ll",
    "b_bol_assem_ll",
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
    "mtm_length",
    "mtm_speed",
    "pix_gal",
    "gal_lat",
    "gal_lon",
    "ecl_lat",
    "ecl_lon",
    "pix_ecl",
    "pix_cel",
    "pix_terr",
    "scan",
    # "adds_per_group_lh",
    # "adds_per_group_ll",
    # "adds_per_group_rh",
    # "adds_per_group_rl",
    "adds_per_group",
    "gain_rh",
    "gain_rl",
    "gain_lh",
    "gain_ll",
    # "gain",
    # "sweeps_lh",
    # "sweeps_ll",
    # "sweeps_rh",
    # "sweeps_rl",
    "sweeps",
    "dwell_stat_a",
    "dwell_stat_b",
    "engstat_spares_1",
    "engstat_spares_2",
    "engstat_spares_3",
    "engstat_spares_4",
    "engstat_spares_5",
    "engstat_spares_6",
    "engstat_spares_7",
    "engstat_spares_8",
    "engstat_spares_9",
    "engstat_spares_10",
    "engstat_spares2_1",
    "engstat_spares2_2",
    "engstat_spares2_3",
    "engstat_spares2_4",
    "engstat_spares2_5",
    "micro_stat_bus_1",
    "micro_stat_bus_2",
    "micro_stat_bus_3",
    "micro_stat_bus_4",
    "grt_addr_a",
    "grt_addr_b",
    "sky_hrn_temp_a",
    "sky_hrn_temp_b",
    "index",
]

# saving to a h5 file
with tb.open_file(g.PREPROCESSED_DATA_PATH_SKY, mode="w") as h5file:
    group = h5file.create_group("/", "df_data", "Sky Data")

    h5file.create_array(group, "gmt", gmt_str_sky)
    h5file.create_array(group, "ifg_lh", np.stack(sky_df["ifg_lh"].values))
    h5file.create_array(group, "ifg_ll", np.stack(sky_df["ifg_ll"].values))
    h5file.create_array(group, "ifg_rh", np.stack(sky_df["ifg_rh"].values))
    h5file.create_array(group, "ifg_rl", np.stack(sky_df["ifg_rl"].values))

    # for loop to save all the other variables
    for variable in sky_variables:
        print(f"Saving {variable}")
        h5file.create_array(group, variable, sky_df[variable].values)

cal_variables = sky_variables + [
    "a_xcal",
    "b_xcal",
    "upmode",
]

print(f"Saving calibration data")

# saving to a h5 file
with tb.open_file(g.PREPROCESSED_DATA_PATH_CAL, mode="w") as h5file:
    group = h5file.create_group("/", "df_data", "Calibration Data")

    h5file.create_array(group, "gmt", gmt_str_cal)
    h5file.create_array(group, "ifg_lh", np.stack(calibration_df["ifg_lh"].values))
    h5file.create_array(group, "ifg_ll", np.stack(calibration_df["ifg_ll"].values))
    h5file.create_array(group, "ifg_rh", np.stack(calibration_df["ifg_rh"].values))
    h5file.create_array(group, "ifg_rl", np.stack(calibration_df["ifg_rl"].values))

    # for loop to save all the other variables
    for variable in cal_variables:
        print(f"Saving {variable}")
        h5file.create_array(group, variable, calibration_df[variable].values)

end_time = time.time()
print(f"Time taken: {(end_time - start_time)/60} minutes")
