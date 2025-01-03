"""
Data processing script has to be re-written again, not matching eng_time/time now,
and instead interpolating over gmt to get the engineering data corresponding to the science data.
"""

import time

import h5py
import healpy as hp
import numpy as np
import pandas as pd
import tables as tb
from scipy import interpolate
from utils.my_utils import (
    clean_variable,
    convert_gain,
    get_temperature_hl,
    parse_date_string,
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

# the following variables i don't think necessarily need to be the same in all channels for the record to be valid
for channel in channels:
    # adds per group
    merged_df = merged_df[
        (merged_df[f"adds_per_group_{channel}"] == 1)
        | (merged_df[f"adds_per_group_{channel}"] == 2)
        | (merged_df[f"adds_per_group_{channel}"] == 3)
        | (merged_df[f"adds_per_group_{channel}"] == 8)
        | (merged_df[f"adds_per_group_{channel}"] == 12)
    ]
    # check for nans in adds_per_group
    if np.isnan(merged_df[f"adds_per_group_{channel}"]).any():
        print(f"There are nans in adds_per_group (checkpoint 3, channel {channel})")
        # how many nans are there?
        print(np.sum(np.isnan(merged_df[f"adds_per_group_{channel}"])))
    else:
        print(f"No nans in adds_per_group (checkpoint 3, channel {channel})")
    # sweeps - TODO: CHECK IF 1 IS A VALID VALUE
    merged_df = merged_df[
        (merged_df[f"sweeps_{channel}"] == 1)
        | (merged_df[f"sweeps_{channel}"] == 4)
        | (merged_df[f"sweeps_{channel}"] == 16)
    ]

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
    }
)

print(f"Columns of df_eng: {df_eng.columns}")

df_eng["ical"] = df_eng.apply(get_temperature_hl, axis=1, args=("ical",))
df_eng = df_eng.drop(
    columns=[
        "a_hi_ical",
        "a_lo_ical",
        "b_hi_ical",
        "b_lo_ical",
    ]
)
df_eng["dihedral"] = df_eng.apply(get_temperature_hl, axis=1, args=("dihedral",))
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

# initialize list for interpolated values
ical_new = []
dihedral_new = []
# also need to interpolate bolometer voltages and commanded biases
bol_volt_rh_new = []
bol_volt_rl_new = []
bol_volt_lh_new = []
bol_volt_ll_new = []
bol_cmd_bias_rh_new = []
bol_cmd_bias_rl_new = []
bol_cmd_bias_lh_new = []
bol_cmd_bias_ll_new = []

# iterate over science times and find engineering indices within 1 minute
for target_time in science_times:
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
            df_eng["ical"][indices],
            fill_value="extrapolate",
        )
        ical_new.append(f(target_time))
        # xcal_pos_series = merged_df[merged_df["gmt"] == target_time]["xcal_pos"]
        # if not xcal_pos_series.empty and xcal_pos_series.iloc[0] == 1: # FIX TO NOT GET XCAL TEMP FOR ALL
        f = interpolate.interp1d(
            engineering_times[indices],
            df_eng["dihedral"][indices],
            fill_value="extrapolate",
        )
        dihedral_new.append(f(target_time))
        f = interpolate.interp1d(
            engineering_times[indices],
            df_eng["bol_volt_rh"][indices],
            fill_value="extrapolate",
        )
        bol_volt_rh_new.append(f(target_time))
        f = interpolate.interp1d(
            engineering_times[indices],
            df_eng["bol_volt_rl"][indices],
            fill_value="extrapolate",
        )
        bol_volt_rl_new.append(f(target_time))
        f = interpolate.interp1d(
            engineering_times[indices],
            df_eng["bol_volt_lh"][indices],
            fill_value="extrapolate",
        )
        bol_volt_lh_new.append(f(target_time))
        f = interpolate.interp1d(
            engineering_times[indices],
            df_eng["bol_volt_ll"][indices],
            fill_value="extrapolate",
        )
        bol_volt_ll_new.append(f(target_time))
        f = interpolate.interp1d(
            engineering_times[indices],
            df_eng["bol_cmd_bias_rh"][indices],
            fill_value="extrapolate",
        )
        bol_cmd_bias_rh_new.append(f(target_time))
        f = interpolate.interp1d(
            engineering_times[indices],
            df_eng["bol_cmd_bias_rl"][indices],
            fill_value="extrapolate",
        )
        bol_cmd_bias_rl_new.append(f(target_time))
        f = interpolate.interp1d(
            engineering_times[indices],
            df_eng["bol_cmd_bias_lh"][indices],
            fill_value="extrapolate",
        )
        bol_cmd_bias_lh_new.append(f(target_time))
        f = interpolate.interp1d(
            engineering_times[indices],
            df_eng["bol_cmd_bias_ll"][indices],
            fill_value="extrapolate",
        )
        bol_cmd_bias_ll_new.append(f(target_time))
        # elif not xcal_pos_series.empty and xcal_pos_series.iloc[0] == 2:
        #     xcal_new.append(np.nan)
    else:
        ical_new.append(np.nan)
        dihedral_new.append(np.nan)
        bol_volt_rh_new.append(np.nan)
        bol_volt_rl_new.append(np.nan)
        bol_volt_lh_new.append(np.nan)
        bol_volt_ll_new.append(np.nan)
        bol_cmd_bias_rh_new.append(np.nan)
        bol_cmd_bias_rl_new.append(np.nan)
        bol_cmd_bias_lh_new.append(np.nan)
        bol_cmd_bias_ll_new.append(np.nan)


merged_df["ical"] = ical_new
merged_df["dihedral"] = dihedral_new
merged_df["bol_volt_rh"] = bol_volt_rh_new
merged_df["bol_volt_rl"] = bol_volt_rl_new
merged_df["bol_volt_lh"] = bol_volt_lh_new
merged_df["bol_volt_ll"] = bol_volt_ll_new
merged_df["bol_cmd_bias_rh"] = bol_cmd_bias_rh_new
merged_df["bol_cmd_bias_rl"] = bol_cmd_bias_rl_new
merged_df["bol_cmd_bias_lh"] = bol_cmd_bias_lh_new
merged_df["bol_cmd_bias_ll"] = bol_cmd_bias_ll_new

# REJECTS
# dihedral temperature must be lower than or equal to 5.5 K
merged_df = merged_df[merged_df["dihedral"] <= 5.5]

# drop rows without ical data
merged_df = merged_df[merged_df["ical"].notna()]

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

df_xcal["xcal"] = df_xcal.apply(get_temperature_hl, axis=1, args=("xcal",))
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

xcal_new = []

# iterate over science times and find engineering indices within 1 minute
for target_time in science_times:
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
            df_xcal["xcal"][indices],
            fill_value="extrapolate",
        )
        xcal_new.append(f(target_time))
    else:
        xcal_new.append(np.nan)

calibration_df["xcal"] = xcal_new

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

# saving to a h5 file
with tb.open_file("./../../data/sky_v1.h5", mode="w") as h5file:
    group = h5file.create_group("/", "df_data", "Sky Data")

    h5file.create_array(group, "gmt", gmt_str_sky)
    h5file.create_array(group, "ifg_lh", np.stack(sky_df["ifg_lh"].values))
    h5file.create_array(group, "ifg_ll", np.stack(sky_df["ifg_ll"].values))
    h5file.create_array(group, "ifg_rh", np.stack(sky_df["ifg_rh"].values))
    h5file.create_array(group, "ifg_rl", np.stack(sky_df["ifg_rl"].values))
    h5file.create_array(group, "ical", sky_df["ical"].tolist())
    h5file.create_array(group, "dihedral", sky_df["dihedral"].tolist())
    h5file.create_array(group, "mtm_length", sky_df["mtm_length"].values)
    h5file.create_array(group, "mtm_speed", sky_df["mtm_speed"].values)
    h5file.create_array(group, "pix_gal", sky_df["pix_gal"].values)
    h5file.create_array(group, "pix_ecl", sky_df["pix_ecl"].values)
    h5file.create_array(group, "pix_cel", sky_df["pix_cel"].values)
    h5file.create_array(group, "pix_terr", sky_df["pix_terr"].values)
    h5file.create_array(group, "adds_per_group_lh", sky_df["adds_per_group_lh"].values)
    h5file.create_array(group, "adds_per_group_ll", sky_df["adds_per_group_ll"].values)
    h5file.create_array(group, "adds_per_group_rh", sky_df["adds_per_group_rh"].values)
    h5file.create_array(group, "adds_per_group_rl", sky_df["adds_per_group_rl"].values)
    h5file.create_array(group, "bol_volt_rh", sky_df["bol_volt_rh"].tolist())
    h5file.create_array(group, "bol_volt_rl", sky_df["bol_volt_rl"].tolist())
    h5file.create_array(group, "bol_volt_lh", sky_df["bol_volt_lh"].tolist())
    h5file.create_array(group, "bol_volt_ll", sky_df["bol_volt_ll"].tolist())
    h5file.create_array(group, "bol_cmd_bias_rh", sky_df["bol_cmd_bias_rh"].tolist())
    h5file.create_array(group, "bol_cmd_bias_rl", sky_df["bol_cmd_bias_rl"].tolist())
    h5file.create_array(group, "bol_cmd_bias_lh", sky_df["bol_cmd_bias_lh"].tolist())
    h5file.create_array(group, "bol_cmd_bias_ll", sky_df["bol_cmd_bias_ll"].tolist())
    h5file.create_array(group, "time", sky_df["time"].values)
    h5file.create_array(group, "gain_rh", sky_df["gain_rh"].tolist())
    h5file.create_array(group, "gain_rl", sky_df["gain_rl"].tolist())
    h5file.create_array(group, "gain_lh", sky_df["gain_lh"].tolist())
    h5file.create_array(group, "gain_ll", sky_df["gain_ll"].tolist())
    h5file.create_array(group, "sweeps_lh", sky_df["sweeps_lh"].values)
    h5file.create_array(group, "sweeps_ll", sky_df["sweeps_ll"].values)
    h5file.create_array(group, "sweeps_rh", sky_df["sweeps_rh"].values)
    h5file.create_array(group, "sweeps_rl", sky_df["sweeps_rl"].values)

# saving to a h5 file
with tb.open_file("./../../data/cal_v1.h5", mode="w") as h5file:
    group = h5file.create_group("/", "df_data", "Calibration Data")

    h5file.create_array(group, "gmt", gmt_str_sky)
    h5file.create_array(group, "ifg_lh", np.stack(calibration_df["ifg_lh"].values))
    h5file.create_array(group, "ifg_ll", np.stack(calibration_df["ifg_ll"].values))
    h5file.create_array(group, "ifg_rh", np.stack(calibration_df["ifg_rh"].values))
    h5file.create_array(group, "ifg_rl", np.stack(calibration_df["ifg_rl"].values))
    h5file.create_array(group, "ical", calibration_df["ical"].tolist())
    h5file.create_array(group, "dihedral", calibration_df["dihedral"].tolist())
    h5file.create_array(group, "mtm_length", calibration_df["mtm_length"].values)
    h5file.create_array(group, "mtm_speed", calibration_df["mtm_speed"].values)
    h5file.create_array(group, "pix_gal", calibration_df["pix_gal"].values)
    h5file.create_array(group, "pix_ecl", calibration_df["pix_ecl"].values)
    h5file.create_array(group, "pix_cel", calibration_df["pix_cel"].values)
    h5file.create_array(group, "pix_terr", calibration_df["pix_terr"].values)
    h5file.create_array(group, "upmode", calibration_df["upmode"].values)
    h5file.create_array(
        group, "adds_per_group_lh", calibration_df["adds_per_group_lh"].values
    )
    h5file.create_array(
        group, "adds_per_group_ll", calibration_df["adds_per_group_ll"].values
    )
    h5file.create_array(
        group, "adds_per_group_rh", calibration_df["adds_per_group_rh"].values
    )
    h5file.create_array(
        group, "adds_per_group_rl", calibration_df["adds_per_group_rl"].values
    )
    h5file.create_array(group, "bol_volt_rh", calibration_df["bol_volt_rh"].tolist())
    h5file.create_array(group, "bol_volt_rl", calibration_df["bol_volt_rl"].tolist())
    h5file.create_array(group, "bol_volt_lh", calibration_df["bol_volt_lh"].tolist())
    h5file.create_array(group, "bol_volt_ll", calibration_df["bol_volt_ll"].tolist())
    h5file.create_array(
        group, "bol_cmd_bias_rh", calibration_df["bol_cmd_bias_rh"].tolist()
    )
    h5file.create_array(
        group, "bol_cmd_bias_rl", calibration_df["bol_cmd_bias_rl"].tolist()
    )
    h5file.create_array(
        group, "bol_cmd_bias_lh", calibration_df["bol_cmd_bias_lh"].tolist()
    )
    h5file.create_array(
        group, "bol_cmd_bias_ll", calibration_df["bol_cmd_bias_ll"].tolist()
    )
    h5file.create_array(group, "time", calibration_df["time"].values)
    h5file.create_array(group, "gain_rh", calibration_df["gain_rh"].tolist())
    h5file.create_array(group, "gain_rl", calibration_df["gain_rl"].tolist())
    h5file.create_array(group, "gain_lh", calibration_df["gain_lh"].tolist())
    h5file.create_array(group, "gain_ll", calibration_df["gain_ll"].tolist())
    h5file.create_array(group, "sweeps_lh", calibration_df["sweeps_lh"].values)
    h5file.create_array(group, "sweeps_ll", calibration_df["sweeps_ll"].values)
    h5file.create_array(group, "sweeps_rh", calibration_df["sweeps_rh"].values)
    h5file.create_array(group, "sweeps_rl", calibration_df["sweeps_rl"].values)

end_time = time.time()
print(f"Time taken: {(end_time - start_time)/60} minutes")
