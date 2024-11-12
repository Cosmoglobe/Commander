"""
Data processing script has to be re-written again, not matching eng_time/time now,
and instead interpolating over gmt to get the engineering data corresponding to the science data.
"""

import h5py
import numpy as np
import pandas as pd
import tables as tb
from scipy import interpolate
import time

from utils.my_utils import (
    parse_date_string,
    clean_variable,
)

# check how much time the script takes to run
start_time = time.time()

# OPENING ORIGINAL DATA FILES
fdq_sdf = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5")
fdq_eng = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5")

# PARSING THE H5 DATA INTO A PANDAS DATAFRAME TO MANIPULATE DATA EASIER
# decoding the gmt data
gmt_lh = np.array(fdq_sdf["fdq_sdf_lh/ct_head/gmt"]).astype(str)  # .astype(int)
gmt_ll = np.array(fdq_sdf["fdq_sdf_ll/ct_head/gmt"]).astype(str)  # .astype(int)
gmt_rh = np.array(fdq_sdf["fdq_sdf_rh/ct_head/gmt"]).astype(str)  # .astype(int)
gmt_rl = np.array(fdq_sdf["fdq_sdf_rl/ct_head/gmt"]).astype(str)  # .astype(int)

gmt_lh_parsed = []
gmt_ll_parsed = []
gmt_rh_parsed = []
gmt_rl_parsed = []

for gmt_nb in gmt_lh:
    gmt_lh_parsed.append(parse_date_string(gmt_nb))
for gmt_nb in gmt_ll:
    gmt_ll_parsed.append(parse_date_string(gmt_nb))
for gmt_nb in gmt_rh:
    # gmt_rh_parsed.append(parse_date_string((gmt_nb.astype(int) - 250).astype(str)))
    gmt_rh_parsed.append(parse_date_string(gmt_nb))
for gmt_nb in gmt_rl:
    # gmt_rl_parsed.append(parse_date_string((gmt_nb.astype(int) - 250).astype(str)))
    gmt_rl_parsed.append(parse_date_string(gmt_nb))

# getting the ifg data
ifg_lh = fdq_sdf["fdq_sdf_lh/ifg_data/ifg"]
ifg_ll = fdq_sdf["fdq_sdf_ll/ifg_data/ifg"]
ifg_rh = fdq_sdf["fdq_sdf_rh/ifg_data/ifg"]
ifg_rl = fdq_sdf["fdq_sdf_rl/ifg_data/ifg"]

# getting the xcal_pos data
xcal_pos_lh = fdq_sdf["fdq_sdf_lh/dq_data/xcal_pos"]
xcal_pos_ll = fdq_sdf["fdq_sdf_ll/dq_data/xcal_pos"]
xcal_pos_rh = fdq_sdf["fdq_sdf_rh/dq_data/xcal_pos"]
xcal_pos_rl = fdq_sdf["fdq_sdf_rl/dq_data/xcal_pos"]

# mtm length
mtm_length_lh = fdq_sdf["fdq_sdf_lh/sci_head/mtm_length"]
mtm_length_ll = fdq_sdf["fdq_sdf_ll/sci_head/mtm_length"]
mtm_length_rh = fdq_sdf["fdq_sdf_rh/sci_head/mtm_length"]
mtm_length_rl = fdq_sdf["fdq_sdf_rl/sci_head/mtm_length"]

# mtm speed
mtm_speed_lh = fdq_sdf["fdq_sdf_lh/sci_head/mtm_speed"]
mtm_speed_ll = fdq_sdf["fdq_sdf_ll/sci_head/mtm_speed"]
mtm_speed_rh = fdq_sdf["fdq_sdf_rh/sci_head/mtm_speed"]
mtm_speed_rl = fdq_sdf["fdq_sdf_rl/sci_head/mtm_speed"]

# fake-it mode
fake_lh = fdq_sdf["fdq_sdf_lh/dq_data/fake"]
fake_ll = fdq_sdf["fdq_sdf_ll/dq_data/fake"]
fake_rh = fdq_sdf["fdq_sdf_rh/dq_data/fake"]
fake_rl = fdq_sdf["fdq_sdf_rl/dq_data/fake"]

# upmode - science data is taken in mode 4
upmode_lh = fdq_sdf["fdq_sdf_lh/sci_head/sc_head1a"]
upmode_ll = fdq_sdf["fdq_sdf_ll/sci_head/sc_head1a"]
upmode_rh = fdq_sdf["fdq_sdf_rh/sci_head/sc_head1a"]
upmode_rl = fdq_sdf["fdq_sdf_rl/sci_head/sc_head1a"]

# adds per group
adds_per_group_lh = fdq_sdf["fdq_sdf_lh/sci_head/sc_head9"]
adds_per_group_ll = fdq_sdf["fdq_sdf_ll/sci_head/sc_head9"]
adds_per_group_rh = fdq_sdf["fdq_sdf_rh/sci_head/sc_head9"]
adds_per_group_rl = fdq_sdf["fdq_sdf_rl/sci_head/sc_head9"]

# getting each channel into its own df
df_lh = pd.DataFrame(
    {
        "gmt": gmt_lh_parsed,
        "ifg": list(ifg_lh),
        "xcal_pos": list(xcal_pos_lh),
        "mtm_length": list(mtm_length_lh),
        "mtm_speed": list(mtm_speed_lh),
        "fake": list(fake_lh),
        "upmode": list(upmode_lh),
        "adds_per_group": list(adds_per_group_lh),
    }
).sort_values("gmt")
df_ll = pd.DataFrame(
    {
        "gmt": gmt_ll_parsed,
        "ifg": list(ifg_ll),
        "xcal_pos": list(xcal_pos_ll),
        "mtm_length": list(mtm_length_ll),
        "mtm_speed": list(mtm_speed_ll),
        "fake": list(fake_ll),
        "upmode": list(upmode_ll),
        "adds_per_group": list(adds_per_group_ll),
    }
).sort_values("gmt")
df_rh = pd.DataFrame(
    {
        "gmt": gmt_rh_parsed,
        "ifg": list(ifg_rh),
        "xcal_pos": list(xcal_pos_rh),
        "mtm_length": list(mtm_length_rh),
        "mtm_speed": list(mtm_speed_rh),
        "fake": list(fake_rh),
        "upmode": list(upmode_rh),
        "adds_per_group": list(adds_per_group_rh),
    }
).sort_values("gmt")
df_rl = pd.DataFrame(
    {
        "gmt": gmt_rl_parsed,
        "ifg": list(ifg_rl),
        "xcal_pos": list(xcal_pos_rl),
        "mtm_length": list(mtm_length_rl),
        "mtm_speed": list(mtm_speed_rl),
        "fake": list(fake_rl),
        "upmode": list(upmode_rl),
        "adds_per_group": list(adds_per_group_rl),
    }
).sort_values("gmt")

# USING THIS FOR MERGE_ASOF
tolerance = pd.Timedelta(seconds=16)
# getting all possible gmts so that we can do an outer join using merge_asof
unified_timestamps = (
    pd.DataFrame(pd.concat([df_lh["gmt"], df_ll["gmt"], df_rh["gmt"], df_rl["gmt"]]))
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

# # outer-join the dataframes on gmt
merged_df = pd.merge_asof(
    non_duplicate_timestamps, df_lh, on="gmt", direction="nearest", tolerance=tolerance
)

merged_df = pd.merge_asof(
    # df_lh,
    merged_df,
    df_ll,
    on="gmt",
    # how="outer",
    direction="nearest",
    suffixes=("_lh", "_ll"),
    tolerance=tolerance,
)
# merged_df = pd.merge(
merged_df = pd.merge_asof(
    merged_df,
    df_rh,
    on="gmt",
    # how="outer",
    direction="nearest",
    suffixes=("_ll", "_rh"),
    tolerance=tolerance,
)
# merged_df = pd.merge(
merged_df = pd.merge_asof(
    merged_df,
    df_rl,
    on="gmt",
    # how="outer",
    direction="nearest",
    suffixes=("_rh", "_rl"),
    tolerance=tolerance,
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

zero_list = [0] * 512
# filling out ifg nans with zeros
for column in ["ifg_lh", "ifg_ll", "ifg_rh", "ifg_rl"]:
    merged_df[column] = merged_df[column].apply(
        lambda x: zero_list if (isinstance(x, float) and np.isnan(x)) else x
    )

# CLEANING XCAL_POS AND CONSTRAINING FOR ONLY 1 AND 2
# making sure xcal_pos is the same within each record
# merged_df["xcal_pos"] = merged_df.apply(clean_xcal_pos, axis=1)
merged_df["xcal_pos"] = merged_df.apply(clean_variable, axis=1, args=("xcal_pos",))
merged_df = merged_df.drop(
    columns=["xcal_pos_lh", "xcal_pos_ll", "xcal_pos_rh", "xcal_pos_rl"]
)
merged_df = merged_df[(merged_df["xcal_pos"] == 1) | (merged_df["xcal_pos"] == 2)]

# same thing but for mtm length and speed
# merged_df["mtm_length"] = merged_df.apply(clean_mtm_length, axis=1)
merged_df["mtm_length"] = merged_df.apply(clean_variable, axis=1, args=("mtm_length",))
# merged_df["mtm_speed"] = merged_df.apply(clean_mtm_speed, axis=1)
merged_df["mtm_speed"] = merged_df.apply(clean_variable, axis=1, args=("mtm_speed",))
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
merged_df["fake"] = merged_df.apply(clean_variable, axis=1, args=("fake",))
# fake-it mode on is 1, so i am guessing that 0 is when it is off
merged_df = merged_df[(merged_df["fake"] == 0)]
merged_df = merged_df.drop(columns=["fake", "fake_lh", "fake_ll", "fake_rh", "fake_rl"])

# clean up for upmode
merged_df["upmode"] = merged_df.apply(clean_variable, axis=1, args=("upmode",))
merged_df = merged_df[(merged_df["upmode"] == 4)]
merged_df = merged_df.drop(
    columns=["upmode", "upmode_lh", "upmode_ll", "upmode_rh", "upmode_rl"]
)

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

# reset index
merged_df = merged_df.reset_index(drop=True)

print(f"Dataframe after cleaning science data: {merged_df.tail()}")

# initializing the xcal temp column so then we only add the relevant xcal values
# merged_df["xcal"][merged_df["xcal_pos"] == 2] = np.nan
# merged_df["xcal"][merged_df["xcal_pos"] == 1] = 0

# ENGINEERING DATA
a_hi_ical = fdq_eng["en_analog/grt/a_hi_ical"]
a_lo_ical = fdq_eng["en_analog/grt/a_lo_ical"]
b_hi_ical = fdq_eng["en_analog/grt/b_hi_ical"]
b_lo_ical = fdq_eng["en_analog/grt/b_lo_ical"]
ical = np.mean([a_hi_ical, a_lo_ical, b_hi_ical, b_lo_ical], axis=0)

a_hi_xcal_cone = np.array(fdq_eng["en_analog/grt/a_hi_xcal_cone"])
a_hi_xcal_tip = np.array(fdq_eng["en_analog/grt/a_hi_xcal_tip"])
a_lo_xcal_cone = np.array(fdq_eng["en_analog/grt/a_lo_xcal_cone"])
a_lo_xcal_tip = np.array(fdq_eng["en_analog/grt/a_lo_xcal_tip"])
b_hi_xcal_cone = np.array(fdq_eng["en_analog/grt/b_hi_xcal_cone"])
b_lo_xcal_cone = np.array(fdq_eng["en_analog/grt/b_lo_xcal_cone"])
xcal = np.mean(
    [
        a_hi_xcal_cone * 0.9 + a_hi_xcal_tip * 0.1,
        a_lo_xcal_cone * 0.9 + a_lo_xcal_tip * 0.1,
        b_hi_xcal_cone,
        b_lo_xcal_cone,
    ],
    axis=0,
)

gmt_eng = np.array(fdq_eng["ct_head/gmt"]).astype(str)
gmt_eng_parsed = []
for gmt in gmt_eng:
    gmt_eng_parsed.append(parse_date_string(gmt))

# make engineeering data df
df_eng = pd.DataFrame(
    {
        "gmt": gmt_eng_parsed,
        "ical": list(ical),
        "xcal": list(xcal),
    }
)

# precompute timestamps for merged_df and df_eng
science_times = merged_df["gmt"].apply(lambda x: x.timestamp()).values
engineering_times = df_eng["gmt"].apply(lambda x: x.timestamp()).values

# initialize list for interpolated values
ical_new = []
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
            df_eng["ical"][indices],
            fill_value="extrapolate",
        )
        ical_new.append(f(target_time))
        # xcal_pos_series = merged_df[merged_df["gmt"] == target_time]["xcal_pos"]
        # if not xcal_pos_series.empty and xcal_pos_series.iloc[0] == 1: # FIX TO NOT GET XCAL TEMP FOR ALL
        f = interpolate.interp1d(
            engineering_times[indices],
            df_eng["xcal"][indices],
            fill_value="extrapolate",
        )
        xcal_new.append(f(target_time))
        # elif not xcal_pos_series.empty and xcal_pos_series.iloc[0] == 2:
        #     xcal_new.append(np.nan)
    else:
        ical_new.append(np.nan)
        xcal_new.append(np.nan)

merged_df["ical"] = ical_new
merged_df["xcal"] = xcal_new

# drop rows without ical data
merged_df = merged_df[merged_df["ical"].notna()]

# set xcal to nan if xcal_pos is 2
# merged_df[merged_df["xcal_pos"] == 2]["xcal"] = np.nan

# converting gmt to string so we can save
gmt_str = merged_df["gmt"].dt.strftime("%Y-%m-%d %H:%M:%S").to_numpy(dtype="S")

print(f"Dataframe after merging engineering data: {merged_df.tail()}")

# saving to a h5 file
with tb.open_file("./../../data/df_v10.h5", mode="w") as h5file:
    group = h5file.create_group("/", "df_data", "Merged Data")

    h5file.create_array(group, "gmt", gmt_str)
    h5file.create_array(group, "ifg_lh", np.stack(merged_df["ifg_lh"].values))
    h5file.create_array(group, "ifg_ll", np.stack(merged_df["ifg_ll"].values))
    h5file.create_array(group, "ifg_rh", np.stack(merged_df["ifg_rh"].values))
    h5file.create_array(group, "ifg_rl", np.stack(merged_df["ifg_rl"].values))
    h5file.create_array(group, "xcal_pos", merged_df["xcal_pos"].values)
    h5file.create_array(group, "ical", merged_df["ical"].tolist())
    h5file.create_array(group, "xcal", merged_df["xcal"].tolist())
    h5file.create_array(group, "mtm_length", merged_df["mtm_length"].values)
    h5file.create_array(group, "mtm_speed", merged_df["mtm_speed"].values)
    # h5file.create_array(group, "fake", merged_df["fake"].values)
    # h5file.create_array(group, "upmode", merged_df["upmode"].values)
    h5file.create_array(group, "adds_per_group", merged_df["adds_per_group"].values)

end_time = time.time()
print(f"Time taken: {(end_time - start_time)/60} minutes")
