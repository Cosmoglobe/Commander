import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import data.utils.my_utils as data_utils
import globals as g


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


if __name__ == "__main__":
    fdq_sdf = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5")
    gmt, midpoint_time, space_time, t, eng_time, df = {}, {}, {}, {}, {}, {}

    for channel in g.CHANNELS.keys():
        print(f"Channel {channel}:")
        gmt[channel] = fdq_sdf[f"/fdq_sdf_{channel}/ct_head/gmt"][:].astype(str)
        print(f"gmt[channel] shape: {gmt[channel].shape}")
        print(f"gmt[channel] sample: {gmt[channel][0]}")
        # parse to datetime
        parse_vectorized = np.vectorize(data_utils.parse_date_string)
        gmt[channel] = parse_vectorized(gmt[channel])
        print(f"gmt: {gmt[channel]}")
        midpoint_time[channel] = binary_to_gmt(
            fdq_sdf[f"/fdq_sdf_{channel}/collect_time/midpoint_time"][:]
        )
        print(f"midpoint_time: {midpoint_time[channel]}")
        # space_time[channel] = fdq_sdf[f"/fdq_sdf_{channel}/ct_head/space_time"][:]
        # print(f"space_time: {space_time[channel]}")
        t[channel] = binary_to_gmt(fdq_sdf[f"/fdq_sdf_{channel}/ct_head/time"][:])
        print(f"time: {t[channel]}")
        eng_time[channel] = fdq_sdf[f"/fdq_sdf_{channel}/dq_data/eng_time"][:]
        print(f"eng_time: {eng_time[channel]}")

        df[channel] = pd.DataFrame(
            {
                "gmt": gmt[channel],
                "midpoint_time": midpoint_time[channel],
                "time": t[channel],
                "eng_time": eng_time[channel],
            }
        )

        # Filter out zero eng_time and add channel suffix to all columns except eng_time
        df[channel] = df[channel][df[channel]["eng_time"] != 0].reset_index(drop=True)
        df[channel] = df[channel].add_suffix(f"_{channel}")
        df[channel].rename(columns={f"eng_time_{channel}": "eng_time"}, inplace=True)

    # match 4 channels together based on eng_time
    merged_df = df["rh"]
    for channel in ["rl", "lh", "ll"]:
        merged_df = pd.merge(merged_df, df[channel], on="eng_time", how="inner")
    print(merged_df)

    print(merged_df.columns)

    for channel in g.CHANNELS.keys():
        print(f"Channel {channel}")
        print(merged_df[f"gmt_{channel}"].values)

    # take trhe difference between each of the channels to see how far they are apart from each other and plot
    merged_df["gmt_diff_rl"] = (
        merged_df["gmt_rl"] - merged_df["gmt_rh"]
    ).dt.total_seconds()
    merged_df["gmt_diff_lh"] = (
        merged_df["gmt_lh"] - merged_df["gmt_rh"]
    ).dt.total_seconds()
    merged_df["gmt_diff_ll"] = (
        merged_df["gmt_ll"] - merged_df["gmt_rh"]
    ).dt.total_seconds()
    plt.plot(merged_df["gmt_diff_rl"], label="rl - rh")
    plt.plot(merged_df["gmt_diff_lh"], label="lh - rh")
    plt.plot(merged_df["gmt_diff_ll"], label="ll - rh")
    plt.legend()
    plt.xlabel("Sample Index")
    plt.ylabel("GMT Difference (seconds)")
    plt.title("GMT Differences Between Channels")
    plt.savefig("eda/output/gmt_differences.png")
    plt.close()

    # now for midpoint_time
    merged_df["midpoint_time_diff_rl"] = (
        merged_df["midpoint_time_rl"] - merged_df["midpoint_time_rh"]
    ).dt.total_seconds()
    merged_df["midpoint_time_diff_lh"] = (
        merged_df["midpoint_time_lh"] - merged_df["midpoint_time_rh"]
    ).dt.total_seconds()
    merged_df["midpoint_time_diff_ll"] = (
        merged_df["midpoint_time_ll"] - merged_df["midpoint_time_rh"]
    ).dt.total_seconds()
    plt.plot(merged_df["midpoint_time_diff_rl"], label="rl - rh")
    plt.plot(merged_df["midpoint_time_diff_lh"], label="lh - rh")
    plt.plot(merged_df["midpoint_time_diff_ll"], label="ll - rh")
    plt.legend()
    plt.xlabel("Sample Index")
    plt.ylabel("Midpoint Time Difference (seconds)")
    plt.title("Midpoint Time Differences Between Channels")
    plt.savefig("eda/output/midpoint_time_differences.png")
    plt.close()

    # now for the binary time
    merged_df["time_diff_rl"] = (
        merged_df["time_rl"] - merged_df["time_rh"]
    ).dt.total_seconds()
    merged_df["time_diff_lh"] = (
        merged_df["time_lh"] - merged_df["time_rh"]
    ).dt.total_seconds()
    merged_df["time_diff_ll"] = (
        merged_df["time_ll"] - merged_df["time_rh"]
    ).dt.total_seconds()
    plt.plot(merged_df["time_diff_rl"], label="rl - rh")
    plt.plot(merged_df["time_diff_lh"], label="lh - rh")
    plt.plot(merged_df["time_diff_ll"], label="ll - rh")
    plt.legend()
    plt.xlabel("Sample Index")
    plt.ylabel("Binary Time Difference (seconds)")
    plt.title("Binary Time Differences Between Channels")
    plt.savefig("eda/output/binary_time_differences.png")
    plt.close()

    # merged_by_gmt = df["rh"]
    # for channel in ["rl", "lh", "ll"]:
    #     merged_by_gmt = pd.merge(
    #         merged_by_gmt,
    #         df[channel],
    #         on=f"gmt",
    #         how="inner",
    #     )
    # print(merged_by_gmt)
