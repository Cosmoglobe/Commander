import h5py
import pandas as pd
import numpy as np

# opening original data files
fdq_sdf = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5")
fdq_eng = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5")

# get the different channels into separate variables
fdq_sdf_lh = fdq_sdf["fdq_sdf_lh"]
fdq_sdf_ll = fdq_sdf["fdq_sdf_ll"]
fdq_sdf_rh = fdq_sdf["fdq_sdf_rh"]
fdq_sdf_rl = fdq_sdf["fdq_sdf_rl"]

# for sci_id_rl in range(len(fdq_sdf_rl["ifg_data/ifg"])):  # rl has the least records
# try:
# sci_id_lh = list(fdq_sdf_lh["dq_data/eng_time"]).index(
#     fdq_sdf_rl["dq_data/eng_time"][sci_id_rl]
# )
# sci_id_ll = list(fdq_sdf_ll["dq_data/eng_time"]).index(
#     fdq_sdf_rl["dq_data/eng_time"][sci_id_rl]
# )
# sci_id_rh = list(fdq_sdf_rh["dq_data/eng_time"]).index(
#     fdq_sdf_rl["dq_data/eng_time"][sci_id_rl]
# )
# eng_id = list(fdq_eng["ct_head/time"]).index(
#     fdq_sdf_rl["dq_data/eng_time"][sci_id_rl]
# )
# print(sci_id_rl, sci_id_lh, sci_id_ll, sci_id_rh, eng_id)

# does each value in fdq_sdf_rl match to any in the other four datasets?
match_lh = np.isin(fdq_sdf_rl["dq_data/eng_time"], fdq_sdf_lh["dq_data/eng_time"])
match_ll = np.isin(fdq_sdf_rl["dq_data/eng_time"], fdq_sdf_ll["dq_data/eng_time"])
match_rh = np.isin(fdq_sdf_rl["dq_data/eng_time"], fdq_sdf_rh["dq_data/eng_time"])
match_eng = np.isin(fdq_sdf_rl["dq_data/eng_time"], fdq_eng["ct_head/time"])

# get the indices of the non-matching values
non_matches_lh = np.where(~match_lh)[0]
non_matches_ll = np.where(~match_ll)[0]
non_matches_rh = np.where(~match_rh)[0]
non_matches_eng = np.where(~match_eng)[0]
# get the indices of the matching values
ids_lh = np.where(match_lh)[0]
ids_ll = np.where(match_ll)[0]
ids_rh = np.where(match_rh)[0]
ids_eng = np.where(match_eng)[0]

# create a list of the sci_ids for each dataset that will have the corresponding indices in each of the datasets in relation to rl
sci_id_lh = []
sci_id_ll = []
sci_id_rh = []
sci_id_eng = []
n_non_matches_lh = 0
n_non_matches_ll = 0
n_non_matches_rh = 0
n_non_matches_eng = 0

# make the lists of sci_ids ordered with the proper indices and -1 for non-matching values
for id in range(len(fdq_sdf_rl["dq_data/eng_time"])):
    if id in non_matches_lh:
        sci_id_lh.append(-1)
        n_non_matches_lh += 1
    else:
        sci_id_lh.append(int(ids_lh[id - n_non_matches_lh]))

    if id in non_matches_ll:
        sci_id_ll.append(-1)
        n_non_matches_ll += 1
    else:
        sci_id_ll.append(int(ids_ll[id - n_non_matches_ll]))

    if id in non_matches_rh:
        sci_id_rh.append(-1)
        n_non_matches_rh += 1
    else:
        sci_id_rh.append(int(ids_rh[id - n_non_matches_rh]))

    if id in non_matches_eng:
        sci_id_eng.append(-1)
        n_non_matches_eng += 1
    else:
        sci_id_eng.append(int(ids_eng[id - n_non_matches_eng]))

id = [-1]
time = []  # change into readable time later - gmt
ifg_lh = []
ifg_ll = []
ifg_rh = []
ifg_rl = []

# check if mode is the same in all 5 arrays
mtm_length = []
mtm_speed = []
# check if coordinates are the same in all 4 arrays, check for where in the sky (galactic center, etc)
latitude = []
longitude = []

# combine temperatures so that they are only one - use their weights, use only xcal and ical only for now, rest is very small
xcal_temp = []
ical_temp = []
skyhorn_temp = []
refhorn_temp = []
bolom1_temp = []
bolom2_temp = []
bolom3_temp = []
bolom4_temp = []
collim_temp = []
dihed_temp = []

# save the proper values for each variable
# for id_rl in range(len(fdq_sdf_rl["dq_data/eng_time"])):
for id_rl in range(10):
    print(id_rl)
    if (
        sci_id_lh[id_rl] == -1
        or sci_id_ll[id_rl] == -1
        or sci_id_rh[id_rl] == -1
        or sci_id_eng[id_rl] == -1
    ):
        continue
    else:
        id_lh = sci_id_lh[id_rl]
        id_ll = sci_id_ll[id_rl]
        id_rh = sci_id_rh[id_rl]
        id_eng = sci_id_eng[id_rl]

        mtm_length_lh = fdq_sdf_lh["sci_head/mtm_length"][id_lh]
        mtm_length_ll = fdq_sdf_ll["sci_head/mtm_length"][id_ll]
        mtm_length_rh = fdq_sdf_rh["sci_head/mtm_length"][id_rh]
        mtm_length_rl = fdq_sdf_rl["sci_head/mtm_length"][id_rl]
        mtm_length_eng_lh = fdq_eng["chan/xmit_mtm_len"][id_eng, 2]
        mtm_length_eng_ll = fdq_eng["chan/xmit_mtm_len"][id_eng, 3]
        mtm_length_eng_rh = fdq_eng["chan/xmit_mtm_len"][id_eng, 0]
        mtm_length_eng_rl = fdq_eng["chan/xmit_mtm_len"][id_eng, 1]

        if (
            mtm_length_lh
            == mtm_length_ll
            == mtm_length_rh
            == mtm_length_rl
            == mtm_length_eng_lh
            == mtm_length_eng_ll
            == mtm_length_eng_rh
            == mtm_length_eng_rl
        ):
            mtm_speed_lh = fdq_sdf_lh["sci_head/mtm_speed"][id_lh]
            mtm_speed_ll = fdq_sdf_ll["sci_head/mtm_speed"][id_ll]
            mtm_speed_rh = fdq_sdf_rh["sci_head/mtm_speed"][id_rh]
            mtm_speed_rl = fdq_sdf_rl["sci_head/mtm_speed"][id_rl]
            mtm_speed_eng_lh = fdq_eng["chan/xmit_mtm_speed"][id_eng, 2]
            mtm_speed_eng_ll = fdq_eng["chan/xmit_mtm_speed"][id_eng, 3]
            mtm_speed_eng_rh = fdq_eng["chan/xmit_mtm_speed"][id_eng, 0]
            mtm_speed_eng_rl = fdq_eng["chan/xmit_mtm_speed"][id_eng, 1]

            if (
                mtm_speed_lh
                == mtm_speed_ll
                == mtm_speed_rh
                == mtm_speed_rl
                == mtm_speed_eng_lh
                == mtm_speed_eng_ll
                == mtm_speed_eng_rh
                == mtm_speed_eng_rl
            ):

                lat_lh = fdq_sdf_lh["attitude/galactic_latitude"][id_lh] / 1e4
                lat_ll = fdq_sdf_ll["attitude/galactic_latitude"][id_ll] / 1e4
                lat_rh = fdq_sdf_rh["attitude/galactic_latitude"][id_rh] / 1e4
                lat_rl = fdq_sdf_rl["attitude/galactic_latitude"][id_rl] / 1e4

                angle_threshold = 8 * np.pi / 180

                if (
                    np.abs(lat_lh - lat_rl) < angle_threshold
                    and np.abs(lat_ll - lat_rl) < angle_threshold
                    and np.abs(lat_rh - lat_rl) < angle_threshold
                ):

                    long_lh = fdq_sdf_lh["attitude/galactic_longitude"][id_lh] / 1e4
                    long_ll = fdq_sdf_ll["attitude/galactic_longitude"][id_ll] / 1e4
                    long_rh = fdq_sdf_rh["attitude/galactic_longitude"][id_rh] / 1e4
                    long_rl = fdq_sdf_rl["attitude/galactic_longitude"][id_rl] / 1e4

                    if (
                        np.abs(long_lh - long_rl) < angle_threshold
                        and np.abs(long_ll - long_rl) < angle_threshold
                        and np.abs(long_rh - long_rl) < angle_threshold
                    ):

                        id.append(id[-1] + 1)
                        time.append(fdq_sdf_rl["ct_head/gmt"][id_rl])
                        ifg_lh.append(fdq_sdf_lh["ifg_data/ifg"][id_lh])
                        ifg_ll.append(fdq_sdf_ll["ifg_data/ifg"][id_ll])
                        ifg_rh.append(fdq_sdf_rh["ifg_data/ifg"][id_rh])
                        ifg_rl.append(fdq_sdf_rl["ifg_data/ifg"][id_rl])
                        mtm_length.append(mtm_length_rl)
                        mtm_speed.append(mtm_speed_rl)

id.pop(0)
print(len(id), len(time), len(ifg_lh), len(ifg_ll), len(ifg_rh), len(ifg_rl))

np.savez(
    "./data/data_ifgs_mtm_pointing.npz",
    id=np.array(id),
    time=np.array(time),
    ifg_lh=np.array(ifg_lh),
    ifg_ll=np.array(ifg_ll),
    ifg_rh=np.array(ifg_rh),
    ifg_rl=np.array(ifg_rl),
    mtm_length=np.array(mtm_length),
    mtm_speed=np.array(mtm_speed),
)

# sci_id_lh = np.where(
#     np.isin(
#         fdq_sdf_rl["dq_data/eng_time"],
#         fdq_sdf_lh["dq_data/eng_time"],
#     )
# )
# sci_id_ll = np.where(
#     np.isin(
#         fdq_sdf_rl["dq_data/eng_time"],
#         fdq_sdf_ll["dq_data/eng_time"],
#     )
# )
# sci_id_rh = np.where(
#     np.isin(
#         fdq_sdf_rl["dq_data/eng_time"],
#         fdq_sdf_rh["dq_data/eng_time"],
#     )
# )
# eng_id = np.where(
#     np.isin(
#         fdq_sdf_rl["dq_data/eng_time"],
#         fdq_eng["ct_head/time"],
#     )
# )

# print(len(sci_id_lh[0]), len(sci_id_ll[0]), len(sci_id_rh[0]), len(eng_id[0]))

# except ValueError:
#     print("ValueError")
#     continue
