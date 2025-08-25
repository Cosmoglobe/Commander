"""
An essential part of the data processing in data_prep.py had to be rewritten
in order to relax the contraint of having a record for all channels. This is that script.
"""

import h5py
import pandas as pd
import numpy as np
import torch

# use_cuda = torch.cuda.is_available()
# device = torch.device("cuda:0" if use_cuda else "cpu")

# opening original data files
fdq_sdf = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_sdf_new.h5")
fdq_eng = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5")

# get the different channels into separate variables
fdq_sdf_lh = fdq_sdf["fdq_sdf_lh"]
fdq_sdf_ll = fdq_sdf["fdq_sdf_ll"]
fdq_sdf_rh = fdq_sdf["fdq_sdf_rh"]
fdq_sdf_rl = fdq_sdf["fdq_sdf_rl"]

# does each value in fdq_sdf_rh match to any in the other four datasets?
sci_id_lh = []
sci_id_ll = []
sci_id_rl = []
sci_id_rh = []
for i in range(len(fdq_eng["ct_head/time"])):
    # for i in range(10):
    print(i)

    match = np.where(fdq_eng["ct_head/time"][i] == fdq_sdf_lh["dq_data/eng_time"])
    if len(match[0]) == 0:
        sci_id_lh.append(-1)
    else:
        sci_id_lh.append(int(match[0]))

    match = np.where(fdq_eng["ct_head/time"][i] == fdq_sdf_ll["dq_data/eng_time"])
    if len(match[0]) == 0:
        sci_id_ll.append(-1)
    else:
        sci_id_ll.append(int(match[0]))

    match = np.where(fdq_eng["ct_head/time"][i] == fdq_sdf_rl["dq_data/eng_time"])
    if len(match[0]) == 0:
        sci_id_rl.append(-1)
    else:
        sci_id_rl.append(int(match[0]))

    match = np.where(fdq_eng["ct_head/time"][i] == fdq_sdf_rh["dq_data/eng_time"])

    if len(match[0]) == 0:
        sci_id_rh.append(-1)
    else:
        sci_id_rh.append(int(match[0]))

id = [-1]
time = []
ifg_lh = []
ifg_ll = []
ifg_rh = []
ifg_rl = []

# check if mode is the same in all 5 arrays
mtm_lengths = []
mtm_speeds = []
# check if coordinates are the same in all 4 arrays, check for where in the sky (galactic center, etc)
latitude = []
longitude = []

# combine temperatures so that they are only one - use their weights, use only xcal and ical only for now, rest is very small
xcal_temp = []
ical_temp = []
# skyhorn_temp = []
# refhorn_temp = []
# bolom1_temp = []
# bolom2_temp = []
# bolom3_temp = []
# bolom4_temp = []
# collim_temp = []
# dihed_temp = []

# save the proper values for each variable
for id_eng in range(len(fdq_eng["ct_head/time"])):
    # for id_eng in range(10):
    print(id_eng)
    if (
        sci_id_rh[id_eng] == -1
        and sci_id_lh[id_eng] == -1
        and sci_id_ll[id_eng] == -1
        and sci_id_rl[id_eng] == -1
    ):
        continue
    else:
        mtm_length_check = True

        if sci_id_lh[id_eng] != -1:
            id_lh = sci_id_lh[id_eng]
            mtm_length_lh = fdq_sdf_lh["sci_head/mtm_length"][id_lh]
            mtm_length_eng_lh = fdq_eng["chan/xmit_mtm_len"][id_eng, 2]

            if mtm_length_eng_lh != mtm_length_lh:
                mtm_length_check = False
                continue

            mtm_length = mtm_length_lh
            lh_exists = True
        else:
            lh_exists = False

        if sci_id_ll[id_eng] != -1:
            id_ll = sci_id_ll[id_eng]
            mtm_length_ll = fdq_sdf_ll["sci_head/mtm_length"][id_ll]
            mtm_length_eng_ll = fdq_eng["chan/xmit_mtm_len"][id_eng, 3]

            if mtm_length_eng_ll != mtm_length_ll:
                mtm_length_check = False
                continue

            if lh_exists and mtm_length_lh != mtm_length_ll:
                mtm_length_check = False
                continue  # technically this continue should not be needed, it only ends the iteration earlier, not needing to finish all the checks

            mtm_length = mtm_length_ll
            ll_exists = True
        else:
            ll_exists = False

        if sci_id_rl[id_eng] != -1:
            id_rl = sci_id_rl[id_eng]
            mtm_length_rl = fdq_sdf_rl["sci_head/mtm_length"][id_rl]
            mtm_length_eng_rl = fdq_eng["chan/xmit_mtm_len"][id_eng, 1]

            if mtm_length_eng_rl != mtm_length_rl:
                mtm_length_check = False
                continue

            if lh_exists and mtm_length_lh != mtm_length_rl:
                mtm_length_check = False
                continue
            if ll_exists and mtm_length_ll != mtm_length_rl:
                mtm_length_check = False
                continue

            mtm_length = mtm_length_rl
            rl_exists = True
        else:
            rl_exists = False

        if sci_id_rh[id_eng] != -1:
            id_rh = sci_id_rh[id_eng]
            mtm_length_rh = fdq_sdf_rh["sci_head/mtm_length"][id_rh]
            mtm_length_eng_rh = fdq_eng["chan/xmit_mtm_len"][id_eng, 0]

            if mtm_length_eng_rh != mtm_length_rh:
                mtm_length_check = False
                continue

            if lh_exists and mtm_length_lh != mtm_length_rh:
                mtm_length_check = False
                continue
            if ll_exists and mtm_length_ll != mtm_length_rh:
                mtm_length_check = False
                continue
            if rl_exists and mtm_length_rl != mtm_length_rh:
                mtm_length_check = False
                continue

            mtm_length = mtm_length_rh
            rh_exists = True
        else:
            rh_exists = False

        if mtm_length_check:
            mtm_speed_check = True

            if lh_exists:
                mtm_speed_lh = fdq_sdf_lh["sci_head/mtm_speed"][id_lh]
                mtm_speed_eng_lh = fdq_eng["chan/xmit_mtm_speed"][id_eng, 2]

                if mtm_speed_lh != mtm_speed_eng_lh:
                    mtm_speed_check = False
                    continue

                mtm_speed = mtm_speed_lh

            if ll_exists:
                mtm_speed_ll = fdq_sdf_ll["sci_head/mtm_speed"][id_ll]
                mtm_speed_eng_ll = fdq_eng["chan/xmit_mtm_speed"][id_eng, 3]

                if mtm_speed_ll != mtm_speed_eng_ll:
                    mtm_speed_check = False
                    continue

                if lh_exists and mtm_speed_lh != mtm_speed_ll:
                    mtm_speed_check = False
                    continue

                mtm_speed = mtm_speed_ll

            if rl_exists:
                mtm_speed_rl = fdq_sdf_rl["sci_head/mtm_speed"][id_rl]
                mtm_speed_eng_rl = fdq_eng["chan/xmit_mtm_speed"][id_eng, 1]

                if mtm_speed_rl != mtm_speed_eng_rl:
                    mtm_speed_check = False
                    continue

                if lh_exists and mtm_speed_lh != mtm_speed_rl:
                    mtm_speed_check = False
                    continue
                if ll_exists and mtm_speed_ll != mtm_speed_rl:
                    mtm_speed_check = False
                    continue

                mtm_speed = mtm_speed_rl

            if rh_exists:
                mtm_speed_rh = fdq_sdf_rh["sci_head/mtm_speed"][id_rh]
                mtm_speed_eng_rh = fdq_eng["chan/xmit_mtm_speed"][id_eng, 0]

                if mtm_speed_rh != mtm_speed_eng_rh:
                    mtm_speed_check = False
                    continue

                if lh_exists and mtm_speed_lh != mtm_speed_rh:
                    mtm_speed_check = False
                    continue
                if ll_exists and mtm_speed_ll != mtm_speed_rh:
                    mtm_speed_check = False
                    continue
                if rl_exists and mtm_speed_rl != mtm_speed_rh:
                    mtm_speed_check = False
                    continue

                mtm_speed = mtm_speed_rh

            if mtm_speed_check:
                # lat_lh = fdq_sdf_lh["attitude/galactic_latitude"][id_lh] / 1e4
                # lat_ll = fdq_sdf_ll["attitude/galactic_latitude"][id_ll] / 1e4
                # lat_rh = fdq_sdf_rh["attitude/galactic_latitude"][id_rh] / 1e4
                # lat_rl = fdq_sdf_rl["attitude/galactic_latitude"][id_rl] / 1e4

                # angle_threshold = 8 * np.pi / 180

                # if (
                #     np.abs(lat_lh - lat_rh) < angle_threshold
                #     and np.abs(lat_ll - lat_rh) < angle_threshold
                #     and np.abs(lat_rl - lat_rh) < angle_threshold
                # ):

                #     long_lh = fdq_sdf_lh["attitude/galactic_longitude"][id_lh] / 1e4
                #     long_ll = fdq_sdf_ll["attitude/galactic_longitude"][id_ll] / 1e4
                #     long_rh = fdq_sdf_rh["attitude/galactic_longitude"][id_rh] / 1e4
                #     long_rl = fdq_sdf_rl["attitude/galactic_longitude"][id_rl] / 1e4

                #     if (
                #         np.abs(long_lh - long_rh) < angle_threshold
                #         and np.abs(long_ll - long_rh) < angle_threshold
                #         and np.abs(long_rl - long_rh) < angle_threshold
                #     ):

                id.append(id[-1] + 1)
                time.append(fdq_sdf_rh["ct_head/gmt"][id_eng])

                if lh_exists:
                    ifg_lh.append(fdq_sdf_lh["ifg_data/ifg"][id_lh])
                else:
                    ifg_lh.append(np.zeros(512))
                if ll_exists:
                    ifg_ll.append(fdq_sdf_ll["ifg_data/ifg"][id_ll])
                else:
                    ifg_ll.append(np.zeros(512))
                if rh_exists:
                    ifg_rh.append(fdq_sdf_rh["ifg_data/ifg"][id_rh])
                else:
                    ifg_rh.append(np.zeros(512))
                if rl_exists:
                    ifg_rl.append(fdq_sdf_rl["ifg_data/ifg"][id_rl])
                else:
                    ifg_rl.append(np.zeros(512))

                mtm_lengths.append(mtm_length)
                mtm_speeds.append(mtm_speed)

id.pop(0)
print(len(id), len(time), len(ifg_lh), len(ifg_ll), len(ifg_rh), len(ifg_rl))

np.savez(
    "./data/data_v1.npz",
    id=np.array(id),
    time=np.array(time),
    ifg_lh=np.array(ifg_lh),
    ifg_ll=np.array(ifg_ll),
    ifg_rh=np.array(ifg_rh),
    ifg_rl=np.array(ifg_rl),
    mtm_length=np.array(mtm_lengths),
    mtm_speed=np.array(mtm_speeds),
)
