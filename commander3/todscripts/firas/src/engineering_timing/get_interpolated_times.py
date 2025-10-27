import astropy.units as u
import h5py
import numpy as np

# DATA_DIR = '/home/dwatts/Commander/commander3/todscripts/firas/src/engineering_timing'
DATA_DIR = "/mn/stornext/d16/cmbco/ola/firas/initial_data/"

plot = True

lens_grt = [
    1,
    1,
    1,
    1,
    1,
    4,
    1,
    4,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    4,
    1,
    4,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    4,
    1,
    4,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    4,
    1,
    4,
    1,
    1,
]

names_en_analog_grt = [
    # grts
    # 'a_lo_grt', 'a_hi_grt', 'b_lo_grt', 'b_hi_grt',
    "a_lo_xcal_tip",
    "a_lo_skyhorn",
    "a_lo_refhorn",
    "a_lo_ical",
    "a_lo_dihedral",
    "a_lo_bol_assem",
    "a_lo_mirror",
    "a_lo_cal_resistors",
    "a_lo_xcal_cone",
    "a_lo_collimator",
    "a_hi_xcal_tip",
    "a_hi_skyhorn",
    "a_hi_refhorn",
    "a_hi_ical",
    "a_hi_dihedral",
    "a_hi_bol_assem",
    "a_hi_mirror",
    "a_hi_cal_resistors",
    "a_hi_xcal_cone",
    "a_hi_collimator",
    "b_lo_xcal_tip",
    "b_lo_skyhorn",
    "b_lo_refhorn",
    "b_lo_ical",
    "b_lo_dihedral",
    "b_lo_bol_assem",
    "b_lo_mirror",
    "b_lo_cal_resistors",
    "b_lo_xcal_cone",
    "b_lo_collimator",
    "b_hi_xcal_tip",
    "b_hi_skyhorn",
    "b_hi_refhorn",
    "b_hi_ical",
    "b_hi_dihedral",
    "b_hi_bol_assem",
    "b_hi_mirror",
    "b_hi_cal_resistors",
    "b_hi_xcal_cone",
    "b_hi_collimator",
]


def get_data():
    data = h5py.File(f"{DATA_DIR}/fdq_eng.h5")
    eng = h5py.File(f"{DATA_DIR}/fdq_eng_new.h5")
    sdf = h5py.File(f"{DATA_DIR}/fdq_sdf_new.h5")

    grt_arr = data["fdq_eng"]["en_analog"]["grt"]

    grts = {}
    offsets = {}

    # "Binary time"
    # ADT Time, in units of 100ns since 1858-11-17
    bin_time = data["fdq_eng"]["ct_head"]["time"] * (100 * u.ns)

    reord = np.argsort(bin_time)

    stat_word_9 = eng["en_stat/stat_word_9"][()]
    filters = stat_word_9 == 16185

    filters = filters[reord]

    ind0 = 0
    print("getting offsets")
    for i in range(len(names_en_analog_grt)):
        grts[names_en_analog_grt[i]] = grt_arr[:, ind0 : ind0 + lens_grt[i]].flatten()
        grts[names_en_analog_grt[i]] = grts[names_en_analog_grt[i]][reord]
        offsets[names_en_analog_grt[i]] = ind0
        ind0 += lens_grt[i]
    time = bin_time[reord]

    inds = (np.arange(len(time)) > 0) & (np.arange(len(time)) < 580_000)
    inds = (np.arange(len(time)) > 100_000) & (np.arange(len(time)) < 125_000)
    time = time.to("s")

    grt_times = {}

    grt_names = [
        "xcal_tip",
        "skyhorn",
        "refhorn",
        "ical",
        "dihedral",
        "mirror",
        "xcal_cone",
        "collimator",
    ]
    offset_ind = [0, 1, 2, 3, 4, 9, 14, 15]

    print("getting the grt times")
    ifg_time_offset = 36 * u.s
    for side in ["a", "b"]:
        for i, grt in enumerate(grt_names):
            not_ok = grts[f"{side}_lo_{grt}"] == -9999
            grts[f"{side}_lo_{grt}"][not_ok] = np.nan
            grts[f"{side}_hi_{grt}"][not_ok] = np.nan
            not_ok = grts[f"{side}_hi_{grt}"] == -9999
            grts[f"{side}_lo_{grt}"][not_ok] = np.nan
            grts[f"{side}_hi_{grt}"][not_ok] = np.nan

            grts[f"{side}_lo_{grt}"][filters] = np.nan
            grts[f"{side}_hi_{grt}"][filters] = np.nan

            # This ordering is an artifact of the time-domain multiplexing, likely wires being plugged in a non-ideal order.
            if (grt == "collimator") or (grt == "xcal_cone"):
                grt_times[f"{side}_lo_{grt}"] = time + 16 * u.s + ifg_time_offset
                grt_times[f"{side}_hi_{grt}"] = time + ifg_time_offset
            else:
                grt_times[f"{side}_lo_{grt}"] = time + ifg_time_offset
                grt_times[f"{side}_hi_{grt}"] = time + 16 * u.s + ifg_time_offset

            grt_times[f"{side}_lo_{grt}"] += offset_ind[i] * u.s
            grt_times[f"{side}_hi_{grt}"] += offset_ind[i] * u.s
    return grts, grt_times


def test_plot():
    import matplotlib.pyplot as plt

    grt = "ical"
    current = "lo"
    side = "a"

    plt.plot(
        grt_times[f"{side}_{current}_{grt}"][inds],
        grts[f"{side}_{current}_{grt}"][inds],
    )

    current = "hi"
    plt.plot(
        grt_times[f"{side}_{current}_{grt}"][inds],
        grts[f"{side}_{current}_{grt}"][inds],
    )

    grt = "ical"
    current = "lo"
    side = "a"

    plt.figure()

    plt.plot(
        grt_times[f"{side}_{current}_{grt}"][inds],
        grts[f"{side}_{current}_{grt}"][inds],
        "o",
    )

    t = grt_times[f"{side}_{current}_{grt}"][inds]
    t_lin = np.linspace(t.min(), t.max(), 10 * len(t))
    plt.plot(t, interpolators[f"{side}_{current}_{grt}"](t))
    plt.show()


def get_interp(grts, grt_times):
    # Get interpolator dictionary
    from scipy.interpolate import interp1d

    interpolators = {}
    for key in grt_times.keys():
        t = grt_times[key]
        T = grts[key]
        f = interp1d(t, T)
        interpolators[key] = f
    return interpolators


if __name__ == "__main__":
    test_plot()
