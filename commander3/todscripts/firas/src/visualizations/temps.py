import h5py
import matplotlib.pyplot as plt

fdq_eng = h5py.File("/mn/stornext/d16/cmbco/ola/firas/initial_data/fdq_eng_new.h5")


def hi_vs_lo_temps(element, side):
    plt.plot(
        fdq_eng[f"en_analog/grt/{side}_lo_{element}"],
        fdq_eng[f"en_analog/grt/{side}_hi_{element}"],
        "o",
        markersize=1,
    )
    plt.xlabel(f"Low {element} temperature (K)")
    plt.ylabel(f"High {element} temperature (K)")
    plt.title(f"High vs Low {element} {side} temperatures")
    plt.savefig(f"../../visualizations/hi_vs_lo_{element}_{side}_temps.png")
    plt.clf()


elements = [
    "bol_assem",
    "collimator",
    "dihedral",
    "ical",
    "refhorn",
    "skyhorn",
    "xcal_cone",
    "xcal_tip",
]
sides = ["a", "b"]

for element in elements:
    for side in sides:
        hi_vs_lo_temps(element, side)
