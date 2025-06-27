"""
This script takes the interferograms into spectra.
"""

import os
import sys
from datetime import datetime

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import globals as g
import my_utils as mu
from utils.config import gen_nyquistl

channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}
modes = {"ss": 0, "lf": 3}  # can change when i have the new cal models

cal_data = h5py.File(
    g.PREPROCESSED_DATA_PATH_CAL,
    "r",
)


print("Getting variables")
variable_names = [
    "gmt",
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
    "mtm_length",
    "mtm_speed",
    "a_bol_assem_rh",
    "b_bol_assem_rh",
    "a_bol_assem_rl",
    "b_bol_assem_rl",
    "a_bol_assem_lh",
    "b_bol_assem_lh",
    "a_bol_assem_ll",
    "b_bol_assem_ll",
    "stat_word_1",
    "stat_word_5",
    "stat_word_9",
    "stat_word_12",
    "stat_word_13",
    "stat_word_16",
    "lvdt_stat_a",
    "lvdt_stat_b",
    "adds_per_group",
    "sweeps",
    "a_xcal",
    "b_xcal",
]
channel_dependent = [
    "ifg",
    "bol_cmd_bias",
    "bol_volt",
    "gain"
]

for channel in channels:
    for variable_name in channel_dependent:
        variable_names = variable_names + [f"{variable_name}_{channel}"]

print("Variable names:", variable_names)

variables = {}
for variable_name in variable_names:
    variables[variable_name] = np.array(cal_data["df_data/" + variable_name][()])


variables["gmt"] = variables["gmt"].astype(str)
variables["gmt"] = np.array(
    [datetime.strptime(gmt, "%Y-%m-%d %H:%M:%S") for gmt in variables["gmt"]]
)

filter_bad = mu.filter_junk(
    variables["stat_word_1"],
    variables["stat_word_5"],
    variables["stat_word_9"],
    variables["stat_word_12"],
    variables["stat_word_13"],
    variables["stat_word_16"],
    variables["lvdt_stat_a"],
    variables["lvdt_stat_b"],
    variables["a_bol_assem_rh"],
    variables["a_bol_assem_rl"],
    variables["a_bol_assem_lh"],
    variables["b_bol_assem_lh"],
    variables["a_bol_assem_ll"],
    variables["b_bol_assem_ll"],
    variables["bol_cmd_bias_lh"],
    variables["bol_cmd_bias_rh"],
)

for variable in variables.keys():
    variables[variable] = variables[variable][filter_bad]

variables["ical"] = (variables["a_ical"] + variables["b_ical"]) / 2
variables["xcal"] = (variables["a_xcal"] + variables["b_xcal"]) / 2
variables["dihedral"] = (variables["a_dihedral"] + variables["b_dihedral"]) / 2
variables["refhorn"] = (variables["a_refhorn"] + variables["b_refhorn"]) / 2
variables["skyhorn"] = (variables["a_skyhorn"] + variables["b_skyhorn"]) / 2
variables["collimator"] = (variables["a_collimator"] + variables["b_collimator"]) / 2
variables["bolometer_rh"] = (
    variables["a_bol_assem_rh"] + variables["b_bol_assem_rh"]
) / 2
variables["bolometer_rl"] = (
    variables["a_bol_assem_rl"] + variables["b_bol_assem_rl"]
) / 2
variables["bolometer_lh"] = (
    variables["a_bol_assem_lh"] + variables["b_bol_assem_lh"]
) / 2
variables["bolometer_ll"] = (
    variables["a_bol_assem_ll"] + variables["b_bol_assem_ll"]
) / 2

# filtering the data based on the mode used
short_filter = variables["mtm_length"] == 0
long_filter = variables["mtm_length"] == 1
slow_filter = variables["mtm_speed"] == 0
fast_filter = variables["mtm_speed"] == 1

filters = {}
filters["ss"] = short_filter & slow_filter
filters["lf"] = long_filter & fast_filter

fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

variablesm = {}
for variable in variables.keys():
    for mode, mode_value in modes.items():
            variablesm[f"{variable}_{mode}"] = variables[variable][filters[mode]]

cal = {}
spec = {}
for mode, mode_value in modes.items():
    for channel, channel_value in channels.items():
        if not(mode == "lf" and (channel == "lh" or channel == "rh")):
            fits_data = fits.open(
                f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
            )
            frec = 4 * (channel_value % 2) + mode_value


            # optical transfer function
            otf = (
                fits_data[1].data["RTRANSFE"][0]
                + 1j * fits_data[1].data["ITRANSFE"][0]
            )
            otf = otf[
                np.abs(otf) > 0
            ]

            apod = fits_data[1].data[
                "APODIZAT"
            ][0]
            etf = (
                fits_data[1].data["RELEX_GA"][0]
                + 1j * fits_data[1].data["IELEX_GA"][0]
            )
            etf = etf[
                np.abs(etf) > 0
            ]

             # bolometer parameters
            S0 = fits_data[1].data[
                "DC_RESPO"
            ][0]
            tau = fits_data[1].data[
                "TIME_CON"
            ][0]
            Jo = fits_data[1].data[
                "BOLPARM8"
            ][0]
            Jg = fits_data[1].data[
                "BOLPARM9"
            ][0]
            T0 = fits_data[1].data[
                "BOLPARM2"
            ][0]
            R0 = fits_data[1].data[
                "BOLPARM_"
            ][0]
            rho = fits_data[1].data[
                "BOLPARM5"
            ][0]
            G1 = fits_data[1].data[
                "BOLPARM3"
            ][0]
            beta = fits_data[1].data[
                "BOLPARM4"
            ][0]
            C3 = fits_data[1].data[
                "BOLPARM7"
            ][0]
            C1 = fits_data[1].data[
                "BOLPARM6"
            ][0]

            print(f"Converting interferograms to spectra for {channel.upper()}{mode.upper()}")
            _, spec[f"spec_{channel}_{mode}"] = mu.ifg_to_spec(
                ifg=variablesm[f"ifg_{channel}_{mode}"],
                mtm_speed=0 if mode[1] == "s" else 1,
                channel=channel_value,
                adds_per_group=variablesm[f"adds_per_group_{mode}"],
                bol_cmd_bias=variablesm[f"bol_cmd_bias_{channel}_{mode}"]
                / 25.5,  # needs this factor to put it into volts (from pipeline)
                bol_volt=variablesm[f"bol_volt_{channel}_{mode}"],
                fnyq_icm=fnyq["icm"][frec],
                otf=otf,
                Jo=Jo,
                Jg=Jg,
                Tbol=variablesm[f"bolometer_{channel}_{mode}"],
                rho=rho,
                R0=R0,
                T0=T0,
                beta=beta,
                G1=G1,
                C3=C3,
                C1=C1,
                gain=variablesm[f"gain_{channel}_{mode}"],
                sweeps=variablesm[f"sweeps_{mode}"],
                apod=apod,
            )

            print("Subtracting known sources from spectra, including XCAL")

            # frequency mapping
            f_ghz = mu.generate_frequencies(channel, mode)


            bb_ical = mu.planck(
                f_ghz,
                variablesm[f"ical_{mode}"],
            )
            ical_emiss = (
                fits_data[1].data["RICAL"][0]
                + 1j * fits_data[1].data["IICAL"][0]
            )
            ical_emiss = ical_emiss[
                np.abs(ical_emiss) > 0
            ]

            bb_xcal = mu.planck(
                f_ghz,
                variablesm[f"xcal_{mode}"],
            )

            # dihedral spectrum
            bb_dihedral = mu.planck(
                f_ghz,
                variablesm[f"dihedral_{mode}"],
            )
            dihedral_emiss = (
                fits_data[1].data["RDIHEDRA"][0]
                + 1j * fits_data[1].data["IDIHEDRA"][0]
            )
            dihedral_emiss = dihedral_emiss[
                np.abs(dihedral_emiss) > 0
            ]

            bb_refhorn = mu.planck(
                f_ghz,
                variablesm[f"refhorn_{mode}"],
            )
            refhorn_emiss = (
                fits_data[1].data["RREFHORN"][0]
                + 1j * fits_data[1].data["IREFHORN"][0]
            )
            refhorn_emiss = refhorn_emiss[
                np.abs(refhorn_emiss) > 0
            ]

            # skyhorn spectrum
            bb_skyhorn = mu.planck(
                f_ghz,
                variablesm[f"skyhorn_{mode}"],
            )
            skyhorn_emiss = (
                fits_data[1].data["RSKYHORN"][0]
                + 1j * fits_data[1].data["ISKYHORN"][0]
            )
            skyhorn_emiss = skyhorn_emiss[
                np.abs(skyhorn_emiss) > 0
            ]

            # collimating mirror
            bb_collimator = mu.planck(
                f_ghz,
                variablesm[f"collimator_{mode}"],
            )
            collimator_emiss = (
                fits_data[1].data["RSTRUCTU"][0]
                + 1j * fits_data[1].data["ISTRUCTU"][0]
            )
            collimator_emiss = collimator_emiss[
                np.abs(collimator_emiss) > 0
            ]

            # bolometer spectrum
            bb_bolometer_rh = mu.planck(
                f_ghz,
                variablesm[f"bolometer_rh_{mode}"],
            )
            bb_bolometer_rl = mu.planck(
                f_ghz,
                variablesm[f"bolometer_rl_{mode}"],
            )
            bb_bolometer_lh = mu.planck(
                f_ghz,
                variablesm[f"bolometer_lh_{mode}"],
            )
            bb_bolometer_ll = mu.planck(
                f_ghz,
                variablesm[f"bolometer_ll_{mode}"],
            )
            bolometer_emiss = (
                fits_data[1].data["RBOLOMET"][0]
                + 1j * fits_data[1].data["IBOLOMET"][0]
            )
            bolometer_emiss = bolometer_emiss[
                np.abs(bolometer_emiss) > 0
            ]

            if mode[1] == "s":
                cutoff = 5
            else:
                cutoff = 7

            cal[f"{channel}_{mode}"] = (
                spec[f"spec_{channel}_{mode}"][
                    :, cutoff : (len(otf) + cutoff)
                ]
                - (
                    (bb_ical * ical_emiss)
                    + (
                        bb_dihedral
                        * dihedral_emiss
                    )
                    + (
                        bb_refhorn
                        * refhorn_emiss
                    )
                    + (
                        bb_skyhorn
                        * skyhorn_emiss
                    )
                    + (
                        bb_collimator
                        * collimator_emiss
                    )
                    + (
                        bb_bolometer_rh
                        * bolometer_emiss
                    )
                    + (
                        bb_bolometer_rl
                        * bolometer_emiss
                    )
                    + (
                        bb_bolometer_lh
                        * bolometer_emiss
                    )
                    + (
                        bb_bolometer_ll
                        * bolometer_emiss
                    )
                )
                / otf[np.newaxis, :] - bb_xcal
            )

print("Saving cal")

# other varibles that i just want to save in the same place
extra_variables = [
    "stat_word_1",
    "stat_word_12",
    "stat_word_4",
    "stat_word_8",
    "power_a_status_a",
    "power_a_status_b",
    "power_b_status_a",
    "power_b_status_b",
    "ref_hrn_temp_a",
    "ref_hrn_temp_b",
    "ext_cal_temp_a",
    "ext_cal_temp_b",
    "hot_spot_cmd_a",
    "hot_spot_cmd_b",
    "int_ref_temp_a",
    "int_ref_temp_b",
    "scan",
]

# save the calibration data
np.savez(
    g.PROCESSED_DATA_PATH_CAL,
    **variablesm,
    **cal,
    **spec
)
