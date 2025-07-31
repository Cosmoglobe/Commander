"""
This script takes the interferograms into spectra.
"""

import os
import sys
from datetime import datetime

import h5py
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
import globals as g
import my_utils as mu
from utils.config import gen_nyquistl

T_CMB = 2.72548  # Fixsen 2009
NSIDE = 32
channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}
modes = {"ss": 0, "lf": 3}  # can change when i have the new cal models

sky_data = h5py.File(
    g.PREPROCESSED_DATA_PATH_SKY,
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
    "pix_gal",
    "pix_terr",
    "pix_ecl",
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
    "gal_lat",
    "gal_lon",
]
channel_dependent = [
    "ifg",
    "bol_cmd_bias",
    "bol_volt",
    "gain",
]

for channel in channels:
    for variable_name in channel_dependent:
        variable_names = variable_names + [f"{variable_name}_{channel}"]

variables = {}
for variable_name in variable_names:
    variables[variable_name] = np.array(sky_data["df_data/" + variable_name][()])

# Get the galactic longitude and latitude into a vector in order to interpolate
gal_vec = mu.tune_pointing(
    gal_lon=variables["gal_lon"],
    gal_lat=variables["gal_lat"],
    gmt=variables["gmt"],
    mtm_length=variables["mtm_length"],
    mtm_speed=variables["mtm_speed"],
    offset=g.OFFSET,
)

# filter out nans that come from gal_vec into all variables
valid_indices = ~np.isnan(gal_vec).any(axis=1)
for variable in variables.keys():
    variables[variable] = variables[variable][valid_indices]
gal_vec = gal_vec[valid_indices]

variables["gal_lon"], variables["gal_lat"] = hp.vec2ang(gal_vec, lonlat=True)
variables["pix_gal"] = hp.ang2pix(
    NSIDE, variables["gal_lon"], variables["gal_lat"], lonlat=True
).astype(int)

# filter out bad data (selected "by eye")
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

variables["pix_gal"] = variables["pix_gal"].astype(int)
variables["pix_terr"] = variables["pix_terr"].astype(int)
variables["pix_ecl"] = variables["pix_ecl"].astype(int)
variables["gmt"] = variables["gmt"].astype(str)
variables["gmt"] = np.array(
    [datetime.strptime(gmt, "%Y-%m-%d %H:%M:%S") for gmt in variables["gmt"]]
)

# constraining to the recs with ical temp at certain bins
variables["ical"] = (
    # variables["a_ical"] * 0.1 + variables["b_ical"] * 0.9
    # )  # using their weights
    (variables["a_ical"] + variables["b_ical"])
    / 2
)
variables["dihedral"] = (variables["a_dihedral"] + variables["b_dihedral"]) / 2
variables["refhorn"] = (variables["a_refhorn"] + variables["b_refhorn"]) / 2
variables["skyhorn"] = (variables["a_skyhorn"] + variables["b_skyhorn"]) / 2
variables["collimator"] = (
    variables["a_collimator"] + variables["b_collimator"]
) / 2
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

del variables["a_ical"], variables["b_ical"], variables["a_dihedral"], variables["b_dihedral"], variables["a_refhorn"], variables["b_refhorn"], variables["a_skyhorn"], variables["b_skyhorn"], variables["a_collimator"], variables["b_collimator"], variables["a_bol_assem_rh"], variables["b_bol_assem_rh"], variables["a_bol_assem_rl"], variables["b_bol_assem_rl"], variables["a_bol_assem_lh"], variables["b_bol_assem_lh"], variables["a_bol_assem_ll"], variables["b_bol_assem_ll"]

# filtering the data based on the mode used
short_filter = variables["mtm_length"] == 0
long_filter = variables["mtm_length"] == 1
slow_filter = variables["mtm_speed"] == 0
fast_filter = variables["mtm_speed"] == 1

filters = {}
filters["ss"] = short_filter & slow_filter
filters["lf"] = long_filter & fast_filter

variablesm = {}
for variable in variables.keys():
    for mode in modes.keys():
        variablesm[f"{variable}_{mode}"] = variables[variable][filters[mode]]

fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)

spec = {}
sky = {}

for channel, channel_value in channels.items():
    for mode, mode_value in modes.items():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            fits_data = fits.open(
                f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
            )

            frec = 4 * (channel_value % 2) + mode_value

            # optical transfer function
            otf = (
                fits_data[1].data["RTRANSFE"][0]
                + 1j * fits_data[1].data["ITRANSFE"][0]
            )
            otf = otf[np.abs(otf) > 0]

            apod = fits_data[1].data[
                "APODIZAT"
            ][0]

            # bolometer parameters
            R0, T0, G1, beta, rho, C1, C3, Jo, Jg = mu.get_bolometer_parameters(channel, mode)

            print(f"Converting interferograms to spectra for {channel.upper()}{mode.upper()}")

            afreq, spec[f"spec_{channel}_{mode}"] = mu.ifg_to_spec(
                ifg=variablesm[f"ifg_{channel}_{mode}"],
                channel=channel,
                mode=mode,
                adds_per_group=variablesm[f"adds_per_group_{mode}"],
                bol_cmd_bias=variablesm[f"bol_cmd_bias_{channel}_{mode}"] / 25.5,  # needs this factor to put it into volts (from pipeline)
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

            print("Subtracting by the known emission sources")

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

            bb_collimator= mu.planck(
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

            print(f"Cutting off the first {cutoff} frequencies")

            # setting the frequency cut-off according to the header of the fits file
            sky[f"{channel}_{mode}"] = (
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
                        bb_bolometer_rh # TODO: keep just the one that is being used to measure?
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
                / otf[np.newaxis, :]
            )

print("Saving sky")

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
    "ext_cal_temp_a",
    "ext_cal_temp_b",
    "grt_addr_a",
    "grt_addr_b",
    "hot_spot_cmd_a",
    "hot_spot_cmd_b",
    "int_ref_temp_a",
    "int_ref_temp_b",
    # "sky_hrn_temp_a",
    # "sky_hrn_temp_b",
    "scan",
]
for variable in extra_variables:
    for mode in modes.keys():
        variablesm[f"{variable}_{mode}"] = np.array(sky_data["df_data/" + variable])[
            valid_indices
        ]
        variablesm[f"{variable}_{mode}"] = variablesm[f"{variable}_{mode}"][filter_bad]
        variablesm[f"{variable}_{mode}"] = variablesm[f"{variable}_{mode}"][
            filters[mode]
        ]

# save the sky
np.savez(
    g.PROCESSED_DATA_PATH,
    **variablesm,
    **sky,
    **spec
)
