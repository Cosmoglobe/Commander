"""
This script takes the interferograms into spectra.
"""

import time

import h5py
import numpy as np
from astropy.io import fits

import globals as g
import utils.my_utils as utils
from calibration import bolometer
from flagging import filter
from pipeline import ifg_spec
from utils.config import gen_nyquistl

# Start timing
script_start_time = time.time()

fnyq = gen_nyquistl(
    "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
)

data = h5py.File(
    g.PREPROCESSED_DATA_PATH_SKY,
    "r+",
)["df_data"]

short_filter = data["mtm_length"][:] == 0
long_filter = data["mtm_length"][:] == 1
slow_filter = data["mtm_speed"][:] == 0
fast_filter = data["mtm_speed"][:] == 1

filters = {}
filters["ss"] = short_filter & slow_filter
filters["lf"] = long_filter & fast_filter

data_mode = {}
print("Filtering variables by mode...")
for mode in g.MODES.keys():
    mode_filter = filters[mode]
    for variable in data.keys():
        data_mode[f"{variable}_{mode}"] = data[variable][mode_filter]

for channel, channel_value in g.CHANNELS.items():
    for mode, mode_value in g.MODES.items():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            # filter out bad data (selected "by eye")
            filter_bad = filter.filter_junk(
                data_mode[f"stat_word_1_{mode}"],
                data_mode[f"stat_word_5_{mode}"],
                data_mode[f"stat_word_9_{mode}"],
                data_mode[f"stat_word_12_{mode}"],
                data_mode[f"stat_word_13_{mode}"],
                data_mode[f"stat_word_16_{mode}"],
                data_mode[f"lvdt_stat_a_{mode}"],
                data_mode[f"lvdt_stat_b_{mode}"],
            )

            sky_data = {}
            for var in data_mode.keys():
                if mode in var:
                    sky_data[var] = data_mode[var][filter_bad]

            frec = 4 * (channel_value % 2) + mode_value

            # bolometer parameters
            R0, T0, G1, beta, rho, C1, C3, Jo, Jg = bolometer.get_bolometer_parameters(
                channel, mode
            )

            print(
                f"Converting interferograms to spectra for {channel.upper()}{mode.upper()}"
            )

            fits_data = fits.open(
                f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
            )
            otf = (
                fits_data[1].data["RTRANSFE"][0] + 1j * fits_data[1].data["ITRANSFE"][0]
            )
            otf = otf[np.abs(otf) > 0]
            apod = fits_data[1].data["APODIZAT"][0]

            sky_data["spec"] = ifg_spec.ifg_to_spec(
                ifg=sky_data[f"ifg_{channel}_{mode}"],
                channel=channel,
                mode=mode,
                adds_per_group=sky_data[f"adds_per_group_{mode}"],
                bol_cmd_bias=sky_data[f"bol_cmd_bias_{channel}_{mode}"]
                / 25.5,  # needs this factor to put it into volts (from pipeline)
                bol_volt=sky_data[f"bol_volt_{channel}_{mode}"],
                fnyq_icm=fnyq["icm"][frec],
                otf=otf,
                Tbol=sky_data[f"bolometer_{channel}_{mode}"],
                gain=sky_data[f"gain_{channel}_{mode}"],
                sweeps=sky_data[f"sweeps_{mode}"],
                apod=apod,
            )

            print("Subtracting by the known emission sources")

            # frequency mapping (computed once per channel/mode)
            f_ghz = utils.generate_frequencies(channel, mode)

            emissivities = {
                "ical": fits_data[1].data["RICAL"][0]
                + 1j * fits_data[1].data["IICAL"][0],
                "dihedral": fits_data[1].data["RDIHEDRA"][0]
                + 1j * fits_data[1].data["IDIHEDRA"][0],
                "refhorn": fits_data[1].data["RREFHORN"][0]
                + 1j * fits_data[1].data["IREFHORN"][0],
                "skyhorn": fits_data[1].data["RSKYHORN"][0]
                + 1j * fits_data[1].data["ISKYHORN"][0],
                "collimator": fits_data[1].data["RSTRUCTU"][0]
                + 1j * fits_data[1].data["ISTRUCTU"][0],
                "bolometer": fits_data[1].data["RBOLOMET"][0]
                + 1j * fits_data[1].data["IBOLOMET"][0],
            }
            for key in emissivities:
                emissivities[key] = emissivities[key][np.abs(emissivities[key]) > 0]

            # Get cached emissivities
            ical_emiss = emissivities["ical"]
            dihedral_emiss = emissivities["dihedral"]
            refhorn_emiss = emissivities["refhorn"]
            skyhorn_emiss = emissivities["skyhorn"]
            collimator_emiss = emissivities["collimator"]
            bolometer_emiss = emissivities["bolometer"]

            # Compute all blackbody spectra at once (vectorized)
            bb_ical = utils.planck(f_ghz, sky_data[f"ical_{mode}"])
            bb_dihedral = utils.planck(f_ghz, sky_data[f"dihedral_{mode}"])
            bb_refhorn = utils.planck(f_ghz, sky_data[f"refhorn_{mode}"])
            bb_skyhorn = utils.planck(f_ghz, sky_data[f"skyhorn_{mode}"])
            bb_collimator = utils.planck(f_ghz, sky_data[f"collimator_{mode}"])
            bb_bolometer = utils.planck(f_ghz, sky_data[f"bolometer_{channel}_{mode}"])
            # bb_bolometer_rl = utils.planck(f_ghz, sky_data["bolometer_rl"])
            # bb_bolometer_lh = utils.planck(f_ghz, sky_data["bolometer_lh"])
            # bb_bolometer_ll = utils.planck(f_ghz, sky_data["bolometer_ll"])

            cutoff = 5 if mode[1] == "s" else 7

            print(f"Cutting off the first {cutoff} frequencies")

            term1 = bb_ical * ical_emiss
            term2 = bb_dihedral * dihedral_emiss
            term3 = bb_refhorn * refhorn_emiss
            term4 = bb_skyhorn * skyhorn_emiss
            term5 = bb_collimator * collimator_emiss
            term6 = bb_bolometer * bolometer_emiss

            total_emission = (
                term1
                + term2
                + term3
                + term4
                + term5
                + term6  # TODO: keep just the one that is being used to measure?
                # + (bb_bolometer_rl * bolometer_emiss)
                # + (bb_bolometer_lh * bolometer_emiss)
                # + (bb_bolometer_ll * bolometer_emiss)
            )

            # Single vectorized operation for sky extraction
            sky_data["sky"] = (
                sky_data["spec"][:, cutoff : (len(otf) + cutoff)]
                - total_emission / otf[np.newaxis, :]
            )

            # check for nans
            if np.isnan(sky_data["sky"]).any():
                print(
                    f"Warning: NaN values found in sky data for {channel.upper()}{mode.upper()}"
                )

            # save
            np.savez(
                f"{g.PROCESSED_DATA_PATH}sky_{channel}_{mode}.npz",
                **sky_data,
            )

# Print timing summary
script_end_time = time.time()
total_time = script_end_time - script_start_time
# saving_time = script_end_time - saving_start_time

print("\n" + "=" * 60)
print("PERFORMANCE SUMMARY")
print("=" * 60)
print(f"Total execution time: {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
# print(f"Data saving time: {saving_time:.2f} seconds")
# print(f"Processing time: {total_time - saving_time:.2f} seconds")
print("=" * 60)
