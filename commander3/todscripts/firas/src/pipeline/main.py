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
from pipeline import ifg_spec
from utils.config import gen_nyquistl

# Start timing
script_start_time = time.time()

# filter out bad data (selected "by eye")
# filter_bad = filter.filter_junk(
#     variables["stat_word_1"],
#     variables["stat_word_5"],
#     variables["stat_word_9"],
#     variables["stat_word_12"],
#     variables["stat_word_13"],
#     variables["stat_word_16"],
#     variables["lvdt_stat_a"],
#     variables["lvdt_stat_b"],
#     variables["a_bol_assem_rh"],
#     variables["a_bol_assem_rl"],
#     variables["a_bol_assem_lh"],
#     variables["b_bol_assem_lh"],
#     variables["a_bol_assem_ll"],
#     variables["b_bol_assem_ll"],
#     variables["bol_cmd_bias_lh"],
#     variables["bol_cmd_bias_rh"],
# )

fnyq = gen_nyquistl(
    "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
)

for channel, channel_value in g.CHANNELS.items():
    sky_data = np.load(
        f"{g.PREPROCESSED_DATA_PATH}/sky_{channel}.npz", allow_pickle=True
    )
    for mode, mode_value in g.MODES.items():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            key = f"{channel}_{mode}"

            length_filter = sky_data["mtm_length"] == (0 if mode[0] == "s" else 1)
            speed_filter = sky_data["mtm_speed"] == (0 if mode[1] == "s" else 1)
            mode_filter = length_filter & speed_filter

            mode_data = {}
            for var in sky_data.files:
                mode_data[var] = sky_data[var][mode_filter]

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

            mode_data["spec"] = ifg_spec.ifg_to_spec(
                ifg=mode_data["ifg"],
                channel=channel,
                mode=mode,
                adds_per_group=mode_data["adds_per_group"],
                bol_cmd_bias=mode_data[f"bol_cmd_bias"]
                / 25.5,  # needs this factor to put it into volts (from pipeline)
                bol_volt=mode_data[f"bol_volt"],
                fnyq_icm=fnyq["icm"][frec],
                otf=otf,
                Tbol=mode_data[f"bolometer"],
                gain=mode_data[f"gain"],
                sweeps=mode_data[f"sweeps"],
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
            bb_ical = utils.planck(f_ghz, mode_data["ical"])
            bb_dihedral = utils.planck(f_ghz, mode_data["dihedral"])
            bb_refhorn = utils.planck(f_ghz, mode_data["refhorn"])
            bb_skyhorn = utils.planck(f_ghz, mode_data["skyhorn"])
            bb_collimator = utils.planck(f_ghz, mode_data["collimator"])
            bb_bolometer = utils.planck(f_ghz, mode_data["bolometer"])
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
            mode_data["sky"] = (
                mode_data["spec"][:, cutoff : (len(otf) + cutoff)] - total_emission
            ) / otf[np.newaxis, :]

            # save
            np.savez(
                f"{g.PROCESSED_DATA_PATH}sky_{channel}_{mode}.npz",
                **mode_data,
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
