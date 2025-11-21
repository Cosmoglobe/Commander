import argparse

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import globals as g
import utils.my_utils as utils
from calibration import bolometer
from pipeline import ifg_spec
from simulations.main import generate_ifg
from utils.config import gen_nyquistl

# argparse add-on to filenames and titles
parser = argparse.ArgumentParser(description="Process some calibration data.")
parser.add_argument(
    "--suffix", type=str, default="", help="Suffix to add to filenames and titles"
)
args = parser.parse_args()

if args.suffix != "":
    suffix = f"_{args.suffix}"
else:
    suffix = ""

suffix = "_otf_3rddeg_ical_2nddeg"

fnyq = gen_nyquistl(
    "../reference/fex_samprate.txt",
    "../reference/fex_nyquist.txt",
    "int",
)

# set all text in figures bigger
plt.rcParams.update(
    {
        "font.size": 16,
        "axes.titlesize": 18,
        "axes.labelsize": 16,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
        "figure.titlesize": 20,
    }
)

print("Loading calibration data...")
for channel in g.CHANNELS:
    data = np.load(f"{g.PREPROCESSED_DATA_PATH}cal_{channel}.npz")

    mtm_length = data["mtm_length"][:]
    mtm_speed = data["mtm_speed"][:]

    for mode in g.MODES:
        if channel[1] == "h" and mode == "lf":
            continue
        print(f"Simulating IFGs for {channel.upper()} {mode.upper()}...")

        fits_data = fits.open(
            f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
        )
        # apod = fits_data[1].data["APODIZAT"][0]
        apod = np.ones(512, dtype=np.float64)  # No apodization for now

        if mode[0] == "s":
            length_filter = mtm_length == 0
        else:
            length_filter = mtm_length == 1
        if mode[1] == "s":
            speed_filter = mtm_speed == 0
        else:
            speed_filter = mtm_speed == 1

        mode_filter = length_filter & speed_filter

        xcal = data["xcal_cone"][mode_filter]  # TODO: update this
        ical = data["ical"][mode_filter]  # - 1.5e-1  # applying -150 mK offset
        dihedral = data["dihedral"][mode_filter]
        refhorn = data["refhorn"][mode_filter]
        skyhorn = data["skyhorn"][mode_filter]
        collimator = data["collimator"][mode_filter]
        bolometer = data["bolometer"][mode_filter]
        temps = np.vstack(
            [xcal, ical, dihedral, refhorn, skyhorn, collimator, bolometer]
        )

        adds_per_group = data["adds_per_group"][mode_filter]
        sweeps = data["sweeps"][mode_filter]
        bol_cmd_bias = data["bol_cmd_bias"][mode_filter]
        bol_volt = data["bol_volt"][mode_filter]
        gain = data["gain"][mode_filter]

        original_ifgs = data[f"ifg"][mode_filter]

        # stuff to process the original IFGs into spectra
        frec = 4 * (g.CHANNELS[channel] % 2) + g.MODES[mode]

        otf = fits_data[1].data["RTRANSFE"][0] + 1j * fits_data[1].data["ITRANSFE"][0]
        # print(f"shape of OTF: {otf.shape}")
        # otf = otf[np.abs(otf) > 0]

        fnyq_icm = fnyq["icm"][frec]

        processed_spectra = ifg_spec.ifg_to_spec(
            original_ifgs,
            channel=channel,
            mode=mode,
            adds_per_group=adds_per_group,
            sweeps=sweeps,
            bol_cmd_bias=bol_cmd_bias / 25.5,
            bol_volt=bol_volt,
            gain=gain,
            fnyq_icm=fnyq_icm,
            otf=otf,
            Tbol=temps[6],
            apod=apod,
        )

        cutoff = 5 if mode[1] == "s" else 7

        f_ghz = utils.generate_frequencies(channel, mode, 257)

        bb_xcal = utils.planck(f_ghz, xcal)
        bb_ical = utils.planck(f_ghz, ical)
        bb_dihedral = utils.planck(f_ghz, dihedral)
        bb_refhorn = utils.planck(f_ghz, refhorn)
        bb_skyhorn = utils.planck(f_ghz, skyhorn)
        bb_collimator = utils.planck(f_ghz, collimator)
        bb_bolometer = utils.planck(f_ghz, bolometer)
        length = len(otf)
        (
            emiss_ical,
            emiss_dihedral,
            emiss_refhorn,
            emiss_skyhorn,
            emiss_bolometer,
            emiss_collimator,
        ) = np.zeros((6, 257), dtype=np.complex128)

        otf_save = otf
        otf = np.zeros(257, dtype=np.complex128)
        otf[cutoff : cutoff + length] = otf_save

        # fit 1st degree polynomial to otf
        otf_fit = np.polyfit(f_ghz, otf, 3)
        otf = np.polyval(otf_fit, f_ghz)
        # otf = np.mean(otf[otf != 0])

        # # otf = np.ones(257, dtype=np.complex128) * 0.01

        emiss_ical[cutoff : cutoff + length] = (
            fits_data[1].data["RICAL"][0] + 1j * fits_data[1].data["IICAL"][0]
        )
        # fit 1st degree polynomial to emiss_ical
        emiss_ical_fit = np.polyfit(f_ghz, emiss_ical, 2)
        emiss_ical = np.polyval(emiss_ical_fit, f_ghz)
        # # emiss_ical = np.mean(emiss_ical[emiss_ical != 0])
        # # emiss_ical = np.ones(257, dtype=np.complex128) * -0.0096
        emiss_dihedral[cutoff : cutoff + length] = (
            fits_data[1].data["RDIHEDRA"][0] + 1j * fits_data[1].data["IDIHEDRA"][0]
        )
        emiss_dihedral = np.mean(emiss_dihedral[emiss_dihedral != 0]) * np.zeros(
            257, dtype=np.complex128
        )
        emiss_refhorn[cutoff : cutoff + length] = (
            fits_data[1].data["RREFHORN"][0] + 1j * fits_data[1].data["IREFHORN"][0]
        )
        emiss_refhorn = np.mean(emiss_refhorn[emiss_refhorn != 0]) * np.zeros(
            257, dtype=np.complex128
        )
        emiss_skyhorn[cutoff : cutoff + length] = (
            fits_data[1].data["RSKYHORN"][0] + 1j * fits_data[1].data["ISKYHORN"][0]
        )
        emiss_skyhorn = np.mean(emiss_skyhorn[emiss_skyhorn != 0]) * np.zeros(
            257, dtype=np.complex128
        )
        emiss_collimator[cutoff : cutoff + length] = (
            fits_data[1].data["RSTRUCTU"][0] + 1j * fits_data[1].data["ISTRUCTU"][0]
        )
        emiss_collimator = np.mean(emiss_collimator[emiss_collimator != 0]) * np.zeros(
            257, dtype=np.complex128
        )
        emiss_bolometer[cutoff : cutoff + length] = (
            fits_data[1].data["RBOLOMET"][0] + 1j * fits_data[1].data["IBOLOMET"][0]
        )
        emiss_bolometer = np.mean(emiss_bolometer[emiss_bolometer != 0]) * np.zeros(
            257, dtype=np.complex128
        )

        # open our generated emissivities
        # data = np.load(f"calibration/output/fitted_emissivities_{channel}_{mode}.npy")

        print(data)

        # print(f"OTF mean: {otf:.6f}")
        # print(f"Emiss ICAL mean: {emiss_ical:.6f}")
        # print(f"Emiss Dihedral mean: {emiss_dihedral:.6f}")
        # print(f"Emiss Refhorn mean: {emiss_refhorn:.6f}")
        # print(f"Emiss Skyhorn mean: {emiss_skyhorn:.6f}")
        # print(f"Emiss Collimator mean: {emiss_collimator:.6f}")
        # print(f"Emiss Bolometer mean: {emiss_bolometer:.6f}")

        # plot emissivities
        fig, ax = plt.subplots(4, 2, figsize=(15, 20), sharex=True, sharey=True)
        fig.suptitle(f"{channel.upper()} {mode.upper()} Emissivities")
        ax.flatten()[1].plot(f_ghz, otf.real, label="Real")
        ax.flatten()[1].plot(f_ghz, otf.imag, label="Imag")
        ax.flatten()[1].set_title("OTF")
        ax.flatten()[1].legend()

        ax.flatten()[2].plot(f_ghz, emiss_ical.real, label="Real")
        ax.flatten()[2].plot(f_ghz, emiss_ical.imag, label="Imag")
        ax.flatten()[2].set_title("ICAL Emissivity")
        ax.flatten()[2].legend()

        ax.flatten()[3].plot(f_ghz, emiss_dihedral.real, label="Real")
        ax.flatten()[3].plot(f_ghz, emiss_dihedral.imag, label="Imag")
        ax.flatten()[3].set_title("Dihedral Emissivity")
        ax.flatten()[3].legend()

        ax.flatten()[4].plot(f_ghz, emiss_refhorn.real, label="Real")
        ax.flatten()[4].plot(f_ghz, emiss_refhorn.imag, label="Imag")
        ax.flatten()[4].set_title("Refhorn Emissivity")
        ax.flatten()[4].legend()

        ax.flatten()[5].plot(f_ghz, emiss_skyhorn.real, label="Real")
        ax.flatten()[5].plot(f_ghz, emiss_skyhorn.imag, label="Imag")
        ax.flatten()[5].set_title("Skyhorn Emissivity")
        ax.flatten()[5].legend()

        ax.flatten()[6].plot(f_ghz, emiss_collimator.real, label="Real")
        ax.flatten()[6].plot(f_ghz, emiss_collimator.imag, label="Imag")
        ax.flatten()[6].set_title("Collimator Emissivity")
        ax.flatten()[6].legend()

        ax.flatten()[7].plot(f_ghz, emiss_bolometer.real, label="Real")
        ax.flatten()[7].plot(f_ghz, emiss_bolometer.imag, label="Imag")
        ax.flatten()[7].set_title("Bolometer Emissivity")
        ax.flatten()[7].legend()
        plt.savefig(f"calibration/output/checks/{channel}_{mode}/emissivities.png")

        # same but all divided by otf
        fig, ax = plt.subplots(3, 2, figsize=(15, 20), sharex=True, sharey=True)
        fig.suptitle(f"{channel.upper()} {mode.upper()} Emissivities / OTF")

        ax.flatten()[0].plot(
            f_ghz,
            emiss_ical.real / otf.real,
            label="Real",
        )
        ax.flatten()[0].plot(
            f_ghz,
            emiss_ical.imag / otf.imag,
            label="Imag",
        )
        ax.flatten()[0].set_title("ICAL")
        ax.flatten()[0].legend()

        ax.flatten()[1].plot(
            f_ghz,
            emiss_dihedral.real / otf.real,
            label="Real",
        )
        ax.flatten()[1].plot(
            f_ghz,
            emiss_dihedral.imag / otf.imag,
            label="Imag",
        )
        ax.flatten()[1].set_title("Dihedral")
        ax.flatten()[1].legend()

        ax.flatten()[2].plot(
            f_ghz,
            emiss_refhorn.real / otf.real,
            label="Real",
        )
        ax.flatten()[2].plot(
            f_ghz,
            emiss_refhorn.imag / otf.imag,
            label="Imag",
        )
        ax.flatten()[2].set_title("Refhorn")
        ax.flatten()[2].legend()

        ax.flatten()[3].plot(
            f_ghz,
            emiss_skyhorn.real / otf.real,
            label="Real",
        )
        ax.flatten()[3].plot(
            f_ghz,
            emiss_skyhorn.imag / otf.imag,
            label="Imag",
        )
        ax.flatten()[3].set_title("Skyhorn")
        ax.flatten()[3].legend()

        ax.flatten()[4].plot(
            f_ghz,
            emiss_collimator.real / otf.real,
            label="Real",
        )
        ax.flatten()[4].plot(
            f_ghz,
            emiss_collimator.imag / otf.imag,
            label="Imag",
        )
        ax.flatten()[4].set_title("Collimator")
        ax.flatten()[4].legend()

        ax.flatten()[5].plot(
            f_ghz,
            emiss_bolometer.real / otf.real,
            label="Real",
        )
        ax.flatten()[5].plot(
            f_ghz,
            emiss_bolometer.imag / otf.imag,
            label="Imag",
        )
        ax.flatten()[5].set_title("Bolometer")
        ax.flatten()[5].legend()
        plt.savefig(
            f"calibration/output/checks/{channel}_{mode}/emissivities_div_otf.png"
        )

        xcal_spectra_processed = (
            processed_spectra
            - (
                bb_ical * emiss_ical
                + bb_dihedral * emiss_dihedral
                + bb_refhorn * emiss_refhorn
                + bb_skyhorn * emiss_skyhorn
                + bb_collimator * emiss_collimator
                + bb_bolometer * emiss_bolometer
            )
            / otf
        )

        ical_spectra_processed = (processed_spectra - bb_xcal) * otf / emiss_ical - (
            bb_dihedral * emiss_dihedral
            + bb_refhorn * emiss_refhorn
            + bb_skyhorn * emiss_skyhorn
            + bb_collimator * emiss_collimator
            + bb_bolometer * emiss_bolometer
        ) / emiss_ical

        dihedral_spectra_processed = (
            processed_spectra - bb_xcal
        ) * otf / emiss_dihedral - (
            bb_ical * emiss_ical
            + bb_refhorn * emiss_refhorn
            + bb_skyhorn * emiss_skyhorn
            + bb_collimator * emiss_collimator
            + bb_bolometer * emiss_bolometer
        ) / emiss_dihedral

        refhorn_spectra_processed = (
            processed_spectra - bb_xcal
        ) * otf / emiss_refhorn - (
            bb_ical * emiss_ical
            + bb_dihedral * emiss_dihedral
            + bb_skyhorn * emiss_skyhorn
            + bb_collimator * emiss_collimator
            + bb_bolometer * emiss_bolometer
        ) / emiss_refhorn

        skyhorn_spectra_processed = (
            processed_spectra - bb_xcal
        ) * otf / emiss_skyhorn - (
            bb_ical * emiss_ical
            + bb_dihedral * emiss_dihedral
            + bb_refhorn * emiss_refhorn
            + bb_collimator * emiss_collimator
            + bb_bolometer * emiss_bolometer
        ) / emiss_skyhorn

        collimator_spectra_processed = (
            processed_spectra - bb_xcal
        ) * otf / emiss_collimator - (
            bb_ical * emiss_ical
            + bb_dihedral * emiss_dihedral
            + bb_refhorn * emiss_refhorn
            + bb_skyhorn * emiss_skyhorn
            + bb_bolometer * emiss_bolometer
        ) / emiss_collimator

        bolometer_spectra_processed = (
            processed_spectra - bb_xcal
        ) * otf / emiss_bolometer - (
            bb_ical * emiss_ical
            + bb_dihedral * emiss_dihedral
            + bb_refhorn * emiss_refhorn
            + bb_skyhorn * emiss_skyhorn
            + bb_collimator * emiss_collimator
        ) / emiss_bolometer

        simulated_ifgs, simulated_spectra = generate_ifg(
            channel=channel,
            mode=mode,
            temps=temps,
            apod=apod,
            adds_per_group=adds_per_group,
            sweeps=sweeps,
            bol_cmd_bias=bol_cmd_bias,
            bol_volt=bol_volt,
            gain=gain,
            emiss_xcal=otf,
            emiss_ical=emiss_ical,
            emiss_dihedral=emiss_dihedral,
            emiss_refhorn=emiss_refhorn,
            emiss_skyhorn=emiss_skyhorn,
            emiss_collimator=emiss_collimator,
            emiss_bolometer=emiss_bolometer,
        )

        # n = np.random.randint(0, simulated_spectra.shape[0])
        # n = 1663
        # n = 1704
        n = 2091

        fig, ax = plt.subplots(4, 2, figsize=(15, 20), sharex=True, sharey=False)
        fig.suptitle(
            f"{channel.upper()} {mode.upper()} Black Body Spectra {n}\nTemps: XCAL={xcal[n]:.2f}, ICAL={ical[n]:.2f}, dihed={dihedral[n]:.2f}, refhorn={refhorn[n]:.2f}, skyhorn={skyhorn[n]:.2f}, collimator={collimator[n]:.2f}, bolometer={bolometer[n]:.2f}"
        )

        # Add secondary x-axis to top subplot showing wavenumber
        ghz_to_icm = g.C / 1e7  # GHz to cm/s
        secax = ax.flatten()[0].secondary_xaxis(
            "top",
            functions=(lambda x: x / ghz_to_icm, lambda x: x * ghz_to_icm),
        )
        secax.set_xlabel("Wavenumber (cm⁻¹)")
        secax = ax.flatten()[1].secondary_xaxis(
            "top", functions=(lambda x: x / ghz_to_icm, lambda x: x * ghz_to_icm)
        )
        secax.set_xlabel("Wavenumber (cm⁻¹)")

        ax.flatten()[0].plot(
            f_ghz,
            # np.abs(processed_spectra[n]),
            processed_spectra[n],
            label="Original Processed",
        )
        # ax.flatten()[0].plot(f_ghz, np.abs(simulated_spectra[n]), label="Sum of BBs")
        ax.flatten()[0].plot(
            f_ghz,
            simulated_spectra[n],
            label="Sum of BBs",
        )
        ax.flatten()[0].set_title(f"Processed Spectra")
        ax.flatten()[0].set_ylabel("MJy/sr")
        ax.flatten()[0].legend()

        ax.flatten()[1].plot(
            f_ghz,
            # np.abs(xcal_spectra_processed[n]),
            xcal_spectra_processed[n],
            label="Original Processed",
        )
        ax.flatten()[1].plot(
            f_ghz,
            bb_xcal[n],
            label="Black Body",
        )
        ax.flatten()[1].set_title(f"XCAL Spectra")
        ax.flatten()[1].legend()

        ax.flatten()[2].plot(
            f_ghz,
            # np.abs(ical_spectra_processed[n]),
            ical_spectra_processed[n],
            label="Original Processed",
        )
        ax.flatten()[2].plot(f_ghz, bb_ical[n], label="Black Body")
        ax.flatten()[2].set_title(f"ICAL Spectra")
        ax.flatten()[2].set_ylabel("MJy/sr")
        ax.flatten()[2].legend()

        ax.flatten()[3].plot(
            # f_ghz, np.abs(dihedral_spectra_processed[n]), label="Original Processed"
            f_ghz,
            dihedral_spectra_processed[n],
            label="Original Processed",
        )
        ax.flatten()[3].plot(f_ghz, bb_dihedral[n], label="Black Body")
        ax.flatten()[3].set_title(f"Dihedral Spectra")
        ax.flatten()[3].legend()

        ax.flatten()[4].plot(
            # f_ghz, np.abs(refhorn_spectra_processed[n]), label="Original Processed"
            f_ghz,
            refhorn_spectra_processed[n],
            label="Original Processed",
        )
        ax.flatten()[4].plot(f_ghz, bb_refhorn[n], label="Black Body")
        ax.flatten()[4].set_title(f"Refhorn Spectra")
        ax.flatten()[4].set_ylabel("MJy/sr")
        ax.flatten()[4].legend()

        ax.flatten()[5].plot(
            # f_ghz, np.abs(skyhorn_spectra_processed[n]), label="Original Processed"
            f_ghz,
            skyhorn_spectra_processed[n],
            label="Original Processed",
        )
        ax.flatten()[5].plot(f_ghz, bb_skyhorn[n], label="Black Body")
        ax.flatten()[5].set_title(f"Skyhorn Spectra")
        ax.flatten()[5].legend()

        ax.flatten()[6].plot(
            # f_ghz, np.abs(collimator_spectra_processed[n]), label="Original Processed"
            f_ghz,
            collimator_spectra_processed[n],
            label="Original Processed",
        )
        ax.flatten()[6].plot(f_ghz, bb_collimator[n], label="Black Body")
        ax.flatten()[6].set_title(f"Collimator Spectra")
        ax.flatten()[6].set_xlabel("Frequency (GHz)")
        ax.flatten()[6].set_ylabel("MJy/sr")
        ax.flatten()[6].legend()

        ax.flatten()[7].plot(
            # f_ghz, np.abs(bolometer_spectra_processed[n]), label="Original Processed"
            f_ghz,
            bolometer_spectra_processed[n],
            label="Original Processed",
        )
        ax.flatten()[7].plot(f_ghz, bb_bolometer[n], label="Black Body")
        ax.flatten()[7].set_title(f"Bolometer Spectra")
        ax.flatten()[7].set_xlabel("Frequency (GHz)")
        ax.flatten()[7].legend()

        fig.tight_layout()

        fig.savefig(
            f"calibration/output/checks/{channel}_{mode}/{n}_01_bb_spectra{suffix}.png",
            bbox_inches="tight",
        )
        plt.close()

        # same plot but with bb * emissivity / otf
        fig, ax = plt.subplots(4, 2, figsize=(15, 20), sharex=True, sharey=False)
        fig.suptitle(
            f"{channel.upper()} {mode.upper()} Black Body Spectra with Emissivities {n}\nTemps: XCAL={xcal[n]:.2f}, ICAL={ical[n]:.2f}, dihed={dihedral[n]:.2f}, refhorn={refhorn[n]:.2f}, skyhorn={skyhorn[n]:.2f}, collimator={collimator[n]:.2f}, bolometer={bolometer[n]:.2f}"
        )

        secax = ax.flatten()[0].secondary_xaxis(
            "top",
            functions=(lambda x: x / ghz_to_icm, lambda x: x * ghz_to_icm),
        )
        secax.set_xlabel("Wavenumber (cm⁻¹)")
        secax = ax.flatten()[1].secondary_xaxis(
            "top", functions=(lambda x: x / ghz_to_icm, lambda x: x * ghz_to_icm)
        )
        secax.set_xlabel("Wavenumber (cm⁻¹)")

        ax.flatten()[0].plot(
            f_ghz,
            processed_spectra[n],
            label="Original Processed",
        )
        ax.flatten()[0].plot(
            f_ghz,
            simulated_spectra[n],
            label="Sum of BBs",
        )
        ax.flatten()[0].set_title(f"Processed Spectra")
        ax.flatten()[0].set_ylabel("MJy/sr")
        ax.flatten()[0].legend()

        ax.flatten()[1].plot(
            f_ghz,
            xcal_spectra_processed[n],
            label="Original Processed",
        )
        ax.flatten()[1].plot(
            f_ghz,
            bb_xcal[n],
            label="BB with Emissivity",
        )
        ax.flatten()[1].set_title(f"XCAL Spectra")
        ax.flatten()[1].legend()

        ax.flatten()[2].plot(
            f_ghz,
            ical_spectra_processed[n] * emiss_ical / otf,
            label="Original Processed",
        )
        ax.flatten()[2].plot(
            f_ghz,
            bb_ical[n] * emiss_ical / otf,
            label="BB with Emissivity",
        )
        ax.flatten()[2].set_title(f"ICAL Spectra")
        ax.flatten()[2].set_ylabel("MJy/sr")
        ax.flatten()[2].legend()

        ax.flatten()[3].plot(
            f_ghz,
            dihedral_spectra_processed[n] * emiss_dihedral / otf,
            label="Original Processed",
        )
        ax.flatten()[3].plot(
            f_ghz,
            bb_dihedral[n] * emiss_dihedral / otf,
            label="BB with Emissivity",
        )
        ax.flatten()[3].set_title(f"Dihedral Spectra")
        ax.flatten()[3].legend()

        ax.flatten()[4].plot(
            f_ghz,
            refhorn_spectra_processed[n] * emiss_refhorn / otf,
            label="Original Processed",
        )
        ax.flatten()[4].plot(
            f_ghz,
            bb_refhorn[n] * emiss_refhorn / otf,
            label="BB with Emissivity",
        )
        ax.flatten()[4].set_title(f"Refhorn Spectra")
        ax.flatten()[4].set_ylabel("MJy/sr")
        ax.flatten()[4].legend()

        ax.flatten()[5].plot(
            f_ghz,
            skyhorn_spectra_processed[n] * emiss_skyhorn / otf,
            label="Original Processed",
        )
        ax.flatten()[5].plot(
            f_ghz,
            bb_skyhorn[n] * emiss_skyhorn / otf,
            label="BB with Emissivity",
        )
        ax.flatten()[5].set_title(f"Skyhorn Spectra")
        ax.flatten()[5].legend()

        ax.flatten()[6].plot(
            f_ghz,
            collimator_spectra_processed[n] * emiss_collimator / otf,
            label="Original Processed",
        )
        ax.flatten()[6].plot(
            f_ghz,
            bb_collimator[n] * emiss_collimator / otf,
            label="BB with Emissivity",
        )
        ax.flatten()[6].set_title(f"Collimator Spectra")
        ax.flatten()[6].set_xlabel("Frequency (GHz)")
        ax.flatten()[6].set_ylabel("MJy/sr")
        ax.flatten()[6].legend()

        ax.flatten()[7].plot(
            f_ghz,
            bolometer_spectra_processed[n] * emiss_bolometer / otf,
            label="Original Processed",
        )
        ax.flatten()[7].plot(
            f_ghz,
            bb_bolometer[n] * emiss_bolometer / otf,
            label="BB with Emissivity",
        )
        ax.flatten()[7].set_title(f"Bolometer Spectra")
        ax.flatten()[7].set_xlabel("Frequency (GHz)")
        ax.flatten()[7].legend()
        fig.tight_layout()
        plt.savefig(
            f"calibration/output/checks/{channel}_{mode}/{n}_01b_bb_spectra_emissivities{suffix}.png",
            bbox_inches="tight",
        )
        plt.close()

        # plot residuals
        fig, ax = plt.subplots(4, 2, figsize=(15, 20), sharex=True, sharey=False)
        fig.suptitle(
            f"{channel.upper()} {mode.upper()} Residuals {n}\nTemps: XCAL={xcal[n]:.2f}, ICAL={ical[n]:.2f}, dihed={dihedral[n]:.2f}, refhorn={refhorn[n]:.2f}, skyhorn={skyhorn[n]:.2f}, collimator={collimator[n]:.2f}, bolometer={bolometer[n]:.2f}"
        )

        secax = ax.flatten()[0].secondary_xaxis(
            "top",
            functions=(lambda x: x / ghz_to_icm, lambda x: x * ghz_to_icm),
        )
        secax.set_xlabel("Wavenumber (cm⁻¹)")
        secax = ax.flatten()[1].secondary_xaxis(
            "top", functions=(lambda x: x / ghz_to_icm, lambda x: x * ghz_to_icm)
        )
        secax.set_xlabel("Wavenumber (cm⁻¹)")

        # ax.flatten()[0].plot(f_ghz, np.abs(processed_spectra[n] - simulated_spectra[n]))
        ax.flatten()[0].plot(f_ghz, processed_spectra[n] - simulated_spectra[n])
        ax.flatten()[0].set_title(f"Processed Spectra")
        ax.flatten()[0].set_ylabel("MJy/sr")

        # ax.flatten()[1].plot(f_ghz, np.abs(xcal_spectra_processed[n] - bb_xcal[n]))
        ax.flatten()[1].plot(f_ghz, xcal_spectra_processed[n] - bb_xcal[n])
        ax.flatten()[1].set_title(f"XCAL Spectra")

        # ax.flatten()[2].plot(f_ghz, np.abs(ical_spectra_processed[n] - bb_ical[n]))
        ax.flatten()[2].plot(f_ghz, ical_spectra_processed[n] - bb_ical[n])
        ax.flatten()[2].set_title(f"ICAL Spectra")
        ax.flatten()[2].set_ylabel("MJy/sr")

        ax.flatten()[3].plot(
            # f_ghz, np.abs(dihedral_spectra_processed[n] - bb_dihedral[n])
            f_ghz,
            dihedral_spectra_processed[n] - bb_dihedral[n],
        )
        ax.flatten()[3].set_title(f"Dihedral Spectra")

        ax.flatten()[4].plot(
            # f_ghz, np.abs(refhorn_spectra_processed[n] - bb_refhorn[n])
            f_ghz,
            refhorn_spectra_processed[n] - bb_refhorn[n],
        )
        ax.flatten()[4].set_title(f"Refhorn Spectra")
        ax.flatten()[4].set_ylabel("MJy/sr")

        ax.flatten()[5].plot(
            # f_ghz, np.abs(skyhorn_spectra_processed[n] - bb_skyhorn[n])
            f_ghz,
            skyhorn_spectra_processed[n] - bb_skyhorn[n],
        )
        ax.flatten()[5].set_title(f"Skyhorn Spectra")

        ax.flatten()[6].plot(
            # f_ghz, np.abs(collimator_spectra_processed[n] - bb_collimator[n])
            f_ghz,
            collimator_spectra_processed[n] - bb_collimator[n],
        )
        ax.flatten()[6].set_title(f"Collimator Spectra")
        ax.flatten()[6].set_xlabel("Frequency (GHz)")
        ax.flatten()[6].set_ylabel("MJy/sr")

        ax.flatten()[7].plot(
            # f_ghz, np.abs(bolometer_spectra_processed[n] - bb_bolometer[n])
            f_ghz,
            bolometer_spectra_processed[n] - bb_bolometer[n],
        )
        ax.flatten()[7].set_title(f"Bolometer Spectra")
        ax.flatten()[7].set_xlabel("Frequency (GHz)")

        fig.tight_layout()

        plt.savefig(
            f"calibration/output/checks/{channel}_{mode}/{n}_02_residuals{suffix}.png",
            bbox_inches="tight",
        )
        plt.close()

        xcal_residuals = xcal_spectra_processed[n] - bb_xcal[n]
        spectra_corrected = processed_spectra[n] - xcal_residuals

        simulated_plus_noise = simulated_spectra[n] + xcal_residuals

        fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        fig.suptitle(
            f"{channel.upper()} {mode.upper()} Spectra {n}\nTemps: XCAL={xcal[n]:.2f}, ICAL={ical[n]:.2f}, dihed={dihedral[n]:.2f}, refhorn={refhorn[n]:.2f}, skyhorn={skyhorn[n]:.2f}, collimator={collimator[n]:.2f}, bolometer={bolometer[n]:.2f}"
        )

        secax = ax.flatten()[0].secondary_xaxis(
            "top",
            functions=(lambda x: x / ghz_to_icm, lambda x: x * ghz_to_icm),
        )
        secax.set_xlabel("Wavenumber (cm⁻¹)")

        ax[0].plot(
            f_ghz,
            # np.abs(spectra_corrected),
            spectra_corrected,
            label="Processed spectrum corrected for XCAL noise",
        )
        ax[0].plot(
            f_ghz,
            # np.abs(simulated_spectra[n]),
            simulated_spectra[n],
            label="Simulated Spectra",
            linestyle="--",
        )
        ax[0].set_ylabel("MJy/sr")
        ax[0].legend()

        ax[1].plot(
            f_ghz,
            # np.abs(simulated_plus_noise),
            simulated_plus_noise,
            label="Simulated spectrum plus XCAL noise",
        )
        ax[1].plot(
            f_ghz,
            # np.abs(processed_spectra[n]),
            processed_spectra[n],
            label="Original Processed Spectra",
            linestyle="--",
        )
        ax[1].set_xlabel("Frequency (GHz)")
        ax[1].set_ylabel("MJy/sr")
        ax[1].legend()
        fig.tight_layout()

        plt.savefig(
            f"calibration/output/checks/{channel}_{mode}/{n}_03_spectra_comparison{suffix}.png",
            bbox_inches="tight",
        )
        plt.close()

        simulated_ifgs_corrected, _ = generate_ifg(
            channel=channel,
            mode=mode,
            temps=temps,
            apod=apod,
            adds_per_group=adds_per_group,
            sweeps=sweeps,
            bol_cmd_bias=bol_cmd_bias,
            bol_volt=bol_volt,
            gain=gain,
            total_spectra=np.nan_to_num(simulated_plus_noise, nan=0.0),
        )

        if channel[0] == "r":
            simulated_ifgs = -simulated_ifgs
            # simulated_ifgs_corrected = -simulated_ifgs_corrected

        for ifg in range(n, n + 1, 1):
            print(f"Plotting IFG {ifg+1} for {channel.upper()} {mode.upper()}...")

            # plot both and residuals
            fig, ax = plt.subplots(4, 1, figsize=(10, 15), sharex=True, sharey=True)
            fig.suptitle(
                f"{channel.upper()} {mode.upper()} IFG {ifg}\nTemps: XCAL={xcal[ifg]:.2f}, ICAL={ical[ifg]:.2f}, dihed={dihedral[ifg]:.2f}, refhorn={refhorn[ifg]:.2f}, skyhorn={skyhorn[ifg]:.2f}, collimator={collimator[ifg]:.2f}, bolometer={bolometer[ifg]:.2f}"
            )

            if mode == "lf":
                scan_length = 7.07  # cm
            else:
                scan_length = 1.76  # cm

            secax = ax[0].secondary_xaxis(
                "top",
                functions=(
                    lambda x: x * scan_length / 512,
                    lambda x: x * 512 / scan_length,
                ),
            )
            secax.set_xlabel("Length (cm)")

            ax[0].plot(original_ifgs[ifg, :] - np.median(original_ifgs[ifg, :]))
            ax[0].axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            ax[0].set_title("Original IFG")

            ax[1].plot(simulated_ifgs[ifg, :])
            ax[1].axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            ax[1].set_title("Simulated IFG")

            ax[2].plot(
                original_ifgs[ifg, :] - np.median(original_ifgs[ifg, :]),
                label="Original IFG",
            )
            ax[2].plot(simulated_ifgs[ifg, :], label="Simulated IFG")
            ax[2].axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            ax[2].set_title("Original vs Simulated IFG")
            ax[2].legend()

            ax[3].plot(
                original_ifgs[ifg, :]
                - np.median(original_ifgs[ifg, :])
                - simulated_ifgs[ifg, :]
            )
            ax[3].plot(
                [
                    np.median(
                        original_ifgs[ifg, :]
                        - np.median(original_ifgs[ifg, :])
                        - simulated_ifgs[ifg, :]
                    )
                ]
                * 512
            )
            ax[3].axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            ax[3].set_title("Residuals (original - simulated)")

            fig.tight_layout()

            plt.savefig(
                f"calibration/output/checks/{channel}_{mode}/{ifg}_04_ifg_residuals{suffix}.png",
                bbox_inches="tight",
            )
            plt.close()

            fig, ax = plt.subplots(4, 1, figsize=(10, 15), sharex=True, sharey=True)
            fig.suptitle(
                f"{channel.upper()} {mode.upper()} Corrected IFG {ifg}\nTemps: XCAL={xcal[ifg]:.2f}, ICAL={ical[ifg]:.2f}, dihed={dihedral[ifg]:.2f}, refhorn={refhorn[ifg]:.2f}, skyhorn={skyhorn[ifg]:.2f}, collimator={collimator[ifg]:.2f}, bolometer={bolometer[ifg]:.2f}"
            )

            secax = ax.flatten()[0].secondary_xaxis(
                "top",
                functions=(
                    lambda x: x * scan_length / 512,
                    lambda x: x * 512 / scan_length,
                ),
            )
            secax.set_xlabel("Length (cm)")

            ax[0].plot(original_ifgs[ifg, :] - np.median(original_ifgs[ifg, :]))
            ax[0].axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            ax[0].set_title("Original IFG")

            ax[1].plot(simulated_ifgs_corrected[ifg, :])
            ax[1].axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            ax[1].set_title("Simulated Corrected IFG")

            ax[2].plot(
                original_ifgs[ifg, :] - np.median(original_ifgs[ifg, :]),
                label="Original IFG",
            )
            ax[2].plot(
                simulated_ifgs_corrected[ifg, :], label="Simulated Corrected IFG"
            )
            ax[2].axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            ax[2].set_title("Original vs Simulated Corrected IFG")
            ax[2].legend()

            ax[3].plot(
                original_ifgs[ifg, :]
                - np.median(original_ifgs[ifg, :])
                - simulated_ifgs_corrected[ifg, :]
            )
            ax[3].plot(
                [
                    np.median(
                        original_ifgs[ifg, :]
                        - np.median(original_ifgs[ifg, :])
                        - simulated_ifgs_corrected[ifg, :]
                    )
                ]
                * 512
            )
            ax[3].axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            ax[3].set_title("Residuals (original - simulated corrected)")

            fig.tight_layout()

            plt.savefig(
                f"calibration/output/checks/{channel}_{mode}/{ifg}_05_ifg_residuals_corrected{suffix}.png",
                bbox_inches="tight",
            )
            plt.close()
