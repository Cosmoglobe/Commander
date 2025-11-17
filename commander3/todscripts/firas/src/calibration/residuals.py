import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import globals as g
import utils.my_utils as utils
from calibration import bolometer
from pipeline import ifg_spec
from simulations.main import generate_ifg
from utils.config import gen_nyquistl

fnyq = gen_nyquistl(
    "../reference/fex_samprate.txt",
    "../reference/fex_nyquist.txt",
    "int",
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
        ical = data["ical"][mode_filter]
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
        print(f"shape of OTF: {otf.shape}")
        # otf = otf[np.abs(otf) > 0]

        processed_spectra = ifg_spec.ifg_to_spec(
            original_ifgs,
            channel=channel,
            mode=mode,
            adds_per_group=adds_per_group,
            sweeps=sweeps,
            bol_cmd_bias=bol_cmd_bias / 25.5,
            bol_volt=bol_volt,
            gain=gain,
            fnyq_icm=fnyq["icm"][frec],
            otf=otf,
            Tbol=temps[6],
            apod=apod,
        )

        cutoff = 5 if mode[1] == "s" else 7

        f_ghz = utils.generate_frequencies(channel, mode, 257)

        bb_ical = utils.planck(f_ghz, ical)
        bb_dihedral = utils.planck(f_ghz, dihedral)
        bb_refhorn = utils.planck(f_ghz, refhorn)
        bb_skyhorn = utils.planck(f_ghz, skyhorn)
        bb_collimator = utils.planck(f_ghz, collimator)
        bb_bolometer = utils.planck(f_ghz, bolometer)
        length = len(otf)
        (
            emiss_xcal,
            emiss_ical,
            emiss_dihedral,
            emiss_refhorn,
            emiss_skyhorn,
            emiss_bolometer,
            emiss_collimator,
        ) = np.zeros((7, 257), dtype=np.complex128)

        emiss_ical[cutoff : cutoff + length] = (
            fits_data[1].data["RICAL"][0] + 1j * fits_data[1].data["IICAL"][0]
        )
        emiss_dihedral[cutoff : cutoff + length] = (
            fits_data[1].data["RDIHEDRA"][0] + 1j * fits_data[1].data["IDIHEDRA"][0]
        )
        emiss_refhorn[cutoff : cutoff + length] = (
            fits_data[1].data["RREFHORN"][0] + 1j * fits_data[1].data["IREFHORN"][0]
        )
        emiss_skyhorn[cutoff : cutoff + length] = (
            fits_data[1].data["RSKYHORN"][0] + 1j * fits_data[1].data["ISKYHORN"][0]
        )
        emiss_collimator[cutoff : cutoff + length] = (
            fits_data[1].data["RSTRUCTU"][0] + 1j * fits_data[1].data["ISTRUCTU"][0]
        )
        emiss_bolometer[cutoff : cutoff + length] = (
            fits_data[1].data["RBOLOMET"][0] + 1j * fits_data[1].data["IBOLOMET"][0]
        )

        total_emission = (
            bb_ical * emiss_ical
            + bb_dihedral * emiss_dihedral
            + bb_refhorn * emiss_refhorn
            + bb_skyhorn * emiss_skyhorn
            + bb_collimator * emiss_collimator
            + bb_bolometer * emiss_bolometer
        )
        xcal_spectra_processed = np.zeros_like(bb_ical)
        xcal_spectra_processed[:, cutoff : (len(otf) + cutoff)] = (
            processed_spectra[:, cutoff : (len(otf) + cutoff)]
            - total_emission[:, cutoff : (len(otf) + cutoff)] / otf
        )

        simulated_ifgs, simulated_spectra, xcal_spectra = generate_ifg(
            channel=channel,
            mode=mode,
            temps=temps,
            real_xcal_spec=xcal_spectra_processed,
            apod=apod,
            adds_per_group=adds_per_group,
            sweeps=sweeps,
            bol_cmd_bias=bol_cmd_bias,
            bol_volt=bol_volt,
            gain=gain,
        )

        if channel[0] == "r":
            simulated_ifgs = -simulated_ifgs
            # simulated_ifgs = simulated_ifgs/3

        # n = np.random.randint(0, simulated_spectra.shape[0])
        # n = 1663
        # n = 1704
        n = 2091
        print(f"peak of simulated spectra: {np.max(np.abs(simulated_spectra[n]))}")

        # np.save(
        #     f"calibration/output/ifgdata_{channel}_{mode}_{n}.npy", simulated_ifgs[n]
        # )

        print(f"Simulated IFGs for {channel.upper()} {mode.upper()}")

        # processed_spectra = data[f"spec_{channel}_{mode}"][:]

        plt.plot(
            np.abs(xcal_spectra[n]),
            label="Simulated XCAL Spectra",
            color="green",
        )
        plt.plot(
            np.abs(xcal_spectra_processed[n]),
            label="Original XCAL Spectra",
            color="red",
        )
        plt.title(f"{channel.upper()} {mode.upper()} XCAL Spectra {n}")
        plt.xlabel("Frequency (GHz)")
        plt.ylabel("MJy/sr")
        plt.legend()
        plt.savefig(
            f"calibration/output/checks/{channel}_{mode}/{n}_01_xcal_spectra.png",
            bbox_inches="tight",
        )
        plt.close()

        plt.plot(np.abs(processed_spectra[n]), label="Original Spectra")
        plt.plot(np.abs(xcal_spectra_processed[n]), label="XCAL Processed Spectra")
        plt.plot(np.abs(xcal_spectra[n]), label="XCAL Simulated Spectra")
        plt.plot(
            np.abs(xcal_spectra_processed[n] - xcal_spectra[n]), label="XCAL Residuals"
        )
        plt.title(f"{channel.upper()} {mode.upper()} Spectra Comparison {n}")
        plt.xlabel("Frequency (GHz)")
        plt.ylabel("MJy/sr")
        plt.legend()
        plt.savefig(
            f"calibration/output/checks/{channel}_{mode}/{n}_02_spectra_calculations.png",
            bbox_inches="tight",
        )
        plt.close()

        plt.plot(
            np.abs(
                processed_spectra[n] - (xcal_spectra_processed[n] - xcal_spectra[n])
            ),
            label="Spectrum without noise",
        )
        plt.plot(np.abs(simulated_spectra[n]), label="Simulated Spectra")
        plt.title(f"{channel.upper()} {mode.upper()} Spectra {n}")
        plt.xlabel("Frequency (GHz)")
        plt.ylabel("MJy/sr")
        plt.legend()
        plt.savefig(
            f"calibration/output/checks/{channel}_{mode}/{n}_03_spectra_comparison.png",
            bbox_inches="tight",
        )
        plt.close()

        print(
            f"original ifg peak: {np.nanmax(np.abs(original_ifgs[n]) - np.median(original_ifgs[n]))}"
        )
        print(f"simulated ifg peak: {np.nanmax(np.abs(simulated_ifgs[n]))}")
        print(
            f"ratio (original/simulated): {(np.max(np.abs(original_ifgs[n])) - np.median(original_ifgs[n]))/ np.nanmax(np.abs(simulated_ifgs[n]))}"
        )

        # for ifg in range(simulated_ifgs.shape[0]):
        for ifg in range(n, n + 1, 1):
            # scale simulated ifg by the value at peak position
            # simulated_ifgs[ifg, :] = simulated_ifgs[ifg, :] * (
            #     (original_ifgs[ifg, peak_positions[f"{channel}_{mode}"]] - np.median(original_ifgs[ifg, :]))
            #     / simulated_ifgs[ifg, peak_positions[f"{channel}_{mode}"]]
            # )

            print(f"Plotting IFG {ifg+1} for {channel.upper()} {mode.upper()}...")

            # plot both and residuals
            plt.subplots(4, 1, figsize=(10, 15), sharex=True, sharey=True)
            plt.suptitle(
                f"{channel.upper()} {mode.upper()} IFG {ifg}\nTemps: XCAL={xcal[ifg]:.2f}, ICAL={ical[ifg]:.2f}, dihed={dihedral[ifg]:.2f}, refhorn={refhorn[ifg]:.2f}, skyhorn={skyhorn[ifg]:.2f}, collimator={collimator[ifg]:.2f}, bolometer={bolometer[ifg]:.2f}"
            )

            plt.subplot(4, 1, 1)
            plt.plot(original_ifgs[ifg, :] - np.median(original_ifgs[ifg, :]))
            plt.axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            plt.title("Original IFG")

            plt.subplot(4, 1, 2)
            plt.plot(simulated_ifgs[ifg, :])
            plt.axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            plt.title("Simulated IFG")

            plt.subplot(4, 1, 3)
            plt.plot(
                original_ifgs[ifg, :] - np.median(original_ifgs[ifg, :]),
                label="Original IFG",
            )
            plt.plot(simulated_ifgs[ifg, :], label="Simulated IFG")
            plt.axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            plt.title("Original vs Simulated IFG")
            plt.legend()

            plt.subplot(4, 1, 4)
            plt.plot(
                original_ifgs[ifg, :]
                - np.median(original_ifgs[ifg, :])
                - simulated_ifgs[ifg, :]
            )
            plt.plot(
                [
                    np.median(
                        original_ifgs[ifg, :]
                        - np.median(original_ifgs[ifg, :])
                        - simulated_ifgs[ifg, :]
                    )
                ]
                * 512
            )
            plt.axvline(
                x=g.PEAK_POSITIONS[f"{channel}_{mode}"],
                color="red",
                linestyle="--",
                label="Peak",
            )
            plt.title("Residuals (original - simulated)")

            plt.savefig(
                f"calibration/output/checks/{channel}_{mode}/{ifg}_04_ifg_residuals.png",
                bbox_inches="tight",
            )
            plt.close()
