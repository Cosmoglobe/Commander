import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib import lines

import globals as g
from pipeline import ifg_spec
from utils.config import gen_nyquistl

fnyq = gen_nyquistl(
    "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
)

for channel in g.CHANNELS:
    data = np.load(f"{g.PREPROCESSED_DATA_PATH}cal_{channel}.npz")

    mtm_length = data["mtm_length"][:]
    mtm_speed = data["mtm_speed"][:]

    for mode in g.MODES:
        if channel[1] == "h" and mode == "lf":
            continue

        fits_data = fits.open(
            f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
        )
        apod = fits_data[1].data["APODIZAT"][0]
        # apod = np.ones(512, dtype=np.float64)  # No apodization for now

        if mode[0] == "s":
            length_filter = mtm_length == 0
        else:
            length_filter = mtm_length == 1
        if mode[1] == "s":
            speed_filter = mtm_speed == 0
        else:
            speed_filter = mtm_speed == 1

        mode_filter = length_filter & speed_filter

        ifg = data["ifg"][mode_filter]
        ifg = ifg - np.median(ifg, axis=1, keepdims=True)
        adds_per_group = data["adds_per_group"][mode_filter]
        bol_cmd_bias = data["bol_cmd_bias"][mode_filter]
        bol_volt = data["bol_volt"][mode_filter]
        gain = data["gain"][mode_filter]
        sweeps = data["sweeps"][mode_filter]
        bolometer = data["bolometer"][mode_filter]

        frec = 4 * (g.CHANNELS[channel] % 2) + g.MODES[mode]

        fits_data = fits.open(
            f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
        )
        otf = fits_data[1].data["RTRANSFE"][0] + 1j * fits_data[1].data["ITRANSFE"][0]
        otf = otf[np.abs(otf) > 0]

        spec = ifg_spec.ifg_to_spec(
            ifg=ifg,
            channel=channel,
            mode=mode,
            adds_per_group=adds_per_group,
            bol_cmd_bias=bol_cmd_bias,
            bol_volt=bol_volt,
            fnyq_icm=fnyq["icm"][frec],
            otf=otf,
            Tbol=bolometer,
            apod=apod,
            gain=gain,
            sweeps=sweeps,
        )

        ifg_replicated = ifg_spec.spec_to_ifg(
            spec=spec,
            channel=channel,
            mode=mode,
            adds_per_group=adds_per_group,
            bol_cmd_bias=bol_cmd_bias,
            bol_volt=bol_volt,
            Tbol=bolometer,
            gain=gain,
            sweeps=sweeps,
            otf=otf,
            apod=apod,
            fnyq_icm=fnyq["icm"][frec],
        )

        # n = np.random.randint(0, ifg.shape[0])
        # n = 1663
        n = 1704
        # n = 2091

        plt.figure(figsize=(10, 6))
        plt.plot(ifg[n].real, label="Original IFG (real)")
        plt.plot(
            ifg_replicated[n].real,
            label="Replicated IFG (real)",
        )
        plt.plot(
            ifg[n].imag, label="Original IFG (imag)", alpha=0.5, linestyle="dashed"
        )
        plt.plot(
            ifg_replicated[n].imag,
            label="Replicated IFG (imag)",
            linestyle="dashed",
            alpha=0.5,
        )
        plt.legend()
        plt.title(f"IFG Replication for {channel.upper()}{mode.upper()} (Sample {n})")
        plt.xlabel("Sample")
        plt.ylabel("Amplitude")
        plt.grid()
        plt.savefig(f"calibration/output/ifg_to_spec_to_ifg/{channel}_{mode}_{n}.png")
        plt.close()
