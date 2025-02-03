import numpy as np
from astropy.io import fits
from my_utils import ifg_to_spec, planck
from utils.config import gen_nyquistl

data_path = "/mn/stornext/d16/cmbco/ola/firas/coadded_interferograms/"

channels = {"rh": 0, "rl": 1, "lh": 2, "ll": 3}
modes = {"ss": 0, "lf": 3}

data = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and channel in ["rh", "lh"]):
            data[f"{channel}_{mode}"] = fits.open(
                data_path
                + f"FIRAS_COADDED_SKY_INTERFEROGRAMS_{channel.upper()}"
                + mode.upper()
                + ".FITS"
            )[1].data

fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)
frec = {}
for channel, channel_value in channels.items():
    for mode, mode_value in modes.items():
        if not (mode == "lf" and channel in ["rh", "lh"]):
            frec[f"{channel}_{mode}"] = 4 * (channel_value % 2) + mode_value

fits_data = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and channel in ["rh", "lh"]):
            fits_data[f"{channel}_{mode}"] = fits.open(
                f"/mn/stornext/d16/cmbco/ola/firas/pub_calibration_model/FIRAS_CALIBRATION_MODEL_{channel.upper()}{mode.upper()}.FITS"
            )

ifg = {}
mtm_speed = {}
channel = 3
bol_cmd_bias = {}
bol_volt = {}
otf = {}
tau = {}
Jo = {}
Jg = {}
Tbol = {}
T0 = {}
R0 = {}
rho = {}
G1 = {}
beta = {}

for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and channel in ["rh", "lh"]):
            ifg[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]["COADDED_"]
            mtm_speed[mode] = 0 if mode == "ss" else 1
            adds_per_group = 3 if mode == "ss" else 8
            bol_cmd_bias[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]["BOLOM_BI"]
            bol_volt[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]["BOLOM_VO"]
            # print(f"{channel}_{mode}")
            # print(f"bol_cmd_bias: {bol_cmd_bias[f'{channel}_{mode}']}")
            # print(f"bol_volt: {bol_volt[f'{channel}_{mode}']}")

            otf[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RTRANSFE"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["ITRANSFE"][0]
            )
            otf[f"{channel}_{mode}"] = otf[f"{channel}_{mode}"][
                np.abs(otf[f"{channel}_{mode}"]) > 0
            ]
            tau[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "TIME_CON"
            ][0]
            Jo[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM8"
            ][0]
            Jg[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM9"
            ][0]
            Tbol[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLOM_B2"
            ][0]
            T0[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM2"
            ][0]
            R0[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM_"
            ][0]
            rho[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM5"
            ][0]
            G1[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM3"
            ][0]
            beta[f"{channel}_{mode}"] = fits_data[f"{channel}_{mode}"][1].data[
                "BOLPARM4"
            ][0]

# print(f"ifg: {ifg}")

spec = {}
for channel, channel_value in channels.items():
    for mode in modes.keys():
        if not (mode == "lf" and channel in ["rh", "lh"]):
            # print(f"{channel}_{mode}")
            spec[f"{channel}_{mode}"] = ifg_to_spec(
                ifg=ifg[f"{channel}_{mode}"],
                mtm_speed=0 if mode == "ss" else 1,
                channel=channel_value,
                adds_per_group=adds_per_group,
                bol_cmd_bias=bol_cmd_bias[f"{channel}_{mode}"],
                bol_volt=bol_volt[f"{channel}_{mode}"],
                fnyq_icm=fnyq["icm"][frec[f"{channel}_{mode}"]],
                fnyq_hz=fnyq["hz"][frec[f"{channel}_{mode}"]],
                otf=otf[f"{channel}_{mode}"],
                Jo=Jo[f"{channel}_{mode}"],
                Jg=Jg[f"{channel}_{mode}"],
                Tbol=Tbol[f"{channel}_{mode}"],
                rho=rho[f"{channel}_{mode}"],
                R0=R0[f"{channel}_{mode}"],
                T0=T0[f"{channel}_{mode}"],
                beta=beta[f"{channel}_{mode}"],
                G1=G1[f"{channel}_{mode}"],
                tau=tau[f"{channel}_{mode}"],
            )
            # print(f"spec: {spec[f'{channel}_{mode}']}")

# print(f"spec: {spec}")

# frequency mapping
nu0 = {"ss": 68.020812, "lf": 23.807283}
dnu = {"ss": 13.604162, "lf": 3.4010405}
nf = {"lh_ss": 210, "ll_lf": 182, "ll_ss": 43, "rh_ss": 210, "rl_lf": 182, "rl_ss": 43}

f_ghz = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and (channel == "lh" or channel == "rh")):
            f_ghz[f"{channel}_{mode}"] = np.linspace(
                nu0[mode],
                nu0[mode] + dnu[mode] * (nf[f"{channel}_{mode}"] - 1),
                nf[f"{channel}_{mode}"],
            )

ical_temp = {}
dihedral_temp = {}
skyhorn_temp = {}
refhorn_temp = {}
mirror_temp = {}
bath_temp = {}
xcal_temp = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and channel in ["rh", "lh"]):
            ical_temp[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]["ICAL_TEM"]
            dihedral_temp[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]["DIHEDRAL"]
            skyhorn_temp[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]["SKYHORN_"]
            refhorn_temp[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]["REFHORN_"]
            mirror_temp[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]["MIRROR_T"]
            bath_temp[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]["BATH_TEM"]
            xcal_temp[f"{channel}_{mode}"] = data[f"{channel}_{mode}"]["XCAL_TEM"]

bb_ical = {}
ical_emiss = {}
bb_dihedral = {}
dihedral_emiss = {}
bb_skyhorn = {}
skyhorn_emiss = {}
bb_refhorn = {}
refhorn_emiss = {}
bb_bol = {}
bol_emiss = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and channel in ["rh", "lh"]):
            bb_ical[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"][: len(spec[f"{channel}_{mode}"])],
                ical_temp[f"{channel}_{mode}"],
            )
            ical_emiss[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RICAL"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["IICAL"][0]
            )
            ical_emiss[f"{channel}_{mode}"] = ical_emiss[f"{channel}_{mode}"][
                np.abs(ical_emiss[f"{channel}_{mode}"]) > 0
            ]

            bb_dihedral[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"][: len(spec[f"{channel}_{mode}"])],
                dihedral_temp[f"{channel}_{mode}"],
            )
            dihedral_emiss[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RDIHEDRA"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["IDIHEDRA"][0]
            )
            dihedral_emiss[f"{channel}_{mode}"] = dihedral_emiss[f"{channel}_{mode}"][
                np.abs(dihedral_emiss[f"{channel}_{mode}"]) > 0
            ]

            bb_skyhorn[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"][: len(spec[f"{channel}_{mode}"])],
                skyhorn_temp[f"{channel}_{mode}"],
            )
            skyhorn_emiss[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RSKYHORN"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["ISKYHORN"][0]
            )
            skyhorn_emiss[f"{channel}_{mode}"] = skyhorn_emiss[f"{channel}_{mode}"][
                np.abs(skyhorn_emiss[f"{channel}_{mode}"]) > 0
            ]

            bb_refhorn[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"][: len(spec[f"{channel}_{mode}"])],
                refhorn_temp[f"{channel}_{mode}"],
            )
            refhorn_emiss[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RREFHORN"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["IREFHORN"][0]
            )
            refhorn_emiss[f"{channel}_{mode}"] = refhorn_emiss[f"{channel}_{mode}"][
                np.abs(refhorn_emiss[f"{channel}_{mode}"]) > 0
            ]

            # bolometer spectrum
            bb_bol[f"{channel}_{mode}"] = planck(
                f_ghz[f"{channel}_{mode}"][: len(spec[f"{channel}_{mode}"])],
                bath_temp[f"{channel}_{mode}"],
            )
            bol_emiss[f"{channel}_{mode}"] = (
                fits_data[f"{channel}_{mode}"][1].data["RBOLOMET"][0]
                + 1j * fits_data[f"{channel}_{mode}"][1].data["IBOLOMET"][0]
            )
            bol_emiss[f"{channel}_{mode}"] = bol_emiss[f"{channel}_{mode}"][
                np.abs(bol_emiss[f"{channel}_{mode}"]) > 0
            ]

sky = {}
for channel in channels.keys():
    for mode in modes.keys():
        if not (mode == "lf" and channel in ["rh", "lh"]):
            # print(f"{channel}_{mode}")
            # print(
            #     f"size of bb_ical: {bb_ical[f'{channel}_{mode}'].shape} and ical_emiss: {ical_emiss[f'{channel}_{mode}'].shape}"
            # )
            sky[f"{channel}_{mode}"] = (
                spec[f"{channel}_{mode}"]
                - (
                    (bb_ical[f"{channel}_{mode}"] * ical_emiss[f"{channel}_{mode}"])
                    + (
                        bb_dihedral[f"{channel}_{mode}"]
                        * dihedral_emiss[f"{channel}_{mode}"]
                    )
                    + (
                        bb_refhorn[f"{channel}_{mode}"]
                        * refhorn_emiss[f"{channel}_{mode}"]
                    )
                    + (
                        bb_skyhorn[f"{channel}_{mode}"]
                        * skyhorn_emiss[f"{channel}_{mode}"]
                    )
                    + (bb_bol[f"{channel}_{mode}"] * bol_emiss[f"{channel}_{mode}"])
                )
                / otf[f"{channel}_{mode}"][np.newaxis, :]
            )

# print(f"sky: {sky}")

# save the sky
np.savez(
    "../../output/data/lambda_sky.npz",
    **sky,
)
