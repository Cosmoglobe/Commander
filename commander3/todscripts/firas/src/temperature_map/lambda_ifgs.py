import numpy as np
from astropy.io import fits
from my_utils import ifg_to_spec, planck
from utils.config import gen_nyquistl

data_path = "/mn/stornext/d16/cmbco/ola/firas/coadded_interferograms/"

# modes = ["ss", "lf"]
modes = {"ss": 0, "lf": 3}

data = {}
for mode in modes:
    data[mode] = fits.open(
        data_path + "FIRAS_COADDED_SKY_INTERFEROGRAMS_LL" + mode.upper() + ".FITS"
    )[1].data

fnyq = gen_nyquistl(
    "../../reference/fex_samprate.txt", "../../reference/fex_nyquist.txt", "int"
)
channel = 3  # LL
frec = {}
for mode, item in modes.items():
    frec[mode] = 4 * (channel % 2) + item

fits_data = {}
for mode in modes:
    fits_data[mode] = fits.open(
        f"/mn/stornext/u3/aimartin/d5/firas-reanalysis/Commander/commander3/todscripts/firas/reference/FIRAS_CALIBRATION_MODEL_LL{mode.upper()}.FITS"
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


for mode in modes:
    ifg[mode] = data[mode]["COADDED_"]
    mtm_speed[mode] = 0 if mode == "ss" else 1
    adds_per_group = 3 if mode == "ss" else 8
    bol_cmd_bias[mode] = data[mode]["BOLOM_BI"]
    bol_volt[mode] = data[mode]["BOLOM_VO"]

    otf[mode] = (
        fits_data[mode][1].data["RTRANSFE"][0]
        + 1j * fits_data[mode][1].data["ITRANSFE"][0]
    )
    otf[mode] = otf[mode][np.abs(otf[mode]) > 0]
    tau[mode] = fits_data[mode][1].data["TIME_CON"][0]
    Jo[mode] = fits_data[mode][1].data["BOLPARM8"][0]
    Jg[mode] = fits_data[mode][1].data["BOLPARM9"][0]
    Tbol[mode] = fits_data[mode][1].data["BOLOM_B2"][0]
    T0[mode] = fits_data[mode][1].data["BOLPARM2"][0]
    R0[mode] = fits_data[mode][1].data["BOLPARM_"][0]
    rho[mode] = fits_data[mode][1].data["BOLPARM5"][0]
    G1[mode] = fits_data[mode][1].data["BOLPARM3"][0]
    beta[mode] = fits_data[mode][1].data["BOLPARM4"][0]

# print(f"ifg: {ifg}")

spec = {}
for mode in modes:
    spec[mode] = ifg_to_spec(
        ifg=ifg[mode],
        mtm_speed=0 if mode == "ss" else 1,
        channel=channel,
        adds_per_group=adds_per_group,
        bol_cmd_bias=bol_cmd_bias[mode],
        bol_volt=bol_volt[mode],
        fnyq_icm=fnyq["icm"][frec[mode]],
        fnyq_hz=fnyq["hz"][frec[mode]],
        otf=otf[mode],
        Jo=Jo[mode],
        Jg=Jg[mode],
        Tbol=Tbol[mode],
        rho=rho[mode],
        R0=R0[mode],
        T0=T0[mode],
        beta=beta[mode],
        G1=G1[mode],
        tau=tau[mode],
    )

# print(f"spec: {spec}")

c = 3e8 * 1e2  # cm/s
f_icm = {}
f_ghz = {}
for mode in modes:
    f_icm[mode] = np.arange(210) * (fnyq["icm"][frec[mode]] / 320) + 2
    f_ghz[mode] = f_icm[mode] * c * 1e-9

ical_temp = {}
dihedral_temp = {}
skyhorn_temp = {}
refhorn_temp = {}
mirror_temp = {}
bath_temp = {}
xcal_temp = {}
for mode in modes:
    ical_temp[mode] = data[mode]["ICAL_TEM"]
    dihedral_temp[mode] = data[mode]["DIHEDRAL"]
    skyhorn_temp[mode] = data[mode]["SKYHORN_"]
    refhorn_temp[mode] = data[mode]["REFHORN_"]
    mirror_temp[mode] = data[mode]["MIRROR_T"]
    bath_temp[mode] = data[mode]["BATH_TEM"]
    xcal_temp[mode] = data[mode]["XCAL_TEM"]

# ical spectrum
bb_ical = {}
ical_emiss = {}
for mode in modes:
    print(f"mode: {mode}")
    bb_ical[mode] = planck(f_ghz[mode][: len(spec[mode])], ical_temp[mode])
    ical_emiss[mode] = (
        fits_data[mode][1].data["RICAL"][0] + 1j * fits_data[mode][1].data["IICAL"][0]
    )
    ical_emiss[mode] = ical_emiss[mode][: len(spec[mode])]

    print(
        f"checking sizes: {len(spec[mode])}, {len(bb_ical[mode])}, {len(ical_emiss[mode])}"
    )

# dihedral spectrum
bb_dihedral = {}
dihedral_emiss = {}
for mode in modes:
    bb_dihedral[mode] = planck(f_ghz[mode][: len(spec[mode])], dihedral_temp[mode])
    dihedral_emiss[mode] = (
        fits_data[mode][1].data["RDIHEDRA"][0]
        + 1j * fits_data[mode][1].data["IDIHEDRA"][0]
    )
    dihedral_emiss[mode] = dihedral_emiss[mode][: len(spec[mode])]

# skyhorn spectrum
bb_skyhorn = {}
skyhorn_emiss = {}
for mode in modes:
    bb_skyhorn[mode] = planck(f_ghz[mode][: len(spec[mode])], skyhorn_temp[mode])
    skyhorn_emiss[mode] = (
        fits_data[mode][1].data["RSKYHORN"][0]
        + 1j * fits_data[mode][1].data["ISKYHORN"][0]
    )
    skyhorn_emiss[mode] = skyhorn_emiss[mode][: len(spec[mode])]

# refhorn spectrum
bb_refhorn = {}
refhorn_emiss = {}
for mode in modes:
    bb_refhorn[mode] = planck(f_ghz[mode][: len(spec[mode])], refhorn_temp[mode])
    refhorn_emiss[mode] = (
        fits_data[mode][1].data["RREFHORN"][0]
        + 1j * fits_data[mode][1].data["IREFHORN"][0]
    )
    refhorn_emiss[mode] = refhorn_emiss[mode][: len(spec[mode])]

# bolometer spectrum
bb_bol = {}
bol_emiss = {}
for mode in modes:
    bb_bol[mode] = planck(f_ghz[mode][: len(spec[mode])], bath_temp[mode])
    bol_emiss[mode] = (
        fits_data[mode][1].data["RBOLOMET"][0]
        + 1j * fits_data[mode][1].data["IBOLOMET"][0]
    )
    bol_emiss[mode] = bol_emiss[mode][: len(spec[mode])]

sky = {}
for mode in modes:
    print(f"mode: {mode}")
    print(
        f"checking sizes: {len(spec[mode])}, {len(bb_ical[mode])}, {len(ical_emiss[mode])}, {len(bb_dihedral[mode])}, {len(dihedral_emiss[mode])}, {len(bb_refhorn[mode])}, {len(refhorn_emiss[mode])}, {len(bb_skyhorn[mode])}, {len(skyhorn_emiss[mode])}, {len(bb_bol[mode])}, {len(bol_emiss[mode])}, {len(otf[mode])}"
    )
    sky[mode] = (
        spec[mode]
        - (
            (bb_ical[mode] * ical_emiss[mode])[:, : len(spec[mode][0])]
            + (bb_dihedral[mode] * dihedral_emiss[mode])[:, : len(spec[mode][0])]
            + (bb_refhorn[mode] * refhorn_emiss[mode])[:, : len(spec[mode][0])]
            + (bb_skyhorn[mode] * skyhorn_emiss[mode])[:, : len(spec[mode][0])]
            + (bb_bol[mode] * bol_emiss[mode])[:, : len(spec[mode][0])]
        )
        / otf[mode][np.newaxis, :]
    )

# print(f"sky: {sky}")

# save the sky
np.savez(
    "../../output/data/lambda_sky.npz",
    **sky,
)
