"""
Script to fit an optical transfer function to the calibration data. For now, it follows the model:
S_sky/XCAL = 1/OTF * 1/ETF * 1/Bol * Y - R,
where Y is the Fourier transformed interferogram and R is defined by
R = 1/OTF * sum over all emitters i (excluding XCAL) of E_i * P(T_i).
"""

import h5py
import numpy as np

import globals as g
from calibration import bolometer
from utils import my_utils as utils

data = h5py.File(g.PREPROCESSED_DATA_PATH_CAL, "r")["df_data"]
print(data.keys())

channels = {"ll": 3}
modes = {"ss": 0}

# TODO: fit for temp weight coefficients too
xcal = (data["a_xcal"][:] + data["b_xcal"][:]) / 2
ical = (data["a_ical"][:] + data["b_ical"][:]) / 2
dihedral = (data["a_dihedral"][:] + data["b_dihedral"][:]) / 2
refhorn = (data["a_refhorn"][:] + data["b_refhorn"][:]) / 2
skyhorn = (data["a_skyhorn"][:] + data["b_skyhorn"][:]) / 2
collimator = (data["a_collimator"][:] + data["b_collimator"][:]) / 2
bolometer_ll = (data["a_bol_assem_ll"][:] + data["b_bol_assem_ll"][:]) / 2
bolometer_lh = (data["a_bol_assem_lh"][:] + data["b_bol_assem_lh"][:]) / 2
bolometer_rl = (data["a_bol_assem_rl"][:] + data["b_bol_assem_rl"][:]) / 2
bolometer_rh = (data["a_bol_assem_rh"][:] + data["b_bol_assem_rh"][:]) / 2

temps = {
    "xcal": xcal,
    "ical": ical,
    "dihedral": dihedral,
    "refhorn": refhorn,
    "skyhorn": skyhorn,
    "collimator": collimator,
    "bolometer_ll": bolometer_ll,
    "bolometer_lh": bolometer_lh,
    "bolometer_rl": bolometer_rl,
    "bolometer_rh": bolometer_rh,
}

adds_per_group = data["adds_per_group"][:]
sweeps = data["sweeps"][:]

for channel in channels:
    ifg = data[f"ifg_{channel}"][:]
    Y = np.fft.rfft(ifg, axis=1)

    bol_cmd_bias = data[f"bol_cmd_bias_{channel}"][:]
    bol_volt = data[f"bol_volt_{channel}"][:]

    for mode in modes:
        R0, T0, G1, beta, rho, C1, C3, Jo, Jg = bolometer.get_bolometer_parameters(
            channel, mode
        )

        # bolometer response function
        S0 = bolometer.calculate_dc_response(
            bol_cmd_bias=bol_cmd_bias / 25.5,
            bol_volt=bol_volt,
            Jo=Jo,
            Jg=Jg,
            Tbol=temps[f"bolometer_{channel}"],
            rho=rho,
            R0=R0,
            T0=T0,
            beta=beta,
            G1=G1,
        )

        tau = bolometer.calculate_time_constant(
            C3=C3,
            Tbol=temps[f"bolometer_{channel}"],
            C1=C1,
            G1=G1,
            beta=beta,
            bol_volt=bol_volt,
            Jo=Jo,
            Jg=Jg,
            bol_cmd_bias=bol_cmd_bias,
            rho=rho,
            T0=T0,
        )

        omega = utils.get_afreq(0 if mode[1] == "s" else 1, channel, 257)

        B = S0[:, np.newaxis] / (1 + 1j * omega[np.newaxis, :] * tau[:, np.newaxis])

        # electronics transfer function
