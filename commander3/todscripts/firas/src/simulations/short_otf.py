import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import globals as g
import utils.my_utils as utils
from utils.config import gen_nyquistl
from utils.frd import elex_transfcnl
from utils.fut import get_recnum
from calibration import bolometer

# TODO: generalise for all channels and modes
channel = "ll"
mode = "ss"

data = np.load(g.PROCESSED_DATA_PATH_CAL)

ifgs_ll = data["ifg_ll_ss"]

n = np.random.randint(0, ifgs_ll.shape[0])

temp_xcal = data["xcal_ss"]
temp_ical = data["ical_ss"]
temp_dihedral = data["dihedral_ss"]
temp_refhorn = data["refhorn_ss"]
temp_skyhorn = data["skyhorn_ss"]
temp_collimator = data["collimator_ss"]
temp_bolometer_ll = data["bolometer_ll_ss"]
temp_bolometer_lh = data["bolometer_lh_ss"]
temp_bolometer_rl = data["bolometer_rl_ss"]
temp_bolometer_rh = data["bolometer_rh_ss"]

frequencies_ll = utils.generate_frequencies("ll", "ss")

bb_xcal = utils.planck(frequencies_ll, np.array(temp_xcal))
bb_ical = utils.planck(frequencies_ll, np.array(temp_ical))
bb_dihedral = utils.planck(frequencies_ll, np.array(temp_dihedral))
bb_refhorn = utils.planck(frequencies_ll, np.array(temp_refhorn))
bb_skyhorn = utils.planck(frequencies_ll, np.array(temp_skyhorn))
bb_collimator = utils.planck(frequencies_ll, np.array(temp_collimator))
bb_bolometer_ll = utils.planck(frequencies_ll, np.array(temp_bolometer_ll))
bb_bolometer_lh = utils.planck(frequencies_ll, np.array(temp_bolometer_lh))
bb_bolometer_rl = utils.planck(frequencies_ll, np.array(temp_bolometer_rl))
bb_bolometer_rh = utils.planck(frequencies_ll, np.array(temp_bolometer_rh))

fits_data = fits.open(f"{g.PUB_MODEL}FIRAS_CALIBRATION_MODEL_LLSS.FITS")
apod = fits_data[1].data["APODIZAT"][0]

emiss_xcal = (
    fits_data[1].data["RTRANSFE"][0][:43] + 1j * fits_data[1].data["ITRANSFE"][0][:43]
)
emiss_ical = (
    fits_data[1].data["RICAL"][0][:43] + 1j * fits_data[1].data["IICAL"][0][:43]
)
emiss_dihedral = (
    fits_data[1].data["RDIHEDRA"][0][:43] + 1j * fits_data[1].data["IDIHEDRA"][0][:43]
)
emiss_refhorn = (
    fits_data[1].data["RREFHORN"][0][:43] + 1j * fits_data[1].data["IREFHORN"][0][:43]
)
emiss_skyhorn = (
    fits_data[1].data["RSKYHORN"][0][:43] + 1j * fits_data[1].data["ISKYHORN"][0][:43]
)
emiss_bolometer = (
    fits_data[1].data["RBOLOMET"][0][:43] + 1j * fits_data[1].data["IBOLOMET"][0][:43]
)
emiss_collimator = (
    fits_data[1].data["RSTRUCTU"][0][:43] + 1j * fits_data[1].data["ISTRUCTU"][0][:43]
)

spec = (
    bb_xcal
    + (
        bb_ical * emiss_ical
        + bb_dihedral * emiss_dihedral
        + bb_refhorn * emiss_refhorn
        + bb_skyhorn * emiss_skyhorn
        + bb_bolometer_ll * emiss_bolometer
        + bb_bolometer_lh * emiss_bolometer
        + bb_bolometer_rl * emiss_bolometer
        + bb_bolometer_rh * emiss_bolometer
        + bb_collimator * emiss_collimator
    )
    / emiss_xcal
)

plt.plot(frequencies_ll, np.abs(spec[n]))
plt.xlabel("Frequency (GHz)")
plt.ylabel("Intensity (MJy/sr)")
# plt.show()

spec = spec / g.FAC_ERG_TO_MJY

frec = 4 * 3 + 0  # LLSS
fnyq_icm = gen_nyquistl(
    "../reference/fex_samprate.txt", "../reference/fex_nyquist.txt", "int"
)["hz"][frec]

spec_norm = fnyq_icm * g.FAC_ETENDU * g.FAC_ADC_SCALE
spec = spec * spec_norm

spec = spec * emiss_xcal

bol_cmd_bias = data["bol_cmd_bias_ll_ss"]
bol_volt = data["bol_volt_ll_ss"]

R0, T0, G1, beta, rho, C1, C3, Jo, Jg = bolometer.get_bolometer_parameters(
    channel, mode
)

S0 = bolometer.calculate_dc_response(
    bol_cmd_bias=bol_cmd_bias,
    bol_volt=bol_volt,
    Tbol=temp_bolometer_ll,
    R0=R0,
    T0=T0,
    G1=G1,
    beta=beta,
    rho=rho,
    Jo=Jo,
    Jg=Jg,
)

tau = bolometer.calculate_time_constant(
    C3=C3,
    Tbol=temp_bolometer_ll,
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

afreq = utils.get_afreq(0, 3, 43)  # LLSS
B = 1 + 1j * tau[:, np.newaxis] * afreq[np.newaxis, :]

# bolometer transfer function
spec = S0[:, np.newaxis] * spec / B

etfl_all = elex_transfcnl(samprate=681.43, nfreq=43)  # LLSS
adds_per_group = data["adds_per_group_ss"]
erecno = get_recnum(0, 3, adds_per_group).astype(np.int32)  # LLSS
etf = etfl_all[erecno, :]

spec = spec / etf

ifg = np.fft.irfft(spec)

# unclean ifg
sweeps = data["sweeps_ss"]
gain = data["gain_ll_ss"]

ifg = ifg * gain[:, np.newaxis] * sweeps[:, np.newaxis]

peak_pos = g.PEAK_POSITIONS[f"{channel}_{mode}"]
ifg = np.roll(ifg, -peak_pos)

fig, ax = plt.subplots(3, 1, figsize=(10, 8))
ax[0].plot(ifg[n])
ax[0].set_title(f"Siutilslated IFG for {channel.upper()}{mode.upper()}")
ax[0].set_xlabel("Sample Index")
ax[0].set_ylabel("Amplitude")
ax[0].axvline(x=peak_pos, color="red", linestyle="--", label="Peak Position")
ax[0].legend()
ax[1].plot(ifgs_ll[n] - np.median(ifgs_ll[n]))
ax[1].set_title(f"Original IFG for {channel.upper()}{mode.upper()}")
ax[1].set_xlabel("Sample Index")
ax[1].set_ylabel("Amplitude")
ax[1].axvline(x=peak_pos, color="red", linestyle="--", label="Peak Position")
ax[1].legend()
ax[2].plot(ifg[n] - (ifgs_ll[n] - np.median(ifgs_ll[n])))
ax[2].set_title(f"Residual IFG for {channel.upper()}{mode.upper()}")
ax[2].set_xlabel("Sample Index")
ax[2].set_ylabel("Amplitude")
ax[2].axvline(x=peak_pos, color="red", linestyle="--", label="Peak Position")
ax[2].legend()
plt.tight_layout()
plt.show()
