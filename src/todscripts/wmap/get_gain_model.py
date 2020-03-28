import numpy as np
from astropy.io import fits
from scipy.io import readsav
from scipy.interpolate import BSpline

import matplotlib.pyplot as plt

'''
The gain model as described by the explanatory supplement requires the total
radiometer bias power, T_FPA, and T_RXB as inputs.

These quantities are indexed by mnemonics, enumerated in table C.2.2, C.2.4, and
C.2.5, respectively. I have listed them in the files t_rxb_mnem.txt,
t_fpa_mnem.txt, and rfbias.txt

AIHK_mnemonic.pro lists the conversion between mnemonic, and the location in the
fits table ("index" and "arr").


AIHK_GetMnemonic.pro returns the physical value associated with a mnemonic.
'''

def get_prt_temp(resist, SerialNum):
    table = readsav('wmap_routines/ref/prt_data/prt_splinedata.xdr')
    splinedata = table['splinedata']
    serial_nos = []
    for b in splinedata['serial']:
        serial_nos.append(bytearray(b).decode('utf8').strip())
    element = np.where(np.array(serial_nos) == SerialNum)[0]
    if len(element) == 0:
        print(f'Serial number {SerialNum} not found, recheck')
        return
    
    
    # spline structure for this PRT
    a = splinedata[element][0]
    # using notation of
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html
    t = a['xknot']
    c = a['bscoef']
    k = a['korder']
    maxind = min(sum(c !=0), sum(t!=0))
    t = t[:maxind]
    c = c[:maxind]
    spl = BSpline(t,c,k)

    temperature = spl(resist)

    return temperature


def G(V, T_RXB, T_FPA, t, pars):
    T0, V0, beta, alpha, m, c = pars
    G = alpha*(V - V0 - beta*(T_RXB-290))/(T_FPA-T0) + m*t + c
    return G


# there are several steps to get from the raw TOD to physical units.

# First read in the data themselves, which are indexed according to the
# mnemonics mentioned before.


# Second, convert the raw value into physical units using the polynomial
# conversion in line 122 of aihk_get_mnemonic.pro.

# Finally, you need to get the serial number of the thermistor used to measure
# the temperature, and get that thermistor's conversion to temperature.

# I'll do this manually as an example, for K1 A side feed.

mnem = 'DFK1AFEEDT'

# index = 18 in fortran indexing, so 17 in idl. Also, AEU1.

# index, arr = aihk_mnemonic(mnem)
index, arr = 17, 2

data = fits.open('wmap.fits')

# Get raw telemtry value
Win = 0
Wflg = 0
if arr == 3:
    Tmp = data[3].data['AEU2'][:, index]
    if (index >= 1) and (index <= 31):
        Wlfg = 1
        Win = data[3].data['AWIN2'][:,index//2]
        # This is some kind of hexadecimal code, but it seems like it's a
        # correction on the order of 30 K. For the sake of my sanity, I'm going
        # to set terms involving Win to zero.
        '''
        if (index % 2) != 1:
            Win = (Win & '00FF'xL)
        else:
            Win = ISHFT((Win & 'FF00'xL), -8)
        '''
        Slp  = 254.968244e-6 
        YInt = 319.5004
        WInt = 129
        WSlp = 256.0
        Rmax = 650.25838
        Rmax = 0
elif arr == 2:
    Tmp = data[3].data['AEU1'][:, index]
    # For AEU1:
    Win = data[3].data['AWIN1'][:,index//2]
    print(Win)
    Slp  = 255.381467e-06
    YInt = 319.5197
    WInt = 129
    WSlp = 256.
    Rmax = 650.58226
    Rmax = 0
else:
    Tmp = data[3].data['PDU'][index]
  


Res = Tmp*Slp + YInt + Rmax*(Win - WInt) / WSlp




# from aihk_mnem2serial_list.pro
SerialNum = 'UG89'
SensorId  = 'IHK1_AC_17'

# from get_prt_temp.pro, converts resistance in ohms to temperature in Kelvin.
# This requires using the prt_splinedata.xdr file in the ref directory.A

Temp = get_prt_temp(Res, SerialNum)

plt.figure()
plt.plot(data[3].data['TIME'], Temp)
plt.xlabel(r'$t$')
plt.ylabel(r'$T$ (Kelvin)')

plt.show()
