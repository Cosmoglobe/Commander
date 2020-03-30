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

mnem_list = np.loadtxt('mnem_indices.txt', delimiter=',', dtype=str)

def aihk_mnemonic(mnem):
    for i in range(len(mnem_list)):
        if mnem in mnem_list[i,0]:
            ind = i
    index, arr = mnem_list[ind,1:]
    # index - 1 because the indexig is from Fortran
    return int(index)-1, int(arr)


def get_resist(data, index, arr):
    Win = 0
    Wflg = 0
    if arr == 3:
        Tmp = data[3].data['AEU2'][:, index]
        if (index >= 1) and (index <= 31):
            Wflg = 1
            Win = data[3].data['AWIN2'][:,index//2]
            # This is some kind of hexadecimal code, but it seems like it's a
            # correction on the order of 30 K. For the sake of my sanity, I'm going
            # to set terms involving Win to zero.
            if index % 2 == 0:
                Win = np.array([w & 0x00FF for w in Win])
            else:
                Win = np.array([(w & 0xFF00)*2**-8 for w in Win])
            Slp  = 254.968244e-6 
            YInt = 319.5004
            WInt = 129
            WSlp = 256.0
            Rmax = 650.25838
            Rmax = 0
    elif arr == 2:
        Tmp = data[3].data['AEU1'][:, index]
        # For AEU1:
        if (index >= 1) and (index <= 31):
            Wflg = 2
            Win = data[3].data['AWIN1'][:,index//2]
            if index % 2 == 0:
                Win = np.array([w & 0x00FF for w in Win])
            else:
                Win = np.array([(w & 0xFF00)*2**-8 for w in Win])
            Slp  = 255.381467e-6
            YInt = 319.5197
            WInt = 129
            WSlp = 256.
            Rmax = 650.58226
    else:
        Tmp = data[3].data['PDU'][index]
    if Wflg == 0:
        return Tmp
    else:
        Res = Tmp*Slp + YInt + Rmax*(Win - WInt) / WSlp

        # from aihk_mnem2serial_list.pro
        SerialNum = 'UG89'
        SensorId  = 'IHK1_AC_17'
        # from get_prt_temp.pro, converts resistance in ohms to temperature in Kelvin.
        # This requires using the prt_splinedata.xdr file in the ref directory.A
        
        Temp = get_prt_temp(Res, SerialNum)
        return Temp

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


def get_val_from_mnem(data, mnem):
    # there are several steps to get from the raw TOD to physical units.
    
    # First read in the data themselves, which are indexed according to the
    # mnemonics mentioned before.
    
    
    # Second, convert the raw value into physical units using the polynomial
    # conversion in line 122 of aihk_get_mnemonic.pro.
    
    # Finally, you need to get the serial number of the thermistor used to measure
    # the temperature, and get that thermistor's conversion to temperature.
    
    # I'll do this manually as an example, for K1 A side feed.
    
    # Thinking of doing something like K1A and then.
    #mnem = get_mnem(band, variable)
    
    # index = 18 in fortran indexing, so 17 in idl. Also, AEU1.
    
    index, arr = aihk_mnemonic(mnem)
    # index, arr = 17, 2
    
    
    # Get temperature from index
    Temp = get_resist(data, index, arr)
    
    coefs_array = np.loadtxt('mnem_coefs.txt', dtype=str, delimiter=',')
    for i in range(len(coefs_array)):
        if mnem in coefs_array[i,0]:
            ind = i
    coefs = coefs_array[ind,1:].astype('float')
    T_pow = np.array([Temp**i for i in range(len(coefs))])
    Val = coefs.dot(T_pow)

    return Val


if __name__ == '__main__':

    # to convert to JD, add 2.4500e06
    data = fits.open('wmap.fits')




    # The only published gain model I've seen in the WMAP paper is Figure 2 of
    # Jarosik et al. 2007, which plots the V223 detector.
    T0 = 5.4309e1
    V0 = -5.4468e-1
    beta = -2.7406e-3
    alpha = -1.6287e1
    m = -3.6946e-7
    c = 2.5173e-1
    pars = np.array([T0, V0, beta, alpha, m ,c])

    mnem_rfb = 'DRV223RFBI14'
    mnem_tfpa = 'DFV22FPATEET'
    mnem_trxb = 'DRV222RXBAMPT'

    T_RXB = get_val_from_mnem(data, mnem_trxb)
    T_FPA = get_val_from_mnem(data, mnem_tfpa)
    RF_bias = get_val_from_mnem(data, mnem_rfb)
    t_JD = data[3].data['TIME'] + 2.45e6
    # 2452131.5 is 0:00 of day 222 of 2001.
    t = t_JD - 2452131.50000

    G_V223 = G(RF_bias, T_RXB, T_FPA, t, pars)

    fig, axes = plt.subplots(sharex=True, nrows=2, ncols=2)
    axs = axes.flatten()
    axs[0].plot(t, T_RXB, '.')
    axs[0].set_title(r'$T_\mathrm{RXB}$')
    axs[1].plot(t, T_FPA, '.')
    axs[1].set_title(r'$T_\mathrm{FPA}$')
    axs[2].plot(t, RF_bias, '.')
    axs[2].set_title('RF Bias')
    axs[3].plot(t, G_V223, '.')
    axs[3].set_title('Gain')
    fig.suptitle('V223')


    # Figure 4 in Hinshaw et al has a couple of figures. K113 and V113;
    # K113
    T0 = 5.3235e1
    V0 = -5.5125e-1
    beta = -7.220e-4
    alpha = 4.4619e1
    m = -2.7274e-7
    c = -6.7556e-1
    pars = np.array([T0, V0, beta, alpha, m ,c])
    mnem_trxb = 'DRK12RXBRIBT'
    mnem_tfpa = 'DFK1AFEEDT'
    mnem_rfb = 'DRK113RFBI0'

    T_RXB = get_val_from_mnem(data, mnem_trxb)
    T_FPA = get_val_from_mnem(data, mnem_tfpa)
    RF_bias = get_val_from_mnem(data, mnem_rfb)
    G_K113 = G(RF_bias, T_RXB, T_FPA, t, pars)

    fig, axes = plt.subplots(sharex=True, nrows=2, ncols=2)
    axs = axes.flatten()
    axs[0].plot(t, T_RXB, '.')
    axs[0].set_title(r'$T_\mathrm{RXB}$')
    axs[1].plot(t, T_FPA, '.')
    axs[1].set_title(r'$T_\mathrm{FPA}$')
    axs[2].plot(t, RF_bias, '.')
    axs[2].set_title('RF Bias')
    axs[3].plot(t, G_K113, '.')
    axs[3].set_title('Gain')
    fig.suptitle('K113')

    #V113
    T0 = 4.4834e1
    V0 = -1.0074e0
    beta = -2.3861e-3
    alpha = -2.1649e1
    m = -2.1019e-7
    c = 4.85e-1
    pars = np.array([T0, V0, beta, alpha, m ,c])
    mnem_trxb = 'DRV111RXBAMPT'
    mnem_tfpa = 'DFV11FPATEET'
    mnem_rfb = 'DRV113RFBI32'

    T_RXB = get_val_from_mnem(data, mnem_trxb)
    T_FPA = get_val_from_mnem(data, mnem_tfpa)
    RF_bias = get_val_from_mnem(data, mnem_rfb)
    G_V113 = G(RF_bias, T_RXB, T_FPA, t, pars)

    fig, axes = plt.subplots(sharex=True, nrows=2, ncols=2)
    axs = axes.flatten()
    axs[0].plot(t, T_RXB, '.')
    axs[0].set_title(r'$T_\mathrm{RXB}$')
    axs[1].plot(t, T_FPA, '.')
    axs[1].set_title(r'$T_\mathrm{FPA}$')
    axs[2].plot(t, RF_bias, '.')
    axs[2].set_title('RF Bias')
    axs[3].plot(t, G_V113, '.')
    axs[3].set_title('Gain')
    fig.suptitle('V113')

    plt.show()
