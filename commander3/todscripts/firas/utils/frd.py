'''FIRAS Reference Data: Collection of functions to generate FIRAS reference data
'''

import numpy as np
import utils.fut as fut

fac_fft_length = [640, 640, 640, 640, 160, 160]
fac_spec_length = [321, 321, 321, 321, 81, 81]
fac_spec_length = [512, 321, 321, 321, 81, 81]
fac_nyq_correct = 1.00159

def apod_fcnl(ictr, good_peak, short, xhsf, xhlf):
    '''This routine generates a 512-point apodization function whose shape depends only on the
    position of the IFG peak. The routine first forms a function in which a symmetric region
    about the peak, equal in length to the short side of the stroke, is "averaged", whereas
    the remaining long part of the stroke is given full weight. The routine then multiplies 
    this function by a smoothing function, which is an inverted octic centered at the peak
    position and tapering to zero at the long end of the long side of the stroke.

    This routine zeros out the first two points of the IFG to compensate for the effects
    of the digital transient response subtraction
    '''

    #Initialize the apodization function
    apod_fcn = np.zeros(512)

    halflim = 256
    twidth = 30
    upperlim = 511
    lowerlim = 2

    if short:
        upperlim = 127
        lowerlim = 2
        halflim = 64
        twidth = 8
        if ictr > 127: good_peak = False
    elif xhsf:
        lowerlim = 4
    elif xhlf:
        upperlim = 508
        lowerlim = 1

    apod_fcn[lowerlim:upperlim+1] = 1.0

    #NOTE: Indices have been checked using IPython to input both the Fortran
    #indices and the Python indices to make sure they are always off by one
    #assuming ictr, lowerlim, upperlim, halflim ar always off by one between
    #Fortran and Python and twidth is the same. It would always be good to check
    #again.

    #Average the apodization function according to the peak position
    if not good_peak:
        #If the peak position is bad, set it to upperlim and do not average the function
        if ictr > upperlim:
            ictr = upperlim
        #x = ictr - 2.0
        x = ictr - 1.0 #since x is a RHS divisor and ictr is off-by-one compared to Fortran
    else:
        #Average the function
        #NOTE: I don't think I need to change the values for x since there are differences
        #in the calculations and x is a RHS divisor independent of the index
        if ictr < halflim:
            x = 1.0 + upperlim - ictr
            for j in range(twidth-1):
                apod_fcn[j + lowerlim] = 0.5 * (1.0 - np.cos((j+1) * np.pi/twidth))
                apod_fcn[j + int(2*ictr) - lowerlim - twidth + 2] = 0.5 * (3.0 - np.cos((j+1)*np.pi/twidth))
            for j in range(int(2*ictr) - lowerlim + 1, upperlim+1):
                apod_fcn[j] = 2.0
        else:
            x = ictr + 1.0 - lowerlim
            for j in range(twidth-1):
                apod_fcn[j + int(2*ictr) - upperlim] = 0.5 * (3.0 + np.cos((j+1)*np.pi/twidth))
                apod_fcn[j + upperlim - twidth + 2] = 0.5 * (1.0 + np.cos((j+1)*np.pi/twidth))
            for j in range(lowerlim, int(2*ictr) - upperlim):
                apod_fcn[j] = 2.0

    #Smooth the apodization function. This is the only place "x" appears
    #in this function
    for j in range(lowerlim, upperlim+1):
        apod_fcn[j] *= (1.0 - ((j-ictr)/x)**4)**2

    return apod_fcn

def apodl():
    '''This function generates the apodization functions for coadded IFGs. The apodization functions are
    512-point floating point arrays that have not been rotated to put the peak in the first position
    '''

    apodl_all = np.zeros([433, 512])
    #Cycle through the instrument states. This is for the "regular" length 512 apodization functions, so
    #set short = 0
    short = 0
    nrec = 0
    for ngroup in range(1, 13):
        for mtm_speed in range(2):
            for mtm_length in range(2):
                for chan in range(1, 3): #low/high (LL/RL, LH/RH)
                    xhsf = (chan - 1) * (1 - mtm_length) * mtm_speed
                    xhlf = (chan - 1) * mtm_length * mtm_speed
                    for upmode in range(1, 3):
                        for linearized in range(2):
                            #For each instrument state, determine the IFG peak position, then generate
                            #the apoziation function
                            ictr, rstatus = fut.default_peak(mtm_speed, mtm_length, chan, ngroup, upmode, linearized)

                            #Peak shift for High Short Fast data
                            if xhsf:
                                if linearized == 0 and upmode == 2:
                                    peak = ictr + 0.5
                                else:
                                    peak = ictr - 0.5
                            else:
                                peak = ictr

                            apod_fcn = apod_fcnl(peak, rstatus, short, xhsf, xhlf)

                            #Write the apodization function to the output array
                            apodl_all[nrec, :] = apod_fcn
                            nrec += 1

    #Cycle through the instrument stats. This is for the short length 128 apodization functions,
    #so set short=1. Manufacture these only for special long fast cases (special short fast data
    #will have been decimated and smoothed to match special long fast)
    xhsf = 0
    xhlf = 0
    short = 1
    for ngroup in range(1, 13):
        mtm_speed = 1 #fast scan speed
        mtm_length = 1 #long scan length (correct for peak retrieval)
        chan = 1 #low channel
        for upmode in range(1, 3):
            for linearized in range(2):
                #For each instrument state, determine the IFG peak position, then
                #generate the apodization function
                ictr, rstatus = fut.default_peak(mtm_speed, mtm_length, chan, ngroup, upmode, linearized)
                peak = ictr
                apod_fcn = apod_fcnl(peak, rstatus, short, xhsf, xhlf)

                apodl_all[nrec, :] = apod_fcn
                nrec += 1

    #Generate an apodization function full of 1's for fakeit data
    apodl_all[nrec, :] = np.ones(512)
    nrec += 1

    if nrec != 433:
        raise ValueError("Number of records in apodl_all is not correct")

    return apodl_all

def elex_transfcnl(samprate, nfreq):
    '''This routine creates an array with the FIRAS post-detector electronics transfer function. The function
    is derived from the analytical expressions of the various analog and digital filters; these are
    (1) a treble boost, (2) a 5-pole Bessel filter, (3) a 6-pole DC-blocking filter, (4) a digital lowpass filter,
    , and (5) smoothing due to data compression.

    Each record is an array of 257 complex numbers representing the complex trasnfer function for a particular
    scna mode at frequencies ranging from DC up to the Nyquist frequency associated with the data compression.

    There are 288 records. These are:
        3 samples rates (fakeit, fast scan, and slow scan)
        4 channels
        2 microprocessor science modes (digital filters on or off)
        12 possible amounts of data compression ("adds per group" = 1 to 12)
    '''
    
    ztrans = np.zeros([288, nfreq], dtype=np.complex128)

    tau = [0.00647, 0.03619, 0.00722, 0.04022] #RH, RL, LH, LL treble boost
    bes3db = 100.0 #Bessel corner freq
    fixed_gain = 31.0 * 1.3823 #Preamp fixed gain
    bessel_gain = 1.2255 * 1.9099 #Bessel DC gain

    #I am assuming the samprate (mission or int) in config is the correct one to use here
    samp_hz = [samprate, samprate, 512.0]

    #Generate electronics transfer function
    dcgain = bessel_gain * fixed_gain

    nrec = 0
    for isampl in range(3): #0, 1 = MTM modes, 2=fakeit
        fsamp = samp_hz[isampl]
        for ichan in range(4):
            tauboost = tau[ichan]
            for micro in range(2): #dig fltr on/off (opposite order than apodization)
                for ncompress in range(1, 13):
                    dfhz = fsamp / ncompress / fac_fft_length[0]
                    ampmax = -1.0
                    for k in range(nfreq):
                        freqhz = k * dfhz
                        zbesl = _bessel(freqhz, bes3db)
                        ztboost = _tboost(freqhz, tauboost)
                        zdcblock = _dcblock(freqhz)
                        zdigfil = _digfltr(freqhz, micro, ichan, fsamp)
                        zsmooth = _compress(freqhz, ncompress, fsamp)

                        zanalog = zbesl * ztboost * zdcblock
                        zdigital = zdigfil * zsmooth
                        zxfer = dcgain * zanalog * zdigital

                        #ampl = np.sqrt(np.real(zxfer * zxfer.conjugate()))
                        #ampl = np.sqrt(np.abs(zxfer)**2)
                        #ampmax = max(ampl, ampmax)
                        #if ampl < ampmax/1000.0:
                        #    ztrans[nrec, k] = 100000.0
                        #else:
                        #    ztrans[nrec, k] = - zxfer
                        ztrans[nrec, k] = - zxfer

                    nrec += 1

    return ztrans

def _bessel(fhz, bes3db):
    zfb = 1j*fhz/bes3db
    zfb2 = zfb * zfb
    zfb3 = zfb * zfb2
    zfb4 = zfb * zfb3
    zfb5 = zfb * zfb4

    zbesl = 1.0/(1.0 + 2.4275*zfb + 2.6189*zfb2 + 1.5894*zfb3 + 0.5511*zfb4 + 0.0892*zfb5)
    
    return zbesl

def _tboost(fhz, tau):
    return 1.0 + 2*np.pi*1j*fhz*tau

def _dcblock(fhz):
    zs = 3.2j * 2 * np.pi * fhz
    return (zs/(1.0 + zs))**5 * zs/(2.0 + zs)

def _digfltr(freqhz, micromode, ichan, samplrate):
    zi = -2j*np.pi
    z = np.exp(zi*freqhz/samplrate)
    zdigfil = 1.0
    if micromode == 1: return zdigfil #digital filter off
    if ichan == 1 or ichan == 3:
        #zdigfil = (1.0 + z*z)**2/(8.0 - 12.5*z*z + 5.0*z**4)/8.0
        zdigfil = ((1.0 + 2.0*z*(1.0 + z + z*z) + z**4) / (8.0 - 8.0*z*z + z**4))/8.0
    else:
        zdigfil = (1.0 + z)**2/(8.0 - 11.0*z + 5.0*z*z)/2.0

    return zdigfil

def _compress(fhz, ncompress, samplrate):
    pif = np.pi * fhz / samplrate
    zsmooth = 1.0
    if pif < 0.001: return zsmooth
    ampl = np.sin(pif*ncompress)/np.sin(pif)
    zphase = pif*(1-ncompress)*1j
    zsmooth = ampl / ncompress * np.exp(zphase)

    return zsmooth

