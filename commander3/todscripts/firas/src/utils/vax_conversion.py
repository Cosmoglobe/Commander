"""
Code from http://yaadc.blogspot.com/2013/01/vax-ffloat-and-dfloat-to-ieee-tfloat.html
Also see
https://idlastro.gsfc.nasa.gov/ftp/obsolete/conv_vax_unix.pro
"""

from __future__ import print_function

import struct
import numpy as np

SIGN_BIT = 0x80000000

VAX_D_SIGN_BIT = SIGN_BIT
VAX_D_EXPONENT_MASK = 0x7F800000
VAX_D_EXPONENT_SIZE = 8
VAX_D_EXPONENT_BIAS = (1 << ( VAX_D_EXPONENT_SIZE - 1))
VAX_D_MANTISSA_MASK = 0x007FFFFF
VAX_D_MANTISSA_SIZE = 23
VAX_D_HIDDEN_BIT    = (1 << VAX_D_MANTISSA_SIZE)

VAX_F_MANTISSA_MASK = 0x007FFFFF
VAX_F_MANTISSA_SIZE = 23
VAX_F_HIDDEN_BIT = (1 << VAX_F_MANTISSA_SIZE)
VAX_F_EXPONENT_SIZE = 8
VAX_F_EXPONENT_BIAS = (1 << (VAX_F_EXPONENT_SIZE - 1))
VAX_F_EXPONENT_MASK = 0x7F800000

IEEE_T_MANTISSA_SIZE = 20
IEEE_T_EXPONENT_SIZE = 11
IEEE_T_EXPONENT_BIAS = ((1 << (IEEE_T_EXPONENT_SIZE - 1)) - 1)

IEEE_S_EXPONENT_SIZE = 8
IEEE_S_EXPONENT_BIAS = ((1 << (IEEE_S_EXPONENT_SIZE - 1)) - 1)
IEEE_S_MANTISSA_SIZE = 23

def bitandnot(v1, v2):
    return (((v1 // 65536) & (0xffff - (v2 // 65536))) * 65536
          + ((v1 % 65536) & (0xffff - (v2 % 65536))))

def vd8tod(stream):
    EXPONENT_RIGHT_SHIFT = (VAX_D_MANTISSA_SIZE - IEEE_T_MANTISSA_SIZE)
    EXPONENT_ADJUSTMENT = (1 + VAX_D_EXPONENT_BIAS - IEEE_T_EXPONENT_BIAS)
    IN_PLACE_EXPONENT_ADJUSTMENT = (EXPONENT_ADJUSTMENT << IEEE_T_MANTISSA_SIZE)

    v1a, v1b, v2a, v2b = struct.unpack('HHHH', stream)

    v1 = v1b + 65536 * v1a
    v2 = v2b + 65536 * v2a

    if (v1 & VAX_D_EXPONENT_MASK) == 0:
        if (v1 & SIGN_BIT) == SIGN_BIT:
            raise
        res = '\x00' * 8
    else:
        i1 = (((v1 & SIGN_BIT) | (bitandnot(v1, SIGN_BIT) >> EXPONENT_RIGHT_SHIFT)) - IN_PLACE_EXPONENT_ADJUSTMENT )
        i2 = (((v1 << (32 - EXPONENT_RIGHT_SHIFT)) | (v2 >> EXPONENT_RIGHT_SHIFT))) & 0xffff

        #res = struct.pack('LL', i2, i1);
        res = struct.pack('II', i2, i1);

    if type(res) == str:
        res = res.encode('utf-8')

    return struct.unpack('d', res)[0];

def vr4tof(stream):
    EXPONENT_ADJUSTMENT = (1 + VAX_F_EXPONENT_BIAS - IEEE_S_EXPONENT_BIAS)
    IN_PLACE_EXPONENT_ADJUSTMENT = (EXPONENT_ADJUSTMENT << IEEE_S_MANTISSA_SIZE)

    v1a, v1b = struct.unpack('HH', stream)
    v1 = v1b + 65536 * v1a

    e = v1 & VAX_F_EXPONENT_MASK

    #if e == 0:
    #    if v1 & SIGN_BIT == SIGN_BIT:
    #        raise
    #    res = '\x00' * 4
    #else:
    e >>= VAX_F_MANTISSA_SIZE
    e -= EXPONENT_ADJUSTMENT

    if e > 0:
        res = v1 - IN_PLACE_EXPONENT_ADJUSTMENT
    else:
        res = (v1 & SIGN_BIT) | ((VAX_F_HIDDEN_BIT | (v1 & VAX_F_MANTISSA_MASK)) >> (1 - e))

    #res = struct.pack('L', res)
    res = struct.pack('I', res)

    return struct.unpack('f', res)[0]

np_vd8tod = np.vectorize(vd8tod)
np_vr4tof = np.vectorize(vr4tof)

if __name__ == "__main__":
    stream = b'\x40\xc1\x00\x00\x00\x00\x00\x00' # -3.0
    print(vd8tod(stream))
    stream = b'\x9b\x41\x33\x33\x00\x00\x00\x00' # approx 4.85
    print(vd8tod(stream))
    print(vr4tof(stream[:4]))
