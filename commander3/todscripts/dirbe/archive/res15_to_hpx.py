import numpy as np

# We need to use SUPER_CENPIX and SUPER_PIXNO to convert res 15 pixels to
# lonlat.

def SUPER_CENPIX(nsuper):
    # c is the output unit vector
    # nsuper is theinput pixel number
    c = np.zeros(3)

    IT28 = 2**28

    ix = np.arange(1024)
    iy = np.arange(1024)

    nface = nsuper/IT28

    n = nsuper % IT28
    i = n % 1024
    n = n // 1024
    j = n % 1024
    k = n // 1024
    jx = 1024*ix[k-1] + 32*ix[j-1] + ix[i-1]
    jy = 1024*iy[k-1] + 32*iy[j-1] + iy[i-1]

    x = (jx - 8191.5)/8192. # database coordinates. Range
    y = (jy - 8191.5)/8192. # Covers the square face


    xi, eta = FORWARD_CUBE(x, y)

    xp, yp = INCUBE(xi, eta)
    xi = xi - (xp-x)
    eta = eta - (yp-y)

    xi, eta = FORWARD_CUBE(x, y)
    xi = xi - (xp-x)
    eta = eta-(yp-y)

    c = XYAXIS(nface, xi, eta)

    return c


def SUPER_PIXNO(c):
    # c is the input unit vector
    # output is nsuper, pixel number from 1 to 6*2**28
    IT14 = 16384
    IT28 = 2**28

    # Create lookup table
    ix = np.zeros(128, dtype='int')
    iy = np.zeros(128, dtype='int')
    for i in range(1, 128+1):
        j = i-1
        k = 0
        ip = 1
        if (j == 0):
            ix[j] = k
            iy[j] = 2*k
        else:
            while (j != 0):
                id = j % 2
                j = j // 2
                k = ip*id + k
                ip = ip*4
    nface, x, y = AXISXY(c)

    i = int(IT14 * x)
    j = int(IT14 * y)
    i = min(IT14-1, i) # Need to roundoff near edge of cube
    j = min(IT14-1, j)



    il = i % 128
    ih = i // 128
    jl = j % 128
    jh = j // 128

    print('nface, IT28, ix[il], iy[jl], IT14, ix[ih], iy[jh]')
    print(nface, IT28, ix[il], iy[jl], IT14, ix[ih], iy[jh])

    nsuper = nface*IT28 + ix[il] + iy[jl] + IT14*(ix[ih] + iy[jh])

    return nsuper




def XYAXIS(nface, xi, eta):
   # Converts face number nface (0-5) and xi, eta (-1 - +1) into unit vector c

   c = np.zeros(3)

   xi1 = max(abs(xi), abs(eta))
   eta1 = min(abs(xi), abs(eta))
   norm = 1./(1+xi1**2+eta1**2)**0.5

   if (nface == 0):
       c[2] = norm
       c[0] = -eta*norm
       c[1] = xi*norm
   elif (nface == 1):
       c[0] = norm
       c[2] = eta*norm
       c[1] = xi*norm
   elif (nface == 2):
       c[1] = norm
       c[2] = eta*norm
       c[0] = -xi*norm
   elif (nface == 3):
       c[0] = -norm
       c[2] = eta*norm
       c[1] = -xi*norm
   elif (nface == 4):
       c[1] = -norm
       c[2] = eta*norm
       c[0] = xi*norm
   elif (nface == 5):
       c[2] = -norm
       c[0] = eta*norm
       c[1] = xi*norm
   else:
       print('Not a valid nface')

   return c


def AXISXY(c):
    # Converts unit vector into nface number and x,y in range of 0-1.

    ac3 = abs(c[2])
    ac2 = abs(c[1])
    ac1 = abs(c[0])

    if ac3 > ac3:
        if ac3 > ac1:
            if c[2] > 0:
                # GOTO 100
                nface = 0
                eta = -c[0]/c[2]
                xi = c[1]/c[2]
            else:
                nface = 5
                eta = -c[0]/c[2]
                xi = -c[1]/c[2]
        else:
            if c[1] > 0:
                nface = 1
                xi = c[1]/c[0]
                eta = c[2]/c[0]
            else:
                nface = 3
                eta = -c[2]/c[0]
                xi = c[1]/c[0]
    else:
        if ac2 > ac1:
            if c[1] > 0:
                nface = 2
                eta = c[2]/c[1]
                xi  -c[0]/c[1]
            else:
                nface = 4
                eta = -c[2]/c[1]
                xi = -c[0]/c[1]
        else:
            if c[0] > 0:
                nface = 1
                xi = c[1]/c[0]
                eta = c[2]/c[0]
            else:
                nface = 3
                eta = -c[2]/c[0]
                xi = c[1]/c[0]

    x, y = INCUBE(xi, eta)
    x = (x+1.)/2
    y = (y+1.)/2
    return nface, x, y



def INCUBE(alpha, beta):
    GSTAR = 1.37484847732
    G=-0.13161671474
    M=0.004869491981
    W1=-0.159596235474
    C00=0.141189631152
    C10 = 0.0809701286525
    C01=-0.281528535557
    C11=0.15384112876
    C20=-0.178251207466
    C02=0.106959469314
    D0=0.0759196200467
    D1=-0.0217762490699
    R0=0.577350269




    aa = alpha**2
    bb = beta**2
    a4 = aa**2
    b4 = bb*2
    onmaa = 1.-aa
    onmbb = 1.-bb

    gstar_1 = 1.-GSTAR
    m_g = M-G
    c_comb = C00+C11*aa+bb

    x = alpha*(
            GSTAR +
            aa * gstar_1 +
            onmaa * (bb * (G + m_g*aa +
                           onmbb * (c_comb + C10*aa + C01*bb +
                                    C20*a4 + C02*b4)) +
                     aa * (W1 - onmaa*(D0 + D1*aa))))

    y = beta*(
            GSTAR +
            bb * gstar_1 +
            onmbb * (aa * (G+m_g*bb +
                           onmaa * (c_comb + C10*bb + C01*aa + 
                                    C20*b4 + C02*a4)) +
                     bb * (W1 - onmbb*(D0 + D1*bb))))
            
    return x,y


def FORWARD_CUBE(x, y):
    # x,y in -1 to +1 are database coordinates
    # xi, eta are in range -1 to + 1, tangent plane coordinates

    # These numbers come from FCFIT.FOR


    P = np.array([
        -0.27292696,    -0.07629969,    -0.02819452,    -0.22797056,
        -0.01471565,    0.27058160,     0.54852384,     0.48051509,
        -0.56800938,    -0.60441560,    -0.62930065,    -1.74114454,
        0.30803317,     1.50880086,     0.93412077,     0.25795794,
        1.71547508,     0.98938102,     -0.93678576,    -1.41601920,
        -0.63915306,    0.02584375,     -0.53022337,    -0.83180469,
        0.08693841,     0.33887446,     0.52032238,     0.14381585])

    xx = x**2
    yy = y**2
    xi = x*(1.+(1.-xx)*(
        P[0] + xx*(P[1] + xx*(P[3]+xx*(P[6]+xx*(P[10]+xx*(P[15]+xx*P[21]))))) +
        yy * (P[2] + xx*(P[4] + xx*(P[7] + xx*(P[11] + xx*(P[16]+xx*P[22])))) +
        yy * (P[5] + xx*(P[8] + xx*(P[12] + xx*(P[17] + xx*P[23]))) +
        yy * (P[9] + xx*(P[13] + xx*(P[18] + xx*P[24])) +
        yy * (P[14] +xx*(P[19] + xx*P[25]) +
        yy * (P[20] + xx*P[26] + yy*P[27])))) )))

    eta = y*(1.+(1.-yy)*(
        P[0] + yy*(P[1] + yy*(P[3]+yy*(P[6]+yy*(P[10]+yy*(P[15]+yy*P[21]))))) +
        xx * (P[2] + yy*(P[4] + yy*(P[7] + yy*(P[11] + yy*(P[16]+yy*P[22])))) +
        xx * (P[5] + yy*(P[8] + yy*(P[12] + yy*(P[17] + yy*P[23]))) +
        xx * (P[9] + yy*(P[13] + yy*(P[18] + yy*P[24])) +
        xx * (P[14] +yy*(P[19] + yy*P[25]) +
        xx * (P[20] + yy*P[26] + xx*P[27])))) )))

    return xi, eta


if __name__ == '__main__':
    c = np.array([1,1,1])/3**0.5
    pix = SUPER_PIXNO(c)
    print(pix)
    print(SUPER_CENPIX(pix))
