#!/usr/bin/env python

import mpmath
DECIMALS = 24
mpmath.mp.dps = DECIMALS

def binomial(n, m):
    if m < 0 or m > n:
        return 0
    else:
        return mpmath.factorial(n) / mpmath.factorial(m) / mpmath.factorial(n-m)

def cart2sph_pure(lx, ly, lz, m):
    '''
    The coefficient that transforms normalized Cartesian spherical harmonic
    GTOs to normalized pure spherical harmonic GTOs.

         gY(l,m) = \sum_{lx+ly+lz=l} gO(lx,ly,lz) * c(l,m,lx,ly,lz)

    gY(l,m) is a pure spherical harmonic GTO and gO(lx,ly,lz) is a Cartesian
    spherical GTO with Condon-Shortley phase convention.

    Ref: 
      H. B. Schlegel and M. J. Frisch, Int. J. Quant. Chem., 54(1955), 83-87. Eq. (15)
    '''
    if (lx+ly+m) % 2:
        return 0

    l = lx + ly + lz
    j = (lx + ly - abs(m)) // 2
    num = mpmath.factorial(lx*2) * mpmath.factorial(ly*2) * mpmath.factorial(lz*2) \
        * mpmath.factorial(l) * mpmath.factorial(l-abs(m))
    div = mpmath.factorial(l*2) * mpmath.factorial(l+abs(m)) \
        * mpmath.factorial(lx) * mpmath.factorial(ly) * mpmath.factorial(lz)
    c0 = mpmath.sqrt(num/div) * (.5**l) / mpmath.factorial(l)

    cp = 0
    for i in range(max(0,j), (l-abs(m))//2+1):
        cp0 = binomial(l, i) * binomial(i, j) \
            * mpmath.factorial(2*l-2*i) / mpmath.factorial(l-abs(m)-2*i)
        if i%2:
            cp -= cp0
        else:
            cp += cp0
    c0 = c0 * cp

    c1 = 0
    for k in range(max(0,(lx-abs(m)+1)//2), min(j,lx//2)+1):
        cp0 = binomial(j, k) * binomial(abs(m), lx-2*k)
        m4 = (abs(m)-lx+2*k) % 4
        if m4 == 0:
            c1 += cp0
        elif m4 == 1:
            c1 += cp0*1j
        elif m4 == 2:
            c1 -= cp0
        else:
            c1 -= cp0*1j

# (-)^m introduced by Condon-Shortley phase convention in Legendre function.
# IJQC_54_83 does not have this phase parameter
    if m >= 0:
        if m%2:
            c =-c0 * c1
        else:
            c = c0 * c1
    else:
# m<0, Y(l,-m) = (-)^m Y(l,m)^*
        c = c0 * c1.conjugate()

    return c


def xyz2sph_real(lx, ly, lz, m):
    '''
    Factor of xyz component for normalized real spherical harmonic functions

         r^l*Y(l,m) = \sum_{lx+ly+lz=l} x^lx*y^lyz^lz * c(l,m,lx,ly,lz)

    Y(l,m) is a real spherical harmonic function with Condon-Shortley phase
    convention
  
    Ref: 
      H. B. Schlegel and M. J. Frisch, Int. J. Quant. Chem., 54(1995), 83-87.
    '''

    if (lx+ly+m) % 2:
        return 0

    l = lx + ly + lz
    j = (lx + ly - abs(m)) // 2
    num = (2*l+1) * mpmath.factorial(l-abs(m))
    div = mpmath.factorial(l+abs(m))
    c0 = mpmath.sqrt(num/(div*4*mpmath.pi)) * (.5**l) / mpmath.factorial(l)

    cp = 0
    for i in range(max(0,j), (l-abs(m))//2+1):
        cp0 = binomial(l,i) * binomial(i,j) \
            * mpmath.factorial(2*l-2*i) / mpmath.factorial(l-abs(m)-2*i)
        if i%2:
            cp -= cp0
        else:
            cp += cp0
    c0 = c0 * cp

    cp = 0
    if m >= 0:
        for k in range(max(0,(lx-abs(m)+1)//2), min(j,lx//2)+1):
            cp0 = binomial(j,k) * binomial(abs(m), lx-2*k)
            r = (abs(m)-lx+2*k) % 4
            if r == 0:
                cp = cp + cp0
            elif r == 2:
                cp = cp - cp0
        if m == 0:
            c = c0 * cp
        else:
            c = mpmath.sqrt(2) * c0 * cp
    else:
        for k in range(max(0,(lx-abs(m)+1)//2), min(j,lx//2)+1):
            cp0 = binomial(j,k) * binomial(abs(m),lx-2*k)
            r = (abs(m)-lx+2*k) % 4
            if r == 1:
                cp = cp - cp0
            elif r == 3:
                cp = cp + cp0
        c = -mpmath.sqrt(2) * c0 * cp

    return c

def xyz2sph_pure(lx, ly, lz, m):
    '''
    Factor of xyz component for normalized pure spherical harmonic functions

         r^l*Y(l,m) = \sum_{lx+ly+lz=l} x^lx*y^lyz^lz * c(l,m,lx,ly,lz)
    Y(l,m) is a pure spherical harmonic function with Condon-Shortley phase
    convention
  
    Ref: 
      H. B. Schlegel and M. J. Frisch, Int. J. Quant. Chem., 54(1995), 83-87.  Eq. (7)
    '''

    if (lx+ly+m)%2:
        return 0

    l = lx + ly + lz
    j = (lx + ly - abs(m)) // 2
    num = (2*l+1) * mpmath.factorial(l-abs(m))
    div = mpmath.factorial(l+abs(m))
    c0 = mpmath.sqrt(num/(div*4*mpmath.pi)) * (.5**l) / mpmath.factorial(l)

    cp = 0
    for i in range(max(0,j), (l-abs(m))//2+1):
        cp0 = binomial(l,i) * binomial(i,j) \
            * mpmath.factorial(2*l-2*i) / mpmath.factorial(l-abs(m)-2*i)
        if i%2:
            cp -= cp0
        else:
            cp += cp0
    c0 = c0 * cp

    c1 = 0
    for k in range(max(0,(lx-abs(m)+1)//2), min(j,lx//2)+1):
        cp0 = binomial(j,k) * binomial(abs(m),lx-2*k)
        r = (abs(m)-lx+2*k) % 4
        if r == 0:
            c1 += cp0
        elif r == 1:
            c1 += cp0*1j
        elif r == 2:
            c1 -= cp0
        else:
            c1 -= cp0*1j

# Condon-Shortley phase convention
    if m >= 0:
        if m%2:
            c =-c0 * c1
        else:
            c = c0 * c1
    else:
        c = c0 * c1.conjugate()

    return c


if __name__ == '__main__':
    l = 6
    for m in range(-l,l+1):
        for lx in range(l, -1, -1):
            for ly in range(l-lx, -1, -1):
                lz = l - lx - ly
                print(xyz2sph_real(lx, ly, lz, m))
