import numpy
import scipy.linalg
import scipy.special

def fmt(t, m, low=None):
#             _ 1           2
#            /     2 m  -t u
# F (t)  =   |    u    e      du,
#  m        _/  0
#
    # fmt1 is alaways more accurate than fmt2 and ~3 times slower than fmt2
    if (t < m + 1.5):
        return fmt1(t, m, low)
    else:
        return fmt2(t, m, low)

def fmt1(t, m, low=None):
#
# F[m] = int u^{2m} e^{-t u^2} du
#      = 1/(2m+1) int e^{-t u^2} d u^{2m+1}
#      = 1/(2m+1) [e^{-t u^2} u^{2m+1}]_0^1 + (2t)/(2m+1) int u^{2m+2} e^{-t u^2} du
#      = 1/(2m+1) e^{-t} + (2t)/(2m+1) F[m+1]
#      = 1/(2m+1) e^{-t} + (2t)/(2m+1)(2m+3) e^{-t} + (2t)^2/(2m+1)(2m+3) F[m+2]
#
    f = numpy.zeros(m+1)
    b = m + 0.5
    e = .5 * numpy.exp(-t)
    x = e
    s = e
    bi = b + 1
    while x > 1e-18:
        x *= t / bi
        s += x
        bi += 1
    f[m] = s / b
    for i in reversed(range(m)):
        b -= 1
        f[i] = (e + t * f[i+1]) / b
    return f

#
# F[m] = int u^{2m} e^{-t u^2} du
#      = -1/2t int u^{2m-1} d e^{-t u^2}
#      = -1/2t [e^{-t u^2} * u^{2m-1}]_0^1 + (2m-1)/2t int u^{2m-2} e^{-t u^2} du
#      = 1/2t (-e^{-t} + (2m-1) F[m-1])
#
def fmt2(t, m, low=None):
    f = numpy.zeros(m+1)
    tt = numpy.sqrt(t)
    f[0] = numpy.pi**.5/2. / tt * scipy.special.erf(tt)
    e = numpy.exp(-t)
    b = .5 / t
    for i in range(m):
        f[i+1] = b * ((2*i+1) * f[i] - e)
    return f

def fmt_erfc(t, m, low=0):
    if m > 3:
        turnover = 4
    else:
        turnover = m + 0.5
    turnover = m + 1.5
    if t < turnover:
        return fmt1_erfc(t, m, low)
    else:
        return fmt2_erfc(t, m, low)

def fmt1_erfc(t, m, low=0):
#
# F[m] = int_s^1 u^{2m} e^{-t u^2} du
#      = 1/(2m+1) int e^{-t u^2} d u^{2m+1}
#      = 1/(2m+1) [e^{-t u^2} u^{2m+1}]_s^1 + (2t)/(2m+1) int u^{2m+2} e^{-t u^2} du
#      = 1/(m+.5) (.5*e^{-t} - .5*e^{-t s^2} s^{2m+1}) + t F[m+1])
#
    f = numpy.zeros(m+1)
    b = m + 0.5
    e = .5 * numpy.exp(-t)
    e1 = .5 * numpy.exp(-t * low**2) * low**(2*m+1)
    x = e
    x1 = e1
    s = e - e1
    bi = b + 1
    while x > 1e-18:
        x *= t / bi
        x1 *= low**2 * t / bi
        s += x - x1
        bi += 1
    f[m] = s / b
    for i in reversed(range(m)):
        b -= 1
        e1 /= low**2
        f[i] = (e - e1 + t * f[i+1]) / b
    return f

#
# F[m] = int_s^1 u^{2m} e^{-t u^2} du
#      = -1/2t int u^{2m-1} d e^{-t u^2}
#      = -1/2t [e^{-t u^2} * u^{2m-1}]_s^1 + (2m-1)/2t int u^{2m-2} e^{-t u^2} du
#      = 1/2t (-e^{-t} + e{-t s^2} s^{2m-1} + (2m-1)/ F[m-1])
# F[0] = int_s^1 e^{-t u^2} du
#
def fmt2_erfc(t, m, low=0):
    f = numpy.zeros(m+1)
    tt = numpy.sqrt(t)
    f[0] = numpy.pi**.5/2. / tt * (scipy.special.erf(tt) - scipy.special.erf(low*tt))
    e = numpy.exp(-t)
    e1 = numpy.exp(-t*low**2) * low
    b = .5 / t
    for i in range(m):
        f[i+1] = b * ((2*i+1) * f[i] - e + e1)
        e1 *= low**2
    return f
########################################


def poly_value1(a, order, x):
    p = a[order]
    for i in range(1, order+1):
        p = p * x + a[order-i]
    return p

def R_dnode(a, rt):
    #TODO: a several step bisection search followed by Newton-Raphson refinement. 
    order = a.size - 1
    accrt = 1e-15

    x1init = 0
    p1init = a[0]
    for m in range(order):
        x0 = x1init
        p0 = p1init
        x1init = rt[m]
        p1init = poly_value1(a, order, x1init)
        if (p0 * p1init > 0):
            raise RuntimeError

        if (x0 <= x1init):
            x1 = x1init
            p1 = p1init
        else:
            x1 = x0
            p1 = p0
            x0 = x1init
            p0 = p1init

        if (p0 == 0):
            rt[m] = x0
            continue
        elif (p1 == 0):
            rt[m] = x1
            continue
        else:
            xi = x0 + (x0 - x1) / (p1 - p0) * p0

        n = 0
        while (x1 > accrt+x0) or (x0 > x1+accrt):
            n += 1
            if (n > 600):
                raise RuntimeError

            pi = poly_value1(a, order, xi)
            if (pi == 0):
                break
            elif (p0 * pi <= 0):
                x1 = xi
                p1 = pi
                xi = x0 * .25 + xi * .75
            else:
                x0 = xi
                p0 = pi
                xi = xi * .75 + x1 * .25

            pi = poly_value1(a, order, xi)
            if (pi == 0):
                break
            elif (p0 * pi <= 0):
                x1 = xi
                p1 = pi
            else:
                x0 = xi
                p0 = pi

            xi = x0 + (x0 - x1) / (p1 - p0) * p0
        rt[m] = xi
    return rt

def R_dsmit(s, n):
    cs = numpy.zeros_like(s)
    for j in range(n):
        fac = s[j,j]
        v = numpy.zeros(j)
        for k in range(j):
            dot = cs[:j+1,k].dot(s[:j+1,j])
            v -= dot * cs[:j,k]
            fac = fac - dot * dot

        fac = fac ** -.5
        cs[j,j] = fac
        cs[:j,j] = fac * v[:j]
    return cs

def rys_root(nroots, x, ff=None):
    if ff is None:
        ff = fmt(x, nroots*2)

    nroots1 = nroots + 1
    s = numpy.zeros((nroots1, nroots1))
    for j in range(nroots1):
        for i in range(nroots1):
            s[i, j] = ff[i + j]

    cs = R_dsmit(s, nroots1)
    rt = numpy.ones(nroots1)

    if nroots == 1:
        rt[0] = ff[1] / ff[0]
    else:
        dum = numpy.sqrt(cs[1,2] * cs[1,2] - 4 * cs[0,2] * cs[2,2])
        rt[0] = .5 * (-cs[1,2] - dum) / cs[2,2]
        rt[1] = .5 * (-cs[1,2] + dum) / cs[2,2]

    for k in range(2, nroots):
        R_dnode(cs[:k+2,k+1], rt)

    roots = numpy.zeros(nroots)
    weights = numpy.zeros(nroots)
    for i in range(nroots):
        root = rt[i]
        dum = 1 / ff[0]
        for j in range(1, nroots):
            poly = poly_value1(cs[:j+1,j], j, root)
            dum += poly * poly
        weights[i] = 1 / dum
        roots[i] = root / (1 - root)
    return roots, weights
