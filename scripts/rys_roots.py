import mpmath
DECIMALS = 100
mpmath.mp.dps = DECIMALS
import numpy as np
from find_polyroots import find_polyroots

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
    half = mpmath.mpf('.5')
    b = m + half
    e = half * mpmath.exp(-t)
    x = e
    s = e
    bi = b + 1
    while x > .1**DECIMALS:
        x *= t / bi
        s += x
        bi += 1
    f = s / b
    out = [f]
    for i in range(m):
        b -= 1
        f = (e + t * f) / b
        out.append(f)
    return np.array(out)[::-1]

#
# F[m] = int u^{2m} e^{-t u^2} du
#      = -1/2t int u^{2m-1} d e^{-t u^2}
#      = -1/2t [e^{-t u^2} * u^{2m-1}]_0^1 + (2m-1)/2t int u^{2m-2} e^{-t u^2} du
#      = 1/2t (-e^{-t} + (2m-1) F[m-1])
#
def fmt2(t, m, low=None):
    half = mpmath.mpf('.5')
    tt = mpmath.sqrt(t)
    f = mpmath.pi**half/2 / tt * mpmath.erf(tt)
    e = mpmath.exp(-t)
    b = half / t
    out = [f]
    for i in range(m):
        f = b * ((2*i+1) * f - e)
        out.append(f)
    return np.array(out)

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
#             _ 1           2
#            /     2 m  -t u
# F (t)  =   |    u    e      du,
#  m        _/  s
#
#
# F[m] = int_s^1 u^{2m} e^{-t u^2} du
#      = 1/(2m+1) int e^{-t u^2} d u^{2m+1}
#      = 1/(2m+1) [e^{-t u^2} u^{2m+1}]_s^1 + (2t)/(2m+1) int u^{2m+2} e^{-t u^2} du
#      = 1/(m+.5) (.5*e^{-t} - .5*e^{-t s^2} s^{2m+1}) + t F[m+1])
#
    half = mpmath.mpf('.5')
    low = mpmath.mpf(low)
    b = m + half
    e = half * mpmath.exp(-t)
    e1 = half * mpmath.exp(-t * low*low) * mpmath.power(low, 2*m+1)
    x = e
    x1 = e1
    s = e - e1
    bi = b
    div = mpmath.mpf(1)
    delta = s
    while abs(delta) > .1**DECIMALS:
        bi += 1
        div *= t / bi
        x1 *= low*low
        delta = (x - x1) * div
        s += delta
    f = s / b
    out = [f]
    for i in range(m):
        b -= 1
        e1 /= low*low
        f = (e - e1 + t * f) / b
        out.append(f)
    return np.array(out)[::-1]

#
# F[m] = int_s^1 u^{2m} e^{-t u^2} du
#      = -1/2t int u^{2m-1} d e^{-t u^2}
#      = -1/2t [e^{-t u^2} * u^{2m-1}]_s^1 + (2m-1)/2t int u^{2m-2} e^{-t u^2} du
#      = 1/2t (-e^{-t} + e{-t s^2} s^{2m-1} + (2m-1)/ F[m-1])
# F[0] = int_s^1 e^{-t u^2} du
#
def fmt2_erfc(t, m, low=0):
    half = mpmath.mpf('.5')
    tt = mpmath.sqrt(t)
    low = mpmath.mpf(low)
    low2 = low * low
    f = mpmath.sqrt(mpmath.pi)/2 / tt * (mpmath.erf(tt) - mpmath.erf(low*tt))
    e = mpmath.exp(-t)
    e1 = mpmath.exp(-t*low2) * low
    b = half / t
    out = [f]
    #.5/t*(1+mpmath.exp(-t*low2)*low/f)
    for i in range(m):
        f = b * ((2*i+1) * f - e + e1)
        e1 *= low2
        out.append(f)
    return np.array(out)
########################################


def poly_value1(a, order, x):
    p = a[order]
    for i in range(1, order+1):
        p = p * x + a[order-i]
    return p

def zeros(shape):
    if isinstance(shape, int):
        size = shape
    else:
        size = np.prod(shape)
    return np.array([mpmath.mpf(0)] * size).reshape(shape)

def R_dsmit(s, n):
    cs = zeros(s.shape)
    for j in range(n):
        fac = s[j,j]
        v = zeros(j)
        for k in range(j):
            dot = cs[k,:j+1].dot(s[j,:j+1])
            v -= dot * cs[k,:j]
            fac = fac - dot * dot

        if fac <= 0:
            raise RuntimeError(f'R_dsmit fac = {fac}')
        fac = 1 / mpmath.sqrt(fac)
        cs[j,j] = fac
        cs[j,:j] = fac * v[:j]
    return cs

def rys_roots_weights_partial(nroots, x, low=None):
    if low is None:
        ff = fmt(x, nroots*2)
        if ff[0] < .1**(DECIMALS/4):
            return zeros(nroots), zeros(nroots)
    else:
        if x * low**2 > DECIMALS*.7:
            return zeros(nroots), zeros(nroots)
        ff = fmt_erfc(x, nroots*2, low)

    nroots1 = nroots + 1
    s = zeros((nroots1, nroots1))
    for j in range(nroots1):
        for i in range(nroots1):
            s[i, j] = ff[i + j]

    cs = R_dsmit(s, nroots1)
    rt = find_polyroots(cs, nroots)

    weights = zeros(nroots)
    for i in range(nroots):
        root = rt[i]
        dum = 1 / ff[0]
        for j in range(1, nroots):
            poly = poly_value1(cs[j,:j+1], j, root)
            dum += poly * poly
        weights[i] = 1 / dum
    return rt[:nroots], weights

def rys_roots_weights(nroots, x, low=None):
    rt, weights = rys_roots_weights_partial(nroots, x, low)
    roots = rt / (1 - rt)
    return roots, weights

# when erf(x**.5) ~= 1
def polyfit_large_x_limits(nroots, x, low=None):
    # polynomial fits for large X
    # roots = rx / (x - rx)
    # weights = mpmath.sqrt(mpmath.pi/4/x) * rat
    #assert x > 1e9
    assert low is None # No large x limits for sr_rys_roots
    r, w = rys_roots_weights_partial(nroots, x, low)
    rx = r * x
    tt = x**.5
    w0 = mpmath.sqrt(mpmath.pi/4/x)
    if low is not None:
        w0 *= mpmath.erfc(low * tt) - mpmath.erfc(tt)
    rat = w / w0 # rat.sum() == 1
    rat[0] = 1 - rat[1:].sum()
    return rx, rat

def polyfit_small_x_limits(nroots, x, low=None):
    # polynomial fits for large X
    # roots = pr1 * x + pr0;
    # weights = pw1 * w + pw0;
    assert x < 1e-12
    x0 = x
    x1 = x / 2**10
    r0, w0 = rys_roots_weights(nroots, x0, low)
    r1, w1 = rys_roots_weights(nroots, x1, low)
    pr1 = (r1 - r0) / (x1 - x0)
    pr0 = r0 - pr1 * x0
    pw1 = (w1 - w0) / (x1 - x0)
    pw0 = w0 - pw1 * x0
    return pr0, pr1, pw0, pw1
