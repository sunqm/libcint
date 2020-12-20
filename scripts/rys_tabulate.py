import os
import mpmath
import numpy as np
from rys_roots import DECIMALS, rys_roots_weights_partial, rys_roots_weights
mpmath.mp.dps = DECIMALS

def chebyshev_roots(n):
    return [mpmath.cos(mpmath.pi * (k + .5) / n) for k in range(n)]

def clenshaw_points(n):
    out = []
    for j in range(n):
        c = [mpmath.cos(mpmath.pi * j * (k + 0.5) / n) for k in range(n)]
        out.append(c)
    return out

ngrids = 14
chebrt = np.array(chebyshev_roots(ngrids))
cs = np.array(clenshaw_points(ngrids))

def get_tabulate_points(tbase, Tp):
    tabulate_points = []
    for r in chebrt:
        if tbase == 0:
            x = (r/2 + Tp) ** 2
        else:
            x = 3 ** (r/2 + Tp - 1)
        tabulate_points.append(x)
    return tabulate_points

def tabulate_erf(nroots, tbase):
    Tmin = mpmath.mpf(tbase)
    Tmax = Tmin + 1
    Tp = (Tmin + Tmax) / 2

    tabulate_points = get_tabulate_points(tbase, Tp)

    fac = mpmath.mpf(2) / ngrids

    rs = []
    ws = []
    for x in tabulate_points:
        r, w = rys_roots_weights_partial(nroots, x)
        rs.append(r)
        ws.append(w)
    rs = np.array(rs)
    ws = np.array(ws)

    # einsum('kn,jk->nj', rs, cs)
    tab_rs = rs.T.dot(cs.T) * fac
    tab_ws = ws.T.dot(cs.T) * fac
    return tab_rs, tab_ws

def tabulate_erfc(nroots, tbase):
    Tmin = mpmath.mpf(tbase)
    Tmax = Tmin + 1
    Tp = (Tmin + Tmax) / 2

    tabulate_points = get_tabulate_points(tbase, Tp)
    boundary_points = [r/2 + 0.5 for r in chebrt]

    fac = mpmath.mpf(2) / ngrids

    rs = []
    ws = []
    for x in tabulate_points:
        for low in boundary_points:
            r, w = rys_roots_weights_partial(nroots, x, low)
            rs.append(r)
            ws.append(w)
    rs = np.array(rs).reshape(ngrids, ngrids, nroots)
    ws = np.array(ws).reshape(ngrids, ngrids, nroots)

    # einsum('lkn,il,jk->nji', rs, cs)
    tab_rs = np.array([cs.dot(rs[:,:,ir]).dot(cs.T).T * fac**2 for ir in range(nroots)])
    tab_ws = np.array([cs.dot(ws[:,:,ir]).dot(cs.T).T * fac**2 for ir in range(nroots)])
    return tab_rs, tab_ws

def clenshaw_d1(x, u, nroots):
    rr = []
    u2 = u * 2
    for i in range(nroots):
        xi = x[i]
        d0 = 0
        d1 = xi[ngrids-1]
        for k in reversed(range(1, ngrids-1)):
            d1, d0 = u2 * d1 - d0 + xi[k], d1
        rr.append(u * d1 - d0 + xi[0] * 0.5)
    return rr


def polyfit_erf(nroots, x):
    t = x
    if t > 19682.99:
        t = 19682.99

    if t > 1.0:
        tt = mpmath.log(t) / mpmath.log(3) + 1.0  # log3(t) + 1
    else:
        tt = mpmath.sqrt(t)

    it = int(tt)
    if isinstance(x, float):
        tt = float(tt) - it
    else:
        tt = tt - it
    tt = 2.0 * tt - 1.0

    tab_rs, tab_ws = tabulate_erf(nroots, it)
    rr = clenshaw_d1(tab_rs.astype(float), tt, nroots)
    ww = clenshaw_d1(tab_ws.astype(float), tt, nroots)
    rr = [r/(1-r) for r in rr]
    return rr, ww

def polyfit_erfc(nroots, x, low):
    t = x
    if t > 19682.99:
        t = 19682.99

    if t > 1.0:
        tt = mpmath.log(t) / mpmath.log(3) + 1.0  # log3(t) + 1
    else:
        tt = mpmath.sqrt(t)

    it = int(tt)
    tt = tt - it
    tt = 2.0 * tt - 1.0  # map [0, 1] to [-1, 1]
    u = low * 2 - 1      # map [0, 1] to [-1, 1]

    tab_rs, tab_ws = tabulate_erfc(nroots, it)
    im = clenshaw_d1(tab_rs.astype(float), u, nroots)
    rr = clenshaw_d1(im, tt, nroots)
    rr = [r/(1-r) for r in rr]

    im = clenshaw_d1(tab_ws.astype(float), u, nroots)
    ww = clenshaw_d1(im, tt, nroots)
    if x * low**2 < DECIMALS*.7:
        factor = mpmath.exp(-x * low**2)
        ww = [w * factor for w in ww]
    return rr, ww

def generate_table(path):
    with open(os.path.join(path, 'roots_x.dat'), 'w') as fx, open(os.path.join(path, 'roots_w.dat'), 'w') as fw:
        for nroots in range(1, 14):
            fx.write('/* root=%d */\n' % nroots)
            fw.write('/* root=%d */\n' % nroots)
            for tbase in range(10):
                print('root %d  tbase %d' % (nroots, tbase))
                tab_rs, tab_ws = tabulate_erfc(nroots, tbase)
                fmt_r = [('    %.17e' % x)[-26:] for x in tab_rs.ravel()]
                fmt_w = [('    %.17e' % x)[-26:] for x in tab_ws.ravel()]
                fx.write(',\n'.join(fmt_r))
                fx.write(',\n')
                fw.write(',\n'.join(fmt_w))
                fw.write(',\n')

if __name__ == '__main__':
    generate_table('.')
