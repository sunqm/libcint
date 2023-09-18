import os
import pickle
from concurrent.futures import ProcessPoolExecutor
import mpmath
import numpy as np
from rys_roots import DECIMALS, rys_roots_weights_partial, rys_roots_weights
mpmath.mp.dps = DECIMALS

def chebyshev_roots(n):
    '''np.cos(np.pi / n * (np.arange(n) + .5))'''
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
TBASE = np.append(np.arange(0, 39, 2.5), np.arange(40, 104, 4))
print(TBASE)

def get_cheb_t_points(tbase):
    if tbase < 1000:
        interval = TBASE[tbase+1] - TBASE[tbase]
        return (chebrt+1) * interval/2 + TBASE[tbase]
    else:
        b1 = float(TBASE[tbase+1])
        b0 = float(TBASE[tbase])
        interval = mpmath.log(b1) - mpmath.log(b0)
        return np.array([mpmath.exp((x+1)*interval/2) * b0 for x in chebrt])

def cheb_t_interval(x, tbase):
    '''map to interval [-1, 1]'''
    if tbase < 1000:
        interval = TBASE[tbase+1] - TBASE[tbase]
        tt = (x - TBASE[tbase]) * 2 / interval - 1
    else:
        b1 = float(TBASE[tbase+1])
        b0 = float(TBASE[tbase])
        interval = mpmath.log(b1) - mpmath.log(b0)
        tt = (mpmath.log(x) - mpmath.log(b0)) * 2 / interval - 1
    return tt

def find_tbase(x):
    return np.searchsorted(TBASE, x) - 1

def tabulate_erf(nroots, tbase):
    tabulate_points = get_cheb_t_points(tbase)

    fac = mpmath.mpf(2) / ngrids

    rs = []
    ws = []
    for x in tabulate_points:
        r, w = rys_roots_weights(nroots, x)
        rs.append(r)
        ws.append(w)
    rs = np.array(rs)
    ws = np.array(ws)

    # einsum('kn,jk->nj', rs, cs)
    tab_rs = rs.T.dot(cs.T) * fac
    tab_ws = ws.T.dot(cs.T) * fac
    return tab_rs, tab_ws

def clenshaw_d1(x, u, nroots):
    rr = []
    u2 = u * 2
    for i in range(nroots):
        xi = x[i]
        ng = len(xi)
        d0 = 0
        d1 = xi[ng-1]
        for k in reversed(range(1, ng-1)):
            d1, d0 = u2 * d1 - d0 + xi[k], d1
        rr.append(u * d1 - d0 + xi[0] * 0.5)
    return np.array(rr)


def polyfit_erf(nroots, x):
    tbase = find_tbase(x)
    tt = cheb_t_interval(x, tbase)
    tab_rs, tab_ws = tabulate_erf(nroots, tbase)
    rr = clenshaw_d1(tab_rs.astype(float), tt, nroots)
    ww = clenshaw_d1(tab_ws.astype(float), tt, nroots)
    #rr = rr/(1-rr)
    return rr, ww

def pkl2table(prefix, pklfile):
    with open(pklfile, 'rb') as f:
        TBASE, rys_tab = pickle.load(f)
    TBASE = TBASE.round(6)
    nt = len(TBASE) - 1
    #nt = find_tbase(81)
    #TBASE = TBASE[:nt+1]
    with open(f'{prefix}_x.dat', 'w') as fx, open(f'{prefix}_w.dat', 'w') as fw:
        fw.write(f'// DATA_TBASE[{len(TBASE)}] = ''{' + (', '.join([str(x) for x in TBASE])) + '};\n')
        fx.write(f'static double DATA_X[] = ''{\n')
        fw.write(f'static double DATA_W[] = ''{\n')
        for i, tab in enumerate(rys_tab):
            nroots = i + 6
            for it in range(nt):
                ttab = tab[it]
                tbase = TBASE[it]
                print(f'root {nroots}  tbase[{it}] {tbase}')
                fx.write(f'/* root={nroots} base[{it}]={tbase} */\n')
                fw.write(f'/* root={nroots} base[{it}]={tbase} */\n')
                tab_rs, tab_ws = ttab
                fmt_r = [('    %.17e' % x)[-26:] for x in tab_rs.ravel()]
                fmt_w = [('    %.17e' % x)[-26:] for x in tab_ws.ravel()]
                fx.write(',\n'.join(fmt_r))
                fx.write(',\n')
                fw.write(',\n'.join(fmt_w))
                fw.write(',\n')
        fx.write('};\n')
        fw.write('};\n')

def generate_table(path):
    with ProcessPoolExecutor(8) as pe:
        tab = []
        nt = len(TBASE) - 1
        for nroots in range(6, 15):
            res = []
            for tbase in range(nt):
                print(f'root {nroots}  tbase {tbase}')
                fut = pe.submit(tabulate_erf, nroots, tbase)
                res.append(fut)
            res = [fut.result() for fut in res]
            res = np.array(res, dtype=float)
            res = res.reshape(nt,2,nroots,ngrids)
            tab.append(res)
    with open(path, 'wb') as f:
        pickle.dump((TBASE, tab), f)

# when erf(x**.5) ~= 1
def polyfit_large_x_limits(nroots, x, low=None):
    # polynomial fits for large X
    # roots = rx / (x - rx)
    # weights = mpmath.sqrt(mpmath.pi/4/x) * rat
    assert low is None # No large x limits for sr_rys_roots
    r, w = rys_roots_weights_partial(nroots, x, low)
    rx = r * x
    w0 = mpmath.sqrt(mpmath.pi/4/x)
    if low is not None:
        w0 *= mpmath.erfc(low * tt) - mpmath.erfc(tt)
    rat = w / w0
    rat[0] = 1 - rat[1:].sum()
    return rx, rat

def polyfit_small_x_limits(nroots, x, low=None):
    # polynomial fits for large X
    # roots = pr1 * x + pr0;
    # weights = pw1 * w + pw0;
    x0 = x
    x1 = x / 2**10
    r0, w0 = rys_roots_weights(nroots, x0, low)
    r1, w1 = rys_roots_weights(nroots, x1, low)
    pr1 = (r1 - r0) / (x1 - x0)
    pr0 = r0 - pr1 * x0
    pw1 = (w1 - w0) / (x1 - x0)
    pw0 = w0 - pw1 * x0
    return pr0, pr1, pw0, pw1

if __name__ == '__main__':
    #generate_table('rys_rw.pkl')
    pkl2table('./rys', 'rys_rw.pkl')
