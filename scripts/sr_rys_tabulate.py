from concurrent.futures import ProcessPoolExecutor
import os
import pickle
import mpmath
import numpy as np
from rys_roots import DECIMALS, rys_roots_weights_partial, rys_roots_weights
from rys_tabulate import clenshaw_d1, chebyshev_roots, clenshaw_points
mpmath.mp.dps = DECIMALS

ngrids = 14
chebrt = np.array(chebyshev_roots(ngrids))
cs = np.array(clenshaw_points(ngrids))

ngrids_u = 14
chebrt_u = np.array(chebyshev_roots(ngrids_u))
cs_u = np.array(clenshaw_points(ngrids_u))

TBASE = np.append(np.arange(0, 12, 4), 12*1.5**np.arange(20)).astype(int)
#TBASE = np.append(np.arange(0, 10, 5), np.arange(14)**2 * 5 + 10)
print(TBASE)
UBASE = np.sort(np.append([0, .2, .3, .4, .5, .6], .1 / 2**np.arange(4)))
print(UBASE)

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

def get_cheb_u_points(ubase):
    if ubase < 30:
        interval = UBASE[ubase+1] - UBASE[ubase]
        return (chebrt_u+1) * interval/2 + UBASE[ubase]
    else:
        b1 = float(UBASE[ubase+1])
        b0 = float(UBASE[ubase])
        interval = mpmath.log(b1) - mpmath.log(b0)
        return np.array([mpmath.exp((x+1)*interval/2) * b0 for x in chebrt_u])

def cheb_u_interval(x, ubase):
    '''map to interval [-1, 1]'''
    if ubase < 30:
        interval = UBASE[ubase+1] - UBASE[ubase]
        uu = (x - UBASE[ubase]) * 2 / interval - 1
    else:
        b1 = float(UBASE[ubase+1])
        b0 = float(UBASE[ubase])
        interval = mpmath.log(b1) - mpmath.log(b0)
        uu = (mpmath.log(x) - mpmath.log(b0)) * 2 / interval - 1
    return uu

def find_ubase(x):
    return np.searchsorted(UBASE, x) - 1

def tabulate_erfc(nroots, tbase, ubase):
    tabulate_points = get_cheb_t_points(tbase)
    boundary_points = get_cheb_u_points(ubase)

    fac = (mpmath.mpf(2) / ngrids) * (mpmath.mpf(2) / ngrids_u)

    rs = []
    ws = []
    for low in boundary_points:
        for x in tabulate_points:
            r, w = rys_roots_weights(nroots, x, low)
            rs.append(r)
            ws.append(w)
    rs = np.array(rs).reshape(ngrids_u, ngrids, nroots)
    ws = np.array(ws).reshape(ngrids_u, ngrids, nroots)

    # einsum('lkn,il,jk->nji', rs, cs)
    tab_rs = np.array([cs_u.dot(rs[:,:,ir]).dot(cs.T) * fac for ir in range(nroots)])
    tab_ws = np.array([cs_u.dot(ws[:,:,ir]).dot(cs.T) * fac for ir in range(nroots)])
    return tab_rs, tab_ws

def polyfit_erfc(nroots, x, low):
    tbase = find_tbase(x)
    ubase = find_ubase(low)
    tt = cheb_t_interval(x, tbase)
    uu = cheb_u_interval(low, ubase)

    tab_rs, tab_ws = tabulate_erfc(nroots, tbase, ubase)
    im = clenshaw_d1(tab_rs.astype(float), uu, nroots)
    rr = clenshaw_d1(im, tt, nroots)

    im = clenshaw_d1(tab_ws.astype(float), uu, nroots)
    ww = clenshaw_d1(im, tt, nroots)
    #if x * low**2 < DECIMALS*.7:
    #    factor = mpmath.exp(-x * low**2)
    #    ww = [w * factor for w in ww]
    return rr, ww

def pkl2table(prefix, pklfile):
    with open(pklfile, 'rb') as f:
        TBASE, UBASE, rys_tab = pickle.load(f)
    UBASE = UBASE.round(6)
    TBASE = TBASE.round(6)
    us = [(0, 5), (0, 5),          (4, len(UBASE)),]
    ts = [(0, 4), (3, len(TBASE)), (0, 11),        ]
    for j in range(3):
        u0, u1 = us[j]
        t0, t1 = ts[j]
        with open(f'{prefix}_part{j}_x.dat', 'w') as fx, open(f'{prefix}_part{j}_w.dat', 'w') as fw:
            fw.write(f'static double SR_DATA{j}_UBASE[{u1-u0}] = ''{' + (', '.join([str(x) for x in UBASE[u0:u1]])) + '};\n')
            fw.write(f'static double SR_DATA{j}_TBASE[{t1-t0}] = ''{' + (', '.join([str(x) for x in TBASE[t0:t1]])) + '};\n')
            #fx.write(f'static double SR_DATA{j}_UBASE[{u1-u0}] = ''{' + (', '.join([str(x) for x in UBASE[u0:u1]])) + '};\n')
            #fx.write(f'static double SR_DATA{j}_TBASE[{t1-t0}] = ''{' + (', '.join([str(x) for x in TBASE[t0:t1]])) + '};\n')
            fx.write(f'static double SR_DATA{j}_X[] = ''{\n')
            fw.write(f'static double SR_DATA{j}_W[] = ''{\n')
            for i, tab in enumerate(rys_tab):
                nroots = i + 1
                #fx.write(f'static const double SR_DATA{j}_X{nroots}[] = ''{\n')
                #fw.write(f'static const double SR_DATA{j}_W{nroots}[] = ''{\n')
                for iu in range(u0, u1-1):
                    utab = tab[iu]
                    ubase = UBASE[iu].round(6)
                    for it in range(t0, t1-1):
                        ttab = utab[it]
                        tbase = TBASE[it].round(6)
                        print(f'root {nroots}  ubase[{iu}] {ubase}  tbase[{it}] {tbase}')
                        fx.write(f'/* root={nroots} ubase[{iu}]={ubase} tbase[{it}]={tbase} */\n')
                        fw.write(f'/* root={nroots} ubase[{iu}]={ubase} tbase[{it}]={tbase} */\n')
                        tab_rs, tab_ws = ttab
                        fmt_r = [('    %.17e' % x)[-26:] for x in tab_rs.ravel()]
                        fmt_w = [('    %.17e' % x)[-26:] for x in tab_ws.ravel()]
                        fx.write(',\n'.join(fmt_r))
                        fx.write(',\n')
                        fw.write(',\n'.join(fmt_w))
                        fw.write(',\n')
                #fx.write('};\n')
                #fw.write('};\n')
            fx.write('};\n')
            fw.write('};\n')

def generate_table(path):
    with ProcessPoolExecutor(8) as pe:
        tab = []
        nt = len(TBASE) - 1
        nu = len(UBASE) - 1
        for nroots in range(1, 6):
            res = []
            for ubase in range(nu):
                for tbase in range(nt):
                    print(f'root {nroots}  ubase {ubase}  tbase {tbase}')
                    fut = pe.submit(tabulate_erfc, nroots, tbase, ubase)
                    res.append(fut)
            res = [fut.result() for fut in res]
            res = np.array(res, dtype=float)
            res = res.reshape(nu,nt,2,nroots,ngrids_u,ngrids)
            tab.append(res)
    with open(path, 'wb') as f:
        pickle.dump((TBASE, UBASE, tab), f)

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
    #generate_table('sr_rys_rw.pkl')
    pkl2table('./sr_roots', 'sr_rys_rw.pkl')
