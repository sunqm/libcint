'''
Accurate for large x (x > 200) and medium boundary (0.1 ... 0.6)
'''
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
#INTERVAL = 5
#TBASE = np.append(0, np.arange(1, 65, INTERVAL)) * .1
TBASE = np.arange(1, 8) ** 2 * .1 + 1.3
UBASE = np.arange(7) * .1
print(TBASE)
print(UBASE)

ngrids_u = 14
chebrt_u = np.array(chebyshev_roots(ngrids_u))
cs_u = np.array(clenshaw_points(ngrids_u))

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
            x0 = (x/low)**2
            r, w = rys_roots_weights(nroots, x0, low)
            rs.append(r)
            ws.append(w)
    rs = np.array(rs).reshape(ngrids_u, ngrids, nroots)
    ws = np.array(ws).reshape(ngrids_u, ngrids, nroots)

    # einsum('lkn,il,jk->nji', rs, cs)
    tab_rs = np.array([cs_u.dot(rs[:,:,ir]).dot(cs.T) * fac for ir in range(nroots)])
    tab_ws = np.array([cs_u.dot(ws[:,:,ir]).dot(cs.T) * fac for ir in range(nroots)])
    return tab_rs, tab_ws

def polyfit_erfc(nroots, x, low):
    x1 = x**.5*low
    tbase = find_tbase(x1)
    ubase = find_ubase(low)
    tt = cheb_t_interval(x1, tbase)
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
    with open(f'{prefix}_x.dat', 'w') as fx, open(f'{prefix}_w.dat', 'w') as fw:
        fw.write(f'static double SR_DATAL_UBASE[{len(UBASE)}] = ''{' + (', '.join([str(x) for x in UBASE])) + '};\n')
        fw.write(f'static double SR_DATAL_TBASE[{len(TBASE)}] = ''{' + (', '.join([str(x) for x in TBASE])) + '};\n')
        #fx.write(f'static double SR_DATAL_UBASE[{len(UBASE)}] = ''{' + (', '.join([str(x) for x in UBASE])) + '};\n')
        #fx.write(f'static double SR_DATAL_TBASE[{len(TBASE)}] = ''{' + (', '.join([str(x) for x in TBASE])) + '};\n')
        fx.write(f'static double SR_DATAL_X[] = ''{\n')
        fw.write(f'static double SR_DATAL_W[] = ''{\n')
        for i, tab in enumerate(rys_tab):
            nroots = i + 1
            #fx.write(f'static const double SR_DATAL_X{nroots}[] = ''{\n')
            #fw.write(f'static const double SR_DATAL_W{nroots}[] = ''{\n')
            for iu, utab in enumerate(tab):
                ubase = UBASE[iu].round(6)
                for it, ttab in enumerate(utab):
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

if __name__ == '__main__':
    #generate_table('sr_rys_rw_xlarge.pkl')
    pkl2table('./sr_roots_part3', 'sr_rys_rw_xlarge.pkl')
