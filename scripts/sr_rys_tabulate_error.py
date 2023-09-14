import pickle
import numpy as np
from scipy.special import erf, erfc
from rys_roots import rys_roots_weights, rys_roots_weights_partial
import sr_rys_tabulate
import sr_rys_tabulate1
import rys_tabulate
from sr_rys_tabulate import clenshaw_d1

np.random.seed(1)
xs = np.sort(np.append(np.random.rand(200) * 1000, 3**((np.random.rand(40)-.3) * 7)))
xs[0] = 1e-15
print(xs)
np.random.seed(4)
lows = np.sort(np.append(np.random.rand(30)*.5, np.random.rand(20)*.1))
lows[0] = 1e-3
print(lows)

fil = 'sr_rys_rw.pkl'
TBASE, UBASE, tab0 = pickle.load(open(fil, 'rb'))
sr_rys_tabulate.TBASE = TBASE
sr_rys_tabulate.UBASE = UBASE

fil = 'sr_rys_rw_xlarge.pkl'
TBASE, UBASE, tab1 = pickle.load(open(fil, 'rb'))
sr_rys_tabulate1.TBASE = TBASE
sr_rys_tabulate1.UBASE = UBASE

nroots = 6
x = 5.14
low = .49
for nroots in range(1, 6):
    for x in xs:
        for low in lows:
            xl = x**.5*low
            if xl >= 6.1:
                continue
            if low < .1 or x < 196:
                it = np.searchsorted(sr_rys_tabulate.TBASE, x) - 1
                iu = np.searchsorted(sr_rys_tabulate.UBASE, low) - 1
                tt = sr_rys_tabulate.cheb_t_interval(x, it)
                uu = sr_rys_tabulate.cheb_u_interval(low, iu)
                tab_rs, tab_ws = tab0[nroots-1][iu][it]
            else:
                it = np.searchsorted(sr_rys_tabulate1.TBASE, xl) - 1
                iu = np.searchsorted(sr_rys_tabulate1.UBASE, low) - 1
                tt = sr_rys_tabulate1.cheb_t_interval(xl, it)
                uu = sr_rys_tabulate1.cheb_u_interval(low, iu)
                tab_rs, tab_ws = tab1[nroots-1][iu][it]
            im = clenshaw_d1(tab_rs.astype(float), uu, nroots)
            rr = np.array(clenshaw_d1(im, tt, nroots), dtype=float)
            im = clenshaw_d1(tab_ws.astype(float), uu, nroots)
            ww = np.array(clenshaw_d1(im, tt, nroots), dtype=float)
            ref = np.array(rys_roots_weights(nroots, x, low), dtype=float)
            rr /= rr + 1
            ref[0] /= ref[0] + 1
            diff1, diff2 = abs((rr - ref[0])/ref[0]).max(), abs((ww - ref[1])).max()
            if diff1 > 1e-12 or diff2 > 1e-13:
                print(nroots, x, low, diff1, diff2)

np.random.seed(1)
xs = np.sort(np.append(np.random.rand(100) * 100, 3**((np.random.rand(30)-.4) * 6)))
xs[0] = 1e-15
fil = 'rys_rw.pkl'
TBASE, tab = pickle.load(open(fil, 'rb'))
rys_tabulate.TBASE = TBASE

for nroots in range(1, 15):
    for x in xs:
        if x >= 100:
            continue

        it = np.searchsorted(rys_tabulate.TBASE, x) - 1
        tt = rys_tabulate.cheb_t_interval(x, it)
        tab_rs, tab_ws = tab[nroots-1][it]
        rr = np.array(clenshaw_d1(tab_rs.astype(float), tt, nroots), dtype=float)
        ww = np.array(clenshaw_d1(tab_ws.astype(float), tt, nroots), dtype=float)
        ref = np.array(rys_roots_weights(nroots, x), dtype=float)
        diff1, diff2 = abs((rr - ref[0])/ref[0]).max(), abs((ww - ref[1])/ref[1]).max()
        if diff1 > 1e-12 or diff2 > 1e-13:
            print(nroots, x, diff1, diff2)
