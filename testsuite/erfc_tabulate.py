import os
import mpmath
import numpy as np
mpmath.mp.dps = 60

#https://www.ams.org/journals/mcom/1981-36-153/S0025-5718-1981-0595058-X/home.html
table_ref = np.array([
mpmath.mpf('0.1177578934567401754080e+01'),
mpmath.mpf('-0.4590054580646477331e-02'),
mpmath.mpf('-0.84249133366517915584e-01'),
mpmath.mpf('0.59209939998191890498e-01'),
mpmath.mpf('-0.26658668435305752277e-01'),
mpmath.mpf('0.9074997670705265094e-02'),
mpmath.mpf('-0.2413163540417608191e-02'),
mpmath.mpf('0.490775836525808632e-03'),
mpmath.mpf('-0.69169733025012064e-04'),
mpmath.mpf('0.4139027986073010e-05'),
mpmath.mpf('0.774038306619849e-06'),
mpmath.mpf('-0.218864010492344e-06'),
mpmath.mpf('0.10764999465671e-07'),
mpmath.mpf('0.4521959811218e-08'),
mpmath.mpf('-0.775440020883e-09'),
mpmath.mpf('-0.63180883409e-10'),
mpmath.mpf('0.28687950109e-10'),
mpmath.mpf('0.194558685e-12'),
mpmath.mpf('-0.965469675e-12'),
mpmath.mpf('0.32525481e-13'),
mpmath.mpf('0.33473119e-13'),
mpmath.mpf('-0.1864563e-14'),
mpmath.mpf('-0.1250795e-14'),
mpmath.mpf('0.74182e-16'),
mpmath.mpf('0.50631e-16'),
mpmath.mpf('-0.2237e-17'),
mpmath.mpf('-0.2187e-17'),
mpmath.mpf('0.27e-19'),
mpmath.mpf('0.97e-19'),
mpmath.mpf('0.3e-20'),
mpmath.mpf('-0.4e-20'),
])

def chebyshev_roots(n):
    return [mpmath.cos(mpmath.pi * (k + .5) / n) for k in range(n)]

def clenshaw_points(n):
    out = []
    for j in range(n):
        c = [mpmath.cos(mpmath.pi * j * (k + 0.5) / n) for k in range(n)]
        out.append(c)
    return out

ngrids = 31
chebrt = np.array(chebyshev_roots(ngrids))
cs = np.array(clenshaw_points(ngrids))

def tabulate_erfc():
    k = mpmath.mpf('3.75')
    values = []
    for t in chebrt:
        x = k*(t+1)/(1-t)
        values.append((1+2*x)*mpmath.exp(x*x)*mpmath.erfc(x))
    values = np.array(values)
    fac = mpmath.mpf(2) / ngrids
    tab = values.dot(cs.T) * fac
    tab[0] /= 2
    return tab

def clenshaw_d1(tab, t):
    d0 = 0
    d1 = tab[ngrids-1]
    for k in reversed(range(1, ngrids-1)):
        d1, d0 = t * 2 * d1 - d0 + tab[k], d1
    return t * d1 - d0 + tab[0]

def erfc(x, tab=None):
    if tab is None:
        tab = tabulate_erfc().astype(float)
    k = 3.75
    t = (x-k)/ (x+k)
    v = clenshaw_d1(tab, t)
    return v / (1 + 2*x) / np.exp(x * x)

if __name__ == '__main__':
    table = tabulate_erfc().astype(float)
    for c in np.arange(-20, 8):
        x = mpmath.mpf('1.5')**c
        ref = mpmath.erfc(x)
        v = erfc(float(x), table)
        print('erfc',1.5**c, 'erfc(x)', v, 'error', float(v / ref - 1))
