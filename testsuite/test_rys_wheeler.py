import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(__file__, '../../scripts')))

import ctypes
import numpy as np
import rys_wheeler
import rys_roots
import rys_tabulate

cint = ctypes.CDLL('../build/libcint.so')
def cint_call(fname, nroots, x, low=None):
    r = np.zeros(nroots)
    w = np.zeros(nroots)
    fun = getattr(cint, fname)
    if low is None:
        fun(ctypes.c_int(nroots), ctypes.c_double(x),
            r.ctypes.data_as(ctypes.c_void_p),
            w.ctypes.data_as(ctypes.c_void_p))
    else:
        fun(ctypes.c_int(nroots), ctypes.c_double(x), ctypes.c_double(low),
            r.ctypes.data_as(ctypes.c_void_p),
            w.ctypes.data_as(ctypes.c_void_p))
    return r, w

def check_rys_wheeler():
    def check(nroots, x):
        r_ref, w_ref = rys_roots.rys_roots_weights(nroots, x)
        r0, w0 = cint_call('CINTrys_schmidt', nroots, x, 0.)
        #r0, w0 = rys_tabulate.polyfit_erf(nroots, x)
        r1, w1 = cint_call('CINTrys_laguerre', nroots, x, 0.)
        r2, w2 = cint_call('CINTrys_jacobi', nroots, x, 0.)
        return np.array([abs(r0 - r_ref).max(),
                         abs(r1 - r_ref).max(),
                         abs(r2 - r_ref).max(),
                         abs(w0 - w_ref).max(),
                         abs(w1 - w_ref).max(),
                         abs(w2 - w_ref).max()]).astype(float)

    failed = False
    es = 2**np.arange(-5, 8, .5)
    for i in range(1, 13):
        for x in es:
            diffs = check(i, x)
            print(i, x, diffs)

def check_rys_wheeler1():
    def check(nroots, x):
        r_ref, w_ref = rys_roots.rys_roots_weights(nroots, x)
        r0, w0 = cint_call('CINTlrys_laguerre', nroots, x, 0.)
        r1, w1 = cint_call('CINTlrys_jacobi', nroots, x, 0.)
        return np.array([abs(r0 - r_ref).max(),
                         abs(r1 - r_ref).max(),
                         abs(w0 - w_ref).max(),
                         abs(w1 - w_ref).max()]).astype(float)

    failed = False
    es = 2**np.arange(-5, 8, .5)
    for i in range(10, 21):
        for x in es:
            diffs = check(i, x)
            print(i, x, diffs)

def check_rys_wheeler2():
    def check(nroots, x):
        r_ref, w_ref = rys_roots.rys_roots_weights(nroots, x)
        r0, w0 = cint_call('CINTqrys_laguerre', nroots, x, 0.)
        r1, w1 = cint_call('CINTqrys_jacobi', nroots, x, 0.)
        return np.array([abs(r0 - r_ref).max(),
                         abs(r1 - r_ref).max(),
                         abs(w0 - w_ref).max(),
                         abs(w1 - w_ref).max()]).astype(float)

    failed = False
    es = 2**np.arange(-5, 8, .5)
    for i in range(13, 25):
        for x in es:
            diffs = check(i, x)
            print(i, x, diffs)


def check_rys_wheeler_sr():
    def check(nroots, x, low):
        r_ref, w_ref = rys_roots.rys_roots_weights(nroots, x, low)
        r0, w0 = cint_call('sr_rys_roots', nroots, x, low)
        #r0, w0 = cint_call('CINTsr_rys_polyfits', nroots, x, low)
        r1, w1 = cint_call('CINTrys_laguerre', nroots, x, low)
        r2, w2 = cint_call('CINTrys_jacobi', nroots, x, low)
        return np.array([abs(r0 - r_ref).max(),
                         abs(r1 - r_ref).max(),
                         abs(r2 - r_ref).max(),
                         abs(w0 - w_ref).max(),
                         abs(w1 - w_ref).max(),
                         abs(w2 - w_ref).max()]).astype(float)

    failed = False
    es = 2**np.arange(-5, 8, .5)
    for i in range(1, 13):
        for low in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            for x in es:
                diffs = check(i, x, low)
                print(i, low, x, diffs)

    # TODO: test the diffcult case. When 40 < low**2*x < 60, rys_wheeler is
    # very inaccurate
    #print(roots_and_weights(7, 77., 0.76))
    #print(roots_and_weights(7, 65., 0.8))

def check_rys_wheeler_sr1():
    def check(nroots, x, low):
        r_ref, w_ref = rys_roots.rys_roots_weights(nroots, x, low)
        r0, w0 = cint_call('CINTlrys_laguerre', nroots, x, low)
        r1, w1 = cint_call('CINTlrys_jacobi', nroots, x, low)
        return np.array([abs(r0 - r_ref).max(),
                         abs(r1 - r_ref).max(),
                         abs(w0 - w_ref).max(),
                         abs(w1 - w_ref).max()]).astype(float)

    failed = False
    es = 2**np.arange(-5, 8, .5)
    for i in range(4, 15):
        for low in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            for x in es:
                diffs = check(i, x, low)
                print(i, low, x, diffs)

def check_rys_wheeler_sr2():
    def check(nroots, x, low):
        r_ref, w_ref = rys_roots.rys_roots_weights(nroots, x, low)
        r0, w0 = cint_call('CINTqrys_laguerre', nroots, x, low)
        r1, w1 = cint_call('CINTqrys_jacobi', nroots, x, low)
        return np.array([abs(r0 - r_ref).max(),
                         abs(r1 - r_ref).max(),
                         abs(w0 - w_ref).max(),
                         abs(w1 - w_ref).max()]).astype(float)

    failed = False
    es = 2**np.arange(-5, 8, .5)
    for i in range(10, 25):
        for low in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            for x in es:
                diffs = check(i, x, low)
                print(i, low, x, diffs)


def test_rys_wheeler1():
    def check(nroots, x, low):
        r_ref, w_ref = rys_roots.rys_roots_weights(nroots, x, low)
        r, w = cint_call('CINTrys_wheeler', nroots, x, low)
        return np.array([abs(r-r_ref).max(), abs(w-w_ref).max()]).astype(float)

    failed = False
    n = 5
    t = 94.1
    low = .8

    diffs = check(n, t, 0)
    failed |= not all(s < 1e-5 for s in diffs)
    diffs = check(n, t, low)
    failed |= not all(s < 1e-5 for s in diffs)

    n = 14
    t = 0.5
    low = .8

    diffs = check(n, t, 0)
    failed |= not all(s < 1e-5 for s in diffs)
    diffs = check(n, t, low)
    failed |= not all(s < 1e-5 for s in diffs)

    es = 2**np.arange(-3, 7, .5)
    for i in range(1, 13):
        for x in es:
            diffs = check(i, x)
            if any(s > 1e-8 for s in diffs):
                print(i, x, diffs)
                failed |= any(s > 1e-5 for s in diffs)

if __name__ == '__main__':
    np.set_printoptions(3)
    check_rys_wheeler()

    #check_rys_wheeler()
    #check_rys_wheeler1()
    #check_rys_wheeler2()
    #check_rys_wheeler_sr()
    #check_rys_wheeler_sr1()
    #check_rys_wheeler_sr2()
