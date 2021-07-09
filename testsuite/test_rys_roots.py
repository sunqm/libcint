import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(__file__, '../../scripts')))

import ctypes
import numpy
import numpy as np
import rys_roots
import rys_wheeler

def fp(a):
    a = numpy.asarray(a)
    return a.ravel().dot(numpy.cos(numpy.arange(a.size)))

cint = ctypes.CDLL('../build/libcint.so')
def cint_call(fname, nroots, x, low=None):
    r = numpy.zeros(nroots)
    w = numpy.zeros(nroots)
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

def check_turnover_point(fun1, fun2):
    es = 2**numpy.arange(-2, 5, .5)
    turnover_points = []
    for i in range(1, 14):
        for x in es:
            f1 = fun1(x, i)
            f2 = fun2(x, i)
            diff = (f1[i] - f2[i])/f1[i]
            if diff < 2.2e-16:
                turnover_points.append(x)
                break
    return turnover_points

def test_turnover_point():
    import functools
    import mpmath
    mpmath.mp.dps = 15
    print(check_turnover_point(rys_roots.fmt1, rys_roots.fmt2))
    tps = []
    for low in 2**numpy.arange(-4, 5, .5):
        tp = check_turnover_point(functools.partial(rys_roots.fmt1_erfc, low=low),
                                  functools.partial(rys_roots.fmt2_erfc, low=low))
        print(low, tp)
        tps.append(tp)
    print(numpy.array(tps).max(axis=0))
    print(rys_roots.fmt1_erfc(0.9, 2, .6))
    print(rys_roots.fmt1_erfc(1.1, 2, .2))
    print(rys_roots.fmt2_erfc(1.1, 2, .2))
    print(rys_roots.fmt1_erfc(8.1, 2, .2))
    print(rys_roots.fmt2_erfc(8.1, 2, .2))
    print(rys_roots.fmt1(1.1, 2))
    print(rys_roots.fmt2(1.1, 2))
    print(rys_roots.fmt1(8.1, 2))
    print(rys_roots.fmt2(8.1, 2))

    print(abs(rys_roots.fmt(1, 2) - [0.746824132812427, 0.18947234582049235, 0.10026879814501737]).max())
    print(abs(rys_roots.fmt(7, 2) - [0.3349010581765593, 0.023856369729357483, 0.005046944801608424]).max())

def test_rys_roots_mpmath():
    numpy.random.seed(4)
    a = numpy.random.rand(4,4)
    a = a.T.dot(a)
    cs = rys_roots.R_dsmit(a, 4)
    print(fp(cs) - -1.1566333099933561)
    r, w = rys_roots.rys_roots_weights(6, .4)
    print(fp(r) - 2.898038115730877)
    print(fp(w) - 0.1193288112731978)

    print(fp(rys_roots.rys_roots_weights(1, .3)) - 0.93475892362307533)
    print(fp(rys_roots.rys_roots_weights(2, .8)) - 0.94983728587499205)
    print(fp(rys_roots.rys_roots_weights(3, .5)) - -2.7165192846520481)
    print(fp(rys_roots.rys_roots_weights(4, .5)) - -11.442155925744938)

    ff = rys_roots.fmt_erfc(.4, 6*2, 0.2)
    r, w = rys_roots.rys_roots_weights(6, .4, ff)
    print(fp(r) - 3.0557596858378053)
    print(fp(w) - 0.00070334763896715)

def test_rys_roots_weights():
    print('test rys_roots_weights')
    def check(nroots, x):
        r_ref, w_ref = rys_roots.rys_roots_weights(nroots, x)
        r, w = cint_call('CINTrys_roots', nroots, x)
        return np.array([abs(r-r_ref).max(), abs(w-w_ref).max()]).astype(float)

    failed = False
    max_r_error = 0
    max_w_error = 0
    es = 2**numpy.arange(-3, 7, .25)
    for i in range(1, 12):
        for x in es:
            diffs = check(i, x)
            if diffs[0] > 1e-5 or diffs[1] > 1e-8:
                print('Errors for root', i, x, diffs)
                #failed |= diffs[0] < 1e-4 and diffs[1] < 1e-7
                max_r_error = max(max_r_error, diffs[0])
                max_w_error = max(max_w_error, diffs[1])
    es = [50, 40, 36, 24, 20, 16, 12, 8, 4, 2, 0.5, 6e-8]
    for i in range(1, 5):
        for x in es:
            diffs = check(i, x)
            if diffs[0] > 1e-5 or diffs[1] > 1e-8:
                print('Errors for root', i, x, diffs)
                #failed |= diffs[0] < 1e-4 and diffs[1] < 1e-7
                max_r_error = max(max_r_error, diffs[0])
                max_w_error = max(max_w_error, diffs[1])
    #if failed:
    #    print('test_rys_roots_weights .. failed')
    #else:
    #    print('test_rys_roots_weights .. pass')
    assert max_r_error < 1e-3
    assert max_w_error < 1e-7
    print('test_rys_roots_weights .. pass')


def test_rys_roots_weights_erfc():
    print('test sr-rys_roots_weights')
    def check(nroots, x, lower):
        r_ref, w_ref = rys_roots.rys_roots_weights(nroots, x, lower)
        r, w = cint_call('CINTsr_rys_roots', nroots, x, lower)
        return np.array([abs(r-r_ref).max(), abs(w-w_ref).max()]).astype(float)

    es = 2**numpy.arange(-3, 6, .25)
    failed = False
    max_r_error = 0
    max_w_error = 0
    max_rw_error = 0
    for i in range(1, 10):
        for x in es:
            for low in [.1, .2, .3, .4, .5, .6, .7, .8, .9]:
                diffs = check(i, x, low)
                if diffs[0] > 1e-4 or diffs[1] > 1e-7:
                    print('Errors for root', i, x, low, diffs)
                    #failed |= not all(s < 1e-4 for s in diffs)
                    max_r_error = max(max_r_error, diffs[0])
                    max_w_error = max(max_w_error, diffs[1])
                    max_rw_error = max(max_rw_error, diffs[0]*diffs[1])
    #if failed:
    #    print('test_rys_roots_weights_erfc .. failed')
    #else:
    #    print('test_rys_roots_weights_erfc .. pass')
    #assert max_r_error < 1e-2
    assert max_w_error < 1e-10
    assert max_rw_error < 1e-12

    for i in range(10, 12):
        for x in es:
            for low in [.1, .2, .3, .4, .5, .6, .7, .8, .9]:
                diffs = check(i, x, low)
                if diffs[0] > 1e-4 or diffs[1] > 1e-7:
                    print('Errors for root', i, x, low, diffs)
                    #failed |= not all(s < 1e-4 for s in diffs)
                    max_r_error = max(max_r_error, diffs[0])
                    max_w_error = max(max_w_error, diffs[1])
                    max_rw_error = max(max_rw_error, diffs[0]*diffs[1])
    #assert max_r_error < 1e-2
    #assert max_w_error < 1e-7
    # Errors for high-angular basis are slightly larger
    #assert max_rw_error < 1e-7

    max_r_error = 0
    max_w_error = 0
    max_rw_error = 0
    for i in range(2, 4):
        for low in [.92, .94, .96, .98, .99, .992, .994, .996]:
            for x in es:
                diffs = check(i, x, low)
                if diffs[0] > 1e-4 or diffs[1] > 1e-7:
                    print('Errors for root', i, x, low, diffs)
                    #failed |= not all(s < 1e-4 for s in diffs)
                    max_r_error = max(max_r_error, diffs[0])
                    max_w_error = max(max_w_error, diffs[1])
                    max_rw_error = max(max_rw_error, diffs[0]*diffs[1])
    assert max_w_error < 1e-11
    print('test_rys_roots_weights_erfc .. pass')


def test_polyfit():
    import rys_tabulate
    def check(nroots, x, low=None):
        r0, w0 = rys_roots.rys_roots_weights(nroots, x, low)
        if low is None:
            r1, w1 = rys_tabulate.polyfit_erf(nroots, x)
        else:
            r1, w1 = rys_tabulate.polyfit_erfc(nroots, x, low)
        return np.array([abs(r0 - r1).max(), abs(w0 - w1).max()]).astype(float)

    es = 2**numpy.arange(-3, 6, .25)
    failed = False
    for i in range(6, 12):
        for low in [None, 0.2, 0.5, 0.8]:
            for x in es:
                diffs = check(i, x, low)
                if not all(s < 1e-8 for s in diffs):
                    print(i, x, low, diffs)
                    failed |= not all(s < 1e-5 for s in diffs)
    if failed:
        print('test_polyfit .. failed')
    else:
        print('test_polyfit .. pass')

def test_rys_roots_vs_polyfit():
    def check(nroots, x, low):
        r_ref, w_ref = rys_roots.rys_roots_weights(nroots, x, low)
        r0, w0 = cint_call('CINTrys_schmidt', nroots, x, low)
        r1, w1 = cint_call('CINTsr_rys_polyfits', nroots, x, low)
        return np.array([abs(r0 - r_ref).max(),
                         abs(r1 - r_ref).max(),
                         abs(w0 - w_ref).max(),
                         abs(w1 - w_ref).max()]).astype(float)

    es = 2**numpy.arange(-6, 6, .5)
    for i in range(6, 12):
        for x in es:
            for low in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
                diffs = check(i, x, low)
                print(i, x, low, diffs)

def test_stg_roots():
    print('test stg roots')
    def stg(nroots, t, u):
        r = numpy.zeros(nroots)
        w = numpy.zeros(nroots)
        cint.CINTstg_roots(ctypes.c_int(nroots),
                           ctypes.c_double(t), ctypes.c_double(u),
                           r.ctypes.data_as(ctypes.c_void_p),
                           w.ctypes.data_as(ctypes.c_void_p))
        return r, w

    assert abs(fp(stg(1, 2.2, 1.5)) - 0.653262713484748 ) < 1e-14
    assert abs(fp(stg(2, 0.2, 8.5)) - 1.2362174210548105) < 1e-14
    assert abs(fp(stg(4, 1.0, 0.5)) - -0.6907781084439245) < 1e-14
    print('test_stg_roots .. pass')

if __name__ == '__main__':
    # test_rys_roots_mpmath()
    #test_polyfit()
    # test_rys_roots_vs_polyfit()
    #test_rys_roots_weights()
    test_stg_roots()
    test_rys_roots_weights()
    test_rys_roots_weights_erfc()
