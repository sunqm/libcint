import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(__file__, '../../scripts')))

import ctypes
import numpy as np
import rys_wheeler
import rys_roots

cint = ctypes.CDLL('../build/libcint.so')

def test_rys_wheeler():
    def call_rw(nroots, t, low):
        r = np.zeros(nroots)
        w = np.zeros(nroots)
        cint.roots_and_weights(ctypes.c_int(nroots),
                               ctypes.c_double(t), ctypes.c_double(low),
                               r.ctypes.data_as(ctypes.c_void_p),
                               w.ctypes.data_as(ctypes.c_void_p))
        return r, w

    n = 3
    t = 6.4
    low = .2

    rref, wref = rys_wheeler.roots_and_weights_partial(n, *rys_wheeler.shifted_jacobi_moments(n*2, t, 0))
    rt, wt = call_rw(n, t, 0)
    #print(rys_roots.rys_roots_weights(n, t))
    print(abs(rref - rt).max(), abs(wref - wt).max())

    rref, wref = rys_wheeler.roots_and_weights_partial(n, *rys_wheeler.shifted_jacobi_moments(n*2, t, low))
    rt, wt = call_rw(n, t, low)
    print(abs(rref - rt).max(), abs(wref - wt).max())
    #print(rys_roots.rys_roots_weights(n, t, low))

    alpha, beta, mus = mp_rys_wheeler.shifted_jacobi_moments(n, t, low)
    alpha = alpha.astype(float)
    beta = beta.astype(float)
    mus = dble_jacobi_mus(n-1, t, low)
    rt, wt = rys_wheeler.roots_and_weights_partial(n//2, alpha, beta, mus)
    rt = rt / (1 - rt)
    print(abs(rref - rt).max(), abs(wref - wt).max())

if __name__ == '__main__':
    # test_rys_wheeler()
    n = 16
    t = .1
    low = .8
    mus_ref = np.array(rys_wheeler.jacobi_mus(n, t, low)).astype(float)
    mus = rys_wheeler.dble_naive_jacobi_mus(n-1, t, low).astype(float)
    print(abs((mus_ref - mus) / mus_ref).max())
