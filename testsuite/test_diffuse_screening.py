'''
Primitive-pair screening must keep significant integrals of diffuse, high angular
momentum shells, such as appear in automatically generated auxiliary basis sets.
'''

import os
import sys
import ctypes
import numpy

_cint = ctypes.CDLL(os.path.abspath(os.path.join(__file__, '../../build/libcint.so')))
_cint.CINTgto_norm.restype = ctypes.c_double
_cint.CINTgto_norm.argtypes = [ctypes.c_int, ctypes.c_double]

PTR_EXPCUTOFF = 0
PTR_ENV_START = 20


def int2e_block(shells, expcutoff, with_opt):
    '''(00|nn) for one-primitive shells (l, exponent, z), each on its own center'''
    env = [0.] * PTR_ENV_START
    env[PTR_EXPCUTOFF] = expcutoff
    atm = []
    bas = []
    for ia, (l, a, z) in enumerate(shells):
        atm += [1, len(env), 1, 0, 0, 0]
        env += [0., 0., z]
        bas += [ia, l, 1, 1, 0, len(env), len(env) + 1, 0]
        env += [a, _cint.CINTgto_norm(l, a)]
    atm = numpy.asarray(atm, dtype=numpy.int32)
    bas = numpy.asarray(bas, dtype=numpy.int32)
    env = numpy.asarray(env, dtype=numpy.double)
    n = len(shells)
    di = 2 * shells[0][0] + 1
    dj = 2 * shells[-1][0] + 1
    buf = numpy.zeros(di * di * dj * dj)
    shls = numpy.asarray([0, 0, n - 1, n - 1], dtype=numpy.int32)
    ptr = lambda x: x.ctypes.data_as(ctypes.c_void_p)
    opt = ctypes.c_void_p()
    if with_opt:
        _cint.int2e_optimizer(ctypes.byref(opt), ptr(atm), ctypes.c_int(n),
                              ptr(bas), ctypes.c_int(n), ptr(env))
    _cint.int2e_sph(ptr(buf), None, ptr(shls), ptr(atm), ctypes.c_int(n),
                    ptr(bas), ctypes.c_int(n), ptr(env), opt, None)
    if with_opt:
        _cint.CINTdel_optimizer(ctypes.byref(opt))
    return buf


failed = False
for with_opt in (False, True):
    for l in (2, 3, 4, 5, 6):
        for a in (1e-1, 1e-2, 1e-3, 2.6e-4, 1e-4):
            for shells in ([(l, a, 0.)], [(l, a, 0.), (l, a, 4.)]):
                ref = int2e_block(shells, 400., with_opt)  # effectively no screening
                v = int2e_block(shells, 0., with_opt)      # default cutoff
                err = abs(v - ref).max()
                if err > 1e-12 * abs(ref).max():
                    failed = True
                    print('FAIL: opt=%s l=%d a=%g ncenter=%d max|ref|=%.6g max|err|=%.6g'
                          % (with_opt, l, a, len(shells), abs(ref).max(), err))
if failed:
    sys.exit(1)
print('pass: int2e_sph screening of diffuse high-l shells')
