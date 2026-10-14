#!/usr/bin/env python
'''
Screening of diffuse primitives of high angular momentum.

When all exponents of normalized primitive shells on one center are scaled by
s, an integral scales exactly as s**p, where p follows from the dimensions of
the operator: p = 1/2 for (ab|cd), -1/4 for (ab|c), -1 for (a|b), 0 for <a|b>,
1/2 for <a|1/r|b> and 1 for <a|-1/2 nabla^2|b>. The range-separated Coulomb
operator erfc(omega r)/r follows the same rule when omega is scaled by
sqrt(s). The integrals at s = 1 therefore serve as the reference for all s,
and a primitive-pair estimate that misses a factor of the exponents shows up
as integrals that are screened to zero for small s.
'''

import os
import sys
import ctypes
import numpy

_cint = numpy.ctypeslib.load_library('libcint', os.path.abspath(os.path.join(__file__, '../../build')))
_cint.CINTgto_norm.restype = ctypes.c_double

PTR_EXPCUTOFF = 0
PTR_RANGE_OMEGA = 8
PTR_ENV_START = 20
ATM_SLOTS = 6
BAS_SLOTS = 8
ANG_OF = 1
NPRIM_OF = 2
NCTR_OF = 3
PTR_EXP = 5
PTR_COEFF = 6

SCALES = [1e-1, 1e-2, 1e-3, 5e-4, 3e-4, 1e-4, 1e-5, 1e-6]

def make_env(shells, scale, omega=0):
    '''shells: list of (l, exponent) on one center at the origin, where
    exponent may be a tuple for a contraction of normalized primitives. The
    exponents are multiplied by scale, and omega by sqrt(scale).'''
    atm = numpy.zeros((1, ATM_SLOTS), dtype=numpy.int32)
    atm[0,0] = 1
    atm[0,1] = PTR_ENV_START
    bas = numpy.zeros((len(shells), BAS_SLOTS), dtype=numpy.int32)
    nprim = [numpy.size(a) for l, a in shells]
    env = numpy.zeros(PTR_ENV_START + 3 + 2*sum(nprim))
    env[PTR_RANGE_OMEGA] = omega * numpy.sqrt(scale)
    off = PTR_ENV_START + 3
    for i, (l, a) in enumerate(shells):
        n = nprim[i]
        a = numpy.ravel(a) * scale
        env[off:off+n] = a
        env[off+n:off+2*n] = [_cint.CINTgto_norm(ctypes.c_int(l), ctypes.c_double(x))
                              for x in a]
        bas[i,ANG_OF] = l
        bas[i,NPRIM_OF] = n
        bas[i,NCTR_OF] = 1
        bas[i,PTR_EXP] = off
        bas[i,PTR_COEFF] = off + n
        off += 2*n
    return atm, bas, env

def eval_intor(intor, shls, atm, bas, env, with_opt):
    nf = 1
    for i in shls:
        nf *= 2 * bas[i,ANG_OF] + 1
    buf = numpy.zeros(nf)
    args = (atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(atm.shape[0]),
            bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(bas.shape[0]),
            env.ctypes.data_as(ctypes.c_void_p))
    opt = ctypes.c_void_p()
    if with_opt:
        getattr(_cint, intor + '_optimizer')(ctypes.byref(opt), *args)
    c_shls = (ctypes.c_int * len(shls))(*shls)
    getattr(_cint, intor + '_sph')(buf.ctypes.data_as(ctypes.c_void_p),
                                   ctypes.c_void_p(), c_shls, *args, opt,
                                   ctypes.c_void_p())
    if with_opt:
        _cint.CINTdel_optimizer(ctypes.byref(opt))
    return buf

def check(intor, power, shells, shls, omega=0, with_opt=False):
    ref = eval_intor(intor, shls, *make_env(shells, 1., omega), with_opt)
    assert abs(ref).max() > 0
    failed = False
    for s in SCALES:
        val = eval_intor(intor, shls, *make_env(shells, s, omega), with_opt)
        expected = ref * s**power
        err = abs(val - expected).max() / abs(expected).max()
        if err > 1e-10:
            print('FAILED %s shells %s shls %s omega %g opt %s scale %g: rel. error %.3e'
                  % (intor, shells, shls, omega, with_opt, s, err))
            failed = True
    return failed

failed = False
for l in (3, 4, 5):
    shells = [(l, 1.)]
    for with_opt in (False, True):
        failed |= check('int2e', .5, shells, (0, 0, 0, 0), with_opt=with_opt)
        failed |= check('int2e', .5, shells, (0, 0, 0, 0), -.5, with_opt)
    # (ll|l) vanishes for odd l; use a d shell for the third index
    failed |= check('int3c2e', -.25, shells + [(2, 1.)], (0, 0, 1))
    failed |= check('int3c2e', -.25, shells + [(2, 1.)], (0, 0, 1), -.5)
    failed |= check('int2c2e', -1., shells, (0, 0))
    failed |= check('int1e_ovlp', 0., shells, (0, 0))
    failed |= check('int1e_kin', 1., shells, (0, 0))
    failed |= check('int1e_nuc', .5, shells, (0, 0))

# shells with different l and different exponents
shells = [(2, 1.), (5, 3.)]
for shls in ((0, 1, 0, 1), (0, 1, 1, 0), (1, 1, 0, 0), (0, 0, 1, 1)):
    for with_opt in (False, True):
        failed |= check('int2e', .5, shells, shls, with_opt=with_opt)
        failed |= check('int2e', .5, shells, shls, -.5, with_opt)
for shls in ((0, 1, 1), (1, 1, 0)):
    failed |= check('int3c2e', -.25, shells, shls)
failed |= check('int1e_ovlp', 0., shells, (1, 1))
failed |= check('int1e_kin', 1., shells, (1, 1))
failed |= check('int1e_nuc', .5, shells, (1, 1))

# contracted shells whose exponents are not sorted in decreasing order
for exps in ((1., 1e4), (1e4, 1., 100.)):
    shells = [(4, exps)]
    for with_opt in (False, True):
        failed |= check('int2e', .5, shells, (0, 0, 0, 0), with_opt=with_opt)
        failed |= check('int2e', .5, shells, (0, 0, 0, 0), -.5, with_opt)
    failed |= check('int3c2e', -.25, shells + [(2, 1.)], (0, 0, 1))
    failed |= check('int1e_nuc', .5, shells, (0, 0))

if failed:
    sys.exit(1)
print('pass')
