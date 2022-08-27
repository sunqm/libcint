#!/usr/bin/env python
# $Id$
# -*- coding: utf-8

'''
test libcint
'''

__author__ = "Qiming Sun <osirpt.sun@gmail.com>"

import sys
import os
import ctypes
import numpy

_cint = numpy.ctypeslib.load_library('libcint', os.path.abspath(os.path.join(__file__, '../../build')))


PTR_LIGHT_SPEED    = 0
PTR_COMMON_ORIG    = 1
PTR_SHIELDING_ORIG = 4
PTR_RINV_ORIG      = 4
PTR_RINV_ZETA      = 7
PTR_ENV_START      = 20

CHARGE_OF  = 0
PTR_COORD  = 1
NUC_MOD_OF = 2
PTR_ZETA   = 3
RAD_GRIDS  = 4
ANG_GRIDS  = 5
ATM_SLOTS  = 6

ATOM_OF   = 0
ANG_OF    = 1
NPRIM_OF  = 2
NCTR_OF   = 3
KAPPA_OF  = 4
PTR_EXP   = 5
PTR_COEFF = 6
BAS_SLOTS = 8

natm = 4
nbas = 0
atm = numpy.zeros((natm,ATM_SLOTS), dtype=numpy.int32)
bas = numpy.zeros((1000,BAS_SLOTS), dtype=numpy.int32)
env = numpy.zeros(10000)
off = PTR_ENV_START
for i in range(natm):
    atm[i, CHARGE_OF] = (i+1)*2
    atm[i, PTR_COORD] = off
    env[off+0] = .2 * (i+1)
    env[off+1] = .3 + (i+1) * .5
    env[off+2] = .1 - (i+1) * .5
    off += 3
off0 = off

# basis with kappa > 0
nh = 0

bas[nh,ATOM_OF ]  = 0
bas[nh,ANG_OF  ]  = 1
bas[nh,KAPPA_OF]  = 1
bas[nh,NPRIM_OF]  = 1
bas[nh,NCTR_OF ]  = 1
bas[nh,PTR_EXP]   = off
env[off+0] = 1
bas[nh,PTR_COEFF] = off + 1
env[off+1] = 1
off += 2
nh += 1

bas[nh,ATOM_OF ]  = 1
bas[nh,ANG_OF  ]  = 2
bas[nh,KAPPA_OF]  = 2
bas[nh,NPRIM_OF]  = 2
bas[nh,NCTR_OF ]  = 2
bas[nh,PTR_EXP]   = off
env[off+0] = 5
env[off+1] = 3
bas[nh,PTR_COEFF] = off + 2
env[off+2] = 1
env[off+3] = 2
env[off+4] = 4
env[off+5] = 1
off += 6
nh += 1

bas[nh,ATOM_OF ]  = 2
bas[nh,ANG_OF  ]  = 3
bas[nh,KAPPA_OF]  = 3
bas[nh,NPRIM_OF]  = 1
bas[nh,NCTR_OF ]  = 1
bas[nh,PTR_EXP ]  = off
env[off+0] = 1
bas[nh,PTR_COEFF] = off + 1
env[off+1] = 1
off += 2
nh += 1

bas[nh,ATOM_OF ]  = 3
bas[nh,ANG_OF  ]  = 4
bas[nh,KAPPA_OF]  = 4
bas[nh,NPRIM_OF]  = 1
bas[nh,NCTR_OF ]  = 1
bas[nh,PTR_EXP ]  = off
env[off+0] = .5
bas[nh,PTR_COEFF] = off + 1
env[off+1] = 1.
off = off + 2
nh += 1

nbas = nh

# basis with kappa < 0
n = off - off0
for i in range(n):
    env[off+i] = env[off0+i]

for i in range(nh):
        bas[i+nh,ATOM_OF ] = bas[i,ATOM_OF ]
        bas[i+nh,ANG_OF  ] = bas[i,ANG_OF  ] - 1
        bas[i+nh,KAPPA_OF] =-bas[i,KAPPA_OF]
        bas[i+nh,NPRIM_OF] = bas[i,NPRIM_OF]
        bas[i+nh,NCTR_OF ] = bas[i,NCTR_OF ]
        bas[i+nh,PTR_EXP ] = bas[i,PTR_EXP ]  + n
        bas[i+nh,PTR_COEFF]= bas[i,PTR_COEFF] + n
        env[bas[i+nh,PTR_COEFF]] /= 2 * env[bas[i,PTR_EXP]]

env[bas[5,PTR_COEFF]+0] = env[bas[1,PTR_COEFF]+0] / (2 * env[bas[1,PTR_EXP]+0])
env[bas[5,PTR_COEFF]+1] = env[bas[1,PTR_COEFF]+1] / (2 * env[bas[1,PTR_EXP]+1])
env[bas[5,PTR_COEFF]+2] = env[bas[1,PTR_COEFF]+2] / (2 * env[bas[1,PTR_EXP]+0])
env[bas[5,PTR_COEFF]+3] = env[bas[1,PTR_COEFF]+3] / (2 * env[bas[1,PTR_EXP]+1])

natm = ctypes.c_int(natm)
nbas = ctypes.c_int(nbas)
c_atm = atm.ctypes.data_as(ctypes.c_void_p)
c_bas = bas.ctypes.data_as(ctypes.c_void_p)
c_env = env.ctypes.data_as(ctypes.c_void_p)

opt = ctypes.POINTER(ctypes.c_void_p)()
_cint.CINTlen_spinor.restype = ctypes.c_int

try:
    import pyscf
except ImportError:
    pyscf = None

def test_int2c1e_sph():
    mol = pyscf.M()
    mol._atm = atm[:natm.value]
    mol._bas = bas[:nbas.value]
    mol._env = env
    coords = mol.atom_coords()
    ao = mol.eval_gto('GTOval_sph', coords)

    fnpp1 = _cint.cint1e_ipiprinv_sph
    fnp1p = _cint.cint1e_iprinvip_sph
    nullptr = ctypes.POINTER(ctypes.c_void_p)()
    def by_pp(shls, shape):
        buf = numpy.empty(shape+(9,), order='F')
        fnpp1(buf.ctypes.data_as(ctypes.c_void_p), (ctypes.c_int*4)(*shls),
              c_atm, natm, c_bas, nbas, c_env, nullptr)
        ref = buf[:,:,0] + buf[:,:,4] + buf[:,:,8]
        fnp1p(buf.ctypes.data_as(ctypes.c_void_p), (ctypes.c_int*4)(*shls),
              c_atm, natm, c_bas, nbas, c_env, nullptr)
        ref+=(buf[:,:,0] + buf[:,:,4] + buf[:,:,8])*2
        shls = (shls[1], shls[0])
        shape = (shape[1], shape[0]) + (9,)
        buf = numpy.empty(shape, order='F')
        fnpp1(buf.ctypes.data_as(ctypes.c_void_p), (ctypes.c_int*4)(*shls),
              c_atm, natm, c_bas, nbas, c_env, nullptr)
        ref+= (buf[:,:,0] + buf[:,:,4] + buf[:,:,8]).transpose(1,0)
        return ref * (-.25/numpy.pi)

    ao_loc = mol.ao_loc_nr()
    for nucid in range(mol.natm):
        mol.set_rinv_orig(coords[nucid])
        for j in range(nbas.value):
            j0 = ao_loc[j]
            j1 = ao_loc[j+1]
            for i in range(j+1):
                di = (bas[i,ANG_OF] * 2 + 1) * bas[i,NCTR_OF]
                dj = (bas[j,ANG_OF] * 2 + 1) * bas[j,NCTR_OF]
                shls = (i, j)
                i0 = ao_loc[i]
                i1 = ao_loc[i+1]
                buf = numpy.einsum('i,j->ij', ao[nucid,i0:i1], ao[nucid,j0:j1])
                ref = by_pp(shls, (di,dj))
                dd = abs(ref - buf).sum()
                if dd > 1e-8:
                    print("* FAIL: cint2c1e", "  shell:", i, j, "err:", dd)
                    return
    print('cint1e_ipiprinv_sph cint1e_iprinvip_sph pass')

def test_int4c1e_sph():
    fnpp1 = _cint.cint2e_ipip1_sph
    fnp1p = _cint.cint2e_ipvip1_sph
    nullptr = ctypes.POINTER(ctypes.c_void_p)()
    def by_pp(shls, shape):
        buf = numpy.empty(shape+(9,), order='F')
        fnpp1(buf.ctypes.data_as(ctypes.c_void_p), (ctypes.c_int*4)(*shls),
              c_atm, natm, c_bas, nbas, c_env, nullptr)
        ref = buf[:,:,:,:,0] + buf[:,:,:,:,4] + buf[:,:,:,:,8]
        fnp1p(buf.ctypes.data_as(ctypes.c_void_p), (ctypes.c_int*4)(*shls),
              c_atm, natm, c_bas, nbas, c_env, nullptr)
        ref+=(buf[:,:,:,:,0] + buf[:,:,:,:,4] + buf[:,:,:,:,8])*2
        shls = (shls[1], shls[0]) + shls[2:]
        shape = (shape[1], shape[0]) + shape[2:] + (9,)
        buf = numpy.empty(shape, order='F')
        fnpp1(buf.ctypes.data_as(ctypes.c_void_p), (ctypes.c_int*4)(*shls),
              c_atm, natm, c_bas, nbas, c_env, nullptr)
        ref+= (buf[:,:,:,:,0] + buf[:,:,:,:,4] + buf[:,:,:,:,8]).transpose(1,0,2,3)
        return ref * (-.25/numpy.pi)

    intor = _cint.cint4c1e_sph
    for l in range(nbas.value):
        for k in range(l+1):
            for j in range(nbas.value):
                for i in range(j+1):
                    di = (bas[i,ANG_OF] * 2 + 1) * bas[i,NCTR_OF]
                    dj = (bas[j,ANG_OF] * 2 + 1) * bas[j,NCTR_OF]
                    dk = (bas[k,ANG_OF] * 2 + 1) * bas[k,NCTR_OF]
                    dl = (bas[l,ANG_OF] * 2 + 1) * bas[l,NCTR_OF]
                    shls = (i, j, k, l)
                    buf = numpy.empty((di,dj,dk,dl), order='F')
                    intor(buf.ctypes.data_as(ctypes.c_void_p), (ctypes.c_int*4)(*shls),
                          c_atm, natm, c_bas, nbas, c_env, nullptr)
                    ref = by_pp(shls, (di,dj,dk,dl))
                    dd = abs(ref - buf).max()
                    if dd > 1e-6:
                        print("* FAIL: cint4c1e", "  shell:", i, j, k, l, "err:", dd)
                        return
    print('cint4c1e_sph pass')

if __name__ == '__main__':
    if pyscf:
        test_int2c1e_sph()
    test_int4c1e_sph()
