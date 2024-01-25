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

_cint = ctypes.CDLL(os.path.abspath(os.path.join(__file__, '../../build/libcint.so')))

PTR_EXPCUTOFF      = 0
PTR_COMMON_ORIG    = 1
PTR_SHIELDING_ORIG = 4
PTR_RINV_ORIG      = 4
PTR_RINV_ZETA      = 7
PTR_RANGE_OMEGA    = 8
PTR_ENV_START      = 20
NGRIDS             = 11
PTR_GRIDS          = 12

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


def test_int1e_grids_sph1(name, vref, dim, place):
    intor = getattr(_cint, name)
    intor.restype = ctypes.c_void_p
    ngrids = 1
    numpy.random.seed(12)
    grids = numpy.random.random((ngrids, 3)) - 5.2
    env_g = numpy.append(env, grids.ravel())
    env_g[NGRIDS] = ngrids
    env_g[PTR_GRIDS] = env.size
    op = numpy.empty(1000000 * dim)
    v1 = 0
    cnt = 0
    for j in range(0,nbas.value*2):
        for i in range(j+1):
            di = (bas[i,ANG_OF] * 2 + 1) * bas[i,NCTR_OF]
            dj = (bas[j,ANG_OF] * 2 + 1) * bas[j,NCTR_OF]
            shls = (ctypes.c_int * 4)(i, j, 0, ngrids)
            intor(op.ctypes.data_as(ctypes.c_void_p),
                  shls, c_atm, natm, c_bas, nbas,
                  env_g.ctypes.data_as(ctypes.c_void_p), opt)
            v1 += abs(op[:ngrids*di*dj*dim]).sum()
            cnt += ngrids*di*dj*dim
            ref = numpy.empty((di, dj), order='F')
            env_g[PTR_RINV_ORIG:PTR_RINV_ORIG+3] = grids[0]
            _cint.cint1e_rinv_sph(ref.ctypes.data_as(ctypes.c_void_p),
                                  shls, c_atm, natm, c_bas, nbas,
                                  env_g.ctypes.data_as(ctypes.c_void_p), opt)
            dat = op[:ref.size].reshape(dj, di).T
            if abs(dat - ref).max() > 1e-12:
                print(i, j, abs(dat - ref).max())
    print('sum', v1, 'diff', v1 - vref)
    if abs(v1 - vref) < 1e-10:
        print('pass')
    else:
        print('failed')

def test_int1e_grids_sph(name, vref, dim, place):
    intor = getattr(_cint, name)
    intor.restype = ctypes.c_void_p
    ngrids = 148
    numpy.random.seed(12)
    grids = numpy.random.random((ngrids, 3)) - 5.2
    env_g = numpy.append(env, grids.ravel())
    env_g[NGRIDS] = ngrids
    env_g[PTR_GRIDS] = env.size
    op = numpy.empty(1000000 * dim)
    v1 = 0
    cnt = 0
    for j in range(nbas.value*2):
        for i in range(j+1):
            di = (bas[i,ANG_OF] * 2 + 1) * bas[i,NCTR_OF]
            dj = (bas[j,ANG_OF] * 2 + 1) * bas[j,NCTR_OF]
            shls = (ctypes.c_int * 4)(i, j, 0, ngrids)
            intor(op.ctypes.data_as(ctypes.c_void_p),
                  shls, c_atm, natm, c_bas, nbas,
                  env_g.ctypes.data_as(ctypes.c_void_p), opt)
            v1 += abs(op[:ngrids*di*dj*dim]).sum()
            cnt += ngrids*di*dj*dim
    print('sum', v1, 'diff', v1 - vref)
    if abs(v1 - vref) < 1e-10:
        print('pass')
    else:
        print('failed')

def test_mol1():
    import time
    import pyscf
    from pyscf import df
    mol0 = pyscf.M(atom='''
C          0.16146       -4.47867        0.00000
H         -0.89492        5.29077        0.00000
H          0.47278        4.59602        0.88488
H          0.47278        4.59602       -0.88488
H         -1.49009        3.01810       -0.87995
            ''', basis='ccpvtz')
    numpy.random.seed(12)
    ngrids = 201
    grids = numpy.random.random((ngrids, 3)) * 12 - 5

    def check(intor, ref):
        j3c = mol.intor(intor, grids=grids)
        diff = abs(j3c - ref).max()
        print(diff)
        return diff > 1e-12

    failed = False
    for omega in (0, 0.1, -0.1):
        for zeta in (0, 10, 1e16):
            print('omega, zeta', omega, zeta)
            if zeta == 0:
                expnt = 1e16
            else:
                expnt = zeta
            mol = mol0.copy()
            mol.omega = omega
            mol.set_rinv_zeta(zeta)
            fmol = pyscf.gto.fakemol_for_charges(grids, expnt)
            ref = df.incore.aux_e2(mol, fmol, intor='int3c2e').transpose(2,0,1)
            failed = check('int1e_grids', ref) or failed

            ref = df.incore.aux_e2(mol, fmol, intor='int3c2e_ip1').transpose(0,3,1,2)
            failed = check('int1e_grids_ip', ref) or failed

            ref = df.incore.aux_e2(mol, fmol, intor='int3c2e_ip1_cart').transpose(0,3,1,2)
            failed = check('int1e_grids_ip_cart', ref) or failed

            ref = df.incore.aux_e2(mol, fmol, intor='int3c2e_ip1_spinor').transpose(0,3,1,2)
            failed = check('int1e_grids_ip_spinor', ref) or failed

            ref = df.r_incore.aux_e2(mol, fmol, intor='int3c2e_spsp1_spinor').transpose(2,0,1)
            failed = check('int1e_grids_spvsp_spinor', ref) or failed
    if failed:
        print('failed')
    else:
        print('pass')

test_int1e_grids_sph1('cint1e_grids_sph', 36.81452996003706, 1, 9)
test_int1e_grids_sph('cint1e_grids_ip_sph', 3279.92861109671, 1, 9)
try:
    test_mol1()
except ImportError:
    pass
