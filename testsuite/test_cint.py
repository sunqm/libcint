#!/usr/bin/env python
# $Id$
# -*- coding: utf-8
from __future__ import print_function

'''
test libcint
'''

__author__ = "Qiming Sun <osirpt.sun@gmail.com>"

import sys
import os
import ctypes
import numpy

_cint = numpy.ctypeslib.load_library('libcint', '.')


PTR_COMMON_ORIG    = 1
PTR_SHIELDING_ORIG = 4
PTR_RINV_ORIG      = 4
PTR_RINV_ZETA      = 7
PTR_RANGE_OMEGA    = 8
PTR_EXPCUTOFF      = 11
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


def close(v1, vref, count, place):
    return round(abs(v1-vref)/count, place) == 0

def test_int1e_sph(name, vref, dim, place):
    intor = getattr(_cint, name)
    intor.restype = ctypes.c_void_p
    op = (ctypes.c_double * (10000 * dim))()
    v1 = 0
    cnt = 0
    for j in range(nbas.value*2):
        for i in range(j+1):
            di = (bas[i,ANG_OF] * 2 + 1) * bas[i,NCTR_OF]
            dj = (bas[j,ANG_OF] * 2 + 1) * bas[j,NCTR_OF]
            shls = (ctypes.c_int * 2)(i, j)
            intor(op, shls, c_atm, natm, c_bas, nbas, c_env)
            v1 += abs(numpy.array(op[:di*dj*dim])).sum()
            cnt += di*dj*dim
    if close(v1, vref, cnt, place):
        print("pass: ", name)
    else:
        print("* FAIL: ", name, ". err:", '%.16g' % abs(v1-vref), "/", vref)

def cdouble_to_cmplx(arr):
    return numpy.array(arr)[0::2] + numpy.array(arr)[1::2] * 1j

def test_int1e_spinor(name, vref, dim, place):
    intor = getattr(_cint, name)
    intor.restype = ctypes.c_void_p
    op = (ctypes.c_double * (20000 * dim))()
    v1 = 0
    cnt = 0
    for j in range(nbas.value*2):
        for i in range(j+1):
            di = _cint.CINTlen_spinor(i, c_bas, nbas) * bas[i,NCTR_OF]
            dj = _cint.CINTlen_spinor(j, c_bas, nbas) * bas[j,NCTR_OF]
            shls = (ctypes.c_int * 2)(i, j)
            intor(op, shls, c_atm, natm, c_bas, nbas, c_env)
            v1 += abs(cdouble_to_cmplx(op[:di*dj*dim*2])).sum()
            cnt += di*dj*dim*2
    if close(v1, vref, cnt, place):
        print("pass: ", name)
    else:
        print("* FAIL: ", name, ". err:", '%.16g' % abs(v1-vref), "/", vref)

def max_loc(arr):
    loc = []
    maxi = arr.argmax()
    n = maxi
    for i in arr.shape:
        loc.append(n % i)
        n /= i
    loc.reverse()
    return maxi, loc

def test_comp1e_spinor(name1, name_ref, shift, dim, place):
    intor     = getattr(_cint, name1)
    intor.restype = ctypes.c_void_p
    intor_ref = getattr(_cint, name_ref)
    intor_ref.restype = ctypes.c_void_p
    op =     (ctypes.c_double * (20000 * dim))()
    op_ref = (ctypes.c_double * (20000 * dim))()

    pfac = 1
    if shift[0] > 0:
        pfac *= -1j
    if shift[1] > 0:
        pfac *= 1j

    for j in range(nbas.value*2 - shift[1]):
        for i in range(min(nbas.value*2-shift[0],j+1)):
            di = _cint.CINTlen_spinor(i, c_bas, nbas) * bas[i,NCTR_OF]
            dj = _cint.CINTlen_spinor(j, c_bas, nbas) * bas[j,NCTR_OF]
            shls = (ctypes.c_int * 2)(i+shift[0], j+shift[1])
            intor(op, shls, c_atm, natm, c_bas, nbas, c_env)
            shls_ref = (ctypes.c_int * 2)(i, j)
            intor_ref(op_ref, shls_ref, c_atm, natm, c_bas, nbas, c_env)
            dd = abs(pfac * cdouble_to_cmplx(op[:di*dj*dim*2]).reshape(di,dj,dim)
                     - cdouble_to_cmplx(op_ref[:di*dj*dim*2]).reshape(di,dj,dim))
            if numpy.round(dd, place).sum():
                maxi = dd.argmax()
                print("* FAIL: ", name1, "/", name_ref, ". shell:", i, j, \
                        "err:", dd.flatten()[maxi], \
                        "/", op_ref[maxi*2]+op_ref[maxi*2+1]*1j)
                return
    print("pass: ", name1, "/", name_ref)

####################
def test_int2e_sph(name, vref, dim, place):
    intor = getattr(_cint, name)
    intor.restype = ctypes.c_void_p
    op = (ctypes.c_double * (1000000 * dim))()
    v1 = 0
    cnt = 0
    for l in range(nbas.value*2):
        for k in range(l+1):
            for j in range(nbas.value*2):
                for i in range(j+1):
                    di = (bas[i,ANG_OF] * 2 + 1) * bas[i,NCTR_OF]
                    dj = (bas[j,ANG_OF] * 2 + 1) * bas[j,NCTR_OF]
                    dk = (bas[k,ANG_OF] * 2 + 1) * bas[k,NCTR_OF]
                    dl = (bas[l,ANG_OF] * 2 + 1) * bas[l,NCTR_OF]
                    shls = (ctypes.c_int * 4)(i, j, k, l)
                    intor(op, shls, c_atm, natm, c_bas, nbas, c_env, opt)
                    v1 += abs(numpy.array(op[:di*dj*dk*dl*dim])).sum()
                    cnt += di*dj*dk*dl*dim
    if close(v1, vref, cnt, place):
        print("pass: ", name)
    else:
        print("* FAIL: ", name, ". err:", '%.16g' % abs(v1-vref), "/", vref)

def test_erf(name, omega, place):
    intor = getattr(_cint, name)
    intor.restype = ctypes.c_void_p
    env_lr = env.copy()
    env_lr[PTR_RANGE_OMEGA] = omega
    env_sr = env.copy()
    env_sr[PTR_RANGE_OMEGA] = -omega
    c_env_lr = env_lr.ctypes.data_as(ctypes.c_void_p)
    c_env_sr = env_sr.ctypes.data_as(ctypes.c_void_p)

    op_lr = numpy.empty(10000)
    op_sr = numpy.empty(10000)
    op = numpy.empty(10000)
    v1 = 0
    cnt = 0
    shls = (ctypes.c_int * 4)(3,3,3,3)
    intor(op.ctypes.data_as(ctypes.c_void_p), shls, c_atm, natm, c_bas, nbas, c_env, opt)
    intor(op_lr.ctypes.data_as(ctypes.c_void_p), shls, c_atm, natm, c_bas, nbas, c_env_lr, opt)
    intor(op_sr.ctypes.data_as(ctypes.c_void_p), shls, c_atm, natm, c_bas, nbas, c_env_sr, opt)

    for l in range(nbas.value):
        for k in range(l+1):
            for j in range(l+1):
                for i in range(j+1):
                    di = (bas[i,ANG_OF] * 2 + 1) * bas[i,NCTR_OF]
                    dj = (bas[j,ANG_OF] * 2 + 1) * bas[j,NCTR_OF]
                    dk = (bas[k,ANG_OF] * 2 + 1) * bas[k,NCTR_OF]
                    dl = (bas[l,ANG_OF] * 2 + 1) * bas[l,NCTR_OF]
                    shls = (ctypes.c_int * 4)(i, j, k, l)
                    intor(op.ctypes.data_as(ctypes.c_void_p), shls, c_atm, natm, c_bas, nbas, c_env, opt)
                    intor(op_lr.ctypes.data_as(ctypes.c_void_p), shls, c_atm, natm, c_bas, nbas, c_env_lr, opt)
                    intor(op_sr.ctypes.data_as(ctypes.c_void_p), shls, c_atm, natm, c_bas, nbas, c_env_sr, opt)
                    cnt = di * dj * dk * dl
                    dd = abs(op_lr[:cnt] + op_sr[:cnt] - op[:cnt])

                    with numpy.errstate(invalid='ignore'):
                        dd[dd > 1e-6] / op[:cnt][dd > 1e-6]
                    if numpy.round(dd.max(), place) > 0:
                        maxi = dd.argmax()
                        print("* FAIL: erf+erfc", name, "omega=", omega,
                              " shell:", i, j, k, l,
                              "erf:", op_lr[maxi], "erfc:", op_sr[maxi],
                              "regular:", op[maxi], "err:", dd[maxi])
                        return
    print("pass: erf+erfc", name, "omega=", omega)

def test_int2e_spinor(name, vref, dim, place):
    intor = getattr(_cint, name)
    intor.restype = ctypes.c_void_p
    op = (ctypes.c_double * (2000000 * dim))()
    v1 = 0
    cnt = 0
    for l in range(nbas.value*2):
        for k in range(l+1):
            for j in range(nbas.value*2):
                for i in range(j+1):
                    di = _cint.CINTlen_spinor(i, c_bas, nbas) * bas[i,NCTR_OF]
                    dj = _cint.CINTlen_spinor(j, c_bas, nbas) * bas[j,NCTR_OF]
                    dk = _cint.CINTlen_spinor(k, c_bas, nbas) * bas[k,NCTR_OF]
                    dl = _cint.CINTlen_spinor(l, c_bas, nbas) * bas[l,NCTR_OF]
                    shls = (ctypes.c_int * 4)(i, j, k, l)
                    intor(op, shls, c_atm, natm, c_bas, nbas, c_env, opt)
                    v1 += abs(cdouble_to_cmplx(op[:di*dj*dk*dl*dim*2])).sum()
                    cnt += di*dj*dk*dl*dim*2
    if close(v1, vref, cnt, place):
        print("pass: ", name)
    else:
        print("* FAIL: ", name, ". err:", '%.16g' % abs(v1-vref), "/", vref)

def test_comp2e_spinor(name1, name_ref, shift, dim, place):
    intor     = getattr(_cint, name1)
    intor.restype = ctypes.c_void_p
    intor_ref = getattr(_cint, name_ref)
    intor_ref.restype = ctypes.c_void_p
    op =     (ctypes.c_double * (2000000 * dim))()
    op_ref = (ctypes.c_double * (2000000 * dim))()

    pfac = 1
    if shift[0] > 0:
        pfac *= -1j
    if shift[1] > 0:
        pfac *= 1j
    if shift[2] > 0:
        pfac *= -1j
    if shift[3] > 0:
        pfac *= 1j

    for l in range(nbas.value*2 - shift[3]):
        for k in range(min(nbas.value*2-shift[0],l+1)):
            for j in range(nbas.value*2 - shift[1]):
                for i in range(min(nbas.value*2-shift[0],j+1)):
                    di = _cint.CINTlen_spinor(i, c_bas, nbas) * bas[i,NCTR_OF]
                    dj = _cint.CINTlen_spinor(j, c_bas, nbas) * bas[j,NCTR_OF]
                    dk = _cint.CINTlen_spinor(k, c_bas, nbas) * bas[k,NCTR_OF]
                    dl = _cint.CINTlen_spinor(l, c_bas, nbas) * bas[l,NCTR_OF]
                    shls = (ctypes.c_int * 4)(i+shift[0], j+shift[1], k+shift[2], l+shift[3])
                    intor(op, shls, c_atm, natm, c_bas, nbas, c_env, opt)
                    shls_ref = (ctypes.c_int * 4)(i, j, k, l)
                    intor_ref(op_ref, shls_ref, c_atm, natm, c_bas, nbas, c_env, opt)
                    dd = abs(pfac * cdouble_to_cmplx(op[:di*dj*dk*dl*dim*2]).reshape(di,dj,dk,dl,dim)
                             - cdouble_to_cmplx(op_ref[:di*dj*dk*dl*dim*2]).reshape(di,dj,dk,dl,dim))
                    if numpy.round(dd, place).sum():
                        maxi = dd.argmax()
                        print("* FAIL: ", name1, "/", name_ref, ". shell:", i, j, k, l, \
                                "err:", dd.flatten()[maxi], \
                                "/", op_ref[maxi*2]+op_ref[maxi*2+1]*1j)
                        return
    print("pass: ", name1, "/", name_ref)



if __name__ == "__main__":
    if "--high-prec" in sys.argv:
        def close(v1, vref, count, place):
            return round(abs(v1-vref), place) == 0

    for f in (('cint1e_ovlp_sph'  , 320.9470780962389, 1, 11),
              ('cint1e_nuc_sph'   , 3664.898206036863, 1, 10),
              ('cint1e_kin_sph'   , 887.2525599069498, 1, 11),
              ('cint1e_ia01p_sph' , 210.475021425001 , 3, 12),
              ('cint1e_cg_irxp_sph', 3464.41761486531, 3, 10),
              ('cint1e_giao_irjxp_sph', 2529.89787038728, 3, 10),
              ('cint1e_igkin_sph' , 107.6417224161130, 3, 11),
              ('cint1e_igovlp_sph', 37.94968860771099, 3, 12),
              ('cint1e_ignuc_sph' , 478.5594827282386, 3, 11),
              ('cint1e_ipovlp_sph', 429.5284222008585, 3, 11),
              ('cint1e_ipkin_sph' , 1307.395170673386, 3, 10),
              ('cint1e_ipnuc_sph' , 8358.422626593954, 3, 10),
              ('cint1e_iprinv_sph', 385.1108471512923, 3, 11),
              ('cint1e_prinvxp_sph',210.475021425001, 3, 11),
              ('cint1e_z_sph'     , 651.8811101988866, 1, 11),
              ('cint1e_zz_sph'    , 1881.075059037941, 1, 10),
              ('cint1e_r_sph'     , 1803.13043674652 , 3, 10),
              ('cint1e_rr_sph'    , 13379.47937680471, 9,  9),
              ('cint1e_r2_sph'    , 5237.899221349136, 1, 10),
             ):
        test_int1e_sph(*f)

    for f in (('cint1e_ovlp'     , 284.4528456839759, 1, 11),
              ('cint1e_nuc'      , 3025.766689838620, 1, 10),
              ('cint1e_gnuc'     , 296.6160944673867, 3, 11),
              ('cint1e_srsr'     , 1430.424389624617, 1, 10),
              ('cint1e_sr'       , 240.4064385362524, 1, 11),
              ('cint1e_srsp'     , 1022.805155947573, 1, 10),
              ('cint1e_spsp'     , 1554.251610462129, 1, 10),
              ('cint1e_sp'       , 265.2448605537422, 1, 11),
              ('cint1e_spspsp'   , 1551.856568558924, 1, 10),
              ('cint1e_spnuc'    , 3905.024862120781, 1, 10),
              ('cint1e_spnucsp'  , 20689.33926165072, 1, 9 ),
              ('cint1e_srnucsr'  , 13408.06084488522, 1, 9 ),
              ('cint1e_cg_sa10sa01', 319.6545034966355, 9, 11),
              ('cint1e_cg_sa10sp'  , 1705.563585675829, 3, 10),
              ('cint1e_cg_sa10nucsp',16506.04502697362, 3, 9 ),
              ('cint1e_giao_sa10sa01' , 358.7833729392868, 9, 11),
              ('cint1e_giao_sa10sp'   , 1070.550400465705, 3, 10),
              ('cint1e_giao_sa10nucsp', 12819.05472701636, 3, 9 ),
              ('cint1e_govlp'    , 23.2674074483772, 3, 12),
              ('cint1e_sa01sp'   , 218.244203172625, 3, 11),
              ('cint1e_spgsp'    , 96.9346217249868, 3, 12),
              ('cint1e_spgnucsp' , 1659.37670007911, 3, 10),
              ('cint1e_spgsa01'  , 37.8884662927634, 9, 12),
              ('cint1e_ipovlp'   , 153.860148521121, 3, 12),
              ('cint1e_ipkin'    , 497.249399637873, 3, 11),
              ('cint1e_ipnuc'    , 4506.61348255897, 3, 10),
              ('cint1e_iprinv'   , 240.036283917245, 3, 11),
              ('cint1e_ipspnucsp', 35059.4071347107, 3,  9),
              ('cint1e_ipsprinvsp',1166.20850563398, 3, 10),
             ):
        test_int1e_spinor(*f)

    for f in (# rys_roots for i,j,k,l=3,3,3,3 has round-off error ~ 1e-5
              ('cint2e_sph'    , 56243.88328768107 , 1, 8 ),
              ('cint2e_ip1_sph', 115489.8643866550 , 3, 8 ),
             ):
        test_int2e_sph(*f)

    test_erf('cint2e_sph', 0.2, 9)
    test_erf('cint2e_sph', 0.5, 9)
    test_erf('cint2e_sph', 0.8, 9)

    if "--quick" not in sys.argv:
        # Four tests marked with "# *" may fail in quadmath mode
        for f in (('cint2e_ip1_sph', 115489.8643866550 , 3, 8 ),
                  ('cint2e_p1vxp1_sph', 89014.88169743448, 3, 9),
                 ):
            test_int2e_sph(*f)
        for f in (('cint2e'             , 37737.11365710611, 1, 8),
                  ('cint2e_spsp1'       , 221528.4764668166, 1, 8),
                  ('cint2e_spsp1spsp2'  , 1391716.876869147, 1, 7),  # *
                  ('cint2e_srsr1'       , 178572.7398308939, 1, 8),
                  ('cint2e_srsr1srsr2'  , 860883.6288270953, 1, 8),  # *
                  ('cint2e_cg_sa10sp1'  , 241519.2143647713, 3, 8),
                  ('cint2e_cg_sa10sp1spsp2'  , 1419443.469767018, 3, 7),  # *
                  ('cint2e_giao_sa10sp1'     , 153861.920807804 , 3, 8),
                  ('cint2e_giao_sa10sp1spsp2', 918284.9464686266, 3, 8),  # *
                  ('cint2e_g1'          , 3755.251591892025, 3, 10),
                  ('cint2e_spgsp1'      , 16626.99103794526, 3, 9 ),
                  ('cint2e_g1spsp2'     , 22186.56654833549, 3, 9 ),
                  ('cint2e_spgsp1spsp2' , 107110.2340526177, 3, 8 ),
                  ('cint2e_ip1'         , 34912.85433806438, 3, 9 ),
                  ('cint2e_ipspsp1'     , 221092.5556043494, 3, 8 ),
                  ('cint2e_ip1spsp2'    , 212447.1029358293, 3, 8 ),
                  ('cint2e_ipspsp1spsp2', 1443972.936563201, 3, 7 ),
                 ):
            test_int2e_spinor(*f)

    test_comp2e_spinor('cint2e_spsp1', 'cint2e', (4,4,0,0), 1, 11)
    test_comp2e_spinor('cint2e_spsp1spsp2', 'cint2e', (4,4,4,4), 1, 11)
    test_comp2e_spinor('cint2e_spsp1spsp2', 'cint2e_spsp1', (0,0,4,4), 1, 11)
    test_comp2e_spinor('cint2e_spgsp1', 'cint2e_g1', (4,4,0,0), 3, 11)
    test_comp2e_spinor('cint2e_g1spsp2', 'cint2e_g1', (0,0,4,4), 3, 11)
    test_comp2e_spinor('cint2e_spgsp1spsp2', 'cint2e_g1', (4,4,4,4), 3, 11)
    test_comp2e_spinor('cint2e_ipspsp1', 'cint2e_ip1', (4,4,0,0), 3, 11)
    test_comp2e_spinor('cint2e_ip1spsp2', 'cint2e_ip1', (0,0,4,4), 3, 11)
    test_comp2e_spinor('cint2e_ipspsp1spsp2', 'cint2e_ip1', (4,4,4,4), 3, 11)

    fz  = getattr(_cint, 'cint1e_z_sph')
    fzz = getattr(_cint, 'cint1e_zz_sph')
    fr  = getattr(_cint, 'cint1e_r_sph')
    fr2 = getattr(_cint, 'cint1e_r2_sph')
    frr = getattr(_cint, 'cint1e_rr_sph')
    v1 = 0
    for j in range(nbas.value*2):
        for i in range(j+1):
            di = (bas[i,ANG_OF] * 2 + 1) * bas[i,NCTR_OF]
            dj = (bas[j,ANG_OF] * 2 + 1) * bas[j,NCTR_OF]
            opz  = numpy.empty((di,dj)  , order='F')
            opzz = numpy.empty((di,dj)  , order='F')
            opr  = numpy.empty((di,dj,3), order='F')
            opr2 = numpy.empty((di,dj)  , order='F')
            oprr = numpy.empty((di,dj,9), order='F')
            shls = (ctypes.c_int * 2)(i, j)
            fz ( opz.ctypes.data_as(ctypes.c_void_p), shls, c_atm, natm, c_bas, nbas, c_env)
            fzz(opzz.ctypes.data_as(ctypes.c_void_p), shls, c_atm, natm, c_bas, nbas, c_env)
            fr ( opr.ctypes.data_as(ctypes.c_void_p), shls, c_atm, natm, c_bas, nbas, c_env)
            fr2(opr2.ctypes.data_as(ctypes.c_void_p), shls, c_atm, natm, c_bas, nbas, c_env)
            frr(oprr.ctypes.data_as(ctypes.c_void_p), shls, c_atm, natm, c_bas, nbas, c_env)
            v1 = abs(opz-opr[:,:,2]).sum()
            v1 += abs(opzz-oprr[:,:,8]).sum()
            v1 += abs(opr2-oprr[:,:,0]-oprr[:,:,4]-oprr[:,:,8]).sum()
            if round(v1/(di*dj), 13):
                print("* FAIL: ", i, j, v1)
