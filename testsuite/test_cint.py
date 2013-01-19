#!/usr/bin/env python
# $Id$
# -*- coding: utf-8

'''
test libcint
'''

__author__ = "Qiming Sun <osirpt.sun@gmail.com>"
__version__ = "$ 0.1 $"

import os
import ctypes
import numpy

alib = os.environ['buildir'] + '/testsuite/.libs/libtestcint.so'
_cint = ctypes.cdll.LoadLibrary(alib)

PTR_LIGHT_SPEED    = 0
PTR_COMMON_ORIG    = 1
PTR_SHIELDING_ORIG = 4
PTR_RINV_ORIG      = 4
PTR_AO_GAUGE       = 7
PTR_ENV_START      = 20

CHARGE_OF  = 0
PTR_COORD  = 1
NUC_MOD_OF = 2
PTR_MASS   = 3
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
GAUGE_OF  = 7
BAS_SLOTS = 8

atm = (ctypes.c_int * 200)()
bas = (ctypes.c_int * 200)()
env = (ctypes.c_double * 400)()
natm = ctypes.c_int()
nbas = ctypes.c_int()
_cint.init_test_env(atm, ctypes.addressof(natm), bas, ctypes.addressof(nbas), env)

_cint.len_spinor.restype = ctypes.c_int

def test_int1e_sph(name, vref, dim, place):
    intor = getattr(_cint, name)
    op = (ctypes.c_double * (10000 * dim))()
    v1 = 0
    for j in range(nbas.value*2):
        for i in range(j+1):
            di = (bas[ANG_OF+BAS_SLOTS*i] * 2 + 1) * bas[NCTR_OF+BAS_SLOTS*i]
            dj = (bas[ANG_OF+BAS_SLOTS*j] * 2 + 1) * bas[NCTR_OF+BAS_SLOTS*j]
            shls = (ctypes.c_int * 2)(i, j)
            intor(op, shls, atm, natm, bas, nbas, env);
            v1 += abs(numpy.array(op[:di*dj*dim])).sum()
    if round(abs(v1-vref), place):
        print "* FAIL: ", name, ". err:", abs(v1-vref), "/", vref
    else:
        print "pass: ", name

def cdouble_to_cmplx(arr):
    return numpy.array(arr)[0::2] + numpy.array(arr)[1::2] * 1j

def test_int1e_spinor(name, vref, dim, place):
    intor = getattr(_cint, name)
    op = (ctypes.c_double * (20000 * dim))()
    v1 = 0
    for j in range(nbas.value*2):
        for i in range(j+1):
            di = _cint.len_spinor(i, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*i]
            dj = _cint.len_spinor(j, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*j]
            shls = (ctypes.c_int * 2)(i, j)
            intor(op, shls, atm, natm, bas, nbas, env);
            v1 += abs(cdouble_to_cmplx(op[:di*dj*dim*2])).sum()
    if round(abs(v1-vref), place):
        print "* FAIL: ", name, ". err:", abs(v1-vref), "/", vref
    else:
        print "pass: ", name

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
    intor_ref = getattr(_cint, name_ref)
    op =     (ctypes.c_double * (20000 * dim))()
    op_ref = (ctypes.c_double * (20000 * dim))()

    pfac = 1
    if shift[0] > 0:
        pfac *= -1j
    if shift[1] > 0:
        pfac *= 1j

    for j in range(nbas.value*2 - shift[1]):
        for i in range(min(nbas.value*2-shift[0],j+1)):
            di = _cint.len_spinor(i, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*i]
            dj = _cint.len_spinor(j, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*j]
            shls = (ctypes.c_int * 2)(i+shift[0], j+shift[1])
            intor(op, shls, atm, natm, bas, nbas, env);
            shls_ref = (ctypes.c_int * 2)(i, j)
            intor_ref(op_ref, shls_ref, atm, natm, bas, nbas, env);
            dd = abs(pfac * cdouble_to_cmplx(op[:di*dj*dim*2]).reshape(di,dj,dim)
                     - cdouble_to_cmplx(op_ref[:di*dj*dim*2]).reshape(di,dj,dim))
            if numpy.round(dd, place).sum():
                maxi = dd.argmax()
                print "* FAIL: ", name1, "/", name_ref, ". shell:", i, j, \
                        "err:", dd.flatten()[maxi], \
                        "/", op_ref[maxi*2]+op_ref[maxi*2+1]*1j
                return
    print "pass: ", name1, "/", name_ref

####################
def test_int2e_sph(name, vref, dim, place):
    intor = getattr(_cint, name)
    op = (ctypes.c_double * (1000000 * dim))()
    v1 = 0
    for l in range(nbas.value*2):
        for k in range(l+1):
            for j in range(nbas.value*2):
                for i in range(j+1):
                    di = (bas[ANG_OF+BAS_SLOTS*i] * 2 + 1) * bas[NCTR_OF+BAS_SLOTS*i]
                    dj = (bas[ANG_OF+BAS_SLOTS*j] * 2 + 1) * bas[NCTR_OF+BAS_SLOTS*j]
                    dk = (bas[ANG_OF+BAS_SLOTS*k] * 2 + 1) * bas[NCTR_OF+BAS_SLOTS*k]
                    dl = (bas[ANG_OF+BAS_SLOTS*l] * 2 + 1) * bas[NCTR_OF+BAS_SLOTS*l]
                    shls = (ctypes.c_int * 4)(i, j, k, l)
                    intor(op, shls, atm, natm, bas, nbas, env);
                    v1 += abs(numpy.array(op[:di*dj*dk*dl*dim])).sum()
    if round(abs(v1-vref), place):
        print "* FAIL: ", name, ". err:", abs(v1-vref), "/", vref
    else:
        print "pass: ", name

def test_int2e_spinor(name, vref, dim, place):
    intor = getattr(_cint, name)
    op = (ctypes.c_double * (2000000 * dim))()
    v1 = 0
    for l in range(nbas.value*2):
        for k in range(l+1):
            for j in range(nbas.value*2):
                for i in range(j+1):
                    di = _cint.len_spinor(i, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*i]
                    dj = _cint.len_spinor(j, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*j]
                    dk = _cint.len_spinor(k, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*k]
                    dl = _cint.len_spinor(l, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*l]
                    shls = (ctypes.c_int * 4)(i, j, k, l)
                    intor(op, shls, atm, natm, bas, nbas, env);
                    v1 += abs(cdouble_to_cmplx(op[:di*dj*dk*dl*dim*2])).sum()
    if round(abs(v1-vref), place):
        print "* FAIL: ", name, ". err:", abs(v1-vref), "/", vref
    else:
        print "pass: ", name

def test_comp2e_spinor(name1, name_ref, shift, dim, place):
    intor     = getattr(_cint, name1)
    intor_ref = getattr(_cint, name_ref)
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
                    di = _cint.len_spinor(i, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*i]
                    dj = _cint.len_spinor(j, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*j]
                    dk = _cint.len_spinor(k, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*k]
                    dl = _cint.len_spinor(l, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*l]
                    shls = (ctypes.c_int * 4)(i+shift[0], j+shift[1], k+shift[2], l+shift[3])
                    intor(op, shls, atm, natm, bas, nbas, env);
                    shls_ref = (ctypes.c_int * 4)(i, j, k, l)
                    intor_ref(op_ref, shls_ref, atm, natm, bas, nbas, env);
                    dd = abs(pfac * cdouble_to_cmplx(op[:di*dj*dk*dl*dim*2]).reshape(di,dj,dk,dl,dim)
                             - cdouble_to_cmplx(op_ref[:di*dj*dk*dl*dim*2]).reshape(di,dj,dk,dl,dim))
                    if numpy.round(dd, place).sum():
                        maxi = dd.argmax()
                        print "* FAIL: ", name1, "/", name_ref, ". shell:", i, j, k, l, \
                                "err:", dd.flatten()[maxi], \
                                "/", op_ref[maxi*2]+op_ref[maxi*2+1]*1j
                        return
    print "pass: ", name1, "/", name_ref



if __name__ == "__main__":
    for f in (('cint1e_ovlp_sph'  , 320.9470780962389, 1, 11),
              ('cint1e_nuc_sph'   , 3664.898206036863, 1, 10),
              ('cint1e_kin_sph'   , 887.2525599069498, 1, 11),
              ('cint1e_ia01p_sph' , 210.475021425001 , 3, 12),
              ('cint1e_irxp_sph'  , 3464.41761486531 , 3, 10),
              ('cint1e_iking_sph' , 105.6144492946662, 3, 11),
              ('cint1e_iovlpg_sph', 37.94968860771099, 3, 12),
              ('cint1e_inucg_sph' , 478.5594827282386, 3, 11),
             ):
        test_int1e_sph(*f)

    for f in (('cint1e_ovlp'     , 284.4528456839759, 1, 11),
              ('cint1e_nuc'      , 3025.766689838620, 1, 10),
              ('cint1e_nucg'     , 296.6160944673867, 3, 11),
              ('cint1e_srsr'     , 1430.424389624617, 1, 10),
              ('cint1e_sr'       , 240.4064385362524, 1, 11),
              ('cint1e_srsp'     , 1022.805155947573, 1, 10),
              ('cint1e_spsp'     , 1554.251610462129, 1, 10),
              ('cint1e_sp'       , 265.2448605537422, 1, 11),
              ('cint1e_spspsp'   , 1551.856568558924, 1, 10),
              ('cint1e_spnuc'    , 3905.024862120781, 1, 10),
              ('cint1e_spnucsp'  , 20689.33926165072, 1, 9 ),
              ('cint1e_srnucsr'  , 13408.06084488522, 1, 9 ),
              ('cint1e_sa10sa01' , 319.6545034966355, 9, 11),
              ('cint1e_ovlpg'    , 23.26740744837723, 3, 12),
              ('cint1e_sa10sp'   , 1705.563585675829, 3, 10),
              ('cint1e_sa10nucsp', 16506.04502697362, 3, 9 ),
              ('cint1e_sa01sp'   , 218.2442031726256, 3, 11),
              ('cint1e_spgsp'    , 96.9346217249868 , 3, 12),
              ('cint1e_spgnucsp' , 1659.37670007911 , 3, 10),
              ('cint1e_spgsa01'  , 37.88846629276348, 9, 12),
             ):
        test_int1e_spinor(*f)

    for f in (('cint2e_sph'    , 56243.883287681034, 1, 8),
              ('cint2e_ig1_sph', 8101.0873343981893, 3, 9),
             ):
        test_int2e_sph(*f)

    for f in (('cint2e'             , 37737.1136571061, 1, 9 ),
              ('cint2e_spsp1'       , 221528.155316353, 1, 8 ),
              ('cint2e_spsp1spsp2'  , 1391536.66523326, 1, 7 ),
              ('cint2e_srsr1'       , 178572.480073500, 1, 8 ),
              ('cint2e_srsr1srsr2'  , 860772.141087694, 1, 8 ),
              ('cint2e_sa10sp1'     , 241519.006395279, 3, 8 ),
              ('cint2e_sa10sp1spsp2', 1419332.59837329, 3, 7 ),
              ('cint2e_g1'          , 3755.25159040905, 3, 10),
              ('cint2e_spgsp1'      , 16626.9223828688, 3, 9 ),
              ('cint2e_g1spsp2'     , 22186.4661168133, 3, 9 ),
              ('cint2e_spgsp1spsp2' , 107096.523556956, 3, 8 ),
             ):
        test_int2e_spinor(*f)
