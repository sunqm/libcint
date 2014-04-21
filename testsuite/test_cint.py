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

if sys.platform == 'darwin':
    alib = os.environ['buildir'] + '/testsuite/.libs/libtestcint.dylib'
else:
    alib = os.environ['buildir'] + '/testsuite/.libs/libtestcint.so'
_cint = ctypes.cdll.LoadLibrary(alib)

PTR_LIGHT_SPEED    = 0
PTR_COMMON_ORIG    = 1
PTR_SHIELDING_ORIG = 4
PTR_RINV_ORIG      = 4
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
BAS_SLOTS = 8

atm = (ctypes.c_int * 200)()
bas = (ctypes.c_int * 200)()
env = (ctypes.c_double * 400)()
natm = ctypes.c_int()
nbas = ctypes.c_int()
_cint.init_test_env(atm, ctypes.addressof(natm), bas, ctypes.addressof(nbas), env)

_cint.CINTlen_spinor.restype = ctypes.c_int

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
        print "* FAIL: ", name, ". err:", '%.16g' % abs(v1-vref), "/", vref
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
            di = _cint.CINTlen_spinor(i, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*i]
            dj = _cint.CINTlen_spinor(j, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*j]
            shls = (ctypes.c_int * 2)(i, j)
            intor(op, shls, atm, natm, bas, nbas, env);
            v1 += abs(cdouble_to_cmplx(op[:di*dj*dim*2])).sum()
    if round(abs(v1-vref), place):
        print "* FAIL: ", name, ". err:", '%.16g' % abs(v1-vref), "/", vref
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
            di = _cint.CINTlen_spinor(i, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*i]
            dj = _cint.CINTlen_spinor(j, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*j]
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
        print "* FAIL: ", name, ". err:", '%.16g' % abs(v1-vref), "/", vref
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
                    di = _cint.CINTlen_spinor(i, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*i]
                    dj = _cint.CINTlen_spinor(j, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*j]
                    dk = _cint.CINTlen_spinor(k, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*k]
                    dl = _cint.CINTlen_spinor(l, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*l]
                    shls = (ctypes.c_int * 4)(i, j, k, l)
                    intor(op, shls, atm, natm, bas, nbas, env);
                    v1 += abs(cdouble_to_cmplx(op[:di*dj*dk*dl*dim*2])).sum()
    if round(abs(v1-vref), place):
        print "* FAIL: ", name, ". err:", '%.16g' % abs(v1-vref), "/", vref
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
                    di = _cint.CINTlen_spinor(i, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*i]
                    dj = _cint.CINTlen_spinor(j, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*j]
                    dk = _cint.CINTlen_spinor(k, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*k]
                    dl = _cint.CINTlen_spinor(l, bas, nbas) * bas[NCTR_OF+BAS_SLOTS*l]
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
              ('cint1e_cg_irxp_sph', 3464.41761486531, 3, 10),
              ('cint1e_giao_irjxp_sph', 2529.89787038728, 3, 10),
              ('cint1e_igkin_sph' , 107.6417224161130, 3, 11),
              ('cint1e_igovlp_sph', 37.94968860771099, 3, 12),
              ('cint1e_ignuc_sph' , 478.5594827282386, 3, 11),
              ('cint1e_ipovlp_sph', 429.5284222008585, 3, 11),
              ('cint1e_ipkin_sph' , 1307.395170673386, 3, 10),
              ('cint1e_ipnuc_sph' , 8358.422626593954, 3, 10),
              ('cint1e_iprinv_sph', 385.1108471512923, 3, 11),
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
              ('cint1e_ipspnucsp', 35059.4071347107, 3, 10),
              ('cint1e_ipsprinvsp',1166.20850563398, 3, 11),
             ):
        test_int1e_spinor(*f)

    for f in (('cint2e_sph'    , 56243.88328820552, 1, 10),
              ('cint2e_ig1_sph', 8101.087330632488, 3, 10),
              ('cint2e_ip1_sph', 115489.8643757788, 3, 9 ),
             ):
        test_int2e_sph(*f)

    for f in (('cint2e'             , 37737.11365741726, 1, 10),
              ('cint2e_spsp1'       , 221528.1553163525, 1, 10),
              ('cint2e_spsp1spsp2'  , 1391536.665233259, 1, 10),
              ('cint2e_srsr1'       , 178572.4800735003, 1, 10),
              ('cint2e_srsr1srsr2'  , 860772.1410876940, 1, 10),
              ('cint2e_cg_sa10sp1'  , 241519.0063952794, 3, 10),
              ('cint2e_cg_sa10sp1spsp2'  , 1419332.598373288, 3, 9 ),
              ('cint2e_giao_sa10sp1'     , 153861.7128168866, 3, 10),
              ('cint2e_giao_sa10sp1spsp2', 918176.1514756213, 3, 10),
              ('cint2e_g1'          , 3755.251590409050, 3, 10),
              ('cint2e_spgsp1'      , 16626.92238286884, 3, 10),
              ('cint2e_g1spsp2'     , 22186.46611681335, 3, 10),
              ('cint2e_spgsp1spsp2' , 107096.5206630970, 3, 10),
              ('cint2e_ip1'         , 34912.85433409394, 3, 10),
              ('cint2e_ipspsp1'     , 221091.4210207796, 3, 10),
              ('cint2e_ip1spsp2'    , 212446.0411584463, 3, 10),
              ('cint2e_ipspsp1spsp2', 1443719.034699620, 3, 9 ),
             ):
        test_int2e_spinor(*f)
