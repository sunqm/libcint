import sys
import os
import ctypes
import numpy

_cint = numpy.ctypeslib.load_library('libcint', os.path.abspath(os.path.join(__file__, '../../build')))


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

def fp(a):
    return numpy.dot(a.ravel(), numpy.cos(numpy.arange(a.size)))

ng = 11
def random_cart(bas_id, seed=12, dtype='d'):
    l = bas[bas_id, ANG_OF]
    ncart = (l + 1) * (l + 2) // 2
    numpy.random.seed(seed)
    cart = numpy.random.random((ng, ncart))
    if dtype == 'z':
        cart = cart + numpy.random.random((ng, ncart)) * 1j
    return numpy.asarray(cart, order='F')

def test_c2s_ket_sph1(name, ref, thr=1e-12):
    v1 = 0
    fn = getattr(_cint, name)
    for j in range(nbas.value*2):
        l = bas[j, ANG_OF]
        dj = (l * 2 + 1)
        cart = random_cart(j)
        sph = numpy.empty((ng, dj), order='F')
        _cint.CINTc2s_ket_sph1(sph.ctypes.data_as(ctypes.c_void_p), cart.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(ng), ctypes.c_int(ng), ctypes.c_int(l))
        v1 += fp(sph)
    assert abs(ref - v1) < thr

def test_c2s_bra_spinor_e1sf(name, ref, thr=1e-12):
    kappa = 0
    v1 = 0
    fn = getattr(_cint, name)
    for j in range(nbas.value*2):
        l = bas[j, ANG_OF]
        dj = (l * 4 + 2)
        cart = numpy.asarray(random_cart(j).T, order='F')
        gsp = numpy.empty((dj, ng, 2), order='F', dtype=numpy.complex128)
        fn(gsp.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(ng), cart.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(kappa), ctypes.c_int(l))
        v1 += fp(gsp)
    assert abs(ref - abs(v1)) < thr

def test_c2s_bra_spinor_sf(name, ref, thr=1e-12):
    kappa = 0
    v1 = 0
    fn = getattr(_cint, name)
    for j in range(nbas.value*2):
        l = bas[j, ANG_OF]
        dj = (l * 4 + 2)
        cart = numpy.asarray(random_cart(j, dtype='z').T, order='F')
        gsp = numpy.empty((dj, ng, 2), order='F', dtype=numpy.complex128)
        fn(gsp.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(ng), cart.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(kappa), ctypes.c_int(l))
        v1 += fp(gsp)
    assert abs(ref - abs(v1)) < thr

def test_c2s_bra_spinor_si(name, ref, thr=1e-12):
    kappa = 0
    v1 = 0
    fn = getattr(_cint, name)
    for j in range(nbas.value*2):
        l = bas[j, ANG_OF]
        dj = (l * 4 + 2)
        cart = numpy.hstack((random_cart(j, dtype='z').T,
                             random_cart(j, seed=13, dtype='z').T))
        cart = numpy.asarray(cart, order='F')
        gsp = numpy.empty((dj, ng), order='F', dtype=numpy.complex128)
        fn(gsp.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(ng), cart.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(kappa), ctypes.c_int(l))
        v1 += fp(gsp)
    assert abs(ref - abs(v1)) < thr

def test_c2s_ket_spinor(name, ref, thr=1e-12):
    kappa = 0
    v1 = 0
    fn = getattr(_cint, name)
    for j in range(nbas.value*2):
        l = bas[j, ANG_OF]
        dj = (l * 4 + 2)
        cart = numpy.hstack((random_cart(j, dtype='z'),
                             random_cart(j, seed=13, dtype='z')))
        cart = numpy.asarray(cart, order='F')
        gsp = numpy.empty((ng, dj), order='F', dtype=numpy.complex128)
        fn(gsp.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(ng), cart.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(kappa), ctypes.c_int(l))
        v1 += fp(gsp)
    assert abs(ref - abs(v1)) < thr

def test_c2s_ket_spinor_sf1(name, ref, thr=1e-12):
    kappa = 0
    v1 = 0
    fn = getattr(_cint, name)
    for j in range(nbas.value*2):
        l = bas[j, ANG_OF]
        dj = (l * 4 + 2)
        cart = random_cart(j)
        gsp = numpy.empty((ng, dj, 2), order='F', dtype=numpy.complex128)
        fn(gsp[:,:,0].ctypes.data_as(ctypes.c_void_p),
           gsp[:,:,1].ctypes.data_as(ctypes.c_void_p),
           cart.ctypes.data_as(ctypes.c_void_p),
           ctypes.c_int(ng), ctypes.c_int(ng), ctypes.c_int(1), ctypes.c_int(kappa), ctypes.c_int(l))
        v1 += fp(gsp)
    assert abs(ref - abs(v1)) < thr

def test_c2s_ket_spinor_si1(name, ref, thr=1e-12):
    kappa = 0
    v1 = 0
    fn = getattr(_cint, name)
    for j in range(nbas.value*2):
        l = bas[j, ANG_OF]
        dj = (l * 4 + 2)
        cart = numpy.hstack([random_cart(j, 12),
                             random_cart(j, 13),
                             random_cart(j, 14),
                             random_cart(j, 15)])
        cart = numpy.asarray(cart, order='F')
        gsp = numpy.empty((ng, dj, 2), order='F', dtype=numpy.complex128)
        fn(gsp[:,:,0].ctypes.data_as(ctypes.c_void_p),
           gsp[:,:,1].ctypes.data_as(ctypes.c_void_p),
           cart.ctypes.data_as(ctypes.c_void_p),
           ctypes.c_int(ng), ctypes.c_int(ng), ctypes.c_int(1), ctypes.c_int(kappa), ctypes.c_int(l))
        v1 += fp(gsp)
    assert abs(ref - abs(v1)) < thr

if __name__ == "__main__":
    test_c2s_ket_sph1('CINTc2s_ket_sph1', 5.597110704570622)
    test_c2s_bra_spinor_e1sf('CINTc2s_bra_spinor_e1sf', 27.53949584857068)
    test_c2s_bra_spinor_sf('CINTc2s_bra_spinor_sf', 39.81581797638495)
    test_c2s_bra_spinor_si('CINTc2s_bra_spinor_si', 12.00244266438091)
    test_c2s_ket_spinor('CINTc2s_ket_spinor', 45.64346086948895)
    test_c2s_ket_spinor('CINTc2s_iket_spinor', 45.64346086948895)

    test_c2s_ket_spinor_sf1('CINTc2s_ket_spinor_sf1', 11.375862317579182)
    test_c2s_ket_spinor_sf1('CINTc2s_iket_spinor_sf1', 11.375862317579182)
    test_c2s_ket_spinor_si1('CINTc2s_ket_spinor_si1', 7.248350912858695)
    test_c2s_ket_spinor_si1('CINTc2s_iket_spinor_si1', 7.248350912858695)
