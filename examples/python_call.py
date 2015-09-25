from math import *
import numpy
import ctypes

# general contracted DZ basis [3s1p/2s1p] for H2
#     exponents    contract-coeff
# S   6.0          0.7               0.4
#     2.0          0.6               0.3
#     0.8          0.5               0.2
# P   0.9          1.

def gto_norm(n, a):
    # normalization factor of function r^n e^{-a r^2}
    s = 2**(2*n+3) * factorial(n+1) * (2*a)**(n+1.5) \
            / (factorial(2*n+2) * sqrt(pi))
    return sqrt(s)

CHARGE_OF  = 1
PTR_COORD  = 2
NUC_MOD_OF = 3
PTR_ZETA   = 4
ATM_SLOTS  = 6

ATOM_OF    = 1
ANG_OF     = 2
NPRIM_OF   = 3
NCTR_OF    = 4
KAPPA_OF   = 5
PTR_EXP    = 6
PTR_COEFF  = 7
BAS_SLOTS  = 8

ptr_env = 20
atm = []
bas = []
env = [0] * ptr_env

i = 0
#           CHARGE_OF,PTR_COORD
atm.append([1,        ptr_env,  0, 0, 0, 0])
#           x  y  z (Bohr)
env.extend([0, 0, -0.8])
ptr_env += 3
atm.append([1,        ptr_env,  0, 0, 0, 0])
env.extend([0, 0,  0.8])

# basis for atom #0
# 3s -> 2s
env.extend([6., 2., .8])
env.extend([.7, .6, .5, .4, .3, .2])
#           ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, KAPPA_OF, PTR_EXP, PTR_COEFF
bas.append([0,       0,      3,        2,       0,        ptr_env, ptr_env+3, 0])
ptr_env += 9
env.extend([.9])
env.extend([1.])
bas.append([0,       1,      1,        1,       0,        ptr_env, ptr_env+1, 0])
ptr_env += 1

# basis functions for atom #1, they are the same to thoes of atom #0
bas.extend(bas[-2:])

# note the integer type
atm = numpy.array(atm,dtype=int32)
bas = numpy.array(bas,dtype=int32)
env = numpy.array(env)

_cint = ctypes.cdll.LoadLibrary('/path/to/libcint.so')
c_natm = ctypes.c_int(atm.size)
c_nbas = ctypes.c_int(bas.size)

_cint.CINTcgto_spheric.restype = ctypes.c_int
di = _cint.CINTcgto_spheric(ctypes.c_int(0), bas.ctype.data)
dj = _cint.CINTcgto_spheric(ctypes.c_int(1), bas.ctype.data)

c_shls = (ctypes.c_int * 2)(0, 1)
buf = numpy.empty((di.value, dj.value, 3))
_cint.cint1e_ipnuc_sph(buf.ctypes.data, c_shls, \
                       atm.ctypes.data, c_natm, \
                       bas.ctypes.data, c_nbas,
                       env.ctypes.data)
