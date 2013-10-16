from math import *
import numpy
import ctypes

def charge_of(atm_id):
    return 1

def coord_of(atm_id):
    return (0,1,2)

def angular_of(bas_id):
    return 1

def nprimitive_of(bas_id):
    return 3

def ncontract_of(bas_id):
    return 2

def kappa_of(bas_id):
    return 0

def coeff_of(bas_id):
    # 3 primitive-GTO -> 2 contracted-GTO
    c = numpy.random.random((3,2)).flatten(order='F')
    n = angular_of(bas_id)
    # absorb the normalization factor of primitive GTO into coefficients
    for i, a in enumerate(exponent_of(bas_id)):
        c[i,:] *= gto_norm(a, n)
    return c

def exponent_of(bas_id):
    # 3 primitive-GTO -> 2 contracted-GTO
    return numpy.random.random(3)

def gto_norm(n, a):
    # normalization factor of function r^n e^{-a r^2}
    s = 2**(2*n+3) * factorial(n+1) * (2*a)**(n+1.5) \
            / (factorial(2*n+2) * sqrt(pi))
    return sqrt(s)

ptr_env = 10
atm = []
bas = []
env = [0] * ptr_env
natm = 2
for atm_id in range(natm):
    env.extend(coord_of(atm_id))
    atm.append([charge_of(atm_id), ptr_env, nuc_mod_of(atm_id), 0, 0, 0])
    ptr_env += 3
    for bas_id in basis_of(atm_id):
        c = coeff_of(bas_id)
        e = exponent_of(bas_id)
        env.extend(e)
        env.extend(c)

        bas.append([atm_id, angular_of(bas_id), \
                    nprimitive_of(bas_id), \
                    ncontract_of(bas_id), \
                    kappa_of(bas_id), \
                    ptr_env, ptr_env+e.size, 0])
        ptr_env += e.size + c.size
atm = numpy.array(atm,dtype=int32)
bas = numpy.array(bas,dtype=int32)
env = numpy.array(env)

_cint = ctypes.cdll.LoadLibrary('/path/to/libcint.so')
c_natm = ctypes.c_int(atm.size)
c_nbas = ctypes.c_int(bas.size)

_cint.cgtos_spheric.restype = ctypes.c_int
di = _cint.cgtos_spheric(ctypes.c_int(0), bas.ctype.data)
dj = _cint.cgtos_spheric(ctypes.c_int(1), bas.ctype.data)

c_shls = (ctypes.c_int * 2)(0, 1)
buf = numpy.empty((di.value, dj.value))
_cint.cint1e_ovlp_sph(buf.ctypes.data, c_shls, \
                      atm.ctypes.data, c_natm, \
                      bas.ctypes.data, c_nbas,
                      env.ctypes.data)
