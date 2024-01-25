
import os
import ctypes
import numpy

_cint = ctypes.CDLL(os.path.abspath(os.path.join(__file__, '../../build/libcint.so')))

from pyscf import gto, lib

mol = gto.M(atom='''H 0 0 0;
            H .2 .5 .8;
            H 1.9 2.1 .1;
            H 2.0 .3 1.4''',
            basis = {'H': gto.basis.parse('''
H    S
      2.3190000              1.
H    P
      1.2000000              1.
H    D
      1.0970000              1.
H    F
      0.7610000              1.
      ''')})

def make_cintopt(atm, bas, env, intor):
    c_atm = numpy.asarray(atm, dtype=numpy.int32, order='C')
    c_bas = numpy.asarray(bas, dtype=numpy.int32, order='C')
    c_env = numpy.asarray(env, dtype=numpy.double, order='C')
    natm = c_atm.shape[0]
    nbas = c_bas.shape[0]
    cintopt = lib.c_null_ptr()
    foptinit = getattr(_cint, intor+'_optimizer')
    foptinit(ctypes.byref(cintopt),
             c_atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(natm),
             c_bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(nbas),
             c_env.ctypes.data_as(ctypes.c_void_p))
    return cintopt

fn1 = getattr(_cint, 'cint2e_sph')
fn2 = getattr(_cint, 'cint2e_yp_sph')
cintopt1 = lib.c_null_ptr()
cintopt1 = make_cintopt(mol._atm, mol._bas, mol._env, 'cint2e_yp_sph')

ao_loc = mol.ao_loc_nr()
mol._env[8] = 1e3
mol._env[9] = 1e-3
mol._env[10] = 1e-4

yp_err = []
for i in range(mol.nbas):
    for j in range(mol.nbas):
        for k in range(mol.nbas):
            for l in range(mol.nbas):
                di = ao_loc[i+1] - ao_loc[i]
                dj = ao_loc[j+1] - ao_loc[j]
                dk = ao_loc[k+1] - ao_loc[k]
                dl = ao_loc[l+1] - ao_loc[l]
                ref = numpy.empty([di,dj,dk,dl])
                buf = numpy.empty([di,dj,dk,dl])
                fn1(ref.ctypes.data_as(ctypes.c_void_p),
                   (ctypes.c_int*4)(i,j,k,l),
                    mol._atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.natm),
                    mol._bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.nbas),
                    mol._env.ctypes.data_as(ctypes.c_void_p), lib.c_null_ptr())
                fn2(buf.ctypes.data_as(ctypes.c_void_p),
                   (ctypes.c_int*4)(i,j,k,l),
                    mol._atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.natm),
                    mol._bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.nbas),
                    mol._env.ctypes.data_as(ctypes.c_void_p), cintopt1)
                err = numpy.linalg.norm(ref-buf)
                yp_err.append(err)
                if err > 1e-3:
                    print('yp', i, j, k, l, numpy.linalg.norm(ref-buf)/numpy.linalg.norm(ref))
                    #exit()
print(max(yp_err))
