
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
   1990.8000000              1.0000000
H    S
      5.0250000              0.2709520              0.2
      1.0130000              0.15                   0.5573680
H    S
     80.8000000              0.0210870             -0.0045400              0.0000000
      3.3190000              0.3461290             -0.1703520              0.0000000
      0.9059000              0.0393780              0.1403820              1.0000000
H    P
      4.1330000              0.0868660              0.0000000
      1.2000000              0.0000000              0.5000000
      0.3827000              0.5010080              1.0000000
H    D
      1.0970000              1.0000000
H    D
      2.1330000              0.1868660              0.0000000
      0.3827000              0.2010080              1.0000000
H    F
      0.7610000              1.0000000        
H    F
      1.1330000              0.3868660              1.0000000
      0.8827000              0.4010080              0.0000000
H    g
      1.1330000              0.3868660              0.0000000
      0.8827000              0.4010080              1.0000000
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

def run(intor, comp=1, suffix='_sph', thr=1e-9):
    if suffix == '_spinor':
        intor3 = 'c%s'%intor
    else:
        intor3 = 'c%s%s'%(intor,suffix)
    intor2 = 'c%s%s'%(intor,suffix)
    print(intor)
    fn1 = getattr(_cint, intor3)
    #fn2 = getattr(_cint4, intor2)
    cintopt = make_cintopt(mol._atm, mol._bas, mol._env, intor)
    #cintopt = lib.c_null_ptr()
    args = (mol._atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.natm),
            mol._bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.nbas),
            mol._env.ctypes.data_as(ctypes.c_void_p), cintopt)
    for i in range(mol.nbas):
        for j in range(mol.nbas):
            for k in range(mol.nbas):
                ref = mol.intor_by_shell(intor2, [i,j,k], comp=comp)
                #fn2(ref.ctypes.data_as(ctypes.c_void_p),
                #   (ctypes.c_int*3)(i,j,k),
                #    mol._atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.natm),
                #    mol._bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.nbas),
                #    mol._env.ctypes.data_as(ctypes.c_void_p), lib.c_null_ptr())
                buf = numpy.empty_like(ref)
                fn1(buf.ctypes.data_as(ctypes.c_void_p),
                   (ctypes.c_int*3)(i,j,k),
                    mol._atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.natm),
                    mol._bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.nbas),
                    mol._env.ctypes.data_as(ctypes.c_void_p), lib.c_null_ptr())
                if abs(ref-buf).max() > thr:
                    print(intor, '| nopt', i, j, k, numpy.linalg.norm(ref-buf))#, ref, buf
                    #exit()
                fn1(buf.ctypes.data_as(ctypes.c_void_p),
                    (ctypes.c_int*3)(i,j,k), *args)
                if abs(ref-buf).max() > thr:
                    print(intor, '|', i, j, k, numpy.linalg.norm(ref-buf))
                    #exit()

run('int3c2e')
run('int3c2e_ip1'               , 3)
run('int3c2e_ip2'               , 3)
run('int3c2e_spsp1'             , suffix='_spinor')
run('int3c2e_ipspsp1'           , 3, suffix='_spinor')
run('int3c2e_spsp1ip2'          , 3, suffix='_spinor')
