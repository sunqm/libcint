
import os
import ctypes
import numpy

_cint = numpy.ctypeslib.load_library('libcint', os.path.abspath(os.path.join(__file__, '../../build')))
#_cint4 = ctypes.cdll.LoadLibrary('libcint.so.4')

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

def run(intor, comp=1, suffix='_sph', thr=1e-7):
    if suffix == '_spinor':
        intor3 = 'c%s'%intor
    else:
        intor3 = 'c%s%s'%(intor,suffix)
    intor2 = 'c%s%s'%(intor,suffix)
    print(intor)
    from pyscf.gto.moleintor import libcgto
    if not hasattr(libcgto, intor2):
        print('skip: %s not provided by the installed pyscf' % intor2)
        return
    fn1 = getattr(_cint, intor3)
    #fn2 = getattr(_cint4, intor2)
    cintopt = make_cintopt(mol._atm, mol._bas, mol._env, intor)
    #cintopt = lib.c_null_ptr()
    args = (mol._atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.natm),
            mol._bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.nbas),
            mol._env.ctypes.data_as(ctypes.c_void_p), cintopt)
    failed = False
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
                if numpy.linalg.norm(ref-buf) > thr:
                    print(intor, '| nopt', i, j, k, numpy.linalg.norm(ref-buf))#, ref, buf
                    failed = True
                fn1(buf.ctypes.data_as(ctypes.c_void_p),
                    (ctypes.c_int*3)(i,j,k), *args)
                if numpy.linalg.norm(ref-buf) > thr:
                    print(intor, '|', i, j, k, numpy.linalg.norm(ref-buf))
                    failed = True
    if failed:
        print('failed')
    else:
        print('pass')

def run_fd(intor, intor0, thr=1e-5, h=1e-5):
    '''Check the bra derivative <nabla i|...> against central differences of
    the base integral intor0.  Only shell triplets whose bra atom carries
    neither j nor k are used, so that -d/dR_i intor0 == intor exactly.
    This does not rely on the integral being available (or correct) in the
    installed pyscf.'''
    print(intor, '(finite differences of %s)' % intor0)
    fn1 = getattr(_cint, 'c%s_sph' % intor)
    fn0 = getattr(_cint, 'c%s_sph' % intor0)
    atm = mol._atm
    bas = mol._bas
    env = mol._env.copy()
    c_atm = atm.ctypes.data_as(ctypes.c_void_p)
    c_bas = bas.ctypes.data_as(ctypes.c_void_p)
    c_env = env.ctypes.data_as(ctypes.c_void_p)
    c_natm = ctypes.c_int(mol.natm)
    c_nbas = ctypes.c_int(mol.nbas)
    PTR_COORD = 1
    failed = False
    # keep the cost moderate: bra shells on atom 0, j on atom 1, k on atoms 2,3
    ibas = [n for n in range(mol.nbas) if bas[n,0] == 0]
    jbas = [n for n in range(mol.nbas) if bas[n,0] == 1]
    kbas = [n for n in range(mol.nbas) if bas[n,0] in (2, 3)]
    for i in ibas:
        di = (bas[i,1]*2+1) * bas[i,3]
        for j in jbas:
            dj = (bas[j,1]*2+1) * bas[j,3]
            for k in kbas:
                dk = (bas[k,1]*2+1) * bas[k,3]
                buf = numpy.empty((3,dk,dj,di))
                fn1(buf.ctypes.data_as(ctypes.c_void_p),
                    (ctypes.c_int*3)(i,j,k),
                    c_atm, c_natm, c_bas, c_nbas, c_env, lib.c_null_ptr())
                bp = numpy.empty((dk,dj,di))
                bm = numpy.empty((dk,dj,di))
                fd = numpy.empty_like(buf)
                ptr = atm[0,PTR_COORD]
                for d in range(3):
                    x0 = env[ptr+d]
                    env[ptr+d] = x0 + h
                    fn0(bp.ctypes.data_as(ctypes.c_void_p),
                        (ctypes.c_int*3)(i,j,k),
                        c_atm, c_natm, c_bas, c_nbas, c_env, lib.c_null_ptr())
                    env[ptr+d] = x0 - h
                    fn0(bm.ctypes.data_as(ctypes.c_void_p),
                        (ctypes.c_int*3)(i,j,k),
                        c_atm, c_natm, c_bas, c_nbas, c_env, lib.c_null_ptr())
                    env[ptr+d] = x0
                    fd[d] = -(bp - bm) / (2*h)
                if numpy.abs(buf-fd).max() > thr:
                    print(intor, '| fd', i, j, k, numpy.abs(buf-fd).max())
                    failed = True
    if failed:
        print('failed')
        exit(1)
    else:
        print('pass')

run('int3c1e')
#run('int3c1e_p2')
run('int3c1e_r2_origk')
run('int3c1e_r4_origk')
run('int3c1e_r6_origk', thr=1e-5)
run('int3c1e_ipip1', comp=9)
run('int3c1e_ip1ip2', comp=9)
run('int3c1e_ipvip1', comp=9)
run('int3c1e_ipip2', comp=9)
run('int3c1e_ip1_r2_origk', comp=3)
run('int3c1e_ip1_r4_origk', comp=3, thr=1e-6)
# int3c1e_ip1_r6_origk is deliberately NOT compared against pyscf: released
# libcint <= 6.1.3 computes its y component from an uninitialized buffer, so
# a pyscf built against it is not a valid reference.  The finite-difference
# checks below are self-contained.
run_fd('int3c1e_ip1_r2_origk', 'int3c1e_r2_origk')
run_fd('int3c1e_ip1_r4_origk', 'int3c1e_r4_origk')
run_fd('int3c1e_ip1_r6_origk', 'int3c1e_r6_origk', thr=2e-4)
