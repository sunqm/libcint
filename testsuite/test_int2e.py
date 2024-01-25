#!/usr/bin/env python

import os
import sys
import ctypes
import numpy

_cint = numpy.ctypeslib.load_library('libcint', os.path.abspath(os.path.join(__file__, '../../build')))
#_cint4 = ctypes.cdll.LoadLibrary('libcint.so.4')

from pyscf import gto, lib

mol = gto.M(atom='''H 0 0 0;
            H .2 .5 .8;
            #H 1.9 2.1 .1;
            #H 2.0 .3 1.4''',
            basis = {'H': gto.basis.parse('''
H    S
   1990.8000000              1.0000000
H    S
      5.0250000              0.2709520              0.2
      1.0130000              0.0                    0.5573680
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
H    F
      0.5827000              1.0000000
H    g
      0.8827000              1.0000000
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
        intor = intor = 'c%s'%intor
    else:
        intor = intor = 'c%s%s'%(intor,suffix)
    print(intor)
    failed = False
    fn1 = getattr(_cint, intor)
    #fn2 = getattr(_cint4, intor)
    cintopt = make_cintopt(mol._atm, mol._bas, mol._env, intor)
    args = (mol._atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.natm),
            mol._bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.nbas),
            mol._env.ctypes.data_as(ctypes.c_void_p), cintopt)
    for i in range(mol.nbas):
        for j in range(mol.nbas):
            for k in range(mol.nbas):
                for l in range(mol.nbas):
                    ref = mol.intor_by_shell(intor, [i,j,k,l], comp=comp)
                    #fn2(ref.ctypes.data_as(ctypes.c_void_p),
                    #   (ctypes.c_int*4)(i,j,k,l),
                    #    mol._atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.natm),
                    #    mol._bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.nbas),
                    #    mol._env.ctypes.data_as(ctypes.c_void_p), lib.c_null_ptr())
                    buf = numpy.empty_like(ref)
                    fn1(buf.ctypes.data_as(ctypes.c_void_p),
                       (ctypes.c_int*4)(i,j,k,l),
                        mol._atm.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.natm),
                        mol._bas.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(mol.nbas),
                        mol._env.ctypes.data_as(ctypes.c_void_p), lib.c_null_ptr())
                    if numpy.linalg.norm(ref-buf) > thr:
                        print(intor, '| nopt', i, j, k, l, numpy.linalg.norm(ref-buf))#, ref, buf
                        failed = True
                        #exit()
                    fn1(buf.ctypes.data_as(ctypes.c_void_p),
                        (ctypes.c_int*4)(i,j,k,l), *args)
                    if numpy.linalg.norm(ref-buf) > thr:
                        print(intor, '|', i, j, k, l, numpy.linalg.norm(ref-buf))
                        failed = True
                        #exit()
    if failed:
        print('failed')
    else:
        print('pass')

run('int2e', thr=1e-12)
run("int2e_ip1", 3, thr=1e-12)
if '--quick' in sys.argv:
    exit()

run('int2e', suffix='_cart')
run("int2e_ig1"               , 3)
run("int2e_p1vxp1"            , 3)
#run("int2e_ig1ig2"            , 9)
run("int2e_spsp1"             , suffix='_spinor')
run("int2e_spsp1spsp2"        , suffix='_spinor', thr=1e-6)
run("int2e_srsr1"             , suffix='_spinor')
run("int2e_srsr1srsr2"        , suffix='_spinor')
run("int2e_cg_sa10sp1"        , 3, suffix='_spinor')
run("int2e_cg_sa10sp1spsp2"   , 3, suffix='_spinor')
run("int2e_giao_sa10sp1"      , 3, suffix='_spinor')
run("int2e_giao_sa10sp1spsp2" , 3, suffix='_spinor')
run("int2e_g1"                , 3, suffix='_spinor')
run("int2e_spgsp1"            , 3, suffix='_spinor')
run("int2e_g1spsp2"           , 3, suffix='_spinor')
run("int2e_spgsp1spsp2"       , 3, suffix='_spinor')
#run("int2e_pp1"               , suffix='_spinor')
#run("int2e_pp2"               , suffix='_spinor')
#run("int2e_pp1pp2"            , suffix='_spinor')
run("int2e_spv1"              , suffix='_spinor')
run("int2e_vsp1"              , suffix='_spinor')
run("int2e_spsp2"             , suffix='_spinor')
run("int2e_spv1spv2"          , suffix='_spinor', thr=1e-6)
run("int2e_vsp1spv2"          , suffix='_spinor', thr=1e-6)
run("int2e_spv1vsp2"          , suffix='_spinor', thr=1e-6)
run("int2e_vsp1vsp2"          , suffix='_spinor', thr=1e-6)
run("int2e_spv1spsp2"         , suffix='_spinor', thr=1e-6)
run("int2e_vsp1spsp2"         , suffix='_spinor', thr=1e-6)
run("int2e_ip1"               , 3, suffix='_spinor')
run("int2e_ipspsp1"           , 3, suffix='_spinor')
run("int2e_ip1spsp2"          , 3, suffix='_spinor')
run("int2e_ipspsp1spsp2"      , 3, suffix='_spinor', thr=1e-5)
run("int2e_ipsrsr1"           , 3, suffix='_spinor')
run("int2e_ip1srsr2"          , 3, suffix='_spinor')
run("int2e_ipsrsr1srsr2"      , 3, suffix='_spinor')
run("int2e_ssp1ssp2"          , suffix='_spinor')
run("int2e_ssp1sps2"          , suffix='_spinor')
run("int2e_sps1ssp2"          , suffix='_spinor')
run("int2e_sps1sps2"          , suffix='_spinor')
run("int2e_cg_ssa10ssp2"      , 3, suffix='_spinor')
run("int2e_giao_ssa10ssp2"    , 3, suffix='_spinor')
run("int2e_gssp1ssp2"         , 3, suffix='_spinor')
run("int2e_gauge_r1_ssp1ssp2" , suffix='_spinor', thr=1e-6)
run("int2e_gauge_r1_ssp1sps2" , suffix='_spinor', thr=1e-6)
run("int2e_gauge_r1_sps1ssp2" , suffix='_spinor', thr=1e-6)
run("int2e_gauge_r1_sps1sps2" , suffix='_spinor', thr=1e-6)
run("int2e_gauge_r2_ssp1ssp2" , suffix='_spinor', thr=1e-6)
run("int2e_gauge_r2_ssp1sps2" , suffix='_spinor', thr=1e-6)
run("int2e_gauge_r2_sps1ssp2" , suffix='_spinor', thr=1e-6)
run("int2e_gauge_r2_sps1sps2" , suffix='_spinor', thr=1e-6)
run("int2e_ipip1"             , 9)
run("int2e_ipvip1"            , 9)
run("int2e_ip1ip2"            , 9)
