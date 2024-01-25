#!/usr/bin/env python
# $Id$
# -*- coding: utf-8

import sys
import os
import ctypes
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(__file__, '../../scripts')))
import cart2sph

_cint = ctypes.CDLL(os.path.abspath(os.path.join(__file__, '../../build/libcint.so')))

pauli = np.array([[[0., 1.],
                   [1., 0.]],  # x
                  [[0.,-1j],
                   [1j, 0.]],  # y
                  [[1., 0.],
                   [0.,-1.]],  # z
                 ])
si = np.vstack([1j * pauli, np.eye(2)[None,:,:]])

np.random.seed(2)
ldc = 2
lds = 3
nctr = 2

class CINTEnvVars(ctypes.Structure):
    _fields_ = [
        ('atm', ctypes.c_void_p),
        ('bas', ctypes.c_void_p),
        ('env', ctypes.c_void_p),
        ('shls', ctypes.c_void_p),
        ('natm', ctypes.c_int),
        ('nbas', ctypes.c_int),
        ('i_l', ctypes.c_int),
        ('j_l', ctypes.c_int),
        ('k_l', ctypes.c_int),
        ('l_l', ctypes.c_int),
        ('nfi', ctypes.c_int),
        ('nfj', ctypes.c_int),
        ('nfk', ctypes.c_int),
        ('nfl', ctypes.c_int),
        ('nf', ctypes.c_int),
        ('_padding', ctypes.c_int),
        ('x_ctr', ctypes.c_int * 4)]


def make_gxyz1(l):
    ncart = (l + 1) * (l + 2) // 2
    nsph = l * 2 + 1
    gx = np.random.rand(nctr, ncart, ldc)
    gy = np.random.rand(nctr, ncart, ldc)
    gz = np.random.rand(nctr, ncart, ldc)
    g1 = np.random.rand(nctr, ncart, ldc)
    gcart = np.array([gx, gy, gz, g1])
    gspa = np.empty((nctr, nsph*2, lds), dtype=np.complex128)
    gspb = np.empty((nctr, nsph*2, lds), dtype=np.complex128)
    gspa[:] = np.nan
    gspb[:] = np.nan
    return gspa, gspb, gcart

def c_c2s(l):
    ua, ub = cart2sph.cart2spinor(l)
    return np.array([ua.astype(np.complex128), ub.astype(np.complex128)])

def check_ket_spinor_si1(l):
    gspa, gspb, gcart = make_gxyz1(l)
    _cint.CINTc2s_ket_spinor_si1(
        gspa.ctypes.data_as(ctypes.c_void_p),
        gspb.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(lds), ctypes.c_int(ldc), ctypes.c_int(nctr),
        ctypes.c_int(0), ctypes.c_int(l))

    c2s = c_c2s(l)
    ref = np.einsum('qxy,qncg,ycj->xnjg', si, gcart, c2s)
    return (abs(gspa[:,:,:ldc] - ref[0]).max() < 1e-12 and
            abs(gspb[:,:,:ldc] - ref[1]).max() < 1e-12)

def check_ket_spinor_sf1(l):
    gspa, gspb, gcart = make_gxyz1(l)
    _cint.CINTc2s_ket_spinor_sf1(
        gspa.ctypes.data_as(ctypes.c_void_p),
        gspb.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(lds), ctypes.c_int(ldc), ctypes.c_int(nctr),
        ctypes.c_int(0), ctypes.c_int(l))

    c2s = c_c2s(l)
    ref = np.einsum('xy,ncg,ycj->xnjg', np.eye(2), gcart[0], c2s)
    return (abs(gspa[:,:,:ldc] - ref[0]).max() < 1e-12 and
            abs(gspb[:,:,:ldc] - ref[1]).max() < 1e-12)

def check_iket_spinor_sf1(l):
    gspa, gspb, gcart = make_gxyz1(l)
    _cint.CINTc2s_iket_spinor_sf1(
        gspa.ctypes.data_as(ctypes.c_void_p),
        gspb.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(lds), ctypes.c_int(ldc), ctypes.c_int(nctr),
        ctypes.c_int(0), ctypes.c_int(l))

    c2s = c_c2s(l)
    ref = np.einsum('xy,ncg,ycj->xnjg', np.eye(2)*1j, gcart[0], c2s)
    return (abs(gspa[:,:,:ldc] - ref[0]).max() < 1e-12 and
            abs(gspb[:,:,:ldc] - ref[1]).max() < 1e-12)

def check_iket_spinor_si1(l):
    gspa, gspb, gcart = make_gxyz1(l)
    _cint.CINTc2s_iket_spinor_si1(
        gspa.ctypes.data_as(ctypes.c_void_p),
        gspb.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(lds), ctypes.c_int(ldc), ctypes.c_int(nctr),
        ctypes.c_int(0), ctypes.c_int(l))

    c2s = c_c2s(l)
    ref = np.einsum('qxy,qncg,ycj->xnjg', si*1j, gcart, c2s)
    return (abs(gspa[:,:,:ldc] - ref[0]).max() < 1e-11 and
            abs(gspb[:,:,:ldc] - ref[1]).max() < 1e-11)

def check_ket_spinor(l):
    ncart = (l + 1) * (l + 2) // 2
    nsph = l * 2 + 1
    nbra = ldc
    gcart = (np.random.rand(ncart*2, nbra) + np.random.rand(ncart*2, nbra) * 1j)
    gsp = np.empty((nsph*2, nbra), dtype=np.complex128)
    _cint.CINTc2s_ket_spinor(
        gsp.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(nbra),
        gcart.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(0), ctypes.c_int(l))

    c2s = c_c2s(l)
    c2s = c2s.reshape(ncart*2, nsph*2)
    ref = np.einsum('cg,cj->jg', gcart, c2s)
    return abs(gsp - ref).max() < 1e-12

def check_iket_spinor(l):
    ncart = (l + 1) * (l + 2) // 2
    nsph = l * 2 + 1
    nbra = ldc
    gcart = (np.random.rand(ncart*2, nbra) + np.random.rand(ncart*2, nbra) * 1j)
    gsp = np.empty((nsph*2, nbra), dtype=np.complex128)
    _cint.CINTc2s_iket_spinor(
        gsp.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(nbra),
        gcart.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(0), ctypes.c_int(l))

    c2s = c_c2s(l)
    c2s = c2s.reshape(ncart*2, nsph*2)
    ref = np.einsum('cg,cj->jg', gcart*1j, c2s)
    return abs(gsp - ref).max() < 1e-12

def check_bra_spinor(l):
    ncart = (l + 1) * (l + 2) // 2
    nsph = l * 2 + 1
    nket = ldc
    gcart = (np.random.rand(ncart*2, nket) + np.random.rand(ncart*2, nket) * 1j)
    gcart = gcart.T.copy('C')
    gsp = np.empty((nket, nsph*2), dtype=np.complex128)
    _cint.CINTc2s_bra_spinor(
        gsp.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(nket),
        gcart.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(0), ctypes.c_int(l))

    c2s = c_c2s(l)
    c2s = c2s.reshape(ncart*2, nsph*2)
    ref = np.einsum('gc,cj->gj', gcart, c2s.conj())
    return abs(gsp - ref).max() < 1e-12

def check_bra_spinor_si(l):
    ncart = (l + 1) * (l + 2) // 2
    nsph = l * 2 + 1
    nket = ldc
    gcart = (np.random.rand(2,nket,ncart) + np.random.rand(2,nket,ncart) * 1j)
    gsp = np.empty((nket, nsph*2), dtype=np.complex128)
    _cint.CINTc2s_bra_spinor_si(
        gsp.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(nket),
        gcart.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(0), ctypes.c_int(l))

    c2s = c_c2s(l)
    ref = np.einsum('xgc,xcj->gj', gcart, c2s.conj())
    return abs(gsp - ref).max() < 1e-12

def check_bra_spinor_sf(l):
    ncart = (l + 1) * (l + 2) // 2
    nsph = l * 2 + 1
    nket = ldc
    gcart = (np.random.rand(ncart, nket) + np.random.rand(ncart, nket) * 1j)
    gcart = gcart.T.copy('C')
    gsp = np.empty((nket*2, nsph*2), dtype=np.complex128)
    _cint.CINTc2s_bra_spinor_sf(
        gsp.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(nket),
        gcart.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(0), ctypes.c_int(l))

    c2s = c_c2s(l)
    ref = np.einsum('xy,gc,xcj->ygj', np.eye(2), gcart, c2s.conj())
    return abs(gsp - ref.reshape(gsp.shape)).max() < 1e-12

def check_bra_spinor_e1sf(l):
    ncart = (l + 1) * (l + 2) // 2
    nsph = l * 2 + 1
    nket = ldc
    gcart = np.random.rand(ncart, nket)
    gcart = gcart.T.copy('C')
    gsp = np.empty((nket*2, nsph*2), dtype=np.complex128)
    _cint.CINTc2s_bra_spinor_e1sf(
        gsp.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(nket),
        gcart.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(0), ctypes.c_int(l))

    c2s = c_c2s(l)
    ref = np.einsum('xy,gc,xcj->ygj', np.eye(2), gcart, c2s.conj())
    return abs(gsp - ref.reshape(gsp.shape)).max() < 1e-12

def make_envs_1e(li, lj, sf=True):
    nfi = (li + 1) * (li + 2) // 2
    nfj = (lj + 1) * (lj + 2) // 2
    nf = nfi * nfj
    di = li * 2 + 1
    dj = lj * 2 + 1
    nfk = nfl = dk = dl = 1
    dims = [di*2*nctr, dj*2*nctr]
    if sf:
        shape = [nfi, nfj, nctr, nctr]
    else:
        shape = [nfi, nfj, nctr, nctr, 4]
    atm = bas = env = np.zeros((2, 8), dtype=np.int32)
    shls = np.array([0, 1], dtype=np.int32)
    x_ctr = (ctypes.c_int * 4)(nctr, nctr, 1, 1)
    dims = np.asarray(dims, dtype=np.int32)
    cache = np.empty(nf*16)
    gcart = np.asarray(np.random.random(shape), order='F')

    envs = CINTEnvVars(atm.ctypes.data_as(ctypes.c_void_p),
                       bas.ctypes.data_as(ctypes.c_void_p),
                       env.ctypes.data_as(ctypes.c_void_p),
                       shls.ctypes.data_as(ctypes.c_void_p),
                       4, 4, li, lj, 0, 0, nfi, nfj, 1, 1, nf,
                       0, x_ctr)
    return envs, gcart, dims, cache

def make_envs_1e_grids(li, lj, sf=True):
    nfi = (li + 1) * (li + 2) // 2
    nfj = (lj + 1) * (lj + 2) // 2
    ngrids = 136
    GRID_BLOCK = 104
    nf = nfi * nfj
    di = li * 2 + 1
    dj = lj * 2 + 1
    nfk = nfl = dk = dl = 1
    dims = [di*2*nctr, dj*2*nctr, ngrids]
    if sf:
        shape = [ngrids, nfi, nfj, nctr, nctr]
    else:
        shape = [ngrids, nfi, nfj, nctr, nctr, 4]
    atm = bas = env = np.zeros((2, 8), dtype=np.int32)
    shls = np.array([0, 1], dtype=np.int32)
    x_ctr = (ctypes.c_int * 4)(nctr, nctr, 1, 1)
    dims = np.asarray(dims, dtype=np.int32)
    cache = np.empty(GRID_BLOCK*nf*16+10)
    gcart0 = np.asarray(np.random.random(shape), order='F')
    gcart1 = []
    if sf:
        for i0 in range(0, ngrids, GRID_BLOCK):
            i1 = min(i0 + GRID_BLOCK, ngrids)
            gcart1.append(gcart0[i0:i1].ravel('F'))
    else:
        for ic in range(4):
            for i0 in range(0, ngrids, GRID_BLOCK):
                i1 = min(i0 + GRID_BLOCK, ngrids)
                gcart1.append(gcart0[i0:i1,:,:,:,:,ic].ravel('F'))
    gcart1 = np.hstack(gcart1)

    envs = CINTEnvVars(atm.ctypes.data_as(ctypes.c_void_p),
                       bas.ctypes.data_as(ctypes.c_void_p),
                       env.ctypes.data_as(ctypes.c_void_p),
                       shls.ctypes.data_as(ctypes.c_void_p),
                       4, 4, li, lj, 0, 0, nfi, nfj, 1, ngrids, nf,
                       0, x_ctr)
    return envs, gcart0, gcart1, dims, cache

def make_envs_2e1(li, lj, lk, ll, sf=True):
    nfi = (li + 1) * (li + 2) // 2
    nfj = (lj + 1) * (lj + 2) // 2
    nfk = (lk + 1) * (lk + 2) // 2
    nfl = (ll + 1) * (ll + 2) // 2
    nf = nfi * nfj * nfk * nfl
    di = li * 2 + 1
    dj = lj * 2 + 1
    dk = lk * 2 + 1
    dl = ll * 2 + 1
    dims = [di*2, nfk, nfl, dj*2, 2, nctr, nctr]
    if sf:
        shape = [nfi, nfk, nfl, nfj, nctr, nctr]
    else:
        shape = [nfi, nfk, nfl, nfj, nctr, nctr, 4]
    atm = bas = env = np.zeros((4, 8), dtype=np.int32)
    shls = np.array([0, 1, 2, 3], dtype=np.int32)
    x_ctr = (ctypes.c_int * 4)(nctr, nctr, 1, 1)
    dims = np.asarray(dims, dtype=np.int32)
    cache = np.empty(nf*64)
    gcart = np.asarray(np.random.random(shape), order='F')

    envs = CINTEnvVars(atm.ctypes.data_as(ctypes.c_void_p),
                       bas.ctypes.data_as(ctypes.c_void_p),
                       env.ctypes.data_as(ctypes.c_void_p),
                       shls.ctypes.data_as(ctypes.c_void_p),
                       4, 4, li, lj, lk, ll, nfi, nfj, nfk, nfl, nf,
                       0, x_ctr)
    return envs, gcart, dims, cache

def make_envs_2e2(li, lj, lk, ll, sf=True):
    nfi = (li + 1) * (li + 2) // 2
    nfj = (lj + 1) * (lj + 2) // 2
    nfk = (lk + 1) * (lk + 2) // 2
    nfl = (ll + 1) * (ll + 2) // 2
    nf = nfi * nfj * nfk * nfl
    di = li * 2 + 1
    dj = lj * 2 + 1
    dk = lk * 2 + 1
    dl = ll * 2 + 1
    dims = [di*2*nctr, dj*2*nctr, dk*2, dl*2]
    if sf:
        shape = [di*2, nfk, nfl, dj*2, 2, nctr, nctr]
    else:
        shape = [di*2, nfk, nfl, dj*2, 2, nctr, nctr, 4]
    atm = bas = env = np.zeros((4, 8), dtype=np.int32)
    shls = np.array([0, 1, 2, 3], dtype=np.int32)
    x_ctr = (ctypes.c_int * 4)(nctr, nctr, 1, 1)
    dims = np.asarray(dims, dtype=np.int32)
    cache = np.empty(nf*64)
    gcart = np.asarray(np.random.random(shape), order='F')

    envs = CINTEnvVars(atm.ctypes.data_as(ctypes.c_void_p),
                       bas.ctypes.data_as(ctypes.c_void_p),
                       env.ctypes.data_as(ctypes.c_void_p),
                       shls.ctypes.data_as(ctypes.c_void_p),
                       4, 4, li, lj, lk, ll, nfi, nfj, nfk, nfl, nf,
                       0, x_ctr)
    return envs, gcart, dims, cache

def make_envs_3c2e1(li, lj, lk, sf=True, ssc=False):
    nfi = (li + 1) * (li + 2) // 2
    nfj = (lj + 1) * (lj + 2) // 2
    nfk = (lk + 1) * (lk + 2) // 2
    nf = nfi * nfj * nfk
    di = li * 2 + 1
    dj = lj * 2 + 1
    dk = lk * 2 + 1
    if ssc:
        dims = [di*2*nctr, dj*2*nctr, nfk]
    else:
        dims = [di*2*nctr, dj*2*nctr, dk]
    if sf:
        shape = [nfi, nfk, nfj, nctr, nctr]
    else:
        shape = [nfi, nfk, nfj, nctr, nctr, 4]
    atm = bas = env = np.zeros((4, 8), dtype=np.int32)
    shls = np.array([0, 1, 2, 3], dtype=np.int32)
    x_ctr = (ctypes.c_int * 4)(nctr, nctr, 1, 1)
    dims = np.asarray(dims, dtype=np.int32)
    cache = np.empty(nf*64)
    gcart = np.asarray(np.random.random(shape), order='F')

    envs = CINTEnvVars(atm.ctypes.data_as(ctypes.c_void_p),
                       bas.ctypes.data_as(ctypes.c_void_p),
                       env.ctypes.data_as(ctypes.c_void_p),
                       shls.ctypes.data_as(ctypes.c_void_p),
                       4, 4, li, lj, lk, 0, nfi, nfj, nfk, 1, nf,
                       0, x_ctr)
    return envs, gcart, dims, cache

def test_c2s_sf_1e(li, lj):
    envs, gcart, dims, cache = make_envs_1e(li, lj)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_sf_1e(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,pqmn,xpi,yqj->minj', np.eye(2), gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_sf_1ei(li, lj):
    envs, gcart, dims, cache = make_envs_1e(li, lj)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_sf_1ei(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,pqmn,xpi,yqj->minj', np.eye(2)*1j, gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_si_1e(li, lj):
    envs, gcart, dims, cache = make_envs_1e(li, lj, sf=False)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_si_1e(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,pqmns,xpi,yqj->minj', si, gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_si_1ei(li, lj):
    envs, gcart, dims, cache = make_envs_1e(li, lj, sf=False)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_si_1ei(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,pqmns,xpi,yqj->minj', si*1j, gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_sf_2e1(li, lj, lk, ll):
    envs, gcart, dims, cache = make_envs_2e1(li, lj, lk, ll)
    out = np.empty(dims, order='F')
    _cint.c2s_sf_2e1(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,pklqmn,xpi,yqj->ikljmn', np.eye(2), gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return (abs(out[:,:,:,:,0] - ref.real).max() < 1e-12 and
            abs(out[:,:,:,:,1] - ref.imag).max() < 1e-12)

def test_c2s_sf_2e1i(li, lj, lk, ll):
    envs, gcart, dims, cache = make_envs_2e1(li, lj, lk, ll)
    out = np.empty(dims, order='F')
    _cint.c2s_sf_2e1i(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,pklqmn,xpi,yqj->ikljmn', np.eye(2)*1j, gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return (abs(out[:,:,:,:,0] - ref.real).max() < 1e-12 and
            abs(out[:,:,:,:,1] - ref.imag).max() < 1e-12)

def test_c2s_sf_2e2(li, lj, lk, ll):
    envs, gcart, dims, cache = make_envs_2e2(li, lj, lk, ll)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_sf_2e2(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,ipqjmn,xpk,yql->minjkl', np.eye(2),
                    gcart[:,:,:,:,0]+gcart[:,:,:,:,1]*1j,
                    c_c2s(lk).conj(), c_c2s(ll))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_sf_2e2i(li, lj, lk, ll):
    envs, gcart, dims, cache = make_envs_2e2(li, lj, lk, ll)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_sf_2e2i(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,ipqjmn,xpk,yql->minjkl', np.eye(2)*1j,
                    gcart[:,:,:,:,0]+gcart[:,:,:,:,1]*1j,
                    c_c2s(lk).conj(), c_c2s(ll))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_si_2e1(li, lj, lk, ll):
    envs, gcart, dims, cache = make_envs_2e1(li, lj, lk, ll, sf=False)
    out = np.empty(dims, order='F')
    _cint.c2s_si_2e1(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,pklqmns,xpi,yqj->ikljmn', si, gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return (abs(out[:,:,:,:,0] - ref.real).max() < 1e-12 and
            abs(out[:,:,:,:,1] - ref.imag).max() < 1e-12)

def test_c2s_si_2e1i(li, lj, lk, ll):
    envs, gcart, dims, cache = make_envs_2e1(li, lj, lk, ll, sf=False)
    out = np.empty(dims, order='F')
    _cint.c2s_si_2e1i(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,pklqmns,xpi,yqj->ikljmn', si*1j, gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return (abs(out[:,:,:,:,0] - ref.real).max() < 1e-12 and
            abs(out[:,:,:,:,1] - ref.imag).max() < 1e-12)

def test_c2s_si_2e2(li, lj, lk, ll):
    envs, gcart, dims, cache = make_envs_2e2(li, lj, lk, ll, sf=False)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_si_2e2(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,ipqjmns,xpk,yql->minjkl', si,
                    gcart[:,:,:,:,0]+gcart[:,:,:,:,1]*1j,
                    c_c2s(lk).conj(), c_c2s(ll))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_si_2e2i(li, lj, lk, ll):
    envs, gcart, dims, cache = make_envs_2e2(li, lj, lk, ll, sf=False)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_si_2e2i(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,ipqjmns,xpk,yql->minjkl', si*1j,
                    gcart[:,:,:,:,0]+gcart[:,:,:,:,1]*1j,
                    c_c2s(lk).conj(), c_c2s(ll))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_sf_3c2e1(li, lj, lk):
    envs, gcart, dims, cache = make_envs_3c2e1(li, lj, lk)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_sf_3c2e1(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,prqmn,xpi,yqj,rk->minjk', np.eye(2), gcart,
                    c_c2s(li).conj(), c_c2s(lj),
                    cart2sph.cart2sph(lk).astype(float))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_sf_3c2e1i(li, lj, lk):
    envs, gcart, dims, cache = make_envs_3c2e1(li, lj, lk)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_sf_3c2e1i(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,prqmn,xpi,yqj,rk->minjk', np.eye(2)*1j, gcart,
                    c_c2s(li).conj(), c_c2s(lj),
                    cart2sph.cart2sph(lk).astype(float))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_si_3c2e1(li, lj, lk):
    envs, gcart, dims, cache = make_envs_3c2e1(li, lj, lk, sf=False)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_si_3c2e1(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,prqmns,xpi,yqj,rk->minjk', si, gcart,
                    c_c2s(li).conj(), c_c2s(lj),
                    cart2sph.cart2sph(lk).astype(float))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_si_3c2e1i(li, lj, lk):
    envs, gcart, dims, cache = make_envs_3c2e1(li, lj, lk, sf=False)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_si_3c2e1i(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,prqmns,xpi,yqj,rk->minjk', si*1j, gcart,
                    c_c2s(li).conj(), c_c2s(lj),
                    cart2sph.cart2sph(lk).astype(float))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_sf_3c2e1_ssc(li, lj, lk):
    envs, gcart, dims, cache = make_envs_3c2e1(li, lj, lk, ssc=True)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_sf_3c2e1_ssc(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,pkqmn,xpi,yqj->minjk', np.eye(2), gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_sf_3c2e1i_ssc(li, lj, lk):
    envs, gcart, dims, cache = make_envs_3c2e1(li, lj, lk, ssc=True)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_sf_3c2e1i_ssc(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,pkqmn,xpi,yqj->minjk', np.eye(2)*1j, gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_si_3c2e1_ssc(li, lj, lk):
    envs, gcart, dims, cache = make_envs_3c2e1(li, lj, lk, sf=False, ssc=True)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_si_3c2e1_ssc(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,pkqmns,xpi,yqj->minjk', si, gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_si_3c2e1i_ssc(li, lj, lk):
    envs, gcart, dims, cache = make_envs_3c2e1(li, lj, lk, sf=False, ssc=True)
    out = np.empty(dims, order='F', dtype=np.complex128)
    _cint.c2s_si_3c2e1i_ssc(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,pkqmns,xpi,yqj->minjk', si*1j, gcart,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_sf_1e_grids(li, lj):
    envs, gcart0, gcart1, dims, cache = make_envs_1e_grids(li, lj)
    out = np.empty(dims[[2,0,1]], order='F', dtype=np.complex128)
    _cint.c2s_sf_1e_grids(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart1.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,gpqmn,xpi,yqj->gminj', np.eye(2), gcart0,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_sf_1e_gridsi(li, lj):
    envs, gcart0, gcart1, dims, cache = make_envs_1e_grids(li, lj)
    out = np.empty(dims[[2,0,1]], order='F', dtype=np.complex128)
    _cint.c2s_sf_1e_gridsi(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart1.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('xy,gpqmn,xpi,yqj->gminj', np.eye(2)*1j, gcart0,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_si_1e_grids(li, lj):
    envs, gcart0, gcart1, dims, cache = make_envs_1e_grids(li, lj, sf=False)
    out = np.empty(dims[[2,0,1]], order='F', dtype=np.complex128)
    _cint.c2s_si_1e_grids(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart1.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,gpqmns,xpi,yqj->gminj', si, gcart0,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

def test_c2s_si_1e_gridsi(li, lj):
    envs, gcart0, gcart1, dims, cache = make_envs_1e_grids(li, lj, sf=False)
    out = np.empty(dims[[2,0,1]], order='F', dtype=np.complex128)
    _cint.c2s_si_1e_gridsi(
        out.ctypes.data_as(ctypes.c_void_p),
        gcart1.ctypes.data_as(ctypes.c_void_p),
        dims.ctypes.data_as(ctypes.c_void_p),
        ctypes.byref(envs),
        cache.ctypes.data_as(ctypes.c_void_p))

    ref = np.einsum('sxy,gpqmns,xpi,yqj->gminj', si*1j, gcart0,
                    c_c2s(li).conj(), c_c2s(lj))
    return abs(out - ref.reshape(out.shape)).max() < 1e-12

if __name__ == "__main__":
    for fn in (
        test_c2s_sf_1e,
        test_c2s_sf_1ei,
        test_c2s_si_1e,
        test_c2s_si_1ei,
    ):
        passed = True
        for li in range(4):
            for lj in range(4):
                if not fn(li, lj):
                    passed = False
        if passed:
            print(f'{fn.__name__}       passed')
        else:
            print(f'{fn.__name__}       failed')

    for fn in (
        test_c2s_sf_2e1,
        test_c2s_sf_2e1i,
        test_c2s_sf_2e2,
        test_c2s_sf_2e2i,
        test_c2s_si_2e1,
        test_c2s_si_2e1i,
        test_c2s_si_2e2,
    ):
        passed = True
        for li in range(3):
            for lj in range(3):
                for lk in range(3):
                    for ll in range(3):
                        if not fn(li, lj, lk, ll):
                            passed = False
        if passed:
            print(f'{fn.__name__}       passed')
        else:
            print(f'{fn.__name__}       failed')

    for fn in (
        test_c2s_sf_3c2e1,
        test_c2s_sf_3c2e1i,
        test_c2s_si_3c2e1,
        test_c2s_si_3c2e1i,
        test_c2s_sf_3c2e1_ssc,
        test_c2s_sf_3c2e1i_ssc,
        test_c2s_si_3c2e1_ssc,
        test_c2s_si_3c2e1i_ssc,
    ):
        passed = True
        for li in range(3):
            for lj in range(3):
                for lk in range(3):
                    if not fn(li, lj, lk):
                        passed = False
        if passed:
            print(f'{fn.__name__}       passed')
        else:
            print(f'{fn.__name__}       failed')

    for fn in (
        test_c2s_sf_1e_grids,
        test_c2s_sf_1e_gridsi,
        test_c2s_si_1e_grids,
        test_c2s_si_1e_gridsi,
    ):
        passed = True
        for li in range(3):
            for lj in range(3):
                if not fn(li, lj):
                    passed = False
        if passed:
            print(f'{fn.__name__}       passed')
        else:
            print(f'{fn.__name__}       failed')

    for fn in (
        check_ket_spinor_si1,
        check_iket_spinor_si1,
        check_ket_spinor_sf1,
        check_iket_spinor_sf1,
        check_ket_spinor,
        check_iket_spinor,
        check_bra_spinor,
        check_bra_spinor_si,
        check_bra_spinor_sf,
        check_bra_spinor_e1sf,
    ):
        passed = True
        for l in range(12):
            if not fn(l):
                passed = False
                print(f'test failed for {fn.__name__} at l={l}')
        if passed:
            print(f'{fn.__name__}       passed')
        else:
            print(f'{fn.__name__}       failed')
