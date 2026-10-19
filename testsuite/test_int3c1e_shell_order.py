'''
int3c1e integrals with a tight shell and a diffuse shell on a distant atom,
and a tight third shell next to the tight one, used to lose all accuracy when
the diffuse shell was the second (j) shell. All shell orders are compared with
high-precision references (Gaussian product theorem for the overlap, Gauss
transform of 1/r for rinv) and with finite differences for iprinv.
'''

import os
import sys
import ctypes
import itertools
import numpy
import mpmath as mp

_cint = ctypes.CDLL(os.path.abspath(os.path.join(__file__, '../../build/libcint.so')))
_cint.CINTgto_norm.restype = ctypes.c_double
_cint.CINTgto_norm.argtypes = [ctypes.c_int, ctypes.c_double]

PTR_RINV_ORIG = 4
PTR_ENV_START = 20
mp.mp.dps = 30


def make_env(shells, rinv_orig):
    '''one-primitive shells (l, exponent, xyz), each on its own center'''
    env = [0.] * PTR_ENV_START
    env[PTR_RINV_ORIG:PTR_RINV_ORIG+3] = rinv_orig
    atm = []
    bas = []
    for ia, (l, a, xyz) in enumerate(shells):
        atm += [1, len(env), 1, 0, 0, 0]
        env += list(xyz)
        bas += [ia, l, 1, 1, 0, len(env), len(env) + 1, 0]
        env += [a, _cint.CINTgto_norm(l, a)]
    return (numpy.asarray(atm, dtype=numpy.int32), numpy.asarray(bas, dtype=numpy.int32),
            numpy.asarray(env, dtype=numpy.double))


def int3c1e(intor, comp, shells, shls, rinv_orig):
    '''(comp, i, j, k) block of the Cartesian integral for the shell order shls'''
    atm, bas, env = make_env(shells, rinv_orig)
    nf = [(shells[s][0] + 1) * (shells[s][0] + 2) // 2 for s in shls]
    buf = numpy.zeros(comp * nf[0] * nf[1] * nf[2])
    ptr = lambda x: x.ctypes.data_as(ctypes.c_void_p)
    n = len(shells)
    getattr(_cint, intor)(ptr(buf), None, ptr(numpy.asarray(shls, dtype=numpy.int32)),
                          ptr(atm), ctypes.c_int(n), ptr(bas), ctypes.c_int(n), ptr(env), None, None)
    return buf.reshape(comp, nf[2], nf[1], nf[0]).transpose(0, 3, 2, 1)


def cart_powers(l):
    return [(lx, ly, l - lx - ly) for lx in range(l, -1, -1) for ly in range(l - lx, -1, -1)]


def reference(shells, rinv_orig=None):
    '''(i, j, k) block of int3c1e_cart, or of int3c1e_rinv_cart if rinv_orig is given.
    The product of the three Gaussians is one Gaussian exp(-s|r-P|^2); every polynomial
    factor is expanded around P (or around the center Q of exp(-s|r-P|^2 - t^2|r-C|^2)
    in the Gauss transform 1/|r-C| = 2/sqrt(pi) int_0^inf exp(-t^2|r-C|^2) dt).'''
    # libcint puts the angular factors of s and p functions into the integral
    fac_sp = {0: 1 / mp.sqrt(4 * mp.pi), 1: mp.sqrt(3 / (4 * mp.pi))}
    ls = [l for l, a, xyz in shells]
    exps = [mp.mpf(a) for l, a, xyz in shells]
    rs = [[mp.mpf(x) for x in xyz] for l, a, xyz in shells]
    s = sum(exps)
    p = [sum(a * r[d] for a, r in zip(exps, rs)) / s for d in range(3)]
    d2 = lambda u, v: sum((x - y)**2 for x, y in zip(u, v))
    pre = mp.exp(-sum(exps[m] * exps[n] * d2(rs[m], rs[n]) for m, n in ((0, 1), (0, 2), (1, 2))) / s)
    for l, a in zip(ls, exps):
        pre *= fac_sp.get(l, 1) * _cint.CINTgto_norm(l, float(a))
    comps = list(itertools.product(*[enumerate(cart_powers(l)) for l in ls]))

    def moments(q, st):
        '''all products of the three polynomial factors integrated with exp(-st|x-q|^2)'''
        out = {}
        for (x, px), (y, py), (z, pz) in comps:
            v = 1
            for d in range(3):
                key = (d, px[d], py[d], pz[d])
                if key not in out:
                    poly = [mp.mpf(1)]
                    for n, r in zip(key[1:], rs):
                        f = [mp.binomial(n, k) * (q[d] - r[d])**(n - k) for k in range(n + 1)]
                        poly = [sum(poly[i] * f[k - i] for i in range(len(poly)) if 0 <= k - i < len(f))
                                for k in range(len(poly) + len(f) - 1)]
                    out[key] = sum(c * mp.gamma(mp.mpf(k + 1) / 2) / st**(mp.mpf(k + 1) / 2)
                                   for k, c in enumerate(poly) if k % 2 == 0)
                v *= out[key]
            out[x, y, z] = v
        return out

    out = numpy.zeros([len(cart_powers(l)) for l in ls])
    if rinv_orig is None:
        m = moments(p, s)
        for (x, _), (y, _), (z, _) in comps:
            out[x, y, z] = pre * m[x, y, z]
        return out
    c = [mp.mpf(x) for x in rinv_orig]
    cache = {}

    def at_t(t):
        if t not in cache:
            st = s + t * t
            q = [(s * p[d] + t * t * c[d]) / st for d in range(3)]
            m = moments(q, st)
            f = mp.exp(-s * t * t / st * d2(p, c))
            cache[t] = {key: f * v for key, v in m.items()}
        return cache[t]
    for (x, _), (y, _), (z, _) in comps:
        out[x, y, z] = pre * 2 / mp.sqrt(mp.pi) * mp.quad(lambda t: at_t(t)[x, y, z], [0, 1, mp.inf])
    return out


def check(name, v, ref, rtol):
    global failed
    err = abs(v - ref).max()
    if not err < rtol * abs(ref).max():
        failed = True
        print('FAIL: %s max|ref|=%.6g max|err|=%.6g' % (name, abs(ref).max(), err))


def in_order(v, perm):
    '''block computed for the shell order perm, transposed to the order (0, 1, 2)'''
    return v.transpose([0] + [perm.index(a) + 1 for a in range(3)])


# tight f, diffuse s/d 8.6 Bohr away, tight i (l = 6) next to the f shell
cases = [
    [(3, 142., (0., 0., 0.)), (0, .08, (0., 0., 8.6)), (6, 222., (0., 0., .01))],
    [(3, 142., (0., 0., 0.)), (2, .08, (.5, -.3, 8.6)), (6, 222., (.02, 0., 0.))],
]
rinv_origs = [(1., -.5, 4.3), (0., 0., 0.), (0., 0., 8.6)]
perms = list(itertools.permutations(range(3)))

failed = False
for shells in cases:
    ls = [s[0] for s in shells]
    ref = reference(shells)
    for perm in perms:
        v = in_order(int3c1e('int3c1e_cart', 1, shells, perm, (0., 0., 0.)), perm)
        check('int3c1e_cart l=%s order=%s' % (ls, perm), v[0], ref, 1e-12)

    for rinv_orig in rinv_origs:
        # the high-precision rinv reference is slow; for the other origins all orders
        # are compared with the first one
        if shells is cases[0] and rinv_orig is rinv_origs[0]:
            ref = reference(shells, rinv_orig)
        else:
            ref = int3c1e('int3c1e_rinv_cart', 1, shells, (0, 1, 2), rinv_orig)[0]
        for perm in perms:
            v = in_order(int3c1e('int3c1e_rinv_cart', 1, shells, perm, rinv_orig), perm)
            check('int3c1e_rinv_cart l=%s rinv_orig=%s order=%s' % (ls, rinv_orig, perm), v[0], ref, 1e-12)

        # iprinv differentiates the first shell: -d/dR_i of rinv
        h = 1e-5
        for perm in perms:
            i = perm[0]
            fd = []
            for d in range(3):
                plus = [list(s) for s in shells]
                minus = [list(s) for s in shells]
                plus[i][2] = numpy.add(shells[i][2], h * numpy.eye(3)[d])
                minus[i][2] = numpy.add(shells[i][2], -h * numpy.eye(3)[d])
                fd.append((int3c1e('int3c1e_rinv_cart', 1, minus, perm, rinv_orig)[0] -
                           int3c1e('int3c1e_rinv_cart', 1, plus, perm, rinv_orig)[0]) / (2 * h))
            v = int3c1e('int3c1e_iprinv_cart', 3, shells, perm, rinv_orig)
            check('int3c1e_iprinv_cart l=%s rinv_orig=%s order=%s' % (ls, rinv_orig, perm), v, numpy.array(fd), 1e-6)
if failed:
    sys.exit(1)
print('pass: int3c1e shell order')
