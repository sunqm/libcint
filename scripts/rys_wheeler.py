import numpy as np
import mpmath
DECIMALS = 60
mpmath.mp.dps = DECIMALS

def jacobi_alpha(n):
    if n < 0:
        return 0
    half = mpmath.mpf('.5')
    return (2*n*(n+half)-mpmath.mpf('.25')) / ((2*n+mpmath.mpf('1.5'))*(2*n-half))

def jacobi_beta(n):
    if n <= 0:
        return 0
    half = mpmath.mpf('.5')
    return (n*n*(n-half)*(n-half)) / ((2*n-half)**2*(2*n-mpmath.mpf('1.5'))*(2*n+half))

_jacobi_coefs_cache = {}
def jacobi_coefs(n):
    if n in _jacobi_coefs_cache:
        return _jacobi_coefs_cache[n]

    if n == 0:
        coefs = np.array([mpmath.mpf(1)])
    elif n == 1:
        coefs = np.array([-jacobi_alpha(0), mpmath.mpf(1)])
    else:
        coefs2 = jacobi_coefs(n - 2)
        coefs1 = jacobi_coefs(n - 1)
        coefs = (-jacobi_alpha(n - 1)) * coefs1
        coefs[1:] += coefs1[:-1]
        coefs = np.append(coefs, coefs1[-1])
        coefs[:n-1] -= jacobi_beta(n - 1) * coefs2
    _jacobi_coefs_cache[n] = coefs
    return coefs

def fmt1_erfc(t, m, low=0):
    t = mpmath.mpf(t)
    low = mpmath.mpf(low)
    half = mpmath.mpf('.5')
    b = (m + half)
    e = half * mpmath.exp(-t)
    e1 = half * mpmath.exp(-t * low*low) * mpmath.power(low, 2*m+1)
    x = e
    x1 = e1
    s = e - e1
    bi = b
    div = 1
    delta = s
    while abs(delta) > 1e-30:
        bi += 1
        div *= t / bi
        x1 *= low*low
        delta = (x - x1) * div
        s += delta
    f = s / b
    out = [f]
    for i in range(m):
        b -= 1
        if low == 0:
            f = (e + t * f) / b
        else:
            e1 /= low*low
            f = (e - e1 + t * f) / b
        out.append(f)
    return np.array(out)[::-1]

def dble_fmt1_erfc(t, m, low=0):
    half = .5
    b = (m + half)
    e = half * np.exp(-t)
    e1 = half * np.exp(-t * low*low) * np.power(low, 2*m+1)
    x = e
    x1 = e1
    s = e - e1
    bi = b
    div = 1
    delta = s
    while abs(delta) > 1e-30:
        bi += 1
        div *= t / bi
        x1 *= low*low
        delta = (x - x1) * div
        s += delta
    f = s / b
    out = [f]
    for i in range(m):
        b -= 1
        if low == 0:
            f = (e + t * f) / b
        else:
            e1 /= low*low
            f = (e - e1 + t * f) / b
        out.append(f)
    return np.array(out)[::-1]

def naive_jacobi_mus(n, t, low):
    fmt = fmt1_erfc(t, n, low)
    mus = []
    for i in range(n):
        coefs = jacobi_coefs(i)
        mus.append((coefs * fmt[:i+1])[::-1].sum())
    return np.array(mus)

def dble_naive_jacobi_mus(n, t, low):
    fmt = dble_fmt1_erfc(t, n, low)
    mus = []
    for i in range(n):
        coefs = jacobi_coefs(i).astype(float)
        order = abs(coefs).argsort()
        mus.append((coefs[order] * fmt[:i+1][order]).sum())
    return np.array(mus)

def jacobi_gs(n, x):
    g0 = mpmath.mpf(0)
    g1 = mpmath.mpf(1)
    gs = [g1]
    for i in range(1, n):
        g2 = (x - jacobi_alpha(i-1)) * g1 - jacobi_beta(i-1) * g0
        gs.append(g2)
        g1, g0 = g2, g1
    return gs

def jacobi_x(n):
    half = mpmath.mpf('.5')
    A = (n + half) / (2 * n + half + 1)
    B = (2*n + half) / (n + 1)
    return A / B

def jacobi_zs(n, low, t):
    x = low ** 2
    gx = np.array(jacobi_gs(n+1, x))
    xs = np.array([jacobi_x(i) for i in range(n)])
    et = mpmath.exp(-t*x)*mpmath.sqrt(x)/(2*t)
    return et * (xs*gx[:-1]-gx[1:])

def jacobi_us(n, low, t):
    half = mpmath.mpf('.5')
    zs = jacobi_zs(n, low, t)
    ns = np.arange(n)
    As = (ns + half) / (2 * ns + half + 1)
    Bs = (2*ns + half) / (ns + 1)
    xs = As / Bs
    cs = xs * 2*ns / (2*ns + 1)

    u0 = ((-1/(2*t) + 3/(4*t*(half+1)))*mpmath.exp(-t) -
          (-low**3/(2*t) + 3*low/(4*t*(half+1)))*mpmath.exp(-t*low**2))
    us = [u0]
    for i in range(1, n):
        u1 = -zs[i] - u0 * cs[i]
        us.append(u1)
        u0 = u1
    return us

# rn = (n+.5)/t + .5 * (2*n+1) / (2*(2*n+1.5)*(2*n-.5))
def jacobi_rn_part(n):
    half = mpmath.mpf('.5')
    return half * (2*n+1) / (2*(2*n+mpmath.mpf('1.5'))*(2*n-half))

def jacobi_sn(n):
    half = mpmath.mpf('.5')
    return mpmath.mpf(n)*(2*n+1)*(2*n-1)*(n-half)/(4*(2*n+half)*(2*n-half)**2*(2*n-mpmath.mpf(1.5)))

def jacobi_mus(n, t, low):
    t = mpmath.mpf(t)
    low = mpmath.mpf(low)
    half = mpmath.mpf('.5')
    n_idx = mpmath.mpf(1) * np.arange(1, n)
    rn = (2*n_idx+1)/(2*t) - (2*n_idx+1)*(half-1) / (2*(2*n_idx+half+1)*(2*n_idx+half-1))
    sn = n_idx*(2*n_idx+1)*(2*n_idx-1)*(n_idx+half-1)/(4*(2*n_idx+half)*(2*n_idx+half-1)**2*(2*n_idx+half-2))
    us = jacobi_us(n, low, t)

    low2 = low * low
    t_inv = half / t
    e0 = mpmath.exp(-t) * t_inv
    et = mpmath.exp(-t * low2) * low * t_inv
    tt = mpmath.sqrt(t)
    alpha0 = jacobi_alpha(0)
    mu0 = mpmath.sqrt(mpmath.pi) / 2 / tt * (mpmath.erfc(low * tt) - mpmath.erfc(tt))
    mu1 = t_inv * mu0 - e0 + et - alpha0 * mu0;

    mus = [mu0, mu1]
    for i in range(1, n-1):
        mu2 = rn[i-1] * mu1 + sn[i-1] * mu0 + us[i-1]
        mus.append(mu2)
        mu1, mu0 = mu2, mu1
    return np.array(mus)

def shifted_jacobi_moments(n, t, low=0):
    if 0:
        mus = jacobi_mus(n, t, low)
    elif 1:
        mus = naive_jacobi_mus(n, t, low)
    else:
        mus = dble_naive_jacobi_mus(n, t, low)
    alphas = np.array([jacobi_alpha(i) for i in range(n)])
    betas = np.array([jacobi_beta(i) for i in range(n)])
    return alphas, betas, mus

def laguerre_mu0(t, low=0):
    if t == 0:
        return 1 - low
    else:
        tt = mpmath.sqrt(t)
        return mpmath.sqrt(mpmath.pi) / 2 / tt * (mpmath.erfc(low * tt) - mpmath.erfc(tt))

def t_scaled_laguerre(T, a, n, x=1):
    l0 = mpmath.mpf(1)
    l1 = mpmath.mpf(x) - (1 + a) / T
    ls = [l0, l1]
    for i in range(2, n):
        l2 = (x - (2*i + a - 1) / T) * l1 - (i - 1) * (i + a - 1) / T**2 * l0
        ls.append(l2)
        l1, l0 = l2, l1
    return np.array(ls)

def laguerre_moments(n, T, low=0):
    half = mpmath.mpf('.5')
    eta = low**2
    moments = (t_scaled_laguerre(T, half, n-1) * mpmath.exp(-T) / (-2 * T) -
               t_scaled_laguerre(T, half, n-1, eta) * mpmath.sqrt(eta) * mpmath.exp(-T * eta) / (-2 * T))
    moments = np.append(laguerre_mu0(T, low), moments)
    idx = mpmath.mpf(1) * np.arange(1, n)
    alpha = (idx * 4 - 3) / (2 * T)
    beta = (idx - 1) * (idx * 2 - 3) / (2 * T**2)
    return alpha, beta, moments

def flocke_jacobi_moments(n, T):
    import scipy.special
    # Flocke's recipe JCP, 131, 064107
    mu1 = 1
    mu2 = 0
    moments = [mu2, mu1]
    additive_for_dp = 20
    # Miller algorithm
    for i in reversed(range(1, n+1+additive_for_dp)):
        r = (2 * i + 1) / (2*T) + (2 * i + 1) / ((4 * i - 1) * (4 * i + 3))
        s = (2 * i * (2 * i + 1) * (2 * i - 1)**2 /
             ((4 * i - 3) * (4 * i + 1) * (4 * i - 1)**2))
        mu0 = (mu2 - r * mu1) / s
        moments.append(mu0)
        mu1, mu2 = mu0, mu1
        #if abs(mu0) > 1e8:
        #    # TODO: scale down moments
        #    raise RuntimeError

    moments = np.array(moments)[::-1]
    if T == 0:
        fmt0 = 1
    else:
        tt = np.sqrt(T)
        fmt0 = np.pi**.5/2. / tt * scipy.special.erf(tt)
    moments = moments[:n] * (fmt0 / mu0)
    idx = np.arange(n - 1)
    alpha = (2 * idx * (idx + .5) - .25) / ((2 * idx + .5)**2 - 1)
    beta = (idx**2 * (idx - .5)**2) / ((2 * idx - .5)**2 * (2*idx-1.5) * (2*idx+.5))
    return alpha, beta, moments

def wheeler(n, alpha, beta, moments):
    sig_m = moments
    sig_0 = moments
    a0 = mpmath.mpf(alpha[0]) + moments[1] / moments[0]
    b0 = mpmath.mpf(0)
    val_a = [a0]
    val_b = [b0]
    for k in range(1, n):
        nc = 2 * n - 2 * k
        sig_k = (sig_0[2:2+nc] - (a0 - alpha[k:k+nc]) * sig_0[1:1+nc] -
                 b0 * sig_m[2:2+nc] + beta[k:k+nc] * sig_0[:nc])
        a1 = alpha[k] - sig_0[1] / sig_0[0] + sig_k[1] / sig_k[0]
        b1 = sig_k[0] / sig_0[0]
        val_a.append(a1)
        val_b.append(b1)
        a0, b0 = a1, b1
        sig_0, sig_m = sig_k, sig_0
    return np.array(val_a), np.array(val_b)

def roots_and_weights_partial(n, alpha, beta, moments):
    a, b = wheeler(n, alpha, beta, moments)
    Tmat = np.diag(a.astype(float))
    idx = np.arange(n - 1)
    Tmat[idx, idx+1] = Tmat[idx+1, idx] = b[1:].astype(float)**.5

    roots, c = np.linalg.eigh(Tmat)
    weights = c[0]**2 * float(moments[0])
    return roots, weights

def roots_and_weights(n, x, low=0):
    x = mpmath.mpf(x)
    low = mpmath.mpf(low)
    if x > 20:
        alpha, beta, moments = laguerre_moments(n*2, x, low)
    #elif low == 0:
    #    alpha, beta, moments = flocke_jacobi_moments(n*2, x)
    else:
        alpha, beta, moments = shifted_jacobi_moments(n*2, x, low)
    roots, weights = roots_and_weights_partial(n, alpha, beta, moments)
    roots = roots / (1 - roots)
    return roots, weights

if __name__ == '__main__':
    #for i in range(64):
    #    mpmath.nprint(jacobi_alpha(i), 36)
    #    mpmath.nprint(jacobi_beta(i), 36)
    #    mpmath.nprint(jacobi_x(i), 17)
    #    mpmath.nprint(2*i / (2*i + 1) * jacobi_x(i), 17)
    #    mpmath.nprint(jacobi_rn_part(i+1), 36)
    #    mpmath.nprint(jacobi_sn(i+1), 36)
    #gs = jacobi_gs(49, 1)
    #for i in gs:
    #    mpmath.nprint(i, 17)
    #for i in range(64):
    #    print(f'// n = {i}')
    #    for c in jacobi_coefs(i):
    #        mpmath.nprint(c, 36)
    #for i in range(64):
    #    cs = jacobi_coefs(i)
    #    print(', '.join([str(i) for i in abs(cs.astype(float)).argsort()]))

    print(roots_and_weights_partial(2, *laguerre_moments(2*2, 5.7, 0)))
    print(roots_and_weights_partial(2, *flocke_jacobi_moments(2*2, 5.7)))
    print(roots_and_weights_partial(2, *laguerre_moments(2*2, 5.7, 0.1)))
    print(roots_and_weights_partial(2, *shifted_jacobi_moments(2*2, 5.7, 0.1)))

    print(roots_and_weights_partial(3, *laguerre_moments(3*2, 11.7, 0)))
    print(roots_and_weights_partial(3, *flocke_jacobi_moments(3*2, 11.7)))
    print(roots_and_weights_partial(3, *laguerre_moments(3*2, 11.7, 0.2)))
    print(roots_and_weights_partial(3, *shifted_jacobi_moments(3*2, 11.7, 0.2)))

    print(roots_and_weights_partial(4, *laguerre_moments(4*2, 2.7, 0)))
    print(roots_and_weights_partial(4, *flocke_jacobi_moments(4*2, 2.7)))
    print(roots_and_weights_partial(4, *laguerre_moments(4*2, 2.7, 0.2)))
    print(roots_and_weights_partial(4, *shifted_jacobi_moments(4*2, 2.7, 0.2)))

    print(roots_and_weights_partial(5, *laguerre_moments(5*2, 1.7, 0)))
    print(roots_and_weights_partial(5, *flocke_jacobi_moments(5*2, 1.7)))
    print(roots_and_weights_partial(5, *laguerre_moments(5*2, 1.7, 0.2)))
    print(roots_and_weights_partial(5, *shifted_jacobi_moments(5*2, 1.7, 0.2)))

    print(roots_and_weights_partial(6, *laguerre_moments(6*2, 1.1, 0)))
    print(roots_and_weights_partial(6, *flocke_jacobi_moments(6*2, 1.1)))
    print(roots_and_weights_partial(6, *laguerre_moments(6*2, 1.1, 0.2)))
    print(roots_and_weights_partial(6, *shifted_jacobi_moments(6*2, 1.1, 0.2)))
