import mpmath
DECIMALS = mpmath.mp.dps
import numpy as np

METHOD = 'eig'

def find_polyroots(cs, nroots):
    if nroots == 1:
        return np.array([-cs[1,0] / cs[1,1]])

    elif nroots == 2:
        dum = mpmath.sqrt(cs[2,1]**2 - 4 * cs[2,0] * cs[2,2])
        rt0 = (-cs[2,1] - dum) / cs[2,2] / 2
        rt1 = (-cs[2,1] + dum) / cs[2,2] / 2
        return np.array([rt0, rt1])

    if METHOD == 'bisection':
        return bisection(cs, nroots)

    A = mpmath.matrix(nroots)
    for m in range(nroots-1):
        A[m+1,m] = mpmath.mpf(1)
    for m in range(nroots):
        A[0,nroots-1-m] = -cs[nroots,m] / cs[nroots,nroots]
    roots = eig(A)
    return np.array(roots[::-1])

def bisection(cs, nroots):
    rt = np.array([mpmath.mpf(1)] * (nroots + 1))
    if nroots == 1:
        #rt[0] = ff[1] / ff[0]
        rt[0] = -cs[1,0] / cs[1,1]
    else:
        dum = mpmath.sqrt(cs[2,1]**2 - 4 * cs[2,0] * cs[2,2])
        rt[0] = (-cs[2,1] - dum) / cs[2,2] / 2
        rt[1] = (-cs[2,1] + dum) / cs[2,2] / 2

    for k in range(2, nroots):
        R_dnode(cs[k+1,:k+2], rt)
    return rt

ACCRT = min(mpmath.power(.1, DECIMALS/2), mpmath.mpf('1e-35'))
def R_dnode(a, rt):
    order = a.size - 1
    quat1 = mpmath.mpf('.25')
    quat3 = mpmath.mpf('.75')

    x1init = mpmath.mpf(0)
    p1init = mpmath.mpf(a[0])
    for m in range(order):
        x0 = x1init
        p0 = p1init
        x1init = mpmath.mpf(rt[m])
        p1init = poly_value1(a, order, x1init)
        if (p0 * p1init > 0):
            raise RuntimeError

        if (x0 <= x1init):
            x1 = x1init
            p1 = p1init
        else:
            x1 = x0
            p1 = p0
            x0 = x1init
            p0 = p1init

        if (p0 == 0):
            rt[m] = x0
            continue
        elif (p1 == 0):
            rt[m] = x1
            continue
        else:
            xi = x0 + (x0 - x1) / (p1 - p0) * p0

        n = 0
        while (x1 > ACCRT+x0) or (x0 > x1+ACCRT):
            n += 1
            if (n > 1000):
                raise RuntimeError

            pi = poly_value1(a, order, xi)
            if (-ACCRT < pi < ACCRT):
                break
            elif (p0 * pi <= 0):
                x1 = xi
                p1 = pi
                xi = x0 * quat1 + xi * quat3
            else:
                x0 = xi
                p0 = pi
                xi = xi * quat3 + x1 * quat1

            pi = poly_value1(a, order, xi)
            if (-ACCRT < pi < ACCRT):
                break
            elif (p0 * pi <= 0):
                x1 = xi
                p1 = pi
            else:
                x0 = xi
                p0 = pi

            xi = x0 + (x0 - x1) / (p1 - p0) * p0
        rt[m] = xi
    return rt

def poly_value1(a, order, x):
    p = a[order]
    for i in range(1, order+1):
        p = p * x + a[order-i]
    return p

def qr_step(n0, n1, A, shift):
    """
    This subroutine executes a single implicitly shifted QR step applied to an
    upper Hessenberg matrix A. Given A and shift as input, first an QR
    decomposition is calculated:

      Q R = A - shift * 1 .

    The output is then following matrix:

      R Q + shift * 1

    parameters:
      n0, n1    (input) Two integers which specify the submatrix A[n0:n1,n0:n1]
                on which this subroutine operators. The subdiagonal elements
                to the left and below this submatrix must be deflated (i.e. zero).
                following restriction is imposed: n1>=n0+2
      A         (input/output) On input, A is an upper Hessenberg matrix.
                On output, A is replaced by "R Q + shift * 1"
      shift     (input) a complex number specifying the shift. idealy close to an
                eigenvalue of the bottemmost part of the submatrix A[n0:n1,n0:n1].

    references:
      Stoer, Bulirsch - Introduction to Numerical Analysis.
      Kresser : Numerical Methods for General and Structured Eigenvalue Problems
    """

    # implicitly shifted and bulge chasing is explained at p.398/399 in "Stoer, Bulirsch - Introduction to Numerical Analysis"
    # for bulge chasing see also "Watkins - The Matrix Eigenvalue Problem" sec.4.5,p.173

    # the Givens rotation we used is determined as follows: let c,s be two complex
    # numbers. then we have following relation:
    #
    #     v = sqrt(|c|^2 + |s|^2)
    #
    #     1/v [ c~  s~]  [c] = [v]
    #         [-s   c ]  [s]   [0]
    #
    # the matrix on the left is our Givens rotation.

    n = A.rows

    # first step

    # calculate givens rotation
    c = A[n0  ,n0] - shift
    s = A[n0+1,n0]

    v = mpmath.sqrt(c**2 + s**2)

    if v == 0:
        v = 1
        c = 1
        s = 0
    else:
        c /= v
        s /= v

    for k in range(n0, n):
        # apply givens rotation from the left
        x = A[n0  ,k]
        y = A[n0+1,k]
        A[n0  ,k] = c * x + s * y
        A[n0+1,k] = c * y - s * x

    for k in range(min(n1, n0+3)):
        # apply givens rotation from the right
        x = A[k,n0  ]
        y = A[k,n0+1]
        A[k,n0  ] = c * x + s * y
        A[k,n0+1] = c * y - s * x

    # chase the bulge

    for j in range(n0, n1 - 2):
        # calculate givens rotation

        c = A[j+1,j]
        s = A[j+2,j]

        v = mpmath.sqrt(c**2 + s**2)
        A[j+1,j] = v
        A[j+2,j] = 0

        if v == 0:
            v = 1
            c = 1
            s = 0
        else:
            c /= v
            s /= v

        for k in range(j+1, n):
            # apply givens rotation from the left
            x = A[j+1,k]
            y = A[j+2,k]
            A[j+1,k] = c * x + s * y
            A[j+2,k] = c * y - s * x

        for k in range(min(n1, j+4)):
            # apply givens rotation from the right
            x = A[k,j+1]
            y = A[k,j+2]
            A[k,j+1] = c * x + s * y
            A[k,j+2] = c * y - s * x

def hessenberg_qr(A, eps, maxits=120):
    n = A.rows
    n0 = 0
    n1 = n
    its = 0
    while 1:
        # kressner p.32 algo 3
        # the active submatrix is A[n0:n1,n0:n1]

        k = n0

        while k + 1 < n1:
            s = abs(A[k,k]) + abs(A[k+1,k+1])
            #if s < eps:
            #    s = 1
            # Ensure relative error converged
            if abs(A[k+1,k]) < eps * s:
                break
            k += 1

        if k + 1 < n1:
            # deflation found at position (k+1, k)

            A[k+1,k] = 0
            n0 = k + 1
            its = 0

            if n0 + 1 >= n1:
                # block of size at most two has converged
                n0 = 0
                n1 = k + 1
                if n1 < 2:
                    # QR algorithm has converged
                    return
        else:
            #    A = [ a b ]       det(x-A)=x*x-x*tr(A)+det(A)
            #        [ c d ]
            #
            # eigenvalues bad:   (tr(A)+sqrt((tr(A))**2-4*det(A)))/2
            #     bad because of cancellation if |c| is small and |a-d| is small, too.
            #
            # eigenvalues good:     (a+d+sqrt((a-d)**2+4*b*c))/2

            t =  A[n1-2,n1-2] + A[n1-1,n1-1]
            s = (A[n1-1,n1-1] - A[n1-2,n1-2]) ** 2 + 4 * A[n1-1,n1-2] * A[n1-2,n1-1]
            if s >= 0:
                s = mpmath.sqrt(s)
                a = (t + s) / 2
                b = (t - s) / 2
                if abs(A[n1-1,n1-1] - a) > abs(A[n1-1,n1-1] - b):
                    shift = b
                else:
                    shift = a
            else:
                if n1 == 2:
                    raise RuntimeError("qr: failed to find real roots")
                shift = t / 2

            its += 1

            qr_step(n0, n1, A, shift)

            if its > maxits:
                raise RuntimeError("qr: failed to converge after %d steps" % its)

def eig(A):
    '''
    Modified based on mpmath.matrices.eigen.eig function.
    This implementation restricts the eigenvalues to real.
    '''
    n = A.rows

    if n == 1:
        return [A[0,0]]

    eps = mpmath.mp.eps / 10
    hessenberg_qr(A, eps)

    E = [A[i,i] for i in range(n)]

    return E
