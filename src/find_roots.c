/*
    Modified based on mpmath.matrices.eigen.eig function.
    This implementation restricts the eigenvalues to real.
*/

#include <stdio.h>
#include <math.h>

#define SQUARE(X) ((X) * (X))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MXRYSROOTS 32

#define POLYNOMIAL_VALUE1(p, a, order, x) \
p = a[order]; \
for (i = 1; i <= order; i++) { \
        p = p * x + a[order-i]; \
}

static int R_dnode(double *a, double *roots, int order)
{
        const double accrt = 1e-15;
        double x0, x1, xi, x1init, p0, p1, pi, p1init;
        int i, m, n;

        x1init = 0;
        p1init = a[0];
        for (m = 0; m < order; ++m) {
                x0 = x1init;
                p0 = p1init;
                x1init = roots[m];
                POLYNOMIAL_VALUE1(p1init, a, order, x1init);

                // When all coefficients a are 0, short-circuit the rest code to
                // ensure the roots from the lower order polynomials are preserved
                if (p1init == 0) {
                        // roots[m] = x1init;
                        continue;
                }
                if (p0 * p1init > 0) {
                        fprintf(stderr, "ROOT NUMBER %d WAS NOT FOUND FOR POLYNOMIAL OF ORDER %d\n",
                                m, order);
                        return 1;
                }
                if (x0 <= x1init) {
                        x1 = x1init;
                        p1 = p1init;
                } else {
                        x1 = x0;
                        p1 = p0;
                        x0 = x1init;
                        p0 = p1init;
                }
                // interpolate/extrapolate between [x0,x1]
                if (p1 == 0) {
                        roots[m] = x1;
                        continue;
                } else if (p0 == 0) {
                        roots[m] = x0;
                        continue;
                } else {
                        xi = x0 + (x0 - x1) / (p1 - p0) * p0;
                }
                n = 0;
                while (fabs(x1 - x0) > x1*accrt) {
                        n++;
                        if (n > 200) {
                                fprintf(stderr, "libcint::rys_roots NO CONV. IN R_dnode\n");
                                return 1;
                        }
                        POLYNOMIAL_VALUE1(pi, a, order, xi);
                        if (pi == 0) {
                                break;
                        } else if (p0 * pi <= 0) {
                                x1 = xi;
                                p1 = pi;
                                xi = x0 * .25 + xi * .75;
                        } else {
                                x0 = xi;
                                p0 = pi;
                                xi = xi * .75 + x1 * .25;
                        }
                        POLYNOMIAL_VALUE1(pi, a, order, xi);
                        if (pi == 0) {
                                break;
                        } else if (p0 * pi <= 0) {
                                x1 = xi;
                                p1 = pi;
                        } else {
                                x0 = xi;
                                p0 = pi;
                        }

                        xi = x0 + (x0 - x1) / (p1 - p0) * p0;
                }
                roots[m] = xi;
        }
        return 0;
}

static void _qr_step(double *A, int nroots, int n0, int n1, double shift)
{
        int m1 = n0 + 1;
        int j, k, m3, j1, j2;
        double c = A[n0*nroots+n0] - shift;
        double s = A[m1*nroots+n0];
        double v = sqrt(c*c + s*s);
        double x, y;

        if (v == 0) {
                v = 1;
                c = 1;
                s = 0;
        }
        v = 1. / v;
        c *= v;
        s *= v;

        for (k = n0; k < nroots; k++) {
                // apply givens rotation from the left
                x = A[n0*nroots+k];
                y = A[m1*nroots+k];
                A[n0*nroots+k] = c * x + s * y;
                A[m1*nroots+k] = c * y - s * x;
        }

        m3 = MIN(n1, n0+3);
        for (k = 0; k < m3; k++) {
                // apply givens rotation from the right
                x = A[k*nroots+n0];
                y = A[k*nroots+m1];
                A[k*nroots+n0] = c * x + s * y;
                A[k*nroots+m1] = c * y - s * x;
        }

        for (j = n0; j < n1 - 2; j++) {
                j1 = j + 1;
                j2 = j + 2;
                // calculate givens rotation
                c = A[j1*nroots+j];
                s = A[j2*nroots+j];
                v = sqrt(c*c + s*s);
                A[j1*nroots+j] = v;
                A[j2*nroots+j] = 0;

                if (v == 0) {
                        v = 1;
                        c = 1;
                        s = 0;
                }
                v = 1. / v;
                c *= v;
                s *= v;

                for (k = j1; k < nroots; k++) {
                        // apply givens rotation from the left
                        x = A[j1*nroots+k];
                        y = A[j2*nroots+k];
                        A[j1*nroots+k] = c * x + s * y;
                        A[j2*nroots+k] = c * y - s * x;
                }
                m3 = MIN(n1, j+4);
                for (k = 0; k < m3; k++) {
                        // apply givens rotation from the right
                        x = A[k*nroots+j1];
                        y = A[k*nroots+j2];
                        A[k*nroots+j1] = c * x + s * y;
                        A[k*nroots+j2] = c * y - s * x;
                }
        }
}

static int _hessenberg_qr(double *A, int nroots)
{
        double eps = 1e-15;
        int maxits = 30;
        int n0 = 0;
        int n1 = nroots;
        int its = 0;
        int k, ic, k1;
        for (ic = 0; ic < nroots*maxits; ic++) {
                k = n0;
                while (k + 1 < n1) {
                        double s = fabs(A[k*nroots+k]) + fabs(A[(k+1)*nroots+k+1]);
                        if (fabs(A[(k+1)*nroots+k]) < eps * s) {
                                break;
                        }
                        k += 1;
                }

                k1 = k + 1;
                if (k1 < n1) {
                        // deflation found at position (k+1, k)
                        A[k1*nroots+k] = 0;
                        n0 = k1;
                        its = 0;

                        if (n0 + 1 >= n1) {
                                // block of size at most two has converged
                                n0 = 0;
                                n1 = k1;
                                if (n1 < 2) {
                                        // QR algorithm has converged
                                        return 0;
                                }
                        }
                } else {
                        int m1 = n1 - 1;
                        int m2 = n1 - 2;
                        double a11 = A[m1*nroots+m1];
                        double a22 = A[m2*nroots+m2];
                        double shift;
                        double t = a11 + a22;
                        double s = SQUARE(a11 - a22);
                        s += 4 * A[m1*nroots+m2] * A[m2*nroots+m1];
                        if (s > 0) {
                                s = sqrt(s);
                                double a = (t + s) * .5;
                                double b = (t - s) * .5;
                                if (fabs(a11 - a) > fabs(a11 - b)) {
                                        shift = b;
                                } else {
                                        shift = a;
                                }
                        } else {
                                if (n1 == 2) {
                                        fprintf(stderr, "hessenberg_qr: failed to find real roots\n");
                                        return 1;
                                }
                                shift = t * .5;
                        }
                        its += 1;
                        _qr_step(A, nroots, n0, n1, shift);
                        if (its > maxits) {
                                fprintf(stderr, "hessenberg_qr: failed to converge after %d steps\n", its);
                                return 1;
                        }
                }
        }
        fprintf(stderr, "hessenberg_qr failed\n");
        return 1;
}

int _CINT_polynomial_roots(double *roots, double *cs, int nroots)
{
        if (nroots == 1) {
                roots[0] = -cs[2] / cs[3];
                return 0;
        } else if (nroots == 2) {
                double dum = sqrt(SQUARE(cs[2*3+1]) - 4*cs[2*3+0]*cs[2*3+2]);
                roots[0] = (-cs[2*3+1] - dum) / cs[2*3+2] / 2;
                roots[1] = (-cs[2*3+1] + dum) / cs[2*3+2] / 2;
                return 0;
        }

        double A[MXRYSROOTS * MXRYSROOTS];
        int nroots1 = nroots + 1;
        // reuse the buffer in coefficients
        int i;
        double fac = -1. / cs[nroots*nroots1+nroots];
        for (i = 0; i < nroots; i++) {
                A[nroots-1-i] = cs[nroots*nroots1+i] * fac;
        }
        for (i = nroots; i < nroots*nroots; i++) {
                A[i] = 0;
        }
        for (i = 0; i < nroots-1; i++) {
                A[(i+1)*nroots+i] = 1.;
        }
        int err = _hessenberg_qr(A, nroots);
        if (err == 0) {
                for (i = 0; i < nroots; i++) {
                        roots[nroots-1-i] = A[i*nroots+i];
                }
        } else {
                int k, order;
                double *a;
                double dum = sqrt(cs[2*nroots1+1] * cs[2*nroots1+1]
                                  - 4 * cs[2*nroots1+0] * cs[2*nroots1+2]);
                roots[0] = .5 * (-cs[2*nroots1+1] - dum) / cs[2*nroots1+2];
                roots[1] = .5 * (-cs[2*nroots1+1] + dum) / cs[2*nroots1+2];
                for (i = 2; i < nroots; i++) {
                        roots[i] = 1;
                }
                for (k = 2; k < nroots; ++k) {
                        order = k + 1;
                        a = cs + order * nroots1;
                        err = R_dnode(a, roots, order);
                        if (err) {
                                break;
                        }
                }
        }
        return err;
}
