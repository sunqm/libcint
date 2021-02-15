/*
 * Copyright (C) 2020-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Eigen solver. Edited based on subroutine DSTEMR from LAPACK 3.9.0
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include "rys_roots.h"

#ifdef LAPACK_FOUND

void dstemr_(
    char const* jobz, char const* range,
    int const* n,
    double* D,
    double* E,
    double const* vl,
    double const* vu, int const* il, int const* iu, int* m,
    double* W,
    double* Z, int const* ldz, int const* nzc, int* ISUPPZ, int* tryrac,
    double* work, int const* lwork,
    int* iwork, int const* liwork,
    int* info);

int _CINTdiagonalize(int n, double *diag, double *diag_off1, double *eig, double *vec)
{
        const char JOBZ = 'V';
        const char RANGE = 'A';
        double VL, VU;
        int IL, IU, M;
        int LDZ = n;
        int NZC = n;
        int ISUPPZ[MXRYSROOTS * 2];
        int TRYRAC = 1;
        double WORK[MXRYSROOTS * 18];
        int LWORK = MXRYSROOTS * 18;
        int IWORK[MXRYSROOTS * 10];
        int LIWORK = MXRYSROOTS * 10;
        int INFO = 0;
        dstemr_(&JOBZ, &RANGE, &n, diag, diag_off1, &VL, &VU, &IL, &IU, &M,
                eig, vec, &LDZ, &NZC, ISUPPZ, &TRYRAC,
                WORK, &LWORK, IWORK, &LIWORK, &INFO);
        return INFO;
}

#else

#define MAX(x, y)       ((x) > (y) ? (x) : (y))
#define MIN(x, y)       ((x) < (y) ? (x) : (y))
#define MAXTRY          6
#define MAXRQITER       6
#define ITMAX           1000
#define FAIL_THRESHOLD  1e16
#define RTOL1           1.5e-8
#define RTOL2           7.45e-11
#define MINGAP          0.001


static int _dlarrk(int n, int iw, double gl, double gu, double *diag, double *e2,
                   double reltol, double *w, double *werr)
{
    int i, it;
    double mid, tmp1, left, right;
    double tnorm;
    int negcnt;
    int info;

    if (n <= 0) {
        return 0;
    }

    tnorm = MAX(fabs(gl), fabs(gu));
    info = -1;
    left = gl - tnorm * 2. * DBL_EPSILON * n;
    right = gu + tnorm * 2. * DBL_EPSILON * n;

    for (it = 0; it < ITMAX; it++) {
        tmp1 = fabs(right - left);
        if (tmp1 <= 0 ||
            tmp1 < reltol * MAX(fabs(right), fabs(left))) {
            info = 0;
            break;
        }

    /*     Count number of negative pivots for mid-point */

        mid = (left + right) * .5;
        negcnt = 0;
        tmp1 = diag[0] - mid;
        if (tmp1 <= 0.) {
            ++negcnt;
        }

        for (i = 1; i < n; ++i) {
            tmp1 = diag[i] - e2[i - 1] / tmp1 - mid;
            if (tmp1 <= 0.) {
                ++negcnt;
            }
        }
        if (negcnt >= iw) {
            right = mid;
        } else {
            left = mid;
        }
    }

/*     Converged or maximum number of iterations reached */

    *w = (left + right) * .5;
    *werr = fabs(right - left) * .5;
    return info;
}


static void _dlarrc(int n, double vl, double vu,
                    double *diag, double *e2, int *lcnt, int *rcnt)
{
    int i;
    int left_count = 0;
    int right_count = 0;
    double tmp, lpivot, rpivot;

    lpivot = diag[0] - vl;
    rpivot = diag[0] - vu;
    if (lpivot <= 0.) {
        ++left_count;
    }
    if (rpivot <= 0.) {
        ++right_count;
    }
    for (i = 0; i < n - 1; ++i) {
        tmp = e2[i];
        lpivot = diag[i + 1] - vl - tmp / lpivot;
        rpivot = diag[i + 1] - vu - tmp / rpivot;
        if (lpivot <= 0.) {
            ++left_count;
        }
        if (rpivot <= 0.) {
            ++right_count;
        }
    }
    *lcnt = left_count;
    *rcnt = right_count;
}

static int _dlasq4(int i0, int n0, int n0init,
                   double *qvecp, double *qvec1p, double *evecp, double *evec1p,
                   double *dmin, double *dn, double *tau)
{
/* an approximation to the smallest eigenvalue using values of d from the previous transform. */
    double a2, b1, b2;
    double gap1, gap2;
    double s = 0.;
    int i;

    if (n0init == n0) {

/*        No eigenvalues deflated. */

        if (dmin[0] == dn[0]) {

/*           Cases 2 and 3. */

            if (dmin[1] == dn[1]) {
                b1 = sqrt(qvecp[n0 - 1] * evecp[n0 - 2]);
                b2 = sqrt(qvecp[n0 - 2] * evecp[n0 - 3]);
                a2 = qvecp[n0 - 2] + evecp[n0 - 2];
                gap2 = dmin[2] - a2 - dmin[2] * .25;
                if (gap2 > b2) {
                    gap1 = a2 - dn[0] - b2 / gap2 * b2;
                } else {
                    gap1 = a2 - dn[0] - (b1 + b2);
                }
                if (gap1 > b1) {
                    s = MAX(dn[0] - b1 / gap1 * b1, dmin[0] * .5);
                } else {
                    s = 0.;
                    if (dn[0] > b1) {
                        s = dn[0] - b1;
                    }
                    if (a2 > b1 + b2) {
                        s = MIN(s, a2 - (b1 + b2));
                    }
                    s = MAX(s, dmin[0] * .333);
                }

            } else {

/*              Case 4. */

                if (evecp[n0 - 2] > qvecp[n0 - 2]) {
                    return 0;
                }

                b2 = evecp[n0 - 2] / qvecp[n0 - 2];
                a2 = b2;

/*              Approximate contribution to norm squared from I < NN-1. */

                for (i = n0 - 3; i >= i0; i--) {
                    b1 = b2;
                    if (evecp[i] > qvecp[i]) {
                        return 0;
                    }
                    b2 *= evecp[i] / qvecp[i];
                    a2 += b2;
                    if (.563 < a2 || MAX(b1, b2) < a2 * .01) {
                        break;
                    }
                }
                a2 *= 1.05;

/*              Rayleigh quotient residual bound. */

                if (a2 < .563) {
                    s = dn[0] * (1. - sqrt(a2)) / (a2 + 1.);
                } else {
                    s = dmin[0] * .25;
                }
            }

        } else if (dmin[0] == dn[1]) {

/*              Case 4. */

            if (evec1p[n0 - 2] > qvec1p[n0 - 1] || evecp[n0 - 3] > qvecp[n0 - 3]) {
                return 0;
            }
            a2 = evec1p[n0 - 2] / qvec1p[n0 - 1];
            b2 = evecp[n0 - 3] / qvecp[n0 - 3];

/*              Approximate contribution to norm squared from I < NN-1. */

            a2 += b2;
            for (i = n0 - 4; i >= i0; i--) {
                if (b2 == 0.) {
                    break;
                }
                b1 = b2;
                if (evecp[i] > qvecp[i]) {
                    return 0;
                }
                b2 *= evecp[i] / qvecp[i];
                a2 += b2;
                if (MAX(b2, b1) * 100. < a2 || .563 < a2) {
                    break;
                }
            }
            a2 *= 1.05;

/*              Rayleigh quotient residual bound. */

            if (a2 < .563) {
                s = dn[1] * (1. - sqrt(a2)) / (a2 + 1.);
            } else {
                s = dmin[0] * .25;
            }

        } else if (dmin[0] == dn[2]) {

/*           Case 5. */

            if (evec1p[n0 - 3] > qvec1p[n0 - 2] || evec1p[n0 - 2] > qvec1p[n0 - 1]) {
                return 0;
            }
            a2 = evec1p[n0 - 3] / qvec1p[n0 - 2] * (evec1p[n0 - 2] / qvec1p[n0 - 1] + 1.);

/*           Approximate contribution to norm squared from I < NN-2. */

            if (n0 - i0 > 3) {
                b2 = evecp[n0 - 4] / qvecp[n0 - 4];
                a2 += b2;
                for (i = n0 - 5; i >= i0; i--) {
                    b1 = b2;
                    if (evecp[i] > qvecp[i]) {
                        return 0;
                    }
                    b2 *= evecp[i] / qvecp[i];
                    a2 += b2;
                    if (.563 < a2 || MAX(b2, b1) < a2 * .01) {
                        break;
                    }
                }
                a2 *= 1.05;
            }

            if (a2 < .563) {
                s = dn[2] * (1. - sqrt(a2)) / (a2 + 1.);
            } else {
                s = dmin[0] * .25;
            }

        } else {  // dmin == 0 or the smallest diag located at the forth or farther place

/*           Case 6, no information to guide us. */

            s = dmin[0] * .25;
        }

    } else if (n0init == n0 + 1) {

/*        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN. */

        if (dmin[1] == dn[1] && dmin[2] == dn[2]) {

            if (evecp[n0 - 2] > qvecp[n0 - 2]) {
                return 0;
            }
            s = dmin[1] * .333;
            b1 = evecp[n0 - 2] / qvecp[n0 - 2];
            b2 = b1;
            if (b2 != 0.) {
                for (i = n0 - 3; i >= i0; i--) {
                    a2 = b1;
                    if (evecp[i] > qvecp[i]) {
                        return 0;
                    }
                    b1 *= evecp[i] / qvecp[i];
                    b2 += b1;
                    if (MAX(b1, a2) * 100. < b2) {
                        break;
                    }
                }
            }
            a2 = dmin[1] / (b2 * 1.05 + 1.);
            b2 = sqrt(b2 * 1.05);
            gap2 = dmin[2] * .5 - a2;
            if (gap2 > 0. && gap2 > b2 * a2) {
                s = MAX(s, a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2));
            } else {
                s = MAX(s, a2 * (1. - b2 * 1.01));
            }
        } else if (dmin[1] == dn[1]) {
            s = dmin[1] * .5;
        } else {
            s = dmin[1] * .25;
        }

    } else if (n0init == n0 + 2) {

/*        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN. */

        if (dmin[2] == dn[2] && evecp[n0 - 2] * 2. < qvecp[n0 - 2]) {
            if (evecp[n0 - 2] > qvecp[n0 - 2]) {
                return 0;
            }
            b1 = evecp[n0 - 2] / qvecp[n0 - 2];
            b2 = b1;
            if (b2 != 0.) {
                for (i = n0 - 2; i > i0; i--) {
                    if (evecp[i - 1] > qvecp[i - 1]) {
                        return 0;
                    }
                    b1 *= evecp[i - 1] / qvecp[i - 1];
                    b2 += b1;
                    if (b1 * 100. < b2) {
                        break;
                    }
                }
            }
            s = dmin[2] * .333;
            a2 = dmin[2] / (b2 * 1.05 + 1.);
            b2 = sqrt(b2 * 1.05);
            gap2 = qvecp[n0 - 2] + evecp[n0 - 3] - sqrt(qvecp[n0 - 3] * evecp[n0 - 3]) - a2;
            if (gap2 > 0. && gap2 > b2 * a2) {
                s = MAX(s, a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2));
            } else {
                s = MAX(s, a2 * (1. - b2 * 1.01));
            }
        } else {
            s = dmin[2] * .25;
        }

    } else if (n0init > n0 + 2) {

/*        Case 12, more than two eigenvalues deflated. No information. */

        s = 0.;
    }

    *tau = s;
    return 0;
}

static void _dlasq5(int i0, int n0,
                    double *qvecp, double *qvec1p, double *evecp, double *evec1p,
                    double tau, double tol, double *dmin, double *dn)
{
    double diag = qvecp[i0] - tau;
    double diag_min = diag;
    double temp;
    int j;

    for (j = i0; j < n0 - 3; j++) {
        qvec1p[j] = diag + evecp[j];
        temp = qvecp[j + 1] / qvec1p[j];
        diag = diag * temp - tau;
        if (diag < tol) {
            diag = 0.;
        }
        diag_min = MIN(diag_min, diag);
        evec1p[j] = evecp[j] * temp;
    }
    dn[2] = diag;

    j = n0 - 3;
    qvec1p[j] = diag + evecp[j];
    temp = qvecp[j + 1] / qvec1p[j];
    evec1p[j] = evecp[j] * temp;
    diag = diag * temp - tau;
    dn[1] = diag;

    j = n0 - 2;
    qvec1p[j] = diag + evecp[j];
    temp = qvecp[j + 1] / qvec1p[j];
    evec1p[j] = evecp[j] * temp;
    diag = diag * temp - tau;
    dn[0] = diag;

    qvec1p[n0 - 1] = diag;

    dmin[2] = diag_min;
    dmin[1] = MIN(dmin[2], dn[1]);
    dmin[0] = MIN(dmin[1], dn[0]);
}


static int _dlasq2(int n, double *work, double *diag, double *diag_off)
{
    int i, j, itry, iwhilb, iter;
    int i0, n0, n1, n2, n0init, nbig;
    double emax, qmin, temp, diag_sum, tol, tol2, s, t;
    double dmin[3] = {0., 0., 0.};
    double dn[3] = {0., 0., 0.};
    double sigma, tau;
    double *qvec, *qvec1, *evec, *evec1;
    double *qvecp, *qvec1p, *evecp, *evec1p, *swap;
    qvec  = work;
    qvec1 = work + n * 1;
    evec  = work + n * 2;
    evec1 = work + n * 3;

    assert(n > 2);

    for (i = 0, j = n - 1; i < n - 1; ++i, --j) {
        temp = fabs(diag[i]);

        qvec1[i] = 0.;
        evec1[i] = 0.;
        evec[j - 1] = diag_off[i] * diag_off[i] * temp;
        qvec[j] = temp;
    }
    //qvec[n - 1] = fabs(diag[n-1]);
    qvec[0] = fabs(diag[n-1]);
    qvec1[n - 1] = 0.;
    evec [n - 1] = 0.;
    evec1[n - 1] = 0.;

/*     Reverse the qd-array, if warranted. */

    if (qvec[0] < qvec[n - 1] * 1.5) {
        for (i = 0, j = n - 1; i < n / 2; i++, --j) {
            temp = qvec[i];
            qvec[i] = qvec[j];
            qvec[j] = temp;

            temp = evec[i];
            evec[i] = evec[j - 1];
            evec[j - 1] = temp;
        }
    }

/*        dqd maps Z to ZZ plus Li's test. */
    diag_sum = qvec[0];
    for (i = 0; i < n - 1; i++) {
            temp = diag_sum + evec[i];
            qvec1[i] = temp;
            evec1[i] = qvec[i + 1] * (evec[i] / temp);
            diag_sum = qvec[i + 1] * (diag_sum / temp);
    }
    qvec1[n-1] = diag_sum;
    diag_sum = qvec1[0];
    for (i = 0; i < n - 1; i++) {
            temp = diag_sum + evec1[i];
            qvec[i] = temp;
            evec[i] = qvec1[i + 1] * (evec1[i] / temp);
            diag_sum = qvec1[i + 1] * (diag_sum / temp);
    }
    qvec[n-1] = diag_sum;

    n0 = n;
    tau = 0.;

    for (itry = 0; n0 > 0; ++itry) {
        if (itry >= n) {
            return 3;
        }

/*        E(N0) holds the value of SIGMA when submatrix in I0:N0 */
/*        splits from the rest of the array, but is negated. */

        sigma = -evec[n0 - 1];
        if (sigma < 0.) {
            return 1;
        }

/*        Find last unreduced submatrix's top index I0, find QMAX and */
/*        EMIN. Find Gershgorin-type bound if Q's much greater than E's. */

        emax = 0.;
        qmin = qvec[n0 - 1];
        for (i = n0 - 1; i > 0; i--) {
            if (evec[i - 1] <= 0.) {
                break;
            }
            if (qmin >= emax * 4.) {
                qmin = MIN(qmin, qvec[i]);
                emax = MAX(emax, evec[i - 1]);
            }
        }

        i0 = i;
        qvecp  = qvec;
        qvec1p = qvec1;
        evecp  = evec;
        evec1p = evec1;

        if (qmin < 0 || emax < 0) {
            fprintf(stderr, "dlasq2: qmin < 0 or emax < 0");
            return 1;
        }

/*        Put -(initial shift) into DMIN. */

        dmin[0] = -MAX(0., qmin - 2. * sqrt(qmin * emax));

/*        Now I0:N0 is unreduced. */

        nbig = (n0 - i0) * 10;
        for (iwhilb = 0; iwhilb < nbig; ++iwhilb) {
/*           While submatrix unfinished take a good dqds step. */

            n0init = n0;
            tol = DBL_EPSILON * 100.;
            tol2 = tol * tol;

            while (n0 > i0) {
                n1 = n0 - 1;
                n2 = n0 - 2;
                if (n1 == i0 || evecp[n2] < tol2 * (sigma + qvecp[n1]) || evec1p[n2] < tol2 * qvecp[n2]) {
                    qvec[n1] = qvecp[n1] + sigma;
                    --n0;
                    continue;
                }
                if (n2 == i0 || evecp[n1 - 2] < tol2 * sigma || evec1p[n1 - 2] < tol2 * qvecp[n1 - 2]) {
                    if (qvecp[n1] > qvecp[n2]) {
                        s = qvecp[n1];
                        qvecp[n1] = qvecp[n2];
                        qvecp[n2] = s;
                    }
                    t = (qvecp[n2] - qvecp[n1] + evecp[n2]) * .5;
                    if (evecp[n2] > qvecp[n1] * tol2 && t != 0.) {
                        s = qvecp[n1] * (evecp[n2] / t);
                        s = qvecp[n1] * (evecp[n2] / (t + sqrt(t * (t + s))));
                        t = qvecp[n2] + (s + evecp[n2]);
                        qvecp[n1] *= qvecp[n2] / t;
                        qvecp[n2] = t;
                    }
                    qvec[n2] = qvecp[n2] + sigma;
                    qvec[n1] = qvecp[n1] + sigma;
                    n0 += -2;
                    continue;
                }
                break;
            }
            if (n0 <= i0) {
                break;
            }

            if (dmin[0] <= 0.) {
                tau = -dmin[0];
            } else {
                _dlasq4(i0, n0, n0init, qvecp, qvec1p, evecp, evec1p, dmin, dn, &tau);
            }

/*     Call dqds until DMIN > 0. */

            tol = DBL_EPSILON * sigma;
            for (iter = 0; iter < 3; iter++) {

                _dlasq5(i0, n0, qvecp, qvec1p, evecp, evec1p, tau, tol, dmin, dn);

                if (dmin[0] >= 0.) {
                    break;
                } else if (dmin[1] > 0.) {
                    tau += dmin[0];
                } else {
                    tau *= .25;
                }
            }
            if (dmin[0] < 0) {
                tau = 0.;
                _dlasq5(i0, n0, qvecp, qvec1p, evecp, evec1p, tau, tol, dmin, dn);
            }

            sigma += tau;

            // swap vec, vec1
            swap = qvecp;
            qvecp = qvec1p;
            qvec1p = swap;
            swap = evecp;
            evecp = evec1p;
            evec1p = swap;
        }

        if (iwhilb == nbig) {
            // TODO: not converged. raise error
            fprintf(stderr, "dlasq2: Maximum number of iterations exceeded");
            return 2;
        }
    }
    return 0;
}

static int _compute_eigenvalues(int n, double *diag, double *diag_off1,
                                double *w, double *werr, double *wgap, double *work)
{
    double gl, gu;
    double eabs, eold, tmp, tmp1, dmax;
    double eps, rtol, rtl;
    double sigma, tau;
    double sgndef, spectral_diameter;
    double isleft, isrght, dpivot;
    int idum, ip, i;
    int lcnt, rcnt, norep;
    int iinfo;

    if (n <= 0) {
        return 0;
    }

    if (n == 1) {
        w[0] = diag[0];
        werr[0] = 0.;
        wgap[0] = 0.;
/*        store the shift for the initial RRR, which is zero in this case */
        diag_off1[0] = 0.;
        return 0;
    }

    eps = 2 * DBL_EPSILON;
    rtl = 2.1e-8;  // sqrt(eps)
    rtol = 16. * DBL_EPSILON;

/*     Init WERR, WGAP. Compute Gerschgorin intervals and spectral diameter. */
/*        Set interval [VL,VU] that contains all eigenvalues */
    gl = diag[0];
    gu = diag[0];
    eold = 0.;
    diag_off1[n-1] = 0.;
    for (i = 0; i < n; ++i) {
        werr[i] = 0.;
        wgap[i] = 0.;
        eabs = fabs(diag_off1[i]);
        tmp1 = eabs + eold;
        gl = MIN(gl, diag[i] - tmp1);
        gu = MAX(gu, diag[i] + tmp1);
        eold = eabs;
    }
/*     Compute spectral diameter. The Gerschgorin bounds give an */
/*     estimate that is wrong by at most a factor of SQRT(2) */
    spectral_diameter = gu - gl;
/* will hold the shift for the initial RRR, for now set it =0 */
    diag_off1[n-1] = 0.;

    if (n == 1) {
        w[0] = diag[0];
        werr[0] = 0.;
        wgap[0] = 0.;
        return 0;
    }

    double *e2 = work;
    work += n;
    for (i = 0; i < n - 1; ++i) {
        e2[i] = diag_off1[i] * diag_off1[i];
    }

/*           Case of DQDS */
/*           Find approximations to the extremal eigenvalues of the block */
    iinfo = _dlarrk(n, 1, gl, gu, diag, e2, rtl, &tmp, &tmp1);
    if (iinfo != 0) {
        return -1;
    }
    isleft = MAX(gl, tmp - tmp1 - eps * 100. * fabs(tmp - tmp1));
    iinfo = _dlarrk(n, n, gl, gu, diag, e2, rtl, &tmp, &tmp1);
    if (iinfo != 0) {
        return -1;
    }
    isrght = MIN(gu, tmp + tmp1 + eps * 100. * fabs(tmp + tmp1));
/*           Improve the estimate of the spectral diameter */
    spectral_diameter = isrght - isleft;
/*        Decide whether the base representation for the current block */
/*        L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I */
/*        should be on the left or the right end of the current block. */
/*        The strategy is to shift to the end which is "more populated" */
/*        Furthermore, decide whether to use DQDS for the computation of */
/*        the eigenvalue approximations at the end of DLARRE or bisection. */
/*        dqds is chosen if all eigenvalues are desired or the number of */
/*        eigenvalues to be computed is large compared to the blocksize. */
/*           If all the eigenvalues have to be computed, we use dqd */
/*           INDL is the local index of the first eigenvalue to compute */

    if (n > 1) {
/*        Compute the negcount at the 1/4 and 3/4 points */
        _dlarrc(n, isleft + spectral_diameter * .25, isrght - spectral_diameter * .25,
                diag, e2, &lcnt, &rcnt);
    }

    if (n == 1) {
        sigma = gl;
        sgndef = 1.;
    } else if (lcnt - 1 >= n - rcnt) {
        sigma = MAX(isleft, gl);
        sgndef = 1.;
    } else {
        sigma = MIN(isrght, gu);
        sgndef = -1.;
    }
/*        An initial SIGMA has been chosen that will be used for computing */
/*        T - SIGMA I = L D L^T */
/*        Define the increment TAU of the shift in case the initial shift */
/*        needs to be refined to obtain a factorization with not too much */
/*        element growth. */
/*           The initial SIGMA was to the outer end of the spectrum */
/*           the matrix is definite and we need not retreat. */
    tau = MAX(spectral_diameter * n, 2. * fabs(sigma)) * eps;

    for (idum = 0; idum < MAXTRY; ++idum) {
/*           Compute L D L^T factorization of tridiagonal matrix T - sigma I. */
/*           Store D in WORK(1:IN), L in WORK(IN+1:2*IN), and reciprocals of */
/*           pivots in WORK(2*IN+1:3*IN) */
        dpivot = diag[0] - sigma;
        work[0] = dpivot;
        dmax = fabs(work[0]);
        for (i = 0; i < n - 1; ++i) {
            work[n * 2 + i] = 1. / work[i];
            tmp = diag_off1[i] * work[n * 2 + i];
            work[n + i] = tmp;
            dpivot = diag[i+1] - sigma - tmp * diag_off1[i];
            work[i+1] = dpivot;
            dmax = MAX(dmax, fabs(dpivot));
        }
/*           check for element growth */
        norep = dmax > spectral_diameter * 64.;
        if (!norep){
/*              Ensure the definiteness of the representation */
/*              All entries of D (of L D L^T) must have the same sign */
            for (i = 0; i < n; ++i) {
                tmp = sgndef * work[i];
                if (tmp < 0.) {
                    norep = 1;
                    break;
                }
            }
        }
        if (norep) {
/*              Note that in the case of IRANGE=ALLRNG, we use the Gerschgorin */
/*              shift which makes the matrix definite. So we should end up */
/*              here really only in the case of IRANGE = VALRNG or INDRNG. */
            if (idum == MAXTRY - 1) {
                if (sgndef == 1.) {
/*                    The fudged Gerschgorin shift should succeed */
                    sigma = gl - spectral_diameter * 2. * eps * n;
                } else {
                    sigma = gu + spectral_diameter * 2. * eps * n;
                }
            } else if (idum == MAXTRY) {
/*        if the program reaches this point, no base representation could be */
/*        found in MAXTRY iterations. */
                return -2;
            } else {
                sigma -= sgndef * tau;
                tau *= 2.;
            }
        } else {
            break;
        }
    }

/*        At this point, we have found an initial base representation */
/*        T - SIGMA I = L D L^T with not too much element growth. */
/*        Store the shift. */
    diag_off1[n-1] = sigma;
/*        Store D and L. */
    for (ip = 0; ip < n; ip++) {
            diag[ip] = work[ip];
    }
    for (ip = 0; ip < n - 1; ip++) {
            diag_off1[ip] = work[n + ip];
    }

/*           Call dqds to get all eigs (and then possibly delete unwanted */
/*           eigenvalues). */
/*           Note that dqds finds the eigenvalues of the L D L^T representation */
/*           of T to high relative accuracy. High relative accuracy */
/*           might be lost when the shift of the RRR is subtracted to obtain */
/*           the eigenvalues of T. However, T is not guaranteed to define its */
/*           eigenvalues to high relative accuracy anyway. */
/*           Set RTOL to the order of the tolerance used in DLASQ2 */
/*           This is an ESTIMATED error, the worst case bound is 4*N*EPS */
/*           which is usually too large and requires unnecessary work to be */
/*           done by bisection when computing the eigenvectors */
    iinfo = _dlasq2(n, work, diag, diag_off1);
    if (iinfo != 0) {
        return -5;
    } else {
        for (i = 0; i < n; ++i) {
            if (work[i] < 0.) {
                fprintf(stderr, "dlarre: negative eigenvalues\n");
                return -6;
            }
        }
    }
    if (sgndef > 0.) {
        for (i = 0; i < n; ++i) {
            w[i] = work[n - 1 - i];
        }
    } else {
        for (i = 0; i < n; ++i) {
            w[i] = -work[i];
        }
    }

    for (i = 0; i < n; ++i) {
        werr[i] = rtol * fabs(w[i]);
    }

/*              compute the right gap between the intervals */
    for (i = 0; i < n - 1; ++i) {
        wgap[i] = MAX(0., w[i + 1] - werr[i + 1] - (w[i] + werr[i]));
    }
    wgap[-1] = MAX(0., w[0] - werr[0] - gl);
    wgap[n-1] = MAX(0., gu - sigma - (w[n-1] + werr[n-1]));

    return 0;
}

static int _dlarrf(int n, double *diag, double *diag_off1, double *ld, int clstrt,
                   double *w, double *wgap, double *werr, double clgapl,
                   double *sigma, double *dplus, double *lplus)
{
    int i, ktry;
    double s, tmp, max1, growthbound, lsigma;

/*     Use a small fudge to make sure that we really shift to the outside */
    lsigma = w[clstrt] - werr[clstrt];
    lsigma -= fabs(lsigma) * 4. * DBL_EPSILON;

    growthbound = diag[0] * 8.;

    for (ktry = 0; ktry < 2; ++ktry) {
/*     Compute the element growth when shifting to both ends of the cluster */
/*     accept the shift if there is no element growth at one of the two ends */
/*     Left end */
        s = -lsigma;
        dplus[0] = diag[0] + s;
        max1 = fabs(dplus[0]);
        for (i = 0; i < n - 1; ++i) {
            tmp = ld[i] / dplus[i];
            lplus[i] = tmp;
            s = s * tmp * diag_off1[i] - lsigma;
            dplus[i + 1] = diag[i + 1] + s;
            max1 = MAX(max1, fabs(dplus[i + 1]));
        }
        *sigma = lsigma;
        if (max1 <= growthbound) {
            return 0;
        }

/*        If we are here, shifts failed the RRR test. */
/*        Back off to the outside */
        lsigma = lsigma - MIN(clgapl * .25, wgap[clstrt]);
    }

    if (max1 > FAIL_THRESHOLD) {
        fprintf(stderr, "dlarrf max1 = %g", max1);
        return 1;
    }
    return 0;
}

static int _dlaneg(int n, double *diag, double *lld, double sigma, int twist_index)
{
    int j;
    int negcnt;
    double p, t;
    double dplus, dminus;

    negcnt = 0;
/*     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T */
    t = -sigma;
    for (j = 0; j < twist_index - 1; ++j) {
        dplus = diag[j] + t;
        if (dplus < 0.) {
            ++negcnt;
        }
        t = t / dplus * lld[j] - sigma;
    }

/*     II) lower part: L D L^T - SIGMA I = U- D- U-^T */
    p = diag[n-1] - sigma;
    for (j = n - 2; j >= twist_index - 1; --j) {
        dminus = lld[j] + p;
        if (dminus < 0.) {
            ++negcnt;
        }
        p = p / dminus * diag[j] - sigma;
    }

/*     III) Twist index */
/*       T was shifted by SIGMA initially. */
    if (t + sigma + p < 0.) {
        ++negcnt;
    }
    return negcnt;
}

/* DLARRB provides limited bisection to locate eigenvalues for more accuracy. */
static int _dlarrb(int n, double *diag, double *lld,
                   int ifirst, int ilast, double rtol1, double rtol2,
                   double *w, double *wgap, double *werr, int twist_index)
{
    int i, iter, negcnt;
    double mid, back, left, right, width;
    double cvrgd = MAX(rtol1*wgap[ifirst],
                       rtol2 * MAX(fabs(w[ifirst]), fabs(w[ilast-1])));

    for (i = ifirst; i < ilast; ++i) {
        if (werr[i] < cvrgd) {
            continue;
        }

/*        Make sure that [LEFT,RIGHT] contains the desired eigenvalue */
/*        Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT */
        left = w[i];
        back = werr[i];
        for (negcnt = ilast; negcnt > i;) {
            left -= back;
            back *= 2.;
            negcnt = _dlaneg(n, diag, lld, left, twist_index);
        }

/*        Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT */
        right = w[i];
        back = werr[i];
        for (negcnt = ifirst; negcnt <= i;) {
            right += back;
            back *= 2.;
            negcnt = _dlaneg(n, diag, lld, right, twist_index);
        }

        for (iter = 0; iter < ITMAX; ++iter) {
/*        Perform one bisection step */
            mid = (left + right) * .5;
            width = right - mid;
            if (width < cvrgd) {
                break;
            }

            negcnt = _dlaneg(n, diag, lld, mid, twist_index);
            if (negcnt <= i) {
                left = mid;
            } else {
                right = mid;
            }
        }
        w[i] = mid;
        werr[i] = width;
    }

    for (i = ifirst; i < ilast-1; ++i) {
        wgap[i] = MAX(0., w[i+1] - werr[i+1] - w[i] - werr[i]);
    }
    return 0;
}

static void _dlar1v(int n, double lambda, double *diag, double *diag_off1,
                    double *ld, double *lld, double gaptol, double *vec, int *negcnt,
                    int *twist_index, double *resid, double *rqcorr, double *work)
{
    int i, r1, r2;
    int neg1, neg2;
    double s, tmp;
    double nrminv, mingma, dplus, dminus;
    double *lplus = work;
    double *uminus = work + n;
    double *work_p = work + n * 2;

    if (*twist_index == -1) {
        r1 = 0;
        r2 = n;
        *twist_index = 0;
    } else {
        r1 = *twist_index;
        r2 = *twist_index + 1;
    }

/*     L D L**T - sigma I. */
/*     Compute the progressive transform (using the differential form) */
/*     until the index R1 */

    neg2 = 0;
    s = diag[n - 1] - lambda;
    work_p[n - 1] = s;
    for (i = n - 2; i >= r1; --i) {
        dminus = lld[i] + s;
        if (dminus < 0.) {
            ++neg2;
        }
        tmp = diag[i] / dminus;
        uminus[i] = diag_off1[i] * tmp;
        s = s * tmp - lambda;
        work_p[i] = s;
    }

/*     Compute the stationary transform (using the differential form) */
/*     until the index R2. */

    neg1 = 0;
    s = -lambda;
    for (i = 0; i < r1; ++i) {
        dplus = diag[i] + s;
        if (dplus < 0.) {
            ++neg1;
        }
        lplus[i] = ld[i] / dplus;
        s = s * lplus[i] * diag_off1[i] - lambda;
    }
    mingma = s + lambda + work_p[r1];
    if (mingma < 0.) {
        ++neg1;
    }

    *negcnt = neg1 + neg2;

/*     Find the index (from R1 to R2) of the largest (in magnitude) */
/*     diagonal element of the inverse */

    for (i = r1; i < r2 - 1; ++i) {
        dplus = diag[i] + s;
        lplus[i] = ld[i] / dplus;
        tmp = s * lplus[i] * diag_off1[i];
        s = tmp - lambda;

        tmp = tmp + work_p[i+1];
        if (fabs(tmp) <= fabs(mingma)) {
            mingma = tmp;
            *twist_index = i + 1;
        }
    }

/*     Compute the FP vector: solve N^T v = e_r */

    vec[*twist_index] = 1.;
    double ztz = 1.;

    for (i = *twist_index - 1; i >= 0; --i) {
        tmp = -(lplus[i] * vec[i + 1]);
        ztz += tmp * tmp;
        vec[i] = tmp;
    }

    for (i = *twist_index; i < n-1; ++i) {
        tmp = -(uminus[i] * vec[i]);
        ztz += tmp * tmp;
        vec[i+1] = tmp;
    }

    tmp = 1. / ztz;
    nrminv = sqrt(tmp);
    for (i = 0; i < n; i++) {
        vec[i] *= nrminv;
    }
    *resid = fabs(mingma) * nrminv;
    *rqcorr = mingma * tmp;
}

/* the eigenvectors of the tridiagonal matrix T = L D LT given L, D and the eigenv alues of L D LT. */
/* TODO: Compute only the fist element of each eigenvectors */
static int _compute_eigenvectors(int n, double *diag, double *diag_off1,
                                 double *w, double *werr, double *wgap,
                                 double *vec, double *work, int *iwork)
{
    int i, j, k, icls, iter, idone, ndepth;
    int ncluster, ncluster1, negcnt;
    int oldfst, oldlst;
    int newfst, newlst;
    int needbs;
    int iinfo;
    double fudge, eps, rqtol, tol, tmp;
    double left, right, gap, bstw, savgap, gaptol;
    double sigma, tau, resid, lambda, bstres;
    double rqcorr, resid_tol, rqcorr_tol;
    double *buf_w, *buf_ld, *buf_lld, *buf_wrk;
    int *twist_indices;
    int *swap, *old_cluster_range, *new_cluster_range;

    if (n <= 0) {
        return 0;
    }

/*     The first N entries of WORK are reserved for the eigenvalues */
    buf_w = work;
    buf_ld = work + n;
    buf_lld = work + n * 2;
    buf_wrk = work + n * 3;
    for (i = 0; i < n * 6; ++i) {
        work[i] = 0.;
    }
    twist_indices = iwork;
    for (i = 0; i < n; ++i) {
        twist_indices[i] = 0;
    }
    old_cluster_range = iwork + n;
    new_cluster_range = iwork + n * 3;

    eps = DBL_EPSILON;
    rqtol = DBL_EPSILON * 2.;
    tol = DBL_EPSILON * 8;

    sigma = diag_off1[n-1];
    if (1 == n) {
        vec[0] = 1.;
        w[0] += sigma;
        return 0;
    }
/*        The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND) */
/*        Note that these can be approximations, in this case, the corresp. */
/*        entries of WERR give the size of the uncertainty interval. */
/*        The eigenvalue approximations will be refined when necessary as */
/*        high relative accuracy is required for the computation of the */
/*        corresponding eigenvectors. */
    for (i = 0; i < n; i++) {
        buf_w[i] = w[i];
        w[i] += sigma;
    }
    ncluster = 1;
    old_cluster_range[0] = 0;
    old_cluster_range[1] = n;
    idone = 0;
    for (ndepth = 0; ndepth < n; ndepth++) {
        if (idone == n) {
            break;
        }

        ncluster1 = ncluster;
        ncluster = 0;
/*           Process the clusters on the current level */
        for (icls = 0; icls < ncluster1; ++icls) {
/*              OLDFST, OLDLST = first, last index of current cluster. */
            oldfst = old_cluster_range[icls * 2    ];
            oldlst = old_cluster_range[icls * 2 + 1];
            if (ndepth > 0) {
/*                 Retrieve relatively robust representation (RRR) of cluster */
/*                 that has been computed at the previous level */
/*                 The RRR is stored in Z and overwritten once the eigenvectors */
/*                 have been computed or when the cluster is refined */
/*                    Get representation from location of the leftmost evalue */
/*                    of the cluster */
                for (i = 0; i < n; i++) {
                    diag[i] = vec[oldfst * n + i];
                }
                for (i = 0; i < n - 1; i++) {
                    diag_off1[i] = vec[(oldfst+1) * n + i];
                }
                sigma = vec[(oldfst+2) * n - 1];
            }
/*              Compute DL and DLL of current RRR */
            for (j = 0; j < n - 1; ++j) {
                tmp = diag[j] * diag_off1[j];
                buf_ld [j] = tmp;
                buf_lld[j] = tmp * diag_off1[j];
            }
/*                 perform limited bisection (if necessary) to get approximate */
/*                 eigenvalues to the precision needed. */
/*
            if (ndepth > 0) {
                _dlarrb(n, diag, buf_lld, oldfst, oldlst, RTOL1, RTOL2, buf_w, wgap, werr, n);
                if (oldfst > 0) {
                    wgap[oldfst-1] = MAX(wgap[oldfst-1], w[oldfst] - werr[oldfst] - w[oldfst-1] - werr[oldfst-1]);
                }
                if (oldlst < n) {
                    wgap[oldlst-1] = MAX(wgap[oldlst-1], w[oldlst] - werr[oldlst] - w[oldlst-1] - werr[oldlst-1]);
                }
                for (j = oldfst; j < oldlst; ++j) {
                    w[j] = buf_w[j] + sigma;
                }
            }
*/

/*              Process the current node. */
            newfst = oldfst;
            for (newlst = oldfst+1; newlst <= oldlst; ++newlst) {
                if (newlst < oldlst && wgap[newlst-1] < MINGAP * fabs(buf_w[newlst-1])) {
                    continue;
                }

                if (newlst - newfst > 1) {
/*                    Compute left- and rightmost eigenvalue of child */
/*                    to high precision in order to shift as close */
/*                    as possible and obtain as large relative gaps */
/*                    as possible */
                    _dlarrb(n, diag, buf_lld, newfst, newfst+1, rqtol, rqtol, buf_w, wgap, werr, n);
                    _dlarrb(n, diag, buf_lld, newlst-1, newlst, rqtol, rqtol, buf_w, wgap, werr, n);

/*                    Compute RRR of child cluster. */
/*                    Note that the new RRR is stored in Z */
                    iinfo = _dlarrf(n, diag, diag_off1, buf_ld, newfst, buf_w,
                                    wgap, werr, wgap[newfst-1],
                                    &tau, vec+newfst*n, vec+(newfst+1)*n);
                    if (iinfo != 0) {
                        return -2;
                    }
/*                       a new RRR for the cluster was found by DLARRF */
/*                       update shift and store it */
                    vec[(newfst+2) * n - 1] = sigma + tau;

                    for (k = newfst; k < newlst; ++k) {
                        fudge = eps * 3. * fabs(buf_w[k]);
                        buf_w[k] -= tau;
                        fudge += eps * 4. * fabs(buf_w[k]);
/*                          Fudge errors */
                        werr[k] += fudge;
/*                          Gaps are not fudged. Provided that WERR is small */
/*                          when eigenvalues are close, a zero gap indicates */
/*                          that a new representation is needed for resolving */
/*                          the cluster. A fudge could lead to a wrong decision */
/*                          of judging eigenvalues 'separated' which in */
/*                          reality are not. This could have a negative impact */
/*                          on the orthogonality of the computed eigenvectors. */
                    }
                    new_cluster_range[ncluster * 2    ] = newfst;
                    new_cluster_range[ncluster * 2 + 1] = newlst;
                    ++ncluster;

                } else {
/*                    Compute eigenvector of singleton */

                    lambda = buf_w[newfst];
                    left = lambda - werr[newfst];
                    right = lambda + werr[newfst];
                    gap = wgap[newfst];
                    if (newfst == 0 || newfst+1 == n) {
                        gaptol = 0;
                    } else {
                        gaptol = gap * eps;
                    }
                    savgap = gap;
                    resid_tol = tol * gap;
                    rqcorr_tol = rqtol * fabs(lambda);

/*                    We want to use the Rayleigh Quotient Correction */
/*                    as often as possible since it converges quadratically */
/*                    when we are close enough to the desired eigenvalue. */
/*                    However, the Rayleigh Quotient can have the wrong sign */
/*                    and lead us away from the desired eigenvalue. In this */
/*                    case, the best we can do is to use bisection. */

                    needbs = 0;
                    bstres = 1e307;
                    bstw = 0;
                    twist_indices[newfst] = -1;
                    for (iter = 0; iter < MAXRQITER; ++iter) {
                        _dlar1v(n, lambda, diag, diag_off1,
                                buf_ld, buf_lld, gaptol, vec+newfst*n,
                                &negcnt, twist_indices+newfst, &resid, &rqcorr, buf_wrk);

                        /* keep track of the best lambda */
                        if (resid < bstres) {
                            bstres = resid;
                            bstw = lambda;
                        }

/*                    Convergence test for Rayleigh-Quotient iteration */

                        if (resid < resid_tol || fabs(rqcorr) < rqcorr_tol) {
                            break;
                        }

/*                       We only use the RQCORR if it improves the */
/*                       the iterate reasonably. */
                        if (lambda + rqcorr > right || lambda + rqcorr < left) {
                            needbs = 1;
                            break;
                        }

/*                       We need to check that the RQCORR update doesn't */
/*                       move the eigenvalue away from the desired one and */
/*                       towards a neighbor. */
                        if (newfst < negcnt) {
/*                          The wanted eigenvalue lies to the left */
                            if (rqcorr > 0) {
                                needbs = 1;
                                break;
                            }
                            right = lambda;
                        } else {
/*                          The wanted eigenvalue lies to the right */
                            if (rqcorr < 0) {
                                needbs = 1;
                                break;
                            }
                            left = lambda;
                        }

                        buf_w[newfst] = (right + left) * .5;
                        lambda += rqcorr;
                        werr[newfst] = (right - left) * .5;

                        if (right - left < rqcorr_tol) {
                            if (bstres < resid) {
                                lambda = bstw;
                                _dlar1v(n, lambda, diag, diag_off1, buf_ld, buf_lld,
                                        gaptol, vec+newfst*n,
                                        &negcnt, twist_indices+newfst, &resid, &rqcorr, buf_wrk);
                            }
                            break;
                        }
                    }

                    if (needbs) {
                        _dlarrb(n, diag, buf_lld, newfst, newfst+1, 0., eps * 2.,
                                buf_w, wgap, werr, twist_indices[newfst] + 1);
                        lambda = buf_w[newfst];

/*                       Reset twist index from inaccurate LAMBDA to */
/*                       force computation of true MINGMA */
                        twist_indices[newfst] = -1;
                        _dlar1v(n, lambda, diag, diag_off1,
                                buf_ld, buf_lld, gaptol, vec+newfst*n,
                                &negcnt, twist_indices+newfst, &resid, &rqcorr, buf_wrk);
                    }

                    w[newfst] = lambda + sigma;

/*                    Recompute the gaps on the left and right */
/*                    But only allow them to become larger and not */
/*                    smaller (which can only happen through "bad" */
/*                    cancellation and doesn't reflect the theory */
/*                    where the initial gaps are underestimated due */
/*                    to WERR being too crude.) */
                    if (newfst > 0) {
                        wgap[newfst-1] = MAX(wgap[newfst-1], w[newfst] - werr[ newfst] - w[newfst-1] - werr[newfst-1]);
                    }
                    if (newfst < n-1) {
                        wgap[newfst] = MAX(savgap, w[newfst+1] - werr[newfst+1] - w[newfst] - werr[newfst]);
                    }

                    ++idone;
                }

                newfst = newlst;
            }
        }

        swap = old_cluster_range;
        old_cluster_range = new_cluster_range;
        new_cluster_range = swap;
    }
    if (idone < n) {
        return -2;
    }
    return 0;
}

static int _dlaev2(double *eig, double *vec, double *diag, double *diag_off1)
{
    double a = diag[0];
    double b = diag_off1[0];
    double c = diag[1];
    double df, cs, ct, tb, sm, tn, rt, tmp;
    double rt1, rt2, cs1, sn1;
    int sgn1, sgn2;

    sm = a + c;
    df = a - c;
    tb = b + b;

    rt = sqrt(tb * tb + df * df);

    if (sm > 0.) {
        rt1 = (sm + rt) * .5;
        sgn1 = 1;
        rt2 = (a * c - b * b) / rt1;
    } else if (sm < 0.) {
        rt1 = (sm - rt) * .5;
        sgn1 = -1;
        rt2 = (a * c - b * b) / rt1;
    } else {
        rt1 = rt * .5;
        rt2 = rt * -.5;
        sgn1 = 1;
    }

/*     Compute the eigenvector */

    if (df >= 0.) {
        cs = df + rt;
        sgn2 = 1;
    } else {
        cs = df - rt;
        sgn2 = -1;
    }

    if (fabs(cs) > fabs(tb)) {
        ct = -tb / cs;
        sn1 = 1. / sqrt(ct * ct + 1.);
        cs1 = ct * sn1;
    } else {
        if (b == 0.) {
            cs1 = 1.;
            sn1 = 0.;
        } else {
            tn = -cs / tb;
            cs1 = 1. / sqrt(tn * tn + 1.);
            sn1 = tn * cs1;
        }
    }

    if (sgn1 == sgn2) {
        tmp = cs1;
        cs1 = -sn1;
        sn1 = tmp;
    }

    eig[0] = rt2;
    eig[1] = rt1;
    vec[0] = -sn1;
    vec[1] = cs1;
    vec[2] = cs1;
    vec[3] = sn1;
    return 0;
}

int _CINTdiagonalize(int n, double *diag, double *diag_off1, double *eig, double *vec)
{
    if (n == 0) {
        return 0;
    } else if (n == 1) {
        eig[0] = diag[0];
        vec[0] = 1.;
        return 0;
    } else if (n == 2) {
        return _dlaev2(eig, vec, diag, diag_off1);
    }

    int iwork[MXRYSROOTS * 5];
    double work[MXRYSROOTS * 9 + 1];
    double *buf_err = work + n;
    double *buf_gp  = work + n * 2 + 1;
    double *buf_wrk = work + n * 3 + 1;
    int info;

    info = _compute_eigenvalues(n, diag, diag_off1, eig, buf_err, buf_gp, buf_wrk);
    if (info == 0) {
            info = _compute_eigenvectors(n, diag, diag_off1, eig, buf_err, buf_gp,
                                         vec, buf_wrk, iwork);
    }
    return info;
}

#endif  // LAPACK_FOUND
