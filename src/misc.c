/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic functions
 */

#include <math.h>
#include <complex.h>
#include "config.h"

void CINTdcmplx_re(const FINT n, double complex *z, const double *re)
{
        FINT i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] + 0 * _Complex_I;
        }
}

void CINTdcmplx_im(const FINT n, double complex *z, const double *im)
{
        FINT i;
        for (i = 0; i < n; i++) {
                z[i] = 0 + im[i] * _Complex_I;
        }
}

void CINTdcmplx_pp(const FINT n, double complex *z,
                   const double *re, const double *im)
{
        FINT i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] + im[i] * _Complex_I;
        }
}
void CINTdcmplx_pn(const FINT n, double complex *z,
                   const double *re, const double *im)
{
        FINT i;
        for (i = 0; i < n; i++) {
                z[i] = re[i] - im[i] * _Complex_I;
        }
}
void CINTdcmplx_np(const FINT n, double complex *z,
                   const double *re, const double *im)
{
        FINT i;
        for (i = 0; i < n; i++) {
                z[i] = -re[i] + im[i] * _Complex_I;
        }
}
void CINTdcmplx_nn(const FINT n, double complex *z,
                   const double *re, const double *im)
{
        FINT i;
        for (i = 0; i < n; i++) {
                z[i] = -re[i] - im[i] * _Complex_I;
        }
}


double CINTsquare_dist(const double *r1, const double *r2)
{
        double r12[3];

        r12[0] = r1[0] - r2[0];
        r12[1] = r1[1] - r2[1];
        r12[2] = r1[2] - r2[2];

        return r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
}

static double _gaussian_int(FINT n, double alpha)
{
        double n1 = (n + 1) * .5;
        return exp(lgamma(n1)) / (2. * pow(alpha, n1));
}

/*
 * Normalized factor for GTO radial part g=r^l e^{-\alpha r^2}
 *
 * \frac{1}{\sqrt{\int g^2 r^2 dr}}
 *   = \sqrt{\frac{2^{2l+3} (l+1)! (2a)^{l+1.5}}{(2l+2)!\sqrt{\pi}}}
 *
 * Ref: H. B. Schlegel and M. J. Frisch, Int. J. Quant.  Chem., 54(1995), 83-87.
 */
double CINTgto_norm(FINT n, double a)
{
        //double nn = pow(2, (2*n+3)) * factorial(n+1) * pow((2*a), (n+1.5)) \
        //        / (factorial(2*n+2) * sqrt(M_PI));
        //return sqrt(nn);
        return 1. / sqrt(_gaussian_int(n*2+2, 2*a));

}
double CINTgto_norm_(FINT *n, double *a)
{
        return CINTgto_norm(*n, *a);
}
