/*
 * File: misc.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic functions
 */

#include <math.h>

void CINTdcmplx_re(const unsigned int n, double *z, const double *re)
{
        unsigned int i;
        for (i = 0; i < n; i++) {
                z[0] = re[i];
                z[1] = 0;
                z += 2;
        }
}

void CINTdcmplx_im(const unsigned int n, double *z, const double *im)
{
        unsigned int i;
        for (i = 0; i < n; i++) {
                z[0] = 0;
                z[1] = im[i];
                z += 2;
        }
}

void CINTdcmplx_pp(const unsigned int n, double *z,
                   const double *re, const double *im)
{
        unsigned int i;
        for (i = 0; i < n; i++) {
                z[0] = re[i];
                z[1] = im[i];
                z += 2;
        }
}
void CINTdcmplx_pn(const unsigned int n, double *z,
                   const double *re, const double *im)
{
        unsigned int i;
        for (i = 0; i < n; i++) {
                z[0] =  re[i];
                z[1] = -im[i];
                z += 2;
        }
}
void CINTdcmplx_np(const unsigned int n, double *z,
                   const double *re, const double *im)
{
        unsigned int i;
        for (i = 0; i < n; i++) {
                z[0] = -re[i];
                z[1] =  im[i];
                z += 2;
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

static int factorial(int n)
{
        int i, fact = 1;
        for (i = 1; i <= n; i++) {
                fact *= i;
        }
        return fact;
}
double CINTgto_norm(int n, double a)
{
        double nn = pow(2, (2*n+3)) * factorial(n+1) * pow((2*a), (n+1.5)) \
                / (factorial(2*n+2) * sqrt(M_PI));
        return sqrt(nn);
}
