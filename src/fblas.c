/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * blas-like functions
 */

#include <string.h>
#include <complex.h>
#include "fblas.h"

#define OF_CMPLX        2

void CINTdset0(const FINT n, double *x)
{
        memset(x, 0, sizeof(double) * n);
}


/*
 * v = a * x + y
 */
void CINTdaxpy2v(const FINT n, const double a,
                 const double *x, const double *y, double *v)
{
        //cblas_dcopy(n, y, 1, v, 1);
        //cblas_daxpy(n, a, x, 1, v, 1);
        FINT i;
        for (i = 0; i < n; i++) {
                v[i] = a * x[i] + y[i];
        }
}

void CINTzmat_transpose(double complex *a_t, const double complex *a,
                        const FINT m, const FINT n)
{
        FINT i, j;

        switch (m) {
        case 2:
                for (i = 0; i < n; i++) {
                        a_t[i  ] = a[2*i+0];
                        a_t[i+n] = a[2*i+1];
                }
                break;
        default:
                switch (n) {
                case 2: for (i = 0; i < m; i++) {
                                a_t[2*i+0] = a[i  ];
                                a_t[2*i+1] = a[i+m];
                        }
                        break;
                default:
                        for (j = 0; j < n; j++) {
                                for (i = 0; i < m; i++) {
                                        a_t[i*n+j] = a[j*m+i];
                                }
                        }
                }
        }
}

void CINTzmat_dagger(double complex *a_t, const double complex *a,
                     const FINT m, const FINT n)
{
        FINT i, j;

        for (j = 0; j < n; j++) {
                for (i = 0; i < m; i++) {
                        a_t[i*n+j] = conj(a[j*m+i]);
                }
        }
}

