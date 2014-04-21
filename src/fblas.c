/*
 * File: fblas.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * blas-like functions
 */

#include <string.h>

#define OF_CMPLX        2

void CINTdset0(const unsigned int n, double *x)
{
        memset(x, 0, sizeof(double) * n);
}


/*
 * v = a * x + y
 */
void CINTdaxpy2v(const unsigned int n, const double a,
                 const double *x, const double *y, double *v)
{
        //cblas_dcopy(n, y, 1, v, 1);
        //cblas_daxpy(n, a, x, 1, v, 1);
        unsigned int i;
        for (i = 0; i < n; i++) {
                v[i] = a * x[i] + y[i];
        }
}


/*
 * a_t[n,m] = transpose of matrix a[m,n]
 */
void CINTdmat_transpose(double *a_t, const double *a,
                        const unsigned int m, const unsigned int n)
{
        unsigned int i, j;
        double *pa_t;

        switch (m) {
                case 1:
                        switch (n) {
                                case 1:
                                        a_t[0] = a[0];
                                        break;
                                case 3:
                                        a_t[0] = a[0];
                                        a_t[1] = a[1];
                                        a_t[2] = a[2];
                                        break;
                                case 5:
                                        a_t[0] = a[0];
                                        a_t[1] = a[1];
                                        a_t[2] = a[2];
                                        a_t[3] = a[3];
                                        a_t[4] = a[4];
                                        break;
                                default:
                                        for (i = 0; i < n; i++)
                                                a_t[i] = a[i];
                        }
                        break;
                case 3:
                        switch (n) {
                                case 1:
                                        a_t[0] = a[0];
                                        a_t[1] = a[1];
                                        a_t[2] = a[2];
                                        break;
                                case 3:
                                        a_t[0] = a[0];
                                        a_t[3] = a[1];
                                        a_t[6] = a[2];
                                        a_t[1] = a[3];
                                        a_t[4] = a[4];
                                        a_t[7] = a[5];
                                        a_t[2] = a[6];
                                        a_t[5] = a[7];
                                        a_t[8] = a[8];
                                        break;
                                case 5:
                                        a_t[ 0] = a[ 0];
                                        a_t[ 5] = a[ 1];
                                        a_t[10] = a[ 2];
                                        a_t[ 1] = a[ 3];
                                        a_t[ 6] = a[ 4];
                                        a_t[11] = a[ 5];
                                        a_t[ 2] = a[ 6];
                                        a_t[ 7] = a[ 7];
                                        a_t[12] = a[ 8];
                                        a_t[ 3] = a[ 9];
                                        a_t[ 8] = a[10];
                                        a_t[13] = a[11];
                                        a_t[ 4] = a[12];
                                        a_t[ 9] = a[13];
                                        a_t[14] = a[14];
                                        break;
                                default:
                                        for (i = 0; i < n; i++) {
                                                a_t[i    ] = a[0];
                                                a_t[i+n  ] = a[1];
                                                a_t[i+n*2] = a[2];
                                                a += 3;
                                        }
                        }
                        break;
                case 5:
                        switch (n) {
                                case 1:
                                        a_t[0] = a[0];
                                        a_t[1] = a[1];
                                        a_t[2] = a[2];
                                        a_t[3] = a[3];
                                        a_t[4] = a[4];
                                        break;
                                case 3:
                                        a_t[ 0] = a[ 0];
                                        a_t[ 3] = a[ 1];
                                        a_t[ 6] = a[ 2];
                                        a_t[ 9] = a[ 3];
                                        a_t[12] = a[ 4];
                                        a_t[ 1] = a[ 5];
                                        a_t[ 4] = a[ 6];
                                        a_t[ 7] = a[ 7];
                                        a_t[10] = a[ 8];
                                        a_t[13] = a[ 9];
                                        a_t[ 2] = a[10];
                                        a_t[ 5] = a[11];
                                        a_t[ 8] = a[12];
                                        a_t[11] = a[13];
                                        a_t[14] = a[14];
                                        break;
                                case 5:
                                        a_t[ 0] = a[ 0];
                                        a_t[ 5] = a[ 1];
                                        a_t[10] = a[ 2];
                                        a_t[15] = a[ 3];
                                        a_t[20] = a[ 4];
                                        a_t[ 1] = a[ 5];
                                        a_t[ 6] = a[ 6];
                                        a_t[11] = a[ 7];
                                        a_t[16] = a[ 8];
                                        a_t[21] = a[ 9];
                                        a_t[ 2] = a[10];
                                        a_t[ 7] = a[11];
                                        a_t[12] = a[12];
                                        a_t[17] = a[13];
                                        a_t[22] = a[14];
                                        a_t[ 3] = a[15];
                                        a_t[ 8] = a[16];
                                        a_t[13] = a[17];
                                        a_t[18] = a[18];
                                        a_t[23] = a[19];
                                        a_t[ 4] = a[20];
                                        a_t[ 9] = a[21];
                                        a_t[14] = a[22];
                                        a_t[19] = a[23];
                                        a_t[24] = a[24];
                                        break;
                                default:
                                        for (i = 0; i < n; i++) {
                                                a_t[i    ] = a[0];
                                                a_t[i+n  ] = a[1];
                                                a_t[i+n*2] = a[2];
                                                a_t[i+n*3] = a[3];
                                                a_t[i+n*4] = a[4];
                                                a += 5;
                                        }
                        }
                        break;
                default:
                        switch (n) {
                                case 1:
                                        for (i = 0; i < m; i++)
                                                a_t[i] = a[i];
                                        break;
                                case 3:
                                        for (i = 0; i < m; i++) {
                                                a_t[0] = a[i    ];
                                                a_t[1] = a[i+m  ];
                                                a_t[2] = a[i+m*2];
                                                a_t += 3;
                                        }
                                        break;
                                case 5:
                                        for (i = 0; i < m; i++) {
                                                a_t[0] = a[i    ];
                                                a_t[1] = a[i+m  ];
                                                a_t[2] = a[i+m*2];
                                                a_t[3] = a[i+m*3];
                                                a_t[4] = a[i+m*4];
                                                a_t += 5;
                                        }
                                        break;
                                default:
                                        for (j = 0; j < n - 1; j += 2) {
                                                pa_t = a_t + j;
                                                for (i = 0; i < m; i++) {
                                                        pa_t[0  ] = a[i    ];
                                                        pa_t[1  ] = a[i+m  ];
                                                        pa_t += n;
                                                }
                                                a += m + m;
                                        }
                                        if (j < n) {
                                                pa_t = a_t + j;
                                                for (i = 0; i < m; i++) {
                                                        *pa_t = *a++;
                                                        pa_t += n;
                                                }
                                        }
                        }
        }
}

void CINTzmat_transpose(double *a_t, const double *a,
                        const unsigned int m, const unsigned int n)
{
        unsigned int i, j;
        double *pa_t;

        if (m == 1) {
                for (i = 0; i < n * OF_CMPLX; i++)
                        a_t[i] = a[i];
        } else if (n == 1) {
                for (i = 0; i < m * OF_CMPLX; i++)
                        a_t[i] = a[i];
        } else {
                for (j = 0; j < n * OF_CMPLX; j += OF_CMPLX) {
                        pa_t = a_t + j;
                        for (i = 0; i < m * OF_CMPLX; i += OF_CMPLX) {
                                pa_t[0] = a[i  ];
                                pa_t[1] = a[i+1];
                                pa_t += n * OF_CMPLX;
                        }
                        a += m * OF_CMPLX;
                }
        }
}

void CINTzmat_dagger(double *a_t, const double *a,
                     const unsigned int m, const unsigned int n)
{
        unsigned int i, j;
        double *pa_t;

        for (j = 0; j < n * OF_CMPLX; j += OF_CMPLX) {
                pa_t = a_t + j;
                for (i = 0; i < m * OF_CMPLX; i += OF_CMPLX) {
                        pa_t[0] = a[i  ];
                        pa_t[1] =-a[i+1];
                        pa_t += n * OF_CMPLX;
                }
                a += m * OF_CMPLX;
        }
}

