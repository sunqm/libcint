/*
 * File: fblas.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * blas-like functions
 */

#include <string.h>
#include <complex.h>

#define OF_CMPLX        2

void CINTdset0(const int n, double *x)
{
        memset(x, 0, sizeof(double) * n);
}


/*
 * v = a * x + y
 */
void CINTdaxpy2v(const int n, const double a,
                 const double *x, const double *y, double *v)
{
        //cblas_dcopy(n, y, 1, v, 1);
        //cblas_daxpy(n, a, x, 1, v, 1);
        int i;
        for (i = 0; i < n; i++) {
                v[i] = a * x[i] + y[i];
        }
}


/*
 * a_t[n,m] = transpose of matrix a[m,n]
 */
void CINTdmat_transpose(double *a_t, const double *a, const int m, const int n)
{
        int i, j;

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
                default: for (i = 0; i < n; i++) a_t[i] = a[i];
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
                                a_t[i    ] = a[3*i+0];
                                a_t[i+n  ] = a[3*i+1];
                                a_t[i+n*2] = a[3*i+2];
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
                                a_t[i    ] = a[5*i+0];
                                a_t[i+n  ] = a[5*i+1];
                                a_t[i+n*2] = a[5*i+2];
                                a_t[i+n*3] = a[5*i+3];
                                a_t[i+n*4] = a[5*i+4];
                        }
                }
                break;
        default:
                switch (n) {
                case 1: for (i = 0; i < m; i++) a_t[i] = a[i]; break;
                case 3: for (i = 0; i < m; i++) {
                                a_t[3*i+0] = a[i    ];
                                a_t[3*i+1] = a[i+m  ];
                                a_t[3*i+2] = a[i+m*2];
                        }
                        break;
                case 5: for (i = 0; i < m; i++) {
                                a_t[5*i+0] = a[i    ];
                                a_t[5*i+1] = a[i+m  ];
                                a_t[5*i+2] = a[i+m*2];
                                a_t[5*i+3] = a[i+m*3];
                                a_t[5*i+4] = a[i+m*4];
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

void CINTzmat_transpose(double complex *a_t, const double complex *a,
                        const int m, const int n)
{
        int i, j;

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
                     const int m, const int n)
{
        int i, j;

        for (j = 0; j < n; j++) {
                for (i = 0; i < m; i++) {
                        a_t[i*n+j] = conj(a[j*m+i]);
                }
        }
}

