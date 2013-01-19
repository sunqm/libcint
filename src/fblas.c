/*
 * File: fblas.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * blas-like functions
 */

#define OF_CMPLX        2

void dset0(const int n, double *x)
{
        int i;
/* for icc
 * #pragma simd
 * #pragma vector aligned
 */
        for (i = 0; i < n; i++)
                x[i] = 0;
}


/*
 * v = a * x + y
 */
void daxpy2v(const int n, const double a, const double *x, const double *y, 
             double *v)
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
void dmat_transpose(double *a_t, double *a, int m, int n)
{
        int i, j;
        double *pa_t;

        if (m == 1) {
                for (i = 0; i < n; i++)
                        a_t[i] = a[i];
        } else if (n == 1) {
                for (i = 0; i < m; i++)
                        a_t[i] = a[i];
        } else {
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
void zmat_transpose(double *a_t, const double *a, int m, int n)
{
        int i, j;
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
void zmat_dagger(double *a_t, const double *a, int m, int n)
{
        int i, j;
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
