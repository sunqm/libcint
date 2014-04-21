/*
 * File: fblas.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * blas interface and blas-like functions
 */

#if defined __cplusplus
extern "C" {
#endif

double dasum_(const unsigned int *n, const double *dx, const unsigned int *incx);
void dscal_(const unsigned int *n, const double *da, double *dx, const unsigned int *incx);
void daxpy_(const unsigned int *n, const double *da, const double *dx,
           const unsigned int *incx, double *dy, const unsigned int *incy);
double ddot_(const unsigned int *n, const double *dx, const unsigned int *incx,
             const double *dy, const unsigned int *incy);
void dcopy_(const unsigned int *n, const double *dx, const unsigned int *incx,
            const double *dy, const unsigned int *incy);
void dgemm_(const char*, const char*,
            const unsigned int*, const unsigned int*, const unsigned int*,
            const double*, const double*, const unsigned int*,
            const double*, const unsigned int*,
            const double*, double*, const unsigned int*);
void dgemv_(const char*, const unsigned int*, const unsigned int*,
            const double*, const double*, const unsigned int*,
            const double*, const unsigned int*,
            const double*, double*, const unsigned int*);
void dger_(const unsigned int *m, const unsigned int *n,
           const double *alpha, const double *x,
           const unsigned int *incx, const double *y, const unsigned int *incy,
           double *a, const unsigned int *lda);
void dsymm_(const char*, const char*, const unsigned int*, const unsigned int*,
            const double*, const double*, const unsigned int*,
            const double*, const unsigned int*,
            const double*, double*, const unsigned int*);

//void dsyrk_
void zgerc_(const unsigned int *m, const unsigned int *n,
            const double *alpha, const double *x, const unsigned int *incx,
            const double *y, const unsigned int *incy,
            double *a, const unsigned int *lda);
void zgemv_(const char*, const unsigned int*, const unsigned int*,
            const double*, const double*, const unsigned int*,
            const double*, const unsigned int*,
            const double*, double*, const unsigned int*);
void zgemm_(const char*, const char*,
            const unsigned int*, const unsigned int*, const unsigned int*,
            const double*, const double*, const unsigned int*,
            const double*, const unsigned int*,
            const double*, double*, const unsigned int*);


void CINTdset0(const unsigned int n, double *x);
void CINTdaxpy2v(const unsigned int n, const double a,
                 const double *x, const double *y, double *v);
void CINTdmat_transpose(double *a_t, const double *a,
                        const unsigned int m, const unsigned int n);
void CINTzmat_transpose(double *a_t, const double *a,
                        const unsigned int m, const unsigned int n);
void CINTzmat_dagger(double *a_c, const double *a,
                     const unsigned int m, const unsigned int n);

#if defined __cplusplus
} // end extern "C"
#endif
