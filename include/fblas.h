/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * blas interface and blas-like functions
 */

#if defined __cplusplus
extern "C" {
#endif
#include <complex.h>
#include "config.h"

double dasum_(const FINT *n, const double *dx, const FINT *incx);
void dscal_(const FINT *n, const double *da, double *dx, const FINT *incx);
void daxpy_(const FINT *n, const double *da, const double *dx,
           const FINT *incx, double *dy, const FINT *incy);
double ddot_(const FINT *n, const double *dx, const FINT *incx,
             const double *dy, const FINT *incy);
void dcopy_(const FINT *n, const double *dx, const FINT *incx,
            const double *dy, const FINT *incy);
void dgemm_(const char*, const char*,
            const FINT*, const FINT*, const FINT*,
            const double*, const double*, const FINT*,
            const double*, const FINT*,
            const double*, double*, const FINT*);
void dgemv_(const char*, const FINT*, const FINT*,
            const double*, const double*, const FINT*,
            const double*, const FINT*,
            const double*, double*, const FINT*);
void dger_(const FINT *m, const FINT *n,
           const double *alpha, const double *x,
           const FINT *incx, const double *y, const FINT *incy,
           double *a, const FINT *lda);
void dsymm_(const char*, const char*, const FINT*, const FINT*,
            const double*, const double*, const FINT*,
            const double*, const FINT*,
            const double*, double*, const FINT*);

//void dsyrk_
void zgerc_(const FINT *m, const FINT *n,
            const double complex *alpha, const double complex *x, const FINT *incx,
            const double complex *y, const FINT *incy,
            double complex *a, const FINT *lda);
void zgemv_(const char*, const FINT*, const FINT*,
            const double complex*, const double complex*, const FINT*,
            const double complex*, const FINT*,
            const double complex*, double complex*, const FINT*);
void zgemm_(const char*, const char*,
            const FINT*, const FINT*, const FINT*,
            const double complex*, const double complex*, const FINT*,
            const double complex*, const FINT*,
            const double complex*, double complex*, const FINT*);


void CINTdset0(const FINT n, double *x);
void CINTdaxpy2v(const FINT n, const double a,
                 const double *x, const double *y, double *v);
void CINTdmat_transpose(double *a_t, const double *a, const FINT m, const FINT n);
void CINTzmat_transpose(double complex *a_t, const double complex *a,
                        const FINT m, const FINT n);
void CINTzmat_dagger(double complex *a_c, const double complex *a,
                     const FINT m, const FINT n);

#if defined __cplusplus
} // end extern "C"
#endif
