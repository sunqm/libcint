/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * blas interface and blas-like functions
 */

#if defined __cplusplus
extern "C" {
#endif
#include <complex.h>
#include "cint.h"
void CINTdset0(FINT n, double *x);
void CINTdaxpy2v(const FINT n, double a, double *x, double *y, double *v);
void CINTdmat_transpose(double *a_t, double *a, FINT m, FINT n);
void CINTdplus_transpose(double *a_t, double *a, FINT m, FINT n);
void CINTzmat_transpose(double complex *a_t, double complex *a, FINT m, FINT n);
void CINTzmat_dagger(double complex *a_c, double complex *a, FINT m, FINT n);

void CINTdgemm_NN(FINT m, FINT n, FINT k,
                  double *a, double *b, double *c);
void CINTdgemm_NN1(FINT m, FINT n, FINT k,
                   double *a, double *b, double *c, FINT ldc);
void CINTdgemm_TN(FINT m, FINT n, FINT k,
                  double *a, double *b, double *c);
void CINTdgemm_NT(FINT m, FINT n, FINT k,
                  double *a, double *b, double *c);
#if defined __cplusplus
} // end extern "C"
#endif
