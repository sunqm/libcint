/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic functions
 */
#include "fblas.h"

void CINTdcmplx_re(const int n, double complex *z, const double *re);
void CINTdcmplx_im(const int n, double complex *z, const double *im);
void CINTdcmplx_pp(const int n, double complex *z, const double *re, const double *im);
void CINTdcmplx_pn(const int n, double complex *z, const double *re, const double *im);
void CINTdcmplx_np(const int n, double complex *z, const double *re, const double *im);

double CINTsquare_dist(const double *r1, const double *r2);

void CINTrys_roots(const int nroots, double x, double *u, double *w);

double CINTgto_norm(int n, double a);

