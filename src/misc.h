/*
 * File: misc.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic functions
 */
#include "fblas.h"

void CINTdcmplx_re(const unsigned int n, double *z, const double *re);
void CINTdcmplx_im(const unsigned int n, double *z, const double *im);
void CINTdcmplx_pp(const unsigned int n, double *z, const double *re, const double *im);
void CINTdcmplx_pn(const unsigned int n, double *z, const double *re, const double *im);
void CINTdcmplx_np(const unsigned int n, double *z, const double *re, const double *im);

double CINTsquare_dist(const double *r1, const double *r2);

void CINTrys_roots(const unsigned int nroots, double x, double *u, double *w);

double CINTgto_norm(int n, double a);

