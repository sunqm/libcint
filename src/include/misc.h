/*
 * File: misc.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic functions
 */
#include "fblas.h"

void dcmplx_re(const int n, double *z, const double *re);
void dcmplx_im(const int n, double *z, const double *im);
void dcmplx_pp(const int n, double *z, const double *re, const double *im);
void dcmplx_pn(const int n, double *z, const double *re, const double *im);
void dcmplx_np(const int n, double *z, const double *re, const double *im);

inline double square_dist(const double *r1, const double *r2);

void rys_roots(const int nroots, double x, double *u, double *w);
