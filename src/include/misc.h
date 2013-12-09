/*
 * File: misc.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic functions
 */
#include "fblas.h"

void dcmplx_re(const unsigned int n, double *z, const double *re);
void dcmplx_im(const unsigned int n, double *z, const double *im);
void dcmplx_pp(const unsigned int n, double *z, const double *re, const double *im);
void dcmplx_pn(const unsigned int n, double *z, const double *re, const double *im);
void dcmplx_np(const unsigned int n, double *z, const double *re, const double *im);

inline double square_dist(const double *r1, const double *r2);

void rys_roots(const unsigned int nroots, double x, double *u, double *w);
