#include "config.h"

void CINTrys_roots(FINT nroots, double x, double *u, double *w);
void CINTstg_roots(FINT nroots, double ta, double ua, double* rr, double* ww);
void CINTerfc_rys_roots(FINT nroots, double x, double lower, double *u, double *w);
void CINTerfc_rys_polyfits(FINT nroots, double x, double lower, double *u, double *w);

void gamma_inc_like(double *f, double t, FINT m);
#ifdef HAVE_QUADMATH_H
void qgamma_inc_like(__float128 *f, __float128 t, FINT m);
#else
void lgamma_inc_like(long double *f, long double t, FINT m);
#endif
void fmt_erfc_like(double *f, double t, double lower, FINT m);

#define MXROOTS   MXRYSROOTS
//#define MXROOTS1  (MXROOTS+1)
#define MXROOTS1  MXROOTS
