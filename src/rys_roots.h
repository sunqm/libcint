#include "config.h"
#include "cint_const.h"

#define MXROOTS   MXRYSROOTS
//#define MXROOTS1  (MXROOTS+1)
#define MXROOTS1  MXROOTS

void CINTrys_roots(FINT nroots, double x, double *u, double *w);
void CINTstg_roots(FINT nroots, double ta, double ua, double* rr, double* ww);
void CINTsr_rys_roots(FINT nroots, double x, double lower, double *u, double *w);
void CINTsr_rys_polyfits(FINT nroots, double x, double lower, double *u, double *w);
int CINTR_droot(FINT nroots, double x, double *roots, double *weights);
#ifdef HAVE_QUADMATH_H
int CINTR_qroot(FINT nroots, double x, double *roots, double *weights);
#else
int CINTR_lroot(FINT nroots, double x, double *roots, double *weights);
#endif

void gamma_inc_like(double *f, double t, FINT m);
void lgamma_inc_like(long double *f, long double t, FINT m);
//void fmt1_gamma_inc_like(double *f, double t, FINT m);
//void fmt1_lgamma_inc_like(long double *f, long double t, FINT m);
void fmt_erfc_like(double *f, double t, double lower, FINT m);
void fmt1_erfc_like(double *f, double t, double lower, FINT m);
void fmt1_lerfc_like(long double *f, long double t, long double lower, FINT m);
#ifdef HAVE_QUADMATH_H
void qgamma_inc_like(__float128 *f, __float128 t, FINT m);
void fmt1_qerfc_like(__float128 *f, __float128 t, __float128 lower, FINT m);
#endif

int CINTrys_wheeler(int n, double x, double lower, double *roots, double *weights);
