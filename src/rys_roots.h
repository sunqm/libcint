#include "cint_config.h"

void CINTrys_roots(int nroots, double x, double *u, double *w);
void CINTsr_rys_roots(int nroots, double x, double lower, double *u, double *w);
void CINTstg_roots(int nroots, double ta, double ua, double* rr, double* ww);
int CINTsr_rys_polyfits(int nroots, double x, double lower, double *u, double *w);

int CINTrys_schmidt(int nroots, double x, double lower, double *roots, double *weights);
int CINTlrys_schmidt(int nroots, double x, double lower, double *roots, double *weights);
int CINTrys_laguerre(int n, double x, double lower, double *roots, double *weights);
int CINTlrys_laguerre(int n, double x, double lower, double *roots, double *weights);
int CINTrys_jacobi(int n, double x, double lower, double *roots, double *weights);
int CINTlrys_jacobi(int n, double x, double lower, double *roots, double *weights);
#ifdef HAVE_QUADMATH_H
int CINTqrys_schmidt(int nroots, double x, double lower, double *roots, double *weights);
int CINTqrys_laguerre(int n, double x, double lower, double *roots, double *weights);
int CINTqrys_jacobi(int n, double x, double lower, double *roots, double *weights);
#else
#define CINTqrys_schmidt        CINTlrys_schmidt
#define CINTqrys_laguerre       CINTlrys_laguerre
#define CINTqrys_jacobi         CINTlrys_jacobi
#endif

void gamma_inc_like(double *f, double t, int m);
void lgamma_inc_like(long double *f, long double t, int m);
//void fmt1_gamma_inc_like(double *f, double t, int m);
//void fmt1_lgamma_inc_like(long double *f, long double t, int m);
void fmt_erfc_like(double *f, double t, double lower, int m);
void fmt1_erfc_like(double *f, double t, double lower, int m);
void fmt_lerfc_like(long double *f, long double t, long double lower, int m);
void fmt1_lerfc_like(long double *f, long double t, long double lower, int m);
#ifdef HAVE_QUADMATH_H
void qgamma_inc_like(__float128 *f, __float128 t, int m);
void fmt_qerfc_like(__float128 *f, __float128 t, __float128 lower, int m);
void fmt1_qerfc_like(__float128 *f, __float128 t, __float128 lower, int m);
#else
#define qgamma_inc_like         lgamma_inc_like
#define fmt_qerfc_like          fmt_lerfc_like
#define fmt1_qerfc_like         fmt1_lerfc_like
#endif

// FIXME:
// short-range Coulomb kernel is numerically very instable when the integrals
// are close to zero (x*lower**2 > 40). Use this cutoff as a temporary solution
// to avoid the numerical issue in sr_rys_roots
#define EXPCUTOFF_SR    40
