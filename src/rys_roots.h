#include "config.h"

void CINTrys_roots(FINT nroots, double x, double *u, double *w);
static void Root123(FINT n, double x, double roots[], double weights[]);
static void Root4(double x, double roots[], double weights[]);
static void Root5(double x, double roots[], double weights[]);
static void R_droot(FINT nroots, double x, double roots[], double weights[]);
static void R_lroot(FINT nroots, double x, double roots[], double weights[]);
#ifdef HAVE_QUADMATH_H
static void R_qroot(FINT nroots, double x, double roots[], double weights[]);
#endif

void CINTstg_roots(int nroots, double ta, double ua, double* rr, double* ww);

