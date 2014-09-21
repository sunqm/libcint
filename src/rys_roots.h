/*
 * File: rys_roots.h
 *
 * Rys quadrature
 
 Code is edited based on
 * PyQuante quantum chemistry program suite http://pyquante.sourceforge.net
 * BDF program package

 */


void CINTrys_roots(int nroots, double x, double *u, double *w);
static void Root123(int n, double x, double roots[], double weights[]);
static void Root4(double x, double roots[], double weights[]);
static void Root5(double x, double roots[], double weights[]);
static void R_droot(int nroots, double x, double roots[], double weights[]);
static void R_qroot(int nroots, double x, double roots[], double weights[]);
