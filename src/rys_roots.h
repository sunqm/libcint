/*
 * File: rys_roots.h
 *
 * Rys quadrature
 
 This code is edited from "crys.c" of PyQuante quantum chemistry program suite.
 http://pyquante.sourceforge.net

 */


void CINTrys_roots(const unsigned int nroots, double x, double *u, double *w);
static void Root123(unsigned int n, double X, double roots[], double weights[]);
static void Root4(double X, double roots[], double weights[]);
static void Root5(double X, double roots[], double weights[]);
static void R_droot(const unsigned int nroots, const double xx, double roots[], double weights[]);
