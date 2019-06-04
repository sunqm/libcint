/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Cartisen GTO to spheric or spinor GTO transformation
 */

/*************************************************
 *
 * transform matrix
 *
 *************************************************/
#include <complex.h>
#include "g1e.h"

void c2s_sph_1e(double *opij, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_sph_2e1(double *fijkl, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_sph_2e2();

void c2s_cart_1e(double *opij, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_cart_2e1(double *fijkl, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_cart_2e2();

void c2s_sf_1e(double complex *opij, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_sf_1ei(double complex *opij, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);

void c2s_si_1e(double complex *opij, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_si_1ei(double complex *opij, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);

void c2s_sf_2e1(double complex *opij, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_sf_2e1i(double complex *opij, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);

void c2s_sf_2e2(double complex *fijkl, double complex *opij, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_sf_2e2i(double complex *fijkl, double complex *opij, FINT *dims, CINTEnvVars *envs, double *cache);

void c2s_si_2e1(double complex *opij, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_si_2e1i(double complex *opij, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);

void c2s_si_2e2(double complex *fijkl, double complex *opij, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_si_2e2i(double complex *fijkl, double complex *opij, FINT *dims, CINTEnvVars *envs, double *cache);

void c2s_sph_3c2e1(double *fijkl, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_cart_3c2e1(double *fijkl, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_sph_3c2e1_ssc(double *fijkl, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);

void c2s_sf_3c2e1(double complex *opijk, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_sf_3c2e1i(double complex *opijk, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_si_3c2e1(double complex *opijk, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_si_3c2e1i(double complex *opijk, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_sf_3c2e1_ssc(double complex *opijk, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_sf_3c2e1i_ssc(double complex *opijk, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_si_3c2e1_ssc(double complex *opijk, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_si_3c2e1i_ssc(double complex *opijk, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);

void c2s_sph_3c1e(double *fijkl, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);
void c2s_cart_3c1e(double *fijkl, double *gctr, FINT *dims, CINTEnvVars *envs, double *cache);

void c2s_dset0(double *out, FINT *dims, FINT *counts);
void c2s_zset0(double complex *out, FINT *dims, FINT *counts);

/*************************************************
 *
 * transform vectors
 *
 *************************************************/
void c2s_sph_vec(double *sph, double *cart, FINT l, FINT nvec);
