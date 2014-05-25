/*
 * File: cart2sph.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * Cartisen GTO to spheric or spinor GTO transformation
 */

/*************************************************
 *
 * transform matrix
 *
 *************************************************/
#include <complex.h>

void c2s_sph_1e(double *opij, const double *gctr,
                const int *shls, const int *bas);

void c2s_sph_2e1(double *fkijl, const double *gctr,
                 const int *shls, const int *bas);
void c2s_sph_2e2();

void c2s_sf_1e(double complex *opij, const double *gctr,
               const int *shls, const int *bas);
void c2s_sf_1ei(double complex *opij, const double *gctr,
                const int *shls, const int *bas);

void c2s_si_1e(double complex *opij, const double *gctr,
               const int *shls, const int *bas);
void c2s_si_1ei(double complex *opij, const double *gctr,
                const int *shls, const int *bas);

void c2s_sf_2e1(double complex *opij, const double *gctr,
                const int *shls, const int *bas);
void c2s_sf_2e1i(double complex *opij, const double *gctr,
                 const int *shls, const int *bas);

void c2s_sf_2e2(double complex *fkijl, const double complex *opij,
                const int *shls, const int *bas);
void c2s_sf_2e2i(double complex *fkijl, const double complex *opij,
                 const int *shls, const int *bas);

void c2s_si_2e1(double complex *opij, const double *gctr,
                const int *shls, const int *bas);
void c2s_si_2e1i(double complex *opij, const double *gctr,
                 const int *shls, const int *bas);

void c2s_si_2e2(double complex *fkijl, const double complex *opij,
                const int *shls, const int *bas);
void c2s_si_2e2i(double complex *fkijl, const double complex *opij,
                 const int *shls, const int *bas);

void c2s_cart_1e(double *opij, const double *gctr,
                 const int *shls, const int *bas);
void c2s_cart_2e1(double *fkijl, const double *gctr,
                  const int *shls, const int *bas);
void c2s_cart_2e2();

