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
void c2s_sph_1e(double *opij, const double *gctr,
                const unsigned int *shls, const int *bas);

void c2s_sph_2e1(double *fkijl, const double *gctr,
                 const unsigned int *shls, const int *bas);
void c2s_sph_2e2();

void c2s_sf_1e(double *opij, const double *gctr,
               const unsigned int *shls, const int *bas);
void c2s_sf_1ei(double *opij, const double *gctr,
                const unsigned int *shls, const int *bas);

void c2s_si_1e(double *opij, const double *gctr,
               const unsigned int *shls, const int *bas);
void c2s_si_1ei(double *opij, const double *gctr,
                const unsigned int *shls, const int *bas);

void c2s_sf_2e1(double *opij, const double *gctr,
                const unsigned int *shls, const int *bas);
void c2s_sf_2e1i(double *opij, const double *gctr,
                 const unsigned int *shls, const int *bas);

void c2s_sf_2e2(double *fkijl, const double *opij,
                const unsigned int *shls, const int *bas);
void c2s_sf_2e2i(double *fkijl, const double *opij,
                 const unsigned int *shls, const int *bas);

void c2s_si_2e1(double *opij, const double *gctr,
                const unsigned int *shls, const int *bas);
void c2s_si_2e1i(double *opij, const double *gctr,
                 const unsigned int *shls, const int *bas);

void c2s_si_2e2(double *fkijl, const double *opij,
                const unsigned int *shls, const int *bas);
void c2s_si_2e2i(double *fkijl, const double *opij,
                 const unsigned int *shls, const int *bas);

void c2s_cart_1e(double *opij, const double *gctr,
                 const unsigned int *shls, const int *bas);
void c2s_cart_2e1(double *fkijl, const double *gctr,
                  const unsigned int *shls, const int *bas);
void c2s_cart_2e2();

