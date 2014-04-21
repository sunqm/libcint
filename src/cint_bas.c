/*
 * File: cint_bas.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO function
 */

#include "cint_bas.h"
/*
#define DEBUG_ON(i, dj, l, k, s) printf("%d : %d, l = %d, kappa = %d, %s\n", \
                                 (i), (dj), (l), (k), (s));
*/
#define DEBUG_ON(...)


/*
 * No. components of a Cartesian GTO, = (l+1)*(l+2)/2
 */
unsigned int CINTlen_cart(const unsigned int l)
{
        switch (l) {
                case 0:
                        return 1;
                case 1:
                        return 3;
                case 2:
                        return 6;
                case 3:
                        return 10;
                case 4:
                        return 15;
                default:
                        return (l + 1) * (l + 2) / 2;
        }
}

unsigned int CINTlen_spinor(const unsigned int bas_id, const int *bas)
{
        if (0 == bas(KAPPA_OF, bas_id)) {
                return 4 * bas(ANG_OF, bas_id) + 2;
        } else if (bas(KAPPA_OF, bas_id) < 0) {
                return 2 * bas(ANG_OF, bas_id) + 2;
        } else {
                return 2 * bas(ANG_OF, bas_id);
        }
}

/* 
 * Num. of contracted cartesian GTO = 2j+1 * n_contraction
 */
unsigned int CINTcgtos_cart(const unsigned int bas_id, const int *bas)
{
        return CINTlen_cart(bas(ANG_OF, bas_id)) * bas(NCTR_OF, bas_id);
}

/* 
 * Num. of contracted spheric GTO = 2j+1 * n_contraction
 */
unsigned int CINTcgtos_spheric(const unsigned int bas_id, const int *bas)
{
        return (bas(ANG_OF, bas_id) * 2 + 1) * bas(NCTR_OF, bas_id);
}

/* 
 * Num. of contracted spinor GTO
 */
unsigned int CINTcgtos_spinor(const unsigned int bas_id, const int *bas)
{
        return CINTlen_spinor(bas_id, bas) * bas(NCTR_OF, bas_id);
}

/*
 * tot. primitive atomic spheric GTOs in a shell
 */
unsigned int CINTtot_pgto_spheric(const int *bas, const int nbas)
{
        unsigned int i;
        unsigned int s = 0;

        for (i = 0; i < nbas; i++) {
                s += (bas(ANG_OF, i) * 2 + 1)
                        * bas(NPRIM_OF, i);
        }
        return s;
}

/*
 * tot. primitive atomic spinors in a shell
 */
unsigned int CINTtot_pgto_spinor(const int *bas, const int nbas)
{
        unsigned int i;
        unsigned int s = 0;

        for (i = 0; i < nbas; i++) {
                s += CINTlen_spinor(i, bas) * bas(NPRIM_OF, i);
        }
        return s;
}

static unsigned int tot_cgto_accum(unsigned int (*f)(), const int *bas, const int nbas)
{
        unsigned int i;
        unsigned int s = 0;

        for (i = 0; i < nbas; i++) {
                s += (*f)(i, bas);
        }
        return s;
}
/*
 * tot. contracted atomic spheric GTOs in a shell
 */
unsigned int CINTtot_cgto_spheric(const int *bas, const int nbas)
{
        return tot_cgto_accum(&CINTcgtos_spheric, bas, nbas);
}

/*
 * tot. contracted atomic spinors in a shell
 */
unsigned int CINTtot_cgto_spinor(const int *bas, const int nbas)
{
        return tot_cgto_accum(&CINTcgtos_spinor, bas, nbas);
}

/*
 * tot. contracted atomic spinors in a shell
 */
unsigned int CINTtot_cgto_cart(const int *bas, const int nbas)
{
        return tot_cgto_accum(&CINTcgtos_cart, bas, nbas);
}

static void shells_cgto_offset(unsigned int (*f)(), int ao_loc[],
                                   const int *bas, const int nbas)
{
        unsigned int i, s;
        for (i = 0, s = 0; i < nbas; i++) {
                ao_loc[i] = s;
                s += (*f)(i, bas);
        }
}
/*
 * offset of each shell for real spheric GTOs
 */
void CINTshells_cart_offset(int ao_loc[], const int *bas, const int nbas)
{
        shells_cgto_offset(&CINTcgtos_cart, ao_loc, bas, nbas);
}

/*
 * offset of each shell for real spheric GTOs
 */
void CINTshells_spheric_offset(int ao_loc[], const int *bas, const int nbas)
{
        shells_cgto_offset(&CINTcgtos_spheric, ao_loc, bas, nbas);
}

/*
 * offset of each shell for AO spinors
 */
void CINTshells_spinor_offset(int ao_loc[], const int *bas, const int nbas)
{
        shells_cgto_offset(&CINTcgtos_spinor, ao_loc, bas, nbas);
}


/*
 * GTO = x^{nx}y^{ny}z^{nz}e^{-ar^2}
 */
void CINTcart_comp(unsigned int *nx, unsigned int *ny, unsigned int *nz,
                   const unsigned int lmax)
{
        unsigned int inc = 0;
        int lx, ly, lz;

        for (lx = lmax; lx >= 0; lx--) {
                for (ly = lmax - lx; ly >= 0; ly--) {
                        lz = lmax - lx - ly;
                        nx[inc] = lx;
                        ny[inc] = ly;
                        nz[inc] = lz;
                        inc++;
                }
        }
}

