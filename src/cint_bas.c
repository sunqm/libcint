/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO function
 */

#include "cint_bas.h"

/*
 * No. components of a Cartesian GTO, = (l+1)*(l+2)/2
 */
FINT CINTlen_cart(const FINT l)
{
        return (l + 1) * (l + 2) / 2;
}

FINT CINTlen_spinor(const FINT bas_id, const FINT *bas)
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
FINT CINTcgtos_cart(const FINT bas_id, const FINT *bas)
{
        FINT l = bas(ANG_OF, bas_id);
        return (l+1)*(l+2)/2 * bas(NCTR_OF, bas_id);
}
FINT CINTcgto_cart(const FINT bas_id, const FINT *bas)
{
        FINT l = bas(ANG_OF, bas_id);
        return (l+1)*(l+2)/2 * bas(NCTR_OF, bas_id);
}

/* 
 * Num. of contracted spheric GTO = 2j+1 * n_contraction
 */
FINT CINTcgtos_spheric(const FINT bas_id, const FINT *bas)
{
        return (bas(ANG_OF, bas_id) * 2 + 1) * bas(NCTR_OF, bas_id);
}
FINT CINTcgto_spheric(const FINT bas_id, const FINT *bas)
{
        return (bas(ANG_OF, bas_id) * 2 + 1) * bas(NCTR_OF, bas_id);
}

/* 
 * Num. of contracted spinor GTO
 */
FINT CINTcgtos_spinor(const FINT bas_id, const FINT *bas)
{
        return CINTlen_spinor(bas_id, bas) * bas(NCTR_OF, bas_id);
}
FINT CINTcgto_spinor(const FINT bas_id, const FINT *bas)
{
        return CINTlen_spinor(bas_id, bas) * bas(NCTR_OF, bas_id);
}

/*
 * tot. primitive atomic spheric GTOs in a shell
 */
FINT CINTtot_pgto_spheric(const FINT *bas, const FINT nbas)
{
        FINT i;
        FINT s = 0;

        for (i = 0; i < nbas; i++) {
                s += (bas(ANG_OF, i) * 2 + 1)
                        * bas(NPRIM_OF, i);
        }
        return s;
}

/*
 * tot. primitive atomic spinors in a shell
 */
FINT CINTtot_pgto_spinor(const FINT *bas, const FINT nbas)
{
        FINT i;
        FINT s = 0;

        for (i = 0; i < nbas; i++) {
                s += CINTlen_spinor(i, bas) * bas(NPRIM_OF, i);
        }
        return s;
}

static FINT tot_cgto_accum(FINT (*f)(), const FINT *bas, const FINT nbas)
{
        FINT i;
        FINT s = 0;

        for (i = 0; i < nbas; i++) {
                s += (*f)(i, bas);
        }
        return s;
}
/*
 * tot. contracted atomic spheric GTOs in a shell
 */
FINT CINTtot_cgto_spheric(const FINT *bas, const FINT nbas)
{
        return tot_cgto_accum(&CINTcgto_spheric, bas, nbas);
}

/*
 * tot. contracted atomic spinors in a shell
 */
FINT CINTtot_cgto_spinor(const FINT *bas, const FINT nbas)
{
        return tot_cgto_accum(&CINTcgto_spinor, bas, nbas);
}

/*
 * tot. contracted atomic spinors in a shell
 */
FINT CINTtot_cgto_cart(const FINT *bas, const FINT nbas)
{
        return tot_cgto_accum(&CINTcgto_cart, bas, nbas);
}

static void shells_cgto_offset(FINT (*f)(), FINT ao_loc[],
                               const FINT *bas, const FINT nbas)
{
        FINT i;
        ao_loc[0] = 0;
        for (i = 1; i < nbas; i++) {
                ao_loc[i] = ao_loc[i-1] + (*f)(i-1, bas);
        }
}
/*
 * offset of each shell for real spheric GTOs
 */
void CINTshells_cart_offset(FINT ao_loc[], const FINT *bas, const FINT nbas)
{
        shells_cgto_offset(&CINTcgto_cart, ao_loc, bas, nbas);
}

/*
 * offset of each shell for real spheric GTOs
 */
void CINTshells_spheric_offset(FINT ao_loc[], const FINT *bas, const FINT nbas)
{
        shells_cgto_offset(&CINTcgto_spheric, ao_loc, bas, nbas);
}

/*
 * offset of each shell for AO spinors
 */
void CINTshells_spinor_offset(FINT ao_loc[], const FINT *bas, const FINT nbas)
{
        shells_cgto_offset(&CINTcgto_spinor, ao_loc, bas, nbas);
}


/*
 * GTO = x^{nx}y^{ny}z^{nz}e^{-ar^2}
 */
void CINTcart_comp(FINT *nx, FINT *ny, FINT *nz, const FINT lmax)
{
        FINT inc = 0;
        FINT lx, ly, lz;

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

