/*
 * File: c2f.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * c to fortran interface
 */

#include <stdlib.h>
#include <math.h>
#include "cint_bas.h"
#include "g1e.h"
#include "misc.h"
#include "c2f.h"


/*
 * * * * * * * * * * * * * * * * * * * * *
 * for cint_bas.c
 */

int len_spinor_(const unsigned int *bas_id, const int *bas)
{
        return len_spinor(*bas_id, bas);
}

int cgtos_cart_(const unsigned int *bas_id, const int *bas)
{
        return cgtos_cart(*bas_id, bas);
}

int cgtos_spheric_(const unsigned int *bas_id, const int *bas)
{
        return cgtos_spheric(*bas_id, bas);
}

int cgtos_spinor_(const unsigned int *bas_id, const int *bas)
{
        return cgtos_spinor(*bas_id, bas);
}

/* 
 * tot. primitive atomic spheric GTOs in a shell
 */
int tot_pgto_spheric_(const int *bas, const int *nbas)
{
        return tot_pgto_spheric(bas, *nbas);
}

/* 
 * tot. primitive atomic spinors in a shell
 */
int tot_pgto_spinor_(const int *bas, const int *nbas)
{
        return tot_pgto_spinor(bas, *nbas);
}

/* 
 * tot. contracted atomic cartesian GTOs in a shell
 */
int tot_cgto_cart_(const int *bas, const int *nbas)
{
        return tot_cgto_cart(bas, *nbas);
}

/* 
 * tot. contracted atomic spheric GTOs in a shell
 */
int tot_cgto_spheric_(const int *bas, const int *nbas)
{
        return tot_cgto_spheric(bas, *nbas);
}

/* 
 * tot. contracted atomic spinors in a shell
 */
int tot_cgto_spinor_(const int *bas, const int *nbas)
{
        return tot_cgto_spinor(bas, *nbas);
}

/* 
 * offset of each shell for cartesian GTOs
 */
void shells_cart_offset_(int ao_loc[], const int *bas, const int *nbas)
{
        shells_cart_offset(ao_loc, bas, *nbas);
}

/* 
 * offset of each shell for real spheric GTOs
 */
void shells_spheric_offset_(int ao_loc[], const int *bas, const int *nbas)
{
        shells_spheric_offset(ao_loc, bas, *nbas);
}

/* 
 * offset of each shell for AO spinors
 */
void shells_spinor_offset_(int ao_loc[], const int *bas, const int *nbas)
{
        shells_spinor_offset(ao_loc, bas, *nbas);
}


/* 
 * GTO = x^{nx}y^{ny}z^{nz}e^{-ar^2}
 */
void cart_comp_(unsigned int *nx, unsigned int *ny, unsigned int *nz,
                const unsigned int *lmax)
{
        cart_comp(nx, ny, nz, *lmax);
}

