/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * c to fortran interface
 */

#include <stdlib.h>
#include <math.h>
#include "cint_bas.h"
#include "g1e.h"
#include "misc.h"
#include "c2f.h"
#include "optimizer.h"


/*
 * * * * * * * * * * * * * * * * * * * * *
 * for cint_bas.c
 */

int cintlen_spinor_(const int *bas_id, const int *bas)
{
        return CINTlen_spinor(*bas_id, bas);
}

int cintcgtos_cart_(const int *bas_id, const int *bas)
{
        return CINTcgto_cart(*bas_id, bas);
}
int cintcgto_cart_(const int *bas_id, const int *bas)
{
        return CINTcgto_cart(*bas_id, bas);
}

int cintcgtos_spheric_(const int *bas_id, const int *bas)
{
        return CINTcgto_spheric(*bas_id, bas);
}
int cintcgto_spheric_(const int *bas_id, const int *bas)
{
        return CINTcgto_spheric(*bas_id, bas);
}

int cintcgtos_spinor_(const int *bas_id, const int *bas)
{
        return CINTcgto_spinor(*bas_id, bas);
}
int cintcgto_spinor_(const int *bas_id, const int *bas)
{
        return CINTcgto_spinor(*bas_id, bas);
}

/* 
 * tot. primitive atomic spheric GTOs in a shell
 */
int cinttot_pgto_spheric_(const int *bas, const int *nbas)
{
        return CINTtot_pgto_spheric(bas, *nbas);
}

/* 
 * tot. primitive atomic spinors in a shell
 */
int cinttot_pgto_spinor_(const int *bas, const int *nbas)
{
        return CINTtot_pgto_spinor(bas, *nbas);
}

/* 
 * tot. contracted atomic cartesian GTOs in a shell
 */
int cinttot_cgto_cart_(const int *bas, const int *nbas)
{
        return CINTtot_cgto_cart(bas, *nbas);
}

/* 
 * tot. contracted atomic spheric GTOs in a shell
 */
int cinttot_cgto_spheric_(const int *bas, const int *nbas)
{
        return CINTtot_cgto_spheric(bas, *nbas);
}

/* 
 * tot. contracted atomic spinors in a shell
 */
int cinttot_cgto_spinor_(const int *bas, const int *nbas)
{
        return CINTtot_cgto_spinor(bas, *nbas);
}

/* 
 * offset of each shell for cartesian GTOs
 */
void cintshells_cart_offset_(int ao_loc[], const int *bas, const int *nbas)
{
        CINTshells_cart_offset(ao_loc, bas, *nbas);
}

/* 
 * offset of each shell for real spheric GTOs
 */
void cintshells_spheric_offset_(int ao_loc[], const int *bas, const int *nbas)
{
        CINTshells_spheric_offset(ao_loc, bas, *nbas);
}

/* 
 * offset of each shell for AO spinors
 */
void cintshells_spinor_offset_(int ao_loc[], const int *bas, const int *nbas)
{
        CINTshells_spinor_offset(ao_loc, bas, *nbas);
}


/* 
 * GTO = x^{nx}y^{ny}z^{nz}e^{-ar^2}
 */
void cintcart_comp_(int *nx, int *ny, int *nz, const int *lmax)
{
        CINTcart_comp(nx, ny, nz, *lmax);
}


double cintgto_norm_(int *n, double *a)
{
        return CINTgto_norm(*n, *a);
}

/*
 * * * * * * * * * * * * * * * * * * * * *
 * let Fortran be able to change CINTOpt
 */
/* in Fortran, pass an integer(8) to hold the pointer of CINTOpt */
//typedef long CINTOptPtrAsInteger8;
void cintinit_2e_optimizer_(CINTOptPtrAsInteger8 *optptr,
                            int *atm, int *natm,
                            int *bas, int *nbas, double *env)
{
        CINTOpt **opt = (CINTOpt **)optptr;
        CINTinit_2e_optimizer(opt, atm, *natm, bas, *nbas, env);
}
void cintinit_optimizer_(CINTOptPtrAsInteger8 *optptr,
                         int *atm, int *natm,
                         int *bas, int *nbas, double *env)
{
        cintinit_2e_optimizer_(optptr, atm, natm, bas, nbas, env);
}
void cintdel_2e_optimizer_(CINTOptPtrAsInteger8 *optptr)
{
        CINTOpt **opt = (CINTOpt **)optptr;
        CINTdel_2e_optimizer(opt);
}
void cintdel_optimizer_(CINTOptPtrAsInteger8 *optptr)
{
        cintdel_2e_optimizer_(optptr);
}
