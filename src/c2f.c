/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * c to fortran interface
 */

#ifdef WITH_FORTRAN
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

FINT cintlen_spinor_(const FINT *bas_id, const FINT *bas)
{
        return CINTlen_spinor(*bas_id, bas);
}

FINT cintcgtos_cart_(const FINT *bas_id, const FINT *bas)
{
        return CINTcgto_cart(*bas_id, bas);
}
FINT cintcgto_cart_(const FINT *bas_id, const FINT *bas)
{
        return CINTcgto_cart(*bas_id, bas);
}

FINT cintcgtos_spheric_(const FINT *bas_id, const FINT *bas)
{
        return CINTcgto_spheric(*bas_id, bas);
}
FINT cintcgto_spheric_(const FINT *bas_id, const FINT *bas)
{
        return CINTcgto_spheric(*bas_id, bas);
}

FINT cintcgtos_spinor_(const FINT *bas_id, const FINT *bas)
{
        return CINTcgto_spinor(*bas_id, bas);
}
FINT cintcgto_spinor_(const FINT *bas_id, const FINT *bas)
{
        return CINTcgto_spinor(*bas_id, bas);
}

/* 
 * tot. primitive atomic spheric GTOs in a shell
 */
FINT cinttot_pgto_spheric_(const FINT *bas, const FINT *nbas)
{
        return CINTtot_pgto_spheric(bas, *nbas);
}

/* 
 * tot. primitive atomic spinors in a shell
 */
FINT cinttot_pgto_spinor_(const FINT *bas, const FINT *nbas)
{
        return CINTtot_pgto_spinor(bas, *nbas);
}

/* 
 * tot. contracted atomic cartesian GTOs in a shell
 */
FINT cinttot_cgto_cart_(const FINT *bas, const FINT *nbas)
{
        return CINTtot_cgto_cart(bas, *nbas);
}

/* 
 * tot. contracted atomic spheric GTOs in a shell
 */
FINT cinttot_cgto_spheric_(const FINT *bas, const FINT *nbas)
{
        return CINTtot_cgto_spheric(bas, *nbas);
}

/* 
 * tot. contracted atomic spinors in a shell
 */
FINT cinttot_cgto_spinor_(const FINT *bas, const FINT *nbas)
{
        return CINTtot_cgto_spinor(bas, *nbas);
}

/* 
 * offset of each shell for cartesian GTOs
 */
void cintshells_cart_offset_(FINT ao_loc[], const FINT *bas, const FINT *nbas)
{
        CINTshells_cart_offset(ao_loc, bas, *nbas);
}

/* 
 * offset of each shell for real spheric GTOs
 */
void cintshells_spheric_offset_(FINT ao_loc[], const FINT *bas, const FINT *nbas)
{
        CINTshells_spheric_offset(ao_loc, bas, *nbas);
}

/* 
 * offset of each shell for AO spinors
 */
void cintshells_spinor_offset_(FINT ao_loc[], const FINT *bas, const FINT *nbas)
{
        CINTshells_spinor_offset(ao_loc, bas, *nbas);
}


double cintgto_norm_(FINT *n, double *a)
{
        return CINTgto_norm(*n, *a);
}

/*
 * * * * * * * * * * * * * * * * * * * * *
 * let Fortran be able to change CINTOpt
 */
/* in Fortran, pass an integer(8) to hold the pointer of CINTOpt */
//typedef int64_t CINTOptPtrAsInteger8;
void cintinit_2e_optimizer_(CINTOpt **opt,
                            FINT *atm, FINT *natm,
                            FINT *bas, FINT *nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, *natm, bas, *nbas, env);
}
void cintinit_optimizer_(CINTOpt **opt,
                         FINT *atm, FINT *natm,
                         FINT *bas, FINT *nbas, double *env)
{
        cintinit_2e_optimizer_(opt, atm, natm, bas, nbas, env);
}
void cintdel_2e_optimizer_(CINTOpt **opt)
{
        CINTdel_2e_optimizer(opt);
}
void cintdel_optimizer_(CINTOpt **opt)
{
        cintdel_2e_optimizer_(opt);
}
#endif
