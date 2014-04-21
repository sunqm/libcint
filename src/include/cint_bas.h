/*
 * File: cint_bas.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO function
 */

#include "cint_const.h"

#if !defined HAVE_DEFINED_CINTOPT_H
#define HAVE_DEFINED_CINTOPT_H
typedef struct {
    unsigned int **index_xyz_array; // ANG_MAX**4 pointers to index_xyz
    unsigned int *ptr_log_coeff;
    double *log_coeff; // -log(c) where c is the largest coeff of a pGTO
    double expcutoff;
} CINTOpt;
#endif

#define bas(SLOT,I)     bas[BAS_SLOTS * (I) + (SLOT)]
#define atm(SLOT,I)     atm[ATM_SLOTS * (I) + (SLOT)]

unsigned int CINTlen_cart(const unsigned int l);
unsigned int CINTlen_spinor(const unsigned int bas_id, const int *bas);

unsigned int CINTcgtos_cart(const unsigned int bas_id, const int *bas);
unsigned int CINTcgtos_spheric(const unsigned int bas_id, const int *bas);
unsigned int CINTcgtos_spinor(const unsigned int bas_id, const int *bas);

unsigned int CINTtot_pgto_spheric(const int *bas, const int nbas);
unsigned int CINTtot_pgto_spinor(const int *bas, const int nbas);

unsigned int CINTtot_cgto_cart(const int *bas, const int nbas);
unsigned int CINTtot_cgto_spheric(const int *bas, const int nbas);
unsigned int CINTtot_cgto_spinor(const int *bas, const int nbas);

void CINTshells_cart_offset(int ao_loc[], const int *bas, const int nbas);
void CINTshells_spheric_offset(int ao_loc[], const int *bas, const int nbas);
void CINTshells_spinor_offset(int ao_loc[], const int *bas, const int nbas);

void CINTcart_comp(unsigned int *nx, unsigned int *ny, unsigned int *nz,
                   const unsigned int lmax);

// misc.c
double CINTgto_norm(int n, double a);
