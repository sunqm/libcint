/*
 * File: cint_bas.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO function
 */

#include "cint_const.h"

#define bas(SLOT,I)     bas[BAS_SLOTS * (I) + (SLOT)]
#define atm(SLOT,I)     atm[ATM_SLOTS * (I) + (SLOT)]

int CINTlen_cart(const int l);
int CINTlen_spinor(const int bas_id, const int *bas);

int CINTcgtos_cart(const int bas_id, const int *bas);
int CINTcgtos_spheric(const int bas_id, const int *bas);
int CINTcgtos_spinor(const int bas_id, const int *bas);
int CINTcgto_cart(const int bas_id, const int *bas);
int CINTcgto_spheric(const int bas_id, const int *bas);
int CINTcgto_spinor(const int bas_id, const int *bas);

int CINTtot_pgto_spheric(const int *bas, const int nbas);
int CINTtot_pgto_spinor(const int *bas, const int nbas);

int CINTtot_cgto_cart(const int *bas, const int nbas);
int CINTtot_cgto_spheric(const int *bas, const int nbas);
int CINTtot_cgto_spinor(const int *bas, const int nbas);

void CINTshells_cart_offset(int ao_loc[], const int *bas, const int nbas);
void CINTshells_spheric_offset(int ao_loc[], const int *bas, const int nbas);
void CINTshells_spinor_offset(int ao_loc[], const int *bas, const int nbas);

void CINTcart_comp(int *nx, int *ny, int *nz, const int lmax);

