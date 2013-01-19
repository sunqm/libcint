/*
 * File: cint_bas.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO function
 */

#include "cint_const.h"

#define bas(SLOT,I)     bas[BAS_SLOTS * (I) + (SLOT)]
#define atm(SLOT,I)     atm[ATM_SLOTS * (I) + (SLOT)]


int len_cart(const int l);
int len_spinor(const int bas_id, const int *bas);

int cgtos_cart(const int bas_id, const int *bas);
int cgtos_spheric(const int bas_id, const int *bas);
int cgtos_spinor(const int bas_id, const int *bas);

int tot_pgto_spheric(const int *bas, const int nbas);
int tot_pgto_spinor(const int *bas, const int nbas);

int tot_cgto_cart(const int *bas, const int nbas);
int tot_cgto_spheric(const int *bas, const int nbas);
int tot_cgto_spinor(const int *bas, const int nbas);

void shells_cart_offset(int ao_loc[], const int *bas, const int nbas);
void shells_spheric_offset(int ao_loc[], const int *bas, const int nbas);
void shells_spinor_offset(int ao_loc[], const int *bas, const int nbas);

void cart_comp(int *nx, int *ny, int *nz, const int lmax);
