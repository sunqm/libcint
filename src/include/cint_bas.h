/*
 * File: cint_bas.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO function
 */

#include "cint_const.h"

#define bas(SLOT,I)     bas[BAS_SLOTS * (I) + (SLOT)]
#define atm(SLOT,I)     atm[ATM_SLOTS * (I) + (SLOT)]

// ref to init_int1e_CintEnvVars, init_int2e_CintEnvVars
typedef struct {
        int natm;
        int nbas;
        const int *atm;
        const int *bas;
        const double *env;
        double *_mem_align_padding1;
        const unsigned int *shls;
        const unsigned int *ng;

        unsigned int i_l;
        unsigned int j_l;
        unsigned int k_l;
        unsigned int l_l;
        unsigned int i_prim;
        unsigned int j_prim;
        unsigned int k_prim;
        unsigned int l_prim;
        unsigned int i_ctr;
        unsigned int j_ctr;
        unsigned int k_ctr;
        unsigned int l_ctr;
        unsigned int nfi;  // number of cartesion components
        unsigned int nfj;
        unsigned int nfk;
        unsigned int nfl;
        unsigned int nf;  // = nfi*nfj*nfk*nfl;
        unsigned int g_size;  // ref to cint2e.c g = malloc(sizeof(double)*g_size)
        unsigned int g_stride_i; // ng[RYS_ROOTS]  shift of (i++,k,l,j)
        unsigned int g_stride_k; // ng[RYS_ROOTS]*ng[0]  shift of (i,k++,l,j)
        unsigned int g_stride_l; // ng[RYS_ROOTS]*ng[0]*ng[1]  shift of (i,k,l++,j)
        unsigned int g_stride_j; // ng[RYS_ROOTS]*ng[0]*ng[1]*ng[2]  shift of (i,k,l,j++)

        const double *ri;
        const double *rj;
        const double *rk;
        const double *rl;
        double ai;
        double aj;
        double ak;
        double al;
} CintEnvVars;


unsigned int len_cart(const unsigned int l);
unsigned int len_spinor(const unsigned int bas_id, const int *bas);

unsigned int cgtos_cart(const unsigned int bas_id, const int *bas);
unsigned int cgtos_spheric(const unsigned int bas_id, const int *bas);
unsigned int cgtos_spinor(const unsigned int bas_id, const int *bas);

unsigned int tot_pgto_spheric(const int *bas, const int nbas);
unsigned int tot_pgto_spinor(const int *bas, const int nbas);

unsigned int tot_cgto_cart(const int *bas, const int nbas);
unsigned int tot_cgto_spheric(const int *bas, const int nbas);
unsigned int tot_cgto_spinor(const int *bas, const int nbas);

void shells_cart_offset(int ao_loc[], const int *bas, const int nbas);
void shells_spheric_offset(int ao_loc[], const int *bas, const int nbas);
void shells_spinor_offset(int ao_loc[], const int *bas, const int nbas);

void cart_comp(unsigned int *nx, unsigned int *ny, unsigned int *nz,
               const unsigned int lmax);
