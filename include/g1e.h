/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include "config.h"

#if !defined HAVE_DEFINED_CINTENVVARS_H
#define HAVE_DEFINED_CINTENVVARS_H
// ref to CINTinit_int1e_EnvVars, CINTinit_int2e_EnvVars
typedef struct {
        FINT *atm;
        FINT *bas;
        double *env;
        FINT *shls;
        FINT natm;
        FINT nbas;

        FINT i_l;
        FINT j_l;
        FINT k_l;
        FINT l_l;
        FINT nfi;  // number of cartesion components
        FINT nfj;
        FINT nfk;
        FINT nfl;
        FINT nf;  // = nfi*nfj*nfk*nfl;
        FINT _padding;
        FINT x_ctr[4];

        FINT gbits;
        FINT ncomp_e1; // = 1 if spin free, = 4 when spin included, it
        FINT ncomp_e2; // corresponds to POSX,POSY,POSZ,POS1, see cint_const.h
        FINT ncomp_tensor; // e.g. = 3 for gradients

        /* values may diff based on the g0_2d4d algorithm */
        FINT li_ceil; // power of x, == i_l if nabla is involved, otherwise == i_l
        FINT lj_ceil;
        FINT lk_ceil;
        FINT ll_ceil;
        FINT g_stride_i; // nrys_roots * shift of (i++,k,l,j)
        FINT g_stride_k; // nrys_roots * shift of (i,k++,l,j)
        FINT g_stride_l; // nrys_roots * shift of (i,k,l++,j)
        FINT g_stride_j; // nrys_roots * shift of (i,k,l,j++)
        FINT nrys_roots;
        FINT g_size;  // ref to cint2e.c g = malloc(sizeof(double)*g_size)

        FINT g2d_ijmax;
        FINT g2d_klmax;
        double common_factor;
        double _padding1;
        double rirj[3]; // diff by sign in different g0_2d4d algorithm
        double rkrl[3];
        double *rx_in_rijrx;
        double *rx_in_rklrx;

        double *ri;
        double *rj;
        double *rk;
        double *rl;

        void (*f_g0_2e)();
        void (*f_g0_2d4d)();
        void (*f_gout)();

        /* values are assigned during calculation */
        FINT *idx;
        double ai;
        double aj;
        double ak;
        double al;
        double rij[3];
        double rijrx[3];
        double aij;
        double rkl[3];
        double rklrx[3];
        double akl;
} CINTEnvVars;
#endif

void CINTinit_int1e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTinit_int3c1e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                              FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);

void CINTg1e_index_xyz(FINT *idx, CINTEnvVars *envs);

void CINTg_ovlp(double *g, double ai, double aj, double fac, CINTEnvVars *envs);

void CINTg_nuc(double *g, double aij, double *rij,
               double *cr, double t2, double fac, CINTEnvVars *envs);

void CINTnabla1i_1e(double *f, double *g,
                    FINT li, FINT lj, FINT lk, CINTEnvVars *envs);

void CINTnabla1j_1e(double *f, double *g,
                    FINT li, FINT lj, FINT lk, CINTEnvVars *envs);

void CINTnabla1k_1e(double *f, double *g,
                    FINT li, FINT lj, FINT lk, CINTEnvVars *envs);

void CINTx1i_1e(double *f, double *g, double ri[3],
                FINT li, FINT lj, FINT lk, CINTEnvVars *envs);

void CINTx1j_1e(double *f, double *g, double rj[3],
                FINT li, FINT lj, FINT lk, CINTEnvVars *envs);

void CINTx1k_1e(double *f, double *g, double rk[3],
                FINT li, FINT lj, FINT lk, CINTEnvVars *envs);

void CINTprim_to_ctr(double *gc, FINT nf, double *gp,
                     FINT inc, FINT nprim,
                     FINT nctr, double *pcoeff);
void CINTprim_to_ctr_0(double *gc, FINT nf, double *gp,
                       FINT nprim, FINT nctr, double *coeff);
void CINTprim_to_ctr_1(double *gc, FINT nf, double *gp,
                       FINT nprim, FINT nctr, double *coeff);
void CINTprim_to_ctr_opt(double *gc, FINT nf, double *gp,
                         double *non0coeff, FINT *non0idx, FINT non0ctr);

double CINTcommon_fac_sp(FINT l);

#define G1E_D_I(f, g, li, lj, lk)   CINTnabla1i_1e(f, g, li, lj, lk, envs)
#define G1E_D_J(f, g, li, lj, lk)   CINTnabla1j_1e(f, g, li, lj, lk, envs)
#define G1E_D_K(f, g, li, lj, lk)   CINTnabla1k_1e(f, g, li, lj, lk, envs)
/* r-R_0, R_0 is (0,0,0) */
#define G1E_R0I(f, g, li, lj, lk)   CINTx1i_1e(f, g, envs->ri, li, lj, lk, envs)
#define G1E_R0J(f, g, li, lj, lk)   CINTx1j_1e(f, g, envs->rj, li, lj, lk, envs)
#define G1E_R0K(f, g, li, lj, lk)   CINTx1k_1e(f, g, envs->rk, li, lj, lk, envs)
/* r-R_C, R_C is common origin */
#define G1E_RCI(f, g, li, lj, lk)   CINTx1i_1e(f, g, dri, li, lj, lk, envs)
#define G1E_RCJ(f, g, li, lj, lk)   CINTx1j_1e(f, g, drj, li, lj, lk, envs)
#define G1E_RCK(f, g, li, lj, lk)   CINTx1k_1e(f, g, drk, li, lj, lk, envs)
/* origin from center of each basis
 * x1[ij]_1e(f, g, ng, li, lj, 0d0) */
#define G1E_R_I(f, g, li, lj, lk)   f = g + envs->g_stride_i
#define G1E_R_J(f, g, li, lj, lk)   f = g + envs->g_stride_j
#define G1E_R_K(f, g, li, lj, lk)   f = g + envs->g_stride_k
