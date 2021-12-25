/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include "cint.h"

void CINTinit_int1e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTinit_int3c1e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                              FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);

void CINTg1e_index_xyz(FINT *idx, CINTEnvVars *envs);

FINT CINTg1e_ovlp(double *g, CINTEnvVars *envs);

FINT CINTg1e_nuc(double *g, CINTEnvVars *envs, FINT nuc_id);

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

double CINTcommon_fac_sp(FINT l);

void CINTprim_to_ctr_0(double *gc, double *gp, double *coeff, size_t nf,
                       FINT nprim, FINT nctr, FINT non0ctr, FINT *sortedidx);
void CINTprim_to_ctr_1(double *gc, double *gp, double *coeff, size_t nf,
                       FINT nprim, FINT nctr, FINT non0ctr, FINT *sortedidx);

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
