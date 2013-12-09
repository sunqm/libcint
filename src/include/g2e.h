/*
 * File: g2e.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * Provide the intermediate variable g(nroots,i,j,k,l,[xyz])
 */

#include "g1e.h"

void g2e_index_xyz(unsigned int *idx, const CintEnvVars *envs);

int init_int2e_CintEnvVars(CintEnvVars *envs, const unsigned int ng[],
                           const unsigned int *shls,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);

void g0_2e(double *g, const RijSets *rij, const RijSets *rkl,
           const double fac, const CintEnvVars *envs);

double g0_2e_ssss(const RijSets *rij, const RijSets *rkl,
                  const double fac, const CintEnvVars *envs);

void nabla1i_2e(double *f, const double *g,
                const unsigned int li, const unsigned int lj,
                const unsigned int lk, const unsigned int ll,
                const CintEnvVars *envs);

void nabla1j_2e(double *f, const double *g,
                const unsigned int li, const unsigned int lj,
                const unsigned int lk, const unsigned int ll,
                const CintEnvVars *envs);

void nabla1k_2e(double *f, const double *g,
                const unsigned int li, const unsigned int lj,
                const unsigned int lk, const unsigned int ll,
                const CintEnvVars *envs);

void nabla1l_2e(double *f, const double *g,
                const unsigned int li, const unsigned int lj,
                const unsigned int lk, const unsigned int ll,
                const CintEnvVars *envs);

void x1i_2e(double *f, const double *g,
            const unsigned int li, const unsigned int lj,
            const unsigned int lk, const unsigned int ll,
            const double *ri, const CintEnvVars *envs);

void x1j_2e(double *f, const double *g,
            const unsigned int li, const unsigned int lj,
            const unsigned int lk, const unsigned int ll,
            const double *rj, const CintEnvVars *envs);

void x1k_2e(double *f, const double *g,
            const unsigned int li, const unsigned int lj,
            const unsigned int lk, const unsigned int ll,
            const double *rk, const CintEnvVars *envs);

void x1l_2e(double *f, const double *g,
            const unsigned int li, const unsigned int lj,
            const unsigned int lk, const unsigned int ll,
            const double *rl, const CintEnvVars *envs);


#define G2E_D_I(f, g, li, lj, lk, ll)   nabla1i_2e(f, g, li, lj, lk, ll, envs)
#define G2E_D_J(f, g, li, lj, lk, ll)   nabla1j_2e(f, g, li, lj, lk, ll, envs)
#define G2E_D_K(f, g, li, lj, lk, ll)   nabla1k_2e(f, g, li, lj, lk, ll, envs)
#define G2E_D_L(f, g, li, lj, lk, ll)   nabla1l_2e(f, g, li, lj, lk, ll, envs)
/* r-R_0, R_0 is (0,0,0) */
#define G2E_R0I(f, g, li, lj, lk, ll)   x1i_2e(f, g, li, lj, lk, ll, ri, envs)
#define G2E_R0J(f, g, li, lj, lk, ll)   x1j_2e(f, g, li, lj, lk, ll, rj, envs)
#define G2E_R0K(f, g, li, lj, lk, ll)   x1k_2e(f, g, li, lj, lk, ll, rk, envs)
#define G2E_R0L(f, g, li, lj, lk, ll)   x1l_2e(f, g, li, lj, lk, ll, rl, envs)
/* r-R_C, R_C is common origin */
#define G2E_RCI(f, g, li, lj, lk, ll)   x1i_2e(f, g, li, lj, lk, ll, dri, envs)
#define G2E_RCJ(f, g, li, lj, lk, ll)   x1j_2e(f, g, li, lj, lk, ll, drj, envs)
#define G2E_RCK(f, g, li, lj, lk, ll)   x1k_2e(f, g, li, lj, lk, ll, drk, envs)
#define G2E_RCL(f, g, li, lj, lk, ll)   x1l_2e(f, g, li, lj, lk, ll, drl, envs)
/* origin from center of each basis
 * x1[ijkl]_2e(f, g, ng, li, lj, lk, ll, 0d0) */
#define G2E_R_I(f, g, li, lj, lk, ll)   f = g + envs->g_stride_i
#define G2E_R_K(f, g, li, lj, lk, ll)   f = g + envs->g_stride_k
#define G2E_R_L(f, g, li, lj, lk, ll)   f = g + envs->g_stride_l
#define G2E_R_J(f, g, li, lj, lk, ll)   f = g + envs->g_stride_j
