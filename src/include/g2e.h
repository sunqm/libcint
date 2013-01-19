/*
 * File: g2e.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * Provide the intermediate variable g(nroots,i,j,k,l,[xyz])
 */

#include "g1e.h"

void g2e_index_xyz(int *idx, const int *ng, const int *shls,
                   const int *bas);

inline int g_size(const int *ng);

inline void extract_dim(const int *ng, int *di, int *dj, int *dk, int *dl);

inline void extract_fg_xyz(double *f, const double *g, const int *ng,
                           const double **pgx, const double **pgy, const double **pgz,
                           double **pfx, double **pfy, double **pfz);

void g0_2e(double *g, const int *ng,
           const double ai, const double aj,
           const double ak, const double al,
           const double *ri, const double *rj,
           const double *rk, const double *rl, const double fac);

void nabla1i_2e(double *f, const double *g, const int *ng,
                const int li, const int lj, const int lk, const int ll,
                const double ai);

void nabla1j_2e(double *f, const double *g, const int *ng,
                const int li, const int lj, const int lk, const int ll,
                const double aj);

void nabla1k_2e(double *f, const double *g, const int *ng,
                const int li, const int lj, const int lk, const int ll,
                const double ak);

void nabla1l_2e(double *f, const double *g, const int *ng,
                const int li, const int lj, const int lk, const int ll,
                const double al);

void x1i_2e(double *f, const double *g, const int *ng,
            const int li, const int lj, const int lk, const int ll,
            const double ri[3]);

void x1j_2e(double *f, const double *g, const int *ng,
            const int li, const int lj, const int lk, const int ll,
            const double rj[3]);

void x1k_2e(double *f, const double *g, const int *ng,
            const int li, const int lj, const int lk, const int ll,
            const double rk[3]);

void x1l_2e(double *f, const double *g, const int *ng,
            const int li, const int lj, const int lk, const int ll,
            const double rl[3]);


#define G2E_D_I(f, g, li, lj, lk, ll)   nabla1i_2e(f, g, ng, li, lj, lk, ll, ai)
#define G2E_D_J(f, g, li, lj, lk, ll)   nabla1j_2e(f, g, ng, li, lj, lk, ll, aj)
#define G2E_D_K(f, g, li, lj, lk, ll)   nabla1k_2e(f, g, ng, li, lj, lk, ll, ak)
#define G2E_D_L(f, g, li, lj, lk, ll)   nabla1l_2e(f, g, ng, li, lj, lk, ll, al)
/* r-R_0, R_0 is (0,0,0) */
#define G2E_R0I(f, g, li, lj, lk, ll)   x1i_2e(f, g, ng, li, lj, lk, ll, ri)
#define G2E_R0J(f, g, li, lj, lk, ll)   x1j_2e(f, g, ng, li, lj, lk, ll, rj)
#define G2E_R0K(f, g, li, lj, lk, ll)   x1k_2e(f, g, ng, li, lj, lk, ll, rk)
#define G2E_R0L(f, g, li, lj, lk, ll)   x1l_2e(f, g, ng, li, lj, lk, ll, rl)
/* r-R_C, R_C is common origin */
#define G2E_RCI(f, g, li, lj, lk, ll)   x1i_2e(f, g, ng, li, lj, lk, ll, dri)
#define G2E_RCJ(f, g, li, lj, lk, ll)   x1j_2e(f, g, ng, li, lj, lk, ll, drj)
#define G2E_RCK(f, g, li, lj, lk, ll)   x1k_2e(f, g, ng, li, lj, lk, ll, drk)
#define G2E_RCL(f, g, li, lj, lk, ll)   x1l_2e(f, g, ng, li, lj, lk, ll, drl)
/* origin from center of each basis
 * x1[ijkl]_2e(f, g, ng, li, lj, lk, ll, 0d0) */
#define G2E_R_I(f, g, li, lj, lk, ll)   f = g + ng[RYS_ROOTS]
#define G2E_R_K(f, g, li, lj, lk, ll)   f = g + ng[RYS_ROOTS] * ng[0]
#define G2E_R_L(f, g, li, lj, lk, ll)   f = g + ng[RYS_ROOTS] * ng[0] * ng[1]
#define G2E_R_J(f, g, li, lj, lk, ll)   f = g + ng[RYS_ROOTS] * ng[0] * ng[1] * ng[2]
