/*
 * File: g1e.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 */

void g1e_index_xyz(int idx[], const int *ng, const int shls[],
                   const int *bas);

void g_ovlp(double *g, const int *ng,
            const double ai, const double aj,
            const double *ri, const double *rj, const double fac);

void g_nuc(double *g, const int *ng,
           const double aij, const double *rij,
           const double *ri, const double *rj,
           const double *cr, const double t2, const double fac);

void nabla1i_1e(double *f, const double *g, const int *ng,
                const int li, const int lj,
                const double ai);

void nabla1j_1e(double *f, const double *g, const int *ng,
                const int li, const int lj,
                const double aj);

void x1i_1e(double *f, const double *g, const int *ng,
            const int li, const int lj,
            const double ri[3]);

void x1j_1e(double *f, const double *g, const int *ng,
            const int li, const int lj,
            const double rj[3]);

void prim_to_ctr(double *gc, const int nf, const double *gp, const int inc,
                 const int shl, const int i_pgto,
                 const int *bas, const double *env);

#define G1E_D_I(f, g, li, lj)   nabla1i_1e(f, g, ng, li, lj, ai)
#define G1E_D_J(f, g, li, lj)   nabla1j_1e(f, g, ng, li, lj, aj)
/* r-R_0, R_0 is (0,0,0) */
#define G1E_R0I(f, g, li, lj)   x1i_1e(f, g, ng, li, lj, ri)
#define G1E_R0J(f, g, li, lj)   x1j_1e(f, g, ng, li, lj, rj)
/* r-R_C, R_C is common origin */
#define G1E_RCI(f, g, li, lj)   x1i_1e(f, g, ng, li, lj, dri)
#define G1E_RCJ(f, g, li, lj)   x1j_1e(f, g, ng, li, lj, drj)
/* origin from center of each basis
 * x1[ij]_1e(f, g, ng, li, lj, 0d0) */
#define G1E_R_I(f, g, li, lj)   f = g + 1
#define G1E_R_J(f, g, li, lj)   f = g + ng[0]
