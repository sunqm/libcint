/*
 * File: g1e.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#if !defined HAVE_DEFINED_CINTENVVARS_H
#define HAVE_DEFINED_CINTENVVARS_H
// ref to CINTinit_int1e_EnvVars, CINTinit_int2e_EnvVars
typedef struct {
        const int *atm;
        const int *bas;
        const double *env;
        const int *shls;
        int natm;
        int nbas;

        int i_l;
        int j_l;
        int k_l;
        int l_l;
        int i_prim;
        int j_prim;
        int k_prim;
        int l_prim;
        int i_ctr;
        int j_ctr;
        int k_ctr;
        int l_ctr;
        int nfi;  // number of cartesion components
        int nfj;
        int nfk;
        int nfl;
        int nf;  // = nfi*nfj*nfk*nfl;
        int _padding1;
        const double *ri;
        const double *rj;
        const double *rk;
        const double *rl;
        double common_factor;

        int gbits;
        int ncomp_e1; // = 1 if spin free, = 4 when spin included, it
        int ncomp_e2; // corresponds to POSX,POSY,POSZ,POS1, see cint_const.h
        int ncomp_tensor; // e.g. = 3 for gradients

        /* values may diff based on the g0_2d4d algorithm */
        int li_ceil; // power of x, == i_l if nabla is involved, otherwise == i_l
        int lj_ceil;
        int lk_ceil;
        int ll_ceil;
        int g_stride_i; // nrys_roots * shift of (i++,k,l,j)
        int g_stride_k; // nrys_roots * shift of (i,k++,l,j)
        int g_stride_l; // nrys_roots * shift of (i,k,l++,j)
        int g_stride_j; // nrys_roots * shift of (i,k,l,j++)
        int nrys_roots;
        int g_size;  // ref to cint2e.c g = malloc(sizeof(double)*g_size)

        int g2d_ijmax;
        int g2d_klmax;
        const double *rx_in_rijrx;
        const double *rx_in_rklrx;
        double rirj[3]; // diff by an sign in different g0_2d4d algorithm
        double rkrl[3];

        void (*f_g0_2d4d)();

        /* */
        void (*f_gout)();

        /* values are assigned during calculation */
        int *idx;
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

void CINTg1e_index_xyz(int idx[], const int *ng,
                       const int shls[], const int *bas);

void CINTg_ovlp(double *g, const int *ng,
                const double ai, const double aj,
                const double *ri, const double *rj, const double fac);

void CINTg_nuc(double *g, const int *ng,
               const double aij, const double *rij,
               const double *ri, const double *rj,
               const double *cr, const double t2, const double fac);

void CINTnabla1i_1e(double *f, const double *g, const int *ng,
                    const int li, const int lj,
                    const double ai);

void CINTnabla1j_1e(double *f, const double *g, const int *ng,
                    const int li, const int lj,
                    const double aj);

void CINTx1i_1e(double *f, const double *g, const int *ng,
                const int li, const int lj,
                const double ri[3]);

void CINTx1j_1e(double *f, const double *g, const int *ng,
                const int li, const int lj,
                const double rj[3]);

void CINTprim_to_ctr(double *gc, const int nf, const double *gp,
                     const int inc, const int nprim,
                     const int nctr, const double *pcoeff);

double CINTcommon_fac_sp(int l);

#define G1E_D_I(f, g, li, lj)   CINTnabla1i_1e(f, g, ng, li, lj, ai)
#define G1E_D_J(f, g, li, lj)   CINTnabla1j_1e(f, g, ng, li, lj, aj)
/* r-R_0, R_0 is (0,0,0) */
#define G1E_R0I(f, g, li, lj)   CINTx1i_1e(f, g, ng, li, lj, ri)
#define G1E_R0J(f, g, li, lj)   CINTx1j_1e(f, g, ng, li, lj, rj)
/* r-R_C, R_C is common origin */
#define G1E_RCI(f, g, li, lj)   CINTx1i_1e(f, g, ng, li, lj, dri)
#define G1E_RCJ(f, g, li, lj)   CINTx1j_1e(f, g, ng, li, lj, drj)
/* origin from center of each basis
 * x1[ij]_1e(f, g, ng, li, lj, 0d0) */
#define G1E_R_I(f, g, li, lj)   f = g + 1
#define G1E_R_J(f, g, li, lj)   f = g + ng[0]
