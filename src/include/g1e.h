/*
 * File: g1e.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#if !defined HAVE_DEFINED_CINTENVVARS_H
#define HAVE_DEFINED_CINTENVVARS_H
// ref to CINTinit_int1e_EnvVars, CINTinit_int2e_EnvVars
typedef struct {
        int natm;
        int nbas;
        const int *atm;
        const int *bas;
        const double *env;
        const unsigned int *shls;

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
        const double *ri;
        const double *rj;
        const double *rk;
        const double *rl;
        double common_factor;

        unsigned int gbits;
        unsigned int ncomp_e1; // = 1 if spin free, = 4 when spin included, it
        unsigned int ncomp_e2; // corresponds to POSX,POSY,POSZ,POS1, see cint_const.h
        unsigned int ncomp_tensor; // e.g. = 3 for gradients

        /* values may diff based on the g0_2d4d algorithm */
        unsigned int nrys_roots;
        unsigned int li_ceil; // power of x, == i_l if nabla is involved, otherwise == i_l
        unsigned int lj_ceil;
        unsigned int lk_ceil;
        unsigned int ll_ceil;
        unsigned int g_stride_i; // nrys_roots * shift of (i++,k,l,j)
        unsigned int g_stride_k; // nrys_roots * shift of (i,k++,l,j)
        unsigned int g_stride_l; // nrys_roots * shift of (i,k,l++,j)
        unsigned int g_stride_j; // nrys_roots * shift of (i,k,l,j++)
        unsigned int g_size;  // ref to cint2e.c g = malloc(sizeof(double)*g_size)

        unsigned int g2d_ijmax;
        unsigned int g2d_klmax;
        const double *rx_in_rijrx;
        const double *rx_in_rklrx;
        double rirj[3]; // diff by an sign in different g0_2d4d algorithm
        double rkrl[3];

        void (*f_g0_2d4d)();

        /* */
        void (*f_gout)();

        /* values are assigned during calculation */
        unsigned int *idx;
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

void CINTg1e_index_xyz(unsigned int idx[], const unsigned int *ng,
                       const unsigned int shls[], const int *bas);

void CINTg_ovlp(double *g, const unsigned int *ng,
                const double ai, const double aj,
                const double *ri, const double *rj, const double fac);

void CINTg_nuc(double *g, const unsigned int *ng,
               const double aij, const double *rij,
               const double *ri, const double *rj,
               const double *cr, const double t2, const double fac);

void CINTnabla1i_1e(double *f, const double *g, const unsigned int *ng,
                    const int li, const int lj,
                    const double ai);

void CINTnabla1j_1e(double *f, const double *g, const unsigned int *ng,
                    const int li, const int lj,
                    const double aj);

void CINTx1i_1e(double *f, const double *g, const unsigned int *ng,
                const int li, const int lj,
                const double ri[3]);

void CINTx1j_1e(double *f, const double *g, const unsigned int *ng,
                const int li, const int lj,
                const double rj[3]);

void CINTprim_to_ctr(double *gc, const unsigned int nf, const double *gp,
                     const unsigned int inc, const unsigned int nprim,
                     const unsigned int nctr, const double *pcoeff);

double CINTcommon_fac_sp(unsigned int l);

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
