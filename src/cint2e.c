/*
 * File: cint2e.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO integrals
 */

#include <stdlib.h>
#include <math.h>
#include "cint_bas.h"
#include "g2e.h"
#include "cint2e.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"

#define MIN(X,Y) (X)<(Y)?(X):(Y)

#define SWITCH_CTR(parentsymb, ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr == 1) {\
                ctrsymb##empty = 0; \
        } else { \
                if (ctrsymb##empty && \
                    (parentsymb##empty || (parentsymb##_ctr>1))) { \
                        prim_to_ctr_0(gctr##ctrsymb, ngp, gp, ctrsymb##_prim, \
                                      ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                        ctrsymb##empty = 0; \
                } else { \
                        prim_to_ctr_1(gctr##ctrsymb, ngp, gp, ctrsymb##_prim, \
                                      ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } \
        }

void prim_to_ctr_0(double *gc, const unsigned int nf, const double *gp,
                   const unsigned int nprim, const unsigned int nctr,
                   const double *coeff)
{
        unsigned int n, i;

        if (nctr > 1) {
                double *pgc = gc + nf;
                const double *pcoeff = coeff + nprim;
                for (i = 1; i < nctr; i++) {
                        if (*pcoeff != 0) {
                                for (n = 0; n < nf; n++) {
                                        pgc[n] = *pcoeff * gp[n];
                                }
                        }
                        pcoeff += nprim;
                        pgc += nf;
                }
        }

        if (*coeff != 0) {
                for (n = 0; n < nf; n++) {
                        gc[n] = *coeff * gp[n];
                }
        }
}

void prim_to_ctr_1(double *gc, const unsigned int nf, const double *gp,
                   const unsigned int nprim, const unsigned int nctr,
                   const double *coeff)
{
        unsigned int n, i;

        for (i = 0; i < nctr; i++) {
                if (*coeff != 0) {
                        for (n = 0; n < nf; n++) {
                                gc[n] += *coeff * gp[n];
                        }
                }
                coeff += nprim;
                gc += nf;
        }
}


// i_ctr = j_ctr = k_ctr = l_ctr = 1;
int cint2e_1111_loop(double *gctr, const double fac,
                     void (*const f_gout)(), CintEnvVars *envs)
{
        const unsigned int *shls  = envs->shls;
        const int *bas   = envs->bas;
        const double *env = envs->env;

        const unsigned int *ng    = envs->ng;
        const unsigned int i_prim = envs->i_prim;
        const unsigned int j_prim = envs->j_prim;
        const unsigned int k_prim = envs->k_prim;
        const unsigned int l_prim = envs->l_prim;
        const unsigned int nf     = envs->nf;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *rk = envs->rk;
        const double *rl = envs->rl;
        const unsigned int i_sh = shls[0];
        const unsigned int j_sh = shls[1];
        const unsigned int k_sh = shls[2];
        const unsigned int l_sh = shls[3];
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        const double *ak = env + bas(PTR_EXP, k_sh);
        const double *al = env + bas(PTR_EXP, l_sh);
        const double *ci = env + bas(PTR_COEFF, i_sh);
        const double *cj = env + bas(PTR_COEFF, j_sh);
        const double *ck = env + bas(PTR_COEFF, k_sh);
        const double *cl = env + bas(PTR_COEFF, l_sh);
        const unsigned int n_comp = ng[POS_E1] * ng[POS_E2] * ng[TENSOR];

        unsigned int ip, jp, kp, lp, n;
        double dist_ij, eij = 0;
        double dist_kl, ekl = 0;
        RijSets rij, rkl;

        rij.r1r2[0] = rj[0] - ri[0];
        rij.r1r2[1] = rj[1] - ri[1];
        rij.r1r2[2] = rj[2] - ri[2];
        dist_ij = square_dist(ri, rj);
        rkl.r1r2[0] = rl[0] - rk[0];
        rkl.r1r2[1] = rl[1] - rk[1];
        rkl.r1r2[2] = rl[2] - rk[2];
        dist_kl = square_dist(rk, rl);

        unsigned int *const idx = malloc(sizeof(unsigned int) * nf * 3);
        g2e_index_xyz(idx, envs);
        int empty = 1;
        double fac1i, fac1j, fac1k;
        int len = envs->g_size * 3 * ((1<<ng[GSHIFT])+1) // (irys,i,j,k,l,coord,0:1);
                + nf * n_comp;
        double *const g = (double *)malloc(sizeof(double) * len);
        double *gout;
        if (n_comp == 1) {
                gout = gctr;
        } else {
                gout = g + envs->g_size * 3 * ((1<<ng[GSHIFT])+1);
        }

        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                for (kp = 0; kp < k_prim; kp++) {
                        envs->ak = ak[kp];
                        rkl.a12 = ak[kp] + al[lp];
                        ekl = dist_kl * ak[kp] * al[lp] / rkl.a12;
                        if (ekl > EXPCUTOFF) {
                                goto k_contracted;
                        }
                        rkl.r12[0] = (ak[kp]*rk[0] + al[lp]*rl[0]) / rkl.a12;
                        rkl.r12[1] = (ak[kp]*rk[1] + al[lp]*rl[1]) / rkl.a12;
                        rkl.r12[2] = (ak[kp]*rk[2] + al[lp]*rl[2]) / rkl.a12;
                        rkl.r12r2[0] = rkl.r12[0] - rl[0];
                        rkl.r12r2[1] = rkl.r12[1] - rl[1];
                        rkl.r12r2[2] = rkl.r12[2] - rl[2];
                        fac1k = (M_PI*M_PI*M_PI)*2/SQRTPI * fac * cl[lp] * ck[kp];

                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++) {
                                        envs->ai = ai[ip];
                                        rij.a12 = ai[ip] + aj[jp];
                                        eij = dist_ij * ai[ip] * aj[jp] / rij.a12;
                                        if (eij > EXPCUTOFF) {
                                                goto i_contracted;
                                        }
                                        rij.r12[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / rij.a12;
                                        rij.r12[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / rij.a12;
                                        rij.r12[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / rij.a12;
                                        fac1i = fac1j*ci[ip]*exp(-(eij+ekl));
                                        if (ng[2]+ng[3] == 2) { // nmax = 0, mmax = 0
                                                if (empty) {
                                                        gout[0] = g0_2e_ssss(&rij, &rkl, fac1i, envs);
                                                } else {
                                                        gout[0]+= g0_2e_ssss(&rij, &rkl, fac1i, envs);
                                                }
                                        } else {
                                                rij.r12r2[0] = rij.r12[0] - rj[0];
                                                rij.r12r2[1] = rij.r12[1] - rj[1];
                                                rij.r12r2[2] = rij.r12[2] - rj[2];
                                                g0_2e(g, &rij, &rkl, fac1i, envs);
                                                (*f_gout)(g, gout, idx, envs, empty);
                                        }
                                        empty = 0;
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim
        if (n_comp > 1 && !empty) {
                for (kp = 0; kp < n_comp; kp++) {
                        for (n = 0, ip = kp; n < nf; n++, ip+=n_comp) {
                                gctr[n] = gout[ip];
                        }
                        gctr += nf;
                }
        }
        free(idx);
        free(g);
        return !empty;
}


// note ng[:4] is ordered as [k,i,j,l]
// add fac*<2e-integrals> to gctr
int cint2e_loop(double *gctr, const double fac,
                void (*const f_gout)(), CintEnvVars *envs)
{
        const unsigned int *shls  = envs->shls;
        const int *bas   = envs->bas;
        const double *env = envs->env;

        const unsigned int *ng    = envs->ng;
        const unsigned int i_prim = envs->i_prim;
        const unsigned int j_prim = envs->j_prim;
        const unsigned int k_prim = envs->k_prim;
        const unsigned int l_prim = envs->l_prim;
        const unsigned int i_ctr  = envs->i_ctr;
        const unsigned int j_ctr  = envs->j_ctr;
        const unsigned int k_ctr  = envs->k_ctr;
        const unsigned int l_ctr  = envs->l_ctr;
        const unsigned int nf     = envs->nf;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *rk = envs->rk;
        const double *rl = envs->rl;
        const unsigned int i_sh = shls[0];
        const unsigned int j_sh = shls[1];
        const unsigned int k_sh = shls[2];
        const unsigned int l_sh = shls[3];
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        const double *ak = env + bas(PTR_EXP, k_sh);
        const double *al = env + bas(PTR_EXP, l_sh);
        const double *ci = env + bas(PTR_COEFF, i_sh);
        const double *cj = env + bas(PTR_COEFF, j_sh);
        const double *ck = env + bas(PTR_COEFF, k_sh);
        const double *cl = env + bas(PTR_COEFF, l_sh);
        const unsigned int nc = i_ctr * j_ctr * k_ctr * l_ctr;
        const unsigned int n_comp = ng[POS_E1] * ng[POS_E2] * ng[TENSOR];

        unsigned int ip, jp, kp, lp, n;
        double dist_ij, eij = 0;
        double dist_kl, ekl = 0;
        RijSets rij, rkl;

        rij.r1r2[0] = rj[0] - ri[0];
        rij.r1r2[1] = rj[1] - ri[1];
        rij.r1r2[2] = rj[2] - ri[2];
        dist_ij = square_dist(ri, rj);
        rkl.r1r2[0] = rl[0] - rk[0];
        rkl.r1r2[1] = rl[1] - rk[1];
        rkl.r1r2[2] = rl[2] - rk[2];
        dist_kl = square_dist(rk, rl);

        unsigned int *const idx = malloc(sizeof(unsigned int) * nf * 3);
        int iempty, jempty, kempty, lempty = 1;
        double fac1i, fac1j, fac1k, fac1l;
        int len = envs->g_size * 3 * ((1<<ng[GSHIFT])+1) // (irys,i,j,k,l,coord,0:1);
                + nf * nc * n_comp // gctrl
                + nf * i_ctr * j_ctr * k_ctr * n_comp // gctrk
                + nf * i_ctr * j_ctr * n_comp // gctrj
                + nf * i_ctr * n_comp // gctri
                + nf * n_comp; // gout
        double *const g = (double *)malloc(sizeof(double) * len);
        double *gout, *gctri, *gctrj, *gctrk, *gctrl;
        double *g1 = g + envs->g_size*3*((1<<ng[GSHIFT])+1);

        if (n_comp == 1) {
                gctrl = gctr;
        } else {
                gctrl = g1;
                g1 += nf * nc * n_comp;
        }
        if (l_ctr == 1) {
                gctrk = gctrl;
        } else {
                gctrk = g1;
                g1 += nf*i_ctr*j_ctr*k_ctr*n_comp;
        }
        if (k_ctr == 1) {
                gctrj = gctrk;
        } else {
                gctrj = g1;
                g1 += nf*i_ctr*j_ctr*n_comp;
        }
        if (j_ctr == 1) {
                gctri = gctrj;
        } else {
                gctri = g1;
                g1 += nf*i_ctr*n_comp;
        }
        if (i_ctr == 1) {
                gout = gctri;
        } else {
                gout = g1;
        }

        g2e_index_xyz(idx, envs);

        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                if (l_ctr == 1) {
                        fac1l = (M_PI*M_PI*M_PI)*2/SQRTPI * fac * cl[lp];
                        kempty = lempty; // same address for gctrk, gctrl
                } else {
                        fac1l = (M_PI*M_PI*M_PI)*2/SQRTPI * fac;
                        kempty = 1;
                }
                for (kp = 0; kp < k_prim; kp++) {
                        envs->ak = ak[kp];
                        rkl.a12 = ak[kp] + al[lp];
                        ekl = dist_kl * ak[kp] * al[lp] / rkl.a12;
                        if (ekl > EXPCUTOFF) {
                                goto k_contracted;
                        }
                        rkl.r12[0] = (ak[kp]*rk[0] + al[lp]*rl[0]) / rkl.a12;
                        rkl.r12[1] = (ak[kp]*rk[1] + al[lp]*rl[1]) / rkl.a12;
                        rkl.r12[2] = (ak[kp]*rk[2] + al[lp]*rl[2]) / rkl.a12;
                        rkl.r12r2[0] = rkl.r12[0] - rl[0];
                        rkl.r12r2[1] = rkl.r12[1] - rl[1];
                        rkl.r12r2[2] = rkl.r12[2] - rl[2];
                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                                jempty = kempty; // same address for gctrj, gctrk
                        } else {
                                fac1k = fac1l;
                                jempty = 1;
                        }

                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                if (j_ctr == 1) {
                                        fac1j = fac1k * cj[jp];
                                        iempty = jempty; // same address for gctri, gctrj
                                } else {
                                        fac1j = fac1k;
                                        iempty = 1;
                                }
                                for (ip = 0; ip < i_prim; ip++) {
                                        envs->ai = ai[ip];
                                        rij.a12 = ai[ip] + aj[jp];
                                        eij = dist_ij * ai[ip] * aj[jp] / rij.a12;
                                        if (eij > EXPCUTOFF) {
                                                goto i_contracted;
                                        }
                                        if (i_ctr == 1) {
                                                fac1i = fac1j*ci[ip]*exp(-(eij+ekl));
                                        } else {
                                                fac1i = fac1j*exp(-(eij+ekl));
                                        }
                                        rij.r12[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / rij.a12;
                                        rij.r12[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / rij.a12;
                                        rij.r12[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / rij.a12;

                                        if (ng[2]+ng[3] == 2) { // nmax = 0, mmax = 0
                                                if (i_ctr > 1 || iempty) {
                                                        *gout = g0_2e_ssss(&rij, &rkl, fac1i, envs);
                                                } else { // if same address for gctri and gout and
                                                         // it has been initialized.
                                                        *gout+= g0_2e_ssss(&rij, &rkl, fac1i, envs);
                                                }
                                        } else {
                                                rij.r12r2[0] = rij.r12[0] - rj[0];
                                                rij.r12r2[1] = rij.r12[1] - rj[1];
                                                rij.r12r2[2] = rij.r12[2] - rj[2];
                                                g0_2e(g, &rij, &rkl, fac1i, envs);
                                                (*f_gout)(g, gout, idx, envs,(i_ctr>1||iempty));
                                        }
                                        SWITCH_CTR(j, i, gout, nf*n_comp);
i_contracted: ;
                                } // end loop i_prim
                                if (!iempty) {
                                        SWITCH_CTR(k,j, gctri, nf*i_ctr*n_comp);
                                }
                        } // end loop j_prim
                        if (!jempty) {
                                SWITCH_CTR(l, k, gctrj, nf*i_ctr*j_ctr*n_comp);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!kempty) {
                        if (l_ctr > 1) {
                                n = nf*i_ctr*j_ctr*k_ctr*n_comp;
                                if (lempty) {
                                        prim_to_ctr_0(gctrl, n, gctrk, l_prim, l_ctr, cl+lp);
                                        lempty = 0;
                                } else {
                                        prim_to_ctr_1(gctrl, n, gctrk, l_prim, l_ctr, cl+lp);
                                }
                        } else {
                                lempty = 0;
                        }
                }
        } // end loop l_prim

        if (n_comp > 1 && !lempty) {
                for (kp = 0; kp < n_comp; kp++) {
                        for (n = 0, ip = kp; n < nf*nc; n++, ip+=n_comp) {
                                gctr[n] = gctrl[ip];
                        }
                        gctr += nf * nc;
                }
        }
        free(idx);
        free(g);
        return !lempty;
}


int cint2e_drv(double *opijkl, unsigned int *ng, const double fac,
               void (*const f_gout)(),
               void (*const f_e1_c2s)(), void (*const f_e2_c2s)(),
               const unsigned int *shls, const int *atm, const int natm,
               const int *bas, const int nbas, const double *env)
{
        ng[RYS_ROOTS] = (ng[2] + ng[3]) / 2;
        CintEnvVars envs;
        init_int2e_CintEnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);

        const unsigned int i_sh = shls[0];
        const unsigned int j_sh = shls[1];
        const unsigned int k_sh = shls[2];
        const unsigned int l_sh = shls[3];
        const unsigned int nc = envs.nf * envs.i_ctr * envs.j_ctr
                               * envs.k_ctr * envs.l_ctr * ng[POS_E1];
        double *const gctr = (double *)malloc(sizeof(double) * nc
                                              * ng[POS_E2] * ng[TENSOR]);
        double *pgctr = gctr;
        double *opij;
        int ip, jp, kp, lp, nop;
        int n, n1, n2;
        int has_value;

        if ((envs.i_ctr+envs.j_ctr+envs.k_ctr+envs.l_ctr) == 4) {
                has_value = cint2e_1111_loop(gctr, fac, f_gout, &envs);
        } else {
                has_value = cint2e_loop(gctr, fac, f_gout, &envs);
        }

        unsigned int (*fcgtos)();
        if (f_e1_c2s == &c2s_sph_2e1) {
                fcgtos = cgtos_spheric;
        } else if (f_e1_c2s == &c2s_cart_2e1) {
                fcgtos = cgtos_cart;
        } else {
                fcgtos = cgtos_spinor;
        }
        ip = (*fcgtos)(i_sh, bas);
        jp = (*fcgtos)(j_sh, bas);
        kp = (*fcgtos)(k_sh, bas);
        lp = (*fcgtos)(l_sh, bas);
        nop = ip * jp * kp * lp;

        if (f_e1_c2s == &c2s_sph_2e1 || f_e1_c2s == &c2s_cart_2e1) {
                if (has_value) {
                        for (n = 0; n < ng[TENSOR]; n++) {
                                (*f_e1_c2s)(opijkl, pgctr, shls, bas);
                                opijkl += nop;
                                pgctr += nc;
                        }
                } else {
                        dset0(nop * ng[TENSOR], opijkl);
                }
        } else {
                if (has_value) {
                        n1 = ip * envs.nfk * envs.k_ctr
                                * envs.nfl * envs.l_ctr * jp * OF_CMPLX;
                        opij = (double *)malloc(sizeof(double)*n1*ng[POS_E2]);
                        for (n = 0; n < ng[TENSOR]; n++) {
                                for (n2 = 0; n2 < ng[POS_E2]; n2++) {
                                        (*f_e1_c2s)(opij + n1 * n2, pgctr, shls, bas);
                                        pgctr += nc;
                                }
                                (*f_e2_c2s)(opijkl, opij, shls, bas);
                                opijkl += nop * OF_CMPLX;
                        }
                        free(opij);
                } else {
                        dset0(nop * OF_CMPLX * ng[TENSOR], opijkl);
                }
        }
        free(gctr);
        return has_value;
}


/*
 * <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
void gout2e(double *g, double *gout, const unsigned int *idx,
            const CintEnvVars *envs, int gout_empty)
{
        const unsigned int *ng = envs->ng;
        unsigned int nf = envs->nf;
        unsigned int i, ix, iy, iz, n;

        if (gout_empty) {
                switch (ng[RYS_ROOTS]) {
                        case 1:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix] * g[iy] * g[iz];
                                }
                                break;
                        case 2:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1];
                                }
                                break;
                        case 3:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2];
                                }
                                break;
                        case 4:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3];
                                }
                                break;
                        case 5:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4];
                                }
                                break;
                        case 6:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5];
                                }
                                break;
                        case 7:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6];
                                }
                                break;
                        default:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = 0;
                                        for (i = 0; i < ng[RYS_ROOTS]; i++)
                                                gout[n] += g[ix+i] * g[iy+i] * g[iz+i];
                                }
                                break;
                } // end switch nroots
        } else { // not flag_acc
                switch (ng[RYS_ROOTS]) {
                        case 1:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] += g[ix] * g[iy] * g[iz];
                                }
                                break;
                        case 2:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1];
                                }
                                break;
                        case 3:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2];
                                }
                                break;
                        case 4:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3];
                                }
                                break;
                        case 5:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4];
                                }
                                break;
                        case 6:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5];
                                }
                                break;
                        case 7:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6];
                                }
                                break;
                        default:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        for (i = 0; i < ng[RYS_ROOTS]; i++)
                                                gout[n] += g[ix+i] * g[iy+i] * g[iz+i];
                                }
                                break;
                } // end switch nroots
        }
}


int cint2e_sph(double *opijkl, const unsigned int *shls,
               const int *atm, const int natm,
               const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int k_sh = shls[2];
        const int l_sh = shls[3];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        const int k_l = bas(ANG_OF, k_sh);
        const int l_l = bas(ANG_OF, l_sh);
        unsigned int ng[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

        // maximum l supported: l=4
        ng[0] = i_l + 1;       // 0:i_l
        ng[1] = k_l + 1;       // 0:k_l
        ng[2] = l_l + ng[1];   // 0:mmax, mmax = lk + ll, 0:l_l
        ng[3] = j_l + ng[0];   // 0:nmax, nmax = li + lj, 0:j_l

        return cint2e_drv(opijkl, ng, 1, &gout2e,
                          &c2s_sph_2e1, NULL,
                          shls, atm, natm, bas, nbas, env);
}

int cint2e_cart(double *opijkl, const unsigned int *shls,
                const int *atm, const int natm,
                const int *bas, const int nbas, const double *env)
{
        const unsigned int i_sh = shls[0];
        const unsigned int j_sh = shls[1];
        const unsigned int k_sh = shls[2];
        const unsigned int l_sh = shls[3];
        const unsigned int i_l = bas(ANG_OF, i_sh);
        const unsigned int j_l = bas(ANG_OF, j_sh);
        const unsigned int k_l = bas(ANG_OF, k_sh);
        const unsigned int l_l = bas(ANG_OF, l_sh);
        unsigned int ng[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

        // maximum l supported: l=4
        ng[0] = i_l + 1;       // 0:i_l
        ng[1] = k_l + 1;       // 0:k_l
        ng[2] = l_l + ng[1];   // 0:mmax, mmax = lk + ll, 0:l_l
        ng[3] = j_l + ng[0];   // 0:nmax, nmax = li + lj, 0:j_l

        return cint2e_drv(opijkl, ng, 1, &gout2e,
                          &c2s_cart_2e1, NULL,
                          shls, atm, natm, bas, nbas, env);
}


/*
 * spinor <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
int cint2e(double *opijkl, const unsigned int *shls,
           const int *atm, const int natm,
           const int *bas, const int nbas, const double *env)
{
        const unsigned int i_sh = shls[0];
        const unsigned int j_sh = shls[1];
        const unsigned int k_sh = shls[2];
        const unsigned int l_sh = shls[3];
        const unsigned int i_l = bas(ANG_OF, i_sh);
        const unsigned int j_l = bas(ANG_OF, j_sh);
        const unsigned int k_l = bas(ANG_OF, k_sh);
        const unsigned int l_l = bas(ANG_OF, l_sh);
        unsigned int ng[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

        // maximum l supported: l=4
        ng[0] = i_l + 1;       // 0:i_l
        ng[1] = k_l + 1;       // 0:k_l
        ng[2] = l_l + ng[1];   // 0:mmax, mmax = lk + ll, 0:l_l
        ng[3] = j_l + ng[0];   // 0:nmax, nmax = li + lj, 0:j_l

        return cint2e_drv(opijkl, ng, 1, &gout2e,
                          &c2s_sf_2e1, &c2s_sf_2e2,
                          shls, atm, natm, bas, nbas, env);
}


/*
 * * * * * * * * * * * * * * * * * * * * *
 * c to fortran interface
 */

C2F_(cint2e_cart);
C2F_(cint2e_sph);
C2F_(cint2e);

