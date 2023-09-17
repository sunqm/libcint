/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 3-center 2-electron integrals
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "optimizer.h"
#include "g2e.h"
#include "cint2e.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"

#define gctrg   gout
#define gctrm   gctr
#define mempty  empty
#define m_ctr   n_comp
#define ALIAS_ADDR_IF_EQUAL(x, y) \
        if (y##_ctr == 1) { \
                gctr##x = gctr##y; \
                x##empty = y##empty; \
        } else { \
                gctr##x = g1; \
                g1 += len##x; \
        }

#define PRIM2CTR0(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } \
        } \
        *ctrsymb##empty = 0

#define TRANSPOSE(a) \
        if (*empty) { \
                CINTdmat_transpose(gctr, a, nf*nc, n_comp); \
        } else { \
                CINTdplus_transpose(gctr, a, nf*nc, n_comp); \
        } \
        *empty = 0;


FINT CINT3c2e_loop_nopt(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
{
        FINT *shls  = envs->shls;
        FINT *bas = envs->bas;
        double *env = envs->env;
        FINT i_sh = shls[0];
        FINT j_sh = shls[1];
        FINT k_sh = shls[2];
        FINT i_ctr = envs->x_ctr[0];
        FINT j_ctr = envs->x_ctr[1];
        FINT k_ctr = envs->x_ctr[2];
        FINT i_prim = bas(NPRIM_OF, i_sh);
        FINT j_prim = bas(NPRIM_OF, j_sh);
        FINT k_prim = bas(NPRIM_OF, k_sh);
        //double *ri = envs->ri;
        //double *rj = envs->rj;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);

        double expcutoff = envs->expcutoff;
        double rr_ij = SQUARE(envs->rirj);
        double *log_maxci, *log_maxcj;
        PairData *pdata_base, *pdata_ij;
        MALLOC_INSTACK(log_maxci, i_prim+j_prim);
        MALLOC_INSTACK(pdata_base, i_prim*j_prim);
        log_maxcj = log_maxci + i_prim;
        CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
        CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
        if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                             log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                             i_prim, j_prim, rr_ij, expcutoff, env)) {
                return 0;
        }

        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        size_t nf = envs->nf;
        double fac1i, fac1j, fac1k;
        FINT ip, jp, kp;
        FINT _empty[4] = {1, 1, 1, 1};
        FINT *iempty = _empty + 0;
        FINT *jempty = _empty + 1;
        FINT *kempty = _empty + 2;
        FINT *gempty = _empty + 3;
        /* COMMON_ENVS_AND_DECLARE end */

        double expij, cutoff;
        double *rij;
        double *rkl = envs->rk;
        double omega = env[PTR_RANGE_OMEGA];
        if (omega < 0 && envs->rys_order > 1) {
                double r_guess = 8.;
                double omega2 = omega * omega;
                int lij = envs->li_ceil + envs->lj_ceil;
                if (lij > 0) {
                        double dist_ij = sqrt(rr_ij);
                        double aij = ai[i_prim-1] + aj[j_prim-1];
                        double theta = omega2 / (omega2 + aij);
                        expcutoff += lij * approx_log(
                                (dist_ij+theta*r_guess+1.)/(dist_ij+1.));
                }
                if (envs->lk_ceil > 0) {
                        double theta = omega2 / (omega2 + ak[k_prim-1]);
                        expcutoff += envs->lk_ceil * approx_log(theta*r_guess+1.);
                }
        }

        FINT *idx;
        MALLOC_INSTACK(idx, nf * 3);
        CINTg2e_index_xyz(idx, envs);

        FINT *non0ctri, *non0ctrj, *non0ctrk;
        FINT *non0idxi, *non0idxj, *non0idxk;
        MALLOC_INSTACK(non0ctri, i_prim+j_prim+k_prim+i_prim*i_ctr+j_prim*j_ctr+k_prim*k_ctr);
        non0ctrj = non0ctri + i_prim;
        non0ctrk = non0ctrj + j_prim;
        non0idxi = non0ctrk + k_prim;
        non0idxj = non0idxi + i_prim*i_ctr;
        non0idxk = non0idxj + j_prim*j_ctr;
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);

        FINT nc = i_ctr * j_ctr * k_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        size_t lenk = nf * nc * n_comp; // gctrk
        size_t lenj = nf * i_ctr * j_ctr * n_comp; // gctrj
        size_t leni = nf * i_ctr * n_comp; // gctri
        size_t len0 = nf * n_comp; // gout
        size_t len = leng + lenk + lenj + leni + len0;
        double *g;
        MALLOC_INSTACK(g, len);  // must be allocated last in this function
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk;

        ALIAS_ADDR_IF_EQUAL(k, m);
        ALIAS_ADDR_IF_EQUAL(j, k);
        ALIAS_ADDR_IF_EQUAL(i, j);
        ALIAS_ADDR_IF_EQUAL(g, i);

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak[0] = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *jempty = 1;
                }

                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj[0] = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                if (pdata_ij->cceij > expcutoff) {
                                        goto i_contracted;
                                }
                                envs->ai[0] = ai[ip];
                                expij = pdata_ij->eij;
                                rij = pdata_ij->rij;
                                cutoff = expcutoff - pdata_ij->cceij;
                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*expij;
                                } else {
                                        fac1i = fac1j*expij;
                                }
                                envs->fac[0] = fac1i;
                                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs)) {
                                        (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                        PRIM2CTR0(i, gout, len0);
                                }
i_contracted: ;
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR0(j, gctri, leni);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR0(k, gctrj, lenj);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                TRANSPOSE(gctrk);
        }
        return !*empty;
}


#define COMMON_ENVS_AND_DECLARE \
        FINT *shls = envs->shls; \
        FINT *bas = envs->bas; \
        double *env = envs->env; \
        FINT i_sh = shls[0]; \
        FINT j_sh = shls[1]; \
        CINTOpt *opt = envs->opt; \
        if (opt->pairdata != NULL && \
            opt->pairdata[i_sh*opt->nbas+j_sh] == NOVALUE) { \
                return 0; \
        } \
        FINT k_sh = shls[2]; \
        FINT i_ctr = envs->x_ctr[0]; \
        FINT j_ctr = envs->x_ctr[1]; \
        FINT k_ctr = envs->x_ctr[2]; \
        FINT i_prim = bas(NPRIM_OF, i_sh); \
        FINT j_prim = bas(NPRIM_OF, j_sh); \
        FINT k_prim = bas(NPRIM_OF, k_sh); \
        double *ai = env + bas(PTR_EXP, i_sh); \
        double *aj = env + bas(PTR_EXP, j_sh); \
        double *ak = env + bas(PTR_EXP, k_sh); \
        double *ci = env + bas(PTR_COEFF, i_sh); \
        double *cj = env + bas(PTR_COEFF, j_sh); \
        double *ck = env + bas(PTR_COEFF, k_sh); \
        double expcutoff = envs->expcutoff; \
        double rr_ij = SQUARE(envs->rirj); \
        PairData *pdata_base, *pdata_ij; \
        if (opt->pairdata != NULL) { \
                pdata_base = opt->pairdata[i_sh*opt->nbas+j_sh]; \
        } else { \
                double *log_maxci = opt->log_max_coeff[i_sh]; \
                double *log_maxcj = opt->log_max_coeff[j_sh]; \
                MALLOC_INSTACK(pdata_base, i_prim*j_prim); \
                if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj, \
                                     log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil, \
                                     i_prim, j_prim, rr_ij, expcutoff, env)) { \
                        return 0; \
                } \
        } \
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor; \
        size_t nf = envs->nf; \
        double fac1i, fac1j, fac1k; \
        FINT ip, jp, kp; \
        FINT _empty[4] = {1, 1, 1, 1}; \
        FINT *iempty = _empty + 0; \
        FINT *jempty = _empty + 1; \
        FINT *kempty = _empty + 2; \
        FINT *gempty = _empty + 3; \
        FINT *non0ctri = opt->non0ctr[i_sh]; \
        FINT *non0ctrj = opt->non0ctr[j_sh]; \
        FINT *non0idxi = opt->sortedidx[i_sh]; \
        FINT *non0idxj = opt->sortedidx[j_sh]; \
        FINT *non0ctrk, *non0idxk; \
        MALLOC_INSTACK(non0ctrk, k_prim+k_prim*k_ctr); \
        non0idxk = non0ctrk + k_prim; \
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr); \
        double expij, cutoff; \
        double *rij; \
        double *rkl = envs->rkl; \
        FINT *idx = opt->index_xyz_array[envs->i_l*LMAX1*LMAX1 \
                                        +envs->j_l*LMAX1+envs->k_l]; \
        if (idx == NULL) { \
                MALLOC_INSTACK(idx, nf * 3); \
                CINTg2e_index_xyz(idx, envs); \
        }

#define ADJUST_CUTOFF      \
        double omega = env[PTR_RANGE_OMEGA]; \
        if (omega < 0 && envs->rys_order > 1) { \
                double r_guess = 8.; \
                double omega2 = omega * omega; \
                int lij = envs->li_ceil + envs->lj_ceil; \
                if (lij > 0) { \
                        double dist_ij = sqrt(rr_ij); \
                        double aij = ai[i_prim-1] + aj[j_prim-1]; \
                        double theta = omega2 / (omega2 + aij); \
                        expcutoff += lij * approx_log( \
                                (dist_ij+theta*r_guess+1.)/(dist_ij+1.)); \
                } \
                if (envs->lk_ceil > 0) { \
                        double theta = omega2 / (omega2 + ak[k_prim-1]); \
                        expcutoff += envs->lk_ceil * approx_log(theta*r_guess+1.); \
                } \
        }

#define SET_RIJ    \
        if (pdata_ij->cceij > expcutoff) { \
                goto i_contracted; \
        } \
        envs->ai[0] = ai[ip]; \
        expij = pdata_ij->eij; \
        rij = pdata_ij->rij;

#define PRIM2CTR(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } \
        } \
        *ctrsymb##empty = 0


// i_ctr = j_ctr = k_ctr = 1;
FINT CINT3c2e_111_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
{
        COMMON_ENVS_AND_DECLARE;
        ADJUST_CUTOFF;
        FINT nc = 1;
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        size_t len0 = envs->nf * n_comp;
        size_t len = leng + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *gout;
        if (n_comp == 1) {
                gout = gctr;
                gempty = empty;
        } else {
                gout = g + leng;
        }

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak[0] = ak[kp];
                fac1k = envs->common_factor * ck[kp];

                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj[0] = aj[jp];
                        fac1j = fac1k * cj[jp];
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                SET_RIJ;
                                cutoff = expcutoff - pdata_ij->cceij;
                                fac1i = fac1j*ci[ip]*expij;
                                envs->fac[0] = fac1i;
                                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs)) {
                                        (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                        *gempty = 0;
                                }
i_contracted: ;
                        } // end loop i_prim
                } // end loop j_prim
        } // end loop k_prim

        if (n_comp > 1 && !*gempty) {
                TRANSPOSE(gout);
        }
        return !*empty;
}

// i_ctr = n; j_ctr = k_ctr = 1;
FINT CINT3c2e_n11_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
{
        COMMON_ENVS_AND_DECLARE;
        ADJUST_CUTOFF;
        FINT nc = i_ctr;
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        size_t leni = nf * i_ctr * n_comp; // gctri
        size_t len0 = nf * n_comp; // gout
        size_t len = leng + leni + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *g1 = g + leng;
        double *gout, *gctri;
        ALIAS_ADDR_IF_EQUAL(i, m);
        gout = g1;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak[0] = ak[kp];
                fac1k = envs->common_factor * ck[kp];
                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj[0] = aj[jp];
                        fac1j = fac1k * cj[jp];
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                SET_RIJ;
                                cutoff = expcutoff - pdata_ij->cceij;
                                fac1i = fac1j*expij;
                                envs->fac[0] = fac1i;
                                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs)) {
                                        (*envs->f_gout)(gout, g, idx, envs, 1);
                                        PRIM2CTR(i, gout, len0);
                                }
i_contracted: ;
                        } // end loop i_prim
                } // end loop j_prim
        } // end loop k_prim

        if (n_comp > 1 && !*iempty) {
                TRANSPOSE(gctri);
        }
        return !*empty;
}

// j_ctr = n; i_ctr = k_ctr = 1;
FINT CINT3c2e_1n1_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
{
        COMMON_ENVS_AND_DECLARE;
        ADJUST_CUTOFF;
        FINT nc = j_ctr;
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        size_t lenj = nf * j_ctr * n_comp; // gctrj
        size_t len0 = nf * n_comp; // gout
        size_t len = leng + lenj + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *g1 = g + leng;
        double *gout, *gctrj;
        ALIAS_ADDR_IF_EQUAL(j, m);
        gout = g1;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak[0] = ak[kp];
                fac1k = envs->common_factor * ck[kp];
                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj[0] = aj[jp];
                        fac1j = fac1k;
                        *iempty = 1;
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                SET_RIJ;
                                cutoff = expcutoff - pdata_ij->cceij;
                                fac1i = fac1j*ci[ip]*expij;
                                envs->fac[0] = fac1i;
                                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs)) {
                                        (*envs->f_gout)(gout, g, idx, envs, *iempty);
                                        *iempty = 0;
                                }
i_contracted: ;
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR(j, gout, len0);
                        }
                } // end loop j_prim
        } // end loop k_prim

        if (n_comp > 1 && !*jempty) {
                TRANSPOSE(gctrj);
        }
        return !*empty;
}


FINT CINT3c2e_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
{
        COMMON_ENVS_AND_DECLARE;
        ADJUST_CUTOFF;
        FINT nc = i_ctr * j_ctr * k_ctr;
        // (irys,i,j,k,coord,0:1); +1 for nabla-r12
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        size_t lenk = nf * nc * n_comp; // gctrk
        size_t lenj = nf * i_ctr * j_ctr * n_comp; // gctrj
        size_t leni = nf * i_ctr * n_comp; // gctri
        size_t len0 = nf * n_comp; // gout
        size_t len = leng + lenk + lenj + leni + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk;

        ALIAS_ADDR_IF_EQUAL(k, m);
        ALIAS_ADDR_IF_EQUAL(j, k);
        ALIAS_ADDR_IF_EQUAL(i, j);
        ALIAS_ADDR_IF_EQUAL(g, i);

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak[0] = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *jempty = 1;
                }
                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj[0] = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                SET_RIJ;
                                cutoff = expcutoff - pdata_ij->cceij;
                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*expij;
                                } else {
                                        fac1i = fac1j*expij;
                                }
                                envs->fac[0] = fac1i;
                                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs)) {
                                        (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                        PRIM2CTR(i, gout, len0);
                                }
i_contracted: ;
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR(j, gctri, leni);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR0(k, gctrj, lenj);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                TRANSPOSE(gctrk);
        }
        return !*empty;
}

static FINT (*CINTf_3c2e_loop[8])(double *, CINTEnvVars *, double *, FINT *) = {
        CINT3c2e_loop,
        CINT3c2e_loop,
        CINT3c2e_loop,
        CINT3c2e_n11_loop,
        CINT3c2e_loop,
        CINT3c2e_1n1_loop,
        CINT3c2e_loop,
        CINT3c2e_111_loop,
};

#define PAIRDATA_NON0IDX_SIZE(ps) \
                FINT *bas = envs->bas; \
                FINT *shls  = envs->shls; \
                FINT i_prim = bas(NPRIM_OF, shls[0]); \
                FINT j_prim = bas(NPRIM_OF, shls[1]); \
                FINT k_prim = bas(NPRIM_OF, shls[2]); \
                FINT ps = (i_prim*j_prim * 5 \
                           + i_prim * x_ctr[0] \
                           + j_prim * x_ctr[1] \
                           + k_prim * x_ctr[2] \
                           +(i_prim+j_prim)*2 + k_prim + envs->nf*3 + 16);

CACHE_SIZE_T CINT3c2e_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                         double *cache, void (*f_e1_c2s)(), FINT is_ssc)
{
        FINT *x_ctr = envs->x_ctr;
        size_t nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        if (out == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                CACHE_SIZE_T leng = envs->g_size*3*((1<<envs->gbits)+1);
                CACHE_SIZE_T len0 = envs->nf*n_comp;
                CACHE_SIZE_T cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                      nc*n_comp+envs->nf*3);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1);
                size_t len0 = envs->nf*n_comp;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                      nc*n_comp+envs->nf*3);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        FINT n;
        FINT empty = 1;
        if (opt != NULL) {
                envs->opt = opt;
                n = ((envs->x_ctr[0]==1) << 2) + ((envs->x_ctr[1]==1) << 1) + (envs->x_ctr[2]==1);
                CINTf_3c2e_loop[n](gctr, envs, cache, &empty);
        } else {
                CINT3c2e_loop_nopt(gctr, envs, cache, &empty);
        }

        FINT counts[4];
        if (f_e1_c2s == &c2s_sph_3c2e1) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->j_l*2+1) * x_ctr[1];
                if (is_ssc) {
                        counts[2] = envs->nfk * x_ctr[2];
                } else {
                        counts[2] = (envs->k_l*2+1) * x_ctr[2];
                }
        } else {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfj * x_ctr[1];
                counts[2] = envs->nfk * x_ctr[2];
        }
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1] * dims[2];
        if (!empty) {
                for (n = 0; n < n_comp; n++) {
                        (*f_e1_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return !empty;
}
CACHE_SIZE_T CINT3c2e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)(), FINT is_ssc)
{
        FINT *x_ctr = envs->x_ctr;
        FINT counts[4];
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        if (is_ssc) {
                counts[2] = envs->nfk * x_ctr[2];
        } else {
                counts[2] = (envs->k_l*2+1) * x_ctr[2];
        }
        counts[3] = 1;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1);
                size_t len0 = envs->nf*n_comp;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                      nc*n_comp + envs->nf*14*OF_CMPLX);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1);
                size_t len0 = envs->nf*n_comp;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                        nc*n_comp + envs->nf*14*OF_CMPLX);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        FINT n;
        FINT empty = 1;
        if (opt != NULL) {
                envs->opt = opt;
                n = ((envs->x_ctr[0]==1) << 2) + ((envs->x_ctr[1]==1) << 1) + (envs->x_ctr[2]==1);
                CINTf_3c2e_loop[n](gctr, envs, cache, &empty);
        } else {
                CINT3c2e_loop_nopt(gctr, envs, cache, &empty);
        }

        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1] * dims[2];
        if (!empty) {
                for (n = 0; n < envs->ncomp_e2 * envs->ncomp_tensor; n++) {
                        (*f_e1_c2s)(out+nout*n, gctr, dims, envs, cache);
                        gctr += nc * envs->ncomp_e1;
                }
        } else {
                for (n = 0; n < envs->ncomp_e2 * envs->ncomp_tensor; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return !empty;
}


CACHE_SIZE_T int3c2e_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_drv(out, dims, &envs, opt, cache, &c2s_sph_3c2e1, 0);
}
void int3c2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_3c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

CACHE_SIZE_T int3c2e_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                 FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_drv(out, dims, &envs, opt, cache, &c2s_cart_3c2e1, 0);
}

CACHE_SIZE_T int3c2e_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                   FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1, 0);
}

CACHE_SIZE_T int3c2e_sph_ssc(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                    FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_drv(out, dims, &envs, opt, cache, &c2s_sph_3c2e1_ssc, 1);
}
CACHE_SIZE_T int3c2e_spinor_ssc(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1_ssc, 1);
}

void CINTgout2e_int3c2e_spsp1(double *g,
double *gout, FINT *idx, CINTEnvVars *envs, FINT gout_empty);
CACHE_SIZE_T int3c2e_spsp1_spinor_ssc(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                             FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {1, 1, 0, 0, 2, 4, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int3c2e_spsp1;
        return CINT3c2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_si_3c2e1_ssc, 1);
}
void int3c2e_ssc_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env)
{
        int3c2e_optimizer(opt, atm, natm, bas, nbas, env);
}



ALL_CINT(int3c2e)
ALL_CINT_FORTRAN_(int3c2e)

