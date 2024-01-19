/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO integrals
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "cint_bas.h"
#include "g1e.h"
#include "g2e.h"
#include "optimizer.h"
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

#define TRANSPOSE(a) \
        if (*empty) { \
                CINTdmat_transpose(gctr, a, nf*nc, n_comp); \
                *empty = 0; \
        } else { \
                CINTdplus_transpose(gctr, a, nf*nc, n_comp); \
        } \

FINT CINT2e_loop_nopt(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
{
        /* COMMON_ENVS_AND_DECLARE */
        FINT *shls  = envs->shls;
        FINT *bas = envs->bas;
        double *env = envs->env;
        FINT i_sh = shls[0];
        FINT j_sh = shls[1];
        FINT k_sh = shls[2];
        FINT l_sh = shls[3];
        FINT i_ctr = envs->x_ctr[0];
        FINT j_ctr = envs->x_ctr[1];
        FINT k_ctr = envs->x_ctr[2];
        FINT l_ctr = envs->x_ctr[3];
        FINT i_prim = bas(NPRIM_OF, i_sh);
        FINT j_prim = bas(NPRIM_OF, j_sh);
        FINT k_prim = bas(NPRIM_OF, k_sh);
        FINT l_prim = bas(NPRIM_OF, l_sh);
        //double *ri = envs->ri;
        //double *rj = envs->rj;
        double *rk = envs->rk;
        double *rl = envs->rl;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *al = env + bas(PTR_EXP, l_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        double *cl = env + bas(PTR_COEFF, l_sh);
        double expcutoff = envs->expcutoff;
        double rr_ij = SQUARE(envs->rirj);
        double rr_kl = SQUARE(envs->rkrl);
        double *log_maxci, *log_maxcj, *log_maxck, *log_maxcl;
        PairData *pdata_base, *pdata_ij;
        MALLOC_INSTACK(log_maxci, i_prim+j_prim+k_prim+l_prim);
        MALLOC_INSTACK(pdata_base, i_prim*j_prim);
        log_maxcj = log_maxci + i_prim;
        log_maxck = log_maxcj + j_prim;
        log_maxcl = log_maxck + k_prim;
        CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
        CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
        if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                             log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                             i_prim, j_prim, rr_ij, expcutoff, env)) {
                return 0;
        }
        CINTOpt_log_max_pgto_coeff(log_maxck, ck, k_prim, k_ctr);
        CINTOpt_log_max_pgto_coeff(log_maxcl, cl, l_prim, l_ctr);

        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        size_t nf = envs->nf;
        double fac1i, fac1j, fac1k, fac1l;
        FINT ip, jp, kp, lp;
        FINT _empty[5] = {1, 1, 1, 1, 1};
        FINT *iempty = _empty + 0;
        FINT *jempty = _empty + 1;
        FINT *kempty = _empty + 2;
        FINT *lempty = _empty + 3;
        FINT *gempty = _empty + 4;
        /* COMMON_ENVS_AND_DECLARE end */

        int lkl = envs->lk_ceil + envs->ll_ceil;
        double akl, ekl, expijkl, ccekl, log_rr_kl, eijcutoff, cutoff;
        double rkl[3];
        double *rij;
        akl = ak[k_prim-1] + al[l_prim-1];
        log_rr_kl = 1.7 - 1.5 * approx_log(akl);
        double omega = env[PTR_RANGE_OMEGA];
        if (omega < 0) {
                // Normally the factor
                //    (aj*d/aij+theta*R)^li * (ai*d/aij+theta*R)^lj * pi^1.5/aij^{(li+lj+3)/2}
                // is a good approximation for polynomial parts in SR-ERIs.
                //    <~ (aj*d/aij+theta*R)^li * (ai*d/aij+theta*R)^lj * (pi/aij)^1.5
                //    <~ (d+theta*R)^li * (d+theta*R)^lj * (pi/aij)^1.5
                if (envs->rys_order > 1) {
                        double r_guess = 8.;
                        double omega2 = omega * omega;
                        int lij = envs->li_ceil + envs->lj_ceil;
                        if (lij > 0) {
                                double aij = ai[i_prim-1] + aj[j_prim-1];
                                double dist_ij = sqrt(rr_ij);
                                double theta = omega2 / (omega2 + aij);
                                expcutoff += lij * approx_log(
                                        (dist_ij+theta*r_guess+1.)/(dist_ij+1.));
                        }
                        if (lkl > 0) {
                                double theta = omega2 / (omega2 + akl);
                                log_rr_kl += lkl * approx_log(
                                        sqrt(rr_kl) + theta*r_guess + 1.);
                        }
                }
        } else {
                if (lkl > 0) {
                        log_rr_kl += lkl * approx_log(sqrt(rr_kl) + 1.);
                }
        }

        FINT *idx;
        MALLOC_INSTACK(idx, nf * 3);
        CINTg2e_index_xyz(idx, envs);

        FINT *non0ctri, *non0ctrj, *non0ctrk, *non0ctrl;
        FINT *non0idxi, *non0idxj, *non0idxk, *non0idxl;
        MALLOC_INSTACK(non0ctri, i_prim+j_prim+k_prim+l_prim+i_prim*i_ctr+j_prim*j_ctr+k_prim*k_ctr+l_prim*l_ctr);
        non0ctrj = non0ctri + i_prim;
        non0ctrk = non0ctrj + j_prim;
        non0ctrl = non0ctrk + k_prim;
        non0idxi = non0ctrl + l_prim;
        non0idxj = non0idxi + i_prim*i_ctr;
        non0idxk = non0idxj + j_prim*j_ctr;
        non0idxl = non0idxk + k_prim*k_ctr;
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
        CINTOpt_non0coeff_byshell(non0idxl, non0ctrl, cl, l_prim, l_ctr);

        FINT nc = i_ctr * j_ctr * k_ctr * l_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        size_t lenl = nf * nc * n_comp; // gctrl
        size_t lenk = nf * i_ctr * j_ctr * k_ctr * n_comp; // gctrk
        size_t lenj = nf * i_ctr * j_ctr * n_comp; // gctrj
        size_t leni = nf * i_ctr * n_comp; // gctri
        size_t len0 = nf * n_comp; // gout
        size_t len = leng + lenl + lenk + lenj + leni + len0;
        double *g;
        MALLOC_INSTACK(g, len);  // must be allocated last in this function
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk, *gctrl;
        ALIAS_ADDR_IF_EQUAL(l, m);
        ALIAS_ADDR_IF_EQUAL(k, l);
        ALIAS_ADDR_IF_EQUAL(j, k);
        ALIAS_ADDR_IF_EQUAL(i, j);
        ALIAS_ADDR_IF_EQUAL(g, i);

        for (lp = 0; lp < l_prim; lp++) {
                envs->al[0] = al[lp];
                if (l_ctr == 1) {
                        fac1l = envs->common_factor * cl[lp];
                } else {
                        fac1l = envs->common_factor;
                        *kempty = 1;
                }
                for (kp = 0; kp < k_prim; kp++) {
                        akl = ak[kp] + al[lp];
                        ekl = rr_kl * ak[kp] * al[lp] / akl;
                        ccekl = ekl - log_rr_kl - log_maxck[kp] - log_maxcl[lp];
                        if (ccekl > expcutoff) {
                                goto k_contracted;
                        }
                        envs->ak[0] = ak[kp];
                        rkl[0] = (ak[kp]*rk[0] + al[lp]*rl[0]) / akl;
                        rkl[1] = (ak[kp]*rk[1] + al[lp]*rl[1]) / akl;
                        rkl[2] = (ak[kp]*rk[2] + al[lp]*rl[2]) / akl;
                        eijcutoff = expcutoff - ccekl;
                        ekl = exp(-ekl);

                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                        } else {
                                fac1k = fac1l;
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
                                        if (pdata_ij->cceij > eijcutoff) {
                                                goto i_contracted;
                                        }
                                        envs->ai[0] = ai[ip];
                                        rij = pdata_ij->rij;
                                        cutoff = eijcutoff - pdata_ij->cceij;
                                        expijkl = pdata_ij->eij * ekl;
                                        if (i_ctr == 1) {
                                                fac1i = fac1j*ci[ip]*expijkl;
                                        } else {
                                                fac1i = fac1j*expijkl;
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
                                PRIM2CTR(k, gctrj, lenj);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR(l, gctrk, lenk);
                }
        } // end loop l_prim

        if (n_comp > 1 && !*lempty) {
                TRANSPOSE(gctrl);
        }
        return !*empty;
}


#define COMMON_ENVS_AND_DECLARE \
        FINT *shls = envs->shls; \
        FINT *bas = envs->bas; \
        double *env = envs->env; \
        FINT i_sh = shls[0]; \
        FINT j_sh = shls[1]; \
        FINT k_sh = shls[2]; \
        FINT l_sh = shls[3]; \
        CINTOpt *opt = envs->opt; \
        if (opt->pairdata != NULL && \
            ((opt->pairdata[i_sh*opt->nbas+j_sh] == NOVALUE) || \
             (opt->pairdata[k_sh*opt->nbas+l_sh] == NOVALUE))) { \
                return 0; \
        } \
        FINT i_ctr = envs->x_ctr[0]; \
        FINT j_ctr = envs->x_ctr[1]; \
        FINT k_ctr = envs->x_ctr[2]; \
        FINT l_ctr = envs->x_ctr[3]; \
        FINT i_prim = bas(NPRIM_OF, i_sh); \
        FINT j_prim = bas(NPRIM_OF, j_sh); \
        FINT k_prim = bas(NPRIM_OF, k_sh); \
        FINT l_prim = bas(NPRIM_OF, l_sh); \
        double *ai = env + bas(PTR_EXP, i_sh); \
        double *aj = env + bas(PTR_EXP, j_sh); \
        double *ak = env + bas(PTR_EXP, k_sh); \
        double *al = env + bas(PTR_EXP, l_sh); \
        double *ci = env + bas(PTR_COEFF, i_sh); \
        double *cj = env + bas(PTR_COEFF, j_sh); \
        double *ck = env + bas(PTR_COEFF, k_sh); \
        double *cl = env + bas(PTR_COEFF, l_sh); \
        double expcutoff = envs->expcutoff; \
        double rr_ij = SQUARE(envs->rirj); \
        double rr_kl = SQUARE(envs->rkrl); \
        PairData *_pdata_ij, *_pdata_kl, *pdata_ij, *pdata_kl; \
        if (opt->pairdata != NULL) { \
                _pdata_ij = opt->pairdata[i_sh*opt->nbas+j_sh]; \
                _pdata_kl = opt->pairdata[k_sh*opt->nbas+l_sh]; \
        } else { \
                double *log_maxci = opt->log_max_coeff[i_sh]; \
                double *log_maxcj = opt->log_max_coeff[j_sh]; \
                MALLOC_INSTACK(_pdata_ij, i_prim*j_prim + k_prim*l_prim); \
                if (CINTset_pairdata(_pdata_ij, ai, aj, envs->ri, envs->rj, \
                                     log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil, \
                                     i_prim, j_prim, rr_ij, expcutoff, env)) { \
                        return 0; \
                } \
                double *log_maxck = opt->log_max_coeff[k_sh]; \
                double *log_maxcl = opt->log_max_coeff[l_sh]; \
                _pdata_kl = _pdata_ij + i_prim*j_prim; \
                if (CINTset_pairdata(_pdata_kl, ak, al, envs->rk, envs->rl, \
                                     log_maxck, log_maxcl, envs->lk_ceil, envs->ll_ceil, \
                                     k_prim, l_prim, rr_kl, expcutoff, env)) { \
                        return 0; \
                } \
        } \
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor; \
        size_t nf = envs->nf; \
        double fac1i, fac1j, fac1k, fac1l; \
        FINT ip, jp, kp, lp; \
        FINT _empty[5] = {1, 1, 1, 1, 1}; \
        FINT *iempty = _empty + 0; \
        FINT *jempty = _empty + 1; \
        FINT *kempty = _empty + 2; \
        FINT *lempty = _empty + 3; \
        FINT *gempty = _empty + 4; \
        FINT *non0ctri = opt->non0ctr[i_sh]; \
        FINT *non0ctrj = opt->non0ctr[j_sh]; \
        FINT *non0ctrk = opt->non0ctr[k_sh]; \
        FINT *non0ctrl = opt->non0ctr[l_sh]; \
        FINT *non0idxi = opt->sortedidx[i_sh]; \
        FINT *non0idxj = opt->sortedidx[j_sh]; \
        FINT *non0idxk = opt->sortedidx[k_sh]; \
        FINT *non0idxl = opt->sortedidx[l_sh]; \
        double expij, expkl, eijcutoff, eklcutoff, cutoff; \
        eklcutoff = expcutoff; \
        double *rij, *rkl; \
        FINT *idx = opt->index_xyz_array[envs->i_l*LMAX1*LMAX1*LMAX1 \
                                        +envs->j_l*LMAX1*LMAX1 \
                                        +envs->k_l*LMAX1 \
                                        +envs->l_l]; \
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
                int lkl = envs->lk_ceil + envs->ll_ceil; \
                if (lij > 0) { \
                        double dist_ij = sqrt(rr_ij); \
                        double aij = ai[i_prim-1] + aj[j_prim-1]; \
                        double theta = omega2 / (omega2 + aij); \
                        expcutoff += lij * approx_log( \
                                (dist_ij+theta*r_guess+1.)/(dist_ij+1.)); \
                } \
                if (lkl > 0) { \
                        double dist_kl = sqrt(rr_kl); \
                        double akl = ak[k_prim-1] + al[l_prim-1]; \
                        double theta = omega2 / (omega2 + akl); \
                        expcutoff += lkl * approx_log( \
                                (dist_kl+theta*r_guess+1.)/(dist_kl+1.)); \
                } \
        }

#define SET_RIJ(I,J)    \
        if (pdata_##I##J->cceij > e##I##J##cutoff) { \
                goto I##_contracted; } \
        envs->a##I[0] = a##I[I##p]; \
        exp##I##J = pdata_##I##J->eij; \
        r##I##J = pdata_##I##J->rij;

// i_ctr = j_ctr = k_ctr = l_ctr = 1;
FINT CINT2e_1111_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
{
        COMMON_ENVS_AND_DECLARE;
        ADJUST_CUTOFF;
        FINT nc = 1;
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        size_t len0 = nf * n_comp;
        size_t len = leng + len0;
        double *gout;
        double *g;
        MALLOC_INSTACK(g, len);
        if (n_comp == 1) {
                gout = gctr;
                gempty = empty;
        } else {
                gout = g + leng;
        }

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al[0] = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];
                        eijcutoff = eklcutoff - pdata_kl->cceij;
                        pdata_ij = _pdata_ij;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj[0] = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        SET_RIJ(i, j);
                                        fac1i = fac1j*ci[ip]*expij*expkl;
                                        envs->fac[0] = fac1i;
                                        cutoff = eijcutoff - pdata_ij->cceij;
                                        if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs)) {
                                                (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                                *gempty = 0;
                                        }
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        if (n_comp > 1 && !*gempty) {
                TRANSPOSE(gout);
        }
        return !*empty;
}

// i_ctr = n; j_ctr = k_ctr = l_ctr = 1;
FINT CINT2e_n111_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
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

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al[0] = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];
                        eijcutoff = eklcutoff - pdata_kl->cceij;
                        pdata_ij = _pdata_ij;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj[0] = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        if (pdata_ij->cceij > eijcutoff) {
                                                goto i_contracted;
                                        }
                                        SET_RIJ(i, j);
                                        cutoff = eijcutoff - pdata_ij->cceij;
                                        fac1i = fac1j*expij*expkl;
                                        envs->fac[0] = fac1i;
                                        if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs)) {
                                                (*envs->f_gout)(gout, g, idx, envs, 1);
                                                PRIM2CTR(i, gout, len0);
                                        }
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        if (n_comp > 1 && !*iempty) {
                TRANSPOSE(gctri);
        }
        return !*empty;
}

// j_ctr = n; i_ctr = k_ctr = l_ctr = 1;
FINT CINT2e_1n11_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
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

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al[0] = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];
                        eijcutoff = eklcutoff - pdata_kl->cceij;
                        pdata_ij = _pdata_ij;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj[0] = aj[jp];
                                fac1j = fac1k;
                                *iempty = 1;
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        SET_RIJ(i, j);
                                        cutoff = eijcutoff - pdata_ij->cceij;
                                        fac1i = fac1j*ci[ip]*expij*expkl;
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
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        if (n_comp > 1 && !*jempty) {
                TRANSPOSE(gctrj);
        }
        return !*empty;
}

// k_ctr = n; i_ctr = j_ctr = l_ctr = 1;
FINT CINT2e_11n1_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
{
        COMMON_ENVS_AND_DECLARE;
        ADJUST_CUTOFF;
        FINT nc = k_ctr;
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        size_t lenk = nf * k_ctr * n_comp; // gctrk
        size_t len0 = nf * n_comp; // gout
        size_t len = leng + lenk + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *g1 = g + leng;
        double *gout, *gctrk;
        ALIAS_ADDR_IF_EQUAL(k, m);
        gout = g1;

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al[0] = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l;
                        eijcutoff = eklcutoff - pdata_kl->cceij;
                        pdata_ij = _pdata_ij;
                        *jempty = 1;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj[0] = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        SET_RIJ(i, j);
                                        cutoff = eijcutoff - pdata_ij->cceij;
                                        fac1i = fac1j*ci[ip]*expij*expkl;
                                        envs->fac[0] = fac1i;
                                        if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs)) {
                                                (*envs->f_gout)(gout, g, idx, envs, *jempty);
                                                *jempty = 0;
                                        }
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR(k, gout, len0);
                        }
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        if (n_comp > 1 && !*kempty) {
                TRANSPOSE(gctrk);
        }
        return !*empty;
}

// l_ctr = n; i_ctr = j_ctr = k_ctr = 1;
FINT CINT2e_111n_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
{
        COMMON_ENVS_AND_DECLARE;
        ADJUST_CUTOFF;
        FINT nc = l_ctr;
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        size_t lenl = nf * l_ctr * n_comp; // gctrl
        size_t len0 = nf * n_comp; // gout
        size_t len = leng + lenl + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *g1 = g + leng;
        double *gout, *gctrl;
        ALIAS_ADDR_IF_EQUAL(l, m);
        gout = g1;

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al[0] = al[lp];
                fac1l = envs->common_factor;
                *kempty = 1;
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];
                        eijcutoff = eklcutoff - pdata_kl->cceij;
                        pdata_ij = _pdata_ij;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj[0] = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        SET_RIJ(i, j);
                                        cutoff = eijcutoff - pdata_ij->cceij;
                                        fac1i = fac1j*ci[ip]*expij*expkl;
                                        envs->fac[0] = fac1i;
                                        if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs)) {
                                                (*envs->f_gout)(gout, g, idx, envs, *kempty);
                                                *kempty = 0;
                                        }
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR(l, gout, len0);
                }
        } // end loop l_prim

        if (n_comp > 1 && !*lempty) {
                TRANSPOSE(gctrl);
        }
        return !*empty;
}


FINT CINT2e_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
{
        COMMON_ENVS_AND_DECLARE;
        ADJUST_CUTOFF;
        FINT nc = i_ctr * j_ctr * k_ctr * l_ctr;
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1); // (irys,i,j,k,l,coord,0:1);
        size_t lenl = nf * nc * n_comp; // gctrl
        size_t lenk = nf * i_ctr * j_ctr * k_ctr * n_comp; // gctrk
        size_t lenj = nf * i_ctr * j_ctr * n_comp; // gctrj
        size_t leni = nf * i_ctr * n_comp; // gctri
        size_t len0 = nf * n_comp; // gout
        size_t len = leng + lenl + lenk + lenj + leni + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk, *gctrl;

        ALIAS_ADDR_IF_EQUAL(l, m);
        ALIAS_ADDR_IF_EQUAL(k, l);
        ALIAS_ADDR_IF_EQUAL(j, k);
        ALIAS_ADDR_IF_EQUAL(i, j);
        ALIAS_ADDR_IF_EQUAL(g, i);

        pdata_kl = _pdata_kl;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al[0] = al[lp];
                if (l_ctr == 1) {
                        fac1l = envs->common_factor * cl[lp];
                } else {
                        fac1l = envs->common_factor;
                        *kempty = 1;
                }
                for (kp = 0; kp < k_prim; kp++, pdata_kl++) {
                        /* SET_RIJ(k, l); */
                        if (pdata_kl->cceij > eklcutoff) {
                                goto k_contracted;
                        }
                        envs->ak[0] = ak[kp];
                        expkl = pdata_kl->eij;
                        rkl = pdata_kl->rij;
                        eijcutoff = eklcutoff - pdata_kl->cceij;
                        /* SET_RIJ(k, l); end */
                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                        } else {
                                fac1k = fac1l;
                                *jempty = 1;
                        }

                        pdata_ij = _pdata_ij;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj[0] = aj[jp];
                                if (j_ctr == 1) {
                                        fac1j = fac1k * cj[jp];
                                } else {
                                        fac1j = fac1k;
                                        *iempty = 1;
                                }
                                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                        /* SET_RIJ(i, j); */
                                        if (pdata_ij->cceij > eijcutoff) {
                                                goto i_contracted;
                                        }
                                        envs->ai[0] = ai[ip];
                                        expij = pdata_ij->eij;
                                        rij = pdata_ij->rij;
                                        /* SET_RIJ(i, j); end */
                                        cutoff = eijcutoff - pdata_ij->cceij;
                                        if (i_ctr == 1) {
                                                fac1i = fac1j*ci[ip] * expij*expkl;
                                        } else {
                                                fac1i = fac1j * expij*expkl;
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
                                PRIM2CTR(k, gctrj, lenj);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR(l, gctrk, lenk);
                }
        } // end loop l_prim

        if (n_comp > 1 && !*lempty) {
                TRANSPOSE(gctrl);
        }
        return !*empty;
}


static FINT (*CINTf_2e_loop[16])(double *, CINTEnvVars *, double *, FINT *) = {
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_n111_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_1n11_loop,
        CINT2e_loop,
        CINT2e_11n1_loop,
        CINT2e_111n_loop,
        CINT2e_1111_loop,
};

#define PAIRDATA_NON0IDX_SIZE(ps) \
                FINT *bas = envs->bas; \
                FINT *shls  = envs->shls; \
                FINT i_prim = bas(NPRIM_OF, shls[0]); \
                FINT j_prim = bas(NPRIM_OF, shls[1]); \
                FINT k_prim = bas(NPRIM_OF, shls[2]); \
                FINT l_prim = bas(NPRIM_OF, shls[3]); \
                size_t ps = ((i_prim*j_prim + k_prim*l_prim) * 5 \
                           + i_prim * x_ctr[0] \
                           + j_prim * x_ctr[1] \
                           + k_prim * x_ctr[2] \
                           + l_prim * x_ctr[3] \
                           +(i_prim+j_prim+k_prim+l_prim)*2 + nf*3);

CACHE_SIZE_T CINT2e_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                      double *cache, void (*f_c2s)())
{
        FINT *x_ctr = envs->x_ctr;
        size_t nf = envs->nf;
        size_t nc = nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1);
                size_t len0 = nf*n_comp;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                        nc*n_comp+nf*4);
#if !defined(I8) && !defined(CACHE_SIZE_I8)
                if (cache_size >= INT32_MAX) {
                        fprintf(stderr, "CINT2e_drv cache_size overflow: "
                                "cache_size %zu > %d, nf %zu, nc %zu, n_comp %d\n",
                                cache_size, INT32_MAX, nf, nc, (int)n_comp);
                        cache_size = 0;
                }
#endif
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1);
                size_t len0 = nf*n_comp;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                        nc*n_comp+nf*4);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        FINT n;
        FINT empty = 1;
        if (opt != NULL) {
                envs->opt = opt;
                n = ((x_ctr[0]==1) << 3) + ((x_ctr[1]==1) << 2)
                  + ((x_ctr[2]==1) << 1) +  (x_ctr[3]==1);
                CINTf_2e_loop[n](gctr, envs, cache, &empty);
        } else {
                CINT2e_loop_nopt(gctr, envs, cache, &empty);
        }

        FINT counts[4];
        if (f_c2s == &c2s_sph_2e1) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->j_l*2+1) * x_ctr[1];
                counts[2] = (envs->k_l*2+1) * x_ctr[2];
                counts[3] = (envs->l_l*2+1) * x_ctr[3];
        } else {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfj * x_ctr[1];
                counts[2] = envs->nfk * x_ctr[2];
                counts[3] = envs->nfl * x_ctr[3];
        }
        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1] * dims[2] * dims[3];
        if (!empty) {
                for (n = 0; n < n_comp; n++) {
                        (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
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
CACHE_SIZE_T CINT2e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                      double *cache, void (*f_e1_c2s)(), void (*f_e2_c2s)())
{
        FINT *shls = envs->shls;
        FINT *bas = envs->bas;
        FINT counts[4];
        counts[0] = CINTcgto_spinor(shls[0], bas);
        counts[1] = CINTcgto_spinor(shls[1], bas);
        counts[2] = CINTcgto_spinor(shls[2], bas);
        counts[3] = CINTcgto_spinor(shls[3], bas);
        FINT *x_ctr = envs->x_ctr;
        size_t nf = envs->nf;
        size_t nc = nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        FINT n1 = counts[0] * envs->nfk * x_ctr[2]
                           * envs->nfl * x_ctr[3] * counts[1];
        if (out == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1);
                size_t len0 = nf*n_comp;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                     nc*n_comp + n1*envs->ncomp_e2*OF_CMPLX
                                     + nf*32*OF_CMPLX);
#if !defined(I8) && !defined(CACHE_SIZE_I8)
                if (cache_size >= INT32_MAX) {
                        fprintf(stderr, "CINT2e_drv cache_size overflow: "
                                "cache_size %zu > %d, nf %zu, nc %zu, n_comp %d\n",
                                cache_size, INT32_MAX, nf, nc, (int)n_comp);
                        cache_size = 0;
                }
#endif
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1);
                size_t len0 = nf*n_comp;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                     nc*n_comp + n1*envs->ncomp_e2*OF_CMPLX
                                     + nf*32*OF_CMPLX);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        FINT n, m;
        FINT empty = 1;
        if (opt != NULL) {
                envs->opt = opt;
                n = ((x_ctr[0]==1) << 3) + ((x_ctr[1]==1) << 2)
                  + ((x_ctr[2]==1) << 1) +  (x_ctr[3]==1);
                CINTf_2e_loop[n](gctr, envs, cache, &empty);
        } else {
                CINT2e_loop_nopt(gctr, envs, cache, &empty);
        }

        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1] * dims[2] * dims[3];
        if (!empty) {
                double complex *opij;
                MALLOC_INSTACK(opij, n1*envs->ncomp_e2);
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        for (m = 0; m < envs->ncomp_e2; m++) {
                                (*f_e1_c2s)(opij+n1*m, gctr, dims, envs, cache);
                                gctr += nc * envs->ncomp_e1;
                        }
                        (*f_e2_c2s)(out+nout*n, opij, dims, envs, cache);
                }
        } else {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return !empty;
}


/*
 * <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
#if __SSE3__
#include "gout2e_simd.c"
#else
void CINTgout2e(double *gout, double *g, FINT *idx,
                CINTEnvVars *envs, FINT gout_empty)
{
        FINT nf = envs->nf;
        FINT i, ix, iy, iz, n;
        double s;

        if (gout_empty) {
                switch (envs->nrys_roots) {
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
                        case 8:
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
                                                + g[ix+6] * g[iy+6] * g[iz+6]
                                                + g[ix+7] * g[iy+7] * g[iz+7];
                                }
                                break;
                        default:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        s = 0;
                                        for (i = 0; i < envs->nrys_roots; i++) {
                                                s += g[ix+i] * g[iy+i] * g[iz+i];
                                        }
                                        gout[n] = s;
                                }
                                break;
                } // end switch nroots
        } else { // not flag_acc
                switch (envs->nrys_roots) {
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
                        case 8:
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
                                                + g[ix+6] * g[iy+6] * g[iz+6]
                                                + g[ix+7] * g[iy+7] * g[iz+7];
                                }
                                break;
                        default:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        s = 0;
                                        for (i = 0; i < envs->nrys_roots; i++) {
                                                s += g[ix+i] * g[iy+i] * g[iz+i];
                                        }
                                        gout[n] += s;
                                }
                                break;
                } // end switch nroots
        }
}
#endif

CACHE_SIZE_T int2e_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
              FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_drv(out, dims, &envs, opt, cache, &c2s_sph_2e1);
}
void int2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                     FINT *bas, FINT nbas, double *env)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

CACHE_SIZE_T int2e_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
               FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_drv(out, dims, &envs, opt, cache, &c2s_cart_2e1);
}

/*
 * spinor <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
CACHE_SIZE_T int2e_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                 FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_spinor_drv(out, dims, &envs, opt, cache,
                                 &c2s_sf_2e1, &c2s_sf_2e2);
}


ALL_CINT(int2e)
ALL_CINT_FORTRAN_(int2e)

