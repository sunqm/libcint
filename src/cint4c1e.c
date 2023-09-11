/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 4-center 1-electron integrals
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "g1e.h"
#include "optimizer.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"

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

void CINTg4c1e_ovlp(double *g, CINTEnvVars *envs, double *cache);
void CINTg4c1e_index_xyz(FINT *idx, CINTEnvVars *envs);
CACHE_SIZE_T CINTinit_int4c1e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                             FINT *atm, FINT natm,
                             FINT *bas, FINT nbas, double *env);
void CINTgout1e(double *gout, double *g, FINT *idx, CINTEnvVars *envs, FINT empty);

CACHE_SIZE_T CINT4c1e_loop_nopt(double *gctr, CINTEnvVars *envs, double *cache)
{
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
        double *ri = envs->ri;
        double *rj = envs->rj;
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
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k, fac1l;
        FINT ip, jp, kp, lp;
        FINT empty[5] = {1, 1, 1, 1, 1};
        FINT *iempty = empty + 0;
        FINT *jempty = empty + 1;
        FINT *kempty = empty + 2;
        FINT *lempty = empty + 3;
        FINT *gempty = empty + 4;
        /* COMMON_ENVS_AND_DECLARE end */

        FINT *idx;
        MALLOC_INSTACK(idx, envs->nf * 3);
        CINTg4c1e_index_xyz(idx, envs);

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
        FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        FINT lenl = envs->nf * nc * n_comp; // gctrl
        FINT lenk = envs->nf * i_ctr * j_ctr * k_ctr * n_comp; // gctrk
        FINT lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        FINT leni = envs->nf * i_ctr * n_comp; // gctri
        FINT len0 = envs->nf * n_comp; // gout
        FINT len = leng + lenl + lenk + lenj + leni + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk, *gctrl;

        if (n_comp == 1) {
                gctrl = gctr;
        } else {
                gctrl = g1;
                g1 += lenl;
        }
        if (l_ctr == 1) {
                gctrk = gctrl;
                kempty = lempty;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctrj = gctrk;
                jempty = kempty;
        } else {
                gctrj = g1;
                g1 += lenj;
        }
        if (j_ctr == 1) {
                gctri = gctrj;
                iempty = jempty;
        } else {
                gctri = g1;
                g1 += leni;
        }
        if (i_ctr == 1) {
                gout = gctri;
                gempty = iempty;
        } else {
                gout = g1;
        }

        double eij, ekl, expijkl, eijkl, a0, aij, akl;
        double rijrkl[3];
        double rij[3];
        double rkl[3];
        double dist_ij = SQUARE(envs->rirj);
        double dist_kl = SQUARE(envs->rkrl);
        double common_fac = envs->common_factor * SQRTPI * M_PI *
                CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l) *
                CINTcommon_fac_sp(envs->k_l) * CINTcommon_fac_sp(envs->l_l);

        *lempty = 1;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al[0] = al[lp];
                if (l_ctr == 1) {
                        fac1l = common_fac * cl[lp];
                } else {
                        fac1l = common_fac;
                        *kempty = 1;
                }
                for (kp = 0; kp < k_prim; kp++) {
                        envs->ak[0] = ak[kp];
                        akl = ak[kp] + al[lp];
                        ekl = dist_kl * ak[kp] * al[lp] / akl;
                        if (ekl > EXPCUTOFF) {
                                goto k_contracted;
                        }
                        rkl[0] = (ak[kp]*rk[0] + al[lp]*rl[0]) / akl;
                        rkl[1] = (ak[kp]*rk[1] + al[lp]*rl[1]) / akl;
                        rkl[2] = (ak[kp]*rk[2] + al[lp]*rl[2]) / akl;
                        envs->rkl[0] = rkl[0];
                        envs->rkl[1] = rkl[1];
                        envs->rkl[2] = rkl[2];
                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                        } else {
                                fac1k = fac1l;
                                *jempty = 1;
                        }

                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj[0] = aj[jp];
                                if (j_ctr == 1) {
                                        fac1j = fac1k * cj[jp];
                                } else {
                                        fac1j = fac1k;
                                        *iempty = 1;
                                }
                                for (ip = 0; ip < i_prim; ip++) {
                                        envs->ai[0] = ai[ip];
                                        aij = ai[ip] + aj[jp];
                                        eij = dist_ij * ai[ip] * aj[jp] / aij;
                                        if (eij+ekl > EXPCUTOFF) {
                                                goto i_contracted;
                                        }
                                        rij[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / aij;
                                        rij[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / aij;
                                        rij[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / aij;
                                        envs->rij[0] = rij[0];
                                        envs->rij[1] = rij[1];
                                        envs->rij[2] = rij[2];
                                        rijrkl[0] = rij[0] - rkl[0];
                                        rijrkl[1] = rij[1] - rkl[1];
                                        rijrkl[2] = rij[2] - rkl[2];
                                        a0 = aij * akl / (aij + akl);
                                        eijkl = a0 * SQUARE(rijrkl);
                                        eijkl += eij + ekl;
                                        if (eijkl > EXPCUTOFF) {
                                                goto i_contracted;
                                        }
                                        expijkl = exp(-eijkl);
                                        if (i_ctr == 1) {
                                                fac1i = fac1j*ci[ip]*expijkl;
                                        } else {
                                                fac1i = fac1j*expijkl;
                                        }
                                        envs->fac[0] = fac1i;
                                        CINTg4c1e_ovlp(g, envs, cache);
                                        (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                        PRIM2CTR0(i, gout, envs->nf*n_comp);
i_contracted: ;
                                } // end loop i_prim
                                if (!*iempty) {
                                        PRIM2CTR0(j, gctri, envs->nf*i_ctr*n_comp);
                                }
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR0(k, gctrj,envs->nf*i_ctr*j_ctr*n_comp);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR0(l, gctrk, envs->nf*i_ctr*j_ctr*k_ctr*n_comp);
                }
        } // end loop l_prim

        if (n_comp > 1 && !*lempty) {
                CINTdmat_transpose(gctr, gctrl, envs->nf*nc, n_comp);
        }
        return !*lempty;
}


#define PAIRDATA_NON0IDX_SIZE(ps) \
                FINT *bas = envs->bas; \
                FINT *shls  = envs->shls; \
                FINT i_prim = bas(NPRIM_OF, shls[0]); \
                FINT j_prim = bas(NPRIM_OF, shls[1]); \
                FINT k_prim = bas(NPRIM_OF, shls[2]); \
                FINT l_prim = bas(NPRIM_OF, shls[3]); \
                FINT ps = (i_prim * x_ctr[0] \
                           + j_prim * x_ctr[1] \
                           + k_prim * x_ctr[2] \
                           + l_prim * x_ctr[3] \
                           + envs->nf*3);

CACHE_SIZE_T CINT4c1e_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                          double *cache, void (*f_c2s)())
{
        FINT *x_ctr = envs->x_ctr;
        size_t nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                CACHE_SIZE_T leng = envs->g_size*3*((1<<envs->gbits)+1);
                CACHE_SIZE_T len0 = envs->nf*n_comp;
                CACHE_SIZE_T cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                      nc*n_comp+envs->nf*32*OF_CMPLX);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1);
                size_t len0 = envs->nf*n_comp;
                size_t cache_size = MAX(leng+len0+nc*n_comp*3 + pdata_size,
                                      nc*n_comp+envs->nf*32*OF_CMPLX);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        FINT n, has_value;
        has_value = CINT4c1e_loop_nopt(gctr, envs, cache);

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
        if (has_value) {
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
        return has_value;
}

CACHE_SIZE_T int4c1e_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                 FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int4c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT4c1e_drv(out, dims, &envs, opt, cache, &c2s_sph_2e1);
}
void int4c1e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}

CACHE_SIZE_T int4c1e_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                  FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int4c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT4c1e_drv(out, dims, &envs, opt, cache, &c2s_cart_2e1);
}

CACHE_SIZE_T int4c1e_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                   FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        fprintf(stderr, "int4c1e_spinor not implemented\n");
        exit(1);
}

ALL_CINT(int4c1e)
ALL_CINT_FORTRAN_(int4c1e)
