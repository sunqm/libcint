/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO integrals
 */

#include <stdlib.h>
#include <math.h>
#include "cint_bas.h"
#include "optimizer.h"
#include "g1e.h"
#include "cint1e.h"
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

static void make_g1e_gout(double *gout, double *g, FINT *idx,
                          CINTEnvVars *envs, FINT empty, FINT int1e_type);

/*
 * 1e GTO integral basic loop for < i|j>, no 1/r
 */
FINT CINT1e_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT int1e_type)
{
        FINT *shls = envs->shls;
        FINT *bas = envs->bas;
        double *env = envs->env;
        FINT i_sh = shls[0];
        FINT j_sh = shls[1];
        FINT i_ctr = envs->x_ctr[0];
        FINT j_ctr = envs->x_ctr[1];
        FINT i_prim = bas(NPRIM_OF, i_sh);
        FINT j_prim = bas(NPRIM_OF, j_sh);
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;

        double expcutoff = envs->expcutoff;
        double *log_maxci, *log_maxcj;
        PairData *pdata_base, *pdata_ij;
        MALLOC_INSTACK(log_maxci, i_prim+j_prim);
        MALLOC_INSTACK(pdata_base, i_prim*j_prim);
        log_maxcj = log_maxci + i_prim;
        CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
        CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
        if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                             log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                             i_prim, j_prim, SQUARE(envs->rirj), expcutoff, env)) {
                return 0;
        }

        double fac1i, fac1j, expij;
        FINT ip, jp;
        FINT empty[4] = {1, 1, 1, 1};
        FINT *gempty = empty + 0;
        FINT *iempty = empty + 1;
        FINT *jempty = empty + 2;
        double *rij;
        FINT *idx;
        MALLOC_INSTACK(idx, envs->nf * 3);
        CINTg1e_index_xyz(idx, envs);

        FINT *non0ctri, *non0ctrj;
        FINT *non0idxi, *non0idxj;
        MALLOC_INSTACK(non0ctri, i_prim+j_prim+i_prim*i_ctr+j_prim*j_ctr);
        non0ctrj = non0ctri + i_prim;
        non0idxi = non0ctrj + j_prim;
        non0idxj = non0idxi + i_prim*i_ctr;
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);

        const FINT nc = i_ctr * j_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenj = envs->nf * nc * n_comp; // gctrj
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenj + leni + len0;
        double *g, *gout, *gctri, *gctrj;
        MALLOC_INSTACK(g, len);  // must be allocated last in this function
        double *g1 = g + leng;
        if (n_comp == 1) {
                gctrj = gctr;
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

        double common_factor = envs->common_factor
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l);

        pdata_ij = pdata_base;
        for (jp = 0; jp < j_prim; jp++) {
                envs->aj[0] = aj[jp];
                if (j_ctr == 1) {
                        fac1j = common_factor * cj[jp];
                } else {
                        fac1j = common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                        if (pdata_ij->cceij > expcutoff) {
                                continue;
                        }
                        envs->ai[0] = ai[ip];
                        expij = pdata_ij->eij;
                        rij = pdata_ij->rij;
                        envs->rij[0] = rij[0];
                        envs->rij[1] = rij[1];
                        envs->rij[2] = rij[2];
                        if (i_ctr == 1) {
                                fac1i = fac1j*ci[ip]*expij;
                        } else {
                                fac1i = fac1j*expij;
                        }
                        envs->fac[0] = fac1i;
                        make_g1e_gout(gout, g, idx, envs, *gempty, int1e_type);
                        PRIM2CTR0(i, gout, envs->nf*n_comp);
                }
                if (!*iempty) {
                        PRIM2CTR0(j, gctri, envs->nf*i_ctr*n_comp);
                }
        }

        if (n_comp > 1 && !*jempty) {
                CINTdmat_transpose(gctr, gctrj, envs->nf*nc, n_comp);
        }
        return !*jempty;
}


CACHE_SIZE_T int1e_cache_size(CINTEnvVars *envs)
{
        FINT *shls = envs->shls;
        FINT *bas = envs->bas;
        FINT i_prim = bas(NPRIM_OF, shls[0]);
        FINT j_prim = bas(NPRIM_OF, shls[1]);
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
        FINT lenj = envs->nf * nc * n_comp;
        FINT leni = envs->nf * x_ctr[0] * n_comp;
        FINT len0 = envs->nf*n_comp;
        FINT pdata_size = (i_prim*j_prim * 5
                           + i_prim * x_ctr[0]
                           + j_prim * x_ctr[1]
                           +(i_prim+j_prim)*2 + envs->nf*3);
        FINT cache_size = MAX(nc*n_comp + leng+lenj+leni+len0 + pdata_size,
                             nc*n_comp + envs->nf*8*OF_CMPLX);
        return cache_size;
}

/*
 * 1e integrals <i|O|j> without 1/r
 */
CACHE_SIZE_T CINT1e_drv(double *out, FINT *dims, CINTEnvVars *envs,
               double *cache, void (*f_c2s)(), FINT int1e_type)
{
        if (out == NULL) {
                return int1e_cache_size(envs);
        }
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *stack = NULL;
        if (cache == NULL) {
                size_t cache_size = int1e_cache_size(envs);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        FINT has_value = CINT1e_loop(gctr, envs, cache, int1e_type);

        FINT counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        if (f_c2s == &c2s_sph_1e) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->j_l*2+1) * x_ctr[1];
        } else if (f_c2s == &c2s_cart_1e) {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfj * x_ctr[1];
        }
        counts[2] = 1;
        counts[3] = 1;
        FINT nout = dims[0] * dims[1];
        FINT n;
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

CACHE_SIZE_T CINT1e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs,
                       double *cache, void (*f_c2s)(), FINT int1e_type)
{
        if (out == NULL) {
                return int1e_cache_size(envs);
        }
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1] * envs->ncomp_e1;
        double *stack = NULL;
        if (cache == NULL) {
                size_t cache_size = int1e_cache_size(envs);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*envs->ncomp_tensor);

        FINT has_value = CINT1e_loop(gctr, envs, cache, int1e_type);

        FINT counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        counts[2] = 1;
        counts[3] = 1;
        FINT nout = dims[0] * dims[1];
        FINT n;
        if (has_value) {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }

        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}

static void make_g1e_gout(double *gout, double *g, FINT *idx,
                          CINTEnvVars *envs, FINT empty, FINT int1e_type)
{
        FINT ia;
        switch (int1e_type) {
        case 0:
                CINTg1e_ovlp(g, envs);
                (*envs->f_gout)(gout, g, idx, envs, empty);
                break;
        case 1:
                CINTg1e_nuc(g, envs, -1);
                (*envs->f_gout)(gout, g, idx, envs, empty);
                break;
        case 2:
                for (ia = 0; ia < envs->natm; ia++) {
                        CINTg1e_nuc(g, envs, ia);
                        (*envs->f_gout)(gout, g, idx, envs, (empty && ia == 0));
                }
                break;
        }
}

void CINTgout1e(double *gout, double *g, FINT *idx, CINTEnvVars *envs, FINT empty)
{
        FINT nf = envs->nf;
        FINT n, ix, iy, iz;
        if (empty) {
                for (n = 0; n < nf; n++) {
                        ix = idx[n*3+0];
                        iy = idx[n*3+1];
                        iz = idx[n*3+2];
                        gout[n] = g[ix] * g[iy] * g[iz];
                }
        } else {
                for (n = 0; n < nf; n++) {
                        ix = idx[n*3+0];
                        iy = idx[n*3+1];
                        iz = idx[n*3+2];
                        gout[n] += g[ix] * g[iy] * g[iz];
                }
        }
}

void CINTgout1e_nuc(double *gout, double *g, FINT *idx, CINTEnvVars *envs, FINT empty)
{
        FINT nf = envs->nf;
        FINT nrys_roots = envs->nrys_roots;
        FINT n, i;
        double *gx, *gy, *gz;
        double s;

        if (empty) {
                for (n = 0; n < nf; n++) {
                        gx = g + idx[n*3+0];
                        gy = g + idx[n*3+1];
                        gz = g + idx[n*3+2];
                        s = 0;
                        for (i = 0; i < nrys_roots; i++) {
                                s += gx[i] * gy[i] * gz[i];
                        }
                        gout[n] = s;
                }
        } else {
                for (n = 0; n < nf; n++) {
                        gx = g + idx[n*3+0];
                        gy = g + idx[n*3+1];
                        gz = g + idx[n*3+2];
                        s = 0;
                        for (i = 0; i < nrys_roots; i++) {
                                s += gx[i] * gy[i] * gz[i];
                        }
                        gout[n] += s;
                }
        }
}

CACHE_SIZE_T int1e_ovlp_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                     FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT1e_drv(out, dims, &envs, cache, &c2s_sph_1e, 0);
}

CACHE_SIZE_T int1e_ovlp_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                     FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT1e_drv(out, dims, &envs, cache, &c2s_cart_1e, 0);
}

CACHE_SIZE_T int1e_ovlp_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                     FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT1e_spinor_drv(out, dims, &envs, cache, &c2s_sf_1e, 0);
}

void int1e_ovlp_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                          FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}

CACHE_SIZE_T int1e_nuc_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                     FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_nuc;
        return CINT1e_drv(out, dims, &envs, cache, &c2s_sph_1e, 2);
}

CACHE_SIZE_T int1e_nuc_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                     FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_nuc;
        return CINT1e_drv(out, dims, &envs, cache, &c2s_cart_1e, 2);
}

CACHE_SIZE_T int1e_nuc_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                     FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_nuc;
        return CINT1e_spinor_drv(out, dims, &envs, cache, &c2s_sf_1e, 2);
}

void int1e_nuc_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                          FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}


ALL_CINT(int1e_ovlp);
ALL_CINT(int1e_nuc);
ALL_CINT_FORTRAN_(int1e_ovlp);
ALL_CINT_FORTRAN_(int1e_nuc);
