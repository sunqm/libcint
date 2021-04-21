/*
 * Copyright (C) 2021  Qiming Sun <osirpt.sun@gmail.com>
 *
 * <i|1/r|j> integrals for multiple grids
 */

#include <stdlib.h>
#include <math.h>
#include "cint_bas.h"
#include "optimizer.h"
#include "g1e.h"
#include "g1e_grids.h"
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

static void _transpose_comps(double *gctr, double *gctrj,
                             FINT bgrids, FINT dij, FINT ngrids, FINT n_comp);

FINT CINT1e_grids_loop(double *gctr, CINTEnvVars *envs, double fac, double *cache)
{
        FINT *shls  = envs->shls;
        FINT *atm = envs->atm;
        FINT *bas = envs->bas;
        double *env = envs->env;
        FINT i_sh = shls[0];
        FINT j_sh = shls[1];
        FINT i_l = envs->i_l;
        FINT j_l = envs->j_l;
        FINT i_ctr = envs->x_ctr[0];
        FINT j_ctr = envs->x_ctr[1];
        FINT i_prim = bas(NPRIM_OF, i_sh);
        FINT j_prim = bas(NPRIM_OF, j_sh);
        FINT nf = envs->nf;
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        FINT ngrids = envs->ngrids;
        double *grids = envs->grids;
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);

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
                             i_prim, j_prim, SQUARE(envs->rirj), expcutoff)) {
                return 0;
        }

        double fac1i, fac1j, expij;
        double *rij;
        FINT ip, jp, i, n, grids_offset, bgrids;
        FINT empty[4] = {1, 1, 1, 1};
        FINT *gempty = empty + 0;
        FINT *iempty = empty + 1;
        FINT *jempty = empty + 2;
        FINT all_empty = 1;

        FINT *idx;
        MALLOC_INSTACK(idx, nf * 3);
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
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenj = GRID_BLKSIZE * nf * nc * n_comp; // gctrj
        const FINT leni = GRID_BLKSIZE * nf * i_ctr * n_comp; // gctri
        const FINT len0 = GRID_BLKSIZE * nf * n_comp; // gout
        const FINT len = leng + lenj + leni + len0;
        double *g;
        MALLOC_ALIGN8_INSTACK(g, len);  // must be allocated last in this function
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj;

        for (grids_offset = 0; grids_offset < ngrids; grids_offset += GRID_BLKSIZE) {
                envs->grids_offset = grids_offset;
                bgrids = MIN(ngrids - grids_offset, GRID_BLKSIZE);
                empty[0] = 1;
                empty[1] = 1;
                empty[2] = 1;
                if (n_comp == 1) {
                        gctrj = gctr + grids_offset * nf * nc;
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
                pdata_ij = pdata_base;
                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = envs->common_factor * cj[jp];
                        } else {
                                fac1j = envs->common_factor;
                                *iempty = 1;
                        }
                        for (ip = 0; ip < i_prim; ip++, pdata_ij++) {
                                if (pdata_ij->cceij > expcutoff) {
                                        continue;
                                }
                                envs->ai = ai[ip];
                                envs->aij = ai[ip] + aj[jp];
                                expij = pdata_ij->eij;
                                rij = pdata_ij->rij;
                                envs->rij = rij;
                                envs->rijrx[0] = rij[0] - envs->rx_in_rijrx[0];
                                envs->rijrx[1] = rij[1] - envs->rx_in_rijrx[1];
                                envs->rijrx[2] = rij[2] - envs->rx_in_rijrx[2];
                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*expij;
                                } else {
                                        fac1i = fac1j*expij;
                                }

                                CINTg0_1e_grids(g, fac1i, envs, cache);
                                (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                PRIM2CTR0(i, gout, bgrids * nf * n_comp);
                        }
                        if (!*iempty) {
                                PRIM2CTR0(j, gctri, bgrids * nf * i_ctr * n_comp);
                        }
                }
                if (n_comp > 1 && !*jempty) {
                        _transpose_comps(gctr+grids_offset*nf*nc, gctrj,
                                         bgrids, nf*nc, ngrids, n_comp);
                }
                all_empty &= *jempty;
        }
        return !all_empty;
}

static void _transpose_comps(double *gctr, double *gctrj,
                             FINT bgrids, FINT dij, FINT ngrids, FINT n_comp)
{
        FINT n, ic, ig;
        double *pgctr, *pgctrj;
        for (ic = 0; ic < n_comp; ic++) {
                for (n = 0; n < dij; n++) {
                        pgctr = gctr + (ic * dij + n) * ngrids;
                        pgctrj = gctrj + (n * n_comp + ic) * bgrids;
                        for (ig = 0; ig < bgrids; ig++) {
                                pgctr[ig] = pgctrj[ig];
                        }
                }
        } 
}

FINT int1e_grids_cache_size(CINTEnvVars *envs)
{
        FINT *x_ctr = envs->x_ctr;
        FINT ngrids = envs->ngrids;
        FINT nroots = envs->nrys_roots;
        FINT nf = envs->nf;
        FINT nc = ngrids * nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
        FINT len0 = GRID_BLKSIZE * nf * n_comp;
        FINT leni = len0 * x_ctr[0];
        FINT lenj = leni * x_ctr[1];
        FINT cache_size = MAX(nc*n_comp + leng + len0 + leni + lenj + GRID_BLKSIZE*nroots*3 + GRID_BLKSIZE*6 + 32,
                              nc*n_comp + GRID_BLKSIZE * nf*8*OF_CMPLX);
        return cache_size + 10000;
}

/*
 * 1e integrals <i|O|j> without 1/r
 */
FINT CINT1e_grids_drv(double *out, FINT *dims, CINTEnvVars *envs,
                      double *cache, void (*f_c2s)())
{
        FINT ngrids = envs->ngrids;
        if (out == NULL) {
                return int1e_grids_cache_size(envs);
        }
        FINT *x_ctr = envs->x_ctr;
        FINT ngrids_nf = envs->ngrids * envs->nf;
        FINT nc = ngrids_nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *stack = NULL;
        if (cache == NULL) {
                FINT cache_size = int1e_grids_cache_size(envs);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_ALIGN8_INSTACK(gctr, nc * n_comp);

        FINT has_value = CINT1e_grids_loop(gctr, envs, 1, cache);

        FINT counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        if (f_c2s == &c2s_sph_1e_grids) {
                counts[0] = envs->ngrids;
                counts[1] = (envs->i_l*2+1) * x_ctr[0];
                counts[2] = (envs->j_l*2+1) * x_ctr[1];
                counts[3] = 1;
        } else if (f_c2s == &c2s_cart_1e_grids) {
                counts[0] = envs->ngrids;
                counts[1] = envs->nfi * x_ctr[0];
                counts[2] = envs->nfj * x_ctr[1];
                counts[3] = 1;
        }
        FINT nout = dims[0] * dims[1] * dims[2];
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

FINT CINT1e_grids_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs,
                             double *cache, void (*f_c2s)())
{
        FINT ngrids = envs->ngrids;
        if (out == NULL) {
                return int1e_grids_cache_size(envs);
        }
        FINT *x_ctr = envs->x_ctr;
        FINT ngrids_nf = envs->ngrids * envs->nf;
        FINT nc = ngrids_nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *stack = NULL;
        if (cache == NULL) {
                FINT cache_size = int1e_grids_cache_size(envs);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_ALIGN8_INSTACK(gctr, nc * n_comp);

        FINT has_value = CINT1e_grids_loop(gctr, envs, 1, cache);

        FINT counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        counts[0] = envs->ngrids;
        counts[1] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[2] = CINTcgto_spinor(envs->shls[1], envs->bas);
        counts[3] = 1;
        FINT nout = dims[0] * dims[1] * dims[2];
        FINT n;
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}

FINT int1e_grids_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                     FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_grids_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_grids;
        return CINT1e_grids_drv(out, dims, &envs, cache, &c2s_sph_1e_grids);
}

void int1e_grids_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}


FINT int1e_grids_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                      FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int1e_grids_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e_grids;
        return CINT1e_grids_drv(out, dims, &envs, cache, &c2s_cart_1e_grids);
}

FINT int1e_grids_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                        FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        return -1;
}



ALL_CINT(int1e_grids)

