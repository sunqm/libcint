/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 2-center 2-electron integrals
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint2e.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"

FINT int1e_cache_size(CINTEnvVars *envs);

#define PRIM2CTR0(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, ngp, gp, ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, ngp, gp, ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } \
        } \
        *ctrsymb##empty = 0


FINT CINT2c2e_loop_nopt(double *gctr, CINTEnvVars *envs, double *cache)
{
        FINT *shls  = envs->shls;
        FINT *bas = envs->bas;
        double *env = envs->env;
        FINT i_sh = shls[0];
        FINT k_sh = shls[1];
        FINT i_ctr = envs->x_ctr[0];
        FINT k_ctr = envs->x_ctr[1];
        FINT i_prim = bas(NPRIM_OF, i_sh);
        FINT k_prim = bas(NPRIM_OF, k_sh);
        double *ai = env + bas(PTR_EXP, i_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        FINT n_comp = envs->ncomp_tensor;
        double fac1i, fac1k;
        FINT ip, kp;
        FINT empty[3] = {1, 1, 1};
        FINT *iempty = empty + 0;
        FINT *kempty = empty + 1;
        FINT *gempty = empty + 2;
        /* COMMON_ENVS_AND_DECLARE end */
        const FINT nc = i_ctr * k_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenk = envs->nf * nc * n_comp; // gctrk
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenk + leni + len0;
        double *g;
        MALLOC_INSTACK(g, double, len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrk;

        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctri = gctrk;
                iempty = kempty;
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

        envs->idx = (FINT *)malloc(sizeof(FINT) * envs->nf * 3);
        CINTg1e_index_xyz(envs->idx, envs);

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp]; // to use CINTg0_2e
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        if (i_ctr == 1) {
                                fac1i = fac1k*ci[ip];
                        } else {
                                fac1i = fac1k;
                        }
                        (*envs->f_g0_2e)(g, fac1i, envs);
                        (*envs->f_gout)(gout, g, envs->idx, envs, *gempty);
                        PRIM2CTR0(i, gout, envs->nf*n_comp);
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR0(k, gctri, envs->nf*i_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        free(envs->idx);
        return !*kempty;
}


#define COMMON_ENVS_AND_DECLARE \
        FINT *shls = envs->shls; \
        FINT *bas = envs->bas; \
        double *env = envs->env; \
        FINT i_sh = shls[0]; \
        FINT k_sh = shls[1]; \
        FINT i_ctr = envs->x_ctr[0]; \
        FINT k_ctr = envs->x_ctr[1]; \
        FINT i_prim = bas(NPRIM_OF, i_sh); \
        FINT k_prim = bas(NPRIM_OF, k_sh); \
        double *ai = env + bas(PTR_EXP, i_sh); \
        double *ak = env + bas(PTR_EXP, k_sh); \
        double *ci = env + bas(PTR_COEFF, i_sh); \
        double *ck = env + bas(PTR_COEFF, k_sh); \
        FINT n_comp = envs->ncomp_tensor; \
        double fac1i, fac1k; \
        FINT ip, kp; \
        FINT empty[3] = {1, 1, 1}; \
        FINT *iempty = empty + 0; \
        FINT *kempty = empty + 1; \
        FINT *gempty = empty + 2;

#define USE_OPT \
        FINT off; \
        const FINT io = opt->prim_offset[i_sh]; \
        const FINT ko = opt->prim_offset[k_sh]; \
        envs->idx = opt->index_xyz_array[envs->i_l*LMAX1+envs->k_l]

#define PRIM2CTR(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, ngp, gp, ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } else { \
                        off = ctrsymb##o + ctrsymb##p; \
                        CINTprim_to_ctr_opt(gctr##ctrsymb, ngp, gp, \
                                            opt->non0coeff[off], \
                                            opt->non0idx[off], \
                                            opt->non0ctr[off]); \
                } \
        } \
        *ctrsymb##empty = 0

// i_ctr = k_ctr = 1;
FINT CINT2c2e_11_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt, double *cache)
{
        COMMON_ENVS_AND_DECLARE;
        const FINT nc = 1;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT len0 = envs->nf * n_comp;
        const FINT len = leng + len0;
        double *g;
        MALLOC_INSTACK(g, double, len);
        double *gout;
        if (n_comp == 1) {
                gout = gctr;
        } else {
                gout = g + leng;
        }

        USE_OPT;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        fac1i = fac1k*ci[ip];
                        (*envs->f_g0_2e)(g, fac1i, envs);
                        (*envs->f_gout)(gout, g, envs->idx, envs, *empty);
                        *empty = 0;
                } // end loop i_prim
        } // end loop k_prim

        if (n_comp > 1 && !*empty) {
                CINTdmat_transpose(gctr, gout, envs->nf*nc, n_comp);
        }
        return !*empty;
}

// i_ctr = n; k_ctr = 1;
FINT CINT2c2e_n1_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt, double *cache)
{
        COMMON_ENVS_AND_DECLARE;

        const FINT nc = i_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + leni + len0;
        double *g;
        MALLOC_INSTACK(g, double, len);
        double *g1 = g + leng;
        double *gout, *gctri;
        if (n_comp == 1) {
                gctri = gctr;
        } else {
                gctri = g1;
                g1 += leni;
        }
        gout = g1;

        USE_OPT;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        fac1i = fac1k;
                        (*envs->f_g0_2e)(g, fac1i, envs);
                        (*envs->f_gout)(gout, g, envs->idx, envs, 1);
                        PRIM2CTR(i, gout, envs->nf*n_comp);
                } // end loop i_prim
        } // end loop k_prim

        if (n_comp > 1 && !*iempty) {
                CINTdmat_transpose(gctr, gctri, envs->nf*nc, n_comp);
        }
        return !*iempty;
}

// k_ctr = n; i_ctr = 1;
FINT CINT2c2e_1n_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt, double *cache)
{
        COMMON_ENVS_AND_DECLARE;

        const FINT nc = k_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenk = envs->nf * k_ctr * n_comp; // gctrk
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenk + len0;
        double *g;
        MALLOC_INSTACK(g, double, len);
        double *g1 = g + leng;
        double *gout, *gctrk;
        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        gout = g1;

        USE_OPT;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor;
                *iempty = 1;
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        fac1i = fac1k*ci[ip];
                        (*envs->f_g0_2e)(g, fac1i, envs);
                        (*envs->f_gout)(gout, g, envs->idx, envs, *iempty);
                        *iempty = 0;
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(k, gout,envs->nf*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}


FINT CINT2c2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt, double *cache)
{
        COMMON_ENVS_AND_DECLARE;
        const FINT nc = i_ctr * k_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenk = envs->nf * nc * n_comp; // gctrk
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenk + leni + len0;
        double *g;
        MALLOC_INSTACK(g, double, len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrk;

        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctri = gctrk;
                iempty = kempty;
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

        /* USE_OPT */
        FINT off;
        const FINT io = opt->prim_offset[i_sh];
        const FINT ko = opt->prim_offset[k_sh];
        envs->idx = opt->index_xyz_array[envs->i_l*LMAX1+envs->k_l];
        /* USE_OPT end */

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        if (i_ctr == 1) {
                                fac1i = fac1k*ci[ip];
                        } else {
                                fac1i = fac1k;
                        }
                        (*envs->f_g0_2e)(g, fac1i, envs);
                        (*envs->f_gout)(gout, g, envs->idx, envs, *gempty);
                        PRIM2CTR(i, gout, envs->nf*n_comp);
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(k, gctri, envs->nf*i_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}

static FINT (*CINTf_2c2e_loop[8])() = {
        CINT2c2e_loop,
        CINT2c2e_n1_loop,
        CINT2c2e_1n_loop,
        CINT2c2e_11_loop,
};

FINT CINT2c2e_cart_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                       double *cache)
{
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = leng+len0+nc*n_comp*3;
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = leng+len0+nc*n_comp*3;
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, double, nc*n_comp);

        FINT n;
        FINT has_value;

        if (opt != NULL) {
                n = ((envs->x_ctr[0]==1) << 1) + (envs->x_ctr[1]==1);
                has_value = CINTf_2c2e_loop[n](gctr, envs, opt, cache);
        } else {
                has_value = CINT2c2e_loop_nopt(gctr, envs, cache);
        }

        FINT counts[4];
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfk * x_ctr[1];
        counts[2] = 1;
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_cart_1e(out+nout*n, gctr+nc*n, dims, envs, cache);
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
FINT CINT2c2e_spheric_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                          double *cache)
{
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = MAX(leng+len0+nc*n_comp*3, nc*n_comp+envs->nf*2);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = MAX(leng+len0+nc*n_comp*3, nc*n_comp+envs->nf*2);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, double, nc*n_comp);

        FINT n;
        FINT has_value;

        if (opt != NULL) {
                n = ((envs->x_ctr[0]==1) << 1) + (envs->x_ctr[1]==1);
                has_value = CINTf_2c2e_loop[n](gctr, envs, opt, cache);
        } else {
                has_value = CINT2c2e_loop_nopt(gctr, envs, cache);
        }

        FINT counts[4];
        counts[0] = (envs->i_l*2+1) * x_ctr[0];
        counts[1] = (envs->k_l*2+1) * x_ctr[1];
        counts[2] = 1;
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_sph_1e(out+nout*n, gctr+nc*n, dims, envs, cache);
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
// (spinor|spinor)
FINT CINT2c2e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)())
{
        if (envs->ncomp_e1 > 1 || envs->ncomp_e2 > 1) {
                fprintf(stderr, "CINT2c2e_spinor_drv not implemented\n");
                exit(1);
        }
        if (out == NULL) {
                return int1e_cache_size(envs);
        }
        FINT *x_ctr = envs->x_ctr;
        FINT counts[4];
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        counts[2] = 1;
        counts[3] = 1;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double *stack = NULL;
        if (cache == NULL) {
                FINT cache_size = int1e_cache_size(envs);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, double, nc*n_comp);

        FINT n, has_value;

        if (opt != NULL) {
                n = ((envs->x_ctr[0]==1) << 1) + (envs->x_ctr[1]==1);
                has_value = CINTf_2c2e_loop[n](gctr, envs, opt, cache);
        } else {
                has_value = CINT2c2e_loop_nopt(gctr, envs, cache);
        }

        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        (*f_e1_c2s)(out+nout*n, gctr, dims, envs, cache);
                        gctr += nc;
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


FINT int2c2e_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_spheric_drv(out, dims, &envs, opt, cache);
}
void int2c2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

FINT int2c2e_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                 FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_cart_drv(out, dims, &envs, opt, cache);
}
 
FINT int2c2e_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                   FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_1e);
}


ALL_CINT(int2c2e)
ALL_CINT_FORTRAN_(int2c2e)

