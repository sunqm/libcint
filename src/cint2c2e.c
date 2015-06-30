/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 2-center 2-electron integrals
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "g2e.h"
#include "g3c2e.h"
#include "optimizer.h"
#include "cint2e.h"
#include "misc.h"
#include "fblas.h"
#include "cart2sph.h"
#include "c2f.h"
#include "g2c2e.h"

#define SQUARE(r)       (r)[0]*(r)[0] + (r)[1]*(r)[1] + (r)[2]*(r)[2]

#define PRIM2CTR0(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, ngp, gp, \
                                          envs->ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, ngp, gp, \
                                          envs->ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } \
        } \
        *ctrsymb##empty = 0


FINT CINT2c2e_loop_nopt(double *gctr, CINTEnvVars *envs)
{
        const FINT *shls  = envs->shls;
        const FINT *bas = envs->bas;
        const double *env = envs->env;
        const FINT i_sh = shls[0];
        const FINT k_sh = shls[1];
        const FINT i_ctr  = envs->i_ctr;
        const FINT k_ctr  = envs->k_ctr;
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *ak = env + bas(PTR_EXP, k_sh);
        const double *ci = env + bas(PTR_COEFF, i_sh);
        const double *ck = env + bas(PTR_COEFF, k_sh);
        const FINT n_comp = envs->ncomp_tensor;
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
        double *const g = (double *)malloc(sizeof(double)*len);
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
        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp]; // to use CINTg0_2e
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        if (i_ctr == 1) {
                                fac1i = fac1k*ci[ip];
                        } else {
                                fac1i = fac1k;
                        }
                        CINT2e_core(gout, g, fac1i, envs, *gempty);
                        PRIM2CTR0(i, gout, envs->nf*n_comp);
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR0(k, gctri, envs->nf*i_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        free(g);
        free(envs->idx);
        return !*kempty;
}


#define COMMON_ENVS_AND_DECLARE \
        const FINT *shls = envs->shls; \
        const FINT *bas = envs->bas; \
        const double *env = envs->env; \
        const FINT i_ctr  = envs->i_ctr; \
        const FINT k_ctr  = envs->k_ctr; \
        const FINT i_sh = shls[0]; \
        const FINT k_sh = shls[1]; \
        const double *ai = env + bas(PTR_EXP, i_sh); \
        const double *ak = env + bas(PTR_EXP, k_sh); \
        const double *ci = env + bas(PTR_COEFF, i_sh); \
        const double *ck = env + bas(PTR_COEFF, k_sh); \
        const FINT n_comp = envs->ncomp_tensor; \
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
        envs->idx = opt->index_xyz_array[envs->i_l*ANG_MAX+envs->k_l]

#define PRIM2CTR(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, ngp, gp, \
                                          envs->ctrsymb##_prim, \
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
FINT CINT2c2e_11_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const FINT nc = 1;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT len0 = envs->nf * n_comp;
        const FINT len = leng + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
        double *gout;
        if (n_comp == 1) {
                gout = gctr;
        } else {
                gout = g + leng;
        }

        USE_OPT;

        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        fac1i = fac1k*ci[ip];
                        CINT2e_core(gout, g, fac1i, envs, *empty);
                        *empty = 0;
                } // end loop i_prim
        } // end loop k_prim

        if (n_comp > 1 && !*empty) {
                CINTdmat_transpose(gctr, gout, envs->nf*nc, n_comp);
        }
        free(g);
        return !*empty;
}

// i_ctr = n; k_ctr = 1;
FINT CINT2c2e_n1_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const FINT nc = i_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + leni + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
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

        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        fac1i = fac1k;
                        CINT2e_core(gout, g, fac1i, envs, 1);
                        PRIM2CTR(i, gout, envs->nf*n_comp);
                } // end loop i_prim
        } // end loop k_prim

        if (n_comp > 1 && !*iempty) {
                CINTdmat_transpose(gctr, gctri, envs->nf*nc, n_comp);
        }
        free(g);
        return !*iempty;
}

// k_ctr = n; i_ctr = 1;
FINT CINT2c2e_1n_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const FINT nc = k_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenk = envs->nf * k_ctr * n_comp; // gctrk
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenk + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
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

        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor;
                *iempty = 1;
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        fac1i = fac1k*ci[ip];
                        CINT2e_core(gout, g, fac1i, envs, *iempty);
                        *iempty = 0;
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(k, gout,envs->nf*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        free(g);
        return !*kempty;
}


FINT CINT2c2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const FINT nc = i_ctr * k_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenk = envs->nf * nc * n_comp; // gctrk
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenk + leni + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
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
        envs->idx = opt->index_xyz_array[envs->i_l*ANG_MAX+envs->k_l];
        /* USE_OPT end */

        *kempty = 1;
        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        if (i_ctr == 1) {
                                fac1i = fac1k*ci[ip];
                        } else {
                                fac1i = fac1k;
                        }
                        CINT2e_core(gout, g, fac1i, envs, *gempty);
                        PRIM2CTR(i, gout, envs->nf*n_comp);
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(k, gctri, envs->nf*i_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        free(g);
        return !*kempty;
}

static FINT (*CINTf_2c2e_loop[8])() = {
        CINT2c2e_loop,
        CINT2c2e_n1_loop,
        CINT2c2e_1n_loop,
        CINT2c2e_11_loop,
};

FINT CINT2c2e_cart_drv(double *opij, CINTEnvVars *envs, const CINTOpt *opt)
{
        const FINT ip = CINTcgto_cart(envs->shls[0], envs->bas);
        const FINT kp = CINTcgto_cart(envs->shls[1], envs->bas);
        const FINT nop = ip * kp;
        const FINT nc = envs->nf * envs->i_ctr * envs->k_ctr;
        double *const gctr = malloc(sizeof(double) * nc * envs->ncomp_tensor);
        double *pgctr = gctr;
        FINT n;
        FINT has_value;

        if (opt) {
                n = ((envs->i_ctr==1) << 1) + (envs->k_ctr==1);
                has_value = CINTf_2c2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT2c2e_loop_nopt(gctr, envs);
        }

        if (has_value) {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        c2s_cart_1e(opij, pgctr, envs);
                        opij += nop;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nop * envs->ncomp_tensor, opij);
        }
        free(gctr);
        return has_value;
}
FINT CINT2c2e_spheric_drv(double *opij, CINTEnvVars *envs, const CINTOpt *opt)
{
        const FINT ip = CINTcgto_spheric(envs->shls[0], envs->bas);
        const FINT kp = CINTcgto_spheric(envs->shls[1], envs->bas);
        const FINT nop = ip * kp;
        const FINT nc = envs->nf * envs->i_ctr * envs->k_ctr;
        double *const gctr = malloc(sizeof(double) * nc * envs->ncomp_tensor);
        double *pgctr = gctr;
        FINT n;
        FINT has_value;

        if (opt) {
                n = ((envs->i_ctr==1) << 1) + (envs->k_ctr==1);
                has_value = CINTf_2c2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT2c2e_loop_nopt(gctr, envs);
        }

        if (has_value) {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        c2s_sph_1e(opij, pgctr, envs);
                        opij += nop;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nop * envs->ncomp_tensor, opij);
        }
        free(gctr);
        return has_value;
}


FINT cint2c2e_sph(double *opij, const FINT *shls,
                 const FINT *atm, const FINT natm,
                 const FINT *bas, const FINT nbas, const double *env,
                 const CINTOpt *opt)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_spheric_drv(opij, &envs, opt);
}
void cint2c2e_sph_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                          const FINT *bas, const FINT nbas, const double *env)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_2cindex_xyz(*opt, ng, atm, natm, bas, nbas, env);
}

FINT cint2c2e_cart(double *opij, const FINT *shls,
                  const FINT *atm, const FINT natm,
                  const FINT *bas, const FINT nbas, const double *env,
                  const CINTOpt *opt)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_cart_drv(opij, &envs, opt);
}
void cint2c2e_cart_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                           const FINT *bas, const FINT nbas, const double *env)
{
        cint2c2e_sph_optimizer(opt, atm, natm, bas, nbas, env);
}


/*
 * * * * * * * * * * * * * * * * * * * * *
 * c to fortran interface
 */

C2Fo_(cint2c2e_cart);
C2Fo_(cint2c2e_sph);
OPTIMIZER2F_(cint2c2e_cart_optimizer);
OPTIMIZER2F_(cint2c2e_sph_optimizer);

