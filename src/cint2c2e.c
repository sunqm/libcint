/*
 * File: cint2e.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO integrals
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
#include "cart2sph.h"
#include "c2f.h"

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
        const FINT j_sh = shls[1];
        const FINT i_ctr  = envs->i_ctr;
        const FINT j_ctr  = envs->j_ctr;
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        const double *ci = env + bas(PTR_COEFF, i_sh);
        const double *cj = env + bas(PTR_COEFF, j_sh);
        double fac1i, fac1j;
        FINT ip, jp;
        FINT empty[3] = {1, 1, 1};
        FINT *iempty = empty + 0;
        FINT *jempty = empty + 1;
        FINT *gempty = empty + 2;
        /* COMMON_ENVS_AND_DECLARE end */
        const FINT nc = i_ctr * j_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenj = envs->nf * nc; // gctrj
        const FINT leni = envs->nf * i_ctr; // gctri
        const FINT len0 = envs->nf; // gout
        const FINT len = leng + lenj + leni + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj;

        gctrj = gctr;
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

        envs->idx = (FINT *)malloc(sizeof(FINT) * envs->nf * 3);
        CINTg1e_index_xyz(envs->idx, envs);

        *jempty = 1;
        for (jp = 0; jp < envs->j_prim; jp++) {
                envs->akl = aj[jp]; // to use CINTg0_2e
                if (j_ctr == 1) {
                        fac1j = envs->common_factor * cj[jp];
                } else {
                        fac1j = envs->common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->aij = ai[ip];
                        if (i_ctr == 1) {
                                fac1i = fac1j*ci[ip];
                        } else {
                                fac1i = fac1j;
                        }
                        CINT2e_core(gout, g, fac1i, envs, *gempty);
                        PRIM2CTR0(i, gout, envs->nf);
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR0(j, gctri, envs->nf*i_ctr);
                }
        } // end loop j_prim
        free(g);
        free(envs->idx);
        return !*jempty;
}


#define COMMON_ENVS_AND_DECLARE \
        const FINT *shls = envs->shls; \
        const FINT *bas = envs->bas; \
        const double *env = envs->env; \
        const FINT i_ctr  = envs->i_ctr; \
        const FINT j_ctr  = envs->j_ctr; \
        const FINT i_sh = shls[0]; \
        const FINT j_sh = shls[1]; \
        const double *ai = env + bas(PTR_EXP, i_sh); \
        const double *aj = env + bas(PTR_EXP, j_sh); \
        const double *ci = env + bas(PTR_COEFF, i_sh); \
        const double *cj = env + bas(PTR_COEFF, j_sh); \
        double fac1i, fac1j; \
        FINT ip, jp; \
        FINT empty[3] = {1, 1, 1}; \
        FINT *iempty = empty + 0; \
        FINT *jempty = empty + 1; \
        FINT *gempty = empty + 2;

#define USE_OPT \
        FINT off; \
        const FINT io = opt->prim_offset[i_sh]; \
        const FINT jo = opt->prim_offset[j_sh]; \
        envs->idx = opt->index_xyz_array[envs->i_l*ANG_MAX+envs->j_l]

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

// i_ctr = j_ctr = 1;
FINT CINT2c2e_11_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT len0 = envs->nf;
        const FINT len = leng + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
        double *gout = gctr;

        USE_OPT;

        for (jp = 0; jp < envs->j_prim; jp++) {
                envs->akl = aj[jp];
                fac1j = envs->common_factor * cj[jp];
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->aij = ai[ip];
                        fac1i = fac1j*ci[ip];
                        CINT2e_core(gout, g, fac1i, envs, *empty);
                        *empty = 0;
                } // end loop i_prim
        } // end loop j_prim

        free(g);
        return !(*empty);
}

// i_ctr = n; j_ctr = 1;
FINT CINT2c2e_n1_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT leni = envs->nf * i_ctr; // gctri
        const FINT len0 = envs->nf; // gout
        const FINT len = leng + leni + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
        double *g1 = g + leng;
        double *gout, *gctri;

        gctri = gctr;

        gout = g1;

        USE_OPT;

        for (jp = 0; jp < envs->j_prim; jp++) {
                envs->akl = aj[jp];
                fac1j = envs->common_factor * cj[jp];
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->aij = ai[ip];
                        fac1i = fac1j;
                        CINT2e_core(gout, g, fac1i, envs, 1);
                        PRIM2CTR(i, gout,envs->nf);
                } // end loop i_prim
        } // end loop j_prim

        free(g);
        return !(*iempty);
}

// j_ctr = n; i_ctr = 1;
FINT CINT2c2e_1n_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenj = envs->nf * j_ctr; // gctrj
        const FINT len0 = envs->nf; // gout
        const FINT len = leng + lenj + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
        double *g1 = g + leng;
        double *gout, *gctrj;

        gctrj = gctr;
        gout = g1;

        USE_OPT;

        for (jp = 0; jp < envs->j_prim; jp++) {
                envs->akl = aj[jp];
                fac1j = envs->common_factor;
                *iempty = 1;
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->aij = ai[ip];
                        fac1i = fac1j*ci[ip];
                        CINT2e_core(gout, g, fac1i, envs, *iempty);
                        *iempty = 0;
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(j, gout,envs->nf);
                }
        } // end loop j_prim

        free(g);
        return !(*jempty);
}


FINT CINT2c2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const FINT nc = i_ctr * j_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenj = envs->nf * nc; // gctrj
        const FINT leni = envs->nf * i_ctr; // gctri
        const FINT len0 = envs->nf; // gout
        const FINT len = leng + lenj + leni + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj;

        gctrj = gctr;
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

        /* USE_OPT */
        FINT off;
        const FINT io = opt->prim_offset[i_sh];
        const FINT jo = opt->prim_offset[j_sh];
        envs->idx = opt->index_xyz_array[envs->i_l*ANG_MAX+envs->j_l];
        /* USE_OPT end */

        *jempty = 1;
        for (jp = 0; jp < envs->j_prim; jp++) {
                envs->akl = aj[jp];
                if (j_ctr == 1) {
                        fac1j = envs->common_factor * cj[jp];
                } else {
                        fac1j = envs->common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->aij = ai[ip];
                        if (i_ctr == 1) {
                                fac1i = fac1j*ci[ip];
                        } else {
                                fac1i = fac1j;
                        }
                        CINT2e_core(gout, g, fac1i, envs, *gempty);
                        PRIM2CTR(i, gout, envs->nf);
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(j, gctri, envs->nf*i_ctr);
                }
        } // end loop j_prim

        free(g);
        return !(*jempty);
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
        const FINT jp = CINTcgto_cart(envs->shls[1], envs->bas);
        const FINT nop = ip * jp;
        const FINT nc = envs->nf * envs->i_ctr * envs->j_ctr;
        double *const gctr = malloc(sizeof(double) * nc * envs->ncomp_tensor);
        double *pgctr = gctr;
        FINT n;
        FINT has_value;

        if (opt) {
                n = ((envs->i_ctr==1) << 1) + (envs->j_ctr==1);
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
        const FINT jp = CINTcgto_spheric(envs->shls[1], envs->bas);
        const FINT nop = ip * jp;
        const FINT nc = envs->nf * envs->i_ctr * envs->j_ctr;
        double *const gctr = malloc(sizeof(double) * nc * envs->ncomp_tensor);
        double *pgctr = gctr;
        FINT n;
        FINT has_value;

        if (opt) {
                n = ((envs->i_ctr==1) << 1) + (envs->j_ctr==1);
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

