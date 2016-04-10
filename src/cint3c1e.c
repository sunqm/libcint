/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 3-center 1-electron integrals
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "g1e.h"
#include "optimizer.h"
#include "misc.h"
#include "fblas.h"
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

void CINTg3c1e_ovlp(double *g, double ai, double aj, double ak,
                    double fac, const CINTEnvVars *envs);
void CINTg3c1e_index_xyz(FINT *idx, const CINTEnvVars *envs);
FINT CINTinit_int3c1e_EnvVars(CINTEnvVars *envs, const FINT *ng, const FINT *shls,
                             const FINT *atm, const FINT natm,
                             const FINT *bas, const FINT nbas, const double *env);


FINT CINT3c1e_loop_nopt(double *gctr, CINTEnvVars *envs)
{
        const FINT *shls  = envs->shls;
        const FINT *bas = envs->bas;
        const double *env = envs->env;
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT k_sh = shls[2];
        const FINT i_ctr  = envs->i_ctr;
        const FINT j_ctr  = envs->j_ctr;
        const FINT k_ctr  = envs->k_ctr;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *rk = envs->rk;
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        const double *ak = env + bas(PTR_EXP, k_sh);
        const double *ci = env + bas(PTR_COEFF, i_sh);
        const double *cj = env + bas(PTR_COEFF, j_sh);
        const double *ck = env + bas(PTR_COEFF, k_sh);
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k;
        FINT ip, jp, kp;
        FINT empty[4] = {1, 1, 1, 1};
        FINT *iempty = empty + 0;
        FINT *jempty = empty + 1;
        FINT *kempty = empty + 2;
        FINT *gempty = empty + 3;
        /* COMMON_ENVS_AND_DECLARE end */
        const FINT nc = i_ctr * j_ctr * k_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenk = envs->nf * nc * n_comp; // gctrk
        const FINT lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenk + lenj + leni + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk;

        if (n_comp == 1) {
                gctrk = gctr;
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

        double eijk, dijk, aijk;
        double aiajrr, aiakrr, ajakrr;
        double rirk[3];
        double rjrk[3];
        rirk[0] = ri[0] - rk[0];
        rirk[1] = ri[1] - rk[1];
        rirk[2] = ri[2] - rk[2];
        rjrk[0] = rj[0] - rk[0];
        rjrk[1] = rj[1] - rk[1];
        rjrk[2] = rj[2] - rk[2];
        const double rr_ij = SQUARE(envs->rirj);
        const double rr_ik = SQUARE(      rirk);
        const double rr_jk = SQUARE(      rjrk);
        envs->idx = (FINT *)malloc(sizeof(FINT) * envs->nf * 3);
        CINTg3c1e_index_xyz(envs->idx, envs);

        *kempty = 1;
        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *jempty = 1;
                }

                for (jp = 0; jp < envs->j_prim; jp++) {
                        envs->aj = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        ajakrr = aj[jp] * ak[kp] * rr_jk;
                        for (ip = 0; ip < envs->i_prim; ip++) {
                                envs->ai = ai[ip];
                                aijk = ai[ip] + aj[jp] + ak[kp];
                                aiakrr = ai[ip] * ak[kp] * rr_ik;
                                aiajrr = ai[ip] * aj[jp] * rr_ij;
                                eijk = (aiajrr+aiakrr+ajakrr) / aijk;
                                if (eijk > EXPCUTOFF) {
                                        continue;
                                }

                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*exp(-eijk);
                                } else {
                                        fac1i = fac1j*exp(-eijk);
                                }
                                dijk = fac1i / (aijk * sqrt(aijk));
                                CINTg3c1e_ovlp(g, ai[ip], aj[jp], ak[kp], dijk, envs);
                                (*envs->f_gout)(g, gout, envs->idx, envs, *gempty);

                                PRIM2CTR0(i, gout, envs->nf*n_comp);
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR0(j, gctri, envs->nf*i_ctr*n_comp);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR0(k, gctrj,envs->nf*i_ctr*j_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        free(g);
        free(envs->idx);
        return !*kempty;
}


FINT CINT3c1e_cart_drv(double *opijk, CINTEnvVars *envs, const CINTOpt *opt)
{
        const FINT nc = envs->nf * envs->i_ctr * envs->j_ctr * envs->k_ctr;
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *const gctr = malloc(sizeof(double) * nc * n_comp);
        FINT n;
        FINT has_value;
        double *pgctr = gctr;

        //if (opt != NULL) {
        //        n = ((envs->i_ctr==1) << 2) + ((envs->j_ctr==1) << 1)
        //          + (envs->k_ctr==1);
        //        //has_value = CINTf_3c1e_loop[n](gctr, envs, opt);
        //} else {
                has_value = CINT3c1e_loop_nopt(gctr, envs);
        //}

        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_cart_3c1e(opijk, pgctr, envs);
                        opijk += nc;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nc * n_comp, opijk);
        }
        free(gctr);
        return has_value;
}
FINT CINT3c1e_spheric_drv(double *opijk, CINTEnvVars *envs, const CINTOpt *opt,
                         void (*const f_e1_c2s)(), FINT is_ssc)
{
        const FINT ip = CINTcgto_spheric(envs->shls[0], envs->bas);
        const FINT jp = CINTcgto_spheric(envs->shls[1], envs->bas);
        const FINT kp = CINTcgto_spheric(envs->shls[2], envs->bas);
        const FINT nop = ip * jp * kp;
        const FINT nc = envs->nf * envs->i_ctr * envs->j_ctr * envs->k_ctr;
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *const gctr = malloc(sizeof(double) * nc * n_comp);
        double *pgctr = gctr;
        FINT n;
        FINT has_value;

        //if (opt != NULL) {
        //        n = ((envs->i_ctr==1) << 2) + ((envs->j_ctr==1) << 1)
        //          + (envs->k_ctr==1);
        //        //has_value = CINTf_3c1e_loop[n](gctr, envs, opt);
        //} else {
                has_value = CINT3c1e_loop_nopt(gctr, envs);
        //}

        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        (*f_e1_c2s)(opijk, pgctr, envs);
                        opijk += nop;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nop * n_comp, opijk);
        }
        free(gctr);
        return has_value;
}

void CINTgout3c1e(double *g, double *gout, const FINT *idx,
                     const CINTEnvVars *envs, FINT gout_empty)
{
        FINT ix, iy, iz, n;

        if (gout_empty) {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        gout[n] = g[ix] * g[iy] * g[iz];
                }
        } else {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        gout[n] += g[ix] * g[iy] * g[iz];
                }
        }
}

FINT cint3c1e_sph(double *opijk, const FINT *shls,
                 const FINT *atm, const FINT natm,
                 const FINT *bas, const FINT nbas, const double *env,
                 const CINTOpt *opt)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_spheric_drv(opijk, &envs, opt, &c2s_sph_3c1e, 0);
}
void cint3c1e_sph_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                          const FINT *bas, const FINT nbas, const double *env)
{
        *opt = NULL;
}

FINT cint3c1e_cart(double *opijk, const FINT *shls,
                  const FINT *atm, const FINT natm,
                  const FINT *bas, const FINT nbas, const double *env,
                  const CINTOpt *opt)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_cart_drv(opijk, &envs, opt);
}
void cint3c1e_cart_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                           const FINT *bas, const FINT nbas, const double *env)
{
        cint3c1e_sph_optimizer(opt, atm, natm, bas, nbas, env);
}

/*
 * * * * * * * * * * * * * * * * * * * * *
 * c to fortran interface
 */

C2Fo_(cint3c1e_cart);
C2Fo_(cint3c1e_sph);
OPTIMIZER2F_(cint3c1e_cart_optimizer);
OPTIMIZER2F_(cint3c1e_sph_optimizer);

