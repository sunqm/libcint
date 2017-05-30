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

#define SQUARE(r)       ((r)[0]*(r)[0] + (r)[1]*(r)[1] + (r)[2]*(r)[2])

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

void CINTg4c1e_ovlp(double *g, double fac, const CINTEnvVars *envs);
void CINTg4c1e_index_xyz(FINT *idx, const CINTEnvVars *envs);
FINT CINTinit_int4c1e_EnvVars(CINTEnvVars *envs, const FINT *ng, const FINT *shls,
                             const FINT *atm, const FINT natm,
                             const FINT *bas, const FINT nbas, const double *env);
void CINTgout3c1e(double *g, double *gout, const FINT *idx,
                  const CINTEnvVars *envs, FINT gout_empty);

FINT CINT4c1e_loop_nopt(double *gctr, CINTEnvVars *envs)
{
        const FINT *shls  = envs->shls;
        const FINT *bas = envs->bas;
        const double *env = envs->env;
        const FINT i_ctr  = envs->i_ctr;
        const FINT j_ctr  = envs->j_ctr;
        const FINT k_ctr  = envs->k_ctr;
        const FINT l_ctr  = envs->l_ctr;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *rk = envs->rk;
        const double *rl = envs->rl;
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT k_sh = shls[2];
        const FINT l_sh = shls[3];
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        const double *ak = env + bas(PTR_EXP, k_sh);
        const double *al = env + bas(PTR_EXP, l_sh);
        const double *ci = env + bas(PTR_COEFF, i_sh);
        const double *cj = env + bas(PTR_COEFF, j_sh);
        const double *ck = env + bas(PTR_COEFF, k_sh);
        const double *cl = env + bas(PTR_COEFF, l_sh);
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k, fac1l;
        FINT ip, jp, kp, lp;
        FINT empty[5] = {1, 1, 1, 1, 1};
        FINT *iempty = empty + 0;
        FINT *jempty = empty + 1;
        FINT *kempty = empty + 2;
        FINT *lempty = empty + 3;
        FINT *gempty = empty + 4;
        /* COMMON_ENVS_AND_DECLARE end */
        const FINT nc = i_ctr * j_ctr * k_ctr * l_ctr;
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenl = envs->nf * nc * n_comp; // gctrl
        const FINT lenk = envs->nf * i_ctr * j_ctr * k_ctr * n_comp; // gctrk
        const FINT lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenl + lenk + lenj + leni + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
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

        double eij, ekl, expijkl, eijkl, a0;
        double rijrkl[3];
        const double dist_ij = SQUARE(envs->rirj);
        const double dist_kl = SQUARE(envs->rkrl);
        envs->idx = (FINT *)malloc(sizeof(FINT) * envs->nf * 3);
        CINTg4c1e_index_xyz(envs->idx, envs);
        double common_fac = envs->common_factor * SQRTPI * M_PI *
                CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l) *
                CINTcommon_fac_sp(envs->k_l) * CINTcommon_fac_sp(envs->l_l);

        *lempty = 1;
        for (lp = 0; lp < envs->l_prim; lp++) {
                envs->al = al[lp];
                if (l_ctr == 1) {
                        fac1l = common_fac * cl[lp];
                } else {
                        fac1l = common_fac;
                        *kempty = 1;
                }
                for (kp = 0; kp < envs->k_prim; kp++) {
                        envs->ak = ak[kp];
                        envs->akl = ak[kp] + al[lp];
                        ekl = dist_kl * ak[kp] * al[lp] / envs->akl;
                        if (ekl > EXPCUTOFF) {
                                goto k_contracted;
                        }
                        envs->rkl[0] = (ak[kp]*rk[0] + al[lp]*rl[0]) / envs->akl;
                        envs->rkl[1] = (ak[kp]*rk[1] + al[lp]*rl[1]) / envs->akl;
                        envs->rkl[2] = (ak[kp]*rk[2] + al[lp]*rl[2]) / envs->akl;
                        envs->rklrx[0] = envs->rkl[0] - envs->rx_in_rklrx[0];
                        envs->rklrx[1] = envs->rkl[1] - envs->rx_in_rklrx[1];
                        envs->rklrx[2] = envs->rkl[2] - envs->rx_in_rklrx[2];
                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                        } else {
                                fac1k = fac1l;
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
                                for (ip = 0; ip < envs->i_prim; ip++) {
                                        envs->ai = ai[ip];
                                        envs->aij = ai[ip] + aj[jp];
                                        eij = dist_ij * ai[ip] * aj[jp] / envs->aij;
                                        if (eij+ekl > EXPCUTOFF) {
                                                goto i_contracted;
                                        }
                                        envs->rij[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / envs->aij;
                                        envs->rij[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / envs->aij;
                                        envs->rij[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / envs->aij;
                                        envs->rijrx[0] = envs->rij[0] - envs->rx_in_rijrx[0];
                                        envs->rijrx[1] = envs->rij[1] - envs->rx_in_rijrx[1];
                                        envs->rijrx[2] = envs->rij[2] - envs->rx_in_rijrx[2];
                                        rijrkl[0] = envs->rij[0] - envs->rkl[0];
                                        rijrkl[1] = envs->rij[1] - envs->rkl[1];
                                        rijrkl[2] = envs->rij[2] - envs->rkl[2];
                                        a0 = envs->aij * envs->akl / (envs->aij + envs->akl);
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
                                        CINTg4c1e_ovlp(g, fac1i, envs);
                                        (*envs->f_gout)(g, gout, envs->idx, envs, *gempty);
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
        free(g);
        free(envs->idx);
        return !*lempty;
}


FINT CINT4c1e_cart_drv(double *out, CINTEnvVars *envs, const CINTOpt *opt)
{
        const FINT nc = envs->nf * envs->i_ctr * envs->j_ctr
                                 * envs->k_ctr * envs->l_ctr;
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double *const gctr = malloc(sizeof(double) * nc * n_comp);
        FINT n;
        FINT has_value;
        double *pgctr = gctr;

        has_value = CINT4c1e_loop_nopt(gctr, envs);

        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_cart_2e1(out, pgctr, envs);
                        out += nc;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nc * n_comp, out);
        }
        free(gctr);
        return has_value;
}
FINT CINT4c1e_spheric_drv(double *out, CINTEnvVars *envs, const CINTOpt *opt)
{
        const FINT ip = CINTcgto_spheric(envs->shls[0], envs->bas);
        const FINT jp = CINTcgto_spheric(envs->shls[1], envs->bas);
        const FINT kp = CINTcgto_spheric(envs->shls[2], envs->bas);
        const FINT lp = CINTcgto_spheric(envs->shls[3], envs->bas);
        const FINT nop = ip * jp * kp * lp;
        const FINT nc = envs->nf * envs->i_ctr * envs->j_ctr
                                 * envs->k_ctr * envs->l_ctr;
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double *const gctr = malloc(sizeof(double) * nc * n_comp);
        double *pgctr = gctr;
        FINT n;
        FINT has_value;

        has_value = CINT4c1e_loop_nopt(gctr, envs);

        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_sph_2e1(out, pgctr, envs);
                        out += nop;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nop * n_comp, out);
        }
        free(gctr);
        return has_value;
}

FINT cint4c1e_ovlp_sph(double *out, const FINT *shls,
                       const FINT *atm, const FINT natm,
                       const FINT *bas, const FINT nbas, const double *env,
                       const CINTOpt *opt)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int4c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT4c1e_spheric_drv(out, &envs, opt);
}
void cint4c1e_ovlp_sph_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                                 const FINT *bas, const FINT nbas, const double *env)
{
        *opt = NULL;
}

FINT cint4c1e_ovlp_cart(double *out, const FINT *shls,
                        const FINT *atm, const FINT natm,
                        const FINT *bas, const FINT nbas, const double *env,
                        const CINTOpt *opt)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int4c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT4c1e_cart_drv(out, &envs, opt);
}
void cint4c1e_ovlp_cart_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                                  const FINT *bas, const FINT nbas, const double *env)
{
        cint4c1e_ovlp_sph_optimizer(opt, atm, natm, bas, nbas, env);
}

/*
 * * * * * * * * * * * * * * * * * * * * *
 * c to fortran interface
 */

C2Fo_(cint4c1e_ovlp_cart);
C2Fo_(cint4c1e_ovlp_sph);
OPTIMIZER2F_(cint4c1e_ovlp_cart_optimizer);
OPTIMIZER2F_(cint4c1e_ovlp_sph_optimizer);

