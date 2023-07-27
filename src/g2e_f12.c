/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Yukawa potential and Slater-type geminal
 */


#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "config.h"
#include "cint_bas.h"
#include "rys_roots.h"
#include "misc.h"
#include "g2e.h"

FINT CINTg0_2e_stg(double *g, double *rij, double *rkl, double cutoff, CINTEnvVars *envs);
FINT CINTg0_2e_yp(double *g, double *rij, double *rkl, double cutoff, CINTEnvVars *envs);
void CINTg0_2e_stg_lj2d4d(double *g, struct _BC *bc, const CINTEnvVars *envs);

void CINTinit_int2e_yp_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                               FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        FINT i_sh = shls[0];
        FINT j_sh = shls[1];
        FINT k_sh = shls[2];
        FINT l_sh = shls[3];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = bas(ANG_OF, l_sh);
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->x_ctr[2] = bas(NCTR_OF, k_sh);
        envs->x_ctr[3] = bas(NCTR_OF, l_sh);
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nfl = (envs->l_l+1)*(envs->l_l+2)/2;
        envs->nf = envs->nfi * envs->nfk * envs->nfl * envs->nfj;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));
        envs->rl = env + atm(PTR_COORD, bas(ATOM_OF, l_sh));

        envs->common_factor = (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l) * CINTcommon_fac_sp(envs->l_l);
        if (env[PTR_EXPCUTOFF] == 0) {
                envs->expcutoff = EXPCUTOFF;
        } else {
                envs->expcutoff = MAX(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
        }

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = envs->l_l + ng[LINC];

        // ceil(L_tot/2) + 1
        FINT nroots = (envs->li_ceil + envs->lj_ceil +
                      envs->lk_ceil + envs->ll_ceil + 3)/2;
        envs->nrys_roots = nroots;
        assert(nroots < MXRYSROOTS);

        FINT dli, dlj, dlk, dll;
        FINT ibase = envs->li_ceil > envs->lj_ceil;
        FINT kbase = envs->lk_ceil > envs->ll_ceil;
        if (kbase) {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_ik2d4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_kj2d4d;
                }
        } else {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_il2d4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_stg_lj2d4d;
                }
        }
        envs->f_g0_2e = &CINTg0_2e_yp;

        if (kbase) {
                dlk = envs->lk_ceil + envs->ll_ceil + 1;
                dll = envs->ll_ceil + 1;
        } else {
                dlk = envs->lk_ceil + 1;
                dll = envs->lk_ceil + envs->ll_ceil + 1;
        }

        if (ibase) {
                dli = envs->li_ceil + envs->lj_ceil + 1;
                dlj = envs->lj_ceil + 1;
        } else {
                dli = envs->li_ceil + 1;
                dlj = envs->li_ceil + envs->lj_ceil + 1;
        }
        envs->g_stride_i = nroots;
        envs->g_stride_k = nroots * dli;
        envs->g_stride_l = nroots * dli * dlk;
        envs->g_stride_j = nroots * dli * dlk * dll;
        envs->g_size     = nroots * dli * dlk * dll * dlj;

        if (kbase) {
                envs->g2d_klmax = envs->g_stride_k;
                envs->rx_in_rklrx = envs->rk;
                envs->rkrl[0] = envs->rk[0] - envs->rl[0];
                envs->rkrl[1] = envs->rk[1] - envs->rl[1];
                envs->rkrl[2] = envs->rk[2] - envs->rl[2];
        } else {
                envs->g2d_klmax = envs->g_stride_l;
                envs->rx_in_rklrx = envs->rl;
                envs->rkrl[0] = envs->rl[0] - envs->rk[0];
                envs->rkrl[1] = envs->rl[1] - envs->rk[1];
                envs->rkrl[2] = envs->rl[2] - envs->rk[2];
        }

        if (ibase) {
                envs->g2d_ijmax = envs->g_stride_i;
                envs->rx_in_rijrx = envs->ri;
                envs->rirj[0] = envs->ri[0] - envs->rj[0];
                envs->rirj[1] = envs->ri[1] - envs->rj[1];
                envs->rirj[2] = envs->ri[2] - envs->rj[2];
        } else {
                envs->g2d_ijmax = envs->g_stride_j;
                envs->rx_in_rijrx = envs->rj;
                envs->rirj[0] = envs->rj[0] - envs->ri[0];
                envs->rirj[1] = envs->rj[1] - envs->ri[1];
                envs->rirj[2] = envs->rj[2] - envs->ri[2];
        }
}

void CINTinit_int2e_stg_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_int2e_yp_EnvVars(envs, ng, shls, atm, natm, bas, nbas, env);
        envs->f_g0_2e = &CINTg0_2e_stg;
}


void CINTg0_2e_stg_lj2d4d(double *g, struct _BC *bc, const CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_lj2d_4d(g, envs);
}

FINT CINTg0_2e_yp(double *g, double *rij, double *rkl, double cutoff, CINTEnvVars *envs)
{
        double aij, akl, a0, a1, fac1, x;
        double ua = 0;
        double rijrkl[3];
        double u[MXRYSROOTS];
        double *w = g + envs->g_size * 2; // ~ gz
        double zeta = envs->env[PTR_F12_ZETA];
        FINT nroots = envs->nrys_roots;
        FINT i;

        aij = envs->ai[0] + envs->aj[0];
        akl = envs->ak[0] + envs->al[0];
        a1 = aij * akl;
        a0 = a1 / (aij + akl);
        //fac1 = sqrt(a0 / (a1 * a1 * a1)) * envs->fac[idsimd];
        fac1 = envs->fac[0] / (sqrt(aij+akl) * a1);

        //:ua[k] = zeta*zeta / a0[k];
        if (zeta > 0) {
                ua = .25 * zeta * zeta / a0;
        }

        rijrkl[0] = rij[0] - rkl[0];
        rijrkl[1] = rij[1] - rkl[1];
        rijrkl[2] = rij[2] - rkl[2];
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);
        if (zeta > 0) {
                CINTstg_roots(nroots, x, ua, u, w);
        } else {
                CINTrys_roots(nroots, x, u, w);
        }

        if (zeta > 0) {
                //:w *= t;
                //:u -> t/(1-t);
                for (i = 0; i < nroots; i++) {
                        w[i] *= u[i];
                        u[i] = u[i] / (1 - u[i]);
                }
        }

        if (envs->g_size == 1) {
                g[0] = 1;
                g[1] = 1;
                g[2] *= fac1;
                return 1;
        }

        double u2, div, tmp1, tmp2, tmp3, tmp4;
        double rijrx[3];
        double rklrx[3];
        rijrx[0] = envs->rij[0] - envs->rx_in_rijrx[0];
        rijrx[1] = envs->rij[1] - envs->rx_in_rijrx[1];
        rijrx[2] = envs->rij[2] - envs->rx_in_rijrx[2];
        rklrx[0] = envs->rkl[0] - envs->rx_in_rklrx[0];
        rklrx[1] = envs->rkl[1] - envs->rx_in_rklrx[1];
        rklrx[2] = envs->rkl[2] - envs->rx_in_rklrx[2];
        struct _BC bc;
        double *c00 = bc.c00;
        double *c0p = bc.c0p;

        for (i = 0; i < nroots; i++, c00+=3, c0p+=3) {
                /*
                 *u(i) = t2/(1-t2)
                 *t2 = u(i)/(1+u(i))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = a0 * u[i];
                div = 1 / (u2 * (aij + akl) + a1);
                tmp1 = u2 * div;
                tmp2 = tmp1 * akl;
                tmp3 = tmp1 * aij;
                tmp4 = .5 * div;
                bc.b00[i] = 0.5 * tmp1;
                bc.b10[i] = bc.b00[i] + tmp4 * akl;
                bc.b01[i] = bc.b00[i] + tmp4 * aij;
                c00[0] = rijrx[0] - tmp2 * rijrkl[0];
                c00[1] = rijrx[1] - tmp2 * rijrkl[1];
                c00[2] = rijrx[2] - tmp2 * rijrkl[2];
                c0p[0] = rklrx[0] + tmp3 * rijrkl[0];
                c0p[1] = rklrx[1] + tmp3 * rijrkl[1];
                c0p[2] = rklrx[2] + tmp3 * rijrkl[2];
                w[i] *= fac1;
        }

        (*envs->f_g0_2d4d)(g, &bc, envs);
        return 1;
}


FINT CINTg0_2e_stg(double *g, double *rij, double *rkl, double cutoff, CINTEnvVars *envs)
{
        double aij, akl, a0, a1, fac1, x;
        double ua = 0;
        double rijrkl[3];
        double u[MXRYSROOTS];
        double *w = g + envs->g_size * 2; // ~ gz
        double zeta = envs->env[PTR_F12_ZETA];
        FINT nroots = envs->nrys_roots;
        FINT i;

        aij = envs->ai[0] + envs->aj[0];
        akl = envs->ak[0] + envs->al[0];
        a1 = aij * akl;
        a0 = a1 / (aij + akl);
        //fac1 = sqrt(a0 / (a1 * a1 * a1)) * fac;
        fac1 = envs->fac[0] / (sqrt(aij+akl) * a1);

        //:ua[k] = zeta*zeta / a0[k];
        if (zeta > 0) {
                ua = .25 * zeta * zeta / a0;
        }

        rijrkl[0] = rij[0] - rkl[0];
        rijrkl[1] = rij[1] - rkl[1];
        rijrkl[2] = rij[2] - rkl[2];
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);
        if (zeta > 0) {
                CINTstg_roots(nroots, x, ua, u, w);
        } else {
                CINTrys_roots(nroots, x, u, w);
        }

        if (zeta > 0) {
                //:w *= (1-t) * 2*ua/zeta;
                //:u -> t/(1-t);
                double ua2 = 2. * ua / zeta;
                for (i = 0; i < nroots; i++) {
                        w[i] *= (1-u[i]) * ua2;
                        u[i] = u[i] / (1 - u[i]);
                }
        }

        if (envs->g_size == 1) {
                g[0] = 1;
                g[1] = 1;
                g[2] *= fac1;
                return 1;
        }

        double u2, div, tmp1, tmp2, tmp3, tmp4;
        double rijrx[3];
        double rklrx[3];
        rijrx[0] = envs->rij[0] - envs->rx_in_rijrx[0];
        rijrx[1] = envs->rij[1] - envs->rx_in_rijrx[1];
        rijrx[2] = envs->rij[2] - envs->rx_in_rijrx[2];
        rklrx[0] = envs->rkl[0] - envs->rx_in_rklrx[0];
        rklrx[1] = envs->rkl[1] - envs->rx_in_rklrx[1];
        rklrx[2] = envs->rkl[2] - envs->rx_in_rklrx[2];
        struct _BC bc;
        double *c00 = bc.c00;
        double *c0p = bc.c0p;

        for (i = 0; i < nroots; i++, c00+=3, c0p+=3) {
                /*
                 *u(i) = t2/(1-t2)
                 *t2 = u(i)/(1+u(i))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = a0 * u[i];
                div = 1 / (u2 * (aij + akl) + a1);
                tmp1 = u2 * div;
                tmp2 = tmp1 * akl;
                tmp3 = tmp1 * aij;
                tmp4 = .5 * div;
                bc.b00[i] = 0.5 * tmp1;
                bc.b10[i] = bc.b00[i] + tmp4 * akl;
                bc.b01[i] = bc.b00[i] + tmp4 * aij;
                c00[0] = rijrx[0] - tmp2 * rijrkl[0];
                c00[1] = rijrx[1] - tmp2 * rijrkl[1];
                c00[2] = rijrx[2] - tmp2 * rijrkl[2];
                c0p[0] = rklrx[0] + tmp3 * rijrkl[0];
                c0p[1] = rklrx[1] + tmp3 * rijrkl[1];
                c0p[2] = rklrx[2] + tmp3 * rijrkl[2];
                w[i] *= fac1;
        }

        (*envs->f_g0_2d4d)(g, &bc, envs);
        return 1;
}
