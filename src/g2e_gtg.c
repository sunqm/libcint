/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "config.h"
#include "cint_bas.h"
#include "misc.h"
#include "g2e.h"

void CINTg0_2e_lj2d4d_regular(double *g, struct _BC *bc, const CINTEnvVars *envs);
FINT CINTg0_2e_gtg(double *g, double *rij, double *rkl, double cutoff, CINTEnvVars *envs);

void CINTinit_int2e_gtg_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT k_sh = shls[2];
        const FINT l_sh = shls[3];
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

        envs->common_factor = (M_PI*M_PI*M_PI)
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
        envs->nrys_roots = 1;

        assert(i_sh < SHLS_MAX);
        assert(j_sh < SHLS_MAX);
        assert(k_sh < SHLS_MAX);
        assert(l_sh < SHLS_MAX);
        assert(envs->i_l < ANG_MAX);
        assert(envs->j_l < ANG_MAX);
        assert(envs->k_l < ANG_MAX);
        assert(envs->l_l < ANG_MAX);
        assert(bas(ATOM_OF,i_sh) >= 0);
        assert(bas(ATOM_OF,j_sh) >= 0);
        assert(bas(ATOM_OF,k_sh) >= 0);
        assert(bas(ATOM_OF,l_sh) >= 0);
        assert(bas(ATOM_OF,i_sh) < natm);
        assert(bas(ATOM_OF,j_sh) < natm);
        assert(bas(ATOM_OF,k_sh) < natm);
        assert(bas(ATOM_OF,l_sh) < natm);
        assert(envs->nrys_roots < MXRYSROOTS);

        FINT dli, dlj, dlk, dll;
        FINT ibase = envs->li_ceil > envs->lj_ceil;
        FINT kbase = envs->lk_ceil > envs->ll_ceil;

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
        envs->g_stride_i = envs->nrys_roots;
        envs->g_stride_k = envs->nrys_roots * dli;
        envs->g_stride_l = envs->nrys_roots * dli * dlk;
        envs->g_stride_j = envs->nrys_roots * dli * dlk * dll;
        envs->g_size     = envs->nrys_roots * dli * dlk * dll * dlj;

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
                        envs->f_g0_2d4d = &CINTg0_2e_lj2d4d_regular;
                }
        }
        envs->f_g0_2e = &CINTg0_2e_gtg;
}


void CINTg0_2e_lj2d4d_regular(double *g, struct _BC *bc, const CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_lj2d_4d(g, envs);
        return;
}

/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
FINT CINTg0_2e_gtg(double *g, double *rij, double *rkl, double cutoff, CINTEnvVars *envs)
{
        double aij = envs->ai[0] + envs->aj[0];
        double akl = envs->ak[0] + envs->al[0];
        double *env = envs->env;
        double zeta = env[PTR_GTG_ZETA];
        double a0, a1, fac1, x, t;
        double *gz = g + envs->g_size * 2;
        double rijrkl[3];
        rijrkl[0] = envs->rij[0] - envs->rkl[0];
        rijrkl[1] = envs->rij[1] - envs->rkl[1];
        rijrkl[2] = envs->rij[2] - envs->rkl[2];

        a1 = aij * akl;
        a0 = a1 / (aij + akl);
        t = zeta / (zeta + a0);
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);
        fac1 = (1-t) / a1;
        gz[0] = fac1*sqrt(fac1) * exp(-t * x) * envs->fac[0];
        if (envs->g_size == 1) {
                g[0] = 1;
                g[1] = 1;
                return 1;
        }

        double div, tmp1, tmp2, tmp3, tmp4;
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
        double *b00 = bc.b00;
        double *b10 = bc.b10;
        double *b01 = bc.b01;

        div = 1 / (zeta * (aij + akl) + a1);
        tmp1 = zeta * div;
        tmp2 = tmp1 * akl;
        tmp3 = tmp1 * aij;
        tmp4 = .5 * div;
        b00[0] = 0.5 * tmp1;
        b10[0] = b00[0] + tmp4 * akl;
        b01[0] = b00[0] + tmp4 * aij;
        c00[0] = rijrx[0] - tmp2 * rijrkl[0];
        c00[1] = rijrx[1] - tmp2 * rijrkl[1];
        c00[2] = rijrx[2] - tmp2 * rijrkl[2];
        c0p[0] = rklrx[0] + tmp3 * rijrkl[0];
        c0p[1] = rklrx[1] + tmp3 * rijrkl[1];
        c0p[2] = rklrx[2] + tmp3 * rijrkl[2];

        (*envs->f_g0_2d4d)(g, &bc, envs);
        return 1;
}

