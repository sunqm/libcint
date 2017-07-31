/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Attenuated coulomb operator exp(-w r_{12}^2) / r_{12}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cint_const.h"
#include "cint_bas.h"
#include "rys_roots.h"
#include "g2e.h"

void CINTg0_2e_coulerf(double *g, double fac, CINTEnvVars *envs);

void CINTinit_int2e_coulerf_EnvVars(CINTEnvVars *envs, int *ng, int *shls,
                                    int *atm, int natm, int *bas, int nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        int i_sh = shls[0];
        int j_sh = shls[1];
        int k_sh = shls[2];
        int l_sh = shls[3];
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

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = envs->l_l + ng[LINC];
        envs->nrys_roots =(envs->li_ceil + envs->lj_ceil
                         + envs->lk_ceil + envs->ll_ceil)/2 + 1;

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

        int dli, dlj, dlk, dll;
        int ibase = envs->li_ceil > envs->lj_ceil;
        int kbase = envs->lk_ceil > envs->ll_ceil;
        if (envs->nrys_roots <= 2) { // use the fully optimized lj_4d algorithm
                ibase = 0;
                kbase = 0;
        }
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
                        envs->f_g0_2d4d = &CINTg0_2e_lj2d4d;
                }
        }
        envs->f_g0_2e = &CINTg0_2e_coulerf;
}

/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
void CINTg0_2e_coulerf(double *g, double fac, CINTEnvVars *envs)
{
        double aij = envs->aij;
        double akl = envs->akl;
        double omega = envs->env[PTR_RANGE_OMEGA];
        double a0, a1, fac1, x;
        double u[MXRYSROOTS];
        double *w = g + envs->g_size * 2; // ~ gz
        double rijrkl[3];
        rijrkl[0] = envs->rij[0] - envs->rkl[0];
        rijrkl[1] = envs->rij[1] - envs->rkl[1];
        rijrkl[2] = envs->rij[2] - envs->rkl[2];

        a1 = aij * akl;
        a0 = a1 / (aij + akl);

        double theta = 0;
        if (omega > 0) {
// For long-range part of range-separated Coulomb operator
                theta = omega * omega / (omega * omega + a0);
                a0 *= theta;
        }

        fac1 = sqrt(a0 / (a1 * a1 * a1)) * fac;
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);
        CINTrys_roots(envs->nrys_roots, x, u, w);
        if (envs->g_size == 1) {
                g[0] = 1;
                g[1] = 1;
                g[2] *= fac1;
                return;
        }

        int irys;
        if (omega > 0) {
                /* u[:] = tau^2 / (1 - tau^2)
                 * transform u[:] to theta^-1 tau^2 / (theta^-1 - tau^2)
                 * so the rest code can be reused.
                 */
                for (irys = 0; irys < envs->nrys_roots; irys++) {
                        u[irys] /= u[irys] + 1 - u[irys] * theta;
                }
        }

        double u2, div, tmp1, tmp2, tmp3, tmp4;
        double *rijrx = envs->rijrx;
        double *rklrx = envs->rklrx;
        struct _BC bc;
        double *c00 = bc.c00;
        double *c0p = bc.c0p;

        for (irys = 0; irys < envs->nrys_roots; irys++, c00+=3, c0p+=3) {
                /*
                 *t2 = u(irys)/(1+u(irys))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = a0 * u[irys];
                div = 1 / (u2 * (aij + akl) + a1);
                tmp1 = u2 * div;
                tmp2 = tmp1 * akl;
                tmp3 = tmp1 * aij;
                tmp4 = .5 * div;
                bc.b00[irys] = 0.5 * tmp1;
                bc.b10[irys] = bc.b00[irys] + tmp4 * akl;
                bc.b01[irys] = bc.b00[irys] + tmp4 * aij;
                c00[0] = rijrx[0] - tmp2 * rijrkl[0];
                c00[1] = rijrx[1] - tmp2 * rijrkl[1];
                c00[2] = rijrx[2] - tmp2 * rijrkl[2];
                c0p[0] = rklrx[0] + tmp3 * rijrkl[0];
                c0p[1] = rklrx[1] + tmp3 * rijrkl[1];
                c0p[2] = rklrx[2] + tmp3 * rijrkl[2];
                w[irys] *= fac1;
        }

        (*envs->f_g0_2d4d)(g, &bc, envs);
}
