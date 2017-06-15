/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cint_bas.h"
#include "misc.h"
#include "g2e.h"

void CINTinit_int2c2e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                              FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        const FINT i_sh = shls[0];
        const FINT k_sh = shls[1];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = 0;
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = 0;
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, k_sh);
        envs->x_ctr[2] = 1;
        envs->x_ctr[3] = 1;
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = 1;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nfl = 1;
        envs->nf = envs->nfi * envs->nfk;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));

        envs->common_factor = (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->k_l);

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = 0;
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = 0;
        envs->nrys_roots =(envs->li_ceil + envs->lk_ceil)/2 + 1;

        FINT dli = envs->li_ceil + 1;
        FINT dlk = envs->lk_ceil + 1;
        envs->g_stride_i = envs->nrys_roots;
        envs->g_stride_k = envs->nrys_roots * dli;
        envs->g_stride_l = envs->g_stride_k;
        envs->g_size     = envs->nrys_roots * dli * dlk;

        envs->aj = 0;
        envs->al = 0;
        envs->rij[0] = envs->ri[0];
        envs->rij[1] = envs->ri[1];
        envs->rij[2] = envs->ri[2];
        envs->rkl[0] = envs->rk[0];
        envs->rkl[1] = envs->rk[1];
        envs->rkl[2] = envs->rk[2];
        envs->g2d_ijmax = envs->g_stride_i;
        envs->g2d_klmax = envs->g_stride_k;
        envs->rkrl[0] = envs->rk[0];
        envs->rkrl[1] = envs->rk[1];
        envs->rkrl[2] = envs->rk[2];
        envs->rklrx[0] = 0;
        envs->rklrx[1] = 0;
        envs->rklrx[2] = 0;
        envs->rirj[0] = envs->ri[0];
        envs->rirj[1] = envs->ri[1];
        envs->rirj[2] = envs->ri[2];
        envs->rijrx[0] = 0;
        envs->rijrx[1] = 0;
        envs->rijrx[2] = 0;
        envs->rx_in_rklrx = envs->rk;
        envs->rx_in_rijrx = envs->ri;

        envs->f_g0_2d4d = &CINTg0_2e_2d;
        envs->f_g0_2e = &CINTg0_2e;

// initialize j_l, j_ctr, nfj because they are used in c2s_sph_1e and
// CINTg1e_index_xyz
        envs->j_l = envs->k_l;
        envs->nfj = envs->nfk;
        envs->g_stride_j = envs->g_stride_k;
}

