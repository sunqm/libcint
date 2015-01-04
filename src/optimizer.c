/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * optimizer for 2e integrals.  Note if CINT2e_drv is only called a few
 * hundred times, this optimizer cannot really speed up the integration. 
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cint_const.h"
#include "cint_bas.h"
#include "g2e.h"
#include "optimizer.h"
#include "misc.h"

#define MAX(X,Y) (X)>(Y)?(X):(Y)

// generate caller to CINTinit_2e_optimizer for each type of function
void CINTinit_2e_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                           const FINT *bas, const FINT nbas, const double *env)
{
        CINTOpt *opt0 = (CINTOpt *)malloc(sizeof(CINTOpt));
        opt0->index_xyz_array = NULL;
        opt0->prim_offset = NULL;
        opt0->non0ctr = NULL;
        opt0->non0idx = NULL;
        opt0->non0coeff = NULL;
        opt0->expij = NULL;
        opt0->rij = NULL;
        opt0->cceij = NULL;
        opt0->tot_prim = 0;
        *opt = opt0;
}
void CINTinit_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                        const FINT *bas, const FINT nbas, const double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
}

void CINTdel_2e_optimizer(CINTOpt **opt)
{
        CINTOpt *opt0 = *opt;
        if (!opt0) { // when opt is created by CINTno_optimizer
                return;
        }

        FINT i;

        if (opt0->index_xyz_array) {
                for (i = 0; i < ANG_MAX*ANG_MAX*ANG_MAX*ANG_MAX; i++) {
                        if (opt0->index_xyz_array[i]) {
                                free(opt0->index_xyz_array[i]);
                        };
                }
                free(opt0->index_xyz_array);
        }

        if (opt0->expij) {
                for (i = 0; i < opt0->tot_prim; i++) {
                        free(opt0->expij[i]);
                        free(opt0->rij[i]);
                        free(opt0->cceij[i]);
                }
                free(opt0->expij);
                free(opt0->rij);
                free(opt0->cceij);
        }

        if (opt0->non0ctr) {
                free(opt0->non0ctr);
                for (i = 0; i < opt0->tot_prim; i++) {
                        free(opt0->non0idx[i]);
                        free(opt0->non0coeff[i]);
                }
                free(opt0->non0idx);
                free(opt0->non0coeff);
        }

        if (opt0->prim_offset) {
                free(opt0->prim_offset);
        }

        opt0->tot_prim = 0;

        free(opt0);
        *opt = NULL;
}
void CINTdel_optimizer(CINTOpt **opt)
{
        CINTdel_2e_optimizer(opt);
}

void CINTno_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                      const FINT *bas, const FINT nbas, const double *env)
{
        *opt = NULL;
}

void CINTuse_all_optimizer(CINTOpt **opt, FINT *ng,
                           const FINT *atm, const FINT natm,
                           const FINT *bas, const FINT nbas, const double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_index_xyz(*opt, ng, atm, natm, bas, nbas, env);
}

/* len(ng) = 8. The first 4 items are the increment adding to envs.li_ceil
 * ... envs.ll_ceil for shell i, j, k, l */
void
CINTOpt_set_index_xyz(CINTOpt *opt, FINT *ng,
                      const FINT *atm, const FINT natm,
                      const FINT *bas, const FINT nbas, const double *env)
{
        FINT i, j, k, l, ptr;
        FINT n = ANG_MAX*ANG_MAX*ANG_MAX*ANG_MAX;
        opt->index_xyz_array = (FINT **)malloc(sizeof(FINT *) * n);
        for (i = 0; i < n; i++) {
                opt->index_xyz_array[i] = NULL;
        }

        FINT max_l = 0;
        for (i = 0; i < nbas; i++) {
                max_l = MAX(max_l, bas(ANG_OF,i));
        }

        FINT fakenbas = max_l+1;
        FINT fakebas[BAS_SLOTS*fakenbas];
        memset(fakebas, 0, sizeof(FINT)*BAS_SLOTS*fakenbas);
        // fakebas only initializes ANG_OF, since the others does not
        // affect index_xyz
        for (i = 0; i <= max_l; i++) {
                fakebas[BAS_SLOTS*i+ANG_OF] = i;
        }

        CINTEnvVars envs;
        FINT shls[4];
        for (i = 0; i <= max_l; i++) {
        for (j = 0; j <= max_l; j++) {
        for (k = 0; k <= max_l; k++) {
        for (l = 0; l <= max_l; l++) {
                shls[0] = i; shls[1] = j; shls[2] = k; shls[3] = l;
                CINTinit_int2e_EnvVars(&envs, ng, shls,
                                       atm, natm, fakebas, fakenbas, env);

                ptr = i*ANG_MAX*ANG_MAX*ANG_MAX
                    + j*ANG_MAX*ANG_MAX
                    + k*ANG_MAX
                    + l;
                opt->index_xyz_array[ptr] = 
                        (FINT *)malloc(sizeof(FINT)*envs.nf*3);
                CINTg2e_index_xyz(opt->index_xyz_array[ptr], &envs);
        } } } }
}


// for the coeffs of the pGTO, find the maximum abs(coeff)
static double max_pgto_coeff(const double *coeff, FINT nprim, FINT nctr,
                             FINT prim_id)
{
        FINT i;
        double maxc = 0;
        for (i = 0; i < nctr; i++) {
                maxc = MAX(maxc, fabs(coeff[i*nprim+prim_id]));
        }
        return maxc;
}

void CINTOpt_setij(CINTOpt *opt, FINT *ng,
                   const FINT *atm, const FINT natm,
                   const FINT *bas, const FINT nbas, const double *env)
{
        FINT i, j, ip, jp, io, jo, off;
        if (!opt->prim_offset) {
                opt->prim_offset = (FINT *)malloc(sizeof(FINT) * nbas);
                opt->tot_prim = 0;
                for (i = 0; i < nbas; i++) {
                        opt->prim_offset[i] = opt->tot_prim;
                        opt->tot_prim += bas(NPRIM_OF, i);
                }
        }

        FINT ik_inc, jl_inc;
        if ((ng[IINC]+ng[JINC]) > (ng[KINC]+ng[LINC])) {
                ik_inc = ng[IINC];
                jl_inc = ng[JINC];
        } else {
                ik_inc = ng[KINC];
                jl_inc = ng[LINC];
        }

        FINT iprim, ictr, jprim, jctr, il, jl;
        double eij, aij, rr, maxci, maxcj, rirj_g4d;
        const double *ai, *aj, *ri, *rj, *ci, *cj;
        double *expij, *rij;
        FINT *cceij;
        opt->expij = (double **)malloc(sizeof(double *) * opt->tot_prim);
        opt->rij = (double **)malloc(sizeof(double *) * opt->tot_prim);
        opt->cceij = (FINT **)malloc(sizeof(FINT *) * opt->tot_prim);
        for (i = 0; i < nbas; i++) {
                ri = env + atm(PTR_COORD,bas(ATOM_OF,i));
                ai = env + bas(PTR_EXP,i);
                io = opt->prim_offset[i];
                iprim = bas(NPRIM_OF,i);
                ictr = bas(NCTR_OF,i);
                ci = env + bas(PTR_COEFF,i);
// For derivative/dipole operator, the l-value in g2e is virtually increased
                il = bas(ANG_OF,i) + ik_inc;
                for (ip = 0; ip < bas(NPRIM_OF,i); ip++) {
                        maxci = max_pgto_coeff(ci, iprim, ictr, ip);
                        maxci = maxci / CINTgto_norm(il, ai[ip]);
                        expij = (double *)malloc(sizeof(double)*opt->tot_prim);
                        rij = (double *)malloc(sizeof(double)*opt->tot_prim*3);
                        cceij = (FINT *)malloc(sizeof(FINT) * opt->tot_prim);
                        opt->expij[io+ip] = expij;
                        opt->rij[io+ip] = rij;
                        opt->cceij[io+ip] = cceij;

                        for (j = 0; j < nbas; j++) {
                                rj = env + atm(PTR_COORD,bas(ATOM_OF,j));
                                aj = env + bas(PTR_EXP,j);
                                jo = opt->prim_offset[j];
                                jprim = bas(NPRIM_OF,j);
                                jctr = bas(NCTR_OF,j);
                                cj = env + bas(PTR_COEFF,j);
                                jl = bas(ANG_OF,j) + jl_inc;
                                rr = (ri[0]-rj[0])*(ri[0]-rj[0])
                                   + (ri[1]-rj[1])*(ri[1]-rj[1])
                                   + (ri[2]-rj[2])*(ri[2]-rj[2]);
                                for (jp = 0; jp < bas(NPRIM_OF,j); jp++) {
                                        maxcj = max_pgto_coeff(cj, jprim, jctr, jp);
                                        maxcj = maxcj/CINTgto_norm(jl, aj[jp]);
                                        aij = ai[ip] + aj[jp];
                                        off = jo + jp;
                                        eij = rr * ai[ip] * aj[jp] / aij;
                                        expij[off] = exp(-eij);
                                        rij[off*3+0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / aij;
                                        rij[off*3+1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / aij;
                                        rij[off*3+2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / aij;
/* estimation of the value, based on g0_2e_2d and g0_xx2d_4d,
 * value ~< exp(-eij)*(il+jl+2)!*(aij/2)^(il+jl)*(ri_or_rj-rij)^(ij+jl)*rirj^max(il,jl)
 *       ~< *       exp(-eij)*(il+jl+2)!*(aij/2)^(il+jl)*rirj^((il+jl)+max(il,jl))
 * But in practice, rirj^((il+jl)/2) is usually large enough to cover all other factors */
/* rr+1 to prevent log() diverge when i,j on same center */
                                        rirj_g4d = pow((rr+1), (il+jl+1)/2);
/*cceij[off] =-log(expij[off]*maxci*maxcj*rirj_g4d);
  when eij is big, expij == 0, singular value in cceij */
                                        cceij[off] = eij - log(maxci*maxcj*rirj_g4d);
                                }
                        }
                }
        }
}

void CINTOpt_set_non0coeff(CINTOpt *opt, const FINT *atm, const FINT natm,
                           const FINT *bas, const FINT nbas, const double *env)
{
        FINT i, j, k, ip, io;
        if (!opt->prim_offset) {
                opt->prim_offset = (FINT *)malloc(sizeof(FINT) * nbas);
                opt->tot_prim = 0;
                for (i = 0; i < nbas; i++) {
                        opt->prim_offset[i] = opt->tot_prim;
                        opt->tot_prim += bas(NPRIM_OF, i);
                }
        }

        FINT iprim, ictr;
        const double *ci;
        double *non0coeff;
        FINT *non0idx;
        opt->non0ctr = (FINT *)malloc(sizeof(FINT) * opt->tot_prim);
        opt->non0idx = (FINT **)malloc(sizeof(FINT *) * opt->tot_prim);
        opt->non0coeff = (double **)malloc(sizeof(double *) * opt->tot_prim);
        for (i = 0; i < nbas; i++) {
                io = opt->prim_offset[i];
                iprim = bas(NPRIM_OF,i);
                ictr = bas(NCTR_OF,i);
                ci = env + bas(PTR_COEFF,i);
                for (ip = 0; ip < bas(NPRIM_OF,i); ip++) {
                        non0idx = (FINT *)malloc(sizeof(FINT) * ictr);
                        non0coeff = (double *)malloc(sizeof(double) * ictr);
                        opt->non0idx[io+ip] = non0idx;
                        opt->non0coeff[io+ip] = non0coeff;

                        for (j = 0, k = 0; j < ictr; j++) {
                                if (ci[iprim*j+ip] != 0) {
                                        non0coeff[k] = ci[iprim*j+ip];
                                        non0idx[k] = j;
                                        k++;
                                }
                        }
                        opt->non0ctr[io+ip] = k;
                }
        }
}

