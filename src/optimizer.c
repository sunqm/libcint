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
#include "g3c1e.h"
#include "optimizer.h"
#include "misc.h"

// generate caller to CINTinit_2e_optimizer for each type of function
void CINTinit_2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env)
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
void CINTinit_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                        FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
}

void CINTdel_2e_optimizer(CINTOpt **opt)
{
        CINTOpt *opt0 = *opt;
        if (opt0 == NULL) { // when opt is created by CINTno_optimizer
                return;
        }

        FINT i;

        if (opt0->index_xyz_array != NULL) {
                free(opt0->index_xyz_array[0]);
                free(opt0->index_xyz_array);
        }

        if (opt0->expij != NULL) {
                for (i = 0; i < opt0->tot_prim; i++) {
                        free(opt0->expij[i]);
                        free(opt0->rij[i]);
                        free(opt0->cceij[i]);
                }
                free(opt0->expij);
                free(opt0->rij);
                free(opt0->cceij);
        }

        if (opt0->non0ctr != NULL) {
                free(opt0->non0ctr);
                for (i = 0; i < opt0->tot_prim; i++) {
                        free(opt0->non0idx[i]);
                        free(opt0->non0coeff[i]);
                }
                free(opt0->non0idx);
                free(opt0->non0coeff);
        }

        if (opt0->prim_offset != NULL) {
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

void CINTno_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                      FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}

static FINT _make_fakebas(FINT *fakebas, FINT *bas, FINT nbas, double *env)
{
        FINT i;
        FINT max_l = 0;
        for (i = 0; i < nbas; i++) {
                max_l = MAX(max_l, bas(ANG_OF,i));
        }

        FINT fakenbas = max_l + 1;
        for (i = 0; i < BAS_SLOTS*fakenbas; i++) {
                fakebas[i] = 0;
        }
        // fakebas only initializes ANG_OF, since the others does not
        // affect index_xyz
        for (i = 0; i <= max_l; i++) {
                fakebas[BAS_SLOTS*i+ANG_OF] = i;
        }
        return max_l;
}
static FINT *_allocate_index_xyz(CINTOpt *opt, FINT max_l, FINT order)
{
        FINT i;
        FINT cumcart = (max_l+1) * (max_l+2) * (max_l+3) / 6;
        FINT ll = max_l + 1;
        FINT cc = cumcart;
        for (i = 1; i < order; i++) {
                ll *= LMAX1;
                cc *= cumcart;
        }
        FINT *buf = malloc(sizeof(FINT) * cc * 3);
        FINT **ppbuf = malloc(sizeof(FINT*) * ll);
        ppbuf[0] = buf;
        for (i = 1; i < ll; i++) {
                ppbuf[i] = NULL;
        }
        opt->index_xyz_array = ppbuf;
        return buf;
}
static void gen_idx(CINTOpt *opt, void (*finit)(), void (*findex_xyz)(),
                    FINT order, FINT max_l, FINT *ng,
                    FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        FINT i, j, k, l, ptr;
        FINT fakebas[BAS_SLOTS*LMAX1];
        FINT max_l1 = _make_fakebas(fakebas, bas, nbas, env);
        if (max_l == 0) {
                max_l = max_l1;
        } else {
                max_l = MIN(max_l, max_l1);
        }
        FINT fakenbas = max_l+1;
        FINT *buf = _allocate_index_xyz(opt, max_l, order);

        CINTEnvVars envs;
        FINT shls[4];
        if (order == 2) {
                for (i = 0; i <= max_l; i++) {
                for (j = 0; j <= max_l; j++) {
                        shls[0] = i; shls[1] = j;
                        (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                        ptr = i*LMAX1 + j;
                        opt->index_xyz_array[ptr] = buf;
                        (*findex_xyz)(opt->index_xyz_array[ptr], &envs);
                        buf += envs.nf * 3;
                } }

        } else if (order == 3) {
                for (i = 0; i <= max_l; i++) {
                for (j = 0; j <= max_l; j++) {
                for (k = 0; k <= max_l; k++) {
                        shls[0] = i; shls[1] = j; shls[2] = k;
                        (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                        ptr = i*LMAX1*LMAX1 + j*LMAX1 + k;
                        opt->index_xyz_array[ptr] = buf;
                        (*findex_xyz)(opt->index_xyz_array[ptr], &envs);
                        buf += envs.nf * 3;
                } } }

        } else {
                for (i = 0; i <= max_l; i++) {
                for (j = 0; j <= max_l; j++) {
                for (k = 0; k <= max_l; k++) {
                for (l = 0; l <= max_l; l++) {
                        shls[0] = i; shls[1] = j; shls[2] = k; shls[3] = l;
                        (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                        ptr = i*LMAX1*LMAX1*LMAX1
                            + j*LMAX1*LMAX1
                            + k*LMAX1
                            + l;
                        opt->index_xyz_array[ptr] = buf;
                        (*findex_xyz)(opt->index_xyz_array[ptr], &envs);
                        buf += envs.nf * 3;
                } } } }
        }
}

void CINTall_1e_optimizer(CINTOpt **opt, FINT *ng,
                          FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int1e_EnvVars, &CINTg1e_index_xyz,
                2, 0, ng, atm, natm, bas, nbas, env);
}

void CINTall_2e_optimizer(CINTOpt **opt, FINT *ng,
                          FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int2e_EnvVars, &CINTg2e_index_xyz,
                4, 0, ng, atm, natm, bas, nbas, env);
}

void CINTall_3c2e_optimizer(CINTOpt **opt, FINT *ng,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int3c2e_EnvVars, &CINTg2e_index_xyz,
                3, 0, ng, atm, natm, bas, nbas, env);
}

void CINTall_2c2e_optimizer(CINTOpt **opt, FINT *ng,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int2c2e_EnvVars, &CINTg1e_index_xyz,
                2, 0, ng, atm, natm, bas, nbas, env);
}

void CINTg3c1e_index_xyz(FINT *idx, const CINTEnvVars *envs);
void CINTall_3c1e_optimizer(CINTOpt **opt, FINT *ng,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int3c1e_EnvVars, &CINTg3c1e_index_xyz,
                3, 0, ng, atm, natm, bas, nbas, env);
}

#ifdef WITH_F12
void CINTinit_int2e_stg_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                           FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTall_2e_stg_optimizer(CINTOpt **opt, FINT *ng,
                              FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int2e_stg_EnvVars, &CINTg2e_index_xyz,
                4, 0, ng, atm, natm, bas, nbas, env);
}
#endif


// for the coeffs of the pGTO, find the maximum abs(coeff)
static double max_pgto_coeff(double *coeff, FINT nprim, FINT nctr,
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
                   FINT *atm, FINT natm,
                   FINT *bas, FINT nbas, double *env)
{
        FINT i, j, ip, jp, io, jo, off;
        if (opt->prim_offset == NULL) {
                opt->prim_offset = (FINT *)malloc(sizeof(FINT) * nbas);
                opt->tot_prim = 0;
                for (i = 0; i < nbas; i++) {
                        opt->prim_offset[i] = opt->tot_prim;
                        opt->tot_prim += bas(NPRIM_OF, i);
                }
        }

        FINT ijkl_inc;
        if ((ng[IINC]+ng[JINC]) > (ng[KINC]+ng[LINC])) {
                ijkl_inc = ng[IINC] + ng[JINC];
        } else {
                ijkl_inc = ng[KINC] + ng[LINC];
        }

        FINT iprim, ictr, jprim, jctr, il, jl;
        double eij, aij, rr, maxci, maxcj, rirj_g4d;
        double *ai, *aj, *ri, *rj, *ci, *cj;
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
                il = bas(ANG_OF,i);
                for (ip = 0; ip < bas(NPRIM_OF,i); ip++) {
                        maxci = max_pgto_coeff(ci, iprim, ictr, ip);
                        maxci = maxci / CINTgto_norm(bas(ANG_OF,i), ai[ip]);
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
                                jl = bas(ANG_OF,j);
                                rr = (ri[0]-rj[0])*(ri[0]-rj[0])
                                   + (ri[1]-rj[1])*(ri[1]-rj[1])
                                   + (ri[2]-rj[2])*(ri[2]-rj[2]);
                                for (jp = 0; jp < bas(NPRIM_OF,j); jp++) {
                                        maxcj = max_pgto_coeff(cj, jprim, jctr, jp);
                                        maxcj = maxcj / CINTgto_norm(bas(ANG_OF,j), aj[jp]);
                                        aij = ai[ip] + aj[jp];
                                        off = jo + jp;
                                        eij = rr * ai[ip] * aj[jp] / aij;
                                        expij[off] = exp(-eij);
                                        rij[off*3+0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / aij;
                                        rij[off*3+1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / aij;
                                        rij[off*3+2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / aij;

        if (maxci*maxcj == 0) {
                cceij[off] = 750;
        } else if (rr > 1e-12) {
/* value estimation based on g0_2e_2d and g0_xx2d_4d,
 * value/exp(-eij) ~< (il+jl+2)!*(aij/2)^(il+jl)*(ri_or_rj-rij)^(ij+jl)*rirj^max(il,jl)
 *                 ~< (il+jl+2)!*(aij/2)^(il+jl)*|rirj|^((il+jl)+max(il,jl))
 * But in practice, |rirj|^((il+jl)/2) is large enough to cover all other factors */
                rirj_g4d = pow(rr+0.5, (il+jl+ijkl_inc+1)/2);
                cceij[off] = eij - log(maxci*maxcj*rirj_g4d);
        } else {
/* If basis on the same center, include the (ss|ss)^{1/2} contribution
 * (ss|ss) = 2\sqrt{aij/pi} */
                cceij[off] = -log(maxci*maxcj) - log(aij)/4;
        }
                                }
                        }
                }
        }
}

void CINTOpt_set_non0coeff(CINTOpt *opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env)
{
        FINT i, j, k, ip, io;
        if (opt->prim_offset == NULL) {
                opt->prim_offset = (FINT *)malloc(sizeof(FINT) * nbas);
                opt->tot_prim = 0;
                for (i = 0; i < nbas; i++) {
                        opt->prim_offset[i] = opt->tot_prim;
                        opt->tot_prim += bas(NPRIM_OF, i);
                }
        }

        FINT iprim, ictr;
        double *ci;
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

void CINTdel_pairdata_optimizer(CINTOpt *cintopt)
{
        if (cintopt != NULL && cintopt->expij != NULL) {
                int i;
                for (i = 0; i < cintopt->tot_prim; i++) {
                        free(cintopt->expij[i]);
                        free(cintopt->rij[i]);
                        free(cintopt->cceij[i]);
                }
                free(cintopt->expij);
                free(cintopt->rij);
                free(cintopt->cceij);
                cintopt->expij = NULL;
                cintopt->rij = NULL;
                cintopt->cceij = NULL;
        }
}
