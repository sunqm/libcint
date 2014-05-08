/*
 * File: optimizer.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * optimizer for 2e integrals.  Note if CINT2e_drv is only called a few
 * hundred times, this optimizer cannot really speed up the integration. 
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cint_bas.h"
#include "g2e.h"
#include "optimizer.h"
#include "misc.h"

#define MAX(X,Y) (X)>(Y)?(X):(Y)

// generate caller to CINTinit_2e_optimizer for each type of function
void CINTinit_2e_optimizer(CINTOpt **opt, const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env)
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
void CINTinit_optimizer(CINTOpt **opt, const int *atm, const int natm,
                        const int *bas, const int nbas, const double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
}

void CINTdel_2e_optimizer(CINTOpt **opt)
{
        CINTOpt *opt0 = *opt;
        if (!opt0) { // when opt is created by CINTno_optimizer
                return;
        }

        int i;

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

void CINTno_optimizer(CINTOpt **opt, const int *atm, const int natm,
                      const int *bas, const int nbas, const double *env)
{
        *opt = NULL;
}

void CINTuse_all_optimizer(CINTOpt **opt, int *ng,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_index_xyz(*opt, ng, atm, natm, bas, nbas, env);
}


/* len(ng) = 9. The first 4 items are the increment adding to envs.li_ceil
 * ... envs.ll_ceil for shell i, j, k, l */
void
CINTOpt_set_index_xyz(CINTOpt *opt, int *ng,
                      const int *atm, const int natm,
                      const int *bas, const int nbas, const double *env)
{
        int i, j, k, l, ptr;
        int n = ANG_MAX*ANG_MAX*ANG_MAX*ANG_MAX;
        opt->index_xyz_array = (int **)malloc(sizeof(int *) * n);
        for (i = 0; i < n; i++) {
                opt->index_xyz_array[i] = NULL;
        }

        int max_l = 0;
        for (i = 0; i < nbas; i++) {
                max_l = MAX(max_l, bas(ANG_OF,i));
        }

        int fakebas[BAS_SLOTS*max_l];
        int fakenbas = max_l+1;
        // fakebas only initializes ANG_OF, since the others does not
        // affect index_xyz
        memset(fakebas, 0, sizeof(int)*BAS_SLOTS*fakenbas);
        for (i = 0; i <= max_l; i++) {
                fakebas[BAS_SLOTS*i+ANG_OF] = i;
        }

        CINTEnvVars envs;
        int shls[4];
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
                        (int *)malloc(sizeof(int)*envs.nf*3);
                CINTg2e_index_xyz(opt->index_xyz_array[ptr], &envs);
        } } } }
}


// for the coeffs of the pGTO, find the maximum abs(coeff)
static double max_pgto_coeff(const double *coeff, int nprim, int nctr,
                             int prim_id)
{
        int i;
        double maxc = 0;
        for (i = 0; i < nctr; i++) {
                maxc = MAX(maxc, fabs(coeff[i*nprim+prim_id]));
        }
        return maxc;
}

void CINTOpt_setij(CINTOpt *opt, const int *atm, const int natm,
                   const int *bas, const int nbas, const double *env)
{
        int i, j, ip, jp, io, jo, off;
        if (!opt->prim_offset) {
                opt->prim_offset = (int *)malloc(sizeof(int) * nbas);
                opt->tot_prim = 0;
                for (i = 0; i < nbas; i++) {
                        opt->prim_offset[i] = opt->tot_prim;
                        opt->tot_prim += bas(NPRIM_OF, i);
                }
        }

        int iprim, ictr, jprim, jctr, il, jl;
        double eij, aij, rr, maxci, maxcj, compensation;
        const double *ai, *aj, *ri, *rj, *ci, *cj;
        double *expij, *rij;
        int *cceij;
        opt->expij = (double **)malloc(sizeof(double *) * opt->tot_prim);
        opt->rij = (double **)malloc(sizeof(double *) * opt->tot_prim);
        opt->cceij = (int **)malloc(sizeof(int *) * opt->tot_prim);
        for (i = 0; i < nbas; i++) {
                ri = env + atm(PTR_COORD,bas(ATOM_OF,i));
                ai = env + bas(PTR_EXP,i);
                io = opt->prim_offset[i];
                iprim = bas(NPRIM_OF,i);
                ictr = bas(NCTR_OF,i);
                ci = env + bas(PTR_COEFF,i);
                il = bas(ANG_OF,i);
                for (ip = 0; ip < bas(NPRIM_OF,i); ip++) {
                        maxci = max_pgto_coeff(ci, iprim, ictr, ip);
                        maxci = maxci / CINTgto_norm(il, ai[ip]);
                        expij = (double *)malloc(sizeof(double)*opt->tot_prim);
                        rij = (double *)malloc(sizeof(double)*opt->tot_prim*3);
                        cceij = (int *)malloc(sizeof(int)*opt->tot_prim);
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
                                        maxcj = maxcj/CINTgto_norm(jl, aj[jp]);
                                        aij = ai[ip] + aj[jp];
                                        off = jo + jp;
                                        eij = rr * ai[ip] * aj[jp] / aij;
                                        expij[off] = exp(-eij);
                                        rij[off*3+0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / aij;
                                        rij[off*3+1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / aij;
                                        rij[off*3+2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / aij;
                                        /* TODO: better screen estimation,
                                         * e.g. based on (ab|cd) ~< (00|00)
                                        if (expij[off] > exp(-EXPCUTOFF/2)) {
                                                cceij[off] = 0;
                                        } else {
                                                //cceij[off] = (eij[off] > EXPCUTOFF);
                                                cceij[off] = (expij[off]*maxci*maxcj < exp(-EXPCUTOFF));
                                        } */
                                        /* in experiment */
                                        compensation = exp(il*2+jl*2);
                                        cceij[off] =-log(expij[off]*maxci*maxcj
                                                         *compensation);
                                }
                        }
                }
        }
}

void CINTOpt_set_non0coeff(CINTOpt *opt, const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env)
{
        int i, j, k, ip, io;
        if (!opt->prim_offset) {
                opt->prim_offset = (int *)malloc(sizeof(int) * nbas);
                opt->tot_prim = 0;
                for (i = 0; i < nbas; i++) {
                        opt->prim_offset[i] = opt->tot_prim;
                        opt->tot_prim += bas(NPRIM_OF, i);
                }
        }

        int iprim, ictr;
        const double *ci;
        double *non0coeff;
        int *non0idx;
        opt->non0ctr = (int *)malloc(sizeof(int) * opt->tot_prim);
        opt->non0idx = (int **)malloc(sizeof(int *) * opt->tot_prim);
        opt->non0coeff = (double **)malloc(sizeof(double *) * opt->tot_prim);
        for (i = 0; i < nbas; i++) {
                io = opt->prim_offset[i];
                iprim = bas(NPRIM_OF,i);
                ictr = bas(NCTR_OF,i);
                ci = env + bas(PTR_COEFF,i);
                for (ip = 0; ip < bas(NPRIM_OF,i); ip++) {
                        non0idx = (int *)malloc(sizeof(int) * ictr);
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

