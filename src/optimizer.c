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
        opt0->eij = NULL;
        opt0->rij = NULL;
        opt0->screenij = NULL;
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

        unsigned int i;

        if (opt0->index_xyz_array) {
                for (i = 0; i < ANG_MAX*ANG_MAX*ANG_MAX*ANG_MAX; i++) {
                        if (opt0->index_xyz_array[i]) {
                                free(opt0->index_xyz_array[i]);
                        };
                }
                free(opt0->index_xyz_array);
        }

        if (opt0->prim_offset) {
                free(opt0->prim_offset);
                for (i = 0; i < opt0->tot_prim; i++) {
                        free(opt0->eij[i]);
                        free(opt0->rij[i]);
                        free(opt0->screenij[i]);
                }
                free(opt0->eij);
                free(opt0->rij);
                free(opt0->screenij);
        }

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

void CINTuse_all_optimizer(CINTOpt **opt, unsigned int *ng,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_index_xyz(*opt, ng, atm, natm, bas, nbas, env);
}


/* len(ng) = 9. The first 4 items are the increment adding to envs.li_ceil
 * ... envs.ll_ceil for shell i, j, k, l */
void
CINTOpt_set_index_xyz(CINTOpt *opt, unsigned int *ng,
                      const int *atm, const int natm,
                      const int *bas, const int nbas, const double *env)
{
        unsigned int i, j, k, l, ptr;
        unsigned int n = ANG_MAX*ANG_MAX*ANG_MAX*ANG_MAX;
        opt->index_xyz_array = (unsigned int **)malloc(sizeof(int *) * n);
        for (i = 0; i < n; i++) {
                opt->index_xyz_array[i] = NULL;
        }

        unsigned int max_l = 0;
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
        unsigned int shls[4];
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
                        (unsigned int *)malloc(sizeof(unsigned int)*envs.nf*3);
                CINTg2e_index_xyz(opt->index_xyz_array[ptr], &envs);
        } } } }
}


// for the coeffs of the pGTO, find the maximum abs(coeff)
static double max_pgto_coeff(const double *coeff, int nprim, int nctr,
                             unsigned int prim_id)
{
        unsigned int i;
        double maxc = 0;
        for (i = 0; i < nctr; i++) {
                maxc = MAX(maxc, fabs(coeff[i*nprim+prim_id]));
        }
        return maxc;
}

void CINTOpt_setij(CINTOpt *opt, const int *atm, const int natm,
                   const int *bas, const int nbas, const double *env)
{
        opt->prim_offset = (unsigned int *)malloc(sizeof(unsigned int) * nbas);
        opt->tot_prim = 0;
        unsigned int i, j, ip, jp, io, jo, off;
        for (i = 0; i < nbas; i++) {
                opt->prim_offset[i] = opt->tot_prim;
                opt->tot_prim += bas(NPRIM_OF, i);
        }

        unsigned int iprim, ictr, jprim, jctr, il, jl;
        double aij, rr, maxci, maxcj, logci, logcj;
        const double *ai, *aj, *ri, *rj, *ci, *cj;
        double *eij, *rij;
        int *screenij;
        opt->eij = (double **)malloc(sizeof(double *) * opt->tot_prim);
        opt->rij = (double **)malloc(sizeof(double *) * opt->tot_prim);
        opt->screenij = (int **)malloc(sizeof(int *) * opt->tot_prim);
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
                        logci = log(maxci/CINTgto_norm(il, ai[ip]));
                        eij = (double *)malloc(sizeof(double)*opt->tot_prim);
                        rij = (double *)malloc(sizeof(double)*opt->tot_prim*3);
                        screenij = (int *)malloc(sizeof(int)*opt->tot_prim);
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
                                        logcj = log(maxcj/CINTgto_norm(jl, aj[jp]));
                                        aij = ai[ip] + aj[jp];
                                        off = jo + jp;
                                        eij[off] = rr * ai[ip] * aj[jp] / aij;
                                        rij[off*3+0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / aij;
                                        rij[off*3+1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / aij;
                                        rij[off*3+2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / aij;
                                        /* TODO: better screen estimation,
                                         * e.g. based on (ab|cd) < (00|00) */
                                        if (eij[off] < EXPCUTOFF/2) {
                                                screenij[off] = 0;
                                        } else {
                                                //screenij[off] = (eij[off] > EXPCUTOFF);
                                                //screenij[off] = (eij[off]-log(maxci)-log(maxcj) > EXPCUTOFF);
                                                screenij[off] = (eij[off]-logci-logcj > EXPCUTOFF);
                                        }
                                }
                        }
                        opt->eij[io+ip] = eij;
                        opt->rij[io+ip] = rij;
                        opt->screenij[io+ip] = screenij;
                }
        }
}
