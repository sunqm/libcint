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

void
CINTOpt_set_log_coeff(CINTOpt *opt, const int *atm, const int natm,
                      const int *bas, const int nbas, const double *env);

// generate caller to CINTinit_2e_optimizer for each type of function
void CINTinit_2e_optimizer(CINTOpt **opt, const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env)
{
        CINTOpt *opt0 = malloc(sizeof(CINTOpt));
        opt0->expcutoff = 46; // ~ 1e-20
        opt0->ptr_log_coeff = NULL;
        opt0->log_coeff = NULL;
        opt0->index_xyz_array = NULL;
        *opt = opt0;
}

void CINTdel_2e_optimizer(CINTOpt *opt)
{
        if (!opt) { // when opt is created by CINTno_optimizer
                return;
        }

        unsigned int i;

        if (opt->ptr_log_coeff) {
                free(opt->ptr_log_coeff);
                free(opt->log_coeff);
        }

        if (opt->index_xyz_array) {
                for (i = 0; i < ANG_MAX*ANG_MAX*ANG_MAX*ANG_MAX; i++) {
                        if (opt->index_xyz_array[i]) {
                                free(opt->index_xyz_array[i]);
                        };
                }
        }

        free(opt);
}

void CINTno_optimizer(CINTOpt **opt, const int *atm, const int natm,
                      const int *bas, const int nbas, const double *env)
{
        *opt = NULL;
}

void CINTcutoff_optimizer(CINTOpt **opt, const int *atm, const int natm,
                          const int *bas, const int nbas, const double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_set_log_coeff(*opt, atm, natm, bas, nbas, env);
}

void CINTuse_all_optimizer(CINTOpt **opt, unsigned int *ng,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env)
{
      CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
      CINTOpt_set_log_coeff(*opt, atm, natm, bas, nbas, env);
      CINTOpt_set_index_xyz(*opt, ng, atm, natm, bas, nbas, env);
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

void
CINTOpt_set_log_coeff(CINTOpt *opt, const int *atm, const int natm,
                      const int *bas, const int nbas, const double *env)
{
        opt->ptr_log_coeff = (unsigned int *)malloc(sizeof(unsigned int)*nbas);
        unsigned int i, j, nprim, nctr;
        unsigned int ptr = 0;
        int l;
        const double *c;
        double maxc, expnt;
        for (i = 0; i < nbas; i++) {
                opt->ptr_log_coeff[i] = ptr;
                ptr += bas(NPRIM_OF,i);
        }

        opt->log_coeff = (double *)malloc(sizeof(double) * ptr);
        for (i = 0; i < nbas; i++) {
                nprim = bas(NPRIM_OF,i);
                nctr = bas(NCTR_OF,i);
                c = env + bas(PTR_COEFF,i);
                l = bas(ANG_OF,i);
                for (j = 0; j < nprim; j++) {
                        maxc = max_pgto_coeff(c, nprim, nctr, j);
                        expnt = env[bas(PTR_EXP,i)+j];
                        opt->log_coeff[i] = -log(maxc/CINTgto_norm(l, expnt));
                }
        }
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

