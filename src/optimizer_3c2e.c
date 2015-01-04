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
#include "cint_const.h"
#include "cint_bas.h"
#include "g3c2e.h"
#include "g2c2e.h"
#include "optimizer.h"
#include "misc.h"

#define MAX(X,Y) (X)>(Y)?(X):(Y)

void CINTOpt_set_3cindex_xyz(CINTOpt *opt, FINT *ng,
                             const FINT *atm, const FINT natm,
                             const FINT *bas, const FINT nbas, const double *env)
{
        FINT i, j, k, ptr;
// index_xyz_array only needs memory ~ ANG_MAX**3
// allocate ANG_MAX**4 to make it work with CINTdel_optimizer
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
        FINT shls[3];
        for (i = 0; i <= max_l; i++) {
        for (j = 0; j <= max_l; j++) {
        for (k = 0; k <= max_l; k++) {
                shls[0] = i; shls[1] = j; shls[2] = k;
                CINTinit_int3c2e_EnvVars(&envs, ng, shls,
                                         atm, natm, fakebas, fakenbas, env);

                ptr = i*ANG_MAX*ANG_MAX + j*ANG_MAX + k;
                opt->index_xyz_array[ptr] =
                        (FINT *)malloc(sizeof(FINT)*envs.nf*3);
                CINTg3c2e_index_xyz(opt->index_xyz_array[ptr], &envs);
        } } }
}

void CINTOpt_set_2cindex_xyz(CINTOpt *opt, FINT *ng,
                             const FINT *atm, const FINT natm,
                             const FINT *bas, const FINT nbas, const double *env)
{
        FINT i, j, ptr;
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
        FINT shls[2];
        for (i = 0; i <= max_l; i++) {
        for (j = 0; j <= max_l; j++) {
                shls[0] = i; shls[1] = j;
                CINTinit_int2c2e_EnvVars(&envs, ng, shls,
                                         atm, natm, fakebas, fakenbas, env);

                ptr = i*ANG_MAX + j;
                opt->index_xyz_array[ptr] =
                        (FINT *)malloc(sizeof(FINT)*envs.nf*3);
                CINTg1e_index_xyz(opt->index_xyz_array[ptr], &envs);
        } }
}

