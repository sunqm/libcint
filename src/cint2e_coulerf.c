/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 */

#include <stdlib.h>
#include "cint_bas.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint2e.h"
#include "cart2sph.h"

FINT CINTinit_int2e_coulerf_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                           FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTgout2e(double *gout, double *g, FINT *idx, CINTEnvVars *envs, FINT gout_empty);

CACHE_SIZE_T int2e_coulerf_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                      FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_coulerf_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_drv(out, dims, &envs, opt, cache, &c2s_sph_2e1);
}
void int2e_coulerf_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                             FINT *bas, FINT nbas, double *env)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

#define ALL_CINT(NAME) \
FINT c##NAME##_sph(double *out, FINT *shls, FINT *atm, FINT natm, \
            FINT *bas, FINT nbas, double *env, CINTOpt *opt) { \
        return NAME##_sph(out, NULL, shls, atm, natm, bas, nbas, env, opt, NULL); \
} \
void c##NAME##_sph_optimizer(CINTOpt **opt, FINT *atm, FINT natm, \
                         FINT *bas, FINT nbas, double *env) { \
        NAME##_optimizer(opt, atm, natm, bas, nbas, env); \
}

ALL_CINT(int2e_coulerf)
