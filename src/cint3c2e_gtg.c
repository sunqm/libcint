/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 */

#include <stdlib.h>
#include "cint_bas.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint2e.h"
#include "cart2sph.h"

CACHE_SIZE_T int3c2e_gtg_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                  FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_gtg_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_drv(out, dims, &envs, opt, cache, &c2s_sph_3c2e1, 0);
}
void int3c2e_gtg_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                         FINT *bas, FINT nbas, double *env)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_3c2e_gtg_optimizer(opt, ng, atm, natm, bas, nbas, env);
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

ALL_CINT(int3c2e_gtg)

