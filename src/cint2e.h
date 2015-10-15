/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include "config.h"

void CINT2e_core(double *gout, double *g, double fac1i,
                 CINTEnvVars *envs, FINT empty);

void CINTgout2e(double *g, double *gout, const FINT *idx,
                const CINTEnvVars *envs, FINT gout_empty);

FINT CINT2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt);

FINT CINT2e_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt,
                void (*const f_e1_c2s)(), void (*const f_e2_c2s)());

FINT CINT2e_cart_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt);
FINT CINT2e_spheric_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt);
FINT CINT2e_spinor_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt,
                       void (*const f_e1_c2s)(), void (*const f_e2_c2s)());

FINT CINT3c2e_cart_drv(double *opijk, CINTEnvVars *envs, const CINTOpt *opt);
FINT CINT3c2e_spheric_drv(double *opijk, CINTEnvVars *envs, const CINTOpt *opt,
                         void (*const f_e1_c2s)(), FINT is_ssc);
FINT CINT3c2e_spinor_drv(double *opijk, CINTEnvVars *envs, const CINTOpt *opt,
                        void (*const f_e1_c2s)(), FINT is_ssc);
FINT CINT2c2e_cart_drv(double *opij, CINTEnvVars *envs, const CINTOpt *opt);
FINT CINT2c2e_spheric_drv(double *opij, CINTEnvVars *envs, const CINTOpt *opt);

