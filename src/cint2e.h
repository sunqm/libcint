/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <complex.h>
#include "g1e.h"
#include "config.h"

void CINTgout2e(double *g, double *gout, FINT *idx,
                CINTEnvVars *envs, FINT gout_empty);

FINT CINT2e_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt, double *cache);

FINT CINT2e_cart_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                    double *cache);
FINT CINT2e_spheric_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                       double *cache);
FINT CINT2e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                      double *cache, void (*f_e1_c2s)(), void (*f_e2_c2s)());

FINT CINT3c2e_cart_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                      double *cache);
FINT CINT3c2e_spheric_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                         double *cache, void (*f_e1_c2s)(), FINT is_ssc);
FINT CINT3c2e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)(), FINT is_ssc);
FINT CINT2c2e_cart_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                      double *cache);
FINT CINT2c2e_spheric_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                      double *cache);
