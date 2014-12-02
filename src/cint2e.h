/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include "config.h"

void CINTprim_to_ctr_0(double *gc, const FINT nf, const double *gp,
                       const FINT nprim, const FINT nctr, const double *coeff);
void CINTprim_to_ctr_1(double *gc, const FINT nf, const double *gp,
                       const FINT nprim, const FINT nctr, const double *coeff);

void CINTgout2e(double *g, double *gout, const FINT *idx,
                const CINTEnvVars *envs, FINT gout_empty);

FINT CINT2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt);

FINT CINT2e_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt,
                void (*const f_e1_c2s)(), void (*const f_e2_c2s)());

FINT CINT2e_cart_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt);
FINT CINT2e_spheric_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt);
FINT CINT2e_spinor_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt,
                       void (*const f_e1_c2s)(), void (*const f_e2_c2s)());
