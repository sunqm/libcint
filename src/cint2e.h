/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

void CINTprim_to_ctr_0(double *gc, const int nf, const double *gp,
                       const int nprim, const int nctr, const double *coeff);
void CINTprim_to_ctr_1(double *gc, const int nf, const double *gp,
                       const int nprim, const int nctr, const double *coeff);

void CINTgout2e(double *g, double *gout, const int *idx,
                const CINTEnvVars *envs, int gout_empty);

int CINT2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt);

int CINT2e_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt,
               void (*const f_e1_c2s)(), void (*const f_e2_c2s)());

int CINT2e_cart_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt);
int CINT2e_spheric_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt);
int CINT2e_spinor_drv(double *opijkl, CINTEnvVars *envs, const CINTOpt *opt,
                      void (*const f_e1_c2s)(), void (*const f_e2_c2s)());
