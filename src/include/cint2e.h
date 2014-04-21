/*
 * File: cint2e.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 */

void CINTprim_to_ctr_0(double *gc, const unsigned int nf, const double *gp,
                       const unsigned int nprim, const unsigned int nctr,
                       const double *coeff);
void CINTprim_to_ctr_1(double *gc, const unsigned int nf, const double *gp,
                       const unsigned int nprim, const unsigned int nctr,
                       const double *coeff);

void CINTgout2e(double *g, double *gout, const unsigned int *idx,
                const CINTEnvVars *envs, int gout_empty);

int CINT2e_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt);

int CINT2e_drv(double *opijkl, CINTEnvVars *envs, CINTOpt *opt,
               void (*const f_e1_c2s)(), void (*const f_e2_c2s)());

int CINT2e_cart_drv(double *opijkl, CINTEnvVars *envs, CINTOpt *opt);
int CINT2e_spheric_drv(double *opijkl, CINTEnvVars *envs, CINTOpt *opt);
int CINT2e_spinor_drv(double *opijkl, CINTEnvVars *envs, CINTOpt *opt,
                      void (*const f_e1_c2s)(), void (*const f_e2_c2s)());
