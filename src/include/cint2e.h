/*
 * File: cint2e.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 */
void prim_to_ctr_0(double *gc, const unsigned int nf, const double *gp,
                   const unsigned int nprim, const unsigned int nctr,
                   const double *coeff);
void prim_to_ctr_1(double *gc, const unsigned int nf, const double *gp,
                   const unsigned int nprim, const unsigned int nctr,
                   const double *coeff);

void gout2e(double *g, double *gout, const unsigned int *idx,
            const CintEnvVars *envs, int gout_empty);

int cint2e_loop(double *gctr, const double fac,
                void (*const f_gout)(), CintEnvVars *envs);

int cint2e_drv(double *opkijl, unsigned int *ng, const double fac,
               void (*const f_gout)(),
               void (*const f_e1_c2s)(), void (*const f_e2_c2s)(),
               const unsigned int *shls, const int *atm, const int natm,
               const int *bas, const int nbas, const double *env);
