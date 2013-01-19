/*
 * File: cint2e.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 */


int cint2e_loop(double *gctr, const int *ng, const double fac,
                void (*const f_gout)(),
                const int *shls, const int *atm, const int *bas, const double *env);

int cint2e_drv(double *opkijl, int *ng, const double fac,
               void (*const f_gout)(), void (*const f_e1_c2s)(),
               void (*const f_e2_c2s)(),
               const int *shls, const int *atm, const int natm,
               const int *bas, const int nbas, const double *env);
