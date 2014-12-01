/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

int CINT1e_loop(double *gctr, const int *ng, double fac,
                void (*const f_gout)(), const int *shls,
                const int *atm, const int *bas, const double *env);

int CINT1e_nuc_loop(double *gctr, const int *ng, double fac,
                void (*const f_gout)(), const int nuc_id,
                const int *shls,
                const int *atm, const int *bas, const double *env);

int CINT1e_drv(double *opij, int *ng, double fac,
               void (*const f_gout)(), void (*const f_c2s)(),
               const int *shls, const int *atm, const int natm,
               const int *bas, const int nbas, const double *env);

int CINT1e_rinv_drv(double *opij, int *ng, double fac,
                    void (*const f_gout)(), void (*const f_c2s)(),
                    const int *shls, const int *atm, const int natm,
                    const int *bas, const int nbas, const double *env);

int CINT1e_nuc_drv(double *opij, int *ng, double fac,
                   void (*const f_gout)(), void (*const f_c2s)(),
                   const int *shls, const int *atm, const int natm,
                   const int *bas, const int nbas, const double *env);

