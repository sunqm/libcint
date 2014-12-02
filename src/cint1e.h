/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

int CINT1e_loop(double *gctr, CINTEnvVars *envs, double fac);

int CINT1e_nuc_loop(double *gctr, CINTEnvVars *envs, double fac, int nuc_id);

int CINT1e_drv(double *opij, CINTEnvVars *envs, double fac,
               void (*const f_c2s)());

int CINT1e_rinv_drv(double *opij, CINTEnvVars *envs, double fac,
                    void (*const f_c2s)());

int CINT1e_nuc_drv(double *opij, CINTEnvVars *envs, double fac,
                    void (*const f_c2s)());

