/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include "config.h"

FINT CINT1e_loop(double *gctr, CINTEnvVars *envs, double fac);

FINT CINT1e_nuc_loop(double *gctr, CINTEnvVars *envs, double fac, FINT nuc_id);

FINT CINT1e_drv(double *opij, CINTEnvVars *envs, double fac,
                void (*const f_c2s)());

FINT CINT1e_rinv_drv(double *opij, CINTEnvVars *envs, double fac,
                     void (*const f_c2s)());

FINT CINT1e_nuc_drv(double *opij, CINTEnvVars *envs, double fac,
                    void (*const f_c2s)());
