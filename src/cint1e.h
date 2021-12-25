/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <complex.h>
#include "config.h"

FINT CINT1e_loop(double *gctr, CINTEnvVars *envs, double *cache, FINT int1e_type);

CACHE_SIZE_T CINT1e_drv(double *out, FINT *dims, CINTEnvVars *envs,
               double *cache, void (*f_c2s)(), FINT int1e_type);

CACHE_SIZE_T CINT1e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs,
                       double *cache, void (*f_c2s)(), FINT int1e_type);

double CINTnuc_mod(double aij, FINT nuc_id, FINT *atm, double *env);

CACHE_SIZE_T int1e_cache_size(CINTEnvVars *envs);

CACHE_SIZE_T CINT3c1e_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                         double *cache, void (*f_e1_c2s)(), FINT int_type, FINT is_ssc);
CACHE_SIZE_T CINT3c1e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)(), FINT int_type, FINT is_ssc);

#define INT1E_TYPE_OVLP 0
#define INT1E_TYPE_RINV 1
#define INT1E_TYPE_NUC  2

CACHE_SIZE_T CINT1e_grids_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs,
                             double *cache, void (*f_c2s)());
CACHE_SIZE_T CINT1e_grids_drv(double *out, FINT *dims, CINTEnvVars *envs,
                      double *cache, void (*f_c2s)());
