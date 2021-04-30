/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#ifdef WITH_FORTRAN
#include "config.h"

#define ALL_CINT_FORTRAN_(NAME) \
FINT c##NAME##_sph_(double *out, FINT *shls, FINT *atm, FINT *natm, \
                    FINT *bas, FINT *nbas, double *env, size_t optptr_as_integer8) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        return NAME##_sph(out, NULL, shls, \
                          atm, *natm, bas, *nbas, env, *opt, NULL); \
} \
void c##NAME##_sph_optimizer_(size_t optptr_as_integer8, FINT *atm, FINT *natm, \
                              FINT *bas, FINT *nbas, double *env) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
} \
FINT c##NAME##_cart_(double *out, FINT *shls, FINT *atm, FINT *natm, \
                     FINT *bas, FINT *nbas, double *env, size_t optptr_as_integer8) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        return NAME##_cart(out, NULL, shls, \
                           atm, *natm, bas, *nbas, env, *opt, NULL); \
} \
void c##NAME##_cart_optimizer_(CINTOpt **opt, FINT *atm, FINT *natm, \
                               FINT *bas, FINT *nbas, double *env) { \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
} \
FINT c##NAME##_(double *out, FINT *shls, FINT *atm, FINT *natm, \
                FINT *bas, FINT *nbas, double *env, size_t optptr_as_integer8) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        return NAME##_spinor((double complex *)out, NULL, shls, \
                             atm, *natm, bas, *nbas, env, *opt, NULL); \
} \
void c##NAME##_optimizer_(size_t optptr_as_integer8, FINT *atm, FINT *natm, \
                         FINT *bas, FINT *nbas, double *env) { \
        CINTOpt **opt = (CINTOpt **)optptr_as_integer8; \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
}

#define ALL_CINT1E_FORTRAN_(NAME) \
FINT c##NAME##_sph_(double *out, FINT *shls, FINT *atm, FINT *natm, \
                    FINT *bas, FINT *nbas, double *env) { \
        return NAME##_sph(out, NULL, shls, atm, *natm, bas, *nbas, env, NULL, NULL); \
} \
FINT c##NAME##_cart_(double *out, FINT *shls, FINT *atm, FINT *natm, \
                     FINT *bas, FINT *nbas, double *env) { \
        return NAME##_cart(out, NULL, shls, \
                           atm, *natm, bas, *nbas, env, NULL, NULL); \
} \
FINT c##NAME##_(double *out, FINT *shls, FINT *atm, FINT *natm, \
                FINT *bas, FINT *nbas, double *env) { \
        return NAME##_spinor((double complex *)out, NULL, shls, \
                             atm, *natm, bas, *nbas, env, NULL, NULL); \
}

#else

#define ALL_CINT_FORTRAN_(NAME)
#define ALL_CINT1E_FORTRAN_(NAME)

#endif
