/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include "config.h"

#define ALL_CINT_FORTRAN_(NAME) \
int c##NAME##_sph_(double *out, int *shls, int *atm, int *natm, \
                    int *bas, int *nbas, double *env, CINTOpt *opt) { \
        return NAME##_sph(out, NULL, shls, \
                          atm, *natm, bas, *nbas, env, opt, NULL); \
} \
void c##NAME##_sph_optimizer_(CINTOpt **opt, int *atm, int *natm, \
                              int *bas, int *nbas, double *env) { \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
} \
int c##NAME##_cart_(double *out, int *shls, int *atm, int *natm, \
                    int *bas, int *nbas, double *env, CINTOpt *opt) { \
        return NAME##_cart(out, NULL, shls, \
                           atm, *natm, bas, *nbas, env, opt, NULL); \
} \
void c##NAME##_cart_optimizer_(CINTOpt **opt, int *atm, int *natm, \
                               int *bas, int *nbas, double *env) { \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
} \
int c##NAME##_(double *out, int *shls, int *atm, int *natm, \
               int *bas, int *nbas, double *env, CINTOpt *opt) { \
        return NAME##_spinor((double complex *)out, NULL, shls, \
                             atm, *natm, bas, *nbas, env, opt, NULL); \
} \
void c##NAME##_optimizer_(CINTOpt **opt, int *atm, int *natm, \
                         int *bas, int *nbas, double *env) { \
        NAME##_optimizer(opt, atm, *natm, bas, *nbas, env); \
}

#define ALL_CINT1E_FORTRAN_(NAME) \
int c##NAME##_sph_(double *out, int *shls, int *atm, int *natm, \
                    int *bas, int *nbas, double *env) { \
        return NAME##_sph(out, NULL, shls, atm, *natm, bas, *nbas, env, NULL, NULL); \
} \
int c##NAME##_cart_(double *out, int *shls, int *atm, int *natm, \
                    int *bas, int *nbas, double *env) { \
        return NAME##_cart(out, NULL, shls, \
                           atm, *natm, bas, *nbas, env, NULL, NULL); \
} \
int c##NAME##_(double *out, int *shls, int *atm, int *natm, \
               int *bas, int *nbas, double *env) { \
        return NAME##_spinor((double complex *)out, NULL, shls, \
                             atm, *natm, bas, *nbas, env, NULL, NULL); \
}
