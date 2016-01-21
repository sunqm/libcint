/*
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * Breit = Gaunt + gauge
 * Gaunt ~ - \sigma1\dot\sigma2/r12
 * gauge ~  1/2 \sigma1\dot\sigma2/r12 - 1/2 (\sigma1\dot r12) (\sigma2\dot r12)/r12^3
 * Breit ~ -1/2 \sigma1\dot\sigma2/r12 - 1/2 (\sigma1\dot r12) (\sigma2\dot r12)/r12^3
 */

#include <stdlib.h>
#include <complex.h>
#include "cint_bas.h"
#include "optimizer.h"

#define DECLARE(X)      FINT X(double complex *opijkl, FINT *shls, \
                              FINT *atm, FINT natm, \
                              FINT *bas, FINT nbas, double *env, CINTOpt *opt)

#define BREIT0(X) \
DECLARE(cint2e_##X); \
DECLARE(cint2e_gauge_r1_##X); \
DECLARE(cint2e_gauge_r2_##X); \
void cint2e_breit_##X##_optimizer(CINTOpt **opt, FINT *atm, FINT natm, \
                                  FINT *bas, FINT nbas, double *env) \
{ \
        *opt = NULL; \
} \
FINT cint2e_breit_##X(double complex *opijkl, FINT *shls, \
                          FINT *atm, FINT natm, \
                          FINT *bas, FINT nbas, double *env, CINTOpt *opt) \
{ \
        FINT has_value = cint2e_##X(opijkl, shls, atm, natm, bas, nbas, env, NULL); \
 \
        const FINT ip = CINTcgto_spinor(shls[0], bas); \
        const FINT jp = CINTcgto_spinor(shls[1], bas); \
        const FINT kp = CINTcgto_spinor(shls[2], bas); \
        const FINT lp = CINTcgto_spinor(shls[3], bas); \
        const FINT nop = ip * jp * kp * lp; \
        double complex *buf = malloc(sizeof(double complex) * nop); \
        FINT i; \
        has_value = (cint2e_gauge_r1_##X(buf, shls, atm, natm, bas, nbas, env, NULL) || \
                     has_value); \
        /* [-1/2 gaunt] - [1/2 xxx*\sigma\dot r1] */ \
        if (has_value) { \
                for (i = 0; i < nop; i++) { \
                        opijkl[i] = -opijkl[i] - buf[i]; \
                } \
        } \
        /* ... [- 1/2 xxx*\sigma\dot(-r2)] */ \
        has_value = (cint2e_gauge_r2_##X(buf, shls, atm, natm, bas, nbas, env, NULL) || \
                     has_value); \
        if (has_value) { \
                for (i = 0; i < nop; i++) { \
                        opijkl[i] = (opijkl[i] + buf[i]) * .5; \
                } \
        } \
        free(buf); \
        return has_value; \
}


BREIT0(ssp1ssp2);
BREIT0(ssp1sps2);
BREIT0(sps1ssp2);
BREIT0(sps1sps2);
