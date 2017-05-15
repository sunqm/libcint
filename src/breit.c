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

#define DECLARE(X)      FINT X(double complex *out, FINT *dims, FINT *shls, \
                              FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, \
                              CINTOpt *opt, double *cache)

#define BREIT0(X, ncomp_tensor) \
DECLARE(int2e_##X##_spinor); \
DECLARE(int2e_gauge_r1_##X##_spinor); \
DECLARE(int2e_gauge_r2_##X##_spinor); \
void int2e_breit_##X##_optimizer(CINTOpt **opt, FINT *atm, FINT natm, \
                                 FINT *bas, FINT nbas, double *env) \
{ \
        *opt = NULL; \
} \
FINT int2e_breit_##X##_spinor(double complex *out, FINT *dims, FINT *shls, \
                             FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, \
                             CINTOpt *opt, double *cache) \
{ \
        _int2e_breit_drv(out, dims, shls, atm, natm, bas, nbas, env, opt, cache, \
                         1, &int2e_##X##_spinor, \
                         &int2e_gauge_r1_##X##_spinor, &int2e_gauge_r2_##X##_spinor); \
} \
FINT cint2e_breit_##X(double complex *out, FINT *shls, \
                      FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env, \
                      CINTOpt *opt) \
{ \
        int2e_breit_##X##_spinor(out, NULL, shls, atm, natm, bas, nbas, env, opt, NULL); \
}

static void _copy_to_out(double complex *out, double complex *in, FINT *dims, FINT *counts)
{
        if (out == in) {
                return;
        }
        FINT ni = dims[0];
        FINT nj = dims[1];
        FINT nk = dims[2];
        FINT nij = ni * nj;
        FINT nijk = nij * nk;
        FINT di = counts[0];
        FINT dj = counts[1];
        FINT dk = counts[2];
        FINT dl = counts[3];
        FINT dij = di * dj;
        FINT dijk = dij * dk;
        FINT i, j, k, l;
        double complex *pin, *pout;
        for (l = 0; l < dl; l++) {
                for (k = 0; k < dk; k++) {
                        pin  = in  + k * dij;
                        pout = out + k * nij;
                        for (j = 0; j < dj; j++) {
                        for (i = 0; i < di; i++) {
                                pout[j*ni+i] = pin[j*di+i];
                        } }
                }
                in  += dijk;
                out += nijk;
        }
}

static FINT _int2e_breit_drv(double complex *out, FINT *dims, FINT *shls,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
                            CINTOpt *opt, double *cache, int ncomp_tensor,
                            int (*f_gaunt)(), int (*f_gauge_r1)(), int (*f_gauge_r2)())
{
        FINT counts[4];
        counts[0] = CINTcgto_spinor(shls[0], bas);
        counts[1] = CINTcgto_spinor(shls[1], bas);
        counts[2] = CINTcgto_spinor(shls[2], bas);
        counts[3] = CINTcgto_spinor(shls[3], bas);
        FINT nop = counts[0] * counts[1] * counts[2] * counts[3] * ncomp_tensor;
        double complex *buf = malloc(sizeof(double complex) * nop*2);
        double complex *buf1;
        if (dims == NULL) {
                dims = counts;
                buf1 = out;
        } else {
                buf1 = buf + nop;
        }

        FINT has_value = (*f_gaunt)(buf1, NULL, shls, atm, natm, bas, nbas, env, NULL, cache);

        FINT i;
        has_value = ((*f_gauge_r1)(buf, NULL, shls, atm, natm, bas, nbas, env, NULL, cache) ||
                     has_value);
        /* [1/2 gaunt] - [1/2 xxx*\sigma1\dot r1] */
        if (has_value) {
                for (i = 0; i < nop; i++) {
                        buf1[i] = -buf1[i] - buf[i];
                }
        }
        /* ... [- 1/2 xxx*\sigma1\dot(-r2)] */
        has_value = ((*f_gauge_r2)(buf, NULL, shls, atm, natm, bas, nbas, env, NULL, cache) ||
                     has_value);
        if (has_value) {
                for (i = 0; i < nop; i++) {
                        buf1[i] = (buf1[i] + buf[i]) * .5;
                }
        }
        _copy_to_out(out, buf1, dims, counts);
        free(buf);
        return has_value;
}


BREIT0(ssp1ssp2, 1);
BREIT0(ssp1sps2, 1);
BREIT0(sps1ssp2, 1);
BREIT0(sps1sps2, 1);
