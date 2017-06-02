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
#include "cart2sph.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint2e.h"
#include "misc.h"
#include "c2f.h"

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
        /* Gaunt */ \
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
        /* [1/2 gaunt] - [1/2 xxx*\sigma1\dot r1] */ \
        if (has_value) { \
                for (i = 0; i < nop; i++) { \
                        opijkl[i] = -opijkl[i] - buf[i]; \
                } \
        } \
        /* ... [- 1/2 xxx*\sigma1\dot(-r2)] */ \
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


/* modfied from
 * '("cint2e_breit_r1p2_sph"  spheric  ( nabla \, r0 \| dot nabla-r12 \| \, nabla ))
 */
static void CINTgout2e_cint2e_breit_r1p2_sph(double *g, double *gout,
        const FINT *idx, const CINTEnvVars *envs, FINT gout_empty) {
        const double *env = envs->env;
        const FINT nf = envs->nf;
        const FINT i_l = envs->i_l;
        const FINT j_l = envs->j_l;
        const FINT k_l = envs->k_l;
        const FINT l_l = envs->l_l;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *rk = envs->rk;
        const double *rl = envs->rl;
        FINT ix, iy, iz, i, n;
        double *g0 = g;
        double *g1 = g0 + envs->g_size * 3;
        double *g2 = g1 + envs->g_size * 3;
        double *g3 = g2 + envs->g_size * 3;
        double *g4 = g3 + envs->g_size * 3;
        double *g5 = g4 + envs->g_size * 3;
        double *g6 = g5 + envs->g_size * 3;
        double *g7 = g6 + envs->g_size * 3;
        double *g8 = g7 + envs->g_size * 3;
        double *g9 = g8 + envs->g_size * 3;
        double *g10 = g9 + envs->g_size * 3;
        double *g11 = g10 + envs->g_size * 3;
        double *g12 = g11 + envs->g_size * 3;
        double *g13 = g12 + envs->g_size * 3;
        double *g14 = g13 + envs->g_size * 3;
        double *g15 = g14 + envs->g_size * 3;
        double *g16 = g15 + envs->g_size * 3;
        double s[9];
        G2E_D_L(g1, g0, i_l+1, j_l+2, k_l+0, l_l+0);
        G2E_R0J(g3, g1, i_l+1, j_l+0, k_l, l_l);
        G2E_D_J(g4, g0, i_l+1, j_l+1, k_l, l_l);
        G2E_D_I(g5, g0, i_l+1, j_l+1, k_l, l_l);
        for (ix = 0; ix < envs->g_size * 3; ix++) {g4[ix] += g5[ix];}
        G2E_D_J(g5, g1, i_l+1, j_l+1, k_l, l_l);
        G2E_D_I(g6, g1, i_l+1, j_l+1, k_l, l_l);
        for (ix = 0; ix < envs->g_size * 3; ix++) {g5[ix] += g6[ix];}
        G2E_R0J(g7, g5, i_l+1, j_l+0, k_l, l_l);
        G2E_D_I(g12, g4, i_l+0, j_l, k_l, l_l);
        G2E_D_I(g15, g7, i_l+0, j_l, k_l, l_l);
        for (n = 0; n < nf; n++, idx+=3) {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                CINTdset0(9, s);
                for (i = 0; i < envs->nrys_roots; i++) {
                        s[0] += g15[ix+i] * g0[iy+i] * g0[iz+i];
                        s[1] += g12[ix+i] * g3[iy+i] * g0[iz+i];
                        s[2] += g12[ix+i] * g0[iy+i] * g3[iz+i];
                        s[3] += g3[ix+i] * g12[iy+i] * g0[iz+i];
                        s[4] += g0[ix+i] * g15[iy+i] * g0[iz+i];
                        s[5] += g0[ix+i] * g12[iy+i] * g3[iz+i];
                        s[6] += g3[ix+i] * g0[iy+i] * g12[iz+i];
                        s[7] += g0[ix+i] * g3[iy+i] * g12[iz+i];
                        s[8] += g0[ix+i] * g0[iy+i] * g15[iz+i];
                }
                if (gout_empty) {
                        gout[n] = s[0] + s[3] + s[6] + s[1] + s[4] + s[7] + s[2] + s[5] + s[8];
                } else {
                        gout[n] += s[0] + s[3] + s[6] + s[1] + s[4] + s[7] + s[2] + s[5] + s[8];
                }
        }
}
void cint2e_breit_r1p2_sph_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                                     const FINT *bas, const FINT nbas, const double *env) {
        FINT ng[] = {1, 2, 0, 1, 4, 1, 1, 1};
        CINTuse_all_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
FINT cint2e_breit_r1p2_sph(double *opijkl, const FINT *shls,
                           const FINT *atm, const FINT natm,
                           const FINT *bas, const FINT nbas, const double *env, CINTOpt *opt) {
        FINT ng[] = {1, 2, 0, 1, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_cint2e_breit_r1p2_sph;
        return CINT2e_spheric_drv(opijkl, &envs, opt);
}
OPTIMIZER2F_(cint2e_breit_r1p2_sph_optimizer);

/* modfied from
 * '("cint2e_breit_r2p2_sph"  spheric  ( nabla \, r0 \| dot nabla-r12 \| \, nabla ))
 */
static void CINTgout2e_cint2e_breit_r2p2_sph(double *g, double *gout,
        const FINT *idx, const CINTEnvVars *envs, FINT gout_empty) {
        const double *env = envs->env;
        const FINT nf = envs->nf;
        const FINT i_l = envs->i_l;
        const FINT j_l = envs->j_l;
        const FINT k_l = envs->k_l;
        const FINT l_l = envs->l_l;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *rk = envs->rk;
        const double *rl = envs->rl;
        FINT ix, iy, iz, i, n;
        double *g0 = g;
        double *g1 = g0 + envs->g_size * 3;
        double *g2 = g1 + envs->g_size * 3;
        double *g3 = g2 + envs->g_size * 3;
        double *g4 = g3 + envs->g_size * 3;
        double *g5 = g4 + envs->g_size * 3;
        double *g6 = g5 + envs->g_size * 3;
        double *g7 = g6 + envs->g_size * 3;
        double *g8 = g7 + envs->g_size * 3;
        double *g9 = g8 + envs->g_size * 3;
        double *g10 = g9 + envs->g_size * 3;
        double *g11 = g10 + envs->g_size * 3;
        double *g12 = g11 + envs->g_size * 3;
        double *g13 = g12 + envs->g_size * 3;
        double *g14 = g13 + envs->g_size * 3;
        double *g15 = g14 + envs->g_size * 3;
        double *g16 = g15 + envs->g_size * 3;
        double s[9];
        G2E_R0L(g2, g0, i_l+1, j_l+1, k_l+0, l_l+1);
        G2E_D_L(g3, g2, i_l+1, j_l+1, k_l+0, l_l+0);
        G2E_D_J(g4, g0, i_l+1, j_l+0, k_l, l_l);
        G2E_D_I(g5, g0, i_l+1, j_l+0, k_l, l_l);
        for (ix = 0; ix < envs->g_size * 3; ix++) {g4[ix] += g5[ix];}
        G2E_D_J(g7, g3, i_l+1, j_l+0, k_l, l_l);
        G2E_D_I(g8, g3, i_l+1, j_l+0, k_l, l_l);
        for (ix = 0; ix < envs->g_size * 3; ix++) {g7[ix] += g8[ix];}
        G2E_D_I(g12, g4, i_l+0, j_l, k_l, l_l);
        G2E_D_I(g15, g7, i_l+0, j_l, k_l, l_l);
        for (n = 0; n < nf; n++, idx+=3) {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                CINTdset0(9, s);
                for (i = 0; i < envs->nrys_roots; i++) {
                        s[0] += g15[ix+i] * g0[iy+i] * g0[iz+i];
                        s[1] += g12[ix+i] * g3[iy+i] * g0[iz+i];
                        s[2] += g12[ix+i] * g0[iy+i] * g3[iz+i];
                        s[3] += g3[ix+i] * g12[iy+i] * g0[iz+i];
                        s[4] += g0[ix+i] * g15[iy+i] * g0[iz+i];
                        s[5] += g0[ix+i] * g12[iy+i] * g3[iz+i];
                        s[6] += g3[ix+i] * g0[iy+i] * g12[iz+i];
                        s[7] += g0[ix+i] * g3[iy+i] * g12[iz+i];
                        s[8] += g0[ix+i] * g0[iy+i] * g15[iz+i];
                }
                if (gout_empty) {
                        gout[n] = s[0] + s[3] + s[6] + s[1] + s[4] + s[7] + s[2] + s[5] + s[8];
                } else {
                        gout[n] += s[0] + s[3] + s[6] + s[1] + s[4] + s[7] + s[2] + s[5] + s[8];
                }
        }
}
void cint2e_breit_r2p2_sph_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                                     const FINT *bas, const FINT nbas, const double *env) {
        FINT ng[] = {1, 1, 0, 2, 4, 1, 1, 1};
        CINTuse_all_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
FINT cint2e_breit_r2p2_sph(double *opijkl, const FINT *shls,
                           const FINT *atm, const FINT natm,
                           const FINT *bas, const FINT nbas, const double *env, CINTOpt *opt) {
        FINT ng[] = {1, 1, 0, 2, 4, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_cint2e_breit_r2p2_sph;
        return CINT2e_spheric_drv(opijkl, &envs, opt);
}
OPTIMIZER2F_(cint2e_breit_r2p2_sph_optimizer);
