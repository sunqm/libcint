/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO integrals
 */

#include <stdlib.h>
#include <math.h>
#include "cint_bas.h"
#include "g1e.h"
#include "cint1e.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"


/*
 * 1e GTO integral basic loop for < i|j>, no 1/r
 */
FINT CINT1e_loop(double *gctr, CINTEnvVars *envs, double fac)
{
        const FINT *shls  = envs->shls;
        const FINT *bas = envs->bas;
        const double *env = envs->env;
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT i_l = envs->i_l;
        const FINT j_l = envs->j_l;
        const FINT i_ctr = envs->i_ctr;
        const FINT j_ctr = envs->j_ctr;
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        const FINT nf = envs->nf;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        const double *ci = env + bas(PTR_COEFF, i_sh);
        const double *cj = env + bas(PTR_COEFF, j_sh);
        FINT ip, jp, n;
        FINT has_value = 0;
        FINT *const idx = malloc(sizeof(FINT) * nf * 3);
        double aij, dij, eij, rrij;
        double *g = malloc(sizeof(double) * envs->g_size * 3
                           * ((1<<envs->gbits)+1)); // +1 as buffer
        double *gout = malloc(sizeof(double) * nf * n_comp);
        double *gctri = malloc(sizeof(double) * nf * i_ctr * n_comp);

        CINTg1e_index_xyz(idx, envs);

        rrij = CINTsquare_dist(ri, rj);
        fac *= SQRTPI * M_PI * CINTcommon_fac_sp(i_l) * CINTcommon_fac_sp(j_l);

        for (jp = 0; jp < envs->j_prim; jp++) {
                envs->aj = aj[jp];
                n = nf * i_ctr * n_comp;
                CINTdset0(n, gctri);
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->ai = ai[ip];
                        aij = ai[ip] + aj[jp];
                        eij = (ai[ip] * aj[jp] / aij) * rrij;
                        if (eij > EXPCUTOFF)
                                continue;
                        has_value = 1;

                        dij = exp(-eij) / (aij * sqrt(aij)) * fac;
                        CINTg_ovlp(g, ai[ip], aj[jp], dij, envs);

                        CINTdset0(nf * n_comp, gout);
                        (*envs->f_gout)(g, gout, idx, envs);

                        n = nf * n_comp;
                        CINTprim_to_ctr(gctri, n, gout, 1, envs->i_prim,
                                        i_ctr, ci+ip);
                }
                n = nf * i_ctr;
                CINTprim_to_ctr(gctr, n, gctri, n_comp, envs->j_prim,
                                j_ctr, cj+jp);
        }
        free(g);
        free(idx);
        free(gout);
        free(gctri);

        return has_value;
}

/*
 * For given charge distribution, calculate temporary parameter tau.
 * The charge parameter zeta is defined as    rho(r) = Norm * exp(-zeta*r^2)
 */
static double CINTnuc_mod(const double aij, const FINT nuc_id,
                          const FINT *atm, const double *env)
{
        double zeta;
        if (nuc_id < 0) {
                zeta = env[PTR_RINV_ZETA];
        } else if (atm(NUC_MOD_OF, nuc_id) == POINT_NUC) {
                return 1;
        } else {
                zeta = env[atm(PTR_ZETA, nuc_id)];
        }

        if (zeta > 0) {
                return sqrt(zeta / (aij + zeta));
        } else {
                return 1;
        }
}

/*
 * 1e GTO integral basic loop for < i|1/r|j>, no 1/r
 * if nuc_id >= 0: nuclear attraction, use nuclear model
 * if nuc_id <  0: 1/r potential, do not use nuclear model
 */
FINT CINT1e_nuc_loop(double *gctr, CINTEnvVars *envs, double fac, FINT nuc_id)
{
        const FINT *shls  = envs->shls;
        const FINT *atm = envs->atm;
        const FINT *bas = envs->bas;
        const double *env = envs->env;
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT i_l = envs->i_l;
        const FINT j_l = envs->j_l;
        const FINT i_ctr = envs->i_ctr;
        const FINT j_ctr = envs->j_ctr;
        const FINT nf = envs->nf;
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        const double *ci = env + bas(PTR_COEFF, i_sh);
        const double *cj = env + bas(PTR_COEFF, j_sh);
        FINT ip, jp, i, n;
        FINT has_value = 0;
        double tau;
        const double *cr;
        double x, u[MXRYSROOTS], w[MXRYSROOTS];
        FINT *const idx = malloc(sizeof(FINT) * nf * 3);
        double rij[3], aij, dij, eij, rrij, t2;
        double *g = malloc(sizeof(double) * envs->g_size * 3
                           * ((1<<envs->gbits)+1)); // +1 as buffer
        double *const gout = malloc(sizeof(double) * nf * n_comp);
        double *const gctri = malloc(sizeof(double) * nf * i_ctr * n_comp);

        if (nuc_id < 0) {
                cr = &env[PTR_RINV_ORIG];
        } else {
                cr = &env[atm(PTR_COORD, nuc_id)];
        }

        CINTg1e_index_xyz(idx, envs);

        rrij = CINTsquare_dist(ri, rj);
        fac *= 2 * M_PI * CINTcommon_fac_sp(i_l) * CINTcommon_fac_sp(j_l);

        for (jp = 0; jp < envs->j_prim; jp++) {
                envs->aj = aj[jp];
                n = nf * i_ctr * n_comp;
                CINTdset0(n, gctri);
                for (ip = 0; ip < envs->i_prim; ip++) {
                        envs->ai = ai[ip];
                        aij = ai[ip] + aj[jp];
                        eij = (ai[ip] * aj[jp] / aij) * rrij;
                        if (eij > EXPCUTOFF)
                                continue;
                        has_value = 1;

                        rij[0] = (ai[ip] * ri[0] + aj[jp] * rj[0]) / aij;
                        rij[1] = (ai[ip] * ri[1] + aj[jp] * rj[1]) / aij;
                        rij[2] = (ai[ip] * ri[2] + aj[jp] * rj[2]) / aij;
                        tau = CINTnuc_mod(aij, nuc_id, atm, env);
                        x = aij * CINTsquare_dist(rij, cr) * tau * tau;
                        CINTrys_roots(envs->nrys_roots, x, u, w);

                        dij = exp(-eij) / aij * fac;
                        CINTdset0(nf * n_comp, gout);
                        for (i = 0; i < envs->nrys_roots; i++) {
                                t2 = u[i] / (1 + u[i]) * tau * tau;
                                CINTg_nuc(g, aij, rij, cr, t2,
                                          dij * w[i] * tau, envs);

                                (*envs->f_gout)(g, gout, idx, envs);
                        }

                        n = nf * n_comp;
                        CINTprim_to_ctr(gctri, n, gout, 1, envs->i_prim,
                                        i_ctr, ci+ip);
                }
                n = nf * i_ctr;
                CINTprim_to_ctr(gctr, n, gctri, n_comp, envs->j_prim,
                                j_ctr, cj+jp);
        }
        free(g);
        free(idx);
        free(gout);
        free(gctri);

        return has_value;
}


/*
 * 1e integrals <i|O|j> without 1/r
 */
FINT CINT1e_drv(double *opij, CINTEnvVars *envs, double fac,
               void (*const f_c2s)())
{
        const FINT *shls  = envs->shls;
        const FINT *bas = envs->bas;
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT i_ctr = envs->i_ctr;
        const FINT j_ctr = envs->j_ctr;
        const FINT nfi = envs->nfi;
        const FINT nfj = envs->nfj;
        const FINT nc = nfi * nfj * i_ctr * j_ctr * envs->ncomp_e1;
        FINT ip, jp, nop;
        FINT n;
        FINT has_value;
        double *gctr = malloc(sizeof(double) * nc * envs->ncomp_tensor);
        double *pgctr = gctr;

        CINTdset0(nc*envs->ncomp_tensor, gctr);
        has_value = CINT1e_loop(gctr, envs, fac);

        if (f_c2s == c2s_sph_1e) {
                ip = CINTcgto_spheric(i_sh, bas);
                jp = CINTcgto_spheric(j_sh, bas);
                nop = ip * jp;
        } else if (f_c2s == c2s_cart_1e) {
                ip = CINTcgto_cart(i_sh, bas);
                jp = CINTcgto_cart(j_sh, bas);
                nop = ip * jp;
        } else {
                ip = CINTcgto_spinor(i_sh, bas);
                jp = CINTcgto_spinor(j_sh, bas);
                nop = ip * jp * OF_CMPLX;
        }

        if (!has_value) {
                CINTdset0(nop * envs->ncomp_tensor, opij);
        } else {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        (*f_c2s)(opij, pgctr, envs);
                        opij += nop;
                        pgctr += nc;
                }
        }
        free(gctr);
        return has_value;
}


/*
 * 1e integrals <i|O|j> with 1/r
 */
FINT CINT1e_rinv_drv(double *opij, CINTEnvVars *envs, double fac,
                    void (*const f_c2s)())
{
        const FINT *shls  = envs->shls;
        const FINT *bas = envs->bas;
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT i_ctr = envs->i_ctr;
        const FINT j_ctr = envs->j_ctr;
        const FINT nfi = envs->nfi;
        const FINT nfj = envs->nfj;
        const FINT nc = nfi * nfj * i_ctr * j_ctr * envs->ncomp_e1;
        FINT ip, jp, nop;
        FINT n;
        FINT has_value;
        double *gctr = malloc(sizeof(double) * nc * envs->ncomp_tensor);
        double *pgctr = gctr;

        CINTdset0(nc*envs->ncomp_tensor, gctr);
        has_value = CINT1e_nuc_loop(gctr, envs, fac, -1);

        if (f_c2s == c2s_sph_1e) {
                ip = CINTcgto_spheric(i_sh, bas);
                jp = CINTcgto_spheric(j_sh, bas);
                nop = ip * jp;
        } else if (f_c2s == c2s_cart_1e) {
                ip = CINTcgto_cart(i_sh, bas);
                jp = CINTcgto_cart(j_sh, bas);
                nop = ip * jp;
        } else {
                ip = CINTcgto_spinor(i_sh, bas);
                jp = CINTcgto_spinor(j_sh, bas);
                nop = ip * jp * OF_CMPLX;
        }

        if (!has_value) {
                CINTdset0(nop * envs->ncomp_tensor, opij);
        } else {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        (*f_c2s)(opij, pgctr, envs);
                        opij += nop;
                        pgctr += nc;
                }
        }
        free(gctr);
        return has_value;
}


/*
 * 1e integrals <i|O|j> with nuclear attraction
 * TODO: add the gaussian nuclear model
 */
FINT CINT1e_nuc_drv(double *opij, CINTEnvVars *envs, double fac,
                   void (*const f_c2s)())
{
        const FINT *shls  = envs->shls;
        const FINT *atm = envs->atm;
        const FINT *bas = envs->bas;
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT i_ctr = envs->i_ctr;
        const FINT j_ctr = envs->j_ctr;
        const FINT nfi = envs->nfi;
        const FINT nfj = envs->nfj;
        const FINT nc = nfi * nfj * i_ctr * j_ctr * envs->ncomp_e1;
        FINT has_value = 0, has_value0;
        FINT ip, jp, nop;
        FINT n;
        double *gctr = malloc(sizeof(double) * nc * envs->ncomp_tensor);
        double *pgctr = gctr;

        CINTdset0(nc * envs->ncomp_tensor, gctr);
        for (n = 0; n < envs->natm; n++) {
                has_value0 = CINT1e_nuc_loop(gctr, envs,
                                             -fabs(atm(CHARGE_OF,n))*fac, n);
                has_value = has_value || has_value0;
        }

        if (f_c2s == c2s_sph_1e) {
                ip = CINTcgto_spheric(i_sh, bas);
                jp = CINTcgto_spheric(j_sh, bas);
                nop = ip * jp;
        } else if (f_c2s == c2s_cart_1e) {
                ip = CINTcgto_cart(i_sh, bas);
                jp = CINTcgto_cart(j_sh, bas);
                nop = ip * jp;
        } else {
                ip = CINTcgto_spinor(i_sh, bas);
                jp = CINTcgto_spinor(j_sh, bas);
                nop = ip * jp * OF_CMPLX;
        }

        if (!has_value) {
                CINTdset0(nop * envs->ncomp_tensor, opij);
        } else {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        (*f_c2s)(opij, pgctr, envs);
                        opij += nop;
                        pgctr += nc;
                }
        }
        free(gctr);
        return has_value;
}


/* template
 *
void gout1e(double *g, double *gout, const FINT *idx, const CINTEnvVars *envs)
{
        FINT nf = envs->nf;
        FINT ix, iy, iz, n;

        for (n = 0; n < nf; n++, idx+=3) {
                ix = idx[0];
                iy = idx[1];
                iz = idx[2];
                gout[n] += g[ix] * g[iy] * g[iz];
        }
}

FINT cint1e_ovlp_sph(double *opij, const FINT *shls,
                    const FINT *atm, const FINT natm,
                    const FINT *bas, const FINT nbas, const double *env)
{
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT i_l = bas(ANG_OF, i_sh);
        const FINT j_l = bas(ANG_OF, j_sh);
        FINT ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &gout1e;

        return CINT1e_drv(opij, envs, 1, &c2s_sph_1e);
}

FINT cint1e_ovlp(double *opij, const FINT *shls,
                const FINT *atm, const FINT natm,
                const FINT *bas, const FINT nbas, const double *env)
{
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT i_l = bas(ANG_OF, i_sh);
        const FINT j_l = bas(ANG_OF, j_sh);
        FINT ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &gout1e;

        return CINT1e_drv(opij, envs, 1, &c2s_sf_1e);
}

FINT cint1e_nuc_sph(double *opij, const FINT *shls,
                   const FINT *atm, const FINT natm,
                   const FINT *bas, const FINT nbas, const double *env)
{
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT i_l = bas(ANG_OF, i_sh);
        const FINT j_l = bas(ANG_OF, j_sh);
        FINT ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVars envs;
        CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &gout1e;

        return CINT1e_nuc_drv(opij, envs, 1, c2s_sph_1e);
}

FINT cint1e_nuc(double *opij, const FINT *shls,
               const FINT *atm, const FINT natm,
               const FINT *bas, const FINT nbas, const double *env)
{
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT i_l = bas(ANG_OF, i_sh);
        const FINT j_l = bas(ANG_OF, j_sh);
        FINT ng[] = {0, 0, 0, 0, 0, 1, 0, 1};

        return CINT1e_nuc_drv(opij, envs, 1, c2s_sf_1e);
}

// cint1e to fortran interface

C2F_(cint1e_ovlp_sph)
C2F_(cint1e_nuc_sph)
C2F_(cint1e_ovlp)
C2F_(cint1e_nuc)

template */
