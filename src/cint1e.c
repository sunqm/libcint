/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO integrals
 */

#include <stdlib.h>
#include <math.h>
#include "cint_bas.h"
#include "optimizer.h"
#include "g1e.h"
#include "cint1e.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"
#include "rys_roots.h"


/*
 * 1e GTO integral basic loop for < i|j>, no 1/r
 */
FINT CINT1e_loop(double *gctr, CINTEnvVars *envs, double *cache)
{
        FINT *shls  = envs->shls;
        FINT *bas = envs->bas;
        double *env = envs->env;
        FINT i_sh = shls[0];
        FINT j_sh = shls[1];
        FINT i_l = envs->i_l;
        FINT j_l = envs->j_l;
        FINT i_ctr = envs->x_ctr[0];
        FINT j_ctr = envs->x_ctr[1];
        FINT i_prim = bas(NPRIM_OF, i_sh);
        FINT j_prim = bas(NPRIM_OF, j_sh);
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        FINT nf = envs->nf;
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        FINT ip, jp, n;
        FINT has_value = 0;
        FINT *idx = malloc(sizeof(FINT) * nf * 3);
        double aij, dij, eij, rrij;
        double *g, *gout, *gctri;
        MALLOC_INSTACK(g, double, envs->g_size * 3 * ((1<<envs->gbits)+1)); // +1 as buffer
        MALLOC_INSTACK(gout, double, nf * n_comp);
        MALLOC_INSTACK(gctri, double, nf * i_ctr * n_comp);
        CINTg1e_index_xyz(idx, envs);

        rrij = CINTsquare_dist(ri, rj);
        double fac = envs->common_factor * CINTcommon_fac_sp(i_l) * CINTcommon_fac_sp(j_l);

        for (jp = 0; jp < j_prim; jp++) {
                envs->aj = aj[jp];
                n = nf * i_ctr * n_comp;
                CINTdset0(n, gctri);
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        aij = ai[ip] + aj[jp];
                        eij = (ai[ip] * aj[jp] / aij) * rrij;
                        if (eij > EXPCUTOFF)
                                continue;
                        has_value = 1;

                        dij = exp(-eij) / (aij * sqrt(aij)) * fac;
                        CINTg_ovlp(g, ai[ip], aj[jp], dij, envs);

                        CINTdset0(nf * n_comp, gout);
                        (*envs->f_gout)(gout, g, idx, envs, 1);

                        n = nf * n_comp;
                        CINTprim_to_ctr(gctri, n, gout, 1, i_prim, i_ctr, ci+ip);
                }
                n = nf * i_ctr;
                CINTprim_to_ctr(gctr, n, gctri, n_comp, j_prim, j_ctr, cj+jp);
        }
        free(idx);
        return has_value;
}

/*
 * For given charge distribution, calculate temporary parameter tau.
 * The charge parameter zeta is defined as    rho(r) = Norm * exp(-zeta*r^2)
 */
double CINTnuc_mod(double aij, FINT nuc_id, FINT *atm, double *env)
{
        double zeta;
        if (nuc_id < 0) {
                zeta = env[PTR_RINV_ZETA];
        } else if (atm(NUC_MOD_OF, nuc_id) == GAUSSIAN_NUC) {
                zeta = env[atm(PTR_ZETA, nuc_id)];
        } else {
                zeta = 0;
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
FINT CINT1e_nuc_loop(double *gctr, CINTEnvVars *envs, double fac, FINT nuc_id, double *cache)
{
        FINT *shls  = envs->shls;
        FINT *atm = envs->atm;
        FINT *bas = envs->bas;
        double *env = envs->env;
        FINT i_sh = shls[0];
        FINT j_sh = shls[1];
        FINT i_l = envs->i_l;
        FINT j_l = envs->j_l;
        FINT i_ctr = envs->x_ctr[0];
        FINT j_ctr = envs->x_ctr[1];
        FINT i_prim = bas(NPRIM_OF, i_sh);
        FINT j_prim = bas(NPRIM_OF, j_sh);
        FINT nf = envs->nf;
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        FINT ip, jp, i, n;
        FINT has_value = 0;
        double tau;
        double *cr;
        double x, u[MXRYSROOTS], w[MXRYSROOTS];
        FINT *idx = malloc(sizeof(FINT) * nf * 3);
        double rij[3], aij, dij, eij, rrij, t2;
        double *g, *gout, *gctri;
        MALLOC_INSTACK(g, double, envs->g_size * 3 * ((1<<envs->gbits)+1)); // +1 as buffer
        MALLOC_INSTACK(gout, double, nf * n_comp);
        MALLOC_INSTACK(gctri, double, nf * i_ctr * n_comp);

        if (nuc_id < 0) {
                cr = &env[PTR_RINV_ORIG];
        } else {
                cr = &env[atm(PTR_COORD, nuc_id)];
        }

        CINTg1e_index_xyz(idx, envs);

        rrij = CINTsquare_dist(ri, rj);
        fac *= envs->common_factor * CINTcommon_fac_sp(i_l) * CINTcommon_fac_sp(j_l);

        for (jp = 0; jp < j_prim; jp++) {
                envs->aj = aj[jp];
                n = nf * i_ctr * n_comp;
                CINTdset0(n, gctri);
                for (ip = 0; ip < i_prim; ip++) {
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

                                (*envs->f_gout)(gout, g, idx, envs, 1);
                        }

                        n = nf * n_comp;
                        CINTprim_to_ctr(gctri, n, gout, 1, i_prim, i_ctr, ci+ip);
                }
                n = nf * i_ctr;
                CINTprim_to_ctr(gctr, n, gctri, n_comp, j_prim, j_ctr, cj+jp);
        }
        free(idx);
        return has_value;
}


FINT int1e_cache_size(CINTEnvVars *envs)
{
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
        FINT len0 = envs->nf*n_comp;
        FINT cache_size = MAX(leng+len0*2+nc*n_comp*3,
                             nc*n_comp + envs->nf*8*OF_CMPLX);
        return cache_size;
}

/*
 * 1e integrals <i|O|j> without 1/r
 */
FINT CINT1e_drv(double *out, FINT *dims, CINTEnvVars *envs,
               double *cache, void (*f_c2s)(), FINT int1e_type)
{
        if (out == NULL) {
                return int1e_cache_size(envs);
        }
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *stack = NULL;
        if (cache == NULL) {
                FINT cache_size = int1e_cache_size(envs);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, double, nc*n_comp);

        FINT nout;
        FINT n;
        FINT has_value = 0;
        FINT *atm = envs->atm;
        double charge_fac;

        CINTdset0(nc*n_comp, gctr);
        for (n = 0; n < nc*n_comp; n++) {
                gctr[n] = 0;
        }
        switch (int1e_type) {
        case INT1E_TYPE_OVLP:
                has_value = CINT1e_loop(gctr, envs, cache);
                break;
        case INT1E_TYPE_RINV:
                has_value = CINT1e_nuc_loop(gctr, envs, 1, -1, cache);
                break;
        default:
                for (n = 0; n < envs->natm; n++) {
                        if (atm(CHARGE_OF,n) != 0) {
                                charge_fac = -abs(atm(CHARGE_OF,n));
                                has_value = CINT1e_nuc_loop(gctr, envs, charge_fac, n, cache)
                                        || has_value;
                        }
                }
        }

        FINT counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        if (f_c2s == &c2s_sph_1e) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->j_l*2+1) * x_ctr[1];
        } else if (f_c2s == &c2s_cart_1e) {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfj * x_ctr[1];
        }
        nout = dims[0] * dims[1];
        for (n = 0; n < n_comp; n++) {
                (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}

FINT CINT1e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs,
                       double *cache, void (*f_c2s)(), FINT int1e_type)
{
        if (out == NULL) {
                return int1e_cache_size(envs);
        }
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1] * envs->ncomp_e1;
        double *stack = NULL;
        if (cache == NULL) {
                FINT cache_size = int1e_cache_size(envs);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, double, nc*envs->ncomp_tensor);

        FINT nout;
        FINT n;
        FINT has_value = 0;
        FINT *atm = envs->atm;
        double charge_fac;

        CINTdset0(nc*envs->ncomp_tensor, gctr);
        switch (int1e_type) {
        case INT1E_TYPE_OVLP:
                has_value = CINT1e_loop(gctr, envs, cache);
                break;
        case INT1E_TYPE_RINV:
                has_value = CINT1e_nuc_loop(gctr, envs, 1, -1, cache);
                break;
        default:
                for (n = 0; n < envs->natm; n++) {
                        if (atm(CHARGE_OF,n) != 0) {
                                charge_fac = -abs(atm(CHARGE_OF,n));
                                has_value = CINT1e_nuc_loop(gctr, envs, charge_fac, n, cache)
                                        || has_value;
                        }
                }
        }

        FINT counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        nout = dims[0] * dims[1];
        for (n = 0; n < envs->ncomp_tensor; n++) {
                (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}

void int1e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                     FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}

