/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 3-center 1-electron integrals
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "g2e.h"
#include "g3c1e.h"
#include "optimizer.h"
#include "cint1e.h"
#include "cint2e.h"
#include "misc.h"
#include "cart2sph.h"
#include "rys_roots.h"
#include "c2f.h"

#define PRIM2CTR0(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, ngp, gp, ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, ngp, gp, ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } \
        } \
        *ctrsymb##empty = 0


FINT CINT3c1e_loop_nopt(double *gctr, CINTEnvVars *envs, double *cache)
{
        FINT *shls  = envs->shls;
        FINT *bas = envs->bas;
        double *env = envs->env;
        FINT i_sh = shls[0];
        FINT j_sh = shls[1];
        FINT k_sh = shls[2];
        FINT i_ctr = envs->x_ctr[0];
        FINT j_ctr = envs->x_ctr[1];
        FINT k_ctr = envs->x_ctr[2];
        FINT i_prim = bas(NPRIM_OF, i_sh);
        FINT j_prim = bas(NPRIM_OF, j_sh);
        FINT k_prim = bas(NPRIM_OF, k_sh);
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *rk = envs->rk;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k;
        FINT ip, jp, kp;
        FINT empty[4] = {1, 1, 1, 1};
        FINT *iempty = empty + 0;
        FINT *jempty = empty + 1;
        FINT *kempty = empty + 2;
        FINT *gempty = empty + 3;
        /* COMMON_ENVS_AND_DECLARE end */
        const FINT nc = i_ctr * j_ctr * k_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenk = envs->nf * nc * n_comp; // gctrk
        const FINT lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenk + lenj + leni + len0;
        double *g;
        MALLOC_INSTACK(g, double, len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk;

        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctrj = gctrk;
                jempty = kempty;
        } else {
                gctrj = g1;
                g1 += lenj;
        }
        if (j_ctr == 1) {
                gctri = gctrj;
                iempty = jempty;
        } else {
                gctri = g1;
                g1 += leni;
        }
        if (i_ctr == 1) {
                gout = gctri;
                gempty = iempty;
        } else {
                gout = g1;
        }

        double eijk, dijk, aijk;
        double aiajrr, aiakrr, ajakrr;
        double rirk[3];
        double rjrk[3];
        rirk[0] = ri[0] - rk[0];
        rirk[1] = ri[1] - rk[1];
        rirk[2] = ri[2] - rk[2];
        rjrk[0] = rj[0] - rk[0];
        rjrk[1] = rj[1] - rk[1];
        rjrk[2] = rj[2] - rk[2];
        const double rr_ij = SQUARE(envs->rirj);
        const double rr_ik = SQUARE(      rirk);
        const double rr_jk = SQUARE(      rjrk);
        envs->idx = malloc(sizeof(FINT) * envs->nf * 3);
        CINTg2e_index_xyz(envs->idx, envs);

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *jempty = 1;
                }

                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        ajakrr = aj[jp] * ak[kp] * rr_jk;
                        for (ip = 0; ip < i_prim; ip++) {
                                envs->ai = ai[ip];
                                aijk = ai[ip] + aj[jp] + ak[kp];
                                aiakrr = ai[ip] * ak[kp] * rr_ik;
                                aiajrr = ai[ip] * aj[jp] * rr_ij;
                                eijk = (aiajrr+aiakrr+ajakrr) / aijk;
                                if (eijk > EXPCUTOFF) {
                                        continue;
                                }

                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*exp(-eijk);
                                } else {
                                        fac1i = fac1j*exp(-eijk);
                                }
                                dijk = fac1i / (aijk * sqrt(aijk));
                                CINTg3c1e_ovlp(g, ai[ip], aj[jp], ak[kp], dijk, envs);
                                (*envs->f_gout)(gout, g, envs->idx, envs, *gempty);

                                PRIM2CTR0(i, gout, envs->nf*n_comp);
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR0(j, gctri, envs->nf*i_ctr*n_comp);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR0(k, gctrj,envs->nf*i_ctr*j_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        free(envs->idx);
        return !*kempty;
}

FINT CINT3c1e_nuc_loop_nopt(double *gctr, CINTEnvVars *envs,
                            double fac, int nuc_id, double *cache)
{
        FINT *shls  = envs->shls;
        FINT *atm = envs->atm;
        FINT *bas = envs->bas;
        double *env = envs->env;
        FINT i_sh = shls[0];
        FINT j_sh = shls[1];
        FINT k_sh = shls[2];
        FINT i_ctr = envs->x_ctr[0];
        FINT j_ctr = envs->x_ctr[1];
        FINT k_ctr = envs->x_ctr[2];
        FINT i_prim = bas(NPRIM_OF, i_sh);
        FINT j_prim = bas(NPRIM_OF, j_sh);
        FINT k_prim = bas(NPRIM_OF, k_sh);
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *rk = envs->rk;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k;
        FINT i, ip, jp, kp;
        FINT empty[4] = {1, 1, 1, 1};
        FINT *iempty = empty + 0;
        FINT *jempty = empty + 1;
        FINT *kempty = empty + 2;
        FINT *gempty = empty + 3;
        int rys_empty;
        double *cr;
        double t2, tau, x, u[MXRYSROOTS], w[MXRYSROOTS];
        /* COMMON_ENVS_AND_DECLARE end */
        const FINT nc = i_ctr * j_ctr * k_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenk = envs->nf * nc * n_comp; // gctrk
        const FINT lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenk + lenj + leni + len0;
        double *g;
        MALLOC_INSTACK(g, double, len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk;

        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctrj = gctrk;
                jempty = kempty;
        } else {
                gctrj = g1;
                g1 += lenj;
        }
        if (j_ctr == 1) {
                gctri = gctrj;
                iempty = jempty;
        } else {
                gctri = g1;
                g1 += leni;
        }
        if (i_ctr == 1) {
                gout = gctri;
                gempty = iempty;
        } else {
                gout = g1;
        }

        if (nuc_id < 0) {
                cr = &env[PTR_RINV_ORIG];
        } else {
                cr = &env[atm(PTR_COORD, nuc_id)];
        }

        double eijk, dijk, aijk;
        double aiajrr, aiakrr, ajakrr;
        double rirk[3];
        double rjrk[3];
        double rijk[3];
        rirk[0] = ri[0] - rk[0];
        rirk[1] = ri[1] - rk[1];
        rirk[2] = ri[2] - rk[2];
        rjrk[0] = rj[0] - rk[0];
        rjrk[1] = rj[1] - rk[1];
        rjrk[2] = rj[2] - rk[2];
        const double rr_ij = SQUARE(envs->rirj);
        const double rr_ik = SQUARE(      rirk);
        const double rr_jk = SQUARE(      rjrk);
        envs->idx = malloc(sizeof(FINT) * envs->nf * 3);
        CINTg2e_index_xyz(envs->idx, envs);
        fac *= envs->common_factor;

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                if (k_ctr == 1) {
                        fac1k = fac * ck[kp];
                } else {
                        fac1k = fac;
                        *jempty = 1;
                }

                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        ajakrr = aj[jp] * ak[kp] * rr_jk;
                        for (ip = 0; ip < i_prim; ip++) {
                                envs->ai = ai[ip];
                                aijk = ai[ip] + aj[jp] + ak[kp];
                                aiakrr = ai[ip] * ak[kp] * rr_ik;
                                aiajrr = ai[ip] * aj[jp] * rr_ij;
                                eijk = (aiajrr+aiakrr+ajakrr) / aijk;
                                if (eijk > EXPCUTOFF) {
                                        continue;
                                }

                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*exp(-eijk);
                                } else {
                                        fac1i = fac1j*exp(-eijk);
                                }
                                dijk = fac1i / aijk;

                                rijk[0] = (ai[ip] * ri[0] + aj[jp] * rj[0] + ak[kp] * rk[0]) / aijk;
                                rijk[1] = (ai[ip] * ri[1] + aj[jp] * rj[1] + ak[kp] * rk[1]) / aijk;
                                rijk[2] = (ai[ip] * ri[2] + aj[jp] * rj[2] + ak[kp] * rk[2]) / aijk;
                                tau = CINTnuc_mod(aijk, nuc_id, atm, env);
                                x = aijk * CINTsquare_dist(rijk, cr) * tau * tau;
                                CINTrys_roots(envs->nrys_roots, x, u, w);
                                rys_empty = *gempty;
                                for (i = 0; i < envs->nrys_roots; i++) {
                                        t2 = u[i] / (1 + u[i]) * tau * tau;
                                        CINTg3c1e_nuc(g, ai[ip], aj[jp], ak[kp], rijk, cr, t2,
                                                      dijk * w[i] * tau, envs);

                                        (*envs->f_gout)(gout, g, envs->idx, envs, rys_empty);
                                        rys_empty = 0;
                                }

                                PRIM2CTR0(i, gout, envs->nf*n_comp);
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR0(j, gctri, envs->nf*i_ctr*n_comp);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR0(k, gctrj,envs->nf*i_ctr*j_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        free(envs->idx);
        return !*kempty;
}


FINT CINT3c1e_cart_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                       double *cache, FINT int_type)
{
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        if (out == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = leng + len0 + nc*n_comp*4;
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = leng + len0 + nc*n_comp*4;
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        FINT n;
        FINT has_value = 0;

        if (int_type == INT1E_TYPE_OVLP) {
                MALLOC_INSTACK(gctr, double, nc*n_comp);
                has_value = CINT3c1e_loop_nopt(gctr, envs, cache);
        } else if (int_type == INT1E_TYPE_RINV) {
                MALLOC_INSTACK(gctr, double, nc*n_comp);
                has_value = CINT3c1e_nuc_loop_nopt(gctr, envs, 1, -1, cache);
        } else {
                FINT *atm = envs->atm;
                FINT i;
                FINT buf_filled;
                double fac;
                double *buf;
                MALLOC_INSTACK(gctr, double, nc*n_comp);
                MALLOC_INSTACK(buf, double, nc*n_comp);
                for (i = 0; i < nc*n_comp; i++) {
                        gctr[i] = 0;
                }
                for (n = 0; n < envs->natm; n++) {
                        if (atm(CHARGE_OF,n) != 0) {
                                fac = -abs(atm(CHARGE_OF,n));
                                buf_filled = CINT3c1e_nuc_loop_nopt(buf, envs, fac, n, cache);
                                if (buf_filled) {
                                        for (i = 0; i < nc*n_comp; i++) {
                                                gctr[i] += buf[i];
                                        }
                                }
                                has_value |= buf_filled;
                        }
                }
                cache = buf; // release memory for c2s function
        }

        FINT counts[4];
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfj * x_ctr[1];
        counts[2] = envs->nfk * x_ctr[2];
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1] * dims[2];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_cart_3c2e1(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}
FINT CINT3c1e_spheric_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                         double *cache, void (*f_e1_c2s)(), FINT int_type, FINT is_ssc)
{
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        if (out == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = MAX(leng+len0+nc*n_comp*4, nc*n_comp+envs->nf*3);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = MAX(leng+len0+nc*n_comp*4, nc*n_comp+envs->nf*3);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        FINT n;
        FINT has_value = 0;

        if (int_type == INT1E_TYPE_OVLP) {
                MALLOC_INSTACK(gctr, double, nc*n_comp);
                has_value = CINT3c1e_loop_nopt(gctr, envs, cache);
        } else if (int_type == INT1E_TYPE_RINV) {
                MALLOC_INSTACK(gctr, double, nc*n_comp);
                has_value = CINT3c1e_nuc_loop_nopt(gctr, envs, 1, -1, cache);
        } else {
                FINT *atm = envs->atm;
                FINT i;
                FINT buf_filled;
                double fac;
                double *buf;
                MALLOC_INSTACK(gctr, double, nc*n_comp);
                MALLOC_INSTACK(buf, double, nc*n_comp);
                for (i = 0; i < nc*n_comp; i++) {
                        gctr[i] = 0;
                }
                for (n = 0; n < envs->natm; n++) {
                        if (atm(CHARGE_OF,n) != 0) {
                                fac = -abs(atm(CHARGE_OF,n));
                                buf_filled = CINT3c1e_nuc_loop_nopt(buf, envs, fac, n, cache);
                                if (buf_filled) {
                                        for (i = 0; i < nc*n_comp; i++) {
                                                gctr[i] += buf[i];
                                        }
                                }
                                has_value |= buf_filled;
                        }
                }
                cache = buf; // release memory for c2s function
        }

        FINT counts[4];
        counts[0] = (envs->i_l*2+1) * x_ctr[0];
        counts[1] = (envs->j_l*2+1) * x_ctr[1];
        if (is_ssc) {
                counts[2] = envs->nfk * x_ctr[2];
        } else {
                counts[2] = (envs->k_l*2+1) * x_ctr[2];
        }
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1] * dims[2];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        (*f_e1_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}

// TODO: ssc type c2s transformation
FINT CINT3c1e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)(), FINT int_type, FINT is_ssc)
{
        fprintf(stderr, "CINT3c1e_spinor_drv not implemented");
        exit(1);
}

void CINTgout3c1e(double *gout, double *g, FINT *idx, CINTEnvVars *envs, FINT gout_empty)
{
        FINT ix, iy, iz, n;

        if (gout_empty) {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        gout[n] = g[ix] * g[iy] * g[iz];
                }
        } else {
                for (n = 0; n < envs->nf; n++, idx+=3) {
                        ix = idx[0];
                        iy = idx[1];
                        iz = idx[2];
                        gout[n] += g[ix] * g[iy] * g[iz];
                }
        }
}

FINT int3c1e_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                 FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_spheric_drv(out, dims, &envs, opt, cache, &c2s_sph_3c1e, 0, 0);
}
void int3c1e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}

FINT int3c1e_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                  FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_cart_drv(out, dims, &envs, opt, cache, 0);
}

FINT int3c1e_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                   FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1, 0, 0);
}

FINT int3c1e_rinv_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                      FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_spheric_drv(out, dims, &envs, opt, cache, &c2s_sph_3c1e, 1, 0);
}
void int3c1e_rinv_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}

FINT int3c1e_rinv_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_cart_drv(out, dims, &envs, opt, cache, 1);
}

FINT int3c1e_rinv_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                         FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c1e;
        return CINT3c1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1, 1, 0);
}


ALL_CINT(int3c1e)
ALL_CINT_FORTRAN_(int3c1e)

