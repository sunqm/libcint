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

#define gctrg   gout
#define gctrm   gctr
#define mempty  empty
#define m_ctr   n_comp
#define ALIAS_ADDR_IF_EQUAL(x, y) \
        if (y##_ctr == 1) { \
                gctr##x = gctr##y; \
                x##empty = y##empty; \
        } else { \
                gctr##x = g1; \
                g1 += len##x; \
        }

#define PRIM2CTR0(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } \
        } \
        *ctrsymb##empty = 0

#define TRANSPOSE(a) \
        if (*empty) { \
                CINTdmat_transpose(gctr, a, nf*nc, n_comp); \
        } else { \
                CINTdplus_transpose(gctr, a, nf*nc, n_comp); \
        } \
        *empty = 0;


FINT CINT3c1e_loop_nopt(double *gctr, CINTEnvVars *envs, double *cache, FINT *empty)
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
        size_t nf = envs->nf;
        double fac1i, fac1j, fac1k;
        FINT ip, jp, kp;
        FINT _empty[4] = {1, 1, 1, 1};
        FINT *iempty = _empty + 0;
        FINT *jempty = _empty + 1;
        FINT *kempty = _empty + 2;
        FINT *gempty = _empty + 3;

        FINT *idx;
        MALLOC_INSTACK(idx, envs->nf * 3);
        CINTg2e_index_xyz(idx, envs);

        FINT *non0ctri, *non0ctrj, *non0ctrk;
        FINT *non0idxi, *non0idxj, *non0idxk;
        MALLOC_INSTACK(non0ctri, i_prim+j_prim+k_prim+i_prim*i_ctr+j_prim*j_ctr+k_prim*k_ctr);
        non0ctrj = non0ctri + i_prim;
        non0ctrk = non0ctrj + j_prim;
        non0idxi = non0ctrk + k_prim;
        non0idxj = non0idxi + i_prim*i_ctr;
        non0idxk = non0idxj + j_prim*j_ctr;
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);

        FINT nc = i_ctr * j_ctr * k_ctr;
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        size_t lenk = envs->nf * nc * n_comp; // gctrk
        size_t lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        size_t leni = envs->nf * i_ctr * n_comp; // gctri
        size_t len0 = envs->nf * n_comp; // gout
        size_t len = leng + lenk + lenj + leni + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk;

        ALIAS_ADDR_IF_EQUAL(k, m);
        ALIAS_ADDR_IF_EQUAL(j, k);
        ALIAS_ADDR_IF_EQUAL(i, j);
        ALIAS_ADDR_IF_EQUAL(g, i);

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
        double rr_ij = SQUARE(envs->rirj);
        double rr_ik = SQUARE(      rirk);
        double rr_jk = SQUARE(      rjrk);

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak[0] = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *jempty = 1;
                }

                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj[0] = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        ajakrr = aj[jp] * ak[kp] * rr_jk;
                        for (ip = 0; ip < i_prim; ip++) {
                                envs->ai[0] = ai[ip];
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
                                envs->fac[0] = dijk;
                                CINTg3c1e_ovlp(g, ai[ip], aj[jp], ak[kp], envs);
                                (*envs->f_gout)(gout, g, idx, envs, *gempty);

                                PRIM2CTR0(i, gout, len0);
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR0(j, gctri, leni);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR0(k, gctrj, lenj);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                TRANSPOSE(gctrk);
        }
        return !*empty;
}

FINT CINT3c1e_nuc_loop_nopt(double *gctr, CINTEnvVars *envs,
                            double fac, FINT nuc_id, double *cache, FINT *empty)
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
        size_t nf = envs->nf;
        double fac1i, fac1j, fac1k;
        FINT i, ip, jp, kp;
        FINT _empty[4] = {1, 1, 1, 1};
        FINT *iempty = _empty + 0;
        FINT *jempty = _empty + 1;
        FINT *kempty = _empty + 2;
        FINT *gempty = _empty + 3;
        FINT rys_empty;

        FINT *idx;
        MALLOC_INSTACK(idx, envs->nf * 3);
        CINTg2e_index_xyz(idx, envs);

        FINT *non0ctri, *non0ctrj, *non0ctrk;
        FINT *non0idxi, *non0idxj, *non0idxk;
        MALLOC_INSTACK(non0ctri, i_prim+j_prim+k_prim+i_prim*i_ctr+j_prim*j_ctr+k_prim*k_ctr);
        non0ctrj = non0ctri + i_prim;
        non0ctrk = non0ctrj + j_prim;
        non0idxi = non0ctrk + k_prim;
        non0idxj = non0idxi + i_prim*i_ctr;
        non0idxk = non0idxj + j_prim*j_ctr;
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
        CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);

        double *cr;
        double t2, tau, x, u[MXRYSROOTS], w[MXRYSROOTS];
        FINT nc = i_ctr * j_ctr * k_ctr;
        size_t leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        size_t lenk = envs->nf * nc * n_comp; // gctrk
        size_t lenj = envs->nf * i_ctr * j_ctr * n_comp; // gctrj
        size_t leni = envs->nf * i_ctr * n_comp; // gctri
        size_t len0 = envs->nf * n_comp; // gout
        size_t len = leng + lenk + lenj + leni + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrj, *gctrk;

        ALIAS_ADDR_IF_EQUAL(k, m);
        ALIAS_ADDR_IF_EQUAL(j, k);
        ALIAS_ADDR_IF_EQUAL(i, j);
        ALIAS_ADDR_IF_EQUAL(g, i);

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
        double rr_ij = SQUARE(envs->rirj);
        double rr_ik = SQUARE(      rirk);
        double rr_jk = SQUARE(      rjrk);
        fac *= envs->common_factor;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak[0] = ak[kp];
                if (k_ctr == 1) {
                        fac1k = fac * ck[kp];
                } else {
                        fac1k = fac;
                        *jempty = 1;
                }

                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj[0] = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        ajakrr = aj[jp] * ak[kp] * rr_jk;
                        for (ip = 0; ip < i_prim; ip++) {
                                envs->ai[0] = ai[ip];
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
                                        envs->fac[0] = dijk * w[i] * tau;
                                        CINTg3c1e_nuc(g, ai[ip], aj[jp], ak[kp], rijk, cr, t2, envs);
                                        (*envs->f_gout)(gout, g, idx, envs, rys_empty);
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
                TRANSPOSE(gctrk);
        }
        return !*empty;
}

#define PAIRDATA_NON0IDX_SIZE(ps) \
                FINT *bas = envs->bas; \
                FINT *shls  = envs->shls; \
                FINT i_prim = bas(NPRIM_OF, shls[0]); \
                FINT j_prim = bas(NPRIM_OF, shls[1]); \
                FINT k_prim = bas(NPRIM_OF, shls[2]); \
                FINT ps = (i_prim * x_ctr[0] \
                           + j_prim * x_ctr[1] \
                           + k_prim * x_ctr[2] \
                           + envs->nf*3);

CACHE_SIZE_T CINT3c1e_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                         double *cache, void (*f_e1_c2s)(), FINT int_type, FINT is_ssc)
{
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        if (out == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1);
                size_t len0 = envs->nf*n_comp;
                FINT cache_size = MAX(leng+len0+nc*n_comp*4 + pdata_size,
                                      nc*n_comp+envs->nf*3);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                PAIRDATA_NON0IDX_SIZE(pdata_size);
                size_t leng = envs->g_size*3*((1<<envs->gbits)+1);
                size_t len0 = envs->nf*n_comp;
                FINT cache_size = MAX(leng+len0+nc*n_comp*4 + pdata_size,
                                      nc*n_comp+envs->nf*3);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        FINT n;
        FINT empty = 1;

        if (int_type == INT1E_TYPE_OVLP) {
                MALLOC_INSTACK(gctr, nc*n_comp);
                CINT3c1e_loop_nopt(gctr, envs, cache, &empty);
        } else if (int_type == INT1E_TYPE_RINV) {
                MALLOC_INSTACK(gctr, nc*n_comp);
                CINT3c1e_nuc_loop_nopt(gctr, envs, 1, -1, cache, &empty);
        } else {
                FINT *atm = envs->atm;
                FINT i;
                double fac;
                double *buf;
                MALLOC_INSTACK(gctr, nc*n_comp);
                MALLOC_INSTACK(buf, nc*n_comp);
                for (n = 0; n < envs->natm; n++) {
                        if (atm(CHARGE_OF,n) != 0) {
                                fac = -abs(atm(CHARGE_OF,n));
                                CINT3c1e_nuc_loop_nopt(buf, envs, fac, n, cache, &empty);
                        }
                }
                cache = buf; // release memory for c2s function
        }

        FINT counts[4];
        if (f_e1_c2s == &c2s_sph_3c1e || f_e1_c2s == &c2s_sph_3c2e1) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->j_l*2+1) * x_ctr[1];
                if (is_ssc) {
                        counts[2] = envs->nfk * x_ctr[2];
                } else {
                        counts[2] = (envs->k_l*2+1) * x_ctr[2];
                }
        } else {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfj * x_ctr[1];
                counts[2] = envs->nfk * x_ctr[2];
        }
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1] * dims[2];
        if (!empty) {
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
        return !empty;
}

// TODO: ssc type c2s transformation
CACHE_SIZE_T CINT3c1e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)(), FINT int_type, FINT is_ssc)
{
        fprintf(stderr, "CINT3c1e_spinor_drv not implemented");
        exit(1);
}

void CINTgout1e(double *gout, double *g, FINT *idx, CINTEnvVars *envs, FINT empty);

CACHE_SIZE_T int3c1e_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                 FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT3c1e_drv(out, dims, &envs, opt, cache, &c2s_sph_3c1e, 0, 0);
}
void int3c1e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}

CACHE_SIZE_T int3c1e_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                  FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT3c1e_drv(out, dims, &envs, opt, cache, &c2s_cart_3c1e, 0, 0);
}

CACHE_SIZE_T int3c1e_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                   FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT3c1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1, 0, 0);
}

CACHE_SIZE_T int3c1e_rinv_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                      FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT3c1e_drv(out, dims, &envs, opt, cache, &c2s_sph_3c1e, 1, 0);
}
void int3c1e_rinv_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}

CACHE_SIZE_T int3c1e_rinv_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT3c1e_drv(out, dims, &envs, opt, cache, &c2s_cart_3c1e, 1, 0);
}

CACHE_SIZE_T int3c1e_rinv_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                         FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e;
        return CINT3c1e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1, 1, 0);
}


ALL_CINT(int3c1e)
ALL_CINT_FORTRAN_(int3c1e)

