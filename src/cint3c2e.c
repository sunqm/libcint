/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 3-center 2-electron integrals
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "optimizer.h"
#include "g2e.h"
#include "cint2e.h"
#include "misc.h"
#include "cart2sph.h"
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


FINT CINT3c2e_loop_nopt(double *gctr, CINTEnvVars *envs, double *cache)
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
        // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
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

        double eij, expij;
        const double dist_ij = SQUARE(envs->rirj);
        envs->idx = malloc(sizeof(FINT) * envs->nf * 3);
        CINTg2e_index_xyz(envs->idx, envs);

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
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
                        for (ip = 0; ip < i_prim; ip++) {
                                envs->ai = ai[ip];
                                envs->aij = ai[ip] + aj[jp];
                                eij = dist_ij * ai[ip] * aj[jp] / envs->aij;
                                if (eij > EXPCUTOFF) {
                                        goto i_contracted;
                                }
                                expij = exp(-eij);
                                envs->rij[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / envs->aij;
                                envs->rij[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / envs->aij;
                                envs->rij[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / envs->aij;
                                envs->rijrx[0] = envs->rij[0] - envs->rx_in_rijrx[0];
                                envs->rijrx[1] = envs->rij[1] - envs->rx_in_rijrx[1];
                                envs->rijrx[2] = envs->rij[2] - envs->rx_in_rijrx[2];
                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*expij;
                                } else {
                                        fac1i = fac1j*expij;
                                }
                                (*envs->f_g0_2e)(g, fac1i, envs);
                                (*envs->f_gout)(gout, g, envs->idx, envs, *gempty);
                                PRIM2CTR0(i, gout, envs->nf*n_comp);
i_contracted: ;
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


#define COMMON_ENVS_AND_DECLARE \
        FINT *shls = envs->shls; \
        FINT *bas = envs->bas; \
        double *env = envs->env; \
        FINT i_sh = shls[0]; \
        FINT j_sh = shls[1]; \
        FINT k_sh = shls[2]; \
        FINT i_ctr = envs->x_ctr[0]; \
        FINT j_ctr = envs->x_ctr[1]; \
        FINT k_ctr = envs->x_ctr[2]; \
        FINT i_prim = bas(NPRIM_OF, i_sh); \
        FINT j_prim = bas(NPRIM_OF, j_sh); \
        FINT k_prim = bas(NPRIM_OF, k_sh); \
        double *ri = envs->ri; \
        double *rj = envs->rj; \
        double *ai = env + bas(PTR_EXP, i_sh); \
        double *aj = env + bas(PTR_EXP, j_sh); \
        double *ak = env + bas(PTR_EXP, k_sh); \
        double *ci = env + bas(PTR_COEFF, i_sh); \
        double *cj = env + bas(PTR_COEFF, j_sh); \
        double *ck = env + bas(PTR_COEFF, k_sh); \
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor; \
        double fac1i, fac1j, fac1k; \
        FINT ip, jp, kp; \
        FINT empty[4] = {1, 1, 1, 1}; \
        FINT *iempty = empty + 0; \
        FINT *jempty = empty + 1; \
        FINT *kempty = empty + 2; \
        FINT *gempty = empty + 3;

#define USE_OPT \
        FINT off; \
        const FINT io = opt->prim_offset[i_sh]; \
        const FINT jo = opt->prim_offset[j_sh]; \
        const FINT ko = opt->prim_offset[k_sh]; \
        double eij, expij; \
        const double dist_ij = SQUARE(envs->rirj); \
        envs->idx = opt->index_xyz_array[envs->i_l*LMAX1*LMAX1 \
                                        +envs->j_l*LMAX1+envs->k_l]

#define SET_RIJ    \
        envs->ai  = ai[ip]; \
        envs->aij = ai[ip] + aj[jp]; \
        eij = dist_ij * ai[ip] * aj[jp] / envs->aij; \
        if (eij > EXPCUTOFF) { \
                goto i_contracted; \
        } \
        expij = exp(-eij); \
        envs->rij[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / envs->aij; \
        envs->rij[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / envs->aij; \
        envs->rij[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / envs->aij; \
        envs->rijrx[0] = envs->rij[0] - envs->rx_in_rijrx[0]; \
        envs->rijrx[1] = envs->rij[1] - envs->rx_in_rijrx[1]; \
        envs->rijrx[2] = envs->rij[2] - envs->rx_in_rijrx[2]

#define PRIM2CTR(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, ngp, gp, ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } else { \
                        off = ctrsymb##o + ctrsymb##p; \
                        CINTprim_to_ctr_opt(gctr##ctrsymb, ngp, gp, \
                                            opt->non0coeff[off], \
                                            opt->non0idx[off], \
                                            opt->non0ctr[off]); \
                } \
        } \
        *ctrsymb##empty = 0


// i_ctr = j_ctr = k_ctr = 1;
FINT CINT3c2e_111_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt, double *cache)
{
        COMMON_ENVS_AND_DECLARE;
        const FINT nc = 1;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT len0 = envs->nf * n_comp;
        const FINT len = leng + len0;
        double *g;
        MALLOC_INSTACK(g, double, len);
        double *gout;
        if (n_comp == 1) {
                gout = gctr;
        } else {
                gout = g + leng;
        }

        USE_OPT;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];

                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj = aj[jp];
                        fac1j = fac1k * cj[jp];
                        for (ip = 0; ip < i_prim; ip++) {
                                SET_RIJ;
                                fac1i = fac1j*ci[ip]*expij;
                                (*envs->f_g0_2e)(g, fac1i, envs);
                                (*envs->f_gout)(gout, g, envs->idx, envs, *empty);
                                *empty = 0;
i_contracted: ;
                        } // end loop i_prim
                } // end loop j_prim
        } // end loop k_prim

        if (n_comp > 1 && !*empty) {
                CINTdmat_transpose(gctr, gout, envs->nf*nc, n_comp);
        }
        return !*empty;
}

// i_ctr = n; j_ctr = k_ctr = 1;
FINT CINT3c2e_n11_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt, double *cache)
{
        COMMON_ENVS_AND_DECLARE;

        const FINT nc = i_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + leni + len0;
        double *g;
        MALLOC_INSTACK(g, double, len);
        double *g1 = g + leng;
        double *gout, *gctri;
        if (n_comp == 1) {
                gctri = gctr;
        } else {
                gctri = g1;
                g1 += leni;
        }
        gout = g1;

        USE_OPT;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];

                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj = aj[jp];
                        fac1j = fac1k * cj[jp];
                        for (ip = 0; ip < i_prim; ip++) {
                                SET_RIJ;
                                fac1i = fac1j*expij;
                                (*envs->f_g0_2e)(g, fac1i, envs);
                                (*envs->f_gout)(gout, g, envs->idx, envs, 1);
                                PRIM2CTR(i, gout,envs->nf*n_comp);
i_contracted: ;
                        } // end loop i_prim
                } // end loop j_prim
        } // end loop k_prim

        if (n_comp > 1 && !*iempty) {
                CINTdmat_transpose(gctr, gctri, envs->nf*nc, n_comp);
        }
        return !*iempty;
}

// j_ctr = n; i_ctr = k_ctr = 1;
FINT CINT3c2e_1n1_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt, double *cache)
{
        COMMON_ENVS_AND_DECLARE;

        const FINT nc = j_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenj = envs->nf * j_ctr * n_comp; // gctrj
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenj + len0;
        double *g;
        MALLOC_INSTACK(g, double, len);
        double *g1 = g + leng;
        double *gout, *gctrj;
        if (n_comp == 1) {
                gctrj = gctr;
        } else {
                gctrj = g1;
                g1 += lenj;
        }
        gout = g1;

        USE_OPT;

        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];

                for (jp = 0; jp < j_prim; jp++) {
                        envs->aj = aj[jp];
                        fac1j = fac1k;
                        *iempty = 1;
                        for (ip = 0; ip < i_prim; ip++) {
                                SET_RIJ;
                                fac1i = fac1j*ci[ip]*expij;
                                (*envs->f_g0_2e)(g, fac1i, envs);
                                (*envs->f_gout)(gout, g, envs->idx, envs, *iempty);
                                *iempty = 0;
i_contracted: ;
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR(j, gout,envs->nf*n_comp);
                        }
                } // end loop j_prim
        } // end loop k_prim

        if (n_comp > 1 && !*jempty) {
                CINTdmat_transpose(gctr, gctrj, envs->nf*nc, n_comp);
        }
        return !*jempty;
}


FINT CINT3c2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt, double *cache)
{
        COMMON_ENVS_AND_DECLARE;
        const FINT nc = i_ctr * j_ctr * k_ctr;
        // (irys,i,j,k,coord,0:1); +1 for nabla-r12
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

        /* USE_OPT */
        FINT off;
        const FINT io = opt->prim_offset[i_sh];
        const FINT jo = opt->prim_offset[j_sh];
        const FINT ko = opt->prim_offset[k_sh];
        double eij, expij;
        const double dist_ij = SQUARE(envs->rirj);
        envs->idx = opt->index_xyz_array[envs->i_l*LMAX1*LMAX1
                                        +envs->j_l*LMAX1+envs->k_l];
        /* USE_OPT end */

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
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
                        for (ip = 0; ip < i_prim; ip++) {
                                /* SET_RIJ; */
                                envs->ai = ai[ip];
                                envs->aij = ai[ip] + aj[jp];
                                eij = dist_ij * ai[ip] * aj[jp] / envs->aij;
                                if (eij > EXPCUTOFF) {
                                        goto i_contracted;
                                }
                                expij = exp(-eij);
                                envs->rij[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / envs->aij;
                                envs->rij[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / envs->aij;
                                envs->rij[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / envs->aij;
                                envs->rijrx[0] = envs->rij[0] - envs->rx_in_rijrx[0];
                                envs->rijrx[1] = envs->rij[1] - envs->rx_in_rijrx[1];
                                envs->rijrx[2] = envs->rij[2] - envs->rx_in_rijrx[2];
                                /* SET_RIJ; end */
                                if (i_ctr == 1) {
                                        fac1i = fac1j*ci[ip]*expij;
                                } else {
                                        fac1i = fac1j*expij;
                                }
                                (*envs->f_g0_2e)(g, fac1i, envs);
                                (*envs->f_gout)(gout, g, envs->idx, envs, *gempty);
                                PRIM2CTR(i, gout, envs->nf*n_comp);
i_contracted: ;
                        } // end loop i_prim
                        if (!*iempty) {
                                PRIM2CTR(j, gctri, envs->nf*i_ctr*n_comp);
                        }
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR(k, gctrj, envs->nf*i_ctr*j_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}

static FINT (*CINTf_3c2e_loop[8])() = {
        CINT3c2e_loop,
        CINT3c2e_loop,
        CINT3c2e_loop,
        CINT3c2e_n11_loop,
        CINT3c2e_loop,
        CINT3c2e_1n1_loop,
        CINT3c2e_loop,
        CINT3c2e_111_loop,
};

FINT CINT3c2e_cart_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                      double *cache)
{
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        if (out == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = leng+len0+nc*n_comp*3;
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = leng+len0+nc*n_comp*3;
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, double, nc*n_comp);

        FINT n;
        FINT has_value;

        if (opt != NULL && opt->expij != NULL) {
                n = ((envs->x_ctr[0]==1) << 2) + ((envs->x_ctr[1]==1) << 1) + (envs->x_ctr[2]==1);
                has_value = CINTf_3c2e_loop[n](gctr, envs, opt, cache);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs, cache);
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
FINT CINT3c2e_spheric_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                         double *cache, void (*f_e1_c2s)(), FINT is_ssc)
{
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        if (out == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = MAX(leng+len0+nc*n_comp*3, nc*n_comp+envs->nf*3);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = MAX(leng+len0+nc*n_comp*3, nc*n_comp+envs->nf*3);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, double, nc*n_comp);

        FINT n;
        FINT has_value;

        if (opt != NULL && opt->expij != NULL) {
                n = ((envs->x_ctr[0]==1) << 2) + ((envs->x_ctr[1]==1) << 1) + (envs->x_ctr[2]==1);
                has_value = CINTf_3c2e_loop[n](gctr, envs, opt, cache);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs, cache);
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
FINT CINT3c2e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)(), FINT is_ssc)
{
        FINT *x_ctr = envs->x_ctr;
        FINT counts[4];
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        if (is_ssc) {
                counts[2] = envs->nfk * x_ctr[2];
        } else {
                counts[2] = (envs->k_l*2+1) * x_ctr[2];
        }
        counts[3] = 1;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = MAX(leng+len0+nc*n_comp*3,
                                     nc*n_comp + envs->nf*14*OF_CMPLX);
                return cache_size;
        }
        double *stack = NULL;
        if (cache == NULL) {
                FINT leng = envs->g_size*3*((1<<envs->gbits)+1);
                FINT len0 = envs->nf*n_comp;
                FINT cache_size = MAX(leng+len0+nc*n_comp*3,
                                     nc*n_comp + envs->nf*14*OF_CMPLX);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, double, nc*n_comp);

        FINT n;
        FINT has_value;

        if (opt != NULL && opt->expij != NULL) {
                n = ((envs->x_ctr[0]==1) << 2) + ((envs->x_ctr[1]==1) << 1) + (envs->x_ctr[2]==1);
                has_value = CINTf_3c2e_loop[n](gctr, envs, opt, cache);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs, cache);
        }

        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1] * dims[2];
        if (has_value) {
                for (n = 0; n < envs->ncomp_e2 * envs->ncomp_tensor; n++) {
                        (*f_e1_c2s)(out+nout*n, gctr, dims, envs, cache);
                        gctr += nc * envs->ncomp_e1;
                }
        } else {
                for (n = 0; n < envs->ncomp_e2 * envs->ncomp_tensor; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}


FINT int3c2e_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spheric_drv(out, dims, &envs, opt, cache, &c2s_sph_3c2e1, 0);
}
void int3c2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_3c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

FINT int3c2e_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                 FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_cart_drv(out, dims, &envs, opt, cache);
}

FINT int3c2e_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                   FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1, 0);
}

FINT int3c2e_sph_ssc(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                    FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spheric_drv(out, dims, &envs, opt, cache, &c2s_sph_3c2e1_ssc, 1);
}
FINT int3c2e_spinor_ssc(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_3c2e1_ssc, 1);
}

void CINTgout2e_int3c2e_spsp1(double *g,
double *gout, FINT *idx, CINTEnvVars *envs, FINT gout_empty);
FINT int3c2e_spsp1_spinor_ssc(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                             FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {1, 1, 0, 0, 2, 4, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e_int3c2e_spsp1;
        return CINT3c2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_si_3c2e1_ssc, 1);
}


ALL_CINT(int3c2e)
ALL_CINT_FORTRAN_(int3c2e)

