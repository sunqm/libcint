/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 3-center 2-electron integrals
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "g3c2e.h"
#include "optimizer.h"
#include "cint2e.h"
#include "misc.h"
#include "fblas.h"
#include "cart2sph.h"
#include "c2f.h"

#define SQUARE(r)       (r)[0]*(r)[0] + (r)[1]*(r)[1] + (r)[2]*(r)[2]

#define PRIM2CTR0(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, ngp, gp, \
                                          envs->ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, ngp, gp, \
                                          envs->ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } \
        } \
        *ctrsymb##empty = 0


FINT CINT3c2e_loop_nopt(double *gctr, CINTEnvVars *envs)
{
        const FINT *shls  = envs->shls;
        const FINT *bas = envs->bas;
        const double *env = envs->env;
        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT k_sh = shls[2];
        const FINT i_ctr  = envs->i_ctr;
        const FINT j_ctr  = envs->j_ctr;
        const FINT k_ctr  = envs->k_ctr;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        const double *ak = env + bas(PTR_EXP, k_sh);
        const double *ci = env + bas(PTR_COEFF, i_sh);
        const double *cj = env + bas(PTR_COEFF, j_sh);
        const double *ck = env + bas(PTR_COEFF, k_sh);
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k;
        FINT ip, jp, kp, n;
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
        double *const g = (double *)malloc(sizeof(double)*len);
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
        envs->idx = (FINT *)malloc(sizeof(FINT) * envs->nf * 3);
        CINTg3c2e_index_xyz(envs->idx, envs);

        *kempty = 1;
        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *jempty = 1;
                }

                for (jp = 0; jp < envs->j_prim; jp++) {
                        envs->aj = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        for (ip = 0; ip < envs->i_prim; ip++) {
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
                                CINT2e_core(gout, g, fac1i, envs, *gempty);
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
        free(g);
        free(envs->idx);
        return !*kempty;
}


#define COMMON_ENVS_AND_DECLARE \
        const FINT *shls = envs->shls; \
        const FINT *bas = envs->bas; \
        const double *env = envs->env; \
        const FINT i_ctr  = envs->i_ctr; \
        const FINT j_ctr  = envs->j_ctr; \
        const FINT k_ctr  = envs->k_ctr; \
        const FINT i_sh = shls[0]; \
        const FINT j_sh = shls[1]; \
        const FINT k_sh = shls[2]; \
        const double *ri = envs->ri; \
        const double *rj = envs->rj; \
        const double *ai = env + bas(PTR_EXP, i_sh); \
        const double *aj = env + bas(PTR_EXP, j_sh); \
        const double *ak = env + bas(PTR_EXP, k_sh); \
        const double *ci = env + bas(PTR_COEFF, i_sh); \
        const double *cj = env + bas(PTR_COEFF, j_sh); \
        const double *ck = env + bas(PTR_COEFF, k_sh); \
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor; \
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
        envs->idx = opt->index_xyz_array[envs->i_l*ANG_MAX*ANG_MAX \
                                        +envs->j_l*ANG_MAX+envs->k_l]

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
                        CINTprim_to_ctr_0(gctr##ctrsymb, ngp, gp, \
                                          envs->ctrsymb##_prim, \
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
FINT CINT3c2e_111_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const FINT nc = 1;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT len0 = envs->nf * n_comp;
        const FINT len = leng + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
        double *gout;
        if (n_comp == 1) {
                gout = gctr;
        } else {
                gout = g + leng;
        }

        USE_OPT;

        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];

                for (jp = 0; jp < envs->j_prim; jp++) {
                        envs->aj = aj[jp];
                        fac1j = fac1k * cj[jp];
                        for (ip = 0; ip < envs->i_prim; ip++) {
                                SET_RIJ;
                                fac1i = fac1j*ci[ip]*expij;
                                CINT2e_core(gout, g, fac1i, envs, *empty);
                                *empty = 0;
i_contracted: ;
                        } // end loop i_prim
                } // end loop j_prim
        } // end loop k_prim

        if (n_comp > 1 && !*empty) {
                CINTdmat_transpose(gctr, gout, envs->nf*nc, n_comp);
        }
        free(g);
        return !*empty;
}

// i_ctr = n; j_ctr = k_ctr = 1;
FINT CINT3c2e_n11_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const FINT nc = i_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + leni + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
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

        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];

                for (jp = 0; jp < envs->j_prim; jp++) {
                        envs->aj = aj[jp];
                        fac1j = fac1k * cj[jp];
                        for (ip = 0; ip < envs->i_prim; ip++) {
                                SET_RIJ;
                                fac1i = fac1j*expij;
                                CINT2e_core(gout, g, fac1i, envs, 1);
                                PRIM2CTR(i, gout,envs->nf*n_comp);
i_contracted: ;
                        } // end loop i_prim
                } // end loop j_prim
        } // end loop k_prim

        if (n_comp > 1 && !*iempty) {
                CINTdmat_transpose(gctr, gctri, envs->nf*nc, n_comp);
        }
        free(g);
        return !*iempty;
}

// j_ctr = n; i_ctr = k_ctr = 1;
FINT CINT3c2e_1n1_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const FINT nc = j_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenj = envs->nf * j_ctr * n_comp; // gctrj
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenj + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
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

        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];

                for (jp = 0; jp < envs->j_prim; jp++) {
                        envs->aj = aj[jp];
                        fac1j = fac1k;
                        *iempty = 1;
                        for (ip = 0; ip < envs->i_prim; ip++) {
                                SET_RIJ;
                                fac1i = fac1j*ci[ip]*expij;
                                CINT2e_core(gout, g, fac1i, envs, *iempty);
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
        free(g);
        return !*jempty;
}

// k_ctr = n; i_ctr = j_ctr = l_ctr = 1;
FINT CINT3c2e_11n_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const FINT nc = k_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenk = envs->nf * k_ctr * n_comp; // gctrk
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenk + len0;
        double *const g = (double *)malloc(sizeof(double)*len);
        double *g1 = g + leng;
        double *gout, *gctrk;
        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        gout = g1;

        USE_OPT;

        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                fac1k = envs->common_factor * ck[kp];
                *jempty = 1;
                for (jp = 0; jp < envs->j_prim; jp++) {
                        envs->aj = aj[jp];
                        fac1j = fac1k * cj[jp];
                        for (ip = 0; ip < envs->i_prim; ip++) {
                                SET_RIJ;
                                fac1i = fac1j*ci[ip]*expij;
                                CINT2e_core(gout, g, fac1i, envs, *jempty);
                                *jempty = 0;
i_contracted: ;
                        } // end loop i_prim
                } // end loop j_prim
                if (!*jempty) {
                        PRIM2CTR(k, gout,envs->nf*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        free(g);
        return !*kempty;
}


FINT CINT3c2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt)
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
        double *const g = (double *)malloc(sizeof(double)*len);
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
        envs->idx = opt->index_xyz_array[envs->i_l*ANG_MAX*ANG_MAX
                                        +envs->j_l*ANG_MAX+envs->k_l];
        /* USE_OPT end */

        *kempty = 1;
        for (kp = 0; kp < envs->k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *jempty = 1;
                }

                for (jp = 0; jp < envs->j_prim; jp++) {
                        envs->aj = aj[jp];
                        if (j_ctr == 1) {
                                fac1j = fac1k * cj[jp];
                        } else {
                                fac1j = fac1k;
                                *iempty = 1;
                        }
                        for (ip = 0; ip < envs->i_prim; ip++) {
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
                                CINT2e_core(gout, g, fac1i, envs, *gempty);
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
        free(g);
        return !*kempty;
}

static FINT (*CINTf_3c2e_loop[8])() = {
        CINT3c2e_loop,
        CINT3c2e_loop,
        CINT3c2e_loop,
        CINT3c2e_n11_loop,
        CINT3c2e_loop,
        CINT3c2e_1n1_loop,
        CINT3c2e_11n_loop,
        CINT3c2e_111_loop,
};

FINT CINT3c2e_cart_drv(double *opijk, CINTEnvVars *envs, const CINTOpt *opt)
{
        const FINT ip = CINTcgto_cart(envs->shls[0], envs->bas);
        const FINT jp = CINTcgto_cart(envs->shls[1], envs->bas);
        const FINT kp = CINTcgto_cart(envs->shls[2], envs->bas);
        const FINT nop = ip * jp * kp;
        const FINT nc = envs->nf * envs->i_ctr * envs->j_ctr * envs->k_ctr;
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *const gctr = malloc(sizeof(double) * nc * n_comp);
        double *pgctr = gctr;
        FINT n;
        FINT has_value;

        if (opt != NULL) {
                n = ((envs->i_ctr==1) << 2) + ((envs->j_ctr==1) << 1)
                  + (envs->k_ctr==1);
                has_value = CINTf_3c2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs);
        }

        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        c2s_cart_3c2e1(opijk, pgctr, envs);
                        opijk += nop;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nop * n_comp, opijk);
        }
        free(gctr);
        return has_value;
}
FINT CINT3c2e_spheric_drv(double *opijk, CINTEnvVars *envs, const CINTOpt *opt,
                         void (*const f_e1_c2s)(), FINT is_ssc)
{
        const FINT ip = CINTcgto_spheric(envs->shls[0], envs->bas);
        const FINT jp = CINTcgto_spheric(envs->shls[1], envs->bas);
        FINT kp;
        if (is_ssc) {
                kp = CINTcgto_cart(envs->shls[2], envs->bas);
        } else {
                kp = CINTcgto_spheric(envs->shls[2], envs->bas);
        }
        const FINT nop = ip * jp * kp;
        const FINT nc = envs->nf * envs->i_ctr * envs->j_ctr * envs->k_ctr;
        const FINT n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *const gctr = malloc(sizeof(double) * nc * n_comp);
        double *pgctr = gctr;
        FINT n;
        FINT has_value;

        if (opt != NULL) {
                n = ((envs->i_ctr==1) << 2) + ((envs->j_ctr==1) << 1)
                  + (envs->k_ctr==1);
                has_value = CINTf_3c2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs);
        }

        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        (*f_e1_c2s)(opijk, pgctr, envs);
                        opijk += nop;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nop * n_comp, opijk);
        }
        free(gctr);
        return has_value;
}
FINT CINT3c2e_spinor_drv(double *opijk, CINTEnvVars *envs, const CINTOpt *opt,
                        void (*const f_e1_c2s)(), FINT is_ssc)
{
        const FINT ip = CINTcgto_spinor(envs->shls[0], envs->bas);
        const FINT jp = CINTcgto_spinor(envs->shls[1], envs->bas);
        FINT kp;
        if (is_ssc) {
                kp = CINTcgto_cart(envs->shls[2], envs->bas);
        } else {
                kp = CINTcgto_spheric(envs->shls[2], envs->bas);
        }
        const FINT nop = ip * jp * kp;
        const FINT nc = envs->nf * envs->i_ctr * envs->j_ctr
                                * envs->k_ctr * envs->ncomp_e1;
        double *gctr = malloc(sizeof(double)*nc*envs->ncomp_tensor);
        double *pgctr = gctr;
        FINT n;
        FINT has_value;

        if (opt != NULL) {
                n = ((envs->i_ctr==1) << 2) + ((envs->j_ctr==1) << 1)
                  + (envs->k_ctr==1);
                has_value = CINTf_3c2e_loop[n](gctr, envs, opt);
        } else {
                has_value = CINT3c2e_loop_nopt(gctr, envs);
        }

        if (has_value) {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        (*f_e1_c2s)(opijk, pgctr, envs);
                        pgctr += nc;
                        opijk += nop * OF_CMPLX;
                }
        } else {
                CINTdset0(nop * OF_CMPLX * envs->ncomp_tensor, opijk);
        }
        free(gctr);
        return has_value;
}


FINT cint3c2e_sph(double *opijk, const FINT *shls,
                 const FINT *atm, const FINT natm,
                 const FINT *bas, const FINT nbas, const double *env,
                 const CINTOpt *opt)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spheric_drv(opijk, &envs, opt, &c2s_sph_3c2e1, 0);
}
void cint3c2e_sph_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                          const FINT *bas, const FINT nbas, const double *env)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_3cindex_xyz(*opt, ng, atm, natm, bas, nbas, env);
}

FINT cint3c2e_cart(double *opijk, const FINT *shls,
                  const FINT *atm, const FINT natm,
                  const FINT *bas, const FINT nbas, const double *env,
                  const CINTOpt *opt)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_cart_drv(opijk, &envs, opt);
}
void cint3c2e_cart_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                           const FINT *bas, const FINT nbas, const double *env)
{
        cint3c2e_sph_optimizer(opt, atm, natm, bas, nbas, env);
}


FINT cint3c2e_spinor(double *opijk, const FINT *shls,
                    const FINT *atm, const FINT natm,
                    const FINT *bas, const FINT nbas, const double *env,
                    const CINTOpt *opt)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spinor_drv(opijk, &envs, opt, &c2s_sf_3c2e1, 0);
}
void cint3c2e_spinor_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                               const FINT *bas, const FINT nbas, const double *env)
{
        cint3c2e_sph_optimizer(opt, atm, natm, bas, nbas, env);
}

FINT cint3c2e_sph_ssc(double *opijk, const FINT *shls,
                     const FINT *atm, const FINT natm,
                     const FINT *bas, const FINT nbas, const double *env,
                     const CINTOpt *opt)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spheric_drv(opijk, &envs, opt, &c2s_sph_3c2e1_ssc, 1);
}
void cint3c2e_sph_ssc_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                                const FINT *bas, const FINT nbas, const double *env)
{
        cint3c2e_sph_optimizer(opt, atm, natm, bas, nbas, env);
}
FINT cint3c2e_spinor_ssc(double *opijk, const FINT *shls,
                        const FINT *atm, const FINT natm,
                        const FINT *bas, const FINT nbas, const double *env,
                        const CINTOpt *opt)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT3c2e_spinor_drv(opijk, &envs, opt, &c2s_sf_3c2e1_ssc, 1);
}
void cint3c2e_spinor_ssc_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                                   const FINT *bas, const FINT nbas, const double *env)
{
        cint3c2e_sph_optimizer(opt, atm, natm, bas, nbas, env);
}

void CINTgout3c2e_cint3c2e_spsp1_spinor(double *g,
double *gout, const FINT *idx, const CINTEnvVars *envs, FINT gout_empty);
void cint3c2e_spsp1_spinor_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                                     const FINT *bas, const FINT nbas, const double *env);
FINT cint3c2e_spsp1_spinor_ssc(double *opijkl, const FINT *shls,
                              const FINT *atm, const FINT natm,
                              const FINT *bas, const FINT nbas, const double *env,
                              CINTOpt *opt)
{
        FINT ng[] = {1, 1, 0, 0, 2, 4, 1, 1};
        CINTEnvVars envs;
        CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout3c2e_cint3c2e_spsp1_spinor;
        return CINT3c2e_spinor_drv(opijkl, &envs, opt, &c2s_si_3c2e1_ssc, 1);
}
void cint3c2e_spsp1_spinor_ssc_optimizer(CINTOpt **opt, const FINT *atm, const FINT natm,
                                         const FINT *bas, const FINT nbas, const double *env) {
        cint3c2e_spsp1_spinor_optimizer(opt, atm, natm, bas, nbas, env);
}


/*
 * * * * * * * * * * * * * * * * * * * * *
 * c to fortran interface
 */

C2Fo_(cint3c2e_cart);
C2Fo_(cint3c2e_sph);
C2Fo_(cint3c2e_sph_ssc);
C2Fo_(cint3c2e_spinor);
C2Fo_(cint3c2e_spinor_ssc);
C2Fo_(cint3c2e_spsp1_spinor_ssc);
OPTIMIZER2F_(cint3c2e_cart_optimizer);
OPTIMIZER2F_(cint3c2e_sph_optimizer);
OPTIMIZER2F_(cint3c2e_spinor_optimizer);
OPTIMIZER2F_(cint3c2e_sph_ssc_optimizer);
OPTIMIZER2F_(cint3c2e_spinor_ssc_optimizer);
OPTIMIZER2F_(cint3c2e_spsp1_spinor_ssc_optimizer);
