/*
 * File: cint2e.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO integrals
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint2e.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"

#define SQUARE(r)       (r)[0]*(r)[0] + (r)[1]*(r)[1] + (r)[2]*(r)[2]

#define COMMON_ENVS_AND_DECLARE \
        const unsigned int *shls = envs->shls; \
        const int *bas = envs->bas; \
        const double *env = envs->env; \
        const unsigned int i_prim = envs->i_prim; \
        const unsigned int j_prim = envs->j_prim; \
        const unsigned int k_prim = envs->k_prim; \
        const unsigned int l_prim = envs->l_prim; \
        const unsigned int i_ctr  = envs->i_ctr; \
        const unsigned int j_ctr  = envs->j_ctr; \
        const unsigned int k_ctr  = envs->k_ctr; \
        const unsigned int l_ctr  = envs->l_ctr; \
        const unsigned int nf     = envs->nf; \
        const double *ri = envs->ri; \
        const double *rj = envs->rj; \
        const double *rk = envs->rk; \
        const double *rl = envs->rl; \
        const unsigned int i_sh = shls[0]; \
        const unsigned int j_sh = shls[1]; \
        const unsigned int k_sh = shls[2]; \
        const unsigned int l_sh = shls[3]; \
        const double *ai = env + bas(PTR_EXP, i_sh); \
        const double *aj = env + bas(PTR_EXP, j_sh); \
        const double *ak = env + bas(PTR_EXP, k_sh); \
        const double *al = env + bas(PTR_EXP, l_sh); \
        const double *ci = env + bas(PTR_COEFF, i_sh); \
        const double *cj = env + bas(PTR_COEFF, j_sh); \
        const double *ck = env + bas(PTR_COEFF, k_sh); \
        const double *cl = env + bas(PTR_COEFF, l_sh); \
        const unsigned int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor; \
        double fac1i, fac1j, fac1k, fac1l; \
        unsigned int ip, jp, kp, lp, n; \
        int empty[5] = {1, 1, 1, 1, 1}; \
        int *iempty = empty + 0; \
        int *jempty = empty + 1; \
        int *kempty = empty + 2; \
        int *lempty = empty + 3; \
        int *gempty = empty + 4;

#define USE_OPT \
        unsigned int io, jo, ko, lo, off; \
        double eij, ekl, expij, expkl, expijkl; \
        double dist_ij, dist_kl; \
        double *prij; \
        if (opt) { \
                io = opt->prim_offset[i_sh]; \
                jo = opt->prim_offset[j_sh]; \
                ko = opt->prim_offset[k_sh]; \
                lo = opt->prim_offset[l_sh]; \
                if (opt->index_xyz_array) { \
                        envs->idx = opt->index_xyz_array[envs->i_l*ANG_MAX*ANG_MAX*ANG_MAX \
                                +envs->j_l*ANG_MAX*ANG_MAX+envs->k_l*ANG_MAX+envs->l_l]; \
                } else { \
                        envs->idx = malloc(sizeof(unsigned int) * nf * 3); \
                        CINTg2e_index_xyz(envs->idx, envs); \
                } \
        } else { \
                dist_ij = SQUARE(envs->rirj); \
                dist_kl = SQUARE(envs->rkrl); \
                envs->idx = malloc(sizeof(unsigned int) * nf * 3); \
                CINTg2e_index_xyz(envs->idx, envs); \
        }

#define SET_RIJ(I,J)    \
        envs->a##I = a##I[I##p]; \
        envs->a##I##J = a##I[I##p] + a##J[J##p]; \
        if (opt) { \
                off = I##o + I##p; \
                if (opt->screenij[J##o+J##p][off]) { \
                        goto I##_contracted; } \
                exp##I##J = opt->expij[J##o+J##p][off]; \
                prij = opt->rij[J##o+J##p]; \
                envs->r##I##J[0] = prij[off*3+0]; \
                envs->r##I##J[1] = prij[off*3+1]; \
                envs->r##I##J[2] = prij[off*3+2]; \
        } else { \
                e##I##J = dist_##I##J * a##I[I##p] * a##J[J##p] / envs->a##I##J; \
                if (e##I##J > EXPCUTOFF) { \
                        goto I##_contracted; } \
                envs->r##I##J[0] = (a##I[I##p]*r##I[0]+a##J[J##p]*r##J[0]) / envs->a##I##J; \
                envs->r##I##J[1] = (a##I[I##p]*r##I[1]+a##J[J##p]*r##J[1]) / envs->a##I##J; \
                envs->r##I##J[2] = (a##I[I##p]*r##I[2]+a##J[J##p]*r##J[2]) / envs->a##I##J; \
        } \
        envs->r##I##J##rx[0] = envs->r##I##J[0] - envs->rx_in_r##I##J##rx[0]; \
        envs->r##I##J##rx[1] = envs->r##I##J[1] - envs->rx_in_r##I##J##rx[1]; \
        envs->r##I##J##rx[2] = envs->r##I##J[2] - envs->rx_in_r##I##J##rx[2];

#define PRIM2CTR(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, ngp, gp, ctrsymb##_prim, \
                                          ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                } else { \
                        if (opt) { \
                                off = ctrsymb##o + ctrsymb##p; \
                                prim_to_ctr_opt(gctr##ctrsymb, ngp, gp, \
                                                opt->non0coeff[off], \
                                                opt->non0idx[off], \
                                                opt->non0ctr[off]); \
                        } else { \
                                CINTprim_to_ctr_1(gctr##ctrsymb, ngp, gp, ctrsymb##_prim, \
                                                  ctrsymb##_ctr, c##ctrsymb+ctrsymb##p); \
                        } \
                } \
        } \
        *ctrsymb##empty = 0

#define COPY_AND_CLOSING(GCTRL, EMPTY) \
        if (n_comp > 1 && !(EMPTY)) { \
                const unsigned int INC1 = 1; \
                double *gctr1, *gctr2, *gctr3; \
                switch (n_comp) { \
                case 3: \
                        gctr1 = gctr  + nf*nc; \
                        gctr2 = gctr1 + nf*nc; \
                        for (n = 0, ip = 0; n < nf*nc; n++, ip+=3) { \
                                gctr [n] = GCTRL[ip+0]; \
                                gctr1[n] = GCTRL[ip+1]; \
                                gctr2[n] = GCTRL[ip+2]; \
                        } \
                        break; \
                default: \
                        for (kp = 0; kp < n_comp-3; kp+=4) { \
                                gctr1 = gctr  + nf*nc; \
                                gctr2 = gctr1 + nf*nc; \
                                gctr3 = gctr2 + nf*nc; \
                                for (n = 0, ip = kp; n < nf*nc; n++,ip+=n_comp) { \
                                        gctr [n] = GCTRL[ip+0]; \
                                        gctr1[n] = GCTRL[ip+1]; \
                                        gctr2[n] = GCTRL[ip+2]; \
                                        gctr3[n] = GCTRL[ip+3]; \
                                } \
                                gctr += nf*nc * 4; \
                        } \
                        for (; kp < n_comp; kp++) { \
                                n = nf*nc; \
                                dcopy_(&n, GCTRL+kp, &n_comp, gctr, &INC1); \
                                gctr += nf*nc; \
                        } \
                } \
        } \
        free(g); \
        if (!opt || !opt->index_xyz_array) { free(envs->idx); } \
        return !(EMPTY);

void CINTprim_to_ctr_0(double *gc, const unsigned int nf, const double *gp,
                       const unsigned int nprim, const unsigned int nctr,
                       const double *coeff)
{
        unsigned int n, i;
        double c0, c1, c2, c3;
        double *p0, *p1, *p2, *p3;
        double non0coeff[32];
        double *non0pgc[32];
        unsigned int ncoeff = 0;

        for (i = 0; i < nctr; i++) {
                if (coeff[nprim*i] != 0) {
                        non0coeff[ncoeff] = coeff[nprim*i];
                        non0pgc[ncoeff] = gc + nf * i;
                        ncoeff++;
                } else { // need to initialize the memory, since += is used later
                        memset(gc+nf*i, 0, sizeof(double)*nf);
                }
        }

        switch (ncoeff) {
                case 1:
                        c0 = non0coeff[0];
                        p0 = non0pgc[0];
                        for (n = 0; n < nf; n++) {
                                p0[n] = c0 * gp[n];
                        }
                        break;
                case 2:
                        c0 = non0coeff[0];
                        c1 = non0coeff[1];
                        p0 = non0pgc[0];
                        p1 = non0pgc[1];
                        for (n = 0; n < nf; n++) {
                                p0[n] = c0 * gp[n];
                                p1[n] = c1 * gp[n];
                        }
                        break;
                case 3:
                        c0 = non0coeff[0];
                        c1 = non0coeff[1];
                        c2 = non0coeff[2];
                        p0 = non0pgc[0];
                        p1 = non0pgc[1];
                        p2 = non0pgc[2];
                        for (n = 0; n < nf; n++) {
                                p0[n] = c0 * gp[n];
                                p1[n] = c1 * gp[n];
                                p2[n] = c2 * gp[n];
                        }
                        break;
                case 4:
                        c0 = non0coeff[0];
                        c1 = non0coeff[1];
                        c2 = non0coeff[2];
                        c3 = non0coeff[3];
                        p0 = non0pgc[0];
                        p1 = non0pgc[1];
                        p2 = non0pgc[2];
                        p3 = non0pgc[3];
                        for (n = 0; n < nf; n++) {
                                p0[n] = c0 * gp[n];
                                p1[n] = c1 * gp[n];
                                p2[n] = c2 * gp[n];
                                p3[n] = c3 * gp[n];
                        }
                        break;
                default:
                        for (i = 0; i < ncoeff; i++) {
                                c1 = non0coeff[i];
                                p1 = non0pgc[i];
                                for (n = 0; n < nf; n++) {
                                        p1[n] = c1 * gp[n];
                                }
                        }
        }
}

void CINTprim_to_ctr_1(double *gc, const unsigned int nf, const double *gp,
                       const unsigned int nprim, const unsigned int nctr,
                       const double *coeff)
{
        unsigned int n, i;
        double c0, c1, c2, c3;
        double *p0, *p1, *p2, *p3;
        double non0coeff[32];
        double *non0pgc[32];
        unsigned int ncoeff = 0;

        for (i = 0; i < nctr; i++) {
                if (coeff[nprim*i] != 0) {
                        non0coeff[ncoeff] = coeff[nprim*i];
                        non0pgc[ncoeff] = gc + nf * i;
                        ncoeff++;
                }
        }

        switch (ncoeff) {
                case 1:
                        c0 = non0coeff[0];
                        p0 = non0pgc[0];
                        for (n = 0; n < nf; n++) {
                                p0[n] += c0 * gp[n];
                        }
                        break;
                case 2:
                        c0 = non0coeff[0];
                        c1 = non0coeff[1];
                        p0 = non0pgc[0];
                        p1 = non0pgc[1];
                        for (n = 0; n < nf; n++) {
                                p0[n] += c0 * gp[n];
                                p1[n] += c1 * gp[n];
                        }
                        break;
                case 3:
                        c0 = non0coeff[0];
                        c1 = non0coeff[1];
                        c2 = non0coeff[2];
                        p0 = non0pgc[0];
                        p1 = non0pgc[1];
                        p2 = non0pgc[2];
                        for (n = 0; n < nf; n++) {
                                p0[n] += c0 * gp[n];
                                p1[n] += c1 * gp[n];
                                p2[n] += c2 * gp[n];
                        }
                        break;
                case 4:
                        c0 = non0coeff[0];
                        c1 = non0coeff[1];
                        c2 = non0coeff[2];
                        c3 = non0coeff[3];
                        p0 = non0pgc[0];
                        p1 = non0pgc[1];
                        p2 = non0pgc[2];
                        p3 = non0pgc[3];
                        for (n = 0; n < nf; n++) {
                                p0[n] += c0 * gp[n];
                                p1[n] += c1 * gp[n];
                                p2[n] += c2 * gp[n];
                                p3[n] += c3 * gp[n];
                        }
                        break;
                default:
                        for (i = 0; i < ncoeff; i++) {
                                c1 = non0coeff[i];
                                p1 = non0pgc[i];
                                for (n = 0; n < nf; n++) {
                                        p1[n] += c1 * gp[n];
                                }
                        }
        }
}

static void prim_to_ctr_opt(double *gc, const unsigned int nf, const double *gp,
                            double *non0coeff, unsigned int *non0idx,
                            unsigned int non0ctr)
{
        unsigned int n, i;
        double c0, c1, c2, c3;
        double *p0, *p1, *p2, *p3;

        switch (non0ctr) {
                case 1:
                        c0 = non0coeff[0];
                        p0 = gc + nf*non0idx[0];
                        for (n = 0; n < nf; n++) {
                                p0[n] += c0 * gp[n];
                        }
                        break;
                case 2:
                        c0 = non0coeff[0];
                        c1 = non0coeff[1];
                        p0 = gc + nf*non0idx[0];
                        p1 = gc + nf*non0idx[1];
                        for (n = 0; n < nf; n++) {
                                p0[n] += c0 * gp[n];
                                p1[n] += c1 * gp[n];
                        }
                        break;
                case 3:
                        c0 = non0coeff[0];
                        c1 = non0coeff[1];
                        c2 = non0coeff[2];
                        p0 = gc + nf*non0idx[0];
                        p1 = gc + nf*non0idx[1];
                        p2 = gc + nf*non0idx[2];
                        for (n = 0; n < nf; n++) {
                                p0[n] += c0 * gp[n];
                                p1[n] += c1 * gp[n];
                                p2[n] += c2 * gp[n];
                        }
                        break;
                case 4:
                        c0 = non0coeff[0];
                        c1 = non0coeff[1];
                        c2 = non0coeff[2];
                        c3 = non0coeff[3];
                        p0 = gc + nf*non0idx[0];
                        p1 = gc + nf*non0idx[1];
                        p2 = gc + nf*non0idx[2];
                        p3 = gc + nf*non0idx[3];
                        for (n = 0; n < nf; n++) {
                                p0[n] += c0 * gp[n];
                                p1[n] += c1 * gp[n];
                                p2[n] += c2 * gp[n];
                                p3[n] += c3 * gp[n];
                        }
                        break;
                default:
                        for (i = 0; i < non0ctr; i++) {
                                c1 = non0coeff[i];
                                p1 = gc + nf*non0idx[i];
                                for (n = 0; n < nf; n++) {
                                        p1[n] += c1 * gp[n];
                                }
                        }
        }
}

inline void CINT2e_core(double *gout, double *g, double fac1i,
                        CINTEnvVars *envs, int empty)
{
        if (envs->g_size == 1) {
                if (empty) {
                        *gout = CINTg0_2e_ssss(fac1i, envs);
                } else { // if same address for gctri and gout and
                         // it has been initialized.
                        *gout+= CINTg0_2e_ssss(fac1i, envs);
                }
        } else {
                CINTg0_2e(g, fac1i, envs);
                (*envs->f_gout)(g, gout, envs->idx, envs, empty);
        }
}


// i_ctr = j_ctr = k_ctr = l_ctr = 1;
int CINT2e_1111_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;
        const unsigned int nc = 1;
        int len = envs->g_size * 3 * ((1<<envs->gbits)+1)
                + nf * n_comp;
        double *const g = (double *)malloc(sizeof(double) * len);
        double *gout;
        if (n_comp == 1) {
                gout = gctr;
        } else {
                gout = g + envs->g_size * 3 * ((1<<envs->gbits)+1);
        }

        USE_OPT;

        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];

                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++) {
                                        SET_RIJ(i, j);
                                        if (opt) {
                                                fac1i = fac1j*ci[ip]*expij*expkl;
                                        } else {
                                                fac1i = fac1j*ci[ip]*exp(-(eij+ekl));
                                        }
                                        CINT2e_core(gout, g, fac1i, envs, *empty);
                                        *empty = 0;
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        COPY_AND_CLOSING(gout, *empty);
}

// i_ctr = n; j_ctr = k_ctr = l_ctr = 1;
int CINT2e_n111_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const unsigned int nc = i_ctr;
        int len = envs->g_size * 3 * ((1<<envs->gbits)+1)
                + nf * i_ctr * n_comp // gctri
                + nf * n_comp; // gout
        double *const g = (double *)malloc(sizeof(double) * len);
        double *g1 = g + envs->g_size*3*((1<<envs->gbits)+1);
        double *gout, *gctri;
        if (n_comp == 1) {
                gctri = gctr;
        } else {
                gctri = g1;
                g1 += nf * nc * n_comp;
        }
        gout = g1;

        USE_OPT;

        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];

                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++) {
                                        SET_RIJ(i, j);
                                        if (opt) {
                                                fac1i = fac1j*expij*expkl;
                                        } else {
                                                fac1i = fac1j*exp(-(eij+ekl));
                                        }
                                        CINT2e_core(gout, g, fac1i, envs, 1);
                                        PRIM2CTR(i, gout, nf*n_comp);
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        COPY_AND_CLOSING(gctri, *iempty);
}

// j_ctr = n; i_ctr = k_ctr = l_ctr = 1;
int CINT2e_1n11_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const unsigned int nc = j_ctr;
        int len = envs->g_size * 3 * ((1<<envs->gbits)+1)
                + nf * j_ctr * n_comp // gctrj
                + nf * n_comp; // gout
        double *const g = (double *)malloc(sizeof(double) * len);
        double *g1 = g + envs->g_size*3*((1<<envs->gbits)+1);
        double *gout, *gctrj;
        if (n_comp == 1) {
                gctrj = gctr;
        } else {
                gctrj = g1;
                g1 += nf * nc * n_comp;
        }
        gout = g1;

        USE_OPT;

        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];

                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k;
                                *iempty = 1;
                                for (ip = 0; ip < i_prim; ip++) {
                                        SET_RIJ(i, j);
                                        if (opt) {
                                                fac1i = fac1j*ci[ip]*expij*expkl;
                                        } else {
                                                fac1i = fac1j*ci[ip]*exp(-(eij+ekl));
                                        }
                                        CINT2e_core(gout, g, fac1i, envs, *iempty);
                                        *iempty = 0;
i_contracted: ;
                                } // end loop i_prim
                                if (!*iempty) {
                                        PRIM2CTR(j, gout, nf*n_comp);
                                }
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        COPY_AND_CLOSING(gctrj, *jempty);
}

// k_ctr = n; i_ctr = j_ctr = l_ctr = 1;
int CINT2e_11n1_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const unsigned int nc = k_ctr;
        int len = envs->g_size * 3 * ((1<<envs->gbits)+1)
                + nf * k_ctr * n_comp // gctrk
                + nf * n_comp; // gout
        double *const g = (double *)malloc(sizeof(double) * len);
        double *g1 = g + envs->g_size*3*((1<<envs->gbits)+1);
        double *gout, *gctrk;
        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += nf * nc * n_comp;
        }
        gout = g1;

        USE_OPT;

        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor * cl[lp];
                for (kp = 0; kp < k_prim; kp++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l;
                        *jempty = 1;
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++) {
                                        SET_RIJ(i, j);
                                        if (opt) {
                                                fac1i = fac1j*ci[ip]*expij*expkl;
                                        } else {
                                                fac1i = fac1j*ci[ip]*exp(-(eij+ekl));
                                        }
                                        CINT2e_core(gout, g, fac1i, envs, *jempty);
                                        *jempty = 0;
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR(k, gout, nf*n_comp);
                        }
k_contracted: ;
                } // end loop k_prim
        } // end loop l_prim

        COPY_AND_CLOSING(gctrk, *kempty);
}

// l_ctr = n; i_ctr = j_ctr = k_ctr = 1;
int CINT2e_111n_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
{
        COMMON_ENVS_AND_DECLARE;

        const unsigned int nc = l_ctr;
        int len = envs->g_size * 3 * ((1<<envs->gbits)+1)
                + nf * l_ctr * n_comp // gctrl
                + nf * n_comp; // gout
        double *const g = (double *)malloc(sizeof(double) * len);
        double *g1 = g + envs->g_size*3*((1<<envs->gbits)+1);
        double *gout, *gctrl;
        if (n_comp == 1) {
                gctrl = gctr;
        } else {
                gctrl = g1;
                g1 += nf * nc * n_comp;
        }
        gout = g1;

        USE_OPT;

        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                fac1l = envs->common_factor;
                *kempty = 1;
                for (kp = 0; kp < k_prim; kp++) {
                        SET_RIJ(k, l);
                        fac1k = fac1l * ck[kp];
                        for (jp = 0; jp < j_prim; jp++) {
                                envs->aj = aj[jp];
                                fac1j = fac1k * cj[jp];
                                for (ip = 0; ip < i_prim; ip++) {
                                        SET_RIJ(i, j);
                                        if (opt) {
                                                fac1i = fac1j*ci[ip]*expij*expkl;
                                        } else {
                                                fac1i = fac1j*ci[ip]*exp(-(eij+ekl));
                                        }
                                        CINT2e_core(gout, g, fac1i, envs, *kempty);
                                        *kempty = 0;
i_contracted: ;
                                } // end loop i_prim
                        } // end loop j_prim
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
                        PRIM2CTR(l, gout, nf*n_comp);
                }
        } // end loop l_prim

        COPY_AND_CLOSING(gctrl, *lempty);
}


int CINT2e_loop(double *gctr, CINTEnvVars *envs, CINTOpt *opt)
{
        /* COMMON_ENVS_AND_DECLARE */
        const unsigned int *shls  = envs->shls;
        const int *bas = envs->bas;
        const double *env = envs->env;
        const unsigned int i_prim = envs->i_prim;
        const unsigned int j_prim = envs->j_prim;
        const unsigned int k_prim = envs->k_prim;
        const unsigned int l_prim = envs->l_prim;
        const unsigned int i_ctr  = envs->i_ctr;
        const unsigned int j_ctr  = envs->j_ctr;
        const unsigned int k_ctr  = envs->k_ctr;
        const unsigned int l_ctr  = envs->l_ctr;
        const unsigned int nf     = envs->nf;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *rk = envs->rk;
        const double *rl = envs->rl;
        const unsigned int i_sh = shls[0];
        const unsigned int j_sh = shls[1];
        const unsigned int k_sh = shls[2];
        const unsigned int l_sh = shls[3];
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        const double *ak = env + bas(PTR_EXP, k_sh);
        const double *al = env + bas(PTR_EXP, l_sh);
        const double *ci = env + bas(PTR_COEFF, i_sh);
        const double *cj = env + bas(PTR_COEFF, j_sh);
        const double *ck = env + bas(PTR_COEFF, k_sh);
        const double *cl = env + bas(PTR_COEFF, l_sh);
        const unsigned int n_comp = envs->ncomp_e1 * envs->ncomp_e2
                                  * envs->ncomp_tensor;
        double fac1i, fac1j, fac1k, fac1l;
        unsigned int ip, jp, kp, lp, n;
        int empty[5] = {1, 1, 1, 1, 1};
        int *iempty = empty + 0;
        int *jempty = empty + 1;
        int *kempty = empty + 2;
        int *lempty = empty + 3;
        int *gempty = empty + 4;
        /* COMMON_ENVS_AND_DECLARE end */
        const unsigned int nc = i_ctr * j_ctr * k_ctr * l_ctr;
        int len = envs->g_size * 3 * ((1<<envs->gbits)+1) // (irys,i,j,k,l,coord,0:1);
                + nf * nc * n_comp // gctrl
                + nf * i_ctr * j_ctr * k_ctr * n_comp // gctrk
                + nf * i_ctr * j_ctr * n_comp // gctrj
                + nf * i_ctr * n_comp // gctri
                + nf * n_comp; // gout
        double *const g = (double *)malloc(sizeof(double) * len);
        double *gout, *gctri, *gctrj, *gctrk, *gctrl;
        double *g1 = g + envs->g_size*3*((1<<envs->gbits)+1);

        if (n_comp == 1) {
                gctrl = gctr;
        } else {
                gctrl = g1;
                g1 += nf * nc * n_comp;
        }
        if (l_ctr == 1) {
                gctrk = gctrl;
                kempty = lempty;
        } else {
                gctrk = g1;
                g1 += nf*i_ctr*j_ctr*k_ctr*n_comp;
        }
        if (k_ctr == 1) {
                gctrj = gctrk;
                jempty = kempty;
        } else {
                gctrj = g1;
                g1 += nf*i_ctr*j_ctr*n_comp;
        }
        if (j_ctr == 1) {
                gctri = gctrj;
                iempty = jempty;
        } else {
                gctri = g1;
                g1 += nf*i_ctr*n_comp;
        }
        if (i_ctr == 1) {
                gout = gctri;
                gempty = iempty;
        } else {
                gout = g1;
        }

        /* USE_OPT */
        unsigned int io, jo, ko, lo, off;
        double eij, ekl, expij, expkl, expijkl;
        double dist_ij, dist_kl;
        double *prij;
        if (opt) {
                io = opt->prim_offset[i_sh];
                jo = opt->prim_offset[j_sh];
                ko = opt->prim_offset[k_sh];
                lo = opt->prim_offset[l_sh];
                if (opt->index_xyz_array) {
                        envs->idx = opt->index_xyz_array[envs->i_l*ANG_MAX*ANG_MAX*ANG_MAX
                                                        +envs->j_l*ANG_MAX*ANG_MAX
                                                        +envs->k_l*ANG_MAX
                                                        +envs->l_l];
                } else {
                        envs->idx = malloc(sizeof(unsigned int) * nf * 3);
                        CINTg2e_index_xyz(envs->idx, envs);
                }
        } else {
                dist_ij = SQUARE(envs->rirj);
                dist_kl = SQUARE(envs->rkrl);
                envs->idx = malloc(sizeof(unsigned int) * nf * 3);
                CINTg2e_index_xyz(envs->idx, envs);
        }
        /* USE_OPT end */

        *lempty = 1;
        for (lp = 0; lp < l_prim; lp++) {
                envs->al = al[lp];
                if (l_ctr == 1) {
                        fac1l = envs->common_factor * cl[lp];
                } else {
                        fac1l = envs->common_factor;
                        *kempty = 1;
                }
                for (kp = 0; kp < k_prim; kp++) {
                        /* SET_RIJ(k, l); */
                        envs->ak = ak[kp];
                        envs->akl = ak[kp] + al[lp];
                        if (opt) {
                                off = ko + kp;
                                if (opt->screenij[lo+lp][off]) {
                                        goto k_contracted;
                                }
                                expkl = opt->expij[lo+lp][off];
                                prij = opt->rij[lo+lp];
                                envs->rkl[0] = prij[off*3+0];
                                envs->rkl[1] = prij[off*3+1];
                                envs->rkl[2] = prij[off*3+2];
                        } else {
                                ekl = dist_kl * ak[kp] * al[lp] / envs->akl;
                                if (ekl > EXPCUTOFF) {
                                        goto k_contracted;
                                }
                                envs->rkl[0] = (ak[kp]*rk[0] + al[lp]*rl[0]) / envs->akl;
                                envs->rkl[1] = (ak[kp]*rk[1] + al[lp]*rl[1]) / envs->akl;
                                envs->rkl[2] = (ak[kp]*rk[2] + al[lp]*rl[2]) / envs->akl;
                        }
                        envs->rklrx[0] = envs->rkl[0] - envs->rx_in_rklrx[0];
                        envs->rklrx[1] = envs->rkl[1] - envs->rx_in_rklrx[1];
                        envs->rklrx[2] = envs->rkl[2] - envs->rx_in_rklrx[2];
                        /* SET_RIJ(k, l); end */
                        if (k_ctr == 1) {
                                fac1k = fac1l * ck[kp];
                        } else {
                                fac1k = fac1l;
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
                                        /* SET_RIJ(i, j); */
                                        envs->ai = ai[ip];
                                        envs->aij = ai[ip] + aj[jp];
                                        if (opt) {
                                                off = io + ip;
                                                if (opt->screenij[jo+jp][off]) {
                                                        goto i_contracted;
                                                }
                                                expij = opt->expij[jo+jp][off];
                                                prij = opt->rij[jo+jp];
                                                envs->rij[0] = prij[off*3+0];
                                                envs->rij[1] = prij[off*3+1];
                                                envs->rij[2] = prij[off*3+2];
                                        } else {
                                                eij = dist_ij * ai[ip] * aj[jp] / envs->aij;
                                                if (eij > EXPCUTOFF) {
                                                        goto i_contracted;
                                                }
                                                envs->rij[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) / envs->aij;
                                                envs->rij[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) / envs->aij;
                                                envs->rij[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) / envs->aij;
                                        }
                                        envs->rijrx[0] = envs->rij[0] - envs->rx_in_rijrx[0];
                                        envs->rijrx[1] = envs->rij[1] - envs->rx_in_rijrx[1];
                                        envs->rijrx[2] = envs->rij[2] - envs->rx_in_rijrx[2];
                                        /* SET_RIJ(i, j); end */
                                        if (opt) {
                                                expijkl = expij*expkl;
                                        } else {
                                                expijkl = exp(-(eij+ekl));
                                        }
                                        if (i_ctr == 1) {
                                                fac1i = fac1j*ci[ip]*expijkl;
                                        } else {
                                                fac1i = fac1j*expijkl;
                                        }
                                        CINT2e_core(gout, g, fac1i, envs, *gempty);
                                        PRIM2CTR(i, gout, nf*n_comp);
i_contracted: ;
                                } // end loop i_prim
                                if (!*iempty) {
                                        PRIM2CTR(j, gctri, nf*i_ctr*n_comp);
                                }
                        } // end loop j_prim
                        if (!*jempty) {
                                PRIM2CTR(k, gctrj, nf*i_ctr*j_ctr*n_comp);
                        }
k_contracted: ;
                } // end loop k_prim
                if (!*kempty) {
//TODO: merge this contraction with COPY_AND_CLOSING for n_comp>1
                        PRIM2CTR(l, gctrk, nf*i_ctr*j_ctr*k_ctr*n_comp);
                }
        } // end loop l_prim

        /* COPY_AND_CLOSING(gctrl, *lempty); */
        if (n_comp > 1 && !*lempty) {
                const unsigned int INC1 = 1;
                double *gctr1, *gctr2, *gctr3;
                switch (n_comp) {
                case 3:
                        gctr1 = gctr  + nf*nc;
                        gctr2 = gctr1 + nf*nc;
                        for (n = 0, ip = 0; n < nf*nc; n++, ip+=3) {
                                gctr [n] = gctrl[ip+0];
                                gctr1[n] = gctrl[ip+1];
                                gctr2[n] = gctrl[ip+2];
                        }
                        break;
                default:
                        for (kp = 0; kp < n_comp-3; kp+=4) {
                                gctr1 = gctr  + nf*nc;
                                gctr2 = gctr1 + nf*nc;
                                gctr3 = gctr2 + nf*nc;
                                for (n = 0, ip = kp; n < nf*nc; n++,ip+=n_comp) {
                                        gctr [n] = gctrl[ip+0];
                                        gctr1[n] = gctrl[ip+1];
                                        gctr2[n] = gctrl[ip+2];
                                        gctr3[n] = gctrl[ip+3];
                                }
                                gctr += nf*nc * 4;
                        }
                        for (; kp < n_comp; kp++) {
                                n = nf*nc;
                                dcopy_(&n, gctrl+kp, &n_comp, gctr, &INC1);
                                gctr += nf*nc;
                        }
                }
        }
        free(g);
        if (!opt || !opt->index_xyz_array) {
                free(envs->idx);
        }
        return !*lempty;
        /* COPY_AND_CLOSING(gctrl, *lempty); end */
}


static int (*CINTf_2e_loop[16])() = {
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_n111_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_loop,
        CINT2e_1n11_loop,
        CINT2e_loop,
        CINT2e_11n1_loop,
        CINT2e_111n_loop,
        CINT2e_1111_loop,
};

int CINT2e_cart_drv(double *opijkl, CINTEnvVars *envs, CINTOpt *opt)
{
        const unsigned int nc = envs->nf * envs->i_ctr * envs->j_ctr
                               * envs->k_ctr * envs->l_ctr * envs->ncomp_e1;
        double *const gctr = (double *)malloc(sizeof(double) * nc
                                              * envs->ncomp_e1 * envs->ncomp_tensor);
        double *pgctr = gctr;
        unsigned int ip, jp, kp, lp, nop;
        unsigned int n;
        int has_value;

        n = ((envs->i_ctr==1) << 3) + ((envs->j_ctr==1) << 2)
          + ((envs->k_ctr==1) << 1) +  (envs->l_ctr==1);
        has_value = CINTf_2e_loop[n](gctr, envs, opt);

        ip = CINTcgtos_cart(envs->shls[0], envs->bas);
        jp = CINTcgtos_cart(envs->shls[1], envs->bas);
        kp = CINTcgtos_cart(envs->shls[2], envs->bas);
        lp = CINTcgtos_cart(envs->shls[3], envs->bas);
        nop = ip * jp * kp * lp;

        if (has_value) {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        c2s_cart_2e1(opijkl, pgctr, envs->shls, envs->bas);
                        opijkl += nop;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nop * envs->ncomp_tensor, opijkl);
        }
        free(gctr);
        return has_value;
}
int CINT2e_spheric_drv(double *opijkl, CINTEnvVars *envs, CINTOpt *opt)
{
        const unsigned int nc = envs->nf * envs->i_ctr * envs->j_ctr
                               * envs->k_ctr * envs->l_ctr * envs->ncomp_e1;
        double *const gctr = (double *)malloc(sizeof(double) * nc
                                              * envs->ncomp_e2 * envs->ncomp_tensor);
        double *pgctr = gctr;
        unsigned int ip, jp, kp, lp, nop;
        unsigned int n;
        int has_value;

        n = ((envs->i_ctr==1) << 3) + ((envs->j_ctr==1) << 2)
          + ((envs->k_ctr==1) << 1) +  (envs->l_ctr==1);
        has_value = CINTf_2e_loop[n](gctr, envs, opt);

        ip = CINTcgtos_spheric(envs->shls[0], envs->bas);
        jp = CINTcgtos_spheric(envs->shls[1], envs->bas);
        kp = CINTcgtos_spheric(envs->shls[2], envs->bas);
        lp = CINTcgtos_spheric(envs->shls[3], envs->bas);
        nop = ip * jp * kp * lp;

        if (has_value) {
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        c2s_sph_2e1(opijkl, pgctr, envs->shls, envs->bas);
                        opijkl += nop;
                        pgctr += nc;
                }
        } else {
                CINTdset0(nop * envs->ncomp_tensor, opijkl);
        }
        free(gctr);
        return has_value;
}
int CINT2e_spinor_drv(double *opijkl, CINTEnvVars *envs, CINTOpt *opt,
                      void (*const f_e1_c2s)(), void (*const f_e2_c2s)())
{
        const unsigned int nc = envs->nf * envs->i_ctr * envs->j_ctr
                               * envs->k_ctr * envs->l_ctr * envs->ncomp_e1;
        double *const gctr = (double *)malloc(sizeof(double) * nc
                                              * envs->ncomp_e2 * envs->ncomp_tensor);
        double *pgctr = gctr;
        double *opij;
        unsigned int ip, jp, kp, lp, nop;
        unsigned int n, n1, n2;
        int has_value;

        n = ((envs->i_ctr==1) << 3) + ((envs->j_ctr==1) << 2)
          + ((envs->k_ctr==1) << 1) +  (envs->l_ctr==1);
        has_value = CINTf_2e_loop[n](gctr, envs, opt);

        ip = CINTcgtos_spinor(envs->shls[0], envs->bas);
        jp = CINTcgtos_spinor(envs->shls[1], envs->bas);
        kp = CINTcgtos_spinor(envs->shls[2], envs->bas);
        lp = CINTcgtos_spinor(envs->shls[3], envs->bas);
        nop = ip * jp * kp * lp;

        if (has_value) {
                n1 = ip * envs->nfk * envs->k_ctr
                        * envs->nfl * envs->l_ctr * jp * OF_CMPLX;
                opij = (double *)malloc(sizeof(double)*n1*envs->ncomp_e2);
                for (n = 0; n < envs->ncomp_tensor; n++) {
                        for (n2 = 0; n2 < envs->ncomp_e2; n2++) {
                                (*f_e1_c2s)(opij+n1*n2, pgctr, envs->shls,
                                            envs->bas);
                                pgctr += nc;
                        }
                        (*f_e2_c2s)(opijkl, opij, envs->shls, envs->bas);
                        opijkl += nop * OF_CMPLX;
                }
                free(opij);
        } else {
                CINTdset0(nop * OF_CMPLX * envs->ncomp_tensor, opijkl);
        }
        free(gctr);
        return has_value;
}


/*
 * <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
void CINTgout2e(double *g, double *gout, const unsigned int *idx,
                const CINTEnvVars *envs, int gout_empty)
{
        unsigned int nf = envs->nf;
        unsigned int i, ix, iy, iz, n;

        if (gout_empty) {
                switch (envs->nrys_roots) {
                        case 1:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix] * g[iy] * g[iz];
                                }
                                break;
                        case 2:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1];
                                }
                                break;
                        case 3:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2];
                                }
                                break;
                        case 4:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3];
                                }
                                break;
                        case 5:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4];
                                }
                                break;
                        case 6:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5];
                                }
                                break;
                        case 7:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6];
                                }
                                break;
                        default:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] = 0;
                                        for (i = 0; i < envs->nrys_roots; i++)
                                                gout[n] += g[ix+i] * g[iy+i] * g[iz+i];
                                }
                                break;
                } // end switch nroots
        } else { // not flag_acc
                switch (envs->nrys_roots) {
                        case 1:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] += g[ix] * g[iy] * g[iz];
                                }
                                break;
                        case 2:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1];
                                }
                                break;
                        case 3:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2];
                                }
                                break;
                        case 4:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3];
                                }
                                break;
                        case 5:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4];
                                }
                                break;
                        case 6:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5];
                                }
                                break;
                        case 7:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        gout[n] +=g[ix  ] * g[iy  ] * g[iz  ]
                                                + g[ix+1] * g[iy+1] * g[iz+1]
                                                + g[ix+2] * g[iy+2] * g[iz+2]
                                                + g[ix+3] * g[iy+3] * g[iz+3]
                                                + g[ix+4] * g[iy+4] * g[iz+4]
                                                + g[ix+5] * g[iy+5] * g[iz+5]
                                                + g[ix+6] * g[iy+6] * g[iz+6];
                                }
                                break;
                        default:
                                for (n = 0; n < nf; n++, idx+=3) {
                                        ix = idx[0];
                                        iy = idx[1];
                                        iz = idx[2];
                                        for (i = 0; i < envs->nrys_roots; i++)
                                                gout[n] += g[ix+i] * g[iy+i] * g[iz+i];
                                }
                                break;
                } // end switch nroots
        }
}


int cint2e_sph(double *opijkl, const unsigned int *shls,
               const int *atm, const int natm,
               const int *bas, const int nbas, const double *env,
               CINTOpt *opt)
{
        unsigned int ng[] = {0, 0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_spheric_drv(opijkl, &envs, opt);
}
void cint2e_sph_optimizer(CINTOpt **opt, const int *atm, const int natm,
                          const int *bas, const int nbas, const double *env)
{
        unsigned int ng[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        CINTuse_all_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

int cint2e_cart(double *opijkl, const unsigned int *shls,
                const int *atm, const int natm,
                const int *bas, const int nbas, const double *env,
                CINTOpt *opt)
{
        unsigned int ng[] = {0, 0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_cart_drv(opijkl, &envs, opt);
}
void cint2e_cart_optimizer(CINTOpt **opt, const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env)
{
        cint2e_sph_optimizer(opt, atm, natm, bas, nbas, env);
}


/*
 * spinor <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
int cint2e(double *opijkl, const unsigned int *shls,
           const int *atm, const int natm,
           const int *bas, const int nbas, const double *env,
           CINTOpt *opt)
{
        unsigned int ng[] = {0, 0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2e_spinor_drv(opijkl, &envs, opt, &c2s_sf_2e1, &c2s_sf_2e2);
}
void cint2e_optimizer(CINTOpt **opt, const int *atm, const int natm,
                      const int *bas, const int nbas, const double *env)
{
        cint2e_sph_optimizer(opt, atm, natm, bas, nbas, env);
}


/*
 * * * * * * * * * * * * * * * * * * * * *
 * c to fortran interface
 */

C2Fo_(cint2e_cart);
C2Fo_(cint2e_sph);
C2Fo_(cint2e);
OPTIMIZER2F_(cint2e_cart_optimizer);
OPTIMIZER2F_(cint2e_sph_optimizer);
OPTIMIZER2F_(cint2e_optimizer);

