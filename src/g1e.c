/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <string.h>
#include <math.h>
#include <assert.h>
#include "cint_bas.h"
#include "misc.h"
#include "g1e.h"


void CINTinit_int1e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                           FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nf = envs->nfi * envs->nfj;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->common_factor = 1;

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->nrys_roots =(envs->li_ceil + envs->lj_ceil)/2 + 1;

        assert(i_sh < SHLS_MAX);
        assert(j_sh < SHLS_MAX);
        assert(envs->i_l < ANG_MAX);
        assert(envs->j_l < ANG_MAX);
        assert(bas(ATOM_OF,i_sh) >= 0);
        assert(bas(ATOM_OF,j_sh) >= 0);
        assert(bas(ATOM_OF,i_sh) < natm);
        assert(bas(ATOM_OF,j_sh) < natm);
        assert(envs->nrys_roots < MXRYSROOTS);

        FINT dli = envs->li_ceil + envs->lj_ceil + 1;
        FINT dlj = envs->lj_ceil + 1;
        envs->g_stride_i = 1;
        envs->g_stride_j = dli;
        envs->g_stride_k = dli * dlj;
        envs->g_size     = dli * dlj;
}

void CINTg1e_index_xyz(FINT *idx, CINTEnvVars *envs)
{
        const FINT i_l = envs->i_l;
        const FINT j_l = envs->j_l;
        const FINT nfi = envs->nfi;
        const FINT nfj = envs->nfj;
        const FINT di = envs->g_stride_i;
        const FINT dj = envs->g_stride_j;
        FINT i, j, n;
        FINT ofx, ofjx;
        FINT ofy, ofjy;
        FINT ofz, ofjz;
        FINT i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
        FINT j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];

        CINTcart_comp(i_nx, i_ny, i_nz, i_l);
        CINTcart_comp(j_nx, j_ny, j_nz, j_l);

        ofx = 0;
        ofy = envs->g_size;
        ofz = envs->g_size * 2;
        n = 0;
        for (j = 0; j < nfj; j++) {
                ofjx = ofx + dj * j_nx[j];
                ofjy = ofy + dj * j_ny[j];
                ofjz = ofz + dj * j_nz[j];
                for (i = 0; i < nfi; i++) {
                        idx[n+0] = ofjx + di * i_nx[i];
                        idx[n+1] = ofjy + di * i_ny[i];
                        idx[n+2] = ofjz + di * i_nz[i];
                        n += 3;
                }
        }
}


void CINTg_ovlp(double *g, double ai, double aj, double fac, CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT dj = envs->g_stride_j;
        const double aij = ai + aj;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        FINT i, j, ptr;
        double rirj[3], ririj[3];
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;

        rirj[0] = ri[0] - rj[0];
        rirj[1] = ri[1] - rj[1];
        rirj[2] = ri[2] - rj[2];
        ririj[0] = ri[0] - (ai * ri[0] + aj * rj[0]) / aij;
        ririj[1] = ri[1] - (ai * ri[1] + aj * rj[1]) / aij;
        ririj[2] = ri[2] - (ai * ri[2] + aj * rj[2]) / aij;

        gx[0] = 1;
        gy[0] = 1;
        gz[0] = SQRTPI * M_PI * fac;
        if (nmax > 0) {
                gx[1] = -ririj[0] * gx[0];
                gy[1] = -ririj[1] * gy[0];
                gz[1] = -ririj[2] * gz[0];
        }

        for (i = 1; i < nmax; i++) {
                gx[i+1] = 0.5 * i / aij * gx[i-1] - ririj[0] * gx[i];
                gy[i+1] = 0.5 * i / aij * gy[i-1] - ririj[1] * gy[i];
                gz[i+1] = 0.5 * i / aij * gz[i-1] - ririj[2] * gz[i];
        }

        for (j = 1; j <= lj; j++) {
                ptr = dj * j;
                for (i = ptr; i <= ptr + nmax - j; i++) {
                        gx[i] = gx[i+1-dj] + rirj[0] * gx[i-dj];
                        gy[i] = gy[i+1-dj] + rirj[1] * gy[i-dj];
                        gz[i] = gz[i+1-dj] + rirj[2] * gz[i-dj];
                }
        }
}

void CINTg_nuc(double *g, double aij, double *rij,
               double *cr, double t2, double fac, CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT dj = envs->g_stride_j;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        FINT i, j, ptr;
        double rir0[3], rirj[3];
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;

        rir0[0] = ri[0] - (rij[0] + t2 * (cr[0] - rij[0]));
        rir0[1] = ri[1] - (rij[1] + t2 * (cr[1] - rij[1]));
        rir0[2] = ri[2] - (rij[2] + t2 * (cr[2] - rij[2]));
        rirj[0] = ri[0] - rj[0];
        rirj[1] = ri[1] - rj[1];
        rirj[2] = ri[2] - rj[2];

        gx[0] = 1;
        gy[0] = 1;
        gz[0] = 2 * M_PI * fac;
        if (nmax > 0) {
                gx[1] = -rir0[0] * gx[0];
                gy[1] = -rir0[1] * gy[0];
                gz[1] = -rir0[2] * gz[0];
        }

        for (i = 1; i < nmax; i++) {
                gx[i+1] = 0.5 * (1 - t2) * i / aij * gx[i-1] - rir0[0] * gx[i];
                gy[i+1] = 0.5 * (1 - t2) * i / aij * gy[i-1] - rir0[1] * gy[i];
                gz[i+1] = 0.5 * (1 - t2) * i / aij * gz[i-1] - rir0[2] * gz[i];
        }

        for (j = 1; j <= lj; j++) {
                ptr = dj * j;
                for (i = ptr; i <= ptr + nmax - j; i++) {
                        gx[i] = gx[i+1-dj] + rirj[0] * gx[i-dj];
                        gy[i] = gy[i+1-dj] + rirj[1] * gy[i-dj];
                        gz[i] = gz[i+1-dj] + rirj[2] * gz[i-dj];
                }
        }
}

void CINTnabla1i_1e(double *f, double *g,
                    FINT li, FINT lj, FINT lk, CINTEnvVars *envs)
{
        const FINT dj = envs->g_stride_j;
        const FINT dk = envs->g_stride_k;
        const double ai2 = -2 * envs->ai;
        FINT i, j, k, ptr;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;

        for (k = 0; k <= lk; k++) {
        for (j = 0; j <= lj; j++) {
                ptr = dj * j + dk * k;
                //f(...,0,...) = -2*ai*g(...,1,...)
                fx[ptr] = ai2 * gx[ptr+1];
                fy[ptr] = ai2 * gy[ptr+1];
                fz[ptr] = ai2 * gz[ptr+1];
                //f(...,i,...) = i*g(...,i-1,...)-2*ai*g(...,i+1,...)
                for (i = 1; i <= li; i++) {
                        fx[ptr+i] = i * gx[ptr+i-1] + ai2 * gx[ptr+i+1];
                        fy[ptr+i] = i * gy[ptr+i-1] + ai2 * gy[ptr+i+1];
                        fz[ptr+i] = i * gz[ptr+i-1] + ai2 * gz[ptr+i+1];
                }
        } }
}

void CINTnabla1j_1e(double *f, double *g,
                    FINT li, FINT lj, FINT lk, CINTEnvVars *envs)
{
        const FINT dj = envs->g_stride_j;
        const FINT dk = envs->g_stride_k;
        const double aj2 = -2 * envs->aj;
        FINT i, j, k, ptr;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;

        for (k = 0; k <= lk; k++) {
                ptr = dk * k;
                //f(...,0,...) = -2*aj*g(...,1,...)
                for (i = ptr; i <= ptr+li; i++) {
                        fx[i] = aj2 * gx[i+dj];
                        fy[i] = aj2 * gy[i+dj];
                        fz[i] = aj2 * gz[i+dj];
                }
                //f(...,j,...) = j*g(...,j-1,...)-2*aj*g(...,j+1,...)
                for (j = 1; j <= lj; j++) {
                        ptr = dj * j + dk * k;
                        for (i = ptr; i <= ptr+li; i++) {
                                fx[i] = j * gx[i-dj] + aj2 * gx[i+dj];
                                fy[i] = j * gy[i-dj] + aj2 * gy[i+dj];
                                fz[i] = j * gz[i-dj] + aj2 * gz[i+dj];
                        }
                }
        }
}

/*
 * ( ij | \nabla k )
 */
void CINTnabla1k_1e(double *f, double *g,
                    FINT li, FINT lj, FINT lk, CINTEnvVars *envs)
{
        const FINT dj = envs->g_stride_j;
        const FINT dk = envs->g_stride_k;
        const double ak2 = -2 * envs->ak;
        FINT i, j, k, ptr;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;

        for (j = 0; j <= lj; j++) {
                ptr = dj * j;
                for (i = ptr; i <= ptr+li; i++) {
                        fx[i] = ak2 * gx[i+dk];
                        fy[i] = ak2 * gy[i+dk];
                        fz[i] = ak2 * gz[i+dk];
                }
        }
        for (k = 1; k <= lk; k++) {
                for (j = 0; j <= lj; j++) {
                        ptr = dj * j + dk * k;
                        for (i = ptr; i <= ptr+li; i++) {
                                fx[i] = k * gx[i-dk] + ak2 * gx[i+dk];
                                fy[i] = k * gy[i-dk] + ak2 * gy[i+dk];
                                fz[i] = k * gz[i-dk] + ak2 * gz[i+dk];
                        }
                }
        }
}


/*
 * ( x^1 i j | k )
 * ri is the shift from the center R_O to the center of |i>
 * r - R_O = (r-R_i) + ri, ri = R_i - R_O
 */
void CINTx1i_1e(double *f, double *g, double ri[3],
                FINT li, FINT lj, FINT lk, CINTEnvVars *envs)
{
        FINT i, j, k, ptr;
        const FINT dj = envs->g_stride_j;
        const FINT dk = envs->g_stride_k;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;

        for (k = 0; k <= lk; k++) {
        for (j = 0; j <= lj; j++) {
                ptr = dj * j + dk * k;
                for (i = ptr; i <= ptr+li; i++) {
                        fx[i] = gx[i+1] + ri[0] * gx[i];
                        fy[i] = gy[i+1] + ri[1] * gy[i];
                        fz[i] = gz[i+1] + ri[2] * gz[i];
                }
        } }
}

void CINTx1j_1e(double *f, double *g, double rj[3],
                FINT li, FINT lj, FINT lk, CINTEnvVars *envs)
{
        FINT i, j, k, ptr;
        const FINT dj = envs->g_stride_j;
        const FINT dk = envs->g_stride_k;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;

        for (k = 0; k <= lk; k++) {
        for (j = 0; j <= lj; j++) {
                ptr = dj * j + dk * k;
                for (i = ptr; i <= ptr+li; i++) {
                        fx[i] = gx[i+dj] + rj[0] * gx[i];
                        fy[i] = gy[i+dj] + rj[1] * gy[i];
                        fz[i] = gz[i+dj] + rj[2] * gz[i];
                }
        } }
}

void CINTx1k_1e(double *f, double *g, double *rk,
                FINT li, FINT lj, FINT lk, CINTEnvVars *envs)
{
        FINT i, j, k, ptr;
        const FINT dj = envs->g_stride_j;
        const FINT dk = envs->g_stride_k;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;

        for (k = 0; k <= lk; k++) {
        for (j = 0; j <= lj; j++) {
                ptr = dj * j + dk * k;
                for (i = ptr; i <= ptr+li; i++) {
                        fx[i] = gx[i+dk] + rk[0] * gx[i];
                        fy[i] = gy[i+dk] + rk[1] * gy[i];
                        fz[i] = gz[i+dk] + rk[2] * gz[i];
                }
        } }
}

/*
 * gc    contracted GTO integral
 * nf    number of primitive integral
 * gp    primitive GTO integral
 * inc   increment of gp
 * shl   nth shell
 * ip    ith-1 primitive GTO
 */
void CINTprim_to_ctr(double *gc, FINT nf, double *gp,
                     FINT inc, FINT nprim, FINT nctr, double *coeff)
{
        FINT n, i, k;
        double *pgc = gc;
        double c;

        for (i = 0; i < inc; i++) {
                //dger(nf, nctr, 1.d0, gp(i+1), inc, env(ptr), nprim, gc(1,i*nctr+1), nf)
                for (n = 0; n < nctr; n++) {
                        c = coeff[nprim*n];
                        if (c != 0) {
                                for (k = 0; k < nf; k++) {
                                        pgc[k] += c * gp[k*inc+i];
                                }
                        }
                        // next cgto block
                        pgc += nf;
                }
        }
}

/* optimized
 * memset(gc, 0, sizeof(double)*nf*nctr);
 * CINTprim_to_ctr(gc, nf, gp, 1, nprim, nprim, nctr, coeff); */
void CINTprim_to_ctr_0(double *gc, FINT nf, double *gp,
                       FINT nprim, FINT nctr, double *coeff)
{
        FINT n, i;
        double c0, c1;
        double *p0, *p1;
        double non0coeff[32];
        double *non0pgc[32];
        FINT ncoeff = 0;

        for (i = 0; i < nctr; i++) {
                if (coeff[nprim*i] != 0) {
                        non0coeff[ncoeff] = coeff[nprim*i];
                        non0pgc[ncoeff] = gc + nf * i;
                        ncoeff++;
                } else { // need to initialize the memory, since += is used later
                        memset(gc+nf*i, 0, sizeof(double)*nf);
                }
        }

        for (i = 0; i < ncoeff-1; i+=2) {
                c0 = non0coeff[i  ];
                c1 = non0coeff[i+1];
                p0 = non0pgc[i  ];
                p1 = non0pgc[i+1];
                for (n = 0; n < nf; n++) {
                        p0[n] = c0 * gp[n];
                        p1[n] = c1 * gp[n];
                }
        }
        if (i < ncoeff) {
                c0 = non0coeff[i];
                p0 = non0pgc[i];
                for (n = 0; n < nf; n++) {
                        p0[n] = c0 * gp[n];
                }
        }
}

/* optimized
 * CINTprim_to_ctr(gc, nf, gp, 1, nprim, nprim, nctr, coeff);
 * with opt->non0coeff, opt->non0idx, opt->non0ctr */
void CINTprim_to_ctr_opt(double *gc, FINT nf, double *gp,
                         double *non0coeff, FINT *non0idx, FINT non0ctr)
{
        FINT n, i;
        double c0, c1;
        double *p0, *p1;

        for (i = 0; i < non0ctr-1; i+=2) {
                c0 = non0coeff[i  ];
                c1 = non0coeff[i+1];
                p0 = gc + nf*non0idx[i  ];
                p1 = gc + nf*non0idx[i+1];
                for (n = 0; n < nf; n++) {
                        p0[n] += c0 * gp[n];
                        p1[n] += c1 * gp[n];
                }
        }
        if (i < non0ctr) {
                c0 = non0coeff[i];
                p0 = gc + nf*non0idx[i];
                for (n = 0; n < nf; n++) {
                        p0[n] += c0 * gp[n];
                }
        }
}

/* optimized
 * CINTprim_to_ctr(gc, nf, gp, 1, nprim, nprim, nctr, coeff); */
void CINTprim_to_ctr_1(double *gc, FINT nf, double *gp,
                       FINT nprim, FINT nctr, double *coeff)
{
        FINT i;
        double non0coeff[32];
        FINT non0idx[32];
        FINT non0ctr = 0;

        for (i = 0; i < nctr; i++) {
                if (coeff[nprim*i] != 0) {
                        non0coeff[non0ctr] = coeff[nprim*i];
                        non0idx[non0ctr] = i;
                        non0ctr++;
                }
        }
        CINTprim_to_ctr_opt(gc, nf, gp, non0coeff, non0idx, non0ctr);
}

/*
 * to optimize memory copy in cart2sph.c, remove the common factor for s
 * and p function in cart2sph
 */
double CINTcommon_fac_sp(FINT l)
{
        switch (l) {
                case 0: return 0.282094791773878143;
                case 1: return 0.488602511902919921;
                default: return 1;
        }
}
