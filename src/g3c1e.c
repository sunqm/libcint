/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <string.h>
#include <math.h>
#include <assert.h>
#include "cint_bas.h"
#include "misc.h"
#include "g1e.h"

void CINTinit_int3c1e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                              FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT k_sh = shls[2];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = 0;
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->x_ctr[2] = bas(NCTR_OF, k_sh);
        envs->x_ctr[3] = 1;
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nfl = 1;
        envs->nf = envs->nfi * envs->nfj * envs->nfk;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = 0;
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = 0;
        envs->nrys_roots =(envs->li_ceil + envs->lj_ceil + envs->lk_ceil)/2 + 1;

        envs->common_factor = SQRTPI * M_PI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l);

        FINT dli = envs->li_ceil + 1;
        FINT dlj = envs->lj_ceil + envs->lk_ceil + 1;
        FINT dlk = envs->lk_ceil + 1;
        envs->g_stride_i = 1;
        envs->g_stride_j = dli;
        envs->g_stride_k = dli * dlj;
        envs->g_stride_l = envs->g_stride_k;
        FINT nmax = envs->li_ceil + dlj;
        envs->g_size     = MAX(dli*dlj*dlk, dli*nmax);

        envs->rirj[0] = envs->ri[0] - envs->rj[0];
        envs->rirj[1] = envs->ri[1] - envs->rj[1];
        envs->rirj[2] = envs->ri[2] - envs->rj[2];
}

void CINTg3c1e_index_xyz(FINT *idx, const CINTEnvVars *envs)
{
        const FINT i_l = envs->i_l;
        const FINT j_l = envs->j_l;
        const FINT k_l = envs->k_l;
        const FINT nfi = envs->nfi;
        const FINT nfj = envs->nfj;
        const FINT nfk = envs->nfk;
        const FINT dj = envs->g_stride_j;
        const FINT dk = envs->g_stride_k;
        FINT i, j, k, n;
        FINT ofx, ofjx, ofkx;
        FINT ofy, ofjy, ofky;
        FINT ofz, ofjz, ofkz;
        FINT i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
        FINT j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];
        FINT k_nx[CART_MAX], k_ny[CART_MAX], k_nz[CART_MAX];

        CINTcart_comp(i_nx, i_ny, i_nz, i_l);
        CINTcart_comp(j_nx, j_ny, j_nz, j_l);
        CINTcart_comp(k_nx, k_ny, k_nz, k_l);

        ofx = 0;
        ofy = envs->g_size;
        ofz = envs->g_size * 2;
        n = 0;
        for (k = 0; k < nfk; k++) {
                ofkx = ofx + dk * k_nx[k];
                ofky = ofy + dk * k_ny[k];
                ofkz = ofz + dk * k_nz[k];
                for (j = 0; j < nfj; j++) {
                        ofjx = ofkx + dj * j_nx[j];
                        ofjy = ofky + dj * j_ny[j];
                        ofjz = ofkz + dj * j_nz[j];
                        for (i = 0; i < nfi; i++) {
                                idx[n+0] = ofjx + i_nx[i];
                                idx[n+1] = ofjy + i_ny[i];
                                idx[n+2] = ofjz + i_nz[i];
                                n += 3;
                        }
                }
        }
}


void CINTg3c1e_ovlp(double *g, double ai, double aj, double ak,
                    double fac, const CINTEnvVars *envs)
{
        const FINT li = envs->li_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT lk = envs->lk_ceil;
        const FINT nmax = li + lj + lk;
        const FINT mmax = lj + lk;
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        gx[0] = 1;
        gy[0] = 1;
        gz[0] = fac;
        if (nmax == 0) {
                return;
        }

        FINT dj = li + 1;
        const FINT dk = envs->g_stride_k;
        const double aijk = ai + aj + ak;
        const double aijk1 = .5 / aijk;
        const double *ri = envs->ri;
        const double *rj = envs->rj;
        const double *rk = envs->rk;
        FINT i, j, k, off;
        const double *rirj = envs->rirj;
        double rjrk[3], rjrijk[3];

        rjrk[0] = rj[0] - rk[0];
        rjrk[1] = rj[1] - rk[1];
        rjrk[2] = rj[2] - rk[2];

        rjrijk[0] = rj[0] - (ai * ri[0] + aj * rj[0] + ak * rk[0]) / aijk;
        rjrijk[1] = rj[1] - (ai * ri[1] + aj * rj[1] + ak * rk[1]) / aijk;
        rjrijk[2] = rj[2] - (ai * ri[2] + aj * rj[2] + ak * rk[2]) / aijk;

        gx[dj] = -rjrijk[0] * gx[0];
        gy[dj] = -rjrijk[1] * gy[0];
        gz[dj] = -rjrijk[2] * gz[0];

        for (j = 1; j < nmax; j++) {
                gx[(j+1)*dj] = aijk1 * j * gx[(j-1)*dj] - rjrijk[0] * gx[j*dj];
                gy[(j+1)*dj] = aijk1 * j * gy[(j-1)*dj] - rjrijk[1] * gy[j*dj];
                gz[(j+1)*dj] = aijk1 * j * gz[(j-1)*dj] - rjrijk[2] * gz[j*dj];
        }

        for (i = 1; i <= li; i++) {
                for (j = 0; j <= nmax-i; j++) { // upper limit lj+lk
                        gx[i+j*dj] = gx[i-1+(j+1)*dj] - rirj[0] * gx[i-1+j*dj];
                        gy[i+j*dj] = gy[i-1+(j+1)*dj] - rirj[1] * gy[i-1+j*dj];
                        gz[i+j*dj] = gz[i-1+(j+1)*dj] - rirj[2] * gz[i-1+j*dj];
                }
        }

        dj = envs->g_stride_j;
        for (k = 1; k <= lk; k++) {
                for (j = 0; j <= mmax-k; j++) {
                        off = k * dk + j * dj;
                        for (i = off; i <= off+li; i++) {
                                gx[i] = gx[i+dj-dk] + rjrk[0] * gx[i-dk];
                                gy[i] = gy[i+dj-dk] + rjrk[1] * gy[i-dk];
                                gz[i] = gz[i+dj-dk] + rjrk[2] * gz[i-dk];
                        }
                }
        }
}

void CINTg3c1e_nuc(double *g, double ai, double aj, double ak, double *rijk,
                   double *cr, double t2, double fac, CINTEnvVars *envs)
{
        const FINT li = envs->li_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT lk = envs->lk_ceil;
        const FINT nmax = li + lj + lk;
        const FINT mmax = lj + lk;
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        gx[0] = 1;
        gy[0] = 1;
        gz[0] = 2/SQRTPI * fac;
        if (nmax == 0) {
                return;
        }

        FINT dj = li + 1;
        const FINT dk = envs->g_stride_k;
        const double aijk = ai + aj + ak;
        const double *rj = envs->rj;
        const double *rk = envs->rk;
        FINT i, j, k, off;
        const double *rirj = envs->rirj;
        double rjrk[3], rjr0[3];

        rjrk[0] = rj[0] - rk[0];
        rjrk[1] = rj[1] - rk[1];
        rjrk[2] = rj[2] - rk[2];

        rjr0[0] = rj[0] - (rijk[0] + t2 * (cr[0] - rijk[0]));
        rjr0[1] = rj[1] - (rijk[1] + t2 * (cr[1] - rijk[1]));
        rjr0[2] = rj[2] - (rijk[2] + t2 * (cr[2] - rijk[2]));

        gx[dj] = -rjr0[0] * gx[0];
        gy[dj] = -rjr0[1] * gy[0];
        gz[dj] = -rjr0[2] * gz[0];

        const double aijk1 = .5 * (1 - t2) / aijk;
        for (j = 1; j < nmax; j++) {
                gx[(j+1)*dj] = aijk1 * j * gx[(j-1)*dj] - rjr0[0] * gx[j*dj];
                gy[(j+1)*dj] = aijk1 * j * gy[(j-1)*dj] - rjr0[1] * gy[j*dj];
                gz[(j+1)*dj] = aijk1 * j * gz[(j-1)*dj] - rjr0[2] * gz[j*dj];
        }

        for (i = 1; i <= li; i++) {
                for (j = 0; j <= nmax-i; j++) { // upper limit lj+lk
                        gx[i+j*dj] = gx[i-1+(j+1)*dj] - rirj[0] * gx[i-1+j*dj];
                        gy[i+j*dj] = gy[i-1+(j+1)*dj] - rirj[1] * gy[i-1+j*dj];
                        gz[i+j*dj] = gz[i-1+(j+1)*dj] - rirj[2] * gz[i-1+j*dj];
                }
        }

        dj = envs->g_stride_j;
        for (k = 1; k <= lk; k++) {
                for (j = 0; j <= mmax-k; j++) {
                        off = k * dk + j * dj;
                        for (i = off; i <= off+li; i++) {
                                gx[i] = gx[i+dj-dk] + rjrk[0] * gx[i-dk];
                                gy[i] = gy[i+dj-dk] + rjrk[1] * gy[i-dk];
                                gz[i] = gz[i+dj-dk] + rjrk[2] * gz[i-dk];
                        }
                }
        }
}

/*
 * ( \nabla i j | k )
 */
void CINTnabla1i_3c1e(double *f, const double *g,
                      const FINT li, const FINT lj, const FINT lk,
                      const CINTEnvVars *envs)
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

/*
 * ( i \nabla j | k )
 */
void CINTnabla1j_3c1e(double *f, const double *g,
                      const FINT li, const FINT lj, const FINT lk,
                      const CINTEnvVars *envs)
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
void CINTnabla1k_3c1e(double *f, const double *g,
                      const FINT li, const FINT lj, const FINT lk,
                      const CINTEnvVars *envs)
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
void CINTx1i_3c1e(double *f, const double *g, const double *ri,
                  const FINT li, const FINT lj, const FINT lk,
                  const CINTEnvVars *envs)
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


/*
 * ( i x^1 j | k )
 */
void CINTx1j_3c1e(double *f, const double *g, const double *rj,
                  const FINT li, const FINT lj, const FINT lk,
                  const CINTEnvVars *envs)
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


/*
 * ( ij | x^1 k )
 */
void CINTx1k_3c1e(double *f, const double *g, const double *rk,
                  const FINT li, const FINT lj, const FINT lk,
                  const CINTEnvVars *envs)
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

