/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cint_const.h"
#include "cint_bas.h"
#include "misc.h"
#include "g3c2e.h"

#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *GX = G; \
        type *GY = G + envs->g_size; \
        type *GZ = G + envs->g_size * 2

void CINTg0_3c2e_kj2d3d(double *g, const CINTEnvVars *envs,struct _BC *bc);
void CINTg0_3c2e_ik2d3d(double *g, const CINTEnvVars *envs,struct _BC *bc);
static void CINTset_g3c2e_params(CINTEnvVars *envs);

/*
 * Note the 3c2e functions takes i,j,k parameters. But we initialize
 * ll_ceil, to reuse g2e_g02d function.
 */

FINT CINTinit_int3c2e_EnvVars(CINTEnvVars *envs, const FINT *ng, const FINT *shls,
                             const FINT *atm, const FINT natm,
                             const FINT *bas, const FINT nbas, const double *env)
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
        envs->i_prim = bas(NPRIM_OF, i_sh);
        envs->j_prim = bas(NPRIM_OF, j_sh);
        envs->k_prim = bas(NPRIM_OF, k_sh);
        envs->i_ctr = bas(NCTR_OF, i_sh);
        envs->j_ctr = bas(NCTR_OF, j_sh);
        envs->k_ctr = bas(NCTR_OF, k_sh);
        envs->nfi = CINTlen_cart(envs->i_l);
        envs->nfj = CINTlen_cart(envs->j_l);
        envs->nfk = CINTlen_cart(envs->k_l);
        envs->nf = envs->nfi * envs->nfk * envs->nfj;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));

        envs->common_factor = (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l);

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = 0; // to reuse CINTg0_2e_2d
        envs->nrys_roots =(envs->li_ceil + envs->lj_ceil
                         + envs->lk_ceil)/2 + 1;

        assert(i_sh < SHLS_MAX);
        assert(j_sh < SHLS_MAX);
        assert(k_sh < SHLS_MAX);
        assert(envs->i_l < ANG_MAX);
        assert(envs->j_l < ANG_MAX);
        assert(envs->k_l < ANG_MAX);
        assert(envs->i_ctr < NCTR_MAX);
        assert(envs->j_ctr < NCTR_MAX);
        assert(envs->k_ctr < NCTR_MAX);
        assert(envs->i_prim < NPRIM_MAX);
        assert(envs->j_prim < NPRIM_MAX);
        assert(envs->k_prim < NPRIM_MAX);
        assert(envs->i_prim >= envs->i_ctr);
        assert(envs->j_prim >= envs->j_ctr);
        assert(envs->k_prim >= envs->k_ctr);
        assert(bas(ATOM_OF,i_sh) >= 0);
        assert(bas(ATOM_OF,j_sh) >= 0);
        assert(bas(ATOM_OF,k_sh) >= 0);
        assert(bas(ATOM_OF,i_sh) < natm);
        assert(bas(ATOM_OF,j_sh) < natm);
        assert(bas(ATOM_OF,k_sh) < natm);
        assert(envs->nrys_roots < MXRYSROOTS);

        CINTset_g3c2e_params(envs);
        return 0;
}

/* set strides and parameters for g0_2d4d algorithm */
static void CINTset_g3c2e_params(CINTEnvVars *envs)
{
        FINT dli, dlj, dlk, dll;
        FINT ibase = envs->li_ceil > envs->lj_ceil;
        FINT kbase = 1;
        if (envs->nrys_roots <= 2) { // use the fully optimized lj_4d algorithm
                ibase = 0;
                kbase = 0;
        }

        if (ibase) {
                dli = envs->li_ceil + envs->lj_ceil + 1;
                dlj = envs->lj_ceil + 1;
        } else {
                dli = envs->li_ceil + 1;
                dlj = envs->li_ceil + envs->lj_ceil + 1;
        }

        if (kbase) {
                dlk = envs->lk_ceil + 1;
                dll = 1;
        } else { // to use _g0_lj_4d_xxxx functions, dll can be > 1
                dlk = envs->lk_ceil + 1;
                dll = dlk;
        }

        envs->g_stride_i = envs->nrys_roots;
        envs->g_stride_k = envs->nrys_roots * dli;
        envs->g_stride_j = envs->nrys_roots * dli * dlk * dll;
        envs->g_size     = envs->nrys_roots * dli * dlk * dll * dlj;

        if (kbase) {
                envs->g2d_klmax = envs->g_stride_k;
                envs->rkrl[0] = envs->rk[0];
                envs->rkrl[1] = envs->rk[1];
                envs->rkrl[2] = envs->rk[2];
                envs->rklrx[0] = 0; // simplify 'envs->rklrx[0] =' in cint2e.c with rl=0, al=0
                envs->rklrx[1] = 0;
                envs->rklrx[2] = 0;
        } else {
                // kbase==0 for calling CINTg0_3c2e_kj2d3d, g2d_klmax will not be used
                envs->g2d_klmax = 0;
                envs->rkrl[0] = -envs->rk[0];
                envs->rkrl[1] = -envs->rk[1];
                envs->rkrl[2] = -envs->rk[2];
                envs->rklrx[0] = envs->rk[0]; // since rl=0, al=0, see cint2e.c
                envs->rklrx[1] = envs->rk[1];
                envs->rklrx[2] = envs->rk[2];
        }
        envs->rkl[0] = envs->rk[0];
        envs->rkl[1] = envs->rk[1];
        envs->rkl[2] = envs->rk[2];

        if (ibase) {
                envs->g2d_ijmax = envs->g_stride_i;
                envs->rx_in_rijrx = envs->ri;
                envs->rirj[0] = envs->ri[0] - envs->rj[0];
                envs->rirj[1] = envs->ri[1] - envs->rj[1];
                envs->rirj[2] = envs->ri[2] - envs->rj[2];
        } else {
                envs->g2d_ijmax = envs->g_stride_j;
                envs->rx_in_rijrx = envs->rj;
                envs->rirj[0] = envs->rj[0] - envs->ri[0];
                envs->rirj[1] = envs->rj[1] - envs->ri[1];
                envs->rirj[2] = envs->rj[2] - envs->ri[2];
        }

        if (ibase) {
                envs->f_g0_2d4d = &CINTg0_3c2e_ik2d3d;
        } else {
                envs->f_g0_2d4d = &CINTg0_3c2e_kj2d3d;
        }
}


void CINTg3c2e_index_xyz(FINT *idx, const CINTEnvVars *envs)
{
        const FINT i_l = envs->i_l;
        const FINT j_l = envs->j_l;
        const FINT k_l = envs->k_l;
        const FINT nfi = envs->nfi;
        const FINT nfj = envs->nfj;
        const FINT nfk = envs->nfk;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dj = envs->g_stride_j;
        FINT i, j, k, n;
        FINT ofx, ofkx;
        FINT ofy, ofky;
        FINT ofz, ofkz;
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
        for (j = 0; j < nfj; j++) {
                for (k = 0; k < nfk; k++) {
                        ofkx = ofx + dj * j_nx[j] + dk * k_nx[k];
                        ofky = ofy + dj * j_ny[j] + dk * k_ny[k];
                        ofkz = ofz + dj * j_nz[j] + dk * k_nz[k];
                        switch (i_l) {
                                case 0:
                                        idx[n+0] = ofkx;
                                        idx[n+1] = ofky;
                                        idx[n+2] = ofkz;
                                        n += 3;
                                        break;
                                case 1:
                                        idx[n+0] = ofkx + di;
                                        idx[n+1] = ofky;
                                        idx[n+2] = ofkz;
                                        idx[n+3] = ofkx;
                                        idx[n+4] = ofky + di;
                                        idx[n+5] = ofkz;
                                        idx[n+6] = ofkx;
                                        idx[n+7] = ofky;
                                        idx[n+8] = ofkz + di;
                                        n += 9;
                                        break;
                                case 2:
                                        idx[n+0 ] = ofkx + di*2;
                                        idx[n+1 ] = ofky;
                                        idx[n+2 ] = ofkz;
                                        idx[n+3 ] = ofkx + di;
                                        idx[n+4 ] = ofky + di;
                                        idx[n+5 ] = ofkz;
                                        idx[n+6 ] = ofkx + di;
                                        idx[n+7 ] = ofky;
                                        idx[n+8 ] = ofkz + di;
                                        idx[n+9 ] = ofkx;
                                        idx[n+10] = ofky + di*2;
                                        idx[n+11] = ofkz;
                                        idx[n+12] = ofkx;
                                        idx[n+13] = ofky + di;
                                        idx[n+14] = ofkz + di;
                                        idx[n+15] = ofkx;
                                        idx[n+16] = ofky;
                                        idx[n+17] = ofkz + di*2;
                                        n += 18;
                                        break;
                                default:
                                        for (i = 0; i < nfi; i++) {
                                                idx[n+0] = ofkx + di * i_nx[i]; //(:,ix,kx,jx,1)
                                                idx[n+1] = ofky + di * i_ny[i]; //(:,iy,ky,jy,2)
                                                idx[n+2] = ofkz + di * i_nz[i]; //(:,iz,kz,jz,3)
                                                n += 3;
                                        } // i
                        }
                } // k
        } // j
}


/*
 * g0[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
/* 2d is based on k,j */
void CINTg0_kj2d_3d(double *g, const CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil;
        const FINT li = envs->li_ceil;
        const FINT nroots = envs->nrys_roots;
        FINT i, j, k, ptr, n;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
        p1x = gx - di;
        p1y = gy - di;
        p1z = gz - di;
        p2x = gx - di + dj;
        p2y = gy - di + dj;
        p2z = gz - di + dj;
        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) {
        for (k = 0; k <= mmax; k++) {
                ptr = j*dj + k*dk + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
                        gx[n] = rirj[0] * p1x[n] + p2x[n];
                        gy[n] = rirj[1] * p1y[n] + p2y[n];
                        gz[n] = rirj[2] * p1z[n] + p2z[n];
                }
        } } }
}

/* to reuse code
void CINTg0_3c2e_lj2d3d(double *g, const CINTEnvVars *envs,struct _BC *bc)
is put in g2e_lj4d.c
*/

void CINTg0_3c2e_ik2d3d(double *g, const CINTEnvVars *envs,struct _BC *bc)
{
        CINTg0_2e_2d(g, bc, envs);

        const FINT lk = envs->lk_ceil;
        const FINT lj = envs->lj_ceil;
        FINT j, k, ptr, n;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        p1x = gx - dj;
        p1y = gy - dj;
        p1z = gz - dj;
        p2x = gx - dj + di;
        p2y = gy - dj + di;
        p2z = gz - dj + di;
        for (j = 1; j <= lj; j++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n++) {
                        gx[n] = rirj[0] * p1x[n] + p2x[n];
                        gy[n] = rirj[1] * p1y[n] + p2y[n];
                        gz[n] = rirj[2] * p1z[n] + p2z[n];
                }
        } }
}


/*
 * ( \nabla i j | k )
 */
void CINTnabla1i_3c2e(double *f, const double *g,
                      const FINT li, const FINT lj, const FINT lk,
                      const CINTEnvVars *envs)
{
        FINT i, j, k, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        const double ai2 = -2 * envs->ai;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx - di;
        const double *p1y = gy - di;
        const double *p1z = gz - di;
        const double *p2x = gx + di;
        const double *p2y = gy + di;
        const double *p2z = gz + di;
        for (j = 0; j <= lj; j++)
        for (k = 0; k <= lk; k++) {
                ptr = dj * j + dk * k;
                //f(...,0,...) = -2*ai*g(...,1,...)
                for (n = ptr; n < ptr+nroots; n++) {
                        fx[n] = ai2 * p2x[n];
                        fy[n] = ai2 * p2y[n];
                        fz[n] = ai2 * p2z[n];
                }
                ptr += di;
                //f(...,i,...) = i*g(...,i-1,...)-2*ai*g(...,i+1,...)
                for (i = 1; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = i*p1x[n] + ai2*p2x[n];
                                fy[n] = i*p1y[n] + ai2*p2y[n];
                                fz[n] = i*p1z[n] + ai2*p2z[n];
                        }
                        ptr += di;
                }
        }
}


/*
 * ( i \nabla j | k )
 */
void CINTnabla1j_3c2e(double *f, const double *g,
                      const FINT li, const FINT lj, const FINT lk,
                      const CINTEnvVars *envs)
{
        FINT i, j, k, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        const double aj2 = -2 * envs->aj;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx - dj;
        const double *p1y = gy - dj;
        const double *p1z = gz - dj;
        const double *p2x = gx + dj;
        const double *p2y = gy + dj;
        const double *p2z = gz + dj;
        //f(...,0,...) = -2*aj*g(...,1,...)
        for (k = 0; k <= lk; k++) {
                ptr = dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = aj2 * p2x[n];
                                fy[n] = aj2 * p2y[n];
                                fz[n] = aj2 * p2z[n];
                        }
                        ptr += di;
                }
        }
        //f(...,j,...) = j*g(...,j-1,...)-2*aj*g(...,j+1,...)
        for (j = 1; j <= lj; j++) {
                for (k = 0; k <= lk; k++) {
                        ptr = dj * j + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
                                        fx[n] = j*p1x[n] + aj2*p2x[n];
                                        fy[n] = j*p1y[n] + aj2*p2y[n];
                                        fz[n] = j*p1z[n] + aj2*p2z[n];
                                }
                                ptr += di;
                        }
                }
        }
}


/*
 * ( ij | \nabla k )
 */
void CINTnabla1k_3c2e(double *f, const double *g,
                      const FINT li, const FINT lj, const FINT lk,
                      const CINTEnvVars *envs)
{
        FINT i, j, k, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        const double ak2 = -2 * envs->ak;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx - dk;
        const double *p1y = gy - dk;
        const double *p1z = gz - dk;
        const double *p2x = gx + dk;
        const double *p2y = gy + dk;
        const double *p2z = gz + dk;
        for (j = 0; j <= lj; j++) {
                ptr = dj * j;
                //f(...,0,...) = -2*ak*g(...,1,...)
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = ak2 * p2x[n];
                                fy[n] = ak2 * p2y[n];
                                fz[n] = ak2 * p2z[n];
                        }
                        ptr += di;
                }
                //f(...,k,...) = k*g(...,k-1,...)-2*ak*g(...,k+1,...)
                for (k = 1; k <= lk; k++) {
                        ptr = dj * j + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
                                        fx[n] = k*p1x[n] + ak2*p2x[n];
                                        fy[n] = k*p1y[n] + ak2*p2y[n];
                                        fz[n] = k*p1z[n] + ak2*p2z[n];
                                }
                                ptr += di;
                        }
                }
        }
}


/*
 * ( x^1 i j | k )
 * ri is the shift from the center R_O to the center of |i>
 * r - R_O = (r-R_i) + ri, ri = R_i - R_O
 */
void CINTx1i_3c2e(double *f, const double *g, const double *ri,
                  const FINT li, const FINT lj, const FINT lk,
                  const CINTEnvVars *envs)
{
        FINT i, j, k, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + di;
        const double *p1y = gy + di;
        const double *p1z = gz + di;
        for (j = 0; j <= lj; j++)
        for (k = 0; k <= lk; k++) {
                //f(...,0:li,...) = g(...,1:li+1,...) + ri(1)*g(...,0:li,...)
                ptr = dj * j + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = p1x[n] + ri[0] * gx[n];
                                fy[n] = p1y[n] + ri[1] * gy[n];
                                fz[n] = p1z[n] + ri[2] * gz[n];
                        }
                        ptr += di;
                }
        }
}


/*
 * ( i x^1 j | k )
 */
void CINTx1j_3c2e(double *f, const double *g, const double *rj,
                  const FINT li, const FINT lj, const FINT lk,
                  const CINTEnvVars *envs)
{
        FINT i, j, k, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + dj;
        const double *p1y = gy + dj;
        const double *p1z = gz + dj;
        for (j = 0; j <= lj; j++)
        for (k = 0; k <= lk; k++) {
                // f(...,0:lj,...) = g(...,1:lj+1,...) + rj(1)*g(...,0:lj,...)
                ptr = dj * j + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = p1x[n] + rj[0] * gx[n];
                                fy[n] = p1y[n] + rj[1] * gy[n];
                                fz[n] = p1z[n] + rj[2] * gz[n];
                        }
                        ptr += di;
                }
        }
}


/*
 * ( ij | x^1 k )
 */
void CINTx1k_3c2e(double *f, const double *g, const double *rk,
                  const FINT li, const FINT lj, const FINT lk,
                  const CINTEnvVars *envs)
{
        FINT i, j, k, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + dk;
        const double *p1y = gy + dk;
        const double *p1z = gz + dk;
        for (j = 0; j <= lj; j++)
        for (k = 0; k <= lk; k++) {
                // f(...,0:lk,...) = g(...,1:lk+1,...) + rk(1)*g(...,0:lk,...)
                ptr = dj * j + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = p1x[n] + rk[0] * gx[n];
                                fy[n] = p1y[n] + rk[1] * gy[n];
                                fz[n] = p1z[n] + rk[2] * gz[n];
                        }
                        ptr += di;
                }
        }
}

