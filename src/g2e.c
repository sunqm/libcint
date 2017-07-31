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
#include "rys_roots.h"
#include "g2e.h"

#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *GX = G; \
        type *GY = G + envs->g_size; \
        type *GZ = G + envs->g_size * 2

void CINTinit_int2e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
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
        const FINT l_sh = shls[3];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = bas(ANG_OF, l_sh);
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->x_ctr[2] = bas(NCTR_OF, k_sh);
        envs->x_ctr[3] = bas(NCTR_OF, l_sh);
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nfl = (envs->l_l+1)*(envs->l_l+2)/2;
        envs->nf = envs->nfi * envs->nfk * envs->nfl * envs->nfj;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));
        envs->rl = env + atm(PTR_COORD, bas(ATOM_OF, l_sh));

        envs->common_factor = (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l) * CINTcommon_fac_sp(envs->l_l);

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = envs->l_l + ng[LINC];
        envs->nrys_roots =(envs->li_ceil + envs->lj_ceil
                         + envs->lk_ceil + envs->ll_ceil)/2 + 1;

        assert(i_sh < SHLS_MAX);
        assert(j_sh < SHLS_MAX);
        assert(k_sh < SHLS_MAX);
        assert(l_sh < SHLS_MAX);
        assert(envs->i_l < ANG_MAX);
        assert(envs->j_l < ANG_MAX);
        assert(envs->k_l < ANG_MAX);
        assert(envs->l_l < ANG_MAX);
        assert(bas(ATOM_OF,i_sh) >= 0);
        assert(bas(ATOM_OF,j_sh) >= 0);
        assert(bas(ATOM_OF,k_sh) >= 0);
        assert(bas(ATOM_OF,l_sh) >= 0);
        assert(bas(ATOM_OF,i_sh) < natm);
        assert(bas(ATOM_OF,j_sh) < natm);
        assert(bas(ATOM_OF,k_sh) < natm);
        assert(bas(ATOM_OF,l_sh) < natm);
        assert(envs->nrys_roots < MXRYSROOTS);

        FINT dli, dlj, dlk, dll;
        FINT ibase = envs->li_ceil > envs->lj_ceil;
        FINT kbase = envs->lk_ceil > envs->ll_ceil;
        if (envs->nrys_roots <= 2) { // use the fully optimized lj_4d algorithm
                ibase = 0;
                kbase = 0;
        }
        if (kbase) {
                dlk = envs->lk_ceil + envs->ll_ceil + 1;
                dll = envs->ll_ceil + 1;
        } else {
                dlk = envs->lk_ceil + 1;
                dll = envs->lk_ceil + envs->ll_ceil + 1;
        }

        if (ibase) {
                dli = envs->li_ceil + envs->lj_ceil + 1;
                dlj = envs->lj_ceil + 1;
        } else {
                dli = envs->li_ceil + 1;
                dlj = envs->li_ceil + envs->lj_ceil + 1;
        }
        envs->g_stride_i = envs->nrys_roots;
        envs->g_stride_k = envs->nrys_roots * dli;
        envs->g_stride_l = envs->nrys_roots * dli * dlk;
        envs->g_stride_j = envs->nrys_roots * dli * dlk * dll;
        envs->g_size     = envs->nrys_roots * dli * dlk * dll * dlj;

        if (kbase) {
                envs->g2d_klmax = envs->g_stride_k;
                envs->rx_in_rklrx = envs->rk;
                envs->rkrl[0] = envs->rk[0] - envs->rl[0];
                envs->rkrl[1] = envs->rk[1] - envs->rl[1];
                envs->rkrl[2] = envs->rk[2] - envs->rl[2];
        } else {
                envs->g2d_klmax = envs->g_stride_l;
                envs->rx_in_rklrx = envs->rl;
                envs->rkrl[0] = envs->rl[0] - envs->rk[0];
                envs->rkrl[1] = envs->rl[1] - envs->rk[1];
                envs->rkrl[2] = envs->rl[2] - envs->rk[2];
        }

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

        if (kbase) {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_ik2d4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_kj2d4d;
                }
        } else {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_il2d4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_lj2d4d;
                }
        }
        envs->f_g0_2e = &CINTg0_2e;
}

void CINTg2e_index_xyz(FINT *idx, const CINTEnvVars *envs)
{
        const FINT i_l = envs->i_l;
        const FINT j_l = envs->j_l;
        const FINT k_l = envs->k_l;
        const FINT l_l = envs->l_l;
        const FINT nfi = envs->nfi;
        const FINT nfj = envs->nfj;
        const FINT nfk = envs->nfk;
        const FINT nfl = envs->nfl;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        FINT i, j, k, l, n;
        FINT ofx, ofkx, oflx;
        FINT ofy, ofky, ofly;
        FINT ofz, ofkz, oflz;
        FINT i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
        FINT j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];
        FINT k_nx[CART_MAX], k_ny[CART_MAX], k_nz[CART_MAX];
        FINT l_nx[CART_MAX], l_ny[CART_MAX], l_nz[CART_MAX];

        CINTcart_comp(i_nx, i_ny, i_nz, i_l);
        CINTcart_comp(j_nx, j_ny, j_nz, j_l);
        CINTcart_comp(k_nx, k_ny, k_nz, k_l);
        CINTcart_comp(l_nx, l_ny, l_nz, l_l);

        ofx = 0;
        ofy = envs->g_size;
        ofz = envs->g_size * 2;
        n = 0;
        for (j = 0; j < nfj; j++) {
                for (l = 0; l < nfl; l++) {
                        oflx = ofx + dj * j_nx[j] + dl * l_nx[l];
                        ofly = ofy + dj * j_ny[j] + dl * l_ny[l];
                        oflz = ofz + dj * j_nz[j] + dl * l_nz[l];
                        for (k = 0; k < nfk; k++) {
                                ofkx = oflx + dk * k_nx[k];
                                ofky = ofly + dk * k_ny[k];
                                ofkz = oflz + dk * k_nz[k];
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
                                                        idx[n+0] = ofkx + di * i_nx[i]; //(:,ix,kx,lx,jx,1)
                                                        idx[n+1] = ofky + di * i_ny[i]; //(:,iy,ky,ly,jy,2)
                                                        idx[n+2] = ofkz + di * i_nz[i]; //(:,iz,kz,lz,jz,3)
                                                        n += 3;
                                                } // i
                                }
                        } // k
                } // l
        } // j
}


/*
 * g(nroots,0:nmax,0:mmax)
 */
void CINTg0_2e_2d(double *g, struct _BC *bc, const CINTEnvVars *envs)
{
        const FINT nroots = envs->nrys_roots;
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        const FINT dm = envs->g2d_klmax;
        const FINT dn = envs->g2d_ijmax;
        FINT i, j, m, n, off;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *c00;
        const double *c0p;
        const double *b01 = bc->b01;
        const double *b00 = bc->b00;
        const double *b10 = bc->b10;
        double *p0x, *p0y, *p0z;
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        for (i = 0; i < nroots; i++) {
                gx[i] = 1;
                gy[i] = 1;
                //gz[i] = w[i];
        }

        if (nmax > 0) {
                p0x = gx + dn;
                p0y = gy + dn;
                p0z = gz + dn;
                p1x = gx - dn;
                p1y = gy - dn;
                p1z = gz - dn;
                // gx(irys,0,1) = c00(irys) * gx(irys,0,0)
                for (c00 = bc->c00, i = 0; i < nroots; i++, c00+=3) {
                        p0x[i] = c00[0] * gx[i];
                        p0y[i] = c00[1] * gy[i];
                        p0z[i] = c00[2] * gz[i];
                }
                // gx(irys,0,n+1) = c00(irys)*gx(irys,0,n)
                // + n*b10(irys)*gx(irys,0,n-1)
                for (n = 1; n < nmax; n++) {
                        off = n * dn;
                        for (c00 = bc->c00, i = 0, j = off;
                             i < nroots; i++, j++, c00+=3) {
                                p0x[j] = c00[0] * gx[j] + n * b10[i] * p1x[j];
                                p0y[j] = c00[1] * gy[j] + n * b10[i] * p1y[j];
                                p0z[j] = c00[2] * gz[j] + n * b10[i] * p1z[j];
                        }
                }
        }

        if (mmax > 0) {
                p0x = gx + dm;
                p0y = gy + dm;
                p0z = gz + dm;
                p1x = gx - dm;
                p1y = gy - dm;
                p1z = gz - dm;
                // gx(irys,1,0) = c0p(irys) * gx(irys,0,0)
                for (c0p = bc->c0p, i = 0; i < nroots; i++, c0p+=3) {
                        p0x[i] = c0p[0] * gx[i];
                        p0y[i] = c0p[1] * gy[i];
                        p0z[i] = c0p[2] * gz[i];
                }
                // gx(irys,m+1,0) = c0p(irys)*gx(irys,m,0)
                // + m*b01(irys)*gx(irys,m-1,0)
                for (m = 1; m < mmax; m++) {
                        off = m * dm;
                        for (c0p = bc->c0p, i = 0, j = off;
                             i < nroots; i++, j++, c0p+=3) {
                                p0x[j] = c0p[0] * gx[j] + m * b01[i] * p1x[j];
                                p0y[j] = c0p[1] * gy[j] + m * b01[i] * p1y[j];
                                p0z[j] = c0p[2] * gz[j] + m * b01[i] * p1z[j];
                        }
                }
        }

        if (nmax > 0 && mmax > 0) {
                p0x = gx + dn;
                p0y = gy + dn;
                p0z = gz + dn;
                p1x = gx - dn;
                p1y = gy - dn;
                p1z = gz - dn;
                p2x = gx - dm;
                p2y = gy - dm;
                p2z = gz - dm;
                // gx(irys,1,1) = c0p(irys)*gx(irys,0,1)
                // + b00(irys)*gx(irys,0,0)
                for (c0p = bc->c0p, i = 0; i < nroots; i++, c0p+=3) {
                        p0x[i+dm] = c0p[0] * p0x[i] + b00[i] * gx[i];
                        p0y[i+dm] = c0p[1] * p0y[i] + b00[i] * gy[i];
                        p0z[i+dm] = c0p[2] * p0z[i] + b00[i] * gz[i];
                }

                // gx(irys,m+1,1) = c0p(irys)*gx(irys,m,1)
                // + m*b01(irys)*gx(irys,m-1,1)
                // + b00(irys)*gx(irys,m,0)
                for (m = 1; m < mmax; m++) {
                        off = m * dm + dn;
                        for (c0p = bc->c0p, i = 0, j = off;
                             i < nroots; i++, j++, c0p+=3) {
                                gx[j+dm] = c0p[0]*gx[j] + m*b01[i]*p2x[j] +b00[i]*p1x[j];
                                gy[j+dm] = c0p[1]*gy[j] + m*b01[i]*p2y[j] +b00[i]*p1y[j];
                                gz[j+dm] = c0p[2]*gz[j] + m*b01[i]*p2z[j] +b00[i]*p1z[j];
                        }
                }

                // gx(irys,m,n+1) = c00(irys)*gx(irys,m,n)
                // + n*b10(irys)*gx(irys,m,n-1)
                // + m*b00(irys)*gx(irys,m-1,n)
                for (m = 1; m <= mmax; m++) {
                        for (n = 1; n < nmax; n++) {
                                off = m * dm + n * dn;
                                for (c00 = bc->c00, i = 0, j = off;
                                     i < nroots; i++, j++, c00+=3) {
                                        p0x[j] = c00[0]*gx[j] +n*b10[i]*p1x[j] + m*b00[i]*p2x[j];
                                        p0y[j] = c00[1]*gy[j] +n*b10[i]*p1y[j] + m*b00[i]*p2y[j];
                                        p0z[j] = c00[2]*gz[j] +n*b10[i]*p1z[j] + m*b00[i]*p2z[j];
                                }
                        }
                }
        }
}


/*
 * g0[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
/* 2d is based on l,j */
void CINTg0_lj2d_4d(double *g, const CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        const FINT li = envs->li_ceil;
        const FINT lk = envs->lk_ceil;
        //const FINT ll = envs->ll_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT nroots = envs->nrys_roots;
        FINT i, j, k, l, ptr, n;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
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
        for (l = 0; l <= mmax; l++) {
                ptr = j*dj + l*dl + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
                        gx[n] = rirj[0] * p1x[n] + p2x[n];
                        gy[n] = rirj[1] * p1y[n] + p2y[n];
                        gz[n] = rirj[2] * p1z[n] + p2z[n];
                }
        } } }

        // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
        p1x = gx - dk;
        p1y = gy - dk;
        p1z = gz - dk;
        p2x = gx - dk + dl;
        p2y = gy - dk + dl;
        p2z = gz - dk + dl;
        for (j = 0; j <= lj; j++) {
        for (k = 1; k <= lk; k++) {
                for (l = 0; l <= mmax-k; l++) {
                        ptr = j*dj + l*dl + k*dk;
                        for (n = ptr; n < ptr+dk; n++) {
                                gx[n] = rkrl[0] * p1x[n] + p2x[n];
                                gy[n] = rkrl[1] * p1y[n] + p2y[n];
                                gz[n] = rkrl[2] * p1z[n] + p2z[n];
                        }
                }
        } }
}
/* 2d is based on k,j */
void CINTg0_kj2d_4d(double *g, const CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        const FINT li = envs->li_ceil;
        //const FINT lk = envs->lk_ceil;
        const FINT ll = envs->ll_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT nroots = envs->nrys_roots;
        FINT i, j, k, l, ptr, n;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
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

        // g(...,k,l,..) = rkrl * g(...,k,l-1,..) + g(...,k+1,l-1,..)
        p1x = gx - dl;
        p1y = gy - dl;
        p1z = gz - dl;
        p2x = gx - dl + dk;
        p2y = gy - dl + dk;
        p2z = gz - dl + dk;
        for (j = 0; j <= lj; j++) {
        for (l = 1; l <= ll; l++) {
        for (k = 0; k <= mmax-l; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk; n++) {
                        gx[n] = rkrl[0] * p1x[n] + p2x[n];
                        gy[n] = rkrl[1] * p1y[n] + p2y[n];
                        gz[n] = rkrl[2] * p1z[n] + p2z[n];
                }
                }
        } }
}
/* 2d is based on i,l */
void CINTg0_il2d_4d(double *g, const CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        //const FINT li = envs->li_ceil;
        const FINT lk = envs->lk_ceil;
        const FINT ll = envs->ll_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT nroots = envs->nrys_roots;
        FINT i, j, k, l, ptr, n;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
        p1x = gx - dk;
        p1y = gy - dk;
        p1z = gz - dk;
        p2x = gx - dk + dl;
        p2y = gy - dk + dl;
        p2z = gz - dk + dl;
        for (k = 1; k <= lk; k++) {
        for (l = 0; l <= mmax-k; l++) {
        for (i = 0; i <= nmax; i++) {
                ptr = l*dl + k*dk + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
                        gx[n] = rkrl[0] * p1x[n] + p2x[n];
                        gy[n] = rkrl[1] * p1y[n] + p2y[n];
                        gz[n] = rkrl[2] * p1z[n] + p2z[n];
                }
        } } }

        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        p1x = gx - dj;
        p1y = gy - dj;
        p1z = gz - dj;
        p2x = gx - dj + di;
        p2y = gy - dj + di;
        p2z = gz - dj + di;
        for (j = 1; j <= lj; j++) {
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n++) {
                        gx[n] = rirj[0] * p1x[n] + p2x[n];
                        gy[n] = rirj[1] * p1y[n] + p2y[n];
                        gz[n] = rirj[2] * p1z[n] + p2z[n];
                }
        } } }
}
/* 2d is based on i,k */
void CINTg0_ik2d_4d(double *g, const CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        //const FINT li = envs->li_ceil;
        const FINT lk = envs->lk_ceil;
        const FINT ll = envs->ll_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT nroots = envs->nrys_roots;
        FINT i, j, k, l, ptr, n;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        // g(...,k,l,..) = rkrl * g(...,k,l-1,..) + g(...,k+1,l-1,..)
        p1x = gx - dl;
        p1y = gy - dl;
        p1z = gz - dl;
        p2x = gx - dl + dk;
        p2y = gy - dl + dk;
        p2z = gz - dl + dk;
        for (l = 1; l <= ll; l++) {
                // (:,i) is full, so loop:k and loop:n can be merged to
                // for(n = l*dl; n < ptr+dl-dk*l; n++)
                for (k = 0; k <= mmax-l; k++) {
                for (i = 0; i <= nmax; i++) {
                        ptr = l*dl + k*dk + i*di;
                        for (n = ptr; n < ptr+nroots; n++) {
                                gx[n] = rkrl[0] * p1x[n] + p2x[n];
                                gy[n] = rkrl[1] * p1y[n] + p2y[n];
                                gz[n] = rkrl[2] * p1z[n] + p2z[n];
                        }
                } }
        }

        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        p1x = gx - dj;
        p1y = gy - dj;
        p1z = gz - dj;
        p2x = gx - dj + di;
        p2y = gy - dj + di;
        p2z = gz - dj + di;
        for (j = 1; j <= lj; j++) {
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n++) {
                        gx[n] = rirj[0] * p1x[n] + p2x[n];
                        gy[n] = rirj[1] * p1y[n] + p2y[n];
                        gz[n] = rirj[2] * p1z[n] + p2z[n];
                }
        } } }
}
/************* some special g0_4d results *************/
/* 4 digits stand for i_ceil, k_ceil, l_ceil, j_ceil */
static inline void _g0_lj_4d_0001(double *g, double *c,
                                  const double *r)
{
        g[0] = 1;
        g[1] = c[0];
        g[2] = 1;
        g[3] = c[1];
        //g[4] = w[0];
        g[5] = c[2] * g[4];
}
static inline void _g0_lj_4d_1000(double *g, double *c,
                                  const double *r)
{
        g[0] = 1;
        g[1] = r[0] + c[0];
        g[4] = 1;
        g[5] = r[1] + c[1];
        //g[8] = w[0];
        g[9] =(r[2] + c[2]) * g[8];
}
static inline void _g0_lj_4d_0002(double *g, double *c, double *b,
                                  const double *r)
{
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = c[0];
        g[3 ] = c[3];
        g[4 ] = c[0] * c[0] + b[0];
        g[5 ] = c[3] * c[3] + b[1];
        g[6 ] = 1;
        g[7 ] = 1;
        g[8 ] = c[1];
        g[9 ] = c[4];
        g[10] = c[1] * c[1] + b[0];
        g[11] = c[4] * c[4] + b[1];
        //g[12] = w[0];
        //g[13] = w[1];
        g[14] = c[2] * g[12];
        g[15] = c[5] * g[13];
        g[16] =(c[2] * c[2] + b[0])* g[12];
        g[17] =(c[5] * c[5] + b[1])* g[13];
}
static inline void _g0_lj_4d_1001(double *g, double *c, double *b,
                                  const double *r)
{
        double rc[] = {r[0]+c[0], r[0]+c[3],
                       r[1]+c[1], r[1]+c[4],
                       r[2]+c[2], r[2]+c[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = c[0];
        g[5 ] = c[3];
        g[6 ] = rc[0] * c[0] + b[0];
        g[7 ] = rc[1] * c[3] + b[1];
        g[12] = 1;
        g[13] = 1;
        g[14] = rc[2];
        g[15] = rc[3];
        g[16] = c[1];
        g[17] = c[4];
        g[18] = rc[2] * c[1] + b[0];
        g[19] = rc[3] * c[4] + b[1];
        //g[24] = w[0];
        //g[25] = w[1];
        g[26] = rc[4] * g[24];
        g[27] = rc[5] * g[25];
        g[28] = c[2] * g[24];
        g[29] = c[5] * g[25];
        g[30] =(rc[4] * c[2] + b[0])* g[24];
        g[31] =(rc[5] * c[5] + b[1])* g[25];
}
static inline void _g0_lj_4d_2000(double *g, double *c, double *b,
                                  const double *r)
{
        double rc[] = {r[0]+c[0], r[0]+c[3],
                       r[1]+c[1], r[1]+c[4],
                       r[2]+c[2], r[2]+c[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = rc[0] * rc[0] + b[0];
        g[5 ] = rc[1] * rc[1] + b[1];
        g[18] = 1;
        g[19] = 1;
        g[20] = rc[2];
        g[21] = rc[3];
        g[22] = rc[2] * rc[2] + b[0];
        g[23] = rc[3] * rc[3] + b[1];
        //g[36] = w[0];
        //g[37] = w[1];
        g[38] = rc[4] * g[36];
        g[39] = rc[5] * g[37];
        g[40] =(rc[4] * rc[4] + b[0])* g[36];
        g[41] =(rc[5] * rc[5] + b[1])* g[37];
}
static inline void _g0_lj_4d_0003(double *g, double *c, double *b,
                                  const double *r)
{
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = c[0];
        g[3 ] = c[3];
        g[4 ] = c[0] * c[0] + b[0];
        g[5 ] = c[3] * c[3] + b[1];
        g[6 ] = c[0] *(c[0] * c[0] + 3 * b[0]);
        g[7 ] = c[3] *(c[3] * c[3] + 3 * b[1]);
        g[8 ] = 1;
        g[9 ] = 1;
        g[10] = c[1];
        g[11] = c[4];
        g[12] = c[1] * c[1] + b[0];
        g[13] = c[4] * c[4] + b[1];
        g[14] = c[1] *(c[1] * c[1] + 3 * b[0]);
        g[15] = c[4] *(c[4] * c[4] + 3 * b[1]);
        //g[16] = w[0];
        //g[17] = w[1];
        g[18] = c[2] * g[16];
        g[19] = c[5] * g[17];
        g[20] =(c[2] * c[2] + b[0])* g[16];
        g[21] =(c[5] * c[5] + b[1])* g[17];
        g[22] =(c[2] * c[2] + 3 * b[0])* c[2] * g[16];
        g[23] =(c[5] * c[5] + 3 * b[1])* c[5] * g[17];
}
static inline void _g0_lj_4d_1002(double *g, double *c, double *b,
                                  const double *r)
{
        double rc[] = {r[0]+c[0], r[0]+c[3],
                       r[1]+c[1], r[1]+c[4],
                       r[2]+c[2], r[2]+c[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = c[0];
        g[5 ] = c[3];
        g[6 ] = rc[0] * c[0] + b[0];
        g[7 ] = rc[1] * c[3] + b[1];
        g[8 ] = c[0] * c[0] + b[0];
        g[9 ] = c[3] * c[3] + b[1];
        g[10] = rc[0]*c[0]*c[0] + b[0]*(rc[0]+2*c[0]);
        g[11] = rc[1]*c[3]*c[3] + b[1]*(rc[1]+2*c[3]);
        g[16] = 1;
        g[17] = 1;
        g[18] = rc[2];
        g[19] = rc[3];
        g[20] = c[1];
        g[21] = c[4];
        g[22] = rc[2] * c[1] + b[0];
        g[23] = rc[3] * c[4] + b[1];
        g[24] = c[1] * c[1] + b[0];
        g[25] = c[4] * c[4] + b[1];
        g[26] = rc[2]*c[1]*c[1] + b[0]*(rc[2]+2*c[1]);
        g[27] = rc[3]*c[4]*c[4] + b[1]*(rc[3]+2*c[4]);
        //g[32] = w[0];
        //g[33] = w[1];
        g[34] = rc[4] * g[32];
        g[35] = rc[5] * g[33];
        g[36] = c[2] * g[32];
        g[37] = c[5] * g[33];
        g[38] =(rc[4] * c[2] + b[0])* g[32];
        g[39] =(rc[5] * c[5] + b[1])* g[33];
        g[40] =(c[2] * c[2] + b[0])* g[32];
        g[41] =(c[5] * c[5] + b[1])* g[33];
        g[42] =(rc[4]*c[2]*c[2]+b[0]*(rc[4]+2*c[2]))*g[32];
        g[43] =(rc[5]*c[5]*c[5]+b[1]*(rc[5]+2*c[5]))*g[33];
}
static inline void _g0_lj_4d_2001(double *g, double *c, double *b,
                                  const double *r)
{
        double rc[] = {r[0]+c[0], r[0]+c[3],
                       r[1]+c[1], r[1]+c[4],
                       r[2]+c[2], r[2]+c[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = rc[0] * rc[0] + b[0];
        g[5 ] = rc[1] * rc[1] + b[1];
        g[6 ] = c[0];
        g[7 ] = c[3];
        g[8 ] = c[0] * rc[0] + b[0];
        g[9 ] = c[3] * rc[1] + b[1];
        g[10] = c[0]*rc[0]*rc[0] + b[0]*(2*rc[0]+c[0]);
        g[11] = c[3]*rc[1]*rc[1] + b[1]*(2*rc[1]+c[3]);
        g[24] = 1;
        g[25] = 1;
        g[26] = rc[2];
        g[27] = rc[3];
        g[28] = rc[2] * rc[2] + b[0];
        g[29] = rc[3] * rc[3] + b[1];
        g[30] = c[1];
        g[31] = c[4];
        g[32] = c[1] * rc[2] + b[0];
        g[33] = c[4] * rc[3] + b[1];
        g[34] = c[1]*rc[2]*rc[2] + b[0]*(2*rc[2]+c[1]);
        g[35] = c[4]*rc[3]*rc[3] + b[1]*(2*rc[3]+c[4]);
        //g[48] = w[0];
        //g[49] = w[1];
        g[50] = rc[4] * g[48];
        g[51] = rc[5] * g[49];
        g[52] =(rc[4] * rc[4] + b[0])* g[48];
        g[53] =(rc[5] * rc[5] + b[1])* g[49];
        g[54] = c[2] * g[48];
        g[55] = c[5] * g[49];
        g[56] =(c[2] * rc[4] + b[0])* g[48];
        g[57] =(c[5] * rc[5] + b[1])* g[49];
        g[58] =(c[2]*rc[4]*rc[4] + b[0]*(2*rc[4]+c[2]))* g[48];
        g[59] =(c[5]*rc[5]*rc[5] + b[1]*(2*rc[5]+c[5]))* g[49];
}
static inline void _g0_lj_4d_3000(double *g, double *c, double *b,
                                  const double *r)
{
        double rc[] = {r[0]+c[0], r[0]+c[3],
                       r[1]+c[1], r[1]+c[4],
                       r[2]+c[2], r[2]+c[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = rc[0] * rc[0] + b[0];
        g[5 ] = rc[1] * rc[1] + b[1];
        g[6 ] = rc[0] *(rc[0] * rc[0] + 3 * b[0]);
        g[7 ] = rc[1] *(rc[1] * rc[1] + 3 * b[1]);
        g[32] = 1;
        g[33] = 1;
        g[34] = rc[2];
        g[35] = rc[3];
        g[36] = rc[2] * rc[2] + b[0];
        g[37] = rc[3] * rc[3] + b[1];
        g[38] = rc[2] *(rc[2] * rc[2] + 3 * b[0]);
        g[39] = rc[3] *(rc[3] * rc[3] + 3 * b[1]);
        //g[64] = w[0];
        //g[65] = w[1];
        g[66] = rc[4] * g[64];
        g[67] = rc[5] * g[65];
        g[68] =(rc[4] * rc[4] + b[0])* g[64];
        g[69] =(rc[5] * rc[5] + b[1])* g[65];
        g[70] =(rc[4] * rc[4] + 3 * b[0])* rc[4] * g[64];
        g[71] =(rc[5] * rc[5] + 3 * b[1])* rc[5] * g[65];
}
static inline void _g0_lj_4d_0011(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp)
{
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = cp[0];
        g[3 ] = cp[3];
        g[4 ] = c0[0];
        g[5 ] = c0[3];
        g[6 ] = cp[0] * c0[0] + b[0];
        g[7 ] = cp[3] * c0[3] + b[1];
        g[8 ] = 1;
        g[9 ] = 1;
        g[10] = cp[1];
        g[11] = cp[4];
        g[12] = c0[1];
        g[13] = c0[4];
        g[14] = cp[1] * c0[1] + b[0];
        g[15] = cp[4] * c0[4] + b[1];
        //g[16] = w[0];
        //g[17] = w[1];
        g[18] = cp[2] * g[16];
        g[19] = cp[5] * g[17];
        g[20] = c0[2] * g[16];
        g[21] = c0[5] * g[17];
        g[22] =(cp[2] * c0[2] + b[0]) * g[16];
        g[23] =(cp[5] * c0[5] + b[1]) * g[17];
}
static inline void _g0_lj_4d_1010(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp)
{
        double rc[] = {r0[0]+c0[0], r0[0]+c0[3],
                       r0[1]+c0[1], r0[1]+c0[4],
                       r0[2]+c0[2], r0[2]+c0[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = cp[0];
        g[5 ] = cp[3];
        g[6 ] = rc[0] * cp[0] + b[0];
        g[7 ] = rc[1] * cp[3] + b[1];
        g[16] = 1;
        g[17] = 1;
        g[18] = rc[2];
        g[19] = rc[3];
        g[20] = cp[1];
        g[21] = cp[4];
        g[22] = rc[2] * cp[1] + b[0];
        g[23] = rc[3] * cp[4] + b[1];
        //g[32] = w[0];
        //g[33] = w[1];
        g[34] = rc[4] * g[32];
        g[35] = rc[5] * g[33];
        g[36] = cp[2] * g[32];
        g[37] = cp[5] * g[33];
        g[38] =(rc[4]*cp[2] + b[0]) * g[32];
        g[39] =(rc[5]*cp[5] + b[1]) * g[33];
}
static inline void _g0_lj_4d_0101(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp)
{
        double rc[] = {rp[0]+cp[0], rp[0]+cp[3],
                       rp[1]+cp[1], rp[1]+cp[4],
                       rp[2]+cp[2], rp[2]+cp[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[8 ] = c0[0];
        g[9 ] = c0[3];
        g[10] = rc[0] * c0[0] + b[0];
        g[11] = rc[1] * c0[3] + b[1];
        g[16] = 1;
        g[17] = 1;
        g[18] = rc[2];
        g[19] = rc[3];
        g[24] = c0[1];
        g[25] = c0[4];
        g[26] = rc[2] * c0[1] + b[0];
        g[27] = rc[3] * c0[4] + b[1];
        //g[32] = w[0];
        //g[33] = w[1];
        g[34] = rc[4] * g[32];
        g[35] = rc[5] * g[33];
        g[40] = c0[2] * g[32];
        g[41] = c0[5] * g[33];
        g[42] =(rc[4]*c0[2] + b[0]) * g[32];
        g[43] =(rc[5]*c0[5] + b[1]) * g[33];
}
static inline void _g0_lj_4d_1100(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp)
{
        double rc0[] = {r0[0]+c0[0], r0[0]+c0[3],
                        r0[1]+c0[1], r0[1]+c0[4],
                        r0[2]+c0[2], r0[2]+c0[5]};
        double rcp[] = {rp[0]+cp[0], rp[0]+cp[3],
                        rp[1]+cp[1], rp[1]+cp[4],
                        rp[2]+cp[2], rp[2]+cp[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rcp[0];
        g[5 ] = rcp[1];
        g[6 ] = rc0[0] * rcp[0] + b[0];
        g[7 ] = rc0[1] * rcp[1] + b[1];
        g[32] = 1;
        g[33] = 1;
        g[34] = rc0[2];
        g[35] = rc0[3];
        g[36] = rcp[2];
        g[37] = rcp[3];
        g[38] = rc0[2] * rcp[2] + b[0];
        g[39] = rc0[3] * rcp[3] + b[1];
        //g[64] = w[0];
        //g[65] = w[1];
        g[66] = rc0[4] * g[64];
        g[67] = rc0[5] * g[65];
        g[68] = rcp[4] * g[64];
        g[69] = rcp[5] * g[65];
        g[70] =(rc0[4]*rcp[4] + b[0]) * g[64];
        g[71] =(rc0[5]*rcp[5] + b[1]) * g[65];
}
static inline void _g0_lj_4d_0021(double *g, double *c0, double *cp,
                                  double *b0, double *b1)
{
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = cp[0];
        g[3 ] = cp[3];
        g[4 ] = cp[0] * cp[0] + b1[0];
        g[5 ] = cp[3] * cp[3] + b1[1];
        g[6 ] = c0[0];
        g[7 ] = c0[3];
        g[8 ] = cp[0] * c0[0] + b0[0];
        g[9 ] = cp[3] * c0[3] + b0[1];
        g[10] = c0[0] * g[4] + 2 * b0[0] * cp[0];
        g[11] = c0[3] * g[5] + 2 * b0[1] * cp[3];
        g[12] = 1;
        g[13] = 1;
        g[14] = cp[1];
        g[15] = cp[4];
        g[16] = cp[1] * cp[1] + b1[0];
        g[17] = cp[4] * cp[4] + b1[1];
        g[18] = c0[1];
        g[19] = c0[4];
        g[20] = cp[1] * c0[1] + b0[0];
        g[21] = cp[4] * c0[4] + b0[1];
        g[22] = c0[1] * g[16] + 2 * b0[0] * cp[1];
        g[23] = c0[4] * g[17] + 2 * b0[1] * cp[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[26] = cp[2] * g[24];
        g[27] = cp[5] * g[25];
        g[28] =(cp[2] * cp[2] + b1[0]) * g[24];
        g[29] =(cp[5] * cp[5] + b1[1]) * g[25];
        g[30] = c0[2] * g[24];
        g[31] = c0[5] * g[25];
        g[32] =(cp[2] * c0[2] + b0[0]) * g[24];
        g[33] =(cp[5] * c0[5] + b0[1]) * g[25];
        g[34] = c0[2] * g[28] + 2 * b0[0] * g[26];
        g[35] = c0[5] * g[29] + 2 * b0[1] * g[27];
}
static inline void _g0_lj_4d_1020(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {r0[0]+c0[0], r0[0]+c0[3],
                       r0[1]+c0[1], r0[1]+c0[4],
                       r0[2]+c0[2], r0[2]+c0[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = cp[0];
        g[5 ] = cp[3];
        g[6 ] = cp[0] * rc[0] + b0[0];
        g[7 ] = cp[3] * rc[1] + b0[1];
        g[8 ] = cp[0] * cp[0] + b1[0];
        g[9 ] = cp[3] * cp[3] + b1[1];
        g[10] = rc[0] * g[8] + 2 * b0[0] * cp[0];
        g[11] = rc[1] * g[9] + 2 * b0[1] * cp[3];
        g[24] = 1;
        g[25] = 1;
        g[26] = rc[2];
        g[27] = rc[3];
        g[28] = cp[1];
        g[29] = cp[4];
        g[30] = cp[1] * rc[2] + b0[0];
        g[31] = cp[4] * rc[3] + b0[1];
        g[32] = cp[1] * cp[1] + b1[0];
        g[33] = cp[4] * cp[4] + b1[1];
        g[34] = rc[2] * g[32] + 2 * b0[0] * cp[1];
        g[35] = rc[3] * g[33] + 2 * b0[1] * cp[4];
        //g[48] = w[0];
        //g[49] = w[1];
        g[50] = rc[4] * g[48];
        g[51] = rc[5] * g[49];
        g[52] = cp[2] * g[48];
        g[53] = cp[5] * g[49];
        g[54] =(cp[2] * rc[4] + b0[0]) * g[48];
        g[55] =(cp[5] * rc[5] + b0[1]) * g[49];
        g[56] =(cp[2] * cp[2] + b1[0]) * g[48];
        g[57] =(cp[5] * cp[5] + b1[1]) * g[49];
        g[58] = rc[4] * g[56] + 2 * b0[0] * g[52];
        g[59] = rc[5] * g[57] + 2 * b0[1] * g[53];
}
static inline void _g0_lj_4d_0111(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {rp[0]+cp[0], rp[0]+cp[3],
                       rp[1]+cp[1], rp[1]+cp[4],
                       rp[2]+cp[2], rp[2]+cp[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = cp[0];
        g[5 ] = cp[3];
        g[6 ] = cp[0] * rc[0] + b1[0];
        g[7 ] = cp[3] * rc[1] + b1[1];
        g[12] = c0[0];
        g[13] = c0[3];
        g[14] = c0[0] * rc[0] + b0[0];
        g[15] = c0[3] * rc[1] + b0[1];
        g[16] = c0[0] * cp[0] + b0[0];
        g[17] = c0[3] * cp[3] + b0[1];
        g[18] = c0[0] * g[6] + b0[0] *(rc[0] + cp[0]);
        g[19] = c0[3] * g[7] + b0[1] *(rc[1] + cp[3]);
        g[24] = 1;
        g[25] = 1;
        g[26] = rc[2];
        g[27] = rc[3];
        g[28] = cp[1];
        g[29] = cp[4];
        g[30] = cp[1] * rc[2] + b1[0];
        g[31] = cp[4] * rc[3] + b1[1];
        g[36] = c0[1];
        g[37] = c0[4];
        g[38] = c0[1] * rc[2] + b0[0];
        g[39] = c0[4] * rc[3] + b0[1];
        g[40] = c0[1] * cp[1] + b0[0];
        g[41] = c0[4] * cp[4] + b0[1];
        g[42] = c0[1] * g[30] + b0[0] *(rc[2] + cp[1]);
        g[43] = c0[4] * g[31] + b0[1] *(rc[3] + cp[4]);
        //g[48] = w[0];
        //g[49] = w[1];
        g[50] = rc[4] * g[48];
        g[51] = rc[5] * g[49];
        g[52] = cp[2] * g[48];
        g[53] = cp[5] * g[49];
        g[54] =(cp[2] * rc[4] + b1[0]) * g[48];
        g[55] =(cp[5] * rc[5] + b1[1]) * g[49];
        g[60] = c0[2] * g[48];
        g[61] = c0[5] * g[49];
        g[62] =(c0[2] * rc[4] + b0[0]) * g[48];
        g[63] =(c0[5] * rc[5] + b0[1]) * g[49];
        g[64] =(c0[2] * cp[2] + b0[0]) * g[48];
        g[65] =(c0[5] * cp[5] + b0[1]) * g[49];
        g[66] = c0[2] * g[54] + b0[0] *(g[50] + g[52]);
        g[67] = c0[5] * g[55] + b0[1] *(g[51] + g[53]);
}
static inline void _g0_lj_4d_1110(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc0[] = {r0[0]+c0[0], r0[0]+c0[3],
                        r0[1]+c0[1], r0[1]+c0[4],
                        r0[2]+c0[2], r0[2]+c0[5]};
        double rcp[] = {rp[0]+cp[0], rp[0]+cp[3],
                        rp[1]+cp[1], rp[1]+cp[4],
                        rp[2]+cp[2], rp[2]+cp[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rcp[0];
        g[5 ] = rcp[1];
        g[6 ] = rcp[0] * rc0[0] + b0[0];
        g[7 ] = rcp[1] * rc0[1] + b0[1];
        g[8 ] = cp[0];
        g[9 ] = cp[3];
        g[10] = cp[0] * rc0[0] + b0[0];
        g[11] = cp[3] * rc0[1] + b0[1];
        g[12] = cp[0] * rcp[0] + b1[0];
        g[13] = cp[3] * rcp[1] + b1[1];
        g[14] = rc0[0] * g[12] + b0[0] *(rcp[0] + cp[0]);
        g[15] = rc0[1] * g[13] + b0[1] *(rcp[1] + cp[3]);
        g[48] = 1;
        g[49] = 1;
        g[50] = rc0[2];
        g[51] = rc0[3];
        g[52] = rcp[2];
        g[53] = rcp[3];
        g[54] = rcp[2] * rc0[2] + b0[0];
        g[55] = rcp[3] * rc0[3] + b0[1];
        g[56] = cp[1];
        g[57] = cp[4];
        g[58] = cp[1] * rc0[2] + b0[0];
        g[59] = cp[4] * rc0[3] + b0[1];
        g[60] = cp[1] * rcp[2] + b1[0];
        g[61] = cp[4] * rcp[3] + b1[1];
        g[62] = rc0[2] * g[60] + b0[0] *(rcp[2] + cp[1]);
        g[63] = rc0[3] * g[61] + b0[1] *(rcp[3] + cp[4]);
        //g[96 ] = w[0];
        //g[97 ] = w[1];
        g[98 ] = rc0[4] * g[96 ];
        g[99 ] = rc0[5] * g[97 ];
        g[100] = rcp[4] * g[96 ];
        g[101] = rcp[5] * g[97 ];
        g[102] =(rcp[4] * rc0[4] + b0[0])* g[96 ];
        g[103] =(rcp[5] * rc0[5] + b0[1])* g[97 ];
        g[104] = cp[2] * g[96 ];
        g[105] = cp[5] * g[97 ];
        g[106] =(cp[2] * rc0[4] + b0[0])* g[96 ];
        g[107] =(cp[5] * rc0[5] + b0[1])* g[97 ];
        g[108] =(cp[2] * rcp[4] + b1[0])* g[96 ];
        g[109] =(cp[5] * rcp[5] + b1[1])* g[97 ];
        g[110] = rc0[4] * g[108] + b0[0] *(g[100] + g[104]);
        g[111] = rc0[5] * g[109] + b0[1] *(g[101] + g[105]);
}
static inline void _g0_lj_4d_0201(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {rp[0]+cp[0], rp[0]+cp[3],
                       rp[1]+cp[1], rp[1]+cp[4],
                       rp[2]+cp[2], rp[2]+cp[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = rc[0] * rc[0] + b1[0];
        g[5 ] = rc[1] * rc[1] + b1[1];
        g[18] = c0[0];
        g[19] = c0[3];
        g[20] = rc[0] * c0[0] + b0[0];
        g[21] = rc[1] * c0[3] + b0[1];
        g[22] = c0[0] * g[4] + 2 * b0[0] * rc[0];
        g[23] = c0[3] * g[5] + 2 * b0[1] * rc[1];
        g[36] = 1;
        g[37] = 1;
        g[38] = rc[2];
        g[39] = rc[3];
        g[40] = rc[2] * rc[2] + b1[0];
        g[41] = rc[3] * rc[3] + b1[1];
        g[54] = c0[1];
        g[55] = c0[4];
        g[56] = rc[2] * c0[1] + b0[0];
        g[57] = rc[3] * c0[4] + b0[1];
        g[58] = c0[1] * g[40] + 2 * b0[0] * rc[2];
        g[59] = c0[4] * g[41] + 2 * b0[1] * rc[3];
        //g[72] = w[0];
        //g[73] = w[1];
        g[74] = rc[4] * g[72];
        g[75] = rc[5] * g[73];
        g[76] =(rc[4] * rc[4] + b1[0])* g[72];
        g[77] =(rc[5] * rc[5] + b1[1])* g[73];
        g[90] = c0[2] * g[72];
        g[91] = c0[5] * g[73];
        g[92] =(rc[4] * c0[2] + b0[0])* g[72];
        g[93] =(rc[5] * c0[5] + b0[1])* g[73];
        g[94] = c0[2] * g[76] + 2 * b0[0] * g[74];
        g[95] = c0[5] * g[77] + 2 * b0[1] * g[75];
}
static inline void _g0_lj_4d_1200(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc0[] = {r0[0]+c0[0], r0[0]+c0[3],
                        r0[1]+c0[1], r0[1]+c0[4],
                        r0[2]+c0[2], r0[2]+c0[5]};
        double rcp[] = {rp[0]+cp[0], rp[0]+cp[3],
                        rp[1]+cp[1], rp[1]+cp[4],
                        rp[2]+cp[2], rp[2]+cp[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rcp[0];
        g[5 ] = rcp[1];
        g[6 ] = rcp[0] * rc0[0] + b0[0];
        g[7 ] = rcp[1] * rc0[1] + b0[1];
        g[8 ] = rcp[0] * rcp[0] + b1[0];
        g[9 ] = rcp[1] * rcp[1] + b1[1];
        g[10] = rc0[0] * g[8] + 2 * b0[0] * rcp[0];
        g[11] = rc0[1] * g[9] + 2 * b0[1] * rcp[1];
        g[72] = 1;
        g[73] = 1;
        g[74] = rc0[2];
        g[75] = rc0[3];
        g[76] = rcp[2];
        g[77] = rcp[3];
        g[78] = rcp[2] * rc0[2] + b0[0];
        g[79] = rcp[3] * rc0[3] + b0[1];
        g[80] = rcp[2] * rcp[2] + b1[0];
        g[81] = rcp[3] * rcp[3] + b1[1];
        g[82] = rc0[2] * g[80] + 2 * b0[0] * rcp[2];
        g[83] = rc0[3] * g[81] + 2 * b0[1] * rcp[3];
        //g[144] = w[0];
        //g[145] = w[1];
        g[146] = rc0[4] * g[144];
        g[147] = rc0[5] * g[145];
        g[148] = rcp[4] * g[144];
        g[149] = rcp[5] * g[145];
        g[150] =(rcp[4] * rc0[4] + b0[0])* g[144];
        g[151] =(rcp[5] * rc0[5] + b0[1])* g[145];
        g[152] =(rcp[4] * rcp[4] + b1[0])* g[144];
        g[153] =(rcp[5] * rcp[5] + b1[1])* g[145];
        g[154] = rc0[4] * g[152] + 2 * b0[0] * g[148];
        g[155] = rc0[5] * g[153] + 2 * b0[1] * g[149];
}
static inline void _g0_lj_4d_0012(double *g, double *c0, double *cp,
                                  double *b0, double *b1)
{
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = cp[0];
        g[3 ] = cp[3];
        g[4 ] = c0[0];
        g[5 ] = c0[3];
        g[6 ] = cp[0] * c0[0] + b0[0];
        g[7 ] = cp[3] * c0[3] + b0[1];
        g[8 ] = c0[0] * c0[0] + b1[0];
        g[9 ] = c0[3] * c0[3] + b1[1];
        g[10] = cp[0] * g[8] + 2 * b0[0] * c0[0];
        g[11] = cp[3] * g[9] + 2 * b0[1] * c0[3];
        g[12] = 1;
        g[13] = 1;
        g[14] = cp[1];
        g[15] = cp[4];
        g[16] = c0[1];
        g[17] = c0[4];
        g[18] = cp[1] * c0[1] + b0[0];
        g[19] = cp[4] * c0[4] + b0[1];
        g[20] = c0[1] * c0[1] + b1[0];
        g[21] = c0[4] * c0[4] + b1[1];
        g[22] = cp[1] * g[20] + 2 * b0[0] * c0[1];
        g[23] = cp[4] * g[21] + 2 * b0[1] * c0[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[26] = cp[2] * g[24];
        g[27] = cp[5] * g[25];
        g[28] = c0[2] * g[24];
        g[29] = c0[5] * g[25];
        g[30] =(cp[2] * c0[2] + b0[0]) * g[24];
        g[31] =(cp[5] * c0[5] + b0[1]) * g[25];
        g[32] =(c0[2] * c0[2] + b1[0]) * g[24];
        g[33] =(c0[5] * c0[5] + b1[1]) * g[25];
        g[34] = cp[2] * g[32] + 2 * b0[0] * g[28];
        g[35] = cp[5] * g[33] + 2 * b0[1] * g[29];
}
static inline void _g0_lj_4d_1011(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {r0[0]+c0[0], r0[0]+c0[3],
                       r0[1]+c0[1], r0[1]+c0[4],
                       r0[2]+c0[2], r0[2]+c0[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = cp[0];
        g[5 ] = cp[3];
        g[6 ] = cp[0] * rc[0] + b0[0];
        g[7 ] = cp[3] * rc[1] + b0[1];
        g[8 ] = c0[0];
        g[9 ] = c0[3];
        g[10] = c0[0] * rc[0] + b1[0];
        g[11] = c0[3] * rc[1] + b1[1];
        g[12] = c0[0] * cp[0] + b0[0];
        g[13] = c0[3] * cp[3] + b0[1];
        g[14] = cp[0] * g[10] + b0[0] *(rc[0] + c0[0]);
        g[15] = cp[3] * g[11] + b0[1] *(rc[1] + c0[3]);
        g[24] = 1;
        g[25] = 1;
        g[26] = rc[2];
        g[27] = rc[3];
        g[28] = cp[1];
        g[29] = cp[4];
        g[30] = cp[1] * rc[2] + b0[0];
        g[31] = cp[4] * rc[3] + b0[1];
        g[32] = c0[1];
        g[33] = c0[4];
        g[34] = c0[1] * rc[2] + b1[0];
        g[35] = c0[4] * rc[3] + b1[1];
        g[36] = c0[1] * cp[1] + b0[0];
        g[37] = c0[4] * cp[4] + b0[1];
        g[38] = cp[1] * g[34] + b0[0] *(rc[2] + c0[1]);
        g[39] = cp[4] * g[35] + b0[1] *(rc[3] + c0[4]);
        //g[48] = w[0];
        //g[49] = w[1];
        g[50] = rc[4] * g[48];
        g[51] = rc[5] * g[49];
        g[52] = cp[2] * g[48];
        g[53] = cp[5] * g[49];
        g[54] =(cp[2] * rc[4] + b0[0])* g[48];
        g[55] =(cp[5] * rc[5] + b0[1])* g[49];
        g[56] = c0[2] * g[48];
        g[57] = c0[5] * g[49];
        g[58] =(c0[2] * rc[4] + b1[0])* g[48];
        g[59] =(c0[5] * rc[5] + b1[1])* g[49];
        g[60] =(c0[2] * cp[2] + b0[0])* g[48];
        g[61] =(c0[5] * cp[5] + b0[1])* g[49];
        g[62] = cp[2] * g[58] + b0[0] *(g[50] + g[56]);
        g[63] = cp[5] * g[59] + b0[1] *(g[51] + g[57]);
}
static inline void _g0_lj_4d_2010(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {r0[0]+c0[0], r0[0]+c0[3],
                       r0[1]+c0[1], r0[1]+c0[4],
                       r0[2]+c0[2], r0[2]+c0[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[4 ] = rc[0] * rc[0] + b1[0];
        g[5 ] = rc[1] * rc[1] + b1[1];
        g[6 ] = cp[0];
        g[7 ] = cp[3];
        g[8 ] = cp[0] * rc[0] + b0[0];
        g[9 ] = cp[3] * rc[1] + b0[1];
        g[10] = cp[0] * g[4] + 2 * b0[0] * rc[0];
        g[11] = cp[3] * g[5] + 2 * b0[1] * rc[1];
        g[36] = 1;
        g[37] = 1;
        g[38] = rc[2];
        g[39] = rc[3];
        g[40] = rc[2] * rc[2] + b1[0];
        g[41] = rc[3] * rc[3] + b1[1];
        g[42] = cp[1];
        g[43] = cp[4];
        g[44] = cp[1] * rc[2] + b0[0];
        g[45] = cp[4] * rc[3] + b0[1];
        g[46] = cp[1] * g[40] + 2 * b0[0] * rc[2];
        g[47] = cp[4] * g[41] + 2 * b0[1] * rc[3];
        //g[72] = w[0];
        //g[73] = w[1];
        g[74] = rc[4] * g[72];
        g[75] = rc[5] * g[73];
        g[76] =(rc[4] * rc[4] + b1[0]) * g[72];
        g[77] =(rc[5] * rc[5] + b1[1]) * g[73];
        g[78] = cp[2] * g[72];
        g[79] = cp[5] * g[73];
        g[80] =(cp[2] * rc[4] + b0[0]) * g[72];
        g[81] =(cp[5] * rc[5] + b0[1]) * g[73];
        g[82] = cp[2] * g[76] + 2 * b0[0] * g[74];
        g[83] = cp[5] * g[77] + 2 * b0[1] * g[75];
}
static inline void _g0_lj_4d_0102(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc[] = {rp[0]+cp[0], rp[0]+cp[3],
                       rp[1]+cp[1], rp[1]+cp[4],
                       rp[2]+cp[2], rp[2]+cp[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc[0];
        g[3 ] = rc[1];
        g[8 ] = c0[0];
        g[9 ] = c0[3];
        g[10] = rc[0] * c0[0] + b0[0];
        g[11] = rc[1] * c0[3] + b0[1];
        g[16] = c0[0] * c0[0] + b1[0];
        g[17] = c0[3] * c0[3] + b1[1];
        g[18] = rc[0] * g[16] + 2 * b0[0] * c0[0];
        g[19] = rc[1] * g[17] + 2 * b0[1] * c0[3];
        g[24] = 1;
        g[25] = 1;
        g[26] = rc[2];
        g[27] = rc[3];
        g[32] = c0[1];
        g[33] = c0[4];
        g[34] = rc[2] * c0[1] + b0[0];
        g[35] = rc[3] * c0[4] + b0[1];
        g[40] = c0[1] * c0[1] + b1[0];
        g[41] = c0[4] * c0[4] + b1[1];
        g[42] = rc[2] * g[40] + 2 * b0[0] * c0[1];
        g[43] = rc[3] * g[41] + 2 * b0[1] * c0[4];
        //g[48] = w[0];
        //g[49] = w[1];
        g[50] = rc[4] * g[48];
        g[51] = rc[5] * g[49];
        g[56] = c0[2] * g[48];
        g[57] = c0[5] * g[49];
        g[58] =(rc[4] * c0[2] + b0[0]) * g[48];
        g[59] =(rc[5] * c0[5] + b0[1]) * g[49];
        g[64] =(c0[2] * c0[2] + b1[0]) * g[48];
        g[65] =(c0[5] * c0[5] + b1[1]) * g[49];
        g[66] = rc[4] * g[64] + 2 * b0[0] * g[56];
        g[67] = rc[5] * g[65] + 2 * b0[1] * g[57];
}
static inline void _g0_lj_4d_1101(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc0[] = {r0[0]+c0[0], r0[0]+c0[3],
                        r0[1]+c0[1], r0[1]+c0[4],
                        r0[2]+c0[2], r0[2]+c0[5]};
        double rcp[] = {rp[0]+cp[0], rp[0]+cp[3],
                        rp[1]+cp[1], rp[1]+cp[4],
                        rp[2]+cp[2], rp[2]+cp[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rcp[0];
        g[5 ] = rcp[1];
        g[6 ] = rcp[0] * rc0[0] + b0[0];
        g[7 ] = rcp[1] * rc0[1] + b0[1];
        g[16] = c0[0];
        g[17] = c0[3];
        g[18] = c0[0] * rc0[0] + b1[0];
        g[19] = c0[3] * rc0[1] + b1[1];
        g[20] = c0[0] * rcp[0] + b0[0];
        g[21] = c0[3] * rcp[1] + b0[1];
        g[22] = rcp[0] * g[18] + b0[0] *(rc0[0] + c0[0]);
        g[23] = rcp[1] * g[19] + b0[1] *(rc0[1] + c0[3]);
        g[48] = 1;
        g[49] = 1;
        g[50] = rc0[2];
        g[51] = rc0[3];
        g[52] = rcp[2];
        g[53] = rcp[3];
        g[54] = rcp[2] * rc0[2] + b0[0];
        g[55] = rcp[3] * rc0[3] + b0[1];
        g[64] = c0[1];
        g[65] = c0[4];
        g[66] = c0[1] * rc0[2] + b1[0];
        g[67] = c0[4] * rc0[3] + b1[1];
        g[68] = c0[1] * rcp[2] + b0[0];
        g[69] = c0[4] * rcp[3] + b0[1];
        g[70] = rcp[2] * g[66] + b0[0] *(rc0[2] + c0[1]);
        g[71] = rcp[3] * g[67] + b0[1] *(rc0[3] + c0[4]);
        //g[96 ] = w[0];
        //g[97 ] = w[1];
        g[98 ] = rc0[4] * g[96];
        g[99 ] = rc0[5] * g[97];
        g[100] = rcp[4] * g[96];
        g[101] = rcp[5] * g[97];
        g[102] =(rcp[4] * rc0[4] + b0[0])* g[96];
        g[103] =(rcp[5] * rc0[5] + b0[1])* g[97];
        g[112] = c0[2] * g[96];
        g[113] = c0[5] * g[97];
        g[114] =(c0[2] * rc0[4] + b1[0])* g[96];
        g[115] =(c0[5] * rc0[5] + b1[1])* g[97];
        g[116] =(c0[2] * rcp[4] + b0[0])* g[96];
        g[117] =(c0[5] * rcp[5] + b0[1])* g[97];
        g[118] = rcp[4] * g[114] + b0[0] *(g[98] + g[112]);
        g[119] = rcp[5] * g[115] + b0[1] *(g[99] + g[113]);
}
static inline void _g0_lj_4d_2100(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp)
{
        double rc0[] = {r0[0]+c0[0], r0[0]+c0[3],
                        r0[1]+c0[1], r0[1]+c0[4],
                        r0[2]+c0[2], r0[2]+c0[5]};
        double rcp[] = {rp[0]+cp[0], rp[0]+cp[3],
                        rp[1]+cp[1], rp[1]+cp[4],
                        rp[2]+cp[2], rp[2]+cp[5]};
        g[0 ] = 1;
        g[1 ] = 1;
        g[2 ] = rc0[0];
        g[3 ] = rc0[1];
        g[4 ] = rc0[0] * rc0[0] + b1[0];
        g[5 ] = rc0[1] * rc0[1] + b1[1];
        g[6 ] = rcp[0];
        g[7 ] = rcp[1];
        g[8 ] = rcp[0] * rc0[0] + b0[0];
        g[9 ] = rcp[1] * rc0[1] + b0[1];
        g[10] = rcp[0] * g[4] + 2 * b0[0] * rc0[0];
        g[11] = rcp[1] * g[5] + 2 * b0[1] * rc0[1];
        g[72] = 1;
        g[73] = 1;
        g[74] = rc0[2];
        g[75] = rc0[3];
        g[76] = rc0[2] * rc0[2] + b1[0];
        g[77] = rc0[3] * rc0[3] + b1[1];
        g[78] = rcp[2];
        g[79] = rcp[3];
        g[80] = rcp[2] * rc0[2] + b0[0];
        g[81] = rcp[3] * rc0[3] + b0[1];
        g[82] = rcp[2] * g[76] + 2 * b0[0] * rc0[2];
        g[83] = rcp[3] * g[77] + 2 * b0[1] * rc0[3];
        //g[144] = w[0];
        //g[145] = w[1];
        g[146] = rc0[4] * g[144];
        g[147] = rc0[5] * g[145];
        g[148] =(rc0[4] * rc0[4] + b1[0])* g[144];
        g[149] =(rc0[5] * rc0[5] + b1[1])* g[145];
        g[150] = rcp[4] * g[144];
        g[151] = rcp[5] * g[145];
        g[152] =(rcp[4] * rc0[4] + b0[0])* g[144];
        g[153] =(rcp[5] * rc0[5] + b0[1])* g[145];
        g[154] = rcp[4] * g[148] + 2 * b0[0] * g[146];
        g[155] = rcp[5] * g[149] + 2 * b0[1] * g[147];
}
/************** end special g0_4d results *************/



void CINTg0_2e_lj2d4d(double *g, struct _BC *bc, const CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        switch (nmax) {
                case 0: switch(mmax) {
                        case 0: goto _g0_4d_default; // ssss
                        case 1: switch (envs->lk_ceil) {
                                case 0: _g0_lj_4d_0001(g, bc->c0p, envs->rkrl); goto normal_end;
                                case 1: _g0_lj_4d_1000(g, bc->c0p, envs->rkrl); goto normal_end;
                                default: goto error; }
                        case 2: switch (envs->lk_ceil) {
                                case 0: _g0_lj_4d_0002(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                case 1: _g0_lj_4d_1001(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                case 2: _g0_lj_4d_2000(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                default: goto error; }
                        case 3: switch (envs->lk_ceil) {
                                case 0: _g0_lj_4d_0003(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                case 1: _g0_lj_4d_1002(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                case 2: _g0_lj_4d_2001(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                case 3: _g0_lj_4d_3000(g, bc->c0p, bc->b01, envs->rkrl); goto normal_end;
                                default: goto error; }
                        default: goto _g0_4d_default; }
                case 1: switch(mmax) {
                        case 0: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0001(g, bc->c00, envs->rirj); goto normal_end;
                                case 1: _g0_lj_4d_1000(g, bc->c00, envs->rirj); goto normal_end;
                                default: goto error; }
                        case 1: switch (envs->lk_ceil) {
                                case 0: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0011(g, bc->c00, bc->c0p, bc->b00, envs->rirj, envs->rkrl); goto normal_end;
                                        case 1: _g0_lj_4d_1010(g, bc->c00, bc->c0p, bc->b00, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                case 1: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0101(g, bc->c00, bc->c0p, bc->b00, envs->rirj, envs->rkrl); goto normal_end;
                                        case 1: _g0_lj_4d_1100(g, bc->c00, bc->c0p, bc->b00, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                default: goto error; }
                        case 2: switch (envs->lk_ceil) {
                                case 0: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0021(g, bc->c00, bc->c0p, bc->b00, bc->b01); goto normal_end;
                                        case 1: _g0_lj_4d_1020(g, bc->c00, bc->c0p, bc->b00, bc->b01, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                case 1: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0111(g, bc->c00, bc->c0p, bc->b00, bc->b01, envs->rirj, envs->rkrl); goto normal_end;
                                        case 1: _g0_lj_4d_1110(g, bc->c00, bc->c0p, bc->b00, bc->b01, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                case 2: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0201(g, bc->c00, bc->c0p, bc->b00, bc->b01, envs->rirj, envs->rkrl); goto normal_end;
                                        case 1: _g0_lj_4d_1200(g, bc->c00, bc->c0p, bc->b00, bc->b01, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                default: goto error; }
                        default: goto _g0_4d_default; }
                case 2: switch(mmax) {
                        case 0: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0002(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                case 1: _g0_lj_4d_1001(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                case 2: _g0_lj_4d_2000(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                default: goto error; }
                        case 1: switch (envs->lk_ceil) {
                                case 0: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0012(g, bc->c00, bc->c0p, bc->b00, bc->b10); goto normal_end;
                                        case 1: _g0_lj_4d_1011(g, bc->c00, bc->c0p, bc->b00, bc->b10, envs->rirj, envs->rkrl); goto normal_end;
                                        case 2: _g0_lj_4d_2010(g, bc->c00, bc->c0p, bc->b00, bc->b10, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                case 1: switch (envs->li_ceil) {
                                        case 0: _g0_lj_4d_0102(g, bc->c00, bc->c0p, bc->b00, bc->b10, envs->rirj, envs->rkrl); goto normal_end;
                                        case 1: _g0_lj_4d_1101(g, bc->c00, bc->c0p, bc->b00, bc->b10, envs->rirj, envs->rkrl); goto normal_end;
                                        case 2: _g0_lj_4d_2100(g, bc->c00, bc->c0p, bc->b00, bc->b10, envs->rirj, envs->rkrl); goto normal_end;
                                        default: goto error; }
                                default: goto error; }
                        default: goto _g0_4d_default; }
                case 3: switch(mmax) {
                        case 0: switch (envs->li_ceil) {
                                case 0: _g0_lj_4d_0003(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                case 1: _g0_lj_4d_1002(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                case 2: _g0_lj_4d_2001(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                case 3: _g0_lj_4d_3000(g, bc->c00, bc->b10, envs->rirj); goto normal_end;
                                default: goto error; }
                        default: goto _g0_4d_default; }
                default:
_g0_4d_default:
                        CINTg0_2e_2d(g, bc, envs);
                        CINTg0_lj2d_4d(g, envs);
        }
normal_end:
        return;
error:
        fprintf(stderr, "Dimension error for CINTg0_2e_lj2d4d: iklj = %d %d %d %d",
               (int)envs->li_ceil, (int)envs->lk_ceil,
               (int)envs->ll_ceil, (int)envs->lj_ceil);
        exit(1);
}

void CINTg0_2e_kj2d4d(double *g, struct _BC *bc, const CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_kj2d_4d(g, envs);
}
void CINTg0_2e_ik2d4d(double *g, struct _BC *bc, const CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_ik2d_4d(g, envs);
}
void CINTg0_2e_il2d4d(double *g, struct _BC *bc, const CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_il2d_4d(g, envs);
}

#ifdef WITH_F12
void CINTg0_2e_stg_lj2d4d(double *g, struct _BC *bc, const CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_lj2d_4d(g, envs);
}
#endif

/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
void CINTg0_2e(double *g, const double fac, const CINTEnvVars *envs)
{
        const double aij = envs->aij;
        const double akl = envs->akl;
        double a0, a1, fac1, x;
        double u[MXRYSROOTS];
        double *w = g + envs->g_size * 2; // ~ gz
        double rijrkl[3];
        rijrkl[0] = envs->rij[0] - envs->rkl[0];
        rijrkl[1] = envs->rij[1] - envs->rkl[1];
        rijrkl[2] = envs->rij[2] - envs->rkl[2];

        a1 = aij * akl;
        a0 = a1 / (aij + akl);

#ifdef WITH_RANGE_COULOMB
        const double omega = envs->env[PTR_RANGE_OMEGA];
        double theta = 0;
        if (omega > 0) {
// For long-range part of range-separated Coulomb operator
                theta = omega * omega / (omega * omega + a0);
                a0 *= theta;
        }
#endif

        fac1 = sqrt(a0 / (a1 * a1 * a1)) * fac;
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);
        CINTrys_roots(envs->nrys_roots, x, u, w);

        FINT irys;
#ifdef WITH_RANGE_COULOMB
        if (omega > 0) {
                /* u[:] = tau^2 / (1 - tau^2)
                 * omega^2u^2 = a0 * tau^2 / (theta^-1 - tau^2)
                 * transform u[:] to theta^-1 tau^2 / (theta^-1 - tau^2)
                 * so the rest code can be reused.
                 */
                for (irys = 0; irys < envs->nrys_roots; irys++) {
                        u[irys] /= u[irys] + 1 - u[irys] * theta;
                }
        }
#endif
        if (envs->g_size == 1) {
                g[0] = 1;
                g[1] = 1;
                g[2] *= fac1;
                return;
        }

        double u2, div, tmp1, tmp2, tmp3, tmp4;
        const double *rijrx = envs->rijrx;
        const double *rklrx = envs->rklrx;
        struct _BC bc;
        double *c00 = bc.c00;
        double *c0p = bc.c0p;

        for (irys = 0; irys < envs->nrys_roots; irys++, c00+=3, c0p+=3)
        {
                /*
                 *u(irys) = t2/(1-t2)
                 *t2 = u(irys)/(1+u(irys))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = a0 * u[irys];
                div = 1 / (u2 * (aij + akl) + a1);
                tmp1 = u2 * div;
                tmp2 = tmp1 * akl;
                tmp3 = tmp1 * aij;
                tmp4 = .5 * div;
                bc.b00[irys] = 0.5 * tmp1;
                bc.b10[irys] = bc.b00[irys] + tmp4 * akl;
                bc.b01[irys] = bc.b00[irys] + tmp4 * aij;
                c00[0] = rijrx[0] - tmp2 * rijrkl[0];
                c00[1] = rijrx[1] - tmp2 * rijrkl[1];
                c00[2] = rijrx[2] - tmp2 * rijrkl[2];
                c0p[0] = rklrx[0] + tmp3 * rijrkl[0];
                c0p[1] = rklrx[1] + tmp3 * rijrkl[1];
                c0p[2] = rklrx[2] + tmp3 * rijrkl[2];
                w[irys] *= fac1;
        }

        (*envs->f_g0_2d4d)(g, &bc, envs);
}

/*
 * ( \nabla i j | kl )
 */
void CINTnabla1i_2e(double *f, const double *g,
                    const FINT li, const FINT lj, const FINT lk, const FINT ll,
                    const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
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
        for (l = 0; l <= ll; l++)
        for (k = 0; k <= lk; k++) {
                ptr = dj * j + dl * l + dk * k;
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
 * ( i \nabla j | kl )
 */
void CINTnabla1j_2e(double *f, const double *g,
                    const FINT li, const FINT lj, const FINT lk, const FINT ll,
                    const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
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
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = dl * l + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = aj2 * p2x[n];
                                fy[n] = aj2 * p2y[n];
                                fz[n] = aj2 * p2z[n];
                        }
                        ptr += di;
                }
        } }
        //f(...,j,...) = j*g(...,j-1,...)-2*aj*g(...,j+1,...)
        for (j = 1; j <= lj; j++) {
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        ptr = dj * j + dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
                                        fx[n] = j*p1x[n] + aj2*p2x[n];
                                        fy[n] = j*p1y[n] + aj2*p2y[n];
                                        fz[n] = j*p1z[n] + aj2*p2z[n];
                                }
                                ptr += di;
                        }
                } }
        }
}


/*
 * ( ij | \nabla k l )
 */
void CINTnabla1k_2e(double *f, const double *g,
                    const FINT li, const FINT lj, const FINT lk, const FINT ll,
                    const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
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
        for (j = 0; j <= lj; j++)
        for (l = 0; l <= ll; l++) {
                ptr = dj * j + dl * l;
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
                        ptr = dj * j + dl * l + dk * k;
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
 * ( ij | k \nabla l )
 */
void CINTnabla1l_2e(double *f, const double *g,
                    const FINT li, const FINT lj, const FINT lk, const FINT ll,
                    const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        const double al2 = -2 * envs->al;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx - dl;
        const double *p1y = gy - dl;
        const double *p1z = gz - dl;
        const double *p2x = gx + dl;
        const double *p2y = gy + dl;
        const double *p2z = gz + dl;
        for (j = 0; j <= lj; j++) {
                //f(...,0,...) = -2*al*g(...,1,...)
                for (k = 0; k <= lk; k++) {
                        ptr = dj * j + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
                                        fx[n] = al2 * p2x[n];
                                        fy[n] = al2 * p2y[n];
                                        fz[n] = al2 * p2z[n];
                                }
                                ptr += di;
                        }
                }
                //f(...,l,...) = l*g(...,l-1,...)-2*al*g(...,l+1,...)
                for (l = 1; l <= ll; l++) {
                        for (k = 0; k <= lk; k++) {
                                ptr = dj * j + dl * l + dk * k;
                                for (i = 0; i <= li; i++, ptr += di) {
                                for (n = ptr; n < ptr+nroots; n++) {
                                        fx[n] = l*p1x[n] + al2*p2x[n];
                                        fy[n] = l*p1y[n] + al2*p2y[n];
                                        fz[n] = l*p1z[n] + al2*p2z[n];
                                } }
                        }
                }
        }
}

/*
 * ( x^1 i j | kl )
 * ri is the shift from the center R_O to the center of |i>
 * r - R_O = (r-R_i) + ri, ri = R_i - R_O
 */
void CINTx1i_2e(double *f, const double *g, const double *ri,
                const FINT li, const FINT lj, const FINT lk, const FINT ll,
                const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + di;
        const double *p1y = gy + di;
        const double *p1z = gz + di;
        for (j = 0; j <= lj; j++)
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                //f(...,0:li,...) = g(...,1:li+1,...) + ri(1)*g(...,0:li,...)
                ptr = dj * j + dl * l + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = p1x[n] + ri[0] * gx[n];
                                fy[n] = p1y[n] + ri[1] * gy[n];
                                fz[n] = p1z[n] + ri[2] * gz[n];
                        }
                        ptr += di;
                }
        } }
}


/*
 * ( i x^1 j | kl )
 */
void CINTx1j_2e(double *f, const double *g, const double *rj,
                const FINT li, const FINT lj, const FINT lk, const FINT ll,
                const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + dj;
        const double *p1y = gy + dj;
        const double *p1z = gz + dj;
        for (j = 0; j <= lj; j++)
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                // f(...,0:lj,...) = g(...,1:lj+1,...) + rj(1)*g(...,0:lj,...)
                ptr = dj * j + dl * l + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = p1x[n] + rj[0] * gx[n];
                                fy[n] = p1y[n] + rj[1] * gy[n];
                                fz[n] = p1z[n] + rj[2] * gz[n];
                        }
                        ptr += di;
                }
        } }
}


/*
 * ( ij | x^1 k l )
 */
void CINTx1k_2e(double *f, const double *g, const double *rk,
                const FINT li, const FINT lj, const FINT lk, const FINT ll,
                const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + dk;
        const double *p1y = gy + dk;
        const double *p1z = gz + dk;
        for (j = 0; j <= lj; j++)
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                // f(...,0:lk,...) = g(...,1:lk+1,...) + rk(1)*g(...,0:lk,...)
                ptr = dj * j + dl * l + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = p1x[n] + rk[0] * gx[n];
                                fy[n] = p1y[n] + rk[1] * gy[n];
                                fz[n] = p1z[n] + rk[2] * gz[n];
                        }
                        ptr += di;
                }
        } }
}


/*
 * ( i j | x^1 kl )
 */
void CINTx1l_2e(double *f, const double *g, const double *rl,
                const FINT li, const FINT lj, const FINT lk, const FINT ll,
                const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + dl;
        const double *p1y = gy + dl;
        const double *p1z = gz + dl;
        for (j = 0; j <= lj; j++)
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                // f(...,0:ll,...) = g(...,1:ll+1,...) + rl(1)*g(...,0:ll,...)
                ptr = dj * j + dl * l + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = p1x[n] + rl[0] * gx[n];
                                fy[n] = p1y[n] + rl[1] * gy[n];
                                fz[n] = p1z[n] + rl[2] * gz[n];
                        }
                        ptr += di;
                }
        } }
}

