/*
 * File: g2e.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cint_bas.h"
#include "misc.h"
#include "g2e.h"

#define MIN(X,Y)        ((X) < (Y) ? (X) : (Y))
#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *GX = G; \
        type *GY = G + envs->g_size; \
        type *GZ = G + envs->g_size * 2

typedef struct {
        double c00[MXRYSROOTS*3];
        double c0p[MXRYSROOTS*3];
        double b01[MXRYSROOTS];
        double b00[MXRYSROOTS];
        double b10[MXRYSROOTS];
} _BC;


int init_int2e_CintEnvVars(CintEnvVars *envs, const unsigned int ng[],
                           const unsigned int *shls,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;
        envs->ng = ng;

        const unsigned int i_sh = shls[0];
        const unsigned int j_sh = shls[1];
        const unsigned int k_sh = shls[2];
        const unsigned int l_sh = shls[3];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = bas(ANG_OF, l_sh);
        envs->i_prim = bas(NPRIM_OF, i_sh);
        envs->j_prim = bas(NPRIM_OF, j_sh);
        envs->k_prim = bas(NPRIM_OF, k_sh);
        envs->l_prim = bas(NPRIM_OF, l_sh);
        envs->i_ctr = bas(NCTR_OF, i_sh);
        envs->j_ctr = bas(NCTR_OF, j_sh);
        envs->k_ctr = bas(NCTR_OF, k_sh);
        envs->l_ctr = bas(NCTR_OF, l_sh);
        envs->nfi = len_cart(envs->i_l);
        envs->nfj = len_cart(envs->j_l);
        envs->nfk = len_cart(envs->k_l);
        envs->nfl = len_cart(envs->l_l);
        envs->nf = envs->nfi * envs->nfk * envs->nfl * envs->nfj;
        envs->g_size = ng[RYS_ROOTS] * ng[0] * ng[1] * ng[2] * ng[3];
        envs->g_stride_i = ng[RYS_ROOTS]; // shift of (i++,k,l,j)
        envs->g_stride_k = ng[RYS_ROOTS] * ng[0]; // shift of (i,k++,l,j)
        envs->g_stride_l = ng[RYS_ROOTS] * ng[0] * ng[1]; // shift of (i,k,l++,j)
        envs->g_stride_j = ng[RYS_ROOTS] * ng[0] * ng[1] * ng[2]; // shift of (i,k,l,j++)
        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));
        envs->rl = env + atm(PTR_COORD, bas(ATOM_OF, l_sh));

        assert(i_sh < SHLS_MAX);
        assert(j_sh < SHLS_MAX);
        assert(k_sh < SHLS_MAX);
        assert(l_sh < SHLS_MAX);
        assert(envs->i_l < ANG_MAX);
        assert(envs->j_l < ANG_MAX);
        assert(envs->k_l < ANG_MAX);
        assert(envs->l_l < ANG_MAX);
        assert(envs->i_ctr < NCTR_MAX);
        assert(envs->j_ctr < NCTR_MAX);
        assert(envs->k_ctr < NCTR_MAX);
        assert(envs->l_ctr < NCTR_MAX);
        assert(envs->i_prim < NPRIM_MAX);
        assert(envs->j_prim < NPRIM_MAX);
        assert(envs->k_prim < NPRIM_MAX);
        assert(envs->l_prim < NPRIM_MAX);
        assert(envs->i_prim >= envs->i_ctr);
        assert(envs->j_prim >= envs->j_ctr);
        assert(envs->k_prim >= envs->k_ctr);
        assert(envs->l_prim >= envs->l_ctr);
        assert(bas(ATOM_OF,i_sh) >= 0);
        assert(bas(ATOM_OF,j_sh) >= 0);
        assert(bas(ATOM_OF,k_sh) >= 0);
        assert(bas(ATOM_OF,l_sh) >= 0);
        assert(bas(ATOM_OF,i_sh) < natm);
        assert(bas(ATOM_OF,j_sh) < natm);
        assert(bas(ATOM_OF,k_sh) < natm);
        assert(bas(ATOM_OF,l_sh) < natm);

        return 0;
}

void g2e_index_xyz(unsigned int *idx, const CintEnvVars *envs)
{
        const unsigned int i_l = envs->i_l;
        const unsigned int j_l = envs->j_l;
        const unsigned int k_l = envs->k_l;
        const unsigned int l_l = envs->l_l;
        const unsigned int nfi = envs->nfi;
        const unsigned int nfj = envs->nfj;
        const unsigned int nfk = envs->nfk;
        const unsigned int nfl = envs->nfl;
        const unsigned int di = envs->g_stride_i;
        const unsigned int dk = envs->g_stride_k;
        const unsigned int dl = envs->g_stride_l;
        const unsigned int dj = envs->g_stride_j;
        unsigned int i, j, k, l, n;
        unsigned int ofx, ofkx, oflx;
        unsigned int ofy, ofky, ofly;
        unsigned int ofz, ofkz, oflz;
        unsigned int i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
        unsigned int j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];
        unsigned int k_nx[CART_MAX], k_ny[CART_MAX], k_nz[CART_MAX];
        unsigned int l_nx[CART_MAX], l_ny[CART_MAX], l_nz[CART_MAX];

        cart_comp(i_nx, i_ny, i_nz, i_l);
        cart_comp(j_nx, j_ny, j_nz, j_l);
        cart_comp(k_nx, k_ny, k_nz, k_l);
        cart_comp(l_nx, l_ny, l_nz, l_l);

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
static void g0_2e_2d(double *g, _BC *bc, const double *w, const double fac,
                     const CintEnvVars *envs)
{
        const unsigned int *ng = envs->ng;
        const unsigned int mmax = ng[2] - 1;
        const unsigned int nmax = ng[3] - 1;
        unsigned int i, m, n, off;
        const unsigned int nroots = ng[RYS_ROOTS];
        const unsigned int dm = envs->g_stride_l;
        const unsigned int dn = envs->g_stride_j;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *c00;
        const double *c0p;
        const double *b01 = bc->b01;
        const double *b00 = bc->b00;
        const double *b10 = bc->b10;

        for (i = 0; i < nroots; i++) {
                gx[i] = 1;
                gy[i] = 1;
                gz[i] = w[i] * fac;
        }

        if (nmax > 0) {
                // gx(irys,0,1) = c00(irys) * gx(irys,0,0)
                for (c00 = bc->c00, i = 0; i < nroots; i++, c00+=3) {
                        gx[i+dn] = c00[0] * gx[i];
                        gy[i+dn] = c00[1] * gy[i];
                        gz[i+dn] = c00[2] * gz[i];
                }
                // gx(irys,0,n+1) = c00(irys)*gx(irys,0,n)
                // + n*b10(irys)*gx(irys,0,n-1)
                for (n = 1; n < nmax; n++) {
                        for (c00 = bc->c00, off = n*dn, i = 0;
                             i < nroots; i++, c00+=3) {
                                gx[off+i+dn] = c00[0] * gx[off+i] + n * b10[i] * gx[off+i-dn];
                                gy[off+i+dn] = c00[1] * gy[off+i] + n * b10[i] * gy[off+i-dn];
                                gz[off+i+dn] = c00[2] * gz[off+i] + n * b10[i] * gz[off+i-dn];
                        }
                }
        }

        if (mmax > 0) {
                // gx(irys,1,0) = c0p(irys) * gx(irys,0,0)
                for (c0p = bc->c0p, i = 0; i < nroots; i++, c0p+=3) {
                        gx[i+dm] = c0p[0] * gx[i];
                        gy[i+dm] = c0p[1] * gy[i];
                        gz[i+dm] = c0p[2] * gz[i];
                }
                // gx(irys,m+1,0) = c0p(irys)*gx(irys,m,0)
                // + m*b01(irys)*gx(irys,m-1,0)
                for (m = 1; m < mmax; m++) {
                        for (c0p = bc->c0p, off = m*dm, i = 0;
                             i < nroots; i++, c0p+=3) {
                                gx[off+i+dm] = c0p[0] * gx[off+i] + m * b01[i] * gx[off+i-dm];
                                gy[off+i+dm] = c0p[1] * gy[off+i] + m * b01[i] * gy[off+i-dm];
                                gz[off+i+dm] = c0p[2] * gz[off+i] + m * b01[i] * gz[off+i-dm];
                        }
                }
        }

        if (nmax > 0 && mmax > 0) {
                // gx(irys,1,1) = c0p(irys)*gx(irys,0,1)
                // + b00(irys)*gx(irys,0,0)
                for (c0p = bc->c0p, i = 0; i < nroots; i++, c0p+=3) {
                        gx[i+dn+dm] = c0p[0] * gx[i+dn] + b00[i] * gx[i];
                        gy[i+dn+dm] = c0p[1] * gy[i+dn] + b00[i] * gy[i];
                        gz[i+dn+dm] = c0p[2] * gz[i+dn] + b00[i] * gz[i];
                }

                // gx(irys,m+1,1) = c0p(irys)*gx(irys,m,1)
                // + m*b01(irys)*gx(irys,m-1,1)
                // + b00(irys)*gx(irys,m,0)
                for (m = 1; m < mmax; m++) {
                        for (c0p = bc->c0p, off = m*dm+dn, i = 0;
                             i < nroots; i++, c0p+=3) {
                                gx[off+i+dm] = c0p[0]*gx[off+i] + m*b01[i]*gx[off+i-dm] + b00[i]*gx[off+i-dn];
                                gy[off+i+dm] = c0p[1]*gy[off+i] + m*b01[i]*gy[off+i-dm] + b00[i]*gy[off+i-dn];
                                gz[off+i+dm] = c0p[2]*gz[off+i] + m*b01[i]*gz[off+i-dm] + b00[i]*gz[off+i-dn];
                        }
                }

                // gx(irys,m,n+1) = c00(irys)*gx(irys,m,n)
                // + n*b10(irys)*gx(irys,m,n-1)
                // + m*b00(irys)*gx(irys,m-1,n)
                for (m = 1; m <= mmax; m++) {
                        for (n = 1; n < nmax; n++) {
                                for (c00 = bc->c00, off = m*dm+n*dn, i = 0;
                                     i < nroots; i++, c00+=3) {
                                        gx[off+i+dn] = c00[0]*gx[off+i] +n*b10[i]*gx[off+i-dn] + m*b00[i]*gx[off+i-dm];
                                        gy[off+i+dn] = c00[1]*gy[off+i] +n*b10[i]*gy[off+i-dn] + m*b00[i]*gy[off+i-dm];
                                        gz[off+i+dn] = c00[2]*gz[off+i] +n*b10[i]*gz[off+i-dn] + m*b00[i]*gz[off+i-dm];
                                }
                        }
                }
        }
}


/*
 * g0[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
void g0_2e_4d(double *g, const double *rirj, const double *rkrl,
              const CintEnvVars *envs)
{
        const unsigned int *ng = envs->ng;
        const unsigned int nmax = ng[3] - 1;
        const unsigned int mmax = ng[2] - 1;
        const unsigned int li = ng[0] - 1;
        const unsigned int lk = ng[1] - 1;
        //const unsigned int ll = mmax - lk;
        const unsigned int lj = nmax - li;
        unsigned int i, j, k, l, ptr, n;
        const unsigned int di = envs->g_stride_i;
        const unsigned int dk = envs->g_stride_k;
        const unsigned int dl = envs->g_stride_l;
        const unsigned int dj = envs->g_stride_j;
        DEF_GXYZ(double, g, gx, gy, gz);

        // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
        for (i = 1; i <= li; i++) {
                for (j = 0; j <= nmax-i; j++) {
                        for (l = 0; l <= mmax; l++) {
                                ptr = j*dj + l*dl + i*di;
                                for (n = ptr; n < ptr+di; n++) {
                                        gx[n] = rirj[0]*gx[n-di] + gx[n-di+dj];
                                        gy[n] = rirj[1]*gy[n-di] + gy[n-di+dj];
                                        gz[n] = rirj[2]*gz[n-di] + gz[n-di+dj];
                                }
                        }
                }
        }

        // g(0,k,l,..) = rkrl * g(0,k-1,l,..) + g(0,k-1,l+1,..)
        for (j = 0; j <= lj; j++) {
                for (k = 1; k <= lk; k++) {
                        for (l = 0; l < mmax-k; l+=2) {
                                ptr = j*dj + l*dl + k*dk;
                                for (n = ptr; n < ptr+dk; n++) {
                                        gx[n] = rkrl[0]*gx[n-dk] + gx[n-dk+dl];
                                        gy[n] = rkrl[1]*gy[n-dk] + gy[n-dk+dl];
                                        gz[n] = rkrl[2]*gz[n-dk] + gz[n-dk+dl];
                                        gx[n+dl] = rkrl[0]*gx[n+dl-dk] + gx[n+dl-dk+dl];
                                        gy[n+dl] = rkrl[1]*gy[n+dl-dk] + gy[n+dl-dk+dl];
                                        gz[n+dl] = rkrl[2]*gz[n+dl-dk] + gz[n+dl-dk+dl];
                                }
                        }
                        if (l <= mmax-k) {
                                ptr = j*dj + l*dl + k*dk;
                                for (n = ptr; n < ptr+dk; n++) {
                                        gx[n] = rkrl[0]*gx[n-dk] + gx[n-dk+dl];
                                        gy[n] = rkrl[1]*gy[n-dk] + gy[n-dk+dl];
                                        gz[n] = rkrl[2]*gz[n-dk] + gz[n-dk+dl];
                                }
                        }
                }
        }
}
/************* some special g0_4d results *************/
static inline void _g0_4d_ng_1112(double *g, double *c,
                                  const double *r, double *w, double fac1)
{
        g[0] = 1;
        g[1] = c[0];
        g[2] = 1;
        g[3] = c[1];
        g[4] = w[0] * fac1;
        g[5] = c[2] * g[4];
}
static inline void _g0_4d_ng_2112(double *g, double *c,
                                  const double *r, double *w, double fac1)
{
        g[0] = 1;
        g[1] = r[0] + c[0];
        g[4] = 1;
        g[5] = r[1] + c[1];
        g[8] = w[0] * fac1;
        g[9] =(r[2] + c[2]) * g[8];
}
static inline void _g0_4d_ng_1113(double *g, double *c, double *b,
                                  const double *r, double *w, double fac1)
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
        g[12] = w[0] * fac1;
        g[13] = w[1] * fac1;
        g[14] = c[2] * g[12];
        g[15] = c[5] * g[13];
        g[16] =(c[2] * c[2] + b[0])* g[12];
        g[17] =(c[5] * c[5] + b[1])* g[13];
}
static inline void _g0_4d_ng_2113(double *g, double *c, double *b,
                                  const double *r, double *w, double fac1)
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
        g[24] = w[0] * fac1;
        g[25] = w[1] * fac1;
        g[26] = rc[4] * g[24];
        g[27] = rc[5] * g[25];
        g[28] = c[2] * g[24];
        g[29] = c[5] * g[25];
        g[30] =(rc[4] * c[2] + b[0])* g[24];
        g[31] =(rc[5] * c[5] + b[1])* g[25];
}
static inline void _g0_4d_ng_3113(double *g, double *c, double *b,
                                  const double *r, double *w, double fac1)
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
        g[36] = w[0] * fac1;
        g[37] = w[1] * fac1;
        g[38] = rc[4] * g[36];
        g[39] = rc[5] * g[37];
        g[40] =(rc[4] * rc[4] + b[0])* g[36];
        g[41] =(rc[5] * rc[5] + b[1])* g[37];
}
static inline void _g0_4d_ng_1114(double *g, double *c, double *b,
                                  const double *r, double *w, double fac1)
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
        g[16] = w[0] * fac1;
        g[17] = w[1] * fac1;
        g[18] = c[2] * g[16];
        g[19] = c[5] * g[17];
        g[20] =(c[2] * c[2] + b[0])* g[16];
        g[21] =(c[5] * c[5] + b[1])* g[17];
        g[22] =(c[2] * c[2] + 3 * b[0])* c[2] * g[16];
        g[23] =(c[5] * c[5] + 3 * b[1])* c[5] * g[17];
}
static inline void _g0_4d_ng_2114(double *g, double *c, double *b,
                                  const double *r, double *w, double fac1)
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
        g[32] = w[0] * fac1;
        g[33] = w[1] * fac1;
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
static inline void _g0_4d_ng_3114(double *g, double *c, double *b,
                                  const double *r, double *w, double fac1)
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
        g[48] = w[0] * fac1;
        g[49] = w[1] * fac1;
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
static inline void _g0_4d_ng_4114(double *g, double *c, double *b,
                                  const double *r, double *w, double fac1)
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
        g[64] = w[0] * fac1;
        g[65] = w[1] * fac1;
        g[66] = rc[4] * g[64];
        g[67] = rc[5] * g[65];
        g[68] =(rc[4] * rc[4] + b[0])* g[64];
        g[69] =(rc[5] * rc[5] + b[1])* g[65];
        g[70] =(rc[4] * rc[4] + 3 * b[0])* rc[4] * g[64];
        g[71] =(rc[5] * rc[5] + 3 * b[1])* rc[5] * g[65];
}
static inline void _g0_4d_ng_1122(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[16] = w[0] * fac1;
        g[17] = w[1] * fac1;
        g[18] = cp[2] * g[16];
        g[19] = cp[5] * g[17];
        g[20] = c0[2] * g[16];
        g[21] = c0[5] * g[17];
        g[22] =(cp[2] * c0[2] + b[0]) * g[16];
        g[23] =(cp[5] * c0[5] + b[1]) * g[17];
}
static inline void _g0_4d_ng_2122(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[32] = w[0] * fac1;
        g[33] = w[1] * fac1;
        g[34] = rc[4] * g[32];
        g[35] = rc[5] * g[33];
        g[36] = cp[2] * g[32];
        g[37] = cp[5] * g[33];
        g[38] =(rc[4]*cp[2] + b[0]) * g[32];
        g[39] =(rc[5]*cp[5] + b[1]) * g[33];
}
static inline void _g0_4d_ng_1222(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[32] = w[0] * fac1;
        g[33] = w[1] * fac1;
        g[34] = rc[4] * g[32];
        g[35] = rc[5] * g[33];
        g[40] = c0[2] * g[32];
        g[41] = c0[5] * g[33];
        g[42] =(rc[4]*c0[2] + b[0]) * g[32];
        g[43] =(rc[5]*c0[5] + b[1]) * g[33];
}
static inline void _g0_4d_ng_2222(double *g, double *c0, double *cp, double *b,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[64] = w[0] * fac1;
        g[65] = w[1] * fac1;
        g[66] = rc0[4] * g[64];
        g[67] = rc0[5] * g[65];
        g[68] = rcp[4] * g[64];
        g[69] = rcp[5] * g[65];
        g[70] =(rc0[4]*rcp[4] + b[0]) * g[64];
        g[71] =(rc0[5]*rcp[5] + b[1]) * g[65];
}
static inline void _g0_4d_ng_1132(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  double *w, double fac1)
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
        g[24] = w[0] * fac1;
        g[25] = w[1] * fac1;
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
static inline void _g0_4d_ng_2132(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[48] = w[0] * fac1;
        g[49] = w[1] * fac1;
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
static inline void _g0_4d_ng_1232(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[48] = w[0] * fac1;
        g[49] = w[1] * fac1;
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
static inline void _g0_4d_ng_2232(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[96 ] = w[0] * fac1;
        g[97 ] = w[1] * fac1;
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
static inline void _g0_4d_ng_1332(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[72] = w[0] * fac1;
        g[73] = w[1] * fac1;
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
static inline void _g0_4d_ng_2332(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[144] = w[0] * fac1;
        g[145] = w[1] * fac1;
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
static inline void _g0_4d_ng_1123(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  double *w, double fac1)
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
        g[24] = w[0] * fac1;
        g[25] = w[1] * fac1;
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
static inline void _g0_4d_ng_2123(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[48] = w[0] * fac1;
        g[49] = w[1] * fac1;
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
static inline void _g0_4d_ng_3123(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[72] = w[0] * fac1;
        g[73] = w[1] * fac1;
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
static inline void _g0_4d_ng_1223(double *g, double *c0, double *cp,
                                  double *b0, double *b1,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[48] = w[0] * fac1;
        g[49] = w[1] * fac1;
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
static inline void _g0_4d_ng_2223(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[96 ] = w[0] * fac1;
        g[97 ] = w[1] * fac1;
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
static inline void _g0_4d_ng_3223(double *g, double *c0, double *cp, 
                                  double *b0, double *b1,
                                  const double *r0, const double *rp,
                                  double *w, double fac1)
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
        g[144] = w[0] * fac1;
        g[145] = w[1] * fac1;
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


/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
void g0_2e(double *g, const RijSets *rij, const RijSets *rkl,
           const double fac, const CintEnvVars *envs)
{
        const unsigned int *ng = envs->ng;
        const double aij = rij->a12;
        const double akl = rkl->a12;
        double a0, a1, fac1, x;
        double u[MXRYSROOTS];
        double w[MXRYSROOTS];
        double rijrkl[3];
        rijrkl[0] = rij->r12[0] - rkl->r12[0];
        rijrkl[1] = rij->r12[1] - rkl->r12[1];
        rijrkl[2] = rij->r12[2] - rkl->r12[2];

        a1 = aij * akl;
        a0 = a1 / (aij + akl);
        fac1 = sqrt(a0 / (a1 * a1 * a1)) * fac;
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);
        rys_roots(ng[RYS_ROOTS], x, u, w);

        unsigned int irys;
        double u2, div;
        const double *rijrj = rij->r12r2;
        const double *rklrl = rkl->r12r2;
        const double *rirj = rij->r1r2;
        const double *rkrl = rkl->r1r2;
        _BC bc;
        double *c00 = bc.c00;
        double *c0p = bc.c0p;

        for (irys = 0; irys < ng[RYS_ROOTS]; irys++, c00+=3, c0p+=3)
        {
                /*
                 *t2 = u(irys)/(1+u(irys))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = a0 * u[irys];
                div = 1 / (u2 * (aij + akl) + aij * akl);
                bc.b00[irys] = 0.5 * u2 * div;
                bc.b10[irys] = 0.5 * (u2 + akl) * div;
                bc.b01[irys] = 0.5 * (u2 + aij) * div;
                c00[0] = rijrj[0] - u2 * akl * rijrkl[0] * div;
                c00[1] = rijrj[1] - u2 * akl * rijrkl[1] * div;
                c00[2] = rijrj[2] - u2 * akl * rijrkl[2] * div;
                c0p[0] = rklrl[0] + u2 * aij * rijrkl[0] * div;
                c0p[1] = rklrl[1] + u2 * aij * rijrkl[1] * div;
                c0p[2] = rklrl[2] + u2 * aij * rijrkl[2] * div;
        }

        switch (ng[3]) {
                case 1: switch(ng[2]) {
                        case 1: goto _g0_4d_default; // ssss
                        case 2: switch (ng[1]) {
                                case 1: _g0_4d_ng_1112(g, bc.c0p, rkrl, w,fac1); goto normal_end;
                                case 2: _g0_4d_ng_2112(g, bc.c0p, rkrl, w,fac1); goto normal_end;
                                default: goto error; }
                        case 3: switch (ng[1]) {
                                case 1: _g0_4d_ng_1113(g, bc.c0p, bc.b01, rkrl, w, fac1); goto normal_end;
                                case 2: _g0_4d_ng_2113(g, bc.c0p, bc.b01, rkrl, w, fac1); goto normal_end;
                                case 3: _g0_4d_ng_3113(g, bc.c0p, bc.b01, rkrl, w, fac1); goto normal_end;
                                default: goto error; }
                        case 4: switch (ng[1]) {
                                case 1: _g0_4d_ng_1114(g, bc.c0p, bc.b01, rkrl, w, fac1); goto normal_end;
                                case 2: _g0_4d_ng_2114(g, bc.c0p, bc.b01, rkrl, w, fac1); goto normal_end;
                                case 3: _g0_4d_ng_3114(g, bc.c0p, bc.b01, rkrl, w, fac1); goto normal_end;
                                case 4: _g0_4d_ng_4114(g, bc.c0p, bc.b01, rkrl, w, fac1); goto normal_end;
                                default: goto error; }
                        default: goto _g0_4d_default; }
                case 2: switch(ng[2]) {
                        case 1: switch (ng[0]) {
                                case 1: _g0_4d_ng_1112(g, bc.c00, rirj, w,fac1); goto normal_end;
                                case 2: _g0_4d_ng_2112(g, bc.c00, rirj, w,fac1); goto normal_end;
                                default: goto error; }
                        case 2: switch (ng[1]) {
                                case 1: switch (ng[0]) {
                                        case 1: _g0_4d_ng_1122(g, bc.c00, bc.c0p, bc.b00, rirj, rkrl, w, fac1); goto normal_end;
                                        case 2: _g0_4d_ng_2122(g, bc.c00, bc.c0p, bc.b00, rirj, rkrl, w, fac1); goto normal_end;
                                        default: goto error; }
                                case 2: switch (ng[0]) {
                                        case 1: _g0_4d_ng_1222(g, bc.c00, bc.c0p, bc.b00, rirj, rkrl, w, fac1); goto normal_end;
                                        case 2: _g0_4d_ng_2222(g, bc.c00, bc.c0p, bc.b00, rirj, rkrl, w, fac1); goto normal_end;
                                        default: goto error; }
                                default: goto error; }
                        case 3: switch (ng[1]) {
                                case 1: switch (ng[0]) {
                                        case 1: _g0_4d_ng_1132(g, bc.c00, bc.c0p, bc.b00, bc.b01, w, fac1); goto normal_end;
                                        case 2: _g0_4d_ng_2132(g, bc.c00, bc.c0p, bc.b00, bc.b01, rirj, rkrl, w, fac1); goto normal_end;
                                        default: goto error; }
                                case 2: switch (ng[0]) {
                                        case 1: _g0_4d_ng_1232(g, bc.c00, bc.c0p, bc.b00, bc.b01, rirj, rkrl, w, fac1); goto normal_end;
                                        case 2: _g0_4d_ng_2232(g, bc.c00, bc.c0p, bc.b00, bc.b01, rirj, rkrl, w, fac1); goto normal_end;
                                        default: goto error; }
                                case 3: switch (ng[0]) {
                                        case 1: _g0_4d_ng_1332(g, bc.c00, bc.c0p, bc.b00, bc.b01, rirj, rkrl, w, fac1); goto normal_end;
                                        case 2: _g0_4d_ng_2332(g, bc.c00, bc.c0p, bc.b00, bc.b01, rirj, rkrl, w, fac1); goto normal_end;
                                        default: goto error; }
                                default: goto error; }
                        default: goto _g0_4d_default; }
                case 3: switch(ng[2]) {
                        case 1: switch (ng[0]) {
                                case 1: _g0_4d_ng_1113(g, bc.c00, bc.b10, rirj, w, fac1); goto normal_end;
                                case 2: _g0_4d_ng_2113(g, bc.c00, bc.b10, rirj, w, fac1); goto normal_end;
                                case 3: _g0_4d_ng_3113(g, bc.c00, bc.b10, rirj, w, fac1); goto normal_end;
                                default: goto error; }
                        case 2: switch (ng[1]) {
                                case 1: switch (ng[0]) {
                                        case 1: _g0_4d_ng_1123(g, bc.c00, bc.c0p, bc.b00, bc.b10, w, fac1); goto normal_end;
                                        case 2: _g0_4d_ng_2123(g, bc.c00, bc.c0p, bc.b00, bc.b10, rirj, rkrl, w, fac1); goto normal_end;
                                        case 3: _g0_4d_ng_3123(g, bc.c00, bc.c0p, bc.b00, bc.b10, rirj, rkrl, w, fac1); goto normal_end;
                                        default: goto error; }
                                case 2: switch (ng[0]) {
                                        case 1: _g0_4d_ng_1223(g, bc.c00, bc.c0p, bc.b00, bc.b10, rirj, rkrl, w, fac1); goto normal_end;
                                        case 2: _g0_4d_ng_2223(g, bc.c00, bc.c0p, bc.b00, bc.b10, rirj, rkrl, w, fac1); goto normal_end;
                                        case 3: _g0_4d_ng_3223(g, bc.c00, bc.c0p, bc.b00, bc.b10, rirj, rkrl, w, fac1); goto normal_end;
                                        default: goto error; }
                                default: goto error; }
                        default: goto _g0_4d_default; }
                case 4: switch(ng[2]) {
                        case 1: switch (ng[0]) {
                                case 1: _g0_4d_ng_1114(g, bc.c00, bc.b10, rirj, w, fac1); goto normal_end;
                                case 2: _g0_4d_ng_2114(g, bc.c00, bc.b10, rirj, w, fac1); goto normal_end;
                                case 3: _g0_4d_ng_3114(g, bc.c00, bc.b10, rirj, w, fac1); goto normal_end;
                                case 4: _g0_4d_ng_4114(g, bc.c00, bc.b10, rirj, w, fac1); goto normal_end;
                                default: goto error; }
                        default: goto _g0_4d_default; }
                default:
_g0_4d_default:
                        g0_2e_2d(g, &bc, w, fac1, envs);
                        g0_2e_4d(g, rij->r1r2, rkl->r1r2, envs);
        }
normal_end:
        return;
error:
        printf("Dimension error: %u %u %u %u\n", ng[0], ng[1], ng[2], ng[3]);
        exit(1);
}


static double rys_root1(double x)
{
        const double pie4 = 7.85398163397448e-01;
        double ww1;
        double f1,e,y;
        if (x < 3.e-7){
                ww1 = 1.0e+00 -x/3.0e+00;
        } else if (x < 1.) {
                f1 = ((((((((-8.36313918003957e-08*x+1.21222603512827e-06 )*x-
                            1.15662609053481e-05 )*x+9.25197374512647e-05 )*x-
                          6.40994113129432e-04 )*x+3.78787044215009e-03 )*x-
                        1.85185172458485e-02 )*x+7.14285713298222e-02 )*x-
                      1.99999999997023e-01 )*x+3.33333333333318e-01;
                ww1 = (x+x)*f1+exp(-x);
        } else if (x < 3.) {
                y = x-2.0e+00;
                f1 = ((((((((((-1.61702782425558e-10*y+1.96215250865776e-09 )*y-
                              2.14234468198419e-08 )*y+2.17216556336318e-07 )*y-
                            1.98850171329371e-06 )*y+1.62429321438911e-05 )*y-
                          1.16740298039895e-04 )*y+7.24888732052332e-04 )*y-
                        3.79490003707156e-03 )*y+1.61723488664661e-02 )*y-
                      5.29428148329736e-02 )*y+1.15702180856167e-01;
                ww1 = (x+x)*f1+exp(-x);
        } else if (x < 5.) {
                y = x-4.0e+00;
                f1 = ((((((((((-2.62453564772299e-11*y+3.24031041623823e-10 )*y-
                              3.614965656163e-09)*y+3.760256799971e-08)*y-
                            3.553558319675e-07)*y+3.022556449731e-06)*y-
                          2.290098979647e-05)*y+1.526537461148e-04)*y-
                        8.81947375894379e-04 )*y+4.33207949514611e-03 )*y-
                      1.75257821619926e-02 )*y+5.28406320615584e-02;
                ww1 = (x+x)*f1+exp(-x);
        } else if (x < 10) {
                e = exp(-x);
                ww1 = (((((( 4.6897511375022e-01/x-6.9955602298985e-01)/x +
                           5.3689283271887e-01)/x-3.2883030418398e-01)/x +
                         2.4645596956002e-01)/x-4.9984072848436e-01)/x -
                       3.1501078774085e-06)*e + sqrt(pie4/x);
        } else if (x < 15) {
                e = exp(-x);
                ww1 = (((-1.8784686463512e-01/x+2.2991849164985e-01)/x -
                        4.9893752514047e-01)/x-2.1916512131607e-05)*e
                        + sqrt(pie4/x);
        } else if (x < 33) {
                e = exp(-x);
                ww1 = (( 1.9623264149430e-01/x-4.9695241464490e-01)/x -
                       6.0156581186481e-05)*e + sqrt(pie4/x);
        } else {
                ww1 = sqrt(pie4/x);
        }
        return ww1;
}
double g0_2e_ssss(const RijSets *rij, const RijSets *rkl,
                  const double fac, const CintEnvVars *envs)
{
        const double aij = rij->a12;
        const double akl = rkl->a12;
        double a0, a1, x;
        double rijrkl[3];
        rijrkl[0] = rij->r12[0] - rkl->r12[0];
        rijrkl[1] = rij->r12[1] - rkl->r12[1];
        rijrkl[2] = rij->r12[2] - rkl->r12[2];

        a1 = aij * akl;
        a0 = a1 / (aij + akl);
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);

        return sqrt(a0 / (a1 * a1 * a1)) * fac * rys_root1(x);
}

/*
 * ( \nabla i j | kl )
 */
void nabla1i_2e(double *f, const double *g,
                const unsigned int li, const unsigned int lj,
                const unsigned int lk, const unsigned int ll,
                const CintEnvVars *envs)
{
        unsigned int i, j, k, l, n, ptr;
        const unsigned int di = envs->g_stride_i;
        const unsigned int dk = envs->g_stride_k;
        const unsigned int dl = envs->g_stride_l;
        const unsigned int dj = envs->g_stride_j;
        const double ai = envs->ai;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++)
                        for (k = 0; k <= lk; k++) {
                                ptr = dj * j + dl * l + dk * k;
                                //f(...,0,...) = -2*ai*g(...,1,...)
                                for (n = ptr; n < ptr + di; n++) {
                                        fx[n] = - 2 * ai * gx[n+di];
                                        fy[n] = - 2 * ai * gy[n+di];
                                        fz[n] = - 2 * ai * gz[n+di];
                                }
                                ptr += di;
                                //f(...,i,...) = i*g(...,i-1,...)-2*ai*g(...,i+1,...)
                                for (i = 1; i <= li; i++) {
                                        for (n = ptr; n < ptr + di; n++) {
                                                fx[n] = i*gx[n-di] - 2*ai*gx[n+di];
                                                fy[n] = i*gy[n-di] - 2*ai*gy[n+di];
                                                fz[n] = i*gz[n-di] - 2*ai*gz[n+di];
                                        }
                                        ptr += di;
                                }
                        }
}


/*
 * ( i \nabla j | kl )
 */
void nabla1j_2e(double *f, const double *g,
                const unsigned int li, const unsigned int lj,
                const unsigned int lk, const unsigned int ll,
                const CintEnvVars *envs)
{
        unsigned int i, j, k, l, ptr;
        const unsigned int di = envs->g_stride_i;
        const unsigned int dk = envs->g_stride_k;
        const unsigned int dl = envs->g_stride_l;
        const unsigned int dj = envs->g_stride_j;
        const double aj = envs->aj;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        //f(...,0,...) = -2*aj*g(...,1,...)
        for (l = 0; l <= ll; l++) {
                ptr = dl * l;
                for (k = 0; k <= lk; k++) {
                        for (i = ptr; i < ptr + di * (li + 1); i++) {
                                fx[i] = - 2 * aj * gx[i+dj];
                                fy[i] = - 2 * aj * gy[i+dj];
                                fz[i] = - 2 * aj * gz[i+dj];
                        }
                        ptr += dk;
                }
        }
        //f(...,j,...) = j*g(...,j-1,...)-2*aj*g(...,j+1,...)
        for (j = 1; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                        ptr = dj * j + dl * l;
                        for (k = 0; k <= lk; k++) {
                                for (i = ptr; i < ptr + di * (li + 1); i++) {
                                        fx[i] = j*gx[i-dj] - 2*aj*gx[i+dj];
                                        fy[i] = j*gy[i-dj] - 2*aj*gy[i+dj];
                                        fz[i] = j*gz[i-dj] - 2*aj*gz[i+dj];
                                }
                                ptr += dk;
                        }
                }
}


/*
 * ( ij | \nabla k l )
 */
void nabla1k_2e(double *f, const double *g,
                const unsigned int li, const unsigned int lj,
                const unsigned int lk, const unsigned int ll,
                const CintEnvVars *envs)
{
        unsigned int i, j, k, l, ptr;
        const unsigned int di = envs->g_stride_i;
        const unsigned int dk = envs->g_stride_k;
        const unsigned int dl = envs->g_stride_l;
        const unsigned int dj = envs->g_stride_j;
        const double ak = envs->ak;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                        ptr = dj * j + dl * l;
                        //f(...,0,...) = -2*ak*g(...,1,...)
                        for (i = ptr; i < ptr + di * (li + 1); i++) {
                                fx[i] = - 2 * ak * gx[i+dk];
                                fy[i] = - 2 * ak * gy[i+dk];
                                fz[i] = - 2 * ak * gz[i+dk];
                        }
                        ptr += dk;
                        //f(...,k,...) = k*g(...,k-1,...)-2*ak*g(...,k+1,...)
                        for (k = 1; k <= lk; k++) {
                                for (i = ptr; i < ptr + di * (li + 1); i++) {
                                        fx[i] = k*gx[i-dk] - 2*ak*gx[i+dk];
                                        fy[i] = k*gy[i-dk] - 2*ak*gy[i+dk];
                                        fz[i] = k*gz[i-dk] - 2*ak*gz[i+dk];
                                }
                                ptr += dk;
                        }
                }
}


/*
 * ( ij | k \nabla l )
 */
void nabla1l_2e(double *f, const double *g,
                const unsigned int li, const unsigned int lj,
                const unsigned int lk, const unsigned int ll,
                const CintEnvVars *envs)
{
        unsigned int i, j, k, l, ptr;
        const unsigned int di = envs->g_stride_i;
        const unsigned int dk = envs->g_stride_k;
        const unsigned int dl = envs->g_stride_l;
        const unsigned int dj = envs->g_stride_j;
        const double al = envs->al;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        for (j = 0; j <= lj; j++) {
                ptr = dj * j;
                //f(...,0,...) = -2*al*g(...,1,...)
                for (k = 0; k <= lk; k++) {
                        for (i = ptr; i < ptr + di * (li + 1); i++) {
                                fx[i] = - 2 * al * gx[i+dl];
                                fy[i] = - 2 * al * gy[i+dl];
                                fz[i] = - 2 * al * gz[i+dl];
                        }
                        ptr += dk;
                }
                //f(...,l,...) = l*g(...,l-1,...)-2*al*g(...,l+1,...)
                for (l = 1; l <= ll; l++) {
                        ptr = dj * j + dl * l;
                        for (k = 0; k <= lk; k++) {
                                for (i = ptr; i < ptr + di * (li + 1); i++) {
                                        fx[i] = l*gx[i-dl] - 2*al*gx[i+dl];
                                        fy[i] = l*gy[i-dl] - 2*al*gy[i+dl];
                                        fz[i] = l*gz[i-dl] - 2*al*gz[i+dl];
                                }
                                ptr += dk;
                        }
                }
        }
}

/*
 * ( x^1 i j | kl )
 * ri is the shift from the center R_O to the center of |i>
 * r - R_O = (r-R_i) + ri, ri = R_i - R_O
 */
void x1i_2e(double *f, const double *g,
            const unsigned int li, const unsigned int lj,
            const unsigned int lk, const unsigned int ll,
            const double *ri, const CintEnvVars *envs)
{
        unsigned int i, j, k, l, ptr;
        const unsigned int di = envs->g_stride_i;
        const unsigned int dk = envs->g_stride_k;
        const unsigned int dl = envs->g_stride_l;
        const unsigned int dj = envs->g_stride_j;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                        ptr = dj * j + dl * l;
                        for (k = 0; k <= lk; k++) {
                                //f(...,0:li,...) = g(...,1:li+1,...) + ri(1)*g(...,0:li,...)
                                for (i = ptr; i < ptr + di * (li + 1); i++) {
                                        fx[i] = gx[i+di] + ri[0] * gx[i];
                                        fy[i] = gy[i+di] + ri[1] * gy[i];
                                        fz[i] = gz[i+di] + ri[2] * gz[i];
                                }
                                ptr += dk;
                        }
                }
}


/*
 * ( i x^1 j | kl )
 */
void x1j_2e(double *f, const double *g,
            const unsigned int li, const unsigned int lj,
            const unsigned int lk, const unsigned int ll,
            const double *rj, const CintEnvVars *envs)
{
        unsigned int i, j, k, l, ptr;
        const unsigned int di = envs->g_stride_i;
        const unsigned int dk = envs->g_stride_k;
        const unsigned int dl = envs->g_stride_l;
        const unsigned int dj = envs->g_stride_j;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                        ptr = dj * j + dl * l;
                        for (k = 0; k <= lk; k++) {
                                // f(...,0:lj,...) = g(...,1:lj+1,...)
                                // + rj(1)*g(...,0:lj,...)
                                for (i = ptr; i < ptr + di * (li + 1); i++) {
                                        fx[i] = gx[i+dj] + rj[0] * gx[i];
                                        fy[i] = gy[i+dj] + rj[1] * gy[i];
                                        fz[i] = gz[i+dj] + rj[2] * gz[i];
                                }
                                ptr += dk;
                        }
                }
}


/*
 * ( ij | x^1 k l )
 */
void x1k_2e(double *f, const double *g,
            const unsigned int li, const unsigned int lj,
            const unsigned int lk, const unsigned int ll,
            const double *rk, const CintEnvVars *envs)
{
        unsigned int i, j, k, l, ptr;
        const unsigned int di = envs->g_stride_i;
        const unsigned int dk = envs->g_stride_k;
        const unsigned int dl = envs->g_stride_l;
        const unsigned int dj = envs->g_stride_j;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                        ptr = dj * j + dl * l;
                        for (k = 0; k <= lk; k++) {
                                // f(...,0:lk,...) = g(...,1:lk+1,...)
                                // + rk(1)*g(...,0:lk,...)
                                for (i = ptr; i < ptr + di * (li + 1); i++) {
                                        fx[i] = gx[i+dk] + rk[0] * gx[i];
                                        fy[i] = gy[i+dk] + rk[1] * gy[i];
                                        fz[i] = gz[i+dk] + rk[2] * gz[i];
                                }
                                ptr += dk;
                        }
                }
}


/*
 * ( i j | x^1 kl )
 */
void x1l_2e(double *f, const double *g,
            const unsigned int li, const unsigned int lj,
            const unsigned int lk, const unsigned int ll,
            const double *rl, const CintEnvVars *envs)
{
        unsigned int i, j, k, l, ptr;
        const unsigned int di = envs->g_stride_i;
        const unsigned int dk = envs->g_stride_k;
        const unsigned int dl = envs->g_stride_l;
        const unsigned int dj = envs->g_stride_j;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        for (j = 0; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                        ptr = dj * j + dl * l;
                        for (k = 0; k <= lk; k++) {
                                // f(...,0:ll,...) = g(...,1:ll+1,...)
                                // + rl(1)*g(...,0:ll,...)
                                for (i = ptr; i < ptr + di * (li + 1); i++) {
                                        fx[i] = gx[i+dl] + rl[0] * gx[i];
                                        fy[i] = gy[i+dl] + rl[1] * gy[i];
                                        fz[i] = gz[i+dl] + rl[2] * gz[i];
                                }
                                ptr += dk;
                        }
                }
}

