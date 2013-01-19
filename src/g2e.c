/*
 * File: g2e.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <math.h>
#include "cint_bas.h"
#include "misc.h"
#include "g2e.h"


inline int g_size(const int *ng)
{
        return ng[RYS_ROOTS] * ng[0] * ng[1] * ng[2] * ng[3];
}

inline void extract_dim(const int *ng, int *di, int *dj, int *dk, int *dl)
{
        *di = ng[RYS_ROOTS]; // shift of (i++,k,l,j)
        *dk = ng[RYS_ROOTS] * ng[0]; // shift of (i,k++,l,j)
        *dl = ng[RYS_ROOTS] * ng[0] * ng[1]; // shift of (i,k,l++,j)
        *dj = ng[RYS_ROOTS] * ng[0] * ng[1] * ng[2]; // shift of (i,k,l,j++)
}


inline void extract_fg_xyz(double *f, const double *g, const int *ng,
                           const double **pgx, const double **pgy, const double **pgz, 
                           double **pfx, double **pfy, double **pfz)
{
        const int len = g_size(ng);
        *pgx = g;
        *pgy = g + len;
        *pgz = g + len + len;
        *pfx = f;
        *pfy = f + len;
        *pfz = f + len + len;
}


void g2e_index_xyz(int *idx, const int *ng, const int *shls,
                   const int *bas)
{
        const int i_l = bas(ANG_OF, shls[0]);
        const int j_l = bas(ANG_OF, shls[1]);
        const int k_l = bas(ANG_OF, shls[2]);
        const int l_l = bas(ANG_OF, shls[3]);
        const int nfi = len_cart(i_l);
        const int nfj = len_cart(j_l);
        const int nfk = len_cart(k_l);
        const int nfl = len_cart(l_l);
        int i, j, k, l, n;
        int ofx, ofkx, oflx;
        int ofy, ofky, ofly;
        int ofz, ofkz, oflz;
        int di, dk, dl, dj;
        int i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
        int j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];
        int k_nx[CART_MAX], k_ny[CART_MAX], k_nz[CART_MAX];
        int l_nx[CART_MAX], l_ny[CART_MAX], l_nz[CART_MAX];
        int *idy = idx + nfi * nfj * nfk * nfl;
        int *idz = idx + nfi * nfj * nfk * nfl * 2;

        cart_comp(i_nx, i_ny, i_nz, i_l);
        cart_comp(j_nx, j_ny, j_nz, j_l);
        cart_comp(k_nx, k_ny, k_nz, k_l);
        cart_comp(l_nx, l_ny, l_nz, l_l);

        ofx = 0;
        ofy = g_size(ng);
        ofz = g_size(ng) * 2;
        extract_dim(ng, &di, &dj, &dk, &dl);
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
                                for (i = 0; i < nfi; i++) {
                                        idx[n] = ofkx + di * i_nx[i]; //(:,ix,kx,lx,jx,1)
                                        idy[n] = ofky + di * i_ny[i]; //(:,iy,ky,ly,jy,2)
                                        idz[n] = ofkz + di * i_nz[i]; //(:,iz,kz,lz,jz,3)
                                        n++;
                                } // i
                        } // k
                } // l
        } // j

}


void g0_2e_bc(const int nroots, const double aij, const double akl,
              const double *rij, const double *rkl,
              const double *ri, const double *rk, const double *u, 
              double *c00, double *c0p, double *b01, double *b00, double *b10)
{
        int irys;
        double u2, div;
        double rijri[3], rijrk[3], rijrkl[3];
        double *c00x = c00;
        double *c00y = c00 + MXRYSROOTS;
        double *c00z = c00 + MXRYSROOTS + MXRYSROOTS;
        double *c0px = c0p;
        double *c0py = c0p + MXRYSROOTS;
        double *c0pz = c0p + MXRYSROOTS + MXRYSROOTS;

        rijri[0] = rij[0] - ri[0];
        rijri[1] = rij[1] - ri[1];
        rijri[2] = rij[2] - ri[2];
        rijrk[0] = rkl[0] - rk[0];
        rijrk[1] = rkl[1] - rk[1];
        rijrk[2] = rkl[2] - rk[2];
        rijrkl[0] = rij[0] - rkl[0];
        rijrkl[1] = rij[1] - rkl[1];
        rijrkl[2] = rij[2] - rkl[2];
        for (irys = 0; irys < nroots; irys++)
        {
                /*
                 *t2 = u(irys)/(1+u(irys))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = aij * akl / (aij + akl) * u[irys];
                div = 1 / (u2 * (aij + akl) + aij * akl);
                b00[irys] = 0.5 * u2 * div;
                b10[irys] = 0.5 * (u2 + akl) * div;
                b01[irys] = 0.5 * (u2 + aij) * div;
                c00x[irys] = rijri[0] - u2 * akl * rijrkl[0] * div;
                c00y[irys] = rijri[1] - u2 * akl * rijrkl[1] * div;
                c00z[irys] = rijri[2] - u2 * akl * rijrkl[2] * div;
                c0px[irys] = rijrk[0] + u2 * aij * rijrkl[0] * div;
                c0py[irys] = rijrk[1] + u2 * aij * rijrkl[1] * div;
                c0pz[irys] = rijrk[2] + u2 * aij * rijrkl[2] * div;
        }
}


/*
 * g(nroots,0:nmax,0:mmax)
 */
void g0_2e_2d(double *g, const int *ng, const double *c00, const double *c0p, 
              const double *b01, const double *b00, const double *b10,
              const double *w, const double fac)
{
        const int nmax = ng[0] - 1;
        const int mmax = ng[1] - 1;
        int i, m, n;
        const int dn = ng[RYS_ROOTS]; // shift of (i++,k,l,j)
        const int dm = ng[RYS_ROOTS] * ng[0]; // shift of (i,k++,l,j)
        double *gx = g;
        double *gy = g + ng[RYS_ROOTS] * ng[0] * ng[1] * ng[2] * ng[3];
        double *gz = g + ng[RYS_ROOTS] * ng[0] * ng[1] * ng[2] * ng[3] * 2;
        double *pgx, *pgy, *pgz;
        const double *c00x = c00;
        const double *c00y = c00 + MXRYSROOTS;
        const double *c00z = c00 + MXRYSROOTS + MXRYSROOTS;
        const double *c0px = c0p;
        const double *c0py = c0p + MXRYSROOTS;
        const double *c0pz = c0p + MXRYSROOTS + MXRYSROOTS;

        for (i = 0; i < dn; i++) {
                gx[i] = 1;
                gy[i] = 1;
                gz[i] = w[i] * fac;
        }

        if (nmax > 0) {
                // gx(irys,1,0) = c00(irys) * gx(irys,0,0)
                for (i = 0; i < dn; i++) {
                        gx[i+dn] = c00x[i] * gx[i];
                        gy[i+dn] = c00y[i] * gy[i];
                        gz[i+dn] = c00z[i] * gz[i];
                }
                pgx = gx + dn;
                pgy = gy + dn;
                pgz = gz + dn;
                // gx(irys,n+1,0) = c00(irys)*gx(irys,n,0)
                // + n*b10(irys)*gx(irys,n-1,0)
                for (n = 1; n < nmax; n++) {
                        for (i = 0; i < dn; i++) {
                                pgx[i+dn] = c00x[i] * pgx[i] + n * b10[i] * pgx[i-dn];
                                pgy[i+dn] = c00y[i] * pgy[i] + n * b10[i] * pgy[i-dn];
                                pgz[i+dn] = c00z[i] * pgz[i] + n * b10[i] * pgz[i-dn];
                        }
                        pgx += dn;
                        pgy += dn;
                        pgz += dn;
                }
        }

        if (mmax > 0) {
                // gx(irys,0,1) = c0p(irys) * gx(irys,0,0)
                for (i = 0; i < dn; i++) {
                        gx[i+dm] = c0px[i] * gx[i];
                        gy[i+dm] = c0py[i] * gy[i];
                        gz[i+dm] = c0pz[i] * gz[i];
                }
                pgx = gx + dm;
                pgy = gy + dm;
                pgz = gz + dm;
                // gx(irys,0,m+1) = c0p(irys)*gx(irys,0,m)
                // + m*b01(irys)*gx(irys,0,m-1)
                for (m = 1; m < mmax; m++) {
                        for (i = 0; i < dn; i++) {
                                pgx[i+dm] = c0px[i] * pgx[i] + m * b01[i] * pgx[i-dm];
                                pgy[i+dm] = c0py[i] * pgy[i] + m * b01[i] * pgy[i-dm];
                                pgz[i+dm] = c0pz[i] * pgz[i] + m * b01[i] * pgz[i-dm];
                        }
                        pgx += dm;
                        pgy += dm;
                        pgz += dm;
                }
        }

        if (nmax > 0 && mmax > 0) {
                // gx(irys,1,1) = c0p(irys)*gx(irys,1,0)
                // + b00(irys)*gx(irys,0,0)
                pgx = gx + dn;
                pgy = gy + dn;
                pgz = gz + dn;
                for (i = 0; i < dn; i++) {
                        pgx[i+dm] = c0px[i] * pgx[i] + b00[i] * pgx[i-dn];
                        pgy[i+dm] = c0py[i] * pgy[i] + b00[i] * pgy[i-dn];
                        pgz[i+dm] = c0pz[i] * pgz[i] + b00[i] * pgz[i-dn];
                }

                pgx = gx + dm + dn;
                pgy = gy + dm + dn;
                pgz = gz + dm + dn;
                // gx(irys,1,m+1) = c0p(irys)*gx(irys,1,m)
                // + m*b01(irys)*gx(irys,1,m-1)
                // + b00(irys)*gx(irys,0,m)
                for (m = 1; m < mmax; m++) {
                        for (i = 0; i < dn; i++) {
                                pgx[i+dm] = c0px[i] * pgx[i] + m * b01[i] * pgx[i-dm] + b00[i] * pgx[i-dn];
                                pgy[i+dm] = c0py[i] * pgy[i] + m * b01[i] * pgy[i-dm] + b00[i] * pgy[i-dn];
                                pgz[i+dm] = c0pz[i] * pgz[i] + m * b01[i] * pgz[i-dm] + b00[i] * pgz[i-dn];
                        }
                        pgx += dm;
                        pgy += dm;
                        pgz += dm;
                }

                // gx(irys,n+1,m) = c00(irys)*gx(irys,n,m)
                // + n*b10(irys)*gx(irys,n-1,m)
                // + m*b00(irys)*gx(irys,n,m-1)
                for (m = 1; m <= mmax; m++) {
                        pgx = gx + dm * m + dn;
                        pgy = gy + dm * m + dn;
                        pgz = gz + dm * m + dn;
                        for (n = 1; n < nmax; n++) {
                                for (i = 0; i < dn; i++) {
                                        pgx[i+dn] = c00x[i] * pgx[i] + n * b10[i] * pgx[i-dn] + m * b00[i] * pgx[i-dm];
                                        pgy[i+dn] = c00y[i] * pgy[i] + n * b10[i] * pgy[i-dn] + m * b00[i] * pgy[i-dm];
                                        pgz[i+dn] = c00z[i] * pgz[i] + n * b10[i] * pgz[i-dn] + m * b00[i] * pgz[i-dm];
                                }
                                pgx += dn;
                                pgy += dn;
                                pgz += dn;
                        }
                }
        }
}


/*
 * g0[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
void g0_2e_4d(double *g, const int *ng, const double *rirj, const double *rkrl)
{
        const int nmax = ng[0] - 1;
        const int mmax = ng[1] - 1;
        const int ll = ng[2] - 1;
        const int lj = ng[3] - 1;
        const int li = nmax - lj;
        const int lk = mmax - ll;
        int i, j, k, l, ptr;
        int di, dk, dl, dj;
        double *gx = g;
        double *gy = g + g_size(ng);
        double *gz = g + g_size(ng) * 2;

        extract_dim(ng, &di, &dj, &dk, &dl);

        for (l = 1; l <= ll; l++) {
                // g(...,k,l) = rkrl * g(...,k,l-1) +  g(...,k+1,l-1)
                //n = di * (nmax + 1) * (mmax - l + 1); // dim of (...,0:mmax-l)
                ptr = dl * l;
                for (i = ptr; i < ptr + dl - dk * l; i++) {
                        gx[i] = rkrl[0] * gx[i-dl] + gx[i-dl+dk];
                        gy[i] = rkrl[1] * gy[i-dl] + gy[i-dl+dk];
                        gz[i] = rkrl[2] * gz[i-dl] + gz[i-dl+dk];
                }
        }

        for (j = 1; j <= lj; j++)
                for (l = 0; l <= ll; l++) {
                        ptr = dj * j + dl * l;
                        for (k = 0; k <= lk; k++) {
                                for (i = ptr; i < ptr + dk - di * j; i++) {
                                        gx[i] = rirj[0] * gx[i-dj] + gx[i-dj+di];
                                        gy[i] = rirj[1] * gy[i-dj] + gy[i-dj+di];
                                        gz[i] = rirj[2] * gz[i-dj] + gz[i-dj+di];
                                }
                                ptr += dk;
                        }
                }
}


/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
void g0_2e(double *g, const int *ng,
           const double ai, const double aj,
           const double ak, const double al,
           const double *ri, const double *rj,
           const double *rk, const double *rl, const double fac)
{
        const double aij = ai + aj;
        const double akl = ak + al;
        double a0, a1, fac1, x;
        double u[MXRYSROOTS];
        double w[MXRYSROOTS];
        double c00[MXRYSROOTS*3];
        double c0p[MXRYSROOTS*3];
        double b01[MXRYSROOTS];
        double b00[MXRYSROOTS];
        double b10[MXRYSROOTS];
        double rij[3], rirj[3];
        double rkl[3], rkrl[3];

        rij[0] = (ai * ri[0] + aj * rj[0]) / aij;
        rij[1] = (ai * ri[1] + aj * rj[1]) / aij;
        rij[2] = (ai * ri[2] + aj * rj[2]) / aij;
        rkl[0] = (ak * rk[0] + al * rl[0]) / akl;
        rkl[1] = (ak * rk[1] + al * rl[1]) / akl;
        rkl[2] = (ak * rk[2] + al * rl[2]) / akl;

        a1 = aij * akl;
        a0 = a1 / (aij + akl);
        fac1 = sqrt(a0 / (a1 * a1 * a1)) * fac
                * (M_PI * M_PI * M_PI) * 2 / SQRTPI;
        x = a0 * square_dist(rij, rkl);
        rys_roots(ng[RYS_ROOTS], x, u, w);

        // recursive_g_r12(g, ng, &ng[RYS_ROOTS], aij, akl, rij, rkl, ri, rj, rk, rl, u, w, &fac1);

        g0_2e_bc(ng[RYS_ROOTS], aij, akl, rij, rkl, ri, rk, u, c00, c0p, b01, b00, b10);
        g0_2e_2d(g, ng, c00, c0p, b01, b00, b10, w, fac1);

        rirj[0] = ri[0] - rj[0];
        rirj[1] = ri[1] - rj[1];
        rirj[2] = ri[2] - rj[2];
        rkrl[0] = rk[0] - rl[0];
        rkrl[1] = rk[1] - rl[1];
        rkrl[2] = rk[2] - rl[2];
        g0_2e_4d(g, ng, rirj, rkrl);
}


/*
 * ( \nabla i j | kl )
 */
void nabla1i_2e(double *f, const double *g, const int *ng,
                const int li, const int lj, const int lk, const int ll,
                const double ai)
{
        int i, j, k, l, n, ptr;
        int di, dk, dl, dj;
        const double *gx, *gy, *gz;
        double *fx, *fy, *fz;

        extract_dim(ng, &di, &dj, &dk, &dl);
        extract_fg_xyz(f, g, ng, &gx, &gy, &gz, &fx, &fy, &fz);

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
                                                fx[n] = i * gx[n-di] - 2 * ai * gx[n+di];
                                                fy[n] = i * gy[n-di] - 2 * ai * gy[n+di];
                                                fz[n] = i * gz[n-di] - 2 * ai * gz[n+di];
                                        }
                                        ptr += di;
                                }
                        }
}


/*
 * ( i \nabla j | kl )
 */
void nabla1j_2e(double *f, const double *g, const int *ng,
                const int li, const int lj, const int lk, const int ll,
                const double aj)
{
        int i, j, k, l, ptr;
        int di, dk, dl, dj;
        const double *gx, *gy, *gz;
        double *fx, *fy, *fz;

        extract_dim(ng, &di, &dj, &dk, &dl);
        extract_fg_xyz(f, g, ng, &gx, &gy, &gz, &fx, &fy, &fz);

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
                                        fx[i] = j * gx[i-dj] - 2 * aj * gx[i+dj];
                                        fy[i] = j * gy[i-dj] - 2 * aj * gy[i+dj];
                                        fz[i] = j * gz[i-dj] - 2 * aj * gz[i+dj];
                                }
                                ptr += dk;
                        }
                }
}


/*
 * ( ij | \nabla k l )
 */
void nabla1k_2e(double *f, const double *g, const int *ng,
                const int li, const int lj, const int lk, const int ll,
                const double ak)
{
        int i, j, k, l, ptr;
        int di, dk, dl, dj;
        const double *gx, *gy, *gz;
        double *fx, *fy, *fz;

        extract_dim(ng, &di, &dj, &dk, &dl);
        extract_fg_xyz(f, g, ng, &gx, &gy, &gz, &fx, &fy, &fz);

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
                                        fx[i] = k * gx[i-dk] - 2 * ak * gx[i+dk];
                                        fy[i] = k * gy[i-dk] - 2 * ak * gy[i+dk];
                                        fz[i] = k * gz[i-dk] - 2 * ak * gz[i+dk];
                                }
                                ptr += dk;
                        }
                }
}


/*
 * ( ij | k \nabla l )
 */
void nabla1l_2e(double *f, const double *g, const int *ng,
                const int li, const int lj, const int lk, const int ll,
                const double al)
{
        int i, j, k, l, ptr;
        int di, dk, dl, dj;
        const double *gx, *gy, *gz;
        double *fx, *fy, *fz;

        extract_dim(ng, &di, &dj, &dk, &dl);
        extract_fg_xyz(f, g, ng, &gx, &gy, &gz, &fx, &fy, &fz);

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
                                        fx[i] = l * gx[i-dl] - 2 * al * gx[i+dl];
                                        fy[i] = l * gy[i-dl] - 2 * al * gy[i+dl];
                                        fz[i] = l * gz[i-dl] - 2 * al * gz[i+dl];
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
void x1i_2e(double *f, const double *g, const int *ng,
            const int li, const int lj, const int lk, const int ll,
            const double ri[3])
{
        int i, j, k, l, ptr;
        int di, dk, dl, dj;
        const double *gx, *gy, *gz;
        double *fx, *fy, *fz;

        extract_dim(ng, &di, &dj, &dk, &dl);
        extract_fg_xyz(f, g, ng, &gx, &gy, &gz, &fx, &fy, &fz);

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
void x1j_2e(double *f, const double *g, const int *ng,
            const int li, const int lj, const int lk, const int ll,
            const double rj[3])
{
        int i, j, k, l, ptr;
        int di, dk, dl, dj;
        const double *gx, *gy, *gz;
        double *fx, *fy, *fz;

        extract_dim(ng, &di, &dj, &dk, &dl);
        extract_fg_xyz(f, g, ng, &gx, &gy, &gz, &fx, &fy, &fz);

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
void x1k_2e(double *f, const double *g, const int *ng,
            const int li, const int lj, const int lk, const int ll,
            const double rk[3])
{
        int i, j, k, l, ptr;
        int di, dk, dl, dj;
        const double *gx, *gy, *gz;
        double *fx, *fy, *fz;

        extract_dim(ng, &di, &dj, &dk, &dl);
        extract_fg_xyz(f, g, ng, &gx, &gy, &gz, &fx, &fy, &fz);

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
void x1l_2e(double *f, const double *g, const int *ng,
            const int li, const int lj, const int lk, const int ll,
            const double rl[3])
{
        int i, j, k, l, ptr;
        int di, dk, dl, dj;
        const double *gx, *gy, *gz;
        double *fx, *fy, *fz;

        extract_dim(ng, &di, &dj, &dk, &dl);
        extract_fg_xyz(f, g, ng, &gx, &gy, &gz, &fx, &fy, &fz);

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
