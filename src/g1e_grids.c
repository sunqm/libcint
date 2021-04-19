/*
 * Copyright (C) 2021  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <string.h>
#include <math.h>
#include <assert.h>
#include "cint_bas.h"
#include "misc.h"
#include "g1e.h"
#include "g1e_grids.h"
#include "rys_roots.h"


void CINTinit_int1e_grids_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                                  FINT *atm, FINT natm,
                                  FINT *bas, FINT nbas, double *env)
{
        CINTinit_int1e_EnvVars(envs, ng, shls, atm, natm, bas, nbas, env);
        FINT ngrids = (FINT)env[NGRIDS];
        double *grids = env + (size_t)env[PTR_GRIDS];

        envs->ngrids = ngrids;
        envs->grids = grids;
        envs->common_factor = 2 * M_PI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l);

        FINT nroots = envs->nrys_roots;
        FINT dli, dlj;
        FINT ibase = 1;//envs->li_ceil > envs->lj_ceil;
        if (ibase) {
                dli = envs->li_ceil + envs->lj_ceil + 1;
                dlj = envs->lj_ceil + 1;
        } else {
                dli = envs->li_ceil + 1;
                dlj = envs->li_ceil + envs->lj_ceil + 1;
        }
        envs->g_stride_i = GRID_BLKSIZE * nroots;
        envs->g_stride_j = GRID_BLKSIZE * nroots * dli;
        envs->g_size     = GRID_BLKSIZE * nroots * dli * dlj;
        envs->f_g0_2e = &CINTg0_1e_grids;

        if (ibase) {
                envs->g2d_ijmax = envs->g_stride_i;
                envs->rx_in_rijrx = envs->ri;
                envs->rirj[0] = envs->ri[0] - envs->rj[0];
                envs->rirj[1] = envs->ri[1] - envs->rj[1];
                envs->rirj[2] = envs->ri[2] - envs->rj[2];
                envs->f_g0_2d4d = &CINTg0_1e_grids_igtj;
        } else {
                envs->g2d_ijmax = envs->g_stride_j;
                envs->rx_in_rijrx = envs->rj;
                envs->rirj[0] = envs->rj[0] - envs->ri[0];
                envs->rirj[1] = envs->rj[1] - envs->ri[1];
                envs->rirj[2] = envs->rj[2] - envs->ri[2];
                envs->f_g0_2d4d = &CINTg0_1e_grids_iltj;
        }
}

int CINTg0_1e_grids(double *g, const double fac, CINTEnvVars *envs,
                    double *cache)
{
        FINT ngrids = envs->ngrids;
        FINT bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        FINT nroots = envs->nrys_roots;
        double *grids = envs->grids + envs->grids_offset * 3;
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        double *w = gz;
        double *rij = envs->rij;
        double ubuf[MXRYSROOTS];
        double wbuf[MXRYSROOTS];
        double *u;
        MALLOC_ALIGN8_INSTACK(u, GRID_BLKSIZE*nroots);
        double aij = envs->aij;
        FINT n, i, j, ig;
        double x, fac1;

        for (i = 0; i < nroots; i++) {
                for (ig = 0; ig < bgrids; ig++) {
                        gx[ig+GRID_BLKSIZE*i] = 1;
                        gy[ig+GRID_BLKSIZE*i] = 1;
                }
        }


#ifdef WITH_RANGE_COULOMB
        const double omega = envs->env[PTR_RANGE_OMEGA];
        double theta, sqrt_theta;

        if (omega == 0) {
                fac1 = fac / aij;
                for (ig = 0; ig < bgrids; ig++) {
                        x = aij * CINTsquare_dist(rij, grids+ig*3);
                        CINTrys_roots(nroots, x, ubuf, wbuf);
                        for (i = 0; i < nroots; i++) {
                                u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1);
                                w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                        }
                }
        } else if (omega < 0) { // short-range part of range-separated Coulomb
                theta = omega * omega / (omega * omega + aij);
                sqrt_theta = sqrt(theta);
                fac1 = fac * sqrt_theta / aij;
                for (ig = 0; ig < bgrids; ig++) {
                        x = aij * CINTsquare_dist(rij, grids+ig*3);
                        if (theta * x > envs->expcutoff) {
                                // very small erfc() leads to ~0 weights
                                for (i = 0; i < nroots; i++) {
                                        u[ig+GRID_BLKSIZE*i] = 0;
                                        w[ig+GRID_BLKSIZE*i] = 0;
                                }
                        } else {
                                CINTsr_rys_roots(nroots, x, sqrt_theta, ubuf, wbuf);
                                for (i = 0; i < nroots; i++) {
                                        u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1);
                                        w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                                }
                        }
                }
        } else {  // long-range part of range-separated Coulomb
                theta = omega * omega / (omega * omega + aij);
                aij *= theta;
                fac1 = fac * sqrt(theta) / aij;
                for (ig = 0; ig < bgrids; ig++) {
                        x = aij * CINTsquare_dist(rij, grids+ig*3);
                        CINTrys_roots(nroots, x, ubuf, wbuf);
                        for (i = 0; i < nroots; i++) {
                                u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1 - ubuf[i] * theta);
                                w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                        }
                }
        }
#else
        fac1 = fac / aij;
        for (ig = 0; ig < bgrids; ig++) {
                x = aij * CINTsquare_dist(rij, grids+ig*3);
                CINTrys_roots(nroots, x, ubuf, wbuf);
                for (i = 0; i < nroots; i++) {
                        u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1);
                        w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                }
        }
#endif
        FINT nmax = envs->li_ceil + envs->lj_ceil;
        if (nmax == 0) {
                return 1;
        }

        FINT lj = envs->lj_ceil;
        FINT di = envs->g2d_ijmax;
        FINT dj = envs->g_stride_j;

        double *rg;
        double *rijrx = envs->rijrx;
        double *p0x, *p0y, *p0z;
        double *p1x, *p1y, *p1z;
        double *p2x, *p2y, *p2z;
        double rijrg[3];
        double *rirg;
        MALLOC_ALIGN8_INSTACK(rirg, GRID_BLKSIZE*3);
        double tmp;
        double aij2 = 0.5 / aij;

        for (n = 0; n < nroots; n++) {
                p0x = gx + GRID_BLKSIZE*n;
                p0y = gy + GRID_BLKSIZE*n;
                p0z = gz + GRID_BLKSIZE*n;
                p1x = p0x + di;
                p1y = p0y + di;
                p1z = p0z + di;
                for (ig = 0; ig < bgrids; ig++) {
                        rg = grids + ig * 3;
                        rijrg[0] = rg[0] - rij[0];
                        rijrg[1] = rg[1] - rij[1];
                        rijrg[2] = rg[2] - rij[2];
                        rirg[ig+GRID_BLKSIZE*0] = rijrx[0] + u[ig+GRID_BLKSIZE*n] * rijrg[0];
                        rirg[ig+GRID_BLKSIZE*1] = rijrx[1] + u[ig+GRID_BLKSIZE*n] * rijrg[1];
                        rirg[ig+GRID_BLKSIZE*2] = rijrx[2] + u[ig+GRID_BLKSIZE*n] * rijrg[2];

                        p1x[ig] = rirg[ig+GRID_BLKSIZE*0] * p0x[ig];
                        p1y[ig] = rirg[ig+GRID_BLKSIZE*1] * p0y[ig];
                        p1z[ig] = rirg[ig+GRID_BLKSIZE*2] * p0z[ig];
                }

                for (i = 1; i < nmax; i++) {
                        p0x = gx + GRID_BLKSIZE*n + i * di;
                        p0y = gy + GRID_BLKSIZE*n + i * di;
                        p0z = gz + GRID_BLKSIZE*n + i * di;
                        p1x = p0x + di;
                        p1y = p0y + di;
                        p1z = p0z + di;
                        p2x = p0x - di;
                        p2y = p0y - di;
                        p2z = p0z - di;
                        for (ig = 0; ig < bgrids; ig++) {
                                tmp = i * aij2 * (1 - u[ig+GRID_BLKSIZE*n]);
                                p1x[ig] = tmp * p2x[ig] + rirg[ig+GRID_BLKSIZE*0] * p0x[ig];
                                p1y[ig] = tmp * p2y[ig] + rirg[ig+GRID_BLKSIZE*1] * p0y[ig];
                                p1z[ig] = tmp * p2z[ig] + rirg[ig+GRID_BLKSIZE*2] * p0z[ig];
                        }
                }
        }
        (*envs->f_g0_2d4d)(g, envs);
        return 1;
}

void CINTg0_1e_grids_igtj(double *g, CINTEnvVars *envs)
{
        FINT ngrids = envs->ngrids;
        FINT bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        FINT nroots = envs->nrys_roots;
        FINT nmax = envs->li_ceil + envs->lj_ceil;
        FINT lj = envs->lj_ceil;
        FINT di = envs->g_stride_i;
        FINT dj = envs->g_stride_j;
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        double *p0x, *p0y, *p0z;
        double *p1x, *p1y, *p1z;
        double *p2x, *p2y, *p2z;
        double *rirj = envs->rirj;
        FINT i, j, ig, n;

        for (j = 1; j <= lj; j++) {
        for (i = 0; i <= nmax - j; i++) {
                p0x = gx + j * dj + i * di;
                p0y = gy + j * dj + i * di;
                p0z = gz + j * dj + i * di;
                p1x = p0x - dj;
                p1y = p0y - dj;
                p1z = p0z - dj;
                p2x = p1x + di;
                p2y = p1y + di;
                p2z = p1z + di;
                for (n = 0; n < nroots; n++) {
                for (ig = 0; ig < bgrids; ig++) {
                        p0x[ig+GRID_BLKSIZE*n] = p2x[ig+GRID_BLKSIZE*n] + rirj[0] * p1x[ig+GRID_BLKSIZE*n];
                        p0y[ig+GRID_BLKSIZE*n] = p2y[ig+GRID_BLKSIZE*n] + rirj[1] * p1y[ig+GRID_BLKSIZE*n];
                        p0z[ig+GRID_BLKSIZE*n] = p2z[ig+GRID_BLKSIZE*n] + rirj[2] * p1z[ig+GRID_BLKSIZE*n];
                } }
        } }
}

void CINTg0_1e_grids_iltj(double *g, CINTEnvVars *envs)
{
        FINT ngrids = envs->ngrids;
        FINT bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        FINT nroots = envs->nrys_roots;
        FINT nmax = envs->li_ceil + envs->lj_ceil;
        FINT li = envs->li_ceil;
        FINT di = envs->g_stride_i;
        FINT dj = envs->g_stride_j;
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        double *p0x, *p0y, *p0z;
        double *p1x, *p1y, *p1z;
        double *p2x, *p2y, *p2z;
        double *rirj = envs->rirj;
        FINT i, j, ig, n;

        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) {
                p0x = gx + j * dj + i * di;
                p0y = gy + j * dj + i * di;
                p0z = gz + j * dj + i * di;
                p1x = p0x - di;
                p1y = p0y - di;
                p1z = p0z - di;
                p2x = p1x + dj;
                p2y = p1y + dj;
                p2z = p1z + dj;
                for (n = 0; n < nroots; n++) {
                for (ig = 0; ig < bgrids; ig++) {
                        p0x[ig+GRID_BLKSIZE*n] = p2x[ig+GRID_BLKSIZE*n] + rirj[0] * p1x[ig+GRID_BLKSIZE*n];
                        p0y[ig+GRID_BLKSIZE*n] = p2y[ig+GRID_BLKSIZE*n] + rirj[1] * p1y[ig+GRID_BLKSIZE*n];
                        p0z[ig+GRID_BLKSIZE*n] = p2z[ig+GRID_BLKSIZE*n] + rirj[2] * p1z[ig+GRID_BLKSIZE*n];
                } }
        } }
}

void CINTgout1e_grids(double *gout, double *g, FINT *idx,
                      CINTEnvVars *envs, FINT gout_empty)
{
        FINT ngrids = envs->ngrids;
        FINT bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        FINT nroots = envs->nrys_roots;
        FINT nf = envs->nf;
        FINT i, n, ig;
        double *gx, *gy, *gz;
        double s[GRID_BLKSIZE];

        if (gout_empty) {
                for (n = 0; n < nf; n++, idx+=3) {
                        gx = g + idx[0];
                        gy = g + idx[1];
                        gz = g + idx[2];
                        for (ig = 0; ig < bgrids; ig++) {
                                s[ig] = 0;
                        }
                        for (i = 0; i < nroots; i++) {
                                for (ig = 0; ig < bgrids; ig++) {
                                        s[ig] += gx[ig+GRID_BLKSIZE*i] * gy[ig+GRID_BLKSIZE*i] * gz[ig+GRID_BLKSIZE*i];
                                }
                        }
                        for (ig = 0; ig < bgrids; ig++) {
                                gout[ig+bgrids*n] = s[ig];
                        }
                }
        } else {
                for (n = 0; n < nf; n++, idx+=3) {
                        gx = g + idx[0];
                        gy = g + idx[1];
                        gz = g + idx[2];
                        for (ig = 0; ig < bgrids; ig++) {
                                s[ig] = 0;
                        }
                        for (i = 0; i < nroots; i++) {
                                for (ig = 0; ig < bgrids; ig++) {
                                        s[ig] += gx[ig+GRID_BLKSIZE*i] * gy[ig+GRID_BLKSIZE*i] * gz[ig+GRID_BLKSIZE*i];
                                }
                        }
                        for (ig = 0; ig < bgrids; ig++) {
                                gout[ig+bgrids*n] += s[ig];
                        }
                }
        }
}

#if 0
void CINTnabla1i_1e_grids(double *f, double *g,
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

void CINTnabla1j_1e_grids(double *f, double *g,
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


void CINTx1i_1e_grids(double *f, double *g, double ri[3],
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

void CINTx1j_1e_grids(double *f, double *g, double rj[3],
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
#endif
