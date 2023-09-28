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
        FINT ngrids = shls[3] - shls[2];
        double *grids = env + (size_t)env[PTR_GRIDS] + shls[2] * 3;

        envs->ngrids = ngrids;
        envs->grids = grids;
        envs->common_factor = 2 * M_PI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l);

        int rys_order = envs->nrys_roots;
        int nroots = rys_order;
        double omega = env[PTR_RANGE_OMEGA];
        if (omega < 0 && rys_order <= 3) {
                nroots *= 2;
        }
        envs->rys_order = rys_order;
        envs->nrys_roots = nroots;

        FINT dli, dlj;
        FINT ibase = envs->li_ceil > envs->lj_ceil;
        if (ibase) {
                dli = envs->li_ceil + envs->lj_ceil + 1;
                dlj = envs->lj_ceil + 1;
                envs->rirj[0] = envs->ri[0] - envs->rj[0];
                envs->rirj[1] = envs->ri[1] - envs->rj[1];
                envs->rirj[2] = envs->ri[2] - envs->rj[2];
        } else {
                dli = envs->li_ceil + 1;
                dlj = envs->li_ceil + envs->lj_ceil + 1;
                envs->rirj[0] = envs->rj[0] - envs->ri[0];
                envs->rirj[1] = envs->rj[1] - envs->ri[1];
                envs->rirj[2] = envs->rj[2] - envs->ri[2];
        }
        envs->g_stride_i = GRID_BLKSIZE * nroots;
        envs->g_stride_j = GRID_BLKSIZE * nroots * dli;
        envs->g_size     = GRID_BLKSIZE * nroots * dli * dlj;
        envs->g_stride_k = envs->g_size;
        envs->g_stride_l = envs->g_size;
}

#define RGSQUARE(r, ig)     (r[ig+GRID_BLKSIZE*0]*r[ig+GRID_BLKSIZE*0] + \
                             r[ig+GRID_BLKSIZE*1]*r[ig+GRID_BLKSIZE*1] + \
                             r[ig+GRID_BLKSIZE*2]*r[ig+GRID_BLKSIZE*2])

FINT CINTg0_1e_grids(double *g, double cutoff,
                     CINTEnvVars *envs, double *cache, double *gridsT)
{
        FINT ngrids = envs->ngrids;
        FINT bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        int nroots = envs->nrys_roots;
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        double *w = gz;
        double *rij = envs->rij;
        double ubuf[MXRYSROOTS];
        double wbuf[MXRYSROOTS];
        double *u;
        MALLOC_ALIGN8_INSTACK(u, GRID_BLKSIZE*nroots);
        double *rijrg;
        MALLOC_ALIGN8_INSTACK(rijrg, GRID_BLKSIZE*3);
        double aij = envs->ai[0] + envs->aj[0];
        FINT n, i, j, ig;
        double x, fac1;

        for (i = 0; i < nroots; i++) {
                for (ig = 0; ig < bgrids; ig++) {
                        gx[ig+GRID_BLKSIZE*i] = 1;
                        gy[ig+GRID_BLKSIZE*i] = 1;
                }
        }
#pragma GCC ivdep
        for (ig = 0; ig < bgrids; ig++) {
                rijrg[ig+GRID_BLKSIZE*0] = gridsT[ig+GRID_BLKSIZE*0] - rij[0];
                rijrg[ig+GRID_BLKSIZE*1] = gridsT[ig+GRID_BLKSIZE*1] - rij[1];
                rijrg[ig+GRID_BLKSIZE*2] = gridsT[ig+GRID_BLKSIZE*2] - rij[2];
        }

        double omega = envs->env[PTR_RANGE_OMEGA];
        double zeta = envs->env[PTR_RINV_ZETA];
        double omega2, theta, sqrt_theta, a0, tau2;

        assert(zeta >= 0);
        if (omega == 0. && zeta == 0.) {
                fac1 = envs->fac[0] / aij;
                for (ig = 0; ig < bgrids; ig++) {
                        x = aij * RGSQUARE(rijrg, ig);
                        CINTrys_roots(nroots, x, ubuf, wbuf);
                        for (i = 0; i < nroots; i++) {
                                // transform to t^2
                                u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1);
                                w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                        }
                }
        } else if (omega < 0.) { // short-range part of range-separated Coulomb
                a0 = aij;
                fac1 = envs->fac[0] / aij;
                if (zeta == 0.) {
                        tau2 = 1.;
                        omega2 = omega * omega;
                        theta = omega2 / (omega2 + aij);
                } else { // zeta > 0.
                        tau2 = zeta / (zeta + aij);
                        a0 *= tau2;
                        fac1 *= sqrt(tau2);
                        omega2 = omega * omega;
                        theta = omega2 / (omega2 + a0);
                }
                sqrt_theta = sqrt(theta);
                // very small erfc() leads to ~0 weights. They can cause
                // numerical issue in sr_rys_roots. Use this cutoff as a
                // temporary solution to avoid numerical issues
                double temp_cutoff = MIN(cutoff, EXPCUTOFF_SR);
                int rorder = envs->rys_order;
                double tau_theta, fac_theta;
                for (ig = 0; ig < bgrids; ig++) {
                        x = a0 * RGSQUARE(rijrg, ig);
                        if (theta * x > temp_cutoff) {
                                // very small erfc() leads to ~0 weights
                                for (i = 0; i < nroots; i++) {
                                        u[ig+GRID_BLKSIZE*i] = 0;
                                        w[ig+GRID_BLKSIZE*i] = 0;
                                }
                        } else if (rorder == nroots) {
                                CINTsr_rys_roots(nroots, x, sqrt_theta, ubuf, wbuf);
                                for (i = 0; i < nroots; i++) {
                                        u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1) * tau2;
                                        w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                                }
                        } else {
                                tau_theta = tau2 * theta;
                                fac_theta = fac1 * -sqrt_theta;
                                CINTrys_roots(rorder, x, ubuf, wbuf);
                                CINTrys_roots(rorder, theta*x, ubuf+rorder, wbuf+rorder);
                                for (i = 0; i < rorder; i++) {
                                        u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1) * tau2;
                                        w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                                        u[ig+GRID_BLKSIZE*(i+rorder)] = ubuf[i+rorder] / (ubuf[i+rorder] + 1) * tau_theta;
                                        w[ig+GRID_BLKSIZE*(i+rorder)] = wbuf[i+rorder] * fac_theta;
                                }
                        }
                }
        } else {
                // * long-range part of range-separated Coulomb
                // * or Gaussian charge model, with rho(r) = Norm exp(-zeta*r^2)
                a0 = aij;
                fac1 = envs->fac[0] / aij;
                if (zeta == 0.) { // omega > 0.
                        omega2 = omega * omega;
                        theta = omega2 / (omega2 + aij);
                        a0 *= theta;
                        fac1 *= sqrt(theta);
                } else if (omega == 0.) { // zeta > 0.
                        theta = zeta / (zeta + aij);
                        a0 *= theta;
                        fac1 *= sqrt(theta);
                } else { // omega > 0. && zeta > 0.
                        omega2 = omega * omega;
                        theta = omega2*zeta / (omega2*zeta + (zeta+omega2)*aij);
                        a0 *= theta;
                        fac1 *= sqrt(theta);
                }
                for (ig = 0; ig < bgrids; ig++) {
                        x = a0 * RGSQUARE(rijrg, ig);
                        CINTrys_roots(nroots, x, ubuf, wbuf);
                        for (i = 0; i < nroots; i++) {
                                // u stores t^2 = tau^2 * theta
                                u[ig+GRID_BLKSIZE*i] = ubuf[i] / (ubuf[i] + 1) * theta;
                                w[ig+GRID_BLKSIZE*i] = wbuf[i] * fac1;
                        }
                }
        }
        FINT nmax = envs->li_ceil + envs->lj_ceil;
        if (nmax == 0) {
                return 1;
        }

        double *rirj = envs->rirj;
        FINT lj, di, dj;
        double *rx;
        if (envs->li_ceil > envs->lj_ceil) {
                //li = envs->li_ceil;
                lj = envs->lj_ceil;
                di = envs->g_stride_i;
                dj = envs->g_stride_j;
                rx = envs->ri;
        } else {
                //li = envs->lj_ceil;
                lj = envs->li_ceil;
                di = envs->g_stride_j;
                dj = envs->g_stride_i;
                rx = envs->rj;
        }
        double rijrx[3];
        rijrx[0] = rij[0] - rx[0];
        rijrx[1] = rij[1] - rx[1];
        rijrx[2] = rij[2] - rx[2];

        double *p0x, *p0y, *p0z;
        double *p1x, *p1y, *p1z;
        double *p2x, *p2y, *p2z;
        double *t2;
        MALLOC_ALIGN8_INSTACK(t2, GRID_BLKSIZE*4);
        double *rirgx = t2 + GRID_BLKSIZE;
        double *rirgy = rirgx + GRID_BLKSIZE;
        double *rirgz = rirgy + GRID_BLKSIZE;
        double aij2 = 0.5 / aij;
        double tx, ty, tz;

        for (n = 0; n < nroots; n++) {
                p0x = gx + GRID_BLKSIZE*n;
                p0y = gy + GRID_BLKSIZE*n;
                p0z = gz + GRID_BLKSIZE*n;
                p1x = p0x + di;
                p1y = p0y + di;
                p1z = p0z + di;
#pragma GCC ivdep
                for (ig = 0; ig < bgrids; ig++) {
                        rirgx[ig] = rijrx[0] + u[ig+GRID_BLKSIZE*n] * rijrg[ig+GRID_BLKSIZE*0];
                        rirgy[ig] = rijrx[1] + u[ig+GRID_BLKSIZE*n] * rijrg[ig+GRID_BLKSIZE*1];
                        rirgz[ig] = rijrx[2] + u[ig+GRID_BLKSIZE*n] * rijrg[ig+GRID_BLKSIZE*2];
                        p1x[ig] = rirgx[ig] * p0x[ig];
                        p1y[ig] = rirgy[ig] * p0y[ig];
                        p1z[ig] = rirgz[ig] * p0z[ig];
                }
                if (nmax > 0) {
                        for (ig = 0; ig < bgrids; ig++) {
                                t2[ig] = aij2 * (1 - u[ig+GRID_BLKSIZE*n]);
                        }
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
#pragma GCC ivdep
                        for (ig = 0; ig < bgrids; ig++) {
                                p1x[ig] = i * t2[ig] * p2x[ig] + rirgx[ig] * p0x[ig];
                                p1y[ig] = i * t2[ig] * p2y[ig] + rirgy[ig] * p0y[ig];
                                p1z[ig] = i * t2[ig] * p2z[ig] + rirgz[ig] * p0z[ig];
                        }
                }
        }

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
#pragma GCC ivdep
                for (ig = 0; ig < bgrids; ig++) {
                        p0x[ig+GRID_BLKSIZE*n] = p2x[ig+GRID_BLKSIZE*n] + rirj[0] * p1x[ig+GRID_BLKSIZE*n];
                        p0y[ig+GRID_BLKSIZE*n] = p2y[ig+GRID_BLKSIZE*n] + rirj[1] * p1y[ig+GRID_BLKSIZE*n];
                        p0z[ig+GRID_BLKSIZE*n] = p2z[ig+GRID_BLKSIZE*n] + rirj[2] * p1z[ig+GRID_BLKSIZE*n];
                } }
        } }
        return 1;
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

void CINTnabla1i_grids(double *f, double *g,
                       FINT li, FINT lj, CINTEnvVars *envs)
{
        FINT ngrids = envs->ngrids;
        FINT bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        FINT nroots = envs->nrys_roots;
        const FINT di = envs->g_stride_i;
        const FINT dj = envs->g_stride_j;
        const double ai2 = -2 * envs->ai[0];
        FINT i, j, n, ig, ptr;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;

        for (j = 0; j <= lj; j++) {
                //f(...,0,...) = -2*ai*g(...,1,...)
                for (n = 0; n < nroots; n++) {
                        ptr = dj * j + n * GRID_BLKSIZE;
#pragma GCC ivdep
                        for (ig = ptr; ig < ptr+bgrids; ig++) {
                                fx[ig] = ai2 * gx[ig+di];
                                fy[ig] = ai2 * gy[ig+di];
                                fz[ig] = ai2 * gz[ig+di];
                        }
                }
                //f(...,i,...) = i*g(...,i-1,...)-2*ai*g(...,i+1,...)
                for (i = 1; i <= li; i++) {
                for (n = 0; n < nroots; n++) {
                        ptr = dj * j + di * i + n * GRID_BLKSIZE;
#pragma GCC ivdep
                        for (ig = ptr; ig < ptr+bgrids; ig++) {
                                fx[ig] = i * gx[ig-di] + ai2 * gx[ig+di];
                                fy[ig] = i * gy[ig-di] + ai2 * gy[ig+di];
                                fz[ig] = i * gz[ig-di] + ai2 * gz[ig+di];
                        }
                } }
        }
}

void CINTnabla1j_grids(double *f, double *g,
                       FINT li, FINT lj, CINTEnvVars *envs)
{
        FINT ngrids = envs->ngrids;
        FINT bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        FINT nroots = envs->nrys_roots;
        const FINT di = envs->g_stride_i;
        const FINT dj = envs->g_stride_j;
        const double aj2 = -2 * envs->aj[0];
        FINT i, j, n, ig, ptr;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;

        //f(...,0,...) = -2*aj*g(...,1,...)
        for (i = 0; i <= li; i++) {
        for (n = 0; n < nroots; n++) {
                ptr = di * i + n * GRID_BLKSIZE;
#pragma GCC ivdep
                for (ig = ptr; ig < ptr+bgrids; ig++) {
                        fx[ig] = aj2 * gx[ig+dj];
                        fy[ig] = aj2 * gy[ig+dj];
                        fz[ig] = aj2 * gz[ig+dj];
                }
        } }
        //f(...,j,...) = j*g(...,j-1,...)-2*aj*g(...,j+1,...)
        for (j = 1; j <= lj; j++) {
        for (i = 0; i <= li; i++) {
        for (n = 0; n < nroots; n++) {
                ptr = dj * j + di * i + n * GRID_BLKSIZE;
#pragma GCC ivdep
                for (ig = ptr; ig < ptr+bgrids; ig++) {
                        fx[ig] = j * gx[ig-dj] + aj2 * gx[ig+dj];
                        fy[ig] = j * gy[ig-dj] + aj2 * gy[ig+dj];
                        fz[ig] = j * gz[ig-dj] + aj2 * gz[ig+dj];
                }
        } } }
}


void CINTx1i_grids(double *f, double *g, double *ri,
                   FINT li, FINT lj, CINTEnvVars *envs)
{
        FINT ngrids = envs->ngrids;
        FINT bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        FINT nroots = envs->nrys_roots;
        FINT i, j, n, ig, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dj = envs->g_stride_j;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;

        for (j = 0; j <= lj; j++) {
        for (i = 0; i <= li; i++) {
        for (n = 0; n < nroots; n++) {
                ptr = dj * j + di * i + n * GRID_BLKSIZE;
#pragma GCC ivdep
                for (ig = ptr; ig < ptr+bgrids; ig++) {
                        fx[ig] = gx[ig+di] + ri[0] * gx[ig];
                        fy[ig] = gy[ig+di] + ri[1] * gy[ig];
                        fz[ig] = gz[ig+di] + ri[2] * gz[ig];
                }
        } } }
}

void CINTx1j_grids(double *f, double *g, double *rj,
                   FINT li, FINT lj, CINTEnvVars *envs)
{
        FINT ngrids = envs->ngrids;
        FINT bgrids = MIN(ngrids - envs->grids_offset, GRID_BLKSIZE);
        FINT nroots = envs->nrys_roots;
        FINT i, j, n, ig, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dj = envs->g_stride_j;
        const double *gx = g;
        const double *gy = g + envs->g_size;
        const double *gz = g + envs->g_size * 2;
        double *fx = f;
        double *fy = f + envs->g_size;
        double *fz = f + envs->g_size * 2;

        for (j = 0; j <= lj; j++) {
        for (i = 0; i <= li; i++) {
        for (n = 0; n < nroots; n++) {
                ptr = dj * j + di * i + n * GRID_BLKSIZE;
#pragma GCC ivdep
                for (ig = ptr; ig < ptr+bgrids; ig++) {
                        fx[ig] = gx[ig+dj] + rj[0] * gx[ig];
                        fy[ig] = gy[ig+dj] + rj[1] * gy[ig];
                        fz[ig] = gz[ig+dj] + rj[2] * gz[ig];
                }
        } } }
}
