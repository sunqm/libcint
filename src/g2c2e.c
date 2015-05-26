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
#include "g2c2e.h"

#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *GX = G; \
        type *GY = G + envs->g_size; \
        type *GZ = G + envs->g_size * 2

void CINTg0_2c2e_ik2d(double *g, const CINTEnvVars *envs,struct _BC *bc);
static void CINTset_g2c2e_params(CINTEnvVars *envs);

FINT CINTinit_int2c2e_EnvVars(CINTEnvVars *envs, const FINT *ng, const FINT *shls,
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
        const FINT k_sh = shls[1];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->i_prim = bas(NPRIM_OF, i_sh);
        envs->k_prim = bas(NPRIM_OF, k_sh);
        envs->i_ctr = bas(NCTR_OF, i_sh);
        envs->k_ctr = bas(NCTR_OF, k_sh);
        envs->nfi = CINTlen_cart(envs->i_l);
        envs->nfk = CINTlen_cart(envs->k_l);
        envs->nf = envs->nfi * envs->nfk;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));

// initialize j_l, j_ctr, nfj because they are used in c2s_sph_1e and
// CINTg1e_index_xyz
        envs->j_l = envs->k_l;
        envs->j_ctr = envs->k_ctr;
        envs->nfj = envs->nfk;

        envs->common_factor = (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->k_l);

        envs->gbits = ng[GSHIFT];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = 0;
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = 0;
        envs->nrys_roots =(envs->li_ceil + envs->lk_ceil)/2 + 1;

        assert(i_sh < SHLS_MAX);
        assert(k_sh < SHLS_MAX);
        assert(envs->i_l < ANG_MAX);
        assert(envs->k_l < ANG_MAX);
        assert(envs->i_ctr < NCTR_MAX);
        assert(envs->k_ctr < NCTR_MAX);
        assert(envs->i_prim < NPRIM_MAX);
        assert(envs->k_prim < NPRIM_MAX);
        assert(envs->i_prim >= envs->i_ctr);
        assert(envs->k_prim >= envs->k_ctr);
        assert(bas(ATOM_OF,i_sh) >= 0);
        assert(bas(ATOM_OF,k_sh) >= 0);
        assert(bas(ATOM_OF,i_sh) < natm);
        assert(bas(ATOM_OF,k_sh) < natm);
        assert(envs->nrys_roots < MXRYSROOTS);

        CINTset_g2c2e_params(envs);
        return 0;
}

/* set strides and parameters for g0_2d4d algorithm */
static void CINTset_g2c2e_params(CINTEnvVars *envs)
{
        FINT dli, dlj, dlk, dll;
        FINT ibase = 1;
        FINT kbase = 1;
        if (envs->nrys_roots <= 2) { // use the fully optimized lj_4d algorithm
                ibase = 0;
                kbase = 0;
        }

        if (ibase) {
                dli = envs->li_ceil + 1;
                dlj = 1;
        } else { // to use _g0_lj_4d_xxxx functions, dll can be > 1
                dli = envs->li_ceil + 1;
                dlj = dli;
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
// initialize g_stride_j to envs->g_stride_k because it's needed in CINTg1e_index_xyz
        envs->g_stride_j = envs->g_stride_k;
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
                // in CINTg0_2c2e_ik2d, g2d_klmax will not be used
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
                envs->rirj[0] = envs->ri[0];
                envs->rirj[1] = envs->ri[1];
                envs->rirj[2] = envs->ri[2];
                envs->rijrx[0] = 0; // since rj=0, aj=0
                envs->rijrx[1] = 0;
                envs->rijrx[2] = 0;
        } else {
                // in CINTg0_2c2e_ik2d, g2d_ijmax will not be used
                envs->g2d_ijmax = 0;
                envs->rirj[0] = -envs->ri[0];
                envs->rirj[1] = -envs->ri[1];
                envs->rirj[2] = -envs->ri[2];
                envs->rijrx[0] = envs->ri[0];
                envs->rijrx[1] = envs->ri[1];
                envs->rijrx[2] = envs->ri[2];
        }
        envs->rij[0] = envs->ri[0];
        envs->rij[1] = envs->ri[1];
        envs->rij[2] = envs->ri[2];

        envs->f_g0_2d4d = &CINTg0_2c2e_ik2d;
}


/*
 * ( \nabla i | k )
 */
void CINTnabla1i_2c2e(double *f, const double *g,
                      const FINT li, const FINT lk, const CINTEnvVars *envs)
{
        FINT i, k, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
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
        for (k = 0; k <= lk; k++) {
                ptr = dk * k;
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
 * ( i | k )
 */
void CINTnabla1k_2c2e(double *f, const double *g,
                      const FINT li, const FINT lk, const CINTEnvVars *envs)
{
        FINT i, k, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
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
        ptr = 0;
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
                ptr = dk * k;
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


/*
 * ( x^1 i | k )
 * ri is the shift from the center R_O to the center of |i>
 * r - R_O = (r-R_i) + ri, ri = R_i - R_O
 */
void CINTx1i_2c2e(double *f, const double *g, const double *ri,
                  const FINT li, const FINT lk, const CINTEnvVars *envs)
{
        FINT i, k, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + di;
        const double *p1y = gy + di;
        const double *p1z = gz + di;
        for (k = 0; k <= lk; k++) {
                //f(...,0:li,...) = g(...,1:li+1,...) + ri(1)*g(...,0:li,...)
                ptr = dk * k;
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
 * ( i | x^1 k )
 */
void CINTx1k_2c2e(double *f, const double *g, const double *rk,
                  const FINT li, const FINT lk, const CINTEnvVars *envs)
{
        FINT i, k, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + dk;
        const double *p1y = gy + dk;
        const double *p1z = gz + dk;
        for (k = 0; k <= lk; k++) {
                // f(...,0:lk,...) = g(...,1:lk+1,...) + rk(1)*g(...,0:lk,...)
                ptr = dk * k;
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

