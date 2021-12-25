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
#include "rys_roots.h"


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
        envs->common_factor = 1;
        if (env[PTR_EXPCUTOFF] == 0) {
                envs->expcutoff = EXPCUTOFF;
        } else {
                envs->expcutoff = MAX(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
        }

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_tensor = ng[TENSOR];
        if (ng[SLOT_RYS_ROOTS] > 0) {
                envs->nrys_roots = ng[SLOT_RYS_ROOTS];
        } else {
                envs->nrys_roots = (envs->li_ceil + envs->lj_ceil)/2 + 1;
        }

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
        envs->g_stride_i = envs->nrys_roots;
        envs->g_stride_j = envs->nrys_roots * dli;
        envs->g_size     = envs->nrys_roots * dli * dlj;
        envs->g_stride_k = envs->g_size;
        envs->g_stride_l = envs->g_size;

        assert(i_sh < SHLS_MAX);
        assert(j_sh < SHLS_MAX);
        assert(envs->i_l < ANG_MAX);
        assert(envs->j_l < ANG_MAX);
        assert(bas(ATOM_OF,i_sh) >= 0);
        assert(bas(ATOM_OF,j_sh) >= 0);
        assert(bas(ATOM_OF,i_sh) < natm);
        assert(bas(ATOM_OF,j_sh) < natm);
        assert(envs->nrys_roots < MXRYSROOTS);
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


FINT CINTg1e_ovlp(double *g, CINTEnvVars *envs)
{
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        double aij = envs->ai[0] + envs->aj[0];

        gx[0] = 1;
        gy[0] = 1;
        gz[0] = envs->fac[0] * SQRTPI*M_PI / (aij * sqrt(aij));

        FINT nmax = envs->li_ceil + envs->lj_ceil;
        if (nmax == 0) {
                return 1;
        }

        double *rij = envs->rij;
        double *rirj = envs->rirj;
        FINT lj, di, dj;
        FINT i, j, n, ptr;
        double *rx;
        if (envs->li_ceil > envs->lj_ceil) {
                // li = envs->li_ceil;
                lj = envs->lj_ceil;
                di = envs->g_stride_i;
                dj = envs->g_stride_j;
                rx = envs->ri;
        } else {
                // li = envs->lj_ceil;
                lj = envs->li_ceil;
                di = envs->g_stride_j;
                dj = envs->g_stride_i;
                rx = envs->rj;
        }
        double rijrx[3];
        rijrx[0] = rij[0] - rx[0];
        rijrx[1] = rij[1] - rx[1];
        rijrx[2] = rij[2] - rx[2];

        gx[di] = rijrx[0] * gx[0];
        gy[di] = rijrx[1] * gy[0];
        gz[di] = rijrx[2] * gz[0];

        double aij2 = .5 / aij;
        for (i = 1; i < nmax; i++) {
                gx[(i+1)*di] = i * aij2 * gx[(i-1)*di] + rijrx[0] * gx[i*di];
                gy[(i+1)*di] = i * aij2 * gy[(i-1)*di] + rijrx[1] * gy[i*di];
                gz[(i+1)*di] = i * aij2 * gz[(i-1)*di] + rijrx[2] * gz[i*di];
        }

        for (j = 1; j <= lj; j++) {
                ptr = dj * j;
                for (i = 0, n = ptr; i <= nmax-j; i++, n+=di) {
                        gx[n] = gx[n+di-dj] + rirj[0] * gx[n-dj];
                        gy[n] = gy[n+di-dj] + rirj[1] * gy[n-dj];
                        gz[n] = gz[n+di-dj] + rirj[2] * gz[n-dj];
                }
        }
        return 1;
}

/*
 * Calculate temporary parameter tau for nuclear charge distribution.
 * The charge parameter zeta is defined as    rho(r) = Norm * exp(-zeta*r^2)
 */
double CINTnuc_mod(double aij, FINT nuc_id, FINT *atm, double *env)
{
        double zeta;
        if (nuc_id < 0) {
                zeta = env[PTR_RINV_ZETA];
        } else if (atm(NUC_MOD_OF, nuc_id) == GAUSSIAN_NUC) {
                zeta = env[atm(PTR_ZETA, nuc_id)];
        } else {
                zeta = 0;
        }

        if (zeta > 0) {
                return sqrt(zeta / (aij + zeta));
        } else {
                return 1;
        }
}

FINT CINTg1e_nuc(double *g, CINTEnvVars *envs, FINT nuc_id)
{
        FINT nrys_roots = envs->nrys_roots;
        FINT *atm = envs->atm;
        double *env = envs->env;
        double *rij = envs->rij;
        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        double u[MXRYSROOTS];
        double *w = gz;
        double *cr;
        FINT i, j, n;
        double crij[3];
        double x, fac1;
        double aij = envs->ai[0] + envs->aj[0];
        double tau = CINTnuc_mod(aij, nuc_id, atm, env);

        if (nuc_id < 0) {
                fac1 = 2*M_PI * envs->fac[0] * tau / aij;
                cr = env + PTR_RINV_ORIG;
        } else if (atm(NUC_MOD_OF, nuc_id) == FRAC_CHARGE_NUC) {
                fac1 = 2*M_PI * -env[atm[PTR_FRAC_CHARGE+nuc_id*ATM_SLOTS]] * envs->fac[0] * tau / aij;
                cr = env + atm(PTR_COORD, nuc_id);
        } else {
                fac1 = 2*M_PI * -fabs(atm[CHARGE_OF+nuc_id*ATM_SLOTS]) * envs->fac[0] * tau / aij;
                cr = env + atm(PTR_COORD, nuc_id);
        }
        crij[0] = cr[0] - rij[0];
        crij[1] = cr[1] - rij[1];
        crij[2] = cr[2] - rij[2];
        x = aij * tau * tau * SQUARE(crij);
        CINTrys_roots(nrys_roots, x, u, w);

        for (i = 0; i < nrys_roots; i++) {
                gx[i] = 1;
                gy[i] = 1;
                gz[i] *= fac1;
        }
        FINT nmax = envs->li_ceil + envs->lj_ceil;
        if (nmax == 0) {
                return 1;
        }

        double *p0x, *p0y, *p0z;
        double *p1x, *p1y, *p1z;
        double *p2x, *p2y, *p2z;
        FINT lj, di, dj;
        double *rx;
        if (envs->li_ceil > envs->lj_ceil) {
                // li = envs->li_ceil;
                lj = envs->lj_ceil;
                di = envs->g_stride_i;
                dj = envs->g_stride_j;
                rx = envs->ri;
        } else {
                // li = envs->lj_ceil;
                lj = envs->li_ceil;
                di = envs->g_stride_j;
                dj = envs->g_stride_i;
                rx = envs->rj;
        }
        double rijrx = rij[0] - rx[0];
        double rijry = rij[1] - rx[1];
        double rijrz = rij[2] - rx[2];
        double aij2 = 0.5 / aij;
        double ru, rt, r0, r1, r2;

        p0x = gx + di;
        p0y = gy + di;
        p0z = gz + di;
        p1x = gx - di;
        p1y = gy - di;
        p1z = gz - di;
        for (n = 0; n < nrys_roots; n++) {
                ru = tau * tau * u[n] / (1 + u[n]);
                rt = aij2 - aij2 * ru;
                r0 = rijrx + ru * crij[0];
                r1 = rijry + ru * crij[1];
                r2 = rijrz + ru * crij[2];

                p0x[n] = r0 * gx[n];
                p0y[n] = r1 * gy[n];
                p0z[n] = r2 * gz[n];
                for (i = 1; i < nmax; i++) {
                        p0x[n+i*di] = i * rt * p1x[n+i*di] + r0 * gx[n+i*di];
                        p0y[n+i*di] = i * rt * p1y[n+i*di] + r1 * gy[n+i*di];
                        p0z[n+i*di] = i * rt * p1z[n+i*di] + r2 * gz[n+i*di];
                }
        }

        double rirjx = envs->rirj[0];
        double rirjy = envs->rirj[1];
        double rirjz = envs->rirj[2];
        for (j = 1; j <= lj; j++) {
                p0x = gx + j * dj;
                p0y = gy + j * dj;
                p0z = gz + j * dj;
                p1x = p0x - dj;
                p1y = p0y - dj;
                p1z = p0z - dj;
                p2x = p1x + di;
                p2y = p1y + di;
                p2z = p1z + di;
                for (i = 0; i <= nmax - j; i++) {
                for (n = 0; n < nrys_roots; n++) {
                        p0x[n+i*di] = p2x[n+i*di] + rirjx * p1x[n+i*di];
                        p0y[n+i*di] = p2y[n+i*di] + rirjy * p1y[n+i*di];
                        p0z[n+i*di] = p2z[n+i*di] + rirjz * p1z[n+i*di];
                } }
        }
        return 1;
}

void CINTnabla1i_1e(double *f, double *g,
                    FINT li, FINT lj, FINT lk, CINTEnvVars *envs)
{
        const FINT dj = envs->g_stride_j;
        const FINT dk = envs->g_stride_k;
        const double ai2 = -2 * envs->ai[0];
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
        const double aj2 = -2 * envs->aj[0];
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
        const double ak2 = -2 * envs->ak[0];
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

void CINTprim_to_ctr_0(double *gc, double *gp, double *coeff, size_t nf,
                       FINT nprim, FINT nctr, FINT non0ctr, FINT *sortedidx)
{
        FINT i;
        size_t n;
        double c0;

        for (i = 0; i < nctr; i++) {
                c0 = coeff[nprim* i];
                for (n = 0; n < nf; n++) {
                        gc[nf*i+n] = c0 * gp[n];
                }
        }
}

void CINTprim_to_ctr_1(double *gc, double *gp, double *coeff, size_t nf,
                       FINT nprim, FINT nctr, FINT non0ctr, FINT *sortedidx)
{
        FINT i, j;
        size_t n;
        double c0;

        for (i = 0; i < non0ctr; i++) {
                c0 = coeff[nprim*sortedidx[i]];
                j = sortedidx[i];
                for (n = 0; n < nf; n++) {
                        gc[nf*j+n] += c0 * gp[n];
                }
        }
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
