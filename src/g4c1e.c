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
#include "g2e.h"

void CINTg0_lj2d_4d(double *g, const CINTEnvVars *envs);
void CINTg0_kj2d_4d(double *g, const CINTEnvVars *envs);
void CINTg0_il2d_4d(double *g, const CINTEnvVars *envs);
void CINTg0_ik2d_4d(double *g, const CINTEnvVars *envs);

FINT CINTinit_int4c1e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
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

        envs->common_factor = 1;
        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = envs->l_l + ng[LINC];
        envs->nrys_roots = 1;
        FINT dli, dlj, dlk, dll;
        FINT ibase = envs->li_ceil > envs->lj_ceil;
        FINT kbase = envs->lk_ceil > envs->ll_ceil;
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
        envs->g_stride_i = 1;
        envs->g_stride_k = dli;
        envs->g_stride_l = dli * dlk;
        envs->g_stride_j = dli * dlk * dll;
        envs->g_size     = dli * dlk * dll * dlj;

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
                        envs->f_g0_2d4d = &CINTg0_ik2d_4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_kj2d_4d;
                }
        } else {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_il2d_4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_lj2d_4d;
                }
        }
        return 0;
}

void CINTg4c1e_index_xyz(FINT *idx, const CINTEnvVars *envs)
{
        CINTg2e_index_xyz(idx, envs);
}


void CINTg4c1e_ovlp(double *g, double fac, const CINTEnvVars *envs)
{
        const FINT li = envs->li_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT lk = envs->lk_ceil;
        const FINT ll = envs->ll_ceil;
        const FINT nmax = li + lj;
        const FINT mmax = lk + ll;
        const FINT db = nmax + mmax + 1;
        const FINT b_size = db*db;
        double buf[3*b_size];
        double *bufx = buf;
        double *bufy = bufx + b_size;
        double *bufz = bufy + b_size;

        double aij = envs->aij;
        double akl = envs->akl;
        const double *rij = envs->rij;
        const double *rkl = envs->rkl;
        const double *ri = envs->rx_in_rijrx;
        const double *rk = envs->rx_in_rklrx;
        double aijkl;
        double r1r12[3];
        double r1r2[3];
        FINT i, j, ptr;

        double *gx = g;
        double *gy = g + envs->g_size;
        double *gz = g + envs->g_size * 2;
        const FINT dn = envs->g2d_ijmax;
        const FINT dm = envs->g2d_klmax;

        aijkl = aij + akl;
        bufx[0] = 1;
        bufy[0] = 1;
        bufz[0] = fac / (aijkl * sqrt(aijkl));

        if (nmax >= mmax) {
                r1r12[0] = ri[0] - (aij*rij[0] + akl*rkl[0]) / aijkl;
                r1r12[1] = ri[1] - (aij*rij[1] + akl*rkl[1]) / aijkl;
                r1r12[2] = ri[2] - (aij*rij[2] + akl*rkl[2]) / aijkl;
                r1r2[0] = ri[0] - rk[0];
                r1r2[1] = ri[1] - rk[1];
                r1r2[2] = ri[2] - rk[2];

                if (nmax+mmax > 0) {
                        bufx[1] = -r1r12[0] * bufx[0];
                        bufy[1] = -r1r12[1] * bufy[0];
                        bufz[1] = -r1r12[2] * bufz[0];
                }

                for (i = 1; i < nmax+mmax; i++) {
                        bufx[i+1] = 0.5 * i / aijkl * bufx[i-1] - r1r12[0] * bufx[i];
                        bufy[i+1] = 0.5 * i / aijkl * bufy[i-1] - r1r12[1] * bufy[i];
                        bufz[i+1] = 0.5 * i / aijkl * bufz[i-1] - r1r12[2] * bufz[i];
                }

                for (j = 1; j <= mmax; j++) {
                        ptr = db * j;
                        for (i = ptr; i <= ptr + nmax+mmax - j; i++) {
                                bufx[i] = bufx[i+1-db] + r1r2[0] * bufx[i-db];
                                bufy[i] = bufy[i+1-db] + r1r2[1] * bufy[i-db];
                                bufz[i] = bufz[i+1-db] + r1r2[2] * bufz[i-db];
                        }
                }

                for (j = 0; j <= mmax; j++) {
                for (i = 0; i <= nmax; i++) {
                        gx[i*dn+j*dm] = bufx[i+j*db];
                        gy[i*dn+j*dm] = bufy[i+j*db];
                        gz[i*dn+j*dm] = bufz[i+j*db];
                } }

        } else {

                r1r12[0] = rk[0] - (aij*rij[0] + akl*rkl[0]) / aijkl;
                r1r12[1] = rk[1] - (aij*rij[1] + akl*rkl[1]) / aijkl;
                r1r12[2] = rk[2] - (aij*rij[2] + akl*rkl[2]) / aijkl;
                r1r2[0] = rk[0] - ri[0];
                r1r2[1] = rk[1] - ri[1];
                r1r2[2] = rk[2] - ri[2];

                if (nmax+mmax > 0) {
                        bufx[1] = -r1r12[0] * bufx[0];
                        bufy[1] = -r1r12[1] * bufy[0];
                        bufz[1] = -r1r12[2] * bufz[0];
                }

                for (i = 1; i < nmax+mmax; i++) {
                        bufx[i+1] = 0.5 * i / aijkl * bufx[i-1] - r1r12[0] * bufx[i];
                        bufy[i+1] = 0.5 * i / aijkl * bufy[i-1] - r1r12[1] * bufy[i];
                        bufz[i+1] = 0.5 * i / aijkl * bufz[i-1] - r1r12[2] * bufz[i];
                }

                for (j = 1; j <= nmax; j++) {
                        ptr = db * j;
                        for (i = ptr; i <= ptr + nmax+mmax - j; i++) {
                                bufx[i] = bufx[i+1-db] + r1r2[0] * bufx[i-db];
                                bufy[i] = bufy[i+1-db] + r1r2[1] * bufy[i-db];
                                bufz[i] = bufz[i+1-db] + r1r2[2] * bufz[i-db];
                        }
                }

                for (j = 0; j <= mmax; j++) {
                for (i = 0; i <= nmax; i++) {
                        gx[i*dn+j*dm] = bufx[i*db+j];
                        gy[i*dn+j*dm] = bufy[i*db+j];
                        gz[i*dn+j*dm] = bufz[i*db+j];
                } }
        }

        (*envs->f_g0_2d4d)(g, envs);
}

