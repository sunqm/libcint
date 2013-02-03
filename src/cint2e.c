/*
 * File: cint2e.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO integrals
 */

#include <stdlib.h>
#include <math.h>
#include "cint_bas.h"
#include "g2e.h"
#include "cint2e.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"


static void _prim_to_ctr(double *gc, const int nf, const double *gp, const int inc, 
                 const int shl, const int i_pgto, 
                 const int *bas, const double *env)
{
        const int INC1 = 1;
        const int off = bas(PTR_COEFF, shl) + i_pgto;
        const double *penv = env + off;

        if (bas(NPRIM_OF, shl) == 1) {
                dscal_(&nf, penv, gc, &INC1);
        } else {
                prim_to_ctr(gc, nf, gp, inc, shl, i_pgto, bas, env);
        }
}

int cint2e_loop(double *gctr, const int *ng, const double fac,
                void (*const f_gout)(),
                const int *shls, const int *atm, const int *bas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int k_sh = shls[2];
        const int l_sh = shls[3];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        const int k_l = bas(ANG_OF, k_sh);
        const int l_l = bas(ANG_OF, l_sh);
        const int i_ctr = bas(NCTR_OF, i_sh);
        const int j_ctr = bas(NCTR_OF, j_sh);
        const int k_ctr = bas(NCTR_OF, k_sh);
        const int l_ctr = bas(NCTR_OF, l_sh);
        const int nfi = len_cart(i_l);
        const int nfj = len_cart(j_l);
        const int nfk = len_cart(k_l);
        const int nfl = len_cart(l_l);
        const int nf = nfi * nfk * nfl * nfj;
        const int n_comp = ng[POS_E1] * ng[POS_E2] * ng[TENSOR];
        const double *ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        const double *rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        const double *rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));
        const double *rl = env + atm(PTR_COORD, bas(ATOM_OF, l_sh));
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        const double *ak = env + bas(PTR_EXP, k_sh);
        const double *al = env + bas(PTR_EXP, l_sh);
        int ip, jp, kp, lp, n;
        int has_value = 0;
        int *const idx = (int *)malloc(sizeof(int) * nf * 3);
        double aij, dij, eij, rrij;
        double akl, dkl, ekl, rrkl;
        double *const g = (double *)malloc(sizeof(double) * g_size(ng) * 3 * ((1<<ng[GSHIFT])+1));  // (irys,i,j,k,l,coord,0:1);
        double *const gctrk = (double *)malloc(sizeof(double) * nf * i_ctr * j_ctr * k_ctr * n_comp);
        double *gout, *gctri, *gctrj;

        if (bas(NPRIM_OF, k_sh) == 1) {
                gctrj = gctrk;
        } else {
                gctrj = (double *)malloc(sizeof(double) * nf * i_ctr * j_ctr * n_comp);
        }
        if (bas(NPRIM_OF, j_sh) == 1) {
                gctri = gctrj;
        } else {
                gctri = (double *)malloc(sizeof(double) * nf * i_ctr * n_comp);
        }
        if (bas(NPRIM_OF, i_sh) == 1) {
                gout = gctri;
        } else {
                gout  = (double *)malloc(sizeof(double) * nf * n_comp);
        }

        g2e_index_xyz(idx, ng, shls, bas);

        rrij = square_dist(ri, rj);
        rrkl = square_dist(rk, rl);

        n = nf * i_ctr * j_ctr * k_ctr * l_ctr * n_comp;
        dset0(n, gctr); 
        for (lp = 0; lp < bas(NPRIM_OF, l_sh); lp++) {
                if (bas(NPRIM_OF, k_sh) > 1) {
                        n = nf * i_ctr * j_ctr * k_ctr * n_comp;
                        dset0(n, gctrk); 
                }
                for (kp = 0; kp < bas(NPRIM_OF, k_sh); kp++) {
                        akl = ak[kp] + al[lp];
                        ekl = (ak[kp] * al[lp] / akl) * rrkl;
                        if (ekl > EXPCUTOFF) {
                                if (bas(NPRIM_OF, k_sh) == 1) {
                                        n = nf * i_ctr * j_ctr * n_comp;
                                        dset0(n, gctrk); 
                                }
                                continue;
                        }
                        dkl = exp(-ekl);

                        if (bas(NPRIM_OF, j_sh) > 1) {
                                n = nf * i_ctr * j_ctr * n_comp;
                                dset0(n, gctrj); 
                        }
                        for (jp = 0; jp < bas(NPRIM_OF, j_sh); jp++) {
                                if (bas(NPRIM_OF, i_sh) > 1) {
                                        n = nf * i_ctr * n_comp;
                                        dset0(n, gctri);
                                }
                                for (ip = 0; ip < bas(NPRIM_OF, i_sh); ip++) {
                                        aij = ai[ip] + aj[jp];
                                        eij = (ai[ip] * aj[jp] / aij) * rrij;
                                        if (eij > EXPCUTOFF || eij+ekl > EXPCUTOFF) {
                                                if (bas(NPRIM_OF, i_sh) == 1) {
                                                        n = nf * n_comp;
                                                        dset0(n, gctri);
                                                }
                                                continue;
                                        }
                                        has_value = 1;
                                        dij = exp(-eij);

                                        g0_2e(g, ng, ai[ip], aj[jp], ak[kp], al[lp],
                                              ri, rj, rk, rl, dij * dkl * fac);

                                        (*f_gout)(g, ng, gout, nf, idx,
                                                  ai[ip], aj[jp], ak[kp], al[lp],
                                                  shls, atm, bas, env);

                                        n = nf * n_comp;
                                        _prim_to_ctr(gctri, n, gout, 1, i_sh, ip, bas, env);
                                }
                                n = nf * i_ctr * n_comp;
                                _prim_to_ctr(gctrj, n, gctri, 1, j_sh, jp, bas, env);
                        }
                        n = nf * i_ctr * j_ctr * n_comp;
                        _prim_to_ctr(gctrk, n, gctrj, 1, k_sh, kp, bas, env);
                }
                n = nf * i_ctr * j_ctr * k_ctr;
                prim_to_ctr(gctr, n, gctrk, n_comp, l_sh, lp, bas, env);
        }
        free(g);
        free(idx);
        if (bas(NPRIM_OF, i_sh) > 1) {
                free(gout);
        }
        if (bas(NPRIM_OF, j_sh) > 1) {
                free(gctri);
        }
        if (bas(NPRIM_OF, k_sh) > 1) {
                free(gctrj);
        }
        free(gctrk);

        return has_value;
}


int cint2e_drv(double *opkijl, int *ng, const double fac,
               void (*const f_gout)(), void (*const f_e1_c2s)(),
               void (*const f_e2_c2s)(),
               const int *shls, const int *atm, const int natm,
               const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int k_sh = shls[2];
        const int l_sh = shls[3];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        const int k_l = bas(ANG_OF, k_sh);
        const int l_l = bas(ANG_OF, l_sh);
        const int i_ctr = bas(NCTR_OF, i_sh);
        const int j_ctr = bas(NCTR_OF, j_sh);
        const int k_ctr = bas(NCTR_OF, k_sh);
        const int l_ctr = bas(NCTR_OF, l_sh);
        const int nfi = len_cart(i_l);
        const int nfj = len_cart(j_l);
        const int nfk = len_cart(k_l);
        const int nfl = len_cart(l_l);
        const int nf = nfi * nfk * nfl * nfj;
        const int nc = nf * i_ctr * j_ctr * k_ctr * l_ctr * ng[POS_E1];
        int ip, jp, kp, lp, nop;
        int n, n1, n2;
        int has_value;
        double *const gctr = (double *)malloc(sizeof(double) * nc * ng[POS_E2] * ng[TENSOR]);
        double *pgctr = gctr;
        double *opij;

        ng[RYS_ROOTS] = (ng[0] + ng[1]) / 2;

        has_value = cint2e_loop(gctr, ng, fac, f_gout,
                                shls, atm, bas, env);

        int (*fcgtos)();
        if (f_e1_c2s == &c2s_sph_2e1) {
                fcgtos = cgtos_spheric;
        } else if (f_e1_c2s == &c2s_cart_2e1) {
                fcgtos = cgtos_cart;
        } else {
                fcgtos = cgtos_spinor;
        }
        ip = (*fcgtos)(i_sh, bas);
        jp = (*fcgtos)(j_sh, bas);
        kp = (*fcgtos)(k_sh, bas);
        lp = (*fcgtos)(l_sh, bas);
        nop = ip * jp * kp * lp;

        if (f_e1_c2s == &c2s_sph_2e1 || f_e1_c2s == &c2s_cart_2e1) {
                if (!has_value) {
                        dset0(nop * ng[TENSOR], opkijl);
                } else {
                        for (n = 0; n < ng[TENSOR]; n++) {
                                (*f_e1_c2s)(opkijl, pgctr, shls, bas);
                                opkijl += nop;
                                pgctr += nc;
                        }
                }
        } else {
                if (!has_value) {
                        dset0(nop * OF_CMPLX * ng[TENSOR], opkijl);
                } else {
                        n1 = ip * nfk * k_ctr * nfl * l_ctr * jp * OF_CMPLX;
                        opij = (double *)malloc(sizeof(double) * n1 * ng[POS_E2]);
                        for (n = 0; n < ng[TENSOR]; n++) {
                                for (n2 = 0; n2 < ng[POS_E2]; n2++) {
                                        (*f_e1_c2s)(opij + n1 * n2, pgctr, shls, bas);
                                        pgctr += nc;
                                }
                                (*f_e2_c2s)(opkijl, opij, shls, bas);
                                opkijl += nop * OF_CMPLX;
                        }
                        free(opij);
                }
        }
        free(gctr);
        return has_value;
}


/*
 * <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
void gout2e(double *g, const int *ng,
            double *gout, const int nf, const int *idx,
            const double ai, const double aj, const double ak, const double al,
            const int *shls, const int *atm, const int *bas, const double *env)
{
        int i, ix, iy, iz, n;
        const int *idy = idx + nf;
        const int *idz = idx + nf * 2;

        switch (ng[RYS_ROOTS]) {
                case 1:
                        for (n = 0; n < nf; n++) {
                                ix = idx[n];
                                iy = idy[n];
                                iz = idz[n];
                                gout[n] = g[ix] * g[iy] * g[iz];
                        }
                        break;
                case 2:
                        for (n = 0; n < nf; n++) {
                                ix = idx[n];
                                iy = idy[n];
                                iz = idz[n];
                                gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                        + g[ix+1] * g[iy+1] * g[iz+1];
                        }
                        break;
                case 3:
                        for (n = 0; n < nf; n++) {
                                ix = idx[n];
                                iy = idy[n];
                                iz = idz[n];
                                gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                        + g[ix+1] * g[iy+1] * g[iz+1]
                                        + g[ix+2] * g[iy+2] * g[iz+2];
                        }
                        break;
                case 4:
                        for (n = 0; n < nf; n++) {
                                ix = idx[n];
                                iy = idy[n];
                                iz = idz[n];
                                gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                        + g[ix+1] * g[iy+1] * g[iz+1]
                                        + g[ix+2] * g[iy+2] * g[iz+2]
                                        + g[ix+3] * g[iy+3] * g[iz+3];
                        }
                        break;
                case 5:
                        for (n = 0; n < nf; n++) {
                                ix = idx[n];
                                iy = idy[n];
                                iz = idz[n];
                                gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                        + g[ix+1] * g[iy+1] * g[iz+1]
                                        + g[ix+2] * g[iy+2] * g[iz+2]
                                        + g[ix+3] * g[iy+3] * g[iz+3]
                                        + g[ix+4] * g[iy+4] * g[iz+4];
                        }
                        break;
                case 6:
                        for (n = 0; n < nf; n++) {
                                ix = idx[n];
                                iy = idy[n];
                                iz = idz[n];
                                gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                        + g[ix+1] * g[iy+1] * g[iz+1]
                                        + g[ix+2] * g[iy+2] * g[iz+2]
                                        + g[ix+3] * g[iy+3] * g[iz+3]
                                        + g[ix+4] * g[iy+4] * g[iz+4]
                                        + g[ix+5] * g[iy+5] * g[iz+5];
                        }
                        break;
                case 7:
                        for (n = 0; n < nf; n++) {
                                ix = idx[n];
                                iy = idy[n];
                                iz = idz[n];
                                gout[n] = g[ix  ] * g[iy  ] * g[iz  ]
                                        + g[ix+1] * g[iy+1] * g[iz+1]
                                        + g[ix+2] * g[iy+2] * g[iz+2]
                                        + g[ix+3] * g[iy+3] * g[iz+3]
                                        + g[ix+4] * g[iy+4] * g[iz+4]
                                        + g[ix+5] * g[iy+5] * g[iz+5]
                                        + g[ix+6] * g[iy+6] * g[iz+6];
                        }
                        break;
                default:
                        for (n = 0; n < nf; n++) {
                                ix = idx[n];
                                iy = idy[n];
                                iz = idz[n];
                                gout[n] = 0;
                                for (i = 0; i < ng[RYS_ROOTS]; i++)
                                        gout[n] += g[ix+i] * g[iy+i] * g[iz+i];
                        }
                        break;
        }
}

int cint2e_sph(double *opkijl, const int *shls,
               const int *atm, const int natm,
               const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int k_sh = shls[2];
        const int l_sh = shls[3];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        const int k_l = bas(ANG_OF, k_sh);
        const int l_l = bas(ANG_OF, l_sh);
        int ng[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

        // maximum l supported: l=4
        ng[3] = j_l     + 1;     // 0:j_l
        ng[2] = l_l     + 1;     // 0:l_l
        ng[1] = k_l     + ng[2]; // 0:mmax, mmax = lk + ll, 0:k_l
        ng[0] = i_l     + ng[3]; // 0:nmax, nmax = li + lj, 0:i_l

        return cint2e_drv(opkijl, ng, 1, &gout2e,
                          &c2s_sph_2e1, NULL,
                          shls, atm, natm, bas, nbas, env);
}

int cint2e_cart(double *opkijl, const int *shls,
                const int *atm, const int natm,
                const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int k_sh = shls[2];
        const int l_sh = shls[3];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        const int k_l = bas(ANG_OF, k_sh);
        const int l_l = bas(ANG_OF, l_sh);
        int ng[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

        // maximum l supported: l=4
        ng[3] = j_l     + 1;     // 0:j_l
        ng[2] = l_l     + 1;     // 0:l_l
        ng[1] = k_l     + ng[2]; // 0:mmax, mmax = lk + ll, 0:k_l
        ng[0] = i_l     + ng[3]; // 0:nmax, nmax = li + lj, 0:i_l

        return cint2e_drv(opkijl, ng, 1, &gout2e,
                          &c2s_cart_2e1, NULL,
                          shls, atm, natm, bas, nbas, env);
}


/*
 * spinor <ki|jl> = (ij|kl); i,j\in electron 1; k,l\in electron 2
 */
int cint2e(double *opkijl, const int *shls,
           const int *atm, const int natm,
           const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int k_sh = shls[2];
        const int l_sh = shls[3];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        const int k_l = bas(ANG_OF, k_sh);
        const int l_l = bas(ANG_OF, l_sh);
        int ng[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

        // maximum l supported: l=4
        ng[3] = j_l     + 1;     // 0:j_l
        ng[2] = l_l     + 1;     // 0:l_l
        ng[1] = k_l     + ng[2]; // 0:mmax, mmax = lk + ll, 0:k_l
        ng[0] = i_l     + ng[3]; // 0:nmax, nmax = li + lj, 0:i_l

        return cint2e_drv(opkijl, ng, 1, &gout2e,
                          &c2s_sf_2e1, &c2s_sf_2e2,
                          shls, atm, natm, bas, nbas, env);
}


/*
 * * * * * * * * * * * * * * * * * * * * *
 * c to fortran interface
 */

C2F_(cint2e_cart)
C2F_(cint2e_sph)
C2F_(cint2e)
