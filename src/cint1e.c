/*
 * File: cint1e.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic cGTO integrals
 */

#include <stdlib.h>
#include <math.h>
#include "cint_bas.h"
#include "g1e.h"
#include "cint1e.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"


/*
 * 1e GTO integral basic loop for < i|j>, no 1/r
 */
int cint1e_loop(double *gctr, const int *ng, const double fac,
                void (*const f_gout)(),
                const int *shls, const int *atm, const int *bas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        const int i_ctr = bas(NCTR_OF, i_sh);
        const int j_ctr = bas(NCTR_OF, j_sh);
        const int nfi = len_cart(i_l);
        const int nfj = len_cart(j_l);
        const int n_comp = ng[POS_E1] * ng[TENSOR];
        const int nf = nfi * nfj;
        const double *ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        const double *rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        int ip, jp, n;
        int g_is_0;
        int *const idx = (int *)malloc(sizeof(int) * nf * 3);
        double aij, dij, eij, rrij;
        double *const g = (double *)malloc(sizeof(double) * ng[0] * ng[1] * 3 * ((1<<ng[GSHIFT])+1)); // +1 as buffer
        double *const gout = (double *)malloc(sizeof(double) * nf * n_comp);
        double *const gctri = (double *)malloc(sizeof(double) * nf * i_ctr * n_comp);

        g1e_index_xyz(idx, ng, shls, bas);

        rrij = square_dist(ri, rj);

        g_is_0 = 1;
        n = nf * i_ctr * j_ctr * n_comp;
        dset0(n, gctr);
        for (jp = 0; jp < bas(NPRIM_OF, j_sh); jp++) {
                n = nf * i_ctr * n_comp;
                dset0(n, gctri);
                for (ip = 0; ip < bas(NPRIM_OF, i_sh); ip++) {
                        aij = ai[ip] + aj[jp];
                        eij = (ai[ip] * aj[jp] / aij) * rrij;
                        if (eij > EXPCUTOFF)
                                continue;
                        g_is_0 = 0;

                        dij = SQRTPI * PI * exp(-eij) / (aij * sqrt(aij)) * fac;
                        g_ovlp(g, ng, ai[ip], aj[jp], ri, rj, dij);

                        dset0(nf * n_comp, gout);
                        (*f_gout)(g, ng, gout, nf, idx, ai[ip], aj[jp],
                                  shls, atm, bas, env);

                        n = nf * n_comp;
                        prim_to_ctr(gctri, n, gout, 1, i_sh, ip, bas, env);
                }
                n = nf * i_ctr;
                prim_to_ctr(gctr, n, gctri, n_comp, j_sh, jp, bas, env);
        }
        free(g);
        free(idx);
        free(gout);
        free(gctri);

        return g_is_0;
}


/*
 * for given nuclear model, calculate temporary parameter tau
 * aij = ai + aj
 */
static double nuc_mod(const double aij, const int nuc_id,
                      const int *atm, const double *env)
{
        const int mass[109] = { \
                1  , 4  , 7  , 9  , 11 , 12 , 14 , 16 , 19  , 20 , \
                23 , 24 , 27 , 28 , 31 , 32 , 35 , 40 , 39  , 40 , \
                45 , 48 , 51 , 52 , 55 , 56 , 59 , 58 , 63  , 64 , \
                69 , 74 , 75 , 80 , 79 , 84 , 85 , 88 , 89  , 90 , \
                93 , 98 , 98 , 102, 103, 106, 107, 114, 115 , 120, \
                121, 130, 127, 132, 133, 138, 139, 140, 141 , 144, \
                145, 152, 153, 158, 159, 162, 162, 168, 169 , 174, \
                175, 180, 181, 184, 187, 192, 193, 195, 197 , 202, \
                205, 208, 209, 209, 210, 222, 223, 226, 227 , 232, \
                231, 238, 237, 244, 243, 247, 247, 251, 252 , 257, \
                258, 259, 262, 261, 262, 263, 262, 265, 266};

        /* screened nuclear potential of Gaussian nuclear model:
         * M. Filatov and D. Cremer, Theor. Chem. Acc. 108, 168 (2002)
         * M. Filatov and D. Cremer, Chem. Phys. Lett. 351, 259 (2002)
        r = (-0.263188 * atm(CHARGE_OF, nuc_id) + 106.016974 \
             + 138.985999 / atm(CHARGE_OF, nuc_id)) \
             / (env[PTR_LIGHT_SPEED] * env[PTR_LIGHT_SPEED]);
        eta = 1 / (r * r); */
        double a = pow(mass[atm(CHARGE_OF, nuc_id)], (double)1 / 3);
        double r = (0.836 * a + 0.570) / 52917.7249;
        double eta = 1.5 / (r * r);

        switch (atm(NUC_MOD_OF, nuc_id)) {
                case POINT_NUC:
                        return 1;
                case GAUSSIAN_NUC:
                        return sqrt(eta / (aij + eta));
                default:
                        return 1;
        }
}
static double no_nuc_mod(const double aij, const int nuc_id,
                         const int *atm, const double *env)
{
        return 1;
}

/*
 * 1e GTO integral basic loop for < i|1/r|j>, no 1/r
 * if nuc_id >= 0: nuclear attraction, use nuclear model
 * if nuc_id <  0: 1/r potential, do not use nuclear model
 */
int cint1e_nuc_loop(double *gctr, const int *ng, const double fac,
                void (*const f_gout)(), const int nuc_id,
                const int *shls, const int *atm, const int *bas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        const int i_ctr = bas(NCTR_OF, i_sh);
        const int j_ctr = bas(NCTR_OF, j_sh);
        const int nfi = len_cart(i_l);
        const int nfj = len_cart(j_l);
        const int nf = nfi * nfj;
        const int n_comp = ng[POS_E1] * ng[TENSOR];
        const double *ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        const double *rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        const double *ai = env + bas(PTR_EXP, i_sh);
        const double *aj = env + bas(PTR_EXP, j_sh);
        int ip, jp, i, n;
        int g_is_0;
        double tau;
        const double *cr;
        double (*f_nuc_mod)();
        double x, u[MXRYSROOTS], w[MXRYSROOTS];
        int *const idx = (int *)malloc(sizeof(int) * nf * 3);
        double rij[3], aij, dij, eij, rrij, t2;
        double *const g = (double *)malloc(sizeof(double) * ng[0] * ng[1] * 3 * ((1<<ng[GSHIFT])+1)); // +1 as buffer
        double *const gout = (double *)malloc(sizeof(double) * nf * n_comp);
        double *const gctri = (double *)malloc(sizeof(double) * nf * i_ctr * n_comp);

        if (nuc_id < 0) {
                cr = &env[PTR_RINV_ORIG];
                f_nuc_mod = no_nuc_mod;
        } else {
                cr = &env[atm(PTR_COORD, nuc_id)],
                f_nuc_mod = nuc_mod;
        }

        g1e_index_xyz(idx, ng, shls, bas);

        rrij = square_dist(ri, rj);

        g_is_0 = 1;
        n = nf * i_ctr * j_ctr * n_comp;
        dset0(n, gctr);
        for (jp = 0; jp < bas(NPRIM_OF, j_sh); jp++) {
                n = nf * i_ctr * n_comp;
                dset0(n, gctri);
                for (ip = 0; ip < bas(NPRIM_OF, i_sh); ip++) {
                        aij = ai[ip] + aj[jp];
                        eij = (ai[ip] * aj[jp] / aij) * rrij;
                        if (eij > EXPCUTOFF)
                                continue;
                        g_is_0 = 0;

                        rij[0] = (ai[ip] * ri[0] + aj[jp] * rj[0]) / aij;
                        rij[1] = (ai[ip] * ri[1] + aj[jp] * rj[1]) / aij;
                        rij[2] = (ai[ip] * ri[2] + aj[jp] * rj[2]) / aij;
                        tau = (*f_nuc_mod)(aij, nuc_id, atm, env);
                        x = aij * square_dist(rij, cr) * tau * tau;
                        rys_roots(ng[RYS_ROOTS], x, u, w);

                        dij = 2 * M_PI * exp(-eij) / aij * fac;
                        dset0(nf * n_comp, gout);
                        for (i = 0; i < ng[RYS_ROOTS]; i++) {
                                t2 = u[i] / (1 + u[i]) * tau * tau;
                                g_nuc(g, ng, aij, rij, ri, rj,
                                      cr, t2, dij * w[i] * tau);

                                (*f_gout)(g, ng, gout, nf, idx, ai[ip], aj[jp],
                                          shls, atm, bas, env);
                        }

                        n = nf * n_comp;
                        prim_to_ctr(gctri, n, gout, 1, i_sh, ip, bas, env);
                }
                n = nf * i_ctr;
                prim_to_ctr(gctr, n, gctri, n_comp, j_sh, jp, bas, env);
        }
        free(g);
        free(idx);
        free(gout);
        free(gctri);

        return g_is_0;
}


/*
 * 1e integrals <i|O|j> without 1/r
 */
int cint1e_drv(double *opij, int *ng, const double fac,
               void (*const f_gout)(), void (*const f_c2s)(),
               const int *shls, const int *atm, const int natm,
               const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        const int i_ctr = bas(NCTR_OF, i_sh);
        const int j_ctr = bas(NCTR_OF, j_sh);
        const int nfi = len_cart(i_l);
        const int nfj = len_cart(j_l);
        const int nc = nfi * nfj * i_ctr * j_ctr * ng[POS_E1];
        int ip, jp, nop;
        int n;
        int g_is_0;
        double *const gctr = (double *)malloc(sizeof(double) * nc * ng[TENSOR]);
        double *pgctr = gctr;

        g_is_0 = cint1e_loop(gctr, ng, fac, f_gout,
                             shls, atm, bas, env);

        if (f_c2s == c2s_sph_1e) {
                ip = cgtos_spheric(i_sh, bas);
                jp = cgtos_spheric(j_sh, bas);
                nop = ip * jp;
        } else if (f_c2s == c2s_cart_1e) {
                ip = cgtos_cart(i_sh, bas);
                jp = cgtos_cart(j_sh, bas);
                nop = ip * jp;
        } else {
                ip = cgtos_spinor(i_sh, bas);
                jp = cgtos_spinor(j_sh, bas);
                nop = ip * jp * OF_CMPLX;
        }

        if (g_is_0) {
                dset0(nop * ng[TENSOR], opij);
        } else {
                for (n = 0; n < ng[TENSOR]; n++) {
                        (*f_c2s)(opij, pgctr, shls, bas);
                        opij += nop;
                        pgctr += nc;
                }
        }
        free(gctr);
        return !g_is_0;
}


/*
 * 1e integrals <i|O|j> with 1/r
 */
int cint1e_rinv_drv(double *opij, int *ng, const double fac,
                    void (*const f_gout)(), void (*const f_c2s)(),
                    const int *shls, const int *atm, const int natm,
                    const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        const int i_ctr = bas(NCTR_OF, i_sh);
        const int j_ctr = bas(NCTR_OF, j_sh);
        const int nfi = len_cart(i_l);
        const int nfj = len_cart(j_l);
        const int nc = nfi * nfj * i_ctr * j_ctr * ng[POS_E1];
        int ip, jp, nop;
        int n;
        int g_is_0;
        double *const gctr = (double *)malloc(sizeof(double) * nc * ng[TENSOR]);
        double *pgctr = gctr;

        ng[RYS_ROOTS] = (ng[0] + 1) / 2; // li + lj + 2
        g_is_0 = cint1e_nuc_loop(gctr, ng, fac, f_gout, -1,
                                 shls, atm, bas, env);

        if (f_c2s == c2s_sph_1e) {
                ip = cgtos_spheric(i_sh, bas);
                jp = cgtos_spheric(j_sh, bas);
                nop = ip * jp;
        } else if (f_c2s == c2s_cart_1e) {
                ip = cgtos_cart(i_sh, bas);
                jp = cgtos_cart(j_sh, bas);
                nop = ip * jp;
        } else {
                ip = cgtos_spinor(i_sh, bas);
                jp = cgtos_spinor(j_sh, bas);
                nop = ip * jp * OF_CMPLX;
        }

        if (g_is_0) {
                dset0(nop * ng[TENSOR], opij);
        } else {
                for (n = 0; n < ng[TENSOR]; n++) {
                        (*f_c2s)(opij, pgctr, shls, bas);
                        opij += nop;
                        pgctr += nc;
                }
        }
        free(gctr);
        return !g_is_0;
}


/*
 * 1e integrals <i|O|j> with nuclear attraction
 * TODO: add the gaussian nuclear model
 */
int cint1e_nuc_drv(double *opij, int *ng, const double fac,
                   void (*const f_gout)(), void (*const f_c2s)(),
                   const int *shls, const int *atm, const int natm,
                   const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        const int i_ctr = bas(NCTR_OF, i_sh);
        const int j_ctr = bas(NCTR_OF, j_sh);
        const int nfi = len_cart(i_l);
        const int nfj = len_cart(j_l);
        const int nc = nfi * nfj * i_ctr * j_ctr * ng[POS_E1];
        int ip, jp, nop;
        int i, n;
        int g_is_0, is_0;
        double *const gctr = (double *)malloc(sizeof(double) * nc * ng[TENSOR]);
        double *const gctr0 = (double *)malloc(sizeof(double) * nc * ng[TENSOR]);
        double *pgctr = gctr;

        ng[RYS_ROOTS] = (ng[0] + 1) / 2; // li + lj + 2
        dset0(nc * ng[TENSOR], gctr);
        g_is_0 = 1;
        for (n = 0; n < natm; n++) {
                is_0 = cint1e_nuc_loop(gctr0, ng, fac, f_gout, n,
                                       shls, atm, bas, env);
                if (!is_0) {
                        g_is_0 = 0;
                        for (i = 0; i < nc * ng[TENSOR]; i++) {
                                gctr[i] += -abs(atm(CHARGE_OF,n)) * gctr0[i];
                        }
                }
        }

        if (f_c2s == c2s_sph_1e) {
                ip = cgtos_spheric(i_sh, bas);
                jp = cgtos_spheric(j_sh, bas);
                nop = ip * jp;
        } else if (f_c2s == c2s_cart_1e) {
                ip = cgtos_cart(i_sh, bas);
                jp = cgtos_cart(j_sh, bas);
                nop = ip * jp;
        } else {
                ip = cgtos_spinor(i_sh, bas);
                jp = cgtos_spinor(j_sh, bas);
                nop = ip * jp * OF_CMPLX;
        }

        if (g_is_0) {
                dset0(nop * ng[TENSOR], opij);
        } else {
                for (n = 0; n < ng[TENSOR]; n++) {
                        (*f_c2s)(opij, pgctr, shls, bas);
                        opij += nop;
                        pgctr += nc;
                }
        }
        free(gctr0);
        free(gctr);
        return !g_is_0;
}


/* template
 *
void gout1e(double *g, const int *ng,
        double *gout, const int nf, const int *idx,
        const double ai, const double aj,
        const int *shls, const int *atm, const int natm,
        const int *bas, const int nbas, const double *env)
{
        int ix, iy, iz, n;
        const int *idy = idx + nf;
        const int *idz = idx + nf * 2;

        for (n = 0; n < nf; n++) {
                ix = idx[n];
                iy = idy[n];
                iz = idz[n];
                gout[n] += g[ix] * g[iy] * g[iz];
        }
}

int cint1e_ovlp_sph(double *opij, const int *shls,
                    const int *atm, const int natm,
                    const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        int ng[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

        ng[1] = j_l     + 1;     // 0:j_l,
        ng[0] = i_l     + ng[1]; // 0:nmax, nmax = li + lj, 0:i_l

        return cint1e_drv(opij, ng, 1, &gout1e, &c2s_sph_1e,
                          shls, atm, natm, bas, nbas, env);
}

int cint1e_ovlp(double *opij, const int *shls,
                const int *atm, const int natm,
                const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        int ng[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

        ng[1] = j_l     + 1;     // 0:j_l,
        ng[0] = i_l     + ng[1]; // 0:nmax, nmax = li + lj, 0:i_l

        return cint1e_drv(opij, ng, 1, &gout1e, &c2s_sf_1e,
                          shls, atm, natm, bas, nbas, env);
}

int cint1e_nuc_sph(double *opij, const int *shls,
                   const int *atm, const int natm,
                   const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        int ng[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

        ng[1] = j_l     + 1;     // 0:j_l,
        ng[0] = i_l     + ng[1]; // 0:nmax, nmax = li + lj, 0:i_l

        return cint1e_nuc_drv(opij, ng, 1, gout1e, c2s_sph_1e,
                              shls, atm, natm, bas, nbas, env);
}

int cint1e_nuc(double *opij, const int *shls,
               const int *atm, const int natm,
               const int *bas, const int nbas, const double *env)
{
        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int i_l = bas(ANG_OF, i_sh);
        const int j_l = bas(ANG_OF, j_sh);
        int ng[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

        ng[1] = j_l     + 1;     // 0:j_l,
        ng[0] = i_l     + ng[1]; // 0:nmax, nmax = li + lj, 0:i_l

        return cint1e_nuc_drv(opij, ng, 1, gout1e, c2s_sf_1e,
                              shls, atm, natm, bas, nbas, env);
}

// cint1e to fortran interface

C2F_(cint1e_ovlp_sph)
C2F_(cint1e_nuc_sph)
C2F_(cint1e_ovlp)
C2F_(cint1e_nuc)

template */
