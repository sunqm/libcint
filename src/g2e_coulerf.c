/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Attenuated coulomb operator erf(-w r_{12}) / r_{12}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "config.h"
#include "cint_bas.h"
#include "rys_roots.h"
#include "g2e.h"

FINT CINTg0_2e_coulerf(double *g, double *rij, double *rkl, double cutoff, CINTEnvVars *envs);

void CINTinit_int2e_coulerf_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                                    FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_int2e_EnvVars(envs, ng, shls, atm, natm, bas, nbas, env);
        envs->f_g0_2e = &CINTg0_2e_coulerf;
}

/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
FINT CINTg0_2e_coulerf(double *g, double *rij, double *rkl, double cutoff, CINTEnvVars *envs)
{
        double aij = envs->ai[0] + envs->aj[0];
        double akl = envs->ak[0] + envs->al[0];
        double omega = envs->env[PTR_RANGE_OMEGA];
        double a0, a1, fac1, x;
        double u[MXRYSROOTS];
        double *w = g + envs->g_size * 2; // ~ gz
        double rijrkl[3];
        rijrkl[0] = rij[0] - rkl[0];
        rijrkl[1] = rij[1] - rkl[1];
        rijrkl[2] = rij[2] - rkl[2];

        a1 = aij * akl;
        a0 = a1 / (aij + akl);

        double theta = 0;
        if (omega > 0) {
// For long-range part of range-separated Coulomb operator
                theta = omega * omega / (omega * omega + a0);
                a0 *= theta;
        }

        fac1 = sqrt(a0 / (a1 * a1 * a1)) * envs->fac[0];
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);
        CINTrys_roots(envs->nrys_roots, x, u, w);
        if (envs->g_size == 1) {
                g[0] = 1;
                g[1] = 1;
                g[2] *= fac1;
                return 1;
        }

        FINT irys;
        if (omega > 0) {
                /* u[:] = tau^2 / (1 - tau^2)
                 * transform u[:] to theta^-1 tau^2 / (theta^-1 - tau^2)
                 * so the rest code can be reused.
                 */
                for (irys = 0; irys < envs->nrys_roots; irys++) {
                        u[irys] /= u[irys] + 1 - u[irys] * theta;
                }
        }

        double u2, div, tmp1, tmp2, tmp3, tmp4;
        double rijrx[3];
        double rklrx[3];
        rijrx[0] = rij[0] - envs->rx_in_rijrx[0];
        rijrx[1] = rij[1] - envs->rx_in_rijrx[1];
        rijrx[2] = rij[2] - envs->rx_in_rijrx[2];
        rklrx[0] = rkl[0] - envs->rx_in_rklrx[0];
        rklrx[1] = rkl[1] - envs->rx_in_rklrx[1];
        rklrx[2] = rkl[2] - envs->rx_in_rklrx[2];
        struct _BC bc;
        double *c00 = bc.c00;
        double *c0p = bc.c0p;

        for (irys = 0; irys < envs->nrys_roots; irys++, c00+=3, c0p+=3) {
                /*
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
        return 1;
}
