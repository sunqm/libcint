/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Attenuated coulomb operator erf(-w r_{12}^2) / r_{12}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cint_const.h"
#include "cint_bas.h"
#include "rys_roots.h"
#include "g2e.h"

void CINTg0_2e_coul_gtg(double *g, double fac, CINTEnvVars *envs);

void CINTinit_int2e_coul_gtg_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                                     FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_int2e_EnvVars(envs, ng, shls, atm, natm, bas, nbas, env);
        envs->f_g0_2e = &CINTg0_2e_coul_gtg;
}

/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
void CINTg0_2e_coul_gtg(double *g, double fac, CINTEnvVars *envs)
{
        double aij = envs->aij;
        double akl = envs->akl;
        double omega = envs->env[PTR_RANGE_OMEGA];
        double a0, a1, fac1, x;
        double u[MXRYSROOTS];
        double *w = g + envs->g_size * 2; // ~ gz
        double rijrkl[3];
        rijrkl[0] = envs->rij[0] - envs->rkl[0];
        rijrkl[1] = envs->rij[1] - envs->rkl[1];
        rijrkl[2] = envs->rij[2] - envs->rkl[2];

        a1 = aij * akl;
        a0 = a1 / (aij + akl);

        fac1 = sqrt(a0 / (a1 * a1 * a1)) * fac;
        x = a0 *(rijrkl[0] * rijrkl[0]
               + rijrkl[1] * rijrkl[1]
               + rijrkl[2] * rijrkl[2]);
        CINTrys_roots(envs->nrys_roots, x, u, w);
        if (envs->g_size == 1) {
                g[0] = 1;
                g[1] = 1;
                g[2] *= fac1;
                return;
        }

        int irys;
        double u2, div, tmp1, tmp2, tmp3, tmp4;
        double *rijrx = envs->rijrx;
        double *rklrx = envs->rklrx;
        struct _BC bc;
        double *c00 = bc.c00;
        double *c0p = bc.c0p;
        double zeta = envs->env[PTR_GTG_ZETA];

        for (irys = 0; irys < envs->nrys_roots; irys++, c00+=3, c0p+=3) {
                /*
                 *u(irys) = t2/(1-t2)
                 *t2 = u(irys)/(1+u(irys))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = a0 * u[irys];
                if (zeta > 0) {
                        u2 += zeta;
                }
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
}
