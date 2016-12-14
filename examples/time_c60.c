/*
 * C60 molecule
 */
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "cint.h"

void run_all(int *atm, int natm, int *bas, int nbas, double *env);

int main()
{
        int natm = 60;
        int nbas = natm*20;
        // ATM_SLOTS = 6; BAS_SLOTS = 8;
        int *atm = malloc(sizeof(int) * natm * ATM_SLOTS);
        int *bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
        double *env = malloc(sizeof(double) * 10000);

        int i, n, off;
        off = PTR_ENV_START; // = 20

        atm(CHARGE_OF, 0)=6; atm(PTR_COORD, 0)=off; env[off+0]=-0.779; env[off+1]= 0.757; env[off+2]= 0.000; off+=3;
        atm(CHARGE_OF, 1)=6; atm(PTR_COORD, 1)=off; env[off+0]= 0.684; env[off+1]= 0.757; env[off+2]= 0.000; off+=3;
        atm(CHARGE_OF, 2)=6; atm(PTR_COORD, 2)=off; env[off+0]= 1.376; env[off+1]= 1.956; env[off+2]= 0.000; off+=3;
        atm(CHARGE_OF, 3)=6; atm(PTR_COORD, 3)=off; env[off+0]= 0.644; env[off+1]= 3.224; env[off+2]=-0.000; off+=3;
        atm(CHARGE_OF, 4)=6; atm(PTR_COORD, 4)=off; env[off+0]=-0.739; env[off+1]= 3.224; env[off+2]= 0.000; off+=3;
        atm(CHARGE_OF, 5)=6; atm(PTR_COORD, 5)=off; env[off+0]=-1.471; env[off+1]= 1.956; env[off+2]= 0.000; off+=3;
        atm(CHARGE_OF, 6)=6; atm(PTR_COORD, 6)=off; env[off+0]=-1.231; env[off+1]=-0.348; env[off+2]=-0.845; off+=3;
        atm(CHARGE_OF, 7)=6; atm(PTR_COORD, 7)=off; env[off+0]=-0.047; env[off+1]=-1.032; env[off+2]=-1.367; off+=3;
        atm(CHARGE_OF, 8)=6; atm(PTR_COORD, 8)=off; env[off+0]= 1.137; env[off+1]=-0.348; env[off+2]=-0.845; off+=3;
        atm(CHARGE_OF, 9)=6; atm(PTR_COORD, 9)=off; env[off+0]= 2.257; env[off+1]=-0.196; env[off+2]=-1.644; off+=3;
        atm(CHARGE_OF,10)=6; atm(PTR_COORD,10)=off; env[off+0]= 2.561; env[off+1]= 2.118; env[off+2]=-0.845; off+=3;
        atm(CHARGE_OF,11)=6; atm(PTR_COORD,11)=off; env[off+0]= 1.376; env[off+1]= 4.169; env[off+2]=-0.845; off+=3;
        atm(CHARGE_OF,12)=6; atm(PTR_COORD,12)=off; env[off+0]= 0.684; env[off+1]= 5.063; env[off+2]=-1.644; off+=3;
        atm(CHARGE_OF,13)=6; atm(PTR_COORD,13)=off; env[off+0]=-0.779; env[off+1]= 5.063; env[off+2]=-1.644; off+=3;
        atm(CHARGE_OF,14)=6; atm(PTR_COORD,14)=off; env[off+0]=-1.471; env[off+1]= 4.169; env[off+2]=-0.845; off+=3;
        atm(CHARGE_OF,15)=6; atm(PTR_COORD,15)=off; env[off+0]=-2.656; env[off+1]= 3.485; env[off+2]=-1.367; off+=3;
        atm(CHARGE_OF,16)=6; atm(PTR_COORD,16)=off; env[off+0]=-2.655; env[off+1]= 2.118; env[off+2]=-0.845; off+=3;
        atm(CHARGE_OF,17)=6; atm(PTR_COORD,17)=off; env[off+0]=-3.083; env[off+1]= 1.071; env[off+2]=-1.644; off+=3;
        atm(CHARGE_OF,18)=6; atm(PTR_COORD,18)=off; env[off+0]=-2.351; env[off+1]=-0.196; env[off+2]=-1.644; off+=3;
        atm(CHARGE_OF,19)=6; atm(PTR_COORD,19)=off; env[off+0]=-0.047; env[off+1]=-1.526; env[off+2]=-2.661; off+=3;
        atm(CHARGE_OF,20)=6; atm(PTR_COORD,20)=off; env[off+0]=-1.231; env[off+1]=-1.365; env[off+2]=-3.506; off+=3;
        atm(CHARGE_OF,21)=6; atm(PTR_COORD,21)=off; env[off+0]=-2.351; env[off+1]=-0.718; env[off+2]=-3.012; off+=3;
        atm(CHARGE_OF,22)=6; atm(PTR_COORD,22)=off; env[off+0]=-3.083; env[off+1]= 0.226; env[off+2]=-3.857; off+=3;
        atm(CHARGE_OF,23)=6; atm(PTR_COORD,23)=off; env[off+0]=-3.536; env[off+1]= 1.332; env[off+2]=-3.012; off+=3;
        atm(CHARGE_OF,24)=6; atm(PTR_COORD,24)=off; env[off+0]=-3.536; env[off+1]= 2.626; env[off+2]=-3.506; off+=3;
        atm(CHARGE_OF,25)=6; atm(PTR_COORD,25)=off; env[off+0]=-3.083; env[off+1]= 3.732; env[off+2]=-2.660; off+=3;
        atm(CHARGE_OF,26)=6; atm(PTR_COORD,26)=off; env[off+0]=-2.352; env[off+1]= 4.677; env[off+2]=-3.506; off+=3;
        atm(CHARGE_OF,27)=6; atm(PTR_COORD,27)=off; env[off+0]=-1.231; env[off+1]= 5.324; env[off+2]=-3.012; off+=3;
        atm(CHARGE_OF,28)=6; atm(PTR_COORD,28)=off; env[off+0]= 2.561; env[off+1]= 3.485; env[off+2]=-1.367; off+=3;
        atm(CHARGE_OF,29)=6; atm(PTR_COORD,29)=off; env[off+0]= 1.136; env[off+1]= 4.307; env[off+2]=-5.673; off+=3;
        atm(CHARGE_OF,30)=6; atm(PTR_COORD,30)=off; env[off+0]= 2.256; env[off+1]= 4.155; env[off+2]=-4.873; off+=3;
        atm(CHARGE_OF,31)=6; atm(PTR_COORD,31)=off; env[off+0]= 2.988; env[off+1]= 2.887; env[off+2]=-4.873; off+=3;
        atm(CHARGE_OF,32)=6; atm(PTR_COORD,32)=off; env[off+0]= 2.561; env[off+1]= 1.841; env[off+2]=-5.673; off+=3;
        atm(CHARGE_OF,33)=6; atm(PTR_COORD,33)=off; env[off+0]= 1.376; env[off+1]= 2.002; env[off+2]=-6.518; off+=3;
        atm(CHARGE_OF,34)=6; atm(PTR_COORD,34)=off; env[off+0]=-0.779; env[off+1]= 3.201; env[off+2]=-6.518; off+=3;
        atm(CHARGE_OF,35)=6; atm(PTR_COORD,35)=off; env[off+0]=-1.231; env[off+1]= 4.307; env[off+2]=-5.673; off+=3;
        atm(CHARGE_OF,36)=6; atm(PTR_COORD,36)=off; env[off+0]=-0.047; env[off+1]= 4.991; env[off+2]=-5.150; off+=3;
        atm(CHARGE_OF,37)=6; atm(PTR_COORD,37)=off; env[off+0]=-0.047; env[off+1]= 5.485; env[off+2]=-3.857; off+=3;
        atm(CHARGE_OF,38)=6; atm(PTR_COORD,38)=off; env[off+0]= 1.136; env[off+1]= 5.324; env[off+2]=-3.012; off+=3;
        atm(CHARGE_OF,39)=6; atm(PTR_COORD,39)=off; env[off+0]= 2.257; env[off+1]= 4.677; env[off+2]=-3.506; off+=3;
        atm(CHARGE_OF,40)=6; atm(PTR_COORD,40)=off; env[off+0]= 3.441; env[off+1]= 2.626; env[off+2]=-3.506; off+=3;
        atm(CHARGE_OF,41)=6; atm(PTR_COORD,41)=off; env[off+0]= 3.441; env[off+1]= 1.332; env[off+2]=-3.012; off+=3;
        atm(CHARGE_OF,42)=6; atm(PTR_COORD,42)=off; env[off+0]= 2.989; env[off+1]= 0.226; env[off+2]=-3.857; off+=3;
        atm(CHARGE_OF,43)=6; atm(PTR_COORD,43)=off; env[off+0]= 2.561; env[off+1]= 0.473; env[off+2]=-5.150; off+=3;
        atm(CHARGE_OF,44)=6; atm(PTR_COORD,44)=off; env[off+0]= 1.376; env[off+1]=-0.210; env[off+2]=-5.673; off+=3;
        atm(CHARGE_OF,45)=6; atm(PTR_COORD,45)=off; env[off+0]= 0.644; env[off+1]= 0.734; env[off+2]=-6.518; off+=3;
        atm(CHARGE_OF,46)=6; atm(PTR_COORD,46)=off; env[off+0]=-0.739; env[off+1]= 0.734; env[off+2]=-6.518; off+=3;
        atm(CHARGE_OF,47)=6; atm(PTR_COORD,47)=off; env[off+0]=-1.471; env[off+1]= 2.002; env[off+2]=-6.518; off+=3;
        atm(CHARGE_OF,48)=6; atm(PTR_COORD,48)=off; env[off+0]=-2.352; env[off+1]= 4.155; env[off+2]=-4.873; off+=3;
        atm(CHARGE_OF,49)=6; atm(PTR_COORD,49)=off; env[off+0]=-3.084; env[off+1]= 2.887; env[off+2]=-4.873; off+=3;
        atm(CHARGE_OF,50)=6; atm(PTR_COORD,50)=off; env[off+0]=-2.656; env[off+1]= 1.840; env[off+2]=-5.673; off+=3;
        atm(CHARGE_OF,51)=6; atm(PTR_COORD,51)=off; env[off+0]=-2.656; env[off+1]= 0.473; env[off+2]=-5.150; off+=3;
        atm(CHARGE_OF,52)=6; atm(PTR_COORD,52)=off; env[off+0]=-1.471; env[off+1]=-0.210; env[off+2]=-5.673; off+=3;
        atm(CHARGE_OF,53)=6; atm(PTR_COORD,53)=off; env[off+0]=-0.779; env[off+1]=-1.104; env[off+2]=-4.873; off+=3;
        atm(CHARGE_OF,54)=6; atm(PTR_COORD,54)=off; env[off+0]= 0.684; env[off+1]=-1.104; env[off+2]=-4.873; off+=3;
        atm(CHARGE_OF,55)=6; atm(PTR_COORD,55)=off; env[off+0]= 1.136; env[off+1]=-1.365; env[off+2]=-3.506; off+=3;
        atm(CHARGE_OF,56)=6; atm(PTR_COORD,56)=off; env[off+0]= 2.257; env[off+1]=-0.718; env[off+2]=-3.012; off+=3;
        atm(CHARGE_OF,57)=6; atm(PTR_COORD,57)=off; env[off+0]= 2.988; env[off+1]= 3.732; env[off+2]=-2.661; off+=3;
        atm(CHARGE_OF,58)=6; atm(PTR_COORD,58)=off; env[off+0]= 2.989; env[off+1]= 1.071; env[off+2]=-1.644; off+=3;
        atm(CHARGE_OF,59)=6; atm(PTR_COORD,59)=off; env[off+0]= 0.684; env[off+1]= 3.201; env[off+2]=-6.518; off+=3;

        // cc-pVDZ
        env[off+ 0] = 6665.0; // s
        env[off+ 1] = 1000.0;
        env[off+ 2] = 228.00;
        env[off+ 3] = 64.710;
        env[off+ 4] = 21.060;
        env[off+ 5] = 7.4950;
        env[off+ 6] = 2.7970;
        env[off+ 7] = 0.5215;
        env[off+ 8] = 0.000692*CINTgto_norm(0,env[off+0]); env[off+16] =-0.000146*CINTgto_norm(0,env[off+0]);
        env[off+ 9] = 0.005329*CINTgto_norm(0,env[off+1]); env[off+17] =-0.001154*CINTgto_norm(0,env[off+1]);
        env[off+10] = 0.027077*CINTgto_norm(0,env[off+2]); env[off+18] =-0.005725*CINTgto_norm(0,env[off+2]);
        env[off+11] = 0.101718*CINTgto_norm(0,env[off+3]); env[off+19] =-0.023312*CINTgto_norm(0,env[off+3]);
        env[off+12] = 0.274740*CINTgto_norm(0,env[off+4]); env[off+20] =-0.063955*CINTgto_norm(0,env[off+4]);
        env[off+13] = 0.448564*CINTgto_norm(0,env[off+5]); env[off+21] =-0.149981*CINTgto_norm(0,env[off+5]);
        env[off+14] = 0.285074*CINTgto_norm(0,env[off+6]); env[off+22] =-0.127262*CINTgto_norm(0,env[off+6]);
        env[off+15] = 0.015204*CINTgto_norm(0,env[off+7]); env[off+23] = 0.544529*CINTgto_norm(0,env[off+7]);
        env[off+24] = 0.1596; // s
        env[off+25] = 1*CINTgto_norm(0,env[off+24]);
        env[off+26] = 9.4390; // p
        env[off+27] = 2.0020;
        env[off+28] = 0.5456;
        env[off+29] = 0.038109*CINTgto_norm(1,env[off+26]);
        env[off+30] = 0.209480*CINTgto_norm(1,env[off+27]);
        env[off+31] = 0.508557*CINTgto_norm(1,env[off+28]);
        env[off+32] = 0.1517; // p
        env[off+33] = 1*CINTgto_norm(1,env[off+32]);
        env[off+34] = 0.55; // d
        env[off+35] = 1*CINTgto_norm(2,env[off+34]);
        for (i = 0, n = 0; i < natm; i++) {
                bas[ATOM_OF  +BAS_SLOTS*n] = i;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 8;
                bas[NCTR_OF  +BAS_SLOTS*n] = 2;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+0;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+8;
                n++;

                bas[ATOM_OF  +BAS_SLOTS*n] = i;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+24;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+25;
                n++;

                bas[ATOM_OF  +BAS_SLOTS*n] = i;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+26;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+29;
                n++;

                bas[ATOM_OF  +BAS_SLOTS*n] = i;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+32;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+33;
                n++;

                bas[ATOM_OF  +BAS_SLOTS*n] = i;
                bas[ANG_OF   +BAS_SLOTS*n] = 2;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+34;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+35;
                n++;
        }
        nbas = n;

        run_all(atm, natm, bas, nbas, env);
        // 6478s on one core of 3.1G I5 CPU
        free(atm);
        free(bas);
        free(env);
}

void run_all(int *atm, int natm, int *bas, int nbas, double *env)
{
        int i, j, k, l, ij, kl;
        int di, dj, dk, dl;
        int kl_max;
        int shls[4];
        double *buf;
        int *ishls = malloc(sizeof(int)*nbas*nbas);
        int *jshls = malloc(sizeof(int)*nbas*nbas);
        for (i = 0, ij = 0; i < nbas; i++) {
                for (j = 0; j <= i; j++, ij++) {
                        ishls[ij] = i;
                        jshls[ij] = j;
                }
        }

        int ncgto = CINTtot_cgto_spheric(bas, nbas);
        printf("\tshells = %d, total cGTO = %d, total pGTO = %d\n",
               nbas, ncgto,
               CINTtot_pgto_spheric(bas, nbas));

        int pct;
        long count;
        double time0, time1 = 0;
        double tt, tot;
        tot = (double)ncgto*ncgto*ncgto*ncgto/8;
        time0 = omp_get_wtime();

        printf("\tcint2e_sph with optimizer: total num ERI = %.2e\n", tot);
        CINTOpt *opt = NULL;
        cint2e_sph_optimizer(&opt, atm, natm, bas, nbas, env);

        pct = 0; count = 0;
#pragma omp parallel default(none) \
        shared(atm, natm, bas, nbas, env, ishls, jshls, opt, time0, pct, count, stdout) \
        private(di, dj, dk, dl, i, j, k, l, ij, kl, kl_max, shls, buf, time1)
#pragma omp for nowait schedule(dynamic, 2)
        for (ij = 0; ij < nbas*(nbas+1)/2; ij++) {
                i = ishls[ij];
                j = jshls[ij];
                di = CINTcgto_spheric(i, bas);
                dj = CINTcgto_spheric(j, bas);
                // when ksh==ish, there exists k<i, so it's possible kl>ij
                kl_max = (i+1)*(i+2)/2;
                for (kl = 0; kl < kl_max; kl++) {
                        k = ishls[kl];
                        l = jshls[kl];
                        dk = CINTcgto_spheric(k, bas);
                        dl = CINTcgto_spheric(l, bas);
                        shls[0] = i;
                        shls[1] = j;
                        shls[2] = k;
                        shls[3] = l;
                        buf = malloc(sizeof(double) * di*dj*dk*dl);
                        cint2e_sph(buf, shls, atm, natm, bas, nbas, env, opt);
                        free(buf);
                }
                count += kl_max;
                if (100l*count/((long)nbas*nbas*(nbas+1)*(nbas+2)/8) > pct) {
                        pct++;
                        time1 = omp_get_wtime();
                        printf("\t%d%%, CPU time = %8.2f\r", pct, time1-time0);
                        fflush(stdout);
                }
        }
        time1 = omp_get_wtime();
        tt = time1-time0;
        printf("\t100%%, CPU time = %8.2f, %8.4f Mflops\n", tt, tot/1e6/tt);

        CINTdel_optimizer(&opt);
        free(ishls);
        free(jshls);
}

