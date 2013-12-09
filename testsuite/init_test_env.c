#include <stdlib.h>
#include <stdio.h>
#include "cint_bas.h"
#include "fblas.h"

// return last offset for env
int init_test_env(int *atm, int *natm, int *bas, int *nbas, double *env)
{
        *natm = 4;
        int i, n, nh, off, off0;

        dset0(PTR_ENV_START, env);
        off = PTR_ENV_START;
        for (i = 0; i < *natm; i++) {
                atm(CHARGE_OF,i) = (i+1) * 2;
                atm(PTR_COORD,i) = off;
                env[off+0] = .2 * (i+1);
                env[off+1] = .3 + (i+1) * .5;
                env[off+2] = .1 - (i+1) * .5;
                off += 3;
        }
        off0 = off;

        // basis with kappa > 0
        nh = 0;

        bas(ATOM_OF ,nh)  = 0;
        bas(ANG_OF  ,nh)  = 1;
        bas(KAPPA_OF,nh)  = 1;
        bas(NPRIM_OF,nh)  = 1;
        bas(NCTR_OF ,nh)  = 1;
        bas(PTR_EXP,nh)   = off;
        env[off+0] = 1;
        bas(PTR_COEFF,nh) = off + 1;
        env[off+1] = 1;
        off += 2;
        nh++;

        bas(ATOM_OF ,nh)  = 1;
        bas(ANG_OF  ,nh)  = 2;
        bas(KAPPA_OF,nh)  = 2;
        bas(NPRIM_OF,nh)  = 2;
        bas(NCTR_OF ,nh)  = 2;
        bas(PTR_EXP,nh)   = off;
        env[off+0] = 5;
        env[off+1] = 3;
        bas(PTR_COEFF,nh) = off + 2;
        env[off+2] = 1;
        env[off+3] = 2;
        env[off+4] = 4;
        env[off+5] = 1;
        off += 6;
        nh++;

        bas(ATOM_OF ,nh)  = 2;
        bas(ANG_OF  ,nh)  = 3;
        bas(KAPPA_OF,nh)  = 3;
        bas(NPRIM_OF,nh)  = 1;
        bas(NCTR_OF ,nh)  = 1;
        bas(PTR_EXP,nh)   = off;
        env[off+0] = 1;
        bas(PTR_COEFF,nh) = off + 1;
        env[off+1] = 1;
        off += 2;
        nh++;

        bas(ATOM_OF ,nh)  = 3;
        bas(ANG_OF  ,nh)  = 4;
        bas(KAPPA_OF,nh)  = 4;
        bas(NPRIM_OF,nh)  = 1;
        bas(NCTR_OF ,nh)  = 1;
        bas(PTR_EXP, nh)  = off;
        env[off+0] = .5;
        bas(PTR_COEFF,nh) = off + 1;
        env[off+1] = 1.;
        off = off + 2;
        nh++;

        *nbas = nh;

        // basis with kappa < 0
        n = off - off0;
        for (i = 0; i < n; i++) {
                env[off+i] = env[off0+i];
        }

        for (i = 0; i < nh; i++) {
                bas(ATOM_OF ,i+nh) = bas(ATOM_OF ,i);
                bas(ANG_OF  ,i+nh) = bas(ANG_OF  ,i) - 1;
                bas(KAPPA_OF,i+nh) =-bas(KAPPA_OF,i);
                bas(NPRIM_OF,i+nh) = bas(NPRIM_OF,i);
                bas(NCTR_OF ,i+nh) = bas(NCTR_OF ,i);
                bas(PTR_EXP ,i+nh) = bas(PTR_EXP ,i)  + n;
                bas(PTR_COEFF,i+nh)= bas(PTR_COEFF,i) + n;
                env[bas(PTR_COEFF,i+nh)] /= 2 * env[bas(PTR_EXP,i)];
        }
        env[bas(PTR_COEFF,5)+0] = env[bas(PTR_COEFF,1)+0] / (2 * env[bas(PTR_EXP,1)+0]);
        env[bas(PTR_COEFF,5)+1] = env[bas(PTR_COEFF,1)+1] / (2 * env[bas(PTR_EXP,1)+1]);
        env[bas(PTR_COEFF,5)+2] = env[bas(PTR_COEFF,1)+2] / (2 * env[bas(PTR_EXP,1)+0]);
        env[bas(PTR_COEFF,5)+3] = env[bas(PTR_COEFF,1)+3] / (2 * env[bas(PTR_EXP,1)+1]);

        return off + n;
}
