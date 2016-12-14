/*
 * ethane molecule
 */
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "cint.h"



void run_all(int *atm, int natm, int *bas, int nbas, double *env);

int cint2e_ip1_sph(double *buf, int *shls,
                   int *atm, int natm, int *bas, int nbas, double *env,
                   CINTOpt *opt);
void cint2e_ip1_sph_optimizer(CINTOpt **opt, int *atm, int natm,
                              int *bas, int nbas, double *env);

int main()
{
        int natm = 8;
        int nbas = natm*20;
        // ATM_SLOTS = 6; BAS_SLOTS = 8;
        int *atm = malloc(sizeof(int) * natm * ATM_SLOTS);
        int *bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
        double *env = malloc(sizeof(double) * 10000);

        int i, j, ia, n, off;
        off = PTR_ENV_START; // = 20

        atm(CHARGE_OF,0)=6; atm(PTR_COORD,0)=off; env[off+0]= 0.000; env[off+1]= 0.000; env[off+2]= 0.769; off+=3;
        atm(CHARGE_OF,1)=1; atm(PTR_COORD,1)=off; env[off+0]= 0.000; env[off+1]= 1.014; env[off+2]= 1.174; off+=3;
        atm(CHARGE_OF,2)=1; atm(PTR_COORD,2)=off; env[off+0]=-0.878; env[off+1]=-0.507; env[off+2]= 1.174; off+=3;
        atm(CHARGE_OF,3)=1; atm(PTR_COORD,3)=off; env[off+0]= 0.878; env[off+1]=-0.507; env[off+2]= 1.174; off+=3;
        atm(CHARGE_OF,4)=6; atm(PTR_COORD,4)=off; env[off+0]= 0.000; env[off+1]= 0.000; env[off+2]=-0.769; off+=3;
        atm(CHARGE_OF,5)=1; atm(PTR_COORD,5)=off; env[off+0]= 0.000; env[off+1]= 1.014; env[off+2]=-1.174; off+=3;
        atm(CHARGE_OF,6)=1; atm(PTR_COORD,6)=off; env[off+0]=-0.878; env[off+1]=-0.507; env[off+2]=-1.174; off+=3;
        atm(CHARGE_OF,7)=1; atm(PTR_COORD,7)=off; env[off+0]= 0.878; env[off+1]=-0.507; env[off+2]=-1.174; off+=3;

        // 6-31G
        env[off+0 ] = 3047.5249; env[off+6 ] = 0.0018347*CINTgto_norm(0,env[off+0 ]);
        env[off+1 ] = 457.36951; env[off+7 ] = 0.0140373*CINTgto_norm(0,env[off+1 ]);
        env[off+2 ] = 103.94869; env[off+8 ] = 0.0688426*CINTgto_norm(0,env[off+2 ]);
        env[off+3 ] = 29.210155; env[off+9 ] = 0.2321844*CINTgto_norm(0,env[off+3 ]);
        env[off+4 ] = 9.2866630; env[off+10] = 0.4679413*CINTgto_norm(0,env[off+4 ]);
        env[off+5 ] = 3.1639270; env[off+11] = 0.3623120*CINTgto_norm(0,env[off+5 ]);
        env[off+12] = 7.8682724; env[off+15] =-0.1193324*CINTgto_norm(0,env[off+12]);
        env[off+13] = 1.8812885; env[off+16] =-0.1608542*CINTgto_norm(0,env[off+13]);
        env[off+14] = 0.5442493; env[off+17] = 1.1434564*CINTgto_norm(0,env[off+14]);
        env[off+18] = 0.1687144; env[off+19] = 1.0000000*CINTgto_norm(0,env[off+18]);
        env[off+20] = 7.8682724; env[off+23] = 0.0689991*CINTgto_norm(1,env[off+20]);
        env[off+21] = 1.8812885; env[off+24] = 0.3164240*CINTgto_norm(1,env[off+21]);
        env[off+22] = 0.5442493; env[off+25] = 0.7443083*CINTgto_norm(1,env[off+22]);
        env[off+26] = 0.1687144; env[off+27] = 1.0000000*CINTgto_norm(1,env[off+26]);
        env[off+28] = 18.731137; env[off+31] = 0.0334946*CINTgto_norm(0,env[off+28]);
        env[off+29] = 2.8253937; env[off+32] = 0.2347269*CINTgto_norm(0,env[off+29]);
        env[off+30] = 0.6401217; env[off+33] = 0.8137573*CINTgto_norm(0,env[off+30]);
        env[off+34] = 0.1612778; env[off+35] = 1.0000000*CINTgto_norm(0,env[off+34]);
        for (i = 0, ia = 0, n = 0; i < 2; i++) {
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 6;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+0;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+6;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+12;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+15;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+18;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+19;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+20;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+23;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+26;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+27;
                n++;
                ia++;
                for (j = 0; j < 3; j++) {
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+28;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+31;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+34;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+35;
                        n++;
                        ia++;
                }
        }
        nbas = n;
        printf("6-31G basis\n");
        run_all(atm, natm, bas, nbas, env);

        // 6-311G**
        env[off+ 0] = 4563.240; env[off+17] = 0.0019666*CINTgto_norm(0,env[off+ 0]);
        env[off+ 1] = 682.0240; env[off+18] = 0.0152306*CINTgto_norm(0,env[off+ 1]);
        env[off+ 2] = 154.9730; env[off+19] = 0.0761269*CINTgto_norm(0,env[off+ 2]);
        env[off+ 3] = 44.45530; env[off+20] = 0.2608010*CINTgto_norm(0,env[off+ 3]);
        env[off+ 4] = 13.02900; env[off+21] = 0.6164620*CINTgto_norm(0,env[off+ 4]);
        env[off+ 5] = 1.827730; env[off+22] = 0.2210060*CINTgto_norm(0,env[off+ 5]);
        env[off+ 6] = 20.96420; env[off+23] = 0.1146600*CINTgto_norm(0,env[off+ 6]);
        env[off+ 7] = 4.803310; env[off+24] = 0.9199990*CINTgto_norm(0,env[off+ 7]);
        env[off+ 8] = 1.459330; env[off+25] = -0.003030*CINTgto_norm(0,env[off+ 8]);
        env[off+ 9] = 0.483456; env[off+26] = 1.0000000*CINTgto_norm(0,env[off+ 9]);
        env[off+10] = 0.145585; env[off+27] = 1.0000000*CINTgto_norm(0,env[off+10]);
        env[off+11] = 20.96420; env[off+28] = 0.0402487*CINTgto_norm(1,env[off+11]);
        env[off+12] = 4.803310; env[off+29] = 0.2375940*CINTgto_norm(1,env[off+12]);
        env[off+13] = 1.459330; env[off+30] = 0.8158540*CINTgto_norm(1,env[off+13]);
        env[off+14] = 0.483456; env[off+31] = 1.0000000*CINTgto_norm(1,env[off+14]);
        env[off+15] = 0.145585; env[off+32] = 1.0000000*CINTgto_norm(1,env[off+15]);
        env[off+16] = 0.626000; env[off+33] = 1.0000000*CINTgto_norm(2,env[off+16]);
        env[off+34] = 33.86500; env[off+40] = 0.0254938*CINTgto_norm(0,env[off+34]);
        env[off+35] = 5.094790; env[off+41] = 0.1903730*CINTgto_norm(0,env[off+35]);
        env[off+36] = 1.158790; env[off+42] = 0.8521610*CINTgto_norm(0,env[off+36]);
        env[off+37] = 0.325840; env[off+43] = 1.0000000*CINTgto_norm(0,env[off+37]);
        env[off+38] = 0.102741; env[off+44] = 1.0000000*CINTgto_norm(0,env[off+38]);
        env[off+39] = 0.750000; env[off+45] = 1.0000000*CINTgto_norm(0,env[off+39]);
        for (i = 0, ia = 0, n = 0; i < 2; i++) {
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 6;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+ 0;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+17;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+ 6;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+23;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+ 9;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+26;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+10;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+27;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+11;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+28;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+14;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+31;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+15;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+32;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 2;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+16;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+33;
                n++;
                ia++;
                for (j = 0; j < 3; j++) {
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+34;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+40;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+37;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+43;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+38;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+44;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 1;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+39;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+45;
                        n++;
                        ia++;
                }
        }
        nbas = n;
        printf("6-311G(dp) basis\n");
        run_all(atm, natm, bas, nbas, env);

        // cc-pVDZ, C
        env[off+ 0] = 6665.0; env[off+ 8]=0.000692*CINTgto_norm(0,env[off+ 0]); env[off+16]=-0.000146*CINTgto_norm(0,env[off+0]);
        env[off+ 1] = 1000.0; env[off+ 9]=0.005329*CINTgto_norm(0,env[off+ 1]); env[off+17]=-0.001154*CINTgto_norm(0,env[off+1]);
        env[off+ 2] = 228.00; env[off+10]=0.027077*CINTgto_norm(0,env[off+ 2]); env[off+18]=-0.005725*CINTgto_norm(0,env[off+2]);
        env[off+ 3] = 64.710; env[off+11]=0.101718*CINTgto_norm(0,env[off+ 3]); env[off+19]=-0.023312*CINTgto_norm(0,env[off+3]);
        env[off+ 4] = 21.060; env[off+12]=0.274740*CINTgto_norm(0,env[off+ 4]); env[off+20]=-0.063955*CINTgto_norm(0,env[off+4]);
        env[off+ 5] = 7.4950; env[off+13]=0.448564*CINTgto_norm(0,env[off+ 5]); env[off+21]=-0.149981*CINTgto_norm(0,env[off+5]);
        env[off+ 6] = 2.7970; env[off+14]=0.285074*CINTgto_norm(0,env[off+ 6]); env[off+22]=-0.127262*CINTgto_norm(0,env[off+6]);
        env[off+ 7] = 0.5215; env[off+15]=0.015204*CINTgto_norm(0,env[off+ 7]); env[off+23]= 0.544529*CINTgto_norm(0,env[off+7]);
        env[off+24] = 0.1596; env[off+25]=1.000000*CINTgto_norm(0,env[off+24]);
        env[off+26] = 9.4390; env[off+29]=0.038109*CINTgto_norm(1,env[off+26]);
        env[off+27] = 2.0020; env[off+30]=0.209480*CINTgto_norm(1,env[off+27]);
        env[off+28] = 0.5456; env[off+31]=0.508557*CINTgto_norm(1,env[off+28]);
        env[off+32] = 0.1517; env[off+33]=1.000000*CINTgto_norm(1,env[off+32]);
        env[off+34] = 0.55  ; env[off+35]=1.000000*CINTgto_norm(2,env[off+34]);
        // H
        env[off+36] = 13.010; env[off+39]=0.019685*CINTgto_norm(0,env[off+36]);
        env[off+37] = 1.9620; env[off+40]=0.137977*CINTgto_norm(0,env[off+37]);
        env[off+38] = 0.4446; env[off+41]=0.478148*CINTgto_norm(0,env[off+38]);
        env[off+42] = 0.1220; env[off+43]=1       *CINTgto_norm(0,env[off+42]);
        env[off+44] = 0.7270; env[off+45]=1       *CINTgto_norm(0,env[off+44]);
        for (i = 0, ia = 0, n = 0; i < 2; i++) {
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 8;
                bas[NCTR_OF  +BAS_SLOTS*n] = 2;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+0;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+8;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+24;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+25;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+26;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+29;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+32;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+33;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 2;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+34;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+35;
                n++;
                ia++;
                for (j = 0; j < 3; j++) {
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+36;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+39;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+42;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+43;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 1;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+44;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+45;
                        n++;
                        ia++;
                }
        }
        nbas = n;
        printf("cc-pVDZ basis\n");
        run_all(atm, natm, bas, nbas, env);

        // cc-pVTZ
        env[off+ 0] = 8236.0; env[off+18]= 0.000531*CINTgto_norm(0,env[off+ 0]); env[off+26]=-0.000113*CINTgto_norm(0,env[off+ 0]);
        env[off+ 1] = 1235.0; env[off+19]= 0.004108*CINTgto_norm(0,env[off+ 1]); env[off+27]=-0.000878*CINTgto_norm(0,env[off+ 1]);
        env[off+ 2] = 280.80; env[off+20]= 0.021087*CINTgto_norm(0,env[off+ 2]); env[off+28]=-0.004540*CINTgto_norm(0,env[off+ 2]);
        env[off+ 3] = 79.270; env[off+21]= 0.081853*CINTgto_norm(0,env[off+ 3]); env[off+29]=-0.018133*CINTgto_norm(0,env[off+ 3]);
        env[off+ 4] = 25.590; env[off+22]= 0.234817*CINTgto_norm(0,env[off+ 4]); env[off+30]=-0.055760*CINTgto_norm(0,env[off+ 4]);
        env[off+ 5] = 8.9970; env[off+23]= 0.434401*CINTgto_norm(0,env[off+ 5]); env[off+31]=-0.126895*CINTgto_norm(0,env[off+ 5]);
        env[off+ 6] = 3.3190; env[off+24]= 0.346129*CINTgto_norm(0,env[off+ 6]); env[off+32]=-0.170352*CINTgto_norm(0,env[off+ 6]);
        env[off+ 7] = 0.3643; env[off+25]=-0.008983*CINTgto_norm(0,env[off+ 7]); env[off+33]= 0.598684*CINTgto_norm(0,env[off+ 7]);
        env[off+ 8] = 0.9059; env[off+34]= 1.000000*CINTgto_norm(0,env[off+ 8]);
        env[off+ 9] = 0.1285; env[off+35]= 1.000000*CINTgto_norm(0,env[off+ 9]);
        env[off+10] = 18.710; env[off+36]= 0.014031*CINTgto_norm(1,env[off+10]);
        env[off+11] = 4.1330; env[off+37]= 0.086866*CINTgto_norm(1,env[off+11]);
        env[off+12] = 1.2000; env[off+38]= 0.290216*CINTgto_norm(1,env[off+12]);
        env[off+13] = 0.3827; env[off+39]= 1.000000*CINTgto_norm(1,env[off+13]);
        env[off+14] = 0.1209; env[off+40]= 1.000000*CINTgto_norm(1,env[off+14]);
        env[off+15] = 1.0970; env[off+41]= 1.000000*CINTgto_norm(2,env[off+15]);
        env[off+16] = 0.3180; env[off+42]= 1.000000*CINTgto_norm(2,env[off+16]);
        env[off+17] = 0.7610; env[off+43]= 1.000000*CINTgto_norm(3,env[off+17]);
        env[off+44] = 33.870; env[off+52]= 0.006068*CINTgto_norm(0,env[off+44]);
        env[off+45] = 5.0950; env[off+53]= 0.045308*CINTgto_norm(0,env[off+45]);
        env[off+46] = 1.1590; env[off+54]= 0.202822*CINTgto_norm(0,env[off+46]);
        env[off+47] = 0.3258; env[off+55]= 1.000000*CINTgto_norm(0,env[off+47]);
        env[off+48] = 0.1027; env[off+56]= 1.000000*CINTgto_norm(0,env[off+48]);
        env[off+49] = 1.4070; env[off+57]= 1.000000*CINTgto_norm(1,env[off+49]);
        env[off+50] = 0.3880; env[off+58]= 1.000000*CINTgto_norm(1,env[off+50]);
        env[off+51] = 1.0570; env[off+59]= 1.000000*CINTgto_norm(2,env[off+51]);
        for (i = 0, ia = 0, n = 0; i < 2; i++) {
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 8;
                bas[NCTR_OF  +BAS_SLOTS*n] = 2;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+ 0;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+18;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+ 8;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+34;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+ 9;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+35;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+10;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+38;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+13;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+39;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+14;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+40;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 2;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+15;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+41;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 2;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+16;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+42;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 3;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+17;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+43;
                n++;
                ia++;
                for (j = 0; j < 3; j++) {
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+44;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+52;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+47;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+55;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+48;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+56;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 1;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+49;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+57;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 1;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+50;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+58;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 2;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+51;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+59;
                        n++;
                        ia++;
                }
        }
        nbas = n;
        printf("cc-pVTZ basis\n");
        run_all(atm, natm, bas, nbas, env);

        env[off+ 0] = 33980.; env[off+24]= 0.000091*CINTgto_norm(0,env[off+ 0]); env[off+33]= -0.000019*CINTgto_norm(0,env[off+0]);
        env[off+ 1] = 5089.0; env[off+25]= 0.000704*CINTgto_norm(0,env[off+ 1]); env[off+34]= -0.000151*CINTgto_norm(0,env[off+1]);
        env[off+ 2] = 1157.0; env[off+26]= 0.003693*CINTgto_norm(0,env[off+ 2]); env[off+35]= -0.000785*CINTgto_norm(0,env[off+2]);
        env[off+ 3] = 326.60; env[off+27]= 0.015360*CINTgto_norm(0,env[off+ 3]); env[off+36]= -0.003324*CINTgto_norm(0,env[off+3]);
        env[off+ 4] = 106.10; env[off+28]= 0.052929*CINTgto_norm(0,env[off+ 4]); env[off+37]= -0.011512*CINTgto_norm(0,env[off+4]);
        env[off+ 5] = 38.110; env[off+29]= 0.147043*CINTgto_norm(0,env[off+ 5]); env[off+38]= -0.034160*CINTgto_norm(0,env[off+5]);
        env[off+ 6] = 14.750; env[off+30]= 0.305631*CINTgto_norm(0,env[off+ 6]); env[off+39]= -0.077173*CINTgto_norm(0,env[off+6]);
        env[off+ 7] = 6.0350; env[off+31]= 0.399345*CINTgto_norm(0,env[off+ 7]); env[off+40]= -0.141493*CINTgto_norm(0,env[off+7]);
        env[off+ 8] = 2.5300; env[off+32]= 0.217051*CINTgto_norm(0,env[off+ 8]); env[off+41]= -0.118019*CINTgto_norm(0,env[off+8]);
        env[off+ 9] = 0.7355; env[off+42]= 1.000000*CINTgto_norm(0,env[off+ 9]);
        env[off+10] = 0.2905; env[off+43]= 1.000000*CINTgto_norm(0,env[off+10]);
        env[off+11] = 0.1111; env[off+44]= 1.000000*CINTgto_norm(0,env[off+11]);
        env[off+12] = 34.510; env[off+45]= 0.005378*CINTgto_norm(1,env[off+12]);
        env[off+13] = 7.9150; env[off+46]= 0.036132*CINTgto_norm(1,env[off+13]);
        env[off+14] = 2.3680; env[off+47]= 0.142493*CINTgto_norm(1,env[off+14]);
        env[off+15] = 0.8132; env[off+48]= 1.000000*CINTgto_norm(1,env[off+15]);
        env[off+16] = 0.2890; env[off+49]= 1.000000*CINTgto_norm(1,env[off+16]);
        env[off+17] = 0.1007; env[off+50]= 1.000000*CINTgto_norm(1,env[off+17]);
        env[off+18] = 1.8480; env[off+51]= 1.000000*CINTgto_norm(2,env[off+18]);
        env[off+19] = 0.6490; env[off+52]= 1.000000*CINTgto_norm(2,env[off+19]);
        env[off+20] = 0.2280; env[off+53]= 1.000000*CINTgto_norm(2,env[off+20]);
        env[off+21] = 1.4190; env[off+54]= 1.000000*CINTgto_norm(3,env[off+21]);
        env[off+22] = 0.4850; env[off+55]= 1.000000*CINTgto_norm(3,env[off+22]);
        env[off+23] = 1.0110; env[off+56]= 1.000000*CINTgto_norm(4,env[off+23]);
        env[off+57] = 82.64; env[off+69] = 0.002006*CINTgto_norm(0,env[off+57]);
        env[off+58] = 12.41; env[off+70] = 0.015343*CINTgto_norm(0,env[off+58]);
        env[off+59] = 2.824; env[off+71] = 0.075579*CINTgto_norm(0,env[off+59]);
        env[off+60] = 0.797; env[off+72] = 1.000000*CINTgto_norm(0,env[off+60]);
        env[off+61] = 0.258; env[off+73] = 1.000000*CINTgto_norm(0,env[off+61]);
        env[off+62] = 0.089; env[off+74] = 1.000000*CINTgto_norm(0,env[off+62]);
        env[off+63] = 2.292; env[off+75] = 1.000000*CINTgto_norm(1,env[off+63]);
        env[off+64] = 0.838; env[off+76] = 1.000000*CINTgto_norm(1,env[off+64]);
        env[off+65] = 0.292; env[off+77] = 1.000000*CINTgto_norm(1,env[off+65]);
        env[off+66] = 2.062; env[off+78] = 1.000000*CINTgto_norm(2,env[off+66]);
        env[off+67] = 0.662; env[off+79] = 1.000000*CINTgto_norm(2,env[off+67]);
        env[off+68] = 1.397; env[off+80] = 1.000000*CINTgto_norm(3,env[off+68]);
        for (i = 0, ia = 0, n = 0; i < 2; i++) {
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 8;
                bas[NCTR_OF  +BAS_SLOTS*n] = 2;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+ 0;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+24;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+ 9;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+42;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+10;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+43;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 0;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+11;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+44;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+12;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+45;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+15;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+48;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+16;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+49;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 1;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+17;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+50;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 2;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+18;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+51;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 2;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+19;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+52;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 2;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+20;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+53;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 3;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+21;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+54;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 3;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+22;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+55;
                n++;
                bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                bas[ANG_OF   +BAS_SLOTS*n] = 4;
                bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                bas[PTR_EXP  +BAS_SLOTS*n] = off+23;
                bas[PTR_COEFF+BAS_SLOTS*n] = off+56;
                n++;
                ia++;
                for (j = 0; j < 3; j++) {
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 3;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+57;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+69;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+60;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+72;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+61;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+73;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 0;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+62;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+74;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 1;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+63;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+75;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 1;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+64;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+76;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 1;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+65;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+77;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 2;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+66;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+78;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 2;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+67;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+79;
                        n++;
                        bas[ATOM_OF  +BAS_SLOTS*n] = ia;
                        bas[ANG_OF   +BAS_SLOTS*n] = 3;
                        bas[NPRIM_OF +BAS_SLOTS*n] = 1;
                        bas[NCTR_OF  +BAS_SLOTS*n] = 1;
                        bas[PTR_EXP  +BAS_SLOTS*n] = off+68;
                        bas[PTR_COEFF+BAS_SLOTS*n] = off+80;
                        n++;
                        ia++;
                }
        }
        nbas = n;
        printf("cc-pVQZ basis\n");
        run_all(atm, natm, bas, nbas, env);

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

        int pct, count;
        double time0, time1 = 0;
        double tt, tot;
        tot = (double)ncgto*ncgto*ncgto*ncgto/8;
        time0 = omp_get_wtime();

        CINTOpt *non_opt = NULL;
        CINTOpt *opt_for_cint2e = NULL;
        cint2e_sph_optimizer(&opt_for_cint2e, atm, natm, bas, nbas, env);
        CINTOpt *opt_for_ip1 = NULL;
        cint2e_ip1_sph_optimizer(&opt_for_ip1, atm, natm, bas, nbas, env);

        printf("\tcint2e_sph without optimizer: total num ERI = %.2e\n", tot);
        pct = 0; count = 0;
#pragma omp parallel default(none) \
        shared(atm, natm, bas, nbas, env, ishls, jshls, non_opt, time0, pct, count, stdout) \
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
                        cint2e_sph(buf, shls, atm, natm, bas, nbas, env, non_opt);
                        free(buf);
                }
                count += kl_max;
                if (100*count/((long)nbas*nbas*(nbas+1)*(nbas+2)/8) > pct) {
                        pct++;
                        time1 = omp_get_wtime();
                        printf("\t%d%%, CPU time = %8.2f\r", pct, time1-time0);
                        fflush(stdout);
                }
        }
        time1 = omp_get_wtime();
        tt = time1-time0;
        printf("\t100%%, CPU time = %8.2f, %8.4f Mflops\n",
               tt, tot/1e6/tt);
        /* */
        time0 = time1;
        printf("\tcint2e_sph with optimizer: total num ERI = %.2e\n", tot);
        pct = 0; count = 0;

#pragma omp parallel default(none) \
        shared(atm, natm, bas, nbas, env, ishls, jshls, opt_for_cint2e, time0, pct, count, stdout) \
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
                        cint2e_sph(buf, shls, atm, natm, bas, nbas, env, opt_for_cint2e);
                        free(buf);
                }
                count += kl_max;
                if (100*count/((long)nbas*nbas*(nbas+1)*(nbas+2)/8) > pct) {
                        pct++;
                        time1 = omp_get_wtime();
                        printf("\t%d%%, CPU time = %8.2f\r", pct, time1-time0);
                        fflush(stdout);
                }
        }
        time1 = omp_get_wtime();
        tt = time1-time0;
        printf("\t100%%, CPU time = %8.2f, %8.4f Mflops\n",
               tt, tot/1e6/tt);

        /* */
        time0 = time1;
        tot = (double)ncgto*ncgto*ncgto*ncgto/2*3;
        printf("\tGradients with optimizer: total num ERI = %.2e\n", tot);
        pct = 0; count = 0;
#pragma omp parallel default(none) \
        shared(atm, natm, bas, nbas, env, ishls, jshls, opt_for_ip1, time0, pct, count, stdout) \
        private(di, dj, dk, dl, i, j, k, l, ij, kl, shls, buf, time1)
#pragma omp for nowait schedule(dynamic, 2)
        for (ij = 0; ij < nbas*nbas; ij++) {
                i = ij / nbas;
                j = ij - nbas*i;
                di = CINTcgto_spheric(i, bas);
                dj = CINTcgto_spheric(j, bas);
                for (kl = 0; kl < nbas*(nbas+1)/2; kl++) {
                        k = ishls[kl];
                        l = jshls[kl];
                        dk = CINTcgto_spheric(k, bas);
                        dl = CINTcgto_spheric(l, bas);
                        shls[0] = i;
                        shls[1] = j;
                        shls[2] = k;
                        shls[3] = l;
                        buf = malloc(sizeof(double) * di*dj*dk*dl*3);
                        cint2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt_for_ip1);
                        free(buf);
                }
                count += nbas*(nbas+1)/2;
                if (100*count/((long)nbas*nbas*nbas*(nbas+1)/2) > pct) {
                        pct++;
                        time1 = omp_get_wtime();
                        printf("\t%d%%, CPU time = %8.2f\r", pct,
                               time1-time0);
                        fflush(stdout);
                }
        }
        time1 = omp_get_wtime();
        tt = time1-time0;
        printf("\t100%%, CPU time = %8.2f, %8.4f Mflops\n", tt, tot/1e6/tt);

        CINTdel_optimizer(&opt_for_cint2e);
        CINTdel_optimizer(&opt_for_ip1);
        free(ishls);
        free(jshls);
}

