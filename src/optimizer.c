/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 * optimizer for 2e integrals.  Note if CINT2e_drv is only called a few
 * hundred times, this optimizer cannot really speed up the integration. 
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "cint_config.h"
#include "cint_bas.h"
#include "g1e.h"
#include "g1e_grids.h"
#include "g2e.h"
#include "g3c1e.h"
#include "optimizer.h"
#include "misc.h"

// generate caller to CINTinit_2e_optimizer for each type of function
void CINTinit_2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env)
{
        CINTOpt *opt0 = (CINTOpt *)malloc(sizeof(CINTOpt));
        opt0->index_xyz_array = NULL;
        opt0->non0ctr = NULL;
        opt0->sortedidx = NULL;
        opt0->nbas = nbas;
        opt0->log_max_coeff = NULL;
        opt0->pairdata = NULL;
        *opt = opt0;
}
void CINTinit_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                        FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
}

void CINTdel_2e_optimizer(CINTOpt **opt)
{
        CINTOpt *opt0 = *opt;
        if (opt0 == NULL) { // when opt is created by CINTno_optimizer
                return;
        }

        if (opt0->index_xyz_array != NULL) {
                free(opt0->index_xyz_array[0]);
                free(opt0->index_xyz_array);
        }

        if (opt0->non0ctr != NULL) {
                free(opt0->sortedidx[0]);
                free(opt0->sortedidx);
                free(opt0->non0ctr[0]);
                free(opt0->non0ctr);
        }

        if (opt0->log_max_coeff != NULL) {
                free(opt0->log_max_coeff[0]);
                free(opt0->log_max_coeff);
        }

        CINTdel_pairdata_optimizer(opt0);

        free(opt0);
        *opt = NULL;
}
void CINTdel_optimizer(CINTOpt **opt)
{
        CINTdel_2e_optimizer(opt);
}

void CINTno_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                      FINT *bas, FINT nbas, double *env)
{
        *opt = NULL;
}

static FINT _make_fakebas(FINT *fakebas, FINT *bas, FINT nbas, double *env)
{
        FINT i;
        FINT max_l = 0;
        for (i = 0; i < nbas; i++) {
                max_l = MAX(max_l, bas(ANG_OF,i));
        }

        FINT fakenbas = max_l + 1;
        for (i = 0; i < BAS_SLOTS*fakenbas; i++) {
                fakebas[i] = 0;
        }
        // fakebas only initializes ANG_OF, since the others does not
        // affect index_xyz
        for (i = 0; i <= max_l; i++) {
                fakebas[BAS_SLOTS*i+ANG_OF] = i;
        }
        return max_l;
}
static FINT *_allocate_index_xyz(CINTOpt *opt, FINT max_l, FINT l_allow, FINT order)
{
        FINT i;
        FINT cumcart = (l_allow+1) * (l_allow+2) * (l_allow+3) / 6;
        size_t ll = max_l + 1;
        size_t cc = cumcart;
        for (i = 1; i < order; i++) {
                ll *= LMAX1;
                cc *= cumcart;
        }
        FINT *buf = malloc(sizeof(FINT) * cc * 3);
        FINT **ppbuf = malloc(sizeof(FINT*) * ll);
        ppbuf[0] = buf;
        for (i = 1; i < ll; i++) {
                ppbuf[i] = NULL;
        }
        opt->index_xyz_array = ppbuf;
        return buf;
}
static void gen_idx(CINTOpt *opt, void (*finit)(), void (*findex_xyz)(),
                    FINT order, FINT l_allow, FINT *ng,
                    FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        FINT i, j, k, l, ptr;
        FINT fakebas[BAS_SLOTS*LMAX1];
        FINT max_l = _make_fakebas(fakebas, bas, nbas, env);
        FINT fakenbas = max_l+1;
        // index_xyz bufsize may blow up for large max_l
        l_allow = MIN(max_l, l_allow);
        FINT *buf = _allocate_index_xyz(opt, max_l, l_allow, order);

        CINTEnvVars envs;
        FINT shls[4] = {0,};
        if (order == 2) {
                for (i = 0; i <= l_allow; i++) {
                for (j = 0; j <= l_allow; j++) {
                        shls[0] = i; shls[1] = j;
                        (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                        ptr = i*LMAX1 + j;
                        opt->index_xyz_array[ptr] = buf;
                        (*findex_xyz)(buf, &envs);
                        buf += envs.nf * 3;
                } }

        } else if (order == 3) {
                for (i = 0; i <= l_allow; i++) {
                for (j = 0; j <= l_allow; j++) {
                for (k = 0; k <= l_allow; k++) {
                        shls[0] = i; shls[1] = j; shls[2] = k;
                        (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                        ptr = i*LMAX1*LMAX1 + j*LMAX1 + k;
                        opt->index_xyz_array[ptr] = buf;
                        (*findex_xyz)(buf, &envs);
                        buf += envs.nf * 3;
                } } }

        } else {
                for (i = 0; i <= l_allow; i++) {
                for (j = 0; j <= l_allow; j++) {
                for (k = 0; k <= l_allow; k++) {
                for (l = 0; l <= l_allow; l++) {
                        shls[0] = i; shls[1] = j; shls[2] = k; shls[3] = l;
                        (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                        ptr = i*LMAX1*LMAX1*LMAX1
                            + j*LMAX1*LMAX1
                            + k*LMAX1
                            + l;
                        opt->index_xyz_array[ptr] = buf;
                        (*findex_xyz)(buf, &envs);
                        buf += envs.nf * 3;
                } } } }
        }
}

void CINTall_1e_optimizer(CINTOpt **opt, FINT *ng,
                          FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_set_log_maxc(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int1e_EnvVars, &CINTg1e_index_xyz,
                2, ANG_MAX, ng, atm, natm, bas, nbas, env);
}

void CINTall_2e_optimizer(CINTOpt **opt, FINT *ng,
                          FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int2e_EnvVars, &CINTg2e_index_xyz,
                4, 6, ng, atm, natm, bas, nbas, env);
}

void CINTall_3c2e_optimizer(CINTOpt **opt, FINT *ng,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int3c2e_EnvVars, &CINTg2e_index_xyz,
                3, 12, ng, atm, natm, bas, nbas, env);
}

void CINTall_2c2e_optimizer(CINTOpt **opt, FINT *ng,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_set_log_maxc(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int2c2e_EnvVars, &CINTg1e_index_xyz,
                2, ANG_MAX, ng, atm, natm, bas, nbas, env);
}

void CINTg3c1e_index_xyz(FINT *idx, const CINTEnvVars *envs);
void CINTall_3c1e_optimizer(CINTOpt **opt, FINT *ng,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int3c1e_EnvVars, &CINTg3c1e_index_xyz,
                3, 12, ng, atm, natm, bas, nbas, env);
}

void CINTall_1e_grids_optimizer(CINTOpt **opt, FINT *ng,
                                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_set_log_maxc(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int1e_grids_EnvVars, &CINTg1e_index_xyz,
                2, ANG_MAX, ng, atm, natm, bas, nbas, env);
}

#ifdef WITH_F12
void CINTinit_int2e_stg_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                           FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTall_2e_stg_optimizer(CINTOpt **opt, FINT *ng,
                              FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int2e_stg_EnvVars, &CINTg2e_index_xyz,
                4, 6, ng, atm, natm, bas, nbas, env);
}
#endif

void CINTOpt_log_max_pgto_coeff(double *log_maxc, double *coeff, FINT nprim, FINT nctr)
{
        FINT i, ip;
        double maxc;
        for (ip = 0; ip < nprim; ip++) {
                maxc = 0;
                for (i = 0; i < nctr; i++) {
                        maxc = MAX(maxc, fabs(coeff[i*nprim+ip]));
                }
                log_maxc[ip] = approx_log(maxc);
        }
}


void CINTOpt_set_log_maxc(CINTOpt *opt, FINT *atm, FINT natm,
                          FINT *bas, FINT nbas, double *env)
{
        FINT i, iprim, ictr;
        double *ci;
        size_t tot_prim = 0;
        for (i = 0; i < nbas; i++) {
                tot_prim += bas(NPRIM_OF, i);
        }
        if (tot_prim == 0) {
                return;
        }

        opt->log_max_coeff = malloc(sizeof(double *) * MAX(nbas, 1));
        double *plog_maxc = malloc(sizeof(double) * tot_prim);
        opt->log_max_coeff[0] = plog_maxc;
        for (i = 0; i < nbas; i++) {
                iprim = bas(NPRIM_OF, i);
                ictr = bas(NCTR_OF, i);
                ci = env + bas(PTR_COEFF, i);
                opt->log_max_coeff[i] = plog_maxc;
                CINTOpt_log_max_pgto_coeff(plog_maxc, ci, iprim, ictr);
                plog_maxc += iprim;
        }
}

FINT CINTset_pairdata(PairData *pairdata, double *ai, double *aj, double *ri, double *rj,
                     double *log_maxci, double *log_maxcj,
                     FINT li_ceil, FINT lj_ceil, FINT iprim, FINT jprim,
                     double rr_ij, double expcutoff, double *env)
{
        FINT ip, jp, n;
        double aij, eij, cceij, wj;
        // Normally
        //    (aj*d/sqrt(aij)+1)^li * (ai*d/sqrt(aij)+1)^lj
        //    * pi^1.5/aij^{(li+lj+3)/2} * exp(-ai*aj/aij*rr_ij)
        // is a good approximation for overlap integrals.
        //    <~ (aj*d/aij+1/sqrt(aij))^li * (ai*d/aij+1/sqrt(aij))^lj * (pi/aij)^1.5
        //    <~ (d+1/sqrt(aij))^(li+lj) * (pi/aij)^1.5
        aij = ai[iprim-1] + aj[jprim-1];
        double log_rr_ij = 1.7 - 1.5 * approx_log(aij);
        int lij = li_ceil + lj_ceil;
        if (lij > 0) {
                double dist_ij = sqrt(rr_ij);
                double omega = env[PTR_RANGE_OMEGA];
                if (omega < 0) {
                        double r_guess = 8.;
                        double omega2 = omega * omega;
                        double theta = omega2 / (omega2 + aij);
                        log_rr_ij += lij * approx_log(dist_ij + theta*r_guess + 1.);
                } else {
                        log_rr_ij += lij * approx_log(dist_ij + 1.);
                }
        }
        PairData *pdata;

        FINT empty = 1;
        for (n = 0, jp = 0; jp < jprim; jp++) {
                for (ip = 0; ip < iprim; ip++, n++) {
                        aij = 1/(ai[ip] + aj[jp]);
                        eij = rr_ij * ai[ip] * aj[jp] * aij;
                        cceij = eij - log_rr_ij - log_maxci[ip] - log_maxcj[jp];
                        pdata = pairdata + n;
                        pdata->cceij = cceij;
                        if (cceij < expcutoff) {
                                empty = 0;
                                wj = aj[jp] * aij;
                                pdata->rij[0] = ri[0] + wj * (rj[0]-ri[0]);
                                pdata->rij[1] = ri[1] + wj * (rj[1]-ri[1]);
                                pdata->rij[2] = ri[2] + wj * (rj[2]-ri[2]);
                                pdata->eij = exp(-eij);
                        } else {
                                pdata->rij[0] = 1e18;
                                pdata->rij[1] = 1e18;
                                pdata->rij[2] = 1e18;
                                pdata->eij = 0;
                        }
                }
        }
        return empty;
}

void CINTOpt_setij(CINTOpt *opt, FINT *ng,
                   FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        FINT i, j, ip, jp;
        FINT iprim, jprim, li, lj;
        double *ai, *aj, *ri, *rj;
        double expcutoff;
        if (env[PTR_EXPCUTOFF] == 0) {
                expcutoff = EXPCUTOFF;
        } else {
                expcutoff = MAX(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
        }

        if (opt->log_max_coeff == NULL) {
                CINTOpt_set_log_maxc(opt, atm, natm, bas, nbas, env);
        }
        double **log_max_coeff = opt->log_max_coeff;
        double *log_maxci, *log_maxcj;

        size_t tot_prim = 0;
        for (i = 0; i < nbas; i++) {
                tot_prim += bas(NPRIM_OF, i);
        }
        if (tot_prim == 0 || tot_prim > MAX_PGTO_FOR_PAIRDATA) {
                return;
        }
        opt->pairdata = malloc(sizeof(PairData *) * MAX(nbas * nbas, 1));
        PairData *pdata = malloc(sizeof(PairData) * tot_prim * tot_prim);
        opt->pairdata[0] = pdata;

        FINT ijkl_inc;
        if ((ng[IINC]+ng[JINC]) > (ng[KINC]+ng[LINC])) {
                ijkl_inc = ng[IINC] + ng[JINC];
        } else {
                ijkl_inc = ng[KINC] + ng[LINC];
        }

        FINT empty;
        double rr;
        PairData *pdata0;
        for (i = 0; i < nbas; i++) {
                ri = env + atm(PTR_COORD,bas(ATOM_OF,i));
                ai = env + bas(PTR_EXP,i);
                iprim = bas(NPRIM_OF,i);
                li = bas(ANG_OF,i);
                log_maxci = log_max_coeff[i];

                for (j = 0; j <= i; j++) {
                        rj = env + atm(PTR_COORD,bas(ATOM_OF,j));
                        aj = env + bas(PTR_EXP,j);
                        jprim = bas(NPRIM_OF,j);
                        lj = bas(ANG_OF,j);
                        log_maxcj = log_max_coeff[j];
                        rr = (ri[0]-rj[0])*(ri[0]-rj[0])
                           + (ri[1]-rj[1])*(ri[1]-rj[1])
                           + (ri[2]-rj[2])*(ri[2]-rj[2]);

                        empty = CINTset_pairdata(pdata, ai, aj, ri, rj, log_maxci, log_maxcj,
                                                 li+ijkl_inc, lj, iprim, jprim, rr, expcutoff, env);
                        if (i == 0 && j == 0) {
                                opt->pairdata[0] = pdata;
                                pdata += iprim * jprim;
                        } else if (!empty) {
                                opt->pairdata[i*nbas+j] = pdata;
                                pdata += iprim * jprim;
                                if (i != j) {
                                        opt->pairdata[j*nbas+i] = pdata;
                                        pdata0 = opt->pairdata[i*nbas+j];
                                        // transpose pairdata
                                        for (ip = 0; ip < iprim; ip++) {
                                        for (jp = 0; jp < jprim; jp++, pdata++) {
                                                memcpy(pdata, pdata0+jp*iprim+ip,
                                                       sizeof(PairData));
                                        } }
                                }
                        } else {
                                opt->pairdata[i*nbas+j] = NOVALUE;
                                opt->pairdata[j*nbas+i] = NOVALUE;
                        }
                }
        }
}

void CINTdel_pairdata_optimizer(CINTOpt *cintopt)
{
        if (cintopt != NULL && cintopt->pairdata != NULL) {
                free(cintopt->pairdata[0]);
                free(cintopt->pairdata);
                cintopt->pairdata = NULL;
        }
}

void CINTOpt_non0coeff_byshell(FINT *sortedidx, FINT *non0ctr, double *ci,
                               FINT iprim, FINT ictr)
{
        FINT ip, j, k, kp;
        FINT zeroidx[ictr];
        for (ip = 0; ip < iprim; ip++) {
                for (j = 0, k = 0, kp = 0; j < ictr; j++) {
                        if (ci[iprim*j+ip] != 0) {
                                sortedidx[k] = j;
                                k++;
                        } else {
                                zeroidx[kp] = j;
                                kp++;
                        }
                }
// Append the index of zero-coeff to sortedidx for function CINTprim_to_ctr_0
                for (j = 0; j < kp; j++) {
                        sortedidx[k+j] = zeroidx[j];
                }
                non0ctr[ip] = k;
                sortedidx += ictr;
        }
}

void CINTOpt_set_non0coeff(CINTOpt *opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env)
{
        FINT i, iprim, ictr;
        double *ci;
        size_t tot_prim = 0;
        size_t tot_prim_ctr = 0;
        for (i = 0; i < nbas; i++) {
                tot_prim += bas(NPRIM_OF, i);
                tot_prim_ctr += bas(NPRIM_OF, i) * bas(NCTR_OF,i);
        }
        if (tot_prim == 0) {
                return;
        }

        opt->non0ctr = malloc(sizeof(FINT *) * MAX(nbas, 1));
        opt->sortedidx = malloc(sizeof(FINT *) * MAX(nbas, 1));
        FINT *pnon0ctr = malloc(sizeof(FINT) * tot_prim);
        FINT *psortedidx = malloc(sizeof(FINT) * tot_prim_ctr);
        opt->non0ctr[0] = pnon0ctr;
        opt->sortedidx[0] = psortedidx;
        for (i = 0; i < nbas; i++) {
                iprim = bas(NPRIM_OF, i);
                ictr = bas(NCTR_OF, i);
                ci = env + bas(PTR_COEFF, i);
                opt->non0ctr[i] = pnon0ctr;
                opt->sortedidx[i] = psortedidx;
                CINTOpt_non0coeff_byshell(psortedidx, pnon0ctr, ci, iprim, ictr);
                pnon0ctr += iprim;
                psortedidx += iprim * ictr;
        }
}
