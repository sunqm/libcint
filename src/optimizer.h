/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 */

#include "cint.h"

#define NOVALUE                 ((void *)0xffffffffffffffffuL)
#define MAX_PGTO_FOR_PAIRDATA   2048

void CINTinit_2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env);
void CINTinit_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                        FINT *bas, FINT nbas, double *env);
void CINTdel_2e_optimizer(CINTOpt **opt);
void CINTdel_optimizer(CINTOpt **opt);
void CINTdel_pairdata_optimizer(CINTOpt *cintopt);
void CINTOpt_log_max_pgto_coeff(double *log_maxc, double *coeff, FINT nprim, FINT nctr);
void CINTOpt_set_log_maxc(CINTOpt *opt, FINT *atm, FINT natm,
                          FINT *bas, FINT nbas, double *env);
void CINTOpt_setij(CINTOpt *opt, FINT *ng,
                   FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTOpt_non0coeff_byshell(FINT *sortedidx, FINT *non0ctr, double *ci,
                               FINT iprim, FINT ictr);
void CINTOpt_set_non0coeff(CINTOpt *opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env);
FINT CINTset_pairdata(PairData *pairdata, double *ai, double *aj, double *ri, double *rj,
                      double *log_maxci, double *log_maxcj,
                      FINT li_ceil, FINT lj_ceil, FINT iprim, FINT jprim,
                      double rr_ij, double expcutoff, double *env);

void CINTOpt_4cindex_xyz(CINTOpt *opt, FINT *ng,
                         FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTOpt_3cindex_xyz(CINTOpt *opt, FINT *ng,
                         FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTOpt_2cindex_xyz(CINTOpt *opt, FINT *ng,
                         FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTOpt_3c1eindex_xyz(CINTOpt *opt, FINT *ng,
                           FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);

// optimizer examples
void CINTno_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                      FINT *bas, FINT nbas, double *env);
void CINTall_1e_optimizer(CINTOpt **opt, FINT *ng,
                          FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTall_2e_optimizer(CINTOpt **opt, FINT *ng,
                          FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTall_3c2e_optimizer(CINTOpt **opt, FINT *ng,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTall_2c2e_optimizer(CINTOpt **opt, FINT *ng,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTall_3c1e_optimizer(CINTOpt **opt, FINT *ng,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
void CINTall_1e_grids_optimizer(CINTOpt **opt, FINT *ng,
                                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);

#ifdef WITH_F12
void CINTall_2e_stg_optimizer(CINTOpt **opt, FINT *ng,
                              FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
#endif

#ifndef HAVE_DEFINED_APPROX_LOG
#define HAVE_DEFINED_APPROX_LOG
#ifdef __X86__
//// little endian on x86
//typedef union {
//    double d;
//    unsigned short s[4];
//} type_IEEE754;
//// ~4 times faster than built-in log
//static inline double approx_log(double x)
//{
//        type_IEEE754 y;
//        y.d = x;
//        return ((y.s[3] >> 4) - 1023 + 1) * 0.693145751953125;
//}
#define approx_log      log
#else
#define approx_log      log
#endif
#endif
