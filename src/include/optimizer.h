/*
 * File: optimizer.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#include "cint_bas.h"

void CINTinit_2e_optimizer(CINTOpt **opt, const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);
void CINTdel_2e_optimizer(CINTOpt *opt);
void CINTOpt_set_log_coeff(CINTOpt *opt, const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);
void CINTOpt_set_index_xyz(CINTOpt *opt, unsigned int *ng,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);

// optimizer examples
void CINTno_optimizer(CINTOpt **opt, const int *atm, const int natm,
                      const int *bas, const int nbas, const double *env);
void CINTcutoff_optimizer(CINTOpt **opt, const int *atm, const int natm,
                          const int *bas, const int nbas, const double *env);
void CINTuse_all_optimizer(CINTOpt **opt, unsigned int *ng,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);

