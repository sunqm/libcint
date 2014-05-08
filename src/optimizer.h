/*
 * File: optimizer.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#if !defined HAVE_DEFINED_CINTOPT_H
#define HAVE_DEFINED_CINTOPT_H
typedef struct {
    int **index_xyz_array; // ANG_MAX**4 pointers to index_xyz
    int *prim_offset;
    int *non0ctr;
    int **non0idx;
    double **non0coeff;
    double **expij;
    double **rij;
    int **cceij;
    int tot_prim;
} CINTOpt;
#endif

void CINTinit_2e_optimizer(CINTOpt **opt, const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);
void CINTinit_optimizer(CINTOpt **opt, const int *atm, const int natm,
                        const int *bas, const int nbas, const double *env);
void CINTdel_2e_optimizer(CINTOpt **opt);
void CINTdel_optimizer(CINTOpt **opt);
void CINTOpt_set_index_xyz(CINTOpt *opt, int *ng,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);
void CINTOpt_setij(CINTOpt *opt, const int *atm, const int natm,
                   const int *bas, const int nbas, const double *env);
void CINTOpt_set_non0coeff(CINTOpt *opt, const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);

// optimizer examples
void CINTno_optimizer(CINTOpt **opt, const int *atm, const int natm,
                      const int *bas, const int nbas, const double *env);
void CINTuse_all_optimizer(CINTOpt **opt, int *ng,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);

