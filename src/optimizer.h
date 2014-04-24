/*
 * File: optimizer.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#if !defined HAVE_DEFINED_CINTOPT_H
#define HAVE_DEFINED_CINTOPT_H
typedef struct {
    unsigned int **index_xyz_array; // ANG_MAX**4 pointers to index_xyz
    unsigned int *prim_offset;
    double **eij;
    double **rij;
    int **screenij;
    unsigned int tot_prim;
} CINTOpt;
#endif

void CINTinit_2e_optimizer(CINTOpt **opt, const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);
void CINTinit_optimizer(CINTOpt **opt, const int *atm, const int natm,
                        const int *bas, const int nbas, const double *env);
void CINTdel_2e_optimizer(CINTOpt **opt);
void CINTdel_optimizer(CINTOpt **opt);
void CINTOpt_set_index_xyz(CINTOpt *opt, unsigned int *ng,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);
void CINTOpt_setij(CINTOpt *opt, const int *atm, const int natm,
                   const int *bas, const int nbas, const double *env);

// optimizer examples
void CINTno_optimizer(CINTOpt **opt, const int *atm, const int natm,
                      const int *bas, const int nbas, const double *env);
void CINTuse_all_optimizer(CINTOpt **opt, unsigned int *ng,
                           const int *atm, const int natm,
                           const int *bas, const int nbas, const double *env);

