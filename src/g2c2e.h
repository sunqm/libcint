#include "g2e.h"

FINT CINTinit_int2c2e_EnvVars(CINTEnvVars *envs, const FINT *ng, const FINT *shls,
                             const FINT *atm, const FINT natm,
                             const FINT *bas, const FINT nbas, const double *env);

void CINTnabla1i_2c2e(double *f, const double *g,
                      const FINT li, const FINT lk, const CINTEnvVars *envs);

void CINTnabla1k_2c2e(double *f, const double *g,
                      const FINT li, const FINT lk, const CINTEnvVars *envs);

void CINTx1i_2c2e(double *f, const double *g, const double *ri,
                  const FINT li, const FINT lk, const CINTEnvVars *envs);

void CINTx1k_2c2e(double *f, const double *g, const double *rk,
                  const FINT li, const FINT lk, const CINTEnvVars *envs);


#define G2C2E_D_I(f, g, li, lk)   CINTnabla1i_2c2e(f, g, li, lk, envs)
#define G2C2E_D_K(f, g, li, lk)   CINTnabla1k_2c2e(f, g, li, lk, envs)
/* r-R_0, R_0 is (0,0,0) */
#define G2C2E_R0I(f, g, li, lk)   CINTx1i_2c2e(f, g, ri, li, lk, envs)
#define G2C2E_R0K(f, g, li, lk)   CINTx1k_2c2e(f, g, rk, li, lk, envs)
/* r-R_C, R_C is common origin */
#define G2C2E_RCI(f, g, li, lk)   CINTx1i_2c2e(f, g, dri, li, lk, envs)
#define G2C2E_RCK(f, g, li, lk)   CINTx1k_2c2e(f, g, drk, li, lk, envs)
/* origin from center of each basis
 * x1[ijk]_2e(f, g, ng, li, lk, 0d0) */
#define G2C2E_R_I(f, g, li, lk)   f = g + envs->g_stride_i
#define G2C2E_R_K(f, g, li, lk)   f = g + envs->g_stride_k
