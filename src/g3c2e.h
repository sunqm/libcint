
FINT CINTinit_int3c2e_EnvVars(CINTEnvVars *envs, const FINT *ng, const FINT *shls,
                             const FINT *atm, const FINT natm,
                             const FINT *bas, const FINT nbas, const double *env);

FINT CINTinit_int2c2e_EnvVars(CINTEnvVars *envs, const FINT *ng, const FINT *shls,
                             const FINT *atm, const FINT natm,
                             const FINT *bas, const FINT nbas, const double *env);

void CINTg3c2e_index_xyz(FINT *idx, const CINTEnvVars *envs);

void CINTnabla1i_3c2e(double *f, const double *g,
                      const FINT li, const FINT lj, const FINT lk,
                      const CINTEnvVars *envs);

void CINTnabla1j_3c2e(double *f, const double *g,
                      const FINT li, const FINT lj, const FINT lk,
                      const CINTEnvVars *envs);

void CINTx1i_3c2e(double *f, const double *g, const double *ri,
                  const FINT li, const FINT lj, const FINT lk,
                  const CINTEnvVars *envs);

void CINTx1j_3c2e(double *f, const double *g, const double *rj,
                  const FINT li, const FINT lj, const FINT lk,
                  const CINTEnvVars *envs);


#define G3C2E_D_I(f, g, li, lj, lk)   CINTnabla1i_3c2e(f, g, li, lj, lk, envs)
#define G3C2E_D_J(f, g, li, lj, lk)   CINTnabla1j_3c2e(f, g, li, lj, lk, envs)
/* r-R_0, R_0 is (0,0,0) */
#define G3C2E_R0I(f, g, li, lj, lk)   CINTx1i_3c2e(f, g, ri, li, lj, lk, envs)
#define G3C2E_R0J(f, g, li, lj, lk)   CINTx1j_3c2e(f, g, rj, li, lj, lk, envs)
/* r-R_C, R_C is common origin */
#define G3C2E_RCI(f, g, li, lj, lk)   CINTx1i_3c2e(f, g, dri, li, lj, lk, envs)
#define G3C2E_RCJ(f, g, li, lj, lk)   CINTx1j_3c2e(f, g, drj, li, lj, lk, envs)
/* origin from center of each basis
 * x1[ijk]_2e(f, g, ng, li, lj, lk, 0d0) */
#define G3C2E_R_I(f, g, li, lj, lk)   f = g + envs->g_stride_i
#define G3C2E_R_J(f, g, li, lj, lk)   f = g + envs->g_stride_j
