void CINTinit_int3c1e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                              FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);

void CINTg3c1e_ovlp(double *g, double ai, double aj, double ak,
                    CINTEnvVars *envs);
void CINTg3c1e_nuc(double *g, double ai, double aj, double ak, double *rijk,
                   double *cr, double t2, CINTEnvVars *envs);

void CINTnabla1i_3c1e(double *f, const double *g,
                      const FINT li, const FINT lj, const FINT lk,
                      const CINTEnvVars *envs);

void CINTnabla1j_3c1e(double *f, const double *g,
                      const FINT li, const FINT lj, const FINT lk,
                      const CINTEnvVars *envs);

void CINTnabla1k_3c1e(double *f, const double *g,
                      const FINT li, const FINT lj, const FINT lk,
                      const CINTEnvVars *envs);
void CINTx1i_3c1e(double *f, const double *g, const double *ri,
                  const FINT li, const FINT lj, const FINT lk,
                  const CINTEnvVars *envs);

void CINTx1j_3c1e(double *f, const double *g, const double *rj,
                  const FINT li, const FINT lj, const FINT lk,
                  const CINTEnvVars *envs);

void CINTx1k_3c1e(double *f, const double *g, const double *rk,
                  const FINT li, const FINT lj, const FINT lk,
                  const CINTEnvVars *envs);


#define G3C1E_D_I(f, g, li, lj, lk)   CINTnabla1i_3c1e(f, g, li, lj, lk, envs)
#define G3C1E_D_J(f, g, li, lj, lk)   CINTnabla1j_3c1e(f, g, li, lj, lk, envs)
#define G3C1E_D_K(f, g, li, lj, lk)   CINTnabla1k_3c1e(f, g, li, lj, lk, envs)
/* r-R_0, R_0 is (0,0,0) */
#define G3C1E_R0I(f, g, li, lj, lk)   CINTx1i_3c1e(f, g, ri, li, lj, lk, envs)
#define G3C1E_R0J(f, g, li, lj, lk)   CINTx1j_3c1e(f, g, rj, li, lj, lk, envs)
#define G3C1E_R0K(f, g, li, lj, lk)   CINTx1k_3c1e(f, g, rk, li, lj, lk, envs)
/* r-R_C, R_C is common origin */
#define G3C1E_RCI(f, g, li, lj, lk)   CINTx1i_3c1e(f, g, dri, li, lj, lk, envs)
#define G3C1E_RCJ(f, g, li, lj, lk)   CINTx1j_3c1e(f, g, drj, li, lj, lk, envs)
#define G3C1E_RCK(f, g, li, lj, lk)   CINTx1k_3c1e(f, g, drk, li, lj, lk, envs)
/* origin from center of each basis */
#define G3C1E_R_I(f, g, li, lj, lk)   f = g + envs->g_stride_i
#define G3C1E_R_J(f, g, li, lj, lk)   f = g + envs->g_stride_j
#define G3C1E_R_K(f, g, li, lj, lk)   f = g + envs->g_stride_k
