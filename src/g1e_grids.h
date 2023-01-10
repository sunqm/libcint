void CINTinit_int1e_grids_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                                  FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);

FINT CINTg0_1e_grids(double *g, double cutoff,
                     CINTEnvVars *envs, double *cache, double *gridsT);

void CINTgout1e_grids(double *gout, double *g, FINT *idx,
                      CINTEnvVars *envs, FINT gout_empty);

void CINTnabla1i_grids(double *f, double *g,
                       FINT li, FINT lj, CINTEnvVars *envs);

void CINTnabla1j_grids(double *f, double *g,
                       FINT li, FINT lj, CINTEnvVars *envs);

void CINTx1i_grids(double *f, double *g, double *ri,
                   FINT li, FINT lj, CINTEnvVars *envs);

void CINTx1j_grids(double *f, double *g, double *rj,
                   FINT li, FINT lj, CINTEnvVars *envs);


#define G1E_GRIDS_D_I(f, g, li, lj)   CINTnabla1i_grids(f, g, li, lj, envs)
#define G1E_GRIDS_D_J(f, g, li, lj)   CINTnabla1j_grids(f, g, li, lj, envs)
/* r-R_0, R_0 is (0,0,0) */
#define G1E_GRIDS_R0I(f, g, li, lj)   CINTx1i_grids(f, g, ri, li, lj, envs)
#define G1E_GRIDS_R0J(f, g, li, lj)   CINTx1j_grids(f, g, rj, li, lj, envs)
/* r-R_C, R_C is common origin */
#define G1E_GRIDS_RCI(f, g, li, lj)   CINTx1i_grids(f, g, dri, li, lj, envs)
#define G1E_GRIDS_RCJ(f, g, li, lj)   CINTx1j_grids(f, g, drj, li, lj, envs)
/* origin from center of each basis */
#define G1E_GRIDS_R_I(f, g, li, lj)   f = g + envs->g_stride_i
#define G1E_GRIDS_R_J(f, g, li, lj)   f = g + envs->g_stride_j
