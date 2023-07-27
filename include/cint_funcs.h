/*
 * Copyright (C) 2019-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Function signature
 */


#include <cint.h>

#if !defined HAVE_DEFINED_CINTINTEGRALFUNCTION
#define HAVE_DEFINED_CINTINTEGRALFUNCTION
typedef void CINTOptimizerFunction(CINTOpt **opt,
                                   FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
typedef CACHE_SIZE_T CINTIntegralFunction(double *out, FINT *dims, FINT *shls,
                                  FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
                                  CINTOpt *opt, double *cache);
#endif

/* Plain ERI (ij|kl) */
extern CINTOptimizerFunction int2e_optimizer;
extern CINTIntegralFunction int2e_cart;
extern CINTIntegralFunction int2e_sph;
extern CINTIntegralFunction int2e_spinor;

/* <i|OVLP |j> */
extern CINTOptimizerFunction int1e_ovlp_optimizer;
extern CINTIntegralFunction int1e_ovlp_cart;
extern CINTIntegralFunction int1e_ovlp_sph;
extern CINTIntegralFunction int1e_ovlp_spinor;

/* <i|NUC |j> */
extern CINTOptimizerFunction int1e_nuc_optimizer;
extern CINTIntegralFunction int1e_nuc_cart;
extern CINTIntegralFunction int1e_nuc_sph;
extern CINTIntegralFunction int1e_nuc_spinor;

/* <i|OVLP |P DOT P j> */
extern CINTOptimizerFunction int1e_kin_optimizer;
extern CINTIntegralFunction int1e_kin_cart;
extern CINTIntegralFunction int1e_kin_sph;
extern CINTIntegralFunction int1e_kin_spinor;

/* <i|NABLA-RINV |CROSS P j> */
extern CINTOptimizerFunction int1e_ia01p_optimizer;
extern CINTIntegralFunction int1e_ia01p_cart;
extern CINTIntegralFunction int1e_ia01p_sph;
extern CINTIntegralFunction int1e_ia01p_spinor;

/* <i|OVLP |R CROSS P j> */
extern CINTOptimizerFunction int1e_giao_irjxp_optimizer;
extern CINTIntegralFunction int1e_giao_irjxp_cart;
extern CINTIntegralFunction int1e_giao_irjxp_sph;
extern CINTIntegralFunction int1e_giao_irjxp_spinor;

/* <i|OVLP |RC CROSS P j> */
extern CINTOptimizerFunction int1e_cg_irxp_optimizer;
extern CINTIntegralFunction int1e_cg_irxp_cart;
extern CINTIntegralFunction int1e_cg_irxp_sph;
extern CINTIntegralFunction int1e_cg_irxp_spinor;

/* <i|NABLA-RINV |R j> */
extern CINTOptimizerFunction int1e_giao_a11part_optimizer;
extern CINTIntegralFunction int1e_giao_a11part_cart;
extern CINTIntegralFunction int1e_giao_a11part_sph;
extern CINTIntegralFunction int1e_giao_a11part_spinor;

/* <i|NABLA-RINV |RC j> */
extern CINTOptimizerFunction int1e_cg_a11part_optimizer;
extern CINTIntegralFunction int1e_cg_a11part_cart;
extern CINTIntegralFunction int1e_cg_a11part_sph;
extern CINTIntegralFunction int1e_cg_a11part_spinor;

/* <G i|NABLA-RINV CROSS P |j> */
extern CINTOptimizerFunction int1e_a01gp_optimizer;
extern CINTIntegralFunction int1e_a01gp_cart;
extern CINTIntegralFunction int1e_a01gp_sph;
extern CINTIntegralFunction int1e_a01gp_spinor;

/* <G i|OVLP |P DOT P j> */
extern CINTOptimizerFunction int1e_igkin_optimizer;
extern CINTIntegralFunction int1e_igkin_cart;
extern CINTIntegralFunction int1e_igkin_sph;
extern CINTIntegralFunction int1e_igkin_spinor;

/* <G i|OVLP |j> */
extern CINTOptimizerFunction int1e_igovlp_optimizer;
extern CINTIntegralFunction int1e_igovlp_cart;
extern CINTIntegralFunction int1e_igovlp_sph;
extern CINTIntegralFunction int1e_igovlp_spinor;

/* <G i|NUC |j> */
extern CINTOptimizerFunction int1e_ignuc_optimizer;
extern CINTIntegralFunction int1e_ignuc_cart;
extern CINTIntegralFunction int1e_ignuc_sph;
extern CINTIntegralFunction int1e_ignuc_spinor;

/* <P* i|NUC DOT P |j> */
extern CINTOptimizerFunction int1e_pnucp_optimizer;
extern CINTIntegralFunction int1e_pnucp_cart;
extern CINTIntegralFunction int1e_pnucp_sph;
extern CINTIntegralFunction int1e_pnucp_spinor;

/* <i|ZC |j> */
extern CINTOptimizerFunction int1e_z_optimizer;
extern CINTIntegralFunction int1e_z_cart;
extern CINTIntegralFunction int1e_z_sph;
extern CINTIntegralFunction int1e_z_spinor;

/* <i|ZC ZC |j> */
extern CINTOptimizerFunction int1e_zz_optimizer;
extern CINTIntegralFunction int1e_zz_cart;
extern CINTIntegralFunction int1e_zz_sph;
extern CINTIntegralFunction int1e_zz_spinor;

/* <i|RC |j> */
extern CINTOptimizerFunction int1e_r_optimizer;
extern CINTIntegralFunction int1e_r_cart;
extern CINTIntegralFunction int1e_r_sph;
extern CINTIntegralFunction int1e_r_spinor;

/* <i|RC DOT RC |j> */
extern CINTOptimizerFunction int1e_r2_optimizer;
extern CINTIntegralFunction int1e_r2_cart;
extern CINTIntegralFunction int1e_r2_sph;
extern CINTIntegralFunction int1e_r2_spinor;

/* <i|RC DOT RC RC DOT RC |j> */
extern CINTOptimizerFunction int1e_r4_optimizer;
extern CINTIntegralFunction int1e_r4_cart;
extern CINTIntegralFunction int1e_r4_sph;
extern CINTIntegralFunction int1e_r4_spinor;

/* <i|RC RC |j> */
extern CINTOptimizerFunction int1e_rr_optimizer;
extern CINTIntegralFunction int1e_rr_cart;
extern CINTIntegralFunction int1e_rr_sph;
extern CINTIntegralFunction int1e_rr_spinor;

/* <i|RC RC RC |j> */
extern CINTOptimizerFunction int1e_rrr_optimizer;
extern CINTIntegralFunction int1e_rrr_cart;
extern CINTIntegralFunction int1e_rrr_sph;
extern CINTIntegralFunction int1e_rrr_spinor;

/* <i|RC RC RC RC |j> */
extern CINTOptimizerFunction int1e_rrrr_optimizer;
extern CINTIntegralFunction int1e_rrrr_cart;
extern CINTIntegralFunction int1e_rrrr_sph;
extern CINTIntegralFunction int1e_rrrr_spinor;

/* <i|Z |j> */
extern CINTOptimizerFunction int1e_z_origj_optimizer;
extern CINTIntegralFunction int1e_z_origj_cart;
extern CINTIntegralFunction int1e_z_origj_sph;
extern CINTIntegralFunction int1e_z_origj_spinor;

/* <i|Z Z |j> */
extern CINTOptimizerFunction int1e_zz_origj_optimizer;
extern CINTIntegralFunction int1e_zz_origj_cart;
extern CINTIntegralFunction int1e_zz_origj_sph;
extern CINTIntegralFunction int1e_zz_origj_spinor;

/* <i|R |j> */
extern CINTOptimizerFunction int1e_r_origj_optimizer;
extern CINTIntegralFunction int1e_r_origj_cart;
extern CINTIntegralFunction int1e_r_origj_sph;
extern CINTIntegralFunction int1e_r_origj_spinor;

/* <i|R R |j> */
extern CINTOptimizerFunction int1e_rr_origj_optimizer;
extern CINTIntegralFunction int1e_rr_origj_cart;
extern CINTIntegralFunction int1e_rr_origj_sph;
extern CINTIntegralFunction int1e_rr_origj_spinor;

/* <i|R DOT R |j> */
extern CINTOptimizerFunction int1e_r2_origj_optimizer;
extern CINTIntegralFunction int1e_r2_origj_cart;
extern CINTIntegralFunction int1e_r2_origj_sph;
extern CINTIntegralFunction int1e_r2_origj_spinor;

/* <i|OVLP |R DOT R R DOT R j> */
extern CINTOptimizerFunction int1e_r4_origj_optimizer;
extern CINTIntegralFunction int1e_r4_origj_cart;
extern CINTIntegralFunction int1e_r4_origj_sph;
extern CINTIntegralFunction int1e_r4_origj_spinor;

/* <P DOT P i|OVLP |P DOT P j> */
extern CINTOptimizerFunction int1e_p4_optimizer;
extern CINTIntegralFunction int1e_p4_cart;
extern CINTIntegralFunction int1e_p4_sph;
extern CINTIntegralFunction int1e_p4_spinor;

/* <P* i|RINV DOT P |j> */
extern CINTOptimizerFunction int1e_prinvp_optimizer;
extern CINTIntegralFunction int1e_prinvp_cart;
extern CINTIntegralFunction int1e_prinvp_sph;
extern CINTIntegralFunction int1e_prinvp_spinor;

/* <P* i|RINV CROSS P |j> */
extern CINTOptimizerFunction int1e_prinvxp_optimizer;
extern CINTIntegralFunction int1e_prinvxp_cart;
extern CINTIntegralFunction int1e_prinvxp_sph;
extern CINTIntegralFunction int1e_prinvxp_spinor;

/* <P* i|NUC CROSS P |j> */
extern CINTOptimizerFunction int1e_pnucxp_optimizer;
extern CINTIntegralFunction int1e_pnucxp_cart;
extern CINTIntegralFunction int1e_pnucxp_sph;
extern CINTIntegralFunction int1e_pnucxp_spinor;

/* <i|RC NABLA |j> */
extern CINTOptimizerFunction int1e_irp_optimizer;
extern CINTIntegralFunction int1e_irp_cart;
extern CINTIntegralFunction int1e_irp_sph;
extern CINTIntegralFunction int1e_irp_spinor;

/* <i|RC RC NABLA |j> */
extern CINTOptimizerFunction int1e_irrp_optimizer;
extern CINTIntegralFunction int1e_irrp_cart;
extern CINTIntegralFunction int1e_irrp_sph;
extern CINTIntegralFunction int1e_irrp_spinor;

/* <i|RC NABLA RC |j> */
extern CINTOptimizerFunction int1e_irpr_optimizer;
extern CINTIntegralFunction int1e_irpr_cart;
extern CINTIntegralFunction int1e_irpr_sph;
extern CINTIntegralFunction int1e_irpr_spinor;

/* <i|G G |j> */
extern CINTOptimizerFunction int1e_ggovlp_optimizer;
extern CINTIntegralFunction int1e_ggovlp_cart;
extern CINTIntegralFunction int1e_ggovlp_sph;
extern CINTIntegralFunction int1e_ggovlp_spinor;

/* <i|G G P DOT P |j> */
extern CINTOptimizerFunction int1e_ggkin_optimizer;
extern CINTIntegralFunction int1e_ggkin_cart;
extern CINTIntegralFunction int1e_ggkin_sph;
extern CINTIntegralFunction int1e_ggkin_spinor;

/* <i|G G NUC |j> */
extern CINTOptimizerFunction int1e_ggnuc_optimizer;
extern CINTIntegralFunction int1e_ggnuc_cart;
extern CINTIntegralFunction int1e_ggnuc_sph;
extern CINTIntegralFunction int1e_ggnuc_spinor;

/* <i|G R CROSS P |j> */
extern CINTOptimizerFunction int1e_grjxp_optimizer;
extern CINTIntegralFunction int1e_grjxp_cart;
extern CINTIntegralFunction int1e_grjxp_sph;
extern CINTIntegralFunction int1e_grjxp_spinor;

/* <i|RINV |j> */
extern CINTOptimizerFunction int1e_rinv_optimizer;
extern CINTIntegralFunction int1e_rinv_cart;
extern CINTIntegralFunction int1e_rinv_sph;
extern CINTIntegralFunction int1e_rinv_spinor;

/* <i|NABLA-RINV |j> */
extern CINTOptimizerFunction int1e_drinv_optimizer;
extern CINTIntegralFunction int1e_drinv_cart;
extern CINTIntegralFunction int1e_drinv_sph;
extern CINTIntegralFunction int1e_drinv_spinor;

/* (G i j|R12 |k l) */
extern CINTOptimizerFunction int2e_ig1_optimizer;
extern CINTIntegralFunction int2e_ig1_cart;
extern CINTIntegralFunction int2e_ig1_sph;
extern CINTIntegralFunction int2e_ig1_spinor;

/* (G G i j|R12 |k l) */
extern CINTOptimizerFunction int2e_gg1_optimizer;
extern CINTIntegralFunction int2e_gg1_cart;
extern CINTIntegralFunction int2e_gg1_sph;
extern CINTIntegralFunction int2e_gg1_spinor;

/* (G i j|R12 |G k l) */
extern CINTOptimizerFunction int2e_g1g2_optimizer;
extern CINTIntegralFunction int2e_g1g2_cart;
extern CINTIntegralFunction int2e_g1g2_sph;
extern CINTIntegralFunction int2e_g1g2_spinor;

/* (P* i CROSS P j|R12 |k l) */
extern CINTOptimizerFunction int2e_p1vxp1_optimizer;
extern CINTIntegralFunction int2e_p1vxp1_cart;
extern CINTIntegralFunction int2e_p1vxp1_sph;
extern CINTIntegralFunction int2e_p1vxp1_spinor;

/* (i RC j|NABLA-R12 |k l) */
extern CINTOptimizerFunction int2e_ip1v_rc1_optimizer;
extern CINTIntegralFunction int2e_ip1v_rc1_cart;
extern CINTIntegralFunction int2e_ip1v_rc1_sph;
extern CINTIntegralFunction int2e_ip1v_rc1_spinor;

/* (i R j|NABLA-R12 |k l) */
extern CINTOptimizerFunction int2e_ip1v_r1_optimizer;
extern CINTIntegralFunction int2e_ip1v_r1_cart;
extern CINTIntegralFunction int2e_ip1v_r1_sph;
extern CINTIntegralFunction int2e_ip1v_r1_spinor;

/* (G i j|NABLA-R12 CROSS P |k l) */
extern CINTOptimizerFunction int2e_ipvg1_xp1_optimizer;
extern CINTIntegralFunction int2e_ipvg1_xp1_cart;
extern CINTIntegralFunction int2e_ipvg1_xp1_sph;
extern CINTIntegralFunction int2e_ipvg1_xp1_spinor;

/* (i j|NABLA-R12 CROSS P |G k l) */
extern CINTOptimizerFunction int2e_ipvg2_xp1_optimizer;
extern CINTIntegralFunction int2e_ipvg2_xp1_cart;
extern CINTIntegralFunction int2e_ipvg2_xp1_sph;
extern CINTIntegralFunction int2e_ipvg2_xp1_spinor;

/* <i|NUC |RC CROSS P j> */
extern CINTOptimizerFunction int1e_inuc_rcxp_optimizer;
extern CINTIntegralFunction int1e_inuc_rcxp_cart;
extern CINTIntegralFunction int1e_inuc_rcxp_sph;
extern CINTIntegralFunction int1e_inuc_rcxp_spinor;

/* <i|NUC |R CROSS P j> */
extern CINTOptimizerFunction int1e_inuc_rxp_optimizer;
extern CINTIntegralFunction int1e_inuc_rxp_cart;
extern CINTIntegralFunction int1e_inuc_rxp_sph;
extern CINTIntegralFunction int1e_inuc_rxp_spinor;

/* <i|OVLP |SIGMA j> */
extern CINTOptimizerFunction int1e_sigma_optimizer;
extern CINTIntegralFunction int1e_sigma_cart;
extern CINTIntegralFunction int1e_sigma_sph;
extern CINTIntegralFunction int1e_sigma_spinor;

/* <SIGMA DOT P i|OVLP |SIGMA SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_spsigmasp_optimizer;
extern CINTIntegralFunction int1e_spsigmasp_cart;
extern CINTIntegralFunction int1e_spsigmasp_sph;
extern CINTIntegralFunction int1e_spsigmasp_spinor;

/* <SIGMA DOT R i|OVLP |SIGMA DOT R j> */
extern CINTOptimizerFunction int1e_srsr_optimizer;
extern CINTIntegralFunction int1e_srsr_cart;
extern CINTIntegralFunction int1e_srsr_sph;
extern CINTIntegralFunction int1e_srsr_spinor;

/* <SIGMA DOT R i|OVLP |j> */
extern CINTOptimizerFunction int1e_sr_optimizer;
extern CINTIntegralFunction int1e_sr_cart;
extern CINTIntegralFunction int1e_sr_sph;
extern CINTIntegralFunction int1e_sr_spinor;

/* <SIGMA DOT R i|OVLP |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_srsp_optimizer;
extern CINTIntegralFunction int1e_srsp_cart;
extern CINTIntegralFunction int1e_srsp_sph;
extern CINTIntegralFunction int1e_srsp_spinor;

/* <SIGMA DOT P i|OVLP |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_spsp_optimizer;
extern CINTIntegralFunction int1e_spsp_cart;
extern CINTIntegralFunction int1e_spsp_sph;
extern CINTIntegralFunction int1e_spsp_spinor;

/* <SIGMA DOT P i|OVLP |j> */
extern CINTOptimizerFunction int1e_sp_optimizer;
extern CINTIntegralFunction int1e_sp_cart;
extern CINTIntegralFunction int1e_sp_sph;
extern CINTIntegralFunction int1e_sp_spinor;

/* <SIGMA DOT P i|NUC |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_spnucsp_optimizer;
extern CINTIntegralFunction int1e_spnucsp_cart;
extern CINTIntegralFunction int1e_spnucsp_sph;
extern CINTIntegralFunction int1e_spnucsp_spinor;

/* <SIGMA DOT P i|RINV |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_sprinvsp_optimizer;
extern CINTIntegralFunction int1e_sprinvsp_cart;
extern CINTIntegralFunction int1e_sprinvsp_sph;
extern CINTIntegralFunction int1e_sprinvsp_spinor;

/* <SIGMA DOT R i|NUC |SIGMA DOT R j> */
extern CINTOptimizerFunction int1e_srnucsr_optimizer;
extern CINTIntegralFunction int1e_srnucsr_cart;
extern CINTIntegralFunction int1e_srnucsr_sph;
extern CINTIntegralFunction int1e_srnucsr_spinor;

/* <SIGMA DOT P i|RC |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_sprsp_optimizer;
extern CINTIntegralFunction int1e_sprsp_cart;
extern CINTIntegralFunction int1e_sprsp_sph;
extern CINTIntegralFunction int1e_sprsp_spinor;

/* <G i|OVLP |j> */
extern CINTOptimizerFunction int1e_govlp_optimizer;
extern CINTIntegralFunction int1e_govlp_cart;
extern CINTIntegralFunction int1e_govlp_sph;
extern CINTIntegralFunction int1e_govlp_spinor;

/* <G i|NUC |j> */
extern CINTOptimizerFunction int1e_gnuc_optimizer;
extern CINTIntegralFunction int1e_gnuc_cart;
extern CINTIntegralFunction int1e_gnuc_sph;
extern CINTIntegralFunction int1e_gnuc_spinor;

/* <SIGMA CROSS RC i|SIGMA CROSS NABLA-RINV |j> */
extern CINTOptimizerFunction int1e_cg_sa10sa01_optimizer;
extern CINTIntegralFunction int1e_cg_sa10sa01_cart;
extern CINTIntegralFunction int1e_cg_sa10sa01_sph;
extern CINTIntegralFunction int1e_cg_sa10sa01_spinor;

/* <RC CROSS SIGMA i|OVLP |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_cg_sa10sp_optimizer;
extern CINTIntegralFunction int1e_cg_sa10sp_cart;
extern CINTIntegralFunction int1e_cg_sa10sp_sph;
extern CINTIntegralFunction int1e_cg_sa10sp_spinor;

/* <RC CROSS SIGMA i|NUC |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_cg_sa10nucsp_optimizer;
extern CINTIntegralFunction int1e_cg_sa10nucsp_cart;
extern CINTIntegralFunction int1e_cg_sa10nucsp_sph;
extern CINTIntegralFunction int1e_cg_sa10nucsp_spinor;

/* <SIGMA CROSS R i|SIGMA CROSS NABLA-RINV |j> */
extern CINTOptimizerFunction int1e_giao_sa10sa01_optimizer;
extern CINTIntegralFunction int1e_giao_sa10sa01_cart;
extern CINTIntegralFunction int1e_giao_sa10sa01_sph;
extern CINTIntegralFunction int1e_giao_sa10sa01_spinor;

/* <R CROSS SIGMA i|OVLP |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_giao_sa10sp_optimizer;
extern CINTIntegralFunction int1e_giao_sa10sp_cart;
extern CINTIntegralFunction int1e_giao_sa10sp_sph;
extern CINTIntegralFunction int1e_giao_sa10sp_spinor;

/* <R CROSS SIGMA i|NUC |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_giao_sa10nucsp_optimizer;
extern CINTIntegralFunction int1e_giao_sa10nucsp_cart;
extern CINTIntegralFunction int1e_giao_sa10nucsp_sph;
extern CINTIntegralFunction int1e_giao_sa10nucsp_spinor;

/* <i|NABLA-RINV CROSS SIGMA |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_sa01sp_optimizer;
extern CINTIntegralFunction int1e_sa01sp_cart;
extern CINTIntegralFunction int1e_sa01sp_sph;
extern CINTIntegralFunction int1e_sa01sp_spinor;

/* <G SIGMA DOT P i|OVLP |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_spgsp_optimizer;
extern CINTIntegralFunction int1e_spgsp_cart;
extern CINTIntegralFunction int1e_spgsp_sph;
extern CINTIntegralFunction int1e_spgsp_spinor;

/* <G SIGMA DOT P i|NUC |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_spgnucsp_optimizer;
extern CINTIntegralFunction int1e_spgnucsp_cart;
extern CINTIntegralFunction int1e_spgnucsp_sph;
extern CINTIntegralFunction int1e_spgnucsp_spinor;

/* <G SIGMA DOT P i|NABLA-RINV CROSS SIGMA |j> */
extern CINTOptimizerFunction int1e_spgsa01_optimizer;
extern CINTIntegralFunction int1e_spgsa01_cart;
extern CINTIntegralFunction int1e_spgsa01_sph;
extern CINTIntegralFunction int1e_spgsa01_spinor;

/* (SIGMA DOT P i SIGMA DOT P j|R12 |k l) */
extern CINTOptimizerFunction int2e_spsp1_optimizer;
extern CINTIntegralFunction int2e_spsp1_cart;
extern CINTIntegralFunction int2e_spsp1_sph;
extern CINTIntegralFunction int2e_spsp1_spinor;

/* (SIGMA DOT P i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_spsp1spsp2_optimizer;
extern CINTIntegralFunction int2e_spsp1spsp2_cart;
extern CINTIntegralFunction int2e_spsp1spsp2_sph;
extern CINTIntegralFunction int2e_spsp1spsp2_spinor;

/* (SIGMA DOT R i SIGMA DOT R j|R12 |k l) */
extern CINTOptimizerFunction int2e_srsr1_optimizer;
extern CINTIntegralFunction int2e_srsr1_cart;
extern CINTIntegralFunction int2e_srsr1_sph;
extern CINTIntegralFunction int2e_srsr1_spinor;

/* (SIGMA DOT R i SIGMA DOT R j|R12 |SIGMA DOT R k SIGMA DOT R l) */
extern CINTOptimizerFunction int2e_srsr1srsr2_optimizer;
extern CINTIntegralFunction int2e_srsr1srsr2_cart;
extern CINTIntegralFunction int2e_srsr1srsr2_sph;
extern CINTIntegralFunction int2e_srsr1srsr2_spinor;

/* (RC CROSS SIGMA i SIGMA DOT P j|R12 |k l) */
extern CINTOptimizerFunction int2e_cg_sa10sp1_optimizer;
extern CINTIntegralFunction int2e_cg_sa10sp1_cart;
extern CINTIntegralFunction int2e_cg_sa10sp1_sph;
extern CINTIntegralFunction int2e_cg_sa10sp1_spinor;

/* (RC CROSS SIGMA i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_cg_sa10sp1spsp2_optimizer;
extern CINTIntegralFunction int2e_cg_sa10sp1spsp2_cart;
extern CINTIntegralFunction int2e_cg_sa10sp1spsp2_sph;
extern CINTIntegralFunction int2e_cg_sa10sp1spsp2_spinor;

/* (R CROSS SIGMA i SIGMA DOT P j|R12 |k l) */
extern CINTOptimizerFunction int2e_giao_sa10sp1_optimizer;
extern CINTIntegralFunction int2e_giao_sa10sp1_cart;
extern CINTIntegralFunction int2e_giao_sa10sp1_sph;
extern CINTIntegralFunction int2e_giao_sa10sp1_spinor;

/* (R CROSS SIGMA i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_giao_sa10sp1spsp2_optimizer;
extern CINTIntegralFunction int2e_giao_sa10sp1spsp2_cart;
extern CINTIntegralFunction int2e_giao_sa10sp1spsp2_sph;
extern CINTIntegralFunction int2e_giao_sa10sp1spsp2_spinor;

/* (G i j|R12 |k l) */
extern CINTOptimizerFunction int2e_g1_optimizer;
extern CINTIntegralFunction int2e_g1_cart;
extern CINTIntegralFunction int2e_g1_sph;
extern CINTIntegralFunction int2e_g1_spinor;

/* (G SIGMA DOT P i SIGMA DOT P j|R12 |k l) */
extern CINTOptimizerFunction int2e_spgsp1_optimizer;
extern CINTIntegralFunction int2e_spgsp1_cart;
extern CINTIntegralFunction int2e_spgsp1_sph;
extern CINTIntegralFunction int2e_spgsp1_spinor;

/* (G i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_g1spsp2_optimizer;
extern CINTIntegralFunction int2e_g1spsp2_cart;
extern CINTIntegralFunction int2e_g1spsp2_sph;
extern CINTIntegralFunction int2e_g1spsp2_spinor;

/* (G SIGMA DOT P i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_spgsp1spsp2_optimizer;
extern CINTIntegralFunction int2e_spgsp1spsp2_cart;
extern CINTIntegralFunction int2e_spgsp1spsp2_sph;
extern CINTIntegralFunction int2e_spgsp1spsp2_spinor;

/* (P* i DOT P j|R12 |k l) */
extern CINTOptimizerFunction int2e_pp1_optimizer;
extern CINTIntegralFunction int2e_pp1_cart;
extern CINTIntegralFunction int2e_pp1_sph;
extern CINTIntegralFunction int2e_pp1_spinor;

/* (i j|R12 |P* k DOT P l) */
extern CINTOptimizerFunction int2e_pp2_optimizer;
extern CINTIntegralFunction int2e_pp2_cart;
extern CINTIntegralFunction int2e_pp2_sph;
extern CINTIntegralFunction int2e_pp2_spinor;

/* (P* i DOT P j|R12 |P* k DOT P l) */
extern CINTOptimizerFunction int2e_pp1pp2_optimizer;
extern CINTIntegralFunction int2e_pp1pp2_cart;
extern CINTIntegralFunction int2e_pp1pp2_sph;
extern CINTIntegralFunction int2e_pp1pp2_spinor;

/* <SIGMA DOT P i|OVLP |SIGMA DOT P SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_spspsp_optimizer;
extern CINTIntegralFunction int1e_spspsp_cart;
extern CINTIntegralFunction int1e_spspsp_sph;
extern CINTIntegralFunction int1e_spspsp_spinor;

/* <SIGMA DOT P i|NUC |j> */
extern CINTOptimizerFunction int1e_spnuc_optimizer;
extern CINTIntegralFunction int1e_spnuc_cart;
extern CINTIntegralFunction int1e_spnuc_sph;
extern CINTIntegralFunction int1e_spnuc_spinor;

/* (SIGMA DOT P i j|R12 |k l) */
extern CINTOptimizerFunction int2e_spv1_optimizer;
extern CINTIntegralFunction int2e_spv1_cart;
extern CINTIntegralFunction int2e_spv1_sph;
extern CINTIntegralFunction int2e_spv1_spinor;

/* (i SIGMA DOT P j|R12 |k l) */
extern CINTOptimizerFunction int2e_vsp1_optimizer;
extern CINTIntegralFunction int2e_vsp1_cart;
extern CINTIntegralFunction int2e_vsp1_sph;
extern CINTIntegralFunction int2e_vsp1_spinor;

/* (i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_spsp2_optimizer;
extern CINTIntegralFunction int2e_spsp2_cart;
extern CINTIntegralFunction int2e_spsp2_sph;
extern CINTIntegralFunction int2e_spsp2_spinor;

/* (SIGMA DOT P i j|R12 |SIGMA DOT P k l) */
extern CINTOptimizerFunction int2e_spv1spv2_optimizer;
extern CINTIntegralFunction int2e_spv1spv2_cart;
extern CINTIntegralFunction int2e_spv1spv2_sph;
extern CINTIntegralFunction int2e_spv1spv2_spinor;

/* (i SIGMA DOT P j|R12 |SIGMA DOT P k l) */
extern CINTOptimizerFunction int2e_vsp1spv2_optimizer;
extern CINTIntegralFunction int2e_vsp1spv2_cart;
extern CINTIntegralFunction int2e_vsp1spv2_sph;
extern CINTIntegralFunction int2e_vsp1spv2_spinor;

/* (SIGMA DOT P i j|R12 |k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_spv1vsp2_optimizer;
extern CINTIntegralFunction int2e_spv1vsp2_cart;
extern CINTIntegralFunction int2e_spv1vsp2_sph;
extern CINTIntegralFunction int2e_spv1vsp2_spinor;

/* (i SIGMA DOT P j|R12 |k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_vsp1vsp2_optimizer;
extern CINTIntegralFunction int2e_vsp1vsp2_cart;
extern CINTIntegralFunction int2e_vsp1vsp2_sph;
extern CINTIntegralFunction int2e_vsp1vsp2_spinor;

/* (SIGMA DOT P i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_spv1spsp2_optimizer;
extern CINTIntegralFunction int2e_spv1spsp2_cart;
extern CINTIntegralFunction int2e_spv1spsp2_sph;
extern CINTIntegralFunction int2e_spv1spsp2_spinor;

/* (i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_vsp1spsp2_optimizer;
extern CINTIntegralFunction int2e_vsp1spsp2_cart;
extern CINTIntegralFunction int2e_vsp1spsp2_sph;
extern CINTIntegralFunction int2e_vsp1spsp2_spinor;

/* <NABLA i|OVLP |j> */
extern CINTOptimizerFunction int1e_ipovlp_optimizer;
extern CINTIntegralFunction int1e_ipovlp_cart;
extern CINTIntegralFunction int1e_ipovlp_sph;
extern CINTIntegralFunction int1e_ipovlp_spinor;

/* <i|OVLP |NABLA j> */
extern CINTOptimizerFunction int1e_ovlpip_optimizer;
extern CINTIntegralFunction int1e_ovlpip_cart;
extern CINTIntegralFunction int1e_ovlpip_sph;
extern CINTIntegralFunction int1e_ovlpip_spinor;

/* <NABLA i|OVLP |P DOT P j> */
extern CINTOptimizerFunction int1e_ipkin_optimizer;
extern CINTIntegralFunction int1e_ipkin_cart;
extern CINTIntegralFunction int1e_ipkin_sph;
extern CINTIntegralFunction int1e_ipkin_spinor;

/* <i|OVLP |P DOT P NABLA j> */
extern CINTOptimizerFunction int1e_kinip_optimizer;
extern CINTIntegralFunction int1e_kinip_cart;
extern CINTIntegralFunction int1e_kinip_sph;
extern CINTIntegralFunction int1e_kinip_spinor;

/* <NABLA i|NUC |j> */
extern CINTOptimizerFunction int1e_ipnuc_optimizer;
extern CINTIntegralFunction int1e_ipnuc_cart;
extern CINTIntegralFunction int1e_ipnuc_sph;
extern CINTIntegralFunction int1e_ipnuc_spinor;

/* <NABLA i|RINV |j> */
extern CINTOptimizerFunction int1e_iprinv_optimizer;
extern CINTIntegralFunction int1e_iprinv_cart;
extern CINTIntegralFunction int1e_iprinv_sph;
extern CINTIntegralFunction int1e_iprinv_spinor;

/* <NABLA SIGMA DOT P i|NUC |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_ipspnucsp_optimizer;
extern CINTIntegralFunction int1e_ipspnucsp_cart;
extern CINTIntegralFunction int1e_ipspnucsp_sph;
extern CINTIntegralFunction int1e_ipspnucsp_spinor;

/* <NABLA SIGMA DOT P i|RINV |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_ipsprinvsp_optimizer;
extern CINTIntegralFunction int1e_ipsprinvsp_cart;
extern CINTIntegralFunction int1e_ipsprinvsp_sph;
extern CINTIntegralFunction int1e_ipsprinvsp_spinor;

/* <P* NABLA i|NUC DOT P |j> */
extern CINTOptimizerFunction int1e_ippnucp_optimizer;
extern CINTIntegralFunction int1e_ippnucp_cart;
extern CINTIntegralFunction int1e_ippnucp_sph;
extern CINTIntegralFunction int1e_ippnucp_spinor;

/* <P* NABLA i|RINV DOT P |j> */
extern CINTOptimizerFunction int1e_ipprinvp_optimizer;
extern CINTIntegralFunction int1e_ipprinvp_cart;
extern CINTIntegralFunction int1e_ipprinvp_sph;
extern CINTIntegralFunction int1e_ipprinvp_spinor;

/* (NABLA i j|R12 |k l) */
extern CINTOptimizerFunction int2e_ip1_optimizer;
extern CINTIntegralFunction int2e_ip1_cart;
extern CINTIntegralFunction int2e_ip1_sph;
extern CINTIntegralFunction int2e_ip1_spinor;

/* (i j|R12 |NABLA k l) */
extern CINTOptimizerFunction int2e_ip2_optimizer;
extern CINTIntegralFunction int2e_ip2_cart;
extern CINTIntegralFunction int2e_ip2_sph;
extern CINTIntegralFunction int2e_ip2_spinor;

/* (NABLA SIGMA DOT P i SIGMA DOT P j|R12 |k l) */
extern CINTOptimizerFunction int2e_ipspsp1_optimizer;
extern CINTIntegralFunction int2e_ipspsp1_cart;
extern CINTIntegralFunction int2e_ipspsp1_sph;
extern CINTIntegralFunction int2e_ipspsp1_spinor;

/* (NABLA i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_ip1spsp2_optimizer;
extern CINTIntegralFunction int2e_ip1spsp2_cart;
extern CINTIntegralFunction int2e_ip1spsp2_sph;
extern CINTIntegralFunction int2e_ip1spsp2_spinor;

/* (NABLA SIGMA DOT P i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_ipspsp1spsp2_optimizer;
extern CINTIntegralFunction int2e_ipspsp1spsp2_cart;
extern CINTIntegralFunction int2e_ipspsp1spsp2_sph;
extern CINTIntegralFunction int2e_ipspsp1spsp2_spinor;

/* (NABLA SIGMA DOT R i SIGMA DOT R j|R12 |k l) */
extern CINTOptimizerFunction int2e_ipsrsr1_optimizer;
extern CINTIntegralFunction int2e_ipsrsr1_cart;
extern CINTIntegralFunction int2e_ipsrsr1_sph;
extern CINTIntegralFunction int2e_ipsrsr1_spinor;

/* (NABLA i j|R12 |SIGMA DOT R k SIGMA DOT R l) */
extern CINTOptimizerFunction int2e_ip1srsr2_optimizer;
extern CINTIntegralFunction int2e_ip1srsr2_cart;
extern CINTIntegralFunction int2e_ip1srsr2_sph;
extern CINTIntegralFunction int2e_ip1srsr2_spinor;

/* (NABLA SIGMA DOT R i SIGMA DOT R j|R12 |SIGMA DOT R k SIGMA DOT R l) */
extern CINTOptimizerFunction int2e_ipsrsr1srsr2_optimizer;
extern CINTIntegralFunction int2e_ipsrsr1srsr2_cart;
extern CINTIntegralFunction int2e_ipsrsr1srsr2_sph;
extern CINTIntegralFunction int2e_ipsrsr1srsr2_spinor;

/* (i SIGMA DOT P j|GAUNT |k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_ssp1ssp2_optimizer;
extern CINTIntegralFunction int2e_ssp1ssp2_cart;
extern CINTIntegralFunction int2e_ssp1ssp2_sph;
extern CINTIntegralFunction int2e_ssp1ssp2_spinor;

/* (i SIGMA DOT P j|GAUNT |SIGMA DOT P k l) */
extern CINTOptimizerFunction int2e_ssp1sps2_optimizer;
extern CINTIntegralFunction int2e_ssp1sps2_cart;
extern CINTIntegralFunction int2e_ssp1sps2_sph;
extern CINTIntegralFunction int2e_ssp1sps2_spinor;

/* (SIGMA DOT P i j|GAUNT |k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_sps1ssp2_optimizer;
extern CINTIntegralFunction int2e_sps1ssp2_cart;
extern CINTIntegralFunction int2e_sps1ssp2_sph;
extern CINTIntegralFunction int2e_sps1ssp2_spinor;

/* (SIGMA DOT P i j|GAUNT |SIGMA DOT P k l) */
extern CINTOptimizerFunction int2e_sps1sps2_optimizer;
extern CINTIntegralFunction int2e_sps1sps2_cart;
extern CINTIntegralFunction int2e_sps1sps2_sph;
extern CINTIntegralFunction int2e_sps1sps2_spinor;

/* (RC CROSS SIGMA i j|GAUNT |k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_cg_ssa10ssp2_optimizer;
extern CINTIntegralFunction int2e_cg_ssa10ssp2_cart;
extern CINTIntegralFunction int2e_cg_ssa10ssp2_sph;
extern CINTIntegralFunction int2e_cg_ssa10ssp2_spinor;

/* (R CROSS SIGMA i j|GAUNT |k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_giao_ssa10ssp2_optimizer;
extern CINTIntegralFunction int2e_giao_ssa10ssp2_cart;
extern CINTIntegralFunction int2e_giao_ssa10ssp2_sph;
extern CINTIntegralFunction int2e_giao_ssa10ssp2_spinor;

/* (G i SIGMA DOT P j|GAUNT |k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_gssp1ssp2_optimizer;
extern CINTIntegralFunction int2e_gssp1ssp2_cart;
extern CINTIntegralFunction int2e_gssp1ssp2_sph;
extern CINTIntegralFunction int2e_gssp1ssp2_spinor;

/* (i R0 SIGMA DOT P j|BREIT-R1 |k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_gauge_r1_ssp1ssp2_optimizer;
extern CINTIntegralFunction int2e_gauge_r1_ssp1ssp2_cart;
extern CINTIntegralFunction int2e_gauge_r1_ssp1ssp2_sph;
extern CINTIntegralFunction int2e_gauge_r1_ssp1ssp2_spinor;

/* (i R0 SIGMA DOT P j|BREIT-R1 |SIGMA DOT P k l) */
extern CINTOptimizerFunction int2e_gauge_r1_ssp1sps2_optimizer;
extern CINTIntegralFunction int2e_gauge_r1_ssp1sps2_cart;
extern CINTIntegralFunction int2e_gauge_r1_ssp1sps2_sph;
extern CINTIntegralFunction int2e_gauge_r1_ssp1sps2_spinor;

/* (SIGMA DOT P i R0 j|BREIT-R1 |k SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_gauge_r1_sps1ssp2_optimizer;
extern CINTIntegralFunction int2e_gauge_r1_sps1ssp2_cart;
extern CINTIntegralFunction int2e_gauge_r1_sps1ssp2_sph;
extern CINTIntegralFunction int2e_gauge_r1_sps1ssp2_spinor;

/* (SIGMA DOT P i R0 j|BREIT-R1 |SIGMA DOT P k l) */
extern CINTOptimizerFunction int2e_gauge_r1_sps1sps2_optimizer;
extern CINTIntegralFunction int2e_gauge_r1_sps1sps2_cart;
extern CINTIntegralFunction int2e_gauge_r1_sps1sps2_sph;
extern CINTIntegralFunction int2e_gauge_r1_sps1sps2_spinor;

/* (i SIGMA DOT P j|BREIT-R2 |k R0 SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_gauge_r2_ssp1ssp2_optimizer;
extern CINTIntegralFunction int2e_gauge_r2_ssp1ssp2_cart;
extern CINTIntegralFunction int2e_gauge_r2_ssp1ssp2_sph;
extern CINTIntegralFunction int2e_gauge_r2_ssp1ssp2_spinor;

/* (i SIGMA DOT P j|BREIT-R2 |SIGMA DOT P k R0 l) */
extern CINTOptimizerFunction int2e_gauge_r2_ssp1sps2_optimizer;
extern CINTIntegralFunction int2e_gauge_r2_ssp1sps2_cart;
extern CINTIntegralFunction int2e_gauge_r2_ssp1sps2_sph;
extern CINTIntegralFunction int2e_gauge_r2_ssp1sps2_spinor;

/* (SIGMA DOT P i j|BREIT-R2 |k R0 SIGMA DOT P l) */
extern CINTOptimizerFunction int2e_gauge_r2_sps1ssp2_optimizer;
extern CINTIntegralFunction int2e_gauge_r2_sps1ssp2_cart;
extern CINTIntegralFunction int2e_gauge_r2_sps1ssp2_sph;
extern CINTIntegralFunction int2e_gauge_r2_sps1ssp2_spinor;

/* (SIGMA DOT P i j|BREIT-R2 |SIGMA DOT P k R0 l) */
extern CINTOptimizerFunction int2e_gauge_r2_sps1sps2_optimizer;
extern CINTIntegralFunction int2e_gauge_r2_sps1sps2_cart;
extern CINTIntegralFunction int2e_gauge_r2_sps1sps2_sph;
extern CINTIntegralFunction int2e_gauge_r2_sps1sps2_spinor;

/* <NABLA NABLA i|OVLP |j> */
extern CINTOptimizerFunction int1e_ipipovlp_optimizer;
extern CINTIntegralFunction int1e_ipipovlp_cart;
extern CINTIntegralFunction int1e_ipipovlp_sph;
extern CINTIntegralFunction int1e_ipipovlp_spinor;

/* <NABLA i|OVLP |NABLA j> */
extern CINTOptimizerFunction int1e_ipovlpip_optimizer;
extern CINTIntegralFunction int1e_ipovlpip_cart;
extern CINTIntegralFunction int1e_ipovlpip_sph;
extern CINTIntegralFunction int1e_ipovlpip_spinor;

/* <NABLA NABLA i|P DOT P |j> */
extern CINTOptimizerFunction int1e_ipipkin_optimizer;
extern CINTIntegralFunction int1e_ipipkin_cart;
extern CINTIntegralFunction int1e_ipipkin_sph;
extern CINTIntegralFunction int1e_ipipkin_spinor;

/* <NABLA i|P DOT P |NABLA j> */
extern CINTOptimizerFunction int1e_ipkinip_optimizer;
extern CINTIntegralFunction int1e_ipkinip_cart;
extern CINTIntegralFunction int1e_ipkinip_sph;
extern CINTIntegralFunction int1e_ipkinip_spinor;

/* <NABLA NABLA i|NUC |j> */
extern CINTOptimizerFunction int1e_ipipnuc_optimizer;
extern CINTIntegralFunction int1e_ipipnuc_cart;
extern CINTIntegralFunction int1e_ipipnuc_sph;
extern CINTIntegralFunction int1e_ipipnuc_spinor;

/* <NABLA i|NUC |NABLA j> */
extern CINTOptimizerFunction int1e_ipnucip_optimizer;
extern CINTIntegralFunction int1e_ipnucip_cart;
extern CINTIntegralFunction int1e_ipnucip_sph;
extern CINTIntegralFunction int1e_ipnucip_spinor;

/* <NABLA NABLA i|RINV |j> */
extern CINTOptimizerFunction int1e_ipiprinv_optimizer;
extern CINTIntegralFunction int1e_ipiprinv_cart;
extern CINTIntegralFunction int1e_ipiprinv_sph;
extern CINTIntegralFunction int1e_ipiprinv_spinor;

/* <NABLA i|RINV |NABLA j> */
extern CINTOptimizerFunction int1e_iprinvip_optimizer;
extern CINTIntegralFunction int1e_iprinvip_cart;
extern CINTIntegralFunction int1e_iprinvip_sph;
extern CINTIntegralFunction int1e_iprinvip_spinor;

/* <NABLA NABLA i|RC |j> */
extern CINTOptimizerFunction int1e_ipipr_optimizer;
extern CINTIntegralFunction int1e_ipipr_cart;
extern CINTIntegralFunction int1e_ipipr_sph;
extern CINTIntegralFunction int1e_ipipr_spinor;

/* <NABLA i|RC |NABLA j> */
extern CINTOptimizerFunction int1e_iprip_optimizer;
extern CINTIntegralFunction int1e_iprip_cart;
extern CINTIntegralFunction int1e_iprip_sph;
extern CINTIntegralFunction int1e_iprip_spinor;

/* (NABLA NABLA i j|R12 |k l) */
extern CINTOptimizerFunction int2e_ipip1_optimizer;
extern CINTIntegralFunction int2e_ipip1_cart;
extern CINTIntegralFunction int2e_ipip1_sph;
extern CINTIntegralFunction int2e_ipip1_spinor;

/* (NABLA i NABLA j|R12 |k l) */
extern CINTOptimizerFunction int2e_ipvip1_optimizer;
extern CINTIntegralFunction int2e_ipvip1_cart;
extern CINTIntegralFunction int2e_ipvip1_sph;
extern CINTIntegralFunction int2e_ipvip1_spinor;

/* (NABLA i j|R12 |NABLA k l) */
extern CINTOptimizerFunction int2e_ip1ip2_optimizer;
extern CINTIntegralFunction int2e_ip1ip2_cart;
extern CINTIntegralFunction int2e_ip1ip2_sph;
extern CINTIntegralFunction int2e_ip1ip2_spinor;

/* <P* NABLA NABLA i|NUC DOT P |j> */
extern CINTOptimizerFunction int1e_ipippnucp_optimizer;
extern CINTIntegralFunction int1e_ipippnucp_cart;
extern CINTIntegralFunction int1e_ipippnucp_sph;
extern CINTIntegralFunction int1e_ipippnucp_spinor;

/* <P* NABLA i|NUC DOT P |NABLA j> */
extern CINTOptimizerFunction int1e_ippnucpip_optimizer;
extern CINTIntegralFunction int1e_ippnucpip_cart;
extern CINTIntegralFunction int1e_ippnucpip_sph;
extern CINTIntegralFunction int1e_ippnucpip_spinor;

/* <P* NABLA NABLA i|RINV DOT P |j> */
extern CINTOptimizerFunction int1e_ipipprinvp_optimizer;
extern CINTIntegralFunction int1e_ipipprinvp_cart;
extern CINTIntegralFunction int1e_ipipprinvp_sph;
extern CINTIntegralFunction int1e_ipipprinvp_spinor;

/* <P* NABLA i|RINV DOT P |NABLA j> */
extern CINTOptimizerFunction int1e_ipprinvpip_optimizer;
extern CINTIntegralFunction int1e_ipprinvpip_cart;
extern CINTIntegralFunction int1e_ipprinvpip_sph;
extern CINTIntegralFunction int1e_ipprinvpip_spinor;

/* <NABLA NABLA SIGMA DOT P i|NUC SIGMA DOT P |j> */
extern CINTOptimizerFunction int1e_ipipspnucsp_optimizer;
extern CINTIntegralFunction int1e_ipipspnucsp_cart;
extern CINTIntegralFunction int1e_ipipspnucsp_sph;
extern CINTIntegralFunction int1e_ipipspnucsp_spinor;

/* <NABLA SIGMA DOT P i|NUC SIGMA DOT P |NABLA j> */
extern CINTOptimizerFunction int1e_ipspnucspip_optimizer;
extern CINTIntegralFunction int1e_ipspnucspip_cart;
extern CINTIntegralFunction int1e_ipspnucspip_sph;
extern CINTIntegralFunction int1e_ipspnucspip_spinor;

/* <NABLA NABLA SIGMA DOT P i|RINV SIGMA DOT P |j> */
extern CINTOptimizerFunction int1e_ipipsprinvsp_optimizer;
extern CINTIntegralFunction int1e_ipipsprinvsp_cart;
extern CINTIntegralFunction int1e_ipipsprinvsp_sph;
extern CINTIntegralFunction int1e_ipipsprinvsp_spinor;

/* <NABLA SIGMA DOT P i|RINV SIGMA DOT P |NABLA j> */
extern CINTOptimizerFunction int1e_ipsprinvspip_optimizer;
extern CINTIntegralFunction int1e_ipsprinvspip_cart;
extern CINTIntegralFunction int1e_ipsprinvspip_sph;
extern CINTIntegralFunction int1e_ipsprinvspip_spinor;

/* (NABLA NABLA i j|R12 |NABLA NABLA k l) */
extern CINTOptimizerFunction int2e_ipip1ipip2_optimizer;
extern CINTIntegralFunction int2e_ipip1ipip2_cart;
extern CINTIntegralFunction int2e_ipip1ipip2_sph;
extern CINTIntegralFunction int2e_ipip1ipip2_spinor;

/* (NABLA i NABLA j|R12 |NABLA k NABLA l) */
extern CINTOptimizerFunction int2e_ipvip1ipvip2_optimizer;
extern CINTIntegralFunction int2e_ipvip1ipvip2_cart;
extern CINTIntegralFunction int2e_ipvip1ipvip2_sph;
extern CINTIntegralFunction int2e_ipvip1ipvip2_spinor;

/* (NABLA i j|R12 |k) */
extern CINTOptimizerFunction int3c2e_ip1_optimizer;
extern CINTIntegralFunction int3c2e_ip1_cart;
extern CINTIntegralFunction int3c2e_ip1_sph;
extern CINTIntegralFunction int3c2e_ip1_spinor;

/* (i j|R12 |NABLA k) */
extern CINTOptimizerFunction int3c2e_ip2_optimizer;
extern CINTIntegralFunction int3c2e_ip2_cart;
extern CINTIntegralFunction int3c2e_ip2_sph;
extern CINTIntegralFunction int3c2e_ip2_spinor;

/* (P* i DOT P j|R12 |k) */
extern CINTOptimizerFunction int3c2e_pvp1_optimizer;
extern CINTIntegralFunction int3c2e_pvp1_cart;
extern CINTIntegralFunction int3c2e_pvp1_sph;
extern CINTIntegralFunction int3c2e_pvp1_spinor;

/* (P* i CROSS P j|R12 |k) */
extern CINTOptimizerFunction int3c2e_pvxp1_optimizer;
extern CINTIntegralFunction int3c2e_pvxp1_cart;
extern CINTIntegralFunction int3c2e_pvxp1_sph;
extern CINTIntegralFunction int3c2e_pvxp1_spinor;

/* (NABLA i |R12 |j) */
extern CINTOptimizerFunction int2c2e_ip1_optimizer;
extern CINTIntegralFunction int2c2e_ip1_cart;
extern CINTIntegralFunction int2c2e_ip1_sph;
extern CINTIntegralFunction int2c2e_ip1_spinor;

/* (i |R12 |NABLA j) */
extern CINTOptimizerFunction int2c2e_ip2_optimizer;
extern CINTIntegralFunction int2c2e_ip2_cart;
extern CINTIntegralFunction int2c2e_ip2_sph;
extern CINTIntegralFunction int2c2e_ip2_spinor;

/* (G i j|R12 |k) */
extern CINTOptimizerFunction int3c2e_ig1_optimizer;
extern CINTIntegralFunction int3c2e_ig1_cart;
extern CINTIntegralFunction int3c2e_ig1_sph;
extern CINTIntegralFunction int3c2e_ig1_spinor;

/* (SIGMA DOT P i SIGMA DOT P j|R12 |k) */
extern CINTOptimizerFunction int3c2e_spsp1_optimizer;
extern CINTIntegralFunction int3c2e_spsp1_cart;
extern CINTIntegralFunction int3c2e_spsp1_sph;
extern CINTIntegralFunction int3c2e_spsp1_spinor;

/* (NABLA SIGMA DOT P i SIGMA DOT P j|R12 |k) */
extern CINTOptimizerFunction int3c2e_ipspsp1_optimizer;
extern CINTIntegralFunction int3c2e_ipspsp1_cart;
extern CINTIntegralFunction int3c2e_ipspsp1_sph;
extern CINTIntegralFunction int3c2e_ipspsp1_spinor;

/* (SIGMA DOT P i SIGMA DOT P j|R12 |NABLA k) */
extern CINTOptimizerFunction int3c2e_spsp1ip2_optimizer;
extern CINTIntegralFunction int3c2e_spsp1ip2_cart;
extern CINTIntegralFunction int3c2e_spsp1ip2_sph;
extern CINTIntegralFunction int3c2e_spsp1ip2_spinor;

/* (NABLA NABLA i j|R12 |k) */
extern CINTOptimizerFunction int3c2e_ipip1_optimizer;
extern CINTIntegralFunction int3c2e_ipip1_cart;
extern CINTIntegralFunction int3c2e_ipip1_sph;
extern CINTIntegralFunction int3c2e_ipip1_spinor;

/* (i j|R12 |NABLA NABLA k) */
extern CINTOptimizerFunction int3c2e_ipip2_optimizer;
extern CINTIntegralFunction int3c2e_ipip2_cart;
extern CINTIntegralFunction int3c2e_ipip2_sph;
extern CINTIntegralFunction int3c2e_ipip2_spinor;

/* (NABLA i NABLA j|R12 |k) */
extern CINTOptimizerFunction int3c2e_ipvip1_optimizer;
extern CINTIntegralFunction int3c2e_ipvip1_cart;
extern CINTIntegralFunction int3c2e_ipvip1_sph;
extern CINTIntegralFunction int3c2e_ipvip1_spinor;

/* (NABLA i j|R12 |NABLA k) */
extern CINTOptimizerFunction int3c2e_ip1ip2_optimizer;
extern CINTIntegralFunction int3c2e_ip1ip2_cart;
extern CINTIntegralFunction int3c2e_ip1ip2_sph;
extern CINTIntegralFunction int3c2e_ip1ip2_spinor;

/* (NABLA NABLA i |R12 |j) */
extern CINTOptimizerFunction int2c2e_ipip1_optimizer;
extern CINTIntegralFunction int2c2e_ipip1_cart;
extern CINTIntegralFunction int2c2e_ipip1_sph;
extern CINTIntegralFunction int2c2e_ipip1_spinor;

/* (NABLA i |R12 |NABLA j) */
extern CINTOptimizerFunction int2c2e_ip1ip2_optimizer;
extern CINTIntegralFunction int2c2e_ip1ip2_cart;
extern CINTIntegralFunction int2c2e_ip1ip2_sph;
extern CINTIntegralFunction int2c2e_ip1ip2_spinor;

/* 3-center 1-electron integral <(i) (j) (P DOT P k)> */
extern CINTOptimizerFunction int3c1e_p2_optimizer;
extern CINTIntegralFunction int3c1e_p2_cart;
extern CINTIntegralFunction int3c1e_p2_sph;
extern CINTIntegralFunction int3c1e_p2_spinor;

/* 3-center 1-electron integral <(P i) (j) (k)> */
extern CINTOptimizerFunction int3c1e_iprinv_optimizer;
extern CINTIntegralFunction int3c1e_iprinv_cart;
extern CINTIntegralFunction int3c1e_iprinv_sph;
extern CINTIntegralFunction int3c1e_iprinv_spinor;

/* 3-center 1-electron integral <(NABLA i) (j) (k)> */
extern CINTOptimizerFunction int3c1e_ip1_optimizer;
extern CINTIntegralFunction int3c1e_ip1_cart;
extern CINTIntegralFunction int3c1e_ip1_sph;
extern CINTIntegralFunction int3c1e_ip1_spinor;

/* <NABLA NABLA NABLA i|NUC |j> */
extern CINTOptimizerFunction int1e_ipipipnuc_optimizer;
extern CINTIntegralFunction int1e_ipipipnuc_cart;
extern CINTIntegralFunction int1e_ipipipnuc_sph;
extern CINTIntegralFunction int1e_ipipipnuc_spinor;

/* <NABLA NABLA NABLA i|RINV |j> */
extern CINTOptimizerFunction int1e_ipipiprinv_optimizer;
extern CINTIntegralFunction int1e_ipipiprinv_cart;
extern CINTIntegralFunction int1e_ipipiprinv_sph;
extern CINTIntegralFunction int1e_ipipiprinv_spinor;

/* <NABLA NABLA i|NUC |NABLA j> */
extern CINTOptimizerFunction int1e_ipipnucip_optimizer;
extern CINTIntegralFunction int1e_ipipnucip_cart;
extern CINTIntegralFunction int1e_ipipnucip_sph;
extern CINTIntegralFunction int1e_ipipnucip_spinor;

/* <NABLA NABLA i|RINV |NABLA j> */
extern CINTOptimizerFunction int1e_ipiprinvip_optimizer;
extern CINTIntegralFunction int1e_ipiprinvip_cart;
extern CINTIntegralFunction int1e_ipiprinvip_sph;
extern CINTIntegralFunction int1e_ipiprinvip_spinor;

/* <NABLA i| 1/r_{grids} |j> */
extern CINTOptimizerFunction int1e_grids_ip_optimizer;
extern CINTIntegralFunction int1e_grids_ip_cart;
extern CINTIntegralFunction int1e_grids_ip_sph;
extern CINTIntegralFunction int1e_grids_ip_spinor;

/* <NABLA i| 1/r_{grids} |NABLA j> */
extern CINTOptimizerFunction int1e_grids_ipvip_optimizer;
extern CINTIntegralFunction int1e_grids_ipvip_cart;
extern CINTIntegralFunction int1e_grids_ipvip_sph;
extern CINTIntegralFunction int1e_grids_ipvip_spinor;

/* <SIGMA DOT P i| 1/r_{grids} |SIGMA DOT P j> */
extern CINTOptimizerFunction int1e_grids_spvsp_optimizer;
extern CINTIntegralFunction int1e_grids_spvsp_cart;
extern CINTIntegralFunction int1e_grids_spvsp_sph;
extern CINTIntegralFunction int1e_grids_spvsp_spinor;

/* <NABLA NABLA i|RINV |NABLA NABLA j> */
extern CINTOptimizerFunction int1e_ipiprinvipip_optimizer;
extern CINTIntegralFunction int1e_ipiprinvipip_cart;
extern CINTIntegralFunction int1e_ipiprinvipip_sph;
extern CINTIntegralFunction int1e_ipiprinvipip_spinor;

/* <NABLA NABLA NABLA i|RINV |NABLA j> */
extern CINTOptimizerFunction int1e_ipipiprinvip_optimizer;
extern CINTIntegralFunction int1e_ipipiprinvip_cart;
extern CINTIntegralFunction int1e_ipipiprinvip_sph;
extern CINTIntegralFunction int1e_ipipiprinvip_spinor;

/* <NABLA NABLA NABLA NABLA i|RINV |j> */
extern CINTOptimizerFunction int1e_ipipipiprinv_optimizer;
extern CINTIntegralFunction int1e_ipipipiprinv_cart;
extern CINTIntegralFunction int1e_ipipipiprinv_sph;
extern CINTIntegralFunction int1e_ipipipiprinv_spinor;
