/*
 * Copyright (C) 2019-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Function signature
 */


#include <cint.h>

typedef void (*CINTOptimizerFunction)(CINTOpt **opt,
                                      FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env);
typedef FINT (*CINTIntegralFunction)(double *out, FINT *dims, FINT *shls,
                                     FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
                                     CINTOpt *opt, double *cache);
/* Plain ERI (ij|kl) */
CINTOptimizerFunction int2e_optimizer;
CINTIntegralFunction int2e_cart;
CINTIntegralFunction int2e_sph;
CINTIntegralFunction int2e_spinor;

/* <i|OVLP |j> */
CINTOptimizerFunction int1e_ovlp_optimizer;
CINTIntegralFunction int1e_ovlp_cart;
CINTIntegralFunction int1e_ovlp_sph;
CINTIntegralFunction int1e_ovlp_spinor;

/* <i|NUC |j> */
CINTOptimizerFunction int1e_nuc_optimizer;
CINTIntegralFunction int1e_nuc_cart;
CINTIntegralFunction int1e_nuc_sph;
CINTIntegralFunction int1e_nuc_spinor;

/* <i|OVLP |P DOT P j> */
CINTOptimizerFunction int1e_kin_optimizer;
CINTIntegralFunction int1e_kin_cart;
CINTIntegralFunction int1e_kin_sph;
CINTIntegralFunction int1e_kin_spinor;

/* <i|NABLA-RINV |CROSS P j> */
CINTOptimizerFunction int1e_ia01p_optimizer;
CINTIntegralFunction int1e_ia01p_cart;
CINTIntegralFunction int1e_ia01p_sph;
CINTIntegralFunction int1e_ia01p_spinor;

/* <i|OVLP |R CROSS P j> */
CINTOptimizerFunction int1e_giao_irjxp_optimizer;
CINTIntegralFunction int1e_giao_irjxp_cart;
CINTIntegralFunction int1e_giao_irjxp_sph;
CINTIntegralFunction int1e_giao_irjxp_spinor;

/* <i|OVLP |RC CROSS P j> */
CINTOptimizerFunction int1e_cg_irxp_optimizer;
CINTIntegralFunction int1e_cg_irxp_cart;
CINTIntegralFunction int1e_cg_irxp_sph;
CINTIntegralFunction int1e_cg_irxp_spinor;

/* <i|NABLA-RINV |R j> */
CINTOptimizerFunction int1e_giao_a11part_optimizer;
CINTIntegralFunction int1e_giao_a11part_cart;
CINTIntegralFunction int1e_giao_a11part_sph;
CINTIntegralFunction int1e_giao_a11part_spinor;

/* <i|NABLA-RINV |RC j> */
CINTOptimizerFunction int1e_cg_a11part_optimizer;
CINTIntegralFunction int1e_cg_a11part_cart;
CINTIntegralFunction int1e_cg_a11part_sph;
CINTIntegralFunction int1e_cg_a11part_spinor;

/* <G i|NABLA-RINV CROSS P |j> */
CINTOptimizerFunction int1e_a01gp_optimizer;
CINTIntegralFunction int1e_a01gp_cart;
CINTIntegralFunction int1e_a01gp_sph;
CINTIntegralFunction int1e_a01gp_spinor;

/* <G i|OVLP |P DOT P j> */
CINTOptimizerFunction int1e_igkin_optimizer;
CINTIntegralFunction int1e_igkin_cart;
CINTIntegralFunction int1e_igkin_sph;
CINTIntegralFunction int1e_igkin_spinor;

/* <G i|OVLP |j> */
CINTOptimizerFunction int1e_igovlp_optimizer;
CINTIntegralFunction int1e_igovlp_cart;
CINTIntegralFunction int1e_igovlp_sph;
CINTIntegralFunction int1e_igovlp_spinor;

/* <G i|NUC |j> */
CINTOptimizerFunction int1e_ignuc_optimizer;
CINTIntegralFunction int1e_ignuc_cart;
CINTIntegralFunction int1e_ignuc_sph;
CINTIntegralFunction int1e_ignuc_spinor;

/* <P* i|NUC DOT P |j> */
CINTOptimizerFunction int1e_pnucp_optimizer;
CINTIntegralFunction int1e_pnucp_cart;
CINTIntegralFunction int1e_pnucp_sph;
CINTIntegralFunction int1e_pnucp_spinor;

/* <i|ZC |j> */
CINTOptimizerFunction int1e_z_optimizer;
CINTIntegralFunction int1e_z_cart;
CINTIntegralFunction int1e_z_sph;
CINTIntegralFunction int1e_z_spinor;

/* <i|ZC ZC |j> */
CINTOptimizerFunction int1e_zz_optimizer;
CINTIntegralFunction int1e_zz_cart;
CINTIntegralFunction int1e_zz_sph;
CINTIntegralFunction int1e_zz_spinor;

/* <i|RC |j> */
CINTOptimizerFunction int1e_r_optimizer;
CINTIntegralFunction int1e_r_cart;
CINTIntegralFunction int1e_r_sph;
CINTIntegralFunction int1e_r_spinor;

/* <i|RC DOT RC |j> */
CINTOptimizerFunction int1e_r2_optimizer;
CINTIntegralFunction int1e_r2_cart;
CINTIntegralFunction int1e_r2_sph;
CINTIntegralFunction int1e_r2_spinor;

/* <i|RC RC |j> */
CINTOptimizerFunction int1e_rr_optimizer;
CINTIntegralFunction int1e_rr_cart;
CINTIntegralFunction int1e_rr_sph;
CINTIntegralFunction int1e_rr_spinor;

/* <i|RC RC RC |j> */
CINTOptimizerFunction int1e_rrr_optimizer;
CINTIntegralFunction int1e_rrr_cart;
CINTIntegralFunction int1e_rrr_sph;
CINTIntegralFunction int1e_rrr_spinor;

/* <i|RC RC RC RC |j> */
CINTOptimizerFunction int1e_rrrr_optimizer;
CINTIntegralFunction int1e_rrrr_cart;
CINTIntegralFunction int1e_rrrr_sph;
CINTIntegralFunction int1e_rrrr_spinor;

/* <i|Z |j> */
CINTOptimizerFunction int1e_z_origj_optimizer;
CINTIntegralFunction int1e_z_origj_cart;
CINTIntegralFunction int1e_z_origj_sph;
CINTIntegralFunction int1e_z_origj_spinor;

/* <i|Z Z |j> */
CINTOptimizerFunction int1e_zz_origj_optimizer;
CINTIntegralFunction int1e_zz_origj_cart;
CINTIntegralFunction int1e_zz_origj_sph;
CINTIntegralFunction int1e_zz_origj_spinor;

/* <i|R |j> */
CINTOptimizerFunction int1e_r_origj_optimizer;
CINTIntegralFunction int1e_r_origj_cart;
CINTIntegralFunction int1e_r_origj_sph;
CINTIntegralFunction int1e_r_origj_spinor;

/* <i|R R |j> */
CINTOptimizerFunction int1e_rr_origj_optimizer;
CINTIntegralFunction int1e_rr_origj_cart;
CINTIntegralFunction int1e_rr_origj_sph;
CINTIntegralFunction int1e_rr_origj_spinor;

/* <i|R DOT R |j> */
CINTOptimizerFunction int1e_r2_origj_optimizer;
CINTIntegralFunction int1e_r2_origj_cart;
CINTIntegralFunction int1e_r2_origj_sph;
CINTIntegralFunction int1e_r2_origj_spinor;

/* <i|OVLP |R DOT R R DOT R j> */
CINTOptimizerFunction int1e_r4_origj_optimizer;
CINTIntegralFunction int1e_r4_origj_cart;
CINTIntegralFunction int1e_r4_origj_sph;
CINTIntegralFunction int1e_r4_origj_spinor;

/* <P DOT P i|OVLP |P DOT P j> */
CINTOptimizerFunction int1e_p4_optimizer;
CINTIntegralFunction int1e_p4_cart;
CINTIntegralFunction int1e_p4_sph;
CINTIntegralFunction int1e_p4_spinor;

/* <P* i|RINV DOT P |j> */
CINTOptimizerFunction int1e_prinvp_optimizer;
CINTIntegralFunction int1e_prinvp_cart;
CINTIntegralFunction int1e_prinvp_sph;
CINTIntegralFunction int1e_prinvp_spinor;

/* <P* i|RINV CROSS P |j> */
CINTOptimizerFunction int1e_prinvxp_optimizer;
CINTIntegralFunction int1e_prinvxp_cart;
CINTIntegralFunction int1e_prinvxp_sph;
CINTIntegralFunction int1e_prinvxp_spinor;

/* <P* i|NUC CROSS P |j> */
CINTOptimizerFunction int1e_pnucxp_optimizer;
CINTIntegralFunction int1e_pnucxp_cart;
CINTIntegralFunction int1e_pnucxp_sph;
CINTIntegralFunction int1e_pnucxp_spinor;

/* <i|RC NABLA |j> */
CINTOptimizerFunction int1e_irp_optimizer;
CINTIntegralFunction int1e_irp_cart;
CINTIntegralFunction int1e_irp_sph;
CINTIntegralFunction int1e_irp_spinor;

/* <i|RC RC NABLA |j> */
CINTOptimizerFunction int1e_irrp_optimizer;
CINTIntegralFunction int1e_irrp_cart;
CINTIntegralFunction int1e_irrp_sph;
CINTIntegralFunction int1e_irrp_spinor;

/* <i|RC NABLA RC |j> */
CINTOptimizerFunction int1e_irpr_optimizer;
CINTIntegralFunction int1e_irpr_cart;
CINTIntegralFunction int1e_irpr_sph;
CINTIntegralFunction int1e_irpr_spinor;

/* <i|G G |j> */
CINTOptimizerFunction int1e_ggovlp_optimizer;
CINTIntegralFunction int1e_ggovlp_cart;
CINTIntegralFunction int1e_ggovlp_sph;
CINTIntegralFunction int1e_ggovlp_spinor;

/* <i|G G NABLA DOT NABLA |j> */
CINTOptimizerFunction int1e_ggkin_optimizer;
CINTIntegralFunction int1e_ggkin_cart;
CINTIntegralFunction int1e_ggkin_sph;
CINTIntegralFunction int1e_ggkin_spinor;

/* <i|G G NUC |j> */
CINTOptimizerFunction int1e_ggnuc_optimizer;
CINTIntegralFunction int1e_ggnuc_cart;
CINTIntegralFunction int1e_ggnuc_sph;
CINTIntegralFunction int1e_ggnuc_spinor;

/* <i|G R CROSS P |j> */
CINTOptimizerFunction int1e_grjxp_optimizer;
CINTIntegralFunction int1e_grjxp_cart;
CINTIntegralFunction int1e_grjxp_sph;
CINTIntegralFunction int1e_grjxp_spinor;

/* <i|RINV |j> */
CINTOptimizerFunction int1e_rinv_optimizer;
CINTIntegralFunction int1e_rinv_cart;
CINTIntegralFunction int1e_rinv_sph;
CINTIntegralFunction int1e_rinv_spinor;

/* <i|NABLA-RINV |j> */
CINTOptimizerFunction int1e_drinv_optimizer;
CINTIntegralFunction int1e_drinv_cart;
CINTIntegralFunction int1e_drinv_sph;
CINTIntegralFunction int1e_drinv_spinor;

/* (G i j|R12 |k l) */
CINTOptimizerFunction int2e_ig1_optimizer;
CINTIntegralFunction int2e_ig1_cart;
CINTIntegralFunction int2e_ig1_sph;
CINTIntegralFunction int2e_ig1_spinor;

/* (G G i j|R12 |k l) */
CINTOptimizerFunction int2e_gg1_optimizer;
CINTIntegralFunction int2e_gg1_cart;
CINTIntegralFunction int2e_gg1_sph;
CINTIntegralFunction int2e_gg1_spinor;

/* (G i j|R12 |G k l) */
CINTOptimizerFunction int2e_g1g2_optimizer;
CINTIntegralFunction int2e_g1g2_cart;
CINTIntegralFunction int2e_g1g2_sph;
CINTIntegralFunction int2e_g1g2_spinor;

/* (P* i CROSS P j|R12 |k l) */
CINTOptimizerFunction int2e_p1vxp1_optimizer;
CINTIntegralFunction int2e_p1vxp1_cart;
CINTIntegralFunction int2e_p1vxp1_sph;
CINTIntegralFunction int2e_p1vxp1_spinor;

/* (i RC j|NABLA-R12 |k l) */
CINTOptimizerFunction int2e_ip1v_rc1_optimizer;
CINTIntegralFunction int2e_ip1v_rc1_cart;
CINTIntegralFunction int2e_ip1v_rc1_sph;
CINTIntegralFunction int2e_ip1v_rc1_spinor;

/* (i R j|NABLA-R12 |k l) */
CINTOptimizerFunction int2e_ip1v_r1_optimizer;
CINTIntegralFunction int2e_ip1v_r1_cart;
CINTIntegralFunction int2e_ip1v_r1_sph;
CINTIntegralFunction int2e_ip1v_r1_spinor;

/* (G i j|NABLA-R12 CROSS P |k l) */
CINTOptimizerFunction int2e_ipvg1_xp1_optimizer;
CINTIntegralFunction int2e_ipvg1_xp1_cart;
CINTIntegralFunction int2e_ipvg1_xp1_sph;
CINTIntegralFunction int2e_ipvg1_xp1_spinor;

/* (i j|NABLA-R12 CROSS P |G k l) */
CINTOptimizerFunction int2e_ipvg2_xp1_optimizer;
CINTIntegralFunction int2e_ipvg2_xp1_cart;
CINTIntegralFunction int2e_ipvg2_xp1_sph;
CINTIntegralFunction int2e_ipvg2_xp1_spinor;

/* <i|NUC |RC CROSS P j> */
CINTOptimizerFunction int1e_inuc_rcxp_optimizer;
CINTIntegralFunction int1e_inuc_rcxp_cart;
CINTIntegralFunction int1e_inuc_rcxp_sph;
CINTIntegralFunction int1e_inuc_rcxp_spinor;

/* <i|NUC |R CROSS P j> */
CINTOptimizerFunction int1e_inuc_rxp_optimizer;
CINTIntegralFunction int1e_inuc_rxp_cart;
CINTIntegralFunction int1e_inuc_rxp_sph;
CINTIntegralFunction int1e_inuc_rxp_spinor;

/* <i|OVLP |SIGMA j> */
CINTOptimizerFunction int1e_sigma_optimizer;
CINTIntegralFunction int1e_sigma_cart;
CINTIntegralFunction int1e_sigma_sph;
CINTIntegralFunction int1e_sigma_spinor;

/* <SIGMA DOT P i|OVLP |SIGMA SIGMA DOT P j> */
CINTOptimizerFunction int1e_spsigmasp_optimizer;
CINTIntegralFunction int1e_spsigmasp_cart;
CINTIntegralFunction int1e_spsigmasp_sph;
CINTIntegralFunction int1e_spsigmasp_spinor;

/* <SIGMA DOT R i|OVLP |SIGMA DOT R j> */
CINTOptimizerFunction int1e_srsr_optimizer;
CINTIntegralFunction int1e_srsr_cart;
CINTIntegralFunction int1e_srsr_sph;
CINTIntegralFunction int1e_srsr_spinor;

/* <SIGMA DOT R i|OVLP |j> */
CINTOptimizerFunction int1e_sr_optimizer;
CINTIntegralFunction int1e_sr_cart;
CINTIntegralFunction int1e_sr_sph;
CINTIntegralFunction int1e_sr_spinor;

/* <SIGMA DOT R i|OVLP |SIGMA DOT P j> */
CINTOptimizerFunction int1e_srsp_optimizer;
CINTIntegralFunction int1e_srsp_cart;
CINTIntegralFunction int1e_srsp_sph;
CINTIntegralFunction int1e_srsp_spinor;

/* <SIGMA DOT P i|OVLP |SIGMA DOT P j> */
CINTOptimizerFunction int1e_spsp_optimizer;
CINTIntegralFunction int1e_spsp_cart;
CINTIntegralFunction int1e_spsp_sph;
CINTIntegralFunction int1e_spsp_spinor;

/* <SIGMA DOT P i|OVLP |j> */
CINTOptimizerFunction int1e_sp_optimizer;
CINTIntegralFunction int1e_sp_cart;
CINTIntegralFunction int1e_sp_sph;
CINTIntegralFunction int1e_sp_spinor;

/* <SIGMA DOT P i|NUC |SIGMA DOT P j> */
CINTOptimizerFunction int1e_spnucsp_optimizer;
CINTIntegralFunction int1e_spnucsp_cart;
CINTIntegralFunction int1e_spnucsp_sph;
CINTIntegralFunction int1e_spnucsp_spinor;

/* <SIGMA DOT P i|RINV |SIGMA DOT P j> */
CINTOptimizerFunction int1e_sprinvsp_optimizer;
CINTIntegralFunction int1e_sprinvsp_cart;
CINTIntegralFunction int1e_sprinvsp_sph;
CINTIntegralFunction int1e_sprinvsp_spinor;

/* <SIGMA DOT R i|NUC |SIGMA DOT R j> */
CINTOptimizerFunction int1e_srnucsr_optimizer;
CINTIntegralFunction int1e_srnucsr_cart;
CINTIntegralFunction int1e_srnucsr_sph;
CINTIntegralFunction int1e_srnucsr_spinor;

/* <SIGMA DOT P i|RC |SIGMA DOT P j> */
CINTOptimizerFunction int1e_sprsp_optimizer;
CINTIntegralFunction int1e_sprsp_cart;
CINTIntegralFunction int1e_sprsp_sph;
CINTIntegralFunction int1e_sprsp_spinor;

/* <G i|OVLP |j> */
CINTOptimizerFunction int1e_govlp_optimizer;
CINTIntegralFunction int1e_govlp_cart;
CINTIntegralFunction int1e_govlp_sph;
CINTIntegralFunction int1e_govlp_spinor;

/* <G i|NUC |j> */
CINTOptimizerFunction int1e_gnuc_optimizer;
CINTIntegralFunction int1e_gnuc_cart;
CINTIntegralFunction int1e_gnuc_sph;
CINTIntegralFunction int1e_gnuc_spinor;

/* <SIGMA CROSS RC i|SIGMA CROSS NABLA-RINV |j> */
CINTOptimizerFunction int1e_cg_sa10sa01_optimizer;
CINTIntegralFunction int1e_cg_sa10sa01_cart;
CINTIntegralFunction int1e_cg_sa10sa01_sph;
CINTIntegralFunction int1e_cg_sa10sa01_spinor;

/* <RC CROSS SIGMA i|OVLP |SIGMA DOT P j> */
CINTOptimizerFunction int1e_cg_sa10sp_optimizer;
CINTIntegralFunction int1e_cg_sa10sp_cart;
CINTIntegralFunction int1e_cg_sa10sp_sph;
CINTIntegralFunction int1e_cg_sa10sp_spinor;

/* <RC CROSS SIGMA i|NUC |SIGMA DOT P j> */
CINTOptimizerFunction int1e_cg_sa10nucsp_optimizer;
CINTIntegralFunction int1e_cg_sa10nucsp_cart;
CINTIntegralFunction int1e_cg_sa10nucsp_sph;
CINTIntegralFunction int1e_cg_sa10nucsp_spinor;

/* <SIGMA CROSS R i|SIGMA CROSS NABLA-RINV |j> */
CINTOptimizerFunction int1e_giao_sa10sa01_optimizer;
CINTIntegralFunction int1e_giao_sa10sa01_cart;
CINTIntegralFunction int1e_giao_sa10sa01_sph;
CINTIntegralFunction int1e_giao_sa10sa01_spinor;

/* <R CROSS SIGMA i|OVLP |SIGMA DOT P j> */
CINTOptimizerFunction int1e_giao_sa10sp_optimizer;
CINTIntegralFunction int1e_giao_sa10sp_cart;
CINTIntegralFunction int1e_giao_sa10sp_sph;
CINTIntegralFunction int1e_giao_sa10sp_spinor;

/* <R CROSS SIGMA i|NUC |SIGMA DOT P j> */
CINTOptimizerFunction int1e_giao_sa10nucsp_optimizer;
CINTIntegralFunction int1e_giao_sa10nucsp_cart;
CINTIntegralFunction int1e_giao_sa10nucsp_sph;
CINTIntegralFunction int1e_giao_sa10nucsp_spinor;

/* <i|NABLA-RINV CROSS SIGMA |SIGMA DOT P j> */
CINTOptimizerFunction int1e_sa01sp_optimizer;
CINTIntegralFunction int1e_sa01sp_cart;
CINTIntegralFunction int1e_sa01sp_sph;
CINTIntegralFunction int1e_sa01sp_spinor;

/* <G SIGMA DOT P i|OVLP |SIGMA DOT P j> */
CINTOptimizerFunction int1e_spgsp_optimizer;
CINTIntegralFunction int1e_spgsp_cart;
CINTIntegralFunction int1e_spgsp_sph;
CINTIntegralFunction int1e_spgsp_spinor;

/* <G SIGMA DOT P i|NUC |SIGMA DOT P j> */
CINTOptimizerFunction int1e_spgnucsp_optimizer;
CINTIntegralFunction int1e_spgnucsp_cart;
CINTIntegralFunction int1e_spgnucsp_sph;
CINTIntegralFunction int1e_spgnucsp_spinor;

/* <G SIGMA DOT P i|NABLA-RINV CROSS SIGMA |j> */
CINTOptimizerFunction int1e_spgsa01_optimizer;
CINTIntegralFunction int1e_spgsa01_cart;
CINTIntegralFunction int1e_spgsa01_sph;
CINTIntegralFunction int1e_spgsa01_spinor;

/* (SIGMA DOT P i SIGMA DOT P j|R12 |k l) */
CINTOptimizerFunction int2e_spsp1_optimizer;
CINTIntegralFunction int2e_spsp1_cart;
CINTIntegralFunction int2e_spsp1_sph;
CINTIntegralFunction int2e_spsp1_spinor;

/* (SIGMA DOT P i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CINTOptimizerFunction int2e_spsp1spsp2_optimizer;
CINTIntegralFunction int2e_spsp1spsp2_cart;
CINTIntegralFunction int2e_spsp1spsp2_sph;
CINTIntegralFunction int2e_spsp1spsp2_spinor;

/* (SIGMA DOT R i SIGMA DOT R j|R12 |k l) */
CINTOptimizerFunction int2e_srsr1_optimizer;
CINTIntegralFunction int2e_srsr1_cart;
CINTIntegralFunction int2e_srsr1_sph;
CINTIntegralFunction int2e_srsr1_spinor;

/* (SIGMA DOT R i SIGMA DOT R j|R12 |SIGMA DOT R k SIGMA DOT R l) */
CINTOptimizerFunction int2e_srsr1srsr2_optimizer;
CINTIntegralFunction int2e_srsr1srsr2_cart;
CINTIntegralFunction int2e_srsr1srsr2_sph;
CINTIntegralFunction int2e_srsr1srsr2_spinor;

/* (RC CROSS SIGMA i SIGMA DOT P j|R12 |k l) */
CINTOptimizerFunction int2e_cg_sa10sp1_optimizer;
CINTIntegralFunction int2e_cg_sa10sp1_cart;
CINTIntegralFunction int2e_cg_sa10sp1_sph;
CINTIntegralFunction int2e_cg_sa10sp1_spinor;

/* (RC CROSS SIGMA i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CINTOptimizerFunction int2e_cg_sa10sp1spsp2_optimizer;
CINTIntegralFunction int2e_cg_sa10sp1spsp2_cart;
CINTIntegralFunction int2e_cg_sa10sp1spsp2_sph;
CINTIntegralFunction int2e_cg_sa10sp1spsp2_spinor;

/* (R CROSS SIGMA i SIGMA DOT P j|R12 |k l) */
CINTOptimizerFunction int2e_giao_sa10sp1_optimizer;
CINTIntegralFunction int2e_giao_sa10sp1_cart;
CINTIntegralFunction int2e_giao_sa10sp1_sph;
CINTIntegralFunction int2e_giao_sa10sp1_spinor;

/* (R CROSS SIGMA i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CINTOptimizerFunction int2e_giao_sa10sp1spsp2_optimizer;
CINTIntegralFunction int2e_giao_sa10sp1spsp2_cart;
CINTIntegralFunction int2e_giao_sa10sp1spsp2_sph;
CINTIntegralFunction int2e_giao_sa10sp1spsp2_spinor;

/* (G i j|R12 |k l) */
CINTOptimizerFunction int2e_g1_optimizer;
CINTIntegralFunction int2e_g1_cart;
CINTIntegralFunction int2e_g1_sph;
CINTIntegralFunction int2e_g1_spinor;

/* (G SIGMA DOT P i SIGMA DOT P j|R12 |k l) */
CINTOptimizerFunction int2e_spgsp1_optimizer;
CINTIntegralFunction int2e_spgsp1_cart;
CINTIntegralFunction int2e_spgsp1_sph;
CINTIntegralFunction int2e_spgsp1_spinor;

/* (G i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CINTOptimizerFunction int2e_g1spsp2_optimizer;
CINTIntegralFunction int2e_g1spsp2_cart;
CINTIntegralFunction int2e_g1spsp2_sph;
CINTIntegralFunction int2e_g1spsp2_spinor;

/* (G SIGMA DOT P i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CINTOptimizerFunction int2e_spgsp1spsp2_optimizer;
CINTIntegralFunction int2e_spgsp1spsp2_cart;
CINTIntegralFunction int2e_spgsp1spsp2_sph;
CINTIntegralFunction int2e_spgsp1spsp2_spinor;

/* (P* i DOT P j|R12 |k l) */
CINTOptimizerFunction int2e_pp1_optimizer;
CINTIntegralFunction int2e_pp1_cart;
CINTIntegralFunction int2e_pp1_sph;
CINTIntegralFunction int2e_pp1_spinor;

/* (i j|R12 |P* k DOT P l) */
CINTOptimizerFunction int2e_pp2_optimizer;
CINTIntegralFunction int2e_pp2_cart;
CINTIntegralFunction int2e_pp2_sph;
CINTIntegralFunction int2e_pp2_spinor;

/* (P* i DOT P j|R12 |P* k DOT P l) */
CINTOptimizerFunction int2e_pp1pp2_optimizer;
CINTIntegralFunction int2e_pp1pp2_cart;
CINTIntegralFunction int2e_pp1pp2_sph;
CINTIntegralFunction int2e_pp1pp2_spinor;

/* <SIGMA DOT P i|OVLP |SIGMA DOT P SIGMA DOT P j> */
CINTOptimizerFunction int1e_spspsp_optimizer;
CINTIntegralFunction int1e_spspsp_cart;
CINTIntegralFunction int1e_spspsp_sph;
CINTIntegralFunction int1e_spspsp_spinor;

/* <SIGMA DOT P i|NUC |j> */
CINTOptimizerFunction int1e_spnuc_optimizer;
CINTIntegralFunction int1e_spnuc_cart;
CINTIntegralFunction int1e_spnuc_sph;
CINTIntegralFunction int1e_spnuc_spinor;

/* (SIGMA DOT P i j|R12 |k l) */
CINTOptimizerFunction int2e_spv1_optimizer;
CINTIntegralFunction int2e_spv1_cart;
CINTIntegralFunction int2e_spv1_sph;
CINTIntegralFunction int2e_spv1_spinor;

/* (i SIGMA DOT P j|R12 |k l) */
CINTOptimizerFunction int2e_vsp1_optimizer;
CINTIntegralFunction int2e_vsp1_cart;
CINTIntegralFunction int2e_vsp1_sph;
CINTIntegralFunction int2e_vsp1_spinor;

/* (i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CINTOptimizerFunction int2e_spsp2_optimizer;
CINTIntegralFunction int2e_spsp2_cart;
CINTIntegralFunction int2e_spsp2_sph;
CINTIntegralFunction int2e_spsp2_spinor;

/* (SIGMA DOT P i j|R12 |SIGMA DOT P k l) */
CINTOptimizerFunction int2e_spv1spv2_optimizer;
CINTIntegralFunction int2e_spv1spv2_cart;
CINTIntegralFunction int2e_spv1spv2_sph;
CINTIntegralFunction int2e_spv1spv2_spinor;

/* (i SIGMA DOT P j|R12 |SIGMA DOT P k l) */
CINTOptimizerFunction int2e_vsp1spv2_optimizer;
CINTIntegralFunction int2e_vsp1spv2_cart;
CINTIntegralFunction int2e_vsp1spv2_sph;
CINTIntegralFunction int2e_vsp1spv2_spinor;

/* (SIGMA DOT P i j|R12 |k SIGMA DOT P l) */
CINTOptimizerFunction int2e_spv1vsp2_optimizer;
CINTIntegralFunction int2e_spv1vsp2_cart;
CINTIntegralFunction int2e_spv1vsp2_sph;
CINTIntegralFunction int2e_spv1vsp2_spinor;

/* (i SIGMA DOT P j|R12 |k SIGMA DOT P l) */
CINTOptimizerFunction int2e_vsp1vsp2_optimizer;
CINTIntegralFunction int2e_vsp1vsp2_cart;
CINTIntegralFunction int2e_vsp1vsp2_sph;
CINTIntegralFunction int2e_vsp1vsp2_spinor;

/* (SIGMA DOT P i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CINTOptimizerFunction int2e_spv1spsp2_optimizer;
CINTIntegralFunction int2e_spv1spsp2_cart;
CINTIntegralFunction int2e_spv1spsp2_sph;
CINTIntegralFunction int2e_spv1spsp2_spinor;

/* (i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CINTOptimizerFunction int2e_vsp1spsp2_optimizer;
CINTIntegralFunction int2e_vsp1spsp2_cart;
CINTIntegralFunction int2e_vsp1spsp2_sph;
CINTIntegralFunction int2e_vsp1spsp2_spinor;

/* <NABLA i|OVLP |j> */
CINTOptimizerFunction int1e_ipovlp_optimizer;
CINTIntegralFunction int1e_ipovlp_cart;
CINTIntegralFunction int1e_ipovlp_sph;
CINTIntegralFunction int1e_ipovlp_spinor;

/* <i|OVLP |NABLA j> */
CINTOptimizerFunction int1e_ovlpip_optimizer;
CINTIntegralFunction int1e_ovlpip_cart;
CINTIntegralFunction int1e_ovlpip_sph;
CINTIntegralFunction int1e_ovlpip_spinor;

/* <NABLA i|OVLP |P DOT P j> */
CINTOptimizerFunction int1e_ipkin_optimizer;
CINTIntegralFunction int1e_ipkin_cart;
CINTIntegralFunction int1e_ipkin_sph;
CINTIntegralFunction int1e_ipkin_spinor;

/* <i|OVLP |P DOT P NABLA j> */
CINTOptimizerFunction int1e_kinip_optimizer;
CINTIntegralFunction int1e_kinip_cart;
CINTIntegralFunction int1e_kinip_sph;
CINTIntegralFunction int1e_kinip_spinor;

/* <NABLA i|NUC |j> */
CINTOptimizerFunction int1e_ipnuc_optimizer;
CINTIntegralFunction int1e_ipnuc_cart;
CINTIntegralFunction int1e_ipnuc_sph;
CINTIntegralFunction int1e_ipnuc_spinor;

/* <NABLA i|RINV |j> */
CINTOptimizerFunction int1e_iprinv_optimizer;
CINTIntegralFunction int1e_iprinv_cart;
CINTIntegralFunction int1e_iprinv_sph;
CINTIntegralFunction int1e_iprinv_spinor;

/* <NABLA SIGMA DOT P i|NUC |SIGMA DOT P j> */
CINTOptimizerFunction int1e_ipspnucsp_optimizer;
CINTIntegralFunction int1e_ipspnucsp_cart;
CINTIntegralFunction int1e_ipspnucsp_sph;
CINTIntegralFunction int1e_ipspnucsp_spinor;

/* <NABLA SIGMA DOT P i|RINV |SIGMA DOT P j> */
CINTOptimizerFunction int1e_ipsprinvsp_optimizer;
CINTIntegralFunction int1e_ipsprinvsp_cart;
CINTIntegralFunction int1e_ipsprinvsp_sph;
CINTIntegralFunction int1e_ipsprinvsp_spinor;

/* <P* NABLA i|NUC DOT P |j> */
CINTOptimizerFunction int1e_ippnucp_optimizer;
CINTIntegralFunction int1e_ippnucp_cart;
CINTIntegralFunction int1e_ippnucp_sph;
CINTIntegralFunction int1e_ippnucp_spinor;

/* <P* NABLA i|RINV DOT P |j> */
CINTOptimizerFunction int1e_ipprinvp_optimizer;
CINTIntegralFunction int1e_ipprinvp_cart;
CINTIntegralFunction int1e_ipprinvp_sph;
CINTIntegralFunction int1e_ipprinvp_spinor;

/* (NABLA i j|R12 |k l) */
CINTOptimizerFunction int2e_ip1_optimizer;
CINTIntegralFunction int2e_ip1_cart;
CINTIntegralFunction int2e_ip1_sph;
CINTIntegralFunction int2e_ip1_spinor;

/* (i j|R12 |NABLA k l) */
CINTOptimizerFunction int2e_ip2_optimizer;
CINTIntegralFunction int2e_ip2_cart;
CINTIntegralFunction int2e_ip2_sph;
CINTIntegralFunction int2e_ip2_spinor;

/* (NABLA SIGMA DOT P i SIGMA DOT P j|R12 |k l) */
CINTOptimizerFunction int2e_ipspsp1_optimizer;
CINTIntegralFunction int2e_ipspsp1_cart;
CINTIntegralFunction int2e_ipspsp1_sph;
CINTIntegralFunction int2e_ipspsp1_spinor;

/* (NABLA i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CINTOptimizerFunction int2e_ip1spsp2_optimizer;
CINTIntegralFunction int2e_ip1spsp2_cart;
CINTIntegralFunction int2e_ip1spsp2_sph;
CINTIntegralFunction int2e_ip1spsp2_spinor;

/* (NABLA SIGMA DOT P i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CINTOptimizerFunction int2e_ipspsp1spsp2_optimizer;
CINTIntegralFunction int2e_ipspsp1spsp2_cart;
CINTIntegralFunction int2e_ipspsp1spsp2_sph;
CINTIntegralFunction int2e_ipspsp1spsp2_spinor;

/* (NABLA SIGMA DOT R i SIGMA DOT R j|R12 |k l) */
CINTOptimizerFunction int2e_ipsrsr1_optimizer;
CINTIntegralFunction int2e_ipsrsr1_cart;
CINTIntegralFunction int2e_ipsrsr1_sph;
CINTIntegralFunction int2e_ipsrsr1_spinor;

/* (NABLA i j|R12 |SIGMA DOT R k SIGMA DOT R l) */
CINTOptimizerFunction int2e_ip1srsr2_optimizer;
CINTIntegralFunction int2e_ip1srsr2_cart;
CINTIntegralFunction int2e_ip1srsr2_sph;
CINTIntegralFunction int2e_ip1srsr2_spinor;

/* (NABLA SIGMA DOT R i SIGMA DOT R j|R12 |SIGMA DOT R k SIGMA DOT R l) */
CINTOptimizerFunction int2e_ipsrsr1srsr2_optimizer;
CINTIntegralFunction int2e_ipsrsr1srsr2_cart;
CINTIntegralFunction int2e_ipsrsr1srsr2_sph;
CINTIntegralFunction int2e_ipsrsr1srsr2_spinor;

/* (i SIGMA DOT P j|GAUNT |k SIGMA DOT P l) */
CINTOptimizerFunction int2e_ssp1ssp2_optimizer;
CINTIntegralFunction int2e_ssp1ssp2_cart;
CINTIntegralFunction int2e_ssp1ssp2_sph;
CINTIntegralFunction int2e_ssp1ssp2_spinor;

/* (i SIGMA DOT P j|GAUNT |SIGMA DOT P k l) */
CINTOptimizerFunction int2e_ssp1sps2_optimizer;
CINTIntegralFunction int2e_ssp1sps2_cart;
CINTIntegralFunction int2e_ssp1sps2_sph;
CINTIntegralFunction int2e_ssp1sps2_spinor;

/* (SIGMA DOT P i j|GAUNT |k SIGMA DOT P l) */
CINTOptimizerFunction int2e_sps1ssp2_optimizer;
CINTIntegralFunction int2e_sps1ssp2_cart;
CINTIntegralFunction int2e_sps1ssp2_sph;
CINTIntegralFunction int2e_sps1ssp2_spinor;

/* (SIGMA DOT P i j|GAUNT |SIGMA DOT P k l) */
CINTOptimizerFunction int2e_sps1sps2_optimizer;
CINTIntegralFunction int2e_sps1sps2_cart;
CINTIntegralFunction int2e_sps1sps2_sph;
CINTIntegralFunction int2e_sps1sps2_spinor;

/* (RC CROSS SIGMA i j|GAUNT |k SIGMA DOT P l) */
CINTOptimizerFunction int2e_cg_ssa10ssp2_optimizer;
CINTIntegralFunction int2e_cg_ssa10ssp2_cart;
CINTIntegralFunction int2e_cg_ssa10ssp2_sph;
CINTIntegralFunction int2e_cg_ssa10ssp2_spinor;

/* (R CROSS SIGMA i j|GAUNT |k SIGMA DOT P l) */
CINTOptimizerFunction int2e_giao_ssa10ssp2_optimizer;
CINTIntegralFunction int2e_giao_ssa10ssp2_cart;
CINTIntegralFunction int2e_giao_ssa10ssp2_sph;
CINTIntegralFunction int2e_giao_ssa10ssp2_spinor;

/* (G i SIGMA DOT P j|GAUNT |k SIGMA DOT P l) */
CINTOptimizerFunction int2e_gssp1ssp2_optimizer;
CINTIntegralFunction int2e_gssp1ssp2_cart;
CINTIntegralFunction int2e_gssp1ssp2_sph;
CINTIntegralFunction int2e_gssp1ssp2_spinor;

/* (i R0 SIGMA DOT P j|BREIT-R1 |k SIGMA DOT P l) */
CINTOptimizerFunction int2e_gauge_r1_ssp1ssp2_optimizer;
CINTIntegralFunction int2e_gauge_r1_ssp1ssp2_cart;
CINTIntegralFunction int2e_gauge_r1_ssp1ssp2_sph;
CINTIntegralFunction int2e_gauge_r1_ssp1ssp2_spinor;

/* (i R0 SIGMA DOT P j|BREIT-R1 |SIGMA DOT P k l) */
CINTOptimizerFunction int2e_gauge_r1_ssp1sps2_optimizer;
CINTIntegralFunction int2e_gauge_r1_ssp1sps2_cart;
CINTIntegralFunction int2e_gauge_r1_ssp1sps2_sph;
CINTIntegralFunction int2e_gauge_r1_ssp1sps2_spinor;

/* (SIGMA DOT P i R0 j|BREIT-R1 |k SIGMA DOT P l) */
CINTOptimizerFunction int2e_gauge_r1_sps1ssp2_optimizer;
CINTIntegralFunction int2e_gauge_r1_sps1ssp2_cart;
CINTIntegralFunction int2e_gauge_r1_sps1ssp2_sph;
CINTIntegralFunction int2e_gauge_r1_sps1ssp2_spinor;

/* (SIGMA DOT P i R0 j|BREIT-R1 |SIGMA DOT P k l) */
CINTOptimizerFunction int2e_gauge_r1_sps1sps2_optimizer;
CINTIntegralFunction int2e_gauge_r1_sps1sps2_cart;
CINTIntegralFunction int2e_gauge_r1_sps1sps2_sph;
CINTIntegralFunction int2e_gauge_r1_sps1sps2_spinor;

/* (i SIGMA DOT P j|BREIT-R2 |k R0 SIGMA DOT P l) */
CINTOptimizerFunction int2e_gauge_r2_ssp1ssp2_optimizer;
CINTIntegralFunction int2e_gauge_r2_ssp1ssp2_cart;
CINTIntegralFunction int2e_gauge_r2_ssp1ssp2_sph;
CINTIntegralFunction int2e_gauge_r2_ssp1ssp2_spinor;

/* (i SIGMA DOT P j|BREIT-R2 |SIGMA DOT P k R0 l) */
CINTOptimizerFunction int2e_gauge_r2_ssp1sps2_optimizer;
CINTIntegralFunction int2e_gauge_r2_ssp1sps2_cart;
CINTIntegralFunction int2e_gauge_r2_ssp1sps2_sph;
CINTIntegralFunction int2e_gauge_r2_ssp1sps2_spinor;

/* (SIGMA DOT P i j|BREIT-R2 |k R0 SIGMA DOT P l) */
CINTOptimizerFunction int2e_gauge_r2_sps1ssp2_optimizer;
CINTIntegralFunction int2e_gauge_r2_sps1ssp2_cart;
CINTIntegralFunction int2e_gauge_r2_sps1ssp2_sph;
CINTIntegralFunction int2e_gauge_r2_sps1ssp2_spinor;

/* (SIGMA DOT P i j|BREIT-R2 |SIGMA DOT P k R0 l) */
CINTOptimizerFunction int2e_gauge_r2_sps1sps2_optimizer;
CINTIntegralFunction int2e_gauge_r2_sps1sps2_cart;
CINTIntegralFunction int2e_gauge_r2_sps1sps2_sph;
CINTIntegralFunction int2e_gauge_r2_sps1sps2_spinor;

/* <NABLA NABLA i|OVLP |j> */
CINTOptimizerFunction int1e_ipipovlp_optimizer;
CINTIntegralFunction int1e_ipipovlp_cart;
CINTIntegralFunction int1e_ipipovlp_sph;
CINTIntegralFunction int1e_ipipovlp_spinor;

/* <NABLA i|OVLP |NABLA j> */
CINTOptimizerFunction int1e_ipovlpip_optimizer;
CINTIntegralFunction int1e_ipovlpip_cart;
CINTIntegralFunction int1e_ipovlpip_sph;
CINTIntegralFunction int1e_ipovlpip_spinor;

/* <NABLA NABLA i|P DOT P |j> */
CINTOptimizerFunction int1e_ipipkin_optimizer;
CINTIntegralFunction int1e_ipipkin_cart;
CINTIntegralFunction int1e_ipipkin_sph;
CINTIntegralFunction int1e_ipipkin_spinor;

/* <NABLA i|P DOT P |NABLA j> */
CINTOptimizerFunction int1e_ipkinip_optimizer;
CINTIntegralFunction int1e_ipkinip_cart;
CINTIntegralFunction int1e_ipkinip_sph;
CINTIntegralFunction int1e_ipkinip_spinor;

/* <NABLA NABLA i|NUC |j> */
CINTOptimizerFunction int1e_ipipnuc_optimizer;
CINTIntegralFunction int1e_ipipnuc_cart;
CINTIntegralFunction int1e_ipipnuc_sph;
CINTIntegralFunction int1e_ipipnuc_spinor;

/* <NABLA i|NUC |NABLA j> */
CINTOptimizerFunction int1e_ipnucip_optimizer;
CINTIntegralFunction int1e_ipnucip_cart;
CINTIntegralFunction int1e_ipnucip_sph;
CINTIntegralFunction int1e_ipnucip_spinor;

/* <NABLA NABLA i|RINV |j> */
CINTOptimizerFunction int1e_ipiprinv_optimizer;
CINTIntegralFunction int1e_ipiprinv_cart;
CINTIntegralFunction int1e_ipiprinv_sph;
CINTIntegralFunction int1e_ipiprinv_spinor;

/* <NABLA i|RINV |NABLA j> */
CINTOptimizerFunction int1e_iprinvip_optimizer;
CINTIntegralFunction int1e_iprinvip_cart;
CINTIntegralFunction int1e_iprinvip_sph;
CINTIntegralFunction int1e_iprinvip_spinor;

/* (NABLA NABLA i j|R12 |k l) */
CINTOptimizerFunction int2e_ipip1_optimizer;
CINTIntegralFunction int2e_ipip1_cart;
CINTIntegralFunction int2e_ipip1_sph;
CINTIntegralFunction int2e_ipip1_spinor;

/* (NABLA i NABLA j|R12 |k l) */
CINTOptimizerFunction int2e_ipvip1_optimizer;
CINTIntegralFunction int2e_ipvip1_cart;
CINTIntegralFunction int2e_ipvip1_sph;
CINTIntegralFunction int2e_ipvip1_spinor;

/* (NABLA i j|R12 |NABLA k l) */
CINTOptimizerFunction int2e_ip1ip2_optimizer;
CINTIntegralFunction int2e_ip1ip2_cart;
CINTIntegralFunction int2e_ip1ip2_sph;
CINTIntegralFunction int2e_ip1ip2_spinor;

/* <P* NABLA NABLA i|NUC DOT P |j> */
CINTOptimizerFunction int1e_ipippnucp_optimizer;
CINTIntegralFunction int1e_ipippnucp_cart;
CINTIntegralFunction int1e_ipippnucp_sph;
CINTIntegralFunction int1e_ipippnucp_spinor;

/* <P* NABLA i|NUC DOT P |NABLA j> */
CINTOptimizerFunction int1e_ippnucpip_optimizer;
CINTIntegralFunction int1e_ippnucpip_cart;
CINTIntegralFunction int1e_ippnucpip_sph;
CINTIntegralFunction int1e_ippnucpip_spinor;

/* <P* NABLA NABLA i|RINV DOT P |j> */
CINTOptimizerFunction int1e_ipipprinvp_optimizer;
CINTIntegralFunction int1e_ipipprinvp_cart;
CINTIntegralFunction int1e_ipipprinvp_sph;
CINTIntegralFunction int1e_ipipprinvp_spinor;

/* <P* NABLA i|RINV DOT P |NABLA j> */
CINTOptimizerFunction int1e_ipprinvpip_optimizer;
CINTIntegralFunction int1e_ipprinvpip_cart;
CINTIntegralFunction int1e_ipprinvpip_sph;
CINTIntegralFunction int1e_ipprinvpip_spinor;

/* <NABLA NABLA SIGMA DOT P i|NUC SIGMA DOT P |j> */
CINTOptimizerFunction int1e_ipipspnucsp_optimizer;
CINTIntegralFunction int1e_ipipspnucsp_cart;
CINTIntegralFunction int1e_ipipspnucsp_sph;
CINTIntegralFunction int1e_ipipspnucsp_spinor;

/* <NABLA SIGMA DOT P i|NUC SIGMA DOT P |NABLA j> */
CINTOptimizerFunction int1e_ipspnucspip_optimizer;
CINTIntegralFunction int1e_ipspnucspip_cart;
CINTIntegralFunction int1e_ipspnucspip_sph;
CINTIntegralFunction int1e_ipspnucspip_spinor;

/* <NABLA NABLA SIGMA DOT P i|RINV SIGMA DOT P |j> */
CINTOptimizerFunction int1e_ipipsprinvsp_optimizer;
CINTIntegralFunction int1e_ipipsprinvsp_cart;
CINTIntegralFunction int1e_ipipsprinvsp_sph;
CINTIntegralFunction int1e_ipipsprinvsp_spinor;

/* <NABLA SIGMA DOT P i|RINV SIGMA DOT P |NABLA j> */
CINTOptimizerFunction int1e_ipsprinvspip_optimizer;
CINTIntegralFunction int1e_ipsprinvspip_cart;
CINTIntegralFunction int1e_ipsprinvspip_sph;
CINTIntegralFunction int1e_ipsprinvspip_spinor;

/* (NABLA NABLA i j|R12 |NABLA NABLA k l) */
CINTOptimizerFunction int2e_ipip1ipip2_optimizer;
CINTIntegralFunction int2e_ipip1ipip2_cart;
CINTIntegralFunction int2e_ipip1ipip2_sph;
CINTIntegralFunction int2e_ipip1ipip2_spinor;

/* (NABLA i NABLA j|R12 |NABLA k NABLA l) */
CINTOptimizerFunction int2e_ipvip1ipvip2_optimizer;
CINTIntegralFunction int2e_ipvip1ipvip2_cart;
CINTIntegralFunction int2e_ipvip1ipvip2_sph;
CINTIntegralFunction int2e_ipvip1ipvip2_spinor;

/* (NABLA i j|R12 |k) */
CINTOptimizerFunction int3c2e_ip1_optimizer;
CINTIntegralFunction int3c2e_ip1_cart;
CINTIntegralFunction int3c2e_ip1_sph;
CINTIntegralFunction int3c2e_ip1_spinor;

/* (i j|R12 |NABLA k) */
CINTOptimizerFunction int3c2e_ip2_optimizer;
CINTIntegralFunction int3c2e_ip2_cart;
CINTIntegralFunction int3c2e_ip2_sph;
CINTIntegralFunction int3c2e_ip2_spinor;

/* (P* i DOT P j|R12 |k) */
CINTOptimizerFunction int3c2e_pvp1_optimizer;
CINTIntegralFunction int3c2e_pvp1_cart;
CINTIntegralFunction int3c2e_pvp1_sph;
CINTIntegralFunction int3c2e_pvp1_spinor;

/* (P* i CROSS P j|R12 |k) */
CINTOptimizerFunction int3c2e_pvxp1_optimizer;
CINTIntegralFunction int3c2e_pvxp1_cart;
CINTIntegralFunction int3c2e_pvxp1_sph;
CINTIntegralFunction int3c2e_pvxp1_spinor;

/* (NABLA i |R12 |j) */
CINTOptimizerFunction int2c2e_ip1_optimizer;
CINTIntegralFunction int2c2e_ip1_cart;
CINTIntegralFunction int2c2e_ip1_sph;
CINTIntegralFunction int2c2e_ip1_spinor;

/* (i |R12 |NABLA j) */
CINTOptimizerFunction int2c2e_ip2_optimizer;
CINTIntegralFunction int2c2e_ip2_cart;
CINTIntegralFunction int2c2e_ip2_sph;
CINTIntegralFunction int2c2e_ip2_spinor;

/* (G i j|R12 |k) */
CINTOptimizerFunction int3c2e_ig1_optimizer;
CINTIntegralFunction int3c2e_ig1_cart;
CINTIntegralFunction int3c2e_ig1_sph;
CINTIntegralFunction int3c2e_ig1_spinor;

/* (SIGMA DOT P i SIGMA DOT P j|R12 |k) */
CINTOptimizerFunction int3c2e_spsp1_optimizer;
CINTIntegralFunction int3c2e_spsp1_cart;
CINTIntegralFunction int3c2e_spsp1_sph;
CINTIntegralFunction int3c2e_spsp1_spinor;

/* (NABLA SIGMA DOT P i SIGMA DOT P j|R12 |k) */
CINTOptimizerFunction int3c2e_ipspsp1_optimizer;
CINTIntegralFunction int3c2e_ipspsp1_cart;
CINTIntegralFunction int3c2e_ipspsp1_sph;
CINTIntegralFunction int3c2e_ipspsp1_spinor;

/* (SIGMA DOT P i SIGMA DOT P j|R12 |NABLA k) */
CINTOptimizerFunction int3c2e_spsp1ip2_optimizer;
CINTIntegralFunction int3c2e_spsp1ip2_cart;
CINTIntegralFunction int3c2e_spsp1ip2_sph;
CINTIntegralFunction int3c2e_spsp1ip2_spinor;

/* (NABLA NABLA i j|R12 |k) */
CINTOptimizerFunction int3c2e_ipip1_optimizer;
CINTIntegralFunction int3c2e_ipip1_cart;
CINTIntegralFunction int3c2e_ipip1_sph;
CINTIntegralFunction int3c2e_ipip1_spinor;

/* (i j|R12 |NABLA NABLA k) */
CINTOptimizerFunction int3c2e_ipip2_optimizer;
CINTIntegralFunction int3c2e_ipip2_cart;
CINTIntegralFunction int3c2e_ipip2_sph;
CINTIntegralFunction int3c2e_ipip2_spinor;

/* (NABLA i NABLA j|R12 |k) */
CINTOptimizerFunction int3c2e_ipvip1_optimizer;
CINTIntegralFunction int3c2e_ipvip1_cart;
CINTIntegralFunction int3c2e_ipvip1_sph;
CINTIntegralFunction int3c2e_ipvip1_spinor;

/* (NABLA i j|R12 |NABLA k) */
CINTOptimizerFunction int3c2e_ip1ip2_optimizer;
CINTIntegralFunction int3c2e_ip1ip2_cart;
CINTIntegralFunction int3c2e_ip1ip2_sph;
CINTIntegralFunction int3c2e_ip1ip2_spinor;

/* (NABLA NABLA i |R12 |j) */
CINTOptimizerFunction int2c2e_ipip1_optimizer;
CINTIntegralFunction int2c2e_ipip1_cart;
CINTIntegralFunction int2c2e_ipip1_sph;
CINTIntegralFunction int2c2e_ipip1_spinor;

/* (NABLA i |R12 |NABLA j) */
CINTOptimizerFunction int2c2e_ip1ip2_optimizer;
CINTIntegralFunction int2c2e_ip1ip2_cart;
CINTIntegralFunction int2c2e_ip1ip2_sph;
CINTIntegralFunction int2c2e_ip1ip2_spinor;

/* 3-center 1-electron integral <(i) (j) (P DOT P k)> */
CINTOptimizerFunction int3c1e_p2_optimizer;
CINTIntegralFunction int3c1e_p2_cart;
CINTIntegralFunction int3c1e_p2_sph;
CINTIntegralFunction int3c1e_p2_spinor;

/* 3-center 1-electron integral <(P i) (j) (k)> */
CINTOptimizerFunction int3c1e_iprinv_optimizer;
CINTIntegralFunction int3c1e_iprinv_cart;
CINTIntegralFunction int3c1e_iprinv_sph;
CINTIntegralFunction int3c1e_iprinv_spinor;

