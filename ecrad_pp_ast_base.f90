MODULE radiation_adding_ica_lw
  CONTAINS
  SUBROUTINE calc_fluxes_no_scattering_lw(ncol_var_0, nlev_var_1, transmittance_var_2, source_up_var_3, source_dn_var_4, emission_surf, albedo_surf, flux_up_var_5, flux_dn_var_6)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ncol_var_0
    INTEGER, INTENT(IN) :: nlev_var_1
    REAL(KIND = 8), INTENT(IN), DIMENSION(ncol_var_0) :: emission_surf, albedo_surf
    REAL(KIND = 8), INTENT(IN), DIMENSION(ncol_var_0, nlev_var_1) :: transmittance_var_2
    REAL(KIND = 8), INTENT(IN), DIMENSION(ncol_var_0, nlev_var_1) :: source_up_var_3, source_dn_var_4
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ncol_var_0, nlev_var_1 + 1) :: flux_up_var_5, flux_dn_var_6
    INTEGER :: jlev_var_7, jcol_var_8
    flux_dn_var_6(:, 1) = 0.0D0
    DO jlev_var_7 = 1, nlev_var_1
      DO jcol_var_8 = 1, ncol_var_0
        flux_dn_var_6(jcol_var_8, jlev_var_7 + 1) = transmittance_var_2(jcol_var_8, jlev_var_7) * flux_dn_var_6(jcol_var_8, jlev_var_7) + source_dn_var_4(jcol_var_8, jlev_var_7)
      END DO
    END DO
    flux_up_var_5(:, nlev_var_1 + 1) = emission_surf + albedo_surf * flux_dn_var_6(:, nlev_var_1 + 1)
    DO jlev_var_7 = nlev_var_1, 1, - 1
      DO jcol_var_8 = 1, ncol_var_0
        flux_up_var_5(jcol_var_8, jlev_var_7) = transmittance_var_2(jcol_var_8, jlev_var_7) * flux_up_var_5(jcol_var_8, jlev_var_7 + 1) + source_up_var_3(jcol_var_8, jlev_var_7)
      END DO
    END DO
  END SUBROUTINE calc_fluxes_no_scattering_lw
END MODULE radiation_adding_ica_lw
MODULE radiation_adding_ica_sw
  CONTAINS
  SUBROUTINE adding_ica_sw(ncol_var_9, nlev_var_10, incoming_toa, albedo_surf_diffuse, albedo_surf_direct, cos_sza_var_11, reflectance_var_12, transmittance_var_13, ref_dir_var_14, trans_dir_diff_var_15, trans_dir_dir_var_16, flux_up_var_17, flux_dn_diffuse_var_18, flux_dn_direct_var_19, albedo_var_20, source_var_21, inv_denominator)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ncol_var_9
    INTEGER, INTENT(IN) :: nlev_var_10
    REAL(KIND = 8), INTENT(IN), DIMENSION(ncol_var_9) :: incoming_toa
    REAL(KIND = 8), INTENT(IN), DIMENSION(ncol_var_9) :: albedo_surf_diffuse, albedo_surf_direct
    REAL(KIND = 8), INTENT(IN) :: cos_sza_var_11
    REAL(KIND = 8), INTENT(IN), DIMENSION(ncol_var_9, nlev_var_10) :: reflectance_var_12, transmittance_var_13
    REAL(KIND = 8), INTENT(IN), DIMENSION(ncol_var_9, nlev_var_10) :: ref_dir_var_14, trans_dir_diff_var_15
    REAL(KIND = 8), INTENT(IN), DIMENSION(ncol_var_9, nlev_var_10) :: trans_dir_dir_var_16
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ncol_var_9, nlev_var_10 + 1) :: flux_up_var_17, flux_dn_diffuse_var_18, flux_dn_direct_var_19
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ncol_var_9, nlev_var_10 + 1) :: albedo_var_20
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ncol_var_9, nlev_var_10 + 1) :: source_var_21
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ncol_var_9, nlev_var_10) :: inv_denominator
    INTEGER :: jlev_var_22, jcol_var_23
    flux_dn_direct_var_19(:, 1) = incoming_toa
    DO jlev_var_22 = 1, nlev_var_10
      flux_dn_direct_var_19(:, jlev_var_22 + 1) = flux_dn_direct_var_19(:, jlev_var_22) * trans_dir_dir_var_16(:, jlev_var_22)
    END DO
    albedo_var_20(:, nlev_var_10 + 1) = albedo_surf_diffuse
    source_var_21(:, nlev_var_10 + 1) = albedo_surf_direct * flux_dn_direct_var_19(:, nlev_var_10 + 1) * cos_sza_var_11
    DO jlev_var_22 = nlev_var_10, 1, - 1
      DO jcol_var_23 = 1, ncol_var_9
        inv_denominator(jcol_var_23, jlev_var_22) = 1.0D0 / (1.0D0 - albedo_var_20(jcol_var_23, jlev_var_22 + 1) * reflectance_var_12(jcol_var_23, jlev_var_22))
        albedo_var_20(jcol_var_23, jlev_var_22) = reflectance_var_12(jcol_var_23, jlev_var_22) + transmittance_var_13(jcol_var_23, jlev_var_22) * transmittance_var_13(jcol_var_23, jlev_var_22) * albedo_var_20(jcol_var_23, jlev_var_22 + 1) * inv_denominator(jcol_var_23, jlev_var_22)
        source_var_21(jcol_var_23, jlev_var_22) = ref_dir_var_14(jcol_var_23, jlev_var_22) * flux_dn_direct_var_19(jcol_var_23, jlev_var_22) + transmittance_var_13(jcol_var_23, jlev_var_22) * (source_var_21(jcol_var_23, jlev_var_22 + 1) + albedo_var_20(jcol_var_23, jlev_var_22 + 1) * trans_dir_diff_var_15(jcol_var_23, jlev_var_22) * flux_dn_direct_var_19(jcol_var_23, jlev_var_22)) * inv_denominator(jcol_var_23, jlev_var_22)
      END DO
    END DO
    flux_dn_diffuse_var_18(:, 1) = 0.0D0
    flux_up_var_17(:, 1) = source_var_21(:, 1)
    DO jlev_var_22 = 1, nlev_var_10
      DO jcol_var_23 = 1, ncol_var_9
        flux_dn_diffuse_var_18(jcol_var_23, jlev_var_22 + 1) = (transmittance_var_13(jcol_var_23, jlev_var_22) * flux_dn_diffuse_var_18(jcol_var_23, jlev_var_22) + reflectance_var_12(jcol_var_23, jlev_var_22) * source_var_21(jcol_var_23, jlev_var_22 + 1) + trans_dir_diff_var_15(jcol_var_23, jlev_var_22) * flux_dn_direct_var_19(jcol_var_23, jlev_var_22)) * inv_denominator(jcol_var_23, jlev_var_22)
        flux_up_var_17(jcol_var_23, jlev_var_22 + 1) = albedo_var_20(jcol_var_23, jlev_var_22 + 1) * flux_dn_diffuse_var_18(jcol_var_23, jlev_var_22 + 1) + source_var_21(jcol_var_23, jlev_var_22 + 1)
        flux_dn_direct_var_19(jcol_var_23, jlev_var_22) = flux_dn_direct_var_19(jcol_var_23, jlev_var_22) * cos_sza_var_11
      END DO
    END DO
    flux_dn_direct_var_19(:, nlev_var_10 + 1) = flux_dn_direct_var_19(:, nlev_var_10 + 1) * cos_sza_var_11
  END SUBROUTINE adding_ica_sw
END MODULE radiation_adding_ica_sw
MODULE radiation_cloud_cover
  CONTAINS
  ELEMENTAL FUNCTION beta2alpha_fn_24(beta_var_25, frac1, frac2)
    IMPLICIT NONE
    REAL(KIND = 8), INTENT(IN) :: beta_var_25, frac1, frac2
    REAL(KIND = 8) :: beta2alpha_fn_24
    REAL(KIND = 8) :: frac_diff
    IF (beta_var_25 < 1.0D0) THEN
      frac_diff = ABS(frac1 - frac2)
      beta2alpha_fn_24 = beta_var_25 + (1.0D0 - beta_var_25) * frac_diff / (frac_diff + 1.0D0 / beta_var_25 - 1.0D0)
    ELSE
      beta2alpha_fn_24 = 1.0D0
    END IF
  END FUNCTION beta2alpha_fn_24
  SUBROUTINE cum_cloud_cover_exp_ran(nlev_var_27, frac_var_28, overlap_param_var_29, cum_cloud_cover_var_30, pair_cloud_cover_var_31, is_beta_overlap)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlev_var_27
    REAL(KIND = 8), INTENT(IN) :: frac_var_28(nlev_var_27)
    REAL(KIND = 8), INTENT(IN) :: overlap_param_var_29(nlev_var_27 - 1)
    LOGICAL, INTENT(IN), OPTIONAL :: is_beta_overlap
    REAL(KIND = 8), INTENT(OUT) :: cum_cloud_cover_var_30(nlev_var_27)
    REAL(KIND = 8), INTENT(OUT) :: pair_cloud_cover_var_31(nlev_var_27 - 1)
    REAL(KIND = 8) :: cum_product
    REAL(KIND = 8) :: overlap_alpha
    LOGICAL :: do_overlap_conversion
    INTEGER :: jlev_var_32
    do_overlap_conversion = is_beta_overlap
    cum_product = 1.0D0 - frac_var_28(1)
    cum_cloud_cover_var_30(1) = frac_var_28(1)
    DO jlev_var_32 = 1, nlev_var_27 - 1
      IF (do_overlap_conversion) THEN
        overlap_alpha = beta2alpha_fn_24(overlap_param_var_29(jlev_var_32), frac_var_28(jlev_var_32), frac_var_28(jlev_var_32 + 1))
      ELSE
        overlap_alpha = overlap_param_var_29(jlev_var_32)
      END IF
      pair_cloud_cover_var_31(jlev_var_32) = overlap_alpha * MAX(frac_var_28(jlev_var_32), frac_var_28(jlev_var_32 + 1)) + (1.0D0 - overlap_alpha) * (frac_var_28(jlev_var_32) + frac_var_28(jlev_var_32 + 1) - frac_var_28(jlev_var_32) * frac_var_28(jlev_var_32 + 1))
      IF (frac_var_28(jlev_var_32) >= 0.9999999999999978D0) THEN
        cum_product = 0.0D0
      ELSE
        cum_product = cum_product * (1.0D0 - pair_cloud_cover_var_31(jlev_var_32)) / (1.0D0 - frac_var_28(jlev_var_32))
      END IF
      cum_cloud_cover_var_30(jlev_var_32 + 1) = 1.0D0 - cum_product
    END DO
  END SUBROUTINE cum_cloud_cover_exp_ran
END MODULE radiation_cloud_cover
MODULE radiation_ice_optics_fu
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE calc_ice_optics_fu_sw(nb_var_33, coeff_var_34, ice_wp_var_35, re_var_36, od_var_37, scat_od_var_38, g_var_39)
    INTEGER, INTENT(IN) :: nb_var_33
    REAL(KIND = 8), INTENT(IN) :: coeff_var_34(:, :)
    REAL(KIND = 8), INTENT(IN) :: ice_wp_var_35
    REAL(KIND = 8), INTENT(IN) :: re_var_36
    REAL(KIND = 8), INTENT(OUT) :: od_var_37(nb_var_33), scat_od_var_38(nb_var_33), g_var_39(nb_var_33)
    REAL(KIND = 8) :: de_um_var_40, inv_de_um_var_41
    REAL(KIND = 8) :: iwp_gm_2_var_42
    INTEGER :: jb_var_43
    de_um_var_40 = MIN(re_var_36, 0.0001D0) * (1539598.4727183152D0)
    inv_de_um_var_41 = 1.0D0 / de_um_var_40
    iwp_gm_2_var_42 = ice_wp_var_35 * 1000.0D0
    DO jb_var_43 = 1, 14
      od_var_37(jb_var_43) = iwp_gm_2_var_42 * (coeff_var_34(jb_var_43, 1) + coeff_var_34(jb_var_43, 2) * inv_de_um_var_41)
      scat_od_var_38(jb_var_43) = od_var_37(jb_var_43) * (1.0D0 - (coeff_var_34(jb_var_43, 3) + de_um_var_40 * (coeff_var_34(jb_var_43, 4) + de_um_var_40 * (coeff_var_34(jb_var_43, 5) + de_um_var_40 * coeff_var_34(jb_var_43, 6)))))
      g_var_39(jb_var_43) = MIN(coeff_var_34(jb_var_43, 7) + de_um_var_40 * (coeff_var_34(jb_var_43, 8) + de_um_var_40 * (coeff_var_34(jb_var_43, 9) + de_um_var_40 * coeff_var_34(jb_var_43, 10))), 0.9999999999999978D0)
    END DO
  END SUBROUTINE calc_ice_optics_fu_sw
  SUBROUTINE calc_ice_optics_fu_lw(nb_var_44, coeff_var_45, ice_wp_var_46, re_var_47, od_var_48, scat_od_var_49, g_var_50)
    INTEGER, INTENT(IN) :: nb_var_44
    REAL(KIND = 8), INTENT(IN) :: coeff_var_45(:, :)
    REAL(KIND = 8), INTENT(IN) :: ice_wp_var_46
    REAL(KIND = 8), INTENT(IN) :: re_var_47
    REAL(KIND = 8), INTENT(OUT) :: od_var_48(nb_var_44), scat_od_var_49(nb_var_44), g_var_50(nb_var_44)
    REAL(KIND = 8) :: de_um_var_51, inv_de_um_var_52
    REAL(KIND = 8) :: iwp_gm_2_var_53
    INTEGER :: jb_var_54
    de_um_var_51 = MIN(re_var_47, 0.0001D0) * (1539598.4727183152D0)
    inv_de_um_var_52 = 1.0D0 / de_um_var_51
    iwp_gm_2_var_53 = ice_wp_var_46 * 1000.0D0
    DO jb_var_54 = 1, 16
      od_var_48(jb_var_54) = iwp_gm_2_var_53 * (coeff_var_45(jb_var_54, 1) + inv_de_um_var_52 * (coeff_var_45(jb_var_54, 2) + inv_de_um_var_52 * coeff_var_45(jb_var_54, 3)))
      scat_od_var_49(jb_var_54) = od_var_48(jb_var_54) - iwp_gm_2_var_53 * inv_de_um_var_52 * (coeff_var_45(jb_var_54, 4) + de_um_var_51 * (coeff_var_45(jb_var_54, 5) + de_um_var_51 * (coeff_var_45(jb_var_54, 6) + de_um_var_51 * coeff_var_45(jb_var_54, 7))))
      g_var_50(jb_var_54) = MIN(coeff_var_45(jb_var_54, 8) + de_um_var_51 * (coeff_var_45(jb_var_54, 9) + de_um_var_51 * (coeff_var_45(jb_var_54, 10) + de_um_var_51 * coeff_var_45(jb_var_54, 11))), 0.9999999999999978D0)
    END DO
  END SUBROUTINE calc_ice_optics_fu_lw
END MODULE radiation_ice_optics_fu
MODULE radiation_ice_optics_yi
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE calc_ice_optics_yi_sw(nb_var_55, coeff_var_56, ice_wp_var_57, re_var_58, od_var_59, scat_od_var_60, g_var_61)
    INTEGER, INTENT(IN) :: nb_var_55
    REAL(KIND = 8), INTENT(IN) :: coeff_var_56(:, :)
    REAL(KIND = 8), INTENT(IN) :: ice_wp_var_57
    REAL(KIND = 8), INTENT(IN) :: re_var_58
    REAL(KIND = 8), INTENT(OUT) :: od_var_59(nb_var_55), scat_od_var_60(nb_var_55), g_var_61(nb_var_55)
    REAL(KIND = 8) :: de_um_var_62
    REAL(KIND = 8) :: iwp_gm_2_var_63
    REAL(KIND = 8) :: wts_1_var_64, wts_2_var_65
    INTEGER(KIND = 4) :: lu_idx_var_66
    de_um_var_62 = re_var_58 * 2000000.0D0
    de_um_var_62 = MAX(de_um_var_62, 10.0D0)
    de_um_var_62 = MIN(de_um_var_62, 119.99D0)
    iwp_gm_2_var_63 = ice_wp_var_57 * 1000.0D0
    lu_idx_var_66 = FLOOR(de_um_var_62 * 0.2D0 - 1.0D0)
    wts_2_var_65 = (de_um_var_62 * 0.2D0 - 1.0D0) - lu_idx_var_66
    wts_1_var_64 = 1.0D0 - wts_2_var_65
    od_var_59 = 0.001D0 * iwp_gm_2_var_63 * (wts_1_var_64 * coeff_var_56(1 : 14, lu_idx_var_66) + wts_2_var_65 * coeff_var_56(1 : 14, lu_idx_var_66 + 1))
    scat_od_var_60 = od_var_59 * (wts_1_var_64 * coeff_var_56(1 : 14, lu_idx_var_66 + 23) + wts_2_var_65 * coeff_var_56(1 : 14, lu_idx_var_66 + 23 + 1))
    g_var_61 = wts_1_var_64 * coeff_var_56(1 : 14, lu_idx_var_66 + 46) + wts_2_var_65 * coeff_var_56(1 : 14, lu_idx_var_66 + 46 + 1)
  END SUBROUTINE calc_ice_optics_yi_sw
  SUBROUTINE calc_ice_optics_yi_lw(nb_var_67, coeff_var_68, ice_wp_var_69, re_var_70, od_var_71, scat_od_var_72, g_var_73)
    INTEGER, INTENT(IN) :: nb_var_67
    REAL(KIND = 8), INTENT(IN) :: coeff_var_68(:, :)
    REAL(KIND = 8), INTENT(IN) :: ice_wp_var_69
    REAL(KIND = 8), INTENT(IN) :: re_var_70
    REAL(KIND = 8), INTENT(OUT) :: od_var_71(nb_var_67), scat_od_var_72(nb_var_67), g_var_73(nb_var_67)
    REAL(KIND = 8) :: de_um_var_74
    REAL(KIND = 8) :: iwp_gm_2_var_75
    REAL(KIND = 8) :: wts_1_var_76, wts_2_var_77
    INTEGER(KIND = 4) :: lu_idx_var_78
    de_um_var_74 = re_var_70 * 2000000.0D0
    de_um_var_74 = MAX(de_um_var_74, 10.0D0)
    de_um_var_74 = MIN(de_um_var_74, 119.99D0)
    iwp_gm_2_var_75 = ice_wp_var_69 * 1000.0D0
    lu_idx_var_78 = FLOOR(de_um_var_74 * 0.2D0 - 1.0D0)
    wts_2_var_77 = (de_um_var_74 * 0.2D0 - 1.0D0) - lu_idx_var_78
    wts_1_var_76 = 1.0D0 - wts_2_var_77
    od_var_71 = 0.001D0 * iwp_gm_2_var_75 * (wts_1_var_76 * coeff_var_68(1 : 16, lu_idx_var_78) + wts_2_var_77 * coeff_var_68(1 : 16, lu_idx_var_78 + 1))
    scat_od_var_72 = od_var_71 * (wts_1_var_76 * coeff_var_68(1 : 16, lu_idx_var_78 + 23) + wts_2_var_77 * coeff_var_68(1 : 16, lu_idx_var_78 + 23 + 1))
    g_var_73 = wts_1_var_76 * coeff_var_68(1 : 16, lu_idx_var_78 + 46) + wts_2_var_77 * coeff_var_68(1 : 16, lu_idx_var_78 + 46 + 1)
  END SUBROUTINE calc_ice_optics_yi_lw
END MODULE radiation_ice_optics_yi
MODULE radiation_liquid_optics_socrates
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE calc_liq_optics_socrates(nb_var_79, coeff_var_80, lwp, re_in, od_var_81, scat_od_var_82, g_var_83)
    INTEGER, INTENT(IN) :: nb_var_79
    REAL(KIND = 8), INTENT(IN) :: coeff_var_80(:, :)
    REAL(KIND = 8), INTENT(IN) :: lwp, re_in
    REAL(KIND = 8), INTENT(OUT) :: od_var_81(nb_var_79), scat_od_var_82(nb_var_79), g_var_83(nb_var_79)
    INTEGER :: jb_var_84
    REAL(KIND = 8) :: re_var_85
    re_var_85 = MAX(1.2000000424450263D-06, MIN(re_in, 4.999999873689376D-05))
    DO jb_var_84 = 1, nb_var_79
      od_var_81(jb_var_84) = lwp * (coeff_var_80(jb_var_84, 1) + re_var_85 * (coeff_var_80(jb_var_84, 2) + re_var_85 * coeff_var_80(jb_var_84, 3))) / (1.0D0 + re_var_85 * (coeff_var_80(jb_var_84, 4) + re_var_85 * (coeff_var_80(jb_var_84, 5) + re_var_85 * coeff_var_80(jb_var_84, 6))))
      scat_od_var_82(jb_var_84) = od_var_81(jb_var_84) * (1.0D0 - (coeff_var_80(jb_var_84, 7) + re_var_85 * (coeff_var_80(jb_var_84, 8) + re_var_85 * coeff_var_80(jb_var_84, 9))) / (1.0D0 + re_var_85 * (coeff_var_80(jb_var_84, 10) + re_var_85 * coeff_var_80(jb_var_84, 11))))
      g_var_83(jb_var_84) = (coeff_var_80(jb_var_84, 12) + re_var_85 * (coeff_var_80(jb_var_84, 13) + re_var_85 * coeff_var_80(jb_var_84, 14))) / (1.0D0 + re_var_85 * (coeff_var_80(jb_var_84, 15) + re_var_85 * coeff_var_80(jb_var_84, 16)))
    END DO
  END SUBROUTINE calc_liq_optics_socrates
END MODULE radiation_liquid_optics_socrates
MODULE radiation_two_stream
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE calc_no_scattering_transmittance_lw(ng_var_86, od_var_87, planck_top, planck_bot, transmittance_var_88, source_up_var_89, source_dn_var_90)
    INTEGER, INTENT(IN) :: ng_var_86
    REAL(KIND = 8), INTENT(IN), DIMENSION(ng_var_86) :: od_var_87
    REAL(KIND = 8), INTENT(IN), DIMENSION(ng_var_86) :: planck_top, planck_bot
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ng_var_86) :: transmittance_var_88
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ng_var_86) :: source_up_var_89, source_dn_var_90
    REAL(KIND = 8) :: coeff_var_91, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot
    INTEGER :: jg_var_92
    transmittance_var_88 = EXP(- 1.66D0 * od_var_87)
    DO jg_var_92 = 1, ng_var_86
      coeff_var_91 = 1.66D0 * od_var_87(jg_var_92)
      IF (od_var_87(jg_var_92) > 0.001D0) THEN
        coeff_var_91 = (planck_bot(jg_var_92) - planck_top(jg_var_92)) / coeff_var_91
        coeff_up_top = coeff_var_91 + planck_top(jg_var_92)
        coeff_up_bot = coeff_var_91 + planck_bot(jg_var_92)
        coeff_dn_top = - coeff_var_91 + planck_top(jg_var_92)
        coeff_dn_bot = - coeff_var_91 + planck_bot(jg_var_92)
        source_up_var_89(jg_var_92) = coeff_up_top - transmittance_var_88(jg_var_92) * coeff_up_bot
        source_dn_var_90(jg_var_92) = coeff_dn_bot - transmittance_var_88(jg_var_92) * coeff_dn_top
      ELSE
        source_up_var_89(jg_var_92) = coeff_var_91 * 0.5D0 * (planck_top(jg_var_92) + planck_bot(jg_var_92))
        source_dn_var_90(jg_var_92) = source_up_var_89(jg_var_92)
      END IF
    END DO
  END SUBROUTINE calc_no_scattering_transmittance_lw
  SUBROUTINE calc_ref_trans_sw(ng_var_93, mu0, od_var_94, ssa_var_95, asymmetry, ref_diff, trans_diff, ref_dir_var_96, trans_dir_diff_var_97, trans_dir_dir_var_98)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng_var_93
    REAL(KIND = 8), INTENT(IN) :: mu0
    REAL(KIND = 8), INTENT(IN), DIMENSION(ng_var_93) :: od_var_94, ssa_var_95, asymmetry
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ng_var_93) :: ref_dir_var_96, trans_dir_diff_var_97
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ng_var_93) :: ref_diff, trans_diff
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ng_var_93) :: trans_dir_dir_var_98
    REAL(KIND = 8), DIMENSION(ng_var_93) :: gamma1, gamma2, gamma3, gamma4
    REAL(KIND = 8), DIMENSION(ng_var_93) :: alpha1, alpha2, k_exponent
    REAL(KIND = 8), DIMENSION(ng_var_93) :: exponential
    REAL(KIND = 8) :: reftrans_factor, factor_var_99
    REAL(KIND = 8) :: exponential2
    REAL(KIND = 8) :: k_mu0, k_gamma3, k_gamma4
    REAL(KIND = 8) :: k_2_exponential, one_minus_kmu0_sqr
    INTEGER :: jg_var_100
    trans_dir_dir_var_98 = MAX(- MAX(od_var_94 * (1.0D0 / mu0), 0.0D0), - 1000.0D0)
    trans_dir_dir_var_98 = EXP(trans_dir_dir_var_98)
    DO jg_var_100 = 1, ng_var_93
      factor_var_99 = 0.75D0 * asymmetry(jg_var_100)
      gamma1(jg_var_100) = 2.0D0 - ssa_var_95(jg_var_100) * (1.25D0 + factor_var_99)
      gamma2(jg_var_100) = ssa_var_95(jg_var_100) * (0.75D0 - factor_var_99)
      gamma3(jg_var_100) = 0.5D0 - mu0 * factor_var_99
      gamma4(jg_var_100) = 1.0D0 - gamma3(jg_var_100)
      alpha1(jg_var_100) = gamma1(jg_var_100) * gamma4(jg_var_100) + gamma2(jg_var_100) * gamma3(jg_var_100)
      alpha2(jg_var_100) = gamma1(jg_var_100) * gamma3(jg_var_100) + gamma2(jg_var_100) * gamma4(jg_var_100)
      k_exponent(jg_var_100) = SQRT(MAX((gamma1(jg_var_100) - gamma2(jg_var_100)) * (gamma1(jg_var_100) + gamma2(jg_var_100)), 1D-12))
    END DO
    exponential = EXP(- k_exponent * od_var_94)
    DO jg_var_100 = 1, ng_var_93
      k_mu0 = k_exponent(jg_var_100) * mu0
      one_minus_kmu0_sqr = 1.0D0 - k_mu0 * k_mu0
      k_gamma3 = k_exponent(jg_var_100) * gamma3(jg_var_100)
      k_gamma4 = k_exponent(jg_var_100) * gamma4(jg_var_100)
      exponential2 = exponential(jg_var_100) * exponential(jg_var_100)
      k_2_exponential = 2.0D0 * k_exponent(jg_var_100) * exponential(jg_var_100)
      reftrans_factor = 1.0D0 / (k_exponent(jg_var_100) + gamma1(jg_var_100) + (k_exponent(jg_var_100) - gamma1(jg_var_100)) * exponential2)
      ref_diff(jg_var_100) = gamma2(jg_var_100) * (1.0D0 - exponential2) * reftrans_factor
      trans_diff(jg_var_100) = MAX(0.0D0, MIN(k_2_exponential * reftrans_factor, 1.0D0 - ref_diff(jg_var_100)))
      reftrans_factor = mu0 * ssa_var_95(jg_var_100) * reftrans_factor / MERGE(one_minus_kmu0_sqr, 2.220446049250313D-16, ABS(one_minus_kmu0_sqr) > 2.220446049250313D-16)
      ref_dir_var_96(jg_var_100) = reftrans_factor * ((1.0D0 - k_mu0) * (alpha2(jg_var_100) + k_gamma3) - (1.0D0 + k_mu0) * (alpha2(jg_var_100) - k_gamma3) * exponential2 - k_2_exponential * (gamma3(jg_var_100) - alpha2(jg_var_100) * mu0) * trans_dir_dir_var_98(jg_var_100))
      trans_dir_diff_var_97(jg_var_100) = reftrans_factor * (k_2_exponential * (gamma4(jg_var_100) + alpha1(jg_var_100) * mu0) - trans_dir_dir_var_98(jg_var_100) * ((1.0D0 + k_mu0) * (alpha1(jg_var_100) + k_gamma4) - (1.0D0 - k_mu0) * (alpha1(jg_var_100) - k_gamma4) * exponential2))
      ref_dir_var_96(jg_var_100) = MAX(0.0D0, MIN(ref_dir_var_96(jg_var_100), mu0 * (1.0D0 - trans_dir_dir_var_98(jg_var_100))))
      trans_dir_diff_var_97(jg_var_100) = MAX(0.0D0, MIN(trans_dir_diff_var_97(jg_var_100), mu0 * (1.0D0 - trans_dir_dir_var_98(jg_var_100)) - ref_dir_var_96(jg_var_100)))
    END DO
  END SUBROUTINE calc_ref_trans_sw
END MODULE radiation_two_stream
MODULE random_numbers_mix
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4), PARAMETER :: jpq = 607
  TYPE :: randomnumberstream
    INTEGER(KIND = 4) :: iused
    INTEGER(KIND = 4) :: inittest
    INTEGER(KIND = 4), DIMENSION(607) :: ix
    REAL(KIND = 8) :: zrm
  END TYPE randomnumberstream
  CONTAINS
  SUBROUTINE initialize_random_numbers(kseed, yd_stream_var_101)
    INTEGER(KIND = 4), INTENT(IN) :: kseed
    TYPE(randomnumberstream), INTENT(INOUT) :: yd_stream_var_101
    INTEGER(KIND = 4) :: idum, jj_var_102, jbit
    REAL(KIND = 8), DIMENSION(999) :: zwarmup
    idum = ABS(IEOR(kseed, 123459876))
    IF (idum == 0) idum = 123459876
    DO jj_var_102 = 1, 64
      IF (BTEST(idum, 31)) THEN
        idum = IBSET(ISHFT(IEOR(idum, 87), 1), 0)
      ELSE
        idum = IBCLR(ISHFT(idum, 1), 0)
      END IF
    END DO
    yd_stream_var_101 % ix(1 : 606) = 0
    yd_stream_var_101 % ix(2) = ISHFT(IBITS(idum, 0, 29), 1)
    yd_stream_var_101 % ix(jpq) = IBITS(idum, 29, BIT_SIZE(idum) + 1 - 30)
    DO jbit = 1, 29
      DO jj_var_102 = 3, 606
        IF (BTEST(idum, 31)) THEN
          idum = IBSET(ISHFT(IEOR(idum, 87), 1), 0)
          yd_stream_var_101 % ix(jj_var_102) = IBSET(yd_stream_var_101 % ix(jj_var_102), jbit)
        ELSE
          idum = IBCLR(ISHFT(idum, 1), 0)
        END IF
      END DO
    END DO
    yd_stream_var_101 % ix(502) = IBSET(yd_stream_var_101 % ix(502), 0)
    yd_stream_var_101 % iused = 607
    yd_stream_var_101 % zrm = 9.313225746154785D-10
    yd_stream_var_101 % inittest = 12345678
    CALL uniform_distribution(zwarmup, yd_stream_var_101)
  END SUBROUTINE initialize_random_numbers
  SUBROUTINE uniform_distribution(px, yd_stream_var_103)
    TYPE(randomnumberstream), INTENT(INOUT) :: yd_stream_var_103
    REAL(KIND = 8), DIMENSION(:), INTENT(OUT) :: px
    INTEGER(KIND = 4) :: jj_var_104, jk_var_105, in_var_106, ifilled
    IF (yd_stream_var_103 % inittest /= 12345678) CALL abor1('uniform_distribution called before initialize_random_numbers')
    in_var_106 = SIZE(px)
    ifilled = 0
    DO jj_var_104 = yd_stream_var_103 % iused + 1, MIN(607, in_var_106 + yd_stream_var_103 % iused)
      px(jj_var_104 - yd_stream_var_103 % iused) = yd_stream_var_103 % ix(jj_var_104) * yd_stream_var_103 % zrm
      ifilled = ifilled + 1
    END DO
    yd_stream_var_103 % iused = yd_stream_var_103 % iused + ifilled
    IF (ifilled == in_var_106) THEN
      RETURN
    END IF
    DO WHILE (ifilled < in_var_106)
      DO jj_var_104 = 1, 273
        yd_stream_var_103 % ix(jj_var_104) = IAND(1073741823, yd_stream_var_103 % ix(jj_var_104) + yd_stream_var_103 % ix(jj_var_104 - 273 + 607))
      END DO
      DO jk_var_105 = 1, 2
        DO jj_var_104 = 274 + (jk_var_105 - 1) * 167, MIN(607, 273 + jk_var_105 * 167)
          yd_stream_var_103 % ix(jj_var_104) = IAND(1073741823, yd_stream_var_103 % ix(jj_var_104) + yd_stream_var_103 % ix(jj_var_104 - 273))
        END DO
      END DO
      yd_stream_var_103 % iused = MIN(607, in_var_106 - ifilled)
      px(ifilled + 1 : ifilled + yd_stream_var_103 % iused) = yd_stream_var_103 % ix(1 : yd_stream_var_103 % iused) * yd_stream_var_103 % zrm
      ifilled = ifilled + yd_stream_var_103 % iused
    END DO
  END SUBROUTINE uniform_distribution
END MODULE random_numbers_mix
MODULE yoerrta1
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_107(10), fracrefb_var_108(10)
  REAL(KIND = 8) :: absa_var_109(65, 10)
  REAL(KIND = 8) :: absb_var_110(235, 10)
  REAL(KIND = 8) :: ka_mn2_var_111(19, 10), kb_mn2(19, 10)
  REAL(KIND = 8) :: selfref_var_112(10, 10), forref_var_113(4, 10)
END MODULE yoerrta1
MODULE yoerrta10
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(6) :: fracrefa_var_114
  REAL(KIND = 8), DIMENSION(6) :: fracrefb_var_115
  REAL(KIND = 8) :: absa_var_116(65, 6)
  REAL(KIND = 8) :: absb_var_117(235, 6)
  REAL(KIND = 8) :: selfref_var_118(10, 6)
  REAL(KIND = 8) :: forref_var_119(4, 6)
END MODULE yoerrta10
MODULE yoerrta11
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(8) :: fracrefa_var_120
  REAL(KIND = 8), DIMENSION(8) :: fracrefb_var_121
  REAL(KIND = 8) :: absa_var_122(65, 8)
  REAL(KIND = 8) :: absb_var_123(235, 8)
  REAL(KIND = 8) :: ka_mo2(19, 8)
  REAL(KIND = 8) :: kb_mo2(19, 8)
  REAL(KIND = 8) :: selfref_var_124(10, 8)
  REAL(KIND = 8) :: forref_var_125(4, 8)
END MODULE yoerrta11
MODULE yoerrta12
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_126(8, 9)
  REAL(KIND = 8) :: absa_var_127(585, 8)
  REAL(KIND = 8) :: selfref_var_128(10, 8)
  REAL(KIND = 8) :: forref_var_129(4, 8)
END MODULE yoerrta12
MODULE yoerrta13
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_130(4, 9)
  REAL(KIND = 8), DIMENSION(4) :: fracrefb_var_131
  REAL(KIND = 8) :: absa_var_132(585, 4)
  REAL(KIND = 8) :: selfref_var_133(10, 4)
  REAL(KIND = 8) :: forref_var_134(4, 4)
  REAL(KIND = 8) :: ka_mco2_var_135(9, 19, 4)
  REAL(KIND = 8) :: ka_mco(9, 19, 4)
  REAL(KIND = 8) :: kb_mo3(19, 4)
END MODULE yoerrta13
MODULE yoerrta14
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(2) :: fracrefa_var_136
  REAL(KIND = 8), DIMENSION(2) :: fracrefb_var_137
  REAL(KIND = 8) :: absa_var_138(65, 2)
  REAL(KIND = 8) :: absb_var_139(235, 2)
  REAL(KIND = 8) :: selfref_var_140(10, 2)
  REAL(KIND = 8) :: forref_var_141(4, 2)
END MODULE yoerrta14
MODULE yoerrta15
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_142(2, 9)
  REAL(KIND = 8) :: absa_var_143(585, 2)
  REAL(KIND = 8) :: ka_mn2_var_144(9, 19, 2)
  REAL(KIND = 8) :: selfref_var_145(10, 2)
  REAL(KIND = 8) :: forref_var_146(4, 2)
END MODULE yoerrta15
MODULE yoerrta16
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_147(2, 9)
  REAL(KIND = 8), DIMENSION(2) :: fracrefb_var_148
  REAL(KIND = 8) :: absa_var_149(585, 2)
  REAL(KIND = 8) :: absb_var_150(235, 2)
  REAL(KIND = 8) :: selfref_var_151(10, 2)
  REAL(KIND = 8) :: forref_var_152(4, 2)
END MODULE yoerrta16
MODULE yoerrta2
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_153(12), fracrefb_var_154(12)
  REAL(KIND = 8) :: absa_var_155(65, 12)
  REAL(KIND = 8) :: absb_var_156(235, 12)
  REAL(KIND = 8) :: selfref_var_157(10, 12), forref_var_158(4, 12)
END MODULE yoerrta2
MODULE yoerrta3
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_159(16, 9), fracrefb_var_160(16, 5)
  REAL(KIND = 8) :: ka_mn2o_var_161(9, 19, 16), kb_mn2o_var_162(5, 19, 16)
  REAL(KIND = 8) :: absa_var_163(585, 16)
  REAL(KIND = 8) :: absb_var_164(1175, 16)
  REAL(KIND = 8) :: selfref_var_165(10, 16)
  REAL(KIND = 8) :: forref_var_166(4, 16)
END MODULE yoerrta3
MODULE yoerrta4
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_167(14, 9), fracrefb_var_168(14, 5)
  REAL(KIND = 8) :: absa_var_169(585, 14)
  REAL(KIND = 8) :: absb_var_170(1175, 14)
  REAL(KIND = 8) :: selfref_var_171(10, 14), forref_var_172(4, 14)
END MODULE yoerrta4
MODULE yoerrta5
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_173(16, 9), fracrefb_var_174(16, 5)
  REAL(KIND = 8), DIMENSION(16) :: ccl4
  REAL(KIND = 8) :: absa_var_175(585, 16)
  REAL(KIND = 8) :: absb_var_176(1175, 16)
  REAL(KIND = 8) :: ka_mo3_var_177(9, 19, 16)
  REAL(KIND = 8) :: selfref_var_178(10, 16)
  REAL(KIND = 8) :: forref_var_179(4, 16)
END MODULE yoerrta5
MODULE yoerrta6
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(8) :: fracrefa_var_180
  REAL(KIND = 8), DIMENSION(8) :: cfc11adj
  REAL(KIND = 8), DIMENSION(8) :: cfc12_var_181
  REAL(KIND = 8) :: absa_var_182(65, 8)
  REAL(KIND = 8) :: selfref_var_183(10, 8)
  REAL(KIND = 8) :: ka_mco2_var_184(19, 8)
  REAL(KIND = 8) :: forref_var_185(4, 8)
END MODULE yoerrta6
MODULE yoerrta7
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_186(12, 9)
  REAL(KIND = 8), DIMENSION(12) :: fracrefb_var_187
  REAL(KIND = 8) :: absa_var_188(585, 12)
  REAL(KIND = 8) :: absb_var_189(235, 12)
  REAL(KIND = 8) :: selfref_var_190(10, 12)
  REAL(KIND = 8) :: ka_mco2_var_191(9, 19, 12)
  REAL(KIND = 8) :: kb_mco2_var_192(19, 12)
  REAL(KIND = 8) :: forref_var_193(4, 12)
END MODULE yoerrta7
MODULE yoerrta8
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(8) :: fracrefa_var_194
  REAL(KIND = 8), DIMENSION(8) :: fracrefb_var_195
  REAL(KIND = 8), DIMENSION(8) :: cfc12_var_196
  REAL(KIND = 8), DIMENSION(8) :: cfc22adj
  REAL(KIND = 8) :: absa_var_197(65, 8)
  REAL(KIND = 8) :: absb_var_198(235, 8)
  REAL(KIND = 8) :: ka_mco2_var_199(19, 8)
  REAL(KIND = 8) :: ka_mn2o_var_200(19, 8)
  REAL(KIND = 8) :: ka_mo3_var_201(19, 8)
  REAL(KIND = 8) :: kb_mco2_var_202(19, 8)
  REAL(KIND = 8) :: kb_mn2o_var_203(19, 8)
  REAL(KIND = 8) :: selfref_var_204(10, 8)
  REAL(KIND = 8) :: forref_var_205(4, 8)
END MODULE yoerrta8
MODULE yoerrta9
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_206(12, 9)
  REAL(KIND = 8), DIMENSION(12) :: fracrefb_var_207
  REAL(KIND = 8) :: absa_var_208(585, 12)
  REAL(KIND = 8) :: absb_var_209(235, 12)
  REAL(KIND = 8) :: ka_mn2o_var_210(9, 19, 12)
  REAL(KIND = 8) :: kb_mn2o_var_211(19, 12)
  REAL(KIND = 8) :: selfref_var_212(10, 12)
  REAL(KIND = 8) :: forref_var_213(4, 12)
END MODULE yoerrta9
MODULE yoerrtm
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4), PARAMETER :: ng3 = 16
  INTEGER(KIND = 4), PARAMETER :: ng4 = 14
  INTEGER(KIND = 4), PARAMETER :: ng5 = 16
  INTEGER(KIND = 4), PARAMETER :: ng7 = 12
  INTEGER(KIND = 4), PARAMETER :: ng9 = 12
  INTEGER(KIND = 4), PARAMETER :: ng12 = 8
  INTEGER(KIND = 4), PARAMETER :: ng13 = 4
  INTEGER(KIND = 4), PARAMETER :: ng15 = 2
  INTEGER(KIND = 4), PARAMETER :: ng16 = 2
END MODULE yoerrtm
MODULE yoerrtrf
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(59) :: preflog_var_214
  REAL(KIND = 8), DIMENSION(59) :: tref_var_215
  REAL(KIND = 8) :: chi_mls(7, 59)
END MODULE yoerrtrf
MODULE yoerrtwn
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4), DIMENSION(16) :: nspa_var_216
  INTEGER(KIND = 4), DIMENSION(16) :: nspb_var_217
  REAL(KIND = 8), DIMENSION(16) :: delwave
  REAL(KIND = 8), DIMENSION(181, 16) :: totplnk
END MODULE yoerrtwn
MODULE yoesrta16
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_218, strrat1
  INTEGER(KIND = 4) :: layreffr_var_219
  REAL(KIND = 8) :: absa_var_220(585, 16)
  REAL(KIND = 8) :: absb_var_221(235, 16)
  REAL(KIND = 8) :: selfrefc_var_222(10, 16), forrefc_var_223(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_224(16)
END MODULE yoesrta16
MODULE yoesrta17
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_225, strrat_var_226
  INTEGER(KIND = 4) :: layreffr_var_227
  REAL(KIND = 8) :: absa_var_228(585, 16)
  REAL(KIND = 8) :: absb_var_229(1175, 16)
  REAL(KIND = 8) :: selfrefc_var_230(10, 16), forrefc_var_231(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_232(16, 5)
END MODULE yoesrta17
MODULE yoesrta18
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_233, strrat_var_234
  INTEGER(KIND = 4) :: layreffr_var_235
  REAL(KIND = 8) :: absa_var_236(585, 16)
  REAL(KIND = 8) :: absb_var_237(235, 16)
  REAL(KIND = 8) :: selfrefc_var_238(10, 16), forrefc_var_239(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_240(16, 9)
END MODULE yoesrta18
MODULE yoesrta19
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_241, strrat_var_242
  INTEGER(KIND = 4) :: layreffr_var_243
  REAL(KIND = 8) :: absa_var_244(585, 16)
  REAL(KIND = 8) :: absb_var_245(235, 16)
  REAL(KIND = 8) :: selfrefc_var_246(10, 16), forrefc_var_247(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_248(16, 9)
END MODULE yoesrta19
MODULE yoesrta20
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_249
  INTEGER(KIND = 4) :: layreffr_var_250
  REAL(KIND = 8) :: absa_var_251(65, 16)
  REAL(KIND = 8) :: absb_var_252(235, 16)
  REAL(KIND = 8) :: selfrefc_var_253(10, 16), forrefc_var_254(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_255(16), absch4c(16)
END MODULE yoesrta20
MODULE yoesrta21
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_256, strrat_var_257
  INTEGER(KIND = 4) :: layreffr_var_258
  REAL(KIND = 8) :: absa_var_259(585, 16)
  REAL(KIND = 8) :: absb_var_260(1175, 16)
  REAL(KIND = 8) :: selfrefc_var_261(10, 16), forrefc_var_262(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_263(16, 9)
END MODULE yoesrta21
MODULE yoesrta22
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_264, strrat_var_265
  INTEGER(KIND = 4) :: layreffr_var_266
  REAL(KIND = 8) :: absa_var_267(585, 16)
  REAL(KIND = 8) :: absb_var_268(235, 16)
  REAL(KIND = 8) :: selfrefc_var_269(10, 16), forrefc_var_270(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_271(16, 9)
END MODULE yoesrta22
MODULE yoesrta23
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: givfac
  INTEGER(KIND = 4) :: layreffr_var_272
  REAL(KIND = 8) :: absa_var_273(65, 16)
  REAL(KIND = 8) :: selfrefc_var_274(10, 16), forrefc_var_275(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_276(16), raylc_var_277(16)
END MODULE yoesrta23
MODULE yoesrta24
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: strrat_var_278
  INTEGER(KIND = 4) :: layreffr_var_279
  REAL(KIND = 8) :: absa_var_280(585, 16)
  REAL(KIND = 8) :: absb_var_281(235, 16)
  REAL(KIND = 8) :: selfrefc_var_282(10, 16), forrefc_var_283(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_284(16, 9)
  REAL(KIND = 8) :: abso3ac_var_285(16), abso3bc_var_286(16), raylac(16, 9), raylbc(16)
END MODULE yoesrta24
MODULE yoesrta25
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4) :: layreffr_var_287
  REAL(KIND = 8) :: absa_var_288(65, 16)
  REAL(KIND = 8) :: sfluxrefc_var_289(16)
  REAL(KIND = 8) :: raylc_var_290(16), abso3ac_var_291(16), abso3bc_var_292(16)
END MODULE yoesrta25
MODULE yoesrta26
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: sfluxrefc_var_293(16), raylc_var_294(16)
END MODULE yoesrta26
MODULE yoesrta27
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: scalekur
  INTEGER(KIND = 4) :: layreffr_var_295
  REAL(KIND = 8) :: absa_var_296(65, 16)
  REAL(KIND = 8) :: absb_var_297(235, 16)
  REAL(KIND = 8) :: sfluxrefc_var_298(16), raylc_var_299(16)
END MODULE yoesrta27
MODULE yoesrta28
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_300, strrat_var_301
  INTEGER(KIND = 4) :: layreffr_var_302
  REAL(KIND = 8) :: absa_var_303(585, 16)
  REAL(KIND = 8) :: absb_var_304(1175, 16)
  REAL(KIND = 8) :: sfluxrefc_var_305(16, 5)
END MODULE yoesrta28
MODULE yoesrta29
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_306
  INTEGER(KIND = 4) :: layreffr_var_307
  REAL(KIND = 8) :: absa_var_308(65, 16)
  REAL(KIND = 8) :: absb_var_309(235, 16)
  REAL(KIND = 8) :: selfrefc_var_310(10, 16), forrefc_var_311(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_312(16), absh2oc(16), absco2c(16)
END MODULE yoesrta29
MODULE yoesrtwn
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4), DIMENSION(16 : 29) :: nspa_var_313
  INTEGER(KIND = 4), DIMENSION(16 : 29) :: nspb_var_314
  REAL(KIND = 8), DIMENSION(59) :: preflog_var_315
  REAL(KIND = 8), DIMENSION(59) :: tref_var_316
  INTEGER(KIND = 4), DIMENSION(14) :: ngc
END MODULE yoesrtwn
MODULE yomlun_ifsaux
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4) :: nulout = 6
  INTEGER(KIND = 4), PARAMETER :: nulerr = 0
END MODULE yomlun_ifsaux
MODULE radiation_io
  USE yomlun_ifsaux, ONLY: nulerr
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE radiation_abort(text)
    USE yomlun_ifsaux, ONLY: nulerr
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: text
    IF (PRESENT(text)) THEN
      WRITE(nulerr, '(a)') text
      ERROR STOP 1
    ELSE
      ERROR STOP 'error in radiation scheme'
    END IF
  END SUBROUTINE radiation_abort
END MODULE radiation_io
MODULE radiation_cloud_optics_data
  IMPLICIT NONE
  TYPE :: cloud_optics_type
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: liq_coeff_lw, liq_coeff_sw, ice_coeff_lw, ice_coeff_sw
  END TYPE cloud_optics_type
  CONTAINS
END MODULE radiation_cloud_optics_data
MODULE radiation_gas
  IMPLICIT NONE
  TYPE :: gas_type
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: mixing_ratio
  END TYPE gas_type
  CONTAINS
  RECURSIVE SUBROUTINE assert_units_gas(this_var_318, iunits, igas, scale_factor, istatus)
    CLASS(gas_type), INTENT(IN) :: this_var_318
    INTEGER, INTENT(IN) :: iunits
    INTEGER, OPTIONAL, INTENT(IN) :: igas
    REAL(KIND = 8), OPTIONAL, INTENT(IN) :: scale_factor
    LOGICAL, OPTIONAL, INTENT(OUT) :: istatus
    REAL(KIND = 8) :: sf
    IF (PRESENT(scale_factor)) THEN
      sf = scale_factor
    ELSE
      sf = 1.0D0
    END IF
    IF (PRESENT(istatus)) THEN
      istatus = .TRUE.
    END IF
  END SUBROUTINE assert_units_gas
END MODULE radiation_gas
MODULE radiation_pdf_sampler
  IMPLICIT NONE
  TYPE :: pdf_sampler_type
    INTEGER :: ncdf, nfsd
    REAL(KIND = 8) :: fsd1, inv_fsd_interval
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: val
  END TYPE pdf_sampler_type
  CONTAINS
  ELEMENTAL SUBROUTINE sample_from_pdf(this_var_320, fsd, cdf, val_var_321)
    CLASS(pdf_sampler_type), INTENT(IN) :: this_var_320
    REAL(KIND = 8), INTENT(IN) :: fsd, cdf
    REAL(KIND = 8), INTENT(OUT) :: val_var_321
    INTEGER :: ifsd, icdf
    REAL(KIND = 8) :: wfsd, wcdf
    wcdf = cdf * (this_var_320 % ncdf - 1) + 1.0D0
    icdf = MAX(1, MIN(INT(wcdf), this_var_320 % ncdf - 1))
    wcdf = MAX(0.0D0, MIN(wcdf - icdf, 1.0D0))
    wfsd = (fsd - this_var_320 % fsd1) * this_var_320 % inv_fsd_interval + 1.0D0
    ifsd = MAX(1, MIN(INT(wfsd), this_var_320 % nfsd - 1))
    wfsd = MAX(0.0D0, MIN(wfsd - ifsd, 1.0D0))
    val_var_321 = (1.0D0 - wcdf) * (1.0D0 - wfsd) * this_var_320 % val(icdf, ifsd) + (1.0D0 - wcdf) * wfsd * this_var_320 % val(icdf, ifsd + 1) + wcdf * (1.0D0 - wfsd) * this_var_320 % val(icdf + 1, ifsd) + wcdf * wfsd * this_var_320 % val(icdf + 1, ifsd + 1)
  END SUBROUTINE sample_from_pdf
END MODULE radiation_pdf_sampler
MODULE radiation_cloud_generator
  CONTAINS
  SUBROUTINE cloud_generator(ng_var_322, nlev_var_323, i_overlap_scheme, iseed_var_324, frac_threshold, frac_var_325, overlap_param_var_326, decorrelation_scaling, fractional_std_var_327, pdf_sampler_var_328, od_scaling_var_329, total_cloud_cover_var_330, use_beta_overlap, use_vectorizable_generator)
    USE radiation_pdf_sampler, ONLY: pdf_sampler_type
    USE random_numbers_mix, ONLY: initialize_random_numbers, randomnumberstream, uniform_distribution
    USE radiation_cloud_cover, ONLY: cum_cloud_cover_exp_ran
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng_var_322
    INTEGER, INTENT(IN) :: nlev_var_323
    INTEGER, INTENT(IN) :: i_overlap_scheme
    INTEGER, INTENT(IN) :: iseed_var_324
    REAL(KIND = 8), INTENT(IN) :: frac_threshold
    REAL(KIND = 8), INTENT(IN) :: frac_var_325(nlev_var_323)
    REAL(KIND = 8), INTENT(IN) :: overlap_param_var_326(nlev_var_323 - 1)
    REAL(KIND = 8), INTENT(IN) :: decorrelation_scaling
    REAL(KIND = 8), INTENT(IN) :: fractional_std_var_327(nlev_var_323)
    TYPE(pdf_sampler_type), INTENT(IN) :: pdf_sampler_var_328
    LOGICAL, INTENT(IN), OPTIONAL :: use_beta_overlap
    LOGICAL, INTENT(IN), OPTIONAL :: use_vectorizable_generator
    REAL(KIND = 8), INTENT(OUT) :: od_scaling_var_329(ng_var_322, nlev_var_323)
    REAL(KIND = 8), INTENT(OUT) :: total_cloud_cover_var_330
    REAL(KIND = 8) :: cum_cloud_cover_var_331(nlev_var_323)
    REAL(KIND = 8) :: trigger
    REAL(KIND = 8) :: rand_top(ng_var_322)
    REAL(KIND = 8) :: overlap_param_inhom_var_332(nlev_var_323 - 1)
    TYPE(randomnumberstream) :: random_stream_var_333
    INTEGER :: ibegin, iend_var_334
    INTEGER :: itrigger_var_335
    INTEGER :: jlev_var_336, jg_var_337
    REAL(KIND = 8), DIMENSION(nlev_var_323 - 1) :: pair_cloud_cover_var_338, overhang_var_339
    LOGICAL :: use_vec_gen
    CALL cum_cloud_cover_exp_ran(nlev_var_323, frac_var_325, overlap_param_var_326, cum_cloud_cover_var_331, pair_cloud_cover_var_338, .FALSE.)
    total_cloud_cover_var_330 = cum_cloud_cover_var_331(nlev_var_323)
    DO jlev_var_336 = 1, nlev_var_323 - 1
      overhang_var_339(jlev_var_336) = cum_cloud_cover_var_331(jlev_var_336 + 1) - cum_cloud_cover_var_331(jlev_var_336)
    END DO
    IF (total_cloud_cover_var_330 < 1D-06) THEN
      total_cloud_cover_var_330 = 0.0D0
    ELSE
      jlev_var_336 = 1
      DO WHILE (frac_var_325(jlev_var_336) <= 0.0D0)
        jlev_var_336 = jlev_var_336 + 1
      END DO
      ibegin = jlev_var_336
      iend_var_334 = jlev_var_336
      DO jlev_var_336 = jlev_var_336 + 1, nlev_var_323
        IF (frac_var_325(jlev_var_336) > 0.0D0) THEN
          iend_var_334 = jlev_var_336
        END IF
      END DO
      overlap_param_inhom_var_332 = overlap_param_var_326
      DO jlev_var_336 = ibegin, iend_var_334 - 1
        IF (overlap_param_var_326(jlev_var_336) > 0.0D0) THEN
          overlap_param_inhom_var_332(jlev_var_336) = overlap_param_var_326(jlev_var_336) ** (2.0D0)
        END IF
      END DO
      od_scaling_var_329 = 0.0D0
      use_vec_gen = .FALSE.
      CALL initialize_random_numbers(iseed_var_324, random_stream_var_333)
      CALL uniform_distribution(rand_top, random_stream_var_333)
      DO jg_var_337 = 1, ng_var_322
        trigger = rand_top(jg_var_337) * total_cloud_cover_var_330
        jlev_var_336 = ibegin
        DO WHILE (trigger > cum_cloud_cover_var_331(jlev_var_336) .AND. jlev_var_336 < iend_var_334)
          jlev_var_336 = jlev_var_336 + 1
        END DO
        itrigger_var_335 = jlev_var_336
        CALL generate_column_exp_ran(ng_var_322, nlev_var_323, jg_var_337, random_stream_var_333, pdf_sampler_var_328, frac_var_325, pair_cloud_cover_var_338, cum_cloud_cover_var_331, overhang_var_339, fractional_std_var_327, overlap_param_inhom_var_332, itrigger_var_335, iend_var_334, od_scaling_var_329)
      END DO
    END IF
  END SUBROUTINE cloud_generator
  SUBROUTINE generate_column_exp_ran(ng_var_340, nlev_var_342, ig_var_341, random_stream_var_343, pdf_sampler_var_344, frac_var_345, pair_cloud_cover_var_348, cum_cloud_cover_var_346, overhang_var_349, fractional_std_var_347, overlap_param_inhom_var_350, itrigger_var_351, iend_var_352, od_scaling_var_353)
    USE random_numbers_mix, ONLY: randomnumberstream, uniform_distribution
    USE radiation_pdf_sampler, ONLY: pdf_sampler_type, sample_from_pdf
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng_var_340, ig_var_341
    INTEGER, INTENT(IN) :: nlev_var_342
    TYPE(randomnumberstream), INTENT(INOUT) :: random_stream_var_343
    TYPE(pdf_sampler_type), INTENT(IN) :: pdf_sampler_var_344
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_342) :: frac_var_345, cum_cloud_cover_var_346, fractional_std_var_347
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_342 - 1) :: pair_cloud_cover_var_348, overhang_var_349
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_342 - 1) :: overlap_param_inhom_var_350
    INTEGER, INTENT(IN) :: itrigger_var_351, iend_var_352
    REAL(KIND = 8), INTENT(INOUT), DIMENSION(ng_var_340, nlev_var_342) :: od_scaling_var_353
    INTEGER :: jlev_var_354, jcloud
    INTEGER :: n_layers_to_scale
    INTEGER :: iy
    LOGICAL :: do_fill_od_scaling
    REAL(KIND = 8) :: rand_cloud(nlev_var_342)
    REAL(KIND = 8) :: rand_inhom1(nlev_var_342), rand_inhom2(nlev_var_342)
    n_layers_to_scale = 1
    iy = 0
    CALL uniform_distribution(rand_cloud(1 : (iend_var_352 + 1 - itrigger_var_351)), random_stream_var_343)
    DO jlev_var_354 = itrigger_var_351 + 1, iend_var_352 + 1
      do_fill_od_scaling = .FALSE.
      IF (jlev_var_354 <= iend_var_352) THEN
        iy = iy + 1
        IF (n_layers_to_scale > 0) THEN
          IF (rand_cloud(iy) * frac_var_345(jlev_var_354 - 1) < frac_var_345(jlev_var_354) + frac_var_345(jlev_var_354 - 1) - pair_cloud_cover_var_348(jlev_var_354 - 1)) THEN
            n_layers_to_scale = n_layers_to_scale + 1
          ELSE
            do_fill_od_scaling = .TRUE.
          END IF
        ELSE
          IF (rand_cloud(iy) * (cum_cloud_cover_var_346(jlev_var_354 - 1) - frac_var_345(jlev_var_354 - 1)) < pair_cloud_cover_var_348(jlev_var_354 - 1) - overhang_var_349(jlev_var_354 - 1) - frac_var_345(jlev_var_354 - 1)) THEN
            n_layers_to_scale = 1
          END IF
        END IF
      ELSE
        do_fill_od_scaling = .TRUE.
      END IF
      IF (do_fill_od_scaling) THEN
        CALL uniform_distribution(rand_inhom1(1 : n_layers_to_scale), random_stream_var_343)
        CALL uniform_distribution(rand_inhom2(1 : n_layers_to_scale), random_stream_var_343)
        DO jcloud = 2, n_layers_to_scale
          IF (rand_inhom2(jcloud) < overlap_param_inhom_var_350(jlev_var_354 - n_layers_to_scale + jcloud - 2)) THEN
            rand_inhom1(jcloud) = rand_inhom1(jcloud - 1)
          END IF
        END DO
        CALL sample_from_pdf(pdf_sampler_var_344, fractional_std_var_347(jlev_var_354 - n_layers_to_scale : jlev_var_354 - 1), rand_inhom1(1 : n_layers_to_scale), od_scaling_var_353(ig_var_341, jlev_var_354 - n_layers_to_scale : jlev_var_354 - 1))
        n_layers_to_scale = 0
      END IF
    END DO
  END SUBROUTINE generate_column_exp_ran
END MODULE radiation_cloud_generator
MODULE radiation_config
  USE radiation_cloud_optics_data, ONLY: cloud_optics_type
  USE radiation_pdf_sampler, ONLY: pdf_sampler_type
  IMPLICIT NONE
  TYPE :: config_type
    INTEGER, ALLOCATABLE, DIMENSION(:) :: i_emiss_from_band_lw
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: sw_albedo_weights
    INTEGER, ALLOCATABLE, DIMENSION(:) :: i_band_from_g_lw
    INTEGER, ALLOCATABLE, DIMENSION(:) :: i_band_from_reordered_g_lw
    INTEGER, ALLOCATABLE, DIMENSION(:) :: i_band_from_reordered_g_sw
    TYPE(cloud_optics_type) :: cloud_optics
    TYPE(pdf_sampler_type) :: pdf_sampler
  END TYPE config_type
  CONTAINS
END MODULE radiation_config
MODULE radiation_aerosol
  IMPLICIT NONE
  TYPE :: aerosol_type
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: od_sw, ssa_sw, g_sw, od_lw, ssa_lw
  END TYPE aerosol_type
  CONTAINS
END MODULE radiation_aerosol
MODULE radiation_flux
  USE radiation_config, ONLY: config_type
  IMPLICIT NONE
  TYPE :: flux_type
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: lw_up, lw_dn, sw_up, sw_dn, sw_dn_direct, lw_up_clear, lw_dn_clear, sw_up_clear, sw_dn_clear, sw_dn_direct_clear
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: lw_dn_surf_g, lw_dn_surf_clear_g, sw_dn_diffuse_surf_g, sw_dn_direct_surf_g, sw_dn_diffuse_surf_clear_g, sw_dn_direct_surf_clear_g
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: sw_dn_surf_band, sw_dn_direct_surf_band, sw_dn_surf_clear_band, sw_dn_direct_surf_clear_band
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: cloud_cover_lw, cloud_cover_sw
  END TYPE flux_type
  CONTAINS
  SUBROUTINE calc_surface_spectral(this_var_361, config_var_362, istartcol_var_363, iendcol_var_364)
    USE radiation_config, ONLY: config_type
    CLASS(flux_type), INTENT(INOUT) :: this_var_361
    TYPE(config_type), INTENT(IN) :: config_var_362
    INTEGER, INTENT(IN) :: istartcol_var_363, iendcol_var_364
    INTEGER :: jcol_var_365
    DO jcol_var_365 = istartcol_var_363, iendcol_var_364
      CALL indexed_sum(this_var_361 % sw_dn_direct_surf_g(:, jcol_var_365), config_var_362 % i_band_from_reordered_g_sw, this_var_361 % sw_dn_direct_surf_band(:, jcol_var_365))
      CALL indexed_sum(this_var_361 % sw_dn_diffuse_surf_g(:, jcol_var_365), config_var_362 % i_band_from_reordered_g_sw, this_var_361 % sw_dn_surf_band(:, jcol_var_365))
      this_var_361 % sw_dn_surf_band(:, jcol_var_365) = this_var_361 % sw_dn_surf_band(:, jcol_var_365) + this_var_361 % sw_dn_direct_surf_band(:, jcol_var_365)
    END DO
    DO jcol_var_365 = istartcol_var_363, iendcol_var_364
      CALL indexed_sum(this_var_361 % sw_dn_direct_surf_clear_g(:, jcol_var_365), config_var_362 % i_band_from_reordered_g_sw, this_var_361 % sw_dn_direct_surf_clear_band(:, jcol_var_365))
      CALL indexed_sum(this_var_361 % sw_dn_diffuse_surf_clear_g(:, jcol_var_365), config_var_362 % i_band_from_reordered_g_sw, this_var_361 % sw_dn_surf_clear_band(:, jcol_var_365))
      this_var_361 % sw_dn_surf_clear_band(:, jcol_var_365) = this_var_361 % sw_dn_surf_clear_band(:, jcol_var_365) + this_var_361 % sw_dn_direct_surf_clear_band(:, jcol_var_365)
    END DO
  END SUBROUTINE calc_surface_spectral
  SUBROUTINE calc_toa_spectral(this_var_366, config_var_367, istartcol_var_368, iendcol_var_369)
    USE radiation_config, ONLY: config_type
    CLASS(flux_type), INTENT(INOUT) :: this_var_366
    TYPE(config_type), INTENT(IN) :: config_var_367
    INTEGER, INTENT(IN) :: istartcol_var_368, iendcol_var_369
  END SUBROUTINE calc_toa_spectral
  PURE SUBROUTINE indexed_sum(source_var_370, ind_var_371, dest)
    REAL(KIND = 8), INTENT(IN) :: source_var_370(:)
    INTEGER, INTENT(IN) :: ind_var_371(:)
    REAL(KIND = 8), INTENT(OUT) :: dest(:)
    INTEGER :: ig_var_372, jg_var_373, istart, iend_var_374
    dest = 0.0
    istart = LBOUND(source_var_370, 1)
    iend_var_374 = UBOUND(source_var_370, 1)
    DO jg_var_373 = istart, iend_var_374
      ig_var_372 = ind_var_371(jg_var_373)
      dest(ig_var_372) = dest(ig_var_372) + source_var_370(jg_var_373)
    END DO
  END SUBROUTINE indexed_sum
END MODULE radiation_flux
MODULE radiation_single_level
  USE radiation_config, ONLY: config_type
  USE yomlun_ifsaux, ONLY: nulerr
  USE radiation_io, ONLY: radiation_abort
  IMPLICIT NONE
  TYPE :: single_level_type
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: cos_sza, skin_temperature
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: sw_albedo, sw_albedo_direct
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: lw_emissivity
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: spectral_solar_scaling
    INTEGER, ALLOCATABLE, DIMENSION(:) :: iseed
  END TYPE single_level_type
  CONTAINS
  SUBROUTINE get_albedos(this_var_378, istartcol_var_380, iendcol_var_381, config_var_379, sw_albedo_direct_var_383, sw_albedo_diffuse_var_384, lw_albedo_var_382)
    USE radiation_config, ONLY: config_type
    USE yomlun_ifsaux, ONLY: nulerr
    USE radiation_io, ONLY: radiation_abort
    CLASS(single_level_type), INTENT(IN) :: this_var_378
    TYPE(config_type), INTENT(IN) :: config_var_379
    INTEGER, INTENT(IN) :: istartcol_var_380, iendcol_var_381
    REAL(KIND = 8), INTENT(OUT), OPTIONAL :: lw_albedo_var_382(140, istartcol_var_380 : iendcol_var_381)
    REAL(KIND = 8), INTENT(OUT), DIMENSION(112, istartcol_var_380 : iendcol_var_381) :: sw_albedo_direct_var_383, sw_albedo_diffuse_var_384
    REAL(KIND = 8) :: sw_albedo_band(istartcol_var_380 : iendcol_var_381, 14)
    INTEGER :: nalbedoband
    INTEGER :: jband_var_385, jalbedoband, jcol_var_386
    nalbedoband = SIZE(config_var_379 % sw_albedo_weights, 1)
    IF (SIZE(this_var_378 % sw_albedo, 2) /= nalbedoband) THEN
      WRITE(nulerr, '(a,i0,a)') '*** error: single_level%sw_albedo does not have the expected ', nalbedoband, ' bands'
      CALL radiation_abort
    END IF
    DO jband_var_385 = 1, 14
      DO jcol_var_386 = istartcol_var_380, iendcol_var_381
        sw_albedo_band(jcol_var_386, jband_var_385) = 0.0D0
      END DO
    END DO
    DO jband_var_385 = 1, 14
      DO jalbedoband = 1, nalbedoband
        IF (config_var_379 % sw_albedo_weights(jalbedoband, jband_var_385) /= 0.0D0) THEN
          DO jcol_var_386 = istartcol_var_380, iendcol_var_381
            sw_albedo_band(jcol_var_386, jband_var_385) = sw_albedo_band(jcol_var_386, jband_var_385) + config_var_379 % sw_albedo_weights(jalbedoband, jband_var_385) * this_var_378 % sw_albedo(jcol_var_386, jalbedoband)
          END DO
        END IF
      END DO
    END DO
    sw_albedo_diffuse_var_384 = TRANSPOSE(sw_albedo_band(istartcol_var_380 : iendcol_var_381, config_var_379 % i_band_from_reordered_g_sw))
    DO jband_var_385 = 1, 14
      DO jcol_var_386 = istartcol_var_380, iendcol_var_381
        sw_albedo_band(jcol_var_386, jband_var_385) = 0.0D0
      END DO
    END DO
    DO jband_var_385 = 1, 14
      DO jalbedoband = 1, nalbedoband
        IF (config_var_379 % sw_albedo_weights(jalbedoband, jband_var_385) /= 0.0D0) THEN
          DO jcol_var_386 = istartcol_var_380, iendcol_var_381
            sw_albedo_band(jcol_var_386, jband_var_385) = sw_albedo_band(jcol_var_386, jband_var_385) + config_var_379 % sw_albedo_weights(jalbedoband, jband_var_385) * this_var_378 % sw_albedo_direct(jcol_var_386, jalbedoband)
          END DO
        END IF
      END DO
    END DO
    sw_albedo_direct_var_383 = TRANSPOSE(sw_albedo_band(istartcol_var_380 : iendcol_var_381, config_var_379 % i_band_from_reordered_g_sw))
    IF (MAXVAL(config_var_379 % i_emiss_from_band_lw) > SIZE(this_var_378 % lw_emissivity, 2)) THEN
      WRITE(nulerr, '(a,i0,a)') '*** error: single_level%lw_emissivity has fewer than required ', MAXVAL(config_var_379 % i_emiss_from_band_lw), ' bands'
      CALL radiation_abort
    END IF
    lw_albedo_var_382 = 1.0D0 - TRANSPOSE(this_var_378 % lw_emissivity(istartcol_var_380 : iendcol_var_381, config_var_379 % i_emiss_from_band_lw(config_var_379 % i_band_from_reordered_g_lw)))
  END SUBROUTINE get_albedos
END MODULE radiation_single_level
MODULE radiation_thermodynamics
  IMPLICIT NONE
  TYPE :: thermodynamics_type
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: pressure_hl, temperature_hl, pressure_fl, temperature_fl
  END TYPE thermodynamics_type
  CONTAINS
END MODULE radiation_thermodynamics
MODULE radiation_aerosol_optics
  USE radiation_config, ONLY: config_type
  USE radiation_thermodynamics, ONLY: thermodynamics_type
  USE radiation_gas, ONLY: gas_type
  USE radiation_aerosol, ONLY: aerosol_type
  USE yomlun_ifsaux, ONLY: nulerr
  USE radiation_io, ONLY: radiation_abort
  IMPLICIT NONE
  CONTAINS
  ELEMENTAL SUBROUTINE delta_eddington_extensive(od_var_389, scat_od_var_390, scat_od_g)
    REAL(KIND = 8), INTENT(INOUT) :: od_var_389, scat_od_var_390, scat_od_g
    REAL(KIND = 8) :: f_var_391, g_var_392
    IF (scat_od_var_390 > 0.0D0) THEN
      g_var_392 = scat_od_g / scat_od_var_390
    ELSE
      g_var_392 = 0.0
    END IF
    f_var_391 = g_var_392 * g_var_392
    od_var_389 = od_var_389 - scat_od_var_390 * f_var_391
    scat_od_var_390 = scat_od_var_390 * (1.0D0 - f_var_391)
    scat_od_g = scat_od_var_390 * g_var_392 / (1.0D0 + g_var_392)
  END SUBROUTINE delta_eddington_extensive
  SUBROUTINE add_aerosol_optics(nlev_var_393, istartcol_var_394, iendcol_var_395, config_var_396, thermodynamics_var_397, gas_var_398, aerosol_var_399, od_lw_var_400, ssa_lw_var_401, g_lw_var_402, od_sw_var_403, ssa_sw_var_404, g_sw_var_405)
    USE radiation_config, ONLY: config_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE radiation_gas, ONLY: gas_type
    USE radiation_aerosol, ONLY: aerosol_type
    INTEGER, INTENT(IN) :: nlev_var_393
    INTEGER, INTENT(IN) :: istartcol_var_394, iendcol_var_395
    TYPE(config_type), INTENT(IN), TARGET :: config_var_396
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics_var_397
    TYPE(gas_type), INTENT(IN) :: gas_var_398
    TYPE(aerosol_type), INTENT(IN) :: aerosol_var_399
    REAL(KIND = 8), DIMENSION(140, nlev_var_393, istartcol_var_394 : iendcol_var_395), INTENT(INOUT) :: od_lw_var_400
    REAL(KIND = 8), DIMENSION(0, nlev_var_393, istartcol_var_394 : iendcol_var_395), INTENT(OUT) :: ssa_lw_var_401, g_lw_var_402
    REAL(KIND = 8), DIMENSION(112, nlev_var_393, istartcol_var_394 : iendcol_var_395), INTENT(INOUT) :: od_sw_var_403, ssa_sw_var_404
    REAL(KIND = 8), DIMENSION(112, nlev_var_393, istartcol_var_394 : iendcol_var_395), INTENT(OUT) :: g_sw_var_405
    CALL add_aerosol_optics_direct(nlev_var_393, istartcol_var_394, iendcol_var_395, config_var_396, aerosol_var_399, od_lw_var_400, ssa_lw_var_401, g_lw_var_402, od_sw_var_403, ssa_sw_var_404, g_sw_var_405)
  END SUBROUTINE add_aerosol_optics
  SUBROUTINE add_aerosol_optics_direct(nlev_var_406, istartcol_var_407, iendcol_var_408, config_var_409, aerosol_var_410, od_lw_var_411, ssa_lw_var_412, g_lw_var_413, od_sw_var_414, ssa_sw_var_415, g_sw_var_416)
    USE radiation_config, ONLY: config_type
    USE radiation_aerosol, ONLY: aerosol_type
    USE yomlun_ifsaux, ONLY: nulerr
    USE radiation_io, ONLY: radiation_abort
    INTEGER, INTENT(IN) :: nlev_var_406
    INTEGER, INTENT(IN) :: istartcol_var_407, iendcol_var_408
    TYPE(config_type), INTENT(IN), TARGET :: config_var_409
    TYPE(aerosol_type), INTENT(IN) :: aerosol_var_410
    REAL(KIND = 8), DIMENSION(140, nlev_var_406, istartcol_var_407 : iendcol_var_408), INTENT(INOUT) :: od_lw_var_411
    REAL(KIND = 8), DIMENSION(0, nlev_var_406, istartcol_var_407 : iendcol_var_408), INTENT(OUT) :: ssa_lw_var_412, g_lw_var_413
    REAL(KIND = 8), DIMENSION(112, nlev_var_406, istartcol_var_407 : iendcol_var_408), INTENT(INOUT) :: od_sw_var_414, ssa_sw_var_415
    REAL(KIND = 8), DIMENSION(112, nlev_var_406, istartcol_var_407 : iendcol_var_408), INTENT(OUT) :: g_sw_var_416
    REAL(KIND = 8) :: local_od, local_scat
    REAL(KIND = 8), DIMENSION(14, nlev_var_406) :: od_sw_aerosol, scat_sw_aerosol, scat_g_sw_aerosol
    REAL(KIND = 8), DIMENSION(16, nlev_var_406) :: od_lw_aerosol
    INTEGER :: jcol_var_417, jlev_var_418, jg_var_419, jb_var_420
    INTEGER :: istartlev_var_421, iendlev_var_422
    INTEGER :: iband_var_423
    IF (UBOUND(aerosol_var_410 % od_sw, 1) /= 14) THEN
      WRITE(nulerr, '(a,i0,a,i0)') '*** error: aerosol%od_sw contains ', UBOUND(aerosol_var_410 % od_sw, 1), ' band, expected ', 14
      CALL radiation_abort
    END IF
    istartlev_var_421 = LBOUND(aerosol_var_410 % od_sw, 2)
    iendlev_var_422 = UBOUND(aerosol_var_410 % od_sw, 2)
    DO jcol_var_417 = istartcol_var_407, iendcol_var_408
      DO jlev_var_418 = 1, nlev_var_406
        DO jg_var_419 = 1, 112
          g_sw_var_416(jg_var_419, jlev_var_418, jcol_var_417) = 0.0D0
        END DO
      END DO
    END DO
    DO jcol_var_417 = istartcol_var_407, iendcol_var_408
      DO jlev_var_418 = istartlev_var_421, iendlev_var_422
        DO jb_var_420 = 1, 14
          od_sw_aerosol(jb_var_420, jlev_var_418) = aerosol_var_410 % od_sw(jb_var_420, jlev_var_418, jcol_var_417)
          scat_sw_aerosol(jb_var_420, jlev_var_418) = aerosol_var_410 % ssa_sw(jb_var_420, jlev_var_418, jcol_var_417) * od_sw_aerosol(jb_var_420, jlev_var_418)
          scat_g_sw_aerosol(jb_var_420, jlev_var_418) = aerosol_var_410 % g_sw(jb_var_420, jlev_var_418, jcol_var_417) * scat_sw_aerosol(jb_var_420, jlev_var_418)
          CALL delta_eddington_extensive(od_sw_aerosol(jb_var_420, jlev_var_418), scat_sw_aerosol(jb_var_420, jlev_var_418), scat_g_sw_aerosol(jb_var_420, jlev_var_418))
        END DO
      END DO
      DO jlev_var_418 = istartlev_var_421, iendlev_var_422
        IF (od_sw_aerosol(1, jlev_var_418) > 0.0D0) THEN
          DO jg_var_419 = 1, 112
            iband_var_423 = config_var_409 % i_band_from_reordered_g_sw(jg_var_419)
            local_od = od_sw_var_414(jg_var_419, jlev_var_418, jcol_var_417) + od_sw_aerosol(iband_var_423, jlev_var_418)
            local_scat = ssa_sw_var_415(jg_var_419, jlev_var_418, jcol_var_417) * od_sw_var_414(jg_var_419, jlev_var_418, jcol_var_417) + scat_sw_aerosol(iband_var_423, jlev_var_418)
            g_sw_var_416(jg_var_419, jlev_var_418, jcol_var_417) = scat_g_sw_aerosol(iband_var_423, jlev_var_418) / local_scat
            local_od = od_sw_var_414(jg_var_419, jlev_var_418, jcol_var_417) + od_sw_aerosol(iband_var_423, jlev_var_418)
            ssa_sw_var_415(jg_var_419, jlev_var_418, jcol_var_417) = local_scat / local_od
            od_sw_var_414(jg_var_419, jlev_var_418, jcol_var_417) = local_od
          END DO
        END IF
      END DO
    END DO
    IF (UBOUND(aerosol_var_410 % od_lw, 1) /= 16) THEN
      WRITE(nulerr, '(a,i0,a,i0)') '*** error: aerosol%od_lw contains ', UBOUND(aerosol_var_410 % od_lw, 1), ' band, expected ', 16
      CALL radiation_abort
    END IF
    istartlev_var_421 = LBOUND(aerosol_var_410 % od_lw, 2)
    iendlev_var_422 = UBOUND(aerosol_var_410 % od_lw, 2)
    DO jcol_var_417 = istartcol_var_407, iendcol_var_408
      DO jlev_var_418 = istartlev_var_421, iendlev_var_422
        DO jb_var_420 = 1, 16
          od_lw_aerosol(jb_var_420, jlev_var_418) = aerosol_var_410 % od_lw(jb_var_420, jlev_var_418, jcol_var_417) * (1.0D0 - aerosol_var_410 % ssa_lw(jb_var_420, jlev_var_418, jcol_var_417))
        END DO
      END DO
      DO jlev_var_418 = istartlev_var_421, iendlev_var_422
        DO jg_var_419 = 1, 140
          od_lw_var_411(jg_var_419, jlev_var_418, jcol_var_417) = od_lw_var_411(jg_var_419, jlev_var_418, jcol_var_417) + od_lw_aerosol(config_var_409 % i_band_from_reordered_g_lw(jg_var_419), jlev_var_418)
        END DO
      END DO
    END DO
  END SUBROUTINE add_aerosol_optics_direct
END MODULE radiation_aerosol_optics
MODULE radiation_cloud
  IMPLICIT NONE
  TYPE :: cloud_type
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: mixing_ratio
    REAL(KIND = 8), POINTER, DIMENSION(:, :) :: q_liq, q_ice, re_liq, re_ice
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: fraction
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: fractional_std
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: overlap_param
  END TYPE cloud_type
  CONTAINS
  SUBROUTINE crop_cloud_fraction(this_var_427, istartcol_var_428, iendcol_var_429, cloud_fraction_threshold, cloud_mixing_ratio_threshold)
    CLASS(cloud_type), INTENT(INOUT) :: this_var_427
    INTEGER, INTENT(IN) :: istartcol_var_428, iendcol_var_429
    INTEGER :: nlev_var_430
    INTEGER :: jcol_var_431, jlev_var_432, jh
    REAL(KIND = 8) :: cloud_fraction_threshold, cloud_mixing_ratio_threshold
    REAL(KIND = 8) :: sum_mixing_ratio(istartcol_var_428 : iendcol_var_429)
    nlev_var_430 = SIZE(this_var_427 % fraction, 2)
    DO jlev_var_432 = 1, nlev_var_430
      DO jcol_var_431 = istartcol_var_428, iendcol_var_429
        sum_mixing_ratio(jcol_var_431) = 0.0D0
      END DO
      DO jh = 1, 2
        DO jcol_var_431 = istartcol_var_428, iendcol_var_429
          sum_mixing_ratio(jcol_var_431) = sum_mixing_ratio(jcol_var_431) + this_var_427 % mixing_ratio(jcol_var_431, jlev_var_432, jh)
        END DO
      END DO
      DO jcol_var_431 = istartcol_var_428, iendcol_var_429
        IF (this_var_427 % fraction(jcol_var_431, jlev_var_432) < cloud_fraction_threshold .OR. sum_mixing_ratio(jcol_var_431) < cloud_mixing_ratio_threshold) THEN
          this_var_427 % fraction(jcol_var_431, jlev_var_432) = 0.0D0
        END IF
      END DO
    END DO
  END SUBROUTINE crop_cloud_fraction
END MODULE radiation_cloud
MODULE radiation_cloud_optics
  USE radiation_config, ONLY: config_type
  USE radiation_thermodynamics, ONLY: thermodynamics_type
  USE radiation_cloud, ONLY: cloud_type
  USE radiation_liquid_optics_socrates, ONLY: calc_liq_optics_socrates
  USE radiation_ice_optics_fu, ONLY: calc_ice_optics_fu_lw, calc_ice_optics_fu_sw
  USE radiation_ice_optics_yi, ONLY: calc_ice_optics_yi_lw, calc_ice_optics_yi_sw
  IMPLICIT NONE
  CONTAINS
  ELEMENTAL SUBROUTINE delta_eddington_scat_od(od_var_433, scat_od_var_434, g_var_435)
    REAL(KIND = 8), INTENT(INOUT) :: od_var_433, scat_od_var_434, g_var_435
    REAL(KIND = 8) :: f_var_436
    f_var_436 = g_var_435 * g_var_435
    od_var_433 = od_var_433 - scat_od_var_434 * f_var_436
    scat_od_var_434 = scat_od_var_434 * (1.0D0 - f_var_436)
    g_var_435 = g_var_435 / (1.0D0 + g_var_435)
  END SUBROUTINE delta_eddington_scat_od
  SUBROUTINE cloud_optics_fn_438(nlev_var_437, istartcol_var_438, iendcol_var_439, config_var_440, thermodynamics_var_441, cloud_var_442, od_lw_cloud_var_443, ssa_lw_cloud_var_444, g_lw_cloud_var_445, od_sw_cloud_var_446, ssa_sw_cloud_var_447, g_sw_cloud_var_448)
    USE radiation_config, ONLY: config_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE radiation_cloud, ONLY: cloud_type
    USE radiation_liquid_optics_socrates, ONLY: calc_liq_optics_socrates
    USE radiation_ice_optics_fu, ONLY: calc_ice_optics_fu_lw, calc_ice_optics_fu_sw
    USE radiation_ice_optics_yi, ONLY: calc_ice_optics_yi_lw, calc_ice_optics_yi_sw
    INTEGER, INTENT(IN) :: nlev_var_437
    INTEGER, INTENT(IN) :: istartcol_var_438, iendcol_var_439
    TYPE(config_type), INTENT(IN), TARGET :: config_var_440
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics_var_441
    TYPE(cloud_type), INTENT(IN) :: cloud_var_442
    REAL(KIND = 8), DIMENSION(16, nlev_var_437, istartcol_var_438 : iendcol_var_439), INTENT(OUT) :: od_lw_cloud_var_443
    REAL(KIND = 8), DIMENSION(0, nlev_var_437, istartcol_var_438 : iendcol_var_439), INTENT(OUT) :: ssa_lw_cloud_var_444, g_lw_cloud_var_445
    REAL(KIND = 8), DIMENSION(14, nlev_var_437, istartcol_var_438 : iendcol_var_439), INTENT(OUT) :: od_sw_cloud_var_446, ssa_sw_cloud_var_447, g_sw_cloud_var_448
    REAL(KIND = 8), DIMENSION(16) :: od_lw_liq, scat_od_lw_liq, g_lw_liq, od_lw_ice, scat_od_lw_ice, g_lw_ice
    REAL(KIND = 8), DIMENSION(14) :: od_sw_liq, scat_od_sw_liq, g_sw_liq, od_sw_ice, scat_od_sw_ice, g_sw_ice
    REAL(KIND = 8) :: lwp_in_cloud, iwp_in_cloud
    REAL(KIND = 8) :: factor_var_449
    INTEGER :: jcol_var_450, jlev_var_451, jb_var_452
    DO jcol_var_450 = istartcol_var_438, iendcol_var_439
      DO jlev_var_451 = 1, nlev_var_437
        DO jb_var_452 = 1, 14
          od_sw_cloud_var_446(jb_var_452, jlev_var_451, jcol_var_450) = 0.0D0
          ssa_sw_cloud_var_447(jb_var_452, jlev_var_451, jcol_var_450) = 0.0D0
          g_sw_cloud_var_448(jb_var_452, jlev_var_451, jcol_var_450) = 0.0D0
        END DO
        DO jb_var_452 = 1, 16
          od_lw_cloud_var_443(jb_var_452, jlev_var_451, jcol_var_450) = 0.0D0
        END DO
        DO jb_var_452 = 1, 0
          ssa_lw_cloud_var_444(jb_var_452, jlev_var_451, jcol_var_450) = 0.0D0
          g_lw_cloud_var_445(jb_var_452, jlev_var_451, jcol_var_450) = 0.0D0
        END DO
      END DO
    END DO
    DO jlev_var_451 = 1, nlev_var_437
      DO jcol_var_450 = istartcol_var_438, iendcol_var_439
        IF (cloud_var_442 % fraction(jcol_var_450, jlev_var_451) > 0.0D0) THEN
          factor_var_449 = (thermodynamics_var_441 % pressure_hl(jcol_var_450, jlev_var_451 + 1) - thermodynamics_var_441 % pressure_hl(jcol_var_450, jlev_var_451)) / (9.80665D0 * cloud_var_442 % fraction(jcol_var_450, jlev_var_451))
          lwp_in_cloud = factor_var_449 * cloud_var_442 % q_liq(jcol_var_450, jlev_var_451)
          iwp_in_cloud = factor_var_449 * cloud_var_442 % q_ice(jcol_var_450, jlev_var_451)
          IF (lwp_in_cloud > 0.0D0) THEN
            CALL calc_liq_optics_socrates(16, config_var_440 % cloud_optics % liq_coeff_lw, lwp_in_cloud, cloud_var_442 % re_liq(jcol_var_450, jlev_var_451), od_lw_liq, scat_od_lw_liq, g_lw_liq)
            CALL calc_liq_optics_socrates(14, config_var_440 % cloud_optics % liq_coeff_sw, lwp_in_cloud, cloud_var_442 % re_liq(jcol_var_450, jlev_var_451), od_sw_liq, scat_od_sw_liq, g_sw_liq)
            CALL delta_eddington_scat_od(od_sw_liq, scat_od_sw_liq, g_sw_liq)
          ELSE
            od_lw_liq = 0.0D0
            scat_od_lw_liq = 0.0D0
            g_lw_liq = 0.0D0
            od_sw_liq = 0.0D0
            scat_od_sw_liq = 0.0D0
            g_sw_liq = 0.0D0
          END IF
          IF (iwp_in_cloud > 0.0D0) THEN
            CALL calc_ice_optics_fu_lw(16, config_var_440 % cloud_optics % ice_coeff_lw, iwp_in_cloud, cloud_var_442 % re_ice(jcol_var_450, jlev_var_451), od_lw_ice, scat_od_lw_ice, g_lw_ice)
            CALL calc_ice_optics_fu_sw(14, config_var_440 % cloud_optics % ice_coeff_sw, iwp_in_cloud, cloud_var_442 % re_ice(jcol_var_450, jlev_var_451), od_sw_ice, scat_od_sw_ice, g_sw_ice)
            CALL calc_ice_optics_yi_lw(16, config_var_440 % cloud_optics % ice_coeff_lw, iwp_in_cloud, cloud_var_442 % re_ice(jcol_var_450, jlev_var_451), od_lw_ice, scat_od_lw_ice, g_lw_ice)
            CALL calc_ice_optics_yi_sw(14, config_var_440 % cloud_optics % ice_coeff_sw, iwp_in_cloud, cloud_var_442 % re_ice(jcol_var_450, jlev_var_451), od_sw_ice, scat_od_sw_ice, g_sw_ice)
            CALL delta_eddington_scat_od(od_sw_ice, scat_od_sw_ice, g_sw_ice)
            CALL delta_eddington_scat_od(od_lw_ice, scat_od_lw_ice, g_lw_ice)
          ELSE
            od_lw_ice = 0.0D0
            scat_od_lw_ice = 0.0D0
            g_lw_ice = 0.0D0
            od_sw_ice = 0.0D0
            scat_od_sw_ice = 0.0D0
            g_sw_ice = 0.0D0
          END IF
          DO jb_var_452 = 1, 16
            od_lw_cloud_var_443(jb_var_452, jlev_var_451, jcol_var_450) = od_lw_liq(jb_var_452) - scat_od_lw_liq(jb_var_452) + od_lw_ice(jb_var_452) - scat_od_lw_ice(jb_var_452)
          END DO
          DO jb_var_452 = 1, 14
            od_sw_cloud_var_446(jb_var_452, jlev_var_451, jcol_var_450) = od_sw_liq(jb_var_452) + od_sw_ice(jb_var_452)
            g_sw_cloud_var_448(jb_var_452, jlev_var_451, jcol_var_450) = (g_sw_liq(jb_var_452) * scat_od_sw_liq(jb_var_452) + g_sw_ice(jb_var_452) * scat_od_sw_ice(jb_var_452)) / (scat_od_sw_liq(jb_var_452) + scat_od_sw_ice(jb_var_452))
            ssa_sw_cloud_var_447(jb_var_452, jlev_var_451, jcol_var_450) = (scat_od_sw_liq(jb_var_452) + scat_od_sw_ice(jb_var_452)) / (od_sw_liq(jb_var_452) + od_sw_ice(jb_var_452))
          END DO
        END IF
      END DO
    END DO
  END SUBROUTINE cloud_optics_fn_438
END MODULE radiation_cloud_optics
MODULE radiation_ifs_rrtm
  USE radiation_config, ONLY: config_type
  USE radiation_single_level, ONLY: single_level_type
  USE radiation_thermodynamics, ONLY: thermodynamics_type
  USE radiation_gas, ONLY: assert_units_gas, gas_type
  USE yoerrtwn, ONLY: delwave, totplnk
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE gas_optics(ncol_var_453, nlev_var_454, istartcol_var_455, iendcol_var_456, config_var_457, single_level_var_458, thermodynamics_var_459, gas_var_460, od_lw_var_462, od_sw_var_463, ssa_sw_var_464, lw_albedo_var_461, planck_hl_var_465, lw_emission_var_466, incoming_sw_var_467)
    USE radiation_config, ONLY: config_type
    USE radiation_single_level, ONLY: single_level_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE radiation_gas, ONLY: assert_units_gas, gas_type
    INTEGER, INTENT(IN) :: ncol_var_453
    INTEGER, INTENT(IN) :: nlev_var_454
    INTEGER, INTENT(IN) :: istartcol_var_455, iendcol_var_456
    TYPE(config_type), INTENT(IN) :: config_var_457
    TYPE(single_level_type), INTENT(IN) :: single_level_var_458
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics_var_459
    TYPE(gas_type), INTENT(IN) :: gas_var_460
    REAL(KIND = 8), DIMENSION(140, istartcol_var_455 : iendcol_var_456), INTENT(IN), OPTIONAL :: lw_albedo_var_461
    REAL(KIND = 8), DIMENSION(140, nlev_var_454, istartcol_var_455 : iendcol_var_456), INTENT(OUT) :: od_lw_var_462
    REAL(KIND = 8), DIMENSION(112, nlev_var_454, istartcol_var_455 : iendcol_var_456), INTENT(OUT) :: od_sw_var_463, ssa_sw_var_464
    REAL(KIND = 8), DIMENSION(140, nlev_var_454 + 1, istartcol_var_455 : iendcol_var_456), INTENT(OUT), OPTIONAL :: planck_hl_var_465
    REAL(KIND = 8), DIMENSION(140, istartcol_var_455 : iendcol_var_456), INTENT(OUT), OPTIONAL :: lw_emission_var_466
    REAL(KIND = 8), DIMENSION(112, istartcol_var_455 : iendcol_var_456), INTENT(OUT), OPTIONAL :: incoming_sw_var_467
    REAL(KIND = 8) :: incoming_sw_scale(istartcol_var_455 : iendcol_var_456)
    REAL(KIND = 8) :: zod_lw(140, nlev_var_454, istartcol_var_455 : iendcol_var_456)
    REAL(KIND = 8) :: zod_sw(istartcol_var_455 : iendcol_var_456, nlev_var_454, 112)
    REAL(KIND = 8) :: zssa_sw(istartcol_var_455 : iendcol_var_456, nlev_var_454, 112)
    REAL(KIND = 8) :: zincsol(istartcol_var_455 : iendcol_var_456, 112)
    REAL(KIND = 8) :: zcolmol(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zcoldry(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zwbrodl(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zcolbrd(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zwkl(istartcol_var_455 : iendcol_var_456, 35, nlev_var_454)
    REAL(KIND = 8) :: zwx(istartcol_var_455 : iendcol_var_456, 4, nlev_var_454)
    REAL(KIND = 8) :: zfluxfac_var_468, zpi
    REAL(KIND = 8) :: ztauaerl(istartcol_var_455 : iendcol_var_456, nlev_var_454, 16)
    REAL(KIND = 8) :: zfac00(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zfac01(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zfac10(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zfac11(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zforfac(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zforfrac(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    INTEGER :: indfor_var_469(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    INTEGER :: indminor_var_470(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zscaleminor(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zscaleminorn2(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zminorfrac(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zrat_h2oco2(istartcol_var_455 : iendcol_var_456, nlev_var_454), zrat_h2oco2_1(istartcol_var_455 : iendcol_var_456, nlev_var_454), zrat_h2oo3(istartcol_var_455 : iendcol_var_456, nlev_var_454), zrat_h2oo3_1(istartcol_var_455 : iendcol_var_456, nlev_var_454), zrat_h2on2o(istartcol_var_455 : iendcol_var_456, nlev_var_454), zrat_h2on2o_1(istartcol_var_455 : iendcol_var_456, nlev_var_454), zrat_h2och4(istartcol_var_455 : iendcol_var_456, nlev_var_454), zrat_h2och4_1(istartcol_var_455 : iendcol_var_456, nlev_var_454), zrat_n2oco2(istartcol_var_455 : iendcol_var_456, nlev_var_454), zrat_n2oco2_1(istartcol_var_455 : iendcol_var_456, nlev_var_454), zrat_o3co2(istartcol_var_455 : iendcol_var_456, nlev_var_454), zrat_o3co2_1(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    INTEGER :: jp_var_471(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    INTEGER :: jt_var_472(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    INTEGER :: jt1_var_473(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zoneminus, zoneminus_array(istartcol_var_455 : iendcol_var_456)
    REAL(KIND = 8) :: zcolh2o(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zcolco2(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zcolo3(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zcoln2o(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zcolch4(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zcolo2(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zco2mult(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    INTEGER :: ilaytrop(istartcol_var_455 : iendcol_var_456)
    INTEGER :: ilayswtch(istartcol_var_455 : iendcol_var_456)
    INTEGER :: ilaylow(istartcol_var_455 : iendcol_var_456)
    REAL(KIND = 8) :: zpavel(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: ztavel(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zpz(istartcol_var_455 : iendcol_var_456, 0 : nlev_var_454)
    REAL(KIND = 8) :: ztz(istartcol_var_455 : iendcol_var_456, 0 : nlev_var_454)
    REAL(KIND = 8) :: zselffac(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zselffrac(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    INTEGER :: indself_var_474(istartcol_var_455 : iendcol_var_456, nlev_var_454)
    REAL(KIND = 8) :: zpfrac(istartcol_var_455 : iendcol_var_456, 140, nlev_var_454)
    INTEGER :: ireflect(istartcol_var_455 : iendcol_var_456)
    REAL(KIND = 8) :: pressure_fl_var_475(ncol_var_453, nlev_var_454), temperature_fl_var_476(ncol_var_453, nlev_var_454)
    INTEGER :: istartlev_var_477, iendlev_var_478
    LOGICAL :: do_sw, do_lw
    INTEGER :: jlev_var_479, jg_var_480, jcol_var_481
    do_sw = .TRUE.
    do_lw = .TRUE.
    iendlev_var_478 = UBOUND(gas_var_460 % mixing_ratio, 2)
    istartlev_var_477 = iendlev_var_478 - nlev_var_454 + 1
    zpi = 1.5707963267948966D0
    zfluxfac_var_468 = 15707.963267948966D0
    zoneminus = 0.999999D0
    DO jcol_var_481 = istartcol_var_455, iendcol_var_456
      zoneminus_array(jcol_var_481) = zoneminus
    END DO
    DO jlev_var_479 = 1, nlev_var_454
      DO jcol_var_481 = istartcol_var_455, iendcol_var_456
        pressure_fl_var_475(jcol_var_481, jlev_var_479) = thermodynamics_var_459 % pressure_fl(jcol_var_481, jlev_var_479)
        temperature_fl_var_476(jcol_var_481, jlev_var_479) = thermodynamics_var_459 % temperature_fl(jcol_var_481, jlev_var_479)
      END DO
    END DO
    CALL assert_units_gas(gas_var_460, 0)
    CALL rrtm_prepare_gases(istartcol_var_455, iendcol_var_456, ncol_var_453, nlev_var_454, thermodynamics_var_459 % pressure_hl(:, istartlev_var_477 : iendlev_var_478 + 1), pressure_fl_var_475, thermodynamics_var_459 % temperature_hl(:, istartlev_var_477 : iendlev_var_478 + 1), temperature_fl_var_476, gas_var_460 % mixing_ratio(:, istartlev_var_477 : iendlev_var_478, 1), gas_var_460 % mixing_ratio(:, istartlev_var_477 : iendlev_var_478, 2), gas_var_460 % mixing_ratio(:, istartlev_var_477 : iendlev_var_478, 6), gas_var_460 % mixing_ratio(:, istartlev_var_477 : iendlev_var_478, 4), gas_var_460 % mixing_ratio(:, istartlev_var_477 : iendlev_var_478, 12), gas_var_460 % mixing_ratio(:, istartlev_var_477 : iendlev_var_478, 8), gas_var_460 % mixing_ratio(:, istartlev_var_477 : iendlev_var_478, 9), gas_var_460 % mixing_ratio(:, istartlev_var_477 : iendlev_var_478, 10), gas_var_460 % mixing_ratio(:, istartlev_var_477 : iendlev_var_478, 11), gas_var_460 % mixing_ratio(:, istartlev_var_477 : iendlev_var_478, 3), zcoldry, zwbrodl, zwkl, zwx, zpavel, ztavel, zpz, ztz, ireflect)
    CALL rrtm_setcoef_140gp(istartcol_var_455, iendcol_var_456, nlev_var_454, zcoldry, zwbrodl, zwkl, zfac00, zfac01, zfac10, zfac11, zforfac, zforfrac, indfor_var_469, jp_var_471, jt_var_472, jt1_var_473, zcolh2o, zcolco2, zcolo3, zcoln2o, zcolch4, zcolo2, zco2mult, zcolbrd, ilaytrop, ilayswtch, ilaylow, zpavel, ztavel, zselffac, zselffrac, indself_var_474, indminor_var_470, zscaleminor, zscaleminorn2, zminorfrac, zrat_h2oco2, zrat_h2oco2_1, zrat_h2oo3, zrat_h2oo3_1, zrat_h2on2o, zrat_h2on2o_1, zrat_h2och4, zrat_h2och4_1, zrat_n2oco2, zrat_n2oco2_1, zrat_o3co2, zrat_o3co2_1)
    DO jg_var_480 = 1, 16
      DO jlev_var_479 = 1, nlev_var_454
        DO jcol_var_481 = istartcol_var_455, iendcol_var_456
          ztauaerl(jcol_var_481, jlev_var_479, jg_var_480) = 0.0D0
        END DO
      END DO
    END DO
    CALL rrtm_gas_optical_depth(istartcol_var_455, iendcol_var_456, nlev_var_454, zod_lw, zpavel, zcoldry, zcolbrd, zwx, ztauaerl, zfac00, zfac01, zfac10, zfac11, zforfac, zforfrac, indfor_var_469, jp_var_471, jt_var_472, jt1_var_473, 0.999999D0, zcolh2o, zcolco2, zcolo3, zcoln2o, zcolch4, zcolo2, zco2mult, ilaytrop, ilayswtch, ilaylow, zselffac, zselffrac, indself_var_474, zpfrac, indminor_var_470, zscaleminor, zscaleminorn2, zminorfrac, zrat_h2oco2, zrat_h2oco2_1, zrat_h2oo3, zrat_h2oo3_1, zrat_h2on2o, zrat_h2on2o_1, zrat_h2och4, zrat_h2och4_1, zrat_n2oco2, zrat_n2oco2_1, zrat_o3co2, zrat_o3co2_1)
    CALL planck_function_atmos(nlev_var_454, istartcol_var_455, iendcol_var_456, config_var_457, thermodynamics_var_459, zpfrac, planck_hl_var_465)
    CALL planck_function_surf(istartcol_var_455, iendcol_var_456, config_var_457, single_level_var_458 % skin_temperature, zpfrac(:, :, 1), lw_emission_var_466)
    DO jcol_var_481 = istartcol_var_455, iendcol_var_456
      DO jg_var_480 = 1, 140
        lw_emission_var_466(jg_var_480, jcol_var_481) = lw_emission_var_466(jg_var_480, jcol_var_481) * (1.0D0 - lw_albedo_var_461(jg_var_480, jcol_var_481))
      END DO
    END DO
    DO jcol_var_481 = istartcol_var_455, iendcol_var_456
      DO jlev_var_479 = 1, nlev_var_454
        DO jg_var_480 = 1, 140
          od_lw_var_462(jg_var_480, jlev_var_479, jcol_var_481) = MAX(1D-15, zod_lw(jg_var_480, nlev_var_454 + 1 - jlev_var_479, jcol_var_481))
        END DO
      END DO
    END DO
    CALL srtm_setcoef(istartcol_var_455, iendcol_var_456, nlev_var_454, zpavel, ztavel, zcoldry, zwkl, ilaytrop, zcolch4, zcolco2, zcolh2o, zcolmol, zcolo2, zcolo3, zforfac, zforfrac, indfor_var_469, zselffac, zselffrac, indself_var_474, zfac00, zfac01, zfac10, zfac11, jp_var_471, jt_var_472, jt1_var_473, single_level_var_458 % cos_sza(istartcol_var_455 : iendcol_var_456))
    DO jg_var_480 = 1, 112
      DO jlev_var_479 = 1, nlev_var_454
        DO jcol_var_481 = istartcol_var_455, iendcol_var_456
          zod_sw(jcol_var_481, jlev_var_479, jg_var_480) = 0.0D0
          zssa_sw(jcol_var_481, jlev_var_479, jg_var_480) = 0.0D0
        END DO
      END DO
    END DO
    DO jg_var_480 = 1, 112
      DO jcol_var_481 = istartcol_var_455, iendcol_var_456
        zincsol(jcol_var_481, jg_var_480) = 0.0D0
      END DO
    END DO
    CALL srtm_gas_optical_depth(istartcol_var_455, iendcol_var_456, nlev_var_454, zoneminus_array, single_level_var_458 % cos_sza(istartcol_var_455 : iendcol_var_456), ilaytrop, zcolch4, zcolco2, zcolh2o, zcolmol, zcolo2, zcolo3, zforfac, zforfrac, indfor_var_469, zselffac, zselffrac, indself_var_474, zfac00, zfac01, zfac10, zfac11, jp_var_471, jt_var_472, jt1_var_473, zod_sw, zssa_sw, zincsol)
    DO jg_var_480 = 1, 112
      DO jcol_var_481 = istartcol_var_455, iendcol_var_456
        zincsol(jcol_var_481, jg_var_480) = zincsol(jcol_var_481, jg_var_480) * single_level_var_458 % spectral_solar_scaling(config_var_457 % i_band_from_reordered_g_sw(jg_var_480))
      END DO
    END DO
    DO jcol_var_481 = istartcol_var_455, iendcol_var_456
      IF (single_level_var_458 % cos_sza(jcol_var_481) > 0.0D0) THEN
        incoming_sw_scale(jcol_var_481) = 1.0D0 / SUM(zincsol(jcol_var_481, :))
      ELSE
        incoming_sw_scale(jcol_var_481) = 1.0D0
      END IF
    END DO
    DO jcol_var_481 = istartcol_var_455, iendcol_var_456
      DO jlev_var_479 = 1, nlev_var_454
        DO jg_var_480 = 1, 112
          od_sw_var_463(jg_var_480, nlev_var_454 + 1 - jlev_var_479, jcol_var_481) = MAX(0.0D0, zod_sw(jcol_var_481, jlev_var_479, jg_var_480))
          ssa_sw_var_464(jg_var_480, nlev_var_454 + 1 - jlev_var_479, jcol_var_481) = zssa_sw(jcol_var_481, jlev_var_479, jg_var_480)
        END DO
      END DO
      DO jg_var_480 = 1, 112
        incoming_sw_var_467(jg_var_480, jcol_var_481) = incoming_sw_scale(jcol_var_481) * zincsol(jcol_var_481, jg_var_480)
      END DO
    END DO
  END SUBROUTINE gas_optics
  SUBROUTINE planck_function_atmos(nlev_var_482, istartcol_var_483, iendcol_var_484, config_var_485, thermodynamics_var_486, pfrac_var_487, planck_hl_var_488)
    USE radiation_config, ONLY: config_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE yoerrtwn, ONLY: delwave, totplnk
    INTEGER, INTENT(IN) :: nlev_var_482
    INTEGER, INTENT(IN) :: istartcol_var_483, iendcol_var_484
    TYPE(config_type), INTENT(IN) :: config_var_485
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics_var_486
    REAL(KIND = 8), INTENT(IN) :: pfrac_var_487(istartcol_var_483 : iendcol_var_484, 140, nlev_var_482)
    REAL(KIND = 8), DIMENSION(140, nlev_var_482 + 1, istartcol_var_483 : iendcol_var_484), INTENT(OUT) :: planck_hl_var_488
    REAL(KIND = 8), DIMENSION(istartcol_var_483 : iendcol_var_484, 16) :: planck_store_var_489
    REAL(KIND = 8), DIMENSION(istartcol_var_483 : iendcol_var_484) :: frac_var_490
    INTEGER, DIMENSION(istartcol_var_483 : iendcol_var_484) :: ind_var_491
    REAL(KIND = 8) :: temperature_var_492
    REAL(KIND = 8) :: factor_var_493, planck_tmp(istartcol_var_483 : iendcol_var_484, 140)
    REAL(KIND = 8) :: zfluxfac_var_494
    INTEGER :: jlev_var_495, jg_var_496, iband_var_497, jband_var_498, jcol_var_499, ilevoffset
    zfluxfac_var_494 = 15707.963267948966D0
    ilevoffset = UBOUND(thermodynamics_var_486 % temperature_hl, 2) - nlev_var_482 - 1
    DO jlev_var_495 = 1, nlev_var_482 + 1
      DO jcol_var_499 = istartcol_var_483, iendcol_var_484
        temperature_var_492 = thermodynamics_var_486 % temperature_hl(jcol_var_499, jlev_var_495 + ilevoffset)
        IF (temperature_var_492 < 339.0D0 .AND. temperature_var_492 >= 160.0D0) THEN
          ind_var_491(jcol_var_499) = INT(temperature_var_492 - 159.0D0)
          frac_var_490(jcol_var_499) = temperature_var_492 - INT(temperature_var_492)
        ELSE IF (temperature_var_492 >= 339.0D0) THEN
          ind_var_491(jcol_var_499) = 180
          frac_var_490(jcol_var_499) = temperature_var_492 - 339.0D0
        ELSE
          ind_var_491(jcol_var_499) = 1
          frac_var_490(jcol_var_499) = 0.0D0
        END IF
      END DO
      DO jband_var_498 = 1, 16
        factor_var_493 = zfluxfac_var_494 * delwave(jband_var_498)
        DO jcol_var_499 = istartcol_var_483, iendcol_var_484
          planck_store_var_489(jcol_var_499, jband_var_498) = factor_var_493 * (totplnk(ind_var_491(jcol_var_499), jband_var_498) + frac_var_490(jcol_var_499) * (totplnk(ind_var_491(jcol_var_499) + 1, jband_var_498) - totplnk(ind_var_491(jcol_var_499), jband_var_498)))
        END DO
      END DO
      IF (jlev_var_495 == 1) THEN
        DO jg_var_496 = 1, 140
          iband_var_497 = config_var_485 % i_band_from_g_lw(jg_var_496)
          DO jcol_var_499 = istartcol_var_483, iendcol_var_484
            planck_hl_var_488(jg_var_496, 1, jcol_var_499) = planck_store_var_489(jcol_var_499, iband_var_497) * pfrac_var_487(jcol_var_499, jg_var_496, nlev_var_482)
          END DO
        END DO
      ELSE
        DO jg_var_496 = 1, 140
          iband_var_497 = config_var_485 % i_band_from_g_lw(jg_var_496)
          DO jcol_var_499 = istartcol_var_483, iendcol_var_484
            planck_tmp(jcol_var_499, jg_var_496) = planck_store_var_489(jcol_var_499, iband_var_497) * pfrac_var_487(jcol_var_499, jg_var_496, nlev_var_482 + 2 - jlev_var_495)
          END DO
        END DO
        DO jcol_var_499 = istartcol_var_483, iendcol_var_484
          DO jg_var_496 = 1, 140
            planck_hl_var_488(jg_var_496, jlev_var_495, jcol_var_499) = planck_tmp(jcol_var_499, jg_var_496)
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE planck_function_atmos
  SUBROUTINE planck_function_surf(istartcol_var_500, iendcol_var_501, config_var_502, temperature_var_503, pfrac_var_504, planck_surf)
    USE radiation_config, ONLY: config_type
    USE yoerrtwn, ONLY: delwave, totplnk
    INTEGER, INTENT(IN) :: istartcol_var_500, iendcol_var_501
    TYPE(config_type), INTENT(IN) :: config_var_502
    REAL(KIND = 8), INTENT(IN) :: temperature_var_503(:)
    REAL(KIND = 8), INTENT(IN) :: pfrac_var_504(istartcol_var_500 : iendcol_var_501, 140)
    REAL(KIND = 8), DIMENSION(140, istartcol_var_500 : iendcol_var_501), INTENT(OUT) :: planck_surf
    REAL(KIND = 8), DIMENSION(istartcol_var_500 : iendcol_var_501, 16) :: planck_store_var_505
    REAL(KIND = 8), DIMENSION(istartcol_var_500 : iendcol_var_501) :: frac_var_506
    INTEGER, DIMENSION(istartcol_var_500 : iendcol_var_501) :: ind_var_507
    REAL(KIND = 8) :: tsurf
    REAL(KIND = 8) :: factor_var_508
    REAL(KIND = 8) :: zfluxfac_var_509
    INTEGER :: jg_var_510, iband_var_511, jband_var_512, jcol_var_513
    zfluxfac_var_509 = 15707.963267948966D0
    DO jcol_var_513 = istartcol_var_500, iendcol_var_501
      tsurf = temperature_var_503(jcol_var_513)
      IF (tsurf < 339.0D0 .AND. tsurf >= 160.0D0) THEN
        ind_var_507(jcol_var_513) = INT(tsurf - 159.0D0)
        frac_var_506(jcol_var_513) = tsurf - INT(tsurf)
      ELSE IF (tsurf >= 339.0D0) THEN
        ind_var_507(jcol_var_513) = 180
        frac_var_506(jcol_var_513) = tsurf - 339.0D0
      ELSE
        ind_var_507(jcol_var_513) = 1
        frac_var_506(jcol_var_513) = 0.0D0
      END IF
    END DO
    DO jband_var_512 = 1, 16
      factor_var_508 = zfluxfac_var_509 * delwave(jband_var_512)
      DO jcol_var_513 = istartcol_var_500, iendcol_var_501
        planck_store_var_505(jcol_var_513, jband_var_512) = factor_var_508 * (totplnk(ind_var_507(jcol_var_513), jband_var_512) + frac_var_506(jcol_var_513) * (totplnk(ind_var_507(jcol_var_513) + 1, jband_var_512) - totplnk(ind_var_507(jcol_var_513), jband_var_512)))
      END DO
    END DO
    DO jg_var_510 = 1, 140
      iband_var_511 = config_var_502 % i_band_from_g_lw(jg_var_510)
      DO jcol_var_513 = istartcol_var_500, iendcol_var_501
        planck_surf(jg_var_510, jcol_var_513) = planck_store_var_505(jcol_var_513, iband_var_511) * pfrac_var_504(jcol_var_513, jg_var_510)
      END DO
    END DO
  END SUBROUTINE planck_function_surf
END MODULE radiation_ifs_rrtm
MODULE radiation_mcica_lw
  CONTAINS
  SUBROUTINE solver_mcica_lw(nlev_var_514, istartcol_var_515, iendcol_var_516, config_var_517, single_level_var_518, cloud_var_519, od_var_520, ssa_var_521, g_var_522, od_cloud_var_523, ssa_cloud_var_524, g_cloud_var_525, planck_hl_var_526, emission, albedo_var_527, flux_var_528)
    USE radiation_config, ONLY: config_type
    USE radiation_single_level, ONLY: single_level_type
    USE radiation_cloud, ONLY: cloud_type
    USE radiation_flux, ONLY: flux_type
    USE radiation_two_stream, ONLY: calc_no_scattering_transmittance_lw
    USE radiation_adding_ica_lw, ONLY: calc_fluxes_no_scattering_lw
    USE radiation_cloud_generator, ONLY: cloud_generator
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlev_var_514
    INTEGER, INTENT(IN) :: istartcol_var_515, iendcol_var_516
    TYPE(config_type), INTENT(IN) :: config_var_517
    TYPE(single_level_type), INTENT(IN) :: single_level_var_518
    TYPE(cloud_type), INTENT(IN) :: cloud_var_519
    REAL(KIND = 8), INTENT(IN), DIMENSION(140, nlev_var_514, istartcol_var_515 : iendcol_var_516) :: od_var_520
    REAL(KIND = 8), INTENT(IN), DIMENSION(0, nlev_var_514, istartcol_var_515 : iendcol_var_516) :: ssa_var_521, g_var_522
    REAL(KIND = 8), INTENT(IN), DIMENSION(16, nlev_var_514, istartcol_var_515 : iendcol_var_516) :: od_cloud_var_523
    REAL(KIND = 8), INTENT(IN), DIMENSION(0, nlev_var_514, istartcol_var_515 : iendcol_var_516) :: ssa_cloud_var_524, g_cloud_var_525
    REAL(KIND = 8), INTENT(IN), DIMENSION(140, nlev_var_514 + 1, istartcol_var_515 : iendcol_var_516) :: planck_hl_var_526
    REAL(KIND = 8), INTENT(IN), DIMENSION(140, istartcol_var_515 : iendcol_var_516) :: emission, albedo_var_527
    TYPE(flux_type), INTENT(INOUT) :: flux_var_528
    REAL(KIND = 8), DIMENSION(140, nlev_var_514) :: ref_clear_var_529, trans_clear_var_530, reflectance_var_531, transmittance_var_532
    REAL(KIND = 8), DIMENSION(140, nlev_var_514) :: source_up_clear, source_dn_clear, source_up_var_533, source_dn_var_534
    REAL(KIND = 8), DIMENSION(140, nlev_var_514 + 1) :: flux_up_var_535, flux_dn_var_536
    REAL(KIND = 8), DIMENSION(140, nlev_var_514 + 1) :: flux_up_clear, flux_dn_clear
    REAL(KIND = 8), DIMENSION(140) :: od_total_var_537, ssa_total_var_538, g_total_var_539
    REAL(KIND = 8), DIMENSION(140, nlev_var_514) :: od_scaling_var_540
    REAL(KIND = 8), DIMENSION(140) :: od_cloud_new_var_541
    REAL(KIND = 8) :: total_cloud_cover_var_542
    LOGICAL :: is_clear_sky_layer(nlev_var_514)
    REAL(KIND = 8) :: sum_up_var_543, sum_dn
    INTEGER :: i_cloud_top
    INTEGER :: ng_var_544
    INTEGER :: jlev_var_545, jcol_var_546, jg_var_547
    ng_var_544 = 140
    DO jcol_var_546 = istartcol_var_515, iendcol_var_516
      CALL calc_no_scattering_transmittance_lw(ng_var_544 * nlev_var_514, od_var_520(:, :, jcol_var_546), planck_hl_var_526(:, 1 : nlev_var_514, jcol_var_546), planck_hl_var_526(:, 2 : nlev_var_514 + 1, jcol_var_546), trans_clear_var_530, source_up_clear, source_dn_clear)
      ref_clear_var_529 = 0.0D0
      CALL calc_fluxes_no_scattering_lw(ng_var_544, nlev_var_514, trans_clear_var_530, source_up_clear, source_dn_clear, emission(:, jcol_var_546), albedo_var_527(:, jcol_var_546), flux_up_clear, flux_dn_clear)
      DO jlev_var_545 = 1, nlev_var_514 + 1
        sum_up_var_543 = 0.0D0
        sum_dn = 0.0D0
        DO jg_var_547 = 1, ng_var_544
          sum_up_var_543 = sum_up_var_543 + flux_up_clear(jg_var_547, jlev_var_545)
          sum_dn = sum_dn + flux_dn_clear(jg_var_547, jlev_var_545)
        END DO
        flux_var_528 % lw_up_clear(jcol_var_546, jlev_var_545) = sum_up_var_543
        flux_var_528 % lw_dn_clear(jcol_var_546, jlev_var_545) = sum_dn
      END DO
      flux_var_528 % lw_dn_surf_clear_g(:, jcol_var_546) = flux_dn_clear(:, nlev_var_514 + 1)
      CALL cloud_generator(ng_var_544, nlev_var_514, 1, single_level_var_518 % iseed(jcol_var_546) + 997, 1D-06, cloud_var_519 % fraction(jcol_var_546, :), cloud_var_519 % overlap_param(jcol_var_546, :), 0.5D0, cloud_var_519 % fractional_std(jcol_var_546, :), config_var_517 % pdf_sampler, od_scaling_var_540, total_cloud_cover_var_542, use_beta_overlap = .FALSE., use_vectorizable_generator = .FALSE.)
      flux_var_528 % cloud_cover_lw(jcol_var_546) = total_cloud_cover_var_542
      IF (total_cloud_cover_var_542 >= 1D-06) THEN
        is_clear_sky_layer = .TRUE.
        i_cloud_top = nlev_var_514 + 1
        DO jlev_var_545 = 1, nlev_var_514
          IF (cloud_var_519 % fraction(jcol_var_546, jlev_var_545) >= 1D-06) THEN
            is_clear_sky_layer(jlev_var_545) = .FALSE.
            IF (i_cloud_top > jlev_var_545) THEN
              i_cloud_top = jlev_var_545
            END IF
            DO jg_var_547 = 1, ng_var_544
              od_cloud_new_var_541(jg_var_547) = od_scaling_var_540(jg_var_547, jlev_var_545) * od_cloud_var_523(config_var_517 % i_band_from_reordered_g_lw(jg_var_547), jlev_var_545, jcol_var_546)
              od_total_var_537(jg_var_547) = od_var_520(jg_var_547, jlev_var_545, jcol_var_546) + od_cloud_new_var_541(jg_var_547)
              ssa_total_var_538(jg_var_547) = 0.0D0
              g_total_var_539(jg_var_547) = 0.0D0
            END DO
            CALL calc_no_scattering_transmittance_lw(ng_var_544, od_total_var_537, planck_hl_var_526(:, jlev_var_545, jcol_var_546), planck_hl_var_526(:, jlev_var_545 + 1, jcol_var_546), transmittance_var_532(:, jlev_var_545), source_up_var_533(:, jlev_var_545), source_dn_var_534(:, jlev_var_545))
          ELSE
            DO jg_var_547 = 1, ng_var_544
              reflectance_var_531(jg_var_547, jlev_var_545) = ref_clear_var_529(jg_var_547, jlev_var_545)
              transmittance_var_532(jg_var_547, jlev_var_545) = trans_clear_var_530(jg_var_547, jlev_var_545)
              source_up_var_533(jg_var_547, jlev_var_545) = source_up_clear(jg_var_547, jlev_var_545)
              source_dn_var_534(jg_var_547, jlev_var_545) = source_dn_clear(jg_var_547, jlev_var_545)
            END DO
          END IF
        END DO
        CALL calc_fluxes_no_scattering_lw(ng_var_544, nlev_var_514, transmittance_var_532, source_up_var_533, source_dn_var_534, emission(:, jcol_var_546), albedo_var_527(:, jcol_var_546), flux_up_var_535, flux_dn_var_536)
        DO jlev_var_545 = 1, nlev_var_514 + 1
          sum_up_var_543 = 0.0D0
          sum_dn = 0.0D0
          DO jg_var_547 = 1, ng_var_544
            sum_up_var_543 = sum_up_var_543 + flux_up_var_535(jg_var_547, jlev_var_545)
            sum_dn = sum_dn + flux_dn_var_536(jg_var_547, jlev_var_545)
          END DO
          flux_var_528 % lw_up(jcol_var_546, jlev_var_545) = sum_up_var_543
          flux_var_528 % lw_dn(jcol_var_546, jlev_var_545) = sum_dn
        END DO
        DO jlev_var_545 = 1, nlev_var_514 + 1
          flux_var_528 % lw_up(jcol_var_546, jlev_var_545) = total_cloud_cover_var_542 * flux_var_528 % lw_up(jcol_var_546, jlev_var_545) + (1.0D0 - total_cloud_cover_var_542) * flux_var_528 % lw_up_clear(jcol_var_546, jlev_var_545)
          flux_var_528 % lw_dn(jcol_var_546, jlev_var_545) = total_cloud_cover_var_542 * flux_var_528 % lw_dn(jcol_var_546, jlev_var_545) + (1.0D0 - total_cloud_cover_var_542) * flux_var_528 % lw_dn_clear(jcol_var_546, jlev_var_545)
        END DO
        flux_var_528 % lw_dn_surf_g(:, jcol_var_546) = total_cloud_cover_var_542 * flux_dn_var_536(:, nlev_var_514 + 1) + (1.0D0 - total_cloud_cover_var_542) * flux_var_528 % lw_dn_surf_clear_g(:, jcol_var_546)
      ELSE
        DO jlev_var_545 = 1, nlev_var_514 + 1
          flux_var_528 % lw_up(jcol_var_546, jlev_var_545) = flux_var_528 % lw_up_clear(jcol_var_546, jlev_var_545)
          flux_var_528 % lw_dn(jcol_var_546, jlev_var_545) = flux_var_528 % lw_dn_clear(jcol_var_546, jlev_var_545)
        END DO
        flux_var_528 % lw_dn_surf_g(:, jcol_var_546) = flux_var_528 % lw_dn_surf_clear_g(:, jcol_var_546)
      END IF
    END DO
  END SUBROUTINE solver_mcica_lw
END MODULE radiation_mcica_lw
MODULE radiation_mcica_sw
  CONTAINS
  SUBROUTINE solver_mcica_sw(nlev_var_548, istartcol_var_549, iendcol_var_550, config_var_551, single_level_var_552, cloud_var_553, od_var_554, ssa_var_555, g_var_556, od_cloud_var_557, ssa_cloud_var_558, g_cloud_var_559, albedo_direct, albedo_diffuse, incoming_sw_var_560, flux_var_561)
    USE radiation_config, ONLY: config_type
    USE radiation_single_level, ONLY: single_level_type
    USE radiation_cloud, ONLY: cloud_type
    USE radiation_flux, ONLY: flux_type
    USE radiation_two_stream, ONLY: calc_ref_trans_sw
    USE radiation_adding_ica_sw, ONLY: adding_ica_sw
    USE radiation_cloud_generator, ONLY: cloud_generator
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlev_var_548
    INTEGER, INTENT(IN) :: istartcol_var_549, iendcol_var_550
    TYPE(config_type), INTENT(IN) :: config_var_551
    TYPE(single_level_type), INTENT(IN) :: single_level_var_552
    TYPE(cloud_type), INTENT(IN) :: cloud_var_553
    REAL(KIND = 8), INTENT(IN), DIMENSION(112, nlev_var_548, istartcol_var_549 : iendcol_var_550) :: od_var_554, ssa_var_555, g_var_556
    REAL(KIND = 8), INTENT(IN), DIMENSION(14, nlev_var_548, istartcol_var_549 : iendcol_var_550) :: od_cloud_var_557, ssa_cloud_var_558, g_cloud_var_559
    REAL(KIND = 8), INTENT(IN), DIMENSION(112, istartcol_var_549 : iendcol_var_550) :: albedo_direct, albedo_diffuse, incoming_sw_var_560
    TYPE(flux_type), INTENT(INOUT) :: flux_var_561
    REAL(KIND = 8) :: cos_sza_var_562
    REAL(KIND = 8), DIMENSION(112, nlev_var_548) :: ref_clear_var_563, trans_clear_var_564, reflectance_var_565, transmittance_var_566
    REAL(KIND = 8), DIMENSION(112, nlev_var_548) :: ref_dir_clear, trans_dir_diff_clear, ref_dir_var_567, trans_dir_diff_var_568
    REAL(KIND = 8), DIMENSION(112, nlev_var_548) :: trans_dir_dir_clear, trans_dir_dir_var_569
    REAL(KIND = 8), DIMENSION(112, nlev_var_548 + 1) :: flux_up_var_570, flux_dn_diffuse_var_571, flux_dn_direct_var_572
    REAL(KIND = 8), DIMENSION(112) :: od_total_var_573, ssa_total_var_574, g_total_var_575
    REAL(KIND = 8) :: scat_od_var_576
    REAL(KIND = 8), DIMENSION(112, nlev_var_548) :: od_scaling_var_577
    REAL(KIND = 8), DIMENSION(112) :: od_cloud_new_var_578
    REAL(KIND = 8), DIMENSION(112, nlev_var_548 + 1) :: tmp_work_albedo, tmp_work_source
    REAL(KIND = 8), DIMENSION(112, nlev_var_548) :: tmp_work_inv_denominator
    REAL(KIND = 8) :: total_cloud_cover_var_579
    REAL(KIND = 8) :: sum_up_var_580, sum_dn_diff, sum_dn_dir
    INTEGER :: ng_var_581
    INTEGER :: jlev_var_582, jcol_var_583, jg_var_584
    ng_var_581 = 112
    DO jcol_var_583 = istartcol_var_549, iendcol_var_550
      IF (single_level_var_552 % cos_sza(jcol_var_583) > 0.0D0) THEN
        cos_sza_var_562 = single_level_var_552 % cos_sza(jcol_var_583)
        CALL calc_ref_trans_sw(ng_var_581 * nlev_var_548, cos_sza_var_562, od_var_554(:, :, jcol_var_583), ssa_var_555(:, :, jcol_var_583), g_var_556(:, :, jcol_var_583), ref_clear_var_563, trans_clear_var_564, ref_dir_clear, trans_dir_diff_clear, trans_dir_dir_clear)
        CALL adding_ica_sw(ng_var_581, nlev_var_548, incoming_sw_var_560(:, jcol_var_583), albedo_diffuse(:, jcol_var_583), albedo_direct(:, jcol_var_583), cos_sza_var_562, ref_clear_var_563, trans_clear_var_564, ref_dir_clear, trans_dir_diff_clear, trans_dir_dir_clear, flux_up_var_570, flux_dn_diffuse_var_571, flux_dn_direct_var_572, albedo_var_20 = tmp_work_albedo, source_var_21 = tmp_work_source, inv_denominator = tmp_work_inv_denominator)
        DO jlev_var_582 = 1, nlev_var_548 + 1
          sum_up_var_580 = 0.0D0
          sum_dn_diff = 0.0D0
          sum_dn_dir = 0.0D0
          DO jg_var_584 = 1, ng_var_581
            sum_up_var_580 = sum_up_var_580 + flux_up_var_570(jg_var_584, jlev_var_582)
            sum_dn_diff = sum_dn_diff + flux_dn_diffuse_var_571(jg_var_584, jlev_var_582)
            sum_dn_dir = sum_dn_dir + flux_dn_direct_var_572(jg_var_584, jlev_var_582)
          END DO
          flux_var_561 % sw_up_clear(jcol_var_583, jlev_var_582) = sum_up_var_580
          flux_var_561 % sw_dn_clear(jcol_var_583, jlev_var_582) = sum_dn_diff + sum_dn_dir
          flux_var_561 % sw_dn_direct_clear(jcol_var_583, jlev_var_582) = sum_dn_dir
        END DO
        DO jg_var_584 = 1, ng_var_581
          flux_var_561 % sw_dn_diffuse_surf_clear_g(jg_var_584, jcol_var_583) = flux_dn_diffuse_var_571(jg_var_584, nlev_var_548 + 1)
          flux_var_561 % sw_dn_direct_surf_clear_g(jg_var_584, jcol_var_583) = flux_dn_direct_var_572(jg_var_584, nlev_var_548 + 1)
        END DO
        CALL cloud_generator(ng_var_581, nlev_var_548, 1, single_level_var_552 % iseed(jcol_var_583), 1D-06, cloud_var_553 % fraction(jcol_var_583, :), cloud_var_553 % overlap_param(jcol_var_583, :), 0.5D0, cloud_var_553 % fractional_std(jcol_var_583, :), config_var_551 % pdf_sampler, od_scaling_var_577, total_cloud_cover_var_579, use_beta_overlap = .FALSE., use_vectorizable_generator = .FALSE.)
        flux_var_561 % cloud_cover_sw(jcol_var_583) = total_cloud_cover_var_579
        IF (total_cloud_cover_var_579 >= 1D-06) THEN
          DO jlev_var_582 = 1, nlev_var_548
            IF (cloud_var_553 % fraction(jcol_var_583, jlev_var_582) >= 1D-06) THEN
              DO jg_var_584 = 1, ng_var_581
                od_cloud_new_var_578(jg_var_584) = od_scaling_var_577(jg_var_584, jlev_var_582) * od_cloud_var_557(config_var_551 % i_band_from_reordered_g_sw(jg_var_584), jlev_var_582, jcol_var_583)
                od_total_var_573(jg_var_584) = od_var_554(jg_var_584, jlev_var_582, jcol_var_583) + od_cloud_new_var_578(jg_var_584)
                ssa_total_var_574(jg_var_584) = 0.0D0
                g_total_var_575(jg_var_584) = 0.0D0
                IF (od_total_var_573(jg_var_584) > 0.0D0) THEN
                  scat_od_var_576 = ssa_var_555(jg_var_584, jlev_var_582, jcol_var_583) * od_var_554(jg_var_584, jlev_var_582, jcol_var_583) + ssa_cloud_var_558(config_var_551 % i_band_from_reordered_g_sw(jg_var_584), jlev_var_582, jcol_var_583) * od_cloud_new_var_578(jg_var_584)
                  ssa_total_var_574(jg_var_584) = scat_od_var_576 / od_total_var_573(jg_var_584)
                  IF (scat_od_var_576 > 0.0D0) THEN
                    g_total_var_575(jg_var_584) = (g_var_556(jg_var_584, jlev_var_582, jcol_var_583) * ssa_var_555(jg_var_584, jlev_var_582, jcol_var_583) * od_var_554(jg_var_584, jlev_var_582, jcol_var_583) + g_cloud_var_559(config_var_551 % i_band_from_reordered_g_sw(jg_var_584), jlev_var_582, jcol_var_583) * ssa_cloud_var_558(config_var_551 % i_band_from_reordered_g_sw(jg_var_584), jlev_var_582, jcol_var_583) * od_cloud_new_var_578(jg_var_584)) / scat_od_var_576
                  END IF
                END IF
              END DO
              CALL calc_ref_trans_sw(ng_var_581, cos_sza_var_562, od_total_var_573, ssa_total_var_574, g_total_var_575, reflectance_var_565(:, jlev_var_582), transmittance_var_566(:, jlev_var_582), ref_dir_var_567(:, jlev_var_582), trans_dir_diff_var_568(:, jlev_var_582), trans_dir_dir_var_569(:, jlev_var_582))
            ELSE
              DO jg_var_584 = 1, ng_var_581
                reflectance_var_565(jg_var_584, jlev_var_582) = ref_clear_var_563(jg_var_584, jlev_var_582)
                transmittance_var_566(jg_var_584, jlev_var_582) = trans_clear_var_564(jg_var_584, jlev_var_582)
                ref_dir_var_567(jg_var_584, jlev_var_582) = ref_dir_clear(jg_var_584, jlev_var_582)
                trans_dir_diff_var_568(jg_var_584, jlev_var_582) = trans_dir_diff_clear(jg_var_584, jlev_var_582)
                trans_dir_dir_var_569(jg_var_584, jlev_var_582) = trans_dir_dir_clear(jg_var_584, jlev_var_582)
              END DO
            END IF
          END DO
          CALL adding_ica_sw(ng_var_581, nlev_var_548, incoming_sw_var_560(:, jcol_var_583), albedo_diffuse(:, jcol_var_583), albedo_direct(:, jcol_var_583), cos_sza_var_562, reflectance_var_565, transmittance_var_566, ref_dir_var_567, trans_dir_diff_var_568, trans_dir_dir_var_569, flux_up_var_570, flux_dn_diffuse_var_571, flux_dn_direct_var_572, albedo_var_20 = tmp_work_albedo, source_var_21 = tmp_work_source, inv_denominator = tmp_work_inv_denominator)
          DO jlev_var_582 = 1, nlev_var_548 + 1
            sum_up_var_580 = 0.0D0
            sum_dn_diff = 0.0D0
            sum_dn_dir = 0.0D0
            DO jg_var_584 = 1, ng_var_581
              sum_up_var_580 = sum_up_var_580 + flux_up_var_570(jg_var_584, jlev_var_582)
              sum_dn_diff = sum_dn_diff + flux_dn_diffuse_var_571(jg_var_584, jlev_var_582)
              sum_dn_dir = sum_dn_dir + flux_dn_direct_var_572(jg_var_584, jlev_var_582)
            END DO
            flux_var_561 % sw_up(jcol_var_583, jlev_var_582) = sum_up_var_580
            flux_var_561 % sw_dn(jcol_var_583, jlev_var_582) = sum_dn_diff + sum_dn_dir
            flux_var_561 % sw_dn_direct(jcol_var_583, jlev_var_582) = sum_dn_dir
          END DO
          DO jlev_var_582 = 1, nlev_var_548 + 1
            flux_var_561 % sw_up(jcol_var_583, jlev_var_582) = total_cloud_cover_var_579 * flux_var_561 % sw_up(jcol_var_583, jlev_var_582) + (1.0D0 - total_cloud_cover_var_579) * flux_var_561 % sw_up_clear(jcol_var_583, jlev_var_582)
            flux_var_561 % sw_dn(jcol_var_583, jlev_var_582) = total_cloud_cover_var_579 * flux_var_561 % sw_dn(jcol_var_583, jlev_var_582) + (1.0D0 - total_cloud_cover_var_579) * flux_var_561 % sw_dn_clear(jcol_var_583, jlev_var_582)
            flux_var_561 % sw_dn_direct(jcol_var_583, jlev_var_582) = total_cloud_cover_var_579 * flux_var_561 % sw_dn_direct(jcol_var_583, jlev_var_582) + (1.0D0 - total_cloud_cover_var_579) * flux_var_561 % sw_dn_direct_clear(jcol_var_583, jlev_var_582)
          END DO
          DO jg_var_584 = 1, ng_var_581
            flux_var_561 % sw_dn_diffuse_surf_g(jg_var_584, jcol_var_583) = flux_dn_diffuse_var_571(jg_var_584, nlev_var_548 + 1)
            flux_var_561 % sw_dn_direct_surf_g(jg_var_584, jcol_var_583) = flux_dn_direct_var_572(jg_var_584, nlev_var_548 + 1)
            flux_var_561 % sw_dn_diffuse_surf_g(jg_var_584, jcol_var_583) = total_cloud_cover_var_579 * flux_var_561 % sw_dn_diffuse_surf_g(jg_var_584, jcol_var_583) + (1.0D0 - total_cloud_cover_var_579) * flux_var_561 % sw_dn_diffuse_surf_clear_g(jg_var_584, jcol_var_583)
            flux_var_561 % sw_dn_direct_surf_g(jg_var_584, jcol_var_583) = total_cloud_cover_var_579 * flux_var_561 % sw_dn_direct_surf_g(jg_var_584, jcol_var_583) + (1.0D0 - total_cloud_cover_var_579) * flux_var_561 % sw_dn_direct_surf_clear_g(jg_var_584, jcol_var_583)
          END DO
        ELSE
          DO jlev_var_582 = 1, nlev_var_548 + 1
            flux_var_561 % sw_up(jcol_var_583, jlev_var_582) = flux_var_561 % sw_up_clear(jcol_var_583, jlev_var_582)
            flux_var_561 % sw_dn(jcol_var_583, jlev_var_582) = flux_var_561 % sw_dn_clear(jcol_var_583, jlev_var_582)
            flux_var_561 % sw_dn_direct(jcol_var_583, jlev_var_582) = flux_var_561 % sw_dn_direct_clear(jcol_var_583, jlev_var_582)
          END DO
          DO jg_var_584 = 1, ng_var_581
            flux_var_561 % sw_dn_diffuse_surf_g(jg_var_584, jcol_var_583) = flux_var_561 % sw_dn_diffuse_surf_clear_g(jg_var_584, jcol_var_583)
            flux_var_561 % sw_dn_direct_surf_g(jg_var_584, jcol_var_583) = flux_var_561 % sw_dn_direct_surf_clear_g(jg_var_584, jcol_var_583)
          END DO
        END IF
      ELSE
        DO jlev_var_582 = 1, nlev_var_548 + 1
          flux_var_561 % sw_up(jcol_var_583, jlev_var_582) = 0.0D0
          flux_var_561 % sw_dn(jcol_var_583, jlev_var_582) = 0.0D0
          flux_var_561 % sw_dn_direct(jcol_var_583, jlev_var_582) = 0.0D0
          flux_var_561 % sw_up_clear(jcol_var_583, jlev_var_582) = 0.0D0
          flux_var_561 % sw_dn_clear(jcol_var_583, jlev_var_582) = 0.0D0
          flux_var_561 % sw_dn_direct_clear(jcol_var_583, jlev_var_582) = 0.0D0
        END DO
        DO jg_var_584 = 1, ng_var_581
          flux_var_561 % sw_dn_diffuse_surf_g(jg_var_584, jcol_var_583) = 0.0D0
          flux_var_561 % sw_dn_direct_surf_g(jg_var_584, jcol_var_583) = 0.0D0
          flux_var_561 % sw_dn_diffuse_surf_clear_g(jg_var_584, jcol_var_583) = 0.0D0
          flux_var_561 % sw_dn_direct_surf_clear_g(jg_var_584, jcol_var_583) = 0.0D0
        END DO
      END IF
    END DO
  END SUBROUTINE solver_mcica_sw
END MODULE radiation_mcica_sw
MODULE radiation_interface
  USE radiation_config, ONLY: config_type
  USE radiation_single_level, ONLY: get_albedos, single_level_type
  USE radiation_thermodynamics, ONLY: thermodynamics_type
  USE radiation_gas, ONLY: gas_type
  USE radiation_cloud, ONLY: cloud_type, crop_cloud_fraction
  USE radiation_aerosol, ONLY: aerosol_type
  USE radiation_flux, ONLY: calc_surface_spectral, calc_toa_spectral, flux_type
  USE radiation_ifs_rrtm, ONLY: gas_optics
  USE radiation_cloud_optics, ONLY: cloud_optics_fn_438
  USE radiation_aerosol_optics, ONLY: add_aerosol_optics
  USE radiation_mcica_lw, ONLY: solver_mcica_lw
  USE radiation_mcica_sw, ONLY: solver_mcica_sw
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE radiation(ncol, nlev, istartcol, iendcol, config, single_level, thermodynamics, gas, cloud, aerosol, flux)
    USE radiation_config, ONLY: config_type
    USE radiation_single_level, ONLY: get_albedos, single_level_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE radiation_gas, ONLY: gas_type
    USE radiation_cloud, ONLY: cloud_type, crop_cloud_fraction
    USE radiation_aerosol, ONLY: aerosol_type
    USE radiation_flux, ONLY: calc_surface_spectral, calc_toa_spectral, flux_type
    USE radiation_ifs_rrtm, ONLY: gas_optics
    USE radiation_cloud_optics, ONLY: cloud_optics_fn_438
    USE radiation_aerosol_optics, ONLY: add_aerosol_optics
    USE radiation_mcica_lw, ONLY: solver_mcica_lw
    USE radiation_mcica_sw, ONLY: solver_mcica_sw
    INTEGER, INTENT(IN) :: ncol
    INTEGER, INTENT(IN) :: nlev
    INTEGER, INTENT(IN) :: istartcol, iendcol
    TYPE(config_type), INTENT(IN) :: config
    TYPE(single_level_type), INTENT(IN) :: single_level
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics
    TYPE(gas_type), INTENT(IN) :: gas
    TYPE(cloud_type), INTENT(INOUT) :: cloud
    TYPE(aerosol_type), INTENT(IN) :: aerosol
    TYPE(flux_type), INTENT(INOUT) :: flux
    REAL(KIND = 8), DIMENSION(140, nlev, istartcol : iendcol) :: od_lw_var_585
    REAL(KIND = 8), DIMENSION(0, nlev, istartcol : iendcol) :: ssa_lw_var_586, g_lw_var_587
    REAL(KIND = 8), DIMENSION(16, nlev, istartcol : iendcol) :: od_lw_cloud_var_588
    REAL(KIND = 8), DIMENSION(0, nlev, istartcol : iendcol) :: ssa_lw_cloud_var_589, g_lw_cloud_var_590
    REAL(KIND = 8), DIMENSION(112, nlev, istartcol : iendcol) :: od_sw_var_591, ssa_sw_var_592, g_sw_var_593
    REAL(KIND = 8), DIMENSION(14, nlev, istartcol : iendcol) :: od_sw_cloud_var_594, ssa_sw_cloud_var_595, g_sw_cloud_var_596
    REAL(KIND = 8), DIMENSION(140, nlev + 1, istartcol : iendcol) :: planck_hl_var_597
    REAL(KIND = 8), DIMENSION(140, istartcol : iendcol) :: lw_emission_var_598
    REAL(KIND = 8), DIMENSION(140, istartcol : iendcol) :: lw_albedo_var_599
    REAL(KIND = 8), DIMENSION(112, istartcol : iendcol) :: sw_albedo_direct_var_600
    REAL(KIND = 8), DIMENSION(112, istartcol : iendcol) :: sw_albedo_diffuse_var_601
    REAL(KIND = 8), DIMENSION(112, istartcol : iendcol) :: incoming_sw_var_602
    CALL global_init_fn
    CALL get_albedos(single_level, istartcol, iendcol, config, sw_albedo_direct_var_600, sw_albedo_diffuse_var_601, lw_albedo_var_599)
    !CALL gas_optics(ncol, nlev, istartcol, iendcol, config, single_level, thermodynamics, gas, od_lw_var_585, od_sw_var_591, ssa_sw_var_592, lw_albedo_var_461 = lw_albedo_var_599, planck_hl_var_465 = planck_hl_var_597, lw_emission_var_466 = lw_emission_var_598, incoming_sw_var_467 = incoming_sw_var_602)
    !CALL crop_cloud_fraction(cloud, istartcol, iendcol, 1D-06, 1D-09)
    !CALL cloud_optics_fn_438(nlev, istartcol, iendcol, config, thermodynamics, cloud, od_lw_cloud_var_588, ssa_lw_cloud_var_589, g_lw_cloud_var_590, od_sw_cloud_var_594, ssa_sw_cloud_var_595, g_sw_cloud_var_596)
    !CALL add_aerosol_optics(nlev, istartcol, iendcol, config, thermodynamics, gas, aerosol, od_lw_var_585, ssa_lw_var_586, g_lw_var_587, od_sw_var_591, ssa_sw_var_592, g_sw_var_593)
    !CALL solver_mcica_lw(nlev, istartcol, iendcol, config, single_level, cloud, od_lw_var_585, ssa_lw_var_586, g_lw_var_587, od_lw_cloud_var_588, ssa_lw_cloud_var_589, g_lw_cloud_var_590, planck_hl_var_597, lw_emission_var_598, lw_albedo_var_599, flux)
    !CALL solver_mcica_sw(nlev, istartcol, iendcol, config, single_level, cloud, od_sw_var_591, ssa_sw_var_592, g_sw_var_593, od_sw_cloud_var_594, ssa_sw_cloud_var_595, g_sw_cloud_var_596, sw_albedo_direct_var_600, sw_albedo_diffuse_var_601, incoming_sw_var_602, flux)
    !CALL calc_surface_spectral(flux, config, istartcol, iendcol)
    !CALL calc_toa_spectral(flux, config, istartcol, iendcol)
  END SUBROUTINE radiation
END MODULE radiation_interface
SUBROUTINE rrtm_taumol16(kidia_var_603, kfdia_var_604, klev_var_605, taug_var_606, p_tauaerl_var_607, fac00_var_608, fac01_var_609, fac10_var_610, fac11_var_611, forfac_var_626, forfrac_var_627, indfor_var_625, jp_var_612, jt_var_613, jt1_var_614, oneminus_var_615, colh2o_var_616, colch4_var_617, laytrop_var_618, selffac_var_619, selffrac_var_620, indself_var_621, fracs_var_622, rat_h2och4_var_623, rat_h2och4_1_var_624)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng16
  USE yoerrta16, ONLY: absa_var_149, absb_var_150, forref_var_152, fracrefa_var_147, fracrefb_var_148, selfref_var_151
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_603
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_604
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_605
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_606(kidia_var_603 : kfdia_var_604, 140, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_607(kidia_var_603 : kfdia_var_604, klev_var_605, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_608(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_609(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_610(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_611(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_612(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_613(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_614(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_615
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_616(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: colch4_var_617(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_618(kidia_var_603 : kfdia_var_604)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_619(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_620(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_621(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_622(kidia_var_603 : kfdia_var_604, 140, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_var_623(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_1_var_624(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_625(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_626(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_627(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4) :: ig_var_628, ind0_var_629, ind1_var_630, inds_var_631, indf_var_632, js_var_633, js1_var_634, jpl_var_635, lay_var_636
  REAL(KIND = 8) :: fac000_var_637, fac100_var_638, fac200_var_639, fac010_var_640, fac110_var_641, fac210_var_642, fac001_var_643, fac101_var_644, fac201_var_645, fac011_var_646, fac111_var_647, fac211_var_648
  REAL(KIND = 8) :: p_var_649, p4_var_650, fk0_var_651, fk1_var_652, fk2_var_653
  REAL(KIND = 8) :: refrat_planck_a_var_654
  REAL(KIND = 8) :: taufor_var_655, tauself_var_656, tau_major_var_657(2), tau_major1_var_658(2)
  REAL(KIND = 8) :: fs_var_659, specmult_var_660, specparm_var_661, speccomb_var_662, fs1_var_663, specmult1_var_664, specparm1_var_665, speccomb1_var_666, fpl_var_667, specmult_planck_var_668, specparm_planck_var_669, speccomb_planck_var_670
  INTEGER(KIND = 4) :: laytrop_min_var_671, laytrop_max_var_672
  INTEGER(KIND = 4) :: ixc_var_673(klev_var_605), ixlow_var_674(kfdia_var_604, klev_var_605), ixhigh_var_675(kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4) :: ich_var_676, icl_var_677, ixc0_var_678, ixp_var_679, jc_var_680, jl_var_681
  laytrop_min_var_671 = MINVAL(laytrop_var_618)
  laytrop_max_var_672 = MAXVAL(laytrop_var_618)
  ixlow_var_674 = 0
  ixhigh_var_675 = 0
  ixc_var_673 = 0
  DO lay_var_636 = laytrop_min_var_671 + 1, laytrop_max_var_672
    icl_var_677 = 0
    ich_var_676 = 0
    DO jc_var_680 = kidia_var_603, kfdia_var_604
      IF (lay_var_636 <= laytrop_var_618(jc_var_680)) THEN
        icl_var_677 = icl_var_677 + 1
        ixlow_var_674(icl_var_677, lay_var_636) = jc_var_680
      ELSE
        ich_var_676 = ich_var_676 + 1
        ixhigh_var_675(ich_var_676, lay_var_636) = jc_var_680
      END IF
    END DO
    ixc_var_673(lay_var_636) = icl_var_677
  END DO
  refrat_planck_a_var_654 = chi_mls(1, 6) / chi_mls(6, 6)
  DO lay_var_636 = 1, laytrop_min_var_671
    DO jl_var_681 = kidia_var_603, kfdia_var_604
      speccomb_var_662 = colh2o_var_616(jl_var_681, lay_var_636) + rat_h2och4_var_623(jl_var_681, lay_var_636) * colch4_var_617(jl_var_681, lay_var_636)
      specparm_var_661 = MIN(colh2o_var_616(jl_var_681, lay_var_636) / speccomb_var_662, oneminus_var_615)
      specmult_var_660 = 8.0D0 * (specparm_var_661)
      js_var_633 = 1 + INT(specmult_var_660)
      fs_var_659 = ((specmult_var_660) - AINT((specmult_var_660)))
      speccomb1_var_666 = colh2o_var_616(jl_var_681, lay_var_636) + rat_h2och4_1_var_624(jl_var_681, lay_var_636) * colch4_var_617(jl_var_681, lay_var_636)
      specparm1_var_665 = MIN(colh2o_var_616(jl_var_681, lay_var_636) / speccomb1_var_666, oneminus_var_615)
      specmult1_var_664 = 8.0D0 * (specparm1_var_665)
      js1_var_634 = 1 + INT(specmult1_var_664)
      fs1_var_663 = ((specmult1_var_664) - AINT((specmult1_var_664)))
      speccomb_planck_var_670 = colh2o_var_616(jl_var_681, lay_var_636) + refrat_planck_a_var_654 * colch4_var_617(jl_var_681, lay_var_636)
      specparm_planck_var_669 = MIN(colh2o_var_616(jl_var_681, lay_var_636) / speccomb_planck_var_670, oneminus_var_615)
      specmult_planck_var_668 = 8.0D0 * specparm_planck_var_669
      jpl_var_635 = 1 + INT(specmult_planck_var_668)
      fpl_var_667 = ((specmult_planck_var_668) - AINT((specmult_planck_var_668)))
      ind0_var_629 = ((jp_var_612(jl_var_681, lay_var_636) - 1) * 5 + (jt_var_613(jl_var_681, lay_var_636) - 1)) * nspa_var_216(16) + js_var_633
      ind1_var_630 = (jp_var_612(jl_var_681, lay_var_636) * 5 + (jt1_var_614(jl_var_681, lay_var_636) - 1)) * nspa_var_216(16) + js1_var_634
      inds_var_631 = indself_var_621(jl_var_681, lay_var_636)
      indf_var_632 = indfor_var_625(jl_var_681, lay_var_636)
      IF (specparm_var_661 .LT. 0.125D0) THEN
        p_var_649 = fs_var_659 - 1.0D0
        p4_var_650 = p_var_649 ** 4
        fk0_var_651 = p4_var_650
        fk1_var_652 = 1.0D0 - p_var_649 - 2.0D0 * p4_var_650
        fk2_var_653 = p_var_649 + p4_var_650
        fac000_var_637 = fk0_var_651 * fac00_var_608(jl_var_681, lay_var_636)
        fac100_var_638 = fk1_var_652 * fac00_var_608(jl_var_681, lay_var_636)
        fac200_var_639 = fk2_var_653 * fac00_var_608(jl_var_681, lay_var_636)
        fac010_var_640 = fk0_var_651 * fac10_var_610(jl_var_681, lay_var_636)
        fac110_var_641 = fk1_var_652 * fac10_var_610(jl_var_681, lay_var_636)
        fac210_var_642 = fk2_var_653 * fac10_var_610(jl_var_681, lay_var_636)
      ELSE IF (specparm_var_661 .GT. 0.875D0) THEN
        p_var_649 = - fs_var_659
        p4_var_650 = p_var_649 ** 4
        fk0_var_651 = p4_var_650
        fk1_var_652 = 1.0D0 - p_var_649 - 2.0D0 * p4_var_650
        fk2_var_653 = p_var_649 + p4_var_650
        fac000_var_637 = fk0_var_651 * fac00_var_608(jl_var_681, lay_var_636)
        fac100_var_638 = fk1_var_652 * fac00_var_608(jl_var_681, lay_var_636)
        fac200_var_639 = fk2_var_653 * fac00_var_608(jl_var_681, lay_var_636)
        fac010_var_640 = fk0_var_651 * fac10_var_610(jl_var_681, lay_var_636)
        fac110_var_641 = fk1_var_652 * fac10_var_610(jl_var_681, lay_var_636)
        fac210_var_642 = fk2_var_653 * fac10_var_610(jl_var_681, lay_var_636)
      ELSE
        fac000_var_637 = (1.0D0 - fs_var_659) * fac00_var_608(jl_var_681, lay_var_636)
        fac010_var_640 = (1.0D0 - fs_var_659) * fac10_var_610(jl_var_681, lay_var_636)
        fac100_var_638 = fs_var_659 * fac00_var_608(jl_var_681, lay_var_636)
        fac110_var_641 = fs_var_659 * fac10_var_610(jl_var_681, lay_var_636)
        fac200_var_639 = 0.0D0
        fac210_var_642 = 0.0D0
      END IF
      IF (specparm1_var_665 .LT. 0.125D0) THEN
        p_var_649 = fs1_var_663 - 1.0D0
        p4_var_650 = p_var_649 ** 4
        fk0_var_651 = p4_var_650
        fk1_var_652 = 1.0D0 - p_var_649 - 2.0D0 * p4_var_650
        fk2_var_653 = p_var_649 + p4_var_650
        fac001_var_643 = fk0_var_651 * fac01_var_609(jl_var_681, lay_var_636)
        fac101_var_644 = fk1_var_652 * fac01_var_609(jl_var_681, lay_var_636)
        fac201_var_645 = fk2_var_653 * fac01_var_609(jl_var_681, lay_var_636)
        fac011_var_646 = fk0_var_651 * fac11_var_611(jl_var_681, lay_var_636)
        fac111_var_647 = fk1_var_652 * fac11_var_611(jl_var_681, lay_var_636)
        fac211_var_648 = fk2_var_653 * fac11_var_611(jl_var_681, lay_var_636)
      ELSE IF (specparm1_var_665 .GT. 0.875D0) THEN
        p_var_649 = - fs1_var_663
        p4_var_650 = p_var_649 ** 4
        fk0_var_651 = p4_var_650
        fk1_var_652 = 1.0D0 - p_var_649 - 2.0D0 * p4_var_650
        fk2_var_653 = p_var_649 + p4_var_650
        fac001_var_643 = fk0_var_651 * fac01_var_609(jl_var_681, lay_var_636)
        fac101_var_644 = fk1_var_652 * fac01_var_609(jl_var_681, lay_var_636)
        fac201_var_645 = fk2_var_653 * fac01_var_609(jl_var_681, lay_var_636)
        fac011_var_646 = fk0_var_651 * fac11_var_611(jl_var_681, lay_var_636)
        fac111_var_647 = fk1_var_652 * fac11_var_611(jl_var_681, lay_var_636)
        fac211_var_648 = fk2_var_653 * fac11_var_611(jl_var_681, lay_var_636)
      ELSE
        fac001_var_643 = (1.0D0 - fs1_var_663) * fac01_var_609(jl_var_681, lay_var_636)
        fac011_var_646 = (1.0D0 - fs1_var_663) * fac11_var_611(jl_var_681, lay_var_636)
        fac101_var_644 = fs1_var_663 * fac01_var_609(jl_var_681, lay_var_636)
        fac111_var_647 = fs1_var_663 * fac11_var_611(jl_var_681, lay_var_636)
        fac201_var_645 = 0.0D0
        fac211_var_648 = 0.0D0
      END IF
      IF (specparm_var_661 .LT. 0.125D0) THEN
        tau_major_var_657(1 : ng16) = speccomb_var_662 * (fac000_var_637 * absa_var_149(ind0_var_629, 1 : 2) + fac100_var_638 * absa_var_149(ind0_var_629 + 1, 1 : 2) + fac200_var_639 * absa_var_149(ind0_var_629 + 2, 1 : 2) + fac010_var_640 * absa_var_149(ind0_var_629 + 9, 1 : 2) + fac110_var_641 * absa_var_149(ind0_var_629 + 10, 1 : 2) + fac210_var_642 * absa_var_149(ind0_var_629 + 11, 1 : 2))
      ELSE IF (specparm_var_661 .GT. 0.875D0) THEN
        tau_major_var_657(1 : ng16) = speccomb_var_662 * (fac200_var_639 * absa_var_149(ind0_var_629 - 1, 1 : 2) + fac100_var_638 * absa_var_149(ind0_var_629, 1 : 2) + fac000_var_637 * absa_var_149(ind0_var_629 + 1, 1 : 2) + fac210_var_642 * absa_var_149(ind0_var_629 + 8, 1 : 2) + fac110_var_641 * absa_var_149(ind0_var_629 + 9, 1 : 2) + fac010_var_640 * absa_var_149(ind0_var_629 + 10, 1 : 2))
      ELSE
        tau_major_var_657(1 : ng16) = speccomb_var_662 * (fac000_var_637 * absa_var_149(ind0_var_629, 1 : 2) + fac100_var_638 * absa_var_149(ind0_var_629 + 1, 1 : 2) + fac010_var_640 * absa_var_149(ind0_var_629 + 9, 1 : 2) + fac110_var_641 * absa_var_149(ind0_var_629 + 10, 1 : 2))
      END IF
      IF (specparm1_var_665 .LT. 0.125D0) THEN
        tau_major1_var_658(1 : ng16) = speccomb1_var_666 * (fac001_var_643 * absa_var_149(ind1_var_630, 1 : 2) + fac101_var_644 * absa_var_149(ind1_var_630 + 1, 1 : 2) + fac201_var_645 * absa_var_149(ind1_var_630 + 2, 1 : 2) + fac011_var_646 * absa_var_149(ind1_var_630 + 9, 1 : 2) + fac111_var_647 * absa_var_149(ind1_var_630 + 10, 1 : 2) + fac211_var_648 * absa_var_149(ind1_var_630 + 11, 1 : 2))
      ELSE IF (specparm1_var_665 .GT. 0.875D0) THEN
        tau_major1_var_658(1 : ng16) = speccomb1_var_666 * (fac201_var_645 * absa_var_149(ind1_var_630 - 1, 1 : 2) + fac101_var_644 * absa_var_149(ind1_var_630, 1 : 2) + fac001_var_643 * absa_var_149(ind1_var_630 + 1, 1 : 2) + fac211_var_648 * absa_var_149(ind1_var_630 + 8, 1 : 2) + fac111_var_647 * absa_var_149(ind1_var_630 + 9, 1 : 2) + fac011_var_646 * absa_var_149(ind1_var_630 + 10, 1 : 2))
      ELSE
        tau_major1_var_658(1 : ng16) = speccomb1_var_666 * (fac001_var_643 * absa_var_149(ind1_var_630, 1 : 2) + fac101_var_644 * absa_var_149(ind1_var_630 + 1, 1 : 2) + fac011_var_646 * absa_var_149(ind1_var_630 + 9, 1 : 2) + fac111_var_647 * absa_var_149(ind1_var_630 + 10, 1 : 2))
      END IF
      DO ig_var_628 = 1, 2
        tauself_var_656 = selffac_var_619(jl_var_681, lay_var_636) * (selfref_var_151(inds_var_631, ig_var_628) + selffrac_var_620(jl_var_681, lay_var_636) * (selfref_var_151(inds_var_631 + 1, ig_var_628) - selfref_var_151(inds_var_631, ig_var_628)))
        taufor_var_655 = forfac_var_626(jl_var_681, lay_var_636) * (forref_var_152(indf_var_632, ig_var_628) + forfrac_var_627(jl_var_681, lay_var_636) * (forref_var_152(indf_var_632 + 1, ig_var_628) - forref_var_152(indf_var_632, ig_var_628)))
        taug_var_606(jl_var_681, 138 + ig_var_628, lay_var_636) = tau_major_var_657(ig_var_628) + tau_major1_var_658(ig_var_628) + tauself_var_656 + taufor_var_655
        fracs_var_622(jl_var_681, 138 + ig_var_628, lay_var_636) = fracrefa_var_147(ig_var_628, jpl_var_635) + fpl_var_667 * (fracrefa_var_147(ig_var_628, jpl_var_635 + 1) - fracrefa_var_147(ig_var_628, jpl_var_635))
      END DO
    END DO
  END DO
  DO lay_var_636 = laytrop_max_var_672 + 1, klev_var_605
    DO jl_var_681 = kidia_var_603, kfdia_var_604
      ind0_var_629 = ((jp_var_612(jl_var_681, lay_var_636) - 13) * 5 + (jt_var_613(jl_var_681, lay_var_636) - 1)) * nspb_var_217(16) + 1
      ind1_var_630 = ((jp_var_612(jl_var_681, lay_var_636) - 12) * 5 + (jt1_var_614(jl_var_681, lay_var_636) - 1)) * nspb_var_217(16) + 1
      DO ig_var_628 = 1, 2
        taug_var_606(jl_var_681, 138 + ig_var_628, lay_var_636) = colch4_var_617(jl_var_681, lay_var_636) * (fac00_var_608(jl_var_681, lay_var_636) * absb_var_150(ind0_var_629, ig_var_628) + fac10_var_610(jl_var_681, lay_var_636) * absb_var_150(ind0_var_629 + 1, ig_var_628) + fac01_var_609(jl_var_681, lay_var_636) * absb_var_150(ind1_var_630, ig_var_628) + fac11_var_611(jl_var_681, lay_var_636) * absb_var_150(ind1_var_630 + 1, ig_var_628))
        fracs_var_622(jl_var_681, 138 + ig_var_628, lay_var_636) = fracrefb_var_148(ig_var_628)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_672 /= laytrop_min_var_671) THEN
    DO lay_var_636 = laytrop_min_var_671 + 1, laytrop_max_var_672
      ixc0_var_678 = ixc_var_673(lay_var_636)
      DO ixp_var_679 = 1, ixc0_var_678
        jl_var_681 = ixlow_var_674(ixp_var_679, lay_var_636)
        speccomb_var_662 = colh2o_var_616(jl_var_681, lay_var_636) + rat_h2och4_var_623(jl_var_681, lay_var_636) * colch4_var_617(jl_var_681, lay_var_636)
        specparm_var_661 = MIN(colh2o_var_616(jl_var_681, lay_var_636) / speccomb_var_662, oneminus_var_615)
        specmult_var_660 = 8.0D0 * (specparm_var_661)
        js_var_633 = 1 + INT(specmult_var_660)
        fs_var_659 = ((specmult_var_660) - AINT((specmult_var_660)))
        speccomb1_var_666 = colh2o_var_616(jl_var_681, lay_var_636) + rat_h2och4_1_var_624(jl_var_681, lay_var_636) * colch4_var_617(jl_var_681, lay_var_636)
        specparm1_var_665 = MIN(colh2o_var_616(jl_var_681, lay_var_636) / speccomb1_var_666, oneminus_var_615)
        specmult1_var_664 = 8.0D0 * (specparm1_var_665)
        js1_var_634 = 1 + INT(specmult1_var_664)
        fs1_var_663 = ((specmult1_var_664) - AINT((specmult1_var_664)))
        speccomb_planck_var_670 = colh2o_var_616(jl_var_681, lay_var_636) + refrat_planck_a_var_654 * colch4_var_617(jl_var_681, lay_var_636)
        specparm_planck_var_669 = MIN(colh2o_var_616(jl_var_681, lay_var_636) / speccomb_planck_var_670, oneminus_var_615)
        specmult_planck_var_668 = 8.0D0 * specparm_planck_var_669
        jpl_var_635 = 1 + INT(specmult_planck_var_668)
        fpl_var_667 = ((specmult_planck_var_668) - AINT((specmult_planck_var_668)))
        ind0_var_629 = ((jp_var_612(jl_var_681, lay_var_636) - 1) * 5 + (jt_var_613(jl_var_681, lay_var_636) - 1)) * nspa_var_216(16) + js_var_633
        ind1_var_630 = (jp_var_612(jl_var_681, lay_var_636) * 5 + (jt1_var_614(jl_var_681, lay_var_636) - 1)) * nspa_var_216(16) + js1_var_634
        inds_var_631 = indself_var_621(jl_var_681, lay_var_636)
        indf_var_632 = indfor_var_625(jl_var_681, lay_var_636)
        IF (specparm_var_661 .LT. 0.125D0) THEN
          p_var_649 = fs_var_659 - 1.0D0
          p4_var_650 = p_var_649 ** 4
          fk0_var_651 = p4_var_650
          fk1_var_652 = 1.0D0 - p_var_649 - 2.0D0 * p4_var_650
          fk2_var_653 = p_var_649 + p4_var_650
          fac000_var_637 = fk0_var_651 * fac00_var_608(jl_var_681, lay_var_636)
          fac100_var_638 = fk1_var_652 * fac00_var_608(jl_var_681, lay_var_636)
          fac200_var_639 = fk2_var_653 * fac00_var_608(jl_var_681, lay_var_636)
          fac010_var_640 = fk0_var_651 * fac10_var_610(jl_var_681, lay_var_636)
          fac110_var_641 = fk1_var_652 * fac10_var_610(jl_var_681, lay_var_636)
          fac210_var_642 = fk2_var_653 * fac10_var_610(jl_var_681, lay_var_636)
        ELSE IF (specparm_var_661 .GT. 0.875D0) THEN
          p_var_649 = - fs_var_659
          p4_var_650 = p_var_649 ** 4
          fk0_var_651 = p4_var_650
          fk1_var_652 = 1.0D0 - p_var_649 - 2.0D0 * p4_var_650
          fk2_var_653 = p_var_649 + p4_var_650
          fac000_var_637 = fk0_var_651 * fac00_var_608(jl_var_681, lay_var_636)
          fac100_var_638 = fk1_var_652 * fac00_var_608(jl_var_681, lay_var_636)
          fac200_var_639 = fk2_var_653 * fac00_var_608(jl_var_681, lay_var_636)
          fac010_var_640 = fk0_var_651 * fac10_var_610(jl_var_681, lay_var_636)
          fac110_var_641 = fk1_var_652 * fac10_var_610(jl_var_681, lay_var_636)
          fac210_var_642 = fk2_var_653 * fac10_var_610(jl_var_681, lay_var_636)
        ELSE
          fac000_var_637 = (1.0D0 - fs_var_659) * fac00_var_608(jl_var_681, lay_var_636)
          fac010_var_640 = (1.0D0 - fs_var_659) * fac10_var_610(jl_var_681, lay_var_636)
          fac100_var_638 = fs_var_659 * fac00_var_608(jl_var_681, lay_var_636)
          fac110_var_641 = fs_var_659 * fac10_var_610(jl_var_681, lay_var_636)
          fac200_var_639 = 0.0D0
          fac210_var_642 = 0.0D0
        END IF
        IF (specparm1_var_665 .LT. 0.125D0) THEN
          p_var_649 = fs1_var_663 - 1.0D0
          p4_var_650 = p_var_649 ** 4
          fk0_var_651 = p4_var_650
          fk1_var_652 = 1.0D0 - p_var_649 - 2.0D0 * p4_var_650
          fk2_var_653 = p_var_649 + p4_var_650
          fac001_var_643 = fk0_var_651 * fac01_var_609(jl_var_681, lay_var_636)
          fac101_var_644 = fk1_var_652 * fac01_var_609(jl_var_681, lay_var_636)
          fac201_var_645 = fk2_var_653 * fac01_var_609(jl_var_681, lay_var_636)
          fac011_var_646 = fk0_var_651 * fac11_var_611(jl_var_681, lay_var_636)
          fac111_var_647 = fk1_var_652 * fac11_var_611(jl_var_681, lay_var_636)
          fac211_var_648 = fk2_var_653 * fac11_var_611(jl_var_681, lay_var_636)
        ELSE IF (specparm1_var_665 .GT. 0.875D0) THEN
          p_var_649 = - fs1_var_663
          p4_var_650 = p_var_649 ** 4
          fk0_var_651 = p4_var_650
          fk1_var_652 = 1.0D0 - p_var_649 - 2.0D0 * p4_var_650
          fk2_var_653 = p_var_649 + p4_var_650
          fac001_var_643 = fk0_var_651 * fac01_var_609(jl_var_681, lay_var_636)
          fac101_var_644 = fk1_var_652 * fac01_var_609(jl_var_681, lay_var_636)
          fac201_var_645 = fk2_var_653 * fac01_var_609(jl_var_681, lay_var_636)
          fac011_var_646 = fk0_var_651 * fac11_var_611(jl_var_681, lay_var_636)
          fac111_var_647 = fk1_var_652 * fac11_var_611(jl_var_681, lay_var_636)
          fac211_var_648 = fk2_var_653 * fac11_var_611(jl_var_681, lay_var_636)
        ELSE
          fac001_var_643 = (1.0D0 - fs1_var_663) * fac01_var_609(jl_var_681, lay_var_636)
          fac011_var_646 = (1.0D0 - fs1_var_663) * fac11_var_611(jl_var_681, lay_var_636)
          fac101_var_644 = fs1_var_663 * fac01_var_609(jl_var_681, lay_var_636)
          fac111_var_647 = fs1_var_663 * fac11_var_611(jl_var_681, lay_var_636)
          fac201_var_645 = 0.0D0
          fac211_var_648 = 0.0D0
        END IF
        IF (specparm_var_661 .LT. 0.125D0) THEN
          tau_major_var_657(1 : ng16) = speccomb_var_662 * (fac000_var_637 * absa_var_149(ind0_var_629, 1 : 2) + fac100_var_638 * absa_var_149(ind0_var_629 + 1, 1 : 2) + fac200_var_639 * absa_var_149(ind0_var_629 + 2, 1 : 2) + fac010_var_640 * absa_var_149(ind0_var_629 + 9, 1 : 2) + fac110_var_641 * absa_var_149(ind0_var_629 + 10, 1 : 2) + fac210_var_642 * absa_var_149(ind0_var_629 + 11, 1 : 2))
        ELSE IF (specparm_var_661 .GT. 0.875D0) THEN
          tau_major_var_657(1 : ng16) = speccomb_var_662 * (fac200_var_639 * absa_var_149(ind0_var_629 - 1, 1 : 2) + fac100_var_638 * absa_var_149(ind0_var_629, 1 : 2) + fac000_var_637 * absa_var_149(ind0_var_629 + 1, 1 : 2) + fac210_var_642 * absa_var_149(ind0_var_629 + 8, 1 : 2) + fac110_var_641 * absa_var_149(ind0_var_629 + 9, 1 : 2) + fac010_var_640 * absa_var_149(ind0_var_629 + 10, 1 : 2))
        ELSE
          tau_major_var_657(1 : ng16) = speccomb_var_662 * (fac000_var_637 * absa_var_149(ind0_var_629, 1 : 2) + fac100_var_638 * absa_var_149(ind0_var_629 + 1, 1 : 2) + fac010_var_640 * absa_var_149(ind0_var_629 + 9, 1 : 2) + fac110_var_641 * absa_var_149(ind0_var_629 + 10, 1 : 2))
        END IF
        IF (specparm1_var_665 .LT. 0.125D0) THEN
          tau_major1_var_658(1 : ng16) = speccomb1_var_666 * (fac001_var_643 * absa_var_149(ind1_var_630, 1 : 2) + fac101_var_644 * absa_var_149(ind1_var_630 + 1, 1 : 2) + fac201_var_645 * absa_var_149(ind1_var_630 + 2, 1 : 2) + fac011_var_646 * absa_var_149(ind1_var_630 + 9, 1 : 2) + fac111_var_647 * absa_var_149(ind1_var_630 + 10, 1 : 2) + fac211_var_648 * absa_var_149(ind1_var_630 + 11, 1 : 2))
        ELSE IF (specparm1_var_665 .GT. 0.875D0) THEN
          tau_major1_var_658(1 : ng16) = speccomb1_var_666 * (fac201_var_645 * absa_var_149(ind1_var_630 - 1, 1 : 2) + fac101_var_644 * absa_var_149(ind1_var_630, 1 : 2) + fac001_var_643 * absa_var_149(ind1_var_630 + 1, 1 : 2) + fac211_var_648 * absa_var_149(ind1_var_630 + 8, 1 : 2) + fac111_var_647 * absa_var_149(ind1_var_630 + 9, 1 : 2) + fac011_var_646 * absa_var_149(ind1_var_630 + 10, 1 : 2))
        ELSE
          tau_major1_var_658(1 : ng16) = speccomb1_var_666 * (fac001_var_643 * absa_var_149(ind1_var_630, 1 : 2) + fac101_var_644 * absa_var_149(ind1_var_630 + 1, 1 : 2) + fac011_var_646 * absa_var_149(ind1_var_630 + 9, 1 : 2) + fac111_var_647 * absa_var_149(ind1_var_630 + 10, 1 : 2))
        END IF
        DO ig_var_628 = 1, 2
          tauself_var_656 = selffac_var_619(jl_var_681, lay_var_636) * (selfref_var_151(inds_var_631, ig_var_628) + selffrac_var_620(jl_var_681, lay_var_636) * (selfref_var_151(inds_var_631 + 1, ig_var_628) - selfref_var_151(inds_var_631, ig_var_628)))
          taufor_var_655 = forfac_var_626(jl_var_681, lay_var_636) * (forref_var_152(indf_var_632, ig_var_628) + forfrac_var_627(jl_var_681, lay_var_636) * (forref_var_152(indf_var_632 + 1, ig_var_628) - forref_var_152(indf_var_632, ig_var_628)))
          taug_var_606(jl_var_681, 138 + ig_var_628, lay_var_636) = tau_major_var_657(ig_var_628) + tau_major1_var_658(ig_var_628) + tauself_var_656 + taufor_var_655
          fracs_var_622(jl_var_681, 138 + ig_var_628, lay_var_636) = fracrefa_var_147(ig_var_628, jpl_var_635) + fpl_var_667 * (fracrefa_var_147(ig_var_628, jpl_var_635 + 1) - fracrefa_var_147(ig_var_628, jpl_var_635))
        END DO
      END DO
      ixc0_var_678 = kfdia_var_604 - kidia_var_603 + 1 - ixc0_var_678
      DO ixp_var_679 = 1, ixc0_var_678
        jl_var_681 = ixhigh_var_675(ixp_var_679, lay_var_636)
        ind0_var_629 = ((jp_var_612(jl_var_681, lay_var_636) - 13) * 5 + (jt_var_613(jl_var_681, lay_var_636) - 1)) * nspb_var_217(16) + 1
        ind1_var_630 = ((jp_var_612(jl_var_681, lay_var_636) - 12) * 5 + (jt1_var_614(jl_var_681, lay_var_636) - 1)) * nspb_var_217(16) + 1
        DO ig_var_628 = 1, 2
          taug_var_606(jl_var_681, 138 + ig_var_628, lay_var_636) = colch4_var_617(jl_var_681, lay_var_636) * (fac00_var_608(jl_var_681, lay_var_636) * absb_var_150(ind0_var_629, ig_var_628) + fac10_var_610(jl_var_681, lay_var_636) * absb_var_150(ind0_var_629 + 1, ig_var_628) + fac01_var_609(jl_var_681, lay_var_636) * absb_var_150(ind1_var_630, ig_var_628) + fac11_var_611(jl_var_681, lay_var_636) * absb_var_150(ind1_var_630 + 1, ig_var_628))
          fracs_var_622(jl_var_681, 138 + ig_var_628, lay_var_636) = fracrefb_var_148(ig_var_628)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol16
SUBROUTINE srtm_taumol20(kidia_var_682, kfdia_var_683, klev_var_684, p_fac00_var_685, p_fac01_var_686, p_fac10_var_687, p_fac11_var_688, k_jp_var_689, k_jt_var_690, k_jt1_var_691, p_colh2o_var_692, p_colch4_var_693, p_colmol_var_694, k_laytrop_var_695, p_selffac_var_696, p_selffrac_var_697, k_indself_var_698, p_forfac_var_699, p_forfrac_var_700, k_indfor_var_701, p_sfluxzen_var_702, p_taug_var_703, p_taur_var_704, prmu0_var_705)
  USE yoesrta20, ONLY: absa_var_251, absb_var_252, absch4c, forrefc_var_254, layreffr_var_250, rayl_var_249, selfrefc_var_253, sfluxrefc_var_255
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_682, kfdia_var_683
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_684
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_685(kidia_var_682 : kfdia_var_683, klev_var_684)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_686(kidia_var_682 : kfdia_var_683, klev_var_684)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_687(kidia_var_682 : kfdia_var_683, klev_var_684)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_688(kidia_var_682 : kfdia_var_683, klev_var_684)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_689(kidia_var_682 : kfdia_var_683, klev_var_684)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_690(kidia_var_682 : kfdia_var_683, klev_var_684)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_691(kidia_var_682 : kfdia_var_683, klev_var_684)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_692(kidia_var_682 : kfdia_var_683, klev_var_684)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_693(kidia_var_682 : kfdia_var_683, klev_var_684)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_694(kidia_var_682 : kfdia_var_683, klev_var_684)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_695(kidia_var_682 : kfdia_var_683)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_696(kidia_var_682 : kfdia_var_683, klev_var_684)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_697(kidia_var_682 : kfdia_var_683, klev_var_684)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_698(kidia_var_682 : kfdia_var_683, klev_var_684)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_699(kidia_var_682 : kfdia_var_683, klev_var_684)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_700(kidia_var_682 : kfdia_var_683, klev_var_684)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_701(kidia_var_682 : kfdia_var_683, klev_var_684)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_702(kidia_var_682 : kfdia_var_683, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_703(kidia_var_682 : kfdia_var_683, klev_var_684, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_704(kidia_var_682 : kfdia_var_683, klev_var_684, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_705(kidia_var_682 : kfdia_var_683)
  INTEGER(KIND = 4) :: ig_var_706, ind0_var_707, ind1_var_708, inds_var_709, indf_var_710, i_lay_var_711, i_laysolfr_var_712(kidia_var_682 : kfdia_var_683), i_nlayers_var_713, iplon_var_714
  INTEGER(KIND = 4) :: laytrop_min_var_715, laytrop_max_var_716
  REAL(KIND = 8) :: z_tauray_var_717
  laytrop_min_var_715 = MINVAL(k_laytrop_var_695(kidia_var_682 : kfdia_var_683))
  laytrop_max_var_716 = MAXVAL(k_laytrop_var_695(kidia_var_682 : kfdia_var_683))
  i_nlayers_var_713 = klev_var_684
  DO iplon_var_714 = kidia_var_682, kfdia_var_683
    i_laysolfr_var_712(iplon_var_714) = k_laytrop_var_695(iplon_var_714)
  END DO
  DO i_lay_var_711 = 1, laytrop_min_var_715
    DO iplon_var_714 = kidia_var_682, kfdia_var_683
      IF (k_jp_var_689(iplon_var_714, i_lay_var_711) < layreffr_var_250 .AND. k_jp_var_689(iplon_var_714, i_lay_var_711 + 1) >= layreffr_var_250) i_laysolfr_var_712(iplon_var_714) = MIN(i_lay_var_711 + 1, k_laytrop_var_695(iplon_var_714))
      ind0_var_707 = ((k_jp_var_689(iplon_var_714, i_lay_var_711) - 1) * 5 + (k_jt_var_690(iplon_var_714, i_lay_var_711) - 1)) * nspa_var_313(20) + 1
      ind1_var_708 = (k_jp_var_689(iplon_var_714, i_lay_var_711) * 5 + (k_jt1_var_691(iplon_var_714, i_lay_var_711) - 1)) * nspa_var_313(20) + 1
      inds_var_709 = k_indself_var_698(iplon_var_714, i_lay_var_711)
      indf_var_710 = k_indfor_var_701(iplon_var_714, i_lay_var_711)
      z_tauray_var_717 = p_colmol_var_694(iplon_var_714, i_lay_var_711) * rayl_var_249
      DO ig_var_706 = 1, 10
        p_taug_var_703(iplon_var_714, i_lay_var_711, ig_var_706) = p_colh2o_var_692(iplon_var_714, i_lay_var_711) * ((p_fac00_var_685(iplon_var_714, i_lay_var_711) * absa_var_251(ind0_var_707, ig_var_706) + p_fac10_var_687(iplon_var_714, i_lay_var_711) * absa_var_251(ind0_var_707 + 1, ig_var_706) + p_fac01_var_686(iplon_var_714, i_lay_var_711) * absa_var_251(ind1_var_708, ig_var_706) + p_fac11_var_688(iplon_var_714, i_lay_var_711) * absa_var_251(ind1_var_708 + 1, ig_var_706)) + p_selffac_var_696(iplon_var_714, i_lay_var_711) * (selfrefc_var_253(inds_var_709, ig_var_706) + p_selffrac_var_697(iplon_var_714, i_lay_var_711) * (selfrefc_var_253(inds_var_709 + 1, ig_var_706) - selfrefc_var_253(inds_var_709, ig_var_706))) + p_forfac_var_699(iplon_var_714, i_lay_var_711) * (forrefc_var_254(indf_var_710, ig_var_706) + p_forfrac_var_700(iplon_var_714, i_lay_var_711) * (forrefc_var_254(indf_var_710 + 1, ig_var_706) - forrefc_var_254(indf_var_710, ig_var_706)))) + p_colch4_var_693(iplon_var_714, i_lay_var_711) * absch4c(ig_var_706)
        p_taur_var_704(iplon_var_714, i_lay_var_711, ig_var_706) = z_tauray_var_717
        IF (i_lay_var_711 == i_laysolfr_var_712(iplon_var_714)) p_sfluxzen_var_702(iplon_var_714, ig_var_706) = sfluxrefc_var_255(ig_var_706)
      END DO
    END DO
  END DO
  DO i_lay_var_711 = laytrop_min_var_715 + 1, laytrop_max_var_716
    DO iplon_var_714 = kidia_var_682, kfdia_var_683
      IF (i_lay_var_711 <= k_laytrop_var_695(iplon_var_714)) THEN
        IF (k_jp_var_689(iplon_var_714, i_lay_var_711) < layreffr_var_250 .AND. k_jp_var_689(iplon_var_714, i_lay_var_711 + 1) >= layreffr_var_250) i_laysolfr_var_712(iplon_var_714) = MIN(i_lay_var_711 + 1, k_laytrop_var_695(iplon_var_714))
        ind0_var_707 = ((k_jp_var_689(iplon_var_714, i_lay_var_711) - 1) * 5 + (k_jt_var_690(iplon_var_714, i_lay_var_711) - 1)) * nspa_var_313(20) + 1
        ind1_var_708 = (k_jp_var_689(iplon_var_714, i_lay_var_711) * 5 + (k_jt1_var_691(iplon_var_714, i_lay_var_711) - 1)) * nspa_var_313(20) + 1
        inds_var_709 = k_indself_var_698(iplon_var_714, i_lay_var_711)
        indf_var_710 = k_indfor_var_701(iplon_var_714, i_lay_var_711)
        z_tauray_var_717 = p_colmol_var_694(iplon_var_714, i_lay_var_711) * rayl_var_249
        DO ig_var_706 = 1, 10
          p_taug_var_703(iplon_var_714, i_lay_var_711, ig_var_706) = p_colh2o_var_692(iplon_var_714, i_lay_var_711) * ((p_fac00_var_685(iplon_var_714, i_lay_var_711) * absa_var_251(ind0_var_707, ig_var_706) + p_fac10_var_687(iplon_var_714, i_lay_var_711) * absa_var_251(ind0_var_707 + 1, ig_var_706) + p_fac01_var_686(iplon_var_714, i_lay_var_711) * absa_var_251(ind1_var_708, ig_var_706) + p_fac11_var_688(iplon_var_714, i_lay_var_711) * absa_var_251(ind1_var_708 + 1, ig_var_706)) + p_selffac_var_696(iplon_var_714, i_lay_var_711) * (selfrefc_var_253(inds_var_709, ig_var_706) + p_selffrac_var_697(iplon_var_714, i_lay_var_711) * (selfrefc_var_253(inds_var_709 + 1, ig_var_706) - selfrefc_var_253(inds_var_709, ig_var_706))) + p_forfac_var_699(iplon_var_714, i_lay_var_711) * (forrefc_var_254(indf_var_710, ig_var_706) + p_forfrac_var_700(iplon_var_714, i_lay_var_711) * (forrefc_var_254(indf_var_710 + 1, ig_var_706) - forrefc_var_254(indf_var_710, ig_var_706)))) + p_colch4_var_693(iplon_var_714, i_lay_var_711) * absch4c(ig_var_706)
          p_taur_var_704(iplon_var_714, i_lay_var_711, ig_var_706) = z_tauray_var_717
          IF (i_lay_var_711 == i_laysolfr_var_712(iplon_var_714)) p_sfluxzen_var_702(iplon_var_714, ig_var_706) = sfluxrefc_var_255(ig_var_706)
        END DO
      ELSE
        ind0_var_707 = ((k_jp_var_689(iplon_var_714, i_lay_var_711) - 13) * 5 + (k_jt_var_690(iplon_var_714, i_lay_var_711) - 1)) * nspb_var_314(20) + 1
        ind1_var_708 = ((k_jp_var_689(iplon_var_714, i_lay_var_711) - 12) * 5 + (k_jt1_var_691(iplon_var_714, i_lay_var_711) - 1)) * nspb_var_314(20) + 1
        indf_var_710 = k_indfor_var_701(iplon_var_714, i_lay_var_711)
        z_tauray_var_717 = p_colmol_var_694(iplon_var_714, i_lay_var_711) * rayl_var_249
        DO ig_var_706 = 1, 10
          p_taug_var_703(iplon_var_714, i_lay_var_711, ig_var_706) = p_colh2o_var_692(iplon_var_714, i_lay_var_711) * (p_fac00_var_685(iplon_var_714, i_lay_var_711) * absb_var_252(ind0_var_707, ig_var_706) + p_fac10_var_687(iplon_var_714, i_lay_var_711) * absb_var_252(ind0_var_707 + 1, ig_var_706) + p_fac01_var_686(iplon_var_714, i_lay_var_711) * absb_var_252(ind1_var_708, ig_var_706) + p_fac11_var_688(iplon_var_714, i_lay_var_711) * absb_var_252(ind1_var_708 + 1, ig_var_706) + p_forfac_var_699(iplon_var_714, i_lay_var_711) * (forrefc_var_254(indf_var_710, ig_var_706) + p_forfrac_var_700(iplon_var_714, i_lay_var_711) * (forrefc_var_254(indf_var_710 + 1, ig_var_706) - forrefc_var_254(indf_var_710, ig_var_706)))) + p_colch4_var_693(iplon_var_714, i_lay_var_711) * absch4c(ig_var_706)
          p_taur_var_704(iplon_var_714, i_lay_var_711, ig_var_706) = z_tauray_var_717
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_711 = laytrop_max_var_716 + 1, i_nlayers_var_713
    DO iplon_var_714 = kidia_var_682, kfdia_var_683
      ind0_var_707 = ((k_jp_var_689(iplon_var_714, i_lay_var_711) - 13) * 5 + (k_jt_var_690(iplon_var_714, i_lay_var_711) - 1)) * nspb_var_314(20) + 1
      ind1_var_708 = ((k_jp_var_689(iplon_var_714, i_lay_var_711) - 12) * 5 + (k_jt1_var_691(iplon_var_714, i_lay_var_711) - 1)) * nspb_var_314(20) + 1
      indf_var_710 = k_indfor_var_701(iplon_var_714, i_lay_var_711)
      z_tauray_var_717 = p_colmol_var_694(iplon_var_714, i_lay_var_711) * rayl_var_249
      DO ig_var_706 = 1, 10
        p_taug_var_703(iplon_var_714, i_lay_var_711, ig_var_706) = p_colh2o_var_692(iplon_var_714, i_lay_var_711) * (p_fac00_var_685(iplon_var_714, i_lay_var_711) * absb_var_252(ind0_var_707, ig_var_706) + p_fac10_var_687(iplon_var_714, i_lay_var_711) * absb_var_252(ind0_var_707 + 1, ig_var_706) + p_fac01_var_686(iplon_var_714, i_lay_var_711) * absb_var_252(ind1_var_708, ig_var_706) + p_fac11_var_688(iplon_var_714, i_lay_var_711) * absb_var_252(ind1_var_708 + 1, ig_var_706) + p_forfac_var_699(iplon_var_714, i_lay_var_711) * (forrefc_var_254(indf_var_710, ig_var_706) + p_forfrac_var_700(iplon_var_714, i_lay_var_711) * (forrefc_var_254(indf_var_710 + 1, ig_var_706) - forrefc_var_254(indf_var_710, ig_var_706)))) + p_colch4_var_693(iplon_var_714, i_lay_var_711) * absch4c(ig_var_706)
        p_taur_var_704(iplon_var_714, i_lay_var_711, ig_var_706) = z_tauray_var_717
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol20
SUBROUTINE srtm_taumol18(kidia_var_718, kfdia_var_719, klev_var_720, p_fac00_var_721, p_fac01_var_722, p_fac10_var_723, p_fac11_var_724, k_jp_var_725, k_jt_var_726, k_jt1_var_727, p_oneminus_var_728, p_colh2o_var_729, p_colch4_var_730, p_colmol_var_731, k_laytrop_var_732, p_selffac_var_733, p_selffrac_var_734, k_indself_var_735, p_forfac_var_736, p_forfrac_var_737, k_indfor_var_738, p_sfluxzen_var_739, p_taug_var_740, p_taur_var_741, prmu0_var_742)
  USE yoesrta18, ONLY: absa_var_236, absb_var_237, forrefc_var_239, layreffr_var_235, rayl_var_233, selfrefc_var_238, sfluxrefc_var_240, strrat_var_234
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_718, kfdia_var_719
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_720
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_721(kidia_var_718 : kfdia_var_719, klev_var_720)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_722(kidia_var_718 : kfdia_var_719, klev_var_720)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_723(kidia_var_718 : kfdia_var_719, klev_var_720)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_724(kidia_var_718 : kfdia_var_719, klev_var_720)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_725(kidia_var_718 : kfdia_var_719, klev_var_720)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_726(kidia_var_718 : kfdia_var_719, klev_var_720)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_727(kidia_var_718 : kfdia_var_719, klev_var_720)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_728(kidia_var_718 : kfdia_var_719)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_729(kidia_var_718 : kfdia_var_719, klev_var_720)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_730(kidia_var_718 : kfdia_var_719, klev_var_720)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_731(kidia_var_718 : kfdia_var_719, klev_var_720)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_732(kidia_var_718 : kfdia_var_719)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_733(kidia_var_718 : kfdia_var_719, klev_var_720)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_734(kidia_var_718 : kfdia_var_719, klev_var_720)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_735(kidia_var_718 : kfdia_var_719, klev_var_720)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_736(kidia_var_718 : kfdia_var_719, klev_var_720)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_737(kidia_var_718 : kfdia_var_719, klev_var_720)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_738(kidia_var_718 : kfdia_var_719, klev_var_720)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_739(kidia_var_718 : kfdia_var_719, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_740(kidia_var_718 : kfdia_var_719, klev_var_720, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_741(kidia_var_718 : kfdia_var_719, klev_var_720, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_742(kidia_var_718 : kfdia_var_719)
  INTEGER(KIND = 4) :: ig_var_743, ind0_var_744, ind1_var_745, inds_var_746, indf_var_747, js_var_748, i_lay_var_749, i_laysolfr_var_750(kidia_var_718 : kfdia_var_719), i_nlayers_var_751, iplon_var_752
  INTEGER(KIND = 4) :: laytrop_min_var_753, laytrop_max_var_754
  REAL(KIND = 8) :: z_fs_var_755, z_speccomb_var_756, z_specmult_var_757, z_specparm_var_758, z_tauray_var_759
  laytrop_min_var_753 = MINVAL(k_laytrop_var_732(kidia_var_718 : kfdia_var_719))
  laytrop_max_var_754 = MAXVAL(k_laytrop_var_732(kidia_var_718 : kfdia_var_719))
  i_nlayers_var_751 = klev_var_720
  DO iplon_var_752 = kidia_var_718, kfdia_var_719
    i_laysolfr_var_750(iplon_var_752) = k_laytrop_var_732(iplon_var_752)
  END DO
  DO i_lay_var_749 = 1, laytrop_min_var_753
    DO iplon_var_752 = kidia_var_718, kfdia_var_719
      IF (k_jp_var_725(iplon_var_752, i_lay_var_749) < layreffr_var_235 .AND. k_jp_var_725(iplon_var_752, i_lay_var_749 + 1) >= layreffr_var_235) i_laysolfr_var_750(iplon_var_752) = MIN(i_lay_var_749 + 1, k_laytrop_var_732(iplon_var_752))
      z_speccomb_var_756 = p_colh2o_var_729(iplon_var_752, i_lay_var_749) + strrat_var_234 * p_colch4_var_730(iplon_var_752, i_lay_var_749)
      z_specparm_var_758 = p_colh2o_var_729(iplon_var_752, i_lay_var_749) / z_speccomb_var_756
      z_specparm_var_758 = MIN(p_oneminus_var_728(iplon_var_752), z_specparm_var_758)
      z_specmult_var_757 = 8.0D0 * (z_specparm_var_758)
      js_var_748 = 1 + INT(z_specmult_var_757)
      z_fs_var_755 = z_specmult_var_757 - AINT(z_specmult_var_757)
      ind0_var_744 = ((k_jp_var_725(iplon_var_752, i_lay_var_749) - 1) * 5 + (k_jt_var_726(iplon_var_752, i_lay_var_749) - 1)) * nspa_var_313(18) + js_var_748
      ind1_var_745 = (k_jp_var_725(iplon_var_752, i_lay_var_749) * 5 + (k_jt1_var_727(iplon_var_752, i_lay_var_749) - 1)) * nspa_var_313(18) + js_var_748
      inds_var_746 = k_indself_var_735(iplon_var_752, i_lay_var_749)
      indf_var_747 = k_indfor_var_738(iplon_var_752, i_lay_var_749)
      z_tauray_var_759 = p_colmol_var_731(iplon_var_752, i_lay_var_749) * rayl_var_233
      DO ig_var_743 = 1, 8
        p_taug_var_740(iplon_var_752, i_lay_var_749, ig_var_743) = z_speccomb_var_756 * ((1.0D0 - z_fs_var_755) * (absa_var_236(ind0_var_744, ig_var_743) * p_fac00_var_721(iplon_var_752, i_lay_var_749) + absa_var_236(ind0_var_744 + 9, ig_var_743) * p_fac10_var_723(iplon_var_752, i_lay_var_749) + absa_var_236(ind1_var_745, ig_var_743) * p_fac01_var_722(iplon_var_752, i_lay_var_749) + absa_var_236(ind1_var_745 + 9, ig_var_743) * p_fac11_var_724(iplon_var_752, i_lay_var_749)) + z_fs_var_755 * (absa_var_236(ind0_var_744 + 1, ig_var_743) * p_fac00_var_721(iplon_var_752, i_lay_var_749) + absa_var_236(ind0_var_744 + 10, ig_var_743) * p_fac10_var_723(iplon_var_752, i_lay_var_749) + absa_var_236(ind1_var_745 + 1, ig_var_743) * p_fac01_var_722(iplon_var_752, i_lay_var_749) + absa_var_236(ind1_var_745 + 10, ig_var_743) * p_fac11_var_724(iplon_var_752, i_lay_var_749))) + p_colh2o_var_729(iplon_var_752, i_lay_var_749) * (p_selffac_var_733(iplon_var_752, i_lay_var_749) * (selfrefc_var_238(inds_var_746, ig_var_743) + p_selffrac_var_734(iplon_var_752, i_lay_var_749) * (selfrefc_var_238(inds_var_746 + 1, ig_var_743) - selfrefc_var_238(inds_var_746, ig_var_743))) + p_forfac_var_736(iplon_var_752, i_lay_var_749) * (forrefc_var_239(indf_var_747, ig_var_743) + p_forfrac_var_737(iplon_var_752, i_lay_var_749) * (forrefc_var_239(indf_var_747 + 1, ig_var_743) - forrefc_var_239(indf_var_747, ig_var_743))))
        IF (i_lay_var_749 == i_laysolfr_var_750(iplon_var_752)) p_sfluxzen_var_739(iplon_var_752, ig_var_743) = sfluxrefc_var_240(ig_var_743, js_var_748) + z_fs_var_755 * (sfluxrefc_var_240(ig_var_743, js_var_748 + 1) - sfluxrefc_var_240(ig_var_743, js_var_748))
        p_taur_var_741(iplon_var_752, i_lay_var_749, ig_var_743) = z_tauray_var_759
      END DO
    END DO
  END DO
  DO i_lay_var_749 = laytrop_min_var_753 + 1, laytrop_max_var_754
    DO iplon_var_752 = kidia_var_718, kfdia_var_719
      IF (i_lay_var_749 <= k_laytrop_var_732(iplon_var_752)) THEN
        IF (k_jp_var_725(iplon_var_752, i_lay_var_749) < layreffr_var_235 .AND. k_jp_var_725(iplon_var_752, i_lay_var_749 + 1) >= layreffr_var_235) i_laysolfr_var_750(iplon_var_752) = MIN(i_lay_var_749 + 1, k_laytrop_var_732(iplon_var_752))
        z_speccomb_var_756 = p_colh2o_var_729(iplon_var_752, i_lay_var_749) + strrat_var_234 * p_colch4_var_730(iplon_var_752, i_lay_var_749)
        z_specparm_var_758 = p_colh2o_var_729(iplon_var_752, i_lay_var_749) / z_speccomb_var_756
        z_specparm_var_758 = MIN(p_oneminus_var_728(iplon_var_752), z_specparm_var_758)
        z_specmult_var_757 = 8.0D0 * (z_specparm_var_758)
        js_var_748 = 1 + INT(z_specmult_var_757)
        z_fs_var_755 = z_specmult_var_757 - AINT(z_specmult_var_757)
        ind0_var_744 = ((k_jp_var_725(iplon_var_752, i_lay_var_749) - 1) * 5 + (k_jt_var_726(iplon_var_752, i_lay_var_749) - 1)) * nspa_var_313(18) + js_var_748
        ind1_var_745 = (k_jp_var_725(iplon_var_752, i_lay_var_749) * 5 + (k_jt1_var_727(iplon_var_752, i_lay_var_749) - 1)) * nspa_var_313(18) + js_var_748
        inds_var_746 = k_indself_var_735(iplon_var_752, i_lay_var_749)
        indf_var_747 = k_indfor_var_738(iplon_var_752, i_lay_var_749)
        z_tauray_var_759 = p_colmol_var_731(iplon_var_752, i_lay_var_749) * rayl_var_233
        DO ig_var_743 = 1, 8
          p_taug_var_740(iplon_var_752, i_lay_var_749, ig_var_743) = z_speccomb_var_756 * ((1.0D0 - z_fs_var_755) * (absa_var_236(ind0_var_744, ig_var_743) * p_fac00_var_721(iplon_var_752, i_lay_var_749) + absa_var_236(ind0_var_744 + 9, ig_var_743) * p_fac10_var_723(iplon_var_752, i_lay_var_749) + absa_var_236(ind1_var_745, ig_var_743) * p_fac01_var_722(iplon_var_752, i_lay_var_749) + absa_var_236(ind1_var_745 + 9, ig_var_743) * p_fac11_var_724(iplon_var_752, i_lay_var_749)) + z_fs_var_755 * (absa_var_236(ind0_var_744 + 1, ig_var_743) * p_fac00_var_721(iplon_var_752, i_lay_var_749) + absa_var_236(ind0_var_744 + 10, ig_var_743) * p_fac10_var_723(iplon_var_752, i_lay_var_749) + absa_var_236(ind1_var_745 + 1, ig_var_743) * p_fac01_var_722(iplon_var_752, i_lay_var_749) + absa_var_236(ind1_var_745 + 10, ig_var_743) * p_fac11_var_724(iplon_var_752, i_lay_var_749))) + p_colh2o_var_729(iplon_var_752, i_lay_var_749) * (p_selffac_var_733(iplon_var_752, i_lay_var_749) * (selfrefc_var_238(inds_var_746, ig_var_743) + p_selffrac_var_734(iplon_var_752, i_lay_var_749) * (selfrefc_var_238(inds_var_746 + 1, ig_var_743) - selfrefc_var_238(inds_var_746, ig_var_743))) + p_forfac_var_736(iplon_var_752, i_lay_var_749) * (forrefc_var_239(indf_var_747, ig_var_743) + p_forfrac_var_737(iplon_var_752, i_lay_var_749) * (forrefc_var_239(indf_var_747 + 1, ig_var_743) - forrefc_var_239(indf_var_747, ig_var_743))))
          IF (i_lay_var_749 == i_laysolfr_var_750(iplon_var_752)) p_sfluxzen_var_739(iplon_var_752, ig_var_743) = sfluxrefc_var_240(ig_var_743, js_var_748) + z_fs_var_755 * (sfluxrefc_var_240(ig_var_743, js_var_748 + 1) - sfluxrefc_var_240(ig_var_743, js_var_748))
          p_taur_var_741(iplon_var_752, i_lay_var_749, ig_var_743) = z_tauray_var_759
        END DO
      ELSE
        ind0_var_744 = ((k_jp_var_725(iplon_var_752, i_lay_var_749) - 13) * 5 + (k_jt_var_726(iplon_var_752, i_lay_var_749) - 1)) * nspb_var_314(18) + 1
        ind1_var_745 = ((k_jp_var_725(iplon_var_752, i_lay_var_749) - 12) * 5 + (k_jt1_var_727(iplon_var_752, i_lay_var_749) - 1)) * nspb_var_314(18) + 1
        z_tauray_var_759 = p_colmol_var_731(iplon_var_752, i_lay_var_749) * rayl_var_233
        DO ig_var_743 = 1, 8
          p_taug_var_740(iplon_var_752, i_lay_var_749, ig_var_743) = p_colch4_var_730(iplon_var_752, i_lay_var_749) * (p_fac00_var_721(iplon_var_752, i_lay_var_749) * absb_var_237(ind0_var_744, ig_var_743) + p_fac10_var_723(iplon_var_752, i_lay_var_749) * absb_var_237(ind0_var_744 + 1, ig_var_743) + p_fac01_var_722(iplon_var_752, i_lay_var_749) * absb_var_237(ind1_var_745, ig_var_743) + p_fac11_var_724(iplon_var_752, i_lay_var_749) * absb_var_237(ind1_var_745 + 1, ig_var_743))
          p_taur_var_741(iplon_var_752, i_lay_var_749, ig_var_743) = z_tauray_var_759
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_749 = laytrop_max_var_754 + 1, i_nlayers_var_751
    DO iplon_var_752 = kidia_var_718, kfdia_var_719
      ind0_var_744 = ((k_jp_var_725(iplon_var_752, i_lay_var_749) - 13) * 5 + (k_jt_var_726(iplon_var_752, i_lay_var_749) - 1)) * nspb_var_314(18) + 1
      ind1_var_745 = ((k_jp_var_725(iplon_var_752, i_lay_var_749) - 12) * 5 + (k_jt1_var_727(iplon_var_752, i_lay_var_749) - 1)) * nspb_var_314(18) + 1
      z_tauray_var_759 = p_colmol_var_731(iplon_var_752, i_lay_var_749) * rayl_var_233
      DO ig_var_743 = 1, 8
        p_taug_var_740(iplon_var_752, i_lay_var_749, ig_var_743) = p_colch4_var_730(iplon_var_752, i_lay_var_749) * (p_fac00_var_721(iplon_var_752, i_lay_var_749) * absb_var_237(ind0_var_744, ig_var_743) + p_fac10_var_723(iplon_var_752, i_lay_var_749) * absb_var_237(ind0_var_744 + 1, ig_var_743) + p_fac01_var_722(iplon_var_752, i_lay_var_749) * absb_var_237(ind1_var_745, ig_var_743) + p_fac11_var_724(iplon_var_752, i_lay_var_749) * absb_var_237(ind1_var_745 + 1, ig_var_743))
        p_taur_var_741(iplon_var_752, i_lay_var_749, ig_var_743) = z_tauray_var_759
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol18
SUBROUTINE srtm_taumol21(kidia_var_760, kfdia_var_761, klev_var_762, p_fac00_var_763, p_fac01_var_764, p_fac10_var_765, p_fac11_var_766, k_jp_var_767, k_jt_var_768, k_jt1_var_769, p_oneminus_var_770, p_colh2o_var_771, p_colco2_var_772, p_colmol_var_773, k_laytrop_var_774, p_selffac_var_775, p_selffrac_var_776, k_indself_var_777, p_forfac_var_778, p_forfrac_var_779, k_indfor_var_780, p_sfluxzen_var_781, p_taug_var_782, p_taur_var_783, prmu0_var_784)
  USE yoesrta21, ONLY: absa_var_259, absb_var_260, forrefc_var_262, layreffr_var_258, rayl_var_256, selfrefc_var_261, sfluxrefc_var_263, strrat_var_257
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_760, kfdia_var_761
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_762
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_763(kidia_var_760 : kfdia_var_761, klev_var_762)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_764(kidia_var_760 : kfdia_var_761, klev_var_762)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_765(kidia_var_760 : kfdia_var_761, klev_var_762)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_766(kidia_var_760 : kfdia_var_761, klev_var_762)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_767(kidia_var_760 : kfdia_var_761, klev_var_762)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_768(kidia_var_760 : kfdia_var_761, klev_var_762)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_769(kidia_var_760 : kfdia_var_761, klev_var_762)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_770(kidia_var_760 : kfdia_var_761)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_771(kidia_var_760 : kfdia_var_761, klev_var_762)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_772(kidia_var_760 : kfdia_var_761, klev_var_762)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_773(kidia_var_760 : kfdia_var_761, klev_var_762)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_774(kidia_var_760 : kfdia_var_761)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_775(kidia_var_760 : kfdia_var_761, klev_var_762)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_776(kidia_var_760 : kfdia_var_761, klev_var_762)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_777(kidia_var_760 : kfdia_var_761, klev_var_762)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_778(kidia_var_760 : kfdia_var_761, klev_var_762)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_779(kidia_var_760 : kfdia_var_761, klev_var_762)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_780(kidia_var_760 : kfdia_var_761, klev_var_762)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_781(kidia_var_760 : kfdia_var_761, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_782(kidia_var_760 : kfdia_var_761, klev_var_762, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_783(kidia_var_760 : kfdia_var_761, klev_var_762, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_784(kidia_var_760 : kfdia_var_761)
  INTEGER(KIND = 4) :: ig_var_785, ind0_var_786, ind1_var_787, inds_var_788, indf_var_789, js_var_790, i_lay_var_791, i_laysolfr_var_792(kidia_var_760 : kfdia_var_761), i_nlayers_var_793, iplon_var_794
  INTEGER(KIND = 4) :: laytrop_min_var_795, laytrop_max_var_796
  REAL(KIND = 8) :: z_fs_var_797, z_speccomb_var_798, z_specmult_var_799, z_specparm_var_800, z_tauray_var_801
  laytrop_min_var_795 = MINVAL(k_laytrop_var_774(kidia_var_760 : kfdia_var_761))
  laytrop_max_var_796 = MAXVAL(k_laytrop_var_774(kidia_var_760 : kfdia_var_761))
  i_nlayers_var_793 = klev_var_762
  DO iplon_var_794 = kidia_var_760, kfdia_var_761
    i_laysolfr_var_792(iplon_var_794) = k_laytrop_var_774(iplon_var_794)
  END DO
  DO i_lay_var_791 = 1, laytrop_min_var_795
    DO iplon_var_794 = kidia_var_760, kfdia_var_761
      IF (k_jp_var_767(iplon_var_794, i_lay_var_791) < layreffr_var_258 .AND. k_jp_var_767(iplon_var_794, i_lay_var_791 + 1) >= layreffr_var_258) i_laysolfr_var_792(iplon_var_794) = MIN(i_lay_var_791 + 1, k_laytrop_var_774(iplon_var_794))
      z_speccomb_var_798 = p_colh2o_var_771(iplon_var_794, i_lay_var_791) + strrat_var_257 * p_colco2_var_772(iplon_var_794, i_lay_var_791)
      z_specparm_var_800 = p_colh2o_var_771(iplon_var_794, i_lay_var_791) / z_speccomb_var_798
      z_specparm_var_800 = MIN(p_oneminus_var_770(iplon_var_794), z_specparm_var_800)
      z_specmult_var_799 = 8.0D0 * (z_specparm_var_800)
      js_var_790 = 1 + INT(z_specmult_var_799)
      z_fs_var_797 = z_specmult_var_799 - AINT(z_specmult_var_799)
      ind0_var_786 = ((k_jp_var_767(iplon_var_794, i_lay_var_791) - 1) * 5 + (k_jt_var_768(iplon_var_794, i_lay_var_791) - 1)) * nspa_var_313(21) + js_var_790
      ind1_var_787 = (k_jp_var_767(iplon_var_794, i_lay_var_791) * 5 + (k_jt1_var_769(iplon_var_794, i_lay_var_791) - 1)) * nspa_var_313(21) + js_var_790
      inds_var_788 = k_indself_var_777(iplon_var_794, i_lay_var_791)
      indf_var_789 = k_indfor_var_780(iplon_var_794, i_lay_var_791)
      z_tauray_var_801 = p_colmol_var_773(iplon_var_794, i_lay_var_791) * rayl_var_256
      DO ig_var_785 = 1, 10
        p_taug_var_782(iplon_var_794, i_lay_var_791, ig_var_785) = z_speccomb_var_798 * ((1.0D0 - z_fs_var_797) * (absa_var_259(ind0_var_786, ig_var_785) * p_fac00_var_763(iplon_var_794, i_lay_var_791) + absa_var_259(ind0_var_786 + 9, ig_var_785) * p_fac10_var_765(iplon_var_794, i_lay_var_791) + absa_var_259(ind1_var_787, ig_var_785) * p_fac01_var_764(iplon_var_794, i_lay_var_791) + absa_var_259(ind1_var_787 + 9, ig_var_785) * p_fac11_var_766(iplon_var_794, i_lay_var_791)) + z_fs_var_797 * (absa_var_259(ind0_var_786 + 1, ig_var_785) * p_fac00_var_763(iplon_var_794, i_lay_var_791) + absa_var_259(ind0_var_786 + 10, ig_var_785) * p_fac10_var_765(iplon_var_794, i_lay_var_791) + absa_var_259(ind1_var_787 + 1, ig_var_785) * p_fac01_var_764(iplon_var_794, i_lay_var_791) + absa_var_259(ind1_var_787 + 10, ig_var_785) * p_fac11_var_766(iplon_var_794, i_lay_var_791))) + p_colh2o_var_771(iplon_var_794, i_lay_var_791) * (p_selffac_var_775(iplon_var_794, i_lay_var_791) * (selfrefc_var_261(inds_var_788, ig_var_785) + p_selffrac_var_776(iplon_var_794, i_lay_var_791) * (selfrefc_var_261(inds_var_788 + 1, ig_var_785) - selfrefc_var_261(inds_var_788, ig_var_785))) + p_forfac_var_778(iplon_var_794, i_lay_var_791) * (forrefc_var_262(indf_var_789, ig_var_785) + p_forfrac_var_779(iplon_var_794, i_lay_var_791) * (forrefc_var_262(indf_var_789 + 1, ig_var_785) - forrefc_var_262(indf_var_789, ig_var_785))))
        IF (i_lay_var_791 == i_laysolfr_var_792(iplon_var_794)) p_sfluxzen_var_781(iplon_var_794, ig_var_785) = sfluxrefc_var_263(ig_var_785, js_var_790) + z_fs_var_797 * (sfluxrefc_var_263(ig_var_785, js_var_790 + 1) - sfluxrefc_var_263(ig_var_785, js_var_790))
        p_taur_var_783(iplon_var_794, i_lay_var_791, ig_var_785) = z_tauray_var_801
      END DO
    END DO
  END DO
  DO i_lay_var_791 = laytrop_min_var_795 + 1, laytrop_max_var_796
    DO iplon_var_794 = kidia_var_760, kfdia_var_761
      IF (i_lay_var_791 <= k_laytrop_var_774(iplon_var_794)) THEN
        IF (k_jp_var_767(iplon_var_794, i_lay_var_791) < layreffr_var_258 .AND. k_jp_var_767(iplon_var_794, i_lay_var_791 + 1) >= layreffr_var_258) i_laysolfr_var_792(iplon_var_794) = MIN(i_lay_var_791 + 1, k_laytrop_var_774(iplon_var_794))
        z_speccomb_var_798 = p_colh2o_var_771(iplon_var_794, i_lay_var_791) + strrat_var_257 * p_colco2_var_772(iplon_var_794, i_lay_var_791)
        z_specparm_var_800 = p_colh2o_var_771(iplon_var_794, i_lay_var_791) / z_speccomb_var_798
        z_specparm_var_800 = MIN(p_oneminus_var_770(iplon_var_794), z_specparm_var_800)
        z_specmult_var_799 = 8.0D0 * (z_specparm_var_800)
        js_var_790 = 1 + INT(z_specmult_var_799)
        z_fs_var_797 = z_specmult_var_799 - AINT(z_specmult_var_799)
        ind0_var_786 = ((k_jp_var_767(iplon_var_794, i_lay_var_791) - 1) * 5 + (k_jt_var_768(iplon_var_794, i_lay_var_791) - 1)) * nspa_var_313(21) + js_var_790
        ind1_var_787 = (k_jp_var_767(iplon_var_794, i_lay_var_791) * 5 + (k_jt1_var_769(iplon_var_794, i_lay_var_791) - 1)) * nspa_var_313(21) + js_var_790
        inds_var_788 = k_indself_var_777(iplon_var_794, i_lay_var_791)
        indf_var_789 = k_indfor_var_780(iplon_var_794, i_lay_var_791)
        z_tauray_var_801 = p_colmol_var_773(iplon_var_794, i_lay_var_791) * rayl_var_256
        DO ig_var_785 = 1, 10
          p_taug_var_782(iplon_var_794, i_lay_var_791, ig_var_785) = z_speccomb_var_798 * ((1.0D0 - z_fs_var_797) * (absa_var_259(ind0_var_786, ig_var_785) * p_fac00_var_763(iplon_var_794, i_lay_var_791) + absa_var_259(ind0_var_786 + 9, ig_var_785) * p_fac10_var_765(iplon_var_794, i_lay_var_791) + absa_var_259(ind1_var_787, ig_var_785) * p_fac01_var_764(iplon_var_794, i_lay_var_791) + absa_var_259(ind1_var_787 + 9, ig_var_785) * p_fac11_var_766(iplon_var_794, i_lay_var_791)) + z_fs_var_797 * (absa_var_259(ind0_var_786 + 1, ig_var_785) * p_fac00_var_763(iplon_var_794, i_lay_var_791) + absa_var_259(ind0_var_786 + 10, ig_var_785) * p_fac10_var_765(iplon_var_794, i_lay_var_791) + absa_var_259(ind1_var_787 + 1, ig_var_785) * p_fac01_var_764(iplon_var_794, i_lay_var_791) + absa_var_259(ind1_var_787 + 10, ig_var_785) * p_fac11_var_766(iplon_var_794, i_lay_var_791))) + p_colh2o_var_771(iplon_var_794, i_lay_var_791) * (p_selffac_var_775(iplon_var_794, i_lay_var_791) * (selfrefc_var_261(inds_var_788, ig_var_785) + p_selffrac_var_776(iplon_var_794, i_lay_var_791) * (selfrefc_var_261(inds_var_788 + 1, ig_var_785) - selfrefc_var_261(inds_var_788, ig_var_785))) + p_forfac_var_778(iplon_var_794, i_lay_var_791) * (forrefc_var_262(indf_var_789, ig_var_785) + p_forfrac_var_779(iplon_var_794, i_lay_var_791) * (forrefc_var_262(indf_var_789 + 1, ig_var_785) - forrefc_var_262(indf_var_789, ig_var_785))))
          IF (i_lay_var_791 == i_laysolfr_var_792(iplon_var_794)) p_sfluxzen_var_781(iplon_var_794, ig_var_785) = sfluxrefc_var_263(ig_var_785, js_var_790) + z_fs_var_797 * (sfluxrefc_var_263(ig_var_785, js_var_790 + 1) - sfluxrefc_var_263(ig_var_785, js_var_790))
          p_taur_var_783(iplon_var_794, i_lay_var_791, ig_var_785) = z_tauray_var_801
        END DO
      ELSE
        z_speccomb_var_798 = p_colh2o_var_771(iplon_var_794, i_lay_var_791) + strrat_var_257 * p_colco2_var_772(iplon_var_794, i_lay_var_791)
        z_specparm_var_800 = p_colh2o_var_771(iplon_var_794, i_lay_var_791) / z_speccomb_var_798
        z_specparm_var_800 = MIN(p_oneminus_var_770(iplon_var_794), z_specparm_var_800)
        z_specmult_var_799 = 4.0D0 * (z_specparm_var_800)
        js_var_790 = 1 + INT(z_specmult_var_799)
        z_fs_var_797 = z_specmult_var_799 - AINT(z_specmult_var_799)
        ind0_var_786 = ((k_jp_var_767(iplon_var_794, i_lay_var_791) - 13) * 5 + (k_jt_var_768(iplon_var_794, i_lay_var_791) - 1)) * nspb_var_314(21) + js_var_790
        ind1_var_787 = ((k_jp_var_767(iplon_var_794, i_lay_var_791) - 12) * 5 + (k_jt1_var_769(iplon_var_794, i_lay_var_791) - 1)) * nspb_var_314(21) + js_var_790
        indf_var_789 = k_indfor_var_780(iplon_var_794, i_lay_var_791)
        z_tauray_var_801 = p_colmol_var_773(iplon_var_794, i_lay_var_791) * rayl_var_256
        DO ig_var_785 = 1, 10
          p_taug_var_782(iplon_var_794, i_lay_var_791, ig_var_785) = z_speccomb_var_798 * ((1.0D0 - z_fs_var_797) * (absb_var_260(ind0_var_786, ig_var_785) * p_fac00_var_763(iplon_var_794, i_lay_var_791) + absb_var_260(ind0_var_786 + 5, ig_var_785) * p_fac10_var_765(iplon_var_794, i_lay_var_791) + absb_var_260(ind1_var_787, ig_var_785) * p_fac01_var_764(iplon_var_794, i_lay_var_791) + absb_var_260(ind1_var_787 + 5, ig_var_785) * p_fac11_var_766(iplon_var_794, i_lay_var_791)) + z_fs_var_797 * (absb_var_260(ind0_var_786 + 1, ig_var_785) * p_fac00_var_763(iplon_var_794, i_lay_var_791) + absb_var_260(ind0_var_786 + 6, ig_var_785) * p_fac10_var_765(iplon_var_794, i_lay_var_791) + absb_var_260(ind1_var_787 + 1, ig_var_785) * p_fac01_var_764(iplon_var_794, i_lay_var_791) + absb_var_260(ind1_var_787 + 6, ig_var_785) * p_fac11_var_766(iplon_var_794, i_lay_var_791))) + p_colh2o_var_771(iplon_var_794, i_lay_var_791) * p_forfac_var_778(iplon_var_794, i_lay_var_791) * (forrefc_var_262(indf_var_789, ig_var_785) + p_forfrac_var_779(iplon_var_794, i_lay_var_791) * (forrefc_var_262(indf_var_789 + 1, ig_var_785) - forrefc_var_262(indf_var_789, ig_var_785)))
          p_taur_var_783(iplon_var_794, i_lay_var_791, ig_var_785) = z_tauray_var_801
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_791 = laytrop_max_var_796 + 1, i_nlayers_var_793
    DO iplon_var_794 = kidia_var_760, kfdia_var_761
      z_speccomb_var_798 = p_colh2o_var_771(iplon_var_794, i_lay_var_791) + strrat_var_257 * p_colco2_var_772(iplon_var_794, i_lay_var_791)
      z_specparm_var_800 = p_colh2o_var_771(iplon_var_794, i_lay_var_791) / z_speccomb_var_798
      z_specparm_var_800 = MIN(p_oneminus_var_770(iplon_var_794), z_specparm_var_800)
      z_specmult_var_799 = 4.0D0 * (z_specparm_var_800)
      js_var_790 = 1 + INT(z_specmult_var_799)
      z_fs_var_797 = z_specmult_var_799 - AINT(z_specmult_var_799)
      ind0_var_786 = ((k_jp_var_767(iplon_var_794, i_lay_var_791) - 13) * 5 + (k_jt_var_768(iplon_var_794, i_lay_var_791) - 1)) * nspb_var_314(21) + js_var_790
      ind1_var_787 = ((k_jp_var_767(iplon_var_794, i_lay_var_791) - 12) * 5 + (k_jt1_var_769(iplon_var_794, i_lay_var_791) - 1)) * nspb_var_314(21) + js_var_790
      indf_var_789 = k_indfor_var_780(iplon_var_794, i_lay_var_791)
      z_tauray_var_801 = p_colmol_var_773(iplon_var_794, i_lay_var_791) * rayl_var_256
      DO ig_var_785 = 1, 10
        p_taug_var_782(iplon_var_794, i_lay_var_791, ig_var_785) = z_speccomb_var_798 * ((1.0D0 - z_fs_var_797) * (absb_var_260(ind0_var_786, ig_var_785) * p_fac00_var_763(iplon_var_794, i_lay_var_791) + absb_var_260(ind0_var_786 + 5, ig_var_785) * p_fac10_var_765(iplon_var_794, i_lay_var_791) + absb_var_260(ind1_var_787, ig_var_785) * p_fac01_var_764(iplon_var_794, i_lay_var_791) + absb_var_260(ind1_var_787 + 5, ig_var_785) * p_fac11_var_766(iplon_var_794, i_lay_var_791)) + z_fs_var_797 * (absb_var_260(ind0_var_786 + 1, ig_var_785) * p_fac00_var_763(iplon_var_794, i_lay_var_791) + absb_var_260(ind0_var_786 + 6, ig_var_785) * p_fac10_var_765(iplon_var_794, i_lay_var_791) + absb_var_260(ind1_var_787 + 1, ig_var_785) * p_fac01_var_764(iplon_var_794, i_lay_var_791) + absb_var_260(ind1_var_787 + 6, ig_var_785) * p_fac11_var_766(iplon_var_794, i_lay_var_791))) + p_colh2o_var_771(iplon_var_794, i_lay_var_791) * p_forfac_var_778(iplon_var_794, i_lay_var_791) * (forrefc_var_262(indf_var_789, ig_var_785) + p_forfrac_var_779(iplon_var_794, i_lay_var_791) * (forrefc_var_262(indf_var_789 + 1, ig_var_785) - forrefc_var_262(indf_var_789, ig_var_785)))
        p_taur_var_783(iplon_var_794, i_lay_var_791, ig_var_785) = z_tauray_var_801
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol21
SUBROUTINE srtm_taumol25(kidia_var_802, kfdia_var_803, klev_var_804, p_fac00_var_805, p_fac01_var_806, p_fac10_var_807, p_fac11_var_808, k_jp_var_809, k_jt_var_810, k_jt1_var_811, p_colh2o_var_812, p_colmol_var_813, p_colo3_var_814, k_laytrop_var_815, p_sfluxzen_var_816, p_taug_var_817, p_taur_var_818, prmu0_var_819)
  USE yoesrta25, ONLY: absa_var_288, abso3ac_var_291, abso3bc_var_292, layreffr_var_287, raylc_var_290, sfluxrefc_var_289
  USE yoesrtwn, ONLY: nspa_var_313
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_802, kfdia_var_803
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_804
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_805(kidia_var_802 : kfdia_var_803, klev_var_804)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_806(kidia_var_802 : kfdia_var_803, klev_var_804)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_807(kidia_var_802 : kfdia_var_803, klev_var_804)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_808(kidia_var_802 : kfdia_var_803, klev_var_804)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_809(kidia_var_802 : kfdia_var_803, klev_var_804)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_810(kidia_var_802 : kfdia_var_803, klev_var_804)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_811(kidia_var_802 : kfdia_var_803, klev_var_804)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_812(kidia_var_802 : kfdia_var_803, klev_var_804)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_813(kidia_var_802 : kfdia_var_803, klev_var_804)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_814(kidia_var_802 : kfdia_var_803, klev_var_804)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_815(kidia_var_802 : kfdia_var_803)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_816(kidia_var_802 : kfdia_var_803, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_817(kidia_var_802 : kfdia_var_803, klev_var_804, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_818(kidia_var_802 : kfdia_var_803, klev_var_804, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_819(kidia_var_802 : kfdia_var_803)
  INTEGER(KIND = 4) :: ig_var_820, ind0_var_821, ind1_var_822, i_lay_var_823, i_laysolfr_var_824(kidia_var_802 : kfdia_var_803), i_nlayers_var_825, iplon_var_826
  INTEGER(KIND = 4) :: laytrop_min_var_827, laytrop_max_var_828
  REAL(KIND = 8) :: z_tauray_var_829
  laytrop_min_var_827 = MINVAL(k_laytrop_var_815(kidia_var_802 : kfdia_var_803))
  laytrop_max_var_828 = MAXVAL(k_laytrop_var_815(kidia_var_802 : kfdia_var_803))
  i_nlayers_var_825 = klev_var_804
  DO iplon_var_826 = kidia_var_802, kfdia_var_803
    i_laysolfr_var_824(iplon_var_826) = k_laytrop_var_815(iplon_var_826)
  END DO
  DO i_lay_var_823 = 1, laytrop_min_var_827
    DO iplon_var_826 = kidia_var_802, kfdia_var_803
      IF (k_jp_var_809(iplon_var_826, i_lay_var_823) < layreffr_var_287 .AND. k_jp_var_809(iplon_var_826, i_lay_var_823 + 1) >= layreffr_var_287) i_laysolfr_var_824(iplon_var_826) = MIN(i_lay_var_823 + 1, k_laytrop_var_815(iplon_var_826))
      ind0_var_821 = ((k_jp_var_809(iplon_var_826, i_lay_var_823) - 1) * 5 + (k_jt_var_810(iplon_var_826, i_lay_var_823) - 1)) * nspa_var_313(25) + 1
      ind1_var_822 = (k_jp_var_809(iplon_var_826, i_lay_var_823) * 5 + (k_jt1_var_811(iplon_var_826, i_lay_var_823) - 1)) * nspa_var_313(25) + 1
      DO ig_var_820 = 1, 6
        z_tauray_var_829 = p_colmol_var_813(iplon_var_826, i_lay_var_823) * raylc_var_290(ig_var_820)
        p_taug_var_817(iplon_var_826, i_lay_var_823, ig_var_820) = p_colh2o_var_812(iplon_var_826, i_lay_var_823) * (p_fac00_var_805(iplon_var_826, i_lay_var_823) * absa_var_288(ind0_var_821, ig_var_820) + p_fac10_var_807(iplon_var_826, i_lay_var_823) * absa_var_288(ind0_var_821 + 1, ig_var_820) + p_fac01_var_806(iplon_var_826, i_lay_var_823) * absa_var_288(ind1_var_822, ig_var_820) + p_fac11_var_808(iplon_var_826, i_lay_var_823) * absa_var_288(ind1_var_822 + 1, ig_var_820)) + p_colo3_var_814(iplon_var_826, i_lay_var_823) * abso3ac_var_291(ig_var_820)
        IF (i_lay_var_823 == i_laysolfr_var_824(iplon_var_826)) p_sfluxzen_var_816(iplon_var_826, ig_var_820) = sfluxrefc_var_289(ig_var_820)
        p_taur_var_818(iplon_var_826, i_lay_var_823, ig_var_820) = z_tauray_var_829
      END DO
    END DO
  END DO
  DO i_lay_var_823 = laytrop_min_var_827 + 1, laytrop_max_var_828
    DO iplon_var_826 = kidia_var_802, kfdia_var_803
      IF (i_lay_var_823 <= k_laytrop_var_815(iplon_var_826)) THEN
        IF (k_jp_var_809(iplon_var_826, i_lay_var_823) < layreffr_var_287 .AND. k_jp_var_809(iplon_var_826, i_lay_var_823 + 1) >= layreffr_var_287) i_laysolfr_var_824(iplon_var_826) = MIN(i_lay_var_823 + 1, k_laytrop_var_815(iplon_var_826))
        ind0_var_821 = ((k_jp_var_809(iplon_var_826, i_lay_var_823) - 1) * 5 + (k_jt_var_810(iplon_var_826, i_lay_var_823) - 1)) * nspa_var_313(25) + 1
        ind1_var_822 = (k_jp_var_809(iplon_var_826, i_lay_var_823) * 5 + (k_jt1_var_811(iplon_var_826, i_lay_var_823) - 1)) * nspa_var_313(25) + 1
        DO ig_var_820 = 1, 6
          z_tauray_var_829 = p_colmol_var_813(iplon_var_826, i_lay_var_823) * raylc_var_290(ig_var_820)
          p_taug_var_817(iplon_var_826, i_lay_var_823, ig_var_820) = p_colh2o_var_812(iplon_var_826, i_lay_var_823) * (p_fac00_var_805(iplon_var_826, i_lay_var_823) * absa_var_288(ind0_var_821, ig_var_820) + p_fac10_var_807(iplon_var_826, i_lay_var_823) * absa_var_288(ind0_var_821 + 1, ig_var_820) + p_fac01_var_806(iplon_var_826, i_lay_var_823) * absa_var_288(ind1_var_822, ig_var_820) + p_fac11_var_808(iplon_var_826, i_lay_var_823) * absa_var_288(ind1_var_822 + 1, ig_var_820)) + p_colo3_var_814(iplon_var_826, i_lay_var_823) * abso3ac_var_291(ig_var_820)
          IF (i_lay_var_823 == i_laysolfr_var_824(iplon_var_826)) p_sfluxzen_var_816(iplon_var_826, ig_var_820) = sfluxrefc_var_289(ig_var_820)
          p_taur_var_818(iplon_var_826, i_lay_var_823, ig_var_820) = z_tauray_var_829
        END DO
      ELSE
        DO ig_var_820 = 1, 6
          z_tauray_var_829 = p_colmol_var_813(iplon_var_826, i_lay_var_823) * raylc_var_290(ig_var_820)
          p_taug_var_817(iplon_var_826, i_lay_var_823, ig_var_820) = p_colo3_var_814(iplon_var_826, i_lay_var_823) * abso3bc_var_292(ig_var_820)
          p_taur_var_818(iplon_var_826, i_lay_var_823, ig_var_820) = z_tauray_var_829
        END DO
      END IF
    END DO
  END DO
  DO ig_var_820 = 1, 6
    DO i_lay_var_823 = laytrop_max_var_828 + 1, i_nlayers_var_825
      DO iplon_var_826 = kidia_var_802, kfdia_var_803
        z_tauray_var_829 = p_colmol_var_813(iplon_var_826, i_lay_var_823) * raylc_var_290(ig_var_820)
        p_taug_var_817(iplon_var_826, i_lay_var_823, ig_var_820) = p_colo3_var_814(iplon_var_826, i_lay_var_823) * abso3bc_var_292(ig_var_820)
        p_taur_var_818(iplon_var_826, i_lay_var_823, ig_var_820) = z_tauray_var_829
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol25
SUBROUTINE srtm_taumol27(kidia_var_830, kfdia_var_831, klev_var_832, p_fac00_var_833, p_fac01_var_834, p_fac10_var_835, p_fac11_var_836, k_jp_var_837, k_jt_var_838, k_jt1_var_839, p_colmol_var_840, p_colo3_var_841, k_laytrop_var_842, p_sfluxzen_var_843, p_taug_var_844, p_taur_var_845, prmu0_var_846)
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  USE yoesrta27, ONLY: absa_var_296, absb_var_297, layreffr_var_295, raylc_var_299, scalekur, sfluxrefc_var_298
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_830, kfdia_var_831
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_832
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_833(kidia_var_830 : kfdia_var_831, klev_var_832)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_834(kidia_var_830 : kfdia_var_831, klev_var_832)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_835(kidia_var_830 : kfdia_var_831, klev_var_832)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_836(kidia_var_830 : kfdia_var_831, klev_var_832)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_837(kidia_var_830 : kfdia_var_831, klev_var_832)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_838(kidia_var_830 : kfdia_var_831, klev_var_832)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_839(kidia_var_830 : kfdia_var_831, klev_var_832)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_840(kidia_var_830 : kfdia_var_831, klev_var_832)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_841(kidia_var_830 : kfdia_var_831, klev_var_832)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_842(kidia_var_830 : kfdia_var_831)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_843(kidia_var_830 : kfdia_var_831, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_844(kidia_var_830 : kfdia_var_831, klev_var_832, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_845(kidia_var_830 : kfdia_var_831, klev_var_832, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_846(kidia_var_830 : kfdia_var_831)
  INTEGER(KIND = 4) :: ig_var_847, ind0_var_848, ind1_var_849, i_lay_var_850, i_laysolfr_var_851(kidia_var_830 : kfdia_var_831), i_nlayers_var_852, iplon_var_853
  INTEGER(KIND = 4) :: laytrop_min_var_854, laytrop_max_var_855
  REAL(KIND = 8) :: z_tauray_var_856
  laytrop_min_var_854 = MINVAL(k_laytrop_var_842(kidia_var_830 : kfdia_var_831))
  laytrop_max_var_855 = MAXVAL(k_laytrop_var_842(kidia_var_830 : kfdia_var_831))
  i_nlayers_var_852 = klev_var_832
  DO iplon_var_853 = kidia_var_830, kfdia_var_831
    i_laysolfr_var_851(iplon_var_853) = i_nlayers_var_852
  END DO
  DO i_lay_var_850 = 1, laytrop_min_var_854
    DO iplon_var_853 = kidia_var_830, kfdia_var_831
      ind0_var_848 = ((k_jp_var_837(iplon_var_853, i_lay_var_850) - 1) * 5 + (k_jt_var_838(iplon_var_853, i_lay_var_850) - 1)) * nspa_var_313(27) + 1
      ind1_var_849 = (k_jp_var_837(iplon_var_853, i_lay_var_850) * 5 + (k_jt1_var_839(iplon_var_853, i_lay_var_850) - 1)) * nspa_var_313(27) + 1
      DO ig_var_847 = 1, 8
        z_tauray_var_856 = p_colmol_var_840(iplon_var_853, i_lay_var_850) * raylc_var_299(ig_var_847)
        p_taug_var_844(iplon_var_853, i_lay_var_850, ig_var_847) = p_colo3_var_841(iplon_var_853, i_lay_var_850) * (p_fac00_var_833(iplon_var_853, i_lay_var_850) * absa_var_296(ind0_var_848, ig_var_847) + p_fac10_var_835(iplon_var_853, i_lay_var_850) * absa_var_296(ind0_var_848 + 1, ig_var_847) + p_fac01_var_834(iplon_var_853, i_lay_var_850) * absa_var_296(ind1_var_849, ig_var_847) + p_fac11_var_836(iplon_var_853, i_lay_var_850) * absa_var_296(ind1_var_849 + 1, ig_var_847))
        p_taur_var_845(iplon_var_853, i_lay_var_850, ig_var_847) = z_tauray_var_856
      END DO
    END DO
  END DO
  DO i_lay_var_850 = laytrop_min_var_854 + 1, laytrop_max_var_855
    DO iplon_var_853 = kidia_var_830, kfdia_var_831
      IF (i_lay_var_850 <= k_laytrop_var_842(iplon_var_853)) THEN
        ind0_var_848 = ((k_jp_var_837(iplon_var_853, i_lay_var_850) - 1) * 5 + (k_jt_var_838(iplon_var_853, i_lay_var_850) - 1)) * nspa_var_313(27) + 1
        ind1_var_849 = (k_jp_var_837(iplon_var_853, i_lay_var_850) * 5 + (k_jt1_var_839(iplon_var_853, i_lay_var_850) - 1)) * nspa_var_313(27) + 1
        DO ig_var_847 = 1, 8
          z_tauray_var_856 = p_colmol_var_840(iplon_var_853, i_lay_var_850) * raylc_var_299(ig_var_847)
          p_taug_var_844(iplon_var_853, i_lay_var_850, ig_var_847) = p_colo3_var_841(iplon_var_853, i_lay_var_850) * (p_fac00_var_833(iplon_var_853, i_lay_var_850) * absa_var_296(ind0_var_848, ig_var_847) + p_fac10_var_835(iplon_var_853, i_lay_var_850) * absa_var_296(ind0_var_848 + 1, ig_var_847) + p_fac01_var_834(iplon_var_853, i_lay_var_850) * absa_var_296(ind1_var_849, ig_var_847) + p_fac11_var_836(iplon_var_853, i_lay_var_850) * absa_var_296(ind1_var_849 + 1, ig_var_847))
          p_taur_var_845(iplon_var_853, i_lay_var_850, ig_var_847) = z_tauray_var_856
        END DO
      ELSE
        IF (k_jp_var_837(iplon_var_853, i_lay_var_850 - 1) < layreffr_var_295 .AND. k_jp_var_837(iplon_var_853, i_lay_var_850) >= layreffr_var_295) i_laysolfr_var_851(iplon_var_853) = i_lay_var_850
        ind0_var_848 = ((k_jp_var_837(iplon_var_853, i_lay_var_850) - 13) * 5 + (k_jt_var_838(iplon_var_853, i_lay_var_850) - 1)) * nspb_var_314(27) + 1
        ind1_var_849 = ((k_jp_var_837(iplon_var_853, i_lay_var_850) - 12) * 5 + (k_jt1_var_839(iplon_var_853, i_lay_var_850) - 1)) * nspb_var_314(27) + 1
        DO ig_var_847 = 1, 8
          z_tauray_var_856 = p_colmol_var_840(iplon_var_853, i_lay_var_850) * raylc_var_299(ig_var_847)
          p_taug_var_844(iplon_var_853, i_lay_var_850, ig_var_847) = p_colo3_var_841(iplon_var_853, i_lay_var_850) * (p_fac00_var_833(iplon_var_853, i_lay_var_850) * absb_var_297(ind0_var_848, ig_var_847) + p_fac10_var_835(iplon_var_853, i_lay_var_850) * absb_var_297(ind0_var_848 + 1, ig_var_847) + p_fac01_var_834(iplon_var_853, i_lay_var_850) * absb_var_297(ind1_var_849, ig_var_847) + p_fac11_var_836(iplon_var_853, i_lay_var_850) * absb_var_297(ind1_var_849 + 1, ig_var_847))
          IF (i_lay_var_850 == i_laysolfr_var_851(iplon_var_853)) p_sfluxzen_var_843(iplon_var_853, ig_var_847) = scalekur * sfluxrefc_var_298(ig_var_847)
          p_taur_var_845(iplon_var_853, i_lay_var_850, ig_var_847) = z_tauray_var_856
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_850 = laytrop_max_var_855 + 1, i_nlayers_var_852
    DO iplon_var_853 = kidia_var_830, kfdia_var_831
      IF (k_jp_var_837(iplon_var_853, i_lay_var_850 - 1) < layreffr_var_295 .AND. k_jp_var_837(iplon_var_853, i_lay_var_850) >= layreffr_var_295) i_laysolfr_var_851(iplon_var_853) = i_lay_var_850
      ind0_var_848 = ((k_jp_var_837(iplon_var_853, i_lay_var_850) - 13) * 5 + (k_jt_var_838(iplon_var_853, i_lay_var_850) - 1)) * nspb_var_314(27) + 1
      ind1_var_849 = ((k_jp_var_837(iplon_var_853, i_lay_var_850) - 12) * 5 + (k_jt1_var_839(iplon_var_853, i_lay_var_850) - 1)) * nspb_var_314(27) + 1
      DO ig_var_847 = 1, 8
        z_tauray_var_856 = p_colmol_var_840(iplon_var_853, i_lay_var_850) * raylc_var_299(ig_var_847)
        p_taug_var_844(iplon_var_853, i_lay_var_850, ig_var_847) = p_colo3_var_841(iplon_var_853, i_lay_var_850) * (p_fac00_var_833(iplon_var_853, i_lay_var_850) * absb_var_297(ind0_var_848, ig_var_847) + p_fac10_var_835(iplon_var_853, i_lay_var_850) * absb_var_297(ind0_var_848 + 1, ig_var_847) + p_fac01_var_834(iplon_var_853, i_lay_var_850) * absb_var_297(ind1_var_849, ig_var_847) + p_fac11_var_836(iplon_var_853, i_lay_var_850) * absb_var_297(ind1_var_849 + 1, ig_var_847))
        IF (i_lay_var_850 == i_laysolfr_var_851(iplon_var_853)) p_sfluxzen_var_843(iplon_var_853, ig_var_847) = scalekur * sfluxrefc_var_298(ig_var_847)
        p_taur_var_845(iplon_var_853, i_lay_var_850, ig_var_847) = z_tauray_var_856
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol27
SUBROUTINE rrtm_prepare_gases(kidia_var_858, kfdia_var_859, klon, klev_var_857, paph, pap, pth, pt, pq, pco2, pch4, pn2o, pno2, pc11, pc12, pc22, pcl4, pozn, pcoldry_var_860, pwbrodl, pwkl_var_861, pwx_var_862, pavel_var_863, ptavel_var_864, pz, ptz, kreflect)
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: klon
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_857
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_858, kfdia_var_859
  REAL(KIND = 8), INTENT(IN) :: paph(klon, klev_var_857 + 1)
  REAL(KIND = 8), INTENT(IN) :: pap(klon, klev_var_857)
  REAL(KIND = 8), INTENT(IN) :: pth(klon, klev_var_857 + 1)
  REAL(KIND = 8), INTENT(IN) :: pt(klon, klev_var_857)
  REAL(KIND = 8), INTENT(IN) :: pq(klon, klev_var_857)
  REAL(KIND = 8), INTENT(IN) :: pco2(klon, klev_var_857)
  REAL(KIND = 8), INTENT(IN) :: pch4(klon, klev_var_857)
  REAL(KIND = 8), INTENT(IN) :: pn2o(klon, klev_var_857)
  REAL(KIND = 8), INTENT(IN) :: pno2(klon, klev_var_857)
  REAL(KIND = 8), INTENT(IN) :: pc11(klon, klev_var_857)
  REAL(KIND = 8), INTENT(IN) :: pc12(klon, klev_var_857)
  REAL(KIND = 8), INTENT(IN) :: pc22(klon, klev_var_857)
  REAL(KIND = 8), INTENT(IN) :: pcl4(klon, klev_var_857)
  REAL(KIND = 8), INTENT(IN) :: pozn(klon, klev_var_857)
  REAL(KIND = 8), INTENT(OUT) :: pcoldry_var_860(kidia_var_858 : kfdia_var_859, klev_var_857)
  REAL(KIND = 8), INTENT(OUT) :: pwbrodl(kidia_var_858 : kfdia_var_859, klev_var_857)
  REAL(KIND = 8), INTENT(OUT) :: pwkl_var_861(kidia_var_858 : kfdia_var_859, 35, klev_var_857)
  REAL(KIND = 8), INTENT(OUT) :: pwx_var_862(kidia_var_858 : kfdia_var_859, 4, klev_var_857)
  REAL(KIND = 8), INTENT(OUT) :: pavel_var_863(kidia_var_858 : kfdia_var_859, klev_var_857)
  REAL(KIND = 8), INTENT(OUT) :: ptavel_var_864(kidia_var_858 : kfdia_var_859, klev_var_857)
  REAL(KIND = 8), INTENT(OUT) :: pz(kidia_var_858 : kfdia_var_859, 0 : klev_var_857)
  REAL(KIND = 8), INTENT(OUT) :: ptz(kidia_var_858 : kfdia_var_859, 0 : klev_var_857)
  INTEGER(KIND = 4), INTENT(OUT) :: kreflect(kidia_var_858 : kfdia_var_859)
  REAL(KIND = 8) :: zamd
  REAL(KIND = 8) :: zamw
  REAL(KIND = 8) :: zamco2
  REAL(KIND = 8) :: zamo
  REAL(KIND = 8) :: zamch4
  REAL(KIND = 8) :: zamn2o
  REAL(KIND = 8) :: zamc11
  REAL(KIND = 8) :: zamc12
  REAL(KIND = 8) :: zamc22
  REAL(KIND = 8) :: zamcl4
  REAL(KIND = 8) :: zavgdro
  REAL(KIND = 8) :: zgravit
  REAL(KIND = 8) :: zsummol
  INTEGER(KIND = 4) :: iatm, jmol, ixmax, j1, j2, jk_var_865, jl_var_866
  REAL(KIND = 8) :: zamm
  zamd = 28.97D0
  zamw = 18.0154D0
  zamco2 = 44.011D0
  zamo = 47.9982D0
  zamch4 = 16.043D0
  zamn2o = 44.013D0
  zamc11 = 137.3686D0
  zamc12 = 120.914D0
  zamc22 = 86.469D0
  zamcl4 = 153.823D0
  zavgdro = 6.02214D+23
  zgravit = 980.665D0
  DO jl_var_866 = kidia_var_858, kfdia_var_859
    kreflect(jl_var_866) = 0
  END DO
  DO j2 = 1, klev_var_857
    DO j1 = 1, 35
      DO jl_var_866 = kidia_var_858, kfdia_var_859
        pwkl_var_861(jl_var_866, j1, j2) = 0.0D0
      END DO
    END DO
  END DO
  iatm = 0
  ixmax = 4
  DO jl_var_866 = kidia_var_858, kfdia_var_859
    pz(jl_var_866, 0) = paph(jl_var_866, klev_var_857 + 1) / 100.0D0
    ptz(jl_var_866, 0) = pth(jl_var_866, klev_var_857 + 1)
  END DO
  DO jk_var_865 = 1, klev_var_857
    DO jl_var_866 = kidia_var_858, kfdia_var_859
      pavel_var_863(jl_var_866, jk_var_865) = pap(jl_var_866, klev_var_857 - jk_var_865 + 1) / 100.0D0
      ptavel_var_864(jl_var_866, jk_var_865) = pt(jl_var_866, klev_var_857 - jk_var_865 + 1)
      pz(jl_var_866, jk_var_865) = paph(jl_var_866, klev_var_857 - jk_var_865 + 1) / 100.0D0
      ptz(jl_var_866, jk_var_865) = pth(jl_var_866, klev_var_857 - jk_var_865 + 1)
      pwkl_var_861(jl_var_866, 1, jk_var_865) = MAX(pq(jl_var_866, klev_var_857 - jk_var_865 + 1), 1E-15) * zamd / zamw
      pwkl_var_861(jl_var_866, 2, jk_var_865) = pco2(jl_var_866, klev_var_857 - jk_var_865 + 1) * zamd / zamco2
      pwkl_var_861(jl_var_866, 3, jk_var_865) = pozn(jl_var_866, klev_var_857 - jk_var_865 + 1) * zamd / zamo
      pwkl_var_861(jl_var_866, 4, jk_var_865) = pn2o(jl_var_866, klev_var_857 - jk_var_865 + 1) * zamd / zamn2o
      pwkl_var_861(jl_var_866, 6, jk_var_865) = pch4(jl_var_866, klev_var_857 - jk_var_865 + 1) * zamd / zamch4
      pwkl_var_861(jl_var_866, 7, jk_var_865) = 0.209488D0
      zamm = (1.0D0 - pwkl_var_861(jl_var_866, 1, jk_var_865)) * zamd + pwkl_var_861(jl_var_866, 1, jk_var_865) * zamw
      pcoldry_var_860(jl_var_866, jk_var_865) = (pz(jl_var_866, jk_var_865 - 1) - pz(jl_var_866, jk_var_865)) * 1000.0D0 * zavgdro / (zgravit * zamm * (1.0D0 + pwkl_var_861(jl_var_866, 1, jk_var_865)))
    END DO
  END DO
  DO j2 = 1, klev_var_857
    DO j1 = 1, 4
      DO jl_var_866 = kidia_var_858, kfdia_var_859
        pwx_var_862(jl_var_866, j1, j2) = 0.0D0
      END DO
    END DO
  END DO
  DO jk_var_865 = 1, klev_var_857
    DO jl_var_866 = kidia_var_858, kfdia_var_859
      pwx_var_862(jl_var_866, 1, jk_var_865) = pcl4(jl_var_866, klev_var_857 - jk_var_865 + 1) * zamd / zamcl4
      pwx_var_862(jl_var_866, 2, jk_var_865) = pc11(jl_var_866, klev_var_857 - jk_var_865 + 1) * zamd / zamc11
      pwx_var_862(jl_var_866, 3, jk_var_865) = pc12(jl_var_866, klev_var_857 - jk_var_865 + 1) * zamd / zamc12
      pwx_var_862(jl_var_866, 4, jk_var_865) = pc22(jl_var_866, klev_var_857 - jk_var_865 + 1) * zamd / zamc22
      pwx_var_862(jl_var_866, 1, jk_var_865) = pcoldry_var_860(jl_var_866, jk_var_865) * pwx_var_862(jl_var_866, 1, jk_var_865) * 1D-20
      pwx_var_862(jl_var_866, 2, jk_var_865) = pcoldry_var_860(jl_var_866, jk_var_865) * pwx_var_862(jl_var_866, 2, jk_var_865) * 1D-20
      pwx_var_862(jl_var_866, 3, jk_var_865) = pcoldry_var_860(jl_var_866, jk_var_865) * pwx_var_862(jl_var_866, 3, jk_var_865) * 1D-20
      pwx_var_862(jl_var_866, 4, jk_var_865) = pcoldry_var_860(jl_var_866, jk_var_865) * pwx_var_862(jl_var_866, 4, jk_var_865) * 1D-20
      zsummol = 0.0D0
      DO jmol = 2, 7
        zsummol = zsummol + pwkl_var_861(jl_var_866, jmol, jk_var_865)
      END DO
      pwbrodl(jl_var_866, jk_var_865) = pcoldry_var_860(jl_var_866, jk_var_865) * (1.0D0 - zsummol)
      DO jmol = 1, 7
        pwkl_var_861(jl_var_866, jmol, jk_var_865) = pcoldry_var_860(jl_var_866, jk_var_865) * pwkl_var_861(jl_var_866, jmol, jk_var_865)
      END DO
    END DO
  END DO
END SUBROUTINE rrtm_prepare_gases
SUBROUTINE rrtm_taumol14(kidia_var_867, kfdia_var_868, klev_var_869, taug_var_870, p_tauaerl_var_871, fac00_var_872, fac01_var_873, fac10_var_874, fac11_var_875, forfac_var_886, forfrac_var_887, indfor_var_885, jp_var_876, jt_var_877, jt1_var_878, colco2_var_879, laytrop_var_880, selffac_var_881, selffrac_var_882, indself_var_883, fracs_var_884)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta14, ONLY: absa_var_138, absb_var_139, forref_var_141, fracrefa_var_136, fracrefb_var_137, selfref_var_140
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_867
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_868
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_869
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_870(kidia_var_867 : kfdia_var_868, 140, klev_var_869)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_871(kidia_var_867 : kfdia_var_868, klev_var_869, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_872(kidia_var_867 : kfdia_var_868, klev_var_869)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_873(kidia_var_867 : kfdia_var_868, klev_var_869)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_874(kidia_var_867 : kfdia_var_868, klev_var_869)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_875(kidia_var_867 : kfdia_var_868, klev_var_869)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_876(kidia_var_867 : kfdia_var_868, klev_var_869)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_877(kidia_var_867 : kfdia_var_868, klev_var_869)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_878(kidia_var_867 : kfdia_var_868, klev_var_869)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_879(kidia_var_867 : kfdia_var_868, klev_var_869)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_880(kidia_var_867 : kfdia_var_868)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_881(kidia_var_867 : kfdia_var_868, klev_var_869)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_882(kidia_var_867 : kfdia_var_868, klev_var_869)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_883(kidia_var_867 : kfdia_var_868, klev_var_869)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_884(kidia_var_867 : kfdia_var_868, 140, klev_var_869)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_885(kidia_var_867 : kfdia_var_868, klev_var_869)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_886(kidia_var_867 : kfdia_var_868, klev_var_869)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_887(kidia_var_867 : kfdia_var_868, klev_var_869)
  INTEGER(KIND = 4) :: ig_var_888, ind0_var_889, ind1_var_890, inds_var_891, indf_var_892, lay_var_893
  REAL(KIND = 8) :: taufor_var_894, tauself_var_895
  INTEGER(KIND = 4) :: laytrop_min_var_896, laytrop_max_var_897
  INTEGER(KIND = 4) :: ixc_var_898(klev_var_869), ixlow_var_899(kfdia_var_868, klev_var_869), ixhigh_var_900(kfdia_var_868, klev_var_869)
  INTEGER(KIND = 4) :: ich_var_901, icl_var_902, ixc0_var_903, ixp_var_904, jc_var_905, jl_var_906
  laytrop_min_var_896 = MINVAL(laytrop_var_880)
  laytrop_max_var_897 = MAXVAL(laytrop_var_880)
  ixlow_var_899 = 0
  ixhigh_var_900 = 0
  ixc_var_898 = 0
  DO lay_var_893 = laytrop_min_var_896 + 1, laytrop_max_var_897
    icl_var_902 = 0
    ich_var_901 = 0
    DO jc_var_905 = kidia_var_867, kfdia_var_868
      IF (lay_var_893 <= laytrop_var_880(jc_var_905)) THEN
        icl_var_902 = icl_var_902 + 1
        ixlow_var_899(icl_var_902, lay_var_893) = jc_var_905
      ELSE
        ich_var_901 = ich_var_901 + 1
        ixhigh_var_900(ich_var_901, lay_var_893) = jc_var_905
      END IF
    END DO
    ixc_var_898(lay_var_893) = icl_var_902
  END DO
  DO lay_var_893 = 1, laytrop_min_var_896
    DO jl_var_906 = kidia_var_867, kfdia_var_868
      ind0_var_889 = ((jp_var_876(jl_var_906, lay_var_893) - 1) * 5 + (jt_var_877(jl_var_906, lay_var_893) - 1)) * nspa_var_216(14) + 1
      ind1_var_890 = (jp_var_876(jl_var_906, lay_var_893) * 5 + (jt1_var_878(jl_var_906, lay_var_893) - 1)) * nspa_var_216(14) + 1
      inds_var_891 = indself_var_883(jl_var_906, lay_var_893)
      indf_var_892 = indfor_var_885(jl_var_906, lay_var_893)
      DO ig_var_888 = 1, 2
        tauself_var_895 = selffac_var_881(jl_var_906, lay_var_893) * (selfref_var_140(inds_var_891, ig_var_888) + selffrac_var_882(jl_var_906, lay_var_893) * (selfref_var_140(inds_var_891 + 1, ig_var_888) - selfref_var_140(inds_var_891, ig_var_888)))
        taufor_var_894 = forfac_var_886(jl_var_906, lay_var_893) * (forref_var_141(indf_var_892, ig_var_888) + forfrac_var_887(jl_var_906, lay_var_893) * (forref_var_141(indf_var_892 + 1, ig_var_888) - forref_var_141(indf_var_892, ig_var_888)))
        taug_var_870(jl_var_906, 134 + ig_var_888, lay_var_893) = colco2_var_879(jl_var_906, lay_var_893) * (fac00_var_872(jl_var_906, lay_var_893) * absa_var_138(ind0_var_889, ig_var_888) + fac10_var_874(jl_var_906, lay_var_893) * absa_var_138(ind0_var_889 + 1, ig_var_888) + fac01_var_873(jl_var_906, lay_var_893) * absa_var_138(ind1_var_890, ig_var_888) + fac11_var_875(jl_var_906, lay_var_893) * absa_var_138(ind1_var_890 + 1, ig_var_888)) + tauself_var_895 + taufor_var_894
        fracs_var_884(jl_var_906, 134 + ig_var_888, lay_var_893) = fracrefa_var_136(ig_var_888)
      END DO
    END DO
  END DO
  DO lay_var_893 = laytrop_max_var_897 + 1, klev_var_869
    DO jl_var_906 = kidia_var_867, kfdia_var_868
      ind0_var_889 = ((jp_var_876(jl_var_906, lay_var_893) - 13) * 5 + (jt_var_877(jl_var_906, lay_var_893) - 1)) * nspb_var_217(14) + 1
      ind1_var_890 = ((jp_var_876(jl_var_906, lay_var_893) - 12) * 5 + (jt1_var_878(jl_var_906, lay_var_893) - 1)) * nspb_var_217(14) + 1
      DO ig_var_888 = 1, 2
        taug_var_870(jl_var_906, 134 + ig_var_888, lay_var_893) = colco2_var_879(jl_var_906, lay_var_893) * (fac00_var_872(jl_var_906, lay_var_893) * absb_var_139(ind0_var_889, ig_var_888) + fac10_var_874(jl_var_906, lay_var_893) * absb_var_139(ind0_var_889 + 1, ig_var_888) + fac01_var_873(jl_var_906, lay_var_893) * absb_var_139(ind1_var_890, ig_var_888) + fac11_var_875(jl_var_906, lay_var_893) * absb_var_139(ind1_var_890 + 1, ig_var_888))
        fracs_var_884(jl_var_906, 134 + ig_var_888, lay_var_893) = fracrefb_var_137(ig_var_888)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_897 /= laytrop_min_var_896) THEN
    DO lay_var_893 = laytrop_min_var_896 + 1, laytrop_max_var_897
      ixc0_var_903 = ixc_var_898(lay_var_893)
      DO ixp_var_904 = 1, ixc0_var_903
        jl_var_906 = ixlow_var_899(ixp_var_904, lay_var_893)
        ind0_var_889 = ((jp_var_876(jl_var_906, lay_var_893) - 1) * 5 + (jt_var_877(jl_var_906, lay_var_893) - 1)) * nspa_var_216(14) + 1
        ind1_var_890 = (jp_var_876(jl_var_906, lay_var_893) * 5 + (jt1_var_878(jl_var_906, lay_var_893) - 1)) * nspa_var_216(14) + 1
        inds_var_891 = indself_var_883(jl_var_906, lay_var_893)
        indf_var_892 = indfor_var_885(jl_var_906, lay_var_893)
        DO ig_var_888 = 1, 2
          tauself_var_895 = selffac_var_881(jl_var_906, lay_var_893) * (selfref_var_140(inds_var_891, ig_var_888) + selffrac_var_882(jl_var_906, lay_var_893) * (selfref_var_140(inds_var_891 + 1, ig_var_888) - selfref_var_140(inds_var_891, ig_var_888)))
          taufor_var_894 = forfac_var_886(jl_var_906, lay_var_893) * (forref_var_141(indf_var_892, ig_var_888) + forfrac_var_887(jl_var_906, lay_var_893) * (forref_var_141(indf_var_892 + 1, ig_var_888) - forref_var_141(indf_var_892, ig_var_888)))
          taug_var_870(jl_var_906, 134 + ig_var_888, lay_var_893) = colco2_var_879(jl_var_906, lay_var_893) * (fac00_var_872(jl_var_906, lay_var_893) * absa_var_138(ind0_var_889, ig_var_888) + fac10_var_874(jl_var_906, lay_var_893) * absa_var_138(ind0_var_889 + 1, ig_var_888) + fac01_var_873(jl_var_906, lay_var_893) * absa_var_138(ind1_var_890, ig_var_888) + fac11_var_875(jl_var_906, lay_var_893) * absa_var_138(ind1_var_890 + 1, ig_var_888)) + tauself_var_895 + taufor_var_894
          fracs_var_884(jl_var_906, 134 + ig_var_888, lay_var_893) = fracrefa_var_136(ig_var_888)
        END DO
      END DO
      ixc0_var_903 = kfdia_var_868 - kidia_var_867 + 1 - ixc0_var_903
      DO ixp_var_904 = 1, ixc0_var_903
        jl_var_906 = ixhigh_var_900(ixp_var_904, lay_var_893)
        ind0_var_889 = ((jp_var_876(jl_var_906, lay_var_893) - 13) * 5 + (jt_var_877(jl_var_906, lay_var_893) - 1)) * nspb_var_217(14) + 1
        ind1_var_890 = ((jp_var_876(jl_var_906, lay_var_893) - 12) * 5 + (jt1_var_878(jl_var_906, lay_var_893) - 1)) * nspb_var_217(14) + 1
        DO ig_var_888 = 1, 2
          taug_var_870(jl_var_906, 134 + ig_var_888, lay_var_893) = colco2_var_879(jl_var_906, lay_var_893) * (fac00_var_872(jl_var_906, lay_var_893) * absb_var_139(ind0_var_889, ig_var_888) + fac10_var_874(jl_var_906, lay_var_893) * absb_var_139(ind0_var_889 + 1, ig_var_888) + fac01_var_873(jl_var_906, lay_var_893) * absb_var_139(ind1_var_890, ig_var_888) + fac11_var_875(jl_var_906, lay_var_893) * absb_var_139(ind1_var_890 + 1, ig_var_888))
          fracs_var_884(jl_var_906, 134 + ig_var_888, lay_var_893) = fracrefb_var_137(ig_var_888)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol14
SUBROUTINE rrtm_gas_optical_depth(kidia_var_907, kfdia_var_908, klev_var_909, pod_var_910, pavel_var_911, pcoldry_var_912, pcolbrd, pwx_var_913, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, poneminus_var_921, pcolh2o_var_922, pcolco2_var_923, pcolo3_var_924, pcoln2o, pcolch4_var_925, pcolo2_var_926, p_co2mult_var_927, klaytrop_var_928, klayswtch, klaylow, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, kindminor, pscaleminor, pscaleminorn2, pminorfrac, prat_h2oco2_var_936, prat_h2oco2_1_var_937, prat_h2oo3_var_938, prat_h2oo3_1_var_939, prat_h2on2o_var_940, prat_h2on2o_1_var_941, prat_h2och4_var_942, prat_h2och4_1_var_943, prat_n2oco2_var_944, prat_n2oco2_1_var_945, prat_o3co2_var_946, prat_o3co2_1_var_947)
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_907
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_908
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_909
  REAL(KIND = 8), INTENT(OUT) :: pod_var_910(140, klev_var_909, kidia_var_907 : kfdia_var_908)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_911(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pcoldry_var_912(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pwx_var_913(kidia_var_907 : kfdia_var_908, 4, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: ptauaerl(kidia_var_907 : kfdia_var_908, klev_var_909, 16)
  REAL(KIND = 8), INTENT(IN) :: pfac00_var_914(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pfac01_var_915(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pfac10_var_916(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pfac11_var_917(kidia_var_907 : kfdia_var_908, klev_var_909)
  INTEGER(KIND = 4), INTENT(IN) :: kjp_var_918(kidia_var_907 : kfdia_var_908, klev_var_909)
  INTEGER(KIND = 4), INTENT(IN) :: kjt_var_919(kidia_var_907 : kfdia_var_908, klev_var_909)
  INTEGER(KIND = 4), INTENT(IN) :: kjt1_var_920(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: poneminus_var_921
  REAL(KIND = 8), INTENT(IN) :: pcolh2o_var_922(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pcolco2_var_923(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pcolo3_var_924(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pcoln2o(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pcolch4_var_925(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pcolo2_var_926(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: p_co2mult_var_927(kidia_var_907 : kfdia_var_908, klev_var_909)
  INTEGER(KIND = 4), INTENT(IN) :: klaytrop_var_928(kidia_var_907 : kfdia_var_908)
  INTEGER(KIND = 4), INTENT(IN) :: klayswtch(kidia_var_907 : kfdia_var_908)
  INTEGER(KIND = 4), INTENT(IN) :: klaylow(kidia_var_907 : kfdia_var_908)
  REAL(KIND = 8), INTENT(IN) :: pselffac_var_929(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pselffrac_var_930(kidia_var_907 : kfdia_var_908, klev_var_909)
  INTEGER(KIND = 4), INTENT(IN) :: kindself_var_931(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(OUT) :: pfrac_var_932(kidia_var_907 : kfdia_var_908, 140, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pforfac_var_933(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pforfrac_var_934(kidia_var_907 : kfdia_var_908, klev_var_909)
  INTEGER(KIND = 4), INTENT(IN) :: kindfor_var_935(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pminorfrac(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pscaleminor(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pscaleminorn2(kidia_var_907 : kfdia_var_908, klev_var_909)
  INTEGER(KIND = 4), INTENT(IN) :: kindminor(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: pcolbrd(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8), INTENT(IN) :: prat_h2oco2_var_936(kidia_var_907 : kfdia_var_908, klev_var_909), prat_h2oco2_1_var_937(kidia_var_907 : kfdia_var_908, klev_var_909), prat_h2oo3_var_938(kidia_var_907 : kfdia_var_908, klev_var_909), prat_h2oo3_1_var_939(kidia_var_907 : kfdia_var_908, klev_var_909), prat_h2on2o_var_940(kidia_var_907 : kfdia_var_908, klev_var_909), prat_h2on2o_1_var_941(kidia_var_907 : kfdia_var_908, klev_var_909), prat_h2och4_var_942(kidia_var_907 : kfdia_var_908, klev_var_909), prat_h2och4_1_var_943(kidia_var_907 : kfdia_var_908, klev_var_909), prat_n2oco2_var_944(kidia_var_907 : kfdia_var_908, klev_var_909), prat_n2oco2_1_var_945(kidia_var_907 : kfdia_var_908, klev_var_909), prat_o3co2_var_946(kidia_var_907 : kfdia_var_908, klev_var_909), prat_o3co2_1_var_947(kidia_var_907 : kfdia_var_908, klev_var_909)
  REAL(KIND = 8) :: ztau(kidia_var_907 : kfdia_var_908, 140, klev_var_909)
  INTEGER(KIND = 4) :: ji, jlev_var_948
  INTEGER(KIND = 4) :: jlon_var_949
  pfrac_var_932(:, :, :) = 0.0D0
  CALL rrtm_taumol1(kidia_var_907, kfdia_var_908, klev_var_909, ztau, pavel_var_911, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, pcolh2o_var_922, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, pminorfrac, kindminor, pscaleminorn2, pcolbrd)
  CALL rrtm_taumol2(kidia_var_907, kfdia_var_908, klev_var_909, ztau, pavel_var_911, pcoldry_var_912, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, pcolh2o_var_922, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932)
  CALL rrtm_taumol3(kidia_var_907, kfdia_var_908, klev_var_909, ztau, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, poneminus_var_921, pcolh2o_var_922, pcolco2_var_923, pcoln2o, pcoldry_var_912, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, prat_h2oco2_var_936, prat_h2oco2_1_var_937, pminorfrac, kindminor)
  CALL rrtm_taumol4(kidia_var_907, kfdia_var_908, klev_var_909, ztau, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, poneminus_var_921, pcolh2o_var_922, pcolco2_var_923, pcolo3_var_924, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, prat_h2oco2_var_936, prat_h2oco2_1_var_937, prat_o3co2_var_946, prat_o3co2_1_var_947)
  CALL rrtm_taumol5(kidia_var_907, kfdia_var_908, klev_var_909, ztau, pwx_var_913, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, poneminus_var_921, pcolh2o_var_922, pcolco2_var_923, pcolo3_var_924, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, prat_h2oco2_var_936, prat_h2oco2_1_var_937, prat_o3co2_var_946, prat_o3co2_1_var_947, pminorfrac, kindminor)
  CALL rrtm_taumol6(kidia_var_907, kfdia_var_908, klev_var_909, ztau, pwx_var_913, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, pcolh2o_var_922, pcolco2_var_923, pcoldry_var_912, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, pminorfrac, kindminor)
  CALL rrtm_taumol7(kidia_var_907, kfdia_var_908, klev_var_909, ztau, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, poneminus_var_921, pcolh2o_var_922, pcolo3_var_924, pcolco2_var_923, pcoldry_var_912, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, prat_h2oo3_var_938, prat_h2oo3_1_var_939, pminorfrac, kindminor)
  CALL rrtm_taumol8(kidia_var_907, kfdia_var_908, klev_var_909, ztau, pwx_var_913, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, pcolh2o_var_922, pcolo3_var_924, pcoln2o, pcolco2_var_923, pcoldry_var_912, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, pminorfrac, kindminor)
  CALL rrtm_taumol9(kidia_var_907, kfdia_var_908, klev_var_909, ztau, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, poneminus_var_921, pcolh2o_var_922, pcoln2o, pcolch4_var_925, pcoldry_var_912, klaytrop_var_928, klayswtch, klaylow, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, prat_h2och4_var_942, prat_h2och4_1_var_943, pminorfrac, kindminor)
  CALL rrtm_taumol10(kidia_var_907, kfdia_var_908, klev_var_909, ztau, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, pcolh2o_var_922, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932)
  CALL rrtm_taumol11(kidia_var_907, kfdia_var_908, klev_var_909, ztau, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, pcolh2o_var_922, pcolo2_var_926, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, pminorfrac, kindminor, pscaleminor)
  CALL rrtm_taumol12(kidia_var_907, kfdia_var_908, klev_var_909, ztau, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, poneminus_var_921, pcolh2o_var_922, pcolco2_var_923, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, prat_h2oco2_var_936, prat_h2oco2_1_var_937)
  CALL rrtm_taumol13(kidia_var_907, kfdia_var_908, klev_var_909, ztau, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, poneminus_var_921, pcolh2o_var_922, pcoln2o, pcolco2_var_923, pcolo3_var_924, pcoldry_var_912, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, prat_h2on2o_var_940, prat_h2on2o_1_var_941, pminorfrac, kindminor)
  CALL rrtm_taumol14(kidia_var_907, kfdia_var_908, klev_var_909, ztau, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, pcolco2_var_923, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932)
  CALL rrtm_taumol15(kidia_var_907, kfdia_var_908, klev_var_909, ztau, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, poneminus_var_921, pcolh2o_var_922, pcolco2_var_923, pcoln2o, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, prat_n2oco2_var_944, prat_n2oco2_1_var_945, pminorfrac, kindminor, pscaleminor, pcolbrd)
  CALL rrtm_taumol16(kidia_var_907, kfdia_var_908, klev_var_909, ztau, ptauaerl, pfac00_var_914, pfac01_var_915, pfac10_var_916, pfac11_var_917, pforfac_var_933, pforfrac_var_934, kindfor_var_935, kjp_var_918, kjt_var_919, kjt1_var_920, poneminus_var_921, pcolh2o_var_922, pcolch4_var_925, klaytrop_var_928, pselffac_var_929, pselffrac_var_930, kindself_var_931, pfrac_var_932, prat_h2och4_var_942, prat_h2och4_1_var_943)
  DO jlev_var_948 = 1, klev_var_909
    DO ji = 1, 140
      DO jlon_var_949 = kidia_var_907, kfdia_var_908
        pod_var_910(ji, jlev_var_948, jlon_var_949) = ztau(jlon_var_949, ji, jlev_var_948)
      END DO
    END DO
  END DO
END SUBROUTINE rrtm_gas_optical_depth
SUBROUTINE rrtm_taumol3(kidia_var_950, kfdia_var_951, klev_var_952, taug_var_953, p_tauaerl_var_954, fac00_var_955, fac01_var_956, fac10_var_957, fac11_var_958, forfac_var_959, forfrac_var_976, indfor_var_975, jp_var_960, jt_var_961, jt1_var_962, oneminus_var_963, colh2o_var_964, colco2_var_965, coln2o_var_966, coldry_var_967, laytrop_var_968, selffac_var_969, selffrac_var_970, indself_var_971, fracs_var_972, rat_h2oco2_var_973, rat_h2oco2_1_var_974, minorfrac_var_977, indminor_var_978)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng3
  USE yoerrta3, ONLY: absa_var_163, absb_var_164, forref_var_166, fracrefa_var_159, fracrefb_var_160, ka_mn2o_var_161, kb_mn2o_var_162, selfref_var_165
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_950
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_951
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_952
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_953(kidia_var_950 : kfdia_var_951, 140, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_954(kidia_var_950 : kfdia_var_951, klev_var_952, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_955(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_956(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_957(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_958(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_959(kidia_var_950 : kfdia_var_951, klev_var_952)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_960(kidia_var_950 : kfdia_var_951, klev_var_952)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_961(kidia_var_950 : kfdia_var_951, klev_var_952)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_962(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_963
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_964(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_965(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_966(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_967(kidia_var_950 : kfdia_var_951, klev_var_952)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_968(kidia_var_950 : kfdia_var_951)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_969(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_970(kidia_var_950 : kfdia_var_951, klev_var_952)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_971(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_972(kidia_var_950 : kfdia_var_951, 140, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_973(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_974(kidia_var_950 : kfdia_var_951, klev_var_952)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_975(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_976(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_977(kidia_var_950 : kfdia_var_951, klev_var_952)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_978(kidia_var_950 : kfdia_var_951, klev_var_952)
  REAL(KIND = 8) :: speccomb_var_979, speccomb1_var_980, speccomb_mn2o_var_981, speccomb_planck_var_982
  REAL(KIND = 8) :: refrat_planck_a_var_983, refrat_planck_b_var_984, refrat_m_a_var_985, refrat_m_b
  INTEGER(KIND = 4) :: ind0_var_986, ind1_var_987, inds_var_988, indf_var_989, indm_var_990
  INTEGER(KIND = 4) :: ig_var_991, js_var_992, lay_var_993, js1_var_994, jmn2o_var_995, jpl_var_996
  REAL(KIND = 8) :: fs_var_997, specmult_var_998, specparm_var_999, fs1_var_1000, specmult1_var_1001, specparm1_var_1002, fmn2o_var_1003, specmult_mn2o_var_1004, specparm_mn2o_var_1005, fpl_var_1006, specmult_planck_var_1007, specparm_planck_var_1008
  REAL(KIND = 8) :: adjfac_var_1009, adjcoln2o_var_1010, ratn2o_var_1011, chi_n2o_var_1012
  REAL(KIND = 8) :: fac000_var_1013, fac100_var_1014, fac200_var_1015, fac010_var_1016, fac110_var_1017, fac210_var_1018, fac001_var_1019, fac101_var_1020, fac201_var_1021, fac011_var_1022, fac111_var_1023, fac211_var_1024
  REAL(KIND = 8) :: p_var_1025, p4_var_1026, fk0_var_1027, fk1_var_1028, fk2_var_1029
  REAL(KIND = 8) :: taufor_var_1030, tauself_var_1031, n2om1_var_1032, n2om2_var_1033, absn2o_var_1034, tau_major_var_1035(16), tau_major1_var_1036(16)
  INTEGER(KIND = 4) :: laytrop_min_var_1037, laytrop_max_var_1038
  INTEGER(KIND = 4) :: ixc_var_1039(klev_var_952), ixlow_var_1040(kfdia_var_951, klev_var_952), ixhigh_var_1041(kfdia_var_951, klev_var_952)
  INTEGER(KIND = 4) :: ich_var_1042, icl_var_1043, ixc0_var_1044, ixp_var_1045, jc_var_1046, jl_var_1047
  laytrop_min_var_1037 = MINVAL(laytrop_var_968)
  laytrop_max_var_1038 = MAXVAL(laytrop_var_968)
  ixlow_var_1040 = 0
  ixhigh_var_1041 = 0
  ixc_var_1039 = 0
  DO lay_var_993 = laytrop_min_var_1037 + 1, laytrop_max_var_1038
    icl_var_1043 = 0
    ich_var_1042 = 0
    DO jc_var_1046 = kidia_var_950, kfdia_var_951
      IF (lay_var_993 <= laytrop_var_968(jc_var_1046)) THEN
        icl_var_1043 = icl_var_1043 + 1
        ixlow_var_1040(icl_var_1043, lay_var_993) = jc_var_1046
      ELSE
        ich_var_1042 = ich_var_1042 + 1
        ixhigh_var_1041(ich_var_1042, lay_var_993) = jc_var_1046
      END IF
    END DO
    ixc_var_1039(lay_var_993) = icl_var_1043
  END DO
  refrat_planck_a_var_983 = chi_mls(1, 9) / chi_mls(2, 9)
  refrat_planck_b_var_984 = chi_mls(1, 13) / chi_mls(2, 13)
  refrat_m_a_var_985 = chi_mls(1, 3) / chi_mls(2, 3)
  refrat_m_b = chi_mls(1, 13) / chi_mls(2, 13)
  DO lay_var_993 = 1, laytrop_min_var_1037
    DO jl_var_1047 = kidia_var_950, kfdia_var_951
      speccomb_var_979 = colh2o_var_964(jl_var_1047, lay_var_993) + rat_h2oco2_var_973(jl_var_1047, lay_var_993) * colco2_var_965(jl_var_1047, lay_var_993)
      specparm_var_999 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_var_979, oneminus_var_963)
      specmult_var_998 = 8.0D0 * (specparm_var_999)
      js_var_992 = 1 + INT(specmult_var_998)
      fs_var_997 = ((specmult_var_998) - AINT((specmult_var_998)))
      speccomb1_var_980 = colh2o_var_964(jl_var_1047, lay_var_993) + rat_h2oco2_1_var_974(jl_var_1047, lay_var_993) * colco2_var_965(jl_var_1047, lay_var_993)
      specparm1_var_1002 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb1_var_980, oneminus_var_963)
      specmult1_var_1001 = 8.0D0 * (specparm1_var_1002)
      js1_var_994 = 1 + INT(specmult1_var_1001)
      fs1_var_1000 = ((specmult1_var_1001) - AINT((specmult1_var_1001)))
      speccomb_mn2o_var_981 = colh2o_var_964(jl_var_1047, lay_var_993) + refrat_m_a_var_985 * colco2_var_965(jl_var_1047, lay_var_993)
      specparm_mn2o_var_1005 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_mn2o_var_981, oneminus_var_963)
      specmult_mn2o_var_1004 = 8.0D0 * specparm_mn2o_var_1005
      jmn2o_var_995 = 1 + INT(specmult_mn2o_var_1004)
      fmn2o_var_1003 = ((specmult_mn2o_var_1004) - AINT((specmult_mn2o_var_1004)))
      chi_n2o_var_1012 = coln2o_var_966(jl_var_1047, lay_var_993) / coldry_var_967(jl_var_1047, lay_var_993)
      ratn2o_var_1011 = 1D+20 * chi_n2o_var_1012 / chi_mls(4, jp_var_960(jl_var_1047, lay_var_993) + 1)
      IF (ratn2o_var_1011 .GT. 1.5D0) THEN
        adjfac_var_1009 = 0.5D0 + (ratn2o_var_1011 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1010 = adjfac_var_1009 * chi_mls(4, jp_var_960(jl_var_1047, lay_var_993) + 1) * coldry_var_967(jl_var_1047, lay_var_993) * 1D-20
      ELSE
        adjcoln2o_var_1010 = coln2o_var_966(jl_var_1047, lay_var_993)
      END IF
      speccomb_planck_var_982 = colh2o_var_964(jl_var_1047, lay_var_993) + refrat_planck_a_var_983 * colco2_var_965(jl_var_1047, lay_var_993)
      specparm_planck_var_1008 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_planck_var_982, oneminus_var_963)
      specmult_planck_var_1007 = 8.0D0 * specparm_planck_var_1008
      jpl_var_996 = 1 + INT(specmult_planck_var_1007)
      fpl_var_1006 = ((specmult_planck_var_1007) - AINT((specmult_planck_var_1007)))
      ind0_var_986 = ((jp_var_960(jl_var_1047, lay_var_993) - 1) * 5 + (jt_var_961(jl_var_1047, lay_var_993) - 1)) * nspa_var_216(3) + js_var_992
      ind1_var_987 = (jp_var_960(jl_var_1047, lay_var_993) * 5 + (jt1_var_962(jl_var_1047, lay_var_993) - 1)) * nspa_var_216(3) + js1_var_994
      inds_var_988 = indself_var_971(jl_var_1047, lay_var_993)
      indf_var_989 = indfor_var_975(jl_var_1047, lay_var_993)
      indm_var_990 = indminor_var_978(jl_var_1047, lay_var_993)
      IF (specparm_var_999 .LT. 0.125D0) THEN
        p_var_1025 = fs_var_997 - 1.0D0
        p4_var_1026 = p_var_1025 ** 4
        fk0_var_1027 = p4_var_1026
        fk1_var_1028 = 1.0D0 - p_var_1025 - 2.0D0 * p4_var_1026
        fk2_var_1029 = p_var_1025 + p4_var_1026
        fac000_var_1013 = fk0_var_1027 * fac00_var_955(jl_var_1047, lay_var_993)
        fac100_var_1014 = fk1_var_1028 * fac00_var_955(jl_var_1047, lay_var_993)
        fac200_var_1015 = fk2_var_1029 * fac00_var_955(jl_var_1047, lay_var_993)
        fac010_var_1016 = fk0_var_1027 * fac10_var_957(jl_var_1047, lay_var_993)
        fac110_var_1017 = fk1_var_1028 * fac10_var_957(jl_var_1047, lay_var_993)
        fac210_var_1018 = fk2_var_1029 * fac10_var_957(jl_var_1047, lay_var_993)
      ELSE IF (specparm_var_999 .GT. 0.875D0) THEN
        p_var_1025 = - fs_var_997
        p4_var_1026 = p_var_1025 ** 4
        fk0_var_1027 = p4_var_1026
        fk1_var_1028 = 1.0D0 - p_var_1025 - 2.0D0 * p4_var_1026
        fk2_var_1029 = p_var_1025 + p4_var_1026
        fac000_var_1013 = fk0_var_1027 * fac00_var_955(jl_var_1047, lay_var_993)
        fac100_var_1014 = fk1_var_1028 * fac00_var_955(jl_var_1047, lay_var_993)
        fac200_var_1015 = fk2_var_1029 * fac00_var_955(jl_var_1047, lay_var_993)
        fac010_var_1016 = fk0_var_1027 * fac10_var_957(jl_var_1047, lay_var_993)
        fac110_var_1017 = fk1_var_1028 * fac10_var_957(jl_var_1047, lay_var_993)
        fac210_var_1018 = fk2_var_1029 * fac10_var_957(jl_var_1047, lay_var_993)
      ELSE
        fac000_var_1013 = (1.0D0 - fs_var_997) * fac00_var_955(jl_var_1047, lay_var_993)
        fac010_var_1016 = (1.0D0 - fs_var_997) * fac10_var_957(jl_var_1047, lay_var_993)
        fac100_var_1014 = fs_var_997 * fac00_var_955(jl_var_1047, lay_var_993)
        fac110_var_1017 = fs_var_997 * fac10_var_957(jl_var_1047, lay_var_993)
        fac200_var_1015 = 0.0D0
        fac210_var_1018 = 0.0D0
      END IF
      IF (specparm1_var_1002 .LT. 0.125D0) THEN
        p_var_1025 = fs1_var_1000 - 1.0D0
        p4_var_1026 = p_var_1025 ** 4
        fk0_var_1027 = p4_var_1026
        fk1_var_1028 = 1.0D0 - p_var_1025 - 2.0D0 * p4_var_1026
        fk2_var_1029 = p_var_1025 + p4_var_1026
        fac001_var_1019 = fk0_var_1027 * fac01_var_956(jl_var_1047, lay_var_993)
        fac101_var_1020 = fk1_var_1028 * fac01_var_956(jl_var_1047, lay_var_993)
        fac201_var_1021 = fk2_var_1029 * fac01_var_956(jl_var_1047, lay_var_993)
        fac011_var_1022 = fk0_var_1027 * fac11_var_958(jl_var_1047, lay_var_993)
        fac111_var_1023 = fk1_var_1028 * fac11_var_958(jl_var_1047, lay_var_993)
        fac211_var_1024 = fk2_var_1029 * fac11_var_958(jl_var_1047, lay_var_993)
      ELSE IF (specparm1_var_1002 .GT. 0.875D0) THEN
        p_var_1025 = - fs1_var_1000
        p4_var_1026 = p_var_1025 ** 4
        fk0_var_1027 = p4_var_1026
        fk1_var_1028 = 1.0D0 - p_var_1025 - 2.0D0 * p4_var_1026
        fk2_var_1029 = p_var_1025 + p4_var_1026
        fac001_var_1019 = fk0_var_1027 * fac01_var_956(jl_var_1047, lay_var_993)
        fac101_var_1020 = fk1_var_1028 * fac01_var_956(jl_var_1047, lay_var_993)
        fac201_var_1021 = fk2_var_1029 * fac01_var_956(jl_var_1047, lay_var_993)
        fac011_var_1022 = fk0_var_1027 * fac11_var_958(jl_var_1047, lay_var_993)
        fac111_var_1023 = fk1_var_1028 * fac11_var_958(jl_var_1047, lay_var_993)
        fac211_var_1024 = fk2_var_1029 * fac11_var_958(jl_var_1047, lay_var_993)
      ELSE
        fac001_var_1019 = (1.0D0 - fs1_var_1000) * fac01_var_956(jl_var_1047, lay_var_993)
        fac011_var_1022 = (1.0D0 - fs1_var_1000) * fac11_var_958(jl_var_1047, lay_var_993)
        fac101_var_1020 = fs1_var_1000 * fac01_var_956(jl_var_1047, lay_var_993)
        fac111_var_1023 = fs1_var_1000 * fac11_var_958(jl_var_1047, lay_var_993)
        fac201_var_1021 = 0.0D0
        fac211_var_1024 = 0.0D0
      END IF
      IF (specparm_var_999 .LT. 0.125D0) THEN
        tau_major_var_1035(1 : ng3) = speccomb_var_979 * (fac000_var_1013 * absa_var_163(ind0_var_986, 1 : 16) + fac100_var_1014 * absa_var_163(ind0_var_986 + 1, 1 : 16) + fac200_var_1015 * absa_var_163(ind0_var_986 + 2, 1 : 16) + fac010_var_1016 * absa_var_163(ind0_var_986 + 9, 1 : 16) + fac110_var_1017 * absa_var_163(ind0_var_986 + 10, 1 : 16) + fac210_var_1018 * absa_var_163(ind0_var_986 + 11, 1 : 16))
      ELSE IF (specparm_var_999 .GT. 0.875D0) THEN
        tau_major_var_1035(1 : ng3) = speccomb_var_979 * (fac200_var_1015 * absa_var_163(ind0_var_986 - 1, 1 : 16) + fac100_var_1014 * absa_var_163(ind0_var_986, 1 : 16) + fac000_var_1013 * absa_var_163(ind0_var_986 + 1, 1 : 16) + fac210_var_1018 * absa_var_163(ind0_var_986 + 8, 1 : 16) + fac110_var_1017 * absa_var_163(ind0_var_986 + 9, 1 : 16) + fac010_var_1016 * absa_var_163(ind0_var_986 + 10, 1 : 16))
      ELSE
        tau_major_var_1035(1 : ng3) = speccomb_var_979 * (fac000_var_1013 * absa_var_163(ind0_var_986, 1 : 16) + fac100_var_1014 * absa_var_163(ind0_var_986 + 1, 1 : 16) + fac010_var_1016 * absa_var_163(ind0_var_986 + 9, 1 : 16) + fac110_var_1017 * absa_var_163(ind0_var_986 + 10, 1 : 16))
      END IF
      IF (specparm1_var_1002 .LT. 0.125D0) THEN
        tau_major1_var_1036(1 : ng3) = speccomb1_var_980 * (fac001_var_1019 * absa_var_163(ind1_var_987, 1 : 16) + fac101_var_1020 * absa_var_163(ind1_var_987 + 1, 1 : 16) + fac201_var_1021 * absa_var_163(ind1_var_987 + 2, 1 : 16) + fac011_var_1022 * absa_var_163(ind1_var_987 + 9, 1 : 16) + fac111_var_1023 * absa_var_163(ind1_var_987 + 10, 1 : 16) + fac211_var_1024 * absa_var_163(ind1_var_987 + 11, 1 : 16))
      ELSE IF (specparm1_var_1002 .GT. 0.875D0) THEN
        tau_major1_var_1036(1 : ng3) = speccomb1_var_980 * (fac201_var_1021 * absa_var_163(ind1_var_987 - 1, 1 : 16) + fac101_var_1020 * absa_var_163(ind1_var_987, 1 : 16) + fac001_var_1019 * absa_var_163(ind1_var_987 + 1, 1 : 16) + fac211_var_1024 * absa_var_163(ind1_var_987 + 8, 1 : 16) + fac111_var_1023 * absa_var_163(ind1_var_987 + 9, 1 : 16) + fac011_var_1022 * absa_var_163(ind1_var_987 + 10, 1 : 16))
      ELSE
        tau_major1_var_1036(1 : ng3) = speccomb1_var_980 * (fac001_var_1019 * absa_var_163(ind1_var_987, 1 : 16) + fac101_var_1020 * absa_var_163(ind1_var_987 + 1, 1 : 16) + fac011_var_1022 * absa_var_163(ind1_var_987 + 9, 1 : 16) + fac111_var_1023 * absa_var_163(ind1_var_987 + 10, 1 : 16))
      END IF
      DO ig_var_991 = 1, 16
        tauself_var_1031 = selffac_var_969(jl_var_1047, lay_var_993) * (selfref_var_165(inds_var_988, ig_var_991) + selffrac_var_970(jl_var_1047, lay_var_993) * (selfref_var_165(inds_var_988 + 1, ig_var_991) - selfref_var_165(inds_var_988, ig_var_991)))
        taufor_var_1030 = forfac_var_959(jl_var_1047, lay_var_993) * (forref_var_166(indf_var_989, ig_var_991) + forfrac_var_976(jl_var_1047, lay_var_993) * (forref_var_166(indf_var_989 + 1, ig_var_991) - forref_var_166(indf_var_989, ig_var_991)))
        n2om1_var_1032 = ka_mn2o_var_161(jmn2o_var_995, indm_var_990, ig_var_991) + fmn2o_var_1003 * (ka_mn2o_var_161(jmn2o_var_995 + 1, indm_var_990, ig_var_991) - ka_mn2o_var_161(jmn2o_var_995, indm_var_990, ig_var_991))
        n2om2_var_1033 = ka_mn2o_var_161(jmn2o_var_995, indm_var_990 + 1, ig_var_991) + fmn2o_var_1003 * (ka_mn2o_var_161(jmn2o_var_995 + 1, indm_var_990 + 1, ig_var_991) - ka_mn2o_var_161(jmn2o_var_995, indm_var_990 + 1, ig_var_991))
        absn2o_var_1034 = n2om1_var_1032 + minorfrac_var_977(jl_var_1047, lay_var_993) * (n2om2_var_1033 - n2om1_var_1032)
        taug_var_953(jl_var_1047, 22 + ig_var_991, lay_var_993) = tau_major_var_1035(ig_var_991) + tau_major1_var_1036(ig_var_991) + tauself_var_1031 + taufor_var_1030 + adjcoln2o_var_1010 * absn2o_var_1034
        fracs_var_972(jl_var_1047, 22 + ig_var_991, lay_var_993) = fracrefa_var_159(ig_var_991, jpl_var_996) + fpl_var_1006 * (fracrefa_var_159(ig_var_991, jpl_var_996 + 1) - fracrefa_var_159(ig_var_991, jpl_var_996))
      END DO
    END DO
  END DO
  DO lay_var_993 = laytrop_max_var_1038 + 1, klev_var_952
    DO jl_var_1047 = kidia_var_950, kfdia_var_951
      speccomb_var_979 = colh2o_var_964(jl_var_1047, lay_var_993) + rat_h2oco2_var_973(jl_var_1047, lay_var_993) * colco2_var_965(jl_var_1047, lay_var_993)
      specparm_var_999 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_var_979, oneminus_var_963)
      specmult_var_998 = 4.0D0 * (specparm_var_999)
      js_var_992 = 1 + INT(specmult_var_998)
      fs_var_997 = ((specmult_var_998) - AINT((specmult_var_998)))
      speccomb1_var_980 = colh2o_var_964(jl_var_1047, lay_var_993) + rat_h2oco2_1_var_974(jl_var_1047, lay_var_993) * colco2_var_965(jl_var_1047, lay_var_993)
      specparm1_var_1002 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb1_var_980, oneminus_var_963)
      specmult1_var_1001 = 4.0D0 * (specparm1_var_1002)
      js1_var_994 = 1 + INT(specmult1_var_1001)
      fs1_var_1000 = ((specmult1_var_1001) - AINT((specmult1_var_1001)))
      fac000_var_1013 = (1.0D0 - fs_var_997) * fac00_var_955(jl_var_1047, lay_var_993)
      fac010_var_1016 = (1.0D0 - fs_var_997) * fac10_var_957(jl_var_1047, lay_var_993)
      fac100_var_1014 = fs_var_997 * fac00_var_955(jl_var_1047, lay_var_993)
      fac110_var_1017 = fs_var_997 * fac10_var_957(jl_var_1047, lay_var_993)
      fac001_var_1019 = (1.0D0 - fs1_var_1000) * fac01_var_956(jl_var_1047, lay_var_993)
      fac011_var_1022 = (1.0D0 - fs1_var_1000) * fac11_var_958(jl_var_1047, lay_var_993)
      fac101_var_1020 = fs1_var_1000 * fac01_var_956(jl_var_1047, lay_var_993)
      fac111_var_1023 = fs1_var_1000 * fac11_var_958(jl_var_1047, lay_var_993)
      speccomb_mn2o_var_981 = colh2o_var_964(jl_var_1047, lay_var_993) + refrat_m_b * colco2_var_965(jl_var_1047, lay_var_993)
      specparm_mn2o_var_1005 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_mn2o_var_981, oneminus_var_963)
      specmult_mn2o_var_1004 = 4.0D0 * specparm_mn2o_var_1005
      jmn2o_var_995 = 1 + INT(specmult_mn2o_var_1004)
      fmn2o_var_1003 = ((specmult_mn2o_var_1004) - AINT((specmult_mn2o_var_1004)))
      chi_n2o_var_1012 = coln2o_var_966(jl_var_1047, lay_var_993) / coldry_var_967(jl_var_1047, lay_var_993)
      ratn2o_var_1011 = 1D+20 * chi_n2o_var_1012 / chi_mls(4, jp_var_960(jl_var_1047, lay_var_993) + 1)
      IF (ratn2o_var_1011 .GT. 1.5D0) THEN
        adjfac_var_1009 = 0.5D0 + (ratn2o_var_1011 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1010 = adjfac_var_1009 * chi_mls(4, jp_var_960(jl_var_1047, lay_var_993) + 1) * coldry_var_967(jl_var_1047, lay_var_993) * 1D-20
      ELSE
        adjcoln2o_var_1010 = coln2o_var_966(jl_var_1047, lay_var_993)
      END IF
      speccomb_planck_var_982 = colh2o_var_964(jl_var_1047, lay_var_993) + refrat_planck_b_var_984 * colco2_var_965(jl_var_1047, lay_var_993)
      specparm_planck_var_1008 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_planck_var_982, oneminus_var_963)
      specmult_planck_var_1007 = 4.0D0 * specparm_planck_var_1008
      jpl_var_996 = 1 + INT(specmult_planck_var_1007)
      fpl_var_1006 = ((specmult_planck_var_1007) - AINT((specmult_planck_var_1007)))
      ind0_var_986 = ((jp_var_960(jl_var_1047, lay_var_993) - 13) * 5 + (jt_var_961(jl_var_1047, lay_var_993) - 1)) * nspb_var_217(3) + js_var_992
      ind1_var_987 = ((jp_var_960(jl_var_1047, lay_var_993) - 12) * 5 + (jt1_var_962(jl_var_1047, lay_var_993) - 1)) * nspb_var_217(3) + js1_var_994
      indf_var_989 = indfor_var_975(jl_var_1047, lay_var_993)
      indm_var_990 = indminor_var_978(jl_var_1047, lay_var_993)
      DO ig_var_991 = 1, 16
        taufor_var_1030 = forfac_var_959(jl_var_1047, lay_var_993) * (forref_var_166(indf_var_989, ig_var_991) + forfrac_var_976(jl_var_1047, lay_var_993) * (forref_var_166(indf_var_989 + 1, ig_var_991) - forref_var_166(indf_var_989, ig_var_991)))
        n2om1_var_1032 = kb_mn2o_var_162(jmn2o_var_995, indm_var_990, ig_var_991) + fmn2o_var_1003 * (kb_mn2o_var_162(jmn2o_var_995 + 1, indm_var_990, ig_var_991) - kb_mn2o_var_162(jmn2o_var_995, indm_var_990, ig_var_991))
        n2om2_var_1033 = kb_mn2o_var_162(jmn2o_var_995, indm_var_990 + 1, ig_var_991) + fmn2o_var_1003 * (kb_mn2o_var_162(jmn2o_var_995 + 1, indm_var_990 + 1, ig_var_991) - kb_mn2o_var_162(jmn2o_var_995, indm_var_990 + 1, ig_var_991))
        absn2o_var_1034 = n2om1_var_1032 + minorfrac_var_977(jl_var_1047, lay_var_993) * (n2om2_var_1033 - n2om1_var_1032)
        taug_var_953(jl_var_1047, 22 + ig_var_991, lay_var_993) = speccomb_var_979 * (fac000_var_1013 * absb_var_164(ind0_var_986, ig_var_991) + fac100_var_1014 * absb_var_164(ind0_var_986 + 1, ig_var_991) + fac010_var_1016 * absb_var_164(ind0_var_986 + 5, ig_var_991) + fac110_var_1017 * absb_var_164(ind0_var_986 + 6, ig_var_991)) + speccomb1_var_980 * (fac001_var_1019 * absb_var_164(ind1_var_987, ig_var_991) + fac101_var_1020 * absb_var_164(ind1_var_987 + 1, ig_var_991) + fac011_var_1022 * absb_var_164(ind1_var_987 + 5, ig_var_991) + fac111_var_1023 * absb_var_164(ind1_var_987 + 6, ig_var_991)) + taufor_var_1030 + adjcoln2o_var_1010 * absn2o_var_1034
        fracs_var_972(jl_var_1047, 22 + ig_var_991, lay_var_993) = fracrefb_var_160(ig_var_991, jpl_var_996) + fpl_var_1006 * (fracrefb_var_160(ig_var_991, jpl_var_996 + 1) - fracrefb_var_160(ig_var_991, jpl_var_996))
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1038 /= laytrop_min_var_1037) THEN
    DO lay_var_993 = laytrop_min_var_1037 + 1, laytrop_max_var_1038
      ixc0_var_1044 = ixc_var_1039(lay_var_993)
      DO ixp_var_1045 = 1, ixc0_var_1044
        jl_var_1047 = ixlow_var_1040(ixp_var_1045, lay_var_993)
        speccomb_var_979 = colh2o_var_964(jl_var_1047, lay_var_993) + rat_h2oco2_var_973(jl_var_1047, lay_var_993) * colco2_var_965(jl_var_1047, lay_var_993)
        specparm_var_999 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_var_979, oneminus_var_963)
        specmult_var_998 = 8.0D0 * (specparm_var_999)
        js_var_992 = 1 + INT(specmult_var_998)
        fs_var_997 = ((specmult_var_998) - AINT((specmult_var_998)))
        speccomb1_var_980 = colh2o_var_964(jl_var_1047, lay_var_993) + rat_h2oco2_1_var_974(jl_var_1047, lay_var_993) * colco2_var_965(jl_var_1047, lay_var_993)
        specparm1_var_1002 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb1_var_980, oneminus_var_963)
        specmult1_var_1001 = 8.0D0 * (specparm1_var_1002)
        js1_var_994 = 1 + INT(specmult1_var_1001)
        fs1_var_1000 = ((specmult1_var_1001) - AINT((specmult1_var_1001)))
        speccomb_mn2o_var_981 = colh2o_var_964(jl_var_1047, lay_var_993) + refrat_m_a_var_985 * colco2_var_965(jl_var_1047, lay_var_993)
        specparm_mn2o_var_1005 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_mn2o_var_981, oneminus_var_963)
        specmult_mn2o_var_1004 = 8.0D0 * specparm_mn2o_var_1005
        jmn2o_var_995 = 1 + INT(specmult_mn2o_var_1004)
        fmn2o_var_1003 = ((specmult_mn2o_var_1004) - AINT((specmult_mn2o_var_1004)))
        chi_n2o_var_1012 = coln2o_var_966(jl_var_1047, lay_var_993) / coldry_var_967(jl_var_1047, lay_var_993)
        ratn2o_var_1011 = 1D+20 * chi_n2o_var_1012 / chi_mls(4, jp_var_960(jl_var_1047, lay_var_993) + 1)
        IF (ratn2o_var_1011 .GT. 1.5D0) THEN
          adjfac_var_1009 = 0.5D0 + (ratn2o_var_1011 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1010 = adjfac_var_1009 * chi_mls(4, jp_var_960(jl_var_1047, lay_var_993) + 1) * coldry_var_967(jl_var_1047, lay_var_993) * 1D-20
        ELSE
          adjcoln2o_var_1010 = coln2o_var_966(jl_var_1047, lay_var_993)
        END IF
        speccomb_planck_var_982 = colh2o_var_964(jl_var_1047, lay_var_993) + refrat_planck_a_var_983 * colco2_var_965(jl_var_1047, lay_var_993)
        specparm_planck_var_1008 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_planck_var_982, oneminus_var_963)
        specmult_planck_var_1007 = 8.0D0 * specparm_planck_var_1008
        jpl_var_996 = 1 + INT(specmult_planck_var_1007)
        fpl_var_1006 = ((specmult_planck_var_1007) - AINT((specmult_planck_var_1007)))
        ind0_var_986 = ((jp_var_960(jl_var_1047, lay_var_993) - 1) * 5 + (jt_var_961(jl_var_1047, lay_var_993) - 1)) * nspa_var_216(3) + js_var_992
        ind1_var_987 = (jp_var_960(jl_var_1047, lay_var_993) * 5 + (jt1_var_962(jl_var_1047, lay_var_993) - 1)) * nspa_var_216(3) + js1_var_994
        inds_var_988 = indself_var_971(jl_var_1047, lay_var_993)
        indf_var_989 = indfor_var_975(jl_var_1047, lay_var_993)
        indm_var_990 = indminor_var_978(jl_var_1047, lay_var_993)
        IF (specparm_var_999 .LT. 0.125D0) THEN
          p_var_1025 = fs_var_997 - 1.0D0
          p4_var_1026 = p_var_1025 ** 4
          fk0_var_1027 = p4_var_1026
          fk1_var_1028 = 1.0D0 - p_var_1025 - 2.0D0 * p4_var_1026
          fk2_var_1029 = p_var_1025 + p4_var_1026
          fac000_var_1013 = fk0_var_1027 * fac00_var_955(jl_var_1047, lay_var_993)
          fac100_var_1014 = fk1_var_1028 * fac00_var_955(jl_var_1047, lay_var_993)
          fac200_var_1015 = fk2_var_1029 * fac00_var_955(jl_var_1047, lay_var_993)
          fac010_var_1016 = fk0_var_1027 * fac10_var_957(jl_var_1047, lay_var_993)
          fac110_var_1017 = fk1_var_1028 * fac10_var_957(jl_var_1047, lay_var_993)
          fac210_var_1018 = fk2_var_1029 * fac10_var_957(jl_var_1047, lay_var_993)
        ELSE IF (specparm_var_999 .GT. 0.875D0) THEN
          p_var_1025 = - fs_var_997
          p4_var_1026 = p_var_1025 ** 4
          fk0_var_1027 = p4_var_1026
          fk1_var_1028 = 1.0D0 - p_var_1025 - 2.0D0 * p4_var_1026
          fk2_var_1029 = p_var_1025 + p4_var_1026
          fac000_var_1013 = fk0_var_1027 * fac00_var_955(jl_var_1047, lay_var_993)
          fac100_var_1014 = fk1_var_1028 * fac00_var_955(jl_var_1047, lay_var_993)
          fac200_var_1015 = fk2_var_1029 * fac00_var_955(jl_var_1047, lay_var_993)
          fac010_var_1016 = fk0_var_1027 * fac10_var_957(jl_var_1047, lay_var_993)
          fac110_var_1017 = fk1_var_1028 * fac10_var_957(jl_var_1047, lay_var_993)
          fac210_var_1018 = fk2_var_1029 * fac10_var_957(jl_var_1047, lay_var_993)
        ELSE
          fac000_var_1013 = (1.0D0 - fs_var_997) * fac00_var_955(jl_var_1047, lay_var_993)
          fac010_var_1016 = (1.0D0 - fs_var_997) * fac10_var_957(jl_var_1047, lay_var_993)
          fac100_var_1014 = fs_var_997 * fac00_var_955(jl_var_1047, lay_var_993)
          fac110_var_1017 = fs_var_997 * fac10_var_957(jl_var_1047, lay_var_993)
          fac200_var_1015 = 0.0D0
          fac210_var_1018 = 0.0D0
        END IF
        IF (specparm1_var_1002 .LT. 0.125D0) THEN
          p_var_1025 = fs1_var_1000 - 1.0D0
          p4_var_1026 = p_var_1025 ** 4
          fk0_var_1027 = p4_var_1026
          fk1_var_1028 = 1.0D0 - p_var_1025 - 2.0D0 * p4_var_1026
          fk2_var_1029 = p_var_1025 + p4_var_1026
          fac001_var_1019 = fk0_var_1027 * fac01_var_956(jl_var_1047, lay_var_993)
          fac101_var_1020 = fk1_var_1028 * fac01_var_956(jl_var_1047, lay_var_993)
          fac201_var_1021 = fk2_var_1029 * fac01_var_956(jl_var_1047, lay_var_993)
          fac011_var_1022 = fk0_var_1027 * fac11_var_958(jl_var_1047, lay_var_993)
          fac111_var_1023 = fk1_var_1028 * fac11_var_958(jl_var_1047, lay_var_993)
          fac211_var_1024 = fk2_var_1029 * fac11_var_958(jl_var_1047, lay_var_993)
        ELSE IF (specparm1_var_1002 .GT. 0.875D0) THEN
          p_var_1025 = - fs1_var_1000
          p4_var_1026 = p_var_1025 ** 4
          fk0_var_1027 = p4_var_1026
          fk1_var_1028 = 1.0D0 - p_var_1025 - 2.0D0 * p4_var_1026
          fk2_var_1029 = p_var_1025 + p4_var_1026
          fac001_var_1019 = fk0_var_1027 * fac01_var_956(jl_var_1047, lay_var_993)
          fac101_var_1020 = fk1_var_1028 * fac01_var_956(jl_var_1047, lay_var_993)
          fac201_var_1021 = fk2_var_1029 * fac01_var_956(jl_var_1047, lay_var_993)
          fac011_var_1022 = fk0_var_1027 * fac11_var_958(jl_var_1047, lay_var_993)
          fac111_var_1023 = fk1_var_1028 * fac11_var_958(jl_var_1047, lay_var_993)
          fac211_var_1024 = fk2_var_1029 * fac11_var_958(jl_var_1047, lay_var_993)
        ELSE
          fac001_var_1019 = (1.0D0 - fs1_var_1000) * fac01_var_956(jl_var_1047, lay_var_993)
          fac011_var_1022 = (1.0D0 - fs1_var_1000) * fac11_var_958(jl_var_1047, lay_var_993)
          fac101_var_1020 = fs1_var_1000 * fac01_var_956(jl_var_1047, lay_var_993)
          fac111_var_1023 = fs1_var_1000 * fac11_var_958(jl_var_1047, lay_var_993)
          fac201_var_1021 = 0.0D0
          fac211_var_1024 = 0.0D0
        END IF
        IF (specparm_var_999 .LT. 0.125D0) THEN
          tau_major_var_1035(1 : ng3) = speccomb_var_979 * (fac000_var_1013 * absa_var_163(ind0_var_986, 1 : 16) + fac100_var_1014 * absa_var_163(ind0_var_986 + 1, 1 : 16) + fac200_var_1015 * absa_var_163(ind0_var_986 + 2, 1 : 16) + fac010_var_1016 * absa_var_163(ind0_var_986 + 9, 1 : 16) + fac110_var_1017 * absa_var_163(ind0_var_986 + 10, 1 : 16) + fac210_var_1018 * absa_var_163(ind0_var_986 + 11, 1 : 16))
        ELSE IF (specparm_var_999 .GT. 0.875D0) THEN
          tau_major_var_1035(1 : ng3) = speccomb_var_979 * (fac200_var_1015 * absa_var_163(ind0_var_986 - 1, 1 : 16) + fac100_var_1014 * absa_var_163(ind0_var_986, 1 : 16) + fac000_var_1013 * absa_var_163(ind0_var_986 + 1, 1 : 16) + fac210_var_1018 * absa_var_163(ind0_var_986 + 8, 1 : 16) + fac110_var_1017 * absa_var_163(ind0_var_986 + 9, 1 : 16) + fac010_var_1016 * absa_var_163(ind0_var_986 + 10, 1 : 16))
        ELSE
          tau_major_var_1035(1 : ng3) = speccomb_var_979 * (fac000_var_1013 * absa_var_163(ind0_var_986, 1 : 16) + fac100_var_1014 * absa_var_163(ind0_var_986 + 1, 1 : 16) + fac010_var_1016 * absa_var_163(ind0_var_986 + 9, 1 : 16) + fac110_var_1017 * absa_var_163(ind0_var_986 + 10, 1 : 16))
        END IF
        IF (specparm1_var_1002 .LT. 0.125D0) THEN
          tau_major1_var_1036(1 : ng3) = speccomb1_var_980 * (fac001_var_1019 * absa_var_163(ind1_var_987, 1 : 16) + fac101_var_1020 * absa_var_163(ind1_var_987 + 1, 1 : 16) + fac201_var_1021 * absa_var_163(ind1_var_987 + 2, 1 : 16) + fac011_var_1022 * absa_var_163(ind1_var_987 + 9, 1 : 16) + fac111_var_1023 * absa_var_163(ind1_var_987 + 10, 1 : 16) + fac211_var_1024 * absa_var_163(ind1_var_987 + 11, 1 : 16))
        ELSE IF (specparm1_var_1002 .GT. 0.875D0) THEN
          tau_major1_var_1036(1 : ng3) = speccomb1_var_980 * (fac201_var_1021 * absa_var_163(ind1_var_987 - 1, 1 : 16) + fac101_var_1020 * absa_var_163(ind1_var_987, 1 : 16) + fac001_var_1019 * absa_var_163(ind1_var_987 + 1, 1 : 16) + fac211_var_1024 * absa_var_163(ind1_var_987 + 8, 1 : 16) + fac111_var_1023 * absa_var_163(ind1_var_987 + 9, 1 : 16) + fac011_var_1022 * absa_var_163(ind1_var_987 + 10, 1 : 16))
        ELSE
          tau_major1_var_1036(1 : ng3) = speccomb1_var_980 * (fac001_var_1019 * absa_var_163(ind1_var_987, 1 : 16) + fac101_var_1020 * absa_var_163(ind1_var_987 + 1, 1 : 16) + fac011_var_1022 * absa_var_163(ind1_var_987 + 9, 1 : 16) + fac111_var_1023 * absa_var_163(ind1_var_987 + 10, 1 : 16))
        END IF
        DO ig_var_991 = 1, 16
          tauself_var_1031 = selffac_var_969(jl_var_1047, lay_var_993) * (selfref_var_165(inds_var_988, ig_var_991) + selffrac_var_970(jl_var_1047, lay_var_993) * (selfref_var_165(inds_var_988 + 1, ig_var_991) - selfref_var_165(inds_var_988, ig_var_991)))
          taufor_var_1030 = forfac_var_959(jl_var_1047, lay_var_993) * (forref_var_166(indf_var_989, ig_var_991) + forfrac_var_976(jl_var_1047, lay_var_993) * (forref_var_166(indf_var_989 + 1, ig_var_991) - forref_var_166(indf_var_989, ig_var_991)))
          n2om1_var_1032 = ka_mn2o_var_161(jmn2o_var_995, indm_var_990, ig_var_991) + fmn2o_var_1003 * (ka_mn2o_var_161(jmn2o_var_995 + 1, indm_var_990, ig_var_991) - ka_mn2o_var_161(jmn2o_var_995, indm_var_990, ig_var_991))
          n2om2_var_1033 = ka_mn2o_var_161(jmn2o_var_995, indm_var_990 + 1, ig_var_991) + fmn2o_var_1003 * (ka_mn2o_var_161(jmn2o_var_995 + 1, indm_var_990 + 1, ig_var_991) - ka_mn2o_var_161(jmn2o_var_995, indm_var_990 + 1, ig_var_991))
          absn2o_var_1034 = n2om1_var_1032 + minorfrac_var_977(jl_var_1047, lay_var_993) * (n2om2_var_1033 - n2om1_var_1032)
          taug_var_953(jl_var_1047, 22 + ig_var_991, lay_var_993) = tau_major_var_1035(ig_var_991) + tau_major1_var_1036(ig_var_991) + tauself_var_1031 + taufor_var_1030 + adjcoln2o_var_1010 * absn2o_var_1034
          fracs_var_972(jl_var_1047, 22 + ig_var_991, lay_var_993) = fracrefa_var_159(ig_var_991, jpl_var_996) + fpl_var_1006 * (fracrefa_var_159(ig_var_991, jpl_var_996 + 1) - fracrefa_var_159(ig_var_991, jpl_var_996))
        END DO
      END DO
      ixc0_var_1044 = kfdia_var_951 - kidia_var_950 + 1 - ixc0_var_1044
      DO ixp_var_1045 = 1, ixc0_var_1044
        jl_var_1047 = ixhigh_var_1041(ixp_var_1045, lay_var_993)
        speccomb_var_979 = colh2o_var_964(jl_var_1047, lay_var_993) + rat_h2oco2_var_973(jl_var_1047, lay_var_993) * colco2_var_965(jl_var_1047, lay_var_993)
        specparm_var_999 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_var_979, oneminus_var_963)
        specmult_var_998 = 4.0D0 * (specparm_var_999)
        js_var_992 = 1 + INT(specmult_var_998)
        fs_var_997 = ((specmult_var_998) - AINT((specmult_var_998)))
        speccomb1_var_980 = colh2o_var_964(jl_var_1047, lay_var_993) + rat_h2oco2_1_var_974(jl_var_1047, lay_var_993) * colco2_var_965(jl_var_1047, lay_var_993)
        specparm1_var_1002 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb1_var_980, oneminus_var_963)
        specmult1_var_1001 = 4.0D0 * (specparm1_var_1002)
        js1_var_994 = 1 + INT(specmult1_var_1001)
        fs1_var_1000 = ((specmult1_var_1001) - AINT((specmult1_var_1001)))
        fac000_var_1013 = (1.0D0 - fs_var_997) * fac00_var_955(jl_var_1047, lay_var_993)
        fac010_var_1016 = (1.0D0 - fs_var_997) * fac10_var_957(jl_var_1047, lay_var_993)
        fac100_var_1014 = fs_var_997 * fac00_var_955(jl_var_1047, lay_var_993)
        fac110_var_1017 = fs_var_997 * fac10_var_957(jl_var_1047, lay_var_993)
        fac001_var_1019 = (1.0D0 - fs1_var_1000) * fac01_var_956(jl_var_1047, lay_var_993)
        fac011_var_1022 = (1.0D0 - fs1_var_1000) * fac11_var_958(jl_var_1047, lay_var_993)
        fac101_var_1020 = fs1_var_1000 * fac01_var_956(jl_var_1047, lay_var_993)
        fac111_var_1023 = fs1_var_1000 * fac11_var_958(jl_var_1047, lay_var_993)
        speccomb_mn2o_var_981 = colh2o_var_964(jl_var_1047, lay_var_993) + refrat_m_b * colco2_var_965(jl_var_1047, lay_var_993)
        specparm_mn2o_var_1005 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_mn2o_var_981, oneminus_var_963)
        specmult_mn2o_var_1004 = 4.0D0 * specparm_mn2o_var_1005
        jmn2o_var_995 = 1 + INT(specmult_mn2o_var_1004)
        fmn2o_var_1003 = ((specmult_mn2o_var_1004) - AINT((specmult_mn2o_var_1004)))
        chi_n2o_var_1012 = coln2o_var_966(jl_var_1047, lay_var_993) / coldry_var_967(jl_var_1047, lay_var_993)
        ratn2o_var_1011 = 1D+20 * chi_n2o_var_1012 / chi_mls(4, jp_var_960(jl_var_1047, lay_var_993) + 1)
        IF (ratn2o_var_1011 .GT. 1.5D0) THEN
          adjfac_var_1009 = 0.5D0 + (ratn2o_var_1011 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1010 = adjfac_var_1009 * chi_mls(4, jp_var_960(jl_var_1047, lay_var_993) + 1) * coldry_var_967(jl_var_1047, lay_var_993) * 1D-20
        ELSE
          adjcoln2o_var_1010 = coln2o_var_966(jl_var_1047, lay_var_993)
        END IF
        speccomb_planck_var_982 = colh2o_var_964(jl_var_1047, lay_var_993) + refrat_planck_b_var_984 * colco2_var_965(jl_var_1047, lay_var_993)
        specparm_planck_var_1008 = MIN(colh2o_var_964(jl_var_1047, lay_var_993) / speccomb_planck_var_982, oneminus_var_963)
        specmult_planck_var_1007 = 4.0D0 * specparm_planck_var_1008
        jpl_var_996 = 1 + INT(specmult_planck_var_1007)
        fpl_var_1006 = ((specmult_planck_var_1007) - AINT((specmult_planck_var_1007)))
        ind0_var_986 = ((jp_var_960(jl_var_1047, lay_var_993) - 13) * 5 + (jt_var_961(jl_var_1047, lay_var_993) - 1)) * nspb_var_217(3) + js_var_992
        ind1_var_987 = ((jp_var_960(jl_var_1047, lay_var_993) - 12) * 5 + (jt1_var_962(jl_var_1047, lay_var_993) - 1)) * nspb_var_217(3) + js1_var_994
        indf_var_989 = indfor_var_975(jl_var_1047, lay_var_993)
        indm_var_990 = indminor_var_978(jl_var_1047, lay_var_993)
        DO ig_var_991 = 1, 16
          taufor_var_1030 = forfac_var_959(jl_var_1047, lay_var_993) * (forref_var_166(indf_var_989, ig_var_991) + forfrac_var_976(jl_var_1047, lay_var_993) * (forref_var_166(indf_var_989 + 1, ig_var_991) - forref_var_166(indf_var_989, ig_var_991)))
          n2om1_var_1032 = kb_mn2o_var_162(jmn2o_var_995, indm_var_990, ig_var_991) + fmn2o_var_1003 * (kb_mn2o_var_162(jmn2o_var_995 + 1, indm_var_990, ig_var_991) - kb_mn2o_var_162(jmn2o_var_995, indm_var_990, ig_var_991))
          n2om2_var_1033 = kb_mn2o_var_162(jmn2o_var_995, indm_var_990 + 1, ig_var_991) + fmn2o_var_1003 * (kb_mn2o_var_162(jmn2o_var_995 + 1, indm_var_990 + 1, ig_var_991) - kb_mn2o_var_162(jmn2o_var_995, indm_var_990 + 1, ig_var_991))
          absn2o_var_1034 = n2om1_var_1032 + minorfrac_var_977(jl_var_1047, lay_var_993) * (n2om2_var_1033 - n2om1_var_1032)
          taug_var_953(jl_var_1047, 22 + ig_var_991, lay_var_993) = speccomb_var_979 * (fac000_var_1013 * absb_var_164(ind0_var_986, ig_var_991) + fac100_var_1014 * absb_var_164(ind0_var_986 + 1, ig_var_991) + fac010_var_1016 * absb_var_164(ind0_var_986 + 5, ig_var_991) + fac110_var_1017 * absb_var_164(ind0_var_986 + 6, ig_var_991)) + speccomb1_var_980 * (fac001_var_1019 * absb_var_164(ind1_var_987, ig_var_991) + fac101_var_1020 * absb_var_164(ind1_var_987 + 1, ig_var_991) + fac011_var_1022 * absb_var_164(ind1_var_987 + 5, ig_var_991) + fac111_var_1023 * absb_var_164(ind1_var_987 + 6, ig_var_991)) + taufor_var_1030 + adjcoln2o_var_1010 * absn2o_var_1034
          fracs_var_972(jl_var_1047, 22 + ig_var_991, lay_var_993) = fracrefb_var_160(ig_var_991, jpl_var_996) + fpl_var_1006 * (fracrefb_var_160(ig_var_991, jpl_var_996 + 1) - fracrefb_var_160(ig_var_991, jpl_var_996))
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol3
SUBROUTINE rrtm_taumol9(kidia_var_1048, kfdia_var_1049, klev_var_1050, taug_var_1051, p_tauaerl_var_1052, fac00_var_1053, fac01_var_1054, fac10_var_1055, fac11_var_1056, forfac_var_1075, forfrac_var_1076, indfor_var_1074, jp_var_1057, jt_var_1058, jt1_var_1059, oneminus_var_1060, colh2o_var_1061, coln2o_var_1062, colch4_var_1063, coldry_var_1064, laytrop_var_1065, k_layswtch_var_1066, k_laylow_var_1067, selffac_var_1068, selffrac_var_1069, indself_var_1070, fracs_var_1071, rat_h2och4_var_1072, rat_h2och4_1_var_1073, minorfrac_var_1077, indminor_var_1078)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng9
  USE yoerrta9, ONLY: absa_var_208, absb_var_209, forref_var_213, fracrefa_var_206, fracrefb_var_207, ka_mn2o_var_210, kb_mn2o_var_211, selfref_var_212
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1048
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1049
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1050
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1051(kidia_var_1048 : kfdia_var_1049, 140, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1052(kidia_var_1048 : kfdia_var_1049, klev_var_1050, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1053(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1054(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1055(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1056(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1057(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1058(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1059(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1060
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1061(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_1062(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: colch4_var_1063(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1064(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1065(kidia_var_1048 : kfdia_var_1049)
  INTEGER(KIND = 4), INTENT(IN) :: k_layswtch_var_1066(kidia_var_1048 : kfdia_var_1049)
  INTEGER(KIND = 4), INTENT(IN) :: k_laylow_var_1067(kidia_var_1048 : kfdia_var_1049)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1068(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1069(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1070(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1071(kidia_var_1048 : kfdia_var_1049, 140, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_var_1072(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_1_var_1073(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1074(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1075(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1076(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1077(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1078(kidia_var_1048 : kfdia_var_1049, klev_var_1050)
  INTEGER(KIND = 4) :: ind0_var_1079, ind1_var_1080, inds_var_1081, indf_var_1082, indm_var_1083
  INTEGER(KIND = 4) :: ig_var_1084, js_var_1085, lay_var_1086, js1_var_1087, jmn2o_var_1088, jpl_var_1089
  REAL(KIND = 8) :: speccomb_var_1090, speccomb1_var_1091, speccomb_mn2o_var_1092, speccomb_planck_var_1093
  REAL(KIND = 8) :: refrat_planck_a_var_1094, refrat_m_a_var_1095
  REAL(KIND = 8) :: fs_var_1096, specmult_var_1097, specparm_var_1098, fs1_var_1099, specmult1_var_1100, specparm1_var_1101, fmn2o_var_1102, specmult_mn2o_var_1103, specparm_mn2o_var_1104, fpl_var_1105, specmult_planck_var_1106, specparm_planck_var_1107
  REAL(KIND = 8) :: adjfac_var_1108, adjcoln2o_var_1109, ratn2o_var_1110, chi_n2o_var_1111
  REAL(KIND = 8) :: fac000_var_1112, fac100_var_1113, fac200_var_1114, fac010_var_1115, fac110_var_1116, fac210_var_1117, fac001_var_1118, fac101_var_1119, fac201_var_1120, fac011_var_1121, fac111_var_1122, fac211_var_1123
  REAL(KIND = 8) :: p_var_1124, p4_var_1125, fk0_var_1126, fk1_var_1127, fk2_var_1128
  REAL(KIND = 8) :: taufor_var_1129, tauself_var_1130, n2om1_var_1131, n2om2_var_1132, absn2o_var_1133, tau_major_var_1134(12), tau_major1_var_1135(12)
  INTEGER(KIND = 4) :: laytrop_min_var_1136, laytrop_max_var_1137
  INTEGER(KIND = 4) :: ixc_var_1138(klev_var_1050), ixlow_var_1139(kfdia_var_1049, klev_var_1050), ixhigh_var_1140(kfdia_var_1049, klev_var_1050)
  INTEGER(KIND = 4) :: ich_var_1141, icl_var_1142, ixc0_var_1143, ixp_var_1144, jc_var_1145, jl_var_1146
  laytrop_min_var_1136 = MINVAL(laytrop_var_1065)
  laytrop_max_var_1137 = MAXVAL(laytrop_var_1065)
  ixlow_var_1139 = 0
  ixhigh_var_1140 = 0
  ixc_var_1138 = 0
  DO lay_var_1086 = laytrop_min_var_1136 + 1, laytrop_max_var_1137
    icl_var_1142 = 0
    ich_var_1141 = 0
    DO jc_var_1145 = kidia_var_1048, kfdia_var_1049
      IF (lay_var_1086 <= laytrop_var_1065(jc_var_1145)) THEN
        icl_var_1142 = icl_var_1142 + 1
        ixlow_var_1139(icl_var_1142, lay_var_1086) = jc_var_1145
      ELSE
        ich_var_1141 = ich_var_1141 + 1
        ixhigh_var_1140(ich_var_1141, lay_var_1086) = jc_var_1145
      END IF
    END DO
    ixc_var_1138(lay_var_1086) = icl_var_1142
  END DO
  refrat_planck_a_var_1094 = chi_mls(1, 9) / chi_mls(6, 9)
  refrat_m_a_var_1095 = chi_mls(1, 3) / chi_mls(6, 3)
  DO lay_var_1086 = 1, laytrop_min_var_1136
    DO jl_var_1146 = kidia_var_1048, kfdia_var_1049
      speccomb_var_1090 = colh2o_var_1061(jl_var_1146, lay_var_1086) + rat_h2och4_var_1072(jl_var_1146, lay_var_1086) * colch4_var_1063(jl_var_1146, lay_var_1086)
      specparm_var_1098 = MIN(colh2o_var_1061(jl_var_1146, lay_var_1086) / speccomb_var_1090, oneminus_var_1060)
      specmult_var_1097 = 8.0D0 * (specparm_var_1098)
      js_var_1085 = 1 + INT(specmult_var_1097)
      fs_var_1096 = ((specmult_var_1097) - AINT((specmult_var_1097)))
      speccomb1_var_1091 = colh2o_var_1061(jl_var_1146, lay_var_1086) + rat_h2och4_1_var_1073(jl_var_1146, lay_var_1086) * colch4_var_1063(jl_var_1146, lay_var_1086)
      specparm1_var_1101 = MIN(colh2o_var_1061(jl_var_1146, lay_var_1086) / speccomb1_var_1091, oneminus_var_1060)
      specmult1_var_1100 = 8.0D0 * (specparm1_var_1101)
      js1_var_1087 = 1 + INT(specmult1_var_1100)
      fs1_var_1099 = ((specmult1_var_1100) - AINT((specmult1_var_1100)))
      speccomb_mn2o_var_1092 = colh2o_var_1061(jl_var_1146, lay_var_1086) + refrat_m_a_var_1095 * colch4_var_1063(jl_var_1146, lay_var_1086)
      specparm_mn2o_var_1104 = MIN(colh2o_var_1061(jl_var_1146, lay_var_1086) / speccomb_mn2o_var_1092, oneminus_var_1060)
      specmult_mn2o_var_1103 = 8.0D0 * specparm_mn2o_var_1104
      jmn2o_var_1088 = 1 + INT(specmult_mn2o_var_1103)
      fmn2o_var_1102 = ((specmult_mn2o_var_1103) - AINT((specmult_mn2o_var_1103)))
      chi_n2o_var_1111 = coln2o_var_1062(jl_var_1146, lay_var_1086) / (coldry_var_1064(jl_var_1146, lay_var_1086))
      ratn2o_var_1110 = 1D+20 * chi_n2o_var_1111 / chi_mls(4, jp_var_1057(jl_var_1146, lay_var_1086) + 1)
      IF (ratn2o_var_1110 .GT. 1.5D0) THEN
        adjfac_var_1108 = 0.5D0 + (ratn2o_var_1110 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1109 = adjfac_var_1108 * chi_mls(4, jp_var_1057(jl_var_1146, lay_var_1086) + 1) * coldry_var_1064(jl_var_1146, lay_var_1086) * 1D-20
      ELSE
        adjcoln2o_var_1109 = coln2o_var_1062(jl_var_1146, lay_var_1086)
      END IF
      speccomb_planck_var_1093 = colh2o_var_1061(jl_var_1146, lay_var_1086) + refrat_planck_a_var_1094 * colch4_var_1063(jl_var_1146, lay_var_1086)
      specparm_planck_var_1107 = MIN(colh2o_var_1061(jl_var_1146, lay_var_1086) / speccomb_planck_var_1093, oneminus_var_1060)
      specmult_planck_var_1106 = 8.0D0 * specparm_planck_var_1107
      jpl_var_1089 = 1 + INT(specmult_planck_var_1106)
      fpl_var_1105 = ((specmult_planck_var_1106) - AINT((specmult_planck_var_1106)))
      ind0_var_1079 = ((jp_var_1057(jl_var_1146, lay_var_1086) - 1) * 5 + (jt_var_1058(jl_var_1146, lay_var_1086) - 1)) * nspa_var_216(9) + js_var_1085
      ind1_var_1080 = (jp_var_1057(jl_var_1146, lay_var_1086) * 5 + (jt1_var_1059(jl_var_1146, lay_var_1086) - 1)) * nspa_var_216(9) + js1_var_1087
      inds_var_1081 = indself_var_1070(jl_var_1146, lay_var_1086)
      indf_var_1082 = indfor_var_1074(jl_var_1146, lay_var_1086)
      indm_var_1083 = indminor_var_1078(jl_var_1146, lay_var_1086)
      IF (specparm_var_1098 .LT. 0.125D0) THEN
        p_var_1124 = fs_var_1096 - 1.0D0
        p4_var_1125 = p_var_1124 ** 4
        fk0_var_1126 = p4_var_1125
        fk1_var_1127 = 1.0D0 - p_var_1124 - 2.0D0 * p4_var_1125
        fk2_var_1128 = p_var_1124 + p4_var_1125
        fac000_var_1112 = fk0_var_1126 * fac00_var_1053(jl_var_1146, lay_var_1086)
        fac100_var_1113 = fk1_var_1127 * fac00_var_1053(jl_var_1146, lay_var_1086)
        fac200_var_1114 = fk2_var_1128 * fac00_var_1053(jl_var_1146, lay_var_1086)
        fac010_var_1115 = fk0_var_1126 * fac10_var_1055(jl_var_1146, lay_var_1086)
        fac110_var_1116 = fk1_var_1127 * fac10_var_1055(jl_var_1146, lay_var_1086)
        fac210_var_1117 = fk2_var_1128 * fac10_var_1055(jl_var_1146, lay_var_1086)
      ELSE IF (specparm_var_1098 .GT. 0.875D0) THEN
        p_var_1124 = - fs_var_1096
        p4_var_1125 = p_var_1124 ** 4
        fk0_var_1126 = p4_var_1125
        fk1_var_1127 = 1.0D0 - p_var_1124 - 2.0D0 * p4_var_1125
        fk2_var_1128 = p_var_1124 + p4_var_1125
        fac000_var_1112 = fk0_var_1126 * fac00_var_1053(jl_var_1146, lay_var_1086)
        fac100_var_1113 = fk1_var_1127 * fac00_var_1053(jl_var_1146, lay_var_1086)
        fac200_var_1114 = fk2_var_1128 * fac00_var_1053(jl_var_1146, lay_var_1086)
        fac010_var_1115 = fk0_var_1126 * fac10_var_1055(jl_var_1146, lay_var_1086)
        fac110_var_1116 = fk1_var_1127 * fac10_var_1055(jl_var_1146, lay_var_1086)
        fac210_var_1117 = fk2_var_1128 * fac10_var_1055(jl_var_1146, lay_var_1086)
      ELSE
        fac000_var_1112 = (1.0D0 - fs_var_1096) * fac00_var_1053(jl_var_1146, lay_var_1086)
        fac010_var_1115 = (1.0D0 - fs_var_1096) * fac10_var_1055(jl_var_1146, lay_var_1086)
        fac100_var_1113 = fs_var_1096 * fac00_var_1053(jl_var_1146, lay_var_1086)
        fac110_var_1116 = fs_var_1096 * fac10_var_1055(jl_var_1146, lay_var_1086)
        fac200_var_1114 = 0.0D0
        fac210_var_1117 = 0.0D0
      END IF
      IF (specparm1_var_1101 .LT. 0.125D0) THEN
        p_var_1124 = fs1_var_1099 - 1.0D0
        p4_var_1125 = p_var_1124 ** 4
        fk0_var_1126 = p4_var_1125
        fk1_var_1127 = 1.0D0 - p_var_1124 - 2.0D0 * p4_var_1125
        fk2_var_1128 = p_var_1124 + p4_var_1125
        fac001_var_1118 = fk0_var_1126 * fac01_var_1054(jl_var_1146, lay_var_1086)
        fac101_var_1119 = fk1_var_1127 * fac01_var_1054(jl_var_1146, lay_var_1086)
        fac201_var_1120 = fk2_var_1128 * fac01_var_1054(jl_var_1146, lay_var_1086)
        fac011_var_1121 = fk0_var_1126 * fac11_var_1056(jl_var_1146, lay_var_1086)
        fac111_var_1122 = fk1_var_1127 * fac11_var_1056(jl_var_1146, lay_var_1086)
        fac211_var_1123 = fk2_var_1128 * fac11_var_1056(jl_var_1146, lay_var_1086)
      ELSE IF (specparm1_var_1101 .GT. 0.875D0) THEN
        p_var_1124 = - fs1_var_1099
        p4_var_1125 = p_var_1124 ** 4
        fk0_var_1126 = p4_var_1125
        fk1_var_1127 = 1.0D0 - p_var_1124 - 2.0D0 * p4_var_1125
        fk2_var_1128 = p_var_1124 + p4_var_1125
        fac001_var_1118 = fk0_var_1126 * fac01_var_1054(jl_var_1146, lay_var_1086)
        fac101_var_1119 = fk1_var_1127 * fac01_var_1054(jl_var_1146, lay_var_1086)
        fac201_var_1120 = fk2_var_1128 * fac01_var_1054(jl_var_1146, lay_var_1086)
        fac011_var_1121 = fk0_var_1126 * fac11_var_1056(jl_var_1146, lay_var_1086)
        fac111_var_1122 = fk1_var_1127 * fac11_var_1056(jl_var_1146, lay_var_1086)
        fac211_var_1123 = fk2_var_1128 * fac11_var_1056(jl_var_1146, lay_var_1086)
      ELSE
        fac001_var_1118 = (1.0D0 - fs1_var_1099) * fac01_var_1054(jl_var_1146, lay_var_1086)
        fac011_var_1121 = (1.0D0 - fs1_var_1099) * fac11_var_1056(jl_var_1146, lay_var_1086)
        fac101_var_1119 = fs1_var_1099 * fac01_var_1054(jl_var_1146, lay_var_1086)
        fac111_var_1122 = fs1_var_1099 * fac11_var_1056(jl_var_1146, lay_var_1086)
        fac201_var_1120 = 0.0D0
        fac211_var_1123 = 0.0D0
      END IF
      IF (specparm_var_1098 .LT. 0.125D0) THEN
        tau_major_var_1134(1 : ng9) = speccomb_var_1090 * (fac000_var_1112 * absa_var_208(ind0_var_1079, 1 : 12) + fac100_var_1113 * absa_var_208(ind0_var_1079 + 1, 1 : 12) + fac200_var_1114 * absa_var_208(ind0_var_1079 + 2, 1 : 12) + fac010_var_1115 * absa_var_208(ind0_var_1079 + 9, 1 : 12) + fac110_var_1116 * absa_var_208(ind0_var_1079 + 10, 1 : 12) + fac210_var_1117 * absa_var_208(ind0_var_1079 + 11, 1 : 12))
      ELSE IF (specparm_var_1098 .GT. 0.875D0) THEN
        tau_major_var_1134(1 : ng9) = speccomb_var_1090 * (fac200_var_1114 * absa_var_208(ind0_var_1079 - 1, 1 : 12) + fac100_var_1113 * absa_var_208(ind0_var_1079, 1 : 12) + fac000_var_1112 * absa_var_208(ind0_var_1079 + 1, 1 : 12) + fac210_var_1117 * absa_var_208(ind0_var_1079 + 8, 1 : 12) + fac110_var_1116 * absa_var_208(ind0_var_1079 + 9, 1 : 12) + fac010_var_1115 * absa_var_208(ind0_var_1079 + 10, 1 : 12))
      ELSE
        tau_major_var_1134(1 : ng9) = speccomb_var_1090 * (fac000_var_1112 * absa_var_208(ind0_var_1079, 1 : 12) + fac100_var_1113 * absa_var_208(ind0_var_1079 + 1, 1 : 12) + fac010_var_1115 * absa_var_208(ind0_var_1079 + 9, 1 : 12) + fac110_var_1116 * absa_var_208(ind0_var_1079 + 10, 1 : 12))
      END IF
      IF (specparm1_var_1101 .LT. 0.125D0) THEN
        tau_major1_var_1135(1 : ng9) = speccomb1_var_1091 * (fac001_var_1118 * absa_var_208(ind1_var_1080, 1 : 12) + fac101_var_1119 * absa_var_208(ind1_var_1080 + 1, 1 : 12) + fac201_var_1120 * absa_var_208(ind1_var_1080 + 2, 1 : 12) + fac011_var_1121 * absa_var_208(ind1_var_1080 + 9, 1 : 12) + fac111_var_1122 * absa_var_208(ind1_var_1080 + 10, 1 : 12) + fac211_var_1123 * absa_var_208(ind1_var_1080 + 11, 1 : 12))
      ELSE IF (specparm1_var_1101 .GT. 0.875D0) THEN
        tau_major1_var_1135(1 : ng9) = speccomb1_var_1091 * (fac201_var_1120 * absa_var_208(ind1_var_1080 - 1, 1 : 12) + fac101_var_1119 * absa_var_208(ind1_var_1080, 1 : 12) + fac001_var_1118 * absa_var_208(ind1_var_1080 + 1, 1 : 12) + fac211_var_1123 * absa_var_208(ind1_var_1080 + 8, 1 : 12) + fac111_var_1122 * absa_var_208(ind1_var_1080 + 9, 1 : 12) + fac011_var_1121 * absa_var_208(ind1_var_1080 + 10, 1 : 12))
      ELSE
        tau_major1_var_1135(1 : ng9) = speccomb1_var_1091 * (fac001_var_1118 * absa_var_208(ind1_var_1080, 1 : 12) + fac101_var_1119 * absa_var_208(ind1_var_1080 + 1, 1 : 12) + fac011_var_1121 * absa_var_208(ind1_var_1080 + 9, 1 : 12) + fac111_var_1122 * absa_var_208(ind1_var_1080 + 10, 1 : 12))
      END IF
      DO ig_var_1084 = 1, 12
        tauself_var_1130 = selffac_var_1068(jl_var_1146, lay_var_1086) * (selfref_var_212(inds_var_1081, ig_var_1084) + selffrac_var_1069(jl_var_1146, lay_var_1086) * (selfref_var_212(inds_var_1081 + 1, ig_var_1084) - selfref_var_212(inds_var_1081, ig_var_1084)))
        taufor_var_1129 = forfac_var_1075(jl_var_1146, lay_var_1086) * (forref_var_213(indf_var_1082, ig_var_1084) + forfrac_var_1076(jl_var_1146, lay_var_1086) * (forref_var_213(indf_var_1082 + 1, ig_var_1084) - forref_var_213(indf_var_1082, ig_var_1084)))
        n2om1_var_1131 = ka_mn2o_var_210(jmn2o_var_1088, indm_var_1083, ig_var_1084) + fmn2o_var_1102 * (ka_mn2o_var_210(jmn2o_var_1088 + 1, indm_var_1083, ig_var_1084) - ka_mn2o_var_210(jmn2o_var_1088, indm_var_1083, ig_var_1084))
        n2om2_var_1132 = ka_mn2o_var_210(jmn2o_var_1088, indm_var_1083 + 1, ig_var_1084) + fmn2o_var_1102 * (ka_mn2o_var_210(jmn2o_var_1088 + 1, indm_var_1083 + 1, ig_var_1084) - ka_mn2o_var_210(jmn2o_var_1088, indm_var_1083 + 1, ig_var_1084))
        absn2o_var_1133 = n2om1_var_1131 + minorfrac_var_1077(jl_var_1146, lay_var_1086) * (n2om2_var_1132 - n2om1_var_1131)
        taug_var_1051(jl_var_1146, 96 + ig_var_1084, lay_var_1086) = tau_major_var_1134(ig_var_1084) + tau_major1_var_1135(ig_var_1084) + tauself_var_1130 + taufor_var_1129 + adjcoln2o_var_1109 * absn2o_var_1133
        fracs_var_1071(jl_var_1146, 96 + ig_var_1084, lay_var_1086) = fracrefa_var_206(ig_var_1084, jpl_var_1089) + fpl_var_1105 * (fracrefa_var_206(ig_var_1084, jpl_var_1089 + 1) - fracrefa_var_206(ig_var_1084, jpl_var_1089))
      END DO
    END DO
  END DO
  DO lay_var_1086 = laytrop_max_var_1137 + 1, klev_var_1050
    DO jl_var_1146 = kidia_var_1048, kfdia_var_1049
      chi_n2o_var_1111 = coln2o_var_1062(jl_var_1146, lay_var_1086) / (coldry_var_1064(jl_var_1146, lay_var_1086))
      ratn2o_var_1110 = 1D+20 * chi_n2o_var_1111 / chi_mls(4, jp_var_1057(jl_var_1146, lay_var_1086) + 1)
      IF (ratn2o_var_1110 .GT. 1.5D0) THEN
        adjfac_var_1108 = 0.5D0 + (ratn2o_var_1110 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1109 = adjfac_var_1108 * chi_mls(4, jp_var_1057(jl_var_1146, lay_var_1086) + 1) * coldry_var_1064(jl_var_1146, lay_var_1086) * 1D-20
      ELSE
        adjcoln2o_var_1109 = coln2o_var_1062(jl_var_1146, lay_var_1086)
      END IF
      ind0_var_1079 = ((jp_var_1057(jl_var_1146, lay_var_1086) - 13) * 5 + (jt_var_1058(jl_var_1146, lay_var_1086) - 1)) * nspb_var_217(9) + 1
      ind1_var_1080 = ((jp_var_1057(jl_var_1146, lay_var_1086) - 12) * 5 + (jt1_var_1059(jl_var_1146, lay_var_1086) - 1)) * nspb_var_217(9) + 1
      indm_var_1083 = indminor_var_1078(jl_var_1146, lay_var_1086)
      DO ig_var_1084 = 1, 12
        absn2o_var_1133 = kb_mn2o_var_211(indm_var_1083, ig_var_1084) + minorfrac_var_1077(jl_var_1146, lay_var_1086) * (kb_mn2o_var_211(indm_var_1083 + 1, ig_var_1084) - kb_mn2o_var_211(indm_var_1083, ig_var_1084))
        taug_var_1051(jl_var_1146, 96 + ig_var_1084, lay_var_1086) = colch4_var_1063(jl_var_1146, lay_var_1086) * (fac00_var_1053(jl_var_1146, lay_var_1086) * absb_var_209(ind0_var_1079, ig_var_1084) + fac10_var_1055(jl_var_1146, lay_var_1086) * absb_var_209(ind0_var_1079 + 1, ig_var_1084) + fac01_var_1054(jl_var_1146, lay_var_1086) * absb_var_209(ind1_var_1080, ig_var_1084) + fac11_var_1056(jl_var_1146, lay_var_1086) * absb_var_209(ind1_var_1080 + 1, ig_var_1084)) + adjcoln2o_var_1109 * absn2o_var_1133
        fracs_var_1071(jl_var_1146, 96 + ig_var_1084, lay_var_1086) = fracrefb_var_207(ig_var_1084)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1137 /= laytrop_min_var_1136) THEN
    DO lay_var_1086 = laytrop_min_var_1136 + 1, laytrop_max_var_1137
      ixc0_var_1143 = ixc_var_1138(lay_var_1086)
      DO ixp_var_1144 = 1, ixc0_var_1143
        jl_var_1146 = ixlow_var_1139(ixp_var_1144, lay_var_1086)
        speccomb_var_1090 = colh2o_var_1061(jl_var_1146, lay_var_1086) + rat_h2och4_var_1072(jl_var_1146, lay_var_1086) * colch4_var_1063(jl_var_1146, lay_var_1086)
        specparm_var_1098 = MIN(colh2o_var_1061(jl_var_1146, lay_var_1086) / speccomb_var_1090, oneminus_var_1060)
        specmult_var_1097 = 8.0D0 * (specparm_var_1098)
        js_var_1085 = 1 + INT(specmult_var_1097)
        fs_var_1096 = ((specmult_var_1097) - AINT((specmult_var_1097)))
        speccomb1_var_1091 = colh2o_var_1061(jl_var_1146, lay_var_1086) + rat_h2och4_1_var_1073(jl_var_1146, lay_var_1086) * colch4_var_1063(jl_var_1146, lay_var_1086)
        specparm1_var_1101 = MIN(colh2o_var_1061(jl_var_1146, lay_var_1086) / speccomb1_var_1091, oneminus_var_1060)
        specmult1_var_1100 = 8.0D0 * (specparm1_var_1101)
        js1_var_1087 = 1 + INT(specmult1_var_1100)
        fs1_var_1099 = ((specmult1_var_1100) - AINT((specmult1_var_1100)))
        speccomb_mn2o_var_1092 = colh2o_var_1061(jl_var_1146, lay_var_1086) + refrat_m_a_var_1095 * colch4_var_1063(jl_var_1146, lay_var_1086)
        specparm_mn2o_var_1104 = MIN(colh2o_var_1061(jl_var_1146, lay_var_1086) / speccomb_mn2o_var_1092, oneminus_var_1060)
        specmult_mn2o_var_1103 = 8.0D0 * specparm_mn2o_var_1104
        jmn2o_var_1088 = 1 + INT(specmult_mn2o_var_1103)
        fmn2o_var_1102 = ((specmult_mn2o_var_1103) - AINT((specmult_mn2o_var_1103)))
        chi_n2o_var_1111 = coln2o_var_1062(jl_var_1146, lay_var_1086) / (coldry_var_1064(jl_var_1146, lay_var_1086))
        ratn2o_var_1110 = 1D+20 * chi_n2o_var_1111 / chi_mls(4, jp_var_1057(jl_var_1146, lay_var_1086) + 1)
        IF (ratn2o_var_1110 .GT. 1.5D0) THEN
          adjfac_var_1108 = 0.5D0 + (ratn2o_var_1110 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1109 = adjfac_var_1108 * chi_mls(4, jp_var_1057(jl_var_1146, lay_var_1086) + 1) * coldry_var_1064(jl_var_1146, lay_var_1086) * 1D-20
        ELSE
          adjcoln2o_var_1109 = coln2o_var_1062(jl_var_1146, lay_var_1086)
        END IF
        speccomb_planck_var_1093 = colh2o_var_1061(jl_var_1146, lay_var_1086) + refrat_planck_a_var_1094 * colch4_var_1063(jl_var_1146, lay_var_1086)
        specparm_planck_var_1107 = MIN(colh2o_var_1061(jl_var_1146, lay_var_1086) / speccomb_planck_var_1093, oneminus_var_1060)
        specmult_planck_var_1106 = 8.0D0 * specparm_planck_var_1107
        jpl_var_1089 = 1 + INT(specmult_planck_var_1106)
        fpl_var_1105 = ((specmult_planck_var_1106) - AINT((specmult_planck_var_1106)))
        ind0_var_1079 = ((jp_var_1057(jl_var_1146, lay_var_1086) - 1) * 5 + (jt_var_1058(jl_var_1146, lay_var_1086) - 1)) * nspa_var_216(9) + js_var_1085
        ind1_var_1080 = (jp_var_1057(jl_var_1146, lay_var_1086) * 5 + (jt1_var_1059(jl_var_1146, lay_var_1086) - 1)) * nspa_var_216(9) + js1_var_1087
        inds_var_1081 = indself_var_1070(jl_var_1146, lay_var_1086)
        indf_var_1082 = indfor_var_1074(jl_var_1146, lay_var_1086)
        indm_var_1083 = indminor_var_1078(jl_var_1146, lay_var_1086)
        IF (specparm_var_1098 .LT. 0.125D0) THEN
          p_var_1124 = fs_var_1096 - 1.0D0
          p4_var_1125 = p_var_1124 ** 4
          fk0_var_1126 = p4_var_1125
          fk1_var_1127 = 1.0D0 - p_var_1124 - 2.0D0 * p4_var_1125
          fk2_var_1128 = p_var_1124 + p4_var_1125
          fac000_var_1112 = fk0_var_1126 * fac00_var_1053(jl_var_1146, lay_var_1086)
          fac100_var_1113 = fk1_var_1127 * fac00_var_1053(jl_var_1146, lay_var_1086)
          fac200_var_1114 = fk2_var_1128 * fac00_var_1053(jl_var_1146, lay_var_1086)
          fac010_var_1115 = fk0_var_1126 * fac10_var_1055(jl_var_1146, lay_var_1086)
          fac110_var_1116 = fk1_var_1127 * fac10_var_1055(jl_var_1146, lay_var_1086)
          fac210_var_1117 = fk2_var_1128 * fac10_var_1055(jl_var_1146, lay_var_1086)
        ELSE IF (specparm_var_1098 .GT. 0.875D0) THEN
          p_var_1124 = - fs_var_1096
          p4_var_1125 = p_var_1124 ** 4
          fk0_var_1126 = p4_var_1125
          fk1_var_1127 = 1.0D0 - p_var_1124 - 2.0D0 * p4_var_1125
          fk2_var_1128 = p_var_1124 + p4_var_1125
          fac000_var_1112 = fk0_var_1126 * fac00_var_1053(jl_var_1146, lay_var_1086)
          fac100_var_1113 = fk1_var_1127 * fac00_var_1053(jl_var_1146, lay_var_1086)
          fac200_var_1114 = fk2_var_1128 * fac00_var_1053(jl_var_1146, lay_var_1086)
          fac010_var_1115 = fk0_var_1126 * fac10_var_1055(jl_var_1146, lay_var_1086)
          fac110_var_1116 = fk1_var_1127 * fac10_var_1055(jl_var_1146, lay_var_1086)
          fac210_var_1117 = fk2_var_1128 * fac10_var_1055(jl_var_1146, lay_var_1086)
        ELSE
          fac000_var_1112 = (1.0D0 - fs_var_1096) * fac00_var_1053(jl_var_1146, lay_var_1086)
          fac010_var_1115 = (1.0D0 - fs_var_1096) * fac10_var_1055(jl_var_1146, lay_var_1086)
          fac100_var_1113 = fs_var_1096 * fac00_var_1053(jl_var_1146, lay_var_1086)
          fac110_var_1116 = fs_var_1096 * fac10_var_1055(jl_var_1146, lay_var_1086)
          fac200_var_1114 = 0.0D0
          fac210_var_1117 = 0.0D0
        END IF
        IF (specparm1_var_1101 .LT. 0.125D0) THEN
          p_var_1124 = fs1_var_1099 - 1.0D0
          p4_var_1125 = p_var_1124 ** 4
          fk0_var_1126 = p4_var_1125
          fk1_var_1127 = 1.0D0 - p_var_1124 - 2.0D0 * p4_var_1125
          fk2_var_1128 = p_var_1124 + p4_var_1125
          fac001_var_1118 = fk0_var_1126 * fac01_var_1054(jl_var_1146, lay_var_1086)
          fac101_var_1119 = fk1_var_1127 * fac01_var_1054(jl_var_1146, lay_var_1086)
          fac201_var_1120 = fk2_var_1128 * fac01_var_1054(jl_var_1146, lay_var_1086)
          fac011_var_1121 = fk0_var_1126 * fac11_var_1056(jl_var_1146, lay_var_1086)
          fac111_var_1122 = fk1_var_1127 * fac11_var_1056(jl_var_1146, lay_var_1086)
          fac211_var_1123 = fk2_var_1128 * fac11_var_1056(jl_var_1146, lay_var_1086)
        ELSE IF (specparm1_var_1101 .GT. 0.875D0) THEN
          p_var_1124 = - fs1_var_1099
          p4_var_1125 = p_var_1124 ** 4
          fk0_var_1126 = p4_var_1125
          fk1_var_1127 = 1.0D0 - p_var_1124 - 2.0D0 * p4_var_1125
          fk2_var_1128 = p_var_1124 + p4_var_1125
          fac001_var_1118 = fk0_var_1126 * fac01_var_1054(jl_var_1146, lay_var_1086)
          fac101_var_1119 = fk1_var_1127 * fac01_var_1054(jl_var_1146, lay_var_1086)
          fac201_var_1120 = fk2_var_1128 * fac01_var_1054(jl_var_1146, lay_var_1086)
          fac011_var_1121 = fk0_var_1126 * fac11_var_1056(jl_var_1146, lay_var_1086)
          fac111_var_1122 = fk1_var_1127 * fac11_var_1056(jl_var_1146, lay_var_1086)
          fac211_var_1123 = fk2_var_1128 * fac11_var_1056(jl_var_1146, lay_var_1086)
        ELSE
          fac001_var_1118 = (1.0D0 - fs1_var_1099) * fac01_var_1054(jl_var_1146, lay_var_1086)
          fac011_var_1121 = (1.0D0 - fs1_var_1099) * fac11_var_1056(jl_var_1146, lay_var_1086)
          fac101_var_1119 = fs1_var_1099 * fac01_var_1054(jl_var_1146, lay_var_1086)
          fac111_var_1122 = fs1_var_1099 * fac11_var_1056(jl_var_1146, lay_var_1086)
          fac201_var_1120 = 0.0D0
          fac211_var_1123 = 0.0D0
        END IF
        IF (specparm_var_1098 .LT. 0.125D0) THEN
          tau_major_var_1134(1 : ng9) = speccomb_var_1090 * (fac000_var_1112 * absa_var_208(ind0_var_1079, 1 : 12) + fac100_var_1113 * absa_var_208(ind0_var_1079 + 1, 1 : 12) + fac200_var_1114 * absa_var_208(ind0_var_1079 + 2, 1 : 12) + fac010_var_1115 * absa_var_208(ind0_var_1079 + 9, 1 : 12) + fac110_var_1116 * absa_var_208(ind0_var_1079 + 10, 1 : 12) + fac210_var_1117 * absa_var_208(ind0_var_1079 + 11, 1 : 12))
        ELSE IF (specparm_var_1098 .GT. 0.875D0) THEN
          tau_major_var_1134(1 : ng9) = speccomb_var_1090 * (fac200_var_1114 * absa_var_208(ind0_var_1079 - 1, 1 : 12) + fac100_var_1113 * absa_var_208(ind0_var_1079, 1 : 12) + fac000_var_1112 * absa_var_208(ind0_var_1079 + 1, 1 : 12) + fac210_var_1117 * absa_var_208(ind0_var_1079 + 8, 1 : 12) + fac110_var_1116 * absa_var_208(ind0_var_1079 + 9, 1 : 12) + fac010_var_1115 * absa_var_208(ind0_var_1079 + 10, 1 : 12))
        ELSE
          tau_major_var_1134(1 : ng9) = speccomb_var_1090 * (fac000_var_1112 * absa_var_208(ind0_var_1079, 1 : 12) + fac100_var_1113 * absa_var_208(ind0_var_1079 + 1, 1 : 12) + fac010_var_1115 * absa_var_208(ind0_var_1079 + 9, 1 : 12) + fac110_var_1116 * absa_var_208(ind0_var_1079 + 10, 1 : 12))
        END IF
        IF (specparm1_var_1101 .LT. 0.125D0) THEN
          tau_major1_var_1135(1 : ng9) = speccomb1_var_1091 * (fac001_var_1118 * absa_var_208(ind1_var_1080, 1 : 12) + fac101_var_1119 * absa_var_208(ind1_var_1080 + 1, 1 : 12) + fac201_var_1120 * absa_var_208(ind1_var_1080 + 2, 1 : 12) + fac011_var_1121 * absa_var_208(ind1_var_1080 + 9, 1 : 12) + fac111_var_1122 * absa_var_208(ind1_var_1080 + 10, 1 : 12) + fac211_var_1123 * absa_var_208(ind1_var_1080 + 11, 1 : 12))
        ELSE IF (specparm1_var_1101 .GT. 0.875D0) THEN
          tau_major1_var_1135(1 : ng9) = speccomb1_var_1091 * (fac201_var_1120 * absa_var_208(ind1_var_1080 - 1, 1 : 12) + fac101_var_1119 * absa_var_208(ind1_var_1080, 1 : 12) + fac001_var_1118 * absa_var_208(ind1_var_1080 + 1, 1 : 12) + fac211_var_1123 * absa_var_208(ind1_var_1080 + 8, 1 : 12) + fac111_var_1122 * absa_var_208(ind1_var_1080 + 9, 1 : 12) + fac011_var_1121 * absa_var_208(ind1_var_1080 + 10, 1 : 12))
        ELSE
          tau_major1_var_1135(1 : ng9) = speccomb1_var_1091 * (fac001_var_1118 * absa_var_208(ind1_var_1080, 1 : 12) + fac101_var_1119 * absa_var_208(ind1_var_1080 + 1, 1 : 12) + fac011_var_1121 * absa_var_208(ind1_var_1080 + 9, 1 : 12) + fac111_var_1122 * absa_var_208(ind1_var_1080 + 10, 1 : 12))
        END IF
        DO ig_var_1084 = 1, 12
          tauself_var_1130 = selffac_var_1068(jl_var_1146, lay_var_1086) * (selfref_var_212(inds_var_1081, ig_var_1084) + selffrac_var_1069(jl_var_1146, lay_var_1086) * (selfref_var_212(inds_var_1081 + 1, ig_var_1084) - selfref_var_212(inds_var_1081, ig_var_1084)))
          taufor_var_1129 = forfac_var_1075(jl_var_1146, lay_var_1086) * (forref_var_213(indf_var_1082, ig_var_1084) + forfrac_var_1076(jl_var_1146, lay_var_1086) * (forref_var_213(indf_var_1082 + 1, ig_var_1084) - forref_var_213(indf_var_1082, ig_var_1084)))
          n2om1_var_1131 = ka_mn2o_var_210(jmn2o_var_1088, indm_var_1083, ig_var_1084) + fmn2o_var_1102 * (ka_mn2o_var_210(jmn2o_var_1088 + 1, indm_var_1083, ig_var_1084) - ka_mn2o_var_210(jmn2o_var_1088, indm_var_1083, ig_var_1084))
          n2om2_var_1132 = ka_mn2o_var_210(jmn2o_var_1088, indm_var_1083 + 1, ig_var_1084) + fmn2o_var_1102 * (ka_mn2o_var_210(jmn2o_var_1088 + 1, indm_var_1083 + 1, ig_var_1084) - ka_mn2o_var_210(jmn2o_var_1088, indm_var_1083 + 1, ig_var_1084))
          absn2o_var_1133 = n2om1_var_1131 + minorfrac_var_1077(jl_var_1146, lay_var_1086) * (n2om2_var_1132 - n2om1_var_1131)
          taug_var_1051(jl_var_1146, 96 + ig_var_1084, lay_var_1086) = tau_major_var_1134(ig_var_1084) + tau_major1_var_1135(ig_var_1084) + tauself_var_1130 + taufor_var_1129 + adjcoln2o_var_1109 * absn2o_var_1133
          fracs_var_1071(jl_var_1146, 96 + ig_var_1084, lay_var_1086) = fracrefa_var_206(ig_var_1084, jpl_var_1089) + fpl_var_1105 * (fracrefa_var_206(ig_var_1084, jpl_var_1089 + 1) - fracrefa_var_206(ig_var_1084, jpl_var_1089))
        END DO
      END DO
      ixc0_var_1143 = kfdia_var_1049 - kidia_var_1048 + 1 - ixc0_var_1143
      DO ixp_var_1144 = 1, ixc0_var_1143
        jl_var_1146 = ixhigh_var_1140(ixp_var_1144, lay_var_1086)
        chi_n2o_var_1111 = coln2o_var_1062(jl_var_1146, lay_var_1086) / (coldry_var_1064(jl_var_1146, lay_var_1086))
        ratn2o_var_1110 = 1D+20 * chi_n2o_var_1111 / chi_mls(4, jp_var_1057(jl_var_1146, lay_var_1086) + 1)
        IF (ratn2o_var_1110 .GT. 1.5D0) THEN
          adjfac_var_1108 = 0.5D0 + (ratn2o_var_1110 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1109 = adjfac_var_1108 * chi_mls(4, jp_var_1057(jl_var_1146, lay_var_1086) + 1) * coldry_var_1064(jl_var_1146, lay_var_1086) * 1D-20
        ELSE
          adjcoln2o_var_1109 = coln2o_var_1062(jl_var_1146, lay_var_1086)
        END IF
        ind0_var_1079 = ((jp_var_1057(jl_var_1146, lay_var_1086) - 13) * 5 + (jt_var_1058(jl_var_1146, lay_var_1086) - 1)) * nspb_var_217(9) + 1
        ind1_var_1080 = ((jp_var_1057(jl_var_1146, lay_var_1086) - 12) * 5 + (jt1_var_1059(jl_var_1146, lay_var_1086) - 1)) * nspb_var_217(9) + 1
        indm_var_1083 = indminor_var_1078(jl_var_1146, lay_var_1086)
        DO ig_var_1084 = 1, 12
          absn2o_var_1133 = kb_mn2o_var_211(indm_var_1083, ig_var_1084) + minorfrac_var_1077(jl_var_1146, lay_var_1086) * (kb_mn2o_var_211(indm_var_1083 + 1, ig_var_1084) - kb_mn2o_var_211(indm_var_1083, ig_var_1084))
          taug_var_1051(jl_var_1146, 96 + ig_var_1084, lay_var_1086) = colch4_var_1063(jl_var_1146, lay_var_1086) * (fac00_var_1053(jl_var_1146, lay_var_1086) * absb_var_209(ind0_var_1079, ig_var_1084) + fac10_var_1055(jl_var_1146, lay_var_1086) * absb_var_209(ind0_var_1079 + 1, ig_var_1084) + fac01_var_1054(jl_var_1146, lay_var_1086) * absb_var_209(ind1_var_1080, ig_var_1084) + fac11_var_1056(jl_var_1146, lay_var_1086) * absb_var_209(ind1_var_1080 + 1, ig_var_1084)) + adjcoln2o_var_1109 * absn2o_var_1133
          fracs_var_1071(jl_var_1146, 96 + ig_var_1084, lay_var_1086) = fracrefb_var_207(ig_var_1084)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol9
SUBROUTINE srtm_taumol24(kidia_var_1147, kfdia_var_1148, klev_var_1149, p_fac00_var_1150, p_fac01_var_1151, p_fac10_var_1152, p_fac11_var_1153, k_jp_var_1154, k_jt_var_1155, k_jt1_var_1156, p_oneminus_var_1157, p_colh2o_var_1158, p_colmol_var_1159, p_colo2_var_1160, p_colo3_var_1161, k_laytrop_var_1162, p_selffac_var_1163, p_selffrac_var_1164, k_indself_var_1165, p_forfac_var_1166, p_forfrac_var_1167, k_indfor_var_1168, p_sfluxzen_var_1169, p_taug_var_1170, p_taur_var_1171, prmu0_var_1172)
  USE yoesrta24, ONLY: absa_var_280, absb_var_281, abso3ac_var_285, abso3bc_var_286, forrefc_var_283, layreffr_var_279, raylac, raylbc, selfrefc_var_282, sfluxrefc_var_284, strrat_var_278
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1147, kfdia_var_1148
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1149
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1150(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1151(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1152(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1153(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1154(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1155(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1156(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1157(kidia_var_1147 : kfdia_var_1148)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1158(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1159(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_1160(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1161(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1162(kidia_var_1147 : kfdia_var_1148)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1163(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1164(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1165(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1166(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1167(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1168(kidia_var_1147 : kfdia_var_1148, klev_var_1149)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1169(kidia_var_1147 : kfdia_var_1148, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1170(kidia_var_1147 : kfdia_var_1148, klev_var_1149, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1171(kidia_var_1147 : kfdia_var_1148, klev_var_1149, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1172(kidia_var_1147 : kfdia_var_1148)
  INTEGER(KIND = 4) :: ig_var_1173, ind0_var_1174, ind1_var_1175, inds_var_1176, indf_var_1177, js_var_1178, i_lay_var_1179, i_laysolfr_var_1180(kidia_var_1147 : kfdia_var_1148), i_nlayers_var_1181, iplon_var_1182
  INTEGER(KIND = 4) :: laytrop_min_var_1183, laytrop_max_var_1184
  REAL(KIND = 8) :: z_fs_var_1185, z_speccomb_var_1186, z_specmult_var_1187, z_specparm_var_1188, z_tauray_var_1189
  laytrop_min_var_1183 = MINVAL(k_laytrop_var_1162(kidia_var_1147 : kfdia_var_1148))
  laytrop_max_var_1184 = MAXVAL(k_laytrop_var_1162(kidia_var_1147 : kfdia_var_1148))
  i_nlayers_var_1181 = klev_var_1149
  DO iplon_var_1182 = kidia_var_1147, kfdia_var_1148
    i_laysolfr_var_1180(iplon_var_1182) = k_laytrop_var_1162(iplon_var_1182)
  END DO
  DO i_lay_var_1179 = 1, laytrop_min_var_1183
    DO iplon_var_1182 = kidia_var_1147, kfdia_var_1148
      IF (k_jp_var_1154(iplon_var_1182, i_lay_var_1179) < layreffr_var_279 .AND. k_jp_var_1154(iplon_var_1182, i_lay_var_1179 + 1) >= layreffr_var_279) i_laysolfr_var_1180(iplon_var_1182) = MIN(i_lay_var_1179 + 1, k_laytrop_var_1162(iplon_var_1182))
      z_speccomb_var_1186 = p_colh2o_var_1158(iplon_var_1182, i_lay_var_1179) + strrat_var_278 * p_colo2_var_1160(iplon_var_1182, i_lay_var_1179)
      z_specparm_var_1188 = p_colh2o_var_1158(iplon_var_1182, i_lay_var_1179) / z_speccomb_var_1186
      z_specparm_var_1188 = MIN(p_oneminus_var_1157(iplon_var_1182), z_specparm_var_1188)
      z_specmult_var_1187 = 8.0D0 * (z_specparm_var_1188)
      js_var_1178 = 1 + INT(z_specmult_var_1187)
      z_fs_var_1185 = z_specmult_var_1187 - AINT(z_specmult_var_1187)
      ind0_var_1174 = ((k_jp_var_1154(iplon_var_1182, i_lay_var_1179) - 1) * 5 + (k_jt_var_1155(iplon_var_1182, i_lay_var_1179) - 1)) * nspa_var_313(24) + js_var_1178
      ind1_var_1175 = (k_jp_var_1154(iplon_var_1182, i_lay_var_1179) * 5 + (k_jt1_var_1156(iplon_var_1182, i_lay_var_1179) - 1)) * nspa_var_313(24) + js_var_1178
      inds_var_1176 = k_indself_var_1165(iplon_var_1182, i_lay_var_1179)
      indf_var_1177 = k_indfor_var_1168(iplon_var_1182, i_lay_var_1179)
      DO ig_var_1173 = 1, 8
        z_tauray_var_1189 = p_colmol_var_1159(iplon_var_1182, i_lay_var_1179) * (raylac(ig_var_1173, js_var_1178) + z_fs_var_1185 * (raylac(ig_var_1173, js_var_1178 + 1) - raylac(ig_var_1173, js_var_1178)))
        p_taug_var_1170(iplon_var_1182, i_lay_var_1179, ig_var_1173) = z_speccomb_var_1186 * ((1.0D0 - z_fs_var_1185) * (absa_var_280(ind0_var_1174, ig_var_1173) * p_fac00_var_1150(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind0_var_1174 + 9, ig_var_1173) * p_fac10_var_1152(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind1_var_1175, ig_var_1173) * p_fac01_var_1151(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind1_var_1175 + 9, ig_var_1173) * p_fac11_var_1153(iplon_var_1182, i_lay_var_1179)) + z_fs_var_1185 * (absa_var_280(ind0_var_1174 + 1, ig_var_1173) * p_fac00_var_1150(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind0_var_1174 + 10, ig_var_1173) * p_fac10_var_1152(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind1_var_1175 + 1, ig_var_1173) * p_fac01_var_1151(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind1_var_1175 + 10, ig_var_1173) * p_fac11_var_1153(iplon_var_1182, i_lay_var_1179))) + p_colo3_var_1161(iplon_var_1182, i_lay_var_1179) * abso3ac_var_285(ig_var_1173) + p_colh2o_var_1158(iplon_var_1182, i_lay_var_1179) * (p_selffac_var_1163(iplon_var_1182, i_lay_var_1179) * (selfrefc_var_282(inds_var_1176, ig_var_1173) + p_selffrac_var_1164(iplon_var_1182, i_lay_var_1179) * (selfrefc_var_282(inds_var_1176 + 1, ig_var_1173) - selfrefc_var_282(inds_var_1176, ig_var_1173))) + p_forfac_var_1166(iplon_var_1182, i_lay_var_1179) * (forrefc_var_283(indf_var_1177, ig_var_1173) + p_forfrac_var_1167(iplon_var_1182, i_lay_var_1179) * (forrefc_var_283(indf_var_1177 + 1, ig_var_1173) - forrefc_var_283(indf_var_1177, ig_var_1173))))
        IF (i_lay_var_1179 == i_laysolfr_var_1180(iplon_var_1182)) p_sfluxzen_var_1169(iplon_var_1182, ig_var_1173) = sfluxrefc_var_284(ig_var_1173, js_var_1178) + z_fs_var_1185 * (sfluxrefc_var_284(ig_var_1173, js_var_1178 + 1) - sfluxrefc_var_284(ig_var_1173, js_var_1178))
        p_taur_var_1171(iplon_var_1182, i_lay_var_1179, ig_var_1173) = z_tauray_var_1189
      END DO
    END DO
  END DO
  DO i_lay_var_1179 = laytrop_min_var_1183 + 1, laytrop_max_var_1184
    DO iplon_var_1182 = kidia_var_1147, kfdia_var_1148
      IF (i_lay_var_1179 <= k_laytrop_var_1162(iplon_var_1182)) THEN
        IF (k_jp_var_1154(iplon_var_1182, i_lay_var_1179) < layreffr_var_279 .AND. k_jp_var_1154(iplon_var_1182, i_lay_var_1179 + 1) >= layreffr_var_279) i_laysolfr_var_1180(iplon_var_1182) = MIN(i_lay_var_1179 + 1, k_laytrop_var_1162(iplon_var_1182))
        z_speccomb_var_1186 = p_colh2o_var_1158(iplon_var_1182, i_lay_var_1179) + strrat_var_278 * p_colo2_var_1160(iplon_var_1182, i_lay_var_1179)
        z_specparm_var_1188 = p_colh2o_var_1158(iplon_var_1182, i_lay_var_1179) / z_speccomb_var_1186
        z_specparm_var_1188 = MIN(p_oneminus_var_1157(iplon_var_1182), z_specparm_var_1188)
        z_specmult_var_1187 = 8.0D0 * (z_specparm_var_1188)
        js_var_1178 = 1 + INT(z_specmult_var_1187)
        z_fs_var_1185 = z_specmult_var_1187 - AINT(z_specmult_var_1187)
        ind0_var_1174 = ((k_jp_var_1154(iplon_var_1182, i_lay_var_1179) - 1) * 5 + (k_jt_var_1155(iplon_var_1182, i_lay_var_1179) - 1)) * nspa_var_313(24) + js_var_1178
        ind1_var_1175 = (k_jp_var_1154(iplon_var_1182, i_lay_var_1179) * 5 + (k_jt1_var_1156(iplon_var_1182, i_lay_var_1179) - 1)) * nspa_var_313(24) + js_var_1178
        inds_var_1176 = k_indself_var_1165(iplon_var_1182, i_lay_var_1179)
        indf_var_1177 = k_indfor_var_1168(iplon_var_1182, i_lay_var_1179)
        DO ig_var_1173 = 1, 8
          z_tauray_var_1189 = p_colmol_var_1159(iplon_var_1182, i_lay_var_1179) * (raylac(ig_var_1173, js_var_1178) + z_fs_var_1185 * (raylac(ig_var_1173, js_var_1178 + 1) - raylac(ig_var_1173, js_var_1178)))
          p_taug_var_1170(iplon_var_1182, i_lay_var_1179, ig_var_1173) = z_speccomb_var_1186 * ((1.0D0 - z_fs_var_1185) * (absa_var_280(ind0_var_1174, ig_var_1173) * p_fac00_var_1150(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind0_var_1174 + 9, ig_var_1173) * p_fac10_var_1152(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind1_var_1175, ig_var_1173) * p_fac01_var_1151(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind1_var_1175 + 9, ig_var_1173) * p_fac11_var_1153(iplon_var_1182, i_lay_var_1179)) + z_fs_var_1185 * (absa_var_280(ind0_var_1174 + 1, ig_var_1173) * p_fac00_var_1150(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind0_var_1174 + 10, ig_var_1173) * p_fac10_var_1152(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind1_var_1175 + 1, ig_var_1173) * p_fac01_var_1151(iplon_var_1182, i_lay_var_1179) + absa_var_280(ind1_var_1175 + 10, ig_var_1173) * p_fac11_var_1153(iplon_var_1182, i_lay_var_1179))) + p_colo3_var_1161(iplon_var_1182, i_lay_var_1179) * abso3ac_var_285(ig_var_1173) + p_colh2o_var_1158(iplon_var_1182, i_lay_var_1179) * (p_selffac_var_1163(iplon_var_1182, i_lay_var_1179) * (selfrefc_var_282(inds_var_1176, ig_var_1173) + p_selffrac_var_1164(iplon_var_1182, i_lay_var_1179) * (selfrefc_var_282(inds_var_1176 + 1, ig_var_1173) - selfrefc_var_282(inds_var_1176, ig_var_1173))) + p_forfac_var_1166(iplon_var_1182, i_lay_var_1179) * (forrefc_var_283(indf_var_1177, ig_var_1173) + p_forfrac_var_1167(iplon_var_1182, i_lay_var_1179) * (forrefc_var_283(indf_var_1177 + 1, ig_var_1173) - forrefc_var_283(indf_var_1177, ig_var_1173))))
          IF (i_lay_var_1179 == i_laysolfr_var_1180(iplon_var_1182)) p_sfluxzen_var_1169(iplon_var_1182, ig_var_1173) = sfluxrefc_var_284(ig_var_1173, js_var_1178) + z_fs_var_1185 * (sfluxrefc_var_284(ig_var_1173, js_var_1178 + 1) - sfluxrefc_var_284(ig_var_1173, js_var_1178))
          p_taur_var_1171(iplon_var_1182, i_lay_var_1179, ig_var_1173) = z_tauray_var_1189
        END DO
      ELSE
        ind0_var_1174 = ((k_jp_var_1154(iplon_var_1182, i_lay_var_1179) - 13) * 5 + (k_jt_var_1155(iplon_var_1182, i_lay_var_1179) - 1)) * nspb_var_314(24) + 1
        ind1_var_1175 = ((k_jp_var_1154(iplon_var_1182, i_lay_var_1179) - 12) * 5 + (k_jt1_var_1156(iplon_var_1182, i_lay_var_1179) - 1)) * nspb_var_314(24) + 1
        DO ig_var_1173 = 1, 8
          z_tauray_var_1189 = p_colmol_var_1159(iplon_var_1182, i_lay_var_1179) * raylbc(ig_var_1173)
          p_taug_var_1170(iplon_var_1182, i_lay_var_1179, ig_var_1173) = p_colo2_var_1160(iplon_var_1182, i_lay_var_1179) * (p_fac00_var_1150(iplon_var_1182, i_lay_var_1179) * absb_var_281(ind0_var_1174, ig_var_1173) + p_fac10_var_1152(iplon_var_1182, i_lay_var_1179) * absb_var_281(ind0_var_1174 + 1, ig_var_1173) + p_fac01_var_1151(iplon_var_1182, i_lay_var_1179) * absb_var_281(ind1_var_1175, ig_var_1173) + p_fac11_var_1153(iplon_var_1182, i_lay_var_1179) * absb_var_281(ind1_var_1175 + 1, ig_var_1173)) + p_colo3_var_1161(iplon_var_1182, i_lay_var_1179) * abso3bc_var_286(ig_var_1173)
          p_taur_var_1171(iplon_var_1182, i_lay_var_1179, ig_var_1173) = z_tauray_var_1189
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1179 = laytrop_max_var_1184 + 1, i_nlayers_var_1181
    DO iplon_var_1182 = kidia_var_1147, kfdia_var_1148
      ind0_var_1174 = ((k_jp_var_1154(iplon_var_1182, i_lay_var_1179) - 13) * 5 + (k_jt_var_1155(iplon_var_1182, i_lay_var_1179) - 1)) * nspb_var_314(24) + 1
      ind1_var_1175 = ((k_jp_var_1154(iplon_var_1182, i_lay_var_1179) - 12) * 5 + (k_jt1_var_1156(iplon_var_1182, i_lay_var_1179) - 1)) * nspb_var_314(24) + 1
      DO ig_var_1173 = 1, 8
        z_tauray_var_1189 = p_colmol_var_1159(iplon_var_1182, i_lay_var_1179) * raylbc(ig_var_1173)
        p_taug_var_1170(iplon_var_1182, i_lay_var_1179, ig_var_1173) = p_colo2_var_1160(iplon_var_1182, i_lay_var_1179) * (p_fac00_var_1150(iplon_var_1182, i_lay_var_1179) * absb_var_281(ind0_var_1174, ig_var_1173) + p_fac10_var_1152(iplon_var_1182, i_lay_var_1179) * absb_var_281(ind0_var_1174 + 1, ig_var_1173) + p_fac01_var_1151(iplon_var_1182, i_lay_var_1179) * absb_var_281(ind1_var_1175, ig_var_1173) + p_fac11_var_1153(iplon_var_1182, i_lay_var_1179) * absb_var_281(ind1_var_1175 + 1, ig_var_1173)) + p_colo3_var_1161(iplon_var_1182, i_lay_var_1179) * abso3bc_var_286(ig_var_1173)
        p_taur_var_1171(iplon_var_1182, i_lay_var_1179, ig_var_1173) = z_tauray_var_1189
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol24
SUBROUTINE srtm_taumol29(kidia_var_1190, kfdia_var_1191, klev_var_1192, p_fac00_var_1193, p_fac01_var_1194, p_fac10_var_1195, p_fac11_var_1196, k_jp_var_1197, k_jt_var_1198, k_jt1_var_1199, p_colh2o_var_1200, p_colco2_var_1201, p_colmol_var_1202, k_laytrop_var_1203, p_selffac_var_1204, p_selffrac_var_1205, k_indself_var_1206, p_forfac_var_1207, p_forfrac_var_1208, k_indfor_var_1209, p_sfluxzen_var_1212, p_taug_var_1213, p_taur_var_1214, prmu0_var_1215)
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  USE yoesrta29, ONLY: absa_var_308, absb_var_309, absco2c, absh2oc, forrefc_var_311, layreffr_var_307, rayl_var_306, selfrefc_var_310, sfluxrefc_var_312
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1190, kfdia_var_1191
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1192
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1193(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1194(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1195(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1196(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1197(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1198(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1199(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1200(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1201(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1202(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1203(kidia_var_1190 : kfdia_var_1191)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1204(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1205(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1206(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1207(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1208(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1209(kidia_var_1190 : kfdia_var_1191, klev_var_1192)
  INTEGER(KIND = 4) :: laytrop_min_var_1210, laytrop_max_var_1211
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1212(kidia_var_1190 : kfdia_var_1191, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1213(kidia_var_1190 : kfdia_var_1191, klev_var_1192, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1214(kidia_var_1190 : kfdia_var_1191, klev_var_1192, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1215(kidia_var_1190 : kfdia_var_1191)
  INTEGER(KIND = 4) :: ig_var_1216, ind0_var_1217, ind1_var_1218, inds_var_1219, indf_var_1220, i_lay_var_1221, i_laysolfr_var_1222(kidia_var_1190 : kfdia_var_1191), i_nlayers_var_1223, iplon_var_1224
  REAL(KIND = 8) :: z_tauray_var_1225
  laytrop_min_var_1210 = MINVAL(k_laytrop_var_1203(kidia_var_1190 : kfdia_var_1191))
  laytrop_max_var_1211 = MAXVAL(k_laytrop_var_1203(kidia_var_1190 : kfdia_var_1191))
  i_nlayers_var_1223 = klev_var_1192
  DO iplon_var_1224 = kidia_var_1190, kfdia_var_1191
    i_laysolfr_var_1222(iplon_var_1224) = i_nlayers_var_1223
  END DO
  DO i_lay_var_1221 = 1, laytrop_min_var_1210
    DO iplon_var_1224 = kidia_var_1190, kfdia_var_1191
      ind0_var_1217 = ((k_jp_var_1197(iplon_var_1224, i_lay_var_1221) - 1) * 5 + (k_jt_var_1198(iplon_var_1224, i_lay_var_1221) - 1)) * nspa_var_313(29) + 1
      ind1_var_1218 = (k_jp_var_1197(iplon_var_1224, i_lay_var_1221) * 5 + (k_jt1_var_1199(iplon_var_1224, i_lay_var_1221) - 1)) * nspa_var_313(29) + 1
      inds_var_1219 = k_indself_var_1206(iplon_var_1224, i_lay_var_1221)
      indf_var_1220 = k_indfor_var_1209(iplon_var_1224, i_lay_var_1221)
      z_tauray_var_1225 = p_colmol_var_1202(iplon_var_1224, i_lay_var_1221) * rayl_var_306
      DO ig_var_1216 = 1, 12
        p_taug_var_1213(iplon_var_1224, i_lay_var_1221, ig_var_1216) = p_colh2o_var_1200(iplon_var_1224, i_lay_var_1221) * ((p_fac00_var_1193(iplon_var_1224, i_lay_var_1221) * absa_var_308(ind0_var_1217, ig_var_1216) + p_fac10_var_1195(iplon_var_1224, i_lay_var_1221) * absa_var_308(ind0_var_1217 + 1, ig_var_1216) + p_fac01_var_1194(iplon_var_1224, i_lay_var_1221) * absa_var_308(ind1_var_1218, ig_var_1216) + p_fac11_var_1196(iplon_var_1224, i_lay_var_1221) * absa_var_308(ind1_var_1218 + 1, ig_var_1216)) + p_selffac_var_1204(iplon_var_1224, i_lay_var_1221) * (selfrefc_var_310(inds_var_1219, ig_var_1216) + p_selffrac_var_1205(iplon_var_1224, i_lay_var_1221) * (selfrefc_var_310(inds_var_1219 + 1, ig_var_1216) - selfrefc_var_310(inds_var_1219, ig_var_1216))) + p_forfac_var_1207(iplon_var_1224, i_lay_var_1221) * (forrefc_var_311(indf_var_1220, ig_var_1216) + p_forfrac_var_1208(iplon_var_1224, i_lay_var_1221) * (forrefc_var_311(indf_var_1220 + 1, ig_var_1216) - forrefc_var_311(indf_var_1220, ig_var_1216)))) + p_colco2_var_1201(iplon_var_1224, i_lay_var_1221) * absco2c(ig_var_1216)
        p_taur_var_1214(iplon_var_1224, i_lay_var_1221, ig_var_1216) = z_tauray_var_1225
      END DO
    END DO
  END DO
  DO i_lay_var_1221 = laytrop_min_var_1210 + 1, laytrop_max_var_1211
    DO iplon_var_1224 = kidia_var_1190, kfdia_var_1191
      IF (i_lay_var_1221 <= k_laytrop_var_1203(iplon_var_1224)) THEN
        ind0_var_1217 = ((k_jp_var_1197(iplon_var_1224, i_lay_var_1221) - 1) * 5 + (k_jt_var_1198(iplon_var_1224, i_lay_var_1221) - 1)) * nspa_var_313(29) + 1
        ind1_var_1218 = (k_jp_var_1197(iplon_var_1224, i_lay_var_1221) * 5 + (k_jt1_var_1199(iplon_var_1224, i_lay_var_1221) - 1)) * nspa_var_313(29) + 1
        inds_var_1219 = k_indself_var_1206(iplon_var_1224, i_lay_var_1221)
        indf_var_1220 = k_indfor_var_1209(iplon_var_1224, i_lay_var_1221)
        z_tauray_var_1225 = p_colmol_var_1202(iplon_var_1224, i_lay_var_1221) * rayl_var_306
        DO ig_var_1216 = 1, 12
          p_taug_var_1213(iplon_var_1224, i_lay_var_1221, ig_var_1216) = p_colh2o_var_1200(iplon_var_1224, i_lay_var_1221) * ((p_fac00_var_1193(iplon_var_1224, i_lay_var_1221) * absa_var_308(ind0_var_1217, ig_var_1216) + p_fac10_var_1195(iplon_var_1224, i_lay_var_1221) * absa_var_308(ind0_var_1217 + 1, ig_var_1216) + p_fac01_var_1194(iplon_var_1224, i_lay_var_1221) * absa_var_308(ind1_var_1218, ig_var_1216) + p_fac11_var_1196(iplon_var_1224, i_lay_var_1221) * absa_var_308(ind1_var_1218 + 1, ig_var_1216)) + p_selffac_var_1204(iplon_var_1224, i_lay_var_1221) * (selfrefc_var_310(inds_var_1219, ig_var_1216) + p_selffrac_var_1205(iplon_var_1224, i_lay_var_1221) * (selfrefc_var_310(inds_var_1219 + 1, ig_var_1216) - selfrefc_var_310(inds_var_1219, ig_var_1216))) + p_forfac_var_1207(iplon_var_1224, i_lay_var_1221) * (forrefc_var_311(indf_var_1220, ig_var_1216) + p_forfrac_var_1208(iplon_var_1224, i_lay_var_1221) * (forrefc_var_311(indf_var_1220 + 1, ig_var_1216) - forrefc_var_311(indf_var_1220, ig_var_1216)))) + p_colco2_var_1201(iplon_var_1224, i_lay_var_1221) * absco2c(ig_var_1216)
          p_taur_var_1214(iplon_var_1224, i_lay_var_1221, ig_var_1216) = z_tauray_var_1225
        END DO
      ELSE
        IF (k_jp_var_1197(iplon_var_1224, i_lay_var_1221 - 1) < layreffr_var_307 .AND. k_jp_var_1197(iplon_var_1224, i_lay_var_1221) >= layreffr_var_307) i_laysolfr_var_1222(iplon_var_1224) = i_lay_var_1221
        ind0_var_1217 = ((k_jp_var_1197(iplon_var_1224, i_lay_var_1221) - 13) * 5 + (k_jt_var_1198(iplon_var_1224, i_lay_var_1221) - 1)) * nspb_var_314(29) + 1
        ind1_var_1218 = ((k_jp_var_1197(iplon_var_1224, i_lay_var_1221) - 12) * 5 + (k_jt1_var_1199(iplon_var_1224, i_lay_var_1221) - 1)) * nspb_var_314(29) + 1
        z_tauray_var_1225 = p_colmol_var_1202(iplon_var_1224, i_lay_var_1221) * rayl_var_306
        DO ig_var_1216 = 1, 12
          p_taug_var_1213(iplon_var_1224, i_lay_var_1221, ig_var_1216) = p_colco2_var_1201(iplon_var_1224, i_lay_var_1221) * (p_fac00_var_1193(iplon_var_1224, i_lay_var_1221) * absb_var_309(ind0_var_1217, ig_var_1216) + p_fac10_var_1195(iplon_var_1224, i_lay_var_1221) * absb_var_309(ind0_var_1217 + 1, ig_var_1216) + p_fac01_var_1194(iplon_var_1224, i_lay_var_1221) * absb_var_309(ind1_var_1218, ig_var_1216) + p_fac11_var_1196(iplon_var_1224, i_lay_var_1221) * absb_var_309(ind1_var_1218 + 1, ig_var_1216)) + p_colh2o_var_1200(iplon_var_1224, i_lay_var_1221) * absh2oc(ig_var_1216)
          IF (i_lay_var_1221 == i_laysolfr_var_1222(iplon_var_1224)) p_sfluxzen_var_1212(iplon_var_1224, ig_var_1216) = sfluxrefc_var_312(ig_var_1216)
          p_taur_var_1214(iplon_var_1224, i_lay_var_1221, ig_var_1216) = z_tauray_var_1225
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1221 = laytrop_max_var_1211 + 1, i_nlayers_var_1223
    DO iplon_var_1224 = kidia_var_1190, kfdia_var_1191
      IF (k_jp_var_1197(iplon_var_1224, i_lay_var_1221 - 1) < layreffr_var_307 .AND. k_jp_var_1197(iplon_var_1224, i_lay_var_1221) >= layreffr_var_307) i_laysolfr_var_1222(iplon_var_1224) = i_lay_var_1221
      ind0_var_1217 = ((k_jp_var_1197(iplon_var_1224, i_lay_var_1221) - 13) * 5 + (k_jt_var_1198(iplon_var_1224, i_lay_var_1221) - 1)) * nspb_var_314(29) + 1
      ind1_var_1218 = ((k_jp_var_1197(iplon_var_1224, i_lay_var_1221) - 12) * 5 + (k_jt1_var_1199(iplon_var_1224, i_lay_var_1221) - 1)) * nspb_var_314(29) + 1
      z_tauray_var_1225 = p_colmol_var_1202(iplon_var_1224, i_lay_var_1221) * rayl_var_306
      DO ig_var_1216 = 1, 12
        p_taug_var_1213(iplon_var_1224, i_lay_var_1221, ig_var_1216) = p_colco2_var_1201(iplon_var_1224, i_lay_var_1221) * (p_fac00_var_1193(iplon_var_1224, i_lay_var_1221) * absb_var_309(ind0_var_1217, ig_var_1216) + p_fac10_var_1195(iplon_var_1224, i_lay_var_1221) * absb_var_309(ind0_var_1217 + 1, ig_var_1216) + p_fac01_var_1194(iplon_var_1224, i_lay_var_1221) * absb_var_309(ind1_var_1218, ig_var_1216) + p_fac11_var_1196(iplon_var_1224, i_lay_var_1221) * absb_var_309(ind1_var_1218 + 1, ig_var_1216)) + p_colh2o_var_1200(iplon_var_1224, i_lay_var_1221) * absh2oc(ig_var_1216)
        IF (i_lay_var_1221 == i_laysolfr_var_1222(iplon_var_1224)) p_sfluxzen_var_1212(iplon_var_1224, ig_var_1216) = sfluxrefc_var_312(ig_var_1216)
        p_taur_var_1214(iplon_var_1224, i_lay_var_1221, ig_var_1216) = z_tauray_var_1225
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol29
SUBROUTINE srtm_taumol28(kidia_var_1226, kfdia_var_1227, klev_var_1228, p_fac00_var_1229, p_fac01_var_1230, p_fac10_var_1231, p_fac11_var_1232, k_jp_var_1233, k_jt_var_1234, k_jt1_var_1235, p_oneminus_var_1236, p_colmol_var_1237, p_colo2_var_1238, p_colo3_var_1239, k_laytrop_var_1240, p_sfluxzen_var_1241, p_taug_var_1242, p_taur_var_1243, prmu0_var_1244)
  USE yoesrta28, ONLY: absa_var_303, absb_var_304, layreffr_var_302, rayl_var_300, sfluxrefc_var_305, strrat_var_301
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1226, kfdia_var_1227
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1228
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1229(kidia_var_1226 : kfdia_var_1227, klev_var_1228)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1230(kidia_var_1226 : kfdia_var_1227, klev_var_1228)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1231(kidia_var_1226 : kfdia_var_1227, klev_var_1228)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1232(kidia_var_1226 : kfdia_var_1227, klev_var_1228)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1233(kidia_var_1226 : kfdia_var_1227, klev_var_1228)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1234(kidia_var_1226 : kfdia_var_1227, klev_var_1228)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1235(kidia_var_1226 : kfdia_var_1227, klev_var_1228)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1236(kidia_var_1226 : kfdia_var_1227)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1237(kidia_var_1226 : kfdia_var_1227, klev_var_1228)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_1238(kidia_var_1226 : kfdia_var_1227, klev_var_1228)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1239(kidia_var_1226 : kfdia_var_1227, klev_var_1228)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1240(kidia_var_1226 : kfdia_var_1227)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1241(kidia_var_1226 : kfdia_var_1227, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1242(kidia_var_1226 : kfdia_var_1227, klev_var_1228, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1243(kidia_var_1226 : kfdia_var_1227, klev_var_1228, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1244(kidia_var_1226 : kfdia_var_1227)
  INTEGER(KIND = 4) :: ig_var_1245, ind0_var_1246, ind1_var_1247, js_var_1248, i_lay_var_1249, i_laysolfr_var_1250(kidia_var_1226 : kfdia_var_1227), i_nlayers_var_1251, iplon_var_1252
  INTEGER(KIND = 4) :: laytrop_min_var_1253, laytrop_max_var_1254
  REAL(KIND = 8) :: z_fs_var_1255, z_speccomb_var_1256, z_specmult_var_1257, z_specparm_var_1258, z_tauray_var_1259
  laytrop_min_var_1253 = MINVAL(k_laytrop_var_1240(kidia_var_1226 : kfdia_var_1227))
  laytrop_max_var_1254 = MAXVAL(k_laytrop_var_1240(kidia_var_1226 : kfdia_var_1227))
  i_nlayers_var_1251 = klev_var_1228
  DO iplon_var_1252 = kidia_var_1226, kfdia_var_1227
    i_laysolfr_var_1250(iplon_var_1252) = i_nlayers_var_1251
  END DO
  DO i_lay_var_1249 = 1, laytrop_min_var_1253
    DO iplon_var_1252 = kidia_var_1226, kfdia_var_1227
      z_speccomb_var_1256 = p_colo3_var_1239(iplon_var_1252, i_lay_var_1249) + strrat_var_301 * p_colo2_var_1238(iplon_var_1252, i_lay_var_1249)
      z_specparm_var_1258 = p_colo3_var_1239(iplon_var_1252, i_lay_var_1249) / z_speccomb_var_1256
      z_specparm_var_1258 = MIN(p_oneminus_var_1236(iplon_var_1252), z_specparm_var_1258)
      z_specmult_var_1257 = 8.0D0 * (z_specparm_var_1258)
      js_var_1248 = 1 + INT(z_specmult_var_1257)
      z_fs_var_1255 = z_specmult_var_1257 - AINT(z_specmult_var_1257)
      ind0_var_1246 = ((k_jp_var_1233(iplon_var_1252, i_lay_var_1249) - 1) * 5 + (k_jt_var_1234(iplon_var_1252, i_lay_var_1249) - 1)) * nspa_var_313(28) + js_var_1248
      ind1_var_1247 = (k_jp_var_1233(iplon_var_1252, i_lay_var_1249) * 5 + (k_jt1_var_1235(iplon_var_1252, i_lay_var_1249) - 1)) * nspa_var_313(28) + js_var_1248
      z_tauray_var_1259 = p_colmol_var_1237(iplon_var_1252, i_lay_var_1249) * rayl_var_300
      DO ig_var_1245 = 1, 6
        p_taug_var_1242(iplon_var_1252, i_lay_var_1249, ig_var_1245) = z_speccomb_var_1256 * ((1.0D0 - z_fs_var_1255) * (absa_var_303(ind0_var_1246, ig_var_1245) * p_fac00_var_1229(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind0_var_1246 + 9, ig_var_1245) * p_fac10_var_1231(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind1_var_1247, ig_var_1245) * p_fac01_var_1230(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind1_var_1247 + 9, ig_var_1245) * p_fac11_var_1232(iplon_var_1252, i_lay_var_1249)) + z_fs_var_1255 * (absa_var_303(ind0_var_1246 + 1, ig_var_1245) * p_fac00_var_1229(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind0_var_1246 + 10, ig_var_1245) * p_fac10_var_1231(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind1_var_1247 + 1, ig_var_1245) * p_fac01_var_1230(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind1_var_1247 + 10, ig_var_1245) * p_fac11_var_1232(iplon_var_1252, i_lay_var_1249)))
        p_taur_var_1243(iplon_var_1252, i_lay_var_1249, ig_var_1245) = z_tauray_var_1259
      END DO
    END DO
  END DO
  DO i_lay_var_1249 = laytrop_min_var_1253 + 1, laytrop_max_var_1254
    DO iplon_var_1252 = kidia_var_1226, kfdia_var_1227
      IF (i_lay_var_1249 <= k_laytrop_var_1240(iplon_var_1252)) THEN
        z_speccomb_var_1256 = p_colo3_var_1239(iplon_var_1252, i_lay_var_1249) + strrat_var_301 * p_colo2_var_1238(iplon_var_1252, i_lay_var_1249)
        z_specparm_var_1258 = p_colo3_var_1239(iplon_var_1252, i_lay_var_1249) / z_speccomb_var_1256
        z_specparm_var_1258 = MIN(p_oneminus_var_1236(iplon_var_1252), z_specparm_var_1258)
        z_specmult_var_1257 = 8.0D0 * (z_specparm_var_1258)
        js_var_1248 = 1 + INT(z_specmult_var_1257)
        z_fs_var_1255 = z_specmult_var_1257 - AINT(z_specmult_var_1257)
        ind0_var_1246 = ((k_jp_var_1233(iplon_var_1252, i_lay_var_1249) - 1) * 5 + (k_jt_var_1234(iplon_var_1252, i_lay_var_1249) - 1)) * nspa_var_313(28) + js_var_1248
        ind1_var_1247 = (k_jp_var_1233(iplon_var_1252, i_lay_var_1249) * 5 + (k_jt1_var_1235(iplon_var_1252, i_lay_var_1249) - 1)) * nspa_var_313(28) + js_var_1248
        z_tauray_var_1259 = p_colmol_var_1237(iplon_var_1252, i_lay_var_1249) * rayl_var_300
        DO ig_var_1245 = 1, 6
          p_taug_var_1242(iplon_var_1252, i_lay_var_1249, ig_var_1245) = z_speccomb_var_1256 * ((1.0D0 - z_fs_var_1255) * (absa_var_303(ind0_var_1246, ig_var_1245) * p_fac00_var_1229(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind0_var_1246 + 9, ig_var_1245) * p_fac10_var_1231(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind1_var_1247, ig_var_1245) * p_fac01_var_1230(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind1_var_1247 + 9, ig_var_1245) * p_fac11_var_1232(iplon_var_1252, i_lay_var_1249)) + z_fs_var_1255 * (absa_var_303(ind0_var_1246 + 1, ig_var_1245) * p_fac00_var_1229(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind0_var_1246 + 10, ig_var_1245) * p_fac10_var_1231(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind1_var_1247 + 1, ig_var_1245) * p_fac01_var_1230(iplon_var_1252, i_lay_var_1249) + absa_var_303(ind1_var_1247 + 10, ig_var_1245) * p_fac11_var_1232(iplon_var_1252, i_lay_var_1249)))
          p_taur_var_1243(iplon_var_1252, i_lay_var_1249, ig_var_1245) = z_tauray_var_1259
        END DO
      ELSE
        IF (k_jp_var_1233(iplon_var_1252, i_lay_var_1249 - 1) < layreffr_var_302 .AND. k_jp_var_1233(iplon_var_1252, i_lay_var_1249) >= layreffr_var_302) i_laysolfr_var_1250(iplon_var_1252) = i_lay_var_1249
        z_speccomb_var_1256 = p_colo3_var_1239(iplon_var_1252, i_lay_var_1249) + strrat_var_301 * p_colo2_var_1238(iplon_var_1252, i_lay_var_1249)
        z_specparm_var_1258 = p_colo3_var_1239(iplon_var_1252, i_lay_var_1249) / z_speccomb_var_1256
        z_specparm_var_1258 = MIN(p_oneminus_var_1236(iplon_var_1252), z_specparm_var_1258)
        z_specmult_var_1257 = 4.0D0 * (z_specparm_var_1258)
        js_var_1248 = 1 + INT(z_specmult_var_1257)
        z_fs_var_1255 = z_specmult_var_1257 - AINT(z_specmult_var_1257)
        ind0_var_1246 = ((k_jp_var_1233(iplon_var_1252, i_lay_var_1249) - 13) * 5 + (k_jt_var_1234(iplon_var_1252, i_lay_var_1249) - 1)) * nspb_var_314(28) + js_var_1248
        ind1_var_1247 = ((k_jp_var_1233(iplon_var_1252, i_lay_var_1249) - 12) * 5 + (k_jt1_var_1235(iplon_var_1252, i_lay_var_1249) - 1)) * nspb_var_314(28) + js_var_1248
        z_tauray_var_1259 = p_colmol_var_1237(iplon_var_1252, i_lay_var_1249) * rayl_var_300
        DO ig_var_1245 = 1, 6
          p_taug_var_1242(iplon_var_1252, i_lay_var_1249, ig_var_1245) = z_speccomb_var_1256 * ((1.0D0 - z_fs_var_1255) * (absb_var_304(ind0_var_1246, ig_var_1245) * p_fac00_var_1229(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind0_var_1246 + 5, ig_var_1245) * p_fac10_var_1231(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind1_var_1247, ig_var_1245) * p_fac01_var_1230(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind1_var_1247 + 5, ig_var_1245) * p_fac11_var_1232(iplon_var_1252, i_lay_var_1249)) + z_fs_var_1255 * (absb_var_304(ind0_var_1246 + 1, ig_var_1245) * p_fac00_var_1229(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind0_var_1246 + 6, ig_var_1245) * p_fac10_var_1231(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind1_var_1247 + 1, ig_var_1245) * p_fac01_var_1230(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind1_var_1247 + 6, ig_var_1245) * p_fac11_var_1232(iplon_var_1252, i_lay_var_1249)))
          IF (i_lay_var_1249 == i_laysolfr_var_1250(iplon_var_1252)) p_sfluxzen_var_1241(iplon_var_1252, ig_var_1245) = sfluxrefc_var_305(ig_var_1245, js_var_1248) + z_fs_var_1255 * (sfluxrefc_var_305(ig_var_1245, js_var_1248 + 1) - sfluxrefc_var_305(ig_var_1245, js_var_1248))
          p_taur_var_1243(iplon_var_1252, i_lay_var_1249, ig_var_1245) = z_tauray_var_1259
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1249 = laytrop_max_var_1254 + 1, i_nlayers_var_1251
    DO iplon_var_1252 = kidia_var_1226, kfdia_var_1227
      IF (k_jp_var_1233(iplon_var_1252, i_lay_var_1249 - 1) < layreffr_var_302 .AND. k_jp_var_1233(iplon_var_1252, i_lay_var_1249) >= layreffr_var_302) i_laysolfr_var_1250(iplon_var_1252) = i_lay_var_1249
      z_speccomb_var_1256 = p_colo3_var_1239(iplon_var_1252, i_lay_var_1249) + strrat_var_301 * p_colo2_var_1238(iplon_var_1252, i_lay_var_1249)
      z_specparm_var_1258 = p_colo3_var_1239(iplon_var_1252, i_lay_var_1249) / z_speccomb_var_1256
      z_specparm_var_1258 = MIN(p_oneminus_var_1236(iplon_var_1252), z_specparm_var_1258)
      z_specmult_var_1257 = 4.0D0 * (z_specparm_var_1258)
      js_var_1248 = 1 + INT(z_specmult_var_1257)
      z_fs_var_1255 = z_specmult_var_1257 - AINT(z_specmult_var_1257)
      ind0_var_1246 = ((k_jp_var_1233(iplon_var_1252, i_lay_var_1249) - 13) * 5 + (k_jt_var_1234(iplon_var_1252, i_lay_var_1249) - 1)) * nspb_var_314(28) + js_var_1248
      ind1_var_1247 = ((k_jp_var_1233(iplon_var_1252, i_lay_var_1249) - 12) * 5 + (k_jt1_var_1235(iplon_var_1252, i_lay_var_1249) - 1)) * nspb_var_314(28) + js_var_1248
      z_tauray_var_1259 = p_colmol_var_1237(iplon_var_1252, i_lay_var_1249) * rayl_var_300
      DO ig_var_1245 = 1, 6
        p_taug_var_1242(iplon_var_1252, i_lay_var_1249, ig_var_1245) = z_speccomb_var_1256 * ((1.0D0 - z_fs_var_1255) * (absb_var_304(ind0_var_1246, ig_var_1245) * p_fac00_var_1229(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind0_var_1246 + 5, ig_var_1245) * p_fac10_var_1231(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind1_var_1247, ig_var_1245) * p_fac01_var_1230(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind1_var_1247 + 5, ig_var_1245) * p_fac11_var_1232(iplon_var_1252, i_lay_var_1249)) + z_fs_var_1255 * (absb_var_304(ind0_var_1246 + 1, ig_var_1245) * p_fac00_var_1229(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind0_var_1246 + 6, ig_var_1245) * p_fac10_var_1231(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind1_var_1247 + 1, ig_var_1245) * p_fac01_var_1230(iplon_var_1252, i_lay_var_1249) + absb_var_304(ind1_var_1247 + 6, ig_var_1245) * p_fac11_var_1232(iplon_var_1252, i_lay_var_1249)))
        IF (i_lay_var_1249 == i_laysolfr_var_1250(iplon_var_1252)) p_sfluxzen_var_1241(iplon_var_1252, ig_var_1245) = sfluxrefc_var_305(ig_var_1245, js_var_1248) + z_fs_var_1255 * (sfluxrefc_var_305(ig_var_1245, js_var_1248 + 1) - sfluxrefc_var_305(ig_var_1245, js_var_1248))
        p_taur_var_1243(iplon_var_1252, i_lay_var_1249, ig_var_1245) = z_tauray_var_1259
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol28
SUBROUTINE srtm_taumol17(kidia_var_1260, kfdia_var_1261, klev_var_1262, p_fac00_var_1263, p_fac01_var_1264, p_fac10_var_1265, p_fac11_var_1266, k_jp_var_1267, k_jt_var_1268, k_jt1_var_1269, p_oneminus_var_1270, p_colh2o_var_1271, p_colco2_var_1272, p_colmol_var_1273, k_laytrop_var_1274, p_selffac_var_1275, p_selffrac_var_1276, k_indself_var_1277, p_forfac_var_1278, p_forfrac_var_1279, k_indfor_var_1280, p_sfluxzen_var_1281, p_taug_var_1282, p_taur_var_1283, prmu0_var_1284)
  USE yoesrta17, ONLY: absa_var_228, absb_var_229, forrefc_var_231, layreffr_var_227, rayl_var_225, selfrefc_var_230, sfluxrefc_var_232, strrat_var_226
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1260, kfdia_var_1261
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1262
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1263(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1264(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1265(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1266(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1267(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1268(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1269(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1270(kidia_var_1260 : kfdia_var_1261)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1271(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1272(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1273(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1274(kidia_var_1260 : kfdia_var_1261)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1275(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1276(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1277(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1278(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1279(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1280(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1281(kidia_var_1260 : kfdia_var_1261, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1282(kidia_var_1260 : kfdia_var_1261, klev_var_1262, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1283(kidia_var_1260 : kfdia_var_1261, klev_var_1262, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1284(kidia_var_1260 : kfdia_var_1261)
  INTEGER(KIND = 4) :: ig_var_1285, ind0_var_1286, ind1_var_1287, inds_var_1288, indf_var_1289, js_var_1290, i_lay_var_1291, i_laysolfr_var_1292(kidia_var_1260 : kfdia_var_1261), i_nlayers_var_1293, iplon_var_1294
  INTEGER(KIND = 4) :: laytrop_min_var_1295, laytrop_max_var_1296
  REAL(KIND = 8) :: z_fs_var_1297, z_speccomb_var_1298, z_specmult_var_1299, z_specparm_var_1300, z_tauray_var_1301
  laytrop_min_var_1295 = MINVAL(k_laytrop_var_1274(kidia_var_1260 : kfdia_var_1261))
  laytrop_max_var_1296 = MAXVAL(k_laytrop_var_1274(kidia_var_1260 : kfdia_var_1261))
  i_nlayers_var_1293 = klev_var_1262
  DO iplon_var_1294 = kidia_var_1260, kfdia_var_1261
    i_laysolfr_var_1292(iplon_var_1294) = i_nlayers_var_1293
  END DO
  DO i_lay_var_1291 = 1, laytrop_min_var_1295
    DO iplon_var_1294 = kidia_var_1260, kfdia_var_1261
      z_speccomb_var_1298 = p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) + strrat_var_226 * p_colco2_var_1272(iplon_var_1294, i_lay_var_1291)
      z_specparm_var_1300 = p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) / z_speccomb_var_1298
      z_specparm_var_1300 = MIN(p_oneminus_var_1270(iplon_var_1294), z_specparm_var_1300)
      z_specmult_var_1299 = 8.0D0 * (z_specparm_var_1300)
      js_var_1290 = 1 + INT(z_specmult_var_1299)
      z_fs_var_1297 = z_specmult_var_1299 - AINT(z_specmult_var_1299)
      ind0_var_1286 = ((k_jp_var_1267(iplon_var_1294, i_lay_var_1291) - 1) * 5 + (k_jt_var_1268(iplon_var_1294, i_lay_var_1291) - 1)) * nspa_var_313(17) + js_var_1290
      ind1_var_1287 = (k_jp_var_1267(iplon_var_1294, i_lay_var_1291) * 5 + (k_jt1_var_1269(iplon_var_1294, i_lay_var_1291) - 1)) * nspa_var_313(17) + js_var_1290
      inds_var_1288 = k_indself_var_1277(iplon_var_1294, i_lay_var_1291)
      indf_var_1289 = k_indfor_var_1280(iplon_var_1294, i_lay_var_1291)
      z_tauray_var_1301 = p_colmol_var_1273(iplon_var_1294, i_lay_var_1291) * rayl_var_225
      DO ig_var_1285 = 1, 12
        p_taug_var_1282(iplon_var_1294, i_lay_var_1291, ig_var_1285) = z_speccomb_var_1298 * ((1.0D0 - z_fs_var_1297) * (absa_var_228(ind0_var_1286, ig_var_1285) * p_fac00_var_1263(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind0_var_1286 + 9, ig_var_1285) * p_fac10_var_1265(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind1_var_1287, ig_var_1285) * p_fac01_var_1264(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind1_var_1287 + 9, ig_var_1285) * p_fac11_var_1266(iplon_var_1294, i_lay_var_1291)) + z_fs_var_1297 * (absa_var_228(ind0_var_1286 + 1, ig_var_1285) * p_fac00_var_1263(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind0_var_1286 + 10, ig_var_1285) * p_fac10_var_1265(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind1_var_1287 + 1, ig_var_1285) * p_fac01_var_1264(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind1_var_1287 + 10, ig_var_1285) * p_fac11_var_1266(iplon_var_1294, i_lay_var_1291))) + p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) * (p_selffac_var_1275(iplon_var_1294, i_lay_var_1291) * (selfrefc_var_230(inds_var_1288, ig_var_1285) + p_selffrac_var_1276(iplon_var_1294, i_lay_var_1291) * (selfrefc_var_230(inds_var_1288 + 1, ig_var_1285) - selfrefc_var_230(inds_var_1288, ig_var_1285))) + p_forfac_var_1278(iplon_var_1294, i_lay_var_1291) * (forrefc_var_231(indf_var_1289, ig_var_1285) + p_forfrac_var_1279(iplon_var_1294, i_lay_var_1291) * (forrefc_var_231(indf_var_1289 + 1, ig_var_1285) - forrefc_var_231(indf_var_1289, ig_var_1285))))
        p_taur_var_1283(iplon_var_1294, i_lay_var_1291, ig_var_1285) = z_tauray_var_1301
      END DO
    END DO
  END DO
  DO i_lay_var_1291 = laytrop_min_var_1295 + 1, laytrop_max_var_1296
    DO iplon_var_1294 = kidia_var_1260, kfdia_var_1261
      IF (i_lay_var_1291 <= k_laytrop_var_1274(iplon_var_1294)) THEN
        z_speccomb_var_1298 = p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) + strrat_var_226 * p_colco2_var_1272(iplon_var_1294, i_lay_var_1291)
        z_specparm_var_1300 = p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) / z_speccomb_var_1298
        z_specparm_var_1300 = MIN(p_oneminus_var_1270(iplon_var_1294), z_specparm_var_1300)
        z_specmult_var_1299 = 8.0D0 * (z_specparm_var_1300)
        js_var_1290 = 1 + INT(z_specmult_var_1299)
        z_fs_var_1297 = z_specmult_var_1299 - AINT(z_specmult_var_1299)
        ind0_var_1286 = ((k_jp_var_1267(iplon_var_1294, i_lay_var_1291) - 1) * 5 + (k_jt_var_1268(iplon_var_1294, i_lay_var_1291) - 1)) * nspa_var_313(17) + js_var_1290
        ind1_var_1287 = (k_jp_var_1267(iplon_var_1294, i_lay_var_1291) * 5 + (k_jt1_var_1269(iplon_var_1294, i_lay_var_1291) - 1)) * nspa_var_313(17) + js_var_1290
        inds_var_1288 = k_indself_var_1277(iplon_var_1294, i_lay_var_1291)
        indf_var_1289 = k_indfor_var_1280(iplon_var_1294, i_lay_var_1291)
        z_tauray_var_1301 = p_colmol_var_1273(iplon_var_1294, i_lay_var_1291) * rayl_var_225
        DO ig_var_1285 = 1, 12
          p_taug_var_1282(iplon_var_1294, i_lay_var_1291, ig_var_1285) = z_speccomb_var_1298 * ((1.0D0 - z_fs_var_1297) * (absa_var_228(ind0_var_1286, ig_var_1285) * p_fac00_var_1263(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind0_var_1286 + 9, ig_var_1285) * p_fac10_var_1265(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind1_var_1287, ig_var_1285) * p_fac01_var_1264(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind1_var_1287 + 9, ig_var_1285) * p_fac11_var_1266(iplon_var_1294, i_lay_var_1291)) + z_fs_var_1297 * (absa_var_228(ind0_var_1286 + 1, ig_var_1285) * p_fac00_var_1263(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind0_var_1286 + 10, ig_var_1285) * p_fac10_var_1265(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind1_var_1287 + 1, ig_var_1285) * p_fac01_var_1264(iplon_var_1294, i_lay_var_1291) + absa_var_228(ind1_var_1287 + 10, ig_var_1285) * p_fac11_var_1266(iplon_var_1294, i_lay_var_1291))) + p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) * (p_selffac_var_1275(iplon_var_1294, i_lay_var_1291) * (selfrefc_var_230(inds_var_1288, ig_var_1285) + p_selffrac_var_1276(iplon_var_1294, i_lay_var_1291) * (selfrefc_var_230(inds_var_1288 + 1, ig_var_1285) - selfrefc_var_230(inds_var_1288, ig_var_1285))) + p_forfac_var_1278(iplon_var_1294, i_lay_var_1291) * (forrefc_var_231(indf_var_1289, ig_var_1285) + p_forfrac_var_1279(iplon_var_1294, i_lay_var_1291) * (forrefc_var_231(indf_var_1289 + 1, ig_var_1285) - forrefc_var_231(indf_var_1289, ig_var_1285))))
          p_taur_var_1283(iplon_var_1294, i_lay_var_1291, ig_var_1285) = z_tauray_var_1301
        END DO
      ELSE
        IF (k_jp_var_1267(iplon_var_1294, i_lay_var_1291 - 1) < layreffr_var_227 .AND. k_jp_var_1267(iplon_var_1294, i_lay_var_1291) >= layreffr_var_227) i_laysolfr_var_1292(iplon_var_1294) = i_lay_var_1291
        z_speccomb_var_1298 = p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) + strrat_var_226 * p_colco2_var_1272(iplon_var_1294, i_lay_var_1291)
        z_specparm_var_1300 = p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) / z_speccomb_var_1298
        z_specparm_var_1300 = MIN(p_oneminus_var_1270(iplon_var_1294), z_specparm_var_1300)
        z_specmult_var_1299 = 4.0D0 * (z_specparm_var_1300)
        js_var_1290 = 1 + INT(z_specmult_var_1299)
        z_fs_var_1297 = z_specmult_var_1299 - AINT(z_specmult_var_1299)
        ind0_var_1286 = ((k_jp_var_1267(iplon_var_1294, i_lay_var_1291) - 13) * 5 + (k_jt_var_1268(iplon_var_1294, i_lay_var_1291) - 1)) * nspb_var_314(17) + js_var_1290
        ind1_var_1287 = ((k_jp_var_1267(iplon_var_1294, i_lay_var_1291) - 12) * 5 + (k_jt1_var_1269(iplon_var_1294, i_lay_var_1291) - 1)) * nspb_var_314(17) + js_var_1290
        indf_var_1289 = k_indfor_var_1280(iplon_var_1294, i_lay_var_1291)
        z_tauray_var_1301 = p_colmol_var_1273(iplon_var_1294, i_lay_var_1291) * rayl_var_225
        DO ig_var_1285 = 1, 12
          p_taug_var_1282(iplon_var_1294, i_lay_var_1291, ig_var_1285) = z_speccomb_var_1298 * ((1.0D0 - z_fs_var_1297) * (absb_var_229(ind0_var_1286, ig_var_1285) * p_fac00_var_1263(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind0_var_1286 + 5, ig_var_1285) * p_fac10_var_1265(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind1_var_1287, ig_var_1285) * p_fac01_var_1264(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind1_var_1287 + 5, ig_var_1285) * p_fac11_var_1266(iplon_var_1294, i_lay_var_1291)) + z_fs_var_1297 * (absb_var_229(ind0_var_1286 + 1, ig_var_1285) * p_fac00_var_1263(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind0_var_1286 + 6, ig_var_1285) * p_fac10_var_1265(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind1_var_1287 + 1, ig_var_1285) * p_fac01_var_1264(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind1_var_1287 + 6, ig_var_1285) * p_fac11_var_1266(iplon_var_1294, i_lay_var_1291))) + p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) * p_forfac_var_1278(iplon_var_1294, i_lay_var_1291) * (forrefc_var_231(indf_var_1289, ig_var_1285) + p_forfrac_var_1279(iplon_var_1294, i_lay_var_1291) * (forrefc_var_231(indf_var_1289 + 1, ig_var_1285) - forrefc_var_231(indf_var_1289, ig_var_1285)))
          IF (i_lay_var_1291 == i_laysolfr_var_1292(iplon_var_1294)) p_sfluxzen_var_1281(iplon_var_1294, ig_var_1285) = sfluxrefc_var_232(ig_var_1285, js_var_1290) + z_fs_var_1297 * (sfluxrefc_var_232(ig_var_1285, js_var_1290 + 1) - sfluxrefc_var_232(ig_var_1285, js_var_1290))
          p_taur_var_1283(iplon_var_1294, i_lay_var_1291, ig_var_1285) = z_tauray_var_1301
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1291 = laytrop_max_var_1296 + 1, i_nlayers_var_1293
    DO iplon_var_1294 = kidia_var_1260, kfdia_var_1261
      IF (k_jp_var_1267(iplon_var_1294, i_lay_var_1291 - 1) < layreffr_var_227 .AND. k_jp_var_1267(iplon_var_1294, i_lay_var_1291) >= layreffr_var_227) i_laysolfr_var_1292(iplon_var_1294) = i_lay_var_1291
      z_speccomb_var_1298 = p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) + strrat_var_226 * p_colco2_var_1272(iplon_var_1294, i_lay_var_1291)
      z_specparm_var_1300 = p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) / z_speccomb_var_1298
      z_specparm_var_1300 = MIN(p_oneminus_var_1270(iplon_var_1294), z_specparm_var_1300)
      z_specmult_var_1299 = 4.0D0 * (z_specparm_var_1300)
      js_var_1290 = 1 + INT(z_specmult_var_1299)
      z_fs_var_1297 = z_specmult_var_1299 - AINT(z_specmult_var_1299)
      ind0_var_1286 = ((k_jp_var_1267(iplon_var_1294, i_lay_var_1291) - 13) * 5 + (k_jt_var_1268(iplon_var_1294, i_lay_var_1291) - 1)) * nspb_var_314(17) + js_var_1290
      ind1_var_1287 = ((k_jp_var_1267(iplon_var_1294, i_lay_var_1291) - 12) * 5 + (k_jt1_var_1269(iplon_var_1294, i_lay_var_1291) - 1)) * nspb_var_314(17) + js_var_1290
      indf_var_1289 = k_indfor_var_1280(iplon_var_1294, i_lay_var_1291)
      z_tauray_var_1301 = p_colmol_var_1273(iplon_var_1294, i_lay_var_1291) * rayl_var_225
      DO ig_var_1285 = 1, 12
        p_taug_var_1282(iplon_var_1294, i_lay_var_1291, ig_var_1285) = z_speccomb_var_1298 * ((1.0D0 - z_fs_var_1297) * (absb_var_229(ind0_var_1286, ig_var_1285) * p_fac00_var_1263(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind0_var_1286 + 5, ig_var_1285) * p_fac10_var_1265(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind1_var_1287, ig_var_1285) * p_fac01_var_1264(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind1_var_1287 + 5, ig_var_1285) * p_fac11_var_1266(iplon_var_1294, i_lay_var_1291)) + z_fs_var_1297 * (absb_var_229(ind0_var_1286 + 1, ig_var_1285) * p_fac00_var_1263(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind0_var_1286 + 6, ig_var_1285) * p_fac10_var_1265(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind1_var_1287 + 1, ig_var_1285) * p_fac01_var_1264(iplon_var_1294, i_lay_var_1291) + absb_var_229(ind1_var_1287 + 6, ig_var_1285) * p_fac11_var_1266(iplon_var_1294, i_lay_var_1291))) + p_colh2o_var_1271(iplon_var_1294, i_lay_var_1291) * p_forfac_var_1278(iplon_var_1294, i_lay_var_1291) * (forrefc_var_231(indf_var_1289, ig_var_1285) + p_forfrac_var_1279(iplon_var_1294, i_lay_var_1291) * (forrefc_var_231(indf_var_1289 + 1, ig_var_1285) - forrefc_var_231(indf_var_1289, ig_var_1285)))
        IF (i_lay_var_1291 == i_laysolfr_var_1292(iplon_var_1294)) p_sfluxzen_var_1281(iplon_var_1294, ig_var_1285) = sfluxrefc_var_232(ig_var_1285, js_var_1290) + z_fs_var_1297 * (sfluxrefc_var_232(ig_var_1285, js_var_1290 + 1) - sfluxrefc_var_232(ig_var_1285, js_var_1290))
        p_taur_var_1283(iplon_var_1294, i_lay_var_1291, ig_var_1285) = z_tauray_var_1301
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol17
SUBROUTINE srtm_setcoef(kidia_var_1302, kfdia_var_1303, klev_var_1304, pavel_var_1305, ptavel_var_1306, pcoldry_var_1307, pwkl_var_1308, klaytrop_var_1309, pcolch4_var_1310, pcolco2_var_1311, pcolh2o_var_1312, pcolmol_var_1313, pcolo2_var_1314, pcolo3_var_1315, pforfac_var_1316, pforfrac_var_1317, kindfor_var_1318, pselffac_var_1319, pselffrac_var_1320, kindself_var_1321, pfac00_var_1322, pfac01_var_1323, pfac10_var_1324, pfac11_var_1325, kjp_var_1326, kjt_var_1327, kjt1_var_1328, prmu0_var_1329)
  USE yoesrtwn, ONLY: preflog_var_315, tref_var_316
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1302, kfdia_var_1303
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1304
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1305(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(IN) :: ptavel_var_1306(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(IN) :: pcoldry_var_1307(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(IN) :: pwkl_var_1308(kidia_var_1302 : kfdia_var_1303, 35, klev_var_1304)
  INTEGER(KIND = 4), INTENT(INOUT) :: klaytrop_var_1309(kidia_var_1302 : kfdia_var_1303)
  REAL(KIND = 8), INTENT(INOUT) :: pcolch4_var_1310(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pcolco2_var_1311(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pcolh2o_var_1312(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pcolmol_var_1313(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pcolo2_var_1314(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pcolo3_var_1315(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pforfac_var_1316(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pforfrac_var_1317(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  INTEGER(KIND = 4), INTENT(INOUT) :: kindfor_var_1318(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pselffac_var_1319(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pselffrac_var_1320(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  INTEGER(KIND = 4), INTENT(INOUT) :: kindself_var_1321(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pfac00_var_1322(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pfac01_var_1323(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pfac10_var_1324(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(INOUT) :: pfac11_var_1325(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjp_var_1326(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjt_var_1327(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjt1_var_1328(kidia_var_1302 : kfdia_var_1303, klev_var_1304)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1329(kidia_var_1302 : kfdia_var_1303)
  INTEGER(KIND = 4) :: i_nlayers_var_1330, jk_var_1331, jl_var_1332, jp1_var_1333
  REAL(KIND = 8) :: z_stpfac_var_1334, z_plog_var_1335
  REAL(KIND = 8) :: z_fp_var_1336, z_ft_var_1337, z_ft1_var_1338, z_water_var_1339, z_scalefac_var_1340
  REAL(KIND = 8) :: z_factor_var_1341, z_co2reg_var_1342, z_compfp_var_1343
  LOGICAL :: goto_0 = .FALSE.
  LOGICAL :: goto_1 = .FALSE.
  z_stpfac_var_1334 = 0.29220138203356366D0
  i_nlayers_var_1330 = klev_var_1304
  DO jk_var_1331 = 1, klev_var_1304
    DO jl_var_1332 = kidia_var_1302, kfdia_var_1303
      pcolmol_var_1313(jl_var_1332, jk_var_1331) = 0.0D0
    END DO
  END DO
  DO jl_var_1332 = kidia_var_1302, kfdia_var_1303
    IF (prmu0_var_1329(jl_var_1332) > 0.0D0) THEN
      klaytrop_var_1309(jl_var_1332) = 0
    END IF
  END DO
  DO jk_var_1331 = 1, i_nlayers_var_1330
    DO jl_var_1332 = kidia_var_1302, kfdia_var_1303
      IF (prmu0_var_1329(jl_var_1332) > 0.0D0) THEN
        z_plog_var_1335 = LOG(pavel_var_1305(jl_var_1332, jk_var_1331))
        kjp_var_1326(jl_var_1332, jk_var_1331) = INT(36.0D0 - 5.0D0 * (z_plog_var_1335 + 0.04D0))
        IF (kjp_var_1326(jl_var_1332, jk_var_1331) < 1) THEN
          kjp_var_1326(jl_var_1332, jk_var_1331) = 1
        ELSE IF (kjp_var_1326(jl_var_1332, jk_var_1331) > 58) THEN
          kjp_var_1326(jl_var_1332, jk_var_1331) = 58
        END IF
        jp1_var_1333 = kjp_var_1326(jl_var_1332, jk_var_1331) + 1
        z_fp_var_1336 = 5.0 * (preflog_var_315(kjp_var_1326(jl_var_1332, jk_var_1331)) - z_plog_var_1335)
        kjt_var_1327(jl_var_1332, jk_var_1331) = INT(3.0 + (ptavel_var_1306(jl_var_1332, jk_var_1331) - tref_var_316(kjp_var_1326(jl_var_1332, jk_var_1331))) / 15.0)
        IF (kjt_var_1327(jl_var_1332, jk_var_1331) < 1) THEN
          kjt_var_1327(jl_var_1332, jk_var_1331) = 1
        ELSE IF (kjt_var_1327(jl_var_1332, jk_var_1331) > 4) THEN
          kjt_var_1327(jl_var_1332, jk_var_1331) = 4
        END IF
        z_ft_var_1337 = ((ptavel_var_1306(jl_var_1332, jk_var_1331) - tref_var_316(kjp_var_1326(jl_var_1332, jk_var_1331))) / 15.0) - REAL(kjt_var_1327(jl_var_1332, jk_var_1331) - 3)
        kjt1_var_1328(jl_var_1332, jk_var_1331) = INT(3.0 + (ptavel_var_1306(jl_var_1332, jk_var_1331) - tref_var_316(jp1_var_1333)) / 15.0)
        IF (kjt1_var_1328(jl_var_1332, jk_var_1331) < 1) THEN
          kjt1_var_1328(jl_var_1332, jk_var_1331) = 1
        ELSE IF (kjt1_var_1328(jl_var_1332, jk_var_1331) > 4) THEN
          kjt1_var_1328(jl_var_1332, jk_var_1331) = 4
        END IF
        z_ft1_var_1338 = ((ptavel_var_1306(jl_var_1332, jk_var_1331) - tref_var_316(jp1_var_1333)) / 15.0) - REAL(kjt1_var_1328(jl_var_1332, jk_var_1331) - 3)
        z_water_var_1339 = pwkl_var_1308(jl_var_1332, 1, jk_var_1331) / pcoldry_var_1307(jl_var_1332, jk_var_1331)
        z_scalefac_var_1340 = pavel_var_1305(jl_var_1332, jk_var_1331) * z_stpfac_var_1334 / ptavel_var_1306(jl_var_1332, jk_var_1331)
        IF (z_plog_var_1335 <= 4.56D0) goto_0 = .TRUE.
        IF (.NOT. (goto_0)) klaytrop_var_1309(jl_var_1332) = klaytrop_var_1309(jl_var_1332) + 1
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pforfac_var_1316(jl_var_1332, jk_var_1331) = z_scalefac_var_1340 / (1.0 + z_water_var_1339)
        IF (.NOT. (goto_0)) z_factor_var_1341 = (332.0 - ptavel_var_1306(jl_var_1332, jk_var_1331)) / 36.0
        IF (.NOT. (goto_0)) kindfor_var_1318(jl_var_1332, jk_var_1331) = MIN(2, MAX(1, INT(z_factor_var_1341)))
        IF (.NOT. (goto_0)) pforfrac_var_1317(jl_var_1332, jk_var_1331) = z_factor_var_1341 - REAL(kindfor_var_1318(jl_var_1332, jk_var_1331))
        IF (.NOT. (goto_0)) pselffac_var_1319(jl_var_1332, jk_var_1331) = z_water_var_1339 * pforfac_var_1316(jl_var_1332, jk_var_1331)
        IF (.NOT. (goto_0)) z_factor_var_1341 = (ptavel_var_1306(jl_var_1332, jk_var_1331) - 188.0) / 7.2
        IF (.NOT. (goto_0)) kindself_var_1321(jl_var_1332, jk_var_1331) = MIN(9, MAX(1, INT(z_factor_var_1341) - 7))
        IF (.NOT. (goto_0)) pselffrac_var_1320(jl_var_1332, jk_var_1331) = z_factor_var_1341 - REAL(kindself_var_1321(jl_var_1332, jk_var_1331) + 7)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolh2o_var_1312(jl_var_1332, jk_var_1331) = 1E-20 * pwkl_var_1308(jl_var_1332, 1, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolco2_var_1311(jl_var_1332, jk_var_1331) = 1E-20 * pwkl_var_1308(jl_var_1332, 2, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolo3_var_1315(jl_var_1332, jk_var_1331) = 1E-20 * pwkl_var_1308(jl_var_1332, 3, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolch4_var_1310(jl_var_1332, jk_var_1331) = 1E-20 * pwkl_var_1308(jl_var_1332, 6, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolo2_var_1314(jl_var_1332, jk_var_1331) = 1E-20 * pwkl_var_1308(jl_var_1332, 7, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolmol_var_1313(jl_var_1332, jk_var_1331) = 1E-20 * pcoldry_var_1307(jl_var_1332, jk_var_1331) + pcolh2o_var_1312(jl_var_1332, jk_var_1331)
        IF (.NOT. (goto_0) .AND. pcolco2_var_1311(jl_var_1332, jk_var_1331) == 0.0) pcolco2_var_1311(jl_var_1332, jk_var_1331) = 1E-32 * pcoldry_var_1307(jl_var_1332, jk_var_1331)
        IF (.NOT. (goto_0) .AND. pcolch4_var_1310(jl_var_1332, jk_var_1331) == 0.0) pcolch4_var_1310(jl_var_1332, jk_var_1331) = 1E-32 * pcoldry_var_1307(jl_var_1332, jk_var_1331)
        IF (.NOT. (goto_0) .AND. pcolo2_var_1314(jl_var_1332, jk_var_1331) == 0.0) pcolo2_var_1314(jl_var_1332, jk_var_1331) = 1E-32 * pcoldry_var_1307(jl_var_1332, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) z_co2reg_var_1342 = 3.55E-24 * pcoldry_var_1307(jl_var_1332, jk_var_1331)
        IF (.NOT. (goto_0)) goto_1 = .TRUE.
5300    CONTINUE
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pforfac_var_1316(jl_var_1332, jk_var_1331) = z_scalefac_var_1340 / (1.0 + z_water_var_1339)
        IF (.NOT. (goto_1)) z_factor_var_1341 = (ptavel_var_1306(jl_var_1332, jk_var_1331) - 188.0) / 36.0
        IF (.NOT. (goto_1)) kindfor_var_1318(jl_var_1332, jk_var_1331) = 3
        IF (.NOT. (goto_1)) pforfrac_var_1317(jl_var_1332, jk_var_1331) = z_factor_var_1341 - 1.0
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolh2o_var_1312(jl_var_1332, jk_var_1331) = 1E-20 * pwkl_var_1308(jl_var_1332, 1, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolco2_var_1311(jl_var_1332, jk_var_1331) = 1E-20 * pwkl_var_1308(jl_var_1332, 2, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolo3_var_1315(jl_var_1332, jk_var_1331) = 1E-20 * pwkl_var_1308(jl_var_1332, 3, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolch4_var_1310(jl_var_1332, jk_var_1331) = 1E-20 * pwkl_var_1308(jl_var_1332, 6, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolo2_var_1314(jl_var_1332, jk_var_1331) = 1E-20 * pwkl_var_1308(jl_var_1332, 7, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolmol_var_1313(jl_var_1332, jk_var_1331) = 1E-20 * pcoldry_var_1307(jl_var_1332, jk_var_1331) + pcolh2o_var_1312(jl_var_1332, jk_var_1331)
        IF (.NOT. (goto_1) .AND. pcolco2_var_1311(jl_var_1332, jk_var_1331) == 0.0) pcolco2_var_1311(jl_var_1332, jk_var_1331) = 1E-32 * pcoldry_var_1307(jl_var_1332, jk_var_1331)
        IF (.NOT. (goto_1) .AND. pcolch4_var_1310(jl_var_1332, jk_var_1331) == 0.0) pcolch4_var_1310(jl_var_1332, jk_var_1331) = 1E-32 * pcoldry_var_1307(jl_var_1332, jk_var_1331)
        IF (.NOT. (goto_1) .AND. pcolo2_var_1314(jl_var_1332, jk_var_1331) == 0.0) pcolo2_var_1314(jl_var_1332, jk_var_1331) = 1E-32 * pcoldry_var_1307(jl_var_1332, jk_var_1331)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) z_co2reg_var_1342 = 3.55E-24 * pcoldry_var_1307(jl_var_1332, jk_var_1331)
        IF (.NOT. (goto_1)) pselffac_var_1319(jl_var_1332, jk_var_1331) = 0.0D0
        IF (.NOT. (goto_1)) pselffrac_var_1320(jl_var_1332, jk_var_1331) = 0.0D0
        IF (.NOT. (goto_1)) kindself_var_1321(jl_var_1332, jk_var_1331) = 0
5400    CONTINUE
        z_compfp_var_1343 = 1.0 - z_fp_var_1336
        pfac10_var_1324(jl_var_1332, jk_var_1331) = z_compfp_var_1343 * z_ft_var_1337
        pfac00_var_1322(jl_var_1332, jk_var_1331) = z_compfp_var_1343 * (1.0 - z_ft_var_1337)
        pfac11_var_1325(jl_var_1332, jk_var_1331) = z_fp_var_1336 * z_ft1_var_1338
        pfac01_var_1323(jl_var_1332, jk_var_1331) = z_fp_var_1336 * (1.0 - z_ft1_var_1338)
9000    FORMAT(1X, 2I3, 3I4, F6.1, 4F7.2, 12E9.2, 2I5)
      END IF
    END DO
  END DO
END SUBROUTINE srtm_setcoef
SUBROUTINE rrtm_taumol5(kidia_var_1344, kfdia_var_1345, klev_var_1346, taug_var_1347, wx_var_1348, p_tauaerl_var_1349, fac00_var_1350, fac01_var_1351, fac10_var_1352, fac11_var_1353, forfac_var_1372, forfrac_var_1371, indfor_var_1370, jp_var_1354, jt_var_1355, jt1_var_1356, oneminus_var_1357, colh2o_var_1358, colco2_var_1359, colo3_var_1360, laytrop_var_1361, selffac_var_1362, selffrac_var_1363, indself_var_1364, fracs_var_1365, rat_h2oco2_var_1366, rat_h2oco2_1_var_1367, rat_o3co2_var_1368, rat_o3co2_1_var_1369, minorfrac_var_1373, indminor_var_1374)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng5
  USE yoerrta5, ONLY: absa_var_175, absb_var_176, ccl4, forref_var_179, fracrefa_var_173, fracrefb_var_174, ka_mo3_var_177, selfref_var_178
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1344
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1345
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1346
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1347(kidia_var_1344 : kfdia_var_1345, 140, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: wx_var_1348(kidia_var_1344 : kfdia_var_1345, 4, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1349(kidia_var_1344 : kfdia_var_1345, klev_var_1346, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1350(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1351(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1352(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1353(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1354(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1355(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1356(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1357
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1358(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1359(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1360(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1361(kidia_var_1344 : kfdia_var_1345)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1362(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1363(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1364(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1365(kidia_var_1344 : kfdia_var_1345, 140, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1366(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1367(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_var_1368(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_1_var_1369(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1370(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1371(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1372(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1373(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1374(kidia_var_1344 : kfdia_var_1345, klev_var_1346)
  REAL(KIND = 8) :: speccomb_var_1375, speccomb1_var_1376, speccomb_mo3, speccomb_planck_var_1377
  INTEGER(KIND = 4) :: ind0_var_1378, ind1_var_1379, inds_var_1380, indf_var_1381, indm_var_1382
  INTEGER(KIND = 4) :: ig_var_1383, js_var_1384, lay_var_1385, js1_var_1386, jpl_var_1387, jmo3
  REAL(KIND = 8) :: refrat_planck_a_var_1388, refrat_planck_b_var_1389, refrat_m_a_var_1390
  REAL(KIND = 8) :: fac000_var_1391, fac100_var_1392, fac200_var_1393, fac010_var_1394, fac110_var_1395, fac210_var_1396, fac001_var_1397, fac101_var_1398, fac201_var_1399, fac011_var_1400, fac111_var_1401, fac211_var_1402
  REAL(KIND = 8) :: p_var_1403, p4_var_1404, fk0_var_1405, fk1_var_1406, fk2_var_1407
  REAL(KIND = 8) :: taufor_var_1408, tauself_var_1409, tau_major_var_1410(16), tau_major1_var_1411(16), o3m1, o3m2, abso3_var_1412
  REAL(KIND = 8) :: fs_var_1413, specmult_var_1414, specparm_var_1415, fs1_var_1416, specmult1_var_1417, specparm1_var_1418, fpl_var_1419, specmult_planck_var_1420, specparm_planck_var_1421, fmo3, specmult_mo3, specparm_mo3
  INTEGER(KIND = 4) :: laytrop_min_var_1422, laytrop_max_var_1423
  INTEGER(KIND = 4) :: ixc_var_1424(klev_var_1346), ixlow_var_1425(kfdia_var_1345, klev_var_1346), ixhigh_var_1426(kfdia_var_1345, klev_var_1346)
  INTEGER(KIND = 4) :: ich_var_1427, icl_var_1428, ixc0_var_1429, ixp_var_1430, jc_var_1431, jl_var_1432
  laytrop_min_var_1422 = MINVAL(laytrop_var_1361)
  laytrop_max_var_1423 = MAXVAL(laytrop_var_1361)
  ixlow_var_1425 = 0
  ixhigh_var_1426 = 0
  ixc_var_1424 = 0
  DO lay_var_1385 = laytrop_min_var_1422 + 1, laytrop_max_var_1423
    icl_var_1428 = 0
    ich_var_1427 = 0
    DO jc_var_1431 = kidia_var_1344, kfdia_var_1345
      IF (lay_var_1385 <= laytrop_var_1361(jc_var_1431)) THEN
        icl_var_1428 = icl_var_1428 + 1
        ixlow_var_1425(icl_var_1428, lay_var_1385) = jc_var_1431
      ELSE
        ich_var_1427 = ich_var_1427 + 1
        ixhigh_var_1426(ich_var_1427, lay_var_1385) = jc_var_1431
      END IF
    END DO
    ixc_var_1424(lay_var_1385) = icl_var_1428
  END DO
  refrat_planck_a_var_1388 = chi_mls(1, 5) / chi_mls(2, 5)
  refrat_planck_b_var_1389 = chi_mls(3, 43) / chi_mls(2, 43)
  refrat_m_a_var_1390 = chi_mls(1, 7) / chi_mls(2, 7)
  DO lay_var_1385 = 1, laytrop_min_var_1422
    DO jl_var_1432 = kidia_var_1344, kfdia_var_1345
      speccomb_var_1375 = colh2o_var_1358(jl_var_1432, lay_var_1385) + rat_h2oco2_var_1366(jl_var_1432, lay_var_1385) * colco2_var_1359(jl_var_1432, lay_var_1385)
      specparm_var_1415 = MIN(colh2o_var_1358(jl_var_1432, lay_var_1385) / speccomb_var_1375, oneminus_var_1357)
      specmult_var_1414 = 8.0D0 * (specparm_var_1415)
      js_var_1384 = 1 + INT(specmult_var_1414)
      fs_var_1413 = ((specmult_var_1414) - AINT((specmult_var_1414)))
      speccomb1_var_1376 = colh2o_var_1358(jl_var_1432, lay_var_1385) + rat_h2oco2_1_var_1367(jl_var_1432, lay_var_1385) * colco2_var_1359(jl_var_1432, lay_var_1385)
      specparm1_var_1418 = MIN(colh2o_var_1358(jl_var_1432, lay_var_1385) / speccomb1_var_1376, oneminus_var_1357)
      specmult1_var_1417 = 8.0D0 * (specparm1_var_1418)
      js1_var_1386 = 1 + INT(specmult1_var_1417)
      fs1_var_1416 = ((specmult1_var_1417) - AINT((specmult1_var_1417)))
      speccomb_mo3 = colh2o_var_1358(jl_var_1432, lay_var_1385) + refrat_m_a_var_1390 * colco2_var_1359(jl_var_1432, lay_var_1385)
      specparm_mo3 = MIN(colh2o_var_1358(jl_var_1432, lay_var_1385) / speccomb_mo3, oneminus_var_1357)
      specmult_mo3 = 8.0D0 * specparm_mo3
      jmo3 = 1 + INT(specmult_mo3)
      fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))
      speccomb_planck_var_1377 = colh2o_var_1358(jl_var_1432, lay_var_1385) + refrat_planck_a_var_1388 * colco2_var_1359(jl_var_1432, lay_var_1385)
      specparm_planck_var_1421 = MIN(colh2o_var_1358(jl_var_1432, lay_var_1385) / speccomb_planck_var_1377, oneminus_var_1357)
      specmult_planck_var_1420 = 8.0D0 * specparm_planck_var_1421
      jpl_var_1387 = 1 + INT(specmult_planck_var_1420)
      fpl_var_1419 = ((specmult_planck_var_1420) - AINT((specmult_planck_var_1420)))
      ind0_var_1378 = ((jp_var_1354(jl_var_1432, lay_var_1385) - 1) * 5 + (jt_var_1355(jl_var_1432, lay_var_1385) - 1)) * nspa_var_216(5) + js_var_1384
      ind1_var_1379 = (jp_var_1354(jl_var_1432, lay_var_1385) * 5 + (jt1_var_1356(jl_var_1432, lay_var_1385) - 1)) * nspa_var_216(5) + js1_var_1386
      inds_var_1380 = indself_var_1364(jl_var_1432, lay_var_1385)
      indf_var_1381 = indfor_var_1370(jl_var_1432, lay_var_1385)
      indm_var_1382 = indminor_var_1374(jl_var_1432, lay_var_1385)
      IF (specparm_var_1415 .LT. 0.125D0) THEN
        p_var_1403 = fs_var_1413 - 1.0D0
        p4_var_1404 = p_var_1403 ** 4
        fk0_var_1405 = p4_var_1404
        fk1_var_1406 = 1.0D0 - p_var_1403 - 2.0D0 * p4_var_1404
        fk2_var_1407 = p_var_1403 + p4_var_1404
        fac000_var_1391 = fk0_var_1405 * fac00_var_1350(jl_var_1432, lay_var_1385)
        fac100_var_1392 = fk1_var_1406 * fac00_var_1350(jl_var_1432, lay_var_1385)
        fac200_var_1393 = fk2_var_1407 * fac00_var_1350(jl_var_1432, lay_var_1385)
        fac010_var_1394 = fk0_var_1405 * fac10_var_1352(jl_var_1432, lay_var_1385)
        fac110_var_1395 = fk1_var_1406 * fac10_var_1352(jl_var_1432, lay_var_1385)
        fac210_var_1396 = fk2_var_1407 * fac10_var_1352(jl_var_1432, lay_var_1385)
      ELSE IF (specparm_var_1415 .GT. 0.875D0) THEN
        p_var_1403 = - fs_var_1413
        p4_var_1404 = p_var_1403 ** 4
        fk0_var_1405 = p4_var_1404
        fk1_var_1406 = 1.0D0 - p_var_1403 - 2.0D0 * p4_var_1404
        fk2_var_1407 = p_var_1403 + p4_var_1404
        fac000_var_1391 = fk0_var_1405 * fac00_var_1350(jl_var_1432, lay_var_1385)
        fac100_var_1392 = fk1_var_1406 * fac00_var_1350(jl_var_1432, lay_var_1385)
        fac200_var_1393 = fk2_var_1407 * fac00_var_1350(jl_var_1432, lay_var_1385)
        fac010_var_1394 = fk0_var_1405 * fac10_var_1352(jl_var_1432, lay_var_1385)
        fac110_var_1395 = fk1_var_1406 * fac10_var_1352(jl_var_1432, lay_var_1385)
        fac210_var_1396 = fk2_var_1407 * fac10_var_1352(jl_var_1432, lay_var_1385)
      ELSE
        fac000_var_1391 = (1.0D0 - fs_var_1413) * fac00_var_1350(jl_var_1432, lay_var_1385)
        fac010_var_1394 = (1.0D0 - fs_var_1413) * fac10_var_1352(jl_var_1432, lay_var_1385)
        fac100_var_1392 = fs_var_1413 * fac00_var_1350(jl_var_1432, lay_var_1385)
        fac110_var_1395 = fs_var_1413 * fac10_var_1352(jl_var_1432, lay_var_1385)
        fac200_var_1393 = 0.0D0
        fac210_var_1396 = 0.0D0
      END IF
      IF (specparm1_var_1418 .LT. 0.125D0) THEN
        p_var_1403 = fs1_var_1416 - 1.0D0
        p4_var_1404 = p_var_1403 ** 4
        fk0_var_1405 = p4_var_1404
        fk1_var_1406 = 1.0D0 - p_var_1403 - 2.0D0 * p4_var_1404
        fk2_var_1407 = p_var_1403 + p4_var_1404
        fac001_var_1397 = fk0_var_1405 * fac01_var_1351(jl_var_1432, lay_var_1385)
        fac101_var_1398 = fk1_var_1406 * fac01_var_1351(jl_var_1432, lay_var_1385)
        fac201_var_1399 = fk2_var_1407 * fac01_var_1351(jl_var_1432, lay_var_1385)
        fac011_var_1400 = fk0_var_1405 * fac11_var_1353(jl_var_1432, lay_var_1385)
        fac111_var_1401 = fk1_var_1406 * fac11_var_1353(jl_var_1432, lay_var_1385)
        fac211_var_1402 = fk2_var_1407 * fac11_var_1353(jl_var_1432, lay_var_1385)
      ELSE IF (specparm1_var_1418 .GT. 0.875D0) THEN
        p_var_1403 = - fs1_var_1416
        p4_var_1404 = p_var_1403 ** 4
        fk0_var_1405 = p4_var_1404
        fk1_var_1406 = 1.0D0 - p_var_1403 - 2.0D0 * p4_var_1404
        fk2_var_1407 = p_var_1403 + p4_var_1404
        fac001_var_1397 = fk0_var_1405 * fac01_var_1351(jl_var_1432, lay_var_1385)
        fac101_var_1398 = fk1_var_1406 * fac01_var_1351(jl_var_1432, lay_var_1385)
        fac201_var_1399 = fk2_var_1407 * fac01_var_1351(jl_var_1432, lay_var_1385)
        fac011_var_1400 = fk0_var_1405 * fac11_var_1353(jl_var_1432, lay_var_1385)
        fac111_var_1401 = fk1_var_1406 * fac11_var_1353(jl_var_1432, lay_var_1385)
        fac211_var_1402 = fk2_var_1407 * fac11_var_1353(jl_var_1432, lay_var_1385)
      ELSE
        fac001_var_1397 = (1.0D0 - fs1_var_1416) * fac01_var_1351(jl_var_1432, lay_var_1385)
        fac011_var_1400 = (1.0D0 - fs1_var_1416) * fac11_var_1353(jl_var_1432, lay_var_1385)
        fac101_var_1398 = fs1_var_1416 * fac01_var_1351(jl_var_1432, lay_var_1385)
        fac111_var_1401 = fs1_var_1416 * fac11_var_1353(jl_var_1432, lay_var_1385)
        fac201_var_1399 = 0.0D0
        fac211_var_1402 = 0.0D0
      END IF
      IF (specparm_var_1415 .LT. 0.125D0) THEN
        tau_major_var_1410(1 : ng5) = speccomb_var_1375 * (fac000_var_1391 * absa_var_175(ind0_var_1378, 1 : 16) + fac100_var_1392 * absa_var_175(ind0_var_1378 + 1, 1 : 16) + fac200_var_1393 * absa_var_175(ind0_var_1378 + 2, 1 : 16) + fac010_var_1394 * absa_var_175(ind0_var_1378 + 9, 1 : 16) + fac110_var_1395 * absa_var_175(ind0_var_1378 + 10, 1 : 16) + fac210_var_1396 * absa_var_175(ind0_var_1378 + 11, 1 : 16))
      ELSE IF (specparm_var_1415 .GT. 0.875D0) THEN
        tau_major_var_1410(1 : ng5) = speccomb_var_1375 * (fac200_var_1393 * absa_var_175(ind0_var_1378 - 1, 1 : 16) + fac100_var_1392 * absa_var_175(ind0_var_1378, 1 : 16) + fac000_var_1391 * absa_var_175(ind0_var_1378 + 1, 1 : 16) + fac210_var_1396 * absa_var_175(ind0_var_1378 + 8, 1 : 16) + fac110_var_1395 * absa_var_175(ind0_var_1378 + 9, 1 : 16) + fac010_var_1394 * absa_var_175(ind0_var_1378 + 10, 1 : 16))
      ELSE
        tau_major_var_1410(1 : ng5) = speccomb_var_1375 * (fac000_var_1391 * absa_var_175(ind0_var_1378, 1 : 16) + fac100_var_1392 * absa_var_175(ind0_var_1378 + 1, 1 : 16) + fac010_var_1394 * absa_var_175(ind0_var_1378 + 9, 1 : 16) + fac110_var_1395 * absa_var_175(ind0_var_1378 + 10, 1 : 16))
      END IF
      IF (specparm1_var_1418 .LT. 0.125D0) THEN
        tau_major1_var_1411(1 : ng5) = speccomb1_var_1376 * (fac001_var_1397 * absa_var_175(ind1_var_1379, 1 : 16) + fac101_var_1398 * absa_var_175(ind1_var_1379 + 1, 1 : 16) + fac201_var_1399 * absa_var_175(ind1_var_1379 + 2, 1 : 16) + fac011_var_1400 * absa_var_175(ind1_var_1379 + 9, 1 : 16) + fac111_var_1401 * absa_var_175(ind1_var_1379 + 10, 1 : 16) + fac211_var_1402 * absa_var_175(ind1_var_1379 + 11, 1 : 16))
      ELSE IF (specparm1_var_1418 .GT. 0.875D0) THEN
        tau_major1_var_1411(1 : ng5) = speccomb1_var_1376 * (fac201_var_1399 * absa_var_175(ind1_var_1379 - 1, 1 : 16) + fac101_var_1398 * absa_var_175(ind1_var_1379, 1 : 16) + fac001_var_1397 * absa_var_175(ind1_var_1379 + 1, 1 : 16) + fac211_var_1402 * absa_var_175(ind1_var_1379 + 8, 1 : 16) + fac111_var_1401 * absa_var_175(ind1_var_1379 + 9, 1 : 16) + fac011_var_1400 * absa_var_175(ind1_var_1379 + 10, 1 : 16))
      ELSE
        tau_major1_var_1411(1 : ng5) = speccomb1_var_1376 * (fac001_var_1397 * absa_var_175(ind1_var_1379, 1 : 16) + fac101_var_1398 * absa_var_175(ind1_var_1379 + 1, 1 : 16) + fac011_var_1400 * absa_var_175(ind1_var_1379 + 9, 1 : 16) + fac111_var_1401 * absa_var_175(ind1_var_1379 + 10, 1 : 16))
      END IF
      DO ig_var_1383 = 1, 16
        tauself_var_1409 = selffac_var_1362(jl_var_1432, lay_var_1385) * (selfref_var_178(inds_var_1380, ig_var_1383) + selffrac_var_1363(jl_var_1432, lay_var_1385) * (selfref_var_178(inds_var_1380 + 1, ig_var_1383) - selfref_var_178(inds_var_1380, ig_var_1383)))
        taufor_var_1408 = forfac_var_1372(jl_var_1432, lay_var_1385) * (forref_var_179(indf_var_1381, ig_var_1383) + forfrac_var_1371(jl_var_1432, lay_var_1385) * (forref_var_179(indf_var_1381 + 1, ig_var_1383) - forref_var_179(indf_var_1381, ig_var_1383)))
        o3m1 = ka_mo3_var_177(jmo3, indm_var_1382, ig_var_1383) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_1382, ig_var_1383) - ka_mo3_var_177(jmo3, indm_var_1382, ig_var_1383))
        o3m2 = ka_mo3_var_177(jmo3, indm_var_1382 + 1, ig_var_1383) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_1382 + 1, ig_var_1383) - ka_mo3_var_177(jmo3, indm_var_1382 + 1, ig_var_1383))
        abso3_var_1412 = o3m1 + minorfrac_var_1373(jl_var_1432, lay_var_1385) * (o3m2 - o3m1)
        taug_var_1347(jl_var_1432, 52 + ig_var_1383, lay_var_1385) = tau_major_var_1410(ig_var_1383) + tau_major1_var_1411(ig_var_1383) + tauself_var_1409 + taufor_var_1408 + abso3_var_1412 * colo3_var_1360(jl_var_1432, lay_var_1385) + wx_var_1348(jl_var_1432, 1, lay_var_1385) * ccl4(ig_var_1383)
        fracs_var_1365(jl_var_1432, 52 + ig_var_1383, lay_var_1385) = fracrefa_var_173(ig_var_1383, jpl_var_1387) + fpl_var_1419 * (fracrefa_var_173(ig_var_1383, jpl_var_1387 + 1) - fracrefa_var_173(ig_var_1383, jpl_var_1387))
      END DO
    END DO
  END DO
  DO lay_var_1385 = laytrop_max_var_1423 + 1, klev_var_1346
    DO jl_var_1432 = kidia_var_1344, kfdia_var_1345
      speccomb_var_1375 = colo3_var_1360(jl_var_1432, lay_var_1385) + rat_o3co2_var_1368(jl_var_1432, lay_var_1385) * colco2_var_1359(jl_var_1432, lay_var_1385)
      specparm_var_1415 = MIN(colo3_var_1360(jl_var_1432, lay_var_1385) / speccomb_var_1375, oneminus_var_1357)
      specmult_var_1414 = 4.0D0 * (specparm_var_1415)
      js_var_1384 = 1 + INT(specmult_var_1414)
      fs_var_1413 = ((specmult_var_1414) - AINT((specmult_var_1414)))
      speccomb1_var_1376 = colo3_var_1360(jl_var_1432, lay_var_1385) + rat_o3co2_1_var_1369(jl_var_1432, lay_var_1385) * colco2_var_1359(jl_var_1432, lay_var_1385)
      specparm1_var_1418 = MIN(colo3_var_1360(jl_var_1432, lay_var_1385) / speccomb1_var_1376, oneminus_var_1357)
      specmult1_var_1417 = 4.0D0 * (specparm1_var_1418)
      js1_var_1386 = 1 + INT(specmult1_var_1417)
      fs1_var_1416 = ((specmult1_var_1417) - AINT((specmult1_var_1417)))
      fac000_var_1391 = (1.0D0 - fs_var_1413) * fac00_var_1350(jl_var_1432, lay_var_1385)
      fac010_var_1394 = (1.0D0 - fs_var_1413) * fac10_var_1352(jl_var_1432, lay_var_1385)
      fac100_var_1392 = fs_var_1413 * fac00_var_1350(jl_var_1432, lay_var_1385)
      fac110_var_1395 = fs_var_1413 * fac10_var_1352(jl_var_1432, lay_var_1385)
      fac001_var_1397 = (1.0D0 - fs1_var_1416) * fac01_var_1351(jl_var_1432, lay_var_1385)
      fac011_var_1400 = (1.0D0 - fs1_var_1416) * fac11_var_1353(jl_var_1432, lay_var_1385)
      fac101_var_1398 = fs1_var_1416 * fac01_var_1351(jl_var_1432, lay_var_1385)
      fac111_var_1401 = fs1_var_1416 * fac11_var_1353(jl_var_1432, lay_var_1385)
      speccomb_planck_var_1377 = colo3_var_1360(jl_var_1432, lay_var_1385) + refrat_planck_b_var_1389 * colco2_var_1359(jl_var_1432, lay_var_1385)
      specparm_planck_var_1421 = MIN(colo3_var_1360(jl_var_1432, lay_var_1385) / speccomb_planck_var_1377, oneminus_var_1357)
      specmult_planck_var_1420 = 4.0D0 * specparm_planck_var_1421
      jpl_var_1387 = 1 + INT(specmult_planck_var_1420)
      fpl_var_1419 = ((specmult_planck_var_1420) - AINT((specmult_planck_var_1420)))
      ind0_var_1378 = ((jp_var_1354(jl_var_1432, lay_var_1385) - 13) * 5 + (jt_var_1355(jl_var_1432, lay_var_1385) - 1)) * nspb_var_217(5) + js_var_1384
      ind1_var_1379 = ((jp_var_1354(jl_var_1432, lay_var_1385) - 12) * 5 + (jt1_var_1356(jl_var_1432, lay_var_1385) - 1)) * nspb_var_217(5) + js1_var_1386
      DO ig_var_1383 = 1, 16
        taug_var_1347(jl_var_1432, 52 + ig_var_1383, lay_var_1385) = speccomb_var_1375 * (fac000_var_1391 * absb_var_176(ind0_var_1378, ig_var_1383) + fac100_var_1392 * absb_var_176(ind0_var_1378 + 1, ig_var_1383) + fac010_var_1394 * absb_var_176(ind0_var_1378 + 5, ig_var_1383) + fac110_var_1395 * absb_var_176(ind0_var_1378 + 6, ig_var_1383)) + speccomb1_var_1376 * (fac001_var_1397 * absb_var_176(ind1_var_1379, ig_var_1383) + fac101_var_1398 * absb_var_176(ind1_var_1379 + 1, ig_var_1383) + fac011_var_1400 * absb_var_176(ind1_var_1379 + 5, ig_var_1383) + fac111_var_1401 * absb_var_176(ind1_var_1379 + 6, ig_var_1383)) + wx_var_1348(jl_var_1432, 1, lay_var_1385) * ccl4(ig_var_1383)
        fracs_var_1365(jl_var_1432, 52 + ig_var_1383, lay_var_1385) = fracrefb_var_174(ig_var_1383, jpl_var_1387) + fpl_var_1419 * (fracrefb_var_174(ig_var_1383, jpl_var_1387 + 1) - fracrefb_var_174(ig_var_1383, jpl_var_1387))
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1423 /= laytrop_min_var_1422) THEN
    DO lay_var_1385 = laytrop_min_var_1422 + 1, laytrop_max_var_1423
      ixc0_var_1429 = ixc_var_1424(lay_var_1385)
      DO ixp_var_1430 = 1, ixc0_var_1429
        jl_var_1432 = ixlow_var_1425(ixp_var_1430, lay_var_1385)
        speccomb_var_1375 = colh2o_var_1358(jl_var_1432, lay_var_1385) + rat_h2oco2_var_1366(jl_var_1432, lay_var_1385) * colco2_var_1359(jl_var_1432, lay_var_1385)
        specparm_var_1415 = MIN(colh2o_var_1358(jl_var_1432, lay_var_1385) / speccomb_var_1375, oneminus_var_1357)
        specmult_var_1414 = 8.0D0 * (specparm_var_1415)
        js_var_1384 = 1 + INT(specmult_var_1414)
        fs_var_1413 = ((specmult_var_1414) - AINT((specmult_var_1414)))
        speccomb1_var_1376 = colh2o_var_1358(jl_var_1432, lay_var_1385) + rat_h2oco2_1_var_1367(jl_var_1432, lay_var_1385) * colco2_var_1359(jl_var_1432, lay_var_1385)
        specparm1_var_1418 = MIN(colh2o_var_1358(jl_var_1432, lay_var_1385) / speccomb1_var_1376, oneminus_var_1357)
        specmult1_var_1417 = 8.0D0 * (specparm1_var_1418)
        js1_var_1386 = 1 + INT(specmult1_var_1417)
        fs1_var_1416 = ((specmult1_var_1417) - AINT((specmult1_var_1417)))
        speccomb_mo3 = colh2o_var_1358(jl_var_1432, lay_var_1385) + refrat_m_a_var_1390 * colco2_var_1359(jl_var_1432, lay_var_1385)
        specparm_mo3 = MIN(colh2o_var_1358(jl_var_1432, lay_var_1385) / speccomb_mo3, oneminus_var_1357)
        specmult_mo3 = 8.0D0 * specparm_mo3
        jmo3 = 1 + INT(specmult_mo3)
        fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))
        speccomb_planck_var_1377 = colh2o_var_1358(jl_var_1432, lay_var_1385) + refrat_planck_a_var_1388 * colco2_var_1359(jl_var_1432, lay_var_1385)
        specparm_planck_var_1421 = MIN(colh2o_var_1358(jl_var_1432, lay_var_1385) / speccomb_planck_var_1377, oneminus_var_1357)
        specmult_planck_var_1420 = 8.0D0 * specparm_planck_var_1421
        jpl_var_1387 = 1 + INT(specmult_planck_var_1420)
        fpl_var_1419 = ((specmult_planck_var_1420) - AINT((specmult_planck_var_1420)))
        ind0_var_1378 = ((jp_var_1354(jl_var_1432, lay_var_1385) - 1) * 5 + (jt_var_1355(jl_var_1432, lay_var_1385) - 1)) * nspa_var_216(5) + js_var_1384
        ind1_var_1379 = (jp_var_1354(jl_var_1432, lay_var_1385) * 5 + (jt1_var_1356(jl_var_1432, lay_var_1385) - 1)) * nspa_var_216(5) + js1_var_1386
        inds_var_1380 = indself_var_1364(jl_var_1432, lay_var_1385)
        indf_var_1381 = indfor_var_1370(jl_var_1432, lay_var_1385)
        indm_var_1382 = indminor_var_1374(jl_var_1432, lay_var_1385)
        IF (specparm_var_1415 .LT. 0.125D0) THEN
          p_var_1403 = fs_var_1413 - 1.0D0
          p4_var_1404 = p_var_1403 ** 4
          fk0_var_1405 = p4_var_1404
          fk1_var_1406 = 1.0D0 - p_var_1403 - 2.0D0 * p4_var_1404
          fk2_var_1407 = p_var_1403 + p4_var_1404
          fac000_var_1391 = fk0_var_1405 * fac00_var_1350(jl_var_1432, lay_var_1385)
          fac100_var_1392 = fk1_var_1406 * fac00_var_1350(jl_var_1432, lay_var_1385)
          fac200_var_1393 = fk2_var_1407 * fac00_var_1350(jl_var_1432, lay_var_1385)
          fac010_var_1394 = fk0_var_1405 * fac10_var_1352(jl_var_1432, lay_var_1385)
          fac110_var_1395 = fk1_var_1406 * fac10_var_1352(jl_var_1432, lay_var_1385)
          fac210_var_1396 = fk2_var_1407 * fac10_var_1352(jl_var_1432, lay_var_1385)
        ELSE IF (specparm_var_1415 .GT. 0.875D0) THEN
          p_var_1403 = - fs_var_1413
          p4_var_1404 = p_var_1403 ** 4
          fk0_var_1405 = p4_var_1404
          fk1_var_1406 = 1.0D0 - p_var_1403 - 2.0D0 * p4_var_1404
          fk2_var_1407 = p_var_1403 + p4_var_1404
          fac000_var_1391 = fk0_var_1405 * fac00_var_1350(jl_var_1432, lay_var_1385)
          fac100_var_1392 = fk1_var_1406 * fac00_var_1350(jl_var_1432, lay_var_1385)
          fac200_var_1393 = fk2_var_1407 * fac00_var_1350(jl_var_1432, lay_var_1385)
          fac010_var_1394 = fk0_var_1405 * fac10_var_1352(jl_var_1432, lay_var_1385)
          fac110_var_1395 = fk1_var_1406 * fac10_var_1352(jl_var_1432, lay_var_1385)
          fac210_var_1396 = fk2_var_1407 * fac10_var_1352(jl_var_1432, lay_var_1385)
        ELSE
          fac000_var_1391 = (1.0D0 - fs_var_1413) * fac00_var_1350(jl_var_1432, lay_var_1385)
          fac010_var_1394 = (1.0D0 - fs_var_1413) * fac10_var_1352(jl_var_1432, lay_var_1385)
          fac100_var_1392 = fs_var_1413 * fac00_var_1350(jl_var_1432, lay_var_1385)
          fac110_var_1395 = fs_var_1413 * fac10_var_1352(jl_var_1432, lay_var_1385)
          fac200_var_1393 = 0.0D0
          fac210_var_1396 = 0.0D0
        END IF
        IF (specparm1_var_1418 .LT. 0.125D0) THEN
          p_var_1403 = fs1_var_1416 - 1.0D0
          p4_var_1404 = p_var_1403 ** 4
          fk0_var_1405 = p4_var_1404
          fk1_var_1406 = 1.0D0 - p_var_1403 - 2.0D0 * p4_var_1404
          fk2_var_1407 = p_var_1403 + p4_var_1404
          fac001_var_1397 = fk0_var_1405 * fac01_var_1351(jl_var_1432, lay_var_1385)
          fac101_var_1398 = fk1_var_1406 * fac01_var_1351(jl_var_1432, lay_var_1385)
          fac201_var_1399 = fk2_var_1407 * fac01_var_1351(jl_var_1432, lay_var_1385)
          fac011_var_1400 = fk0_var_1405 * fac11_var_1353(jl_var_1432, lay_var_1385)
          fac111_var_1401 = fk1_var_1406 * fac11_var_1353(jl_var_1432, lay_var_1385)
          fac211_var_1402 = fk2_var_1407 * fac11_var_1353(jl_var_1432, lay_var_1385)
        ELSE IF (specparm1_var_1418 .GT. 0.875D0) THEN
          p_var_1403 = - fs1_var_1416
          p4_var_1404 = p_var_1403 ** 4
          fk0_var_1405 = p4_var_1404
          fk1_var_1406 = 1.0D0 - p_var_1403 - 2.0D0 * p4_var_1404
          fk2_var_1407 = p_var_1403 + p4_var_1404
          fac001_var_1397 = fk0_var_1405 * fac01_var_1351(jl_var_1432, lay_var_1385)
          fac101_var_1398 = fk1_var_1406 * fac01_var_1351(jl_var_1432, lay_var_1385)
          fac201_var_1399 = fk2_var_1407 * fac01_var_1351(jl_var_1432, lay_var_1385)
          fac011_var_1400 = fk0_var_1405 * fac11_var_1353(jl_var_1432, lay_var_1385)
          fac111_var_1401 = fk1_var_1406 * fac11_var_1353(jl_var_1432, lay_var_1385)
          fac211_var_1402 = fk2_var_1407 * fac11_var_1353(jl_var_1432, lay_var_1385)
        ELSE
          fac001_var_1397 = (1.0D0 - fs1_var_1416) * fac01_var_1351(jl_var_1432, lay_var_1385)
          fac011_var_1400 = (1.0D0 - fs1_var_1416) * fac11_var_1353(jl_var_1432, lay_var_1385)
          fac101_var_1398 = fs1_var_1416 * fac01_var_1351(jl_var_1432, lay_var_1385)
          fac111_var_1401 = fs1_var_1416 * fac11_var_1353(jl_var_1432, lay_var_1385)
          fac201_var_1399 = 0.0D0
          fac211_var_1402 = 0.0D0
        END IF
        IF (specparm_var_1415 .LT. 0.125D0) THEN
          tau_major_var_1410(1 : ng5) = speccomb_var_1375 * (fac000_var_1391 * absa_var_175(ind0_var_1378, 1 : 16) + fac100_var_1392 * absa_var_175(ind0_var_1378 + 1, 1 : 16) + fac200_var_1393 * absa_var_175(ind0_var_1378 + 2, 1 : 16) + fac010_var_1394 * absa_var_175(ind0_var_1378 + 9, 1 : 16) + fac110_var_1395 * absa_var_175(ind0_var_1378 + 10, 1 : 16) + fac210_var_1396 * absa_var_175(ind0_var_1378 + 11, 1 : 16))
        ELSE IF (specparm_var_1415 .GT. 0.875D0) THEN
          tau_major_var_1410(1 : ng5) = speccomb_var_1375 * (fac200_var_1393 * absa_var_175(ind0_var_1378 - 1, 1 : 16) + fac100_var_1392 * absa_var_175(ind0_var_1378, 1 : 16) + fac000_var_1391 * absa_var_175(ind0_var_1378 + 1, 1 : 16) + fac210_var_1396 * absa_var_175(ind0_var_1378 + 8, 1 : 16) + fac110_var_1395 * absa_var_175(ind0_var_1378 + 9, 1 : 16) + fac010_var_1394 * absa_var_175(ind0_var_1378 + 10, 1 : 16))
        ELSE
          tau_major_var_1410(1 : ng5) = speccomb_var_1375 * (fac000_var_1391 * absa_var_175(ind0_var_1378, 1 : 16) + fac100_var_1392 * absa_var_175(ind0_var_1378 + 1, 1 : 16) + fac010_var_1394 * absa_var_175(ind0_var_1378 + 9, 1 : 16) + fac110_var_1395 * absa_var_175(ind0_var_1378 + 10, 1 : 16))
        END IF
        IF (specparm1_var_1418 .LT. 0.125D0) THEN
          tau_major1_var_1411(1 : ng5) = speccomb1_var_1376 * (fac001_var_1397 * absa_var_175(ind1_var_1379, 1 : 16) + fac101_var_1398 * absa_var_175(ind1_var_1379 + 1, 1 : 16) + fac201_var_1399 * absa_var_175(ind1_var_1379 + 2, 1 : 16) + fac011_var_1400 * absa_var_175(ind1_var_1379 + 9, 1 : 16) + fac111_var_1401 * absa_var_175(ind1_var_1379 + 10, 1 : 16) + fac211_var_1402 * absa_var_175(ind1_var_1379 + 11, 1 : 16))
        ELSE IF (specparm1_var_1418 .GT. 0.875D0) THEN
          tau_major1_var_1411(1 : ng5) = speccomb1_var_1376 * (fac201_var_1399 * absa_var_175(ind1_var_1379 - 1, 1 : 16) + fac101_var_1398 * absa_var_175(ind1_var_1379, 1 : 16) + fac001_var_1397 * absa_var_175(ind1_var_1379 + 1, 1 : 16) + fac211_var_1402 * absa_var_175(ind1_var_1379 + 8, 1 : 16) + fac111_var_1401 * absa_var_175(ind1_var_1379 + 9, 1 : 16) + fac011_var_1400 * absa_var_175(ind1_var_1379 + 10, 1 : 16))
        ELSE
          tau_major1_var_1411(1 : ng5) = speccomb1_var_1376 * (fac001_var_1397 * absa_var_175(ind1_var_1379, 1 : 16) + fac101_var_1398 * absa_var_175(ind1_var_1379 + 1, 1 : 16) + fac011_var_1400 * absa_var_175(ind1_var_1379 + 9, 1 : 16) + fac111_var_1401 * absa_var_175(ind1_var_1379 + 10, 1 : 16))
        END IF
        DO ig_var_1383 = 1, 16
          tauself_var_1409 = selffac_var_1362(jl_var_1432, lay_var_1385) * (selfref_var_178(inds_var_1380, ig_var_1383) + selffrac_var_1363(jl_var_1432, lay_var_1385) * (selfref_var_178(inds_var_1380 + 1, ig_var_1383) - selfref_var_178(inds_var_1380, ig_var_1383)))
          taufor_var_1408 = forfac_var_1372(jl_var_1432, lay_var_1385) * (forref_var_179(indf_var_1381, ig_var_1383) + forfrac_var_1371(jl_var_1432, lay_var_1385) * (forref_var_179(indf_var_1381 + 1, ig_var_1383) - forref_var_179(indf_var_1381, ig_var_1383)))
          o3m1 = ka_mo3_var_177(jmo3, indm_var_1382, ig_var_1383) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_1382, ig_var_1383) - ka_mo3_var_177(jmo3, indm_var_1382, ig_var_1383))
          o3m2 = ka_mo3_var_177(jmo3, indm_var_1382 + 1, ig_var_1383) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_1382 + 1, ig_var_1383) - ka_mo3_var_177(jmo3, indm_var_1382 + 1, ig_var_1383))
          abso3_var_1412 = o3m1 + minorfrac_var_1373(jl_var_1432, lay_var_1385) * (o3m2 - o3m1)
          taug_var_1347(jl_var_1432, 52 + ig_var_1383, lay_var_1385) = tau_major_var_1410(ig_var_1383) + tau_major1_var_1411(ig_var_1383) + tauself_var_1409 + taufor_var_1408 + abso3_var_1412 * colo3_var_1360(jl_var_1432, lay_var_1385) + wx_var_1348(jl_var_1432, 1, lay_var_1385) * ccl4(ig_var_1383)
          fracs_var_1365(jl_var_1432, 52 + ig_var_1383, lay_var_1385) = fracrefa_var_173(ig_var_1383, jpl_var_1387) + fpl_var_1419 * (fracrefa_var_173(ig_var_1383, jpl_var_1387 + 1) - fracrefa_var_173(ig_var_1383, jpl_var_1387))
        END DO
      END DO
      ixc0_var_1429 = kfdia_var_1345 - kidia_var_1344 + 1 - ixc0_var_1429
      DO ixp_var_1430 = 1, ixc0_var_1429
        jl_var_1432 = ixhigh_var_1426(ixp_var_1430, lay_var_1385)
        speccomb_var_1375 = colo3_var_1360(jl_var_1432, lay_var_1385) + rat_o3co2_var_1368(jl_var_1432, lay_var_1385) * colco2_var_1359(jl_var_1432, lay_var_1385)
        specparm_var_1415 = MIN(colo3_var_1360(jl_var_1432, lay_var_1385) / speccomb_var_1375, oneminus_var_1357)
        specmult_var_1414 = 4.0D0 * (specparm_var_1415)
        js_var_1384 = 1 + INT(specmult_var_1414)
        fs_var_1413 = ((specmult_var_1414) - AINT((specmult_var_1414)))
        speccomb1_var_1376 = colo3_var_1360(jl_var_1432, lay_var_1385) + rat_o3co2_1_var_1369(jl_var_1432, lay_var_1385) * colco2_var_1359(jl_var_1432, lay_var_1385)
        specparm1_var_1418 = MIN(colo3_var_1360(jl_var_1432, lay_var_1385) / speccomb1_var_1376, oneminus_var_1357)
        specmult1_var_1417 = 4.0D0 * (specparm1_var_1418)
        js1_var_1386 = 1 + INT(specmult1_var_1417)
        fs1_var_1416 = ((specmult1_var_1417) - AINT((specmult1_var_1417)))
        fac000_var_1391 = (1.0D0 - fs_var_1413) * fac00_var_1350(jl_var_1432, lay_var_1385)
        fac010_var_1394 = (1.0D0 - fs_var_1413) * fac10_var_1352(jl_var_1432, lay_var_1385)
        fac100_var_1392 = fs_var_1413 * fac00_var_1350(jl_var_1432, lay_var_1385)
        fac110_var_1395 = fs_var_1413 * fac10_var_1352(jl_var_1432, lay_var_1385)
        fac001_var_1397 = (1.0D0 - fs1_var_1416) * fac01_var_1351(jl_var_1432, lay_var_1385)
        fac011_var_1400 = (1.0D0 - fs1_var_1416) * fac11_var_1353(jl_var_1432, lay_var_1385)
        fac101_var_1398 = fs1_var_1416 * fac01_var_1351(jl_var_1432, lay_var_1385)
        fac111_var_1401 = fs1_var_1416 * fac11_var_1353(jl_var_1432, lay_var_1385)
        speccomb_planck_var_1377 = colo3_var_1360(jl_var_1432, lay_var_1385) + refrat_planck_b_var_1389 * colco2_var_1359(jl_var_1432, lay_var_1385)
        specparm_planck_var_1421 = MIN(colo3_var_1360(jl_var_1432, lay_var_1385) / speccomb_planck_var_1377, oneminus_var_1357)
        specmult_planck_var_1420 = 4.0D0 * specparm_planck_var_1421
        jpl_var_1387 = 1 + INT(specmult_planck_var_1420)
        fpl_var_1419 = ((specmult_planck_var_1420) - AINT((specmult_planck_var_1420)))
        ind0_var_1378 = ((jp_var_1354(jl_var_1432, lay_var_1385) - 13) * 5 + (jt_var_1355(jl_var_1432, lay_var_1385) - 1)) * nspb_var_217(5) + js_var_1384
        ind1_var_1379 = ((jp_var_1354(jl_var_1432, lay_var_1385) - 12) * 5 + (jt1_var_1356(jl_var_1432, lay_var_1385) - 1)) * nspb_var_217(5) + js1_var_1386
        DO ig_var_1383 = 1, 16
          taug_var_1347(jl_var_1432, 52 + ig_var_1383, lay_var_1385) = speccomb_var_1375 * (fac000_var_1391 * absb_var_176(ind0_var_1378, ig_var_1383) + fac100_var_1392 * absb_var_176(ind0_var_1378 + 1, ig_var_1383) + fac010_var_1394 * absb_var_176(ind0_var_1378 + 5, ig_var_1383) + fac110_var_1395 * absb_var_176(ind0_var_1378 + 6, ig_var_1383)) + speccomb1_var_1376 * (fac001_var_1397 * absb_var_176(ind1_var_1379, ig_var_1383) + fac101_var_1398 * absb_var_176(ind1_var_1379 + 1, ig_var_1383) + fac011_var_1400 * absb_var_176(ind1_var_1379 + 5, ig_var_1383) + fac111_var_1401 * absb_var_176(ind1_var_1379 + 6, ig_var_1383)) + wx_var_1348(jl_var_1432, 1, lay_var_1385) * ccl4(ig_var_1383)
          fracs_var_1365(jl_var_1432, 52 + ig_var_1383, lay_var_1385) = fracrefb_var_174(ig_var_1383, jpl_var_1387) + fpl_var_1419 * (fracrefb_var_174(ig_var_1383, jpl_var_1387 + 1) - fracrefb_var_174(ig_var_1383, jpl_var_1387))
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol5
SUBROUTINE srtm_taumol19(kidia_var_1433, kfdia_var_1434, klev_var_1435, p_fac00_var_1436, p_fac01_var_1437, p_fac10_var_1438, p_fac11_var_1439, k_jp_var_1440, k_jt_var_1441, k_jt1_var_1442, p_oneminus_var_1443, p_colh2o_var_1444, p_colco2_var_1445, p_colmol_var_1446, k_laytrop_var_1447, p_selffac_var_1448, p_selffrac_var_1449, k_indself_var_1450, p_forfac_var_1451, p_forfrac_var_1452, k_indfor_var_1453, p_sfluxzen_var_1454, p_taug_var_1455, p_taur_var_1456, prmu0_var_1457)
  USE yoesrta19, ONLY: absa_var_244, absb_var_245, forrefc_var_247, layreffr_var_243, rayl_var_241, selfrefc_var_246, sfluxrefc_var_248, strrat_var_242
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1433, kfdia_var_1434
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1435
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1436(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1437(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1438(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1439(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1440(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1441(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1442(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1443(kidia_var_1433 : kfdia_var_1434)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1444(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1445(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1446(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1447(kidia_var_1433 : kfdia_var_1434)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1448(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1449(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1450(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1451(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1452(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1453(kidia_var_1433 : kfdia_var_1434, klev_var_1435)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1454(kidia_var_1433 : kfdia_var_1434, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1455(kidia_var_1433 : kfdia_var_1434, klev_var_1435, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1456(kidia_var_1433 : kfdia_var_1434, klev_var_1435, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1457(kidia_var_1433 : kfdia_var_1434)
  INTEGER(KIND = 4) :: ig_var_1458, ind0_var_1459, ind1_var_1460, inds_var_1461, indf_var_1462, js_var_1463, i_lay_var_1464, i_laysolfr_var_1465(kidia_var_1433 : kfdia_var_1434), i_nlayers_var_1466, iplon_var_1467
  INTEGER(KIND = 4) :: laytrop_min_var_1468, laytrop_max_var_1469
  REAL(KIND = 8) :: z_fs_var_1470, z_speccomb_var_1471, z_specmult_var_1472, z_specparm_var_1473, z_tauray_var_1474
  laytrop_min_var_1468 = MINVAL(k_laytrop_var_1447(kidia_var_1433 : kfdia_var_1434))
  laytrop_max_var_1469 = MAXVAL(k_laytrop_var_1447(kidia_var_1433 : kfdia_var_1434))
  i_nlayers_var_1466 = klev_var_1435
  DO iplon_var_1467 = kidia_var_1433, kfdia_var_1434
    i_laysolfr_var_1465(iplon_var_1467) = k_laytrop_var_1447(iplon_var_1467)
  END DO
  DO i_lay_var_1464 = 1, laytrop_min_var_1468
    DO iplon_var_1467 = kidia_var_1433, kfdia_var_1434
      IF (k_jp_var_1440(iplon_var_1467, i_lay_var_1464) < layreffr_var_243 .AND. k_jp_var_1440(iplon_var_1467, i_lay_var_1464 + 1) >= layreffr_var_243) i_laysolfr_var_1465(iplon_var_1467) = MIN(i_lay_var_1464 + 1, k_laytrop_var_1447(iplon_var_1467))
      z_speccomb_var_1471 = p_colh2o_var_1444(iplon_var_1467, i_lay_var_1464) + strrat_var_242 * p_colco2_var_1445(iplon_var_1467, i_lay_var_1464)
      z_specparm_var_1473 = p_colh2o_var_1444(iplon_var_1467, i_lay_var_1464) / z_speccomb_var_1471
      z_specparm_var_1473 = MIN(p_oneminus_var_1443(iplon_var_1467), z_specparm_var_1473)
      z_specmult_var_1472 = 8.0D0 * (z_specparm_var_1473)
      js_var_1463 = 1 + INT(z_specmult_var_1472)
      z_fs_var_1470 = z_specmult_var_1472 - AINT(z_specmult_var_1472)
      ind0_var_1459 = ((k_jp_var_1440(iplon_var_1467, i_lay_var_1464) - 1) * 5 + (k_jt_var_1441(iplon_var_1467, i_lay_var_1464) - 1)) * nspa_var_313(19) + js_var_1463
      ind1_var_1460 = (k_jp_var_1440(iplon_var_1467, i_lay_var_1464) * 5 + (k_jt1_var_1442(iplon_var_1467, i_lay_var_1464) - 1)) * nspa_var_313(19) + js_var_1463
      inds_var_1461 = k_indself_var_1450(iplon_var_1467, i_lay_var_1464)
      indf_var_1462 = k_indfor_var_1453(iplon_var_1467, i_lay_var_1464)
      z_tauray_var_1474 = p_colmol_var_1446(iplon_var_1467, i_lay_var_1464) * rayl_var_241
      DO ig_var_1458 = 1, 8
        p_taug_var_1455(iplon_var_1467, i_lay_var_1464, ig_var_1458) = z_speccomb_var_1471 * ((1.0D0 - z_fs_var_1470) * (absa_var_244(ind0_var_1459, ig_var_1458) * p_fac00_var_1436(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind0_var_1459 + 9, ig_var_1458) * p_fac10_var_1438(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind1_var_1460, ig_var_1458) * p_fac01_var_1437(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind1_var_1460 + 9, ig_var_1458) * p_fac11_var_1439(iplon_var_1467, i_lay_var_1464)) + z_fs_var_1470 * (absa_var_244(ind0_var_1459 + 1, ig_var_1458) * p_fac00_var_1436(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind0_var_1459 + 10, ig_var_1458) * p_fac10_var_1438(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind1_var_1460 + 1, ig_var_1458) * p_fac01_var_1437(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind1_var_1460 + 10, ig_var_1458) * p_fac11_var_1439(iplon_var_1467, i_lay_var_1464))) + p_colh2o_var_1444(iplon_var_1467, i_lay_var_1464) * (p_selffac_var_1448(iplon_var_1467, i_lay_var_1464) * (selfrefc_var_246(inds_var_1461, ig_var_1458) + p_selffrac_var_1449(iplon_var_1467, i_lay_var_1464) * (selfrefc_var_246(inds_var_1461 + 1, ig_var_1458) - selfrefc_var_246(inds_var_1461, ig_var_1458))) + p_forfac_var_1451(iplon_var_1467, i_lay_var_1464) * (forrefc_var_247(indf_var_1462, ig_var_1458) + p_forfrac_var_1452(iplon_var_1467, i_lay_var_1464) * (forrefc_var_247(indf_var_1462 + 1, ig_var_1458) - forrefc_var_247(indf_var_1462, ig_var_1458))))
        IF (i_lay_var_1464 == i_laysolfr_var_1465(iplon_var_1467)) p_sfluxzen_var_1454(iplon_var_1467, ig_var_1458) = sfluxrefc_var_248(ig_var_1458, js_var_1463) + z_fs_var_1470 * (sfluxrefc_var_248(ig_var_1458, js_var_1463 + 1) - sfluxrefc_var_248(ig_var_1458, js_var_1463))
        p_taur_var_1456(iplon_var_1467, i_lay_var_1464, ig_var_1458) = z_tauray_var_1474
      END DO
    END DO
  END DO
  DO i_lay_var_1464 = laytrop_min_var_1468 + 1, laytrop_max_var_1469
    DO iplon_var_1467 = kidia_var_1433, kfdia_var_1434
      IF (i_lay_var_1464 <= k_laytrop_var_1447(iplon_var_1467)) THEN
        IF (k_jp_var_1440(iplon_var_1467, i_lay_var_1464) < layreffr_var_243 .AND. k_jp_var_1440(iplon_var_1467, i_lay_var_1464 + 1) >= layreffr_var_243) i_laysolfr_var_1465(iplon_var_1467) = MIN(i_lay_var_1464 + 1, k_laytrop_var_1447(iplon_var_1467))
        z_speccomb_var_1471 = p_colh2o_var_1444(iplon_var_1467, i_lay_var_1464) + strrat_var_242 * p_colco2_var_1445(iplon_var_1467, i_lay_var_1464)
        z_specparm_var_1473 = p_colh2o_var_1444(iplon_var_1467, i_lay_var_1464) / z_speccomb_var_1471
        z_specparm_var_1473 = MIN(p_oneminus_var_1443(iplon_var_1467), z_specparm_var_1473)
        z_specmult_var_1472 = 8.0D0 * (z_specparm_var_1473)
        js_var_1463 = 1 + INT(z_specmult_var_1472)
        z_fs_var_1470 = z_specmult_var_1472 - AINT(z_specmult_var_1472)
        ind0_var_1459 = ((k_jp_var_1440(iplon_var_1467, i_lay_var_1464) - 1) * 5 + (k_jt_var_1441(iplon_var_1467, i_lay_var_1464) - 1)) * nspa_var_313(19) + js_var_1463
        ind1_var_1460 = (k_jp_var_1440(iplon_var_1467, i_lay_var_1464) * 5 + (k_jt1_var_1442(iplon_var_1467, i_lay_var_1464) - 1)) * nspa_var_313(19) + js_var_1463
        inds_var_1461 = k_indself_var_1450(iplon_var_1467, i_lay_var_1464)
        indf_var_1462 = k_indfor_var_1453(iplon_var_1467, i_lay_var_1464)
        z_tauray_var_1474 = p_colmol_var_1446(iplon_var_1467, i_lay_var_1464) * rayl_var_241
        DO ig_var_1458 = 1, 8
          p_taug_var_1455(iplon_var_1467, i_lay_var_1464, ig_var_1458) = z_speccomb_var_1471 * ((1.0D0 - z_fs_var_1470) * (absa_var_244(ind0_var_1459, ig_var_1458) * p_fac00_var_1436(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind0_var_1459 + 9, ig_var_1458) * p_fac10_var_1438(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind1_var_1460, ig_var_1458) * p_fac01_var_1437(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind1_var_1460 + 9, ig_var_1458) * p_fac11_var_1439(iplon_var_1467, i_lay_var_1464)) + z_fs_var_1470 * (absa_var_244(ind0_var_1459 + 1, ig_var_1458) * p_fac00_var_1436(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind0_var_1459 + 10, ig_var_1458) * p_fac10_var_1438(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind1_var_1460 + 1, ig_var_1458) * p_fac01_var_1437(iplon_var_1467, i_lay_var_1464) + absa_var_244(ind1_var_1460 + 10, ig_var_1458) * p_fac11_var_1439(iplon_var_1467, i_lay_var_1464))) + p_colh2o_var_1444(iplon_var_1467, i_lay_var_1464) * (p_selffac_var_1448(iplon_var_1467, i_lay_var_1464) * (selfrefc_var_246(inds_var_1461, ig_var_1458) + p_selffrac_var_1449(iplon_var_1467, i_lay_var_1464) * (selfrefc_var_246(inds_var_1461 + 1, ig_var_1458) - selfrefc_var_246(inds_var_1461, ig_var_1458))) + p_forfac_var_1451(iplon_var_1467, i_lay_var_1464) * (forrefc_var_247(indf_var_1462, ig_var_1458) + p_forfrac_var_1452(iplon_var_1467, i_lay_var_1464) * (forrefc_var_247(indf_var_1462 + 1, ig_var_1458) - forrefc_var_247(indf_var_1462, ig_var_1458))))
          IF (i_lay_var_1464 == i_laysolfr_var_1465(iplon_var_1467)) p_sfluxzen_var_1454(iplon_var_1467, ig_var_1458) = sfluxrefc_var_248(ig_var_1458, js_var_1463) + z_fs_var_1470 * (sfluxrefc_var_248(ig_var_1458, js_var_1463 + 1) - sfluxrefc_var_248(ig_var_1458, js_var_1463))
          p_taur_var_1456(iplon_var_1467, i_lay_var_1464, ig_var_1458) = z_tauray_var_1474
        END DO
      ELSE
        ind0_var_1459 = ((k_jp_var_1440(iplon_var_1467, i_lay_var_1464) - 13) * 5 + (k_jt_var_1441(iplon_var_1467, i_lay_var_1464) - 1)) * nspb_var_314(19) + 1
        ind1_var_1460 = ((k_jp_var_1440(iplon_var_1467, i_lay_var_1464) - 12) * 5 + (k_jt1_var_1442(iplon_var_1467, i_lay_var_1464) - 1)) * nspb_var_314(19) + 1
        z_tauray_var_1474 = p_colmol_var_1446(iplon_var_1467, i_lay_var_1464) * rayl_var_241
        DO ig_var_1458 = 1, 8
          p_taug_var_1455(iplon_var_1467, i_lay_var_1464, ig_var_1458) = p_colco2_var_1445(iplon_var_1467, i_lay_var_1464) * (p_fac00_var_1436(iplon_var_1467, i_lay_var_1464) * absb_var_245(ind0_var_1459, ig_var_1458) + p_fac10_var_1438(iplon_var_1467, i_lay_var_1464) * absb_var_245(ind0_var_1459 + 1, ig_var_1458) + p_fac01_var_1437(iplon_var_1467, i_lay_var_1464) * absb_var_245(ind1_var_1460, ig_var_1458) + p_fac11_var_1439(iplon_var_1467, i_lay_var_1464) * absb_var_245(ind1_var_1460 + 1, ig_var_1458))
          p_taur_var_1456(iplon_var_1467, i_lay_var_1464, ig_var_1458) = z_tauray_var_1474
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1464 = laytrop_max_var_1469 + 1, i_nlayers_var_1466
    DO iplon_var_1467 = kidia_var_1433, kfdia_var_1434
      ind0_var_1459 = ((k_jp_var_1440(iplon_var_1467, i_lay_var_1464) - 13) * 5 + (k_jt_var_1441(iplon_var_1467, i_lay_var_1464) - 1)) * nspb_var_314(19) + 1
      ind1_var_1460 = ((k_jp_var_1440(iplon_var_1467, i_lay_var_1464) - 12) * 5 + (k_jt1_var_1442(iplon_var_1467, i_lay_var_1464) - 1)) * nspb_var_314(19) + 1
      z_tauray_var_1474 = p_colmol_var_1446(iplon_var_1467, i_lay_var_1464) * rayl_var_241
      DO ig_var_1458 = 1, 8
        p_taug_var_1455(iplon_var_1467, i_lay_var_1464, ig_var_1458) = p_colco2_var_1445(iplon_var_1467, i_lay_var_1464) * (p_fac00_var_1436(iplon_var_1467, i_lay_var_1464) * absb_var_245(ind0_var_1459, ig_var_1458) + p_fac10_var_1438(iplon_var_1467, i_lay_var_1464) * absb_var_245(ind0_var_1459 + 1, ig_var_1458) + p_fac01_var_1437(iplon_var_1467, i_lay_var_1464) * absb_var_245(ind1_var_1460, ig_var_1458) + p_fac11_var_1439(iplon_var_1467, i_lay_var_1464) * absb_var_245(ind1_var_1460 + 1, ig_var_1458))
        p_taur_var_1456(iplon_var_1467, i_lay_var_1464, ig_var_1458) = z_tauray_var_1474
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol19
SUBROUTINE rrtm_taumol15(kidia_var_1475, kfdia_var_1476, klev_var_1477, taug_var_1478, p_tauaerl_var_1479, fac00_var_1480, fac01_var_1481, fac10_var_1482, fac11_var_1483, forfac_var_1497, forfrac_var_1498, indfor_var_1496, jp_var_1484, jt_var_1485, jt1_var_1486, oneminus_var_1487, colh2o_var_1488, colco2_var_1489, coln2o_var_1490, laytrop_var_1491, selffac_var_1492, selffrac_var_1493, indself_var_1494, fracs_var_1495, rat_n2oco2, rat_n2oco2_1, minorfrac_var_1499, indminor_var_1500, scaleminor_var_1501, colbrd_var_1502)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrtm, ONLY: ng15
  USE yoerrta15, ONLY: absa_var_143, forref_var_146, fracrefa_var_142, ka_mn2_var_144, selfref_var_145
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1475
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1476
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1477
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1478(kidia_var_1475 : kfdia_var_1476, 140, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1479(kidia_var_1475 : kfdia_var_1476, klev_var_1477, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1480(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1481(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1482(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1483(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1484(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1485(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1486(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1487
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1488(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1489(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_1490(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1491(kidia_var_1475 : kfdia_var_1476)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1492(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1493(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1494(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1495(kidia_var_1475 : kfdia_var_1476, 140, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: rat_n2oco2(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: rat_n2oco2_1(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1496(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1497(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1498(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1499(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1500(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: scaleminor_var_1501(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  REAL(KIND = 8), INTENT(IN) :: colbrd_var_1502(kidia_var_1475 : kfdia_var_1476, klev_var_1477)
  INTEGER(KIND = 4) :: ig_var_1503, ind0_var_1504, ind1_var_1505, inds_var_1506, indf_var_1507, indm_var_1508, js_var_1509, js1_var_1510, jpl_var_1511, jmn2, lay_var_1512
  REAL(KIND = 8) :: refrat_planck_a_var_1513, refrat_m_a_var_1514
  REAL(KIND = 8) :: taufor_var_1515, tauself_var_1516, tau_major_var_1517(2), tau_major1_var_1518(2), n2m1, n2m2, taun2_var_1519, scalen2_var_1520
  REAL(KIND = 8) :: fac000_var_1521, fac100_var_1522, fac200_var_1523, fac010_var_1524, fac110_var_1525, fac210_var_1526, fac001_var_1527, fac101_var_1528, fac201_var_1529, fac011_var_1530, fac111_var_1531, fac211_var_1532
  REAL(KIND = 8) :: p_var_1533, p4_var_1534, fk0_var_1535, fk1_var_1536, fk2_var_1537
  REAL(KIND = 8) :: fs_var_1538, specmult_var_1539, specparm_var_1540, speccomb_var_1541, fs1_var_1542, specmult1_var_1543, specparm1_var_1544, speccomb1_var_1545, fmn2, specmult_mn2, specparm_mn2, speccomb_mn2, fpl_var_1546, specmult_planck_var_1547, specparm_planck_var_1548, speccomb_planck_var_1549
  INTEGER(KIND = 4) :: laytrop_min_var_1550, laytrop_max_var_1551
  INTEGER(KIND = 4) :: ixc_var_1552(klev_var_1477), ixlow_var_1553(kfdia_var_1476, klev_var_1477), ixhigh_var_1554(kfdia_var_1476, klev_var_1477)
  INTEGER(KIND = 4) :: ich_var_1555, icl_var_1556, ixc0_var_1557, ixp_var_1558, jc_var_1559, jl_var_1560
  laytrop_min_var_1550 = MINVAL(laytrop_var_1491)
  laytrop_max_var_1551 = MAXVAL(laytrop_var_1491)
  ixlow_var_1553 = 0
  ixhigh_var_1554 = 0
  ixc_var_1552 = 0
  DO lay_var_1512 = laytrop_min_var_1550 + 1, laytrop_max_var_1551
    icl_var_1556 = 0
    ich_var_1555 = 0
    DO jc_var_1559 = kidia_var_1475, kfdia_var_1476
      IF (lay_var_1512 <= laytrop_var_1491(jc_var_1559)) THEN
        icl_var_1556 = icl_var_1556 + 1
        ixlow_var_1553(icl_var_1556, lay_var_1512) = jc_var_1559
      ELSE
        ich_var_1555 = ich_var_1555 + 1
        ixhigh_var_1554(ich_var_1555, lay_var_1512) = jc_var_1559
      END IF
    END DO
    ixc_var_1552(lay_var_1512) = icl_var_1556
  END DO
  refrat_planck_a_var_1513 = chi_mls(4, 1) / chi_mls(2, 1)
  refrat_m_a_var_1514 = chi_mls(4, 1) / chi_mls(2, 1)
  DO lay_var_1512 = 1, laytrop_min_var_1550
    DO jl_var_1560 = kidia_var_1475, kfdia_var_1476
      speccomb_var_1541 = coln2o_var_1490(jl_var_1560, lay_var_1512) + rat_n2oco2(jl_var_1560, lay_var_1512) * colco2_var_1489(jl_var_1560, lay_var_1512)
      specparm_var_1540 = MIN(coln2o_var_1490(jl_var_1560, lay_var_1512) / speccomb_var_1541, oneminus_var_1487)
      specmult_var_1539 = 8.0D0 * (specparm_var_1540)
      js_var_1509 = 1 + INT(specmult_var_1539)
      fs_var_1538 = ((specmult_var_1539) - AINT((specmult_var_1539)))
      speccomb1_var_1545 = coln2o_var_1490(jl_var_1560, lay_var_1512) + rat_n2oco2_1(jl_var_1560, lay_var_1512) * colco2_var_1489(jl_var_1560, lay_var_1512)
      specparm1_var_1544 = MIN(coln2o_var_1490(jl_var_1560, lay_var_1512) / speccomb1_var_1545, oneminus_var_1487)
      specmult1_var_1543 = 8.0D0 * (specparm1_var_1544)
      js1_var_1510 = 1 + INT(specmult1_var_1543)
      fs1_var_1542 = ((specmult1_var_1543) - AINT((specmult1_var_1543)))
      speccomb_mn2 = coln2o_var_1490(jl_var_1560, lay_var_1512) + refrat_m_a_var_1514 * colco2_var_1489(jl_var_1560, lay_var_1512)
      specparm_mn2 = MIN(coln2o_var_1490(jl_var_1560, lay_var_1512) / speccomb_mn2, oneminus_var_1487)
      specmult_mn2 = 8.0D0 * specparm_mn2
      jmn2 = 1 + INT(specmult_mn2)
      fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))
      speccomb_planck_var_1549 = coln2o_var_1490(jl_var_1560, lay_var_1512) + refrat_planck_a_var_1513 * colco2_var_1489(jl_var_1560, lay_var_1512)
      specparm_planck_var_1548 = MIN(coln2o_var_1490(jl_var_1560, lay_var_1512) / speccomb_planck_var_1549, oneminus_var_1487)
      specmult_planck_var_1547 = 8.0D0 * specparm_planck_var_1548
      jpl_var_1511 = 1 + INT(specmult_planck_var_1547)
      fpl_var_1546 = ((specmult_planck_var_1547) - AINT((specmult_planck_var_1547)))
      ind0_var_1504 = ((jp_var_1484(jl_var_1560, lay_var_1512) - 1) * 5 + (jt_var_1485(jl_var_1560, lay_var_1512) - 1)) * nspa_var_216(15) + js_var_1509
      ind1_var_1505 = (jp_var_1484(jl_var_1560, lay_var_1512) * 5 + (jt1_var_1486(jl_var_1560, lay_var_1512) - 1)) * nspa_var_216(15) + js1_var_1510
      inds_var_1506 = indself_var_1494(jl_var_1560, lay_var_1512)
      indf_var_1507 = indfor_var_1496(jl_var_1560, lay_var_1512)
      indm_var_1508 = indminor_var_1500(jl_var_1560, lay_var_1512)
      scalen2_var_1520 = colbrd_var_1502(jl_var_1560, lay_var_1512) * scaleminor_var_1501(jl_var_1560, lay_var_1512)
      IF (specparm_var_1540 .LT. 0.125D0) THEN
        p_var_1533 = fs_var_1538 - 1.0D0
        p4_var_1534 = p_var_1533 ** 4
        fk0_var_1535 = p4_var_1534
        fk1_var_1536 = 1.0D0 - p_var_1533 - 2.0D0 * p4_var_1534
        fk2_var_1537 = p_var_1533 + p4_var_1534
        fac000_var_1521 = fk0_var_1535 * fac00_var_1480(jl_var_1560, lay_var_1512)
        fac100_var_1522 = fk1_var_1536 * fac00_var_1480(jl_var_1560, lay_var_1512)
        fac200_var_1523 = fk2_var_1537 * fac00_var_1480(jl_var_1560, lay_var_1512)
        fac010_var_1524 = fk0_var_1535 * fac10_var_1482(jl_var_1560, lay_var_1512)
        fac110_var_1525 = fk1_var_1536 * fac10_var_1482(jl_var_1560, lay_var_1512)
        fac210_var_1526 = fk2_var_1537 * fac10_var_1482(jl_var_1560, lay_var_1512)
      ELSE IF (specparm_var_1540 .GT. 0.875D0) THEN
        p_var_1533 = - fs_var_1538
        p4_var_1534 = p_var_1533 ** 4
        fk0_var_1535 = p4_var_1534
        fk1_var_1536 = 1.0D0 - p_var_1533 - 2.0D0 * p4_var_1534
        fk2_var_1537 = p_var_1533 + p4_var_1534
        fac000_var_1521 = fk0_var_1535 * fac00_var_1480(jl_var_1560, lay_var_1512)
        fac100_var_1522 = fk1_var_1536 * fac00_var_1480(jl_var_1560, lay_var_1512)
        fac200_var_1523 = fk2_var_1537 * fac00_var_1480(jl_var_1560, lay_var_1512)
        fac010_var_1524 = fk0_var_1535 * fac10_var_1482(jl_var_1560, lay_var_1512)
        fac110_var_1525 = fk1_var_1536 * fac10_var_1482(jl_var_1560, lay_var_1512)
        fac210_var_1526 = fk2_var_1537 * fac10_var_1482(jl_var_1560, lay_var_1512)
      ELSE
        fac000_var_1521 = (1.0D0 - fs_var_1538) * fac00_var_1480(jl_var_1560, lay_var_1512)
        fac010_var_1524 = (1.0D0 - fs_var_1538) * fac10_var_1482(jl_var_1560, lay_var_1512)
        fac100_var_1522 = fs_var_1538 * fac00_var_1480(jl_var_1560, lay_var_1512)
        fac110_var_1525 = fs_var_1538 * fac10_var_1482(jl_var_1560, lay_var_1512)
        fac200_var_1523 = 0.0D0
        fac210_var_1526 = 0.0D0
      END IF
      IF (specparm1_var_1544 .LT. 0.125D0) THEN
        p_var_1533 = fs1_var_1542 - 1.0D0
        p4_var_1534 = p_var_1533 ** 4
        fk0_var_1535 = p4_var_1534
        fk1_var_1536 = 1.0D0 - p_var_1533 - 2.0D0 * p4_var_1534
        fk2_var_1537 = p_var_1533 + p4_var_1534
        fac001_var_1527 = fk0_var_1535 * fac01_var_1481(jl_var_1560, lay_var_1512)
        fac101_var_1528 = fk1_var_1536 * fac01_var_1481(jl_var_1560, lay_var_1512)
        fac201_var_1529 = fk2_var_1537 * fac01_var_1481(jl_var_1560, lay_var_1512)
        fac011_var_1530 = fk0_var_1535 * fac11_var_1483(jl_var_1560, lay_var_1512)
        fac111_var_1531 = fk1_var_1536 * fac11_var_1483(jl_var_1560, lay_var_1512)
        fac211_var_1532 = fk2_var_1537 * fac11_var_1483(jl_var_1560, lay_var_1512)
      ELSE IF (specparm1_var_1544 .GT. 0.875D0) THEN
        p_var_1533 = - fs1_var_1542
        p4_var_1534 = p_var_1533 ** 4
        fk0_var_1535 = p4_var_1534
        fk1_var_1536 = 1.0D0 - p_var_1533 - 2.0D0 * p4_var_1534
        fk2_var_1537 = p_var_1533 + p4_var_1534
        fac001_var_1527 = fk0_var_1535 * fac01_var_1481(jl_var_1560, lay_var_1512)
        fac101_var_1528 = fk1_var_1536 * fac01_var_1481(jl_var_1560, lay_var_1512)
        fac201_var_1529 = fk2_var_1537 * fac01_var_1481(jl_var_1560, lay_var_1512)
        fac011_var_1530 = fk0_var_1535 * fac11_var_1483(jl_var_1560, lay_var_1512)
        fac111_var_1531 = fk1_var_1536 * fac11_var_1483(jl_var_1560, lay_var_1512)
        fac211_var_1532 = fk2_var_1537 * fac11_var_1483(jl_var_1560, lay_var_1512)
      ELSE
        fac001_var_1527 = (1.0D0 - fs1_var_1542) * fac01_var_1481(jl_var_1560, lay_var_1512)
        fac011_var_1530 = (1.0D0 - fs1_var_1542) * fac11_var_1483(jl_var_1560, lay_var_1512)
        fac101_var_1528 = fs1_var_1542 * fac01_var_1481(jl_var_1560, lay_var_1512)
        fac111_var_1531 = fs1_var_1542 * fac11_var_1483(jl_var_1560, lay_var_1512)
        fac201_var_1529 = 0.0D0
        fac211_var_1532 = 0.0D0
      END IF
      IF (specparm_var_1540 .LT. 0.125D0) THEN
        tau_major_var_1517(1 : ng15) = speccomb_var_1541 * (fac000_var_1521 * absa_var_143(ind0_var_1504, 1 : 2) + fac100_var_1522 * absa_var_143(ind0_var_1504 + 1, 1 : 2) + fac200_var_1523 * absa_var_143(ind0_var_1504 + 2, 1 : 2) + fac010_var_1524 * absa_var_143(ind0_var_1504 + 9, 1 : 2) + fac110_var_1525 * absa_var_143(ind0_var_1504 + 10, 1 : 2) + fac210_var_1526 * absa_var_143(ind0_var_1504 + 11, 1 : 2))
      ELSE IF (specparm_var_1540 .GT. 0.875D0) THEN
        tau_major_var_1517(1 : ng15) = speccomb_var_1541 * (fac200_var_1523 * absa_var_143(ind0_var_1504 - 1, 1 : 2) + fac100_var_1522 * absa_var_143(ind0_var_1504, 1 : 2) + fac000_var_1521 * absa_var_143(ind0_var_1504 + 1, 1 : 2) + fac210_var_1526 * absa_var_143(ind0_var_1504 + 8, 1 : 2) + fac110_var_1525 * absa_var_143(ind0_var_1504 + 9, 1 : 2) + fac010_var_1524 * absa_var_143(ind0_var_1504 + 10, 1 : 2))
      ELSE
        tau_major_var_1517(1 : ng15) = speccomb_var_1541 * (fac000_var_1521 * absa_var_143(ind0_var_1504, 1 : 2) + fac100_var_1522 * absa_var_143(ind0_var_1504 + 1, 1 : 2) + fac010_var_1524 * absa_var_143(ind0_var_1504 + 9, 1 : 2) + fac110_var_1525 * absa_var_143(ind0_var_1504 + 10, 1 : 2))
      END IF
      IF (specparm1_var_1544 .LT. 0.125D0) THEN
        tau_major1_var_1518(1 : ng15) = speccomb1_var_1545 * (fac001_var_1527 * absa_var_143(ind1_var_1505, 1 : 2) + fac101_var_1528 * absa_var_143(ind1_var_1505 + 1, 1 : 2) + fac201_var_1529 * absa_var_143(ind1_var_1505 + 2, 1 : 2) + fac011_var_1530 * absa_var_143(ind1_var_1505 + 9, 1 : 2) + fac111_var_1531 * absa_var_143(ind1_var_1505 + 10, 1 : 2) + fac211_var_1532 * absa_var_143(ind1_var_1505 + 11, 1 : 2))
      ELSE IF (specparm1_var_1544 .GT. 0.875D0) THEN
        tau_major1_var_1518(1 : ng15) = speccomb1_var_1545 * (fac201_var_1529 * absa_var_143(ind1_var_1505 - 1, 1 : 2) + fac101_var_1528 * absa_var_143(ind1_var_1505, 1 : 2) + fac001_var_1527 * absa_var_143(ind1_var_1505 + 1, 1 : 2) + fac211_var_1532 * absa_var_143(ind1_var_1505 + 8, 1 : 2) + fac111_var_1531 * absa_var_143(ind1_var_1505 + 9, 1 : 2) + fac011_var_1530 * absa_var_143(ind1_var_1505 + 10, 1 : 2))
      ELSE
        tau_major1_var_1518(1 : ng15) = speccomb1_var_1545 * (fac001_var_1527 * absa_var_143(ind1_var_1505, 1 : 2) + fac101_var_1528 * absa_var_143(ind1_var_1505 + 1, 1 : 2) + fac011_var_1530 * absa_var_143(ind1_var_1505 + 9, 1 : 2) + fac111_var_1531 * absa_var_143(ind1_var_1505 + 10, 1 : 2))
      END IF
      DO ig_var_1503 = 1, 2
        tauself_var_1516 = selffac_var_1492(jl_var_1560, lay_var_1512) * (selfref_var_145(inds_var_1506, ig_var_1503) + selffrac_var_1493(jl_var_1560, lay_var_1512) * (selfref_var_145(inds_var_1506 + 1, ig_var_1503) - selfref_var_145(inds_var_1506, ig_var_1503)))
        taufor_var_1515 = forfac_var_1497(jl_var_1560, lay_var_1512) * (forref_var_146(indf_var_1507, ig_var_1503) + forfrac_var_1498(jl_var_1560, lay_var_1512) * (forref_var_146(indf_var_1507 + 1, ig_var_1503) - forref_var_146(indf_var_1507, ig_var_1503)))
        n2m1 = ka_mn2_var_144(jmn2, indm_var_1508, ig_var_1503) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_1508, ig_var_1503) - ka_mn2_var_144(jmn2, indm_var_1508, ig_var_1503))
        n2m2 = ka_mn2_var_144(jmn2, indm_var_1508 + 1, ig_var_1503) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_1508 + 1, ig_var_1503) - ka_mn2_var_144(jmn2, indm_var_1508 + 1, ig_var_1503))
        taun2_var_1519 = scalen2_var_1520 * (n2m1 + minorfrac_var_1499(jl_var_1560, lay_var_1512) * (n2m2 - n2m1))
        taug_var_1478(jl_var_1560, 136 + ig_var_1503, lay_var_1512) = tau_major_var_1517(ig_var_1503) + tau_major1_var_1518(ig_var_1503) + tauself_var_1516 + taufor_var_1515 + taun2_var_1519
        fracs_var_1495(jl_var_1560, 136 + ig_var_1503, lay_var_1512) = fracrefa_var_142(ig_var_1503, jpl_var_1511) + fpl_var_1546 * (fracrefa_var_142(ig_var_1503, jpl_var_1511 + 1) - fracrefa_var_142(ig_var_1503, jpl_var_1511))
      END DO
    END DO
  END DO
  DO ig_var_1503 = 1, 2
    DO lay_var_1512 = laytrop_max_var_1551 + 1, klev_var_1477
      DO jl_var_1560 = kidia_var_1475, kfdia_var_1476
        taug_var_1478(jl_var_1560, 136 + ig_var_1503, lay_var_1512) = 0.0D0
        fracs_var_1495(jl_var_1560, 136 + ig_var_1503, lay_var_1512) = 0.0D0
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1551 /= laytrop_min_var_1550) THEN
    DO lay_var_1512 = laytrop_min_var_1550 + 1, laytrop_max_var_1551
      ixc0_var_1557 = ixc_var_1552(lay_var_1512)
      DO ixp_var_1558 = 1, ixc0_var_1557
        jl_var_1560 = ixlow_var_1553(ixp_var_1558, lay_var_1512)
        speccomb_var_1541 = coln2o_var_1490(jl_var_1560, lay_var_1512) + rat_n2oco2(jl_var_1560, lay_var_1512) * colco2_var_1489(jl_var_1560, lay_var_1512)
        specparm_var_1540 = MIN(coln2o_var_1490(jl_var_1560, lay_var_1512) / speccomb_var_1541, oneminus_var_1487)
        specmult_var_1539 = 8.0D0 * (specparm_var_1540)
        js_var_1509 = 1 + INT(specmult_var_1539)
        fs_var_1538 = ((specmult_var_1539) - AINT((specmult_var_1539)))
        speccomb1_var_1545 = coln2o_var_1490(jl_var_1560, lay_var_1512) + rat_n2oco2_1(jl_var_1560, lay_var_1512) * colco2_var_1489(jl_var_1560, lay_var_1512)
        specparm1_var_1544 = MIN(coln2o_var_1490(jl_var_1560, lay_var_1512) / speccomb1_var_1545, oneminus_var_1487)
        specmult1_var_1543 = 8.0D0 * (specparm1_var_1544)
        js1_var_1510 = 1 + INT(specmult1_var_1543)
        fs1_var_1542 = ((specmult1_var_1543) - AINT((specmult1_var_1543)))
        speccomb_mn2 = coln2o_var_1490(jl_var_1560, lay_var_1512) + refrat_m_a_var_1514 * colco2_var_1489(jl_var_1560, lay_var_1512)
        specparm_mn2 = MIN(coln2o_var_1490(jl_var_1560, lay_var_1512) / speccomb_mn2, oneminus_var_1487)
        specmult_mn2 = 8.0D0 * specparm_mn2
        jmn2 = 1 + INT(specmult_mn2)
        fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))
        speccomb_planck_var_1549 = coln2o_var_1490(jl_var_1560, lay_var_1512) + refrat_planck_a_var_1513 * colco2_var_1489(jl_var_1560, lay_var_1512)
        specparm_planck_var_1548 = MIN(coln2o_var_1490(jl_var_1560, lay_var_1512) / speccomb_planck_var_1549, oneminus_var_1487)
        specmult_planck_var_1547 = 8.0D0 * specparm_planck_var_1548
        jpl_var_1511 = 1 + INT(specmult_planck_var_1547)
        fpl_var_1546 = ((specmult_planck_var_1547) - AINT((specmult_planck_var_1547)))
        ind0_var_1504 = ((jp_var_1484(jl_var_1560, lay_var_1512) - 1) * 5 + (jt_var_1485(jl_var_1560, lay_var_1512) - 1)) * nspa_var_216(15) + js_var_1509
        ind1_var_1505 = (jp_var_1484(jl_var_1560, lay_var_1512) * 5 + (jt1_var_1486(jl_var_1560, lay_var_1512) - 1)) * nspa_var_216(15) + js1_var_1510
        inds_var_1506 = indself_var_1494(jl_var_1560, lay_var_1512)
        indf_var_1507 = indfor_var_1496(jl_var_1560, lay_var_1512)
        indm_var_1508 = indminor_var_1500(jl_var_1560, lay_var_1512)
        scalen2_var_1520 = colbrd_var_1502(jl_var_1560, lay_var_1512) * scaleminor_var_1501(jl_var_1560, lay_var_1512)
        IF (specparm_var_1540 .LT. 0.125D0) THEN
          p_var_1533 = fs_var_1538 - 1.0D0
          p4_var_1534 = p_var_1533 ** 4
          fk0_var_1535 = p4_var_1534
          fk1_var_1536 = 1.0D0 - p_var_1533 - 2.0D0 * p4_var_1534
          fk2_var_1537 = p_var_1533 + p4_var_1534
          fac000_var_1521 = fk0_var_1535 * fac00_var_1480(jl_var_1560, lay_var_1512)
          fac100_var_1522 = fk1_var_1536 * fac00_var_1480(jl_var_1560, lay_var_1512)
          fac200_var_1523 = fk2_var_1537 * fac00_var_1480(jl_var_1560, lay_var_1512)
          fac010_var_1524 = fk0_var_1535 * fac10_var_1482(jl_var_1560, lay_var_1512)
          fac110_var_1525 = fk1_var_1536 * fac10_var_1482(jl_var_1560, lay_var_1512)
          fac210_var_1526 = fk2_var_1537 * fac10_var_1482(jl_var_1560, lay_var_1512)
        ELSE IF (specparm_var_1540 .GT. 0.875D0) THEN
          p_var_1533 = - fs_var_1538
          p4_var_1534 = p_var_1533 ** 4
          fk0_var_1535 = p4_var_1534
          fk1_var_1536 = 1.0D0 - p_var_1533 - 2.0D0 * p4_var_1534
          fk2_var_1537 = p_var_1533 + p4_var_1534
          fac000_var_1521 = fk0_var_1535 * fac00_var_1480(jl_var_1560, lay_var_1512)
          fac100_var_1522 = fk1_var_1536 * fac00_var_1480(jl_var_1560, lay_var_1512)
          fac200_var_1523 = fk2_var_1537 * fac00_var_1480(jl_var_1560, lay_var_1512)
          fac010_var_1524 = fk0_var_1535 * fac10_var_1482(jl_var_1560, lay_var_1512)
          fac110_var_1525 = fk1_var_1536 * fac10_var_1482(jl_var_1560, lay_var_1512)
          fac210_var_1526 = fk2_var_1537 * fac10_var_1482(jl_var_1560, lay_var_1512)
        ELSE
          fac000_var_1521 = (1.0D0 - fs_var_1538) * fac00_var_1480(jl_var_1560, lay_var_1512)
          fac010_var_1524 = (1.0D0 - fs_var_1538) * fac10_var_1482(jl_var_1560, lay_var_1512)
          fac100_var_1522 = fs_var_1538 * fac00_var_1480(jl_var_1560, lay_var_1512)
          fac110_var_1525 = fs_var_1538 * fac10_var_1482(jl_var_1560, lay_var_1512)
          fac200_var_1523 = 0.0D0
          fac210_var_1526 = 0.0D0
        END IF
        IF (specparm1_var_1544 .LT. 0.125D0) THEN
          p_var_1533 = fs1_var_1542 - 1.0D0
          p4_var_1534 = p_var_1533 ** 4
          fk0_var_1535 = p4_var_1534
          fk1_var_1536 = 1.0D0 - p_var_1533 - 2.0D0 * p4_var_1534
          fk2_var_1537 = p_var_1533 + p4_var_1534
          fac001_var_1527 = fk0_var_1535 * fac01_var_1481(jl_var_1560, lay_var_1512)
          fac101_var_1528 = fk1_var_1536 * fac01_var_1481(jl_var_1560, lay_var_1512)
          fac201_var_1529 = fk2_var_1537 * fac01_var_1481(jl_var_1560, lay_var_1512)
          fac011_var_1530 = fk0_var_1535 * fac11_var_1483(jl_var_1560, lay_var_1512)
          fac111_var_1531 = fk1_var_1536 * fac11_var_1483(jl_var_1560, lay_var_1512)
          fac211_var_1532 = fk2_var_1537 * fac11_var_1483(jl_var_1560, lay_var_1512)
        ELSE IF (specparm1_var_1544 .GT. 0.875D0) THEN
          p_var_1533 = - fs1_var_1542
          p4_var_1534 = p_var_1533 ** 4
          fk0_var_1535 = p4_var_1534
          fk1_var_1536 = 1.0D0 - p_var_1533 - 2.0D0 * p4_var_1534
          fk2_var_1537 = p_var_1533 + p4_var_1534
          fac001_var_1527 = fk0_var_1535 * fac01_var_1481(jl_var_1560, lay_var_1512)
          fac101_var_1528 = fk1_var_1536 * fac01_var_1481(jl_var_1560, lay_var_1512)
          fac201_var_1529 = fk2_var_1537 * fac01_var_1481(jl_var_1560, lay_var_1512)
          fac011_var_1530 = fk0_var_1535 * fac11_var_1483(jl_var_1560, lay_var_1512)
          fac111_var_1531 = fk1_var_1536 * fac11_var_1483(jl_var_1560, lay_var_1512)
          fac211_var_1532 = fk2_var_1537 * fac11_var_1483(jl_var_1560, lay_var_1512)
        ELSE
          fac001_var_1527 = (1.0D0 - fs1_var_1542) * fac01_var_1481(jl_var_1560, lay_var_1512)
          fac011_var_1530 = (1.0D0 - fs1_var_1542) * fac11_var_1483(jl_var_1560, lay_var_1512)
          fac101_var_1528 = fs1_var_1542 * fac01_var_1481(jl_var_1560, lay_var_1512)
          fac111_var_1531 = fs1_var_1542 * fac11_var_1483(jl_var_1560, lay_var_1512)
          fac201_var_1529 = 0.0D0
          fac211_var_1532 = 0.0D0
        END IF
        IF (specparm_var_1540 .LT. 0.125D0) THEN
          tau_major_var_1517(1 : ng15) = speccomb_var_1541 * (fac000_var_1521 * absa_var_143(ind0_var_1504, 1 : 2) + fac100_var_1522 * absa_var_143(ind0_var_1504 + 1, 1 : 2) + fac200_var_1523 * absa_var_143(ind0_var_1504 + 2, 1 : 2) + fac010_var_1524 * absa_var_143(ind0_var_1504 + 9, 1 : 2) + fac110_var_1525 * absa_var_143(ind0_var_1504 + 10, 1 : 2) + fac210_var_1526 * absa_var_143(ind0_var_1504 + 11, 1 : 2))
        ELSE IF (specparm_var_1540 .GT. 0.875D0) THEN
          tau_major_var_1517(1 : ng15) = speccomb_var_1541 * (fac200_var_1523 * absa_var_143(ind0_var_1504 - 1, 1 : 2) + fac100_var_1522 * absa_var_143(ind0_var_1504, 1 : 2) + fac000_var_1521 * absa_var_143(ind0_var_1504 + 1, 1 : 2) + fac210_var_1526 * absa_var_143(ind0_var_1504 + 8, 1 : 2) + fac110_var_1525 * absa_var_143(ind0_var_1504 + 9, 1 : 2) + fac010_var_1524 * absa_var_143(ind0_var_1504 + 10, 1 : 2))
        ELSE
          tau_major_var_1517(1 : ng15) = speccomb_var_1541 * (fac000_var_1521 * absa_var_143(ind0_var_1504, 1 : 2) + fac100_var_1522 * absa_var_143(ind0_var_1504 + 1, 1 : 2) + fac010_var_1524 * absa_var_143(ind0_var_1504 + 9, 1 : 2) + fac110_var_1525 * absa_var_143(ind0_var_1504 + 10, 1 : 2))
        END IF
        IF (specparm1_var_1544 .LT. 0.125D0) THEN
          tau_major1_var_1518(1 : ng15) = speccomb1_var_1545 * (fac001_var_1527 * absa_var_143(ind1_var_1505, 1 : 2) + fac101_var_1528 * absa_var_143(ind1_var_1505 + 1, 1 : 2) + fac201_var_1529 * absa_var_143(ind1_var_1505 + 2, 1 : 2) + fac011_var_1530 * absa_var_143(ind1_var_1505 + 9, 1 : 2) + fac111_var_1531 * absa_var_143(ind1_var_1505 + 10, 1 : 2) + fac211_var_1532 * absa_var_143(ind1_var_1505 + 11, 1 : 2))
        ELSE IF (specparm1_var_1544 .GT. 0.875D0) THEN
          tau_major1_var_1518(1 : ng15) = speccomb1_var_1545 * (fac201_var_1529 * absa_var_143(ind1_var_1505 - 1, 1 : 2) + fac101_var_1528 * absa_var_143(ind1_var_1505, 1 : 2) + fac001_var_1527 * absa_var_143(ind1_var_1505 + 1, 1 : 2) + fac211_var_1532 * absa_var_143(ind1_var_1505 + 8, 1 : 2) + fac111_var_1531 * absa_var_143(ind1_var_1505 + 9, 1 : 2) + fac011_var_1530 * absa_var_143(ind1_var_1505 + 10, 1 : 2))
        ELSE
          tau_major1_var_1518(1 : ng15) = speccomb1_var_1545 * (fac001_var_1527 * absa_var_143(ind1_var_1505, 1 : 2) + fac101_var_1528 * absa_var_143(ind1_var_1505 + 1, 1 : 2) + fac011_var_1530 * absa_var_143(ind1_var_1505 + 9, 1 : 2) + fac111_var_1531 * absa_var_143(ind1_var_1505 + 10, 1 : 2))
        END IF
        DO ig_var_1503 = 1, 2
          tauself_var_1516 = selffac_var_1492(jl_var_1560, lay_var_1512) * (selfref_var_145(inds_var_1506, ig_var_1503) + selffrac_var_1493(jl_var_1560, lay_var_1512) * (selfref_var_145(inds_var_1506 + 1, ig_var_1503) - selfref_var_145(inds_var_1506, ig_var_1503)))
          taufor_var_1515 = forfac_var_1497(jl_var_1560, lay_var_1512) * (forref_var_146(indf_var_1507, ig_var_1503) + forfrac_var_1498(jl_var_1560, lay_var_1512) * (forref_var_146(indf_var_1507 + 1, ig_var_1503) - forref_var_146(indf_var_1507, ig_var_1503)))
          n2m1 = ka_mn2_var_144(jmn2, indm_var_1508, ig_var_1503) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_1508, ig_var_1503) - ka_mn2_var_144(jmn2, indm_var_1508, ig_var_1503))
          n2m2 = ka_mn2_var_144(jmn2, indm_var_1508 + 1, ig_var_1503) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_1508 + 1, ig_var_1503) - ka_mn2_var_144(jmn2, indm_var_1508 + 1, ig_var_1503))
          taun2_var_1519 = scalen2_var_1520 * (n2m1 + minorfrac_var_1499(jl_var_1560, lay_var_1512) * (n2m2 - n2m1))
          taug_var_1478(jl_var_1560, 136 + ig_var_1503, lay_var_1512) = tau_major_var_1517(ig_var_1503) + tau_major1_var_1518(ig_var_1503) + tauself_var_1516 + taufor_var_1515 + taun2_var_1519
          fracs_var_1495(jl_var_1560, 136 + ig_var_1503, lay_var_1512) = fracrefa_var_142(ig_var_1503, jpl_var_1511) + fpl_var_1546 * (fracrefa_var_142(ig_var_1503, jpl_var_1511 + 1) - fracrefa_var_142(ig_var_1503, jpl_var_1511))
        END DO
      END DO
      ixc0_var_1557 = kfdia_var_1476 - kidia_var_1475 + 1 - ixc0_var_1557
      DO ig_var_1503 = 1, 2
        DO ixp_var_1558 = 1, ixc0_var_1557
          jl_var_1560 = ixhigh_var_1554(ixp_var_1558, lay_var_1512)
          taug_var_1478(jl_var_1560, 136 + ig_var_1503, lay_var_1512) = 0.0D0
          fracs_var_1495(jl_var_1560, 136 + ig_var_1503, lay_var_1512) = 0.0D0
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol15
SUBROUTINE rrtm_taumol13(kidia_var_1561, kfdia_var_1562, klev_var_1563, taug_var_1564, p_tauaerl_var_1565, fac00_var_1566, fac01_var_1567, fac10_var_1568, fac11_var_1569, forfac_var_1585, forfrac_var_1586, indfor_var_1584, jp_var_1570, jt_var_1571, jt1_var_1572, oneminus_var_1573, colh2o_var_1574, coln2o_var_1575, colco2_var_1576, colo3_var_1577, coldry_var_1578, laytrop_var_1579, selffac_var_1580, selffrac_var_1581, indself_var_1582, fracs_var_1583, rat_h2on2o, rat_h2on2o_1, minorfrac_var_1587, indminor_var_1588)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrtm, ONLY: ng13
  USE yoerrta13, ONLY: absa_var_132, forref_var_134, fracrefa_var_130, fracrefb_var_131, ka_mco, ka_mco2_var_135, kb_mo3, selfref_var_133
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1561
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1562
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1563
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1564(kidia_var_1561 : kfdia_var_1562, 140, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1565(kidia_var_1561 : kfdia_var_1562, klev_var_1563, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1566(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1567(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1568(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1569(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1570(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1571(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1572(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1573
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1574(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_1575(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1576(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1577(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1578(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1579(kidia_var_1561 : kfdia_var_1562)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1580(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1581(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1582(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1583(kidia_var_1561 : kfdia_var_1562, 140, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: rat_h2on2o(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: rat_h2on2o_1(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1584(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1585(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1586(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1587(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1588(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  REAL(KIND = 8) :: speccomb_var_1589, speccomb1_var_1590, speccomb_planck_var_1591, speccomb_mco2_var_1592, speccomb_mco
  INTEGER(KIND = 4) :: ind0_var_1593, ind1_var_1594, inds_var_1595, indf_var_1596, indm_var_1597
  INTEGER(KIND = 4) :: ig_var_1598, js_var_1599, lay_var_1600, js1_var_1601, jpl_var_1602, jmco2_var_1603, jmco
  REAL(KIND = 8) :: refrat_planck_a_var_1604, refrat_m_a_var_1605, refrat_m_a3
  REAL(KIND = 8) :: fac000_var_1606, fac100_var_1607, fac200_var_1608, fac010_var_1609, fac110_var_1610, fac210_var_1611, fac001_var_1612, fac101_var_1613, fac201_var_1614, fac011_var_1615, fac111_var_1616, fac211_var_1617
  REAL(KIND = 8) :: p_var_1618, p4_var_1619, fk0_var_1620, fk1_var_1621, fk2_var_1622
  REAL(KIND = 8) :: taufor_var_1623, tauself_var_1624, tau_major_var_1625(4), tau_major1_var_1626(4), co2m1_var_1627, co2m2_var_1628, absco2_var_1629
  REAL(KIND = 8) :: com1, com2, absco, abso3_var_1630
  REAL(KIND = 8) :: chi_co2_var_1631, ratco2_var_1632, adjfac_var_1633, adjcolco2_var_1634
  REAL(KIND = 8) :: fs_var_1635, specmult_var_1636, specparm_var_1637, fs1_var_1638, specmult1_var_1639, specparm1_var_1640, fmco2_var_1641, specmult_mco2_var_1642, specparm_mco2_var_1643, fmco, specmult_mco, specparm_mco, fpl_var_1644, specmult_planck_var_1645, specparm_planck_var_1646
  REAL(KIND = 8) :: colco(kidia_var_1561 : kfdia_var_1562, klev_var_1563)
  INTEGER(KIND = 4) :: laytrop_min_var_1647, laytrop_max_var_1648
  INTEGER(KIND = 4) :: ixc_var_1649(klev_var_1563), ixlow_var_1650(kfdia_var_1562, klev_var_1563), ixhigh_var_1651(kfdia_var_1562, klev_var_1563)
  INTEGER(KIND = 4) :: ich_var_1652, icl_var_1653, ixc0_var_1654, ixp_var_1655, jc_var_1656, jl_var_1657
  DO lay_var_1600 = 1, klev_var_1563
    DO jc_var_1656 = kidia_var_1561, kfdia_var_1562
      colco(jc_var_1656, lay_var_1600) = 0.0D0
    END DO
  END DO
  laytrop_min_var_1647 = MINVAL(laytrop_var_1579)
  laytrop_max_var_1648 = MAXVAL(laytrop_var_1579)
  ixlow_var_1650 = 0
  ixhigh_var_1651 = 0
  ixc_var_1649 = 0
  DO lay_var_1600 = laytrop_min_var_1647 + 1, laytrop_max_var_1648
    icl_var_1653 = 0
    ich_var_1652 = 0
    DO jc_var_1656 = kidia_var_1561, kfdia_var_1562
      IF (lay_var_1600 <= laytrop_var_1579(jc_var_1656)) THEN
        icl_var_1653 = icl_var_1653 + 1
        ixlow_var_1650(icl_var_1653, lay_var_1600) = jc_var_1656
      ELSE
        ich_var_1652 = ich_var_1652 + 1
        ixhigh_var_1651(ich_var_1652, lay_var_1600) = jc_var_1656
      END IF
    END DO
    ixc_var_1649(lay_var_1600) = icl_var_1653
  END DO
  refrat_planck_a_var_1604 = chi_mls(1, 5) / chi_mls(4, 5)
  refrat_m_a_var_1605 = chi_mls(1, 1) / chi_mls(4, 1)
  refrat_m_a3 = chi_mls(1, 3) / chi_mls(4, 3)
  DO lay_var_1600 = 1, laytrop_min_var_1647
    DO jl_var_1657 = kidia_var_1561, kfdia_var_1562
      speccomb_var_1589 = colh2o_var_1574(jl_var_1657, lay_var_1600) + rat_h2on2o(jl_var_1657, lay_var_1600) * coln2o_var_1575(jl_var_1657, lay_var_1600)
      specparm_var_1637 = MIN(colh2o_var_1574(jl_var_1657, lay_var_1600) / speccomb_var_1589, oneminus_var_1573)
      specmult_var_1636 = 8.0D0 * (specparm_var_1637)
      js_var_1599 = 1 + INT(specmult_var_1636)
      fs_var_1635 = ((specmult_var_1636) - AINT((specmult_var_1636)))
      speccomb1_var_1590 = colh2o_var_1574(jl_var_1657, lay_var_1600) + rat_h2on2o_1(jl_var_1657, lay_var_1600) * coln2o_var_1575(jl_var_1657, lay_var_1600)
      specparm1_var_1640 = MIN(colh2o_var_1574(jl_var_1657, lay_var_1600) / speccomb1_var_1590, oneminus_var_1573)
      specmult1_var_1639 = 8.0D0 * (specparm1_var_1640)
      js1_var_1601 = 1 + INT(specmult1_var_1639)
      fs1_var_1638 = ((specmult1_var_1639) - AINT((specmult1_var_1639)))
      speccomb_mco2_var_1592 = colh2o_var_1574(jl_var_1657, lay_var_1600) + refrat_m_a_var_1605 * coln2o_var_1575(jl_var_1657, lay_var_1600)
      specparm_mco2_var_1643 = MIN(colh2o_var_1574(jl_var_1657, lay_var_1600) / speccomb_mco2_var_1592, oneminus_var_1573)
      specmult_mco2_var_1642 = 8.0D0 * specparm_mco2_var_1643
      jmco2_var_1603 = 1 + INT(specmult_mco2_var_1642)
      fmco2_var_1641 = ((specmult_mco2_var_1642) - AINT((specmult_mco2_var_1642)))
      chi_co2_var_1631 = colco2_var_1576(jl_var_1657, lay_var_1600) / (coldry_var_1578(jl_var_1657, lay_var_1600))
      ratco2_var_1632 = 1D+20 * chi_co2_var_1631 / 0.000355D0
      IF (ratco2_var_1632 .GT. 3.0D0) THEN
        adjfac_var_1633 = 2.0D0 + (ratco2_var_1632 - 2.0D0) ** 0.68D0
        adjcolco2_var_1634 = adjfac_var_1633 * 0.000355D0 * coldry_var_1578(jl_var_1657, lay_var_1600) * 1D-20
      ELSE
        adjcolco2_var_1634 = colco2_var_1576(jl_var_1657, lay_var_1600)
      END IF
      speccomb_mco = colh2o_var_1574(jl_var_1657, lay_var_1600) + refrat_m_a3 * coln2o_var_1575(jl_var_1657, lay_var_1600)
      specparm_mco = MIN(colh2o_var_1574(jl_var_1657, lay_var_1600) / speccomb_mco, oneminus_var_1573)
      specmult_mco = 8.0D0 * specparm_mco
      jmco = 1 + INT(specmult_mco)
      fmco = ((specmult_mco) - AINT((specmult_mco)))
      speccomb_planck_var_1591 = colh2o_var_1574(jl_var_1657, lay_var_1600) + refrat_planck_a_var_1604 * coln2o_var_1575(jl_var_1657, lay_var_1600)
      specparm_planck_var_1646 = MIN(colh2o_var_1574(jl_var_1657, lay_var_1600) / speccomb_planck_var_1591, oneminus_var_1573)
      specmult_planck_var_1645 = 8.0D0 * specparm_planck_var_1646
      jpl_var_1602 = 1 + INT(specmult_planck_var_1645)
      fpl_var_1644 = ((specmult_planck_var_1645) - AINT((specmult_planck_var_1645)))
      ind0_var_1593 = ((jp_var_1570(jl_var_1657, lay_var_1600) - 1) * 5 + (jt_var_1571(jl_var_1657, lay_var_1600) - 1)) * nspa_var_216(13) + js_var_1599
      ind1_var_1594 = (jp_var_1570(jl_var_1657, lay_var_1600) * 5 + (jt1_var_1572(jl_var_1657, lay_var_1600) - 1)) * nspa_var_216(13) + js1_var_1601
      inds_var_1595 = indself_var_1582(jl_var_1657, lay_var_1600)
      indf_var_1596 = indfor_var_1584(jl_var_1657, lay_var_1600)
      indm_var_1597 = indminor_var_1588(jl_var_1657, lay_var_1600)
      IF (specparm_var_1637 .LT. 0.125D0) THEN
        p_var_1618 = fs_var_1635 - 1.0D0
        p4_var_1619 = p_var_1618 ** 4
        fk0_var_1620 = p4_var_1619
        fk1_var_1621 = 1.0D0 - p_var_1618 - 2.0D0 * p4_var_1619
        fk2_var_1622 = p_var_1618 + p4_var_1619
        fac000_var_1606 = fk0_var_1620 * fac00_var_1566(jl_var_1657, lay_var_1600)
        fac100_var_1607 = fk1_var_1621 * fac00_var_1566(jl_var_1657, lay_var_1600)
        fac200_var_1608 = fk2_var_1622 * fac00_var_1566(jl_var_1657, lay_var_1600)
        fac010_var_1609 = fk0_var_1620 * fac10_var_1568(jl_var_1657, lay_var_1600)
        fac110_var_1610 = fk1_var_1621 * fac10_var_1568(jl_var_1657, lay_var_1600)
        fac210_var_1611 = fk2_var_1622 * fac10_var_1568(jl_var_1657, lay_var_1600)
      ELSE IF (specparm_var_1637 .GT. 0.875D0) THEN
        p_var_1618 = - fs_var_1635
        p4_var_1619 = p_var_1618 ** 4
        fk0_var_1620 = p4_var_1619
        fk1_var_1621 = 1.0D0 - p_var_1618 - 2.0D0 * p4_var_1619
        fk2_var_1622 = p_var_1618 + p4_var_1619
        fac000_var_1606 = fk0_var_1620 * fac00_var_1566(jl_var_1657, lay_var_1600)
        fac100_var_1607 = fk1_var_1621 * fac00_var_1566(jl_var_1657, lay_var_1600)
        fac200_var_1608 = fk2_var_1622 * fac00_var_1566(jl_var_1657, lay_var_1600)
        fac010_var_1609 = fk0_var_1620 * fac10_var_1568(jl_var_1657, lay_var_1600)
        fac110_var_1610 = fk1_var_1621 * fac10_var_1568(jl_var_1657, lay_var_1600)
        fac210_var_1611 = fk2_var_1622 * fac10_var_1568(jl_var_1657, lay_var_1600)
      ELSE
        fac000_var_1606 = (1.0D0 - fs_var_1635) * fac00_var_1566(jl_var_1657, lay_var_1600)
        fac010_var_1609 = (1.0D0 - fs_var_1635) * fac10_var_1568(jl_var_1657, lay_var_1600)
        fac100_var_1607 = fs_var_1635 * fac00_var_1566(jl_var_1657, lay_var_1600)
        fac110_var_1610 = fs_var_1635 * fac10_var_1568(jl_var_1657, lay_var_1600)
        fac200_var_1608 = 0.0D0
        fac210_var_1611 = 0.0D0
      END IF
      IF (specparm1_var_1640 .LT. 0.125D0) THEN
        p_var_1618 = fs1_var_1638 - 1.0D0
        p4_var_1619 = p_var_1618 ** 4
        fk0_var_1620 = p4_var_1619
        fk1_var_1621 = 1.0D0 - p_var_1618 - 2.0D0 * p4_var_1619
        fk2_var_1622 = p_var_1618 + p4_var_1619
        fac001_var_1612 = fk0_var_1620 * fac01_var_1567(jl_var_1657, lay_var_1600)
        fac101_var_1613 = fk1_var_1621 * fac01_var_1567(jl_var_1657, lay_var_1600)
        fac201_var_1614 = fk2_var_1622 * fac01_var_1567(jl_var_1657, lay_var_1600)
        fac011_var_1615 = fk0_var_1620 * fac11_var_1569(jl_var_1657, lay_var_1600)
        fac111_var_1616 = fk1_var_1621 * fac11_var_1569(jl_var_1657, lay_var_1600)
        fac211_var_1617 = fk2_var_1622 * fac11_var_1569(jl_var_1657, lay_var_1600)
      ELSE IF (specparm1_var_1640 .GT. 0.875D0) THEN
        p_var_1618 = - fs1_var_1638
        p4_var_1619 = p_var_1618 ** 4
        fk0_var_1620 = p4_var_1619
        fk1_var_1621 = 1.0D0 - p_var_1618 - 2.0D0 * p4_var_1619
        fk2_var_1622 = p_var_1618 + p4_var_1619
        fac001_var_1612 = fk0_var_1620 * fac01_var_1567(jl_var_1657, lay_var_1600)
        fac101_var_1613 = fk1_var_1621 * fac01_var_1567(jl_var_1657, lay_var_1600)
        fac201_var_1614 = fk2_var_1622 * fac01_var_1567(jl_var_1657, lay_var_1600)
        fac011_var_1615 = fk0_var_1620 * fac11_var_1569(jl_var_1657, lay_var_1600)
        fac111_var_1616 = fk1_var_1621 * fac11_var_1569(jl_var_1657, lay_var_1600)
        fac211_var_1617 = fk2_var_1622 * fac11_var_1569(jl_var_1657, lay_var_1600)
      ELSE
        fac001_var_1612 = (1.0D0 - fs1_var_1638) * fac01_var_1567(jl_var_1657, lay_var_1600)
        fac011_var_1615 = (1.0D0 - fs1_var_1638) * fac11_var_1569(jl_var_1657, lay_var_1600)
        fac101_var_1613 = fs1_var_1638 * fac01_var_1567(jl_var_1657, lay_var_1600)
        fac111_var_1616 = fs1_var_1638 * fac11_var_1569(jl_var_1657, lay_var_1600)
        fac201_var_1614 = 0.0D0
        fac211_var_1617 = 0.0D0
      END IF
      IF (specparm_var_1637 .LT. 0.125D0) THEN
        tau_major_var_1625(1 : ng13) = speccomb_var_1589 * (fac000_var_1606 * absa_var_132(ind0_var_1593, 1 : 4) + fac100_var_1607 * absa_var_132(ind0_var_1593 + 1, 1 : 4) + fac200_var_1608 * absa_var_132(ind0_var_1593 + 2, 1 : 4) + fac010_var_1609 * absa_var_132(ind0_var_1593 + 9, 1 : 4) + fac110_var_1610 * absa_var_132(ind0_var_1593 + 10, 1 : 4) + fac210_var_1611 * absa_var_132(ind0_var_1593 + 11, 1 : 4))
      ELSE IF (specparm_var_1637 .GT. 0.875D0) THEN
        tau_major_var_1625(1 : ng13) = speccomb_var_1589 * (fac200_var_1608 * absa_var_132(ind0_var_1593 - 1, 1 : 4) + fac100_var_1607 * absa_var_132(ind0_var_1593, 1 : 4) + fac000_var_1606 * absa_var_132(ind0_var_1593 + 1, 1 : 4) + fac210_var_1611 * absa_var_132(ind0_var_1593 + 8, 1 : 4) + fac110_var_1610 * absa_var_132(ind0_var_1593 + 9, 1 : 4) + fac010_var_1609 * absa_var_132(ind0_var_1593 + 10, 1 : 4))
      ELSE
        tau_major_var_1625(1 : ng13) = speccomb_var_1589 * (fac000_var_1606 * absa_var_132(ind0_var_1593, 1 : 4) + fac100_var_1607 * absa_var_132(ind0_var_1593 + 1, 1 : 4) + fac010_var_1609 * absa_var_132(ind0_var_1593 + 9, 1 : 4) + fac110_var_1610 * absa_var_132(ind0_var_1593 + 10, 1 : 4))
      END IF
      IF (specparm1_var_1640 .LT. 0.125D0) THEN
        tau_major1_var_1626(1 : ng13) = speccomb1_var_1590 * (fac001_var_1612 * absa_var_132(ind1_var_1594, 1 : 4) + fac101_var_1613 * absa_var_132(ind1_var_1594 + 1, 1 : 4) + fac201_var_1614 * absa_var_132(ind1_var_1594 + 2, 1 : 4) + fac011_var_1615 * absa_var_132(ind1_var_1594 + 9, 1 : 4) + fac111_var_1616 * absa_var_132(ind1_var_1594 + 10, 1 : 4) + fac211_var_1617 * absa_var_132(ind1_var_1594 + 11, 1 : 4))
      ELSE IF (specparm1_var_1640 .GT. 0.875D0) THEN
        tau_major1_var_1626(1 : ng13) = speccomb1_var_1590 * (fac201_var_1614 * absa_var_132(ind1_var_1594 - 1, 1 : 4) + fac101_var_1613 * absa_var_132(ind1_var_1594, 1 : 4) + fac001_var_1612 * absa_var_132(ind1_var_1594 + 1, 1 : 4) + fac211_var_1617 * absa_var_132(ind1_var_1594 + 8, 1 : 4) + fac111_var_1616 * absa_var_132(ind1_var_1594 + 9, 1 : 4) + fac011_var_1615 * absa_var_132(ind1_var_1594 + 10, 1 : 4))
      ELSE
        tau_major1_var_1626(1 : ng13) = speccomb1_var_1590 * (fac001_var_1612 * absa_var_132(ind1_var_1594, 1 : 4) + fac101_var_1613 * absa_var_132(ind1_var_1594 + 1, 1 : 4) + fac011_var_1615 * absa_var_132(ind1_var_1594 + 9, 1 : 4) + fac111_var_1616 * absa_var_132(ind1_var_1594 + 10, 1 : 4))
      END IF
      DO ig_var_1598 = 1, 4
        tauself_var_1624 = selffac_var_1580(jl_var_1657, lay_var_1600) * (selfref_var_133(inds_var_1595, ig_var_1598) + selffrac_var_1581(jl_var_1657, lay_var_1600) * (selfref_var_133(inds_var_1595 + 1, ig_var_1598) - selfref_var_133(inds_var_1595, ig_var_1598)))
        taufor_var_1623 = forfac_var_1585(jl_var_1657, lay_var_1600) * (forref_var_134(indf_var_1596, ig_var_1598) + forfrac_var_1586(jl_var_1657, lay_var_1600) * (forref_var_134(indf_var_1596 + 1, ig_var_1598) - forref_var_134(indf_var_1596, ig_var_1598)))
        co2m1_var_1627 = ka_mco2_var_135(jmco2_var_1603, indm_var_1597, ig_var_1598) + fmco2_var_1641 * (ka_mco2_var_135(jmco2_var_1603 + 1, indm_var_1597, ig_var_1598) - ka_mco2_var_135(jmco2_var_1603, indm_var_1597, ig_var_1598))
        co2m2_var_1628 = ka_mco2_var_135(jmco2_var_1603, indm_var_1597 + 1, ig_var_1598) + fmco2_var_1641 * (ka_mco2_var_135(jmco2_var_1603 + 1, indm_var_1597 + 1, ig_var_1598) - ka_mco2_var_135(jmco2_var_1603, indm_var_1597 + 1, ig_var_1598))
        absco2_var_1629 = co2m1_var_1627 + minorfrac_var_1587(jl_var_1657, lay_var_1600) * (co2m2_var_1628 - co2m1_var_1627)
        com1 = ka_mco(jmco, indm_var_1597, ig_var_1598) + fmco * (ka_mco(jmco + 1, indm_var_1597, ig_var_1598) - ka_mco(jmco, indm_var_1597, ig_var_1598))
        com2 = ka_mco(jmco, indm_var_1597 + 1, ig_var_1598) + fmco * (ka_mco(jmco + 1, indm_var_1597 + 1, ig_var_1598) - ka_mco(jmco, indm_var_1597 + 1, ig_var_1598))
        absco = com1 + minorfrac_var_1587(jl_var_1657, lay_var_1600) * (com2 - com1)
        taug_var_1564(jl_var_1657, 130 + ig_var_1598, lay_var_1600) = tau_major_var_1625(ig_var_1598) + tau_major1_var_1626(ig_var_1598) + tauself_var_1624 + taufor_var_1623 + adjcolco2_var_1634 * absco2_var_1629 + colco(jl_var_1657, lay_var_1600) * absco
        fracs_var_1583(jl_var_1657, 130 + ig_var_1598, lay_var_1600) = fracrefa_var_130(ig_var_1598, jpl_var_1602) + fpl_var_1644 * (fracrefa_var_130(ig_var_1598, jpl_var_1602 + 1) - fracrefa_var_130(ig_var_1598, jpl_var_1602))
      END DO
    END DO
  END DO
  DO lay_var_1600 = laytrop_max_var_1648 + 1, klev_var_1563
    DO jl_var_1657 = kidia_var_1561, kfdia_var_1562
      indm_var_1597 = indminor_var_1588(jl_var_1657, lay_var_1600)
      DO ig_var_1598 = 1, 4
        abso3_var_1630 = kb_mo3(indm_var_1597, ig_var_1598) + minorfrac_var_1587(jl_var_1657, lay_var_1600) * (kb_mo3(indm_var_1597 + 1, ig_var_1598) - kb_mo3(indm_var_1597, ig_var_1598))
        taug_var_1564(jl_var_1657, 130 + ig_var_1598, lay_var_1600) = colo3_var_1577(jl_var_1657, lay_var_1600) * abso3_var_1630
        fracs_var_1583(jl_var_1657, 130 + ig_var_1598, lay_var_1600) = fracrefb_var_131(ig_var_1598)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1648 /= laytrop_min_var_1647) THEN
    DO lay_var_1600 = laytrop_min_var_1647 + 1, laytrop_max_var_1648
      ixc0_var_1654 = ixc_var_1649(lay_var_1600)
      DO ixp_var_1655 = 1, ixc0_var_1654
        jl_var_1657 = ixlow_var_1650(ixp_var_1655, lay_var_1600)
        speccomb_var_1589 = colh2o_var_1574(jl_var_1657, lay_var_1600) + rat_h2on2o(jl_var_1657, lay_var_1600) * coln2o_var_1575(jl_var_1657, lay_var_1600)
        specparm_var_1637 = MIN(colh2o_var_1574(jl_var_1657, lay_var_1600) / speccomb_var_1589, oneminus_var_1573)
        specmult_var_1636 = 8.0D0 * (specparm_var_1637)
        js_var_1599 = 1 + INT(specmult_var_1636)
        fs_var_1635 = ((specmult_var_1636) - AINT((specmult_var_1636)))
        speccomb1_var_1590 = colh2o_var_1574(jl_var_1657, lay_var_1600) + rat_h2on2o_1(jl_var_1657, lay_var_1600) * coln2o_var_1575(jl_var_1657, lay_var_1600)
        specparm1_var_1640 = MIN(colh2o_var_1574(jl_var_1657, lay_var_1600) / speccomb1_var_1590, oneminus_var_1573)
        specmult1_var_1639 = 8.0D0 * (specparm1_var_1640)
        js1_var_1601 = 1 + INT(specmult1_var_1639)
        fs1_var_1638 = ((specmult1_var_1639) - AINT((specmult1_var_1639)))
        speccomb_mco2_var_1592 = colh2o_var_1574(jl_var_1657, lay_var_1600) + refrat_m_a_var_1605 * coln2o_var_1575(jl_var_1657, lay_var_1600)
        specparm_mco2_var_1643 = MIN(colh2o_var_1574(jl_var_1657, lay_var_1600) / speccomb_mco2_var_1592, oneminus_var_1573)
        specmult_mco2_var_1642 = 8.0D0 * specparm_mco2_var_1643
        jmco2_var_1603 = 1 + INT(specmult_mco2_var_1642)
        fmco2_var_1641 = ((specmult_mco2_var_1642) - AINT((specmult_mco2_var_1642)))
        chi_co2_var_1631 = colco2_var_1576(jl_var_1657, lay_var_1600) / (coldry_var_1578(jl_var_1657, lay_var_1600))
        ratco2_var_1632 = 1D+20 * chi_co2_var_1631 / 0.000355D0
        IF (ratco2_var_1632 .GT. 3.0D0) THEN
          adjfac_var_1633 = 2.0D0 + (ratco2_var_1632 - 2.0D0) ** 0.68D0
          adjcolco2_var_1634 = adjfac_var_1633 * 0.000355D0 * coldry_var_1578(jl_var_1657, lay_var_1600) * 1D-20
        ELSE
          adjcolco2_var_1634 = colco2_var_1576(jl_var_1657, lay_var_1600)
        END IF
        speccomb_mco = colh2o_var_1574(jl_var_1657, lay_var_1600) + refrat_m_a3 * coln2o_var_1575(jl_var_1657, lay_var_1600)
        specparm_mco = MIN(colh2o_var_1574(jl_var_1657, lay_var_1600) / speccomb_mco, oneminus_var_1573)
        specmult_mco = 8.0D0 * specparm_mco
        jmco = 1 + INT(specmult_mco)
        fmco = ((specmult_mco) - AINT((specmult_mco)))
        speccomb_planck_var_1591 = colh2o_var_1574(jl_var_1657, lay_var_1600) + refrat_planck_a_var_1604 * coln2o_var_1575(jl_var_1657, lay_var_1600)
        specparm_planck_var_1646 = MIN(colh2o_var_1574(jl_var_1657, lay_var_1600) / speccomb_planck_var_1591, oneminus_var_1573)
        specmult_planck_var_1645 = 8.0D0 * specparm_planck_var_1646
        jpl_var_1602 = 1 + INT(specmult_planck_var_1645)
        fpl_var_1644 = ((specmult_planck_var_1645) - AINT((specmult_planck_var_1645)))
        ind0_var_1593 = ((jp_var_1570(jl_var_1657, lay_var_1600) - 1) * 5 + (jt_var_1571(jl_var_1657, lay_var_1600) - 1)) * nspa_var_216(13) + js_var_1599
        ind1_var_1594 = (jp_var_1570(jl_var_1657, lay_var_1600) * 5 + (jt1_var_1572(jl_var_1657, lay_var_1600) - 1)) * nspa_var_216(13) + js1_var_1601
        inds_var_1595 = indself_var_1582(jl_var_1657, lay_var_1600)
        indf_var_1596 = indfor_var_1584(jl_var_1657, lay_var_1600)
        indm_var_1597 = indminor_var_1588(jl_var_1657, lay_var_1600)
        IF (specparm_var_1637 .LT. 0.125D0) THEN
          p_var_1618 = fs_var_1635 - 1.0D0
          p4_var_1619 = p_var_1618 ** 4
          fk0_var_1620 = p4_var_1619
          fk1_var_1621 = 1.0D0 - p_var_1618 - 2.0D0 * p4_var_1619
          fk2_var_1622 = p_var_1618 + p4_var_1619
          fac000_var_1606 = fk0_var_1620 * fac00_var_1566(jl_var_1657, lay_var_1600)
          fac100_var_1607 = fk1_var_1621 * fac00_var_1566(jl_var_1657, lay_var_1600)
          fac200_var_1608 = fk2_var_1622 * fac00_var_1566(jl_var_1657, lay_var_1600)
          fac010_var_1609 = fk0_var_1620 * fac10_var_1568(jl_var_1657, lay_var_1600)
          fac110_var_1610 = fk1_var_1621 * fac10_var_1568(jl_var_1657, lay_var_1600)
          fac210_var_1611 = fk2_var_1622 * fac10_var_1568(jl_var_1657, lay_var_1600)
        ELSE IF (specparm_var_1637 .GT. 0.875D0) THEN
          p_var_1618 = - fs_var_1635
          p4_var_1619 = p_var_1618 ** 4
          fk0_var_1620 = p4_var_1619
          fk1_var_1621 = 1.0D0 - p_var_1618 - 2.0D0 * p4_var_1619
          fk2_var_1622 = p_var_1618 + p4_var_1619
          fac000_var_1606 = fk0_var_1620 * fac00_var_1566(jl_var_1657, lay_var_1600)
          fac100_var_1607 = fk1_var_1621 * fac00_var_1566(jl_var_1657, lay_var_1600)
          fac200_var_1608 = fk2_var_1622 * fac00_var_1566(jl_var_1657, lay_var_1600)
          fac010_var_1609 = fk0_var_1620 * fac10_var_1568(jl_var_1657, lay_var_1600)
          fac110_var_1610 = fk1_var_1621 * fac10_var_1568(jl_var_1657, lay_var_1600)
          fac210_var_1611 = fk2_var_1622 * fac10_var_1568(jl_var_1657, lay_var_1600)
        ELSE
          fac000_var_1606 = (1.0D0 - fs_var_1635) * fac00_var_1566(jl_var_1657, lay_var_1600)
          fac010_var_1609 = (1.0D0 - fs_var_1635) * fac10_var_1568(jl_var_1657, lay_var_1600)
          fac100_var_1607 = fs_var_1635 * fac00_var_1566(jl_var_1657, lay_var_1600)
          fac110_var_1610 = fs_var_1635 * fac10_var_1568(jl_var_1657, lay_var_1600)
          fac200_var_1608 = 0.0D0
          fac210_var_1611 = 0.0D0
        END IF
        IF (specparm1_var_1640 .LT. 0.125D0) THEN
          p_var_1618 = fs1_var_1638 - 1.0D0
          p4_var_1619 = p_var_1618 ** 4
          fk0_var_1620 = p4_var_1619
          fk1_var_1621 = 1.0D0 - p_var_1618 - 2.0D0 * p4_var_1619
          fk2_var_1622 = p_var_1618 + p4_var_1619
          fac001_var_1612 = fk0_var_1620 * fac01_var_1567(jl_var_1657, lay_var_1600)
          fac101_var_1613 = fk1_var_1621 * fac01_var_1567(jl_var_1657, lay_var_1600)
          fac201_var_1614 = fk2_var_1622 * fac01_var_1567(jl_var_1657, lay_var_1600)
          fac011_var_1615 = fk0_var_1620 * fac11_var_1569(jl_var_1657, lay_var_1600)
          fac111_var_1616 = fk1_var_1621 * fac11_var_1569(jl_var_1657, lay_var_1600)
          fac211_var_1617 = fk2_var_1622 * fac11_var_1569(jl_var_1657, lay_var_1600)
        ELSE IF (specparm1_var_1640 .GT. 0.875D0) THEN
          p_var_1618 = - fs1_var_1638
          p4_var_1619 = p_var_1618 ** 4
          fk0_var_1620 = p4_var_1619
          fk1_var_1621 = 1.0D0 - p_var_1618 - 2.0D0 * p4_var_1619
          fk2_var_1622 = p_var_1618 + p4_var_1619
          fac001_var_1612 = fk0_var_1620 * fac01_var_1567(jl_var_1657, lay_var_1600)
          fac101_var_1613 = fk1_var_1621 * fac01_var_1567(jl_var_1657, lay_var_1600)
          fac201_var_1614 = fk2_var_1622 * fac01_var_1567(jl_var_1657, lay_var_1600)
          fac011_var_1615 = fk0_var_1620 * fac11_var_1569(jl_var_1657, lay_var_1600)
          fac111_var_1616 = fk1_var_1621 * fac11_var_1569(jl_var_1657, lay_var_1600)
          fac211_var_1617 = fk2_var_1622 * fac11_var_1569(jl_var_1657, lay_var_1600)
        ELSE
          fac001_var_1612 = (1.0D0 - fs1_var_1638) * fac01_var_1567(jl_var_1657, lay_var_1600)
          fac011_var_1615 = (1.0D0 - fs1_var_1638) * fac11_var_1569(jl_var_1657, lay_var_1600)
          fac101_var_1613 = fs1_var_1638 * fac01_var_1567(jl_var_1657, lay_var_1600)
          fac111_var_1616 = fs1_var_1638 * fac11_var_1569(jl_var_1657, lay_var_1600)
          fac201_var_1614 = 0.0D0
          fac211_var_1617 = 0.0D0
        END IF
        IF (specparm_var_1637 .LT. 0.125D0) THEN
          tau_major_var_1625(1 : ng13) = speccomb_var_1589 * (fac000_var_1606 * absa_var_132(ind0_var_1593, 1 : 4) + fac100_var_1607 * absa_var_132(ind0_var_1593 + 1, 1 : 4) + fac200_var_1608 * absa_var_132(ind0_var_1593 + 2, 1 : 4) + fac010_var_1609 * absa_var_132(ind0_var_1593 + 9, 1 : 4) + fac110_var_1610 * absa_var_132(ind0_var_1593 + 10, 1 : 4) + fac210_var_1611 * absa_var_132(ind0_var_1593 + 11, 1 : 4))
        ELSE IF (specparm_var_1637 .GT. 0.875D0) THEN
          tau_major_var_1625(1 : ng13) = speccomb_var_1589 * (fac200_var_1608 * absa_var_132(ind0_var_1593 - 1, 1 : 4) + fac100_var_1607 * absa_var_132(ind0_var_1593, 1 : 4) + fac000_var_1606 * absa_var_132(ind0_var_1593 + 1, 1 : 4) + fac210_var_1611 * absa_var_132(ind0_var_1593 + 8, 1 : 4) + fac110_var_1610 * absa_var_132(ind0_var_1593 + 9, 1 : 4) + fac010_var_1609 * absa_var_132(ind0_var_1593 + 10, 1 : 4))
        ELSE
          tau_major_var_1625(1 : ng13) = speccomb_var_1589 * (fac000_var_1606 * absa_var_132(ind0_var_1593, 1 : 4) + fac100_var_1607 * absa_var_132(ind0_var_1593 + 1, 1 : 4) + fac010_var_1609 * absa_var_132(ind0_var_1593 + 9, 1 : 4) + fac110_var_1610 * absa_var_132(ind0_var_1593 + 10, 1 : 4))
        END IF
        IF (specparm1_var_1640 .LT. 0.125D0) THEN
          tau_major1_var_1626(1 : ng13) = speccomb1_var_1590 * (fac001_var_1612 * absa_var_132(ind1_var_1594, 1 : 4) + fac101_var_1613 * absa_var_132(ind1_var_1594 + 1, 1 : 4) + fac201_var_1614 * absa_var_132(ind1_var_1594 + 2, 1 : 4) + fac011_var_1615 * absa_var_132(ind1_var_1594 + 9, 1 : 4) + fac111_var_1616 * absa_var_132(ind1_var_1594 + 10, 1 : 4) + fac211_var_1617 * absa_var_132(ind1_var_1594 + 11, 1 : 4))
        ELSE IF (specparm1_var_1640 .GT. 0.875D0) THEN
          tau_major1_var_1626(1 : ng13) = speccomb1_var_1590 * (fac201_var_1614 * absa_var_132(ind1_var_1594 - 1, 1 : 4) + fac101_var_1613 * absa_var_132(ind1_var_1594, 1 : 4) + fac001_var_1612 * absa_var_132(ind1_var_1594 + 1, 1 : 4) + fac211_var_1617 * absa_var_132(ind1_var_1594 + 8, 1 : 4) + fac111_var_1616 * absa_var_132(ind1_var_1594 + 9, 1 : 4) + fac011_var_1615 * absa_var_132(ind1_var_1594 + 10, 1 : 4))
        ELSE
          tau_major1_var_1626(1 : ng13) = speccomb1_var_1590 * (fac001_var_1612 * absa_var_132(ind1_var_1594, 1 : 4) + fac101_var_1613 * absa_var_132(ind1_var_1594 + 1, 1 : 4) + fac011_var_1615 * absa_var_132(ind1_var_1594 + 9, 1 : 4) + fac111_var_1616 * absa_var_132(ind1_var_1594 + 10, 1 : 4))
        END IF
        DO ig_var_1598 = 1, 4
          tauself_var_1624 = selffac_var_1580(jl_var_1657, lay_var_1600) * (selfref_var_133(inds_var_1595, ig_var_1598) + selffrac_var_1581(jl_var_1657, lay_var_1600) * (selfref_var_133(inds_var_1595 + 1, ig_var_1598) - selfref_var_133(inds_var_1595, ig_var_1598)))
          taufor_var_1623 = forfac_var_1585(jl_var_1657, lay_var_1600) * (forref_var_134(indf_var_1596, ig_var_1598) + forfrac_var_1586(jl_var_1657, lay_var_1600) * (forref_var_134(indf_var_1596 + 1, ig_var_1598) - forref_var_134(indf_var_1596, ig_var_1598)))
          co2m1_var_1627 = ka_mco2_var_135(jmco2_var_1603, indm_var_1597, ig_var_1598) + fmco2_var_1641 * (ka_mco2_var_135(jmco2_var_1603 + 1, indm_var_1597, ig_var_1598) - ka_mco2_var_135(jmco2_var_1603, indm_var_1597, ig_var_1598))
          co2m2_var_1628 = ka_mco2_var_135(jmco2_var_1603, indm_var_1597 + 1, ig_var_1598) + fmco2_var_1641 * (ka_mco2_var_135(jmco2_var_1603 + 1, indm_var_1597 + 1, ig_var_1598) - ka_mco2_var_135(jmco2_var_1603, indm_var_1597 + 1, ig_var_1598))
          absco2_var_1629 = co2m1_var_1627 + minorfrac_var_1587(jl_var_1657, lay_var_1600) * (co2m2_var_1628 - co2m1_var_1627)
          com1 = ka_mco(jmco, indm_var_1597, ig_var_1598) + fmco * (ka_mco(jmco + 1, indm_var_1597, ig_var_1598) - ka_mco(jmco, indm_var_1597, ig_var_1598))
          com2 = ka_mco(jmco, indm_var_1597 + 1, ig_var_1598) + fmco * (ka_mco(jmco + 1, indm_var_1597 + 1, ig_var_1598) - ka_mco(jmco, indm_var_1597 + 1, ig_var_1598))
          absco = com1 + minorfrac_var_1587(jl_var_1657, lay_var_1600) * (com2 - com1)
          taug_var_1564(jl_var_1657, 130 + ig_var_1598, lay_var_1600) = tau_major_var_1625(ig_var_1598) + tau_major1_var_1626(ig_var_1598) + tauself_var_1624 + taufor_var_1623 + adjcolco2_var_1634 * absco2_var_1629 + colco(jl_var_1657, lay_var_1600) * absco
          fracs_var_1583(jl_var_1657, 130 + ig_var_1598, lay_var_1600) = fracrefa_var_130(ig_var_1598, jpl_var_1602) + fpl_var_1644 * (fracrefa_var_130(ig_var_1598, jpl_var_1602 + 1) - fracrefa_var_130(ig_var_1598, jpl_var_1602))
        END DO
      END DO
      ixc0_var_1654 = kfdia_var_1562 - kidia_var_1561 + 1 - ixc0_var_1654
      DO ixp_var_1655 = 1, ixc0_var_1654
        jl_var_1657 = ixhigh_var_1651(ixp_var_1655, lay_var_1600)
        indm_var_1597 = indminor_var_1588(jl_var_1657, lay_var_1600)
        DO ig_var_1598 = 1, 4
          abso3_var_1630 = kb_mo3(indm_var_1597, ig_var_1598) + minorfrac_var_1587(jl_var_1657, lay_var_1600) * (kb_mo3(indm_var_1597 + 1, ig_var_1598) - kb_mo3(indm_var_1597, ig_var_1598))
          taug_var_1564(jl_var_1657, 130 + ig_var_1598, lay_var_1600) = colo3_var_1577(jl_var_1657, lay_var_1600) * abso3_var_1630
          fracs_var_1583(jl_var_1657, 130 + ig_var_1598, lay_var_1600) = fracrefb_var_131(ig_var_1598)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol13
SUBROUTINE rrtm_taumol8(kidia_var_1658, kfdia_var_1659, klev_var_1660, taug_var_1661, wx_var_1662, p_tauaerl_var_1663, fac00_var_1664, fac01_var_1665, fac10_var_1666, fac11_var_1667, forfac_var_1683, forfrac_var_1682, indfor_var_1681, jp_var_1668, jt_var_1669, jt1_var_1670, colh2o_var_1671, colo3_var_1672, coln2o_var_1673, colco2_var_1674, coldry_var_1675, laytrop_var_1676, selffac_var_1677, selffrac_var_1678, indself_var_1679, fracs_var_1680, minorfrac_var_1684, indminor_var_1685)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta8, ONLY: absa_var_197, absb_var_198, cfc12_var_196, cfc22adj, forref_var_205, fracrefa_var_194, fracrefb_var_195, ka_mco2_var_199, ka_mn2o_var_200, ka_mo3_var_201, kb_mco2_var_202, kb_mn2o_var_203, selfref_var_204
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1658
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1659
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1660
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1661(kidia_var_1658 : kfdia_var_1659, 140, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: wx_var_1662(kidia_var_1658 : kfdia_var_1659, 4, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1663(kidia_var_1658 : kfdia_var_1659, klev_var_1660, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1664(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1665(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1666(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1667(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1668(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1669(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1670(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1671(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1672(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_1673(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1674(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1675(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1676(kidia_var_1658 : kfdia_var_1659)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1677(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1678(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1679(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1680(kidia_var_1658 : kfdia_var_1659, 140, klev_var_1660)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1681(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1682(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1683(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1684(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1685(kidia_var_1658 : kfdia_var_1659, klev_var_1660)
  INTEGER(KIND = 4) :: ind0_var_1686, ind1_var_1687, inds_var_1688, indf_var_1689, indm_var_1690
  INTEGER(KIND = 4) :: ig_var_1691, lay_var_1692
  REAL(KIND = 8) :: chi_co2_var_1693, ratco2_var_1694, adjfac_var_1695, adjcolco2_var_1696
  REAL(KIND = 8) :: taufor_var_1697, tauself_var_1698, abso3_var_1699, absco2_var_1700, absn2o_var_1701
  INTEGER(KIND = 4) :: laytrop_min_var_1702, laytrop_max_var_1703
  INTEGER(KIND = 4) :: ixc_var_1704(klev_var_1660), ixlow_var_1705(kfdia_var_1659, klev_var_1660), ixhigh_var_1706(kfdia_var_1659, klev_var_1660)
  INTEGER(KIND = 4) :: ich_var_1707, icl_var_1708, ixc0_var_1709, ixp_var_1710, jc_var_1711, jl_var_1712
  laytrop_min_var_1702 = MINVAL(laytrop_var_1676)
  laytrop_max_var_1703 = MAXVAL(laytrop_var_1676)
  ixlow_var_1705 = 0
  ixhigh_var_1706 = 0
  ixc_var_1704 = 0
  DO lay_var_1692 = laytrop_min_var_1702 + 1, laytrop_max_var_1703
    icl_var_1708 = 0
    ich_var_1707 = 0
    DO jc_var_1711 = kidia_var_1658, kfdia_var_1659
      IF (lay_var_1692 <= laytrop_var_1676(jc_var_1711)) THEN
        icl_var_1708 = icl_var_1708 + 1
        ixlow_var_1705(icl_var_1708, lay_var_1692) = jc_var_1711
      ELSE
        ich_var_1707 = ich_var_1707 + 1
        ixhigh_var_1706(ich_var_1707, lay_var_1692) = jc_var_1711
      END IF
    END DO
    ixc_var_1704(lay_var_1692) = icl_var_1708
  END DO
  DO lay_var_1692 = 1, laytrop_min_var_1702
    DO jl_var_1712 = kidia_var_1658, kfdia_var_1659
      chi_co2_var_1693 = colco2_var_1674(jl_var_1712, lay_var_1692) / (coldry_var_1675(jl_var_1712, lay_var_1692))
      ratco2_var_1694 = 1D+20 * chi_co2_var_1693 / chi_mls(2, jp_var_1668(jl_var_1712, lay_var_1692) + 1)
      IF (ratco2_var_1694 .GT. 3.0D0) THEN
        adjfac_var_1695 = 2.0D0 + (ratco2_var_1694 - 2.0D0) ** 0.65D0
        adjcolco2_var_1696 = adjfac_var_1695 * chi_mls(2, jp_var_1668(jl_var_1712, lay_var_1692) + 1) * coldry_var_1675(jl_var_1712, lay_var_1692) * 1D-20
      ELSE
        adjcolco2_var_1696 = colco2_var_1674(jl_var_1712, lay_var_1692)
      END IF
      ind0_var_1686 = ((jp_var_1668(jl_var_1712, lay_var_1692) - 1) * 5 + (jt_var_1669(jl_var_1712, lay_var_1692) - 1)) * nspa_var_216(8) + 1
      ind1_var_1687 = (jp_var_1668(jl_var_1712, lay_var_1692) * 5 + (jt1_var_1670(jl_var_1712, lay_var_1692) - 1)) * nspa_var_216(8) + 1
      inds_var_1688 = indself_var_1679(jl_var_1712, lay_var_1692)
      indf_var_1689 = indfor_var_1681(jl_var_1712, lay_var_1692)
      indm_var_1690 = indminor_var_1685(jl_var_1712, lay_var_1692)
      DO ig_var_1691 = 1, 8
        tauself_var_1698 = selffac_var_1677(jl_var_1712, lay_var_1692) * (selfref_var_204(inds_var_1688, ig_var_1691) + selffrac_var_1678(jl_var_1712, lay_var_1692) * (selfref_var_204(inds_var_1688 + 1, ig_var_1691) - selfref_var_204(inds_var_1688, ig_var_1691)))
        taufor_var_1697 = forfac_var_1683(jl_var_1712, lay_var_1692) * (forref_var_205(indf_var_1689, ig_var_1691) + forfrac_var_1682(jl_var_1712, lay_var_1692) * (forref_var_205(indf_var_1689 + 1, ig_var_1691) - forref_var_205(indf_var_1689, ig_var_1691)))
        absco2_var_1700 = (ka_mco2_var_199(indm_var_1690, ig_var_1691) + minorfrac_var_1684(jl_var_1712, lay_var_1692) * (ka_mco2_var_199(indm_var_1690 + 1, ig_var_1691) - ka_mco2_var_199(indm_var_1690, ig_var_1691)))
        abso3_var_1699 = (ka_mo3_var_201(indm_var_1690, ig_var_1691) + minorfrac_var_1684(jl_var_1712, lay_var_1692) * (ka_mo3_var_201(indm_var_1690 + 1, ig_var_1691) - ka_mo3_var_201(indm_var_1690, ig_var_1691)))
        absn2o_var_1701 = (ka_mn2o_var_200(indm_var_1690, ig_var_1691) + minorfrac_var_1684(jl_var_1712, lay_var_1692) * (ka_mn2o_var_200(indm_var_1690 + 1, ig_var_1691) - ka_mn2o_var_200(indm_var_1690, ig_var_1691)))
        taug_var_1661(jl_var_1712, 88 + ig_var_1691, lay_var_1692) = colh2o_var_1671(jl_var_1712, lay_var_1692) * (fac00_var_1664(jl_var_1712, lay_var_1692) * absa_var_197(ind0_var_1686, ig_var_1691) + fac10_var_1666(jl_var_1712, lay_var_1692) * absa_var_197(ind0_var_1686 + 1, ig_var_1691) + fac01_var_1665(jl_var_1712, lay_var_1692) * absa_var_197(ind1_var_1687, ig_var_1691) + fac11_var_1667(jl_var_1712, lay_var_1692) * absa_var_197(ind1_var_1687 + 1, ig_var_1691)) + tauself_var_1698 + taufor_var_1697 + adjcolco2_var_1696 * absco2_var_1700 + colo3_var_1672(jl_var_1712, lay_var_1692) * abso3_var_1699 + coln2o_var_1673(jl_var_1712, lay_var_1692) * absn2o_var_1701 + wx_var_1662(jl_var_1712, 3, lay_var_1692) * cfc12_var_196(ig_var_1691) + wx_var_1662(jl_var_1712, 4, lay_var_1692) * cfc22adj(ig_var_1691)
        fracs_var_1680(jl_var_1712, 88 + ig_var_1691, lay_var_1692) = fracrefa_var_194(ig_var_1691)
      END DO
    END DO
  END DO
  DO lay_var_1692 = laytrop_max_var_1703 + 1, klev_var_1660
    DO jl_var_1712 = kidia_var_1658, kfdia_var_1659
      chi_co2_var_1693 = colco2_var_1674(jl_var_1712, lay_var_1692) / coldry_var_1675(jl_var_1712, lay_var_1692)
      ratco2_var_1694 = 1D+20 * chi_co2_var_1693 / chi_mls(2, jp_var_1668(jl_var_1712, lay_var_1692) + 1)
      IF (ratco2_var_1694 .GT. 3.0D0) THEN
        adjfac_var_1695 = 2.0D0 + (ratco2_var_1694 - 2.0D0) ** 0.65D0
        adjcolco2_var_1696 = adjfac_var_1695 * chi_mls(2, jp_var_1668(jl_var_1712, lay_var_1692) + 1) * coldry_var_1675(jl_var_1712, lay_var_1692) * 1D-20
      ELSE
        adjcolco2_var_1696 = colco2_var_1674(jl_var_1712, lay_var_1692)
      END IF
      ind0_var_1686 = ((jp_var_1668(jl_var_1712, lay_var_1692) - 13) * 5 + (jt_var_1669(jl_var_1712, lay_var_1692) - 1)) * nspb_var_217(8) + 1
      ind1_var_1687 = ((jp_var_1668(jl_var_1712, lay_var_1692) - 12) * 5 + (jt1_var_1670(jl_var_1712, lay_var_1692) - 1)) * nspb_var_217(8) + 1
      indm_var_1690 = indminor_var_1685(jl_var_1712, lay_var_1692)
      DO ig_var_1691 = 1, 8
        absco2_var_1700 = (kb_mco2_var_202(indm_var_1690, ig_var_1691) + minorfrac_var_1684(jl_var_1712, lay_var_1692) * (kb_mco2_var_202(indm_var_1690 + 1, ig_var_1691) - kb_mco2_var_202(indm_var_1690, ig_var_1691)))
        absn2o_var_1701 = (kb_mn2o_var_203(indm_var_1690, ig_var_1691) + minorfrac_var_1684(jl_var_1712, lay_var_1692) * (kb_mn2o_var_203(indm_var_1690 + 1, ig_var_1691) - kb_mn2o_var_203(indm_var_1690, ig_var_1691)))
        taug_var_1661(jl_var_1712, 88 + ig_var_1691, lay_var_1692) = colo3_var_1672(jl_var_1712, lay_var_1692) * (fac00_var_1664(jl_var_1712, lay_var_1692) * absb_var_198(ind0_var_1686, ig_var_1691) + fac10_var_1666(jl_var_1712, lay_var_1692) * absb_var_198(ind0_var_1686 + 1, ig_var_1691) + fac01_var_1665(jl_var_1712, lay_var_1692) * absb_var_198(ind1_var_1687, ig_var_1691) + fac11_var_1667(jl_var_1712, lay_var_1692) * absb_var_198(ind1_var_1687 + 1, ig_var_1691)) + adjcolco2_var_1696 * absco2_var_1700 + coln2o_var_1673(jl_var_1712, lay_var_1692) * absn2o_var_1701 + wx_var_1662(jl_var_1712, 3, lay_var_1692) * cfc12_var_196(ig_var_1691) + wx_var_1662(jl_var_1712, 4, lay_var_1692) * cfc22adj(ig_var_1691)
        fracs_var_1680(jl_var_1712, 88 + ig_var_1691, lay_var_1692) = fracrefb_var_195(ig_var_1691)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1703 /= laytrop_min_var_1702) THEN
    DO lay_var_1692 = laytrop_min_var_1702 + 1, laytrop_max_var_1703
      ixc0_var_1709 = ixc_var_1704(lay_var_1692)
      DO ixp_var_1710 = 1, ixc0_var_1709
        jl_var_1712 = ixlow_var_1705(ixp_var_1710, lay_var_1692)
        chi_co2_var_1693 = colco2_var_1674(jl_var_1712, lay_var_1692) / (coldry_var_1675(jl_var_1712, lay_var_1692))
        ratco2_var_1694 = 1D+20 * chi_co2_var_1693 / chi_mls(2, jp_var_1668(jl_var_1712, lay_var_1692) + 1)
        IF (ratco2_var_1694 .GT. 3.0D0) THEN
          adjfac_var_1695 = 2.0D0 + (ratco2_var_1694 - 2.0D0) ** 0.65D0
          adjcolco2_var_1696 = adjfac_var_1695 * chi_mls(2, jp_var_1668(jl_var_1712, lay_var_1692) + 1) * coldry_var_1675(jl_var_1712, lay_var_1692) * 1D-20
        ELSE
          adjcolco2_var_1696 = colco2_var_1674(jl_var_1712, lay_var_1692)
        END IF
        ind0_var_1686 = ((jp_var_1668(jl_var_1712, lay_var_1692) - 1) * 5 + (jt_var_1669(jl_var_1712, lay_var_1692) - 1)) * nspa_var_216(8) + 1
        ind1_var_1687 = (jp_var_1668(jl_var_1712, lay_var_1692) * 5 + (jt1_var_1670(jl_var_1712, lay_var_1692) - 1)) * nspa_var_216(8) + 1
        inds_var_1688 = indself_var_1679(jl_var_1712, lay_var_1692)
        indf_var_1689 = indfor_var_1681(jl_var_1712, lay_var_1692)
        indm_var_1690 = indminor_var_1685(jl_var_1712, lay_var_1692)
        DO ig_var_1691 = 1, 8
          tauself_var_1698 = selffac_var_1677(jl_var_1712, lay_var_1692) * (selfref_var_204(inds_var_1688, ig_var_1691) + selffrac_var_1678(jl_var_1712, lay_var_1692) * (selfref_var_204(inds_var_1688 + 1, ig_var_1691) - selfref_var_204(inds_var_1688, ig_var_1691)))
          taufor_var_1697 = forfac_var_1683(jl_var_1712, lay_var_1692) * (forref_var_205(indf_var_1689, ig_var_1691) + forfrac_var_1682(jl_var_1712, lay_var_1692) * (forref_var_205(indf_var_1689 + 1, ig_var_1691) - forref_var_205(indf_var_1689, ig_var_1691)))
          absco2_var_1700 = (ka_mco2_var_199(indm_var_1690, ig_var_1691) + minorfrac_var_1684(jl_var_1712, lay_var_1692) * (ka_mco2_var_199(indm_var_1690 + 1, ig_var_1691) - ka_mco2_var_199(indm_var_1690, ig_var_1691)))
          abso3_var_1699 = (ka_mo3_var_201(indm_var_1690, ig_var_1691) + minorfrac_var_1684(jl_var_1712, lay_var_1692) * (ka_mo3_var_201(indm_var_1690 + 1, ig_var_1691) - ka_mo3_var_201(indm_var_1690, ig_var_1691)))
          absn2o_var_1701 = (ka_mn2o_var_200(indm_var_1690, ig_var_1691) + minorfrac_var_1684(jl_var_1712, lay_var_1692) * (ka_mn2o_var_200(indm_var_1690 + 1, ig_var_1691) - ka_mn2o_var_200(indm_var_1690, ig_var_1691)))
          taug_var_1661(jl_var_1712, 88 + ig_var_1691, lay_var_1692) = colh2o_var_1671(jl_var_1712, lay_var_1692) * (fac00_var_1664(jl_var_1712, lay_var_1692) * absa_var_197(ind0_var_1686, ig_var_1691) + fac10_var_1666(jl_var_1712, lay_var_1692) * absa_var_197(ind0_var_1686 + 1, ig_var_1691) + fac01_var_1665(jl_var_1712, lay_var_1692) * absa_var_197(ind1_var_1687, ig_var_1691) + fac11_var_1667(jl_var_1712, lay_var_1692) * absa_var_197(ind1_var_1687 + 1, ig_var_1691)) + tauself_var_1698 + taufor_var_1697 + adjcolco2_var_1696 * absco2_var_1700 + colo3_var_1672(jl_var_1712, lay_var_1692) * abso3_var_1699 + coln2o_var_1673(jl_var_1712, lay_var_1692) * absn2o_var_1701 + wx_var_1662(jl_var_1712, 3, lay_var_1692) * cfc12_var_196(ig_var_1691) + wx_var_1662(jl_var_1712, 4, lay_var_1692) * cfc22adj(ig_var_1691)
          fracs_var_1680(jl_var_1712, 88 + ig_var_1691, lay_var_1692) = fracrefa_var_194(ig_var_1691)
        END DO
      END DO
      ixc0_var_1709 = kfdia_var_1659 - kidia_var_1658 + 1 - ixc0_var_1709
      DO ixp_var_1710 = 1, ixc0_var_1709
        jl_var_1712 = ixhigh_var_1706(ixp_var_1710, lay_var_1692)
        chi_co2_var_1693 = colco2_var_1674(jl_var_1712, lay_var_1692) / coldry_var_1675(jl_var_1712, lay_var_1692)
        ratco2_var_1694 = 1D+20 * chi_co2_var_1693 / chi_mls(2, jp_var_1668(jl_var_1712, lay_var_1692) + 1)
        IF (ratco2_var_1694 .GT. 3.0D0) THEN
          adjfac_var_1695 = 2.0D0 + (ratco2_var_1694 - 2.0D0) ** 0.65D0
          adjcolco2_var_1696 = adjfac_var_1695 * chi_mls(2, jp_var_1668(jl_var_1712, lay_var_1692) + 1) * coldry_var_1675(jl_var_1712, lay_var_1692) * 1D-20
        ELSE
          adjcolco2_var_1696 = colco2_var_1674(jl_var_1712, lay_var_1692)
        END IF
        ind0_var_1686 = ((jp_var_1668(jl_var_1712, lay_var_1692) - 13) * 5 + (jt_var_1669(jl_var_1712, lay_var_1692) - 1)) * nspb_var_217(8) + 1
        ind1_var_1687 = ((jp_var_1668(jl_var_1712, lay_var_1692) - 12) * 5 + (jt1_var_1670(jl_var_1712, lay_var_1692) - 1)) * nspb_var_217(8) + 1
        indm_var_1690 = indminor_var_1685(jl_var_1712, lay_var_1692)
        DO ig_var_1691 = 1, 8
          absco2_var_1700 = (kb_mco2_var_202(indm_var_1690, ig_var_1691) + minorfrac_var_1684(jl_var_1712, lay_var_1692) * (kb_mco2_var_202(indm_var_1690 + 1, ig_var_1691) - kb_mco2_var_202(indm_var_1690, ig_var_1691)))
          absn2o_var_1701 = (kb_mn2o_var_203(indm_var_1690, ig_var_1691) + minorfrac_var_1684(jl_var_1712, lay_var_1692) * (kb_mn2o_var_203(indm_var_1690 + 1, ig_var_1691) - kb_mn2o_var_203(indm_var_1690, ig_var_1691)))
          taug_var_1661(jl_var_1712, 88 + ig_var_1691, lay_var_1692) = colo3_var_1672(jl_var_1712, lay_var_1692) * (fac00_var_1664(jl_var_1712, lay_var_1692) * absb_var_198(ind0_var_1686, ig_var_1691) + fac10_var_1666(jl_var_1712, lay_var_1692) * absb_var_198(ind0_var_1686 + 1, ig_var_1691) + fac01_var_1665(jl_var_1712, lay_var_1692) * absb_var_198(ind1_var_1687, ig_var_1691) + fac11_var_1667(jl_var_1712, lay_var_1692) * absb_var_198(ind1_var_1687 + 1, ig_var_1691)) + adjcolco2_var_1696 * absco2_var_1700 + coln2o_var_1673(jl_var_1712, lay_var_1692) * absn2o_var_1701 + wx_var_1662(jl_var_1712, 3, lay_var_1692) * cfc12_var_196(ig_var_1691) + wx_var_1662(jl_var_1712, 4, lay_var_1692) * cfc22adj(ig_var_1691)
          fracs_var_1680(jl_var_1712, 88 + ig_var_1691, lay_var_1692) = fracrefb_var_195(ig_var_1691)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol8
SUBROUTINE rrtm_taumol1(kidia_var_1713, kfdia_var_1714, klev_var_1715, taug_var_1717, pavel_var_1716, p_tauaerl_var_1718, fac00_var_1719, fac01_var_1720, fac10_var_1721, fac11_var_1722, forfac_var_1723, forfrac_var_1724, indfor_var_1735, jp_var_1725, jt_var_1726, jt1_var_1727, colh2o_var_1728, laytrop_var_1729, selffac_var_1730, selffrac_var_1731, indself_var_1733, fracs_var_1734, minorfrac_var_1732, indminor_var_1736, scaleminorn2, colbrd_var_1737)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta1, ONLY: absa_var_109, absb_var_110, forref_var_113, fracrefa_var_107, fracrefb_var_108, ka_mn2_var_111, kb_mn2, selfref_var_112
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1713
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1714
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1715
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1716(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1717(kidia_var_1713 : kfdia_var_1714, 140, klev_var_1715)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1718(kidia_var_1713 : kfdia_var_1714, klev_var_1715, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1719(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1720(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1721(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1722(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1723(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1724(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1725(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1726(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1727(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1728(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1729(kidia_var_1713 : kfdia_var_1714)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1730(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1731(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1732(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1733(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1734(kidia_var_1713 : kfdia_var_1714, 140, klev_var_1715)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1735(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1736(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(IN) :: scaleminorn2(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  REAL(KIND = 8), INTENT(IN) :: colbrd_var_1737(kidia_var_1713 : kfdia_var_1714, klev_var_1715)
  INTEGER(KIND = 4) :: ind0_var_1738, ind1_var_1739, inds_var_1740
  INTEGER(KIND = 4) :: indf_var_1741, indm_var_1742
  INTEGER(KIND = 4) :: ig_var_1743, lay_var_1744
  REAL(KIND = 8) :: taufor_var_1745, tauself_var_1746, corradj_var_1747, pp_var_1748, scalen2_var_1749, taun2_var_1750
  INTEGER(KIND = 4) :: laytrop_min_var_1751, laytrop_max_var_1752
  INTEGER(KIND = 4) :: ixc_var_1753(klev_var_1715), ixlow_var_1754(kfdia_var_1714, klev_var_1715), ixhigh_var_1755(kfdia_var_1714, klev_var_1715)
  INTEGER(KIND = 4) :: ich_var_1756, icl_var_1757, ixc0_var_1758, ixp_var_1759, jc_var_1760, jl_var_1761
  laytrop_min_var_1751 = MINVAL(laytrop_var_1729)
  laytrop_max_var_1752 = MAXVAL(laytrop_var_1729)
  ixlow_var_1754 = 0
  ixhigh_var_1755 = 0
  ixc_var_1753 = 0
  DO lay_var_1744 = laytrop_min_var_1751 + 1, laytrop_max_var_1752
    icl_var_1757 = 0
    ich_var_1756 = 0
    DO jc_var_1760 = kidia_var_1713, kfdia_var_1714
      IF (lay_var_1744 <= laytrop_var_1729(jc_var_1760)) THEN
        icl_var_1757 = icl_var_1757 + 1
        ixlow_var_1754(icl_var_1757, lay_var_1744) = jc_var_1760
      ELSE
        ich_var_1756 = ich_var_1756 + 1
        ixhigh_var_1755(ich_var_1756, lay_var_1744) = jc_var_1760
      END IF
    END DO
    ixc_var_1753(lay_var_1744) = icl_var_1757
  END DO
  DO lay_var_1744 = 1, laytrop_min_var_1751
    DO jl_var_1761 = kidia_var_1713, kfdia_var_1714
      ind0_var_1738 = ((jp_var_1725(jl_var_1761, lay_var_1744) - 1) * 5 + (jt_var_1726(jl_var_1761, lay_var_1744) - 1)) * nspa_var_216(1) + 1
      ind1_var_1739 = (jp_var_1725(jl_var_1761, lay_var_1744) * 5 + (jt1_var_1727(jl_var_1761, lay_var_1744) - 1)) * nspa_var_216(1) + 1
      inds_var_1740 = indself_var_1733(jl_var_1761, lay_var_1744)
      indf_var_1741 = indfor_var_1735(jl_var_1761, lay_var_1744)
      indm_var_1742 = indminor_var_1736(jl_var_1761, lay_var_1744)
      pp_var_1748 = pavel_var_1716(jl_var_1761, lay_var_1744)
      corradj_var_1747 = 1.0D0
      IF (pp_var_1748 .LT. 250.0D0) THEN
        corradj_var_1747 = 1.0D0 - 0.15D0 * (250.0D0 - pp_var_1748) / 154.4D0
      END IF
      scalen2_var_1749 = colbrd_var_1737(jl_var_1761, lay_var_1744) * scaleminorn2(jl_var_1761, lay_var_1744)
      DO ig_var_1743 = 1, 10
        tauself_var_1746 = selffac_var_1730(jl_var_1761, lay_var_1744) * (selfref_var_112(inds_var_1740, ig_var_1743) + selffrac_var_1731(jl_var_1761, lay_var_1744) * (selfref_var_112(inds_var_1740 + 1, ig_var_1743) - selfref_var_112(inds_var_1740, ig_var_1743)))
        taufor_var_1745 = forfac_var_1723(jl_var_1761, lay_var_1744) * (forref_var_113(indf_var_1741, ig_var_1743) + forfrac_var_1724(jl_var_1761, lay_var_1744) * (forref_var_113(indf_var_1741 + 1, ig_var_1743) - forref_var_113(indf_var_1741, ig_var_1743)))
        taun2_var_1750 = scalen2_var_1749 * (ka_mn2_var_111(indm_var_1742, ig_var_1743) + minorfrac_var_1732(jl_var_1761, lay_var_1744) * (ka_mn2_var_111(indm_var_1742 + 1, ig_var_1743) - ka_mn2_var_111(indm_var_1742, ig_var_1743)))
        taug_var_1717(jl_var_1761, ig_var_1743, lay_var_1744) = corradj_var_1747 * (colh2o_var_1728(jl_var_1761, lay_var_1744) * (fac00_var_1719(jl_var_1761, lay_var_1744) * absa_var_109(ind0_var_1738, ig_var_1743) + fac10_var_1721(jl_var_1761, lay_var_1744) * absa_var_109(ind0_var_1738 + 1, ig_var_1743) + fac01_var_1720(jl_var_1761, lay_var_1744) * absa_var_109(ind1_var_1739, ig_var_1743) + fac11_var_1722(jl_var_1761, lay_var_1744) * absa_var_109(ind1_var_1739 + 1, ig_var_1743)) + tauself_var_1746 + taufor_var_1745 + taun2_var_1750)
        fracs_var_1734(jl_var_1761, ig_var_1743, lay_var_1744) = fracrefa_var_107(ig_var_1743)
      END DO
    END DO
  END DO
  DO lay_var_1744 = laytrop_max_var_1752 + 1, klev_var_1715
    DO jl_var_1761 = kidia_var_1713, kfdia_var_1714
      ind0_var_1738 = ((jp_var_1725(jl_var_1761, lay_var_1744) - 13) * 5 + (jt_var_1726(jl_var_1761, lay_var_1744) - 1)) * nspb_var_217(1) + 1
      ind1_var_1739 = ((jp_var_1725(jl_var_1761, lay_var_1744) - 12) * 5 + (jt1_var_1727(jl_var_1761, lay_var_1744) - 1)) * nspb_var_217(1) + 1
      indf_var_1741 = indfor_var_1735(jl_var_1761, lay_var_1744)
      indm_var_1742 = indminor_var_1736(jl_var_1761, lay_var_1744)
      pp_var_1748 = pavel_var_1716(jl_var_1761, lay_var_1744)
      corradj_var_1747 = 1.0D0 - 0.15D0 * (pp_var_1748 / 95.6D0)
      scalen2_var_1749 = colbrd_var_1737(jl_var_1761, lay_var_1744) * scaleminorn2(jl_var_1761, lay_var_1744)
      DO ig_var_1743 = 1, 10
        taufor_var_1745 = forfac_var_1723(jl_var_1761, lay_var_1744) * (forref_var_113(indf_var_1741, ig_var_1743) + forfrac_var_1724(jl_var_1761, lay_var_1744) * (forref_var_113(indf_var_1741 + 1, ig_var_1743) - forref_var_113(indf_var_1741, ig_var_1743)))
        taun2_var_1750 = scalen2_var_1749 * (kb_mn2(indm_var_1742, ig_var_1743) + minorfrac_var_1732(jl_var_1761, lay_var_1744) * (kb_mn2(indm_var_1742 + 1, ig_var_1743) - kb_mn2(indm_var_1742, ig_var_1743)))
        taug_var_1717(jl_var_1761, ig_var_1743, lay_var_1744) = corradj_var_1747 * (colh2o_var_1728(jl_var_1761, lay_var_1744) * (fac00_var_1719(jl_var_1761, lay_var_1744) * absb_var_110(ind0_var_1738, ig_var_1743) + fac10_var_1721(jl_var_1761, lay_var_1744) * absb_var_110(ind0_var_1738 + 1, ig_var_1743) + fac01_var_1720(jl_var_1761, lay_var_1744) * absb_var_110(ind1_var_1739, ig_var_1743) + fac11_var_1722(jl_var_1761, lay_var_1744) * absb_var_110(ind1_var_1739 + 1, ig_var_1743)) + taufor_var_1745 + taun2_var_1750)
        fracs_var_1734(jl_var_1761, ig_var_1743, lay_var_1744) = fracrefb_var_108(ig_var_1743)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1752 /= laytrop_min_var_1751) THEN
    DO lay_var_1744 = laytrop_min_var_1751 + 1, laytrop_max_var_1752
      ixc0_var_1758 = ixc_var_1753(lay_var_1744)
      DO ixp_var_1759 = 1, ixc0_var_1758
        jl_var_1761 = ixlow_var_1754(ixp_var_1759, lay_var_1744)
        ind0_var_1738 = ((jp_var_1725(jl_var_1761, lay_var_1744) - 1) * 5 + (jt_var_1726(jl_var_1761, lay_var_1744) - 1)) * nspa_var_216(1) + 1
        ind1_var_1739 = (jp_var_1725(jl_var_1761, lay_var_1744) * 5 + (jt1_var_1727(jl_var_1761, lay_var_1744) - 1)) * nspa_var_216(1) + 1
        inds_var_1740 = indself_var_1733(jl_var_1761, lay_var_1744)
        indf_var_1741 = indfor_var_1735(jl_var_1761, lay_var_1744)
        indm_var_1742 = indminor_var_1736(jl_var_1761, lay_var_1744)
        pp_var_1748 = pavel_var_1716(jl_var_1761, lay_var_1744)
        corradj_var_1747 = 1.0D0
        IF (pp_var_1748 .LT. 250.0D0) THEN
          corradj_var_1747 = 1.0D0 - 0.15D0 * (250.0D0 - pp_var_1748) / 154.4D0
        END IF
        scalen2_var_1749 = colbrd_var_1737(jl_var_1761, lay_var_1744) * scaleminorn2(jl_var_1761, lay_var_1744)
        DO ig_var_1743 = 1, 10
          tauself_var_1746 = selffac_var_1730(jl_var_1761, lay_var_1744) * (selfref_var_112(inds_var_1740, ig_var_1743) + selffrac_var_1731(jl_var_1761, lay_var_1744) * (selfref_var_112(inds_var_1740 + 1, ig_var_1743) - selfref_var_112(inds_var_1740, ig_var_1743)))
          taufor_var_1745 = forfac_var_1723(jl_var_1761, lay_var_1744) * (forref_var_113(indf_var_1741, ig_var_1743) + forfrac_var_1724(jl_var_1761, lay_var_1744) * (forref_var_113(indf_var_1741 + 1, ig_var_1743) - forref_var_113(indf_var_1741, ig_var_1743)))
          taun2_var_1750 = scalen2_var_1749 * (ka_mn2_var_111(indm_var_1742, ig_var_1743) + minorfrac_var_1732(jl_var_1761, lay_var_1744) * (ka_mn2_var_111(indm_var_1742 + 1, ig_var_1743) - ka_mn2_var_111(indm_var_1742, ig_var_1743)))
          taug_var_1717(jl_var_1761, ig_var_1743, lay_var_1744) = corradj_var_1747 * (colh2o_var_1728(jl_var_1761, lay_var_1744) * (fac00_var_1719(jl_var_1761, lay_var_1744) * absa_var_109(ind0_var_1738, ig_var_1743) + fac10_var_1721(jl_var_1761, lay_var_1744) * absa_var_109(ind0_var_1738 + 1, ig_var_1743) + fac01_var_1720(jl_var_1761, lay_var_1744) * absa_var_109(ind1_var_1739, ig_var_1743) + fac11_var_1722(jl_var_1761, lay_var_1744) * absa_var_109(ind1_var_1739 + 1, ig_var_1743)) + tauself_var_1746 + taufor_var_1745 + taun2_var_1750)
          fracs_var_1734(jl_var_1761, ig_var_1743, lay_var_1744) = fracrefa_var_107(ig_var_1743)
        END DO
      END DO
      ixc0_var_1758 = kfdia_var_1714 - kidia_var_1713 + 1 - ixc0_var_1758
      DO ixp_var_1759 = 1, ixc0_var_1758
        jl_var_1761 = ixhigh_var_1755(ixp_var_1759, lay_var_1744)
        ind0_var_1738 = ((jp_var_1725(jl_var_1761, lay_var_1744) - 13) * 5 + (jt_var_1726(jl_var_1761, lay_var_1744) - 1)) * nspb_var_217(1) + 1
        ind1_var_1739 = ((jp_var_1725(jl_var_1761, lay_var_1744) - 12) * 5 + (jt1_var_1727(jl_var_1761, lay_var_1744) - 1)) * nspb_var_217(1) + 1
        indf_var_1741 = indfor_var_1735(jl_var_1761, lay_var_1744)
        indm_var_1742 = indminor_var_1736(jl_var_1761, lay_var_1744)
        pp_var_1748 = pavel_var_1716(jl_var_1761, lay_var_1744)
        corradj_var_1747 = 1.0D0 - 0.15D0 * (pp_var_1748 / 95.6D0)
        scalen2_var_1749 = colbrd_var_1737(jl_var_1761, lay_var_1744) * scaleminorn2(jl_var_1761, lay_var_1744)
        DO ig_var_1743 = 1, 10
          taufor_var_1745 = forfac_var_1723(jl_var_1761, lay_var_1744) * (forref_var_113(indf_var_1741, ig_var_1743) + forfrac_var_1724(jl_var_1761, lay_var_1744) * (forref_var_113(indf_var_1741 + 1, ig_var_1743) - forref_var_113(indf_var_1741, ig_var_1743)))
          taun2_var_1750 = scalen2_var_1749 * (kb_mn2(indm_var_1742, ig_var_1743) + minorfrac_var_1732(jl_var_1761, lay_var_1744) * (kb_mn2(indm_var_1742 + 1, ig_var_1743) - kb_mn2(indm_var_1742, ig_var_1743)))
          taug_var_1717(jl_var_1761, ig_var_1743, lay_var_1744) = corradj_var_1747 * (colh2o_var_1728(jl_var_1761, lay_var_1744) * (fac00_var_1719(jl_var_1761, lay_var_1744) * absb_var_110(ind0_var_1738, ig_var_1743) + fac10_var_1721(jl_var_1761, lay_var_1744) * absb_var_110(ind0_var_1738 + 1, ig_var_1743) + fac01_var_1720(jl_var_1761, lay_var_1744) * absb_var_110(ind1_var_1739, ig_var_1743) + fac11_var_1722(jl_var_1761, lay_var_1744) * absb_var_110(ind1_var_1739 + 1, ig_var_1743)) + taufor_var_1745 + taun2_var_1750)
          fracs_var_1734(jl_var_1761, ig_var_1743, lay_var_1744) = fracrefb_var_108(ig_var_1743)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol1
SUBROUTINE srtm_gas_optical_depth(kidia_var_1762, kfdia_var_1763, klev_var_1764, poneminus_var_1765, prmu0_var_1766, klaytrop_var_1767, pcolch4_var_1768, pcolco2_var_1769, pcolh2o_var_1770, pcolmol_var_1771, pcolo2_var_1772, pcolo3_var_1773, pforfac_var_1774, pforfrac_var_1775, kindfor_var_1776, pselffac_var_1777, pselffrac_var_1778, kindself_var_1779, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, pod_var_1787, pssa, pincsol)
  USE yoesrtwn, ONLY: ngc
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1762, kfdia_var_1763
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1764
  REAL(KIND = 8), INTENT(IN) :: poneminus_var_1765(kidia_var_1762 : kfdia_var_1763)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1766(kidia_var_1762 : kfdia_var_1763)
  INTEGER(KIND = 4), INTENT(IN) :: klaytrop_var_1767(kidia_var_1762 : kfdia_var_1763)
  REAL(KIND = 8), INTENT(IN) :: pcolch4_var_1768(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pcolco2_var_1769(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pcolh2o_var_1770(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pcolmol_var_1771(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pcolo2_var_1772(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pcolo3_var_1773(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pforfac_var_1774(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pforfrac_var_1775(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  INTEGER(KIND = 4), INTENT(IN) :: kindfor_var_1776(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pselffac_var_1777(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pselffrac_var_1778(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  INTEGER(KIND = 4), INTENT(IN) :: kindself_var_1779(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pfac00_var_1780(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pfac01_var_1781(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pfac10_var_1782(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(IN) :: pfac11_var_1783(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  INTEGER(KIND = 4), INTENT(IN) :: kjp_var_1784(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  INTEGER(KIND = 4), INTENT(IN) :: kjt_var_1785(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  INTEGER(KIND = 4), INTENT(IN) :: kjt1_var_1786(kidia_var_1762 : kfdia_var_1763, klev_var_1764)
  REAL(KIND = 8), INTENT(INOUT) :: pod_var_1787(kidia_var_1762 : kfdia_var_1763, klev_var_1764, 112)
  REAL(KIND = 8), INTENT(INOUT) :: pssa(kidia_var_1762 : kfdia_var_1763, klev_var_1764, 112)
  REAL(KIND = 8), INTENT(INOUT) :: pincsol(kidia_var_1762 : kfdia_var_1763, 112)
  INTEGER(KIND = 4) :: ib1, ib2, ibm, igt, iw(kidia_var_1762 : kfdia_var_1763), jb_var_1788, jg_var_1789, jk_var_1790, jl_var_1791, icount
  REAL(KIND = 8) :: ztaug(kidia_var_1762 : kfdia_var_1763, klev_var_1764, 16)
  REAL(KIND = 8) :: ztaur(kidia_var_1762 : kfdia_var_1763, klev_var_1764, 16)
  REAL(KIND = 8) :: zsflxzen(kidia_var_1762 : kfdia_var_1763, 16)
  ib1 = 16
  ib2 = 29
  icount = 0
  DO jl_var_1791 = kidia_var_1762, kfdia_var_1763
    IF (prmu0_var_1766(jl_var_1791) > 0.0D0) THEN
      icount = icount + 1
      iw(jl_var_1791) = 0
    END IF
  END DO
  IF (icount /= 0) THEN
    DO jb_var_1788 = ib1, ib2
      ibm = jb_var_1788 - 15
      igt = ngc(ibm)
      IF (jb_var_1788 == 16) THEN
        CALL srtm_taumol16(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, poneminus_var_1765, pcolh2o_var_1770, pcolch4_var_1768, pcolmol_var_1771, klaytrop_var_1767, pselffac_var_1777, pselffrac_var_1778, kindself_var_1779, pforfac_var_1774, pforfrac_var_1775, kindfor_var_1776, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 17) THEN
        CALL srtm_taumol17(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, poneminus_var_1765, pcolh2o_var_1770, pcolco2_var_1769, pcolmol_var_1771, klaytrop_var_1767, pselffac_var_1777, pselffrac_var_1778, kindself_var_1779, pforfac_var_1774, pforfrac_var_1775, kindfor_var_1776, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 18) THEN
        CALL srtm_taumol18(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, poneminus_var_1765, pcolh2o_var_1770, pcolch4_var_1768, pcolmol_var_1771, klaytrop_var_1767, pselffac_var_1777, pselffrac_var_1778, kindself_var_1779, pforfac_var_1774, pforfrac_var_1775, kindfor_var_1776, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 19) THEN
        CALL srtm_taumol19(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, poneminus_var_1765, pcolh2o_var_1770, pcolco2_var_1769, pcolmol_var_1771, klaytrop_var_1767, pselffac_var_1777, pselffrac_var_1778, kindself_var_1779, pforfac_var_1774, pforfrac_var_1775, kindfor_var_1776, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 20) THEN
        CALL srtm_taumol20(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, pcolh2o_var_1770, pcolch4_var_1768, pcolmol_var_1771, klaytrop_var_1767, pselffac_var_1777, pselffrac_var_1778, kindself_var_1779, pforfac_var_1774, pforfrac_var_1775, kindfor_var_1776, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 21) THEN
        CALL srtm_taumol21(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, poneminus_var_1765, pcolh2o_var_1770, pcolco2_var_1769, pcolmol_var_1771, klaytrop_var_1767, pselffac_var_1777, pselffrac_var_1778, kindself_var_1779, pforfac_var_1774, pforfrac_var_1775, kindfor_var_1776, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 22) THEN
        CALL srtm_taumol22(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, poneminus_var_1765, pcolh2o_var_1770, pcolmol_var_1771, pcolo2_var_1772, klaytrop_var_1767, pselffac_var_1777, pselffrac_var_1778, kindself_var_1779, pforfac_var_1774, pforfrac_var_1775, kindfor_var_1776, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 23) THEN
        CALL srtm_taumol23(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, pcolh2o_var_1770, pcolmol_var_1771, klaytrop_var_1767, pselffac_var_1777, pselffrac_var_1778, kindself_var_1779, pforfac_var_1774, pforfrac_var_1775, kindfor_var_1776, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 24) THEN
        CALL srtm_taumol24(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, poneminus_var_1765, pcolh2o_var_1770, pcolmol_var_1771, pcolo2_var_1772, pcolo3_var_1773, klaytrop_var_1767, pselffac_var_1777, pselffrac_var_1778, kindself_var_1779, pforfac_var_1774, pforfrac_var_1775, kindfor_var_1776, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 25) THEN
        CALL srtm_taumol25(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, pcolh2o_var_1770, pcolmol_var_1771, pcolo3_var_1773, klaytrop_var_1767, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 26) THEN
        CALL srtm_taumol26(kidia_var_1762, kfdia_var_1763, klev_var_1764, pcolmol_var_1771, klaytrop_var_1767, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 27) THEN
        CALL srtm_taumol27(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, pcolmol_var_1771, pcolo3_var_1773, klaytrop_var_1767, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 28) THEN
        CALL srtm_taumol28(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, poneminus_var_1765, pcolmol_var_1771, pcolo2_var_1772, pcolo3_var_1773, klaytrop_var_1767, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      ELSE IF (jb_var_1788 == 29) THEN
        CALL srtm_taumol29(kidia_var_1762, kfdia_var_1763, klev_var_1764, pfac00_var_1780, pfac01_var_1781, pfac10_var_1782, pfac11_var_1783, kjp_var_1784, kjt_var_1785, kjt1_var_1786, pcolh2o_var_1770, pcolco2_var_1769, pcolmol_var_1771, klaytrop_var_1767, pselffac_var_1777, pselffrac_var_1778, kindself_var_1779, pforfac_var_1774, pforfrac_var_1775, kindfor_var_1776, zsflxzen, ztaug, ztaur, prmu0_var_1766)
      END IF
      DO jg_var_1789 = 1, igt
        DO jl_var_1791 = kidia_var_1762, kfdia_var_1763
          IF (prmu0_var_1766(jl_var_1791) > 0.0D0) THEN
            iw(jl_var_1791) = iw(jl_var_1791) + 1
            pincsol(jl_var_1791, iw(jl_var_1791)) = zsflxzen(jl_var_1791, jg_var_1789)
          END IF
        END DO
        DO jk_var_1790 = 1, klev_var_1764
          DO jl_var_1791 = kidia_var_1762, kfdia_var_1763
            IF (prmu0_var_1766(jl_var_1791) > 0.0D0) THEN
              pod_var_1787(jl_var_1791, jk_var_1790, iw(jl_var_1791)) = ztaur(jl_var_1791, jk_var_1790, jg_var_1789) + ztaug(jl_var_1791, jk_var_1790, jg_var_1789)
              pssa(jl_var_1791, jk_var_1790, iw(jl_var_1791)) = ztaur(jl_var_1791, jk_var_1790, jg_var_1789) / pod_var_1787(jl_var_1791, jk_var_1790, iw(jl_var_1791))
            END IF
          END DO
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE srtm_gas_optical_depth
SUBROUTINE rrtm_setcoef_140gp(kidia_var_1792, kfdia_var_1793, klev_var_1794, p_coldry, p_wbroad, p_wkl, p_fac00_var_1795, p_fac01_var_1796, p_fac10_var_1797, p_fac11_var_1798, p_forfac_var_1799, p_forfrac_var_1800, k_indfor_var_1817, k_jp_var_1801, k_jt_var_1802, k_jt1_var_1803, p_colh2o_var_1804, p_colco2_var_1805, p_colo3_var_1806, p_coln2o, p_colch4_var_1807, p_colo2_var_1808, p_co2mult_var_1809, p_colbrd, k_laytrop_var_1810, k_layswtch_var_1811, k_laylow_var_1812, pavel_var_1813, p_tavel, p_selffac_var_1814, p_selffrac_var_1815, k_indself_var_1816, k_indminor, p_scaleminor, p_scaleminorn2, p_minorfrac, prat_h2oco2_var_1818, prat_h2oco2_1_var_1819, prat_h2oo3_var_1820, prat_h2oo3_1_var_1821, prat_h2on2o_var_1822, prat_h2on2o_1_var_1823, prat_h2och4_var_1824, prat_h2och4_1_var_1825, prat_n2oco2_var_1826, prat_n2oco2_1_var_1827, prat_o3co2_var_1828, prat_o3co2_1_var_1829)
  USE yoerrtrf, ONLY: chi_mls, preflog_var_214, tref_var_215
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1792
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1793
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1794
  REAL(KIND = 8), INTENT(IN) :: p_coldry(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(IN) :: p_wbroad(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_colbrd(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(IN) :: p_wkl(kidia_var_1792 : kfdia_var_1793, 35, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_fac00_var_1795(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_fac01_var_1796(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_fac10_var_1797(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_fac11_var_1798(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_forfac_var_1799(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_forfrac_var_1800(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jp_var_1801(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jt_var_1802(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jt1_var_1803(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_colh2o_var_1804(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_colco2_var_1805(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_colo3_var_1806(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_coln2o(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_colch4_var_1807(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_colo2_var_1808(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_co2mult_var_1809(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  INTEGER(KIND = 4), INTENT(OUT) :: k_laytrop_var_1810(kidia_var_1792 : kfdia_var_1793)
  INTEGER(KIND = 4), INTENT(OUT) :: k_layswtch_var_1811(kidia_var_1792 : kfdia_var_1793)
  INTEGER(KIND = 4), INTENT(OUT) :: k_laylow_var_1812(kidia_var_1792 : kfdia_var_1793)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1813(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(IN) :: p_tavel(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_selffac_var_1814(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_selffrac_var_1815(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indself_var_1816(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indfor_var_1817(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indminor(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_scaleminor(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_scaleminorn2(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: p_minorfrac(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  REAL(KIND = 8), INTENT(OUT) :: prat_h2oco2_var_1818(kidia_var_1792 : kfdia_var_1793, klev_var_1794), prat_h2oco2_1_var_1819(kidia_var_1792 : kfdia_var_1793, klev_var_1794), prat_h2oo3_var_1820(kidia_var_1792 : kfdia_var_1793, klev_var_1794), prat_h2oo3_1_var_1821(kidia_var_1792 : kfdia_var_1793, klev_var_1794), prat_h2on2o_var_1822(kidia_var_1792 : kfdia_var_1793, klev_var_1794), prat_h2on2o_1_var_1823(kidia_var_1792 : kfdia_var_1793, klev_var_1794), prat_h2och4_var_1824(kidia_var_1792 : kfdia_var_1793, klev_var_1794), prat_h2och4_1_var_1825(kidia_var_1792 : kfdia_var_1793, klev_var_1794), prat_n2oco2_var_1826(kidia_var_1792 : kfdia_var_1793, klev_var_1794), prat_n2oco2_1_var_1827(kidia_var_1792 : kfdia_var_1793, klev_var_1794), prat_o3co2_var_1828(kidia_var_1792 : kfdia_var_1793, klev_var_1794), prat_o3co2_1_var_1829(kidia_var_1792 : kfdia_var_1793, klev_var_1794)
  INTEGER(KIND = 4) :: jp1_var_1830, jlay
  INTEGER(KIND = 4) :: jlon_var_1831
  REAL(KIND = 8) :: z_co2reg_var_1832, z_compfp_var_1833, z_factor_var_1834, z_fp_var_1835, z_ft_var_1836, z_ft1_var_1837, z_plog_var_1838, z_scalefac_var_1839, z_stpfac_var_1840, z_water_var_1841
  DO jlon_var_1831 = kidia_var_1792, kfdia_var_1793
    z_stpfac_var_1840 = 0.29220138203356366D0
    k_laytrop_var_1810(jlon_var_1831) = 0
    k_layswtch_var_1811(jlon_var_1831) = 0
    k_laylow_var_1812(jlon_var_1831) = 0
    DO jlay = 1, klev_var_1794
      z_plog_var_1838 = LOG(pavel_var_1813(jlon_var_1831, jlay))
      k_jp_var_1801(jlon_var_1831, jlay) = INT(36.0D0 - 5 * (z_plog_var_1838 + 0.04D0))
      IF (k_jp_var_1801(jlon_var_1831, jlay) < 1) THEN
        k_jp_var_1801(jlon_var_1831, jlay) = 1
      ELSE IF (k_jp_var_1801(jlon_var_1831, jlay) > 58) THEN
        k_jp_var_1801(jlon_var_1831, jlay) = 58
      END IF
      jp1_var_1830 = k_jp_var_1801(jlon_var_1831, jlay) + 1
      z_fp_var_1835 = 5.0D0 * (preflog_var_214(k_jp_var_1801(jlon_var_1831, jlay)) - z_plog_var_1838)
      z_fp_var_1835 = MAX(- 1.0D0, MIN(1.0D0, z_fp_var_1835))
      k_jt_var_1802(jlon_var_1831, jlay) = INT(3.0D0 + (p_tavel(jlon_var_1831, jlay) - tref_var_215(k_jp_var_1801(jlon_var_1831, jlay))) / 15.0D0)
      IF (k_jt_var_1802(jlon_var_1831, jlay) < 1) THEN
        k_jt_var_1802(jlon_var_1831, jlay) = 1
      ELSE IF (k_jt_var_1802(jlon_var_1831, jlay) > 4) THEN
        k_jt_var_1802(jlon_var_1831, jlay) = 4
      END IF
      z_ft_var_1836 = ((p_tavel(jlon_var_1831, jlay) - tref_var_215(k_jp_var_1801(jlon_var_1831, jlay))) / 15.0D0) - REAL(k_jt_var_1802(jlon_var_1831, jlay) - 3)
      k_jt1_var_1803(jlon_var_1831, jlay) = INT(3.0D0 + (p_tavel(jlon_var_1831, jlay) - tref_var_215(jp1_var_1830)) / 15.0D0)
      IF (k_jt1_var_1803(jlon_var_1831, jlay) < 1) THEN
        k_jt1_var_1803(jlon_var_1831, jlay) = 1
      ELSE IF (k_jt1_var_1803(jlon_var_1831, jlay) > 4) THEN
        k_jt1_var_1803(jlon_var_1831, jlay) = 4
      END IF
      z_ft1_var_1837 = ((p_tavel(jlon_var_1831, jlay) - tref_var_215(jp1_var_1830)) / 15.0D0) - REAL(k_jt1_var_1803(jlon_var_1831, jlay) - 3)
      z_water_var_1841 = p_wkl(jlon_var_1831, 1, jlay) / p_coldry(jlon_var_1831, jlay)
      z_scalefac_var_1839 = pavel_var_1813(jlon_var_1831, jlay) * z_stpfac_var_1840 / p_tavel(jlon_var_1831, jlay)
      IF (z_plog_var_1838 > 4.56D0) THEN
        k_laytrop_var_1810(jlon_var_1831) = k_laytrop_var_1810(jlon_var_1831) + 1
        p_forfac_var_1799(jlon_var_1831, jlay) = z_scalefac_var_1839 / (1.0D0 + z_water_var_1841)
        z_factor_var_1834 = (332.0D0 - p_tavel(jlon_var_1831, jlay)) / 36.0D0
        k_indfor_var_1817(jlon_var_1831, jlay) = MIN(2, MAX(1, INT(z_factor_var_1834)))
        p_forfrac_var_1800(jlon_var_1831, jlay) = z_factor_var_1834 - REAL(k_indfor_var_1817(jlon_var_1831, jlay))
        p_selffac_var_1814(jlon_var_1831, jlay) = z_water_var_1841 * p_forfac_var_1799(jlon_var_1831, jlay)
        z_factor_var_1834 = (p_tavel(jlon_var_1831, jlay) - 188.0D0) / 7.2D0
        k_indself_var_1816(jlon_var_1831, jlay) = MIN(9, MAX(1, INT(z_factor_var_1834) - 7))
        p_selffrac_var_1815(jlon_var_1831, jlay) = z_factor_var_1834 - REAL(k_indself_var_1816(jlon_var_1831, jlay) + 7)
        p_scaleminor(jlon_var_1831, jlay) = pavel_var_1813(jlon_var_1831, jlay) / p_tavel(jlon_var_1831, jlay)
        p_scaleminorn2(jlon_var_1831, jlay) = (pavel_var_1813(jlon_var_1831, jlay) / p_tavel(jlon_var_1831, jlay)) * (p_wbroad(jlon_var_1831, jlay) / (p_coldry(jlon_var_1831, jlay) + p_wkl(jlon_var_1831, 1, jlay)))
        z_factor_var_1834 = (p_tavel(jlon_var_1831, jlay) - 180.8D0) / 7.2D0
        k_indminor(jlon_var_1831, jlay) = MIN(18, MAX(1, INT(z_factor_var_1834)))
        p_minorfrac(jlon_var_1831, jlay) = z_factor_var_1834 - REAL(k_indminor(jlon_var_1831, jlay))
        prat_h2oco2_var_1818(jlon_var_1831, jlay) = chi_mls(1, k_jp_var_1801(jlon_var_1831, jlay)) / chi_mls(2, k_jp_var_1801(jlon_var_1831, jlay))
        prat_h2oco2_1_var_1819(jlon_var_1831, jlay) = chi_mls(1, k_jp_var_1801(jlon_var_1831, jlay) + 1) / chi_mls(2, k_jp_var_1801(jlon_var_1831, jlay) + 1)
        prat_h2oo3_var_1820(jlon_var_1831, jlay) = chi_mls(1, k_jp_var_1801(jlon_var_1831, jlay)) / chi_mls(3, k_jp_var_1801(jlon_var_1831, jlay))
        prat_h2oo3_1_var_1821(jlon_var_1831, jlay) = chi_mls(1, k_jp_var_1801(jlon_var_1831, jlay) + 1) / chi_mls(3, k_jp_var_1801(jlon_var_1831, jlay) + 1)
        prat_h2on2o_var_1822(jlon_var_1831, jlay) = chi_mls(1, k_jp_var_1801(jlon_var_1831, jlay)) / chi_mls(4, k_jp_var_1801(jlon_var_1831, jlay))
        prat_h2on2o_1_var_1823(jlon_var_1831, jlay) = chi_mls(1, k_jp_var_1801(jlon_var_1831, jlay) + 1) / chi_mls(4, k_jp_var_1801(jlon_var_1831, jlay) + 1)
        prat_h2och4_var_1824(jlon_var_1831, jlay) = chi_mls(1, k_jp_var_1801(jlon_var_1831, jlay)) / chi_mls(6, k_jp_var_1801(jlon_var_1831, jlay))
        prat_h2och4_1_var_1825(jlon_var_1831, jlay) = chi_mls(1, k_jp_var_1801(jlon_var_1831, jlay) + 1) / chi_mls(6, k_jp_var_1801(jlon_var_1831, jlay) + 1)
        prat_n2oco2_var_1826(jlon_var_1831, jlay) = chi_mls(4, k_jp_var_1801(jlon_var_1831, jlay)) / chi_mls(2, k_jp_var_1801(jlon_var_1831, jlay))
        prat_n2oco2_1_var_1827(jlon_var_1831, jlay) = chi_mls(4, k_jp_var_1801(jlon_var_1831, jlay) + 1) / chi_mls(2, k_jp_var_1801(jlon_var_1831, jlay) + 1)
        p_colh2o_var_1804(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 1, jlay)
        p_colco2_var_1805(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 2, jlay)
        p_colo3_var_1806(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 3, jlay)
        p_coln2o(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 4, jlay)
        p_colch4_var_1807(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 6, jlay)
        p_colo2_var_1808(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 7, jlay)
        p_colbrd(jlon_var_1831, jlay) = 1D-20 * p_wbroad(jlon_var_1831, jlay)
        IF (p_colco2_var_1805(jlon_var_1831, jlay) == 0.0D0) p_colco2_var_1805(jlon_var_1831, jlay) = 1D-32 * p_coldry(jlon_var_1831, jlay)
        IF (p_coln2o(jlon_var_1831, jlay) == 0.0D0) p_coln2o(jlon_var_1831, jlay) = 1D-32 * p_coldry(jlon_var_1831, jlay)
        IF (p_colch4_var_1807(jlon_var_1831, jlay) == 0.0D0) p_colch4_var_1807(jlon_var_1831, jlay) = 1D-32 * p_coldry(jlon_var_1831, jlay)
        z_co2reg_var_1832 = 3.55D-24 * p_coldry(jlon_var_1831, jlay)
        p_co2mult_var_1809(jlon_var_1831, jlay) = (p_colco2_var_1805(jlon_var_1831, jlay) - z_co2reg_var_1832) * 272.63D0 * EXP(- 1919.4D0 / p_tavel(jlon_var_1831, jlay)) / (0.00087604D0 * p_tavel(jlon_var_1831, jlay))
      ELSE
        p_forfac_var_1799(jlon_var_1831, jlay) = z_scalefac_var_1839 / (1.0D0 + z_water_var_1841)
        z_factor_var_1834 = (p_tavel(jlon_var_1831, jlay) - 188.0D0) / 36.0D0
        k_indfor_var_1817(jlon_var_1831, jlay) = 3
        p_forfrac_var_1800(jlon_var_1831, jlay) = z_factor_var_1834 - 1.0D0
        p_selffac_var_1814(jlon_var_1831, jlay) = z_water_var_1841 * p_forfac_var_1799(jlon_var_1831, jlay)
        p_scaleminor(jlon_var_1831, jlay) = pavel_var_1813(jlon_var_1831, jlay) / p_tavel(jlon_var_1831, jlay)
        p_scaleminorn2(jlon_var_1831, jlay) = (pavel_var_1813(jlon_var_1831, jlay) / p_tavel(jlon_var_1831, jlay)) * (p_wbroad(jlon_var_1831, jlay) / (p_coldry(jlon_var_1831, jlay) + p_wkl(jlon_var_1831, 1, jlay)))
        z_factor_var_1834 = (p_tavel(jlon_var_1831, jlay) - 180.8D0) / 7.2D0
        k_indminor(jlon_var_1831, jlay) = MIN(18, MAX(1, INT(z_factor_var_1834)))
        p_minorfrac(jlon_var_1831, jlay) = z_factor_var_1834 - REAL(k_indminor(jlon_var_1831, jlay))
        prat_h2oco2_var_1818(jlon_var_1831, jlay) = chi_mls(1, k_jp_var_1801(jlon_var_1831, jlay)) / chi_mls(2, k_jp_var_1801(jlon_var_1831, jlay))
        prat_h2oco2_1_var_1819(jlon_var_1831, jlay) = chi_mls(1, k_jp_var_1801(jlon_var_1831, jlay) + 1) / chi_mls(2, k_jp_var_1801(jlon_var_1831, jlay) + 1)
        prat_o3co2_var_1828(jlon_var_1831, jlay) = chi_mls(3, k_jp_var_1801(jlon_var_1831, jlay)) / chi_mls(2, k_jp_var_1801(jlon_var_1831, jlay))
        prat_o3co2_1_var_1829(jlon_var_1831, jlay) = chi_mls(3, k_jp_var_1801(jlon_var_1831, jlay) + 1) / chi_mls(2, k_jp_var_1801(jlon_var_1831, jlay) + 1)
        p_colh2o_var_1804(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 1, jlay)
        p_colco2_var_1805(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 2, jlay)
        p_colo3_var_1806(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 3, jlay)
        p_coln2o(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 4, jlay)
        p_colch4_var_1807(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 6, jlay)
        p_colo2_var_1808(jlon_var_1831, jlay) = 1D-20 * p_wkl(jlon_var_1831, 7, jlay)
        p_colbrd(jlon_var_1831, jlay) = 1D-20 * p_wbroad(jlon_var_1831, jlay)
        IF (p_colco2_var_1805(jlon_var_1831, jlay) == 0.0D0) p_colco2_var_1805(jlon_var_1831, jlay) = 1D-32 * p_coldry(jlon_var_1831, jlay)
        IF (p_coln2o(jlon_var_1831, jlay) == 0.0D0) p_coln2o(jlon_var_1831, jlay) = 1D-32 * p_coldry(jlon_var_1831, jlay)
        IF (p_colch4_var_1807(jlon_var_1831, jlay) == 0.0D0) p_colch4_var_1807(jlon_var_1831, jlay) = 1D-32 * p_coldry(jlon_var_1831, jlay)
        z_co2reg_var_1832 = 3.55D-24 * p_coldry(jlon_var_1831, jlay)
        p_co2mult_var_1809(jlon_var_1831, jlay) = (p_colco2_var_1805(jlon_var_1831, jlay) - z_co2reg_var_1832) * 272.63D0 * EXP(- 1919.4D0 / p_tavel(jlon_var_1831, jlay)) / (0.00087604D0 * p_tavel(jlon_var_1831, jlay))
      END IF
      z_compfp_var_1833 = 1.0D0 - z_fp_var_1835
      p_fac10_var_1797(jlon_var_1831, jlay) = z_compfp_var_1833 * z_ft_var_1836
      p_fac00_var_1795(jlon_var_1831, jlay) = z_compfp_var_1833 * (1.0D0 - z_ft_var_1836)
      p_fac11_var_1798(jlon_var_1831, jlay) = z_fp_var_1835 * z_ft1_var_1837
      p_fac01_var_1796(jlon_var_1831, jlay) = z_fp_var_1835 * (1.0D0 - z_ft1_var_1837)
      p_selffac_var_1814(jlon_var_1831, jlay) = p_colh2o_var_1804(jlon_var_1831, jlay) * p_selffac_var_1814(jlon_var_1831, jlay)
      p_forfac_var_1799(jlon_var_1831, jlay) = p_colh2o_var_1804(jlon_var_1831, jlay) * p_forfac_var_1799(jlon_var_1831, jlay)
    END DO
    IF (k_laylow_var_1812(jlon_var_1831) == 0) k_laylow_var_1812(jlon_var_1831) = 1
  END DO
END SUBROUTINE rrtm_setcoef_140gp
SUBROUTINE rrtm_taumol12(kidia_var_1842, kfdia_var_1843, klev_var_1844, taug_var_1845, p_tauaerl_var_1846, fac00_var_1847, fac01_var_1848, fac10_var_1849, fac11_var_1850, forfac_var_1866, forfrac_var_1865, indfor_var_1864, jp_var_1851, jt_var_1852, jt1_var_1853, oneminus_var_1854, colh2o_var_1855, colco2_var_1856, laytrop_var_1857, selffac_var_1858, selffrac_var_1859, indself_var_1860, fracs_var_1861, rat_h2oco2_var_1862, rat_h2oco2_1_var_1863)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrtm, ONLY: ng12
  USE yoerrta12, ONLY: absa_var_127, forref_var_129, fracrefa_var_126, selfref_var_128
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1842
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1843
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1844
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1845(kidia_var_1842 : kfdia_var_1843, 140, klev_var_1844)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1846(kidia_var_1842 : kfdia_var_1843, klev_var_1844, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1847(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1848(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1849(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1850(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1851(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1852(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1853(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1854
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1855(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1856(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1857(kidia_var_1842 : kfdia_var_1843)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1858(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1859(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1860(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1861(kidia_var_1842 : kfdia_var_1843, 140, klev_var_1844)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1862(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1863(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1864(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1865(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1866(kidia_var_1842 : kfdia_var_1843, klev_var_1844)
  REAL(KIND = 8) :: speccomb_var_1867, speccomb1_var_1868, speccomb_planck_var_1869
  INTEGER(KIND = 4) :: ind0_var_1870, ind1_var_1871, inds_var_1872, indf_var_1873
  INTEGER(KIND = 4) :: ig_var_1874, js_var_1875, lay_var_1876, js1_var_1877, jpl_var_1878
  REAL(KIND = 8) :: fs_var_1879, specmult_var_1880, specparm_var_1881, fs1_var_1882, specmult1_var_1883, specparm1_var_1884, fpl_var_1885, specmult_planck_var_1886, specparm_planck_var_1887
  REAL(KIND = 8) :: fac000_var_1888, fac100_var_1889, fac200_var_1890, fac010_var_1891, fac110_var_1892, fac210_var_1893, fac001_var_1894, fac101_var_1895, fac201_var_1896, fac011_var_1897, fac111_var_1898, fac211_var_1899
  REAL(KIND = 8) :: p_var_1900, p4_var_1901, fk0_var_1902, fk1_var_1903, fk2_var_1904
  REAL(KIND = 8) :: taufor_var_1905, tauself_var_1906, tau_major_var_1907(8), tau_major1_var_1908(8)
  REAL(KIND = 8) :: refrat_planck_a_var_1909
  INTEGER(KIND = 4) :: laytrop_min_var_1910, laytrop_max_var_1911
  INTEGER(KIND = 4) :: ixc_var_1912(klev_var_1844), ixlow_var_1913(kfdia_var_1843, klev_var_1844), ixhigh_var_1914(kfdia_var_1843, klev_var_1844)
  INTEGER(KIND = 4) :: ich_var_1915, icl_var_1916, ixc0_var_1917, ixp_var_1918, jc_var_1919, jl_var_1920
  laytrop_min_var_1910 = MINVAL(laytrop_var_1857)
  laytrop_max_var_1911 = MAXVAL(laytrop_var_1857)
  ixlow_var_1913 = 0
  ixhigh_var_1914 = 0
  ixc_var_1912 = 0
  DO lay_var_1876 = laytrop_min_var_1910 + 1, laytrop_max_var_1911
    icl_var_1916 = 0
    ich_var_1915 = 0
    DO jc_var_1919 = kidia_var_1842, kfdia_var_1843
      IF (lay_var_1876 <= laytrop_var_1857(jc_var_1919)) THEN
        icl_var_1916 = icl_var_1916 + 1
        ixlow_var_1913(icl_var_1916, lay_var_1876) = jc_var_1919
      ELSE
        ich_var_1915 = ich_var_1915 + 1
        ixhigh_var_1914(ich_var_1915, lay_var_1876) = jc_var_1919
      END IF
    END DO
    ixc_var_1912(lay_var_1876) = icl_var_1916
  END DO
  refrat_planck_a_var_1909 = chi_mls(1, 10) / chi_mls(2, 10)
  DO lay_var_1876 = 1, laytrop_min_var_1910
    DO jl_var_1920 = kidia_var_1842, kfdia_var_1843
      speccomb_var_1867 = colh2o_var_1855(jl_var_1920, lay_var_1876) + rat_h2oco2_var_1862(jl_var_1920, lay_var_1876) * colco2_var_1856(jl_var_1920, lay_var_1876)
      specparm_var_1881 = MIN(colh2o_var_1855(jl_var_1920, lay_var_1876) / speccomb_var_1867, oneminus_var_1854)
      specmult_var_1880 = 8.0D0 * (specparm_var_1881)
      js_var_1875 = 1 + INT(specmult_var_1880)
      fs_var_1879 = ((specmult_var_1880) - AINT((specmult_var_1880)))
      speccomb1_var_1868 = colh2o_var_1855(jl_var_1920, lay_var_1876) + rat_h2oco2_1_var_1863(jl_var_1920, lay_var_1876) * colco2_var_1856(jl_var_1920, lay_var_1876)
      specparm1_var_1884 = MIN(colh2o_var_1855(jl_var_1920, lay_var_1876) / speccomb1_var_1868, oneminus_var_1854)
      specmult1_var_1883 = 8.0D0 * (specparm1_var_1884)
      js1_var_1877 = 1 + INT(specmult1_var_1883)
      fs1_var_1882 = ((specmult1_var_1883) - AINT((specmult1_var_1883)))
      speccomb_planck_var_1869 = colh2o_var_1855(jl_var_1920, lay_var_1876) + refrat_planck_a_var_1909 * colco2_var_1856(jl_var_1920, lay_var_1876)
      specparm_planck_var_1887 = MIN(colh2o_var_1855(jl_var_1920, lay_var_1876) / speccomb_planck_var_1869, oneminus_var_1854)
      specmult_planck_var_1886 = 8.0D0 * specparm_planck_var_1887
      jpl_var_1878 = 1 + INT(specmult_planck_var_1886)
      fpl_var_1885 = ((specmult_planck_var_1886) - AINT((specmult_planck_var_1886)))
      ind0_var_1870 = ((jp_var_1851(jl_var_1920, lay_var_1876) - 1) * 5 + (jt_var_1852(jl_var_1920, lay_var_1876) - 1)) * nspa_var_216(12) + js_var_1875
      ind1_var_1871 = (jp_var_1851(jl_var_1920, lay_var_1876) * 5 + (jt1_var_1853(jl_var_1920, lay_var_1876) - 1)) * nspa_var_216(12) + js1_var_1877
      inds_var_1872 = indself_var_1860(jl_var_1920, lay_var_1876)
      indf_var_1873 = indfor_var_1864(jl_var_1920, lay_var_1876)
      IF (specparm_var_1881 .LT. 0.125D0) THEN
        p_var_1900 = fs_var_1879 - 1.0D0
        p4_var_1901 = p_var_1900 ** 4
        fk0_var_1902 = p4_var_1901
        fk1_var_1903 = 1.0D0 - p_var_1900 - 2.0D0 * p4_var_1901
        fk2_var_1904 = p_var_1900 + p4_var_1901
        fac000_var_1888 = fk0_var_1902 * fac00_var_1847(jl_var_1920, lay_var_1876)
        fac100_var_1889 = fk1_var_1903 * fac00_var_1847(jl_var_1920, lay_var_1876)
        fac200_var_1890 = fk2_var_1904 * fac00_var_1847(jl_var_1920, lay_var_1876)
        fac010_var_1891 = fk0_var_1902 * fac10_var_1849(jl_var_1920, lay_var_1876)
        fac110_var_1892 = fk1_var_1903 * fac10_var_1849(jl_var_1920, lay_var_1876)
        fac210_var_1893 = fk2_var_1904 * fac10_var_1849(jl_var_1920, lay_var_1876)
      ELSE IF (specparm_var_1881 .GT. 0.875D0) THEN
        p_var_1900 = - fs_var_1879
        p4_var_1901 = p_var_1900 ** 4
        fk0_var_1902 = p4_var_1901
        fk1_var_1903 = 1.0D0 - p_var_1900 - 2.0D0 * p4_var_1901
        fk2_var_1904 = p_var_1900 + p4_var_1901
        fac000_var_1888 = fk0_var_1902 * fac00_var_1847(jl_var_1920, lay_var_1876)
        fac100_var_1889 = fk1_var_1903 * fac00_var_1847(jl_var_1920, lay_var_1876)
        fac200_var_1890 = fk2_var_1904 * fac00_var_1847(jl_var_1920, lay_var_1876)
        fac010_var_1891 = fk0_var_1902 * fac10_var_1849(jl_var_1920, lay_var_1876)
        fac110_var_1892 = fk1_var_1903 * fac10_var_1849(jl_var_1920, lay_var_1876)
        fac210_var_1893 = fk2_var_1904 * fac10_var_1849(jl_var_1920, lay_var_1876)
      ELSE
        fac000_var_1888 = (1.0D0 - fs_var_1879) * fac00_var_1847(jl_var_1920, lay_var_1876)
        fac010_var_1891 = (1.0D0 - fs_var_1879) * fac10_var_1849(jl_var_1920, lay_var_1876)
        fac100_var_1889 = fs_var_1879 * fac00_var_1847(jl_var_1920, lay_var_1876)
        fac110_var_1892 = fs_var_1879 * fac10_var_1849(jl_var_1920, lay_var_1876)
        fac200_var_1890 = 0.0D0
        fac210_var_1893 = 0.0D0
      END IF
      IF (specparm1_var_1884 .LT. 0.125D0) THEN
        p_var_1900 = fs1_var_1882 - 1.0D0
        p4_var_1901 = p_var_1900 ** 4
        fk0_var_1902 = p4_var_1901
        fk1_var_1903 = 1.0D0 - p_var_1900 - 2.0D0 * p4_var_1901
        fk2_var_1904 = p_var_1900 + p4_var_1901
        fac001_var_1894 = fk0_var_1902 * fac01_var_1848(jl_var_1920, lay_var_1876)
        fac101_var_1895 = fk1_var_1903 * fac01_var_1848(jl_var_1920, lay_var_1876)
        fac201_var_1896 = fk2_var_1904 * fac01_var_1848(jl_var_1920, lay_var_1876)
        fac011_var_1897 = fk0_var_1902 * fac11_var_1850(jl_var_1920, lay_var_1876)
        fac111_var_1898 = fk1_var_1903 * fac11_var_1850(jl_var_1920, lay_var_1876)
        fac211_var_1899 = fk2_var_1904 * fac11_var_1850(jl_var_1920, lay_var_1876)
      ELSE IF (specparm1_var_1884 .GT. 0.875D0) THEN
        p_var_1900 = - fs1_var_1882
        p4_var_1901 = p_var_1900 ** 4
        fk0_var_1902 = p4_var_1901
        fk1_var_1903 = 1.0D0 - p_var_1900 - 2.0D0 * p4_var_1901
        fk2_var_1904 = p_var_1900 + p4_var_1901
        fac001_var_1894 = fk0_var_1902 * fac01_var_1848(jl_var_1920, lay_var_1876)
        fac101_var_1895 = fk1_var_1903 * fac01_var_1848(jl_var_1920, lay_var_1876)
        fac201_var_1896 = fk2_var_1904 * fac01_var_1848(jl_var_1920, lay_var_1876)
        fac011_var_1897 = fk0_var_1902 * fac11_var_1850(jl_var_1920, lay_var_1876)
        fac111_var_1898 = fk1_var_1903 * fac11_var_1850(jl_var_1920, lay_var_1876)
        fac211_var_1899 = fk2_var_1904 * fac11_var_1850(jl_var_1920, lay_var_1876)
      ELSE
        fac001_var_1894 = (1.0D0 - fs1_var_1882) * fac01_var_1848(jl_var_1920, lay_var_1876)
        fac011_var_1897 = (1.0D0 - fs1_var_1882) * fac11_var_1850(jl_var_1920, lay_var_1876)
        fac101_var_1895 = fs1_var_1882 * fac01_var_1848(jl_var_1920, lay_var_1876)
        fac111_var_1898 = fs1_var_1882 * fac11_var_1850(jl_var_1920, lay_var_1876)
        fac201_var_1896 = 0.0D0
        fac211_var_1899 = 0.0D0
      END IF
      IF (specparm_var_1881 .LT. 0.125D0) THEN
        tau_major_var_1907(1 : ng12) = speccomb_var_1867 * (fac000_var_1888 * absa_var_127(ind0_var_1870, 1 : 8) + fac100_var_1889 * absa_var_127(ind0_var_1870 + 1, 1 : 8) + fac200_var_1890 * absa_var_127(ind0_var_1870 + 2, 1 : 8) + fac010_var_1891 * absa_var_127(ind0_var_1870 + 9, 1 : 8) + fac110_var_1892 * absa_var_127(ind0_var_1870 + 10, 1 : 8) + fac210_var_1893 * absa_var_127(ind0_var_1870 + 11, 1 : 8))
      ELSE IF (specparm_var_1881 .GT. 0.875D0) THEN
        tau_major_var_1907(1 : ng12) = speccomb_var_1867 * (fac200_var_1890 * absa_var_127(ind0_var_1870 - 1, 1 : 8) + fac100_var_1889 * absa_var_127(ind0_var_1870, 1 : 8) + fac000_var_1888 * absa_var_127(ind0_var_1870 + 1, 1 : 8) + fac210_var_1893 * absa_var_127(ind0_var_1870 + 8, 1 : 8) + fac110_var_1892 * absa_var_127(ind0_var_1870 + 9, 1 : 8) + fac010_var_1891 * absa_var_127(ind0_var_1870 + 10, 1 : 8))
      ELSE
        tau_major_var_1907(1 : ng12) = speccomb_var_1867 * (fac000_var_1888 * absa_var_127(ind0_var_1870, 1 : 8) + fac100_var_1889 * absa_var_127(ind0_var_1870 + 1, 1 : 8) + fac010_var_1891 * absa_var_127(ind0_var_1870 + 9, 1 : 8) + fac110_var_1892 * absa_var_127(ind0_var_1870 + 10, 1 : 8))
      END IF
      IF (specparm1_var_1884 .LT. 0.125D0) THEN
        tau_major1_var_1908(1 : ng12) = speccomb1_var_1868 * (fac001_var_1894 * absa_var_127(ind1_var_1871, 1 : 8) + fac101_var_1895 * absa_var_127(ind1_var_1871 + 1, 1 : 8) + fac201_var_1896 * absa_var_127(ind1_var_1871 + 2, 1 : 8) + fac011_var_1897 * absa_var_127(ind1_var_1871 + 9, 1 : 8) + fac111_var_1898 * absa_var_127(ind1_var_1871 + 10, 1 : 8) + fac211_var_1899 * absa_var_127(ind1_var_1871 + 11, 1 : 8))
      ELSE IF (specparm1_var_1884 .GT. 0.875D0) THEN
        tau_major1_var_1908(1 : ng12) = speccomb1_var_1868 * (fac201_var_1896 * absa_var_127(ind1_var_1871 - 1, 1 : 8) + fac101_var_1895 * absa_var_127(ind1_var_1871, 1 : 8) + fac001_var_1894 * absa_var_127(ind1_var_1871 + 1, 1 : 8) + fac211_var_1899 * absa_var_127(ind1_var_1871 + 8, 1 : 8) + fac111_var_1898 * absa_var_127(ind1_var_1871 + 9, 1 : 8) + fac011_var_1897 * absa_var_127(ind1_var_1871 + 10, 1 : 8))
      ELSE
        tau_major1_var_1908(1 : ng12) = speccomb1_var_1868 * (fac001_var_1894 * absa_var_127(ind1_var_1871, 1 : 8) + fac101_var_1895 * absa_var_127(ind1_var_1871 + 1, 1 : 8) + fac011_var_1897 * absa_var_127(ind1_var_1871 + 9, 1 : 8) + fac111_var_1898 * absa_var_127(ind1_var_1871 + 10, 1 : 8))
      END IF
      DO ig_var_1874 = 1, 8
        tauself_var_1906 = selffac_var_1858(jl_var_1920, lay_var_1876) * (selfref_var_128(inds_var_1872, ig_var_1874) + selffrac_var_1859(jl_var_1920, lay_var_1876) * (selfref_var_128(inds_var_1872 + 1, ig_var_1874) - selfref_var_128(inds_var_1872, ig_var_1874)))
        taufor_var_1905 = forfac_var_1866(jl_var_1920, lay_var_1876) * (forref_var_129(indf_var_1873, ig_var_1874) + forfrac_var_1865(jl_var_1920, lay_var_1876) * (forref_var_129(indf_var_1873 + 1, ig_var_1874) - forref_var_129(indf_var_1873, ig_var_1874)))
        taug_var_1845(jl_var_1920, 122 + ig_var_1874, lay_var_1876) = tau_major_var_1907(ig_var_1874) + tau_major1_var_1908(ig_var_1874) + tauself_var_1906 + taufor_var_1905
        fracs_var_1861(jl_var_1920, 122 + ig_var_1874, lay_var_1876) = fracrefa_var_126(ig_var_1874, jpl_var_1878) + fpl_var_1885 * (fracrefa_var_126(ig_var_1874, jpl_var_1878 + 1) - fracrefa_var_126(ig_var_1874, jpl_var_1878))
      END DO
    END DO
  END DO
  DO ig_var_1874 = 1, 8
    DO lay_var_1876 = laytrop_max_var_1911 + 1, klev_var_1844
      DO jl_var_1920 = kidia_var_1842, kfdia_var_1843
        taug_var_1845(jl_var_1920, 122 + ig_var_1874, lay_var_1876) = 0.0D0
        fracs_var_1861(jl_var_1920, 122 + ig_var_1874, lay_var_1876) = 0.0D0
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1911 /= laytrop_min_var_1910) THEN
    DO lay_var_1876 = laytrop_min_var_1910 + 1, laytrop_max_var_1911
      ixc0_var_1917 = ixc_var_1912(lay_var_1876)
      DO ixp_var_1918 = 1, ixc0_var_1917
        jl_var_1920 = ixlow_var_1913(ixp_var_1918, lay_var_1876)
        speccomb_var_1867 = colh2o_var_1855(jl_var_1920, lay_var_1876) + rat_h2oco2_var_1862(jl_var_1920, lay_var_1876) * colco2_var_1856(jl_var_1920, lay_var_1876)
        specparm_var_1881 = MIN(colh2o_var_1855(jl_var_1920, lay_var_1876) / speccomb_var_1867, oneminus_var_1854)
        specmult_var_1880 = 8.0D0 * (specparm_var_1881)
        js_var_1875 = 1 + INT(specmult_var_1880)
        fs_var_1879 = ((specmult_var_1880) - AINT((specmult_var_1880)))
        speccomb1_var_1868 = colh2o_var_1855(jl_var_1920, lay_var_1876) + rat_h2oco2_1_var_1863(jl_var_1920, lay_var_1876) * colco2_var_1856(jl_var_1920, lay_var_1876)
        specparm1_var_1884 = MIN(colh2o_var_1855(jl_var_1920, lay_var_1876) / speccomb1_var_1868, oneminus_var_1854)
        specmult1_var_1883 = 8.0D0 * (specparm1_var_1884)
        js1_var_1877 = 1 + INT(specmult1_var_1883)
        fs1_var_1882 = ((specmult1_var_1883) - AINT((specmult1_var_1883)))
        speccomb_planck_var_1869 = colh2o_var_1855(jl_var_1920, lay_var_1876) + refrat_planck_a_var_1909 * colco2_var_1856(jl_var_1920, lay_var_1876)
        specparm_planck_var_1887 = MIN(colh2o_var_1855(jl_var_1920, lay_var_1876) / speccomb_planck_var_1869, oneminus_var_1854)
        specmult_planck_var_1886 = 8.0D0 * specparm_planck_var_1887
        jpl_var_1878 = 1 + INT(specmult_planck_var_1886)
        fpl_var_1885 = ((specmult_planck_var_1886) - AINT((specmult_planck_var_1886)))
        ind0_var_1870 = ((jp_var_1851(jl_var_1920, lay_var_1876) - 1) * 5 + (jt_var_1852(jl_var_1920, lay_var_1876) - 1)) * nspa_var_216(12) + js_var_1875
        ind1_var_1871 = (jp_var_1851(jl_var_1920, lay_var_1876) * 5 + (jt1_var_1853(jl_var_1920, lay_var_1876) - 1)) * nspa_var_216(12) + js1_var_1877
        inds_var_1872 = indself_var_1860(jl_var_1920, lay_var_1876)
        indf_var_1873 = indfor_var_1864(jl_var_1920, lay_var_1876)
        IF (specparm_var_1881 .LT. 0.125D0) THEN
          p_var_1900 = fs_var_1879 - 1.0D0
          p4_var_1901 = p_var_1900 ** 4
          fk0_var_1902 = p4_var_1901
          fk1_var_1903 = 1.0D0 - p_var_1900 - 2.0D0 * p4_var_1901
          fk2_var_1904 = p_var_1900 + p4_var_1901
          fac000_var_1888 = fk0_var_1902 * fac00_var_1847(jl_var_1920, lay_var_1876)
          fac100_var_1889 = fk1_var_1903 * fac00_var_1847(jl_var_1920, lay_var_1876)
          fac200_var_1890 = fk2_var_1904 * fac00_var_1847(jl_var_1920, lay_var_1876)
          fac010_var_1891 = fk0_var_1902 * fac10_var_1849(jl_var_1920, lay_var_1876)
          fac110_var_1892 = fk1_var_1903 * fac10_var_1849(jl_var_1920, lay_var_1876)
          fac210_var_1893 = fk2_var_1904 * fac10_var_1849(jl_var_1920, lay_var_1876)
        ELSE IF (specparm_var_1881 .GT. 0.875D0) THEN
          p_var_1900 = - fs_var_1879
          p4_var_1901 = p_var_1900 ** 4
          fk0_var_1902 = p4_var_1901
          fk1_var_1903 = 1.0D0 - p_var_1900 - 2.0D0 * p4_var_1901
          fk2_var_1904 = p_var_1900 + p4_var_1901
          fac000_var_1888 = fk0_var_1902 * fac00_var_1847(jl_var_1920, lay_var_1876)
          fac100_var_1889 = fk1_var_1903 * fac00_var_1847(jl_var_1920, lay_var_1876)
          fac200_var_1890 = fk2_var_1904 * fac00_var_1847(jl_var_1920, lay_var_1876)
          fac010_var_1891 = fk0_var_1902 * fac10_var_1849(jl_var_1920, lay_var_1876)
          fac110_var_1892 = fk1_var_1903 * fac10_var_1849(jl_var_1920, lay_var_1876)
          fac210_var_1893 = fk2_var_1904 * fac10_var_1849(jl_var_1920, lay_var_1876)
        ELSE
          fac000_var_1888 = (1.0D0 - fs_var_1879) * fac00_var_1847(jl_var_1920, lay_var_1876)
          fac010_var_1891 = (1.0D0 - fs_var_1879) * fac10_var_1849(jl_var_1920, lay_var_1876)
          fac100_var_1889 = fs_var_1879 * fac00_var_1847(jl_var_1920, lay_var_1876)
          fac110_var_1892 = fs_var_1879 * fac10_var_1849(jl_var_1920, lay_var_1876)
          fac200_var_1890 = 0.0D0
          fac210_var_1893 = 0.0D0
        END IF
        IF (specparm1_var_1884 .LT. 0.125D0) THEN
          p_var_1900 = fs1_var_1882 - 1.0D0
          p4_var_1901 = p_var_1900 ** 4
          fk0_var_1902 = p4_var_1901
          fk1_var_1903 = 1.0D0 - p_var_1900 - 2.0D0 * p4_var_1901
          fk2_var_1904 = p_var_1900 + p4_var_1901
          fac001_var_1894 = fk0_var_1902 * fac01_var_1848(jl_var_1920, lay_var_1876)
          fac101_var_1895 = fk1_var_1903 * fac01_var_1848(jl_var_1920, lay_var_1876)
          fac201_var_1896 = fk2_var_1904 * fac01_var_1848(jl_var_1920, lay_var_1876)
          fac011_var_1897 = fk0_var_1902 * fac11_var_1850(jl_var_1920, lay_var_1876)
          fac111_var_1898 = fk1_var_1903 * fac11_var_1850(jl_var_1920, lay_var_1876)
          fac211_var_1899 = fk2_var_1904 * fac11_var_1850(jl_var_1920, lay_var_1876)
        ELSE IF (specparm1_var_1884 .GT. 0.875D0) THEN
          p_var_1900 = - fs1_var_1882
          p4_var_1901 = p_var_1900 ** 4
          fk0_var_1902 = p4_var_1901
          fk1_var_1903 = 1.0D0 - p_var_1900 - 2.0D0 * p4_var_1901
          fk2_var_1904 = p_var_1900 + p4_var_1901
          fac001_var_1894 = fk0_var_1902 * fac01_var_1848(jl_var_1920, lay_var_1876)
          fac101_var_1895 = fk1_var_1903 * fac01_var_1848(jl_var_1920, lay_var_1876)
          fac201_var_1896 = fk2_var_1904 * fac01_var_1848(jl_var_1920, lay_var_1876)
          fac011_var_1897 = fk0_var_1902 * fac11_var_1850(jl_var_1920, lay_var_1876)
          fac111_var_1898 = fk1_var_1903 * fac11_var_1850(jl_var_1920, lay_var_1876)
          fac211_var_1899 = fk2_var_1904 * fac11_var_1850(jl_var_1920, lay_var_1876)
        ELSE
          fac001_var_1894 = (1.0D0 - fs1_var_1882) * fac01_var_1848(jl_var_1920, lay_var_1876)
          fac011_var_1897 = (1.0D0 - fs1_var_1882) * fac11_var_1850(jl_var_1920, lay_var_1876)
          fac101_var_1895 = fs1_var_1882 * fac01_var_1848(jl_var_1920, lay_var_1876)
          fac111_var_1898 = fs1_var_1882 * fac11_var_1850(jl_var_1920, lay_var_1876)
          fac201_var_1896 = 0.0D0
          fac211_var_1899 = 0.0D0
        END IF
        IF (specparm_var_1881 .LT. 0.125D0) THEN
          tau_major_var_1907(1 : ng12) = speccomb_var_1867 * (fac000_var_1888 * absa_var_127(ind0_var_1870, 1 : 8) + fac100_var_1889 * absa_var_127(ind0_var_1870 + 1, 1 : 8) + fac200_var_1890 * absa_var_127(ind0_var_1870 + 2, 1 : 8) + fac010_var_1891 * absa_var_127(ind0_var_1870 + 9, 1 : 8) + fac110_var_1892 * absa_var_127(ind0_var_1870 + 10, 1 : 8) + fac210_var_1893 * absa_var_127(ind0_var_1870 + 11, 1 : 8))
        ELSE IF (specparm_var_1881 .GT. 0.875D0) THEN
          tau_major_var_1907(1 : ng12) = speccomb_var_1867 * (fac200_var_1890 * absa_var_127(ind0_var_1870 - 1, 1 : 8) + fac100_var_1889 * absa_var_127(ind0_var_1870, 1 : 8) + fac000_var_1888 * absa_var_127(ind0_var_1870 + 1, 1 : 8) + fac210_var_1893 * absa_var_127(ind0_var_1870 + 8, 1 : 8) + fac110_var_1892 * absa_var_127(ind0_var_1870 + 9, 1 : 8) + fac010_var_1891 * absa_var_127(ind0_var_1870 + 10, 1 : 8))
        ELSE
          tau_major_var_1907(1 : ng12) = speccomb_var_1867 * (fac000_var_1888 * absa_var_127(ind0_var_1870, 1 : 8) + fac100_var_1889 * absa_var_127(ind0_var_1870 + 1, 1 : 8) + fac010_var_1891 * absa_var_127(ind0_var_1870 + 9, 1 : 8) + fac110_var_1892 * absa_var_127(ind0_var_1870 + 10, 1 : 8))
        END IF
        IF (specparm1_var_1884 .LT. 0.125D0) THEN
          tau_major1_var_1908(1 : ng12) = speccomb1_var_1868 * (fac001_var_1894 * absa_var_127(ind1_var_1871, 1 : 8) + fac101_var_1895 * absa_var_127(ind1_var_1871 + 1, 1 : 8) + fac201_var_1896 * absa_var_127(ind1_var_1871 + 2, 1 : 8) + fac011_var_1897 * absa_var_127(ind1_var_1871 + 9, 1 : 8) + fac111_var_1898 * absa_var_127(ind1_var_1871 + 10, 1 : 8) + fac211_var_1899 * absa_var_127(ind1_var_1871 + 11, 1 : 8))
        ELSE IF (specparm1_var_1884 .GT. 0.875D0) THEN
          tau_major1_var_1908(1 : ng12) = speccomb1_var_1868 * (fac201_var_1896 * absa_var_127(ind1_var_1871 - 1, 1 : 8) + fac101_var_1895 * absa_var_127(ind1_var_1871, 1 : 8) + fac001_var_1894 * absa_var_127(ind1_var_1871 + 1, 1 : 8) + fac211_var_1899 * absa_var_127(ind1_var_1871 + 8, 1 : 8) + fac111_var_1898 * absa_var_127(ind1_var_1871 + 9, 1 : 8) + fac011_var_1897 * absa_var_127(ind1_var_1871 + 10, 1 : 8))
        ELSE
          tau_major1_var_1908(1 : ng12) = speccomb1_var_1868 * (fac001_var_1894 * absa_var_127(ind1_var_1871, 1 : 8) + fac101_var_1895 * absa_var_127(ind1_var_1871 + 1, 1 : 8) + fac011_var_1897 * absa_var_127(ind1_var_1871 + 9, 1 : 8) + fac111_var_1898 * absa_var_127(ind1_var_1871 + 10, 1 : 8))
        END IF
        DO ig_var_1874 = 1, 8
          tauself_var_1906 = selffac_var_1858(jl_var_1920, lay_var_1876) * (selfref_var_128(inds_var_1872, ig_var_1874) + selffrac_var_1859(jl_var_1920, lay_var_1876) * (selfref_var_128(inds_var_1872 + 1, ig_var_1874) - selfref_var_128(inds_var_1872, ig_var_1874)))
          taufor_var_1905 = forfac_var_1866(jl_var_1920, lay_var_1876) * (forref_var_129(indf_var_1873, ig_var_1874) + forfrac_var_1865(jl_var_1920, lay_var_1876) * (forref_var_129(indf_var_1873 + 1, ig_var_1874) - forref_var_129(indf_var_1873, ig_var_1874)))
          taug_var_1845(jl_var_1920, 122 + ig_var_1874, lay_var_1876) = tau_major_var_1907(ig_var_1874) + tau_major1_var_1908(ig_var_1874) + tauself_var_1906 + taufor_var_1905
          fracs_var_1861(jl_var_1920, 122 + ig_var_1874, lay_var_1876) = fracrefa_var_126(ig_var_1874, jpl_var_1878) + fpl_var_1885 * (fracrefa_var_126(ig_var_1874, jpl_var_1878 + 1) - fracrefa_var_126(ig_var_1874, jpl_var_1878))
        END DO
      END DO
      ixc0_var_1917 = kfdia_var_1843 - kidia_var_1842 + 1 - ixc0_var_1917
      DO ig_var_1874 = 1, 8
        DO ixp_var_1918 = 1, ixc0_var_1917
          jl_var_1920 = ixhigh_var_1914(ixp_var_1918, lay_var_1876)
          taug_var_1845(jl_var_1920, 122 + ig_var_1874, lay_var_1876) = 0.0D0
          fracs_var_1861(jl_var_1920, 122 + ig_var_1874, lay_var_1876) = 0.0D0
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol12
SUBROUTINE srtm_taumol26(kidia_var_1921, kfdia_var_1922, klev_var_1923, p_colmol_var_1924, k_laytrop_var_1925, p_sfluxzen_var_1926, p_taug_var_1927, p_taur_var_1928, prmu0_var_1929)
  USE yoesrta26, ONLY: raylc_var_294, sfluxrefc_var_293
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1921, kfdia_var_1922
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1923
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1924(kidia_var_1921 : kfdia_var_1922, klev_var_1923)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1925(kidia_var_1921 : kfdia_var_1922)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1926(kidia_var_1921 : kfdia_var_1922, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1927(kidia_var_1921 : kfdia_var_1922, klev_var_1923, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1928(kidia_var_1921 : kfdia_var_1922, klev_var_1923, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1929(kidia_var_1921 : kfdia_var_1922)
  INTEGER(KIND = 4) :: ig_var_1930, i_lay_var_1931, i_laysolfr_var_1932(kidia_var_1921 : kfdia_var_1922), i_nlayers_var_1933, iplon_var_1934
  INTEGER(KIND = 4) :: laytrop_min_var_1935, laytrop_max_var_1936
  laytrop_min_var_1935 = MINVAL(k_laytrop_var_1925(kidia_var_1921 : kfdia_var_1922))
  laytrop_max_var_1936 = MAXVAL(k_laytrop_var_1925(kidia_var_1921 : kfdia_var_1922))
  i_nlayers_var_1933 = klev_var_1923
  DO iplon_var_1934 = kidia_var_1921, kfdia_var_1922
    i_laysolfr_var_1932(iplon_var_1934) = k_laytrop_var_1925(iplon_var_1934)
  END DO
  DO i_lay_var_1931 = 1, laytrop_min_var_1935
    DO iplon_var_1934 = kidia_var_1921, kfdia_var_1922
      DO ig_var_1930 = 1, 6
        IF (i_lay_var_1931 == i_laysolfr_var_1932(iplon_var_1934)) p_sfluxzen_var_1926(iplon_var_1934, ig_var_1930) = sfluxrefc_var_293(ig_var_1930)
        p_taug_var_1927(iplon_var_1934, i_lay_var_1931, ig_var_1930) = 0.0D0
        p_taur_var_1928(iplon_var_1934, i_lay_var_1931, ig_var_1930) = p_colmol_var_1924(iplon_var_1934, i_lay_var_1931) * raylc_var_294(ig_var_1930)
      END DO
    END DO
  END DO
  DO i_lay_var_1931 = laytrop_min_var_1935 + 1, laytrop_max_var_1936
    DO iplon_var_1934 = kidia_var_1921, kfdia_var_1922
      IF (i_lay_var_1931 <= k_laytrop_var_1925(iplon_var_1934)) THEN
        DO ig_var_1930 = 1, 6
          IF (i_lay_var_1931 == i_laysolfr_var_1932(iplon_var_1934)) p_sfluxzen_var_1926(iplon_var_1934, ig_var_1930) = sfluxrefc_var_293(ig_var_1930)
          p_taug_var_1927(iplon_var_1934, i_lay_var_1931, ig_var_1930) = 0.0D0
          p_taur_var_1928(iplon_var_1934, i_lay_var_1931, ig_var_1930) = p_colmol_var_1924(iplon_var_1934, i_lay_var_1931) * raylc_var_294(ig_var_1930)
        END DO
      ELSE
        DO ig_var_1930 = 1, 6
          p_taug_var_1927(iplon_var_1934, i_lay_var_1931, ig_var_1930) = 0.0D0
          p_taur_var_1928(iplon_var_1934, i_lay_var_1931, ig_var_1930) = p_colmol_var_1924(iplon_var_1934, i_lay_var_1931) * raylc_var_294(ig_var_1930)
        END DO
      END IF
    END DO
  END DO
  DO ig_var_1930 = 1, 6
    DO i_lay_var_1931 = laytrop_max_var_1936 + 1, i_nlayers_var_1933
      DO iplon_var_1934 = kidia_var_1921, kfdia_var_1922
        p_taug_var_1927(iplon_var_1934, i_lay_var_1931, ig_var_1930) = 0.0D0
        p_taur_var_1928(iplon_var_1934, i_lay_var_1931, ig_var_1930) = p_colmol_var_1924(iplon_var_1934, i_lay_var_1931) * raylc_var_294(ig_var_1930)
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol26
SUBROUTINE rrtm_taumol2(kidia_var_1937, kfdia_var_1938, klev_var_1939, taug_var_1940, pavel_var_1941, coldry_var_1942, p_tauaerl_var_1943, fac00_var_1944, fac01_var_1945, fac10_var_1946, fac11_var_1947, forfac_var_1949, forfrac_var_1948, indfor_var_1958, jp_var_1950, jt_var_1951, jt1_var_1952, colh2o_var_1953, laytrop_var_1954, selffac_var_1955, selffrac_var_1956, indself_var_1957, fracs_var_1959)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta2, ONLY: absa_var_155, absb_var_156, forref_var_158, fracrefa_var_153, fracrefb_var_154, selfref_var_157
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1937
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1938
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1939
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1940(kidia_var_1937 : kfdia_var_1938, 140, klev_var_1939)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1941(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1942(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1943(kidia_var_1937 : kfdia_var_1938, klev_var_1939, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1944(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1945(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1946(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1947(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1948(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1949(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1950(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1951(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1952(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1953(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1954(kidia_var_1937 : kfdia_var_1938)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1955(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1956(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1957(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1958(kidia_var_1937 : kfdia_var_1938, klev_var_1939)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1959(kidia_var_1937 : kfdia_var_1938, 140, klev_var_1939)
  INTEGER(KIND = 4) :: ind0_var_1960, ind1_var_1961, inds_var_1962, indf_var_1963
  INTEGER(KIND = 4) :: ig_var_1964, lay_var_1965
  REAL(KIND = 8) :: taufor_var_1966, tauself_var_1967, corradj_var_1968, pp_var_1969
  INTEGER(KIND = 4) :: laytrop_min_var_1970, laytrop_max_var_1971
  INTEGER(KIND = 4) :: ixc_var_1972(klev_var_1939), ixlow_var_1973(kfdia_var_1938, klev_var_1939), ixhigh_var_1974(kfdia_var_1938, klev_var_1939)
  INTEGER(KIND = 4) :: ich_var_1975, icl_var_1976, ixc0_var_1977, ixp_var_1978, jc_var_1979, jl_var_1980
  laytrop_min_var_1970 = MINVAL(laytrop_var_1954)
  laytrop_max_var_1971 = MAXVAL(laytrop_var_1954)
  ixlow_var_1973 = 0
  ixhigh_var_1974 = 0
  ixc_var_1972 = 0
  DO lay_var_1965 = laytrop_min_var_1970 + 1, laytrop_max_var_1971
    icl_var_1976 = 0
    ich_var_1975 = 0
    DO jc_var_1979 = kidia_var_1937, kfdia_var_1938
      IF (lay_var_1965 <= laytrop_var_1954(jc_var_1979)) THEN
        icl_var_1976 = icl_var_1976 + 1
        ixlow_var_1973(icl_var_1976, lay_var_1965) = jc_var_1979
      ELSE
        ich_var_1975 = ich_var_1975 + 1
        ixhigh_var_1974(ich_var_1975, lay_var_1965) = jc_var_1979
      END IF
    END DO
    ixc_var_1972(lay_var_1965) = icl_var_1976
  END DO
  DO lay_var_1965 = 1, laytrop_min_var_1970
    DO jl_var_1980 = kidia_var_1937, kfdia_var_1938
      ind0_var_1960 = ((jp_var_1950(jl_var_1980, lay_var_1965) - 1) * 5 + (jt_var_1951(jl_var_1980, lay_var_1965) - 1)) * nspa_var_216(2) + 1
      ind1_var_1961 = (jp_var_1950(jl_var_1980, lay_var_1965) * 5 + (jt1_var_1952(jl_var_1980, lay_var_1965) - 1)) * nspa_var_216(2) + 1
      inds_var_1962 = indself_var_1957(jl_var_1980, lay_var_1965)
      indf_var_1963 = indfor_var_1958(jl_var_1980, lay_var_1965)
      pp_var_1969 = pavel_var_1941(jl_var_1980, lay_var_1965)
      corradj_var_1968 = 1.0D0 - 0.05D0 * (pp_var_1969 - 100.0D0) / 900.0D0
      DO ig_var_1964 = 1, 12
        tauself_var_1967 = selffac_var_1955(jl_var_1980, lay_var_1965) * (selfref_var_157(inds_var_1962, ig_var_1964) + selffrac_var_1956(jl_var_1980, lay_var_1965) * (selfref_var_157(inds_var_1962 + 1, ig_var_1964) - selfref_var_157(inds_var_1962, ig_var_1964)))
        taufor_var_1966 = forfac_var_1949(jl_var_1980, lay_var_1965) * (forref_var_158(indf_var_1963, ig_var_1964) + forfrac_var_1948(jl_var_1980, lay_var_1965) * (forref_var_158(indf_var_1963 + 1, ig_var_1964) - forref_var_158(indf_var_1963, ig_var_1964)))
        taug_var_1940(jl_var_1980, 10 + ig_var_1964, lay_var_1965) = corradj_var_1968 * (colh2o_var_1953(jl_var_1980, lay_var_1965) * (fac00_var_1944(jl_var_1980, lay_var_1965) * absa_var_155(ind0_var_1960, ig_var_1964) + fac10_var_1946(jl_var_1980, lay_var_1965) * absa_var_155(ind0_var_1960 + 1, ig_var_1964) + fac01_var_1945(jl_var_1980, lay_var_1965) * absa_var_155(ind1_var_1961, ig_var_1964) + fac11_var_1947(jl_var_1980, lay_var_1965) * absa_var_155(ind1_var_1961 + 1, ig_var_1964)) + tauself_var_1967 + taufor_var_1966)
        fracs_var_1959(jl_var_1980, 10 + ig_var_1964, lay_var_1965) = fracrefa_var_153(ig_var_1964)
      END DO
    END DO
  END DO
  DO lay_var_1965 = laytrop_max_var_1971 + 1, klev_var_1939
    DO jl_var_1980 = kidia_var_1937, kfdia_var_1938
      ind0_var_1960 = ((jp_var_1950(jl_var_1980, lay_var_1965) - 13) * 5 + (jt_var_1951(jl_var_1980, lay_var_1965) - 1)) * nspb_var_217(2) + 1
      ind1_var_1961 = ((jp_var_1950(jl_var_1980, lay_var_1965) - 12) * 5 + (jt1_var_1952(jl_var_1980, lay_var_1965) - 1)) * nspb_var_217(2) + 1
      indf_var_1963 = indfor_var_1958(jl_var_1980, lay_var_1965)
      DO ig_var_1964 = 1, 12
        taufor_var_1966 = forfac_var_1949(jl_var_1980, lay_var_1965) * (forref_var_158(indf_var_1963, ig_var_1964) + forfrac_var_1948(jl_var_1980, lay_var_1965) * (forref_var_158(indf_var_1963 + 1, ig_var_1964) - forref_var_158(indf_var_1963, ig_var_1964)))
        taug_var_1940(jl_var_1980, 10 + ig_var_1964, lay_var_1965) = colh2o_var_1953(jl_var_1980, lay_var_1965) * (fac00_var_1944(jl_var_1980, lay_var_1965) * absb_var_156(ind0_var_1960, ig_var_1964) + fac10_var_1946(jl_var_1980, lay_var_1965) * absb_var_156(ind0_var_1960 + 1, ig_var_1964) + fac01_var_1945(jl_var_1980, lay_var_1965) * absb_var_156(ind1_var_1961, ig_var_1964) + fac11_var_1947(jl_var_1980, lay_var_1965) * absb_var_156(ind1_var_1961 + 1, ig_var_1964)) + taufor_var_1966
        fracs_var_1959(jl_var_1980, 10 + ig_var_1964, lay_var_1965) = fracrefb_var_154(ig_var_1964)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1971 /= laytrop_min_var_1970) THEN
    DO lay_var_1965 = laytrop_min_var_1970 + 1, laytrop_max_var_1971
      ixc0_var_1977 = ixc_var_1972(lay_var_1965)
      DO ixp_var_1978 = 1, ixc0_var_1977
        jl_var_1980 = ixlow_var_1973(ixp_var_1978, lay_var_1965)
        ind0_var_1960 = ((jp_var_1950(jl_var_1980, lay_var_1965) - 1) * 5 + (jt_var_1951(jl_var_1980, lay_var_1965) - 1)) * nspa_var_216(2) + 1
        ind1_var_1961 = (jp_var_1950(jl_var_1980, lay_var_1965) * 5 + (jt1_var_1952(jl_var_1980, lay_var_1965) - 1)) * nspa_var_216(2) + 1
        inds_var_1962 = indself_var_1957(jl_var_1980, lay_var_1965)
        indf_var_1963 = indfor_var_1958(jl_var_1980, lay_var_1965)
        pp_var_1969 = pavel_var_1941(jl_var_1980, lay_var_1965)
        corradj_var_1968 = 1.0D0 - 0.05D0 * (pp_var_1969 - 100.0D0) / 900.0D0
        DO ig_var_1964 = 1, 12
          tauself_var_1967 = selffac_var_1955(jl_var_1980, lay_var_1965) * (selfref_var_157(inds_var_1962, ig_var_1964) + selffrac_var_1956(jl_var_1980, lay_var_1965) * (selfref_var_157(inds_var_1962 + 1, ig_var_1964) - selfref_var_157(inds_var_1962, ig_var_1964)))
          taufor_var_1966 = forfac_var_1949(jl_var_1980, lay_var_1965) * (forref_var_158(indf_var_1963, ig_var_1964) + forfrac_var_1948(jl_var_1980, lay_var_1965) * (forref_var_158(indf_var_1963 + 1, ig_var_1964) - forref_var_158(indf_var_1963, ig_var_1964)))
          taug_var_1940(jl_var_1980, 10 + ig_var_1964, lay_var_1965) = corradj_var_1968 * (colh2o_var_1953(jl_var_1980, lay_var_1965) * (fac00_var_1944(jl_var_1980, lay_var_1965) * absa_var_155(ind0_var_1960, ig_var_1964) + fac10_var_1946(jl_var_1980, lay_var_1965) * absa_var_155(ind0_var_1960 + 1, ig_var_1964) + fac01_var_1945(jl_var_1980, lay_var_1965) * absa_var_155(ind1_var_1961, ig_var_1964) + fac11_var_1947(jl_var_1980, lay_var_1965) * absa_var_155(ind1_var_1961 + 1, ig_var_1964)) + tauself_var_1967 + taufor_var_1966)
          fracs_var_1959(jl_var_1980, 10 + ig_var_1964, lay_var_1965) = fracrefa_var_153(ig_var_1964)
        END DO
      END DO
      ixc0_var_1977 = kfdia_var_1938 - kidia_var_1937 + 1 - ixc0_var_1977
      DO ixp_var_1978 = 1, ixc0_var_1977
        jl_var_1980 = ixhigh_var_1974(ixp_var_1978, lay_var_1965)
        ind0_var_1960 = ((jp_var_1950(jl_var_1980, lay_var_1965) - 13) * 5 + (jt_var_1951(jl_var_1980, lay_var_1965) - 1)) * nspb_var_217(2) + 1
        ind1_var_1961 = ((jp_var_1950(jl_var_1980, lay_var_1965) - 12) * 5 + (jt1_var_1952(jl_var_1980, lay_var_1965) - 1)) * nspb_var_217(2) + 1
        indf_var_1963 = indfor_var_1958(jl_var_1980, lay_var_1965)
        DO ig_var_1964 = 1, 12
          taufor_var_1966 = forfac_var_1949(jl_var_1980, lay_var_1965) * (forref_var_158(indf_var_1963, ig_var_1964) + forfrac_var_1948(jl_var_1980, lay_var_1965) * (forref_var_158(indf_var_1963 + 1, ig_var_1964) - forref_var_158(indf_var_1963, ig_var_1964)))
          taug_var_1940(jl_var_1980, 10 + ig_var_1964, lay_var_1965) = colh2o_var_1953(jl_var_1980, lay_var_1965) * (fac00_var_1944(jl_var_1980, lay_var_1965) * absb_var_156(ind0_var_1960, ig_var_1964) + fac10_var_1946(jl_var_1980, lay_var_1965) * absb_var_156(ind0_var_1960 + 1, ig_var_1964) + fac01_var_1945(jl_var_1980, lay_var_1965) * absb_var_156(ind1_var_1961, ig_var_1964) + fac11_var_1947(jl_var_1980, lay_var_1965) * absb_var_156(ind1_var_1961 + 1, ig_var_1964)) + taufor_var_1966
          fracs_var_1959(jl_var_1980, 10 + ig_var_1964, lay_var_1965) = fracrefb_var_154(ig_var_1964)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol2
SUBROUTINE srtm_taumol22(kidia_var_1981, kfdia_var_1982, klev_var_1983, p_fac00_var_1984, p_fac01_var_1985, p_fac10_var_1986, p_fac11_var_1987, k_jp_var_1988, k_jt_var_1989, k_jt1_var_1990, p_oneminus_var_1991, p_colh2o_var_1992, p_colmol_var_1993, p_colo2_var_1994, k_laytrop_var_1995, p_selffac_var_1996, p_selffrac_var_1997, k_indself_var_1998, p_forfac_var_1999, p_forfrac_var_2000, k_indfor_var_2001, p_sfluxzen_var_2002, p_taug_var_2003, p_taur_var_2004, prmu0_var_2005)
  USE yoesrta22, ONLY: absa_var_267, absb_var_268, forrefc_var_270, layreffr_var_266, rayl_var_264, selfrefc_var_269, sfluxrefc_var_271, strrat_var_265
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1981, kfdia_var_1982
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1983
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1984(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1985(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1986(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1987(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1988(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1989(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1990(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1991(kidia_var_1981 : kfdia_var_1982)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1992(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1993(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_1994(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1995(kidia_var_1981 : kfdia_var_1982)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1996(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1997(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1998(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1999(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_2000(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_2001(kidia_var_1981 : kfdia_var_1982, klev_var_1983)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_2002(kidia_var_1981 : kfdia_var_1982, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_2003(kidia_var_1981 : kfdia_var_1982, klev_var_1983, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_2004(kidia_var_1981 : kfdia_var_1982, klev_var_1983, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_2005(kidia_var_1981 : kfdia_var_1982)
  INTEGER(KIND = 4) :: ig_var_2006, ind0_var_2007, ind1_var_2008, inds_var_2009, indf_var_2010, js_var_2011, i_lay_var_2012, i_laysolfr_var_2013(kidia_var_1981 : kfdia_var_1982), i_nlayers_var_2014, iplon_var_2015
  INTEGER(KIND = 4) :: laytrop_min_var_2016, laytrop_max_var_2017
  REAL(KIND = 8) :: z_fs_var_2018, z_speccomb_var_2019, z_specmult_var_2020, z_specparm_var_2021, z_tauray_var_2022, z_o2adj, z_o2cont
  laytrop_min_var_2016 = MINVAL(k_laytrop_var_1995(kidia_var_1981 : kfdia_var_1982))
  laytrop_max_var_2017 = MAXVAL(k_laytrop_var_1995(kidia_var_1981 : kfdia_var_1982))
  i_nlayers_var_2014 = klev_var_1983
  z_o2adj = 1.6D0
  DO iplon_var_2015 = kidia_var_1981, kfdia_var_1982
    i_laysolfr_var_2013(iplon_var_2015) = k_laytrop_var_1995(iplon_var_2015)
  END DO
  DO i_lay_var_2012 = 1, laytrop_min_var_2016
    DO iplon_var_2015 = kidia_var_1981, kfdia_var_1982
      IF (k_jp_var_1988(iplon_var_2015, i_lay_var_2012) < layreffr_var_266 .AND. k_jp_var_1988(iplon_var_2015, i_lay_var_2012 + 1) >= layreffr_var_266) i_laysolfr_var_2013(iplon_var_2015) = MIN(i_lay_var_2012 + 1, k_laytrop_var_1995(iplon_var_2015))
      z_o2cont = 0.000435D0 * p_colo2_var_1994(iplon_var_2015, i_lay_var_2012) / (700.0D0)
      z_speccomb_var_2019 = p_colh2o_var_1992(iplon_var_2015, i_lay_var_2012) + z_o2adj * strrat_var_265 * p_colo2_var_1994(iplon_var_2015, i_lay_var_2012)
      z_specparm_var_2021 = p_colh2o_var_1992(iplon_var_2015, i_lay_var_2012) / z_speccomb_var_2019
      z_specparm_var_2021 = MIN(p_oneminus_var_1991(iplon_var_2015), z_specparm_var_2021)
      z_specmult_var_2020 = 8.0D0 * (z_specparm_var_2021)
      js_var_2011 = 1 + INT(z_specmult_var_2020)
      z_fs_var_2018 = z_specmult_var_2020 - AINT(z_specmult_var_2020)
      ind0_var_2007 = ((k_jp_var_1988(iplon_var_2015, i_lay_var_2012) - 1) * 5 + (k_jt_var_1989(iplon_var_2015, i_lay_var_2012) - 1)) * nspa_var_313(22) + js_var_2011
      ind1_var_2008 = (k_jp_var_1988(iplon_var_2015, i_lay_var_2012) * 5 + (k_jt1_var_1990(iplon_var_2015, i_lay_var_2012) - 1)) * nspa_var_313(22) + js_var_2011
      inds_var_2009 = k_indself_var_1998(iplon_var_2015, i_lay_var_2012)
      indf_var_2010 = k_indfor_var_2001(iplon_var_2015, i_lay_var_2012)
      z_tauray_var_2022 = p_colmol_var_1993(iplon_var_2015, i_lay_var_2012) * rayl_var_264
      DO ig_var_2006 = 1, 2
        p_taug_var_2003(iplon_var_2015, i_lay_var_2012, ig_var_2006) = z_speccomb_var_2019 * ((1.0D0 - z_fs_var_2018) * (absa_var_267(ind0_var_2007, ig_var_2006) * p_fac00_var_1984(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind0_var_2007 + 9, ig_var_2006) * p_fac10_var_1986(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind1_var_2008, ig_var_2006) * p_fac01_var_1985(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind1_var_2008 + 9, ig_var_2006) * p_fac11_var_1987(iplon_var_2015, i_lay_var_2012)) + z_fs_var_2018 * (absa_var_267(ind0_var_2007 + 1, ig_var_2006) * p_fac00_var_1984(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind0_var_2007 + 10, ig_var_2006) * p_fac10_var_1986(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind1_var_2008 + 1, ig_var_2006) * p_fac01_var_1985(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind1_var_2008 + 10, ig_var_2006) * p_fac11_var_1987(iplon_var_2015, i_lay_var_2012))) + p_colh2o_var_1992(iplon_var_2015, i_lay_var_2012) * (p_selffac_var_1996(iplon_var_2015, i_lay_var_2012) * (selfrefc_var_269(inds_var_2009, ig_var_2006) + p_selffrac_var_1997(iplon_var_2015, i_lay_var_2012) * (selfrefc_var_269(inds_var_2009 + 1, ig_var_2006) - selfrefc_var_269(inds_var_2009, ig_var_2006))) + p_forfac_var_1999(iplon_var_2015, i_lay_var_2012) * (forrefc_var_270(indf_var_2010, ig_var_2006) + p_forfrac_var_2000(iplon_var_2015, i_lay_var_2012) * (forrefc_var_270(indf_var_2010 + 1, ig_var_2006) - forrefc_var_270(indf_var_2010, ig_var_2006)))) + z_o2cont
        IF (i_lay_var_2012 == i_laysolfr_var_2013(iplon_var_2015)) p_sfluxzen_var_2002(iplon_var_2015, ig_var_2006) = sfluxrefc_var_271(ig_var_2006, js_var_2011) + z_fs_var_2018 * (sfluxrefc_var_271(ig_var_2006, js_var_2011 + 1) - sfluxrefc_var_271(ig_var_2006, js_var_2011))
        p_taur_var_2004(iplon_var_2015, i_lay_var_2012, ig_var_2006) = z_tauray_var_2022
      END DO
    END DO
  END DO
  DO i_lay_var_2012 = laytrop_min_var_2016 + 1, laytrop_max_var_2017
    DO iplon_var_2015 = kidia_var_1981, kfdia_var_1982
      IF (i_lay_var_2012 <= k_laytrop_var_1995(iplon_var_2015)) THEN
        IF (k_jp_var_1988(iplon_var_2015, i_lay_var_2012) < layreffr_var_266 .AND. k_jp_var_1988(iplon_var_2015, i_lay_var_2012 + 1) >= layreffr_var_266) i_laysolfr_var_2013(iplon_var_2015) = MIN(i_lay_var_2012 + 1, k_laytrop_var_1995(iplon_var_2015))
        z_o2cont = 0.000435D0 * p_colo2_var_1994(iplon_var_2015, i_lay_var_2012) / (700.0D0)
        z_speccomb_var_2019 = p_colh2o_var_1992(iplon_var_2015, i_lay_var_2012) + z_o2adj * strrat_var_265 * p_colo2_var_1994(iplon_var_2015, i_lay_var_2012)
        z_specparm_var_2021 = p_colh2o_var_1992(iplon_var_2015, i_lay_var_2012) / z_speccomb_var_2019
        z_specparm_var_2021 = MIN(p_oneminus_var_1991(iplon_var_2015), z_specparm_var_2021)
        z_specmult_var_2020 = 8.0D0 * (z_specparm_var_2021)
        js_var_2011 = 1 + INT(z_specmult_var_2020)
        z_fs_var_2018 = z_specmult_var_2020 - AINT(z_specmult_var_2020)
        ind0_var_2007 = ((k_jp_var_1988(iplon_var_2015, i_lay_var_2012) - 1) * 5 + (k_jt_var_1989(iplon_var_2015, i_lay_var_2012) - 1)) * nspa_var_313(22) + js_var_2011
        ind1_var_2008 = (k_jp_var_1988(iplon_var_2015, i_lay_var_2012) * 5 + (k_jt1_var_1990(iplon_var_2015, i_lay_var_2012) - 1)) * nspa_var_313(22) + js_var_2011
        inds_var_2009 = k_indself_var_1998(iplon_var_2015, i_lay_var_2012)
        indf_var_2010 = k_indfor_var_2001(iplon_var_2015, i_lay_var_2012)
        z_tauray_var_2022 = p_colmol_var_1993(iplon_var_2015, i_lay_var_2012) * rayl_var_264
        DO ig_var_2006 = 1, 2
          p_taug_var_2003(iplon_var_2015, i_lay_var_2012, ig_var_2006) = z_speccomb_var_2019 * ((1.0D0 - z_fs_var_2018) * (absa_var_267(ind0_var_2007, ig_var_2006) * p_fac00_var_1984(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind0_var_2007 + 9, ig_var_2006) * p_fac10_var_1986(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind1_var_2008, ig_var_2006) * p_fac01_var_1985(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind1_var_2008 + 9, ig_var_2006) * p_fac11_var_1987(iplon_var_2015, i_lay_var_2012)) + z_fs_var_2018 * (absa_var_267(ind0_var_2007 + 1, ig_var_2006) * p_fac00_var_1984(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind0_var_2007 + 10, ig_var_2006) * p_fac10_var_1986(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind1_var_2008 + 1, ig_var_2006) * p_fac01_var_1985(iplon_var_2015, i_lay_var_2012) + absa_var_267(ind1_var_2008 + 10, ig_var_2006) * p_fac11_var_1987(iplon_var_2015, i_lay_var_2012))) + p_colh2o_var_1992(iplon_var_2015, i_lay_var_2012) * (p_selffac_var_1996(iplon_var_2015, i_lay_var_2012) * (selfrefc_var_269(inds_var_2009, ig_var_2006) + p_selffrac_var_1997(iplon_var_2015, i_lay_var_2012) * (selfrefc_var_269(inds_var_2009 + 1, ig_var_2006) - selfrefc_var_269(inds_var_2009, ig_var_2006))) + p_forfac_var_1999(iplon_var_2015, i_lay_var_2012) * (forrefc_var_270(indf_var_2010, ig_var_2006) + p_forfrac_var_2000(iplon_var_2015, i_lay_var_2012) * (forrefc_var_270(indf_var_2010 + 1, ig_var_2006) - forrefc_var_270(indf_var_2010, ig_var_2006)))) + z_o2cont
          IF (i_lay_var_2012 == i_laysolfr_var_2013(iplon_var_2015)) p_sfluxzen_var_2002(iplon_var_2015, ig_var_2006) = sfluxrefc_var_271(ig_var_2006, js_var_2011) + z_fs_var_2018 * (sfluxrefc_var_271(ig_var_2006, js_var_2011 + 1) - sfluxrefc_var_271(ig_var_2006, js_var_2011))
          p_taur_var_2004(iplon_var_2015, i_lay_var_2012, ig_var_2006) = z_tauray_var_2022
        END DO
      ELSE
        z_o2cont = 0.000435D0 * p_colo2_var_1994(iplon_var_2015, i_lay_var_2012) / (700.0D0)
        ind0_var_2007 = ((k_jp_var_1988(iplon_var_2015, i_lay_var_2012) - 13) * 5 + (k_jt_var_1989(iplon_var_2015, i_lay_var_2012) - 1)) * nspb_var_314(22) + 1
        ind1_var_2008 = ((k_jp_var_1988(iplon_var_2015, i_lay_var_2012) - 12) * 5 + (k_jt1_var_1990(iplon_var_2015, i_lay_var_2012) - 1)) * nspb_var_314(22) + 1
        z_tauray_var_2022 = p_colmol_var_1993(iplon_var_2015, i_lay_var_2012) * rayl_var_264
        DO ig_var_2006 = 1, 2
          p_taug_var_2003(iplon_var_2015, i_lay_var_2012, ig_var_2006) = p_colo2_var_1994(iplon_var_2015, i_lay_var_2012) * z_o2adj * (p_fac00_var_1984(iplon_var_2015, i_lay_var_2012) * absb_var_268(ind0_var_2007, ig_var_2006) + p_fac10_var_1986(iplon_var_2015, i_lay_var_2012) * absb_var_268(ind0_var_2007 + 1, ig_var_2006) + p_fac01_var_1985(iplon_var_2015, i_lay_var_2012) * absb_var_268(ind1_var_2008, ig_var_2006) + p_fac11_var_1987(iplon_var_2015, i_lay_var_2012) * absb_var_268(ind1_var_2008 + 1, ig_var_2006)) + z_o2cont
          p_taur_var_2004(iplon_var_2015, i_lay_var_2012, ig_var_2006) = z_tauray_var_2022
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_2012 = laytrop_max_var_2017 + 1, i_nlayers_var_2014
    DO iplon_var_2015 = kidia_var_1981, kfdia_var_1982
      z_o2cont = 0.000435D0 * p_colo2_var_1994(iplon_var_2015, i_lay_var_2012) / (700.0D0)
      ind0_var_2007 = ((k_jp_var_1988(iplon_var_2015, i_lay_var_2012) - 13) * 5 + (k_jt_var_1989(iplon_var_2015, i_lay_var_2012) - 1)) * nspb_var_314(22) + 1
      ind1_var_2008 = ((k_jp_var_1988(iplon_var_2015, i_lay_var_2012) - 12) * 5 + (k_jt1_var_1990(iplon_var_2015, i_lay_var_2012) - 1)) * nspb_var_314(22) + 1
      z_tauray_var_2022 = p_colmol_var_1993(iplon_var_2015, i_lay_var_2012) * rayl_var_264
      DO ig_var_2006 = 1, 2
        p_taug_var_2003(iplon_var_2015, i_lay_var_2012, ig_var_2006) = p_colo2_var_1994(iplon_var_2015, i_lay_var_2012) * z_o2adj * (p_fac00_var_1984(iplon_var_2015, i_lay_var_2012) * absb_var_268(ind0_var_2007, ig_var_2006) + p_fac10_var_1986(iplon_var_2015, i_lay_var_2012) * absb_var_268(ind0_var_2007 + 1, ig_var_2006) + p_fac01_var_1985(iplon_var_2015, i_lay_var_2012) * absb_var_268(ind1_var_2008, ig_var_2006) + p_fac11_var_1987(iplon_var_2015, i_lay_var_2012) * absb_var_268(ind1_var_2008 + 1, ig_var_2006)) + z_o2cont
        p_taur_var_2004(iplon_var_2015, i_lay_var_2012, ig_var_2006) = z_tauray_var_2022
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol22
SUBROUTINE rrtm_taumol6(kidia_var_2023, kfdia_var_2024, klev_var_2025, taug_var_2026, wx_var_2027, p_tauaerl_var_2028, fac00_var_2029, fac01_var_2030, fac10_var_2031, fac11_var_2032, forfac_var_2045, forfrac_var_2046, indfor_var_2044, jp_var_2033, jt_var_2034, jt1_var_2035, colh2o_var_2036, colco2_var_2037, coldry_var_2038, laytrop_var_2039, selffac_var_2040, selffrac_var_2041, indself_var_2042, fracs_var_2043, minorfrac_var_2047, indminor_var_2048)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrta6, ONLY: absa_var_182, cfc11adj, cfc12_var_181, forref_var_185, fracrefa_var_180, ka_mco2_var_184, selfref_var_183
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2023
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2024
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2025
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2026(kidia_var_2023 : kfdia_var_2024, 140, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: wx_var_2027(kidia_var_2023 : kfdia_var_2024, 4, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2028(kidia_var_2023 : kfdia_var_2024, klev_var_2025, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2029(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2030(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2031(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2032(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2033(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2034(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2035(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2036(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2037(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_2038(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2039(kidia_var_2023 : kfdia_var_2024)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2040(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2041(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2042(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2043(kidia_var_2023 : kfdia_var_2024, 140, klev_var_2025)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2044(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2045(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2046(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2047(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2048(kidia_var_2023 : kfdia_var_2024, klev_var_2025)
  INTEGER(KIND = 4) :: ind0_var_2049, ind1_var_2050, inds_var_2051, indf_var_2052, indm_var_2053
  INTEGER(KIND = 4) :: ig_var_2054, lay_var_2055
  REAL(KIND = 8) :: adjfac_var_2056, adjcolco2_var_2057, ratco2_var_2058, chi_co2_var_2059
  REAL(KIND = 8) :: taufor_var_2060, tauself_var_2061, absco2_var_2062
  INTEGER(KIND = 4) :: laytrop_min_var_2063, laytrop_max_var_2064
  INTEGER(KIND = 4) :: ixc_var_2065(klev_var_2025), ixlow_var_2066(kfdia_var_2024, klev_var_2025), ixhigh_var_2067(kfdia_var_2024, klev_var_2025)
  INTEGER(KIND = 4) :: ich_var_2068, icl_var_2069, ixc0_var_2070, ixp_var_2071, jc_var_2072, jl_var_2073
  laytrop_min_var_2063 = MINVAL(laytrop_var_2039)
  laytrop_max_var_2064 = MAXVAL(laytrop_var_2039)
  ixlow_var_2066 = 0
  ixhigh_var_2067 = 0
  ixc_var_2065 = 0
  DO lay_var_2055 = laytrop_min_var_2063 + 1, laytrop_max_var_2064
    icl_var_2069 = 0
    ich_var_2068 = 0
    DO jc_var_2072 = kidia_var_2023, kfdia_var_2024
      IF (lay_var_2055 <= laytrop_var_2039(jc_var_2072)) THEN
        icl_var_2069 = icl_var_2069 + 1
        ixlow_var_2066(icl_var_2069, lay_var_2055) = jc_var_2072
      ELSE
        ich_var_2068 = ich_var_2068 + 1
        ixhigh_var_2067(ich_var_2068, lay_var_2055) = jc_var_2072
      END IF
    END DO
    ixc_var_2065(lay_var_2055) = icl_var_2069
  END DO
  DO lay_var_2055 = 1, laytrop_min_var_2063
    DO jl_var_2073 = kidia_var_2023, kfdia_var_2024
      chi_co2_var_2059 = colco2_var_2037(jl_var_2073, lay_var_2055) / (coldry_var_2038(jl_var_2073, lay_var_2055))
      ratco2_var_2058 = 1D+20 * chi_co2_var_2059 / chi_mls(2, jp_var_2033(jl_var_2073, lay_var_2055) + 1)
      IF (ratco2_var_2058 .GT. 3.0D0) THEN
        adjfac_var_2056 = 2.0D0 + (ratco2_var_2058 - 2.0D0) ** 0.77D0
        adjcolco2_var_2057 = adjfac_var_2056 * chi_mls(2, jp_var_2033(jl_var_2073, lay_var_2055) + 1) * coldry_var_2038(jl_var_2073, lay_var_2055) * 1D-20
      ELSE
        adjcolco2_var_2057 = colco2_var_2037(jl_var_2073, lay_var_2055)
      END IF
      ind0_var_2049 = ((jp_var_2033(jl_var_2073, lay_var_2055) - 1) * 5 + (jt_var_2034(jl_var_2073, lay_var_2055) - 1)) * nspa_var_216(6) + 1
      ind1_var_2050 = (jp_var_2033(jl_var_2073, lay_var_2055) * 5 + (jt1_var_2035(jl_var_2073, lay_var_2055) - 1)) * nspa_var_216(6) + 1
      inds_var_2051 = indself_var_2042(jl_var_2073, lay_var_2055)
      indf_var_2052 = indfor_var_2044(jl_var_2073, lay_var_2055)
      indm_var_2053 = indminor_var_2048(jl_var_2073, lay_var_2055)
      DO ig_var_2054 = 1, 8
        tauself_var_2061 = selffac_var_2040(jl_var_2073, lay_var_2055) * (selfref_var_183(inds_var_2051, ig_var_2054) + selffrac_var_2041(jl_var_2073, lay_var_2055) * (selfref_var_183(inds_var_2051 + 1, ig_var_2054) - selfref_var_183(inds_var_2051, ig_var_2054)))
        taufor_var_2060 = forfac_var_2045(jl_var_2073, lay_var_2055) * (forref_var_185(indf_var_2052, ig_var_2054) + forfrac_var_2046(jl_var_2073, lay_var_2055) * (forref_var_185(indf_var_2052 + 1, ig_var_2054) - forref_var_185(indf_var_2052, ig_var_2054)))
        absco2_var_2062 = (ka_mco2_var_184(indm_var_2053, ig_var_2054) + minorfrac_var_2047(jl_var_2073, lay_var_2055) * (ka_mco2_var_184(indm_var_2053 + 1, ig_var_2054) - ka_mco2_var_184(indm_var_2053, ig_var_2054)))
        taug_var_2026(jl_var_2073, 68 + ig_var_2054, lay_var_2055) = colh2o_var_2036(jl_var_2073, lay_var_2055) * (fac00_var_2029(jl_var_2073, lay_var_2055) * absa_var_182(ind0_var_2049, ig_var_2054) + fac10_var_2031(jl_var_2073, lay_var_2055) * absa_var_182(ind0_var_2049 + 1, ig_var_2054) + fac01_var_2030(jl_var_2073, lay_var_2055) * absa_var_182(ind1_var_2050, ig_var_2054) + fac11_var_2032(jl_var_2073, lay_var_2055) * absa_var_182(ind1_var_2050 + 1, ig_var_2054)) + tauself_var_2061 + taufor_var_2060 + adjcolco2_var_2057 * absco2_var_2062 + wx_var_2027(jl_var_2073, 2, lay_var_2055) * cfc11adj(ig_var_2054) + wx_var_2027(jl_var_2073, 3, lay_var_2055) * cfc12_var_181(ig_var_2054)
        fracs_var_2043(jl_var_2073, 68 + ig_var_2054, lay_var_2055) = fracrefa_var_180(ig_var_2054)
      END DO
    END DO
  END DO
  DO ig_var_2054 = 1, 8
    DO lay_var_2055 = laytrop_max_var_2064 + 1, klev_var_2025
      DO jl_var_2073 = kidia_var_2023, kfdia_var_2024
        taug_var_2026(jl_var_2073, 68 + ig_var_2054, lay_var_2055) = 0.0D0 + wx_var_2027(jl_var_2073, 2, lay_var_2055) * cfc11adj(ig_var_2054) + wx_var_2027(jl_var_2073, 3, lay_var_2055) * cfc12_var_181(ig_var_2054)
        fracs_var_2043(jl_var_2073, 68 + ig_var_2054, lay_var_2055) = fracrefa_var_180(ig_var_2054)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2064 /= laytrop_min_var_2063) THEN
    DO lay_var_2055 = laytrop_min_var_2063 + 1, laytrop_max_var_2064
      ixc0_var_2070 = ixc_var_2065(lay_var_2055)
      DO ixp_var_2071 = 1, ixc0_var_2070
        jl_var_2073 = ixlow_var_2066(ixp_var_2071, lay_var_2055)
        chi_co2_var_2059 = colco2_var_2037(jl_var_2073, lay_var_2055) / (coldry_var_2038(jl_var_2073, lay_var_2055))
        ratco2_var_2058 = 1D+20 * chi_co2_var_2059 / chi_mls(2, jp_var_2033(jl_var_2073, lay_var_2055) + 1)
        IF (ratco2_var_2058 .GT. 3.0D0) THEN
          adjfac_var_2056 = 2.0D0 + (ratco2_var_2058 - 2.0D0) ** 0.77D0
          adjcolco2_var_2057 = adjfac_var_2056 * chi_mls(2, jp_var_2033(jl_var_2073, lay_var_2055) + 1) * coldry_var_2038(jl_var_2073, lay_var_2055) * 1D-20
        ELSE
          adjcolco2_var_2057 = colco2_var_2037(jl_var_2073, lay_var_2055)
        END IF
        ind0_var_2049 = ((jp_var_2033(jl_var_2073, lay_var_2055) - 1) * 5 + (jt_var_2034(jl_var_2073, lay_var_2055) - 1)) * nspa_var_216(6) + 1
        ind1_var_2050 = (jp_var_2033(jl_var_2073, lay_var_2055) * 5 + (jt1_var_2035(jl_var_2073, lay_var_2055) - 1)) * nspa_var_216(6) + 1
        inds_var_2051 = indself_var_2042(jl_var_2073, lay_var_2055)
        indf_var_2052 = indfor_var_2044(jl_var_2073, lay_var_2055)
        indm_var_2053 = indminor_var_2048(jl_var_2073, lay_var_2055)
        DO ig_var_2054 = 1, 8
          tauself_var_2061 = selffac_var_2040(jl_var_2073, lay_var_2055) * (selfref_var_183(inds_var_2051, ig_var_2054) + selffrac_var_2041(jl_var_2073, lay_var_2055) * (selfref_var_183(inds_var_2051 + 1, ig_var_2054) - selfref_var_183(inds_var_2051, ig_var_2054)))
          taufor_var_2060 = forfac_var_2045(jl_var_2073, lay_var_2055) * (forref_var_185(indf_var_2052, ig_var_2054) + forfrac_var_2046(jl_var_2073, lay_var_2055) * (forref_var_185(indf_var_2052 + 1, ig_var_2054) - forref_var_185(indf_var_2052, ig_var_2054)))
          absco2_var_2062 = (ka_mco2_var_184(indm_var_2053, ig_var_2054) + minorfrac_var_2047(jl_var_2073, lay_var_2055) * (ka_mco2_var_184(indm_var_2053 + 1, ig_var_2054) - ka_mco2_var_184(indm_var_2053, ig_var_2054)))
          taug_var_2026(jl_var_2073, 68 + ig_var_2054, lay_var_2055) = colh2o_var_2036(jl_var_2073, lay_var_2055) * (fac00_var_2029(jl_var_2073, lay_var_2055) * absa_var_182(ind0_var_2049, ig_var_2054) + fac10_var_2031(jl_var_2073, lay_var_2055) * absa_var_182(ind0_var_2049 + 1, ig_var_2054) + fac01_var_2030(jl_var_2073, lay_var_2055) * absa_var_182(ind1_var_2050, ig_var_2054) + fac11_var_2032(jl_var_2073, lay_var_2055) * absa_var_182(ind1_var_2050 + 1, ig_var_2054)) + tauself_var_2061 + taufor_var_2060 + adjcolco2_var_2057 * absco2_var_2062 + wx_var_2027(jl_var_2073, 2, lay_var_2055) * cfc11adj(ig_var_2054) + wx_var_2027(jl_var_2073, 3, lay_var_2055) * cfc12_var_181(ig_var_2054)
          fracs_var_2043(jl_var_2073, 68 + ig_var_2054, lay_var_2055) = fracrefa_var_180(ig_var_2054)
        END DO
      END DO
      ixc0_var_2070 = kfdia_var_2024 - kidia_var_2023 + 1 - ixc0_var_2070
      DO ig_var_2054 = 1, 8
        DO ixp_var_2071 = 1, ixc0_var_2070
          jl_var_2073 = ixhigh_var_2067(ixp_var_2071, lay_var_2055)
          taug_var_2026(jl_var_2073, 68 + ig_var_2054, lay_var_2055) = 0.0D0 + wx_var_2027(jl_var_2073, 2, lay_var_2055) * cfc11adj(ig_var_2054) + wx_var_2027(jl_var_2073, 3, lay_var_2055) * cfc12_var_181(ig_var_2054)
          fracs_var_2043(jl_var_2073, 68 + ig_var_2054, lay_var_2055) = fracrefa_var_180(ig_var_2054)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol6
SUBROUTINE rrtm_taumol7(kidia_var_2074, kfdia_var_2075, klev_var_2076, taug_var_2077, p_tauaerl_var_2078, fac00_var_2079, fac01_var_2080, fac10_var_2081, fac11_var_2082, forfac_var_2098, forfrac_var_2097, indfor_var_2096, jp_var_2083, jt_var_2084, jt1_var_2085, oneminus_var_2086, colh2o_var_2087, colo3_var_2088, colco2_var_2089, coldry_var_2090, laytrop_var_2091, selffac_var_2092, selffrac_var_2093, indself_var_2094, fracs_var_2095, rat_h2oo3, rat_h2oo3_1, minorfrac_var_2099, indminor_var_2100)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng7
  USE yoerrta7, ONLY: absa_var_188, absb_var_189, forref_var_193, fracrefa_var_186, fracrefb_var_187, ka_mco2_var_191, kb_mco2_var_192, selfref_var_190
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2074
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2075
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2076
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2077(kidia_var_2074 : kfdia_var_2075, 140, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2078(kidia_var_2074 : kfdia_var_2075, klev_var_2076, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2079(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2080(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2081(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2082(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2083(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2084(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2085(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2086
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2087(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_2088(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2089(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_2090(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2091(kidia_var_2074 : kfdia_var_2075)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2092(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2093(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2094(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2095(kidia_var_2074 : kfdia_var_2075, 140, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oo3(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oo3_1(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2096(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2097(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2098(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2099(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2100(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8) :: speccomb_var_2101, speccomb1_var_2102, speccomb_mco2_var_2103, speccomb_planck_var_2104
  INTEGER(KIND = 4) :: ind0_var_2105, ind1_var_2106, inds_var_2107, indf_var_2108, indm_var_2109
  INTEGER(KIND = 4) :: ig_var_2110, js_var_2111, lay_var_2112, js1_var_2113, jpl_var_2114, jmco2_var_2115
  REAL(KIND = 8) :: refrat_planck_a_var_2116, refrat_m_a_var_2117
  REAL(KIND = 8) :: chi_co2_var_2118, ratco2_var_2119, adjfac_var_2120, adjcolco2_var_2121
  REAL(KIND = 8) :: fac000_var_2122, fac100_var_2123, fac200_var_2124, fac010_var_2125, fac110_var_2126, fac210_var_2127, fac001_var_2128, fac101_var_2129, fac201_var_2130, fac011_var_2131, fac111_var_2132, fac211_var_2133
  REAL(KIND = 8) :: p_var_2134, p4_var_2135, fk0_var_2136, fk1_var_2137, fk2_var_2138
  REAL(KIND = 8) :: taufor_var_2139, tauself_var_2140, tau_major_var_2141(12), tau_major1_var_2142(12), co2m1_var_2143, co2m2_var_2144, absco2_var_2145
  REAL(KIND = 8) :: fs_var_2146, specmult_var_2147, specparm_var_2148, fs1_var_2149, specmult1_var_2150, specparm1_var_2151, fpl_var_2152, specmult_planck_var_2153, specparm_planck_var_2154, fmco2_var_2155, specmult_mco2_var_2156, specparm_mco2_var_2157
  INTEGER(KIND = 4) :: laytrop_min_var_2158, laytrop_max_var_2159
  INTEGER(KIND = 4) :: ixc_var_2160(klev_var_2076), ixlow_var_2161(kfdia_var_2075, klev_var_2076), ixhigh_var_2162(kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4) :: ich_var_2163, icl_var_2164, ixc0_var_2165, ixp_var_2166, jc_var_2167, jl_var_2168
  laytrop_min_var_2158 = MINVAL(laytrop_var_2091)
  laytrop_max_var_2159 = MAXVAL(laytrop_var_2091)
  ixlow_var_2161 = 0
  ixhigh_var_2162 = 0
  ixc_var_2160 = 0
  DO lay_var_2112 = laytrop_min_var_2158 + 1, laytrop_max_var_2159
    icl_var_2164 = 0
    ich_var_2163 = 0
    DO jc_var_2167 = kidia_var_2074, kfdia_var_2075
      IF (lay_var_2112 <= laytrop_var_2091(jc_var_2167)) THEN
        icl_var_2164 = icl_var_2164 + 1
        ixlow_var_2161(icl_var_2164, lay_var_2112) = jc_var_2167
      ELSE
        ich_var_2163 = ich_var_2163 + 1
        ixhigh_var_2162(ich_var_2163, lay_var_2112) = jc_var_2167
      END IF
    END DO
    ixc_var_2160(lay_var_2112) = icl_var_2164
  END DO
  refrat_planck_a_var_2116 = chi_mls(1, 3) / chi_mls(3, 3)
  refrat_m_a_var_2117 = chi_mls(1, 3) / chi_mls(3, 3)
  DO lay_var_2112 = 1, laytrop_min_var_2158
    DO jl_var_2168 = kidia_var_2074, kfdia_var_2075
      speccomb_var_2101 = colh2o_var_2087(jl_var_2168, lay_var_2112) + rat_h2oo3(jl_var_2168, lay_var_2112) * colo3_var_2088(jl_var_2168, lay_var_2112)
      specparm_var_2148 = MIN(colh2o_var_2087(jl_var_2168, lay_var_2112) / speccomb_var_2101, oneminus_var_2086)
      specmult_var_2147 = 8.0D0 * (specparm_var_2148)
      js_var_2111 = 1 + INT(specmult_var_2147)
      fs_var_2146 = ((specmult_var_2147) - AINT((specmult_var_2147)))
      speccomb1_var_2102 = colh2o_var_2087(jl_var_2168, lay_var_2112) + rat_h2oo3_1(jl_var_2168, lay_var_2112) * colo3_var_2088(jl_var_2168, lay_var_2112)
      specparm1_var_2151 = MIN(colh2o_var_2087(jl_var_2168, lay_var_2112) / speccomb1_var_2102, oneminus_var_2086)
      specmult1_var_2150 = 8.0D0 * (specparm1_var_2151)
      js1_var_2113 = 1 + INT(specmult1_var_2150)
      fs1_var_2149 = ((specmult1_var_2150) - AINT((specmult1_var_2150)))
      speccomb_mco2_var_2103 = colh2o_var_2087(jl_var_2168, lay_var_2112) + refrat_m_a_var_2117 * colo3_var_2088(jl_var_2168, lay_var_2112)
      specparm_mco2_var_2157 = MIN(colh2o_var_2087(jl_var_2168, lay_var_2112) / speccomb_mco2_var_2103, oneminus_var_2086)
      specmult_mco2_var_2156 = 8.0D0 * specparm_mco2_var_2157
      jmco2_var_2115 = 1 + INT(specmult_mco2_var_2156)
      fmco2_var_2155 = ((specmult_mco2_var_2156) - AINT((specmult_mco2_var_2156)))
      chi_co2_var_2118 = colco2_var_2089(jl_var_2168, lay_var_2112) / (coldry_var_2090(jl_var_2168, lay_var_2112))
      ratco2_var_2119 = 1D+20 * chi_co2_var_2118 / chi_mls(2, jp_var_2083(jl_var_2168, lay_var_2112) + 1)
      IF (ratco2_var_2119 .GT. 3.0D0) THEN
        adjfac_var_2120 = 3.0D0 + (ratco2_var_2119 - 3.0D0) ** 0.79D0
        adjcolco2_var_2121 = adjfac_var_2120 * chi_mls(2, jp_var_2083(jl_var_2168, lay_var_2112) + 1) * coldry_var_2090(jl_var_2168, lay_var_2112) * 1D-20
      ELSE
        adjcolco2_var_2121 = colco2_var_2089(jl_var_2168, lay_var_2112)
      END IF
      speccomb_planck_var_2104 = colh2o_var_2087(jl_var_2168, lay_var_2112) + refrat_planck_a_var_2116 * colo3_var_2088(jl_var_2168, lay_var_2112)
      specparm_planck_var_2154 = MIN(colh2o_var_2087(jl_var_2168, lay_var_2112) / speccomb_planck_var_2104, oneminus_var_2086)
      specmult_planck_var_2153 = 8.0D0 * specparm_planck_var_2154
      jpl_var_2114 = 1 + INT(specmult_planck_var_2153)
      fpl_var_2152 = ((specmult_planck_var_2153) - AINT((specmult_planck_var_2153)))
      ind0_var_2105 = ((jp_var_2083(jl_var_2168, lay_var_2112) - 1) * 5 + (jt_var_2084(jl_var_2168, lay_var_2112) - 1)) * nspa_var_216(7) + js_var_2111
      ind1_var_2106 = (jp_var_2083(jl_var_2168, lay_var_2112) * 5 + (jt1_var_2085(jl_var_2168, lay_var_2112) - 1)) * nspa_var_216(7) + js1_var_2113
      inds_var_2107 = indself_var_2094(jl_var_2168, lay_var_2112)
      indf_var_2108 = indfor_var_2096(jl_var_2168, lay_var_2112)
      indm_var_2109 = indminor_var_2100(jl_var_2168, lay_var_2112)
      IF (specparm_var_2148 .LT. 0.125D0) THEN
        p_var_2134 = fs_var_2146 - 1.0D0
        p4_var_2135 = p_var_2134 ** 4
        fk0_var_2136 = p4_var_2135
        fk1_var_2137 = 1.0D0 - p_var_2134 - 2.0D0 * p4_var_2135
        fk2_var_2138 = p_var_2134 + p4_var_2135
        fac000_var_2122 = fk0_var_2136 * fac00_var_2079(jl_var_2168, lay_var_2112)
        fac100_var_2123 = fk1_var_2137 * fac00_var_2079(jl_var_2168, lay_var_2112)
        fac200_var_2124 = fk2_var_2138 * fac00_var_2079(jl_var_2168, lay_var_2112)
        fac010_var_2125 = fk0_var_2136 * fac10_var_2081(jl_var_2168, lay_var_2112)
        fac110_var_2126 = fk1_var_2137 * fac10_var_2081(jl_var_2168, lay_var_2112)
        fac210_var_2127 = fk2_var_2138 * fac10_var_2081(jl_var_2168, lay_var_2112)
      ELSE IF (specparm_var_2148 .GT. 0.875D0) THEN
        p_var_2134 = - fs_var_2146
        p4_var_2135 = p_var_2134 ** 4
        fk0_var_2136 = p4_var_2135
        fk1_var_2137 = 1.0D0 - p_var_2134 - 2.0D0 * p4_var_2135
        fk2_var_2138 = p_var_2134 + p4_var_2135
        fac000_var_2122 = fk0_var_2136 * fac00_var_2079(jl_var_2168, lay_var_2112)
        fac100_var_2123 = fk1_var_2137 * fac00_var_2079(jl_var_2168, lay_var_2112)
        fac200_var_2124 = fk2_var_2138 * fac00_var_2079(jl_var_2168, lay_var_2112)
        fac010_var_2125 = fk0_var_2136 * fac10_var_2081(jl_var_2168, lay_var_2112)
        fac110_var_2126 = fk1_var_2137 * fac10_var_2081(jl_var_2168, lay_var_2112)
        fac210_var_2127 = fk2_var_2138 * fac10_var_2081(jl_var_2168, lay_var_2112)
      ELSE
        fac000_var_2122 = (1.0D0 - fs_var_2146) * fac00_var_2079(jl_var_2168, lay_var_2112)
        fac010_var_2125 = (1.0D0 - fs_var_2146) * fac10_var_2081(jl_var_2168, lay_var_2112)
        fac100_var_2123 = fs_var_2146 * fac00_var_2079(jl_var_2168, lay_var_2112)
        fac110_var_2126 = fs_var_2146 * fac10_var_2081(jl_var_2168, lay_var_2112)
        fac200_var_2124 = 0.0D0
        fac210_var_2127 = 0.0D0
      END IF
      IF (specparm1_var_2151 .LT. 0.125D0) THEN
        p_var_2134 = fs1_var_2149 - 1.0D0
        p4_var_2135 = p_var_2134 ** 4
        fk0_var_2136 = p4_var_2135
        fk1_var_2137 = 1.0D0 - p_var_2134 - 2.0D0 * p4_var_2135
        fk2_var_2138 = p_var_2134 + p4_var_2135
        fac001_var_2128 = fk0_var_2136 * fac01_var_2080(jl_var_2168, lay_var_2112)
        fac101_var_2129 = fk1_var_2137 * fac01_var_2080(jl_var_2168, lay_var_2112)
        fac201_var_2130 = fk2_var_2138 * fac01_var_2080(jl_var_2168, lay_var_2112)
        fac011_var_2131 = fk0_var_2136 * fac11_var_2082(jl_var_2168, lay_var_2112)
        fac111_var_2132 = fk1_var_2137 * fac11_var_2082(jl_var_2168, lay_var_2112)
        fac211_var_2133 = fk2_var_2138 * fac11_var_2082(jl_var_2168, lay_var_2112)
      ELSE IF (specparm1_var_2151 .GT. 0.875D0) THEN
        p_var_2134 = - fs1_var_2149
        p4_var_2135 = p_var_2134 ** 4
        fk0_var_2136 = p4_var_2135
        fk1_var_2137 = 1.0D0 - p_var_2134 - 2.0D0 * p4_var_2135
        fk2_var_2138 = p_var_2134 + p4_var_2135
        fac001_var_2128 = fk0_var_2136 * fac01_var_2080(jl_var_2168, lay_var_2112)
        fac101_var_2129 = fk1_var_2137 * fac01_var_2080(jl_var_2168, lay_var_2112)
        fac201_var_2130 = fk2_var_2138 * fac01_var_2080(jl_var_2168, lay_var_2112)
        fac011_var_2131 = fk0_var_2136 * fac11_var_2082(jl_var_2168, lay_var_2112)
        fac111_var_2132 = fk1_var_2137 * fac11_var_2082(jl_var_2168, lay_var_2112)
        fac211_var_2133 = fk2_var_2138 * fac11_var_2082(jl_var_2168, lay_var_2112)
      ELSE
        fac001_var_2128 = (1.0D0 - fs1_var_2149) * fac01_var_2080(jl_var_2168, lay_var_2112)
        fac011_var_2131 = (1.0D0 - fs1_var_2149) * fac11_var_2082(jl_var_2168, lay_var_2112)
        fac101_var_2129 = fs1_var_2149 * fac01_var_2080(jl_var_2168, lay_var_2112)
        fac111_var_2132 = fs1_var_2149 * fac11_var_2082(jl_var_2168, lay_var_2112)
        fac201_var_2130 = 0.0D0
        fac211_var_2133 = 0.0D0
      END IF
      IF (specparm_var_2148 .LT. 0.125D0) THEN
        tau_major_var_2141(1 : ng7) = speccomb_var_2101 * (fac000_var_2122 * absa_var_188(ind0_var_2105, 1 : 12) + fac100_var_2123 * absa_var_188(ind0_var_2105 + 1, 1 : 12) + fac200_var_2124 * absa_var_188(ind0_var_2105 + 2, 1 : 12) + fac010_var_2125 * absa_var_188(ind0_var_2105 + 9, 1 : 12) + fac110_var_2126 * absa_var_188(ind0_var_2105 + 10, 1 : 12) + fac210_var_2127 * absa_var_188(ind0_var_2105 + 11, 1 : 12))
      ELSE IF (specparm_var_2148 .GT. 0.875D0) THEN
        tau_major_var_2141(1 : ng7) = speccomb_var_2101 * (fac200_var_2124 * absa_var_188(ind0_var_2105 - 1, 1 : 12) + fac100_var_2123 * absa_var_188(ind0_var_2105, 1 : 12) + fac000_var_2122 * absa_var_188(ind0_var_2105 + 1, 1 : 12) + fac210_var_2127 * absa_var_188(ind0_var_2105 + 8, 1 : 12) + fac110_var_2126 * absa_var_188(ind0_var_2105 + 9, 1 : 12) + fac010_var_2125 * absa_var_188(ind0_var_2105 + 10, 1 : 12))
      ELSE
        tau_major_var_2141(1 : ng7) = speccomb_var_2101 * (fac000_var_2122 * absa_var_188(ind0_var_2105, 1 : 12) + fac100_var_2123 * absa_var_188(ind0_var_2105 + 1, 1 : 12) + fac010_var_2125 * absa_var_188(ind0_var_2105 + 9, 1 : 12) + fac110_var_2126 * absa_var_188(ind0_var_2105 + 10, 1 : 12))
      END IF
      IF (specparm1_var_2151 .LT. 0.125D0) THEN
        tau_major1_var_2142(1 : ng7) = speccomb1_var_2102 * (fac001_var_2128 * absa_var_188(ind1_var_2106, 1 : 12) + fac101_var_2129 * absa_var_188(ind1_var_2106 + 1, 1 : 12) + fac201_var_2130 * absa_var_188(ind1_var_2106 + 2, 1 : 12) + fac011_var_2131 * absa_var_188(ind1_var_2106 + 9, 1 : 12) + fac111_var_2132 * absa_var_188(ind1_var_2106 + 10, 1 : 12) + fac211_var_2133 * absa_var_188(ind1_var_2106 + 11, 1 : 12))
      ELSE IF (specparm1_var_2151 .GT. 0.875D0) THEN
        tau_major1_var_2142(1 : ng7) = speccomb1_var_2102 * (fac201_var_2130 * absa_var_188(ind1_var_2106 - 1, 1 : 12) + fac101_var_2129 * absa_var_188(ind1_var_2106, 1 : 12) + fac001_var_2128 * absa_var_188(ind1_var_2106 + 1, 1 : 12) + fac211_var_2133 * absa_var_188(ind1_var_2106 + 8, 1 : 12) + fac111_var_2132 * absa_var_188(ind1_var_2106 + 9, 1 : 12) + fac011_var_2131 * absa_var_188(ind1_var_2106 + 10, 1 : 12))
      ELSE
        tau_major1_var_2142(1 : ng7) = speccomb1_var_2102 * (fac001_var_2128 * absa_var_188(ind1_var_2106, 1 : 12) + fac101_var_2129 * absa_var_188(ind1_var_2106 + 1, 1 : 12) + fac011_var_2131 * absa_var_188(ind1_var_2106 + 9, 1 : 12) + fac111_var_2132 * absa_var_188(ind1_var_2106 + 10, 1 : 12))
      END IF
      DO ig_var_2110 = 1, 12
        tauself_var_2140 = selffac_var_2092(jl_var_2168, lay_var_2112) * (selfref_var_190(inds_var_2107, ig_var_2110) + selffrac_var_2093(jl_var_2168, lay_var_2112) * (selfref_var_190(inds_var_2107 + 1, ig_var_2110) - selfref_var_190(inds_var_2107, ig_var_2110)))
        taufor_var_2139 = forfac_var_2098(jl_var_2168, lay_var_2112) * (forref_var_193(indf_var_2108, ig_var_2110) + forfrac_var_2097(jl_var_2168, lay_var_2112) * (forref_var_193(indf_var_2108 + 1, ig_var_2110) - forref_var_193(indf_var_2108, ig_var_2110)))
        co2m1_var_2143 = ka_mco2_var_191(jmco2_var_2115, indm_var_2109, ig_var_2110) + fmco2_var_2155 * (ka_mco2_var_191(jmco2_var_2115 + 1, indm_var_2109, ig_var_2110) - ka_mco2_var_191(jmco2_var_2115, indm_var_2109, ig_var_2110))
        co2m2_var_2144 = ka_mco2_var_191(jmco2_var_2115, indm_var_2109 + 1, ig_var_2110) + fmco2_var_2155 * (ka_mco2_var_191(jmco2_var_2115 + 1, indm_var_2109 + 1, ig_var_2110) - ka_mco2_var_191(jmco2_var_2115, indm_var_2109 + 1, ig_var_2110))
        absco2_var_2145 = co2m1_var_2143 + minorfrac_var_2099(jl_var_2168, lay_var_2112) * (co2m2_var_2144 - co2m1_var_2143)
        taug_var_2077(jl_var_2168, 76 + ig_var_2110, lay_var_2112) = tau_major_var_2141(ig_var_2110) + tau_major1_var_2142(ig_var_2110) + tauself_var_2140 + taufor_var_2139 + adjcolco2_var_2121 * absco2_var_2145
        fracs_var_2095(jl_var_2168, 76 + ig_var_2110, lay_var_2112) = fracrefa_var_186(ig_var_2110, jpl_var_2114) + fpl_var_2152 * (fracrefa_var_186(ig_var_2110, jpl_var_2114 + 1) - fracrefa_var_186(ig_var_2110, jpl_var_2114))
      END DO
    END DO
  END DO
  DO lay_var_2112 = laytrop_max_var_2159 + 1, klev_var_2076
    DO jl_var_2168 = kidia_var_2074, kfdia_var_2075
      chi_co2_var_2118 = colco2_var_2089(jl_var_2168, lay_var_2112) / (coldry_var_2090(jl_var_2168, lay_var_2112))
      ratco2_var_2119 = 1D+20 * chi_co2_var_2118 / chi_mls(2, jp_var_2083(jl_var_2168, lay_var_2112) + 1)
      IF (ratco2_var_2119 .GT. 3.0D0) THEN
        adjfac_var_2120 = 2.0D0 + (ratco2_var_2119 - 2.0D0) ** 0.79D0
        adjcolco2_var_2121 = adjfac_var_2120 * chi_mls(2, jp_var_2083(jl_var_2168, lay_var_2112) + 1) * coldry_var_2090(jl_var_2168, lay_var_2112) * 1D-20
      ELSE
        adjcolco2_var_2121 = colco2_var_2089(jl_var_2168, lay_var_2112)
      END IF
      ind0_var_2105 = ((jp_var_2083(jl_var_2168, lay_var_2112) - 13) * 5 + (jt_var_2084(jl_var_2168, lay_var_2112) - 1)) * nspb_var_217(7) + 1
      ind1_var_2106 = ((jp_var_2083(jl_var_2168, lay_var_2112) - 12) * 5 + (jt1_var_2085(jl_var_2168, lay_var_2112) - 1)) * nspb_var_217(7) + 1
      indm_var_2109 = indminor_var_2100(jl_var_2168, lay_var_2112)
      DO ig_var_2110 = 1, 12
        absco2_var_2145 = kb_mco2_var_192(indm_var_2109, ig_var_2110) + minorfrac_var_2099(jl_var_2168, lay_var_2112) * (kb_mco2_var_192(indm_var_2109 + 1, ig_var_2110) - kb_mco2_var_192(indm_var_2109, ig_var_2110))
        taug_var_2077(jl_var_2168, 76 + ig_var_2110, lay_var_2112) = colo3_var_2088(jl_var_2168, lay_var_2112) * (fac00_var_2079(jl_var_2168, lay_var_2112) * absb_var_189(ind0_var_2105, ig_var_2110) + fac10_var_2081(jl_var_2168, lay_var_2112) * absb_var_189(ind0_var_2105 + 1, ig_var_2110) + fac01_var_2080(jl_var_2168, lay_var_2112) * absb_var_189(ind1_var_2106, ig_var_2110) + fac11_var_2082(jl_var_2168, lay_var_2112) * absb_var_189(ind1_var_2106 + 1, ig_var_2110)) + adjcolco2_var_2121 * absco2_var_2145
        fracs_var_2095(jl_var_2168, 76 + ig_var_2110, lay_var_2112) = fracrefb_var_187(ig_var_2110)
      END DO
    END DO
  END DO
  DO lay_var_2112 = laytrop_max_var_2159 + 1, klev_var_2076
    DO jl_var_2168 = kidia_var_2074, kfdia_var_2075
      taug_var_2077(jl_var_2168, 82, lay_var_2112) = taug_var_2077(jl_var_2168, 82, lay_var_2112) * 0.92D0
      taug_var_2077(jl_var_2168, 83, lay_var_2112) = taug_var_2077(jl_var_2168, 83, lay_var_2112) * 0.88D0
      taug_var_2077(jl_var_2168, 84, lay_var_2112) = taug_var_2077(jl_var_2168, 84, lay_var_2112) * 1.07D0
      taug_var_2077(jl_var_2168, 85, lay_var_2112) = taug_var_2077(jl_var_2168, 85, lay_var_2112) * 1.1D0
      taug_var_2077(jl_var_2168, 86, lay_var_2112) = taug_var_2077(jl_var_2168, 86, lay_var_2112) * 0.99D0
      taug_var_2077(jl_var_2168, 87, lay_var_2112) = taug_var_2077(jl_var_2168, 87, lay_var_2112) * 0.855D0
    END DO
  END DO
  IF (laytrop_max_var_2159 /= laytrop_min_var_2158) THEN
    DO lay_var_2112 = laytrop_min_var_2158 + 1, laytrop_max_var_2159
      ixc0_var_2165 = ixc_var_2160(lay_var_2112)
      DO ixp_var_2166 = 1, ixc0_var_2165
        jl_var_2168 = ixlow_var_2161(ixp_var_2166, lay_var_2112)
        speccomb_var_2101 = colh2o_var_2087(jl_var_2168, lay_var_2112) + rat_h2oo3(jl_var_2168, lay_var_2112) * colo3_var_2088(jl_var_2168, lay_var_2112)
        specparm_var_2148 = MIN(colh2o_var_2087(jl_var_2168, lay_var_2112) / speccomb_var_2101, oneminus_var_2086)
        specmult_var_2147 = 8.0D0 * (specparm_var_2148)
        js_var_2111 = 1 + INT(specmult_var_2147)
        fs_var_2146 = ((specmult_var_2147) - AINT((specmult_var_2147)))
        speccomb1_var_2102 = colh2o_var_2087(jl_var_2168, lay_var_2112) + rat_h2oo3_1(jl_var_2168, lay_var_2112) * colo3_var_2088(jl_var_2168, lay_var_2112)
        specparm1_var_2151 = MIN(colh2o_var_2087(jl_var_2168, lay_var_2112) / speccomb1_var_2102, oneminus_var_2086)
        specmult1_var_2150 = 8.0D0 * (specparm1_var_2151)
        js1_var_2113 = 1 + INT(specmult1_var_2150)
        fs1_var_2149 = ((specmult1_var_2150) - AINT((specmult1_var_2150)))
        speccomb_mco2_var_2103 = colh2o_var_2087(jl_var_2168, lay_var_2112) + refrat_m_a_var_2117 * colo3_var_2088(jl_var_2168, lay_var_2112)
        specparm_mco2_var_2157 = MIN(colh2o_var_2087(jl_var_2168, lay_var_2112) / speccomb_mco2_var_2103, oneminus_var_2086)
        specmult_mco2_var_2156 = 8.0D0 * specparm_mco2_var_2157
        jmco2_var_2115 = 1 + INT(specmult_mco2_var_2156)
        fmco2_var_2155 = ((specmult_mco2_var_2156) - AINT((specmult_mco2_var_2156)))
        chi_co2_var_2118 = colco2_var_2089(jl_var_2168, lay_var_2112) / (coldry_var_2090(jl_var_2168, lay_var_2112))
        ratco2_var_2119 = 1D+20 * chi_co2_var_2118 / chi_mls(2, jp_var_2083(jl_var_2168, lay_var_2112) + 1)
        IF (ratco2_var_2119 .GT. 3.0D0) THEN
          adjfac_var_2120 = 3.0D0 + (ratco2_var_2119 - 3.0D0) ** 0.79D0
          adjcolco2_var_2121 = adjfac_var_2120 * chi_mls(2, jp_var_2083(jl_var_2168, lay_var_2112) + 1) * coldry_var_2090(jl_var_2168, lay_var_2112) * 1D-20
        ELSE
          adjcolco2_var_2121 = colco2_var_2089(jl_var_2168, lay_var_2112)
        END IF
        speccomb_planck_var_2104 = colh2o_var_2087(jl_var_2168, lay_var_2112) + refrat_planck_a_var_2116 * colo3_var_2088(jl_var_2168, lay_var_2112)
        specparm_planck_var_2154 = MIN(colh2o_var_2087(jl_var_2168, lay_var_2112) / speccomb_planck_var_2104, oneminus_var_2086)
        specmult_planck_var_2153 = 8.0D0 * specparm_planck_var_2154
        jpl_var_2114 = 1 + INT(specmult_planck_var_2153)
        fpl_var_2152 = ((specmult_planck_var_2153) - AINT((specmult_planck_var_2153)))
        ind0_var_2105 = ((jp_var_2083(jl_var_2168, lay_var_2112) - 1) * 5 + (jt_var_2084(jl_var_2168, lay_var_2112) - 1)) * nspa_var_216(7) + js_var_2111
        ind1_var_2106 = (jp_var_2083(jl_var_2168, lay_var_2112) * 5 + (jt1_var_2085(jl_var_2168, lay_var_2112) - 1)) * nspa_var_216(7) + js1_var_2113
        inds_var_2107 = indself_var_2094(jl_var_2168, lay_var_2112)
        indf_var_2108 = indfor_var_2096(jl_var_2168, lay_var_2112)
        indm_var_2109 = indminor_var_2100(jl_var_2168, lay_var_2112)
        IF (specparm_var_2148 .LT. 0.125D0) THEN
          p_var_2134 = fs_var_2146 - 1.0D0
          p4_var_2135 = p_var_2134 ** 4
          fk0_var_2136 = p4_var_2135
          fk1_var_2137 = 1.0D0 - p_var_2134 - 2.0D0 * p4_var_2135
          fk2_var_2138 = p_var_2134 + p4_var_2135
          fac000_var_2122 = fk0_var_2136 * fac00_var_2079(jl_var_2168, lay_var_2112)
          fac100_var_2123 = fk1_var_2137 * fac00_var_2079(jl_var_2168, lay_var_2112)
          fac200_var_2124 = fk2_var_2138 * fac00_var_2079(jl_var_2168, lay_var_2112)
          fac010_var_2125 = fk0_var_2136 * fac10_var_2081(jl_var_2168, lay_var_2112)
          fac110_var_2126 = fk1_var_2137 * fac10_var_2081(jl_var_2168, lay_var_2112)
          fac210_var_2127 = fk2_var_2138 * fac10_var_2081(jl_var_2168, lay_var_2112)
        ELSE IF (specparm_var_2148 .GT. 0.875D0) THEN
          p_var_2134 = - fs_var_2146
          p4_var_2135 = p_var_2134 ** 4
          fk0_var_2136 = p4_var_2135
          fk1_var_2137 = 1.0D0 - p_var_2134 - 2.0D0 * p4_var_2135
          fk2_var_2138 = p_var_2134 + p4_var_2135
          fac000_var_2122 = fk0_var_2136 * fac00_var_2079(jl_var_2168, lay_var_2112)
          fac100_var_2123 = fk1_var_2137 * fac00_var_2079(jl_var_2168, lay_var_2112)
          fac200_var_2124 = fk2_var_2138 * fac00_var_2079(jl_var_2168, lay_var_2112)
          fac010_var_2125 = fk0_var_2136 * fac10_var_2081(jl_var_2168, lay_var_2112)
          fac110_var_2126 = fk1_var_2137 * fac10_var_2081(jl_var_2168, lay_var_2112)
          fac210_var_2127 = fk2_var_2138 * fac10_var_2081(jl_var_2168, lay_var_2112)
        ELSE
          fac000_var_2122 = (1.0D0 - fs_var_2146) * fac00_var_2079(jl_var_2168, lay_var_2112)
          fac010_var_2125 = (1.0D0 - fs_var_2146) * fac10_var_2081(jl_var_2168, lay_var_2112)
          fac100_var_2123 = fs_var_2146 * fac00_var_2079(jl_var_2168, lay_var_2112)
          fac110_var_2126 = fs_var_2146 * fac10_var_2081(jl_var_2168, lay_var_2112)
          fac200_var_2124 = 0.0D0
          fac210_var_2127 = 0.0D0
        END IF
        IF (specparm1_var_2151 .LT. 0.125D0) THEN
          p_var_2134 = fs1_var_2149 - 1.0D0
          p4_var_2135 = p_var_2134 ** 4
          fk0_var_2136 = p4_var_2135
          fk1_var_2137 = 1.0D0 - p_var_2134 - 2.0D0 * p4_var_2135
          fk2_var_2138 = p_var_2134 + p4_var_2135
          fac001_var_2128 = fk0_var_2136 * fac01_var_2080(jl_var_2168, lay_var_2112)
          fac101_var_2129 = fk1_var_2137 * fac01_var_2080(jl_var_2168, lay_var_2112)
          fac201_var_2130 = fk2_var_2138 * fac01_var_2080(jl_var_2168, lay_var_2112)
          fac011_var_2131 = fk0_var_2136 * fac11_var_2082(jl_var_2168, lay_var_2112)
          fac111_var_2132 = fk1_var_2137 * fac11_var_2082(jl_var_2168, lay_var_2112)
          fac211_var_2133 = fk2_var_2138 * fac11_var_2082(jl_var_2168, lay_var_2112)
        ELSE IF (specparm1_var_2151 .GT. 0.875D0) THEN
          p_var_2134 = - fs1_var_2149
          p4_var_2135 = p_var_2134 ** 4
          fk0_var_2136 = p4_var_2135
          fk1_var_2137 = 1.0D0 - p_var_2134 - 2.0D0 * p4_var_2135
          fk2_var_2138 = p_var_2134 + p4_var_2135
          fac001_var_2128 = fk0_var_2136 * fac01_var_2080(jl_var_2168, lay_var_2112)
          fac101_var_2129 = fk1_var_2137 * fac01_var_2080(jl_var_2168, lay_var_2112)
          fac201_var_2130 = fk2_var_2138 * fac01_var_2080(jl_var_2168, lay_var_2112)
          fac011_var_2131 = fk0_var_2136 * fac11_var_2082(jl_var_2168, lay_var_2112)
          fac111_var_2132 = fk1_var_2137 * fac11_var_2082(jl_var_2168, lay_var_2112)
          fac211_var_2133 = fk2_var_2138 * fac11_var_2082(jl_var_2168, lay_var_2112)
        ELSE
          fac001_var_2128 = (1.0D0 - fs1_var_2149) * fac01_var_2080(jl_var_2168, lay_var_2112)
          fac011_var_2131 = (1.0D0 - fs1_var_2149) * fac11_var_2082(jl_var_2168, lay_var_2112)
          fac101_var_2129 = fs1_var_2149 * fac01_var_2080(jl_var_2168, lay_var_2112)
          fac111_var_2132 = fs1_var_2149 * fac11_var_2082(jl_var_2168, lay_var_2112)
          fac201_var_2130 = 0.0D0
          fac211_var_2133 = 0.0D0
        END IF
        IF (specparm_var_2148 .LT. 0.125D0) THEN
          tau_major_var_2141(1 : ng7) = speccomb_var_2101 * (fac000_var_2122 * absa_var_188(ind0_var_2105, 1 : 12) + fac100_var_2123 * absa_var_188(ind0_var_2105 + 1, 1 : 12) + fac200_var_2124 * absa_var_188(ind0_var_2105 + 2, 1 : 12) + fac010_var_2125 * absa_var_188(ind0_var_2105 + 9, 1 : 12) + fac110_var_2126 * absa_var_188(ind0_var_2105 + 10, 1 : 12) + fac210_var_2127 * absa_var_188(ind0_var_2105 + 11, 1 : 12))
        ELSE IF (specparm_var_2148 .GT. 0.875D0) THEN
          tau_major_var_2141(1 : ng7) = speccomb_var_2101 * (fac200_var_2124 * absa_var_188(ind0_var_2105 - 1, 1 : 12) + fac100_var_2123 * absa_var_188(ind0_var_2105, 1 : 12) + fac000_var_2122 * absa_var_188(ind0_var_2105 + 1, 1 : 12) + fac210_var_2127 * absa_var_188(ind0_var_2105 + 8, 1 : 12) + fac110_var_2126 * absa_var_188(ind0_var_2105 + 9, 1 : 12) + fac010_var_2125 * absa_var_188(ind0_var_2105 + 10, 1 : 12))
        ELSE
          tau_major_var_2141(1 : ng7) = speccomb_var_2101 * (fac000_var_2122 * absa_var_188(ind0_var_2105, 1 : 12) + fac100_var_2123 * absa_var_188(ind0_var_2105 + 1, 1 : 12) + fac010_var_2125 * absa_var_188(ind0_var_2105 + 9, 1 : 12) + fac110_var_2126 * absa_var_188(ind0_var_2105 + 10, 1 : 12))
        END IF
        IF (specparm1_var_2151 .LT. 0.125D0) THEN
          tau_major1_var_2142(1 : ng7) = speccomb1_var_2102 * (fac001_var_2128 * absa_var_188(ind1_var_2106, 1 : 12) + fac101_var_2129 * absa_var_188(ind1_var_2106 + 1, 1 : 12) + fac201_var_2130 * absa_var_188(ind1_var_2106 + 2, 1 : 12) + fac011_var_2131 * absa_var_188(ind1_var_2106 + 9, 1 : 12) + fac111_var_2132 * absa_var_188(ind1_var_2106 + 10, 1 : 12) + fac211_var_2133 * absa_var_188(ind1_var_2106 + 11, 1 : 12))
        ELSE IF (specparm1_var_2151 .GT. 0.875D0) THEN
          tau_major1_var_2142(1 : ng7) = speccomb1_var_2102 * (fac201_var_2130 * absa_var_188(ind1_var_2106 - 1, 1 : 12) + fac101_var_2129 * absa_var_188(ind1_var_2106, 1 : 12) + fac001_var_2128 * absa_var_188(ind1_var_2106 + 1, 1 : 12) + fac211_var_2133 * absa_var_188(ind1_var_2106 + 8, 1 : 12) + fac111_var_2132 * absa_var_188(ind1_var_2106 + 9, 1 : 12) + fac011_var_2131 * absa_var_188(ind1_var_2106 + 10, 1 : 12))
        ELSE
          tau_major1_var_2142(1 : ng7) = speccomb1_var_2102 * (fac001_var_2128 * absa_var_188(ind1_var_2106, 1 : 12) + fac101_var_2129 * absa_var_188(ind1_var_2106 + 1, 1 : 12) + fac011_var_2131 * absa_var_188(ind1_var_2106 + 9, 1 : 12) + fac111_var_2132 * absa_var_188(ind1_var_2106 + 10, 1 : 12))
        END IF
        DO ig_var_2110 = 1, 12
          tauself_var_2140 = selffac_var_2092(jl_var_2168, lay_var_2112) * (selfref_var_190(inds_var_2107, ig_var_2110) + selffrac_var_2093(jl_var_2168, lay_var_2112) * (selfref_var_190(inds_var_2107 + 1, ig_var_2110) - selfref_var_190(inds_var_2107, ig_var_2110)))
          taufor_var_2139 = forfac_var_2098(jl_var_2168, lay_var_2112) * (forref_var_193(indf_var_2108, ig_var_2110) + forfrac_var_2097(jl_var_2168, lay_var_2112) * (forref_var_193(indf_var_2108 + 1, ig_var_2110) - forref_var_193(indf_var_2108, ig_var_2110)))
          co2m1_var_2143 = ka_mco2_var_191(jmco2_var_2115, indm_var_2109, ig_var_2110) + fmco2_var_2155 * (ka_mco2_var_191(jmco2_var_2115 + 1, indm_var_2109, ig_var_2110) - ka_mco2_var_191(jmco2_var_2115, indm_var_2109, ig_var_2110))
          co2m2_var_2144 = ka_mco2_var_191(jmco2_var_2115, indm_var_2109 + 1, ig_var_2110) + fmco2_var_2155 * (ka_mco2_var_191(jmco2_var_2115 + 1, indm_var_2109 + 1, ig_var_2110) - ka_mco2_var_191(jmco2_var_2115, indm_var_2109 + 1, ig_var_2110))
          absco2_var_2145 = co2m1_var_2143 + minorfrac_var_2099(jl_var_2168, lay_var_2112) * (co2m2_var_2144 - co2m1_var_2143)
          taug_var_2077(jl_var_2168, 76 + ig_var_2110, lay_var_2112) = tau_major_var_2141(ig_var_2110) + tau_major1_var_2142(ig_var_2110) + tauself_var_2140 + taufor_var_2139 + adjcolco2_var_2121 * absco2_var_2145
          fracs_var_2095(jl_var_2168, 76 + ig_var_2110, lay_var_2112) = fracrefa_var_186(ig_var_2110, jpl_var_2114) + fpl_var_2152 * (fracrefa_var_186(ig_var_2110, jpl_var_2114 + 1) - fracrefa_var_186(ig_var_2110, jpl_var_2114))
        END DO
      END DO
      ixc0_var_2165 = kfdia_var_2075 - kidia_var_2074 + 1 - ixc0_var_2165
      DO ixp_var_2166 = 1, ixc0_var_2165
        jl_var_2168 = ixhigh_var_2162(ixp_var_2166, lay_var_2112)
        chi_co2_var_2118 = colco2_var_2089(jl_var_2168, lay_var_2112) / (coldry_var_2090(jl_var_2168, lay_var_2112))
        ratco2_var_2119 = 1D+20 * chi_co2_var_2118 / chi_mls(2, jp_var_2083(jl_var_2168, lay_var_2112) + 1)
        IF (ratco2_var_2119 .GT. 3.0D0) THEN
          adjfac_var_2120 = 2.0D0 + (ratco2_var_2119 - 2.0D0) ** 0.79D0
          adjcolco2_var_2121 = adjfac_var_2120 * chi_mls(2, jp_var_2083(jl_var_2168, lay_var_2112) + 1) * coldry_var_2090(jl_var_2168, lay_var_2112) * 1D-20
        ELSE
          adjcolco2_var_2121 = colco2_var_2089(jl_var_2168, lay_var_2112)
        END IF
        ind0_var_2105 = ((jp_var_2083(jl_var_2168, lay_var_2112) - 13) * 5 + (jt_var_2084(jl_var_2168, lay_var_2112) - 1)) * nspb_var_217(7) + 1
        ind1_var_2106 = ((jp_var_2083(jl_var_2168, lay_var_2112) - 12) * 5 + (jt1_var_2085(jl_var_2168, lay_var_2112) - 1)) * nspb_var_217(7) + 1
        indm_var_2109 = indminor_var_2100(jl_var_2168, lay_var_2112)
        DO ig_var_2110 = 1, 12
          absco2_var_2145 = kb_mco2_var_192(indm_var_2109, ig_var_2110) + minorfrac_var_2099(jl_var_2168, lay_var_2112) * (kb_mco2_var_192(indm_var_2109 + 1, ig_var_2110) - kb_mco2_var_192(indm_var_2109, ig_var_2110))
          taug_var_2077(jl_var_2168, 76 + ig_var_2110, lay_var_2112) = colo3_var_2088(jl_var_2168, lay_var_2112) * (fac00_var_2079(jl_var_2168, lay_var_2112) * absb_var_189(ind0_var_2105, ig_var_2110) + fac10_var_2081(jl_var_2168, lay_var_2112) * absb_var_189(ind0_var_2105 + 1, ig_var_2110) + fac01_var_2080(jl_var_2168, lay_var_2112) * absb_var_189(ind1_var_2106, ig_var_2110) + fac11_var_2082(jl_var_2168, lay_var_2112) * absb_var_189(ind1_var_2106 + 1, ig_var_2110)) + adjcolco2_var_2121 * absco2_var_2145
          fracs_var_2095(jl_var_2168, 76 + ig_var_2110, lay_var_2112) = fracrefb_var_187(ig_var_2110)
        END DO
      END DO
      DO ixp_var_2166 = 1, ixc0_var_2165
        jl_var_2168 = ixhigh_var_2162(ixp_var_2166, lay_var_2112)
        taug_var_2077(jl_var_2168, 82, lay_var_2112) = taug_var_2077(jl_var_2168, 82, lay_var_2112) * 0.92D0
        taug_var_2077(jl_var_2168, 83, lay_var_2112) = taug_var_2077(jl_var_2168, 83, lay_var_2112) * 0.88D0
        taug_var_2077(jl_var_2168, 84, lay_var_2112) = taug_var_2077(jl_var_2168, 84, lay_var_2112) * 1.07D0
        taug_var_2077(jl_var_2168, 85, lay_var_2112) = taug_var_2077(jl_var_2168, 85, lay_var_2112) * 1.1D0
        taug_var_2077(jl_var_2168, 86, lay_var_2112) = taug_var_2077(jl_var_2168, 86, lay_var_2112) * 0.99D0
        taug_var_2077(jl_var_2168, 87, lay_var_2112) = taug_var_2077(jl_var_2168, 87, lay_var_2112) * 0.855D0
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol7
SUBROUTINE rrtm_taumol11(kidia_var_2169, kfdia_var_2170, klev_var_2171, taug_var_2172, p_tauaerl_var_2173, fac00_var_2174, fac01_var_2175, fac10_var_2176, fac11_var_2177, forfac_var_2188, forfrac_var_2189, indfor_var_2187, jp_var_2178, jt_var_2179, jt1_var_2180, colh2o_var_2181, colo2, laytrop_var_2182, selffac_var_2183, selffrac_var_2184, indself_var_2185, fracs_var_2186, minorfrac_var_2190, indminor_var_2191, scaleminor_var_2192)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta11, ONLY: absa_var_122, absb_var_123, forref_var_125, fracrefa_var_120, fracrefb_var_121, ka_mo2, kb_mo2, selfref_var_124
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2169
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2170
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2171
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2172(kidia_var_2169 : kfdia_var_2170, 140, klev_var_2171)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2173(kidia_var_2169 : kfdia_var_2170, klev_var_2171, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2174(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2175(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2176(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2177(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2178(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2179(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2180(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2181(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  REAL(KIND = 8), INTENT(IN) :: colo2(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2182(kidia_var_2169 : kfdia_var_2170)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2183(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2184(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2185(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2186(kidia_var_2169 : kfdia_var_2170, 140, klev_var_2171)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2187(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2188(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2189(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2190(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2191(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  REAL(KIND = 8), INTENT(IN) :: scaleminor_var_2192(kidia_var_2169 : kfdia_var_2170, klev_var_2171)
  INTEGER(KIND = 4) :: ind0_var_2193, ind1_var_2194
  INTEGER(KIND = 4) :: inds_var_2195, indf_var_2196, indm_var_2197
  INTEGER(KIND = 4) :: ig_var_2198, lay_var_2199
  REAL(KIND = 8) :: taufor_var_2200, tauself_var_2201, scaleo2, tauo2
  INTEGER(KIND = 4) :: laytrop_min_var_2202, laytrop_max_var_2203
  INTEGER(KIND = 4) :: ixc_var_2204(klev_var_2171), ixlow_var_2205(kfdia_var_2170, klev_var_2171), ixhigh_var_2206(kfdia_var_2170, klev_var_2171)
  INTEGER(KIND = 4) :: ich_var_2207, icl_var_2208, ixc0_var_2209, ixp_var_2210, jc_var_2211, jl_var_2212
  laytrop_min_var_2202 = MINVAL(laytrop_var_2182)
  laytrop_max_var_2203 = MAXVAL(laytrop_var_2182)
  ixlow_var_2205 = 0
  ixhigh_var_2206 = 0
  ixc_var_2204 = 0
  DO lay_var_2199 = laytrop_min_var_2202 + 1, laytrop_max_var_2203
    icl_var_2208 = 0
    ich_var_2207 = 0
    DO jc_var_2211 = kidia_var_2169, kfdia_var_2170
      IF (lay_var_2199 <= laytrop_var_2182(jc_var_2211)) THEN
        icl_var_2208 = icl_var_2208 + 1
        ixlow_var_2205(icl_var_2208, lay_var_2199) = jc_var_2211
      ELSE
        ich_var_2207 = ich_var_2207 + 1
        ixhigh_var_2206(ich_var_2207, lay_var_2199) = jc_var_2211
      END IF
    END DO
    ixc_var_2204(lay_var_2199) = icl_var_2208
  END DO
  DO lay_var_2199 = 1, laytrop_min_var_2202
    DO jl_var_2212 = kidia_var_2169, kfdia_var_2170
      ind0_var_2193 = ((jp_var_2178(jl_var_2212, lay_var_2199) - 1) * 5 + (jt_var_2179(jl_var_2212, lay_var_2199) - 1)) * nspa_var_216(11) + 1
      ind1_var_2194 = (jp_var_2178(jl_var_2212, lay_var_2199) * 5 + (jt1_var_2180(jl_var_2212, lay_var_2199) - 1)) * nspa_var_216(11) + 1
      inds_var_2195 = indself_var_2185(jl_var_2212, lay_var_2199)
      indf_var_2196 = indfor_var_2187(jl_var_2212, lay_var_2199)
      indm_var_2197 = indminor_var_2191(jl_var_2212, lay_var_2199)
      scaleo2 = colo2(jl_var_2212, lay_var_2199) * scaleminor_var_2192(jl_var_2212, lay_var_2199)
      DO ig_var_2198 = 1, 8
        tauself_var_2201 = selffac_var_2183(jl_var_2212, lay_var_2199) * (selfref_var_124(inds_var_2195, ig_var_2198) + selffrac_var_2184(jl_var_2212, lay_var_2199) * (selfref_var_124(inds_var_2195 + 1, ig_var_2198) - selfref_var_124(inds_var_2195, ig_var_2198)))
        taufor_var_2200 = forfac_var_2188(jl_var_2212, lay_var_2199) * (forref_var_125(indf_var_2196, ig_var_2198) + forfrac_var_2189(jl_var_2212, lay_var_2199) * (forref_var_125(indf_var_2196 + 1, ig_var_2198) - forref_var_125(indf_var_2196, ig_var_2198)))
        tauo2 = scaleo2 * (ka_mo2(indm_var_2197, ig_var_2198) + minorfrac_var_2190(jl_var_2212, lay_var_2199) * (ka_mo2(indm_var_2197 + 1, ig_var_2198) - ka_mo2(indm_var_2197, ig_var_2198)))
        taug_var_2172(jl_var_2212, 114 + ig_var_2198, lay_var_2199) = colh2o_var_2181(jl_var_2212, lay_var_2199) * (fac00_var_2174(jl_var_2212, lay_var_2199) * absa_var_122(ind0_var_2193, ig_var_2198) + fac10_var_2176(jl_var_2212, lay_var_2199) * absa_var_122(ind0_var_2193 + 1, ig_var_2198) + fac01_var_2175(jl_var_2212, lay_var_2199) * absa_var_122(ind1_var_2194, ig_var_2198) + fac11_var_2177(jl_var_2212, lay_var_2199) * absa_var_122(ind1_var_2194 + 1, ig_var_2198)) + tauself_var_2201 + taufor_var_2200 + tauo2
        fracs_var_2186(jl_var_2212, 114 + ig_var_2198, lay_var_2199) = fracrefa_var_120(ig_var_2198)
      END DO
    END DO
  END DO
  DO lay_var_2199 = laytrop_max_var_2203 + 1, klev_var_2171
    DO jl_var_2212 = kidia_var_2169, kfdia_var_2170
      ind0_var_2193 = ((jp_var_2178(jl_var_2212, lay_var_2199) - 13) * 5 + (jt_var_2179(jl_var_2212, lay_var_2199) - 1)) * nspb_var_217(11) + 1
      ind1_var_2194 = ((jp_var_2178(jl_var_2212, lay_var_2199) - 12) * 5 + (jt1_var_2180(jl_var_2212, lay_var_2199) - 1)) * nspb_var_217(11) + 1
      indf_var_2196 = indfor_var_2187(jl_var_2212, lay_var_2199)
      indm_var_2197 = indminor_var_2191(jl_var_2212, lay_var_2199)
      scaleo2 = colo2(jl_var_2212, lay_var_2199) * scaleminor_var_2192(jl_var_2212, lay_var_2199)
      DO ig_var_2198 = 1, 8
        taufor_var_2200 = forfac_var_2188(jl_var_2212, lay_var_2199) * (forref_var_125(indf_var_2196, ig_var_2198) + forfrac_var_2189(jl_var_2212, lay_var_2199) * (forref_var_125(indf_var_2196 + 1, ig_var_2198) - forref_var_125(indf_var_2196, ig_var_2198)))
        tauo2 = scaleo2 * (kb_mo2(indm_var_2197, ig_var_2198) + minorfrac_var_2190(jl_var_2212, lay_var_2199) * (kb_mo2(indm_var_2197 + 1, ig_var_2198) - kb_mo2(indm_var_2197, ig_var_2198)))
        taug_var_2172(jl_var_2212, 114 + ig_var_2198, lay_var_2199) = colh2o_var_2181(jl_var_2212, lay_var_2199) * (fac00_var_2174(jl_var_2212, lay_var_2199) * absb_var_123(ind0_var_2193, ig_var_2198) + fac10_var_2176(jl_var_2212, lay_var_2199) * absb_var_123(ind0_var_2193 + 1, ig_var_2198) + fac01_var_2175(jl_var_2212, lay_var_2199) * absb_var_123(ind1_var_2194, ig_var_2198) + fac11_var_2177(jl_var_2212, lay_var_2199) * absb_var_123(ind1_var_2194 + 1, ig_var_2198)) + taufor_var_2200 + tauo2
        fracs_var_2186(jl_var_2212, 114 + ig_var_2198, lay_var_2199) = fracrefb_var_121(ig_var_2198)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2203 /= laytrop_min_var_2202) THEN
    DO lay_var_2199 = laytrop_min_var_2202 + 1, laytrop_max_var_2203
      ixc0_var_2209 = ixc_var_2204(lay_var_2199)
      DO ixp_var_2210 = 1, ixc0_var_2209
        jl_var_2212 = ixlow_var_2205(ixp_var_2210, lay_var_2199)
        ind0_var_2193 = ((jp_var_2178(jl_var_2212, lay_var_2199) - 1) * 5 + (jt_var_2179(jl_var_2212, lay_var_2199) - 1)) * nspa_var_216(11) + 1
        ind1_var_2194 = (jp_var_2178(jl_var_2212, lay_var_2199) * 5 + (jt1_var_2180(jl_var_2212, lay_var_2199) - 1)) * nspa_var_216(11) + 1
        inds_var_2195 = indself_var_2185(jl_var_2212, lay_var_2199)
        indf_var_2196 = indfor_var_2187(jl_var_2212, lay_var_2199)
        indm_var_2197 = indminor_var_2191(jl_var_2212, lay_var_2199)
        scaleo2 = colo2(jl_var_2212, lay_var_2199) * scaleminor_var_2192(jl_var_2212, lay_var_2199)
        DO ig_var_2198 = 1, 8
          tauself_var_2201 = selffac_var_2183(jl_var_2212, lay_var_2199) * (selfref_var_124(inds_var_2195, ig_var_2198) + selffrac_var_2184(jl_var_2212, lay_var_2199) * (selfref_var_124(inds_var_2195 + 1, ig_var_2198) - selfref_var_124(inds_var_2195, ig_var_2198)))
          taufor_var_2200 = forfac_var_2188(jl_var_2212, lay_var_2199) * (forref_var_125(indf_var_2196, ig_var_2198) + forfrac_var_2189(jl_var_2212, lay_var_2199) * (forref_var_125(indf_var_2196 + 1, ig_var_2198) - forref_var_125(indf_var_2196, ig_var_2198)))
          tauo2 = scaleo2 * (ka_mo2(indm_var_2197, ig_var_2198) + minorfrac_var_2190(jl_var_2212, lay_var_2199) * (ka_mo2(indm_var_2197 + 1, ig_var_2198) - ka_mo2(indm_var_2197, ig_var_2198)))
          taug_var_2172(jl_var_2212, 114 + ig_var_2198, lay_var_2199) = colh2o_var_2181(jl_var_2212, lay_var_2199) * (fac00_var_2174(jl_var_2212, lay_var_2199) * absa_var_122(ind0_var_2193, ig_var_2198) + fac10_var_2176(jl_var_2212, lay_var_2199) * absa_var_122(ind0_var_2193 + 1, ig_var_2198) + fac01_var_2175(jl_var_2212, lay_var_2199) * absa_var_122(ind1_var_2194, ig_var_2198) + fac11_var_2177(jl_var_2212, lay_var_2199) * absa_var_122(ind1_var_2194 + 1, ig_var_2198)) + tauself_var_2201 + taufor_var_2200 + tauo2
          fracs_var_2186(jl_var_2212, 114 + ig_var_2198, lay_var_2199) = fracrefa_var_120(ig_var_2198)
        END DO
      END DO
      ixc0_var_2209 = kfdia_var_2170 - kidia_var_2169 + 1 - ixc0_var_2209
      DO ixp_var_2210 = 1, ixc0_var_2209
        jl_var_2212 = ixhigh_var_2206(ixp_var_2210, lay_var_2199)
        ind0_var_2193 = ((jp_var_2178(jl_var_2212, lay_var_2199) - 13) * 5 + (jt_var_2179(jl_var_2212, lay_var_2199) - 1)) * nspb_var_217(11) + 1
        ind1_var_2194 = ((jp_var_2178(jl_var_2212, lay_var_2199) - 12) * 5 + (jt1_var_2180(jl_var_2212, lay_var_2199) - 1)) * nspb_var_217(11) + 1
        indf_var_2196 = indfor_var_2187(jl_var_2212, lay_var_2199)
        indm_var_2197 = indminor_var_2191(jl_var_2212, lay_var_2199)
        scaleo2 = colo2(jl_var_2212, lay_var_2199) * scaleminor_var_2192(jl_var_2212, lay_var_2199)
        DO ig_var_2198 = 1, 8
          taufor_var_2200 = forfac_var_2188(jl_var_2212, lay_var_2199) * (forref_var_125(indf_var_2196, ig_var_2198) + forfrac_var_2189(jl_var_2212, lay_var_2199) * (forref_var_125(indf_var_2196 + 1, ig_var_2198) - forref_var_125(indf_var_2196, ig_var_2198)))
          tauo2 = scaleo2 * (kb_mo2(indm_var_2197, ig_var_2198) + minorfrac_var_2190(jl_var_2212, lay_var_2199) * (kb_mo2(indm_var_2197 + 1, ig_var_2198) - kb_mo2(indm_var_2197, ig_var_2198)))
          taug_var_2172(jl_var_2212, 114 + ig_var_2198, lay_var_2199) = colh2o_var_2181(jl_var_2212, lay_var_2199) * (fac00_var_2174(jl_var_2212, lay_var_2199) * absb_var_123(ind0_var_2193, ig_var_2198) + fac10_var_2176(jl_var_2212, lay_var_2199) * absb_var_123(ind0_var_2193 + 1, ig_var_2198) + fac01_var_2175(jl_var_2212, lay_var_2199) * absb_var_123(ind1_var_2194, ig_var_2198) + fac11_var_2177(jl_var_2212, lay_var_2199) * absb_var_123(ind1_var_2194 + 1, ig_var_2198)) + taufor_var_2200 + tauo2
          fracs_var_2186(jl_var_2212, 114 + ig_var_2198, lay_var_2199) = fracrefb_var_121(ig_var_2198)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol11
SUBROUTINE rrtm_taumol10(kidia_var_2213, kfdia_var_2214, klev_var_2215, taug_var_2216, p_tauaerl_var_2217, fac00_var_2218, fac01_var_2219, fac10_var_2220, fac11_var_2221, forfac_var_2233, forfrac_var_2232, indfor_var_2231, jp_var_2222, jt_var_2223, jt1_var_2224, colh2o_var_2225, laytrop_var_2226, selffac_var_2228, selffrac_var_2229, indself_var_2230, fracs_var_2227)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta10, ONLY: absa_var_116, absb_var_117, forref_var_119, fracrefa_var_114, fracrefb_var_115, selfref_var_118
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2213
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2214
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2215
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2216(kidia_var_2213 : kfdia_var_2214, 140, klev_var_2215)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2217(kidia_var_2213 : kfdia_var_2214, klev_var_2215, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2218(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2219(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2220(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2221(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2222(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2223(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2224(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2225(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2226(kidia_var_2213 : kfdia_var_2214)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2227(kidia_var_2213 : kfdia_var_2214, 140, klev_var_2215)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2228(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2229(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2230(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2231(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2232(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2233(kidia_var_2213 : kfdia_var_2214, klev_var_2215)
  INTEGER(KIND = 4) :: ind0_var_2234, ind1_var_2235
  INTEGER(KIND = 4) :: inds_var_2236, indf_var_2237
  INTEGER(KIND = 4) :: ig_var_2238, lay_var_2239
  REAL(KIND = 8) :: taufor_var_2240, tauself_var_2241
  INTEGER(KIND = 4) :: laytrop_min_var_2242, laytrop_max_var_2243
  INTEGER(KIND = 4) :: ixc_var_2244(klev_var_2215), ixlow_var_2245(kfdia_var_2214, klev_var_2215), ixhigh_var_2246(kfdia_var_2214, klev_var_2215)
  INTEGER(KIND = 4) :: ich_var_2247, icl_var_2248, ixc0_var_2249, ixp_var_2250, jc_var_2251, jl_var_2252
  laytrop_min_var_2242 = MINVAL(laytrop_var_2226)
  laytrop_max_var_2243 = MAXVAL(laytrop_var_2226)
  ixlow_var_2245 = 0
  ixhigh_var_2246 = 0
  ixc_var_2244 = 0
  DO lay_var_2239 = laytrop_min_var_2242 + 1, laytrop_max_var_2243
    icl_var_2248 = 0
    ich_var_2247 = 0
    DO jc_var_2251 = kidia_var_2213, kfdia_var_2214
      IF (lay_var_2239 <= laytrop_var_2226(jc_var_2251)) THEN
        icl_var_2248 = icl_var_2248 + 1
        ixlow_var_2245(icl_var_2248, lay_var_2239) = jc_var_2251
      ELSE
        ich_var_2247 = ich_var_2247 + 1
        ixhigh_var_2246(ich_var_2247, lay_var_2239) = jc_var_2251
      END IF
    END DO
    ixc_var_2244(lay_var_2239) = icl_var_2248
  END DO
  DO lay_var_2239 = 1, laytrop_min_var_2242
    DO jl_var_2252 = kidia_var_2213, kfdia_var_2214
      ind0_var_2234 = ((jp_var_2222(jl_var_2252, lay_var_2239) - 1) * 5 + (jt_var_2223(jl_var_2252, lay_var_2239) - 1)) * nspa_var_216(10) + 1
      ind1_var_2235 = (jp_var_2222(jl_var_2252, lay_var_2239) * 5 + (jt1_var_2224(jl_var_2252, lay_var_2239) - 1)) * nspa_var_216(10) + 1
      inds_var_2236 = indself_var_2230(jl_var_2252, lay_var_2239)
      indf_var_2237 = indfor_var_2231(jl_var_2252, lay_var_2239)
      DO ig_var_2238 = 1, 6
        tauself_var_2241 = selffac_var_2228(jl_var_2252, lay_var_2239) * (selfref_var_118(inds_var_2236, ig_var_2238) + selffrac_var_2229(jl_var_2252, lay_var_2239) * (selfref_var_118(inds_var_2236 + 1, ig_var_2238) - selfref_var_118(inds_var_2236, ig_var_2238)))
        taufor_var_2240 = forfac_var_2233(jl_var_2252, lay_var_2239) * (forref_var_119(indf_var_2237, ig_var_2238) + forfrac_var_2232(jl_var_2252, lay_var_2239) * (forref_var_119(indf_var_2237 + 1, ig_var_2238) - forref_var_119(indf_var_2237, ig_var_2238)))
        taug_var_2216(jl_var_2252, 108 + ig_var_2238, lay_var_2239) = colh2o_var_2225(jl_var_2252, lay_var_2239) * (fac00_var_2218(jl_var_2252, lay_var_2239) * absa_var_116(ind0_var_2234, ig_var_2238) + fac10_var_2220(jl_var_2252, lay_var_2239) * absa_var_116(ind0_var_2234 + 1, ig_var_2238) + fac01_var_2219(jl_var_2252, lay_var_2239) * absa_var_116(ind1_var_2235, ig_var_2238) + fac11_var_2221(jl_var_2252, lay_var_2239) * absa_var_116(ind1_var_2235 + 1, ig_var_2238)) + tauself_var_2241 + taufor_var_2240
        fracs_var_2227(jl_var_2252, 108 + ig_var_2238, lay_var_2239) = fracrefa_var_114(ig_var_2238)
      END DO
    END DO
  END DO
  DO lay_var_2239 = laytrop_max_var_2243 + 1, klev_var_2215
    DO jl_var_2252 = kidia_var_2213, kfdia_var_2214
      ind0_var_2234 = ((jp_var_2222(jl_var_2252, lay_var_2239) - 13) * 5 + (jt_var_2223(jl_var_2252, lay_var_2239) - 1)) * nspb_var_217(10) + 1
      ind1_var_2235 = ((jp_var_2222(jl_var_2252, lay_var_2239) - 12) * 5 + (jt1_var_2224(jl_var_2252, lay_var_2239) - 1)) * nspb_var_217(10) + 1
      indf_var_2237 = indfor_var_2231(jl_var_2252, lay_var_2239)
      DO ig_var_2238 = 1, 6
        taufor_var_2240 = forfac_var_2233(jl_var_2252, lay_var_2239) * (forref_var_119(indf_var_2237, ig_var_2238) + forfrac_var_2232(jl_var_2252, lay_var_2239) * (forref_var_119(indf_var_2237 + 1, ig_var_2238) - forref_var_119(indf_var_2237, ig_var_2238)))
        taug_var_2216(jl_var_2252, 108 + ig_var_2238, lay_var_2239) = colh2o_var_2225(jl_var_2252, lay_var_2239) * (fac00_var_2218(jl_var_2252, lay_var_2239) * absb_var_117(ind0_var_2234, ig_var_2238) + fac10_var_2220(jl_var_2252, lay_var_2239) * absb_var_117(ind0_var_2234 + 1, ig_var_2238) + fac01_var_2219(jl_var_2252, lay_var_2239) * absb_var_117(ind1_var_2235, ig_var_2238) + fac11_var_2221(jl_var_2252, lay_var_2239) * absb_var_117(ind1_var_2235 + 1, ig_var_2238)) + taufor_var_2240
        fracs_var_2227(jl_var_2252, 108 + ig_var_2238, lay_var_2239) = fracrefb_var_115(ig_var_2238)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2243 /= laytrop_min_var_2242) THEN
    DO lay_var_2239 = laytrop_min_var_2242 + 1, laytrop_max_var_2243
      ixc0_var_2249 = ixc_var_2244(lay_var_2239)
      DO ixp_var_2250 = 1, ixc0_var_2249
        jl_var_2252 = ixlow_var_2245(ixp_var_2250, lay_var_2239)
        ind0_var_2234 = ((jp_var_2222(jl_var_2252, lay_var_2239) - 1) * 5 + (jt_var_2223(jl_var_2252, lay_var_2239) - 1)) * nspa_var_216(10) + 1
        ind1_var_2235 = (jp_var_2222(jl_var_2252, lay_var_2239) * 5 + (jt1_var_2224(jl_var_2252, lay_var_2239) - 1)) * nspa_var_216(10) + 1
        inds_var_2236 = indself_var_2230(jl_var_2252, lay_var_2239)
        indf_var_2237 = indfor_var_2231(jl_var_2252, lay_var_2239)
        DO ig_var_2238 = 1, 6
          tauself_var_2241 = selffac_var_2228(jl_var_2252, lay_var_2239) * (selfref_var_118(inds_var_2236, ig_var_2238) + selffrac_var_2229(jl_var_2252, lay_var_2239) * (selfref_var_118(inds_var_2236 + 1, ig_var_2238) - selfref_var_118(inds_var_2236, ig_var_2238)))
          taufor_var_2240 = forfac_var_2233(jl_var_2252, lay_var_2239) * (forref_var_119(indf_var_2237, ig_var_2238) + forfrac_var_2232(jl_var_2252, lay_var_2239) * (forref_var_119(indf_var_2237 + 1, ig_var_2238) - forref_var_119(indf_var_2237, ig_var_2238)))
          taug_var_2216(jl_var_2252, 108 + ig_var_2238, lay_var_2239) = colh2o_var_2225(jl_var_2252, lay_var_2239) * (fac00_var_2218(jl_var_2252, lay_var_2239) * absa_var_116(ind0_var_2234, ig_var_2238) + fac10_var_2220(jl_var_2252, lay_var_2239) * absa_var_116(ind0_var_2234 + 1, ig_var_2238) + fac01_var_2219(jl_var_2252, lay_var_2239) * absa_var_116(ind1_var_2235, ig_var_2238) + fac11_var_2221(jl_var_2252, lay_var_2239) * absa_var_116(ind1_var_2235 + 1, ig_var_2238)) + tauself_var_2241 + taufor_var_2240
          fracs_var_2227(jl_var_2252, 108 + ig_var_2238, lay_var_2239) = fracrefa_var_114(ig_var_2238)
        END DO
      END DO
      ixc0_var_2249 = kfdia_var_2214 - kidia_var_2213 + 1 - ixc0_var_2249
      DO ixp_var_2250 = 1, ixc0_var_2249
        jl_var_2252 = ixhigh_var_2246(ixp_var_2250, lay_var_2239)
        ind0_var_2234 = ((jp_var_2222(jl_var_2252, lay_var_2239) - 13) * 5 + (jt_var_2223(jl_var_2252, lay_var_2239) - 1)) * nspb_var_217(10) + 1
        ind1_var_2235 = ((jp_var_2222(jl_var_2252, lay_var_2239) - 12) * 5 + (jt1_var_2224(jl_var_2252, lay_var_2239) - 1)) * nspb_var_217(10) + 1
        indf_var_2237 = indfor_var_2231(jl_var_2252, lay_var_2239)
        DO ig_var_2238 = 1, 6
          taufor_var_2240 = forfac_var_2233(jl_var_2252, lay_var_2239) * (forref_var_119(indf_var_2237, ig_var_2238) + forfrac_var_2232(jl_var_2252, lay_var_2239) * (forref_var_119(indf_var_2237 + 1, ig_var_2238) - forref_var_119(indf_var_2237, ig_var_2238)))
          taug_var_2216(jl_var_2252, 108 + ig_var_2238, lay_var_2239) = colh2o_var_2225(jl_var_2252, lay_var_2239) * (fac00_var_2218(jl_var_2252, lay_var_2239) * absb_var_117(ind0_var_2234, ig_var_2238) + fac10_var_2220(jl_var_2252, lay_var_2239) * absb_var_117(ind0_var_2234 + 1, ig_var_2238) + fac01_var_2219(jl_var_2252, lay_var_2239) * absb_var_117(ind1_var_2235, ig_var_2238) + fac11_var_2221(jl_var_2252, lay_var_2239) * absb_var_117(ind1_var_2235 + 1, ig_var_2238)) + taufor_var_2240
          fracs_var_2227(jl_var_2252, 108 + ig_var_2238, lay_var_2239) = fracrefb_var_115(ig_var_2238)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol10
SUBROUTINE srtm_taumol23(kidia_var_2253, kfdia_var_2254, klev_var_2255, p_fac00_var_2256, p_fac01_var_2257, p_fac10_var_2258, p_fac11_var_2259, k_jp_var_2260, k_jt_var_2261, k_jt1_var_2262, p_colh2o_var_2263, p_colmol_var_2264, k_laytrop_var_2265, p_selffac_var_2266, p_selffrac_var_2267, k_indself_var_2268, p_forfac_var_2269, p_forfrac_var_2270, k_indfor_var_2271, p_sfluxzen_var_2272, p_taug_var_2273, p_taur_var_2274, prmu0_var_2275)
  USE yoesrta23, ONLY: absa_var_273, forrefc_var_275, givfac, layreffr_var_272, raylc_var_277, selfrefc_var_274, sfluxrefc_var_276
  USE yoesrtwn, ONLY: nspa_var_313
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2253, kfdia_var_2254
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2255
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_2256(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_2257(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_2258(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_2259(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_2260(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_2261(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_2262(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_2263(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_2264(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_2265(kidia_var_2253 : kfdia_var_2254)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_2266(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_2267(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_2268(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_2269(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_2270(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_2271(kidia_var_2253 : kfdia_var_2254, klev_var_2255)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_2272(kidia_var_2253 : kfdia_var_2254, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_2273(kidia_var_2253 : kfdia_var_2254, klev_var_2255, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_2274(kidia_var_2253 : kfdia_var_2254, klev_var_2255, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_2275(kidia_var_2253 : kfdia_var_2254)
  INTEGER(KIND = 4) :: ig_var_2276, ind0_var_2277, ind1_var_2278, inds_var_2279, indf_var_2280, i_lay_var_2281, i_laysolfr_var_2282(kidia_var_2253 : kfdia_var_2254), i_nlayers_var_2283, iplon_var_2284
  INTEGER(KIND = 4) :: laytrop_min_var_2285, laytrop_max_var_2286
  REAL(KIND = 8) :: z_tauray_var_2287
  laytrop_min_var_2285 = MINVAL(k_laytrop_var_2265(kidia_var_2253 : kfdia_var_2254))
  laytrop_max_var_2286 = MAXVAL(k_laytrop_var_2265(kidia_var_2253 : kfdia_var_2254))
  i_nlayers_var_2283 = klev_var_2255
  DO iplon_var_2284 = kidia_var_2253, kfdia_var_2254
    i_laysolfr_var_2282(iplon_var_2284) = k_laytrop_var_2265(iplon_var_2284)
  END DO
  DO i_lay_var_2281 = 1, laytrop_min_var_2285
    DO iplon_var_2284 = kidia_var_2253, kfdia_var_2254
      IF (k_jp_var_2260(iplon_var_2284, i_lay_var_2281) < layreffr_var_272 .AND. k_jp_var_2260(iplon_var_2284, i_lay_var_2281 + 1) >= layreffr_var_272) i_laysolfr_var_2282(iplon_var_2284) = MIN(i_lay_var_2281 + 1, k_laytrop_var_2265(iplon_var_2284))
      ind0_var_2277 = ((k_jp_var_2260(iplon_var_2284, i_lay_var_2281) - 1) * 5 + (k_jt_var_2261(iplon_var_2284, i_lay_var_2281) - 1)) * nspa_var_313(23) + 1
      ind1_var_2278 = (k_jp_var_2260(iplon_var_2284, i_lay_var_2281) * 5 + (k_jt1_var_2262(iplon_var_2284, i_lay_var_2281) - 1)) * nspa_var_313(23) + 1
      inds_var_2279 = k_indself_var_2268(iplon_var_2284, i_lay_var_2281)
      indf_var_2280 = k_indfor_var_2271(iplon_var_2284, i_lay_var_2281)
      DO ig_var_2276 = 1, 10
        z_tauray_var_2287 = p_colmol_var_2264(iplon_var_2284, i_lay_var_2281) * raylc_var_277(ig_var_2276)
        p_taug_var_2273(iplon_var_2284, i_lay_var_2281, ig_var_2276) = p_colh2o_var_2263(iplon_var_2284, i_lay_var_2281) * (givfac * (p_fac00_var_2256(iplon_var_2284, i_lay_var_2281) * absa_var_273(ind0_var_2277, ig_var_2276) + p_fac10_var_2258(iplon_var_2284, i_lay_var_2281) * absa_var_273(ind0_var_2277 + 1, ig_var_2276) + p_fac01_var_2257(iplon_var_2284, i_lay_var_2281) * absa_var_273(ind1_var_2278, ig_var_2276) + p_fac11_var_2259(iplon_var_2284, i_lay_var_2281) * absa_var_273(ind1_var_2278 + 1, ig_var_2276)) + p_selffac_var_2266(iplon_var_2284, i_lay_var_2281) * (selfrefc_var_274(inds_var_2279, ig_var_2276) + p_selffrac_var_2267(iplon_var_2284, i_lay_var_2281) * (selfrefc_var_274(inds_var_2279 + 1, ig_var_2276) - selfrefc_var_274(inds_var_2279, ig_var_2276))) + p_forfac_var_2269(iplon_var_2284, i_lay_var_2281) * (forrefc_var_275(indf_var_2280, ig_var_2276) + p_forfrac_var_2270(iplon_var_2284, i_lay_var_2281) * (forrefc_var_275(indf_var_2280 + 1, ig_var_2276) - forrefc_var_275(indf_var_2280, ig_var_2276))))
        IF (i_lay_var_2281 == i_laysolfr_var_2282(iplon_var_2284)) p_sfluxzen_var_2272(iplon_var_2284, ig_var_2276) = sfluxrefc_var_276(ig_var_2276)
        p_taur_var_2274(iplon_var_2284, i_lay_var_2281, ig_var_2276) = z_tauray_var_2287
      END DO
    END DO
  END DO
  DO i_lay_var_2281 = laytrop_min_var_2285 + 1, laytrop_max_var_2286
    DO iplon_var_2284 = kidia_var_2253, kfdia_var_2254
      IF (i_lay_var_2281 <= k_laytrop_var_2265(iplon_var_2284)) THEN
        IF (k_jp_var_2260(iplon_var_2284, i_lay_var_2281) < layreffr_var_272 .AND. k_jp_var_2260(iplon_var_2284, i_lay_var_2281 + 1) >= layreffr_var_272) i_laysolfr_var_2282(iplon_var_2284) = MIN(i_lay_var_2281 + 1, k_laytrop_var_2265(iplon_var_2284))
        ind0_var_2277 = ((k_jp_var_2260(iplon_var_2284, i_lay_var_2281) - 1) * 5 + (k_jt_var_2261(iplon_var_2284, i_lay_var_2281) - 1)) * nspa_var_313(23) + 1
        ind1_var_2278 = (k_jp_var_2260(iplon_var_2284, i_lay_var_2281) * 5 + (k_jt1_var_2262(iplon_var_2284, i_lay_var_2281) - 1)) * nspa_var_313(23) + 1
        inds_var_2279 = k_indself_var_2268(iplon_var_2284, i_lay_var_2281)
        indf_var_2280 = k_indfor_var_2271(iplon_var_2284, i_lay_var_2281)
        DO ig_var_2276 = 1, 10
          z_tauray_var_2287 = p_colmol_var_2264(iplon_var_2284, i_lay_var_2281) * raylc_var_277(ig_var_2276)
          p_taug_var_2273(iplon_var_2284, i_lay_var_2281, ig_var_2276) = p_colh2o_var_2263(iplon_var_2284, i_lay_var_2281) * (givfac * (p_fac00_var_2256(iplon_var_2284, i_lay_var_2281) * absa_var_273(ind0_var_2277, ig_var_2276) + p_fac10_var_2258(iplon_var_2284, i_lay_var_2281) * absa_var_273(ind0_var_2277 + 1, ig_var_2276) + p_fac01_var_2257(iplon_var_2284, i_lay_var_2281) * absa_var_273(ind1_var_2278, ig_var_2276) + p_fac11_var_2259(iplon_var_2284, i_lay_var_2281) * absa_var_273(ind1_var_2278 + 1, ig_var_2276)) + p_selffac_var_2266(iplon_var_2284, i_lay_var_2281) * (selfrefc_var_274(inds_var_2279, ig_var_2276) + p_selffrac_var_2267(iplon_var_2284, i_lay_var_2281) * (selfrefc_var_274(inds_var_2279 + 1, ig_var_2276) - selfrefc_var_274(inds_var_2279, ig_var_2276))) + p_forfac_var_2269(iplon_var_2284, i_lay_var_2281) * (forrefc_var_275(indf_var_2280, ig_var_2276) + p_forfrac_var_2270(iplon_var_2284, i_lay_var_2281) * (forrefc_var_275(indf_var_2280 + 1, ig_var_2276) - forrefc_var_275(indf_var_2280, ig_var_2276))))
          IF (i_lay_var_2281 == i_laysolfr_var_2282(iplon_var_2284)) p_sfluxzen_var_2272(iplon_var_2284, ig_var_2276) = sfluxrefc_var_276(ig_var_2276)
          p_taur_var_2274(iplon_var_2284, i_lay_var_2281, ig_var_2276) = z_tauray_var_2287
        END DO
      ELSE
        DO ig_var_2276 = 1, 10
          p_taug_var_2273(iplon_var_2284, i_lay_var_2281, ig_var_2276) = 0.0D0
          p_taur_var_2274(iplon_var_2284, i_lay_var_2281, ig_var_2276) = p_colmol_var_2264(iplon_var_2284, i_lay_var_2281) * raylc_var_277(ig_var_2276)
        END DO
      END IF
    END DO
  END DO
  DO ig_var_2276 = 1, 10
    DO i_lay_var_2281 = laytrop_max_var_2286 + 1, i_nlayers_var_2283
      DO iplon_var_2284 = kidia_var_2253, kfdia_var_2254
        p_taug_var_2273(iplon_var_2284, i_lay_var_2281, ig_var_2276) = 0.0D0
        p_taur_var_2274(iplon_var_2284, i_lay_var_2281, ig_var_2276) = p_colmol_var_2264(iplon_var_2284, i_lay_var_2281) * raylc_var_277(ig_var_2276)
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol23
SUBROUTINE rrtm_taumol4(kidia_var_2288, kfdia_var_2289, klev_var_2290, taug_var_2291, p_tauaerl_var_2292, fac00_var_2293, fac01_var_2294, fac10_var_2295, fac11_var_2296, forfac_var_2314, forfrac_var_2315, indfor_var_2313, jp_var_2297, jt_var_2298, jt1_var_2299, oneminus_var_2300, colh2o_var_2301, colco2_var_2302, colo3_var_2303, laytrop_var_2304, selffac_var_2305, selffrac_var_2306, indself_var_2307, fracs_var_2308, rat_h2oco2_var_2309, rat_h2oco2_1_var_2310, rat_o3co2_var_2311, rat_o3co2_1_var_2312)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng4
  USE yoerrta4, ONLY: absa_var_169, absb_var_170, forref_var_172, fracrefa_var_167, fracrefb_var_168, selfref_var_171
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2288
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2289
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2290
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2291(kidia_var_2288 : kfdia_var_2289, 140, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2292(kidia_var_2288 : kfdia_var_2289, klev_var_2290, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2293(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2294(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2295(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2296(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2297(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2298(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2299(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2300
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2301(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2302(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_2303(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2304(kidia_var_2288 : kfdia_var_2289)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2305(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2306(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2307(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2308(kidia_var_2288 : kfdia_var_2289, 140, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_2309(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_2310(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_var_2311(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_1_var_2312(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2313(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2314(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2315(kidia_var_2288 : kfdia_var_2289, klev_var_2290)
  REAL(KIND = 8) :: speccomb_var_2316, speccomb1_var_2317, speccomb_planck_var_2318
  INTEGER(KIND = 4) :: ind0_var_2319, ind1_var_2320, inds_var_2321, indf_var_2322
  INTEGER(KIND = 4) :: ig_var_2323, js_var_2324, lay_var_2325, js1_var_2326, jpl_var_2327
  REAL(KIND = 8) :: refrat_planck_a_var_2328, refrat_planck_b_var_2329
  REAL(KIND = 8) :: fac000_var_2330, fac100_var_2331, fac200_var_2332, fac010_var_2333, fac110_var_2334, fac210_var_2335, fac001_var_2336, fac101_var_2337, fac201_var_2338, fac011_var_2339, fac111_var_2340, fac211_var_2341
  REAL(KIND = 8) :: p_var_2342, p4_var_2343, fk0_var_2344, fk1_var_2345, fk2_var_2346
  REAL(KIND = 8) :: taufor_var_2347, tauself_var_2348, tau_major_var_2349(14), tau_major1_var_2350(14)
  REAL(KIND = 8) :: fs_var_2351, specmult_var_2352, specparm_var_2353, fs1_var_2354, specmult1_var_2355, specparm1_var_2356, fpl_var_2357, specmult_planck_var_2358, specparm_planck_var_2359
  INTEGER(KIND = 4) :: laytrop_min_var_2360, laytrop_max_var_2361
  INTEGER(KIND = 4) :: ixc_var_2362(klev_var_2290), ixlow_var_2363(kfdia_var_2289, klev_var_2290), ixhigh_var_2364(kfdia_var_2289, klev_var_2290)
  INTEGER(KIND = 4) :: ich_var_2365, icl_var_2366, ixc0_var_2367, ixp_var_2368, jc_var_2369, jl_var_2370
  laytrop_min_var_2360 = MINVAL(laytrop_var_2304)
  laytrop_max_var_2361 = MAXVAL(laytrop_var_2304)
  ixlow_var_2363 = 0
  ixhigh_var_2364 = 0
  ixc_var_2362 = 0
  DO lay_var_2325 = laytrop_min_var_2360 + 1, laytrop_max_var_2361
    icl_var_2366 = 0
    ich_var_2365 = 0
    DO jc_var_2369 = kidia_var_2288, kfdia_var_2289
      IF (lay_var_2325 <= laytrop_var_2304(jc_var_2369)) THEN
        icl_var_2366 = icl_var_2366 + 1
        ixlow_var_2363(icl_var_2366, lay_var_2325) = jc_var_2369
      ELSE
        ich_var_2365 = ich_var_2365 + 1
        ixhigh_var_2364(ich_var_2365, lay_var_2325) = jc_var_2369
      END IF
    END DO
    ixc_var_2362(lay_var_2325) = icl_var_2366
  END DO
  refrat_planck_a_var_2328 = chi_mls(1, 11) / chi_mls(2, 11)
  refrat_planck_b_var_2329 = chi_mls(3, 13) / chi_mls(2, 13)
  DO lay_var_2325 = 1, laytrop_min_var_2360
    DO jl_var_2370 = kidia_var_2288, kfdia_var_2289
      speccomb_var_2316 = colh2o_var_2301(jl_var_2370, lay_var_2325) + rat_h2oco2_var_2309(jl_var_2370, lay_var_2325) * colco2_var_2302(jl_var_2370, lay_var_2325)
      specparm_var_2353 = MIN(colh2o_var_2301(jl_var_2370, lay_var_2325) / speccomb_var_2316, oneminus_var_2300)
      specmult_var_2352 = 8.0D0 * (specparm_var_2353)
      js_var_2324 = 1 + INT(specmult_var_2352)
      fs_var_2351 = ((specmult_var_2352) - AINT((specmult_var_2352)))
      speccomb1_var_2317 = colh2o_var_2301(jl_var_2370, lay_var_2325) + rat_h2oco2_1_var_2310(jl_var_2370, lay_var_2325) * colco2_var_2302(jl_var_2370, lay_var_2325)
      specparm1_var_2356 = MIN(colh2o_var_2301(jl_var_2370, lay_var_2325) / speccomb1_var_2317, oneminus_var_2300)
      specmult1_var_2355 = 8.0D0 * (specparm1_var_2356)
      js1_var_2326 = 1 + INT(specmult1_var_2355)
      fs1_var_2354 = ((specmult1_var_2355) - AINT((specmult1_var_2355)))
      speccomb_planck_var_2318 = colh2o_var_2301(jl_var_2370, lay_var_2325) + refrat_planck_a_var_2328 * colco2_var_2302(jl_var_2370, lay_var_2325)
      specparm_planck_var_2359 = MIN(colh2o_var_2301(jl_var_2370, lay_var_2325) / speccomb_planck_var_2318, oneminus_var_2300)
      specmult_planck_var_2358 = 8.0D0 * specparm_planck_var_2359
      jpl_var_2327 = 1 + INT(specmult_planck_var_2358)
      fpl_var_2357 = ((specmult_planck_var_2358) - AINT((specmult_planck_var_2358)))
      ind0_var_2319 = ((jp_var_2297(jl_var_2370, lay_var_2325) - 1) * 5 + (jt_var_2298(jl_var_2370, lay_var_2325) - 1)) * nspa_var_216(4) + js_var_2324
      ind1_var_2320 = (jp_var_2297(jl_var_2370, lay_var_2325) * 5 + (jt1_var_2299(jl_var_2370, lay_var_2325) - 1)) * nspa_var_216(4) + js1_var_2326
      inds_var_2321 = indself_var_2307(jl_var_2370, lay_var_2325)
      indf_var_2322 = indfor_var_2313(jl_var_2370, lay_var_2325)
      IF (specparm_var_2353 .LT. 0.125D0) THEN
        p_var_2342 = fs_var_2351 - 1.0D0
        p4_var_2343 = p_var_2342 ** 4
        fk0_var_2344 = p4_var_2343
        fk1_var_2345 = 1.0D0 - p_var_2342 - 2.0D0 * p4_var_2343
        fk2_var_2346 = p_var_2342 + p4_var_2343
        fac000_var_2330 = fk0_var_2344 * fac00_var_2293(jl_var_2370, lay_var_2325)
        fac100_var_2331 = fk1_var_2345 * fac00_var_2293(jl_var_2370, lay_var_2325)
        fac200_var_2332 = fk2_var_2346 * fac00_var_2293(jl_var_2370, lay_var_2325)
        fac010_var_2333 = fk0_var_2344 * fac10_var_2295(jl_var_2370, lay_var_2325)
        fac110_var_2334 = fk1_var_2345 * fac10_var_2295(jl_var_2370, lay_var_2325)
        fac210_var_2335 = fk2_var_2346 * fac10_var_2295(jl_var_2370, lay_var_2325)
      ELSE IF (specparm_var_2353 .GT. 0.875D0) THEN
        p_var_2342 = - fs_var_2351
        p4_var_2343 = p_var_2342 ** 4
        fk0_var_2344 = p4_var_2343
        fk1_var_2345 = 1.0D0 - p_var_2342 - 2.0D0 * p4_var_2343
        fk2_var_2346 = p_var_2342 + p4_var_2343
        fac000_var_2330 = fk0_var_2344 * fac00_var_2293(jl_var_2370, lay_var_2325)
        fac100_var_2331 = fk1_var_2345 * fac00_var_2293(jl_var_2370, lay_var_2325)
        fac200_var_2332 = fk2_var_2346 * fac00_var_2293(jl_var_2370, lay_var_2325)
        fac010_var_2333 = fk0_var_2344 * fac10_var_2295(jl_var_2370, lay_var_2325)
        fac110_var_2334 = fk1_var_2345 * fac10_var_2295(jl_var_2370, lay_var_2325)
        fac210_var_2335 = fk2_var_2346 * fac10_var_2295(jl_var_2370, lay_var_2325)
      ELSE
        fac000_var_2330 = (1.0D0 - fs_var_2351) * fac00_var_2293(jl_var_2370, lay_var_2325)
        fac010_var_2333 = (1.0D0 - fs_var_2351) * fac10_var_2295(jl_var_2370, lay_var_2325)
        fac100_var_2331 = fs_var_2351 * fac00_var_2293(jl_var_2370, lay_var_2325)
        fac110_var_2334 = fs_var_2351 * fac10_var_2295(jl_var_2370, lay_var_2325)
        fac200_var_2332 = 0.0D0
        fac210_var_2335 = 0.0D0
      END IF
      IF (specparm1_var_2356 .LT. 0.125D0) THEN
        p_var_2342 = fs1_var_2354 - 1.0D0
        p4_var_2343 = p_var_2342 ** 4
        fk0_var_2344 = p4_var_2343
        fk1_var_2345 = 1.0D0 - p_var_2342 - 2.0D0 * p4_var_2343
        fk2_var_2346 = p_var_2342 + p4_var_2343
        fac001_var_2336 = fk0_var_2344 * fac01_var_2294(jl_var_2370, lay_var_2325)
        fac101_var_2337 = fk1_var_2345 * fac01_var_2294(jl_var_2370, lay_var_2325)
        fac201_var_2338 = fk2_var_2346 * fac01_var_2294(jl_var_2370, lay_var_2325)
        fac011_var_2339 = fk0_var_2344 * fac11_var_2296(jl_var_2370, lay_var_2325)
        fac111_var_2340 = fk1_var_2345 * fac11_var_2296(jl_var_2370, lay_var_2325)
        fac211_var_2341 = fk2_var_2346 * fac11_var_2296(jl_var_2370, lay_var_2325)
      ELSE IF (specparm1_var_2356 .GT. 0.875D0) THEN
        p_var_2342 = - fs1_var_2354
        p4_var_2343 = p_var_2342 ** 4
        fk0_var_2344 = p4_var_2343
        fk1_var_2345 = 1.0D0 - p_var_2342 - 2.0D0 * p4_var_2343
        fk2_var_2346 = p_var_2342 + p4_var_2343
        fac001_var_2336 = fk0_var_2344 * fac01_var_2294(jl_var_2370, lay_var_2325)
        fac101_var_2337 = fk1_var_2345 * fac01_var_2294(jl_var_2370, lay_var_2325)
        fac201_var_2338 = fk2_var_2346 * fac01_var_2294(jl_var_2370, lay_var_2325)
        fac011_var_2339 = fk0_var_2344 * fac11_var_2296(jl_var_2370, lay_var_2325)
        fac111_var_2340 = fk1_var_2345 * fac11_var_2296(jl_var_2370, lay_var_2325)
        fac211_var_2341 = fk2_var_2346 * fac11_var_2296(jl_var_2370, lay_var_2325)
      ELSE
        fac001_var_2336 = (1.0D0 - fs1_var_2354) * fac01_var_2294(jl_var_2370, lay_var_2325)
        fac011_var_2339 = (1.0D0 - fs1_var_2354) * fac11_var_2296(jl_var_2370, lay_var_2325)
        fac101_var_2337 = fs1_var_2354 * fac01_var_2294(jl_var_2370, lay_var_2325)
        fac111_var_2340 = fs1_var_2354 * fac11_var_2296(jl_var_2370, lay_var_2325)
        fac201_var_2338 = 0.0D0
        fac211_var_2341 = 0.0D0
      END IF
      IF (specparm_var_2353 .LT. 0.125D0) THEN
        tau_major_var_2349(1 : ng4) = speccomb_var_2316 * (fac000_var_2330 * absa_var_169(ind0_var_2319, 1 : 14) + fac100_var_2331 * absa_var_169(ind0_var_2319 + 1, 1 : 14) + fac200_var_2332 * absa_var_169(ind0_var_2319 + 2, 1 : 14) + fac010_var_2333 * absa_var_169(ind0_var_2319 + 9, 1 : 14) + fac110_var_2334 * absa_var_169(ind0_var_2319 + 10, 1 : 14) + fac210_var_2335 * absa_var_169(ind0_var_2319 + 11, 1 : 14))
      ELSE IF (specparm_var_2353 .GT. 0.875D0) THEN
        tau_major_var_2349(1 : ng4) = speccomb_var_2316 * (fac200_var_2332 * absa_var_169(ind0_var_2319 - 1, 1 : 14) + fac100_var_2331 * absa_var_169(ind0_var_2319, 1 : 14) + fac000_var_2330 * absa_var_169(ind0_var_2319 + 1, 1 : 14) + fac210_var_2335 * absa_var_169(ind0_var_2319 + 8, 1 : 14) + fac110_var_2334 * absa_var_169(ind0_var_2319 + 9, 1 : 14) + fac010_var_2333 * absa_var_169(ind0_var_2319 + 10, 1 : 14))
      ELSE
        tau_major_var_2349(1 : ng4) = speccomb_var_2316 * (fac000_var_2330 * absa_var_169(ind0_var_2319, 1 : 14) + fac100_var_2331 * absa_var_169(ind0_var_2319 + 1, 1 : 14) + fac010_var_2333 * absa_var_169(ind0_var_2319 + 9, 1 : 14) + fac110_var_2334 * absa_var_169(ind0_var_2319 + 10, 1 : 14))
      END IF
      IF (specparm1_var_2356 .LT. 0.125D0) THEN
        tau_major1_var_2350(1 : ng4) = speccomb1_var_2317 * (fac001_var_2336 * absa_var_169(ind1_var_2320, 1 : 14) + fac101_var_2337 * absa_var_169(ind1_var_2320 + 1, 1 : 14) + fac201_var_2338 * absa_var_169(ind1_var_2320 + 2, 1 : 14) + fac011_var_2339 * absa_var_169(ind1_var_2320 + 9, 1 : 14) + fac111_var_2340 * absa_var_169(ind1_var_2320 + 10, 1 : 14) + fac211_var_2341 * absa_var_169(ind1_var_2320 + 11, 1 : 14))
      ELSE IF (specparm1_var_2356 .GT. 0.875D0) THEN
        tau_major1_var_2350(1 : ng4) = speccomb1_var_2317 * (fac201_var_2338 * absa_var_169(ind1_var_2320 - 1, 1 : 14) + fac101_var_2337 * absa_var_169(ind1_var_2320, 1 : 14) + fac001_var_2336 * absa_var_169(ind1_var_2320 + 1, 1 : 14) + fac211_var_2341 * absa_var_169(ind1_var_2320 + 8, 1 : 14) + fac111_var_2340 * absa_var_169(ind1_var_2320 + 9, 1 : 14) + fac011_var_2339 * absa_var_169(ind1_var_2320 + 10, 1 : 14))
      ELSE
        tau_major1_var_2350(1 : ng4) = speccomb1_var_2317 * (fac001_var_2336 * absa_var_169(ind1_var_2320, 1 : 14) + fac101_var_2337 * absa_var_169(ind1_var_2320 + 1, 1 : 14) + fac011_var_2339 * absa_var_169(ind1_var_2320 + 9, 1 : 14) + fac111_var_2340 * absa_var_169(ind1_var_2320 + 10, 1 : 14))
      END IF
      DO ig_var_2323 = 1, 14
        tauself_var_2348 = selffac_var_2305(jl_var_2370, lay_var_2325) * (selfref_var_171(inds_var_2321, ig_var_2323) + selffrac_var_2306(jl_var_2370, lay_var_2325) * (selfref_var_171(inds_var_2321 + 1, ig_var_2323) - selfref_var_171(inds_var_2321, ig_var_2323)))
        taufor_var_2347 = forfac_var_2314(jl_var_2370, lay_var_2325) * (forref_var_172(indf_var_2322, ig_var_2323) + forfrac_var_2315(jl_var_2370, lay_var_2325) * (forref_var_172(indf_var_2322 + 1, ig_var_2323) - forref_var_172(indf_var_2322, ig_var_2323)))
        taug_var_2291(jl_var_2370, 38 + ig_var_2323, lay_var_2325) = tau_major_var_2349(ig_var_2323) + tau_major1_var_2350(ig_var_2323) + tauself_var_2348 + taufor_var_2347
        fracs_var_2308(jl_var_2370, 38 + ig_var_2323, lay_var_2325) = fracrefa_var_167(ig_var_2323, jpl_var_2327) + fpl_var_2357 * (fracrefa_var_167(ig_var_2323, jpl_var_2327 + 1) - fracrefa_var_167(ig_var_2323, jpl_var_2327))
      END DO
    END DO
  END DO
  DO lay_var_2325 = laytrop_max_var_2361 + 1, klev_var_2290
    DO jl_var_2370 = kidia_var_2288, kfdia_var_2289
      speccomb_var_2316 = colo3_var_2303(jl_var_2370, lay_var_2325) + rat_o3co2_var_2311(jl_var_2370, lay_var_2325) * colco2_var_2302(jl_var_2370, lay_var_2325)
      specparm_var_2353 = MIN(colo3_var_2303(jl_var_2370, lay_var_2325) / speccomb_var_2316, oneminus_var_2300)
      specmult_var_2352 = 4.0D0 * (specparm_var_2353)
      js_var_2324 = 1 + INT(specmult_var_2352)
      fs_var_2351 = ((specmult_var_2352) - AINT((specmult_var_2352)))
      speccomb1_var_2317 = colo3_var_2303(jl_var_2370, lay_var_2325) + rat_o3co2_1_var_2312(jl_var_2370, lay_var_2325) * colco2_var_2302(jl_var_2370, lay_var_2325)
      specparm1_var_2356 = MIN(colo3_var_2303(jl_var_2370, lay_var_2325) / speccomb1_var_2317, oneminus_var_2300)
      specmult1_var_2355 = 4.0D0 * (specparm1_var_2356)
      js1_var_2326 = 1 + INT(specmult1_var_2355)
      fs1_var_2354 = ((specmult1_var_2355) - AINT((specmult1_var_2355)))
      fac000_var_2330 = (1.0D0 - fs_var_2351) * fac00_var_2293(jl_var_2370, lay_var_2325)
      fac010_var_2333 = (1.0D0 - fs_var_2351) * fac10_var_2295(jl_var_2370, lay_var_2325)
      fac100_var_2331 = fs_var_2351 * fac00_var_2293(jl_var_2370, lay_var_2325)
      fac110_var_2334 = fs_var_2351 * fac10_var_2295(jl_var_2370, lay_var_2325)
      fac001_var_2336 = (1.0D0 - fs1_var_2354) * fac01_var_2294(jl_var_2370, lay_var_2325)
      fac011_var_2339 = (1.0D0 - fs1_var_2354) * fac11_var_2296(jl_var_2370, lay_var_2325)
      fac101_var_2337 = fs1_var_2354 * fac01_var_2294(jl_var_2370, lay_var_2325)
      fac111_var_2340 = fs1_var_2354 * fac11_var_2296(jl_var_2370, lay_var_2325)
      speccomb_planck_var_2318 = colo3_var_2303(jl_var_2370, lay_var_2325) + refrat_planck_b_var_2329 * colco2_var_2302(jl_var_2370, lay_var_2325)
      specparm_planck_var_2359 = MIN(colo3_var_2303(jl_var_2370, lay_var_2325) / speccomb_planck_var_2318, oneminus_var_2300)
      specmult_planck_var_2358 = 4.0D0 * specparm_planck_var_2359
      jpl_var_2327 = 1 + INT(specmult_planck_var_2358)
      fpl_var_2357 = ((specmult_planck_var_2358) - AINT((specmult_planck_var_2358)))
      ind0_var_2319 = ((jp_var_2297(jl_var_2370, lay_var_2325) - 13) * 5 + (jt_var_2298(jl_var_2370, lay_var_2325) - 1)) * nspb_var_217(4) + js_var_2324
      ind1_var_2320 = ((jp_var_2297(jl_var_2370, lay_var_2325) - 12) * 5 + (jt1_var_2299(jl_var_2370, lay_var_2325) - 1)) * nspb_var_217(4) + js1_var_2326
      DO ig_var_2323 = 1, 14
        taug_var_2291(jl_var_2370, 38 + ig_var_2323, lay_var_2325) = speccomb_var_2316 * (fac000_var_2330 * absb_var_170(ind0_var_2319, ig_var_2323) + fac100_var_2331 * absb_var_170(ind0_var_2319 + 1, ig_var_2323) + fac010_var_2333 * absb_var_170(ind0_var_2319 + 5, ig_var_2323) + fac110_var_2334 * absb_var_170(ind0_var_2319 + 6, ig_var_2323)) + speccomb1_var_2317 * (fac001_var_2336 * absb_var_170(ind1_var_2320, ig_var_2323) + fac101_var_2337 * absb_var_170(ind1_var_2320 + 1, ig_var_2323) + fac011_var_2339 * absb_var_170(ind1_var_2320 + 5, ig_var_2323) + fac111_var_2340 * absb_var_170(ind1_var_2320 + 6, ig_var_2323))
        fracs_var_2308(jl_var_2370, 38 + ig_var_2323, lay_var_2325) = fracrefb_var_168(ig_var_2323, jpl_var_2327) + fpl_var_2357 * (fracrefb_var_168(ig_var_2323, jpl_var_2327 + 1) - fracrefb_var_168(ig_var_2323, jpl_var_2327))
      END DO
    END DO
  END DO
  DO lay_var_2325 = laytrop_max_var_2361 + 1, klev_var_2290
    DO jl_var_2370 = kidia_var_2288, kfdia_var_2289
      taug_var_2291(jl_var_2370, 46, lay_var_2325) = taug_var_2291(jl_var_2370, 46, lay_var_2325) * 0.92D0
      taug_var_2291(jl_var_2370, 47, lay_var_2325) = taug_var_2291(jl_var_2370, 47, lay_var_2325) * 0.88D0
      taug_var_2291(jl_var_2370, 48, lay_var_2325) = taug_var_2291(jl_var_2370, 48, lay_var_2325) * 1.07D0
      taug_var_2291(jl_var_2370, 49, lay_var_2325) = taug_var_2291(jl_var_2370, 49, lay_var_2325) * 1.1D0
      taug_var_2291(jl_var_2370, 50, lay_var_2325) = taug_var_2291(jl_var_2370, 50, lay_var_2325) * 0.99D0
      taug_var_2291(jl_var_2370, 51, lay_var_2325) = taug_var_2291(jl_var_2370, 51, lay_var_2325) * 0.88D0
      taug_var_2291(jl_var_2370, 52, lay_var_2325) = taug_var_2291(jl_var_2370, 52, lay_var_2325) * 0.943D0
    END DO
  END DO
  IF (laytrop_max_var_2361 /= laytrop_min_var_2360) THEN
    DO lay_var_2325 = laytrop_min_var_2360 + 1, laytrop_max_var_2361
      ixc0_var_2367 = ixc_var_2362(lay_var_2325)
      DO ixp_var_2368 = 1, ixc0_var_2367
        jl_var_2370 = ixlow_var_2363(ixp_var_2368, lay_var_2325)
        speccomb_var_2316 = colh2o_var_2301(jl_var_2370, lay_var_2325) + rat_h2oco2_var_2309(jl_var_2370, lay_var_2325) * colco2_var_2302(jl_var_2370, lay_var_2325)
        specparm_var_2353 = MIN(colh2o_var_2301(jl_var_2370, lay_var_2325) / speccomb_var_2316, oneminus_var_2300)
        specmult_var_2352 = 8.0D0 * (specparm_var_2353)
        js_var_2324 = 1 + INT(specmult_var_2352)
        fs_var_2351 = ((specmult_var_2352) - AINT((specmult_var_2352)))
        speccomb1_var_2317 = colh2o_var_2301(jl_var_2370, lay_var_2325) + rat_h2oco2_1_var_2310(jl_var_2370, lay_var_2325) * colco2_var_2302(jl_var_2370, lay_var_2325)
        specparm1_var_2356 = MIN(colh2o_var_2301(jl_var_2370, lay_var_2325) / speccomb1_var_2317, oneminus_var_2300)
        specmult1_var_2355 = 8.0D0 * (specparm1_var_2356)
        js1_var_2326 = 1 + INT(specmult1_var_2355)
        fs1_var_2354 = ((specmult1_var_2355) - AINT((specmult1_var_2355)))
        speccomb_planck_var_2318 = colh2o_var_2301(jl_var_2370, lay_var_2325) + refrat_planck_a_var_2328 * colco2_var_2302(jl_var_2370, lay_var_2325)
        specparm_planck_var_2359 = MIN(colh2o_var_2301(jl_var_2370, lay_var_2325) / speccomb_planck_var_2318, oneminus_var_2300)
        specmult_planck_var_2358 = 8.0D0 * specparm_planck_var_2359
        jpl_var_2327 = 1 + INT(specmult_planck_var_2358)
        fpl_var_2357 = ((specmult_planck_var_2358) - AINT((specmult_planck_var_2358)))
        ind0_var_2319 = ((jp_var_2297(jl_var_2370, lay_var_2325) - 1) * 5 + (jt_var_2298(jl_var_2370, lay_var_2325) - 1)) * nspa_var_216(4) + js_var_2324
        ind1_var_2320 = (jp_var_2297(jl_var_2370, lay_var_2325) * 5 + (jt1_var_2299(jl_var_2370, lay_var_2325) - 1)) * nspa_var_216(4) + js1_var_2326
        inds_var_2321 = indself_var_2307(jl_var_2370, lay_var_2325)
        indf_var_2322 = indfor_var_2313(jl_var_2370, lay_var_2325)
        IF (specparm_var_2353 .LT. 0.125D0) THEN
          p_var_2342 = fs_var_2351 - 1.0D0
          p4_var_2343 = p_var_2342 ** 4
          fk0_var_2344 = p4_var_2343
          fk1_var_2345 = 1.0D0 - p_var_2342 - 2.0D0 * p4_var_2343
          fk2_var_2346 = p_var_2342 + p4_var_2343
          fac000_var_2330 = fk0_var_2344 * fac00_var_2293(jl_var_2370, lay_var_2325)
          fac100_var_2331 = fk1_var_2345 * fac00_var_2293(jl_var_2370, lay_var_2325)
          fac200_var_2332 = fk2_var_2346 * fac00_var_2293(jl_var_2370, lay_var_2325)
          fac010_var_2333 = fk0_var_2344 * fac10_var_2295(jl_var_2370, lay_var_2325)
          fac110_var_2334 = fk1_var_2345 * fac10_var_2295(jl_var_2370, lay_var_2325)
          fac210_var_2335 = fk2_var_2346 * fac10_var_2295(jl_var_2370, lay_var_2325)
        ELSE IF (specparm_var_2353 .GT. 0.875D0) THEN
          p_var_2342 = - fs_var_2351
          p4_var_2343 = p_var_2342 ** 4
          fk0_var_2344 = p4_var_2343
          fk1_var_2345 = 1.0D0 - p_var_2342 - 2.0D0 * p4_var_2343
          fk2_var_2346 = p_var_2342 + p4_var_2343
          fac000_var_2330 = fk0_var_2344 * fac00_var_2293(jl_var_2370, lay_var_2325)
          fac100_var_2331 = fk1_var_2345 * fac00_var_2293(jl_var_2370, lay_var_2325)
          fac200_var_2332 = fk2_var_2346 * fac00_var_2293(jl_var_2370, lay_var_2325)
          fac010_var_2333 = fk0_var_2344 * fac10_var_2295(jl_var_2370, lay_var_2325)
          fac110_var_2334 = fk1_var_2345 * fac10_var_2295(jl_var_2370, lay_var_2325)
          fac210_var_2335 = fk2_var_2346 * fac10_var_2295(jl_var_2370, lay_var_2325)
        ELSE
          fac000_var_2330 = (1.0D0 - fs_var_2351) * fac00_var_2293(jl_var_2370, lay_var_2325)
          fac010_var_2333 = (1.0D0 - fs_var_2351) * fac10_var_2295(jl_var_2370, lay_var_2325)
          fac100_var_2331 = fs_var_2351 * fac00_var_2293(jl_var_2370, lay_var_2325)
          fac110_var_2334 = fs_var_2351 * fac10_var_2295(jl_var_2370, lay_var_2325)
          fac200_var_2332 = 0.0D0
          fac210_var_2335 = 0.0D0
        END IF
        IF (specparm1_var_2356 .LT. 0.125D0) THEN
          p_var_2342 = fs1_var_2354 - 1.0D0
          p4_var_2343 = p_var_2342 ** 4
          fk0_var_2344 = p4_var_2343
          fk1_var_2345 = 1.0D0 - p_var_2342 - 2.0D0 * p4_var_2343
          fk2_var_2346 = p_var_2342 + p4_var_2343
          fac001_var_2336 = fk0_var_2344 * fac01_var_2294(jl_var_2370, lay_var_2325)
          fac101_var_2337 = fk1_var_2345 * fac01_var_2294(jl_var_2370, lay_var_2325)
          fac201_var_2338 = fk2_var_2346 * fac01_var_2294(jl_var_2370, lay_var_2325)
          fac011_var_2339 = fk0_var_2344 * fac11_var_2296(jl_var_2370, lay_var_2325)
          fac111_var_2340 = fk1_var_2345 * fac11_var_2296(jl_var_2370, lay_var_2325)
          fac211_var_2341 = fk2_var_2346 * fac11_var_2296(jl_var_2370, lay_var_2325)
        ELSE IF (specparm1_var_2356 .GT. 0.875D0) THEN
          p_var_2342 = - fs1_var_2354
          p4_var_2343 = p_var_2342 ** 4
          fk0_var_2344 = p4_var_2343
          fk1_var_2345 = 1.0D0 - p_var_2342 - 2.0D0 * p4_var_2343
          fk2_var_2346 = p_var_2342 + p4_var_2343
          fac001_var_2336 = fk0_var_2344 * fac01_var_2294(jl_var_2370, lay_var_2325)
          fac101_var_2337 = fk1_var_2345 * fac01_var_2294(jl_var_2370, lay_var_2325)
          fac201_var_2338 = fk2_var_2346 * fac01_var_2294(jl_var_2370, lay_var_2325)
          fac011_var_2339 = fk0_var_2344 * fac11_var_2296(jl_var_2370, lay_var_2325)
          fac111_var_2340 = fk1_var_2345 * fac11_var_2296(jl_var_2370, lay_var_2325)
          fac211_var_2341 = fk2_var_2346 * fac11_var_2296(jl_var_2370, lay_var_2325)
        ELSE
          fac001_var_2336 = (1.0D0 - fs1_var_2354) * fac01_var_2294(jl_var_2370, lay_var_2325)
          fac011_var_2339 = (1.0D0 - fs1_var_2354) * fac11_var_2296(jl_var_2370, lay_var_2325)
          fac101_var_2337 = fs1_var_2354 * fac01_var_2294(jl_var_2370, lay_var_2325)
          fac111_var_2340 = fs1_var_2354 * fac11_var_2296(jl_var_2370, lay_var_2325)
          fac201_var_2338 = 0.0D0
          fac211_var_2341 = 0.0D0
        END IF
        IF (specparm_var_2353 .LT. 0.125D0) THEN
          tau_major_var_2349(1 : ng4) = speccomb_var_2316 * (fac000_var_2330 * absa_var_169(ind0_var_2319, 1 : 14) + fac100_var_2331 * absa_var_169(ind0_var_2319 + 1, 1 : 14) + fac200_var_2332 * absa_var_169(ind0_var_2319 + 2, 1 : 14) + fac010_var_2333 * absa_var_169(ind0_var_2319 + 9, 1 : 14) + fac110_var_2334 * absa_var_169(ind0_var_2319 + 10, 1 : 14) + fac210_var_2335 * absa_var_169(ind0_var_2319 + 11, 1 : 14))
        ELSE IF (specparm_var_2353 .GT. 0.875D0) THEN
          tau_major_var_2349(1 : ng4) = speccomb_var_2316 * (fac200_var_2332 * absa_var_169(ind0_var_2319 - 1, 1 : 14) + fac100_var_2331 * absa_var_169(ind0_var_2319, 1 : 14) + fac000_var_2330 * absa_var_169(ind0_var_2319 + 1, 1 : 14) + fac210_var_2335 * absa_var_169(ind0_var_2319 + 8, 1 : 14) + fac110_var_2334 * absa_var_169(ind0_var_2319 + 9, 1 : 14) + fac010_var_2333 * absa_var_169(ind0_var_2319 + 10, 1 : 14))
        ELSE
          tau_major_var_2349(1 : ng4) = speccomb_var_2316 * (fac000_var_2330 * absa_var_169(ind0_var_2319, 1 : 14) + fac100_var_2331 * absa_var_169(ind0_var_2319 + 1, 1 : 14) + fac010_var_2333 * absa_var_169(ind0_var_2319 + 9, 1 : 14) + fac110_var_2334 * absa_var_169(ind0_var_2319 + 10, 1 : 14))
        END IF
        IF (specparm1_var_2356 .LT. 0.125D0) THEN
          tau_major1_var_2350(1 : ng4) = speccomb1_var_2317 * (fac001_var_2336 * absa_var_169(ind1_var_2320, 1 : 14) + fac101_var_2337 * absa_var_169(ind1_var_2320 + 1, 1 : 14) + fac201_var_2338 * absa_var_169(ind1_var_2320 + 2, 1 : 14) + fac011_var_2339 * absa_var_169(ind1_var_2320 + 9, 1 : 14) + fac111_var_2340 * absa_var_169(ind1_var_2320 + 10, 1 : 14) + fac211_var_2341 * absa_var_169(ind1_var_2320 + 11, 1 : 14))
        ELSE IF (specparm1_var_2356 .GT. 0.875D0) THEN
          tau_major1_var_2350(1 : ng4) = speccomb1_var_2317 * (fac201_var_2338 * absa_var_169(ind1_var_2320 - 1, 1 : 14) + fac101_var_2337 * absa_var_169(ind1_var_2320, 1 : 14) + fac001_var_2336 * absa_var_169(ind1_var_2320 + 1, 1 : 14) + fac211_var_2341 * absa_var_169(ind1_var_2320 + 8, 1 : 14) + fac111_var_2340 * absa_var_169(ind1_var_2320 + 9, 1 : 14) + fac011_var_2339 * absa_var_169(ind1_var_2320 + 10, 1 : 14))
        ELSE
          tau_major1_var_2350(1 : ng4) = speccomb1_var_2317 * (fac001_var_2336 * absa_var_169(ind1_var_2320, 1 : 14) + fac101_var_2337 * absa_var_169(ind1_var_2320 + 1, 1 : 14) + fac011_var_2339 * absa_var_169(ind1_var_2320 + 9, 1 : 14) + fac111_var_2340 * absa_var_169(ind1_var_2320 + 10, 1 : 14))
        END IF
        DO ig_var_2323 = 1, 14
          tauself_var_2348 = selffac_var_2305(jl_var_2370, lay_var_2325) * (selfref_var_171(inds_var_2321, ig_var_2323) + selffrac_var_2306(jl_var_2370, lay_var_2325) * (selfref_var_171(inds_var_2321 + 1, ig_var_2323) - selfref_var_171(inds_var_2321, ig_var_2323)))
          taufor_var_2347 = forfac_var_2314(jl_var_2370, lay_var_2325) * (forref_var_172(indf_var_2322, ig_var_2323) + forfrac_var_2315(jl_var_2370, lay_var_2325) * (forref_var_172(indf_var_2322 + 1, ig_var_2323) - forref_var_172(indf_var_2322, ig_var_2323)))
          taug_var_2291(jl_var_2370, 38 + ig_var_2323, lay_var_2325) = tau_major_var_2349(ig_var_2323) + tau_major1_var_2350(ig_var_2323) + tauself_var_2348 + taufor_var_2347
          fracs_var_2308(jl_var_2370, 38 + ig_var_2323, lay_var_2325) = fracrefa_var_167(ig_var_2323, jpl_var_2327) + fpl_var_2357 * (fracrefa_var_167(ig_var_2323, jpl_var_2327 + 1) - fracrefa_var_167(ig_var_2323, jpl_var_2327))
        END DO
      END DO
      ixc0_var_2367 = kfdia_var_2289 - kidia_var_2288 + 1 - ixc0_var_2367
      DO ixp_var_2368 = 1, ixc0_var_2367
        jl_var_2370 = ixhigh_var_2364(ixp_var_2368, lay_var_2325)
        speccomb_var_2316 = colo3_var_2303(jl_var_2370, lay_var_2325) + rat_o3co2_var_2311(jl_var_2370, lay_var_2325) * colco2_var_2302(jl_var_2370, lay_var_2325)
        specparm_var_2353 = MIN(colo3_var_2303(jl_var_2370, lay_var_2325) / speccomb_var_2316, oneminus_var_2300)
        specmult_var_2352 = 4.0D0 * (specparm_var_2353)
        js_var_2324 = 1 + INT(specmult_var_2352)
        fs_var_2351 = ((specmult_var_2352) - AINT((specmult_var_2352)))
        speccomb1_var_2317 = colo3_var_2303(jl_var_2370, lay_var_2325) + rat_o3co2_1_var_2312(jl_var_2370, lay_var_2325) * colco2_var_2302(jl_var_2370, lay_var_2325)
        specparm1_var_2356 = MIN(colo3_var_2303(jl_var_2370, lay_var_2325) / speccomb1_var_2317, oneminus_var_2300)
        specmult1_var_2355 = 4.0D0 * (specparm1_var_2356)
        js1_var_2326 = 1 + INT(specmult1_var_2355)
        fs1_var_2354 = ((specmult1_var_2355) - AINT((specmult1_var_2355)))
        fac000_var_2330 = (1.0D0 - fs_var_2351) * fac00_var_2293(jl_var_2370, lay_var_2325)
        fac010_var_2333 = (1.0D0 - fs_var_2351) * fac10_var_2295(jl_var_2370, lay_var_2325)
        fac100_var_2331 = fs_var_2351 * fac00_var_2293(jl_var_2370, lay_var_2325)
        fac110_var_2334 = fs_var_2351 * fac10_var_2295(jl_var_2370, lay_var_2325)
        fac001_var_2336 = (1.0D0 - fs1_var_2354) * fac01_var_2294(jl_var_2370, lay_var_2325)
        fac011_var_2339 = (1.0D0 - fs1_var_2354) * fac11_var_2296(jl_var_2370, lay_var_2325)
        fac101_var_2337 = fs1_var_2354 * fac01_var_2294(jl_var_2370, lay_var_2325)
        fac111_var_2340 = fs1_var_2354 * fac11_var_2296(jl_var_2370, lay_var_2325)
        speccomb_planck_var_2318 = colo3_var_2303(jl_var_2370, lay_var_2325) + refrat_planck_b_var_2329 * colco2_var_2302(jl_var_2370, lay_var_2325)
        specparm_planck_var_2359 = MIN(colo3_var_2303(jl_var_2370, lay_var_2325) / speccomb_planck_var_2318, oneminus_var_2300)
        specmult_planck_var_2358 = 4.0D0 * specparm_planck_var_2359
        jpl_var_2327 = 1 + INT(specmult_planck_var_2358)
        fpl_var_2357 = ((specmult_planck_var_2358) - AINT((specmult_planck_var_2358)))
        ind0_var_2319 = ((jp_var_2297(jl_var_2370, lay_var_2325) - 13) * 5 + (jt_var_2298(jl_var_2370, lay_var_2325) - 1)) * nspb_var_217(4) + js_var_2324
        ind1_var_2320 = ((jp_var_2297(jl_var_2370, lay_var_2325) - 12) * 5 + (jt1_var_2299(jl_var_2370, lay_var_2325) - 1)) * nspb_var_217(4) + js1_var_2326
        DO ig_var_2323 = 1, 14
          taug_var_2291(jl_var_2370, 38 + ig_var_2323, lay_var_2325) = speccomb_var_2316 * (fac000_var_2330 * absb_var_170(ind0_var_2319, ig_var_2323) + fac100_var_2331 * absb_var_170(ind0_var_2319 + 1, ig_var_2323) + fac010_var_2333 * absb_var_170(ind0_var_2319 + 5, ig_var_2323) + fac110_var_2334 * absb_var_170(ind0_var_2319 + 6, ig_var_2323)) + speccomb1_var_2317 * (fac001_var_2336 * absb_var_170(ind1_var_2320, ig_var_2323) + fac101_var_2337 * absb_var_170(ind1_var_2320 + 1, ig_var_2323) + fac011_var_2339 * absb_var_170(ind1_var_2320 + 5, ig_var_2323) + fac111_var_2340 * absb_var_170(ind1_var_2320 + 6, ig_var_2323))
          fracs_var_2308(jl_var_2370, 38 + ig_var_2323, lay_var_2325) = fracrefb_var_168(ig_var_2323, jpl_var_2327) + fpl_var_2357 * (fracrefb_var_168(ig_var_2323, jpl_var_2327 + 1) - fracrefb_var_168(ig_var_2323, jpl_var_2327))
        END DO
      END DO
      DO ixp_var_2368 = 1, ixc0_var_2367
        jl_var_2370 = ixhigh_var_2364(ixp_var_2368, lay_var_2325)
        taug_var_2291(jl_var_2370, 46, lay_var_2325) = taug_var_2291(jl_var_2370, 46, lay_var_2325) * 0.92D0
        taug_var_2291(jl_var_2370, 47, lay_var_2325) = taug_var_2291(jl_var_2370, 47, lay_var_2325) * 0.88D0
        taug_var_2291(jl_var_2370, 48, lay_var_2325) = taug_var_2291(jl_var_2370, 48, lay_var_2325) * 1.07D0
        taug_var_2291(jl_var_2370, 49, lay_var_2325) = taug_var_2291(jl_var_2370, 49, lay_var_2325) * 1.1D0
        taug_var_2291(jl_var_2370, 50, lay_var_2325) = taug_var_2291(jl_var_2370, 50, lay_var_2325) * 0.99D0
        taug_var_2291(jl_var_2370, 51, lay_var_2325) = taug_var_2291(jl_var_2370, 51, lay_var_2325) * 0.88D0
        taug_var_2291(jl_var_2370, 52, lay_var_2325) = taug_var_2291(jl_var_2370, 52, lay_var_2325) * 0.943D0
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol4
SUBROUTINE srtm_taumol16(kidia_var_2371, kfdia_var_2372, klev_var_2373, p_fac00_var_2374, p_fac01_var_2375, p_fac10_var_2376, p_fac11_var_2377, k_jp_var_2378, k_jt_var_2379, k_jt1_var_2380, p_oneminus_var_2381, p_colh2o_var_2382, p_colch4_var_2383, p_colmol_var_2384, k_laytrop_var_2385, p_selffac_var_2386, p_selffrac_var_2387, k_indself_var_2388, p_forfac_var_2389, p_forfrac_var_2390, k_indfor_var_2391, p_sfluxzen_var_2392, p_taug_var_2393, p_taur_var_2394, prmu0_var_2395)
  USE yoesrta16, ONLY: absa_var_220, absb_var_221, forrefc_var_223, layreffr_var_219, rayl_var_218, selfrefc_var_222, sfluxrefc_var_224, strrat1
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2371, kfdia_var_2372
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2373
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_2374(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_2375(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_2376(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_2377(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_2378(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_2379(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_2380(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_2381(kidia_var_2371 : kfdia_var_2372)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_2382(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_2383(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_2384(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_2385(kidia_var_2371 : kfdia_var_2372)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_2386(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_2387(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_2388(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_2389(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_2390(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_2391(kidia_var_2371 : kfdia_var_2372, klev_var_2373)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_2392(kidia_var_2371 : kfdia_var_2372, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_2393(kidia_var_2371 : kfdia_var_2372, klev_var_2373, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_2394(kidia_var_2371 : kfdia_var_2372, klev_var_2373, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_2395(kidia_var_2371 : kfdia_var_2372)
  INTEGER(KIND = 4) :: ig_var_2396, ind0_var_2397, ind1_var_2398, inds_var_2399, indf_var_2400, js_var_2401, i_lay_var_2402, i_laysolfr_var_2403(kidia_var_2371 : kfdia_var_2372), i_nlayers_var_2404, iplon_var_2405
  INTEGER(KIND = 4) :: laytrop_min_var_2406, laytrop_max_var_2407
  REAL(KIND = 8) :: z_fs_var_2408, z_speccomb_var_2409, z_specmult_var_2410, z_specparm_var_2411, z_tauray_var_2412
  laytrop_min_var_2406 = MINVAL(k_laytrop_var_2385(kidia_var_2371 : kfdia_var_2372))
  laytrop_max_var_2407 = MAXVAL(k_laytrop_var_2385(kidia_var_2371 : kfdia_var_2372))
  i_nlayers_var_2404 = klev_var_2373
  DO iplon_var_2405 = kidia_var_2371, kfdia_var_2372
    i_laysolfr_var_2403(iplon_var_2405) = i_nlayers_var_2404
  END DO
  DO i_lay_var_2402 = 1, laytrop_min_var_2406
    DO iplon_var_2405 = kidia_var_2371, kfdia_var_2372
      z_speccomb_var_2409 = p_colh2o_var_2382(iplon_var_2405, i_lay_var_2402) + strrat1 * p_colch4_var_2383(iplon_var_2405, i_lay_var_2402)
      z_specparm_var_2411 = p_colh2o_var_2382(iplon_var_2405, i_lay_var_2402) / z_speccomb_var_2409
      z_specparm_var_2411 = MIN(p_oneminus_var_2381(iplon_var_2405), z_specparm_var_2411)
      z_specmult_var_2410 = 8.0D0 * (z_specparm_var_2411)
      js_var_2401 = 1 + INT(z_specmult_var_2410)
      z_fs_var_2408 = z_specmult_var_2410 - AINT(z_specmult_var_2410)
      ind0_var_2397 = ((k_jp_var_2378(iplon_var_2405, i_lay_var_2402) - 1) * 5 + (k_jt_var_2379(iplon_var_2405, i_lay_var_2402) - 1)) * nspa_var_313(16) + js_var_2401
      ind1_var_2398 = (k_jp_var_2378(iplon_var_2405, i_lay_var_2402) * 5 + (k_jt1_var_2380(iplon_var_2405, i_lay_var_2402) - 1)) * nspa_var_313(16) + js_var_2401
      inds_var_2399 = k_indself_var_2388(iplon_var_2405, i_lay_var_2402)
      indf_var_2400 = k_indfor_var_2391(iplon_var_2405, i_lay_var_2402)
      z_tauray_var_2412 = p_colmol_var_2384(iplon_var_2405, i_lay_var_2402) * rayl_var_218
      DO ig_var_2396 = 1, 6
        p_taug_var_2393(iplon_var_2405, i_lay_var_2402, ig_var_2396) = z_speccomb_var_2409 * ((1.0D0 - z_fs_var_2408) * (absa_var_220(ind0_var_2397, ig_var_2396) * p_fac00_var_2374(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind0_var_2397 + 9, ig_var_2396) * p_fac10_var_2376(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind1_var_2398, ig_var_2396) * p_fac01_var_2375(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind1_var_2398 + 9, ig_var_2396) * p_fac11_var_2377(iplon_var_2405, i_lay_var_2402)) + z_fs_var_2408 * (absa_var_220(ind0_var_2397 + 1, ig_var_2396) * p_fac00_var_2374(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind0_var_2397 + 10, ig_var_2396) * p_fac10_var_2376(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind1_var_2398 + 1, ig_var_2396) * p_fac01_var_2375(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind1_var_2398 + 10, ig_var_2396) * p_fac11_var_2377(iplon_var_2405, i_lay_var_2402))) + p_colh2o_var_2382(iplon_var_2405, i_lay_var_2402) * (p_selffac_var_2386(iplon_var_2405, i_lay_var_2402) * (selfrefc_var_222(inds_var_2399, ig_var_2396) + p_selffrac_var_2387(iplon_var_2405, i_lay_var_2402) * (selfrefc_var_222(inds_var_2399 + 1, ig_var_2396) - selfrefc_var_222(inds_var_2399, ig_var_2396))) + p_forfac_var_2389(iplon_var_2405, i_lay_var_2402) * (forrefc_var_223(indf_var_2400, ig_var_2396) + p_forfrac_var_2390(iplon_var_2405, i_lay_var_2402) * (forrefc_var_223(indf_var_2400 + 1, ig_var_2396) - forrefc_var_223(indf_var_2400, ig_var_2396))))
        p_taur_var_2394(iplon_var_2405, i_lay_var_2402, ig_var_2396) = z_tauray_var_2412
      END DO
    END DO
  END DO
  DO i_lay_var_2402 = laytrop_min_var_2406 + 1, laytrop_max_var_2407
    DO iplon_var_2405 = kidia_var_2371, kfdia_var_2372
      IF (i_lay_var_2402 <= k_laytrop_var_2385(iplon_var_2405)) THEN
        z_speccomb_var_2409 = p_colh2o_var_2382(iplon_var_2405, i_lay_var_2402) + strrat1 * p_colch4_var_2383(iplon_var_2405, i_lay_var_2402)
        z_specparm_var_2411 = p_colh2o_var_2382(iplon_var_2405, i_lay_var_2402) / z_speccomb_var_2409
        z_specparm_var_2411 = MIN(p_oneminus_var_2381(iplon_var_2405), z_specparm_var_2411)
        z_specmult_var_2410 = 8.0D0 * (z_specparm_var_2411)
        js_var_2401 = 1 + INT(z_specmult_var_2410)
        z_fs_var_2408 = z_specmult_var_2410 - AINT(z_specmult_var_2410)
        ind0_var_2397 = ((k_jp_var_2378(iplon_var_2405, i_lay_var_2402) - 1) * 5 + (k_jt_var_2379(iplon_var_2405, i_lay_var_2402) - 1)) * nspa_var_313(16) + js_var_2401
        ind1_var_2398 = (k_jp_var_2378(iplon_var_2405, i_lay_var_2402) * 5 + (k_jt1_var_2380(iplon_var_2405, i_lay_var_2402) - 1)) * nspa_var_313(16) + js_var_2401
        inds_var_2399 = k_indself_var_2388(iplon_var_2405, i_lay_var_2402)
        indf_var_2400 = k_indfor_var_2391(iplon_var_2405, i_lay_var_2402)
        z_tauray_var_2412 = p_colmol_var_2384(iplon_var_2405, i_lay_var_2402) * rayl_var_218
        DO ig_var_2396 = 1, 6
          p_taug_var_2393(iplon_var_2405, i_lay_var_2402, ig_var_2396) = z_speccomb_var_2409 * ((1.0D0 - z_fs_var_2408) * (absa_var_220(ind0_var_2397, ig_var_2396) * p_fac00_var_2374(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind0_var_2397 + 9, ig_var_2396) * p_fac10_var_2376(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind1_var_2398, ig_var_2396) * p_fac01_var_2375(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind1_var_2398 + 9, ig_var_2396) * p_fac11_var_2377(iplon_var_2405, i_lay_var_2402)) + z_fs_var_2408 * (absa_var_220(ind0_var_2397 + 1, ig_var_2396) * p_fac00_var_2374(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind0_var_2397 + 10, ig_var_2396) * p_fac10_var_2376(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind1_var_2398 + 1, ig_var_2396) * p_fac01_var_2375(iplon_var_2405, i_lay_var_2402) + absa_var_220(ind1_var_2398 + 10, ig_var_2396) * p_fac11_var_2377(iplon_var_2405, i_lay_var_2402))) + p_colh2o_var_2382(iplon_var_2405, i_lay_var_2402) * (p_selffac_var_2386(iplon_var_2405, i_lay_var_2402) * (selfrefc_var_222(inds_var_2399, ig_var_2396) + p_selffrac_var_2387(iplon_var_2405, i_lay_var_2402) * (selfrefc_var_222(inds_var_2399 + 1, ig_var_2396) - selfrefc_var_222(inds_var_2399, ig_var_2396))) + p_forfac_var_2389(iplon_var_2405, i_lay_var_2402) * (forrefc_var_223(indf_var_2400, ig_var_2396) + p_forfrac_var_2390(iplon_var_2405, i_lay_var_2402) * (forrefc_var_223(indf_var_2400 + 1, ig_var_2396) - forrefc_var_223(indf_var_2400, ig_var_2396))))
          p_taur_var_2394(iplon_var_2405, i_lay_var_2402, ig_var_2396) = z_tauray_var_2412
        END DO
      ELSE
        IF (k_jp_var_2378(iplon_var_2405, i_lay_var_2402 - 1) < layreffr_var_219 .AND. k_jp_var_2378(iplon_var_2405, i_lay_var_2402) >= layreffr_var_219) i_laysolfr_var_2403(iplon_var_2405) = i_lay_var_2402
        ind0_var_2397 = ((k_jp_var_2378(iplon_var_2405, i_lay_var_2402) - 13) * 5 + (k_jt_var_2379(iplon_var_2405, i_lay_var_2402) - 1)) * nspb_var_314(16) + 1
        ind1_var_2398 = ((k_jp_var_2378(iplon_var_2405, i_lay_var_2402) - 12) * 5 + (k_jt1_var_2380(iplon_var_2405, i_lay_var_2402) - 1)) * nspb_var_314(16) + 1
        z_tauray_var_2412 = p_colmol_var_2384(iplon_var_2405, i_lay_var_2402) * rayl_var_218
        DO ig_var_2396 = 1, 6
          p_taug_var_2393(iplon_var_2405, i_lay_var_2402, ig_var_2396) = p_colch4_var_2383(iplon_var_2405, i_lay_var_2402) * (p_fac00_var_2374(iplon_var_2405, i_lay_var_2402) * absb_var_221(ind0_var_2397, ig_var_2396) + p_fac10_var_2376(iplon_var_2405, i_lay_var_2402) * absb_var_221(ind0_var_2397 + 1, ig_var_2396) + p_fac01_var_2375(iplon_var_2405, i_lay_var_2402) * absb_var_221(ind1_var_2398, ig_var_2396) + p_fac11_var_2377(iplon_var_2405, i_lay_var_2402) * absb_var_221(ind1_var_2398 + 1, ig_var_2396))
          IF (i_lay_var_2402 == i_laysolfr_var_2403(iplon_var_2405)) p_sfluxzen_var_2392(iplon_var_2405, ig_var_2396) = sfluxrefc_var_224(ig_var_2396)
          p_taur_var_2394(iplon_var_2405, i_lay_var_2402, ig_var_2396) = z_tauray_var_2412
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_2402 = laytrop_max_var_2407 + 1, i_nlayers_var_2404
    DO iplon_var_2405 = kidia_var_2371, kfdia_var_2372
      IF (k_jp_var_2378(iplon_var_2405, i_lay_var_2402 - 1) < layreffr_var_219 .AND. k_jp_var_2378(iplon_var_2405, i_lay_var_2402) >= layreffr_var_219) i_laysolfr_var_2403(iplon_var_2405) = i_lay_var_2402
      ind0_var_2397 = ((k_jp_var_2378(iplon_var_2405, i_lay_var_2402) - 13) * 5 + (k_jt_var_2379(iplon_var_2405, i_lay_var_2402) - 1)) * nspb_var_314(16) + 1
      ind1_var_2398 = ((k_jp_var_2378(iplon_var_2405, i_lay_var_2402) - 12) * 5 + (k_jt1_var_2380(iplon_var_2405, i_lay_var_2402) - 1)) * nspb_var_314(16) + 1
      z_tauray_var_2412 = p_colmol_var_2384(iplon_var_2405, i_lay_var_2402) * rayl_var_218
      DO ig_var_2396 = 1, 6
        p_taug_var_2393(iplon_var_2405, i_lay_var_2402, ig_var_2396) = p_colch4_var_2383(iplon_var_2405, i_lay_var_2402) * (p_fac00_var_2374(iplon_var_2405, i_lay_var_2402) * absb_var_221(ind0_var_2397, ig_var_2396) + p_fac10_var_2376(iplon_var_2405, i_lay_var_2402) * absb_var_221(ind0_var_2397 + 1, ig_var_2396) + p_fac01_var_2375(iplon_var_2405, i_lay_var_2402) * absb_var_221(ind1_var_2398, ig_var_2396) + p_fac11_var_2377(iplon_var_2405, i_lay_var_2402) * absb_var_221(ind1_var_2398 + 1, ig_var_2396))
        IF (i_lay_var_2402 == i_laysolfr_var_2403(iplon_var_2405)) p_sfluxzen_var_2392(iplon_var_2405, ig_var_2396) = sfluxrefc_var_224(ig_var_2396)
        p_taur_var_2394(iplon_var_2405, i_lay_var_2402, ig_var_2396) = z_tauray_var_2412
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol16
SUBROUTINE abor1(cdtext)
  USE yomlun_ifsaux, ONLY: nulerr, nulout
  IMPLICIT NONE
  CHARACTER(LEN = *), INTENT(IN) :: cdtext
  IF (nulout >= 0) WRITE(nulout, '(1x,a)') cdtext
  WRITE(nulerr, '(1x,a,a)') 'abort! ', cdtext
  IF (nulout >= 0) THEN
    IF (nulout /= 0 .AND. nulout /= 6) CLOSE(UNIT = nulout)
  END IF
  ERROR STOP 1
END SUBROUTINE abor1
SUBROUTINE global_init_fn
  USE yomlun_ifsaux, ONLY: nulout
  IMPLICIT NONE
  nulout = 6
END SUBROUTINE global_init_fn