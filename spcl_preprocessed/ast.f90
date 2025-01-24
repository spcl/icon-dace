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
    PRIVATE
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
    WRITE(nulerr, '(a)') text
    ERROR STOP 1
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
    sf = 1.0D0
    istatus = .TRUE.
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
    CALL gas_optics(ncol, nlev, istartcol, iendcol, config, single_level, thermodynamics, gas, od_lw_var_585, od_sw_var_591, ssa_sw_var_592, lw_albedo_var_461 = lw_albedo_var_599, planck_hl_var_465 = planck_hl_var_597, lw_emission_var_466 = lw_emission_var_598, incoming_sw_var_467 = incoming_sw_var_602)
    CALL crop_cloud_fraction(cloud, istartcol, iendcol, 1D-06, 1D-09)
    CALL cloud_optics_fn_438(nlev, istartcol, iendcol, config, thermodynamics, cloud, od_lw_cloud_var_588, ssa_lw_cloud_var_589, g_lw_cloud_var_590, od_sw_cloud_var_594, ssa_sw_cloud_var_595, g_sw_cloud_var_596)
    CALL add_aerosol_optics(nlev, istartcol, iendcol, config, thermodynamics, gas, aerosol, od_lw_var_585, ssa_lw_var_586, g_lw_var_587, od_sw_var_591, ssa_sw_var_592, g_sw_var_593)
    CALL solver_mcica_lw(nlev, istartcol, iendcol, config, single_level, cloud, od_lw_var_585, ssa_lw_var_586, g_lw_var_587, od_lw_cloud_var_588, ssa_lw_cloud_var_589, g_lw_cloud_var_590, planck_hl_var_597, lw_emission_var_598, lw_albedo_var_599, flux)
    CALL solver_mcica_sw(nlev, istartcol, iendcol, config, single_level, cloud, od_sw_var_591, ssa_sw_var_592, g_sw_var_593, od_sw_cloud_var_594, ssa_sw_cloud_var_595, g_sw_cloud_var_596, sw_albedo_direct_var_600, sw_albedo_diffuse_var_601, incoming_sw_var_602, flux)
    CALL calc_surface_spectral(flux, config, istartcol, iendcol)
    CALL calc_toa_spectral(flux, config, istartcol, iendcol)
  END SUBROUTINE radiation
END MODULE radiation_interface
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
SUBROUTINE srtm_setcoef(kidia_var_603, kfdia_var_604, klev_var_605, pavel_var_606, ptavel_var_607, pcoldry_var_608, pwkl_var_609, klaytrop_var_610, pcolch4_var_611, pcolco2_var_612, pcolh2o_var_613, pcolmol_var_614, pcolo2_var_615, pcolo3_var_616, pforfac_var_617, pforfrac_var_618, kindfor_var_619, pselffac_var_620, pselffrac_var_621, kindself_var_622, pfac00_var_623, pfac01_var_624, pfac10_var_625, pfac11_var_626, kjp_var_627, kjt_var_628, kjt1_var_629, prmu0_var_630)
  USE yoesrtwn, ONLY: preflog_var_315, tref_var_316
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_603, kfdia_var_604
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_605
  REAL(KIND = 8), INTENT(IN) :: pavel_var_606(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: ptavel_var_607(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: pcoldry_var_608(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: pwkl_var_609(kidia_var_603 : kfdia_var_604, 35, klev_var_605)
  INTEGER(KIND = 4), INTENT(INOUT) :: klaytrop_var_610(kidia_var_603 : kfdia_var_604)
  REAL(KIND = 8), INTENT(INOUT) :: pcolch4_var_611(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pcolco2_var_612(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pcolh2o_var_613(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pcolmol_var_614(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pcolo2_var_615(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pcolo3_var_616(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pforfac_var_617(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pforfrac_var_618(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(INOUT) :: kindfor_var_619(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pselffac_var_620(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pselffrac_var_621(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(INOUT) :: kindself_var_622(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pfac00_var_623(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pfac01_var_624(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pfac10_var_625(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(INOUT) :: pfac11_var_626(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjp_var_627(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjt_var_628(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjt1_var_629(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_630(kidia_var_603 : kfdia_var_604)
  INTEGER(KIND = 4) :: i_nlayers_var_631, jk_var_632, jl_var_633, jp1_var_634
  REAL(KIND = 8) :: z_stpfac_var_635, z_plog_var_636
  REAL(KIND = 8) :: z_fp_var_637, z_ft_var_638, z_ft1_var_639, z_water_var_640, z_scalefac_var_641
  REAL(KIND = 8) :: z_factor_var_642, z_co2reg_var_643, z_compfp_var_644
  z_stpfac_var_635 = 0.29220138203356366D0
  i_nlayers_var_631 = klev_var_605
  DO jk_var_632 = 1, klev_var_605
    DO jl_var_633 = kidia_var_603, kfdia_var_604
      pcolmol_var_614(jl_var_633, jk_var_632) = 0.0D0
    END DO
  END DO
  DO jl_var_633 = kidia_var_603, kfdia_var_604
    IF (prmu0_var_630(jl_var_633) > 0.0D0) THEN
      klaytrop_var_610(jl_var_633) = 0
    END IF
  END DO
  DO jk_var_632 = 1, i_nlayers_var_631
    DO jl_var_633 = kidia_var_603, kfdia_var_604
      IF (prmu0_var_630(jl_var_633) > 0.0D0) THEN
        z_plog_var_636 = LOG(pavel_var_606(jl_var_633, jk_var_632))
        kjp_var_627(jl_var_633, jk_var_632) = INT(36.0D0 - 5.0D0 * (z_plog_var_636 + 0.04D0))
        IF (kjp_var_627(jl_var_633, jk_var_632) < 1) THEN
          kjp_var_627(jl_var_633, jk_var_632) = 1
        ELSE IF (kjp_var_627(jl_var_633, jk_var_632) > 58) THEN
          kjp_var_627(jl_var_633, jk_var_632) = 58
        END IF
        jp1_var_634 = kjp_var_627(jl_var_633, jk_var_632) + 1
        z_fp_var_637 = 5.0 * (preflog_var_315(kjp_var_627(jl_var_633, jk_var_632)) - z_plog_var_636)
        kjt_var_628(jl_var_633, jk_var_632) = INT(3.0 + (ptavel_var_607(jl_var_633, jk_var_632) - tref_var_316(kjp_var_627(jl_var_633, jk_var_632))) / 15.0)
        IF (kjt_var_628(jl_var_633, jk_var_632) < 1) THEN
          kjt_var_628(jl_var_633, jk_var_632) = 1
        ELSE IF (kjt_var_628(jl_var_633, jk_var_632) > 4) THEN
          kjt_var_628(jl_var_633, jk_var_632) = 4
        END IF
        z_ft_var_638 = ((ptavel_var_607(jl_var_633, jk_var_632) - tref_var_316(kjp_var_627(jl_var_633, jk_var_632))) / 15.0) - REAL(kjt_var_628(jl_var_633, jk_var_632) - 3)
        kjt1_var_629(jl_var_633, jk_var_632) = INT(3.0 + (ptavel_var_607(jl_var_633, jk_var_632) - tref_var_316(jp1_var_634)) / 15.0)
        IF (kjt1_var_629(jl_var_633, jk_var_632) < 1) THEN
          kjt1_var_629(jl_var_633, jk_var_632) = 1
        ELSE IF (kjt1_var_629(jl_var_633, jk_var_632) > 4) THEN
          kjt1_var_629(jl_var_633, jk_var_632) = 4
        END IF
        z_ft1_var_639 = ((ptavel_var_607(jl_var_633, jk_var_632) - tref_var_316(jp1_var_634)) / 15.0) - REAL(kjt1_var_629(jl_var_633, jk_var_632) - 3)
        z_water_var_640 = pwkl_var_609(jl_var_633, 1, jk_var_632) / pcoldry_var_608(jl_var_633, jk_var_632)
        z_scalefac_var_641 = pavel_var_606(jl_var_633, jk_var_632) * z_stpfac_var_635 / ptavel_var_607(jl_var_633, jk_var_632)
        IF (z_plog_var_636 <= 4.56D0) GO TO 5300
        klaytrop_var_610(jl_var_633) = klaytrop_var_610(jl_var_633) + 1
        pforfac_var_617(jl_var_633, jk_var_632) = z_scalefac_var_641 / (1.0 + z_water_var_640)
        z_factor_var_642 = (332.0 - ptavel_var_607(jl_var_633, jk_var_632)) / 36.0
        kindfor_var_619(jl_var_633, jk_var_632) = MIN(2, MAX(1, INT(z_factor_var_642)))
        pforfrac_var_618(jl_var_633, jk_var_632) = z_factor_var_642 - REAL(kindfor_var_619(jl_var_633, jk_var_632))
        pselffac_var_620(jl_var_633, jk_var_632) = z_water_var_640 * pforfac_var_617(jl_var_633, jk_var_632)
        z_factor_var_642 = (ptavel_var_607(jl_var_633, jk_var_632) - 188.0) / 7.2
        kindself_var_622(jl_var_633, jk_var_632) = MIN(9, MAX(1, INT(z_factor_var_642) - 7))
        pselffrac_var_621(jl_var_633, jk_var_632) = z_factor_var_642 - REAL(kindself_var_622(jl_var_633, jk_var_632) + 7)
        pcolh2o_var_613(jl_var_633, jk_var_632) = 1E-20 * pwkl_var_609(jl_var_633, 1, jk_var_632)
        pcolco2_var_612(jl_var_633, jk_var_632) = 1E-20 * pwkl_var_609(jl_var_633, 2, jk_var_632)
        pcolo3_var_616(jl_var_633, jk_var_632) = 1E-20 * pwkl_var_609(jl_var_633, 3, jk_var_632)
        pcolch4_var_611(jl_var_633, jk_var_632) = 1E-20 * pwkl_var_609(jl_var_633, 6, jk_var_632)
        pcolo2_var_615(jl_var_633, jk_var_632) = 1E-20 * pwkl_var_609(jl_var_633, 7, jk_var_632)
        pcolmol_var_614(jl_var_633, jk_var_632) = 1E-20 * pcoldry_var_608(jl_var_633, jk_var_632) + pcolh2o_var_613(jl_var_633, jk_var_632)
        IF (pcolco2_var_612(jl_var_633, jk_var_632) == 0.0) pcolco2_var_612(jl_var_633, jk_var_632) = 1E-32 * pcoldry_var_608(jl_var_633, jk_var_632)
        IF (pcolch4_var_611(jl_var_633, jk_var_632) == 0.0) pcolch4_var_611(jl_var_633, jk_var_632) = 1E-32 * pcoldry_var_608(jl_var_633, jk_var_632)
        IF (pcolo2_var_615(jl_var_633, jk_var_632) == 0.0) pcolo2_var_615(jl_var_633, jk_var_632) = 1E-32 * pcoldry_var_608(jl_var_633, jk_var_632)
        z_co2reg_var_643 = 3.55E-24 * pcoldry_var_608(jl_var_633, jk_var_632)
        GO TO 5400
5300    CONTINUE
        pforfac_var_617(jl_var_633, jk_var_632) = z_scalefac_var_641 / (1.0 + z_water_var_640)
        z_factor_var_642 = (ptavel_var_607(jl_var_633, jk_var_632) - 188.0) / 36.0
        kindfor_var_619(jl_var_633, jk_var_632) = 3
        pforfrac_var_618(jl_var_633, jk_var_632) = z_factor_var_642 - 1.0
        pcolh2o_var_613(jl_var_633, jk_var_632) = 1E-20 * pwkl_var_609(jl_var_633, 1, jk_var_632)
        pcolco2_var_612(jl_var_633, jk_var_632) = 1E-20 * pwkl_var_609(jl_var_633, 2, jk_var_632)
        pcolo3_var_616(jl_var_633, jk_var_632) = 1E-20 * pwkl_var_609(jl_var_633, 3, jk_var_632)
        pcolch4_var_611(jl_var_633, jk_var_632) = 1E-20 * pwkl_var_609(jl_var_633, 6, jk_var_632)
        pcolo2_var_615(jl_var_633, jk_var_632) = 1E-20 * pwkl_var_609(jl_var_633, 7, jk_var_632)
        pcolmol_var_614(jl_var_633, jk_var_632) = 1E-20 * pcoldry_var_608(jl_var_633, jk_var_632) + pcolh2o_var_613(jl_var_633, jk_var_632)
        IF (pcolco2_var_612(jl_var_633, jk_var_632) == 0.0) pcolco2_var_612(jl_var_633, jk_var_632) = 1E-32 * pcoldry_var_608(jl_var_633, jk_var_632)
        IF (pcolch4_var_611(jl_var_633, jk_var_632) == 0.0) pcolch4_var_611(jl_var_633, jk_var_632) = 1E-32 * pcoldry_var_608(jl_var_633, jk_var_632)
        IF (pcolo2_var_615(jl_var_633, jk_var_632) == 0.0) pcolo2_var_615(jl_var_633, jk_var_632) = 1E-32 * pcoldry_var_608(jl_var_633, jk_var_632)
        z_co2reg_var_643 = 3.55E-24 * pcoldry_var_608(jl_var_633, jk_var_632)
        pselffac_var_620(jl_var_633, jk_var_632) = 0.0D0
        pselffrac_var_621(jl_var_633, jk_var_632) = 0.0D0
        kindself_var_622(jl_var_633, jk_var_632) = 0
5400    CONTINUE
        z_compfp_var_644 = 1.0 - z_fp_var_637
        pfac10_var_625(jl_var_633, jk_var_632) = z_compfp_var_644 * z_ft_var_638
        pfac00_var_623(jl_var_633, jk_var_632) = z_compfp_var_644 * (1.0 - z_ft_var_638)
        pfac11_var_626(jl_var_633, jk_var_632) = z_fp_var_637 * z_ft1_var_639
        pfac01_var_624(jl_var_633, jk_var_632) = z_fp_var_637 * (1.0 - z_ft1_var_639)
9000    FORMAT(1X, 2I3, 3I4, F6.1, 4F7.2, 12E9.2, 2I5)
      END IF
    END DO
  END DO
END SUBROUTINE srtm_setcoef
SUBROUTINE rrtm_gas_optical_depth(kidia_var_645, kfdia_var_646, klev_var_647, pod_var_648, pavel_var_649, pcoldry_var_650, pcolbrd, pwx_var_651, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, poneminus_var_659, pcolh2o_var_660, pcolco2_var_661, pcolo3_var_662, pcoln2o, pcolch4_var_663, pcolo2_var_664, p_co2mult_var_665, klaytrop_var_666, klayswtch, klaylow, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, kindminor, pscaleminor, pscaleminorn2, pminorfrac, prat_h2oco2_var_674, prat_h2oco2_1_var_675, prat_h2oo3_var_676, prat_h2oo3_1_var_677, prat_h2on2o_var_678, prat_h2on2o_1_var_679, prat_h2och4_var_680, prat_h2och4_1_var_681, prat_n2oco2_var_682, prat_n2oco2_1_var_683, prat_o3co2_var_684, prat_o3co2_1_var_685)
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_645
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_646
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_647
  REAL(KIND = 8), INTENT(OUT) :: pod_var_648(140, klev_var_647, kidia_var_645 : kfdia_var_646)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_649(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pcoldry_var_650(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pwx_var_651(kidia_var_645 : kfdia_var_646, 4, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: ptauaerl(kidia_var_645 : kfdia_var_646, klev_var_647, 16)
  REAL(KIND = 8), INTENT(IN) :: pfac00_var_652(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pfac01_var_653(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pfac10_var_654(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pfac11_var_655(kidia_var_645 : kfdia_var_646, klev_var_647)
  INTEGER(KIND = 4), INTENT(IN) :: kjp_var_656(kidia_var_645 : kfdia_var_646, klev_var_647)
  INTEGER(KIND = 4), INTENT(IN) :: kjt_var_657(kidia_var_645 : kfdia_var_646, klev_var_647)
  INTEGER(KIND = 4), INTENT(IN) :: kjt1_var_658(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: poneminus_var_659
  REAL(KIND = 8), INTENT(IN) :: pcolh2o_var_660(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pcolco2_var_661(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pcolo3_var_662(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pcoln2o(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pcolch4_var_663(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pcolo2_var_664(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: p_co2mult_var_665(kidia_var_645 : kfdia_var_646, klev_var_647)
  INTEGER(KIND = 4), INTENT(IN) :: klaytrop_var_666(kidia_var_645 : kfdia_var_646)
  INTEGER(KIND = 4), INTENT(IN) :: klayswtch(kidia_var_645 : kfdia_var_646)
  INTEGER(KIND = 4), INTENT(IN) :: klaylow(kidia_var_645 : kfdia_var_646)
  REAL(KIND = 8), INTENT(IN) :: pselffac_var_667(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pselffrac_var_668(kidia_var_645 : kfdia_var_646, klev_var_647)
  INTEGER(KIND = 4), INTENT(IN) :: kindself_var_669(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(OUT) :: pfrac_var_670(kidia_var_645 : kfdia_var_646, 140, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pforfac_var_671(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pforfrac_var_672(kidia_var_645 : kfdia_var_646, klev_var_647)
  INTEGER(KIND = 4), INTENT(IN) :: kindfor_var_673(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pminorfrac(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pscaleminor(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pscaleminorn2(kidia_var_645 : kfdia_var_646, klev_var_647)
  INTEGER(KIND = 4), INTENT(IN) :: kindminor(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: pcolbrd(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8), INTENT(IN) :: prat_h2oco2_var_674(kidia_var_645 : kfdia_var_646, klev_var_647), prat_h2oco2_1_var_675(kidia_var_645 : kfdia_var_646, klev_var_647), prat_h2oo3_var_676(kidia_var_645 : kfdia_var_646, klev_var_647), prat_h2oo3_1_var_677(kidia_var_645 : kfdia_var_646, klev_var_647), prat_h2on2o_var_678(kidia_var_645 : kfdia_var_646, klev_var_647), prat_h2on2o_1_var_679(kidia_var_645 : kfdia_var_646, klev_var_647), prat_h2och4_var_680(kidia_var_645 : kfdia_var_646, klev_var_647), prat_h2och4_1_var_681(kidia_var_645 : kfdia_var_646, klev_var_647), prat_n2oco2_var_682(kidia_var_645 : kfdia_var_646, klev_var_647), prat_n2oco2_1_var_683(kidia_var_645 : kfdia_var_646, klev_var_647), prat_o3co2_var_684(kidia_var_645 : kfdia_var_646, klev_var_647), prat_o3co2_1_var_685(kidia_var_645 : kfdia_var_646, klev_var_647)
  REAL(KIND = 8) :: ztau(kidia_var_645 : kfdia_var_646, 140, klev_var_647)
  INTEGER(KIND = 4) :: ji, jlev_var_686
  INTEGER(KIND = 4) :: jlon_var_687
  pfrac_var_670(:, :, :) = 0.0D0
  CALL rrtm_taumol1(kidia_var_645, kfdia_var_646, klev_var_647, ztau, pavel_var_649, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, pcolh2o_var_660, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, pminorfrac, kindminor, pscaleminorn2, pcolbrd)
  CALL rrtm_taumol2(kidia_var_645, kfdia_var_646, klev_var_647, ztau, pavel_var_649, pcoldry_var_650, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, pcolh2o_var_660, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670)
  CALL rrtm_taumol3(kidia_var_645, kfdia_var_646, klev_var_647, ztau, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, poneminus_var_659, pcolh2o_var_660, pcolco2_var_661, pcoln2o, pcoldry_var_650, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, prat_h2oco2_var_674, prat_h2oco2_1_var_675, pminorfrac, kindminor)
  CALL rrtm_taumol4(kidia_var_645, kfdia_var_646, klev_var_647, ztau, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, poneminus_var_659, pcolh2o_var_660, pcolco2_var_661, pcolo3_var_662, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, prat_h2oco2_var_674, prat_h2oco2_1_var_675, prat_o3co2_var_684, prat_o3co2_1_var_685)
  CALL rrtm_taumol5(kidia_var_645, kfdia_var_646, klev_var_647, ztau, pwx_var_651, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, poneminus_var_659, pcolh2o_var_660, pcolco2_var_661, pcolo3_var_662, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, prat_h2oco2_var_674, prat_h2oco2_1_var_675, prat_o3co2_var_684, prat_o3co2_1_var_685, pminorfrac, kindminor)
  CALL rrtm_taumol6(kidia_var_645, kfdia_var_646, klev_var_647, ztau, pwx_var_651, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, pcolh2o_var_660, pcolco2_var_661, pcoldry_var_650, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, pminorfrac, kindminor)
  CALL rrtm_taumol7(kidia_var_645, kfdia_var_646, klev_var_647, ztau, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, poneminus_var_659, pcolh2o_var_660, pcolo3_var_662, pcolco2_var_661, pcoldry_var_650, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, prat_h2oo3_var_676, prat_h2oo3_1_var_677, pminorfrac, kindminor)
  CALL rrtm_taumol8(kidia_var_645, kfdia_var_646, klev_var_647, ztau, pwx_var_651, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, pcolh2o_var_660, pcolo3_var_662, pcoln2o, pcolco2_var_661, pcoldry_var_650, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, pminorfrac, kindminor)
  CALL rrtm_taumol9(kidia_var_645, kfdia_var_646, klev_var_647, ztau, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, poneminus_var_659, pcolh2o_var_660, pcoln2o, pcolch4_var_663, pcoldry_var_650, klaytrop_var_666, klayswtch, klaylow, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, prat_h2och4_var_680, prat_h2och4_1_var_681, pminorfrac, kindminor)
  CALL rrtm_taumol10(kidia_var_645, kfdia_var_646, klev_var_647, ztau, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, pcolh2o_var_660, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670)
  CALL rrtm_taumol11(kidia_var_645, kfdia_var_646, klev_var_647, ztau, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, pcolh2o_var_660, pcolo2_var_664, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, pminorfrac, kindminor, pscaleminor)
  CALL rrtm_taumol12(kidia_var_645, kfdia_var_646, klev_var_647, ztau, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, poneminus_var_659, pcolh2o_var_660, pcolco2_var_661, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, prat_h2oco2_var_674, prat_h2oco2_1_var_675)
  CALL rrtm_taumol13(kidia_var_645, kfdia_var_646, klev_var_647, ztau, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, poneminus_var_659, pcolh2o_var_660, pcoln2o, pcolco2_var_661, pcolo3_var_662, pcoldry_var_650, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, prat_h2on2o_var_678, prat_h2on2o_1_var_679, pminorfrac, kindminor)
  CALL rrtm_taumol14(kidia_var_645, kfdia_var_646, klev_var_647, ztau, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, pcolco2_var_661, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670)
  CALL rrtm_taumol15(kidia_var_645, kfdia_var_646, klev_var_647, ztau, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, poneminus_var_659, pcolh2o_var_660, pcolco2_var_661, pcoln2o, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, prat_n2oco2_var_682, prat_n2oco2_1_var_683, pminorfrac, kindminor, pscaleminor, pcolbrd)
  CALL rrtm_taumol16(kidia_var_645, kfdia_var_646, klev_var_647, ztau, ptauaerl, pfac00_var_652, pfac01_var_653, pfac10_var_654, pfac11_var_655, pforfac_var_671, pforfrac_var_672, kindfor_var_673, kjp_var_656, kjt_var_657, kjt1_var_658, poneminus_var_659, pcolh2o_var_660, pcolch4_var_663, klaytrop_var_666, pselffac_var_667, pselffrac_var_668, kindself_var_669, pfrac_var_670, prat_h2och4_var_680, prat_h2och4_1_var_681)
  DO jlev_var_686 = 1, klev_var_647
    DO ji = 1, 140
      DO jlon_var_687 = kidia_var_645, kfdia_var_646
        pod_var_648(ji, jlev_var_686, jlon_var_687) = ztau(jlon_var_687, ji, jlev_var_686)
      END DO
    END DO
  END DO
END SUBROUTINE rrtm_gas_optical_depth
SUBROUTINE rrtm_taumol9(kidia_var_688, kfdia_var_689, klev_var_690, taug_var_691, p_tauaerl_var_692, fac00_var_693, fac01_var_694, fac10_var_695, fac11_var_696, forfac_var_715, forfrac_var_716, indfor_var_714, jp_var_697, jt_var_698, jt1_var_699, oneminus_var_700, colh2o_var_701, coln2o_var_702, colch4_var_703, coldry_var_704, laytrop_var_705, k_layswtch_var_706, k_laylow_var_707, selffac_var_708, selffrac_var_709, indself_var_710, fracs_var_711, rat_h2och4_var_712, rat_h2och4_1_var_713, minorfrac_var_717, indminor_var_718)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng9
  USE yoerrta9, ONLY: absa_var_208, absb_var_209, forref_var_213, fracrefa_var_206, fracrefb_var_207, ka_mn2o_var_210, kb_mn2o_var_211, selfref_var_212
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_688
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_689
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_690
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_691(kidia_var_688 : kfdia_var_689, 140, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_692(kidia_var_688 : kfdia_var_689, klev_var_690, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_693(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_694(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_695(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_696(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_697(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_698(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_699(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_700
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_701(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_702(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: colch4_var_703(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_704(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_705(kidia_var_688 : kfdia_var_689)
  INTEGER(KIND = 4), INTENT(IN) :: k_layswtch_var_706(kidia_var_688 : kfdia_var_689)
  INTEGER(KIND = 4), INTENT(IN) :: k_laylow_var_707(kidia_var_688 : kfdia_var_689)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_708(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_709(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_710(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_711(kidia_var_688 : kfdia_var_689, 140, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_var_712(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_1_var_713(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_714(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_715(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_716(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_717(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_718(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4) :: ind0_var_719, ind1_var_720, inds_var_721, indf_var_722, indm_var_723
  INTEGER(KIND = 4) :: ig_var_724, js_var_725, lay_var_726, js1_var_727, jmn2o_var_728, jpl_var_729
  REAL(KIND = 8) :: speccomb_var_730, speccomb1_var_731, speccomb_mn2o_var_732, speccomb_planck_var_733
  REAL(KIND = 8) :: refrat_planck_a_var_734, refrat_m_a_var_735
  REAL(KIND = 8) :: fs_var_736, specmult_var_737, specparm_var_738, fs1_var_739, specmult1_var_740, specparm1_var_741, fmn2o_var_742, specmult_mn2o_var_743, specparm_mn2o_var_744, fpl_var_745, specmult_planck_var_746, specparm_planck_var_747
  REAL(KIND = 8) :: adjfac_var_748, adjcoln2o_var_749, ratn2o_var_750, chi_n2o_var_751
  REAL(KIND = 8) :: fac000_var_752, fac100_var_753, fac200_var_754, fac010_var_755, fac110_var_756, fac210_var_757, fac001_var_758, fac101_var_759, fac201_var_760, fac011_var_761, fac111_var_762, fac211_var_763
  REAL(KIND = 8) :: p_var_764, p4_var_765, fk0_var_766, fk1_var_767, fk2_var_768
  REAL(KIND = 8) :: taufor_var_769, tauself_var_770, n2om1_var_771, n2om2_var_772, absn2o_var_773, tau_major_var_774(12), tau_major1_var_775(12)
  INTEGER(KIND = 4) :: laytrop_min_var_776, laytrop_max_var_777
  INTEGER(KIND = 4) :: ixc_var_778(klev_var_690), ixlow_var_779(kfdia_var_689, klev_var_690), ixhigh_var_780(kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4) :: ich_var_781, icl_var_782, ixc0_var_783, ixp_var_784, jc_var_785, jl_var_786
  laytrop_min_var_776 = MINVAL(laytrop_var_705)
  laytrop_max_var_777 = MAXVAL(laytrop_var_705)
  ixlow_var_779 = 0
  ixhigh_var_780 = 0
  ixc_var_778 = 0
  DO lay_var_726 = laytrop_min_var_776 + 1, laytrop_max_var_777
    icl_var_782 = 0
    ich_var_781 = 0
    DO jc_var_785 = kidia_var_688, kfdia_var_689
      IF (lay_var_726 <= laytrop_var_705(jc_var_785)) THEN
        icl_var_782 = icl_var_782 + 1
        ixlow_var_779(icl_var_782, lay_var_726) = jc_var_785
      ELSE
        ich_var_781 = ich_var_781 + 1
        ixhigh_var_780(ich_var_781, lay_var_726) = jc_var_785
      END IF
    END DO
    ixc_var_778(lay_var_726) = icl_var_782
  END DO
  refrat_planck_a_var_734 = chi_mls(1, 9) / chi_mls(6, 9)
  refrat_m_a_var_735 = chi_mls(1, 3) / chi_mls(6, 3)
  DO lay_var_726 = 1, laytrop_min_var_776
    DO jl_var_786 = kidia_var_688, kfdia_var_689
      speccomb_var_730 = colh2o_var_701(jl_var_786, lay_var_726) + rat_h2och4_var_712(jl_var_786, lay_var_726) * colch4_var_703(jl_var_786, lay_var_726)
      specparm_var_738 = MIN(colh2o_var_701(jl_var_786, lay_var_726) / speccomb_var_730, oneminus_var_700)
      specmult_var_737 = 8.0D0 * (specparm_var_738)
      js_var_725 = 1 + INT(specmult_var_737)
      fs_var_736 = ((specmult_var_737) - AINT((specmult_var_737)))
      speccomb1_var_731 = colh2o_var_701(jl_var_786, lay_var_726) + rat_h2och4_1_var_713(jl_var_786, lay_var_726) * colch4_var_703(jl_var_786, lay_var_726)
      specparm1_var_741 = MIN(colh2o_var_701(jl_var_786, lay_var_726) / speccomb1_var_731, oneminus_var_700)
      specmult1_var_740 = 8.0D0 * (specparm1_var_741)
      js1_var_727 = 1 + INT(specmult1_var_740)
      fs1_var_739 = ((specmult1_var_740) - AINT((specmult1_var_740)))
      speccomb_mn2o_var_732 = colh2o_var_701(jl_var_786, lay_var_726) + refrat_m_a_var_735 * colch4_var_703(jl_var_786, lay_var_726)
      specparm_mn2o_var_744 = MIN(colh2o_var_701(jl_var_786, lay_var_726) / speccomb_mn2o_var_732, oneminus_var_700)
      specmult_mn2o_var_743 = 8.0D0 * specparm_mn2o_var_744
      jmn2o_var_728 = 1 + INT(specmult_mn2o_var_743)
      fmn2o_var_742 = ((specmult_mn2o_var_743) - AINT((specmult_mn2o_var_743)))
      chi_n2o_var_751 = coln2o_var_702(jl_var_786, lay_var_726) / (coldry_var_704(jl_var_786, lay_var_726))
      ratn2o_var_750 = 1D+20 * chi_n2o_var_751 / chi_mls(4, jp_var_697(jl_var_786, lay_var_726) + 1)
      IF (ratn2o_var_750 .GT. 1.5D0) THEN
        adjfac_var_748 = 0.5D0 + (ratn2o_var_750 - 0.5D0) ** 0.65D0
        adjcoln2o_var_749 = adjfac_var_748 * chi_mls(4, jp_var_697(jl_var_786, lay_var_726) + 1) * coldry_var_704(jl_var_786, lay_var_726) * 1D-20
      ELSE
        adjcoln2o_var_749 = coln2o_var_702(jl_var_786, lay_var_726)
      END IF
      speccomb_planck_var_733 = colh2o_var_701(jl_var_786, lay_var_726) + refrat_planck_a_var_734 * colch4_var_703(jl_var_786, lay_var_726)
      specparm_planck_var_747 = MIN(colh2o_var_701(jl_var_786, lay_var_726) / speccomb_planck_var_733, oneminus_var_700)
      specmult_planck_var_746 = 8.0D0 * specparm_planck_var_747
      jpl_var_729 = 1 + INT(specmult_planck_var_746)
      fpl_var_745 = ((specmult_planck_var_746) - AINT((specmult_planck_var_746)))
      ind0_var_719 = ((jp_var_697(jl_var_786, lay_var_726) - 1) * 5 + (jt_var_698(jl_var_786, lay_var_726) - 1)) * nspa_var_216(9) + js_var_725
      ind1_var_720 = (jp_var_697(jl_var_786, lay_var_726) * 5 + (jt1_var_699(jl_var_786, lay_var_726) - 1)) * nspa_var_216(9) + js1_var_727
      inds_var_721 = indself_var_710(jl_var_786, lay_var_726)
      indf_var_722 = indfor_var_714(jl_var_786, lay_var_726)
      indm_var_723 = indminor_var_718(jl_var_786, lay_var_726)
      IF (specparm_var_738 .LT. 0.125D0) THEN
        p_var_764 = fs_var_736 - 1.0D0
        p4_var_765 = p_var_764 ** 4
        fk0_var_766 = p4_var_765
        fk1_var_767 = 1.0D0 - p_var_764 - 2.0D0 * p4_var_765
        fk2_var_768 = p_var_764 + p4_var_765
        fac000_var_752 = fk0_var_766 * fac00_var_693(jl_var_786, lay_var_726)
        fac100_var_753 = fk1_var_767 * fac00_var_693(jl_var_786, lay_var_726)
        fac200_var_754 = fk2_var_768 * fac00_var_693(jl_var_786, lay_var_726)
        fac010_var_755 = fk0_var_766 * fac10_var_695(jl_var_786, lay_var_726)
        fac110_var_756 = fk1_var_767 * fac10_var_695(jl_var_786, lay_var_726)
        fac210_var_757 = fk2_var_768 * fac10_var_695(jl_var_786, lay_var_726)
      ELSE IF (specparm_var_738 .GT. 0.875D0) THEN
        p_var_764 = - fs_var_736
        p4_var_765 = p_var_764 ** 4
        fk0_var_766 = p4_var_765
        fk1_var_767 = 1.0D0 - p_var_764 - 2.0D0 * p4_var_765
        fk2_var_768 = p_var_764 + p4_var_765
        fac000_var_752 = fk0_var_766 * fac00_var_693(jl_var_786, lay_var_726)
        fac100_var_753 = fk1_var_767 * fac00_var_693(jl_var_786, lay_var_726)
        fac200_var_754 = fk2_var_768 * fac00_var_693(jl_var_786, lay_var_726)
        fac010_var_755 = fk0_var_766 * fac10_var_695(jl_var_786, lay_var_726)
        fac110_var_756 = fk1_var_767 * fac10_var_695(jl_var_786, lay_var_726)
        fac210_var_757 = fk2_var_768 * fac10_var_695(jl_var_786, lay_var_726)
      ELSE
        fac000_var_752 = (1.0D0 - fs_var_736) * fac00_var_693(jl_var_786, lay_var_726)
        fac010_var_755 = (1.0D0 - fs_var_736) * fac10_var_695(jl_var_786, lay_var_726)
        fac100_var_753 = fs_var_736 * fac00_var_693(jl_var_786, lay_var_726)
        fac110_var_756 = fs_var_736 * fac10_var_695(jl_var_786, lay_var_726)
        fac200_var_754 = 0.0D0
        fac210_var_757 = 0.0D0
      END IF
      IF (specparm1_var_741 .LT. 0.125D0) THEN
        p_var_764 = fs1_var_739 - 1.0D0
        p4_var_765 = p_var_764 ** 4
        fk0_var_766 = p4_var_765
        fk1_var_767 = 1.0D0 - p_var_764 - 2.0D0 * p4_var_765
        fk2_var_768 = p_var_764 + p4_var_765
        fac001_var_758 = fk0_var_766 * fac01_var_694(jl_var_786, lay_var_726)
        fac101_var_759 = fk1_var_767 * fac01_var_694(jl_var_786, lay_var_726)
        fac201_var_760 = fk2_var_768 * fac01_var_694(jl_var_786, lay_var_726)
        fac011_var_761 = fk0_var_766 * fac11_var_696(jl_var_786, lay_var_726)
        fac111_var_762 = fk1_var_767 * fac11_var_696(jl_var_786, lay_var_726)
        fac211_var_763 = fk2_var_768 * fac11_var_696(jl_var_786, lay_var_726)
      ELSE IF (specparm1_var_741 .GT. 0.875D0) THEN
        p_var_764 = - fs1_var_739
        p4_var_765 = p_var_764 ** 4
        fk0_var_766 = p4_var_765
        fk1_var_767 = 1.0D0 - p_var_764 - 2.0D0 * p4_var_765
        fk2_var_768 = p_var_764 + p4_var_765
        fac001_var_758 = fk0_var_766 * fac01_var_694(jl_var_786, lay_var_726)
        fac101_var_759 = fk1_var_767 * fac01_var_694(jl_var_786, lay_var_726)
        fac201_var_760 = fk2_var_768 * fac01_var_694(jl_var_786, lay_var_726)
        fac011_var_761 = fk0_var_766 * fac11_var_696(jl_var_786, lay_var_726)
        fac111_var_762 = fk1_var_767 * fac11_var_696(jl_var_786, lay_var_726)
        fac211_var_763 = fk2_var_768 * fac11_var_696(jl_var_786, lay_var_726)
      ELSE
        fac001_var_758 = (1.0D0 - fs1_var_739) * fac01_var_694(jl_var_786, lay_var_726)
        fac011_var_761 = (1.0D0 - fs1_var_739) * fac11_var_696(jl_var_786, lay_var_726)
        fac101_var_759 = fs1_var_739 * fac01_var_694(jl_var_786, lay_var_726)
        fac111_var_762 = fs1_var_739 * fac11_var_696(jl_var_786, lay_var_726)
        fac201_var_760 = 0.0D0
        fac211_var_763 = 0.0D0
      END IF
      IF (specparm_var_738 .LT. 0.125D0) THEN
        tau_major_var_774(1 : ng9) = speccomb_var_730 * (fac000_var_752 * absa_var_208(ind0_var_719, 1 : 12) + fac100_var_753 * absa_var_208(ind0_var_719 + 1, 1 : 12) + fac200_var_754 * absa_var_208(ind0_var_719 + 2, 1 : 12) + fac010_var_755 * absa_var_208(ind0_var_719 + 9, 1 : 12) + fac110_var_756 * absa_var_208(ind0_var_719 + 10, 1 : 12) + fac210_var_757 * absa_var_208(ind0_var_719 + 11, 1 : 12))
      ELSE IF (specparm_var_738 .GT. 0.875D0) THEN
        tau_major_var_774(1 : ng9) = speccomb_var_730 * (fac200_var_754 * absa_var_208(ind0_var_719 - 1, 1 : 12) + fac100_var_753 * absa_var_208(ind0_var_719, 1 : 12) + fac000_var_752 * absa_var_208(ind0_var_719 + 1, 1 : 12) + fac210_var_757 * absa_var_208(ind0_var_719 + 8, 1 : 12) + fac110_var_756 * absa_var_208(ind0_var_719 + 9, 1 : 12) + fac010_var_755 * absa_var_208(ind0_var_719 + 10, 1 : 12))
      ELSE
        tau_major_var_774(1 : ng9) = speccomb_var_730 * (fac000_var_752 * absa_var_208(ind0_var_719, 1 : 12) + fac100_var_753 * absa_var_208(ind0_var_719 + 1, 1 : 12) + fac010_var_755 * absa_var_208(ind0_var_719 + 9, 1 : 12) + fac110_var_756 * absa_var_208(ind0_var_719 + 10, 1 : 12))
      END IF
      IF (specparm1_var_741 .LT. 0.125D0) THEN
        tau_major1_var_775(1 : ng9) = speccomb1_var_731 * (fac001_var_758 * absa_var_208(ind1_var_720, 1 : 12) + fac101_var_759 * absa_var_208(ind1_var_720 + 1, 1 : 12) + fac201_var_760 * absa_var_208(ind1_var_720 + 2, 1 : 12) + fac011_var_761 * absa_var_208(ind1_var_720 + 9, 1 : 12) + fac111_var_762 * absa_var_208(ind1_var_720 + 10, 1 : 12) + fac211_var_763 * absa_var_208(ind1_var_720 + 11, 1 : 12))
      ELSE IF (specparm1_var_741 .GT. 0.875D0) THEN
        tau_major1_var_775(1 : ng9) = speccomb1_var_731 * (fac201_var_760 * absa_var_208(ind1_var_720 - 1, 1 : 12) + fac101_var_759 * absa_var_208(ind1_var_720, 1 : 12) + fac001_var_758 * absa_var_208(ind1_var_720 + 1, 1 : 12) + fac211_var_763 * absa_var_208(ind1_var_720 + 8, 1 : 12) + fac111_var_762 * absa_var_208(ind1_var_720 + 9, 1 : 12) + fac011_var_761 * absa_var_208(ind1_var_720 + 10, 1 : 12))
      ELSE
        tau_major1_var_775(1 : ng9) = speccomb1_var_731 * (fac001_var_758 * absa_var_208(ind1_var_720, 1 : 12) + fac101_var_759 * absa_var_208(ind1_var_720 + 1, 1 : 12) + fac011_var_761 * absa_var_208(ind1_var_720 + 9, 1 : 12) + fac111_var_762 * absa_var_208(ind1_var_720 + 10, 1 : 12))
      END IF
      DO ig_var_724 = 1, 12
        tauself_var_770 = selffac_var_708(jl_var_786, lay_var_726) * (selfref_var_212(inds_var_721, ig_var_724) + selffrac_var_709(jl_var_786, lay_var_726) * (selfref_var_212(inds_var_721 + 1, ig_var_724) - selfref_var_212(inds_var_721, ig_var_724)))
        taufor_var_769 = forfac_var_715(jl_var_786, lay_var_726) * (forref_var_213(indf_var_722, ig_var_724) + forfrac_var_716(jl_var_786, lay_var_726) * (forref_var_213(indf_var_722 + 1, ig_var_724) - forref_var_213(indf_var_722, ig_var_724)))
        n2om1_var_771 = ka_mn2o_var_210(jmn2o_var_728, indm_var_723, ig_var_724) + fmn2o_var_742 * (ka_mn2o_var_210(jmn2o_var_728 + 1, indm_var_723, ig_var_724) - ka_mn2o_var_210(jmn2o_var_728, indm_var_723, ig_var_724))
        n2om2_var_772 = ka_mn2o_var_210(jmn2o_var_728, indm_var_723 + 1, ig_var_724) + fmn2o_var_742 * (ka_mn2o_var_210(jmn2o_var_728 + 1, indm_var_723 + 1, ig_var_724) - ka_mn2o_var_210(jmn2o_var_728, indm_var_723 + 1, ig_var_724))
        absn2o_var_773 = n2om1_var_771 + minorfrac_var_717(jl_var_786, lay_var_726) * (n2om2_var_772 - n2om1_var_771)
        taug_var_691(jl_var_786, 96 + ig_var_724, lay_var_726) = tau_major_var_774(ig_var_724) + tau_major1_var_775(ig_var_724) + tauself_var_770 + taufor_var_769 + adjcoln2o_var_749 * absn2o_var_773
        fracs_var_711(jl_var_786, 96 + ig_var_724, lay_var_726) = fracrefa_var_206(ig_var_724, jpl_var_729) + fpl_var_745 * (fracrefa_var_206(ig_var_724, jpl_var_729 + 1) - fracrefa_var_206(ig_var_724, jpl_var_729))
      END DO
    END DO
  END DO
  DO lay_var_726 = laytrop_max_var_777 + 1, klev_var_690
    DO jl_var_786 = kidia_var_688, kfdia_var_689
      chi_n2o_var_751 = coln2o_var_702(jl_var_786, lay_var_726) / (coldry_var_704(jl_var_786, lay_var_726))
      ratn2o_var_750 = 1D+20 * chi_n2o_var_751 / chi_mls(4, jp_var_697(jl_var_786, lay_var_726) + 1)
      IF (ratn2o_var_750 .GT. 1.5D0) THEN
        adjfac_var_748 = 0.5D0 + (ratn2o_var_750 - 0.5D0) ** 0.65D0
        adjcoln2o_var_749 = adjfac_var_748 * chi_mls(4, jp_var_697(jl_var_786, lay_var_726) + 1) * coldry_var_704(jl_var_786, lay_var_726) * 1D-20
      ELSE
        adjcoln2o_var_749 = coln2o_var_702(jl_var_786, lay_var_726)
      END IF
      ind0_var_719 = ((jp_var_697(jl_var_786, lay_var_726) - 13) * 5 + (jt_var_698(jl_var_786, lay_var_726) - 1)) * nspb_var_217(9) + 1
      ind1_var_720 = ((jp_var_697(jl_var_786, lay_var_726) - 12) * 5 + (jt1_var_699(jl_var_786, lay_var_726) - 1)) * nspb_var_217(9) + 1
      indm_var_723 = indminor_var_718(jl_var_786, lay_var_726)
      DO ig_var_724 = 1, 12
        absn2o_var_773 = kb_mn2o_var_211(indm_var_723, ig_var_724) + minorfrac_var_717(jl_var_786, lay_var_726) * (kb_mn2o_var_211(indm_var_723 + 1, ig_var_724) - kb_mn2o_var_211(indm_var_723, ig_var_724))
        taug_var_691(jl_var_786, 96 + ig_var_724, lay_var_726) = colch4_var_703(jl_var_786, lay_var_726) * (fac00_var_693(jl_var_786, lay_var_726) * absb_var_209(ind0_var_719, ig_var_724) + fac10_var_695(jl_var_786, lay_var_726) * absb_var_209(ind0_var_719 + 1, ig_var_724) + fac01_var_694(jl_var_786, lay_var_726) * absb_var_209(ind1_var_720, ig_var_724) + fac11_var_696(jl_var_786, lay_var_726) * absb_var_209(ind1_var_720 + 1, ig_var_724)) + adjcoln2o_var_749 * absn2o_var_773
        fracs_var_711(jl_var_786, 96 + ig_var_724, lay_var_726) = fracrefb_var_207(ig_var_724)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_777 /= laytrop_min_var_776) THEN
    DO lay_var_726 = laytrop_min_var_776 + 1, laytrop_max_var_777
      ixc0_var_783 = ixc_var_778(lay_var_726)
      DO ixp_var_784 = 1, ixc0_var_783
        jl_var_786 = ixlow_var_779(ixp_var_784, lay_var_726)
        speccomb_var_730 = colh2o_var_701(jl_var_786, lay_var_726) + rat_h2och4_var_712(jl_var_786, lay_var_726) * colch4_var_703(jl_var_786, lay_var_726)
        specparm_var_738 = MIN(colh2o_var_701(jl_var_786, lay_var_726) / speccomb_var_730, oneminus_var_700)
        specmult_var_737 = 8.0D0 * (specparm_var_738)
        js_var_725 = 1 + INT(specmult_var_737)
        fs_var_736 = ((specmult_var_737) - AINT((specmult_var_737)))
        speccomb1_var_731 = colh2o_var_701(jl_var_786, lay_var_726) + rat_h2och4_1_var_713(jl_var_786, lay_var_726) * colch4_var_703(jl_var_786, lay_var_726)
        specparm1_var_741 = MIN(colh2o_var_701(jl_var_786, lay_var_726) / speccomb1_var_731, oneminus_var_700)
        specmult1_var_740 = 8.0D0 * (specparm1_var_741)
        js1_var_727 = 1 + INT(specmult1_var_740)
        fs1_var_739 = ((specmult1_var_740) - AINT((specmult1_var_740)))
        speccomb_mn2o_var_732 = colh2o_var_701(jl_var_786, lay_var_726) + refrat_m_a_var_735 * colch4_var_703(jl_var_786, lay_var_726)
        specparm_mn2o_var_744 = MIN(colh2o_var_701(jl_var_786, lay_var_726) / speccomb_mn2o_var_732, oneminus_var_700)
        specmult_mn2o_var_743 = 8.0D0 * specparm_mn2o_var_744
        jmn2o_var_728 = 1 + INT(specmult_mn2o_var_743)
        fmn2o_var_742 = ((specmult_mn2o_var_743) - AINT((specmult_mn2o_var_743)))
        chi_n2o_var_751 = coln2o_var_702(jl_var_786, lay_var_726) / (coldry_var_704(jl_var_786, lay_var_726))
        ratn2o_var_750 = 1D+20 * chi_n2o_var_751 / chi_mls(4, jp_var_697(jl_var_786, lay_var_726) + 1)
        IF (ratn2o_var_750 .GT. 1.5D0) THEN
          adjfac_var_748 = 0.5D0 + (ratn2o_var_750 - 0.5D0) ** 0.65D0
          adjcoln2o_var_749 = adjfac_var_748 * chi_mls(4, jp_var_697(jl_var_786, lay_var_726) + 1) * coldry_var_704(jl_var_786, lay_var_726) * 1D-20
        ELSE
          adjcoln2o_var_749 = coln2o_var_702(jl_var_786, lay_var_726)
        END IF
        speccomb_planck_var_733 = colh2o_var_701(jl_var_786, lay_var_726) + refrat_planck_a_var_734 * colch4_var_703(jl_var_786, lay_var_726)
        specparm_planck_var_747 = MIN(colh2o_var_701(jl_var_786, lay_var_726) / speccomb_planck_var_733, oneminus_var_700)
        specmult_planck_var_746 = 8.0D0 * specparm_planck_var_747
        jpl_var_729 = 1 + INT(specmult_planck_var_746)
        fpl_var_745 = ((specmult_planck_var_746) - AINT((specmult_planck_var_746)))
        ind0_var_719 = ((jp_var_697(jl_var_786, lay_var_726) - 1) * 5 + (jt_var_698(jl_var_786, lay_var_726) - 1)) * nspa_var_216(9) + js_var_725
        ind1_var_720 = (jp_var_697(jl_var_786, lay_var_726) * 5 + (jt1_var_699(jl_var_786, lay_var_726) - 1)) * nspa_var_216(9) + js1_var_727
        inds_var_721 = indself_var_710(jl_var_786, lay_var_726)
        indf_var_722 = indfor_var_714(jl_var_786, lay_var_726)
        indm_var_723 = indminor_var_718(jl_var_786, lay_var_726)
        IF (specparm_var_738 .LT. 0.125D0) THEN
          p_var_764 = fs_var_736 - 1.0D0
          p4_var_765 = p_var_764 ** 4
          fk0_var_766 = p4_var_765
          fk1_var_767 = 1.0D0 - p_var_764 - 2.0D0 * p4_var_765
          fk2_var_768 = p_var_764 + p4_var_765
          fac000_var_752 = fk0_var_766 * fac00_var_693(jl_var_786, lay_var_726)
          fac100_var_753 = fk1_var_767 * fac00_var_693(jl_var_786, lay_var_726)
          fac200_var_754 = fk2_var_768 * fac00_var_693(jl_var_786, lay_var_726)
          fac010_var_755 = fk0_var_766 * fac10_var_695(jl_var_786, lay_var_726)
          fac110_var_756 = fk1_var_767 * fac10_var_695(jl_var_786, lay_var_726)
          fac210_var_757 = fk2_var_768 * fac10_var_695(jl_var_786, lay_var_726)
        ELSE IF (specparm_var_738 .GT. 0.875D0) THEN
          p_var_764 = - fs_var_736
          p4_var_765 = p_var_764 ** 4
          fk0_var_766 = p4_var_765
          fk1_var_767 = 1.0D0 - p_var_764 - 2.0D0 * p4_var_765
          fk2_var_768 = p_var_764 + p4_var_765
          fac000_var_752 = fk0_var_766 * fac00_var_693(jl_var_786, lay_var_726)
          fac100_var_753 = fk1_var_767 * fac00_var_693(jl_var_786, lay_var_726)
          fac200_var_754 = fk2_var_768 * fac00_var_693(jl_var_786, lay_var_726)
          fac010_var_755 = fk0_var_766 * fac10_var_695(jl_var_786, lay_var_726)
          fac110_var_756 = fk1_var_767 * fac10_var_695(jl_var_786, lay_var_726)
          fac210_var_757 = fk2_var_768 * fac10_var_695(jl_var_786, lay_var_726)
        ELSE
          fac000_var_752 = (1.0D0 - fs_var_736) * fac00_var_693(jl_var_786, lay_var_726)
          fac010_var_755 = (1.0D0 - fs_var_736) * fac10_var_695(jl_var_786, lay_var_726)
          fac100_var_753 = fs_var_736 * fac00_var_693(jl_var_786, lay_var_726)
          fac110_var_756 = fs_var_736 * fac10_var_695(jl_var_786, lay_var_726)
          fac200_var_754 = 0.0D0
          fac210_var_757 = 0.0D0
        END IF
        IF (specparm1_var_741 .LT. 0.125D0) THEN
          p_var_764 = fs1_var_739 - 1.0D0
          p4_var_765 = p_var_764 ** 4
          fk0_var_766 = p4_var_765
          fk1_var_767 = 1.0D0 - p_var_764 - 2.0D0 * p4_var_765
          fk2_var_768 = p_var_764 + p4_var_765
          fac001_var_758 = fk0_var_766 * fac01_var_694(jl_var_786, lay_var_726)
          fac101_var_759 = fk1_var_767 * fac01_var_694(jl_var_786, lay_var_726)
          fac201_var_760 = fk2_var_768 * fac01_var_694(jl_var_786, lay_var_726)
          fac011_var_761 = fk0_var_766 * fac11_var_696(jl_var_786, lay_var_726)
          fac111_var_762 = fk1_var_767 * fac11_var_696(jl_var_786, lay_var_726)
          fac211_var_763 = fk2_var_768 * fac11_var_696(jl_var_786, lay_var_726)
        ELSE IF (specparm1_var_741 .GT. 0.875D0) THEN
          p_var_764 = - fs1_var_739
          p4_var_765 = p_var_764 ** 4
          fk0_var_766 = p4_var_765
          fk1_var_767 = 1.0D0 - p_var_764 - 2.0D0 * p4_var_765
          fk2_var_768 = p_var_764 + p4_var_765
          fac001_var_758 = fk0_var_766 * fac01_var_694(jl_var_786, lay_var_726)
          fac101_var_759 = fk1_var_767 * fac01_var_694(jl_var_786, lay_var_726)
          fac201_var_760 = fk2_var_768 * fac01_var_694(jl_var_786, lay_var_726)
          fac011_var_761 = fk0_var_766 * fac11_var_696(jl_var_786, lay_var_726)
          fac111_var_762 = fk1_var_767 * fac11_var_696(jl_var_786, lay_var_726)
          fac211_var_763 = fk2_var_768 * fac11_var_696(jl_var_786, lay_var_726)
        ELSE
          fac001_var_758 = (1.0D0 - fs1_var_739) * fac01_var_694(jl_var_786, lay_var_726)
          fac011_var_761 = (1.0D0 - fs1_var_739) * fac11_var_696(jl_var_786, lay_var_726)
          fac101_var_759 = fs1_var_739 * fac01_var_694(jl_var_786, lay_var_726)
          fac111_var_762 = fs1_var_739 * fac11_var_696(jl_var_786, lay_var_726)
          fac201_var_760 = 0.0D0
          fac211_var_763 = 0.0D0
        END IF
        IF (specparm_var_738 .LT. 0.125D0) THEN
          tau_major_var_774(1 : ng9) = speccomb_var_730 * (fac000_var_752 * absa_var_208(ind0_var_719, 1 : 12) + fac100_var_753 * absa_var_208(ind0_var_719 + 1, 1 : 12) + fac200_var_754 * absa_var_208(ind0_var_719 + 2, 1 : 12) + fac010_var_755 * absa_var_208(ind0_var_719 + 9, 1 : 12) + fac110_var_756 * absa_var_208(ind0_var_719 + 10, 1 : 12) + fac210_var_757 * absa_var_208(ind0_var_719 + 11, 1 : 12))
        ELSE IF (specparm_var_738 .GT. 0.875D0) THEN
          tau_major_var_774(1 : ng9) = speccomb_var_730 * (fac200_var_754 * absa_var_208(ind0_var_719 - 1, 1 : 12) + fac100_var_753 * absa_var_208(ind0_var_719, 1 : 12) + fac000_var_752 * absa_var_208(ind0_var_719 + 1, 1 : 12) + fac210_var_757 * absa_var_208(ind0_var_719 + 8, 1 : 12) + fac110_var_756 * absa_var_208(ind0_var_719 + 9, 1 : 12) + fac010_var_755 * absa_var_208(ind0_var_719 + 10, 1 : 12))
        ELSE
          tau_major_var_774(1 : ng9) = speccomb_var_730 * (fac000_var_752 * absa_var_208(ind0_var_719, 1 : 12) + fac100_var_753 * absa_var_208(ind0_var_719 + 1, 1 : 12) + fac010_var_755 * absa_var_208(ind0_var_719 + 9, 1 : 12) + fac110_var_756 * absa_var_208(ind0_var_719 + 10, 1 : 12))
        END IF
        IF (specparm1_var_741 .LT. 0.125D0) THEN
          tau_major1_var_775(1 : ng9) = speccomb1_var_731 * (fac001_var_758 * absa_var_208(ind1_var_720, 1 : 12) + fac101_var_759 * absa_var_208(ind1_var_720 + 1, 1 : 12) + fac201_var_760 * absa_var_208(ind1_var_720 + 2, 1 : 12) + fac011_var_761 * absa_var_208(ind1_var_720 + 9, 1 : 12) + fac111_var_762 * absa_var_208(ind1_var_720 + 10, 1 : 12) + fac211_var_763 * absa_var_208(ind1_var_720 + 11, 1 : 12))
        ELSE IF (specparm1_var_741 .GT. 0.875D0) THEN
          tau_major1_var_775(1 : ng9) = speccomb1_var_731 * (fac201_var_760 * absa_var_208(ind1_var_720 - 1, 1 : 12) + fac101_var_759 * absa_var_208(ind1_var_720, 1 : 12) + fac001_var_758 * absa_var_208(ind1_var_720 + 1, 1 : 12) + fac211_var_763 * absa_var_208(ind1_var_720 + 8, 1 : 12) + fac111_var_762 * absa_var_208(ind1_var_720 + 9, 1 : 12) + fac011_var_761 * absa_var_208(ind1_var_720 + 10, 1 : 12))
        ELSE
          tau_major1_var_775(1 : ng9) = speccomb1_var_731 * (fac001_var_758 * absa_var_208(ind1_var_720, 1 : 12) + fac101_var_759 * absa_var_208(ind1_var_720 + 1, 1 : 12) + fac011_var_761 * absa_var_208(ind1_var_720 + 9, 1 : 12) + fac111_var_762 * absa_var_208(ind1_var_720 + 10, 1 : 12))
        END IF
        DO ig_var_724 = 1, 12
          tauself_var_770 = selffac_var_708(jl_var_786, lay_var_726) * (selfref_var_212(inds_var_721, ig_var_724) + selffrac_var_709(jl_var_786, lay_var_726) * (selfref_var_212(inds_var_721 + 1, ig_var_724) - selfref_var_212(inds_var_721, ig_var_724)))
          taufor_var_769 = forfac_var_715(jl_var_786, lay_var_726) * (forref_var_213(indf_var_722, ig_var_724) + forfrac_var_716(jl_var_786, lay_var_726) * (forref_var_213(indf_var_722 + 1, ig_var_724) - forref_var_213(indf_var_722, ig_var_724)))
          n2om1_var_771 = ka_mn2o_var_210(jmn2o_var_728, indm_var_723, ig_var_724) + fmn2o_var_742 * (ka_mn2o_var_210(jmn2o_var_728 + 1, indm_var_723, ig_var_724) - ka_mn2o_var_210(jmn2o_var_728, indm_var_723, ig_var_724))
          n2om2_var_772 = ka_mn2o_var_210(jmn2o_var_728, indm_var_723 + 1, ig_var_724) + fmn2o_var_742 * (ka_mn2o_var_210(jmn2o_var_728 + 1, indm_var_723 + 1, ig_var_724) - ka_mn2o_var_210(jmn2o_var_728, indm_var_723 + 1, ig_var_724))
          absn2o_var_773 = n2om1_var_771 + minorfrac_var_717(jl_var_786, lay_var_726) * (n2om2_var_772 - n2om1_var_771)
          taug_var_691(jl_var_786, 96 + ig_var_724, lay_var_726) = tau_major_var_774(ig_var_724) + tau_major1_var_775(ig_var_724) + tauself_var_770 + taufor_var_769 + adjcoln2o_var_749 * absn2o_var_773
          fracs_var_711(jl_var_786, 96 + ig_var_724, lay_var_726) = fracrefa_var_206(ig_var_724, jpl_var_729) + fpl_var_745 * (fracrefa_var_206(ig_var_724, jpl_var_729 + 1) - fracrefa_var_206(ig_var_724, jpl_var_729))
        END DO
      END DO
      ixc0_var_783 = kfdia_var_689 - kidia_var_688 + 1 - ixc0_var_783
      DO ixp_var_784 = 1, ixc0_var_783
        jl_var_786 = ixhigh_var_780(ixp_var_784, lay_var_726)
        chi_n2o_var_751 = coln2o_var_702(jl_var_786, lay_var_726) / (coldry_var_704(jl_var_786, lay_var_726))
        ratn2o_var_750 = 1D+20 * chi_n2o_var_751 / chi_mls(4, jp_var_697(jl_var_786, lay_var_726) + 1)
        IF (ratn2o_var_750 .GT. 1.5D0) THEN
          adjfac_var_748 = 0.5D0 + (ratn2o_var_750 - 0.5D0) ** 0.65D0
          adjcoln2o_var_749 = adjfac_var_748 * chi_mls(4, jp_var_697(jl_var_786, lay_var_726) + 1) * coldry_var_704(jl_var_786, lay_var_726) * 1D-20
        ELSE
          adjcoln2o_var_749 = coln2o_var_702(jl_var_786, lay_var_726)
        END IF
        ind0_var_719 = ((jp_var_697(jl_var_786, lay_var_726) - 13) * 5 + (jt_var_698(jl_var_786, lay_var_726) - 1)) * nspb_var_217(9) + 1
        ind1_var_720 = ((jp_var_697(jl_var_786, lay_var_726) - 12) * 5 + (jt1_var_699(jl_var_786, lay_var_726) - 1)) * nspb_var_217(9) + 1
        indm_var_723 = indminor_var_718(jl_var_786, lay_var_726)
        DO ig_var_724 = 1, 12
          absn2o_var_773 = kb_mn2o_var_211(indm_var_723, ig_var_724) + minorfrac_var_717(jl_var_786, lay_var_726) * (kb_mn2o_var_211(indm_var_723 + 1, ig_var_724) - kb_mn2o_var_211(indm_var_723, ig_var_724))
          taug_var_691(jl_var_786, 96 + ig_var_724, lay_var_726) = colch4_var_703(jl_var_786, lay_var_726) * (fac00_var_693(jl_var_786, lay_var_726) * absb_var_209(ind0_var_719, ig_var_724) + fac10_var_695(jl_var_786, lay_var_726) * absb_var_209(ind0_var_719 + 1, ig_var_724) + fac01_var_694(jl_var_786, lay_var_726) * absb_var_209(ind1_var_720, ig_var_724) + fac11_var_696(jl_var_786, lay_var_726) * absb_var_209(ind1_var_720 + 1, ig_var_724)) + adjcoln2o_var_749 * absn2o_var_773
          fracs_var_711(jl_var_786, 96 + ig_var_724, lay_var_726) = fracrefb_var_207(ig_var_724)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol9
SUBROUTINE rrtm_taumol11(kidia_var_787, kfdia_var_788, klev_var_789, taug_var_790, p_tauaerl_var_791, fac00_var_792, fac01_var_793, fac10_var_794, fac11_var_795, forfac_var_806, forfrac_var_807, indfor_var_805, jp_var_796, jt_var_797, jt1_var_798, colh2o_var_799, colo2, laytrop_var_800, selffac_var_801, selffrac_var_802, indself_var_803, fracs_var_804, minorfrac_var_808, indminor_var_809, scaleminor_var_810)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta11, ONLY: absa_var_122, absb_var_123, forref_var_125, fracrefa_var_120, fracrefb_var_121, ka_mo2, kb_mo2, selfref_var_124
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_787
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_788
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_789
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_790(kidia_var_787 : kfdia_var_788, 140, klev_var_789)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_791(kidia_var_787 : kfdia_var_788, klev_var_789, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_792(kidia_var_787 : kfdia_var_788, klev_var_789)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_793(kidia_var_787 : kfdia_var_788, klev_var_789)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_794(kidia_var_787 : kfdia_var_788, klev_var_789)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_795(kidia_var_787 : kfdia_var_788, klev_var_789)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_796(kidia_var_787 : kfdia_var_788, klev_var_789)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_797(kidia_var_787 : kfdia_var_788, klev_var_789)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_798(kidia_var_787 : kfdia_var_788, klev_var_789)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_799(kidia_var_787 : kfdia_var_788, klev_var_789)
  REAL(KIND = 8), INTENT(IN) :: colo2(kidia_var_787 : kfdia_var_788, klev_var_789)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_800(kidia_var_787 : kfdia_var_788)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_801(kidia_var_787 : kfdia_var_788, klev_var_789)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_802(kidia_var_787 : kfdia_var_788, klev_var_789)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_803(kidia_var_787 : kfdia_var_788, klev_var_789)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_804(kidia_var_787 : kfdia_var_788, 140, klev_var_789)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_805(kidia_var_787 : kfdia_var_788, klev_var_789)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_806(kidia_var_787 : kfdia_var_788, klev_var_789)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_807(kidia_var_787 : kfdia_var_788, klev_var_789)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_808(kidia_var_787 : kfdia_var_788, klev_var_789)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_809(kidia_var_787 : kfdia_var_788, klev_var_789)
  REAL(KIND = 8), INTENT(IN) :: scaleminor_var_810(kidia_var_787 : kfdia_var_788, klev_var_789)
  INTEGER(KIND = 4) :: ind0_var_811, ind1_var_812
  INTEGER(KIND = 4) :: inds_var_813, indf_var_814, indm_var_815
  INTEGER(KIND = 4) :: ig_var_816, lay_var_817
  REAL(KIND = 8) :: taufor_var_818, tauself_var_819, scaleo2, tauo2
  INTEGER(KIND = 4) :: laytrop_min_var_820, laytrop_max_var_821
  INTEGER(KIND = 4) :: ixc_var_822(klev_var_789), ixlow_var_823(kfdia_var_788, klev_var_789), ixhigh_var_824(kfdia_var_788, klev_var_789)
  INTEGER(KIND = 4) :: ich_var_825, icl_var_826, ixc0_var_827, ixp_var_828, jc_var_829, jl_var_830
  laytrop_min_var_820 = MINVAL(laytrop_var_800)
  laytrop_max_var_821 = MAXVAL(laytrop_var_800)
  ixlow_var_823 = 0
  ixhigh_var_824 = 0
  ixc_var_822 = 0
  DO lay_var_817 = laytrop_min_var_820 + 1, laytrop_max_var_821
    icl_var_826 = 0
    ich_var_825 = 0
    DO jc_var_829 = kidia_var_787, kfdia_var_788
      IF (lay_var_817 <= laytrop_var_800(jc_var_829)) THEN
        icl_var_826 = icl_var_826 + 1
        ixlow_var_823(icl_var_826, lay_var_817) = jc_var_829
      ELSE
        ich_var_825 = ich_var_825 + 1
        ixhigh_var_824(ich_var_825, lay_var_817) = jc_var_829
      END IF
    END DO
    ixc_var_822(lay_var_817) = icl_var_826
  END DO
  DO lay_var_817 = 1, laytrop_min_var_820
    DO jl_var_830 = kidia_var_787, kfdia_var_788
      ind0_var_811 = ((jp_var_796(jl_var_830, lay_var_817) - 1) * 5 + (jt_var_797(jl_var_830, lay_var_817) - 1)) * nspa_var_216(11) + 1
      ind1_var_812 = (jp_var_796(jl_var_830, lay_var_817) * 5 + (jt1_var_798(jl_var_830, lay_var_817) - 1)) * nspa_var_216(11) + 1
      inds_var_813 = indself_var_803(jl_var_830, lay_var_817)
      indf_var_814 = indfor_var_805(jl_var_830, lay_var_817)
      indm_var_815 = indminor_var_809(jl_var_830, lay_var_817)
      scaleo2 = colo2(jl_var_830, lay_var_817) * scaleminor_var_810(jl_var_830, lay_var_817)
      DO ig_var_816 = 1, 8
        tauself_var_819 = selffac_var_801(jl_var_830, lay_var_817) * (selfref_var_124(inds_var_813, ig_var_816) + selffrac_var_802(jl_var_830, lay_var_817) * (selfref_var_124(inds_var_813 + 1, ig_var_816) - selfref_var_124(inds_var_813, ig_var_816)))
        taufor_var_818 = forfac_var_806(jl_var_830, lay_var_817) * (forref_var_125(indf_var_814, ig_var_816) + forfrac_var_807(jl_var_830, lay_var_817) * (forref_var_125(indf_var_814 + 1, ig_var_816) - forref_var_125(indf_var_814, ig_var_816)))
        tauo2 = scaleo2 * (ka_mo2(indm_var_815, ig_var_816) + minorfrac_var_808(jl_var_830, lay_var_817) * (ka_mo2(indm_var_815 + 1, ig_var_816) - ka_mo2(indm_var_815, ig_var_816)))
        taug_var_790(jl_var_830, 114 + ig_var_816, lay_var_817) = colh2o_var_799(jl_var_830, lay_var_817) * (fac00_var_792(jl_var_830, lay_var_817) * absa_var_122(ind0_var_811, ig_var_816) + fac10_var_794(jl_var_830, lay_var_817) * absa_var_122(ind0_var_811 + 1, ig_var_816) + fac01_var_793(jl_var_830, lay_var_817) * absa_var_122(ind1_var_812, ig_var_816) + fac11_var_795(jl_var_830, lay_var_817) * absa_var_122(ind1_var_812 + 1, ig_var_816)) + tauself_var_819 + taufor_var_818 + tauo2
        fracs_var_804(jl_var_830, 114 + ig_var_816, lay_var_817) = fracrefa_var_120(ig_var_816)
      END DO
    END DO
  END DO
  DO lay_var_817 = laytrop_max_var_821 + 1, klev_var_789
    DO jl_var_830 = kidia_var_787, kfdia_var_788
      ind0_var_811 = ((jp_var_796(jl_var_830, lay_var_817) - 13) * 5 + (jt_var_797(jl_var_830, lay_var_817) - 1)) * nspb_var_217(11) + 1
      ind1_var_812 = ((jp_var_796(jl_var_830, lay_var_817) - 12) * 5 + (jt1_var_798(jl_var_830, lay_var_817) - 1)) * nspb_var_217(11) + 1
      indf_var_814 = indfor_var_805(jl_var_830, lay_var_817)
      indm_var_815 = indminor_var_809(jl_var_830, lay_var_817)
      scaleo2 = colo2(jl_var_830, lay_var_817) * scaleminor_var_810(jl_var_830, lay_var_817)
      DO ig_var_816 = 1, 8
        taufor_var_818 = forfac_var_806(jl_var_830, lay_var_817) * (forref_var_125(indf_var_814, ig_var_816) + forfrac_var_807(jl_var_830, lay_var_817) * (forref_var_125(indf_var_814 + 1, ig_var_816) - forref_var_125(indf_var_814, ig_var_816)))
        tauo2 = scaleo2 * (kb_mo2(indm_var_815, ig_var_816) + minorfrac_var_808(jl_var_830, lay_var_817) * (kb_mo2(indm_var_815 + 1, ig_var_816) - kb_mo2(indm_var_815, ig_var_816)))
        taug_var_790(jl_var_830, 114 + ig_var_816, lay_var_817) = colh2o_var_799(jl_var_830, lay_var_817) * (fac00_var_792(jl_var_830, lay_var_817) * absb_var_123(ind0_var_811, ig_var_816) + fac10_var_794(jl_var_830, lay_var_817) * absb_var_123(ind0_var_811 + 1, ig_var_816) + fac01_var_793(jl_var_830, lay_var_817) * absb_var_123(ind1_var_812, ig_var_816) + fac11_var_795(jl_var_830, lay_var_817) * absb_var_123(ind1_var_812 + 1, ig_var_816)) + taufor_var_818 + tauo2
        fracs_var_804(jl_var_830, 114 + ig_var_816, lay_var_817) = fracrefb_var_121(ig_var_816)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_821 /= laytrop_min_var_820) THEN
    DO lay_var_817 = laytrop_min_var_820 + 1, laytrop_max_var_821
      ixc0_var_827 = ixc_var_822(lay_var_817)
      DO ixp_var_828 = 1, ixc0_var_827
        jl_var_830 = ixlow_var_823(ixp_var_828, lay_var_817)
        ind0_var_811 = ((jp_var_796(jl_var_830, lay_var_817) - 1) * 5 + (jt_var_797(jl_var_830, lay_var_817) - 1)) * nspa_var_216(11) + 1
        ind1_var_812 = (jp_var_796(jl_var_830, lay_var_817) * 5 + (jt1_var_798(jl_var_830, lay_var_817) - 1)) * nspa_var_216(11) + 1
        inds_var_813 = indself_var_803(jl_var_830, lay_var_817)
        indf_var_814 = indfor_var_805(jl_var_830, lay_var_817)
        indm_var_815 = indminor_var_809(jl_var_830, lay_var_817)
        scaleo2 = colo2(jl_var_830, lay_var_817) * scaleminor_var_810(jl_var_830, lay_var_817)
        DO ig_var_816 = 1, 8
          tauself_var_819 = selffac_var_801(jl_var_830, lay_var_817) * (selfref_var_124(inds_var_813, ig_var_816) + selffrac_var_802(jl_var_830, lay_var_817) * (selfref_var_124(inds_var_813 + 1, ig_var_816) - selfref_var_124(inds_var_813, ig_var_816)))
          taufor_var_818 = forfac_var_806(jl_var_830, lay_var_817) * (forref_var_125(indf_var_814, ig_var_816) + forfrac_var_807(jl_var_830, lay_var_817) * (forref_var_125(indf_var_814 + 1, ig_var_816) - forref_var_125(indf_var_814, ig_var_816)))
          tauo2 = scaleo2 * (ka_mo2(indm_var_815, ig_var_816) + minorfrac_var_808(jl_var_830, lay_var_817) * (ka_mo2(indm_var_815 + 1, ig_var_816) - ka_mo2(indm_var_815, ig_var_816)))
          taug_var_790(jl_var_830, 114 + ig_var_816, lay_var_817) = colh2o_var_799(jl_var_830, lay_var_817) * (fac00_var_792(jl_var_830, lay_var_817) * absa_var_122(ind0_var_811, ig_var_816) + fac10_var_794(jl_var_830, lay_var_817) * absa_var_122(ind0_var_811 + 1, ig_var_816) + fac01_var_793(jl_var_830, lay_var_817) * absa_var_122(ind1_var_812, ig_var_816) + fac11_var_795(jl_var_830, lay_var_817) * absa_var_122(ind1_var_812 + 1, ig_var_816)) + tauself_var_819 + taufor_var_818 + tauo2
          fracs_var_804(jl_var_830, 114 + ig_var_816, lay_var_817) = fracrefa_var_120(ig_var_816)
        END DO
      END DO
      ixc0_var_827 = kfdia_var_788 - kidia_var_787 + 1 - ixc0_var_827
      DO ixp_var_828 = 1, ixc0_var_827
        jl_var_830 = ixhigh_var_824(ixp_var_828, lay_var_817)
        ind0_var_811 = ((jp_var_796(jl_var_830, lay_var_817) - 13) * 5 + (jt_var_797(jl_var_830, lay_var_817) - 1)) * nspb_var_217(11) + 1
        ind1_var_812 = ((jp_var_796(jl_var_830, lay_var_817) - 12) * 5 + (jt1_var_798(jl_var_830, lay_var_817) - 1)) * nspb_var_217(11) + 1
        indf_var_814 = indfor_var_805(jl_var_830, lay_var_817)
        indm_var_815 = indminor_var_809(jl_var_830, lay_var_817)
        scaleo2 = colo2(jl_var_830, lay_var_817) * scaleminor_var_810(jl_var_830, lay_var_817)
        DO ig_var_816 = 1, 8
          taufor_var_818 = forfac_var_806(jl_var_830, lay_var_817) * (forref_var_125(indf_var_814, ig_var_816) + forfrac_var_807(jl_var_830, lay_var_817) * (forref_var_125(indf_var_814 + 1, ig_var_816) - forref_var_125(indf_var_814, ig_var_816)))
          tauo2 = scaleo2 * (kb_mo2(indm_var_815, ig_var_816) + minorfrac_var_808(jl_var_830, lay_var_817) * (kb_mo2(indm_var_815 + 1, ig_var_816) - kb_mo2(indm_var_815, ig_var_816)))
          taug_var_790(jl_var_830, 114 + ig_var_816, lay_var_817) = colh2o_var_799(jl_var_830, lay_var_817) * (fac00_var_792(jl_var_830, lay_var_817) * absb_var_123(ind0_var_811, ig_var_816) + fac10_var_794(jl_var_830, lay_var_817) * absb_var_123(ind0_var_811 + 1, ig_var_816) + fac01_var_793(jl_var_830, lay_var_817) * absb_var_123(ind1_var_812, ig_var_816) + fac11_var_795(jl_var_830, lay_var_817) * absb_var_123(ind1_var_812 + 1, ig_var_816)) + taufor_var_818 + tauo2
          fracs_var_804(jl_var_830, 114 + ig_var_816, lay_var_817) = fracrefb_var_121(ig_var_816)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol11
SUBROUTINE rrtm_taumol10(kidia_var_831, kfdia_var_832, klev_var_833, taug_var_834, p_tauaerl_var_835, fac00_var_836, fac01_var_837, fac10_var_838, fac11_var_839, forfac_var_851, forfrac_var_850, indfor_var_849, jp_var_840, jt_var_841, jt1_var_842, colh2o_var_843, laytrop_var_844, selffac_var_846, selffrac_var_847, indself_var_848, fracs_var_845)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta10, ONLY: absa_var_116, absb_var_117, forref_var_119, fracrefa_var_114, fracrefb_var_115, selfref_var_118
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_831
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_832
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_833
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_834(kidia_var_831 : kfdia_var_832, 140, klev_var_833)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_835(kidia_var_831 : kfdia_var_832, klev_var_833, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_836(kidia_var_831 : kfdia_var_832, klev_var_833)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_837(kidia_var_831 : kfdia_var_832, klev_var_833)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_838(kidia_var_831 : kfdia_var_832, klev_var_833)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_839(kidia_var_831 : kfdia_var_832, klev_var_833)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_840(kidia_var_831 : kfdia_var_832, klev_var_833)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_841(kidia_var_831 : kfdia_var_832, klev_var_833)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_842(kidia_var_831 : kfdia_var_832, klev_var_833)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_843(kidia_var_831 : kfdia_var_832, klev_var_833)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_844(kidia_var_831 : kfdia_var_832)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_845(kidia_var_831 : kfdia_var_832, 140, klev_var_833)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_846(kidia_var_831 : kfdia_var_832, klev_var_833)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_847(kidia_var_831 : kfdia_var_832, klev_var_833)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_848(kidia_var_831 : kfdia_var_832, klev_var_833)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_849(kidia_var_831 : kfdia_var_832, klev_var_833)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_850(kidia_var_831 : kfdia_var_832, klev_var_833)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_851(kidia_var_831 : kfdia_var_832, klev_var_833)
  INTEGER(KIND = 4) :: ind0_var_852, ind1_var_853
  INTEGER(KIND = 4) :: inds_var_854, indf_var_855
  INTEGER(KIND = 4) :: ig_var_856, lay_var_857
  REAL(KIND = 8) :: taufor_var_858, tauself_var_859
  INTEGER(KIND = 4) :: laytrop_min_var_860, laytrop_max_var_861
  INTEGER(KIND = 4) :: ixc_var_862(klev_var_833), ixlow_var_863(kfdia_var_832, klev_var_833), ixhigh_var_864(kfdia_var_832, klev_var_833)
  INTEGER(KIND = 4) :: ich_var_865, icl_var_866, ixc0_var_867, ixp_var_868, jc_var_869, jl_var_870
  laytrop_min_var_860 = MINVAL(laytrop_var_844)
  laytrop_max_var_861 = MAXVAL(laytrop_var_844)
  ixlow_var_863 = 0
  ixhigh_var_864 = 0
  ixc_var_862 = 0
  DO lay_var_857 = laytrop_min_var_860 + 1, laytrop_max_var_861
    icl_var_866 = 0
    ich_var_865 = 0
    DO jc_var_869 = kidia_var_831, kfdia_var_832
      IF (lay_var_857 <= laytrop_var_844(jc_var_869)) THEN
        icl_var_866 = icl_var_866 + 1
        ixlow_var_863(icl_var_866, lay_var_857) = jc_var_869
      ELSE
        ich_var_865 = ich_var_865 + 1
        ixhigh_var_864(ich_var_865, lay_var_857) = jc_var_869
      END IF
    END DO
    ixc_var_862(lay_var_857) = icl_var_866
  END DO
  DO lay_var_857 = 1, laytrop_min_var_860
    DO jl_var_870 = kidia_var_831, kfdia_var_832
      ind0_var_852 = ((jp_var_840(jl_var_870, lay_var_857) - 1) * 5 + (jt_var_841(jl_var_870, lay_var_857) - 1)) * nspa_var_216(10) + 1
      ind1_var_853 = (jp_var_840(jl_var_870, lay_var_857) * 5 + (jt1_var_842(jl_var_870, lay_var_857) - 1)) * nspa_var_216(10) + 1
      inds_var_854 = indself_var_848(jl_var_870, lay_var_857)
      indf_var_855 = indfor_var_849(jl_var_870, lay_var_857)
      DO ig_var_856 = 1, 6
        tauself_var_859 = selffac_var_846(jl_var_870, lay_var_857) * (selfref_var_118(inds_var_854, ig_var_856) + selffrac_var_847(jl_var_870, lay_var_857) * (selfref_var_118(inds_var_854 + 1, ig_var_856) - selfref_var_118(inds_var_854, ig_var_856)))
        taufor_var_858 = forfac_var_851(jl_var_870, lay_var_857) * (forref_var_119(indf_var_855, ig_var_856) + forfrac_var_850(jl_var_870, lay_var_857) * (forref_var_119(indf_var_855 + 1, ig_var_856) - forref_var_119(indf_var_855, ig_var_856)))
        taug_var_834(jl_var_870, 108 + ig_var_856, lay_var_857) = colh2o_var_843(jl_var_870, lay_var_857) * (fac00_var_836(jl_var_870, lay_var_857) * absa_var_116(ind0_var_852, ig_var_856) + fac10_var_838(jl_var_870, lay_var_857) * absa_var_116(ind0_var_852 + 1, ig_var_856) + fac01_var_837(jl_var_870, lay_var_857) * absa_var_116(ind1_var_853, ig_var_856) + fac11_var_839(jl_var_870, lay_var_857) * absa_var_116(ind1_var_853 + 1, ig_var_856)) + tauself_var_859 + taufor_var_858
        fracs_var_845(jl_var_870, 108 + ig_var_856, lay_var_857) = fracrefa_var_114(ig_var_856)
      END DO
    END DO
  END DO
  DO lay_var_857 = laytrop_max_var_861 + 1, klev_var_833
    DO jl_var_870 = kidia_var_831, kfdia_var_832
      ind0_var_852 = ((jp_var_840(jl_var_870, lay_var_857) - 13) * 5 + (jt_var_841(jl_var_870, lay_var_857) - 1)) * nspb_var_217(10) + 1
      ind1_var_853 = ((jp_var_840(jl_var_870, lay_var_857) - 12) * 5 + (jt1_var_842(jl_var_870, lay_var_857) - 1)) * nspb_var_217(10) + 1
      indf_var_855 = indfor_var_849(jl_var_870, lay_var_857)
      DO ig_var_856 = 1, 6
        taufor_var_858 = forfac_var_851(jl_var_870, lay_var_857) * (forref_var_119(indf_var_855, ig_var_856) + forfrac_var_850(jl_var_870, lay_var_857) * (forref_var_119(indf_var_855 + 1, ig_var_856) - forref_var_119(indf_var_855, ig_var_856)))
        taug_var_834(jl_var_870, 108 + ig_var_856, lay_var_857) = colh2o_var_843(jl_var_870, lay_var_857) * (fac00_var_836(jl_var_870, lay_var_857) * absb_var_117(ind0_var_852, ig_var_856) + fac10_var_838(jl_var_870, lay_var_857) * absb_var_117(ind0_var_852 + 1, ig_var_856) + fac01_var_837(jl_var_870, lay_var_857) * absb_var_117(ind1_var_853, ig_var_856) + fac11_var_839(jl_var_870, lay_var_857) * absb_var_117(ind1_var_853 + 1, ig_var_856)) + taufor_var_858
        fracs_var_845(jl_var_870, 108 + ig_var_856, lay_var_857) = fracrefb_var_115(ig_var_856)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_861 /= laytrop_min_var_860) THEN
    DO lay_var_857 = laytrop_min_var_860 + 1, laytrop_max_var_861
      ixc0_var_867 = ixc_var_862(lay_var_857)
      DO ixp_var_868 = 1, ixc0_var_867
        jl_var_870 = ixlow_var_863(ixp_var_868, lay_var_857)
        ind0_var_852 = ((jp_var_840(jl_var_870, lay_var_857) - 1) * 5 + (jt_var_841(jl_var_870, lay_var_857) - 1)) * nspa_var_216(10) + 1
        ind1_var_853 = (jp_var_840(jl_var_870, lay_var_857) * 5 + (jt1_var_842(jl_var_870, lay_var_857) - 1)) * nspa_var_216(10) + 1
        inds_var_854 = indself_var_848(jl_var_870, lay_var_857)
        indf_var_855 = indfor_var_849(jl_var_870, lay_var_857)
        DO ig_var_856 = 1, 6
          tauself_var_859 = selffac_var_846(jl_var_870, lay_var_857) * (selfref_var_118(inds_var_854, ig_var_856) + selffrac_var_847(jl_var_870, lay_var_857) * (selfref_var_118(inds_var_854 + 1, ig_var_856) - selfref_var_118(inds_var_854, ig_var_856)))
          taufor_var_858 = forfac_var_851(jl_var_870, lay_var_857) * (forref_var_119(indf_var_855, ig_var_856) + forfrac_var_850(jl_var_870, lay_var_857) * (forref_var_119(indf_var_855 + 1, ig_var_856) - forref_var_119(indf_var_855, ig_var_856)))
          taug_var_834(jl_var_870, 108 + ig_var_856, lay_var_857) = colh2o_var_843(jl_var_870, lay_var_857) * (fac00_var_836(jl_var_870, lay_var_857) * absa_var_116(ind0_var_852, ig_var_856) + fac10_var_838(jl_var_870, lay_var_857) * absa_var_116(ind0_var_852 + 1, ig_var_856) + fac01_var_837(jl_var_870, lay_var_857) * absa_var_116(ind1_var_853, ig_var_856) + fac11_var_839(jl_var_870, lay_var_857) * absa_var_116(ind1_var_853 + 1, ig_var_856)) + tauself_var_859 + taufor_var_858
          fracs_var_845(jl_var_870, 108 + ig_var_856, lay_var_857) = fracrefa_var_114(ig_var_856)
        END DO
      END DO
      ixc0_var_867 = kfdia_var_832 - kidia_var_831 + 1 - ixc0_var_867
      DO ixp_var_868 = 1, ixc0_var_867
        jl_var_870 = ixhigh_var_864(ixp_var_868, lay_var_857)
        ind0_var_852 = ((jp_var_840(jl_var_870, lay_var_857) - 13) * 5 + (jt_var_841(jl_var_870, lay_var_857) - 1)) * nspb_var_217(10) + 1
        ind1_var_853 = ((jp_var_840(jl_var_870, lay_var_857) - 12) * 5 + (jt1_var_842(jl_var_870, lay_var_857) - 1)) * nspb_var_217(10) + 1
        indf_var_855 = indfor_var_849(jl_var_870, lay_var_857)
        DO ig_var_856 = 1, 6
          taufor_var_858 = forfac_var_851(jl_var_870, lay_var_857) * (forref_var_119(indf_var_855, ig_var_856) + forfrac_var_850(jl_var_870, lay_var_857) * (forref_var_119(indf_var_855 + 1, ig_var_856) - forref_var_119(indf_var_855, ig_var_856)))
          taug_var_834(jl_var_870, 108 + ig_var_856, lay_var_857) = colh2o_var_843(jl_var_870, lay_var_857) * (fac00_var_836(jl_var_870, lay_var_857) * absb_var_117(ind0_var_852, ig_var_856) + fac10_var_838(jl_var_870, lay_var_857) * absb_var_117(ind0_var_852 + 1, ig_var_856) + fac01_var_837(jl_var_870, lay_var_857) * absb_var_117(ind1_var_853, ig_var_856) + fac11_var_839(jl_var_870, lay_var_857) * absb_var_117(ind1_var_853 + 1, ig_var_856)) + taufor_var_858
          fracs_var_845(jl_var_870, 108 + ig_var_856, lay_var_857) = fracrefb_var_115(ig_var_856)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol10
SUBROUTINE rrtm_taumol8(kidia_var_871, kfdia_var_872, klev_var_873, taug_var_874, wx_var_875, p_tauaerl_var_876, fac00_var_877, fac01_var_878, fac10_var_879, fac11_var_880, forfac_var_896, forfrac_var_895, indfor_var_894, jp_var_881, jt_var_882, jt1_var_883, colh2o_var_884, colo3_var_885, coln2o_var_886, colco2_var_887, coldry_var_888, laytrop_var_889, selffac_var_890, selffrac_var_891, indself_var_892, fracs_var_893, minorfrac_var_897, indminor_var_898)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta8, ONLY: absa_var_197, absb_var_198, cfc12_var_196, cfc22adj, forref_var_205, fracrefa_var_194, fracrefb_var_195, ka_mco2_var_199, ka_mn2o_var_200, ka_mo3_var_201, kb_mco2_var_202, kb_mn2o_var_203, selfref_var_204
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_871
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_872
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_873
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_874(kidia_var_871 : kfdia_var_872, 140, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: wx_var_875(kidia_var_871 : kfdia_var_872, 4, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_876(kidia_var_871 : kfdia_var_872, klev_var_873, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_877(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_878(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_879(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_880(kidia_var_871 : kfdia_var_872, klev_var_873)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_881(kidia_var_871 : kfdia_var_872, klev_var_873)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_882(kidia_var_871 : kfdia_var_872, klev_var_873)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_883(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_884(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_885(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_886(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_887(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_888(kidia_var_871 : kfdia_var_872, klev_var_873)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_889(kidia_var_871 : kfdia_var_872)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_890(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_891(kidia_var_871 : kfdia_var_872, klev_var_873)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_892(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_893(kidia_var_871 : kfdia_var_872, 140, klev_var_873)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_894(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_895(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_896(kidia_var_871 : kfdia_var_872, klev_var_873)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_897(kidia_var_871 : kfdia_var_872, klev_var_873)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_898(kidia_var_871 : kfdia_var_872, klev_var_873)
  INTEGER(KIND = 4) :: ind0_var_899, ind1_var_900, inds_var_901, indf_var_902, indm_var_903
  INTEGER(KIND = 4) :: ig_var_904, lay_var_905
  REAL(KIND = 8) :: chi_co2_var_906, ratco2_var_907, adjfac_var_908, adjcolco2_var_909
  REAL(KIND = 8) :: taufor_var_910, tauself_var_911, abso3_var_912, absco2_var_913, absn2o_var_914
  INTEGER(KIND = 4) :: laytrop_min_var_915, laytrop_max_var_916
  INTEGER(KIND = 4) :: ixc_var_917(klev_var_873), ixlow_var_918(kfdia_var_872, klev_var_873), ixhigh_var_919(kfdia_var_872, klev_var_873)
  INTEGER(KIND = 4) :: ich_var_920, icl_var_921, ixc0_var_922, ixp_var_923, jc_var_924, jl_var_925
  laytrop_min_var_915 = MINVAL(laytrop_var_889)
  laytrop_max_var_916 = MAXVAL(laytrop_var_889)
  ixlow_var_918 = 0
  ixhigh_var_919 = 0
  ixc_var_917 = 0
  DO lay_var_905 = laytrop_min_var_915 + 1, laytrop_max_var_916
    icl_var_921 = 0
    ich_var_920 = 0
    DO jc_var_924 = kidia_var_871, kfdia_var_872
      IF (lay_var_905 <= laytrop_var_889(jc_var_924)) THEN
        icl_var_921 = icl_var_921 + 1
        ixlow_var_918(icl_var_921, lay_var_905) = jc_var_924
      ELSE
        ich_var_920 = ich_var_920 + 1
        ixhigh_var_919(ich_var_920, lay_var_905) = jc_var_924
      END IF
    END DO
    ixc_var_917(lay_var_905) = icl_var_921
  END DO
  DO lay_var_905 = 1, laytrop_min_var_915
    DO jl_var_925 = kidia_var_871, kfdia_var_872
      chi_co2_var_906 = colco2_var_887(jl_var_925, lay_var_905) / (coldry_var_888(jl_var_925, lay_var_905))
      ratco2_var_907 = 1D+20 * chi_co2_var_906 / chi_mls(2, jp_var_881(jl_var_925, lay_var_905) + 1)
      IF (ratco2_var_907 .GT. 3.0D0) THEN
        adjfac_var_908 = 2.0D0 + (ratco2_var_907 - 2.0D0) ** 0.65D0
        adjcolco2_var_909 = adjfac_var_908 * chi_mls(2, jp_var_881(jl_var_925, lay_var_905) + 1) * coldry_var_888(jl_var_925, lay_var_905) * 1D-20
      ELSE
        adjcolco2_var_909 = colco2_var_887(jl_var_925, lay_var_905)
      END IF
      ind0_var_899 = ((jp_var_881(jl_var_925, lay_var_905) - 1) * 5 + (jt_var_882(jl_var_925, lay_var_905) - 1)) * nspa_var_216(8) + 1
      ind1_var_900 = (jp_var_881(jl_var_925, lay_var_905) * 5 + (jt1_var_883(jl_var_925, lay_var_905) - 1)) * nspa_var_216(8) + 1
      inds_var_901 = indself_var_892(jl_var_925, lay_var_905)
      indf_var_902 = indfor_var_894(jl_var_925, lay_var_905)
      indm_var_903 = indminor_var_898(jl_var_925, lay_var_905)
      DO ig_var_904 = 1, 8
        tauself_var_911 = selffac_var_890(jl_var_925, lay_var_905) * (selfref_var_204(inds_var_901, ig_var_904) + selffrac_var_891(jl_var_925, lay_var_905) * (selfref_var_204(inds_var_901 + 1, ig_var_904) - selfref_var_204(inds_var_901, ig_var_904)))
        taufor_var_910 = forfac_var_896(jl_var_925, lay_var_905) * (forref_var_205(indf_var_902, ig_var_904) + forfrac_var_895(jl_var_925, lay_var_905) * (forref_var_205(indf_var_902 + 1, ig_var_904) - forref_var_205(indf_var_902, ig_var_904)))
        absco2_var_913 = (ka_mco2_var_199(indm_var_903, ig_var_904) + minorfrac_var_897(jl_var_925, lay_var_905) * (ka_mco2_var_199(indm_var_903 + 1, ig_var_904) - ka_mco2_var_199(indm_var_903, ig_var_904)))
        abso3_var_912 = (ka_mo3_var_201(indm_var_903, ig_var_904) + minorfrac_var_897(jl_var_925, lay_var_905) * (ka_mo3_var_201(indm_var_903 + 1, ig_var_904) - ka_mo3_var_201(indm_var_903, ig_var_904)))
        absn2o_var_914 = (ka_mn2o_var_200(indm_var_903, ig_var_904) + minorfrac_var_897(jl_var_925, lay_var_905) * (ka_mn2o_var_200(indm_var_903 + 1, ig_var_904) - ka_mn2o_var_200(indm_var_903, ig_var_904)))
        taug_var_874(jl_var_925, 88 + ig_var_904, lay_var_905) = colh2o_var_884(jl_var_925, lay_var_905) * (fac00_var_877(jl_var_925, lay_var_905) * absa_var_197(ind0_var_899, ig_var_904) + fac10_var_879(jl_var_925, lay_var_905) * absa_var_197(ind0_var_899 + 1, ig_var_904) + fac01_var_878(jl_var_925, lay_var_905) * absa_var_197(ind1_var_900, ig_var_904) + fac11_var_880(jl_var_925, lay_var_905) * absa_var_197(ind1_var_900 + 1, ig_var_904)) + tauself_var_911 + taufor_var_910 + adjcolco2_var_909 * absco2_var_913 + colo3_var_885(jl_var_925, lay_var_905) * abso3_var_912 + coln2o_var_886(jl_var_925, lay_var_905) * absn2o_var_914 + wx_var_875(jl_var_925, 3, lay_var_905) * cfc12_var_196(ig_var_904) + wx_var_875(jl_var_925, 4, lay_var_905) * cfc22adj(ig_var_904)
        fracs_var_893(jl_var_925, 88 + ig_var_904, lay_var_905) = fracrefa_var_194(ig_var_904)
      END DO
    END DO
  END DO
  DO lay_var_905 = laytrop_max_var_916 + 1, klev_var_873
    DO jl_var_925 = kidia_var_871, kfdia_var_872
      chi_co2_var_906 = colco2_var_887(jl_var_925, lay_var_905) / coldry_var_888(jl_var_925, lay_var_905)
      ratco2_var_907 = 1D+20 * chi_co2_var_906 / chi_mls(2, jp_var_881(jl_var_925, lay_var_905) + 1)
      IF (ratco2_var_907 .GT. 3.0D0) THEN
        adjfac_var_908 = 2.0D0 + (ratco2_var_907 - 2.0D0) ** 0.65D0
        adjcolco2_var_909 = adjfac_var_908 * chi_mls(2, jp_var_881(jl_var_925, lay_var_905) + 1) * coldry_var_888(jl_var_925, lay_var_905) * 1D-20
      ELSE
        adjcolco2_var_909 = colco2_var_887(jl_var_925, lay_var_905)
      END IF
      ind0_var_899 = ((jp_var_881(jl_var_925, lay_var_905) - 13) * 5 + (jt_var_882(jl_var_925, lay_var_905) - 1)) * nspb_var_217(8) + 1
      ind1_var_900 = ((jp_var_881(jl_var_925, lay_var_905) - 12) * 5 + (jt1_var_883(jl_var_925, lay_var_905) - 1)) * nspb_var_217(8) + 1
      indm_var_903 = indminor_var_898(jl_var_925, lay_var_905)
      DO ig_var_904 = 1, 8
        absco2_var_913 = (kb_mco2_var_202(indm_var_903, ig_var_904) + minorfrac_var_897(jl_var_925, lay_var_905) * (kb_mco2_var_202(indm_var_903 + 1, ig_var_904) - kb_mco2_var_202(indm_var_903, ig_var_904)))
        absn2o_var_914 = (kb_mn2o_var_203(indm_var_903, ig_var_904) + minorfrac_var_897(jl_var_925, lay_var_905) * (kb_mn2o_var_203(indm_var_903 + 1, ig_var_904) - kb_mn2o_var_203(indm_var_903, ig_var_904)))
        taug_var_874(jl_var_925, 88 + ig_var_904, lay_var_905) = colo3_var_885(jl_var_925, lay_var_905) * (fac00_var_877(jl_var_925, lay_var_905) * absb_var_198(ind0_var_899, ig_var_904) + fac10_var_879(jl_var_925, lay_var_905) * absb_var_198(ind0_var_899 + 1, ig_var_904) + fac01_var_878(jl_var_925, lay_var_905) * absb_var_198(ind1_var_900, ig_var_904) + fac11_var_880(jl_var_925, lay_var_905) * absb_var_198(ind1_var_900 + 1, ig_var_904)) + adjcolco2_var_909 * absco2_var_913 + coln2o_var_886(jl_var_925, lay_var_905) * absn2o_var_914 + wx_var_875(jl_var_925, 3, lay_var_905) * cfc12_var_196(ig_var_904) + wx_var_875(jl_var_925, 4, lay_var_905) * cfc22adj(ig_var_904)
        fracs_var_893(jl_var_925, 88 + ig_var_904, lay_var_905) = fracrefb_var_195(ig_var_904)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_916 /= laytrop_min_var_915) THEN
    DO lay_var_905 = laytrop_min_var_915 + 1, laytrop_max_var_916
      ixc0_var_922 = ixc_var_917(lay_var_905)
      DO ixp_var_923 = 1, ixc0_var_922
        jl_var_925 = ixlow_var_918(ixp_var_923, lay_var_905)
        chi_co2_var_906 = colco2_var_887(jl_var_925, lay_var_905) / (coldry_var_888(jl_var_925, lay_var_905))
        ratco2_var_907 = 1D+20 * chi_co2_var_906 / chi_mls(2, jp_var_881(jl_var_925, lay_var_905) + 1)
        IF (ratco2_var_907 .GT. 3.0D0) THEN
          adjfac_var_908 = 2.0D0 + (ratco2_var_907 - 2.0D0) ** 0.65D0
          adjcolco2_var_909 = adjfac_var_908 * chi_mls(2, jp_var_881(jl_var_925, lay_var_905) + 1) * coldry_var_888(jl_var_925, lay_var_905) * 1D-20
        ELSE
          adjcolco2_var_909 = colco2_var_887(jl_var_925, lay_var_905)
        END IF
        ind0_var_899 = ((jp_var_881(jl_var_925, lay_var_905) - 1) * 5 + (jt_var_882(jl_var_925, lay_var_905) - 1)) * nspa_var_216(8) + 1
        ind1_var_900 = (jp_var_881(jl_var_925, lay_var_905) * 5 + (jt1_var_883(jl_var_925, lay_var_905) - 1)) * nspa_var_216(8) + 1
        inds_var_901 = indself_var_892(jl_var_925, lay_var_905)
        indf_var_902 = indfor_var_894(jl_var_925, lay_var_905)
        indm_var_903 = indminor_var_898(jl_var_925, lay_var_905)
        DO ig_var_904 = 1, 8
          tauself_var_911 = selffac_var_890(jl_var_925, lay_var_905) * (selfref_var_204(inds_var_901, ig_var_904) + selffrac_var_891(jl_var_925, lay_var_905) * (selfref_var_204(inds_var_901 + 1, ig_var_904) - selfref_var_204(inds_var_901, ig_var_904)))
          taufor_var_910 = forfac_var_896(jl_var_925, lay_var_905) * (forref_var_205(indf_var_902, ig_var_904) + forfrac_var_895(jl_var_925, lay_var_905) * (forref_var_205(indf_var_902 + 1, ig_var_904) - forref_var_205(indf_var_902, ig_var_904)))
          absco2_var_913 = (ka_mco2_var_199(indm_var_903, ig_var_904) + minorfrac_var_897(jl_var_925, lay_var_905) * (ka_mco2_var_199(indm_var_903 + 1, ig_var_904) - ka_mco2_var_199(indm_var_903, ig_var_904)))
          abso3_var_912 = (ka_mo3_var_201(indm_var_903, ig_var_904) + minorfrac_var_897(jl_var_925, lay_var_905) * (ka_mo3_var_201(indm_var_903 + 1, ig_var_904) - ka_mo3_var_201(indm_var_903, ig_var_904)))
          absn2o_var_914 = (ka_mn2o_var_200(indm_var_903, ig_var_904) + minorfrac_var_897(jl_var_925, lay_var_905) * (ka_mn2o_var_200(indm_var_903 + 1, ig_var_904) - ka_mn2o_var_200(indm_var_903, ig_var_904)))
          taug_var_874(jl_var_925, 88 + ig_var_904, lay_var_905) = colh2o_var_884(jl_var_925, lay_var_905) * (fac00_var_877(jl_var_925, lay_var_905) * absa_var_197(ind0_var_899, ig_var_904) + fac10_var_879(jl_var_925, lay_var_905) * absa_var_197(ind0_var_899 + 1, ig_var_904) + fac01_var_878(jl_var_925, lay_var_905) * absa_var_197(ind1_var_900, ig_var_904) + fac11_var_880(jl_var_925, lay_var_905) * absa_var_197(ind1_var_900 + 1, ig_var_904)) + tauself_var_911 + taufor_var_910 + adjcolco2_var_909 * absco2_var_913 + colo3_var_885(jl_var_925, lay_var_905) * abso3_var_912 + coln2o_var_886(jl_var_925, lay_var_905) * absn2o_var_914 + wx_var_875(jl_var_925, 3, lay_var_905) * cfc12_var_196(ig_var_904) + wx_var_875(jl_var_925, 4, lay_var_905) * cfc22adj(ig_var_904)
          fracs_var_893(jl_var_925, 88 + ig_var_904, lay_var_905) = fracrefa_var_194(ig_var_904)
        END DO
      END DO
      ixc0_var_922 = kfdia_var_872 - kidia_var_871 + 1 - ixc0_var_922
      DO ixp_var_923 = 1, ixc0_var_922
        jl_var_925 = ixhigh_var_919(ixp_var_923, lay_var_905)
        chi_co2_var_906 = colco2_var_887(jl_var_925, lay_var_905) / coldry_var_888(jl_var_925, lay_var_905)
        ratco2_var_907 = 1D+20 * chi_co2_var_906 / chi_mls(2, jp_var_881(jl_var_925, lay_var_905) + 1)
        IF (ratco2_var_907 .GT. 3.0D0) THEN
          adjfac_var_908 = 2.0D0 + (ratco2_var_907 - 2.0D0) ** 0.65D0
          adjcolco2_var_909 = adjfac_var_908 * chi_mls(2, jp_var_881(jl_var_925, lay_var_905) + 1) * coldry_var_888(jl_var_925, lay_var_905) * 1D-20
        ELSE
          adjcolco2_var_909 = colco2_var_887(jl_var_925, lay_var_905)
        END IF
        ind0_var_899 = ((jp_var_881(jl_var_925, lay_var_905) - 13) * 5 + (jt_var_882(jl_var_925, lay_var_905) - 1)) * nspb_var_217(8) + 1
        ind1_var_900 = ((jp_var_881(jl_var_925, lay_var_905) - 12) * 5 + (jt1_var_883(jl_var_925, lay_var_905) - 1)) * nspb_var_217(8) + 1
        indm_var_903 = indminor_var_898(jl_var_925, lay_var_905)
        DO ig_var_904 = 1, 8
          absco2_var_913 = (kb_mco2_var_202(indm_var_903, ig_var_904) + minorfrac_var_897(jl_var_925, lay_var_905) * (kb_mco2_var_202(indm_var_903 + 1, ig_var_904) - kb_mco2_var_202(indm_var_903, ig_var_904)))
          absn2o_var_914 = (kb_mn2o_var_203(indm_var_903, ig_var_904) + minorfrac_var_897(jl_var_925, lay_var_905) * (kb_mn2o_var_203(indm_var_903 + 1, ig_var_904) - kb_mn2o_var_203(indm_var_903, ig_var_904)))
          taug_var_874(jl_var_925, 88 + ig_var_904, lay_var_905) = colo3_var_885(jl_var_925, lay_var_905) * (fac00_var_877(jl_var_925, lay_var_905) * absb_var_198(ind0_var_899, ig_var_904) + fac10_var_879(jl_var_925, lay_var_905) * absb_var_198(ind0_var_899 + 1, ig_var_904) + fac01_var_878(jl_var_925, lay_var_905) * absb_var_198(ind1_var_900, ig_var_904) + fac11_var_880(jl_var_925, lay_var_905) * absb_var_198(ind1_var_900 + 1, ig_var_904)) + adjcolco2_var_909 * absco2_var_913 + coln2o_var_886(jl_var_925, lay_var_905) * absn2o_var_914 + wx_var_875(jl_var_925, 3, lay_var_905) * cfc12_var_196(ig_var_904) + wx_var_875(jl_var_925, 4, lay_var_905) * cfc22adj(ig_var_904)
          fracs_var_893(jl_var_925, 88 + ig_var_904, lay_var_905) = fracrefb_var_195(ig_var_904)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol8
SUBROUTINE rrtm_prepare_gases(kidia_var_927, kfdia_var_928, klon, klev_var_926, paph, pap, pth, pt, pq, pco2, pch4, pn2o, pno2, pc11, pc12, pc22, pcl4, pozn, pcoldry_var_929, pwbrodl, pwkl_var_930, pwx_var_931, pavel_var_932, ptavel_var_933, pz, ptz, kreflect)
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: klon
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_926
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_927, kfdia_var_928
  REAL(KIND = 8), INTENT(IN) :: paph(klon, klev_var_926 + 1)
  REAL(KIND = 8), INTENT(IN) :: pap(klon, klev_var_926)
  REAL(KIND = 8), INTENT(IN) :: pth(klon, klev_var_926 + 1)
  REAL(KIND = 8), INTENT(IN) :: pt(klon, klev_var_926)
  REAL(KIND = 8), INTENT(IN) :: pq(klon, klev_var_926)
  REAL(KIND = 8), INTENT(IN) :: pco2(klon, klev_var_926)
  REAL(KIND = 8), INTENT(IN) :: pch4(klon, klev_var_926)
  REAL(KIND = 8), INTENT(IN) :: pn2o(klon, klev_var_926)
  REAL(KIND = 8), INTENT(IN) :: pno2(klon, klev_var_926)
  REAL(KIND = 8), INTENT(IN) :: pc11(klon, klev_var_926)
  REAL(KIND = 8), INTENT(IN) :: pc12(klon, klev_var_926)
  REAL(KIND = 8), INTENT(IN) :: pc22(klon, klev_var_926)
  REAL(KIND = 8), INTENT(IN) :: pcl4(klon, klev_var_926)
  REAL(KIND = 8), INTENT(IN) :: pozn(klon, klev_var_926)
  REAL(KIND = 8), INTENT(OUT) :: pcoldry_var_929(kidia_var_927 : kfdia_var_928, klev_var_926)
  REAL(KIND = 8), INTENT(OUT) :: pwbrodl(kidia_var_927 : kfdia_var_928, klev_var_926)
  REAL(KIND = 8), INTENT(OUT) :: pwkl_var_930(kidia_var_927 : kfdia_var_928, 35, klev_var_926)
  REAL(KIND = 8), INTENT(OUT) :: pwx_var_931(kidia_var_927 : kfdia_var_928, 4, klev_var_926)
  REAL(KIND = 8), INTENT(OUT) :: pavel_var_932(kidia_var_927 : kfdia_var_928, klev_var_926)
  REAL(KIND = 8), INTENT(OUT) :: ptavel_var_933(kidia_var_927 : kfdia_var_928, klev_var_926)
  REAL(KIND = 8), INTENT(OUT) :: pz(kidia_var_927 : kfdia_var_928, 0 : klev_var_926)
  REAL(KIND = 8), INTENT(OUT) :: ptz(kidia_var_927 : kfdia_var_928, 0 : klev_var_926)
  INTEGER(KIND = 4), INTENT(OUT) :: kreflect(kidia_var_927 : kfdia_var_928)
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
  INTEGER(KIND = 4) :: iatm, jmol, ixmax, j1, j2, jk_var_934, jl_var_935
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
  DO jl_var_935 = kidia_var_927, kfdia_var_928
    kreflect(jl_var_935) = 0
  END DO
  DO j2 = 1, klev_var_926
    DO j1 = 1, 35
      DO jl_var_935 = kidia_var_927, kfdia_var_928
        pwkl_var_930(jl_var_935, j1, j2) = 0.0D0
      END DO
    END DO
  END DO
  iatm = 0
  ixmax = 4
  DO jl_var_935 = kidia_var_927, kfdia_var_928
    pz(jl_var_935, 0) = paph(jl_var_935, klev_var_926 + 1) / 100.0D0
    ptz(jl_var_935, 0) = pth(jl_var_935, klev_var_926 + 1)
  END DO
  DO jk_var_934 = 1, klev_var_926
    DO jl_var_935 = kidia_var_927, kfdia_var_928
      pavel_var_932(jl_var_935, jk_var_934) = pap(jl_var_935, klev_var_926 - jk_var_934 + 1) / 100.0D0
      ptavel_var_933(jl_var_935, jk_var_934) = pt(jl_var_935, klev_var_926 - jk_var_934 + 1)
      pz(jl_var_935, jk_var_934) = paph(jl_var_935, klev_var_926 - jk_var_934 + 1) / 100.0D0
      ptz(jl_var_935, jk_var_934) = pth(jl_var_935, klev_var_926 - jk_var_934 + 1)
      pwkl_var_930(jl_var_935, 1, jk_var_934) = MAX(pq(jl_var_935, klev_var_926 - jk_var_934 + 1), 1E-15) * zamd / zamw
      pwkl_var_930(jl_var_935, 2, jk_var_934) = pco2(jl_var_935, klev_var_926 - jk_var_934 + 1) * zamd / zamco2
      pwkl_var_930(jl_var_935, 3, jk_var_934) = pozn(jl_var_935, klev_var_926 - jk_var_934 + 1) * zamd / zamo
      pwkl_var_930(jl_var_935, 4, jk_var_934) = pn2o(jl_var_935, klev_var_926 - jk_var_934 + 1) * zamd / zamn2o
      pwkl_var_930(jl_var_935, 6, jk_var_934) = pch4(jl_var_935, klev_var_926 - jk_var_934 + 1) * zamd / zamch4
      pwkl_var_930(jl_var_935, 7, jk_var_934) = 0.209488D0
      zamm = (1.0D0 - pwkl_var_930(jl_var_935, 1, jk_var_934)) * zamd + pwkl_var_930(jl_var_935, 1, jk_var_934) * zamw
      pcoldry_var_929(jl_var_935, jk_var_934) = (pz(jl_var_935, jk_var_934 - 1) - pz(jl_var_935, jk_var_934)) * 1000.0D0 * zavgdro / (zgravit * zamm * (1.0D0 + pwkl_var_930(jl_var_935, 1, jk_var_934)))
    END DO
  END DO
  DO j2 = 1, klev_var_926
    DO j1 = 1, 4
      DO jl_var_935 = kidia_var_927, kfdia_var_928
        pwx_var_931(jl_var_935, j1, j2) = 0.0D0
      END DO
    END DO
  END DO
  DO jk_var_934 = 1, klev_var_926
    DO jl_var_935 = kidia_var_927, kfdia_var_928
      pwx_var_931(jl_var_935, 1, jk_var_934) = pcl4(jl_var_935, klev_var_926 - jk_var_934 + 1) * zamd / zamcl4
      pwx_var_931(jl_var_935, 2, jk_var_934) = pc11(jl_var_935, klev_var_926 - jk_var_934 + 1) * zamd / zamc11
      pwx_var_931(jl_var_935, 3, jk_var_934) = pc12(jl_var_935, klev_var_926 - jk_var_934 + 1) * zamd / zamc12
      pwx_var_931(jl_var_935, 4, jk_var_934) = pc22(jl_var_935, klev_var_926 - jk_var_934 + 1) * zamd / zamc22
      pwx_var_931(jl_var_935, 1, jk_var_934) = pcoldry_var_929(jl_var_935, jk_var_934) * pwx_var_931(jl_var_935, 1, jk_var_934) * 1D-20
      pwx_var_931(jl_var_935, 2, jk_var_934) = pcoldry_var_929(jl_var_935, jk_var_934) * pwx_var_931(jl_var_935, 2, jk_var_934) * 1D-20
      pwx_var_931(jl_var_935, 3, jk_var_934) = pcoldry_var_929(jl_var_935, jk_var_934) * pwx_var_931(jl_var_935, 3, jk_var_934) * 1D-20
      pwx_var_931(jl_var_935, 4, jk_var_934) = pcoldry_var_929(jl_var_935, jk_var_934) * pwx_var_931(jl_var_935, 4, jk_var_934) * 1D-20
      zsummol = 0.0D0
      DO jmol = 2, 7
        zsummol = zsummol + pwkl_var_930(jl_var_935, jmol, jk_var_934)
      END DO
      pwbrodl(jl_var_935, jk_var_934) = pcoldry_var_929(jl_var_935, jk_var_934) * (1.0D0 - zsummol)
      DO jmol = 1, 7
        pwkl_var_930(jl_var_935, jmol, jk_var_934) = pcoldry_var_929(jl_var_935, jk_var_934) * pwkl_var_930(jl_var_935, jmol, jk_var_934)
      END DO
    END DO
  END DO
END SUBROUTINE rrtm_prepare_gases
SUBROUTINE srtm_gas_optical_depth(kidia_var_936, kfdia_var_937, klev_var_938, poneminus_var_939, prmu0_var_940, klaytrop_var_941, pcolch4_var_942, pcolco2_var_943, pcolh2o_var_944, pcolmol_var_945, pcolo2_var_946, pcolo3_var_947, pforfac_var_948, pforfrac_var_949, kindfor_var_950, pselffac_var_951, pselffrac_var_952, kindself_var_953, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, pod_var_961, pssa, pincsol)
  USE yoesrtwn, ONLY: ngc
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_936, kfdia_var_937
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_938
  REAL(KIND = 8), INTENT(IN) :: poneminus_var_939(kidia_var_936 : kfdia_var_937)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_940(kidia_var_936 : kfdia_var_937)
  INTEGER(KIND = 4), INTENT(IN) :: klaytrop_var_941(kidia_var_936 : kfdia_var_937)
  REAL(KIND = 8), INTENT(IN) :: pcolch4_var_942(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pcolco2_var_943(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pcolh2o_var_944(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pcolmol_var_945(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pcolo2_var_946(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pcolo3_var_947(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pforfac_var_948(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pforfrac_var_949(kidia_var_936 : kfdia_var_937, klev_var_938)
  INTEGER(KIND = 4), INTENT(IN) :: kindfor_var_950(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pselffac_var_951(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pselffrac_var_952(kidia_var_936 : kfdia_var_937, klev_var_938)
  INTEGER(KIND = 4), INTENT(IN) :: kindself_var_953(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pfac00_var_954(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pfac01_var_955(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pfac10_var_956(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(IN) :: pfac11_var_957(kidia_var_936 : kfdia_var_937, klev_var_938)
  INTEGER(KIND = 4), INTENT(IN) :: kjp_var_958(kidia_var_936 : kfdia_var_937, klev_var_938)
  INTEGER(KIND = 4), INTENT(IN) :: kjt_var_959(kidia_var_936 : kfdia_var_937, klev_var_938)
  INTEGER(KIND = 4), INTENT(IN) :: kjt1_var_960(kidia_var_936 : kfdia_var_937, klev_var_938)
  REAL(KIND = 8), INTENT(INOUT) :: pod_var_961(kidia_var_936 : kfdia_var_937, klev_var_938, 112)
  REAL(KIND = 8), INTENT(INOUT) :: pssa(kidia_var_936 : kfdia_var_937, klev_var_938, 112)
  REAL(KIND = 8), INTENT(INOUT) :: pincsol(kidia_var_936 : kfdia_var_937, 112)
  INTEGER(KIND = 4) :: ib1, ib2, ibm, igt, iw(kidia_var_936 : kfdia_var_937), jb_var_962, jg_var_963, jk_var_964, jl_var_965, icount
  REAL(KIND = 8) :: ztaug(kidia_var_936 : kfdia_var_937, klev_var_938, 16)
  REAL(KIND = 8) :: ztaur(kidia_var_936 : kfdia_var_937, klev_var_938, 16)
  REAL(KIND = 8) :: zsflxzen(kidia_var_936 : kfdia_var_937, 16)
  ib1 = 16
  ib2 = 29
  icount = 0
  DO jl_var_965 = kidia_var_936, kfdia_var_937
    IF (prmu0_var_940(jl_var_965) > 0.0D0) THEN
      icount = icount + 1
      iw(jl_var_965) = 0
    END IF
  END DO
  IF (icount /= 0) THEN
    DO jb_var_962 = ib1, ib2
      ibm = jb_var_962 - 15
      igt = ngc(ibm)
      IF (jb_var_962 == 16) THEN
        CALL srtm_taumol16(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, poneminus_var_939, pcolh2o_var_944, pcolch4_var_942, pcolmol_var_945, klaytrop_var_941, pselffac_var_951, pselffrac_var_952, kindself_var_953, pforfac_var_948, pforfrac_var_949, kindfor_var_950, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 17) THEN
        CALL srtm_taumol17(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, poneminus_var_939, pcolh2o_var_944, pcolco2_var_943, pcolmol_var_945, klaytrop_var_941, pselffac_var_951, pselffrac_var_952, kindself_var_953, pforfac_var_948, pforfrac_var_949, kindfor_var_950, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 18) THEN
        CALL srtm_taumol18(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, poneminus_var_939, pcolh2o_var_944, pcolch4_var_942, pcolmol_var_945, klaytrop_var_941, pselffac_var_951, pselffrac_var_952, kindself_var_953, pforfac_var_948, pforfrac_var_949, kindfor_var_950, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 19) THEN
        CALL srtm_taumol19(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, poneminus_var_939, pcolh2o_var_944, pcolco2_var_943, pcolmol_var_945, klaytrop_var_941, pselffac_var_951, pselffrac_var_952, kindself_var_953, pforfac_var_948, pforfrac_var_949, kindfor_var_950, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 20) THEN
        CALL srtm_taumol20(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, pcolh2o_var_944, pcolch4_var_942, pcolmol_var_945, klaytrop_var_941, pselffac_var_951, pselffrac_var_952, kindself_var_953, pforfac_var_948, pforfrac_var_949, kindfor_var_950, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 21) THEN
        CALL srtm_taumol21(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, poneminus_var_939, pcolh2o_var_944, pcolco2_var_943, pcolmol_var_945, klaytrop_var_941, pselffac_var_951, pselffrac_var_952, kindself_var_953, pforfac_var_948, pforfrac_var_949, kindfor_var_950, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 22) THEN
        CALL srtm_taumol22(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, poneminus_var_939, pcolh2o_var_944, pcolmol_var_945, pcolo2_var_946, klaytrop_var_941, pselffac_var_951, pselffrac_var_952, kindself_var_953, pforfac_var_948, pforfrac_var_949, kindfor_var_950, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 23) THEN
        CALL srtm_taumol23(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, pcolh2o_var_944, pcolmol_var_945, klaytrop_var_941, pselffac_var_951, pselffrac_var_952, kindself_var_953, pforfac_var_948, pforfrac_var_949, kindfor_var_950, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 24) THEN
        CALL srtm_taumol24(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, poneminus_var_939, pcolh2o_var_944, pcolmol_var_945, pcolo2_var_946, pcolo3_var_947, klaytrop_var_941, pselffac_var_951, pselffrac_var_952, kindself_var_953, pforfac_var_948, pforfrac_var_949, kindfor_var_950, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 25) THEN
        CALL srtm_taumol25(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, pcolh2o_var_944, pcolmol_var_945, pcolo3_var_947, klaytrop_var_941, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 26) THEN
        CALL srtm_taumol26(kidia_var_936, kfdia_var_937, klev_var_938, pcolmol_var_945, klaytrop_var_941, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 27) THEN
        CALL srtm_taumol27(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, pcolmol_var_945, pcolo3_var_947, klaytrop_var_941, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 28) THEN
        CALL srtm_taumol28(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, poneminus_var_939, pcolmol_var_945, pcolo2_var_946, pcolo3_var_947, klaytrop_var_941, zsflxzen, ztaug, ztaur, prmu0_var_940)
      ELSE IF (jb_var_962 == 29) THEN
        CALL srtm_taumol29(kidia_var_936, kfdia_var_937, klev_var_938, pfac00_var_954, pfac01_var_955, pfac10_var_956, pfac11_var_957, kjp_var_958, kjt_var_959, kjt1_var_960, pcolh2o_var_944, pcolco2_var_943, pcolmol_var_945, klaytrop_var_941, pselffac_var_951, pselffrac_var_952, kindself_var_953, pforfac_var_948, pforfrac_var_949, kindfor_var_950, zsflxzen, ztaug, ztaur, prmu0_var_940)
      END IF
      DO jg_var_963 = 1, igt
        DO jl_var_965 = kidia_var_936, kfdia_var_937
          IF (prmu0_var_940(jl_var_965) > 0.0D0) THEN
            iw(jl_var_965) = iw(jl_var_965) + 1
            pincsol(jl_var_965, iw(jl_var_965)) = zsflxzen(jl_var_965, jg_var_963)
          END IF
        END DO
        DO jk_var_964 = 1, klev_var_938
          DO jl_var_965 = kidia_var_936, kfdia_var_937
            IF (prmu0_var_940(jl_var_965) > 0.0D0) THEN
              pod_var_961(jl_var_965, jk_var_964, iw(jl_var_965)) = ztaur(jl_var_965, jk_var_964, jg_var_963) + ztaug(jl_var_965, jk_var_964, jg_var_963)
              pssa(jl_var_965, jk_var_964, iw(jl_var_965)) = ztaur(jl_var_965, jk_var_964, jg_var_963) / pod_var_961(jl_var_965, jk_var_964, iw(jl_var_965))
            END IF
          END DO
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE srtm_gas_optical_depth
SUBROUTINE srtm_taumol28(kidia_var_966, kfdia_var_967, klev_var_968, p_fac00_var_969, p_fac01_var_970, p_fac10_var_971, p_fac11_var_972, k_jp_var_973, k_jt_var_974, k_jt1_var_975, p_oneminus_var_976, p_colmol_var_977, p_colo2_var_978, p_colo3_var_979, k_laytrop_var_980, p_sfluxzen_var_981, p_taug_var_982, p_taur_var_983, prmu0_var_984)
  USE yoesrta28, ONLY: absa_var_303, absb_var_304, layreffr_var_302, rayl_var_300, sfluxrefc_var_305, strrat_var_301
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_966, kfdia_var_967
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_968
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_969(kidia_var_966 : kfdia_var_967, klev_var_968)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_970(kidia_var_966 : kfdia_var_967, klev_var_968)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_971(kidia_var_966 : kfdia_var_967, klev_var_968)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_972(kidia_var_966 : kfdia_var_967, klev_var_968)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_973(kidia_var_966 : kfdia_var_967, klev_var_968)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_974(kidia_var_966 : kfdia_var_967, klev_var_968)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_975(kidia_var_966 : kfdia_var_967, klev_var_968)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_976(kidia_var_966 : kfdia_var_967)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_977(kidia_var_966 : kfdia_var_967, klev_var_968)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_978(kidia_var_966 : kfdia_var_967, klev_var_968)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_979(kidia_var_966 : kfdia_var_967, klev_var_968)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_980(kidia_var_966 : kfdia_var_967)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_981(kidia_var_966 : kfdia_var_967, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_982(kidia_var_966 : kfdia_var_967, klev_var_968, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_983(kidia_var_966 : kfdia_var_967, klev_var_968, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_984(kidia_var_966 : kfdia_var_967)
  INTEGER(KIND = 4) :: ig_var_985, ind0_var_986, ind1_var_987, js_var_988, i_lay_var_989, i_laysolfr_var_990(kidia_var_966 : kfdia_var_967), i_nlayers_var_991, iplon_var_992
  INTEGER(KIND = 4) :: laytrop_min_var_993, laytrop_max_var_994
  REAL(KIND = 8) :: z_fs_var_995, z_speccomb_var_996, z_specmult_var_997, z_specparm_var_998, z_tauray_var_999
  laytrop_min_var_993 = MINVAL(k_laytrop_var_980(kidia_var_966 : kfdia_var_967))
  laytrop_max_var_994 = MAXVAL(k_laytrop_var_980(kidia_var_966 : kfdia_var_967))
  i_nlayers_var_991 = klev_var_968
  DO iplon_var_992 = kidia_var_966, kfdia_var_967
    i_laysolfr_var_990(iplon_var_992) = i_nlayers_var_991
  END DO
  DO i_lay_var_989 = 1, laytrop_min_var_993
    DO iplon_var_992 = kidia_var_966, kfdia_var_967
      z_speccomb_var_996 = p_colo3_var_979(iplon_var_992, i_lay_var_989) + strrat_var_301 * p_colo2_var_978(iplon_var_992, i_lay_var_989)
      z_specparm_var_998 = p_colo3_var_979(iplon_var_992, i_lay_var_989) / z_speccomb_var_996
      z_specparm_var_998 = MIN(p_oneminus_var_976(iplon_var_992), z_specparm_var_998)
      z_specmult_var_997 = 8.0D0 * (z_specparm_var_998)
      js_var_988 = 1 + INT(z_specmult_var_997)
      z_fs_var_995 = z_specmult_var_997 - AINT(z_specmult_var_997)
      ind0_var_986 = ((k_jp_var_973(iplon_var_992, i_lay_var_989) - 1) * 5 + (k_jt_var_974(iplon_var_992, i_lay_var_989) - 1)) * nspa_var_313(28) + js_var_988
      ind1_var_987 = (k_jp_var_973(iplon_var_992, i_lay_var_989) * 5 + (k_jt1_var_975(iplon_var_992, i_lay_var_989) - 1)) * nspa_var_313(28) + js_var_988
      z_tauray_var_999 = p_colmol_var_977(iplon_var_992, i_lay_var_989) * rayl_var_300
      DO ig_var_985 = 1, 6
        p_taug_var_982(iplon_var_992, i_lay_var_989, ig_var_985) = z_speccomb_var_996 * ((1.0D0 - z_fs_var_995) * (absa_var_303(ind0_var_986, ig_var_985) * p_fac00_var_969(iplon_var_992, i_lay_var_989) + absa_var_303(ind0_var_986 + 9, ig_var_985) * p_fac10_var_971(iplon_var_992, i_lay_var_989) + absa_var_303(ind1_var_987, ig_var_985) * p_fac01_var_970(iplon_var_992, i_lay_var_989) + absa_var_303(ind1_var_987 + 9, ig_var_985) * p_fac11_var_972(iplon_var_992, i_lay_var_989)) + z_fs_var_995 * (absa_var_303(ind0_var_986 + 1, ig_var_985) * p_fac00_var_969(iplon_var_992, i_lay_var_989) + absa_var_303(ind0_var_986 + 10, ig_var_985) * p_fac10_var_971(iplon_var_992, i_lay_var_989) + absa_var_303(ind1_var_987 + 1, ig_var_985) * p_fac01_var_970(iplon_var_992, i_lay_var_989) + absa_var_303(ind1_var_987 + 10, ig_var_985) * p_fac11_var_972(iplon_var_992, i_lay_var_989)))
        p_taur_var_983(iplon_var_992, i_lay_var_989, ig_var_985) = z_tauray_var_999
      END DO
    END DO
  END DO
  DO i_lay_var_989 = laytrop_min_var_993 + 1, laytrop_max_var_994
    DO iplon_var_992 = kidia_var_966, kfdia_var_967
      IF (i_lay_var_989 <= k_laytrop_var_980(iplon_var_992)) THEN
        z_speccomb_var_996 = p_colo3_var_979(iplon_var_992, i_lay_var_989) + strrat_var_301 * p_colo2_var_978(iplon_var_992, i_lay_var_989)
        z_specparm_var_998 = p_colo3_var_979(iplon_var_992, i_lay_var_989) / z_speccomb_var_996
        z_specparm_var_998 = MIN(p_oneminus_var_976(iplon_var_992), z_specparm_var_998)
        z_specmult_var_997 = 8.0D0 * (z_specparm_var_998)
        js_var_988 = 1 + INT(z_specmult_var_997)
        z_fs_var_995 = z_specmult_var_997 - AINT(z_specmult_var_997)
        ind0_var_986 = ((k_jp_var_973(iplon_var_992, i_lay_var_989) - 1) * 5 + (k_jt_var_974(iplon_var_992, i_lay_var_989) - 1)) * nspa_var_313(28) + js_var_988
        ind1_var_987 = (k_jp_var_973(iplon_var_992, i_lay_var_989) * 5 + (k_jt1_var_975(iplon_var_992, i_lay_var_989) - 1)) * nspa_var_313(28) + js_var_988
        z_tauray_var_999 = p_colmol_var_977(iplon_var_992, i_lay_var_989) * rayl_var_300
        DO ig_var_985 = 1, 6
          p_taug_var_982(iplon_var_992, i_lay_var_989, ig_var_985) = z_speccomb_var_996 * ((1.0D0 - z_fs_var_995) * (absa_var_303(ind0_var_986, ig_var_985) * p_fac00_var_969(iplon_var_992, i_lay_var_989) + absa_var_303(ind0_var_986 + 9, ig_var_985) * p_fac10_var_971(iplon_var_992, i_lay_var_989) + absa_var_303(ind1_var_987, ig_var_985) * p_fac01_var_970(iplon_var_992, i_lay_var_989) + absa_var_303(ind1_var_987 + 9, ig_var_985) * p_fac11_var_972(iplon_var_992, i_lay_var_989)) + z_fs_var_995 * (absa_var_303(ind0_var_986 + 1, ig_var_985) * p_fac00_var_969(iplon_var_992, i_lay_var_989) + absa_var_303(ind0_var_986 + 10, ig_var_985) * p_fac10_var_971(iplon_var_992, i_lay_var_989) + absa_var_303(ind1_var_987 + 1, ig_var_985) * p_fac01_var_970(iplon_var_992, i_lay_var_989) + absa_var_303(ind1_var_987 + 10, ig_var_985) * p_fac11_var_972(iplon_var_992, i_lay_var_989)))
          p_taur_var_983(iplon_var_992, i_lay_var_989, ig_var_985) = z_tauray_var_999
        END DO
      ELSE
        IF (k_jp_var_973(iplon_var_992, i_lay_var_989 - 1) < layreffr_var_302 .AND. k_jp_var_973(iplon_var_992, i_lay_var_989) >= layreffr_var_302) i_laysolfr_var_990(iplon_var_992) = i_lay_var_989
        z_speccomb_var_996 = p_colo3_var_979(iplon_var_992, i_lay_var_989) + strrat_var_301 * p_colo2_var_978(iplon_var_992, i_lay_var_989)
        z_specparm_var_998 = p_colo3_var_979(iplon_var_992, i_lay_var_989) / z_speccomb_var_996
        z_specparm_var_998 = MIN(p_oneminus_var_976(iplon_var_992), z_specparm_var_998)
        z_specmult_var_997 = 4.0D0 * (z_specparm_var_998)
        js_var_988 = 1 + INT(z_specmult_var_997)
        z_fs_var_995 = z_specmult_var_997 - AINT(z_specmult_var_997)
        ind0_var_986 = ((k_jp_var_973(iplon_var_992, i_lay_var_989) - 13) * 5 + (k_jt_var_974(iplon_var_992, i_lay_var_989) - 1)) * nspb_var_314(28) + js_var_988
        ind1_var_987 = ((k_jp_var_973(iplon_var_992, i_lay_var_989) - 12) * 5 + (k_jt1_var_975(iplon_var_992, i_lay_var_989) - 1)) * nspb_var_314(28) + js_var_988
        z_tauray_var_999 = p_colmol_var_977(iplon_var_992, i_lay_var_989) * rayl_var_300
        DO ig_var_985 = 1, 6
          p_taug_var_982(iplon_var_992, i_lay_var_989, ig_var_985) = z_speccomb_var_996 * ((1.0D0 - z_fs_var_995) * (absb_var_304(ind0_var_986, ig_var_985) * p_fac00_var_969(iplon_var_992, i_lay_var_989) + absb_var_304(ind0_var_986 + 5, ig_var_985) * p_fac10_var_971(iplon_var_992, i_lay_var_989) + absb_var_304(ind1_var_987, ig_var_985) * p_fac01_var_970(iplon_var_992, i_lay_var_989) + absb_var_304(ind1_var_987 + 5, ig_var_985) * p_fac11_var_972(iplon_var_992, i_lay_var_989)) + z_fs_var_995 * (absb_var_304(ind0_var_986 + 1, ig_var_985) * p_fac00_var_969(iplon_var_992, i_lay_var_989) + absb_var_304(ind0_var_986 + 6, ig_var_985) * p_fac10_var_971(iplon_var_992, i_lay_var_989) + absb_var_304(ind1_var_987 + 1, ig_var_985) * p_fac01_var_970(iplon_var_992, i_lay_var_989) + absb_var_304(ind1_var_987 + 6, ig_var_985) * p_fac11_var_972(iplon_var_992, i_lay_var_989)))
          IF (i_lay_var_989 == i_laysolfr_var_990(iplon_var_992)) p_sfluxzen_var_981(iplon_var_992, ig_var_985) = sfluxrefc_var_305(ig_var_985, js_var_988) + z_fs_var_995 * (sfluxrefc_var_305(ig_var_985, js_var_988 + 1) - sfluxrefc_var_305(ig_var_985, js_var_988))
          p_taur_var_983(iplon_var_992, i_lay_var_989, ig_var_985) = z_tauray_var_999
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_989 = laytrop_max_var_994 + 1, i_nlayers_var_991
    DO iplon_var_992 = kidia_var_966, kfdia_var_967
      IF (k_jp_var_973(iplon_var_992, i_lay_var_989 - 1) < layreffr_var_302 .AND. k_jp_var_973(iplon_var_992, i_lay_var_989) >= layreffr_var_302) i_laysolfr_var_990(iplon_var_992) = i_lay_var_989
      z_speccomb_var_996 = p_colo3_var_979(iplon_var_992, i_lay_var_989) + strrat_var_301 * p_colo2_var_978(iplon_var_992, i_lay_var_989)
      z_specparm_var_998 = p_colo3_var_979(iplon_var_992, i_lay_var_989) / z_speccomb_var_996
      z_specparm_var_998 = MIN(p_oneminus_var_976(iplon_var_992), z_specparm_var_998)
      z_specmult_var_997 = 4.0D0 * (z_specparm_var_998)
      js_var_988 = 1 + INT(z_specmult_var_997)
      z_fs_var_995 = z_specmult_var_997 - AINT(z_specmult_var_997)
      ind0_var_986 = ((k_jp_var_973(iplon_var_992, i_lay_var_989) - 13) * 5 + (k_jt_var_974(iplon_var_992, i_lay_var_989) - 1)) * nspb_var_314(28) + js_var_988
      ind1_var_987 = ((k_jp_var_973(iplon_var_992, i_lay_var_989) - 12) * 5 + (k_jt1_var_975(iplon_var_992, i_lay_var_989) - 1)) * nspb_var_314(28) + js_var_988
      z_tauray_var_999 = p_colmol_var_977(iplon_var_992, i_lay_var_989) * rayl_var_300
      DO ig_var_985 = 1, 6
        p_taug_var_982(iplon_var_992, i_lay_var_989, ig_var_985) = z_speccomb_var_996 * ((1.0D0 - z_fs_var_995) * (absb_var_304(ind0_var_986, ig_var_985) * p_fac00_var_969(iplon_var_992, i_lay_var_989) + absb_var_304(ind0_var_986 + 5, ig_var_985) * p_fac10_var_971(iplon_var_992, i_lay_var_989) + absb_var_304(ind1_var_987, ig_var_985) * p_fac01_var_970(iplon_var_992, i_lay_var_989) + absb_var_304(ind1_var_987 + 5, ig_var_985) * p_fac11_var_972(iplon_var_992, i_lay_var_989)) + z_fs_var_995 * (absb_var_304(ind0_var_986 + 1, ig_var_985) * p_fac00_var_969(iplon_var_992, i_lay_var_989) + absb_var_304(ind0_var_986 + 6, ig_var_985) * p_fac10_var_971(iplon_var_992, i_lay_var_989) + absb_var_304(ind1_var_987 + 1, ig_var_985) * p_fac01_var_970(iplon_var_992, i_lay_var_989) + absb_var_304(ind1_var_987 + 6, ig_var_985) * p_fac11_var_972(iplon_var_992, i_lay_var_989)))
        IF (i_lay_var_989 == i_laysolfr_var_990(iplon_var_992)) p_sfluxzen_var_981(iplon_var_992, ig_var_985) = sfluxrefc_var_305(ig_var_985, js_var_988) + z_fs_var_995 * (sfluxrefc_var_305(ig_var_985, js_var_988 + 1) - sfluxrefc_var_305(ig_var_985, js_var_988))
        p_taur_var_983(iplon_var_992, i_lay_var_989, ig_var_985) = z_tauray_var_999
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol28
SUBROUTINE srtm_taumol29(kidia_var_1000, kfdia_var_1001, klev_var_1002, p_fac00_var_1003, p_fac01_var_1004, p_fac10_var_1005, p_fac11_var_1006, k_jp_var_1007, k_jt_var_1008, k_jt1_var_1009, p_colh2o_var_1010, p_colco2_var_1011, p_colmol_var_1012, k_laytrop_var_1013, p_selffac_var_1014, p_selffrac_var_1015, k_indself_var_1016, p_forfac_var_1017, p_forfrac_var_1018, k_indfor_var_1019, p_sfluxzen_var_1022, p_taug_var_1023, p_taur_var_1024, prmu0_var_1025)
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  USE yoesrta29, ONLY: absa_var_308, absb_var_309, absco2c, absh2oc, forrefc_var_311, layreffr_var_307, rayl_var_306, selfrefc_var_310, sfluxrefc_var_312
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1000, kfdia_var_1001
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1002
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1003(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1004(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1005(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1006(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1007(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1008(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1009(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1010(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1011(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1012(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1013(kidia_var_1000 : kfdia_var_1001)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1014(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1015(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1016(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1017(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1018(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1019(kidia_var_1000 : kfdia_var_1001, klev_var_1002)
  INTEGER(KIND = 4) :: laytrop_min_var_1020, laytrop_max_var_1021
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1022(kidia_var_1000 : kfdia_var_1001, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1023(kidia_var_1000 : kfdia_var_1001, klev_var_1002, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1024(kidia_var_1000 : kfdia_var_1001, klev_var_1002, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1025(kidia_var_1000 : kfdia_var_1001)
  INTEGER(KIND = 4) :: ig_var_1026, ind0_var_1027, ind1_var_1028, inds_var_1029, indf_var_1030, i_lay_var_1031, i_laysolfr_var_1032(kidia_var_1000 : kfdia_var_1001), i_nlayers_var_1033, iplon_var_1034
  REAL(KIND = 8) :: z_tauray_var_1035
  laytrop_min_var_1020 = MINVAL(k_laytrop_var_1013(kidia_var_1000 : kfdia_var_1001))
  laytrop_max_var_1021 = MAXVAL(k_laytrop_var_1013(kidia_var_1000 : kfdia_var_1001))
  i_nlayers_var_1033 = klev_var_1002
  DO iplon_var_1034 = kidia_var_1000, kfdia_var_1001
    i_laysolfr_var_1032(iplon_var_1034) = i_nlayers_var_1033
  END DO
  DO i_lay_var_1031 = 1, laytrop_min_var_1020
    DO iplon_var_1034 = kidia_var_1000, kfdia_var_1001
      ind0_var_1027 = ((k_jp_var_1007(iplon_var_1034, i_lay_var_1031) - 1) * 5 + (k_jt_var_1008(iplon_var_1034, i_lay_var_1031) - 1)) * nspa_var_313(29) + 1
      ind1_var_1028 = (k_jp_var_1007(iplon_var_1034, i_lay_var_1031) * 5 + (k_jt1_var_1009(iplon_var_1034, i_lay_var_1031) - 1)) * nspa_var_313(29) + 1
      inds_var_1029 = k_indself_var_1016(iplon_var_1034, i_lay_var_1031)
      indf_var_1030 = k_indfor_var_1019(iplon_var_1034, i_lay_var_1031)
      z_tauray_var_1035 = p_colmol_var_1012(iplon_var_1034, i_lay_var_1031) * rayl_var_306
      DO ig_var_1026 = 1, 12
        p_taug_var_1023(iplon_var_1034, i_lay_var_1031, ig_var_1026) = p_colh2o_var_1010(iplon_var_1034, i_lay_var_1031) * ((p_fac00_var_1003(iplon_var_1034, i_lay_var_1031) * absa_var_308(ind0_var_1027, ig_var_1026) + p_fac10_var_1005(iplon_var_1034, i_lay_var_1031) * absa_var_308(ind0_var_1027 + 1, ig_var_1026) + p_fac01_var_1004(iplon_var_1034, i_lay_var_1031) * absa_var_308(ind1_var_1028, ig_var_1026) + p_fac11_var_1006(iplon_var_1034, i_lay_var_1031) * absa_var_308(ind1_var_1028 + 1, ig_var_1026)) + p_selffac_var_1014(iplon_var_1034, i_lay_var_1031) * (selfrefc_var_310(inds_var_1029, ig_var_1026) + p_selffrac_var_1015(iplon_var_1034, i_lay_var_1031) * (selfrefc_var_310(inds_var_1029 + 1, ig_var_1026) - selfrefc_var_310(inds_var_1029, ig_var_1026))) + p_forfac_var_1017(iplon_var_1034, i_lay_var_1031) * (forrefc_var_311(indf_var_1030, ig_var_1026) + p_forfrac_var_1018(iplon_var_1034, i_lay_var_1031) * (forrefc_var_311(indf_var_1030 + 1, ig_var_1026) - forrefc_var_311(indf_var_1030, ig_var_1026)))) + p_colco2_var_1011(iplon_var_1034, i_lay_var_1031) * absco2c(ig_var_1026)
        p_taur_var_1024(iplon_var_1034, i_lay_var_1031, ig_var_1026) = z_tauray_var_1035
      END DO
    END DO
  END DO
  DO i_lay_var_1031 = laytrop_min_var_1020 + 1, laytrop_max_var_1021
    DO iplon_var_1034 = kidia_var_1000, kfdia_var_1001
      IF (i_lay_var_1031 <= k_laytrop_var_1013(iplon_var_1034)) THEN
        ind0_var_1027 = ((k_jp_var_1007(iplon_var_1034, i_lay_var_1031) - 1) * 5 + (k_jt_var_1008(iplon_var_1034, i_lay_var_1031) - 1)) * nspa_var_313(29) + 1
        ind1_var_1028 = (k_jp_var_1007(iplon_var_1034, i_lay_var_1031) * 5 + (k_jt1_var_1009(iplon_var_1034, i_lay_var_1031) - 1)) * nspa_var_313(29) + 1
        inds_var_1029 = k_indself_var_1016(iplon_var_1034, i_lay_var_1031)
        indf_var_1030 = k_indfor_var_1019(iplon_var_1034, i_lay_var_1031)
        z_tauray_var_1035 = p_colmol_var_1012(iplon_var_1034, i_lay_var_1031) * rayl_var_306
        DO ig_var_1026 = 1, 12
          p_taug_var_1023(iplon_var_1034, i_lay_var_1031, ig_var_1026) = p_colh2o_var_1010(iplon_var_1034, i_lay_var_1031) * ((p_fac00_var_1003(iplon_var_1034, i_lay_var_1031) * absa_var_308(ind0_var_1027, ig_var_1026) + p_fac10_var_1005(iplon_var_1034, i_lay_var_1031) * absa_var_308(ind0_var_1027 + 1, ig_var_1026) + p_fac01_var_1004(iplon_var_1034, i_lay_var_1031) * absa_var_308(ind1_var_1028, ig_var_1026) + p_fac11_var_1006(iplon_var_1034, i_lay_var_1031) * absa_var_308(ind1_var_1028 + 1, ig_var_1026)) + p_selffac_var_1014(iplon_var_1034, i_lay_var_1031) * (selfrefc_var_310(inds_var_1029, ig_var_1026) + p_selffrac_var_1015(iplon_var_1034, i_lay_var_1031) * (selfrefc_var_310(inds_var_1029 + 1, ig_var_1026) - selfrefc_var_310(inds_var_1029, ig_var_1026))) + p_forfac_var_1017(iplon_var_1034, i_lay_var_1031) * (forrefc_var_311(indf_var_1030, ig_var_1026) + p_forfrac_var_1018(iplon_var_1034, i_lay_var_1031) * (forrefc_var_311(indf_var_1030 + 1, ig_var_1026) - forrefc_var_311(indf_var_1030, ig_var_1026)))) + p_colco2_var_1011(iplon_var_1034, i_lay_var_1031) * absco2c(ig_var_1026)
          p_taur_var_1024(iplon_var_1034, i_lay_var_1031, ig_var_1026) = z_tauray_var_1035
        END DO
      ELSE
        IF (k_jp_var_1007(iplon_var_1034, i_lay_var_1031 - 1) < layreffr_var_307 .AND. k_jp_var_1007(iplon_var_1034, i_lay_var_1031) >= layreffr_var_307) i_laysolfr_var_1032(iplon_var_1034) = i_lay_var_1031
        ind0_var_1027 = ((k_jp_var_1007(iplon_var_1034, i_lay_var_1031) - 13) * 5 + (k_jt_var_1008(iplon_var_1034, i_lay_var_1031) - 1)) * nspb_var_314(29) + 1
        ind1_var_1028 = ((k_jp_var_1007(iplon_var_1034, i_lay_var_1031) - 12) * 5 + (k_jt1_var_1009(iplon_var_1034, i_lay_var_1031) - 1)) * nspb_var_314(29) + 1
        z_tauray_var_1035 = p_colmol_var_1012(iplon_var_1034, i_lay_var_1031) * rayl_var_306
        DO ig_var_1026 = 1, 12
          p_taug_var_1023(iplon_var_1034, i_lay_var_1031, ig_var_1026) = p_colco2_var_1011(iplon_var_1034, i_lay_var_1031) * (p_fac00_var_1003(iplon_var_1034, i_lay_var_1031) * absb_var_309(ind0_var_1027, ig_var_1026) + p_fac10_var_1005(iplon_var_1034, i_lay_var_1031) * absb_var_309(ind0_var_1027 + 1, ig_var_1026) + p_fac01_var_1004(iplon_var_1034, i_lay_var_1031) * absb_var_309(ind1_var_1028, ig_var_1026) + p_fac11_var_1006(iplon_var_1034, i_lay_var_1031) * absb_var_309(ind1_var_1028 + 1, ig_var_1026)) + p_colh2o_var_1010(iplon_var_1034, i_lay_var_1031) * absh2oc(ig_var_1026)
          IF (i_lay_var_1031 == i_laysolfr_var_1032(iplon_var_1034)) p_sfluxzen_var_1022(iplon_var_1034, ig_var_1026) = sfluxrefc_var_312(ig_var_1026)
          p_taur_var_1024(iplon_var_1034, i_lay_var_1031, ig_var_1026) = z_tauray_var_1035
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1031 = laytrop_max_var_1021 + 1, i_nlayers_var_1033
    DO iplon_var_1034 = kidia_var_1000, kfdia_var_1001
      IF (k_jp_var_1007(iplon_var_1034, i_lay_var_1031 - 1) < layreffr_var_307 .AND. k_jp_var_1007(iplon_var_1034, i_lay_var_1031) >= layreffr_var_307) i_laysolfr_var_1032(iplon_var_1034) = i_lay_var_1031
      ind0_var_1027 = ((k_jp_var_1007(iplon_var_1034, i_lay_var_1031) - 13) * 5 + (k_jt_var_1008(iplon_var_1034, i_lay_var_1031) - 1)) * nspb_var_314(29) + 1
      ind1_var_1028 = ((k_jp_var_1007(iplon_var_1034, i_lay_var_1031) - 12) * 5 + (k_jt1_var_1009(iplon_var_1034, i_lay_var_1031) - 1)) * nspb_var_314(29) + 1
      z_tauray_var_1035 = p_colmol_var_1012(iplon_var_1034, i_lay_var_1031) * rayl_var_306
      DO ig_var_1026 = 1, 12
        p_taug_var_1023(iplon_var_1034, i_lay_var_1031, ig_var_1026) = p_colco2_var_1011(iplon_var_1034, i_lay_var_1031) * (p_fac00_var_1003(iplon_var_1034, i_lay_var_1031) * absb_var_309(ind0_var_1027, ig_var_1026) + p_fac10_var_1005(iplon_var_1034, i_lay_var_1031) * absb_var_309(ind0_var_1027 + 1, ig_var_1026) + p_fac01_var_1004(iplon_var_1034, i_lay_var_1031) * absb_var_309(ind1_var_1028, ig_var_1026) + p_fac11_var_1006(iplon_var_1034, i_lay_var_1031) * absb_var_309(ind1_var_1028 + 1, ig_var_1026)) + p_colh2o_var_1010(iplon_var_1034, i_lay_var_1031) * absh2oc(ig_var_1026)
        IF (i_lay_var_1031 == i_laysolfr_var_1032(iplon_var_1034)) p_sfluxzen_var_1022(iplon_var_1034, ig_var_1026) = sfluxrefc_var_312(ig_var_1026)
        p_taur_var_1024(iplon_var_1034, i_lay_var_1031, ig_var_1026) = z_tauray_var_1035
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol29
SUBROUTINE srtm_taumol17(kidia_var_1036, kfdia_var_1037, klev_var_1038, p_fac00_var_1039, p_fac01_var_1040, p_fac10_var_1041, p_fac11_var_1042, k_jp_var_1043, k_jt_var_1044, k_jt1_var_1045, p_oneminus_var_1046, p_colh2o_var_1047, p_colco2_var_1048, p_colmol_var_1049, k_laytrop_var_1050, p_selffac_var_1051, p_selffrac_var_1052, k_indself_var_1053, p_forfac_var_1054, p_forfrac_var_1055, k_indfor_var_1056, p_sfluxzen_var_1057, p_taug_var_1058, p_taur_var_1059, prmu0_var_1060)
  USE yoesrta17, ONLY: absa_var_228, absb_var_229, forrefc_var_231, layreffr_var_227, rayl_var_225, selfrefc_var_230, sfluxrefc_var_232, strrat_var_226
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1036, kfdia_var_1037
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1038
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1039(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1040(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1041(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1042(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1043(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1044(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1045(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1046(kidia_var_1036 : kfdia_var_1037)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1047(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1048(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1049(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1050(kidia_var_1036 : kfdia_var_1037)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1051(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1052(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1053(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1054(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1055(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1056(kidia_var_1036 : kfdia_var_1037, klev_var_1038)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1057(kidia_var_1036 : kfdia_var_1037, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1058(kidia_var_1036 : kfdia_var_1037, klev_var_1038, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1059(kidia_var_1036 : kfdia_var_1037, klev_var_1038, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1060(kidia_var_1036 : kfdia_var_1037)
  INTEGER(KIND = 4) :: ig_var_1061, ind0_var_1062, ind1_var_1063, inds_var_1064, indf_var_1065, js_var_1066, i_lay_var_1067, i_laysolfr_var_1068(kidia_var_1036 : kfdia_var_1037), i_nlayers_var_1069, iplon_var_1070
  INTEGER(KIND = 4) :: laytrop_min_var_1071, laytrop_max_var_1072
  REAL(KIND = 8) :: z_fs_var_1073, z_speccomb_var_1074, z_specmult_var_1075, z_specparm_var_1076, z_tauray_var_1077
  laytrop_min_var_1071 = MINVAL(k_laytrop_var_1050(kidia_var_1036 : kfdia_var_1037))
  laytrop_max_var_1072 = MAXVAL(k_laytrop_var_1050(kidia_var_1036 : kfdia_var_1037))
  i_nlayers_var_1069 = klev_var_1038
  DO iplon_var_1070 = kidia_var_1036, kfdia_var_1037
    i_laysolfr_var_1068(iplon_var_1070) = i_nlayers_var_1069
  END DO
  DO i_lay_var_1067 = 1, laytrop_min_var_1071
    DO iplon_var_1070 = kidia_var_1036, kfdia_var_1037
      z_speccomb_var_1074 = p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) + strrat_var_226 * p_colco2_var_1048(iplon_var_1070, i_lay_var_1067)
      z_specparm_var_1076 = p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) / z_speccomb_var_1074
      z_specparm_var_1076 = MIN(p_oneminus_var_1046(iplon_var_1070), z_specparm_var_1076)
      z_specmult_var_1075 = 8.0D0 * (z_specparm_var_1076)
      js_var_1066 = 1 + INT(z_specmult_var_1075)
      z_fs_var_1073 = z_specmult_var_1075 - AINT(z_specmult_var_1075)
      ind0_var_1062 = ((k_jp_var_1043(iplon_var_1070, i_lay_var_1067) - 1) * 5 + (k_jt_var_1044(iplon_var_1070, i_lay_var_1067) - 1)) * nspa_var_313(17) + js_var_1066
      ind1_var_1063 = (k_jp_var_1043(iplon_var_1070, i_lay_var_1067) * 5 + (k_jt1_var_1045(iplon_var_1070, i_lay_var_1067) - 1)) * nspa_var_313(17) + js_var_1066
      inds_var_1064 = k_indself_var_1053(iplon_var_1070, i_lay_var_1067)
      indf_var_1065 = k_indfor_var_1056(iplon_var_1070, i_lay_var_1067)
      z_tauray_var_1077 = p_colmol_var_1049(iplon_var_1070, i_lay_var_1067) * rayl_var_225
      DO ig_var_1061 = 1, 12
        p_taug_var_1058(iplon_var_1070, i_lay_var_1067, ig_var_1061) = z_speccomb_var_1074 * ((1.0D0 - z_fs_var_1073) * (absa_var_228(ind0_var_1062, ig_var_1061) * p_fac00_var_1039(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind0_var_1062 + 9, ig_var_1061) * p_fac10_var_1041(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind1_var_1063, ig_var_1061) * p_fac01_var_1040(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind1_var_1063 + 9, ig_var_1061) * p_fac11_var_1042(iplon_var_1070, i_lay_var_1067)) + z_fs_var_1073 * (absa_var_228(ind0_var_1062 + 1, ig_var_1061) * p_fac00_var_1039(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind0_var_1062 + 10, ig_var_1061) * p_fac10_var_1041(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind1_var_1063 + 1, ig_var_1061) * p_fac01_var_1040(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind1_var_1063 + 10, ig_var_1061) * p_fac11_var_1042(iplon_var_1070, i_lay_var_1067))) + p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) * (p_selffac_var_1051(iplon_var_1070, i_lay_var_1067) * (selfrefc_var_230(inds_var_1064, ig_var_1061) + p_selffrac_var_1052(iplon_var_1070, i_lay_var_1067) * (selfrefc_var_230(inds_var_1064 + 1, ig_var_1061) - selfrefc_var_230(inds_var_1064, ig_var_1061))) + p_forfac_var_1054(iplon_var_1070, i_lay_var_1067) * (forrefc_var_231(indf_var_1065, ig_var_1061) + p_forfrac_var_1055(iplon_var_1070, i_lay_var_1067) * (forrefc_var_231(indf_var_1065 + 1, ig_var_1061) - forrefc_var_231(indf_var_1065, ig_var_1061))))
        p_taur_var_1059(iplon_var_1070, i_lay_var_1067, ig_var_1061) = z_tauray_var_1077
      END DO
    END DO
  END DO
  DO i_lay_var_1067 = laytrop_min_var_1071 + 1, laytrop_max_var_1072
    DO iplon_var_1070 = kidia_var_1036, kfdia_var_1037
      IF (i_lay_var_1067 <= k_laytrop_var_1050(iplon_var_1070)) THEN
        z_speccomb_var_1074 = p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) + strrat_var_226 * p_colco2_var_1048(iplon_var_1070, i_lay_var_1067)
        z_specparm_var_1076 = p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) / z_speccomb_var_1074
        z_specparm_var_1076 = MIN(p_oneminus_var_1046(iplon_var_1070), z_specparm_var_1076)
        z_specmult_var_1075 = 8.0D0 * (z_specparm_var_1076)
        js_var_1066 = 1 + INT(z_specmult_var_1075)
        z_fs_var_1073 = z_specmult_var_1075 - AINT(z_specmult_var_1075)
        ind0_var_1062 = ((k_jp_var_1043(iplon_var_1070, i_lay_var_1067) - 1) * 5 + (k_jt_var_1044(iplon_var_1070, i_lay_var_1067) - 1)) * nspa_var_313(17) + js_var_1066
        ind1_var_1063 = (k_jp_var_1043(iplon_var_1070, i_lay_var_1067) * 5 + (k_jt1_var_1045(iplon_var_1070, i_lay_var_1067) - 1)) * nspa_var_313(17) + js_var_1066
        inds_var_1064 = k_indself_var_1053(iplon_var_1070, i_lay_var_1067)
        indf_var_1065 = k_indfor_var_1056(iplon_var_1070, i_lay_var_1067)
        z_tauray_var_1077 = p_colmol_var_1049(iplon_var_1070, i_lay_var_1067) * rayl_var_225
        DO ig_var_1061 = 1, 12
          p_taug_var_1058(iplon_var_1070, i_lay_var_1067, ig_var_1061) = z_speccomb_var_1074 * ((1.0D0 - z_fs_var_1073) * (absa_var_228(ind0_var_1062, ig_var_1061) * p_fac00_var_1039(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind0_var_1062 + 9, ig_var_1061) * p_fac10_var_1041(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind1_var_1063, ig_var_1061) * p_fac01_var_1040(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind1_var_1063 + 9, ig_var_1061) * p_fac11_var_1042(iplon_var_1070, i_lay_var_1067)) + z_fs_var_1073 * (absa_var_228(ind0_var_1062 + 1, ig_var_1061) * p_fac00_var_1039(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind0_var_1062 + 10, ig_var_1061) * p_fac10_var_1041(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind1_var_1063 + 1, ig_var_1061) * p_fac01_var_1040(iplon_var_1070, i_lay_var_1067) + absa_var_228(ind1_var_1063 + 10, ig_var_1061) * p_fac11_var_1042(iplon_var_1070, i_lay_var_1067))) + p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) * (p_selffac_var_1051(iplon_var_1070, i_lay_var_1067) * (selfrefc_var_230(inds_var_1064, ig_var_1061) + p_selffrac_var_1052(iplon_var_1070, i_lay_var_1067) * (selfrefc_var_230(inds_var_1064 + 1, ig_var_1061) - selfrefc_var_230(inds_var_1064, ig_var_1061))) + p_forfac_var_1054(iplon_var_1070, i_lay_var_1067) * (forrefc_var_231(indf_var_1065, ig_var_1061) + p_forfrac_var_1055(iplon_var_1070, i_lay_var_1067) * (forrefc_var_231(indf_var_1065 + 1, ig_var_1061) - forrefc_var_231(indf_var_1065, ig_var_1061))))
          p_taur_var_1059(iplon_var_1070, i_lay_var_1067, ig_var_1061) = z_tauray_var_1077
        END DO
      ELSE
        IF (k_jp_var_1043(iplon_var_1070, i_lay_var_1067 - 1) < layreffr_var_227 .AND. k_jp_var_1043(iplon_var_1070, i_lay_var_1067) >= layreffr_var_227) i_laysolfr_var_1068(iplon_var_1070) = i_lay_var_1067
        z_speccomb_var_1074 = p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) + strrat_var_226 * p_colco2_var_1048(iplon_var_1070, i_lay_var_1067)
        z_specparm_var_1076 = p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) / z_speccomb_var_1074
        z_specparm_var_1076 = MIN(p_oneminus_var_1046(iplon_var_1070), z_specparm_var_1076)
        z_specmult_var_1075 = 4.0D0 * (z_specparm_var_1076)
        js_var_1066 = 1 + INT(z_specmult_var_1075)
        z_fs_var_1073 = z_specmult_var_1075 - AINT(z_specmult_var_1075)
        ind0_var_1062 = ((k_jp_var_1043(iplon_var_1070, i_lay_var_1067) - 13) * 5 + (k_jt_var_1044(iplon_var_1070, i_lay_var_1067) - 1)) * nspb_var_314(17) + js_var_1066
        ind1_var_1063 = ((k_jp_var_1043(iplon_var_1070, i_lay_var_1067) - 12) * 5 + (k_jt1_var_1045(iplon_var_1070, i_lay_var_1067) - 1)) * nspb_var_314(17) + js_var_1066
        indf_var_1065 = k_indfor_var_1056(iplon_var_1070, i_lay_var_1067)
        z_tauray_var_1077 = p_colmol_var_1049(iplon_var_1070, i_lay_var_1067) * rayl_var_225
        DO ig_var_1061 = 1, 12
          p_taug_var_1058(iplon_var_1070, i_lay_var_1067, ig_var_1061) = z_speccomb_var_1074 * ((1.0D0 - z_fs_var_1073) * (absb_var_229(ind0_var_1062, ig_var_1061) * p_fac00_var_1039(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind0_var_1062 + 5, ig_var_1061) * p_fac10_var_1041(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind1_var_1063, ig_var_1061) * p_fac01_var_1040(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind1_var_1063 + 5, ig_var_1061) * p_fac11_var_1042(iplon_var_1070, i_lay_var_1067)) + z_fs_var_1073 * (absb_var_229(ind0_var_1062 + 1, ig_var_1061) * p_fac00_var_1039(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind0_var_1062 + 6, ig_var_1061) * p_fac10_var_1041(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind1_var_1063 + 1, ig_var_1061) * p_fac01_var_1040(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind1_var_1063 + 6, ig_var_1061) * p_fac11_var_1042(iplon_var_1070, i_lay_var_1067))) + p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) * p_forfac_var_1054(iplon_var_1070, i_lay_var_1067) * (forrefc_var_231(indf_var_1065, ig_var_1061) + p_forfrac_var_1055(iplon_var_1070, i_lay_var_1067) * (forrefc_var_231(indf_var_1065 + 1, ig_var_1061) - forrefc_var_231(indf_var_1065, ig_var_1061)))
          IF (i_lay_var_1067 == i_laysolfr_var_1068(iplon_var_1070)) p_sfluxzen_var_1057(iplon_var_1070, ig_var_1061) = sfluxrefc_var_232(ig_var_1061, js_var_1066) + z_fs_var_1073 * (sfluxrefc_var_232(ig_var_1061, js_var_1066 + 1) - sfluxrefc_var_232(ig_var_1061, js_var_1066))
          p_taur_var_1059(iplon_var_1070, i_lay_var_1067, ig_var_1061) = z_tauray_var_1077
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1067 = laytrop_max_var_1072 + 1, i_nlayers_var_1069
    DO iplon_var_1070 = kidia_var_1036, kfdia_var_1037
      IF (k_jp_var_1043(iplon_var_1070, i_lay_var_1067 - 1) < layreffr_var_227 .AND. k_jp_var_1043(iplon_var_1070, i_lay_var_1067) >= layreffr_var_227) i_laysolfr_var_1068(iplon_var_1070) = i_lay_var_1067
      z_speccomb_var_1074 = p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) + strrat_var_226 * p_colco2_var_1048(iplon_var_1070, i_lay_var_1067)
      z_specparm_var_1076 = p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) / z_speccomb_var_1074
      z_specparm_var_1076 = MIN(p_oneminus_var_1046(iplon_var_1070), z_specparm_var_1076)
      z_specmult_var_1075 = 4.0D0 * (z_specparm_var_1076)
      js_var_1066 = 1 + INT(z_specmult_var_1075)
      z_fs_var_1073 = z_specmult_var_1075 - AINT(z_specmult_var_1075)
      ind0_var_1062 = ((k_jp_var_1043(iplon_var_1070, i_lay_var_1067) - 13) * 5 + (k_jt_var_1044(iplon_var_1070, i_lay_var_1067) - 1)) * nspb_var_314(17) + js_var_1066
      ind1_var_1063 = ((k_jp_var_1043(iplon_var_1070, i_lay_var_1067) - 12) * 5 + (k_jt1_var_1045(iplon_var_1070, i_lay_var_1067) - 1)) * nspb_var_314(17) + js_var_1066
      indf_var_1065 = k_indfor_var_1056(iplon_var_1070, i_lay_var_1067)
      z_tauray_var_1077 = p_colmol_var_1049(iplon_var_1070, i_lay_var_1067) * rayl_var_225
      DO ig_var_1061 = 1, 12
        p_taug_var_1058(iplon_var_1070, i_lay_var_1067, ig_var_1061) = z_speccomb_var_1074 * ((1.0D0 - z_fs_var_1073) * (absb_var_229(ind0_var_1062, ig_var_1061) * p_fac00_var_1039(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind0_var_1062 + 5, ig_var_1061) * p_fac10_var_1041(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind1_var_1063, ig_var_1061) * p_fac01_var_1040(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind1_var_1063 + 5, ig_var_1061) * p_fac11_var_1042(iplon_var_1070, i_lay_var_1067)) + z_fs_var_1073 * (absb_var_229(ind0_var_1062 + 1, ig_var_1061) * p_fac00_var_1039(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind0_var_1062 + 6, ig_var_1061) * p_fac10_var_1041(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind1_var_1063 + 1, ig_var_1061) * p_fac01_var_1040(iplon_var_1070, i_lay_var_1067) + absb_var_229(ind1_var_1063 + 6, ig_var_1061) * p_fac11_var_1042(iplon_var_1070, i_lay_var_1067))) + p_colh2o_var_1047(iplon_var_1070, i_lay_var_1067) * p_forfac_var_1054(iplon_var_1070, i_lay_var_1067) * (forrefc_var_231(indf_var_1065, ig_var_1061) + p_forfrac_var_1055(iplon_var_1070, i_lay_var_1067) * (forrefc_var_231(indf_var_1065 + 1, ig_var_1061) - forrefc_var_231(indf_var_1065, ig_var_1061)))
        IF (i_lay_var_1067 == i_laysolfr_var_1068(iplon_var_1070)) p_sfluxzen_var_1057(iplon_var_1070, ig_var_1061) = sfluxrefc_var_232(ig_var_1061, js_var_1066) + z_fs_var_1073 * (sfluxrefc_var_232(ig_var_1061, js_var_1066 + 1) - sfluxrefc_var_232(ig_var_1061, js_var_1066))
        p_taur_var_1059(iplon_var_1070, i_lay_var_1067, ig_var_1061) = z_tauray_var_1077
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol17
SUBROUTINE rrtm_setcoef_140gp(kidia_var_1078, kfdia_var_1079, klev_var_1080, p_coldry, p_wbroad, p_wkl, p_fac00_var_1081, p_fac01_var_1082, p_fac10_var_1083, p_fac11_var_1084, p_forfac_var_1085, p_forfrac_var_1086, k_indfor_var_1103, k_jp_var_1087, k_jt_var_1088, k_jt1_var_1089, p_colh2o_var_1090, p_colco2_var_1091, p_colo3_var_1092, p_coln2o, p_colch4_var_1093, p_colo2_var_1094, p_co2mult_var_1095, p_colbrd, k_laytrop_var_1096, k_layswtch_var_1097, k_laylow_var_1098, pavel_var_1099, p_tavel, p_selffac_var_1100, p_selffrac_var_1101, k_indself_var_1102, k_indminor, p_scaleminor, p_scaleminorn2, p_minorfrac, prat_h2oco2_var_1104, prat_h2oco2_1_var_1105, prat_h2oo3_var_1106, prat_h2oo3_1_var_1107, prat_h2on2o_var_1108, prat_h2on2o_1_var_1109, prat_h2och4_var_1110, prat_h2och4_1_var_1111, prat_n2oco2_var_1112, prat_n2oco2_1_var_1113, prat_o3co2_var_1114, prat_o3co2_1_var_1115)
  USE yoerrtrf, ONLY: chi_mls, preflog_var_214, tref_var_215
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1078
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1079
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1080
  REAL(KIND = 8), INTENT(IN) :: p_coldry(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(IN) :: p_wbroad(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_colbrd(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(IN) :: p_wkl(kidia_var_1078 : kfdia_var_1079, 35, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_fac00_var_1081(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_fac01_var_1082(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_fac10_var_1083(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_fac11_var_1084(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_forfac_var_1085(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_forfrac_var_1086(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jp_var_1087(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jt_var_1088(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jt1_var_1089(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_colh2o_var_1090(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_colco2_var_1091(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_colo3_var_1092(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_coln2o(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_colch4_var_1093(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_colo2_var_1094(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_co2mult_var_1095(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  INTEGER(KIND = 4), INTENT(OUT) :: k_laytrop_var_1096(kidia_var_1078 : kfdia_var_1079)
  INTEGER(KIND = 4), INTENT(OUT) :: k_layswtch_var_1097(kidia_var_1078 : kfdia_var_1079)
  INTEGER(KIND = 4), INTENT(OUT) :: k_laylow_var_1098(kidia_var_1078 : kfdia_var_1079)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1099(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(IN) :: p_tavel(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_selffac_var_1100(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_selffrac_var_1101(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indself_var_1102(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indfor_var_1103(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indminor(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_scaleminor(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_scaleminorn2(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: p_minorfrac(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  REAL(KIND = 8), INTENT(OUT) :: prat_h2oco2_var_1104(kidia_var_1078 : kfdia_var_1079, klev_var_1080), prat_h2oco2_1_var_1105(kidia_var_1078 : kfdia_var_1079, klev_var_1080), prat_h2oo3_var_1106(kidia_var_1078 : kfdia_var_1079, klev_var_1080), prat_h2oo3_1_var_1107(kidia_var_1078 : kfdia_var_1079, klev_var_1080), prat_h2on2o_var_1108(kidia_var_1078 : kfdia_var_1079, klev_var_1080), prat_h2on2o_1_var_1109(kidia_var_1078 : kfdia_var_1079, klev_var_1080), prat_h2och4_var_1110(kidia_var_1078 : kfdia_var_1079, klev_var_1080), prat_h2och4_1_var_1111(kidia_var_1078 : kfdia_var_1079, klev_var_1080), prat_n2oco2_var_1112(kidia_var_1078 : kfdia_var_1079, klev_var_1080), prat_n2oco2_1_var_1113(kidia_var_1078 : kfdia_var_1079, klev_var_1080), prat_o3co2_var_1114(kidia_var_1078 : kfdia_var_1079, klev_var_1080), prat_o3co2_1_var_1115(kidia_var_1078 : kfdia_var_1079, klev_var_1080)
  INTEGER(KIND = 4) :: jp1_var_1116, jlay
  INTEGER(KIND = 4) :: jlon_var_1117
  REAL(KIND = 8) :: z_co2reg_var_1118, z_compfp_var_1119, z_factor_var_1120, z_fp_var_1121, z_ft_var_1122, z_ft1_var_1123, z_plog_var_1124, z_scalefac_var_1125, z_stpfac_var_1126, z_water_var_1127
  DO jlon_var_1117 = kidia_var_1078, kfdia_var_1079
    z_stpfac_var_1126 = 0.29220138203356366D0
    k_laytrop_var_1096(jlon_var_1117) = 0
    k_layswtch_var_1097(jlon_var_1117) = 0
    k_laylow_var_1098(jlon_var_1117) = 0
    DO jlay = 1, klev_var_1080
      z_plog_var_1124 = LOG(pavel_var_1099(jlon_var_1117, jlay))
      k_jp_var_1087(jlon_var_1117, jlay) = INT(36.0D0 - 5 * (z_plog_var_1124 + 0.04D0))
      IF (k_jp_var_1087(jlon_var_1117, jlay) < 1) THEN
        k_jp_var_1087(jlon_var_1117, jlay) = 1
      ELSE IF (k_jp_var_1087(jlon_var_1117, jlay) > 58) THEN
        k_jp_var_1087(jlon_var_1117, jlay) = 58
      END IF
      jp1_var_1116 = k_jp_var_1087(jlon_var_1117, jlay) + 1
      z_fp_var_1121 = 5.0D0 * (preflog_var_214(k_jp_var_1087(jlon_var_1117, jlay)) - z_plog_var_1124)
      z_fp_var_1121 = MAX(- 1.0D0, MIN(1.0D0, z_fp_var_1121))
      k_jt_var_1088(jlon_var_1117, jlay) = INT(3.0D0 + (p_tavel(jlon_var_1117, jlay) - tref_var_215(k_jp_var_1087(jlon_var_1117, jlay))) / 15.0D0)
      IF (k_jt_var_1088(jlon_var_1117, jlay) < 1) THEN
        k_jt_var_1088(jlon_var_1117, jlay) = 1
      ELSE IF (k_jt_var_1088(jlon_var_1117, jlay) > 4) THEN
        k_jt_var_1088(jlon_var_1117, jlay) = 4
      END IF
      z_ft_var_1122 = ((p_tavel(jlon_var_1117, jlay) - tref_var_215(k_jp_var_1087(jlon_var_1117, jlay))) / 15.0D0) - REAL(k_jt_var_1088(jlon_var_1117, jlay) - 3)
      k_jt1_var_1089(jlon_var_1117, jlay) = INT(3.0D0 + (p_tavel(jlon_var_1117, jlay) - tref_var_215(jp1_var_1116)) / 15.0D0)
      IF (k_jt1_var_1089(jlon_var_1117, jlay) < 1) THEN
        k_jt1_var_1089(jlon_var_1117, jlay) = 1
      ELSE IF (k_jt1_var_1089(jlon_var_1117, jlay) > 4) THEN
        k_jt1_var_1089(jlon_var_1117, jlay) = 4
      END IF
      z_ft1_var_1123 = ((p_tavel(jlon_var_1117, jlay) - tref_var_215(jp1_var_1116)) / 15.0D0) - REAL(k_jt1_var_1089(jlon_var_1117, jlay) - 3)
      z_water_var_1127 = p_wkl(jlon_var_1117, 1, jlay) / p_coldry(jlon_var_1117, jlay)
      z_scalefac_var_1125 = pavel_var_1099(jlon_var_1117, jlay) * z_stpfac_var_1126 / p_tavel(jlon_var_1117, jlay)
      IF (z_plog_var_1124 > 4.56D0) THEN
        k_laytrop_var_1096(jlon_var_1117) = k_laytrop_var_1096(jlon_var_1117) + 1
        p_forfac_var_1085(jlon_var_1117, jlay) = z_scalefac_var_1125 / (1.0D0 + z_water_var_1127)
        z_factor_var_1120 = (332.0D0 - p_tavel(jlon_var_1117, jlay)) / 36.0D0
        k_indfor_var_1103(jlon_var_1117, jlay) = MIN(2, MAX(1, INT(z_factor_var_1120)))
        p_forfrac_var_1086(jlon_var_1117, jlay) = z_factor_var_1120 - REAL(k_indfor_var_1103(jlon_var_1117, jlay))
        p_selffac_var_1100(jlon_var_1117, jlay) = z_water_var_1127 * p_forfac_var_1085(jlon_var_1117, jlay)
        z_factor_var_1120 = (p_tavel(jlon_var_1117, jlay) - 188.0D0) / 7.2D0
        k_indself_var_1102(jlon_var_1117, jlay) = MIN(9, MAX(1, INT(z_factor_var_1120) - 7))
        p_selffrac_var_1101(jlon_var_1117, jlay) = z_factor_var_1120 - REAL(k_indself_var_1102(jlon_var_1117, jlay) + 7)
        p_scaleminor(jlon_var_1117, jlay) = pavel_var_1099(jlon_var_1117, jlay) / p_tavel(jlon_var_1117, jlay)
        p_scaleminorn2(jlon_var_1117, jlay) = (pavel_var_1099(jlon_var_1117, jlay) / p_tavel(jlon_var_1117, jlay)) * (p_wbroad(jlon_var_1117, jlay) / (p_coldry(jlon_var_1117, jlay) + p_wkl(jlon_var_1117, 1, jlay)))
        z_factor_var_1120 = (p_tavel(jlon_var_1117, jlay) - 180.8D0) / 7.2D0
        k_indminor(jlon_var_1117, jlay) = MIN(18, MAX(1, INT(z_factor_var_1120)))
        p_minorfrac(jlon_var_1117, jlay) = z_factor_var_1120 - REAL(k_indminor(jlon_var_1117, jlay))
        prat_h2oco2_var_1104(jlon_var_1117, jlay) = chi_mls(1, k_jp_var_1087(jlon_var_1117, jlay)) / chi_mls(2, k_jp_var_1087(jlon_var_1117, jlay))
        prat_h2oco2_1_var_1105(jlon_var_1117, jlay) = chi_mls(1, k_jp_var_1087(jlon_var_1117, jlay) + 1) / chi_mls(2, k_jp_var_1087(jlon_var_1117, jlay) + 1)
        prat_h2oo3_var_1106(jlon_var_1117, jlay) = chi_mls(1, k_jp_var_1087(jlon_var_1117, jlay)) / chi_mls(3, k_jp_var_1087(jlon_var_1117, jlay))
        prat_h2oo3_1_var_1107(jlon_var_1117, jlay) = chi_mls(1, k_jp_var_1087(jlon_var_1117, jlay) + 1) / chi_mls(3, k_jp_var_1087(jlon_var_1117, jlay) + 1)
        prat_h2on2o_var_1108(jlon_var_1117, jlay) = chi_mls(1, k_jp_var_1087(jlon_var_1117, jlay)) / chi_mls(4, k_jp_var_1087(jlon_var_1117, jlay))
        prat_h2on2o_1_var_1109(jlon_var_1117, jlay) = chi_mls(1, k_jp_var_1087(jlon_var_1117, jlay) + 1) / chi_mls(4, k_jp_var_1087(jlon_var_1117, jlay) + 1)
        prat_h2och4_var_1110(jlon_var_1117, jlay) = chi_mls(1, k_jp_var_1087(jlon_var_1117, jlay)) / chi_mls(6, k_jp_var_1087(jlon_var_1117, jlay))
        prat_h2och4_1_var_1111(jlon_var_1117, jlay) = chi_mls(1, k_jp_var_1087(jlon_var_1117, jlay) + 1) / chi_mls(6, k_jp_var_1087(jlon_var_1117, jlay) + 1)
        prat_n2oco2_var_1112(jlon_var_1117, jlay) = chi_mls(4, k_jp_var_1087(jlon_var_1117, jlay)) / chi_mls(2, k_jp_var_1087(jlon_var_1117, jlay))
        prat_n2oco2_1_var_1113(jlon_var_1117, jlay) = chi_mls(4, k_jp_var_1087(jlon_var_1117, jlay) + 1) / chi_mls(2, k_jp_var_1087(jlon_var_1117, jlay) + 1)
        p_colh2o_var_1090(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 1, jlay)
        p_colco2_var_1091(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 2, jlay)
        p_colo3_var_1092(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 3, jlay)
        p_coln2o(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 4, jlay)
        p_colch4_var_1093(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 6, jlay)
        p_colo2_var_1094(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 7, jlay)
        p_colbrd(jlon_var_1117, jlay) = 1D-20 * p_wbroad(jlon_var_1117, jlay)
        IF (p_colco2_var_1091(jlon_var_1117, jlay) == 0.0D0) p_colco2_var_1091(jlon_var_1117, jlay) = 1D-32 * p_coldry(jlon_var_1117, jlay)
        IF (p_coln2o(jlon_var_1117, jlay) == 0.0D0) p_coln2o(jlon_var_1117, jlay) = 1D-32 * p_coldry(jlon_var_1117, jlay)
        IF (p_colch4_var_1093(jlon_var_1117, jlay) == 0.0D0) p_colch4_var_1093(jlon_var_1117, jlay) = 1D-32 * p_coldry(jlon_var_1117, jlay)
        z_co2reg_var_1118 = 3.55D-24 * p_coldry(jlon_var_1117, jlay)
        p_co2mult_var_1095(jlon_var_1117, jlay) = (p_colco2_var_1091(jlon_var_1117, jlay) - z_co2reg_var_1118) * 272.63D0 * EXP(- 1919.4D0 / p_tavel(jlon_var_1117, jlay)) / (0.00087604D0 * p_tavel(jlon_var_1117, jlay))
      ELSE
        p_forfac_var_1085(jlon_var_1117, jlay) = z_scalefac_var_1125 / (1.0D0 + z_water_var_1127)
        z_factor_var_1120 = (p_tavel(jlon_var_1117, jlay) - 188.0D0) / 36.0D0
        k_indfor_var_1103(jlon_var_1117, jlay) = 3
        p_forfrac_var_1086(jlon_var_1117, jlay) = z_factor_var_1120 - 1.0D0
        p_selffac_var_1100(jlon_var_1117, jlay) = z_water_var_1127 * p_forfac_var_1085(jlon_var_1117, jlay)
        p_scaleminor(jlon_var_1117, jlay) = pavel_var_1099(jlon_var_1117, jlay) / p_tavel(jlon_var_1117, jlay)
        p_scaleminorn2(jlon_var_1117, jlay) = (pavel_var_1099(jlon_var_1117, jlay) / p_tavel(jlon_var_1117, jlay)) * (p_wbroad(jlon_var_1117, jlay) / (p_coldry(jlon_var_1117, jlay) + p_wkl(jlon_var_1117, 1, jlay)))
        z_factor_var_1120 = (p_tavel(jlon_var_1117, jlay) - 180.8D0) / 7.2D0
        k_indminor(jlon_var_1117, jlay) = MIN(18, MAX(1, INT(z_factor_var_1120)))
        p_minorfrac(jlon_var_1117, jlay) = z_factor_var_1120 - REAL(k_indminor(jlon_var_1117, jlay))
        prat_h2oco2_var_1104(jlon_var_1117, jlay) = chi_mls(1, k_jp_var_1087(jlon_var_1117, jlay)) / chi_mls(2, k_jp_var_1087(jlon_var_1117, jlay))
        prat_h2oco2_1_var_1105(jlon_var_1117, jlay) = chi_mls(1, k_jp_var_1087(jlon_var_1117, jlay) + 1) / chi_mls(2, k_jp_var_1087(jlon_var_1117, jlay) + 1)
        prat_o3co2_var_1114(jlon_var_1117, jlay) = chi_mls(3, k_jp_var_1087(jlon_var_1117, jlay)) / chi_mls(2, k_jp_var_1087(jlon_var_1117, jlay))
        prat_o3co2_1_var_1115(jlon_var_1117, jlay) = chi_mls(3, k_jp_var_1087(jlon_var_1117, jlay) + 1) / chi_mls(2, k_jp_var_1087(jlon_var_1117, jlay) + 1)
        p_colh2o_var_1090(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 1, jlay)
        p_colco2_var_1091(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 2, jlay)
        p_colo3_var_1092(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 3, jlay)
        p_coln2o(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 4, jlay)
        p_colch4_var_1093(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 6, jlay)
        p_colo2_var_1094(jlon_var_1117, jlay) = 1D-20 * p_wkl(jlon_var_1117, 7, jlay)
        p_colbrd(jlon_var_1117, jlay) = 1D-20 * p_wbroad(jlon_var_1117, jlay)
        IF (p_colco2_var_1091(jlon_var_1117, jlay) == 0.0D0) p_colco2_var_1091(jlon_var_1117, jlay) = 1D-32 * p_coldry(jlon_var_1117, jlay)
        IF (p_coln2o(jlon_var_1117, jlay) == 0.0D0) p_coln2o(jlon_var_1117, jlay) = 1D-32 * p_coldry(jlon_var_1117, jlay)
        IF (p_colch4_var_1093(jlon_var_1117, jlay) == 0.0D0) p_colch4_var_1093(jlon_var_1117, jlay) = 1D-32 * p_coldry(jlon_var_1117, jlay)
        z_co2reg_var_1118 = 3.55D-24 * p_coldry(jlon_var_1117, jlay)
        p_co2mult_var_1095(jlon_var_1117, jlay) = (p_colco2_var_1091(jlon_var_1117, jlay) - z_co2reg_var_1118) * 272.63D0 * EXP(- 1919.4D0 / p_tavel(jlon_var_1117, jlay)) / (0.00087604D0 * p_tavel(jlon_var_1117, jlay))
      END IF
      z_compfp_var_1119 = 1.0D0 - z_fp_var_1121
      p_fac10_var_1083(jlon_var_1117, jlay) = z_compfp_var_1119 * z_ft_var_1122
      p_fac00_var_1081(jlon_var_1117, jlay) = z_compfp_var_1119 * (1.0D0 - z_ft_var_1122)
      p_fac11_var_1084(jlon_var_1117, jlay) = z_fp_var_1121 * z_ft1_var_1123
      p_fac01_var_1082(jlon_var_1117, jlay) = z_fp_var_1121 * (1.0D0 - z_ft1_var_1123)
      p_selffac_var_1100(jlon_var_1117, jlay) = p_colh2o_var_1090(jlon_var_1117, jlay) * p_selffac_var_1100(jlon_var_1117, jlay)
      p_forfac_var_1085(jlon_var_1117, jlay) = p_colh2o_var_1090(jlon_var_1117, jlay) * p_forfac_var_1085(jlon_var_1117, jlay)
    END DO
    IF (k_laylow_var_1098(jlon_var_1117) == 0) k_laylow_var_1098(jlon_var_1117) = 1
  END DO
END SUBROUTINE rrtm_setcoef_140gp
SUBROUTINE srtm_taumol16(kidia_var_1128, kfdia_var_1129, klev_var_1130, p_fac00_var_1131, p_fac01_var_1132, p_fac10_var_1133, p_fac11_var_1134, k_jp_var_1135, k_jt_var_1136, k_jt1_var_1137, p_oneminus_var_1138, p_colh2o_var_1139, p_colch4_var_1140, p_colmol_var_1141, k_laytrop_var_1142, p_selffac_var_1143, p_selffrac_var_1144, k_indself_var_1145, p_forfac_var_1146, p_forfrac_var_1147, k_indfor_var_1148, p_sfluxzen_var_1149, p_taug_var_1150, p_taur_var_1151, prmu0_var_1152)
  USE yoesrta16, ONLY: absa_var_220, absb_var_221, forrefc_var_223, layreffr_var_219, rayl_var_218, selfrefc_var_222, sfluxrefc_var_224, strrat1
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1128, kfdia_var_1129
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1130
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1131(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1132(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1133(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1134(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1135(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1136(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1137(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1138(kidia_var_1128 : kfdia_var_1129)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1139(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_1140(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1141(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1142(kidia_var_1128 : kfdia_var_1129)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1143(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1144(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1145(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1146(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1147(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1148(kidia_var_1128 : kfdia_var_1129, klev_var_1130)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1149(kidia_var_1128 : kfdia_var_1129, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1150(kidia_var_1128 : kfdia_var_1129, klev_var_1130, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1151(kidia_var_1128 : kfdia_var_1129, klev_var_1130, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1152(kidia_var_1128 : kfdia_var_1129)
  INTEGER(KIND = 4) :: ig_var_1153, ind0_var_1154, ind1_var_1155, inds_var_1156, indf_var_1157, js_var_1158, i_lay_var_1159, i_laysolfr_var_1160(kidia_var_1128 : kfdia_var_1129), i_nlayers_var_1161, iplon_var_1162
  INTEGER(KIND = 4) :: laytrop_min_var_1163, laytrop_max_var_1164
  REAL(KIND = 8) :: z_fs_var_1165, z_speccomb_var_1166, z_specmult_var_1167, z_specparm_var_1168, z_tauray_var_1169
  laytrop_min_var_1163 = MINVAL(k_laytrop_var_1142(kidia_var_1128 : kfdia_var_1129))
  laytrop_max_var_1164 = MAXVAL(k_laytrop_var_1142(kidia_var_1128 : kfdia_var_1129))
  i_nlayers_var_1161 = klev_var_1130
  DO iplon_var_1162 = kidia_var_1128, kfdia_var_1129
    i_laysolfr_var_1160(iplon_var_1162) = i_nlayers_var_1161
  END DO
  DO i_lay_var_1159 = 1, laytrop_min_var_1163
    DO iplon_var_1162 = kidia_var_1128, kfdia_var_1129
      z_speccomb_var_1166 = p_colh2o_var_1139(iplon_var_1162, i_lay_var_1159) + strrat1 * p_colch4_var_1140(iplon_var_1162, i_lay_var_1159)
      z_specparm_var_1168 = p_colh2o_var_1139(iplon_var_1162, i_lay_var_1159) / z_speccomb_var_1166
      z_specparm_var_1168 = MIN(p_oneminus_var_1138(iplon_var_1162), z_specparm_var_1168)
      z_specmult_var_1167 = 8.0D0 * (z_specparm_var_1168)
      js_var_1158 = 1 + INT(z_specmult_var_1167)
      z_fs_var_1165 = z_specmult_var_1167 - AINT(z_specmult_var_1167)
      ind0_var_1154 = ((k_jp_var_1135(iplon_var_1162, i_lay_var_1159) - 1) * 5 + (k_jt_var_1136(iplon_var_1162, i_lay_var_1159) - 1)) * nspa_var_313(16) + js_var_1158
      ind1_var_1155 = (k_jp_var_1135(iplon_var_1162, i_lay_var_1159) * 5 + (k_jt1_var_1137(iplon_var_1162, i_lay_var_1159) - 1)) * nspa_var_313(16) + js_var_1158
      inds_var_1156 = k_indself_var_1145(iplon_var_1162, i_lay_var_1159)
      indf_var_1157 = k_indfor_var_1148(iplon_var_1162, i_lay_var_1159)
      z_tauray_var_1169 = p_colmol_var_1141(iplon_var_1162, i_lay_var_1159) * rayl_var_218
      DO ig_var_1153 = 1, 6
        p_taug_var_1150(iplon_var_1162, i_lay_var_1159, ig_var_1153) = z_speccomb_var_1166 * ((1.0D0 - z_fs_var_1165) * (absa_var_220(ind0_var_1154, ig_var_1153) * p_fac00_var_1131(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind0_var_1154 + 9, ig_var_1153) * p_fac10_var_1133(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind1_var_1155, ig_var_1153) * p_fac01_var_1132(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind1_var_1155 + 9, ig_var_1153) * p_fac11_var_1134(iplon_var_1162, i_lay_var_1159)) + z_fs_var_1165 * (absa_var_220(ind0_var_1154 + 1, ig_var_1153) * p_fac00_var_1131(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind0_var_1154 + 10, ig_var_1153) * p_fac10_var_1133(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind1_var_1155 + 1, ig_var_1153) * p_fac01_var_1132(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind1_var_1155 + 10, ig_var_1153) * p_fac11_var_1134(iplon_var_1162, i_lay_var_1159))) + p_colh2o_var_1139(iplon_var_1162, i_lay_var_1159) * (p_selffac_var_1143(iplon_var_1162, i_lay_var_1159) * (selfrefc_var_222(inds_var_1156, ig_var_1153) + p_selffrac_var_1144(iplon_var_1162, i_lay_var_1159) * (selfrefc_var_222(inds_var_1156 + 1, ig_var_1153) - selfrefc_var_222(inds_var_1156, ig_var_1153))) + p_forfac_var_1146(iplon_var_1162, i_lay_var_1159) * (forrefc_var_223(indf_var_1157, ig_var_1153) + p_forfrac_var_1147(iplon_var_1162, i_lay_var_1159) * (forrefc_var_223(indf_var_1157 + 1, ig_var_1153) - forrefc_var_223(indf_var_1157, ig_var_1153))))
        p_taur_var_1151(iplon_var_1162, i_lay_var_1159, ig_var_1153) = z_tauray_var_1169
      END DO
    END DO
  END DO
  DO i_lay_var_1159 = laytrop_min_var_1163 + 1, laytrop_max_var_1164
    DO iplon_var_1162 = kidia_var_1128, kfdia_var_1129
      IF (i_lay_var_1159 <= k_laytrop_var_1142(iplon_var_1162)) THEN
        z_speccomb_var_1166 = p_colh2o_var_1139(iplon_var_1162, i_lay_var_1159) + strrat1 * p_colch4_var_1140(iplon_var_1162, i_lay_var_1159)
        z_specparm_var_1168 = p_colh2o_var_1139(iplon_var_1162, i_lay_var_1159) / z_speccomb_var_1166
        z_specparm_var_1168 = MIN(p_oneminus_var_1138(iplon_var_1162), z_specparm_var_1168)
        z_specmult_var_1167 = 8.0D0 * (z_specparm_var_1168)
        js_var_1158 = 1 + INT(z_specmult_var_1167)
        z_fs_var_1165 = z_specmult_var_1167 - AINT(z_specmult_var_1167)
        ind0_var_1154 = ((k_jp_var_1135(iplon_var_1162, i_lay_var_1159) - 1) * 5 + (k_jt_var_1136(iplon_var_1162, i_lay_var_1159) - 1)) * nspa_var_313(16) + js_var_1158
        ind1_var_1155 = (k_jp_var_1135(iplon_var_1162, i_lay_var_1159) * 5 + (k_jt1_var_1137(iplon_var_1162, i_lay_var_1159) - 1)) * nspa_var_313(16) + js_var_1158
        inds_var_1156 = k_indself_var_1145(iplon_var_1162, i_lay_var_1159)
        indf_var_1157 = k_indfor_var_1148(iplon_var_1162, i_lay_var_1159)
        z_tauray_var_1169 = p_colmol_var_1141(iplon_var_1162, i_lay_var_1159) * rayl_var_218
        DO ig_var_1153 = 1, 6
          p_taug_var_1150(iplon_var_1162, i_lay_var_1159, ig_var_1153) = z_speccomb_var_1166 * ((1.0D0 - z_fs_var_1165) * (absa_var_220(ind0_var_1154, ig_var_1153) * p_fac00_var_1131(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind0_var_1154 + 9, ig_var_1153) * p_fac10_var_1133(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind1_var_1155, ig_var_1153) * p_fac01_var_1132(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind1_var_1155 + 9, ig_var_1153) * p_fac11_var_1134(iplon_var_1162, i_lay_var_1159)) + z_fs_var_1165 * (absa_var_220(ind0_var_1154 + 1, ig_var_1153) * p_fac00_var_1131(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind0_var_1154 + 10, ig_var_1153) * p_fac10_var_1133(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind1_var_1155 + 1, ig_var_1153) * p_fac01_var_1132(iplon_var_1162, i_lay_var_1159) + absa_var_220(ind1_var_1155 + 10, ig_var_1153) * p_fac11_var_1134(iplon_var_1162, i_lay_var_1159))) + p_colh2o_var_1139(iplon_var_1162, i_lay_var_1159) * (p_selffac_var_1143(iplon_var_1162, i_lay_var_1159) * (selfrefc_var_222(inds_var_1156, ig_var_1153) + p_selffrac_var_1144(iplon_var_1162, i_lay_var_1159) * (selfrefc_var_222(inds_var_1156 + 1, ig_var_1153) - selfrefc_var_222(inds_var_1156, ig_var_1153))) + p_forfac_var_1146(iplon_var_1162, i_lay_var_1159) * (forrefc_var_223(indf_var_1157, ig_var_1153) + p_forfrac_var_1147(iplon_var_1162, i_lay_var_1159) * (forrefc_var_223(indf_var_1157 + 1, ig_var_1153) - forrefc_var_223(indf_var_1157, ig_var_1153))))
          p_taur_var_1151(iplon_var_1162, i_lay_var_1159, ig_var_1153) = z_tauray_var_1169
        END DO
      ELSE
        IF (k_jp_var_1135(iplon_var_1162, i_lay_var_1159 - 1) < layreffr_var_219 .AND. k_jp_var_1135(iplon_var_1162, i_lay_var_1159) >= layreffr_var_219) i_laysolfr_var_1160(iplon_var_1162) = i_lay_var_1159
        ind0_var_1154 = ((k_jp_var_1135(iplon_var_1162, i_lay_var_1159) - 13) * 5 + (k_jt_var_1136(iplon_var_1162, i_lay_var_1159) - 1)) * nspb_var_314(16) + 1
        ind1_var_1155 = ((k_jp_var_1135(iplon_var_1162, i_lay_var_1159) - 12) * 5 + (k_jt1_var_1137(iplon_var_1162, i_lay_var_1159) - 1)) * nspb_var_314(16) + 1
        z_tauray_var_1169 = p_colmol_var_1141(iplon_var_1162, i_lay_var_1159) * rayl_var_218
        DO ig_var_1153 = 1, 6
          p_taug_var_1150(iplon_var_1162, i_lay_var_1159, ig_var_1153) = p_colch4_var_1140(iplon_var_1162, i_lay_var_1159) * (p_fac00_var_1131(iplon_var_1162, i_lay_var_1159) * absb_var_221(ind0_var_1154, ig_var_1153) + p_fac10_var_1133(iplon_var_1162, i_lay_var_1159) * absb_var_221(ind0_var_1154 + 1, ig_var_1153) + p_fac01_var_1132(iplon_var_1162, i_lay_var_1159) * absb_var_221(ind1_var_1155, ig_var_1153) + p_fac11_var_1134(iplon_var_1162, i_lay_var_1159) * absb_var_221(ind1_var_1155 + 1, ig_var_1153))
          IF (i_lay_var_1159 == i_laysolfr_var_1160(iplon_var_1162)) p_sfluxzen_var_1149(iplon_var_1162, ig_var_1153) = sfluxrefc_var_224(ig_var_1153)
          p_taur_var_1151(iplon_var_1162, i_lay_var_1159, ig_var_1153) = z_tauray_var_1169
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1159 = laytrop_max_var_1164 + 1, i_nlayers_var_1161
    DO iplon_var_1162 = kidia_var_1128, kfdia_var_1129
      IF (k_jp_var_1135(iplon_var_1162, i_lay_var_1159 - 1) < layreffr_var_219 .AND. k_jp_var_1135(iplon_var_1162, i_lay_var_1159) >= layreffr_var_219) i_laysolfr_var_1160(iplon_var_1162) = i_lay_var_1159
      ind0_var_1154 = ((k_jp_var_1135(iplon_var_1162, i_lay_var_1159) - 13) * 5 + (k_jt_var_1136(iplon_var_1162, i_lay_var_1159) - 1)) * nspb_var_314(16) + 1
      ind1_var_1155 = ((k_jp_var_1135(iplon_var_1162, i_lay_var_1159) - 12) * 5 + (k_jt1_var_1137(iplon_var_1162, i_lay_var_1159) - 1)) * nspb_var_314(16) + 1
      z_tauray_var_1169 = p_colmol_var_1141(iplon_var_1162, i_lay_var_1159) * rayl_var_218
      DO ig_var_1153 = 1, 6
        p_taug_var_1150(iplon_var_1162, i_lay_var_1159, ig_var_1153) = p_colch4_var_1140(iplon_var_1162, i_lay_var_1159) * (p_fac00_var_1131(iplon_var_1162, i_lay_var_1159) * absb_var_221(ind0_var_1154, ig_var_1153) + p_fac10_var_1133(iplon_var_1162, i_lay_var_1159) * absb_var_221(ind0_var_1154 + 1, ig_var_1153) + p_fac01_var_1132(iplon_var_1162, i_lay_var_1159) * absb_var_221(ind1_var_1155, ig_var_1153) + p_fac11_var_1134(iplon_var_1162, i_lay_var_1159) * absb_var_221(ind1_var_1155 + 1, ig_var_1153))
        IF (i_lay_var_1159 == i_laysolfr_var_1160(iplon_var_1162)) p_sfluxzen_var_1149(iplon_var_1162, ig_var_1153) = sfluxrefc_var_224(ig_var_1153)
        p_taur_var_1151(iplon_var_1162, i_lay_var_1159, ig_var_1153) = z_tauray_var_1169
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol16
SUBROUTINE srtm_taumol27(kidia_var_1170, kfdia_var_1171, klev_var_1172, p_fac00_var_1173, p_fac01_var_1174, p_fac10_var_1175, p_fac11_var_1176, k_jp_var_1177, k_jt_var_1178, k_jt1_var_1179, p_colmol_var_1180, p_colo3_var_1181, k_laytrop_var_1182, p_sfluxzen_var_1183, p_taug_var_1184, p_taur_var_1185, prmu0_var_1186)
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  USE yoesrta27, ONLY: absa_var_296, absb_var_297, layreffr_var_295, raylc_var_299, scalekur, sfluxrefc_var_298
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1170, kfdia_var_1171
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1172
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1173(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1174(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1175(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1176(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1177(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1178(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1179(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1180(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1181(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1182(kidia_var_1170 : kfdia_var_1171)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1183(kidia_var_1170 : kfdia_var_1171, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1184(kidia_var_1170 : kfdia_var_1171, klev_var_1172, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1185(kidia_var_1170 : kfdia_var_1171, klev_var_1172, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1186(kidia_var_1170 : kfdia_var_1171)
  INTEGER(KIND = 4) :: ig_var_1187, ind0_var_1188, ind1_var_1189, i_lay_var_1190, i_laysolfr_var_1191(kidia_var_1170 : kfdia_var_1171), i_nlayers_var_1192, iplon_var_1193
  INTEGER(KIND = 4) :: laytrop_min_var_1194, laytrop_max_var_1195
  REAL(KIND = 8) :: z_tauray_var_1196
  laytrop_min_var_1194 = MINVAL(k_laytrop_var_1182(kidia_var_1170 : kfdia_var_1171))
  laytrop_max_var_1195 = MAXVAL(k_laytrop_var_1182(kidia_var_1170 : kfdia_var_1171))
  i_nlayers_var_1192 = klev_var_1172
  DO iplon_var_1193 = kidia_var_1170, kfdia_var_1171
    i_laysolfr_var_1191(iplon_var_1193) = i_nlayers_var_1192
  END DO
  DO i_lay_var_1190 = 1, laytrop_min_var_1194
    DO iplon_var_1193 = kidia_var_1170, kfdia_var_1171
      ind0_var_1188 = ((k_jp_var_1177(iplon_var_1193, i_lay_var_1190) - 1) * 5 + (k_jt_var_1178(iplon_var_1193, i_lay_var_1190) - 1)) * nspa_var_313(27) + 1
      ind1_var_1189 = (k_jp_var_1177(iplon_var_1193, i_lay_var_1190) * 5 + (k_jt1_var_1179(iplon_var_1193, i_lay_var_1190) - 1)) * nspa_var_313(27) + 1
      DO ig_var_1187 = 1, 8
        z_tauray_var_1196 = p_colmol_var_1180(iplon_var_1193, i_lay_var_1190) * raylc_var_299(ig_var_1187)
        p_taug_var_1184(iplon_var_1193, i_lay_var_1190, ig_var_1187) = p_colo3_var_1181(iplon_var_1193, i_lay_var_1190) * (p_fac00_var_1173(iplon_var_1193, i_lay_var_1190) * absa_var_296(ind0_var_1188, ig_var_1187) + p_fac10_var_1175(iplon_var_1193, i_lay_var_1190) * absa_var_296(ind0_var_1188 + 1, ig_var_1187) + p_fac01_var_1174(iplon_var_1193, i_lay_var_1190) * absa_var_296(ind1_var_1189, ig_var_1187) + p_fac11_var_1176(iplon_var_1193, i_lay_var_1190) * absa_var_296(ind1_var_1189 + 1, ig_var_1187))
        p_taur_var_1185(iplon_var_1193, i_lay_var_1190, ig_var_1187) = z_tauray_var_1196
      END DO
    END DO
  END DO
  DO i_lay_var_1190 = laytrop_min_var_1194 + 1, laytrop_max_var_1195
    DO iplon_var_1193 = kidia_var_1170, kfdia_var_1171
      IF (i_lay_var_1190 <= k_laytrop_var_1182(iplon_var_1193)) THEN
        ind0_var_1188 = ((k_jp_var_1177(iplon_var_1193, i_lay_var_1190) - 1) * 5 + (k_jt_var_1178(iplon_var_1193, i_lay_var_1190) - 1)) * nspa_var_313(27) + 1
        ind1_var_1189 = (k_jp_var_1177(iplon_var_1193, i_lay_var_1190) * 5 + (k_jt1_var_1179(iplon_var_1193, i_lay_var_1190) - 1)) * nspa_var_313(27) + 1
        DO ig_var_1187 = 1, 8
          z_tauray_var_1196 = p_colmol_var_1180(iplon_var_1193, i_lay_var_1190) * raylc_var_299(ig_var_1187)
          p_taug_var_1184(iplon_var_1193, i_lay_var_1190, ig_var_1187) = p_colo3_var_1181(iplon_var_1193, i_lay_var_1190) * (p_fac00_var_1173(iplon_var_1193, i_lay_var_1190) * absa_var_296(ind0_var_1188, ig_var_1187) + p_fac10_var_1175(iplon_var_1193, i_lay_var_1190) * absa_var_296(ind0_var_1188 + 1, ig_var_1187) + p_fac01_var_1174(iplon_var_1193, i_lay_var_1190) * absa_var_296(ind1_var_1189, ig_var_1187) + p_fac11_var_1176(iplon_var_1193, i_lay_var_1190) * absa_var_296(ind1_var_1189 + 1, ig_var_1187))
          p_taur_var_1185(iplon_var_1193, i_lay_var_1190, ig_var_1187) = z_tauray_var_1196
        END DO
      ELSE
        IF (k_jp_var_1177(iplon_var_1193, i_lay_var_1190 - 1) < layreffr_var_295 .AND. k_jp_var_1177(iplon_var_1193, i_lay_var_1190) >= layreffr_var_295) i_laysolfr_var_1191(iplon_var_1193) = i_lay_var_1190
        ind0_var_1188 = ((k_jp_var_1177(iplon_var_1193, i_lay_var_1190) - 13) * 5 + (k_jt_var_1178(iplon_var_1193, i_lay_var_1190) - 1)) * nspb_var_314(27) + 1
        ind1_var_1189 = ((k_jp_var_1177(iplon_var_1193, i_lay_var_1190) - 12) * 5 + (k_jt1_var_1179(iplon_var_1193, i_lay_var_1190) - 1)) * nspb_var_314(27) + 1
        DO ig_var_1187 = 1, 8
          z_tauray_var_1196 = p_colmol_var_1180(iplon_var_1193, i_lay_var_1190) * raylc_var_299(ig_var_1187)
          p_taug_var_1184(iplon_var_1193, i_lay_var_1190, ig_var_1187) = p_colo3_var_1181(iplon_var_1193, i_lay_var_1190) * (p_fac00_var_1173(iplon_var_1193, i_lay_var_1190) * absb_var_297(ind0_var_1188, ig_var_1187) + p_fac10_var_1175(iplon_var_1193, i_lay_var_1190) * absb_var_297(ind0_var_1188 + 1, ig_var_1187) + p_fac01_var_1174(iplon_var_1193, i_lay_var_1190) * absb_var_297(ind1_var_1189, ig_var_1187) + p_fac11_var_1176(iplon_var_1193, i_lay_var_1190) * absb_var_297(ind1_var_1189 + 1, ig_var_1187))
          IF (i_lay_var_1190 == i_laysolfr_var_1191(iplon_var_1193)) p_sfluxzen_var_1183(iplon_var_1193, ig_var_1187) = scalekur * sfluxrefc_var_298(ig_var_1187)
          p_taur_var_1185(iplon_var_1193, i_lay_var_1190, ig_var_1187) = z_tauray_var_1196
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1190 = laytrop_max_var_1195 + 1, i_nlayers_var_1192
    DO iplon_var_1193 = kidia_var_1170, kfdia_var_1171
      IF (k_jp_var_1177(iplon_var_1193, i_lay_var_1190 - 1) < layreffr_var_295 .AND. k_jp_var_1177(iplon_var_1193, i_lay_var_1190) >= layreffr_var_295) i_laysolfr_var_1191(iplon_var_1193) = i_lay_var_1190
      ind0_var_1188 = ((k_jp_var_1177(iplon_var_1193, i_lay_var_1190) - 13) * 5 + (k_jt_var_1178(iplon_var_1193, i_lay_var_1190) - 1)) * nspb_var_314(27) + 1
      ind1_var_1189 = ((k_jp_var_1177(iplon_var_1193, i_lay_var_1190) - 12) * 5 + (k_jt1_var_1179(iplon_var_1193, i_lay_var_1190) - 1)) * nspb_var_314(27) + 1
      DO ig_var_1187 = 1, 8
        z_tauray_var_1196 = p_colmol_var_1180(iplon_var_1193, i_lay_var_1190) * raylc_var_299(ig_var_1187)
        p_taug_var_1184(iplon_var_1193, i_lay_var_1190, ig_var_1187) = p_colo3_var_1181(iplon_var_1193, i_lay_var_1190) * (p_fac00_var_1173(iplon_var_1193, i_lay_var_1190) * absb_var_297(ind0_var_1188, ig_var_1187) + p_fac10_var_1175(iplon_var_1193, i_lay_var_1190) * absb_var_297(ind0_var_1188 + 1, ig_var_1187) + p_fac01_var_1174(iplon_var_1193, i_lay_var_1190) * absb_var_297(ind1_var_1189, ig_var_1187) + p_fac11_var_1176(iplon_var_1193, i_lay_var_1190) * absb_var_297(ind1_var_1189 + 1, ig_var_1187))
        IF (i_lay_var_1190 == i_laysolfr_var_1191(iplon_var_1193)) p_sfluxzen_var_1183(iplon_var_1193, ig_var_1187) = scalekur * sfluxrefc_var_298(ig_var_1187)
        p_taur_var_1185(iplon_var_1193, i_lay_var_1190, ig_var_1187) = z_tauray_var_1196
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol27
SUBROUTINE srtm_taumol26(kidia_var_1197, kfdia_var_1198, klev_var_1199, p_colmol_var_1200, k_laytrop_var_1201, p_sfluxzen_var_1202, p_taug_var_1203, p_taur_var_1204, prmu0_var_1205)
  USE yoesrta26, ONLY: raylc_var_294, sfluxrefc_var_293
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1197, kfdia_var_1198
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1199
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1200(kidia_var_1197 : kfdia_var_1198, klev_var_1199)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1201(kidia_var_1197 : kfdia_var_1198)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1202(kidia_var_1197 : kfdia_var_1198, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1203(kidia_var_1197 : kfdia_var_1198, klev_var_1199, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1204(kidia_var_1197 : kfdia_var_1198, klev_var_1199, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1205(kidia_var_1197 : kfdia_var_1198)
  INTEGER(KIND = 4) :: ig_var_1206, i_lay_var_1207, i_laysolfr_var_1208(kidia_var_1197 : kfdia_var_1198), i_nlayers_var_1209, iplon_var_1210
  INTEGER(KIND = 4) :: laytrop_min_var_1211, laytrop_max_var_1212
  laytrop_min_var_1211 = MINVAL(k_laytrop_var_1201(kidia_var_1197 : kfdia_var_1198))
  laytrop_max_var_1212 = MAXVAL(k_laytrop_var_1201(kidia_var_1197 : kfdia_var_1198))
  i_nlayers_var_1209 = klev_var_1199
  DO iplon_var_1210 = kidia_var_1197, kfdia_var_1198
    i_laysolfr_var_1208(iplon_var_1210) = k_laytrop_var_1201(iplon_var_1210)
  END DO
  DO i_lay_var_1207 = 1, laytrop_min_var_1211
    DO iplon_var_1210 = kidia_var_1197, kfdia_var_1198
      DO ig_var_1206 = 1, 6
        IF (i_lay_var_1207 == i_laysolfr_var_1208(iplon_var_1210)) p_sfluxzen_var_1202(iplon_var_1210, ig_var_1206) = sfluxrefc_var_293(ig_var_1206)
        p_taug_var_1203(iplon_var_1210, i_lay_var_1207, ig_var_1206) = 0.0D0
        p_taur_var_1204(iplon_var_1210, i_lay_var_1207, ig_var_1206) = p_colmol_var_1200(iplon_var_1210, i_lay_var_1207) * raylc_var_294(ig_var_1206)
      END DO
    END DO
  END DO
  DO i_lay_var_1207 = laytrop_min_var_1211 + 1, laytrop_max_var_1212
    DO iplon_var_1210 = kidia_var_1197, kfdia_var_1198
      IF (i_lay_var_1207 <= k_laytrop_var_1201(iplon_var_1210)) THEN
        DO ig_var_1206 = 1, 6
          IF (i_lay_var_1207 == i_laysolfr_var_1208(iplon_var_1210)) p_sfluxzen_var_1202(iplon_var_1210, ig_var_1206) = sfluxrefc_var_293(ig_var_1206)
          p_taug_var_1203(iplon_var_1210, i_lay_var_1207, ig_var_1206) = 0.0D0
          p_taur_var_1204(iplon_var_1210, i_lay_var_1207, ig_var_1206) = p_colmol_var_1200(iplon_var_1210, i_lay_var_1207) * raylc_var_294(ig_var_1206)
        END DO
      ELSE
        DO ig_var_1206 = 1, 6
          p_taug_var_1203(iplon_var_1210, i_lay_var_1207, ig_var_1206) = 0.0D0
          p_taur_var_1204(iplon_var_1210, i_lay_var_1207, ig_var_1206) = p_colmol_var_1200(iplon_var_1210, i_lay_var_1207) * raylc_var_294(ig_var_1206)
        END DO
      END IF
    END DO
  END DO
  DO ig_var_1206 = 1, 6
    DO i_lay_var_1207 = laytrop_max_var_1212 + 1, i_nlayers_var_1209
      DO iplon_var_1210 = kidia_var_1197, kfdia_var_1198
        p_taug_var_1203(iplon_var_1210, i_lay_var_1207, ig_var_1206) = 0.0D0
        p_taur_var_1204(iplon_var_1210, i_lay_var_1207, ig_var_1206) = p_colmol_var_1200(iplon_var_1210, i_lay_var_1207) * raylc_var_294(ig_var_1206)
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol26
SUBROUTINE srtm_taumol18(kidia_var_1213, kfdia_var_1214, klev_var_1215, p_fac00_var_1216, p_fac01_var_1217, p_fac10_var_1218, p_fac11_var_1219, k_jp_var_1220, k_jt_var_1221, k_jt1_var_1222, p_oneminus_var_1223, p_colh2o_var_1224, p_colch4_var_1225, p_colmol_var_1226, k_laytrop_var_1227, p_selffac_var_1228, p_selffrac_var_1229, k_indself_var_1230, p_forfac_var_1231, p_forfrac_var_1232, k_indfor_var_1233, p_sfluxzen_var_1234, p_taug_var_1235, p_taur_var_1236, prmu0_var_1237)
  USE yoesrta18, ONLY: absa_var_236, absb_var_237, forrefc_var_239, layreffr_var_235, rayl_var_233, selfrefc_var_238, sfluxrefc_var_240, strrat_var_234
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1213, kfdia_var_1214
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1215
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1216(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1217(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1218(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1219(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1220(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1221(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1222(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1223(kidia_var_1213 : kfdia_var_1214)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1224(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_1225(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1226(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1227(kidia_var_1213 : kfdia_var_1214)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1228(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1229(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1230(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1231(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1232(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1233(kidia_var_1213 : kfdia_var_1214, klev_var_1215)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1234(kidia_var_1213 : kfdia_var_1214, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1235(kidia_var_1213 : kfdia_var_1214, klev_var_1215, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1236(kidia_var_1213 : kfdia_var_1214, klev_var_1215, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1237(kidia_var_1213 : kfdia_var_1214)
  INTEGER(KIND = 4) :: ig_var_1238, ind0_var_1239, ind1_var_1240, inds_var_1241, indf_var_1242, js_var_1243, i_lay_var_1244, i_laysolfr_var_1245(kidia_var_1213 : kfdia_var_1214), i_nlayers_var_1246, iplon_var_1247
  INTEGER(KIND = 4) :: laytrop_min_var_1248, laytrop_max_var_1249
  REAL(KIND = 8) :: z_fs_var_1250, z_speccomb_var_1251, z_specmult_var_1252, z_specparm_var_1253, z_tauray_var_1254
  laytrop_min_var_1248 = MINVAL(k_laytrop_var_1227(kidia_var_1213 : kfdia_var_1214))
  laytrop_max_var_1249 = MAXVAL(k_laytrop_var_1227(kidia_var_1213 : kfdia_var_1214))
  i_nlayers_var_1246 = klev_var_1215
  DO iplon_var_1247 = kidia_var_1213, kfdia_var_1214
    i_laysolfr_var_1245(iplon_var_1247) = k_laytrop_var_1227(iplon_var_1247)
  END DO
  DO i_lay_var_1244 = 1, laytrop_min_var_1248
    DO iplon_var_1247 = kidia_var_1213, kfdia_var_1214
      IF (k_jp_var_1220(iplon_var_1247, i_lay_var_1244) < layreffr_var_235 .AND. k_jp_var_1220(iplon_var_1247, i_lay_var_1244 + 1) >= layreffr_var_235) i_laysolfr_var_1245(iplon_var_1247) = MIN(i_lay_var_1244 + 1, k_laytrop_var_1227(iplon_var_1247))
      z_speccomb_var_1251 = p_colh2o_var_1224(iplon_var_1247, i_lay_var_1244) + strrat_var_234 * p_colch4_var_1225(iplon_var_1247, i_lay_var_1244)
      z_specparm_var_1253 = p_colh2o_var_1224(iplon_var_1247, i_lay_var_1244) / z_speccomb_var_1251
      z_specparm_var_1253 = MIN(p_oneminus_var_1223(iplon_var_1247), z_specparm_var_1253)
      z_specmult_var_1252 = 8.0D0 * (z_specparm_var_1253)
      js_var_1243 = 1 + INT(z_specmult_var_1252)
      z_fs_var_1250 = z_specmult_var_1252 - AINT(z_specmult_var_1252)
      ind0_var_1239 = ((k_jp_var_1220(iplon_var_1247, i_lay_var_1244) - 1) * 5 + (k_jt_var_1221(iplon_var_1247, i_lay_var_1244) - 1)) * nspa_var_313(18) + js_var_1243
      ind1_var_1240 = (k_jp_var_1220(iplon_var_1247, i_lay_var_1244) * 5 + (k_jt1_var_1222(iplon_var_1247, i_lay_var_1244) - 1)) * nspa_var_313(18) + js_var_1243
      inds_var_1241 = k_indself_var_1230(iplon_var_1247, i_lay_var_1244)
      indf_var_1242 = k_indfor_var_1233(iplon_var_1247, i_lay_var_1244)
      z_tauray_var_1254 = p_colmol_var_1226(iplon_var_1247, i_lay_var_1244) * rayl_var_233
      DO ig_var_1238 = 1, 8
        p_taug_var_1235(iplon_var_1247, i_lay_var_1244, ig_var_1238) = z_speccomb_var_1251 * ((1.0D0 - z_fs_var_1250) * (absa_var_236(ind0_var_1239, ig_var_1238) * p_fac00_var_1216(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind0_var_1239 + 9, ig_var_1238) * p_fac10_var_1218(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind1_var_1240, ig_var_1238) * p_fac01_var_1217(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind1_var_1240 + 9, ig_var_1238) * p_fac11_var_1219(iplon_var_1247, i_lay_var_1244)) + z_fs_var_1250 * (absa_var_236(ind0_var_1239 + 1, ig_var_1238) * p_fac00_var_1216(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind0_var_1239 + 10, ig_var_1238) * p_fac10_var_1218(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind1_var_1240 + 1, ig_var_1238) * p_fac01_var_1217(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind1_var_1240 + 10, ig_var_1238) * p_fac11_var_1219(iplon_var_1247, i_lay_var_1244))) + p_colh2o_var_1224(iplon_var_1247, i_lay_var_1244) * (p_selffac_var_1228(iplon_var_1247, i_lay_var_1244) * (selfrefc_var_238(inds_var_1241, ig_var_1238) + p_selffrac_var_1229(iplon_var_1247, i_lay_var_1244) * (selfrefc_var_238(inds_var_1241 + 1, ig_var_1238) - selfrefc_var_238(inds_var_1241, ig_var_1238))) + p_forfac_var_1231(iplon_var_1247, i_lay_var_1244) * (forrefc_var_239(indf_var_1242, ig_var_1238) + p_forfrac_var_1232(iplon_var_1247, i_lay_var_1244) * (forrefc_var_239(indf_var_1242 + 1, ig_var_1238) - forrefc_var_239(indf_var_1242, ig_var_1238))))
        IF (i_lay_var_1244 == i_laysolfr_var_1245(iplon_var_1247)) p_sfluxzen_var_1234(iplon_var_1247, ig_var_1238) = sfluxrefc_var_240(ig_var_1238, js_var_1243) + z_fs_var_1250 * (sfluxrefc_var_240(ig_var_1238, js_var_1243 + 1) - sfluxrefc_var_240(ig_var_1238, js_var_1243))
        p_taur_var_1236(iplon_var_1247, i_lay_var_1244, ig_var_1238) = z_tauray_var_1254
      END DO
    END DO
  END DO
  DO i_lay_var_1244 = laytrop_min_var_1248 + 1, laytrop_max_var_1249
    DO iplon_var_1247 = kidia_var_1213, kfdia_var_1214
      IF (i_lay_var_1244 <= k_laytrop_var_1227(iplon_var_1247)) THEN
        IF (k_jp_var_1220(iplon_var_1247, i_lay_var_1244) < layreffr_var_235 .AND. k_jp_var_1220(iplon_var_1247, i_lay_var_1244 + 1) >= layreffr_var_235) i_laysolfr_var_1245(iplon_var_1247) = MIN(i_lay_var_1244 + 1, k_laytrop_var_1227(iplon_var_1247))
        z_speccomb_var_1251 = p_colh2o_var_1224(iplon_var_1247, i_lay_var_1244) + strrat_var_234 * p_colch4_var_1225(iplon_var_1247, i_lay_var_1244)
        z_specparm_var_1253 = p_colh2o_var_1224(iplon_var_1247, i_lay_var_1244) / z_speccomb_var_1251
        z_specparm_var_1253 = MIN(p_oneminus_var_1223(iplon_var_1247), z_specparm_var_1253)
        z_specmult_var_1252 = 8.0D0 * (z_specparm_var_1253)
        js_var_1243 = 1 + INT(z_specmult_var_1252)
        z_fs_var_1250 = z_specmult_var_1252 - AINT(z_specmult_var_1252)
        ind0_var_1239 = ((k_jp_var_1220(iplon_var_1247, i_lay_var_1244) - 1) * 5 + (k_jt_var_1221(iplon_var_1247, i_lay_var_1244) - 1)) * nspa_var_313(18) + js_var_1243
        ind1_var_1240 = (k_jp_var_1220(iplon_var_1247, i_lay_var_1244) * 5 + (k_jt1_var_1222(iplon_var_1247, i_lay_var_1244) - 1)) * nspa_var_313(18) + js_var_1243
        inds_var_1241 = k_indself_var_1230(iplon_var_1247, i_lay_var_1244)
        indf_var_1242 = k_indfor_var_1233(iplon_var_1247, i_lay_var_1244)
        z_tauray_var_1254 = p_colmol_var_1226(iplon_var_1247, i_lay_var_1244) * rayl_var_233
        DO ig_var_1238 = 1, 8
          p_taug_var_1235(iplon_var_1247, i_lay_var_1244, ig_var_1238) = z_speccomb_var_1251 * ((1.0D0 - z_fs_var_1250) * (absa_var_236(ind0_var_1239, ig_var_1238) * p_fac00_var_1216(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind0_var_1239 + 9, ig_var_1238) * p_fac10_var_1218(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind1_var_1240, ig_var_1238) * p_fac01_var_1217(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind1_var_1240 + 9, ig_var_1238) * p_fac11_var_1219(iplon_var_1247, i_lay_var_1244)) + z_fs_var_1250 * (absa_var_236(ind0_var_1239 + 1, ig_var_1238) * p_fac00_var_1216(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind0_var_1239 + 10, ig_var_1238) * p_fac10_var_1218(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind1_var_1240 + 1, ig_var_1238) * p_fac01_var_1217(iplon_var_1247, i_lay_var_1244) + absa_var_236(ind1_var_1240 + 10, ig_var_1238) * p_fac11_var_1219(iplon_var_1247, i_lay_var_1244))) + p_colh2o_var_1224(iplon_var_1247, i_lay_var_1244) * (p_selffac_var_1228(iplon_var_1247, i_lay_var_1244) * (selfrefc_var_238(inds_var_1241, ig_var_1238) + p_selffrac_var_1229(iplon_var_1247, i_lay_var_1244) * (selfrefc_var_238(inds_var_1241 + 1, ig_var_1238) - selfrefc_var_238(inds_var_1241, ig_var_1238))) + p_forfac_var_1231(iplon_var_1247, i_lay_var_1244) * (forrefc_var_239(indf_var_1242, ig_var_1238) + p_forfrac_var_1232(iplon_var_1247, i_lay_var_1244) * (forrefc_var_239(indf_var_1242 + 1, ig_var_1238) - forrefc_var_239(indf_var_1242, ig_var_1238))))
          IF (i_lay_var_1244 == i_laysolfr_var_1245(iplon_var_1247)) p_sfluxzen_var_1234(iplon_var_1247, ig_var_1238) = sfluxrefc_var_240(ig_var_1238, js_var_1243) + z_fs_var_1250 * (sfluxrefc_var_240(ig_var_1238, js_var_1243 + 1) - sfluxrefc_var_240(ig_var_1238, js_var_1243))
          p_taur_var_1236(iplon_var_1247, i_lay_var_1244, ig_var_1238) = z_tauray_var_1254
        END DO
      ELSE
        ind0_var_1239 = ((k_jp_var_1220(iplon_var_1247, i_lay_var_1244) - 13) * 5 + (k_jt_var_1221(iplon_var_1247, i_lay_var_1244) - 1)) * nspb_var_314(18) + 1
        ind1_var_1240 = ((k_jp_var_1220(iplon_var_1247, i_lay_var_1244) - 12) * 5 + (k_jt1_var_1222(iplon_var_1247, i_lay_var_1244) - 1)) * nspb_var_314(18) + 1
        z_tauray_var_1254 = p_colmol_var_1226(iplon_var_1247, i_lay_var_1244) * rayl_var_233
        DO ig_var_1238 = 1, 8
          p_taug_var_1235(iplon_var_1247, i_lay_var_1244, ig_var_1238) = p_colch4_var_1225(iplon_var_1247, i_lay_var_1244) * (p_fac00_var_1216(iplon_var_1247, i_lay_var_1244) * absb_var_237(ind0_var_1239, ig_var_1238) + p_fac10_var_1218(iplon_var_1247, i_lay_var_1244) * absb_var_237(ind0_var_1239 + 1, ig_var_1238) + p_fac01_var_1217(iplon_var_1247, i_lay_var_1244) * absb_var_237(ind1_var_1240, ig_var_1238) + p_fac11_var_1219(iplon_var_1247, i_lay_var_1244) * absb_var_237(ind1_var_1240 + 1, ig_var_1238))
          p_taur_var_1236(iplon_var_1247, i_lay_var_1244, ig_var_1238) = z_tauray_var_1254
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1244 = laytrop_max_var_1249 + 1, i_nlayers_var_1246
    DO iplon_var_1247 = kidia_var_1213, kfdia_var_1214
      ind0_var_1239 = ((k_jp_var_1220(iplon_var_1247, i_lay_var_1244) - 13) * 5 + (k_jt_var_1221(iplon_var_1247, i_lay_var_1244) - 1)) * nspb_var_314(18) + 1
      ind1_var_1240 = ((k_jp_var_1220(iplon_var_1247, i_lay_var_1244) - 12) * 5 + (k_jt1_var_1222(iplon_var_1247, i_lay_var_1244) - 1)) * nspb_var_314(18) + 1
      z_tauray_var_1254 = p_colmol_var_1226(iplon_var_1247, i_lay_var_1244) * rayl_var_233
      DO ig_var_1238 = 1, 8
        p_taug_var_1235(iplon_var_1247, i_lay_var_1244, ig_var_1238) = p_colch4_var_1225(iplon_var_1247, i_lay_var_1244) * (p_fac00_var_1216(iplon_var_1247, i_lay_var_1244) * absb_var_237(ind0_var_1239, ig_var_1238) + p_fac10_var_1218(iplon_var_1247, i_lay_var_1244) * absb_var_237(ind0_var_1239 + 1, ig_var_1238) + p_fac01_var_1217(iplon_var_1247, i_lay_var_1244) * absb_var_237(ind1_var_1240, ig_var_1238) + p_fac11_var_1219(iplon_var_1247, i_lay_var_1244) * absb_var_237(ind1_var_1240 + 1, ig_var_1238))
        p_taur_var_1236(iplon_var_1247, i_lay_var_1244, ig_var_1238) = z_tauray_var_1254
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol18
SUBROUTINE srtm_taumol24(kidia_var_1255, kfdia_var_1256, klev_var_1257, p_fac00_var_1258, p_fac01_var_1259, p_fac10_var_1260, p_fac11_var_1261, k_jp_var_1262, k_jt_var_1263, k_jt1_var_1264, p_oneminus_var_1265, p_colh2o_var_1266, p_colmol_var_1267, p_colo2_var_1268, p_colo3_var_1269, k_laytrop_var_1270, p_selffac_var_1271, p_selffrac_var_1272, k_indself_var_1273, p_forfac_var_1274, p_forfrac_var_1275, k_indfor_var_1276, p_sfluxzen_var_1277, p_taug_var_1278, p_taur_var_1279, prmu0_var_1280)
  USE yoesrta24, ONLY: absa_var_280, absb_var_281, abso3ac_var_285, abso3bc_var_286, forrefc_var_283, layreffr_var_279, raylac, raylbc, selfrefc_var_282, sfluxrefc_var_284, strrat_var_278
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1255, kfdia_var_1256
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1257
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1258(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1259(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1260(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1261(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1262(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1263(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1264(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1265(kidia_var_1255 : kfdia_var_1256)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1266(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1267(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_1268(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1269(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1270(kidia_var_1255 : kfdia_var_1256)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1271(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1272(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1273(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1274(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1275(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1276(kidia_var_1255 : kfdia_var_1256, klev_var_1257)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1277(kidia_var_1255 : kfdia_var_1256, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1278(kidia_var_1255 : kfdia_var_1256, klev_var_1257, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1279(kidia_var_1255 : kfdia_var_1256, klev_var_1257, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1280(kidia_var_1255 : kfdia_var_1256)
  INTEGER(KIND = 4) :: ig_var_1281, ind0_var_1282, ind1_var_1283, inds_var_1284, indf_var_1285, js_var_1286, i_lay_var_1287, i_laysolfr_var_1288(kidia_var_1255 : kfdia_var_1256), i_nlayers_var_1289, iplon_var_1290
  INTEGER(KIND = 4) :: laytrop_min_var_1291, laytrop_max_var_1292
  REAL(KIND = 8) :: z_fs_var_1293, z_speccomb_var_1294, z_specmult_var_1295, z_specparm_var_1296, z_tauray_var_1297
  laytrop_min_var_1291 = MINVAL(k_laytrop_var_1270(kidia_var_1255 : kfdia_var_1256))
  laytrop_max_var_1292 = MAXVAL(k_laytrop_var_1270(kidia_var_1255 : kfdia_var_1256))
  i_nlayers_var_1289 = klev_var_1257
  DO iplon_var_1290 = kidia_var_1255, kfdia_var_1256
    i_laysolfr_var_1288(iplon_var_1290) = k_laytrop_var_1270(iplon_var_1290)
  END DO
  DO i_lay_var_1287 = 1, laytrop_min_var_1291
    DO iplon_var_1290 = kidia_var_1255, kfdia_var_1256
      IF (k_jp_var_1262(iplon_var_1290, i_lay_var_1287) < layreffr_var_279 .AND. k_jp_var_1262(iplon_var_1290, i_lay_var_1287 + 1) >= layreffr_var_279) i_laysolfr_var_1288(iplon_var_1290) = MIN(i_lay_var_1287 + 1, k_laytrop_var_1270(iplon_var_1290))
      z_speccomb_var_1294 = p_colh2o_var_1266(iplon_var_1290, i_lay_var_1287) + strrat_var_278 * p_colo2_var_1268(iplon_var_1290, i_lay_var_1287)
      z_specparm_var_1296 = p_colh2o_var_1266(iplon_var_1290, i_lay_var_1287) / z_speccomb_var_1294
      z_specparm_var_1296 = MIN(p_oneminus_var_1265(iplon_var_1290), z_specparm_var_1296)
      z_specmult_var_1295 = 8.0D0 * (z_specparm_var_1296)
      js_var_1286 = 1 + INT(z_specmult_var_1295)
      z_fs_var_1293 = z_specmult_var_1295 - AINT(z_specmult_var_1295)
      ind0_var_1282 = ((k_jp_var_1262(iplon_var_1290, i_lay_var_1287) - 1) * 5 + (k_jt_var_1263(iplon_var_1290, i_lay_var_1287) - 1)) * nspa_var_313(24) + js_var_1286
      ind1_var_1283 = (k_jp_var_1262(iplon_var_1290, i_lay_var_1287) * 5 + (k_jt1_var_1264(iplon_var_1290, i_lay_var_1287) - 1)) * nspa_var_313(24) + js_var_1286
      inds_var_1284 = k_indself_var_1273(iplon_var_1290, i_lay_var_1287)
      indf_var_1285 = k_indfor_var_1276(iplon_var_1290, i_lay_var_1287)
      DO ig_var_1281 = 1, 8
        z_tauray_var_1297 = p_colmol_var_1267(iplon_var_1290, i_lay_var_1287) * (raylac(ig_var_1281, js_var_1286) + z_fs_var_1293 * (raylac(ig_var_1281, js_var_1286 + 1) - raylac(ig_var_1281, js_var_1286)))
        p_taug_var_1278(iplon_var_1290, i_lay_var_1287, ig_var_1281) = z_speccomb_var_1294 * ((1.0D0 - z_fs_var_1293) * (absa_var_280(ind0_var_1282, ig_var_1281) * p_fac00_var_1258(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind0_var_1282 + 9, ig_var_1281) * p_fac10_var_1260(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind1_var_1283, ig_var_1281) * p_fac01_var_1259(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind1_var_1283 + 9, ig_var_1281) * p_fac11_var_1261(iplon_var_1290, i_lay_var_1287)) + z_fs_var_1293 * (absa_var_280(ind0_var_1282 + 1, ig_var_1281) * p_fac00_var_1258(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind0_var_1282 + 10, ig_var_1281) * p_fac10_var_1260(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind1_var_1283 + 1, ig_var_1281) * p_fac01_var_1259(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind1_var_1283 + 10, ig_var_1281) * p_fac11_var_1261(iplon_var_1290, i_lay_var_1287))) + p_colo3_var_1269(iplon_var_1290, i_lay_var_1287) * abso3ac_var_285(ig_var_1281) + p_colh2o_var_1266(iplon_var_1290, i_lay_var_1287) * (p_selffac_var_1271(iplon_var_1290, i_lay_var_1287) * (selfrefc_var_282(inds_var_1284, ig_var_1281) + p_selffrac_var_1272(iplon_var_1290, i_lay_var_1287) * (selfrefc_var_282(inds_var_1284 + 1, ig_var_1281) - selfrefc_var_282(inds_var_1284, ig_var_1281))) + p_forfac_var_1274(iplon_var_1290, i_lay_var_1287) * (forrefc_var_283(indf_var_1285, ig_var_1281) + p_forfrac_var_1275(iplon_var_1290, i_lay_var_1287) * (forrefc_var_283(indf_var_1285 + 1, ig_var_1281) - forrefc_var_283(indf_var_1285, ig_var_1281))))
        IF (i_lay_var_1287 == i_laysolfr_var_1288(iplon_var_1290)) p_sfluxzen_var_1277(iplon_var_1290, ig_var_1281) = sfluxrefc_var_284(ig_var_1281, js_var_1286) + z_fs_var_1293 * (sfluxrefc_var_284(ig_var_1281, js_var_1286 + 1) - sfluxrefc_var_284(ig_var_1281, js_var_1286))
        p_taur_var_1279(iplon_var_1290, i_lay_var_1287, ig_var_1281) = z_tauray_var_1297
      END DO
    END DO
  END DO
  DO i_lay_var_1287 = laytrop_min_var_1291 + 1, laytrop_max_var_1292
    DO iplon_var_1290 = kidia_var_1255, kfdia_var_1256
      IF (i_lay_var_1287 <= k_laytrop_var_1270(iplon_var_1290)) THEN
        IF (k_jp_var_1262(iplon_var_1290, i_lay_var_1287) < layreffr_var_279 .AND. k_jp_var_1262(iplon_var_1290, i_lay_var_1287 + 1) >= layreffr_var_279) i_laysolfr_var_1288(iplon_var_1290) = MIN(i_lay_var_1287 + 1, k_laytrop_var_1270(iplon_var_1290))
        z_speccomb_var_1294 = p_colh2o_var_1266(iplon_var_1290, i_lay_var_1287) + strrat_var_278 * p_colo2_var_1268(iplon_var_1290, i_lay_var_1287)
        z_specparm_var_1296 = p_colh2o_var_1266(iplon_var_1290, i_lay_var_1287) / z_speccomb_var_1294
        z_specparm_var_1296 = MIN(p_oneminus_var_1265(iplon_var_1290), z_specparm_var_1296)
        z_specmult_var_1295 = 8.0D0 * (z_specparm_var_1296)
        js_var_1286 = 1 + INT(z_specmult_var_1295)
        z_fs_var_1293 = z_specmult_var_1295 - AINT(z_specmult_var_1295)
        ind0_var_1282 = ((k_jp_var_1262(iplon_var_1290, i_lay_var_1287) - 1) * 5 + (k_jt_var_1263(iplon_var_1290, i_lay_var_1287) - 1)) * nspa_var_313(24) + js_var_1286
        ind1_var_1283 = (k_jp_var_1262(iplon_var_1290, i_lay_var_1287) * 5 + (k_jt1_var_1264(iplon_var_1290, i_lay_var_1287) - 1)) * nspa_var_313(24) + js_var_1286
        inds_var_1284 = k_indself_var_1273(iplon_var_1290, i_lay_var_1287)
        indf_var_1285 = k_indfor_var_1276(iplon_var_1290, i_lay_var_1287)
        DO ig_var_1281 = 1, 8
          z_tauray_var_1297 = p_colmol_var_1267(iplon_var_1290, i_lay_var_1287) * (raylac(ig_var_1281, js_var_1286) + z_fs_var_1293 * (raylac(ig_var_1281, js_var_1286 + 1) - raylac(ig_var_1281, js_var_1286)))
          p_taug_var_1278(iplon_var_1290, i_lay_var_1287, ig_var_1281) = z_speccomb_var_1294 * ((1.0D0 - z_fs_var_1293) * (absa_var_280(ind0_var_1282, ig_var_1281) * p_fac00_var_1258(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind0_var_1282 + 9, ig_var_1281) * p_fac10_var_1260(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind1_var_1283, ig_var_1281) * p_fac01_var_1259(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind1_var_1283 + 9, ig_var_1281) * p_fac11_var_1261(iplon_var_1290, i_lay_var_1287)) + z_fs_var_1293 * (absa_var_280(ind0_var_1282 + 1, ig_var_1281) * p_fac00_var_1258(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind0_var_1282 + 10, ig_var_1281) * p_fac10_var_1260(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind1_var_1283 + 1, ig_var_1281) * p_fac01_var_1259(iplon_var_1290, i_lay_var_1287) + absa_var_280(ind1_var_1283 + 10, ig_var_1281) * p_fac11_var_1261(iplon_var_1290, i_lay_var_1287))) + p_colo3_var_1269(iplon_var_1290, i_lay_var_1287) * abso3ac_var_285(ig_var_1281) + p_colh2o_var_1266(iplon_var_1290, i_lay_var_1287) * (p_selffac_var_1271(iplon_var_1290, i_lay_var_1287) * (selfrefc_var_282(inds_var_1284, ig_var_1281) + p_selffrac_var_1272(iplon_var_1290, i_lay_var_1287) * (selfrefc_var_282(inds_var_1284 + 1, ig_var_1281) - selfrefc_var_282(inds_var_1284, ig_var_1281))) + p_forfac_var_1274(iplon_var_1290, i_lay_var_1287) * (forrefc_var_283(indf_var_1285, ig_var_1281) + p_forfrac_var_1275(iplon_var_1290, i_lay_var_1287) * (forrefc_var_283(indf_var_1285 + 1, ig_var_1281) - forrefc_var_283(indf_var_1285, ig_var_1281))))
          IF (i_lay_var_1287 == i_laysolfr_var_1288(iplon_var_1290)) p_sfluxzen_var_1277(iplon_var_1290, ig_var_1281) = sfluxrefc_var_284(ig_var_1281, js_var_1286) + z_fs_var_1293 * (sfluxrefc_var_284(ig_var_1281, js_var_1286 + 1) - sfluxrefc_var_284(ig_var_1281, js_var_1286))
          p_taur_var_1279(iplon_var_1290, i_lay_var_1287, ig_var_1281) = z_tauray_var_1297
        END DO
      ELSE
        ind0_var_1282 = ((k_jp_var_1262(iplon_var_1290, i_lay_var_1287) - 13) * 5 + (k_jt_var_1263(iplon_var_1290, i_lay_var_1287) - 1)) * nspb_var_314(24) + 1
        ind1_var_1283 = ((k_jp_var_1262(iplon_var_1290, i_lay_var_1287) - 12) * 5 + (k_jt1_var_1264(iplon_var_1290, i_lay_var_1287) - 1)) * nspb_var_314(24) + 1
        DO ig_var_1281 = 1, 8
          z_tauray_var_1297 = p_colmol_var_1267(iplon_var_1290, i_lay_var_1287) * raylbc(ig_var_1281)
          p_taug_var_1278(iplon_var_1290, i_lay_var_1287, ig_var_1281) = p_colo2_var_1268(iplon_var_1290, i_lay_var_1287) * (p_fac00_var_1258(iplon_var_1290, i_lay_var_1287) * absb_var_281(ind0_var_1282, ig_var_1281) + p_fac10_var_1260(iplon_var_1290, i_lay_var_1287) * absb_var_281(ind0_var_1282 + 1, ig_var_1281) + p_fac01_var_1259(iplon_var_1290, i_lay_var_1287) * absb_var_281(ind1_var_1283, ig_var_1281) + p_fac11_var_1261(iplon_var_1290, i_lay_var_1287) * absb_var_281(ind1_var_1283 + 1, ig_var_1281)) + p_colo3_var_1269(iplon_var_1290, i_lay_var_1287) * abso3bc_var_286(ig_var_1281)
          p_taur_var_1279(iplon_var_1290, i_lay_var_1287, ig_var_1281) = z_tauray_var_1297
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1287 = laytrop_max_var_1292 + 1, i_nlayers_var_1289
    DO iplon_var_1290 = kidia_var_1255, kfdia_var_1256
      ind0_var_1282 = ((k_jp_var_1262(iplon_var_1290, i_lay_var_1287) - 13) * 5 + (k_jt_var_1263(iplon_var_1290, i_lay_var_1287) - 1)) * nspb_var_314(24) + 1
      ind1_var_1283 = ((k_jp_var_1262(iplon_var_1290, i_lay_var_1287) - 12) * 5 + (k_jt1_var_1264(iplon_var_1290, i_lay_var_1287) - 1)) * nspb_var_314(24) + 1
      DO ig_var_1281 = 1, 8
        z_tauray_var_1297 = p_colmol_var_1267(iplon_var_1290, i_lay_var_1287) * raylbc(ig_var_1281)
        p_taug_var_1278(iplon_var_1290, i_lay_var_1287, ig_var_1281) = p_colo2_var_1268(iplon_var_1290, i_lay_var_1287) * (p_fac00_var_1258(iplon_var_1290, i_lay_var_1287) * absb_var_281(ind0_var_1282, ig_var_1281) + p_fac10_var_1260(iplon_var_1290, i_lay_var_1287) * absb_var_281(ind0_var_1282 + 1, ig_var_1281) + p_fac01_var_1259(iplon_var_1290, i_lay_var_1287) * absb_var_281(ind1_var_1283, ig_var_1281) + p_fac11_var_1261(iplon_var_1290, i_lay_var_1287) * absb_var_281(ind1_var_1283 + 1, ig_var_1281)) + p_colo3_var_1269(iplon_var_1290, i_lay_var_1287) * abso3bc_var_286(ig_var_1281)
        p_taur_var_1279(iplon_var_1290, i_lay_var_1287, ig_var_1281) = z_tauray_var_1297
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol24
SUBROUTINE srtm_taumol25(kidia_var_1298, kfdia_var_1299, klev_var_1300, p_fac00_var_1301, p_fac01_var_1302, p_fac10_var_1303, p_fac11_var_1304, k_jp_var_1305, k_jt_var_1306, k_jt1_var_1307, p_colh2o_var_1308, p_colmol_var_1309, p_colo3_var_1310, k_laytrop_var_1311, p_sfluxzen_var_1312, p_taug_var_1313, p_taur_var_1314, prmu0_var_1315)
  USE yoesrta25, ONLY: absa_var_288, abso3ac_var_291, abso3bc_var_292, layreffr_var_287, raylc_var_290, sfluxrefc_var_289
  USE yoesrtwn, ONLY: nspa_var_313
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1298, kfdia_var_1299
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1300
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1301(kidia_var_1298 : kfdia_var_1299, klev_var_1300)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1302(kidia_var_1298 : kfdia_var_1299, klev_var_1300)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1303(kidia_var_1298 : kfdia_var_1299, klev_var_1300)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1304(kidia_var_1298 : kfdia_var_1299, klev_var_1300)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1305(kidia_var_1298 : kfdia_var_1299, klev_var_1300)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1306(kidia_var_1298 : kfdia_var_1299, klev_var_1300)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1307(kidia_var_1298 : kfdia_var_1299, klev_var_1300)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1308(kidia_var_1298 : kfdia_var_1299, klev_var_1300)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1309(kidia_var_1298 : kfdia_var_1299, klev_var_1300)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1310(kidia_var_1298 : kfdia_var_1299, klev_var_1300)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1311(kidia_var_1298 : kfdia_var_1299)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1312(kidia_var_1298 : kfdia_var_1299, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1313(kidia_var_1298 : kfdia_var_1299, klev_var_1300, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1314(kidia_var_1298 : kfdia_var_1299, klev_var_1300, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1315(kidia_var_1298 : kfdia_var_1299)
  INTEGER(KIND = 4) :: ig_var_1316, ind0_var_1317, ind1_var_1318, i_lay_var_1319, i_laysolfr_var_1320(kidia_var_1298 : kfdia_var_1299), i_nlayers_var_1321, iplon_var_1322
  INTEGER(KIND = 4) :: laytrop_min_var_1323, laytrop_max_var_1324
  REAL(KIND = 8) :: z_tauray_var_1325
  laytrop_min_var_1323 = MINVAL(k_laytrop_var_1311(kidia_var_1298 : kfdia_var_1299))
  laytrop_max_var_1324 = MAXVAL(k_laytrop_var_1311(kidia_var_1298 : kfdia_var_1299))
  i_nlayers_var_1321 = klev_var_1300
  DO iplon_var_1322 = kidia_var_1298, kfdia_var_1299
    i_laysolfr_var_1320(iplon_var_1322) = k_laytrop_var_1311(iplon_var_1322)
  END DO
  DO i_lay_var_1319 = 1, laytrop_min_var_1323
    DO iplon_var_1322 = kidia_var_1298, kfdia_var_1299
      IF (k_jp_var_1305(iplon_var_1322, i_lay_var_1319) < layreffr_var_287 .AND. k_jp_var_1305(iplon_var_1322, i_lay_var_1319 + 1) >= layreffr_var_287) i_laysolfr_var_1320(iplon_var_1322) = MIN(i_lay_var_1319 + 1, k_laytrop_var_1311(iplon_var_1322))
      ind0_var_1317 = ((k_jp_var_1305(iplon_var_1322, i_lay_var_1319) - 1) * 5 + (k_jt_var_1306(iplon_var_1322, i_lay_var_1319) - 1)) * nspa_var_313(25) + 1
      ind1_var_1318 = (k_jp_var_1305(iplon_var_1322, i_lay_var_1319) * 5 + (k_jt1_var_1307(iplon_var_1322, i_lay_var_1319) - 1)) * nspa_var_313(25) + 1
      DO ig_var_1316 = 1, 6
        z_tauray_var_1325 = p_colmol_var_1309(iplon_var_1322, i_lay_var_1319) * raylc_var_290(ig_var_1316)
        p_taug_var_1313(iplon_var_1322, i_lay_var_1319, ig_var_1316) = p_colh2o_var_1308(iplon_var_1322, i_lay_var_1319) * (p_fac00_var_1301(iplon_var_1322, i_lay_var_1319) * absa_var_288(ind0_var_1317, ig_var_1316) + p_fac10_var_1303(iplon_var_1322, i_lay_var_1319) * absa_var_288(ind0_var_1317 + 1, ig_var_1316) + p_fac01_var_1302(iplon_var_1322, i_lay_var_1319) * absa_var_288(ind1_var_1318, ig_var_1316) + p_fac11_var_1304(iplon_var_1322, i_lay_var_1319) * absa_var_288(ind1_var_1318 + 1, ig_var_1316)) + p_colo3_var_1310(iplon_var_1322, i_lay_var_1319) * abso3ac_var_291(ig_var_1316)
        IF (i_lay_var_1319 == i_laysolfr_var_1320(iplon_var_1322)) p_sfluxzen_var_1312(iplon_var_1322, ig_var_1316) = sfluxrefc_var_289(ig_var_1316)
        p_taur_var_1314(iplon_var_1322, i_lay_var_1319, ig_var_1316) = z_tauray_var_1325
      END DO
    END DO
  END DO
  DO i_lay_var_1319 = laytrop_min_var_1323 + 1, laytrop_max_var_1324
    DO iplon_var_1322 = kidia_var_1298, kfdia_var_1299
      IF (i_lay_var_1319 <= k_laytrop_var_1311(iplon_var_1322)) THEN
        IF (k_jp_var_1305(iplon_var_1322, i_lay_var_1319) < layreffr_var_287 .AND. k_jp_var_1305(iplon_var_1322, i_lay_var_1319 + 1) >= layreffr_var_287) i_laysolfr_var_1320(iplon_var_1322) = MIN(i_lay_var_1319 + 1, k_laytrop_var_1311(iplon_var_1322))
        ind0_var_1317 = ((k_jp_var_1305(iplon_var_1322, i_lay_var_1319) - 1) * 5 + (k_jt_var_1306(iplon_var_1322, i_lay_var_1319) - 1)) * nspa_var_313(25) + 1
        ind1_var_1318 = (k_jp_var_1305(iplon_var_1322, i_lay_var_1319) * 5 + (k_jt1_var_1307(iplon_var_1322, i_lay_var_1319) - 1)) * nspa_var_313(25) + 1
        DO ig_var_1316 = 1, 6
          z_tauray_var_1325 = p_colmol_var_1309(iplon_var_1322, i_lay_var_1319) * raylc_var_290(ig_var_1316)
          p_taug_var_1313(iplon_var_1322, i_lay_var_1319, ig_var_1316) = p_colh2o_var_1308(iplon_var_1322, i_lay_var_1319) * (p_fac00_var_1301(iplon_var_1322, i_lay_var_1319) * absa_var_288(ind0_var_1317, ig_var_1316) + p_fac10_var_1303(iplon_var_1322, i_lay_var_1319) * absa_var_288(ind0_var_1317 + 1, ig_var_1316) + p_fac01_var_1302(iplon_var_1322, i_lay_var_1319) * absa_var_288(ind1_var_1318, ig_var_1316) + p_fac11_var_1304(iplon_var_1322, i_lay_var_1319) * absa_var_288(ind1_var_1318 + 1, ig_var_1316)) + p_colo3_var_1310(iplon_var_1322, i_lay_var_1319) * abso3ac_var_291(ig_var_1316)
          IF (i_lay_var_1319 == i_laysolfr_var_1320(iplon_var_1322)) p_sfluxzen_var_1312(iplon_var_1322, ig_var_1316) = sfluxrefc_var_289(ig_var_1316)
          p_taur_var_1314(iplon_var_1322, i_lay_var_1319, ig_var_1316) = z_tauray_var_1325
        END DO
      ELSE
        DO ig_var_1316 = 1, 6
          z_tauray_var_1325 = p_colmol_var_1309(iplon_var_1322, i_lay_var_1319) * raylc_var_290(ig_var_1316)
          p_taug_var_1313(iplon_var_1322, i_lay_var_1319, ig_var_1316) = p_colo3_var_1310(iplon_var_1322, i_lay_var_1319) * abso3bc_var_292(ig_var_1316)
          p_taur_var_1314(iplon_var_1322, i_lay_var_1319, ig_var_1316) = z_tauray_var_1325
        END DO
      END IF
    END DO
  END DO
  DO ig_var_1316 = 1, 6
    DO i_lay_var_1319 = laytrop_max_var_1324 + 1, i_nlayers_var_1321
      DO iplon_var_1322 = kidia_var_1298, kfdia_var_1299
        z_tauray_var_1325 = p_colmol_var_1309(iplon_var_1322, i_lay_var_1319) * raylc_var_290(ig_var_1316)
        p_taug_var_1313(iplon_var_1322, i_lay_var_1319, ig_var_1316) = p_colo3_var_1310(iplon_var_1322, i_lay_var_1319) * abso3bc_var_292(ig_var_1316)
        p_taur_var_1314(iplon_var_1322, i_lay_var_1319, ig_var_1316) = z_tauray_var_1325
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol25
SUBROUTINE srtm_taumol19(kidia_var_1326, kfdia_var_1327, klev_var_1328, p_fac00_var_1329, p_fac01_var_1330, p_fac10_var_1331, p_fac11_var_1332, k_jp_var_1333, k_jt_var_1334, k_jt1_var_1335, p_oneminus_var_1336, p_colh2o_var_1337, p_colco2_var_1338, p_colmol_var_1339, k_laytrop_var_1340, p_selffac_var_1341, p_selffrac_var_1342, k_indself_var_1343, p_forfac_var_1344, p_forfrac_var_1345, k_indfor_var_1346, p_sfluxzen_var_1347, p_taug_var_1348, p_taur_var_1349, prmu0_var_1350)
  USE yoesrta19, ONLY: absa_var_244, absb_var_245, forrefc_var_247, layreffr_var_243, rayl_var_241, selfrefc_var_246, sfluxrefc_var_248, strrat_var_242
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1326, kfdia_var_1327
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1328
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1329(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1330(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1331(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1332(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1333(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1334(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1335(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1336(kidia_var_1326 : kfdia_var_1327)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1337(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1338(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1339(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1340(kidia_var_1326 : kfdia_var_1327)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1341(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1342(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1343(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1344(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1345(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1346(kidia_var_1326 : kfdia_var_1327, klev_var_1328)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1347(kidia_var_1326 : kfdia_var_1327, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1348(kidia_var_1326 : kfdia_var_1327, klev_var_1328, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1349(kidia_var_1326 : kfdia_var_1327, klev_var_1328, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1350(kidia_var_1326 : kfdia_var_1327)
  INTEGER(KIND = 4) :: ig_var_1351, ind0_var_1352, ind1_var_1353, inds_var_1354, indf_var_1355, js_var_1356, i_lay_var_1357, i_laysolfr_var_1358(kidia_var_1326 : kfdia_var_1327), i_nlayers_var_1359, iplon_var_1360
  INTEGER(KIND = 4) :: laytrop_min_var_1361, laytrop_max_var_1362
  REAL(KIND = 8) :: z_fs_var_1363, z_speccomb_var_1364, z_specmult_var_1365, z_specparm_var_1366, z_tauray_var_1367
  laytrop_min_var_1361 = MINVAL(k_laytrop_var_1340(kidia_var_1326 : kfdia_var_1327))
  laytrop_max_var_1362 = MAXVAL(k_laytrop_var_1340(kidia_var_1326 : kfdia_var_1327))
  i_nlayers_var_1359 = klev_var_1328
  DO iplon_var_1360 = kidia_var_1326, kfdia_var_1327
    i_laysolfr_var_1358(iplon_var_1360) = k_laytrop_var_1340(iplon_var_1360)
  END DO
  DO i_lay_var_1357 = 1, laytrop_min_var_1361
    DO iplon_var_1360 = kidia_var_1326, kfdia_var_1327
      IF (k_jp_var_1333(iplon_var_1360, i_lay_var_1357) < layreffr_var_243 .AND. k_jp_var_1333(iplon_var_1360, i_lay_var_1357 + 1) >= layreffr_var_243) i_laysolfr_var_1358(iplon_var_1360) = MIN(i_lay_var_1357 + 1, k_laytrop_var_1340(iplon_var_1360))
      z_speccomb_var_1364 = p_colh2o_var_1337(iplon_var_1360, i_lay_var_1357) + strrat_var_242 * p_colco2_var_1338(iplon_var_1360, i_lay_var_1357)
      z_specparm_var_1366 = p_colh2o_var_1337(iplon_var_1360, i_lay_var_1357) / z_speccomb_var_1364
      z_specparm_var_1366 = MIN(p_oneminus_var_1336(iplon_var_1360), z_specparm_var_1366)
      z_specmult_var_1365 = 8.0D0 * (z_specparm_var_1366)
      js_var_1356 = 1 + INT(z_specmult_var_1365)
      z_fs_var_1363 = z_specmult_var_1365 - AINT(z_specmult_var_1365)
      ind0_var_1352 = ((k_jp_var_1333(iplon_var_1360, i_lay_var_1357) - 1) * 5 + (k_jt_var_1334(iplon_var_1360, i_lay_var_1357) - 1)) * nspa_var_313(19) + js_var_1356
      ind1_var_1353 = (k_jp_var_1333(iplon_var_1360, i_lay_var_1357) * 5 + (k_jt1_var_1335(iplon_var_1360, i_lay_var_1357) - 1)) * nspa_var_313(19) + js_var_1356
      inds_var_1354 = k_indself_var_1343(iplon_var_1360, i_lay_var_1357)
      indf_var_1355 = k_indfor_var_1346(iplon_var_1360, i_lay_var_1357)
      z_tauray_var_1367 = p_colmol_var_1339(iplon_var_1360, i_lay_var_1357) * rayl_var_241
      DO ig_var_1351 = 1, 8
        p_taug_var_1348(iplon_var_1360, i_lay_var_1357, ig_var_1351) = z_speccomb_var_1364 * ((1.0D0 - z_fs_var_1363) * (absa_var_244(ind0_var_1352, ig_var_1351) * p_fac00_var_1329(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind0_var_1352 + 9, ig_var_1351) * p_fac10_var_1331(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind1_var_1353, ig_var_1351) * p_fac01_var_1330(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind1_var_1353 + 9, ig_var_1351) * p_fac11_var_1332(iplon_var_1360, i_lay_var_1357)) + z_fs_var_1363 * (absa_var_244(ind0_var_1352 + 1, ig_var_1351) * p_fac00_var_1329(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind0_var_1352 + 10, ig_var_1351) * p_fac10_var_1331(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind1_var_1353 + 1, ig_var_1351) * p_fac01_var_1330(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind1_var_1353 + 10, ig_var_1351) * p_fac11_var_1332(iplon_var_1360, i_lay_var_1357))) + p_colh2o_var_1337(iplon_var_1360, i_lay_var_1357) * (p_selffac_var_1341(iplon_var_1360, i_lay_var_1357) * (selfrefc_var_246(inds_var_1354, ig_var_1351) + p_selffrac_var_1342(iplon_var_1360, i_lay_var_1357) * (selfrefc_var_246(inds_var_1354 + 1, ig_var_1351) - selfrefc_var_246(inds_var_1354, ig_var_1351))) + p_forfac_var_1344(iplon_var_1360, i_lay_var_1357) * (forrefc_var_247(indf_var_1355, ig_var_1351) + p_forfrac_var_1345(iplon_var_1360, i_lay_var_1357) * (forrefc_var_247(indf_var_1355 + 1, ig_var_1351) - forrefc_var_247(indf_var_1355, ig_var_1351))))
        IF (i_lay_var_1357 == i_laysolfr_var_1358(iplon_var_1360)) p_sfluxzen_var_1347(iplon_var_1360, ig_var_1351) = sfluxrefc_var_248(ig_var_1351, js_var_1356) + z_fs_var_1363 * (sfluxrefc_var_248(ig_var_1351, js_var_1356 + 1) - sfluxrefc_var_248(ig_var_1351, js_var_1356))
        p_taur_var_1349(iplon_var_1360, i_lay_var_1357, ig_var_1351) = z_tauray_var_1367
      END DO
    END DO
  END DO
  DO i_lay_var_1357 = laytrop_min_var_1361 + 1, laytrop_max_var_1362
    DO iplon_var_1360 = kidia_var_1326, kfdia_var_1327
      IF (i_lay_var_1357 <= k_laytrop_var_1340(iplon_var_1360)) THEN
        IF (k_jp_var_1333(iplon_var_1360, i_lay_var_1357) < layreffr_var_243 .AND. k_jp_var_1333(iplon_var_1360, i_lay_var_1357 + 1) >= layreffr_var_243) i_laysolfr_var_1358(iplon_var_1360) = MIN(i_lay_var_1357 + 1, k_laytrop_var_1340(iplon_var_1360))
        z_speccomb_var_1364 = p_colh2o_var_1337(iplon_var_1360, i_lay_var_1357) + strrat_var_242 * p_colco2_var_1338(iplon_var_1360, i_lay_var_1357)
        z_specparm_var_1366 = p_colh2o_var_1337(iplon_var_1360, i_lay_var_1357) / z_speccomb_var_1364
        z_specparm_var_1366 = MIN(p_oneminus_var_1336(iplon_var_1360), z_specparm_var_1366)
        z_specmult_var_1365 = 8.0D0 * (z_specparm_var_1366)
        js_var_1356 = 1 + INT(z_specmult_var_1365)
        z_fs_var_1363 = z_specmult_var_1365 - AINT(z_specmult_var_1365)
        ind0_var_1352 = ((k_jp_var_1333(iplon_var_1360, i_lay_var_1357) - 1) * 5 + (k_jt_var_1334(iplon_var_1360, i_lay_var_1357) - 1)) * nspa_var_313(19) + js_var_1356
        ind1_var_1353 = (k_jp_var_1333(iplon_var_1360, i_lay_var_1357) * 5 + (k_jt1_var_1335(iplon_var_1360, i_lay_var_1357) - 1)) * nspa_var_313(19) + js_var_1356
        inds_var_1354 = k_indself_var_1343(iplon_var_1360, i_lay_var_1357)
        indf_var_1355 = k_indfor_var_1346(iplon_var_1360, i_lay_var_1357)
        z_tauray_var_1367 = p_colmol_var_1339(iplon_var_1360, i_lay_var_1357) * rayl_var_241
        DO ig_var_1351 = 1, 8
          p_taug_var_1348(iplon_var_1360, i_lay_var_1357, ig_var_1351) = z_speccomb_var_1364 * ((1.0D0 - z_fs_var_1363) * (absa_var_244(ind0_var_1352, ig_var_1351) * p_fac00_var_1329(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind0_var_1352 + 9, ig_var_1351) * p_fac10_var_1331(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind1_var_1353, ig_var_1351) * p_fac01_var_1330(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind1_var_1353 + 9, ig_var_1351) * p_fac11_var_1332(iplon_var_1360, i_lay_var_1357)) + z_fs_var_1363 * (absa_var_244(ind0_var_1352 + 1, ig_var_1351) * p_fac00_var_1329(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind0_var_1352 + 10, ig_var_1351) * p_fac10_var_1331(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind1_var_1353 + 1, ig_var_1351) * p_fac01_var_1330(iplon_var_1360, i_lay_var_1357) + absa_var_244(ind1_var_1353 + 10, ig_var_1351) * p_fac11_var_1332(iplon_var_1360, i_lay_var_1357))) + p_colh2o_var_1337(iplon_var_1360, i_lay_var_1357) * (p_selffac_var_1341(iplon_var_1360, i_lay_var_1357) * (selfrefc_var_246(inds_var_1354, ig_var_1351) + p_selffrac_var_1342(iplon_var_1360, i_lay_var_1357) * (selfrefc_var_246(inds_var_1354 + 1, ig_var_1351) - selfrefc_var_246(inds_var_1354, ig_var_1351))) + p_forfac_var_1344(iplon_var_1360, i_lay_var_1357) * (forrefc_var_247(indf_var_1355, ig_var_1351) + p_forfrac_var_1345(iplon_var_1360, i_lay_var_1357) * (forrefc_var_247(indf_var_1355 + 1, ig_var_1351) - forrefc_var_247(indf_var_1355, ig_var_1351))))
          IF (i_lay_var_1357 == i_laysolfr_var_1358(iplon_var_1360)) p_sfluxzen_var_1347(iplon_var_1360, ig_var_1351) = sfluxrefc_var_248(ig_var_1351, js_var_1356) + z_fs_var_1363 * (sfluxrefc_var_248(ig_var_1351, js_var_1356 + 1) - sfluxrefc_var_248(ig_var_1351, js_var_1356))
          p_taur_var_1349(iplon_var_1360, i_lay_var_1357, ig_var_1351) = z_tauray_var_1367
        END DO
      ELSE
        ind0_var_1352 = ((k_jp_var_1333(iplon_var_1360, i_lay_var_1357) - 13) * 5 + (k_jt_var_1334(iplon_var_1360, i_lay_var_1357) - 1)) * nspb_var_314(19) + 1
        ind1_var_1353 = ((k_jp_var_1333(iplon_var_1360, i_lay_var_1357) - 12) * 5 + (k_jt1_var_1335(iplon_var_1360, i_lay_var_1357) - 1)) * nspb_var_314(19) + 1
        z_tauray_var_1367 = p_colmol_var_1339(iplon_var_1360, i_lay_var_1357) * rayl_var_241
        DO ig_var_1351 = 1, 8
          p_taug_var_1348(iplon_var_1360, i_lay_var_1357, ig_var_1351) = p_colco2_var_1338(iplon_var_1360, i_lay_var_1357) * (p_fac00_var_1329(iplon_var_1360, i_lay_var_1357) * absb_var_245(ind0_var_1352, ig_var_1351) + p_fac10_var_1331(iplon_var_1360, i_lay_var_1357) * absb_var_245(ind0_var_1352 + 1, ig_var_1351) + p_fac01_var_1330(iplon_var_1360, i_lay_var_1357) * absb_var_245(ind1_var_1353, ig_var_1351) + p_fac11_var_1332(iplon_var_1360, i_lay_var_1357) * absb_var_245(ind1_var_1353 + 1, ig_var_1351))
          p_taur_var_1349(iplon_var_1360, i_lay_var_1357, ig_var_1351) = z_tauray_var_1367
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1357 = laytrop_max_var_1362 + 1, i_nlayers_var_1359
    DO iplon_var_1360 = kidia_var_1326, kfdia_var_1327
      ind0_var_1352 = ((k_jp_var_1333(iplon_var_1360, i_lay_var_1357) - 13) * 5 + (k_jt_var_1334(iplon_var_1360, i_lay_var_1357) - 1)) * nspb_var_314(19) + 1
      ind1_var_1353 = ((k_jp_var_1333(iplon_var_1360, i_lay_var_1357) - 12) * 5 + (k_jt1_var_1335(iplon_var_1360, i_lay_var_1357) - 1)) * nspb_var_314(19) + 1
      z_tauray_var_1367 = p_colmol_var_1339(iplon_var_1360, i_lay_var_1357) * rayl_var_241
      DO ig_var_1351 = 1, 8
        p_taug_var_1348(iplon_var_1360, i_lay_var_1357, ig_var_1351) = p_colco2_var_1338(iplon_var_1360, i_lay_var_1357) * (p_fac00_var_1329(iplon_var_1360, i_lay_var_1357) * absb_var_245(ind0_var_1352, ig_var_1351) + p_fac10_var_1331(iplon_var_1360, i_lay_var_1357) * absb_var_245(ind0_var_1352 + 1, ig_var_1351) + p_fac01_var_1330(iplon_var_1360, i_lay_var_1357) * absb_var_245(ind1_var_1353, ig_var_1351) + p_fac11_var_1332(iplon_var_1360, i_lay_var_1357) * absb_var_245(ind1_var_1353 + 1, ig_var_1351))
        p_taur_var_1349(iplon_var_1360, i_lay_var_1357, ig_var_1351) = z_tauray_var_1367
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol19
SUBROUTINE srtm_taumol21(kidia_var_1368, kfdia_var_1369, klev_var_1370, p_fac00_var_1371, p_fac01_var_1372, p_fac10_var_1373, p_fac11_var_1374, k_jp_var_1375, k_jt_var_1376, k_jt1_var_1377, p_oneminus_var_1378, p_colh2o_var_1379, p_colco2_var_1380, p_colmol_var_1381, k_laytrop_var_1382, p_selffac_var_1383, p_selffrac_var_1384, k_indself_var_1385, p_forfac_var_1386, p_forfrac_var_1387, k_indfor_var_1388, p_sfluxzen_var_1389, p_taug_var_1390, p_taur_var_1391, prmu0_var_1392)
  USE yoesrta21, ONLY: absa_var_259, absb_var_260, forrefc_var_262, layreffr_var_258, rayl_var_256, selfrefc_var_261, sfluxrefc_var_263, strrat_var_257
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1368, kfdia_var_1369
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1370
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1371(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1372(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1373(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1374(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1375(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1376(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1377(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1378(kidia_var_1368 : kfdia_var_1369)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1379(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1380(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1381(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1382(kidia_var_1368 : kfdia_var_1369)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1383(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1384(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1385(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1386(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1387(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1388(kidia_var_1368 : kfdia_var_1369, klev_var_1370)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1389(kidia_var_1368 : kfdia_var_1369, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1390(kidia_var_1368 : kfdia_var_1369, klev_var_1370, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1391(kidia_var_1368 : kfdia_var_1369, klev_var_1370, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1392(kidia_var_1368 : kfdia_var_1369)
  INTEGER(KIND = 4) :: ig_var_1393, ind0_var_1394, ind1_var_1395, inds_var_1396, indf_var_1397, js_var_1398, i_lay_var_1399, i_laysolfr_var_1400(kidia_var_1368 : kfdia_var_1369), i_nlayers_var_1401, iplon_var_1402
  INTEGER(KIND = 4) :: laytrop_min_var_1403, laytrop_max_var_1404
  REAL(KIND = 8) :: z_fs_var_1405, z_speccomb_var_1406, z_specmult_var_1407, z_specparm_var_1408, z_tauray_var_1409
  laytrop_min_var_1403 = MINVAL(k_laytrop_var_1382(kidia_var_1368 : kfdia_var_1369))
  laytrop_max_var_1404 = MAXVAL(k_laytrop_var_1382(kidia_var_1368 : kfdia_var_1369))
  i_nlayers_var_1401 = klev_var_1370
  DO iplon_var_1402 = kidia_var_1368, kfdia_var_1369
    i_laysolfr_var_1400(iplon_var_1402) = k_laytrop_var_1382(iplon_var_1402)
  END DO
  DO i_lay_var_1399 = 1, laytrop_min_var_1403
    DO iplon_var_1402 = kidia_var_1368, kfdia_var_1369
      IF (k_jp_var_1375(iplon_var_1402, i_lay_var_1399) < layreffr_var_258 .AND. k_jp_var_1375(iplon_var_1402, i_lay_var_1399 + 1) >= layreffr_var_258) i_laysolfr_var_1400(iplon_var_1402) = MIN(i_lay_var_1399 + 1, k_laytrop_var_1382(iplon_var_1402))
      z_speccomb_var_1406 = p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) + strrat_var_257 * p_colco2_var_1380(iplon_var_1402, i_lay_var_1399)
      z_specparm_var_1408 = p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) / z_speccomb_var_1406
      z_specparm_var_1408 = MIN(p_oneminus_var_1378(iplon_var_1402), z_specparm_var_1408)
      z_specmult_var_1407 = 8.0D0 * (z_specparm_var_1408)
      js_var_1398 = 1 + INT(z_specmult_var_1407)
      z_fs_var_1405 = z_specmult_var_1407 - AINT(z_specmult_var_1407)
      ind0_var_1394 = ((k_jp_var_1375(iplon_var_1402, i_lay_var_1399) - 1) * 5 + (k_jt_var_1376(iplon_var_1402, i_lay_var_1399) - 1)) * nspa_var_313(21) + js_var_1398
      ind1_var_1395 = (k_jp_var_1375(iplon_var_1402, i_lay_var_1399) * 5 + (k_jt1_var_1377(iplon_var_1402, i_lay_var_1399) - 1)) * nspa_var_313(21) + js_var_1398
      inds_var_1396 = k_indself_var_1385(iplon_var_1402, i_lay_var_1399)
      indf_var_1397 = k_indfor_var_1388(iplon_var_1402, i_lay_var_1399)
      z_tauray_var_1409 = p_colmol_var_1381(iplon_var_1402, i_lay_var_1399) * rayl_var_256
      DO ig_var_1393 = 1, 10
        p_taug_var_1390(iplon_var_1402, i_lay_var_1399, ig_var_1393) = z_speccomb_var_1406 * ((1.0D0 - z_fs_var_1405) * (absa_var_259(ind0_var_1394, ig_var_1393) * p_fac00_var_1371(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind0_var_1394 + 9, ig_var_1393) * p_fac10_var_1373(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind1_var_1395, ig_var_1393) * p_fac01_var_1372(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind1_var_1395 + 9, ig_var_1393) * p_fac11_var_1374(iplon_var_1402, i_lay_var_1399)) + z_fs_var_1405 * (absa_var_259(ind0_var_1394 + 1, ig_var_1393) * p_fac00_var_1371(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind0_var_1394 + 10, ig_var_1393) * p_fac10_var_1373(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind1_var_1395 + 1, ig_var_1393) * p_fac01_var_1372(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind1_var_1395 + 10, ig_var_1393) * p_fac11_var_1374(iplon_var_1402, i_lay_var_1399))) + p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) * (p_selffac_var_1383(iplon_var_1402, i_lay_var_1399) * (selfrefc_var_261(inds_var_1396, ig_var_1393) + p_selffrac_var_1384(iplon_var_1402, i_lay_var_1399) * (selfrefc_var_261(inds_var_1396 + 1, ig_var_1393) - selfrefc_var_261(inds_var_1396, ig_var_1393))) + p_forfac_var_1386(iplon_var_1402, i_lay_var_1399) * (forrefc_var_262(indf_var_1397, ig_var_1393) + p_forfrac_var_1387(iplon_var_1402, i_lay_var_1399) * (forrefc_var_262(indf_var_1397 + 1, ig_var_1393) - forrefc_var_262(indf_var_1397, ig_var_1393))))
        IF (i_lay_var_1399 == i_laysolfr_var_1400(iplon_var_1402)) p_sfluxzen_var_1389(iplon_var_1402, ig_var_1393) = sfluxrefc_var_263(ig_var_1393, js_var_1398) + z_fs_var_1405 * (sfluxrefc_var_263(ig_var_1393, js_var_1398 + 1) - sfluxrefc_var_263(ig_var_1393, js_var_1398))
        p_taur_var_1391(iplon_var_1402, i_lay_var_1399, ig_var_1393) = z_tauray_var_1409
      END DO
    END DO
  END DO
  DO i_lay_var_1399 = laytrop_min_var_1403 + 1, laytrop_max_var_1404
    DO iplon_var_1402 = kidia_var_1368, kfdia_var_1369
      IF (i_lay_var_1399 <= k_laytrop_var_1382(iplon_var_1402)) THEN
        IF (k_jp_var_1375(iplon_var_1402, i_lay_var_1399) < layreffr_var_258 .AND. k_jp_var_1375(iplon_var_1402, i_lay_var_1399 + 1) >= layreffr_var_258) i_laysolfr_var_1400(iplon_var_1402) = MIN(i_lay_var_1399 + 1, k_laytrop_var_1382(iplon_var_1402))
        z_speccomb_var_1406 = p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) + strrat_var_257 * p_colco2_var_1380(iplon_var_1402, i_lay_var_1399)
        z_specparm_var_1408 = p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) / z_speccomb_var_1406
        z_specparm_var_1408 = MIN(p_oneminus_var_1378(iplon_var_1402), z_specparm_var_1408)
        z_specmult_var_1407 = 8.0D0 * (z_specparm_var_1408)
        js_var_1398 = 1 + INT(z_specmult_var_1407)
        z_fs_var_1405 = z_specmult_var_1407 - AINT(z_specmult_var_1407)
        ind0_var_1394 = ((k_jp_var_1375(iplon_var_1402, i_lay_var_1399) - 1) * 5 + (k_jt_var_1376(iplon_var_1402, i_lay_var_1399) - 1)) * nspa_var_313(21) + js_var_1398
        ind1_var_1395 = (k_jp_var_1375(iplon_var_1402, i_lay_var_1399) * 5 + (k_jt1_var_1377(iplon_var_1402, i_lay_var_1399) - 1)) * nspa_var_313(21) + js_var_1398
        inds_var_1396 = k_indself_var_1385(iplon_var_1402, i_lay_var_1399)
        indf_var_1397 = k_indfor_var_1388(iplon_var_1402, i_lay_var_1399)
        z_tauray_var_1409 = p_colmol_var_1381(iplon_var_1402, i_lay_var_1399) * rayl_var_256
        DO ig_var_1393 = 1, 10
          p_taug_var_1390(iplon_var_1402, i_lay_var_1399, ig_var_1393) = z_speccomb_var_1406 * ((1.0D0 - z_fs_var_1405) * (absa_var_259(ind0_var_1394, ig_var_1393) * p_fac00_var_1371(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind0_var_1394 + 9, ig_var_1393) * p_fac10_var_1373(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind1_var_1395, ig_var_1393) * p_fac01_var_1372(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind1_var_1395 + 9, ig_var_1393) * p_fac11_var_1374(iplon_var_1402, i_lay_var_1399)) + z_fs_var_1405 * (absa_var_259(ind0_var_1394 + 1, ig_var_1393) * p_fac00_var_1371(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind0_var_1394 + 10, ig_var_1393) * p_fac10_var_1373(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind1_var_1395 + 1, ig_var_1393) * p_fac01_var_1372(iplon_var_1402, i_lay_var_1399) + absa_var_259(ind1_var_1395 + 10, ig_var_1393) * p_fac11_var_1374(iplon_var_1402, i_lay_var_1399))) + p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) * (p_selffac_var_1383(iplon_var_1402, i_lay_var_1399) * (selfrefc_var_261(inds_var_1396, ig_var_1393) + p_selffrac_var_1384(iplon_var_1402, i_lay_var_1399) * (selfrefc_var_261(inds_var_1396 + 1, ig_var_1393) - selfrefc_var_261(inds_var_1396, ig_var_1393))) + p_forfac_var_1386(iplon_var_1402, i_lay_var_1399) * (forrefc_var_262(indf_var_1397, ig_var_1393) + p_forfrac_var_1387(iplon_var_1402, i_lay_var_1399) * (forrefc_var_262(indf_var_1397 + 1, ig_var_1393) - forrefc_var_262(indf_var_1397, ig_var_1393))))
          IF (i_lay_var_1399 == i_laysolfr_var_1400(iplon_var_1402)) p_sfluxzen_var_1389(iplon_var_1402, ig_var_1393) = sfluxrefc_var_263(ig_var_1393, js_var_1398) + z_fs_var_1405 * (sfluxrefc_var_263(ig_var_1393, js_var_1398 + 1) - sfluxrefc_var_263(ig_var_1393, js_var_1398))
          p_taur_var_1391(iplon_var_1402, i_lay_var_1399, ig_var_1393) = z_tauray_var_1409
        END DO
      ELSE
        z_speccomb_var_1406 = p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) + strrat_var_257 * p_colco2_var_1380(iplon_var_1402, i_lay_var_1399)
        z_specparm_var_1408 = p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) / z_speccomb_var_1406
        z_specparm_var_1408 = MIN(p_oneminus_var_1378(iplon_var_1402), z_specparm_var_1408)
        z_specmult_var_1407 = 4.0D0 * (z_specparm_var_1408)
        js_var_1398 = 1 + INT(z_specmult_var_1407)
        z_fs_var_1405 = z_specmult_var_1407 - AINT(z_specmult_var_1407)
        ind0_var_1394 = ((k_jp_var_1375(iplon_var_1402, i_lay_var_1399) - 13) * 5 + (k_jt_var_1376(iplon_var_1402, i_lay_var_1399) - 1)) * nspb_var_314(21) + js_var_1398
        ind1_var_1395 = ((k_jp_var_1375(iplon_var_1402, i_lay_var_1399) - 12) * 5 + (k_jt1_var_1377(iplon_var_1402, i_lay_var_1399) - 1)) * nspb_var_314(21) + js_var_1398
        indf_var_1397 = k_indfor_var_1388(iplon_var_1402, i_lay_var_1399)
        z_tauray_var_1409 = p_colmol_var_1381(iplon_var_1402, i_lay_var_1399) * rayl_var_256
        DO ig_var_1393 = 1, 10
          p_taug_var_1390(iplon_var_1402, i_lay_var_1399, ig_var_1393) = z_speccomb_var_1406 * ((1.0D0 - z_fs_var_1405) * (absb_var_260(ind0_var_1394, ig_var_1393) * p_fac00_var_1371(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind0_var_1394 + 5, ig_var_1393) * p_fac10_var_1373(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind1_var_1395, ig_var_1393) * p_fac01_var_1372(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind1_var_1395 + 5, ig_var_1393) * p_fac11_var_1374(iplon_var_1402, i_lay_var_1399)) + z_fs_var_1405 * (absb_var_260(ind0_var_1394 + 1, ig_var_1393) * p_fac00_var_1371(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind0_var_1394 + 6, ig_var_1393) * p_fac10_var_1373(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind1_var_1395 + 1, ig_var_1393) * p_fac01_var_1372(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind1_var_1395 + 6, ig_var_1393) * p_fac11_var_1374(iplon_var_1402, i_lay_var_1399))) + p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) * p_forfac_var_1386(iplon_var_1402, i_lay_var_1399) * (forrefc_var_262(indf_var_1397, ig_var_1393) + p_forfrac_var_1387(iplon_var_1402, i_lay_var_1399) * (forrefc_var_262(indf_var_1397 + 1, ig_var_1393) - forrefc_var_262(indf_var_1397, ig_var_1393)))
          p_taur_var_1391(iplon_var_1402, i_lay_var_1399, ig_var_1393) = z_tauray_var_1409
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1399 = laytrop_max_var_1404 + 1, i_nlayers_var_1401
    DO iplon_var_1402 = kidia_var_1368, kfdia_var_1369
      z_speccomb_var_1406 = p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) + strrat_var_257 * p_colco2_var_1380(iplon_var_1402, i_lay_var_1399)
      z_specparm_var_1408 = p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) / z_speccomb_var_1406
      z_specparm_var_1408 = MIN(p_oneminus_var_1378(iplon_var_1402), z_specparm_var_1408)
      z_specmult_var_1407 = 4.0D0 * (z_specparm_var_1408)
      js_var_1398 = 1 + INT(z_specmult_var_1407)
      z_fs_var_1405 = z_specmult_var_1407 - AINT(z_specmult_var_1407)
      ind0_var_1394 = ((k_jp_var_1375(iplon_var_1402, i_lay_var_1399) - 13) * 5 + (k_jt_var_1376(iplon_var_1402, i_lay_var_1399) - 1)) * nspb_var_314(21) + js_var_1398
      ind1_var_1395 = ((k_jp_var_1375(iplon_var_1402, i_lay_var_1399) - 12) * 5 + (k_jt1_var_1377(iplon_var_1402, i_lay_var_1399) - 1)) * nspb_var_314(21) + js_var_1398
      indf_var_1397 = k_indfor_var_1388(iplon_var_1402, i_lay_var_1399)
      z_tauray_var_1409 = p_colmol_var_1381(iplon_var_1402, i_lay_var_1399) * rayl_var_256
      DO ig_var_1393 = 1, 10
        p_taug_var_1390(iplon_var_1402, i_lay_var_1399, ig_var_1393) = z_speccomb_var_1406 * ((1.0D0 - z_fs_var_1405) * (absb_var_260(ind0_var_1394, ig_var_1393) * p_fac00_var_1371(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind0_var_1394 + 5, ig_var_1393) * p_fac10_var_1373(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind1_var_1395, ig_var_1393) * p_fac01_var_1372(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind1_var_1395 + 5, ig_var_1393) * p_fac11_var_1374(iplon_var_1402, i_lay_var_1399)) + z_fs_var_1405 * (absb_var_260(ind0_var_1394 + 1, ig_var_1393) * p_fac00_var_1371(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind0_var_1394 + 6, ig_var_1393) * p_fac10_var_1373(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind1_var_1395 + 1, ig_var_1393) * p_fac01_var_1372(iplon_var_1402, i_lay_var_1399) + absb_var_260(ind1_var_1395 + 6, ig_var_1393) * p_fac11_var_1374(iplon_var_1402, i_lay_var_1399))) + p_colh2o_var_1379(iplon_var_1402, i_lay_var_1399) * p_forfac_var_1386(iplon_var_1402, i_lay_var_1399) * (forrefc_var_262(indf_var_1397, ig_var_1393) + p_forfrac_var_1387(iplon_var_1402, i_lay_var_1399) * (forrefc_var_262(indf_var_1397 + 1, ig_var_1393) - forrefc_var_262(indf_var_1397, ig_var_1393)))
        p_taur_var_1391(iplon_var_1402, i_lay_var_1399, ig_var_1393) = z_tauray_var_1409
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol21
SUBROUTINE srtm_taumol20(kidia_var_1410, kfdia_var_1411, klev_var_1412, p_fac00_var_1413, p_fac01_var_1414, p_fac10_var_1415, p_fac11_var_1416, k_jp_var_1417, k_jt_var_1418, k_jt1_var_1419, p_colh2o_var_1420, p_colch4_var_1421, p_colmol_var_1422, k_laytrop_var_1423, p_selffac_var_1424, p_selffrac_var_1425, k_indself_var_1426, p_forfac_var_1427, p_forfrac_var_1428, k_indfor_var_1429, p_sfluxzen_var_1430, p_taug_var_1431, p_taur_var_1432, prmu0_var_1433)
  USE yoesrta20, ONLY: absa_var_251, absb_var_252, absch4c, forrefc_var_254, layreffr_var_250, rayl_var_249, selfrefc_var_253, sfluxrefc_var_255
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1410, kfdia_var_1411
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1412
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1413(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1414(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1415(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1416(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1417(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1418(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1419(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1420(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_1421(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1422(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1423(kidia_var_1410 : kfdia_var_1411)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1424(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1425(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1426(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1427(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1428(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1429(kidia_var_1410 : kfdia_var_1411, klev_var_1412)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1430(kidia_var_1410 : kfdia_var_1411, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1431(kidia_var_1410 : kfdia_var_1411, klev_var_1412, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1432(kidia_var_1410 : kfdia_var_1411, klev_var_1412, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1433(kidia_var_1410 : kfdia_var_1411)
  INTEGER(KIND = 4) :: ig_var_1434, ind0_var_1435, ind1_var_1436, inds_var_1437, indf_var_1438, i_lay_var_1439, i_laysolfr_var_1440(kidia_var_1410 : kfdia_var_1411), i_nlayers_var_1441, iplon_var_1442
  INTEGER(KIND = 4) :: laytrop_min_var_1443, laytrop_max_var_1444
  REAL(KIND = 8) :: z_tauray_var_1445
  laytrop_min_var_1443 = MINVAL(k_laytrop_var_1423(kidia_var_1410 : kfdia_var_1411))
  laytrop_max_var_1444 = MAXVAL(k_laytrop_var_1423(kidia_var_1410 : kfdia_var_1411))
  i_nlayers_var_1441 = klev_var_1412
  DO iplon_var_1442 = kidia_var_1410, kfdia_var_1411
    i_laysolfr_var_1440(iplon_var_1442) = k_laytrop_var_1423(iplon_var_1442)
  END DO
  DO i_lay_var_1439 = 1, laytrop_min_var_1443
    DO iplon_var_1442 = kidia_var_1410, kfdia_var_1411
      IF (k_jp_var_1417(iplon_var_1442, i_lay_var_1439) < layreffr_var_250 .AND. k_jp_var_1417(iplon_var_1442, i_lay_var_1439 + 1) >= layreffr_var_250) i_laysolfr_var_1440(iplon_var_1442) = MIN(i_lay_var_1439 + 1, k_laytrop_var_1423(iplon_var_1442))
      ind0_var_1435 = ((k_jp_var_1417(iplon_var_1442, i_lay_var_1439) - 1) * 5 + (k_jt_var_1418(iplon_var_1442, i_lay_var_1439) - 1)) * nspa_var_313(20) + 1
      ind1_var_1436 = (k_jp_var_1417(iplon_var_1442, i_lay_var_1439) * 5 + (k_jt1_var_1419(iplon_var_1442, i_lay_var_1439) - 1)) * nspa_var_313(20) + 1
      inds_var_1437 = k_indself_var_1426(iplon_var_1442, i_lay_var_1439)
      indf_var_1438 = k_indfor_var_1429(iplon_var_1442, i_lay_var_1439)
      z_tauray_var_1445 = p_colmol_var_1422(iplon_var_1442, i_lay_var_1439) * rayl_var_249
      DO ig_var_1434 = 1, 10
        p_taug_var_1431(iplon_var_1442, i_lay_var_1439, ig_var_1434) = p_colh2o_var_1420(iplon_var_1442, i_lay_var_1439) * ((p_fac00_var_1413(iplon_var_1442, i_lay_var_1439) * absa_var_251(ind0_var_1435, ig_var_1434) + p_fac10_var_1415(iplon_var_1442, i_lay_var_1439) * absa_var_251(ind0_var_1435 + 1, ig_var_1434) + p_fac01_var_1414(iplon_var_1442, i_lay_var_1439) * absa_var_251(ind1_var_1436, ig_var_1434) + p_fac11_var_1416(iplon_var_1442, i_lay_var_1439) * absa_var_251(ind1_var_1436 + 1, ig_var_1434)) + p_selffac_var_1424(iplon_var_1442, i_lay_var_1439) * (selfrefc_var_253(inds_var_1437, ig_var_1434) + p_selffrac_var_1425(iplon_var_1442, i_lay_var_1439) * (selfrefc_var_253(inds_var_1437 + 1, ig_var_1434) - selfrefc_var_253(inds_var_1437, ig_var_1434))) + p_forfac_var_1427(iplon_var_1442, i_lay_var_1439) * (forrefc_var_254(indf_var_1438, ig_var_1434) + p_forfrac_var_1428(iplon_var_1442, i_lay_var_1439) * (forrefc_var_254(indf_var_1438 + 1, ig_var_1434) - forrefc_var_254(indf_var_1438, ig_var_1434)))) + p_colch4_var_1421(iplon_var_1442, i_lay_var_1439) * absch4c(ig_var_1434)
        p_taur_var_1432(iplon_var_1442, i_lay_var_1439, ig_var_1434) = z_tauray_var_1445
        IF (i_lay_var_1439 == i_laysolfr_var_1440(iplon_var_1442)) p_sfluxzen_var_1430(iplon_var_1442, ig_var_1434) = sfluxrefc_var_255(ig_var_1434)
      END DO
    END DO
  END DO
  DO i_lay_var_1439 = laytrop_min_var_1443 + 1, laytrop_max_var_1444
    DO iplon_var_1442 = kidia_var_1410, kfdia_var_1411
      IF (i_lay_var_1439 <= k_laytrop_var_1423(iplon_var_1442)) THEN
        IF (k_jp_var_1417(iplon_var_1442, i_lay_var_1439) < layreffr_var_250 .AND. k_jp_var_1417(iplon_var_1442, i_lay_var_1439 + 1) >= layreffr_var_250) i_laysolfr_var_1440(iplon_var_1442) = MIN(i_lay_var_1439 + 1, k_laytrop_var_1423(iplon_var_1442))
        ind0_var_1435 = ((k_jp_var_1417(iplon_var_1442, i_lay_var_1439) - 1) * 5 + (k_jt_var_1418(iplon_var_1442, i_lay_var_1439) - 1)) * nspa_var_313(20) + 1
        ind1_var_1436 = (k_jp_var_1417(iplon_var_1442, i_lay_var_1439) * 5 + (k_jt1_var_1419(iplon_var_1442, i_lay_var_1439) - 1)) * nspa_var_313(20) + 1
        inds_var_1437 = k_indself_var_1426(iplon_var_1442, i_lay_var_1439)
        indf_var_1438 = k_indfor_var_1429(iplon_var_1442, i_lay_var_1439)
        z_tauray_var_1445 = p_colmol_var_1422(iplon_var_1442, i_lay_var_1439) * rayl_var_249
        DO ig_var_1434 = 1, 10
          p_taug_var_1431(iplon_var_1442, i_lay_var_1439, ig_var_1434) = p_colh2o_var_1420(iplon_var_1442, i_lay_var_1439) * ((p_fac00_var_1413(iplon_var_1442, i_lay_var_1439) * absa_var_251(ind0_var_1435, ig_var_1434) + p_fac10_var_1415(iplon_var_1442, i_lay_var_1439) * absa_var_251(ind0_var_1435 + 1, ig_var_1434) + p_fac01_var_1414(iplon_var_1442, i_lay_var_1439) * absa_var_251(ind1_var_1436, ig_var_1434) + p_fac11_var_1416(iplon_var_1442, i_lay_var_1439) * absa_var_251(ind1_var_1436 + 1, ig_var_1434)) + p_selffac_var_1424(iplon_var_1442, i_lay_var_1439) * (selfrefc_var_253(inds_var_1437, ig_var_1434) + p_selffrac_var_1425(iplon_var_1442, i_lay_var_1439) * (selfrefc_var_253(inds_var_1437 + 1, ig_var_1434) - selfrefc_var_253(inds_var_1437, ig_var_1434))) + p_forfac_var_1427(iplon_var_1442, i_lay_var_1439) * (forrefc_var_254(indf_var_1438, ig_var_1434) + p_forfrac_var_1428(iplon_var_1442, i_lay_var_1439) * (forrefc_var_254(indf_var_1438 + 1, ig_var_1434) - forrefc_var_254(indf_var_1438, ig_var_1434)))) + p_colch4_var_1421(iplon_var_1442, i_lay_var_1439) * absch4c(ig_var_1434)
          p_taur_var_1432(iplon_var_1442, i_lay_var_1439, ig_var_1434) = z_tauray_var_1445
          IF (i_lay_var_1439 == i_laysolfr_var_1440(iplon_var_1442)) p_sfluxzen_var_1430(iplon_var_1442, ig_var_1434) = sfluxrefc_var_255(ig_var_1434)
        END DO
      ELSE
        ind0_var_1435 = ((k_jp_var_1417(iplon_var_1442, i_lay_var_1439) - 13) * 5 + (k_jt_var_1418(iplon_var_1442, i_lay_var_1439) - 1)) * nspb_var_314(20) + 1
        ind1_var_1436 = ((k_jp_var_1417(iplon_var_1442, i_lay_var_1439) - 12) * 5 + (k_jt1_var_1419(iplon_var_1442, i_lay_var_1439) - 1)) * nspb_var_314(20) + 1
        indf_var_1438 = k_indfor_var_1429(iplon_var_1442, i_lay_var_1439)
        z_tauray_var_1445 = p_colmol_var_1422(iplon_var_1442, i_lay_var_1439) * rayl_var_249
        DO ig_var_1434 = 1, 10
          p_taug_var_1431(iplon_var_1442, i_lay_var_1439, ig_var_1434) = p_colh2o_var_1420(iplon_var_1442, i_lay_var_1439) * (p_fac00_var_1413(iplon_var_1442, i_lay_var_1439) * absb_var_252(ind0_var_1435, ig_var_1434) + p_fac10_var_1415(iplon_var_1442, i_lay_var_1439) * absb_var_252(ind0_var_1435 + 1, ig_var_1434) + p_fac01_var_1414(iplon_var_1442, i_lay_var_1439) * absb_var_252(ind1_var_1436, ig_var_1434) + p_fac11_var_1416(iplon_var_1442, i_lay_var_1439) * absb_var_252(ind1_var_1436 + 1, ig_var_1434) + p_forfac_var_1427(iplon_var_1442, i_lay_var_1439) * (forrefc_var_254(indf_var_1438, ig_var_1434) + p_forfrac_var_1428(iplon_var_1442, i_lay_var_1439) * (forrefc_var_254(indf_var_1438 + 1, ig_var_1434) - forrefc_var_254(indf_var_1438, ig_var_1434)))) + p_colch4_var_1421(iplon_var_1442, i_lay_var_1439) * absch4c(ig_var_1434)
          p_taur_var_1432(iplon_var_1442, i_lay_var_1439, ig_var_1434) = z_tauray_var_1445
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1439 = laytrop_max_var_1444 + 1, i_nlayers_var_1441
    DO iplon_var_1442 = kidia_var_1410, kfdia_var_1411
      ind0_var_1435 = ((k_jp_var_1417(iplon_var_1442, i_lay_var_1439) - 13) * 5 + (k_jt_var_1418(iplon_var_1442, i_lay_var_1439) - 1)) * nspb_var_314(20) + 1
      ind1_var_1436 = ((k_jp_var_1417(iplon_var_1442, i_lay_var_1439) - 12) * 5 + (k_jt1_var_1419(iplon_var_1442, i_lay_var_1439) - 1)) * nspb_var_314(20) + 1
      indf_var_1438 = k_indfor_var_1429(iplon_var_1442, i_lay_var_1439)
      z_tauray_var_1445 = p_colmol_var_1422(iplon_var_1442, i_lay_var_1439) * rayl_var_249
      DO ig_var_1434 = 1, 10
        p_taug_var_1431(iplon_var_1442, i_lay_var_1439, ig_var_1434) = p_colh2o_var_1420(iplon_var_1442, i_lay_var_1439) * (p_fac00_var_1413(iplon_var_1442, i_lay_var_1439) * absb_var_252(ind0_var_1435, ig_var_1434) + p_fac10_var_1415(iplon_var_1442, i_lay_var_1439) * absb_var_252(ind0_var_1435 + 1, ig_var_1434) + p_fac01_var_1414(iplon_var_1442, i_lay_var_1439) * absb_var_252(ind1_var_1436, ig_var_1434) + p_fac11_var_1416(iplon_var_1442, i_lay_var_1439) * absb_var_252(ind1_var_1436 + 1, ig_var_1434) + p_forfac_var_1427(iplon_var_1442, i_lay_var_1439) * (forrefc_var_254(indf_var_1438, ig_var_1434) + p_forfrac_var_1428(iplon_var_1442, i_lay_var_1439) * (forrefc_var_254(indf_var_1438 + 1, ig_var_1434) - forrefc_var_254(indf_var_1438, ig_var_1434)))) + p_colch4_var_1421(iplon_var_1442, i_lay_var_1439) * absch4c(ig_var_1434)
        p_taur_var_1432(iplon_var_1442, i_lay_var_1439, ig_var_1434) = z_tauray_var_1445
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol20
SUBROUTINE srtm_taumol22(kidia_var_1446, kfdia_var_1447, klev_var_1448, p_fac00_var_1449, p_fac01_var_1450, p_fac10_var_1451, p_fac11_var_1452, k_jp_var_1453, k_jt_var_1454, k_jt1_var_1455, p_oneminus_var_1456, p_colh2o_var_1457, p_colmol_var_1458, p_colo2_var_1459, k_laytrop_var_1460, p_selffac_var_1461, p_selffrac_var_1462, k_indself_var_1463, p_forfac_var_1464, p_forfrac_var_1465, k_indfor_var_1466, p_sfluxzen_var_1467, p_taug_var_1468, p_taur_var_1469, prmu0_var_1470)
  USE yoesrta22, ONLY: absa_var_267, absb_var_268, forrefc_var_270, layreffr_var_266, rayl_var_264, selfrefc_var_269, sfluxrefc_var_271, strrat_var_265
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1446, kfdia_var_1447
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1448
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1449(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1450(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1451(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1452(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1453(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1454(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1455(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1456(kidia_var_1446 : kfdia_var_1447)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1457(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1458(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_1459(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1460(kidia_var_1446 : kfdia_var_1447)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1461(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1462(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1463(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1464(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1465(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1466(kidia_var_1446 : kfdia_var_1447, klev_var_1448)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1467(kidia_var_1446 : kfdia_var_1447, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1468(kidia_var_1446 : kfdia_var_1447, klev_var_1448, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1469(kidia_var_1446 : kfdia_var_1447, klev_var_1448, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1470(kidia_var_1446 : kfdia_var_1447)
  INTEGER(KIND = 4) :: ig_var_1471, ind0_var_1472, ind1_var_1473, inds_var_1474, indf_var_1475, js_var_1476, i_lay_var_1477, i_laysolfr_var_1478(kidia_var_1446 : kfdia_var_1447), i_nlayers_var_1479, iplon_var_1480
  INTEGER(KIND = 4) :: laytrop_min_var_1481, laytrop_max_var_1482
  REAL(KIND = 8) :: z_fs_var_1483, z_speccomb_var_1484, z_specmult_var_1485, z_specparm_var_1486, z_tauray_var_1487, z_o2adj, z_o2cont
  laytrop_min_var_1481 = MINVAL(k_laytrop_var_1460(kidia_var_1446 : kfdia_var_1447))
  laytrop_max_var_1482 = MAXVAL(k_laytrop_var_1460(kidia_var_1446 : kfdia_var_1447))
  i_nlayers_var_1479 = klev_var_1448
  z_o2adj = 1.6D0
  DO iplon_var_1480 = kidia_var_1446, kfdia_var_1447
    i_laysolfr_var_1478(iplon_var_1480) = k_laytrop_var_1460(iplon_var_1480)
  END DO
  DO i_lay_var_1477 = 1, laytrop_min_var_1481
    DO iplon_var_1480 = kidia_var_1446, kfdia_var_1447
      IF (k_jp_var_1453(iplon_var_1480, i_lay_var_1477) < layreffr_var_266 .AND. k_jp_var_1453(iplon_var_1480, i_lay_var_1477 + 1) >= layreffr_var_266) i_laysolfr_var_1478(iplon_var_1480) = MIN(i_lay_var_1477 + 1, k_laytrop_var_1460(iplon_var_1480))
      z_o2cont = 0.000435D0 * p_colo2_var_1459(iplon_var_1480, i_lay_var_1477) / (700.0D0)
      z_speccomb_var_1484 = p_colh2o_var_1457(iplon_var_1480, i_lay_var_1477) + z_o2adj * strrat_var_265 * p_colo2_var_1459(iplon_var_1480, i_lay_var_1477)
      z_specparm_var_1486 = p_colh2o_var_1457(iplon_var_1480, i_lay_var_1477) / z_speccomb_var_1484
      z_specparm_var_1486 = MIN(p_oneminus_var_1456(iplon_var_1480), z_specparm_var_1486)
      z_specmult_var_1485 = 8.0D0 * (z_specparm_var_1486)
      js_var_1476 = 1 + INT(z_specmult_var_1485)
      z_fs_var_1483 = z_specmult_var_1485 - AINT(z_specmult_var_1485)
      ind0_var_1472 = ((k_jp_var_1453(iplon_var_1480, i_lay_var_1477) - 1) * 5 + (k_jt_var_1454(iplon_var_1480, i_lay_var_1477) - 1)) * nspa_var_313(22) + js_var_1476
      ind1_var_1473 = (k_jp_var_1453(iplon_var_1480, i_lay_var_1477) * 5 + (k_jt1_var_1455(iplon_var_1480, i_lay_var_1477) - 1)) * nspa_var_313(22) + js_var_1476
      inds_var_1474 = k_indself_var_1463(iplon_var_1480, i_lay_var_1477)
      indf_var_1475 = k_indfor_var_1466(iplon_var_1480, i_lay_var_1477)
      z_tauray_var_1487 = p_colmol_var_1458(iplon_var_1480, i_lay_var_1477) * rayl_var_264
      DO ig_var_1471 = 1, 2
        p_taug_var_1468(iplon_var_1480, i_lay_var_1477, ig_var_1471) = z_speccomb_var_1484 * ((1.0D0 - z_fs_var_1483) * (absa_var_267(ind0_var_1472, ig_var_1471) * p_fac00_var_1449(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind0_var_1472 + 9, ig_var_1471) * p_fac10_var_1451(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind1_var_1473, ig_var_1471) * p_fac01_var_1450(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind1_var_1473 + 9, ig_var_1471) * p_fac11_var_1452(iplon_var_1480, i_lay_var_1477)) + z_fs_var_1483 * (absa_var_267(ind0_var_1472 + 1, ig_var_1471) * p_fac00_var_1449(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind0_var_1472 + 10, ig_var_1471) * p_fac10_var_1451(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind1_var_1473 + 1, ig_var_1471) * p_fac01_var_1450(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind1_var_1473 + 10, ig_var_1471) * p_fac11_var_1452(iplon_var_1480, i_lay_var_1477))) + p_colh2o_var_1457(iplon_var_1480, i_lay_var_1477) * (p_selffac_var_1461(iplon_var_1480, i_lay_var_1477) * (selfrefc_var_269(inds_var_1474, ig_var_1471) + p_selffrac_var_1462(iplon_var_1480, i_lay_var_1477) * (selfrefc_var_269(inds_var_1474 + 1, ig_var_1471) - selfrefc_var_269(inds_var_1474, ig_var_1471))) + p_forfac_var_1464(iplon_var_1480, i_lay_var_1477) * (forrefc_var_270(indf_var_1475, ig_var_1471) + p_forfrac_var_1465(iplon_var_1480, i_lay_var_1477) * (forrefc_var_270(indf_var_1475 + 1, ig_var_1471) - forrefc_var_270(indf_var_1475, ig_var_1471)))) + z_o2cont
        IF (i_lay_var_1477 == i_laysolfr_var_1478(iplon_var_1480)) p_sfluxzen_var_1467(iplon_var_1480, ig_var_1471) = sfluxrefc_var_271(ig_var_1471, js_var_1476) + z_fs_var_1483 * (sfluxrefc_var_271(ig_var_1471, js_var_1476 + 1) - sfluxrefc_var_271(ig_var_1471, js_var_1476))
        p_taur_var_1469(iplon_var_1480, i_lay_var_1477, ig_var_1471) = z_tauray_var_1487
      END DO
    END DO
  END DO
  DO i_lay_var_1477 = laytrop_min_var_1481 + 1, laytrop_max_var_1482
    DO iplon_var_1480 = kidia_var_1446, kfdia_var_1447
      IF (i_lay_var_1477 <= k_laytrop_var_1460(iplon_var_1480)) THEN
        IF (k_jp_var_1453(iplon_var_1480, i_lay_var_1477) < layreffr_var_266 .AND. k_jp_var_1453(iplon_var_1480, i_lay_var_1477 + 1) >= layreffr_var_266) i_laysolfr_var_1478(iplon_var_1480) = MIN(i_lay_var_1477 + 1, k_laytrop_var_1460(iplon_var_1480))
        z_o2cont = 0.000435D0 * p_colo2_var_1459(iplon_var_1480, i_lay_var_1477) / (700.0D0)
        z_speccomb_var_1484 = p_colh2o_var_1457(iplon_var_1480, i_lay_var_1477) + z_o2adj * strrat_var_265 * p_colo2_var_1459(iplon_var_1480, i_lay_var_1477)
        z_specparm_var_1486 = p_colh2o_var_1457(iplon_var_1480, i_lay_var_1477) / z_speccomb_var_1484
        z_specparm_var_1486 = MIN(p_oneminus_var_1456(iplon_var_1480), z_specparm_var_1486)
        z_specmult_var_1485 = 8.0D0 * (z_specparm_var_1486)
        js_var_1476 = 1 + INT(z_specmult_var_1485)
        z_fs_var_1483 = z_specmult_var_1485 - AINT(z_specmult_var_1485)
        ind0_var_1472 = ((k_jp_var_1453(iplon_var_1480, i_lay_var_1477) - 1) * 5 + (k_jt_var_1454(iplon_var_1480, i_lay_var_1477) - 1)) * nspa_var_313(22) + js_var_1476
        ind1_var_1473 = (k_jp_var_1453(iplon_var_1480, i_lay_var_1477) * 5 + (k_jt1_var_1455(iplon_var_1480, i_lay_var_1477) - 1)) * nspa_var_313(22) + js_var_1476
        inds_var_1474 = k_indself_var_1463(iplon_var_1480, i_lay_var_1477)
        indf_var_1475 = k_indfor_var_1466(iplon_var_1480, i_lay_var_1477)
        z_tauray_var_1487 = p_colmol_var_1458(iplon_var_1480, i_lay_var_1477) * rayl_var_264
        DO ig_var_1471 = 1, 2
          p_taug_var_1468(iplon_var_1480, i_lay_var_1477, ig_var_1471) = z_speccomb_var_1484 * ((1.0D0 - z_fs_var_1483) * (absa_var_267(ind0_var_1472, ig_var_1471) * p_fac00_var_1449(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind0_var_1472 + 9, ig_var_1471) * p_fac10_var_1451(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind1_var_1473, ig_var_1471) * p_fac01_var_1450(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind1_var_1473 + 9, ig_var_1471) * p_fac11_var_1452(iplon_var_1480, i_lay_var_1477)) + z_fs_var_1483 * (absa_var_267(ind0_var_1472 + 1, ig_var_1471) * p_fac00_var_1449(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind0_var_1472 + 10, ig_var_1471) * p_fac10_var_1451(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind1_var_1473 + 1, ig_var_1471) * p_fac01_var_1450(iplon_var_1480, i_lay_var_1477) + absa_var_267(ind1_var_1473 + 10, ig_var_1471) * p_fac11_var_1452(iplon_var_1480, i_lay_var_1477))) + p_colh2o_var_1457(iplon_var_1480, i_lay_var_1477) * (p_selffac_var_1461(iplon_var_1480, i_lay_var_1477) * (selfrefc_var_269(inds_var_1474, ig_var_1471) + p_selffrac_var_1462(iplon_var_1480, i_lay_var_1477) * (selfrefc_var_269(inds_var_1474 + 1, ig_var_1471) - selfrefc_var_269(inds_var_1474, ig_var_1471))) + p_forfac_var_1464(iplon_var_1480, i_lay_var_1477) * (forrefc_var_270(indf_var_1475, ig_var_1471) + p_forfrac_var_1465(iplon_var_1480, i_lay_var_1477) * (forrefc_var_270(indf_var_1475 + 1, ig_var_1471) - forrefc_var_270(indf_var_1475, ig_var_1471)))) + z_o2cont
          IF (i_lay_var_1477 == i_laysolfr_var_1478(iplon_var_1480)) p_sfluxzen_var_1467(iplon_var_1480, ig_var_1471) = sfluxrefc_var_271(ig_var_1471, js_var_1476) + z_fs_var_1483 * (sfluxrefc_var_271(ig_var_1471, js_var_1476 + 1) - sfluxrefc_var_271(ig_var_1471, js_var_1476))
          p_taur_var_1469(iplon_var_1480, i_lay_var_1477, ig_var_1471) = z_tauray_var_1487
        END DO
      ELSE
        z_o2cont = 0.000435D0 * p_colo2_var_1459(iplon_var_1480, i_lay_var_1477) / (700.0D0)
        ind0_var_1472 = ((k_jp_var_1453(iplon_var_1480, i_lay_var_1477) - 13) * 5 + (k_jt_var_1454(iplon_var_1480, i_lay_var_1477) - 1)) * nspb_var_314(22) + 1
        ind1_var_1473 = ((k_jp_var_1453(iplon_var_1480, i_lay_var_1477) - 12) * 5 + (k_jt1_var_1455(iplon_var_1480, i_lay_var_1477) - 1)) * nspb_var_314(22) + 1
        z_tauray_var_1487 = p_colmol_var_1458(iplon_var_1480, i_lay_var_1477) * rayl_var_264
        DO ig_var_1471 = 1, 2
          p_taug_var_1468(iplon_var_1480, i_lay_var_1477, ig_var_1471) = p_colo2_var_1459(iplon_var_1480, i_lay_var_1477) * z_o2adj * (p_fac00_var_1449(iplon_var_1480, i_lay_var_1477) * absb_var_268(ind0_var_1472, ig_var_1471) + p_fac10_var_1451(iplon_var_1480, i_lay_var_1477) * absb_var_268(ind0_var_1472 + 1, ig_var_1471) + p_fac01_var_1450(iplon_var_1480, i_lay_var_1477) * absb_var_268(ind1_var_1473, ig_var_1471) + p_fac11_var_1452(iplon_var_1480, i_lay_var_1477) * absb_var_268(ind1_var_1473 + 1, ig_var_1471)) + z_o2cont
          p_taur_var_1469(iplon_var_1480, i_lay_var_1477, ig_var_1471) = z_tauray_var_1487
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1477 = laytrop_max_var_1482 + 1, i_nlayers_var_1479
    DO iplon_var_1480 = kidia_var_1446, kfdia_var_1447
      z_o2cont = 0.000435D0 * p_colo2_var_1459(iplon_var_1480, i_lay_var_1477) / (700.0D0)
      ind0_var_1472 = ((k_jp_var_1453(iplon_var_1480, i_lay_var_1477) - 13) * 5 + (k_jt_var_1454(iplon_var_1480, i_lay_var_1477) - 1)) * nspb_var_314(22) + 1
      ind1_var_1473 = ((k_jp_var_1453(iplon_var_1480, i_lay_var_1477) - 12) * 5 + (k_jt1_var_1455(iplon_var_1480, i_lay_var_1477) - 1)) * nspb_var_314(22) + 1
      z_tauray_var_1487 = p_colmol_var_1458(iplon_var_1480, i_lay_var_1477) * rayl_var_264
      DO ig_var_1471 = 1, 2
        p_taug_var_1468(iplon_var_1480, i_lay_var_1477, ig_var_1471) = p_colo2_var_1459(iplon_var_1480, i_lay_var_1477) * z_o2adj * (p_fac00_var_1449(iplon_var_1480, i_lay_var_1477) * absb_var_268(ind0_var_1472, ig_var_1471) + p_fac10_var_1451(iplon_var_1480, i_lay_var_1477) * absb_var_268(ind0_var_1472 + 1, ig_var_1471) + p_fac01_var_1450(iplon_var_1480, i_lay_var_1477) * absb_var_268(ind1_var_1473, ig_var_1471) + p_fac11_var_1452(iplon_var_1480, i_lay_var_1477) * absb_var_268(ind1_var_1473 + 1, ig_var_1471)) + z_o2cont
        p_taur_var_1469(iplon_var_1480, i_lay_var_1477, ig_var_1471) = z_tauray_var_1487
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol22
SUBROUTINE srtm_taumol23(kidia_var_1488, kfdia_var_1489, klev_var_1490, p_fac00_var_1491, p_fac01_var_1492, p_fac10_var_1493, p_fac11_var_1494, k_jp_var_1495, k_jt_var_1496, k_jt1_var_1497, p_colh2o_var_1498, p_colmol_var_1499, k_laytrop_var_1500, p_selffac_var_1501, p_selffrac_var_1502, k_indself_var_1503, p_forfac_var_1504, p_forfrac_var_1505, k_indfor_var_1506, p_sfluxzen_var_1507, p_taug_var_1508, p_taur_var_1509, prmu0_var_1510)
  USE yoesrta23, ONLY: absa_var_273, forrefc_var_275, givfac, layreffr_var_272, raylc_var_277, selfrefc_var_274, sfluxrefc_var_276
  USE yoesrtwn, ONLY: nspa_var_313
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1488, kfdia_var_1489
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1490
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1491(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1492(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1493(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1494(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1495(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1496(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1497(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1498(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1499(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1500(kidia_var_1488 : kfdia_var_1489)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1501(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1502(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1503(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1504(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1505(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1506(kidia_var_1488 : kfdia_var_1489, klev_var_1490)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1507(kidia_var_1488 : kfdia_var_1489, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1508(kidia_var_1488 : kfdia_var_1489, klev_var_1490, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1509(kidia_var_1488 : kfdia_var_1489, klev_var_1490, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1510(kidia_var_1488 : kfdia_var_1489)
  INTEGER(KIND = 4) :: ig_var_1511, ind0_var_1512, ind1_var_1513, inds_var_1514, indf_var_1515, i_lay_var_1516, i_laysolfr_var_1517(kidia_var_1488 : kfdia_var_1489), i_nlayers_var_1518, iplon_var_1519
  INTEGER(KIND = 4) :: laytrop_min_var_1520, laytrop_max_var_1521
  REAL(KIND = 8) :: z_tauray_var_1522
  laytrop_min_var_1520 = MINVAL(k_laytrop_var_1500(kidia_var_1488 : kfdia_var_1489))
  laytrop_max_var_1521 = MAXVAL(k_laytrop_var_1500(kidia_var_1488 : kfdia_var_1489))
  i_nlayers_var_1518 = klev_var_1490
  DO iplon_var_1519 = kidia_var_1488, kfdia_var_1489
    i_laysolfr_var_1517(iplon_var_1519) = k_laytrop_var_1500(iplon_var_1519)
  END DO
  DO i_lay_var_1516 = 1, laytrop_min_var_1520
    DO iplon_var_1519 = kidia_var_1488, kfdia_var_1489
      IF (k_jp_var_1495(iplon_var_1519, i_lay_var_1516) < layreffr_var_272 .AND. k_jp_var_1495(iplon_var_1519, i_lay_var_1516 + 1) >= layreffr_var_272) i_laysolfr_var_1517(iplon_var_1519) = MIN(i_lay_var_1516 + 1, k_laytrop_var_1500(iplon_var_1519))
      ind0_var_1512 = ((k_jp_var_1495(iplon_var_1519, i_lay_var_1516) - 1) * 5 + (k_jt_var_1496(iplon_var_1519, i_lay_var_1516) - 1)) * nspa_var_313(23) + 1
      ind1_var_1513 = (k_jp_var_1495(iplon_var_1519, i_lay_var_1516) * 5 + (k_jt1_var_1497(iplon_var_1519, i_lay_var_1516) - 1)) * nspa_var_313(23) + 1
      inds_var_1514 = k_indself_var_1503(iplon_var_1519, i_lay_var_1516)
      indf_var_1515 = k_indfor_var_1506(iplon_var_1519, i_lay_var_1516)
      DO ig_var_1511 = 1, 10
        z_tauray_var_1522 = p_colmol_var_1499(iplon_var_1519, i_lay_var_1516) * raylc_var_277(ig_var_1511)
        p_taug_var_1508(iplon_var_1519, i_lay_var_1516, ig_var_1511) = p_colh2o_var_1498(iplon_var_1519, i_lay_var_1516) * (givfac * (p_fac00_var_1491(iplon_var_1519, i_lay_var_1516) * absa_var_273(ind0_var_1512, ig_var_1511) + p_fac10_var_1493(iplon_var_1519, i_lay_var_1516) * absa_var_273(ind0_var_1512 + 1, ig_var_1511) + p_fac01_var_1492(iplon_var_1519, i_lay_var_1516) * absa_var_273(ind1_var_1513, ig_var_1511) + p_fac11_var_1494(iplon_var_1519, i_lay_var_1516) * absa_var_273(ind1_var_1513 + 1, ig_var_1511)) + p_selffac_var_1501(iplon_var_1519, i_lay_var_1516) * (selfrefc_var_274(inds_var_1514, ig_var_1511) + p_selffrac_var_1502(iplon_var_1519, i_lay_var_1516) * (selfrefc_var_274(inds_var_1514 + 1, ig_var_1511) - selfrefc_var_274(inds_var_1514, ig_var_1511))) + p_forfac_var_1504(iplon_var_1519, i_lay_var_1516) * (forrefc_var_275(indf_var_1515, ig_var_1511) + p_forfrac_var_1505(iplon_var_1519, i_lay_var_1516) * (forrefc_var_275(indf_var_1515 + 1, ig_var_1511) - forrefc_var_275(indf_var_1515, ig_var_1511))))
        IF (i_lay_var_1516 == i_laysolfr_var_1517(iplon_var_1519)) p_sfluxzen_var_1507(iplon_var_1519, ig_var_1511) = sfluxrefc_var_276(ig_var_1511)
        p_taur_var_1509(iplon_var_1519, i_lay_var_1516, ig_var_1511) = z_tauray_var_1522
      END DO
    END DO
  END DO
  DO i_lay_var_1516 = laytrop_min_var_1520 + 1, laytrop_max_var_1521
    DO iplon_var_1519 = kidia_var_1488, kfdia_var_1489
      IF (i_lay_var_1516 <= k_laytrop_var_1500(iplon_var_1519)) THEN
        IF (k_jp_var_1495(iplon_var_1519, i_lay_var_1516) < layreffr_var_272 .AND. k_jp_var_1495(iplon_var_1519, i_lay_var_1516 + 1) >= layreffr_var_272) i_laysolfr_var_1517(iplon_var_1519) = MIN(i_lay_var_1516 + 1, k_laytrop_var_1500(iplon_var_1519))
        ind0_var_1512 = ((k_jp_var_1495(iplon_var_1519, i_lay_var_1516) - 1) * 5 + (k_jt_var_1496(iplon_var_1519, i_lay_var_1516) - 1)) * nspa_var_313(23) + 1
        ind1_var_1513 = (k_jp_var_1495(iplon_var_1519, i_lay_var_1516) * 5 + (k_jt1_var_1497(iplon_var_1519, i_lay_var_1516) - 1)) * nspa_var_313(23) + 1
        inds_var_1514 = k_indself_var_1503(iplon_var_1519, i_lay_var_1516)
        indf_var_1515 = k_indfor_var_1506(iplon_var_1519, i_lay_var_1516)
        DO ig_var_1511 = 1, 10
          z_tauray_var_1522 = p_colmol_var_1499(iplon_var_1519, i_lay_var_1516) * raylc_var_277(ig_var_1511)
          p_taug_var_1508(iplon_var_1519, i_lay_var_1516, ig_var_1511) = p_colh2o_var_1498(iplon_var_1519, i_lay_var_1516) * (givfac * (p_fac00_var_1491(iplon_var_1519, i_lay_var_1516) * absa_var_273(ind0_var_1512, ig_var_1511) + p_fac10_var_1493(iplon_var_1519, i_lay_var_1516) * absa_var_273(ind0_var_1512 + 1, ig_var_1511) + p_fac01_var_1492(iplon_var_1519, i_lay_var_1516) * absa_var_273(ind1_var_1513, ig_var_1511) + p_fac11_var_1494(iplon_var_1519, i_lay_var_1516) * absa_var_273(ind1_var_1513 + 1, ig_var_1511)) + p_selffac_var_1501(iplon_var_1519, i_lay_var_1516) * (selfrefc_var_274(inds_var_1514, ig_var_1511) + p_selffrac_var_1502(iplon_var_1519, i_lay_var_1516) * (selfrefc_var_274(inds_var_1514 + 1, ig_var_1511) - selfrefc_var_274(inds_var_1514, ig_var_1511))) + p_forfac_var_1504(iplon_var_1519, i_lay_var_1516) * (forrefc_var_275(indf_var_1515, ig_var_1511) + p_forfrac_var_1505(iplon_var_1519, i_lay_var_1516) * (forrefc_var_275(indf_var_1515 + 1, ig_var_1511) - forrefc_var_275(indf_var_1515, ig_var_1511))))
          IF (i_lay_var_1516 == i_laysolfr_var_1517(iplon_var_1519)) p_sfluxzen_var_1507(iplon_var_1519, ig_var_1511) = sfluxrefc_var_276(ig_var_1511)
          p_taur_var_1509(iplon_var_1519, i_lay_var_1516, ig_var_1511) = z_tauray_var_1522
        END DO
      ELSE
        DO ig_var_1511 = 1, 10
          p_taug_var_1508(iplon_var_1519, i_lay_var_1516, ig_var_1511) = 0.0D0
          p_taur_var_1509(iplon_var_1519, i_lay_var_1516, ig_var_1511) = p_colmol_var_1499(iplon_var_1519, i_lay_var_1516) * raylc_var_277(ig_var_1511)
        END DO
      END IF
    END DO
  END DO
  DO ig_var_1511 = 1, 10
    DO i_lay_var_1516 = laytrop_max_var_1521 + 1, i_nlayers_var_1518
      DO iplon_var_1519 = kidia_var_1488, kfdia_var_1489
        p_taug_var_1508(iplon_var_1519, i_lay_var_1516, ig_var_1511) = 0.0D0
        p_taur_var_1509(iplon_var_1519, i_lay_var_1516, ig_var_1511) = p_colmol_var_1499(iplon_var_1519, i_lay_var_1516) * raylc_var_277(ig_var_1511)
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol23
SUBROUTINE rrtm_taumol5(kidia_var_1523, kfdia_var_1524, klev_var_1525, taug_var_1526, wx_var_1527, p_tauaerl_var_1528, fac00_var_1529, fac01_var_1530, fac10_var_1531, fac11_var_1532, forfac_var_1551, forfrac_var_1550, indfor_var_1549, jp_var_1533, jt_var_1534, jt1_var_1535, oneminus_var_1536, colh2o_var_1537, colco2_var_1538, colo3_var_1539, laytrop_var_1540, selffac_var_1541, selffrac_var_1542, indself_var_1543, fracs_var_1544, rat_h2oco2_var_1545, rat_h2oco2_1_var_1546, rat_o3co2_var_1547, rat_o3co2_1_var_1548, minorfrac_var_1552, indminor_var_1553)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng5
  USE yoerrta5, ONLY: absa_var_175, absb_var_176, ccl4, forref_var_179, fracrefa_var_173, fracrefb_var_174, ka_mo3_var_177, selfref_var_178
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1523
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1524
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1525
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1526(kidia_var_1523 : kfdia_var_1524, 140, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: wx_var_1527(kidia_var_1523 : kfdia_var_1524, 4, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1528(kidia_var_1523 : kfdia_var_1524, klev_var_1525, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1529(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1530(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1531(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1532(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1533(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1534(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1535(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1536
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1537(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1538(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1539(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1540(kidia_var_1523 : kfdia_var_1524)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1541(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1542(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1543(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1544(kidia_var_1523 : kfdia_var_1524, 140, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1545(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1546(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_var_1547(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_1_var_1548(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1549(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1550(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1551(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1552(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1553(kidia_var_1523 : kfdia_var_1524, klev_var_1525)
  REAL(KIND = 8) :: speccomb_var_1554, speccomb1_var_1555, speccomb_mo3, speccomb_planck_var_1556
  INTEGER(KIND = 4) :: ind0_var_1557, ind1_var_1558, inds_var_1559, indf_var_1560, indm_var_1561
  INTEGER(KIND = 4) :: ig_var_1562, js_var_1563, lay_var_1564, js1_var_1565, jpl_var_1566, jmo3
  REAL(KIND = 8) :: refrat_planck_a_var_1567, refrat_planck_b_var_1568, refrat_m_a_var_1569
  REAL(KIND = 8) :: fac000_var_1570, fac100_var_1571, fac200_var_1572, fac010_var_1573, fac110_var_1574, fac210_var_1575, fac001_var_1576, fac101_var_1577, fac201_var_1578, fac011_var_1579, fac111_var_1580, fac211_var_1581
  REAL(KIND = 8) :: p_var_1582, p4_var_1583, fk0_var_1584, fk1_var_1585, fk2_var_1586
  REAL(KIND = 8) :: taufor_var_1587, tauself_var_1588, tau_major_var_1589(16), tau_major1_var_1590(16), o3m1, o3m2, abso3_var_1591
  REAL(KIND = 8) :: fs_var_1592, specmult_var_1593, specparm_var_1594, fs1_var_1595, specmult1_var_1596, specparm1_var_1597, fpl_var_1598, specmult_planck_var_1599, specparm_planck_var_1600, fmo3, specmult_mo3, specparm_mo3
  INTEGER(KIND = 4) :: laytrop_min_var_1601, laytrop_max_var_1602
  INTEGER(KIND = 4) :: ixc_var_1603(klev_var_1525), ixlow_var_1604(kfdia_var_1524, klev_var_1525), ixhigh_var_1605(kfdia_var_1524, klev_var_1525)
  INTEGER(KIND = 4) :: ich_var_1606, icl_var_1607, ixc0_var_1608, ixp_var_1609, jc_var_1610, jl_var_1611
  laytrop_min_var_1601 = MINVAL(laytrop_var_1540)
  laytrop_max_var_1602 = MAXVAL(laytrop_var_1540)
  ixlow_var_1604 = 0
  ixhigh_var_1605 = 0
  ixc_var_1603 = 0
  DO lay_var_1564 = laytrop_min_var_1601 + 1, laytrop_max_var_1602
    icl_var_1607 = 0
    ich_var_1606 = 0
    DO jc_var_1610 = kidia_var_1523, kfdia_var_1524
      IF (lay_var_1564 <= laytrop_var_1540(jc_var_1610)) THEN
        icl_var_1607 = icl_var_1607 + 1
        ixlow_var_1604(icl_var_1607, lay_var_1564) = jc_var_1610
      ELSE
        ich_var_1606 = ich_var_1606 + 1
        ixhigh_var_1605(ich_var_1606, lay_var_1564) = jc_var_1610
      END IF
    END DO
    ixc_var_1603(lay_var_1564) = icl_var_1607
  END DO
  refrat_planck_a_var_1567 = chi_mls(1, 5) / chi_mls(2, 5)
  refrat_planck_b_var_1568 = chi_mls(3, 43) / chi_mls(2, 43)
  refrat_m_a_var_1569 = chi_mls(1, 7) / chi_mls(2, 7)
  DO lay_var_1564 = 1, laytrop_min_var_1601
    DO jl_var_1611 = kidia_var_1523, kfdia_var_1524
      speccomb_var_1554 = colh2o_var_1537(jl_var_1611, lay_var_1564) + rat_h2oco2_var_1545(jl_var_1611, lay_var_1564) * colco2_var_1538(jl_var_1611, lay_var_1564)
      specparm_var_1594 = MIN(colh2o_var_1537(jl_var_1611, lay_var_1564) / speccomb_var_1554, oneminus_var_1536)
      specmult_var_1593 = 8.0D0 * (specparm_var_1594)
      js_var_1563 = 1 + INT(specmult_var_1593)
      fs_var_1592 = ((specmult_var_1593) - AINT((specmult_var_1593)))
      speccomb1_var_1555 = colh2o_var_1537(jl_var_1611, lay_var_1564) + rat_h2oco2_1_var_1546(jl_var_1611, lay_var_1564) * colco2_var_1538(jl_var_1611, lay_var_1564)
      specparm1_var_1597 = MIN(colh2o_var_1537(jl_var_1611, lay_var_1564) / speccomb1_var_1555, oneminus_var_1536)
      specmult1_var_1596 = 8.0D0 * (specparm1_var_1597)
      js1_var_1565 = 1 + INT(specmult1_var_1596)
      fs1_var_1595 = ((specmult1_var_1596) - AINT((specmult1_var_1596)))
      speccomb_mo3 = colh2o_var_1537(jl_var_1611, lay_var_1564) + refrat_m_a_var_1569 * colco2_var_1538(jl_var_1611, lay_var_1564)
      specparm_mo3 = MIN(colh2o_var_1537(jl_var_1611, lay_var_1564) / speccomb_mo3, oneminus_var_1536)
      specmult_mo3 = 8.0D0 * specparm_mo3
      jmo3 = 1 + INT(specmult_mo3)
      fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))
      speccomb_planck_var_1556 = colh2o_var_1537(jl_var_1611, lay_var_1564) + refrat_planck_a_var_1567 * colco2_var_1538(jl_var_1611, lay_var_1564)
      specparm_planck_var_1600 = MIN(colh2o_var_1537(jl_var_1611, lay_var_1564) / speccomb_planck_var_1556, oneminus_var_1536)
      specmult_planck_var_1599 = 8.0D0 * specparm_planck_var_1600
      jpl_var_1566 = 1 + INT(specmult_planck_var_1599)
      fpl_var_1598 = ((specmult_planck_var_1599) - AINT((specmult_planck_var_1599)))
      ind0_var_1557 = ((jp_var_1533(jl_var_1611, lay_var_1564) - 1) * 5 + (jt_var_1534(jl_var_1611, lay_var_1564) - 1)) * nspa_var_216(5) + js_var_1563
      ind1_var_1558 = (jp_var_1533(jl_var_1611, lay_var_1564) * 5 + (jt1_var_1535(jl_var_1611, lay_var_1564) - 1)) * nspa_var_216(5) + js1_var_1565
      inds_var_1559 = indself_var_1543(jl_var_1611, lay_var_1564)
      indf_var_1560 = indfor_var_1549(jl_var_1611, lay_var_1564)
      indm_var_1561 = indminor_var_1553(jl_var_1611, lay_var_1564)
      IF (specparm_var_1594 .LT. 0.125D0) THEN
        p_var_1582 = fs_var_1592 - 1.0D0
        p4_var_1583 = p_var_1582 ** 4
        fk0_var_1584 = p4_var_1583
        fk1_var_1585 = 1.0D0 - p_var_1582 - 2.0D0 * p4_var_1583
        fk2_var_1586 = p_var_1582 + p4_var_1583
        fac000_var_1570 = fk0_var_1584 * fac00_var_1529(jl_var_1611, lay_var_1564)
        fac100_var_1571 = fk1_var_1585 * fac00_var_1529(jl_var_1611, lay_var_1564)
        fac200_var_1572 = fk2_var_1586 * fac00_var_1529(jl_var_1611, lay_var_1564)
        fac010_var_1573 = fk0_var_1584 * fac10_var_1531(jl_var_1611, lay_var_1564)
        fac110_var_1574 = fk1_var_1585 * fac10_var_1531(jl_var_1611, lay_var_1564)
        fac210_var_1575 = fk2_var_1586 * fac10_var_1531(jl_var_1611, lay_var_1564)
      ELSE IF (specparm_var_1594 .GT. 0.875D0) THEN
        p_var_1582 = - fs_var_1592
        p4_var_1583 = p_var_1582 ** 4
        fk0_var_1584 = p4_var_1583
        fk1_var_1585 = 1.0D0 - p_var_1582 - 2.0D0 * p4_var_1583
        fk2_var_1586 = p_var_1582 + p4_var_1583
        fac000_var_1570 = fk0_var_1584 * fac00_var_1529(jl_var_1611, lay_var_1564)
        fac100_var_1571 = fk1_var_1585 * fac00_var_1529(jl_var_1611, lay_var_1564)
        fac200_var_1572 = fk2_var_1586 * fac00_var_1529(jl_var_1611, lay_var_1564)
        fac010_var_1573 = fk0_var_1584 * fac10_var_1531(jl_var_1611, lay_var_1564)
        fac110_var_1574 = fk1_var_1585 * fac10_var_1531(jl_var_1611, lay_var_1564)
        fac210_var_1575 = fk2_var_1586 * fac10_var_1531(jl_var_1611, lay_var_1564)
      ELSE
        fac000_var_1570 = (1.0D0 - fs_var_1592) * fac00_var_1529(jl_var_1611, lay_var_1564)
        fac010_var_1573 = (1.0D0 - fs_var_1592) * fac10_var_1531(jl_var_1611, lay_var_1564)
        fac100_var_1571 = fs_var_1592 * fac00_var_1529(jl_var_1611, lay_var_1564)
        fac110_var_1574 = fs_var_1592 * fac10_var_1531(jl_var_1611, lay_var_1564)
        fac200_var_1572 = 0.0D0
        fac210_var_1575 = 0.0D0
      END IF
      IF (specparm1_var_1597 .LT. 0.125D0) THEN
        p_var_1582 = fs1_var_1595 - 1.0D0
        p4_var_1583 = p_var_1582 ** 4
        fk0_var_1584 = p4_var_1583
        fk1_var_1585 = 1.0D0 - p_var_1582 - 2.0D0 * p4_var_1583
        fk2_var_1586 = p_var_1582 + p4_var_1583
        fac001_var_1576 = fk0_var_1584 * fac01_var_1530(jl_var_1611, lay_var_1564)
        fac101_var_1577 = fk1_var_1585 * fac01_var_1530(jl_var_1611, lay_var_1564)
        fac201_var_1578 = fk2_var_1586 * fac01_var_1530(jl_var_1611, lay_var_1564)
        fac011_var_1579 = fk0_var_1584 * fac11_var_1532(jl_var_1611, lay_var_1564)
        fac111_var_1580 = fk1_var_1585 * fac11_var_1532(jl_var_1611, lay_var_1564)
        fac211_var_1581 = fk2_var_1586 * fac11_var_1532(jl_var_1611, lay_var_1564)
      ELSE IF (specparm1_var_1597 .GT. 0.875D0) THEN
        p_var_1582 = - fs1_var_1595
        p4_var_1583 = p_var_1582 ** 4
        fk0_var_1584 = p4_var_1583
        fk1_var_1585 = 1.0D0 - p_var_1582 - 2.0D0 * p4_var_1583
        fk2_var_1586 = p_var_1582 + p4_var_1583
        fac001_var_1576 = fk0_var_1584 * fac01_var_1530(jl_var_1611, lay_var_1564)
        fac101_var_1577 = fk1_var_1585 * fac01_var_1530(jl_var_1611, lay_var_1564)
        fac201_var_1578 = fk2_var_1586 * fac01_var_1530(jl_var_1611, lay_var_1564)
        fac011_var_1579 = fk0_var_1584 * fac11_var_1532(jl_var_1611, lay_var_1564)
        fac111_var_1580 = fk1_var_1585 * fac11_var_1532(jl_var_1611, lay_var_1564)
        fac211_var_1581 = fk2_var_1586 * fac11_var_1532(jl_var_1611, lay_var_1564)
      ELSE
        fac001_var_1576 = (1.0D0 - fs1_var_1595) * fac01_var_1530(jl_var_1611, lay_var_1564)
        fac011_var_1579 = (1.0D0 - fs1_var_1595) * fac11_var_1532(jl_var_1611, lay_var_1564)
        fac101_var_1577 = fs1_var_1595 * fac01_var_1530(jl_var_1611, lay_var_1564)
        fac111_var_1580 = fs1_var_1595 * fac11_var_1532(jl_var_1611, lay_var_1564)
        fac201_var_1578 = 0.0D0
        fac211_var_1581 = 0.0D0
      END IF
      IF (specparm_var_1594 .LT. 0.125D0) THEN
        tau_major_var_1589(1 : ng5) = speccomb_var_1554 * (fac000_var_1570 * absa_var_175(ind0_var_1557, 1 : 16) + fac100_var_1571 * absa_var_175(ind0_var_1557 + 1, 1 : 16) + fac200_var_1572 * absa_var_175(ind0_var_1557 + 2, 1 : 16) + fac010_var_1573 * absa_var_175(ind0_var_1557 + 9, 1 : 16) + fac110_var_1574 * absa_var_175(ind0_var_1557 + 10, 1 : 16) + fac210_var_1575 * absa_var_175(ind0_var_1557 + 11, 1 : 16))
      ELSE IF (specparm_var_1594 .GT. 0.875D0) THEN
        tau_major_var_1589(1 : ng5) = speccomb_var_1554 * (fac200_var_1572 * absa_var_175(ind0_var_1557 - 1, 1 : 16) + fac100_var_1571 * absa_var_175(ind0_var_1557, 1 : 16) + fac000_var_1570 * absa_var_175(ind0_var_1557 + 1, 1 : 16) + fac210_var_1575 * absa_var_175(ind0_var_1557 + 8, 1 : 16) + fac110_var_1574 * absa_var_175(ind0_var_1557 + 9, 1 : 16) + fac010_var_1573 * absa_var_175(ind0_var_1557 + 10, 1 : 16))
      ELSE
        tau_major_var_1589(1 : ng5) = speccomb_var_1554 * (fac000_var_1570 * absa_var_175(ind0_var_1557, 1 : 16) + fac100_var_1571 * absa_var_175(ind0_var_1557 + 1, 1 : 16) + fac010_var_1573 * absa_var_175(ind0_var_1557 + 9, 1 : 16) + fac110_var_1574 * absa_var_175(ind0_var_1557 + 10, 1 : 16))
      END IF
      IF (specparm1_var_1597 .LT. 0.125D0) THEN
        tau_major1_var_1590(1 : ng5) = speccomb1_var_1555 * (fac001_var_1576 * absa_var_175(ind1_var_1558, 1 : 16) + fac101_var_1577 * absa_var_175(ind1_var_1558 + 1, 1 : 16) + fac201_var_1578 * absa_var_175(ind1_var_1558 + 2, 1 : 16) + fac011_var_1579 * absa_var_175(ind1_var_1558 + 9, 1 : 16) + fac111_var_1580 * absa_var_175(ind1_var_1558 + 10, 1 : 16) + fac211_var_1581 * absa_var_175(ind1_var_1558 + 11, 1 : 16))
      ELSE IF (specparm1_var_1597 .GT. 0.875D0) THEN
        tau_major1_var_1590(1 : ng5) = speccomb1_var_1555 * (fac201_var_1578 * absa_var_175(ind1_var_1558 - 1, 1 : 16) + fac101_var_1577 * absa_var_175(ind1_var_1558, 1 : 16) + fac001_var_1576 * absa_var_175(ind1_var_1558 + 1, 1 : 16) + fac211_var_1581 * absa_var_175(ind1_var_1558 + 8, 1 : 16) + fac111_var_1580 * absa_var_175(ind1_var_1558 + 9, 1 : 16) + fac011_var_1579 * absa_var_175(ind1_var_1558 + 10, 1 : 16))
      ELSE
        tau_major1_var_1590(1 : ng5) = speccomb1_var_1555 * (fac001_var_1576 * absa_var_175(ind1_var_1558, 1 : 16) + fac101_var_1577 * absa_var_175(ind1_var_1558 + 1, 1 : 16) + fac011_var_1579 * absa_var_175(ind1_var_1558 + 9, 1 : 16) + fac111_var_1580 * absa_var_175(ind1_var_1558 + 10, 1 : 16))
      END IF
      DO ig_var_1562 = 1, 16
        tauself_var_1588 = selffac_var_1541(jl_var_1611, lay_var_1564) * (selfref_var_178(inds_var_1559, ig_var_1562) + selffrac_var_1542(jl_var_1611, lay_var_1564) * (selfref_var_178(inds_var_1559 + 1, ig_var_1562) - selfref_var_178(inds_var_1559, ig_var_1562)))
        taufor_var_1587 = forfac_var_1551(jl_var_1611, lay_var_1564) * (forref_var_179(indf_var_1560, ig_var_1562) + forfrac_var_1550(jl_var_1611, lay_var_1564) * (forref_var_179(indf_var_1560 + 1, ig_var_1562) - forref_var_179(indf_var_1560, ig_var_1562)))
        o3m1 = ka_mo3_var_177(jmo3, indm_var_1561, ig_var_1562) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_1561, ig_var_1562) - ka_mo3_var_177(jmo3, indm_var_1561, ig_var_1562))
        o3m2 = ka_mo3_var_177(jmo3, indm_var_1561 + 1, ig_var_1562) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_1561 + 1, ig_var_1562) - ka_mo3_var_177(jmo3, indm_var_1561 + 1, ig_var_1562))
        abso3_var_1591 = o3m1 + minorfrac_var_1552(jl_var_1611, lay_var_1564) * (o3m2 - o3m1)
        taug_var_1526(jl_var_1611, 52 + ig_var_1562, lay_var_1564) = tau_major_var_1589(ig_var_1562) + tau_major1_var_1590(ig_var_1562) + tauself_var_1588 + taufor_var_1587 + abso3_var_1591 * colo3_var_1539(jl_var_1611, lay_var_1564) + wx_var_1527(jl_var_1611, 1, lay_var_1564) * ccl4(ig_var_1562)
        fracs_var_1544(jl_var_1611, 52 + ig_var_1562, lay_var_1564) = fracrefa_var_173(ig_var_1562, jpl_var_1566) + fpl_var_1598 * (fracrefa_var_173(ig_var_1562, jpl_var_1566 + 1) - fracrefa_var_173(ig_var_1562, jpl_var_1566))
      END DO
    END DO
  END DO
  DO lay_var_1564 = laytrop_max_var_1602 + 1, klev_var_1525
    DO jl_var_1611 = kidia_var_1523, kfdia_var_1524
      speccomb_var_1554 = colo3_var_1539(jl_var_1611, lay_var_1564) + rat_o3co2_var_1547(jl_var_1611, lay_var_1564) * colco2_var_1538(jl_var_1611, lay_var_1564)
      specparm_var_1594 = MIN(colo3_var_1539(jl_var_1611, lay_var_1564) / speccomb_var_1554, oneminus_var_1536)
      specmult_var_1593 = 4.0D0 * (specparm_var_1594)
      js_var_1563 = 1 + INT(specmult_var_1593)
      fs_var_1592 = ((specmult_var_1593) - AINT((specmult_var_1593)))
      speccomb1_var_1555 = colo3_var_1539(jl_var_1611, lay_var_1564) + rat_o3co2_1_var_1548(jl_var_1611, lay_var_1564) * colco2_var_1538(jl_var_1611, lay_var_1564)
      specparm1_var_1597 = MIN(colo3_var_1539(jl_var_1611, lay_var_1564) / speccomb1_var_1555, oneminus_var_1536)
      specmult1_var_1596 = 4.0D0 * (specparm1_var_1597)
      js1_var_1565 = 1 + INT(specmult1_var_1596)
      fs1_var_1595 = ((specmult1_var_1596) - AINT((specmult1_var_1596)))
      fac000_var_1570 = (1.0D0 - fs_var_1592) * fac00_var_1529(jl_var_1611, lay_var_1564)
      fac010_var_1573 = (1.0D0 - fs_var_1592) * fac10_var_1531(jl_var_1611, lay_var_1564)
      fac100_var_1571 = fs_var_1592 * fac00_var_1529(jl_var_1611, lay_var_1564)
      fac110_var_1574 = fs_var_1592 * fac10_var_1531(jl_var_1611, lay_var_1564)
      fac001_var_1576 = (1.0D0 - fs1_var_1595) * fac01_var_1530(jl_var_1611, lay_var_1564)
      fac011_var_1579 = (1.0D0 - fs1_var_1595) * fac11_var_1532(jl_var_1611, lay_var_1564)
      fac101_var_1577 = fs1_var_1595 * fac01_var_1530(jl_var_1611, lay_var_1564)
      fac111_var_1580 = fs1_var_1595 * fac11_var_1532(jl_var_1611, lay_var_1564)
      speccomb_planck_var_1556 = colo3_var_1539(jl_var_1611, lay_var_1564) + refrat_planck_b_var_1568 * colco2_var_1538(jl_var_1611, lay_var_1564)
      specparm_planck_var_1600 = MIN(colo3_var_1539(jl_var_1611, lay_var_1564) / speccomb_planck_var_1556, oneminus_var_1536)
      specmult_planck_var_1599 = 4.0D0 * specparm_planck_var_1600
      jpl_var_1566 = 1 + INT(specmult_planck_var_1599)
      fpl_var_1598 = ((specmult_planck_var_1599) - AINT((specmult_planck_var_1599)))
      ind0_var_1557 = ((jp_var_1533(jl_var_1611, lay_var_1564) - 13) * 5 + (jt_var_1534(jl_var_1611, lay_var_1564) - 1)) * nspb_var_217(5) + js_var_1563
      ind1_var_1558 = ((jp_var_1533(jl_var_1611, lay_var_1564) - 12) * 5 + (jt1_var_1535(jl_var_1611, lay_var_1564) - 1)) * nspb_var_217(5) + js1_var_1565
      DO ig_var_1562 = 1, 16
        taug_var_1526(jl_var_1611, 52 + ig_var_1562, lay_var_1564) = speccomb_var_1554 * (fac000_var_1570 * absb_var_176(ind0_var_1557, ig_var_1562) + fac100_var_1571 * absb_var_176(ind0_var_1557 + 1, ig_var_1562) + fac010_var_1573 * absb_var_176(ind0_var_1557 + 5, ig_var_1562) + fac110_var_1574 * absb_var_176(ind0_var_1557 + 6, ig_var_1562)) + speccomb1_var_1555 * (fac001_var_1576 * absb_var_176(ind1_var_1558, ig_var_1562) + fac101_var_1577 * absb_var_176(ind1_var_1558 + 1, ig_var_1562) + fac011_var_1579 * absb_var_176(ind1_var_1558 + 5, ig_var_1562) + fac111_var_1580 * absb_var_176(ind1_var_1558 + 6, ig_var_1562)) + wx_var_1527(jl_var_1611, 1, lay_var_1564) * ccl4(ig_var_1562)
        fracs_var_1544(jl_var_1611, 52 + ig_var_1562, lay_var_1564) = fracrefb_var_174(ig_var_1562, jpl_var_1566) + fpl_var_1598 * (fracrefb_var_174(ig_var_1562, jpl_var_1566 + 1) - fracrefb_var_174(ig_var_1562, jpl_var_1566))
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1602 /= laytrop_min_var_1601) THEN
    DO lay_var_1564 = laytrop_min_var_1601 + 1, laytrop_max_var_1602
      ixc0_var_1608 = ixc_var_1603(lay_var_1564)
      DO ixp_var_1609 = 1, ixc0_var_1608
        jl_var_1611 = ixlow_var_1604(ixp_var_1609, lay_var_1564)
        speccomb_var_1554 = colh2o_var_1537(jl_var_1611, lay_var_1564) + rat_h2oco2_var_1545(jl_var_1611, lay_var_1564) * colco2_var_1538(jl_var_1611, lay_var_1564)
        specparm_var_1594 = MIN(colh2o_var_1537(jl_var_1611, lay_var_1564) / speccomb_var_1554, oneminus_var_1536)
        specmult_var_1593 = 8.0D0 * (specparm_var_1594)
        js_var_1563 = 1 + INT(specmult_var_1593)
        fs_var_1592 = ((specmult_var_1593) - AINT((specmult_var_1593)))
        speccomb1_var_1555 = colh2o_var_1537(jl_var_1611, lay_var_1564) + rat_h2oco2_1_var_1546(jl_var_1611, lay_var_1564) * colco2_var_1538(jl_var_1611, lay_var_1564)
        specparm1_var_1597 = MIN(colh2o_var_1537(jl_var_1611, lay_var_1564) / speccomb1_var_1555, oneminus_var_1536)
        specmult1_var_1596 = 8.0D0 * (specparm1_var_1597)
        js1_var_1565 = 1 + INT(specmult1_var_1596)
        fs1_var_1595 = ((specmult1_var_1596) - AINT((specmult1_var_1596)))
        speccomb_mo3 = colh2o_var_1537(jl_var_1611, lay_var_1564) + refrat_m_a_var_1569 * colco2_var_1538(jl_var_1611, lay_var_1564)
        specparm_mo3 = MIN(colh2o_var_1537(jl_var_1611, lay_var_1564) / speccomb_mo3, oneminus_var_1536)
        specmult_mo3 = 8.0D0 * specparm_mo3
        jmo3 = 1 + INT(specmult_mo3)
        fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))
        speccomb_planck_var_1556 = colh2o_var_1537(jl_var_1611, lay_var_1564) + refrat_planck_a_var_1567 * colco2_var_1538(jl_var_1611, lay_var_1564)
        specparm_planck_var_1600 = MIN(colh2o_var_1537(jl_var_1611, lay_var_1564) / speccomb_planck_var_1556, oneminus_var_1536)
        specmult_planck_var_1599 = 8.0D0 * specparm_planck_var_1600
        jpl_var_1566 = 1 + INT(specmult_planck_var_1599)
        fpl_var_1598 = ((specmult_planck_var_1599) - AINT((specmult_planck_var_1599)))
        ind0_var_1557 = ((jp_var_1533(jl_var_1611, lay_var_1564) - 1) * 5 + (jt_var_1534(jl_var_1611, lay_var_1564) - 1)) * nspa_var_216(5) + js_var_1563
        ind1_var_1558 = (jp_var_1533(jl_var_1611, lay_var_1564) * 5 + (jt1_var_1535(jl_var_1611, lay_var_1564) - 1)) * nspa_var_216(5) + js1_var_1565
        inds_var_1559 = indself_var_1543(jl_var_1611, lay_var_1564)
        indf_var_1560 = indfor_var_1549(jl_var_1611, lay_var_1564)
        indm_var_1561 = indminor_var_1553(jl_var_1611, lay_var_1564)
        IF (specparm_var_1594 .LT. 0.125D0) THEN
          p_var_1582 = fs_var_1592 - 1.0D0
          p4_var_1583 = p_var_1582 ** 4
          fk0_var_1584 = p4_var_1583
          fk1_var_1585 = 1.0D0 - p_var_1582 - 2.0D0 * p4_var_1583
          fk2_var_1586 = p_var_1582 + p4_var_1583
          fac000_var_1570 = fk0_var_1584 * fac00_var_1529(jl_var_1611, lay_var_1564)
          fac100_var_1571 = fk1_var_1585 * fac00_var_1529(jl_var_1611, lay_var_1564)
          fac200_var_1572 = fk2_var_1586 * fac00_var_1529(jl_var_1611, lay_var_1564)
          fac010_var_1573 = fk0_var_1584 * fac10_var_1531(jl_var_1611, lay_var_1564)
          fac110_var_1574 = fk1_var_1585 * fac10_var_1531(jl_var_1611, lay_var_1564)
          fac210_var_1575 = fk2_var_1586 * fac10_var_1531(jl_var_1611, lay_var_1564)
        ELSE IF (specparm_var_1594 .GT. 0.875D0) THEN
          p_var_1582 = - fs_var_1592
          p4_var_1583 = p_var_1582 ** 4
          fk0_var_1584 = p4_var_1583
          fk1_var_1585 = 1.0D0 - p_var_1582 - 2.0D0 * p4_var_1583
          fk2_var_1586 = p_var_1582 + p4_var_1583
          fac000_var_1570 = fk0_var_1584 * fac00_var_1529(jl_var_1611, lay_var_1564)
          fac100_var_1571 = fk1_var_1585 * fac00_var_1529(jl_var_1611, lay_var_1564)
          fac200_var_1572 = fk2_var_1586 * fac00_var_1529(jl_var_1611, lay_var_1564)
          fac010_var_1573 = fk0_var_1584 * fac10_var_1531(jl_var_1611, lay_var_1564)
          fac110_var_1574 = fk1_var_1585 * fac10_var_1531(jl_var_1611, lay_var_1564)
          fac210_var_1575 = fk2_var_1586 * fac10_var_1531(jl_var_1611, lay_var_1564)
        ELSE
          fac000_var_1570 = (1.0D0 - fs_var_1592) * fac00_var_1529(jl_var_1611, lay_var_1564)
          fac010_var_1573 = (1.0D0 - fs_var_1592) * fac10_var_1531(jl_var_1611, lay_var_1564)
          fac100_var_1571 = fs_var_1592 * fac00_var_1529(jl_var_1611, lay_var_1564)
          fac110_var_1574 = fs_var_1592 * fac10_var_1531(jl_var_1611, lay_var_1564)
          fac200_var_1572 = 0.0D0
          fac210_var_1575 = 0.0D0
        END IF
        IF (specparm1_var_1597 .LT. 0.125D0) THEN
          p_var_1582 = fs1_var_1595 - 1.0D0
          p4_var_1583 = p_var_1582 ** 4
          fk0_var_1584 = p4_var_1583
          fk1_var_1585 = 1.0D0 - p_var_1582 - 2.0D0 * p4_var_1583
          fk2_var_1586 = p_var_1582 + p4_var_1583
          fac001_var_1576 = fk0_var_1584 * fac01_var_1530(jl_var_1611, lay_var_1564)
          fac101_var_1577 = fk1_var_1585 * fac01_var_1530(jl_var_1611, lay_var_1564)
          fac201_var_1578 = fk2_var_1586 * fac01_var_1530(jl_var_1611, lay_var_1564)
          fac011_var_1579 = fk0_var_1584 * fac11_var_1532(jl_var_1611, lay_var_1564)
          fac111_var_1580 = fk1_var_1585 * fac11_var_1532(jl_var_1611, lay_var_1564)
          fac211_var_1581 = fk2_var_1586 * fac11_var_1532(jl_var_1611, lay_var_1564)
        ELSE IF (specparm1_var_1597 .GT. 0.875D0) THEN
          p_var_1582 = - fs1_var_1595
          p4_var_1583 = p_var_1582 ** 4
          fk0_var_1584 = p4_var_1583
          fk1_var_1585 = 1.0D0 - p_var_1582 - 2.0D0 * p4_var_1583
          fk2_var_1586 = p_var_1582 + p4_var_1583
          fac001_var_1576 = fk0_var_1584 * fac01_var_1530(jl_var_1611, lay_var_1564)
          fac101_var_1577 = fk1_var_1585 * fac01_var_1530(jl_var_1611, lay_var_1564)
          fac201_var_1578 = fk2_var_1586 * fac01_var_1530(jl_var_1611, lay_var_1564)
          fac011_var_1579 = fk0_var_1584 * fac11_var_1532(jl_var_1611, lay_var_1564)
          fac111_var_1580 = fk1_var_1585 * fac11_var_1532(jl_var_1611, lay_var_1564)
          fac211_var_1581 = fk2_var_1586 * fac11_var_1532(jl_var_1611, lay_var_1564)
        ELSE
          fac001_var_1576 = (1.0D0 - fs1_var_1595) * fac01_var_1530(jl_var_1611, lay_var_1564)
          fac011_var_1579 = (1.0D0 - fs1_var_1595) * fac11_var_1532(jl_var_1611, lay_var_1564)
          fac101_var_1577 = fs1_var_1595 * fac01_var_1530(jl_var_1611, lay_var_1564)
          fac111_var_1580 = fs1_var_1595 * fac11_var_1532(jl_var_1611, lay_var_1564)
          fac201_var_1578 = 0.0D0
          fac211_var_1581 = 0.0D0
        END IF
        IF (specparm_var_1594 .LT. 0.125D0) THEN
          tau_major_var_1589(1 : ng5) = speccomb_var_1554 * (fac000_var_1570 * absa_var_175(ind0_var_1557, 1 : 16) + fac100_var_1571 * absa_var_175(ind0_var_1557 + 1, 1 : 16) + fac200_var_1572 * absa_var_175(ind0_var_1557 + 2, 1 : 16) + fac010_var_1573 * absa_var_175(ind0_var_1557 + 9, 1 : 16) + fac110_var_1574 * absa_var_175(ind0_var_1557 + 10, 1 : 16) + fac210_var_1575 * absa_var_175(ind0_var_1557 + 11, 1 : 16))
        ELSE IF (specparm_var_1594 .GT. 0.875D0) THEN
          tau_major_var_1589(1 : ng5) = speccomb_var_1554 * (fac200_var_1572 * absa_var_175(ind0_var_1557 - 1, 1 : 16) + fac100_var_1571 * absa_var_175(ind0_var_1557, 1 : 16) + fac000_var_1570 * absa_var_175(ind0_var_1557 + 1, 1 : 16) + fac210_var_1575 * absa_var_175(ind0_var_1557 + 8, 1 : 16) + fac110_var_1574 * absa_var_175(ind0_var_1557 + 9, 1 : 16) + fac010_var_1573 * absa_var_175(ind0_var_1557 + 10, 1 : 16))
        ELSE
          tau_major_var_1589(1 : ng5) = speccomb_var_1554 * (fac000_var_1570 * absa_var_175(ind0_var_1557, 1 : 16) + fac100_var_1571 * absa_var_175(ind0_var_1557 + 1, 1 : 16) + fac010_var_1573 * absa_var_175(ind0_var_1557 + 9, 1 : 16) + fac110_var_1574 * absa_var_175(ind0_var_1557 + 10, 1 : 16))
        END IF
        IF (specparm1_var_1597 .LT. 0.125D0) THEN
          tau_major1_var_1590(1 : ng5) = speccomb1_var_1555 * (fac001_var_1576 * absa_var_175(ind1_var_1558, 1 : 16) + fac101_var_1577 * absa_var_175(ind1_var_1558 + 1, 1 : 16) + fac201_var_1578 * absa_var_175(ind1_var_1558 + 2, 1 : 16) + fac011_var_1579 * absa_var_175(ind1_var_1558 + 9, 1 : 16) + fac111_var_1580 * absa_var_175(ind1_var_1558 + 10, 1 : 16) + fac211_var_1581 * absa_var_175(ind1_var_1558 + 11, 1 : 16))
        ELSE IF (specparm1_var_1597 .GT. 0.875D0) THEN
          tau_major1_var_1590(1 : ng5) = speccomb1_var_1555 * (fac201_var_1578 * absa_var_175(ind1_var_1558 - 1, 1 : 16) + fac101_var_1577 * absa_var_175(ind1_var_1558, 1 : 16) + fac001_var_1576 * absa_var_175(ind1_var_1558 + 1, 1 : 16) + fac211_var_1581 * absa_var_175(ind1_var_1558 + 8, 1 : 16) + fac111_var_1580 * absa_var_175(ind1_var_1558 + 9, 1 : 16) + fac011_var_1579 * absa_var_175(ind1_var_1558 + 10, 1 : 16))
        ELSE
          tau_major1_var_1590(1 : ng5) = speccomb1_var_1555 * (fac001_var_1576 * absa_var_175(ind1_var_1558, 1 : 16) + fac101_var_1577 * absa_var_175(ind1_var_1558 + 1, 1 : 16) + fac011_var_1579 * absa_var_175(ind1_var_1558 + 9, 1 : 16) + fac111_var_1580 * absa_var_175(ind1_var_1558 + 10, 1 : 16))
        END IF
        DO ig_var_1562 = 1, 16
          tauself_var_1588 = selffac_var_1541(jl_var_1611, lay_var_1564) * (selfref_var_178(inds_var_1559, ig_var_1562) + selffrac_var_1542(jl_var_1611, lay_var_1564) * (selfref_var_178(inds_var_1559 + 1, ig_var_1562) - selfref_var_178(inds_var_1559, ig_var_1562)))
          taufor_var_1587 = forfac_var_1551(jl_var_1611, lay_var_1564) * (forref_var_179(indf_var_1560, ig_var_1562) + forfrac_var_1550(jl_var_1611, lay_var_1564) * (forref_var_179(indf_var_1560 + 1, ig_var_1562) - forref_var_179(indf_var_1560, ig_var_1562)))
          o3m1 = ka_mo3_var_177(jmo3, indm_var_1561, ig_var_1562) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_1561, ig_var_1562) - ka_mo3_var_177(jmo3, indm_var_1561, ig_var_1562))
          o3m2 = ka_mo3_var_177(jmo3, indm_var_1561 + 1, ig_var_1562) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_1561 + 1, ig_var_1562) - ka_mo3_var_177(jmo3, indm_var_1561 + 1, ig_var_1562))
          abso3_var_1591 = o3m1 + minorfrac_var_1552(jl_var_1611, lay_var_1564) * (o3m2 - o3m1)
          taug_var_1526(jl_var_1611, 52 + ig_var_1562, lay_var_1564) = tau_major_var_1589(ig_var_1562) + tau_major1_var_1590(ig_var_1562) + tauself_var_1588 + taufor_var_1587 + abso3_var_1591 * colo3_var_1539(jl_var_1611, lay_var_1564) + wx_var_1527(jl_var_1611, 1, lay_var_1564) * ccl4(ig_var_1562)
          fracs_var_1544(jl_var_1611, 52 + ig_var_1562, lay_var_1564) = fracrefa_var_173(ig_var_1562, jpl_var_1566) + fpl_var_1598 * (fracrefa_var_173(ig_var_1562, jpl_var_1566 + 1) - fracrefa_var_173(ig_var_1562, jpl_var_1566))
        END DO
      END DO
      ixc0_var_1608 = kfdia_var_1524 - kidia_var_1523 + 1 - ixc0_var_1608
      DO ixp_var_1609 = 1, ixc0_var_1608
        jl_var_1611 = ixhigh_var_1605(ixp_var_1609, lay_var_1564)
        speccomb_var_1554 = colo3_var_1539(jl_var_1611, lay_var_1564) + rat_o3co2_var_1547(jl_var_1611, lay_var_1564) * colco2_var_1538(jl_var_1611, lay_var_1564)
        specparm_var_1594 = MIN(colo3_var_1539(jl_var_1611, lay_var_1564) / speccomb_var_1554, oneminus_var_1536)
        specmult_var_1593 = 4.0D0 * (specparm_var_1594)
        js_var_1563 = 1 + INT(specmult_var_1593)
        fs_var_1592 = ((specmult_var_1593) - AINT((specmult_var_1593)))
        speccomb1_var_1555 = colo3_var_1539(jl_var_1611, lay_var_1564) + rat_o3co2_1_var_1548(jl_var_1611, lay_var_1564) * colco2_var_1538(jl_var_1611, lay_var_1564)
        specparm1_var_1597 = MIN(colo3_var_1539(jl_var_1611, lay_var_1564) / speccomb1_var_1555, oneminus_var_1536)
        specmult1_var_1596 = 4.0D0 * (specparm1_var_1597)
        js1_var_1565 = 1 + INT(specmult1_var_1596)
        fs1_var_1595 = ((specmult1_var_1596) - AINT((specmult1_var_1596)))
        fac000_var_1570 = (1.0D0 - fs_var_1592) * fac00_var_1529(jl_var_1611, lay_var_1564)
        fac010_var_1573 = (1.0D0 - fs_var_1592) * fac10_var_1531(jl_var_1611, lay_var_1564)
        fac100_var_1571 = fs_var_1592 * fac00_var_1529(jl_var_1611, lay_var_1564)
        fac110_var_1574 = fs_var_1592 * fac10_var_1531(jl_var_1611, lay_var_1564)
        fac001_var_1576 = (1.0D0 - fs1_var_1595) * fac01_var_1530(jl_var_1611, lay_var_1564)
        fac011_var_1579 = (1.0D0 - fs1_var_1595) * fac11_var_1532(jl_var_1611, lay_var_1564)
        fac101_var_1577 = fs1_var_1595 * fac01_var_1530(jl_var_1611, lay_var_1564)
        fac111_var_1580 = fs1_var_1595 * fac11_var_1532(jl_var_1611, lay_var_1564)
        speccomb_planck_var_1556 = colo3_var_1539(jl_var_1611, lay_var_1564) + refrat_planck_b_var_1568 * colco2_var_1538(jl_var_1611, lay_var_1564)
        specparm_planck_var_1600 = MIN(colo3_var_1539(jl_var_1611, lay_var_1564) / speccomb_planck_var_1556, oneminus_var_1536)
        specmult_planck_var_1599 = 4.0D0 * specparm_planck_var_1600
        jpl_var_1566 = 1 + INT(specmult_planck_var_1599)
        fpl_var_1598 = ((specmult_planck_var_1599) - AINT((specmult_planck_var_1599)))
        ind0_var_1557 = ((jp_var_1533(jl_var_1611, lay_var_1564) - 13) * 5 + (jt_var_1534(jl_var_1611, lay_var_1564) - 1)) * nspb_var_217(5) + js_var_1563
        ind1_var_1558 = ((jp_var_1533(jl_var_1611, lay_var_1564) - 12) * 5 + (jt1_var_1535(jl_var_1611, lay_var_1564) - 1)) * nspb_var_217(5) + js1_var_1565
        DO ig_var_1562 = 1, 16
          taug_var_1526(jl_var_1611, 52 + ig_var_1562, lay_var_1564) = speccomb_var_1554 * (fac000_var_1570 * absb_var_176(ind0_var_1557, ig_var_1562) + fac100_var_1571 * absb_var_176(ind0_var_1557 + 1, ig_var_1562) + fac010_var_1573 * absb_var_176(ind0_var_1557 + 5, ig_var_1562) + fac110_var_1574 * absb_var_176(ind0_var_1557 + 6, ig_var_1562)) + speccomb1_var_1555 * (fac001_var_1576 * absb_var_176(ind1_var_1558, ig_var_1562) + fac101_var_1577 * absb_var_176(ind1_var_1558 + 1, ig_var_1562) + fac011_var_1579 * absb_var_176(ind1_var_1558 + 5, ig_var_1562) + fac111_var_1580 * absb_var_176(ind1_var_1558 + 6, ig_var_1562)) + wx_var_1527(jl_var_1611, 1, lay_var_1564) * ccl4(ig_var_1562)
          fracs_var_1544(jl_var_1611, 52 + ig_var_1562, lay_var_1564) = fracrefb_var_174(ig_var_1562, jpl_var_1566) + fpl_var_1598 * (fracrefb_var_174(ig_var_1562, jpl_var_1566 + 1) - fracrefb_var_174(ig_var_1562, jpl_var_1566))
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol5
SUBROUTINE rrtm_taumol4(kidia_var_1612, kfdia_var_1613, klev_var_1614, taug_var_1615, p_tauaerl_var_1616, fac00_var_1617, fac01_var_1618, fac10_var_1619, fac11_var_1620, forfac_var_1638, forfrac_var_1639, indfor_var_1637, jp_var_1621, jt_var_1622, jt1_var_1623, oneminus_var_1624, colh2o_var_1625, colco2_var_1626, colo3_var_1627, laytrop_var_1628, selffac_var_1629, selffrac_var_1630, indself_var_1631, fracs_var_1632, rat_h2oco2_var_1633, rat_h2oco2_1_var_1634, rat_o3co2_var_1635, rat_o3co2_1_var_1636)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng4
  USE yoerrta4, ONLY: absa_var_169, absb_var_170, forref_var_172, fracrefa_var_167, fracrefb_var_168, selfref_var_171
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1612
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1613
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1614
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1615(kidia_var_1612 : kfdia_var_1613, 140, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1616(kidia_var_1612 : kfdia_var_1613, klev_var_1614, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1617(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1618(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1619(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1620(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1621(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1622(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1623(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1624
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1625(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1626(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1627(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1628(kidia_var_1612 : kfdia_var_1613)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1629(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1630(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1631(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1632(kidia_var_1612 : kfdia_var_1613, 140, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1633(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1634(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_var_1635(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_1_var_1636(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1637(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1638(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1639(kidia_var_1612 : kfdia_var_1613, klev_var_1614)
  REAL(KIND = 8) :: speccomb_var_1640, speccomb1_var_1641, speccomb_planck_var_1642
  INTEGER(KIND = 4) :: ind0_var_1643, ind1_var_1644, inds_var_1645, indf_var_1646
  INTEGER(KIND = 4) :: ig_var_1647, js_var_1648, lay_var_1649, js1_var_1650, jpl_var_1651
  REAL(KIND = 8) :: refrat_planck_a_var_1652, refrat_planck_b_var_1653
  REAL(KIND = 8) :: fac000_var_1654, fac100_var_1655, fac200_var_1656, fac010_var_1657, fac110_var_1658, fac210_var_1659, fac001_var_1660, fac101_var_1661, fac201_var_1662, fac011_var_1663, fac111_var_1664, fac211_var_1665
  REAL(KIND = 8) :: p_var_1666, p4_var_1667, fk0_var_1668, fk1_var_1669, fk2_var_1670
  REAL(KIND = 8) :: taufor_var_1671, tauself_var_1672, tau_major_var_1673(14), tau_major1_var_1674(14)
  REAL(KIND = 8) :: fs_var_1675, specmult_var_1676, specparm_var_1677, fs1_var_1678, specmult1_var_1679, specparm1_var_1680, fpl_var_1681, specmult_planck_var_1682, specparm_planck_var_1683
  INTEGER(KIND = 4) :: laytrop_min_var_1684, laytrop_max_var_1685
  INTEGER(KIND = 4) :: ixc_var_1686(klev_var_1614), ixlow_var_1687(kfdia_var_1613, klev_var_1614), ixhigh_var_1688(kfdia_var_1613, klev_var_1614)
  INTEGER(KIND = 4) :: ich_var_1689, icl_var_1690, ixc0_var_1691, ixp_var_1692, jc_var_1693, jl_var_1694
  laytrop_min_var_1684 = MINVAL(laytrop_var_1628)
  laytrop_max_var_1685 = MAXVAL(laytrop_var_1628)
  ixlow_var_1687 = 0
  ixhigh_var_1688 = 0
  ixc_var_1686 = 0
  DO lay_var_1649 = laytrop_min_var_1684 + 1, laytrop_max_var_1685
    icl_var_1690 = 0
    ich_var_1689 = 0
    DO jc_var_1693 = kidia_var_1612, kfdia_var_1613
      IF (lay_var_1649 <= laytrop_var_1628(jc_var_1693)) THEN
        icl_var_1690 = icl_var_1690 + 1
        ixlow_var_1687(icl_var_1690, lay_var_1649) = jc_var_1693
      ELSE
        ich_var_1689 = ich_var_1689 + 1
        ixhigh_var_1688(ich_var_1689, lay_var_1649) = jc_var_1693
      END IF
    END DO
    ixc_var_1686(lay_var_1649) = icl_var_1690
  END DO
  refrat_planck_a_var_1652 = chi_mls(1, 11) / chi_mls(2, 11)
  refrat_planck_b_var_1653 = chi_mls(3, 13) / chi_mls(2, 13)
  DO lay_var_1649 = 1, laytrop_min_var_1684
    DO jl_var_1694 = kidia_var_1612, kfdia_var_1613
      speccomb_var_1640 = colh2o_var_1625(jl_var_1694, lay_var_1649) + rat_h2oco2_var_1633(jl_var_1694, lay_var_1649) * colco2_var_1626(jl_var_1694, lay_var_1649)
      specparm_var_1677 = MIN(colh2o_var_1625(jl_var_1694, lay_var_1649) / speccomb_var_1640, oneminus_var_1624)
      specmult_var_1676 = 8.0D0 * (specparm_var_1677)
      js_var_1648 = 1 + INT(specmult_var_1676)
      fs_var_1675 = ((specmult_var_1676) - AINT((specmult_var_1676)))
      speccomb1_var_1641 = colh2o_var_1625(jl_var_1694, lay_var_1649) + rat_h2oco2_1_var_1634(jl_var_1694, lay_var_1649) * colco2_var_1626(jl_var_1694, lay_var_1649)
      specparm1_var_1680 = MIN(colh2o_var_1625(jl_var_1694, lay_var_1649) / speccomb1_var_1641, oneminus_var_1624)
      specmult1_var_1679 = 8.0D0 * (specparm1_var_1680)
      js1_var_1650 = 1 + INT(specmult1_var_1679)
      fs1_var_1678 = ((specmult1_var_1679) - AINT((specmult1_var_1679)))
      speccomb_planck_var_1642 = colh2o_var_1625(jl_var_1694, lay_var_1649) + refrat_planck_a_var_1652 * colco2_var_1626(jl_var_1694, lay_var_1649)
      specparm_planck_var_1683 = MIN(colh2o_var_1625(jl_var_1694, lay_var_1649) / speccomb_planck_var_1642, oneminus_var_1624)
      specmult_planck_var_1682 = 8.0D0 * specparm_planck_var_1683
      jpl_var_1651 = 1 + INT(specmult_planck_var_1682)
      fpl_var_1681 = ((specmult_planck_var_1682) - AINT((specmult_planck_var_1682)))
      ind0_var_1643 = ((jp_var_1621(jl_var_1694, lay_var_1649) - 1) * 5 + (jt_var_1622(jl_var_1694, lay_var_1649) - 1)) * nspa_var_216(4) + js_var_1648
      ind1_var_1644 = (jp_var_1621(jl_var_1694, lay_var_1649) * 5 + (jt1_var_1623(jl_var_1694, lay_var_1649) - 1)) * nspa_var_216(4) + js1_var_1650
      inds_var_1645 = indself_var_1631(jl_var_1694, lay_var_1649)
      indf_var_1646 = indfor_var_1637(jl_var_1694, lay_var_1649)
      IF (specparm_var_1677 .LT. 0.125D0) THEN
        p_var_1666 = fs_var_1675 - 1.0D0
        p4_var_1667 = p_var_1666 ** 4
        fk0_var_1668 = p4_var_1667
        fk1_var_1669 = 1.0D0 - p_var_1666 - 2.0D0 * p4_var_1667
        fk2_var_1670 = p_var_1666 + p4_var_1667
        fac000_var_1654 = fk0_var_1668 * fac00_var_1617(jl_var_1694, lay_var_1649)
        fac100_var_1655 = fk1_var_1669 * fac00_var_1617(jl_var_1694, lay_var_1649)
        fac200_var_1656 = fk2_var_1670 * fac00_var_1617(jl_var_1694, lay_var_1649)
        fac010_var_1657 = fk0_var_1668 * fac10_var_1619(jl_var_1694, lay_var_1649)
        fac110_var_1658 = fk1_var_1669 * fac10_var_1619(jl_var_1694, lay_var_1649)
        fac210_var_1659 = fk2_var_1670 * fac10_var_1619(jl_var_1694, lay_var_1649)
      ELSE IF (specparm_var_1677 .GT. 0.875D0) THEN
        p_var_1666 = - fs_var_1675
        p4_var_1667 = p_var_1666 ** 4
        fk0_var_1668 = p4_var_1667
        fk1_var_1669 = 1.0D0 - p_var_1666 - 2.0D0 * p4_var_1667
        fk2_var_1670 = p_var_1666 + p4_var_1667
        fac000_var_1654 = fk0_var_1668 * fac00_var_1617(jl_var_1694, lay_var_1649)
        fac100_var_1655 = fk1_var_1669 * fac00_var_1617(jl_var_1694, lay_var_1649)
        fac200_var_1656 = fk2_var_1670 * fac00_var_1617(jl_var_1694, lay_var_1649)
        fac010_var_1657 = fk0_var_1668 * fac10_var_1619(jl_var_1694, lay_var_1649)
        fac110_var_1658 = fk1_var_1669 * fac10_var_1619(jl_var_1694, lay_var_1649)
        fac210_var_1659 = fk2_var_1670 * fac10_var_1619(jl_var_1694, lay_var_1649)
      ELSE
        fac000_var_1654 = (1.0D0 - fs_var_1675) * fac00_var_1617(jl_var_1694, lay_var_1649)
        fac010_var_1657 = (1.0D0 - fs_var_1675) * fac10_var_1619(jl_var_1694, lay_var_1649)
        fac100_var_1655 = fs_var_1675 * fac00_var_1617(jl_var_1694, lay_var_1649)
        fac110_var_1658 = fs_var_1675 * fac10_var_1619(jl_var_1694, lay_var_1649)
        fac200_var_1656 = 0.0D0
        fac210_var_1659 = 0.0D0
      END IF
      IF (specparm1_var_1680 .LT. 0.125D0) THEN
        p_var_1666 = fs1_var_1678 - 1.0D0
        p4_var_1667 = p_var_1666 ** 4
        fk0_var_1668 = p4_var_1667
        fk1_var_1669 = 1.0D0 - p_var_1666 - 2.0D0 * p4_var_1667
        fk2_var_1670 = p_var_1666 + p4_var_1667
        fac001_var_1660 = fk0_var_1668 * fac01_var_1618(jl_var_1694, lay_var_1649)
        fac101_var_1661 = fk1_var_1669 * fac01_var_1618(jl_var_1694, lay_var_1649)
        fac201_var_1662 = fk2_var_1670 * fac01_var_1618(jl_var_1694, lay_var_1649)
        fac011_var_1663 = fk0_var_1668 * fac11_var_1620(jl_var_1694, lay_var_1649)
        fac111_var_1664 = fk1_var_1669 * fac11_var_1620(jl_var_1694, lay_var_1649)
        fac211_var_1665 = fk2_var_1670 * fac11_var_1620(jl_var_1694, lay_var_1649)
      ELSE IF (specparm1_var_1680 .GT. 0.875D0) THEN
        p_var_1666 = - fs1_var_1678
        p4_var_1667 = p_var_1666 ** 4
        fk0_var_1668 = p4_var_1667
        fk1_var_1669 = 1.0D0 - p_var_1666 - 2.0D0 * p4_var_1667
        fk2_var_1670 = p_var_1666 + p4_var_1667
        fac001_var_1660 = fk0_var_1668 * fac01_var_1618(jl_var_1694, lay_var_1649)
        fac101_var_1661 = fk1_var_1669 * fac01_var_1618(jl_var_1694, lay_var_1649)
        fac201_var_1662 = fk2_var_1670 * fac01_var_1618(jl_var_1694, lay_var_1649)
        fac011_var_1663 = fk0_var_1668 * fac11_var_1620(jl_var_1694, lay_var_1649)
        fac111_var_1664 = fk1_var_1669 * fac11_var_1620(jl_var_1694, lay_var_1649)
        fac211_var_1665 = fk2_var_1670 * fac11_var_1620(jl_var_1694, lay_var_1649)
      ELSE
        fac001_var_1660 = (1.0D0 - fs1_var_1678) * fac01_var_1618(jl_var_1694, lay_var_1649)
        fac011_var_1663 = (1.0D0 - fs1_var_1678) * fac11_var_1620(jl_var_1694, lay_var_1649)
        fac101_var_1661 = fs1_var_1678 * fac01_var_1618(jl_var_1694, lay_var_1649)
        fac111_var_1664 = fs1_var_1678 * fac11_var_1620(jl_var_1694, lay_var_1649)
        fac201_var_1662 = 0.0D0
        fac211_var_1665 = 0.0D0
      END IF
      IF (specparm_var_1677 .LT. 0.125D0) THEN
        tau_major_var_1673(1 : ng4) = speccomb_var_1640 * (fac000_var_1654 * absa_var_169(ind0_var_1643, 1 : 14) + fac100_var_1655 * absa_var_169(ind0_var_1643 + 1, 1 : 14) + fac200_var_1656 * absa_var_169(ind0_var_1643 + 2, 1 : 14) + fac010_var_1657 * absa_var_169(ind0_var_1643 + 9, 1 : 14) + fac110_var_1658 * absa_var_169(ind0_var_1643 + 10, 1 : 14) + fac210_var_1659 * absa_var_169(ind0_var_1643 + 11, 1 : 14))
      ELSE IF (specparm_var_1677 .GT. 0.875D0) THEN
        tau_major_var_1673(1 : ng4) = speccomb_var_1640 * (fac200_var_1656 * absa_var_169(ind0_var_1643 - 1, 1 : 14) + fac100_var_1655 * absa_var_169(ind0_var_1643, 1 : 14) + fac000_var_1654 * absa_var_169(ind0_var_1643 + 1, 1 : 14) + fac210_var_1659 * absa_var_169(ind0_var_1643 + 8, 1 : 14) + fac110_var_1658 * absa_var_169(ind0_var_1643 + 9, 1 : 14) + fac010_var_1657 * absa_var_169(ind0_var_1643 + 10, 1 : 14))
      ELSE
        tau_major_var_1673(1 : ng4) = speccomb_var_1640 * (fac000_var_1654 * absa_var_169(ind0_var_1643, 1 : 14) + fac100_var_1655 * absa_var_169(ind0_var_1643 + 1, 1 : 14) + fac010_var_1657 * absa_var_169(ind0_var_1643 + 9, 1 : 14) + fac110_var_1658 * absa_var_169(ind0_var_1643 + 10, 1 : 14))
      END IF
      IF (specparm1_var_1680 .LT. 0.125D0) THEN
        tau_major1_var_1674(1 : ng4) = speccomb1_var_1641 * (fac001_var_1660 * absa_var_169(ind1_var_1644, 1 : 14) + fac101_var_1661 * absa_var_169(ind1_var_1644 + 1, 1 : 14) + fac201_var_1662 * absa_var_169(ind1_var_1644 + 2, 1 : 14) + fac011_var_1663 * absa_var_169(ind1_var_1644 + 9, 1 : 14) + fac111_var_1664 * absa_var_169(ind1_var_1644 + 10, 1 : 14) + fac211_var_1665 * absa_var_169(ind1_var_1644 + 11, 1 : 14))
      ELSE IF (specparm1_var_1680 .GT. 0.875D0) THEN
        tau_major1_var_1674(1 : ng4) = speccomb1_var_1641 * (fac201_var_1662 * absa_var_169(ind1_var_1644 - 1, 1 : 14) + fac101_var_1661 * absa_var_169(ind1_var_1644, 1 : 14) + fac001_var_1660 * absa_var_169(ind1_var_1644 + 1, 1 : 14) + fac211_var_1665 * absa_var_169(ind1_var_1644 + 8, 1 : 14) + fac111_var_1664 * absa_var_169(ind1_var_1644 + 9, 1 : 14) + fac011_var_1663 * absa_var_169(ind1_var_1644 + 10, 1 : 14))
      ELSE
        tau_major1_var_1674(1 : ng4) = speccomb1_var_1641 * (fac001_var_1660 * absa_var_169(ind1_var_1644, 1 : 14) + fac101_var_1661 * absa_var_169(ind1_var_1644 + 1, 1 : 14) + fac011_var_1663 * absa_var_169(ind1_var_1644 + 9, 1 : 14) + fac111_var_1664 * absa_var_169(ind1_var_1644 + 10, 1 : 14))
      END IF
      DO ig_var_1647 = 1, 14
        tauself_var_1672 = selffac_var_1629(jl_var_1694, lay_var_1649) * (selfref_var_171(inds_var_1645, ig_var_1647) + selffrac_var_1630(jl_var_1694, lay_var_1649) * (selfref_var_171(inds_var_1645 + 1, ig_var_1647) - selfref_var_171(inds_var_1645, ig_var_1647)))
        taufor_var_1671 = forfac_var_1638(jl_var_1694, lay_var_1649) * (forref_var_172(indf_var_1646, ig_var_1647) + forfrac_var_1639(jl_var_1694, lay_var_1649) * (forref_var_172(indf_var_1646 + 1, ig_var_1647) - forref_var_172(indf_var_1646, ig_var_1647)))
        taug_var_1615(jl_var_1694, 38 + ig_var_1647, lay_var_1649) = tau_major_var_1673(ig_var_1647) + tau_major1_var_1674(ig_var_1647) + tauself_var_1672 + taufor_var_1671
        fracs_var_1632(jl_var_1694, 38 + ig_var_1647, lay_var_1649) = fracrefa_var_167(ig_var_1647, jpl_var_1651) + fpl_var_1681 * (fracrefa_var_167(ig_var_1647, jpl_var_1651 + 1) - fracrefa_var_167(ig_var_1647, jpl_var_1651))
      END DO
    END DO
  END DO
  DO lay_var_1649 = laytrop_max_var_1685 + 1, klev_var_1614
    DO jl_var_1694 = kidia_var_1612, kfdia_var_1613
      speccomb_var_1640 = colo3_var_1627(jl_var_1694, lay_var_1649) + rat_o3co2_var_1635(jl_var_1694, lay_var_1649) * colco2_var_1626(jl_var_1694, lay_var_1649)
      specparm_var_1677 = MIN(colo3_var_1627(jl_var_1694, lay_var_1649) / speccomb_var_1640, oneminus_var_1624)
      specmult_var_1676 = 4.0D0 * (specparm_var_1677)
      js_var_1648 = 1 + INT(specmult_var_1676)
      fs_var_1675 = ((specmult_var_1676) - AINT((specmult_var_1676)))
      speccomb1_var_1641 = colo3_var_1627(jl_var_1694, lay_var_1649) + rat_o3co2_1_var_1636(jl_var_1694, lay_var_1649) * colco2_var_1626(jl_var_1694, lay_var_1649)
      specparm1_var_1680 = MIN(colo3_var_1627(jl_var_1694, lay_var_1649) / speccomb1_var_1641, oneminus_var_1624)
      specmult1_var_1679 = 4.0D0 * (specparm1_var_1680)
      js1_var_1650 = 1 + INT(specmult1_var_1679)
      fs1_var_1678 = ((specmult1_var_1679) - AINT((specmult1_var_1679)))
      fac000_var_1654 = (1.0D0 - fs_var_1675) * fac00_var_1617(jl_var_1694, lay_var_1649)
      fac010_var_1657 = (1.0D0 - fs_var_1675) * fac10_var_1619(jl_var_1694, lay_var_1649)
      fac100_var_1655 = fs_var_1675 * fac00_var_1617(jl_var_1694, lay_var_1649)
      fac110_var_1658 = fs_var_1675 * fac10_var_1619(jl_var_1694, lay_var_1649)
      fac001_var_1660 = (1.0D0 - fs1_var_1678) * fac01_var_1618(jl_var_1694, lay_var_1649)
      fac011_var_1663 = (1.0D0 - fs1_var_1678) * fac11_var_1620(jl_var_1694, lay_var_1649)
      fac101_var_1661 = fs1_var_1678 * fac01_var_1618(jl_var_1694, lay_var_1649)
      fac111_var_1664 = fs1_var_1678 * fac11_var_1620(jl_var_1694, lay_var_1649)
      speccomb_planck_var_1642 = colo3_var_1627(jl_var_1694, lay_var_1649) + refrat_planck_b_var_1653 * colco2_var_1626(jl_var_1694, lay_var_1649)
      specparm_planck_var_1683 = MIN(colo3_var_1627(jl_var_1694, lay_var_1649) / speccomb_planck_var_1642, oneminus_var_1624)
      specmult_planck_var_1682 = 4.0D0 * specparm_planck_var_1683
      jpl_var_1651 = 1 + INT(specmult_planck_var_1682)
      fpl_var_1681 = ((specmult_planck_var_1682) - AINT((specmult_planck_var_1682)))
      ind0_var_1643 = ((jp_var_1621(jl_var_1694, lay_var_1649) - 13) * 5 + (jt_var_1622(jl_var_1694, lay_var_1649) - 1)) * nspb_var_217(4) + js_var_1648
      ind1_var_1644 = ((jp_var_1621(jl_var_1694, lay_var_1649) - 12) * 5 + (jt1_var_1623(jl_var_1694, lay_var_1649) - 1)) * nspb_var_217(4) + js1_var_1650
      DO ig_var_1647 = 1, 14
        taug_var_1615(jl_var_1694, 38 + ig_var_1647, lay_var_1649) = speccomb_var_1640 * (fac000_var_1654 * absb_var_170(ind0_var_1643, ig_var_1647) + fac100_var_1655 * absb_var_170(ind0_var_1643 + 1, ig_var_1647) + fac010_var_1657 * absb_var_170(ind0_var_1643 + 5, ig_var_1647) + fac110_var_1658 * absb_var_170(ind0_var_1643 + 6, ig_var_1647)) + speccomb1_var_1641 * (fac001_var_1660 * absb_var_170(ind1_var_1644, ig_var_1647) + fac101_var_1661 * absb_var_170(ind1_var_1644 + 1, ig_var_1647) + fac011_var_1663 * absb_var_170(ind1_var_1644 + 5, ig_var_1647) + fac111_var_1664 * absb_var_170(ind1_var_1644 + 6, ig_var_1647))
        fracs_var_1632(jl_var_1694, 38 + ig_var_1647, lay_var_1649) = fracrefb_var_168(ig_var_1647, jpl_var_1651) + fpl_var_1681 * (fracrefb_var_168(ig_var_1647, jpl_var_1651 + 1) - fracrefb_var_168(ig_var_1647, jpl_var_1651))
      END DO
    END DO
  END DO
  DO lay_var_1649 = laytrop_max_var_1685 + 1, klev_var_1614
    DO jl_var_1694 = kidia_var_1612, kfdia_var_1613
      taug_var_1615(jl_var_1694, 46, lay_var_1649) = taug_var_1615(jl_var_1694, 46, lay_var_1649) * 0.92D0
      taug_var_1615(jl_var_1694, 47, lay_var_1649) = taug_var_1615(jl_var_1694, 47, lay_var_1649) * 0.88D0
      taug_var_1615(jl_var_1694, 48, lay_var_1649) = taug_var_1615(jl_var_1694, 48, lay_var_1649) * 1.07D0
      taug_var_1615(jl_var_1694, 49, lay_var_1649) = taug_var_1615(jl_var_1694, 49, lay_var_1649) * 1.1D0
      taug_var_1615(jl_var_1694, 50, lay_var_1649) = taug_var_1615(jl_var_1694, 50, lay_var_1649) * 0.99D0
      taug_var_1615(jl_var_1694, 51, lay_var_1649) = taug_var_1615(jl_var_1694, 51, lay_var_1649) * 0.88D0
      taug_var_1615(jl_var_1694, 52, lay_var_1649) = taug_var_1615(jl_var_1694, 52, lay_var_1649) * 0.943D0
    END DO
  END DO
  IF (laytrop_max_var_1685 /= laytrop_min_var_1684) THEN
    DO lay_var_1649 = laytrop_min_var_1684 + 1, laytrop_max_var_1685
      ixc0_var_1691 = ixc_var_1686(lay_var_1649)
      DO ixp_var_1692 = 1, ixc0_var_1691
        jl_var_1694 = ixlow_var_1687(ixp_var_1692, lay_var_1649)
        speccomb_var_1640 = colh2o_var_1625(jl_var_1694, lay_var_1649) + rat_h2oco2_var_1633(jl_var_1694, lay_var_1649) * colco2_var_1626(jl_var_1694, lay_var_1649)
        specparm_var_1677 = MIN(colh2o_var_1625(jl_var_1694, lay_var_1649) / speccomb_var_1640, oneminus_var_1624)
        specmult_var_1676 = 8.0D0 * (specparm_var_1677)
        js_var_1648 = 1 + INT(specmult_var_1676)
        fs_var_1675 = ((specmult_var_1676) - AINT((specmult_var_1676)))
        speccomb1_var_1641 = colh2o_var_1625(jl_var_1694, lay_var_1649) + rat_h2oco2_1_var_1634(jl_var_1694, lay_var_1649) * colco2_var_1626(jl_var_1694, lay_var_1649)
        specparm1_var_1680 = MIN(colh2o_var_1625(jl_var_1694, lay_var_1649) / speccomb1_var_1641, oneminus_var_1624)
        specmult1_var_1679 = 8.0D0 * (specparm1_var_1680)
        js1_var_1650 = 1 + INT(specmult1_var_1679)
        fs1_var_1678 = ((specmult1_var_1679) - AINT((specmult1_var_1679)))
        speccomb_planck_var_1642 = colh2o_var_1625(jl_var_1694, lay_var_1649) + refrat_planck_a_var_1652 * colco2_var_1626(jl_var_1694, lay_var_1649)
        specparm_planck_var_1683 = MIN(colh2o_var_1625(jl_var_1694, lay_var_1649) / speccomb_planck_var_1642, oneminus_var_1624)
        specmult_planck_var_1682 = 8.0D0 * specparm_planck_var_1683
        jpl_var_1651 = 1 + INT(specmult_planck_var_1682)
        fpl_var_1681 = ((specmult_planck_var_1682) - AINT((specmult_planck_var_1682)))
        ind0_var_1643 = ((jp_var_1621(jl_var_1694, lay_var_1649) - 1) * 5 + (jt_var_1622(jl_var_1694, lay_var_1649) - 1)) * nspa_var_216(4) + js_var_1648
        ind1_var_1644 = (jp_var_1621(jl_var_1694, lay_var_1649) * 5 + (jt1_var_1623(jl_var_1694, lay_var_1649) - 1)) * nspa_var_216(4) + js1_var_1650
        inds_var_1645 = indself_var_1631(jl_var_1694, lay_var_1649)
        indf_var_1646 = indfor_var_1637(jl_var_1694, lay_var_1649)
        IF (specparm_var_1677 .LT. 0.125D0) THEN
          p_var_1666 = fs_var_1675 - 1.0D0
          p4_var_1667 = p_var_1666 ** 4
          fk0_var_1668 = p4_var_1667
          fk1_var_1669 = 1.0D0 - p_var_1666 - 2.0D0 * p4_var_1667
          fk2_var_1670 = p_var_1666 + p4_var_1667
          fac000_var_1654 = fk0_var_1668 * fac00_var_1617(jl_var_1694, lay_var_1649)
          fac100_var_1655 = fk1_var_1669 * fac00_var_1617(jl_var_1694, lay_var_1649)
          fac200_var_1656 = fk2_var_1670 * fac00_var_1617(jl_var_1694, lay_var_1649)
          fac010_var_1657 = fk0_var_1668 * fac10_var_1619(jl_var_1694, lay_var_1649)
          fac110_var_1658 = fk1_var_1669 * fac10_var_1619(jl_var_1694, lay_var_1649)
          fac210_var_1659 = fk2_var_1670 * fac10_var_1619(jl_var_1694, lay_var_1649)
        ELSE IF (specparm_var_1677 .GT. 0.875D0) THEN
          p_var_1666 = - fs_var_1675
          p4_var_1667 = p_var_1666 ** 4
          fk0_var_1668 = p4_var_1667
          fk1_var_1669 = 1.0D0 - p_var_1666 - 2.0D0 * p4_var_1667
          fk2_var_1670 = p_var_1666 + p4_var_1667
          fac000_var_1654 = fk0_var_1668 * fac00_var_1617(jl_var_1694, lay_var_1649)
          fac100_var_1655 = fk1_var_1669 * fac00_var_1617(jl_var_1694, lay_var_1649)
          fac200_var_1656 = fk2_var_1670 * fac00_var_1617(jl_var_1694, lay_var_1649)
          fac010_var_1657 = fk0_var_1668 * fac10_var_1619(jl_var_1694, lay_var_1649)
          fac110_var_1658 = fk1_var_1669 * fac10_var_1619(jl_var_1694, lay_var_1649)
          fac210_var_1659 = fk2_var_1670 * fac10_var_1619(jl_var_1694, lay_var_1649)
        ELSE
          fac000_var_1654 = (1.0D0 - fs_var_1675) * fac00_var_1617(jl_var_1694, lay_var_1649)
          fac010_var_1657 = (1.0D0 - fs_var_1675) * fac10_var_1619(jl_var_1694, lay_var_1649)
          fac100_var_1655 = fs_var_1675 * fac00_var_1617(jl_var_1694, lay_var_1649)
          fac110_var_1658 = fs_var_1675 * fac10_var_1619(jl_var_1694, lay_var_1649)
          fac200_var_1656 = 0.0D0
          fac210_var_1659 = 0.0D0
        END IF
        IF (specparm1_var_1680 .LT. 0.125D0) THEN
          p_var_1666 = fs1_var_1678 - 1.0D0
          p4_var_1667 = p_var_1666 ** 4
          fk0_var_1668 = p4_var_1667
          fk1_var_1669 = 1.0D0 - p_var_1666 - 2.0D0 * p4_var_1667
          fk2_var_1670 = p_var_1666 + p4_var_1667
          fac001_var_1660 = fk0_var_1668 * fac01_var_1618(jl_var_1694, lay_var_1649)
          fac101_var_1661 = fk1_var_1669 * fac01_var_1618(jl_var_1694, lay_var_1649)
          fac201_var_1662 = fk2_var_1670 * fac01_var_1618(jl_var_1694, lay_var_1649)
          fac011_var_1663 = fk0_var_1668 * fac11_var_1620(jl_var_1694, lay_var_1649)
          fac111_var_1664 = fk1_var_1669 * fac11_var_1620(jl_var_1694, lay_var_1649)
          fac211_var_1665 = fk2_var_1670 * fac11_var_1620(jl_var_1694, lay_var_1649)
        ELSE IF (specparm1_var_1680 .GT. 0.875D0) THEN
          p_var_1666 = - fs1_var_1678
          p4_var_1667 = p_var_1666 ** 4
          fk0_var_1668 = p4_var_1667
          fk1_var_1669 = 1.0D0 - p_var_1666 - 2.0D0 * p4_var_1667
          fk2_var_1670 = p_var_1666 + p4_var_1667
          fac001_var_1660 = fk0_var_1668 * fac01_var_1618(jl_var_1694, lay_var_1649)
          fac101_var_1661 = fk1_var_1669 * fac01_var_1618(jl_var_1694, lay_var_1649)
          fac201_var_1662 = fk2_var_1670 * fac01_var_1618(jl_var_1694, lay_var_1649)
          fac011_var_1663 = fk0_var_1668 * fac11_var_1620(jl_var_1694, lay_var_1649)
          fac111_var_1664 = fk1_var_1669 * fac11_var_1620(jl_var_1694, lay_var_1649)
          fac211_var_1665 = fk2_var_1670 * fac11_var_1620(jl_var_1694, lay_var_1649)
        ELSE
          fac001_var_1660 = (1.0D0 - fs1_var_1678) * fac01_var_1618(jl_var_1694, lay_var_1649)
          fac011_var_1663 = (1.0D0 - fs1_var_1678) * fac11_var_1620(jl_var_1694, lay_var_1649)
          fac101_var_1661 = fs1_var_1678 * fac01_var_1618(jl_var_1694, lay_var_1649)
          fac111_var_1664 = fs1_var_1678 * fac11_var_1620(jl_var_1694, lay_var_1649)
          fac201_var_1662 = 0.0D0
          fac211_var_1665 = 0.0D0
        END IF
        IF (specparm_var_1677 .LT. 0.125D0) THEN
          tau_major_var_1673(1 : ng4) = speccomb_var_1640 * (fac000_var_1654 * absa_var_169(ind0_var_1643, 1 : 14) + fac100_var_1655 * absa_var_169(ind0_var_1643 + 1, 1 : 14) + fac200_var_1656 * absa_var_169(ind0_var_1643 + 2, 1 : 14) + fac010_var_1657 * absa_var_169(ind0_var_1643 + 9, 1 : 14) + fac110_var_1658 * absa_var_169(ind0_var_1643 + 10, 1 : 14) + fac210_var_1659 * absa_var_169(ind0_var_1643 + 11, 1 : 14))
        ELSE IF (specparm_var_1677 .GT. 0.875D0) THEN
          tau_major_var_1673(1 : ng4) = speccomb_var_1640 * (fac200_var_1656 * absa_var_169(ind0_var_1643 - 1, 1 : 14) + fac100_var_1655 * absa_var_169(ind0_var_1643, 1 : 14) + fac000_var_1654 * absa_var_169(ind0_var_1643 + 1, 1 : 14) + fac210_var_1659 * absa_var_169(ind0_var_1643 + 8, 1 : 14) + fac110_var_1658 * absa_var_169(ind0_var_1643 + 9, 1 : 14) + fac010_var_1657 * absa_var_169(ind0_var_1643 + 10, 1 : 14))
        ELSE
          tau_major_var_1673(1 : ng4) = speccomb_var_1640 * (fac000_var_1654 * absa_var_169(ind0_var_1643, 1 : 14) + fac100_var_1655 * absa_var_169(ind0_var_1643 + 1, 1 : 14) + fac010_var_1657 * absa_var_169(ind0_var_1643 + 9, 1 : 14) + fac110_var_1658 * absa_var_169(ind0_var_1643 + 10, 1 : 14))
        END IF
        IF (specparm1_var_1680 .LT. 0.125D0) THEN
          tau_major1_var_1674(1 : ng4) = speccomb1_var_1641 * (fac001_var_1660 * absa_var_169(ind1_var_1644, 1 : 14) + fac101_var_1661 * absa_var_169(ind1_var_1644 + 1, 1 : 14) + fac201_var_1662 * absa_var_169(ind1_var_1644 + 2, 1 : 14) + fac011_var_1663 * absa_var_169(ind1_var_1644 + 9, 1 : 14) + fac111_var_1664 * absa_var_169(ind1_var_1644 + 10, 1 : 14) + fac211_var_1665 * absa_var_169(ind1_var_1644 + 11, 1 : 14))
        ELSE IF (specparm1_var_1680 .GT. 0.875D0) THEN
          tau_major1_var_1674(1 : ng4) = speccomb1_var_1641 * (fac201_var_1662 * absa_var_169(ind1_var_1644 - 1, 1 : 14) + fac101_var_1661 * absa_var_169(ind1_var_1644, 1 : 14) + fac001_var_1660 * absa_var_169(ind1_var_1644 + 1, 1 : 14) + fac211_var_1665 * absa_var_169(ind1_var_1644 + 8, 1 : 14) + fac111_var_1664 * absa_var_169(ind1_var_1644 + 9, 1 : 14) + fac011_var_1663 * absa_var_169(ind1_var_1644 + 10, 1 : 14))
        ELSE
          tau_major1_var_1674(1 : ng4) = speccomb1_var_1641 * (fac001_var_1660 * absa_var_169(ind1_var_1644, 1 : 14) + fac101_var_1661 * absa_var_169(ind1_var_1644 + 1, 1 : 14) + fac011_var_1663 * absa_var_169(ind1_var_1644 + 9, 1 : 14) + fac111_var_1664 * absa_var_169(ind1_var_1644 + 10, 1 : 14))
        END IF
        DO ig_var_1647 = 1, 14
          tauself_var_1672 = selffac_var_1629(jl_var_1694, lay_var_1649) * (selfref_var_171(inds_var_1645, ig_var_1647) + selffrac_var_1630(jl_var_1694, lay_var_1649) * (selfref_var_171(inds_var_1645 + 1, ig_var_1647) - selfref_var_171(inds_var_1645, ig_var_1647)))
          taufor_var_1671 = forfac_var_1638(jl_var_1694, lay_var_1649) * (forref_var_172(indf_var_1646, ig_var_1647) + forfrac_var_1639(jl_var_1694, lay_var_1649) * (forref_var_172(indf_var_1646 + 1, ig_var_1647) - forref_var_172(indf_var_1646, ig_var_1647)))
          taug_var_1615(jl_var_1694, 38 + ig_var_1647, lay_var_1649) = tau_major_var_1673(ig_var_1647) + tau_major1_var_1674(ig_var_1647) + tauself_var_1672 + taufor_var_1671
          fracs_var_1632(jl_var_1694, 38 + ig_var_1647, lay_var_1649) = fracrefa_var_167(ig_var_1647, jpl_var_1651) + fpl_var_1681 * (fracrefa_var_167(ig_var_1647, jpl_var_1651 + 1) - fracrefa_var_167(ig_var_1647, jpl_var_1651))
        END DO
      END DO
      ixc0_var_1691 = kfdia_var_1613 - kidia_var_1612 + 1 - ixc0_var_1691
      DO ixp_var_1692 = 1, ixc0_var_1691
        jl_var_1694 = ixhigh_var_1688(ixp_var_1692, lay_var_1649)
        speccomb_var_1640 = colo3_var_1627(jl_var_1694, lay_var_1649) + rat_o3co2_var_1635(jl_var_1694, lay_var_1649) * colco2_var_1626(jl_var_1694, lay_var_1649)
        specparm_var_1677 = MIN(colo3_var_1627(jl_var_1694, lay_var_1649) / speccomb_var_1640, oneminus_var_1624)
        specmult_var_1676 = 4.0D0 * (specparm_var_1677)
        js_var_1648 = 1 + INT(specmult_var_1676)
        fs_var_1675 = ((specmult_var_1676) - AINT((specmult_var_1676)))
        speccomb1_var_1641 = colo3_var_1627(jl_var_1694, lay_var_1649) + rat_o3co2_1_var_1636(jl_var_1694, lay_var_1649) * colco2_var_1626(jl_var_1694, lay_var_1649)
        specparm1_var_1680 = MIN(colo3_var_1627(jl_var_1694, lay_var_1649) / speccomb1_var_1641, oneminus_var_1624)
        specmult1_var_1679 = 4.0D0 * (specparm1_var_1680)
        js1_var_1650 = 1 + INT(specmult1_var_1679)
        fs1_var_1678 = ((specmult1_var_1679) - AINT((specmult1_var_1679)))
        fac000_var_1654 = (1.0D0 - fs_var_1675) * fac00_var_1617(jl_var_1694, lay_var_1649)
        fac010_var_1657 = (1.0D0 - fs_var_1675) * fac10_var_1619(jl_var_1694, lay_var_1649)
        fac100_var_1655 = fs_var_1675 * fac00_var_1617(jl_var_1694, lay_var_1649)
        fac110_var_1658 = fs_var_1675 * fac10_var_1619(jl_var_1694, lay_var_1649)
        fac001_var_1660 = (1.0D0 - fs1_var_1678) * fac01_var_1618(jl_var_1694, lay_var_1649)
        fac011_var_1663 = (1.0D0 - fs1_var_1678) * fac11_var_1620(jl_var_1694, lay_var_1649)
        fac101_var_1661 = fs1_var_1678 * fac01_var_1618(jl_var_1694, lay_var_1649)
        fac111_var_1664 = fs1_var_1678 * fac11_var_1620(jl_var_1694, lay_var_1649)
        speccomb_planck_var_1642 = colo3_var_1627(jl_var_1694, lay_var_1649) + refrat_planck_b_var_1653 * colco2_var_1626(jl_var_1694, lay_var_1649)
        specparm_planck_var_1683 = MIN(colo3_var_1627(jl_var_1694, lay_var_1649) / speccomb_planck_var_1642, oneminus_var_1624)
        specmult_planck_var_1682 = 4.0D0 * specparm_planck_var_1683
        jpl_var_1651 = 1 + INT(specmult_planck_var_1682)
        fpl_var_1681 = ((specmult_planck_var_1682) - AINT((specmult_planck_var_1682)))
        ind0_var_1643 = ((jp_var_1621(jl_var_1694, lay_var_1649) - 13) * 5 + (jt_var_1622(jl_var_1694, lay_var_1649) - 1)) * nspb_var_217(4) + js_var_1648
        ind1_var_1644 = ((jp_var_1621(jl_var_1694, lay_var_1649) - 12) * 5 + (jt1_var_1623(jl_var_1694, lay_var_1649) - 1)) * nspb_var_217(4) + js1_var_1650
        DO ig_var_1647 = 1, 14
          taug_var_1615(jl_var_1694, 38 + ig_var_1647, lay_var_1649) = speccomb_var_1640 * (fac000_var_1654 * absb_var_170(ind0_var_1643, ig_var_1647) + fac100_var_1655 * absb_var_170(ind0_var_1643 + 1, ig_var_1647) + fac010_var_1657 * absb_var_170(ind0_var_1643 + 5, ig_var_1647) + fac110_var_1658 * absb_var_170(ind0_var_1643 + 6, ig_var_1647)) + speccomb1_var_1641 * (fac001_var_1660 * absb_var_170(ind1_var_1644, ig_var_1647) + fac101_var_1661 * absb_var_170(ind1_var_1644 + 1, ig_var_1647) + fac011_var_1663 * absb_var_170(ind1_var_1644 + 5, ig_var_1647) + fac111_var_1664 * absb_var_170(ind1_var_1644 + 6, ig_var_1647))
          fracs_var_1632(jl_var_1694, 38 + ig_var_1647, lay_var_1649) = fracrefb_var_168(ig_var_1647, jpl_var_1651) + fpl_var_1681 * (fracrefb_var_168(ig_var_1647, jpl_var_1651 + 1) - fracrefb_var_168(ig_var_1647, jpl_var_1651))
        END DO
      END DO
      DO ixp_var_1692 = 1, ixc0_var_1691
        jl_var_1694 = ixhigh_var_1688(ixp_var_1692, lay_var_1649)
        taug_var_1615(jl_var_1694, 46, lay_var_1649) = taug_var_1615(jl_var_1694, 46, lay_var_1649) * 0.92D0
        taug_var_1615(jl_var_1694, 47, lay_var_1649) = taug_var_1615(jl_var_1694, 47, lay_var_1649) * 0.88D0
        taug_var_1615(jl_var_1694, 48, lay_var_1649) = taug_var_1615(jl_var_1694, 48, lay_var_1649) * 1.07D0
        taug_var_1615(jl_var_1694, 49, lay_var_1649) = taug_var_1615(jl_var_1694, 49, lay_var_1649) * 1.1D0
        taug_var_1615(jl_var_1694, 50, lay_var_1649) = taug_var_1615(jl_var_1694, 50, lay_var_1649) * 0.99D0
        taug_var_1615(jl_var_1694, 51, lay_var_1649) = taug_var_1615(jl_var_1694, 51, lay_var_1649) * 0.88D0
        taug_var_1615(jl_var_1694, 52, lay_var_1649) = taug_var_1615(jl_var_1694, 52, lay_var_1649) * 0.943D0
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol4
SUBROUTINE rrtm_taumol6(kidia_var_1695, kfdia_var_1696, klev_var_1697, taug_var_1698, wx_var_1699, p_tauaerl_var_1700, fac00_var_1701, fac01_var_1702, fac10_var_1703, fac11_var_1704, forfac_var_1717, forfrac_var_1718, indfor_var_1716, jp_var_1705, jt_var_1706, jt1_var_1707, colh2o_var_1708, colco2_var_1709, coldry_var_1710, laytrop_var_1711, selffac_var_1712, selffrac_var_1713, indself_var_1714, fracs_var_1715, minorfrac_var_1719, indminor_var_1720)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrta6, ONLY: absa_var_182, cfc11adj, cfc12_var_181, forref_var_185, fracrefa_var_180, ka_mco2_var_184, selfref_var_183
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1695
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1696
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1697
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1698(kidia_var_1695 : kfdia_var_1696, 140, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: wx_var_1699(kidia_var_1695 : kfdia_var_1696, 4, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1700(kidia_var_1695 : kfdia_var_1696, klev_var_1697, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1701(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1702(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1703(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1704(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1705(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1706(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1707(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1708(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1709(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1710(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1711(kidia_var_1695 : kfdia_var_1696)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1712(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1713(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1714(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1715(kidia_var_1695 : kfdia_var_1696, 140, klev_var_1697)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1716(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1717(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1718(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1719(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1720(kidia_var_1695 : kfdia_var_1696, klev_var_1697)
  INTEGER(KIND = 4) :: ind0_var_1721, ind1_var_1722, inds_var_1723, indf_var_1724, indm_var_1725
  INTEGER(KIND = 4) :: ig_var_1726, lay_var_1727
  REAL(KIND = 8) :: adjfac_var_1728, adjcolco2_var_1729, ratco2_var_1730, chi_co2_var_1731
  REAL(KIND = 8) :: taufor_var_1732, tauself_var_1733, absco2_var_1734
  INTEGER(KIND = 4) :: laytrop_min_var_1735, laytrop_max_var_1736
  INTEGER(KIND = 4) :: ixc_var_1737(klev_var_1697), ixlow_var_1738(kfdia_var_1696, klev_var_1697), ixhigh_var_1739(kfdia_var_1696, klev_var_1697)
  INTEGER(KIND = 4) :: ich_var_1740, icl_var_1741, ixc0_var_1742, ixp_var_1743, jc_var_1744, jl_var_1745
  laytrop_min_var_1735 = MINVAL(laytrop_var_1711)
  laytrop_max_var_1736 = MAXVAL(laytrop_var_1711)
  ixlow_var_1738 = 0
  ixhigh_var_1739 = 0
  ixc_var_1737 = 0
  DO lay_var_1727 = laytrop_min_var_1735 + 1, laytrop_max_var_1736
    icl_var_1741 = 0
    ich_var_1740 = 0
    DO jc_var_1744 = kidia_var_1695, kfdia_var_1696
      IF (lay_var_1727 <= laytrop_var_1711(jc_var_1744)) THEN
        icl_var_1741 = icl_var_1741 + 1
        ixlow_var_1738(icl_var_1741, lay_var_1727) = jc_var_1744
      ELSE
        ich_var_1740 = ich_var_1740 + 1
        ixhigh_var_1739(ich_var_1740, lay_var_1727) = jc_var_1744
      END IF
    END DO
    ixc_var_1737(lay_var_1727) = icl_var_1741
  END DO
  DO lay_var_1727 = 1, laytrop_min_var_1735
    DO jl_var_1745 = kidia_var_1695, kfdia_var_1696
      chi_co2_var_1731 = colco2_var_1709(jl_var_1745, lay_var_1727) / (coldry_var_1710(jl_var_1745, lay_var_1727))
      ratco2_var_1730 = 1D+20 * chi_co2_var_1731 / chi_mls(2, jp_var_1705(jl_var_1745, lay_var_1727) + 1)
      IF (ratco2_var_1730 .GT. 3.0D0) THEN
        adjfac_var_1728 = 2.0D0 + (ratco2_var_1730 - 2.0D0) ** 0.77D0
        adjcolco2_var_1729 = adjfac_var_1728 * chi_mls(2, jp_var_1705(jl_var_1745, lay_var_1727) + 1) * coldry_var_1710(jl_var_1745, lay_var_1727) * 1D-20
      ELSE
        adjcolco2_var_1729 = colco2_var_1709(jl_var_1745, lay_var_1727)
      END IF
      ind0_var_1721 = ((jp_var_1705(jl_var_1745, lay_var_1727) - 1) * 5 + (jt_var_1706(jl_var_1745, lay_var_1727) - 1)) * nspa_var_216(6) + 1
      ind1_var_1722 = (jp_var_1705(jl_var_1745, lay_var_1727) * 5 + (jt1_var_1707(jl_var_1745, lay_var_1727) - 1)) * nspa_var_216(6) + 1
      inds_var_1723 = indself_var_1714(jl_var_1745, lay_var_1727)
      indf_var_1724 = indfor_var_1716(jl_var_1745, lay_var_1727)
      indm_var_1725 = indminor_var_1720(jl_var_1745, lay_var_1727)
      DO ig_var_1726 = 1, 8
        tauself_var_1733 = selffac_var_1712(jl_var_1745, lay_var_1727) * (selfref_var_183(inds_var_1723, ig_var_1726) + selffrac_var_1713(jl_var_1745, lay_var_1727) * (selfref_var_183(inds_var_1723 + 1, ig_var_1726) - selfref_var_183(inds_var_1723, ig_var_1726)))
        taufor_var_1732 = forfac_var_1717(jl_var_1745, lay_var_1727) * (forref_var_185(indf_var_1724, ig_var_1726) + forfrac_var_1718(jl_var_1745, lay_var_1727) * (forref_var_185(indf_var_1724 + 1, ig_var_1726) - forref_var_185(indf_var_1724, ig_var_1726)))
        absco2_var_1734 = (ka_mco2_var_184(indm_var_1725, ig_var_1726) + minorfrac_var_1719(jl_var_1745, lay_var_1727) * (ka_mco2_var_184(indm_var_1725 + 1, ig_var_1726) - ka_mco2_var_184(indm_var_1725, ig_var_1726)))
        taug_var_1698(jl_var_1745, 68 + ig_var_1726, lay_var_1727) = colh2o_var_1708(jl_var_1745, lay_var_1727) * (fac00_var_1701(jl_var_1745, lay_var_1727) * absa_var_182(ind0_var_1721, ig_var_1726) + fac10_var_1703(jl_var_1745, lay_var_1727) * absa_var_182(ind0_var_1721 + 1, ig_var_1726) + fac01_var_1702(jl_var_1745, lay_var_1727) * absa_var_182(ind1_var_1722, ig_var_1726) + fac11_var_1704(jl_var_1745, lay_var_1727) * absa_var_182(ind1_var_1722 + 1, ig_var_1726)) + tauself_var_1733 + taufor_var_1732 + adjcolco2_var_1729 * absco2_var_1734 + wx_var_1699(jl_var_1745, 2, lay_var_1727) * cfc11adj(ig_var_1726) + wx_var_1699(jl_var_1745, 3, lay_var_1727) * cfc12_var_181(ig_var_1726)
        fracs_var_1715(jl_var_1745, 68 + ig_var_1726, lay_var_1727) = fracrefa_var_180(ig_var_1726)
      END DO
    END DO
  END DO
  DO ig_var_1726 = 1, 8
    DO lay_var_1727 = laytrop_max_var_1736 + 1, klev_var_1697
      DO jl_var_1745 = kidia_var_1695, kfdia_var_1696
        taug_var_1698(jl_var_1745, 68 + ig_var_1726, lay_var_1727) = 0.0D0 + wx_var_1699(jl_var_1745, 2, lay_var_1727) * cfc11adj(ig_var_1726) + wx_var_1699(jl_var_1745, 3, lay_var_1727) * cfc12_var_181(ig_var_1726)
        fracs_var_1715(jl_var_1745, 68 + ig_var_1726, lay_var_1727) = fracrefa_var_180(ig_var_1726)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1736 /= laytrop_min_var_1735) THEN
    DO lay_var_1727 = laytrop_min_var_1735 + 1, laytrop_max_var_1736
      ixc0_var_1742 = ixc_var_1737(lay_var_1727)
      DO ixp_var_1743 = 1, ixc0_var_1742
        jl_var_1745 = ixlow_var_1738(ixp_var_1743, lay_var_1727)
        chi_co2_var_1731 = colco2_var_1709(jl_var_1745, lay_var_1727) / (coldry_var_1710(jl_var_1745, lay_var_1727))
        ratco2_var_1730 = 1D+20 * chi_co2_var_1731 / chi_mls(2, jp_var_1705(jl_var_1745, lay_var_1727) + 1)
        IF (ratco2_var_1730 .GT. 3.0D0) THEN
          adjfac_var_1728 = 2.0D0 + (ratco2_var_1730 - 2.0D0) ** 0.77D0
          adjcolco2_var_1729 = adjfac_var_1728 * chi_mls(2, jp_var_1705(jl_var_1745, lay_var_1727) + 1) * coldry_var_1710(jl_var_1745, lay_var_1727) * 1D-20
        ELSE
          adjcolco2_var_1729 = colco2_var_1709(jl_var_1745, lay_var_1727)
        END IF
        ind0_var_1721 = ((jp_var_1705(jl_var_1745, lay_var_1727) - 1) * 5 + (jt_var_1706(jl_var_1745, lay_var_1727) - 1)) * nspa_var_216(6) + 1
        ind1_var_1722 = (jp_var_1705(jl_var_1745, lay_var_1727) * 5 + (jt1_var_1707(jl_var_1745, lay_var_1727) - 1)) * nspa_var_216(6) + 1
        inds_var_1723 = indself_var_1714(jl_var_1745, lay_var_1727)
        indf_var_1724 = indfor_var_1716(jl_var_1745, lay_var_1727)
        indm_var_1725 = indminor_var_1720(jl_var_1745, lay_var_1727)
        DO ig_var_1726 = 1, 8
          tauself_var_1733 = selffac_var_1712(jl_var_1745, lay_var_1727) * (selfref_var_183(inds_var_1723, ig_var_1726) + selffrac_var_1713(jl_var_1745, lay_var_1727) * (selfref_var_183(inds_var_1723 + 1, ig_var_1726) - selfref_var_183(inds_var_1723, ig_var_1726)))
          taufor_var_1732 = forfac_var_1717(jl_var_1745, lay_var_1727) * (forref_var_185(indf_var_1724, ig_var_1726) + forfrac_var_1718(jl_var_1745, lay_var_1727) * (forref_var_185(indf_var_1724 + 1, ig_var_1726) - forref_var_185(indf_var_1724, ig_var_1726)))
          absco2_var_1734 = (ka_mco2_var_184(indm_var_1725, ig_var_1726) + minorfrac_var_1719(jl_var_1745, lay_var_1727) * (ka_mco2_var_184(indm_var_1725 + 1, ig_var_1726) - ka_mco2_var_184(indm_var_1725, ig_var_1726)))
          taug_var_1698(jl_var_1745, 68 + ig_var_1726, lay_var_1727) = colh2o_var_1708(jl_var_1745, lay_var_1727) * (fac00_var_1701(jl_var_1745, lay_var_1727) * absa_var_182(ind0_var_1721, ig_var_1726) + fac10_var_1703(jl_var_1745, lay_var_1727) * absa_var_182(ind0_var_1721 + 1, ig_var_1726) + fac01_var_1702(jl_var_1745, lay_var_1727) * absa_var_182(ind1_var_1722, ig_var_1726) + fac11_var_1704(jl_var_1745, lay_var_1727) * absa_var_182(ind1_var_1722 + 1, ig_var_1726)) + tauself_var_1733 + taufor_var_1732 + adjcolco2_var_1729 * absco2_var_1734 + wx_var_1699(jl_var_1745, 2, lay_var_1727) * cfc11adj(ig_var_1726) + wx_var_1699(jl_var_1745, 3, lay_var_1727) * cfc12_var_181(ig_var_1726)
          fracs_var_1715(jl_var_1745, 68 + ig_var_1726, lay_var_1727) = fracrefa_var_180(ig_var_1726)
        END DO
      END DO
      ixc0_var_1742 = kfdia_var_1696 - kidia_var_1695 + 1 - ixc0_var_1742
      DO ig_var_1726 = 1, 8
        DO ixp_var_1743 = 1, ixc0_var_1742
          jl_var_1745 = ixhigh_var_1739(ixp_var_1743, lay_var_1727)
          taug_var_1698(jl_var_1745, 68 + ig_var_1726, lay_var_1727) = 0.0D0 + wx_var_1699(jl_var_1745, 2, lay_var_1727) * cfc11adj(ig_var_1726) + wx_var_1699(jl_var_1745, 3, lay_var_1727) * cfc12_var_181(ig_var_1726)
          fracs_var_1715(jl_var_1745, 68 + ig_var_1726, lay_var_1727) = fracrefa_var_180(ig_var_1726)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol6
SUBROUTINE rrtm_taumol7(kidia_var_1746, kfdia_var_1747, klev_var_1748, taug_var_1749, p_tauaerl_var_1750, fac00_var_1751, fac01_var_1752, fac10_var_1753, fac11_var_1754, forfac_var_1770, forfrac_var_1769, indfor_var_1768, jp_var_1755, jt_var_1756, jt1_var_1757, oneminus_var_1758, colh2o_var_1759, colo3_var_1760, colco2_var_1761, coldry_var_1762, laytrop_var_1763, selffac_var_1764, selffrac_var_1765, indself_var_1766, fracs_var_1767, rat_h2oo3, rat_h2oo3_1, minorfrac_var_1771, indminor_var_1772)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng7
  USE yoerrta7, ONLY: absa_var_188, absb_var_189, forref_var_193, fracrefa_var_186, fracrefb_var_187, ka_mco2_var_191, kb_mco2_var_192, selfref_var_190
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1746
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1747
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1748
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1749(kidia_var_1746 : kfdia_var_1747, 140, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1750(kidia_var_1746 : kfdia_var_1747, klev_var_1748, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1751(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1752(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1753(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1754(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1755(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1756(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1757(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1758
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1759(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1760(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1761(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1762(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1763(kidia_var_1746 : kfdia_var_1747)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1764(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1765(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1766(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1767(kidia_var_1746 : kfdia_var_1747, 140, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oo3(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oo3_1(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1768(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1769(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1770(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1771(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1772(kidia_var_1746 : kfdia_var_1747, klev_var_1748)
  REAL(KIND = 8) :: speccomb_var_1773, speccomb1_var_1774, speccomb_mco2_var_1775, speccomb_planck_var_1776
  INTEGER(KIND = 4) :: ind0_var_1777, ind1_var_1778, inds_var_1779, indf_var_1780, indm_var_1781
  INTEGER(KIND = 4) :: ig_var_1782, js_var_1783, lay_var_1784, js1_var_1785, jpl_var_1786, jmco2_var_1787
  REAL(KIND = 8) :: refrat_planck_a_var_1788, refrat_m_a_var_1789
  REAL(KIND = 8) :: chi_co2_var_1790, ratco2_var_1791, adjfac_var_1792, adjcolco2_var_1793
  REAL(KIND = 8) :: fac000_var_1794, fac100_var_1795, fac200_var_1796, fac010_var_1797, fac110_var_1798, fac210_var_1799, fac001_var_1800, fac101_var_1801, fac201_var_1802, fac011_var_1803, fac111_var_1804, fac211_var_1805
  REAL(KIND = 8) :: p_var_1806, p4_var_1807, fk0_var_1808, fk1_var_1809, fk2_var_1810
  REAL(KIND = 8) :: taufor_var_1811, tauself_var_1812, tau_major_var_1813(12), tau_major1_var_1814(12), co2m1_var_1815, co2m2_var_1816, absco2_var_1817
  REAL(KIND = 8) :: fs_var_1818, specmult_var_1819, specparm_var_1820, fs1_var_1821, specmult1_var_1822, specparm1_var_1823, fpl_var_1824, specmult_planck_var_1825, specparm_planck_var_1826, fmco2_var_1827, specmult_mco2_var_1828, specparm_mco2_var_1829
  INTEGER(KIND = 4) :: laytrop_min_var_1830, laytrop_max_var_1831
  INTEGER(KIND = 4) :: ixc_var_1832(klev_var_1748), ixlow_var_1833(kfdia_var_1747, klev_var_1748), ixhigh_var_1834(kfdia_var_1747, klev_var_1748)
  INTEGER(KIND = 4) :: ich_var_1835, icl_var_1836, ixc0_var_1837, ixp_var_1838, jc_var_1839, jl_var_1840
  laytrop_min_var_1830 = MINVAL(laytrop_var_1763)
  laytrop_max_var_1831 = MAXVAL(laytrop_var_1763)
  ixlow_var_1833 = 0
  ixhigh_var_1834 = 0
  ixc_var_1832 = 0
  DO lay_var_1784 = laytrop_min_var_1830 + 1, laytrop_max_var_1831
    icl_var_1836 = 0
    ich_var_1835 = 0
    DO jc_var_1839 = kidia_var_1746, kfdia_var_1747
      IF (lay_var_1784 <= laytrop_var_1763(jc_var_1839)) THEN
        icl_var_1836 = icl_var_1836 + 1
        ixlow_var_1833(icl_var_1836, lay_var_1784) = jc_var_1839
      ELSE
        ich_var_1835 = ich_var_1835 + 1
        ixhigh_var_1834(ich_var_1835, lay_var_1784) = jc_var_1839
      END IF
    END DO
    ixc_var_1832(lay_var_1784) = icl_var_1836
  END DO
  refrat_planck_a_var_1788 = chi_mls(1, 3) / chi_mls(3, 3)
  refrat_m_a_var_1789 = chi_mls(1, 3) / chi_mls(3, 3)
  DO lay_var_1784 = 1, laytrop_min_var_1830
    DO jl_var_1840 = kidia_var_1746, kfdia_var_1747
      speccomb_var_1773 = colh2o_var_1759(jl_var_1840, lay_var_1784) + rat_h2oo3(jl_var_1840, lay_var_1784) * colo3_var_1760(jl_var_1840, lay_var_1784)
      specparm_var_1820 = MIN(colh2o_var_1759(jl_var_1840, lay_var_1784) / speccomb_var_1773, oneminus_var_1758)
      specmult_var_1819 = 8.0D0 * (specparm_var_1820)
      js_var_1783 = 1 + INT(specmult_var_1819)
      fs_var_1818 = ((specmult_var_1819) - AINT((specmult_var_1819)))
      speccomb1_var_1774 = colh2o_var_1759(jl_var_1840, lay_var_1784) + rat_h2oo3_1(jl_var_1840, lay_var_1784) * colo3_var_1760(jl_var_1840, lay_var_1784)
      specparm1_var_1823 = MIN(colh2o_var_1759(jl_var_1840, lay_var_1784) / speccomb1_var_1774, oneminus_var_1758)
      specmult1_var_1822 = 8.0D0 * (specparm1_var_1823)
      js1_var_1785 = 1 + INT(specmult1_var_1822)
      fs1_var_1821 = ((specmult1_var_1822) - AINT((specmult1_var_1822)))
      speccomb_mco2_var_1775 = colh2o_var_1759(jl_var_1840, lay_var_1784) + refrat_m_a_var_1789 * colo3_var_1760(jl_var_1840, lay_var_1784)
      specparm_mco2_var_1829 = MIN(colh2o_var_1759(jl_var_1840, lay_var_1784) / speccomb_mco2_var_1775, oneminus_var_1758)
      specmult_mco2_var_1828 = 8.0D0 * specparm_mco2_var_1829
      jmco2_var_1787 = 1 + INT(specmult_mco2_var_1828)
      fmco2_var_1827 = ((specmult_mco2_var_1828) - AINT((specmult_mco2_var_1828)))
      chi_co2_var_1790 = colco2_var_1761(jl_var_1840, lay_var_1784) / (coldry_var_1762(jl_var_1840, lay_var_1784))
      ratco2_var_1791 = 1D+20 * chi_co2_var_1790 / chi_mls(2, jp_var_1755(jl_var_1840, lay_var_1784) + 1)
      IF (ratco2_var_1791 .GT. 3.0D0) THEN
        adjfac_var_1792 = 3.0D0 + (ratco2_var_1791 - 3.0D0) ** 0.79D0
        adjcolco2_var_1793 = adjfac_var_1792 * chi_mls(2, jp_var_1755(jl_var_1840, lay_var_1784) + 1) * coldry_var_1762(jl_var_1840, lay_var_1784) * 1D-20
      ELSE
        adjcolco2_var_1793 = colco2_var_1761(jl_var_1840, lay_var_1784)
      END IF
      speccomb_planck_var_1776 = colh2o_var_1759(jl_var_1840, lay_var_1784) + refrat_planck_a_var_1788 * colo3_var_1760(jl_var_1840, lay_var_1784)
      specparm_planck_var_1826 = MIN(colh2o_var_1759(jl_var_1840, lay_var_1784) / speccomb_planck_var_1776, oneminus_var_1758)
      specmult_planck_var_1825 = 8.0D0 * specparm_planck_var_1826
      jpl_var_1786 = 1 + INT(specmult_planck_var_1825)
      fpl_var_1824 = ((specmult_planck_var_1825) - AINT((specmult_planck_var_1825)))
      ind0_var_1777 = ((jp_var_1755(jl_var_1840, lay_var_1784) - 1) * 5 + (jt_var_1756(jl_var_1840, lay_var_1784) - 1)) * nspa_var_216(7) + js_var_1783
      ind1_var_1778 = (jp_var_1755(jl_var_1840, lay_var_1784) * 5 + (jt1_var_1757(jl_var_1840, lay_var_1784) - 1)) * nspa_var_216(7) + js1_var_1785
      inds_var_1779 = indself_var_1766(jl_var_1840, lay_var_1784)
      indf_var_1780 = indfor_var_1768(jl_var_1840, lay_var_1784)
      indm_var_1781 = indminor_var_1772(jl_var_1840, lay_var_1784)
      IF (specparm_var_1820 .LT. 0.125D0) THEN
        p_var_1806 = fs_var_1818 - 1.0D0
        p4_var_1807 = p_var_1806 ** 4
        fk0_var_1808 = p4_var_1807
        fk1_var_1809 = 1.0D0 - p_var_1806 - 2.0D0 * p4_var_1807
        fk2_var_1810 = p_var_1806 + p4_var_1807
        fac000_var_1794 = fk0_var_1808 * fac00_var_1751(jl_var_1840, lay_var_1784)
        fac100_var_1795 = fk1_var_1809 * fac00_var_1751(jl_var_1840, lay_var_1784)
        fac200_var_1796 = fk2_var_1810 * fac00_var_1751(jl_var_1840, lay_var_1784)
        fac010_var_1797 = fk0_var_1808 * fac10_var_1753(jl_var_1840, lay_var_1784)
        fac110_var_1798 = fk1_var_1809 * fac10_var_1753(jl_var_1840, lay_var_1784)
        fac210_var_1799 = fk2_var_1810 * fac10_var_1753(jl_var_1840, lay_var_1784)
      ELSE IF (specparm_var_1820 .GT. 0.875D0) THEN
        p_var_1806 = - fs_var_1818
        p4_var_1807 = p_var_1806 ** 4
        fk0_var_1808 = p4_var_1807
        fk1_var_1809 = 1.0D0 - p_var_1806 - 2.0D0 * p4_var_1807
        fk2_var_1810 = p_var_1806 + p4_var_1807
        fac000_var_1794 = fk0_var_1808 * fac00_var_1751(jl_var_1840, lay_var_1784)
        fac100_var_1795 = fk1_var_1809 * fac00_var_1751(jl_var_1840, lay_var_1784)
        fac200_var_1796 = fk2_var_1810 * fac00_var_1751(jl_var_1840, lay_var_1784)
        fac010_var_1797 = fk0_var_1808 * fac10_var_1753(jl_var_1840, lay_var_1784)
        fac110_var_1798 = fk1_var_1809 * fac10_var_1753(jl_var_1840, lay_var_1784)
        fac210_var_1799 = fk2_var_1810 * fac10_var_1753(jl_var_1840, lay_var_1784)
      ELSE
        fac000_var_1794 = (1.0D0 - fs_var_1818) * fac00_var_1751(jl_var_1840, lay_var_1784)
        fac010_var_1797 = (1.0D0 - fs_var_1818) * fac10_var_1753(jl_var_1840, lay_var_1784)
        fac100_var_1795 = fs_var_1818 * fac00_var_1751(jl_var_1840, lay_var_1784)
        fac110_var_1798 = fs_var_1818 * fac10_var_1753(jl_var_1840, lay_var_1784)
        fac200_var_1796 = 0.0D0
        fac210_var_1799 = 0.0D0
      END IF
      IF (specparm1_var_1823 .LT. 0.125D0) THEN
        p_var_1806 = fs1_var_1821 - 1.0D0
        p4_var_1807 = p_var_1806 ** 4
        fk0_var_1808 = p4_var_1807
        fk1_var_1809 = 1.0D0 - p_var_1806 - 2.0D0 * p4_var_1807
        fk2_var_1810 = p_var_1806 + p4_var_1807
        fac001_var_1800 = fk0_var_1808 * fac01_var_1752(jl_var_1840, lay_var_1784)
        fac101_var_1801 = fk1_var_1809 * fac01_var_1752(jl_var_1840, lay_var_1784)
        fac201_var_1802 = fk2_var_1810 * fac01_var_1752(jl_var_1840, lay_var_1784)
        fac011_var_1803 = fk0_var_1808 * fac11_var_1754(jl_var_1840, lay_var_1784)
        fac111_var_1804 = fk1_var_1809 * fac11_var_1754(jl_var_1840, lay_var_1784)
        fac211_var_1805 = fk2_var_1810 * fac11_var_1754(jl_var_1840, lay_var_1784)
      ELSE IF (specparm1_var_1823 .GT. 0.875D0) THEN
        p_var_1806 = - fs1_var_1821
        p4_var_1807 = p_var_1806 ** 4
        fk0_var_1808 = p4_var_1807
        fk1_var_1809 = 1.0D0 - p_var_1806 - 2.0D0 * p4_var_1807
        fk2_var_1810 = p_var_1806 + p4_var_1807
        fac001_var_1800 = fk0_var_1808 * fac01_var_1752(jl_var_1840, lay_var_1784)
        fac101_var_1801 = fk1_var_1809 * fac01_var_1752(jl_var_1840, lay_var_1784)
        fac201_var_1802 = fk2_var_1810 * fac01_var_1752(jl_var_1840, lay_var_1784)
        fac011_var_1803 = fk0_var_1808 * fac11_var_1754(jl_var_1840, lay_var_1784)
        fac111_var_1804 = fk1_var_1809 * fac11_var_1754(jl_var_1840, lay_var_1784)
        fac211_var_1805 = fk2_var_1810 * fac11_var_1754(jl_var_1840, lay_var_1784)
      ELSE
        fac001_var_1800 = (1.0D0 - fs1_var_1821) * fac01_var_1752(jl_var_1840, lay_var_1784)
        fac011_var_1803 = (1.0D0 - fs1_var_1821) * fac11_var_1754(jl_var_1840, lay_var_1784)
        fac101_var_1801 = fs1_var_1821 * fac01_var_1752(jl_var_1840, lay_var_1784)
        fac111_var_1804 = fs1_var_1821 * fac11_var_1754(jl_var_1840, lay_var_1784)
        fac201_var_1802 = 0.0D0
        fac211_var_1805 = 0.0D0
      END IF
      IF (specparm_var_1820 .LT. 0.125D0) THEN
        tau_major_var_1813(1 : ng7) = speccomb_var_1773 * (fac000_var_1794 * absa_var_188(ind0_var_1777, 1 : 12) + fac100_var_1795 * absa_var_188(ind0_var_1777 + 1, 1 : 12) + fac200_var_1796 * absa_var_188(ind0_var_1777 + 2, 1 : 12) + fac010_var_1797 * absa_var_188(ind0_var_1777 + 9, 1 : 12) + fac110_var_1798 * absa_var_188(ind0_var_1777 + 10, 1 : 12) + fac210_var_1799 * absa_var_188(ind0_var_1777 + 11, 1 : 12))
      ELSE IF (specparm_var_1820 .GT. 0.875D0) THEN
        tau_major_var_1813(1 : ng7) = speccomb_var_1773 * (fac200_var_1796 * absa_var_188(ind0_var_1777 - 1, 1 : 12) + fac100_var_1795 * absa_var_188(ind0_var_1777, 1 : 12) + fac000_var_1794 * absa_var_188(ind0_var_1777 + 1, 1 : 12) + fac210_var_1799 * absa_var_188(ind0_var_1777 + 8, 1 : 12) + fac110_var_1798 * absa_var_188(ind0_var_1777 + 9, 1 : 12) + fac010_var_1797 * absa_var_188(ind0_var_1777 + 10, 1 : 12))
      ELSE
        tau_major_var_1813(1 : ng7) = speccomb_var_1773 * (fac000_var_1794 * absa_var_188(ind0_var_1777, 1 : 12) + fac100_var_1795 * absa_var_188(ind0_var_1777 + 1, 1 : 12) + fac010_var_1797 * absa_var_188(ind0_var_1777 + 9, 1 : 12) + fac110_var_1798 * absa_var_188(ind0_var_1777 + 10, 1 : 12))
      END IF
      IF (specparm1_var_1823 .LT. 0.125D0) THEN
        tau_major1_var_1814(1 : ng7) = speccomb1_var_1774 * (fac001_var_1800 * absa_var_188(ind1_var_1778, 1 : 12) + fac101_var_1801 * absa_var_188(ind1_var_1778 + 1, 1 : 12) + fac201_var_1802 * absa_var_188(ind1_var_1778 + 2, 1 : 12) + fac011_var_1803 * absa_var_188(ind1_var_1778 + 9, 1 : 12) + fac111_var_1804 * absa_var_188(ind1_var_1778 + 10, 1 : 12) + fac211_var_1805 * absa_var_188(ind1_var_1778 + 11, 1 : 12))
      ELSE IF (specparm1_var_1823 .GT. 0.875D0) THEN
        tau_major1_var_1814(1 : ng7) = speccomb1_var_1774 * (fac201_var_1802 * absa_var_188(ind1_var_1778 - 1, 1 : 12) + fac101_var_1801 * absa_var_188(ind1_var_1778, 1 : 12) + fac001_var_1800 * absa_var_188(ind1_var_1778 + 1, 1 : 12) + fac211_var_1805 * absa_var_188(ind1_var_1778 + 8, 1 : 12) + fac111_var_1804 * absa_var_188(ind1_var_1778 + 9, 1 : 12) + fac011_var_1803 * absa_var_188(ind1_var_1778 + 10, 1 : 12))
      ELSE
        tau_major1_var_1814(1 : ng7) = speccomb1_var_1774 * (fac001_var_1800 * absa_var_188(ind1_var_1778, 1 : 12) + fac101_var_1801 * absa_var_188(ind1_var_1778 + 1, 1 : 12) + fac011_var_1803 * absa_var_188(ind1_var_1778 + 9, 1 : 12) + fac111_var_1804 * absa_var_188(ind1_var_1778 + 10, 1 : 12))
      END IF
      DO ig_var_1782 = 1, 12
        tauself_var_1812 = selffac_var_1764(jl_var_1840, lay_var_1784) * (selfref_var_190(inds_var_1779, ig_var_1782) + selffrac_var_1765(jl_var_1840, lay_var_1784) * (selfref_var_190(inds_var_1779 + 1, ig_var_1782) - selfref_var_190(inds_var_1779, ig_var_1782)))
        taufor_var_1811 = forfac_var_1770(jl_var_1840, lay_var_1784) * (forref_var_193(indf_var_1780, ig_var_1782) + forfrac_var_1769(jl_var_1840, lay_var_1784) * (forref_var_193(indf_var_1780 + 1, ig_var_1782) - forref_var_193(indf_var_1780, ig_var_1782)))
        co2m1_var_1815 = ka_mco2_var_191(jmco2_var_1787, indm_var_1781, ig_var_1782) + fmco2_var_1827 * (ka_mco2_var_191(jmco2_var_1787 + 1, indm_var_1781, ig_var_1782) - ka_mco2_var_191(jmco2_var_1787, indm_var_1781, ig_var_1782))
        co2m2_var_1816 = ka_mco2_var_191(jmco2_var_1787, indm_var_1781 + 1, ig_var_1782) + fmco2_var_1827 * (ka_mco2_var_191(jmco2_var_1787 + 1, indm_var_1781 + 1, ig_var_1782) - ka_mco2_var_191(jmco2_var_1787, indm_var_1781 + 1, ig_var_1782))
        absco2_var_1817 = co2m1_var_1815 + minorfrac_var_1771(jl_var_1840, lay_var_1784) * (co2m2_var_1816 - co2m1_var_1815)
        taug_var_1749(jl_var_1840, 76 + ig_var_1782, lay_var_1784) = tau_major_var_1813(ig_var_1782) + tau_major1_var_1814(ig_var_1782) + tauself_var_1812 + taufor_var_1811 + adjcolco2_var_1793 * absco2_var_1817
        fracs_var_1767(jl_var_1840, 76 + ig_var_1782, lay_var_1784) = fracrefa_var_186(ig_var_1782, jpl_var_1786) + fpl_var_1824 * (fracrefa_var_186(ig_var_1782, jpl_var_1786 + 1) - fracrefa_var_186(ig_var_1782, jpl_var_1786))
      END DO
    END DO
  END DO
  DO lay_var_1784 = laytrop_max_var_1831 + 1, klev_var_1748
    DO jl_var_1840 = kidia_var_1746, kfdia_var_1747
      chi_co2_var_1790 = colco2_var_1761(jl_var_1840, lay_var_1784) / (coldry_var_1762(jl_var_1840, lay_var_1784))
      ratco2_var_1791 = 1D+20 * chi_co2_var_1790 / chi_mls(2, jp_var_1755(jl_var_1840, lay_var_1784) + 1)
      IF (ratco2_var_1791 .GT. 3.0D0) THEN
        adjfac_var_1792 = 2.0D0 + (ratco2_var_1791 - 2.0D0) ** 0.79D0
        adjcolco2_var_1793 = adjfac_var_1792 * chi_mls(2, jp_var_1755(jl_var_1840, lay_var_1784) + 1) * coldry_var_1762(jl_var_1840, lay_var_1784) * 1D-20
      ELSE
        adjcolco2_var_1793 = colco2_var_1761(jl_var_1840, lay_var_1784)
      END IF
      ind0_var_1777 = ((jp_var_1755(jl_var_1840, lay_var_1784) - 13) * 5 + (jt_var_1756(jl_var_1840, lay_var_1784) - 1)) * nspb_var_217(7) + 1
      ind1_var_1778 = ((jp_var_1755(jl_var_1840, lay_var_1784) - 12) * 5 + (jt1_var_1757(jl_var_1840, lay_var_1784) - 1)) * nspb_var_217(7) + 1
      indm_var_1781 = indminor_var_1772(jl_var_1840, lay_var_1784)
      DO ig_var_1782 = 1, 12
        absco2_var_1817 = kb_mco2_var_192(indm_var_1781, ig_var_1782) + minorfrac_var_1771(jl_var_1840, lay_var_1784) * (kb_mco2_var_192(indm_var_1781 + 1, ig_var_1782) - kb_mco2_var_192(indm_var_1781, ig_var_1782))
        taug_var_1749(jl_var_1840, 76 + ig_var_1782, lay_var_1784) = colo3_var_1760(jl_var_1840, lay_var_1784) * (fac00_var_1751(jl_var_1840, lay_var_1784) * absb_var_189(ind0_var_1777, ig_var_1782) + fac10_var_1753(jl_var_1840, lay_var_1784) * absb_var_189(ind0_var_1777 + 1, ig_var_1782) + fac01_var_1752(jl_var_1840, lay_var_1784) * absb_var_189(ind1_var_1778, ig_var_1782) + fac11_var_1754(jl_var_1840, lay_var_1784) * absb_var_189(ind1_var_1778 + 1, ig_var_1782)) + adjcolco2_var_1793 * absco2_var_1817
        fracs_var_1767(jl_var_1840, 76 + ig_var_1782, lay_var_1784) = fracrefb_var_187(ig_var_1782)
      END DO
    END DO
  END DO
  DO lay_var_1784 = laytrop_max_var_1831 + 1, klev_var_1748
    DO jl_var_1840 = kidia_var_1746, kfdia_var_1747
      taug_var_1749(jl_var_1840, 82, lay_var_1784) = taug_var_1749(jl_var_1840, 82, lay_var_1784) * 0.92D0
      taug_var_1749(jl_var_1840, 83, lay_var_1784) = taug_var_1749(jl_var_1840, 83, lay_var_1784) * 0.88D0
      taug_var_1749(jl_var_1840, 84, lay_var_1784) = taug_var_1749(jl_var_1840, 84, lay_var_1784) * 1.07D0
      taug_var_1749(jl_var_1840, 85, lay_var_1784) = taug_var_1749(jl_var_1840, 85, lay_var_1784) * 1.1D0
      taug_var_1749(jl_var_1840, 86, lay_var_1784) = taug_var_1749(jl_var_1840, 86, lay_var_1784) * 0.99D0
      taug_var_1749(jl_var_1840, 87, lay_var_1784) = taug_var_1749(jl_var_1840, 87, lay_var_1784) * 0.855D0
    END DO
  END DO
  IF (laytrop_max_var_1831 /= laytrop_min_var_1830) THEN
    DO lay_var_1784 = laytrop_min_var_1830 + 1, laytrop_max_var_1831
      ixc0_var_1837 = ixc_var_1832(lay_var_1784)
      DO ixp_var_1838 = 1, ixc0_var_1837
        jl_var_1840 = ixlow_var_1833(ixp_var_1838, lay_var_1784)
        speccomb_var_1773 = colh2o_var_1759(jl_var_1840, lay_var_1784) + rat_h2oo3(jl_var_1840, lay_var_1784) * colo3_var_1760(jl_var_1840, lay_var_1784)
        specparm_var_1820 = MIN(colh2o_var_1759(jl_var_1840, lay_var_1784) / speccomb_var_1773, oneminus_var_1758)
        specmult_var_1819 = 8.0D0 * (specparm_var_1820)
        js_var_1783 = 1 + INT(specmult_var_1819)
        fs_var_1818 = ((specmult_var_1819) - AINT((specmult_var_1819)))
        speccomb1_var_1774 = colh2o_var_1759(jl_var_1840, lay_var_1784) + rat_h2oo3_1(jl_var_1840, lay_var_1784) * colo3_var_1760(jl_var_1840, lay_var_1784)
        specparm1_var_1823 = MIN(colh2o_var_1759(jl_var_1840, lay_var_1784) / speccomb1_var_1774, oneminus_var_1758)
        specmult1_var_1822 = 8.0D0 * (specparm1_var_1823)
        js1_var_1785 = 1 + INT(specmult1_var_1822)
        fs1_var_1821 = ((specmult1_var_1822) - AINT((specmult1_var_1822)))
        speccomb_mco2_var_1775 = colh2o_var_1759(jl_var_1840, lay_var_1784) + refrat_m_a_var_1789 * colo3_var_1760(jl_var_1840, lay_var_1784)
        specparm_mco2_var_1829 = MIN(colh2o_var_1759(jl_var_1840, lay_var_1784) / speccomb_mco2_var_1775, oneminus_var_1758)
        specmult_mco2_var_1828 = 8.0D0 * specparm_mco2_var_1829
        jmco2_var_1787 = 1 + INT(specmult_mco2_var_1828)
        fmco2_var_1827 = ((specmult_mco2_var_1828) - AINT((specmult_mco2_var_1828)))
        chi_co2_var_1790 = colco2_var_1761(jl_var_1840, lay_var_1784) / (coldry_var_1762(jl_var_1840, lay_var_1784))
        ratco2_var_1791 = 1D+20 * chi_co2_var_1790 / chi_mls(2, jp_var_1755(jl_var_1840, lay_var_1784) + 1)
        IF (ratco2_var_1791 .GT. 3.0D0) THEN
          adjfac_var_1792 = 3.0D0 + (ratco2_var_1791 - 3.0D0) ** 0.79D0
          adjcolco2_var_1793 = adjfac_var_1792 * chi_mls(2, jp_var_1755(jl_var_1840, lay_var_1784) + 1) * coldry_var_1762(jl_var_1840, lay_var_1784) * 1D-20
        ELSE
          adjcolco2_var_1793 = colco2_var_1761(jl_var_1840, lay_var_1784)
        END IF
        speccomb_planck_var_1776 = colh2o_var_1759(jl_var_1840, lay_var_1784) + refrat_planck_a_var_1788 * colo3_var_1760(jl_var_1840, lay_var_1784)
        specparm_planck_var_1826 = MIN(colh2o_var_1759(jl_var_1840, lay_var_1784) / speccomb_planck_var_1776, oneminus_var_1758)
        specmult_planck_var_1825 = 8.0D0 * specparm_planck_var_1826
        jpl_var_1786 = 1 + INT(specmult_planck_var_1825)
        fpl_var_1824 = ((specmult_planck_var_1825) - AINT((specmult_planck_var_1825)))
        ind0_var_1777 = ((jp_var_1755(jl_var_1840, lay_var_1784) - 1) * 5 + (jt_var_1756(jl_var_1840, lay_var_1784) - 1)) * nspa_var_216(7) + js_var_1783
        ind1_var_1778 = (jp_var_1755(jl_var_1840, lay_var_1784) * 5 + (jt1_var_1757(jl_var_1840, lay_var_1784) - 1)) * nspa_var_216(7) + js1_var_1785
        inds_var_1779 = indself_var_1766(jl_var_1840, lay_var_1784)
        indf_var_1780 = indfor_var_1768(jl_var_1840, lay_var_1784)
        indm_var_1781 = indminor_var_1772(jl_var_1840, lay_var_1784)
        IF (specparm_var_1820 .LT. 0.125D0) THEN
          p_var_1806 = fs_var_1818 - 1.0D0
          p4_var_1807 = p_var_1806 ** 4
          fk0_var_1808 = p4_var_1807
          fk1_var_1809 = 1.0D0 - p_var_1806 - 2.0D0 * p4_var_1807
          fk2_var_1810 = p_var_1806 + p4_var_1807
          fac000_var_1794 = fk0_var_1808 * fac00_var_1751(jl_var_1840, lay_var_1784)
          fac100_var_1795 = fk1_var_1809 * fac00_var_1751(jl_var_1840, lay_var_1784)
          fac200_var_1796 = fk2_var_1810 * fac00_var_1751(jl_var_1840, lay_var_1784)
          fac010_var_1797 = fk0_var_1808 * fac10_var_1753(jl_var_1840, lay_var_1784)
          fac110_var_1798 = fk1_var_1809 * fac10_var_1753(jl_var_1840, lay_var_1784)
          fac210_var_1799 = fk2_var_1810 * fac10_var_1753(jl_var_1840, lay_var_1784)
        ELSE IF (specparm_var_1820 .GT. 0.875D0) THEN
          p_var_1806 = - fs_var_1818
          p4_var_1807 = p_var_1806 ** 4
          fk0_var_1808 = p4_var_1807
          fk1_var_1809 = 1.0D0 - p_var_1806 - 2.0D0 * p4_var_1807
          fk2_var_1810 = p_var_1806 + p4_var_1807
          fac000_var_1794 = fk0_var_1808 * fac00_var_1751(jl_var_1840, lay_var_1784)
          fac100_var_1795 = fk1_var_1809 * fac00_var_1751(jl_var_1840, lay_var_1784)
          fac200_var_1796 = fk2_var_1810 * fac00_var_1751(jl_var_1840, lay_var_1784)
          fac010_var_1797 = fk0_var_1808 * fac10_var_1753(jl_var_1840, lay_var_1784)
          fac110_var_1798 = fk1_var_1809 * fac10_var_1753(jl_var_1840, lay_var_1784)
          fac210_var_1799 = fk2_var_1810 * fac10_var_1753(jl_var_1840, lay_var_1784)
        ELSE
          fac000_var_1794 = (1.0D0 - fs_var_1818) * fac00_var_1751(jl_var_1840, lay_var_1784)
          fac010_var_1797 = (1.0D0 - fs_var_1818) * fac10_var_1753(jl_var_1840, lay_var_1784)
          fac100_var_1795 = fs_var_1818 * fac00_var_1751(jl_var_1840, lay_var_1784)
          fac110_var_1798 = fs_var_1818 * fac10_var_1753(jl_var_1840, lay_var_1784)
          fac200_var_1796 = 0.0D0
          fac210_var_1799 = 0.0D0
        END IF
        IF (specparm1_var_1823 .LT. 0.125D0) THEN
          p_var_1806 = fs1_var_1821 - 1.0D0
          p4_var_1807 = p_var_1806 ** 4
          fk0_var_1808 = p4_var_1807
          fk1_var_1809 = 1.0D0 - p_var_1806 - 2.0D0 * p4_var_1807
          fk2_var_1810 = p_var_1806 + p4_var_1807
          fac001_var_1800 = fk0_var_1808 * fac01_var_1752(jl_var_1840, lay_var_1784)
          fac101_var_1801 = fk1_var_1809 * fac01_var_1752(jl_var_1840, lay_var_1784)
          fac201_var_1802 = fk2_var_1810 * fac01_var_1752(jl_var_1840, lay_var_1784)
          fac011_var_1803 = fk0_var_1808 * fac11_var_1754(jl_var_1840, lay_var_1784)
          fac111_var_1804 = fk1_var_1809 * fac11_var_1754(jl_var_1840, lay_var_1784)
          fac211_var_1805 = fk2_var_1810 * fac11_var_1754(jl_var_1840, lay_var_1784)
        ELSE IF (specparm1_var_1823 .GT. 0.875D0) THEN
          p_var_1806 = - fs1_var_1821
          p4_var_1807 = p_var_1806 ** 4
          fk0_var_1808 = p4_var_1807
          fk1_var_1809 = 1.0D0 - p_var_1806 - 2.0D0 * p4_var_1807
          fk2_var_1810 = p_var_1806 + p4_var_1807
          fac001_var_1800 = fk0_var_1808 * fac01_var_1752(jl_var_1840, lay_var_1784)
          fac101_var_1801 = fk1_var_1809 * fac01_var_1752(jl_var_1840, lay_var_1784)
          fac201_var_1802 = fk2_var_1810 * fac01_var_1752(jl_var_1840, lay_var_1784)
          fac011_var_1803 = fk0_var_1808 * fac11_var_1754(jl_var_1840, lay_var_1784)
          fac111_var_1804 = fk1_var_1809 * fac11_var_1754(jl_var_1840, lay_var_1784)
          fac211_var_1805 = fk2_var_1810 * fac11_var_1754(jl_var_1840, lay_var_1784)
        ELSE
          fac001_var_1800 = (1.0D0 - fs1_var_1821) * fac01_var_1752(jl_var_1840, lay_var_1784)
          fac011_var_1803 = (1.0D0 - fs1_var_1821) * fac11_var_1754(jl_var_1840, lay_var_1784)
          fac101_var_1801 = fs1_var_1821 * fac01_var_1752(jl_var_1840, lay_var_1784)
          fac111_var_1804 = fs1_var_1821 * fac11_var_1754(jl_var_1840, lay_var_1784)
          fac201_var_1802 = 0.0D0
          fac211_var_1805 = 0.0D0
        END IF
        IF (specparm_var_1820 .LT. 0.125D0) THEN
          tau_major_var_1813(1 : ng7) = speccomb_var_1773 * (fac000_var_1794 * absa_var_188(ind0_var_1777, 1 : 12) + fac100_var_1795 * absa_var_188(ind0_var_1777 + 1, 1 : 12) + fac200_var_1796 * absa_var_188(ind0_var_1777 + 2, 1 : 12) + fac010_var_1797 * absa_var_188(ind0_var_1777 + 9, 1 : 12) + fac110_var_1798 * absa_var_188(ind0_var_1777 + 10, 1 : 12) + fac210_var_1799 * absa_var_188(ind0_var_1777 + 11, 1 : 12))
        ELSE IF (specparm_var_1820 .GT. 0.875D0) THEN
          tau_major_var_1813(1 : ng7) = speccomb_var_1773 * (fac200_var_1796 * absa_var_188(ind0_var_1777 - 1, 1 : 12) + fac100_var_1795 * absa_var_188(ind0_var_1777, 1 : 12) + fac000_var_1794 * absa_var_188(ind0_var_1777 + 1, 1 : 12) + fac210_var_1799 * absa_var_188(ind0_var_1777 + 8, 1 : 12) + fac110_var_1798 * absa_var_188(ind0_var_1777 + 9, 1 : 12) + fac010_var_1797 * absa_var_188(ind0_var_1777 + 10, 1 : 12))
        ELSE
          tau_major_var_1813(1 : ng7) = speccomb_var_1773 * (fac000_var_1794 * absa_var_188(ind0_var_1777, 1 : 12) + fac100_var_1795 * absa_var_188(ind0_var_1777 + 1, 1 : 12) + fac010_var_1797 * absa_var_188(ind0_var_1777 + 9, 1 : 12) + fac110_var_1798 * absa_var_188(ind0_var_1777 + 10, 1 : 12))
        END IF
        IF (specparm1_var_1823 .LT. 0.125D0) THEN
          tau_major1_var_1814(1 : ng7) = speccomb1_var_1774 * (fac001_var_1800 * absa_var_188(ind1_var_1778, 1 : 12) + fac101_var_1801 * absa_var_188(ind1_var_1778 + 1, 1 : 12) + fac201_var_1802 * absa_var_188(ind1_var_1778 + 2, 1 : 12) + fac011_var_1803 * absa_var_188(ind1_var_1778 + 9, 1 : 12) + fac111_var_1804 * absa_var_188(ind1_var_1778 + 10, 1 : 12) + fac211_var_1805 * absa_var_188(ind1_var_1778 + 11, 1 : 12))
        ELSE IF (specparm1_var_1823 .GT. 0.875D0) THEN
          tau_major1_var_1814(1 : ng7) = speccomb1_var_1774 * (fac201_var_1802 * absa_var_188(ind1_var_1778 - 1, 1 : 12) + fac101_var_1801 * absa_var_188(ind1_var_1778, 1 : 12) + fac001_var_1800 * absa_var_188(ind1_var_1778 + 1, 1 : 12) + fac211_var_1805 * absa_var_188(ind1_var_1778 + 8, 1 : 12) + fac111_var_1804 * absa_var_188(ind1_var_1778 + 9, 1 : 12) + fac011_var_1803 * absa_var_188(ind1_var_1778 + 10, 1 : 12))
        ELSE
          tau_major1_var_1814(1 : ng7) = speccomb1_var_1774 * (fac001_var_1800 * absa_var_188(ind1_var_1778, 1 : 12) + fac101_var_1801 * absa_var_188(ind1_var_1778 + 1, 1 : 12) + fac011_var_1803 * absa_var_188(ind1_var_1778 + 9, 1 : 12) + fac111_var_1804 * absa_var_188(ind1_var_1778 + 10, 1 : 12))
        END IF
        DO ig_var_1782 = 1, 12
          tauself_var_1812 = selffac_var_1764(jl_var_1840, lay_var_1784) * (selfref_var_190(inds_var_1779, ig_var_1782) + selffrac_var_1765(jl_var_1840, lay_var_1784) * (selfref_var_190(inds_var_1779 + 1, ig_var_1782) - selfref_var_190(inds_var_1779, ig_var_1782)))
          taufor_var_1811 = forfac_var_1770(jl_var_1840, lay_var_1784) * (forref_var_193(indf_var_1780, ig_var_1782) + forfrac_var_1769(jl_var_1840, lay_var_1784) * (forref_var_193(indf_var_1780 + 1, ig_var_1782) - forref_var_193(indf_var_1780, ig_var_1782)))
          co2m1_var_1815 = ka_mco2_var_191(jmco2_var_1787, indm_var_1781, ig_var_1782) + fmco2_var_1827 * (ka_mco2_var_191(jmco2_var_1787 + 1, indm_var_1781, ig_var_1782) - ka_mco2_var_191(jmco2_var_1787, indm_var_1781, ig_var_1782))
          co2m2_var_1816 = ka_mco2_var_191(jmco2_var_1787, indm_var_1781 + 1, ig_var_1782) + fmco2_var_1827 * (ka_mco2_var_191(jmco2_var_1787 + 1, indm_var_1781 + 1, ig_var_1782) - ka_mco2_var_191(jmco2_var_1787, indm_var_1781 + 1, ig_var_1782))
          absco2_var_1817 = co2m1_var_1815 + minorfrac_var_1771(jl_var_1840, lay_var_1784) * (co2m2_var_1816 - co2m1_var_1815)
          taug_var_1749(jl_var_1840, 76 + ig_var_1782, lay_var_1784) = tau_major_var_1813(ig_var_1782) + tau_major1_var_1814(ig_var_1782) + tauself_var_1812 + taufor_var_1811 + adjcolco2_var_1793 * absco2_var_1817
          fracs_var_1767(jl_var_1840, 76 + ig_var_1782, lay_var_1784) = fracrefa_var_186(ig_var_1782, jpl_var_1786) + fpl_var_1824 * (fracrefa_var_186(ig_var_1782, jpl_var_1786 + 1) - fracrefa_var_186(ig_var_1782, jpl_var_1786))
        END DO
      END DO
      ixc0_var_1837 = kfdia_var_1747 - kidia_var_1746 + 1 - ixc0_var_1837
      DO ixp_var_1838 = 1, ixc0_var_1837
        jl_var_1840 = ixhigh_var_1834(ixp_var_1838, lay_var_1784)
        chi_co2_var_1790 = colco2_var_1761(jl_var_1840, lay_var_1784) / (coldry_var_1762(jl_var_1840, lay_var_1784))
        ratco2_var_1791 = 1D+20 * chi_co2_var_1790 / chi_mls(2, jp_var_1755(jl_var_1840, lay_var_1784) + 1)
        IF (ratco2_var_1791 .GT. 3.0D0) THEN
          adjfac_var_1792 = 2.0D0 + (ratco2_var_1791 - 2.0D0) ** 0.79D0
          adjcolco2_var_1793 = adjfac_var_1792 * chi_mls(2, jp_var_1755(jl_var_1840, lay_var_1784) + 1) * coldry_var_1762(jl_var_1840, lay_var_1784) * 1D-20
        ELSE
          adjcolco2_var_1793 = colco2_var_1761(jl_var_1840, lay_var_1784)
        END IF
        ind0_var_1777 = ((jp_var_1755(jl_var_1840, lay_var_1784) - 13) * 5 + (jt_var_1756(jl_var_1840, lay_var_1784) - 1)) * nspb_var_217(7) + 1
        ind1_var_1778 = ((jp_var_1755(jl_var_1840, lay_var_1784) - 12) * 5 + (jt1_var_1757(jl_var_1840, lay_var_1784) - 1)) * nspb_var_217(7) + 1
        indm_var_1781 = indminor_var_1772(jl_var_1840, lay_var_1784)
        DO ig_var_1782 = 1, 12
          absco2_var_1817 = kb_mco2_var_192(indm_var_1781, ig_var_1782) + minorfrac_var_1771(jl_var_1840, lay_var_1784) * (kb_mco2_var_192(indm_var_1781 + 1, ig_var_1782) - kb_mco2_var_192(indm_var_1781, ig_var_1782))
          taug_var_1749(jl_var_1840, 76 + ig_var_1782, lay_var_1784) = colo3_var_1760(jl_var_1840, lay_var_1784) * (fac00_var_1751(jl_var_1840, lay_var_1784) * absb_var_189(ind0_var_1777, ig_var_1782) + fac10_var_1753(jl_var_1840, lay_var_1784) * absb_var_189(ind0_var_1777 + 1, ig_var_1782) + fac01_var_1752(jl_var_1840, lay_var_1784) * absb_var_189(ind1_var_1778, ig_var_1782) + fac11_var_1754(jl_var_1840, lay_var_1784) * absb_var_189(ind1_var_1778 + 1, ig_var_1782)) + adjcolco2_var_1793 * absco2_var_1817
          fracs_var_1767(jl_var_1840, 76 + ig_var_1782, lay_var_1784) = fracrefb_var_187(ig_var_1782)
        END DO
      END DO
      DO ixp_var_1838 = 1, ixc0_var_1837
        jl_var_1840 = ixhigh_var_1834(ixp_var_1838, lay_var_1784)
        taug_var_1749(jl_var_1840, 82, lay_var_1784) = taug_var_1749(jl_var_1840, 82, lay_var_1784) * 0.92D0
        taug_var_1749(jl_var_1840, 83, lay_var_1784) = taug_var_1749(jl_var_1840, 83, lay_var_1784) * 0.88D0
        taug_var_1749(jl_var_1840, 84, lay_var_1784) = taug_var_1749(jl_var_1840, 84, lay_var_1784) * 1.07D0
        taug_var_1749(jl_var_1840, 85, lay_var_1784) = taug_var_1749(jl_var_1840, 85, lay_var_1784) * 1.1D0
        taug_var_1749(jl_var_1840, 86, lay_var_1784) = taug_var_1749(jl_var_1840, 86, lay_var_1784) * 0.99D0
        taug_var_1749(jl_var_1840, 87, lay_var_1784) = taug_var_1749(jl_var_1840, 87, lay_var_1784) * 0.855D0
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol7
SUBROUTINE rrtm_taumol3(kidia_var_1841, kfdia_var_1842, klev_var_1843, taug_var_1844, p_tauaerl_var_1845, fac00_var_1846, fac01_var_1847, fac10_var_1848, fac11_var_1849, forfac_var_1850, forfrac_var_1867, indfor_var_1866, jp_var_1851, jt_var_1852, jt1_var_1853, oneminus_var_1854, colh2o_var_1855, colco2_var_1856, coln2o_var_1857, coldry_var_1858, laytrop_var_1859, selffac_var_1860, selffrac_var_1861, indself_var_1862, fracs_var_1863, rat_h2oco2_var_1864, rat_h2oco2_1_var_1865, minorfrac_var_1868, indminor_var_1869)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng3
  USE yoerrta3, ONLY: absa_var_163, absb_var_164, forref_var_166, fracrefa_var_159, fracrefb_var_160, ka_mn2o_var_161, kb_mn2o_var_162, selfref_var_165
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1841
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1842
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1843
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1844(kidia_var_1841 : kfdia_var_1842, 140, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1845(kidia_var_1841 : kfdia_var_1842, klev_var_1843, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1846(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1847(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1848(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1849(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1850(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1851(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1852(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1853(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1854
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1855(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1856(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_1857(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1858(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1859(kidia_var_1841 : kfdia_var_1842)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1860(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1861(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1862(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1863(kidia_var_1841 : kfdia_var_1842, 140, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1864(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1865(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1866(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1867(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1868(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1869(kidia_var_1841 : kfdia_var_1842, klev_var_1843)
  REAL(KIND = 8) :: speccomb_var_1870, speccomb1_var_1871, speccomb_mn2o_var_1872, speccomb_planck_var_1873
  REAL(KIND = 8) :: refrat_planck_a_var_1874, refrat_planck_b_var_1875, refrat_m_a_var_1876, refrat_m_b
  INTEGER(KIND = 4) :: ind0_var_1877, ind1_var_1878, inds_var_1879, indf_var_1880, indm_var_1881
  INTEGER(KIND = 4) :: ig_var_1882, js_var_1883, lay_var_1884, js1_var_1885, jmn2o_var_1886, jpl_var_1887
  REAL(KIND = 8) :: fs_var_1888, specmult_var_1889, specparm_var_1890, fs1_var_1891, specmult1_var_1892, specparm1_var_1893, fmn2o_var_1894, specmult_mn2o_var_1895, specparm_mn2o_var_1896, fpl_var_1897, specmult_planck_var_1898, specparm_planck_var_1899
  REAL(KIND = 8) :: adjfac_var_1900, adjcoln2o_var_1901, ratn2o_var_1902, chi_n2o_var_1903
  REAL(KIND = 8) :: fac000_var_1904, fac100_var_1905, fac200_var_1906, fac010_var_1907, fac110_var_1908, fac210_var_1909, fac001_var_1910, fac101_var_1911, fac201_var_1912, fac011_var_1913, fac111_var_1914, fac211_var_1915
  REAL(KIND = 8) :: p_var_1916, p4_var_1917, fk0_var_1918, fk1_var_1919, fk2_var_1920
  REAL(KIND = 8) :: taufor_var_1921, tauself_var_1922, n2om1_var_1923, n2om2_var_1924, absn2o_var_1925, tau_major_var_1926(16), tau_major1_var_1927(16)
  INTEGER(KIND = 4) :: laytrop_min_var_1928, laytrop_max_var_1929
  INTEGER(KIND = 4) :: ixc_var_1930(klev_var_1843), ixlow_var_1931(kfdia_var_1842, klev_var_1843), ixhigh_var_1932(kfdia_var_1842, klev_var_1843)
  INTEGER(KIND = 4) :: ich_var_1933, icl_var_1934, ixc0_var_1935, ixp_var_1936, jc_var_1937, jl_var_1938
  laytrop_min_var_1928 = MINVAL(laytrop_var_1859)
  laytrop_max_var_1929 = MAXVAL(laytrop_var_1859)
  ixlow_var_1931 = 0
  ixhigh_var_1932 = 0
  ixc_var_1930 = 0
  DO lay_var_1884 = laytrop_min_var_1928 + 1, laytrop_max_var_1929
    icl_var_1934 = 0
    ich_var_1933 = 0
    DO jc_var_1937 = kidia_var_1841, kfdia_var_1842
      IF (lay_var_1884 <= laytrop_var_1859(jc_var_1937)) THEN
        icl_var_1934 = icl_var_1934 + 1
        ixlow_var_1931(icl_var_1934, lay_var_1884) = jc_var_1937
      ELSE
        ich_var_1933 = ich_var_1933 + 1
        ixhigh_var_1932(ich_var_1933, lay_var_1884) = jc_var_1937
      END IF
    END DO
    ixc_var_1930(lay_var_1884) = icl_var_1934
  END DO
  refrat_planck_a_var_1874 = chi_mls(1, 9) / chi_mls(2, 9)
  refrat_planck_b_var_1875 = chi_mls(1, 13) / chi_mls(2, 13)
  refrat_m_a_var_1876 = chi_mls(1, 3) / chi_mls(2, 3)
  refrat_m_b = chi_mls(1, 13) / chi_mls(2, 13)
  DO lay_var_1884 = 1, laytrop_min_var_1928
    DO jl_var_1938 = kidia_var_1841, kfdia_var_1842
      speccomb_var_1870 = colh2o_var_1855(jl_var_1938, lay_var_1884) + rat_h2oco2_var_1864(jl_var_1938, lay_var_1884) * colco2_var_1856(jl_var_1938, lay_var_1884)
      specparm_var_1890 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_var_1870, oneminus_var_1854)
      specmult_var_1889 = 8.0D0 * (specparm_var_1890)
      js_var_1883 = 1 + INT(specmult_var_1889)
      fs_var_1888 = ((specmult_var_1889) - AINT((specmult_var_1889)))
      speccomb1_var_1871 = colh2o_var_1855(jl_var_1938, lay_var_1884) + rat_h2oco2_1_var_1865(jl_var_1938, lay_var_1884) * colco2_var_1856(jl_var_1938, lay_var_1884)
      specparm1_var_1893 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb1_var_1871, oneminus_var_1854)
      specmult1_var_1892 = 8.0D0 * (specparm1_var_1893)
      js1_var_1885 = 1 + INT(specmult1_var_1892)
      fs1_var_1891 = ((specmult1_var_1892) - AINT((specmult1_var_1892)))
      speccomb_mn2o_var_1872 = colh2o_var_1855(jl_var_1938, lay_var_1884) + refrat_m_a_var_1876 * colco2_var_1856(jl_var_1938, lay_var_1884)
      specparm_mn2o_var_1896 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_mn2o_var_1872, oneminus_var_1854)
      specmult_mn2o_var_1895 = 8.0D0 * specparm_mn2o_var_1896
      jmn2o_var_1886 = 1 + INT(specmult_mn2o_var_1895)
      fmn2o_var_1894 = ((specmult_mn2o_var_1895) - AINT((specmult_mn2o_var_1895)))
      chi_n2o_var_1903 = coln2o_var_1857(jl_var_1938, lay_var_1884) / coldry_var_1858(jl_var_1938, lay_var_1884)
      ratn2o_var_1902 = 1D+20 * chi_n2o_var_1903 / chi_mls(4, jp_var_1851(jl_var_1938, lay_var_1884) + 1)
      IF (ratn2o_var_1902 .GT. 1.5D0) THEN
        adjfac_var_1900 = 0.5D0 + (ratn2o_var_1902 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1901 = adjfac_var_1900 * chi_mls(4, jp_var_1851(jl_var_1938, lay_var_1884) + 1) * coldry_var_1858(jl_var_1938, lay_var_1884) * 1D-20
      ELSE
        adjcoln2o_var_1901 = coln2o_var_1857(jl_var_1938, lay_var_1884)
      END IF
      speccomb_planck_var_1873 = colh2o_var_1855(jl_var_1938, lay_var_1884) + refrat_planck_a_var_1874 * colco2_var_1856(jl_var_1938, lay_var_1884)
      specparm_planck_var_1899 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_planck_var_1873, oneminus_var_1854)
      specmult_planck_var_1898 = 8.0D0 * specparm_planck_var_1899
      jpl_var_1887 = 1 + INT(specmult_planck_var_1898)
      fpl_var_1897 = ((specmult_planck_var_1898) - AINT((specmult_planck_var_1898)))
      ind0_var_1877 = ((jp_var_1851(jl_var_1938, lay_var_1884) - 1) * 5 + (jt_var_1852(jl_var_1938, lay_var_1884) - 1)) * nspa_var_216(3) + js_var_1883
      ind1_var_1878 = (jp_var_1851(jl_var_1938, lay_var_1884) * 5 + (jt1_var_1853(jl_var_1938, lay_var_1884) - 1)) * nspa_var_216(3) + js1_var_1885
      inds_var_1879 = indself_var_1862(jl_var_1938, lay_var_1884)
      indf_var_1880 = indfor_var_1866(jl_var_1938, lay_var_1884)
      indm_var_1881 = indminor_var_1869(jl_var_1938, lay_var_1884)
      IF (specparm_var_1890 .LT. 0.125D0) THEN
        p_var_1916 = fs_var_1888 - 1.0D0
        p4_var_1917 = p_var_1916 ** 4
        fk0_var_1918 = p4_var_1917
        fk1_var_1919 = 1.0D0 - p_var_1916 - 2.0D0 * p4_var_1917
        fk2_var_1920 = p_var_1916 + p4_var_1917
        fac000_var_1904 = fk0_var_1918 * fac00_var_1846(jl_var_1938, lay_var_1884)
        fac100_var_1905 = fk1_var_1919 * fac00_var_1846(jl_var_1938, lay_var_1884)
        fac200_var_1906 = fk2_var_1920 * fac00_var_1846(jl_var_1938, lay_var_1884)
        fac010_var_1907 = fk0_var_1918 * fac10_var_1848(jl_var_1938, lay_var_1884)
        fac110_var_1908 = fk1_var_1919 * fac10_var_1848(jl_var_1938, lay_var_1884)
        fac210_var_1909 = fk2_var_1920 * fac10_var_1848(jl_var_1938, lay_var_1884)
      ELSE IF (specparm_var_1890 .GT. 0.875D0) THEN
        p_var_1916 = - fs_var_1888
        p4_var_1917 = p_var_1916 ** 4
        fk0_var_1918 = p4_var_1917
        fk1_var_1919 = 1.0D0 - p_var_1916 - 2.0D0 * p4_var_1917
        fk2_var_1920 = p_var_1916 + p4_var_1917
        fac000_var_1904 = fk0_var_1918 * fac00_var_1846(jl_var_1938, lay_var_1884)
        fac100_var_1905 = fk1_var_1919 * fac00_var_1846(jl_var_1938, lay_var_1884)
        fac200_var_1906 = fk2_var_1920 * fac00_var_1846(jl_var_1938, lay_var_1884)
        fac010_var_1907 = fk0_var_1918 * fac10_var_1848(jl_var_1938, lay_var_1884)
        fac110_var_1908 = fk1_var_1919 * fac10_var_1848(jl_var_1938, lay_var_1884)
        fac210_var_1909 = fk2_var_1920 * fac10_var_1848(jl_var_1938, lay_var_1884)
      ELSE
        fac000_var_1904 = (1.0D0 - fs_var_1888) * fac00_var_1846(jl_var_1938, lay_var_1884)
        fac010_var_1907 = (1.0D0 - fs_var_1888) * fac10_var_1848(jl_var_1938, lay_var_1884)
        fac100_var_1905 = fs_var_1888 * fac00_var_1846(jl_var_1938, lay_var_1884)
        fac110_var_1908 = fs_var_1888 * fac10_var_1848(jl_var_1938, lay_var_1884)
        fac200_var_1906 = 0.0D0
        fac210_var_1909 = 0.0D0
      END IF
      IF (specparm1_var_1893 .LT. 0.125D0) THEN
        p_var_1916 = fs1_var_1891 - 1.0D0
        p4_var_1917 = p_var_1916 ** 4
        fk0_var_1918 = p4_var_1917
        fk1_var_1919 = 1.0D0 - p_var_1916 - 2.0D0 * p4_var_1917
        fk2_var_1920 = p_var_1916 + p4_var_1917
        fac001_var_1910 = fk0_var_1918 * fac01_var_1847(jl_var_1938, lay_var_1884)
        fac101_var_1911 = fk1_var_1919 * fac01_var_1847(jl_var_1938, lay_var_1884)
        fac201_var_1912 = fk2_var_1920 * fac01_var_1847(jl_var_1938, lay_var_1884)
        fac011_var_1913 = fk0_var_1918 * fac11_var_1849(jl_var_1938, lay_var_1884)
        fac111_var_1914 = fk1_var_1919 * fac11_var_1849(jl_var_1938, lay_var_1884)
        fac211_var_1915 = fk2_var_1920 * fac11_var_1849(jl_var_1938, lay_var_1884)
      ELSE IF (specparm1_var_1893 .GT. 0.875D0) THEN
        p_var_1916 = - fs1_var_1891
        p4_var_1917 = p_var_1916 ** 4
        fk0_var_1918 = p4_var_1917
        fk1_var_1919 = 1.0D0 - p_var_1916 - 2.0D0 * p4_var_1917
        fk2_var_1920 = p_var_1916 + p4_var_1917
        fac001_var_1910 = fk0_var_1918 * fac01_var_1847(jl_var_1938, lay_var_1884)
        fac101_var_1911 = fk1_var_1919 * fac01_var_1847(jl_var_1938, lay_var_1884)
        fac201_var_1912 = fk2_var_1920 * fac01_var_1847(jl_var_1938, lay_var_1884)
        fac011_var_1913 = fk0_var_1918 * fac11_var_1849(jl_var_1938, lay_var_1884)
        fac111_var_1914 = fk1_var_1919 * fac11_var_1849(jl_var_1938, lay_var_1884)
        fac211_var_1915 = fk2_var_1920 * fac11_var_1849(jl_var_1938, lay_var_1884)
      ELSE
        fac001_var_1910 = (1.0D0 - fs1_var_1891) * fac01_var_1847(jl_var_1938, lay_var_1884)
        fac011_var_1913 = (1.0D0 - fs1_var_1891) * fac11_var_1849(jl_var_1938, lay_var_1884)
        fac101_var_1911 = fs1_var_1891 * fac01_var_1847(jl_var_1938, lay_var_1884)
        fac111_var_1914 = fs1_var_1891 * fac11_var_1849(jl_var_1938, lay_var_1884)
        fac201_var_1912 = 0.0D0
        fac211_var_1915 = 0.0D0
      END IF
      IF (specparm_var_1890 .LT. 0.125D0) THEN
        tau_major_var_1926(1 : ng3) = speccomb_var_1870 * (fac000_var_1904 * absa_var_163(ind0_var_1877, 1 : 16) + fac100_var_1905 * absa_var_163(ind0_var_1877 + 1, 1 : 16) + fac200_var_1906 * absa_var_163(ind0_var_1877 + 2, 1 : 16) + fac010_var_1907 * absa_var_163(ind0_var_1877 + 9, 1 : 16) + fac110_var_1908 * absa_var_163(ind0_var_1877 + 10, 1 : 16) + fac210_var_1909 * absa_var_163(ind0_var_1877 + 11, 1 : 16))
      ELSE IF (specparm_var_1890 .GT. 0.875D0) THEN
        tau_major_var_1926(1 : ng3) = speccomb_var_1870 * (fac200_var_1906 * absa_var_163(ind0_var_1877 - 1, 1 : 16) + fac100_var_1905 * absa_var_163(ind0_var_1877, 1 : 16) + fac000_var_1904 * absa_var_163(ind0_var_1877 + 1, 1 : 16) + fac210_var_1909 * absa_var_163(ind0_var_1877 + 8, 1 : 16) + fac110_var_1908 * absa_var_163(ind0_var_1877 + 9, 1 : 16) + fac010_var_1907 * absa_var_163(ind0_var_1877 + 10, 1 : 16))
      ELSE
        tau_major_var_1926(1 : ng3) = speccomb_var_1870 * (fac000_var_1904 * absa_var_163(ind0_var_1877, 1 : 16) + fac100_var_1905 * absa_var_163(ind0_var_1877 + 1, 1 : 16) + fac010_var_1907 * absa_var_163(ind0_var_1877 + 9, 1 : 16) + fac110_var_1908 * absa_var_163(ind0_var_1877 + 10, 1 : 16))
      END IF
      IF (specparm1_var_1893 .LT. 0.125D0) THEN
        tau_major1_var_1927(1 : ng3) = speccomb1_var_1871 * (fac001_var_1910 * absa_var_163(ind1_var_1878, 1 : 16) + fac101_var_1911 * absa_var_163(ind1_var_1878 + 1, 1 : 16) + fac201_var_1912 * absa_var_163(ind1_var_1878 + 2, 1 : 16) + fac011_var_1913 * absa_var_163(ind1_var_1878 + 9, 1 : 16) + fac111_var_1914 * absa_var_163(ind1_var_1878 + 10, 1 : 16) + fac211_var_1915 * absa_var_163(ind1_var_1878 + 11, 1 : 16))
      ELSE IF (specparm1_var_1893 .GT. 0.875D0) THEN
        tau_major1_var_1927(1 : ng3) = speccomb1_var_1871 * (fac201_var_1912 * absa_var_163(ind1_var_1878 - 1, 1 : 16) + fac101_var_1911 * absa_var_163(ind1_var_1878, 1 : 16) + fac001_var_1910 * absa_var_163(ind1_var_1878 + 1, 1 : 16) + fac211_var_1915 * absa_var_163(ind1_var_1878 + 8, 1 : 16) + fac111_var_1914 * absa_var_163(ind1_var_1878 + 9, 1 : 16) + fac011_var_1913 * absa_var_163(ind1_var_1878 + 10, 1 : 16))
      ELSE
        tau_major1_var_1927(1 : ng3) = speccomb1_var_1871 * (fac001_var_1910 * absa_var_163(ind1_var_1878, 1 : 16) + fac101_var_1911 * absa_var_163(ind1_var_1878 + 1, 1 : 16) + fac011_var_1913 * absa_var_163(ind1_var_1878 + 9, 1 : 16) + fac111_var_1914 * absa_var_163(ind1_var_1878 + 10, 1 : 16))
      END IF
      DO ig_var_1882 = 1, 16
        tauself_var_1922 = selffac_var_1860(jl_var_1938, lay_var_1884) * (selfref_var_165(inds_var_1879, ig_var_1882) + selffrac_var_1861(jl_var_1938, lay_var_1884) * (selfref_var_165(inds_var_1879 + 1, ig_var_1882) - selfref_var_165(inds_var_1879, ig_var_1882)))
        taufor_var_1921 = forfac_var_1850(jl_var_1938, lay_var_1884) * (forref_var_166(indf_var_1880, ig_var_1882) + forfrac_var_1867(jl_var_1938, lay_var_1884) * (forref_var_166(indf_var_1880 + 1, ig_var_1882) - forref_var_166(indf_var_1880, ig_var_1882)))
        n2om1_var_1923 = ka_mn2o_var_161(jmn2o_var_1886, indm_var_1881, ig_var_1882) + fmn2o_var_1894 * (ka_mn2o_var_161(jmn2o_var_1886 + 1, indm_var_1881, ig_var_1882) - ka_mn2o_var_161(jmn2o_var_1886, indm_var_1881, ig_var_1882))
        n2om2_var_1924 = ka_mn2o_var_161(jmn2o_var_1886, indm_var_1881 + 1, ig_var_1882) + fmn2o_var_1894 * (ka_mn2o_var_161(jmn2o_var_1886 + 1, indm_var_1881 + 1, ig_var_1882) - ka_mn2o_var_161(jmn2o_var_1886, indm_var_1881 + 1, ig_var_1882))
        absn2o_var_1925 = n2om1_var_1923 + minorfrac_var_1868(jl_var_1938, lay_var_1884) * (n2om2_var_1924 - n2om1_var_1923)
        taug_var_1844(jl_var_1938, 22 + ig_var_1882, lay_var_1884) = tau_major_var_1926(ig_var_1882) + tau_major1_var_1927(ig_var_1882) + tauself_var_1922 + taufor_var_1921 + adjcoln2o_var_1901 * absn2o_var_1925
        fracs_var_1863(jl_var_1938, 22 + ig_var_1882, lay_var_1884) = fracrefa_var_159(ig_var_1882, jpl_var_1887) + fpl_var_1897 * (fracrefa_var_159(ig_var_1882, jpl_var_1887 + 1) - fracrefa_var_159(ig_var_1882, jpl_var_1887))
      END DO
    END DO
  END DO
  DO lay_var_1884 = laytrop_max_var_1929 + 1, klev_var_1843
    DO jl_var_1938 = kidia_var_1841, kfdia_var_1842
      speccomb_var_1870 = colh2o_var_1855(jl_var_1938, lay_var_1884) + rat_h2oco2_var_1864(jl_var_1938, lay_var_1884) * colco2_var_1856(jl_var_1938, lay_var_1884)
      specparm_var_1890 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_var_1870, oneminus_var_1854)
      specmult_var_1889 = 4.0D0 * (specparm_var_1890)
      js_var_1883 = 1 + INT(specmult_var_1889)
      fs_var_1888 = ((specmult_var_1889) - AINT((specmult_var_1889)))
      speccomb1_var_1871 = colh2o_var_1855(jl_var_1938, lay_var_1884) + rat_h2oco2_1_var_1865(jl_var_1938, lay_var_1884) * colco2_var_1856(jl_var_1938, lay_var_1884)
      specparm1_var_1893 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb1_var_1871, oneminus_var_1854)
      specmult1_var_1892 = 4.0D0 * (specparm1_var_1893)
      js1_var_1885 = 1 + INT(specmult1_var_1892)
      fs1_var_1891 = ((specmult1_var_1892) - AINT((specmult1_var_1892)))
      fac000_var_1904 = (1.0D0 - fs_var_1888) * fac00_var_1846(jl_var_1938, lay_var_1884)
      fac010_var_1907 = (1.0D0 - fs_var_1888) * fac10_var_1848(jl_var_1938, lay_var_1884)
      fac100_var_1905 = fs_var_1888 * fac00_var_1846(jl_var_1938, lay_var_1884)
      fac110_var_1908 = fs_var_1888 * fac10_var_1848(jl_var_1938, lay_var_1884)
      fac001_var_1910 = (1.0D0 - fs1_var_1891) * fac01_var_1847(jl_var_1938, lay_var_1884)
      fac011_var_1913 = (1.0D0 - fs1_var_1891) * fac11_var_1849(jl_var_1938, lay_var_1884)
      fac101_var_1911 = fs1_var_1891 * fac01_var_1847(jl_var_1938, lay_var_1884)
      fac111_var_1914 = fs1_var_1891 * fac11_var_1849(jl_var_1938, lay_var_1884)
      speccomb_mn2o_var_1872 = colh2o_var_1855(jl_var_1938, lay_var_1884) + refrat_m_b * colco2_var_1856(jl_var_1938, lay_var_1884)
      specparm_mn2o_var_1896 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_mn2o_var_1872, oneminus_var_1854)
      specmult_mn2o_var_1895 = 4.0D0 * specparm_mn2o_var_1896
      jmn2o_var_1886 = 1 + INT(specmult_mn2o_var_1895)
      fmn2o_var_1894 = ((specmult_mn2o_var_1895) - AINT((specmult_mn2o_var_1895)))
      chi_n2o_var_1903 = coln2o_var_1857(jl_var_1938, lay_var_1884) / coldry_var_1858(jl_var_1938, lay_var_1884)
      ratn2o_var_1902 = 1D+20 * chi_n2o_var_1903 / chi_mls(4, jp_var_1851(jl_var_1938, lay_var_1884) + 1)
      IF (ratn2o_var_1902 .GT. 1.5D0) THEN
        adjfac_var_1900 = 0.5D0 + (ratn2o_var_1902 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1901 = adjfac_var_1900 * chi_mls(4, jp_var_1851(jl_var_1938, lay_var_1884) + 1) * coldry_var_1858(jl_var_1938, lay_var_1884) * 1D-20
      ELSE
        adjcoln2o_var_1901 = coln2o_var_1857(jl_var_1938, lay_var_1884)
      END IF
      speccomb_planck_var_1873 = colh2o_var_1855(jl_var_1938, lay_var_1884) + refrat_planck_b_var_1875 * colco2_var_1856(jl_var_1938, lay_var_1884)
      specparm_planck_var_1899 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_planck_var_1873, oneminus_var_1854)
      specmult_planck_var_1898 = 4.0D0 * specparm_planck_var_1899
      jpl_var_1887 = 1 + INT(specmult_planck_var_1898)
      fpl_var_1897 = ((specmult_planck_var_1898) - AINT((specmult_planck_var_1898)))
      ind0_var_1877 = ((jp_var_1851(jl_var_1938, lay_var_1884) - 13) * 5 + (jt_var_1852(jl_var_1938, lay_var_1884) - 1)) * nspb_var_217(3) + js_var_1883
      ind1_var_1878 = ((jp_var_1851(jl_var_1938, lay_var_1884) - 12) * 5 + (jt1_var_1853(jl_var_1938, lay_var_1884) - 1)) * nspb_var_217(3) + js1_var_1885
      indf_var_1880 = indfor_var_1866(jl_var_1938, lay_var_1884)
      indm_var_1881 = indminor_var_1869(jl_var_1938, lay_var_1884)
      DO ig_var_1882 = 1, 16
        taufor_var_1921 = forfac_var_1850(jl_var_1938, lay_var_1884) * (forref_var_166(indf_var_1880, ig_var_1882) + forfrac_var_1867(jl_var_1938, lay_var_1884) * (forref_var_166(indf_var_1880 + 1, ig_var_1882) - forref_var_166(indf_var_1880, ig_var_1882)))
        n2om1_var_1923 = kb_mn2o_var_162(jmn2o_var_1886, indm_var_1881, ig_var_1882) + fmn2o_var_1894 * (kb_mn2o_var_162(jmn2o_var_1886 + 1, indm_var_1881, ig_var_1882) - kb_mn2o_var_162(jmn2o_var_1886, indm_var_1881, ig_var_1882))
        n2om2_var_1924 = kb_mn2o_var_162(jmn2o_var_1886, indm_var_1881 + 1, ig_var_1882) + fmn2o_var_1894 * (kb_mn2o_var_162(jmn2o_var_1886 + 1, indm_var_1881 + 1, ig_var_1882) - kb_mn2o_var_162(jmn2o_var_1886, indm_var_1881 + 1, ig_var_1882))
        absn2o_var_1925 = n2om1_var_1923 + minorfrac_var_1868(jl_var_1938, lay_var_1884) * (n2om2_var_1924 - n2om1_var_1923)
        taug_var_1844(jl_var_1938, 22 + ig_var_1882, lay_var_1884) = speccomb_var_1870 * (fac000_var_1904 * absb_var_164(ind0_var_1877, ig_var_1882) + fac100_var_1905 * absb_var_164(ind0_var_1877 + 1, ig_var_1882) + fac010_var_1907 * absb_var_164(ind0_var_1877 + 5, ig_var_1882) + fac110_var_1908 * absb_var_164(ind0_var_1877 + 6, ig_var_1882)) + speccomb1_var_1871 * (fac001_var_1910 * absb_var_164(ind1_var_1878, ig_var_1882) + fac101_var_1911 * absb_var_164(ind1_var_1878 + 1, ig_var_1882) + fac011_var_1913 * absb_var_164(ind1_var_1878 + 5, ig_var_1882) + fac111_var_1914 * absb_var_164(ind1_var_1878 + 6, ig_var_1882)) + taufor_var_1921 + adjcoln2o_var_1901 * absn2o_var_1925
        fracs_var_1863(jl_var_1938, 22 + ig_var_1882, lay_var_1884) = fracrefb_var_160(ig_var_1882, jpl_var_1887) + fpl_var_1897 * (fracrefb_var_160(ig_var_1882, jpl_var_1887 + 1) - fracrefb_var_160(ig_var_1882, jpl_var_1887))
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1929 /= laytrop_min_var_1928) THEN
    DO lay_var_1884 = laytrop_min_var_1928 + 1, laytrop_max_var_1929
      ixc0_var_1935 = ixc_var_1930(lay_var_1884)
      DO ixp_var_1936 = 1, ixc0_var_1935
        jl_var_1938 = ixlow_var_1931(ixp_var_1936, lay_var_1884)
        speccomb_var_1870 = colh2o_var_1855(jl_var_1938, lay_var_1884) + rat_h2oco2_var_1864(jl_var_1938, lay_var_1884) * colco2_var_1856(jl_var_1938, lay_var_1884)
        specparm_var_1890 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_var_1870, oneminus_var_1854)
        specmult_var_1889 = 8.0D0 * (specparm_var_1890)
        js_var_1883 = 1 + INT(specmult_var_1889)
        fs_var_1888 = ((specmult_var_1889) - AINT((specmult_var_1889)))
        speccomb1_var_1871 = colh2o_var_1855(jl_var_1938, lay_var_1884) + rat_h2oco2_1_var_1865(jl_var_1938, lay_var_1884) * colco2_var_1856(jl_var_1938, lay_var_1884)
        specparm1_var_1893 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb1_var_1871, oneminus_var_1854)
        specmult1_var_1892 = 8.0D0 * (specparm1_var_1893)
        js1_var_1885 = 1 + INT(specmult1_var_1892)
        fs1_var_1891 = ((specmult1_var_1892) - AINT((specmult1_var_1892)))
        speccomb_mn2o_var_1872 = colh2o_var_1855(jl_var_1938, lay_var_1884) + refrat_m_a_var_1876 * colco2_var_1856(jl_var_1938, lay_var_1884)
        specparm_mn2o_var_1896 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_mn2o_var_1872, oneminus_var_1854)
        specmult_mn2o_var_1895 = 8.0D0 * specparm_mn2o_var_1896
        jmn2o_var_1886 = 1 + INT(specmult_mn2o_var_1895)
        fmn2o_var_1894 = ((specmult_mn2o_var_1895) - AINT((specmult_mn2o_var_1895)))
        chi_n2o_var_1903 = coln2o_var_1857(jl_var_1938, lay_var_1884) / coldry_var_1858(jl_var_1938, lay_var_1884)
        ratn2o_var_1902 = 1D+20 * chi_n2o_var_1903 / chi_mls(4, jp_var_1851(jl_var_1938, lay_var_1884) + 1)
        IF (ratn2o_var_1902 .GT. 1.5D0) THEN
          adjfac_var_1900 = 0.5D0 + (ratn2o_var_1902 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1901 = adjfac_var_1900 * chi_mls(4, jp_var_1851(jl_var_1938, lay_var_1884) + 1) * coldry_var_1858(jl_var_1938, lay_var_1884) * 1D-20
        ELSE
          adjcoln2o_var_1901 = coln2o_var_1857(jl_var_1938, lay_var_1884)
        END IF
        speccomb_planck_var_1873 = colh2o_var_1855(jl_var_1938, lay_var_1884) + refrat_planck_a_var_1874 * colco2_var_1856(jl_var_1938, lay_var_1884)
        specparm_planck_var_1899 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_planck_var_1873, oneminus_var_1854)
        specmult_planck_var_1898 = 8.0D0 * specparm_planck_var_1899
        jpl_var_1887 = 1 + INT(specmult_planck_var_1898)
        fpl_var_1897 = ((specmult_planck_var_1898) - AINT((specmult_planck_var_1898)))
        ind0_var_1877 = ((jp_var_1851(jl_var_1938, lay_var_1884) - 1) * 5 + (jt_var_1852(jl_var_1938, lay_var_1884) - 1)) * nspa_var_216(3) + js_var_1883
        ind1_var_1878 = (jp_var_1851(jl_var_1938, lay_var_1884) * 5 + (jt1_var_1853(jl_var_1938, lay_var_1884) - 1)) * nspa_var_216(3) + js1_var_1885
        inds_var_1879 = indself_var_1862(jl_var_1938, lay_var_1884)
        indf_var_1880 = indfor_var_1866(jl_var_1938, lay_var_1884)
        indm_var_1881 = indminor_var_1869(jl_var_1938, lay_var_1884)
        IF (specparm_var_1890 .LT. 0.125D0) THEN
          p_var_1916 = fs_var_1888 - 1.0D0
          p4_var_1917 = p_var_1916 ** 4
          fk0_var_1918 = p4_var_1917
          fk1_var_1919 = 1.0D0 - p_var_1916 - 2.0D0 * p4_var_1917
          fk2_var_1920 = p_var_1916 + p4_var_1917
          fac000_var_1904 = fk0_var_1918 * fac00_var_1846(jl_var_1938, lay_var_1884)
          fac100_var_1905 = fk1_var_1919 * fac00_var_1846(jl_var_1938, lay_var_1884)
          fac200_var_1906 = fk2_var_1920 * fac00_var_1846(jl_var_1938, lay_var_1884)
          fac010_var_1907 = fk0_var_1918 * fac10_var_1848(jl_var_1938, lay_var_1884)
          fac110_var_1908 = fk1_var_1919 * fac10_var_1848(jl_var_1938, lay_var_1884)
          fac210_var_1909 = fk2_var_1920 * fac10_var_1848(jl_var_1938, lay_var_1884)
        ELSE IF (specparm_var_1890 .GT. 0.875D0) THEN
          p_var_1916 = - fs_var_1888
          p4_var_1917 = p_var_1916 ** 4
          fk0_var_1918 = p4_var_1917
          fk1_var_1919 = 1.0D0 - p_var_1916 - 2.0D0 * p4_var_1917
          fk2_var_1920 = p_var_1916 + p4_var_1917
          fac000_var_1904 = fk0_var_1918 * fac00_var_1846(jl_var_1938, lay_var_1884)
          fac100_var_1905 = fk1_var_1919 * fac00_var_1846(jl_var_1938, lay_var_1884)
          fac200_var_1906 = fk2_var_1920 * fac00_var_1846(jl_var_1938, lay_var_1884)
          fac010_var_1907 = fk0_var_1918 * fac10_var_1848(jl_var_1938, lay_var_1884)
          fac110_var_1908 = fk1_var_1919 * fac10_var_1848(jl_var_1938, lay_var_1884)
          fac210_var_1909 = fk2_var_1920 * fac10_var_1848(jl_var_1938, lay_var_1884)
        ELSE
          fac000_var_1904 = (1.0D0 - fs_var_1888) * fac00_var_1846(jl_var_1938, lay_var_1884)
          fac010_var_1907 = (1.0D0 - fs_var_1888) * fac10_var_1848(jl_var_1938, lay_var_1884)
          fac100_var_1905 = fs_var_1888 * fac00_var_1846(jl_var_1938, lay_var_1884)
          fac110_var_1908 = fs_var_1888 * fac10_var_1848(jl_var_1938, lay_var_1884)
          fac200_var_1906 = 0.0D0
          fac210_var_1909 = 0.0D0
        END IF
        IF (specparm1_var_1893 .LT. 0.125D0) THEN
          p_var_1916 = fs1_var_1891 - 1.0D0
          p4_var_1917 = p_var_1916 ** 4
          fk0_var_1918 = p4_var_1917
          fk1_var_1919 = 1.0D0 - p_var_1916 - 2.0D0 * p4_var_1917
          fk2_var_1920 = p_var_1916 + p4_var_1917
          fac001_var_1910 = fk0_var_1918 * fac01_var_1847(jl_var_1938, lay_var_1884)
          fac101_var_1911 = fk1_var_1919 * fac01_var_1847(jl_var_1938, lay_var_1884)
          fac201_var_1912 = fk2_var_1920 * fac01_var_1847(jl_var_1938, lay_var_1884)
          fac011_var_1913 = fk0_var_1918 * fac11_var_1849(jl_var_1938, lay_var_1884)
          fac111_var_1914 = fk1_var_1919 * fac11_var_1849(jl_var_1938, lay_var_1884)
          fac211_var_1915 = fk2_var_1920 * fac11_var_1849(jl_var_1938, lay_var_1884)
        ELSE IF (specparm1_var_1893 .GT. 0.875D0) THEN
          p_var_1916 = - fs1_var_1891
          p4_var_1917 = p_var_1916 ** 4
          fk0_var_1918 = p4_var_1917
          fk1_var_1919 = 1.0D0 - p_var_1916 - 2.0D0 * p4_var_1917
          fk2_var_1920 = p_var_1916 + p4_var_1917
          fac001_var_1910 = fk0_var_1918 * fac01_var_1847(jl_var_1938, lay_var_1884)
          fac101_var_1911 = fk1_var_1919 * fac01_var_1847(jl_var_1938, lay_var_1884)
          fac201_var_1912 = fk2_var_1920 * fac01_var_1847(jl_var_1938, lay_var_1884)
          fac011_var_1913 = fk0_var_1918 * fac11_var_1849(jl_var_1938, lay_var_1884)
          fac111_var_1914 = fk1_var_1919 * fac11_var_1849(jl_var_1938, lay_var_1884)
          fac211_var_1915 = fk2_var_1920 * fac11_var_1849(jl_var_1938, lay_var_1884)
        ELSE
          fac001_var_1910 = (1.0D0 - fs1_var_1891) * fac01_var_1847(jl_var_1938, lay_var_1884)
          fac011_var_1913 = (1.0D0 - fs1_var_1891) * fac11_var_1849(jl_var_1938, lay_var_1884)
          fac101_var_1911 = fs1_var_1891 * fac01_var_1847(jl_var_1938, lay_var_1884)
          fac111_var_1914 = fs1_var_1891 * fac11_var_1849(jl_var_1938, lay_var_1884)
          fac201_var_1912 = 0.0D0
          fac211_var_1915 = 0.0D0
        END IF
        IF (specparm_var_1890 .LT. 0.125D0) THEN
          tau_major_var_1926(1 : ng3) = speccomb_var_1870 * (fac000_var_1904 * absa_var_163(ind0_var_1877, 1 : 16) + fac100_var_1905 * absa_var_163(ind0_var_1877 + 1, 1 : 16) + fac200_var_1906 * absa_var_163(ind0_var_1877 + 2, 1 : 16) + fac010_var_1907 * absa_var_163(ind0_var_1877 + 9, 1 : 16) + fac110_var_1908 * absa_var_163(ind0_var_1877 + 10, 1 : 16) + fac210_var_1909 * absa_var_163(ind0_var_1877 + 11, 1 : 16))
        ELSE IF (specparm_var_1890 .GT. 0.875D0) THEN
          tau_major_var_1926(1 : ng3) = speccomb_var_1870 * (fac200_var_1906 * absa_var_163(ind0_var_1877 - 1, 1 : 16) + fac100_var_1905 * absa_var_163(ind0_var_1877, 1 : 16) + fac000_var_1904 * absa_var_163(ind0_var_1877 + 1, 1 : 16) + fac210_var_1909 * absa_var_163(ind0_var_1877 + 8, 1 : 16) + fac110_var_1908 * absa_var_163(ind0_var_1877 + 9, 1 : 16) + fac010_var_1907 * absa_var_163(ind0_var_1877 + 10, 1 : 16))
        ELSE
          tau_major_var_1926(1 : ng3) = speccomb_var_1870 * (fac000_var_1904 * absa_var_163(ind0_var_1877, 1 : 16) + fac100_var_1905 * absa_var_163(ind0_var_1877 + 1, 1 : 16) + fac010_var_1907 * absa_var_163(ind0_var_1877 + 9, 1 : 16) + fac110_var_1908 * absa_var_163(ind0_var_1877 + 10, 1 : 16))
        END IF
        IF (specparm1_var_1893 .LT. 0.125D0) THEN
          tau_major1_var_1927(1 : ng3) = speccomb1_var_1871 * (fac001_var_1910 * absa_var_163(ind1_var_1878, 1 : 16) + fac101_var_1911 * absa_var_163(ind1_var_1878 + 1, 1 : 16) + fac201_var_1912 * absa_var_163(ind1_var_1878 + 2, 1 : 16) + fac011_var_1913 * absa_var_163(ind1_var_1878 + 9, 1 : 16) + fac111_var_1914 * absa_var_163(ind1_var_1878 + 10, 1 : 16) + fac211_var_1915 * absa_var_163(ind1_var_1878 + 11, 1 : 16))
        ELSE IF (specparm1_var_1893 .GT. 0.875D0) THEN
          tau_major1_var_1927(1 : ng3) = speccomb1_var_1871 * (fac201_var_1912 * absa_var_163(ind1_var_1878 - 1, 1 : 16) + fac101_var_1911 * absa_var_163(ind1_var_1878, 1 : 16) + fac001_var_1910 * absa_var_163(ind1_var_1878 + 1, 1 : 16) + fac211_var_1915 * absa_var_163(ind1_var_1878 + 8, 1 : 16) + fac111_var_1914 * absa_var_163(ind1_var_1878 + 9, 1 : 16) + fac011_var_1913 * absa_var_163(ind1_var_1878 + 10, 1 : 16))
        ELSE
          tau_major1_var_1927(1 : ng3) = speccomb1_var_1871 * (fac001_var_1910 * absa_var_163(ind1_var_1878, 1 : 16) + fac101_var_1911 * absa_var_163(ind1_var_1878 + 1, 1 : 16) + fac011_var_1913 * absa_var_163(ind1_var_1878 + 9, 1 : 16) + fac111_var_1914 * absa_var_163(ind1_var_1878 + 10, 1 : 16))
        END IF
        DO ig_var_1882 = 1, 16
          tauself_var_1922 = selffac_var_1860(jl_var_1938, lay_var_1884) * (selfref_var_165(inds_var_1879, ig_var_1882) + selffrac_var_1861(jl_var_1938, lay_var_1884) * (selfref_var_165(inds_var_1879 + 1, ig_var_1882) - selfref_var_165(inds_var_1879, ig_var_1882)))
          taufor_var_1921 = forfac_var_1850(jl_var_1938, lay_var_1884) * (forref_var_166(indf_var_1880, ig_var_1882) + forfrac_var_1867(jl_var_1938, lay_var_1884) * (forref_var_166(indf_var_1880 + 1, ig_var_1882) - forref_var_166(indf_var_1880, ig_var_1882)))
          n2om1_var_1923 = ka_mn2o_var_161(jmn2o_var_1886, indm_var_1881, ig_var_1882) + fmn2o_var_1894 * (ka_mn2o_var_161(jmn2o_var_1886 + 1, indm_var_1881, ig_var_1882) - ka_mn2o_var_161(jmn2o_var_1886, indm_var_1881, ig_var_1882))
          n2om2_var_1924 = ka_mn2o_var_161(jmn2o_var_1886, indm_var_1881 + 1, ig_var_1882) + fmn2o_var_1894 * (ka_mn2o_var_161(jmn2o_var_1886 + 1, indm_var_1881 + 1, ig_var_1882) - ka_mn2o_var_161(jmn2o_var_1886, indm_var_1881 + 1, ig_var_1882))
          absn2o_var_1925 = n2om1_var_1923 + minorfrac_var_1868(jl_var_1938, lay_var_1884) * (n2om2_var_1924 - n2om1_var_1923)
          taug_var_1844(jl_var_1938, 22 + ig_var_1882, lay_var_1884) = tau_major_var_1926(ig_var_1882) + tau_major1_var_1927(ig_var_1882) + tauself_var_1922 + taufor_var_1921 + adjcoln2o_var_1901 * absn2o_var_1925
          fracs_var_1863(jl_var_1938, 22 + ig_var_1882, lay_var_1884) = fracrefa_var_159(ig_var_1882, jpl_var_1887) + fpl_var_1897 * (fracrefa_var_159(ig_var_1882, jpl_var_1887 + 1) - fracrefa_var_159(ig_var_1882, jpl_var_1887))
        END DO
      END DO
      ixc0_var_1935 = kfdia_var_1842 - kidia_var_1841 + 1 - ixc0_var_1935
      DO ixp_var_1936 = 1, ixc0_var_1935
        jl_var_1938 = ixhigh_var_1932(ixp_var_1936, lay_var_1884)
        speccomb_var_1870 = colh2o_var_1855(jl_var_1938, lay_var_1884) + rat_h2oco2_var_1864(jl_var_1938, lay_var_1884) * colco2_var_1856(jl_var_1938, lay_var_1884)
        specparm_var_1890 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_var_1870, oneminus_var_1854)
        specmult_var_1889 = 4.0D0 * (specparm_var_1890)
        js_var_1883 = 1 + INT(specmult_var_1889)
        fs_var_1888 = ((specmult_var_1889) - AINT((specmult_var_1889)))
        speccomb1_var_1871 = colh2o_var_1855(jl_var_1938, lay_var_1884) + rat_h2oco2_1_var_1865(jl_var_1938, lay_var_1884) * colco2_var_1856(jl_var_1938, lay_var_1884)
        specparm1_var_1893 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb1_var_1871, oneminus_var_1854)
        specmult1_var_1892 = 4.0D0 * (specparm1_var_1893)
        js1_var_1885 = 1 + INT(specmult1_var_1892)
        fs1_var_1891 = ((specmult1_var_1892) - AINT((specmult1_var_1892)))
        fac000_var_1904 = (1.0D0 - fs_var_1888) * fac00_var_1846(jl_var_1938, lay_var_1884)
        fac010_var_1907 = (1.0D0 - fs_var_1888) * fac10_var_1848(jl_var_1938, lay_var_1884)
        fac100_var_1905 = fs_var_1888 * fac00_var_1846(jl_var_1938, lay_var_1884)
        fac110_var_1908 = fs_var_1888 * fac10_var_1848(jl_var_1938, lay_var_1884)
        fac001_var_1910 = (1.0D0 - fs1_var_1891) * fac01_var_1847(jl_var_1938, lay_var_1884)
        fac011_var_1913 = (1.0D0 - fs1_var_1891) * fac11_var_1849(jl_var_1938, lay_var_1884)
        fac101_var_1911 = fs1_var_1891 * fac01_var_1847(jl_var_1938, lay_var_1884)
        fac111_var_1914 = fs1_var_1891 * fac11_var_1849(jl_var_1938, lay_var_1884)
        speccomb_mn2o_var_1872 = colh2o_var_1855(jl_var_1938, lay_var_1884) + refrat_m_b * colco2_var_1856(jl_var_1938, lay_var_1884)
        specparm_mn2o_var_1896 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_mn2o_var_1872, oneminus_var_1854)
        specmult_mn2o_var_1895 = 4.0D0 * specparm_mn2o_var_1896
        jmn2o_var_1886 = 1 + INT(specmult_mn2o_var_1895)
        fmn2o_var_1894 = ((specmult_mn2o_var_1895) - AINT((specmult_mn2o_var_1895)))
        chi_n2o_var_1903 = coln2o_var_1857(jl_var_1938, lay_var_1884) / coldry_var_1858(jl_var_1938, lay_var_1884)
        ratn2o_var_1902 = 1D+20 * chi_n2o_var_1903 / chi_mls(4, jp_var_1851(jl_var_1938, lay_var_1884) + 1)
        IF (ratn2o_var_1902 .GT. 1.5D0) THEN
          adjfac_var_1900 = 0.5D0 + (ratn2o_var_1902 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1901 = adjfac_var_1900 * chi_mls(4, jp_var_1851(jl_var_1938, lay_var_1884) + 1) * coldry_var_1858(jl_var_1938, lay_var_1884) * 1D-20
        ELSE
          adjcoln2o_var_1901 = coln2o_var_1857(jl_var_1938, lay_var_1884)
        END IF
        speccomb_planck_var_1873 = colh2o_var_1855(jl_var_1938, lay_var_1884) + refrat_planck_b_var_1875 * colco2_var_1856(jl_var_1938, lay_var_1884)
        specparm_planck_var_1899 = MIN(colh2o_var_1855(jl_var_1938, lay_var_1884) / speccomb_planck_var_1873, oneminus_var_1854)
        specmult_planck_var_1898 = 4.0D0 * specparm_planck_var_1899
        jpl_var_1887 = 1 + INT(specmult_planck_var_1898)
        fpl_var_1897 = ((specmult_planck_var_1898) - AINT((specmult_planck_var_1898)))
        ind0_var_1877 = ((jp_var_1851(jl_var_1938, lay_var_1884) - 13) * 5 + (jt_var_1852(jl_var_1938, lay_var_1884) - 1)) * nspb_var_217(3) + js_var_1883
        ind1_var_1878 = ((jp_var_1851(jl_var_1938, lay_var_1884) - 12) * 5 + (jt1_var_1853(jl_var_1938, lay_var_1884) - 1)) * nspb_var_217(3) + js1_var_1885
        indf_var_1880 = indfor_var_1866(jl_var_1938, lay_var_1884)
        indm_var_1881 = indminor_var_1869(jl_var_1938, lay_var_1884)
        DO ig_var_1882 = 1, 16
          taufor_var_1921 = forfac_var_1850(jl_var_1938, lay_var_1884) * (forref_var_166(indf_var_1880, ig_var_1882) + forfrac_var_1867(jl_var_1938, lay_var_1884) * (forref_var_166(indf_var_1880 + 1, ig_var_1882) - forref_var_166(indf_var_1880, ig_var_1882)))
          n2om1_var_1923 = kb_mn2o_var_162(jmn2o_var_1886, indm_var_1881, ig_var_1882) + fmn2o_var_1894 * (kb_mn2o_var_162(jmn2o_var_1886 + 1, indm_var_1881, ig_var_1882) - kb_mn2o_var_162(jmn2o_var_1886, indm_var_1881, ig_var_1882))
          n2om2_var_1924 = kb_mn2o_var_162(jmn2o_var_1886, indm_var_1881 + 1, ig_var_1882) + fmn2o_var_1894 * (kb_mn2o_var_162(jmn2o_var_1886 + 1, indm_var_1881 + 1, ig_var_1882) - kb_mn2o_var_162(jmn2o_var_1886, indm_var_1881 + 1, ig_var_1882))
          absn2o_var_1925 = n2om1_var_1923 + minorfrac_var_1868(jl_var_1938, lay_var_1884) * (n2om2_var_1924 - n2om1_var_1923)
          taug_var_1844(jl_var_1938, 22 + ig_var_1882, lay_var_1884) = speccomb_var_1870 * (fac000_var_1904 * absb_var_164(ind0_var_1877, ig_var_1882) + fac100_var_1905 * absb_var_164(ind0_var_1877 + 1, ig_var_1882) + fac010_var_1907 * absb_var_164(ind0_var_1877 + 5, ig_var_1882) + fac110_var_1908 * absb_var_164(ind0_var_1877 + 6, ig_var_1882)) + speccomb1_var_1871 * (fac001_var_1910 * absb_var_164(ind1_var_1878, ig_var_1882) + fac101_var_1911 * absb_var_164(ind1_var_1878 + 1, ig_var_1882) + fac011_var_1913 * absb_var_164(ind1_var_1878 + 5, ig_var_1882) + fac111_var_1914 * absb_var_164(ind1_var_1878 + 6, ig_var_1882)) + taufor_var_1921 + adjcoln2o_var_1901 * absn2o_var_1925
          fracs_var_1863(jl_var_1938, 22 + ig_var_1882, lay_var_1884) = fracrefb_var_160(ig_var_1882, jpl_var_1887) + fpl_var_1897 * (fracrefb_var_160(ig_var_1882, jpl_var_1887 + 1) - fracrefb_var_160(ig_var_1882, jpl_var_1887))
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol3
SUBROUTINE rrtm_taumol2(kidia_var_1939, kfdia_var_1940, klev_var_1941, taug_var_1942, pavel_var_1943, coldry_var_1944, p_tauaerl_var_1945, fac00_var_1946, fac01_var_1947, fac10_var_1948, fac11_var_1949, forfac_var_1951, forfrac_var_1950, indfor_var_1960, jp_var_1952, jt_var_1953, jt1_var_1954, colh2o_var_1955, laytrop_var_1956, selffac_var_1957, selffrac_var_1958, indself_var_1959, fracs_var_1961)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta2, ONLY: absa_var_155, absb_var_156, forref_var_158, fracrefa_var_153, fracrefb_var_154, selfref_var_157
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1939
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1940
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1941
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1942(kidia_var_1939 : kfdia_var_1940, 140, klev_var_1941)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1943(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1944(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1945(kidia_var_1939 : kfdia_var_1940, klev_var_1941, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1946(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1947(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1948(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1949(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1950(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1951(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1952(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1953(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1954(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1955(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1956(kidia_var_1939 : kfdia_var_1940)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1957(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1958(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1959(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1960(kidia_var_1939 : kfdia_var_1940, klev_var_1941)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1961(kidia_var_1939 : kfdia_var_1940, 140, klev_var_1941)
  INTEGER(KIND = 4) :: ind0_var_1962, ind1_var_1963, inds_var_1964, indf_var_1965
  INTEGER(KIND = 4) :: ig_var_1966, lay_var_1967
  REAL(KIND = 8) :: taufor_var_1968, tauself_var_1969, corradj_var_1970, pp_var_1971
  INTEGER(KIND = 4) :: laytrop_min_var_1972, laytrop_max_var_1973
  INTEGER(KIND = 4) :: ixc_var_1974(klev_var_1941), ixlow_var_1975(kfdia_var_1940, klev_var_1941), ixhigh_var_1976(kfdia_var_1940, klev_var_1941)
  INTEGER(KIND = 4) :: ich_var_1977, icl_var_1978, ixc0_var_1979, ixp_var_1980, jc_var_1981, jl_var_1982
  laytrop_min_var_1972 = MINVAL(laytrop_var_1956)
  laytrop_max_var_1973 = MAXVAL(laytrop_var_1956)
  ixlow_var_1975 = 0
  ixhigh_var_1976 = 0
  ixc_var_1974 = 0
  DO lay_var_1967 = laytrop_min_var_1972 + 1, laytrop_max_var_1973
    icl_var_1978 = 0
    ich_var_1977 = 0
    DO jc_var_1981 = kidia_var_1939, kfdia_var_1940
      IF (lay_var_1967 <= laytrop_var_1956(jc_var_1981)) THEN
        icl_var_1978 = icl_var_1978 + 1
        ixlow_var_1975(icl_var_1978, lay_var_1967) = jc_var_1981
      ELSE
        ich_var_1977 = ich_var_1977 + 1
        ixhigh_var_1976(ich_var_1977, lay_var_1967) = jc_var_1981
      END IF
    END DO
    ixc_var_1974(lay_var_1967) = icl_var_1978
  END DO
  DO lay_var_1967 = 1, laytrop_min_var_1972
    DO jl_var_1982 = kidia_var_1939, kfdia_var_1940
      ind0_var_1962 = ((jp_var_1952(jl_var_1982, lay_var_1967) - 1) * 5 + (jt_var_1953(jl_var_1982, lay_var_1967) - 1)) * nspa_var_216(2) + 1
      ind1_var_1963 = (jp_var_1952(jl_var_1982, lay_var_1967) * 5 + (jt1_var_1954(jl_var_1982, lay_var_1967) - 1)) * nspa_var_216(2) + 1
      inds_var_1964 = indself_var_1959(jl_var_1982, lay_var_1967)
      indf_var_1965 = indfor_var_1960(jl_var_1982, lay_var_1967)
      pp_var_1971 = pavel_var_1943(jl_var_1982, lay_var_1967)
      corradj_var_1970 = 1.0D0 - 0.05D0 * (pp_var_1971 - 100.0D0) / 900.0D0
      DO ig_var_1966 = 1, 12
        tauself_var_1969 = selffac_var_1957(jl_var_1982, lay_var_1967) * (selfref_var_157(inds_var_1964, ig_var_1966) + selffrac_var_1958(jl_var_1982, lay_var_1967) * (selfref_var_157(inds_var_1964 + 1, ig_var_1966) - selfref_var_157(inds_var_1964, ig_var_1966)))
        taufor_var_1968 = forfac_var_1951(jl_var_1982, lay_var_1967) * (forref_var_158(indf_var_1965, ig_var_1966) + forfrac_var_1950(jl_var_1982, lay_var_1967) * (forref_var_158(indf_var_1965 + 1, ig_var_1966) - forref_var_158(indf_var_1965, ig_var_1966)))
        taug_var_1942(jl_var_1982, 10 + ig_var_1966, lay_var_1967) = corradj_var_1970 * (colh2o_var_1955(jl_var_1982, lay_var_1967) * (fac00_var_1946(jl_var_1982, lay_var_1967) * absa_var_155(ind0_var_1962, ig_var_1966) + fac10_var_1948(jl_var_1982, lay_var_1967) * absa_var_155(ind0_var_1962 + 1, ig_var_1966) + fac01_var_1947(jl_var_1982, lay_var_1967) * absa_var_155(ind1_var_1963, ig_var_1966) + fac11_var_1949(jl_var_1982, lay_var_1967) * absa_var_155(ind1_var_1963 + 1, ig_var_1966)) + tauself_var_1969 + taufor_var_1968)
        fracs_var_1961(jl_var_1982, 10 + ig_var_1966, lay_var_1967) = fracrefa_var_153(ig_var_1966)
      END DO
    END DO
  END DO
  DO lay_var_1967 = laytrop_max_var_1973 + 1, klev_var_1941
    DO jl_var_1982 = kidia_var_1939, kfdia_var_1940
      ind0_var_1962 = ((jp_var_1952(jl_var_1982, lay_var_1967) - 13) * 5 + (jt_var_1953(jl_var_1982, lay_var_1967) - 1)) * nspb_var_217(2) + 1
      ind1_var_1963 = ((jp_var_1952(jl_var_1982, lay_var_1967) - 12) * 5 + (jt1_var_1954(jl_var_1982, lay_var_1967) - 1)) * nspb_var_217(2) + 1
      indf_var_1965 = indfor_var_1960(jl_var_1982, lay_var_1967)
      DO ig_var_1966 = 1, 12
        taufor_var_1968 = forfac_var_1951(jl_var_1982, lay_var_1967) * (forref_var_158(indf_var_1965, ig_var_1966) + forfrac_var_1950(jl_var_1982, lay_var_1967) * (forref_var_158(indf_var_1965 + 1, ig_var_1966) - forref_var_158(indf_var_1965, ig_var_1966)))
        taug_var_1942(jl_var_1982, 10 + ig_var_1966, lay_var_1967) = colh2o_var_1955(jl_var_1982, lay_var_1967) * (fac00_var_1946(jl_var_1982, lay_var_1967) * absb_var_156(ind0_var_1962, ig_var_1966) + fac10_var_1948(jl_var_1982, lay_var_1967) * absb_var_156(ind0_var_1962 + 1, ig_var_1966) + fac01_var_1947(jl_var_1982, lay_var_1967) * absb_var_156(ind1_var_1963, ig_var_1966) + fac11_var_1949(jl_var_1982, lay_var_1967) * absb_var_156(ind1_var_1963 + 1, ig_var_1966)) + taufor_var_1968
        fracs_var_1961(jl_var_1982, 10 + ig_var_1966, lay_var_1967) = fracrefb_var_154(ig_var_1966)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1973 /= laytrop_min_var_1972) THEN
    DO lay_var_1967 = laytrop_min_var_1972 + 1, laytrop_max_var_1973
      ixc0_var_1979 = ixc_var_1974(lay_var_1967)
      DO ixp_var_1980 = 1, ixc0_var_1979
        jl_var_1982 = ixlow_var_1975(ixp_var_1980, lay_var_1967)
        ind0_var_1962 = ((jp_var_1952(jl_var_1982, lay_var_1967) - 1) * 5 + (jt_var_1953(jl_var_1982, lay_var_1967) - 1)) * nspa_var_216(2) + 1
        ind1_var_1963 = (jp_var_1952(jl_var_1982, lay_var_1967) * 5 + (jt1_var_1954(jl_var_1982, lay_var_1967) - 1)) * nspa_var_216(2) + 1
        inds_var_1964 = indself_var_1959(jl_var_1982, lay_var_1967)
        indf_var_1965 = indfor_var_1960(jl_var_1982, lay_var_1967)
        pp_var_1971 = pavel_var_1943(jl_var_1982, lay_var_1967)
        corradj_var_1970 = 1.0D0 - 0.05D0 * (pp_var_1971 - 100.0D0) / 900.0D0
        DO ig_var_1966 = 1, 12
          tauself_var_1969 = selffac_var_1957(jl_var_1982, lay_var_1967) * (selfref_var_157(inds_var_1964, ig_var_1966) + selffrac_var_1958(jl_var_1982, lay_var_1967) * (selfref_var_157(inds_var_1964 + 1, ig_var_1966) - selfref_var_157(inds_var_1964, ig_var_1966)))
          taufor_var_1968 = forfac_var_1951(jl_var_1982, lay_var_1967) * (forref_var_158(indf_var_1965, ig_var_1966) + forfrac_var_1950(jl_var_1982, lay_var_1967) * (forref_var_158(indf_var_1965 + 1, ig_var_1966) - forref_var_158(indf_var_1965, ig_var_1966)))
          taug_var_1942(jl_var_1982, 10 + ig_var_1966, lay_var_1967) = corradj_var_1970 * (colh2o_var_1955(jl_var_1982, lay_var_1967) * (fac00_var_1946(jl_var_1982, lay_var_1967) * absa_var_155(ind0_var_1962, ig_var_1966) + fac10_var_1948(jl_var_1982, lay_var_1967) * absa_var_155(ind0_var_1962 + 1, ig_var_1966) + fac01_var_1947(jl_var_1982, lay_var_1967) * absa_var_155(ind1_var_1963, ig_var_1966) + fac11_var_1949(jl_var_1982, lay_var_1967) * absa_var_155(ind1_var_1963 + 1, ig_var_1966)) + tauself_var_1969 + taufor_var_1968)
          fracs_var_1961(jl_var_1982, 10 + ig_var_1966, lay_var_1967) = fracrefa_var_153(ig_var_1966)
        END DO
      END DO
      ixc0_var_1979 = kfdia_var_1940 - kidia_var_1939 + 1 - ixc0_var_1979
      DO ixp_var_1980 = 1, ixc0_var_1979
        jl_var_1982 = ixhigh_var_1976(ixp_var_1980, lay_var_1967)
        ind0_var_1962 = ((jp_var_1952(jl_var_1982, lay_var_1967) - 13) * 5 + (jt_var_1953(jl_var_1982, lay_var_1967) - 1)) * nspb_var_217(2) + 1
        ind1_var_1963 = ((jp_var_1952(jl_var_1982, lay_var_1967) - 12) * 5 + (jt1_var_1954(jl_var_1982, lay_var_1967) - 1)) * nspb_var_217(2) + 1
        indf_var_1965 = indfor_var_1960(jl_var_1982, lay_var_1967)
        DO ig_var_1966 = 1, 12
          taufor_var_1968 = forfac_var_1951(jl_var_1982, lay_var_1967) * (forref_var_158(indf_var_1965, ig_var_1966) + forfrac_var_1950(jl_var_1982, lay_var_1967) * (forref_var_158(indf_var_1965 + 1, ig_var_1966) - forref_var_158(indf_var_1965, ig_var_1966)))
          taug_var_1942(jl_var_1982, 10 + ig_var_1966, lay_var_1967) = colh2o_var_1955(jl_var_1982, lay_var_1967) * (fac00_var_1946(jl_var_1982, lay_var_1967) * absb_var_156(ind0_var_1962, ig_var_1966) + fac10_var_1948(jl_var_1982, lay_var_1967) * absb_var_156(ind0_var_1962 + 1, ig_var_1966) + fac01_var_1947(jl_var_1982, lay_var_1967) * absb_var_156(ind1_var_1963, ig_var_1966) + fac11_var_1949(jl_var_1982, lay_var_1967) * absb_var_156(ind1_var_1963 + 1, ig_var_1966)) + taufor_var_1968
          fracs_var_1961(jl_var_1982, 10 + ig_var_1966, lay_var_1967) = fracrefb_var_154(ig_var_1966)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol2
SUBROUTINE rrtm_taumol1(kidia_var_1983, kfdia_var_1984, klev_var_1985, taug_var_1987, pavel_var_1986, p_tauaerl_var_1988, fac00_var_1989, fac01_var_1990, fac10_var_1991, fac11_var_1992, forfac_var_1993, forfrac_var_1994, indfor_var_2005, jp_var_1995, jt_var_1996, jt1_var_1997, colh2o_var_1998, laytrop_var_1999, selffac_var_2000, selffrac_var_2001, indself_var_2003, fracs_var_2004, minorfrac_var_2002, indminor_var_2006, scaleminorn2, colbrd_var_2007)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta1, ONLY: absa_var_109, absb_var_110, forref_var_113, fracrefa_var_107, fracrefb_var_108, ka_mn2_var_111, kb_mn2, selfref_var_112
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1983
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1984
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1985
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1986(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1987(kidia_var_1983 : kfdia_var_1984, 140, klev_var_1985)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1988(kidia_var_1983 : kfdia_var_1984, klev_var_1985, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1989(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1990(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1991(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1992(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1993(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1994(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1995(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1996(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1997(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1998(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1999(kidia_var_1983 : kfdia_var_1984)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2000(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2001(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2002(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2003(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2004(kidia_var_1983 : kfdia_var_1984, 140, klev_var_1985)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2005(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2006(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(IN) :: scaleminorn2(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  REAL(KIND = 8), INTENT(IN) :: colbrd_var_2007(kidia_var_1983 : kfdia_var_1984, klev_var_1985)
  INTEGER(KIND = 4) :: ind0_var_2008, ind1_var_2009, inds_var_2010
  INTEGER(KIND = 4) :: indf_var_2011, indm_var_2012
  INTEGER(KIND = 4) :: ig_var_2013, lay_var_2014
  REAL(KIND = 8) :: taufor_var_2015, tauself_var_2016, corradj_var_2017, pp_var_2018, scalen2_var_2019, taun2_var_2020
  INTEGER(KIND = 4) :: laytrop_min_var_2021, laytrop_max_var_2022
  INTEGER(KIND = 4) :: ixc_var_2023(klev_var_1985), ixlow_var_2024(kfdia_var_1984, klev_var_1985), ixhigh_var_2025(kfdia_var_1984, klev_var_1985)
  INTEGER(KIND = 4) :: ich_var_2026, icl_var_2027, ixc0_var_2028, ixp_var_2029, jc_var_2030, jl_var_2031
  laytrop_min_var_2021 = MINVAL(laytrop_var_1999)
  laytrop_max_var_2022 = MAXVAL(laytrop_var_1999)
  ixlow_var_2024 = 0
  ixhigh_var_2025 = 0
  ixc_var_2023 = 0
  DO lay_var_2014 = laytrop_min_var_2021 + 1, laytrop_max_var_2022
    icl_var_2027 = 0
    ich_var_2026 = 0
    DO jc_var_2030 = kidia_var_1983, kfdia_var_1984
      IF (lay_var_2014 <= laytrop_var_1999(jc_var_2030)) THEN
        icl_var_2027 = icl_var_2027 + 1
        ixlow_var_2024(icl_var_2027, lay_var_2014) = jc_var_2030
      ELSE
        ich_var_2026 = ich_var_2026 + 1
        ixhigh_var_2025(ich_var_2026, lay_var_2014) = jc_var_2030
      END IF
    END DO
    ixc_var_2023(lay_var_2014) = icl_var_2027
  END DO
  DO lay_var_2014 = 1, laytrop_min_var_2021
    DO jl_var_2031 = kidia_var_1983, kfdia_var_1984
      ind0_var_2008 = ((jp_var_1995(jl_var_2031, lay_var_2014) - 1) * 5 + (jt_var_1996(jl_var_2031, lay_var_2014) - 1)) * nspa_var_216(1) + 1
      ind1_var_2009 = (jp_var_1995(jl_var_2031, lay_var_2014) * 5 + (jt1_var_1997(jl_var_2031, lay_var_2014) - 1)) * nspa_var_216(1) + 1
      inds_var_2010 = indself_var_2003(jl_var_2031, lay_var_2014)
      indf_var_2011 = indfor_var_2005(jl_var_2031, lay_var_2014)
      indm_var_2012 = indminor_var_2006(jl_var_2031, lay_var_2014)
      pp_var_2018 = pavel_var_1986(jl_var_2031, lay_var_2014)
      corradj_var_2017 = 1.0D0
      IF (pp_var_2018 .LT. 250.0D0) THEN
        corradj_var_2017 = 1.0D0 - 0.15D0 * (250.0D0 - pp_var_2018) / 154.4D0
      END IF
      scalen2_var_2019 = colbrd_var_2007(jl_var_2031, lay_var_2014) * scaleminorn2(jl_var_2031, lay_var_2014)
      DO ig_var_2013 = 1, 10
        tauself_var_2016 = selffac_var_2000(jl_var_2031, lay_var_2014) * (selfref_var_112(inds_var_2010, ig_var_2013) + selffrac_var_2001(jl_var_2031, lay_var_2014) * (selfref_var_112(inds_var_2010 + 1, ig_var_2013) - selfref_var_112(inds_var_2010, ig_var_2013)))
        taufor_var_2015 = forfac_var_1993(jl_var_2031, lay_var_2014) * (forref_var_113(indf_var_2011, ig_var_2013) + forfrac_var_1994(jl_var_2031, lay_var_2014) * (forref_var_113(indf_var_2011 + 1, ig_var_2013) - forref_var_113(indf_var_2011, ig_var_2013)))
        taun2_var_2020 = scalen2_var_2019 * (ka_mn2_var_111(indm_var_2012, ig_var_2013) + minorfrac_var_2002(jl_var_2031, lay_var_2014) * (ka_mn2_var_111(indm_var_2012 + 1, ig_var_2013) - ka_mn2_var_111(indm_var_2012, ig_var_2013)))
        taug_var_1987(jl_var_2031, ig_var_2013, lay_var_2014) = corradj_var_2017 * (colh2o_var_1998(jl_var_2031, lay_var_2014) * (fac00_var_1989(jl_var_2031, lay_var_2014) * absa_var_109(ind0_var_2008, ig_var_2013) + fac10_var_1991(jl_var_2031, lay_var_2014) * absa_var_109(ind0_var_2008 + 1, ig_var_2013) + fac01_var_1990(jl_var_2031, lay_var_2014) * absa_var_109(ind1_var_2009, ig_var_2013) + fac11_var_1992(jl_var_2031, lay_var_2014) * absa_var_109(ind1_var_2009 + 1, ig_var_2013)) + tauself_var_2016 + taufor_var_2015 + taun2_var_2020)
        fracs_var_2004(jl_var_2031, ig_var_2013, lay_var_2014) = fracrefa_var_107(ig_var_2013)
      END DO
    END DO
  END DO
  DO lay_var_2014 = laytrop_max_var_2022 + 1, klev_var_1985
    DO jl_var_2031 = kidia_var_1983, kfdia_var_1984
      ind0_var_2008 = ((jp_var_1995(jl_var_2031, lay_var_2014) - 13) * 5 + (jt_var_1996(jl_var_2031, lay_var_2014) - 1)) * nspb_var_217(1) + 1
      ind1_var_2009 = ((jp_var_1995(jl_var_2031, lay_var_2014) - 12) * 5 + (jt1_var_1997(jl_var_2031, lay_var_2014) - 1)) * nspb_var_217(1) + 1
      indf_var_2011 = indfor_var_2005(jl_var_2031, lay_var_2014)
      indm_var_2012 = indminor_var_2006(jl_var_2031, lay_var_2014)
      pp_var_2018 = pavel_var_1986(jl_var_2031, lay_var_2014)
      corradj_var_2017 = 1.0D0 - 0.15D0 * (pp_var_2018 / 95.6D0)
      scalen2_var_2019 = colbrd_var_2007(jl_var_2031, lay_var_2014) * scaleminorn2(jl_var_2031, lay_var_2014)
      DO ig_var_2013 = 1, 10
        taufor_var_2015 = forfac_var_1993(jl_var_2031, lay_var_2014) * (forref_var_113(indf_var_2011, ig_var_2013) + forfrac_var_1994(jl_var_2031, lay_var_2014) * (forref_var_113(indf_var_2011 + 1, ig_var_2013) - forref_var_113(indf_var_2011, ig_var_2013)))
        taun2_var_2020 = scalen2_var_2019 * (kb_mn2(indm_var_2012, ig_var_2013) + minorfrac_var_2002(jl_var_2031, lay_var_2014) * (kb_mn2(indm_var_2012 + 1, ig_var_2013) - kb_mn2(indm_var_2012, ig_var_2013)))
        taug_var_1987(jl_var_2031, ig_var_2013, lay_var_2014) = corradj_var_2017 * (colh2o_var_1998(jl_var_2031, lay_var_2014) * (fac00_var_1989(jl_var_2031, lay_var_2014) * absb_var_110(ind0_var_2008, ig_var_2013) + fac10_var_1991(jl_var_2031, lay_var_2014) * absb_var_110(ind0_var_2008 + 1, ig_var_2013) + fac01_var_1990(jl_var_2031, lay_var_2014) * absb_var_110(ind1_var_2009, ig_var_2013) + fac11_var_1992(jl_var_2031, lay_var_2014) * absb_var_110(ind1_var_2009 + 1, ig_var_2013)) + taufor_var_2015 + taun2_var_2020)
        fracs_var_2004(jl_var_2031, ig_var_2013, lay_var_2014) = fracrefb_var_108(ig_var_2013)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2022 /= laytrop_min_var_2021) THEN
    DO lay_var_2014 = laytrop_min_var_2021 + 1, laytrop_max_var_2022
      ixc0_var_2028 = ixc_var_2023(lay_var_2014)
      DO ixp_var_2029 = 1, ixc0_var_2028
        jl_var_2031 = ixlow_var_2024(ixp_var_2029, lay_var_2014)
        ind0_var_2008 = ((jp_var_1995(jl_var_2031, lay_var_2014) - 1) * 5 + (jt_var_1996(jl_var_2031, lay_var_2014) - 1)) * nspa_var_216(1) + 1
        ind1_var_2009 = (jp_var_1995(jl_var_2031, lay_var_2014) * 5 + (jt1_var_1997(jl_var_2031, lay_var_2014) - 1)) * nspa_var_216(1) + 1
        inds_var_2010 = indself_var_2003(jl_var_2031, lay_var_2014)
        indf_var_2011 = indfor_var_2005(jl_var_2031, lay_var_2014)
        indm_var_2012 = indminor_var_2006(jl_var_2031, lay_var_2014)
        pp_var_2018 = pavel_var_1986(jl_var_2031, lay_var_2014)
        corradj_var_2017 = 1.0D0
        IF (pp_var_2018 .LT. 250.0D0) THEN
          corradj_var_2017 = 1.0D0 - 0.15D0 * (250.0D0 - pp_var_2018) / 154.4D0
        END IF
        scalen2_var_2019 = colbrd_var_2007(jl_var_2031, lay_var_2014) * scaleminorn2(jl_var_2031, lay_var_2014)
        DO ig_var_2013 = 1, 10
          tauself_var_2016 = selffac_var_2000(jl_var_2031, lay_var_2014) * (selfref_var_112(inds_var_2010, ig_var_2013) + selffrac_var_2001(jl_var_2031, lay_var_2014) * (selfref_var_112(inds_var_2010 + 1, ig_var_2013) - selfref_var_112(inds_var_2010, ig_var_2013)))
          taufor_var_2015 = forfac_var_1993(jl_var_2031, lay_var_2014) * (forref_var_113(indf_var_2011, ig_var_2013) + forfrac_var_1994(jl_var_2031, lay_var_2014) * (forref_var_113(indf_var_2011 + 1, ig_var_2013) - forref_var_113(indf_var_2011, ig_var_2013)))
          taun2_var_2020 = scalen2_var_2019 * (ka_mn2_var_111(indm_var_2012, ig_var_2013) + minorfrac_var_2002(jl_var_2031, lay_var_2014) * (ka_mn2_var_111(indm_var_2012 + 1, ig_var_2013) - ka_mn2_var_111(indm_var_2012, ig_var_2013)))
          taug_var_1987(jl_var_2031, ig_var_2013, lay_var_2014) = corradj_var_2017 * (colh2o_var_1998(jl_var_2031, lay_var_2014) * (fac00_var_1989(jl_var_2031, lay_var_2014) * absa_var_109(ind0_var_2008, ig_var_2013) + fac10_var_1991(jl_var_2031, lay_var_2014) * absa_var_109(ind0_var_2008 + 1, ig_var_2013) + fac01_var_1990(jl_var_2031, lay_var_2014) * absa_var_109(ind1_var_2009, ig_var_2013) + fac11_var_1992(jl_var_2031, lay_var_2014) * absa_var_109(ind1_var_2009 + 1, ig_var_2013)) + tauself_var_2016 + taufor_var_2015 + taun2_var_2020)
          fracs_var_2004(jl_var_2031, ig_var_2013, lay_var_2014) = fracrefa_var_107(ig_var_2013)
        END DO
      END DO
      ixc0_var_2028 = kfdia_var_1984 - kidia_var_1983 + 1 - ixc0_var_2028
      DO ixp_var_2029 = 1, ixc0_var_2028
        jl_var_2031 = ixhigh_var_2025(ixp_var_2029, lay_var_2014)
        ind0_var_2008 = ((jp_var_1995(jl_var_2031, lay_var_2014) - 13) * 5 + (jt_var_1996(jl_var_2031, lay_var_2014) - 1)) * nspb_var_217(1) + 1
        ind1_var_2009 = ((jp_var_1995(jl_var_2031, lay_var_2014) - 12) * 5 + (jt1_var_1997(jl_var_2031, lay_var_2014) - 1)) * nspb_var_217(1) + 1
        indf_var_2011 = indfor_var_2005(jl_var_2031, lay_var_2014)
        indm_var_2012 = indminor_var_2006(jl_var_2031, lay_var_2014)
        pp_var_2018 = pavel_var_1986(jl_var_2031, lay_var_2014)
        corradj_var_2017 = 1.0D0 - 0.15D0 * (pp_var_2018 / 95.6D0)
        scalen2_var_2019 = colbrd_var_2007(jl_var_2031, lay_var_2014) * scaleminorn2(jl_var_2031, lay_var_2014)
        DO ig_var_2013 = 1, 10
          taufor_var_2015 = forfac_var_1993(jl_var_2031, lay_var_2014) * (forref_var_113(indf_var_2011, ig_var_2013) + forfrac_var_1994(jl_var_2031, lay_var_2014) * (forref_var_113(indf_var_2011 + 1, ig_var_2013) - forref_var_113(indf_var_2011, ig_var_2013)))
          taun2_var_2020 = scalen2_var_2019 * (kb_mn2(indm_var_2012, ig_var_2013) + minorfrac_var_2002(jl_var_2031, lay_var_2014) * (kb_mn2(indm_var_2012 + 1, ig_var_2013) - kb_mn2(indm_var_2012, ig_var_2013)))
          taug_var_1987(jl_var_2031, ig_var_2013, lay_var_2014) = corradj_var_2017 * (colh2o_var_1998(jl_var_2031, lay_var_2014) * (fac00_var_1989(jl_var_2031, lay_var_2014) * absb_var_110(ind0_var_2008, ig_var_2013) + fac10_var_1991(jl_var_2031, lay_var_2014) * absb_var_110(ind0_var_2008 + 1, ig_var_2013) + fac01_var_1990(jl_var_2031, lay_var_2014) * absb_var_110(ind1_var_2009, ig_var_2013) + fac11_var_1992(jl_var_2031, lay_var_2014) * absb_var_110(ind1_var_2009 + 1, ig_var_2013)) + taufor_var_2015 + taun2_var_2020)
          fracs_var_2004(jl_var_2031, ig_var_2013, lay_var_2014) = fracrefb_var_108(ig_var_2013)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol1
SUBROUTINE rrtm_taumol14(kidia_var_2032, kfdia_var_2033, klev_var_2034, taug_var_2035, p_tauaerl_var_2036, fac00_var_2037, fac01_var_2038, fac10_var_2039, fac11_var_2040, forfac_var_2051, forfrac_var_2052, indfor_var_2050, jp_var_2041, jt_var_2042, jt1_var_2043, colco2_var_2044, laytrop_var_2045, selffac_var_2046, selffrac_var_2047, indself_var_2048, fracs_var_2049)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta14, ONLY: absa_var_138, absb_var_139, forref_var_141, fracrefa_var_136, fracrefb_var_137, selfref_var_140
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2032
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2033
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2034
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2035(kidia_var_2032 : kfdia_var_2033, 140, klev_var_2034)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2036(kidia_var_2032 : kfdia_var_2033, klev_var_2034, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2037(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2038(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2039(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2040(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2041(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2042(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2043(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2044(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2045(kidia_var_2032 : kfdia_var_2033)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2046(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2047(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2048(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2049(kidia_var_2032 : kfdia_var_2033, 140, klev_var_2034)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2050(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2051(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2052(kidia_var_2032 : kfdia_var_2033, klev_var_2034)
  INTEGER(KIND = 4) :: ig_var_2053, ind0_var_2054, ind1_var_2055, inds_var_2056, indf_var_2057, lay_var_2058
  REAL(KIND = 8) :: taufor_var_2059, tauself_var_2060
  INTEGER(KIND = 4) :: laytrop_min_var_2061, laytrop_max_var_2062
  INTEGER(KIND = 4) :: ixc_var_2063(klev_var_2034), ixlow_var_2064(kfdia_var_2033, klev_var_2034), ixhigh_var_2065(kfdia_var_2033, klev_var_2034)
  INTEGER(KIND = 4) :: ich_var_2066, icl_var_2067, ixc0_var_2068, ixp_var_2069, jc_var_2070, jl_var_2071
  laytrop_min_var_2061 = MINVAL(laytrop_var_2045)
  laytrop_max_var_2062 = MAXVAL(laytrop_var_2045)
  ixlow_var_2064 = 0
  ixhigh_var_2065 = 0
  ixc_var_2063 = 0
  DO lay_var_2058 = laytrop_min_var_2061 + 1, laytrop_max_var_2062
    icl_var_2067 = 0
    ich_var_2066 = 0
    DO jc_var_2070 = kidia_var_2032, kfdia_var_2033
      IF (lay_var_2058 <= laytrop_var_2045(jc_var_2070)) THEN
        icl_var_2067 = icl_var_2067 + 1
        ixlow_var_2064(icl_var_2067, lay_var_2058) = jc_var_2070
      ELSE
        ich_var_2066 = ich_var_2066 + 1
        ixhigh_var_2065(ich_var_2066, lay_var_2058) = jc_var_2070
      END IF
    END DO
    ixc_var_2063(lay_var_2058) = icl_var_2067
  END DO
  DO lay_var_2058 = 1, laytrop_min_var_2061
    DO jl_var_2071 = kidia_var_2032, kfdia_var_2033
      ind0_var_2054 = ((jp_var_2041(jl_var_2071, lay_var_2058) - 1) * 5 + (jt_var_2042(jl_var_2071, lay_var_2058) - 1)) * nspa_var_216(14) + 1
      ind1_var_2055 = (jp_var_2041(jl_var_2071, lay_var_2058) * 5 + (jt1_var_2043(jl_var_2071, lay_var_2058) - 1)) * nspa_var_216(14) + 1
      inds_var_2056 = indself_var_2048(jl_var_2071, lay_var_2058)
      indf_var_2057 = indfor_var_2050(jl_var_2071, lay_var_2058)
      DO ig_var_2053 = 1, 2
        tauself_var_2060 = selffac_var_2046(jl_var_2071, lay_var_2058) * (selfref_var_140(inds_var_2056, ig_var_2053) + selffrac_var_2047(jl_var_2071, lay_var_2058) * (selfref_var_140(inds_var_2056 + 1, ig_var_2053) - selfref_var_140(inds_var_2056, ig_var_2053)))
        taufor_var_2059 = forfac_var_2051(jl_var_2071, lay_var_2058) * (forref_var_141(indf_var_2057, ig_var_2053) + forfrac_var_2052(jl_var_2071, lay_var_2058) * (forref_var_141(indf_var_2057 + 1, ig_var_2053) - forref_var_141(indf_var_2057, ig_var_2053)))
        taug_var_2035(jl_var_2071, 134 + ig_var_2053, lay_var_2058) = colco2_var_2044(jl_var_2071, lay_var_2058) * (fac00_var_2037(jl_var_2071, lay_var_2058) * absa_var_138(ind0_var_2054, ig_var_2053) + fac10_var_2039(jl_var_2071, lay_var_2058) * absa_var_138(ind0_var_2054 + 1, ig_var_2053) + fac01_var_2038(jl_var_2071, lay_var_2058) * absa_var_138(ind1_var_2055, ig_var_2053) + fac11_var_2040(jl_var_2071, lay_var_2058) * absa_var_138(ind1_var_2055 + 1, ig_var_2053)) + tauself_var_2060 + taufor_var_2059
        fracs_var_2049(jl_var_2071, 134 + ig_var_2053, lay_var_2058) = fracrefa_var_136(ig_var_2053)
      END DO
    END DO
  END DO
  DO lay_var_2058 = laytrop_max_var_2062 + 1, klev_var_2034
    DO jl_var_2071 = kidia_var_2032, kfdia_var_2033
      ind0_var_2054 = ((jp_var_2041(jl_var_2071, lay_var_2058) - 13) * 5 + (jt_var_2042(jl_var_2071, lay_var_2058) - 1)) * nspb_var_217(14) + 1
      ind1_var_2055 = ((jp_var_2041(jl_var_2071, lay_var_2058) - 12) * 5 + (jt1_var_2043(jl_var_2071, lay_var_2058) - 1)) * nspb_var_217(14) + 1
      DO ig_var_2053 = 1, 2
        taug_var_2035(jl_var_2071, 134 + ig_var_2053, lay_var_2058) = colco2_var_2044(jl_var_2071, lay_var_2058) * (fac00_var_2037(jl_var_2071, lay_var_2058) * absb_var_139(ind0_var_2054, ig_var_2053) + fac10_var_2039(jl_var_2071, lay_var_2058) * absb_var_139(ind0_var_2054 + 1, ig_var_2053) + fac01_var_2038(jl_var_2071, lay_var_2058) * absb_var_139(ind1_var_2055, ig_var_2053) + fac11_var_2040(jl_var_2071, lay_var_2058) * absb_var_139(ind1_var_2055 + 1, ig_var_2053))
        fracs_var_2049(jl_var_2071, 134 + ig_var_2053, lay_var_2058) = fracrefb_var_137(ig_var_2053)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2062 /= laytrop_min_var_2061) THEN
    DO lay_var_2058 = laytrop_min_var_2061 + 1, laytrop_max_var_2062
      ixc0_var_2068 = ixc_var_2063(lay_var_2058)
      DO ixp_var_2069 = 1, ixc0_var_2068
        jl_var_2071 = ixlow_var_2064(ixp_var_2069, lay_var_2058)
        ind0_var_2054 = ((jp_var_2041(jl_var_2071, lay_var_2058) - 1) * 5 + (jt_var_2042(jl_var_2071, lay_var_2058) - 1)) * nspa_var_216(14) + 1
        ind1_var_2055 = (jp_var_2041(jl_var_2071, lay_var_2058) * 5 + (jt1_var_2043(jl_var_2071, lay_var_2058) - 1)) * nspa_var_216(14) + 1
        inds_var_2056 = indself_var_2048(jl_var_2071, lay_var_2058)
        indf_var_2057 = indfor_var_2050(jl_var_2071, lay_var_2058)
        DO ig_var_2053 = 1, 2
          tauself_var_2060 = selffac_var_2046(jl_var_2071, lay_var_2058) * (selfref_var_140(inds_var_2056, ig_var_2053) + selffrac_var_2047(jl_var_2071, lay_var_2058) * (selfref_var_140(inds_var_2056 + 1, ig_var_2053) - selfref_var_140(inds_var_2056, ig_var_2053)))
          taufor_var_2059 = forfac_var_2051(jl_var_2071, lay_var_2058) * (forref_var_141(indf_var_2057, ig_var_2053) + forfrac_var_2052(jl_var_2071, lay_var_2058) * (forref_var_141(indf_var_2057 + 1, ig_var_2053) - forref_var_141(indf_var_2057, ig_var_2053)))
          taug_var_2035(jl_var_2071, 134 + ig_var_2053, lay_var_2058) = colco2_var_2044(jl_var_2071, lay_var_2058) * (fac00_var_2037(jl_var_2071, lay_var_2058) * absa_var_138(ind0_var_2054, ig_var_2053) + fac10_var_2039(jl_var_2071, lay_var_2058) * absa_var_138(ind0_var_2054 + 1, ig_var_2053) + fac01_var_2038(jl_var_2071, lay_var_2058) * absa_var_138(ind1_var_2055, ig_var_2053) + fac11_var_2040(jl_var_2071, lay_var_2058) * absa_var_138(ind1_var_2055 + 1, ig_var_2053)) + tauself_var_2060 + taufor_var_2059
          fracs_var_2049(jl_var_2071, 134 + ig_var_2053, lay_var_2058) = fracrefa_var_136(ig_var_2053)
        END DO
      END DO
      ixc0_var_2068 = kfdia_var_2033 - kidia_var_2032 + 1 - ixc0_var_2068
      DO ixp_var_2069 = 1, ixc0_var_2068
        jl_var_2071 = ixhigh_var_2065(ixp_var_2069, lay_var_2058)
        ind0_var_2054 = ((jp_var_2041(jl_var_2071, lay_var_2058) - 13) * 5 + (jt_var_2042(jl_var_2071, lay_var_2058) - 1)) * nspb_var_217(14) + 1
        ind1_var_2055 = ((jp_var_2041(jl_var_2071, lay_var_2058) - 12) * 5 + (jt1_var_2043(jl_var_2071, lay_var_2058) - 1)) * nspb_var_217(14) + 1
        DO ig_var_2053 = 1, 2
          taug_var_2035(jl_var_2071, 134 + ig_var_2053, lay_var_2058) = colco2_var_2044(jl_var_2071, lay_var_2058) * (fac00_var_2037(jl_var_2071, lay_var_2058) * absb_var_139(ind0_var_2054, ig_var_2053) + fac10_var_2039(jl_var_2071, lay_var_2058) * absb_var_139(ind0_var_2054 + 1, ig_var_2053) + fac01_var_2038(jl_var_2071, lay_var_2058) * absb_var_139(ind1_var_2055, ig_var_2053) + fac11_var_2040(jl_var_2071, lay_var_2058) * absb_var_139(ind1_var_2055 + 1, ig_var_2053))
          fracs_var_2049(jl_var_2071, 134 + ig_var_2053, lay_var_2058) = fracrefb_var_137(ig_var_2053)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol14
SUBROUTINE rrtm_taumol15(kidia_var_2072, kfdia_var_2073, klev_var_2074, taug_var_2075, p_tauaerl_var_2076, fac00_var_2077, fac01_var_2078, fac10_var_2079, fac11_var_2080, forfac_var_2094, forfrac_var_2095, indfor_var_2093, jp_var_2081, jt_var_2082, jt1_var_2083, oneminus_var_2084, colh2o_var_2085, colco2_var_2086, coln2o_var_2087, laytrop_var_2088, selffac_var_2089, selffrac_var_2090, indself_var_2091, fracs_var_2092, rat_n2oco2, rat_n2oco2_1, minorfrac_var_2096, indminor_var_2097, scaleminor_var_2098, colbrd_var_2099)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrtm, ONLY: ng15
  USE yoerrta15, ONLY: absa_var_143, forref_var_146, fracrefa_var_142, ka_mn2_var_144, selfref_var_145
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2072
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2073
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2074
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2075(kidia_var_2072 : kfdia_var_2073, 140, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2076(kidia_var_2072 : kfdia_var_2073, klev_var_2074, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2077(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2078(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2079(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2080(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2081(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2082(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2083(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2084
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2085(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2086(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_2087(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2088(kidia_var_2072 : kfdia_var_2073)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2089(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2090(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2091(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2092(kidia_var_2072 : kfdia_var_2073, 140, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: rat_n2oco2(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: rat_n2oco2_1(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2093(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2094(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2095(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2096(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2097(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: scaleminor_var_2098(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  REAL(KIND = 8), INTENT(IN) :: colbrd_var_2099(kidia_var_2072 : kfdia_var_2073, klev_var_2074)
  INTEGER(KIND = 4) :: ig_var_2100, ind0_var_2101, ind1_var_2102, inds_var_2103, indf_var_2104, indm_var_2105, js_var_2106, js1_var_2107, jpl_var_2108, jmn2, lay_var_2109
  REAL(KIND = 8) :: refrat_planck_a_var_2110, refrat_m_a_var_2111
  REAL(KIND = 8) :: taufor_var_2112, tauself_var_2113, tau_major_var_2114(2), tau_major1_var_2115(2), n2m1, n2m2, taun2_var_2116, scalen2_var_2117
  REAL(KIND = 8) :: fac000_var_2118, fac100_var_2119, fac200_var_2120, fac010_var_2121, fac110_var_2122, fac210_var_2123, fac001_var_2124, fac101_var_2125, fac201_var_2126, fac011_var_2127, fac111_var_2128, fac211_var_2129
  REAL(KIND = 8) :: p_var_2130, p4_var_2131, fk0_var_2132, fk1_var_2133, fk2_var_2134
  REAL(KIND = 8) :: fs_var_2135, specmult_var_2136, specparm_var_2137, speccomb_var_2138, fs1_var_2139, specmult1_var_2140, specparm1_var_2141, speccomb1_var_2142, fmn2, specmult_mn2, specparm_mn2, speccomb_mn2, fpl_var_2143, specmult_planck_var_2144, specparm_planck_var_2145, speccomb_planck_var_2146
  INTEGER(KIND = 4) :: laytrop_min_var_2147, laytrop_max_var_2148
  INTEGER(KIND = 4) :: ixc_var_2149(klev_var_2074), ixlow_var_2150(kfdia_var_2073, klev_var_2074), ixhigh_var_2151(kfdia_var_2073, klev_var_2074)
  INTEGER(KIND = 4) :: ich_var_2152, icl_var_2153, ixc0_var_2154, ixp_var_2155, jc_var_2156, jl_var_2157
  laytrop_min_var_2147 = MINVAL(laytrop_var_2088)
  laytrop_max_var_2148 = MAXVAL(laytrop_var_2088)
  ixlow_var_2150 = 0
  ixhigh_var_2151 = 0
  ixc_var_2149 = 0
  DO lay_var_2109 = laytrop_min_var_2147 + 1, laytrop_max_var_2148
    icl_var_2153 = 0
    ich_var_2152 = 0
    DO jc_var_2156 = kidia_var_2072, kfdia_var_2073
      IF (lay_var_2109 <= laytrop_var_2088(jc_var_2156)) THEN
        icl_var_2153 = icl_var_2153 + 1
        ixlow_var_2150(icl_var_2153, lay_var_2109) = jc_var_2156
      ELSE
        ich_var_2152 = ich_var_2152 + 1
        ixhigh_var_2151(ich_var_2152, lay_var_2109) = jc_var_2156
      END IF
    END DO
    ixc_var_2149(lay_var_2109) = icl_var_2153
  END DO
  refrat_planck_a_var_2110 = chi_mls(4, 1) / chi_mls(2, 1)
  refrat_m_a_var_2111 = chi_mls(4, 1) / chi_mls(2, 1)
  DO lay_var_2109 = 1, laytrop_min_var_2147
    DO jl_var_2157 = kidia_var_2072, kfdia_var_2073
      speccomb_var_2138 = coln2o_var_2087(jl_var_2157, lay_var_2109) + rat_n2oco2(jl_var_2157, lay_var_2109) * colco2_var_2086(jl_var_2157, lay_var_2109)
      specparm_var_2137 = MIN(coln2o_var_2087(jl_var_2157, lay_var_2109) / speccomb_var_2138, oneminus_var_2084)
      specmult_var_2136 = 8.0D0 * (specparm_var_2137)
      js_var_2106 = 1 + INT(specmult_var_2136)
      fs_var_2135 = ((specmult_var_2136) - AINT((specmult_var_2136)))
      speccomb1_var_2142 = coln2o_var_2087(jl_var_2157, lay_var_2109) + rat_n2oco2_1(jl_var_2157, lay_var_2109) * colco2_var_2086(jl_var_2157, lay_var_2109)
      specparm1_var_2141 = MIN(coln2o_var_2087(jl_var_2157, lay_var_2109) / speccomb1_var_2142, oneminus_var_2084)
      specmult1_var_2140 = 8.0D0 * (specparm1_var_2141)
      js1_var_2107 = 1 + INT(specmult1_var_2140)
      fs1_var_2139 = ((specmult1_var_2140) - AINT((specmult1_var_2140)))
      speccomb_mn2 = coln2o_var_2087(jl_var_2157, lay_var_2109) + refrat_m_a_var_2111 * colco2_var_2086(jl_var_2157, lay_var_2109)
      specparm_mn2 = MIN(coln2o_var_2087(jl_var_2157, lay_var_2109) / speccomb_mn2, oneminus_var_2084)
      specmult_mn2 = 8.0D0 * specparm_mn2
      jmn2 = 1 + INT(specmult_mn2)
      fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))
      speccomb_planck_var_2146 = coln2o_var_2087(jl_var_2157, lay_var_2109) + refrat_planck_a_var_2110 * colco2_var_2086(jl_var_2157, lay_var_2109)
      specparm_planck_var_2145 = MIN(coln2o_var_2087(jl_var_2157, lay_var_2109) / speccomb_planck_var_2146, oneminus_var_2084)
      specmult_planck_var_2144 = 8.0D0 * specparm_planck_var_2145
      jpl_var_2108 = 1 + INT(specmult_planck_var_2144)
      fpl_var_2143 = ((specmult_planck_var_2144) - AINT((specmult_planck_var_2144)))
      ind0_var_2101 = ((jp_var_2081(jl_var_2157, lay_var_2109) - 1) * 5 + (jt_var_2082(jl_var_2157, lay_var_2109) - 1)) * nspa_var_216(15) + js_var_2106
      ind1_var_2102 = (jp_var_2081(jl_var_2157, lay_var_2109) * 5 + (jt1_var_2083(jl_var_2157, lay_var_2109) - 1)) * nspa_var_216(15) + js1_var_2107
      inds_var_2103 = indself_var_2091(jl_var_2157, lay_var_2109)
      indf_var_2104 = indfor_var_2093(jl_var_2157, lay_var_2109)
      indm_var_2105 = indminor_var_2097(jl_var_2157, lay_var_2109)
      scalen2_var_2117 = colbrd_var_2099(jl_var_2157, lay_var_2109) * scaleminor_var_2098(jl_var_2157, lay_var_2109)
      IF (specparm_var_2137 .LT. 0.125D0) THEN
        p_var_2130 = fs_var_2135 - 1.0D0
        p4_var_2131 = p_var_2130 ** 4
        fk0_var_2132 = p4_var_2131
        fk1_var_2133 = 1.0D0 - p_var_2130 - 2.0D0 * p4_var_2131
        fk2_var_2134 = p_var_2130 + p4_var_2131
        fac000_var_2118 = fk0_var_2132 * fac00_var_2077(jl_var_2157, lay_var_2109)
        fac100_var_2119 = fk1_var_2133 * fac00_var_2077(jl_var_2157, lay_var_2109)
        fac200_var_2120 = fk2_var_2134 * fac00_var_2077(jl_var_2157, lay_var_2109)
        fac010_var_2121 = fk0_var_2132 * fac10_var_2079(jl_var_2157, lay_var_2109)
        fac110_var_2122 = fk1_var_2133 * fac10_var_2079(jl_var_2157, lay_var_2109)
        fac210_var_2123 = fk2_var_2134 * fac10_var_2079(jl_var_2157, lay_var_2109)
      ELSE IF (specparm_var_2137 .GT. 0.875D0) THEN
        p_var_2130 = - fs_var_2135
        p4_var_2131 = p_var_2130 ** 4
        fk0_var_2132 = p4_var_2131
        fk1_var_2133 = 1.0D0 - p_var_2130 - 2.0D0 * p4_var_2131
        fk2_var_2134 = p_var_2130 + p4_var_2131
        fac000_var_2118 = fk0_var_2132 * fac00_var_2077(jl_var_2157, lay_var_2109)
        fac100_var_2119 = fk1_var_2133 * fac00_var_2077(jl_var_2157, lay_var_2109)
        fac200_var_2120 = fk2_var_2134 * fac00_var_2077(jl_var_2157, lay_var_2109)
        fac010_var_2121 = fk0_var_2132 * fac10_var_2079(jl_var_2157, lay_var_2109)
        fac110_var_2122 = fk1_var_2133 * fac10_var_2079(jl_var_2157, lay_var_2109)
        fac210_var_2123 = fk2_var_2134 * fac10_var_2079(jl_var_2157, lay_var_2109)
      ELSE
        fac000_var_2118 = (1.0D0 - fs_var_2135) * fac00_var_2077(jl_var_2157, lay_var_2109)
        fac010_var_2121 = (1.0D0 - fs_var_2135) * fac10_var_2079(jl_var_2157, lay_var_2109)
        fac100_var_2119 = fs_var_2135 * fac00_var_2077(jl_var_2157, lay_var_2109)
        fac110_var_2122 = fs_var_2135 * fac10_var_2079(jl_var_2157, lay_var_2109)
        fac200_var_2120 = 0.0D0
        fac210_var_2123 = 0.0D0
      END IF
      IF (specparm1_var_2141 .LT. 0.125D0) THEN
        p_var_2130 = fs1_var_2139 - 1.0D0
        p4_var_2131 = p_var_2130 ** 4
        fk0_var_2132 = p4_var_2131
        fk1_var_2133 = 1.0D0 - p_var_2130 - 2.0D0 * p4_var_2131
        fk2_var_2134 = p_var_2130 + p4_var_2131
        fac001_var_2124 = fk0_var_2132 * fac01_var_2078(jl_var_2157, lay_var_2109)
        fac101_var_2125 = fk1_var_2133 * fac01_var_2078(jl_var_2157, lay_var_2109)
        fac201_var_2126 = fk2_var_2134 * fac01_var_2078(jl_var_2157, lay_var_2109)
        fac011_var_2127 = fk0_var_2132 * fac11_var_2080(jl_var_2157, lay_var_2109)
        fac111_var_2128 = fk1_var_2133 * fac11_var_2080(jl_var_2157, lay_var_2109)
        fac211_var_2129 = fk2_var_2134 * fac11_var_2080(jl_var_2157, lay_var_2109)
      ELSE IF (specparm1_var_2141 .GT. 0.875D0) THEN
        p_var_2130 = - fs1_var_2139
        p4_var_2131 = p_var_2130 ** 4
        fk0_var_2132 = p4_var_2131
        fk1_var_2133 = 1.0D0 - p_var_2130 - 2.0D0 * p4_var_2131
        fk2_var_2134 = p_var_2130 + p4_var_2131
        fac001_var_2124 = fk0_var_2132 * fac01_var_2078(jl_var_2157, lay_var_2109)
        fac101_var_2125 = fk1_var_2133 * fac01_var_2078(jl_var_2157, lay_var_2109)
        fac201_var_2126 = fk2_var_2134 * fac01_var_2078(jl_var_2157, lay_var_2109)
        fac011_var_2127 = fk0_var_2132 * fac11_var_2080(jl_var_2157, lay_var_2109)
        fac111_var_2128 = fk1_var_2133 * fac11_var_2080(jl_var_2157, lay_var_2109)
        fac211_var_2129 = fk2_var_2134 * fac11_var_2080(jl_var_2157, lay_var_2109)
      ELSE
        fac001_var_2124 = (1.0D0 - fs1_var_2139) * fac01_var_2078(jl_var_2157, lay_var_2109)
        fac011_var_2127 = (1.0D0 - fs1_var_2139) * fac11_var_2080(jl_var_2157, lay_var_2109)
        fac101_var_2125 = fs1_var_2139 * fac01_var_2078(jl_var_2157, lay_var_2109)
        fac111_var_2128 = fs1_var_2139 * fac11_var_2080(jl_var_2157, lay_var_2109)
        fac201_var_2126 = 0.0D0
        fac211_var_2129 = 0.0D0
      END IF
      IF (specparm_var_2137 .LT. 0.125D0) THEN
        tau_major_var_2114(1 : ng15) = speccomb_var_2138 * (fac000_var_2118 * absa_var_143(ind0_var_2101, 1 : 2) + fac100_var_2119 * absa_var_143(ind0_var_2101 + 1, 1 : 2) + fac200_var_2120 * absa_var_143(ind0_var_2101 + 2, 1 : 2) + fac010_var_2121 * absa_var_143(ind0_var_2101 + 9, 1 : 2) + fac110_var_2122 * absa_var_143(ind0_var_2101 + 10, 1 : 2) + fac210_var_2123 * absa_var_143(ind0_var_2101 + 11, 1 : 2))
      ELSE IF (specparm_var_2137 .GT. 0.875D0) THEN
        tau_major_var_2114(1 : ng15) = speccomb_var_2138 * (fac200_var_2120 * absa_var_143(ind0_var_2101 - 1, 1 : 2) + fac100_var_2119 * absa_var_143(ind0_var_2101, 1 : 2) + fac000_var_2118 * absa_var_143(ind0_var_2101 + 1, 1 : 2) + fac210_var_2123 * absa_var_143(ind0_var_2101 + 8, 1 : 2) + fac110_var_2122 * absa_var_143(ind0_var_2101 + 9, 1 : 2) + fac010_var_2121 * absa_var_143(ind0_var_2101 + 10, 1 : 2))
      ELSE
        tau_major_var_2114(1 : ng15) = speccomb_var_2138 * (fac000_var_2118 * absa_var_143(ind0_var_2101, 1 : 2) + fac100_var_2119 * absa_var_143(ind0_var_2101 + 1, 1 : 2) + fac010_var_2121 * absa_var_143(ind0_var_2101 + 9, 1 : 2) + fac110_var_2122 * absa_var_143(ind0_var_2101 + 10, 1 : 2))
      END IF
      IF (specparm1_var_2141 .LT. 0.125D0) THEN
        tau_major1_var_2115(1 : ng15) = speccomb1_var_2142 * (fac001_var_2124 * absa_var_143(ind1_var_2102, 1 : 2) + fac101_var_2125 * absa_var_143(ind1_var_2102 + 1, 1 : 2) + fac201_var_2126 * absa_var_143(ind1_var_2102 + 2, 1 : 2) + fac011_var_2127 * absa_var_143(ind1_var_2102 + 9, 1 : 2) + fac111_var_2128 * absa_var_143(ind1_var_2102 + 10, 1 : 2) + fac211_var_2129 * absa_var_143(ind1_var_2102 + 11, 1 : 2))
      ELSE IF (specparm1_var_2141 .GT. 0.875D0) THEN
        tau_major1_var_2115(1 : ng15) = speccomb1_var_2142 * (fac201_var_2126 * absa_var_143(ind1_var_2102 - 1, 1 : 2) + fac101_var_2125 * absa_var_143(ind1_var_2102, 1 : 2) + fac001_var_2124 * absa_var_143(ind1_var_2102 + 1, 1 : 2) + fac211_var_2129 * absa_var_143(ind1_var_2102 + 8, 1 : 2) + fac111_var_2128 * absa_var_143(ind1_var_2102 + 9, 1 : 2) + fac011_var_2127 * absa_var_143(ind1_var_2102 + 10, 1 : 2))
      ELSE
        tau_major1_var_2115(1 : ng15) = speccomb1_var_2142 * (fac001_var_2124 * absa_var_143(ind1_var_2102, 1 : 2) + fac101_var_2125 * absa_var_143(ind1_var_2102 + 1, 1 : 2) + fac011_var_2127 * absa_var_143(ind1_var_2102 + 9, 1 : 2) + fac111_var_2128 * absa_var_143(ind1_var_2102 + 10, 1 : 2))
      END IF
      DO ig_var_2100 = 1, 2
        tauself_var_2113 = selffac_var_2089(jl_var_2157, lay_var_2109) * (selfref_var_145(inds_var_2103, ig_var_2100) + selffrac_var_2090(jl_var_2157, lay_var_2109) * (selfref_var_145(inds_var_2103 + 1, ig_var_2100) - selfref_var_145(inds_var_2103, ig_var_2100)))
        taufor_var_2112 = forfac_var_2094(jl_var_2157, lay_var_2109) * (forref_var_146(indf_var_2104, ig_var_2100) + forfrac_var_2095(jl_var_2157, lay_var_2109) * (forref_var_146(indf_var_2104 + 1, ig_var_2100) - forref_var_146(indf_var_2104, ig_var_2100)))
        n2m1 = ka_mn2_var_144(jmn2, indm_var_2105, ig_var_2100) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_2105, ig_var_2100) - ka_mn2_var_144(jmn2, indm_var_2105, ig_var_2100))
        n2m2 = ka_mn2_var_144(jmn2, indm_var_2105 + 1, ig_var_2100) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_2105 + 1, ig_var_2100) - ka_mn2_var_144(jmn2, indm_var_2105 + 1, ig_var_2100))
        taun2_var_2116 = scalen2_var_2117 * (n2m1 + minorfrac_var_2096(jl_var_2157, lay_var_2109) * (n2m2 - n2m1))
        taug_var_2075(jl_var_2157, 136 + ig_var_2100, lay_var_2109) = tau_major_var_2114(ig_var_2100) + tau_major1_var_2115(ig_var_2100) + tauself_var_2113 + taufor_var_2112 + taun2_var_2116
        fracs_var_2092(jl_var_2157, 136 + ig_var_2100, lay_var_2109) = fracrefa_var_142(ig_var_2100, jpl_var_2108) + fpl_var_2143 * (fracrefa_var_142(ig_var_2100, jpl_var_2108 + 1) - fracrefa_var_142(ig_var_2100, jpl_var_2108))
      END DO
    END DO
  END DO
  DO ig_var_2100 = 1, 2
    DO lay_var_2109 = laytrop_max_var_2148 + 1, klev_var_2074
      DO jl_var_2157 = kidia_var_2072, kfdia_var_2073
        taug_var_2075(jl_var_2157, 136 + ig_var_2100, lay_var_2109) = 0.0D0
        fracs_var_2092(jl_var_2157, 136 + ig_var_2100, lay_var_2109) = 0.0D0
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2148 /= laytrop_min_var_2147) THEN
    DO lay_var_2109 = laytrop_min_var_2147 + 1, laytrop_max_var_2148
      ixc0_var_2154 = ixc_var_2149(lay_var_2109)
      DO ixp_var_2155 = 1, ixc0_var_2154
        jl_var_2157 = ixlow_var_2150(ixp_var_2155, lay_var_2109)
        speccomb_var_2138 = coln2o_var_2087(jl_var_2157, lay_var_2109) + rat_n2oco2(jl_var_2157, lay_var_2109) * colco2_var_2086(jl_var_2157, lay_var_2109)
        specparm_var_2137 = MIN(coln2o_var_2087(jl_var_2157, lay_var_2109) / speccomb_var_2138, oneminus_var_2084)
        specmult_var_2136 = 8.0D0 * (specparm_var_2137)
        js_var_2106 = 1 + INT(specmult_var_2136)
        fs_var_2135 = ((specmult_var_2136) - AINT((specmult_var_2136)))
        speccomb1_var_2142 = coln2o_var_2087(jl_var_2157, lay_var_2109) + rat_n2oco2_1(jl_var_2157, lay_var_2109) * colco2_var_2086(jl_var_2157, lay_var_2109)
        specparm1_var_2141 = MIN(coln2o_var_2087(jl_var_2157, lay_var_2109) / speccomb1_var_2142, oneminus_var_2084)
        specmult1_var_2140 = 8.0D0 * (specparm1_var_2141)
        js1_var_2107 = 1 + INT(specmult1_var_2140)
        fs1_var_2139 = ((specmult1_var_2140) - AINT((specmult1_var_2140)))
        speccomb_mn2 = coln2o_var_2087(jl_var_2157, lay_var_2109) + refrat_m_a_var_2111 * colco2_var_2086(jl_var_2157, lay_var_2109)
        specparm_mn2 = MIN(coln2o_var_2087(jl_var_2157, lay_var_2109) / speccomb_mn2, oneminus_var_2084)
        specmult_mn2 = 8.0D0 * specparm_mn2
        jmn2 = 1 + INT(specmult_mn2)
        fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))
        speccomb_planck_var_2146 = coln2o_var_2087(jl_var_2157, lay_var_2109) + refrat_planck_a_var_2110 * colco2_var_2086(jl_var_2157, lay_var_2109)
        specparm_planck_var_2145 = MIN(coln2o_var_2087(jl_var_2157, lay_var_2109) / speccomb_planck_var_2146, oneminus_var_2084)
        specmult_planck_var_2144 = 8.0D0 * specparm_planck_var_2145
        jpl_var_2108 = 1 + INT(specmult_planck_var_2144)
        fpl_var_2143 = ((specmult_planck_var_2144) - AINT((specmult_planck_var_2144)))
        ind0_var_2101 = ((jp_var_2081(jl_var_2157, lay_var_2109) - 1) * 5 + (jt_var_2082(jl_var_2157, lay_var_2109) - 1)) * nspa_var_216(15) + js_var_2106
        ind1_var_2102 = (jp_var_2081(jl_var_2157, lay_var_2109) * 5 + (jt1_var_2083(jl_var_2157, lay_var_2109) - 1)) * nspa_var_216(15) + js1_var_2107
        inds_var_2103 = indself_var_2091(jl_var_2157, lay_var_2109)
        indf_var_2104 = indfor_var_2093(jl_var_2157, lay_var_2109)
        indm_var_2105 = indminor_var_2097(jl_var_2157, lay_var_2109)
        scalen2_var_2117 = colbrd_var_2099(jl_var_2157, lay_var_2109) * scaleminor_var_2098(jl_var_2157, lay_var_2109)
        IF (specparm_var_2137 .LT. 0.125D0) THEN
          p_var_2130 = fs_var_2135 - 1.0D0
          p4_var_2131 = p_var_2130 ** 4
          fk0_var_2132 = p4_var_2131
          fk1_var_2133 = 1.0D0 - p_var_2130 - 2.0D0 * p4_var_2131
          fk2_var_2134 = p_var_2130 + p4_var_2131
          fac000_var_2118 = fk0_var_2132 * fac00_var_2077(jl_var_2157, lay_var_2109)
          fac100_var_2119 = fk1_var_2133 * fac00_var_2077(jl_var_2157, lay_var_2109)
          fac200_var_2120 = fk2_var_2134 * fac00_var_2077(jl_var_2157, lay_var_2109)
          fac010_var_2121 = fk0_var_2132 * fac10_var_2079(jl_var_2157, lay_var_2109)
          fac110_var_2122 = fk1_var_2133 * fac10_var_2079(jl_var_2157, lay_var_2109)
          fac210_var_2123 = fk2_var_2134 * fac10_var_2079(jl_var_2157, lay_var_2109)
        ELSE IF (specparm_var_2137 .GT. 0.875D0) THEN
          p_var_2130 = - fs_var_2135
          p4_var_2131 = p_var_2130 ** 4
          fk0_var_2132 = p4_var_2131
          fk1_var_2133 = 1.0D0 - p_var_2130 - 2.0D0 * p4_var_2131
          fk2_var_2134 = p_var_2130 + p4_var_2131
          fac000_var_2118 = fk0_var_2132 * fac00_var_2077(jl_var_2157, lay_var_2109)
          fac100_var_2119 = fk1_var_2133 * fac00_var_2077(jl_var_2157, lay_var_2109)
          fac200_var_2120 = fk2_var_2134 * fac00_var_2077(jl_var_2157, lay_var_2109)
          fac010_var_2121 = fk0_var_2132 * fac10_var_2079(jl_var_2157, lay_var_2109)
          fac110_var_2122 = fk1_var_2133 * fac10_var_2079(jl_var_2157, lay_var_2109)
          fac210_var_2123 = fk2_var_2134 * fac10_var_2079(jl_var_2157, lay_var_2109)
        ELSE
          fac000_var_2118 = (1.0D0 - fs_var_2135) * fac00_var_2077(jl_var_2157, lay_var_2109)
          fac010_var_2121 = (1.0D0 - fs_var_2135) * fac10_var_2079(jl_var_2157, lay_var_2109)
          fac100_var_2119 = fs_var_2135 * fac00_var_2077(jl_var_2157, lay_var_2109)
          fac110_var_2122 = fs_var_2135 * fac10_var_2079(jl_var_2157, lay_var_2109)
          fac200_var_2120 = 0.0D0
          fac210_var_2123 = 0.0D0
        END IF
        IF (specparm1_var_2141 .LT. 0.125D0) THEN
          p_var_2130 = fs1_var_2139 - 1.0D0
          p4_var_2131 = p_var_2130 ** 4
          fk0_var_2132 = p4_var_2131
          fk1_var_2133 = 1.0D0 - p_var_2130 - 2.0D0 * p4_var_2131
          fk2_var_2134 = p_var_2130 + p4_var_2131
          fac001_var_2124 = fk0_var_2132 * fac01_var_2078(jl_var_2157, lay_var_2109)
          fac101_var_2125 = fk1_var_2133 * fac01_var_2078(jl_var_2157, lay_var_2109)
          fac201_var_2126 = fk2_var_2134 * fac01_var_2078(jl_var_2157, lay_var_2109)
          fac011_var_2127 = fk0_var_2132 * fac11_var_2080(jl_var_2157, lay_var_2109)
          fac111_var_2128 = fk1_var_2133 * fac11_var_2080(jl_var_2157, lay_var_2109)
          fac211_var_2129 = fk2_var_2134 * fac11_var_2080(jl_var_2157, lay_var_2109)
        ELSE IF (specparm1_var_2141 .GT. 0.875D0) THEN
          p_var_2130 = - fs1_var_2139
          p4_var_2131 = p_var_2130 ** 4
          fk0_var_2132 = p4_var_2131
          fk1_var_2133 = 1.0D0 - p_var_2130 - 2.0D0 * p4_var_2131
          fk2_var_2134 = p_var_2130 + p4_var_2131
          fac001_var_2124 = fk0_var_2132 * fac01_var_2078(jl_var_2157, lay_var_2109)
          fac101_var_2125 = fk1_var_2133 * fac01_var_2078(jl_var_2157, lay_var_2109)
          fac201_var_2126 = fk2_var_2134 * fac01_var_2078(jl_var_2157, lay_var_2109)
          fac011_var_2127 = fk0_var_2132 * fac11_var_2080(jl_var_2157, lay_var_2109)
          fac111_var_2128 = fk1_var_2133 * fac11_var_2080(jl_var_2157, lay_var_2109)
          fac211_var_2129 = fk2_var_2134 * fac11_var_2080(jl_var_2157, lay_var_2109)
        ELSE
          fac001_var_2124 = (1.0D0 - fs1_var_2139) * fac01_var_2078(jl_var_2157, lay_var_2109)
          fac011_var_2127 = (1.0D0 - fs1_var_2139) * fac11_var_2080(jl_var_2157, lay_var_2109)
          fac101_var_2125 = fs1_var_2139 * fac01_var_2078(jl_var_2157, lay_var_2109)
          fac111_var_2128 = fs1_var_2139 * fac11_var_2080(jl_var_2157, lay_var_2109)
          fac201_var_2126 = 0.0D0
          fac211_var_2129 = 0.0D0
        END IF
        IF (specparm_var_2137 .LT. 0.125D0) THEN
          tau_major_var_2114(1 : ng15) = speccomb_var_2138 * (fac000_var_2118 * absa_var_143(ind0_var_2101, 1 : 2) + fac100_var_2119 * absa_var_143(ind0_var_2101 + 1, 1 : 2) + fac200_var_2120 * absa_var_143(ind0_var_2101 + 2, 1 : 2) + fac010_var_2121 * absa_var_143(ind0_var_2101 + 9, 1 : 2) + fac110_var_2122 * absa_var_143(ind0_var_2101 + 10, 1 : 2) + fac210_var_2123 * absa_var_143(ind0_var_2101 + 11, 1 : 2))
        ELSE IF (specparm_var_2137 .GT. 0.875D0) THEN
          tau_major_var_2114(1 : ng15) = speccomb_var_2138 * (fac200_var_2120 * absa_var_143(ind0_var_2101 - 1, 1 : 2) + fac100_var_2119 * absa_var_143(ind0_var_2101, 1 : 2) + fac000_var_2118 * absa_var_143(ind0_var_2101 + 1, 1 : 2) + fac210_var_2123 * absa_var_143(ind0_var_2101 + 8, 1 : 2) + fac110_var_2122 * absa_var_143(ind0_var_2101 + 9, 1 : 2) + fac010_var_2121 * absa_var_143(ind0_var_2101 + 10, 1 : 2))
        ELSE
          tau_major_var_2114(1 : ng15) = speccomb_var_2138 * (fac000_var_2118 * absa_var_143(ind0_var_2101, 1 : 2) + fac100_var_2119 * absa_var_143(ind0_var_2101 + 1, 1 : 2) + fac010_var_2121 * absa_var_143(ind0_var_2101 + 9, 1 : 2) + fac110_var_2122 * absa_var_143(ind0_var_2101 + 10, 1 : 2))
        END IF
        IF (specparm1_var_2141 .LT. 0.125D0) THEN
          tau_major1_var_2115(1 : ng15) = speccomb1_var_2142 * (fac001_var_2124 * absa_var_143(ind1_var_2102, 1 : 2) + fac101_var_2125 * absa_var_143(ind1_var_2102 + 1, 1 : 2) + fac201_var_2126 * absa_var_143(ind1_var_2102 + 2, 1 : 2) + fac011_var_2127 * absa_var_143(ind1_var_2102 + 9, 1 : 2) + fac111_var_2128 * absa_var_143(ind1_var_2102 + 10, 1 : 2) + fac211_var_2129 * absa_var_143(ind1_var_2102 + 11, 1 : 2))
        ELSE IF (specparm1_var_2141 .GT. 0.875D0) THEN
          tau_major1_var_2115(1 : ng15) = speccomb1_var_2142 * (fac201_var_2126 * absa_var_143(ind1_var_2102 - 1, 1 : 2) + fac101_var_2125 * absa_var_143(ind1_var_2102, 1 : 2) + fac001_var_2124 * absa_var_143(ind1_var_2102 + 1, 1 : 2) + fac211_var_2129 * absa_var_143(ind1_var_2102 + 8, 1 : 2) + fac111_var_2128 * absa_var_143(ind1_var_2102 + 9, 1 : 2) + fac011_var_2127 * absa_var_143(ind1_var_2102 + 10, 1 : 2))
        ELSE
          tau_major1_var_2115(1 : ng15) = speccomb1_var_2142 * (fac001_var_2124 * absa_var_143(ind1_var_2102, 1 : 2) + fac101_var_2125 * absa_var_143(ind1_var_2102 + 1, 1 : 2) + fac011_var_2127 * absa_var_143(ind1_var_2102 + 9, 1 : 2) + fac111_var_2128 * absa_var_143(ind1_var_2102 + 10, 1 : 2))
        END IF
        DO ig_var_2100 = 1, 2
          tauself_var_2113 = selffac_var_2089(jl_var_2157, lay_var_2109) * (selfref_var_145(inds_var_2103, ig_var_2100) + selffrac_var_2090(jl_var_2157, lay_var_2109) * (selfref_var_145(inds_var_2103 + 1, ig_var_2100) - selfref_var_145(inds_var_2103, ig_var_2100)))
          taufor_var_2112 = forfac_var_2094(jl_var_2157, lay_var_2109) * (forref_var_146(indf_var_2104, ig_var_2100) + forfrac_var_2095(jl_var_2157, lay_var_2109) * (forref_var_146(indf_var_2104 + 1, ig_var_2100) - forref_var_146(indf_var_2104, ig_var_2100)))
          n2m1 = ka_mn2_var_144(jmn2, indm_var_2105, ig_var_2100) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_2105, ig_var_2100) - ka_mn2_var_144(jmn2, indm_var_2105, ig_var_2100))
          n2m2 = ka_mn2_var_144(jmn2, indm_var_2105 + 1, ig_var_2100) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_2105 + 1, ig_var_2100) - ka_mn2_var_144(jmn2, indm_var_2105 + 1, ig_var_2100))
          taun2_var_2116 = scalen2_var_2117 * (n2m1 + minorfrac_var_2096(jl_var_2157, lay_var_2109) * (n2m2 - n2m1))
          taug_var_2075(jl_var_2157, 136 + ig_var_2100, lay_var_2109) = tau_major_var_2114(ig_var_2100) + tau_major1_var_2115(ig_var_2100) + tauself_var_2113 + taufor_var_2112 + taun2_var_2116
          fracs_var_2092(jl_var_2157, 136 + ig_var_2100, lay_var_2109) = fracrefa_var_142(ig_var_2100, jpl_var_2108) + fpl_var_2143 * (fracrefa_var_142(ig_var_2100, jpl_var_2108 + 1) - fracrefa_var_142(ig_var_2100, jpl_var_2108))
        END DO
      END DO
      ixc0_var_2154 = kfdia_var_2073 - kidia_var_2072 + 1 - ixc0_var_2154
      DO ig_var_2100 = 1, 2
        DO ixp_var_2155 = 1, ixc0_var_2154
          jl_var_2157 = ixhigh_var_2151(ixp_var_2155, lay_var_2109)
          taug_var_2075(jl_var_2157, 136 + ig_var_2100, lay_var_2109) = 0.0D0
          fracs_var_2092(jl_var_2157, 136 + ig_var_2100, lay_var_2109) = 0.0D0
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol15
SUBROUTINE rrtm_taumol16(kidia_var_2158, kfdia_var_2159, klev_var_2160, taug_var_2161, p_tauaerl_var_2162, fac00_var_2163, fac01_var_2164, fac10_var_2165, fac11_var_2166, forfac_var_2181, forfrac_var_2182, indfor_var_2180, jp_var_2167, jt_var_2168, jt1_var_2169, oneminus_var_2170, colh2o_var_2171, colch4_var_2172, laytrop_var_2173, selffac_var_2174, selffrac_var_2175, indself_var_2176, fracs_var_2177, rat_h2och4_var_2178, rat_h2och4_1_var_2179)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng16
  USE yoerrta16, ONLY: absa_var_149, absb_var_150, forref_var_152, fracrefa_var_147, fracrefb_var_148, selfref_var_151
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2158
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2159
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2160
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2161(kidia_var_2158 : kfdia_var_2159, 140, klev_var_2160)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2162(kidia_var_2158 : kfdia_var_2159, klev_var_2160, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2163(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2164(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2165(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2166(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2167(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2168(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2169(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2170
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2171(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  REAL(KIND = 8), INTENT(IN) :: colch4_var_2172(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2173(kidia_var_2158 : kfdia_var_2159)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2174(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2175(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2176(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2177(kidia_var_2158 : kfdia_var_2159, 140, klev_var_2160)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_var_2178(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_1_var_2179(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2180(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2181(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2182(kidia_var_2158 : kfdia_var_2159, klev_var_2160)
  INTEGER(KIND = 4) :: ig_var_2183, ind0_var_2184, ind1_var_2185, inds_var_2186, indf_var_2187, js_var_2188, js1_var_2189, jpl_var_2190, lay_var_2191
  REAL(KIND = 8) :: fac000_var_2192, fac100_var_2193, fac200_var_2194, fac010_var_2195, fac110_var_2196, fac210_var_2197, fac001_var_2198, fac101_var_2199, fac201_var_2200, fac011_var_2201, fac111_var_2202, fac211_var_2203
  REAL(KIND = 8) :: p_var_2204, p4_var_2205, fk0_var_2206, fk1_var_2207, fk2_var_2208
  REAL(KIND = 8) :: refrat_planck_a_var_2209
  REAL(KIND = 8) :: taufor_var_2210, tauself_var_2211, tau_major_var_2212(2), tau_major1_var_2213(2)
  REAL(KIND = 8) :: fs_var_2214, specmult_var_2215, specparm_var_2216, speccomb_var_2217, fs1_var_2218, specmult1_var_2219, specparm1_var_2220, speccomb1_var_2221, fpl_var_2222, specmult_planck_var_2223, specparm_planck_var_2224, speccomb_planck_var_2225
  INTEGER(KIND = 4) :: laytrop_min_var_2226, laytrop_max_var_2227
  INTEGER(KIND = 4) :: ixc_var_2228(klev_var_2160), ixlow_var_2229(kfdia_var_2159, klev_var_2160), ixhigh_var_2230(kfdia_var_2159, klev_var_2160)
  INTEGER(KIND = 4) :: ich_var_2231, icl_var_2232, ixc0_var_2233, ixp_var_2234, jc_var_2235, jl_var_2236
  laytrop_min_var_2226 = MINVAL(laytrop_var_2173)
  laytrop_max_var_2227 = MAXVAL(laytrop_var_2173)
  ixlow_var_2229 = 0
  ixhigh_var_2230 = 0
  ixc_var_2228 = 0
  DO lay_var_2191 = laytrop_min_var_2226 + 1, laytrop_max_var_2227
    icl_var_2232 = 0
    ich_var_2231 = 0
    DO jc_var_2235 = kidia_var_2158, kfdia_var_2159
      IF (lay_var_2191 <= laytrop_var_2173(jc_var_2235)) THEN
        icl_var_2232 = icl_var_2232 + 1
        ixlow_var_2229(icl_var_2232, lay_var_2191) = jc_var_2235
      ELSE
        ich_var_2231 = ich_var_2231 + 1
        ixhigh_var_2230(ich_var_2231, lay_var_2191) = jc_var_2235
      END IF
    END DO
    ixc_var_2228(lay_var_2191) = icl_var_2232
  END DO
  refrat_planck_a_var_2209 = chi_mls(1, 6) / chi_mls(6, 6)
  DO lay_var_2191 = 1, laytrop_min_var_2226
    DO jl_var_2236 = kidia_var_2158, kfdia_var_2159
      speccomb_var_2217 = colh2o_var_2171(jl_var_2236, lay_var_2191) + rat_h2och4_var_2178(jl_var_2236, lay_var_2191) * colch4_var_2172(jl_var_2236, lay_var_2191)
      specparm_var_2216 = MIN(colh2o_var_2171(jl_var_2236, lay_var_2191) / speccomb_var_2217, oneminus_var_2170)
      specmult_var_2215 = 8.0D0 * (specparm_var_2216)
      js_var_2188 = 1 + INT(specmult_var_2215)
      fs_var_2214 = ((specmult_var_2215) - AINT((specmult_var_2215)))
      speccomb1_var_2221 = colh2o_var_2171(jl_var_2236, lay_var_2191) + rat_h2och4_1_var_2179(jl_var_2236, lay_var_2191) * colch4_var_2172(jl_var_2236, lay_var_2191)
      specparm1_var_2220 = MIN(colh2o_var_2171(jl_var_2236, lay_var_2191) / speccomb1_var_2221, oneminus_var_2170)
      specmult1_var_2219 = 8.0D0 * (specparm1_var_2220)
      js1_var_2189 = 1 + INT(specmult1_var_2219)
      fs1_var_2218 = ((specmult1_var_2219) - AINT((specmult1_var_2219)))
      speccomb_planck_var_2225 = colh2o_var_2171(jl_var_2236, lay_var_2191) + refrat_planck_a_var_2209 * colch4_var_2172(jl_var_2236, lay_var_2191)
      specparm_planck_var_2224 = MIN(colh2o_var_2171(jl_var_2236, lay_var_2191) / speccomb_planck_var_2225, oneminus_var_2170)
      specmult_planck_var_2223 = 8.0D0 * specparm_planck_var_2224
      jpl_var_2190 = 1 + INT(specmult_planck_var_2223)
      fpl_var_2222 = ((specmult_planck_var_2223) - AINT((specmult_planck_var_2223)))
      ind0_var_2184 = ((jp_var_2167(jl_var_2236, lay_var_2191) - 1) * 5 + (jt_var_2168(jl_var_2236, lay_var_2191) - 1)) * nspa_var_216(16) + js_var_2188
      ind1_var_2185 = (jp_var_2167(jl_var_2236, lay_var_2191) * 5 + (jt1_var_2169(jl_var_2236, lay_var_2191) - 1)) * nspa_var_216(16) + js1_var_2189
      inds_var_2186 = indself_var_2176(jl_var_2236, lay_var_2191)
      indf_var_2187 = indfor_var_2180(jl_var_2236, lay_var_2191)
      IF (specparm_var_2216 .LT. 0.125D0) THEN
        p_var_2204 = fs_var_2214 - 1.0D0
        p4_var_2205 = p_var_2204 ** 4
        fk0_var_2206 = p4_var_2205
        fk1_var_2207 = 1.0D0 - p_var_2204 - 2.0D0 * p4_var_2205
        fk2_var_2208 = p_var_2204 + p4_var_2205
        fac000_var_2192 = fk0_var_2206 * fac00_var_2163(jl_var_2236, lay_var_2191)
        fac100_var_2193 = fk1_var_2207 * fac00_var_2163(jl_var_2236, lay_var_2191)
        fac200_var_2194 = fk2_var_2208 * fac00_var_2163(jl_var_2236, lay_var_2191)
        fac010_var_2195 = fk0_var_2206 * fac10_var_2165(jl_var_2236, lay_var_2191)
        fac110_var_2196 = fk1_var_2207 * fac10_var_2165(jl_var_2236, lay_var_2191)
        fac210_var_2197 = fk2_var_2208 * fac10_var_2165(jl_var_2236, lay_var_2191)
      ELSE IF (specparm_var_2216 .GT. 0.875D0) THEN
        p_var_2204 = - fs_var_2214
        p4_var_2205 = p_var_2204 ** 4
        fk0_var_2206 = p4_var_2205
        fk1_var_2207 = 1.0D0 - p_var_2204 - 2.0D0 * p4_var_2205
        fk2_var_2208 = p_var_2204 + p4_var_2205
        fac000_var_2192 = fk0_var_2206 * fac00_var_2163(jl_var_2236, lay_var_2191)
        fac100_var_2193 = fk1_var_2207 * fac00_var_2163(jl_var_2236, lay_var_2191)
        fac200_var_2194 = fk2_var_2208 * fac00_var_2163(jl_var_2236, lay_var_2191)
        fac010_var_2195 = fk0_var_2206 * fac10_var_2165(jl_var_2236, lay_var_2191)
        fac110_var_2196 = fk1_var_2207 * fac10_var_2165(jl_var_2236, lay_var_2191)
        fac210_var_2197 = fk2_var_2208 * fac10_var_2165(jl_var_2236, lay_var_2191)
      ELSE
        fac000_var_2192 = (1.0D0 - fs_var_2214) * fac00_var_2163(jl_var_2236, lay_var_2191)
        fac010_var_2195 = (1.0D0 - fs_var_2214) * fac10_var_2165(jl_var_2236, lay_var_2191)
        fac100_var_2193 = fs_var_2214 * fac00_var_2163(jl_var_2236, lay_var_2191)
        fac110_var_2196 = fs_var_2214 * fac10_var_2165(jl_var_2236, lay_var_2191)
        fac200_var_2194 = 0.0D0
        fac210_var_2197 = 0.0D0
      END IF
      IF (specparm1_var_2220 .LT. 0.125D0) THEN
        p_var_2204 = fs1_var_2218 - 1.0D0
        p4_var_2205 = p_var_2204 ** 4
        fk0_var_2206 = p4_var_2205
        fk1_var_2207 = 1.0D0 - p_var_2204 - 2.0D0 * p4_var_2205
        fk2_var_2208 = p_var_2204 + p4_var_2205
        fac001_var_2198 = fk0_var_2206 * fac01_var_2164(jl_var_2236, lay_var_2191)
        fac101_var_2199 = fk1_var_2207 * fac01_var_2164(jl_var_2236, lay_var_2191)
        fac201_var_2200 = fk2_var_2208 * fac01_var_2164(jl_var_2236, lay_var_2191)
        fac011_var_2201 = fk0_var_2206 * fac11_var_2166(jl_var_2236, lay_var_2191)
        fac111_var_2202 = fk1_var_2207 * fac11_var_2166(jl_var_2236, lay_var_2191)
        fac211_var_2203 = fk2_var_2208 * fac11_var_2166(jl_var_2236, lay_var_2191)
      ELSE IF (specparm1_var_2220 .GT. 0.875D0) THEN
        p_var_2204 = - fs1_var_2218
        p4_var_2205 = p_var_2204 ** 4
        fk0_var_2206 = p4_var_2205
        fk1_var_2207 = 1.0D0 - p_var_2204 - 2.0D0 * p4_var_2205
        fk2_var_2208 = p_var_2204 + p4_var_2205
        fac001_var_2198 = fk0_var_2206 * fac01_var_2164(jl_var_2236, lay_var_2191)
        fac101_var_2199 = fk1_var_2207 * fac01_var_2164(jl_var_2236, lay_var_2191)
        fac201_var_2200 = fk2_var_2208 * fac01_var_2164(jl_var_2236, lay_var_2191)
        fac011_var_2201 = fk0_var_2206 * fac11_var_2166(jl_var_2236, lay_var_2191)
        fac111_var_2202 = fk1_var_2207 * fac11_var_2166(jl_var_2236, lay_var_2191)
        fac211_var_2203 = fk2_var_2208 * fac11_var_2166(jl_var_2236, lay_var_2191)
      ELSE
        fac001_var_2198 = (1.0D0 - fs1_var_2218) * fac01_var_2164(jl_var_2236, lay_var_2191)
        fac011_var_2201 = (1.0D0 - fs1_var_2218) * fac11_var_2166(jl_var_2236, lay_var_2191)
        fac101_var_2199 = fs1_var_2218 * fac01_var_2164(jl_var_2236, lay_var_2191)
        fac111_var_2202 = fs1_var_2218 * fac11_var_2166(jl_var_2236, lay_var_2191)
        fac201_var_2200 = 0.0D0
        fac211_var_2203 = 0.0D0
      END IF
      IF (specparm_var_2216 .LT. 0.125D0) THEN
        tau_major_var_2212(1 : ng16) = speccomb_var_2217 * (fac000_var_2192 * absa_var_149(ind0_var_2184, 1 : 2) + fac100_var_2193 * absa_var_149(ind0_var_2184 + 1, 1 : 2) + fac200_var_2194 * absa_var_149(ind0_var_2184 + 2, 1 : 2) + fac010_var_2195 * absa_var_149(ind0_var_2184 + 9, 1 : 2) + fac110_var_2196 * absa_var_149(ind0_var_2184 + 10, 1 : 2) + fac210_var_2197 * absa_var_149(ind0_var_2184 + 11, 1 : 2))
      ELSE IF (specparm_var_2216 .GT. 0.875D0) THEN
        tau_major_var_2212(1 : ng16) = speccomb_var_2217 * (fac200_var_2194 * absa_var_149(ind0_var_2184 - 1, 1 : 2) + fac100_var_2193 * absa_var_149(ind0_var_2184, 1 : 2) + fac000_var_2192 * absa_var_149(ind0_var_2184 + 1, 1 : 2) + fac210_var_2197 * absa_var_149(ind0_var_2184 + 8, 1 : 2) + fac110_var_2196 * absa_var_149(ind0_var_2184 + 9, 1 : 2) + fac010_var_2195 * absa_var_149(ind0_var_2184 + 10, 1 : 2))
      ELSE
        tau_major_var_2212(1 : ng16) = speccomb_var_2217 * (fac000_var_2192 * absa_var_149(ind0_var_2184, 1 : 2) + fac100_var_2193 * absa_var_149(ind0_var_2184 + 1, 1 : 2) + fac010_var_2195 * absa_var_149(ind0_var_2184 + 9, 1 : 2) + fac110_var_2196 * absa_var_149(ind0_var_2184 + 10, 1 : 2))
      END IF
      IF (specparm1_var_2220 .LT. 0.125D0) THEN
        tau_major1_var_2213(1 : ng16) = speccomb1_var_2221 * (fac001_var_2198 * absa_var_149(ind1_var_2185, 1 : 2) + fac101_var_2199 * absa_var_149(ind1_var_2185 + 1, 1 : 2) + fac201_var_2200 * absa_var_149(ind1_var_2185 + 2, 1 : 2) + fac011_var_2201 * absa_var_149(ind1_var_2185 + 9, 1 : 2) + fac111_var_2202 * absa_var_149(ind1_var_2185 + 10, 1 : 2) + fac211_var_2203 * absa_var_149(ind1_var_2185 + 11, 1 : 2))
      ELSE IF (specparm1_var_2220 .GT. 0.875D0) THEN
        tau_major1_var_2213(1 : ng16) = speccomb1_var_2221 * (fac201_var_2200 * absa_var_149(ind1_var_2185 - 1, 1 : 2) + fac101_var_2199 * absa_var_149(ind1_var_2185, 1 : 2) + fac001_var_2198 * absa_var_149(ind1_var_2185 + 1, 1 : 2) + fac211_var_2203 * absa_var_149(ind1_var_2185 + 8, 1 : 2) + fac111_var_2202 * absa_var_149(ind1_var_2185 + 9, 1 : 2) + fac011_var_2201 * absa_var_149(ind1_var_2185 + 10, 1 : 2))
      ELSE
        tau_major1_var_2213(1 : ng16) = speccomb1_var_2221 * (fac001_var_2198 * absa_var_149(ind1_var_2185, 1 : 2) + fac101_var_2199 * absa_var_149(ind1_var_2185 + 1, 1 : 2) + fac011_var_2201 * absa_var_149(ind1_var_2185 + 9, 1 : 2) + fac111_var_2202 * absa_var_149(ind1_var_2185 + 10, 1 : 2))
      END IF
      DO ig_var_2183 = 1, 2
        tauself_var_2211 = selffac_var_2174(jl_var_2236, lay_var_2191) * (selfref_var_151(inds_var_2186, ig_var_2183) + selffrac_var_2175(jl_var_2236, lay_var_2191) * (selfref_var_151(inds_var_2186 + 1, ig_var_2183) - selfref_var_151(inds_var_2186, ig_var_2183)))
        taufor_var_2210 = forfac_var_2181(jl_var_2236, lay_var_2191) * (forref_var_152(indf_var_2187, ig_var_2183) + forfrac_var_2182(jl_var_2236, lay_var_2191) * (forref_var_152(indf_var_2187 + 1, ig_var_2183) - forref_var_152(indf_var_2187, ig_var_2183)))
        taug_var_2161(jl_var_2236, 138 + ig_var_2183, lay_var_2191) = tau_major_var_2212(ig_var_2183) + tau_major1_var_2213(ig_var_2183) + tauself_var_2211 + taufor_var_2210
        fracs_var_2177(jl_var_2236, 138 + ig_var_2183, lay_var_2191) = fracrefa_var_147(ig_var_2183, jpl_var_2190) + fpl_var_2222 * (fracrefa_var_147(ig_var_2183, jpl_var_2190 + 1) - fracrefa_var_147(ig_var_2183, jpl_var_2190))
      END DO
    END DO
  END DO
  DO lay_var_2191 = laytrop_max_var_2227 + 1, klev_var_2160
    DO jl_var_2236 = kidia_var_2158, kfdia_var_2159
      ind0_var_2184 = ((jp_var_2167(jl_var_2236, lay_var_2191) - 13) * 5 + (jt_var_2168(jl_var_2236, lay_var_2191) - 1)) * nspb_var_217(16) + 1
      ind1_var_2185 = ((jp_var_2167(jl_var_2236, lay_var_2191) - 12) * 5 + (jt1_var_2169(jl_var_2236, lay_var_2191) - 1)) * nspb_var_217(16) + 1
      DO ig_var_2183 = 1, 2
        taug_var_2161(jl_var_2236, 138 + ig_var_2183, lay_var_2191) = colch4_var_2172(jl_var_2236, lay_var_2191) * (fac00_var_2163(jl_var_2236, lay_var_2191) * absb_var_150(ind0_var_2184, ig_var_2183) + fac10_var_2165(jl_var_2236, lay_var_2191) * absb_var_150(ind0_var_2184 + 1, ig_var_2183) + fac01_var_2164(jl_var_2236, lay_var_2191) * absb_var_150(ind1_var_2185, ig_var_2183) + fac11_var_2166(jl_var_2236, lay_var_2191) * absb_var_150(ind1_var_2185 + 1, ig_var_2183))
        fracs_var_2177(jl_var_2236, 138 + ig_var_2183, lay_var_2191) = fracrefb_var_148(ig_var_2183)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2227 /= laytrop_min_var_2226) THEN
    DO lay_var_2191 = laytrop_min_var_2226 + 1, laytrop_max_var_2227
      ixc0_var_2233 = ixc_var_2228(lay_var_2191)
      DO ixp_var_2234 = 1, ixc0_var_2233
        jl_var_2236 = ixlow_var_2229(ixp_var_2234, lay_var_2191)
        speccomb_var_2217 = colh2o_var_2171(jl_var_2236, lay_var_2191) + rat_h2och4_var_2178(jl_var_2236, lay_var_2191) * colch4_var_2172(jl_var_2236, lay_var_2191)
        specparm_var_2216 = MIN(colh2o_var_2171(jl_var_2236, lay_var_2191) / speccomb_var_2217, oneminus_var_2170)
        specmult_var_2215 = 8.0D0 * (specparm_var_2216)
        js_var_2188 = 1 + INT(specmult_var_2215)
        fs_var_2214 = ((specmult_var_2215) - AINT((specmult_var_2215)))
        speccomb1_var_2221 = colh2o_var_2171(jl_var_2236, lay_var_2191) + rat_h2och4_1_var_2179(jl_var_2236, lay_var_2191) * colch4_var_2172(jl_var_2236, lay_var_2191)
        specparm1_var_2220 = MIN(colh2o_var_2171(jl_var_2236, lay_var_2191) / speccomb1_var_2221, oneminus_var_2170)
        specmult1_var_2219 = 8.0D0 * (specparm1_var_2220)
        js1_var_2189 = 1 + INT(specmult1_var_2219)
        fs1_var_2218 = ((specmult1_var_2219) - AINT((specmult1_var_2219)))
        speccomb_planck_var_2225 = colh2o_var_2171(jl_var_2236, lay_var_2191) + refrat_planck_a_var_2209 * colch4_var_2172(jl_var_2236, lay_var_2191)
        specparm_planck_var_2224 = MIN(colh2o_var_2171(jl_var_2236, lay_var_2191) / speccomb_planck_var_2225, oneminus_var_2170)
        specmult_planck_var_2223 = 8.0D0 * specparm_planck_var_2224
        jpl_var_2190 = 1 + INT(specmult_planck_var_2223)
        fpl_var_2222 = ((specmult_planck_var_2223) - AINT((specmult_planck_var_2223)))
        ind0_var_2184 = ((jp_var_2167(jl_var_2236, lay_var_2191) - 1) * 5 + (jt_var_2168(jl_var_2236, lay_var_2191) - 1)) * nspa_var_216(16) + js_var_2188
        ind1_var_2185 = (jp_var_2167(jl_var_2236, lay_var_2191) * 5 + (jt1_var_2169(jl_var_2236, lay_var_2191) - 1)) * nspa_var_216(16) + js1_var_2189
        inds_var_2186 = indself_var_2176(jl_var_2236, lay_var_2191)
        indf_var_2187 = indfor_var_2180(jl_var_2236, lay_var_2191)
        IF (specparm_var_2216 .LT. 0.125D0) THEN
          p_var_2204 = fs_var_2214 - 1.0D0
          p4_var_2205 = p_var_2204 ** 4
          fk0_var_2206 = p4_var_2205
          fk1_var_2207 = 1.0D0 - p_var_2204 - 2.0D0 * p4_var_2205
          fk2_var_2208 = p_var_2204 + p4_var_2205
          fac000_var_2192 = fk0_var_2206 * fac00_var_2163(jl_var_2236, lay_var_2191)
          fac100_var_2193 = fk1_var_2207 * fac00_var_2163(jl_var_2236, lay_var_2191)
          fac200_var_2194 = fk2_var_2208 * fac00_var_2163(jl_var_2236, lay_var_2191)
          fac010_var_2195 = fk0_var_2206 * fac10_var_2165(jl_var_2236, lay_var_2191)
          fac110_var_2196 = fk1_var_2207 * fac10_var_2165(jl_var_2236, lay_var_2191)
          fac210_var_2197 = fk2_var_2208 * fac10_var_2165(jl_var_2236, lay_var_2191)
        ELSE IF (specparm_var_2216 .GT. 0.875D0) THEN
          p_var_2204 = - fs_var_2214
          p4_var_2205 = p_var_2204 ** 4
          fk0_var_2206 = p4_var_2205
          fk1_var_2207 = 1.0D0 - p_var_2204 - 2.0D0 * p4_var_2205
          fk2_var_2208 = p_var_2204 + p4_var_2205
          fac000_var_2192 = fk0_var_2206 * fac00_var_2163(jl_var_2236, lay_var_2191)
          fac100_var_2193 = fk1_var_2207 * fac00_var_2163(jl_var_2236, lay_var_2191)
          fac200_var_2194 = fk2_var_2208 * fac00_var_2163(jl_var_2236, lay_var_2191)
          fac010_var_2195 = fk0_var_2206 * fac10_var_2165(jl_var_2236, lay_var_2191)
          fac110_var_2196 = fk1_var_2207 * fac10_var_2165(jl_var_2236, lay_var_2191)
          fac210_var_2197 = fk2_var_2208 * fac10_var_2165(jl_var_2236, lay_var_2191)
        ELSE
          fac000_var_2192 = (1.0D0 - fs_var_2214) * fac00_var_2163(jl_var_2236, lay_var_2191)
          fac010_var_2195 = (1.0D0 - fs_var_2214) * fac10_var_2165(jl_var_2236, lay_var_2191)
          fac100_var_2193 = fs_var_2214 * fac00_var_2163(jl_var_2236, lay_var_2191)
          fac110_var_2196 = fs_var_2214 * fac10_var_2165(jl_var_2236, lay_var_2191)
          fac200_var_2194 = 0.0D0
          fac210_var_2197 = 0.0D0
        END IF
        IF (specparm1_var_2220 .LT. 0.125D0) THEN
          p_var_2204 = fs1_var_2218 - 1.0D0
          p4_var_2205 = p_var_2204 ** 4
          fk0_var_2206 = p4_var_2205
          fk1_var_2207 = 1.0D0 - p_var_2204 - 2.0D0 * p4_var_2205
          fk2_var_2208 = p_var_2204 + p4_var_2205
          fac001_var_2198 = fk0_var_2206 * fac01_var_2164(jl_var_2236, lay_var_2191)
          fac101_var_2199 = fk1_var_2207 * fac01_var_2164(jl_var_2236, lay_var_2191)
          fac201_var_2200 = fk2_var_2208 * fac01_var_2164(jl_var_2236, lay_var_2191)
          fac011_var_2201 = fk0_var_2206 * fac11_var_2166(jl_var_2236, lay_var_2191)
          fac111_var_2202 = fk1_var_2207 * fac11_var_2166(jl_var_2236, lay_var_2191)
          fac211_var_2203 = fk2_var_2208 * fac11_var_2166(jl_var_2236, lay_var_2191)
        ELSE IF (specparm1_var_2220 .GT. 0.875D0) THEN
          p_var_2204 = - fs1_var_2218
          p4_var_2205 = p_var_2204 ** 4
          fk0_var_2206 = p4_var_2205
          fk1_var_2207 = 1.0D0 - p_var_2204 - 2.0D0 * p4_var_2205
          fk2_var_2208 = p_var_2204 + p4_var_2205
          fac001_var_2198 = fk0_var_2206 * fac01_var_2164(jl_var_2236, lay_var_2191)
          fac101_var_2199 = fk1_var_2207 * fac01_var_2164(jl_var_2236, lay_var_2191)
          fac201_var_2200 = fk2_var_2208 * fac01_var_2164(jl_var_2236, lay_var_2191)
          fac011_var_2201 = fk0_var_2206 * fac11_var_2166(jl_var_2236, lay_var_2191)
          fac111_var_2202 = fk1_var_2207 * fac11_var_2166(jl_var_2236, lay_var_2191)
          fac211_var_2203 = fk2_var_2208 * fac11_var_2166(jl_var_2236, lay_var_2191)
        ELSE
          fac001_var_2198 = (1.0D0 - fs1_var_2218) * fac01_var_2164(jl_var_2236, lay_var_2191)
          fac011_var_2201 = (1.0D0 - fs1_var_2218) * fac11_var_2166(jl_var_2236, lay_var_2191)
          fac101_var_2199 = fs1_var_2218 * fac01_var_2164(jl_var_2236, lay_var_2191)
          fac111_var_2202 = fs1_var_2218 * fac11_var_2166(jl_var_2236, lay_var_2191)
          fac201_var_2200 = 0.0D0
          fac211_var_2203 = 0.0D0
        END IF
        IF (specparm_var_2216 .LT. 0.125D0) THEN
          tau_major_var_2212(1 : ng16) = speccomb_var_2217 * (fac000_var_2192 * absa_var_149(ind0_var_2184, 1 : 2) + fac100_var_2193 * absa_var_149(ind0_var_2184 + 1, 1 : 2) + fac200_var_2194 * absa_var_149(ind0_var_2184 + 2, 1 : 2) + fac010_var_2195 * absa_var_149(ind0_var_2184 + 9, 1 : 2) + fac110_var_2196 * absa_var_149(ind0_var_2184 + 10, 1 : 2) + fac210_var_2197 * absa_var_149(ind0_var_2184 + 11, 1 : 2))
        ELSE IF (specparm_var_2216 .GT. 0.875D0) THEN
          tau_major_var_2212(1 : ng16) = speccomb_var_2217 * (fac200_var_2194 * absa_var_149(ind0_var_2184 - 1, 1 : 2) + fac100_var_2193 * absa_var_149(ind0_var_2184, 1 : 2) + fac000_var_2192 * absa_var_149(ind0_var_2184 + 1, 1 : 2) + fac210_var_2197 * absa_var_149(ind0_var_2184 + 8, 1 : 2) + fac110_var_2196 * absa_var_149(ind0_var_2184 + 9, 1 : 2) + fac010_var_2195 * absa_var_149(ind0_var_2184 + 10, 1 : 2))
        ELSE
          tau_major_var_2212(1 : ng16) = speccomb_var_2217 * (fac000_var_2192 * absa_var_149(ind0_var_2184, 1 : 2) + fac100_var_2193 * absa_var_149(ind0_var_2184 + 1, 1 : 2) + fac010_var_2195 * absa_var_149(ind0_var_2184 + 9, 1 : 2) + fac110_var_2196 * absa_var_149(ind0_var_2184 + 10, 1 : 2))
        END IF
        IF (specparm1_var_2220 .LT. 0.125D0) THEN
          tau_major1_var_2213(1 : ng16) = speccomb1_var_2221 * (fac001_var_2198 * absa_var_149(ind1_var_2185, 1 : 2) + fac101_var_2199 * absa_var_149(ind1_var_2185 + 1, 1 : 2) + fac201_var_2200 * absa_var_149(ind1_var_2185 + 2, 1 : 2) + fac011_var_2201 * absa_var_149(ind1_var_2185 + 9, 1 : 2) + fac111_var_2202 * absa_var_149(ind1_var_2185 + 10, 1 : 2) + fac211_var_2203 * absa_var_149(ind1_var_2185 + 11, 1 : 2))
        ELSE IF (specparm1_var_2220 .GT. 0.875D0) THEN
          tau_major1_var_2213(1 : ng16) = speccomb1_var_2221 * (fac201_var_2200 * absa_var_149(ind1_var_2185 - 1, 1 : 2) + fac101_var_2199 * absa_var_149(ind1_var_2185, 1 : 2) + fac001_var_2198 * absa_var_149(ind1_var_2185 + 1, 1 : 2) + fac211_var_2203 * absa_var_149(ind1_var_2185 + 8, 1 : 2) + fac111_var_2202 * absa_var_149(ind1_var_2185 + 9, 1 : 2) + fac011_var_2201 * absa_var_149(ind1_var_2185 + 10, 1 : 2))
        ELSE
          tau_major1_var_2213(1 : ng16) = speccomb1_var_2221 * (fac001_var_2198 * absa_var_149(ind1_var_2185, 1 : 2) + fac101_var_2199 * absa_var_149(ind1_var_2185 + 1, 1 : 2) + fac011_var_2201 * absa_var_149(ind1_var_2185 + 9, 1 : 2) + fac111_var_2202 * absa_var_149(ind1_var_2185 + 10, 1 : 2))
        END IF
        DO ig_var_2183 = 1, 2
          tauself_var_2211 = selffac_var_2174(jl_var_2236, lay_var_2191) * (selfref_var_151(inds_var_2186, ig_var_2183) + selffrac_var_2175(jl_var_2236, lay_var_2191) * (selfref_var_151(inds_var_2186 + 1, ig_var_2183) - selfref_var_151(inds_var_2186, ig_var_2183)))
          taufor_var_2210 = forfac_var_2181(jl_var_2236, lay_var_2191) * (forref_var_152(indf_var_2187, ig_var_2183) + forfrac_var_2182(jl_var_2236, lay_var_2191) * (forref_var_152(indf_var_2187 + 1, ig_var_2183) - forref_var_152(indf_var_2187, ig_var_2183)))
          taug_var_2161(jl_var_2236, 138 + ig_var_2183, lay_var_2191) = tau_major_var_2212(ig_var_2183) + tau_major1_var_2213(ig_var_2183) + tauself_var_2211 + taufor_var_2210
          fracs_var_2177(jl_var_2236, 138 + ig_var_2183, lay_var_2191) = fracrefa_var_147(ig_var_2183, jpl_var_2190) + fpl_var_2222 * (fracrefa_var_147(ig_var_2183, jpl_var_2190 + 1) - fracrefa_var_147(ig_var_2183, jpl_var_2190))
        END DO
      END DO
      ixc0_var_2233 = kfdia_var_2159 - kidia_var_2158 + 1 - ixc0_var_2233
      DO ixp_var_2234 = 1, ixc0_var_2233
        jl_var_2236 = ixhigh_var_2230(ixp_var_2234, lay_var_2191)
        ind0_var_2184 = ((jp_var_2167(jl_var_2236, lay_var_2191) - 13) * 5 + (jt_var_2168(jl_var_2236, lay_var_2191) - 1)) * nspb_var_217(16) + 1
        ind1_var_2185 = ((jp_var_2167(jl_var_2236, lay_var_2191) - 12) * 5 + (jt1_var_2169(jl_var_2236, lay_var_2191) - 1)) * nspb_var_217(16) + 1
        DO ig_var_2183 = 1, 2
          taug_var_2161(jl_var_2236, 138 + ig_var_2183, lay_var_2191) = colch4_var_2172(jl_var_2236, lay_var_2191) * (fac00_var_2163(jl_var_2236, lay_var_2191) * absb_var_150(ind0_var_2184, ig_var_2183) + fac10_var_2165(jl_var_2236, lay_var_2191) * absb_var_150(ind0_var_2184 + 1, ig_var_2183) + fac01_var_2164(jl_var_2236, lay_var_2191) * absb_var_150(ind1_var_2185, ig_var_2183) + fac11_var_2166(jl_var_2236, lay_var_2191) * absb_var_150(ind1_var_2185 + 1, ig_var_2183))
          fracs_var_2177(jl_var_2236, 138 + ig_var_2183, lay_var_2191) = fracrefb_var_148(ig_var_2183)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol16
SUBROUTINE rrtm_taumol12(kidia_var_2237, kfdia_var_2238, klev_var_2239, taug_var_2240, p_tauaerl_var_2241, fac00_var_2242, fac01_var_2243, fac10_var_2244, fac11_var_2245, forfac_var_2261, forfrac_var_2260, indfor_var_2259, jp_var_2246, jt_var_2247, jt1_var_2248, oneminus_var_2249, colh2o_var_2250, colco2_var_2251, laytrop_var_2252, selffac_var_2253, selffrac_var_2254, indself_var_2255, fracs_var_2256, rat_h2oco2_var_2257, rat_h2oco2_1_var_2258)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrtm, ONLY: ng12
  USE yoerrta12, ONLY: absa_var_127, forref_var_129, fracrefa_var_126, selfref_var_128
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2237
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2238
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2239
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2240(kidia_var_2237 : kfdia_var_2238, 140, klev_var_2239)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2241(kidia_var_2237 : kfdia_var_2238, klev_var_2239, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2242(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2243(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2244(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2245(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2246(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2247(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2248(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2249
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2250(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2251(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2252(kidia_var_2237 : kfdia_var_2238)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2253(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2254(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2255(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2256(kidia_var_2237 : kfdia_var_2238, 140, klev_var_2239)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_2257(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_2258(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2259(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2260(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2261(kidia_var_2237 : kfdia_var_2238, klev_var_2239)
  REAL(KIND = 8) :: speccomb_var_2262, speccomb1_var_2263, speccomb_planck_var_2264
  INTEGER(KIND = 4) :: ind0_var_2265, ind1_var_2266, inds_var_2267, indf_var_2268
  INTEGER(KIND = 4) :: ig_var_2269, js_var_2270, lay_var_2271, js1_var_2272, jpl_var_2273
  REAL(KIND = 8) :: fs_var_2274, specmult_var_2275, specparm_var_2276, fs1_var_2277, specmult1_var_2278, specparm1_var_2279, fpl_var_2280, specmult_planck_var_2281, specparm_planck_var_2282
  REAL(KIND = 8) :: fac000_var_2283, fac100_var_2284, fac200_var_2285, fac010_var_2286, fac110_var_2287, fac210_var_2288, fac001_var_2289, fac101_var_2290, fac201_var_2291, fac011_var_2292, fac111_var_2293, fac211_var_2294
  REAL(KIND = 8) :: p_var_2295, p4_var_2296, fk0_var_2297, fk1_var_2298, fk2_var_2299
  REAL(KIND = 8) :: taufor_var_2300, tauself_var_2301, tau_major_var_2302(8), tau_major1_var_2303(8)
  REAL(KIND = 8) :: refrat_planck_a_var_2304
  INTEGER(KIND = 4) :: laytrop_min_var_2305, laytrop_max_var_2306
  INTEGER(KIND = 4) :: ixc_var_2307(klev_var_2239), ixlow_var_2308(kfdia_var_2238, klev_var_2239), ixhigh_var_2309(kfdia_var_2238, klev_var_2239)
  INTEGER(KIND = 4) :: ich_var_2310, icl_var_2311, ixc0_var_2312, ixp_var_2313, jc_var_2314, jl_var_2315
  laytrop_min_var_2305 = MINVAL(laytrop_var_2252)
  laytrop_max_var_2306 = MAXVAL(laytrop_var_2252)
  ixlow_var_2308 = 0
  ixhigh_var_2309 = 0
  ixc_var_2307 = 0
  DO lay_var_2271 = laytrop_min_var_2305 + 1, laytrop_max_var_2306
    icl_var_2311 = 0
    ich_var_2310 = 0
    DO jc_var_2314 = kidia_var_2237, kfdia_var_2238
      IF (lay_var_2271 <= laytrop_var_2252(jc_var_2314)) THEN
        icl_var_2311 = icl_var_2311 + 1
        ixlow_var_2308(icl_var_2311, lay_var_2271) = jc_var_2314
      ELSE
        ich_var_2310 = ich_var_2310 + 1
        ixhigh_var_2309(ich_var_2310, lay_var_2271) = jc_var_2314
      END IF
    END DO
    ixc_var_2307(lay_var_2271) = icl_var_2311
  END DO
  refrat_planck_a_var_2304 = chi_mls(1, 10) / chi_mls(2, 10)
  DO lay_var_2271 = 1, laytrop_min_var_2305
    DO jl_var_2315 = kidia_var_2237, kfdia_var_2238
      speccomb_var_2262 = colh2o_var_2250(jl_var_2315, lay_var_2271) + rat_h2oco2_var_2257(jl_var_2315, lay_var_2271) * colco2_var_2251(jl_var_2315, lay_var_2271)
      specparm_var_2276 = MIN(colh2o_var_2250(jl_var_2315, lay_var_2271) / speccomb_var_2262, oneminus_var_2249)
      specmult_var_2275 = 8.0D0 * (specparm_var_2276)
      js_var_2270 = 1 + INT(specmult_var_2275)
      fs_var_2274 = ((specmult_var_2275) - AINT((specmult_var_2275)))
      speccomb1_var_2263 = colh2o_var_2250(jl_var_2315, lay_var_2271) + rat_h2oco2_1_var_2258(jl_var_2315, lay_var_2271) * colco2_var_2251(jl_var_2315, lay_var_2271)
      specparm1_var_2279 = MIN(colh2o_var_2250(jl_var_2315, lay_var_2271) / speccomb1_var_2263, oneminus_var_2249)
      specmult1_var_2278 = 8.0D0 * (specparm1_var_2279)
      js1_var_2272 = 1 + INT(specmult1_var_2278)
      fs1_var_2277 = ((specmult1_var_2278) - AINT((specmult1_var_2278)))
      speccomb_planck_var_2264 = colh2o_var_2250(jl_var_2315, lay_var_2271) + refrat_planck_a_var_2304 * colco2_var_2251(jl_var_2315, lay_var_2271)
      specparm_planck_var_2282 = MIN(colh2o_var_2250(jl_var_2315, lay_var_2271) / speccomb_planck_var_2264, oneminus_var_2249)
      specmult_planck_var_2281 = 8.0D0 * specparm_planck_var_2282
      jpl_var_2273 = 1 + INT(specmult_planck_var_2281)
      fpl_var_2280 = ((specmult_planck_var_2281) - AINT((specmult_planck_var_2281)))
      ind0_var_2265 = ((jp_var_2246(jl_var_2315, lay_var_2271) - 1) * 5 + (jt_var_2247(jl_var_2315, lay_var_2271) - 1)) * nspa_var_216(12) + js_var_2270
      ind1_var_2266 = (jp_var_2246(jl_var_2315, lay_var_2271) * 5 + (jt1_var_2248(jl_var_2315, lay_var_2271) - 1)) * nspa_var_216(12) + js1_var_2272
      inds_var_2267 = indself_var_2255(jl_var_2315, lay_var_2271)
      indf_var_2268 = indfor_var_2259(jl_var_2315, lay_var_2271)
      IF (specparm_var_2276 .LT. 0.125D0) THEN
        p_var_2295 = fs_var_2274 - 1.0D0
        p4_var_2296 = p_var_2295 ** 4
        fk0_var_2297 = p4_var_2296
        fk1_var_2298 = 1.0D0 - p_var_2295 - 2.0D0 * p4_var_2296
        fk2_var_2299 = p_var_2295 + p4_var_2296
        fac000_var_2283 = fk0_var_2297 * fac00_var_2242(jl_var_2315, lay_var_2271)
        fac100_var_2284 = fk1_var_2298 * fac00_var_2242(jl_var_2315, lay_var_2271)
        fac200_var_2285 = fk2_var_2299 * fac00_var_2242(jl_var_2315, lay_var_2271)
        fac010_var_2286 = fk0_var_2297 * fac10_var_2244(jl_var_2315, lay_var_2271)
        fac110_var_2287 = fk1_var_2298 * fac10_var_2244(jl_var_2315, lay_var_2271)
        fac210_var_2288 = fk2_var_2299 * fac10_var_2244(jl_var_2315, lay_var_2271)
      ELSE IF (specparm_var_2276 .GT. 0.875D0) THEN
        p_var_2295 = - fs_var_2274
        p4_var_2296 = p_var_2295 ** 4
        fk0_var_2297 = p4_var_2296
        fk1_var_2298 = 1.0D0 - p_var_2295 - 2.0D0 * p4_var_2296
        fk2_var_2299 = p_var_2295 + p4_var_2296
        fac000_var_2283 = fk0_var_2297 * fac00_var_2242(jl_var_2315, lay_var_2271)
        fac100_var_2284 = fk1_var_2298 * fac00_var_2242(jl_var_2315, lay_var_2271)
        fac200_var_2285 = fk2_var_2299 * fac00_var_2242(jl_var_2315, lay_var_2271)
        fac010_var_2286 = fk0_var_2297 * fac10_var_2244(jl_var_2315, lay_var_2271)
        fac110_var_2287 = fk1_var_2298 * fac10_var_2244(jl_var_2315, lay_var_2271)
        fac210_var_2288 = fk2_var_2299 * fac10_var_2244(jl_var_2315, lay_var_2271)
      ELSE
        fac000_var_2283 = (1.0D0 - fs_var_2274) * fac00_var_2242(jl_var_2315, lay_var_2271)
        fac010_var_2286 = (1.0D0 - fs_var_2274) * fac10_var_2244(jl_var_2315, lay_var_2271)
        fac100_var_2284 = fs_var_2274 * fac00_var_2242(jl_var_2315, lay_var_2271)
        fac110_var_2287 = fs_var_2274 * fac10_var_2244(jl_var_2315, lay_var_2271)
        fac200_var_2285 = 0.0D0
        fac210_var_2288 = 0.0D0
      END IF
      IF (specparm1_var_2279 .LT. 0.125D0) THEN
        p_var_2295 = fs1_var_2277 - 1.0D0
        p4_var_2296 = p_var_2295 ** 4
        fk0_var_2297 = p4_var_2296
        fk1_var_2298 = 1.0D0 - p_var_2295 - 2.0D0 * p4_var_2296
        fk2_var_2299 = p_var_2295 + p4_var_2296
        fac001_var_2289 = fk0_var_2297 * fac01_var_2243(jl_var_2315, lay_var_2271)
        fac101_var_2290 = fk1_var_2298 * fac01_var_2243(jl_var_2315, lay_var_2271)
        fac201_var_2291 = fk2_var_2299 * fac01_var_2243(jl_var_2315, lay_var_2271)
        fac011_var_2292 = fk0_var_2297 * fac11_var_2245(jl_var_2315, lay_var_2271)
        fac111_var_2293 = fk1_var_2298 * fac11_var_2245(jl_var_2315, lay_var_2271)
        fac211_var_2294 = fk2_var_2299 * fac11_var_2245(jl_var_2315, lay_var_2271)
      ELSE IF (specparm1_var_2279 .GT. 0.875D0) THEN
        p_var_2295 = - fs1_var_2277
        p4_var_2296 = p_var_2295 ** 4
        fk0_var_2297 = p4_var_2296
        fk1_var_2298 = 1.0D0 - p_var_2295 - 2.0D0 * p4_var_2296
        fk2_var_2299 = p_var_2295 + p4_var_2296
        fac001_var_2289 = fk0_var_2297 * fac01_var_2243(jl_var_2315, lay_var_2271)
        fac101_var_2290 = fk1_var_2298 * fac01_var_2243(jl_var_2315, lay_var_2271)
        fac201_var_2291 = fk2_var_2299 * fac01_var_2243(jl_var_2315, lay_var_2271)
        fac011_var_2292 = fk0_var_2297 * fac11_var_2245(jl_var_2315, lay_var_2271)
        fac111_var_2293 = fk1_var_2298 * fac11_var_2245(jl_var_2315, lay_var_2271)
        fac211_var_2294 = fk2_var_2299 * fac11_var_2245(jl_var_2315, lay_var_2271)
      ELSE
        fac001_var_2289 = (1.0D0 - fs1_var_2277) * fac01_var_2243(jl_var_2315, lay_var_2271)
        fac011_var_2292 = (1.0D0 - fs1_var_2277) * fac11_var_2245(jl_var_2315, lay_var_2271)
        fac101_var_2290 = fs1_var_2277 * fac01_var_2243(jl_var_2315, lay_var_2271)
        fac111_var_2293 = fs1_var_2277 * fac11_var_2245(jl_var_2315, lay_var_2271)
        fac201_var_2291 = 0.0D0
        fac211_var_2294 = 0.0D0
      END IF
      IF (specparm_var_2276 .LT. 0.125D0) THEN
        tau_major_var_2302(1 : ng12) = speccomb_var_2262 * (fac000_var_2283 * absa_var_127(ind0_var_2265, 1 : 8) + fac100_var_2284 * absa_var_127(ind0_var_2265 + 1, 1 : 8) + fac200_var_2285 * absa_var_127(ind0_var_2265 + 2, 1 : 8) + fac010_var_2286 * absa_var_127(ind0_var_2265 + 9, 1 : 8) + fac110_var_2287 * absa_var_127(ind0_var_2265 + 10, 1 : 8) + fac210_var_2288 * absa_var_127(ind0_var_2265 + 11, 1 : 8))
      ELSE IF (specparm_var_2276 .GT. 0.875D0) THEN
        tau_major_var_2302(1 : ng12) = speccomb_var_2262 * (fac200_var_2285 * absa_var_127(ind0_var_2265 - 1, 1 : 8) + fac100_var_2284 * absa_var_127(ind0_var_2265, 1 : 8) + fac000_var_2283 * absa_var_127(ind0_var_2265 + 1, 1 : 8) + fac210_var_2288 * absa_var_127(ind0_var_2265 + 8, 1 : 8) + fac110_var_2287 * absa_var_127(ind0_var_2265 + 9, 1 : 8) + fac010_var_2286 * absa_var_127(ind0_var_2265 + 10, 1 : 8))
      ELSE
        tau_major_var_2302(1 : ng12) = speccomb_var_2262 * (fac000_var_2283 * absa_var_127(ind0_var_2265, 1 : 8) + fac100_var_2284 * absa_var_127(ind0_var_2265 + 1, 1 : 8) + fac010_var_2286 * absa_var_127(ind0_var_2265 + 9, 1 : 8) + fac110_var_2287 * absa_var_127(ind0_var_2265 + 10, 1 : 8))
      END IF
      IF (specparm1_var_2279 .LT. 0.125D0) THEN
        tau_major1_var_2303(1 : ng12) = speccomb1_var_2263 * (fac001_var_2289 * absa_var_127(ind1_var_2266, 1 : 8) + fac101_var_2290 * absa_var_127(ind1_var_2266 + 1, 1 : 8) + fac201_var_2291 * absa_var_127(ind1_var_2266 + 2, 1 : 8) + fac011_var_2292 * absa_var_127(ind1_var_2266 + 9, 1 : 8) + fac111_var_2293 * absa_var_127(ind1_var_2266 + 10, 1 : 8) + fac211_var_2294 * absa_var_127(ind1_var_2266 + 11, 1 : 8))
      ELSE IF (specparm1_var_2279 .GT. 0.875D0) THEN
        tau_major1_var_2303(1 : ng12) = speccomb1_var_2263 * (fac201_var_2291 * absa_var_127(ind1_var_2266 - 1, 1 : 8) + fac101_var_2290 * absa_var_127(ind1_var_2266, 1 : 8) + fac001_var_2289 * absa_var_127(ind1_var_2266 + 1, 1 : 8) + fac211_var_2294 * absa_var_127(ind1_var_2266 + 8, 1 : 8) + fac111_var_2293 * absa_var_127(ind1_var_2266 + 9, 1 : 8) + fac011_var_2292 * absa_var_127(ind1_var_2266 + 10, 1 : 8))
      ELSE
        tau_major1_var_2303(1 : ng12) = speccomb1_var_2263 * (fac001_var_2289 * absa_var_127(ind1_var_2266, 1 : 8) + fac101_var_2290 * absa_var_127(ind1_var_2266 + 1, 1 : 8) + fac011_var_2292 * absa_var_127(ind1_var_2266 + 9, 1 : 8) + fac111_var_2293 * absa_var_127(ind1_var_2266 + 10, 1 : 8))
      END IF
      DO ig_var_2269 = 1, 8
        tauself_var_2301 = selffac_var_2253(jl_var_2315, lay_var_2271) * (selfref_var_128(inds_var_2267, ig_var_2269) + selffrac_var_2254(jl_var_2315, lay_var_2271) * (selfref_var_128(inds_var_2267 + 1, ig_var_2269) - selfref_var_128(inds_var_2267, ig_var_2269)))
        taufor_var_2300 = forfac_var_2261(jl_var_2315, lay_var_2271) * (forref_var_129(indf_var_2268, ig_var_2269) + forfrac_var_2260(jl_var_2315, lay_var_2271) * (forref_var_129(indf_var_2268 + 1, ig_var_2269) - forref_var_129(indf_var_2268, ig_var_2269)))
        taug_var_2240(jl_var_2315, 122 + ig_var_2269, lay_var_2271) = tau_major_var_2302(ig_var_2269) + tau_major1_var_2303(ig_var_2269) + tauself_var_2301 + taufor_var_2300
        fracs_var_2256(jl_var_2315, 122 + ig_var_2269, lay_var_2271) = fracrefa_var_126(ig_var_2269, jpl_var_2273) + fpl_var_2280 * (fracrefa_var_126(ig_var_2269, jpl_var_2273 + 1) - fracrefa_var_126(ig_var_2269, jpl_var_2273))
      END DO
    END DO
  END DO
  DO ig_var_2269 = 1, 8
    DO lay_var_2271 = laytrop_max_var_2306 + 1, klev_var_2239
      DO jl_var_2315 = kidia_var_2237, kfdia_var_2238
        taug_var_2240(jl_var_2315, 122 + ig_var_2269, lay_var_2271) = 0.0D0
        fracs_var_2256(jl_var_2315, 122 + ig_var_2269, lay_var_2271) = 0.0D0
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2306 /= laytrop_min_var_2305) THEN
    DO lay_var_2271 = laytrop_min_var_2305 + 1, laytrop_max_var_2306
      ixc0_var_2312 = ixc_var_2307(lay_var_2271)
      DO ixp_var_2313 = 1, ixc0_var_2312
        jl_var_2315 = ixlow_var_2308(ixp_var_2313, lay_var_2271)
        speccomb_var_2262 = colh2o_var_2250(jl_var_2315, lay_var_2271) + rat_h2oco2_var_2257(jl_var_2315, lay_var_2271) * colco2_var_2251(jl_var_2315, lay_var_2271)
        specparm_var_2276 = MIN(colh2o_var_2250(jl_var_2315, lay_var_2271) / speccomb_var_2262, oneminus_var_2249)
        specmult_var_2275 = 8.0D0 * (specparm_var_2276)
        js_var_2270 = 1 + INT(specmult_var_2275)
        fs_var_2274 = ((specmult_var_2275) - AINT((specmult_var_2275)))
        speccomb1_var_2263 = colh2o_var_2250(jl_var_2315, lay_var_2271) + rat_h2oco2_1_var_2258(jl_var_2315, lay_var_2271) * colco2_var_2251(jl_var_2315, lay_var_2271)
        specparm1_var_2279 = MIN(colh2o_var_2250(jl_var_2315, lay_var_2271) / speccomb1_var_2263, oneminus_var_2249)
        specmult1_var_2278 = 8.0D0 * (specparm1_var_2279)
        js1_var_2272 = 1 + INT(specmult1_var_2278)
        fs1_var_2277 = ((specmult1_var_2278) - AINT((specmult1_var_2278)))
        speccomb_planck_var_2264 = colh2o_var_2250(jl_var_2315, lay_var_2271) + refrat_planck_a_var_2304 * colco2_var_2251(jl_var_2315, lay_var_2271)
        specparm_planck_var_2282 = MIN(colh2o_var_2250(jl_var_2315, lay_var_2271) / speccomb_planck_var_2264, oneminus_var_2249)
        specmult_planck_var_2281 = 8.0D0 * specparm_planck_var_2282
        jpl_var_2273 = 1 + INT(specmult_planck_var_2281)
        fpl_var_2280 = ((specmult_planck_var_2281) - AINT((specmult_planck_var_2281)))
        ind0_var_2265 = ((jp_var_2246(jl_var_2315, lay_var_2271) - 1) * 5 + (jt_var_2247(jl_var_2315, lay_var_2271) - 1)) * nspa_var_216(12) + js_var_2270
        ind1_var_2266 = (jp_var_2246(jl_var_2315, lay_var_2271) * 5 + (jt1_var_2248(jl_var_2315, lay_var_2271) - 1)) * nspa_var_216(12) + js1_var_2272
        inds_var_2267 = indself_var_2255(jl_var_2315, lay_var_2271)
        indf_var_2268 = indfor_var_2259(jl_var_2315, lay_var_2271)
        IF (specparm_var_2276 .LT. 0.125D0) THEN
          p_var_2295 = fs_var_2274 - 1.0D0
          p4_var_2296 = p_var_2295 ** 4
          fk0_var_2297 = p4_var_2296
          fk1_var_2298 = 1.0D0 - p_var_2295 - 2.0D0 * p4_var_2296
          fk2_var_2299 = p_var_2295 + p4_var_2296
          fac000_var_2283 = fk0_var_2297 * fac00_var_2242(jl_var_2315, lay_var_2271)
          fac100_var_2284 = fk1_var_2298 * fac00_var_2242(jl_var_2315, lay_var_2271)
          fac200_var_2285 = fk2_var_2299 * fac00_var_2242(jl_var_2315, lay_var_2271)
          fac010_var_2286 = fk0_var_2297 * fac10_var_2244(jl_var_2315, lay_var_2271)
          fac110_var_2287 = fk1_var_2298 * fac10_var_2244(jl_var_2315, lay_var_2271)
          fac210_var_2288 = fk2_var_2299 * fac10_var_2244(jl_var_2315, lay_var_2271)
        ELSE IF (specparm_var_2276 .GT. 0.875D0) THEN
          p_var_2295 = - fs_var_2274
          p4_var_2296 = p_var_2295 ** 4
          fk0_var_2297 = p4_var_2296
          fk1_var_2298 = 1.0D0 - p_var_2295 - 2.0D0 * p4_var_2296
          fk2_var_2299 = p_var_2295 + p4_var_2296
          fac000_var_2283 = fk0_var_2297 * fac00_var_2242(jl_var_2315, lay_var_2271)
          fac100_var_2284 = fk1_var_2298 * fac00_var_2242(jl_var_2315, lay_var_2271)
          fac200_var_2285 = fk2_var_2299 * fac00_var_2242(jl_var_2315, lay_var_2271)
          fac010_var_2286 = fk0_var_2297 * fac10_var_2244(jl_var_2315, lay_var_2271)
          fac110_var_2287 = fk1_var_2298 * fac10_var_2244(jl_var_2315, lay_var_2271)
          fac210_var_2288 = fk2_var_2299 * fac10_var_2244(jl_var_2315, lay_var_2271)
        ELSE
          fac000_var_2283 = (1.0D0 - fs_var_2274) * fac00_var_2242(jl_var_2315, lay_var_2271)
          fac010_var_2286 = (1.0D0 - fs_var_2274) * fac10_var_2244(jl_var_2315, lay_var_2271)
          fac100_var_2284 = fs_var_2274 * fac00_var_2242(jl_var_2315, lay_var_2271)
          fac110_var_2287 = fs_var_2274 * fac10_var_2244(jl_var_2315, lay_var_2271)
          fac200_var_2285 = 0.0D0
          fac210_var_2288 = 0.0D0
        END IF
        IF (specparm1_var_2279 .LT. 0.125D0) THEN
          p_var_2295 = fs1_var_2277 - 1.0D0
          p4_var_2296 = p_var_2295 ** 4
          fk0_var_2297 = p4_var_2296
          fk1_var_2298 = 1.0D0 - p_var_2295 - 2.0D0 * p4_var_2296
          fk2_var_2299 = p_var_2295 + p4_var_2296
          fac001_var_2289 = fk0_var_2297 * fac01_var_2243(jl_var_2315, lay_var_2271)
          fac101_var_2290 = fk1_var_2298 * fac01_var_2243(jl_var_2315, lay_var_2271)
          fac201_var_2291 = fk2_var_2299 * fac01_var_2243(jl_var_2315, lay_var_2271)
          fac011_var_2292 = fk0_var_2297 * fac11_var_2245(jl_var_2315, lay_var_2271)
          fac111_var_2293 = fk1_var_2298 * fac11_var_2245(jl_var_2315, lay_var_2271)
          fac211_var_2294 = fk2_var_2299 * fac11_var_2245(jl_var_2315, lay_var_2271)
        ELSE IF (specparm1_var_2279 .GT. 0.875D0) THEN
          p_var_2295 = - fs1_var_2277
          p4_var_2296 = p_var_2295 ** 4
          fk0_var_2297 = p4_var_2296
          fk1_var_2298 = 1.0D0 - p_var_2295 - 2.0D0 * p4_var_2296
          fk2_var_2299 = p_var_2295 + p4_var_2296
          fac001_var_2289 = fk0_var_2297 * fac01_var_2243(jl_var_2315, lay_var_2271)
          fac101_var_2290 = fk1_var_2298 * fac01_var_2243(jl_var_2315, lay_var_2271)
          fac201_var_2291 = fk2_var_2299 * fac01_var_2243(jl_var_2315, lay_var_2271)
          fac011_var_2292 = fk0_var_2297 * fac11_var_2245(jl_var_2315, lay_var_2271)
          fac111_var_2293 = fk1_var_2298 * fac11_var_2245(jl_var_2315, lay_var_2271)
          fac211_var_2294 = fk2_var_2299 * fac11_var_2245(jl_var_2315, lay_var_2271)
        ELSE
          fac001_var_2289 = (1.0D0 - fs1_var_2277) * fac01_var_2243(jl_var_2315, lay_var_2271)
          fac011_var_2292 = (1.0D0 - fs1_var_2277) * fac11_var_2245(jl_var_2315, lay_var_2271)
          fac101_var_2290 = fs1_var_2277 * fac01_var_2243(jl_var_2315, lay_var_2271)
          fac111_var_2293 = fs1_var_2277 * fac11_var_2245(jl_var_2315, lay_var_2271)
          fac201_var_2291 = 0.0D0
          fac211_var_2294 = 0.0D0
        END IF
        IF (specparm_var_2276 .LT. 0.125D0) THEN
          tau_major_var_2302(1 : ng12) = speccomb_var_2262 * (fac000_var_2283 * absa_var_127(ind0_var_2265, 1 : 8) + fac100_var_2284 * absa_var_127(ind0_var_2265 + 1, 1 : 8) + fac200_var_2285 * absa_var_127(ind0_var_2265 + 2, 1 : 8) + fac010_var_2286 * absa_var_127(ind0_var_2265 + 9, 1 : 8) + fac110_var_2287 * absa_var_127(ind0_var_2265 + 10, 1 : 8) + fac210_var_2288 * absa_var_127(ind0_var_2265 + 11, 1 : 8))
        ELSE IF (specparm_var_2276 .GT. 0.875D0) THEN
          tau_major_var_2302(1 : ng12) = speccomb_var_2262 * (fac200_var_2285 * absa_var_127(ind0_var_2265 - 1, 1 : 8) + fac100_var_2284 * absa_var_127(ind0_var_2265, 1 : 8) + fac000_var_2283 * absa_var_127(ind0_var_2265 + 1, 1 : 8) + fac210_var_2288 * absa_var_127(ind0_var_2265 + 8, 1 : 8) + fac110_var_2287 * absa_var_127(ind0_var_2265 + 9, 1 : 8) + fac010_var_2286 * absa_var_127(ind0_var_2265 + 10, 1 : 8))
        ELSE
          tau_major_var_2302(1 : ng12) = speccomb_var_2262 * (fac000_var_2283 * absa_var_127(ind0_var_2265, 1 : 8) + fac100_var_2284 * absa_var_127(ind0_var_2265 + 1, 1 : 8) + fac010_var_2286 * absa_var_127(ind0_var_2265 + 9, 1 : 8) + fac110_var_2287 * absa_var_127(ind0_var_2265 + 10, 1 : 8))
        END IF
        IF (specparm1_var_2279 .LT. 0.125D0) THEN
          tau_major1_var_2303(1 : ng12) = speccomb1_var_2263 * (fac001_var_2289 * absa_var_127(ind1_var_2266, 1 : 8) + fac101_var_2290 * absa_var_127(ind1_var_2266 + 1, 1 : 8) + fac201_var_2291 * absa_var_127(ind1_var_2266 + 2, 1 : 8) + fac011_var_2292 * absa_var_127(ind1_var_2266 + 9, 1 : 8) + fac111_var_2293 * absa_var_127(ind1_var_2266 + 10, 1 : 8) + fac211_var_2294 * absa_var_127(ind1_var_2266 + 11, 1 : 8))
        ELSE IF (specparm1_var_2279 .GT. 0.875D0) THEN
          tau_major1_var_2303(1 : ng12) = speccomb1_var_2263 * (fac201_var_2291 * absa_var_127(ind1_var_2266 - 1, 1 : 8) + fac101_var_2290 * absa_var_127(ind1_var_2266, 1 : 8) + fac001_var_2289 * absa_var_127(ind1_var_2266 + 1, 1 : 8) + fac211_var_2294 * absa_var_127(ind1_var_2266 + 8, 1 : 8) + fac111_var_2293 * absa_var_127(ind1_var_2266 + 9, 1 : 8) + fac011_var_2292 * absa_var_127(ind1_var_2266 + 10, 1 : 8))
        ELSE
          tau_major1_var_2303(1 : ng12) = speccomb1_var_2263 * (fac001_var_2289 * absa_var_127(ind1_var_2266, 1 : 8) + fac101_var_2290 * absa_var_127(ind1_var_2266 + 1, 1 : 8) + fac011_var_2292 * absa_var_127(ind1_var_2266 + 9, 1 : 8) + fac111_var_2293 * absa_var_127(ind1_var_2266 + 10, 1 : 8))
        END IF
        DO ig_var_2269 = 1, 8
          tauself_var_2301 = selffac_var_2253(jl_var_2315, lay_var_2271) * (selfref_var_128(inds_var_2267, ig_var_2269) + selffrac_var_2254(jl_var_2315, lay_var_2271) * (selfref_var_128(inds_var_2267 + 1, ig_var_2269) - selfref_var_128(inds_var_2267, ig_var_2269)))
          taufor_var_2300 = forfac_var_2261(jl_var_2315, lay_var_2271) * (forref_var_129(indf_var_2268, ig_var_2269) + forfrac_var_2260(jl_var_2315, lay_var_2271) * (forref_var_129(indf_var_2268 + 1, ig_var_2269) - forref_var_129(indf_var_2268, ig_var_2269)))
          taug_var_2240(jl_var_2315, 122 + ig_var_2269, lay_var_2271) = tau_major_var_2302(ig_var_2269) + tau_major1_var_2303(ig_var_2269) + tauself_var_2301 + taufor_var_2300
          fracs_var_2256(jl_var_2315, 122 + ig_var_2269, lay_var_2271) = fracrefa_var_126(ig_var_2269, jpl_var_2273) + fpl_var_2280 * (fracrefa_var_126(ig_var_2269, jpl_var_2273 + 1) - fracrefa_var_126(ig_var_2269, jpl_var_2273))
        END DO
      END DO
      ixc0_var_2312 = kfdia_var_2238 - kidia_var_2237 + 1 - ixc0_var_2312
      DO ig_var_2269 = 1, 8
        DO ixp_var_2313 = 1, ixc0_var_2312
          jl_var_2315 = ixhigh_var_2309(ixp_var_2313, lay_var_2271)
          taug_var_2240(jl_var_2315, 122 + ig_var_2269, lay_var_2271) = 0.0D0
          fracs_var_2256(jl_var_2315, 122 + ig_var_2269, lay_var_2271) = 0.0D0
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol12
SUBROUTINE rrtm_taumol13(kidia_var_2316, kfdia_var_2317, klev_var_2318, taug_var_2319, p_tauaerl_var_2320, fac00_var_2321, fac01_var_2322, fac10_var_2323, fac11_var_2324, forfac_var_2340, forfrac_var_2341, indfor_var_2339, jp_var_2325, jt_var_2326, jt1_var_2327, oneminus_var_2328, colh2o_var_2329, coln2o_var_2330, colco2_var_2331, colo3_var_2332, coldry_var_2333, laytrop_var_2334, selffac_var_2335, selffrac_var_2336, indself_var_2337, fracs_var_2338, rat_h2on2o, rat_h2on2o_1, minorfrac_var_2342, indminor_var_2343)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrtm, ONLY: ng13
  USE yoerrta13, ONLY: absa_var_132, forref_var_134, fracrefa_var_130, fracrefb_var_131, ka_mco, ka_mco2_var_135, kb_mo3, selfref_var_133
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2316
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2317
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2318
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2319(kidia_var_2316 : kfdia_var_2317, 140, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2320(kidia_var_2316 : kfdia_var_2317, klev_var_2318, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2321(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2322(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2323(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2324(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2325(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2326(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2327(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2328
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2329(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_2330(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2331(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_2332(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_2333(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2334(kidia_var_2316 : kfdia_var_2317)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2335(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2336(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2337(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2338(kidia_var_2316 : kfdia_var_2317, 140, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: rat_h2on2o(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: rat_h2on2o_1(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2339(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2340(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2341(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2342(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2343(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  REAL(KIND = 8) :: speccomb_var_2344, speccomb1_var_2345, speccomb_planck_var_2346, speccomb_mco2_var_2347, speccomb_mco
  INTEGER(KIND = 4) :: ind0_var_2348, ind1_var_2349, inds_var_2350, indf_var_2351, indm_var_2352
  INTEGER(KIND = 4) :: ig_var_2353, js_var_2354, lay_var_2355, js1_var_2356, jpl_var_2357, jmco2_var_2358, jmco
  REAL(KIND = 8) :: refrat_planck_a_var_2359, refrat_m_a_var_2360, refrat_m_a3
  REAL(KIND = 8) :: fac000_var_2361, fac100_var_2362, fac200_var_2363, fac010_var_2364, fac110_var_2365, fac210_var_2366, fac001_var_2367, fac101_var_2368, fac201_var_2369, fac011_var_2370, fac111_var_2371, fac211_var_2372
  REAL(KIND = 8) :: p_var_2373, p4_var_2374, fk0_var_2375, fk1_var_2376, fk2_var_2377
  REAL(KIND = 8) :: taufor_var_2378, tauself_var_2379, tau_major_var_2380(4), tau_major1_var_2381(4), co2m1_var_2382, co2m2_var_2383, absco2_var_2384
  REAL(KIND = 8) :: com1, com2, absco, abso3_var_2385
  REAL(KIND = 8) :: chi_co2_var_2386, ratco2_var_2387, adjfac_var_2388, adjcolco2_var_2389
  REAL(KIND = 8) :: fs_var_2390, specmult_var_2391, specparm_var_2392, fs1_var_2393, specmult1_var_2394, specparm1_var_2395, fmco2_var_2396, specmult_mco2_var_2397, specparm_mco2_var_2398, fmco, specmult_mco, specparm_mco, fpl_var_2399, specmult_planck_var_2400, specparm_planck_var_2401
  REAL(KIND = 8) :: colco(kidia_var_2316 : kfdia_var_2317, klev_var_2318)
  INTEGER(KIND = 4) :: laytrop_min_var_2402, laytrop_max_var_2403
  INTEGER(KIND = 4) :: ixc_var_2404(klev_var_2318), ixlow_var_2405(kfdia_var_2317, klev_var_2318), ixhigh_var_2406(kfdia_var_2317, klev_var_2318)
  INTEGER(KIND = 4) :: ich_var_2407, icl_var_2408, ixc0_var_2409, ixp_var_2410, jc_var_2411, jl_var_2412
  DO lay_var_2355 = 1, klev_var_2318
    DO jc_var_2411 = kidia_var_2316, kfdia_var_2317
      colco(jc_var_2411, lay_var_2355) = 0.0D0
    END DO
  END DO
  laytrop_min_var_2402 = MINVAL(laytrop_var_2334)
  laytrop_max_var_2403 = MAXVAL(laytrop_var_2334)
  ixlow_var_2405 = 0
  ixhigh_var_2406 = 0
  ixc_var_2404 = 0
  DO lay_var_2355 = laytrop_min_var_2402 + 1, laytrop_max_var_2403
    icl_var_2408 = 0
    ich_var_2407 = 0
    DO jc_var_2411 = kidia_var_2316, kfdia_var_2317
      IF (lay_var_2355 <= laytrop_var_2334(jc_var_2411)) THEN
        icl_var_2408 = icl_var_2408 + 1
        ixlow_var_2405(icl_var_2408, lay_var_2355) = jc_var_2411
      ELSE
        ich_var_2407 = ich_var_2407 + 1
        ixhigh_var_2406(ich_var_2407, lay_var_2355) = jc_var_2411
      END IF
    END DO
    ixc_var_2404(lay_var_2355) = icl_var_2408
  END DO
  refrat_planck_a_var_2359 = chi_mls(1, 5) / chi_mls(4, 5)
  refrat_m_a_var_2360 = chi_mls(1, 1) / chi_mls(4, 1)
  refrat_m_a3 = chi_mls(1, 3) / chi_mls(4, 3)
  DO lay_var_2355 = 1, laytrop_min_var_2402
    DO jl_var_2412 = kidia_var_2316, kfdia_var_2317
      speccomb_var_2344 = colh2o_var_2329(jl_var_2412, lay_var_2355) + rat_h2on2o(jl_var_2412, lay_var_2355) * coln2o_var_2330(jl_var_2412, lay_var_2355)
      specparm_var_2392 = MIN(colh2o_var_2329(jl_var_2412, lay_var_2355) / speccomb_var_2344, oneminus_var_2328)
      specmult_var_2391 = 8.0D0 * (specparm_var_2392)
      js_var_2354 = 1 + INT(specmult_var_2391)
      fs_var_2390 = ((specmult_var_2391) - AINT((specmult_var_2391)))
      speccomb1_var_2345 = colh2o_var_2329(jl_var_2412, lay_var_2355) + rat_h2on2o_1(jl_var_2412, lay_var_2355) * coln2o_var_2330(jl_var_2412, lay_var_2355)
      specparm1_var_2395 = MIN(colh2o_var_2329(jl_var_2412, lay_var_2355) / speccomb1_var_2345, oneminus_var_2328)
      specmult1_var_2394 = 8.0D0 * (specparm1_var_2395)
      js1_var_2356 = 1 + INT(specmult1_var_2394)
      fs1_var_2393 = ((specmult1_var_2394) - AINT((specmult1_var_2394)))
      speccomb_mco2_var_2347 = colh2o_var_2329(jl_var_2412, lay_var_2355) + refrat_m_a_var_2360 * coln2o_var_2330(jl_var_2412, lay_var_2355)
      specparm_mco2_var_2398 = MIN(colh2o_var_2329(jl_var_2412, lay_var_2355) / speccomb_mco2_var_2347, oneminus_var_2328)
      specmult_mco2_var_2397 = 8.0D0 * specparm_mco2_var_2398
      jmco2_var_2358 = 1 + INT(specmult_mco2_var_2397)
      fmco2_var_2396 = ((specmult_mco2_var_2397) - AINT((specmult_mco2_var_2397)))
      chi_co2_var_2386 = colco2_var_2331(jl_var_2412, lay_var_2355) / (coldry_var_2333(jl_var_2412, lay_var_2355))
      ratco2_var_2387 = 1D+20 * chi_co2_var_2386 / 0.000355D0
      IF (ratco2_var_2387 .GT. 3.0D0) THEN
        adjfac_var_2388 = 2.0D0 + (ratco2_var_2387 - 2.0D0) ** 0.68D0
        adjcolco2_var_2389 = adjfac_var_2388 * 0.000355D0 * coldry_var_2333(jl_var_2412, lay_var_2355) * 1D-20
      ELSE
        adjcolco2_var_2389 = colco2_var_2331(jl_var_2412, lay_var_2355)
      END IF
      speccomb_mco = colh2o_var_2329(jl_var_2412, lay_var_2355) + refrat_m_a3 * coln2o_var_2330(jl_var_2412, lay_var_2355)
      specparm_mco = MIN(colh2o_var_2329(jl_var_2412, lay_var_2355) / speccomb_mco, oneminus_var_2328)
      specmult_mco = 8.0D0 * specparm_mco
      jmco = 1 + INT(specmult_mco)
      fmco = ((specmult_mco) - AINT((specmult_mco)))
      speccomb_planck_var_2346 = colh2o_var_2329(jl_var_2412, lay_var_2355) + refrat_planck_a_var_2359 * coln2o_var_2330(jl_var_2412, lay_var_2355)
      specparm_planck_var_2401 = MIN(colh2o_var_2329(jl_var_2412, lay_var_2355) / speccomb_planck_var_2346, oneminus_var_2328)
      specmult_planck_var_2400 = 8.0D0 * specparm_planck_var_2401
      jpl_var_2357 = 1 + INT(specmult_planck_var_2400)
      fpl_var_2399 = ((specmult_planck_var_2400) - AINT((specmult_planck_var_2400)))
      ind0_var_2348 = ((jp_var_2325(jl_var_2412, lay_var_2355) - 1) * 5 + (jt_var_2326(jl_var_2412, lay_var_2355) - 1)) * nspa_var_216(13) + js_var_2354
      ind1_var_2349 = (jp_var_2325(jl_var_2412, lay_var_2355) * 5 + (jt1_var_2327(jl_var_2412, lay_var_2355) - 1)) * nspa_var_216(13) + js1_var_2356
      inds_var_2350 = indself_var_2337(jl_var_2412, lay_var_2355)
      indf_var_2351 = indfor_var_2339(jl_var_2412, lay_var_2355)
      indm_var_2352 = indminor_var_2343(jl_var_2412, lay_var_2355)
      IF (specparm_var_2392 .LT. 0.125D0) THEN
        p_var_2373 = fs_var_2390 - 1.0D0
        p4_var_2374 = p_var_2373 ** 4
        fk0_var_2375 = p4_var_2374
        fk1_var_2376 = 1.0D0 - p_var_2373 - 2.0D0 * p4_var_2374
        fk2_var_2377 = p_var_2373 + p4_var_2374
        fac000_var_2361 = fk0_var_2375 * fac00_var_2321(jl_var_2412, lay_var_2355)
        fac100_var_2362 = fk1_var_2376 * fac00_var_2321(jl_var_2412, lay_var_2355)
        fac200_var_2363 = fk2_var_2377 * fac00_var_2321(jl_var_2412, lay_var_2355)
        fac010_var_2364 = fk0_var_2375 * fac10_var_2323(jl_var_2412, lay_var_2355)
        fac110_var_2365 = fk1_var_2376 * fac10_var_2323(jl_var_2412, lay_var_2355)
        fac210_var_2366 = fk2_var_2377 * fac10_var_2323(jl_var_2412, lay_var_2355)
      ELSE IF (specparm_var_2392 .GT. 0.875D0) THEN
        p_var_2373 = - fs_var_2390
        p4_var_2374 = p_var_2373 ** 4
        fk0_var_2375 = p4_var_2374
        fk1_var_2376 = 1.0D0 - p_var_2373 - 2.0D0 * p4_var_2374
        fk2_var_2377 = p_var_2373 + p4_var_2374
        fac000_var_2361 = fk0_var_2375 * fac00_var_2321(jl_var_2412, lay_var_2355)
        fac100_var_2362 = fk1_var_2376 * fac00_var_2321(jl_var_2412, lay_var_2355)
        fac200_var_2363 = fk2_var_2377 * fac00_var_2321(jl_var_2412, lay_var_2355)
        fac010_var_2364 = fk0_var_2375 * fac10_var_2323(jl_var_2412, lay_var_2355)
        fac110_var_2365 = fk1_var_2376 * fac10_var_2323(jl_var_2412, lay_var_2355)
        fac210_var_2366 = fk2_var_2377 * fac10_var_2323(jl_var_2412, lay_var_2355)
      ELSE
        fac000_var_2361 = (1.0D0 - fs_var_2390) * fac00_var_2321(jl_var_2412, lay_var_2355)
        fac010_var_2364 = (1.0D0 - fs_var_2390) * fac10_var_2323(jl_var_2412, lay_var_2355)
        fac100_var_2362 = fs_var_2390 * fac00_var_2321(jl_var_2412, lay_var_2355)
        fac110_var_2365 = fs_var_2390 * fac10_var_2323(jl_var_2412, lay_var_2355)
        fac200_var_2363 = 0.0D0
        fac210_var_2366 = 0.0D0
      END IF
      IF (specparm1_var_2395 .LT. 0.125D0) THEN
        p_var_2373 = fs1_var_2393 - 1.0D0
        p4_var_2374 = p_var_2373 ** 4
        fk0_var_2375 = p4_var_2374
        fk1_var_2376 = 1.0D0 - p_var_2373 - 2.0D0 * p4_var_2374
        fk2_var_2377 = p_var_2373 + p4_var_2374
        fac001_var_2367 = fk0_var_2375 * fac01_var_2322(jl_var_2412, lay_var_2355)
        fac101_var_2368 = fk1_var_2376 * fac01_var_2322(jl_var_2412, lay_var_2355)
        fac201_var_2369 = fk2_var_2377 * fac01_var_2322(jl_var_2412, lay_var_2355)
        fac011_var_2370 = fk0_var_2375 * fac11_var_2324(jl_var_2412, lay_var_2355)
        fac111_var_2371 = fk1_var_2376 * fac11_var_2324(jl_var_2412, lay_var_2355)
        fac211_var_2372 = fk2_var_2377 * fac11_var_2324(jl_var_2412, lay_var_2355)
      ELSE IF (specparm1_var_2395 .GT. 0.875D0) THEN
        p_var_2373 = - fs1_var_2393
        p4_var_2374 = p_var_2373 ** 4
        fk0_var_2375 = p4_var_2374
        fk1_var_2376 = 1.0D0 - p_var_2373 - 2.0D0 * p4_var_2374
        fk2_var_2377 = p_var_2373 + p4_var_2374
        fac001_var_2367 = fk0_var_2375 * fac01_var_2322(jl_var_2412, lay_var_2355)
        fac101_var_2368 = fk1_var_2376 * fac01_var_2322(jl_var_2412, lay_var_2355)
        fac201_var_2369 = fk2_var_2377 * fac01_var_2322(jl_var_2412, lay_var_2355)
        fac011_var_2370 = fk0_var_2375 * fac11_var_2324(jl_var_2412, lay_var_2355)
        fac111_var_2371 = fk1_var_2376 * fac11_var_2324(jl_var_2412, lay_var_2355)
        fac211_var_2372 = fk2_var_2377 * fac11_var_2324(jl_var_2412, lay_var_2355)
      ELSE
        fac001_var_2367 = (1.0D0 - fs1_var_2393) * fac01_var_2322(jl_var_2412, lay_var_2355)
        fac011_var_2370 = (1.0D0 - fs1_var_2393) * fac11_var_2324(jl_var_2412, lay_var_2355)
        fac101_var_2368 = fs1_var_2393 * fac01_var_2322(jl_var_2412, lay_var_2355)
        fac111_var_2371 = fs1_var_2393 * fac11_var_2324(jl_var_2412, lay_var_2355)
        fac201_var_2369 = 0.0D0
        fac211_var_2372 = 0.0D0
      END IF
      IF (specparm_var_2392 .LT. 0.125D0) THEN
        tau_major_var_2380(1 : ng13) = speccomb_var_2344 * (fac000_var_2361 * absa_var_132(ind0_var_2348, 1 : 4) + fac100_var_2362 * absa_var_132(ind0_var_2348 + 1, 1 : 4) + fac200_var_2363 * absa_var_132(ind0_var_2348 + 2, 1 : 4) + fac010_var_2364 * absa_var_132(ind0_var_2348 + 9, 1 : 4) + fac110_var_2365 * absa_var_132(ind0_var_2348 + 10, 1 : 4) + fac210_var_2366 * absa_var_132(ind0_var_2348 + 11, 1 : 4))
      ELSE IF (specparm_var_2392 .GT. 0.875D0) THEN
        tau_major_var_2380(1 : ng13) = speccomb_var_2344 * (fac200_var_2363 * absa_var_132(ind0_var_2348 - 1, 1 : 4) + fac100_var_2362 * absa_var_132(ind0_var_2348, 1 : 4) + fac000_var_2361 * absa_var_132(ind0_var_2348 + 1, 1 : 4) + fac210_var_2366 * absa_var_132(ind0_var_2348 + 8, 1 : 4) + fac110_var_2365 * absa_var_132(ind0_var_2348 + 9, 1 : 4) + fac010_var_2364 * absa_var_132(ind0_var_2348 + 10, 1 : 4))
      ELSE
        tau_major_var_2380(1 : ng13) = speccomb_var_2344 * (fac000_var_2361 * absa_var_132(ind0_var_2348, 1 : 4) + fac100_var_2362 * absa_var_132(ind0_var_2348 + 1, 1 : 4) + fac010_var_2364 * absa_var_132(ind0_var_2348 + 9, 1 : 4) + fac110_var_2365 * absa_var_132(ind0_var_2348 + 10, 1 : 4))
      END IF
      IF (specparm1_var_2395 .LT. 0.125D0) THEN
        tau_major1_var_2381(1 : ng13) = speccomb1_var_2345 * (fac001_var_2367 * absa_var_132(ind1_var_2349, 1 : 4) + fac101_var_2368 * absa_var_132(ind1_var_2349 + 1, 1 : 4) + fac201_var_2369 * absa_var_132(ind1_var_2349 + 2, 1 : 4) + fac011_var_2370 * absa_var_132(ind1_var_2349 + 9, 1 : 4) + fac111_var_2371 * absa_var_132(ind1_var_2349 + 10, 1 : 4) + fac211_var_2372 * absa_var_132(ind1_var_2349 + 11, 1 : 4))
      ELSE IF (specparm1_var_2395 .GT. 0.875D0) THEN
        tau_major1_var_2381(1 : ng13) = speccomb1_var_2345 * (fac201_var_2369 * absa_var_132(ind1_var_2349 - 1, 1 : 4) + fac101_var_2368 * absa_var_132(ind1_var_2349, 1 : 4) + fac001_var_2367 * absa_var_132(ind1_var_2349 + 1, 1 : 4) + fac211_var_2372 * absa_var_132(ind1_var_2349 + 8, 1 : 4) + fac111_var_2371 * absa_var_132(ind1_var_2349 + 9, 1 : 4) + fac011_var_2370 * absa_var_132(ind1_var_2349 + 10, 1 : 4))
      ELSE
        tau_major1_var_2381(1 : ng13) = speccomb1_var_2345 * (fac001_var_2367 * absa_var_132(ind1_var_2349, 1 : 4) + fac101_var_2368 * absa_var_132(ind1_var_2349 + 1, 1 : 4) + fac011_var_2370 * absa_var_132(ind1_var_2349 + 9, 1 : 4) + fac111_var_2371 * absa_var_132(ind1_var_2349 + 10, 1 : 4))
      END IF
      DO ig_var_2353 = 1, 4
        tauself_var_2379 = selffac_var_2335(jl_var_2412, lay_var_2355) * (selfref_var_133(inds_var_2350, ig_var_2353) + selffrac_var_2336(jl_var_2412, lay_var_2355) * (selfref_var_133(inds_var_2350 + 1, ig_var_2353) - selfref_var_133(inds_var_2350, ig_var_2353)))
        taufor_var_2378 = forfac_var_2340(jl_var_2412, lay_var_2355) * (forref_var_134(indf_var_2351, ig_var_2353) + forfrac_var_2341(jl_var_2412, lay_var_2355) * (forref_var_134(indf_var_2351 + 1, ig_var_2353) - forref_var_134(indf_var_2351, ig_var_2353)))
        co2m1_var_2382 = ka_mco2_var_135(jmco2_var_2358, indm_var_2352, ig_var_2353) + fmco2_var_2396 * (ka_mco2_var_135(jmco2_var_2358 + 1, indm_var_2352, ig_var_2353) - ka_mco2_var_135(jmco2_var_2358, indm_var_2352, ig_var_2353))
        co2m2_var_2383 = ka_mco2_var_135(jmco2_var_2358, indm_var_2352 + 1, ig_var_2353) + fmco2_var_2396 * (ka_mco2_var_135(jmco2_var_2358 + 1, indm_var_2352 + 1, ig_var_2353) - ka_mco2_var_135(jmco2_var_2358, indm_var_2352 + 1, ig_var_2353))
        absco2_var_2384 = co2m1_var_2382 + minorfrac_var_2342(jl_var_2412, lay_var_2355) * (co2m2_var_2383 - co2m1_var_2382)
        com1 = ka_mco(jmco, indm_var_2352, ig_var_2353) + fmco * (ka_mco(jmco + 1, indm_var_2352, ig_var_2353) - ka_mco(jmco, indm_var_2352, ig_var_2353))
        com2 = ka_mco(jmco, indm_var_2352 + 1, ig_var_2353) + fmco * (ka_mco(jmco + 1, indm_var_2352 + 1, ig_var_2353) - ka_mco(jmco, indm_var_2352 + 1, ig_var_2353))
        absco = com1 + minorfrac_var_2342(jl_var_2412, lay_var_2355) * (com2 - com1)
        taug_var_2319(jl_var_2412, 130 + ig_var_2353, lay_var_2355) = tau_major_var_2380(ig_var_2353) + tau_major1_var_2381(ig_var_2353) + tauself_var_2379 + taufor_var_2378 + adjcolco2_var_2389 * absco2_var_2384 + colco(jl_var_2412, lay_var_2355) * absco
        fracs_var_2338(jl_var_2412, 130 + ig_var_2353, lay_var_2355) = fracrefa_var_130(ig_var_2353, jpl_var_2357) + fpl_var_2399 * (fracrefa_var_130(ig_var_2353, jpl_var_2357 + 1) - fracrefa_var_130(ig_var_2353, jpl_var_2357))
      END DO
    END DO
  END DO
  DO lay_var_2355 = laytrop_max_var_2403 + 1, klev_var_2318
    DO jl_var_2412 = kidia_var_2316, kfdia_var_2317
      indm_var_2352 = indminor_var_2343(jl_var_2412, lay_var_2355)
      DO ig_var_2353 = 1, 4
        abso3_var_2385 = kb_mo3(indm_var_2352, ig_var_2353) + minorfrac_var_2342(jl_var_2412, lay_var_2355) * (kb_mo3(indm_var_2352 + 1, ig_var_2353) - kb_mo3(indm_var_2352, ig_var_2353))
        taug_var_2319(jl_var_2412, 130 + ig_var_2353, lay_var_2355) = colo3_var_2332(jl_var_2412, lay_var_2355) * abso3_var_2385
        fracs_var_2338(jl_var_2412, 130 + ig_var_2353, lay_var_2355) = fracrefb_var_131(ig_var_2353)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2403 /= laytrop_min_var_2402) THEN
    DO lay_var_2355 = laytrop_min_var_2402 + 1, laytrop_max_var_2403
      ixc0_var_2409 = ixc_var_2404(lay_var_2355)
      DO ixp_var_2410 = 1, ixc0_var_2409
        jl_var_2412 = ixlow_var_2405(ixp_var_2410, lay_var_2355)
        speccomb_var_2344 = colh2o_var_2329(jl_var_2412, lay_var_2355) + rat_h2on2o(jl_var_2412, lay_var_2355) * coln2o_var_2330(jl_var_2412, lay_var_2355)
        specparm_var_2392 = MIN(colh2o_var_2329(jl_var_2412, lay_var_2355) / speccomb_var_2344, oneminus_var_2328)
        specmult_var_2391 = 8.0D0 * (specparm_var_2392)
        js_var_2354 = 1 + INT(specmult_var_2391)
        fs_var_2390 = ((specmult_var_2391) - AINT((specmult_var_2391)))
        speccomb1_var_2345 = colh2o_var_2329(jl_var_2412, lay_var_2355) + rat_h2on2o_1(jl_var_2412, lay_var_2355) * coln2o_var_2330(jl_var_2412, lay_var_2355)
        specparm1_var_2395 = MIN(colh2o_var_2329(jl_var_2412, lay_var_2355) / speccomb1_var_2345, oneminus_var_2328)
        specmult1_var_2394 = 8.0D0 * (specparm1_var_2395)
        js1_var_2356 = 1 + INT(specmult1_var_2394)
        fs1_var_2393 = ((specmult1_var_2394) - AINT((specmult1_var_2394)))
        speccomb_mco2_var_2347 = colh2o_var_2329(jl_var_2412, lay_var_2355) + refrat_m_a_var_2360 * coln2o_var_2330(jl_var_2412, lay_var_2355)
        specparm_mco2_var_2398 = MIN(colh2o_var_2329(jl_var_2412, lay_var_2355) / speccomb_mco2_var_2347, oneminus_var_2328)
        specmult_mco2_var_2397 = 8.0D0 * specparm_mco2_var_2398
        jmco2_var_2358 = 1 + INT(specmult_mco2_var_2397)
        fmco2_var_2396 = ((specmult_mco2_var_2397) - AINT((specmult_mco2_var_2397)))
        chi_co2_var_2386 = colco2_var_2331(jl_var_2412, lay_var_2355) / (coldry_var_2333(jl_var_2412, lay_var_2355))
        ratco2_var_2387 = 1D+20 * chi_co2_var_2386 / 0.000355D0
        IF (ratco2_var_2387 .GT. 3.0D0) THEN
          adjfac_var_2388 = 2.0D0 + (ratco2_var_2387 - 2.0D0) ** 0.68D0
          adjcolco2_var_2389 = adjfac_var_2388 * 0.000355D0 * coldry_var_2333(jl_var_2412, lay_var_2355) * 1D-20
        ELSE
          adjcolco2_var_2389 = colco2_var_2331(jl_var_2412, lay_var_2355)
        END IF
        speccomb_mco = colh2o_var_2329(jl_var_2412, lay_var_2355) + refrat_m_a3 * coln2o_var_2330(jl_var_2412, lay_var_2355)
        specparm_mco = MIN(colh2o_var_2329(jl_var_2412, lay_var_2355) / speccomb_mco, oneminus_var_2328)
        specmult_mco = 8.0D0 * specparm_mco
        jmco = 1 + INT(specmult_mco)
        fmco = ((specmult_mco) - AINT((specmult_mco)))
        speccomb_planck_var_2346 = colh2o_var_2329(jl_var_2412, lay_var_2355) + refrat_planck_a_var_2359 * coln2o_var_2330(jl_var_2412, lay_var_2355)
        specparm_planck_var_2401 = MIN(colh2o_var_2329(jl_var_2412, lay_var_2355) / speccomb_planck_var_2346, oneminus_var_2328)
        specmult_planck_var_2400 = 8.0D0 * specparm_planck_var_2401
        jpl_var_2357 = 1 + INT(specmult_planck_var_2400)
        fpl_var_2399 = ((specmult_planck_var_2400) - AINT((specmult_planck_var_2400)))
        ind0_var_2348 = ((jp_var_2325(jl_var_2412, lay_var_2355) - 1) * 5 + (jt_var_2326(jl_var_2412, lay_var_2355) - 1)) * nspa_var_216(13) + js_var_2354
        ind1_var_2349 = (jp_var_2325(jl_var_2412, lay_var_2355) * 5 + (jt1_var_2327(jl_var_2412, lay_var_2355) - 1)) * nspa_var_216(13) + js1_var_2356
        inds_var_2350 = indself_var_2337(jl_var_2412, lay_var_2355)
        indf_var_2351 = indfor_var_2339(jl_var_2412, lay_var_2355)
        indm_var_2352 = indminor_var_2343(jl_var_2412, lay_var_2355)
        IF (specparm_var_2392 .LT. 0.125D0) THEN
          p_var_2373 = fs_var_2390 - 1.0D0
          p4_var_2374 = p_var_2373 ** 4
          fk0_var_2375 = p4_var_2374
          fk1_var_2376 = 1.0D0 - p_var_2373 - 2.0D0 * p4_var_2374
          fk2_var_2377 = p_var_2373 + p4_var_2374
          fac000_var_2361 = fk0_var_2375 * fac00_var_2321(jl_var_2412, lay_var_2355)
          fac100_var_2362 = fk1_var_2376 * fac00_var_2321(jl_var_2412, lay_var_2355)
          fac200_var_2363 = fk2_var_2377 * fac00_var_2321(jl_var_2412, lay_var_2355)
          fac010_var_2364 = fk0_var_2375 * fac10_var_2323(jl_var_2412, lay_var_2355)
          fac110_var_2365 = fk1_var_2376 * fac10_var_2323(jl_var_2412, lay_var_2355)
          fac210_var_2366 = fk2_var_2377 * fac10_var_2323(jl_var_2412, lay_var_2355)
        ELSE IF (specparm_var_2392 .GT. 0.875D0) THEN
          p_var_2373 = - fs_var_2390
          p4_var_2374 = p_var_2373 ** 4
          fk0_var_2375 = p4_var_2374
          fk1_var_2376 = 1.0D0 - p_var_2373 - 2.0D0 * p4_var_2374
          fk2_var_2377 = p_var_2373 + p4_var_2374
          fac000_var_2361 = fk0_var_2375 * fac00_var_2321(jl_var_2412, lay_var_2355)
          fac100_var_2362 = fk1_var_2376 * fac00_var_2321(jl_var_2412, lay_var_2355)
          fac200_var_2363 = fk2_var_2377 * fac00_var_2321(jl_var_2412, lay_var_2355)
          fac010_var_2364 = fk0_var_2375 * fac10_var_2323(jl_var_2412, lay_var_2355)
          fac110_var_2365 = fk1_var_2376 * fac10_var_2323(jl_var_2412, lay_var_2355)
          fac210_var_2366 = fk2_var_2377 * fac10_var_2323(jl_var_2412, lay_var_2355)
        ELSE
          fac000_var_2361 = (1.0D0 - fs_var_2390) * fac00_var_2321(jl_var_2412, lay_var_2355)
          fac010_var_2364 = (1.0D0 - fs_var_2390) * fac10_var_2323(jl_var_2412, lay_var_2355)
          fac100_var_2362 = fs_var_2390 * fac00_var_2321(jl_var_2412, lay_var_2355)
          fac110_var_2365 = fs_var_2390 * fac10_var_2323(jl_var_2412, lay_var_2355)
          fac200_var_2363 = 0.0D0
          fac210_var_2366 = 0.0D0
        END IF
        IF (specparm1_var_2395 .LT. 0.125D0) THEN
          p_var_2373 = fs1_var_2393 - 1.0D0
          p4_var_2374 = p_var_2373 ** 4
          fk0_var_2375 = p4_var_2374
          fk1_var_2376 = 1.0D0 - p_var_2373 - 2.0D0 * p4_var_2374
          fk2_var_2377 = p_var_2373 + p4_var_2374
          fac001_var_2367 = fk0_var_2375 * fac01_var_2322(jl_var_2412, lay_var_2355)
          fac101_var_2368 = fk1_var_2376 * fac01_var_2322(jl_var_2412, lay_var_2355)
          fac201_var_2369 = fk2_var_2377 * fac01_var_2322(jl_var_2412, lay_var_2355)
          fac011_var_2370 = fk0_var_2375 * fac11_var_2324(jl_var_2412, lay_var_2355)
          fac111_var_2371 = fk1_var_2376 * fac11_var_2324(jl_var_2412, lay_var_2355)
          fac211_var_2372 = fk2_var_2377 * fac11_var_2324(jl_var_2412, lay_var_2355)
        ELSE IF (specparm1_var_2395 .GT. 0.875D0) THEN
          p_var_2373 = - fs1_var_2393
          p4_var_2374 = p_var_2373 ** 4
          fk0_var_2375 = p4_var_2374
          fk1_var_2376 = 1.0D0 - p_var_2373 - 2.0D0 * p4_var_2374
          fk2_var_2377 = p_var_2373 + p4_var_2374
          fac001_var_2367 = fk0_var_2375 * fac01_var_2322(jl_var_2412, lay_var_2355)
          fac101_var_2368 = fk1_var_2376 * fac01_var_2322(jl_var_2412, lay_var_2355)
          fac201_var_2369 = fk2_var_2377 * fac01_var_2322(jl_var_2412, lay_var_2355)
          fac011_var_2370 = fk0_var_2375 * fac11_var_2324(jl_var_2412, lay_var_2355)
          fac111_var_2371 = fk1_var_2376 * fac11_var_2324(jl_var_2412, lay_var_2355)
          fac211_var_2372 = fk2_var_2377 * fac11_var_2324(jl_var_2412, lay_var_2355)
        ELSE
          fac001_var_2367 = (1.0D0 - fs1_var_2393) * fac01_var_2322(jl_var_2412, lay_var_2355)
          fac011_var_2370 = (1.0D0 - fs1_var_2393) * fac11_var_2324(jl_var_2412, lay_var_2355)
          fac101_var_2368 = fs1_var_2393 * fac01_var_2322(jl_var_2412, lay_var_2355)
          fac111_var_2371 = fs1_var_2393 * fac11_var_2324(jl_var_2412, lay_var_2355)
          fac201_var_2369 = 0.0D0
          fac211_var_2372 = 0.0D0
        END IF
        IF (specparm_var_2392 .LT. 0.125D0) THEN
          tau_major_var_2380(1 : ng13) = speccomb_var_2344 * (fac000_var_2361 * absa_var_132(ind0_var_2348, 1 : 4) + fac100_var_2362 * absa_var_132(ind0_var_2348 + 1, 1 : 4) + fac200_var_2363 * absa_var_132(ind0_var_2348 + 2, 1 : 4) + fac010_var_2364 * absa_var_132(ind0_var_2348 + 9, 1 : 4) + fac110_var_2365 * absa_var_132(ind0_var_2348 + 10, 1 : 4) + fac210_var_2366 * absa_var_132(ind0_var_2348 + 11, 1 : 4))
        ELSE IF (specparm_var_2392 .GT. 0.875D0) THEN
          tau_major_var_2380(1 : ng13) = speccomb_var_2344 * (fac200_var_2363 * absa_var_132(ind0_var_2348 - 1, 1 : 4) + fac100_var_2362 * absa_var_132(ind0_var_2348, 1 : 4) + fac000_var_2361 * absa_var_132(ind0_var_2348 + 1, 1 : 4) + fac210_var_2366 * absa_var_132(ind0_var_2348 + 8, 1 : 4) + fac110_var_2365 * absa_var_132(ind0_var_2348 + 9, 1 : 4) + fac010_var_2364 * absa_var_132(ind0_var_2348 + 10, 1 : 4))
        ELSE
          tau_major_var_2380(1 : ng13) = speccomb_var_2344 * (fac000_var_2361 * absa_var_132(ind0_var_2348, 1 : 4) + fac100_var_2362 * absa_var_132(ind0_var_2348 + 1, 1 : 4) + fac010_var_2364 * absa_var_132(ind0_var_2348 + 9, 1 : 4) + fac110_var_2365 * absa_var_132(ind0_var_2348 + 10, 1 : 4))
        END IF
        IF (specparm1_var_2395 .LT. 0.125D0) THEN
          tau_major1_var_2381(1 : ng13) = speccomb1_var_2345 * (fac001_var_2367 * absa_var_132(ind1_var_2349, 1 : 4) + fac101_var_2368 * absa_var_132(ind1_var_2349 + 1, 1 : 4) + fac201_var_2369 * absa_var_132(ind1_var_2349 + 2, 1 : 4) + fac011_var_2370 * absa_var_132(ind1_var_2349 + 9, 1 : 4) + fac111_var_2371 * absa_var_132(ind1_var_2349 + 10, 1 : 4) + fac211_var_2372 * absa_var_132(ind1_var_2349 + 11, 1 : 4))
        ELSE IF (specparm1_var_2395 .GT. 0.875D0) THEN
          tau_major1_var_2381(1 : ng13) = speccomb1_var_2345 * (fac201_var_2369 * absa_var_132(ind1_var_2349 - 1, 1 : 4) + fac101_var_2368 * absa_var_132(ind1_var_2349, 1 : 4) + fac001_var_2367 * absa_var_132(ind1_var_2349 + 1, 1 : 4) + fac211_var_2372 * absa_var_132(ind1_var_2349 + 8, 1 : 4) + fac111_var_2371 * absa_var_132(ind1_var_2349 + 9, 1 : 4) + fac011_var_2370 * absa_var_132(ind1_var_2349 + 10, 1 : 4))
        ELSE
          tau_major1_var_2381(1 : ng13) = speccomb1_var_2345 * (fac001_var_2367 * absa_var_132(ind1_var_2349, 1 : 4) + fac101_var_2368 * absa_var_132(ind1_var_2349 + 1, 1 : 4) + fac011_var_2370 * absa_var_132(ind1_var_2349 + 9, 1 : 4) + fac111_var_2371 * absa_var_132(ind1_var_2349 + 10, 1 : 4))
        END IF
        DO ig_var_2353 = 1, 4
          tauself_var_2379 = selffac_var_2335(jl_var_2412, lay_var_2355) * (selfref_var_133(inds_var_2350, ig_var_2353) + selffrac_var_2336(jl_var_2412, lay_var_2355) * (selfref_var_133(inds_var_2350 + 1, ig_var_2353) - selfref_var_133(inds_var_2350, ig_var_2353)))
          taufor_var_2378 = forfac_var_2340(jl_var_2412, lay_var_2355) * (forref_var_134(indf_var_2351, ig_var_2353) + forfrac_var_2341(jl_var_2412, lay_var_2355) * (forref_var_134(indf_var_2351 + 1, ig_var_2353) - forref_var_134(indf_var_2351, ig_var_2353)))
          co2m1_var_2382 = ka_mco2_var_135(jmco2_var_2358, indm_var_2352, ig_var_2353) + fmco2_var_2396 * (ka_mco2_var_135(jmco2_var_2358 + 1, indm_var_2352, ig_var_2353) - ka_mco2_var_135(jmco2_var_2358, indm_var_2352, ig_var_2353))
          co2m2_var_2383 = ka_mco2_var_135(jmco2_var_2358, indm_var_2352 + 1, ig_var_2353) + fmco2_var_2396 * (ka_mco2_var_135(jmco2_var_2358 + 1, indm_var_2352 + 1, ig_var_2353) - ka_mco2_var_135(jmco2_var_2358, indm_var_2352 + 1, ig_var_2353))
          absco2_var_2384 = co2m1_var_2382 + minorfrac_var_2342(jl_var_2412, lay_var_2355) * (co2m2_var_2383 - co2m1_var_2382)
          com1 = ka_mco(jmco, indm_var_2352, ig_var_2353) + fmco * (ka_mco(jmco + 1, indm_var_2352, ig_var_2353) - ka_mco(jmco, indm_var_2352, ig_var_2353))
          com2 = ka_mco(jmco, indm_var_2352 + 1, ig_var_2353) + fmco * (ka_mco(jmco + 1, indm_var_2352 + 1, ig_var_2353) - ka_mco(jmco, indm_var_2352 + 1, ig_var_2353))
          absco = com1 + minorfrac_var_2342(jl_var_2412, lay_var_2355) * (com2 - com1)
          taug_var_2319(jl_var_2412, 130 + ig_var_2353, lay_var_2355) = tau_major_var_2380(ig_var_2353) + tau_major1_var_2381(ig_var_2353) + tauself_var_2379 + taufor_var_2378 + adjcolco2_var_2389 * absco2_var_2384 + colco(jl_var_2412, lay_var_2355) * absco
          fracs_var_2338(jl_var_2412, 130 + ig_var_2353, lay_var_2355) = fracrefa_var_130(ig_var_2353, jpl_var_2357) + fpl_var_2399 * (fracrefa_var_130(ig_var_2353, jpl_var_2357 + 1) - fracrefa_var_130(ig_var_2353, jpl_var_2357))
        END DO
      END DO
      ixc0_var_2409 = kfdia_var_2317 - kidia_var_2316 + 1 - ixc0_var_2409
      DO ixp_var_2410 = 1, ixc0_var_2409
        jl_var_2412 = ixhigh_var_2406(ixp_var_2410, lay_var_2355)
        indm_var_2352 = indminor_var_2343(jl_var_2412, lay_var_2355)
        DO ig_var_2353 = 1, 4
          abso3_var_2385 = kb_mo3(indm_var_2352, ig_var_2353) + minorfrac_var_2342(jl_var_2412, lay_var_2355) * (kb_mo3(indm_var_2352 + 1, ig_var_2353) - kb_mo3(indm_var_2352, ig_var_2353))
          taug_var_2319(jl_var_2412, 130 + ig_var_2353, lay_var_2355) = colo3_var_2332(jl_var_2412, lay_var_2355) * abso3_var_2385
          fracs_var_2338(jl_var_2412, 130 + ig_var_2353, lay_var_2355) = fracrefb_var_131(ig_var_2353)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol13
SUBROUTINE global_init_fn
  USE yomlun_ifsaux, ONLY: nulout
  IMPLICIT NONE
  nulout = 6
END SUBROUTINE global_init_fn