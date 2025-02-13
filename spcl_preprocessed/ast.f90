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
    do_overlap_conversion = .FALSE.
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
  INTEGER(KIND = 4), PARAMETER :: nulout = 6
  INTEGER(KIND = 4), PARAMETER :: nulerr = 0
END MODULE yomlun_ifsaux
MODULE radiation_io
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE radiation_abort(text)
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: text
    ERROR STOP 'error in radiation scheme'
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
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: mixing_ratio, effective_radius
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
          lwp_in_cloud = factor_var_449 * cloud_var_442 % mixing_ratio(jcol_var_450, jlev_var_451, 1)
          iwp_in_cloud = factor_var_449 * cloud_var_442 % mixing_ratio(jcol_var_450, jlev_var_451, 2)
          IF (lwp_in_cloud > 0.0D0) THEN
            CALL calc_liq_optics_socrates(16, config_var_440 % cloud_optics % liq_coeff_lw, lwp_in_cloud, cloud_var_442 % effective_radius(jcol_var_450, jlev_var_451, 1), od_lw_liq, scat_od_lw_liq, g_lw_liq)
            CALL calc_liq_optics_socrates(14, config_var_440 % cloud_optics % liq_coeff_sw, lwp_in_cloud, cloud_var_442 % effective_radius(jcol_var_450, jlev_var_451, 1), od_sw_liq, scat_od_sw_liq, g_sw_liq)
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
            CALL calc_ice_optics_fu_lw(16, config_var_440 % cloud_optics % ice_coeff_lw, iwp_in_cloud, cloud_var_442 % effective_radius(jcol_var_450, jlev_var_451, 2), od_lw_ice, scat_od_lw_ice, g_lw_ice)
            CALL calc_ice_optics_fu_sw(14, config_var_440 % cloud_optics % ice_coeff_sw, iwp_in_cloud, cloud_var_442 % effective_radius(jcol_var_450, jlev_var_451, 2), od_sw_ice, scat_od_sw_ice, g_sw_ice)
            CALL calc_ice_optics_yi_lw(16, config_var_440 % cloud_optics % ice_coeff_lw, iwp_in_cloud, cloud_var_442 % effective_radius(jcol_var_450, jlev_var_451, 2), od_lw_ice, scat_od_lw_ice, g_lw_ice)
            CALL calc_ice_optics_yi_sw(14, config_var_440 % cloud_optics % ice_coeff_sw, iwp_in_cloud, cloud_var_442 % effective_radius(jcol_var_450, jlev_var_451, 2), od_sw_ice, scat_od_sw_ice, g_sw_ice)
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
  WRITE(nulout, '(1x,a)') cdtext
  WRITE(nulerr, '(1x,a,a)') 'abort! ', cdtext
  ERROR STOP 1
END SUBROUTINE abor1
SUBROUTINE srtm_taumol27(kidia_var_603, kfdia_var_604, klev_var_605, p_fac00_var_606, p_fac01_var_607, p_fac10_var_608, p_fac11_var_609, k_jp_var_610, k_jt_var_611, k_jt1_var_612, p_colmol_var_613, p_colo3_var_614, k_laytrop_var_615, p_sfluxzen_var_616, p_taug_var_617, p_taur_var_618, prmu0_var_619)
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  USE yoesrta27, ONLY: absa_var_296, absb_var_297, layreffr_var_295, raylc_var_299, scalekur, sfluxrefc_var_298
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_603, kfdia_var_604
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_605
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_606(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_607(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_608(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_609(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_610(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_611(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_612(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_613(kidia_var_603 : kfdia_var_604, klev_var_605)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_614(kidia_var_603 : kfdia_var_604, klev_var_605)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_615(kidia_var_603 : kfdia_var_604)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_616(kidia_var_603 : kfdia_var_604, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_617(kidia_var_603 : kfdia_var_604, klev_var_605, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_618(kidia_var_603 : kfdia_var_604, klev_var_605, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_619(kidia_var_603 : kfdia_var_604)
  INTEGER(KIND = 4) :: ig_var_620, ind0_var_621, ind1_var_622, i_lay_var_623, i_laysolfr_var_624(kidia_var_603 : kfdia_var_604), i_nlayers_var_625, iplon_var_626
  INTEGER(KIND = 4) :: laytrop_min_var_627, laytrop_max_var_628
  REAL(KIND = 8) :: z_tauray_var_629
  laytrop_min_var_627 = MINVAL(k_laytrop_var_615(kidia_var_603 : kfdia_var_604))
  laytrop_max_var_628 = MAXVAL(k_laytrop_var_615(kidia_var_603 : kfdia_var_604))
  i_nlayers_var_625 = klev_var_605
  DO iplon_var_626 = kidia_var_603, kfdia_var_604
    i_laysolfr_var_624(iplon_var_626) = i_nlayers_var_625
  END DO
  DO i_lay_var_623 = 1, laytrop_min_var_627
    DO iplon_var_626 = kidia_var_603, kfdia_var_604
      ind0_var_621 = ((k_jp_var_610(iplon_var_626, i_lay_var_623) - 1) * 5 + (k_jt_var_611(iplon_var_626, i_lay_var_623) - 1)) * nspa_var_313(27) + 1
      ind1_var_622 = (k_jp_var_610(iplon_var_626, i_lay_var_623) * 5 + (k_jt1_var_612(iplon_var_626, i_lay_var_623) - 1)) * nspa_var_313(27) + 1
      DO ig_var_620 = 1, 8
        z_tauray_var_629 = p_colmol_var_613(iplon_var_626, i_lay_var_623) * raylc_var_299(ig_var_620)
        p_taug_var_617(iplon_var_626, i_lay_var_623, ig_var_620) = p_colo3_var_614(iplon_var_626, i_lay_var_623) * (p_fac00_var_606(iplon_var_626, i_lay_var_623) * absa_var_296(ind0_var_621, ig_var_620) + p_fac10_var_608(iplon_var_626, i_lay_var_623) * absa_var_296(ind0_var_621 + 1, ig_var_620) + p_fac01_var_607(iplon_var_626, i_lay_var_623) * absa_var_296(ind1_var_622, ig_var_620) + p_fac11_var_609(iplon_var_626, i_lay_var_623) * absa_var_296(ind1_var_622 + 1, ig_var_620))
        p_taur_var_618(iplon_var_626, i_lay_var_623, ig_var_620) = z_tauray_var_629
      END DO
    END DO
  END DO
  DO i_lay_var_623 = laytrop_min_var_627 + 1, laytrop_max_var_628
    DO iplon_var_626 = kidia_var_603, kfdia_var_604
      IF (i_lay_var_623 <= k_laytrop_var_615(iplon_var_626)) THEN
        ind0_var_621 = ((k_jp_var_610(iplon_var_626, i_lay_var_623) - 1) * 5 + (k_jt_var_611(iplon_var_626, i_lay_var_623) - 1)) * nspa_var_313(27) + 1
        ind1_var_622 = (k_jp_var_610(iplon_var_626, i_lay_var_623) * 5 + (k_jt1_var_612(iplon_var_626, i_lay_var_623) - 1)) * nspa_var_313(27) + 1
        DO ig_var_620 = 1, 8
          z_tauray_var_629 = p_colmol_var_613(iplon_var_626, i_lay_var_623) * raylc_var_299(ig_var_620)
          p_taug_var_617(iplon_var_626, i_lay_var_623, ig_var_620) = p_colo3_var_614(iplon_var_626, i_lay_var_623) * (p_fac00_var_606(iplon_var_626, i_lay_var_623) * absa_var_296(ind0_var_621, ig_var_620) + p_fac10_var_608(iplon_var_626, i_lay_var_623) * absa_var_296(ind0_var_621 + 1, ig_var_620) + p_fac01_var_607(iplon_var_626, i_lay_var_623) * absa_var_296(ind1_var_622, ig_var_620) + p_fac11_var_609(iplon_var_626, i_lay_var_623) * absa_var_296(ind1_var_622 + 1, ig_var_620))
          p_taur_var_618(iplon_var_626, i_lay_var_623, ig_var_620) = z_tauray_var_629
        END DO
      ELSE
        IF (k_jp_var_610(iplon_var_626, i_lay_var_623 - 1) < layreffr_var_295 .AND. k_jp_var_610(iplon_var_626, i_lay_var_623) >= layreffr_var_295) i_laysolfr_var_624(iplon_var_626) = i_lay_var_623
        ind0_var_621 = ((k_jp_var_610(iplon_var_626, i_lay_var_623) - 13) * 5 + (k_jt_var_611(iplon_var_626, i_lay_var_623) - 1)) * nspb_var_314(27) + 1
        ind1_var_622 = ((k_jp_var_610(iplon_var_626, i_lay_var_623) - 12) * 5 + (k_jt1_var_612(iplon_var_626, i_lay_var_623) - 1)) * nspb_var_314(27) + 1
        DO ig_var_620 = 1, 8
          z_tauray_var_629 = p_colmol_var_613(iplon_var_626, i_lay_var_623) * raylc_var_299(ig_var_620)
          p_taug_var_617(iplon_var_626, i_lay_var_623, ig_var_620) = p_colo3_var_614(iplon_var_626, i_lay_var_623) * (p_fac00_var_606(iplon_var_626, i_lay_var_623) * absb_var_297(ind0_var_621, ig_var_620) + p_fac10_var_608(iplon_var_626, i_lay_var_623) * absb_var_297(ind0_var_621 + 1, ig_var_620) + p_fac01_var_607(iplon_var_626, i_lay_var_623) * absb_var_297(ind1_var_622, ig_var_620) + p_fac11_var_609(iplon_var_626, i_lay_var_623) * absb_var_297(ind1_var_622 + 1, ig_var_620))
          IF (i_lay_var_623 == i_laysolfr_var_624(iplon_var_626)) p_sfluxzen_var_616(iplon_var_626, ig_var_620) = scalekur * sfluxrefc_var_298(ig_var_620)
          p_taur_var_618(iplon_var_626, i_lay_var_623, ig_var_620) = z_tauray_var_629
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_623 = laytrop_max_var_628 + 1, i_nlayers_var_625
    DO iplon_var_626 = kidia_var_603, kfdia_var_604
      IF (k_jp_var_610(iplon_var_626, i_lay_var_623 - 1) < layreffr_var_295 .AND. k_jp_var_610(iplon_var_626, i_lay_var_623) >= layreffr_var_295) i_laysolfr_var_624(iplon_var_626) = i_lay_var_623
      ind0_var_621 = ((k_jp_var_610(iplon_var_626, i_lay_var_623) - 13) * 5 + (k_jt_var_611(iplon_var_626, i_lay_var_623) - 1)) * nspb_var_314(27) + 1
      ind1_var_622 = ((k_jp_var_610(iplon_var_626, i_lay_var_623) - 12) * 5 + (k_jt1_var_612(iplon_var_626, i_lay_var_623) - 1)) * nspb_var_314(27) + 1
      DO ig_var_620 = 1, 8
        z_tauray_var_629 = p_colmol_var_613(iplon_var_626, i_lay_var_623) * raylc_var_299(ig_var_620)
        p_taug_var_617(iplon_var_626, i_lay_var_623, ig_var_620) = p_colo3_var_614(iplon_var_626, i_lay_var_623) * (p_fac00_var_606(iplon_var_626, i_lay_var_623) * absb_var_297(ind0_var_621, ig_var_620) + p_fac10_var_608(iplon_var_626, i_lay_var_623) * absb_var_297(ind0_var_621 + 1, ig_var_620) + p_fac01_var_607(iplon_var_626, i_lay_var_623) * absb_var_297(ind1_var_622, ig_var_620) + p_fac11_var_609(iplon_var_626, i_lay_var_623) * absb_var_297(ind1_var_622 + 1, ig_var_620))
        IF (i_lay_var_623 == i_laysolfr_var_624(iplon_var_626)) p_sfluxzen_var_616(iplon_var_626, ig_var_620) = scalekur * sfluxrefc_var_298(ig_var_620)
        p_taur_var_618(iplon_var_626, i_lay_var_623, ig_var_620) = z_tauray_var_629
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol27
SUBROUTINE srtm_taumol26(kidia_var_630, kfdia_var_631, klev_var_632, p_colmol_var_633, k_laytrop_var_634, p_sfluxzen_var_635, p_taug_var_636, p_taur_var_637, prmu0_var_638)
  USE yoesrta26, ONLY: raylc_var_294, sfluxrefc_var_293
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_630, kfdia_var_631
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_632
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_633(kidia_var_630 : kfdia_var_631, klev_var_632)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_634(kidia_var_630 : kfdia_var_631)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_635(kidia_var_630 : kfdia_var_631, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_636(kidia_var_630 : kfdia_var_631, klev_var_632, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_637(kidia_var_630 : kfdia_var_631, klev_var_632, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_638(kidia_var_630 : kfdia_var_631)
  INTEGER(KIND = 4) :: ig_var_639, i_lay_var_640, i_laysolfr_var_641(kidia_var_630 : kfdia_var_631), i_nlayers_var_642, iplon_var_643
  INTEGER(KIND = 4) :: laytrop_min_var_644, laytrop_max_var_645
  laytrop_min_var_644 = MINVAL(k_laytrop_var_634(kidia_var_630 : kfdia_var_631))
  laytrop_max_var_645 = MAXVAL(k_laytrop_var_634(kidia_var_630 : kfdia_var_631))
  i_nlayers_var_642 = klev_var_632
  DO iplon_var_643 = kidia_var_630, kfdia_var_631
    i_laysolfr_var_641(iplon_var_643) = k_laytrop_var_634(iplon_var_643)
  END DO
  DO i_lay_var_640 = 1, laytrop_min_var_644
    DO iplon_var_643 = kidia_var_630, kfdia_var_631
      DO ig_var_639 = 1, 6
        IF (i_lay_var_640 == i_laysolfr_var_641(iplon_var_643)) p_sfluxzen_var_635(iplon_var_643, ig_var_639) = sfluxrefc_var_293(ig_var_639)
        p_taug_var_636(iplon_var_643, i_lay_var_640, ig_var_639) = 0.0D0
        p_taur_var_637(iplon_var_643, i_lay_var_640, ig_var_639) = p_colmol_var_633(iplon_var_643, i_lay_var_640) * raylc_var_294(ig_var_639)
      END DO
    END DO
  END DO
  DO i_lay_var_640 = laytrop_min_var_644 + 1, laytrop_max_var_645
    DO iplon_var_643 = kidia_var_630, kfdia_var_631
      IF (i_lay_var_640 <= k_laytrop_var_634(iplon_var_643)) THEN
        DO ig_var_639 = 1, 6
          IF (i_lay_var_640 == i_laysolfr_var_641(iplon_var_643)) p_sfluxzen_var_635(iplon_var_643, ig_var_639) = sfluxrefc_var_293(ig_var_639)
          p_taug_var_636(iplon_var_643, i_lay_var_640, ig_var_639) = 0.0D0
          p_taur_var_637(iplon_var_643, i_lay_var_640, ig_var_639) = p_colmol_var_633(iplon_var_643, i_lay_var_640) * raylc_var_294(ig_var_639)
        END DO
      ELSE
        DO ig_var_639 = 1, 6
          p_taug_var_636(iplon_var_643, i_lay_var_640, ig_var_639) = 0.0D0
          p_taur_var_637(iplon_var_643, i_lay_var_640, ig_var_639) = p_colmol_var_633(iplon_var_643, i_lay_var_640) * raylc_var_294(ig_var_639)
        END DO
      END IF
    END DO
  END DO
  DO ig_var_639 = 1, 6
    DO i_lay_var_640 = laytrop_max_var_645 + 1, i_nlayers_var_642
      DO iplon_var_643 = kidia_var_630, kfdia_var_631
        p_taug_var_636(iplon_var_643, i_lay_var_640, ig_var_639) = 0.0D0
        p_taur_var_637(iplon_var_643, i_lay_var_640, ig_var_639) = p_colmol_var_633(iplon_var_643, i_lay_var_640) * raylc_var_294(ig_var_639)
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol26
SUBROUTINE srtm_taumol18(kidia_var_646, kfdia_var_647, klev_var_648, p_fac00_var_649, p_fac01_var_650, p_fac10_var_651, p_fac11_var_652, k_jp_var_653, k_jt_var_654, k_jt1_var_655, p_oneminus_var_656, p_colh2o_var_657, p_colch4_var_658, p_colmol_var_659, k_laytrop_var_660, p_selffac_var_661, p_selffrac_var_662, k_indself_var_663, p_forfac_var_664, p_forfrac_var_665, k_indfor_var_666, p_sfluxzen_var_667, p_taug_var_668, p_taur_var_669, prmu0_var_670)
  USE yoesrta18, ONLY: absa_var_236, absb_var_237, forrefc_var_239, layreffr_var_235, rayl_var_233, selfrefc_var_238, sfluxrefc_var_240, strrat_var_234
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_646, kfdia_var_647
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_648
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_649(kidia_var_646 : kfdia_var_647, klev_var_648)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_650(kidia_var_646 : kfdia_var_647, klev_var_648)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_651(kidia_var_646 : kfdia_var_647, klev_var_648)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_652(kidia_var_646 : kfdia_var_647, klev_var_648)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_653(kidia_var_646 : kfdia_var_647, klev_var_648)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_654(kidia_var_646 : kfdia_var_647, klev_var_648)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_655(kidia_var_646 : kfdia_var_647, klev_var_648)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_656(kidia_var_646 : kfdia_var_647)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_657(kidia_var_646 : kfdia_var_647, klev_var_648)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_658(kidia_var_646 : kfdia_var_647, klev_var_648)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_659(kidia_var_646 : kfdia_var_647, klev_var_648)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_660(kidia_var_646 : kfdia_var_647)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_661(kidia_var_646 : kfdia_var_647, klev_var_648)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_662(kidia_var_646 : kfdia_var_647, klev_var_648)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_663(kidia_var_646 : kfdia_var_647, klev_var_648)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_664(kidia_var_646 : kfdia_var_647, klev_var_648)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_665(kidia_var_646 : kfdia_var_647, klev_var_648)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_666(kidia_var_646 : kfdia_var_647, klev_var_648)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_667(kidia_var_646 : kfdia_var_647, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_668(kidia_var_646 : kfdia_var_647, klev_var_648, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_669(kidia_var_646 : kfdia_var_647, klev_var_648, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_670(kidia_var_646 : kfdia_var_647)
  INTEGER(KIND = 4) :: ig_var_671, ind0_var_672, ind1_var_673, inds_var_674, indf_var_675, js_var_676, i_lay_var_677, i_laysolfr_var_678(kidia_var_646 : kfdia_var_647), i_nlayers_var_679, iplon_var_680
  INTEGER(KIND = 4) :: laytrop_min_var_681, laytrop_max_var_682
  REAL(KIND = 8) :: z_fs_var_683, z_speccomb_var_684, z_specmult_var_685, z_specparm_var_686, z_tauray_var_687
  laytrop_min_var_681 = MINVAL(k_laytrop_var_660(kidia_var_646 : kfdia_var_647))
  laytrop_max_var_682 = MAXVAL(k_laytrop_var_660(kidia_var_646 : kfdia_var_647))
  i_nlayers_var_679 = klev_var_648
  DO iplon_var_680 = kidia_var_646, kfdia_var_647
    i_laysolfr_var_678(iplon_var_680) = k_laytrop_var_660(iplon_var_680)
  END DO
  DO i_lay_var_677 = 1, laytrop_min_var_681
    DO iplon_var_680 = kidia_var_646, kfdia_var_647
      IF (k_jp_var_653(iplon_var_680, i_lay_var_677) < layreffr_var_235 .AND. k_jp_var_653(iplon_var_680, i_lay_var_677 + 1) >= layreffr_var_235) i_laysolfr_var_678(iplon_var_680) = MIN(i_lay_var_677 + 1, k_laytrop_var_660(iplon_var_680))
      z_speccomb_var_684 = p_colh2o_var_657(iplon_var_680, i_lay_var_677) + strrat_var_234 * p_colch4_var_658(iplon_var_680, i_lay_var_677)
      z_specparm_var_686 = p_colh2o_var_657(iplon_var_680, i_lay_var_677) / z_speccomb_var_684
      z_specparm_var_686 = MIN(p_oneminus_var_656(iplon_var_680), z_specparm_var_686)
      z_specmult_var_685 = 8.0D0 * (z_specparm_var_686)
      js_var_676 = 1 + INT(z_specmult_var_685)
      z_fs_var_683 = z_specmult_var_685 - AINT(z_specmult_var_685)
      ind0_var_672 = ((k_jp_var_653(iplon_var_680, i_lay_var_677) - 1) * 5 + (k_jt_var_654(iplon_var_680, i_lay_var_677) - 1)) * nspa_var_313(18) + js_var_676
      ind1_var_673 = (k_jp_var_653(iplon_var_680, i_lay_var_677) * 5 + (k_jt1_var_655(iplon_var_680, i_lay_var_677) - 1)) * nspa_var_313(18) + js_var_676
      inds_var_674 = k_indself_var_663(iplon_var_680, i_lay_var_677)
      indf_var_675 = k_indfor_var_666(iplon_var_680, i_lay_var_677)
      z_tauray_var_687 = p_colmol_var_659(iplon_var_680, i_lay_var_677) * rayl_var_233
      DO ig_var_671 = 1, 8
        p_taug_var_668(iplon_var_680, i_lay_var_677, ig_var_671) = z_speccomb_var_684 * ((1.0D0 - z_fs_var_683) * (absa_var_236(ind0_var_672, ig_var_671) * p_fac00_var_649(iplon_var_680, i_lay_var_677) + absa_var_236(ind0_var_672 + 9, ig_var_671) * p_fac10_var_651(iplon_var_680, i_lay_var_677) + absa_var_236(ind1_var_673, ig_var_671) * p_fac01_var_650(iplon_var_680, i_lay_var_677) + absa_var_236(ind1_var_673 + 9, ig_var_671) * p_fac11_var_652(iplon_var_680, i_lay_var_677)) + z_fs_var_683 * (absa_var_236(ind0_var_672 + 1, ig_var_671) * p_fac00_var_649(iplon_var_680, i_lay_var_677) + absa_var_236(ind0_var_672 + 10, ig_var_671) * p_fac10_var_651(iplon_var_680, i_lay_var_677) + absa_var_236(ind1_var_673 + 1, ig_var_671) * p_fac01_var_650(iplon_var_680, i_lay_var_677) + absa_var_236(ind1_var_673 + 10, ig_var_671) * p_fac11_var_652(iplon_var_680, i_lay_var_677))) + p_colh2o_var_657(iplon_var_680, i_lay_var_677) * (p_selffac_var_661(iplon_var_680, i_lay_var_677) * (selfrefc_var_238(inds_var_674, ig_var_671) + p_selffrac_var_662(iplon_var_680, i_lay_var_677) * (selfrefc_var_238(inds_var_674 + 1, ig_var_671) - selfrefc_var_238(inds_var_674, ig_var_671))) + p_forfac_var_664(iplon_var_680, i_lay_var_677) * (forrefc_var_239(indf_var_675, ig_var_671) + p_forfrac_var_665(iplon_var_680, i_lay_var_677) * (forrefc_var_239(indf_var_675 + 1, ig_var_671) - forrefc_var_239(indf_var_675, ig_var_671))))
        IF (i_lay_var_677 == i_laysolfr_var_678(iplon_var_680)) p_sfluxzen_var_667(iplon_var_680, ig_var_671) = sfluxrefc_var_240(ig_var_671, js_var_676) + z_fs_var_683 * (sfluxrefc_var_240(ig_var_671, js_var_676 + 1) - sfluxrefc_var_240(ig_var_671, js_var_676))
        p_taur_var_669(iplon_var_680, i_lay_var_677, ig_var_671) = z_tauray_var_687
      END DO
    END DO
  END DO
  DO i_lay_var_677 = laytrop_min_var_681 + 1, laytrop_max_var_682
    DO iplon_var_680 = kidia_var_646, kfdia_var_647
      IF (i_lay_var_677 <= k_laytrop_var_660(iplon_var_680)) THEN
        IF (k_jp_var_653(iplon_var_680, i_lay_var_677) < layreffr_var_235 .AND. k_jp_var_653(iplon_var_680, i_lay_var_677 + 1) >= layreffr_var_235) i_laysolfr_var_678(iplon_var_680) = MIN(i_lay_var_677 + 1, k_laytrop_var_660(iplon_var_680))
        z_speccomb_var_684 = p_colh2o_var_657(iplon_var_680, i_lay_var_677) + strrat_var_234 * p_colch4_var_658(iplon_var_680, i_lay_var_677)
        z_specparm_var_686 = p_colh2o_var_657(iplon_var_680, i_lay_var_677) / z_speccomb_var_684
        z_specparm_var_686 = MIN(p_oneminus_var_656(iplon_var_680), z_specparm_var_686)
        z_specmult_var_685 = 8.0D0 * (z_specparm_var_686)
        js_var_676 = 1 + INT(z_specmult_var_685)
        z_fs_var_683 = z_specmult_var_685 - AINT(z_specmult_var_685)
        ind0_var_672 = ((k_jp_var_653(iplon_var_680, i_lay_var_677) - 1) * 5 + (k_jt_var_654(iplon_var_680, i_lay_var_677) - 1)) * nspa_var_313(18) + js_var_676
        ind1_var_673 = (k_jp_var_653(iplon_var_680, i_lay_var_677) * 5 + (k_jt1_var_655(iplon_var_680, i_lay_var_677) - 1)) * nspa_var_313(18) + js_var_676
        inds_var_674 = k_indself_var_663(iplon_var_680, i_lay_var_677)
        indf_var_675 = k_indfor_var_666(iplon_var_680, i_lay_var_677)
        z_tauray_var_687 = p_colmol_var_659(iplon_var_680, i_lay_var_677) * rayl_var_233
        DO ig_var_671 = 1, 8
          p_taug_var_668(iplon_var_680, i_lay_var_677, ig_var_671) = z_speccomb_var_684 * ((1.0D0 - z_fs_var_683) * (absa_var_236(ind0_var_672, ig_var_671) * p_fac00_var_649(iplon_var_680, i_lay_var_677) + absa_var_236(ind0_var_672 + 9, ig_var_671) * p_fac10_var_651(iplon_var_680, i_lay_var_677) + absa_var_236(ind1_var_673, ig_var_671) * p_fac01_var_650(iplon_var_680, i_lay_var_677) + absa_var_236(ind1_var_673 + 9, ig_var_671) * p_fac11_var_652(iplon_var_680, i_lay_var_677)) + z_fs_var_683 * (absa_var_236(ind0_var_672 + 1, ig_var_671) * p_fac00_var_649(iplon_var_680, i_lay_var_677) + absa_var_236(ind0_var_672 + 10, ig_var_671) * p_fac10_var_651(iplon_var_680, i_lay_var_677) + absa_var_236(ind1_var_673 + 1, ig_var_671) * p_fac01_var_650(iplon_var_680, i_lay_var_677) + absa_var_236(ind1_var_673 + 10, ig_var_671) * p_fac11_var_652(iplon_var_680, i_lay_var_677))) + p_colh2o_var_657(iplon_var_680, i_lay_var_677) * (p_selffac_var_661(iplon_var_680, i_lay_var_677) * (selfrefc_var_238(inds_var_674, ig_var_671) + p_selffrac_var_662(iplon_var_680, i_lay_var_677) * (selfrefc_var_238(inds_var_674 + 1, ig_var_671) - selfrefc_var_238(inds_var_674, ig_var_671))) + p_forfac_var_664(iplon_var_680, i_lay_var_677) * (forrefc_var_239(indf_var_675, ig_var_671) + p_forfrac_var_665(iplon_var_680, i_lay_var_677) * (forrefc_var_239(indf_var_675 + 1, ig_var_671) - forrefc_var_239(indf_var_675, ig_var_671))))
          IF (i_lay_var_677 == i_laysolfr_var_678(iplon_var_680)) p_sfluxzen_var_667(iplon_var_680, ig_var_671) = sfluxrefc_var_240(ig_var_671, js_var_676) + z_fs_var_683 * (sfluxrefc_var_240(ig_var_671, js_var_676 + 1) - sfluxrefc_var_240(ig_var_671, js_var_676))
          p_taur_var_669(iplon_var_680, i_lay_var_677, ig_var_671) = z_tauray_var_687
        END DO
      ELSE
        ind0_var_672 = ((k_jp_var_653(iplon_var_680, i_lay_var_677) - 13) * 5 + (k_jt_var_654(iplon_var_680, i_lay_var_677) - 1)) * nspb_var_314(18) + 1
        ind1_var_673 = ((k_jp_var_653(iplon_var_680, i_lay_var_677) - 12) * 5 + (k_jt1_var_655(iplon_var_680, i_lay_var_677) - 1)) * nspb_var_314(18) + 1
        z_tauray_var_687 = p_colmol_var_659(iplon_var_680, i_lay_var_677) * rayl_var_233
        DO ig_var_671 = 1, 8
          p_taug_var_668(iplon_var_680, i_lay_var_677, ig_var_671) = p_colch4_var_658(iplon_var_680, i_lay_var_677) * (p_fac00_var_649(iplon_var_680, i_lay_var_677) * absb_var_237(ind0_var_672, ig_var_671) + p_fac10_var_651(iplon_var_680, i_lay_var_677) * absb_var_237(ind0_var_672 + 1, ig_var_671) + p_fac01_var_650(iplon_var_680, i_lay_var_677) * absb_var_237(ind1_var_673, ig_var_671) + p_fac11_var_652(iplon_var_680, i_lay_var_677) * absb_var_237(ind1_var_673 + 1, ig_var_671))
          p_taur_var_669(iplon_var_680, i_lay_var_677, ig_var_671) = z_tauray_var_687
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_677 = laytrop_max_var_682 + 1, i_nlayers_var_679
    DO iplon_var_680 = kidia_var_646, kfdia_var_647
      ind0_var_672 = ((k_jp_var_653(iplon_var_680, i_lay_var_677) - 13) * 5 + (k_jt_var_654(iplon_var_680, i_lay_var_677) - 1)) * nspb_var_314(18) + 1
      ind1_var_673 = ((k_jp_var_653(iplon_var_680, i_lay_var_677) - 12) * 5 + (k_jt1_var_655(iplon_var_680, i_lay_var_677) - 1)) * nspb_var_314(18) + 1
      z_tauray_var_687 = p_colmol_var_659(iplon_var_680, i_lay_var_677) * rayl_var_233
      DO ig_var_671 = 1, 8
        p_taug_var_668(iplon_var_680, i_lay_var_677, ig_var_671) = p_colch4_var_658(iplon_var_680, i_lay_var_677) * (p_fac00_var_649(iplon_var_680, i_lay_var_677) * absb_var_237(ind0_var_672, ig_var_671) + p_fac10_var_651(iplon_var_680, i_lay_var_677) * absb_var_237(ind0_var_672 + 1, ig_var_671) + p_fac01_var_650(iplon_var_680, i_lay_var_677) * absb_var_237(ind1_var_673, ig_var_671) + p_fac11_var_652(iplon_var_680, i_lay_var_677) * absb_var_237(ind1_var_673 + 1, ig_var_671))
        p_taur_var_669(iplon_var_680, i_lay_var_677, ig_var_671) = z_tauray_var_687
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol18
SUBROUTINE srtm_taumol24(kidia_var_688, kfdia_var_689, klev_var_690, p_fac00_var_691, p_fac01_var_692, p_fac10_var_693, p_fac11_var_694, k_jp_var_695, k_jt_var_696, k_jt1_var_697, p_oneminus_var_698, p_colh2o_var_699, p_colmol_var_700, p_colo2_var_701, p_colo3_var_702, k_laytrop_var_703, p_selffac_var_704, p_selffrac_var_705, k_indself_var_706, p_forfac_var_707, p_forfrac_var_708, k_indfor_var_709, p_sfluxzen_var_710, p_taug_var_711, p_taur_var_712, prmu0_var_713)
  USE yoesrta24, ONLY: absa_var_280, absb_var_281, abso3ac_var_285, abso3bc_var_286, forrefc_var_283, layreffr_var_279, raylac, raylbc, selfrefc_var_282, sfluxrefc_var_284, strrat_var_278
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_688, kfdia_var_689
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_690
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_691(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_692(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_693(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_694(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_695(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_696(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_697(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_698(kidia_var_688 : kfdia_var_689)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_699(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_700(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_701(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_702(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_703(kidia_var_688 : kfdia_var_689)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_704(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_705(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_706(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_707(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_708(kidia_var_688 : kfdia_var_689, klev_var_690)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_709(kidia_var_688 : kfdia_var_689, klev_var_690)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_710(kidia_var_688 : kfdia_var_689, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_711(kidia_var_688 : kfdia_var_689, klev_var_690, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_712(kidia_var_688 : kfdia_var_689, klev_var_690, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_713(kidia_var_688 : kfdia_var_689)
  INTEGER(KIND = 4) :: ig_var_714, ind0_var_715, ind1_var_716, inds_var_717, indf_var_718, js_var_719, i_lay_var_720, i_laysolfr_var_721(kidia_var_688 : kfdia_var_689), i_nlayers_var_722, iplon_var_723
  INTEGER(KIND = 4) :: laytrop_min_var_724, laytrop_max_var_725
  REAL(KIND = 8) :: z_fs_var_726, z_speccomb_var_727, z_specmult_var_728, z_specparm_var_729, z_tauray_var_730
  laytrop_min_var_724 = MINVAL(k_laytrop_var_703(kidia_var_688 : kfdia_var_689))
  laytrop_max_var_725 = MAXVAL(k_laytrop_var_703(kidia_var_688 : kfdia_var_689))
  i_nlayers_var_722 = klev_var_690
  DO iplon_var_723 = kidia_var_688, kfdia_var_689
    i_laysolfr_var_721(iplon_var_723) = k_laytrop_var_703(iplon_var_723)
  END DO
  DO i_lay_var_720 = 1, laytrop_min_var_724
    DO iplon_var_723 = kidia_var_688, kfdia_var_689
      IF (k_jp_var_695(iplon_var_723, i_lay_var_720) < layreffr_var_279 .AND. k_jp_var_695(iplon_var_723, i_lay_var_720 + 1) >= layreffr_var_279) i_laysolfr_var_721(iplon_var_723) = MIN(i_lay_var_720 + 1, k_laytrop_var_703(iplon_var_723))
      z_speccomb_var_727 = p_colh2o_var_699(iplon_var_723, i_lay_var_720) + strrat_var_278 * p_colo2_var_701(iplon_var_723, i_lay_var_720)
      z_specparm_var_729 = p_colh2o_var_699(iplon_var_723, i_lay_var_720) / z_speccomb_var_727
      z_specparm_var_729 = MIN(p_oneminus_var_698(iplon_var_723), z_specparm_var_729)
      z_specmult_var_728 = 8.0D0 * (z_specparm_var_729)
      js_var_719 = 1 + INT(z_specmult_var_728)
      z_fs_var_726 = z_specmult_var_728 - AINT(z_specmult_var_728)
      ind0_var_715 = ((k_jp_var_695(iplon_var_723, i_lay_var_720) - 1) * 5 + (k_jt_var_696(iplon_var_723, i_lay_var_720) - 1)) * nspa_var_313(24) + js_var_719
      ind1_var_716 = (k_jp_var_695(iplon_var_723, i_lay_var_720) * 5 + (k_jt1_var_697(iplon_var_723, i_lay_var_720) - 1)) * nspa_var_313(24) + js_var_719
      inds_var_717 = k_indself_var_706(iplon_var_723, i_lay_var_720)
      indf_var_718 = k_indfor_var_709(iplon_var_723, i_lay_var_720)
      DO ig_var_714 = 1, 8
        z_tauray_var_730 = p_colmol_var_700(iplon_var_723, i_lay_var_720) * (raylac(ig_var_714, js_var_719) + z_fs_var_726 * (raylac(ig_var_714, js_var_719 + 1) - raylac(ig_var_714, js_var_719)))
        p_taug_var_711(iplon_var_723, i_lay_var_720, ig_var_714) = z_speccomb_var_727 * ((1.0D0 - z_fs_var_726) * (absa_var_280(ind0_var_715, ig_var_714) * p_fac00_var_691(iplon_var_723, i_lay_var_720) + absa_var_280(ind0_var_715 + 9, ig_var_714) * p_fac10_var_693(iplon_var_723, i_lay_var_720) + absa_var_280(ind1_var_716, ig_var_714) * p_fac01_var_692(iplon_var_723, i_lay_var_720) + absa_var_280(ind1_var_716 + 9, ig_var_714) * p_fac11_var_694(iplon_var_723, i_lay_var_720)) + z_fs_var_726 * (absa_var_280(ind0_var_715 + 1, ig_var_714) * p_fac00_var_691(iplon_var_723, i_lay_var_720) + absa_var_280(ind0_var_715 + 10, ig_var_714) * p_fac10_var_693(iplon_var_723, i_lay_var_720) + absa_var_280(ind1_var_716 + 1, ig_var_714) * p_fac01_var_692(iplon_var_723, i_lay_var_720) + absa_var_280(ind1_var_716 + 10, ig_var_714) * p_fac11_var_694(iplon_var_723, i_lay_var_720))) + p_colo3_var_702(iplon_var_723, i_lay_var_720) * abso3ac_var_285(ig_var_714) + p_colh2o_var_699(iplon_var_723, i_lay_var_720) * (p_selffac_var_704(iplon_var_723, i_lay_var_720) * (selfrefc_var_282(inds_var_717, ig_var_714) + p_selffrac_var_705(iplon_var_723, i_lay_var_720) * (selfrefc_var_282(inds_var_717 + 1, ig_var_714) - selfrefc_var_282(inds_var_717, ig_var_714))) + p_forfac_var_707(iplon_var_723, i_lay_var_720) * (forrefc_var_283(indf_var_718, ig_var_714) + p_forfrac_var_708(iplon_var_723, i_lay_var_720) * (forrefc_var_283(indf_var_718 + 1, ig_var_714) - forrefc_var_283(indf_var_718, ig_var_714))))
        IF (i_lay_var_720 == i_laysolfr_var_721(iplon_var_723)) p_sfluxzen_var_710(iplon_var_723, ig_var_714) = sfluxrefc_var_284(ig_var_714, js_var_719) + z_fs_var_726 * (sfluxrefc_var_284(ig_var_714, js_var_719 + 1) - sfluxrefc_var_284(ig_var_714, js_var_719))
        p_taur_var_712(iplon_var_723, i_lay_var_720, ig_var_714) = z_tauray_var_730
      END DO
    END DO
  END DO
  DO i_lay_var_720 = laytrop_min_var_724 + 1, laytrop_max_var_725
    DO iplon_var_723 = kidia_var_688, kfdia_var_689
      IF (i_lay_var_720 <= k_laytrop_var_703(iplon_var_723)) THEN
        IF (k_jp_var_695(iplon_var_723, i_lay_var_720) < layreffr_var_279 .AND. k_jp_var_695(iplon_var_723, i_lay_var_720 + 1) >= layreffr_var_279) i_laysolfr_var_721(iplon_var_723) = MIN(i_lay_var_720 + 1, k_laytrop_var_703(iplon_var_723))
        z_speccomb_var_727 = p_colh2o_var_699(iplon_var_723, i_lay_var_720) + strrat_var_278 * p_colo2_var_701(iplon_var_723, i_lay_var_720)
        z_specparm_var_729 = p_colh2o_var_699(iplon_var_723, i_lay_var_720) / z_speccomb_var_727
        z_specparm_var_729 = MIN(p_oneminus_var_698(iplon_var_723), z_specparm_var_729)
        z_specmult_var_728 = 8.0D0 * (z_specparm_var_729)
        js_var_719 = 1 + INT(z_specmult_var_728)
        z_fs_var_726 = z_specmult_var_728 - AINT(z_specmult_var_728)
        ind0_var_715 = ((k_jp_var_695(iplon_var_723, i_lay_var_720) - 1) * 5 + (k_jt_var_696(iplon_var_723, i_lay_var_720) - 1)) * nspa_var_313(24) + js_var_719
        ind1_var_716 = (k_jp_var_695(iplon_var_723, i_lay_var_720) * 5 + (k_jt1_var_697(iplon_var_723, i_lay_var_720) - 1)) * nspa_var_313(24) + js_var_719
        inds_var_717 = k_indself_var_706(iplon_var_723, i_lay_var_720)
        indf_var_718 = k_indfor_var_709(iplon_var_723, i_lay_var_720)
        DO ig_var_714 = 1, 8
          z_tauray_var_730 = p_colmol_var_700(iplon_var_723, i_lay_var_720) * (raylac(ig_var_714, js_var_719) + z_fs_var_726 * (raylac(ig_var_714, js_var_719 + 1) - raylac(ig_var_714, js_var_719)))
          p_taug_var_711(iplon_var_723, i_lay_var_720, ig_var_714) = z_speccomb_var_727 * ((1.0D0 - z_fs_var_726) * (absa_var_280(ind0_var_715, ig_var_714) * p_fac00_var_691(iplon_var_723, i_lay_var_720) + absa_var_280(ind0_var_715 + 9, ig_var_714) * p_fac10_var_693(iplon_var_723, i_lay_var_720) + absa_var_280(ind1_var_716, ig_var_714) * p_fac01_var_692(iplon_var_723, i_lay_var_720) + absa_var_280(ind1_var_716 + 9, ig_var_714) * p_fac11_var_694(iplon_var_723, i_lay_var_720)) + z_fs_var_726 * (absa_var_280(ind0_var_715 + 1, ig_var_714) * p_fac00_var_691(iplon_var_723, i_lay_var_720) + absa_var_280(ind0_var_715 + 10, ig_var_714) * p_fac10_var_693(iplon_var_723, i_lay_var_720) + absa_var_280(ind1_var_716 + 1, ig_var_714) * p_fac01_var_692(iplon_var_723, i_lay_var_720) + absa_var_280(ind1_var_716 + 10, ig_var_714) * p_fac11_var_694(iplon_var_723, i_lay_var_720))) + p_colo3_var_702(iplon_var_723, i_lay_var_720) * abso3ac_var_285(ig_var_714) + p_colh2o_var_699(iplon_var_723, i_lay_var_720) * (p_selffac_var_704(iplon_var_723, i_lay_var_720) * (selfrefc_var_282(inds_var_717, ig_var_714) + p_selffrac_var_705(iplon_var_723, i_lay_var_720) * (selfrefc_var_282(inds_var_717 + 1, ig_var_714) - selfrefc_var_282(inds_var_717, ig_var_714))) + p_forfac_var_707(iplon_var_723, i_lay_var_720) * (forrefc_var_283(indf_var_718, ig_var_714) + p_forfrac_var_708(iplon_var_723, i_lay_var_720) * (forrefc_var_283(indf_var_718 + 1, ig_var_714) - forrefc_var_283(indf_var_718, ig_var_714))))
          IF (i_lay_var_720 == i_laysolfr_var_721(iplon_var_723)) p_sfluxzen_var_710(iplon_var_723, ig_var_714) = sfluxrefc_var_284(ig_var_714, js_var_719) + z_fs_var_726 * (sfluxrefc_var_284(ig_var_714, js_var_719 + 1) - sfluxrefc_var_284(ig_var_714, js_var_719))
          p_taur_var_712(iplon_var_723, i_lay_var_720, ig_var_714) = z_tauray_var_730
        END DO
      ELSE
        ind0_var_715 = ((k_jp_var_695(iplon_var_723, i_lay_var_720) - 13) * 5 + (k_jt_var_696(iplon_var_723, i_lay_var_720) - 1)) * nspb_var_314(24) + 1
        ind1_var_716 = ((k_jp_var_695(iplon_var_723, i_lay_var_720) - 12) * 5 + (k_jt1_var_697(iplon_var_723, i_lay_var_720) - 1)) * nspb_var_314(24) + 1
        DO ig_var_714 = 1, 8
          z_tauray_var_730 = p_colmol_var_700(iplon_var_723, i_lay_var_720) * raylbc(ig_var_714)
          p_taug_var_711(iplon_var_723, i_lay_var_720, ig_var_714) = p_colo2_var_701(iplon_var_723, i_lay_var_720) * (p_fac00_var_691(iplon_var_723, i_lay_var_720) * absb_var_281(ind0_var_715, ig_var_714) + p_fac10_var_693(iplon_var_723, i_lay_var_720) * absb_var_281(ind0_var_715 + 1, ig_var_714) + p_fac01_var_692(iplon_var_723, i_lay_var_720) * absb_var_281(ind1_var_716, ig_var_714) + p_fac11_var_694(iplon_var_723, i_lay_var_720) * absb_var_281(ind1_var_716 + 1, ig_var_714)) + p_colo3_var_702(iplon_var_723, i_lay_var_720) * abso3bc_var_286(ig_var_714)
          p_taur_var_712(iplon_var_723, i_lay_var_720, ig_var_714) = z_tauray_var_730
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_720 = laytrop_max_var_725 + 1, i_nlayers_var_722
    DO iplon_var_723 = kidia_var_688, kfdia_var_689
      ind0_var_715 = ((k_jp_var_695(iplon_var_723, i_lay_var_720) - 13) * 5 + (k_jt_var_696(iplon_var_723, i_lay_var_720) - 1)) * nspb_var_314(24) + 1
      ind1_var_716 = ((k_jp_var_695(iplon_var_723, i_lay_var_720) - 12) * 5 + (k_jt1_var_697(iplon_var_723, i_lay_var_720) - 1)) * nspb_var_314(24) + 1
      DO ig_var_714 = 1, 8
        z_tauray_var_730 = p_colmol_var_700(iplon_var_723, i_lay_var_720) * raylbc(ig_var_714)
        p_taug_var_711(iplon_var_723, i_lay_var_720, ig_var_714) = p_colo2_var_701(iplon_var_723, i_lay_var_720) * (p_fac00_var_691(iplon_var_723, i_lay_var_720) * absb_var_281(ind0_var_715, ig_var_714) + p_fac10_var_693(iplon_var_723, i_lay_var_720) * absb_var_281(ind0_var_715 + 1, ig_var_714) + p_fac01_var_692(iplon_var_723, i_lay_var_720) * absb_var_281(ind1_var_716, ig_var_714) + p_fac11_var_694(iplon_var_723, i_lay_var_720) * absb_var_281(ind1_var_716 + 1, ig_var_714)) + p_colo3_var_702(iplon_var_723, i_lay_var_720) * abso3bc_var_286(ig_var_714)
        p_taur_var_712(iplon_var_723, i_lay_var_720, ig_var_714) = z_tauray_var_730
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol24
SUBROUTINE srtm_taumol25(kidia_var_731, kfdia_var_732, klev_var_733, p_fac00_var_734, p_fac01_var_735, p_fac10_var_736, p_fac11_var_737, k_jp_var_738, k_jt_var_739, k_jt1_var_740, p_colh2o_var_741, p_colmol_var_742, p_colo3_var_743, k_laytrop_var_744, p_sfluxzen_var_745, p_taug_var_746, p_taur_var_747, prmu0_var_748)
  USE yoesrta25, ONLY: absa_var_288, abso3ac_var_291, abso3bc_var_292, layreffr_var_287, raylc_var_290, sfluxrefc_var_289
  USE yoesrtwn, ONLY: nspa_var_313
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_731, kfdia_var_732
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_733
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_734(kidia_var_731 : kfdia_var_732, klev_var_733)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_735(kidia_var_731 : kfdia_var_732, klev_var_733)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_736(kidia_var_731 : kfdia_var_732, klev_var_733)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_737(kidia_var_731 : kfdia_var_732, klev_var_733)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_738(kidia_var_731 : kfdia_var_732, klev_var_733)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_739(kidia_var_731 : kfdia_var_732, klev_var_733)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_740(kidia_var_731 : kfdia_var_732, klev_var_733)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_741(kidia_var_731 : kfdia_var_732, klev_var_733)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_742(kidia_var_731 : kfdia_var_732, klev_var_733)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_743(kidia_var_731 : kfdia_var_732, klev_var_733)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_744(kidia_var_731 : kfdia_var_732)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_745(kidia_var_731 : kfdia_var_732, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_746(kidia_var_731 : kfdia_var_732, klev_var_733, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_747(kidia_var_731 : kfdia_var_732, klev_var_733, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_748(kidia_var_731 : kfdia_var_732)
  INTEGER(KIND = 4) :: ig_var_749, ind0_var_750, ind1_var_751, i_lay_var_752, i_laysolfr_var_753(kidia_var_731 : kfdia_var_732), i_nlayers_var_754, iplon_var_755
  INTEGER(KIND = 4) :: laytrop_min_var_756, laytrop_max_var_757
  REAL(KIND = 8) :: z_tauray_var_758
  laytrop_min_var_756 = MINVAL(k_laytrop_var_744(kidia_var_731 : kfdia_var_732))
  laytrop_max_var_757 = MAXVAL(k_laytrop_var_744(kidia_var_731 : kfdia_var_732))
  i_nlayers_var_754 = klev_var_733
  DO iplon_var_755 = kidia_var_731, kfdia_var_732
    i_laysolfr_var_753(iplon_var_755) = k_laytrop_var_744(iplon_var_755)
  END DO
  DO i_lay_var_752 = 1, laytrop_min_var_756
    DO iplon_var_755 = kidia_var_731, kfdia_var_732
      IF (k_jp_var_738(iplon_var_755, i_lay_var_752) < layreffr_var_287 .AND. k_jp_var_738(iplon_var_755, i_lay_var_752 + 1) >= layreffr_var_287) i_laysolfr_var_753(iplon_var_755) = MIN(i_lay_var_752 + 1, k_laytrop_var_744(iplon_var_755))
      ind0_var_750 = ((k_jp_var_738(iplon_var_755, i_lay_var_752) - 1) * 5 + (k_jt_var_739(iplon_var_755, i_lay_var_752) - 1)) * nspa_var_313(25) + 1
      ind1_var_751 = (k_jp_var_738(iplon_var_755, i_lay_var_752) * 5 + (k_jt1_var_740(iplon_var_755, i_lay_var_752) - 1)) * nspa_var_313(25) + 1
      DO ig_var_749 = 1, 6
        z_tauray_var_758 = p_colmol_var_742(iplon_var_755, i_lay_var_752) * raylc_var_290(ig_var_749)
        p_taug_var_746(iplon_var_755, i_lay_var_752, ig_var_749) = p_colh2o_var_741(iplon_var_755, i_lay_var_752) * (p_fac00_var_734(iplon_var_755, i_lay_var_752) * absa_var_288(ind0_var_750, ig_var_749) + p_fac10_var_736(iplon_var_755, i_lay_var_752) * absa_var_288(ind0_var_750 + 1, ig_var_749) + p_fac01_var_735(iplon_var_755, i_lay_var_752) * absa_var_288(ind1_var_751, ig_var_749) + p_fac11_var_737(iplon_var_755, i_lay_var_752) * absa_var_288(ind1_var_751 + 1, ig_var_749)) + p_colo3_var_743(iplon_var_755, i_lay_var_752) * abso3ac_var_291(ig_var_749)
        IF (i_lay_var_752 == i_laysolfr_var_753(iplon_var_755)) p_sfluxzen_var_745(iplon_var_755, ig_var_749) = sfluxrefc_var_289(ig_var_749)
        p_taur_var_747(iplon_var_755, i_lay_var_752, ig_var_749) = z_tauray_var_758
      END DO
    END DO
  END DO
  DO i_lay_var_752 = laytrop_min_var_756 + 1, laytrop_max_var_757
    DO iplon_var_755 = kidia_var_731, kfdia_var_732
      IF (i_lay_var_752 <= k_laytrop_var_744(iplon_var_755)) THEN
        IF (k_jp_var_738(iplon_var_755, i_lay_var_752) < layreffr_var_287 .AND. k_jp_var_738(iplon_var_755, i_lay_var_752 + 1) >= layreffr_var_287) i_laysolfr_var_753(iplon_var_755) = MIN(i_lay_var_752 + 1, k_laytrop_var_744(iplon_var_755))
        ind0_var_750 = ((k_jp_var_738(iplon_var_755, i_lay_var_752) - 1) * 5 + (k_jt_var_739(iplon_var_755, i_lay_var_752) - 1)) * nspa_var_313(25) + 1
        ind1_var_751 = (k_jp_var_738(iplon_var_755, i_lay_var_752) * 5 + (k_jt1_var_740(iplon_var_755, i_lay_var_752) - 1)) * nspa_var_313(25) + 1
        DO ig_var_749 = 1, 6
          z_tauray_var_758 = p_colmol_var_742(iplon_var_755, i_lay_var_752) * raylc_var_290(ig_var_749)
          p_taug_var_746(iplon_var_755, i_lay_var_752, ig_var_749) = p_colh2o_var_741(iplon_var_755, i_lay_var_752) * (p_fac00_var_734(iplon_var_755, i_lay_var_752) * absa_var_288(ind0_var_750, ig_var_749) + p_fac10_var_736(iplon_var_755, i_lay_var_752) * absa_var_288(ind0_var_750 + 1, ig_var_749) + p_fac01_var_735(iplon_var_755, i_lay_var_752) * absa_var_288(ind1_var_751, ig_var_749) + p_fac11_var_737(iplon_var_755, i_lay_var_752) * absa_var_288(ind1_var_751 + 1, ig_var_749)) + p_colo3_var_743(iplon_var_755, i_lay_var_752) * abso3ac_var_291(ig_var_749)
          IF (i_lay_var_752 == i_laysolfr_var_753(iplon_var_755)) p_sfluxzen_var_745(iplon_var_755, ig_var_749) = sfluxrefc_var_289(ig_var_749)
          p_taur_var_747(iplon_var_755, i_lay_var_752, ig_var_749) = z_tauray_var_758
        END DO
      ELSE
        DO ig_var_749 = 1, 6
          z_tauray_var_758 = p_colmol_var_742(iplon_var_755, i_lay_var_752) * raylc_var_290(ig_var_749)
          p_taug_var_746(iplon_var_755, i_lay_var_752, ig_var_749) = p_colo3_var_743(iplon_var_755, i_lay_var_752) * abso3bc_var_292(ig_var_749)
          p_taur_var_747(iplon_var_755, i_lay_var_752, ig_var_749) = z_tauray_var_758
        END DO
      END IF
    END DO
  END DO
  DO ig_var_749 = 1, 6
    DO i_lay_var_752 = laytrop_max_var_757 + 1, i_nlayers_var_754
      DO iplon_var_755 = kidia_var_731, kfdia_var_732
        z_tauray_var_758 = p_colmol_var_742(iplon_var_755, i_lay_var_752) * raylc_var_290(ig_var_749)
        p_taug_var_746(iplon_var_755, i_lay_var_752, ig_var_749) = p_colo3_var_743(iplon_var_755, i_lay_var_752) * abso3bc_var_292(ig_var_749)
        p_taur_var_747(iplon_var_755, i_lay_var_752, ig_var_749) = z_tauray_var_758
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol25
SUBROUTINE srtm_taumol19(kidia_var_759, kfdia_var_760, klev_var_761, p_fac00_var_762, p_fac01_var_763, p_fac10_var_764, p_fac11_var_765, k_jp_var_766, k_jt_var_767, k_jt1_var_768, p_oneminus_var_769, p_colh2o_var_770, p_colco2_var_771, p_colmol_var_772, k_laytrop_var_773, p_selffac_var_774, p_selffrac_var_775, k_indself_var_776, p_forfac_var_777, p_forfrac_var_778, k_indfor_var_779, p_sfluxzen_var_780, p_taug_var_781, p_taur_var_782, prmu0_var_783)
  USE yoesrta19, ONLY: absa_var_244, absb_var_245, forrefc_var_247, layreffr_var_243, rayl_var_241, selfrefc_var_246, sfluxrefc_var_248, strrat_var_242
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_759, kfdia_var_760
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_761
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_762(kidia_var_759 : kfdia_var_760, klev_var_761)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_763(kidia_var_759 : kfdia_var_760, klev_var_761)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_764(kidia_var_759 : kfdia_var_760, klev_var_761)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_765(kidia_var_759 : kfdia_var_760, klev_var_761)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_766(kidia_var_759 : kfdia_var_760, klev_var_761)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_767(kidia_var_759 : kfdia_var_760, klev_var_761)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_768(kidia_var_759 : kfdia_var_760, klev_var_761)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_769(kidia_var_759 : kfdia_var_760)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_770(kidia_var_759 : kfdia_var_760, klev_var_761)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_771(kidia_var_759 : kfdia_var_760, klev_var_761)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_772(kidia_var_759 : kfdia_var_760, klev_var_761)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_773(kidia_var_759 : kfdia_var_760)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_774(kidia_var_759 : kfdia_var_760, klev_var_761)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_775(kidia_var_759 : kfdia_var_760, klev_var_761)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_776(kidia_var_759 : kfdia_var_760, klev_var_761)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_777(kidia_var_759 : kfdia_var_760, klev_var_761)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_778(kidia_var_759 : kfdia_var_760, klev_var_761)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_779(kidia_var_759 : kfdia_var_760, klev_var_761)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_780(kidia_var_759 : kfdia_var_760, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_781(kidia_var_759 : kfdia_var_760, klev_var_761, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_782(kidia_var_759 : kfdia_var_760, klev_var_761, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_783(kidia_var_759 : kfdia_var_760)
  INTEGER(KIND = 4) :: ig_var_784, ind0_var_785, ind1_var_786, inds_var_787, indf_var_788, js_var_789, i_lay_var_790, i_laysolfr_var_791(kidia_var_759 : kfdia_var_760), i_nlayers_var_792, iplon_var_793
  INTEGER(KIND = 4) :: laytrop_min_var_794, laytrop_max_var_795
  REAL(KIND = 8) :: z_fs_var_796, z_speccomb_var_797, z_specmult_var_798, z_specparm_var_799, z_tauray_var_800
  laytrop_min_var_794 = MINVAL(k_laytrop_var_773(kidia_var_759 : kfdia_var_760))
  laytrop_max_var_795 = MAXVAL(k_laytrop_var_773(kidia_var_759 : kfdia_var_760))
  i_nlayers_var_792 = klev_var_761
  DO iplon_var_793 = kidia_var_759, kfdia_var_760
    i_laysolfr_var_791(iplon_var_793) = k_laytrop_var_773(iplon_var_793)
  END DO
  DO i_lay_var_790 = 1, laytrop_min_var_794
    DO iplon_var_793 = kidia_var_759, kfdia_var_760
      IF (k_jp_var_766(iplon_var_793, i_lay_var_790) < layreffr_var_243 .AND. k_jp_var_766(iplon_var_793, i_lay_var_790 + 1) >= layreffr_var_243) i_laysolfr_var_791(iplon_var_793) = MIN(i_lay_var_790 + 1, k_laytrop_var_773(iplon_var_793))
      z_speccomb_var_797 = p_colh2o_var_770(iplon_var_793, i_lay_var_790) + strrat_var_242 * p_colco2_var_771(iplon_var_793, i_lay_var_790)
      z_specparm_var_799 = p_colh2o_var_770(iplon_var_793, i_lay_var_790) / z_speccomb_var_797
      z_specparm_var_799 = MIN(p_oneminus_var_769(iplon_var_793), z_specparm_var_799)
      z_specmult_var_798 = 8.0D0 * (z_specparm_var_799)
      js_var_789 = 1 + INT(z_specmult_var_798)
      z_fs_var_796 = z_specmult_var_798 - AINT(z_specmult_var_798)
      ind0_var_785 = ((k_jp_var_766(iplon_var_793, i_lay_var_790) - 1) * 5 + (k_jt_var_767(iplon_var_793, i_lay_var_790) - 1)) * nspa_var_313(19) + js_var_789
      ind1_var_786 = (k_jp_var_766(iplon_var_793, i_lay_var_790) * 5 + (k_jt1_var_768(iplon_var_793, i_lay_var_790) - 1)) * nspa_var_313(19) + js_var_789
      inds_var_787 = k_indself_var_776(iplon_var_793, i_lay_var_790)
      indf_var_788 = k_indfor_var_779(iplon_var_793, i_lay_var_790)
      z_tauray_var_800 = p_colmol_var_772(iplon_var_793, i_lay_var_790) * rayl_var_241
      DO ig_var_784 = 1, 8
        p_taug_var_781(iplon_var_793, i_lay_var_790, ig_var_784) = z_speccomb_var_797 * ((1.0D0 - z_fs_var_796) * (absa_var_244(ind0_var_785, ig_var_784) * p_fac00_var_762(iplon_var_793, i_lay_var_790) + absa_var_244(ind0_var_785 + 9, ig_var_784) * p_fac10_var_764(iplon_var_793, i_lay_var_790) + absa_var_244(ind1_var_786, ig_var_784) * p_fac01_var_763(iplon_var_793, i_lay_var_790) + absa_var_244(ind1_var_786 + 9, ig_var_784) * p_fac11_var_765(iplon_var_793, i_lay_var_790)) + z_fs_var_796 * (absa_var_244(ind0_var_785 + 1, ig_var_784) * p_fac00_var_762(iplon_var_793, i_lay_var_790) + absa_var_244(ind0_var_785 + 10, ig_var_784) * p_fac10_var_764(iplon_var_793, i_lay_var_790) + absa_var_244(ind1_var_786 + 1, ig_var_784) * p_fac01_var_763(iplon_var_793, i_lay_var_790) + absa_var_244(ind1_var_786 + 10, ig_var_784) * p_fac11_var_765(iplon_var_793, i_lay_var_790))) + p_colh2o_var_770(iplon_var_793, i_lay_var_790) * (p_selffac_var_774(iplon_var_793, i_lay_var_790) * (selfrefc_var_246(inds_var_787, ig_var_784) + p_selffrac_var_775(iplon_var_793, i_lay_var_790) * (selfrefc_var_246(inds_var_787 + 1, ig_var_784) - selfrefc_var_246(inds_var_787, ig_var_784))) + p_forfac_var_777(iplon_var_793, i_lay_var_790) * (forrefc_var_247(indf_var_788, ig_var_784) + p_forfrac_var_778(iplon_var_793, i_lay_var_790) * (forrefc_var_247(indf_var_788 + 1, ig_var_784) - forrefc_var_247(indf_var_788, ig_var_784))))
        IF (i_lay_var_790 == i_laysolfr_var_791(iplon_var_793)) p_sfluxzen_var_780(iplon_var_793, ig_var_784) = sfluxrefc_var_248(ig_var_784, js_var_789) + z_fs_var_796 * (sfluxrefc_var_248(ig_var_784, js_var_789 + 1) - sfluxrefc_var_248(ig_var_784, js_var_789))
        p_taur_var_782(iplon_var_793, i_lay_var_790, ig_var_784) = z_tauray_var_800
      END DO
    END DO
  END DO
  DO i_lay_var_790 = laytrop_min_var_794 + 1, laytrop_max_var_795
    DO iplon_var_793 = kidia_var_759, kfdia_var_760
      IF (i_lay_var_790 <= k_laytrop_var_773(iplon_var_793)) THEN
        IF (k_jp_var_766(iplon_var_793, i_lay_var_790) < layreffr_var_243 .AND. k_jp_var_766(iplon_var_793, i_lay_var_790 + 1) >= layreffr_var_243) i_laysolfr_var_791(iplon_var_793) = MIN(i_lay_var_790 + 1, k_laytrop_var_773(iplon_var_793))
        z_speccomb_var_797 = p_colh2o_var_770(iplon_var_793, i_lay_var_790) + strrat_var_242 * p_colco2_var_771(iplon_var_793, i_lay_var_790)
        z_specparm_var_799 = p_colh2o_var_770(iplon_var_793, i_lay_var_790) / z_speccomb_var_797
        z_specparm_var_799 = MIN(p_oneminus_var_769(iplon_var_793), z_specparm_var_799)
        z_specmult_var_798 = 8.0D0 * (z_specparm_var_799)
        js_var_789 = 1 + INT(z_specmult_var_798)
        z_fs_var_796 = z_specmult_var_798 - AINT(z_specmult_var_798)
        ind0_var_785 = ((k_jp_var_766(iplon_var_793, i_lay_var_790) - 1) * 5 + (k_jt_var_767(iplon_var_793, i_lay_var_790) - 1)) * nspa_var_313(19) + js_var_789
        ind1_var_786 = (k_jp_var_766(iplon_var_793, i_lay_var_790) * 5 + (k_jt1_var_768(iplon_var_793, i_lay_var_790) - 1)) * nspa_var_313(19) + js_var_789
        inds_var_787 = k_indself_var_776(iplon_var_793, i_lay_var_790)
        indf_var_788 = k_indfor_var_779(iplon_var_793, i_lay_var_790)
        z_tauray_var_800 = p_colmol_var_772(iplon_var_793, i_lay_var_790) * rayl_var_241
        DO ig_var_784 = 1, 8
          p_taug_var_781(iplon_var_793, i_lay_var_790, ig_var_784) = z_speccomb_var_797 * ((1.0D0 - z_fs_var_796) * (absa_var_244(ind0_var_785, ig_var_784) * p_fac00_var_762(iplon_var_793, i_lay_var_790) + absa_var_244(ind0_var_785 + 9, ig_var_784) * p_fac10_var_764(iplon_var_793, i_lay_var_790) + absa_var_244(ind1_var_786, ig_var_784) * p_fac01_var_763(iplon_var_793, i_lay_var_790) + absa_var_244(ind1_var_786 + 9, ig_var_784) * p_fac11_var_765(iplon_var_793, i_lay_var_790)) + z_fs_var_796 * (absa_var_244(ind0_var_785 + 1, ig_var_784) * p_fac00_var_762(iplon_var_793, i_lay_var_790) + absa_var_244(ind0_var_785 + 10, ig_var_784) * p_fac10_var_764(iplon_var_793, i_lay_var_790) + absa_var_244(ind1_var_786 + 1, ig_var_784) * p_fac01_var_763(iplon_var_793, i_lay_var_790) + absa_var_244(ind1_var_786 + 10, ig_var_784) * p_fac11_var_765(iplon_var_793, i_lay_var_790))) + p_colh2o_var_770(iplon_var_793, i_lay_var_790) * (p_selffac_var_774(iplon_var_793, i_lay_var_790) * (selfrefc_var_246(inds_var_787, ig_var_784) + p_selffrac_var_775(iplon_var_793, i_lay_var_790) * (selfrefc_var_246(inds_var_787 + 1, ig_var_784) - selfrefc_var_246(inds_var_787, ig_var_784))) + p_forfac_var_777(iplon_var_793, i_lay_var_790) * (forrefc_var_247(indf_var_788, ig_var_784) + p_forfrac_var_778(iplon_var_793, i_lay_var_790) * (forrefc_var_247(indf_var_788 + 1, ig_var_784) - forrefc_var_247(indf_var_788, ig_var_784))))
          IF (i_lay_var_790 == i_laysolfr_var_791(iplon_var_793)) p_sfluxzen_var_780(iplon_var_793, ig_var_784) = sfluxrefc_var_248(ig_var_784, js_var_789) + z_fs_var_796 * (sfluxrefc_var_248(ig_var_784, js_var_789 + 1) - sfluxrefc_var_248(ig_var_784, js_var_789))
          p_taur_var_782(iplon_var_793, i_lay_var_790, ig_var_784) = z_tauray_var_800
        END DO
      ELSE
        ind0_var_785 = ((k_jp_var_766(iplon_var_793, i_lay_var_790) - 13) * 5 + (k_jt_var_767(iplon_var_793, i_lay_var_790) - 1)) * nspb_var_314(19) + 1
        ind1_var_786 = ((k_jp_var_766(iplon_var_793, i_lay_var_790) - 12) * 5 + (k_jt1_var_768(iplon_var_793, i_lay_var_790) - 1)) * nspb_var_314(19) + 1
        z_tauray_var_800 = p_colmol_var_772(iplon_var_793, i_lay_var_790) * rayl_var_241
        DO ig_var_784 = 1, 8
          p_taug_var_781(iplon_var_793, i_lay_var_790, ig_var_784) = p_colco2_var_771(iplon_var_793, i_lay_var_790) * (p_fac00_var_762(iplon_var_793, i_lay_var_790) * absb_var_245(ind0_var_785, ig_var_784) + p_fac10_var_764(iplon_var_793, i_lay_var_790) * absb_var_245(ind0_var_785 + 1, ig_var_784) + p_fac01_var_763(iplon_var_793, i_lay_var_790) * absb_var_245(ind1_var_786, ig_var_784) + p_fac11_var_765(iplon_var_793, i_lay_var_790) * absb_var_245(ind1_var_786 + 1, ig_var_784))
          p_taur_var_782(iplon_var_793, i_lay_var_790, ig_var_784) = z_tauray_var_800
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_790 = laytrop_max_var_795 + 1, i_nlayers_var_792
    DO iplon_var_793 = kidia_var_759, kfdia_var_760
      ind0_var_785 = ((k_jp_var_766(iplon_var_793, i_lay_var_790) - 13) * 5 + (k_jt_var_767(iplon_var_793, i_lay_var_790) - 1)) * nspb_var_314(19) + 1
      ind1_var_786 = ((k_jp_var_766(iplon_var_793, i_lay_var_790) - 12) * 5 + (k_jt1_var_768(iplon_var_793, i_lay_var_790) - 1)) * nspb_var_314(19) + 1
      z_tauray_var_800 = p_colmol_var_772(iplon_var_793, i_lay_var_790) * rayl_var_241
      DO ig_var_784 = 1, 8
        p_taug_var_781(iplon_var_793, i_lay_var_790, ig_var_784) = p_colco2_var_771(iplon_var_793, i_lay_var_790) * (p_fac00_var_762(iplon_var_793, i_lay_var_790) * absb_var_245(ind0_var_785, ig_var_784) + p_fac10_var_764(iplon_var_793, i_lay_var_790) * absb_var_245(ind0_var_785 + 1, ig_var_784) + p_fac01_var_763(iplon_var_793, i_lay_var_790) * absb_var_245(ind1_var_786, ig_var_784) + p_fac11_var_765(iplon_var_793, i_lay_var_790) * absb_var_245(ind1_var_786 + 1, ig_var_784))
        p_taur_var_782(iplon_var_793, i_lay_var_790, ig_var_784) = z_tauray_var_800
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol19
SUBROUTINE srtm_taumol21(kidia_var_801, kfdia_var_802, klev_var_803, p_fac00_var_804, p_fac01_var_805, p_fac10_var_806, p_fac11_var_807, k_jp_var_808, k_jt_var_809, k_jt1_var_810, p_oneminus_var_811, p_colh2o_var_812, p_colco2_var_813, p_colmol_var_814, k_laytrop_var_815, p_selffac_var_816, p_selffrac_var_817, k_indself_var_818, p_forfac_var_819, p_forfrac_var_820, k_indfor_var_821, p_sfluxzen_var_822, p_taug_var_823, p_taur_var_824, prmu0_var_825)
  USE yoesrta21, ONLY: absa_var_259, absb_var_260, forrefc_var_262, layreffr_var_258, rayl_var_256, selfrefc_var_261, sfluxrefc_var_263, strrat_var_257
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_801, kfdia_var_802
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_803
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_804(kidia_var_801 : kfdia_var_802, klev_var_803)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_805(kidia_var_801 : kfdia_var_802, klev_var_803)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_806(kidia_var_801 : kfdia_var_802, klev_var_803)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_807(kidia_var_801 : kfdia_var_802, klev_var_803)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_808(kidia_var_801 : kfdia_var_802, klev_var_803)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_809(kidia_var_801 : kfdia_var_802, klev_var_803)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_810(kidia_var_801 : kfdia_var_802, klev_var_803)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_811(kidia_var_801 : kfdia_var_802)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_812(kidia_var_801 : kfdia_var_802, klev_var_803)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_813(kidia_var_801 : kfdia_var_802, klev_var_803)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_814(kidia_var_801 : kfdia_var_802, klev_var_803)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_815(kidia_var_801 : kfdia_var_802)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_816(kidia_var_801 : kfdia_var_802, klev_var_803)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_817(kidia_var_801 : kfdia_var_802, klev_var_803)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_818(kidia_var_801 : kfdia_var_802, klev_var_803)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_819(kidia_var_801 : kfdia_var_802, klev_var_803)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_820(kidia_var_801 : kfdia_var_802, klev_var_803)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_821(kidia_var_801 : kfdia_var_802, klev_var_803)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_822(kidia_var_801 : kfdia_var_802, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_823(kidia_var_801 : kfdia_var_802, klev_var_803, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_824(kidia_var_801 : kfdia_var_802, klev_var_803, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_825(kidia_var_801 : kfdia_var_802)
  INTEGER(KIND = 4) :: ig_var_826, ind0_var_827, ind1_var_828, inds_var_829, indf_var_830, js_var_831, i_lay_var_832, i_laysolfr_var_833(kidia_var_801 : kfdia_var_802), i_nlayers_var_834, iplon_var_835
  INTEGER(KIND = 4) :: laytrop_min_var_836, laytrop_max_var_837
  REAL(KIND = 8) :: z_fs_var_838, z_speccomb_var_839, z_specmult_var_840, z_specparm_var_841, z_tauray_var_842
  laytrop_min_var_836 = MINVAL(k_laytrop_var_815(kidia_var_801 : kfdia_var_802))
  laytrop_max_var_837 = MAXVAL(k_laytrop_var_815(kidia_var_801 : kfdia_var_802))
  i_nlayers_var_834 = klev_var_803
  DO iplon_var_835 = kidia_var_801, kfdia_var_802
    i_laysolfr_var_833(iplon_var_835) = k_laytrop_var_815(iplon_var_835)
  END DO
  DO i_lay_var_832 = 1, laytrop_min_var_836
    DO iplon_var_835 = kidia_var_801, kfdia_var_802
      IF (k_jp_var_808(iplon_var_835, i_lay_var_832) < layreffr_var_258 .AND. k_jp_var_808(iplon_var_835, i_lay_var_832 + 1) >= layreffr_var_258) i_laysolfr_var_833(iplon_var_835) = MIN(i_lay_var_832 + 1, k_laytrop_var_815(iplon_var_835))
      z_speccomb_var_839 = p_colh2o_var_812(iplon_var_835, i_lay_var_832) + strrat_var_257 * p_colco2_var_813(iplon_var_835, i_lay_var_832)
      z_specparm_var_841 = p_colh2o_var_812(iplon_var_835, i_lay_var_832) / z_speccomb_var_839
      z_specparm_var_841 = MIN(p_oneminus_var_811(iplon_var_835), z_specparm_var_841)
      z_specmult_var_840 = 8.0D0 * (z_specparm_var_841)
      js_var_831 = 1 + INT(z_specmult_var_840)
      z_fs_var_838 = z_specmult_var_840 - AINT(z_specmult_var_840)
      ind0_var_827 = ((k_jp_var_808(iplon_var_835, i_lay_var_832) - 1) * 5 + (k_jt_var_809(iplon_var_835, i_lay_var_832) - 1)) * nspa_var_313(21) + js_var_831
      ind1_var_828 = (k_jp_var_808(iplon_var_835, i_lay_var_832) * 5 + (k_jt1_var_810(iplon_var_835, i_lay_var_832) - 1)) * nspa_var_313(21) + js_var_831
      inds_var_829 = k_indself_var_818(iplon_var_835, i_lay_var_832)
      indf_var_830 = k_indfor_var_821(iplon_var_835, i_lay_var_832)
      z_tauray_var_842 = p_colmol_var_814(iplon_var_835, i_lay_var_832) * rayl_var_256
      DO ig_var_826 = 1, 10
        p_taug_var_823(iplon_var_835, i_lay_var_832, ig_var_826) = z_speccomb_var_839 * ((1.0D0 - z_fs_var_838) * (absa_var_259(ind0_var_827, ig_var_826) * p_fac00_var_804(iplon_var_835, i_lay_var_832) + absa_var_259(ind0_var_827 + 9, ig_var_826) * p_fac10_var_806(iplon_var_835, i_lay_var_832) + absa_var_259(ind1_var_828, ig_var_826) * p_fac01_var_805(iplon_var_835, i_lay_var_832) + absa_var_259(ind1_var_828 + 9, ig_var_826) * p_fac11_var_807(iplon_var_835, i_lay_var_832)) + z_fs_var_838 * (absa_var_259(ind0_var_827 + 1, ig_var_826) * p_fac00_var_804(iplon_var_835, i_lay_var_832) + absa_var_259(ind0_var_827 + 10, ig_var_826) * p_fac10_var_806(iplon_var_835, i_lay_var_832) + absa_var_259(ind1_var_828 + 1, ig_var_826) * p_fac01_var_805(iplon_var_835, i_lay_var_832) + absa_var_259(ind1_var_828 + 10, ig_var_826) * p_fac11_var_807(iplon_var_835, i_lay_var_832))) + p_colh2o_var_812(iplon_var_835, i_lay_var_832) * (p_selffac_var_816(iplon_var_835, i_lay_var_832) * (selfrefc_var_261(inds_var_829, ig_var_826) + p_selffrac_var_817(iplon_var_835, i_lay_var_832) * (selfrefc_var_261(inds_var_829 + 1, ig_var_826) - selfrefc_var_261(inds_var_829, ig_var_826))) + p_forfac_var_819(iplon_var_835, i_lay_var_832) * (forrefc_var_262(indf_var_830, ig_var_826) + p_forfrac_var_820(iplon_var_835, i_lay_var_832) * (forrefc_var_262(indf_var_830 + 1, ig_var_826) - forrefc_var_262(indf_var_830, ig_var_826))))
        IF (i_lay_var_832 == i_laysolfr_var_833(iplon_var_835)) p_sfluxzen_var_822(iplon_var_835, ig_var_826) = sfluxrefc_var_263(ig_var_826, js_var_831) + z_fs_var_838 * (sfluxrefc_var_263(ig_var_826, js_var_831 + 1) - sfluxrefc_var_263(ig_var_826, js_var_831))
        p_taur_var_824(iplon_var_835, i_lay_var_832, ig_var_826) = z_tauray_var_842
      END DO
    END DO
  END DO
  DO i_lay_var_832 = laytrop_min_var_836 + 1, laytrop_max_var_837
    DO iplon_var_835 = kidia_var_801, kfdia_var_802
      IF (i_lay_var_832 <= k_laytrop_var_815(iplon_var_835)) THEN
        IF (k_jp_var_808(iplon_var_835, i_lay_var_832) < layreffr_var_258 .AND. k_jp_var_808(iplon_var_835, i_lay_var_832 + 1) >= layreffr_var_258) i_laysolfr_var_833(iplon_var_835) = MIN(i_lay_var_832 + 1, k_laytrop_var_815(iplon_var_835))
        z_speccomb_var_839 = p_colh2o_var_812(iplon_var_835, i_lay_var_832) + strrat_var_257 * p_colco2_var_813(iplon_var_835, i_lay_var_832)
        z_specparm_var_841 = p_colh2o_var_812(iplon_var_835, i_lay_var_832) / z_speccomb_var_839
        z_specparm_var_841 = MIN(p_oneminus_var_811(iplon_var_835), z_specparm_var_841)
        z_specmult_var_840 = 8.0D0 * (z_specparm_var_841)
        js_var_831 = 1 + INT(z_specmult_var_840)
        z_fs_var_838 = z_specmult_var_840 - AINT(z_specmult_var_840)
        ind0_var_827 = ((k_jp_var_808(iplon_var_835, i_lay_var_832) - 1) * 5 + (k_jt_var_809(iplon_var_835, i_lay_var_832) - 1)) * nspa_var_313(21) + js_var_831
        ind1_var_828 = (k_jp_var_808(iplon_var_835, i_lay_var_832) * 5 + (k_jt1_var_810(iplon_var_835, i_lay_var_832) - 1)) * nspa_var_313(21) + js_var_831
        inds_var_829 = k_indself_var_818(iplon_var_835, i_lay_var_832)
        indf_var_830 = k_indfor_var_821(iplon_var_835, i_lay_var_832)
        z_tauray_var_842 = p_colmol_var_814(iplon_var_835, i_lay_var_832) * rayl_var_256
        DO ig_var_826 = 1, 10
          p_taug_var_823(iplon_var_835, i_lay_var_832, ig_var_826) = z_speccomb_var_839 * ((1.0D0 - z_fs_var_838) * (absa_var_259(ind0_var_827, ig_var_826) * p_fac00_var_804(iplon_var_835, i_lay_var_832) + absa_var_259(ind0_var_827 + 9, ig_var_826) * p_fac10_var_806(iplon_var_835, i_lay_var_832) + absa_var_259(ind1_var_828, ig_var_826) * p_fac01_var_805(iplon_var_835, i_lay_var_832) + absa_var_259(ind1_var_828 + 9, ig_var_826) * p_fac11_var_807(iplon_var_835, i_lay_var_832)) + z_fs_var_838 * (absa_var_259(ind0_var_827 + 1, ig_var_826) * p_fac00_var_804(iplon_var_835, i_lay_var_832) + absa_var_259(ind0_var_827 + 10, ig_var_826) * p_fac10_var_806(iplon_var_835, i_lay_var_832) + absa_var_259(ind1_var_828 + 1, ig_var_826) * p_fac01_var_805(iplon_var_835, i_lay_var_832) + absa_var_259(ind1_var_828 + 10, ig_var_826) * p_fac11_var_807(iplon_var_835, i_lay_var_832))) + p_colh2o_var_812(iplon_var_835, i_lay_var_832) * (p_selffac_var_816(iplon_var_835, i_lay_var_832) * (selfrefc_var_261(inds_var_829, ig_var_826) + p_selffrac_var_817(iplon_var_835, i_lay_var_832) * (selfrefc_var_261(inds_var_829 + 1, ig_var_826) - selfrefc_var_261(inds_var_829, ig_var_826))) + p_forfac_var_819(iplon_var_835, i_lay_var_832) * (forrefc_var_262(indf_var_830, ig_var_826) + p_forfrac_var_820(iplon_var_835, i_lay_var_832) * (forrefc_var_262(indf_var_830 + 1, ig_var_826) - forrefc_var_262(indf_var_830, ig_var_826))))
          IF (i_lay_var_832 == i_laysolfr_var_833(iplon_var_835)) p_sfluxzen_var_822(iplon_var_835, ig_var_826) = sfluxrefc_var_263(ig_var_826, js_var_831) + z_fs_var_838 * (sfluxrefc_var_263(ig_var_826, js_var_831 + 1) - sfluxrefc_var_263(ig_var_826, js_var_831))
          p_taur_var_824(iplon_var_835, i_lay_var_832, ig_var_826) = z_tauray_var_842
        END DO
      ELSE
        z_speccomb_var_839 = p_colh2o_var_812(iplon_var_835, i_lay_var_832) + strrat_var_257 * p_colco2_var_813(iplon_var_835, i_lay_var_832)
        z_specparm_var_841 = p_colh2o_var_812(iplon_var_835, i_lay_var_832) / z_speccomb_var_839
        z_specparm_var_841 = MIN(p_oneminus_var_811(iplon_var_835), z_specparm_var_841)
        z_specmult_var_840 = 4.0D0 * (z_specparm_var_841)
        js_var_831 = 1 + INT(z_specmult_var_840)
        z_fs_var_838 = z_specmult_var_840 - AINT(z_specmult_var_840)
        ind0_var_827 = ((k_jp_var_808(iplon_var_835, i_lay_var_832) - 13) * 5 + (k_jt_var_809(iplon_var_835, i_lay_var_832) - 1)) * nspb_var_314(21) + js_var_831
        ind1_var_828 = ((k_jp_var_808(iplon_var_835, i_lay_var_832) - 12) * 5 + (k_jt1_var_810(iplon_var_835, i_lay_var_832) - 1)) * nspb_var_314(21) + js_var_831
        indf_var_830 = k_indfor_var_821(iplon_var_835, i_lay_var_832)
        z_tauray_var_842 = p_colmol_var_814(iplon_var_835, i_lay_var_832) * rayl_var_256
        DO ig_var_826 = 1, 10
          p_taug_var_823(iplon_var_835, i_lay_var_832, ig_var_826) = z_speccomb_var_839 * ((1.0D0 - z_fs_var_838) * (absb_var_260(ind0_var_827, ig_var_826) * p_fac00_var_804(iplon_var_835, i_lay_var_832) + absb_var_260(ind0_var_827 + 5, ig_var_826) * p_fac10_var_806(iplon_var_835, i_lay_var_832) + absb_var_260(ind1_var_828, ig_var_826) * p_fac01_var_805(iplon_var_835, i_lay_var_832) + absb_var_260(ind1_var_828 + 5, ig_var_826) * p_fac11_var_807(iplon_var_835, i_lay_var_832)) + z_fs_var_838 * (absb_var_260(ind0_var_827 + 1, ig_var_826) * p_fac00_var_804(iplon_var_835, i_lay_var_832) + absb_var_260(ind0_var_827 + 6, ig_var_826) * p_fac10_var_806(iplon_var_835, i_lay_var_832) + absb_var_260(ind1_var_828 + 1, ig_var_826) * p_fac01_var_805(iplon_var_835, i_lay_var_832) + absb_var_260(ind1_var_828 + 6, ig_var_826) * p_fac11_var_807(iplon_var_835, i_lay_var_832))) + p_colh2o_var_812(iplon_var_835, i_lay_var_832) * p_forfac_var_819(iplon_var_835, i_lay_var_832) * (forrefc_var_262(indf_var_830, ig_var_826) + p_forfrac_var_820(iplon_var_835, i_lay_var_832) * (forrefc_var_262(indf_var_830 + 1, ig_var_826) - forrefc_var_262(indf_var_830, ig_var_826)))
          p_taur_var_824(iplon_var_835, i_lay_var_832, ig_var_826) = z_tauray_var_842
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_832 = laytrop_max_var_837 + 1, i_nlayers_var_834
    DO iplon_var_835 = kidia_var_801, kfdia_var_802
      z_speccomb_var_839 = p_colh2o_var_812(iplon_var_835, i_lay_var_832) + strrat_var_257 * p_colco2_var_813(iplon_var_835, i_lay_var_832)
      z_specparm_var_841 = p_colh2o_var_812(iplon_var_835, i_lay_var_832) / z_speccomb_var_839
      z_specparm_var_841 = MIN(p_oneminus_var_811(iplon_var_835), z_specparm_var_841)
      z_specmult_var_840 = 4.0D0 * (z_specparm_var_841)
      js_var_831 = 1 + INT(z_specmult_var_840)
      z_fs_var_838 = z_specmult_var_840 - AINT(z_specmult_var_840)
      ind0_var_827 = ((k_jp_var_808(iplon_var_835, i_lay_var_832) - 13) * 5 + (k_jt_var_809(iplon_var_835, i_lay_var_832) - 1)) * nspb_var_314(21) + js_var_831
      ind1_var_828 = ((k_jp_var_808(iplon_var_835, i_lay_var_832) - 12) * 5 + (k_jt1_var_810(iplon_var_835, i_lay_var_832) - 1)) * nspb_var_314(21) + js_var_831
      indf_var_830 = k_indfor_var_821(iplon_var_835, i_lay_var_832)
      z_tauray_var_842 = p_colmol_var_814(iplon_var_835, i_lay_var_832) * rayl_var_256
      DO ig_var_826 = 1, 10
        p_taug_var_823(iplon_var_835, i_lay_var_832, ig_var_826) = z_speccomb_var_839 * ((1.0D0 - z_fs_var_838) * (absb_var_260(ind0_var_827, ig_var_826) * p_fac00_var_804(iplon_var_835, i_lay_var_832) + absb_var_260(ind0_var_827 + 5, ig_var_826) * p_fac10_var_806(iplon_var_835, i_lay_var_832) + absb_var_260(ind1_var_828, ig_var_826) * p_fac01_var_805(iplon_var_835, i_lay_var_832) + absb_var_260(ind1_var_828 + 5, ig_var_826) * p_fac11_var_807(iplon_var_835, i_lay_var_832)) + z_fs_var_838 * (absb_var_260(ind0_var_827 + 1, ig_var_826) * p_fac00_var_804(iplon_var_835, i_lay_var_832) + absb_var_260(ind0_var_827 + 6, ig_var_826) * p_fac10_var_806(iplon_var_835, i_lay_var_832) + absb_var_260(ind1_var_828 + 1, ig_var_826) * p_fac01_var_805(iplon_var_835, i_lay_var_832) + absb_var_260(ind1_var_828 + 6, ig_var_826) * p_fac11_var_807(iplon_var_835, i_lay_var_832))) + p_colh2o_var_812(iplon_var_835, i_lay_var_832) * p_forfac_var_819(iplon_var_835, i_lay_var_832) * (forrefc_var_262(indf_var_830, ig_var_826) + p_forfrac_var_820(iplon_var_835, i_lay_var_832) * (forrefc_var_262(indf_var_830 + 1, ig_var_826) - forrefc_var_262(indf_var_830, ig_var_826)))
        p_taur_var_824(iplon_var_835, i_lay_var_832, ig_var_826) = z_tauray_var_842
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol21
SUBROUTINE srtm_taumol20(kidia_var_843, kfdia_var_844, klev_var_845, p_fac00_var_846, p_fac01_var_847, p_fac10_var_848, p_fac11_var_849, k_jp_var_850, k_jt_var_851, k_jt1_var_852, p_colh2o_var_853, p_colch4_var_854, p_colmol_var_855, k_laytrop_var_856, p_selffac_var_857, p_selffrac_var_858, k_indself_var_859, p_forfac_var_860, p_forfrac_var_861, k_indfor_var_862, p_sfluxzen_var_863, p_taug_var_864, p_taur_var_865, prmu0_var_866)
  USE yoesrta20, ONLY: absa_var_251, absb_var_252, absch4c, forrefc_var_254, layreffr_var_250, rayl_var_249, selfrefc_var_253, sfluxrefc_var_255
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_843, kfdia_var_844
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_845
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_846(kidia_var_843 : kfdia_var_844, klev_var_845)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_847(kidia_var_843 : kfdia_var_844, klev_var_845)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_848(kidia_var_843 : kfdia_var_844, klev_var_845)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_849(kidia_var_843 : kfdia_var_844, klev_var_845)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_850(kidia_var_843 : kfdia_var_844, klev_var_845)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_851(kidia_var_843 : kfdia_var_844, klev_var_845)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_852(kidia_var_843 : kfdia_var_844, klev_var_845)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_853(kidia_var_843 : kfdia_var_844, klev_var_845)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_854(kidia_var_843 : kfdia_var_844, klev_var_845)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_855(kidia_var_843 : kfdia_var_844, klev_var_845)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_856(kidia_var_843 : kfdia_var_844)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_857(kidia_var_843 : kfdia_var_844, klev_var_845)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_858(kidia_var_843 : kfdia_var_844, klev_var_845)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_859(kidia_var_843 : kfdia_var_844, klev_var_845)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_860(kidia_var_843 : kfdia_var_844, klev_var_845)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_861(kidia_var_843 : kfdia_var_844, klev_var_845)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_862(kidia_var_843 : kfdia_var_844, klev_var_845)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_863(kidia_var_843 : kfdia_var_844, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_864(kidia_var_843 : kfdia_var_844, klev_var_845, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_865(kidia_var_843 : kfdia_var_844, klev_var_845, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_866(kidia_var_843 : kfdia_var_844)
  INTEGER(KIND = 4) :: ig_var_867, ind0_var_868, ind1_var_869, inds_var_870, indf_var_871, i_lay_var_872, i_laysolfr_var_873(kidia_var_843 : kfdia_var_844), i_nlayers_var_874, iplon_var_875
  INTEGER(KIND = 4) :: laytrop_min_var_876, laytrop_max_var_877
  REAL(KIND = 8) :: z_tauray_var_878
  laytrop_min_var_876 = MINVAL(k_laytrop_var_856(kidia_var_843 : kfdia_var_844))
  laytrop_max_var_877 = MAXVAL(k_laytrop_var_856(kidia_var_843 : kfdia_var_844))
  i_nlayers_var_874 = klev_var_845
  DO iplon_var_875 = kidia_var_843, kfdia_var_844
    i_laysolfr_var_873(iplon_var_875) = k_laytrop_var_856(iplon_var_875)
  END DO
  DO i_lay_var_872 = 1, laytrop_min_var_876
    DO iplon_var_875 = kidia_var_843, kfdia_var_844
      IF (k_jp_var_850(iplon_var_875, i_lay_var_872) < layreffr_var_250 .AND. k_jp_var_850(iplon_var_875, i_lay_var_872 + 1) >= layreffr_var_250) i_laysolfr_var_873(iplon_var_875) = MIN(i_lay_var_872 + 1, k_laytrop_var_856(iplon_var_875))
      ind0_var_868 = ((k_jp_var_850(iplon_var_875, i_lay_var_872) - 1) * 5 + (k_jt_var_851(iplon_var_875, i_lay_var_872) - 1)) * nspa_var_313(20) + 1
      ind1_var_869 = (k_jp_var_850(iplon_var_875, i_lay_var_872) * 5 + (k_jt1_var_852(iplon_var_875, i_lay_var_872) - 1)) * nspa_var_313(20) + 1
      inds_var_870 = k_indself_var_859(iplon_var_875, i_lay_var_872)
      indf_var_871 = k_indfor_var_862(iplon_var_875, i_lay_var_872)
      z_tauray_var_878 = p_colmol_var_855(iplon_var_875, i_lay_var_872) * rayl_var_249
      DO ig_var_867 = 1, 10
        p_taug_var_864(iplon_var_875, i_lay_var_872, ig_var_867) = p_colh2o_var_853(iplon_var_875, i_lay_var_872) * ((p_fac00_var_846(iplon_var_875, i_lay_var_872) * absa_var_251(ind0_var_868, ig_var_867) + p_fac10_var_848(iplon_var_875, i_lay_var_872) * absa_var_251(ind0_var_868 + 1, ig_var_867) + p_fac01_var_847(iplon_var_875, i_lay_var_872) * absa_var_251(ind1_var_869, ig_var_867) + p_fac11_var_849(iplon_var_875, i_lay_var_872) * absa_var_251(ind1_var_869 + 1, ig_var_867)) + p_selffac_var_857(iplon_var_875, i_lay_var_872) * (selfrefc_var_253(inds_var_870, ig_var_867) + p_selffrac_var_858(iplon_var_875, i_lay_var_872) * (selfrefc_var_253(inds_var_870 + 1, ig_var_867) - selfrefc_var_253(inds_var_870, ig_var_867))) + p_forfac_var_860(iplon_var_875, i_lay_var_872) * (forrefc_var_254(indf_var_871, ig_var_867) + p_forfrac_var_861(iplon_var_875, i_lay_var_872) * (forrefc_var_254(indf_var_871 + 1, ig_var_867) - forrefc_var_254(indf_var_871, ig_var_867)))) + p_colch4_var_854(iplon_var_875, i_lay_var_872) * absch4c(ig_var_867)
        p_taur_var_865(iplon_var_875, i_lay_var_872, ig_var_867) = z_tauray_var_878
        IF (i_lay_var_872 == i_laysolfr_var_873(iplon_var_875)) p_sfluxzen_var_863(iplon_var_875, ig_var_867) = sfluxrefc_var_255(ig_var_867)
      END DO
    END DO
  END DO
  DO i_lay_var_872 = laytrop_min_var_876 + 1, laytrop_max_var_877
    DO iplon_var_875 = kidia_var_843, kfdia_var_844
      IF (i_lay_var_872 <= k_laytrop_var_856(iplon_var_875)) THEN
        IF (k_jp_var_850(iplon_var_875, i_lay_var_872) < layreffr_var_250 .AND. k_jp_var_850(iplon_var_875, i_lay_var_872 + 1) >= layreffr_var_250) i_laysolfr_var_873(iplon_var_875) = MIN(i_lay_var_872 + 1, k_laytrop_var_856(iplon_var_875))
        ind0_var_868 = ((k_jp_var_850(iplon_var_875, i_lay_var_872) - 1) * 5 + (k_jt_var_851(iplon_var_875, i_lay_var_872) - 1)) * nspa_var_313(20) + 1
        ind1_var_869 = (k_jp_var_850(iplon_var_875, i_lay_var_872) * 5 + (k_jt1_var_852(iplon_var_875, i_lay_var_872) - 1)) * nspa_var_313(20) + 1
        inds_var_870 = k_indself_var_859(iplon_var_875, i_lay_var_872)
        indf_var_871 = k_indfor_var_862(iplon_var_875, i_lay_var_872)
        z_tauray_var_878 = p_colmol_var_855(iplon_var_875, i_lay_var_872) * rayl_var_249
        DO ig_var_867 = 1, 10
          p_taug_var_864(iplon_var_875, i_lay_var_872, ig_var_867) = p_colh2o_var_853(iplon_var_875, i_lay_var_872) * ((p_fac00_var_846(iplon_var_875, i_lay_var_872) * absa_var_251(ind0_var_868, ig_var_867) + p_fac10_var_848(iplon_var_875, i_lay_var_872) * absa_var_251(ind0_var_868 + 1, ig_var_867) + p_fac01_var_847(iplon_var_875, i_lay_var_872) * absa_var_251(ind1_var_869, ig_var_867) + p_fac11_var_849(iplon_var_875, i_lay_var_872) * absa_var_251(ind1_var_869 + 1, ig_var_867)) + p_selffac_var_857(iplon_var_875, i_lay_var_872) * (selfrefc_var_253(inds_var_870, ig_var_867) + p_selffrac_var_858(iplon_var_875, i_lay_var_872) * (selfrefc_var_253(inds_var_870 + 1, ig_var_867) - selfrefc_var_253(inds_var_870, ig_var_867))) + p_forfac_var_860(iplon_var_875, i_lay_var_872) * (forrefc_var_254(indf_var_871, ig_var_867) + p_forfrac_var_861(iplon_var_875, i_lay_var_872) * (forrefc_var_254(indf_var_871 + 1, ig_var_867) - forrefc_var_254(indf_var_871, ig_var_867)))) + p_colch4_var_854(iplon_var_875, i_lay_var_872) * absch4c(ig_var_867)
          p_taur_var_865(iplon_var_875, i_lay_var_872, ig_var_867) = z_tauray_var_878
          IF (i_lay_var_872 == i_laysolfr_var_873(iplon_var_875)) p_sfluxzen_var_863(iplon_var_875, ig_var_867) = sfluxrefc_var_255(ig_var_867)
        END DO
      ELSE
        ind0_var_868 = ((k_jp_var_850(iplon_var_875, i_lay_var_872) - 13) * 5 + (k_jt_var_851(iplon_var_875, i_lay_var_872) - 1)) * nspb_var_314(20) + 1
        ind1_var_869 = ((k_jp_var_850(iplon_var_875, i_lay_var_872) - 12) * 5 + (k_jt1_var_852(iplon_var_875, i_lay_var_872) - 1)) * nspb_var_314(20) + 1
        indf_var_871 = k_indfor_var_862(iplon_var_875, i_lay_var_872)
        z_tauray_var_878 = p_colmol_var_855(iplon_var_875, i_lay_var_872) * rayl_var_249
        DO ig_var_867 = 1, 10
          p_taug_var_864(iplon_var_875, i_lay_var_872, ig_var_867) = p_colh2o_var_853(iplon_var_875, i_lay_var_872) * (p_fac00_var_846(iplon_var_875, i_lay_var_872) * absb_var_252(ind0_var_868, ig_var_867) + p_fac10_var_848(iplon_var_875, i_lay_var_872) * absb_var_252(ind0_var_868 + 1, ig_var_867) + p_fac01_var_847(iplon_var_875, i_lay_var_872) * absb_var_252(ind1_var_869, ig_var_867) + p_fac11_var_849(iplon_var_875, i_lay_var_872) * absb_var_252(ind1_var_869 + 1, ig_var_867) + p_forfac_var_860(iplon_var_875, i_lay_var_872) * (forrefc_var_254(indf_var_871, ig_var_867) + p_forfrac_var_861(iplon_var_875, i_lay_var_872) * (forrefc_var_254(indf_var_871 + 1, ig_var_867) - forrefc_var_254(indf_var_871, ig_var_867)))) + p_colch4_var_854(iplon_var_875, i_lay_var_872) * absch4c(ig_var_867)
          p_taur_var_865(iplon_var_875, i_lay_var_872, ig_var_867) = z_tauray_var_878
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_872 = laytrop_max_var_877 + 1, i_nlayers_var_874
    DO iplon_var_875 = kidia_var_843, kfdia_var_844
      ind0_var_868 = ((k_jp_var_850(iplon_var_875, i_lay_var_872) - 13) * 5 + (k_jt_var_851(iplon_var_875, i_lay_var_872) - 1)) * nspb_var_314(20) + 1
      ind1_var_869 = ((k_jp_var_850(iplon_var_875, i_lay_var_872) - 12) * 5 + (k_jt1_var_852(iplon_var_875, i_lay_var_872) - 1)) * nspb_var_314(20) + 1
      indf_var_871 = k_indfor_var_862(iplon_var_875, i_lay_var_872)
      z_tauray_var_878 = p_colmol_var_855(iplon_var_875, i_lay_var_872) * rayl_var_249
      DO ig_var_867 = 1, 10
        p_taug_var_864(iplon_var_875, i_lay_var_872, ig_var_867) = p_colh2o_var_853(iplon_var_875, i_lay_var_872) * (p_fac00_var_846(iplon_var_875, i_lay_var_872) * absb_var_252(ind0_var_868, ig_var_867) + p_fac10_var_848(iplon_var_875, i_lay_var_872) * absb_var_252(ind0_var_868 + 1, ig_var_867) + p_fac01_var_847(iplon_var_875, i_lay_var_872) * absb_var_252(ind1_var_869, ig_var_867) + p_fac11_var_849(iplon_var_875, i_lay_var_872) * absb_var_252(ind1_var_869 + 1, ig_var_867) + p_forfac_var_860(iplon_var_875, i_lay_var_872) * (forrefc_var_254(indf_var_871, ig_var_867) + p_forfrac_var_861(iplon_var_875, i_lay_var_872) * (forrefc_var_254(indf_var_871 + 1, ig_var_867) - forrefc_var_254(indf_var_871, ig_var_867)))) + p_colch4_var_854(iplon_var_875, i_lay_var_872) * absch4c(ig_var_867)
        p_taur_var_865(iplon_var_875, i_lay_var_872, ig_var_867) = z_tauray_var_878
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol20
SUBROUTINE srtm_taumol22(kidia_var_879, kfdia_var_880, klev_var_881, p_fac00_var_882, p_fac01_var_883, p_fac10_var_884, p_fac11_var_885, k_jp_var_886, k_jt_var_887, k_jt1_var_888, p_oneminus_var_889, p_colh2o_var_890, p_colmol_var_891, p_colo2_var_892, k_laytrop_var_893, p_selffac_var_894, p_selffrac_var_895, k_indself_var_896, p_forfac_var_897, p_forfrac_var_898, k_indfor_var_899, p_sfluxzen_var_900, p_taug_var_901, p_taur_var_902, prmu0_var_903)
  USE yoesrta22, ONLY: absa_var_267, absb_var_268, forrefc_var_270, layreffr_var_266, rayl_var_264, selfrefc_var_269, sfluxrefc_var_271, strrat_var_265
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_879, kfdia_var_880
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_881
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_882(kidia_var_879 : kfdia_var_880, klev_var_881)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_883(kidia_var_879 : kfdia_var_880, klev_var_881)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_884(kidia_var_879 : kfdia_var_880, klev_var_881)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_885(kidia_var_879 : kfdia_var_880, klev_var_881)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_886(kidia_var_879 : kfdia_var_880, klev_var_881)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_887(kidia_var_879 : kfdia_var_880, klev_var_881)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_888(kidia_var_879 : kfdia_var_880, klev_var_881)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_889(kidia_var_879 : kfdia_var_880)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_890(kidia_var_879 : kfdia_var_880, klev_var_881)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_891(kidia_var_879 : kfdia_var_880, klev_var_881)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_892(kidia_var_879 : kfdia_var_880, klev_var_881)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_893(kidia_var_879 : kfdia_var_880)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_894(kidia_var_879 : kfdia_var_880, klev_var_881)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_895(kidia_var_879 : kfdia_var_880, klev_var_881)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_896(kidia_var_879 : kfdia_var_880, klev_var_881)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_897(kidia_var_879 : kfdia_var_880, klev_var_881)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_898(kidia_var_879 : kfdia_var_880, klev_var_881)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_899(kidia_var_879 : kfdia_var_880, klev_var_881)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_900(kidia_var_879 : kfdia_var_880, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_901(kidia_var_879 : kfdia_var_880, klev_var_881, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_902(kidia_var_879 : kfdia_var_880, klev_var_881, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_903(kidia_var_879 : kfdia_var_880)
  INTEGER(KIND = 4) :: ig_var_904, ind0_var_905, ind1_var_906, inds_var_907, indf_var_908, js_var_909, i_lay_var_910, i_laysolfr_var_911(kidia_var_879 : kfdia_var_880), i_nlayers_var_912, iplon_var_913
  INTEGER(KIND = 4) :: laytrop_min_var_914, laytrop_max_var_915
  REAL(KIND = 8) :: z_fs_var_916, z_speccomb_var_917, z_specmult_var_918, z_specparm_var_919, z_tauray_var_920, z_o2adj, z_o2cont
  laytrop_min_var_914 = MINVAL(k_laytrop_var_893(kidia_var_879 : kfdia_var_880))
  laytrop_max_var_915 = MAXVAL(k_laytrop_var_893(kidia_var_879 : kfdia_var_880))
  i_nlayers_var_912 = klev_var_881
  z_o2adj = 1.6D0
  DO iplon_var_913 = kidia_var_879, kfdia_var_880
    i_laysolfr_var_911(iplon_var_913) = k_laytrop_var_893(iplon_var_913)
  END DO
  DO i_lay_var_910 = 1, laytrop_min_var_914
    DO iplon_var_913 = kidia_var_879, kfdia_var_880
      IF (k_jp_var_886(iplon_var_913, i_lay_var_910) < layreffr_var_266 .AND. k_jp_var_886(iplon_var_913, i_lay_var_910 + 1) >= layreffr_var_266) i_laysolfr_var_911(iplon_var_913) = MIN(i_lay_var_910 + 1, k_laytrop_var_893(iplon_var_913))
      z_o2cont = 0.000435D0 * p_colo2_var_892(iplon_var_913, i_lay_var_910) / (700.0D0)
      z_speccomb_var_917 = p_colh2o_var_890(iplon_var_913, i_lay_var_910) + z_o2adj * strrat_var_265 * p_colo2_var_892(iplon_var_913, i_lay_var_910)
      z_specparm_var_919 = p_colh2o_var_890(iplon_var_913, i_lay_var_910) / z_speccomb_var_917
      z_specparm_var_919 = MIN(p_oneminus_var_889(iplon_var_913), z_specparm_var_919)
      z_specmult_var_918 = 8.0D0 * (z_specparm_var_919)
      js_var_909 = 1 + INT(z_specmult_var_918)
      z_fs_var_916 = z_specmult_var_918 - AINT(z_specmult_var_918)
      ind0_var_905 = ((k_jp_var_886(iplon_var_913, i_lay_var_910) - 1) * 5 + (k_jt_var_887(iplon_var_913, i_lay_var_910) - 1)) * nspa_var_313(22) + js_var_909
      ind1_var_906 = (k_jp_var_886(iplon_var_913, i_lay_var_910) * 5 + (k_jt1_var_888(iplon_var_913, i_lay_var_910) - 1)) * nspa_var_313(22) + js_var_909
      inds_var_907 = k_indself_var_896(iplon_var_913, i_lay_var_910)
      indf_var_908 = k_indfor_var_899(iplon_var_913, i_lay_var_910)
      z_tauray_var_920 = p_colmol_var_891(iplon_var_913, i_lay_var_910) * rayl_var_264
      DO ig_var_904 = 1, 2
        p_taug_var_901(iplon_var_913, i_lay_var_910, ig_var_904) = z_speccomb_var_917 * ((1.0D0 - z_fs_var_916) * (absa_var_267(ind0_var_905, ig_var_904) * p_fac00_var_882(iplon_var_913, i_lay_var_910) + absa_var_267(ind0_var_905 + 9, ig_var_904) * p_fac10_var_884(iplon_var_913, i_lay_var_910) + absa_var_267(ind1_var_906, ig_var_904) * p_fac01_var_883(iplon_var_913, i_lay_var_910) + absa_var_267(ind1_var_906 + 9, ig_var_904) * p_fac11_var_885(iplon_var_913, i_lay_var_910)) + z_fs_var_916 * (absa_var_267(ind0_var_905 + 1, ig_var_904) * p_fac00_var_882(iplon_var_913, i_lay_var_910) + absa_var_267(ind0_var_905 + 10, ig_var_904) * p_fac10_var_884(iplon_var_913, i_lay_var_910) + absa_var_267(ind1_var_906 + 1, ig_var_904) * p_fac01_var_883(iplon_var_913, i_lay_var_910) + absa_var_267(ind1_var_906 + 10, ig_var_904) * p_fac11_var_885(iplon_var_913, i_lay_var_910))) + p_colh2o_var_890(iplon_var_913, i_lay_var_910) * (p_selffac_var_894(iplon_var_913, i_lay_var_910) * (selfrefc_var_269(inds_var_907, ig_var_904) + p_selffrac_var_895(iplon_var_913, i_lay_var_910) * (selfrefc_var_269(inds_var_907 + 1, ig_var_904) - selfrefc_var_269(inds_var_907, ig_var_904))) + p_forfac_var_897(iplon_var_913, i_lay_var_910) * (forrefc_var_270(indf_var_908, ig_var_904) + p_forfrac_var_898(iplon_var_913, i_lay_var_910) * (forrefc_var_270(indf_var_908 + 1, ig_var_904) - forrefc_var_270(indf_var_908, ig_var_904)))) + z_o2cont
        IF (i_lay_var_910 == i_laysolfr_var_911(iplon_var_913)) p_sfluxzen_var_900(iplon_var_913, ig_var_904) = sfluxrefc_var_271(ig_var_904, js_var_909) + z_fs_var_916 * (sfluxrefc_var_271(ig_var_904, js_var_909 + 1) - sfluxrefc_var_271(ig_var_904, js_var_909))
        p_taur_var_902(iplon_var_913, i_lay_var_910, ig_var_904) = z_tauray_var_920
      END DO
    END DO
  END DO
  DO i_lay_var_910 = laytrop_min_var_914 + 1, laytrop_max_var_915
    DO iplon_var_913 = kidia_var_879, kfdia_var_880
      IF (i_lay_var_910 <= k_laytrop_var_893(iplon_var_913)) THEN
        IF (k_jp_var_886(iplon_var_913, i_lay_var_910) < layreffr_var_266 .AND. k_jp_var_886(iplon_var_913, i_lay_var_910 + 1) >= layreffr_var_266) i_laysolfr_var_911(iplon_var_913) = MIN(i_lay_var_910 + 1, k_laytrop_var_893(iplon_var_913))
        z_o2cont = 0.000435D0 * p_colo2_var_892(iplon_var_913, i_lay_var_910) / (700.0D0)
        z_speccomb_var_917 = p_colh2o_var_890(iplon_var_913, i_lay_var_910) + z_o2adj * strrat_var_265 * p_colo2_var_892(iplon_var_913, i_lay_var_910)
        z_specparm_var_919 = p_colh2o_var_890(iplon_var_913, i_lay_var_910) / z_speccomb_var_917
        z_specparm_var_919 = MIN(p_oneminus_var_889(iplon_var_913), z_specparm_var_919)
        z_specmult_var_918 = 8.0D0 * (z_specparm_var_919)
        js_var_909 = 1 + INT(z_specmult_var_918)
        z_fs_var_916 = z_specmult_var_918 - AINT(z_specmult_var_918)
        ind0_var_905 = ((k_jp_var_886(iplon_var_913, i_lay_var_910) - 1) * 5 + (k_jt_var_887(iplon_var_913, i_lay_var_910) - 1)) * nspa_var_313(22) + js_var_909
        ind1_var_906 = (k_jp_var_886(iplon_var_913, i_lay_var_910) * 5 + (k_jt1_var_888(iplon_var_913, i_lay_var_910) - 1)) * nspa_var_313(22) + js_var_909
        inds_var_907 = k_indself_var_896(iplon_var_913, i_lay_var_910)
        indf_var_908 = k_indfor_var_899(iplon_var_913, i_lay_var_910)
        z_tauray_var_920 = p_colmol_var_891(iplon_var_913, i_lay_var_910) * rayl_var_264
        DO ig_var_904 = 1, 2
          p_taug_var_901(iplon_var_913, i_lay_var_910, ig_var_904) = z_speccomb_var_917 * ((1.0D0 - z_fs_var_916) * (absa_var_267(ind0_var_905, ig_var_904) * p_fac00_var_882(iplon_var_913, i_lay_var_910) + absa_var_267(ind0_var_905 + 9, ig_var_904) * p_fac10_var_884(iplon_var_913, i_lay_var_910) + absa_var_267(ind1_var_906, ig_var_904) * p_fac01_var_883(iplon_var_913, i_lay_var_910) + absa_var_267(ind1_var_906 + 9, ig_var_904) * p_fac11_var_885(iplon_var_913, i_lay_var_910)) + z_fs_var_916 * (absa_var_267(ind0_var_905 + 1, ig_var_904) * p_fac00_var_882(iplon_var_913, i_lay_var_910) + absa_var_267(ind0_var_905 + 10, ig_var_904) * p_fac10_var_884(iplon_var_913, i_lay_var_910) + absa_var_267(ind1_var_906 + 1, ig_var_904) * p_fac01_var_883(iplon_var_913, i_lay_var_910) + absa_var_267(ind1_var_906 + 10, ig_var_904) * p_fac11_var_885(iplon_var_913, i_lay_var_910))) + p_colh2o_var_890(iplon_var_913, i_lay_var_910) * (p_selffac_var_894(iplon_var_913, i_lay_var_910) * (selfrefc_var_269(inds_var_907, ig_var_904) + p_selffrac_var_895(iplon_var_913, i_lay_var_910) * (selfrefc_var_269(inds_var_907 + 1, ig_var_904) - selfrefc_var_269(inds_var_907, ig_var_904))) + p_forfac_var_897(iplon_var_913, i_lay_var_910) * (forrefc_var_270(indf_var_908, ig_var_904) + p_forfrac_var_898(iplon_var_913, i_lay_var_910) * (forrefc_var_270(indf_var_908 + 1, ig_var_904) - forrefc_var_270(indf_var_908, ig_var_904)))) + z_o2cont
          IF (i_lay_var_910 == i_laysolfr_var_911(iplon_var_913)) p_sfluxzen_var_900(iplon_var_913, ig_var_904) = sfluxrefc_var_271(ig_var_904, js_var_909) + z_fs_var_916 * (sfluxrefc_var_271(ig_var_904, js_var_909 + 1) - sfluxrefc_var_271(ig_var_904, js_var_909))
          p_taur_var_902(iplon_var_913, i_lay_var_910, ig_var_904) = z_tauray_var_920
        END DO
      ELSE
        z_o2cont = 0.000435D0 * p_colo2_var_892(iplon_var_913, i_lay_var_910) / (700.0D0)
        ind0_var_905 = ((k_jp_var_886(iplon_var_913, i_lay_var_910) - 13) * 5 + (k_jt_var_887(iplon_var_913, i_lay_var_910) - 1)) * nspb_var_314(22) + 1
        ind1_var_906 = ((k_jp_var_886(iplon_var_913, i_lay_var_910) - 12) * 5 + (k_jt1_var_888(iplon_var_913, i_lay_var_910) - 1)) * nspb_var_314(22) + 1
        z_tauray_var_920 = p_colmol_var_891(iplon_var_913, i_lay_var_910) * rayl_var_264
        DO ig_var_904 = 1, 2
          p_taug_var_901(iplon_var_913, i_lay_var_910, ig_var_904) = p_colo2_var_892(iplon_var_913, i_lay_var_910) * z_o2adj * (p_fac00_var_882(iplon_var_913, i_lay_var_910) * absb_var_268(ind0_var_905, ig_var_904) + p_fac10_var_884(iplon_var_913, i_lay_var_910) * absb_var_268(ind0_var_905 + 1, ig_var_904) + p_fac01_var_883(iplon_var_913, i_lay_var_910) * absb_var_268(ind1_var_906, ig_var_904) + p_fac11_var_885(iplon_var_913, i_lay_var_910) * absb_var_268(ind1_var_906 + 1, ig_var_904)) + z_o2cont
          p_taur_var_902(iplon_var_913, i_lay_var_910, ig_var_904) = z_tauray_var_920
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_910 = laytrop_max_var_915 + 1, i_nlayers_var_912
    DO iplon_var_913 = kidia_var_879, kfdia_var_880
      z_o2cont = 0.000435D0 * p_colo2_var_892(iplon_var_913, i_lay_var_910) / (700.0D0)
      ind0_var_905 = ((k_jp_var_886(iplon_var_913, i_lay_var_910) - 13) * 5 + (k_jt_var_887(iplon_var_913, i_lay_var_910) - 1)) * nspb_var_314(22) + 1
      ind1_var_906 = ((k_jp_var_886(iplon_var_913, i_lay_var_910) - 12) * 5 + (k_jt1_var_888(iplon_var_913, i_lay_var_910) - 1)) * nspb_var_314(22) + 1
      z_tauray_var_920 = p_colmol_var_891(iplon_var_913, i_lay_var_910) * rayl_var_264
      DO ig_var_904 = 1, 2
        p_taug_var_901(iplon_var_913, i_lay_var_910, ig_var_904) = p_colo2_var_892(iplon_var_913, i_lay_var_910) * z_o2adj * (p_fac00_var_882(iplon_var_913, i_lay_var_910) * absb_var_268(ind0_var_905, ig_var_904) + p_fac10_var_884(iplon_var_913, i_lay_var_910) * absb_var_268(ind0_var_905 + 1, ig_var_904) + p_fac01_var_883(iplon_var_913, i_lay_var_910) * absb_var_268(ind1_var_906, ig_var_904) + p_fac11_var_885(iplon_var_913, i_lay_var_910) * absb_var_268(ind1_var_906 + 1, ig_var_904)) + z_o2cont
        p_taur_var_902(iplon_var_913, i_lay_var_910, ig_var_904) = z_tauray_var_920
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol22
SUBROUTINE srtm_taumol23(kidia_var_921, kfdia_var_922, klev_var_923, p_fac00_var_924, p_fac01_var_925, p_fac10_var_926, p_fac11_var_927, k_jp_var_928, k_jt_var_929, k_jt1_var_930, p_colh2o_var_931, p_colmol_var_932, k_laytrop_var_933, p_selffac_var_934, p_selffrac_var_935, k_indself_var_936, p_forfac_var_937, p_forfrac_var_938, k_indfor_var_939, p_sfluxzen_var_940, p_taug_var_941, p_taur_var_942, prmu0_var_943)
  USE yoesrta23, ONLY: absa_var_273, forrefc_var_275, givfac, layreffr_var_272, raylc_var_277, selfrefc_var_274, sfluxrefc_var_276
  USE yoesrtwn, ONLY: nspa_var_313
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_921, kfdia_var_922
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_923
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_924(kidia_var_921 : kfdia_var_922, klev_var_923)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_925(kidia_var_921 : kfdia_var_922, klev_var_923)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_926(kidia_var_921 : kfdia_var_922, klev_var_923)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_927(kidia_var_921 : kfdia_var_922, klev_var_923)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_928(kidia_var_921 : kfdia_var_922, klev_var_923)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_929(kidia_var_921 : kfdia_var_922, klev_var_923)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_930(kidia_var_921 : kfdia_var_922, klev_var_923)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_931(kidia_var_921 : kfdia_var_922, klev_var_923)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_932(kidia_var_921 : kfdia_var_922, klev_var_923)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_933(kidia_var_921 : kfdia_var_922)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_934(kidia_var_921 : kfdia_var_922, klev_var_923)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_935(kidia_var_921 : kfdia_var_922, klev_var_923)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_936(kidia_var_921 : kfdia_var_922, klev_var_923)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_937(kidia_var_921 : kfdia_var_922, klev_var_923)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_938(kidia_var_921 : kfdia_var_922, klev_var_923)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_939(kidia_var_921 : kfdia_var_922, klev_var_923)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_940(kidia_var_921 : kfdia_var_922, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_941(kidia_var_921 : kfdia_var_922, klev_var_923, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_942(kidia_var_921 : kfdia_var_922, klev_var_923, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_943(kidia_var_921 : kfdia_var_922)
  INTEGER(KIND = 4) :: ig_var_944, ind0_var_945, ind1_var_946, inds_var_947, indf_var_948, i_lay_var_949, i_laysolfr_var_950(kidia_var_921 : kfdia_var_922), i_nlayers_var_951, iplon_var_952
  INTEGER(KIND = 4) :: laytrop_min_var_953, laytrop_max_var_954
  REAL(KIND = 8) :: z_tauray_var_955
  laytrop_min_var_953 = MINVAL(k_laytrop_var_933(kidia_var_921 : kfdia_var_922))
  laytrop_max_var_954 = MAXVAL(k_laytrop_var_933(kidia_var_921 : kfdia_var_922))
  i_nlayers_var_951 = klev_var_923
  DO iplon_var_952 = kidia_var_921, kfdia_var_922
    i_laysolfr_var_950(iplon_var_952) = k_laytrop_var_933(iplon_var_952)
  END DO
  DO i_lay_var_949 = 1, laytrop_min_var_953
    DO iplon_var_952 = kidia_var_921, kfdia_var_922
      IF (k_jp_var_928(iplon_var_952, i_lay_var_949) < layreffr_var_272 .AND. k_jp_var_928(iplon_var_952, i_lay_var_949 + 1) >= layreffr_var_272) i_laysolfr_var_950(iplon_var_952) = MIN(i_lay_var_949 + 1, k_laytrop_var_933(iplon_var_952))
      ind0_var_945 = ((k_jp_var_928(iplon_var_952, i_lay_var_949) - 1) * 5 + (k_jt_var_929(iplon_var_952, i_lay_var_949) - 1)) * nspa_var_313(23) + 1
      ind1_var_946 = (k_jp_var_928(iplon_var_952, i_lay_var_949) * 5 + (k_jt1_var_930(iplon_var_952, i_lay_var_949) - 1)) * nspa_var_313(23) + 1
      inds_var_947 = k_indself_var_936(iplon_var_952, i_lay_var_949)
      indf_var_948 = k_indfor_var_939(iplon_var_952, i_lay_var_949)
      DO ig_var_944 = 1, 10
        z_tauray_var_955 = p_colmol_var_932(iplon_var_952, i_lay_var_949) * raylc_var_277(ig_var_944)
        p_taug_var_941(iplon_var_952, i_lay_var_949, ig_var_944) = p_colh2o_var_931(iplon_var_952, i_lay_var_949) * (givfac * (p_fac00_var_924(iplon_var_952, i_lay_var_949) * absa_var_273(ind0_var_945, ig_var_944) + p_fac10_var_926(iplon_var_952, i_lay_var_949) * absa_var_273(ind0_var_945 + 1, ig_var_944) + p_fac01_var_925(iplon_var_952, i_lay_var_949) * absa_var_273(ind1_var_946, ig_var_944) + p_fac11_var_927(iplon_var_952, i_lay_var_949) * absa_var_273(ind1_var_946 + 1, ig_var_944)) + p_selffac_var_934(iplon_var_952, i_lay_var_949) * (selfrefc_var_274(inds_var_947, ig_var_944) + p_selffrac_var_935(iplon_var_952, i_lay_var_949) * (selfrefc_var_274(inds_var_947 + 1, ig_var_944) - selfrefc_var_274(inds_var_947, ig_var_944))) + p_forfac_var_937(iplon_var_952, i_lay_var_949) * (forrefc_var_275(indf_var_948, ig_var_944) + p_forfrac_var_938(iplon_var_952, i_lay_var_949) * (forrefc_var_275(indf_var_948 + 1, ig_var_944) - forrefc_var_275(indf_var_948, ig_var_944))))
        IF (i_lay_var_949 == i_laysolfr_var_950(iplon_var_952)) p_sfluxzen_var_940(iplon_var_952, ig_var_944) = sfluxrefc_var_276(ig_var_944)
        p_taur_var_942(iplon_var_952, i_lay_var_949, ig_var_944) = z_tauray_var_955
      END DO
    END DO
  END DO
  DO i_lay_var_949 = laytrop_min_var_953 + 1, laytrop_max_var_954
    DO iplon_var_952 = kidia_var_921, kfdia_var_922
      IF (i_lay_var_949 <= k_laytrop_var_933(iplon_var_952)) THEN
        IF (k_jp_var_928(iplon_var_952, i_lay_var_949) < layreffr_var_272 .AND. k_jp_var_928(iplon_var_952, i_lay_var_949 + 1) >= layreffr_var_272) i_laysolfr_var_950(iplon_var_952) = MIN(i_lay_var_949 + 1, k_laytrop_var_933(iplon_var_952))
        ind0_var_945 = ((k_jp_var_928(iplon_var_952, i_lay_var_949) - 1) * 5 + (k_jt_var_929(iplon_var_952, i_lay_var_949) - 1)) * nspa_var_313(23) + 1
        ind1_var_946 = (k_jp_var_928(iplon_var_952, i_lay_var_949) * 5 + (k_jt1_var_930(iplon_var_952, i_lay_var_949) - 1)) * nspa_var_313(23) + 1
        inds_var_947 = k_indself_var_936(iplon_var_952, i_lay_var_949)
        indf_var_948 = k_indfor_var_939(iplon_var_952, i_lay_var_949)
        DO ig_var_944 = 1, 10
          z_tauray_var_955 = p_colmol_var_932(iplon_var_952, i_lay_var_949) * raylc_var_277(ig_var_944)
          p_taug_var_941(iplon_var_952, i_lay_var_949, ig_var_944) = p_colh2o_var_931(iplon_var_952, i_lay_var_949) * (givfac * (p_fac00_var_924(iplon_var_952, i_lay_var_949) * absa_var_273(ind0_var_945, ig_var_944) + p_fac10_var_926(iplon_var_952, i_lay_var_949) * absa_var_273(ind0_var_945 + 1, ig_var_944) + p_fac01_var_925(iplon_var_952, i_lay_var_949) * absa_var_273(ind1_var_946, ig_var_944) + p_fac11_var_927(iplon_var_952, i_lay_var_949) * absa_var_273(ind1_var_946 + 1, ig_var_944)) + p_selffac_var_934(iplon_var_952, i_lay_var_949) * (selfrefc_var_274(inds_var_947, ig_var_944) + p_selffrac_var_935(iplon_var_952, i_lay_var_949) * (selfrefc_var_274(inds_var_947 + 1, ig_var_944) - selfrefc_var_274(inds_var_947, ig_var_944))) + p_forfac_var_937(iplon_var_952, i_lay_var_949) * (forrefc_var_275(indf_var_948, ig_var_944) + p_forfrac_var_938(iplon_var_952, i_lay_var_949) * (forrefc_var_275(indf_var_948 + 1, ig_var_944) - forrefc_var_275(indf_var_948, ig_var_944))))
          IF (i_lay_var_949 == i_laysolfr_var_950(iplon_var_952)) p_sfluxzen_var_940(iplon_var_952, ig_var_944) = sfluxrefc_var_276(ig_var_944)
          p_taur_var_942(iplon_var_952, i_lay_var_949, ig_var_944) = z_tauray_var_955
        END DO
      ELSE
        DO ig_var_944 = 1, 10
          p_taug_var_941(iplon_var_952, i_lay_var_949, ig_var_944) = 0.0D0
          p_taur_var_942(iplon_var_952, i_lay_var_949, ig_var_944) = p_colmol_var_932(iplon_var_952, i_lay_var_949) * raylc_var_277(ig_var_944)
        END DO
      END IF
    END DO
  END DO
  DO ig_var_944 = 1, 10
    DO i_lay_var_949 = laytrop_max_var_954 + 1, i_nlayers_var_951
      DO iplon_var_952 = kidia_var_921, kfdia_var_922
        p_taug_var_941(iplon_var_952, i_lay_var_949, ig_var_944) = 0.0D0
        p_taur_var_942(iplon_var_952, i_lay_var_949, ig_var_944) = p_colmol_var_932(iplon_var_952, i_lay_var_949) * raylc_var_277(ig_var_944)
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol23
SUBROUTINE rrtm_taumol5(kidia_var_956, kfdia_var_957, klev_var_958, taug_var_959, wx_var_960, p_tauaerl_var_961, fac00_var_962, fac01_var_963, fac10_var_964, fac11_var_965, forfac_var_984, forfrac_var_983, indfor_var_982, jp_var_966, jt_var_967, jt1_var_968, oneminus_var_969, colh2o_var_970, colco2_var_971, colo3_var_972, laytrop_var_973, selffac_var_974, selffrac_var_975, indself_var_976, fracs_var_977, rat_h2oco2_var_978, rat_h2oco2_1_var_979, rat_o3co2_var_980, rat_o3co2_1_var_981, minorfrac_var_985, indminor_var_986)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng5
  USE yoerrta5, ONLY: absa_var_175, absb_var_176, ccl4, forref_var_179, fracrefa_var_173, fracrefb_var_174, ka_mo3_var_177, selfref_var_178
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_956
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_957
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_958
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_959(kidia_var_956 : kfdia_var_957, 140, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: wx_var_960(kidia_var_956 : kfdia_var_957, 4, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_961(kidia_var_956 : kfdia_var_957, klev_var_958, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_962(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_963(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_964(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_965(kidia_var_956 : kfdia_var_957, klev_var_958)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_966(kidia_var_956 : kfdia_var_957, klev_var_958)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_967(kidia_var_956 : kfdia_var_957, klev_var_958)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_968(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_969
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_970(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_971(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_972(kidia_var_956 : kfdia_var_957, klev_var_958)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_973(kidia_var_956 : kfdia_var_957)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_974(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_975(kidia_var_956 : kfdia_var_957, klev_var_958)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_976(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_977(kidia_var_956 : kfdia_var_957, 140, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_978(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_979(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_var_980(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_1_var_981(kidia_var_956 : kfdia_var_957, klev_var_958)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_982(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_983(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_984(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_985(kidia_var_956 : kfdia_var_957, klev_var_958)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_986(kidia_var_956 : kfdia_var_957, klev_var_958)
  REAL(KIND = 8) :: speccomb_var_987, speccomb1_var_988, speccomb_mo3, speccomb_planck_var_989
  INTEGER(KIND = 4) :: ind0_var_990, ind1_var_991, inds_var_992, indf_var_993, indm_var_994
  INTEGER(KIND = 4) :: ig_var_995, js_var_996, lay_var_997, js1_var_998, jpl_var_999, jmo3
  REAL(KIND = 8) :: refrat_planck_a_var_1000, refrat_planck_b_var_1001, refrat_m_a_var_1002
  REAL(KIND = 8) :: fac000_var_1003, fac100_var_1004, fac200_var_1005, fac010_var_1006, fac110_var_1007, fac210_var_1008, fac001_var_1009, fac101_var_1010, fac201_var_1011, fac011_var_1012, fac111_var_1013, fac211_var_1014
  REAL(KIND = 8) :: p_var_1015, p4_var_1016, fk0_var_1017, fk1_var_1018, fk2_var_1019
  REAL(KIND = 8) :: taufor_var_1020, tauself_var_1021, tau_major_var_1022(16), tau_major1_var_1023(16), o3m1, o3m2, abso3_var_1024
  REAL(KIND = 8) :: fs_var_1025, specmult_var_1026, specparm_var_1027, fs1_var_1028, specmult1_var_1029, specparm1_var_1030, fpl_var_1031, specmult_planck_var_1032, specparm_planck_var_1033, fmo3, specmult_mo3, specparm_mo3
  INTEGER(KIND = 4) :: laytrop_min_var_1034, laytrop_max_var_1035
  INTEGER(KIND = 4) :: ixc_var_1036(klev_var_958), ixlow_var_1037(kfdia_var_957, klev_var_958), ixhigh_var_1038(kfdia_var_957, klev_var_958)
  INTEGER(KIND = 4) :: ich_var_1039, icl_var_1040, ixc0_var_1041, ixp_var_1042, jc_var_1043, jl_var_1044
  laytrop_min_var_1034 = MINVAL(laytrop_var_973)
  laytrop_max_var_1035 = MAXVAL(laytrop_var_973)
  ixlow_var_1037 = 0
  ixhigh_var_1038 = 0
  ixc_var_1036 = 0
  DO lay_var_997 = laytrop_min_var_1034 + 1, laytrop_max_var_1035
    icl_var_1040 = 0
    ich_var_1039 = 0
    DO jc_var_1043 = kidia_var_956, kfdia_var_957
      IF (lay_var_997 <= laytrop_var_973(jc_var_1043)) THEN
        icl_var_1040 = icl_var_1040 + 1
        ixlow_var_1037(icl_var_1040, lay_var_997) = jc_var_1043
      ELSE
        ich_var_1039 = ich_var_1039 + 1
        ixhigh_var_1038(ich_var_1039, lay_var_997) = jc_var_1043
      END IF
    END DO
    ixc_var_1036(lay_var_997) = icl_var_1040
  END DO
  refrat_planck_a_var_1000 = chi_mls(1, 5) / chi_mls(2, 5)
  refrat_planck_b_var_1001 = chi_mls(3, 43) / chi_mls(2, 43)
  refrat_m_a_var_1002 = chi_mls(1, 7) / chi_mls(2, 7)
  DO lay_var_997 = 1, laytrop_min_var_1034
    DO jl_var_1044 = kidia_var_956, kfdia_var_957
      speccomb_var_987 = colh2o_var_970(jl_var_1044, lay_var_997) + rat_h2oco2_var_978(jl_var_1044, lay_var_997) * colco2_var_971(jl_var_1044, lay_var_997)
      specparm_var_1027 = MIN(colh2o_var_970(jl_var_1044, lay_var_997) / speccomb_var_987, 0.999999D0)
      specmult_var_1026 = 8.0D0 * (specparm_var_1027)
      js_var_996 = 1 + INT(specmult_var_1026)
      fs_var_1025 = ((specmult_var_1026) - AINT((specmult_var_1026)))
      speccomb1_var_988 = colh2o_var_970(jl_var_1044, lay_var_997) + rat_h2oco2_1_var_979(jl_var_1044, lay_var_997) * colco2_var_971(jl_var_1044, lay_var_997)
      specparm1_var_1030 = MIN(colh2o_var_970(jl_var_1044, lay_var_997) / speccomb1_var_988, 0.999999D0)
      specmult1_var_1029 = 8.0D0 * (specparm1_var_1030)
      js1_var_998 = 1 + INT(specmult1_var_1029)
      fs1_var_1028 = ((specmult1_var_1029) - AINT((specmult1_var_1029)))
      speccomb_mo3 = colh2o_var_970(jl_var_1044, lay_var_997) + refrat_m_a_var_1002 * colco2_var_971(jl_var_1044, lay_var_997)
      specparm_mo3 = MIN(colh2o_var_970(jl_var_1044, lay_var_997) / speccomb_mo3, 0.999999D0)
      specmult_mo3 = 8.0D0 * specparm_mo3
      jmo3 = 1 + INT(specmult_mo3)
      fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))
      speccomb_planck_var_989 = colh2o_var_970(jl_var_1044, lay_var_997) + refrat_planck_a_var_1000 * colco2_var_971(jl_var_1044, lay_var_997)
      specparm_planck_var_1033 = MIN(colh2o_var_970(jl_var_1044, lay_var_997) / speccomb_planck_var_989, 0.999999D0)
      specmult_planck_var_1032 = 8.0D0 * specparm_planck_var_1033
      jpl_var_999 = 1 + INT(specmult_planck_var_1032)
      fpl_var_1031 = ((specmult_planck_var_1032) - AINT((specmult_planck_var_1032)))
      ind0_var_990 = ((jp_var_966(jl_var_1044, lay_var_997) - 1) * 5 + (jt_var_967(jl_var_1044, lay_var_997) - 1)) * nspa_var_216(5) + js_var_996
      ind1_var_991 = (jp_var_966(jl_var_1044, lay_var_997) * 5 + (jt1_var_968(jl_var_1044, lay_var_997) - 1)) * nspa_var_216(5) + js1_var_998
      inds_var_992 = indself_var_976(jl_var_1044, lay_var_997)
      indf_var_993 = indfor_var_982(jl_var_1044, lay_var_997)
      indm_var_994 = indminor_var_986(jl_var_1044, lay_var_997)
      IF (specparm_var_1027 .LT. 0.125D0) THEN
        p_var_1015 = fs_var_1025 - 1.0D0
        p4_var_1016 = p_var_1015 ** 4
        fk0_var_1017 = p4_var_1016
        fk1_var_1018 = 1.0D0 - p_var_1015 - 2.0D0 * p4_var_1016
        fk2_var_1019 = p_var_1015 + p4_var_1016
        fac000_var_1003 = fk0_var_1017 * fac00_var_962(jl_var_1044, lay_var_997)
        fac100_var_1004 = fk1_var_1018 * fac00_var_962(jl_var_1044, lay_var_997)
        fac200_var_1005 = fk2_var_1019 * fac00_var_962(jl_var_1044, lay_var_997)
        fac010_var_1006 = fk0_var_1017 * fac10_var_964(jl_var_1044, lay_var_997)
        fac110_var_1007 = fk1_var_1018 * fac10_var_964(jl_var_1044, lay_var_997)
        fac210_var_1008 = fk2_var_1019 * fac10_var_964(jl_var_1044, lay_var_997)
      ELSE IF (specparm_var_1027 .GT. 0.875D0) THEN
        p_var_1015 = - fs_var_1025
        p4_var_1016 = p_var_1015 ** 4
        fk0_var_1017 = p4_var_1016
        fk1_var_1018 = 1.0D0 - p_var_1015 - 2.0D0 * p4_var_1016
        fk2_var_1019 = p_var_1015 + p4_var_1016
        fac000_var_1003 = fk0_var_1017 * fac00_var_962(jl_var_1044, lay_var_997)
        fac100_var_1004 = fk1_var_1018 * fac00_var_962(jl_var_1044, lay_var_997)
        fac200_var_1005 = fk2_var_1019 * fac00_var_962(jl_var_1044, lay_var_997)
        fac010_var_1006 = fk0_var_1017 * fac10_var_964(jl_var_1044, lay_var_997)
        fac110_var_1007 = fk1_var_1018 * fac10_var_964(jl_var_1044, lay_var_997)
        fac210_var_1008 = fk2_var_1019 * fac10_var_964(jl_var_1044, lay_var_997)
      ELSE
        fac000_var_1003 = (1.0D0 - fs_var_1025) * fac00_var_962(jl_var_1044, lay_var_997)
        fac010_var_1006 = (1.0D0 - fs_var_1025) * fac10_var_964(jl_var_1044, lay_var_997)
        fac100_var_1004 = fs_var_1025 * fac00_var_962(jl_var_1044, lay_var_997)
        fac110_var_1007 = fs_var_1025 * fac10_var_964(jl_var_1044, lay_var_997)
        fac200_var_1005 = 0.0D0
        fac210_var_1008 = 0.0D0
      END IF
      IF (specparm1_var_1030 .LT. 0.125D0) THEN
        p_var_1015 = fs1_var_1028 - 1.0D0
        p4_var_1016 = p_var_1015 ** 4
        fk0_var_1017 = p4_var_1016
        fk1_var_1018 = 1.0D0 - p_var_1015 - 2.0D0 * p4_var_1016
        fk2_var_1019 = p_var_1015 + p4_var_1016
        fac001_var_1009 = fk0_var_1017 * fac01_var_963(jl_var_1044, lay_var_997)
        fac101_var_1010 = fk1_var_1018 * fac01_var_963(jl_var_1044, lay_var_997)
        fac201_var_1011 = fk2_var_1019 * fac01_var_963(jl_var_1044, lay_var_997)
        fac011_var_1012 = fk0_var_1017 * fac11_var_965(jl_var_1044, lay_var_997)
        fac111_var_1013 = fk1_var_1018 * fac11_var_965(jl_var_1044, lay_var_997)
        fac211_var_1014 = fk2_var_1019 * fac11_var_965(jl_var_1044, lay_var_997)
      ELSE IF (specparm1_var_1030 .GT. 0.875D0) THEN
        p_var_1015 = - fs1_var_1028
        p4_var_1016 = p_var_1015 ** 4
        fk0_var_1017 = p4_var_1016
        fk1_var_1018 = 1.0D0 - p_var_1015 - 2.0D0 * p4_var_1016
        fk2_var_1019 = p_var_1015 + p4_var_1016
        fac001_var_1009 = fk0_var_1017 * fac01_var_963(jl_var_1044, lay_var_997)
        fac101_var_1010 = fk1_var_1018 * fac01_var_963(jl_var_1044, lay_var_997)
        fac201_var_1011 = fk2_var_1019 * fac01_var_963(jl_var_1044, lay_var_997)
        fac011_var_1012 = fk0_var_1017 * fac11_var_965(jl_var_1044, lay_var_997)
        fac111_var_1013 = fk1_var_1018 * fac11_var_965(jl_var_1044, lay_var_997)
        fac211_var_1014 = fk2_var_1019 * fac11_var_965(jl_var_1044, lay_var_997)
      ELSE
        fac001_var_1009 = (1.0D0 - fs1_var_1028) * fac01_var_963(jl_var_1044, lay_var_997)
        fac011_var_1012 = (1.0D0 - fs1_var_1028) * fac11_var_965(jl_var_1044, lay_var_997)
        fac101_var_1010 = fs1_var_1028 * fac01_var_963(jl_var_1044, lay_var_997)
        fac111_var_1013 = fs1_var_1028 * fac11_var_965(jl_var_1044, lay_var_997)
        fac201_var_1011 = 0.0D0
        fac211_var_1014 = 0.0D0
      END IF
      IF (specparm_var_1027 .LT. 0.125D0) THEN
        tau_major_var_1022(1 : ng5) = speccomb_var_987 * (fac000_var_1003 * absa_var_175(ind0_var_990, 1 : 16) + fac100_var_1004 * absa_var_175(ind0_var_990 + 1, 1 : 16) + fac200_var_1005 * absa_var_175(ind0_var_990 + 2, 1 : 16) + fac010_var_1006 * absa_var_175(ind0_var_990 + 9, 1 : 16) + fac110_var_1007 * absa_var_175(ind0_var_990 + 10, 1 : 16) + fac210_var_1008 * absa_var_175(ind0_var_990 + 11, 1 : 16))
      ELSE IF (specparm_var_1027 .GT. 0.875D0) THEN
        tau_major_var_1022(1 : ng5) = speccomb_var_987 * (fac200_var_1005 * absa_var_175(ind0_var_990 - 1, 1 : 16) + fac100_var_1004 * absa_var_175(ind0_var_990, 1 : 16) + fac000_var_1003 * absa_var_175(ind0_var_990 + 1, 1 : 16) + fac210_var_1008 * absa_var_175(ind0_var_990 + 8, 1 : 16) + fac110_var_1007 * absa_var_175(ind0_var_990 + 9, 1 : 16) + fac010_var_1006 * absa_var_175(ind0_var_990 + 10, 1 : 16))
      ELSE
        tau_major_var_1022(1 : ng5) = speccomb_var_987 * (fac000_var_1003 * absa_var_175(ind0_var_990, 1 : 16) + fac100_var_1004 * absa_var_175(ind0_var_990 + 1, 1 : 16) + fac010_var_1006 * absa_var_175(ind0_var_990 + 9, 1 : 16) + fac110_var_1007 * absa_var_175(ind0_var_990 + 10, 1 : 16))
      END IF
      IF (specparm1_var_1030 .LT. 0.125D0) THEN
        tau_major1_var_1023(1 : ng5) = speccomb1_var_988 * (fac001_var_1009 * absa_var_175(ind1_var_991, 1 : 16) + fac101_var_1010 * absa_var_175(ind1_var_991 + 1, 1 : 16) + fac201_var_1011 * absa_var_175(ind1_var_991 + 2, 1 : 16) + fac011_var_1012 * absa_var_175(ind1_var_991 + 9, 1 : 16) + fac111_var_1013 * absa_var_175(ind1_var_991 + 10, 1 : 16) + fac211_var_1014 * absa_var_175(ind1_var_991 + 11, 1 : 16))
      ELSE IF (specparm1_var_1030 .GT. 0.875D0) THEN
        tau_major1_var_1023(1 : ng5) = speccomb1_var_988 * (fac201_var_1011 * absa_var_175(ind1_var_991 - 1, 1 : 16) + fac101_var_1010 * absa_var_175(ind1_var_991, 1 : 16) + fac001_var_1009 * absa_var_175(ind1_var_991 + 1, 1 : 16) + fac211_var_1014 * absa_var_175(ind1_var_991 + 8, 1 : 16) + fac111_var_1013 * absa_var_175(ind1_var_991 + 9, 1 : 16) + fac011_var_1012 * absa_var_175(ind1_var_991 + 10, 1 : 16))
      ELSE
        tau_major1_var_1023(1 : ng5) = speccomb1_var_988 * (fac001_var_1009 * absa_var_175(ind1_var_991, 1 : 16) + fac101_var_1010 * absa_var_175(ind1_var_991 + 1, 1 : 16) + fac011_var_1012 * absa_var_175(ind1_var_991 + 9, 1 : 16) + fac111_var_1013 * absa_var_175(ind1_var_991 + 10, 1 : 16))
      END IF
      DO ig_var_995 = 1, 16
        tauself_var_1021 = selffac_var_974(jl_var_1044, lay_var_997) * (selfref_var_178(inds_var_992, ig_var_995) + selffrac_var_975(jl_var_1044, lay_var_997) * (selfref_var_178(inds_var_992 + 1, ig_var_995) - selfref_var_178(inds_var_992, ig_var_995)))
        taufor_var_1020 = forfac_var_984(jl_var_1044, lay_var_997) * (forref_var_179(indf_var_993, ig_var_995) + forfrac_var_983(jl_var_1044, lay_var_997) * (forref_var_179(indf_var_993 + 1, ig_var_995) - forref_var_179(indf_var_993, ig_var_995)))
        o3m1 = ka_mo3_var_177(jmo3, indm_var_994, ig_var_995) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_994, ig_var_995) - ka_mo3_var_177(jmo3, indm_var_994, ig_var_995))
        o3m2 = ka_mo3_var_177(jmo3, indm_var_994 + 1, ig_var_995) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_994 + 1, ig_var_995) - ka_mo3_var_177(jmo3, indm_var_994 + 1, ig_var_995))
        abso3_var_1024 = o3m1 + minorfrac_var_985(jl_var_1044, lay_var_997) * (o3m2 - o3m1)
        taug_var_959(jl_var_1044, 52 + ig_var_995, lay_var_997) = tau_major_var_1022(ig_var_995) + tau_major1_var_1023(ig_var_995) + tauself_var_1021 + taufor_var_1020 + abso3_var_1024 * colo3_var_972(jl_var_1044, lay_var_997) + wx_var_960(jl_var_1044, 1, lay_var_997) * ccl4(ig_var_995)
        fracs_var_977(jl_var_1044, 52 + ig_var_995, lay_var_997) = fracrefa_var_173(ig_var_995, jpl_var_999) + fpl_var_1031 * (fracrefa_var_173(ig_var_995, jpl_var_999 + 1) - fracrefa_var_173(ig_var_995, jpl_var_999))
      END DO
    END DO
  END DO
  DO lay_var_997 = laytrop_max_var_1035 + 1, klev_var_958
    DO jl_var_1044 = kidia_var_956, kfdia_var_957
      speccomb_var_987 = colo3_var_972(jl_var_1044, lay_var_997) + rat_o3co2_var_980(jl_var_1044, lay_var_997) * colco2_var_971(jl_var_1044, lay_var_997)
      specparm_var_1027 = MIN(colo3_var_972(jl_var_1044, lay_var_997) / speccomb_var_987, 0.999999D0)
      specmult_var_1026 = 4.0D0 * (specparm_var_1027)
      js_var_996 = 1 + INT(specmult_var_1026)
      fs_var_1025 = ((specmult_var_1026) - AINT((specmult_var_1026)))
      speccomb1_var_988 = colo3_var_972(jl_var_1044, lay_var_997) + rat_o3co2_1_var_981(jl_var_1044, lay_var_997) * colco2_var_971(jl_var_1044, lay_var_997)
      specparm1_var_1030 = MIN(colo3_var_972(jl_var_1044, lay_var_997) / speccomb1_var_988, 0.999999D0)
      specmult1_var_1029 = 4.0D0 * (specparm1_var_1030)
      js1_var_998 = 1 + INT(specmult1_var_1029)
      fs1_var_1028 = ((specmult1_var_1029) - AINT((specmult1_var_1029)))
      fac000_var_1003 = (1.0D0 - fs_var_1025) * fac00_var_962(jl_var_1044, lay_var_997)
      fac010_var_1006 = (1.0D0 - fs_var_1025) * fac10_var_964(jl_var_1044, lay_var_997)
      fac100_var_1004 = fs_var_1025 * fac00_var_962(jl_var_1044, lay_var_997)
      fac110_var_1007 = fs_var_1025 * fac10_var_964(jl_var_1044, lay_var_997)
      fac001_var_1009 = (1.0D0 - fs1_var_1028) * fac01_var_963(jl_var_1044, lay_var_997)
      fac011_var_1012 = (1.0D0 - fs1_var_1028) * fac11_var_965(jl_var_1044, lay_var_997)
      fac101_var_1010 = fs1_var_1028 * fac01_var_963(jl_var_1044, lay_var_997)
      fac111_var_1013 = fs1_var_1028 * fac11_var_965(jl_var_1044, lay_var_997)
      speccomb_planck_var_989 = colo3_var_972(jl_var_1044, lay_var_997) + refrat_planck_b_var_1001 * colco2_var_971(jl_var_1044, lay_var_997)
      specparm_planck_var_1033 = MIN(colo3_var_972(jl_var_1044, lay_var_997) / speccomb_planck_var_989, 0.999999D0)
      specmult_planck_var_1032 = 4.0D0 * specparm_planck_var_1033
      jpl_var_999 = 1 + INT(specmult_planck_var_1032)
      fpl_var_1031 = ((specmult_planck_var_1032) - AINT((specmult_planck_var_1032)))
      ind0_var_990 = ((jp_var_966(jl_var_1044, lay_var_997) - 13) * 5 + (jt_var_967(jl_var_1044, lay_var_997) - 1)) * nspb_var_217(5) + js_var_996
      ind1_var_991 = ((jp_var_966(jl_var_1044, lay_var_997) - 12) * 5 + (jt1_var_968(jl_var_1044, lay_var_997) - 1)) * nspb_var_217(5) + js1_var_998
      DO ig_var_995 = 1, 16
        taug_var_959(jl_var_1044, 52 + ig_var_995, lay_var_997) = speccomb_var_987 * (fac000_var_1003 * absb_var_176(ind0_var_990, ig_var_995) + fac100_var_1004 * absb_var_176(ind0_var_990 + 1, ig_var_995) + fac010_var_1006 * absb_var_176(ind0_var_990 + 5, ig_var_995) + fac110_var_1007 * absb_var_176(ind0_var_990 + 6, ig_var_995)) + speccomb1_var_988 * (fac001_var_1009 * absb_var_176(ind1_var_991, ig_var_995) + fac101_var_1010 * absb_var_176(ind1_var_991 + 1, ig_var_995) + fac011_var_1012 * absb_var_176(ind1_var_991 + 5, ig_var_995) + fac111_var_1013 * absb_var_176(ind1_var_991 + 6, ig_var_995)) + wx_var_960(jl_var_1044, 1, lay_var_997) * ccl4(ig_var_995)
        fracs_var_977(jl_var_1044, 52 + ig_var_995, lay_var_997) = fracrefb_var_174(ig_var_995, jpl_var_999) + fpl_var_1031 * (fracrefb_var_174(ig_var_995, jpl_var_999 + 1) - fracrefb_var_174(ig_var_995, jpl_var_999))
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1035 /= laytrop_min_var_1034) THEN
    DO lay_var_997 = laytrop_min_var_1034 + 1, laytrop_max_var_1035
      ixc0_var_1041 = ixc_var_1036(lay_var_997)
      DO ixp_var_1042 = 1, ixc0_var_1041
        jl_var_1044 = ixlow_var_1037(ixp_var_1042, lay_var_997)
        speccomb_var_987 = colh2o_var_970(jl_var_1044, lay_var_997) + rat_h2oco2_var_978(jl_var_1044, lay_var_997) * colco2_var_971(jl_var_1044, lay_var_997)
        specparm_var_1027 = MIN(colh2o_var_970(jl_var_1044, lay_var_997) / speccomb_var_987, 0.999999D0)
        specmult_var_1026 = 8.0D0 * (specparm_var_1027)
        js_var_996 = 1 + INT(specmult_var_1026)
        fs_var_1025 = ((specmult_var_1026) - AINT((specmult_var_1026)))
        speccomb1_var_988 = colh2o_var_970(jl_var_1044, lay_var_997) + rat_h2oco2_1_var_979(jl_var_1044, lay_var_997) * colco2_var_971(jl_var_1044, lay_var_997)
        specparm1_var_1030 = MIN(colh2o_var_970(jl_var_1044, lay_var_997) / speccomb1_var_988, 0.999999D0)
        specmult1_var_1029 = 8.0D0 * (specparm1_var_1030)
        js1_var_998 = 1 + INT(specmult1_var_1029)
        fs1_var_1028 = ((specmult1_var_1029) - AINT((specmult1_var_1029)))
        speccomb_mo3 = colh2o_var_970(jl_var_1044, lay_var_997) + refrat_m_a_var_1002 * colco2_var_971(jl_var_1044, lay_var_997)
        specparm_mo3 = MIN(colh2o_var_970(jl_var_1044, lay_var_997) / speccomb_mo3, 0.999999D0)
        specmult_mo3 = 8.0D0 * specparm_mo3
        jmo3 = 1 + INT(specmult_mo3)
        fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))
        speccomb_planck_var_989 = colh2o_var_970(jl_var_1044, lay_var_997) + refrat_planck_a_var_1000 * colco2_var_971(jl_var_1044, lay_var_997)
        specparm_planck_var_1033 = MIN(colh2o_var_970(jl_var_1044, lay_var_997) / speccomb_planck_var_989, 0.999999D0)
        specmult_planck_var_1032 = 8.0D0 * specparm_planck_var_1033
        jpl_var_999 = 1 + INT(specmult_planck_var_1032)
        fpl_var_1031 = ((specmult_planck_var_1032) - AINT((specmult_planck_var_1032)))
        ind0_var_990 = ((jp_var_966(jl_var_1044, lay_var_997) - 1) * 5 + (jt_var_967(jl_var_1044, lay_var_997) - 1)) * nspa_var_216(5) + js_var_996
        ind1_var_991 = (jp_var_966(jl_var_1044, lay_var_997) * 5 + (jt1_var_968(jl_var_1044, lay_var_997) - 1)) * nspa_var_216(5) + js1_var_998
        inds_var_992 = indself_var_976(jl_var_1044, lay_var_997)
        indf_var_993 = indfor_var_982(jl_var_1044, lay_var_997)
        indm_var_994 = indminor_var_986(jl_var_1044, lay_var_997)
        IF (specparm_var_1027 .LT. 0.125D0) THEN
          p_var_1015 = fs_var_1025 - 1.0D0
          p4_var_1016 = p_var_1015 ** 4
          fk0_var_1017 = p4_var_1016
          fk1_var_1018 = 1.0D0 - p_var_1015 - 2.0D0 * p4_var_1016
          fk2_var_1019 = p_var_1015 + p4_var_1016
          fac000_var_1003 = fk0_var_1017 * fac00_var_962(jl_var_1044, lay_var_997)
          fac100_var_1004 = fk1_var_1018 * fac00_var_962(jl_var_1044, lay_var_997)
          fac200_var_1005 = fk2_var_1019 * fac00_var_962(jl_var_1044, lay_var_997)
          fac010_var_1006 = fk0_var_1017 * fac10_var_964(jl_var_1044, lay_var_997)
          fac110_var_1007 = fk1_var_1018 * fac10_var_964(jl_var_1044, lay_var_997)
          fac210_var_1008 = fk2_var_1019 * fac10_var_964(jl_var_1044, lay_var_997)
        ELSE IF (specparm_var_1027 .GT. 0.875D0) THEN
          p_var_1015 = - fs_var_1025
          p4_var_1016 = p_var_1015 ** 4
          fk0_var_1017 = p4_var_1016
          fk1_var_1018 = 1.0D0 - p_var_1015 - 2.0D0 * p4_var_1016
          fk2_var_1019 = p_var_1015 + p4_var_1016
          fac000_var_1003 = fk0_var_1017 * fac00_var_962(jl_var_1044, lay_var_997)
          fac100_var_1004 = fk1_var_1018 * fac00_var_962(jl_var_1044, lay_var_997)
          fac200_var_1005 = fk2_var_1019 * fac00_var_962(jl_var_1044, lay_var_997)
          fac010_var_1006 = fk0_var_1017 * fac10_var_964(jl_var_1044, lay_var_997)
          fac110_var_1007 = fk1_var_1018 * fac10_var_964(jl_var_1044, lay_var_997)
          fac210_var_1008 = fk2_var_1019 * fac10_var_964(jl_var_1044, lay_var_997)
        ELSE
          fac000_var_1003 = (1.0D0 - fs_var_1025) * fac00_var_962(jl_var_1044, lay_var_997)
          fac010_var_1006 = (1.0D0 - fs_var_1025) * fac10_var_964(jl_var_1044, lay_var_997)
          fac100_var_1004 = fs_var_1025 * fac00_var_962(jl_var_1044, lay_var_997)
          fac110_var_1007 = fs_var_1025 * fac10_var_964(jl_var_1044, lay_var_997)
          fac200_var_1005 = 0.0D0
          fac210_var_1008 = 0.0D0
        END IF
        IF (specparm1_var_1030 .LT. 0.125D0) THEN
          p_var_1015 = fs1_var_1028 - 1.0D0
          p4_var_1016 = p_var_1015 ** 4
          fk0_var_1017 = p4_var_1016
          fk1_var_1018 = 1.0D0 - p_var_1015 - 2.0D0 * p4_var_1016
          fk2_var_1019 = p_var_1015 + p4_var_1016
          fac001_var_1009 = fk0_var_1017 * fac01_var_963(jl_var_1044, lay_var_997)
          fac101_var_1010 = fk1_var_1018 * fac01_var_963(jl_var_1044, lay_var_997)
          fac201_var_1011 = fk2_var_1019 * fac01_var_963(jl_var_1044, lay_var_997)
          fac011_var_1012 = fk0_var_1017 * fac11_var_965(jl_var_1044, lay_var_997)
          fac111_var_1013 = fk1_var_1018 * fac11_var_965(jl_var_1044, lay_var_997)
          fac211_var_1014 = fk2_var_1019 * fac11_var_965(jl_var_1044, lay_var_997)
        ELSE IF (specparm1_var_1030 .GT. 0.875D0) THEN
          p_var_1015 = - fs1_var_1028
          p4_var_1016 = p_var_1015 ** 4
          fk0_var_1017 = p4_var_1016
          fk1_var_1018 = 1.0D0 - p_var_1015 - 2.0D0 * p4_var_1016
          fk2_var_1019 = p_var_1015 + p4_var_1016
          fac001_var_1009 = fk0_var_1017 * fac01_var_963(jl_var_1044, lay_var_997)
          fac101_var_1010 = fk1_var_1018 * fac01_var_963(jl_var_1044, lay_var_997)
          fac201_var_1011 = fk2_var_1019 * fac01_var_963(jl_var_1044, lay_var_997)
          fac011_var_1012 = fk0_var_1017 * fac11_var_965(jl_var_1044, lay_var_997)
          fac111_var_1013 = fk1_var_1018 * fac11_var_965(jl_var_1044, lay_var_997)
          fac211_var_1014 = fk2_var_1019 * fac11_var_965(jl_var_1044, lay_var_997)
        ELSE
          fac001_var_1009 = (1.0D0 - fs1_var_1028) * fac01_var_963(jl_var_1044, lay_var_997)
          fac011_var_1012 = (1.0D0 - fs1_var_1028) * fac11_var_965(jl_var_1044, lay_var_997)
          fac101_var_1010 = fs1_var_1028 * fac01_var_963(jl_var_1044, lay_var_997)
          fac111_var_1013 = fs1_var_1028 * fac11_var_965(jl_var_1044, lay_var_997)
          fac201_var_1011 = 0.0D0
          fac211_var_1014 = 0.0D0
        END IF
        IF (specparm_var_1027 .LT. 0.125D0) THEN
          tau_major_var_1022(1 : ng5) = speccomb_var_987 * (fac000_var_1003 * absa_var_175(ind0_var_990, 1 : 16) + fac100_var_1004 * absa_var_175(ind0_var_990 + 1, 1 : 16) + fac200_var_1005 * absa_var_175(ind0_var_990 + 2, 1 : 16) + fac010_var_1006 * absa_var_175(ind0_var_990 + 9, 1 : 16) + fac110_var_1007 * absa_var_175(ind0_var_990 + 10, 1 : 16) + fac210_var_1008 * absa_var_175(ind0_var_990 + 11, 1 : 16))
        ELSE IF (specparm_var_1027 .GT. 0.875D0) THEN
          tau_major_var_1022(1 : ng5) = speccomb_var_987 * (fac200_var_1005 * absa_var_175(ind0_var_990 - 1, 1 : 16) + fac100_var_1004 * absa_var_175(ind0_var_990, 1 : 16) + fac000_var_1003 * absa_var_175(ind0_var_990 + 1, 1 : 16) + fac210_var_1008 * absa_var_175(ind0_var_990 + 8, 1 : 16) + fac110_var_1007 * absa_var_175(ind0_var_990 + 9, 1 : 16) + fac010_var_1006 * absa_var_175(ind0_var_990 + 10, 1 : 16))
        ELSE
          tau_major_var_1022(1 : ng5) = speccomb_var_987 * (fac000_var_1003 * absa_var_175(ind0_var_990, 1 : 16) + fac100_var_1004 * absa_var_175(ind0_var_990 + 1, 1 : 16) + fac010_var_1006 * absa_var_175(ind0_var_990 + 9, 1 : 16) + fac110_var_1007 * absa_var_175(ind0_var_990 + 10, 1 : 16))
        END IF
        IF (specparm1_var_1030 .LT. 0.125D0) THEN
          tau_major1_var_1023(1 : ng5) = speccomb1_var_988 * (fac001_var_1009 * absa_var_175(ind1_var_991, 1 : 16) + fac101_var_1010 * absa_var_175(ind1_var_991 + 1, 1 : 16) + fac201_var_1011 * absa_var_175(ind1_var_991 + 2, 1 : 16) + fac011_var_1012 * absa_var_175(ind1_var_991 + 9, 1 : 16) + fac111_var_1013 * absa_var_175(ind1_var_991 + 10, 1 : 16) + fac211_var_1014 * absa_var_175(ind1_var_991 + 11, 1 : 16))
        ELSE IF (specparm1_var_1030 .GT. 0.875D0) THEN
          tau_major1_var_1023(1 : ng5) = speccomb1_var_988 * (fac201_var_1011 * absa_var_175(ind1_var_991 - 1, 1 : 16) + fac101_var_1010 * absa_var_175(ind1_var_991, 1 : 16) + fac001_var_1009 * absa_var_175(ind1_var_991 + 1, 1 : 16) + fac211_var_1014 * absa_var_175(ind1_var_991 + 8, 1 : 16) + fac111_var_1013 * absa_var_175(ind1_var_991 + 9, 1 : 16) + fac011_var_1012 * absa_var_175(ind1_var_991 + 10, 1 : 16))
        ELSE
          tau_major1_var_1023(1 : ng5) = speccomb1_var_988 * (fac001_var_1009 * absa_var_175(ind1_var_991, 1 : 16) + fac101_var_1010 * absa_var_175(ind1_var_991 + 1, 1 : 16) + fac011_var_1012 * absa_var_175(ind1_var_991 + 9, 1 : 16) + fac111_var_1013 * absa_var_175(ind1_var_991 + 10, 1 : 16))
        END IF
        DO ig_var_995 = 1, 16
          tauself_var_1021 = selffac_var_974(jl_var_1044, lay_var_997) * (selfref_var_178(inds_var_992, ig_var_995) + selffrac_var_975(jl_var_1044, lay_var_997) * (selfref_var_178(inds_var_992 + 1, ig_var_995) - selfref_var_178(inds_var_992, ig_var_995)))
          taufor_var_1020 = forfac_var_984(jl_var_1044, lay_var_997) * (forref_var_179(indf_var_993, ig_var_995) + forfrac_var_983(jl_var_1044, lay_var_997) * (forref_var_179(indf_var_993 + 1, ig_var_995) - forref_var_179(indf_var_993, ig_var_995)))
          o3m1 = ka_mo3_var_177(jmo3, indm_var_994, ig_var_995) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_994, ig_var_995) - ka_mo3_var_177(jmo3, indm_var_994, ig_var_995))
          o3m2 = ka_mo3_var_177(jmo3, indm_var_994 + 1, ig_var_995) + fmo3 * (ka_mo3_var_177(jmo3 + 1, indm_var_994 + 1, ig_var_995) - ka_mo3_var_177(jmo3, indm_var_994 + 1, ig_var_995))
          abso3_var_1024 = o3m1 + minorfrac_var_985(jl_var_1044, lay_var_997) * (o3m2 - o3m1)
          taug_var_959(jl_var_1044, 52 + ig_var_995, lay_var_997) = tau_major_var_1022(ig_var_995) + tau_major1_var_1023(ig_var_995) + tauself_var_1021 + taufor_var_1020 + abso3_var_1024 * colo3_var_972(jl_var_1044, lay_var_997) + wx_var_960(jl_var_1044, 1, lay_var_997) * ccl4(ig_var_995)
          fracs_var_977(jl_var_1044, 52 + ig_var_995, lay_var_997) = fracrefa_var_173(ig_var_995, jpl_var_999) + fpl_var_1031 * (fracrefa_var_173(ig_var_995, jpl_var_999 + 1) - fracrefa_var_173(ig_var_995, jpl_var_999))
        END DO
      END DO
      ixc0_var_1041 = kfdia_var_957 - kidia_var_956 + 1 - ixc0_var_1041
      DO ixp_var_1042 = 1, ixc0_var_1041
        jl_var_1044 = ixhigh_var_1038(ixp_var_1042, lay_var_997)
        speccomb_var_987 = colo3_var_972(jl_var_1044, lay_var_997) + rat_o3co2_var_980(jl_var_1044, lay_var_997) * colco2_var_971(jl_var_1044, lay_var_997)
        specparm_var_1027 = MIN(colo3_var_972(jl_var_1044, lay_var_997) / speccomb_var_987, 0.999999D0)
        specmult_var_1026 = 4.0D0 * (specparm_var_1027)
        js_var_996 = 1 + INT(specmult_var_1026)
        fs_var_1025 = ((specmult_var_1026) - AINT((specmult_var_1026)))
        speccomb1_var_988 = colo3_var_972(jl_var_1044, lay_var_997) + rat_o3co2_1_var_981(jl_var_1044, lay_var_997) * colco2_var_971(jl_var_1044, lay_var_997)
        specparm1_var_1030 = MIN(colo3_var_972(jl_var_1044, lay_var_997) / speccomb1_var_988, 0.999999D0)
        specmult1_var_1029 = 4.0D0 * (specparm1_var_1030)
        js1_var_998 = 1 + INT(specmult1_var_1029)
        fs1_var_1028 = ((specmult1_var_1029) - AINT((specmult1_var_1029)))
        fac000_var_1003 = (1.0D0 - fs_var_1025) * fac00_var_962(jl_var_1044, lay_var_997)
        fac010_var_1006 = (1.0D0 - fs_var_1025) * fac10_var_964(jl_var_1044, lay_var_997)
        fac100_var_1004 = fs_var_1025 * fac00_var_962(jl_var_1044, lay_var_997)
        fac110_var_1007 = fs_var_1025 * fac10_var_964(jl_var_1044, lay_var_997)
        fac001_var_1009 = (1.0D0 - fs1_var_1028) * fac01_var_963(jl_var_1044, lay_var_997)
        fac011_var_1012 = (1.0D0 - fs1_var_1028) * fac11_var_965(jl_var_1044, lay_var_997)
        fac101_var_1010 = fs1_var_1028 * fac01_var_963(jl_var_1044, lay_var_997)
        fac111_var_1013 = fs1_var_1028 * fac11_var_965(jl_var_1044, lay_var_997)
        speccomb_planck_var_989 = colo3_var_972(jl_var_1044, lay_var_997) + refrat_planck_b_var_1001 * colco2_var_971(jl_var_1044, lay_var_997)
        specparm_planck_var_1033 = MIN(colo3_var_972(jl_var_1044, lay_var_997) / speccomb_planck_var_989, 0.999999D0)
        specmult_planck_var_1032 = 4.0D0 * specparm_planck_var_1033
        jpl_var_999 = 1 + INT(specmult_planck_var_1032)
        fpl_var_1031 = ((specmult_planck_var_1032) - AINT((specmult_planck_var_1032)))
        ind0_var_990 = ((jp_var_966(jl_var_1044, lay_var_997) - 13) * 5 + (jt_var_967(jl_var_1044, lay_var_997) - 1)) * nspb_var_217(5) + js_var_996
        ind1_var_991 = ((jp_var_966(jl_var_1044, lay_var_997) - 12) * 5 + (jt1_var_968(jl_var_1044, lay_var_997) - 1)) * nspb_var_217(5) + js1_var_998
        DO ig_var_995 = 1, 16
          taug_var_959(jl_var_1044, 52 + ig_var_995, lay_var_997) = speccomb_var_987 * (fac000_var_1003 * absb_var_176(ind0_var_990, ig_var_995) + fac100_var_1004 * absb_var_176(ind0_var_990 + 1, ig_var_995) + fac010_var_1006 * absb_var_176(ind0_var_990 + 5, ig_var_995) + fac110_var_1007 * absb_var_176(ind0_var_990 + 6, ig_var_995)) + speccomb1_var_988 * (fac001_var_1009 * absb_var_176(ind1_var_991, ig_var_995) + fac101_var_1010 * absb_var_176(ind1_var_991 + 1, ig_var_995) + fac011_var_1012 * absb_var_176(ind1_var_991 + 5, ig_var_995) + fac111_var_1013 * absb_var_176(ind1_var_991 + 6, ig_var_995)) + wx_var_960(jl_var_1044, 1, lay_var_997) * ccl4(ig_var_995)
          fracs_var_977(jl_var_1044, 52 + ig_var_995, lay_var_997) = fracrefb_var_174(ig_var_995, jpl_var_999) + fpl_var_1031 * (fracrefb_var_174(ig_var_995, jpl_var_999 + 1) - fracrefb_var_174(ig_var_995, jpl_var_999))
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol5
SUBROUTINE srtm_setcoef(kidia_var_1045, kfdia_var_1046, klev_var_1047, pavel_var_1048, ptavel_var_1049, pcoldry_var_1050, pwkl_var_1051, klaytrop_var_1052, pcolch4_var_1053, pcolco2_var_1054, pcolh2o_var_1055, pcolmol_var_1056, pcolo2_var_1057, pcolo3_var_1058, pforfac_var_1059, pforfrac_var_1060, kindfor_var_1061, pselffac_var_1062, pselffrac_var_1063, kindself_var_1064, pfac00_var_1065, pfac01_var_1066, pfac10_var_1067, pfac11_var_1068, kjp_var_1069, kjt_var_1070, kjt1_var_1071, prmu0_var_1072)
  USE yoesrtwn, ONLY: preflog_var_315, tref_var_316
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1045, kfdia_var_1046
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1047
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1048(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(IN) :: ptavel_var_1049(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(IN) :: pcoldry_var_1050(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(IN) :: pwkl_var_1051(kidia_var_1045 : kfdia_var_1046, 35, klev_var_1047)
  INTEGER(KIND = 4), INTENT(INOUT) :: klaytrop_var_1052(kidia_var_1045 : kfdia_var_1046)
  REAL(KIND = 8), INTENT(INOUT) :: pcolch4_var_1053(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pcolco2_var_1054(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pcolh2o_var_1055(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pcolmol_var_1056(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pcolo2_var_1057(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pcolo3_var_1058(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pforfac_var_1059(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pforfrac_var_1060(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  INTEGER(KIND = 4), INTENT(INOUT) :: kindfor_var_1061(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pselffac_var_1062(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pselffrac_var_1063(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  INTEGER(KIND = 4), INTENT(INOUT) :: kindself_var_1064(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pfac00_var_1065(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pfac01_var_1066(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pfac10_var_1067(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(INOUT) :: pfac11_var_1068(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjp_var_1069(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjt_var_1070(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjt1_var_1071(kidia_var_1045 : kfdia_var_1046, klev_var_1047)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1072(kidia_var_1045 : kfdia_var_1046)
  INTEGER(KIND = 4) :: i_nlayers_var_1073, jk_var_1074, jl_var_1075, jp1_var_1076
  REAL(KIND = 8) :: z_stpfac_var_1077, z_plog_var_1078
  REAL(KIND = 8) :: z_fp_var_1079, z_ft_var_1080, z_ft1_var_1081, z_water_var_1082, z_scalefac_var_1083
  REAL(KIND = 8) :: z_factor_var_1084, z_co2reg_var_1085, z_compfp_var_1086
  LOGICAL :: goto_0 = .FALSE.
  LOGICAL :: goto_1 = .FALSE.
  z_stpfac_var_1077 = 0.29220138203356366D0
  i_nlayers_var_1073 = klev_var_1047
  DO jk_var_1074 = 1, klev_var_1047
    DO jl_var_1075 = kidia_var_1045, kfdia_var_1046
      pcolmol_var_1056(jl_var_1075, jk_var_1074) = 0.0D0
    END DO
  END DO
  DO jl_var_1075 = kidia_var_1045, kfdia_var_1046
    IF (prmu0_var_1072(jl_var_1075) > 0.0D0) THEN
      klaytrop_var_1052(jl_var_1075) = 0
    END IF
  END DO
  DO jk_var_1074 = 1, i_nlayers_var_1073
    DO jl_var_1075 = kidia_var_1045, kfdia_var_1046
      IF (prmu0_var_1072(jl_var_1075) > 0.0D0) THEN
        z_plog_var_1078 = LOG(pavel_var_1048(jl_var_1075, jk_var_1074))
        kjp_var_1069(jl_var_1075, jk_var_1074) = INT(36.0D0 - 5.0D0 * (z_plog_var_1078 + 0.04D0))
        IF (kjp_var_1069(jl_var_1075, jk_var_1074) < 1) THEN
          kjp_var_1069(jl_var_1075, jk_var_1074) = 1
        ELSE IF (kjp_var_1069(jl_var_1075, jk_var_1074) > 58) THEN
          kjp_var_1069(jl_var_1075, jk_var_1074) = 58
        END IF
        jp1_var_1076 = kjp_var_1069(jl_var_1075, jk_var_1074) + 1
        z_fp_var_1079 = 5.0 * (preflog_var_315(kjp_var_1069(jl_var_1075, jk_var_1074)) - z_plog_var_1078)
        kjt_var_1070(jl_var_1075, jk_var_1074) = INT(3.0 + (ptavel_var_1049(jl_var_1075, jk_var_1074) - tref_var_316(kjp_var_1069(jl_var_1075, jk_var_1074))) / 15.0)
        IF (kjt_var_1070(jl_var_1075, jk_var_1074) < 1) THEN
          kjt_var_1070(jl_var_1075, jk_var_1074) = 1
        ELSE IF (kjt_var_1070(jl_var_1075, jk_var_1074) > 4) THEN
          kjt_var_1070(jl_var_1075, jk_var_1074) = 4
        END IF
        z_ft_var_1080 = ((ptavel_var_1049(jl_var_1075, jk_var_1074) - tref_var_316(kjp_var_1069(jl_var_1075, jk_var_1074))) / 15.0) - REAL(kjt_var_1070(jl_var_1075, jk_var_1074) - 3)
        kjt1_var_1071(jl_var_1075, jk_var_1074) = INT(3.0 + (ptavel_var_1049(jl_var_1075, jk_var_1074) - tref_var_316(jp1_var_1076)) / 15.0)
        IF (kjt1_var_1071(jl_var_1075, jk_var_1074) < 1) THEN
          kjt1_var_1071(jl_var_1075, jk_var_1074) = 1
        ELSE IF (kjt1_var_1071(jl_var_1075, jk_var_1074) > 4) THEN
          kjt1_var_1071(jl_var_1075, jk_var_1074) = 4
        END IF
        z_ft1_var_1081 = ((ptavel_var_1049(jl_var_1075, jk_var_1074) - tref_var_316(jp1_var_1076)) / 15.0) - REAL(kjt1_var_1071(jl_var_1075, jk_var_1074) - 3)
        z_water_var_1082 = pwkl_var_1051(jl_var_1075, 1, jk_var_1074) / pcoldry_var_1050(jl_var_1075, jk_var_1074)
        z_scalefac_var_1083 = pavel_var_1048(jl_var_1075, jk_var_1074) * z_stpfac_var_1077 / ptavel_var_1049(jl_var_1075, jk_var_1074)
        IF (z_plog_var_1078 <= 4.56D0) goto_0 = .TRUE.
        IF (.NOT. (goto_0)) klaytrop_var_1052(jl_var_1075) = klaytrop_var_1052(jl_var_1075) + 1
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pforfac_var_1059(jl_var_1075, jk_var_1074) = z_scalefac_var_1083 / (1.0 + z_water_var_1082)
        IF (.NOT. (goto_0)) z_factor_var_1084 = (332.0 - ptavel_var_1049(jl_var_1075, jk_var_1074)) / 36.0
        IF (.NOT. (goto_0)) kindfor_var_1061(jl_var_1075, jk_var_1074) = MIN(2, MAX(1, INT(z_factor_var_1084)))
        IF (.NOT. (goto_0)) pforfrac_var_1060(jl_var_1075, jk_var_1074) = z_factor_var_1084 - REAL(kindfor_var_1061(jl_var_1075, jk_var_1074))
        IF (.NOT. (goto_0)) pselffac_var_1062(jl_var_1075, jk_var_1074) = z_water_var_1082 * pforfac_var_1059(jl_var_1075, jk_var_1074)
        IF (.NOT. (goto_0)) z_factor_var_1084 = (ptavel_var_1049(jl_var_1075, jk_var_1074) - 188.0) / 7.2
        IF (.NOT. (goto_0)) kindself_var_1064(jl_var_1075, jk_var_1074) = MIN(9, MAX(1, INT(z_factor_var_1084) - 7))
        IF (.NOT. (goto_0)) pselffrac_var_1063(jl_var_1075, jk_var_1074) = z_factor_var_1084 - REAL(kindself_var_1064(jl_var_1075, jk_var_1074) + 7)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolh2o_var_1055(jl_var_1075, jk_var_1074) = 1E-20 * pwkl_var_1051(jl_var_1075, 1, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolco2_var_1054(jl_var_1075, jk_var_1074) = 1E-20 * pwkl_var_1051(jl_var_1075, 2, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolo3_var_1058(jl_var_1075, jk_var_1074) = 1E-20 * pwkl_var_1051(jl_var_1075, 3, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolch4_var_1053(jl_var_1075, jk_var_1074) = 1E-20 * pwkl_var_1051(jl_var_1075, 6, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolo2_var_1057(jl_var_1075, jk_var_1074) = 1E-20 * pwkl_var_1051(jl_var_1075, 7, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolmol_var_1056(jl_var_1075, jk_var_1074) = 1E-20 * pcoldry_var_1050(jl_var_1075, jk_var_1074) + pcolh2o_var_1055(jl_var_1075, jk_var_1074)
        IF (.NOT. (goto_0) .AND. pcolco2_var_1054(jl_var_1075, jk_var_1074) == 0.0) pcolco2_var_1054(jl_var_1075, jk_var_1074) = 1E-32 * pcoldry_var_1050(jl_var_1075, jk_var_1074)
        IF (.NOT. (goto_0) .AND. pcolch4_var_1053(jl_var_1075, jk_var_1074) == 0.0) pcolch4_var_1053(jl_var_1075, jk_var_1074) = 1E-32 * pcoldry_var_1050(jl_var_1075, jk_var_1074)
        IF (.NOT. (goto_0) .AND. pcolo2_var_1057(jl_var_1075, jk_var_1074) == 0.0) pcolo2_var_1057(jl_var_1075, jk_var_1074) = 1E-32 * pcoldry_var_1050(jl_var_1075, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) z_co2reg_var_1085 = 3.55E-24 * pcoldry_var_1050(jl_var_1075, jk_var_1074)
        IF (.NOT. (goto_0)) goto_1 = .TRUE.
5300    CONTINUE
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pforfac_var_1059(jl_var_1075, jk_var_1074) = z_scalefac_var_1083 / (1.0 + z_water_var_1082)
        IF (.NOT. (goto_1)) z_factor_var_1084 = (ptavel_var_1049(jl_var_1075, jk_var_1074) - 188.0) / 36.0
        IF (.NOT. (goto_1)) kindfor_var_1061(jl_var_1075, jk_var_1074) = 3
        IF (.NOT. (goto_1)) pforfrac_var_1060(jl_var_1075, jk_var_1074) = z_factor_var_1084 - 1.0
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolh2o_var_1055(jl_var_1075, jk_var_1074) = 1E-20 * pwkl_var_1051(jl_var_1075, 1, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolco2_var_1054(jl_var_1075, jk_var_1074) = 1E-20 * pwkl_var_1051(jl_var_1075, 2, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolo3_var_1058(jl_var_1075, jk_var_1074) = 1E-20 * pwkl_var_1051(jl_var_1075, 3, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolch4_var_1053(jl_var_1075, jk_var_1074) = 1E-20 * pwkl_var_1051(jl_var_1075, 6, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolo2_var_1057(jl_var_1075, jk_var_1074) = 1E-20 * pwkl_var_1051(jl_var_1075, 7, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) pcolmol_var_1056(jl_var_1075, jk_var_1074) = 1E-20 * pcoldry_var_1050(jl_var_1075, jk_var_1074) + pcolh2o_var_1055(jl_var_1075, jk_var_1074)
        IF (.NOT. (goto_1) .AND. pcolco2_var_1054(jl_var_1075, jk_var_1074) == 0.0) pcolco2_var_1054(jl_var_1075, jk_var_1074) = 1E-32 * pcoldry_var_1050(jl_var_1075, jk_var_1074)
        IF (.NOT. (goto_1) .AND. pcolch4_var_1053(jl_var_1075, jk_var_1074) == 0.0) pcolch4_var_1053(jl_var_1075, jk_var_1074) = 1E-32 * pcoldry_var_1050(jl_var_1075, jk_var_1074)
        IF (.NOT. (goto_1) .AND. pcolo2_var_1057(jl_var_1075, jk_var_1074) == 0.0) pcolo2_var_1057(jl_var_1075, jk_var_1074) = 1E-32 * pcoldry_var_1050(jl_var_1075, jk_var_1074)
        IF (.NOT. (goto_1) .AND. .NOT. (goto_0)) z_co2reg_var_1085 = 3.55E-24 * pcoldry_var_1050(jl_var_1075, jk_var_1074)
        IF (.NOT. (goto_1)) pselffac_var_1062(jl_var_1075, jk_var_1074) = 0.0D0
        IF (.NOT. (goto_1)) pselffrac_var_1063(jl_var_1075, jk_var_1074) = 0.0D0
        IF (.NOT. (goto_1)) kindself_var_1064(jl_var_1075, jk_var_1074) = 0
5400    CONTINUE
        z_compfp_var_1086 = 1.0 - z_fp_var_1079
        pfac10_var_1067(jl_var_1075, jk_var_1074) = z_compfp_var_1086 * z_ft_var_1080
        pfac00_var_1065(jl_var_1075, jk_var_1074) = z_compfp_var_1086 * (1.0 - z_ft_var_1080)
        pfac11_var_1068(jl_var_1075, jk_var_1074) = z_fp_var_1079 * z_ft1_var_1081
        pfac01_var_1066(jl_var_1075, jk_var_1074) = z_fp_var_1079 * (1.0 - z_ft1_var_1081)
9000    FORMAT(1X, 2I3, 3I4, F6.1, 4F7.2, 12E9.2, 2I5)
      END IF
    END DO
  END DO
END SUBROUTINE srtm_setcoef
SUBROUTINE rrtm_taumol4(kidia_var_1087, kfdia_var_1088, klev_var_1089, taug_var_1090, p_tauaerl_var_1091, fac00_var_1092, fac01_var_1093, fac10_var_1094, fac11_var_1095, forfac_var_1113, forfrac_var_1114, indfor_var_1112, jp_var_1096, jt_var_1097, jt1_var_1098, oneminus_var_1099, colh2o_var_1100, colco2_var_1101, colo3_var_1102, laytrop_var_1103, selffac_var_1104, selffrac_var_1105, indself_var_1106, fracs_var_1107, rat_h2oco2_var_1108, rat_h2oco2_1_var_1109, rat_o3co2_var_1110, rat_o3co2_1_var_1111)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng4
  USE yoerrta4, ONLY: absa_var_169, absb_var_170, forref_var_172, fracrefa_var_167, fracrefb_var_168, selfref_var_171
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1087
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1088
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1089
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1090(kidia_var_1087 : kfdia_var_1088, 140, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1091(kidia_var_1087 : kfdia_var_1088, klev_var_1089, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1092(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1093(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1094(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1095(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1096(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1097(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1098(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1099
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1100(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1101(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1102(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1103(kidia_var_1087 : kfdia_var_1088)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1104(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1105(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1106(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1107(kidia_var_1087 : kfdia_var_1088, 140, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1108(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1109(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_var_1110(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_1_var_1111(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1112(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1113(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1114(kidia_var_1087 : kfdia_var_1088, klev_var_1089)
  REAL(KIND = 8) :: speccomb_var_1115, speccomb1_var_1116, speccomb_planck_var_1117
  INTEGER(KIND = 4) :: ind0_var_1118, ind1_var_1119, inds_var_1120, indf_var_1121
  INTEGER(KIND = 4) :: ig_var_1122, js_var_1123, lay_var_1124, js1_var_1125, jpl_var_1126
  REAL(KIND = 8) :: refrat_planck_a_var_1127, refrat_planck_b_var_1128
  REAL(KIND = 8) :: fac000_var_1129, fac100_var_1130, fac200_var_1131, fac010_var_1132, fac110_var_1133, fac210_var_1134, fac001_var_1135, fac101_var_1136, fac201_var_1137, fac011_var_1138, fac111_var_1139, fac211_var_1140
  REAL(KIND = 8) :: p_var_1141, p4_var_1142, fk0_var_1143, fk1_var_1144, fk2_var_1145
  REAL(KIND = 8) :: taufor_var_1146, tauself_var_1147, tau_major_var_1148(14), tau_major1_var_1149(14)
  REAL(KIND = 8) :: fs_var_1150, specmult_var_1151, specparm_var_1152, fs1_var_1153, specmult1_var_1154, specparm1_var_1155, fpl_var_1156, specmult_planck_var_1157, specparm_planck_var_1158
  INTEGER(KIND = 4) :: laytrop_min_var_1159, laytrop_max_var_1160
  INTEGER(KIND = 4) :: ixc_var_1161(klev_var_1089), ixlow_var_1162(kfdia_var_1088, klev_var_1089), ixhigh_var_1163(kfdia_var_1088, klev_var_1089)
  INTEGER(KIND = 4) :: ich_var_1164, icl_var_1165, ixc0_var_1166, ixp_var_1167, jc_var_1168, jl_var_1169
  laytrop_min_var_1159 = MINVAL(laytrop_var_1103)
  laytrop_max_var_1160 = MAXVAL(laytrop_var_1103)
  ixlow_var_1162 = 0
  ixhigh_var_1163 = 0
  ixc_var_1161 = 0
  DO lay_var_1124 = laytrop_min_var_1159 + 1, laytrop_max_var_1160
    icl_var_1165 = 0
    ich_var_1164 = 0
    DO jc_var_1168 = kidia_var_1087, kfdia_var_1088
      IF (lay_var_1124 <= laytrop_var_1103(jc_var_1168)) THEN
        icl_var_1165 = icl_var_1165 + 1
        ixlow_var_1162(icl_var_1165, lay_var_1124) = jc_var_1168
      ELSE
        ich_var_1164 = ich_var_1164 + 1
        ixhigh_var_1163(ich_var_1164, lay_var_1124) = jc_var_1168
      END IF
    END DO
    ixc_var_1161(lay_var_1124) = icl_var_1165
  END DO
  refrat_planck_a_var_1127 = chi_mls(1, 11) / chi_mls(2, 11)
  refrat_planck_b_var_1128 = chi_mls(3, 13) / chi_mls(2, 13)
  DO lay_var_1124 = 1, laytrop_min_var_1159
    DO jl_var_1169 = kidia_var_1087, kfdia_var_1088
      speccomb_var_1115 = colh2o_var_1100(jl_var_1169, lay_var_1124) + rat_h2oco2_var_1108(jl_var_1169, lay_var_1124) * colco2_var_1101(jl_var_1169, lay_var_1124)
      specparm_var_1152 = MIN(colh2o_var_1100(jl_var_1169, lay_var_1124) / speccomb_var_1115, 0.999999D0)
      specmult_var_1151 = 8.0D0 * (specparm_var_1152)
      js_var_1123 = 1 + INT(specmult_var_1151)
      fs_var_1150 = ((specmult_var_1151) - AINT((specmult_var_1151)))
      speccomb1_var_1116 = colh2o_var_1100(jl_var_1169, lay_var_1124) + rat_h2oco2_1_var_1109(jl_var_1169, lay_var_1124) * colco2_var_1101(jl_var_1169, lay_var_1124)
      specparm1_var_1155 = MIN(colh2o_var_1100(jl_var_1169, lay_var_1124) / speccomb1_var_1116, 0.999999D0)
      specmult1_var_1154 = 8.0D0 * (specparm1_var_1155)
      js1_var_1125 = 1 + INT(specmult1_var_1154)
      fs1_var_1153 = ((specmult1_var_1154) - AINT((specmult1_var_1154)))
      speccomb_planck_var_1117 = colh2o_var_1100(jl_var_1169, lay_var_1124) + refrat_planck_a_var_1127 * colco2_var_1101(jl_var_1169, lay_var_1124)
      specparm_planck_var_1158 = MIN(colh2o_var_1100(jl_var_1169, lay_var_1124) / speccomb_planck_var_1117, 0.999999D0)
      specmult_planck_var_1157 = 8.0D0 * specparm_planck_var_1158
      jpl_var_1126 = 1 + INT(specmult_planck_var_1157)
      fpl_var_1156 = ((specmult_planck_var_1157) - AINT((specmult_planck_var_1157)))
      ind0_var_1118 = ((jp_var_1096(jl_var_1169, lay_var_1124) - 1) * 5 + (jt_var_1097(jl_var_1169, lay_var_1124) - 1)) * nspa_var_216(4) + js_var_1123
      ind1_var_1119 = (jp_var_1096(jl_var_1169, lay_var_1124) * 5 + (jt1_var_1098(jl_var_1169, lay_var_1124) - 1)) * nspa_var_216(4) + js1_var_1125
      inds_var_1120 = indself_var_1106(jl_var_1169, lay_var_1124)
      indf_var_1121 = indfor_var_1112(jl_var_1169, lay_var_1124)
      IF (specparm_var_1152 .LT. 0.125D0) THEN
        p_var_1141 = fs_var_1150 - 1.0D0
        p4_var_1142 = p_var_1141 ** 4
        fk0_var_1143 = p4_var_1142
        fk1_var_1144 = 1.0D0 - p_var_1141 - 2.0D0 * p4_var_1142
        fk2_var_1145 = p_var_1141 + p4_var_1142
        fac000_var_1129 = fk0_var_1143 * fac00_var_1092(jl_var_1169, lay_var_1124)
        fac100_var_1130 = fk1_var_1144 * fac00_var_1092(jl_var_1169, lay_var_1124)
        fac200_var_1131 = fk2_var_1145 * fac00_var_1092(jl_var_1169, lay_var_1124)
        fac010_var_1132 = fk0_var_1143 * fac10_var_1094(jl_var_1169, lay_var_1124)
        fac110_var_1133 = fk1_var_1144 * fac10_var_1094(jl_var_1169, lay_var_1124)
        fac210_var_1134 = fk2_var_1145 * fac10_var_1094(jl_var_1169, lay_var_1124)
      ELSE IF (specparm_var_1152 .GT. 0.875D0) THEN
        p_var_1141 = - fs_var_1150
        p4_var_1142 = p_var_1141 ** 4
        fk0_var_1143 = p4_var_1142
        fk1_var_1144 = 1.0D0 - p_var_1141 - 2.0D0 * p4_var_1142
        fk2_var_1145 = p_var_1141 + p4_var_1142
        fac000_var_1129 = fk0_var_1143 * fac00_var_1092(jl_var_1169, lay_var_1124)
        fac100_var_1130 = fk1_var_1144 * fac00_var_1092(jl_var_1169, lay_var_1124)
        fac200_var_1131 = fk2_var_1145 * fac00_var_1092(jl_var_1169, lay_var_1124)
        fac010_var_1132 = fk0_var_1143 * fac10_var_1094(jl_var_1169, lay_var_1124)
        fac110_var_1133 = fk1_var_1144 * fac10_var_1094(jl_var_1169, lay_var_1124)
        fac210_var_1134 = fk2_var_1145 * fac10_var_1094(jl_var_1169, lay_var_1124)
      ELSE
        fac000_var_1129 = (1.0D0 - fs_var_1150) * fac00_var_1092(jl_var_1169, lay_var_1124)
        fac010_var_1132 = (1.0D0 - fs_var_1150) * fac10_var_1094(jl_var_1169, lay_var_1124)
        fac100_var_1130 = fs_var_1150 * fac00_var_1092(jl_var_1169, lay_var_1124)
        fac110_var_1133 = fs_var_1150 * fac10_var_1094(jl_var_1169, lay_var_1124)
        fac200_var_1131 = 0.0D0
        fac210_var_1134 = 0.0D0
      END IF
      IF (specparm1_var_1155 .LT. 0.125D0) THEN
        p_var_1141 = fs1_var_1153 - 1.0D0
        p4_var_1142 = p_var_1141 ** 4
        fk0_var_1143 = p4_var_1142
        fk1_var_1144 = 1.0D0 - p_var_1141 - 2.0D0 * p4_var_1142
        fk2_var_1145 = p_var_1141 + p4_var_1142
        fac001_var_1135 = fk0_var_1143 * fac01_var_1093(jl_var_1169, lay_var_1124)
        fac101_var_1136 = fk1_var_1144 * fac01_var_1093(jl_var_1169, lay_var_1124)
        fac201_var_1137 = fk2_var_1145 * fac01_var_1093(jl_var_1169, lay_var_1124)
        fac011_var_1138 = fk0_var_1143 * fac11_var_1095(jl_var_1169, lay_var_1124)
        fac111_var_1139 = fk1_var_1144 * fac11_var_1095(jl_var_1169, lay_var_1124)
        fac211_var_1140 = fk2_var_1145 * fac11_var_1095(jl_var_1169, lay_var_1124)
      ELSE IF (specparm1_var_1155 .GT. 0.875D0) THEN
        p_var_1141 = - fs1_var_1153
        p4_var_1142 = p_var_1141 ** 4
        fk0_var_1143 = p4_var_1142
        fk1_var_1144 = 1.0D0 - p_var_1141 - 2.0D0 * p4_var_1142
        fk2_var_1145 = p_var_1141 + p4_var_1142
        fac001_var_1135 = fk0_var_1143 * fac01_var_1093(jl_var_1169, lay_var_1124)
        fac101_var_1136 = fk1_var_1144 * fac01_var_1093(jl_var_1169, lay_var_1124)
        fac201_var_1137 = fk2_var_1145 * fac01_var_1093(jl_var_1169, lay_var_1124)
        fac011_var_1138 = fk0_var_1143 * fac11_var_1095(jl_var_1169, lay_var_1124)
        fac111_var_1139 = fk1_var_1144 * fac11_var_1095(jl_var_1169, lay_var_1124)
        fac211_var_1140 = fk2_var_1145 * fac11_var_1095(jl_var_1169, lay_var_1124)
      ELSE
        fac001_var_1135 = (1.0D0 - fs1_var_1153) * fac01_var_1093(jl_var_1169, lay_var_1124)
        fac011_var_1138 = (1.0D0 - fs1_var_1153) * fac11_var_1095(jl_var_1169, lay_var_1124)
        fac101_var_1136 = fs1_var_1153 * fac01_var_1093(jl_var_1169, lay_var_1124)
        fac111_var_1139 = fs1_var_1153 * fac11_var_1095(jl_var_1169, lay_var_1124)
        fac201_var_1137 = 0.0D0
        fac211_var_1140 = 0.0D0
      END IF
      IF (specparm_var_1152 .LT. 0.125D0) THEN
        tau_major_var_1148(1 : ng4) = speccomb_var_1115 * (fac000_var_1129 * absa_var_169(ind0_var_1118, 1 : 14) + fac100_var_1130 * absa_var_169(ind0_var_1118 + 1, 1 : 14) + fac200_var_1131 * absa_var_169(ind0_var_1118 + 2, 1 : 14) + fac010_var_1132 * absa_var_169(ind0_var_1118 + 9, 1 : 14) + fac110_var_1133 * absa_var_169(ind0_var_1118 + 10, 1 : 14) + fac210_var_1134 * absa_var_169(ind0_var_1118 + 11, 1 : 14))
      ELSE IF (specparm_var_1152 .GT. 0.875D0) THEN
        tau_major_var_1148(1 : ng4) = speccomb_var_1115 * (fac200_var_1131 * absa_var_169(ind0_var_1118 - 1, 1 : 14) + fac100_var_1130 * absa_var_169(ind0_var_1118, 1 : 14) + fac000_var_1129 * absa_var_169(ind0_var_1118 + 1, 1 : 14) + fac210_var_1134 * absa_var_169(ind0_var_1118 + 8, 1 : 14) + fac110_var_1133 * absa_var_169(ind0_var_1118 + 9, 1 : 14) + fac010_var_1132 * absa_var_169(ind0_var_1118 + 10, 1 : 14))
      ELSE
        tau_major_var_1148(1 : ng4) = speccomb_var_1115 * (fac000_var_1129 * absa_var_169(ind0_var_1118, 1 : 14) + fac100_var_1130 * absa_var_169(ind0_var_1118 + 1, 1 : 14) + fac010_var_1132 * absa_var_169(ind0_var_1118 + 9, 1 : 14) + fac110_var_1133 * absa_var_169(ind0_var_1118 + 10, 1 : 14))
      END IF
      IF (specparm1_var_1155 .LT. 0.125D0) THEN
        tau_major1_var_1149(1 : ng4) = speccomb1_var_1116 * (fac001_var_1135 * absa_var_169(ind1_var_1119, 1 : 14) + fac101_var_1136 * absa_var_169(ind1_var_1119 + 1, 1 : 14) + fac201_var_1137 * absa_var_169(ind1_var_1119 + 2, 1 : 14) + fac011_var_1138 * absa_var_169(ind1_var_1119 + 9, 1 : 14) + fac111_var_1139 * absa_var_169(ind1_var_1119 + 10, 1 : 14) + fac211_var_1140 * absa_var_169(ind1_var_1119 + 11, 1 : 14))
      ELSE IF (specparm1_var_1155 .GT. 0.875D0) THEN
        tau_major1_var_1149(1 : ng4) = speccomb1_var_1116 * (fac201_var_1137 * absa_var_169(ind1_var_1119 - 1, 1 : 14) + fac101_var_1136 * absa_var_169(ind1_var_1119, 1 : 14) + fac001_var_1135 * absa_var_169(ind1_var_1119 + 1, 1 : 14) + fac211_var_1140 * absa_var_169(ind1_var_1119 + 8, 1 : 14) + fac111_var_1139 * absa_var_169(ind1_var_1119 + 9, 1 : 14) + fac011_var_1138 * absa_var_169(ind1_var_1119 + 10, 1 : 14))
      ELSE
        tau_major1_var_1149(1 : ng4) = speccomb1_var_1116 * (fac001_var_1135 * absa_var_169(ind1_var_1119, 1 : 14) + fac101_var_1136 * absa_var_169(ind1_var_1119 + 1, 1 : 14) + fac011_var_1138 * absa_var_169(ind1_var_1119 + 9, 1 : 14) + fac111_var_1139 * absa_var_169(ind1_var_1119 + 10, 1 : 14))
      END IF
      DO ig_var_1122 = 1, 14
        tauself_var_1147 = selffac_var_1104(jl_var_1169, lay_var_1124) * (selfref_var_171(inds_var_1120, ig_var_1122) + selffrac_var_1105(jl_var_1169, lay_var_1124) * (selfref_var_171(inds_var_1120 + 1, ig_var_1122) - selfref_var_171(inds_var_1120, ig_var_1122)))
        taufor_var_1146 = forfac_var_1113(jl_var_1169, lay_var_1124) * (forref_var_172(indf_var_1121, ig_var_1122) + forfrac_var_1114(jl_var_1169, lay_var_1124) * (forref_var_172(indf_var_1121 + 1, ig_var_1122) - forref_var_172(indf_var_1121, ig_var_1122)))
        taug_var_1090(jl_var_1169, 38 + ig_var_1122, lay_var_1124) = tau_major_var_1148(ig_var_1122) + tau_major1_var_1149(ig_var_1122) + tauself_var_1147 + taufor_var_1146
        fracs_var_1107(jl_var_1169, 38 + ig_var_1122, lay_var_1124) = fracrefa_var_167(ig_var_1122, jpl_var_1126) + fpl_var_1156 * (fracrefa_var_167(ig_var_1122, jpl_var_1126 + 1) - fracrefa_var_167(ig_var_1122, jpl_var_1126))
      END DO
    END DO
  END DO
  DO lay_var_1124 = laytrop_max_var_1160 + 1, klev_var_1089
    DO jl_var_1169 = kidia_var_1087, kfdia_var_1088
      speccomb_var_1115 = colo3_var_1102(jl_var_1169, lay_var_1124) + rat_o3co2_var_1110(jl_var_1169, lay_var_1124) * colco2_var_1101(jl_var_1169, lay_var_1124)
      specparm_var_1152 = MIN(colo3_var_1102(jl_var_1169, lay_var_1124) / speccomb_var_1115, 0.999999D0)
      specmult_var_1151 = 4.0D0 * (specparm_var_1152)
      js_var_1123 = 1 + INT(specmult_var_1151)
      fs_var_1150 = ((specmult_var_1151) - AINT((specmult_var_1151)))
      speccomb1_var_1116 = colo3_var_1102(jl_var_1169, lay_var_1124) + rat_o3co2_1_var_1111(jl_var_1169, lay_var_1124) * colco2_var_1101(jl_var_1169, lay_var_1124)
      specparm1_var_1155 = MIN(colo3_var_1102(jl_var_1169, lay_var_1124) / speccomb1_var_1116, 0.999999D0)
      specmult1_var_1154 = 4.0D0 * (specparm1_var_1155)
      js1_var_1125 = 1 + INT(specmult1_var_1154)
      fs1_var_1153 = ((specmult1_var_1154) - AINT((specmult1_var_1154)))
      fac000_var_1129 = (1.0D0 - fs_var_1150) * fac00_var_1092(jl_var_1169, lay_var_1124)
      fac010_var_1132 = (1.0D0 - fs_var_1150) * fac10_var_1094(jl_var_1169, lay_var_1124)
      fac100_var_1130 = fs_var_1150 * fac00_var_1092(jl_var_1169, lay_var_1124)
      fac110_var_1133 = fs_var_1150 * fac10_var_1094(jl_var_1169, lay_var_1124)
      fac001_var_1135 = (1.0D0 - fs1_var_1153) * fac01_var_1093(jl_var_1169, lay_var_1124)
      fac011_var_1138 = (1.0D0 - fs1_var_1153) * fac11_var_1095(jl_var_1169, lay_var_1124)
      fac101_var_1136 = fs1_var_1153 * fac01_var_1093(jl_var_1169, lay_var_1124)
      fac111_var_1139 = fs1_var_1153 * fac11_var_1095(jl_var_1169, lay_var_1124)
      speccomb_planck_var_1117 = colo3_var_1102(jl_var_1169, lay_var_1124) + refrat_planck_b_var_1128 * colco2_var_1101(jl_var_1169, lay_var_1124)
      specparm_planck_var_1158 = MIN(colo3_var_1102(jl_var_1169, lay_var_1124) / speccomb_planck_var_1117, 0.999999D0)
      specmult_planck_var_1157 = 4.0D0 * specparm_planck_var_1158
      jpl_var_1126 = 1 + INT(specmult_planck_var_1157)
      fpl_var_1156 = ((specmult_planck_var_1157) - AINT((specmult_planck_var_1157)))
      ind0_var_1118 = ((jp_var_1096(jl_var_1169, lay_var_1124) - 13) * 5 + (jt_var_1097(jl_var_1169, lay_var_1124) - 1)) * nspb_var_217(4) + js_var_1123
      ind1_var_1119 = ((jp_var_1096(jl_var_1169, lay_var_1124) - 12) * 5 + (jt1_var_1098(jl_var_1169, lay_var_1124) - 1)) * nspb_var_217(4) + js1_var_1125
      DO ig_var_1122 = 1, 14
        taug_var_1090(jl_var_1169, 38 + ig_var_1122, lay_var_1124) = speccomb_var_1115 * (fac000_var_1129 * absb_var_170(ind0_var_1118, ig_var_1122) + fac100_var_1130 * absb_var_170(ind0_var_1118 + 1, ig_var_1122) + fac010_var_1132 * absb_var_170(ind0_var_1118 + 5, ig_var_1122) + fac110_var_1133 * absb_var_170(ind0_var_1118 + 6, ig_var_1122)) + speccomb1_var_1116 * (fac001_var_1135 * absb_var_170(ind1_var_1119, ig_var_1122) + fac101_var_1136 * absb_var_170(ind1_var_1119 + 1, ig_var_1122) + fac011_var_1138 * absb_var_170(ind1_var_1119 + 5, ig_var_1122) + fac111_var_1139 * absb_var_170(ind1_var_1119 + 6, ig_var_1122))
        fracs_var_1107(jl_var_1169, 38 + ig_var_1122, lay_var_1124) = fracrefb_var_168(ig_var_1122, jpl_var_1126) + fpl_var_1156 * (fracrefb_var_168(ig_var_1122, jpl_var_1126 + 1) - fracrefb_var_168(ig_var_1122, jpl_var_1126))
      END DO
    END DO
  END DO
  DO lay_var_1124 = laytrop_max_var_1160 + 1, klev_var_1089
    DO jl_var_1169 = kidia_var_1087, kfdia_var_1088
      taug_var_1090(jl_var_1169, 46, lay_var_1124) = taug_var_1090(jl_var_1169, 46, lay_var_1124) * 0.92D0
      taug_var_1090(jl_var_1169, 47, lay_var_1124) = taug_var_1090(jl_var_1169, 47, lay_var_1124) * 0.88D0
      taug_var_1090(jl_var_1169, 48, lay_var_1124) = taug_var_1090(jl_var_1169, 48, lay_var_1124) * 1.07D0
      taug_var_1090(jl_var_1169, 49, lay_var_1124) = taug_var_1090(jl_var_1169, 49, lay_var_1124) * 1.1D0
      taug_var_1090(jl_var_1169, 50, lay_var_1124) = taug_var_1090(jl_var_1169, 50, lay_var_1124) * 0.99D0
      taug_var_1090(jl_var_1169, 51, lay_var_1124) = taug_var_1090(jl_var_1169, 51, lay_var_1124) * 0.88D0
      taug_var_1090(jl_var_1169, 52, lay_var_1124) = taug_var_1090(jl_var_1169, 52, lay_var_1124) * 0.943D0
    END DO
  END DO
  IF (laytrop_max_var_1160 /= laytrop_min_var_1159) THEN
    DO lay_var_1124 = laytrop_min_var_1159 + 1, laytrop_max_var_1160
      ixc0_var_1166 = ixc_var_1161(lay_var_1124)
      DO ixp_var_1167 = 1, ixc0_var_1166
        jl_var_1169 = ixlow_var_1162(ixp_var_1167, lay_var_1124)
        speccomb_var_1115 = colh2o_var_1100(jl_var_1169, lay_var_1124) + rat_h2oco2_var_1108(jl_var_1169, lay_var_1124) * colco2_var_1101(jl_var_1169, lay_var_1124)
        specparm_var_1152 = MIN(colh2o_var_1100(jl_var_1169, lay_var_1124) / speccomb_var_1115, 0.999999D0)
        specmult_var_1151 = 8.0D0 * (specparm_var_1152)
        js_var_1123 = 1 + INT(specmult_var_1151)
        fs_var_1150 = ((specmult_var_1151) - AINT((specmult_var_1151)))
        speccomb1_var_1116 = colh2o_var_1100(jl_var_1169, lay_var_1124) + rat_h2oco2_1_var_1109(jl_var_1169, lay_var_1124) * colco2_var_1101(jl_var_1169, lay_var_1124)
        specparm1_var_1155 = MIN(colh2o_var_1100(jl_var_1169, lay_var_1124) / speccomb1_var_1116, 0.999999D0)
        specmult1_var_1154 = 8.0D0 * (specparm1_var_1155)
        js1_var_1125 = 1 + INT(specmult1_var_1154)
        fs1_var_1153 = ((specmult1_var_1154) - AINT((specmult1_var_1154)))
        speccomb_planck_var_1117 = colh2o_var_1100(jl_var_1169, lay_var_1124) + refrat_planck_a_var_1127 * colco2_var_1101(jl_var_1169, lay_var_1124)
        specparm_planck_var_1158 = MIN(colh2o_var_1100(jl_var_1169, lay_var_1124) / speccomb_planck_var_1117, 0.999999D0)
        specmult_planck_var_1157 = 8.0D0 * specparm_planck_var_1158
        jpl_var_1126 = 1 + INT(specmult_planck_var_1157)
        fpl_var_1156 = ((specmult_planck_var_1157) - AINT((specmult_planck_var_1157)))
        ind0_var_1118 = ((jp_var_1096(jl_var_1169, lay_var_1124) - 1) * 5 + (jt_var_1097(jl_var_1169, lay_var_1124) - 1)) * nspa_var_216(4) + js_var_1123
        ind1_var_1119 = (jp_var_1096(jl_var_1169, lay_var_1124) * 5 + (jt1_var_1098(jl_var_1169, lay_var_1124) - 1)) * nspa_var_216(4) + js1_var_1125
        inds_var_1120 = indself_var_1106(jl_var_1169, lay_var_1124)
        indf_var_1121 = indfor_var_1112(jl_var_1169, lay_var_1124)
        IF (specparm_var_1152 .LT. 0.125D0) THEN
          p_var_1141 = fs_var_1150 - 1.0D0
          p4_var_1142 = p_var_1141 ** 4
          fk0_var_1143 = p4_var_1142
          fk1_var_1144 = 1.0D0 - p_var_1141 - 2.0D0 * p4_var_1142
          fk2_var_1145 = p_var_1141 + p4_var_1142
          fac000_var_1129 = fk0_var_1143 * fac00_var_1092(jl_var_1169, lay_var_1124)
          fac100_var_1130 = fk1_var_1144 * fac00_var_1092(jl_var_1169, lay_var_1124)
          fac200_var_1131 = fk2_var_1145 * fac00_var_1092(jl_var_1169, lay_var_1124)
          fac010_var_1132 = fk0_var_1143 * fac10_var_1094(jl_var_1169, lay_var_1124)
          fac110_var_1133 = fk1_var_1144 * fac10_var_1094(jl_var_1169, lay_var_1124)
          fac210_var_1134 = fk2_var_1145 * fac10_var_1094(jl_var_1169, lay_var_1124)
        ELSE IF (specparm_var_1152 .GT. 0.875D0) THEN
          p_var_1141 = - fs_var_1150
          p4_var_1142 = p_var_1141 ** 4
          fk0_var_1143 = p4_var_1142
          fk1_var_1144 = 1.0D0 - p_var_1141 - 2.0D0 * p4_var_1142
          fk2_var_1145 = p_var_1141 + p4_var_1142
          fac000_var_1129 = fk0_var_1143 * fac00_var_1092(jl_var_1169, lay_var_1124)
          fac100_var_1130 = fk1_var_1144 * fac00_var_1092(jl_var_1169, lay_var_1124)
          fac200_var_1131 = fk2_var_1145 * fac00_var_1092(jl_var_1169, lay_var_1124)
          fac010_var_1132 = fk0_var_1143 * fac10_var_1094(jl_var_1169, lay_var_1124)
          fac110_var_1133 = fk1_var_1144 * fac10_var_1094(jl_var_1169, lay_var_1124)
          fac210_var_1134 = fk2_var_1145 * fac10_var_1094(jl_var_1169, lay_var_1124)
        ELSE
          fac000_var_1129 = (1.0D0 - fs_var_1150) * fac00_var_1092(jl_var_1169, lay_var_1124)
          fac010_var_1132 = (1.0D0 - fs_var_1150) * fac10_var_1094(jl_var_1169, lay_var_1124)
          fac100_var_1130 = fs_var_1150 * fac00_var_1092(jl_var_1169, lay_var_1124)
          fac110_var_1133 = fs_var_1150 * fac10_var_1094(jl_var_1169, lay_var_1124)
          fac200_var_1131 = 0.0D0
          fac210_var_1134 = 0.0D0
        END IF
        IF (specparm1_var_1155 .LT. 0.125D0) THEN
          p_var_1141 = fs1_var_1153 - 1.0D0
          p4_var_1142 = p_var_1141 ** 4
          fk0_var_1143 = p4_var_1142
          fk1_var_1144 = 1.0D0 - p_var_1141 - 2.0D0 * p4_var_1142
          fk2_var_1145 = p_var_1141 + p4_var_1142
          fac001_var_1135 = fk0_var_1143 * fac01_var_1093(jl_var_1169, lay_var_1124)
          fac101_var_1136 = fk1_var_1144 * fac01_var_1093(jl_var_1169, lay_var_1124)
          fac201_var_1137 = fk2_var_1145 * fac01_var_1093(jl_var_1169, lay_var_1124)
          fac011_var_1138 = fk0_var_1143 * fac11_var_1095(jl_var_1169, lay_var_1124)
          fac111_var_1139 = fk1_var_1144 * fac11_var_1095(jl_var_1169, lay_var_1124)
          fac211_var_1140 = fk2_var_1145 * fac11_var_1095(jl_var_1169, lay_var_1124)
        ELSE IF (specparm1_var_1155 .GT. 0.875D0) THEN
          p_var_1141 = - fs1_var_1153
          p4_var_1142 = p_var_1141 ** 4
          fk0_var_1143 = p4_var_1142
          fk1_var_1144 = 1.0D0 - p_var_1141 - 2.0D0 * p4_var_1142
          fk2_var_1145 = p_var_1141 + p4_var_1142
          fac001_var_1135 = fk0_var_1143 * fac01_var_1093(jl_var_1169, lay_var_1124)
          fac101_var_1136 = fk1_var_1144 * fac01_var_1093(jl_var_1169, lay_var_1124)
          fac201_var_1137 = fk2_var_1145 * fac01_var_1093(jl_var_1169, lay_var_1124)
          fac011_var_1138 = fk0_var_1143 * fac11_var_1095(jl_var_1169, lay_var_1124)
          fac111_var_1139 = fk1_var_1144 * fac11_var_1095(jl_var_1169, lay_var_1124)
          fac211_var_1140 = fk2_var_1145 * fac11_var_1095(jl_var_1169, lay_var_1124)
        ELSE
          fac001_var_1135 = (1.0D0 - fs1_var_1153) * fac01_var_1093(jl_var_1169, lay_var_1124)
          fac011_var_1138 = (1.0D0 - fs1_var_1153) * fac11_var_1095(jl_var_1169, lay_var_1124)
          fac101_var_1136 = fs1_var_1153 * fac01_var_1093(jl_var_1169, lay_var_1124)
          fac111_var_1139 = fs1_var_1153 * fac11_var_1095(jl_var_1169, lay_var_1124)
          fac201_var_1137 = 0.0D0
          fac211_var_1140 = 0.0D0
        END IF
        IF (specparm_var_1152 .LT. 0.125D0) THEN
          tau_major_var_1148(1 : ng4) = speccomb_var_1115 * (fac000_var_1129 * absa_var_169(ind0_var_1118, 1 : 14) + fac100_var_1130 * absa_var_169(ind0_var_1118 + 1, 1 : 14) + fac200_var_1131 * absa_var_169(ind0_var_1118 + 2, 1 : 14) + fac010_var_1132 * absa_var_169(ind0_var_1118 + 9, 1 : 14) + fac110_var_1133 * absa_var_169(ind0_var_1118 + 10, 1 : 14) + fac210_var_1134 * absa_var_169(ind0_var_1118 + 11, 1 : 14))
        ELSE IF (specparm_var_1152 .GT. 0.875D0) THEN
          tau_major_var_1148(1 : ng4) = speccomb_var_1115 * (fac200_var_1131 * absa_var_169(ind0_var_1118 - 1, 1 : 14) + fac100_var_1130 * absa_var_169(ind0_var_1118, 1 : 14) + fac000_var_1129 * absa_var_169(ind0_var_1118 + 1, 1 : 14) + fac210_var_1134 * absa_var_169(ind0_var_1118 + 8, 1 : 14) + fac110_var_1133 * absa_var_169(ind0_var_1118 + 9, 1 : 14) + fac010_var_1132 * absa_var_169(ind0_var_1118 + 10, 1 : 14))
        ELSE
          tau_major_var_1148(1 : ng4) = speccomb_var_1115 * (fac000_var_1129 * absa_var_169(ind0_var_1118, 1 : 14) + fac100_var_1130 * absa_var_169(ind0_var_1118 + 1, 1 : 14) + fac010_var_1132 * absa_var_169(ind0_var_1118 + 9, 1 : 14) + fac110_var_1133 * absa_var_169(ind0_var_1118 + 10, 1 : 14))
        END IF
        IF (specparm1_var_1155 .LT. 0.125D0) THEN
          tau_major1_var_1149(1 : ng4) = speccomb1_var_1116 * (fac001_var_1135 * absa_var_169(ind1_var_1119, 1 : 14) + fac101_var_1136 * absa_var_169(ind1_var_1119 + 1, 1 : 14) + fac201_var_1137 * absa_var_169(ind1_var_1119 + 2, 1 : 14) + fac011_var_1138 * absa_var_169(ind1_var_1119 + 9, 1 : 14) + fac111_var_1139 * absa_var_169(ind1_var_1119 + 10, 1 : 14) + fac211_var_1140 * absa_var_169(ind1_var_1119 + 11, 1 : 14))
        ELSE IF (specparm1_var_1155 .GT. 0.875D0) THEN
          tau_major1_var_1149(1 : ng4) = speccomb1_var_1116 * (fac201_var_1137 * absa_var_169(ind1_var_1119 - 1, 1 : 14) + fac101_var_1136 * absa_var_169(ind1_var_1119, 1 : 14) + fac001_var_1135 * absa_var_169(ind1_var_1119 + 1, 1 : 14) + fac211_var_1140 * absa_var_169(ind1_var_1119 + 8, 1 : 14) + fac111_var_1139 * absa_var_169(ind1_var_1119 + 9, 1 : 14) + fac011_var_1138 * absa_var_169(ind1_var_1119 + 10, 1 : 14))
        ELSE
          tau_major1_var_1149(1 : ng4) = speccomb1_var_1116 * (fac001_var_1135 * absa_var_169(ind1_var_1119, 1 : 14) + fac101_var_1136 * absa_var_169(ind1_var_1119 + 1, 1 : 14) + fac011_var_1138 * absa_var_169(ind1_var_1119 + 9, 1 : 14) + fac111_var_1139 * absa_var_169(ind1_var_1119 + 10, 1 : 14))
        END IF
        DO ig_var_1122 = 1, 14
          tauself_var_1147 = selffac_var_1104(jl_var_1169, lay_var_1124) * (selfref_var_171(inds_var_1120, ig_var_1122) + selffrac_var_1105(jl_var_1169, lay_var_1124) * (selfref_var_171(inds_var_1120 + 1, ig_var_1122) - selfref_var_171(inds_var_1120, ig_var_1122)))
          taufor_var_1146 = forfac_var_1113(jl_var_1169, lay_var_1124) * (forref_var_172(indf_var_1121, ig_var_1122) + forfrac_var_1114(jl_var_1169, lay_var_1124) * (forref_var_172(indf_var_1121 + 1, ig_var_1122) - forref_var_172(indf_var_1121, ig_var_1122)))
          taug_var_1090(jl_var_1169, 38 + ig_var_1122, lay_var_1124) = tau_major_var_1148(ig_var_1122) + tau_major1_var_1149(ig_var_1122) + tauself_var_1147 + taufor_var_1146
          fracs_var_1107(jl_var_1169, 38 + ig_var_1122, lay_var_1124) = fracrefa_var_167(ig_var_1122, jpl_var_1126) + fpl_var_1156 * (fracrefa_var_167(ig_var_1122, jpl_var_1126 + 1) - fracrefa_var_167(ig_var_1122, jpl_var_1126))
        END DO
      END DO
      ixc0_var_1166 = kfdia_var_1088 - kidia_var_1087 + 1 - ixc0_var_1166
      DO ixp_var_1167 = 1, ixc0_var_1166
        jl_var_1169 = ixhigh_var_1163(ixp_var_1167, lay_var_1124)
        speccomb_var_1115 = colo3_var_1102(jl_var_1169, lay_var_1124) + rat_o3co2_var_1110(jl_var_1169, lay_var_1124) * colco2_var_1101(jl_var_1169, lay_var_1124)
        specparm_var_1152 = MIN(colo3_var_1102(jl_var_1169, lay_var_1124) / speccomb_var_1115, 0.999999D0)
        specmult_var_1151 = 4.0D0 * (specparm_var_1152)
        js_var_1123 = 1 + INT(specmult_var_1151)
        fs_var_1150 = ((specmult_var_1151) - AINT((specmult_var_1151)))
        speccomb1_var_1116 = colo3_var_1102(jl_var_1169, lay_var_1124) + rat_o3co2_1_var_1111(jl_var_1169, lay_var_1124) * colco2_var_1101(jl_var_1169, lay_var_1124)
        specparm1_var_1155 = MIN(colo3_var_1102(jl_var_1169, lay_var_1124) / speccomb1_var_1116, 0.999999D0)
        specmult1_var_1154 = 4.0D0 * (specparm1_var_1155)
        js1_var_1125 = 1 + INT(specmult1_var_1154)
        fs1_var_1153 = ((specmult1_var_1154) - AINT((specmult1_var_1154)))
        fac000_var_1129 = (1.0D0 - fs_var_1150) * fac00_var_1092(jl_var_1169, lay_var_1124)
        fac010_var_1132 = (1.0D0 - fs_var_1150) * fac10_var_1094(jl_var_1169, lay_var_1124)
        fac100_var_1130 = fs_var_1150 * fac00_var_1092(jl_var_1169, lay_var_1124)
        fac110_var_1133 = fs_var_1150 * fac10_var_1094(jl_var_1169, lay_var_1124)
        fac001_var_1135 = (1.0D0 - fs1_var_1153) * fac01_var_1093(jl_var_1169, lay_var_1124)
        fac011_var_1138 = (1.0D0 - fs1_var_1153) * fac11_var_1095(jl_var_1169, lay_var_1124)
        fac101_var_1136 = fs1_var_1153 * fac01_var_1093(jl_var_1169, lay_var_1124)
        fac111_var_1139 = fs1_var_1153 * fac11_var_1095(jl_var_1169, lay_var_1124)
        speccomb_planck_var_1117 = colo3_var_1102(jl_var_1169, lay_var_1124) + refrat_planck_b_var_1128 * colco2_var_1101(jl_var_1169, lay_var_1124)
        specparm_planck_var_1158 = MIN(colo3_var_1102(jl_var_1169, lay_var_1124) / speccomb_planck_var_1117, 0.999999D0)
        specmult_planck_var_1157 = 4.0D0 * specparm_planck_var_1158
        jpl_var_1126 = 1 + INT(specmult_planck_var_1157)
        fpl_var_1156 = ((specmult_planck_var_1157) - AINT((specmult_planck_var_1157)))
        ind0_var_1118 = ((jp_var_1096(jl_var_1169, lay_var_1124) - 13) * 5 + (jt_var_1097(jl_var_1169, lay_var_1124) - 1)) * nspb_var_217(4) + js_var_1123
        ind1_var_1119 = ((jp_var_1096(jl_var_1169, lay_var_1124) - 12) * 5 + (jt1_var_1098(jl_var_1169, lay_var_1124) - 1)) * nspb_var_217(4) + js1_var_1125
        DO ig_var_1122 = 1, 14
          taug_var_1090(jl_var_1169, 38 + ig_var_1122, lay_var_1124) = speccomb_var_1115 * (fac000_var_1129 * absb_var_170(ind0_var_1118, ig_var_1122) + fac100_var_1130 * absb_var_170(ind0_var_1118 + 1, ig_var_1122) + fac010_var_1132 * absb_var_170(ind0_var_1118 + 5, ig_var_1122) + fac110_var_1133 * absb_var_170(ind0_var_1118 + 6, ig_var_1122)) + speccomb1_var_1116 * (fac001_var_1135 * absb_var_170(ind1_var_1119, ig_var_1122) + fac101_var_1136 * absb_var_170(ind1_var_1119 + 1, ig_var_1122) + fac011_var_1138 * absb_var_170(ind1_var_1119 + 5, ig_var_1122) + fac111_var_1139 * absb_var_170(ind1_var_1119 + 6, ig_var_1122))
          fracs_var_1107(jl_var_1169, 38 + ig_var_1122, lay_var_1124) = fracrefb_var_168(ig_var_1122, jpl_var_1126) + fpl_var_1156 * (fracrefb_var_168(ig_var_1122, jpl_var_1126 + 1) - fracrefb_var_168(ig_var_1122, jpl_var_1126))
        END DO
      END DO
      DO ixp_var_1167 = 1, ixc0_var_1166
        jl_var_1169 = ixhigh_var_1163(ixp_var_1167, lay_var_1124)
        taug_var_1090(jl_var_1169, 46, lay_var_1124) = taug_var_1090(jl_var_1169, 46, lay_var_1124) * 0.92D0
        taug_var_1090(jl_var_1169, 47, lay_var_1124) = taug_var_1090(jl_var_1169, 47, lay_var_1124) * 0.88D0
        taug_var_1090(jl_var_1169, 48, lay_var_1124) = taug_var_1090(jl_var_1169, 48, lay_var_1124) * 1.07D0
        taug_var_1090(jl_var_1169, 49, lay_var_1124) = taug_var_1090(jl_var_1169, 49, lay_var_1124) * 1.1D0
        taug_var_1090(jl_var_1169, 50, lay_var_1124) = taug_var_1090(jl_var_1169, 50, lay_var_1124) * 0.99D0
        taug_var_1090(jl_var_1169, 51, lay_var_1124) = taug_var_1090(jl_var_1169, 51, lay_var_1124) * 0.88D0
        taug_var_1090(jl_var_1169, 52, lay_var_1124) = taug_var_1090(jl_var_1169, 52, lay_var_1124) * 0.943D0
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol4
SUBROUTINE rrtm_taumol6(kidia_var_1170, kfdia_var_1171, klev_var_1172, taug_var_1173, wx_var_1174, p_tauaerl_var_1175, fac00_var_1176, fac01_var_1177, fac10_var_1178, fac11_var_1179, forfac_var_1192, forfrac_var_1193, indfor_var_1191, jp_var_1180, jt_var_1181, jt1_var_1182, colh2o_var_1183, colco2_var_1184, coldry_var_1185, laytrop_var_1186, selffac_var_1187, selffrac_var_1188, indself_var_1189, fracs_var_1190, minorfrac_var_1194, indminor_var_1195)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrta6, ONLY: absa_var_182, cfc11adj, cfc12_var_181, forref_var_185, fracrefa_var_180, ka_mco2_var_184, selfref_var_183
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1170
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1171
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1172
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1173(kidia_var_1170 : kfdia_var_1171, 140, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: wx_var_1174(kidia_var_1170 : kfdia_var_1171, 4, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1175(kidia_var_1170 : kfdia_var_1171, klev_var_1172, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1176(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1177(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1178(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1179(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1180(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1181(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1182(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1183(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1184(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1185(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1186(kidia_var_1170 : kfdia_var_1171)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1187(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1188(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1189(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1190(kidia_var_1170 : kfdia_var_1171, 140, klev_var_1172)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1191(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1192(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1193(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1194(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1195(kidia_var_1170 : kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4) :: ind0_var_1196, ind1_var_1197, inds_var_1198, indf_var_1199, indm_var_1200
  INTEGER(KIND = 4) :: ig_var_1201, lay_var_1202
  REAL(KIND = 8) :: adjfac_var_1203, adjcolco2_var_1204, ratco2_var_1205, chi_co2_var_1206
  REAL(KIND = 8) :: taufor_var_1207, tauself_var_1208, absco2_var_1209
  INTEGER(KIND = 4) :: laytrop_min_var_1210, laytrop_max_var_1211
  INTEGER(KIND = 4) :: ixc_var_1212(klev_var_1172), ixlow_var_1213(kfdia_var_1171, klev_var_1172), ixhigh_var_1214(kfdia_var_1171, klev_var_1172)
  INTEGER(KIND = 4) :: ich_var_1215, icl_var_1216, ixc0_var_1217, ixp_var_1218, jc_var_1219, jl_var_1220
  laytrop_min_var_1210 = MINVAL(laytrop_var_1186)
  laytrop_max_var_1211 = MAXVAL(laytrop_var_1186)
  ixlow_var_1213 = 0
  ixhigh_var_1214 = 0
  ixc_var_1212 = 0
  DO lay_var_1202 = laytrop_min_var_1210 + 1, laytrop_max_var_1211
    icl_var_1216 = 0
    ich_var_1215 = 0
    DO jc_var_1219 = kidia_var_1170, kfdia_var_1171
      IF (lay_var_1202 <= laytrop_var_1186(jc_var_1219)) THEN
        icl_var_1216 = icl_var_1216 + 1
        ixlow_var_1213(icl_var_1216, lay_var_1202) = jc_var_1219
      ELSE
        ich_var_1215 = ich_var_1215 + 1
        ixhigh_var_1214(ich_var_1215, lay_var_1202) = jc_var_1219
      END IF
    END DO
    ixc_var_1212(lay_var_1202) = icl_var_1216
  END DO
  DO lay_var_1202 = 1, laytrop_min_var_1210
    DO jl_var_1220 = kidia_var_1170, kfdia_var_1171
      chi_co2_var_1206 = colco2_var_1184(jl_var_1220, lay_var_1202) / (coldry_var_1185(jl_var_1220, lay_var_1202))
      ratco2_var_1205 = 1D+20 * chi_co2_var_1206 / chi_mls(2, jp_var_1180(jl_var_1220, lay_var_1202) + 1)
      IF (ratco2_var_1205 .GT. 3.0D0) THEN
        adjfac_var_1203 = 2.0D0 + (ratco2_var_1205 - 2.0D0) ** 0.77D0
        adjcolco2_var_1204 = adjfac_var_1203 * chi_mls(2, jp_var_1180(jl_var_1220, lay_var_1202) + 1) * coldry_var_1185(jl_var_1220, lay_var_1202) * 1D-20
      ELSE
        adjcolco2_var_1204 = colco2_var_1184(jl_var_1220, lay_var_1202)
      END IF
      ind0_var_1196 = ((jp_var_1180(jl_var_1220, lay_var_1202) - 1) * 5 + (jt_var_1181(jl_var_1220, lay_var_1202) - 1)) * nspa_var_216(6) + 1
      ind1_var_1197 = (jp_var_1180(jl_var_1220, lay_var_1202) * 5 + (jt1_var_1182(jl_var_1220, lay_var_1202) - 1)) * nspa_var_216(6) + 1
      inds_var_1198 = indself_var_1189(jl_var_1220, lay_var_1202)
      indf_var_1199 = indfor_var_1191(jl_var_1220, lay_var_1202)
      indm_var_1200 = indminor_var_1195(jl_var_1220, lay_var_1202)
      DO ig_var_1201 = 1, 8
        tauself_var_1208 = selffac_var_1187(jl_var_1220, lay_var_1202) * (selfref_var_183(inds_var_1198, ig_var_1201) + selffrac_var_1188(jl_var_1220, lay_var_1202) * (selfref_var_183(inds_var_1198 + 1, ig_var_1201) - selfref_var_183(inds_var_1198, ig_var_1201)))
        taufor_var_1207 = forfac_var_1192(jl_var_1220, lay_var_1202) * (forref_var_185(indf_var_1199, ig_var_1201) + forfrac_var_1193(jl_var_1220, lay_var_1202) * (forref_var_185(indf_var_1199 + 1, ig_var_1201) - forref_var_185(indf_var_1199, ig_var_1201)))
        absco2_var_1209 = (ka_mco2_var_184(indm_var_1200, ig_var_1201) + minorfrac_var_1194(jl_var_1220, lay_var_1202) * (ka_mco2_var_184(indm_var_1200 + 1, ig_var_1201) - ka_mco2_var_184(indm_var_1200, ig_var_1201)))
        taug_var_1173(jl_var_1220, 68 + ig_var_1201, lay_var_1202) = colh2o_var_1183(jl_var_1220, lay_var_1202) * (fac00_var_1176(jl_var_1220, lay_var_1202) * absa_var_182(ind0_var_1196, ig_var_1201) + fac10_var_1178(jl_var_1220, lay_var_1202) * absa_var_182(ind0_var_1196 + 1, ig_var_1201) + fac01_var_1177(jl_var_1220, lay_var_1202) * absa_var_182(ind1_var_1197, ig_var_1201) + fac11_var_1179(jl_var_1220, lay_var_1202) * absa_var_182(ind1_var_1197 + 1, ig_var_1201)) + tauself_var_1208 + taufor_var_1207 + adjcolco2_var_1204 * absco2_var_1209 + wx_var_1174(jl_var_1220, 2, lay_var_1202) * cfc11adj(ig_var_1201) + wx_var_1174(jl_var_1220, 3, lay_var_1202) * cfc12_var_181(ig_var_1201)
        fracs_var_1190(jl_var_1220, 68 + ig_var_1201, lay_var_1202) = fracrefa_var_180(ig_var_1201)
      END DO
    END DO
  END DO
  DO ig_var_1201 = 1, 8
    DO lay_var_1202 = laytrop_max_var_1211 + 1, klev_var_1172
      DO jl_var_1220 = kidia_var_1170, kfdia_var_1171
        taug_var_1173(jl_var_1220, 68 + ig_var_1201, lay_var_1202) = 0.0D0 + wx_var_1174(jl_var_1220, 2, lay_var_1202) * cfc11adj(ig_var_1201) + wx_var_1174(jl_var_1220, 3, lay_var_1202) * cfc12_var_181(ig_var_1201)
        fracs_var_1190(jl_var_1220, 68 + ig_var_1201, lay_var_1202) = fracrefa_var_180(ig_var_1201)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1211 /= laytrop_min_var_1210) THEN
    DO lay_var_1202 = laytrop_min_var_1210 + 1, laytrop_max_var_1211
      ixc0_var_1217 = ixc_var_1212(lay_var_1202)
      DO ixp_var_1218 = 1, ixc0_var_1217
        jl_var_1220 = ixlow_var_1213(ixp_var_1218, lay_var_1202)
        chi_co2_var_1206 = colco2_var_1184(jl_var_1220, lay_var_1202) / (coldry_var_1185(jl_var_1220, lay_var_1202))
        ratco2_var_1205 = 1D+20 * chi_co2_var_1206 / chi_mls(2, jp_var_1180(jl_var_1220, lay_var_1202) + 1)
        IF (ratco2_var_1205 .GT. 3.0D0) THEN
          adjfac_var_1203 = 2.0D0 + (ratco2_var_1205 - 2.0D0) ** 0.77D0
          adjcolco2_var_1204 = adjfac_var_1203 * chi_mls(2, jp_var_1180(jl_var_1220, lay_var_1202) + 1) * coldry_var_1185(jl_var_1220, lay_var_1202) * 1D-20
        ELSE
          adjcolco2_var_1204 = colco2_var_1184(jl_var_1220, lay_var_1202)
        END IF
        ind0_var_1196 = ((jp_var_1180(jl_var_1220, lay_var_1202) - 1) * 5 + (jt_var_1181(jl_var_1220, lay_var_1202) - 1)) * nspa_var_216(6) + 1
        ind1_var_1197 = (jp_var_1180(jl_var_1220, lay_var_1202) * 5 + (jt1_var_1182(jl_var_1220, lay_var_1202) - 1)) * nspa_var_216(6) + 1
        inds_var_1198 = indself_var_1189(jl_var_1220, lay_var_1202)
        indf_var_1199 = indfor_var_1191(jl_var_1220, lay_var_1202)
        indm_var_1200 = indminor_var_1195(jl_var_1220, lay_var_1202)
        DO ig_var_1201 = 1, 8
          tauself_var_1208 = selffac_var_1187(jl_var_1220, lay_var_1202) * (selfref_var_183(inds_var_1198, ig_var_1201) + selffrac_var_1188(jl_var_1220, lay_var_1202) * (selfref_var_183(inds_var_1198 + 1, ig_var_1201) - selfref_var_183(inds_var_1198, ig_var_1201)))
          taufor_var_1207 = forfac_var_1192(jl_var_1220, lay_var_1202) * (forref_var_185(indf_var_1199, ig_var_1201) + forfrac_var_1193(jl_var_1220, lay_var_1202) * (forref_var_185(indf_var_1199 + 1, ig_var_1201) - forref_var_185(indf_var_1199, ig_var_1201)))
          absco2_var_1209 = (ka_mco2_var_184(indm_var_1200, ig_var_1201) + minorfrac_var_1194(jl_var_1220, lay_var_1202) * (ka_mco2_var_184(indm_var_1200 + 1, ig_var_1201) - ka_mco2_var_184(indm_var_1200, ig_var_1201)))
          taug_var_1173(jl_var_1220, 68 + ig_var_1201, lay_var_1202) = colh2o_var_1183(jl_var_1220, lay_var_1202) * (fac00_var_1176(jl_var_1220, lay_var_1202) * absa_var_182(ind0_var_1196, ig_var_1201) + fac10_var_1178(jl_var_1220, lay_var_1202) * absa_var_182(ind0_var_1196 + 1, ig_var_1201) + fac01_var_1177(jl_var_1220, lay_var_1202) * absa_var_182(ind1_var_1197, ig_var_1201) + fac11_var_1179(jl_var_1220, lay_var_1202) * absa_var_182(ind1_var_1197 + 1, ig_var_1201)) + tauself_var_1208 + taufor_var_1207 + adjcolco2_var_1204 * absco2_var_1209 + wx_var_1174(jl_var_1220, 2, lay_var_1202) * cfc11adj(ig_var_1201) + wx_var_1174(jl_var_1220, 3, lay_var_1202) * cfc12_var_181(ig_var_1201)
          fracs_var_1190(jl_var_1220, 68 + ig_var_1201, lay_var_1202) = fracrefa_var_180(ig_var_1201)
        END DO
      END DO
      ixc0_var_1217 = kfdia_var_1171 - kidia_var_1170 + 1 - ixc0_var_1217
      DO ig_var_1201 = 1, 8
        DO ixp_var_1218 = 1, ixc0_var_1217
          jl_var_1220 = ixhigh_var_1214(ixp_var_1218, lay_var_1202)
          taug_var_1173(jl_var_1220, 68 + ig_var_1201, lay_var_1202) = 0.0D0 + wx_var_1174(jl_var_1220, 2, lay_var_1202) * cfc11adj(ig_var_1201) + wx_var_1174(jl_var_1220, 3, lay_var_1202) * cfc12_var_181(ig_var_1201)
          fracs_var_1190(jl_var_1220, 68 + ig_var_1201, lay_var_1202) = fracrefa_var_180(ig_var_1201)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol6
SUBROUTINE rrtm_taumol7(kidia_var_1221, kfdia_var_1222, klev_var_1223, taug_var_1224, p_tauaerl_var_1225, fac00_var_1226, fac01_var_1227, fac10_var_1228, fac11_var_1229, forfac_var_1245, forfrac_var_1244, indfor_var_1243, jp_var_1230, jt_var_1231, jt1_var_1232, oneminus_var_1233, colh2o_var_1234, colo3_var_1235, colco2_var_1236, coldry_var_1237, laytrop_var_1238, selffac_var_1239, selffrac_var_1240, indself_var_1241, fracs_var_1242, rat_h2oo3, rat_h2oo3_1, minorfrac_var_1246, indminor_var_1247)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng7
  USE yoerrta7, ONLY: absa_var_188, absb_var_189, forref_var_193, fracrefa_var_186, fracrefb_var_187, ka_mco2_var_191, kb_mco2_var_192, selfref_var_190
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1221
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1222
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1223
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1224(kidia_var_1221 : kfdia_var_1222, 140, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1225(kidia_var_1221 : kfdia_var_1222, klev_var_1223, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1226(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1227(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1228(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1229(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1230(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1231(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1232(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1233
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1234(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1235(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1236(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1237(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1238(kidia_var_1221 : kfdia_var_1222)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1239(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1240(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1241(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1242(kidia_var_1221 : kfdia_var_1222, 140, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oo3(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oo3_1(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1243(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1244(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1245(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1246(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1247(kidia_var_1221 : kfdia_var_1222, klev_var_1223)
  REAL(KIND = 8) :: speccomb_var_1248, speccomb1_var_1249, speccomb_mco2_var_1250, speccomb_planck_var_1251
  INTEGER(KIND = 4) :: ind0_var_1252, ind1_var_1253, inds_var_1254, indf_var_1255, indm_var_1256
  INTEGER(KIND = 4) :: ig_var_1257, js_var_1258, lay_var_1259, js1_var_1260, jpl_var_1261, jmco2_var_1262
  REAL(KIND = 8) :: refrat_planck_a_var_1263, refrat_m_a_var_1264
  REAL(KIND = 8) :: chi_co2_var_1265, ratco2_var_1266, adjfac_var_1267, adjcolco2_var_1268
  REAL(KIND = 8) :: fac000_var_1269, fac100_var_1270, fac200_var_1271, fac010_var_1272, fac110_var_1273, fac210_var_1274, fac001_var_1275, fac101_var_1276, fac201_var_1277, fac011_var_1278, fac111_var_1279, fac211_var_1280
  REAL(KIND = 8) :: p_var_1281, p4_var_1282, fk0_var_1283, fk1_var_1284, fk2_var_1285
  REAL(KIND = 8) :: taufor_var_1286, tauself_var_1287, tau_major_var_1288(12), tau_major1_var_1289(12), co2m1_var_1290, co2m2_var_1291, absco2_var_1292
  REAL(KIND = 8) :: fs_var_1293, specmult_var_1294, specparm_var_1295, fs1_var_1296, specmult1_var_1297, specparm1_var_1298, fpl_var_1299, specmult_planck_var_1300, specparm_planck_var_1301, fmco2_var_1302, specmult_mco2_var_1303, specparm_mco2_var_1304
  INTEGER(KIND = 4) :: laytrop_min_var_1305, laytrop_max_var_1306
  INTEGER(KIND = 4) :: ixc_var_1307(klev_var_1223), ixlow_var_1308(kfdia_var_1222, klev_var_1223), ixhigh_var_1309(kfdia_var_1222, klev_var_1223)
  INTEGER(KIND = 4) :: ich_var_1310, icl_var_1311, ixc0_var_1312, ixp_var_1313, jc_var_1314, jl_var_1315
  laytrop_min_var_1305 = MINVAL(laytrop_var_1238)
  laytrop_max_var_1306 = MAXVAL(laytrop_var_1238)
  ixlow_var_1308 = 0
  ixhigh_var_1309 = 0
  ixc_var_1307 = 0
  DO lay_var_1259 = laytrop_min_var_1305 + 1, laytrop_max_var_1306
    icl_var_1311 = 0
    ich_var_1310 = 0
    DO jc_var_1314 = kidia_var_1221, kfdia_var_1222
      IF (lay_var_1259 <= laytrop_var_1238(jc_var_1314)) THEN
        icl_var_1311 = icl_var_1311 + 1
        ixlow_var_1308(icl_var_1311, lay_var_1259) = jc_var_1314
      ELSE
        ich_var_1310 = ich_var_1310 + 1
        ixhigh_var_1309(ich_var_1310, lay_var_1259) = jc_var_1314
      END IF
    END DO
    ixc_var_1307(lay_var_1259) = icl_var_1311
  END DO
  refrat_planck_a_var_1263 = chi_mls(1, 3) / chi_mls(3, 3)
  refrat_m_a_var_1264 = chi_mls(1, 3) / chi_mls(3, 3)
  DO lay_var_1259 = 1, laytrop_min_var_1305
    DO jl_var_1315 = kidia_var_1221, kfdia_var_1222
      speccomb_var_1248 = colh2o_var_1234(jl_var_1315, lay_var_1259) + rat_h2oo3(jl_var_1315, lay_var_1259) * colo3_var_1235(jl_var_1315, lay_var_1259)
      specparm_var_1295 = MIN(colh2o_var_1234(jl_var_1315, lay_var_1259) / speccomb_var_1248, 0.999999D0)
      specmult_var_1294 = 8.0D0 * (specparm_var_1295)
      js_var_1258 = 1 + INT(specmult_var_1294)
      fs_var_1293 = ((specmult_var_1294) - AINT((specmult_var_1294)))
      speccomb1_var_1249 = colh2o_var_1234(jl_var_1315, lay_var_1259) + rat_h2oo3_1(jl_var_1315, lay_var_1259) * colo3_var_1235(jl_var_1315, lay_var_1259)
      specparm1_var_1298 = MIN(colh2o_var_1234(jl_var_1315, lay_var_1259) / speccomb1_var_1249, 0.999999D0)
      specmult1_var_1297 = 8.0D0 * (specparm1_var_1298)
      js1_var_1260 = 1 + INT(specmult1_var_1297)
      fs1_var_1296 = ((specmult1_var_1297) - AINT((specmult1_var_1297)))
      speccomb_mco2_var_1250 = colh2o_var_1234(jl_var_1315, lay_var_1259) + refrat_m_a_var_1264 * colo3_var_1235(jl_var_1315, lay_var_1259)
      specparm_mco2_var_1304 = MIN(colh2o_var_1234(jl_var_1315, lay_var_1259) / speccomb_mco2_var_1250, 0.999999D0)
      specmult_mco2_var_1303 = 8.0D0 * specparm_mco2_var_1304
      jmco2_var_1262 = 1 + INT(specmult_mco2_var_1303)
      fmco2_var_1302 = ((specmult_mco2_var_1303) - AINT((specmult_mco2_var_1303)))
      chi_co2_var_1265 = colco2_var_1236(jl_var_1315, lay_var_1259) / (coldry_var_1237(jl_var_1315, lay_var_1259))
      ratco2_var_1266 = 1D+20 * chi_co2_var_1265 / chi_mls(2, jp_var_1230(jl_var_1315, lay_var_1259) + 1)
      IF (ratco2_var_1266 .GT. 3.0D0) THEN
        adjfac_var_1267 = 3.0D0 + (ratco2_var_1266 - 3.0D0) ** 0.79D0
        adjcolco2_var_1268 = adjfac_var_1267 * chi_mls(2, jp_var_1230(jl_var_1315, lay_var_1259) + 1) * coldry_var_1237(jl_var_1315, lay_var_1259) * 1D-20
      ELSE
        adjcolco2_var_1268 = colco2_var_1236(jl_var_1315, lay_var_1259)
      END IF
      speccomb_planck_var_1251 = colh2o_var_1234(jl_var_1315, lay_var_1259) + refrat_planck_a_var_1263 * colo3_var_1235(jl_var_1315, lay_var_1259)
      specparm_planck_var_1301 = MIN(colh2o_var_1234(jl_var_1315, lay_var_1259) / speccomb_planck_var_1251, 0.999999D0)
      specmult_planck_var_1300 = 8.0D0 * specparm_planck_var_1301
      jpl_var_1261 = 1 + INT(specmult_planck_var_1300)
      fpl_var_1299 = ((specmult_planck_var_1300) - AINT((specmult_planck_var_1300)))
      ind0_var_1252 = ((jp_var_1230(jl_var_1315, lay_var_1259) - 1) * 5 + (jt_var_1231(jl_var_1315, lay_var_1259) - 1)) * nspa_var_216(7) + js_var_1258
      ind1_var_1253 = (jp_var_1230(jl_var_1315, lay_var_1259) * 5 + (jt1_var_1232(jl_var_1315, lay_var_1259) - 1)) * nspa_var_216(7) + js1_var_1260
      inds_var_1254 = indself_var_1241(jl_var_1315, lay_var_1259)
      indf_var_1255 = indfor_var_1243(jl_var_1315, lay_var_1259)
      indm_var_1256 = indminor_var_1247(jl_var_1315, lay_var_1259)
      IF (specparm_var_1295 .LT. 0.125D0) THEN
        p_var_1281 = fs_var_1293 - 1.0D0
        p4_var_1282 = p_var_1281 ** 4
        fk0_var_1283 = p4_var_1282
        fk1_var_1284 = 1.0D0 - p_var_1281 - 2.0D0 * p4_var_1282
        fk2_var_1285 = p_var_1281 + p4_var_1282
        fac000_var_1269 = fk0_var_1283 * fac00_var_1226(jl_var_1315, lay_var_1259)
        fac100_var_1270 = fk1_var_1284 * fac00_var_1226(jl_var_1315, lay_var_1259)
        fac200_var_1271 = fk2_var_1285 * fac00_var_1226(jl_var_1315, lay_var_1259)
        fac010_var_1272 = fk0_var_1283 * fac10_var_1228(jl_var_1315, lay_var_1259)
        fac110_var_1273 = fk1_var_1284 * fac10_var_1228(jl_var_1315, lay_var_1259)
        fac210_var_1274 = fk2_var_1285 * fac10_var_1228(jl_var_1315, lay_var_1259)
      ELSE IF (specparm_var_1295 .GT. 0.875D0) THEN
        p_var_1281 = - fs_var_1293
        p4_var_1282 = p_var_1281 ** 4
        fk0_var_1283 = p4_var_1282
        fk1_var_1284 = 1.0D0 - p_var_1281 - 2.0D0 * p4_var_1282
        fk2_var_1285 = p_var_1281 + p4_var_1282
        fac000_var_1269 = fk0_var_1283 * fac00_var_1226(jl_var_1315, lay_var_1259)
        fac100_var_1270 = fk1_var_1284 * fac00_var_1226(jl_var_1315, lay_var_1259)
        fac200_var_1271 = fk2_var_1285 * fac00_var_1226(jl_var_1315, lay_var_1259)
        fac010_var_1272 = fk0_var_1283 * fac10_var_1228(jl_var_1315, lay_var_1259)
        fac110_var_1273 = fk1_var_1284 * fac10_var_1228(jl_var_1315, lay_var_1259)
        fac210_var_1274 = fk2_var_1285 * fac10_var_1228(jl_var_1315, lay_var_1259)
      ELSE
        fac000_var_1269 = (1.0D0 - fs_var_1293) * fac00_var_1226(jl_var_1315, lay_var_1259)
        fac010_var_1272 = (1.0D0 - fs_var_1293) * fac10_var_1228(jl_var_1315, lay_var_1259)
        fac100_var_1270 = fs_var_1293 * fac00_var_1226(jl_var_1315, lay_var_1259)
        fac110_var_1273 = fs_var_1293 * fac10_var_1228(jl_var_1315, lay_var_1259)
        fac200_var_1271 = 0.0D0
        fac210_var_1274 = 0.0D0
      END IF
      IF (specparm1_var_1298 .LT. 0.125D0) THEN
        p_var_1281 = fs1_var_1296 - 1.0D0
        p4_var_1282 = p_var_1281 ** 4
        fk0_var_1283 = p4_var_1282
        fk1_var_1284 = 1.0D0 - p_var_1281 - 2.0D0 * p4_var_1282
        fk2_var_1285 = p_var_1281 + p4_var_1282
        fac001_var_1275 = fk0_var_1283 * fac01_var_1227(jl_var_1315, lay_var_1259)
        fac101_var_1276 = fk1_var_1284 * fac01_var_1227(jl_var_1315, lay_var_1259)
        fac201_var_1277 = fk2_var_1285 * fac01_var_1227(jl_var_1315, lay_var_1259)
        fac011_var_1278 = fk0_var_1283 * fac11_var_1229(jl_var_1315, lay_var_1259)
        fac111_var_1279 = fk1_var_1284 * fac11_var_1229(jl_var_1315, lay_var_1259)
        fac211_var_1280 = fk2_var_1285 * fac11_var_1229(jl_var_1315, lay_var_1259)
      ELSE IF (specparm1_var_1298 .GT. 0.875D0) THEN
        p_var_1281 = - fs1_var_1296
        p4_var_1282 = p_var_1281 ** 4
        fk0_var_1283 = p4_var_1282
        fk1_var_1284 = 1.0D0 - p_var_1281 - 2.0D0 * p4_var_1282
        fk2_var_1285 = p_var_1281 + p4_var_1282
        fac001_var_1275 = fk0_var_1283 * fac01_var_1227(jl_var_1315, lay_var_1259)
        fac101_var_1276 = fk1_var_1284 * fac01_var_1227(jl_var_1315, lay_var_1259)
        fac201_var_1277 = fk2_var_1285 * fac01_var_1227(jl_var_1315, lay_var_1259)
        fac011_var_1278 = fk0_var_1283 * fac11_var_1229(jl_var_1315, lay_var_1259)
        fac111_var_1279 = fk1_var_1284 * fac11_var_1229(jl_var_1315, lay_var_1259)
        fac211_var_1280 = fk2_var_1285 * fac11_var_1229(jl_var_1315, lay_var_1259)
      ELSE
        fac001_var_1275 = (1.0D0 - fs1_var_1296) * fac01_var_1227(jl_var_1315, lay_var_1259)
        fac011_var_1278 = (1.0D0 - fs1_var_1296) * fac11_var_1229(jl_var_1315, lay_var_1259)
        fac101_var_1276 = fs1_var_1296 * fac01_var_1227(jl_var_1315, lay_var_1259)
        fac111_var_1279 = fs1_var_1296 * fac11_var_1229(jl_var_1315, lay_var_1259)
        fac201_var_1277 = 0.0D0
        fac211_var_1280 = 0.0D0
      END IF
      IF (specparm_var_1295 .LT. 0.125D0) THEN
        tau_major_var_1288(1 : ng7) = speccomb_var_1248 * (fac000_var_1269 * absa_var_188(ind0_var_1252, 1 : 12) + fac100_var_1270 * absa_var_188(ind0_var_1252 + 1, 1 : 12) + fac200_var_1271 * absa_var_188(ind0_var_1252 + 2, 1 : 12) + fac010_var_1272 * absa_var_188(ind0_var_1252 + 9, 1 : 12) + fac110_var_1273 * absa_var_188(ind0_var_1252 + 10, 1 : 12) + fac210_var_1274 * absa_var_188(ind0_var_1252 + 11, 1 : 12))
      ELSE IF (specparm_var_1295 .GT. 0.875D0) THEN
        tau_major_var_1288(1 : ng7) = speccomb_var_1248 * (fac200_var_1271 * absa_var_188(ind0_var_1252 - 1, 1 : 12) + fac100_var_1270 * absa_var_188(ind0_var_1252, 1 : 12) + fac000_var_1269 * absa_var_188(ind0_var_1252 + 1, 1 : 12) + fac210_var_1274 * absa_var_188(ind0_var_1252 + 8, 1 : 12) + fac110_var_1273 * absa_var_188(ind0_var_1252 + 9, 1 : 12) + fac010_var_1272 * absa_var_188(ind0_var_1252 + 10, 1 : 12))
      ELSE
        tau_major_var_1288(1 : ng7) = speccomb_var_1248 * (fac000_var_1269 * absa_var_188(ind0_var_1252, 1 : 12) + fac100_var_1270 * absa_var_188(ind0_var_1252 + 1, 1 : 12) + fac010_var_1272 * absa_var_188(ind0_var_1252 + 9, 1 : 12) + fac110_var_1273 * absa_var_188(ind0_var_1252 + 10, 1 : 12))
      END IF
      IF (specparm1_var_1298 .LT. 0.125D0) THEN
        tau_major1_var_1289(1 : ng7) = speccomb1_var_1249 * (fac001_var_1275 * absa_var_188(ind1_var_1253, 1 : 12) + fac101_var_1276 * absa_var_188(ind1_var_1253 + 1, 1 : 12) + fac201_var_1277 * absa_var_188(ind1_var_1253 + 2, 1 : 12) + fac011_var_1278 * absa_var_188(ind1_var_1253 + 9, 1 : 12) + fac111_var_1279 * absa_var_188(ind1_var_1253 + 10, 1 : 12) + fac211_var_1280 * absa_var_188(ind1_var_1253 + 11, 1 : 12))
      ELSE IF (specparm1_var_1298 .GT. 0.875D0) THEN
        tau_major1_var_1289(1 : ng7) = speccomb1_var_1249 * (fac201_var_1277 * absa_var_188(ind1_var_1253 - 1, 1 : 12) + fac101_var_1276 * absa_var_188(ind1_var_1253, 1 : 12) + fac001_var_1275 * absa_var_188(ind1_var_1253 + 1, 1 : 12) + fac211_var_1280 * absa_var_188(ind1_var_1253 + 8, 1 : 12) + fac111_var_1279 * absa_var_188(ind1_var_1253 + 9, 1 : 12) + fac011_var_1278 * absa_var_188(ind1_var_1253 + 10, 1 : 12))
      ELSE
        tau_major1_var_1289(1 : ng7) = speccomb1_var_1249 * (fac001_var_1275 * absa_var_188(ind1_var_1253, 1 : 12) + fac101_var_1276 * absa_var_188(ind1_var_1253 + 1, 1 : 12) + fac011_var_1278 * absa_var_188(ind1_var_1253 + 9, 1 : 12) + fac111_var_1279 * absa_var_188(ind1_var_1253 + 10, 1 : 12))
      END IF
      DO ig_var_1257 = 1, 12
        tauself_var_1287 = selffac_var_1239(jl_var_1315, lay_var_1259) * (selfref_var_190(inds_var_1254, ig_var_1257) + selffrac_var_1240(jl_var_1315, lay_var_1259) * (selfref_var_190(inds_var_1254 + 1, ig_var_1257) - selfref_var_190(inds_var_1254, ig_var_1257)))
        taufor_var_1286 = forfac_var_1245(jl_var_1315, lay_var_1259) * (forref_var_193(indf_var_1255, ig_var_1257) + forfrac_var_1244(jl_var_1315, lay_var_1259) * (forref_var_193(indf_var_1255 + 1, ig_var_1257) - forref_var_193(indf_var_1255, ig_var_1257)))
        co2m1_var_1290 = ka_mco2_var_191(jmco2_var_1262, indm_var_1256, ig_var_1257) + fmco2_var_1302 * (ka_mco2_var_191(jmco2_var_1262 + 1, indm_var_1256, ig_var_1257) - ka_mco2_var_191(jmco2_var_1262, indm_var_1256, ig_var_1257))
        co2m2_var_1291 = ka_mco2_var_191(jmco2_var_1262, indm_var_1256 + 1, ig_var_1257) + fmco2_var_1302 * (ka_mco2_var_191(jmco2_var_1262 + 1, indm_var_1256 + 1, ig_var_1257) - ka_mco2_var_191(jmco2_var_1262, indm_var_1256 + 1, ig_var_1257))
        absco2_var_1292 = co2m1_var_1290 + minorfrac_var_1246(jl_var_1315, lay_var_1259) * (co2m2_var_1291 - co2m1_var_1290)
        taug_var_1224(jl_var_1315, 76 + ig_var_1257, lay_var_1259) = tau_major_var_1288(ig_var_1257) + tau_major1_var_1289(ig_var_1257) + tauself_var_1287 + taufor_var_1286 + adjcolco2_var_1268 * absco2_var_1292
        fracs_var_1242(jl_var_1315, 76 + ig_var_1257, lay_var_1259) = fracrefa_var_186(ig_var_1257, jpl_var_1261) + fpl_var_1299 * (fracrefa_var_186(ig_var_1257, jpl_var_1261 + 1) - fracrefa_var_186(ig_var_1257, jpl_var_1261))
      END DO
    END DO
  END DO
  DO lay_var_1259 = laytrop_max_var_1306 + 1, klev_var_1223
    DO jl_var_1315 = kidia_var_1221, kfdia_var_1222
      chi_co2_var_1265 = colco2_var_1236(jl_var_1315, lay_var_1259) / (coldry_var_1237(jl_var_1315, lay_var_1259))
      ratco2_var_1266 = 1D+20 * chi_co2_var_1265 / chi_mls(2, jp_var_1230(jl_var_1315, lay_var_1259) + 1)
      IF (ratco2_var_1266 .GT. 3.0D0) THEN
        adjfac_var_1267 = 2.0D0 + (ratco2_var_1266 - 2.0D0) ** 0.79D0
        adjcolco2_var_1268 = adjfac_var_1267 * chi_mls(2, jp_var_1230(jl_var_1315, lay_var_1259) + 1) * coldry_var_1237(jl_var_1315, lay_var_1259) * 1D-20
      ELSE
        adjcolco2_var_1268 = colco2_var_1236(jl_var_1315, lay_var_1259)
      END IF
      ind0_var_1252 = ((jp_var_1230(jl_var_1315, lay_var_1259) - 13) * 5 + (jt_var_1231(jl_var_1315, lay_var_1259) - 1)) * nspb_var_217(7) + 1
      ind1_var_1253 = ((jp_var_1230(jl_var_1315, lay_var_1259) - 12) * 5 + (jt1_var_1232(jl_var_1315, lay_var_1259) - 1)) * nspb_var_217(7) + 1
      indm_var_1256 = indminor_var_1247(jl_var_1315, lay_var_1259)
      DO ig_var_1257 = 1, 12
        absco2_var_1292 = kb_mco2_var_192(indm_var_1256, ig_var_1257) + minorfrac_var_1246(jl_var_1315, lay_var_1259) * (kb_mco2_var_192(indm_var_1256 + 1, ig_var_1257) - kb_mco2_var_192(indm_var_1256, ig_var_1257))
        taug_var_1224(jl_var_1315, 76 + ig_var_1257, lay_var_1259) = colo3_var_1235(jl_var_1315, lay_var_1259) * (fac00_var_1226(jl_var_1315, lay_var_1259) * absb_var_189(ind0_var_1252, ig_var_1257) + fac10_var_1228(jl_var_1315, lay_var_1259) * absb_var_189(ind0_var_1252 + 1, ig_var_1257) + fac01_var_1227(jl_var_1315, lay_var_1259) * absb_var_189(ind1_var_1253, ig_var_1257) + fac11_var_1229(jl_var_1315, lay_var_1259) * absb_var_189(ind1_var_1253 + 1, ig_var_1257)) + adjcolco2_var_1268 * absco2_var_1292
        fracs_var_1242(jl_var_1315, 76 + ig_var_1257, lay_var_1259) = fracrefb_var_187(ig_var_1257)
      END DO
    END DO
  END DO
  DO lay_var_1259 = laytrop_max_var_1306 + 1, klev_var_1223
    DO jl_var_1315 = kidia_var_1221, kfdia_var_1222
      taug_var_1224(jl_var_1315, 82, lay_var_1259) = taug_var_1224(jl_var_1315, 82, lay_var_1259) * 0.92D0
      taug_var_1224(jl_var_1315, 83, lay_var_1259) = taug_var_1224(jl_var_1315, 83, lay_var_1259) * 0.88D0
      taug_var_1224(jl_var_1315, 84, lay_var_1259) = taug_var_1224(jl_var_1315, 84, lay_var_1259) * 1.07D0
      taug_var_1224(jl_var_1315, 85, lay_var_1259) = taug_var_1224(jl_var_1315, 85, lay_var_1259) * 1.1D0
      taug_var_1224(jl_var_1315, 86, lay_var_1259) = taug_var_1224(jl_var_1315, 86, lay_var_1259) * 0.99D0
      taug_var_1224(jl_var_1315, 87, lay_var_1259) = taug_var_1224(jl_var_1315, 87, lay_var_1259) * 0.855D0
    END DO
  END DO
  IF (laytrop_max_var_1306 /= laytrop_min_var_1305) THEN
    DO lay_var_1259 = laytrop_min_var_1305 + 1, laytrop_max_var_1306
      ixc0_var_1312 = ixc_var_1307(lay_var_1259)
      DO ixp_var_1313 = 1, ixc0_var_1312
        jl_var_1315 = ixlow_var_1308(ixp_var_1313, lay_var_1259)
        speccomb_var_1248 = colh2o_var_1234(jl_var_1315, lay_var_1259) + rat_h2oo3(jl_var_1315, lay_var_1259) * colo3_var_1235(jl_var_1315, lay_var_1259)
        specparm_var_1295 = MIN(colh2o_var_1234(jl_var_1315, lay_var_1259) / speccomb_var_1248, 0.999999D0)
        specmult_var_1294 = 8.0D0 * (specparm_var_1295)
        js_var_1258 = 1 + INT(specmult_var_1294)
        fs_var_1293 = ((specmult_var_1294) - AINT((specmult_var_1294)))
        speccomb1_var_1249 = colh2o_var_1234(jl_var_1315, lay_var_1259) + rat_h2oo3_1(jl_var_1315, lay_var_1259) * colo3_var_1235(jl_var_1315, lay_var_1259)
        specparm1_var_1298 = MIN(colh2o_var_1234(jl_var_1315, lay_var_1259) / speccomb1_var_1249, 0.999999D0)
        specmult1_var_1297 = 8.0D0 * (specparm1_var_1298)
        js1_var_1260 = 1 + INT(specmult1_var_1297)
        fs1_var_1296 = ((specmult1_var_1297) - AINT((specmult1_var_1297)))
        speccomb_mco2_var_1250 = colh2o_var_1234(jl_var_1315, lay_var_1259) + refrat_m_a_var_1264 * colo3_var_1235(jl_var_1315, lay_var_1259)
        specparm_mco2_var_1304 = MIN(colh2o_var_1234(jl_var_1315, lay_var_1259) / speccomb_mco2_var_1250, 0.999999D0)
        specmult_mco2_var_1303 = 8.0D0 * specparm_mco2_var_1304
        jmco2_var_1262 = 1 + INT(specmult_mco2_var_1303)
        fmco2_var_1302 = ((specmult_mco2_var_1303) - AINT((specmult_mco2_var_1303)))
        chi_co2_var_1265 = colco2_var_1236(jl_var_1315, lay_var_1259) / (coldry_var_1237(jl_var_1315, lay_var_1259))
        ratco2_var_1266 = 1D+20 * chi_co2_var_1265 / chi_mls(2, jp_var_1230(jl_var_1315, lay_var_1259) + 1)
        IF (ratco2_var_1266 .GT. 3.0D0) THEN
          adjfac_var_1267 = 3.0D0 + (ratco2_var_1266 - 3.0D0) ** 0.79D0
          adjcolco2_var_1268 = adjfac_var_1267 * chi_mls(2, jp_var_1230(jl_var_1315, lay_var_1259) + 1) * coldry_var_1237(jl_var_1315, lay_var_1259) * 1D-20
        ELSE
          adjcolco2_var_1268 = colco2_var_1236(jl_var_1315, lay_var_1259)
        END IF
        speccomb_planck_var_1251 = colh2o_var_1234(jl_var_1315, lay_var_1259) + refrat_planck_a_var_1263 * colo3_var_1235(jl_var_1315, lay_var_1259)
        specparm_planck_var_1301 = MIN(colh2o_var_1234(jl_var_1315, lay_var_1259) / speccomb_planck_var_1251, 0.999999D0)
        specmult_planck_var_1300 = 8.0D0 * specparm_planck_var_1301
        jpl_var_1261 = 1 + INT(specmult_planck_var_1300)
        fpl_var_1299 = ((specmult_planck_var_1300) - AINT((specmult_planck_var_1300)))
        ind0_var_1252 = ((jp_var_1230(jl_var_1315, lay_var_1259) - 1) * 5 + (jt_var_1231(jl_var_1315, lay_var_1259) - 1)) * nspa_var_216(7) + js_var_1258
        ind1_var_1253 = (jp_var_1230(jl_var_1315, lay_var_1259) * 5 + (jt1_var_1232(jl_var_1315, lay_var_1259) - 1)) * nspa_var_216(7) + js1_var_1260
        inds_var_1254 = indself_var_1241(jl_var_1315, lay_var_1259)
        indf_var_1255 = indfor_var_1243(jl_var_1315, lay_var_1259)
        indm_var_1256 = indminor_var_1247(jl_var_1315, lay_var_1259)
        IF (specparm_var_1295 .LT. 0.125D0) THEN
          p_var_1281 = fs_var_1293 - 1.0D0
          p4_var_1282 = p_var_1281 ** 4
          fk0_var_1283 = p4_var_1282
          fk1_var_1284 = 1.0D0 - p_var_1281 - 2.0D0 * p4_var_1282
          fk2_var_1285 = p_var_1281 + p4_var_1282
          fac000_var_1269 = fk0_var_1283 * fac00_var_1226(jl_var_1315, lay_var_1259)
          fac100_var_1270 = fk1_var_1284 * fac00_var_1226(jl_var_1315, lay_var_1259)
          fac200_var_1271 = fk2_var_1285 * fac00_var_1226(jl_var_1315, lay_var_1259)
          fac010_var_1272 = fk0_var_1283 * fac10_var_1228(jl_var_1315, lay_var_1259)
          fac110_var_1273 = fk1_var_1284 * fac10_var_1228(jl_var_1315, lay_var_1259)
          fac210_var_1274 = fk2_var_1285 * fac10_var_1228(jl_var_1315, lay_var_1259)
        ELSE IF (specparm_var_1295 .GT. 0.875D0) THEN
          p_var_1281 = - fs_var_1293
          p4_var_1282 = p_var_1281 ** 4
          fk0_var_1283 = p4_var_1282
          fk1_var_1284 = 1.0D0 - p_var_1281 - 2.0D0 * p4_var_1282
          fk2_var_1285 = p_var_1281 + p4_var_1282
          fac000_var_1269 = fk0_var_1283 * fac00_var_1226(jl_var_1315, lay_var_1259)
          fac100_var_1270 = fk1_var_1284 * fac00_var_1226(jl_var_1315, lay_var_1259)
          fac200_var_1271 = fk2_var_1285 * fac00_var_1226(jl_var_1315, lay_var_1259)
          fac010_var_1272 = fk0_var_1283 * fac10_var_1228(jl_var_1315, lay_var_1259)
          fac110_var_1273 = fk1_var_1284 * fac10_var_1228(jl_var_1315, lay_var_1259)
          fac210_var_1274 = fk2_var_1285 * fac10_var_1228(jl_var_1315, lay_var_1259)
        ELSE
          fac000_var_1269 = (1.0D0 - fs_var_1293) * fac00_var_1226(jl_var_1315, lay_var_1259)
          fac010_var_1272 = (1.0D0 - fs_var_1293) * fac10_var_1228(jl_var_1315, lay_var_1259)
          fac100_var_1270 = fs_var_1293 * fac00_var_1226(jl_var_1315, lay_var_1259)
          fac110_var_1273 = fs_var_1293 * fac10_var_1228(jl_var_1315, lay_var_1259)
          fac200_var_1271 = 0.0D0
          fac210_var_1274 = 0.0D0
        END IF
        IF (specparm1_var_1298 .LT. 0.125D0) THEN
          p_var_1281 = fs1_var_1296 - 1.0D0
          p4_var_1282 = p_var_1281 ** 4
          fk0_var_1283 = p4_var_1282
          fk1_var_1284 = 1.0D0 - p_var_1281 - 2.0D0 * p4_var_1282
          fk2_var_1285 = p_var_1281 + p4_var_1282
          fac001_var_1275 = fk0_var_1283 * fac01_var_1227(jl_var_1315, lay_var_1259)
          fac101_var_1276 = fk1_var_1284 * fac01_var_1227(jl_var_1315, lay_var_1259)
          fac201_var_1277 = fk2_var_1285 * fac01_var_1227(jl_var_1315, lay_var_1259)
          fac011_var_1278 = fk0_var_1283 * fac11_var_1229(jl_var_1315, lay_var_1259)
          fac111_var_1279 = fk1_var_1284 * fac11_var_1229(jl_var_1315, lay_var_1259)
          fac211_var_1280 = fk2_var_1285 * fac11_var_1229(jl_var_1315, lay_var_1259)
        ELSE IF (specparm1_var_1298 .GT. 0.875D0) THEN
          p_var_1281 = - fs1_var_1296
          p4_var_1282 = p_var_1281 ** 4
          fk0_var_1283 = p4_var_1282
          fk1_var_1284 = 1.0D0 - p_var_1281 - 2.0D0 * p4_var_1282
          fk2_var_1285 = p_var_1281 + p4_var_1282
          fac001_var_1275 = fk0_var_1283 * fac01_var_1227(jl_var_1315, lay_var_1259)
          fac101_var_1276 = fk1_var_1284 * fac01_var_1227(jl_var_1315, lay_var_1259)
          fac201_var_1277 = fk2_var_1285 * fac01_var_1227(jl_var_1315, lay_var_1259)
          fac011_var_1278 = fk0_var_1283 * fac11_var_1229(jl_var_1315, lay_var_1259)
          fac111_var_1279 = fk1_var_1284 * fac11_var_1229(jl_var_1315, lay_var_1259)
          fac211_var_1280 = fk2_var_1285 * fac11_var_1229(jl_var_1315, lay_var_1259)
        ELSE
          fac001_var_1275 = (1.0D0 - fs1_var_1296) * fac01_var_1227(jl_var_1315, lay_var_1259)
          fac011_var_1278 = (1.0D0 - fs1_var_1296) * fac11_var_1229(jl_var_1315, lay_var_1259)
          fac101_var_1276 = fs1_var_1296 * fac01_var_1227(jl_var_1315, lay_var_1259)
          fac111_var_1279 = fs1_var_1296 * fac11_var_1229(jl_var_1315, lay_var_1259)
          fac201_var_1277 = 0.0D0
          fac211_var_1280 = 0.0D0
        END IF
        IF (specparm_var_1295 .LT. 0.125D0) THEN
          tau_major_var_1288(1 : ng7) = speccomb_var_1248 * (fac000_var_1269 * absa_var_188(ind0_var_1252, 1 : 12) + fac100_var_1270 * absa_var_188(ind0_var_1252 + 1, 1 : 12) + fac200_var_1271 * absa_var_188(ind0_var_1252 + 2, 1 : 12) + fac010_var_1272 * absa_var_188(ind0_var_1252 + 9, 1 : 12) + fac110_var_1273 * absa_var_188(ind0_var_1252 + 10, 1 : 12) + fac210_var_1274 * absa_var_188(ind0_var_1252 + 11, 1 : 12))
        ELSE IF (specparm_var_1295 .GT. 0.875D0) THEN
          tau_major_var_1288(1 : ng7) = speccomb_var_1248 * (fac200_var_1271 * absa_var_188(ind0_var_1252 - 1, 1 : 12) + fac100_var_1270 * absa_var_188(ind0_var_1252, 1 : 12) + fac000_var_1269 * absa_var_188(ind0_var_1252 + 1, 1 : 12) + fac210_var_1274 * absa_var_188(ind0_var_1252 + 8, 1 : 12) + fac110_var_1273 * absa_var_188(ind0_var_1252 + 9, 1 : 12) + fac010_var_1272 * absa_var_188(ind0_var_1252 + 10, 1 : 12))
        ELSE
          tau_major_var_1288(1 : ng7) = speccomb_var_1248 * (fac000_var_1269 * absa_var_188(ind0_var_1252, 1 : 12) + fac100_var_1270 * absa_var_188(ind0_var_1252 + 1, 1 : 12) + fac010_var_1272 * absa_var_188(ind0_var_1252 + 9, 1 : 12) + fac110_var_1273 * absa_var_188(ind0_var_1252 + 10, 1 : 12))
        END IF
        IF (specparm1_var_1298 .LT. 0.125D0) THEN
          tau_major1_var_1289(1 : ng7) = speccomb1_var_1249 * (fac001_var_1275 * absa_var_188(ind1_var_1253, 1 : 12) + fac101_var_1276 * absa_var_188(ind1_var_1253 + 1, 1 : 12) + fac201_var_1277 * absa_var_188(ind1_var_1253 + 2, 1 : 12) + fac011_var_1278 * absa_var_188(ind1_var_1253 + 9, 1 : 12) + fac111_var_1279 * absa_var_188(ind1_var_1253 + 10, 1 : 12) + fac211_var_1280 * absa_var_188(ind1_var_1253 + 11, 1 : 12))
        ELSE IF (specparm1_var_1298 .GT. 0.875D0) THEN
          tau_major1_var_1289(1 : ng7) = speccomb1_var_1249 * (fac201_var_1277 * absa_var_188(ind1_var_1253 - 1, 1 : 12) + fac101_var_1276 * absa_var_188(ind1_var_1253, 1 : 12) + fac001_var_1275 * absa_var_188(ind1_var_1253 + 1, 1 : 12) + fac211_var_1280 * absa_var_188(ind1_var_1253 + 8, 1 : 12) + fac111_var_1279 * absa_var_188(ind1_var_1253 + 9, 1 : 12) + fac011_var_1278 * absa_var_188(ind1_var_1253 + 10, 1 : 12))
        ELSE
          tau_major1_var_1289(1 : ng7) = speccomb1_var_1249 * (fac001_var_1275 * absa_var_188(ind1_var_1253, 1 : 12) + fac101_var_1276 * absa_var_188(ind1_var_1253 + 1, 1 : 12) + fac011_var_1278 * absa_var_188(ind1_var_1253 + 9, 1 : 12) + fac111_var_1279 * absa_var_188(ind1_var_1253 + 10, 1 : 12))
        END IF
        DO ig_var_1257 = 1, 12
          tauself_var_1287 = selffac_var_1239(jl_var_1315, lay_var_1259) * (selfref_var_190(inds_var_1254, ig_var_1257) + selffrac_var_1240(jl_var_1315, lay_var_1259) * (selfref_var_190(inds_var_1254 + 1, ig_var_1257) - selfref_var_190(inds_var_1254, ig_var_1257)))
          taufor_var_1286 = forfac_var_1245(jl_var_1315, lay_var_1259) * (forref_var_193(indf_var_1255, ig_var_1257) + forfrac_var_1244(jl_var_1315, lay_var_1259) * (forref_var_193(indf_var_1255 + 1, ig_var_1257) - forref_var_193(indf_var_1255, ig_var_1257)))
          co2m1_var_1290 = ka_mco2_var_191(jmco2_var_1262, indm_var_1256, ig_var_1257) + fmco2_var_1302 * (ka_mco2_var_191(jmco2_var_1262 + 1, indm_var_1256, ig_var_1257) - ka_mco2_var_191(jmco2_var_1262, indm_var_1256, ig_var_1257))
          co2m2_var_1291 = ka_mco2_var_191(jmco2_var_1262, indm_var_1256 + 1, ig_var_1257) + fmco2_var_1302 * (ka_mco2_var_191(jmco2_var_1262 + 1, indm_var_1256 + 1, ig_var_1257) - ka_mco2_var_191(jmco2_var_1262, indm_var_1256 + 1, ig_var_1257))
          absco2_var_1292 = co2m1_var_1290 + minorfrac_var_1246(jl_var_1315, lay_var_1259) * (co2m2_var_1291 - co2m1_var_1290)
          taug_var_1224(jl_var_1315, 76 + ig_var_1257, lay_var_1259) = tau_major_var_1288(ig_var_1257) + tau_major1_var_1289(ig_var_1257) + tauself_var_1287 + taufor_var_1286 + adjcolco2_var_1268 * absco2_var_1292
          fracs_var_1242(jl_var_1315, 76 + ig_var_1257, lay_var_1259) = fracrefa_var_186(ig_var_1257, jpl_var_1261) + fpl_var_1299 * (fracrefa_var_186(ig_var_1257, jpl_var_1261 + 1) - fracrefa_var_186(ig_var_1257, jpl_var_1261))
        END DO
      END DO
      ixc0_var_1312 = kfdia_var_1222 - kidia_var_1221 + 1 - ixc0_var_1312
      DO ixp_var_1313 = 1, ixc0_var_1312
        jl_var_1315 = ixhigh_var_1309(ixp_var_1313, lay_var_1259)
        chi_co2_var_1265 = colco2_var_1236(jl_var_1315, lay_var_1259) / (coldry_var_1237(jl_var_1315, lay_var_1259))
        ratco2_var_1266 = 1D+20 * chi_co2_var_1265 / chi_mls(2, jp_var_1230(jl_var_1315, lay_var_1259) + 1)
        IF (ratco2_var_1266 .GT. 3.0D0) THEN
          adjfac_var_1267 = 2.0D0 + (ratco2_var_1266 - 2.0D0) ** 0.79D0
          adjcolco2_var_1268 = adjfac_var_1267 * chi_mls(2, jp_var_1230(jl_var_1315, lay_var_1259) + 1) * coldry_var_1237(jl_var_1315, lay_var_1259) * 1D-20
        ELSE
          adjcolco2_var_1268 = colco2_var_1236(jl_var_1315, lay_var_1259)
        END IF
        ind0_var_1252 = ((jp_var_1230(jl_var_1315, lay_var_1259) - 13) * 5 + (jt_var_1231(jl_var_1315, lay_var_1259) - 1)) * nspb_var_217(7) + 1
        ind1_var_1253 = ((jp_var_1230(jl_var_1315, lay_var_1259) - 12) * 5 + (jt1_var_1232(jl_var_1315, lay_var_1259) - 1)) * nspb_var_217(7) + 1
        indm_var_1256 = indminor_var_1247(jl_var_1315, lay_var_1259)
        DO ig_var_1257 = 1, 12
          absco2_var_1292 = kb_mco2_var_192(indm_var_1256, ig_var_1257) + minorfrac_var_1246(jl_var_1315, lay_var_1259) * (kb_mco2_var_192(indm_var_1256 + 1, ig_var_1257) - kb_mco2_var_192(indm_var_1256, ig_var_1257))
          taug_var_1224(jl_var_1315, 76 + ig_var_1257, lay_var_1259) = colo3_var_1235(jl_var_1315, lay_var_1259) * (fac00_var_1226(jl_var_1315, lay_var_1259) * absb_var_189(ind0_var_1252, ig_var_1257) + fac10_var_1228(jl_var_1315, lay_var_1259) * absb_var_189(ind0_var_1252 + 1, ig_var_1257) + fac01_var_1227(jl_var_1315, lay_var_1259) * absb_var_189(ind1_var_1253, ig_var_1257) + fac11_var_1229(jl_var_1315, lay_var_1259) * absb_var_189(ind1_var_1253 + 1, ig_var_1257)) + adjcolco2_var_1268 * absco2_var_1292
          fracs_var_1242(jl_var_1315, 76 + ig_var_1257, lay_var_1259) = fracrefb_var_187(ig_var_1257)
        END DO
      END DO
      DO ixp_var_1313 = 1, ixc0_var_1312
        jl_var_1315 = ixhigh_var_1309(ixp_var_1313, lay_var_1259)
        taug_var_1224(jl_var_1315, 82, lay_var_1259) = taug_var_1224(jl_var_1315, 82, lay_var_1259) * 0.92D0
        taug_var_1224(jl_var_1315, 83, lay_var_1259) = taug_var_1224(jl_var_1315, 83, lay_var_1259) * 0.88D0
        taug_var_1224(jl_var_1315, 84, lay_var_1259) = taug_var_1224(jl_var_1315, 84, lay_var_1259) * 1.07D0
        taug_var_1224(jl_var_1315, 85, lay_var_1259) = taug_var_1224(jl_var_1315, 85, lay_var_1259) * 1.1D0
        taug_var_1224(jl_var_1315, 86, lay_var_1259) = taug_var_1224(jl_var_1315, 86, lay_var_1259) * 0.99D0
        taug_var_1224(jl_var_1315, 87, lay_var_1259) = taug_var_1224(jl_var_1315, 87, lay_var_1259) * 0.855D0
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol7
SUBROUTINE rrtm_taumol3(kidia_var_1316, kfdia_var_1317, klev_var_1318, taug_var_1319, p_tauaerl_var_1320, fac00_var_1321, fac01_var_1322, fac10_var_1323, fac11_var_1324, forfac_var_1325, forfrac_var_1342, indfor_var_1341, jp_var_1326, jt_var_1327, jt1_var_1328, oneminus_var_1329, colh2o_var_1330, colco2_var_1331, coln2o_var_1332, coldry_var_1333, laytrop_var_1334, selffac_var_1335, selffrac_var_1336, indself_var_1337, fracs_var_1338, rat_h2oco2_var_1339, rat_h2oco2_1_var_1340, minorfrac_var_1343, indminor_var_1344)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng3
  USE yoerrta3, ONLY: absa_var_163, absb_var_164, forref_var_166, fracrefa_var_159, fracrefb_var_160, ka_mn2o_var_161, kb_mn2o_var_162, selfref_var_165
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1316
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1317
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1318
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1319(kidia_var_1316 : kfdia_var_1317, 140, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1320(kidia_var_1316 : kfdia_var_1317, klev_var_1318, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1321(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1322(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1323(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1324(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1325(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1326(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1327(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1328(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1329
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1330(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1331(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_1332(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1333(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1334(kidia_var_1316 : kfdia_var_1317)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1335(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1336(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1337(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1338(kidia_var_1316 : kfdia_var_1317, 140, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1339(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1340(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1341(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1342(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1343(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1344(kidia_var_1316 : kfdia_var_1317, klev_var_1318)
  REAL(KIND = 8) :: speccomb_var_1345, speccomb1_var_1346, speccomb_mn2o_var_1347, speccomb_planck_var_1348
  REAL(KIND = 8) :: refrat_planck_a_var_1349, refrat_planck_b_var_1350, refrat_m_a_var_1351, refrat_m_b
  INTEGER(KIND = 4) :: ind0_var_1352, ind1_var_1353, inds_var_1354, indf_var_1355, indm_var_1356
  INTEGER(KIND = 4) :: ig_var_1357, js_var_1358, lay_var_1359, js1_var_1360, jmn2o_var_1361, jpl_var_1362
  REAL(KIND = 8) :: fs_var_1363, specmult_var_1364, specparm_var_1365, fs1_var_1366, specmult1_var_1367, specparm1_var_1368, fmn2o_var_1369, specmult_mn2o_var_1370, specparm_mn2o_var_1371, fpl_var_1372, specmult_planck_var_1373, specparm_planck_var_1374
  REAL(KIND = 8) :: adjfac_var_1375, adjcoln2o_var_1376, ratn2o_var_1377, chi_n2o_var_1378
  REAL(KIND = 8) :: fac000_var_1379, fac100_var_1380, fac200_var_1381, fac010_var_1382, fac110_var_1383, fac210_var_1384, fac001_var_1385, fac101_var_1386, fac201_var_1387, fac011_var_1388, fac111_var_1389, fac211_var_1390
  REAL(KIND = 8) :: p_var_1391, p4_var_1392, fk0_var_1393, fk1_var_1394, fk2_var_1395
  REAL(KIND = 8) :: taufor_var_1396, tauself_var_1397, n2om1_var_1398, n2om2_var_1399, absn2o_var_1400, tau_major_var_1401(16), tau_major1_var_1402(16)
  INTEGER(KIND = 4) :: laytrop_min_var_1403, laytrop_max_var_1404
  INTEGER(KIND = 4) :: ixc_var_1405(klev_var_1318), ixlow_var_1406(kfdia_var_1317, klev_var_1318), ixhigh_var_1407(kfdia_var_1317, klev_var_1318)
  INTEGER(KIND = 4) :: ich_var_1408, icl_var_1409, ixc0_var_1410, ixp_var_1411, jc_var_1412, jl_var_1413
  laytrop_min_var_1403 = MINVAL(laytrop_var_1334)
  laytrop_max_var_1404 = MAXVAL(laytrop_var_1334)
  ixlow_var_1406 = 0
  ixhigh_var_1407 = 0
  ixc_var_1405 = 0
  DO lay_var_1359 = laytrop_min_var_1403 + 1, laytrop_max_var_1404
    icl_var_1409 = 0
    ich_var_1408 = 0
    DO jc_var_1412 = kidia_var_1316, kfdia_var_1317
      IF (lay_var_1359 <= laytrop_var_1334(jc_var_1412)) THEN
        icl_var_1409 = icl_var_1409 + 1
        ixlow_var_1406(icl_var_1409, lay_var_1359) = jc_var_1412
      ELSE
        ich_var_1408 = ich_var_1408 + 1
        ixhigh_var_1407(ich_var_1408, lay_var_1359) = jc_var_1412
      END IF
    END DO
    ixc_var_1405(lay_var_1359) = icl_var_1409
  END DO
  refrat_planck_a_var_1349 = chi_mls(1, 9) / chi_mls(2, 9)
  refrat_planck_b_var_1350 = chi_mls(1, 13) / chi_mls(2, 13)
  refrat_m_a_var_1351 = chi_mls(1, 3) / chi_mls(2, 3)
  refrat_m_b = chi_mls(1, 13) / chi_mls(2, 13)
  DO lay_var_1359 = 1, laytrop_min_var_1403
    DO jl_var_1413 = kidia_var_1316, kfdia_var_1317
      speccomb_var_1345 = colh2o_var_1330(jl_var_1413, lay_var_1359) + rat_h2oco2_var_1339(jl_var_1413, lay_var_1359) * colco2_var_1331(jl_var_1413, lay_var_1359)
      specparm_var_1365 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_var_1345, 0.999999D0)
      specmult_var_1364 = 8.0D0 * (specparm_var_1365)
      js_var_1358 = 1 + INT(specmult_var_1364)
      fs_var_1363 = ((specmult_var_1364) - AINT((specmult_var_1364)))
      speccomb1_var_1346 = colh2o_var_1330(jl_var_1413, lay_var_1359) + rat_h2oco2_1_var_1340(jl_var_1413, lay_var_1359) * colco2_var_1331(jl_var_1413, lay_var_1359)
      specparm1_var_1368 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb1_var_1346, 0.999999D0)
      specmult1_var_1367 = 8.0D0 * (specparm1_var_1368)
      js1_var_1360 = 1 + INT(specmult1_var_1367)
      fs1_var_1366 = ((specmult1_var_1367) - AINT((specmult1_var_1367)))
      speccomb_mn2o_var_1347 = colh2o_var_1330(jl_var_1413, lay_var_1359) + refrat_m_a_var_1351 * colco2_var_1331(jl_var_1413, lay_var_1359)
      specparm_mn2o_var_1371 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_mn2o_var_1347, 0.999999D0)
      specmult_mn2o_var_1370 = 8.0D0 * specparm_mn2o_var_1371
      jmn2o_var_1361 = 1 + INT(specmult_mn2o_var_1370)
      fmn2o_var_1369 = ((specmult_mn2o_var_1370) - AINT((specmult_mn2o_var_1370)))
      chi_n2o_var_1378 = coln2o_var_1332(jl_var_1413, lay_var_1359) / coldry_var_1333(jl_var_1413, lay_var_1359)
      ratn2o_var_1377 = 1D+20 * chi_n2o_var_1378 / chi_mls(4, jp_var_1326(jl_var_1413, lay_var_1359) + 1)
      IF (ratn2o_var_1377 .GT. 1.5D0) THEN
        adjfac_var_1375 = 0.5D0 + (ratn2o_var_1377 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1376 = adjfac_var_1375 * chi_mls(4, jp_var_1326(jl_var_1413, lay_var_1359) + 1) * coldry_var_1333(jl_var_1413, lay_var_1359) * 1D-20
      ELSE
        adjcoln2o_var_1376 = coln2o_var_1332(jl_var_1413, lay_var_1359)
      END IF
      speccomb_planck_var_1348 = colh2o_var_1330(jl_var_1413, lay_var_1359) + refrat_planck_a_var_1349 * colco2_var_1331(jl_var_1413, lay_var_1359)
      specparm_planck_var_1374 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_planck_var_1348, 0.999999D0)
      specmult_planck_var_1373 = 8.0D0 * specparm_planck_var_1374
      jpl_var_1362 = 1 + INT(specmult_planck_var_1373)
      fpl_var_1372 = ((specmult_planck_var_1373) - AINT((specmult_planck_var_1373)))
      ind0_var_1352 = ((jp_var_1326(jl_var_1413, lay_var_1359) - 1) * 5 + (jt_var_1327(jl_var_1413, lay_var_1359) - 1)) * nspa_var_216(3) + js_var_1358
      ind1_var_1353 = (jp_var_1326(jl_var_1413, lay_var_1359) * 5 + (jt1_var_1328(jl_var_1413, lay_var_1359) - 1)) * nspa_var_216(3) + js1_var_1360
      inds_var_1354 = indself_var_1337(jl_var_1413, lay_var_1359)
      indf_var_1355 = indfor_var_1341(jl_var_1413, lay_var_1359)
      indm_var_1356 = indminor_var_1344(jl_var_1413, lay_var_1359)
      IF (specparm_var_1365 .LT. 0.125D0) THEN
        p_var_1391 = fs_var_1363 - 1.0D0
        p4_var_1392 = p_var_1391 ** 4
        fk0_var_1393 = p4_var_1392
        fk1_var_1394 = 1.0D0 - p_var_1391 - 2.0D0 * p4_var_1392
        fk2_var_1395 = p_var_1391 + p4_var_1392
        fac000_var_1379 = fk0_var_1393 * fac00_var_1321(jl_var_1413, lay_var_1359)
        fac100_var_1380 = fk1_var_1394 * fac00_var_1321(jl_var_1413, lay_var_1359)
        fac200_var_1381 = fk2_var_1395 * fac00_var_1321(jl_var_1413, lay_var_1359)
        fac010_var_1382 = fk0_var_1393 * fac10_var_1323(jl_var_1413, lay_var_1359)
        fac110_var_1383 = fk1_var_1394 * fac10_var_1323(jl_var_1413, lay_var_1359)
        fac210_var_1384 = fk2_var_1395 * fac10_var_1323(jl_var_1413, lay_var_1359)
      ELSE IF (specparm_var_1365 .GT. 0.875D0) THEN
        p_var_1391 = - fs_var_1363
        p4_var_1392 = p_var_1391 ** 4
        fk0_var_1393 = p4_var_1392
        fk1_var_1394 = 1.0D0 - p_var_1391 - 2.0D0 * p4_var_1392
        fk2_var_1395 = p_var_1391 + p4_var_1392
        fac000_var_1379 = fk0_var_1393 * fac00_var_1321(jl_var_1413, lay_var_1359)
        fac100_var_1380 = fk1_var_1394 * fac00_var_1321(jl_var_1413, lay_var_1359)
        fac200_var_1381 = fk2_var_1395 * fac00_var_1321(jl_var_1413, lay_var_1359)
        fac010_var_1382 = fk0_var_1393 * fac10_var_1323(jl_var_1413, lay_var_1359)
        fac110_var_1383 = fk1_var_1394 * fac10_var_1323(jl_var_1413, lay_var_1359)
        fac210_var_1384 = fk2_var_1395 * fac10_var_1323(jl_var_1413, lay_var_1359)
      ELSE
        fac000_var_1379 = (1.0D0 - fs_var_1363) * fac00_var_1321(jl_var_1413, lay_var_1359)
        fac010_var_1382 = (1.0D0 - fs_var_1363) * fac10_var_1323(jl_var_1413, lay_var_1359)
        fac100_var_1380 = fs_var_1363 * fac00_var_1321(jl_var_1413, lay_var_1359)
        fac110_var_1383 = fs_var_1363 * fac10_var_1323(jl_var_1413, lay_var_1359)
        fac200_var_1381 = 0.0D0
        fac210_var_1384 = 0.0D0
      END IF
      IF (specparm1_var_1368 .LT. 0.125D0) THEN
        p_var_1391 = fs1_var_1366 - 1.0D0
        p4_var_1392 = p_var_1391 ** 4
        fk0_var_1393 = p4_var_1392
        fk1_var_1394 = 1.0D0 - p_var_1391 - 2.0D0 * p4_var_1392
        fk2_var_1395 = p_var_1391 + p4_var_1392
        fac001_var_1385 = fk0_var_1393 * fac01_var_1322(jl_var_1413, lay_var_1359)
        fac101_var_1386 = fk1_var_1394 * fac01_var_1322(jl_var_1413, lay_var_1359)
        fac201_var_1387 = fk2_var_1395 * fac01_var_1322(jl_var_1413, lay_var_1359)
        fac011_var_1388 = fk0_var_1393 * fac11_var_1324(jl_var_1413, lay_var_1359)
        fac111_var_1389 = fk1_var_1394 * fac11_var_1324(jl_var_1413, lay_var_1359)
        fac211_var_1390 = fk2_var_1395 * fac11_var_1324(jl_var_1413, lay_var_1359)
      ELSE IF (specparm1_var_1368 .GT. 0.875D0) THEN
        p_var_1391 = - fs1_var_1366
        p4_var_1392 = p_var_1391 ** 4
        fk0_var_1393 = p4_var_1392
        fk1_var_1394 = 1.0D0 - p_var_1391 - 2.0D0 * p4_var_1392
        fk2_var_1395 = p_var_1391 + p4_var_1392
        fac001_var_1385 = fk0_var_1393 * fac01_var_1322(jl_var_1413, lay_var_1359)
        fac101_var_1386 = fk1_var_1394 * fac01_var_1322(jl_var_1413, lay_var_1359)
        fac201_var_1387 = fk2_var_1395 * fac01_var_1322(jl_var_1413, lay_var_1359)
        fac011_var_1388 = fk0_var_1393 * fac11_var_1324(jl_var_1413, lay_var_1359)
        fac111_var_1389 = fk1_var_1394 * fac11_var_1324(jl_var_1413, lay_var_1359)
        fac211_var_1390 = fk2_var_1395 * fac11_var_1324(jl_var_1413, lay_var_1359)
      ELSE
        fac001_var_1385 = (1.0D0 - fs1_var_1366) * fac01_var_1322(jl_var_1413, lay_var_1359)
        fac011_var_1388 = (1.0D0 - fs1_var_1366) * fac11_var_1324(jl_var_1413, lay_var_1359)
        fac101_var_1386 = fs1_var_1366 * fac01_var_1322(jl_var_1413, lay_var_1359)
        fac111_var_1389 = fs1_var_1366 * fac11_var_1324(jl_var_1413, lay_var_1359)
        fac201_var_1387 = 0.0D0
        fac211_var_1390 = 0.0D0
      END IF
      IF (specparm_var_1365 .LT. 0.125D0) THEN
        tau_major_var_1401(1 : ng3) = speccomb_var_1345 * (fac000_var_1379 * absa_var_163(ind0_var_1352, 1 : 16) + fac100_var_1380 * absa_var_163(ind0_var_1352 + 1, 1 : 16) + fac200_var_1381 * absa_var_163(ind0_var_1352 + 2, 1 : 16) + fac010_var_1382 * absa_var_163(ind0_var_1352 + 9, 1 : 16) + fac110_var_1383 * absa_var_163(ind0_var_1352 + 10, 1 : 16) + fac210_var_1384 * absa_var_163(ind0_var_1352 + 11, 1 : 16))
      ELSE IF (specparm_var_1365 .GT. 0.875D0) THEN
        tau_major_var_1401(1 : ng3) = speccomb_var_1345 * (fac200_var_1381 * absa_var_163(ind0_var_1352 - 1, 1 : 16) + fac100_var_1380 * absa_var_163(ind0_var_1352, 1 : 16) + fac000_var_1379 * absa_var_163(ind0_var_1352 + 1, 1 : 16) + fac210_var_1384 * absa_var_163(ind0_var_1352 + 8, 1 : 16) + fac110_var_1383 * absa_var_163(ind0_var_1352 + 9, 1 : 16) + fac010_var_1382 * absa_var_163(ind0_var_1352 + 10, 1 : 16))
      ELSE
        tau_major_var_1401(1 : ng3) = speccomb_var_1345 * (fac000_var_1379 * absa_var_163(ind0_var_1352, 1 : 16) + fac100_var_1380 * absa_var_163(ind0_var_1352 + 1, 1 : 16) + fac010_var_1382 * absa_var_163(ind0_var_1352 + 9, 1 : 16) + fac110_var_1383 * absa_var_163(ind0_var_1352 + 10, 1 : 16))
      END IF
      IF (specparm1_var_1368 .LT. 0.125D0) THEN
        tau_major1_var_1402(1 : ng3) = speccomb1_var_1346 * (fac001_var_1385 * absa_var_163(ind1_var_1353, 1 : 16) + fac101_var_1386 * absa_var_163(ind1_var_1353 + 1, 1 : 16) + fac201_var_1387 * absa_var_163(ind1_var_1353 + 2, 1 : 16) + fac011_var_1388 * absa_var_163(ind1_var_1353 + 9, 1 : 16) + fac111_var_1389 * absa_var_163(ind1_var_1353 + 10, 1 : 16) + fac211_var_1390 * absa_var_163(ind1_var_1353 + 11, 1 : 16))
      ELSE IF (specparm1_var_1368 .GT. 0.875D0) THEN
        tau_major1_var_1402(1 : ng3) = speccomb1_var_1346 * (fac201_var_1387 * absa_var_163(ind1_var_1353 - 1, 1 : 16) + fac101_var_1386 * absa_var_163(ind1_var_1353, 1 : 16) + fac001_var_1385 * absa_var_163(ind1_var_1353 + 1, 1 : 16) + fac211_var_1390 * absa_var_163(ind1_var_1353 + 8, 1 : 16) + fac111_var_1389 * absa_var_163(ind1_var_1353 + 9, 1 : 16) + fac011_var_1388 * absa_var_163(ind1_var_1353 + 10, 1 : 16))
      ELSE
        tau_major1_var_1402(1 : ng3) = speccomb1_var_1346 * (fac001_var_1385 * absa_var_163(ind1_var_1353, 1 : 16) + fac101_var_1386 * absa_var_163(ind1_var_1353 + 1, 1 : 16) + fac011_var_1388 * absa_var_163(ind1_var_1353 + 9, 1 : 16) + fac111_var_1389 * absa_var_163(ind1_var_1353 + 10, 1 : 16))
      END IF
      DO ig_var_1357 = 1, 16
        tauself_var_1397 = selffac_var_1335(jl_var_1413, lay_var_1359) * (selfref_var_165(inds_var_1354, ig_var_1357) + selffrac_var_1336(jl_var_1413, lay_var_1359) * (selfref_var_165(inds_var_1354 + 1, ig_var_1357) - selfref_var_165(inds_var_1354, ig_var_1357)))
        taufor_var_1396 = forfac_var_1325(jl_var_1413, lay_var_1359) * (forref_var_166(indf_var_1355, ig_var_1357) + forfrac_var_1342(jl_var_1413, lay_var_1359) * (forref_var_166(indf_var_1355 + 1, ig_var_1357) - forref_var_166(indf_var_1355, ig_var_1357)))
        n2om1_var_1398 = ka_mn2o_var_161(jmn2o_var_1361, indm_var_1356, ig_var_1357) + fmn2o_var_1369 * (ka_mn2o_var_161(jmn2o_var_1361 + 1, indm_var_1356, ig_var_1357) - ka_mn2o_var_161(jmn2o_var_1361, indm_var_1356, ig_var_1357))
        n2om2_var_1399 = ka_mn2o_var_161(jmn2o_var_1361, indm_var_1356 + 1, ig_var_1357) + fmn2o_var_1369 * (ka_mn2o_var_161(jmn2o_var_1361 + 1, indm_var_1356 + 1, ig_var_1357) - ka_mn2o_var_161(jmn2o_var_1361, indm_var_1356 + 1, ig_var_1357))
        absn2o_var_1400 = n2om1_var_1398 + minorfrac_var_1343(jl_var_1413, lay_var_1359) * (n2om2_var_1399 - n2om1_var_1398)
        taug_var_1319(jl_var_1413, 22 + ig_var_1357, lay_var_1359) = tau_major_var_1401(ig_var_1357) + tau_major1_var_1402(ig_var_1357) + tauself_var_1397 + taufor_var_1396 + adjcoln2o_var_1376 * absn2o_var_1400
        fracs_var_1338(jl_var_1413, 22 + ig_var_1357, lay_var_1359) = fracrefa_var_159(ig_var_1357, jpl_var_1362) + fpl_var_1372 * (fracrefa_var_159(ig_var_1357, jpl_var_1362 + 1) - fracrefa_var_159(ig_var_1357, jpl_var_1362))
      END DO
    END DO
  END DO
  DO lay_var_1359 = laytrop_max_var_1404 + 1, klev_var_1318
    DO jl_var_1413 = kidia_var_1316, kfdia_var_1317
      speccomb_var_1345 = colh2o_var_1330(jl_var_1413, lay_var_1359) + rat_h2oco2_var_1339(jl_var_1413, lay_var_1359) * colco2_var_1331(jl_var_1413, lay_var_1359)
      specparm_var_1365 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_var_1345, 0.999999D0)
      specmult_var_1364 = 4.0D0 * (specparm_var_1365)
      js_var_1358 = 1 + INT(specmult_var_1364)
      fs_var_1363 = ((specmult_var_1364) - AINT((specmult_var_1364)))
      speccomb1_var_1346 = colh2o_var_1330(jl_var_1413, lay_var_1359) + rat_h2oco2_1_var_1340(jl_var_1413, lay_var_1359) * colco2_var_1331(jl_var_1413, lay_var_1359)
      specparm1_var_1368 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb1_var_1346, 0.999999D0)
      specmult1_var_1367 = 4.0D0 * (specparm1_var_1368)
      js1_var_1360 = 1 + INT(specmult1_var_1367)
      fs1_var_1366 = ((specmult1_var_1367) - AINT((specmult1_var_1367)))
      fac000_var_1379 = (1.0D0 - fs_var_1363) * fac00_var_1321(jl_var_1413, lay_var_1359)
      fac010_var_1382 = (1.0D0 - fs_var_1363) * fac10_var_1323(jl_var_1413, lay_var_1359)
      fac100_var_1380 = fs_var_1363 * fac00_var_1321(jl_var_1413, lay_var_1359)
      fac110_var_1383 = fs_var_1363 * fac10_var_1323(jl_var_1413, lay_var_1359)
      fac001_var_1385 = (1.0D0 - fs1_var_1366) * fac01_var_1322(jl_var_1413, lay_var_1359)
      fac011_var_1388 = (1.0D0 - fs1_var_1366) * fac11_var_1324(jl_var_1413, lay_var_1359)
      fac101_var_1386 = fs1_var_1366 * fac01_var_1322(jl_var_1413, lay_var_1359)
      fac111_var_1389 = fs1_var_1366 * fac11_var_1324(jl_var_1413, lay_var_1359)
      speccomb_mn2o_var_1347 = colh2o_var_1330(jl_var_1413, lay_var_1359) + refrat_m_b * colco2_var_1331(jl_var_1413, lay_var_1359)
      specparm_mn2o_var_1371 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_mn2o_var_1347, 0.999999D0)
      specmult_mn2o_var_1370 = 4.0D0 * specparm_mn2o_var_1371
      jmn2o_var_1361 = 1 + INT(specmult_mn2o_var_1370)
      fmn2o_var_1369 = ((specmult_mn2o_var_1370) - AINT((specmult_mn2o_var_1370)))
      chi_n2o_var_1378 = coln2o_var_1332(jl_var_1413, lay_var_1359) / coldry_var_1333(jl_var_1413, lay_var_1359)
      ratn2o_var_1377 = 1D+20 * chi_n2o_var_1378 / chi_mls(4, jp_var_1326(jl_var_1413, lay_var_1359) + 1)
      IF (ratn2o_var_1377 .GT. 1.5D0) THEN
        adjfac_var_1375 = 0.5D0 + (ratn2o_var_1377 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1376 = adjfac_var_1375 * chi_mls(4, jp_var_1326(jl_var_1413, lay_var_1359) + 1) * coldry_var_1333(jl_var_1413, lay_var_1359) * 1D-20
      ELSE
        adjcoln2o_var_1376 = coln2o_var_1332(jl_var_1413, lay_var_1359)
      END IF
      speccomb_planck_var_1348 = colh2o_var_1330(jl_var_1413, lay_var_1359) + refrat_planck_b_var_1350 * colco2_var_1331(jl_var_1413, lay_var_1359)
      specparm_planck_var_1374 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_planck_var_1348, 0.999999D0)
      specmult_planck_var_1373 = 4.0D0 * specparm_planck_var_1374
      jpl_var_1362 = 1 + INT(specmult_planck_var_1373)
      fpl_var_1372 = ((specmult_planck_var_1373) - AINT((specmult_planck_var_1373)))
      ind0_var_1352 = ((jp_var_1326(jl_var_1413, lay_var_1359) - 13) * 5 + (jt_var_1327(jl_var_1413, lay_var_1359) - 1)) * nspb_var_217(3) + js_var_1358
      ind1_var_1353 = ((jp_var_1326(jl_var_1413, lay_var_1359) - 12) * 5 + (jt1_var_1328(jl_var_1413, lay_var_1359) - 1)) * nspb_var_217(3) + js1_var_1360
      indf_var_1355 = indfor_var_1341(jl_var_1413, lay_var_1359)
      indm_var_1356 = indminor_var_1344(jl_var_1413, lay_var_1359)
      DO ig_var_1357 = 1, 16
        taufor_var_1396 = forfac_var_1325(jl_var_1413, lay_var_1359) * (forref_var_166(indf_var_1355, ig_var_1357) + forfrac_var_1342(jl_var_1413, lay_var_1359) * (forref_var_166(indf_var_1355 + 1, ig_var_1357) - forref_var_166(indf_var_1355, ig_var_1357)))
        n2om1_var_1398 = kb_mn2o_var_162(jmn2o_var_1361, indm_var_1356, ig_var_1357) + fmn2o_var_1369 * (kb_mn2o_var_162(jmn2o_var_1361 + 1, indm_var_1356, ig_var_1357) - kb_mn2o_var_162(jmn2o_var_1361, indm_var_1356, ig_var_1357))
        n2om2_var_1399 = kb_mn2o_var_162(jmn2o_var_1361, indm_var_1356 + 1, ig_var_1357) + fmn2o_var_1369 * (kb_mn2o_var_162(jmn2o_var_1361 + 1, indm_var_1356 + 1, ig_var_1357) - kb_mn2o_var_162(jmn2o_var_1361, indm_var_1356 + 1, ig_var_1357))
        absn2o_var_1400 = n2om1_var_1398 + minorfrac_var_1343(jl_var_1413, lay_var_1359) * (n2om2_var_1399 - n2om1_var_1398)
        taug_var_1319(jl_var_1413, 22 + ig_var_1357, lay_var_1359) = speccomb_var_1345 * (fac000_var_1379 * absb_var_164(ind0_var_1352, ig_var_1357) + fac100_var_1380 * absb_var_164(ind0_var_1352 + 1, ig_var_1357) + fac010_var_1382 * absb_var_164(ind0_var_1352 + 5, ig_var_1357) + fac110_var_1383 * absb_var_164(ind0_var_1352 + 6, ig_var_1357)) + speccomb1_var_1346 * (fac001_var_1385 * absb_var_164(ind1_var_1353, ig_var_1357) + fac101_var_1386 * absb_var_164(ind1_var_1353 + 1, ig_var_1357) + fac011_var_1388 * absb_var_164(ind1_var_1353 + 5, ig_var_1357) + fac111_var_1389 * absb_var_164(ind1_var_1353 + 6, ig_var_1357)) + taufor_var_1396 + adjcoln2o_var_1376 * absn2o_var_1400
        fracs_var_1338(jl_var_1413, 22 + ig_var_1357, lay_var_1359) = fracrefb_var_160(ig_var_1357, jpl_var_1362) + fpl_var_1372 * (fracrefb_var_160(ig_var_1357, jpl_var_1362 + 1) - fracrefb_var_160(ig_var_1357, jpl_var_1362))
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1404 /= laytrop_min_var_1403) THEN
    DO lay_var_1359 = laytrop_min_var_1403 + 1, laytrop_max_var_1404
      ixc0_var_1410 = ixc_var_1405(lay_var_1359)
      DO ixp_var_1411 = 1, ixc0_var_1410
        jl_var_1413 = ixlow_var_1406(ixp_var_1411, lay_var_1359)
        speccomb_var_1345 = colh2o_var_1330(jl_var_1413, lay_var_1359) + rat_h2oco2_var_1339(jl_var_1413, lay_var_1359) * colco2_var_1331(jl_var_1413, lay_var_1359)
        specparm_var_1365 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_var_1345, 0.999999D0)
        specmult_var_1364 = 8.0D0 * (specparm_var_1365)
        js_var_1358 = 1 + INT(specmult_var_1364)
        fs_var_1363 = ((specmult_var_1364) - AINT((specmult_var_1364)))
        speccomb1_var_1346 = colh2o_var_1330(jl_var_1413, lay_var_1359) + rat_h2oco2_1_var_1340(jl_var_1413, lay_var_1359) * colco2_var_1331(jl_var_1413, lay_var_1359)
        specparm1_var_1368 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb1_var_1346, 0.999999D0)
        specmult1_var_1367 = 8.0D0 * (specparm1_var_1368)
        js1_var_1360 = 1 + INT(specmult1_var_1367)
        fs1_var_1366 = ((specmult1_var_1367) - AINT((specmult1_var_1367)))
        speccomb_mn2o_var_1347 = colh2o_var_1330(jl_var_1413, lay_var_1359) + refrat_m_a_var_1351 * colco2_var_1331(jl_var_1413, lay_var_1359)
        specparm_mn2o_var_1371 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_mn2o_var_1347, 0.999999D0)
        specmult_mn2o_var_1370 = 8.0D0 * specparm_mn2o_var_1371
        jmn2o_var_1361 = 1 + INT(specmult_mn2o_var_1370)
        fmn2o_var_1369 = ((specmult_mn2o_var_1370) - AINT((specmult_mn2o_var_1370)))
        chi_n2o_var_1378 = coln2o_var_1332(jl_var_1413, lay_var_1359) / coldry_var_1333(jl_var_1413, lay_var_1359)
        ratn2o_var_1377 = 1D+20 * chi_n2o_var_1378 / chi_mls(4, jp_var_1326(jl_var_1413, lay_var_1359) + 1)
        IF (ratn2o_var_1377 .GT. 1.5D0) THEN
          adjfac_var_1375 = 0.5D0 + (ratn2o_var_1377 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1376 = adjfac_var_1375 * chi_mls(4, jp_var_1326(jl_var_1413, lay_var_1359) + 1) * coldry_var_1333(jl_var_1413, lay_var_1359) * 1D-20
        ELSE
          adjcoln2o_var_1376 = coln2o_var_1332(jl_var_1413, lay_var_1359)
        END IF
        speccomb_planck_var_1348 = colh2o_var_1330(jl_var_1413, lay_var_1359) + refrat_planck_a_var_1349 * colco2_var_1331(jl_var_1413, lay_var_1359)
        specparm_planck_var_1374 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_planck_var_1348, 0.999999D0)
        specmult_planck_var_1373 = 8.0D0 * specparm_planck_var_1374
        jpl_var_1362 = 1 + INT(specmult_planck_var_1373)
        fpl_var_1372 = ((specmult_planck_var_1373) - AINT((specmult_planck_var_1373)))
        ind0_var_1352 = ((jp_var_1326(jl_var_1413, lay_var_1359) - 1) * 5 + (jt_var_1327(jl_var_1413, lay_var_1359) - 1)) * nspa_var_216(3) + js_var_1358
        ind1_var_1353 = (jp_var_1326(jl_var_1413, lay_var_1359) * 5 + (jt1_var_1328(jl_var_1413, lay_var_1359) - 1)) * nspa_var_216(3) + js1_var_1360
        inds_var_1354 = indself_var_1337(jl_var_1413, lay_var_1359)
        indf_var_1355 = indfor_var_1341(jl_var_1413, lay_var_1359)
        indm_var_1356 = indminor_var_1344(jl_var_1413, lay_var_1359)
        IF (specparm_var_1365 .LT. 0.125D0) THEN
          p_var_1391 = fs_var_1363 - 1.0D0
          p4_var_1392 = p_var_1391 ** 4
          fk0_var_1393 = p4_var_1392
          fk1_var_1394 = 1.0D0 - p_var_1391 - 2.0D0 * p4_var_1392
          fk2_var_1395 = p_var_1391 + p4_var_1392
          fac000_var_1379 = fk0_var_1393 * fac00_var_1321(jl_var_1413, lay_var_1359)
          fac100_var_1380 = fk1_var_1394 * fac00_var_1321(jl_var_1413, lay_var_1359)
          fac200_var_1381 = fk2_var_1395 * fac00_var_1321(jl_var_1413, lay_var_1359)
          fac010_var_1382 = fk0_var_1393 * fac10_var_1323(jl_var_1413, lay_var_1359)
          fac110_var_1383 = fk1_var_1394 * fac10_var_1323(jl_var_1413, lay_var_1359)
          fac210_var_1384 = fk2_var_1395 * fac10_var_1323(jl_var_1413, lay_var_1359)
        ELSE IF (specparm_var_1365 .GT. 0.875D0) THEN
          p_var_1391 = - fs_var_1363
          p4_var_1392 = p_var_1391 ** 4
          fk0_var_1393 = p4_var_1392
          fk1_var_1394 = 1.0D0 - p_var_1391 - 2.0D0 * p4_var_1392
          fk2_var_1395 = p_var_1391 + p4_var_1392
          fac000_var_1379 = fk0_var_1393 * fac00_var_1321(jl_var_1413, lay_var_1359)
          fac100_var_1380 = fk1_var_1394 * fac00_var_1321(jl_var_1413, lay_var_1359)
          fac200_var_1381 = fk2_var_1395 * fac00_var_1321(jl_var_1413, lay_var_1359)
          fac010_var_1382 = fk0_var_1393 * fac10_var_1323(jl_var_1413, lay_var_1359)
          fac110_var_1383 = fk1_var_1394 * fac10_var_1323(jl_var_1413, lay_var_1359)
          fac210_var_1384 = fk2_var_1395 * fac10_var_1323(jl_var_1413, lay_var_1359)
        ELSE
          fac000_var_1379 = (1.0D0 - fs_var_1363) * fac00_var_1321(jl_var_1413, lay_var_1359)
          fac010_var_1382 = (1.0D0 - fs_var_1363) * fac10_var_1323(jl_var_1413, lay_var_1359)
          fac100_var_1380 = fs_var_1363 * fac00_var_1321(jl_var_1413, lay_var_1359)
          fac110_var_1383 = fs_var_1363 * fac10_var_1323(jl_var_1413, lay_var_1359)
          fac200_var_1381 = 0.0D0
          fac210_var_1384 = 0.0D0
        END IF
        IF (specparm1_var_1368 .LT. 0.125D0) THEN
          p_var_1391 = fs1_var_1366 - 1.0D0
          p4_var_1392 = p_var_1391 ** 4
          fk0_var_1393 = p4_var_1392
          fk1_var_1394 = 1.0D0 - p_var_1391 - 2.0D0 * p4_var_1392
          fk2_var_1395 = p_var_1391 + p4_var_1392
          fac001_var_1385 = fk0_var_1393 * fac01_var_1322(jl_var_1413, lay_var_1359)
          fac101_var_1386 = fk1_var_1394 * fac01_var_1322(jl_var_1413, lay_var_1359)
          fac201_var_1387 = fk2_var_1395 * fac01_var_1322(jl_var_1413, lay_var_1359)
          fac011_var_1388 = fk0_var_1393 * fac11_var_1324(jl_var_1413, lay_var_1359)
          fac111_var_1389 = fk1_var_1394 * fac11_var_1324(jl_var_1413, lay_var_1359)
          fac211_var_1390 = fk2_var_1395 * fac11_var_1324(jl_var_1413, lay_var_1359)
        ELSE IF (specparm1_var_1368 .GT. 0.875D0) THEN
          p_var_1391 = - fs1_var_1366
          p4_var_1392 = p_var_1391 ** 4
          fk0_var_1393 = p4_var_1392
          fk1_var_1394 = 1.0D0 - p_var_1391 - 2.0D0 * p4_var_1392
          fk2_var_1395 = p_var_1391 + p4_var_1392
          fac001_var_1385 = fk0_var_1393 * fac01_var_1322(jl_var_1413, lay_var_1359)
          fac101_var_1386 = fk1_var_1394 * fac01_var_1322(jl_var_1413, lay_var_1359)
          fac201_var_1387 = fk2_var_1395 * fac01_var_1322(jl_var_1413, lay_var_1359)
          fac011_var_1388 = fk0_var_1393 * fac11_var_1324(jl_var_1413, lay_var_1359)
          fac111_var_1389 = fk1_var_1394 * fac11_var_1324(jl_var_1413, lay_var_1359)
          fac211_var_1390 = fk2_var_1395 * fac11_var_1324(jl_var_1413, lay_var_1359)
        ELSE
          fac001_var_1385 = (1.0D0 - fs1_var_1366) * fac01_var_1322(jl_var_1413, lay_var_1359)
          fac011_var_1388 = (1.0D0 - fs1_var_1366) * fac11_var_1324(jl_var_1413, lay_var_1359)
          fac101_var_1386 = fs1_var_1366 * fac01_var_1322(jl_var_1413, lay_var_1359)
          fac111_var_1389 = fs1_var_1366 * fac11_var_1324(jl_var_1413, lay_var_1359)
          fac201_var_1387 = 0.0D0
          fac211_var_1390 = 0.0D0
        END IF
        IF (specparm_var_1365 .LT. 0.125D0) THEN
          tau_major_var_1401(1 : ng3) = speccomb_var_1345 * (fac000_var_1379 * absa_var_163(ind0_var_1352, 1 : 16) + fac100_var_1380 * absa_var_163(ind0_var_1352 + 1, 1 : 16) + fac200_var_1381 * absa_var_163(ind0_var_1352 + 2, 1 : 16) + fac010_var_1382 * absa_var_163(ind0_var_1352 + 9, 1 : 16) + fac110_var_1383 * absa_var_163(ind0_var_1352 + 10, 1 : 16) + fac210_var_1384 * absa_var_163(ind0_var_1352 + 11, 1 : 16))
        ELSE IF (specparm_var_1365 .GT. 0.875D0) THEN
          tau_major_var_1401(1 : ng3) = speccomb_var_1345 * (fac200_var_1381 * absa_var_163(ind0_var_1352 - 1, 1 : 16) + fac100_var_1380 * absa_var_163(ind0_var_1352, 1 : 16) + fac000_var_1379 * absa_var_163(ind0_var_1352 + 1, 1 : 16) + fac210_var_1384 * absa_var_163(ind0_var_1352 + 8, 1 : 16) + fac110_var_1383 * absa_var_163(ind0_var_1352 + 9, 1 : 16) + fac010_var_1382 * absa_var_163(ind0_var_1352 + 10, 1 : 16))
        ELSE
          tau_major_var_1401(1 : ng3) = speccomb_var_1345 * (fac000_var_1379 * absa_var_163(ind0_var_1352, 1 : 16) + fac100_var_1380 * absa_var_163(ind0_var_1352 + 1, 1 : 16) + fac010_var_1382 * absa_var_163(ind0_var_1352 + 9, 1 : 16) + fac110_var_1383 * absa_var_163(ind0_var_1352 + 10, 1 : 16))
        END IF
        IF (specparm1_var_1368 .LT. 0.125D0) THEN
          tau_major1_var_1402(1 : ng3) = speccomb1_var_1346 * (fac001_var_1385 * absa_var_163(ind1_var_1353, 1 : 16) + fac101_var_1386 * absa_var_163(ind1_var_1353 + 1, 1 : 16) + fac201_var_1387 * absa_var_163(ind1_var_1353 + 2, 1 : 16) + fac011_var_1388 * absa_var_163(ind1_var_1353 + 9, 1 : 16) + fac111_var_1389 * absa_var_163(ind1_var_1353 + 10, 1 : 16) + fac211_var_1390 * absa_var_163(ind1_var_1353 + 11, 1 : 16))
        ELSE IF (specparm1_var_1368 .GT. 0.875D0) THEN
          tau_major1_var_1402(1 : ng3) = speccomb1_var_1346 * (fac201_var_1387 * absa_var_163(ind1_var_1353 - 1, 1 : 16) + fac101_var_1386 * absa_var_163(ind1_var_1353, 1 : 16) + fac001_var_1385 * absa_var_163(ind1_var_1353 + 1, 1 : 16) + fac211_var_1390 * absa_var_163(ind1_var_1353 + 8, 1 : 16) + fac111_var_1389 * absa_var_163(ind1_var_1353 + 9, 1 : 16) + fac011_var_1388 * absa_var_163(ind1_var_1353 + 10, 1 : 16))
        ELSE
          tau_major1_var_1402(1 : ng3) = speccomb1_var_1346 * (fac001_var_1385 * absa_var_163(ind1_var_1353, 1 : 16) + fac101_var_1386 * absa_var_163(ind1_var_1353 + 1, 1 : 16) + fac011_var_1388 * absa_var_163(ind1_var_1353 + 9, 1 : 16) + fac111_var_1389 * absa_var_163(ind1_var_1353 + 10, 1 : 16))
        END IF
        DO ig_var_1357 = 1, 16
          tauself_var_1397 = selffac_var_1335(jl_var_1413, lay_var_1359) * (selfref_var_165(inds_var_1354, ig_var_1357) + selffrac_var_1336(jl_var_1413, lay_var_1359) * (selfref_var_165(inds_var_1354 + 1, ig_var_1357) - selfref_var_165(inds_var_1354, ig_var_1357)))
          taufor_var_1396 = forfac_var_1325(jl_var_1413, lay_var_1359) * (forref_var_166(indf_var_1355, ig_var_1357) + forfrac_var_1342(jl_var_1413, lay_var_1359) * (forref_var_166(indf_var_1355 + 1, ig_var_1357) - forref_var_166(indf_var_1355, ig_var_1357)))
          n2om1_var_1398 = ka_mn2o_var_161(jmn2o_var_1361, indm_var_1356, ig_var_1357) + fmn2o_var_1369 * (ka_mn2o_var_161(jmn2o_var_1361 + 1, indm_var_1356, ig_var_1357) - ka_mn2o_var_161(jmn2o_var_1361, indm_var_1356, ig_var_1357))
          n2om2_var_1399 = ka_mn2o_var_161(jmn2o_var_1361, indm_var_1356 + 1, ig_var_1357) + fmn2o_var_1369 * (ka_mn2o_var_161(jmn2o_var_1361 + 1, indm_var_1356 + 1, ig_var_1357) - ka_mn2o_var_161(jmn2o_var_1361, indm_var_1356 + 1, ig_var_1357))
          absn2o_var_1400 = n2om1_var_1398 + minorfrac_var_1343(jl_var_1413, lay_var_1359) * (n2om2_var_1399 - n2om1_var_1398)
          taug_var_1319(jl_var_1413, 22 + ig_var_1357, lay_var_1359) = tau_major_var_1401(ig_var_1357) + tau_major1_var_1402(ig_var_1357) + tauself_var_1397 + taufor_var_1396 + adjcoln2o_var_1376 * absn2o_var_1400
          fracs_var_1338(jl_var_1413, 22 + ig_var_1357, lay_var_1359) = fracrefa_var_159(ig_var_1357, jpl_var_1362) + fpl_var_1372 * (fracrefa_var_159(ig_var_1357, jpl_var_1362 + 1) - fracrefa_var_159(ig_var_1357, jpl_var_1362))
        END DO
      END DO
      ixc0_var_1410 = kfdia_var_1317 - kidia_var_1316 + 1 - ixc0_var_1410
      DO ixp_var_1411 = 1, ixc0_var_1410
        jl_var_1413 = ixhigh_var_1407(ixp_var_1411, lay_var_1359)
        speccomb_var_1345 = colh2o_var_1330(jl_var_1413, lay_var_1359) + rat_h2oco2_var_1339(jl_var_1413, lay_var_1359) * colco2_var_1331(jl_var_1413, lay_var_1359)
        specparm_var_1365 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_var_1345, 0.999999D0)
        specmult_var_1364 = 4.0D0 * (specparm_var_1365)
        js_var_1358 = 1 + INT(specmult_var_1364)
        fs_var_1363 = ((specmult_var_1364) - AINT((specmult_var_1364)))
        speccomb1_var_1346 = colh2o_var_1330(jl_var_1413, lay_var_1359) + rat_h2oco2_1_var_1340(jl_var_1413, lay_var_1359) * colco2_var_1331(jl_var_1413, lay_var_1359)
        specparm1_var_1368 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb1_var_1346, 0.999999D0)
        specmult1_var_1367 = 4.0D0 * (specparm1_var_1368)
        js1_var_1360 = 1 + INT(specmult1_var_1367)
        fs1_var_1366 = ((specmult1_var_1367) - AINT((specmult1_var_1367)))
        fac000_var_1379 = (1.0D0 - fs_var_1363) * fac00_var_1321(jl_var_1413, lay_var_1359)
        fac010_var_1382 = (1.0D0 - fs_var_1363) * fac10_var_1323(jl_var_1413, lay_var_1359)
        fac100_var_1380 = fs_var_1363 * fac00_var_1321(jl_var_1413, lay_var_1359)
        fac110_var_1383 = fs_var_1363 * fac10_var_1323(jl_var_1413, lay_var_1359)
        fac001_var_1385 = (1.0D0 - fs1_var_1366) * fac01_var_1322(jl_var_1413, lay_var_1359)
        fac011_var_1388 = (1.0D0 - fs1_var_1366) * fac11_var_1324(jl_var_1413, lay_var_1359)
        fac101_var_1386 = fs1_var_1366 * fac01_var_1322(jl_var_1413, lay_var_1359)
        fac111_var_1389 = fs1_var_1366 * fac11_var_1324(jl_var_1413, lay_var_1359)
        speccomb_mn2o_var_1347 = colh2o_var_1330(jl_var_1413, lay_var_1359) + refrat_m_b * colco2_var_1331(jl_var_1413, lay_var_1359)
        specparm_mn2o_var_1371 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_mn2o_var_1347, 0.999999D0)
        specmult_mn2o_var_1370 = 4.0D0 * specparm_mn2o_var_1371
        jmn2o_var_1361 = 1 + INT(specmult_mn2o_var_1370)
        fmn2o_var_1369 = ((specmult_mn2o_var_1370) - AINT((specmult_mn2o_var_1370)))
        chi_n2o_var_1378 = coln2o_var_1332(jl_var_1413, lay_var_1359) / coldry_var_1333(jl_var_1413, lay_var_1359)
        ratn2o_var_1377 = 1D+20 * chi_n2o_var_1378 / chi_mls(4, jp_var_1326(jl_var_1413, lay_var_1359) + 1)
        IF (ratn2o_var_1377 .GT. 1.5D0) THEN
          adjfac_var_1375 = 0.5D0 + (ratn2o_var_1377 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1376 = adjfac_var_1375 * chi_mls(4, jp_var_1326(jl_var_1413, lay_var_1359) + 1) * coldry_var_1333(jl_var_1413, lay_var_1359) * 1D-20
        ELSE
          adjcoln2o_var_1376 = coln2o_var_1332(jl_var_1413, lay_var_1359)
        END IF
        speccomb_planck_var_1348 = colh2o_var_1330(jl_var_1413, lay_var_1359) + refrat_planck_b_var_1350 * colco2_var_1331(jl_var_1413, lay_var_1359)
        specparm_planck_var_1374 = MIN(colh2o_var_1330(jl_var_1413, lay_var_1359) / speccomb_planck_var_1348, 0.999999D0)
        specmult_planck_var_1373 = 4.0D0 * specparm_planck_var_1374
        jpl_var_1362 = 1 + INT(specmult_planck_var_1373)
        fpl_var_1372 = ((specmult_planck_var_1373) - AINT((specmult_planck_var_1373)))
        ind0_var_1352 = ((jp_var_1326(jl_var_1413, lay_var_1359) - 13) * 5 + (jt_var_1327(jl_var_1413, lay_var_1359) - 1)) * nspb_var_217(3) + js_var_1358
        ind1_var_1353 = ((jp_var_1326(jl_var_1413, lay_var_1359) - 12) * 5 + (jt1_var_1328(jl_var_1413, lay_var_1359) - 1)) * nspb_var_217(3) + js1_var_1360
        indf_var_1355 = indfor_var_1341(jl_var_1413, lay_var_1359)
        indm_var_1356 = indminor_var_1344(jl_var_1413, lay_var_1359)
        DO ig_var_1357 = 1, 16
          taufor_var_1396 = forfac_var_1325(jl_var_1413, lay_var_1359) * (forref_var_166(indf_var_1355, ig_var_1357) + forfrac_var_1342(jl_var_1413, lay_var_1359) * (forref_var_166(indf_var_1355 + 1, ig_var_1357) - forref_var_166(indf_var_1355, ig_var_1357)))
          n2om1_var_1398 = kb_mn2o_var_162(jmn2o_var_1361, indm_var_1356, ig_var_1357) + fmn2o_var_1369 * (kb_mn2o_var_162(jmn2o_var_1361 + 1, indm_var_1356, ig_var_1357) - kb_mn2o_var_162(jmn2o_var_1361, indm_var_1356, ig_var_1357))
          n2om2_var_1399 = kb_mn2o_var_162(jmn2o_var_1361, indm_var_1356 + 1, ig_var_1357) + fmn2o_var_1369 * (kb_mn2o_var_162(jmn2o_var_1361 + 1, indm_var_1356 + 1, ig_var_1357) - kb_mn2o_var_162(jmn2o_var_1361, indm_var_1356 + 1, ig_var_1357))
          absn2o_var_1400 = n2om1_var_1398 + minorfrac_var_1343(jl_var_1413, lay_var_1359) * (n2om2_var_1399 - n2om1_var_1398)
          taug_var_1319(jl_var_1413, 22 + ig_var_1357, lay_var_1359) = speccomb_var_1345 * (fac000_var_1379 * absb_var_164(ind0_var_1352, ig_var_1357) + fac100_var_1380 * absb_var_164(ind0_var_1352 + 1, ig_var_1357) + fac010_var_1382 * absb_var_164(ind0_var_1352 + 5, ig_var_1357) + fac110_var_1383 * absb_var_164(ind0_var_1352 + 6, ig_var_1357)) + speccomb1_var_1346 * (fac001_var_1385 * absb_var_164(ind1_var_1353, ig_var_1357) + fac101_var_1386 * absb_var_164(ind1_var_1353 + 1, ig_var_1357) + fac011_var_1388 * absb_var_164(ind1_var_1353 + 5, ig_var_1357) + fac111_var_1389 * absb_var_164(ind1_var_1353 + 6, ig_var_1357)) + taufor_var_1396 + adjcoln2o_var_1376 * absn2o_var_1400
          fracs_var_1338(jl_var_1413, 22 + ig_var_1357, lay_var_1359) = fracrefb_var_160(ig_var_1357, jpl_var_1362) + fpl_var_1372 * (fracrefb_var_160(ig_var_1357, jpl_var_1362 + 1) - fracrefb_var_160(ig_var_1357, jpl_var_1362))
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol3
SUBROUTINE rrtm_taumol2(kidia_var_1414, kfdia_var_1415, klev_var_1416, taug_var_1417, pavel_var_1418, coldry_var_1419, p_tauaerl_var_1420, fac00_var_1421, fac01_var_1422, fac10_var_1423, fac11_var_1424, forfac_var_1426, forfrac_var_1425, indfor_var_1435, jp_var_1427, jt_var_1428, jt1_var_1429, colh2o_var_1430, laytrop_var_1431, selffac_var_1432, selffrac_var_1433, indself_var_1434, fracs_var_1436)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta2, ONLY: absa_var_155, absb_var_156, forref_var_158, fracrefa_var_153, fracrefb_var_154, selfref_var_157
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1414
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1415
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1416
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1417(kidia_var_1414 : kfdia_var_1415, 140, klev_var_1416)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1418(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1419(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1420(kidia_var_1414 : kfdia_var_1415, klev_var_1416, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1421(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1422(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1423(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1424(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1425(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1426(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1427(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1428(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1429(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1430(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1431(kidia_var_1414 : kfdia_var_1415)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1432(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1433(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1434(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1435(kidia_var_1414 : kfdia_var_1415, klev_var_1416)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1436(kidia_var_1414 : kfdia_var_1415, 140, klev_var_1416)
  INTEGER(KIND = 4) :: ind0_var_1437, ind1_var_1438, inds_var_1439, indf_var_1440
  INTEGER(KIND = 4) :: ig_var_1441, lay_var_1442
  REAL(KIND = 8) :: taufor_var_1443, tauself_var_1444, corradj_var_1445, pp_var_1446
  INTEGER(KIND = 4) :: laytrop_min_var_1447, laytrop_max_var_1448
  INTEGER(KIND = 4) :: ixc_var_1449(klev_var_1416), ixlow_var_1450(kfdia_var_1415, klev_var_1416), ixhigh_var_1451(kfdia_var_1415, klev_var_1416)
  INTEGER(KIND = 4) :: ich_var_1452, icl_var_1453, ixc0_var_1454, ixp_var_1455, jc_var_1456, jl_var_1457
  laytrop_min_var_1447 = MINVAL(laytrop_var_1431)
  laytrop_max_var_1448 = MAXVAL(laytrop_var_1431)
  ixlow_var_1450 = 0
  ixhigh_var_1451 = 0
  ixc_var_1449 = 0
  DO lay_var_1442 = laytrop_min_var_1447 + 1, laytrop_max_var_1448
    icl_var_1453 = 0
    ich_var_1452 = 0
    DO jc_var_1456 = kidia_var_1414, kfdia_var_1415
      IF (lay_var_1442 <= laytrop_var_1431(jc_var_1456)) THEN
        icl_var_1453 = icl_var_1453 + 1
        ixlow_var_1450(icl_var_1453, lay_var_1442) = jc_var_1456
      ELSE
        ich_var_1452 = ich_var_1452 + 1
        ixhigh_var_1451(ich_var_1452, lay_var_1442) = jc_var_1456
      END IF
    END DO
    ixc_var_1449(lay_var_1442) = icl_var_1453
  END DO
  DO lay_var_1442 = 1, laytrop_min_var_1447
    DO jl_var_1457 = kidia_var_1414, kfdia_var_1415
      ind0_var_1437 = ((jp_var_1427(jl_var_1457, lay_var_1442) - 1) * 5 + (jt_var_1428(jl_var_1457, lay_var_1442) - 1)) * nspa_var_216(2) + 1
      ind1_var_1438 = (jp_var_1427(jl_var_1457, lay_var_1442) * 5 + (jt1_var_1429(jl_var_1457, lay_var_1442) - 1)) * nspa_var_216(2) + 1
      inds_var_1439 = indself_var_1434(jl_var_1457, lay_var_1442)
      indf_var_1440 = indfor_var_1435(jl_var_1457, lay_var_1442)
      pp_var_1446 = pavel_var_1418(jl_var_1457, lay_var_1442)
      corradj_var_1445 = 1.0D0 - 0.05D0 * (pp_var_1446 - 100.0D0) / 900.0D0
      DO ig_var_1441 = 1, 12
        tauself_var_1444 = selffac_var_1432(jl_var_1457, lay_var_1442) * (selfref_var_157(inds_var_1439, ig_var_1441) + selffrac_var_1433(jl_var_1457, lay_var_1442) * (selfref_var_157(inds_var_1439 + 1, ig_var_1441) - selfref_var_157(inds_var_1439, ig_var_1441)))
        taufor_var_1443 = forfac_var_1426(jl_var_1457, lay_var_1442) * (forref_var_158(indf_var_1440, ig_var_1441) + forfrac_var_1425(jl_var_1457, lay_var_1442) * (forref_var_158(indf_var_1440 + 1, ig_var_1441) - forref_var_158(indf_var_1440, ig_var_1441)))
        taug_var_1417(jl_var_1457, 10 + ig_var_1441, lay_var_1442) = corradj_var_1445 * (colh2o_var_1430(jl_var_1457, lay_var_1442) * (fac00_var_1421(jl_var_1457, lay_var_1442) * absa_var_155(ind0_var_1437, ig_var_1441) + fac10_var_1423(jl_var_1457, lay_var_1442) * absa_var_155(ind0_var_1437 + 1, ig_var_1441) + fac01_var_1422(jl_var_1457, lay_var_1442) * absa_var_155(ind1_var_1438, ig_var_1441) + fac11_var_1424(jl_var_1457, lay_var_1442) * absa_var_155(ind1_var_1438 + 1, ig_var_1441)) + tauself_var_1444 + taufor_var_1443)
        fracs_var_1436(jl_var_1457, 10 + ig_var_1441, lay_var_1442) = fracrefa_var_153(ig_var_1441)
      END DO
    END DO
  END DO
  DO lay_var_1442 = laytrop_max_var_1448 + 1, klev_var_1416
    DO jl_var_1457 = kidia_var_1414, kfdia_var_1415
      ind0_var_1437 = ((jp_var_1427(jl_var_1457, lay_var_1442) - 13) * 5 + (jt_var_1428(jl_var_1457, lay_var_1442) - 1)) * nspb_var_217(2) + 1
      ind1_var_1438 = ((jp_var_1427(jl_var_1457, lay_var_1442) - 12) * 5 + (jt1_var_1429(jl_var_1457, lay_var_1442) - 1)) * nspb_var_217(2) + 1
      indf_var_1440 = indfor_var_1435(jl_var_1457, lay_var_1442)
      DO ig_var_1441 = 1, 12
        taufor_var_1443 = forfac_var_1426(jl_var_1457, lay_var_1442) * (forref_var_158(indf_var_1440, ig_var_1441) + forfrac_var_1425(jl_var_1457, lay_var_1442) * (forref_var_158(indf_var_1440 + 1, ig_var_1441) - forref_var_158(indf_var_1440, ig_var_1441)))
        taug_var_1417(jl_var_1457, 10 + ig_var_1441, lay_var_1442) = colh2o_var_1430(jl_var_1457, lay_var_1442) * (fac00_var_1421(jl_var_1457, lay_var_1442) * absb_var_156(ind0_var_1437, ig_var_1441) + fac10_var_1423(jl_var_1457, lay_var_1442) * absb_var_156(ind0_var_1437 + 1, ig_var_1441) + fac01_var_1422(jl_var_1457, lay_var_1442) * absb_var_156(ind1_var_1438, ig_var_1441) + fac11_var_1424(jl_var_1457, lay_var_1442) * absb_var_156(ind1_var_1438 + 1, ig_var_1441)) + taufor_var_1443
        fracs_var_1436(jl_var_1457, 10 + ig_var_1441, lay_var_1442) = fracrefb_var_154(ig_var_1441)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1448 /= laytrop_min_var_1447) THEN
    DO lay_var_1442 = laytrop_min_var_1447 + 1, laytrop_max_var_1448
      ixc0_var_1454 = ixc_var_1449(lay_var_1442)
      DO ixp_var_1455 = 1, ixc0_var_1454
        jl_var_1457 = ixlow_var_1450(ixp_var_1455, lay_var_1442)
        ind0_var_1437 = ((jp_var_1427(jl_var_1457, lay_var_1442) - 1) * 5 + (jt_var_1428(jl_var_1457, lay_var_1442) - 1)) * nspa_var_216(2) + 1
        ind1_var_1438 = (jp_var_1427(jl_var_1457, lay_var_1442) * 5 + (jt1_var_1429(jl_var_1457, lay_var_1442) - 1)) * nspa_var_216(2) + 1
        inds_var_1439 = indself_var_1434(jl_var_1457, lay_var_1442)
        indf_var_1440 = indfor_var_1435(jl_var_1457, lay_var_1442)
        pp_var_1446 = pavel_var_1418(jl_var_1457, lay_var_1442)
        corradj_var_1445 = 1.0D0 - 0.05D0 * (pp_var_1446 - 100.0D0) / 900.0D0
        DO ig_var_1441 = 1, 12
          tauself_var_1444 = selffac_var_1432(jl_var_1457, lay_var_1442) * (selfref_var_157(inds_var_1439, ig_var_1441) + selffrac_var_1433(jl_var_1457, lay_var_1442) * (selfref_var_157(inds_var_1439 + 1, ig_var_1441) - selfref_var_157(inds_var_1439, ig_var_1441)))
          taufor_var_1443 = forfac_var_1426(jl_var_1457, lay_var_1442) * (forref_var_158(indf_var_1440, ig_var_1441) + forfrac_var_1425(jl_var_1457, lay_var_1442) * (forref_var_158(indf_var_1440 + 1, ig_var_1441) - forref_var_158(indf_var_1440, ig_var_1441)))
          taug_var_1417(jl_var_1457, 10 + ig_var_1441, lay_var_1442) = corradj_var_1445 * (colh2o_var_1430(jl_var_1457, lay_var_1442) * (fac00_var_1421(jl_var_1457, lay_var_1442) * absa_var_155(ind0_var_1437, ig_var_1441) + fac10_var_1423(jl_var_1457, lay_var_1442) * absa_var_155(ind0_var_1437 + 1, ig_var_1441) + fac01_var_1422(jl_var_1457, lay_var_1442) * absa_var_155(ind1_var_1438, ig_var_1441) + fac11_var_1424(jl_var_1457, lay_var_1442) * absa_var_155(ind1_var_1438 + 1, ig_var_1441)) + tauself_var_1444 + taufor_var_1443)
          fracs_var_1436(jl_var_1457, 10 + ig_var_1441, lay_var_1442) = fracrefa_var_153(ig_var_1441)
        END DO
      END DO
      ixc0_var_1454 = kfdia_var_1415 - kidia_var_1414 + 1 - ixc0_var_1454
      DO ixp_var_1455 = 1, ixc0_var_1454
        jl_var_1457 = ixhigh_var_1451(ixp_var_1455, lay_var_1442)
        ind0_var_1437 = ((jp_var_1427(jl_var_1457, lay_var_1442) - 13) * 5 + (jt_var_1428(jl_var_1457, lay_var_1442) - 1)) * nspb_var_217(2) + 1
        ind1_var_1438 = ((jp_var_1427(jl_var_1457, lay_var_1442) - 12) * 5 + (jt1_var_1429(jl_var_1457, lay_var_1442) - 1)) * nspb_var_217(2) + 1
        indf_var_1440 = indfor_var_1435(jl_var_1457, lay_var_1442)
        DO ig_var_1441 = 1, 12
          taufor_var_1443 = forfac_var_1426(jl_var_1457, lay_var_1442) * (forref_var_158(indf_var_1440, ig_var_1441) + forfrac_var_1425(jl_var_1457, lay_var_1442) * (forref_var_158(indf_var_1440 + 1, ig_var_1441) - forref_var_158(indf_var_1440, ig_var_1441)))
          taug_var_1417(jl_var_1457, 10 + ig_var_1441, lay_var_1442) = colh2o_var_1430(jl_var_1457, lay_var_1442) * (fac00_var_1421(jl_var_1457, lay_var_1442) * absb_var_156(ind0_var_1437, ig_var_1441) + fac10_var_1423(jl_var_1457, lay_var_1442) * absb_var_156(ind0_var_1437 + 1, ig_var_1441) + fac01_var_1422(jl_var_1457, lay_var_1442) * absb_var_156(ind1_var_1438, ig_var_1441) + fac11_var_1424(jl_var_1457, lay_var_1442) * absb_var_156(ind1_var_1438 + 1, ig_var_1441)) + taufor_var_1443
          fracs_var_1436(jl_var_1457, 10 + ig_var_1441, lay_var_1442) = fracrefb_var_154(ig_var_1441)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol2
SUBROUTINE rrtm_taumol1(kidia_var_1458, kfdia_var_1459, klev_var_1460, taug_var_1462, pavel_var_1461, p_tauaerl_var_1463, fac00_var_1464, fac01_var_1465, fac10_var_1466, fac11_var_1467, forfac_var_1468, forfrac_var_1469, indfor_var_1480, jp_var_1470, jt_var_1471, jt1_var_1472, colh2o_var_1473, laytrop_var_1474, selffac_var_1475, selffrac_var_1476, indself_var_1478, fracs_var_1479, minorfrac_var_1477, indminor_var_1481, scaleminorn2, colbrd_var_1482)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta1, ONLY: absa_var_109, absb_var_110, forref_var_113, fracrefa_var_107, fracrefb_var_108, ka_mn2_var_111, kb_mn2, selfref_var_112
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1458
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1459
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1460
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1461(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1462(kidia_var_1458 : kfdia_var_1459, 140, klev_var_1460)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1463(kidia_var_1458 : kfdia_var_1459, klev_var_1460, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1464(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1465(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1466(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1467(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1468(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1469(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1470(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1471(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1472(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1473(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1474(kidia_var_1458 : kfdia_var_1459)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1475(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1476(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1477(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1478(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1479(kidia_var_1458 : kfdia_var_1459, 140, klev_var_1460)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1480(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1481(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(IN) :: scaleminorn2(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  REAL(KIND = 8), INTENT(IN) :: colbrd_var_1482(kidia_var_1458 : kfdia_var_1459, klev_var_1460)
  INTEGER(KIND = 4) :: ind0_var_1483, ind1_var_1484, inds_var_1485
  INTEGER(KIND = 4) :: indf_var_1486, indm_var_1487
  INTEGER(KIND = 4) :: ig_var_1488, lay_var_1489
  REAL(KIND = 8) :: taufor_var_1490, tauself_var_1491, corradj_var_1492, pp_var_1493, scalen2_var_1494, taun2_var_1495
  INTEGER(KIND = 4) :: laytrop_min_var_1496, laytrop_max_var_1497
  INTEGER(KIND = 4) :: ixc_var_1498(klev_var_1460), ixlow_var_1499(kfdia_var_1459, klev_var_1460), ixhigh_var_1500(kfdia_var_1459, klev_var_1460)
  INTEGER(KIND = 4) :: ich_var_1501, icl_var_1502, ixc0_var_1503, ixp_var_1504, jc_var_1505, jl_var_1506
  laytrop_min_var_1496 = MINVAL(laytrop_var_1474)
  laytrop_max_var_1497 = MAXVAL(laytrop_var_1474)
  ixlow_var_1499 = 0
  ixhigh_var_1500 = 0
  ixc_var_1498 = 0
  DO lay_var_1489 = laytrop_min_var_1496 + 1, laytrop_max_var_1497
    icl_var_1502 = 0
    ich_var_1501 = 0
    DO jc_var_1505 = kidia_var_1458, kfdia_var_1459
      IF (lay_var_1489 <= laytrop_var_1474(jc_var_1505)) THEN
        icl_var_1502 = icl_var_1502 + 1
        ixlow_var_1499(icl_var_1502, lay_var_1489) = jc_var_1505
      ELSE
        ich_var_1501 = ich_var_1501 + 1
        ixhigh_var_1500(ich_var_1501, lay_var_1489) = jc_var_1505
      END IF
    END DO
    ixc_var_1498(lay_var_1489) = icl_var_1502
  END DO
  DO lay_var_1489 = 1, laytrop_min_var_1496
    DO jl_var_1506 = kidia_var_1458, kfdia_var_1459
      ind0_var_1483 = ((jp_var_1470(jl_var_1506, lay_var_1489) - 1) * 5 + (jt_var_1471(jl_var_1506, lay_var_1489) - 1)) * nspa_var_216(1) + 1
      ind1_var_1484 = (jp_var_1470(jl_var_1506, lay_var_1489) * 5 + (jt1_var_1472(jl_var_1506, lay_var_1489) - 1)) * nspa_var_216(1) + 1
      inds_var_1485 = indself_var_1478(jl_var_1506, lay_var_1489)
      indf_var_1486 = indfor_var_1480(jl_var_1506, lay_var_1489)
      indm_var_1487 = indminor_var_1481(jl_var_1506, lay_var_1489)
      pp_var_1493 = pavel_var_1461(jl_var_1506, lay_var_1489)
      corradj_var_1492 = 1.0D0
      IF (pp_var_1493 .LT. 250.0D0) THEN
        corradj_var_1492 = 1.0D0 - 0.15D0 * (250.0D0 - pp_var_1493) / 154.4D0
      END IF
      scalen2_var_1494 = colbrd_var_1482(jl_var_1506, lay_var_1489) * scaleminorn2(jl_var_1506, lay_var_1489)
      DO ig_var_1488 = 1, 10
        tauself_var_1491 = selffac_var_1475(jl_var_1506, lay_var_1489) * (selfref_var_112(inds_var_1485, ig_var_1488) + selffrac_var_1476(jl_var_1506, lay_var_1489) * (selfref_var_112(inds_var_1485 + 1, ig_var_1488) - selfref_var_112(inds_var_1485, ig_var_1488)))
        taufor_var_1490 = forfac_var_1468(jl_var_1506, lay_var_1489) * (forref_var_113(indf_var_1486, ig_var_1488) + forfrac_var_1469(jl_var_1506, lay_var_1489) * (forref_var_113(indf_var_1486 + 1, ig_var_1488) - forref_var_113(indf_var_1486, ig_var_1488)))
        taun2_var_1495 = scalen2_var_1494 * (ka_mn2_var_111(indm_var_1487, ig_var_1488) + minorfrac_var_1477(jl_var_1506, lay_var_1489) * (ka_mn2_var_111(indm_var_1487 + 1, ig_var_1488) - ka_mn2_var_111(indm_var_1487, ig_var_1488)))
        taug_var_1462(jl_var_1506, ig_var_1488, lay_var_1489) = corradj_var_1492 * (colh2o_var_1473(jl_var_1506, lay_var_1489) * (fac00_var_1464(jl_var_1506, lay_var_1489) * absa_var_109(ind0_var_1483, ig_var_1488) + fac10_var_1466(jl_var_1506, lay_var_1489) * absa_var_109(ind0_var_1483 + 1, ig_var_1488) + fac01_var_1465(jl_var_1506, lay_var_1489) * absa_var_109(ind1_var_1484, ig_var_1488) + fac11_var_1467(jl_var_1506, lay_var_1489) * absa_var_109(ind1_var_1484 + 1, ig_var_1488)) + tauself_var_1491 + taufor_var_1490 + taun2_var_1495)
        fracs_var_1479(jl_var_1506, ig_var_1488, lay_var_1489) = fracrefa_var_107(ig_var_1488)
      END DO
    END DO
  END DO
  DO lay_var_1489 = laytrop_max_var_1497 + 1, klev_var_1460
    DO jl_var_1506 = kidia_var_1458, kfdia_var_1459
      ind0_var_1483 = ((jp_var_1470(jl_var_1506, lay_var_1489) - 13) * 5 + (jt_var_1471(jl_var_1506, lay_var_1489) - 1)) * nspb_var_217(1) + 1
      ind1_var_1484 = ((jp_var_1470(jl_var_1506, lay_var_1489) - 12) * 5 + (jt1_var_1472(jl_var_1506, lay_var_1489) - 1)) * nspb_var_217(1) + 1
      indf_var_1486 = indfor_var_1480(jl_var_1506, lay_var_1489)
      indm_var_1487 = indminor_var_1481(jl_var_1506, lay_var_1489)
      pp_var_1493 = pavel_var_1461(jl_var_1506, lay_var_1489)
      corradj_var_1492 = 1.0D0 - 0.15D0 * (pp_var_1493 / 95.6D0)
      scalen2_var_1494 = colbrd_var_1482(jl_var_1506, lay_var_1489) * scaleminorn2(jl_var_1506, lay_var_1489)
      DO ig_var_1488 = 1, 10
        taufor_var_1490 = forfac_var_1468(jl_var_1506, lay_var_1489) * (forref_var_113(indf_var_1486, ig_var_1488) + forfrac_var_1469(jl_var_1506, lay_var_1489) * (forref_var_113(indf_var_1486 + 1, ig_var_1488) - forref_var_113(indf_var_1486, ig_var_1488)))
        taun2_var_1495 = scalen2_var_1494 * (kb_mn2(indm_var_1487, ig_var_1488) + minorfrac_var_1477(jl_var_1506, lay_var_1489) * (kb_mn2(indm_var_1487 + 1, ig_var_1488) - kb_mn2(indm_var_1487, ig_var_1488)))
        taug_var_1462(jl_var_1506, ig_var_1488, lay_var_1489) = corradj_var_1492 * (colh2o_var_1473(jl_var_1506, lay_var_1489) * (fac00_var_1464(jl_var_1506, lay_var_1489) * absb_var_110(ind0_var_1483, ig_var_1488) + fac10_var_1466(jl_var_1506, lay_var_1489) * absb_var_110(ind0_var_1483 + 1, ig_var_1488) + fac01_var_1465(jl_var_1506, lay_var_1489) * absb_var_110(ind1_var_1484, ig_var_1488) + fac11_var_1467(jl_var_1506, lay_var_1489) * absb_var_110(ind1_var_1484 + 1, ig_var_1488)) + taufor_var_1490 + taun2_var_1495)
        fracs_var_1479(jl_var_1506, ig_var_1488, lay_var_1489) = fracrefb_var_108(ig_var_1488)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1497 /= laytrop_min_var_1496) THEN
    DO lay_var_1489 = laytrop_min_var_1496 + 1, laytrop_max_var_1497
      ixc0_var_1503 = ixc_var_1498(lay_var_1489)
      DO ixp_var_1504 = 1, ixc0_var_1503
        jl_var_1506 = ixlow_var_1499(ixp_var_1504, lay_var_1489)
        ind0_var_1483 = ((jp_var_1470(jl_var_1506, lay_var_1489) - 1) * 5 + (jt_var_1471(jl_var_1506, lay_var_1489) - 1)) * nspa_var_216(1) + 1
        ind1_var_1484 = (jp_var_1470(jl_var_1506, lay_var_1489) * 5 + (jt1_var_1472(jl_var_1506, lay_var_1489) - 1)) * nspa_var_216(1) + 1
        inds_var_1485 = indself_var_1478(jl_var_1506, lay_var_1489)
        indf_var_1486 = indfor_var_1480(jl_var_1506, lay_var_1489)
        indm_var_1487 = indminor_var_1481(jl_var_1506, lay_var_1489)
        pp_var_1493 = pavel_var_1461(jl_var_1506, lay_var_1489)
        corradj_var_1492 = 1.0D0
        IF (pp_var_1493 .LT. 250.0D0) THEN
          corradj_var_1492 = 1.0D0 - 0.15D0 * (250.0D0 - pp_var_1493) / 154.4D0
        END IF
        scalen2_var_1494 = colbrd_var_1482(jl_var_1506, lay_var_1489) * scaleminorn2(jl_var_1506, lay_var_1489)
        DO ig_var_1488 = 1, 10
          tauself_var_1491 = selffac_var_1475(jl_var_1506, lay_var_1489) * (selfref_var_112(inds_var_1485, ig_var_1488) + selffrac_var_1476(jl_var_1506, lay_var_1489) * (selfref_var_112(inds_var_1485 + 1, ig_var_1488) - selfref_var_112(inds_var_1485, ig_var_1488)))
          taufor_var_1490 = forfac_var_1468(jl_var_1506, lay_var_1489) * (forref_var_113(indf_var_1486, ig_var_1488) + forfrac_var_1469(jl_var_1506, lay_var_1489) * (forref_var_113(indf_var_1486 + 1, ig_var_1488) - forref_var_113(indf_var_1486, ig_var_1488)))
          taun2_var_1495 = scalen2_var_1494 * (ka_mn2_var_111(indm_var_1487, ig_var_1488) + minorfrac_var_1477(jl_var_1506, lay_var_1489) * (ka_mn2_var_111(indm_var_1487 + 1, ig_var_1488) - ka_mn2_var_111(indm_var_1487, ig_var_1488)))
          taug_var_1462(jl_var_1506, ig_var_1488, lay_var_1489) = corradj_var_1492 * (colh2o_var_1473(jl_var_1506, lay_var_1489) * (fac00_var_1464(jl_var_1506, lay_var_1489) * absa_var_109(ind0_var_1483, ig_var_1488) + fac10_var_1466(jl_var_1506, lay_var_1489) * absa_var_109(ind0_var_1483 + 1, ig_var_1488) + fac01_var_1465(jl_var_1506, lay_var_1489) * absa_var_109(ind1_var_1484, ig_var_1488) + fac11_var_1467(jl_var_1506, lay_var_1489) * absa_var_109(ind1_var_1484 + 1, ig_var_1488)) + tauself_var_1491 + taufor_var_1490 + taun2_var_1495)
          fracs_var_1479(jl_var_1506, ig_var_1488, lay_var_1489) = fracrefa_var_107(ig_var_1488)
        END DO
      END DO
      ixc0_var_1503 = kfdia_var_1459 - kidia_var_1458 + 1 - ixc0_var_1503
      DO ixp_var_1504 = 1, ixc0_var_1503
        jl_var_1506 = ixhigh_var_1500(ixp_var_1504, lay_var_1489)
        ind0_var_1483 = ((jp_var_1470(jl_var_1506, lay_var_1489) - 13) * 5 + (jt_var_1471(jl_var_1506, lay_var_1489) - 1)) * nspb_var_217(1) + 1
        ind1_var_1484 = ((jp_var_1470(jl_var_1506, lay_var_1489) - 12) * 5 + (jt1_var_1472(jl_var_1506, lay_var_1489) - 1)) * nspb_var_217(1) + 1
        indf_var_1486 = indfor_var_1480(jl_var_1506, lay_var_1489)
        indm_var_1487 = indminor_var_1481(jl_var_1506, lay_var_1489)
        pp_var_1493 = pavel_var_1461(jl_var_1506, lay_var_1489)
        corradj_var_1492 = 1.0D0 - 0.15D0 * (pp_var_1493 / 95.6D0)
        scalen2_var_1494 = colbrd_var_1482(jl_var_1506, lay_var_1489) * scaleminorn2(jl_var_1506, lay_var_1489)
        DO ig_var_1488 = 1, 10
          taufor_var_1490 = forfac_var_1468(jl_var_1506, lay_var_1489) * (forref_var_113(indf_var_1486, ig_var_1488) + forfrac_var_1469(jl_var_1506, lay_var_1489) * (forref_var_113(indf_var_1486 + 1, ig_var_1488) - forref_var_113(indf_var_1486, ig_var_1488)))
          taun2_var_1495 = scalen2_var_1494 * (kb_mn2(indm_var_1487, ig_var_1488) + minorfrac_var_1477(jl_var_1506, lay_var_1489) * (kb_mn2(indm_var_1487 + 1, ig_var_1488) - kb_mn2(indm_var_1487, ig_var_1488)))
          taug_var_1462(jl_var_1506, ig_var_1488, lay_var_1489) = corradj_var_1492 * (colh2o_var_1473(jl_var_1506, lay_var_1489) * (fac00_var_1464(jl_var_1506, lay_var_1489) * absb_var_110(ind0_var_1483, ig_var_1488) + fac10_var_1466(jl_var_1506, lay_var_1489) * absb_var_110(ind0_var_1483 + 1, ig_var_1488) + fac01_var_1465(jl_var_1506, lay_var_1489) * absb_var_110(ind1_var_1484, ig_var_1488) + fac11_var_1467(jl_var_1506, lay_var_1489) * absb_var_110(ind1_var_1484 + 1, ig_var_1488)) + taufor_var_1490 + taun2_var_1495)
          fracs_var_1479(jl_var_1506, ig_var_1488, lay_var_1489) = fracrefb_var_108(ig_var_1488)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol1
SUBROUTINE rrtm_taumol14(kidia_var_1507, kfdia_var_1508, klev_var_1509, taug_var_1510, p_tauaerl_var_1511, fac00_var_1512, fac01_var_1513, fac10_var_1514, fac11_var_1515, forfac_var_1526, forfrac_var_1527, indfor_var_1525, jp_var_1516, jt_var_1517, jt1_var_1518, colco2_var_1519, laytrop_var_1520, selffac_var_1521, selffrac_var_1522, indself_var_1523, fracs_var_1524)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta14, ONLY: absa_var_138, absb_var_139, forref_var_141, fracrefa_var_136, fracrefb_var_137, selfref_var_140
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1507
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1508
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1509
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1510(kidia_var_1507 : kfdia_var_1508, 140, klev_var_1509)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1511(kidia_var_1507 : kfdia_var_1508, klev_var_1509, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1512(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1513(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1514(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1515(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1516(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1517(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1518(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1519(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1520(kidia_var_1507 : kfdia_var_1508)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1521(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1522(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1523(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1524(kidia_var_1507 : kfdia_var_1508, 140, klev_var_1509)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1525(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1526(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1527(kidia_var_1507 : kfdia_var_1508, klev_var_1509)
  INTEGER(KIND = 4) :: ig_var_1528, ind0_var_1529, ind1_var_1530, inds_var_1531, indf_var_1532, lay_var_1533
  REAL(KIND = 8) :: taufor_var_1534, tauself_var_1535
  INTEGER(KIND = 4) :: laytrop_min_var_1536, laytrop_max_var_1537
  INTEGER(KIND = 4) :: ixc_var_1538(klev_var_1509), ixlow_var_1539(kfdia_var_1508, klev_var_1509), ixhigh_var_1540(kfdia_var_1508, klev_var_1509)
  INTEGER(KIND = 4) :: ich_var_1541, icl_var_1542, ixc0_var_1543, ixp_var_1544, jc_var_1545, jl_var_1546
  laytrop_min_var_1536 = MINVAL(laytrop_var_1520)
  laytrop_max_var_1537 = MAXVAL(laytrop_var_1520)
  ixlow_var_1539 = 0
  ixhigh_var_1540 = 0
  ixc_var_1538 = 0
  DO lay_var_1533 = laytrop_min_var_1536 + 1, laytrop_max_var_1537
    icl_var_1542 = 0
    ich_var_1541 = 0
    DO jc_var_1545 = kidia_var_1507, kfdia_var_1508
      IF (lay_var_1533 <= laytrop_var_1520(jc_var_1545)) THEN
        icl_var_1542 = icl_var_1542 + 1
        ixlow_var_1539(icl_var_1542, lay_var_1533) = jc_var_1545
      ELSE
        ich_var_1541 = ich_var_1541 + 1
        ixhigh_var_1540(ich_var_1541, lay_var_1533) = jc_var_1545
      END IF
    END DO
    ixc_var_1538(lay_var_1533) = icl_var_1542
  END DO
  DO lay_var_1533 = 1, laytrop_min_var_1536
    DO jl_var_1546 = kidia_var_1507, kfdia_var_1508
      ind0_var_1529 = ((jp_var_1516(jl_var_1546, lay_var_1533) - 1) * 5 + (jt_var_1517(jl_var_1546, lay_var_1533) - 1)) * nspa_var_216(14) + 1
      ind1_var_1530 = (jp_var_1516(jl_var_1546, lay_var_1533) * 5 + (jt1_var_1518(jl_var_1546, lay_var_1533) - 1)) * nspa_var_216(14) + 1
      inds_var_1531 = indself_var_1523(jl_var_1546, lay_var_1533)
      indf_var_1532 = indfor_var_1525(jl_var_1546, lay_var_1533)
      DO ig_var_1528 = 1, 2
        tauself_var_1535 = selffac_var_1521(jl_var_1546, lay_var_1533) * (selfref_var_140(inds_var_1531, ig_var_1528) + selffrac_var_1522(jl_var_1546, lay_var_1533) * (selfref_var_140(inds_var_1531 + 1, ig_var_1528) - selfref_var_140(inds_var_1531, ig_var_1528)))
        taufor_var_1534 = forfac_var_1526(jl_var_1546, lay_var_1533) * (forref_var_141(indf_var_1532, ig_var_1528) + forfrac_var_1527(jl_var_1546, lay_var_1533) * (forref_var_141(indf_var_1532 + 1, ig_var_1528) - forref_var_141(indf_var_1532, ig_var_1528)))
        taug_var_1510(jl_var_1546, 134 + ig_var_1528, lay_var_1533) = colco2_var_1519(jl_var_1546, lay_var_1533) * (fac00_var_1512(jl_var_1546, lay_var_1533) * absa_var_138(ind0_var_1529, ig_var_1528) + fac10_var_1514(jl_var_1546, lay_var_1533) * absa_var_138(ind0_var_1529 + 1, ig_var_1528) + fac01_var_1513(jl_var_1546, lay_var_1533) * absa_var_138(ind1_var_1530, ig_var_1528) + fac11_var_1515(jl_var_1546, lay_var_1533) * absa_var_138(ind1_var_1530 + 1, ig_var_1528)) + tauself_var_1535 + taufor_var_1534
        fracs_var_1524(jl_var_1546, 134 + ig_var_1528, lay_var_1533) = fracrefa_var_136(ig_var_1528)
      END DO
    END DO
  END DO
  DO lay_var_1533 = laytrop_max_var_1537 + 1, klev_var_1509
    DO jl_var_1546 = kidia_var_1507, kfdia_var_1508
      ind0_var_1529 = ((jp_var_1516(jl_var_1546, lay_var_1533) - 13) * 5 + (jt_var_1517(jl_var_1546, lay_var_1533) - 1)) * nspb_var_217(14) + 1
      ind1_var_1530 = ((jp_var_1516(jl_var_1546, lay_var_1533) - 12) * 5 + (jt1_var_1518(jl_var_1546, lay_var_1533) - 1)) * nspb_var_217(14) + 1
      DO ig_var_1528 = 1, 2
        taug_var_1510(jl_var_1546, 134 + ig_var_1528, lay_var_1533) = colco2_var_1519(jl_var_1546, lay_var_1533) * (fac00_var_1512(jl_var_1546, lay_var_1533) * absb_var_139(ind0_var_1529, ig_var_1528) + fac10_var_1514(jl_var_1546, lay_var_1533) * absb_var_139(ind0_var_1529 + 1, ig_var_1528) + fac01_var_1513(jl_var_1546, lay_var_1533) * absb_var_139(ind1_var_1530, ig_var_1528) + fac11_var_1515(jl_var_1546, lay_var_1533) * absb_var_139(ind1_var_1530 + 1, ig_var_1528))
        fracs_var_1524(jl_var_1546, 134 + ig_var_1528, lay_var_1533) = fracrefb_var_137(ig_var_1528)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1537 /= laytrop_min_var_1536) THEN
    DO lay_var_1533 = laytrop_min_var_1536 + 1, laytrop_max_var_1537
      ixc0_var_1543 = ixc_var_1538(lay_var_1533)
      DO ixp_var_1544 = 1, ixc0_var_1543
        jl_var_1546 = ixlow_var_1539(ixp_var_1544, lay_var_1533)
        ind0_var_1529 = ((jp_var_1516(jl_var_1546, lay_var_1533) - 1) * 5 + (jt_var_1517(jl_var_1546, lay_var_1533) - 1)) * nspa_var_216(14) + 1
        ind1_var_1530 = (jp_var_1516(jl_var_1546, lay_var_1533) * 5 + (jt1_var_1518(jl_var_1546, lay_var_1533) - 1)) * nspa_var_216(14) + 1
        inds_var_1531 = indself_var_1523(jl_var_1546, lay_var_1533)
        indf_var_1532 = indfor_var_1525(jl_var_1546, lay_var_1533)
        DO ig_var_1528 = 1, 2
          tauself_var_1535 = selffac_var_1521(jl_var_1546, lay_var_1533) * (selfref_var_140(inds_var_1531, ig_var_1528) + selffrac_var_1522(jl_var_1546, lay_var_1533) * (selfref_var_140(inds_var_1531 + 1, ig_var_1528) - selfref_var_140(inds_var_1531, ig_var_1528)))
          taufor_var_1534 = forfac_var_1526(jl_var_1546, lay_var_1533) * (forref_var_141(indf_var_1532, ig_var_1528) + forfrac_var_1527(jl_var_1546, lay_var_1533) * (forref_var_141(indf_var_1532 + 1, ig_var_1528) - forref_var_141(indf_var_1532, ig_var_1528)))
          taug_var_1510(jl_var_1546, 134 + ig_var_1528, lay_var_1533) = colco2_var_1519(jl_var_1546, lay_var_1533) * (fac00_var_1512(jl_var_1546, lay_var_1533) * absa_var_138(ind0_var_1529, ig_var_1528) + fac10_var_1514(jl_var_1546, lay_var_1533) * absa_var_138(ind0_var_1529 + 1, ig_var_1528) + fac01_var_1513(jl_var_1546, lay_var_1533) * absa_var_138(ind1_var_1530, ig_var_1528) + fac11_var_1515(jl_var_1546, lay_var_1533) * absa_var_138(ind1_var_1530 + 1, ig_var_1528)) + tauself_var_1535 + taufor_var_1534
          fracs_var_1524(jl_var_1546, 134 + ig_var_1528, lay_var_1533) = fracrefa_var_136(ig_var_1528)
        END DO
      END DO
      ixc0_var_1543 = kfdia_var_1508 - kidia_var_1507 + 1 - ixc0_var_1543
      DO ixp_var_1544 = 1, ixc0_var_1543
        jl_var_1546 = ixhigh_var_1540(ixp_var_1544, lay_var_1533)
        ind0_var_1529 = ((jp_var_1516(jl_var_1546, lay_var_1533) - 13) * 5 + (jt_var_1517(jl_var_1546, lay_var_1533) - 1)) * nspb_var_217(14) + 1
        ind1_var_1530 = ((jp_var_1516(jl_var_1546, lay_var_1533) - 12) * 5 + (jt1_var_1518(jl_var_1546, lay_var_1533) - 1)) * nspb_var_217(14) + 1
        DO ig_var_1528 = 1, 2
          taug_var_1510(jl_var_1546, 134 + ig_var_1528, lay_var_1533) = colco2_var_1519(jl_var_1546, lay_var_1533) * (fac00_var_1512(jl_var_1546, lay_var_1533) * absb_var_139(ind0_var_1529, ig_var_1528) + fac10_var_1514(jl_var_1546, lay_var_1533) * absb_var_139(ind0_var_1529 + 1, ig_var_1528) + fac01_var_1513(jl_var_1546, lay_var_1533) * absb_var_139(ind1_var_1530, ig_var_1528) + fac11_var_1515(jl_var_1546, lay_var_1533) * absb_var_139(ind1_var_1530 + 1, ig_var_1528))
          fracs_var_1524(jl_var_1546, 134 + ig_var_1528, lay_var_1533) = fracrefb_var_137(ig_var_1528)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol14
SUBROUTINE rrtm_taumol15(kidia_var_1547, kfdia_var_1548, klev_var_1549, taug_var_1550, p_tauaerl_var_1551, fac00_var_1552, fac01_var_1553, fac10_var_1554, fac11_var_1555, forfac_var_1569, forfrac_var_1570, indfor_var_1568, jp_var_1556, jt_var_1557, jt1_var_1558, oneminus_var_1559, colh2o_var_1560, colco2_var_1561, coln2o_var_1562, laytrop_var_1563, selffac_var_1564, selffrac_var_1565, indself_var_1566, fracs_var_1567, rat_n2oco2, rat_n2oco2_1, minorfrac_var_1571, indminor_var_1572, scaleminor_var_1573, colbrd_var_1574)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrtm, ONLY: ng15
  USE yoerrta15, ONLY: absa_var_143, forref_var_146, fracrefa_var_142, ka_mn2_var_144, selfref_var_145
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1547
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1548
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1549
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1550(kidia_var_1547 : kfdia_var_1548, 140, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1551(kidia_var_1547 : kfdia_var_1548, klev_var_1549, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1552(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1553(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1554(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1555(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1556(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1557(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1558(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1559
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1560(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1561(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_1562(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1563(kidia_var_1547 : kfdia_var_1548)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1564(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1565(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1566(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1567(kidia_var_1547 : kfdia_var_1548, 140, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: rat_n2oco2(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: rat_n2oco2_1(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1568(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1569(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1570(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1571(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1572(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: scaleminor_var_1573(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  REAL(KIND = 8), INTENT(IN) :: colbrd_var_1574(kidia_var_1547 : kfdia_var_1548, klev_var_1549)
  INTEGER(KIND = 4) :: ig_var_1575, ind0_var_1576, ind1_var_1577, inds_var_1578, indf_var_1579, indm_var_1580, js_var_1581, js1_var_1582, jpl_var_1583, jmn2, lay_var_1584
  REAL(KIND = 8) :: refrat_planck_a_var_1585, refrat_m_a_var_1586
  REAL(KIND = 8) :: taufor_var_1587, tauself_var_1588, tau_major_var_1589(2), tau_major1_var_1590(2), n2m1, n2m2, taun2_var_1591, scalen2_var_1592
  REAL(KIND = 8) :: fac000_var_1593, fac100_var_1594, fac200_var_1595, fac010_var_1596, fac110_var_1597, fac210_var_1598, fac001_var_1599, fac101_var_1600, fac201_var_1601, fac011_var_1602, fac111_var_1603, fac211_var_1604
  REAL(KIND = 8) :: p_var_1605, p4_var_1606, fk0_var_1607, fk1_var_1608, fk2_var_1609
  REAL(KIND = 8) :: fs_var_1610, specmult_var_1611, specparm_var_1612, speccomb_var_1613, fs1_var_1614, specmult1_var_1615, specparm1_var_1616, speccomb1_var_1617, fmn2, specmult_mn2, specparm_mn2, speccomb_mn2, fpl_var_1618, specmult_planck_var_1619, specparm_planck_var_1620, speccomb_planck_var_1621
  INTEGER(KIND = 4) :: laytrop_min_var_1622, laytrop_max_var_1623
  INTEGER(KIND = 4) :: ixc_var_1624(klev_var_1549), ixlow_var_1625(kfdia_var_1548, klev_var_1549), ixhigh_var_1626(kfdia_var_1548, klev_var_1549)
  INTEGER(KIND = 4) :: ich_var_1627, icl_var_1628, ixc0_var_1629, ixp_var_1630, jc_var_1631, jl_var_1632
  laytrop_min_var_1622 = MINVAL(laytrop_var_1563)
  laytrop_max_var_1623 = MAXVAL(laytrop_var_1563)
  ixlow_var_1625 = 0
  ixhigh_var_1626 = 0
  ixc_var_1624 = 0
  DO lay_var_1584 = laytrop_min_var_1622 + 1, laytrop_max_var_1623
    icl_var_1628 = 0
    ich_var_1627 = 0
    DO jc_var_1631 = kidia_var_1547, kfdia_var_1548
      IF (lay_var_1584 <= laytrop_var_1563(jc_var_1631)) THEN
        icl_var_1628 = icl_var_1628 + 1
        ixlow_var_1625(icl_var_1628, lay_var_1584) = jc_var_1631
      ELSE
        ich_var_1627 = ich_var_1627 + 1
        ixhigh_var_1626(ich_var_1627, lay_var_1584) = jc_var_1631
      END IF
    END DO
    ixc_var_1624(lay_var_1584) = icl_var_1628
  END DO
  refrat_planck_a_var_1585 = chi_mls(4, 1) / chi_mls(2, 1)
  refrat_m_a_var_1586 = chi_mls(4, 1) / chi_mls(2, 1)
  DO lay_var_1584 = 1, laytrop_min_var_1622
    DO jl_var_1632 = kidia_var_1547, kfdia_var_1548
      speccomb_var_1613 = coln2o_var_1562(jl_var_1632, lay_var_1584) + rat_n2oco2(jl_var_1632, lay_var_1584) * colco2_var_1561(jl_var_1632, lay_var_1584)
      specparm_var_1612 = MIN(coln2o_var_1562(jl_var_1632, lay_var_1584) / speccomb_var_1613, 0.999999D0)
      specmult_var_1611 = 8.0D0 * (specparm_var_1612)
      js_var_1581 = 1 + INT(specmult_var_1611)
      fs_var_1610 = ((specmult_var_1611) - AINT((specmult_var_1611)))
      speccomb1_var_1617 = coln2o_var_1562(jl_var_1632, lay_var_1584) + rat_n2oco2_1(jl_var_1632, lay_var_1584) * colco2_var_1561(jl_var_1632, lay_var_1584)
      specparm1_var_1616 = MIN(coln2o_var_1562(jl_var_1632, lay_var_1584) / speccomb1_var_1617, 0.999999D0)
      specmult1_var_1615 = 8.0D0 * (specparm1_var_1616)
      js1_var_1582 = 1 + INT(specmult1_var_1615)
      fs1_var_1614 = ((specmult1_var_1615) - AINT((specmult1_var_1615)))
      speccomb_mn2 = coln2o_var_1562(jl_var_1632, lay_var_1584) + refrat_m_a_var_1586 * colco2_var_1561(jl_var_1632, lay_var_1584)
      specparm_mn2 = MIN(coln2o_var_1562(jl_var_1632, lay_var_1584) / speccomb_mn2, 0.999999D0)
      specmult_mn2 = 8.0D0 * specparm_mn2
      jmn2 = 1 + INT(specmult_mn2)
      fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))
      speccomb_planck_var_1621 = coln2o_var_1562(jl_var_1632, lay_var_1584) + refrat_planck_a_var_1585 * colco2_var_1561(jl_var_1632, lay_var_1584)
      specparm_planck_var_1620 = MIN(coln2o_var_1562(jl_var_1632, lay_var_1584) / speccomb_planck_var_1621, 0.999999D0)
      specmult_planck_var_1619 = 8.0D0 * specparm_planck_var_1620
      jpl_var_1583 = 1 + INT(specmult_planck_var_1619)
      fpl_var_1618 = ((specmult_planck_var_1619) - AINT((specmult_planck_var_1619)))
      ind0_var_1576 = ((jp_var_1556(jl_var_1632, lay_var_1584) - 1) * 5 + (jt_var_1557(jl_var_1632, lay_var_1584) - 1)) * nspa_var_216(15) + js_var_1581
      ind1_var_1577 = (jp_var_1556(jl_var_1632, lay_var_1584) * 5 + (jt1_var_1558(jl_var_1632, lay_var_1584) - 1)) * nspa_var_216(15) + js1_var_1582
      inds_var_1578 = indself_var_1566(jl_var_1632, lay_var_1584)
      indf_var_1579 = indfor_var_1568(jl_var_1632, lay_var_1584)
      indm_var_1580 = indminor_var_1572(jl_var_1632, lay_var_1584)
      scalen2_var_1592 = colbrd_var_1574(jl_var_1632, lay_var_1584) * scaleminor_var_1573(jl_var_1632, lay_var_1584)
      IF (specparm_var_1612 .LT. 0.125D0) THEN
        p_var_1605 = fs_var_1610 - 1.0D0
        p4_var_1606 = p_var_1605 ** 4
        fk0_var_1607 = p4_var_1606
        fk1_var_1608 = 1.0D0 - p_var_1605 - 2.0D0 * p4_var_1606
        fk2_var_1609 = p_var_1605 + p4_var_1606
        fac000_var_1593 = fk0_var_1607 * fac00_var_1552(jl_var_1632, lay_var_1584)
        fac100_var_1594 = fk1_var_1608 * fac00_var_1552(jl_var_1632, lay_var_1584)
        fac200_var_1595 = fk2_var_1609 * fac00_var_1552(jl_var_1632, lay_var_1584)
        fac010_var_1596 = fk0_var_1607 * fac10_var_1554(jl_var_1632, lay_var_1584)
        fac110_var_1597 = fk1_var_1608 * fac10_var_1554(jl_var_1632, lay_var_1584)
        fac210_var_1598 = fk2_var_1609 * fac10_var_1554(jl_var_1632, lay_var_1584)
      ELSE IF (specparm_var_1612 .GT. 0.875D0) THEN
        p_var_1605 = - fs_var_1610
        p4_var_1606 = p_var_1605 ** 4
        fk0_var_1607 = p4_var_1606
        fk1_var_1608 = 1.0D0 - p_var_1605 - 2.0D0 * p4_var_1606
        fk2_var_1609 = p_var_1605 + p4_var_1606
        fac000_var_1593 = fk0_var_1607 * fac00_var_1552(jl_var_1632, lay_var_1584)
        fac100_var_1594 = fk1_var_1608 * fac00_var_1552(jl_var_1632, lay_var_1584)
        fac200_var_1595 = fk2_var_1609 * fac00_var_1552(jl_var_1632, lay_var_1584)
        fac010_var_1596 = fk0_var_1607 * fac10_var_1554(jl_var_1632, lay_var_1584)
        fac110_var_1597 = fk1_var_1608 * fac10_var_1554(jl_var_1632, lay_var_1584)
        fac210_var_1598 = fk2_var_1609 * fac10_var_1554(jl_var_1632, lay_var_1584)
      ELSE
        fac000_var_1593 = (1.0D0 - fs_var_1610) * fac00_var_1552(jl_var_1632, lay_var_1584)
        fac010_var_1596 = (1.0D0 - fs_var_1610) * fac10_var_1554(jl_var_1632, lay_var_1584)
        fac100_var_1594 = fs_var_1610 * fac00_var_1552(jl_var_1632, lay_var_1584)
        fac110_var_1597 = fs_var_1610 * fac10_var_1554(jl_var_1632, lay_var_1584)
        fac200_var_1595 = 0.0D0
        fac210_var_1598 = 0.0D0
      END IF
      IF (specparm1_var_1616 .LT. 0.125D0) THEN
        p_var_1605 = fs1_var_1614 - 1.0D0
        p4_var_1606 = p_var_1605 ** 4
        fk0_var_1607 = p4_var_1606
        fk1_var_1608 = 1.0D0 - p_var_1605 - 2.0D0 * p4_var_1606
        fk2_var_1609 = p_var_1605 + p4_var_1606
        fac001_var_1599 = fk0_var_1607 * fac01_var_1553(jl_var_1632, lay_var_1584)
        fac101_var_1600 = fk1_var_1608 * fac01_var_1553(jl_var_1632, lay_var_1584)
        fac201_var_1601 = fk2_var_1609 * fac01_var_1553(jl_var_1632, lay_var_1584)
        fac011_var_1602 = fk0_var_1607 * fac11_var_1555(jl_var_1632, lay_var_1584)
        fac111_var_1603 = fk1_var_1608 * fac11_var_1555(jl_var_1632, lay_var_1584)
        fac211_var_1604 = fk2_var_1609 * fac11_var_1555(jl_var_1632, lay_var_1584)
      ELSE IF (specparm1_var_1616 .GT. 0.875D0) THEN
        p_var_1605 = - fs1_var_1614
        p4_var_1606 = p_var_1605 ** 4
        fk0_var_1607 = p4_var_1606
        fk1_var_1608 = 1.0D0 - p_var_1605 - 2.0D0 * p4_var_1606
        fk2_var_1609 = p_var_1605 + p4_var_1606
        fac001_var_1599 = fk0_var_1607 * fac01_var_1553(jl_var_1632, lay_var_1584)
        fac101_var_1600 = fk1_var_1608 * fac01_var_1553(jl_var_1632, lay_var_1584)
        fac201_var_1601 = fk2_var_1609 * fac01_var_1553(jl_var_1632, lay_var_1584)
        fac011_var_1602 = fk0_var_1607 * fac11_var_1555(jl_var_1632, lay_var_1584)
        fac111_var_1603 = fk1_var_1608 * fac11_var_1555(jl_var_1632, lay_var_1584)
        fac211_var_1604 = fk2_var_1609 * fac11_var_1555(jl_var_1632, lay_var_1584)
      ELSE
        fac001_var_1599 = (1.0D0 - fs1_var_1614) * fac01_var_1553(jl_var_1632, lay_var_1584)
        fac011_var_1602 = (1.0D0 - fs1_var_1614) * fac11_var_1555(jl_var_1632, lay_var_1584)
        fac101_var_1600 = fs1_var_1614 * fac01_var_1553(jl_var_1632, lay_var_1584)
        fac111_var_1603 = fs1_var_1614 * fac11_var_1555(jl_var_1632, lay_var_1584)
        fac201_var_1601 = 0.0D0
        fac211_var_1604 = 0.0D0
      END IF
      IF (specparm_var_1612 .LT. 0.125D0) THEN
        tau_major_var_1589(1 : ng15) = speccomb_var_1613 * (fac000_var_1593 * absa_var_143(ind0_var_1576, 1 : 2) + fac100_var_1594 * absa_var_143(ind0_var_1576 + 1, 1 : 2) + fac200_var_1595 * absa_var_143(ind0_var_1576 + 2, 1 : 2) + fac010_var_1596 * absa_var_143(ind0_var_1576 + 9, 1 : 2) + fac110_var_1597 * absa_var_143(ind0_var_1576 + 10, 1 : 2) + fac210_var_1598 * absa_var_143(ind0_var_1576 + 11, 1 : 2))
      ELSE IF (specparm_var_1612 .GT. 0.875D0) THEN
        tau_major_var_1589(1 : ng15) = speccomb_var_1613 * (fac200_var_1595 * absa_var_143(ind0_var_1576 - 1, 1 : 2) + fac100_var_1594 * absa_var_143(ind0_var_1576, 1 : 2) + fac000_var_1593 * absa_var_143(ind0_var_1576 + 1, 1 : 2) + fac210_var_1598 * absa_var_143(ind0_var_1576 + 8, 1 : 2) + fac110_var_1597 * absa_var_143(ind0_var_1576 + 9, 1 : 2) + fac010_var_1596 * absa_var_143(ind0_var_1576 + 10, 1 : 2))
      ELSE
        tau_major_var_1589(1 : ng15) = speccomb_var_1613 * (fac000_var_1593 * absa_var_143(ind0_var_1576, 1 : 2) + fac100_var_1594 * absa_var_143(ind0_var_1576 + 1, 1 : 2) + fac010_var_1596 * absa_var_143(ind0_var_1576 + 9, 1 : 2) + fac110_var_1597 * absa_var_143(ind0_var_1576 + 10, 1 : 2))
      END IF
      IF (specparm1_var_1616 .LT. 0.125D0) THEN
        tau_major1_var_1590(1 : ng15) = speccomb1_var_1617 * (fac001_var_1599 * absa_var_143(ind1_var_1577, 1 : 2) + fac101_var_1600 * absa_var_143(ind1_var_1577 + 1, 1 : 2) + fac201_var_1601 * absa_var_143(ind1_var_1577 + 2, 1 : 2) + fac011_var_1602 * absa_var_143(ind1_var_1577 + 9, 1 : 2) + fac111_var_1603 * absa_var_143(ind1_var_1577 + 10, 1 : 2) + fac211_var_1604 * absa_var_143(ind1_var_1577 + 11, 1 : 2))
      ELSE IF (specparm1_var_1616 .GT. 0.875D0) THEN
        tau_major1_var_1590(1 : ng15) = speccomb1_var_1617 * (fac201_var_1601 * absa_var_143(ind1_var_1577 - 1, 1 : 2) + fac101_var_1600 * absa_var_143(ind1_var_1577, 1 : 2) + fac001_var_1599 * absa_var_143(ind1_var_1577 + 1, 1 : 2) + fac211_var_1604 * absa_var_143(ind1_var_1577 + 8, 1 : 2) + fac111_var_1603 * absa_var_143(ind1_var_1577 + 9, 1 : 2) + fac011_var_1602 * absa_var_143(ind1_var_1577 + 10, 1 : 2))
      ELSE
        tau_major1_var_1590(1 : ng15) = speccomb1_var_1617 * (fac001_var_1599 * absa_var_143(ind1_var_1577, 1 : 2) + fac101_var_1600 * absa_var_143(ind1_var_1577 + 1, 1 : 2) + fac011_var_1602 * absa_var_143(ind1_var_1577 + 9, 1 : 2) + fac111_var_1603 * absa_var_143(ind1_var_1577 + 10, 1 : 2))
      END IF
      DO ig_var_1575 = 1, 2
        tauself_var_1588 = selffac_var_1564(jl_var_1632, lay_var_1584) * (selfref_var_145(inds_var_1578, ig_var_1575) + selffrac_var_1565(jl_var_1632, lay_var_1584) * (selfref_var_145(inds_var_1578 + 1, ig_var_1575) - selfref_var_145(inds_var_1578, ig_var_1575)))
        taufor_var_1587 = forfac_var_1569(jl_var_1632, lay_var_1584) * (forref_var_146(indf_var_1579, ig_var_1575) + forfrac_var_1570(jl_var_1632, lay_var_1584) * (forref_var_146(indf_var_1579 + 1, ig_var_1575) - forref_var_146(indf_var_1579, ig_var_1575)))
        n2m1 = ka_mn2_var_144(jmn2, indm_var_1580, ig_var_1575) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_1580, ig_var_1575) - ka_mn2_var_144(jmn2, indm_var_1580, ig_var_1575))
        n2m2 = ka_mn2_var_144(jmn2, indm_var_1580 + 1, ig_var_1575) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_1580 + 1, ig_var_1575) - ka_mn2_var_144(jmn2, indm_var_1580 + 1, ig_var_1575))
        taun2_var_1591 = scalen2_var_1592 * (n2m1 + minorfrac_var_1571(jl_var_1632, lay_var_1584) * (n2m2 - n2m1))
        taug_var_1550(jl_var_1632, 136 + ig_var_1575, lay_var_1584) = tau_major_var_1589(ig_var_1575) + tau_major1_var_1590(ig_var_1575) + tauself_var_1588 + taufor_var_1587 + taun2_var_1591
        fracs_var_1567(jl_var_1632, 136 + ig_var_1575, lay_var_1584) = fracrefa_var_142(ig_var_1575, jpl_var_1583) + fpl_var_1618 * (fracrefa_var_142(ig_var_1575, jpl_var_1583 + 1) - fracrefa_var_142(ig_var_1575, jpl_var_1583))
      END DO
    END DO
  END DO
  DO ig_var_1575 = 1, 2
    DO lay_var_1584 = laytrop_max_var_1623 + 1, klev_var_1549
      DO jl_var_1632 = kidia_var_1547, kfdia_var_1548
        taug_var_1550(jl_var_1632, 136 + ig_var_1575, lay_var_1584) = 0.0D0
        fracs_var_1567(jl_var_1632, 136 + ig_var_1575, lay_var_1584) = 0.0D0
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1623 /= laytrop_min_var_1622) THEN
    DO lay_var_1584 = laytrop_min_var_1622 + 1, laytrop_max_var_1623
      ixc0_var_1629 = ixc_var_1624(lay_var_1584)
      DO ixp_var_1630 = 1, ixc0_var_1629
        jl_var_1632 = ixlow_var_1625(ixp_var_1630, lay_var_1584)
        speccomb_var_1613 = coln2o_var_1562(jl_var_1632, lay_var_1584) + rat_n2oco2(jl_var_1632, lay_var_1584) * colco2_var_1561(jl_var_1632, lay_var_1584)
        specparm_var_1612 = MIN(coln2o_var_1562(jl_var_1632, lay_var_1584) / speccomb_var_1613, 0.999999D0)
        specmult_var_1611 = 8.0D0 * (specparm_var_1612)
        js_var_1581 = 1 + INT(specmult_var_1611)
        fs_var_1610 = ((specmult_var_1611) - AINT((specmult_var_1611)))
        speccomb1_var_1617 = coln2o_var_1562(jl_var_1632, lay_var_1584) + rat_n2oco2_1(jl_var_1632, lay_var_1584) * colco2_var_1561(jl_var_1632, lay_var_1584)
        specparm1_var_1616 = MIN(coln2o_var_1562(jl_var_1632, lay_var_1584) / speccomb1_var_1617, 0.999999D0)
        specmult1_var_1615 = 8.0D0 * (specparm1_var_1616)
        js1_var_1582 = 1 + INT(specmult1_var_1615)
        fs1_var_1614 = ((specmult1_var_1615) - AINT((specmult1_var_1615)))
        speccomb_mn2 = coln2o_var_1562(jl_var_1632, lay_var_1584) + refrat_m_a_var_1586 * colco2_var_1561(jl_var_1632, lay_var_1584)
        specparm_mn2 = MIN(coln2o_var_1562(jl_var_1632, lay_var_1584) / speccomb_mn2, 0.999999D0)
        specmult_mn2 = 8.0D0 * specparm_mn2
        jmn2 = 1 + INT(specmult_mn2)
        fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))
        speccomb_planck_var_1621 = coln2o_var_1562(jl_var_1632, lay_var_1584) + refrat_planck_a_var_1585 * colco2_var_1561(jl_var_1632, lay_var_1584)
        specparm_planck_var_1620 = MIN(coln2o_var_1562(jl_var_1632, lay_var_1584) / speccomb_planck_var_1621, 0.999999D0)
        specmult_planck_var_1619 = 8.0D0 * specparm_planck_var_1620
        jpl_var_1583 = 1 + INT(specmult_planck_var_1619)
        fpl_var_1618 = ((specmult_planck_var_1619) - AINT((specmult_planck_var_1619)))
        ind0_var_1576 = ((jp_var_1556(jl_var_1632, lay_var_1584) - 1) * 5 + (jt_var_1557(jl_var_1632, lay_var_1584) - 1)) * nspa_var_216(15) + js_var_1581
        ind1_var_1577 = (jp_var_1556(jl_var_1632, lay_var_1584) * 5 + (jt1_var_1558(jl_var_1632, lay_var_1584) - 1)) * nspa_var_216(15) + js1_var_1582
        inds_var_1578 = indself_var_1566(jl_var_1632, lay_var_1584)
        indf_var_1579 = indfor_var_1568(jl_var_1632, lay_var_1584)
        indm_var_1580 = indminor_var_1572(jl_var_1632, lay_var_1584)
        scalen2_var_1592 = colbrd_var_1574(jl_var_1632, lay_var_1584) * scaleminor_var_1573(jl_var_1632, lay_var_1584)
        IF (specparm_var_1612 .LT. 0.125D0) THEN
          p_var_1605 = fs_var_1610 - 1.0D0
          p4_var_1606 = p_var_1605 ** 4
          fk0_var_1607 = p4_var_1606
          fk1_var_1608 = 1.0D0 - p_var_1605 - 2.0D0 * p4_var_1606
          fk2_var_1609 = p_var_1605 + p4_var_1606
          fac000_var_1593 = fk0_var_1607 * fac00_var_1552(jl_var_1632, lay_var_1584)
          fac100_var_1594 = fk1_var_1608 * fac00_var_1552(jl_var_1632, lay_var_1584)
          fac200_var_1595 = fk2_var_1609 * fac00_var_1552(jl_var_1632, lay_var_1584)
          fac010_var_1596 = fk0_var_1607 * fac10_var_1554(jl_var_1632, lay_var_1584)
          fac110_var_1597 = fk1_var_1608 * fac10_var_1554(jl_var_1632, lay_var_1584)
          fac210_var_1598 = fk2_var_1609 * fac10_var_1554(jl_var_1632, lay_var_1584)
        ELSE IF (specparm_var_1612 .GT. 0.875D0) THEN
          p_var_1605 = - fs_var_1610
          p4_var_1606 = p_var_1605 ** 4
          fk0_var_1607 = p4_var_1606
          fk1_var_1608 = 1.0D0 - p_var_1605 - 2.0D0 * p4_var_1606
          fk2_var_1609 = p_var_1605 + p4_var_1606
          fac000_var_1593 = fk0_var_1607 * fac00_var_1552(jl_var_1632, lay_var_1584)
          fac100_var_1594 = fk1_var_1608 * fac00_var_1552(jl_var_1632, lay_var_1584)
          fac200_var_1595 = fk2_var_1609 * fac00_var_1552(jl_var_1632, lay_var_1584)
          fac010_var_1596 = fk0_var_1607 * fac10_var_1554(jl_var_1632, lay_var_1584)
          fac110_var_1597 = fk1_var_1608 * fac10_var_1554(jl_var_1632, lay_var_1584)
          fac210_var_1598 = fk2_var_1609 * fac10_var_1554(jl_var_1632, lay_var_1584)
        ELSE
          fac000_var_1593 = (1.0D0 - fs_var_1610) * fac00_var_1552(jl_var_1632, lay_var_1584)
          fac010_var_1596 = (1.0D0 - fs_var_1610) * fac10_var_1554(jl_var_1632, lay_var_1584)
          fac100_var_1594 = fs_var_1610 * fac00_var_1552(jl_var_1632, lay_var_1584)
          fac110_var_1597 = fs_var_1610 * fac10_var_1554(jl_var_1632, lay_var_1584)
          fac200_var_1595 = 0.0D0
          fac210_var_1598 = 0.0D0
        END IF
        IF (specparm1_var_1616 .LT. 0.125D0) THEN
          p_var_1605 = fs1_var_1614 - 1.0D0
          p4_var_1606 = p_var_1605 ** 4
          fk0_var_1607 = p4_var_1606
          fk1_var_1608 = 1.0D0 - p_var_1605 - 2.0D0 * p4_var_1606
          fk2_var_1609 = p_var_1605 + p4_var_1606
          fac001_var_1599 = fk0_var_1607 * fac01_var_1553(jl_var_1632, lay_var_1584)
          fac101_var_1600 = fk1_var_1608 * fac01_var_1553(jl_var_1632, lay_var_1584)
          fac201_var_1601 = fk2_var_1609 * fac01_var_1553(jl_var_1632, lay_var_1584)
          fac011_var_1602 = fk0_var_1607 * fac11_var_1555(jl_var_1632, lay_var_1584)
          fac111_var_1603 = fk1_var_1608 * fac11_var_1555(jl_var_1632, lay_var_1584)
          fac211_var_1604 = fk2_var_1609 * fac11_var_1555(jl_var_1632, lay_var_1584)
        ELSE IF (specparm1_var_1616 .GT. 0.875D0) THEN
          p_var_1605 = - fs1_var_1614
          p4_var_1606 = p_var_1605 ** 4
          fk0_var_1607 = p4_var_1606
          fk1_var_1608 = 1.0D0 - p_var_1605 - 2.0D0 * p4_var_1606
          fk2_var_1609 = p_var_1605 + p4_var_1606
          fac001_var_1599 = fk0_var_1607 * fac01_var_1553(jl_var_1632, lay_var_1584)
          fac101_var_1600 = fk1_var_1608 * fac01_var_1553(jl_var_1632, lay_var_1584)
          fac201_var_1601 = fk2_var_1609 * fac01_var_1553(jl_var_1632, lay_var_1584)
          fac011_var_1602 = fk0_var_1607 * fac11_var_1555(jl_var_1632, lay_var_1584)
          fac111_var_1603 = fk1_var_1608 * fac11_var_1555(jl_var_1632, lay_var_1584)
          fac211_var_1604 = fk2_var_1609 * fac11_var_1555(jl_var_1632, lay_var_1584)
        ELSE
          fac001_var_1599 = (1.0D0 - fs1_var_1614) * fac01_var_1553(jl_var_1632, lay_var_1584)
          fac011_var_1602 = (1.0D0 - fs1_var_1614) * fac11_var_1555(jl_var_1632, lay_var_1584)
          fac101_var_1600 = fs1_var_1614 * fac01_var_1553(jl_var_1632, lay_var_1584)
          fac111_var_1603 = fs1_var_1614 * fac11_var_1555(jl_var_1632, lay_var_1584)
          fac201_var_1601 = 0.0D0
          fac211_var_1604 = 0.0D0
        END IF
        IF (specparm_var_1612 .LT. 0.125D0) THEN
          tau_major_var_1589(1 : ng15) = speccomb_var_1613 * (fac000_var_1593 * absa_var_143(ind0_var_1576, 1 : 2) + fac100_var_1594 * absa_var_143(ind0_var_1576 + 1, 1 : 2) + fac200_var_1595 * absa_var_143(ind0_var_1576 + 2, 1 : 2) + fac010_var_1596 * absa_var_143(ind0_var_1576 + 9, 1 : 2) + fac110_var_1597 * absa_var_143(ind0_var_1576 + 10, 1 : 2) + fac210_var_1598 * absa_var_143(ind0_var_1576 + 11, 1 : 2))
        ELSE IF (specparm_var_1612 .GT. 0.875D0) THEN
          tau_major_var_1589(1 : ng15) = speccomb_var_1613 * (fac200_var_1595 * absa_var_143(ind0_var_1576 - 1, 1 : 2) + fac100_var_1594 * absa_var_143(ind0_var_1576, 1 : 2) + fac000_var_1593 * absa_var_143(ind0_var_1576 + 1, 1 : 2) + fac210_var_1598 * absa_var_143(ind0_var_1576 + 8, 1 : 2) + fac110_var_1597 * absa_var_143(ind0_var_1576 + 9, 1 : 2) + fac010_var_1596 * absa_var_143(ind0_var_1576 + 10, 1 : 2))
        ELSE
          tau_major_var_1589(1 : ng15) = speccomb_var_1613 * (fac000_var_1593 * absa_var_143(ind0_var_1576, 1 : 2) + fac100_var_1594 * absa_var_143(ind0_var_1576 + 1, 1 : 2) + fac010_var_1596 * absa_var_143(ind0_var_1576 + 9, 1 : 2) + fac110_var_1597 * absa_var_143(ind0_var_1576 + 10, 1 : 2))
        END IF
        IF (specparm1_var_1616 .LT. 0.125D0) THEN
          tau_major1_var_1590(1 : ng15) = speccomb1_var_1617 * (fac001_var_1599 * absa_var_143(ind1_var_1577, 1 : 2) + fac101_var_1600 * absa_var_143(ind1_var_1577 + 1, 1 : 2) + fac201_var_1601 * absa_var_143(ind1_var_1577 + 2, 1 : 2) + fac011_var_1602 * absa_var_143(ind1_var_1577 + 9, 1 : 2) + fac111_var_1603 * absa_var_143(ind1_var_1577 + 10, 1 : 2) + fac211_var_1604 * absa_var_143(ind1_var_1577 + 11, 1 : 2))
        ELSE IF (specparm1_var_1616 .GT. 0.875D0) THEN
          tau_major1_var_1590(1 : ng15) = speccomb1_var_1617 * (fac201_var_1601 * absa_var_143(ind1_var_1577 - 1, 1 : 2) + fac101_var_1600 * absa_var_143(ind1_var_1577, 1 : 2) + fac001_var_1599 * absa_var_143(ind1_var_1577 + 1, 1 : 2) + fac211_var_1604 * absa_var_143(ind1_var_1577 + 8, 1 : 2) + fac111_var_1603 * absa_var_143(ind1_var_1577 + 9, 1 : 2) + fac011_var_1602 * absa_var_143(ind1_var_1577 + 10, 1 : 2))
        ELSE
          tau_major1_var_1590(1 : ng15) = speccomb1_var_1617 * (fac001_var_1599 * absa_var_143(ind1_var_1577, 1 : 2) + fac101_var_1600 * absa_var_143(ind1_var_1577 + 1, 1 : 2) + fac011_var_1602 * absa_var_143(ind1_var_1577 + 9, 1 : 2) + fac111_var_1603 * absa_var_143(ind1_var_1577 + 10, 1 : 2))
        END IF
        DO ig_var_1575 = 1, 2
          tauself_var_1588 = selffac_var_1564(jl_var_1632, lay_var_1584) * (selfref_var_145(inds_var_1578, ig_var_1575) + selffrac_var_1565(jl_var_1632, lay_var_1584) * (selfref_var_145(inds_var_1578 + 1, ig_var_1575) - selfref_var_145(inds_var_1578, ig_var_1575)))
          taufor_var_1587 = forfac_var_1569(jl_var_1632, lay_var_1584) * (forref_var_146(indf_var_1579, ig_var_1575) + forfrac_var_1570(jl_var_1632, lay_var_1584) * (forref_var_146(indf_var_1579 + 1, ig_var_1575) - forref_var_146(indf_var_1579, ig_var_1575)))
          n2m1 = ka_mn2_var_144(jmn2, indm_var_1580, ig_var_1575) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_1580, ig_var_1575) - ka_mn2_var_144(jmn2, indm_var_1580, ig_var_1575))
          n2m2 = ka_mn2_var_144(jmn2, indm_var_1580 + 1, ig_var_1575) + fmn2 * (ka_mn2_var_144(jmn2 + 1, indm_var_1580 + 1, ig_var_1575) - ka_mn2_var_144(jmn2, indm_var_1580 + 1, ig_var_1575))
          taun2_var_1591 = scalen2_var_1592 * (n2m1 + minorfrac_var_1571(jl_var_1632, lay_var_1584) * (n2m2 - n2m1))
          taug_var_1550(jl_var_1632, 136 + ig_var_1575, lay_var_1584) = tau_major_var_1589(ig_var_1575) + tau_major1_var_1590(ig_var_1575) + tauself_var_1588 + taufor_var_1587 + taun2_var_1591
          fracs_var_1567(jl_var_1632, 136 + ig_var_1575, lay_var_1584) = fracrefa_var_142(ig_var_1575, jpl_var_1583) + fpl_var_1618 * (fracrefa_var_142(ig_var_1575, jpl_var_1583 + 1) - fracrefa_var_142(ig_var_1575, jpl_var_1583))
        END DO
      END DO
      ixc0_var_1629 = kfdia_var_1548 - kidia_var_1547 + 1 - ixc0_var_1629
      DO ig_var_1575 = 1, 2
        DO ixp_var_1630 = 1, ixc0_var_1629
          jl_var_1632 = ixhigh_var_1626(ixp_var_1630, lay_var_1584)
          taug_var_1550(jl_var_1632, 136 + ig_var_1575, lay_var_1584) = 0.0D0
          fracs_var_1567(jl_var_1632, 136 + ig_var_1575, lay_var_1584) = 0.0D0
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol15
SUBROUTINE rrtm_taumol16(kidia_var_1633, kfdia_var_1634, klev_var_1635, taug_var_1636, p_tauaerl_var_1637, fac00_var_1638, fac01_var_1639, fac10_var_1640, fac11_var_1641, forfac_var_1656, forfrac_var_1657, indfor_var_1655, jp_var_1642, jt_var_1643, jt1_var_1644, oneminus_var_1645, colh2o_var_1646, colch4_var_1647, laytrop_var_1648, selffac_var_1649, selffrac_var_1650, indself_var_1651, fracs_var_1652, rat_h2och4_var_1653, rat_h2och4_1_var_1654)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng16
  USE yoerrta16, ONLY: absa_var_149, absb_var_150, forref_var_152, fracrefa_var_147, fracrefb_var_148, selfref_var_151
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1633
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1634
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1635
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1636(kidia_var_1633 : kfdia_var_1634, 140, klev_var_1635)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1637(kidia_var_1633 : kfdia_var_1634, klev_var_1635, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1638(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1639(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1640(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1641(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1642(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1643(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1644(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1645
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1646(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  REAL(KIND = 8), INTENT(IN) :: colch4_var_1647(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1648(kidia_var_1633 : kfdia_var_1634)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1649(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1650(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1651(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1652(kidia_var_1633 : kfdia_var_1634, 140, klev_var_1635)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_var_1653(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_1_var_1654(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1655(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1656(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1657(kidia_var_1633 : kfdia_var_1634, klev_var_1635)
  INTEGER(KIND = 4) :: ig_var_1658, ind0_var_1659, ind1_var_1660, inds_var_1661, indf_var_1662, js_var_1663, js1_var_1664, jpl_var_1665, lay_var_1666
  REAL(KIND = 8) :: fac000_var_1667, fac100_var_1668, fac200_var_1669, fac010_var_1670, fac110_var_1671, fac210_var_1672, fac001_var_1673, fac101_var_1674, fac201_var_1675, fac011_var_1676, fac111_var_1677, fac211_var_1678
  REAL(KIND = 8) :: p_var_1679, p4_var_1680, fk0_var_1681, fk1_var_1682, fk2_var_1683
  REAL(KIND = 8) :: refrat_planck_a_var_1684
  REAL(KIND = 8) :: taufor_var_1685, tauself_var_1686, tau_major_var_1687(2), tau_major1_var_1688(2)
  REAL(KIND = 8) :: fs_var_1689, specmult_var_1690, specparm_var_1691, speccomb_var_1692, fs1_var_1693, specmult1_var_1694, specparm1_var_1695, speccomb1_var_1696, fpl_var_1697, specmult_planck_var_1698, specparm_planck_var_1699, speccomb_planck_var_1700
  INTEGER(KIND = 4) :: laytrop_min_var_1701, laytrop_max_var_1702
  INTEGER(KIND = 4) :: ixc_var_1703(klev_var_1635), ixlow_var_1704(kfdia_var_1634, klev_var_1635), ixhigh_var_1705(kfdia_var_1634, klev_var_1635)
  INTEGER(KIND = 4) :: ich_var_1706, icl_var_1707, ixc0_var_1708, ixp_var_1709, jc_var_1710, jl_var_1711
  laytrop_min_var_1701 = MINVAL(laytrop_var_1648)
  laytrop_max_var_1702 = MAXVAL(laytrop_var_1648)
  ixlow_var_1704 = 0
  ixhigh_var_1705 = 0
  ixc_var_1703 = 0
  DO lay_var_1666 = laytrop_min_var_1701 + 1, laytrop_max_var_1702
    icl_var_1707 = 0
    ich_var_1706 = 0
    DO jc_var_1710 = kidia_var_1633, kfdia_var_1634
      IF (lay_var_1666 <= laytrop_var_1648(jc_var_1710)) THEN
        icl_var_1707 = icl_var_1707 + 1
        ixlow_var_1704(icl_var_1707, lay_var_1666) = jc_var_1710
      ELSE
        ich_var_1706 = ich_var_1706 + 1
        ixhigh_var_1705(ich_var_1706, lay_var_1666) = jc_var_1710
      END IF
    END DO
    ixc_var_1703(lay_var_1666) = icl_var_1707
  END DO
  refrat_planck_a_var_1684 = chi_mls(1, 6) / chi_mls(6, 6)
  DO lay_var_1666 = 1, laytrop_min_var_1701
    DO jl_var_1711 = kidia_var_1633, kfdia_var_1634
      speccomb_var_1692 = colh2o_var_1646(jl_var_1711, lay_var_1666) + rat_h2och4_var_1653(jl_var_1711, lay_var_1666) * colch4_var_1647(jl_var_1711, lay_var_1666)
      specparm_var_1691 = MIN(colh2o_var_1646(jl_var_1711, lay_var_1666) / speccomb_var_1692, 0.999999D0)
      specmult_var_1690 = 8.0D0 * (specparm_var_1691)
      js_var_1663 = 1 + INT(specmult_var_1690)
      fs_var_1689 = ((specmult_var_1690) - AINT((specmult_var_1690)))
      speccomb1_var_1696 = colh2o_var_1646(jl_var_1711, lay_var_1666) + rat_h2och4_1_var_1654(jl_var_1711, lay_var_1666) * colch4_var_1647(jl_var_1711, lay_var_1666)
      specparm1_var_1695 = MIN(colh2o_var_1646(jl_var_1711, lay_var_1666) / speccomb1_var_1696, 0.999999D0)
      specmult1_var_1694 = 8.0D0 * (specparm1_var_1695)
      js1_var_1664 = 1 + INT(specmult1_var_1694)
      fs1_var_1693 = ((specmult1_var_1694) - AINT((specmult1_var_1694)))
      speccomb_planck_var_1700 = colh2o_var_1646(jl_var_1711, lay_var_1666) + refrat_planck_a_var_1684 * colch4_var_1647(jl_var_1711, lay_var_1666)
      specparm_planck_var_1699 = MIN(colh2o_var_1646(jl_var_1711, lay_var_1666) / speccomb_planck_var_1700, 0.999999D0)
      specmult_planck_var_1698 = 8.0D0 * specparm_planck_var_1699
      jpl_var_1665 = 1 + INT(specmult_planck_var_1698)
      fpl_var_1697 = ((specmult_planck_var_1698) - AINT((specmult_planck_var_1698)))
      ind0_var_1659 = ((jp_var_1642(jl_var_1711, lay_var_1666) - 1) * 5 + (jt_var_1643(jl_var_1711, lay_var_1666) - 1)) * nspa_var_216(16) + js_var_1663
      ind1_var_1660 = (jp_var_1642(jl_var_1711, lay_var_1666) * 5 + (jt1_var_1644(jl_var_1711, lay_var_1666) - 1)) * nspa_var_216(16) + js1_var_1664
      inds_var_1661 = indself_var_1651(jl_var_1711, lay_var_1666)
      indf_var_1662 = indfor_var_1655(jl_var_1711, lay_var_1666)
      IF (specparm_var_1691 .LT. 0.125D0) THEN
        p_var_1679 = fs_var_1689 - 1.0D0
        p4_var_1680 = p_var_1679 ** 4
        fk0_var_1681 = p4_var_1680
        fk1_var_1682 = 1.0D0 - p_var_1679 - 2.0D0 * p4_var_1680
        fk2_var_1683 = p_var_1679 + p4_var_1680
        fac000_var_1667 = fk0_var_1681 * fac00_var_1638(jl_var_1711, lay_var_1666)
        fac100_var_1668 = fk1_var_1682 * fac00_var_1638(jl_var_1711, lay_var_1666)
        fac200_var_1669 = fk2_var_1683 * fac00_var_1638(jl_var_1711, lay_var_1666)
        fac010_var_1670 = fk0_var_1681 * fac10_var_1640(jl_var_1711, lay_var_1666)
        fac110_var_1671 = fk1_var_1682 * fac10_var_1640(jl_var_1711, lay_var_1666)
        fac210_var_1672 = fk2_var_1683 * fac10_var_1640(jl_var_1711, lay_var_1666)
      ELSE IF (specparm_var_1691 .GT. 0.875D0) THEN
        p_var_1679 = - fs_var_1689
        p4_var_1680 = p_var_1679 ** 4
        fk0_var_1681 = p4_var_1680
        fk1_var_1682 = 1.0D0 - p_var_1679 - 2.0D0 * p4_var_1680
        fk2_var_1683 = p_var_1679 + p4_var_1680
        fac000_var_1667 = fk0_var_1681 * fac00_var_1638(jl_var_1711, lay_var_1666)
        fac100_var_1668 = fk1_var_1682 * fac00_var_1638(jl_var_1711, lay_var_1666)
        fac200_var_1669 = fk2_var_1683 * fac00_var_1638(jl_var_1711, lay_var_1666)
        fac010_var_1670 = fk0_var_1681 * fac10_var_1640(jl_var_1711, lay_var_1666)
        fac110_var_1671 = fk1_var_1682 * fac10_var_1640(jl_var_1711, lay_var_1666)
        fac210_var_1672 = fk2_var_1683 * fac10_var_1640(jl_var_1711, lay_var_1666)
      ELSE
        fac000_var_1667 = (1.0D0 - fs_var_1689) * fac00_var_1638(jl_var_1711, lay_var_1666)
        fac010_var_1670 = (1.0D0 - fs_var_1689) * fac10_var_1640(jl_var_1711, lay_var_1666)
        fac100_var_1668 = fs_var_1689 * fac00_var_1638(jl_var_1711, lay_var_1666)
        fac110_var_1671 = fs_var_1689 * fac10_var_1640(jl_var_1711, lay_var_1666)
        fac200_var_1669 = 0.0D0
        fac210_var_1672 = 0.0D0
      END IF
      IF (specparm1_var_1695 .LT. 0.125D0) THEN
        p_var_1679 = fs1_var_1693 - 1.0D0
        p4_var_1680 = p_var_1679 ** 4
        fk0_var_1681 = p4_var_1680
        fk1_var_1682 = 1.0D0 - p_var_1679 - 2.0D0 * p4_var_1680
        fk2_var_1683 = p_var_1679 + p4_var_1680
        fac001_var_1673 = fk0_var_1681 * fac01_var_1639(jl_var_1711, lay_var_1666)
        fac101_var_1674 = fk1_var_1682 * fac01_var_1639(jl_var_1711, lay_var_1666)
        fac201_var_1675 = fk2_var_1683 * fac01_var_1639(jl_var_1711, lay_var_1666)
        fac011_var_1676 = fk0_var_1681 * fac11_var_1641(jl_var_1711, lay_var_1666)
        fac111_var_1677 = fk1_var_1682 * fac11_var_1641(jl_var_1711, lay_var_1666)
        fac211_var_1678 = fk2_var_1683 * fac11_var_1641(jl_var_1711, lay_var_1666)
      ELSE IF (specparm1_var_1695 .GT. 0.875D0) THEN
        p_var_1679 = - fs1_var_1693
        p4_var_1680 = p_var_1679 ** 4
        fk0_var_1681 = p4_var_1680
        fk1_var_1682 = 1.0D0 - p_var_1679 - 2.0D0 * p4_var_1680
        fk2_var_1683 = p_var_1679 + p4_var_1680
        fac001_var_1673 = fk0_var_1681 * fac01_var_1639(jl_var_1711, lay_var_1666)
        fac101_var_1674 = fk1_var_1682 * fac01_var_1639(jl_var_1711, lay_var_1666)
        fac201_var_1675 = fk2_var_1683 * fac01_var_1639(jl_var_1711, lay_var_1666)
        fac011_var_1676 = fk0_var_1681 * fac11_var_1641(jl_var_1711, lay_var_1666)
        fac111_var_1677 = fk1_var_1682 * fac11_var_1641(jl_var_1711, lay_var_1666)
        fac211_var_1678 = fk2_var_1683 * fac11_var_1641(jl_var_1711, lay_var_1666)
      ELSE
        fac001_var_1673 = (1.0D0 - fs1_var_1693) * fac01_var_1639(jl_var_1711, lay_var_1666)
        fac011_var_1676 = (1.0D0 - fs1_var_1693) * fac11_var_1641(jl_var_1711, lay_var_1666)
        fac101_var_1674 = fs1_var_1693 * fac01_var_1639(jl_var_1711, lay_var_1666)
        fac111_var_1677 = fs1_var_1693 * fac11_var_1641(jl_var_1711, lay_var_1666)
        fac201_var_1675 = 0.0D0
        fac211_var_1678 = 0.0D0
      END IF
      IF (specparm_var_1691 .LT. 0.125D0) THEN
        tau_major_var_1687(1 : ng16) = speccomb_var_1692 * (fac000_var_1667 * absa_var_149(ind0_var_1659, 1 : 2) + fac100_var_1668 * absa_var_149(ind0_var_1659 + 1, 1 : 2) + fac200_var_1669 * absa_var_149(ind0_var_1659 + 2, 1 : 2) + fac010_var_1670 * absa_var_149(ind0_var_1659 + 9, 1 : 2) + fac110_var_1671 * absa_var_149(ind0_var_1659 + 10, 1 : 2) + fac210_var_1672 * absa_var_149(ind0_var_1659 + 11, 1 : 2))
      ELSE IF (specparm_var_1691 .GT. 0.875D0) THEN
        tau_major_var_1687(1 : ng16) = speccomb_var_1692 * (fac200_var_1669 * absa_var_149(ind0_var_1659 - 1, 1 : 2) + fac100_var_1668 * absa_var_149(ind0_var_1659, 1 : 2) + fac000_var_1667 * absa_var_149(ind0_var_1659 + 1, 1 : 2) + fac210_var_1672 * absa_var_149(ind0_var_1659 + 8, 1 : 2) + fac110_var_1671 * absa_var_149(ind0_var_1659 + 9, 1 : 2) + fac010_var_1670 * absa_var_149(ind0_var_1659 + 10, 1 : 2))
      ELSE
        tau_major_var_1687(1 : ng16) = speccomb_var_1692 * (fac000_var_1667 * absa_var_149(ind0_var_1659, 1 : 2) + fac100_var_1668 * absa_var_149(ind0_var_1659 + 1, 1 : 2) + fac010_var_1670 * absa_var_149(ind0_var_1659 + 9, 1 : 2) + fac110_var_1671 * absa_var_149(ind0_var_1659 + 10, 1 : 2))
      END IF
      IF (specparm1_var_1695 .LT. 0.125D0) THEN
        tau_major1_var_1688(1 : ng16) = speccomb1_var_1696 * (fac001_var_1673 * absa_var_149(ind1_var_1660, 1 : 2) + fac101_var_1674 * absa_var_149(ind1_var_1660 + 1, 1 : 2) + fac201_var_1675 * absa_var_149(ind1_var_1660 + 2, 1 : 2) + fac011_var_1676 * absa_var_149(ind1_var_1660 + 9, 1 : 2) + fac111_var_1677 * absa_var_149(ind1_var_1660 + 10, 1 : 2) + fac211_var_1678 * absa_var_149(ind1_var_1660 + 11, 1 : 2))
      ELSE IF (specparm1_var_1695 .GT. 0.875D0) THEN
        tau_major1_var_1688(1 : ng16) = speccomb1_var_1696 * (fac201_var_1675 * absa_var_149(ind1_var_1660 - 1, 1 : 2) + fac101_var_1674 * absa_var_149(ind1_var_1660, 1 : 2) + fac001_var_1673 * absa_var_149(ind1_var_1660 + 1, 1 : 2) + fac211_var_1678 * absa_var_149(ind1_var_1660 + 8, 1 : 2) + fac111_var_1677 * absa_var_149(ind1_var_1660 + 9, 1 : 2) + fac011_var_1676 * absa_var_149(ind1_var_1660 + 10, 1 : 2))
      ELSE
        tau_major1_var_1688(1 : ng16) = speccomb1_var_1696 * (fac001_var_1673 * absa_var_149(ind1_var_1660, 1 : 2) + fac101_var_1674 * absa_var_149(ind1_var_1660 + 1, 1 : 2) + fac011_var_1676 * absa_var_149(ind1_var_1660 + 9, 1 : 2) + fac111_var_1677 * absa_var_149(ind1_var_1660 + 10, 1 : 2))
      END IF
      DO ig_var_1658 = 1, 2
        tauself_var_1686 = selffac_var_1649(jl_var_1711, lay_var_1666) * (selfref_var_151(inds_var_1661, ig_var_1658) + selffrac_var_1650(jl_var_1711, lay_var_1666) * (selfref_var_151(inds_var_1661 + 1, ig_var_1658) - selfref_var_151(inds_var_1661, ig_var_1658)))
        taufor_var_1685 = forfac_var_1656(jl_var_1711, lay_var_1666) * (forref_var_152(indf_var_1662, ig_var_1658) + forfrac_var_1657(jl_var_1711, lay_var_1666) * (forref_var_152(indf_var_1662 + 1, ig_var_1658) - forref_var_152(indf_var_1662, ig_var_1658)))
        taug_var_1636(jl_var_1711, 138 + ig_var_1658, lay_var_1666) = tau_major_var_1687(ig_var_1658) + tau_major1_var_1688(ig_var_1658) + tauself_var_1686 + taufor_var_1685
        fracs_var_1652(jl_var_1711, 138 + ig_var_1658, lay_var_1666) = fracrefa_var_147(ig_var_1658, jpl_var_1665) + fpl_var_1697 * (fracrefa_var_147(ig_var_1658, jpl_var_1665 + 1) - fracrefa_var_147(ig_var_1658, jpl_var_1665))
      END DO
    END DO
  END DO
  DO lay_var_1666 = laytrop_max_var_1702 + 1, klev_var_1635
    DO jl_var_1711 = kidia_var_1633, kfdia_var_1634
      ind0_var_1659 = ((jp_var_1642(jl_var_1711, lay_var_1666) - 13) * 5 + (jt_var_1643(jl_var_1711, lay_var_1666) - 1)) * nspb_var_217(16) + 1
      ind1_var_1660 = ((jp_var_1642(jl_var_1711, lay_var_1666) - 12) * 5 + (jt1_var_1644(jl_var_1711, lay_var_1666) - 1)) * nspb_var_217(16) + 1
      DO ig_var_1658 = 1, 2
        taug_var_1636(jl_var_1711, 138 + ig_var_1658, lay_var_1666) = colch4_var_1647(jl_var_1711, lay_var_1666) * (fac00_var_1638(jl_var_1711, lay_var_1666) * absb_var_150(ind0_var_1659, ig_var_1658) + fac10_var_1640(jl_var_1711, lay_var_1666) * absb_var_150(ind0_var_1659 + 1, ig_var_1658) + fac01_var_1639(jl_var_1711, lay_var_1666) * absb_var_150(ind1_var_1660, ig_var_1658) + fac11_var_1641(jl_var_1711, lay_var_1666) * absb_var_150(ind1_var_1660 + 1, ig_var_1658))
        fracs_var_1652(jl_var_1711, 138 + ig_var_1658, lay_var_1666) = fracrefb_var_148(ig_var_1658)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1702 /= laytrop_min_var_1701) THEN
    DO lay_var_1666 = laytrop_min_var_1701 + 1, laytrop_max_var_1702
      ixc0_var_1708 = ixc_var_1703(lay_var_1666)
      DO ixp_var_1709 = 1, ixc0_var_1708
        jl_var_1711 = ixlow_var_1704(ixp_var_1709, lay_var_1666)
        speccomb_var_1692 = colh2o_var_1646(jl_var_1711, lay_var_1666) + rat_h2och4_var_1653(jl_var_1711, lay_var_1666) * colch4_var_1647(jl_var_1711, lay_var_1666)
        specparm_var_1691 = MIN(colh2o_var_1646(jl_var_1711, lay_var_1666) / speccomb_var_1692, 0.999999D0)
        specmult_var_1690 = 8.0D0 * (specparm_var_1691)
        js_var_1663 = 1 + INT(specmult_var_1690)
        fs_var_1689 = ((specmult_var_1690) - AINT((specmult_var_1690)))
        speccomb1_var_1696 = colh2o_var_1646(jl_var_1711, lay_var_1666) + rat_h2och4_1_var_1654(jl_var_1711, lay_var_1666) * colch4_var_1647(jl_var_1711, lay_var_1666)
        specparm1_var_1695 = MIN(colh2o_var_1646(jl_var_1711, lay_var_1666) / speccomb1_var_1696, 0.999999D0)
        specmult1_var_1694 = 8.0D0 * (specparm1_var_1695)
        js1_var_1664 = 1 + INT(specmult1_var_1694)
        fs1_var_1693 = ((specmult1_var_1694) - AINT((specmult1_var_1694)))
        speccomb_planck_var_1700 = colh2o_var_1646(jl_var_1711, lay_var_1666) + refrat_planck_a_var_1684 * colch4_var_1647(jl_var_1711, lay_var_1666)
        specparm_planck_var_1699 = MIN(colh2o_var_1646(jl_var_1711, lay_var_1666) / speccomb_planck_var_1700, 0.999999D0)
        specmult_planck_var_1698 = 8.0D0 * specparm_planck_var_1699
        jpl_var_1665 = 1 + INT(specmult_planck_var_1698)
        fpl_var_1697 = ((specmult_planck_var_1698) - AINT((specmult_planck_var_1698)))
        ind0_var_1659 = ((jp_var_1642(jl_var_1711, lay_var_1666) - 1) * 5 + (jt_var_1643(jl_var_1711, lay_var_1666) - 1)) * nspa_var_216(16) + js_var_1663
        ind1_var_1660 = (jp_var_1642(jl_var_1711, lay_var_1666) * 5 + (jt1_var_1644(jl_var_1711, lay_var_1666) - 1)) * nspa_var_216(16) + js1_var_1664
        inds_var_1661 = indself_var_1651(jl_var_1711, lay_var_1666)
        indf_var_1662 = indfor_var_1655(jl_var_1711, lay_var_1666)
        IF (specparm_var_1691 .LT. 0.125D0) THEN
          p_var_1679 = fs_var_1689 - 1.0D0
          p4_var_1680 = p_var_1679 ** 4
          fk0_var_1681 = p4_var_1680
          fk1_var_1682 = 1.0D0 - p_var_1679 - 2.0D0 * p4_var_1680
          fk2_var_1683 = p_var_1679 + p4_var_1680
          fac000_var_1667 = fk0_var_1681 * fac00_var_1638(jl_var_1711, lay_var_1666)
          fac100_var_1668 = fk1_var_1682 * fac00_var_1638(jl_var_1711, lay_var_1666)
          fac200_var_1669 = fk2_var_1683 * fac00_var_1638(jl_var_1711, lay_var_1666)
          fac010_var_1670 = fk0_var_1681 * fac10_var_1640(jl_var_1711, lay_var_1666)
          fac110_var_1671 = fk1_var_1682 * fac10_var_1640(jl_var_1711, lay_var_1666)
          fac210_var_1672 = fk2_var_1683 * fac10_var_1640(jl_var_1711, lay_var_1666)
        ELSE IF (specparm_var_1691 .GT. 0.875D0) THEN
          p_var_1679 = - fs_var_1689
          p4_var_1680 = p_var_1679 ** 4
          fk0_var_1681 = p4_var_1680
          fk1_var_1682 = 1.0D0 - p_var_1679 - 2.0D0 * p4_var_1680
          fk2_var_1683 = p_var_1679 + p4_var_1680
          fac000_var_1667 = fk0_var_1681 * fac00_var_1638(jl_var_1711, lay_var_1666)
          fac100_var_1668 = fk1_var_1682 * fac00_var_1638(jl_var_1711, lay_var_1666)
          fac200_var_1669 = fk2_var_1683 * fac00_var_1638(jl_var_1711, lay_var_1666)
          fac010_var_1670 = fk0_var_1681 * fac10_var_1640(jl_var_1711, lay_var_1666)
          fac110_var_1671 = fk1_var_1682 * fac10_var_1640(jl_var_1711, lay_var_1666)
          fac210_var_1672 = fk2_var_1683 * fac10_var_1640(jl_var_1711, lay_var_1666)
        ELSE
          fac000_var_1667 = (1.0D0 - fs_var_1689) * fac00_var_1638(jl_var_1711, lay_var_1666)
          fac010_var_1670 = (1.0D0 - fs_var_1689) * fac10_var_1640(jl_var_1711, lay_var_1666)
          fac100_var_1668 = fs_var_1689 * fac00_var_1638(jl_var_1711, lay_var_1666)
          fac110_var_1671 = fs_var_1689 * fac10_var_1640(jl_var_1711, lay_var_1666)
          fac200_var_1669 = 0.0D0
          fac210_var_1672 = 0.0D0
        END IF
        IF (specparm1_var_1695 .LT. 0.125D0) THEN
          p_var_1679 = fs1_var_1693 - 1.0D0
          p4_var_1680 = p_var_1679 ** 4
          fk0_var_1681 = p4_var_1680
          fk1_var_1682 = 1.0D0 - p_var_1679 - 2.0D0 * p4_var_1680
          fk2_var_1683 = p_var_1679 + p4_var_1680
          fac001_var_1673 = fk0_var_1681 * fac01_var_1639(jl_var_1711, lay_var_1666)
          fac101_var_1674 = fk1_var_1682 * fac01_var_1639(jl_var_1711, lay_var_1666)
          fac201_var_1675 = fk2_var_1683 * fac01_var_1639(jl_var_1711, lay_var_1666)
          fac011_var_1676 = fk0_var_1681 * fac11_var_1641(jl_var_1711, lay_var_1666)
          fac111_var_1677 = fk1_var_1682 * fac11_var_1641(jl_var_1711, lay_var_1666)
          fac211_var_1678 = fk2_var_1683 * fac11_var_1641(jl_var_1711, lay_var_1666)
        ELSE IF (specparm1_var_1695 .GT. 0.875D0) THEN
          p_var_1679 = - fs1_var_1693
          p4_var_1680 = p_var_1679 ** 4
          fk0_var_1681 = p4_var_1680
          fk1_var_1682 = 1.0D0 - p_var_1679 - 2.0D0 * p4_var_1680
          fk2_var_1683 = p_var_1679 + p4_var_1680
          fac001_var_1673 = fk0_var_1681 * fac01_var_1639(jl_var_1711, lay_var_1666)
          fac101_var_1674 = fk1_var_1682 * fac01_var_1639(jl_var_1711, lay_var_1666)
          fac201_var_1675 = fk2_var_1683 * fac01_var_1639(jl_var_1711, lay_var_1666)
          fac011_var_1676 = fk0_var_1681 * fac11_var_1641(jl_var_1711, lay_var_1666)
          fac111_var_1677 = fk1_var_1682 * fac11_var_1641(jl_var_1711, lay_var_1666)
          fac211_var_1678 = fk2_var_1683 * fac11_var_1641(jl_var_1711, lay_var_1666)
        ELSE
          fac001_var_1673 = (1.0D0 - fs1_var_1693) * fac01_var_1639(jl_var_1711, lay_var_1666)
          fac011_var_1676 = (1.0D0 - fs1_var_1693) * fac11_var_1641(jl_var_1711, lay_var_1666)
          fac101_var_1674 = fs1_var_1693 * fac01_var_1639(jl_var_1711, lay_var_1666)
          fac111_var_1677 = fs1_var_1693 * fac11_var_1641(jl_var_1711, lay_var_1666)
          fac201_var_1675 = 0.0D0
          fac211_var_1678 = 0.0D0
        END IF
        IF (specparm_var_1691 .LT. 0.125D0) THEN
          tau_major_var_1687(1 : ng16) = speccomb_var_1692 * (fac000_var_1667 * absa_var_149(ind0_var_1659, 1 : 2) + fac100_var_1668 * absa_var_149(ind0_var_1659 + 1, 1 : 2) + fac200_var_1669 * absa_var_149(ind0_var_1659 + 2, 1 : 2) + fac010_var_1670 * absa_var_149(ind0_var_1659 + 9, 1 : 2) + fac110_var_1671 * absa_var_149(ind0_var_1659 + 10, 1 : 2) + fac210_var_1672 * absa_var_149(ind0_var_1659 + 11, 1 : 2))
        ELSE IF (specparm_var_1691 .GT. 0.875D0) THEN
          tau_major_var_1687(1 : ng16) = speccomb_var_1692 * (fac200_var_1669 * absa_var_149(ind0_var_1659 - 1, 1 : 2) + fac100_var_1668 * absa_var_149(ind0_var_1659, 1 : 2) + fac000_var_1667 * absa_var_149(ind0_var_1659 + 1, 1 : 2) + fac210_var_1672 * absa_var_149(ind0_var_1659 + 8, 1 : 2) + fac110_var_1671 * absa_var_149(ind0_var_1659 + 9, 1 : 2) + fac010_var_1670 * absa_var_149(ind0_var_1659 + 10, 1 : 2))
        ELSE
          tau_major_var_1687(1 : ng16) = speccomb_var_1692 * (fac000_var_1667 * absa_var_149(ind0_var_1659, 1 : 2) + fac100_var_1668 * absa_var_149(ind0_var_1659 + 1, 1 : 2) + fac010_var_1670 * absa_var_149(ind0_var_1659 + 9, 1 : 2) + fac110_var_1671 * absa_var_149(ind0_var_1659 + 10, 1 : 2))
        END IF
        IF (specparm1_var_1695 .LT. 0.125D0) THEN
          tau_major1_var_1688(1 : ng16) = speccomb1_var_1696 * (fac001_var_1673 * absa_var_149(ind1_var_1660, 1 : 2) + fac101_var_1674 * absa_var_149(ind1_var_1660 + 1, 1 : 2) + fac201_var_1675 * absa_var_149(ind1_var_1660 + 2, 1 : 2) + fac011_var_1676 * absa_var_149(ind1_var_1660 + 9, 1 : 2) + fac111_var_1677 * absa_var_149(ind1_var_1660 + 10, 1 : 2) + fac211_var_1678 * absa_var_149(ind1_var_1660 + 11, 1 : 2))
        ELSE IF (specparm1_var_1695 .GT. 0.875D0) THEN
          tau_major1_var_1688(1 : ng16) = speccomb1_var_1696 * (fac201_var_1675 * absa_var_149(ind1_var_1660 - 1, 1 : 2) + fac101_var_1674 * absa_var_149(ind1_var_1660, 1 : 2) + fac001_var_1673 * absa_var_149(ind1_var_1660 + 1, 1 : 2) + fac211_var_1678 * absa_var_149(ind1_var_1660 + 8, 1 : 2) + fac111_var_1677 * absa_var_149(ind1_var_1660 + 9, 1 : 2) + fac011_var_1676 * absa_var_149(ind1_var_1660 + 10, 1 : 2))
        ELSE
          tau_major1_var_1688(1 : ng16) = speccomb1_var_1696 * (fac001_var_1673 * absa_var_149(ind1_var_1660, 1 : 2) + fac101_var_1674 * absa_var_149(ind1_var_1660 + 1, 1 : 2) + fac011_var_1676 * absa_var_149(ind1_var_1660 + 9, 1 : 2) + fac111_var_1677 * absa_var_149(ind1_var_1660 + 10, 1 : 2))
        END IF
        DO ig_var_1658 = 1, 2
          tauself_var_1686 = selffac_var_1649(jl_var_1711, lay_var_1666) * (selfref_var_151(inds_var_1661, ig_var_1658) + selffrac_var_1650(jl_var_1711, lay_var_1666) * (selfref_var_151(inds_var_1661 + 1, ig_var_1658) - selfref_var_151(inds_var_1661, ig_var_1658)))
          taufor_var_1685 = forfac_var_1656(jl_var_1711, lay_var_1666) * (forref_var_152(indf_var_1662, ig_var_1658) + forfrac_var_1657(jl_var_1711, lay_var_1666) * (forref_var_152(indf_var_1662 + 1, ig_var_1658) - forref_var_152(indf_var_1662, ig_var_1658)))
          taug_var_1636(jl_var_1711, 138 + ig_var_1658, lay_var_1666) = tau_major_var_1687(ig_var_1658) + tau_major1_var_1688(ig_var_1658) + tauself_var_1686 + taufor_var_1685
          fracs_var_1652(jl_var_1711, 138 + ig_var_1658, lay_var_1666) = fracrefa_var_147(ig_var_1658, jpl_var_1665) + fpl_var_1697 * (fracrefa_var_147(ig_var_1658, jpl_var_1665 + 1) - fracrefa_var_147(ig_var_1658, jpl_var_1665))
        END DO
      END DO
      ixc0_var_1708 = kfdia_var_1634 - kidia_var_1633 + 1 - ixc0_var_1708
      DO ixp_var_1709 = 1, ixc0_var_1708
        jl_var_1711 = ixhigh_var_1705(ixp_var_1709, lay_var_1666)
        ind0_var_1659 = ((jp_var_1642(jl_var_1711, lay_var_1666) - 13) * 5 + (jt_var_1643(jl_var_1711, lay_var_1666) - 1)) * nspb_var_217(16) + 1
        ind1_var_1660 = ((jp_var_1642(jl_var_1711, lay_var_1666) - 12) * 5 + (jt1_var_1644(jl_var_1711, lay_var_1666) - 1)) * nspb_var_217(16) + 1
        DO ig_var_1658 = 1, 2
          taug_var_1636(jl_var_1711, 138 + ig_var_1658, lay_var_1666) = colch4_var_1647(jl_var_1711, lay_var_1666) * (fac00_var_1638(jl_var_1711, lay_var_1666) * absb_var_150(ind0_var_1659, ig_var_1658) + fac10_var_1640(jl_var_1711, lay_var_1666) * absb_var_150(ind0_var_1659 + 1, ig_var_1658) + fac01_var_1639(jl_var_1711, lay_var_1666) * absb_var_150(ind1_var_1660, ig_var_1658) + fac11_var_1641(jl_var_1711, lay_var_1666) * absb_var_150(ind1_var_1660 + 1, ig_var_1658))
          fracs_var_1652(jl_var_1711, 138 + ig_var_1658, lay_var_1666) = fracrefb_var_148(ig_var_1658)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol16
SUBROUTINE rrtm_taumol12(kidia_var_1712, kfdia_var_1713, klev_var_1714, taug_var_1715, p_tauaerl_var_1716, fac00_var_1717, fac01_var_1718, fac10_var_1719, fac11_var_1720, forfac_var_1736, forfrac_var_1735, indfor_var_1734, jp_var_1721, jt_var_1722, jt1_var_1723, oneminus_var_1724, colh2o_var_1725, colco2_var_1726, laytrop_var_1727, selffac_var_1728, selffrac_var_1729, indself_var_1730, fracs_var_1731, rat_h2oco2_var_1732, rat_h2oco2_1_var_1733)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrtm, ONLY: ng12
  USE yoerrta12, ONLY: absa_var_127, forref_var_129, fracrefa_var_126, selfref_var_128
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1712
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1713
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1714
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1715(kidia_var_1712 : kfdia_var_1713, 140, klev_var_1714)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1716(kidia_var_1712 : kfdia_var_1713, klev_var_1714, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1717(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1718(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1719(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1720(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1721(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1722(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1723(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1724
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1725(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1726(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1727(kidia_var_1712 : kfdia_var_1713)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1728(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1729(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1730(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1731(kidia_var_1712 : kfdia_var_1713, 140, klev_var_1714)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1732(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1733(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1734(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1735(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1736(kidia_var_1712 : kfdia_var_1713, klev_var_1714)
  REAL(KIND = 8) :: speccomb_var_1737, speccomb1_var_1738, speccomb_planck_var_1739
  INTEGER(KIND = 4) :: ind0_var_1740, ind1_var_1741, inds_var_1742, indf_var_1743
  INTEGER(KIND = 4) :: ig_var_1744, js_var_1745, lay_var_1746, js1_var_1747, jpl_var_1748
  REAL(KIND = 8) :: fs_var_1749, specmult_var_1750, specparm_var_1751, fs1_var_1752, specmult1_var_1753, specparm1_var_1754, fpl_var_1755, specmult_planck_var_1756, specparm_planck_var_1757
  REAL(KIND = 8) :: fac000_var_1758, fac100_var_1759, fac200_var_1760, fac010_var_1761, fac110_var_1762, fac210_var_1763, fac001_var_1764, fac101_var_1765, fac201_var_1766, fac011_var_1767, fac111_var_1768, fac211_var_1769
  REAL(KIND = 8) :: p_var_1770, p4_var_1771, fk0_var_1772, fk1_var_1773, fk2_var_1774
  REAL(KIND = 8) :: taufor_var_1775, tauself_var_1776, tau_major_var_1777(8), tau_major1_var_1778(8)
  REAL(KIND = 8) :: refrat_planck_a_var_1779
  INTEGER(KIND = 4) :: laytrop_min_var_1780, laytrop_max_var_1781
  INTEGER(KIND = 4) :: ixc_var_1782(klev_var_1714), ixlow_var_1783(kfdia_var_1713, klev_var_1714), ixhigh_var_1784(kfdia_var_1713, klev_var_1714)
  INTEGER(KIND = 4) :: ich_var_1785, icl_var_1786, ixc0_var_1787, ixp_var_1788, jc_var_1789, jl_var_1790
  laytrop_min_var_1780 = MINVAL(laytrop_var_1727)
  laytrop_max_var_1781 = MAXVAL(laytrop_var_1727)
  ixlow_var_1783 = 0
  ixhigh_var_1784 = 0
  ixc_var_1782 = 0
  DO lay_var_1746 = laytrop_min_var_1780 + 1, laytrop_max_var_1781
    icl_var_1786 = 0
    ich_var_1785 = 0
    DO jc_var_1789 = kidia_var_1712, kfdia_var_1713
      IF (lay_var_1746 <= laytrop_var_1727(jc_var_1789)) THEN
        icl_var_1786 = icl_var_1786 + 1
        ixlow_var_1783(icl_var_1786, lay_var_1746) = jc_var_1789
      ELSE
        ich_var_1785 = ich_var_1785 + 1
        ixhigh_var_1784(ich_var_1785, lay_var_1746) = jc_var_1789
      END IF
    END DO
    ixc_var_1782(lay_var_1746) = icl_var_1786
  END DO
  refrat_planck_a_var_1779 = chi_mls(1, 10) / chi_mls(2, 10)
  DO lay_var_1746 = 1, laytrop_min_var_1780
    DO jl_var_1790 = kidia_var_1712, kfdia_var_1713
      speccomb_var_1737 = colh2o_var_1725(jl_var_1790, lay_var_1746) + rat_h2oco2_var_1732(jl_var_1790, lay_var_1746) * colco2_var_1726(jl_var_1790, lay_var_1746)
      specparm_var_1751 = MIN(colh2o_var_1725(jl_var_1790, lay_var_1746) / speccomb_var_1737, 0.999999D0)
      specmult_var_1750 = 8.0D0 * (specparm_var_1751)
      js_var_1745 = 1 + INT(specmult_var_1750)
      fs_var_1749 = ((specmult_var_1750) - AINT((specmult_var_1750)))
      speccomb1_var_1738 = colh2o_var_1725(jl_var_1790, lay_var_1746) + rat_h2oco2_1_var_1733(jl_var_1790, lay_var_1746) * colco2_var_1726(jl_var_1790, lay_var_1746)
      specparm1_var_1754 = MIN(colh2o_var_1725(jl_var_1790, lay_var_1746) / speccomb1_var_1738, 0.999999D0)
      specmult1_var_1753 = 8.0D0 * (specparm1_var_1754)
      js1_var_1747 = 1 + INT(specmult1_var_1753)
      fs1_var_1752 = ((specmult1_var_1753) - AINT((specmult1_var_1753)))
      speccomb_planck_var_1739 = colh2o_var_1725(jl_var_1790, lay_var_1746) + refrat_planck_a_var_1779 * colco2_var_1726(jl_var_1790, lay_var_1746)
      specparm_planck_var_1757 = MIN(colh2o_var_1725(jl_var_1790, lay_var_1746) / speccomb_planck_var_1739, 0.999999D0)
      specmult_planck_var_1756 = 8.0D0 * specparm_planck_var_1757
      jpl_var_1748 = 1 + INT(specmult_planck_var_1756)
      fpl_var_1755 = ((specmult_planck_var_1756) - AINT((specmult_planck_var_1756)))
      ind0_var_1740 = ((jp_var_1721(jl_var_1790, lay_var_1746) - 1) * 5 + (jt_var_1722(jl_var_1790, lay_var_1746) - 1)) * nspa_var_216(12) + js_var_1745
      ind1_var_1741 = (jp_var_1721(jl_var_1790, lay_var_1746) * 5 + (jt1_var_1723(jl_var_1790, lay_var_1746) - 1)) * nspa_var_216(12) + js1_var_1747
      inds_var_1742 = indself_var_1730(jl_var_1790, lay_var_1746)
      indf_var_1743 = indfor_var_1734(jl_var_1790, lay_var_1746)
      IF (specparm_var_1751 .LT. 0.125D0) THEN
        p_var_1770 = fs_var_1749 - 1.0D0
        p4_var_1771 = p_var_1770 ** 4
        fk0_var_1772 = p4_var_1771
        fk1_var_1773 = 1.0D0 - p_var_1770 - 2.0D0 * p4_var_1771
        fk2_var_1774 = p_var_1770 + p4_var_1771
        fac000_var_1758 = fk0_var_1772 * fac00_var_1717(jl_var_1790, lay_var_1746)
        fac100_var_1759 = fk1_var_1773 * fac00_var_1717(jl_var_1790, lay_var_1746)
        fac200_var_1760 = fk2_var_1774 * fac00_var_1717(jl_var_1790, lay_var_1746)
        fac010_var_1761 = fk0_var_1772 * fac10_var_1719(jl_var_1790, lay_var_1746)
        fac110_var_1762 = fk1_var_1773 * fac10_var_1719(jl_var_1790, lay_var_1746)
        fac210_var_1763 = fk2_var_1774 * fac10_var_1719(jl_var_1790, lay_var_1746)
      ELSE IF (specparm_var_1751 .GT. 0.875D0) THEN
        p_var_1770 = - fs_var_1749
        p4_var_1771 = p_var_1770 ** 4
        fk0_var_1772 = p4_var_1771
        fk1_var_1773 = 1.0D0 - p_var_1770 - 2.0D0 * p4_var_1771
        fk2_var_1774 = p_var_1770 + p4_var_1771
        fac000_var_1758 = fk0_var_1772 * fac00_var_1717(jl_var_1790, lay_var_1746)
        fac100_var_1759 = fk1_var_1773 * fac00_var_1717(jl_var_1790, lay_var_1746)
        fac200_var_1760 = fk2_var_1774 * fac00_var_1717(jl_var_1790, lay_var_1746)
        fac010_var_1761 = fk0_var_1772 * fac10_var_1719(jl_var_1790, lay_var_1746)
        fac110_var_1762 = fk1_var_1773 * fac10_var_1719(jl_var_1790, lay_var_1746)
        fac210_var_1763 = fk2_var_1774 * fac10_var_1719(jl_var_1790, lay_var_1746)
      ELSE
        fac000_var_1758 = (1.0D0 - fs_var_1749) * fac00_var_1717(jl_var_1790, lay_var_1746)
        fac010_var_1761 = (1.0D0 - fs_var_1749) * fac10_var_1719(jl_var_1790, lay_var_1746)
        fac100_var_1759 = fs_var_1749 * fac00_var_1717(jl_var_1790, lay_var_1746)
        fac110_var_1762 = fs_var_1749 * fac10_var_1719(jl_var_1790, lay_var_1746)
        fac200_var_1760 = 0.0D0
        fac210_var_1763 = 0.0D0
      END IF
      IF (specparm1_var_1754 .LT. 0.125D0) THEN
        p_var_1770 = fs1_var_1752 - 1.0D0
        p4_var_1771 = p_var_1770 ** 4
        fk0_var_1772 = p4_var_1771
        fk1_var_1773 = 1.0D0 - p_var_1770 - 2.0D0 * p4_var_1771
        fk2_var_1774 = p_var_1770 + p4_var_1771
        fac001_var_1764 = fk0_var_1772 * fac01_var_1718(jl_var_1790, lay_var_1746)
        fac101_var_1765 = fk1_var_1773 * fac01_var_1718(jl_var_1790, lay_var_1746)
        fac201_var_1766 = fk2_var_1774 * fac01_var_1718(jl_var_1790, lay_var_1746)
        fac011_var_1767 = fk0_var_1772 * fac11_var_1720(jl_var_1790, lay_var_1746)
        fac111_var_1768 = fk1_var_1773 * fac11_var_1720(jl_var_1790, lay_var_1746)
        fac211_var_1769 = fk2_var_1774 * fac11_var_1720(jl_var_1790, lay_var_1746)
      ELSE IF (specparm1_var_1754 .GT. 0.875D0) THEN
        p_var_1770 = - fs1_var_1752
        p4_var_1771 = p_var_1770 ** 4
        fk0_var_1772 = p4_var_1771
        fk1_var_1773 = 1.0D0 - p_var_1770 - 2.0D0 * p4_var_1771
        fk2_var_1774 = p_var_1770 + p4_var_1771
        fac001_var_1764 = fk0_var_1772 * fac01_var_1718(jl_var_1790, lay_var_1746)
        fac101_var_1765 = fk1_var_1773 * fac01_var_1718(jl_var_1790, lay_var_1746)
        fac201_var_1766 = fk2_var_1774 * fac01_var_1718(jl_var_1790, lay_var_1746)
        fac011_var_1767 = fk0_var_1772 * fac11_var_1720(jl_var_1790, lay_var_1746)
        fac111_var_1768 = fk1_var_1773 * fac11_var_1720(jl_var_1790, lay_var_1746)
        fac211_var_1769 = fk2_var_1774 * fac11_var_1720(jl_var_1790, lay_var_1746)
      ELSE
        fac001_var_1764 = (1.0D0 - fs1_var_1752) * fac01_var_1718(jl_var_1790, lay_var_1746)
        fac011_var_1767 = (1.0D0 - fs1_var_1752) * fac11_var_1720(jl_var_1790, lay_var_1746)
        fac101_var_1765 = fs1_var_1752 * fac01_var_1718(jl_var_1790, lay_var_1746)
        fac111_var_1768 = fs1_var_1752 * fac11_var_1720(jl_var_1790, lay_var_1746)
        fac201_var_1766 = 0.0D0
        fac211_var_1769 = 0.0D0
      END IF
      IF (specparm_var_1751 .LT. 0.125D0) THEN
        tau_major_var_1777(1 : ng12) = speccomb_var_1737 * (fac000_var_1758 * absa_var_127(ind0_var_1740, 1 : 8) + fac100_var_1759 * absa_var_127(ind0_var_1740 + 1, 1 : 8) + fac200_var_1760 * absa_var_127(ind0_var_1740 + 2, 1 : 8) + fac010_var_1761 * absa_var_127(ind0_var_1740 + 9, 1 : 8) + fac110_var_1762 * absa_var_127(ind0_var_1740 + 10, 1 : 8) + fac210_var_1763 * absa_var_127(ind0_var_1740 + 11, 1 : 8))
      ELSE IF (specparm_var_1751 .GT. 0.875D0) THEN
        tau_major_var_1777(1 : ng12) = speccomb_var_1737 * (fac200_var_1760 * absa_var_127(ind0_var_1740 - 1, 1 : 8) + fac100_var_1759 * absa_var_127(ind0_var_1740, 1 : 8) + fac000_var_1758 * absa_var_127(ind0_var_1740 + 1, 1 : 8) + fac210_var_1763 * absa_var_127(ind0_var_1740 + 8, 1 : 8) + fac110_var_1762 * absa_var_127(ind0_var_1740 + 9, 1 : 8) + fac010_var_1761 * absa_var_127(ind0_var_1740 + 10, 1 : 8))
      ELSE
        tau_major_var_1777(1 : ng12) = speccomb_var_1737 * (fac000_var_1758 * absa_var_127(ind0_var_1740, 1 : 8) + fac100_var_1759 * absa_var_127(ind0_var_1740 + 1, 1 : 8) + fac010_var_1761 * absa_var_127(ind0_var_1740 + 9, 1 : 8) + fac110_var_1762 * absa_var_127(ind0_var_1740 + 10, 1 : 8))
      END IF
      IF (specparm1_var_1754 .LT. 0.125D0) THEN
        tau_major1_var_1778(1 : ng12) = speccomb1_var_1738 * (fac001_var_1764 * absa_var_127(ind1_var_1741, 1 : 8) + fac101_var_1765 * absa_var_127(ind1_var_1741 + 1, 1 : 8) + fac201_var_1766 * absa_var_127(ind1_var_1741 + 2, 1 : 8) + fac011_var_1767 * absa_var_127(ind1_var_1741 + 9, 1 : 8) + fac111_var_1768 * absa_var_127(ind1_var_1741 + 10, 1 : 8) + fac211_var_1769 * absa_var_127(ind1_var_1741 + 11, 1 : 8))
      ELSE IF (specparm1_var_1754 .GT. 0.875D0) THEN
        tau_major1_var_1778(1 : ng12) = speccomb1_var_1738 * (fac201_var_1766 * absa_var_127(ind1_var_1741 - 1, 1 : 8) + fac101_var_1765 * absa_var_127(ind1_var_1741, 1 : 8) + fac001_var_1764 * absa_var_127(ind1_var_1741 + 1, 1 : 8) + fac211_var_1769 * absa_var_127(ind1_var_1741 + 8, 1 : 8) + fac111_var_1768 * absa_var_127(ind1_var_1741 + 9, 1 : 8) + fac011_var_1767 * absa_var_127(ind1_var_1741 + 10, 1 : 8))
      ELSE
        tau_major1_var_1778(1 : ng12) = speccomb1_var_1738 * (fac001_var_1764 * absa_var_127(ind1_var_1741, 1 : 8) + fac101_var_1765 * absa_var_127(ind1_var_1741 + 1, 1 : 8) + fac011_var_1767 * absa_var_127(ind1_var_1741 + 9, 1 : 8) + fac111_var_1768 * absa_var_127(ind1_var_1741 + 10, 1 : 8))
      END IF
      DO ig_var_1744 = 1, 8
        tauself_var_1776 = selffac_var_1728(jl_var_1790, lay_var_1746) * (selfref_var_128(inds_var_1742, ig_var_1744) + selffrac_var_1729(jl_var_1790, lay_var_1746) * (selfref_var_128(inds_var_1742 + 1, ig_var_1744) - selfref_var_128(inds_var_1742, ig_var_1744)))
        taufor_var_1775 = forfac_var_1736(jl_var_1790, lay_var_1746) * (forref_var_129(indf_var_1743, ig_var_1744) + forfrac_var_1735(jl_var_1790, lay_var_1746) * (forref_var_129(indf_var_1743 + 1, ig_var_1744) - forref_var_129(indf_var_1743, ig_var_1744)))
        taug_var_1715(jl_var_1790, 122 + ig_var_1744, lay_var_1746) = tau_major_var_1777(ig_var_1744) + tau_major1_var_1778(ig_var_1744) + tauself_var_1776 + taufor_var_1775
        fracs_var_1731(jl_var_1790, 122 + ig_var_1744, lay_var_1746) = fracrefa_var_126(ig_var_1744, jpl_var_1748) + fpl_var_1755 * (fracrefa_var_126(ig_var_1744, jpl_var_1748 + 1) - fracrefa_var_126(ig_var_1744, jpl_var_1748))
      END DO
    END DO
  END DO
  DO ig_var_1744 = 1, 8
    DO lay_var_1746 = laytrop_max_var_1781 + 1, klev_var_1714
      DO jl_var_1790 = kidia_var_1712, kfdia_var_1713
        taug_var_1715(jl_var_1790, 122 + ig_var_1744, lay_var_1746) = 0.0D0
        fracs_var_1731(jl_var_1790, 122 + ig_var_1744, lay_var_1746) = 0.0D0
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1781 /= laytrop_min_var_1780) THEN
    DO lay_var_1746 = laytrop_min_var_1780 + 1, laytrop_max_var_1781
      ixc0_var_1787 = ixc_var_1782(lay_var_1746)
      DO ixp_var_1788 = 1, ixc0_var_1787
        jl_var_1790 = ixlow_var_1783(ixp_var_1788, lay_var_1746)
        speccomb_var_1737 = colh2o_var_1725(jl_var_1790, lay_var_1746) + rat_h2oco2_var_1732(jl_var_1790, lay_var_1746) * colco2_var_1726(jl_var_1790, lay_var_1746)
        specparm_var_1751 = MIN(colh2o_var_1725(jl_var_1790, lay_var_1746) / speccomb_var_1737, 0.999999D0)
        specmult_var_1750 = 8.0D0 * (specparm_var_1751)
        js_var_1745 = 1 + INT(specmult_var_1750)
        fs_var_1749 = ((specmult_var_1750) - AINT((specmult_var_1750)))
        speccomb1_var_1738 = colh2o_var_1725(jl_var_1790, lay_var_1746) + rat_h2oco2_1_var_1733(jl_var_1790, lay_var_1746) * colco2_var_1726(jl_var_1790, lay_var_1746)
        specparm1_var_1754 = MIN(colh2o_var_1725(jl_var_1790, lay_var_1746) / speccomb1_var_1738, 0.999999D0)
        specmult1_var_1753 = 8.0D0 * (specparm1_var_1754)
        js1_var_1747 = 1 + INT(specmult1_var_1753)
        fs1_var_1752 = ((specmult1_var_1753) - AINT((specmult1_var_1753)))
        speccomb_planck_var_1739 = colh2o_var_1725(jl_var_1790, lay_var_1746) + refrat_planck_a_var_1779 * colco2_var_1726(jl_var_1790, lay_var_1746)
        specparm_planck_var_1757 = MIN(colh2o_var_1725(jl_var_1790, lay_var_1746) / speccomb_planck_var_1739, 0.999999D0)
        specmult_planck_var_1756 = 8.0D0 * specparm_planck_var_1757
        jpl_var_1748 = 1 + INT(specmult_planck_var_1756)
        fpl_var_1755 = ((specmult_planck_var_1756) - AINT((specmult_planck_var_1756)))
        ind0_var_1740 = ((jp_var_1721(jl_var_1790, lay_var_1746) - 1) * 5 + (jt_var_1722(jl_var_1790, lay_var_1746) - 1)) * nspa_var_216(12) + js_var_1745
        ind1_var_1741 = (jp_var_1721(jl_var_1790, lay_var_1746) * 5 + (jt1_var_1723(jl_var_1790, lay_var_1746) - 1)) * nspa_var_216(12) + js1_var_1747
        inds_var_1742 = indself_var_1730(jl_var_1790, lay_var_1746)
        indf_var_1743 = indfor_var_1734(jl_var_1790, lay_var_1746)
        IF (specparm_var_1751 .LT. 0.125D0) THEN
          p_var_1770 = fs_var_1749 - 1.0D0
          p4_var_1771 = p_var_1770 ** 4
          fk0_var_1772 = p4_var_1771
          fk1_var_1773 = 1.0D0 - p_var_1770 - 2.0D0 * p4_var_1771
          fk2_var_1774 = p_var_1770 + p4_var_1771
          fac000_var_1758 = fk0_var_1772 * fac00_var_1717(jl_var_1790, lay_var_1746)
          fac100_var_1759 = fk1_var_1773 * fac00_var_1717(jl_var_1790, lay_var_1746)
          fac200_var_1760 = fk2_var_1774 * fac00_var_1717(jl_var_1790, lay_var_1746)
          fac010_var_1761 = fk0_var_1772 * fac10_var_1719(jl_var_1790, lay_var_1746)
          fac110_var_1762 = fk1_var_1773 * fac10_var_1719(jl_var_1790, lay_var_1746)
          fac210_var_1763 = fk2_var_1774 * fac10_var_1719(jl_var_1790, lay_var_1746)
        ELSE IF (specparm_var_1751 .GT. 0.875D0) THEN
          p_var_1770 = - fs_var_1749
          p4_var_1771 = p_var_1770 ** 4
          fk0_var_1772 = p4_var_1771
          fk1_var_1773 = 1.0D0 - p_var_1770 - 2.0D0 * p4_var_1771
          fk2_var_1774 = p_var_1770 + p4_var_1771
          fac000_var_1758 = fk0_var_1772 * fac00_var_1717(jl_var_1790, lay_var_1746)
          fac100_var_1759 = fk1_var_1773 * fac00_var_1717(jl_var_1790, lay_var_1746)
          fac200_var_1760 = fk2_var_1774 * fac00_var_1717(jl_var_1790, lay_var_1746)
          fac010_var_1761 = fk0_var_1772 * fac10_var_1719(jl_var_1790, lay_var_1746)
          fac110_var_1762 = fk1_var_1773 * fac10_var_1719(jl_var_1790, lay_var_1746)
          fac210_var_1763 = fk2_var_1774 * fac10_var_1719(jl_var_1790, lay_var_1746)
        ELSE
          fac000_var_1758 = (1.0D0 - fs_var_1749) * fac00_var_1717(jl_var_1790, lay_var_1746)
          fac010_var_1761 = (1.0D0 - fs_var_1749) * fac10_var_1719(jl_var_1790, lay_var_1746)
          fac100_var_1759 = fs_var_1749 * fac00_var_1717(jl_var_1790, lay_var_1746)
          fac110_var_1762 = fs_var_1749 * fac10_var_1719(jl_var_1790, lay_var_1746)
          fac200_var_1760 = 0.0D0
          fac210_var_1763 = 0.0D0
        END IF
        IF (specparm1_var_1754 .LT. 0.125D0) THEN
          p_var_1770 = fs1_var_1752 - 1.0D0
          p4_var_1771 = p_var_1770 ** 4
          fk0_var_1772 = p4_var_1771
          fk1_var_1773 = 1.0D0 - p_var_1770 - 2.0D0 * p4_var_1771
          fk2_var_1774 = p_var_1770 + p4_var_1771
          fac001_var_1764 = fk0_var_1772 * fac01_var_1718(jl_var_1790, lay_var_1746)
          fac101_var_1765 = fk1_var_1773 * fac01_var_1718(jl_var_1790, lay_var_1746)
          fac201_var_1766 = fk2_var_1774 * fac01_var_1718(jl_var_1790, lay_var_1746)
          fac011_var_1767 = fk0_var_1772 * fac11_var_1720(jl_var_1790, lay_var_1746)
          fac111_var_1768 = fk1_var_1773 * fac11_var_1720(jl_var_1790, lay_var_1746)
          fac211_var_1769 = fk2_var_1774 * fac11_var_1720(jl_var_1790, lay_var_1746)
        ELSE IF (specparm1_var_1754 .GT. 0.875D0) THEN
          p_var_1770 = - fs1_var_1752
          p4_var_1771 = p_var_1770 ** 4
          fk0_var_1772 = p4_var_1771
          fk1_var_1773 = 1.0D0 - p_var_1770 - 2.0D0 * p4_var_1771
          fk2_var_1774 = p_var_1770 + p4_var_1771
          fac001_var_1764 = fk0_var_1772 * fac01_var_1718(jl_var_1790, lay_var_1746)
          fac101_var_1765 = fk1_var_1773 * fac01_var_1718(jl_var_1790, lay_var_1746)
          fac201_var_1766 = fk2_var_1774 * fac01_var_1718(jl_var_1790, lay_var_1746)
          fac011_var_1767 = fk0_var_1772 * fac11_var_1720(jl_var_1790, lay_var_1746)
          fac111_var_1768 = fk1_var_1773 * fac11_var_1720(jl_var_1790, lay_var_1746)
          fac211_var_1769 = fk2_var_1774 * fac11_var_1720(jl_var_1790, lay_var_1746)
        ELSE
          fac001_var_1764 = (1.0D0 - fs1_var_1752) * fac01_var_1718(jl_var_1790, lay_var_1746)
          fac011_var_1767 = (1.0D0 - fs1_var_1752) * fac11_var_1720(jl_var_1790, lay_var_1746)
          fac101_var_1765 = fs1_var_1752 * fac01_var_1718(jl_var_1790, lay_var_1746)
          fac111_var_1768 = fs1_var_1752 * fac11_var_1720(jl_var_1790, lay_var_1746)
          fac201_var_1766 = 0.0D0
          fac211_var_1769 = 0.0D0
        END IF
        IF (specparm_var_1751 .LT. 0.125D0) THEN
          tau_major_var_1777(1 : ng12) = speccomb_var_1737 * (fac000_var_1758 * absa_var_127(ind0_var_1740, 1 : 8) + fac100_var_1759 * absa_var_127(ind0_var_1740 + 1, 1 : 8) + fac200_var_1760 * absa_var_127(ind0_var_1740 + 2, 1 : 8) + fac010_var_1761 * absa_var_127(ind0_var_1740 + 9, 1 : 8) + fac110_var_1762 * absa_var_127(ind0_var_1740 + 10, 1 : 8) + fac210_var_1763 * absa_var_127(ind0_var_1740 + 11, 1 : 8))
        ELSE IF (specparm_var_1751 .GT. 0.875D0) THEN
          tau_major_var_1777(1 : ng12) = speccomb_var_1737 * (fac200_var_1760 * absa_var_127(ind0_var_1740 - 1, 1 : 8) + fac100_var_1759 * absa_var_127(ind0_var_1740, 1 : 8) + fac000_var_1758 * absa_var_127(ind0_var_1740 + 1, 1 : 8) + fac210_var_1763 * absa_var_127(ind0_var_1740 + 8, 1 : 8) + fac110_var_1762 * absa_var_127(ind0_var_1740 + 9, 1 : 8) + fac010_var_1761 * absa_var_127(ind0_var_1740 + 10, 1 : 8))
        ELSE
          tau_major_var_1777(1 : ng12) = speccomb_var_1737 * (fac000_var_1758 * absa_var_127(ind0_var_1740, 1 : 8) + fac100_var_1759 * absa_var_127(ind0_var_1740 + 1, 1 : 8) + fac010_var_1761 * absa_var_127(ind0_var_1740 + 9, 1 : 8) + fac110_var_1762 * absa_var_127(ind0_var_1740 + 10, 1 : 8))
        END IF
        IF (specparm1_var_1754 .LT. 0.125D0) THEN
          tau_major1_var_1778(1 : ng12) = speccomb1_var_1738 * (fac001_var_1764 * absa_var_127(ind1_var_1741, 1 : 8) + fac101_var_1765 * absa_var_127(ind1_var_1741 + 1, 1 : 8) + fac201_var_1766 * absa_var_127(ind1_var_1741 + 2, 1 : 8) + fac011_var_1767 * absa_var_127(ind1_var_1741 + 9, 1 : 8) + fac111_var_1768 * absa_var_127(ind1_var_1741 + 10, 1 : 8) + fac211_var_1769 * absa_var_127(ind1_var_1741 + 11, 1 : 8))
        ELSE IF (specparm1_var_1754 .GT. 0.875D0) THEN
          tau_major1_var_1778(1 : ng12) = speccomb1_var_1738 * (fac201_var_1766 * absa_var_127(ind1_var_1741 - 1, 1 : 8) + fac101_var_1765 * absa_var_127(ind1_var_1741, 1 : 8) + fac001_var_1764 * absa_var_127(ind1_var_1741 + 1, 1 : 8) + fac211_var_1769 * absa_var_127(ind1_var_1741 + 8, 1 : 8) + fac111_var_1768 * absa_var_127(ind1_var_1741 + 9, 1 : 8) + fac011_var_1767 * absa_var_127(ind1_var_1741 + 10, 1 : 8))
        ELSE
          tau_major1_var_1778(1 : ng12) = speccomb1_var_1738 * (fac001_var_1764 * absa_var_127(ind1_var_1741, 1 : 8) + fac101_var_1765 * absa_var_127(ind1_var_1741 + 1, 1 : 8) + fac011_var_1767 * absa_var_127(ind1_var_1741 + 9, 1 : 8) + fac111_var_1768 * absa_var_127(ind1_var_1741 + 10, 1 : 8))
        END IF
        DO ig_var_1744 = 1, 8
          tauself_var_1776 = selffac_var_1728(jl_var_1790, lay_var_1746) * (selfref_var_128(inds_var_1742, ig_var_1744) + selffrac_var_1729(jl_var_1790, lay_var_1746) * (selfref_var_128(inds_var_1742 + 1, ig_var_1744) - selfref_var_128(inds_var_1742, ig_var_1744)))
          taufor_var_1775 = forfac_var_1736(jl_var_1790, lay_var_1746) * (forref_var_129(indf_var_1743, ig_var_1744) + forfrac_var_1735(jl_var_1790, lay_var_1746) * (forref_var_129(indf_var_1743 + 1, ig_var_1744) - forref_var_129(indf_var_1743, ig_var_1744)))
          taug_var_1715(jl_var_1790, 122 + ig_var_1744, lay_var_1746) = tau_major_var_1777(ig_var_1744) + tau_major1_var_1778(ig_var_1744) + tauself_var_1776 + taufor_var_1775
          fracs_var_1731(jl_var_1790, 122 + ig_var_1744, lay_var_1746) = fracrefa_var_126(ig_var_1744, jpl_var_1748) + fpl_var_1755 * (fracrefa_var_126(ig_var_1744, jpl_var_1748 + 1) - fracrefa_var_126(ig_var_1744, jpl_var_1748))
        END DO
      END DO
      ixc0_var_1787 = kfdia_var_1713 - kidia_var_1712 + 1 - ixc0_var_1787
      DO ig_var_1744 = 1, 8
        DO ixp_var_1788 = 1, ixc0_var_1787
          jl_var_1790 = ixhigh_var_1784(ixp_var_1788, lay_var_1746)
          taug_var_1715(jl_var_1790, 122 + ig_var_1744, lay_var_1746) = 0.0D0
          fracs_var_1731(jl_var_1790, 122 + ig_var_1744, lay_var_1746) = 0.0D0
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol12
SUBROUTINE rrtm_taumol13(kidia_var_1791, kfdia_var_1792, klev_var_1793, taug_var_1794, p_tauaerl_var_1795, fac00_var_1796, fac01_var_1797, fac10_var_1798, fac11_var_1799, forfac_var_1815, forfrac_var_1816, indfor_var_1814, jp_var_1800, jt_var_1801, jt1_var_1802, oneminus_var_1803, colh2o_var_1804, coln2o_var_1805, colco2_var_1806, colo3_var_1807, coldry_var_1808, laytrop_var_1809, selffac_var_1810, selffrac_var_1811, indself_var_1812, fracs_var_1813, rat_h2on2o, rat_h2on2o_1, minorfrac_var_1817, indminor_var_1818)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216
  USE yoerrtm, ONLY: ng13
  USE yoerrta13, ONLY: absa_var_132, forref_var_134, fracrefa_var_130, fracrefb_var_131, ka_mco, ka_mco2_var_135, kb_mo3, selfref_var_133
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1791
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1792
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1793
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1794(kidia_var_1791 : kfdia_var_1792, 140, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1795(kidia_var_1791 : kfdia_var_1792, klev_var_1793, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1796(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1797(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1798(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1799(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1800(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1801(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1802(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1803
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1804(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_1805(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1806(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1807(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1808(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1809(kidia_var_1791 : kfdia_var_1792)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1810(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1811(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1812(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1813(kidia_var_1791 : kfdia_var_1792, 140, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: rat_h2on2o(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: rat_h2on2o_1(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1814(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1815(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1816(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1817(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1818(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  REAL(KIND = 8) :: speccomb_var_1819, speccomb1_var_1820, speccomb_planck_var_1821, speccomb_mco2_var_1822, speccomb_mco
  INTEGER(KIND = 4) :: ind0_var_1823, ind1_var_1824, inds_var_1825, indf_var_1826, indm_var_1827
  INTEGER(KIND = 4) :: ig_var_1828, js_var_1829, lay_var_1830, js1_var_1831, jpl_var_1832, jmco2_var_1833, jmco
  REAL(KIND = 8) :: refrat_planck_a_var_1834, refrat_m_a_var_1835, refrat_m_a3
  REAL(KIND = 8) :: fac000_var_1836, fac100_var_1837, fac200_var_1838, fac010_var_1839, fac110_var_1840, fac210_var_1841, fac001_var_1842, fac101_var_1843, fac201_var_1844, fac011_var_1845, fac111_var_1846, fac211_var_1847
  REAL(KIND = 8) :: p_var_1848, p4_var_1849, fk0_var_1850, fk1_var_1851, fk2_var_1852
  REAL(KIND = 8) :: taufor_var_1853, tauself_var_1854, tau_major_var_1855(4), tau_major1_var_1856(4), co2m1_var_1857, co2m2_var_1858, absco2_var_1859
  REAL(KIND = 8) :: com1, com2, absco, abso3_var_1860
  REAL(KIND = 8) :: chi_co2_var_1861, ratco2_var_1862, adjfac_var_1863, adjcolco2_var_1864
  REAL(KIND = 8) :: fs_var_1865, specmult_var_1866, specparm_var_1867, fs1_var_1868, specmult1_var_1869, specparm1_var_1870, fmco2_var_1871, specmult_mco2_var_1872, specparm_mco2_var_1873, fmco, specmult_mco, specparm_mco, fpl_var_1874, specmult_planck_var_1875, specparm_planck_var_1876
  REAL(KIND = 8) :: colco(kidia_var_1791 : kfdia_var_1792, klev_var_1793)
  INTEGER(KIND = 4) :: laytrop_min_var_1877, laytrop_max_var_1878
  INTEGER(KIND = 4) :: ixc_var_1879(klev_var_1793), ixlow_var_1880(kfdia_var_1792, klev_var_1793), ixhigh_var_1881(kfdia_var_1792, klev_var_1793)
  INTEGER(KIND = 4) :: ich_var_1882, icl_var_1883, ixc0_var_1884, ixp_var_1885, jc_var_1886, jl_var_1887
  DO lay_var_1830 = 1, klev_var_1793
    DO jc_var_1886 = kidia_var_1791, kfdia_var_1792
      colco(jc_var_1886, lay_var_1830) = 0.0D0
    END DO
  END DO
  laytrop_min_var_1877 = MINVAL(laytrop_var_1809)
  laytrop_max_var_1878 = MAXVAL(laytrop_var_1809)
  ixlow_var_1880 = 0
  ixhigh_var_1881 = 0
  ixc_var_1879 = 0
  DO lay_var_1830 = laytrop_min_var_1877 + 1, laytrop_max_var_1878
    icl_var_1883 = 0
    ich_var_1882 = 0
    DO jc_var_1886 = kidia_var_1791, kfdia_var_1792
      IF (lay_var_1830 <= laytrop_var_1809(jc_var_1886)) THEN
        icl_var_1883 = icl_var_1883 + 1
        ixlow_var_1880(icl_var_1883, lay_var_1830) = jc_var_1886
      ELSE
        ich_var_1882 = ich_var_1882 + 1
        ixhigh_var_1881(ich_var_1882, lay_var_1830) = jc_var_1886
      END IF
    END DO
    ixc_var_1879(lay_var_1830) = icl_var_1883
  END DO
  refrat_planck_a_var_1834 = chi_mls(1, 5) / chi_mls(4, 5)
  refrat_m_a_var_1835 = chi_mls(1, 1) / chi_mls(4, 1)
  refrat_m_a3 = chi_mls(1, 3) / chi_mls(4, 3)
  DO lay_var_1830 = 1, laytrop_min_var_1877
    DO jl_var_1887 = kidia_var_1791, kfdia_var_1792
      speccomb_var_1819 = colh2o_var_1804(jl_var_1887, lay_var_1830) + rat_h2on2o(jl_var_1887, lay_var_1830) * coln2o_var_1805(jl_var_1887, lay_var_1830)
      specparm_var_1867 = MIN(colh2o_var_1804(jl_var_1887, lay_var_1830) / speccomb_var_1819, 0.999999D0)
      specmult_var_1866 = 8.0D0 * (specparm_var_1867)
      js_var_1829 = 1 + INT(specmult_var_1866)
      fs_var_1865 = ((specmult_var_1866) - AINT((specmult_var_1866)))
      speccomb1_var_1820 = colh2o_var_1804(jl_var_1887, lay_var_1830) + rat_h2on2o_1(jl_var_1887, lay_var_1830) * coln2o_var_1805(jl_var_1887, lay_var_1830)
      specparm1_var_1870 = MIN(colh2o_var_1804(jl_var_1887, lay_var_1830) / speccomb1_var_1820, 0.999999D0)
      specmult1_var_1869 = 8.0D0 * (specparm1_var_1870)
      js1_var_1831 = 1 + INT(specmult1_var_1869)
      fs1_var_1868 = ((specmult1_var_1869) - AINT((specmult1_var_1869)))
      speccomb_mco2_var_1822 = colh2o_var_1804(jl_var_1887, lay_var_1830) + refrat_m_a_var_1835 * coln2o_var_1805(jl_var_1887, lay_var_1830)
      specparm_mco2_var_1873 = MIN(colh2o_var_1804(jl_var_1887, lay_var_1830) / speccomb_mco2_var_1822, 0.999999D0)
      specmult_mco2_var_1872 = 8.0D0 * specparm_mco2_var_1873
      jmco2_var_1833 = 1 + INT(specmult_mco2_var_1872)
      fmco2_var_1871 = ((specmult_mco2_var_1872) - AINT((specmult_mco2_var_1872)))
      chi_co2_var_1861 = colco2_var_1806(jl_var_1887, lay_var_1830) / (coldry_var_1808(jl_var_1887, lay_var_1830))
      ratco2_var_1862 = 1D+20 * chi_co2_var_1861 / 0.000355D0
      IF (ratco2_var_1862 .GT. 3.0D0) THEN
        adjfac_var_1863 = 2.0D0 + (ratco2_var_1862 - 2.0D0) ** 0.68D0
        adjcolco2_var_1864 = adjfac_var_1863 * 0.000355D0 * coldry_var_1808(jl_var_1887, lay_var_1830) * 1D-20
      ELSE
        adjcolco2_var_1864 = colco2_var_1806(jl_var_1887, lay_var_1830)
      END IF
      speccomb_mco = colh2o_var_1804(jl_var_1887, lay_var_1830) + refrat_m_a3 * coln2o_var_1805(jl_var_1887, lay_var_1830)
      specparm_mco = MIN(colh2o_var_1804(jl_var_1887, lay_var_1830) / speccomb_mco, 0.999999D0)
      specmult_mco = 8.0D0 * specparm_mco
      jmco = 1 + INT(specmult_mco)
      fmco = ((specmult_mco) - AINT((specmult_mco)))
      speccomb_planck_var_1821 = colh2o_var_1804(jl_var_1887, lay_var_1830) + refrat_planck_a_var_1834 * coln2o_var_1805(jl_var_1887, lay_var_1830)
      specparm_planck_var_1876 = MIN(colh2o_var_1804(jl_var_1887, lay_var_1830) / speccomb_planck_var_1821, 0.999999D0)
      specmult_planck_var_1875 = 8.0D0 * specparm_planck_var_1876
      jpl_var_1832 = 1 + INT(specmult_planck_var_1875)
      fpl_var_1874 = ((specmult_planck_var_1875) - AINT((specmult_planck_var_1875)))
      ind0_var_1823 = ((jp_var_1800(jl_var_1887, lay_var_1830) - 1) * 5 + (jt_var_1801(jl_var_1887, lay_var_1830) - 1)) * nspa_var_216(13) + js_var_1829
      ind1_var_1824 = (jp_var_1800(jl_var_1887, lay_var_1830) * 5 + (jt1_var_1802(jl_var_1887, lay_var_1830) - 1)) * nspa_var_216(13) + js1_var_1831
      inds_var_1825 = indself_var_1812(jl_var_1887, lay_var_1830)
      indf_var_1826 = indfor_var_1814(jl_var_1887, lay_var_1830)
      indm_var_1827 = indminor_var_1818(jl_var_1887, lay_var_1830)
      IF (specparm_var_1867 .LT. 0.125D0) THEN
        p_var_1848 = fs_var_1865 - 1.0D0
        p4_var_1849 = p_var_1848 ** 4
        fk0_var_1850 = p4_var_1849
        fk1_var_1851 = 1.0D0 - p_var_1848 - 2.0D0 * p4_var_1849
        fk2_var_1852 = p_var_1848 + p4_var_1849
        fac000_var_1836 = fk0_var_1850 * fac00_var_1796(jl_var_1887, lay_var_1830)
        fac100_var_1837 = fk1_var_1851 * fac00_var_1796(jl_var_1887, lay_var_1830)
        fac200_var_1838 = fk2_var_1852 * fac00_var_1796(jl_var_1887, lay_var_1830)
        fac010_var_1839 = fk0_var_1850 * fac10_var_1798(jl_var_1887, lay_var_1830)
        fac110_var_1840 = fk1_var_1851 * fac10_var_1798(jl_var_1887, lay_var_1830)
        fac210_var_1841 = fk2_var_1852 * fac10_var_1798(jl_var_1887, lay_var_1830)
      ELSE IF (specparm_var_1867 .GT. 0.875D0) THEN
        p_var_1848 = - fs_var_1865
        p4_var_1849 = p_var_1848 ** 4
        fk0_var_1850 = p4_var_1849
        fk1_var_1851 = 1.0D0 - p_var_1848 - 2.0D0 * p4_var_1849
        fk2_var_1852 = p_var_1848 + p4_var_1849
        fac000_var_1836 = fk0_var_1850 * fac00_var_1796(jl_var_1887, lay_var_1830)
        fac100_var_1837 = fk1_var_1851 * fac00_var_1796(jl_var_1887, lay_var_1830)
        fac200_var_1838 = fk2_var_1852 * fac00_var_1796(jl_var_1887, lay_var_1830)
        fac010_var_1839 = fk0_var_1850 * fac10_var_1798(jl_var_1887, lay_var_1830)
        fac110_var_1840 = fk1_var_1851 * fac10_var_1798(jl_var_1887, lay_var_1830)
        fac210_var_1841 = fk2_var_1852 * fac10_var_1798(jl_var_1887, lay_var_1830)
      ELSE
        fac000_var_1836 = (1.0D0 - fs_var_1865) * fac00_var_1796(jl_var_1887, lay_var_1830)
        fac010_var_1839 = (1.0D0 - fs_var_1865) * fac10_var_1798(jl_var_1887, lay_var_1830)
        fac100_var_1837 = fs_var_1865 * fac00_var_1796(jl_var_1887, lay_var_1830)
        fac110_var_1840 = fs_var_1865 * fac10_var_1798(jl_var_1887, lay_var_1830)
        fac200_var_1838 = 0.0D0
        fac210_var_1841 = 0.0D0
      END IF
      IF (specparm1_var_1870 .LT. 0.125D0) THEN
        p_var_1848 = fs1_var_1868 - 1.0D0
        p4_var_1849 = p_var_1848 ** 4
        fk0_var_1850 = p4_var_1849
        fk1_var_1851 = 1.0D0 - p_var_1848 - 2.0D0 * p4_var_1849
        fk2_var_1852 = p_var_1848 + p4_var_1849
        fac001_var_1842 = fk0_var_1850 * fac01_var_1797(jl_var_1887, lay_var_1830)
        fac101_var_1843 = fk1_var_1851 * fac01_var_1797(jl_var_1887, lay_var_1830)
        fac201_var_1844 = fk2_var_1852 * fac01_var_1797(jl_var_1887, lay_var_1830)
        fac011_var_1845 = fk0_var_1850 * fac11_var_1799(jl_var_1887, lay_var_1830)
        fac111_var_1846 = fk1_var_1851 * fac11_var_1799(jl_var_1887, lay_var_1830)
        fac211_var_1847 = fk2_var_1852 * fac11_var_1799(jl_var_1887, lay_var_1830)
      ELSE IF (specparm1_var_1870 .GT. 0.875D0) THEN
        p_var_1848 = - fs1_var_1868
        p4_var_1849 = p_var_1848 ** 4
        fk0_var_1850 = p4_var_1849
        fk1_var_1851 = 1.0D0 - p_var_1848 - 2.0D0 * p4_var_1849
        fk2_var_1852 = p_var_1848 + p4_var_1849
        fac001_var_1842 = fk0_var_1850 * fac01_var_1797(jl_var_1887, lay_var_1830)
        fac101_var_1843 = fk1_var_1851 * fac01_var_1797(jl_var_1887, lay_var_1830)
        fac201_var_1844 = fk2_var_1852 * fac01_var_1797(jl_var_1887, lay_var_1830)
        fac011_var_1845 = fk0_var_1850 * fac11_var_1799(jl_var_1887, lay_var_1830)
        fac111_var_1846 = fk1_var_1851 * fac11_var_1799(jl_var_1887, lay_var_1830)
        fac211_var_1847 = fk2_var_1852 * fac11_var_1799(jl_var_1887, lay_var_1830)
      ELSE
        fac001_var_1842 = (1.0D0 - fs1_var_1868) * fac01_var_1797(jl_var_1887, lay_var_1830)
        fac011_var_1845 = (1.0D0 - fs1_var_1868) * fac11_var_1799(jl_var_1887, lay_var_1830)
        fac101_var_1843 = fs1_var_1868 * fac01_var_1797(jl_var_1887, lay_var_1830)
        fac111_var_1846 = fs1_var_1868 * fac11_var_1799(jl_var_1887, lay_var_1830)
        fac201_var_1844 = 0.0D0
        fac211_var_1847 = 0.0D0
      END IF
      IF (specparm_var_1867 .LT. 0.125D0) THEN
        tau_major_var_1855(1 : ng13) = speccomb_var_1819 * (fac000_var_1836 * absa_var_132(ind0_var_1823, 1 : 4) + fac100_var_1837 * absa_var_132(ind0_var_1823 + 1, 1 : 4) + fac200_var_1838 * absa_var_132(ind0_var_1823 + 2, 1 : 4) + fac010_var_1839 * absa_var_132(ind0_var_1823 + 9, 1 : 4) + fac110_var_1840 * absa_var_132(ind0_var_1823 + 10, 1 : 4) + fac210_var_1841 * absa_var_132(ind0_var_1823 + 11, 1 : 4))
      ELSE IF (specparm_var_1867 .GT. 0.875D0) THEN
        tau_major_var_1855(1 : ng13) = speccomb_var_1819 * (fac200_var_1838 * absa_var_132(ind0_var_1823 - 1, 1 : 4) + fac100_var_1837 * absa_var_132(ind0_var_1823, 1 : 4) + fac000_var_1836 * absa_var_132(ind0_var_1823 + 1, 1 : 4) + fac210_var_1841 * absa_var_132(ind0_var_1823 + 8, 1 : 4) + fac110_var_1840 * absa_var_132(ind0_var_1823 + 9, 1 : 4) + fac010_var_1839 * absa_var_132(ind0_var_1823 + 10, 1 : 4))
      ELSE
        tau_major_var_1855(1 : ng13) = speccomb_var_1819 * (fac000_var_1836 * absa_var_132(ind0_var_1823, 1 : 4) + fac100_var_1837 * absa_var_132(ind0_var_1823 + 1, 1 : 4) + fac010_var_1839 * absa_var_132(ind0_var_1823 + 9, 1 : 4) + fac110_var_1840 * absa_var_132(ind0_var_1823 + 10, 1 : 4))
      END IF
      IF (specparm1_var_1870 .LT. 0.125D0) THEN
        tau_major1_var_1856(1 : ng13) = speccomb1_var_1820 * (fac001_var_1842 * absa_var_132(ind1_var_1824, 1 : 4) + fac101_var_1843 * absa_var_132(ind1_var_1824 + 1, 1 : 4) + fac201_var_1844 * absa_var_132(ind1_var_1824 + 2, 1 : 4) + fac011_var_1845 * absa_var_132(ind1_var_1824 + 9, 1 : 4) + fac111_var_1846 * absa_var_132(ind1_var_1824 + 10, 1 : 4) + fac211_var_1847 * absa_var_132(ind1_var_1824 + 11, 1 : 4))
      ELSE IF (specparm1_var_1870 .GT. 0.875D0) THEN
        tau_major1_var_1856(1 : ng13) = speccomb1_var_1820 * (fac201_var_1844 * absa_var_132(ind1_var_1824 - 1, 1 : 4) + fac101_var_1843 * absa_var_132(ind1_var_1824, 1 : 4) + fac001_var_1842 * absa_var_132(ind1_var_1824 + 1, 1 : 4) + fac211_var_1847 * absa_var_132(ind1_var_1824 + 8, 1 : 4) + fac111_var_1846 * absa_var_132(ind1_var_1824 + 9, 1 : 4) + fac011_var_1845 * absa_var_132(ind1_var_1824 + 10, 1 : 4))
      ELSE
        tau_major1_var_1856(1 : ng13) = speccomb1_var_1820 * (fac001_var_1842 * absa_var_132(ind1_var_1824, 1 : 4) + fac101_var_1843 * absa_var_132(ind1_var_1824 + 1, 1 : 4) + fac011_var_1845 * absa_var_132(ind1_var_1824 + 9, 1 : 4) + fac111_var_1846 * absa_var_132(ind1_var_1824 + 10, 1 : 4))
      END IF
      DO ig_var_1828 = 1, 4
        tauself_var_1854 = selffac_var_1810(jl_var_1887, lay_var_1830) * (selfref_var_133(inds_var_1825, ig_var_1828) + selffrac_var_1811(jl_var_1887, lay_var_1830) * (selfref_var_133(inds_var_1825 + 1, ig_var_1828) - selfref_var_133(inds_var_1825, ig_var_1828)))
        taufor_var_1853 = forfac_var_1815(jl_var_1887, lay_var_1830) * (forref_var_134(indf_var_1826, ig_var_1828) + forfrac_var_1816(jl_var_1887, lay_var_1830) * (forref_var_134(indf_var_1826 + 1, ig_var_1828) - forref_var_134(indf_var_1826, ig_var_1828)))
        co2m1_var_1857 = ka_mco2_var_135(jmco2_var_1833, indm_var_1827, ig_var_1828) + fmco2_var_1871 * (ka_mco2_var_135(jmco2_var_1833 + 1, indm_var_1827, ig_var_1828) - ka_mco2_var_135(jmco2_var_1833, indm_var_1827, ig_var_1828))
        co2m2_var_1858 = ka_mco2_var_135(jmco2_var_1833, indm_var_1827 + 1, ig_var_1828) + fmco2_var_1871 * (ka_mco2_var_135(jmco2_var_1833 + 1, indm_var_1827 + 1, ig_var_1828) - ka_mco2_var_135(jmco2_var_1833, indm_var_1827 + 1, ig_var_1828))
        absco2_var_1859 = co2m1_var_1857 + minorfrac_var_1817(jl_var_1887, lay_var_1830) * (co2m2_var_1858 - co2m1_var_1857)
        com1 = ka_mco(jmco, indm_var_1827, ig_var_1828) + fmco * (ka_mco(jmco + 1, indm_var_1827, ig_var_1828) - ka_mco(jmco, indm_var_1827, ig_var_1828))
        com2 = ka_mco(jmco, indm_var_1827 + 1, ig_var_1828) + fmco * (ka_mco(jmco + 1, indm_var_1827 + 1, ig_var_1828) - ka_mco(jmco, indm_var_1827 + 1, ig_var_1828))
        absco = com1 + minorfrac_var_1817(jl_var_1887, lay_var_1830) * (com2 - com1)
        taug_var_1794(jl_var_1887, 130 + ig_var_1828, lay_var_1830) = tau_major_var_1855(ig_var_1828) + tau_major1_var_1856(ig_var_1828) + tauself_var_1854 + taufor_var_1853 + adjcolco2_var_1864 * absco2_var_1859 + colco(jl_var_1887, lay_var_1830) * absco
        fracs_var_1813(jl_var_1887, 130 + ig_var_1828, lay_var_1830) = fracrefa_var_130(ig_var_1828, jpl_var_1832) + fpl_var_1874 * (fracrefa_var_130(ig_var_1828, jpl_var_1832 + 1) - fracrefa_var_130(ig_var_1828, jpl_var_1832))
      END DO
    END DO
  END DO
  DO lay_var_1830 = laytrop_max_var_1878 + 1, klev_var_1793
    DO jl_var_1887 = kidia_var_1791, kfdia_var_1792
      indm_var_1827 = indminor_var_1818(jl_var_1887, lay_var_1830)
      DO ig_var_1828 = 1, 4
        abso3_var_1860 = kb_mo3(indm_var_1827, ig_var_1828) + minorfrac_var_1817(jl_var_1887, lay_var_1830) * (kb_mo3(indm_var_1827 + 1, ig_var_1828) - kb_mo3(indm_var_1827, ig_var_1828))
        taug_var_1794(jl_var_1887, 130 + ig_var_1828, lay_var_1830) = colo3_var_1807(jl_var_1887, lay_var_1830) * abso3_var_1860
        fracs_var_1813(jl_var_1887, 130 + ig_var_1828, lay_var_1830) = fracrefb_var_131(ig_var_1828)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1878 /= laytrop_min_var_1877) THEN
    DO lay_var_1830 = laytrop_min_var_1877 + 1, laytrop_max_var_1878
      ixc0_var_1884 = ixc_var_1879(lay_var_1830)
      DO ixp_var_1885 = 1, ixc0_var_1884
        jl_var_1887 = ixlow_var_1880(ixp_var_1885, lay_var_1830)
        speccomb_var_1819 = colh2o_var_1804(jl_var_1887, lay_var_1830) + rat_h2on2o(jl_var_1887, lay_var_1830) * coln2o_var_1805(jl_var_1887, lay_var_1830)
        specparm_var_1867 = MIN(colh2o_var_1804(jl_var_1887, lay_var_1830) / speccomb_var_1819, 0.999999D0)
        specmult_var_1866 = 8.0D0 * (specparm_var_1867)
        js_var_1829 = 1 + INT(specmult_var_1866)
        fs_var_1865 = ((specmult_var_1866) - AINT((specmult_var_1866)))
        speccomb1_var_1820 = colh2o_var_1804(jl_var_1887, lay_var_1830) + rat_h2on2o_1(jl_var_1887, lay_var_1830) * coln2o_var_1805(jl_var_1887, lay_var_1830)
        specparm1_var_1870 = MIN(colh2o_var_1804(jl_var_1887, lay_var_1830) / speccomb1_var_1820, 0.999999D0)
        specmult1_var_1869 = 8.0D0 * (specparm1_var_1870)
        js1_var_1831 = 1 + INT(specmult1_var_1869)
        fs1_var_1868 = ((specmult1_var_1869) - AINT((specmult1_var_1869)))
        speccomb_mco2_var_1822 = colh2o_var_1804(jl_var_1887, lay_var_1830) + refrat_m_a_var_1835 * coln2o_var_1805(jl_var_1887, lay_var_1830)
        specparm_mco2_var_1873 = MIN(colh2o_var_1804(jl_var_1887, lay_var_1830) / speccomb_mco2_var_1822, 0.999999D0)
        specmult_mco2_var_1872 = 8.0D0 * specparm_mco2_var_1873
        jmco2_var_1833 = 1 + INT(specmult_mco2_var_1872)
        fmco2_var_1871 = ((specmult_mco2_var_1872) - AINT((specmult_mco2_var_1872)))
        chi_co2_var_1861 = colco2_var_1806(jl_var_1887, lay_var_1830) / (coldry_var_1808(jl_var_1887, lay_var_1830))
        ratco2_var_1862 = 1D+20 * chi_co2_var_1861 / 0.000355D0
        IF (ratco2_var_1862 .GT. 3.0D0) THEN
          adjfac_var_1863 = 2.0D0 + (ratco2_var_1862 - 2.0D0) ** 0.68D0
          adjcolco2_var_1864 = adjfac_var_1863 * 0.000355D0 * coldry_var_1808(jl_var_1887, lay_var_1830) * 1D-20
        ELSE
          adjcolco2_var_1864 = colco2_var_1806(jl_var_1887, lay_var_1830)
        END IF
        speccomb_mco = colh2o_var_1804(jl_var_1887, lay_var_1830) + refrat_m_a3 * coln2o_var_1805(jl_var_1887, lay_var_1830)
        specparm_mco = MIN(colh2o_var_1804(jl_var_1887, lay_var_1830) / speccomb_mco, 0.999999D0)
        specmult_mco = 8.0D0 * specparm_mco
        jmco = 1 + INT(specmult_mco)
        fmco = ((specmult_mco) - AINT((specmult_mco)))
        speccomb_planck_var_1821 = colh2o_var_1804(jl_var_1887, lay_var_1830) + refrat_planck_a_var_1834 * coln2o_var_1805(jl_var_1887, lay_var_1830)
        specparm_planck_var_1876 = MIN(colh2o_var_1804(jl_var_1887, lay_var_1830) / speccomb_planck_var_1821, 0.999999D0)
        specmult_planck_var_1875 = 8.0D0 * specparm_planck_var_1876
        jpl_var_1832 = 1 + INT(specmult_planck_var_1875)
        fpl_var_1874 = ((specmult_planck_var_1875) - AINT((specmult_planck_var_1875)))
        ind0_var_1823 = ((jp_var_1800(jl_var_1887, lay_var_1830) - 1) * 5 + (jt_var_1801(jl_var_1887, lay_var_1830) - 1)) * nspa_var_216(13) + js_var_1829
        ind1_var_1824 = (jp_var_1800(jl_var_1887, lay_var_1830) * 5 + (jt1_var_1802(jl_var_1887, lay_var_1830) - 1)) * nspa_var_216(13) + js1_var_1831
        inds_var_1825 = indself_var_1812(jl_var_1887, lay_var_1830)
        indf_var_1826 = indfor_var_1814(jl_var_1887, lay_var_1830)
        indm_var_1827 = indminor_var_1818(jl_var_1887, lay_var_1830)
        IF (specparm_var_1867 .LT. 0.125D0) THEN
          p_var_1848 = fs_var_1865 - 1.0D0
          p4_var_1849 = p_var_1848 ** 4
          fk0_var_1850 = p4_var_1849
          fk1_var_1851 = 1.0D0 - p_var_1848 - 2.0D0 * p4_var_1849
          fk2_var_1852 = p_var_1848 + p4_var_1849
          fac000_var_1836 = fk0_var_1850 * fac00_var_1796(jl_var_1887, lay_var_1830)
          fac100_var_1837 = fk1_var_1851 * fac00_var_1796(jl_var_1887, lay_var_1830)
          fac200_var_1838 = fk2_var_1852 * fac00_var_1796(jl_var_1887, lay_var_1830)
          fac010_var_1839 = fk0_var_1850 * fac10_var_1798(jl_var_1887, lay_var_1830)
          fac110_var_1840 = fk1_var_1851 * fac10_var_1798(jl_var_1887, lay_var_1830)
          fac210_var_1841 = fk2_var_1852 * fac10_var_1798(jl_var_1887, lay_var_1830)
        ELSE IF (specparm_var_1867 .GT. 0.875D0) THEN
          p_var_1848 = - fs_var_1865
          p4_var_1849 = p_var_1848 ** 4
          fk0_var_1850 = p4_var_1849
          fk1_var_1851 = 1.0D0 - p_var_1848 - 2.0D0 * p4_var_1849
          fk2_var_1852 = p_var_1848 + p4_var_1849
          fac000_var_1836 = fk0_var_1850 * fac00_var_1796(jl_var_1887, lay_var_1830)
          fac100_var_1837 = fk1_var_1851 * fac00_var_1796(jl_var_1887, lay_var_1830)
          fac200_var_1838 = fk2_var_1852 * fac00_var_1796(jl_var_1887, lay_var_1830)
          fac010_var_1839 = fk0_var_1850 * fac10_var_1798(jl_var_1887, lay_var_1830)
          fac110_var_1840 = fk1_var_1851 * fac10_var_1798(jl_var_1887, lay_var_1830)
          fac210_var_1841 = fk2_var_1852 * fac10_var_1798(jl_var_1887, lay_var_1830)
        ELSE
          fac000_var_1836 = (1.0D0 - fs_var_1865) * fac00_var_1796(jl_var_1887, lay_var_1830)
          fac010_var_1839 = (1.0D0 - fs_var_1865) * fac10_var_1798(jl_var_1887, lay_var_1830)
          fac100_var_1837 = fs_var_1865 * fac00_var_1796(jl_var_1887, lay_var_1830)
          fac110_var_1840 = fs_var_1865 * fac10_var_1798(jl_var_1887, lay_var_1830)
          fac200_var_1838 = 0.0D0
          fac210_var_1841 = 0.0D0
        END IF
        IF (specparm1_var_1870 .LT. 0.125D0) THEN
          p_var_1848 = fs1_var_1868 - 1.0D0
          p4_var_1849 = p_var_1848 ** 4
          fk0_var_1850 = p4_var_1849
          fk1_var_1851 = 1.0D0 - p_var_1848 - 2.0D0 * p4_var_1849
          fk2_var_1852 = p_var_1848 + p4_var_1849
          fac001_var_1842 = fk0_var_1850 * fac01_var_1797(jl_var_1887, lay_var_1830)
          fac101_var_1843 = fk1_var_1851 * fac01_var_1797(jl_var_1887, lay_var_1830)
          fac201_var_1844 = fk2_var_1852 * fac01_var_1797(jl_var_1887, lay_var_1830)
          fac011_var_1845 = fk0_var_1850 * fac11_var_1799(jl_var_1887, lay_var_1830)
          fac111_var_1846 = fk1_var_1851 * fac11_var_1799(jl_var_1887, lay_var_1830)
          fac211_var_1847 = fk2_var_1852 * fac11_var_1799(jl_var_1887, lay_var_1830)
        ELSE IF (specparm1_var_1870 .GT. 0.875D0) THEN
          p_var_1848 = - fs1_var_1868
          p4_var_1849 = p_var_1848 ** 4
          fk0_var_1850 = p4_var_1849
          fk1_var_1851 = 1.0D0 - p_var_1848 - 2.0D0 * p4_var_1849
          fk2_var_1852 = p_var_1848 + p4_var_1849
          fac001_var_1842 = fk0_var_1850 * fac01_var_1797(jl_var_1887, lay_var_1830)
          fac101_var_1843 = fk1_var_1851 * fac01_var_1797(jl_var_1887, lay_var_1830)
          fac201_var_1844 = fk2_var_1852 * fac01_var_1797(jl_var_1887, lay_var_1830)
          fac011_var_1845 = fk0_var_1850 * fac11_var_1799(jl_var_1887, lay_var_1830)
          fac111_var_1846 = fk1_var_1851 * fac11_var_1799(jl_var_1887, lay_var_1830)
          fac211_var_1847 = fk2_var_1852 * fac11_var_1799(jl_var_1887, lay_var_1830)
        ELSE
          fac001_var_1842 = (1.0D0 - fs1_var_1868) * fac01_var_1797(jl_var_1887, lay_var_1830)
          fac011_var_1845 = (1.0D0 - fs1_var_1868) * fac11_var_1799(jl_var_1887, lay_var_1830)
          fac101_var_1843 = fs1_var_1868 * fac01_var_1797(jl_var_1887, lay_var_1830)
          fac111_var_1846 = fs1_var_1868 * fac11_var_1799(jl_var_1887, lay_var_1830)
          fac201_var_1844 = 0.0D0
          fac211_var_1847 = 0.0D0
        END IF
        IF (specparm_var_1867 .LT. 0.125D0) THEN
          tau_major_var_1855(1 : ng13) = speccomb_var_1819 * (fac000_var_1836 * absa_var_132(ind0_var_1823, 1 : 4) + fac100_var_1837 * absa_var_132(ind0_var_1823 + 1, 1 : 4) + fac200_var_1838 * absa_var_132(ind0_var_1823 + 2, 1 : 4) + fac010_var_1839 * absa_var_132(ind0_var_1823 + 9, 1 : 4) + fac110_var_1840 * absa_var_132(ind0_var_1823 + 10, 1 : 4) + fac210_var_1841 * absa_var_132(ind0_var_1823 + 11, 1 : 4))
        ELSE IF (specparm_var_1867 .GT. 0.875D0) THEN
          tau_major_var_1855(1 : ng13) = speccomb_var_1819 * (fac200_var_1838 * absa_var_132(ind0_var_1823 - 1, 1 : 4) + fac100_var_1837 * absa_var_132(ind0_var_1823, 1 : 4) + fac000_var_1836 * absa_var_132(ind0_var_1823 + 1, 1 : 4) + fac210_var_1841 * absa_var_132(ind0_var_1823 + 8, 1 : 4) + fac110_var_1840 * absa_var_132(ind0_var_1823 + 9, 1 : 4) + fac010_var_1839 * absa_var_132(ind0_var_1823 + 10, 1 : 4))
        ELSE
          tau_major_var_1855(1 : ng13) = speccomb_var_1819 * (fac000_var_1836 * absa_var_132(ind0_var_1823, 1 : 4) + fac100_var_1837 * absa_var_132(ind0_var_1823 + 1, 1 : 4) + fac010_var_1839 * absa_var_132(ind0_var_1823 + 9, 1 : 4) + fac110_var_1840 * absa_var_132(ind0_var_1823 + 10, 1 : 4))
        END IF
        IF (specparm1_var_1870 .LT. 0.125D0) THEN
          tau_major1_var_1856(1 : ng13) = speccomb1_var_1820 * (fac001_var_1842 * absa_var_132(ind1_var_1824, 1 : 4) + fac101_var_1843 * absa_var_132(ind1_var_1824 + 1, 1 : 4) + fac201_var_1844 * absa_var_132(ind1_var_1824 + 2, 1 : 4) + fac011_var_1845 * absa_var_132(ind1_var_1824 + 9, 1 : 4) + fac111_var_1846 * absa_var_132(ind1_var_1824 + 10, 1 : 4) + fac211_var_1847 * absa_var_132(ind1_var_1824 + 11, 1 : 4))
        ELSE IF (specparm1_var_1870 .GT. 0.875D0) THEN
          tau_major1_var_1856(1 : ng13) = speccomb1_var_1820 * (fac201_var_1844 * absa_var_132(ind1_var_1824 - 1, 1 : 4) + fac101_var_1843 * absa_var_132(ind1_var_1824, 1 : 4) + fac001_var_1842 * absa_var_132(ind1_var_1824 + 1, 1 : 4) + fac211_var_1847 * absa_var_132(ind1_var_1824 + 8, 1 : 4) + fac111_var_1846 * absa_var_132(ind1_var_1824 + 9, 1 : 4) + fac011_var_1845 * absa_var_132(ind1_var_1824 + 10, 1 : 4))
        ELSE
          tau_major1_var_1856(1 : ng13) = speccomb1_var_1820 * (fac001_var_1842 * absa_var_132(ind1_var_1824, 1 : 4) + fac101_var_1843 * absa_var_132(ind1_var_1824 + 1, 1 : 4) + fac011_var_1845 * absa_var_132(ind1_var_1824 + 9, 1 : 4) + fac111_var_1846 * absa_var_132(ind1_var_1824 + 10, 1 : 4))
        END IF
        DO ig_var_1828 = 1, 4
          tauself_var_1854 = selffac_var_1810(jl_var_1887, lay_var_1830) * (selfref_var_133(inds_var_1825, ig_var_1828) + selffrac_var_1811(jl_var_1887, lay_var_1830) * (selfref_var_133(inds_var_1825 + 1, ig_var_1828) - selfref_var_133(inds_var_1825, ig_var_1828)))
          taufor_var_1853 = forfac_var_1815(jl_var_1887, lay_var_1830) * (forref_var_134(indf_var_1826, ig_var_1828) + forfrac_var_1816(jl_var_1887, lay_var_1830) * (forref_var_134(indf_var_1826 + 1, ig_var_1828) - forref_var_134(indf_var_1826, ig_var_1828)))
          co2m1_var_1857 = ka_mco2_var_135(jmco2_var_1833, indm_var_1827, ig_var_1828) + fmco2_var_1871 * (ka_mco2_var_135(jmco2_var_1833 + 1, indm_var_1827, ig_var_1828) - ka_mco2_var_135(jmco2_var_1833, indm_var_1827, ig_var_1828))
          co2m2_var_1858 = ka_mco2_var_135(jmco2_var_1833, indm_var_1827 + 1, ig_var_1828) + fmco2_var_1871 * (ka_mco2_var_135(jmco2_var_1833 + 1, indm_var_1827 + 1, ig_var_1828) - ka_mco2_var_135(jmco2_var_1833, indm_var_1827 + 1, ig_var_1828))
          absco2_var_1859 = co2m1_var_1857 + minorfrac_var_1817(jl_var_1887, lay_var_1830) * (co2m2_var_1858 - co2m1_var_1857)
          com1 = ka_mco(jmco, indm_var_1827, ig_var_1828) + fmco * (ka_mco(jmco + 1, indm_var_1827, ig_var_1828) - ka_mco(jmco, indm_var_1827, ig_var_1828))
          com2 = ka_mco(jmco, indm_var_1827 + 1, ig_var_1828) + fmco * (ka_mco(jmco + 1, indm_var_1827 + 1, ig_var_1828) - ka_mco(jmco, indm_var_1827 + 1, ig_var_1828))
          absco = com1 + minorfrac_var_1817(jl_var_1887, lay_var_1830) * (com2 - com1)
          taug_var_1794(jl_var_1887, 130 + ig_var_1828, lay_var_1830) = tau_major_var_1855(ig_var_1828) + tau_major1_var_1856(ig_var_1828) + tauself_var_1854 + taufor_var_1853 + adjcolco2_var_1864 * absco2_var_1859 + colco(jl_var_1887, lay_var_1830) * absco
          fracs_var_1813(jl_var_1887, 130 + ig_var_1828, lay_var_1830) = fracrefa_var_130(ig_var_1828, jpl_var_1832) + fpl_var_1874 * (fracrefa_var_130(ig_var_1828, jpl_var_1832 + 1) - fracrefa_var_130(ig_var_1828, jpl_var_1832))
        END DO
      END DO
      ixc0_var_1884 = kfdia_var_1792 - kidia_var_1791 + 1 - ixc0_var_1884
      DO ixp_var_1885 = 1, ixc0_var_1884
        jl_var_1887 = ixhigh_var_1881(ixp_var_1885, lay_var_1830)
        indm_var_1827 = indminor_var_1818(jl_var_1887, lay_var_1830)
        DO ig_var_1828 = 1, 4
          abso3_var_1860 = kb_mo3(indm_var_1827, ig_var_1828) + minorfrac_var_1817(jl_var_1887, lay_var_1830) * (kb_mo3(indm_var_1827 + 1, ig_var_1828) - kb_mo3(indm_var_1827, ig_var_1828))
          taug_var_1794(jl_var_1887, 130 + ig_var_1828, lay_var_1830) = colo3_var_1807(jl_var_1887, lay_var_1830) * abso3_var_1860
          fracs_var_1813(jl_var_1887, 130 + ig_var_1828, lay_var_1830) = fracrefb_var_131(ig_var_1828)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol13
SUBROUTINE rrtm_gas_optical_depth(kidia_var_1888, kfdia_var_1889, klev_var_1890, pod_var_1891, pavel_var_1892, pcoldry_var_1893, pcolbrd, pwx_var_1894, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, poneminus_var_1902, pcolh2o_var_1903, pcolco2_var_1904, pcolo3_var_1905, pcoln2o, pcolch4_var_1906, pcolo2_var_1907, p_co2mult_var_1908, klaytrop_var_1909, klayswtch, klaylow, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, kindminor, pscaleminor, pscaleminorn2, pminorfrac, prat_h2oco2_var_1917, prat_h2oco2_1_var_1918, prat_h2oo3_var_1919, prat_h2oo3_1_var_1920, prat_h2on2o_var_1921, prat_h2on2o_1_var_1922, prat_h2och4_var_1923, prat_h2och4_1_var_1924, prat_n2oco2_var_1925, prat_n2oco2_1_var_1926, prat_o3co2_var_1927, prat_o3co2_1_var_1928)
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1888
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1889
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1890
  REAL(KIND = 8), INTENT(OUT) :: pod_var_1891(140, klev_var_1890, kidia_var_1888 : kfdia_var_1889)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1892(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pcoldry_var_1893(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pwx_var_1894(kidia_var_1888 : kfdia_var_1889, 4, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: ptauaerl(kidia_var_1888 : kfdia_var_1889, klev_var_1890, 16)
  REAL(KIND = 8), INTENT(IN) :: pfac00_var_1895(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pfac01_var_1896(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pfac10_var_1897(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pfac11_var_1898(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  INTEGER(KIND = 4), INTENT(IN) :: kjp_var_1899(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  INTEGER(KIND = 4), INTENT(IN) :: kjt_var_1900(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  INTEGER(KIND = 4), INTENT(IN) :: kjt1_var_1901(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: poneminus_var_1902
  REAL(KIND = 8), INTENT(IN) :: pcolh2o_var_1903(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pcolco2_var_1904(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pcolo3_var_1905(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pcoln2o(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pcolch4_var_1906(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pcolo2_var_1907(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: p_co2mult_var_1908(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  INTEGER(KIND = 4), INTENT(IN) :: klaytrop_var_1909(kidia_var_1888 : kfdia_var_1889)
  INTEGER(KIND = 4), INTENT(IN) :: klayswtch(kidia_var_1888 : kfdia_var_1889)
  INTEGER(KIND = 4), INTENT(IN) :: klaylow(kidia_var_1888 : kfdia_var_1889)
  REAL(KIND = 8), INTENT(IN) :: pselffac_var_1910(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pselffrac_var_1911(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  INTEGER(KIND = 4), INTENT(IN) :: kindself_var_1912(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(OUT) :: pfrac_var_1913(kidia_var_1888 : kfdia_var_1889, 140, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pforfac_var_1914(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pforfrac_var_1915(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  INTEGER(KIND = 4), INTENT(IN) :: kindfor_var_1916(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pminorfrac(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pscaleminor(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pscaleminorn2(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  INTEGER(KIND = 4), INTENT(IN) :: kindminor(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: pcolbrd(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8), INTENT(IN) :: prat_h2oco2_var_1917(kidia_var_1888 : kfdia_var_1889, klev_var_1890), prat_h2oco2_1_var_1918(kidia_var_1888 : kfdia_var_1889, klev_var_1890), prat_h2oo3_var_1919(kidia_var_1888 : kfdia_var_1889, klev_var_1890), prat_h2oo3_1_var_1920(kidia_var_1888 : kfdia_var_1889, klev_var_1890), prat_h2on2o_var_1921(kidia_var_1888 : kfdia_var_1889, klev_var_1890), prat_h2on2o_1_var_1922(kidia_var_1888 : kfdia_var_1889, klev_var_1890), prat_h2och4_var_1923(kidia_var_1888 : kfdia_var_1889, klev_var_1890), prat_h2och4_1_var_1924(kidia_var_1888 : kfdia_var_1889, klev_var_1890), prat_n2oco2_var_1925(kidia_var_1888 : kfdia_var_1889, klev_var_1890), prat_n2oco2_1_var_1926(kidia_var_1888 : kfdia_var_1889, klev_var_1890), prat_o3co2_var_1927(kidia_var_1888 : kfdia_var_1889, klev_var_1890), prat_o3co2_1_var_1928(kidia_var_1888 : kfdia_var_1889, klev_var_1890)
  REAL(KIND = 8) :: ztau(kidia_var_1888 : kfdia_var_1889, 140, klev_var_1890)
  INTEGER(KIND = 4) :: ji, jlev_var_1929
  INTEGER(KIND = 4) :: jlon_var_1930
  pfrac_var_1913(:, :, :) = 0.0D0
  CALL rrtm_taumol1(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, pavel_var_1892, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, pcolh2o_var_1903, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, pminorfrac, kindminor, pscaleminorn2, pcolbrd)
  CALL rrtm_taumol2(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, pavel_var_1892, pcoldry_var_1893, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, pcolh2o_var_1903, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913)
  CALL rrtm_taumol3(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, 0.999999D0, pcolh2o_var_1903, pcolco2_var_1904, pcoln2o, pcoldry_var_1893, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, prat_h2oco2_var_1917, prat_h2oco2_1_var_1918, pminorfrac, kindminor)
  CALL rrtm_taumol4(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, 0.999999D0, pcolh2o_var_1903, pcolco2_var_1904, pcolo3_var_1905, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, prat_h2oco2_var_1917, prat_h2oco2_1_var_1918, prat_o3co2_var_1927, prat_o3co2_1_var_1928)
  CALL rrtm_taumol5(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, pwx_var_1894, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, 0.999999D0, pcolh2o_var_1903, pcolco2_var_1904, pcolo3_var_1905, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, prat_h2oco2_var_1917, prat_h2oco2_1_var_1918, prat_o3co2_var_1927, prat_o3co2_1_var_1928, pminorfrac, kindminor)
  CALL rrtm_taumol6(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, pwx_var_1894, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, pcolh2o_var_1903, pcolco2_var_1904, pcoldry_var_1893, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, pminorfrac, kindminor)
  CALL rrtm_taumol7(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, 0.999999D0, pcolh2o_var_1903, pcolo3_var_1905, pcolco2_var_1904, pcoldry_var_1893, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, prat_h2oo3_var_1919, prat_h2oo3_1_var_1920, pminorfrac, kindminor)
  CALL rrtm_taumol8(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, pwx_var_1894, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, pcolh2o_var_1903, pcolo3_var_1905, pcoln2o, pcolco2_var_1904, pcoldry_var_1893, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, pminorfrac, kindminor)
  CALL rrtm_taumol9(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, 0.999999D0, pcolh2o_var_1903, pcoln2o, pcolch4_var_1906, pcoldry_var_1893, klaytrop_var_1909, klayswtch, klaylow, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, prat_h2och4_var_1923, prat_h2och4_1_var_1924, pminorfrac, kindminor)
  CALL rrtm_taumol10(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, pcolh2o_var_1903, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913)
  CALL rrtm_taumol11(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, pcolh2o_var_1903, pcolo2_var_1907, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, pminorfrac, kindminor, pscaleminor)
  CALL rrtm_taumol12(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, 0.999999D0, pcolh2o_var_1903, pcolco2_var_1904, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, prat_h2oco2_var_1917, prat_h2oco2_1_var_1918)
  CALL rrtm_taumol13(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, 0.999999D0, pcolh2o_var_1903, pcoln2o, pcolco2_var_1904, pcolo3_var_1905, pcoldry_var_1893, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, prat_h2on2o_var_1921, prat_h2on2o_1_var_1922, pminorfrac, kindminor)
  CALL rrtm_taumol14(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, pcolco2_var_1904, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913)
  CALL rrtm_taumol15(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, 0.999999D0, pcolh2o_var_1903, pcolco2_var_1904, pcoln2o, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, prat_n2oco2_var_1925, prat_n2oco2_1_var_1926, pminorfrac, kindminor, pscaleminor, pcolbrd)
  CALL rrtm_taumol16(kidia_var_1888, kfdia_var_1889, klev_var_1890, ztau, ptauaerl, pfac00_var_1895, pfac01_var_1896, pfac10_var_1897, pfac11_var_1898, pforfac_var_1914, pforfrac_var_1915, kindfor_var_1916, kjp_var_1899, kjt_var_1900, kjt1_var_1901, 0.999999D0, pcolh2o_var_1903, pcolch4_var_1906, klaytrop_var_1909, pselffac_var_1910, pselffrac_var_1911, kindself_var_1912, pfrac_var_1913, prat_h2och4_var_1923, prat_h2och4_1_var_1924)
  DO jlev_var_1929 = 1, klev_var_1890
    DO ji = 1, 140
      DO jlon_var_1930 = kidia_var_1888, kfdia_var_1889
        pod_var_1891(ji, jlev_var_1929, jlon_var_1930) = ztau(jlon_var_1930, ji, jlev_var_1929)
      END DO
    END DO
  END DO
END SUBROUTINE rrtm_gas_optical_depth
SUBROUTINE rrtm_taumol9(kidia_var_1931, kfdia_var_1932, klev_var_1933, taug_var_1934, p_tauaerl_var_1935, fac00_var_1936, fac01_var_1937, fac10_var_1938, fac11_var_1939, forfac_var_1958, forfrac_var_1959, indfor_var_1957, jp_var_1940, jt_var_1941, jt1_var_1942, oneminus_var_1943, colh2o_var_1944, coln2o_var_1945, colch4_var_1946, coldry_var_1947, laytrop_var_1948, k_layswtch_var_1949, k_laylow_var_1950, selffac_var_1951, selffrac_var_1952, indself_var_1953, fracs_var_1954, rat_h2och4_var_1955, rat_h2och4_1_var_1956, minorfrac_var_1960, indminor_var_1961)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrtm, ONLY: ng9
  USE yoerrta9, ONLY: absa_var_208, absb_var_209, forref_var_213, fracrefa_var_206, fracrefb_var_207, ka_mn2o_var_210, kb_mn2o_var_211, selfref_var_212
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1931
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1932
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1933
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1934(kidia_var_1931 : kfdia_var_1932, 140, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1935(kidia_var_1931 : kfdia_var_1932, klev_var_1933, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1936(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1937(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1938(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1939(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1940(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1941(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1942(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1943
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1944(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_1945(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: colch4_var_1946(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1947(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1948(kidia_var_1931 : kfdia_var_1932)
  INTEGER(KIND = 4), INTENT(IN) :: k_layswtch_var_1949(kidia_var_1931 : kfdia_var_1932)
  INTEGER(KIND = 4), INTENT(IN) :: k_laylow_var_1950(kidia_var_1931 : kfdia_var_1932)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1951(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1952(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1953(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1954(kidia_var_1931 : kfdia_var_1932, 140, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_var_1955(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_1_var_1956(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1957(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1958(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1959(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1960(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1961(kidia_var_1931 : kfdia_var_1932, klev_var_1933)
  INTEGER(KIND = 4) :: ind0_var_1962, ind1_var_1963, inds_var_1964, indf_var_1965, indm_var_1966
  INTEGER(KIND = 4) :: ig_var_1967, js_var_1968, lay_var_1969, js1_var_1970, jmn2o_var_1971, jpl_var_1972
  REAL(KIND = 8) :: speccomb_var_1973, speccomb1_var_1974, speccomb_mn2o_var_1975, speccomb_planck_var_1976
  REAL(KIND = 8) :: refrat_planck_a_var_1977, refrat_m_a_var_1978
  REAL(KIND = 8) :: fs_var_1979, specmult_var_1980, specparm_var_1981, fs1_var_1982, specmult1_var_1983, specparm1_var_1984, fmn2o_var_1985, specmult_mn2o_var_1986, specparm_mn2o_var_1987, fpl_var_1988, specmult_planck_var_1989, specparm_planck_var_1990
  REAL(KIND = 8) :: adjfac_var_1991, adjcoln2o_var_1992, ratn2o_var_1993, chi_n2o_var_1994
  REAL(KIND = 8) :: fac000_var_1995, fac100_var_1996, fac200_var_1997, fac010_var_1998, fac110_var_1999, fac210_var_2000, fac001_var_2001, fac101_var_2002, fac201_var_2003, fac011_var_2004, fac111_var_2005, fac211_var_2006
  REAL(KIND = 8) :: p_var_2007, p4_var_2008, fk0_var_2009, fk1_var_2010, fk2_var_2011
  REAL(KIND = 8) :: taufor_var_2012, tauself_var_2013, n2om1_var_2014, n2om2_var_2015, absn2o_var_2016, tau_major_var_2017(12), tau_major1_var_2018(12)
  INTEGER(KIND = 4) :: laytrop_min_var_2019, laytrop_max_var_2020
  INTEGER(KIND = 4) :: ixc_var_2021(klev_var_1933), ixlow_var_2022(kfdia_var_1932, klev_var_1933), ixhigh_var_2023(kfdia_var_1932, klev_var_1933)
  INTEGER(KIND = 4) :: ich_var_2024, icl_var_2025, ixc0_var_2026, ixp_var_2027, jc_var_2028, jl_var_2029
  laytrop_min_var_2019 = MINVAL(laytrop_var_1948)
  laytrop_max_var_2020 = MAXVAL(laytrop_var_1948)
  ixlow_var_2022 = 0
  ixhigh_var_2023 = 0
  ixc_var_2021 = 0
  DO lay_var_1969 = laytrop_min_var_2019 + 1, laytrop_max_var_2020
    icl_var_2025 = 0
    ich_var_2024 = 0
    DO jc_var_2028 = kidia_var_1931, kfdia_var_1932
      IF (lay_var_1969 <= laytrop_var_1948(jc_var_2028)) THEN
        icl_var_2025 = icl_var_2025 + 1
        ixlow_var_2022(icl_var_2025, lay_var_1969) = jc_var_2028
      ELSE
        ich_var_2024 = ich_var_2024 + 1
        ixhigh_var_2023(ich_var_2024, lay_var_1969) = jc_var_2028
      END IF
    END DO
    ixc_var_2021(lay_var_1969) = icl_var_2025
  END DO
  refrat_planck_a_var_1977 = chi_mls(1, 9) / chi_mls(6, 9)
  refrat_m_a_var_1978 = chi_mls(1, 3) / chi_mls(6, 3)
  DO lay_var_1969 = 1, laytrop_min_var_2019
    DO jl_var_2029 = kidia_var_1931, kfdia_var_1932
      speccomb_var_1973 = colh2o_var_1944(jl_var_2029, lay_var_1969) + rat_h2och4_var_1955(jl_var_2029, lay_var_1969) * colch4_var_1946(jl_var_2029, lay_var_1969)
      specparm_var_1981 = MIN(colh2o_var_1944(jl_var_2029, lay_var_1969) / speccomb_var_1973, 0.999999D0)
      specmult_var_1980 = 8.0D0 * (specparm_var_1981)
      js_var_1968 = 1 + INT(specmult_var_1980)
      fs_var_1979 = ((specmult_var_1980) - AINT((specmult_var_1980)))
      speccomb1_var_1974 = colh2o_var_1944(jl_var_2029, lay_var_1969) + rat_h2och4_1_var_1956(jl_var_2029, lay_var_1969) * colch4_var_1946(jl_var_2029, lay_var_1969)
      specparm1_var_1984 = MIN(colh2o_var_1944(jl_var_2029, lay_var_1969) / speccomb1_var_1974, 0.999999D0)
      specmult1_var_1983 = 8.0D0 * (specparm1_var_1984)
      js1_var_1970 = 1 + INT(specmult1_var_1983)
      fs1_var_1982 = ((specmult1_var_1983) - AINT((specmult1_var_1983)))
      speccomb_mn2o_var_1975 = colh2o_var_1944(jl_var_2029, lay_var_1969) + refrat_m_a_var_1978 * colch4_var_1946(jl_var_2029, lay_var_1969)
      specparm_mn2o_var_1987 = MIN(colh2o_var_1944(jl_var_2029, lay_var_1969) / speccomb_mn2o_var_1975, 0.999999D0)
      specmult_mn2o_var_1986 = 8.0D0 * specparm_mn2o_var_1987
      jmn2o_var_1971 = 1 + INT(specmult_mn2o_var_1986)
      fmn2o_var_1985 = ((specmult_mn2o_var_1986) - AINT((specmult_mn2o_var_1986)))
      chi_n2o_var_1994 = coln2o_var_1945(jl_var_2029, lay_var_1969) / (coldry_var_1947(jl_var_2029, lay_var_1969))
      ratn2o_var_1993 = 1D+20 * chi_n2o_var_1994 / chi_mls(4, jp_var_1940(jl_var_2029, lay_var_1969) + 1)
      IF (ratn2o_var_1993 .GT. 1.5D0) THEN
        adjfac_var_1991 = 0.5D0 + (ratn2o_var_1993 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1992 = adjfac_var_1991 * chi_mls(4, jp_var_1940(jl_var_2029, lay_var_1969) + 1) * coldry_var_1947(jl_var_2029, lay_var_1969) * 1D-20
      ELSE
        adjcoln2o_var_1992 = coln2o_var_1945(jl_var_2029, lay_var_1969)
      END IF
      speccomb_planck_var_1976 = colh2o_var_1944(jl_var_2029, lay_var_1969) + refrat_planck_a_var_1977 * colch4_var_1946(jl_var_2029, lay_var_1969)
      specparm_planck_var_1990 = MIN(colh2o_var_1944(jl_var_2029, lay_var_1969) / speccomb_planck_var_1976, 0.999999D0)
      specmult_planck_var_1989 = 8.0D0 * specparm_planck_var_1990
      jpl_var_1972 = 1 + INT(specmult_planck_var_1989)
      fpl_var_1988 = ((specmult_planck_var_1989) - AINT((specmult_planck_var_1989)))
      ind0_var_1962 = ((jp_var_1940(jl_var_2029, lay_var_1969) - 1) * 5 + (jt_var_1941(jl_var_2029, lay_var_1969) - 1)) * nspa_var_216(9) + js_var_1968
      ind1_var_1963 = (jp_var_1940(jl_var_2029, lay_var_1969) * 5 + (jt1_var_1942(jl_var_2029, lay_var_1969) - 1)) * nspa_var_216(9) + js1_var_1970
      inds_var_1964 = indself_var_1953(jl_var_2029, lay_var_1969)
      indf_var_1965 = indfor_var_1957(jl_var_2029, lay_var_1969)
      indm_var_1966 = indminor_var_1961(jl_var_2029, lay_var_1969)
      IF (specparm_var_1981 .LT. 0.125D0) THEN
        p_var_2007 = fs_var_1979 - 1.0D0
        p4_var_2008 = p_var_2007 ** 4
        fk0_var_2009 = p4_var_2008
        fk1_var_2010 = 1.0D0 - p_var_2007 - 2.0D0 * p4_var_2008
        fk2_var_2011 = p_var_2007 + p4_var_2008
        fac000_var_1995 = fk0_var_2009 * fac00_var_1936(jl_var_2029, lay_var_1969)
        fac100_var_1996 = fk1_var_2010 * fac00_var_1936(jl_var_2029, lay_var_1969)
        fac200_var_1997 = fk2_var_2011 * fac00_var_1936(jl_var_2029, lay_var_1969)
        fac010_var_1998 = fk0_var_2009 * fac10_var_1938(jl_var_2029, lay_var_1969)
        fac110_var_1999 = fk1_var_2010 * fac10_var_1938(jl_var_2029, lay_var_1969)
        fac210_var_2000 = fk2_var_2011 * fac10_var_1938(jl_var_2029, lay_var_1969)
      ELSE IF (specparm_var_1981 .GT. 0.875D0) THEN
        p_var_2007 = - fs_var_1979
        p4_var_2008 = p_var_2007 ** 4
        fk0_var_2009 = p4_var_2008
        fk1_var_2010 = 1.0D0 - p_var_2007 - 2.0D0 * p4_var_2008
        fk2_var_2011 = p_var_2007 + p4_var_2008
        fac000_var_1995 = fk0_var_2009 * fac00_var_1936(jl_var_2029, lay_var_1969)
        fac100_var_1996 = fk1_var_2010 * fac00_var_1936(jl_var_2029, lay_var_1969)
        fac200_var_1997 = fk2_var_2011 * fac00_var_1936(jl_var_2029, lay_var_1969)
        fac010_var_1998 = fk0_var_2009 * fac10_var_1938(jl_var_2029, lay_var_1969)
        fac110_var_1999 = fk1_var_2010 * fac10_var_1938(jl_var_2029, lay_var_1969)
        fac210_var_2000 = fk2_var_2011 * fac10_var_1938(jl_var_2029, lay_var_1969)
      ELSE
        fac000_var_1995 = (1.0D0 - fs_var_1979) * fac00_var_1936(jl_var_2029, lay_var_1969)
        fac010_var_1998 = (1.0D0 - fs_var_1979) * fac10_var_1938(jl_var_2029, lay_var_1969)
        fac100_var_1996 = fs_var_1979 * fac00_var_1936(jl_var_2029, lay_var_1969)
        fac110_var_1999 = fs_var_1979 * fac10_var_1938(jl_var_2029, lay_var_1969)
        fac200_var_1997 = 0.0D0
        fac210_var_2000 = 0.0D0
      END IF
      IF (specparm1_var_1984 .LT. 0.125D0) THEN
        p_var_2007 = fs1_var_1982 - 1.0D0
        p4_var_2008 = p_var_2007 ** 4
        fk0_var_2009 = p4_var_2008
        fk1_var_2010 = 1.0D0 - p_var_2007 - 2.0D0 * p4_var_2008
        fk2_var_2011 = p_var_2007 + p4_var_2008
        fac001_var_2001 = fk0_var_2009 * fac01_var_1937(jl_var_2029, lay_var_1969)
        fac101_var_2002 = fk1_var_2010 * fac01_var_1937(jl_var_2029, lay_var_1969)
        fac201_var_2003 = fk2_var_2011 * fac01_var_1937(jl_var_2029, lay_var_1969)
        fac011_var_2004 = fk0_var_2009 * fac11_var_1939(jl_var_2029, lay_var_1969)
        fac111_var_2005 = fk1_var_2010 * fac11_var_1939(jl_var_2029, lay_var_1969)
        fac211_var_2006 = fk2_var_2011 * fac11_var_1939(jl_var_2029, lay_var_1969)
      ELSE IF (specparm1_var_1984 .GT. 0.875D0) THEN
        p_var_2007 = - fs1_var_1982
        p4_var_2008 = p_var_2007 ** 4
        fk0_var_2009 = p4_var_2008
        fk1_var_2010 = 1.0D0 - p_var_2007 - 2.0D0 * p4_var_2008
        fk2_var_2011 = p_var_2007 + p4_var_2008
        fac001_var_2001 = fk0_var_2009 * fac01_var_1937(jl_var_2029, lay_var_1969)
        fac101_var_2002 = fk1_var_2010 * fac01_var_1937(jl_var_2029, lay_var_1969)
        fac201_var_2003 = fk2_var_2011 * fac01_var_1937(jl_var_2029, lay_var_1969)
        fac011_var_2004 = fk0_var_2009 * fac11_var_1939(jl_var_2029, lay_var_1969)
        fac111_var_2005 = fk1_var_2010 * fac11_var_1939(jl_var_2029, lay_var_1969)
        fac211_var_2006 = fk2_var_2011 * fac11_var_1939(jl_var_2029, lay_var_1969)
      ELSE
        fac001_var_2001 = (1.0D0 - fs1_var_1982) * fac01_var_1937(jl_var_2029, lay_var_1969)
        fac011_var_2004 = (1.0D0 - fs1_var_1982) * fac11_var_1939(jl_var_2029, lay_var_1969)
        fac101_var_2002 = fs1_var_1982 * fac01_var_1937(jl_var_2029, lay_var_1969)
        fac111_var_2005 = fs1_var_1982 * fac11_var_1939(jl_var_2029, lay_var_1969)
        fac201_var_2003 = 0.0D0
        fac211_var_2006 = 0.0D0
      END IF
      IF (specparm_var_1981 .LT. 0.125D0) THEN
        tau_major_var_2017(1 : ng9) = speccomb_var_1973 * (fac000_var_1995 * absa_var_208(ind0_var_1962, 1 : 12) + fac100_var_1996 * absa_var_208(ind0_var_1962 + 1, 1 : 12) + fac200_var_1997 * absa_var_208(ind0_var_1962 + 2, 1 : 12) + fac010_var_1998 * absa_var_208(ind0_var_1962 + 9, 1 : 12) + fac110_var_1999 * absa_var_208(ind0_var_1962 + 10, 1 : 12) + fac210_var_2000 * absa_var_208(ind0_var_1962 + 11, 1 : 12))
      ELSE IF (specparm_var_1981 .GT. 0.875D0) THEN
        tau_major_var_2017(1 : ng9) = speccomb_var_1973 * (fac200_var_1997 * absa_var_208(ind0_var_1962 - 1, 1 : 12) + fac100_var_1996 * absa_var_208(ind0_var_1962, 1 : 12) + fac000_var_1995 * absa_var_208(ind0_var_1962 + 1, 1 : 12) + fac210_var_2000 * absa_var_208(ind0_var_1962 + 8, 1 : 12) + fac110_var_1999 * absa_var_208(ind0_var_1962 + 9, 1 : 12) + fac010_var_1998 * absa_var_208(ind0_var_1962 + 10, 1 : 12))
      ELSE
        tau_major_var_2017(1 : ng9) = speccomb_var_1973 * (fac000_var_1995 * absa_var_208(ind0_var_1962, 1 : 12) + fac100_var_1996 * absa_var_208(ind0_var_1962 + 1, 1 : 12) + fac010_var_1998 * absa_var_208(ind0_var_1962 + 9, 1 : 12) + fac110_var_1999 * absa_var_208(ind0_var_1962 + 10, 1 : 12))
      END IF
      IF (specparm1_var_1984 .LT. 0.125D0) THEN
        tau_major1_var_2018(1 : ng9) = speccomb1_var_1974 * (fac001_var_2001 * absa_var_208(ind1_var_1963, 1 : 12) + fac101_var_2002 * absa_var_208(ind1_var_1963 + 1, 1 : 12) + fac201_var_2003 * absa_var_208(ind1_var_1963 + 2, 1 : 12) + fac011_var_2004 * absa_var_208(ind1_var_1963 + 9, 1 : 12) + fac111_var_2005 * absa_var_208(ind1_var_1963 + 10, 1 : 12) + fac211_var_2006 * absa_var_208(ind1_var_1963 + 11, 1 : 12))
      ELSE IF (specparm1_var_1984 .GT. 0.875D0) THEN
        tau_major1_var_2018(1 : ng9) = speccomb1_var_1974 * (fac201_var_2003 * absa_var_208(ind1_var_1963 - 1, 1 : 12) + fac101_var_2002 * absa_var_208(ind1_var_1963, 1 : 12) + fac001_var_2001 * absa_var_208(ind1_var_1963 + 1, 1 : 12) + fac211_var_2006 * absa_var_208(ind1_var_1963 + 8, 1 : 12) + fac111_var_2005 * absa_var_208(ind1_var_1963 + 9, 1 : 12) + fac011_var_2004 * absa_var_208(ind1_var_1963 + 10, 1 : 12))
      ELSE
        tau_major1_var_2018(1 : ng9) = speccomb1_var_1974 * (fac001_var_2001 * absa_var_208(ind1_var_1963, 1 : 12) + fac101_var_2002 * absa_var_208(ind1_var_1963 + 1, 1 : 12) + fac011_var_2004 * absa_var_208(ind1_var_1963 + 9, 1 : 12) + fac111_var_2005 * absa_var_208(ind1_var_1963 + 10, 1 : 12))
      END IF
      DO ig_var_1967 = 1, 12
        tauself_var_2013 = selffac_var_1951(jl_var_2029, lay_var_1969) * (selfref_var_212(inds_var_1964, ig_var_1967) + selffrac_var_1952(jl_var_2029, lay_var_1969) * (selfref_var_212(inds_var_1964 + 1, ig_var_1967) - selfref_var_212(inds_var_1964, ig_var_1967)))
        taufor_var_2012 = forfac_var_1958(jl_var_2029, lay_var_1969) * (forref_var_213(indf_var_1965, ig_var_1967) + forfrac_var_1959(jl_var_2029, lay_var_1969) * (forref_var_213(indf_var_1965 + 1, ig_var_1967) - forref_var_213(indf_var_1965, ig_var_1967)))
        n2om1_var_2014 = ka_mn2o_var_210(jmn2o_var_1971, indm_var_1966, ig_var_1967) + fmn2o_var_1985 * (ka_mn2o_var_210(jmn2o_var_1971 + 1, indm_var_1966, ig_var_1967) - ka_mn2o_var_210(jmn2o_var_1971, indm_var_1966, ig_var_1967))
        n2om2_var_2015 = ka_mn2o_var_210(jmn2o_var_1971, indm_var_1966 + 1, ig_var_1967) + fmn2o_var_1985 * (ka_mn2o_var_210(jmn2o_var_1971 + 1, indm_var_1966 + 1, ig_var_1967) - ka_mn2o_var_210(jmn2o_var_1971, indm_var_1966 + 1, ig_var_1967))
        absn2o_var_2016 = n2om1_var_2014 + minorfrac_var_1960(jl_var_2029, lay_var_1969) * (n2om2_var_2015 - n2om1_var_2014)
        taug_var_1934(jl_var_2029, 96 + ig_var_1967, lay_var_1969) = tau_major_var_2017(ig_var_1967) + tau_major1_var_2018(ig_var_1967) + tauself_var_2013 + taufor_var_2012 + adjcoln2o_var_1992 * absn2o_var_2016
        fracs_var_1954(jl_var_2029, 96 + ig_var_1967, lay_var_1969) = fracrefa_var_206(ig_var_1967, jpl_var_1972) + fpl_var_1988 * (fracrefa_var_206(ig_var_1967, jpl_var_1972 + 1) - fracrefa_var_206(ig_var_1967, jpl_var_1972))
      END DO
    END DO
  END DO
  DO lay_var_1969 = laytrop_max_var_2020 + 1, klev_var_1933
    DO jl_var_2029 = kidia_var_1931, kfdia_var_1932
      chi_n2o_var_1994 = coln2o_var_1945(jl_var_2029, lay_var_1969) / (coldry_var_1947(jl_var_2029, lay_var_1969))
      ratn2o_var_1993 = 1D+20 * chi_n2o_var_1994 / chi_mls(4, jp_var_1940(jl_var_2029, lay_var_1969) + 1)
      IF (ratn2o_var_1993 .GT. 1.5D0) THEN
        adjfac_var_1991 = 0.5D0 + (ratn2o_var_1993 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1992 = adjfac_var_1991 * chi_mls(4, jp_var_1940(jl_var_2029, lay_var_1969) + 1) * coldry_var_1947(jl_var_2029, lay_var_1969) * 1D-20
      ELSE
        adjcoln2o_var_1992 = coln2o_var_1945(jl_var_2029, lay_var_1969)
      END IF
      ind0_var_1962 = ((jp_var_1940(jl_var_2029, lay_var_1969) - 13) * 5 + (jt_var_1941(jl_var_2029, lay_var_1969) - 1)) * nspb_var_217(9) + 1
      ind1_var_1963 = ((jp_var_1940(jl_var_2029, lay_var_1969) - 12) * 5 + (jt1_var_1942(jl_var_2029, lay_var_1969) - 1)) * nspb_var_217(9) + 1
      indm_var_1966 = indminor_var_1961(jl_var_2029, lay_var_1969)
      DO ig_var_1967 = 1, 12
        absn2o_var_2016 = kb_mn2o_var_211(indm_var_1966, ig_var_1967) + minorfrac_var_1960(jl_var_2029, lay_var_1969) * (kb_mn2o_var_211(indm_var_1966 + 1, ig_var_1967) - kb_mn2o_var_211(indm_var_1966, ig_var_1967))
        taug_var_1934(jl_var_2029, 96 + ig_var_1967, lay_var_1969) = colch4_var_1946(jl_var_2029, lay_var_1969) * (fac00_var_1936(jl_var_2029, lay_var_1969) * absb_var_209(ind0_var_1962, ig_var_1967) + fac10_var_1938(jl_var_2029, lay_var_1969) * absb_var_209(ind0_var_1962 + 1, ig_var_1967) + fac01_var_1937(jl_var_2029, lay_var_1969) * absb_var_209(ind1_var_1963, ig_var_1967) + fac11_var_1939(jl_var_2029, lay_var_1969) * absb_var_209(ind1_var_1963 + 1, ig_var_1967)) + adjcoln2o_var_1992 * absn2o_var_2016
        fracs_var_1954(jl_var_2029, 96 + ig_var_1967, lay_var_1969) = fracrefb_var_207(ig_var_1967)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2020 /= laytrop_min_var_2019) THEN
    DO lay_var_1969 = laytrop_min_var_2019 + 1, laytrop_max_var_2020
      ixc0_var_2026 = ixc_var_2021(lay_var_1969)
      DO ixp_var_2027 = 1, ixc0_var_2026
        jl_var_2029 = ixlow_var_2022(ixp_var_2027, lay_var_1969)
        speccomb_var_1973 = colh2o_var_1944(jl_var_2029, lay_var_1969) + rat_h2och4_var_1955(jl_var_2029, lay_var_1969) * colch4_var_1946(jl_var_2029, lay_var_1969)
        specparm_var_1981 = MIN(colh2o_var_1944(jl_var_2029, lay_var_1969) / speccomb_var_1973, 0.999999D0)
        specmult_var_1980 = 8.0D0 * (specparm_var_1981)
        js_var_1968 = 1 + INT(specmult_var_1980)
        fs_var_1979 = ((specmult_var_1980) - AINT((specmult_var_1980)))
        speccomb1_var_1974 = colh2o_var_1944(jl_var_2029, lay_var_1969) + rat_h2och4_1_var_1956(jl_var_2029, lay_var_1969) * colch4_var_1946(jl_var_2029, lay_var_1969)
        specparm1_var_1984 = MIN(colh2o_var_1944(jl_var_2029, lay_var_1969) / speccomb1_var_1974, 0.999999D0)
        specmult1_var_1983 = 8.0D0 * (specparm1_var_1984)
        js1_var_1970 = 1 + INT(specmult1_var_1983)
        fs1_var_1982 = ((specmult1_var_1983) - AINT((specmult1_var_1983)))
        speccomb_mn2o_var_1975 = colh2o_var_1944(jl_var_2029, lay_var_1969) + refrat_m_a_var_1978 * colch4_var_1946(jl_var_2029, lay_var_1969)
        specparm_mn2o_var_1987 = MIN(colh2o_var_1944(jl_var_2029, lay_var_1969) / speccomb_mn2o_var_1975, 0.999999D0)
        specmult_mn2o_var_1986 = 8.0D0 * specparm_mn2o_var_1987
        jmn2o_var_1971 = 1 + INT(specmult_mn2o_var_1986)
        fmn2o_var_1985 = ((specmult_mn2o_var_1986) - AINT((specmult_mn2o_var_1986)))
        chi_n2o_var_1994 = coln2o_var_1945(jl_var_2029, lay_var_1969) / (coldry_var_1947(jl_var_2029, lay_var_1969))
        ratn2o_var_1993 = 1D+20 * chi_n2o_var_1994 / chi_mls(4, jp_var_1940(jl_var_2029, lay_var_1969) + 1)
        IF (ratn2o_var_1993 .GT. 1.5D0) THEN
          adjfac_var_1991 = 0.5D0 + (ratn2o_var_1993 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1992 = adjfac_var_1991 * chi_mls(4, jp_var_1940(jl_var_2029, lay_var_1969) + 1) * coldry_var_1947(jl_var_2029, lay_var_1969) * 1D-20
        ELSE
          adjcoln2o_var_1992 = coln2o_var_1945(jl_var_2029, lay_var_1969)
        END IF
        speccomb_planck_var_1976 = colh2o_var_1944(jl_var_2029, lay_var_1969) + refrat_planck_a_var_1977 * colch4_var_1946(jl_var_2029, lay_var_1969)
        specparm_planck_var_1990 = MIN(colh2o_var_1944(jl_var_2029, lay_var_1969) / speccomb_planck_var_1976, 0.999999D0)
        specmult_planck_var_1989 = 8.0D0 * specparm_planck_var_1990
        jpl_var_1972 = 1 + INT(specmult_planck_var_1989)
        fpl_var_1988 = ((specmult_planck_var_1989) - AINT((specmult_planck_var_1989)))
        ind0_var_1962 = ((jp_var_1940(jl_var_2029, lay_var_1969) - 1) * 5 + (jt_var_1941(jl_var_2029, lay_var_1969) - 1)) * nspa_var_216(9) + js_var_1968
        ind1_var_1963 = (jp_var_1940(jl_var_2029, lay_var_1969) * 5 + (jt1_var_1942(jl_var_2029, lay_var_1969) - 1)) * nspa_var_216(9) + js1_var_1970
        inds_var_1964 = indself_var_1953(jl_var_2029, lay_var_1969)
        indf_var_1965 = indfor_var_1957(jl_var_2029, lay_var_1969)
        indm_var_1966 = indminor_var_1961(jl_var_2029, lay_var_1969)
        IF (specparm_var_1981 .LT. 0.125D0) THEN
          p_var_2007 = fs_var_1979 - 1.0D0
          p4_var_2008 = p_var_2007 ** 4
          fk0_var_2009 = p4_var_2008
          fk1_var_2010 = 1.0D0 - p_var_2007 - 2.0D0 * p4_var_2008
          fk2_var_2011 = p_var_2007 + p4_var_2008
          fac000_var_1995 = fk0_var_2009 * fac00_var_1936(jl_var_2029, lay_var_1969)
          fac100_var_1996 = fk1_var_2010 * fac00_var_1936(jl_var_2029, lay_var_1969)
          fac200_var_1997 = fk2_var_2011 * fac00_var_1936(jl_var_2029, lay_var_1969)
          fac010_var_1998 = fk0_var_2009 * fac10_var_1938(jl_var_2029, lay_var_1969)
          fac110_var_1999 = fk1_var_2010 * fac10_var_1938(jl_var_2029, lay_var_1969)
          fac210_var_2000 = fk2_var_2011 * fac10_var_1938(jl_var_2029, lay_var_1969)
        ELSE IF (specparm_var_1981 .GT. 0.875D0) THEN
          p_var_2007 = - fs_var_1979
          p4_var_2008 = p_var_2007 ** 4
          fk0_var_2009 = p4_var_2008
          fk1_var_2010 = 1.0D0 - p_var_2007 - 2.0D0 * p4_var_2008
          fk2_var_2011 = p_var_2007 + p4_var_2008
          fac000_var_1995 = fk0_var_2009 * fac00_var_1936(jl_var_2029, lay_var_1969)
          fac100_var_1996 = fk1_var_2010 * fac00_var_1936(jl_var_2029, lay_var_1969)
          fac200_var_1997 = fk2_var_2011 * fac00_var_1936(jl_var_2029, lay_var_1969)
          fac010_var_1998 = fk0_var_2009 * fac10_var_1938(jl_var_2029, lay_var_1969)
          fac110_var_1999 = fk1_var_2010 * fac10_var_1938(jl_var_2029, lay_var_1969)
          fac210_var_2000 = fk2_var_2011 * fac10_var_1938(jl_var_2029, lay_var_1969)
        ELSE
          fac000_var_1995 = (1.0D0 - fs_var_1979) * fac00_var_1936(jl_var_2029, lay_var_1969)
          fac010_var_1998 = (1.0D0 - fs_var_1979) * fac10_var_1938(jl_var_2029, lay_var_1969)
          fac100_var_1996 = fs_var_1979 * fac00_var_1936(jl_var_2029, lay_var_1969)
          fac110_var_1999 = fs_var_1979 * fac10_var_1938(jl_var_2029, lay_var_1969)
          fac200_var_1997 = 0.0D0
          fac210_var_2000 = 0.0D0
        END IF
        IF (specparm1_var_1984 .LT. 0.125D0) THEN
          p_var_2007 = fs1_var_1982 - 1.0D0
          p4_var_2008 = p_var_2007 ** 4
          fk0_var_2009 = p4_var_2008
          fk1_var_2010 = 1.0D0 - p_var_2007 - 2.0D0 * p4_var_2008
          fk2_var_2011 = p_var_2007 + p4_var_2008
          fac001_var_2001 = fk0_var_2009 * fac01_var_1937(jl_var_2029, lay_var_1969)
          fac101_var_2002 = fk1_var_2010 * fac01_var_1937(jl_var_2029, lay_var_1969)
          fac201_var_2003 = fk2_var_2011 * fac01_var_1937(jl_var_2029, lay_var_1969)
          fac011_var_2004 = fk0_var_2009 * fac11_var_1939(jl_var_2029, lay_var_1969)
          fac111_var_2005 = fk1_var_2010 * fac11_var_1939(jl_var_2029, lay_var_1969)
          fac211_var_2006 = fk2_var_2011 * fac11_var_1939(jl_var_2029, lay_var_1969)
        ELSE IF (specparm1_var_1984 .GT. 0.875D0) THEN
          p_var_2007 = - fs1_var_1982
          p4_var_2008 = p_var_2007 ** 4
          fk0_var_2009 = p4_var_2008
          fk1_var_2010 = 1.0D0 - p_var_2007 - 2.0D0 * p4_var_2008
          fk2_var_2011 = p_var_2007 + p4_var_2008
          fac001_var_2001 = fk0_var_2009 * fac01_var_1937(jl_var_2029, lay_var_1969)
          fac101_var_2002 = fk1_var_2010 * fac01_var_1937(jl_var_2029, lay_var_1969)
          fac201_var_2003 = fk2_var_2011 * fac01_var_1937(jl_var_2029, lay_var_1969)
          fac011_var_2004 = fk0_var_2009 * fac11_var_1939(jl_var_2029, lay_var_1969)
          fac111_var_2005 = fk1_var_2010 * fac11_var_1939(jl_var_2029, lay_var_1969)
          fac211_var_2006 = fk2_var_2011 * fac11_var_1939(jl_var_2029, lay_var_1969)
        ELSE
          fac001_var_2001 = (1.0D0 - fs1_var_1982) * fac01_var_1937(jl_var_2029, lay_var_1969)
          fac011_var_2004 = (1.0D0 - fs1_var_1982) * fac11_var_1939(jl_var_2029, lay_var_1969)
          fac101_var_2002 = fs1_var_1982 * fac01_var_1937(jl_var_2029, lay_var_1969)
          fac111_var_2005 = fs1_var_1982 * fac11_var_1939(jl_var_2029, lay_var_1969)
          fac201_var_2003 = 0.0D0
          fac211_var_2006 = 0.0D0
        END IF
        IF (specparm_var_1981 .LT. 0.125D0) THEN
          tau_major_var_2017(1 : ng9) = speccomb_var_1973 * (fac000_var_1995 * absa_var_208(ind0_var_1962, 1 : 12) + fac100_var_1996 * absa_var_208(ind0_var_1962 + 1, 1 : 12) + fac200_var_1997 * absa_var_208(ind0_var_1962 + 2, 1 : 12) + fac010_var_1998 * absa_var_208(ind0_var_1962 + 9, 1 : 12) + fac110_var_1999 * absa_var_208(ind0_var_1962 + 10, 1 : 12) + fac210_var_2000 * absa_var_208(ind0_var_1962 + 11, 1 : 12))
        ELSE IF (specparm_var_1981 .GT. 0.875D0) THEN
          tau_major_var_2017(1 : ng9) = speccomb_var_1973 * (fac200_var_1997 * absa_var_208(ind0_var_1962 - 1, 1 : 12) + fac100_var_1996 * absa_var_208(ind0_var_1962, 1 : 12) + fac000_var_1995 * absa_var_208(ind0_var_1962 + 1, 1 : 12) + fac210_var_2000 * absa_var_208(ind0_var_1962 + 8, 1 : 12) + fac110_var_1999 * absa_var_208(ind0_var_1962 + 9, 1 : 12) + fac010_var_1998 * absa_var_208(ind0_var_1962 + 10, 1 : 12))
        ELSE
          tau_major_var_2017(1 : ng9) = speccomb_var_1973 * (fac000_var_1995 * absa_var_208(ind0_var_1962, 1 : 12) + fac100_var_1996 * absa_var_208(ind0_var_1962 + 1, 1 : 12) + fac010_var_1998 * absa_var_208(ind0_var_1962 + 9, 1 : 12) + fac110_var_1999 * absa_var_208(ind0_var_1962 + 10, 1 : 12))
        END IF
        IF (specparm1_var_1984 .LT. 0.125D0) THEN
          tau_major1_var_2018(1 : ng9) = speccomb1_var_1974 * (fac001_var_2001 * absa_var_208(ind1_var_1963, 1 : 12) + fac101_var_2002 * absa_var_208(ind1_var_1963 + 1, 1 : 12) + fac201_var_2003 * absa_var_208(ind1_var_1963 + 2, 1 : 12) + fac011_var_2004 * absa_var_208(ind1_var_1963 + 9, 1 : 12) + fac111_var_2005 * absa_var_208(ind1_var_1963 + 10, 1 : 12) + fac211_var_2006 * absa_var_208(ind1_var_1963 + 11, 1 : 12))
        ELSE IF (specparm1_var_1984 .GT. 0.875D0) THEN
          tau_major1_var_2018(1 : ng9) = speccomb1_var_1974 * (fac201_var_2003 * absa_var_208(ind1_var_1963 - 1, 1 : 12) + fac101_var_2002 * absa_var_208(ind1_var_1963, 1 : 12) + fac001_var_2001 * absa_var_208(ind1_var_1963 + 1, 1 : 12) + fac211_var_2006 * absa_var_208(ind1_var_1963 + 8, 1 : 12) + fac111_var_2005 * absa_var_208(ind1_var_1963 + 9, 1 : 12) + fac011_var_2004 * absa_var_208(ind1_var_1963 + 10, 1 : 12))
        ELSE
          tau_major1_var_2018(1 : ng9) = speccomb1_var_1974 * (fac001_var_2001 * absa_var_208(ind1_var_1963, 1 : 12) + fac101_var_2002 * absa_var_208(ind1_var_1963 + 1, 1 : 12) + fac011_var_2004 * absa_var_208(ind1_var_1963 + 9, 1 : 12) + fac111_var_2005 * absa_var_208(ind1_var_1963 + 10, 1 : 12))
        END IF
        DO ig_var_1967 = 1, 12
          tauself_var_2013 = selffac_var_1951(jl_var_2029, lay_var_1969) * (selfref_var_212(inds_var_1964, ig_var_1967) + selffrac_var_1952(jl_var_2029, lay_var_1969) * (selfref_var_212(inds_var_1964 + 1, ig_var_1967) - selfref_var_212(inds_var_1964, ig_var_1967)))
          taufor_var_2012 = forfac_var_1958(jl_var_2029, lay_var_1969) * (forref_var_213(indf_var_1965, ig_var_1967) + forfrac_var_1959(jl_var_2029, lay_var_1969) * (forref_var_213(indf_var_1965 + 1, ig_var_1967) - forref_var_213(indf_var_1965, ig_var_1967)))
          n2om1_var_2014 = ka_mn2o_var_210(jmn2o_var_1971, indm_var_1966, ig_var_1967) + fmn2o_var_1985 * (ka_mn2o_var_210(jmn2o_var_1971 + 1, indm_var_1966, ig_var_1967) - ka_mn2o_var_210(jmn2o_var_1971, indm_var_1966, ig_var_1967))
          n2om2_var_2015 = ka_mn2o_var_210(jmn2o_var_1971, indm_var_1966 + 1, ig_var_1967) + fmn2o_var_1985 * (ka_mn2o_var_210(jmn2o_var_1971 + 1, indm_var_1966 + 1, ig_var_1967) - ka_mn2o_var_210(jmn2o_var_1971, indm_var_1966 + 1, ig_var_1967))
          absn2o_var_2016 = n2om1_var_2014 + minorfrac_var_1960(jl_var_2029, lay_var_1969) * (n2om2_var_2015 - n2om1_var_2014)
          taug_var_1934(jl_var_2029, 96 + ig_var_1967, lay_var_1969) = tau_major_var_2017(ig_var_1967) + tau_major1_var_2018(ig_var_1967) + tauself_var_2013 + taufor_var_2012 + adjcoln2o_var_1992 * absn2o_var_2016
          fracs_var_1954(jl_var_2029, 96 + ig_var_1967, lay_var_1969) = fracrefa_var_206(ig_var_1967, jpl_var_1972) + fpl_var_1988 * (fracrefa_var_206(ig_var_1967, jpl_var_1972 + 1) - fracrefa_var_206(ig_var_1967, jpl_var_1972))
        END DO
      END DO
      ixc0_var_2026 = kfdia_var_1932 - kidia_var_1931 + 1 - ixc0_var_2026
      DO ixp_var_2027 = 1, ixc0_var_2026
        jl_var_2029 = ixhigh_var_2023(ixp_var_2027, lay_var_1969)
        chi_n2o_var_1994 = coln2o_var_1945(jl_var_2029, lay_var_1969) / (coldry_var_1947(jl_var_2029, lay_var_1969))
        ratn2o_var_1993 = 1D+20 * chi_n2o_var_1994 / chi_mls(4, jp_var_1940(jl_var_2029, lay_var_1969) + 1)
        IF (ratn2o_var_1993 .GT. 1.5D0) THEN
          adjfac_var_1991 = 0.5D0 + (ratn2o_var_1993 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1992 = adjfac_var_1991 * chi_mls(4, jp_var_1940(jl_var_2029, lay_var_1969) + 1) * coldry_var_1947(jl_var_2029, lay_var_1969) * 1D-20
        ELSE
          adjcoln2o_var_1992 = coln2o_var_1945(jl_var_2029, lay_var_1969)
        END IF
        ind0_var_1962 = ((jp_var_1940(jl_var_2029, lay_var_1969) - 13) * 5 + (jt_var_1941(jl_var_2029, lay_var_1969) - 1)) * nspb_var_217(9) + 1
        ind1_var_1963 = ((jp_var_1940(jl_var_2029, lay_var_1969) - 12) * 5 + (jt1_var_1942(jl_var_2029, lay_var_1969) - 1)) * nspb_var_217(9) + 1
        indm_var_1966 = indminor_var_1961(jl_var_2029, lay_var_1969)
        DO ig_var_1967 = 1, 12
          absn2o_var_2016 = kb_mn2o_var_211(indm_var_1966, ig_var_1967) + minorfrac_var_1960(jl_var_2029, lay_var_1969) * (kb_mn2o_var_211(indm_var_1966 + 1, ig_var_1967) - kb_mn2o_var_211(indm_var_1966, ig_var_1967))
          taug_var_1934(jl_var_2029, 96 + ig_var_1967, lay_var_1969) = colch4_var_1946(jl_var_2029, lay_var_1969) * (fac00_var_1936(jl_var_2029, lay_var_1969) * absb_var_209(ind0_var_1962, ig_var_1967) + fac10_var_1938(jl_var_2029, lay_var_1969) * absb_var_209(ind0_var_1962 + 1, ig_var_1967) + fac01_var_1937(jl_var_2029, lay_var_1969) * absb_var_209(ind1_var_1963, ig_var_1967) + fac11_var_1939(jl_var_2029, lay_var_1969) * absb_var_209(ind1_var_1963 + 1, ig_var_1967)) + adjcoln2o_var_1992 * absn2o_var_2016
          fracs_var_1954(jl_var_2029, 96 + ig_var_1967, lay_var_1969) = fracrefb_var_207(ig_var_1967)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol9
SUBROUTINE rrtm_taumol11(kidia_var_2030, kfdia_var_2031, klev_var_2032, taug_var_2033, p_tauaerl_var_2034, fac00_var_2035, fac01_var_2036, fac10_var_2037, fac11_var_2038, forfac_var_2049, forfrac_var_2050, indfor_var_2048, jp_var_2039, jt_var_2040, jt1_var_2041, colh2o_var_2042, colo2, laytrop_var_2043, selffac_var_2044, selffrac_var_2045, indself_var_2046, fracs_var_2047, minorfrac_var_2051, indminor_var_2052, scaleminor_var_2053)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta11, ONLY: absa_var_122, absb_var_123, forref_var_125, fracrefa_var_120, fracrefb_var_121, ka_mo2, kb_mo2, selfref_var_124
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2030
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2031
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2032
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2033(kidia_var_2030 : kfdia_var_2031, 140, klev_var_2032)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2034(kidia_var_2030 : kfdia_var_2031, klev_var_2032, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2035(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2036(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2037(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2038(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2039(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2040(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2041(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2042(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  REAL(KIND = 8), INTENT(IN) :: colo2(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2043(kidia_var_2030 : kfdia_var_2031)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2044(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2045(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2046(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2047(kidia_var_2030 : kfdia_var_2031, 140, klev_var_2032)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2048(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2049(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2050(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2051(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2052(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  REAL(KIND = 8), INTENT(IN) :: scaleminor_var_2053(kidia_var_2030 : kfdia_var_2031, klev_var_2032)
  INTEGER(KIND = 4) :: ind0_var_2054, ind1_var_2055
  INTEGER(KIND = 4) :: inds_var_2056, indf_var_2057, indm_var_2058
  INTEGER(KIND = 4) :: ig_var_2059, lay_var_2060
  REAL(KIND = 8) :: taufor_var_2061, tauself_var_2062, scaleo2, tauo2
  INTEGER(KIND = 4) :: laytrop_min_var_2063, laytrop_max_var_2064
  INTEGER(KIND = 4) :: ixc_var_2065(klev_var_2032), ixlow_var_2066(kfdia_var_2031, klev_var_2032), ixhigh_var_2067(kfdia_var_2031, klev_var_2032)
  INTEGER(KIND = 4) :: ich_var_2068, icl_var_2069, ixc0_var_2070, ixp_var_2071, jc_var_2072, jl_var_2073
  laytrop_min_var_2063 = MINVAL(laytrop_var_2043)
  laytrop_max_var_2064 = MAXVAL(laytrop_var_2043)
  ixlow_var_2066 = 0
  ixhigh_var_2067 = 0
  ixc_var_2065 = 0
  DO lay_var_2060 = laytrop_min_var_2063 + 1, laytrop_max_var_2064
    icl_var_2069 = 0
    ich_var_2068 = 0
    DO jc_var_2072 = kidia_var_2030, kfdia_var_2031
      IF (lay_var_2060 <= laytrop_var_2043(jc_var_2072)) THEN
        icl_var_2069 = icl_var_2069 + 1
        ixlow_var_2066(icl_var_2069, lay_var_2060) = jc_var_2072
      ELSE
        ich_var_2068 = ich_var_2068 + 1
        ixhigh_var_2067(ich_var_2068, lay_var_2060) = jc_var_2072
      END IF
    END DO
    ixc_var_2065(lay_var_2060) = icl_var_2069
  END DO
  DO lay_var_2060 = 1, laytrop_min_var_2063
    DO jl_var_2073 = kidia_var_2030, kfdia_var_2031
      ind0_var_2054 = ((jp_var_2039(jl_var_2073, lay_var_2060) - 1) * 5 + (jt_var_2040(jl_var_2073, lay_var_2060) - 1)) * nspa_var_216(11) + 1
      ind1_var_2055 = (jp_var_2039(jl_var_2073, lay_var_2060) * 5 + (jt1_var_2041(jl_var_2073, lay_var_2060) - 1)) * nspa_var_216(11) + 1
      inds_var_2056 = indself_var_2046(jl_var_2073, lay_var_2060)
      indf_var_2057 = indfor_var_2048(jl_var_2073, lay_var_2060)
      indm_var_2058 = indminor_var_2052(jl_var_2073, lay_var_2060)
      scaleo2 = colo2(jl_var_2073, lay_var_2060) * scaleminor_var_2053(jl_var_2073, lay_var_2060)
      DO ig_var_2059 = 1, 8
        tauself_var_2062 = selffac_var_2044(jl_var_2073, lay_var_2060) * (selfref_var_124(inds_var_2056, ig_var_2059) + selffrac_var_2045(jl_var_2073, lay_var_2060) * (selfref_var_124(inds_var_2056 + 1, ig_var_2059) - selfref_var_124(inds_var_2056, ig_var_2059)))
        taufor_var_2061 = forfac_var_2049(jl_var_2073, lay_var_2060) * (forref_var_125(indf_var_2057, ig_var_2059) + forfrac_var_2050(jl_var_2073, lay_var_2060) * (forref_var_125(indf_var_2057 + 1, ig_var_2059) - forref_var_125(indf_var_2057, ig_var_2059)))
        tauo2 = scaleo2 * (ka_mo2(indm_var_2058, ig_var_2059) + minorfrac_var_2051(jl_var_2073, lay_var_2060) * (ka_mo2(indm_var_2058 + 1, ig_var_2059) - ka_mo2(indm_var_2058, ig_var_2059)))
        taug_var_2033(jl_var_2073, 114 + ig_var_2059, lay_var_2060) = colh2o_var_2042(jl_var_2073, lay_var_2060) * (fac00_var_2035(jl_var_2073, lay_var_2060) * absa_var_122(ind0_var_2054, ig_var_2059) + fac10_var_2037(jl_var_2073, lay_var_2060) * absa_var_122(ind0_var_2054 + 1, ig_var_2059) + fac01_var_2036(jl_var_2073, lay_var_2060) * absa_var_122(ind1_var_2055, ig_var_2059) + fac11_var_2038(jl_var_2073, lay_var_2060) * absa_var_122(ind1_var_2055 + 1, ig_var_2059)) + tauself_var_2062 + taufor_var_2061 + tauo2
        fracs_var_2047(jl_var_2073, 114 + ig_var_2059, lay_var_2060) = fracrefa_var_120(ig_var_2059)
      END DO
    END DO
  END DO
  DO lay_var_2060 = laytrop_max_var_2064 + 1, klev_var_2032
    DO jl_var_2073 = kidia_var_2030, kfdia_var_2031
      ind0_var_2054 = ((jp_var_2039(jl_var_2073, lay_var_2060) - 13) * 5 + (jt_var_2040(jl_var_2073, lay_var_2060) - 1)) * nspb_var_217(11) + 1
      ind1_var_2055 = ((jp_var_2039(jl_var_2073, lay_var_2060) - 12) * 5 + (jt1_var_2041(jl_var_2073, lay_var_2060) - 1)) * nspb_var_217(11) + 1
      indf_var_2057 = indfor_var_2048(jl_var_2073, lay_var_2060)
      indm_var_2058 = indminor_var_2052(jl_var_2073, lay_var_2060)
      scaleo2 = colo2(jl_var_2073, lay_var_2060) * scaleminor_var_2053(jl_var_2073, lay_var_2060)
      DO ig_var_2059 = 1, 8
        taufor_var_2061 = forfac_var_2049(jl_var_2073, lay_var_2060) * (forref_var_125(indf_var_2057, ig_var_2059) + forfrac_var_2050(jl_var_2073, lay_var_2060) * (forref_var_125(indf_var_2057 + 1, ig_var_2059) - forref_var_125(indf_var_2057, ig_var_2059)))
        tauo2 = scaleo2 * (kb_mo2(indm_var_2058, ig_var_2059) + minorfrac_var_2051(jl_var_2073, lay_var_2060) * (kb_mo2(indm_var_2058 + 1, ig_var_2059) - kb_mo2(indm_var_2058, ig_var_2059)))
        taug_var_2033(jl_var_2073, 114 + ig_var_2059, lay_var_2060) = colh2o_var_2042(jl_var_2073, lay_var_2060) * (fac00_var_2035(jl_var_2073, lay_var_2060) * absb_var_123(ind0_var_2054, ig_var_2059) + fac10_var_2037(jl_var_2073, lay_var_2060) * absb_var_123(ind0_var_2054 + 1, ig_var_2059) + fac01_var_2036(jl_var_2073, lay_var_2060) * absb_var_123(ind1_var_2055, ig_var_2059) + fac11_var_2038(jl_var_2073, lay_var_2060) * absb_var_123(ind1_var_2055 + 1, ig_var_2059)) + taufor_var_2061 + tauo2
        fracs_var_2047(jl_var_2073, 114 + ig_var_2059, lay_var_2060) = fracrefb_var_121(ig_var_2059)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2064 /= laytrop_min_var_2063) THEN
    DO lay_var_2060 = laytrop_min_var_2063 + 1, laytrop_max_var_2064
      ixc0_var_2070 = ixc_var_2065(lay_var_2060)
      DO ixp_var_2071 = 1, ixc0_var_2070
        jl_var_2073 = ixlow_var_2066(ixp_var_2071, lay_var_2060)
        ind0_var_2054 = ((jp_var_2039(jl_var_2073, lay_var_2060) - 1) * 5 + (jt_var_2040(jl_var_2073, lay_var_2060) - 1)) * nspa_var_216(11) + 1
        ind1_var_2055 = (jp_var_2039(jl_var_2073, lay_var_2060) * 5 + (jt1_var_2041(jl_var_2073, lay_var_2060) - 1)) * nspa_var_216(11) + 1
        inds_var_2056 = indself_var_2046(jl_var_2073, lay_var_2060)
        indf_var_2057 = indfor_var_2048(jl_var_2073, lay_var_2060)
        indm_var_2058 = indminor_var_2052(jl_var_2073, lay_var_2060)
        scaleo2 = colo2(jl_var_2073, lay_var_2060) * scaleminor_var_2053(jl_var_2073, lay_var_2060)
        DO ig_var_2059 = 1, 8
          tauself_var_2062 = selffac_var_2044(jl_var_2073, lay_var_2060) * (selfref_var_124(inds_var_2056, ig_var_2059) + selffrac_var_2045(jl_var_2073, lay_var_2060) * (selfref_var_124(inds_var_2056 + 1, ig_var_2059) - selfref_var_124(inds_var_2056, ig_var_2059)))
          taufor_var_2061 = forfac_var_2049(jl_var_2073, lay_var_2060) * (forref_var_125(indf_var_2057, ig_var_2059) + forfrac_var_2050(jl_var_2073, lay_var_2060) * (forref_var_125(indf_var_2057 + 1, ig_var_2059) - forref_var_125(indf_var_2057, ig_var_2059)))
          tauo2 = scaleo2 * (ka_mo2(indm_var_2058, ig_var_2059) + minorfrac_var_2051(jl_var_2073, lay_var_2060) * (ka_mo2(indm_var_2058 + 1, ig_var_2059) - ka_mo2(indm_var_2058, ig_var_2059)))
          taug_var_2033(jl_var_2073, 114 + ig_var_2059, lay_var_2060) = colh2o_var_2042(jl_var_2073, lay_var_2060) * (fac00_var_2035(jl_var_2073, lay_var_2060) * absa_var_122(ind0_var_2054, ig_var_2059) + fac10_var_2037(jl_var_2073, lay_var_2060) * absa_var_122(ind0_var_2054 + 1, ig_var_2059) + fac01_var_2036(jl_var_2073, lay_var_2060) * absa_var_122(ind1_var_2055, ig_var_2059) + fac11_var_2038(jl_var_2073, lay_var_2060) * absa_var_122(ind1_var_2055 + 1, ig_var_2059)) + tauself_var_2062 + taufor_var_2061 + tauo2
          fracs_var_2047(jl_var_2073, 114 + ig_var_2059, lay_var_2060) = fracrefa_var_120(ig_var_2059)
        END DO
      END DO
      ixc0_var_2070 = kfdia_var_2031 - kidia_var_2030 + 1 - ixc0_var_2070
      DO ixp_var_2071 = 1, ixc0_var_2070
        jl_var_2073 = ixhigh_var_2067(ixp_var_2071, lay_var_2060)
        ind0_var_2054 = ((jp_var_2039(jl_var_2073, lay_var_2060) - 13) * 5 + (jt_var_2040(jl_var_2073, lay_var_2060) - 1)) * nspb_var_217(11) + 1
        ind1_var_2055 = ((jp_var_2039(jl_var_2073, lay_var_2060) - 12) * 5 + (jt1_var_2041(jl_var_2073, lay_var_2060) - 1)) * nspb_var_217(11) + 1
        indf_var_2057 = indfor_var_2048(jl_var_2073, lay_var_2060)
        indm_var_2058 = indminor_var_2052(jl_var_2073, lay_var_2060)
        scaleo2 = colo2(jl_var_2073, lay_var_2060) * scaleminor_var_2053(jl_var_2073, lay_var_2060)
        DO ig_var_2059 = 1, 8
          taufor_var_2061 = forfac_var_2049(jl_var_2073, lay_var_2060) * (forref_var_125(indf_var_2057, ig_var_2059) + forfrac_var_2050(jl_var_2073, lay_var_2060) * (forref_var_125(indf_var_2057 + 1, ig_var_2059) - forref_var_125(indf_var_2057, ig_var_2059)))
          tauo2 = scaleo2 * (kb_mo2(indm_var_2058, ig_var_2059) + minorfrac_var_2051(jl_var_2073, lay_var_2060) * (kb_mo2(indm_var_2058 + 1, ig_var_2059) - kb_mo2(indm_var_2058, ig_var_2059)))
          taug_var_2033(jl_var_2073, 114 + ig_var_2059, lay_var_2060) = colh2o_var_2042(jl_var_2073, lay_var_2060) * (fac00_var_2035(jl_var_2073, lay_var_2060) * absb_var_123(ind0_var_2054, ig_var_2059) + fac10_var_2037(jl_var_2073, lay_var_2060) * absb_var_123(ind0_var_2054 + 1, ig_var_2059) + fac01_var_2036(jl_var_2073, lay_var_2060) * absb_var_123(ind1_var_2055, ig_var_2059) + fac11_var_2038(jl_var_2073, lay_var_2060) * absb_var_123(ind1_var_2055 + 1, ig_var_2059)) + taufor_var_2061 + tauo2
          fracs_var_2047(jl_var_2073, 114 + ig_var_2059, lay_var_2060) = fracrefb_var_121(ig_var_2059)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol11
SUBROUTINE rrtm_taumol10(kidia_var_2074, kfdia_var_2075, klev_var_2076, taug_var_2077, p_tauaerl_var_2078, fac00_var_2079, fac01_var_2080, fac10_var_2081, fac11_var_2082, forfac_var_2094, forfrac_var_2093, indfor_var_2092, jp_var_2083, jt_var_2084, jt1_var_2085, colh2o_var_2086, laytrop_var_2087, selffac_var_2089, selffrac_var_2090, indself_var_2091, fracs_var_2088)
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta10, ONLY: absa_var_116, absb_var_117, forref_var_119, fracrefa_var_114, fracrefb_var_115, selfref_var_118
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
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2086(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2087(kidia_var_2074 : kfdia_var_2075)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2088(kidia_var_2074 : kfdia_var_2075, 140, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2089(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2090(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2091(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2092(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2093(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2094(kidia_var_2074 : kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4) :: ind0_var_2095, ind1_var_2096
  INTEGER(KIND = 4) :: inds_var_2097, indf_var_2098
  INTEGER(KIND = 4) :: ig_var_2099, lay_var_2100
  REAL(KIND = 8) :: taufor_var_2101, tauself_var_2102
  INTEGER(KIND = 4) :: laytrop_min_var_2103, laytrop_max_var_2104
  INTEGER(KIND = 4) :: ixc_var_2105(klev_var_2076), ixlow_var_2106(kfdia_var_2075, klev_var_2076), ixhigh_var_2107(kfdia_var_2075, klev_var_2076)
  INTEGER(KIND = 4) :: ich_var_2108, icl_var_2109, ixc0_var_2110, ixp_var_2111, jc_var_2112, jl_var_2113
  laytrop_min_var_2103 = MINVAL(laytrop_var_2087)
  laytrop_max_var_2104 = MAXVAL(laytrop_var_2087)
  ixlow_var_2106 = 0
  ixhigh_var_2107 = 0
  ixc_var_2105 = 0
  DO lay_var_2100 = laytrop_min_var_2103 + 1, laytrop_max_var_2104
    icl_var_2109 = 0
    ich_var_2108 = 0
    DO jc_var_2112 = kidia_var_2074, kfdia_var_2075
      IF (lay_var_2100 <= laytrop_var_2087(jc_var_2112)) THEN
        icl_var_2109 = icl_var_2109 + 1
        ixlow_var_2106(icl_var_2109, lay_var_2100) = jc_var_2112
      ELSE
        ich_var_2108 = ich_var_2108 + 1
        ixhigh_var_2107(ich_var_2108, lay_var_2100) = jc_var_2112
      END IF
    END DO
    ixc_var_2105(lay_var_2100) = icl_var_2109
  END DO
  DO lay_var_2100 = 1, laytrop_min_var_2103
    DO jl_var_2113 = kidia_var_2074, kfdia_var_2075
      ind0_var_2095 = ((jp_var_2083(jl_var_2113, lay_var_2100) - 1) * 5 + (jt_var_2084(jl_var_2113, lay_var_2100) - 1)) * nspa_var_216(10) + 1
      ind1_var_2096 = (jp_var_2083(jl_var_2113, lay_var_2100) * 5 + (jt1_var_2085(jl_var_2113, lay_var_2100) - 1)) * nspa_var_216(10) + 1
      inds_var_2097 = indself_var_2091(jl_var_2113, lay_var_2100)
      indf_var_2098 = indfor_var_2092(jl_var_2113, lay_var_2100)
      DO ig_var_2099 = 1, 6
        tauself_var_2102 = selffac_var_2089(jl_var_2113, lay_var_2100) * (selfref_var_118(inds_var_2097, ig_var_2099) + selffrac_var_2090(jl_var_2113, lay_var_2100) * (selfref_var_118(inds_var_2097 + 1, ig_var_2099) - selfref_var_118(inds_var_2097, ig_var_2099)))
        taufor_var_2101 = forfac_var_2094(jl_var_2113, lay_var_2100) * (forref_var_119(indf_var_2098, ig_var_2099) + forfrac_var_2093(jl_var_2113, lay_var_2100) * (forref_var_119(indf_var_2098 + 1, ig_var_2099) - forref_var_119(indf_var_2098, ig_var_2099)))
        taug_var_2077(jl_var_2113, 108 + ig_var_2099, lay_var_2100) = colh2o_var_2086(jl_var_2113, lay_var_2100) * (fac00_var_2079(jl_var_2113, lay_var_2100) * absa_var_116(ind0_var_2095, ig_var_2099) + fac10_var_2081(jl_var_2113, lay_var_2100) * absa_var_116(ind0_var_2095 + 1, ig_var_2099) + fac01_var_2080(jl_var_2113, lay_var_2100) * absa_var_116(ind1_var_2096, ig_var_2099) + fac11_var_2082(jl_var_2113, lay_var_2100) * absa_var_116(ind1_var_2096 + 1, ig_var_2099)) + tauself_var_2102 + taufor_var_2101
        fracs_var_2088(jl_var_2113, 108 + ig_var_2099, lay_var_2100) = fracrefa_var_114(ig_var_2099)
      END DO
    END DO
  END DO
  DO lay_var_2100 = laytrop_max_var_2104 + 1, klev_var_2076
    DO jl_var_2113 = kidia_var_2074, kfdia_var_2075
      ind0_var_2095 = ((jp_var_2083(jl_var_2113, lay_var_2100) - 13) * 5 + (jt_var_2084(jl_var_2113, lay_var_2100) - 1)) * nspb_var_217(10) + 1
      ind1_var_2096 = ((jp_var_2083(jl_var_2113, lay_var_2100) - 12) * 5 + (jt1_var_2085(jl_var_2113, lay_var_2100) - 1)) * nspb_var_217(10) + 1
      indf_var_2098 = indfor_var_2092(jl_var_2113, lay_var_2100)
      DO ig_var_2099 = 1, 6
        taufor_var_2101 = forfac_var_2094(jl_var_2113, lay_var_2100) * (forref_var_119(indf_var_2098, ig_var_2099) + forfrac_var_2093(jl_var_2113, lay_var_2100) * (forref_var_119(indf_var_2098 + 1, ig_var_2099) - forref_var_119(indf_var_2098, ig_var_2099)))
        taug_var_2077(jl_var_2113, 108 + ig_var_2099, lay_var_2100) = colh2o_var_2086(jl_var_2113, lay_var_2100) * (fac00_var_2079(jl_var_2113, lay_var_2100) * absb_var_117(ind0_var_2095, ig_var_2099) + fac10_var_2081(jl_var_2113, lay_var_2100) * absb_var_117(ind0_var_2095 + 1, ig_var_2099) + fac01_var_2080(jl_var_2113, lay_var_2100) * absb_var_117(ind1_var_2096, ig_var_2099) + fac11_var_2082(jl_var_2113, lay_var_2100) * absb_var_117(ind1_var_2096 + 1, ig_var_2099)) + taufor_var_2101
        fracs_var_2088(jl_var_2113, 108 + ig_var_2099, lay_var_2100) = fracrefb_var_115(ig_var_2099)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2104 /= laytrop_min_var_2103) THEN
    DO lay_var_2100 = laytrop_min_var_2103 + 1, laytrop_max_var_2104
      ixc0_var_2110 = ixc_var_2105(lay_var_2100)
      DO ixp_var_2111 = 1, ixc0_var_2110
        jl_var_2113 = ixlow_var_2106(ixp_var_2111, lay_var_2100)
        ind0_var_2095 = ((jp_var_2083(jl_var_2113, lay_var_2100) - 1) * 5 + (jt_var_2084(jl_var_2113, lay_var_2100) - 1)) * nspa_var_216(10) + 1
        ind1_var_2096 = (jp_var_2083(jl_var_2113, lay_var_2100) * 5 + (jt1_var_2085(jl_var_2113, lay_var_2100) - 1)) * nspa_var_216(10) + 1
        inds_var_2097 = indself_var_2091(jl_var_2113, lay_var_2100)
        indf_var_2098 = indfor_var_2092(jl_var_2113, lay_var_2100)
        DO ig_var_2099 = 1, 6
          tauself_var_2102 = selffac_var_2089(jl_var_2113, lay_var_2100) * (selfref_var_118(inds_var_2097, ig_var_2099) + selffrac_var_2090(jl_var_2113, lay_var_2100) * (selfref_var_118(inds_var_2097 + 1, ig_var_2099) - selfref_var_118(inds_var_2097, ig_var_2099)))
          taufor_var_2101 = forfac_var_2094(jl_var_2113, lay_var_2100) * (forref_var_119(indf_var_2098, ig_var_2099) + forfrac_var_2093(jl_var_2113, lay_var_2100) * (forref_var_119(indf_var_2098 + 1, ig_var_2099) - forref_var_119(indf_var_2098, ig_var_2099)))
          taug_var_2077(jl_var_2113, 108 + ig_var_2099, lay_var_2100) = colh2o_var_2086(jl_var_2113, lay_var_2100) * (fac00_var_2079(jl_var_2113, lay_var_2100) * absa_var_116(ind0_var_2095, ig_var_2099) + fac10_var_2081(jl_var_2113, lay_var_2100) * absa_var_116(ind0_var_2095 + 1, ig_var_2099) + fac01_var_2080(jl_var_2113, lay_var_2100) * absa_var_116(ind1_var_2096, ig_var_2099) + fac11_var_2082(jl_var_2113, lay_var_2100) * absa_var_116(ind1_var_2096 + 1, ig_var_2099)) + tauself_var_2102 + taufor_var_2101
          fracs_var_2088(jl_var_2113, 108 + ig_var_2099, lay_var_2100) = fracrefa_var_114(ig_var_2099)
        END DO
      END DO
      ixc0_var_2110 = kfdia_var_2075 - kidia_var_2074 + 1 - ixc0_var_2110
      DO ixp_var_2111 = 1, ixc0_var_2110
        jl_var_2113 = ixhigh_var_2107(ixp_var_2111, lay_var_2100)
        ind0_var_2095 = ((jp_var_2083(jl_var_2113, lay_var_2100) - 13) * 5 + (jt_var_2084(jl_var_2113, lay_var_2100) - 1)) * nspb_var_217(10) + 1
        ind1_var_2096 = ((jp_var_2083(jl_var_2113, lay_var_2100) - 12) * 5 + (jt1_var_2085(jl_var_2113, lay_var_2100) - 1)) * nspb_var_217(10) + 1
        indf_var_2098 = indfor_var_2092(jl_var_2113, lay_var_2100)
        DO ig_var_2099 = 1, 6
          taufor_var_2101 = forfac_var_2094(jl_var_2113, lay_var_2100) * (forref_var_119(indf_var_2098, ig_var_2099) + forfrac_var_2093(jl_var_2113, lay_var_2100) * (forref_var_119(indf_var_2098 + 1, ig_var_2099) - forref_var_119(indf_var_2098, ig_var_2099)))
          taug_var_2077(jl_var_2113, 108 + ig_var_2099, lay_var_2100) = colh2o_var_2086(jl_var_2113, lay_var_2100) * (fac00_var_2079(jl_var_2113, lay_var_2100) * absb_var_117(ind0_var_2095, ig_var_2099) + fac10_var_2081(jl_var_2113, lay_var_2100) * absb_var_117(ind0_var_2095 + 1, ig_var_2099) + fac01_var_2080(jl_var_2113, lay_var_2100) * absb_var_117(ind1_var_2096, ig_var_2099) + fac11_var_2082(jl_var_2113, lay_var_2100) * absb_var_117(ind1_var_2096 + 1, ig_var_2099)) + taufor_var_2101
          fracs_var_2088(jl_var_2113, 108 + ig_var_2099, lay_var_2100) = fracrefb_var_115(ig_var_2099)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol10
SUBROUTINE rrtm_taumol8(kidia_var_2114, kfdia_var_2115, klev_var_2116, taug_var_2117, wx_var_2118, p_tauaerl_var_2119, fac00_var_2120, fac01_var_2121, fac10_var_2122, fac11_var_2123, forfac_var_2139, forfrac_var_2138, indfor_var_2137, jp_var_2124, jt_var_2125, jt1_var_2126, colh2o_var_2127, colo3_var_2128, coln2o_var_2129, colco2_var_2130, coldry_var_2131, laytrop_var_2132, selffac_var_2133, selffrac_var_2134, indself_var_2135, fracs_var_2136, minorfrac_var_2140, indminor_var_2141)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_216, nspb_var_217
  USE yoerrta8, ONLY: absa_var_197, absb_var_198, cfc12_var_196, cfc22adj, forref_var_205, fracrefa_var_194, fracrefb_var_195, ka_mco2_var_199, ka_mn2o_var_200, ka_mo3_var_201, kb_mco2_var_202, kb_mn2o_var_203, selfref_var_204
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2114
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2115
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2116
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2117(kidia_var_2114 : kfdia_var_2115, 140, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: wx_var_2118(kidia_var_2114 : kfdia_var_2115, 4, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2119(kidia_var_2114 : kfdia_var_2115, klev_var_2116, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2120(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2121(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2122(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2123(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2124(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2125(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2126(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2127(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_2128(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_2129(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2130(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_2131(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2132(kidia_var_2114 : kfdia_var_2115)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2133(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2134(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2135(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2136(kidia_var_2114 : kfdia_var_2115, 140, klev_var_2116)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2137(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2138(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2139(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2140(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2141(kidia_var_2114 : kfdia_var_2115, klev_var_2116)
  INTEGER(KIND = 4) :: ind0_var_2142, ind1_var_2143, inds_var_2144, indf_var_2145, indm_var_2146
  INTEGER(KIND = 4) :: ig_var_2147, lay_var_2148
  REAL(KIND = 8) :: chi_co2_var_2149, ratco2_var_2150, adjfac_var_2151, adjcolco2_var_2152
  REAL(KIND = 8) :: taufor_var_2153, tauself_var_2154, abso3_var_2155, absco2_var_2156, absn2o_var_2157
  INTEGER(KIND = 4) :: laytrop_min_var_2158, laytrop_max_var_2159
  INTEGER(KIND = 4) :: ixc_var_2160(klev_var_2116), ixlow_var_2161(kfdia_var_2115, klev_var_2116), ixhigh_var_2162(kfdia_var_2115, klev_var_2116)
  INTEGER(KIND = 4) :: ich_var_2163, icl_var_2164, ixc0_var_2165, ixp_var_2166, jc_var_2167, jl_var_2168
  laytrop_min_var_2158 = MINVAL(laytrop_var_2132)
  laytrop_max_var_2159 = MAXVAL(laytrop_var_2132)
  ixlow_var_2161 = 0
  ixhigh_var_2162 = 0
  ixc_var_2160 = 0
  DO lay_var_2148 = laytrop_min_var_2158 + 1, laytrop_max_var_2159
    icl_var_2164 = 0
    ich_var_2163 = 0
    DO jc_var_2167 = kidia_var_2114, kfdia_var_2115
      IF (lay_var_2148 <= laytrop_var_2132(jc_var_2167)) THEN
        icl_var_2164 = icl_var_2164 + 1
        ixlow_var_2161(icl_var_2164, lay_var_2148) = jc_var_2167
      ELSE
        ich_var_2163 = ich_var_2163 + 1
        ixhigh_var_2162(ich_var_2163, lay_var_2148) = jc_var_2167
      END IF
    END DO
    ixc_var_2160(lay_var_2148) = icl_var_2164
  END DO
  DO lay_var_2148 = 1, laytrop_min_var_2158
    DO jl_var_2168 = kidia_var_2114, kfdia_var_2115
      chi_co2_var_2149 = colco2_var_2130(jl_var_2168, lay_var_2148) / (coldry_var_2131(jl_var_2168, lay_var_2148))
      ratco2_var_2150 = 1D+20 * chi_co2_var_2149 / chi_mls(2, jp_var_2124(jl_var_2168, lay_var_2148) + 1)
      IF (ratco2_var_2150 .GT. 3.0D0) THEN
        adjfac_var_2151 = 2.0D0 + (ratco2_var_2150 - 2.0D0) ** 0.65D0
        adjcolco2_var_2152 = adjfac_var_2151 * chi_mls(2, jp_var_2124(jl_var_2168, lay_var_2148) + 1) * coldry_var_2131(jl_var_2168, lay_var_2148) * 1D-20
      ELSE
        adjcolco2_var_2152 = colco2_var_2130(jl_var_2168, lay_var_2148)
      END IF
      ind0_var_2142 = ((jp_var_2124(jl_var_2168, lay_var_2148) - 1) * 5 + (jt_var_2125(jl_var_2168, lay_var_2148) - 1)) * nspa_var_216(8) + 1
      ind1_var_2143 = (jp_var_2124(jl_var_2168, lay_var_2148) * 5 + (jt1_var_2126(jl_var_2168, lay_var_2148) - 1)) * nspa_var_216(8) + 1
      inds_var_2144 = indself_var_2135(jl_var_2168, lay_var_2148)
      indf_var_2145 = indfor_var_2137(jl_var_2168, lay_var_2148)
      indm_var_2146 = indminor_var_2141(jl_var_2168, lay_var_2148)
      DO ig_var_2147 = 1, 8
        tauself_var_2154 = selffac_var_2133(jl_var_2168, lay_var_2148) * (selfref_var_204(inds_var_2144, ig_var_2147) + selffrac_var_2134(jl_var_2168, lay_var_2148) * (selfref_var_204(inds_var_2144 + 1, ig_var_2147) - selfref_var_204(inds_var_2144, ig_var_2147)))
        taufor_var_2153 = forfac_var_2139(jl_var_2168, lay_var_2148) * (forref_var_205(indf_var_2145, ig_var_2147) + forfrac_var_2138(jl_var_2168, lay_var_2148) * (forref_var_205(indf_var_2145 + 1, ig_var_2147) - forref_var_205(indf_var_2145, ig_var_2147)))
        absco2_var_2156 = (ka_mco2_var_199(indm_var_2146, ig_var_2147) + minorfrac_var_2140(jl_var_2168, lay_var_2148) * (ka_mco2_var_199(indm_var_2146 + 1, ig_var_2147) - ka_mco2_var_199(indm_var_2146, ig_var_2147)))
        abso3_var_2155 = (ka_mo3_var_201(indm_var_2146, ig_var_2147) + minorfrac_var_2140(jl_var_2168, lay_var_2148) * (ka_mo3_var_201(indm_var_2146 + 1, ig_var_2147) - ka_mo3_var_201(indm_var_2146, ig_var_2147)))
        absn2o_var_2157 = (ka_mn2o_var_200(indm_var_2146, ig_var_2147) + minorfrac_var_2140(jl_var_2168, lay_var_2148) * (ka_mn2o_var_200(indm_var_2146 + 1, ig_var_2147) - ka_mn2o_var_200(indm_var_2146, ig_var_2147)))
        taug_var_2117(jl_var_2168, 88 + ig_var_2147, lay_var_2148) = colh2o_var_2127(jl_var_2168, lay_var_2148) * (fac00_var_2120(jl_var_2168, lay_var_2148) * absa_var_197(ind0_var_2142, ig_var_2147) + fac10_var_2122(jl_var_2168, lay_var_2148) * absa_var_197(ind0_var_2142 + 1, ig_var_2147) + fac01_var_2121(jl_var_2168, lay_var_2148) * absa_var_197(ind1_var_2143, ig_var_2147) + fac11_var_2123(jl_var_2168, lay_var_2148) * absa_var_197(ind1_var_2143 + 1, ig_var_2147)) + tauself_var_2154 + taufor_var_2153 + adjcolco2_var_2152 * absco2_var_2156 + colo3_var_2128(jl_var_2168, lay_var_2148) * abso3_var_2155 + coln2o_var_2129(jl_var_2168, lay_var_2148) * absn2o_var_2157 + wx_var_2118(jl_var_2168, 3, lay_var_2148) * cfc12_var_196(ig_var_2147) + wx_var_2118(jl_var_2168, 4, lay_var_2148) * cfc22adj(ig_var_2147)
        fracs_var_2136(jl_var_2168, 88 + ig_var_2147, lay_var_2148) = fracrefa_var_194(ig_var_2147)
      END DO
    END DO
  END DO
  DO lay_var_2148 = laytrop_max_var_2159 + 1, klev_var_2116
    DO jl_var_2168 = kidia_var_2114, kfdia_var_2115
      chi_co2_var_2149 = colco2_var_2130(jl_var_2168, lay_var_2148) / coldry_var_2131(jl_var_2168, lay_var_2148)
      ratco2_var_2150 = 1D+20 * chi_co2_var_2149 / chi_mls(2, jp_var_2124(jl_var_2168, lay_var_2148) + 1)
      IF (ratco2_var_2150 .GT. 3.0D0) THEN
        adjfac_var_2151 = 2.0D0 + (ratco2_var_2150 - 2.0D0) ** 0.65D0
        adjcolco2_var_2152 = adjfac_var_2151 * chi_mls(2, jp_var_2124(jl_var_2168, lay_var_2148) + 1) * coldry_var_2131(jl_var_2168, lay_var_2148) * 1D-20
      ELSE
        adjcolco2_var_2152 = colco2_var_2130(jl_var_2168, lay_var_2148)
      END IF
      ind0_var_2142 = ((jp_var_2124(jl_var_2168, lay_var_2148) - 13) * 5 + (jt_var_2125(jl_var_2168, lay_var_2148) - 1)) * nspb_var_217(8) + 1
      ind1_var_2143 = ((jp_var_2124(jl_var_2168, lay_var_2148) - 12) * 5 + (jt1_var_2126(jl_var_2168, lay_var_2148) - 1)) * nspb_var_217(8) + 1
      indm_var_2146 = indminor_var_2141(jl_var_2168, lay_var_2148)
      DO ig_var_2147 = 1, 8
        absco2_var_2156 = (kb_mco2_var_202(indm_var_2146, ig_var_2147) + minorfrac_var_2140(jl_var_2168, lay_var_2148) * (kb_mco2_var_202(indm_var_2146 + 1, ig_var_2147) - kb_mco2_var_202(indm_var_2146, ig_var_2147)))
        absn2o_var_2157 = (kb_mn2o_var_203(indm_var_2146, ig_var_2147) + minorfrac_var_2140(jl_var_2168, lay_var_2148) * (kb_mn2o_var_203(indm_var_2146 + 1, ig_var_2147) - kb_mn2o_var_203(indm_var_2146, ig_var_2147)))
        taug_var_2117(jl_var_2168, 88 + ig_var_2147, lay_var_2148) = colo3_var_2128(jl_var_2168, lay_var_2148) * (fac00_var_2120(jl_var_2168, lay_var_2148) * absb_var_198(ind0_var_2142, ig_var_2147) + fac10_var_2122(jl_var_2168, lay_var_2148) * absb_var_198(ind0_var_2142 + 1, ig_var_2147) + fac01_var_2121(jl_var_2168, lay_var_2148) * absb_var_198(ind1_var_2143, ig_var_2147) + fac11_var_2123(jl_var_2168, lay_var_2148) * absb_var_198(ind1_var_2143 + 1, ig_var_2147)) + adjcolco2_var_2152 * absco2_var_2156 + coln2o_var_2129(jl_var_2168, lay_var_2148) * absn2o_var_2157 + wx_var_2118(jl_var_2168, 3, lay_var_2148) * cfc12_var_196(ig_var_2147) + wx_var_2118(jl_var_2168, 4, lay_var_2148) * cfc22adj(ig_var_2147)
        fracs_var_2136(jl_var_2168, 88 + ig_var_2147, lay_var_2148) = fracrefb_var_195(ig_var_2147)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2159 /= laytrop_min_var_2158) THEN
    DO lay_var_2148 = laytrop_min_var_2158 + 1, laytrop_max_var_2159
      ixc0_var_2165 = ixc_var_2160(lay_var_2148)
      DO ixp_var_2166 = 1, ixc0_var_2165
        jl_var_2168 = ixlow_var_2161(ixp_var_2166, lay_var_2148)
        chi_co2_var_2149 = colco2_var_2130(jl_var_2168, lay_var_2148) / (coldry_var_2131(jl_var_2168, lay_var_2148))
        ratco2_var_2150 = 1D+20 * chi_co2_var_2149 / chi_mls(2, jp_var_2124(jl_var_2168, lay_var_2148) + 1)
        IF (ratco2_var_2150 .GT. 3.0D0) THEN
          adjfac_var_2151 = 2.0D0 + (ratco2_var_2150 - 2.0D0) ** 0.65D0
          adjcolco2_var_2152 = adjfac_var_2151 * chi_mls(2, jp_var_2124(jl_var_2168, lay_var_2148) + 1) * coldry_var_2131(jl_var_2168, lay_var_2148) * 1D-20
        ELSE
          adjcolco2_var_2152 = colco2_var_2130(jl_var_2168, lay_var_2148)
        END IF
        ind0_var_2142 = ((jp_var_2124(jl_var_2168, lay_var_2148) - 1) * 5 + (jt_var_2125(jl_var_2168, lay_var_2148) - 1)) * nspa_var_216(8) + 1
        ind1_var_2143 = (jp_var_2124(jl_var_2168, lay_var_2148) * 5 + (jt1_var_2126(jl_var_2168, lay_var_2148) - 1)) * nspa_var_216(8) + 1
        inds_var_2144 = indself_var_2135(jl_var_2168, lay_var_2148)
        indf_var_2145 = indfor_var_2137(jl_var_2168, lay_var_2148)
        indm_var_2146 = indminor_var_2141(jl_var_2168, lay_var_2148)
        DO ig_var_2147 = 1, 8
          tauself_var_2154 = selffac_var_2133(jl_var_2168, lay_var_2148) * (selfref_var_204(inds_var_2144, ig_var_2147) + selffrac_var_2134(jl_var_2168, lay_var_2148) * (selfref_var_204(inds_var_2144 + 1, ig_var_2147) - selfref_var_204(inds_var_2144, ig_var_2147)))
          taufor_var_2153 = forfac_var_2139(jl_var_2168, lay_var_2148) * (forref_var_205(indf_var_2145, ig_var_2147) + forfrac_var_2138(jl_var_2168, lay_var_2148) * (forref_var_205(indf_var_2145 + 1, ig_var_2147) - forref_var_205(indf_var_2145, ig_var_2147)))
          absco2_var_2156 = (ka_mco2_var_199(indm_var_2146, ig_var_2147) + minorfrac_var_2140(jl_var_2168, lay_var_2148) * (ka_mco2_var_199(indm_var_2146 + 1, ig_var_2147) - ka_mco2_var_199(indm_var_2146, ig_var_2147)))
          abso3_var_2155 = (ka_mo3_var_201(indm_var_2146, ig_var_2147) + minorfrac_var_2140(jl_var_2168, lay_var_2148) * (ka_mo3_var_201(indm_var_2146 + 1, ig_var_2147) - ka_mo3_var_201(indm_var_2146, ig_var_2147)))
          absn2o_var_2157 = (ka_mn2o_var_200(indm_var_2146, ig_var_2147) + minorfrac_var_2140(jl_var_2168, lay_var_2148) * (ka_mn2o_var_200(indm_var_2146 + 1, ig_var_2147) - ka_mn2o_var_200(indm_var_2146, ig_var_2147)))
          taug_var_2117(jl_var_2168, 88 + ig_var_2147, lay_var_2148) = colh2o_var_2127(jl_var_2168, lay_var_2148) * (fac00_var_2120(jl_var_2168, lay_var_2148) * absa_var_197(ind0_var_2142, ig_var_2147) + fac10_var_2122(jl_var_2168, lay_var_2148) * absa_var_197(ind0_var_2142 + 1, ig_var_2147) + fac01_var_2121(jl_var_2168, lay_var_2148) * absa_var_197(ind1_var_2143, ig_var_2147) + fac11_var_2123(jl_var_2168, lay_var_2148) * absa_var_197(ind1_var_2143 + 1, ig_var_2147)) + tauself_var_2154 + taufor_var_2153 + adjcolco2_var_2152 * absco2_var_2156 + colo3_var_2128(jl_var_2168, lay_var_2148) * abso3_var_2155 + coln2o_var_2129(jl_var_2168, lay_var_2148) * absn2o_var_2157 + wx_var_2118(jl_var_2168, 3, lay_var_2148) * cfc12_var_196(ig_var_2147) + wx_var_2118(jl_var_2168, 4, lay_var_2148) * cfc22adj(ig_var_2147)
          fracs_var_2136(jl_var_2168, 88 + ig_var_2147, lay_var_2148) = fracrefa_var_194(ig_var_2147)
        END DO
      END DO
      ixc0_var_2165 = kfdia_var_2115 - kidia_var_2114 + 1 - ixc0_var_2165
      DO ixp_var_2166 = 1, ixc0_var_2165
        jl_var_2168 = ixhigh_var_2162(ixp_var_2166, lay_var_2148)
        chi_co2_var_2149 = colco2_var_2130(jl_var_2168, lay_var_2148) / coldry_var_2131(jl_var_2168, lay_var_2148)
        ratco2_var_2150 = 1D+20 * chi_co2_var_2149 / chi_mls(2, jp_var_2124(jl_var_2168, lay_var_2148) + 1)
        IF (ratco2_var_2150 .GT. 3.0D0) THEN
          adjfac_var_2151 = 2.0D0 + (ratco2_var_2150 - 2.0D0) ** 0.65D0
          adjcolco2_var_2152 = adjfac_var_2151 * chi_mls(2, jp_var_2124(jl_var_2168, lay_var_2148) + 1) * coldry_var_2131(jl_var_2168, lay_var_2148) * 1D-20
        ELSE
          adjcolco2_var_2152 = colco2_var_2130(jl_var_2168, lay_var_2148)
        END IF
        ind0_var_2142 = ((jp_var_2124(jl_var_2168, lay_var_2148) - 13) * 5 + (jt_var_2125(jl_var_2168, lay_var_2148) - 1)) * nspb_var_217(8) + 1
        ind1_var_2143 = ((jp_var_2124(jl_var_2168, lay_var_2148) - 12) * 5 + (jt1_var_2126(jl_var_2168, lay_var_2148) - 1)) * nspb_var_217(8) + 1
        indm_var_2146 = indminor_var_2141(jl_var_2168, lay_var_2148)
        DO ig_var_2147 = 1, 8
          absco2_var_2156 = (kb_mco2_var_202(indm_var_2146, ig_var_2147) + minorfrac_var_2140(jl_var_2168, lay_var_2148) * (kb_mco2_var_202(indm_var_2146 + 1, ig_var_2147) - kb_mco2_var_202(indm_var_2146, ig_var_2147)))
          absn2o_var_2157 = (kb_mn2o_var_203(indm_var_2146, ig_var_2147) + minorfrac_var_2140(jl_var_2168, lay_var_2148) * (kb_mn2o_var_203(indm_var_2146 + 1, ig_var_2147) - kb_mn2o_var_203(indm_var_2146, ig_var_2147)))
          taug_var_2117(jl_var_2168, 88 + ig_var_2147, lay_var_2148) = colo3_var_2128(jl_var_2168, lay_var_2148) * (fac00_var_2120(jl_var_2168, lay_var_2148) * absb_var_198(ind0_var_2142, ig_var_2147) + fac10_var_2122(jl_var_2168, lay_var_2148) * absb_var_198(ind0_var_2142 + 1, ig_var_2147) + fac01_var_2121(jl_var_2168, lay_var_2148) * absb_var_198(ind1_var_2143, ig_var_2147) + fac11_var_2123(jl_var_2168, lay_var_2148) * absb_var_198(ind1_var_2143 + 1, ig_var_2147)) + adjcolco2_var_2152 * absco2_var_2156 + coln2o_var_2129(jl_var_2168, lay_var_2148) * absn2o_var_2157 + wx_var_2118(jl_var_2168, 3, lay_var_2148) * cfc12_var_196(ig_var_2147) + wx_var_2118(jl_var_2168, 4, lay_var_2148) * cfc22adj(ig_var_2147)
          fracs_var_2136(jl_var_2168, 88 + ig_var_2147, lay_var_2148) = fracrefb_var_195(ig_var_2147)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol8
SUBROUTINE rrtm_prepare_gases(kidia_var_2170, kfdia_var_2171, klon, klev_var_2169, paph, pap, pth, pt, pq, pco2, pch4, pn2o, pno2, pc11, pc12, pc22, pcl4, pozn, pcoldry_var_2172, pwbrodl, pwkl_var_2173, pwx_var_2174, pavel_var_2175, ptavel_var_2176, pz, ptz, kreflect)
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: klon
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2169
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2170, kfdia_var_2171
  REAL(KIND = 8), INTENT(IN) :: paph(klon, klev_var_2169 + 1)
  REAL(KIND = 8), INTENT(IN) :: pap(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(IN) :: pth(klon, klev_var_2169 + 1)
  REAL(KIND = 8), INTENT(IN) :: pt(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(IN) :: pq(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(IN) :: pco2(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(IN) :: pch4(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(IN) :: pn2o(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(IN) :: pno2(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(IN) :: pc11(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(IN) :: pc12(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(IN) :: pc22(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(IN) :: pcl4(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(IN) :: pozn(klon, klev_var_2169)
  REAL(KIND = 8), INTENT(OUT) :: pcoldry_var_2172(kidia_var_2170 : kfdia_var_2171, klev_var_2169)
  REAL(KIND = 8), INTENT(OUT) :: pwbrodl(kidia_var_2170 : kfdia_var_2171, klev_var_2169)
  REAL(KIND = 8), INTENT(OUT) :: pwkl_var_2173(kidia_var_2170 : kfdia_var_2171, 35, klev_var_2169)
  REAL(KIND = 8), INTENT(OUT) :: pwx_var_2174(kidia_var_2170 : kfdia_var_2171, 4, klev_var_2169)
  REAL(KIND = 8), INTENT(OUT) :: pavel_var_2175(kidia_var_2170 : kfdia_var_2171, klev_var_2169)
  REAL(KIND = 8), INTENT(OUT) :: ptavel_var_2176(kidia_var_2170 : kfdia_var_2171, klev_var_2169)
  REAL(KIND = 8), INTENT(OUT) :: pz(kidia_var_2170 : kfdia_var_2171, 0 : klev_var_2169)
  REAL(KIND = 8), INTENT(OUT) :: ptz(kidia_var_2170 : kfdia_var_2171, 0 : klev_var_2169)
  INTEGER(KIND = 4), INTENT(OUT) :: kreflect(kidia_var_2170 : kfdia_var_2171)
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
  INTEGER(KIND = 4) :: iatm, jmol, ixmax, j1, j2, jk_var_2177, jl_var_2178
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
  DO jl_var_2178 = kidia_var_2170, kfdia_var_2171
    kreflect(jl_var_2178) = 0
  END DO
  DO j2 = 1, klev_var_2169
    DO j1 = 1, 35
      DO jl_var_2178 = kidia_var_2170, kfdia_var_2171
        pwkl_var_2173(jl_var_2178, j1, j2) = 0.0D0
      END DO
    END DO
  END DO
  iatm = 0
  ixmax = 4
  DO jl_var_2178 = kidia_var_2170, kfdia_var_2171
    pz(jl_var_2178, 0) = paph(jl_var_2178, klev_var_2169 + 1) / 100.0D0
    ptz(jl_var_2178, 0) = pth(jl_var_2178, klev_var_2169 + 1)
  END DO
  DO jk_var_2177 = 1, klev_var_2169
    DO jl_var_2178 = kidia_var_2170, kfdia_var_2171
      pavel_var_2175(jl_var_2178, jk_var_2177) = pap(jl_var_2178, klev_var_2169 - jk_var_2177 + 1) / 100.0D0
      ptavel_var_2176(jl_var_2178, jk_var_2177) = pt(jl_var_2178, klev_var_2169 - jk_var_2177 + 1)
      pz(jl_var_2178, jk_var_2177) = paph(jl_var_2178, klev_var_2169 - jk_var_2177 + 1) / 100.0D0
      ptz(jl_var_2178, jk_var_2177) = pth(jl_var_2178, klev_var_2169 - jk_var_2177 + 1)
      pwkl_var_2173(jl_var_2178, 1, jk_var_2177) = MAX(pq(jl_var_2178, klev_var_2169 - jk_var_2177 + 1), 1E-15) * zamd / zamw
      pwkl_var_2173(jl_var_2178, 2, jk_var_2177) = pco2(jl_var_2178, klev_var_2169 - jk_var_2177 + 1) * zamd / zamco2
      pwkl_var_2173(jl_var_2178, 3, jk_var_2177) = pozn(jl_var_2178, klev_var_2169 - jk_var_2177 + 1) * zamd / zamo
      pwkl_var_2173(jl_var_2178, 4, jk_var_2177) = pn2o(jl_var_2178, klev_var_2169 - jk_var_2177 + 1) * zamd / zamn2o
      pwkl_var_2173(jl_var_2178, 6, jk_var_2177) = pch4(jl_var_2178, klev_var_2169 - jk_var_2177 + 1) * zamd / zamch4
      pwkl_var_2173(jl_var_2178, 7, jk_var_2177) = 0.209488D0
      zamm = (1.0D0 - pwkl_var_2173(jl_var_2178, 1, jk_var_2177)) * zamd + pwkl_var_2173(jl_var_2178, 1, jk_var_2177) * zamw
      pcoldry_var_2172(jl_var_2178, jk_var_2177) = (pz(jl_var_2178, jk_var_2177 - 1) - pz(jl_var_2178, jk_var_2177)) * 1000.0D0 * zavgdro / (zgravit * zamm * (1.0D0 + pwkl_var_2173(jl_var_2178, 1, jk_var_2177)))
    END DO
  END DO
  DO j2 = 1, klev_var_2169
    DO j1 = 1, 4
      DO jl_var_2178 = kidia_var_2170, kfdia_var_2171
        pwx_var_2174(jl_var_2178, j1, j2) = 0.0D0
      END DO
    END DO
  END DO
  DO jk_var_2177 = 1, klev_var_2169
    DO jl_var_2178 = kidia_var_2170, kfdia_var_2171
      pwx_var_2174(jl_var_2178, 1, jk_var_2177) = pcl4(jl_var_2178, klev_var_2169 - jk_var_2177 + 1) * zamd / zamcl4
      pwx_var_2174(jl_var_2178, 2, jk_var_2177) = pc11(jl_var_2178, klev_var_2169 - jk_var_2177 + 1) * zamd / zamc11
      pwx_var_2174(jl_var_2178, 3, jk_var_2177) = pc12(jl_var_2178, klev_var_2169 - jk_var_2177 + 1) * zamd / zamc12
      pwx_var_2174(jl_var_2178, 4, jk_var_2177) = pc22(jl_var_2178, klev_var_2169 - jk_var_2177 + 1) * zamd / zamc22
      pwx_var_2174(jl_var_2178, 1, jk_var_2177) = pcoldry_var_2172(jl_var_2178, jk_var_2177) * pwx_var_2174(jl_var_2178, 1, jk_var_2177) * 1D-20
      pwx_var_2174(jl_var_2178, 2, jk_var_2177) = pcoldry_var_2172(jl_var_2178, jk_var_2177) * pwx_var_2174(jl_var_2178, 2, jk_var_2177) * 1D-20
      pwx_var_2174(jl_var_2178, 3, jk_var_2177) = pcoldry_var_2172(jl_var_2178, jk_var_2177) * pwx_var_2174(jl_var_2178, 3, jk_var_2177) * 1D-20
      pwx_var_2174(jl_var_2178, 4, jk_var_2177) = pcoldry_var_2172(jl_var_2178, jk_var_2177) * pwx_var_2174(jl_var_2178, 4, jk_var_2177) * 1D-20
      zsummol = 0.0D0
      DO jmol = 2, 7
        zsummol = zsummol + pwkl_var_2173(jl_var_2178, jmol, jk_var_2177)
      END DO
      pwbrodl(jl_var_2178, jk_var_2177) = pcoldry_var_2172(jl_var_2178, jk_var_2177) * (1.0D0 - zsummol)
      DO jmol = 1, 7
        pwkl_var_2173(jl_var_2178, jmol, jk_var_2177) = pcoldry_var_2172(jl_var_2178, jk_var_2177) * pwkl_var_2173(jl_var_2178, jmol, jk_var_2177)
      END DO
    END DO
  END DO
END SUBROUTINE rrtm_prepare_gases
SUBROUTINE srtm_gas_optical_depth(kidia_var_2179, kfdia_var_2180, klev_var_2181, poneminus_var_2182, prmu0_var_2183, klaytrop_var_2184, pcolch4_var_2185, pcolco2_var_2186, pcolh2o_var_2187, pcolmol_var_2188, pcolo2_var_2189, pcolo3_var_2190, pforfac_var_2191, pforfrac_var_2192, kindfor_var_2193, pselffac_var_2194, pselffrac_var_2195, kindself_var_2196, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, pod_var_2204, pssa, pincsol)
  USE yoesrtwn, ONLY: ngc
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2179, kfdia_var_2180
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2181
  REAL(KIND = 8), INTENT(IN) :: poneminus_var_2182(kidia_var_2179 : kfdia_var_2180)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_2183(kidia_var_2179 : kfdia_var_2180)
  INTEGER(KIND = 4), INTENT(IN) :: klaytrop_var_2184(kidia_var_2179 : kfdia_var_2180)
  REAL(KIND = 8), INTENT(IN) :: pcolch4_var_2185(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pcolco2_var_2186(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pcolh2o_var_2187(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pcolmol_var_2188(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pcolo2_var_2189(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pcolo3_var_2190(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pforfac_var_2191(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pforfrac_var_2192(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  INTEGER(KIND = 4), INTENT(IN) :: kindfor_var_2193(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pselffac_var_2194(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pselffrac_var_2195(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  INTEGER(KIND = 4), INTENT(IN) :: kindself_var_2196(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pfac00_var_2197(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pfac01_var_2198(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pfac10_var_2199(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(IN) :: pfac11_var_2200(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  INTEGER(KIND = 4), INTENT(IN) :: kjp_var_2201(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  INTEGER(KIND = 4), INTENT(IN) :: kjt_var_2202(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  INTEGER(KIND = 4), INTENT(IN) :: kjt1_var_2203(kidia_var_2179 : kfdia_var_2180, klev_var_2181)
  REAL(KIND = 8), INTENT(INOUT) :: pod_var_2204(kidia_var_2179 : kfdia_var_2180, klev_var_2181, 112)
  REAL(KIND = 8), INTENT(INOUT) :: pssa(kidia_var_2179 : kfdia_var_2180, klev_var_2181, 112)
  REAL(KIND = 8), INTENT(INOUT) :: pincsol(kidia_var_2179 : kfdia_var_2180, 112)
  INTEGER(KIND = 4) :: ib1, ib2, ibm, igt, iw(kidia_var_2179 : kfdia_var_2180), jb_var_2205, jg_var_2206, jk_var_2207, jl_var_2208, icount
  REAL(KIND = 8) :: ztaug(kidia_var_2179 : kfdia_var_2180, klev_var_2181, 16)
  REAL(KIND = 8) :: ztaur(kidia_var_2179 : kfdia_var_2180, klev_var_2181, 16)
  REAL(KIND = 8) :: zsflxzen(kidia_var_2179 : kfdia_var_2180, 16)
  ib1 = 16
  ib2 = 29
  icount = 0
  DO jl_var_2208 = kidia_var_2179, kfdia_var_2180
    IF (prmu0_var_2183(jl_var_2208) > 0.0D0) THEN
      icount = icount + 1
      iw(jl_var_2208) = 0
    END IF
  END DO
  IF (icount /= 0) THEN
    DO jb_var_2205 = ib1, ib2
      ibm = jb_var_2205 - 15
      igt = ngc(ibm)
      IF (jb_var_2205 == 16) THEN
        CALL srtm_taumol16(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, poneminus_var_2182, pcolh2o_var_2187, pcolch4_var_2185, pcolmol_var_2188, klaytrop_var_2184, pselffac_var_2194, pselffrac_var_2195, kindself_var_2196, pforfac_var_2191, pforfrac_var_2192, kindfor_var_2193, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 17) THEN
        CALL srtm_taumol17(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, poneminus_var_2182, pcolh2o_var_2187, pcolco2_var_2186, pcolmol_var_2188, klaytrop_var_2184, pselffac_var_2194, pselffrac_var_2195, kindself_var_2196, pforfac_var_2191, pforfrac_var_2192, kindfor_var_2193, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 18) THEN
        CALL srtm_taumol18(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, poneminus_var_2182, pcolh2o_var_2187, pcolch4_var_2185, pcolmol_var_2188, klaytrop_var_2184, pselffac_var_2194, pselffrac_var_2195, kindself_var_2196, pforfac_var_2191, pforfrac_var_2192, kindfor_var_2193, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 19) THEN
        CALL srtm_taumol19(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, poneminus_var_2182, pcolh2o_var_2187, pcolco2_var_2186, pcolmol_var_2188, klaytrop_var_2184, pselffac_var_2194, pselffrac_var_2195, kindself_var_2196, pforfac_var_2191, pforfrac_var_2192, kindfor_var_2193, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 20) THEN
        CALL srtm_taumol20(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, pcolh2o_var_2187, pcolch4_var_2185, pcolmol_var_2188, klaytrop_var_2184, pselffac_var_2194, pselffrac_var_2195, kindself_var_2196, pforfac_var_2191, pforfrac_var_2192, kindfor_var_2193, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 21) THEN
        CALL srtm_taumol21(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, poneminus_var_2182, pcolh2o_var_2187, pcolco2_var_2186, pcolmol_var_2188, klaytrop_var_2184, pselffac_var_2194, pselffrac_var_2195, kindself_var_2196, pforfac_var_2191, pforfrac_var_2192, kindfor_var_2193, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 22) THEN
        CALL srtm_taumol22(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, poneminus_var_2182, pcolh2o_var_2187, pcolmol_var_2188, pcolo2_var_2189, klaytrop_var_2184, pselffac_var_2194, pselffrac_var_2195, kindself_var_2196, pforfac_var_2191, pforfrac_var_2192, kindfor_var_2193, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 23) THEN
        CALL srtm_taumol23(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, pcolh2o_var_2187, pcolmol_var_2188, klaytrop_var_2184, pselffac_var_2194, pselffrac_var_2195, kindself_var_2196, pforfac_var_2191, pforfrac_var_2192, kindfor_var_2193, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 24) THEN
        CALL srtm_taumol24(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, poneminus_var_2182, pcolh2o_var_2187, pcolmol_var_2188, pcolo2_var_2189, pcolo3_var_2190, klaytrop_var_2184, pselffac_var_2194, pselffrac_var_2195, kindself_var_2196, pforfac_var_2191, pforfrac_var_2192, kindfor_var_2193, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 25) THEN
        CALL srtm_taumol25(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, pcolh2o_var_2187, pcolmol_var_2188, pcolo3_var_2190, klaytrop_var_2184, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 26) THEN
        CALL srtm_taumol26(kidia_var_2179, kfdia_var_2180, klev_var_2181, pcolmol_var_2188, klaytrop_var_2184, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 27) THEN
        CALL srtm_taumol27(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, pcolmol_var_2188, pcolo3_var_2190, klaytrop_var_2184, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 28) THEN
        CALL srtm_taumol28(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, poneminus_var_2182, pcolmol_var_2188, pcolo2_var_2189, pcolo3_var_2190, klaytrop_var_2184, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      ELSE IF (jb_var_2205 == 29) THEN
        CALL srtm_taumol29(kidia_var_2179, kfdia_var_2180, klev_var_2181, pfac00_var_2197, pfac01_var_2198, pfac10_var_2199, pfac11_var_2200, kjp_var_2201, kjt_var_2202, kjt1_var_2203, pcolh2o_var_2187, pcolco2_var_2186, pcolmol_var_2188, klaytrop_var_2184, pselffac_var_2194, pselffrac_var_2195, kindself_var_2196, pforfac_var_2191, pforfrac_var_2192, kindfor_var_2193, zsflxzen, ztaug, ztaur, prmu0_var_2183)
      END IF
      DO jg_var_2206 = 1, igt
        DO jl_var_2208 = kidia_var_2179, kfdia_var_2180
          IF (prmu0_var_2183(jl_var_2208) > 0.0D0) THEN
            iw(jl_var_2208) = iw(jl_var_2208) + 1
            pincsol(jl_var_2208, iw(jl_var_2208)) = zsflxzen(jl_var_2208, jg_var_2206)
          END IF
        END DO
        DO jk_var_2207 = 1, klev_var_2181
          DO jl_var_2208 = kidia_var_2179, kfdia_var_2180
            IF (prmu0_var_2183(jl_var_2208) > 0.0D0) THEN
              pod_var_2204(jl_var_2208, jk_var_2207, iw(jl_var_2208)) = ztaur(jl_var_2208, jk_var_2207, jg_var_2206) + ztaug(jl_var_2208, jk_var_2207, jg_var_2206)
              pssa(jl_var_2208, jk_var_2207, iw(jl_var_2208)) = ztaur(jl_var_2208, jk_var_2207, jg_var_2206) / pod_var_2204(jl_var_2208, jk_var_2207, iw(jl_var_2208))
            END IF
          END DO
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE srtm_gas_optical_depth
SUBROUTINE srtm_taumol28(kidia_var_2209, kfdia_var_2210, klev_var_2211, p_fac00_var_2212, p_fac01_var_2213, p_fac10_var_2214, p_fac11_var_2215, k_jp_var_2216, k_jt_var_2217, k_jt1_var_2218, p_oneminus_var_2219, p_colmol_var_2220, p_colo2_var_2221, p_colo3_var_2222, k_laytrop_var_2223, p_sfluxzen_var_2224, p_taug_var_2225, p_taur_var_2226, prmu0_var_2227)
  USE yoesrta28, ONLY: absa_var_303, absb_var_304, layreffr_var_302, rayl_var_300, sfluxrefc_var_305, strrat_var_301
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2209, kfdia_var_2210
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2211
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_2212(kidia_var_2209 : kfdia_var_2210, klev_var_2211)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_2213(kidia_var_2209 : kfdia_var_2210, klev_var_2211)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_2214(kidia_var_2209 : kfdia_var_2210, klev_var_2211)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_2215(kidia_var_2209 : kfdia_var_2210, klev_var_2211)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_2216(kidia_var_2209 : kfdia_var_2210, klev_var_2211)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_2217(kidia_var_2209 : kfdia_var_2210, klev_var_2211)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_2218(kidia_var_2209 : kfdia_var_2210, klev_var_2211)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_2219(kidia_var_2209 : kfdia_var_2210)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_2220(kidia_var_2209 : kfdia_var_2210, klev_var_2211)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_2221(kidia_var_2209 : kfdia_var_2210, klev_var_2211)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_2222(kidia_var_2209 : kfdia_var_2210, klev_var_2211)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_2223(kidia_var_2209 : kfdia_var_2210)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_2224(kidia_var_2209 : kfdia_var_2210, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_2225(kidia_var_2209 : kfdia_var_2210, klev_var_2211, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_2226(kidia_var_2209 : kfdia_var_2210, klev_var_2211, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_2227(kidia_var_2209 : kfdia_var_2210)
  INTEGER(KIND = 4) :: ig_var_2228, ind0_var_2229, ind1_var_2230, js_var_2231, i_lay_var_2232, i_laysolfr_var_2233(kidia_var_2209 : kfdia_var_2210), i_nlayers_var_2234, iplon_var_2235
  INTEGER(KIND = 4) :: laytrop_min_var_2236, laytrop_max_var_2237
  REAL(KIND = 8) :: z_fs_var_2238, z_speccomb_var_2239, z_specmult_var_2240, z_specparm_var_2241, z_tauray_var_2242
  laytrop_min_var_2236 = MINVAL(k_laytrop_var_2223(kidia_var_2209 : kfdia_var_2210))
  laytrop_max_var_2237 = MAXVAL(k_laytrop_var_2223(kidia_var_2209 : kfdia_var_2210))
  i_nlayers_var_2234 = klev_var_2211
  DO iplon_var_2235 = kidia_var_2209, kfdia_var_2210
    i_laysolfr_var_2233(iplon_var_2235) = i_nlayers_var_2234
  END DO
  DO i_lay_var_2232 = 1, laytrop_min_var_2236
    DO iplon_var_2235 = kidia_var_2209, kfdia_var_2210
      z_speccomb_var_2239 = p_colo3_var_2222(iplon_var_2235, i_lay_var_2232) + strrat_var_301 * p_colo2_var_2221(iplon_var_2235, i_lay_var_2232)
      z_specparm_var_2241 = p_colo3_var_2222(iplon_var_2235, i_lay_var_2232) / z_speccomb_var_2239
      z_specparm_var_2241 = MIN(p_oneminus_var_2219(iplon_var_2235), z_specparm_var_2241)
      z_specmult_var_2240 = 8.0D0 * (z_specparm_var_2241)
      js_var_2231 = 1 + INT(z_specmult_var_2240)
      z_fs_var_2238 = z_specmult_var_2240 - AINT(z_specmult_var_2240)
      ind0_var_2229 = ((k_jp_var_2216(iplon_var_2235, i_lay_var_2232) - 1) * 5 + (k_jt_var_2217(iplon_var_2235, i_lay_var_2232) - 1)) * nspa_var_313(28) + js_var_2231
      ind1_var_2230 = (k_jp_var_2216(iplon_var_2235, i_lay_var_2232) * 5 + (k_jt1_var_2218(iplon_var_2235, i_lay_var_2232) - 1)) * nspa_var_313(28) + js_var_2231
      z_tauray_var_2242 = p_colmol_var_2220(iplon_var_2235, i_lay_var_2232) * rayl_var_300
      DO ig_var_2228 = 1, 6
        p_taug_var_2225(iplon_var_2235, i_lay_var_2232, ig_var_2228) = z_speccomb_var_2239 * ((1.0D0 - z_fs_var_2238) * (absa_var_303(ind0_var_2229, ig_var_2228) * p_fac00_var_2212(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind0_var_2229 + 9, ig_var_2228) * p_fac10_var_2214(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind1_var_2230, ig_var_2228) * p_fac01_var_2213(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind1_var_2230 + 9, ig_var_2228) * p_fac11_var_2215(iplon_var_2235, i_lay_var_2232)) + z_fs_var_2238 * (absa_var_303(ind0_var_2229 + 1, ig_var_2228) * p_fac00_var_2212(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind0_var_2229 + 10, ig_var_2228) * p_fac10_var_2214(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind1_var_2230 + 1, ig_var_2228) * p_fac01_var_2213(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind1_var_2230 + 10, ig_var_2228) * p_fac11_var_2215(iplon_var_2235, i_lay_var_2232)))
        p_taur_var_2226(iplon_var_2235, i_lay_var_2232, ig_var_2228) = z_tauray_var_2242
      END DO
    END DO
  END DO
  DO i_lay_var_2232 = laytrop_min_var_2236 + 1, laytrop_max_var_2237
    DO iplon_var_2235 = kidia_var_2209, kfdia_var_2210
      IF (i_lay_var_2232 <= k_laytrop_var_2223(iplon_var_2235)) THEN
        z_speccomb_var_2239 = p_colo3_var_2222(iplon_var_2235, i_lay_var_2232) + strrat_var_301 * p_colo2_var_2221(iplon_var_2235, i_lay_var_2232)
        z_specparm_var_2241 = p_colo3_var_2222(iplon_var_2235, i_lay_var_2232) / z_speccomb_var_2239
        z_specparm_var_2241 = MIN(p_oneminus_var_2219(iplon_var_2235), z_specparm_var_2241)
        z_specmult_var_2240 = 8.0D0 * (z_specparm_var_2241)
        js_var_2231 = 1 + INT(z_specmult_var_2240)
        z_fs_var_2238 = z_specmult_var_2240 - AINT(z_specmult_var_2240)
        ind0_var_2229 = ((k_jp_var_2216(iplon_var_2235, i_lay_var_2232) - 1) * 5 + (k_jt_var_2217(iplon_var_2235, i_lay_var_2232) - 1)) * nspa_var_313(28) + js_var_2231
        ind1_var_2230 = (k_jp_var_2216(iplon_var_2235, i_lay_var_2232) * 5 + (k_jt1_var_2218(iplon_var_2235, i_lay_var_2232) - 1)) * nspa_var_313(28) + js_var_2231
        z_tauray_var_2242 = p_colmol_var_2220(iplon_var_2235, i_lay_var_2232) * rayl_var_300
        DO ig_var_2228 = 1, 6
          p_taug_var_2225(iplon_var_2235, i_lay_var_2232, ig_var_2228) = z_speccomb_var_2239 * ((1.0D0 - z_fs_var_2238) * (absa_var_303(ind0_var_2229, ig_var_2228) * p_fac00_var_2212(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind0_var_2229 + 9, ig_var_2228) * p_fac10_var_2214(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind1_var_2230, ig_var_2228) * p_fac01_var_2213(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind1_var_2230 + 9, ig_var_2228) * p_fac11_var_2215(iplon_var_2235, i_lay_var_2232)) + z_fs_var_2238 * (absa_var_303(ind0_var_2229 + 1, ig_var_2228) * p_fac00_var_2212(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind0_var_2229 + 10, ig_var_2228) * p_fac10_var_2214(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind1_var_2230 + 1, ig_var_2228) * p_fac01_var_2213(iplon_var_2235, i_lay_var_2232) + absa_var_303(ind1_var_2230 + 10, ig_var_2228) * p_fac11_var_2215(iplon_var_2235, i_lay_var_2232)))
          p_taur_var_2226(iplon_var_2235, i_lay_var_2232, ig_var_2228) = z_tauray_var_2242
        END DO
      ELSE
        IF (k_jp_var_2216(iplon_var_2235, i_lay_var_2232 - 1) < layreffr_var_302 .AND. k_jp_var_2216(iplon_var_2235, i_lay_var_2232) >= layreffr_var_302) i_laysolfr_var_2233(iplon_var_2235) = i_lay_var_2232
        z_speccomb_var_2239 = p_colo3_var_2222(iplon_var_2235, i_lay_var_2232) + strrat_var_301 * p_colo2_var_2221(iplon_var_2235, i_lay_var_2232)
        z_specparm_var_2241 = p_colo3_var_2222(iplon_var_2235, i_lay_var_2232) / z_speccomb_var_2239
        z_specparm_var_2241 = MIN(p_oneminus_var_2219(iplon_var_2235), z_specparm_var_2241)
        z_specmult_var_2240 = 4.0D0 * (z_specparm_var_2241)
        js_var_2231 = 1 + INT(z_specmult_var_2240)
        z_fs_var_2238 = z_specmult_var_2240 - AINT(z_specmult_var_2240)
        ind0_var_2229 = ((k_jp_var_2216(iplon_var_2235, i_lay_var_2232) - 13) * 5 + (k_jt_var_2217(iplon_var_2235, i_lay_var_2232) - 1)) * nspb_var_314(28) + js_var_2231
        ind1_var_2230 = ((k_jp_var_2216(iplon_var_2235, i_lay_var_2232) - 12) * 5 + (k_jt1_var_2218(iplon_var_2235, i_lay_var_2232) - 1)) * nspb_var_314(28) + js_var_2231
        z_tauray_var_2242 = p_colmol_var_2220(iplon_var_2235, i_lay_var_2232) * rayl_var_300
        DO ig_var_2228 = 1, 6
          p_taug_var_2225(iplon_var_2235, i_lay_var_2232, ig_var_2228) = z_speccomb_var_2239 * ((1.0D0 - z_fs_var_2238) * (absb_var_304(ind0_var_2229, ig_var_2228) * p_fac00_var_2212(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind0_var_2229 + 5, ig_var_2228) * p_fac10_var_2214(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind1_var_2230, ig_var_2228) * p_fac01_var_2213(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind1_var_2230 + 5, ig_var_2228) * p_fac11_var_2215(iplon_var_2235, i_lay_var_2232)) + z_fs_var_2238 * (absb_var_304(ind0_var_2229 + 1, ig_var_2228) * p_fac00_var_2212(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind0_var_2229 + 6, ig_var_2228) * p_fac10_var_2214(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind1_var_2230 + 1, ig_var_2228) * p_fac01_var_2213(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind1_var_2230 + 6, ig_var_2228) * p_fac11_var_2215(iplon_var_2235, i_lay_var_2232)))
          IF (i_lay_var_2232 == i_laysolfr_var_2233(iplon_var_2235)) p_sfluxzen_var_2224(iplon_var_2235, ig_var_2228) = sfluxrefc_var_305(ig_var_2228, js_var_2231) + z_fs_var_2238 * (sfluxrefc_var_305(ig_var_2228, js_var_2231 + 1) - sfluxrefc_var_305(ig_var_2228, js_var_2231))
          p_taur_var_2226(iplon_var_2235, i_lay_var_2232, ig_var_2228) = z_tauray_var_2242
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_2232 = laytrop_max_var_2237 + 1, i_nlayers_var_2234
    DO iplon_var_2235 = kidia_var_2209, kfdia_var_2210
      IF (k_jp_var_2216(iplon_var_2235, i_lay_var_2232 - 1) < layreffr_var_302 .AND. k_jp_var_2216(iplon_var_2235, i_lay_var_2232) >= layreffr_var_302) i_laysolfr_var_2233(iplon_var_2235) = i_lay_var_2232
      z_speccomb_var_2239 = p_colo3_var_2222(iplon_var_2235, i_lay_var_2232) + strrat_var_301 * p_colo2_var_2221(iplon_var_2235, i_lay_var_2232)
      z_specparm_var_2241 = p_colo3_var_2222(iplon_var_2235, i_lay_var_2232) / z_speccomb_var_2239
      z_specparm_var_2241 = MIN(p_oneminus_var_2219(iplon_var_2235), z_specparm_var_2241)
      z_specmult_var_2240 = 4.0D0 * (z_specparm_var_2241)
      js_var_2231 = 1 + INT(z_specmult_var_2240)
      z_fs_var_2238 = z_specmult_var_2240 - AINT(z_specmult_var_2240)
      ind0_var_2229 = ((k_jp_var_2216(iplon_var_2235, i_lay_var_2232) - 13) * 5 + (k_jt_var_2217(iplon_var_2235, i_lay_var_2232) - 1)) * nspb_var_314(28) + js_var_2231
      ind1_var_2230 = ((k_jp_var_2216(iplon_var_2235, i_lay_var_2232) - 12) * 5 + (k_jt1_var_2218(iplon_var_2235, i_lay_var_2232) - 1)) * nspb_var_314(28) + js_var_2231
      z_tauray_var_2242 = p_colmol_var_2220(iplon_var_2235, i_lay_var_2232) * rayl_var_300
      DO ig_var_2228 = 1, 6
        p_taug_var_2225(iplon_var_2235, i_lay_var_2232, ig_var_2228) = z_speccomb_var_2239 * ((1.0D0 - z_fs_var_2238) * (absb_var_304(ind0_var_2229, ig_var_2228) * p_fac00_var_2212(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind0_var_2229 + 5, ig_var_2228) * p_fac10_var_2214(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind1_var_2230, ig_var_2228) * p_fac01_var_2213(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind1_var_2230 + 5, ig_var_2228) * p_fac11_var_2215(iplon_var_2235, i_lay_var_2232)) + z_fs_var_2238 * (absb_var_304(ind0_var_2229 + 1, ig_var_2228) * p_fac00_var_2212(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind0_var_2229 + 6, ig_var_2228) * p_fac10_var_2214(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind1_var_2230 + 1, ig_var_2228) * p_fac01_var_2213(iplon_var_2235, i_lay_var_2232) + absb_var_304(ind1_var_2230 + 6, ig_var_2228) * p_fac11_var_2215(iplon_var_2235, i_lay_var_2232)))
        IF (i_lay_var_2232 == i_laysolfr_var_2233(iplon_var_2235)) p_sfluxzen_var_2224(iplon_var_2235, ig_var_2228) = sfluxrefc_var_305(ig_var_2228, js_var_2231) + z_fs_var_2238 * (sfluxrefc_var_305(ig_var_2228, js_var_2231 + 1) - sfluxrefc_var_305(ig_var_2228, js_var_2231))
        p_taur_var_2226(iplon_var_2235, i_lay_var_2232, ig_var_2228) = z_tauray_var_2242
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol28
SUBROUTINE srtm_taumol29(kidia_var_2243, kfdia_var_2244, klev_var_2245, p_fac00_var_2246, p_fac01_var_2247, p_fac10_var_2248, p_fac11_var_2249, k_jp_var_2250, k_jt_var_2251, k_jt1_var_2252, p_colh2o_var_2253, p_colco2_var_2254, p_colmol_var_2255, k_laytrop_var_2256, p_selffac_var_2257, p_selffrac_var_2258, k_indself_var_2259, p_forfac_var_2260, p_forfrac_var_2261, k_indfor_var_2262, p_sfluxzen_var_2265, p_taug_var_2266, p_taur_var_2267, prmu0_var_2268)
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  USE yoesrta29, ONLY: absa_var_308, absb_var_309, absco2c, absh2oc, forrefc_var_311, layreffr_var_307, rayl_var_306, selfrefc_var_310, sfluxrefc_var_312
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2243, kfdia_var_2244
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2245
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_2246(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_2247(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_2248(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_2249(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_2250(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_2251(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_2252(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_2253(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_2254(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_2255(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_2256(kidia_var_2243 : kfdia_var_2244)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_2257(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_2258(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_2259(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_2260(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_2261(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_2262(kidia_var_2243 : kfdia_var_2244, klev_var_2245)
  INTEGER(KIND = 4) :: laytrop_min_var_2263, laytrop_max_var_2264
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_2265(kidia_var_2243 : kfdia_var_2244, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_2266(kidia_var_2243 : kfdia_var_2244, klev_var_2245, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_2267(kidia_var_2243 : kfdia_var_2244, klev_var_2245, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_2268(kidia_var_2243 : kfdia_var_2244)
  INTEGER(KIND = 4) :: ig_var_2269, ind0_var_2270, ind1_var_2271, inds_var_2272, indf_var_2273, i_lay_var_2274, i_laysolfr_var_2275(kidia_var_2243 : kfdia_var_2244), i_nlayers_var_2276, iplon_var_2277
  REAL(KIND = 8) :: z_tauray_var_2278
  laytrop_min_var_2263 = MINVAL(k_laytrop_var_2256(kidia_var_2243 : kfdia_var_2244))
  laytrop_max_var_2264 = MAXVAL(k_laytrop_var_2256(kidia_var_2243 : kfdia_var_2244))
  i_nlayers_var_2276 = klev_var_2245
  DO iplon_var_2277 = kidia_var_2243, kfdia_var_2244
    i_laysolfr_var_2275(iplon_var_2277) = i_nlayers_var_2276
  END DO
  DO i_lay_var_2274 = 1, laytrop_min_var_2263
    DO iplon_var_2277 = kidia_var_2243, kfdia_var_2244
      ind0_var_2270 = ((k_jp_var_2250(iplon_var_2277, i_lay_var_2274) - 1) * 5 + (k_jt_var_2251(iplon_var_2277, i_lay_var_2274) - 1)) * nspa_var_313(29) + 1
      ind1_var_2271 = (k_jp_var_2250(iplon_var_2277, i_lay_var_2274) * 5 + (k_jt1_var_2252(iplon_var_2277, i_lay_var_2274) - 1)) * nspa_var_313(29) + 1
      inds_var_2272 = k_indself_var_2259(iplon_var_2277, i_lay_var_2274)
      indf_var_2273 = k_indfor_var_2262(iplon_var_2277, i_lay_var_2274)
      z_tauray_var_2278 = p_colmol_var_2255(iplon_var_2277, i_lay_var_2274) * rayl_var_306
      DO ig_var_2269 = 1, 12
        p_taug_var_2266(iplon_var_2277, i_lay_var_2274, ig_var_2269) = p_colh2o_var_2253(iplon_var_2277, i_lay_var_2274) * ((p_fac00_var_2246(iplon_var_2277, i_lay_var_2274) * absa_var_308(ind0_var_2270, ig_var_2269) + p_fac10_var_2248(iplon_var_2277, i_lay_var_2274) * absa_var_308(ind0_var_2270 + 1, ig_var_2269) + p_fac01_var_2247(iplon_var_2277, i_lay_var_2274) * absa_var_308(ind1_var_2271, ig_var_2269) + p_fac11_var_2249(iplon_var_2277, i_lay_var_2274) * absa_var_308(ind1_var_2271 + 1, ig_var_2269)) + p_selffac_var_2257(iplon_var_2277, i_lay_var_2274) * (selfrefc_var_310(inds_var_2272, ig_var_2269) + p_selffrac_var_2258(iplon_var_2277, i_lay_var_2274) * (selfrefc_var_310(inds_var_2272 + 1, ig_var_2269) - selfrefc_var_310(inds_var_2272, ig_var_2269))) + p_forfac_var_2260(iplon_var_2277, i_lay_var_2274) * (forrefc_var_311(indf_var_2273, ig_var_2269) + p_forfrac_var_2261(iplon_var_2277, i_lay_var_2274) * (forrefc_var_311(indf_var_2273 + 1, ig_var_2269) - forrefc_var_311(indf_var_2273, ig_var_2269)))) + p_colco2_var_2254(iplon_var_2277, i_lay_var_2274) * absco2c(ig_var_2269)
        p_taur_var_2267(iplon_var_2277, i_lay_var_2274, ig_var_2269) = z_tauray_var_2278
      END DO
    END DO
  END DO
  DO i_lay_var_2274 = laytrop_min_var_2263 + 1, laytrop_max_var_2264
    DO iplon_var_2277 = kidia_var_2243, kfdia_var_2244
      IF (i_lay_var_2274 <= k_laytrop_var_2256(iplon_var_2277)) THEN
        ind0_var_2270 = ((k_jp_var_2250(iplon_var_2277, i_lay_var_2274) - 1) * 5 + (k_jt_var_2251(iplon_var_2277, i_lay_var_2274) - 1)) * nspa_var_313(29) + 1
        ind1_var_2271 = (k_jp_var_2250(iplon_var_2277, i_lay_var_2274) * 5 + (k_jt1_var_2252(iplon_var_2277, i_lay_var_2274) - 1)) * nspa_var_313(29) + 1
        inds_var_2272 = k_indself_var_2259(iplon_var_2277, i_lay_var_2274)
        indf_var_2273 = k_indfor_var_2262(iplon_var_2277, i_lay_var_2274)
        z_tauray_var_2278 = p_colmol_var_2255(iplon_var_2277, i_lay_var_2274) * rayl_var_306
        DO ig_var_2269 = 1, 12
          p_taug_var_2266(iplon_var_2277, i_lay_var_2274, ig_var_2269) = p_colh2o_var_2253(iplon_var_2277, i_lay_var_2274) * ((p_fac00_var_2246(iplon_var_2277, i_lay_var_2274) * absa_var_308(ind0_var_2270, ig_var_2269) + p_fac10_var_2248(iplon_var_2277, i_lay_var_2274) * absa_var_308(ind0_var_2270 + 1, ig_var_2269) + p_fac01_var_2247(iplon_var_2277, i_lay_var_2274) * absa_var_308(ind1_var_2271, ig_var_2269) + p_fac11_var_2249(iplon_var_2277, i_lay_var_2274) * absa_var_308(ind1_var_2271 + 1, ig_var_2269)) + p_selffac_var_2257(iplon_var_2277, i_lay_var_2274) * (selfrefc_var_310(inds_var_2272, ig_var_2269) + p_selffrac_var_2258(iplon_var_2277, i_lay_var_2274) * (selfrefc_var_310(inds_var_2272 + 1, ig_var_2269) - selfrefc_var_310(inds_var_2272, ig_var_2269))) + p_forfac_var_2260(iplon_var_2277, i_lay_var_2274) * (forrefc_var_311(indf_var_2273, ig_var_2269) + p_forfrac_var_2261(iplon_var_2277, i_lay_var_2274) * (forrefc_var_311(indf_var_2273 + 1, ig_var_2269) - forrefc_var_311(indf_var_2273, ig_var_2269)))) + p_colco2_var_2254(iplon_var_2277, i_lay_var_2274) * absco2c(ig_var_2269)
          p_taur_var_2267(iplon_var_2277, i_lay_var_2274, ig_var_2269) = z_tauray_var_2278
        END DO
      ELSE
        IF (k_jp_var_2250(iplon_var_2277, i_lay_var_2274 - 1) < layreffr_var_307 .AND. k_jp_var_2250(iplon_var_2277, i_lay_var_2274) >= layreffr_var_307) i_laysolfr_var_2275(iplon_var_2277) = i_lay_var_2274
        ind0_var_2270 = ((k_jp_var_2250(iplon_var_2277, i_lay_var_2274) - 13) * 5 + (k_jt_var_2251(iplon_var_2277, i_lay_var_2274) - 1)) * nspb_var_314(29) + 1
        ind1_var_2271 = ((k_jp_var_2250(iplon_var_2277, i_lay_var_2274) - 12) * 5 + (k_jt1_var_2252(iplon_var_2277, i_lay_var_2274) - 1)) * nspb_var_314(29) + 1
        z_tauray_var_2278 = p_colmol_var_2255(iplon_var_2277, i_lay_var_2274) * rayl_var_306
        DO ig_var_2269 = 1, 12
          p_taug_var_2266(iplon_var_2277, i_lay_var_2274, ig_var_2269) = p_colco2_var_2254(iplon_var_2277, i_lay_var_2274) * (p_fac00_var_2246(iplon_var_2277, i_lay_var_2274) * absb_var_309(ind0_var_2270, ig_var_2269) + p_fac10_var_2248(iplon_var_2277, i_lay_var_2274) * absb_var_309(ind0_var_2270 + 1, ig_var_2269) + p_fac01_var_2247(iplon_var_2277, i_lay_var_2274) * absb_var_309(ind1_var_2271, ig_var_2269) + p_fac11_var_2249(iplon_var_2277, i_lay_var_2274) * absb_var_309(ind1_var_2271 + 1, ig_var_2269)) + p_colh2o_var_2253(iplon_var_2277, i_lay_var_2274) * absh2oc(ig_var_2269)
          IF (i_lay_var_2274 == i_laysolfr_var_2275(iplon_var_2277)) p_sfluxzen_var_2265(iplon_var_2277, ig_var_2269) = sfluxrefc_var_312(ig_var_2269)
          p_taur_var_2267(iplon_var_2277, i_lay_var_2274, ig_var_2269) = z_tauray_var_2278
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_2274 = laytrop_max_var_2264 + 1, i_nlayers_var_2276
    DO iplon_var_2277 = kidia_var_2243, kfdia_var_2244
      IF (k_jp_var_2250(iplon_var_2277, i_lay_var_2274 - 1) < layreffr_var_307 .AND. k_jp_var_2250(iplon_var_2277, i_lay_var_2274) >= layreffr_var_307) i_laysolfr_var_2275(iplon_var_2277) = i_lay_var_2274
      ind0_var_2270 = ((k_jp_var_2250(iplon_var_2277, i_lay_var_2274) - 13) * 5 + (k_jt_var_2251(iplon_var_2277, i_lay_var_2274) - 1)) * nspb_var_314(29) + 1
      ind1_var_2271 = ((k_jp_var_2250(iplon_var_2277, i_lay_var_2274) - 12) * 5 + (k_jt1_var_2252(iplon_var_2277, i_lay_var_2274) - 1)) * nspb_var_314(29) + 1
      z_tauray_var_2278 = p_colmol_var_2255(iplon_var_2277, i_lay_var_2274) * rayl_var_306
      DO ig_var_2269 = 1, 12
        p_taug_var_2266(iplon_var_2277, i_lay_var_2274, ig_var_2269) = p_colco2_var_2254(iplon_var_2277, i_lay_var_2274) * (p_fac00_var_2246(iplon_var_2277, i_lay_var_2274) * absb_var_309(ind0_var_2270, ig_var_2269) + p_fac10_var_2248(iplon_var_2277, i_lay_var_2274) * absb_var_309(ind0_var_2270 + 1, ig_var_2269) + p_fac01_var_2247(iplon_var_2277, i_lay_var_2274) * absb_var_309(ind1_var_2271, ig_var_2269) + p_fac11_var_2249(iplon_var_2277, i_lay_var_2274) * absb_var_309(ind1_var_2271 + 1, ig_var_2269)) + p_colh2o_var_2253(iplon_var_2277, i_lay_var_2274) * absh2oc(ig_var_2269)
        IF (i_lay_var_2274 == i_laysolfr_var_2275(iplon_var_2277)) p_sfluxzen_var_2265(iplon_var_2277, ig_var_2269) = sfluxrefc_var_312(ig_var_2269)
        p_taur_var_2267(iplon_var_2277, i_lay_var_2274, ig_var_2269) = z_tauray_var_2278
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol29
SUBROUTINE srtm_taumol17(kidia_var_2279, kfdia_var_2280, klev_var_2281, p_fac00_var_2282, p_fac01_var_2283, p_fac10_var_2284, p_fac11_var_2285, k_jp_var_2286, k_jt_var_2287, k_jt1_var_2288, p_oneminus_var_2289, p_colh2o_var_2290, p_colco2_var_2291, p_colmol_var_2292, k_laytrop_var_2293, p_selffac_var_2294, p_selffrac_var_2295, k_indself_var_2296, p_forfac_var_2297, p_forfrac_var_2298, k_indfor_var_2299, p_sfluxzen_var_2300, p_taug_var_2301, p_taur_var_2302, prmu0_var_2303)
  USE yoesrta17, ONLY: absa_var_228, absb_var_229, forrefc_var_231, layreffr_var_227, rayl_var_225, selfrefc_var_230, sfluxrefc_var_232, strrat_var_226
  USE yoesrtwn, ONLY: nspa_var_313, nspb_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2279, kfdia_var_2280
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2281
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_2282(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_2283(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_2284(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_2285(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_2286(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_2287(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_2288(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_2289(kidia_var_2279 : kfdia_var_2280)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_2290(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_2291(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_2292(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_2293(kidia_var_2279 : kfdia_var_2280)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_2294(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_2295(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_2296(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_2297(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_2298(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_2299(kidia_var_2279 : kfdia_var_2280, klev_var_2281)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_2300(kidia_var_2279 : kfdia_var_2280, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_2301(kidia_var_2279 : kfdia_var_2280, klev_var_2281, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_2302(kidia_var_2279 : kfdia_var_2280, klev_var_2281, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_2303(kidia_var_2279 : kfdia_var_2280)
  INTEGER(KIND = 4) :: ig_var_2304, ind0_var_2305, ind1_var_2306, inds_var_2307, indf_var_2308, js_var_2309, i_lay_var_2310, i_laysolfr_var_2311(kidia_var_2279 : kfdia_var_2280), i_nlayers_var_2312, iplon_var_2313
  INTEGER(KIND = 4) :: laytrop_min_var_2314, laytrop_max_var_2315
  REAL(KIND = 8) :: z_fs_var_2316, z_speccomb_var_2317, z_specmult_var_2318, z_specparm_var_2319, z_tauray_var_2320
  laytrop_min_var_2314 = MINVAL(k_laytrop_var_2293(kidia_var_2279 : kfdia_var_2280))
  laytrop_max_var_2315 = MAXVAL(k_laytrop_var_2293(kidia_var_2279 : kfdia_var_2280))
  i_nlayers_var_2312 = klev_var_2281
  DO iplon_var_2313 = kidia_var_2279, kfdia_var_2280
    i_laysolfr_var_2311(iplon_var_2313) = i_nlayers_var_2312
  END DO
  DO i_lay_var_2310 = 1, laytrop_min_var_2314
    DO iplon_var_2313 = kidia_var_2279, kfdia_var_2280
      z_speccomb_var_2317 = p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) + strrat_var_226 * p_colco2_var_2291(iplon_var_2313, i_lay_var_2310)
      z_specparm_var_2319 = p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) / z_speccomb_var_2317
      z_specparm_var_2319 = MIN(p_oneminus_var_2289(iplon_var_2313), z_specparm_var_2319)
      z_specmult_var_2318 = 8.0D0 * (z_specparm_var_2319)
      js_var_2309 = 1 + INT(z_specmult_var_2318)
      z_fs_var_2316 = z_specmult_var_2318 - AINT(z_specmult_var_2318)
      ind0_var_2305 = ((k_jp_var_2286(iplon_var_2313, i_lay_var_2310) - 1) * 5 + (k_jt_var_2287(iplon_var_2313, i_lay_var_2310) - 1)) * nspa_var_313(17) + js_var_2309
      ind1_var_2306 = (k_jp_var_2286(iplon_var_2313, i_lay_var_2310) * 5 + (k_jt1_var_2288(iplon_var_2313, i_lay_var_2310) - 1)) * nspa_var_313(17) + js_var_2309
      inds_var_2307 = k_indself_var_2296(iplon_var_2313, i_lay_var_2310)
      indf_var_2308 = k_indfor_var_2299(iplon_var_2313, i_lay_var_2310)
      z_tauray_var_2320 = p_colmol_var_2292(iplon_var_2313, i_lay_var_2310) * rayl_var_225
      DO ig_var_2304 = 1, 12
        p_taug_var_2301(iplon_var_2313, i_lay_var_2310, ig_var_2304) = z_speccomb_var_2317 * ((1.0D0 - z_fs_var_2316) * (absa_var_228(ind0_var_2305, ig_var_2304) * p_fac00_var_2282(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind0_var_2305 + 9, ig_var_2304) * p_fac10_var_2284(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind1_var_2306, ig_var_2304) * p_fac01_var_2283(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind1_var_2306 + 9, ig_var_2304) * p_fac11_var_2285(iplon_var_2313, i_lay_var_2310)) + z_fs_var_2316 * (absa_var_228(ind0_var_2305 + 1, ig_var_2304) * p_fac00_var_2282(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind0_var_2305 + 10, ig_var_2304) * p_fac10_var_2284(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind1_var_2306 + 1, ig_var_2304) * p_fac01_var_2283(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind1_var_2306 + 10, ig_var_2304) * p_fac11_var_2285(iplon_var_2313, i_lay_var_2310))) + p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) * (p_selffac_var_2294(iplon_var_2313, i_lay_var_2310) * (selfrefc_var_230(inds_var_2307, ig_var_2304) + p_selffrac_var_2295(iplon_var_2313, i_lay_var_2310) * (selfrefc_var_230(inds_var_2307 + 1, ig_var_2304) - selfrefc_var_230(inds_var_2307, ig_var_2304))) + p_forfac_var_2297(iplon_var_2313, i_lay_var_2310) * (forrefc_var_231(indf_var_2308, ig_var_2304) + p_forfrac_var_2298(iplon_var_2313, i_lay_var_2310) * (forrefc_var_231(indf_var_2308 + 1, ig_var_2304) - forrefc_var_231(indf_var_2308, ig_var_2304))))
        p_taur_var_2302(iplon_var_2313, i_lay_var_2310, ig_var_2304) = z_tauray_var_2320
      END DO
    END DO
  END DO
  DO i_lay_var_2310 = laytrop_min_var_2314 + 1, laytrop_max_var_2315
    DO iplon_var_2313 = kidia_var_2279, kfdia_var_2280
      IF (i_lay_var_2310 <= k_laytrop_var_2293(iplon_var_2313)) THEN
        z_speccomb_var_2317 = p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) + strrat_var_226 * p_colco2_var_2291(iplon_var_2313, i_lay_var_2310)
        z_specparm_var_2319 = p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) / z_speccomb_var_2317
        z_specparm_var_2319 = MIN(p_oneminus_var_2289(iplon_var_2313), z_specparm_var_2319)
        z_specmult_var_2318 = 8.0D0 * (z_specparm_var_2319)
        js_var_2309 = 1 + INT(z_specmult_var_2318)
        z_fs_var_2316 = z_specmult_var_2318 - AINT(z_specmult_var_2318)
        ind0_var_2305 = ((k_jp_var_2286(iplon_var_2313, i_lay_var_2310) - 1) * 5 + (k_jt_var_2287(iplon_var_2313, i_lay_var_2310) - 1)) * nspa_var_313(17) + js_var_2309
        ind1_var_2306 = (k_jp_var_2286(iplon_var_2313, i_lay_var_2310) * 5 + (k_jt1_var_2288(iplon_var_2313, i_lay_var_2310) - 1)) * nspa_var_313(17) + js_var_2309
        inds_var_2307 = k_indself_var_2296(iplon_var_2313, i_lay_var_2310)
        indf_var_2308 = k_indfor_var_2299(iplon_var_2313, i_lay_var_2310)
        z_tauray_var_2320 = p_colmol_var_2292(iplon_var_2313, i_lay_var_2310) * rayl_var_225
        DO ig_var_2304 = 1, 12
          p_taug_var_2301(iplon_var_2313, i_lay_var_2310, ig_var_2304) = z_speccomb_var_2317 * ((1.0D0 - z_fs_var_2316) * (absa_var_228(ind0_var_2305, ig_var_2304) * p_fac00_var_2282(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind0_var_2305 + 9, ig_var_2304) * p_fac10_var_2284(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind1_var_2306, ig_var_2304) * p_fac01_var_2283(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind1_var_2306 + 9, ig_var_2304) * p_fac11_var_2285(iplon_var_2313, i_lay_var_2310)) + z_fs_var_2316 * (absa_var_228(ind0_var_2305 + 1, ig_var_2304) * p_fac00_var_2282(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind0_var_2305 + 10, ig_var_2304) * p_fac10_var_2284(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind1_var_2306 + 1, ig_var_2304) * p_fac01_var_2283(iplon_var_2313, i_lay_var_2310) + absa_var_228(ind1_var_2306 + 10, ig_var_2304) * p_fac11_var_2285(iplon_var_2313, i_lay_var_2310))) + p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) * (p_selffac_var_2294(iplon_var_2313, i_lay_var_2310) * (selfrefc_var_230(inds_var_2307, ig_var_2304) + p_selffrac_var_2295(iplon_var_2313, i_lay_var_2310) * (selfrefc_var_230(inds_var_2307 + 1, ig_var_2304) - selfrefc_var_230(inds_var_2307, ig_var_2304))) + p_forfac_var_2297(iplon_var_2313, i_lay_var_2310) * (forrefc_var_231(indf_var_2308, ig_var_2304) + p_forfrac_var_2298(iplon_var_2313, i_lay_var_2310) * (forrefc_var_231(indf_var_2308 + 1, ig_var_2304) - forrefc_var_231(indf_var_2308, ig_var_2304))))
          p_taur_var_2302(iplon_var_2313, i_lay_var_2310, ig_var_2304) = z_tauray_var_2320
        END DO
      ELSE
        IF (k_jp_var_2286(iplon_var_2313, i_lay_var_2310 - 1) < layreffr_var_227 .AND. k_jp_var_2286(iplon_var_2313, i_lay_var_2310) >= layreffr_var_227) i_laysolfr_var_2311(iplon_var_2313) = i_lay_var_2310
        z_speccomb_var_2317 = p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) + strrat_var_226 * p_colco2_var_2291(iplon_var_2313, i_lay_var_2310)
        z_specparm_var_2319 = p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) / z_speccomb_var_2317
        z_specparm_var_2319 = MIN(p_oneminus_var_2289(iplon_var_2313), z_specparm_var_2319)
        z_specmult_var_2318 = 4.0D0 * (z_specparm_var_2319)
        js_var_2309 = 1 + INT(z_specmult_var_2318)
        z_fs_var_2316 = z_specmult_var_2318 - AINT(z_specmult_var_2318)
        ind0_var_2305 = ((k_jp_var_2286(iplon_var_2313, i_lay_var_2310) - 13) * 5 + (k_jt_var_2287(iplon_var_2313, i_lay_var_2310) - 1)) * nspb_var_314(17) + js_var_2309
        ind1_var_2306 = ((k_jp_var_2286(iplon_var_2313, i_lay_var_2310) - 12) * 5 + (k_jt1_var_2288(iplon_var_2313, i_lay_var_2310) - 1)) * nspb_var_314(17) + js_var_2309
        indf_var_2308 = k_indfor_var_2299(iplon_var_2313, i_lay_var_2310)
        z_tauray_var_2320 = p_colmol_var_2292(iplon_var_2313, i_lay_var_2310) * rayl_var_225
        DO ig_var_2304 = 1, 12
          p_taug_var_2301(iplon_var_2313, i_lay_var_2310, ig_var_2304) = z_speccomb_var_2317 * ((1.0D0 - z_fs_var_2316) * (absb_var_229(ind0_var_2305, ig_var_2304) * p_fac00_var_2282(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind0_var_2305 + 5, ig_var_2304) * p_fac10_var_2284(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind1_var_2306, ig_var_2304) * p_fac01_var_2283(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind1_var_2306 + 5, ig_var_2304) * p_fac11_var_2285(iplon_var_2313, i_lay_var_2310)) + z_fs_var_2316 * (absb_var_229(ind0_var_2305 + 1, ig_var_2304) * p_fac00_var_2282(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind0_var_2305 + 6, ig_var_2304) * p_fac10_var_2284(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind1_var_2306 + 1, ig_var_2304) * p_fac01_var_2283(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind1_var_2306 + 6, ig_var_2304) * p_fac11_var_2285(iplon_var_2313, i_lay_var_2310))) + p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) * p_forfac_var_2297(iplon_var_2313, i_lay_var_2310) * (forrefc_var_231(indf_var_2308, ig_var_2304) + p_forfrac_var_2298(iplon_var_2313, i_lay_var_2310) * (forrefc_var_231(indf_var_2308 + 1, ig_var_2304) - forrefc_var_231(indf_var_2308, ig_var_2304)))
          IF (i_lay_var_2310 == i_laysolfr_var_2311(iplon_var_2313)) p_sfluxzen_var_2300(iplon_var_2313, ig_var_2304) = sfluxrefc_var_232(ig_var_2304, js_var_2309) + z_fs_var_2316 * (sfluxrefc_var_232(ig_var_2304, js_var_2309 + 1) - sfluxrefc_var_232(ig_var_2304, js_var_2309))
          p_taur_var_2302(iplon_var_2313, i_lay_var_2310, ig_var_2304) = z_tauray_var_2320
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_2310 = laytrop_max_var_2315 + 1, i_nlayers_var_2312
    DO iplon_var_2313 = kidia_var_2279, kfdia_var_2280
      IF (k_jp_var_2286(iplon_var_2313, i_lay_var_2310 - 1) < layreffr_var_227 .AND. k_jp_var_2286(iplon_var_2313, i_lay_var_2310) >= layreffr_var_227) i_laysolfr_var_2311(iplon_var_2313) = i_lay_var_2310
      z_speccomb_var_2317 = p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) + strrat_var_226 * p_colco2_var_2291(iplon_var_2313, i_lay_var_2310)
      z_specparm_var_2319 = p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) / z_speccomb_var_2317
      z_specparm_var_2319 = MIN(p_oneminus_var_2289(iplon_var_2313), z_specparm_var_2319)
      z_specmult_var_2318 = 4.0D0 * (z_specparm_var_2319)
      js_var_2309 = 1 + INT(z_specmult_var_2318)
      z_fs_var_2316 = z_specmult_var_2318 - AINT(z_specmult_var_2318)
      ind0_var_2305 = ((k_jp_var_2286(iplon_var_2313, i_lay_var_2310) - 13) * 5 + (k_jt_var_2287(iplon_var_2313, i_lay_var_2310) - 1)) * nspb_var_314(17) + js_var_2309
      ind1_var_2306 = ((k_jp_var_2286(iplon_var_2313, i_lay_var_2310) - 12) * 5 + (k_jt1_var_2288(iplon_var_2313, i_lay_var_2310) - 1)) * nspb_var_314(17) + js_var_2309
      indf_var_2308 = k_indfor_var_2299(iplon_var_2313, i_lay_var_2310)
      z_tauray_var_2320 = p_colmol_var_2292(iplon_var_2313, i_lay_var_2310) * rayl_var_225
      DO ig_var_2304 = 1, 12
        p_taug_var_2301(iplon_var_2313, i_lay_var_2310, ig_var_2304) = z_speccomb_var_2317 * ((1.0D0 - z_fs_var_2316) * (absb_var_229(ind0_var_2305, ig_var_2304) * p_fac00_var_2282(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind0_var_2305 + 5, ig_var_2304) * p_fac10_var_2284(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind1_var_2306, ig_var_2304) * p_fac01_var_2283(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind1_var_2306 + 5, ig_var_2304) * p_fac11_var_2285(iplon_var_2313, i_lay_var_2310)) + z_fs_var_2316 * (absb_var_229(ind0_var_2305 + 1, ig_var_2304) * p_fac00_var_2282(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind0_var_2305 + 6, ig_var_2304) * p_fac10_var_2284(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind1_var_2306 + 1, ig_var_2304) * p_fac01_var_2283(iplon_var_2313, i_lay_var_2310) + absb_var_229(ind1_var_2306 + 6, ig_var_2304) * p_fac11_var_2285(iplon_var_2313, i_lay_var_2310))) + p_colh2o_var_2290(iplon_var_2313, i_lay_var_2310) * p_forfac_var_2297(iplon_var_2313, i_lay_var_2310) * (forrefc_var_231(indf_var_2308, ig_var_2304) + p_forfrac_var_2298(iplon_var_2313, i_lay_var_2310) * (forrefc_var_231(indf_var_2308 + 1, ig_var_2304) - forrefc_var_231(indf_var_2308, ig_var_2304)))
        IF (i_lay_var_2310 == i_laysolfr_var_2311(iplon_var_2313)) p_sfluxzen_var_2300(iplon_var_2313, ig_var_2304) = sfluxrefc_var_232(ig_var_2304, js_var_2309) + z_fs_var_2316 * (sfluxrefc_var_232(ig_var_2304, js_var_2309 + 1) - sfluxrefc_var_232(ig_var_2304, js_var_2309))
        p_taur_var_2302(iplon_var_2313, i_lay_var_2310, ig_var_2304) = z_tauray_var_2320
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol17
SUBROUTINE rrtm_setcoef_140gp(kidia_var_2321, kfdia_var_2322, klev_var_2323, p_coldry, p_wbroad, p_wkl, p_fac00_var_2324, p_fac01_var_2325, p_fac10_var_2326, p_fac11_var_2327, p_forfac_var_2328, p_forfrac_var_2329, k_indfor_var_2346, k_jp_var_2330, k_jt_var_2331, k_jt1_var_2332, p_colh2o_var_2333, p_colco2_var_2334, p_colo3_var_2335, p_coln2o, p_colch4_var_2336, p_colo2_var_2337, p_co2mult_var_2338, p_colbrd, k_laytrop_var_2339, k_layswtch_var_2340, k_laylow_var_2341, pavel_var_2342, p_tavel, p_selffac_var_2343, p_selffrac_var_2344, k_indself_var_2345, k_indminor, p_scaleminor, p_scaleminorn2, p_minorfrac, prat_h2oco2_var_2347, prat_h2oco2_1_var_2348, prat_h2oo3_var_2349, prat_h2oo3_1_var_2350, prat_h2on2o_var_2351, prat_h2on2o_1_var_2352, prat_h2och4_var_2353, prat_h2och4_1_var_2354, prat_n2oco2_var_2355, prat_n2oco2_1_var_2356, prat_o3co2_var_2357, prat_o3co2_1_var_2358)
  USE yoerrtrf, ONLY: chi_mls, preflog_var_214, tref_var_215
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2321
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2322
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2323
  REAL(KIND = 8), INTENT(IN) :: p_coldry(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(IN) :: p_wbroad(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_colbrd(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(IN) :: p_wkl(kidia_var_2321 : kfdia_var_2322, 35, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_fac00_var_2324(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_fac01_var_2325(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_fac10_var_2326(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_fac11_var_2327(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_forfac_var_2328(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_forfrac_var_2329(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jp_var_2330(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jt_var_2331(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jt1_var_2332(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_colh2o_var_2333(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_colco2_var_2334(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_colo3_var_2335(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_coln2o(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_colch4_var_2336(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_colo2_var_2337(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_co2mult_var_2338(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  INTEGER(KIND = 4), INTENT(OUT) :: k_laytrop_var_2339(kidia_var_2321 : kfdia_var_2322)
  INTEGER(KIND = 4), INTENT(OUT) :: k_layswtch_var_2340(kidia_var_2321 : kfdia_var_2322)
  INTEGER(KIND = 4), INTENT(OUT) :: k_laylow_var_2341(kidia_var_2321 : kfdia_var_2322)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_2342(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(IN) :: p_tavel(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_selffac_var_2343(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_selffrac_var_2344(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indself_var_2345(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indfor_var_2346(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indminor(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_scaleminor(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_scaleminorn2(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: p_minorfrac(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  REAL(KIND = 8), INTENT(OUT) :: prat_h2oco2_var_2347(kidia_var_2321 : kfdia_var_2322, klev_var_2323), prat_h2oco2_1_var_2348(kidia_var_2321 : kfdia_var_2322, klev_var_2323), prat_h2oo3_var_2349(kidia_var_2321 : kfdia_var_2322, klev_var_2323), prat_h2oo3_1_var_2350(kidia_var_2321 : kfdia_var_2322, klev_var_2323), prat_h2on2o_var_2351(kidia_var_2321 : kfdia_var_2322, klev_var_2323), prat_h2on2o_1_var_2352(kidia_var_2321 : kfdia_var_2322, klev_var_2323), prat_h2och4_var_2353(kidia_var_2321 : kfdia_var_2322, klev_var_2323), prat_h2och4_1_var_2354(kidia_var_2321 : kfdia_var_2322, klev_var_2323), prat_n2oco2_var_2355(kidia_var_2321 : kfdia_var_2322, klev_var_2323), prat_n2oco2_1_var_2356(kidia_var_2321 : kfdia_var_2322, klev_var_2323), prat_o3co2_var_2357(kidia_var_2321 : kfdia_var_2322, klev_var_2323), prat_o3co2_1_var_2358(kidia_var_2321 : kfdia_var_2322, klev_var_2323)
  INTEGER(KIND = 4) :: jp1_var_2359, jlay
  INTEGER(KIND = 4) :: jlon_var_2360
  REAL(KIND = 8) :: z_co2reg_var_2361, z_compfp_var_2362, z_factor_var_2363, z_fp_var_2364, z_ft_var_2365, z_ft1_var_2366, z_plog_var_2367, z_scalefac_var_2368, z_stpfac_var_2369, z_water_var_2370
  DO jlon_var_2360 = kidia_var_2321, kfdia_var_2322
    z_stpfac_var_2369 = 0.29220138203356366D0
    k_laytrop_var_2339(jlon_var_2360) = 0
    k_layswtch_var_2340(jlon_var_2360) = 0
    k_laylow_var_2341(jlon_var_2360) = 0
    DO jlay = 1, klev_var_2323
      z_plog_var_2367 = LOG(pavel_var_2342(jlon_var_2360, jlay))
      k_jp_var_2330(jlon_var_2360, jlay) = INT(36.0D0 - 5 * (z_plog_var_2367 + 0.04D0))
      IF (k_jp_var_2330(jlon_var_2360, jlay) < 1) THEN
        k_jp_var_2330(jlon_var_2360, jlay) = 1
      ELSE IF (k_jp_var_2330(jlon_var_2360, jlay) > 58) THEN
        k_jp_var_2330(jlon_var_2360, jlay) = 58
      END IF
      jp1_var_2359 = k_jp_var_2330(jlon_var_2360, jlay) + 1
      z_fp_var_2364 = 5.0D0 * (preflog_var_214(k_jp_var_2330(jlon_var_2360, jlay)) - z_plog_var_2367)
      z_fp_var_2364 = MAX(- 1.0D0, MIN(1.0D0, z_fp_var_2364))
      k_jt_var_2331(jlon_var_2360, jlay) = INT(3.0D0 + (p_tavel(jlon_var_2360, jlay) - tref_var_215(k_jp_var_2330(jlon_var_2360, jlay))) / 15.0D0)
      IF (k_jt_var_2331(jlon_var_2360, jlay) < 1) THEN
        k_jt_var_2331(jlon_var_2360, jlay) = 1
      ELSE IF (k_jt_var_2331(jlon_var_2360, jlay) > 4) THEN
        k_jt_var_2331(jlon_var_2360, jlay) = 4
      END IF
      z_ft_var_2365 = ((p_tavel(jlon_var_2360, jlay) - tref_var_215(k_jp_var_2330(jlon_var_2360, jlay))) / 15.0D0) - REAL(k_jt_var_2331(jlon_var_2360, jlay) - 3)
      k_jt1_var_2332(jlon_var_2360, jlay) = INT(3.0D0 + (p_tavel(jlon_var_2360, jlay) - tref_var_215(jp1_var_2359)) / 15.0D0)
      IF (k_jt1_var_2332(jlon_var_2360, jlay) < 1) THEN
        k_jt1_var_2332(jlon_var_2360, jlay) = 1
      ELSE IF (k_jt1_var_2332(jlon_var_2360, jlay) > 4) THEN
        k_jt1_var_2332(jlon_var_2360, jlay) = 4
      END IF
      z_ft1_var_2366 = ((p_tavel(jlon_var_2360, jlay) - tref_var_215(jp1_var_2359)) / 15.0D0) - REAL(k_jt1_var_2332(jlon_var_2360, jlay) - 3)
      z_water_var_2370 = p_wkl(jlon_var_2360, 1, jlay) / p_coldry(jlon_var_2360, jlay)
      z_scalefac_var_2368 = pavel_var_2342(jlon_var_2360, jlay) * z_stpfac_var_2369 / p_tavel(jlon_var_2360, jlay)
      IF (z_plog_var_2367 > 4.56D0) THEN
        k_laytrop_var_2339(jlon_var_2360) = k_laytrop_var_2339(jlon_var_2360) + 1
        p_forfac_var_2328(jlon_var_2360, jlay) = z_scalefac_var_2368 / (1.0D0 + z_water_var_2370)
        z_factor_var_2363 = (332.0D0 - p_tavel(jlon_var_2360, jlay)) / 36.0D0
        k_indfor_var_2346(jlon_var_2360, jlay) = MIN(2, MAX(1, INT(z_factor_var_2363)))
        p_forfrac_var_2329(jlon_var_2360, jlay) = z_factor_var_2363 - REAL(k_indfor_var_2346(jlon_var_2360, jlay))
        p_selffac_var_2343(jlon_var_2360, jlay) = z_water_var_2370 * p_forfac_var_2328(jlon_var_2360, jlay)
        z_factor_var_2363 = (p_tavel(jlon_var_2360, jlay) - 188.0D0) / 7.2D0
        k_indself_var_2345(jlon_var_2360, jlay) = MIN(9, MAX(1, INT(z_factor_var_2363) - 7))
        p_selffrac_var_2344(jlon_var_2360, jlay) = z_factor_var_2363 - REAL(k_indself_var_2345(jlon_var_2360, jlay) + 7)
        p_scaleminor(jlon_var_2360, jlay) = pavel_var_2342(jlon_var_2360, jlay) / p_tavel(jlon_var_2360, jlay)
        p_scaleminorn2(jlon_var_2360, jlay) = (pavel_var_2342(jlon_var_2360, jlay) / p_tavel(jlon_var_2360, jlay)) * (p_wbroad(jlon_var_2360, jlay) / (p_coldry(jlon_var_2360, jlay) + p_wkl(jlon_var_2360, 1, jlay)))
        z_factor_var_2363 = (p_tavel(jlon_var_2360, jlay) - 180.8D0) / 7.2D0
        k_indminor(jlon_var_2360, jlay) = MIN(18, MAX(1, INT(z_factor_var_2363)))
        p_minorfrac(jlon_var_2360, jlay) = z_factor_var_2363 - REAL(k_indminor(jlon_var_2360, jlay))
        prat_h2oco2_var_2347(jlon_var_2360, jlay) = chi_mls(1, k_jp_var_2330(jlon_var_2360, jlay)) / chi_mls(2, k_jp_var_2330(jlon_var_2360, jlay))
        prat_h2oco2_1_var_2348(jlon_var_2360, jlay) = chi_mls(1, k_jp_var_2330(jlon_var_2360, jlay) + 1) / chi_mls(2, k_jp_var_2330(jlon_var_2360, jlay) + 1)
        prat_h2oo3_var_2349(jlon_var_2360, jlay) = chi_mls(1, k_jp_var_2330(jlon_var_2360, jlay)) / chi_mls(3, k_jp_var_2330(jlon_var_2360, jlay))
        prat_h2oo3_1_var_2350(jlon_var_2360, jlay) = chi_mls(1, k_jp_var_2330(jlon_var_2360, jlay) + 1) / chi_mls(3, k_jp_var_2330(jlon_var_2360, jlay) + 1)
        prat_h2on2o_var_2351(jlon_var_2360, jlay) = chi_mls(1, k_jp_var_2330(jlon_var_2360, jlay)) / chi_mls(4, k_jp_var_2330(jlon_var_2360, jlay))
        prat_h2on2o_1_var_2352(jlon_var_2360, jlay) = chi_mls(1, k_jp_var_2330(jlon_var_2360, jlay) + 1) / chi_mls(4, k_jp_var_2330(jlon_var_2360, jlay) + 1)
        prat_h2och4_var_2353(jlon_var_2360, jlay) = chi_mls(1, k_jp_var_2330(jlon_var_2360, jlay)) / chi_mls(6, k_jp_var_2330(jlon_var_2360, jlay))
        prat_h2och4_1_var_2354(jlon_var_2360, jlay) = chi_mls(1, k_jp_var_2330(jlon_var_2360, jlay) + 1) / chi_mls(6, k_jp_var_2330(jlon_var_2360, jlay) + 1)
        prat_n2oco2_var_2355(jlon_var_2360, jlay) = chi_mls(4, k_jp_var_2330(jlon_var_2360, jlay)) / chi_mls(2, k_jp_var_2330(jlon_var_2360, jlay))
        prat_n2oco2_1_var_2356(jlon_var_2360, jlay) = chi_mls(4, k_jp_var_2330(jlon_var_2360, jlay) + 1) / chi_mls(2, k_jp_var_2330(jlon_var_2360, jlay) + 1)
        p_colh2o_var_2333(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 1, jlay)
        p_colco2_var_2334(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 2, jlay)
        p_colo3_var_2335(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 3, jlay)
        p_coln2o(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 4, jlay)
        p_colch4_var_2336(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 6, jlay)
        p_colo2_var_2337(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 7, jlay)
        p_colbrd(jlon_var_2360, jlay) = 1D-20 * p_wbroad(jlon_var_2360, jlay)
        IF (p_colco2_var_2334(jlon_var_2360, jlay) == 0.0D0) p_colco2_var_2334(jlon_var_2360, jlay) = 1D-32 * p_coldry(jlon_var_2360, jlay)
        IF (p_coln2o(jlon_var_2360, jlay) == 0.0D0) p_coln2o(jlon_var_2360, jlay) = 1D-32 * p_coldry(jlon_var_2360, jlay)
        IF (p_colch4_var_2336(jlon_var_2360, jlay) == 0.0D0) p_colch4_var_2336(jlon_var_2360, jlay) = 1D-32 * p_coldry(jlon_var_2360, jlay)
        z_co2reg_var_2361 = 3.55D-24 * p_coldry(jlon_var_2360, jlay)
        p_co2mult_var_2338(jlon_var_2360, jlay) = (p_colco2_var_2334(jlon_var_2360, jlay) - z_co2reg_var_2361) * 272.63D0 * EXP(- 1919.4D0 / p_tavel(jlon_var_2360, jlay)) / (0.00087604D0 * p_tavel(jlon_var_2360, jlay))
      ELSE
        p_forfac_var_2328(jlon_var_2360, jlay) = z_scalefac_var_2368 / (1.0D0 + z_water_var_2370)
        z_factor_var_2363 = (p_tavel(jlon_var_2360, jlay) - 188.0D0) / 36.0D0
        k_indfor_var_2346(jlon_var_2360, jlay) = 3
        p_forfrac_var_2329(jlon_var_2360, jlay) = z_factor_var_2363 - 1.0D0
        p_selffac_var_2343(jlon_var_2360, jlay) = z_water_var_2370 * p_forfac_var_2328(jlon_var_2360, jlay)
        p_scaleminor(jlon_var_2360, jlay) = pavel_var_2342(jlon_var_2360, jlay) / p_tavel(jlon_var_2360, jlay)
        p_scaleminorn2(jlon_var_2360, jlay) = (pavel_var_2342(jlon_var_2360, jlay) / p_tavel(jlon_var_2360, jlay)) * (p_wbroad(jlon_var_2360, jlay) / (p_coldry(jlon_var_2360, jlay) + p_wkl(jlon_var_2360, 1, jlay)))
        z_factor_var_2363 = (p_tavel(jlon_var_2360, jlay) - 180.8D0) / 7.2D0
        k_indminor(jlon_var_2360, jlay) = MIN(18, MAX(1, INT(z_factor_var_2363)))
        p_minorfrac(jlon_var_2360, jlay) = z_factor_var_2363 - REAL(k_indminor(jlon_var_2360, jlay))
        prat_h2oco2_var_2347(jlon_var_2360, jlay) = chi_mls(1, k_jp_var_2330(jlon_var_2360, jlay)) / chi_mls(2, k_jp_var_2330(jlon_var_2360, jlay))
        prat_h2oco2_1_var_2348(jlon_var_2360, jlay) = chi_mls(1, k_jp_var_2330(jlon_var_2360, jlay) + 1) / chi_mls(2, k_jp_var_2330(jlon_var_2360, jlay) + 1)
        prat_o3co2_var_2357(jlon_var_2360, jlay) = chi_mls(3, k_jp_var_2330(jlon_var_2360, jlay)) / chi_mls(2, k_jp_var_2330(jlon_var_2360, jlay))
        prat_o3co2_1_var_2358(jlon_var_2360, jlay) = chi_mls(3, k_jp_var_2330(jlon_var_2360, jlay) + 1) / chi_mls(2, k_jp_var_2330(jlon_var_2360, jlay) + 1)
        p_colh2o_var_2333(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 1, jlay)
        p_colco2_var_2334(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 2, jlay)
        p_colo3_var_2335(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 3, jlay)
        p_coln2o(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 4, jlay)
        p_colch4_var_2336(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 6, jlay)
        p_colo2_var_2337(jlon_var_2360, jlay) = 1D-20 * p_wkl(jlon_var_2360, 7, jlay)
        p_colbrd(jlon_var_2360, jlay) = 1D-20 * p_wbroad(jlon_var_2360, jlay)
        IF (p_colco2_var_2334(jlon_var_2360, jlay) == 0.0D0) p_colco2_var_2334(jlon_var_2360, jlay) = 1D-32 * p_coldry(jlon_var_2360, jlay)
        IF (p_coln2o(jlon_var_2360, jlay) == 0.0D0) p_coln2o(jlon_var_2360, jlay) = 1D-32 * p_coldry(jlon_var_2360, jlay)
        IF (p_colch4_var_2336(jlon_var_2360, jlay) == 0.0D0) p_colch4_var_2336(jlon_var_2360, jlay) = 1D-32 * p_coldry(jlon_var_2360, jlay)
        z_co2reg_var_2361 = 3.55D-24 * p_coldry(jlon_var_2360, jlay)
        p_co2mult_var_2338(jlon_var_2360, jlay) = (p_colco2_var_2334(jlon_var_2360, jlay) - z_co2reg_var_2361) * 272.63D0 * EXP(- 1919.4D0 / p_tavel(jlon_var_2360, jlay)) / (0.00087604D0 * p_tavel(jlon_var_2360, jlay))
      END IF
      z_compfp_var_2362 = 1.0D0 - z_fp_var_2364
      p_fac10_var_2326(jlon_var_2360, jlay) = z_compfp_var_2362 * z_ft_var_2365
      p_fac00_var_2324(jlon_var_2360, jlay) = z_compfp_var_2362 * (1.0D0 - z_ft_var_2365)
      p_fac11_var_2327(jlon_var_2360, jlay) = z_fp_var_2364 * z_ft1_var_2366
      p_fac01_var_2325(jlon_var_2360, jlay) = z_fp_var_2364 * (1.0D0 - z_ft1_var_2366)
      p_selffac_var_2343(jlon_var_2360, jlay) = p_colh2o_var_2333(jlon_var_2360, jlay) * p_selffac_var_2343(jlon_var_2360, jlay)
      p_forfac_var_2328(jlon_var_2360, jlay) = p_colh2o_var_2333(jlon_var_2360, jlay) * p_forfac_var_2328(jlon_var_2360, jlay)
    END DO
    IF (k_laylow_var_2341(jlon_var_2360) == 0) k_laylow_var_2341(jlon_var_2360) = 1
  END DO
END SUBROUTINE rrtm_setcoef_140gp
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