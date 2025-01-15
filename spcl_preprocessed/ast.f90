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
  ELEMENTAL FUNCTION beta2alpha_fn_24(beta, frac1, frac2)
    IMPLICIT NONE
    REAL(KIND = 8), INTENT(IN) :: beta, frac1, frac2
    REAL(KIND = 8) :: beta2alpha_fn_24
    REAL(KIND = 8) :: frac_diff
    IF (beta < 1.0D0) THEN
      frac_diff = ABS(frac1 - frac2)
      beta2alpha_fn_24 = beta + (1.0D0 - beta) * frac_diff / (frac_diff + 1.0D0 / beta - 1.0D0)
    ELSE
      beta2alpha_fn_24 = 1.0D0
    END IF
  END FUNCTION beta2alpha_fn_24
  SUBROUTINE cum_cloud_cover_exp_ran(nlev_var_26, frac_var_27, overlap_param_var_28, cum_cloud_cover_var_29, pair_cloud_cover_var_30, is_beta_overlap)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlev_var_26
    REAL(KIND = 8), INTENT(IN) :: frac_var_27(nlev_var_26)
    REAL(KIND = 8), INTENT(IN) :: overlap_param_var_28(nlev_var_26 - 1)
    LOGICAL, INTENT(IN), OPTIONAL :: is_beta_overlap
    REAL(KIND = 8), INTENT(OUT) :: cum_cloud_cover_var_29(nlev_var_26)
    REAL(KIND = 8), INTENT(OUT) :: pair_cloud_cover_var_30(nlev_var_26 - 1)
    REAL(KIND = 8) :: cum_product
    REAL(KIND = 8) :: overlap_alpha
    LOGICAL :: do_overlap_conversion
    INTEGER :: jlev_var_31
    do_overlap_conversion = is_beta_overlap
    cum_product = 1.0D0 - frac_var_27(1)
    cum_cloud_cover_var_29(1) = frac_var_27(1)
    DO jlev_var_31 = 1, nlev_var_26 - 1
      IF (do_overlap_conversion) THEN
        overlap_alpha = beta2alpha_fn_24(overlap_param_var_28(jlev_var_31), frac_var_27(jlev_var_31), frac_var_27(jlev_var_31 + 1))
      ELSE
        overlap_alpha = overlap_param_var_28(jlev_var_31)
      END IF
      pair_cloud_cover_var_30(jlev_var_31) = overlap_alpha * MAX(frac_var_27(jlev_var_31), frac_var_27(jlev_var_31 + 1)) + (1.0D0 - overlap_alpha) * (frac_var_27(jlev_var_31) + frac_var_27(jlev_var_31 + 1) - frac_var_27(jlev_var_31) * frac_var_27(jlev_var_31 + 1))
      IF (frac_var_27(jlev_var_31) >= 0.9999999999999978D0) THEN
        cum_product = 0.0D0
      ELSE
        cum_product = cum_product * (1.0D0 - pair_cloud_cover_var_30(jlev_var_31)) / (1.0D0 - frac_var_27(jlev_var_31))
      END IF
      cum_cloud_cover_var_29(jlev_var_31 + 1) = 1.0D0 - cum_product
    END DO
  END SUBROUTINE cum_cloud_cover_exp_ran
END MODULE radiation_cloud_cover
MODULE radiation_ice_optics_fu
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE calc_ice_optics_fu_sw(nb_var_32, coeff_var_33, ice_wp_var_34, re_var_35, od_var_36, scat_od_var_37, g_var_38)
    INTEGER, INTENT(IN) :: nb_var_32
    REAL(KIND = 8), INTENT(IN) :: coeff_var_33(:, :)
    REAL(KIND = 8), INTENT(IN) :: ice_wp_var_34
    REAL(KIND = 8), INTENT(IN) :: re_var_35
    REAL(KIND = 8), INTENT(OUT) :: od_var_36(nb_var_32), scat_od_var_37(nb_var_32), g_var_38(nb_var_32)
    REAL(KIND = 8) :: de_um_var_39, inv_de_um_var_40
    REAL(KIND = 8) :: iwp_gm_2_var_41
    INTEGER :: jb_var_42
    de_um_var_39 = MIN(re_var_35, 0.0001D0) * (1539598.4727183152D0)
    inv_de_um_var_40 = 1.0D0 / de_um_var_39
    iwp_gm_2_var_41 = ice_wp_var_34 * 1000.0D0
    DO jb_var_42 = 1, 14
      od_var_36(jb_var_42) = iwp_gm_2_var_41 * (coeff_var_33(jb_var_42, 1) + coeff_var_33(jb_var_42, 2) * inv_de_um_var_40)
      scat_od_var_37(jb_var_42) = od_var_36(jb_var_42) * (1.0D0 - (coeff_var_33(jb_var_42, 3) + de_um_var_39 * (coeff_var_33(jb_var_42, 4) + de_um_var_39 * (coeff_var_33(jb_var_42, 5) + de_um_var_39 * coeff_var_33(jb_var_42, 6)))))
      g_var_38(jb_var_42) = MIN(coeff_var_33(jb_var_42, 7) + de_um_var_39 * (coeff_var_33(jb_var_42, 8) + de_um_var_39 * (coeff_var_33(jb_var_42, 9) + de_um_var_39 * coeff_var_33(jb_var_42, 10))), 0.9999999999999978D0)
    END DO
  END SUBROUTINE calc_ice_optics_fu_sw
  SUBROUTINE calc_ice_optics_fu_lw(nb_var_43, coeff_var_44, ice_wp_var_45, re_var_46, od_var_47, scat_od_var_48, g_var_49)
    INTEGER, INTENT(IN) :: nb_var_43
    REAL(KIND = 8), INTENT(IN) :: coeff_var_44(:, :)
    REAL(KIND = 8), INTENT(IN) :: ice_wp_var_45
    REAL(KIND = 8), INTENT(IN) :: re_var_46
    REAL(KIND = 8), INTENT(OUT) :: od_var_47(nb_var_43), scat_od_var_48(nb_var_43), g_var_49(nb_var_43)
    REAL(KIND = 8) :: de_um_var_50, inv_de_um_var_51
    REAL(KIND = 8) :: iwp_gm_2_var_52
    INTEGER :: jb_var_53
    de_um_var_50 = MIN(re_var_46, 0.0001D0) * (1539598.4727183152D0)
    inv_de_um_var_51 = 1.0D0 / de_um_var_50
    iwp_gm_2_var_52 = ice_wp_var_45 * 1000.0D0
    DO jb_var_53 = 1, 16
      od_var_47(jb_var_53) = iwp_gm_2_var_52 * (coeff_var_44(jb_var_53, 1) + inv_de_um_var_51 * (coeff_var_44(jb_var_53, 2) + inv_de_um_var_51 * coeff_var_44(jb_var_53, 3)))
      scat_od_var_48(jb_var_53) = od_var_47(jb_var_53) - iwp_gm_2_var_52 * inv_de_um_var_51 * (coeff_var_44(jb_var_53, 4) + de_um_var_50 * (coeff_var_44(jb_var_53, 5) + de_um_var_50 * (coeff_var_44(jb_var_53, 6) + de_um_var_50 * coeff_var_44(jb_var_53, 7))))
      g_var_49(jb_var_53) = MIN(coeff_var_44(jb_var_53, 8) + de_um_var_50 * (coeff_var_44(jb_var_53, 9) + de_um_var_50 * (coeff_var_44(jb_var_53, 10) + de_um_var_50 * coeff_var_44(jb_var_53, 11))), 0.9999999999999978D0)
    END DO
  END SUBROUTINE calc_ice_optics_fu_lw
END MODULE radiation_ice_optics_fu
MODULE radiation_ice_optics_yi
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE calc_ice_optics_yi_sw(nb_var_54, coeff_var_55, ice_wp_var_56, re_var_57, od_var_58, scat_od_var_59, g_var_60)
    INTEGER, INTENT(IN) :: nb_var_54
    REAL(KIND = 8), INTENT(IN) :: coeff_var_55(:, :)
    REAL(KIND = 8), INTENT(IN) :: ice_wp_var_56
    REAL(KIND = 8), INTENT(IN) :: re_var_57
    REAL(KIND = 8), INTENT(OUT) :: od_var_58(nb_var_54), scat_od_var_59(nb_var_54), g_var_60(nb_var_54)
    REAL(KIND = 8) :: de_um_var_61
    REAL(KIND = 8) :: iwp_gm_2_var_62
    REAL(KIND = 8) :: wts_1_var_63, wts_2_var_64
    INTEGER(KIND = 4) :: lu_idx_var_65
    de_um_var_61 = re_var_57 * 2000000.0D0
    de_um_var_61 = MAX(de_um_var_61, 10.0D0)
    de_um_var_61 = MIN(de_um_var_61, 119.99D0)
    iwp_gm_2_var_62 = ice_wp_var_56 * 1000.0D0
    lu_idx_var_65 = FLOOR(de_um_var_61 * 0.2D0 - 1.0D0)
    wts_2_var_64 = (de_um_var_61 * 0.2D0 - 1.0D0) - lu_idx_var_65
    wts_1_var_63 = 1.0D0 - wts_2_var_64
    od_var_58 = 0.001D0 * iwp_gm_2_var_62 * (wts_1_var_63 * coeff_var_55(1 : 14, lu_idx_var_65) + wts_2_var_64 * coeff_var_55(1 : 14, lu_idx_var_65 + 1))
    scat_od_var_59 = od_var_58 * (wts_1_var_63 * coeff_var_55(1 : 14, lu_idx_var_65 + 23) + wts_2_var_64 * coeff_var_55(1 : 14, lu_idx_var_65 + 23 + 1))
    g_var_60 = wts_1_var_63 * coeff_var_55(1 : 14, lu_idx_var_65 + 46) + wts_2_var_64 * coeff_var_55(1 : 14, lu_idx_var_65 + 46 + 1)
  END SUBROUTINE calc_ice_optics_yi_sw
  SUBROUTINE calc_ice_optics_yi_lw(nb_var_66, coeff_var_67, ice_wp_var_68, re_var_69, od_var_70, scat_od_var_71, g_var_72)
    INTEGER, INTENT(IN) :: nb_var_66
    REAL(KIND = 8), INTENT(IN) :: coeff_var_67(:, :)
    REAL(KIND = 8), INTENT(IN) :: ice_wp_var_68
    REAL(KIND = 8), INTENT(IN) :: re_var_69
    REAL(KIND = 8), INTENT(OUT) :: od_var_70(nb_var_66), scat_od_var_71(nb_var_66), g_var_72(nb_var_66)
    REAL(KIND = 8) :: de_um_var_73
    REAL(KIND = 8) :: iwp_gm_2_var_74
    REAL(KIND = 8) :: wts_1_var_75, wts_2_var_76
    INTEGER(KIND = 4) :: lu_idx_var_77
    de_um_var_73 = re_var_69 * 2000000.0D0
    de_um_var_73 = MAX(de_um_var_73, 10.0D0)
    de_um_var_73 = MIN(de_um_var_73, 119.99D0)
    iwp_gm_2_var_74 = ice_wp_var_68 * 1000.0D0
    lu_idx_var_77 = FLOOR(de_um_var_73 * 0.2D0 - 1.0D0)
    wts_2_var_76 = (de_um_var_73 * 0.2D0 - 1.0D0) - lu_idx_var_77
    wts_1_var_75 = 1.0D0 - wts_2_var_76
    od_var_70 = 0.001D0 * iwp_gm_2_var_74 * (wts_1_var_75 * coeff_var_67(1 : 16, lu_idx_var_77) + wts_2_var_76 * coeff_var_67(1 : 16, lu_idx_var_77 + 1))
    scat_od_var_71 = od_var_70 * (wts_1_var_75 * coeff_var_67(1 : 16, lu_idx_var_77 + 23) + wts_2_var_76 * coeff_var_67(1 : 16, lu_idx_var_77 + 23 + 1))
    g_var_72 = wts_1_var_75 * coeff_var_67(1 : 16, lu_idx_var_77 + 46) + wts_2_var_76 * coeff_var_67(1 : 16, lu_idx_var_77 + 46 + 1)
  END SUBROUTINE calc_ice_optics_yi_lw
END MODULE radiation_ice_optics_yi
MODULE radiation_liquid_optics_socrates
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE calc_liq_optics_socrates(nb_var_78, coeff_var_79, lwp, re_in, od_var_80, scat_od_var_81, g_var_82)
    INTEGER, INTENT(IN) :: nb_var_78
    REAL(KIND = 8), INTENT(IN) :: coeff_var_79(:, :)
    REAL(KIND = 8), INTENT(IN) :: lwp, re_in
    REAL(KIND = 8), INTENT(OUT) :: od_var_80(nb_var_78), scat_od_var_81(nb_var_78), g_var_82(nb_var_78)
    INTEGER :: jb_var_83
    REAL(KIND = 8) :: re_var_84
    re_var_84 = MAX(1.2000000424450263D-06, MIN(re_in, 4.999999873689376D-05))
    DO jb_var_83 = 1, nb_var_78
      od_var_80(jb_var_83) = lwp * (coeff_var_79(jb_var_83, 1) + re_var_84 * (coeff_var_79(jb_var_83, 2) + re_var_84 * coeff_var_79(jb_var_83, 3))) / (1.0D0 + re_var_84 * (coeff_var_79(jb_var_83, 4) + re_var_84 * (coeff_var_79(jb_var_83, 5) + re_var_84 * coeff_var_79(jb_var_83, 6))))
      scat_od_var_81(jb_var_83) = od_var_80(jb_var_83) * (1.0D0 - (coeff_var_79(jb_var_83, 7) + re_var_84 * (coeff_var_79(jb_var_83, 8) + re_var_84 * coeff_var_79(jb_var_83, 9))) / (1.0D0 + re_var_84 * (coeff_var_79(jb_var_83, 10) + re_var_84 * coeff_var_79(jb_var_83, 11))))
      g_var_82(jb_var_83) = (coeff_var_79(jb_var_83, 12) + re_var_84 * (coeff_var_79(jb_var_83, 13) + re_var_84 * coeff_var_79(jb_var_83, 14))) / (1.0D0 + re_var_84 * (coeff_var_79(jb_var_83, 15) + re_var_84 * coeff_var_79(jb_var_83, 16)))
    END DO
  END SUBROUTINE calc_liq_optics_socrates
END MODULE radiation_liquid_optics_socrates
MODULE radiation_random_numbers
  IMPLICIT NONE
  TYPE :: rng_type
    INTEGER(KIND = 4) :: itype = 0
    REAL(KIND = 8) :: istate(512)
    INTEGER(KIND = 4) :: nmaxstreams = 512
    INTEGER(KIND = 4) :: iseed
  END TYPE rng_type
  CONTAINS
  SUBROUTINE initialize(this_var_88, itype_var_89, iseed_var_90, nmaxstreams_var_91)
    CLASS(rng_type), INTENT(INOUT) :: this_var_88
    INTEGER(KIND = 4), INTENT(IN), OPTIONAL :: itype_var_89
    INTEGER(KIND = 4), INTENT(IN), OPTIONAL :: iseed_var_90
    INTEGER(KIND = 4), INTENT(IN), OPTIONAL :: nmaxstreams_var_91
    INTEGER :: jstr
    REAL(KIND = 8) :: rseed
    this_var_88 % itype = 1
    this_var_88 % iseed = iseed_var_90
    this_var_88 % nmaxstreams = nmaxstreams_var_91
    rseed = REAL(ABS(this_var_88 % iseed), 8)
    DO jstr = 1, this_var_88 % nmaxstreams
      this_var_88 % istate(jstr) = NINT(MOD(rseed * jstr * (1.0D0 - 0.05D0 * jstr + 0.005D0 * jstr ** 2) * 16807.0D0, 2147483647.0D0), kind = 8)
    END DO
    DO jstr = 1, this_var_88 % nmaxstreams
      this_var_88 % istate(jstr) = MOD(48271.0D0 * this_var_88 % istate(jstr), 2147483647.0D0)
    END DO
  END SUBROUTINE initialize
  SUBROUTINE uniform_distribution_1d(this_var_92, randnum_var_93)
    CLASS(rng_type), INTENT(INOUT) :: this_var_92
    REAL(KIND = 8), INTENT(OUT) :: randnum_var_93(:)
    INTEGER :: imax_var_94, i_var_95
    IF (this_var_92 % itype == 1) THEN
      imax_var_94 = MIN(this_var_92 % nmaxstreams, SIZE(randnum_var_93))
      DO i_var_95 = 1, imax_var_94
        this_var_92 % istate(i_var_95) = MOD(48271.0D0 * this_var_92 % istate(i_var_95), 2147483647.0D0)
        randnum_var_93(i_var_95) = 4.656612875245797D-10 * this_var_92 % istate(i_var_95)
      END DO
    ELSE
      CALL RANDOM_NUMBER(randnum_var_93)
    END IF
  END SUBROUTINE uniform_distribution_1d
  SUBROUTINE uniform_distribution_2d(this_var_96, randnum_var_97)
    CLASS(rng_type), INTENT(INOUT) :: this_var_96
    REAL(KIND = 8), INTENT(OUT) :: randnum_var_97(:, :)
    INTEGER :: imax_var_98, jblock_var_99, i_var_100
    IF (this_var_96 % itype == 1) THEN
      imax_var_98 = MIN(this_var_96 % nmaxstreams, SIZE(randnum_var_97, 1))
      DO jblock_var_99 = 1, SIZE(randnum_var_97, 2)
        DO i_var_100 = 1, imax_var_98
          this_var_96 % istate(i_var_100) = MOD(48271.0D0 * this_var_96 % istate(i_var_100), 2147483647.0D0)
          randnum_var_97(i_var_100, jblock_var_99) = 4.656612875245797D-10 * this_var_96 % istate(i_var_100)
        END DO
      END DO
    ELSE
      CALL RANDOM_NUMBER(randnum_var_97)
    END IF
  END SUBROUTINE uniform_distribution_2d
  SUBROUTINE uniform_distribution_2d_masked(this_var_101, randnum_var_102, mask_var_103)
    CLASS(rng_type), INTENT(INOUT) :: this_var_101
    REAL(KIND = 8), INTENT(INOUT) :: randnum_var_102(:, :)
    LOGICAL, INTENT(IN) :: mask_var_103(:)
    INTEGER :: imax_var_104, jblock_var_105, i_var_106
    IF (this_var_101 % itype == 1) THEN
      imax_var_104 = MIN(this_var_101 % nmaxstreams, SIZE(randnum_var_102, 1))
      DO jblock_var_105 = 1, SIZE(randnum_var_102, 2)
        IF (mask_var_103(jblock_var_105)) THEN
          DO i_var_106 = 1, imax_var_104
            this_var_101 % istate(i_var_106) = MOD(48271.0D0 * this_var_101 % istate(i_var_106), 2147483647.0D0)
            randnum_var_102(i_var_106, jblock_var_105) = 4.656612875245797D-10 * this_var_101 % istate(i_var_106)
          END DO
        END IF
      END DO
    ELSE
      DO jblock_var_105 = 1, SIZE(randnum_var_102, 2)
        IF (mask_var_103(jblock_var_105)) THEN
          CALL RANDOM_NUMBER(randnum_var_102(:, jblock_var_105))
        END IF
      END DO
    END IF
  END SUBROUTINE uniform_distribution_2d_masked
END MODULE radiation_random_numbers
MODULE radiation_two_stream
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE calc_no_scattering_transmittance_lw(ng_var_107, od_var_108, planck_top, planck_bot, transmittance_var_109, source_up_var_110, source_dn_var_111)
    INTEGER, INTENT(IN) :: ng_var_107
    REAL(KIND = 8), INTENT(IN), DIMENSION(ng_var_107) :: od_var_108
    REAL(KIND = 8), INTENT(IN), DIMENSION(ng_var_107) :: planck_top, planck_bot
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ng_var_107) :: transmittance_var_109
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ng_var_107) :: source_up_var_110, source_dn_var_111
    REAL(KIND = 8) :: coeff_var_112, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot
    INTEGER :: jg_var_113
    transmittance_var_109 = EXP(- 1.66D0 * od_var_108)
    DO jg_var_113 = 1, ng_var_107
      coeff_var_112 = 1.66D0 * od_var_108(jg_var_113)
      IF (od_var_108(jg_var_113) > 0.001D0) THEN
        coeff_var_112 = (planck_bot(jg_var_113) - planck_top(jg_var_113)) / coeff_var_112
        coeff_up_top = coeff_var_112 + planck_top(jg_var_113)
        coeff_up_bot = coeff_var_112 + planck_bot(jg_var_113)
        coeff_dn_top = - coeff_var_112 + planck_top(jg_var_113)
        coeff_dn_bot = - coeff_var_112 + planck_bot(jg_var_113)
        source_up_var_110(jg_var_113) = coeff_up_top - transmittance_var_109(jg_var_113) * coeff_up_bot
        source_dn_var_111(jg_var_113) = coeff_dn_bot - transmittance_var_109(jg_var_113) * coeff_dn_top
      ELSE
        source_up_var_110(jg_var_113) = coeff_var_112 * 0.5D0 * (planck_top(jg_var_113) + planck_bot(jg_var_113))
        source_dn_var_111(jg_var_113) = source_up_var_110(jg_var_113)
      END IF
    END DO
  END SUBROUTINE calc_no_scattering_transmittance_lw
  SUBROUTINE calc_ref_trans_sw(ng_var_114, mu0, od_var_115, ssa_var_116, asymmetry, ref_diff, trans_diff, ref_dir_var_117, trans_dir_diff_var_118, trans_dir_dir_var_119)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng_var_114
    REAL(KIND = 8), INTENT(IN) :: mu0
    REAL(KIND = 8), INTENT(IN), DIMENSION(ng_var_114) :: od_var_115, ssa_var_116, asymmetry
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ng_var_114) :: ref_dir_var_117, trans_dir_diff_var_118
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ng_var_114) :: ref_diff, trans_diff
    REAL(KIND = 8), INTENT(OUT), DIMENSION(ng_var_114) :: trans_dir_dir_var_119
    REAL(KIND = 8), DIMENSION(ng_var_114) :: gamma1, gamma2, gamma3, gamma4
    REAL(KIND = 8), DIMENSION(ng_var_114) :: alpha1, alpha2, k_exponent
    REAL(KIND = 8), DIMENSION(ng_var_114) :: exponential
    REAL(KIND = 8) :: reftrans_factor, factor_var_120
    REAL(KIND = 8) :: exponential2
    REAL(KIND = 8) :: k_mu0, k_gamma3, k_gamma4
    REAL(KIND = 8) :: k_2_exponential, one_minus_kmu0_sqr
    INTEGER :: jg_var_121
    trans_dir_dir_var_119 = MAX(- MAX(od_var_115 * (1.0D0 / mu0), 0.0D0), - 1000.0D0)
    trans_dir_dir_var_119 = EXP(trans_dir_dir_var_119)
    DO jg_var_121 = 1, ng_var_114
      factor_var_120 = 0.75D0 * asymmetry(jg_var_121)
      gamma1(jg_var_121) = 2.0D0 - ssa_var_116(jg_var_121) * (1.25D0 + factor_var_120)
      gamma2(jg_var_121) = ssa_var_116(jg_var_121) * (0.75D0 - factor_var_120)
      gamma3(jg_var_121) = 0.5D0 - mu0 * factor_var_120
      gamma4(jg_var_121) = 1.0D0 - gamma3(jg_var_121)
      alpha1(jg_var_121) = gamma1(jg_var_121) * gamma4(jg_var_121) + gamma2(jg_var_121) * gamma3(jg_var_121)
      alpha2(jg_var_121) = gamma1(jg_var_121) * gamma3(jg_var_121) + gamma2(jg_var_121) * gamma4(jg_var_121)
      k_exponent(jg_var_121) = SQRT(MAX((gamma1(jg_var_121) - gamma2(jg_var_121)) * (gamma1(jg_var_121) + gamma2(jg_var_121)), 1D-12))
    END DO
    exponential = EXP(- k_exponent * od_var_115)
    DO jg_var_121 = 1, ng_var_114
      k_mu0 = k_exponent(jg_var_121) * mu0
      one_minus_kmu0_sqr = 1.0D0 - k_mu0 * k_mu0
      k_gamma3 = k_exponent(jg_var_121) * gamma3(jg_var_121)
      k_gamma4 = k_exponent(jg_var_121) * gamma4(jg_var_121)
      exponential2 = exponential(jg_var_121) * exponential(jg_var_121)
      k_2_exponential = 2.0D0 * k_exponent(jg_var_121) * exponential(jg_var_121)
      reftrans_factor = 1.0D0 / (k_exponent(jg_var_121) + gamma1(jg_var_121) + (k_exponent(jg_var_121) - gamma1(jg_var_121)) * exponential2)
      ref_diff(jg_var_121) = gamma2(jg_var_121) * (1.0D0 - exponential2) * reftrans_factor
      trans_diff(jg_var_121) = MAX(0.0D0, MIN(k_2_exponential * reftrans_factor, 1.0D0 - ref_diff(jg_var_121)))
      reftrans_factor = mu0 * ssa_var_116(jg_var_121) * reftrans_factor / MERGE(one_minus_kmu0_sqr, 2.220446049250313D-16, ABS(one_minus_kmu0_sqr) > 2.220446049250313D-16)
      ref_dir_var_117(jg_var_121) = reftrans_factor * ((1.0D0 - k_mu0) * (alpha2(jg_var_121) + k_gamma3) - (1.0D0 + k_mu0) * (alpha2(jg_var_121) - k_gamma3) * exponential2 - k_2_exponential * (gamma3(jg_var_121) - alpha2(jg_var_121) * mu0) * trans_dir_dir_var_119(jg_var_121))
      trans_dir_diff_var_118(jg_var_121) = reftrans_factor * (k_2_exponential * (gamma4(jg_var_121) + alpha1(jg_var_121) * mu0) - trans_dir_dir_var_119(jg_var_121) * ((1.0D0 + k_mu0) * (alpha1(jg_var_121) + k_gamma4) - (1.0D0 - k_mu0) * (alpha1(jg_var_121) - k_gamma4) * exponential2))
      ref_dir_var_117(jg_var_121) = MAX(0.0D0, MIN(ref_dir_var_117(jg_var_121), mu0 * (1.0D0 - trans_dir_dir_var_119(jg_var_121))))
      trans_dir_diff_var_118(jg_var_121) = MAX(0.0D0, MIN(trans_dir_diff_var_118(jg_var_121), mu0 * (1.0D0 - trans_dir_dir_var_119(jg_var_121)) - ref_dir_var_117(jg_var_121)))
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
  SUBROUTINE initialize_random_numbers(kseed, yd_stream_var_122)
    INTEGER(KIND = 4), INTENT(IN) :: kseed
    TYPE(randomnumberstream), INTENT(INOUT) :: yd_stream_var_122
    INTEGER(KIND = 4) :: idum, jj_var_123, jbit
    REAL(KIND = 8), DIMENSION(999) :: zwarmup
    idum = ABS(IEOR(kseed, 123459876))
    IF (idum == 0) idum = 123459876
    DO jj_var_123 = 1, 64
      IF (BTEST(idum, 31)) THEN
        idum = IBSET(ISHFT(IEOR(idum, 87), 1), 0)
      ELSE
        idum = IBCLR(ISHFT(idum, 1), 0)
      END IF
    END DO
    yd_stream_var_122 % ix(1 : 606) = 0
    yd_stream_var_122 % ix(2) = ISHFT(IBITS(idum, 0, 29), 1)
    yd_stream_var_122 % ix(jpq) = IBITS(idum, 29, BIT_SIZE(idum) + 1 - 30)
    DO jbit = 1, 29
      DO jj_var_123 = 3, 606
        IF (BTEST(idum, 31)) THEN
          idum = IBSET(ISHFT(IEOR(idum, 87), 1), 0)
          yd_stream_var_122 % ix(jj_var_123) = IBSET(yd_stream_var_122 % ix(jj_var_123), jbit)
        ELSE
          idum = IBCLR(ISHFT(idum, 1), 0)
        END IF
      END DO
    END DO
    yd_stream_var_122 % ix(502) = IBSET(yd_stream_var_122 % ix(502), 0)
    yd_stream_var_122 % iused = 607
    yd_stream_var_122 % zrm = 9.313225746154785D-10
    yd_stream_var_122 % inittest = 12345678
    CALL uniform_distribution(zwarmup, yd_stream_var_122)
  END SUBROUTINE initialize_random_numbers
  SUBROUTINE uniform_distribution(px, yd_stream_var_124)
    TYPE(randomnumberstream), INTENT(INOUT) :: yd_stream_var_124
    REAL(KIND = 8), DIMENSION(:), INTENT(OUT) :: px
    INTEGER(KIND = 4) :: jj_var_125, jk_var_126, in, ifilled
    IF (yd_stream_var_124 % inittest /= 12345678) CALL abor1('uniform_distribution called before initialize_random_numbers')
    in = SIZE(px)
    ifilled = 0
    DO jj_var_125 = yd_stream_var_124 % iused + 1, MIN(607, in + yd_stream_var_124 % iused)
      px(jj_var_125 - yd_stream_var_124 % iused) = yd_stream_var_124 % ix(jj_var_125) * yd_stream_var_124 % zrm
      ifilled = 1
    END DO
    yd_stream_var_124 % iused = yd_stream_var_124 % iused + 1
    IF (1 == in) THEN
      RETURN
    END IF
    DO WHILE (ifilled < in)
      DO jj_var_125 = 1, 273
        yd_stream_var_124 % ix(jj_var_125) = IAND(1073741823, yd_stream_var_124 % ix(jj_var_125) + yd_stream_var_124 % ix(jj_var_125 - 273 + 607))
      END DO
      DO jk_var_126 = 1, 2
        DO jj_var_125 = 274 + (jk_var_126 - 1) * 167, MIN(607, 273 + jk_var_126 * 167)
          yd_stream_var_124 % ix(jj_var_125) = IAND(1073741823, yd_stream_var_124 % ix(jj_var_125) + yd_stream_var_124 % ix(jj_var_125 - 273))
        END DO
      END DO
      yd_stream_var_124 % iused = MIN(607, in - 1)
      px(ifilled + 1 : ifilled + yd_stream_var_124 % iused) = yd_stream_var_124 % ix(1 : yd_stream_var_124 % iused) * yd_stream_var_124 % zrm
      ifilled = 1 + yd_stream_var_124 % iused
    END DO
  END SUBROUTINE uniform_distribution
END MODULE random_numbers_mix
MODULE yoerrta1
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_127(10), fracrefb_var_128(10)
  REAL(KIND = 8) :: absa_var_129(65, 10)
  REAL(KIND = 8) :: absb_var_130(235, 10)
  REAL(KIND = 8) :: ka_mn2_var_131(19, 10), kb_mn2(19, 10)
  REAL(KIND = 8) :: selfref_var_132(10, 10), forref_var_133(4, 10)
END MODULE yoerrta1
MODULE yoerrta10
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(6) :: fracrefa_var_134
  REAL(KIND = 8), DIMENSION(6) :: fracrefb_var_135
  REAL(KIND = 8) :: absa_var_136(65, 6)
  REAL(KIND = 8) :: absb_var_137(235, 6)
  REAL(KIND = 8) :: selfref_var_138(10, 6)
  REAL(KIND = 8) :: forref_var_139(4, 6)
END MODULE yoerrta10
MODULE yoerrta11
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(8) :: fracrefa_var_140
  REAL(KIND = 8), DIMENSION(8) :: fracrefb_var_141
  REAL(KIND = 8) :: absa_var_142(65, 8)
  REAL(KIND = 8) :: absb_var_143(235, 8)
  REAL(KIND = 8) :: ka_mo2(19, 8)
  REAL(KIND = 8) :: kb_mo2(19, 8)
  REAL(KIND = 8) :: selfref_var_144(10, 8)
  REAL(KIND = 8) :: forref_var_145(4, 8)
END MODULE yoerrta11
MODULE yoerrta12
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_146(8, 9)
  REAL(KIND = 8) :: absa_var_147(585, 8)
  REAL(KIND = 8) :: selfref_var_148(10, 8)
  REAL(KIND = 8) :: forref_var_149(4, 8)
END MODULE yoerrta12
MODULE yoerrta13
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_150(4, 9)
  REAL(KIND = 8), DIMENSION(4) :: fracrefb_var_151
  REAL(KIND = 8) :: absa_var_152(585, 4)
  REAL(KIND = 8) :: selfref_var_153(10, 4)
  REAL(KIND = 8) :: forref_var_154(4, 4)
  REAL(KIND = 8) :: ka_mco2_var_155(9, 19, 4)
  REAL(KIND = 8) :: ka_mco(9, 19, 4)
  REAL(KIND = 8) :: kb_mo3(19, 4)
END MODULE yoerrta13
MODULE yoerrta14
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(2) :: fracrefa_var_156
  REAL(KIND = 8), DIMENSION(2) :: fracrefb_var_157
  REAL(KIND = 8) :: absa_var_158(65, 2)
  REAL(KIND = 8) :: absb_var_159(235, 2)
  REAL(KIND = 8) :: selfref_var_160(10, 2)
  REAL(KIND = 8) :: forref_var_161(4, 2)
END MODULE yoerrta14
MODULE yoerrta15
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_162(2, 9)
  REAL(KIND = 8) :: absa_var_163(585, 2)
  REAL(KIND = 8) :: ka_mn2_var_164(9, 19, 2)
  REAL(KIND = 8) :: selfref_var_165(10, 2)
  REAL(KIND = 8) :: forref_var_166(4, 2)
END MODULE yoerrta15
MODULE yoerrta16
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_167(2, 9)
  REAL(KIND = 8), DIMENSION(2) :: fracrefb_var_168
  REAL(KIND = 8) :: absa_var_169(585, 2)
  REAL(KIND = 8) :: absb_var_170(235, 2)
  REAL(KIND = 8) :: selfref_var_171(10, 2)
  REAL(KIND = 8) :: forref_var_172(4, 2)
END MODULE yoerrta16
MODULE yoerrta2
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_173(12), fracrefb_var_174(12)
  REAL(KIND = 8) :: absa_var_175(65, 12)
  REAL(KIND = 8) :: absb_var_176(235, 12)
  REAL(KIND = 8) :: selfref_var_177(10, 12), forref_var_178(4, 12)
END MODULE yoerrta2
MODULE yoerrta3
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_179(16, 9), fracrefb_var_180(16, 5)
  REAL(KIND = 8) :: ka_mn2o_var_181(9, 19, 16), kb_mn2o_var_182(5, 19, 16)
  REAL(KIND = 8) :: absa_var_183(585, 16)
  REAL(KIND = 8) :: absb_var_184(1175, 16)
  REAL(KIND = 8) :: selfref_var_185(10, 16)
  REAL(KIND = 8) :: forref_var_186(4, 16)
END MODULE yoerrta3
MODULE yoerrta4
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_187(14, 9), fracrefb_var_188(14, 5)
  REAL(KIND = 8) :: absa_var_189(585, 14)
  REAL(KIND = 8) :: absb_var_190(1175, 14)
  REAL(KIND = 8) :: selfref_var_191(10, 14), forref_var_192(4, 14)
END MODULE yoerrta4
MODULE yoerrta5
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_193(16, 9), fracrefb_var_194(16, 5)
  REAL(KIND = 8), DIMENSION(16) :: ccl4
  REAL(KIND = 8) :: absa_var_195(585, 16)
  REAL(KIND = 8) :: absb_var_196(1175, 16)
  REAL(KIND = 8) :: ka_mo3_var_197(9, 19, 16)
  REAL(KIND = 8) :: selfref_var_198(10, 16)
  REAL(KIND = 8) :: forref_var_199(4, 16)
END MODULE yoerrta5
MODULE yoerrta6
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(8) :: fracrefa_var_200
  REAL(KIND = 8), DIMENSION(8) :: cfc11adj
  REAL(KIND = 8), DIMENSION(8) :: cfc12_var_201
  REAL(KIND = 8) :: absa_var_202(65, 8)
  REAL(KIND = 8) :: selfref_var_203(10, 8)
  REAL(KIND = 8) :: ka_mco2_var_204(19, 8)
  REAL(KIND = 8) :: forref_var_205(4, 8)
END MODULE yoerrta6
MODULE yoerrta7
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_206(12, 9)
  REAL(KIND = 8), DIMENSION(12) :: fracrefb_var_207
  REAL(KIND = 8) :: absa_var_208(585, 12)
  REAL(KIND = 8) :: absb_var_209(235, 12)
  REAL(KIND = 8) :: selfref_var_210(10, 12)
  REAL(KIND = 8) :: ka_mco2_var_211(9, 19, 12)
  REAL(KIND = 8) :: kb_mco2_var_212(19, 12)
  REAL(KIND = 8) :: forref_var_213(4, 12)
END MODULE yoerrta7
MODULE yoerrta8
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(8) :: fracrefa_var_214
  REAL(KIND = 8), DIMENSION(8) :: fracrefb_var_215
  REAL(KIND = 8), DIMENSION(8) :: cfc12_var_216
  REAL(KIND = 8), DIMENSION(8) :: cfc22adj
  REAL(KIND = 8) :: absa_var_217(65, 8)
  REAL(KIND = 8) :: absb_var_218(235, 8)
  REAL(KIND = 8) :: ka_mco2_var_219(19, 8)
  REAL(KIND = 8) :: ka_mn2o_var_220(19, 8)
  REAL(KIND = 8) :: ka_mo3_var_221(19, 8)
  REAL(KIND = 8) :: kb_mco2_var_222(19, 8)
  REAL(KIND = 8) :: kb_mn2o_var_223(19, 8)
  REAL(KIND = 8) :: selfref_var_224(10, 8)
  REAL(KIND = 8) :: forref_var_225(4, 8)
END MODULE yoerrta8
MODULE yoerrta9
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_226(12, 9)
  REAL(KIND = 8), DIMENSION(12) :: fracrefb_var_227
  REAL(KIND = 8) :: absa_var_228(585, 12)
  REAL(KIND = 8) :: absb_var_229(235, 12)
  REAL(KIND = 8) :: ka_mn2o_var_230(9, 19, 12)
  REAL(KIND = 8) :: kb_mn2o_var_231(19, 12)
  REAL(KIND = 8) :: selfref_var_232(10, 12)
  REAL(KIND = 8) :: forref_var_233(4, 12)
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
  REAL(KIND = 8), DIMENSION(59) :: preflog_var_234
  REAL(KIND = 8), DIMENSION(59) :: tref_var_235
  REAL(KIND = 8) :: chi_mls(7, 59)
END MODULE yoerrtrf
MODULE yoerrtwn
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4), DIMENSION(16) :: nspa_var_236
  INTEGER(KIND = 4), DIMENSION(16) :: nspb_var_237
  REAL(KIND = 8), DIMENSION(16) :: delwave
  REAL(KIND = 8), DIMENSION(181, 16) :: totplnk
END MODULE yoerrtwn
MODULE yoesrta16
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_238, strrat1
  INTEGER(KIND = 4) :: layreffr_var_239
  REAL(KIND = 8) :: absa_var_240(585, 16)
  REAL(KIND = 8) :: absb_var_241(235, 16)
  REAL(KIND = 8) :: selfrefc_var_242(10, 16), forrefc_var_243(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_244(16)
END MODULE yoesrta16
MODULE yoesrta17
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_245, strrat_var_246
  INTEGER(KIND = 4) :: layreffr_var_247
  REAL(KIND = 8) :: absa_var_248(585, 16)
  REAL(KIND = 8) :: absb_var_249(1175, 16)
  REAL(KIND = 8) :: selfrefc_var_250(10, 16), forrefc_var_251(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_252(16, 5)
END MODULE yoesrta17
MODULE yoesrta18
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_253, strrat_var_254
  INTEGER(KIND = 4) :: layreffr_var_255
  REAL(KIND = 8) :: absa_var_256(585, 16)
  REAL(KIND = 8) :: absb_var_257(235, 16)
  REAL(KIND = 8) :: selfrefc_var_258(10, 16), forrefc_var_259(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_260(16, 9)
END MODULE yoesrta18
MODULE yoesrta19
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_261, strrat_var_262
  INTEGER(KIND = 4) :: layreffr_var_263
  REAL(KIND = 8) :: absa_var_264(585, 16)
  REAL(KIND = 8) :: absb_var_265(235, 16)
  REAL(KIND = 8) :: selfrefc_var_266(10, 16), forrefc_var_267(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_268(16, 9)
END MODULE yoesrta19
MODULE yoesrta20
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_269
  INTEGER(KIND = 4) :: layreffr_var_270
  REAL(KIND = 8) :: absa_var_271(65, 16)
  REAL(KIND = 8) :: absb_var_272(235, 16)
  REAL(KIND = 8) :: selfrefc_var_273(10, 16), forrefc_var_274(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_275(16), absch4c(16)
END MODULE yoesrta20
MODULE yoesrta21
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_276, strrat_var_277
  INTEGER(KIND = 4) :: layreffr_var_278
  REAL(KIND = 8) :: absa_var_279(585, 16)
  REAL(KIND = 8) :: absb_var_280(1175, 16)
  REAL(KIND = 8) :: selfrefc_var_281(10, 16), forrefc_var_282(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_283(16, 9)
END MODULE yoesrta21
MODULE yoesrta22
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_284, strrat_var_285
  INTEGER(KIND = 4) :: layreffr_var_286
  REAL(KIND = 8) :: absa_var_287(585, 16)
  REAL(KIND = 8) :: absb_var_288(235, 16)
  REAL(KIND = 8) :: selfrefc_var_289(10, 16), forrefc_var_290(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_291(16, 9)
END MODULE yoesrta22
MODULE yoesrta23
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: givfac
  INTEGER(KIND = 4) :: layreffr_var_292
  REAL(KIND = 8) :: absa_var_293(65, 16)
  REAL(KIND = 8) :: selfrefc_var_294(10, 16), forrefc_var_295(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_296(16), raylc_var_297(16)
END MODULE yoesrta23
MODULE yoesrta24
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: strrat_var_298
  INTEGER(KIND = 4) :: layreffr_var_299
  REAL(KIND = 8) :: absa_var_300(585, 16)
  REAL(KIND = 8) :: absb_var_301(235, 16)
  REAL(KIND = 8) :: selfrefc_var_302(10, 16), forrefc_var_303(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_304(16, 9)
  REAL(KIND = 8) :: abso3ac_var_305(16), abso3bc_var_306(16), raylac(16, 9), raylbc(16)
END MODULE yoesrta24
MODULE yoesrta25
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4) :: layreffr_var_307
  REAL(KIND = 8) :: absa_var_308(65, 16)
  REAL(KIND = 8) :: sfluxrefc_var_309(16)
  REAL(KIND = 8) :: raylc_var_310(16), abso3ac_var_311(16), abso3bc_var_312(16)
END MODULE yoesrta25
MODULE yoesrta26
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: sfluxrefc_var_313(16), raylc_var_314(16)
END MODULE yoesrta26
MODULE yoesrta27
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: scalekur
  INTEGER(KIND = 4) :: layreffr_var_315
  REAL(KIND = 8) :: absa_var_316(65, 16)
  REAL(KIND = 8) :: absb_var_317(235, 16)
  REAL(KIND = 8) :: sfluxrefc_var_318(16), raylc_var_319(16)
END MODULE yoesrta27
MODULE yoesrta28
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_320, strrat_var_321
  INTEGER(KIND = 4) :: layreffr_var_322
  REAL(KIND = 8) :: absa_var_323(585, 16)
  REAL(KIND = 8) :: absb_var_324(1175, 16)
  REAL(KIND = 8) :: sfluxrefc_var_325(16, 5)
END MODULE yoesrta28
MODULE yoesrta29
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_326
  INTEGER(KIND = 4) :: layreffr_var_327
  REAL(KIND = 8) :: absa_var_328(65, 16)
  REAL(KIND = 8) :: absb_var_329(235, 16)
  REAL(KIND = 8) :: selfrefc_var_330(10, 16), forrefc_var_331(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_332(16), absh2oc(16), absco2c(16)
END MODULE yoesrta29
MODULE yoesrtwn
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4), DIMENSION(16 : 29) :: nspa_var_333
  INTEGER(KIND = 4), DIMENSION(16 : 29) :: nspb_var_334
  REAL(KIND = 8), DIMENSION(59) :: preflog_var_335
  REAL(KIND = 8), DIMENSION(59) :: tref_var_336
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
  RECURSIVE SUBROUTINE assert_units_gas(this_var_338, iunits, igas, scale_factor, istatus)
    CLASS(gas_type), INTENT(IN) :: this_var_338
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
  ELEMENTAL SUBROUTINE sample_from_pdf(this_var_340, fsd_var_341, cdf_var_342, val_var_343)
    CLASS(pdf_sampler_type), INTENT(IN) :: this_var_340
    REAL(KIND = 8), INTENT(IN) :: fsd_var_341, cdf_var_342
    REAL(KIND = 8), INTENT(OUT) :: val_var_343
    INTEGER :: ifsd_var_344, icdf_var_345
    REAL(KIND = 8) :: wfsd_var_346, wcdf_var_347
    wcdf_var_347 = cdf_var_342 * (this_var_340 % ncdf - 1) + 1.0D0
    icdf_var_345 = MAX(1, MIN(INT(wcdf_var_347), this_var_340 % ncdf - 1))
    wcdf_var_347 = MAX(0.0D0, MIN(wcdf_var_347 - icdf_var_345, 1.0D0))
    wfsd_var_346 = (fsd_var_341 - this_var_340 % fsd1) * this_var_340 % inv_fsd_interval + 1.0D0
    ifsd_var_344 = MAX(1, MIN(INT(wfsd_var_346), this_var_340 % nfsd - 1))
    wfsd_var_346 = MAX(0.0D0, MIN(wfsd_var_346 - ifsd_var_344, 1.0D0))
    val_var_343 = (1.0D0 - wcdf_var_347) * (1.0D0 - wfsd_var_346) * this_var_340 % val(icdf_var_345, ifsd_var_344) + (1.0D0 - wcdf_var_347) * wfsd_var_346 * this_var_340 % val(icdf_var_345, ifsd_var_344 + 1) + wcdf_var_347 * (1.0D0 - wfsd_var_346) * this_var_340 % val(icdf_var_345 + 1, ifsd_var_344) + wcdf_var_347 * wfsd_var_346 * this_var_340 % val(icdf_var_345 + 1, ifsd_var_344 + 1)
  END SUBROUTINE sample_from_pdf
  SUBROUTINE sample_from_pdf_masked_block(this_var_348, nz, ng_var_349, fsd_var_350, cdf_var_351, val_var_352, mask_var_353)
    CLASS(pdf_sampler_type), INTENT(IN) :: this_var_348
    INTEGER, INTENT(IN) :: nz, ng_var_349
    REAL(KIND = 8), INTENT(IN) :: fsd_var_350(nz), cdf_var_351(ng_var_349, nz)
    REAL(KIND = 8), INTENT(OUT) :: val_var_352(:, :)
    LOGICAL, INTENT(IN), OPTIONAL :: mask_var_353(nz)
    INTEGER :: jz, jg_var_354
    INTEGER :: ifsd_var_355, icdf_var_356
    REAL(KIND = 8) :: wfsd_var_357, wcdf_var_358
    DO jz = 1, nz
      IF (mask_var_353(jz)) THEN
        DO jg_var_354 = 1, ng_var_349
          IF (cdf_var_351(jg_var_354, jz) > 0.0D0) THEN
            wcdf_var_358 = cdf_var_351(jg_var_354, jz) * (this_var_348 % ncdf - 1) + 1.0D0
            icdf_var_356 = MAX(1, MIN(INT(wcdf_var_358), this_var_348 % ncdf - 1))
            wcdf_var_358 = MAX(0.0D0, MIN(wcdf_var_358 - icdf_var_356, 1.0D0))
            wfsd_var_357 = (fsd_var_350(jz) - this_var_348 % fsd1) * this_var_348 % inv_fsd_interval + 1.0D0
            ifsd_var_355 = MAX(1, MIN(INT(wfsd_var_357), this_var_348 % nfsd - 1))
            wfsd_var_357 = MAX(0.0D0, MIN(wfsd_var_357 - ifsd_var_355, 1.0D0))
            val_var_352(jg_var_354, jz) = (1.0D0 - wcdf_var_358) * (1.0D0 - wfsd_var_357) * this_var_348 % val(icdf_var_356, ifsd_var_355) + (1.0D0 - wcdf_var_358) * wfsd_var_357 * this_var_348 % val(icdf_var_356, ifsd_var_355 + 1) + wcdf_var_358 * (1.0D0 - wfsd_var_357) * this_var_348 % val(icdf_var_356 + 1, ifsd_var_355) + wcdf_var_358 * wfsd_var_357 * this_var_348 % val(icdf_var_356 + 1, ifsd_var_355 + 1)
          ELSE
            val_var_352(jg_var_354, jz) = 0.0D0
          END IF
        END DO
      END IF
    END DO
  END SUBROUTINE sample_from_pdf_masked_block
END MODULE radiation_pdf_sampler
MODULE radiation_cloud_generator
  CONTAINS
  SUBROUTINE cloud_generator(ng_var_359, nlev_var_360, i_overlap_scheme, iseed_var_361, frac_threshold_var_362, frac_var_363, overlap_param_var_364, decorrelation_scaling, fractional_std_var_365, pdf_sampler_var_366, od_scaling_var_367, total_cloud_cover_var_368, use_beta_overlap, use_vectorizable_generator)
    USE radiation_pdf_sampler, ONLY: pdf_sampler_type
    USE random_numbers_mix, ONLY: initialize_random_numbers, randomnumberstream, uniform_distribution
    USE radiation_cloud_cover, ONLY: cum_cloud_cover_exp_ran
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng_var_359
    INTEGER, INTENT(IN) :: nlev_var_360
    INTEGER, INTENT(IN) :: i_overlap_scheme
    INTEGER, INTENT(IN) :: iseed_var_361
    REAL(KIND = 8), INTENT(IN) :: frac_threshold_var_362
    REAL(KIND = 8), INTENT(IN) :: frac_var_363(nlev_var_360)
    REAL(KIND = 8), INTENT(IN) :: overlap_param_var_364(nlev_var_360 - 1)
    REAL(KIND = 8), INTENT(IN) :: decorrelation_scaling
    REAL(KIND = 8), INTENT(IN) :: fractional_std_var_365(nlev_var_360)
    TYPE(pdf_sampler_type), INTENT(IN) :: pdf_sampler_var_366
    LOGICAL, INTENT(IN), OPTIONAL :: use_beta_overlap
    LOGICAL, INTENT(IN), OPTIONAL :: use_vectorizable_generator
    REAL(KIND = 8), INTENT(OUT) :: od_scaling_var_367(ng_var_359, nlev_var_360)
    REAL(KIND = 8), INTENT(OUT) :: total_cloud_cover_var_368
    REAL(KIND = 8) :: cum_cloud_cover_var_369(nlev_var_360)
    REAL(KIND = 8) :: trigger_var_370
    REAL(KIND = 8) :: rand_top(ng_var_359)
    REAL(KIND = 8) :: overlap_param_inhom_var_371(nlev_var_360 - 1)
    TYPE(randomnumberstream) :: random_stream_var_372
    INTEGER :: ibegin_var_373, iend_var_374
    INTEGER :: itrigger_var_375
    INTEGER :: jlev_var_376, jg_var_377
    REAL(KIND = 8), DIMENSION(nlev_var_360 - 1) :: pair_cloud_cover_var_378, overhang_var_379
    LOGICAL :: use_vec_gen
    CALL cum_cloud_cover_exp_ran(nlev_var_360, frac_var_363, overlap_param_var_364, cum_cloud_cover_var_369, pair_cloud_cover_var_378, .FALSE.)
    total_cloud_cover_var_368 = cum_cloud_cover_var_369(nlev_var_360)
    DO jlev_var_376 = 1, nlev_var_360 - 1
      overhang_var_379(jlev_var_376) = cum_cloud_cover_var_369(jlev_var_376 + 1) - cum_cloud_cover_var_369(jlev_var_376)
    END DO
    IF (total_cloud_cover_var_368 < 1D-06) THEN
      total_cloud_cover_var_368 = 0.0D0
    ELSE
      jlev_var_376 = 1
      DO WHILE (frac_var_363(jlev_var_376) <= 0.0D0)
        jlev_var_376 = jlev_var_376 + 1
      END DO
      ibegin_var_373 = jlev_var_376
      iend_var_374 = jlev_var_376
      DO jlev_var_376 = jlev_var_376 + 1, nlev_var_360
        IF (frac_var_363(jlev_var_376) > 0.0D0) THEN
          iend_var_374 = jlev_var_376
        END IF
      END DO
      overlap_param_inhom_var_371 = overlap_param_var_364
      DO jlev_var_376 = ibegin_var_373, iend_var_374 - 1
        IF (overlap_param_var_364(jlev_var_376) > 0.0D0) THEN
          overlap_param_inhom_var_371(jlev_var_376) = overlap_param_var_364(jlev_var_376) ** (2.0D0)
        END IF
      END DO
      od_scaling_var_367 = 0.0D0
      use_vec_gen = .FALSE.
      IF (.NOT. use_vec_gen) THEN
        CALL initialize_random_numbers(iseed_var_361, random_stream_var_372)
        CALL uniform_distribution(rand_top, random_stream_var_372)
        DO jg_var_377 = 1, ng_var_359
          trigger_var_370 = rand_top(jg_var_377) * total_cloud_cover_var_368
          jlev_var_376 = ibegin_var_373
          DO WHILE (trigger_var_370 > cum_cloud_cover_var_369(jlev_var_376) .AND. jlev_var_376 < iend_var_374)
            jlev_var_376 = jlev_var_376 + 1
          END DO
          itrigger_var_375 = jlev_var_376
          CALL generate_column_exp_ran(ng_var_359, nlev_var_360, jg_var_377, random_stream_var_372, pdf_sampler_var_366, frac_var_363, pair_cloud_cover_var_378, cum_cloud_cover_var_369, overhang_var_379, fractional_std_var_365, overlap_param_inhom_var_371, itrigger_var_375, iend_var_374, od_scaling_var_367)
        END DO
      ELSE
        CALL generate_columns_exp_ran(ng_var_359, nlev_var_360, iseed_var_361, pdf_sampler_var_366, total_cloud_cover_var_368, 1D-06, frac_var_363, pair_cloud_cover_var_378, cum_cloud_cover_var_369, overhang_var_379, fractional_std_var_365, overlap_param_inhom_var_371, ibegin_var_373, iend_var_374, od_scaling_var_367)
      END IF
    END IF
  END SUBROUTINE cloud_generator
  SUBROUTINE generate_column_exp_ran(ng_var_380, nlev_var_382, ig_var_381, random_stream_var_383, pdf_sampler_var_384, frac_var_385, pair_cloud_cover_var_388, cum_cloud_cover_var_386, overhang_var_389, fractional_std_var_387, overlap_param_inhom_var_390, itrigger_var_391, iend_var_392, od_scaling_var_393)
    USE random_numbers_mix, ONLY: randomnumberstream, uniform_distribution
    USE radiation_pdf_sampler, ONLY: pdf_sampler_type, sample_from_pdf
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng_var_380, ig_var_381
    INTEGER, INTENT(IN) :: nlev_var_382
    TYPE(randomnumberstream), INTENT(INOUT) :: random_stream_var_383
    TYPE(pdf_sampler_type), INTENT(IN) :: pdf_sampler_var_384
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_382) :: frac_var_385, cum_cloud_cover_var_386, fractional_std_var_387
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_382 - 1) :: pair_cloud_cover_var_388, overhang_var_389
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_382 - 1) :: overlap_param_inhom_var_390
    INTEGER, INTENT(IN) :: itrigger_var_391, iend_var_392
    REAL(KIND = 8), INTENT(INOUT), DIMENSION(ng_var_380, nlev_var_382) :: od_scaling_var_393
    INTEGER :: jlev_var_394, jcloud
    INTEGER :: n_layers_to_scale
    INTEGER :: iy
    LOGICAL :: do_fill_od_scaling
    REAL(KIND = 8) :: rand_cloud_var_395(nlev_var_382)
    REAL(KIND = 8) :: rand_inhom1(nlev_var_382), rand_inhom2_var_396(nlev_var_382)
    n_layers_to_scale = 1
    iy = 0
    CALL uniform_distribution(rand_cloud_var_395(1 : (iend_var_392 + 1 - itrigger_var_391)), random_stream_var_383)
    DO jlev_var_394 = itrigger_var_391 + 1, iend_var_392 + 1
      do_fill_od_scaling = .FALSE.
      IF (jlev_var_394 <= iend_var_392) THEN
        iy = iy + 1
        IF (n_layers_to_scale > 0) THEN
          IF (rand_cloud_var_395(iy) * frac_var_385(jlev_var_394 - 1) < frac_var_385(jlev_var_394) + frac_var_385(jlev_var_394 - 1) - pair_cloud_cover_var_388(jlev_var_394 - 1)) THEN
            n_layers_to_scale = n_layers_to_scale + 1
          ELSE
            do_fill_od_scaling = .TRUE.
          END IF
        ELSE
          IF (rand_cloud_var_395(iy) * (cum_cloud_cover_var_386(jlev_var_394 - 1) - frac_var_385(jlev_var_394 - 1)) < pair_cloud_cover_var_388(jlev_var_394 - 1) - overhang_var_389(jlev_var_394 - 1) - frac_var_385(jlev_var_394 - 1)) THEN
            n_layers_to_scale = 1
          END IF
        END IF
      ELSE
        do_fill_od_scaling = .TRUE.
      END IF
      IF (do_fill_od_scaling) THEN
        CALL uniform_distribution(rand_inhom1(1 : n_layers_to_scale), random_stream_var_383)
        CALL uniform_distribution(rand_inhom2_var_396(1 : n_layers_to_scale), random_stream_var_383)
        DO jcloud = 2, n_layers_to_scale
          IF (rand_inhom2_var_396(jcloud) < overlap_param_inhom_var_390(jlev_var_394 - n_layers_to_scale + jcloud - 2)) THEN
            rand_inhom1(jcloud) = rand_inhom1(jcloud - 1)
          END IF
        END DO
        CALL sample_from_pdf(pdf_sampler_var_384, fractional_std_var_387(jlev_var_394 - n_layers_to_scale : jlev_var_394 - 1), rand_inhom1(1 : n_layers_to_scale), od_scaling_var_393(ig_var_381, jlev_var_394 - n_layers_to_scale : jlev_var_394 - 1))
        n_layers_to_scale = 0
      END IF
    END DO
  END SUBROUTINE generate_column_exp_ran
  SUBROUTINE generate_columns_exp_ran(ng_var_397, nlev_var_398, iseed_var_399, pdf_sampler_var_400, total_cloud_cover_var_401, frac_threshold_var_402, frac_var_403, pair_cloud_cover_var_406, cum_cloud_cover_var_404, overhang_var_407, fractional_std_var_405, overlap_param_inhom_var_408, ibegin_var_409, iend_var_410, od_scaling_var_411)
    USE radiation_random_numbers, ONLY: initialize, rng_type, uniform_distribution_1d, uniform_distribution_2d, uniform_distribution_2d_masked
    USE radiation_pdf_sampler, ONLY: pdf_sampler_type, sample_from_pdf_masked_block
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng_var_397
    INTEGER, INTENT(IN) :: nlev_var_398
    INTEGER, INTENT(IN) :: iseed_var_399
    TYPE(rng_type) :: random_number_generator
    TYPE(pdf_sampler_type), INTENT(IN) :: pdf_sampler_var_400
    REAL(KIND = 8), INTENT(IN) :: total_cloud_cover_var_401
    REAL(KIND = 8), INTENT(IN) :: frac_threshold_var_402
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_398) :: frac_var_403, cum_cloud_cover_var_404, fractional_std_var_405
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_398 - 1) :: pair_cloud_cover_var_406, overhang_var_407
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_398 - 1) :: overlap_param_inhom_var_408
    INTEGER, INTENT(INOUT) :: ibegin_var_409, iend_var_410
    REAL(KIND = 8), INTENT(INOUT), DIMENSION(ng_var_397, nlev_var_398) :: od_scaling_var_411
    INTEGER :: jlev_var_412, jg_var_413
    REAL(KIND = 8) :: rand_cloud_var_414(ng_var_397, ibegin_var_409 : iend_var_410)
    REAL(KIND = 8) :: rand_inhom(ng_var_397, ibegin_var_409 - 1 : iend_var_410), rand_inhom2_var_415(ng_var_397, ibegin_var_409 : iend_var_410)
    LOGICAL :: is_any_cloud(ibegin_var_409 : iend_var_410)
    REAL(KIND = 8) :: trigger_var_416(ng_var_397)
    LOGICAL :: is_cloud(ng_var_397)
    LOGICAL :: prev_cloud(ng_var_397)
    LOGICAL :: first_cloud(ng_var_397)
    LOGICAL :: found_cloud(ng_var_397)
    is_any_cloud = (frac_var_403(ibegin_var_409 : iend_var_410) >= frac_threshold_var_402)
    CALL initialize(random_number_generator, 1, iseed_var_90 = iseed_var_399, nmaxstreams_var_91 = ng_var_397)
    CALL uniform_distribution_1d(random_number_generator, trigger_var_416)
    CALL uniform_distribution_2d_masked(random_number_generator, rand_cloud_var_414, is_any_cloud)
    CALL uniform_distribution_2d(random_number_generator, rand_inhom)
    CALL uniform_distribution_2d_masked(random_number_generator, rand_inhom2_var_415, is_any_cloud)
    trigger_var_416 = trigger_var_416 * total_cloud_cover_var_401
    found_cloud = .FALSE.
    is_cloud = .FALSE.
    first_cloud = .FALSE.
    DO jlev_var_412 = ibegin_var_409, iend_var_410
      IF (is_any_cloud(jlev_var_412)) THEN
        DO jg_var_413 = 1, ng_var_397
          prev_cloud(jg_var_413) = is_cloud(jg_var_413)
          first_cloud(jg_var_413) = (trigger_var_416(jg_var_413) <= cum_cloud_cover_var_404(jlev_var_412) .AND. .NOT. found_cloud(jg_var_413))
          found_cloud(jg_var_413) = found_cloud(jg_var_413) .OR. first_cloud(jg_var_413)
          is_cloud(jg_var_413) = first_cloud(jg_var_413) .OR. found_cloud(jg_var_413) .AND. MERGE(rand_cloud_var_414(jg_var_413, jlev_var_412) * frac_var_403(jlev_var_412 - 1) < frac_var_403(jlev_var_412) + frac_var_403(jlev_var_412 - 1) - pair_cloud_cover_var_406(jlev_var_412 - 1), rand_cloud_var_414(jg_var_413, jlev_var_412) * (cum_cloud_cover_var_404(jlev_var_412 - 1) - frac_var_403(jlev_var_412 - 1)) < pair_cloud_cover_var_406(jlev_var_412 - 1) - overhang_var_407(jlev_var_412 - 1) - frac_var_403(jlev_var_412 - 1), prev_cloud(jg_var_413))
          rand_inhom(jg_var_413, jlev_var_412) = MERGE(MERGE(rand_inhom(jg_var_413, jlev_var_412 - 1), rand_inhom(jg_var_413, jlev_var_412), rand_inhom2_var_415(jg_var_413, jlev_var_412) < overlap_param_inhom_var_408(jlev_var_412 - 1) .AND. prev_cloud(jg_var_413)), 0.0D0, is_cloud(jg_var_413))
        END DO
      ELSE
        is_cloud = .FALSE.
      END IF
    END DO
    CALL sample_from_pdf_masked_block(pdf_sampler_var_400, iend_var_410 - ibegin_var_409 + 1, ng_var_397, fractional_std_var_405(ibegin_var_409 : iend_var_410), rand_inhom(:, ibegin_var_409 : iend_var_410), od_scaling_var_411(:, ibegin_var_409 : iend_var_410), is_any_cloud)
  END SUBROUTINE generate_columns_exp_ran
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
  SUBROUTINE calc_surface_spectral(this_var_423, config, istartcol_var_424, iendcol_var_425)
    USE radiation_config, ONLY: config_type
    CLASS(flux_type), INTENT(INOUT) :: this_var_423
    TYPE(config_type), INTENT(IN) :: config
    INTEGER, INTENT(IN) :: istartcol_var_424, iendcol_var_425
    INTEGER :: jcol_var_426
    DO jcol_var_426 = istartcol_var_424, iendcol_var_425
      CALL indexed_sum(this_var_423 % sw_dn_direct_surf_g(:, jcol_var_426), config % i_band_from_reordered_g_sw, this_var_423 % sw_dn_direct_surf_band(:, jcol_var_426))
      CALL indexed_sum(this_var_423 % sw_dn_diffuse_surf_g(:, jcol_var_426), config % i_band_from_reordered_g_sw, this_var_423 % sw_dn_surf_band(:, jcol_var_426))
      this_var_423 % sw_dn_surf_band(:, jcol_var_426) = this_var_423 % sw_dn_surf_band(:, jcol_var_426) + this_var_423 % sw_dn_direct_surf_band(:, jcol_var_426)
    END DO
    DO jcol_var_426 = istartcol_var_424, iendcol_var_425
      CALL indexed_sum(this_var_423 % sw_dn_direct_surf_clear_g(:, jcol_var_426), config % i_band_from_reordered_g_sw, this_var_423 % sw_dn_direct_surf_clear_band(:, jcol_var_426))
      CALL indexed_sum(this_var_423 % sw_dn_diffuse_surf_clear_g(:, jcol_var_426), config % i_band_from_reordered_g_sw, this_var_423 % sw_dn_surf_clear_band(:, jcol_var_426))
      this_var_423 % sw_dn_surf_clear_band(:, jcol_var_426) = this_var_423 % sw_dn_surf_clear_band(:, jcol_var_426) + this_var_423 % sw_dn_direct_surf_clear_band(:, jcol_var_426)
    END DO
  END SUBROUTINE calc_surface_spectral
  SUBROUTINE calc_toa_spectral(this_var_427, config, istartcol_var_428, iendcol_var_429)
    USE radiation_config, ONLY: config_type
    CLASS(flux_type), INTENT(INOUT) :: this_var_427
    TYPE(config_type), INTENT(IN) :: config
    INTEGER, INTENT(IN) :: istartcol_var_428, iendcol_var_429
  END SUBROUTINE calc_toa_spectral
  PURE SUBROUTINE indexed_sum(source_var_430, ind_var_431, dest)
    REAL(KIND = 8), INTENT(IN) :: source_var_430(:)
    INTEGER, INTENT(IN) :: ind_var_431(:)
    REAL(KIND = 8), INTENT(OUT) :: dest(:)
    INTEGER :: ig_var_432, jg_var_433, istart, iend_var_434
    dest = 0.0
    istart = LBOUND(source_var_430, 1)
    iend_var_434 = UBOUND(source_var_430, 1)
    DO jg_var_433 = istart, iend_var_434
      ig_var_432 = ind_var_431(jg_var_433)
      dest(ig_var_432) = dest(ig_var_432) + source_var_430(jg_var_433)
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
  SUBROUTINE get_albedos(this_var_438, istartcol_var_439, iendcol_var_440, config, sw_albedo_direct_var_442, sw_albedo_diffuse_var_443, lw_albedo_var_441)
    USE radiation_config, ONLY: config_type
    USE yomlun_ifsaux, ONLY: nulerr
    USE radiation_io, ONLY: radiation_abort
    CLASS(single_level_type), INTENT(IN) :: this_var_438
    TYPE(config_type), INTENT(IN) :: config
    INTEGER, INTENT(IN) :: istartcol_var_439, iendcol_var_440
    REAL(KIND = 8), INTENT(OUT), OPTIONAL :: lw_albedo_var_441(140, istartcol_var_439 : iendcol_var_440)
    REAL(KIND = 8), INTENT(OUT), DIMENSION(112, istartcol_var_439 : iendcol_var_440) :: sw_albedo_direct_var_442, sw_albedo_diffuse_var_443
    REAL(KIND = 8) :: sw_albedo_band(istartcol_var_439 : iendcol_var_440, 14)
    INTEGER :: nalbedoband
    INTEGER :: jband_var_444, jalbedoband, jcol_var_445
    nalbedoband = SIZE(config % sw_albedo_weights, 1)
    IF (SIZE(this_var_438 % sw_albedo, 2) /= nalbedoband) THEN
      WRITE(nulerr, '(a,i0,a)') '*** error: single_level%sw_albedo does not have the expected ', nalbedoband, ' bands'
      CALL radiation_abort
    END IF
    DO jband_var_444 = 1, 14
      DO jcol_var_445 = istartcol_var_439, iendcol_var_440
        sw_albedo_band(jcol_var_445, jband_var_444) = 0.0D0
      END DO
    END DO
    DO jband_var_444 = 1, 14
      DO jalbedoband = 1, nalbedoband
        IF (config % sw_albedo_weights(jalbedoband, jband_var_444) /= 0.0D0) THEN
          DO jcol_var_445 = istartcol_var_439, iendcol_var_440
            sw_albedo_band(jcol_var_445, jband_var_444) = sw_albedo_band(jcol_var_445, jband_var_444) + config % sw_albedo_weights(jalbedoband, jband_var_444) * this_var_438 % sw_albedo(jcol_var_445, jalbedoband)
          END DO
        END IF
      END DO
    END DO
    sw_albedo_diffuse_var_443 = TRANSPOSE(sw_albedo_band(istartcol_var_439 : iendcol_var_440, config % i_band_from_reordered_g_sw))
    DO jband_var_444 = 1, 14
      DO jcol_var_445 = istartcol_var_439, iendcol_var_440
        sw_albedo_band(jcol_var_445, jband_var_444) = 0.0D0
      END DO
    END DO
    DO jband_var_444 = 1, 14
      DO jalbedoband = 1, nalbedoband
        IF (config % sw_albedo_weights(jalbedoband, jband_var_444) /= 0.0D0) THEN
          DO jcol_var_445 = istartcol_var_439, iendcol_var_440
            sw_albedo_band(jcol_var_445, jband_var_444) = sw_albedo_band(jcol_var_445, jband_var_444) + config % sw_albedo_weights(jalbedoband, jband_var_444) * this_var_438 % sw_albedo_direct(jcol_var_445, jalbedoband)
          END DO
        END IF
      END DO
    END DO
    sw_albedo_direct_var_442 = TRANSPOSE(sw_albedo_band(istartcol_var_439 : iendcol_var_440, config % i_band_from_reordered_g_sw))
    IF (MAXVAL(config % i_emiss_from_band_lw) > SIZE(this_var_438 % lw_emissivity, 2)) THEN
      WRITE(nulerr, '(a,i0,a)') '*** error: single_level%lw_emissivity has fewer than required ', MAXVAL(config % i_emiss_from_band_lw), ' bands'
      CALL radiation_abort
    END IF
    lw_albedo_var_441 = 1.0D0 - TRANSPOSE(this_var_438 % lw_emissivity(istartcol_var_439 : iendcol_var_440, config % i_emiss_from_band_lw(config % i_band_from_reordered_g_lw)))
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
  ELEMENTAL SUBROUTINE delta_eddington_extensive(od_var_448, scat_od_var_449, scat_od_g)
    REAL(KIND = 8), INTENT(INOUT) :: od_var_448, scat_od_var_449, scat_od_g
    REAL(KIND = 8) :: f_var_450, g_var_451
    IF (scat_od_var_449 > 0.0D0) THEN
      g_var_451 = scat_od_g / scat_od_var_449
    ELSE
      g_var_451 = 0.0
    END IF
    f_var_450 = g_var_451 * g_var_451
    od_var_448 = od_var_448 - scat_od_var_449 * f_var_450
    scat_od_var_449 = scat_od_var_449 * (1.0D0 - f_var_450)
    scat_od_g = scat_od_var_449 * g_var_451 / (1.0D0 + g_var_451)
  END SUBROUTINE delta_eddington_extensive
  SUBROUTINE add_aerosol_optics(nlev_var_452, istartcol_var_453, iendcol_var_454, config, thermodynamics, gas, aerosol, od_lw_var_455, ssa_lw_var_456, g_lw_var_457, od_sw_var_458, ssa_sw_var_459, g_sw_var_460)
    USE radiation_config, ONLY: config_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE radiation_gas, ONLY: gas_type
    USE radiation_aerosol, ONLY: aerosol_type
    INTEGER, INTENT(IN) :: nlev_var_452
    INTEGER, INTENT(IN) :: istartcol_var_453, iendcol_var_454
    TYPE(config_type), INTENT(IN), TARGET :: config
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics
    TYPE(gas_type), INTENT(IN) :: gas
    TYPE(aerosol_type), INTENT(IN) :: aerosol
    REAL(KIND = 8), DIMENSION(140, nlev_var_452, istartcol_var_453 : iendcol_var_454), INTENT(INOUT) :: od_lw_var_455
    REAL(KIND = 8), DIMENSION(0, nlev_var_452, istartcol_var_453 : iendcol_var_454), INTENT(OUT) :: ssa_lw_var_456, g_lw_var_457
    REAL(KIND = 8), DIMENSION(112, nlev_var_452, istartcol_var_453 : iendcol_var_454), INTENT(INOUT) :: od_sw_var_458, ssa_sw_var_459
    REAL(KIND = 8), DIMENSION(112, nlev_var_452, istartcol_var_453 : iendcol_var_454), INTENT(OUT) :: g_sw_var_460
    CALL add_aerosol_optics_direct(nlev_var_452, istartcol_var_453, iendcol_var_454, config, aerosol, od_lw_var_455, ssa_lw_var_456, g_lw_var_457, od_sw_var_458, ssa_sw_var_459, g_sw_var_460)
  END SUBROUTINE add_aerosol_optics
  SUBROUTINE add_aerosol_optics_direct(nlev_var_461, istartcol_var_462, iendcol_var_463, config, aerosol, od_lw_var_464, ssa_lw_var_465, g_lw_var_466, od_sw_var_467, ssa_sw_var_468, g_sw_var_469)
    USE radiation_config, ONLY: config_type
    USE radiation_aerosol, ONLY: aerosol_type
    USE yomlun_ifsaux, ONLY: nulerr
    USE radiation_io, ONLY: radiation_abort
    INTEGER, INTENT(IN) :: nlev_var_461
    INTEGER, INTENT(IN) :: istartcol_var_462, iendcol_var_463
    TYPE(config_type), INTENT(IN), TARGET :: config
    TYPE(aerosol_type), INTENT(IN) :: aerosol
    REAL(KIND = 8), DIMENSION(140, nlev_var_461, istartcol_var_462 : iendcol_var_463), INTENT(INOUT) :: od_lw_var_464
    REAL(KIND = 8), DIMENSION(0, nlev_var_461, istartcol_var_462 : iendcol_var_463), INTENT(OUT) :: ssa_lw_var_465, g_lw_var_466
    REAL(KIND = 8), DIMENSION(112, nlev_var_461, istartcol_var_462 : iendcol_var_463), INTENT(INOUT) :: od_sw_var_467, ssa_sw_var_468
    REAL(KIND = 8), DIMENSION(112, nlev_var_461, istartcol_var_462 : iendcol_var_463), INTENT(OUT) :: g_sw_var_469
    REAL(KIND = 8) :: local_od, local_scat
    REAL(KIND = 8), DIMENSION(14, nlev_var_461) :: od_sw_aerosol, scat_sw_aerosol, scat_g_sw_aerosol
    REAL(KIND = 8), DIMENSION(16, nlev_var_461) :: od_lw_aerosol
    INTEGER :: jcol_var_470, jlev_var_471, jg_var_472, jb_var_473
    INTEGER :: istartlev_var_474, iendlev_var_475
    INTEGER :: iband_var_476
    IF (UBOUND(aerosol % od_sw, 1) /= 14) THEN
      WRITE(nulerr, '(a,i0,a,i0)') '*** error: aerosol%od_sw contains ', UBOUND(aerosol % od_sw, 1), ' band, expected ', 14
      CALL radiation_abort
    END IF
    istartlev_var_474 = LBOUND(aerosol % od_sw, 2)
    iendlev_var_475 = UBOUND(aerosol % od_sw, 2)
    DO jcol_var_470 = istartcol_var_462, iendcol_var_463
      DO jlev_var_471 = 1, nlev_var_461
        DO jg_var_472 = 1, 112
          g_sw_var_469(jg_var_472, jlev_var_471, jcol_var_470) = 0.0D0
        END DO
      END DO
    END DO
    DO jcol_var_470 = istartcol_var_462, iendcol_var_463
      DO jlev_var_471 = istartlev_var_474, iendlev_var_475
        DO jb_var_473 = 1, 14
          od_sw_aerosol(jb_var_473, jlev_var_471) = aerosol % od_sw(jb_var_473, jlev_var_471, jcol_var_470)
          scat_sw_aerosol(jb_var_473, jlev_var_471) = aerosol % ssa_sw(jb_var_473, jlev_var_471, jcol_var_470) * od_sw_aerosol(jb_var_473, jlev_var_471)
          scat_g_sw_aerosol(jb_var_473, jlev_var_471) = aerosol % g_sw(jb_var_473, jlev_var_471, jcol_var_470) * scat_sw_aerosol(jb_var_473, jlev_var_471)
          CALL delta_eddington_extensive(od_sw_aerosol(jb_var_473, jlev_var_471), scat_sw_aerosol(jb_var_473, jlev_var_471), scat_g_sw_aerosol(jb_var_473, jlev_var_471))
        END DO
      END DO
      DO jlev_var_471 = istartlev_var_474, iendlev_var_475
        IF (od_sw_aerosol(1, jlev_var_471) > 0.0D0) THEN
          DO jg_var_472 = 1, 112
            iband_var_476 = config % i_band_from_reordered_g_sw(jg_var_472)
            local_od = od_sw_var_467(jg_var_472, jlev_var_471, jcol_var_470) + od_sw_aerosol(iband_var_476, jlev_var_471)
            local_scat = ssa_sw_var_468(jg_var_472, jlev_var_471, jcol_var_470) * od_sw_var_467(jg_var_472, jlev_var_471, jcol_var_470) + scat_sw_aerosol(iband_var_476, jlev_var_471)
            g_sw_var_469(jg_var_472, jlev_var_471, jcol_var_470) = scat_g_sw_aerosol(iband_var_476, jlev_var_471) / local_scat
            local_od = od_sw_var_467(jg_var_472, jlev_var_471, jcol_var_470) + od_sw_aerosol(iband_var_476, jlev_var_471)
            ssa_sw_var_468(jg_var_472, jlev_var_471, jcol_var_470) = local_scat / local_od
            od_sw_var_467(jg_var_472, jlev_var_471, jcol_var_470) = local_od
          END DO
        END IF
      END DO
    END DO
    IF (UBOUND(aerosol % od_lw, 1) /= 16) THEN
      WRITE(nulerr, '(a,i0,a,i0)') '*** error: aerosol%od_lw contains ', UBOUND(aerosol % od_lw, 1), ' band, expected ', 16
      CALL radiation_abort
    END IF
    istartlev_var_474 = LBOUND(aerosol % od_lw, 2)
    iendlev_var_475 = UBOUND(aerosol % od_lw, 2)
    DO jcol_var_470 = istartcol_var_462, iendcol_var_463
      DO jlev_var_471 = istartlev_var_474, iendlev_var_475
        DO jb_var_473 = 1, 16
          od_lw_aerosol(jb_var_473, jlev_var_471) = aerosol % od_lw(jb_var_473, jlev_var_471, jcol_var_470) * (1.0D0 - aerosol % ssa_lw(jb_var_473, jlev_var_471, jcol_var_470))
        END DO
      END DO
      DO jlev_var_471 = istartlev_var_474, iendlev_var_475
        DO jg_var_472 = 1, 140
          od_lw_var_464(jg_var_472, jlev_var_471, jcol_var_470) = od_lw_var_464(jg_var_472, jlev_var_471, jcol_var_470) + od_lw_aerosol(config % i_band_from_reordered_g_lw(jg_var_472), jlev_var_471)
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
  SUBROUTINE crop_cloud_fraction(this_var_480, istartcol_var_481, iendcol_var_482, cloud_fraction_threshold, cloud_mixing_ratio_threshold)
    CLASS(cloud_type), INTENT(INOUT) :: this_var_480
    INTEGER, INTENT(IN) :: istartcol_var_481, iendcol_var_482
    INTEGER :: nlev_var_483
    INTEGER :: jcol_var_484, jlev_var_485, jh
    REAL(KIND = 8) :: cloud_fraction_threshold, cloud_mixing_ratio_threshold
    REAL(KIND = 8) :: sum_mixing_ratio(istartcol_var_481 : iendcol_var_482)
    nlev_var_483 = SIZE(this_var_480 % fraction, 2)
    DO jlev_var_485 = 1, nlev_var_483
      DO jcol_var_484 = istartcol_var_481, iendcol_var_482
        sum_mixing_ratio(jcol_var_484) = 0.0D0
      END DO
      DO jh = 1, 2
        DO jcol_var_484 = istartcol_var_481, iendcol_var_482
          sum_mixing_ratio(jcol_var_484) = sum_mixing_ratio(jcol_var_484) + this_var_480 % mixing_ratio(jcol_var_484, jlev_var_485, jh)
        END DO
      END DO
      DO jcol_var_484 = istartcol_var_481, iendcol_var_482
        IF (this_var_480 % fraction(jcol_var_484, jlev_var_485) < cloud_fraction_threshold .OR. sum_mixing_ratio(jcol_var_484) < cloud_mixing_ratio_threshold) THEN
          this_var_480 % fraction(jcol_var_484, jlev_var_485) = 0.0D0
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
  ELEMENTAL SUBROUTINE delta_eddington_scat_od(od_var_486, scat_od_var_487, g_var_488)
    REAL(KIND = 8), INTENT(INOUT) :: od_var_486, scat_od_var_487, g_var_488
    REAL(KIND = 8) :: f_var_489
    f_var_489 = g_var_488 * g_var_488
    od_var_486 = od_var_486 - scat_od_var_487 * f_var_489
    scat_od_var_487 = scat_od_var_487 * (1.0D0 - f_var_489)
    g_var_488 = g_var_488 / (1.0D0 + g_var_488)
  END SUBROUTINE delta_eddington_scat_od
  SUBROUTINE cloud_optics_fn_500(nlev_var_490, istartcol_var_491, iendcol_var_492, config, thermodynamics, cloud, od_lw_cloud_var_493, ssa_lw_cloud_var_494, g_lw_cloud_var_495, od_sw_cloud_var_496, ssa_sw_cloud_var_497, g_sw_cloud_var_498)
    USE radiation_config, ONLY: config_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE radiation_cloud, ONLY: cloud_type
    USE radiation_liquid_optics_socrates, ONLY: calc_liq_optics_socrates
    USE radiation_ice_optics_fu, ONLY: calc_ice_optics_fu_lw, calc_ice_optics_fu_sw
    USE radiation_ice_optics_yi, ONLY: calc_ice_optics_yi_lw, calc_ice_optics_yi_sw
    INTEGER, INTENT(IN) :: nlev_var_490
    INTEGER, INTENT(IN) :: istartcol_var_491, iendcol_var_492
    TYPE(config_type), INTENT(IN), TARGET :: config
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics
    TYPE(cloud_type), INTENT(IN) :: cloud
    REAL(KIND = 8), DIMENSION(16, nlev_var_490, istartcol_var_491 : iendcol_var_492), INTENT(OUT) :: od_lw_cloud_var_493
    REAL(KIND = 8), DIMENSION(0, nlev_var_490, istartcol_var_491 : iendcol_var_492), INTENT(OUT) :: ssa_lw_cloud_var_494, g_lw_cloud_var_495
    REAL(KIND = 8), DIMENSION(14, nlev_var_490, istartcol_var_491 : iendcol_var_492), INTENT(OUT) :: od_sw_cloud_var_496, ssa_sw_cloud_var_497, g_sw_cloud_var_498
    REAL(KIND = 8), DIMENSION(16) :: od_lw_liq, scat_od_lw_liq, g_lw_liq, od_lw_ice, scat_od_lw_ice, g_lw_ice
    REAL(KIND = 8), DIMENSION(14) :: od_sw_liq, scat_od_sw_liq, g_sw_liq, od_sw_ice, scat_od_sw_ice, g_sw_ice
    REAL(KIND = 8) :: lwp_in_cloud, iwp_in_cloud
    REAL(KIND = 8) :: factor_var_499
    INTEGER :: jcol_var_500, jlev_var_501, jb_var_502
    DO jcol_var_500 = istartcol_var_491, iendcol_var_492
      DO jlev_var_501 = 1, nlev_var_490
        DO jb_var_502 = 1, 14
          od_sw_cloud_var_496(jb_var_502, jlev_var_501, jcol_var_500) = 0.0D0
          ssa_sw_cloud_var_497(jb_var_502, jlev_var_501, jcol_var_500) = 0.0D0
          g_sw_cloud_var_498(jb_var_502, jlev_var_501, jcol_var_500) = 0.0D0
        END DO
        DO jb_var_502 = 1, 16
          od_lw_cloud_var_493(jb_var_502, jlev_var_501, jcol_var_500) = 0.0D0
        END DO
        DO jb_var_502 = 1, 0
          ssa_lw_cloud_var_494(jb_var_502, jlev_var_501, jcol_var_500) = 0.0D0
          g_lw_cloud_var_495(jb_var_502, jlev_var_501, jcol_var_500) = 0.0D0
        END DO
      END DO
    END DO
    DO jlev_var_501 = 1, nlev_var_490
      DO jcol_var_500 = istartcol_var_491, iendcol_var_492
        IF (cloud % fraction(jcol_var_500, jlev_var_501) > 0.0D0) THEN
          factor_var_499 = (thermodynamics % pressure_hl(jcol_var_500, jlev_var_501 + 1) - thermodynamics % pressure_hl(jcol_var_500, jlev_var_501)) / (9.80665D0 * cloud % fraction(jcol_var_500, jlev_var_501))
          lwp_in_cloud = factor_var_499 * cloud % q_liq(jcol_var_500, jlev_var_501)
          iwp_in_cloud = factor_var_499 * cloud % q_ice(jcol_var_500, jlev_var_501)
          IF (lwp_in_cloud > 0.0D0) THEN
            CALL calc_liq_optics_socrates(16, config % cloud_optics % liq_coeff_lw, lwp_in_cloud, cloud % re_liq(jcol_var_500, jlev_var_501), od_lw_liq, scat_od_lw_liq, g_lw_liq)
            CALL calc_liq_optics_socrates(14, config % cloud_optics % liq_coeff_sw, lwp_in_cloud, cloud % re_liq(jcol_var_500, jlev_var_501), od_sw_liq, scat_od_sw_liq, g_sw_liq)
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
            CALL calc_ice_optics_fu_lw(16, config % cloud_optics % ice_coeff_lw, iwp_in_cloud, cloud % re_ice(jcol_var_500, jlev_var_501), od_lw_ice, scat_od_lw_ice, g_lw_ice)
            CALL calc_ice_optics_fu_sw(14, config % cloud_optics % ice_coeff_sw, iwp_in_cloud, cloud % re_ice(jcol_var_500, jlev_var_501), od_sw_ice, scat_od_sw_ice, g_sw_ice)
            CALL calc_ice_optics_yi_lw(16, config % cloud_optics % ice_coeff_lw, iwp_in_cloud, cloud % re_ice(jcol_var_500, jlev_var_501), od_lw_ice, scat_od_lw_ice, g_lw_ice)
            CALL calc_ice_optics_yi_sw(14, config % cloud_optics % ice_coeff_sw, iwp_in_cloud, cloud % re_ice(jcol_var_500, jlev_var_501), od_sw_ice, scat_od_sw_ice, g_sw_ice)
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
          DO jb_var_502 = 1, 16
            od_lw_cloud_var_493(jb_var_502, jlev_var_501, jcol_var_500) = od_lw_liq(jb_var_502) - scat_od_lw_liq(jb_var_502) + od_lw_ice(jb_var_502) - scat_od_lw_ice(jb_var_502)
          END DO
          DO jb_var_502 = 1, 14
            od_sw_cloud_var_496(jb_var_502, jlev_var_501, jcol_var_500) = od_sw_liq(jb_var_502) + od_sw_ice(jb_var_502)
            g_sw_cloud_var_498(jb_var_502, jlev_var_501, jcol_var_500) = (g_sw_liq(jb_var_502) * scat_od_sw_liq(jb_var_502) + g_sw_ice(jb_var_502) * scat_od_sw_ice(jb_var_502)) / (scat_od_sw_liq(jb_var_502) + scat_od_sw_ice(jb_var_502))
            ssa_sw_cloud_var_497(jb_var_502, jlev_var_501, jcol_var_500) = (scat_od_sw_liq(jb_var_502) + scat_od_sw_ice(jb_var_502)) / (od_sw_liq(jb_var_502) + od_sw_ice(jb_var_502))
          END DO
        END IF
      END DO
    END DO
  END SUBROUTINE cloud_optics_fn_500
END MODULE radiation_cloud_optics
MODULE radiation_ifs_rrtm
  USE radiation_config, ONLY: config_type
  USE radiation_single_level, ONLY: single_level_type
  USE radiation_thermodynamics, ONLY: thermodynamics_type
  USE radiation_gas, ONLY: assert_units_gas, gas_type
  USE yoerrtwn, ONLY: delwave, totplnk
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE gas_optics(ncol_var_503, nlev_var_504, istartcol_var_505, iendcol_var_506, config, single_level, thermodynamics, gas, od_lw_var_508, od_sw_var_509, ssa_sw_var_510, lw_albedo_var_507, planck_hl_var_511, lw_emission_var_512, incoming_sw_var_513)
    USE radiation_config, ONLY: config_type
    USE radiation_single_level, ONLY: single_level_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE radiation_gas, ONLY: assert_units_gas, gas_type
    INTEGER, INTENT(IN) :: ncol_var_503
    INTEGER, INTENT(IN) :: nlev_var_504
    INTEGER, INTENT(IN) :: istartcol_var_505, iendcol_var_506
    TYPE(config_type), INTENT(IN) :: config
    TYPE(single_level_type), INTENT(IN) :: single_level
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics
    TYPE(gas_type), INTENT(IN) :: gas
    REAL(KIND = 8), DIMENSION(140, istartcol_var_505 : iendcol_var_506), INTENT(IN), OPTIONAL :: lw_albedo_var_507
    REAL(KIND = 8), DIMENSION(140, nlev_var_504, istartcol_var_505 : iendcol_var_506), INTENT(OUT) :: od_lw_var_508
    REAL(KIND = 8), DIMENSION(112, nlev_var_504, istartcol_var_505 : iendcol_var_506), INTENT(OUT) :: od_sw_var_509, ssa_sw_var_510
    REAL(KIND = 8), DIMENSION(140, nlev_var_504 + 1, istartcol_var_505 : iendcol_var_506), INTENT(OUT), OPTIONAL :: planck_hl_var_511
    REAL(KIND = 8), DIMENSION(140, istartcol_var_505 : iendcol_var_506), INTENT(OUT), OPTIONAL :: lw_emission_var_512
    REAL(KIND = 8), DIMENSION(112, istartcol_var_505 : iendcol_var_506), INTENT(OUT), OPTIONAL :: incoming_sw_var_513
    REAL(KIND = 8) :: incoming_sw_scale(istartcol_var_505 : iendcol_var_506)
    REAL(KIND = 8) :: zod_lw(140, nlev_var_504, istartcol_var_505 : iendcol_var_506)
    REAL(KIND = 8) :: zod_sw(istartcol_var_505 : iendcol_var_506, nlev_var_504, 112)
    REAL(KIND = 8) :: zssa_sw(istartcol_var_505 : iendcol_var_506, nlev_var_504, 112)
    REAL(KIND = 8) :: zincsol(istartcol_var_505 : iendcol_var_506, 112)
    REAL(KIND = 8) :: zcolmol(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zcoldry(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zwbrodl(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zcolbrd(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zwkl(istartcol_var_505 : iendcol_var_506, 35, nlev_var_504)
    REAL(KIND = 8) :: zwx(istartcol_var_505 : iendcol_var_506, 4, nlev_var_504)
    REAL(KIND = 8) :: zfluxfac_var_514, zpi
    REAL(KIND = 8) :: ztauaerl(istartcol_var_505 : iendcol_var_506, nlev_var_504, 16)
    REAL(KIND = 8) :: zfac00(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zfac01(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zfac10(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zfac11(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zforfac(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zforfrac(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    INTEGER :: indfor_var_515(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    INTEGER :: indminor_var_516(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zscaleminor(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zscaleminorn2(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zminorfrac(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zrat_h2oco2(istartcol_var_505 : iendcol_var_506, nlev_var_504), zrat_h2oco2_1(istartcol_var_505 : iendcol_var_506, nlev_var_504), zrat_h2oo3(istartcol_var_505 : iendcol_var_506, nlev_var_504), zrat_h2oo3_1(istartcol_var_505 : iendcol_var_506, nlev_var_504), zrat_h2on2o(istartcol_var_505 : iendcol_var_506, nlev_var_504), zrat_h2on2o_1(istartcol_var_505 : iendcol_var_506, nlev_var_504), zrat_h2och4(istartcol_var_505 : iendcol_var_506, nlev_var_504), zrat_h2och4_1(istartcol_var_505 : iendcol_var_506, nlev_var_504), zrat_n2oco2(istartcol_var_505 : iendcol_var_506, nlev_var_504), zrat_n2oco2_1(istartcol_var_505 : iendcol_var_506, nlev_var_504), zrat_o3co2(istartcol_var_505 : iendcol_var_506, nlev_var_504), zrat_o3co2_1(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    INTEGER :: jp_var_517(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    INTEGER :: jt_var_518(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    INTEGER :: jt1_var_519(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zoneminus, zoneminus_array(istartcol_var_505 : iendcol_var_506)
    REAL(KIND = 8) :: zcolh2o(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zcolco2(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zcolo3(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zcoln2o(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zcolch4(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zcolo2(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zco2mult(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    INTEGER :: ilaytrop(istartcol_var_505 : iendcol_var_506)
    INTEGER :: ilayswtch(istartcol_var_505 : iendcol_var_506)
    INTEGER :: ilaylow(istartcol_var_505 : iendcol_var_506)
    REAL(KIND = 8) :: zpavel(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: ztavel(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zpz(istartcol_var_505 : iendcol_var_506, 0 : nlev_var_504)
    REAL(KIND = 8) :: ztz(istartcol_var_505 : iendcol_var_506, 0 : nlev_var_504)
    REAL(KIND = 8) :: zselffac(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zselffrac(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    INTEGER :: indself_var_520(istartcol_var_505 : iendcol_var_506, nlev_var_504)
    REAL(KIND = 8) :: zpfrac(istartcol_var_505 : iendcol_var_506, 140, nlev_var_504)
    INTEGER :: ireflect(istartcol_var_505 : iendcol_var_506)
    REAL(KIND = 8) :: pressure_fl_var_521(ncol_var_503, nlev_var_504), temperature_fl_var_522(ncol_var_503, nlev_var_504)
    INTEGER :: istartlev_var_523, iendlev_var_524
    LOGICAL :: do_sw, do_lw
    INTEGER :: jlev_var_525, jg_var_526, jcol_var_527
    do_sw = .TRUE.
    do_lw = .TRUE.
    iendlev_var_524 = UBOUND(gas % mixing_ratio, 2)
    istartlev_var_523 = iendlev_var_524 - nlev_var_504 + 1
    zpi = 2.0D0 * ASIN(1.0D0)
    zfluxfac_var_514 = zpi * 10000.0
    zoneminus = 0.999999D0
    DO jcol_var_527 = istartcol_var_505, iendcol_var_506
      zoneminus_array(jcol_var_527) = 0.999999D0
    END DO
    DO jlev_var_525 = 1, nlev_var_504
      DO jcol_var_527 = istartcol_var_505, iendcol_var_506
        pressure_fl_var_521(jcol_var_527, jlev_var_525) = thermodynamics % pressure_fl(jcol_var_527, jlev_var_525)
        temperature_fl_var_522(jcol_var_527, jlev_var_525) = thermodynamics % temperature_fl(jcol_var_527, jlev_var_525)
      END DO
    END DO
    CALL assert_units_gas(gas, 0)
    CALL rrtm_prepare_gases(istartcol_var_505, iendcol_var_506, ncol_var_503, nlev_var_504, thermodynamics % pressure_hl(:, istartlev_var_523 : iendlev_var_524 + 1), pressure_fl_var_521, thermodynamics % temperature_hl(:, istartlev_var_523 : iendlev_var_524 + 1), temperature_fl_var_522, gas % mixing_ratio(:, istartlev_var_523 : iendlev_var_524, 1), gas % mixing_ratio(:, istartlev_var_523 : iendlev_var_524, 2), gas % mixing_ratio(:, istartlev_var_523 : iendlev_var_524, 6), gas % mixing_ratio(:, istartlev_var_523 : iendlev_var_524, 4), gas % mixing_ratio(:, istartlev_var_523 : iendlev_var_524, 12), gas % mixing_ratio(:, istartlev_var_523 : iendlev_var_524, 8), gas % mixing_ratio(:, istartlev_var_523 : iendlev_var_524, 9), gas % mixing_ratio(:, istartlev_var_523 : iendlev_var_524, 10), gas % mixing_ratio(:, istartlev_var_523 : iendlev_var_524, 11), gas % mixing_ratio(:, istartlev_var_523 : iendlev_var_524, 3), zcoldry, zwbrodl, zwkl, zwx, zpavel, ztavel, zpz, ztz, ireflect)
    CALL rrtm_setcoef_140gp(istartcol_var_505, iendcol_var_506, nlev_var_504, zcoldry, zwbrodl, zwkl, zfac00, zfac01, zfac10, zfac11, zforfac, zforfrac, indfor_var_515, jp_var_517, jt_var_518, jt1_var_519, zcolh2o, zcolco2, zcolo3, zcoln2o, zcolch4, zcolo2, zco2mult, zcolbrd, ilaytrop, ilayswtch, ilaylow, zpavel, ztavel, zselffac, zselffrac, indself_var_520, indminor_var_516, zscaleminor, zscaleminorn2, zminorfrac, zrat_h2oco2, zrat_h2oco2_1, zrat_h2oo3, zrat_h2oo3_1, zrat_h2on2o, zrat_h2on2o_1, zrat_h2och4, zrat_h2och4_1, zrat_n2oco2, zrat_n2oco2_1, zrat_o3co2, zrat_o3co2_1)
    DO jg_var_526 = 1, 16
      DO jlev_var_525 = 1, nlev_var_504
        DO jcol_var_527 = istartcol_var_505, iendcol_var_506
          ztauaerl(jcol_var_527, jlev_var_525, jg_var_526) = 0.0D0
        END DO
      END DO
    END DO
    CALL rrtm_gas_optical_depth(istartcol_var_505, iendcol_var_506, nlev_var_504, zod_lw, zpavel, zcoldry, zcolbrd, zwx, ztauaerl, zfac00, zfac01, zfac10, zfac11, zforfac, zforfrac, indfor_var_515, jp_var_517, jt_var_518, jt1_var_519, zoneminus, zcolh2o, zcolco2, zcolo3, zcoln2o, zcolch4, zcolo2, zco2mult, ilaytrop, ilayswtch, ilaylow, zselffac, zselffrac, indself_var_520, zpfrac, indminor_var_516, zscaleminor, zscaleminorn2, zminorfrac, zrat_h2oco2, zrat_h2oco2_1, zrat_h2oo3, zrat_h2oo3_1, zrat_h2on2o, zrat_h2on2o_1, zrat_h2och4, zrat_h2och4_1, zrat_n2oco2, zrat_n2oco2_1, zrat_o3co2, zrat_o3co2_1)
    CALL planck_function_atmos(nlev_var_504, istartcol_var_505, iendcol_var_506, config, thermodynamics, zpfrac, planck_hl_var_511)
    CALL planck_function_surf(istartcol_var_505, iendcol_var_506, config, single_level % skin_temperature, zpfrac(:, :, 1), lw_emission_var_512)
    DO jcol_var_527 = istartcol_var_505, iendcol_var_506
      DO jg_var_526 = 1, 140
        lw_emission_var_512(jg_var_526, jcol_var_527) = lw_emission_var_512(jg_var_526, jcol_var_527) * (1.0D0 - lw_albedo_var_507(jg_var_526, jcol_var_527))
      END DO
    END DO
    DO jcol_var_527 = istartcol_var_505, iendcol_var_506
      DO jlev_var_525 = 1, nlev_var_504
        DO jg_var_526 = 1, 140
          od_lw_var_508(jg_var_526, jlev_var_525, jcol_var_527) = MAX(1D-15, zod_lw(jg_var_526, nlev_var_504 + 1 - jlev_var_525, jcol_var_527))
        END DO
      END DO
    END DO
    CALL srtm_setcoef(istartcol_var_505, iendcol_var_506, nlev_var_504, zpavel, ztavel, zcoldry, zwkl, ilaytrop, zcolch4, zcolco2, zcolh2o, zcolmol, zcolo2, zcolo3, zforfac, zforfrac, indfor_var_515, zselffac, zselffrac, indself_var_520, zfac00, zfac01, zfac10, zfac11, jp_var_517, jt_var_518, jt1_var_519, single_level % cos_sza(istartcol_var_505 : iendcol_var_506))
    DO jg_var_526 = 1, 112
      DO jlev_var_525 = 1, nlev_var_504
        DO jcol_var_527 = istartcol_var_505, iendcol_var_506
          zod_sw(jcol_var_527, jlev_var_525, jg_var_526) = 0.0D0
          zssa_sw(jcol_var_527, jlev_var_525, jg_var_526) = 0.0D0
        END DO
      END DO
    END DO
    DO jg_var_526 = 1, 112
      DO jcol_var_527 = istartcol_var_505, iendcol_var_506
        zincsol(jcol_var_527, jg_var_526) = 0.0D0
      END DO
    END DO
    CALL srtm_gas_optical_depth(istartcol_var_505, iendcol_var_506, nlev_var_504, zoneminus_array, single_level % cos_sza(istartcol_var_505 : iendcol_var_506), ilaytrop, zcolch4, zcolco2, zcolh2o, zcolmol, zcolo2, zcolo3, zforfac, zforfrac, indfor_var_515, zselffac, zselffrac, indself_var_520, zfac00, zfac01, zfac10, zfac11, jp_var_517, jt_var_518, jt1_var_519, zod_sw, zssa_sw, zincsol)
    DO jg_var_526 = 1, 112
      DO jcol_var_527 = istartcol_var_505, iendcol_var_506
        zincsol(jcol_var_527, jg_var_526) = zincsol(jcol_var_527, jg_var_526) * single_level % spectral_solar_scaling(config % i_band_from_reordered_g_sw(jg_var_526))
      END DO
    END DO
    DO jcol_var_527 = istartcol_var_505, iendcol_var_506
      IF (single_level % cos_sza(jcol_var_527) > 0.0D0) THEN
        incoming_sw_scale(jcol_var_527) = 1.0D0 / SUM(zincsol(jcol_var_527, :))
      ELSE
        incoming_sw_scale(jcol_var_527) = 1.0D0
      END IF
    END DO
    DO jcol_var_527 = istartcol_var_505, iendcol_var_506
      DO jlev_var_525 = 1, nlev_var_504
        DO jg_var_526 = 1, 112
          od_sw_var_509(jg_var_526, nlev_var_504 + 1 - jlev_var_525, jcol_var_527) = MAX(0.0D0, zod_sw(jcol_var_527, jlev_var_525, jg_var_526))
          ssa_sw_var_510(jg_var_526, nlev_var_504 + 1 - jlev_var_525, jcol_var_527) = zssa_sw(jcol_var_527, jlev_var_525, jg_var_526)
        END DO
      END DO
      DO jg_var_526 = 1, 112
        incoming_sw_var_513(jg_var_526, jcol_var_527) = incoming_sw_scale(jcol_var_527) * zincsol(jcol_var_527, jg_var_526)
      END DO
    END DO
  END SUBROUTINE gas_optics
  SUBROUTINE planck_function_atmos(nlev_var_528, istartcol_var_529, iendcol_var_530, config, thermodynamics, pfrac_var_531, planck_hl_var_532)
    USE radiation_config, ONLY: config_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE yoerrtwn, ONLY: delwave, totplnk
    INTEGER, INTENT(IN) :: nlev_var_528
    INTEGER, INTENT(IN) :: istartcol_var_529, iendcol_var_530
    TYPE(config_type), INTENT(IN) :: config
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics
    REAL(KIND = 8), INTENT(IN) :: pfrac_var_531(istartcol_var_529 : iendcol_var_530, 140, nlev_var_528)
    REAL(KIND = 8), DIMENSION(140, nlev_var_528 + 1, istartcol_var_529 : iendcol_var_530), INTENT(OUT) :: planck_hl_var_532
    REAL(KIND = 8), DIMENSION(istartcol_var_529 : iendcol_var_530, 16) :: planck_store_var_533
    REAL(KIND = 8), DIMENSION(istartcol_var_529 : iendcol_var_530) :: frac_var_534
    INTEGER, DIMENSION(istartcol_var_529 : iendcol_var_530) :: ind_var_535
    REAL(KIND = 8) :: temperature_var_536
    REAL(KIND = 8) :: factor_var_537, planck_tmp(istartcol_var_529 : iendcol_var_530, 140)
    REAL(KIND = 8) :: zfluxfac_var_538
    INTEGER :: jlev_var_539, jg_var_540, iband_var_541, jband_var_542, jcol_var_543, ilevoffset
    zfluxfac_var_538 = 2.0D0 * ASIN(1.0D0) * 10000.0D0
    ilevoffset = UBOUND(thermodynamics % temperature_hl, 2) - nlev_var_528 - 1
    DO jlev_var_539 = 1, nlev_var_528 + 1
      DO jcol_var_543 = istartcol_var_529, iendcol_var_530
        temperature_var_536 = thermodynamics % temperature_hl(jcol_var_543, jlev_var_539 + ilevoffset)
        IF (temperature_var_536 < 339.0D0 .AND. temperature_var_536 >= 160.0D0) THEN
          ind_var_535(jcol_var_543) = INT(temperature_var_536 - 159.0D0)
          frac_var_534(jcol_var_543) = temperature_var_536 - INT(temperature_var_536)
        ELSE IF (temperature_var_536 >= 339.0D0) THEN
          ind_var_535(jcol_var_543) = 180
          frac_var_534(jcol_var_543) = temperature_var_536 - 339.0D0
        ELSE
          ind_var_535(jcol_var_543) = 1
          frac_var_534(jcol_var_543) = 0.0D0
        END IF
      END DO
      DO jband_var_542 = 1, 16
        factor_var_537 = zfluxfac_var_538 * delwave(jband_var_542)
        DO jcol_var_543 = istartcol_var_529, iendcol_var_530
          planck_store_var_533(jcol_var_543, jband_var_542) = factor_var_537 * (totplnk(ind_var_535(jcol_var_543), jband_var_542) + frac_var_534(jcol_var_543) * (totplnk(ind_var_535(jcol_var_543) + 1, jband_var_542) - totplnk(ind_var_535(jcol_var_543), jband_var_542)))
        END DO
      END DO
      IF (jlev_var_539 == 1) THEN
        DO jg_var_540 = 1, 140
          iband_var_541 = config % i_band_from_g_lw(jg_var_540)
          DO jcol_var_543 = istartcol_var_529, iendcol_var_530
            planck_hl_var_532(jg_var_540, 1, jcol_var_543) = planck_store_var_533(jcol_var_543, iband_var_541) * pfrac_var_531(jcol_var_543, jg_var_540, nlev_var_528)
          END DO
        END DO
      ELSE
        DO jg_var_540 = 1, 140
          iband_var_541 = config % i_band_from_g_lw(jg_var_540)
          DO jcol_var_543 = istartcol_var_529, iendcol_var_530
            planck_tmp(jcol_var_543, jg_var_540) = planck_store_var_533(jcol_var_543, iband_var_541) * pfrac_var_531(jcol_var_543, jg_var_540, nlev_var_528 + 2 - jlev_var_539)
          END DO
        END DO
        DO jcol_var_543 = istartcol_var_529, iendcol_var_530
          DO jg_var_540 = 1, 140
            planck_hl_var_532(jg_var_540, jlev_var_539, jcol_var_543) = planck_tmp(jcol_var_543, jg_var_540)
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE planck_function_atmos
  SUBROUTINE planck_function_surf(istartcol_var_544, iendcol_var_545, config, temperature_var_546, pfrac_var_547, planck_surf)
    USE radiation_config, ONLY: config_type
    USE yoerrtwn, ONLY: delwave, totplnk
    INTEGER, INTENT(IN) :: istartcol_var_544, iendcol_var_545
    TYPE(config_type), INTENT(IN) :: config
    REAL(KIND = 8), INTENT(IN) :: temperature_var_546(:)
    REAL(KIND = 8), INTENT(IN) :: pfrac_var_547(istartcol_var_544 : iendcol_var_545, 140)
    REAL(KIND = 8), DIMENSION(140, istartcol_var_544 : iendcol_var_545), INTENT(OUT) :: planck_surf
    REAL(KIND = 8), DIMENSION(istartcol_var_544 : iendcol_var_545, 16) :: planck_store_var_548
    REAL(KIND = 8), DIMENSION(istartcol_var_544 : iendcol_var_545) :: frac_var_549
    INTEGER, DIMENSION(istartcol_var_544 : iendcol_var_545) :: ind_var_550
    REAL(KIND = 8) :: tsurf
    REAL(KIND = 8) :: factor_var_551
    REAL(KIND = 8) :: zfluxfac_var_552
    INTEGER :: jg_var_553, iband_var_554, jband_var_555, jcol_var_556
    zfluxfac_var_552 = 2.0D0 * ASIN(1.0D0) * 10000.0D0
    DO jcol_var_556 = istartcol_var_544, iendcol_var_545
      tsurf = temperature_var_546(jcol_var_556)
      IF (tsurf < 339.0D0 .AND. tsurf >= 160.0D0) THEN
        ind_var_550(jcol_var_556) = INT(tsurf - 159.0D0)
        frac_var_549(jcol_var_556) = tsurf - INT(tsurf)
      ELSE IF (tsurf >= 339.0D0) THEN
        ind_var_550(jcol_var_556) = 180
        frac_var_549(jcol_var_556) = tsurf - 339.0D0
      ELSE
        ind_var_550(jcol_var_556) = 1
        frac_var_549(jcol_var_556) = 0.0D0
      END IF
    END DO
    DO jband_var_555 = 1, 16
      factor_var_551 = zfluxfac_var_552 * delwave(jband_var_555)
      DO jcol_var_556 = istartcol_var_544, iendcol_var_545
        planck_store_var_548(jcol_var_556, jband_var_555) = factor_var_551 * (totplnk(ind_var_550(jcol_var_556), jband_var_555) + frac_var_549(jcol_var_556) * (totplnk(ind_var_550(jcol_var_556) + 1, jband_var_555) - totplnk(ind_var_550(jcol_var_556), jband_var_555)))
      END DO
    END DO
    DO jg_var_553 = 1, 140
      iband_var_554 = config % i_band_from_g_lw(jg_var_553)
      DO jcol_var_556 = istartcol_var_544, iendcol_var_545
        planck_surf(jg_var_553, jcol_var_556) = planck_store_var_548(jcol_var_556, iband_var_554) * pfrac_var_547(jcol_var_556, jg_var_553)
      END DO
    END DO
  END SUBROUTINE planck_function_surf
END MODULE radiation_ifs_rrtm
MODULE radiation_mcica_lw
  CONTAINS
  SUBROUTINE solver_mcica_lw(nlev_var_557, istartcol_var_558, iendcol_var_559, config, single_level, cloud, od_var_560, ssa_var_561, g_var_562, od_cloud_var_563, ssa_cloud_var_564, g_cloud_var_565, planck_hl_var_566, emission, albedo_var_567, flux)
    USE radiation_config, ONLY: config_type
    USE radiation_single_level, ONLY: single_level_type
    USE radiation_cloud, ONLY: cloud_type
    USE radiation_flux, ONLY: flux_type
    USE radiation_two_stream, ONLY: calc_no_scattering_transmittance_lw
    USE radiation_adding_ica_lw, ONLY: calc_fluxes_no_scattering_lw
    USE radiation_cloud_generator, ONLY: cloud_generator
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlev_var_557
    INTEGER, INTENT(IN) :: istartcol_var_558, iendcol_var_559
    TYPE(config_type), INTENT(IN) :: config
    TYPE(single_level_type), INTENT(IN) :: single_level
    TYPE(cloud_type), INTENT(IN) :: cloud
    REAL(KIND = 8), INTENT(IN), DIMENSION(140, nlev_var_557, istartcol_var_558 : iendcol_var_559) :: od_var_560
    REAL(KIND = 8), INTENT(IN), DIMENSION(0, nlev_var_557, istartcol_var_558 : iendcol_var_559) :: ssa_var_561, g_var_562
    REAL(KIND = 8), INTENT(IN), DIMENSION(16, nlev_var_557, istartcol_var_558 : iendcol_var_559) :: od_cloud_var_563
    REAL(KIND = 8), INTENT(IN), DIMENSION(0, nlev_var_557, istartcol_var_558 : iendcol_var_559) :: ssa_cloud_var_564, g_cloud_var_565
    REAL(KIND = 8), INTENT(IN), DIMENSION(140, nlev_var_557 + 1, istartcol_var_558 : iendcol_var_559) :: planck_hl_var_566
    REAL(KIND = 8), INTENT(IN), DIMENSION(140, istartcol_var_558 : iendcol_var_559) :: emission, albedo_var_567
    TYPE(flux_type), INTENT(INOUT) :: flux
    REAL(KIND = 8), DIMENSION(140, nlev_var_557) :: ref_clear_var_568, trans_clear_var_569, reflectance_var_570, transmittance_var_571
    REAL(KIND = 8), DIMENSION(140, nlev_var_557) :: source_up_clear, source_dn_clear, source_up_var_572, source_dn_var_573
    REAL(KIND = 8), DIMENSION(140, nlev_var_557 + 1) :: flux_up_var_574, flux_dn_var_575
    REAL(KIND = 8), DIMENSION(140, nlev_var_557 + 1) :: flux_up_clear, flux_dn_clear
    REAL(KIND = 8), DIMENSION(140) :: od_total_var_576, ssa_total_var_577, g_total_var_578
    REAL(KIND = 8), DIMENSION(140, nlev_var_557) :: od_scaling_var_579
    REAL(KIND = 8), DIMENSION(140) :: od_cloud_new_var_580
    REAL(KIND = 8) :: total_cloud_cover_var_581
    LOGICAL :: is_clear_sky_layer(nlev_var_557)
    REAL(KIND = 8) :: sum_up_var_582, sum_dn
    INTEGER :: i_cloud_top
    INTEGER :: ng_var_583
    INTEGER :: jlev_var_584, jcol_var_585, jg_var_586
    ng_var_583 = 140
    DO jcol_var_585 = istartcol_var_558, iendcol_var_559
      CALL calc_no_scattering_transmittance_lw(140 * nlev_var_557, od_var_560(:, :, jcol_var_585), planck_hl_var_566(:, 1 : nlev_var_557, jcol_var_585), planck_hl_var_566(:, 2 : nlev_var_557 + 1, jcol_var_585), trans_clear_var_569, source_up_clear, source_dn_clear)
      ref_clear_var_568 = 0.0D0
      CALL calc_fluxes_no_scattering_lw(140, nlev_var_557, trans_clear_var_569, source_up_clear, source_dn_clear, emission(:, jcol_var_585), albedo_var_567(:, jcol_var_585), flux_up_clear, flux_dn_clear)
      DO jlev_var_584 = 1, nlev_var_557 + 1
        sum_up_var_582 = 0.0D0
        sum_dn = 0.0D0
        DO jg_var_586 = 1, ng_var_583
          sum_up_var_582 = 0.0D0 + flux_up_clear(jg_var_586, jlev_var_584)
          sum_dn = 0.0D0 + flux_dn_clear(jg_var_586, jlev_var_584)
        END DO
        flux % lw_up_clear(jcol_var_585, jlev_var_584) = sum_up_var_582
        flux % lw_dn_clear(jcol_var_585, jlev_var_584) = sum_dn
      END DO
      flux % lw_dn_surf_clear_g(:, jcol_var_585) = flux_dn_clear(:, nlev_var_557 + 1)
      CALL cloud_generator(140, nlev_var_557, 1, single_level % iseed(jcol_var_585) + 997, 1D-06, cloud % fraction(jcol_var_585, :), cloud % overlap_param(jcol_var_585, :), 0.5D0, cloud % fractional_std(jcol_var_585, :), config % pdf_sampler, od_scaling_var_579, total_cloud_cover_var_581, use_beta_overlap = .FALSE., use_vectorizable_generator = .FALSE.)
      flux % cloud_cover_lw(jcol_var_585) = total_cloud_cover_var_581
      IF (total_cloud_cover_var_581 >= 1D-06) THEN
        is_clear_sky_layer = .TRUE.
        i_cloud_top = nlev_var_557 + 1
        DO jlev_var_584 = 1, nlev_var_557
          IF (cloud % fraction(jcol_var_585, jlev_var_584) >= 1D-06) THEN
            is_clear_sky_layer(jlev_var_584) = .FALSE.
            IF (i_cloud_top > jlev_var_584) THEN
              i_cloud_top = jlev_var_584
            END IF
            DO jg_var_586 = 1, ng_var_583
              od_cloud_new_var_580(jg_var_586) = od_scaling_var_579(jg_var_586, jlev_var_584) * od_cloud_var_563(config % i_band_from_reordered_g_lw(jg_var_586), jlev_var_584, jcol_var_585)
              od_total_var_576(jg_var_586) = od_var_560(jg_var_586, jlev_var_584, jcol_var_585) + od_cloud_new_var_580(jg_var_586)
              ssa_total_var_577(jg_var_586) = 0.0D0
              g_total_var_578(jg_var_586) = 0.0D0
            END DO
            CALL calc_no_scattering_transmittance_lw(ng_var_583, od_total_var_576, planck_hl_var_566(:, jlev_var_584, jcol_var_585), planck_hl_var_566(:, jlev_var_584 + 1, jcol_var_585), transmittance_var_571(:, jlev_var_584), source_up_var_572(:, jlev_var_584), source_dn_var_573(:, jlev_var_584))
          ELSE
            DO jg_var_586 = 1, ng_var_583
              reflectance_var_570(jg_var_586, jlev_var_584) = ref_clear_var_568(jg_var_586, jlev_var_584)
              transmittance_var_571(jg_var_586, jlev_var_584) = trans_clear_var_569(jg_var_586, jlev_var_584)
              source_up_var_572(jg_var_586, jlev_var_584) = source_up_clear(jg_var_586, jlev_var_584)
              source_dn_var_573(jg_var_586, jlev_var_584) = source_dn_clear(jg_var_586, jlev_var_584)
            END DO
          END IF
        END DO
        CALL calc_fluxes_no_scattering_lw(ng_var_583, nlev_var_557, transmittance_var_571, source_up_var_572, source_dn_var_573, emission(:, jcol_var_585), albedo_var_567(:, jcol_var_585), flux_up_var_574, flux_dn_var_575)
        DO jlev_var_584 = 1, nlev_var_557 + 1
          sum_up_var_582 = 0.0D0
          sum_dn = 0.0D0
          DO jg_var_586 = 1, ng_var_583
            sum_up_var_582 = 0.0D0 + flux_up_var_574(jg_var_586, jlev_var_584)
            sum_dn = 0.0D0 + flux_dn_var_575(jg_var_586, jlev_var_584)
          END DO
          flux % lw_up(jcol_var_585, jlev_var_584) = sum_up_var_582
          flux % lw_dn(jcol_var_585, jlev_var_584) = sum_dn
        END DO
        DO jlev_var_584 = 1, nlev_var_557 + 1
          flux % lw_up(jcol_var_585, jlev_var_584) = total_cloud_cover_var_581 * flux % lw_up(jcol_var_585, jlev_var_584) + (1.0D0 - total_cloud_cover_var_581) * flux % lw_up_clear(jcol_var_585, jlev_var_584)
          flux % lw_dn(jcol_var_585, jlev_var_584) = total_cloud_cover_var_581 * flux % lw_dn(jcol_var_585, jlev_var_584) + (1.0D0 - total_cloud_cover_var_581) * flux % lw_dn_clear(jcol_var_585, jlev_var_584)
        END DO
        flux % lw_dn_surf_g(:, jcol_var_585) = total_cloud_cover_var_581 * flux_dn_var_575(:, nlev_var_557 + 1) + (1.0D0 - total_cloud_cover_var_581) * flux % lw_dn_surf_clear_g(:, jcol_var_585)
      ELSE
        DO jlev_var_584 = 1, nlev_var_557 + 1
          flux % lw_up(jcol_var_585, jlev_var_584) = flux % lw_up_clear(jcol_var_585, jlev_var_584)
          flux % lw_dn(jcol_var_585, jlev_var_584) = flux % lw_dn_clear(jcol_var_585, jlev_var_584)
        END DO
        flux % lw_dn_surf_g(:, jcol_var_585) = flux % lw_dn_surf_clear_g(:, jcol_var_585)
      END IF
    END DO
  END SUBROUTINE solver_mcica_lw
END MODULE radiation_mcica_lw
MODULE radiation_mcica_sw
  CONTAINS
  SUBROUTINE solver_mcica_sw(nlev_var_587, istartcol_var_588, iendcol_var_589, config, single_level, cloud, od_var_590, ssa_var_591, g_var_592, od_cloud_var_593, ssa_cloud_var_594, g_cloud_var_595, albedo_direct, albedo_diffuse, incoming_sw_var_596, flux)
    USE radiation_config, ONLY: config_type
    USE radiation_single_level, ONLY: single_level_type
    USE radiation_cloud, ONLY: cloud_type
    USE radiation_flux, ONLY: flux_type
    USE radiation_two_stream, ONLY: calc_ref_trans_sw
    USE radiation_adding_ica_sw, ONLY: adding_ica_sw
    USE radiation_cloud_generator, ONLY: cloud_generator
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlev_var_587
    INTEGER, INTENT(IN) :: istartcol_var_588, iendcol_var_589
    TYPE(config_type), INTENT(IN) :: config
    TYPE(single_level_type), INTENT(IN) :: single_level
    TYPE(cloud_type), INTENT(IN) :: cloud
    REAL(KIND = 8), INTENT(IN), DIMENSION(112, nlev_var_587, istartcol_var_588 : iendcol_var_589) :: od_var_590, ssa_var_591, g_var_592
    REAL(KIND = 8), INTENT(IN), DIMENSION(14, nlev_var_587, istartcol_var_588 : iendcol_var_589) :: od_cloud_var_593, ssa_cloud_var_594, g_cloud_var_595
    REAL(KIND = 8), INTENT(IN), DIMENSION(112, istartcol_var_588 : iendcol_var_589) :: albedo_direct, albedo_diffuse, incoming_sw_var_596
    TYPE(flux_type), INTENT(INOUT) :: flux
    REAL(KIND = 8) :: cos_sza_var_597
    REAL(KIND = 8), DIMENSION(112, nlev_var_587) :: ref_clear_var_598, trans_clear_var_599, reflectance_var_600, transmittance_var_601
    REAL(KIND = 8), DIMENSION(112, nlev_var_587) :: ref_dir_clear, trans_dir_diff_clear, ref_dir_var_602, trans_dir_diff_var_603
    REAL(KIND = 8), DIMENSION(112, nlev_var_587) :: trans_dir_dir_clear, trans_dir_dir_var_604
    REAL(KIND = 8), DIMENSION(112, nlev_var_587 + 1) :: flux_up_var_605, flux_dn_diffuse_var_606, flux_dn_direct_var_607
    REAL(KIND = 8), DIMENSION(112) :: od_total_var_608, ssa_total_var_609, g_total_var_610
    REAL(KIND = 8) :: scat_od_var_611
    REAL(KIND = 8), DIMENSION(112, nlev_var_587) :: od_scaling_var_612
    REAL(KIND = 8), DIMENSION(112) :: od_cloud_new_var_613
    REAL(KIND = 8), DIMENSION(112, nlev_var_587 + 1) :: tmp_work_albedo, tmp_work_source
    REAL(KIND = 8), DIMENSION(112, nlev_var_587) :: tmp_work_inv_denominator
    REAL(KIND = 8) :: total_cloud_cover_var_614
    REAL(KIND = 8) :: sum_up_var_615, sum_dn_diff, sum_dn_dir
    INTEGER :: ng_var_616
    INTEGER :: jlev_var_617, jcol_var_618, jg_var_619
    ng_var_616 = 112
    DO jcol_var_618 = istartcol_var_588, iendcol_var_589
      IF (single_level % cos_sza(jcol_var_618) > 0.0D0) THEN
        cos_sza_var_597 = single_level % cos_sza(jcol_var_618)
        CALL calc_ref_trans_sw(ng_var_616 * nlev_var_587, cos_sza_var_597, od_var_590(:, :, jcol_var_618), ssa_var_591(:, :, jcol_var_618), g_var_592(:, :, jcol_var_618), ref_clear_var_598, trans_clear_var_599, ref_dir_clear, trans_dir_diff_clear, trans_dir_dir_clear)
        CALL adding_ica_sw(ng_var_616, nlev_var_587, incoming_sw_var_596(:, jcol_var_618), albedo_diffuse(:, jcol_var_618), albedo_direct(:, jcol_var_618), cos_sza_var_597, ref_clear_var_598, trans_clear_var_599, ref_dir_clear, trans_dir_diff_clear, trans_dir_dir_clear, flux_up_var_605, flux_dn_diffuse_var_606, flux_dn_direct_var_607, albedo_var_20 = tmp_work_albedo, source_var_21 = tmp_work_source, inv_denominator = tmp_work_inv_denominator)
        DO jlev_var_617 = 1, nlev_var_587 + 1
          sum_up_var_615 = 0.0D0
          sum_dn_diff = 0.0D0
          sum_dn_dir = 0.0D0
          DO jg_var_619 = 1, ng_var_616
            sum_up_var_615 = 0.0D0 + flux_up_var_605(jg_var_619, jlev_var_617)
            sum_dn_diff = 0.0D0 + flux_dn_diffuse_var_606(jg_var_619, jlev_var_617)
            sum_dn_dir = 0.0D0 + flux_dn_direct_var_607(jg_var_619, jlev_var_617)
          END DO
          flux % sw_up_clear(jcol_var_618, jlev_var_617) = sum_up_var_615
          flux % sw_dn_clear(jcol_var_618, jlev_var_617) = sum_dn_diff + sum_dn_dir
          flux % sw_dn_direct_clear(jcol_var_618, jlev_var_617) = sum_dn_dir
        END DO
        DO jg_var_619 = 1, ng_var_616
          flux % sw_dn_diffuse_surf_clear_g(jg_var_619, jcol_var_618) = flux_dn_diffuse_var_606(jg_var_619, nlev_var_587 + 1)
          flux % sw_dn_direct_surf_clear_g(jg_var_619, jcol_var_618) = flux_dn_direct_var_607(jg_var_619, nlev_var_587 + 1)
        END DO
        CALL cloud_generator(ng_var_616, nlev_var_587, 1, single_level % iseed(jcol_var_618), 1D-06, cloud % fraction(jcol_var_618, :), cloud % overlap_param(jcol_var_618, :), 0.5D0, cloud % fractional_std(jcol_var_618, :), config % pdf_sampler, od_scaling_var_612, total_cloud_cover_var_614, use_beta_overlap = .FALSE., use_vectorizable_generator = .FALSE.)
        flux % cloud_cover_sw(jcol_var_618) = total_cloud_cover_var_614
        IF (total_cloud_cover_var_614 >= 1D-06) THEN
          DO jlev_var_617 = 1, nlev_var_587
            IF (cloud % fraction(jcol_var_618, jlev_var_617) >= 1D-06) THEN
              DO jg_var_619 = 1, ng_var_616
                od_cloud_new_var_613(jg_var_619) = od_scaling_var_612(jg_var_619, jlev_var_617) * od_cloud_var_593(config % i_band_from_reordered_g_sw(jg_var_619), jlev_var_617, jcol_var_618)
                od_total_var_608(jg_var_619) = od_var_590(jg_var_619, jlev_var_617, jcol_var_618) + od_cloud_new_var_613(jg_var_619)
                ssa_total_var_609(jg_var_619) = 0.0D0
                g_total_var_610(jg_var_619) = 0.0D0
                IF (od_total_var_608(jg_var_619) > 0.0D0) THEN
                  scat_od_var_611 = ssa_var_591(jg_var_619, jlev_var_617, jcol_var_618) * od_var_590(jg_var_619, jlev_var_617, jcol_var_618) + ssa_cloud_var_594(config % i_band_from_reordered_g_sw(jg_var_619), jlev_var_617, jcol_var_618) * od_cloud_new_var_613(jg_var_619)
                  ssa_total_var_609(jg_var_619) = scat_od_var_611 / od_total_var_608(jg_var_619)
                  IF (scat_od_var_611 > 0.0D0) THEN
                    g_total_var_610(jg_var_619) = (g_var_592(jg_var_619, jlev_var_617, jcol_var_618) * ssa_var_591(jg_var_619, jlev_var_617, jcol_var_618) * od_var_590(jg_var_619, jlev_var_617, jcol_var_618) + g_cloud_var_595(config % i_band_from_reordered_g_sw(jg_var_619), jlev_var_617, jcol_var_618) * ssa_cloud_var_594(config % i_band_from_reordered_g_sw(jg_var_619), jlev_var_617, jcol_var_618) * od_cloud_new_var_613(jg_var_619)) / scat_od_var_611
                  END IF
                END IF
              END DO
              CALL calc_ref_trans_sw(ng_var_616, cos_sza_var_597, od_total_var_608, ssa_total_var_609, g_total_var_610, reflectance_var_600(:, jlev_var_617), transmittance_var_601(:, jlev_var_617), ref_dir_var_602(:, jlev_var_617), trans_dir_diff_var_603(:, jlev_var_617), trans_dir_dir_var_604(:, jlev_var_617))
            ELSE
              DO jg_var_619 = 1, ng_var_616
                reflectance_var_600(jg_var_619, jlev_var_617) = ref_clear_var_598(jg_var_619, jlev_var_617)
                transmittance_var_601(jg_var_619, jlev_var_617) = trans_clear_var_599(jg_var_619, jlev_var_617)
                ref_dir_var_602(jg_var_619, jlev_var_617) = ref_dir_clear(jg_var_619, jlev_var_617)
                trans_dir_diff_var_603(jg_var_619, jlev_var_617) = trans_dir_diff_clear(jg_var_619, jlev_var_617)
                trans_dir_dir_var_604(jg_var_619, jlev_var_617) = trans_dir_dir_clear(jg_var_619, jlev_var_617)
              END DO
            END IF
          END DO
          CALL adding_ica_sw(ng_var_616, nlev_var_587, incoming_sw_var_596(:, jcol_var_618), albedo_diffuse(:, jcol_var_618), albedo_direct(:, jcol_var_618), cos_sza_var_597, reflectance_var_600, transmittance_var_601, ref_dir_var_602, trans_dir_diff_var_603, trans_dir_dir_var_604, flux_up_var_605, flux_dn_diffuse_var_606, flux_dn_direct_var_607, albedo_var_20 = tmp_work_albedo, source_var_21 = tmp_work_source, inv_denominator = tmp_work_inv_denominator)
          DO jlev_var_617 = 1, nlev_var_587 + 1
            sum_up_var_615 = 0.0D0
            sum_dn_diff = 0.0D0
            sum_dn_dir = 0.0D0
            DO jg_var_619 = 1, ng_var_616
              sum_up_var_615 = 0.0D0 + flux_up_var_605(jg_var_619, jlev_var_617)
              sum_dn_diff = 0.0D0 + flux_dn_diffuse_var_606(jg_var_619, jlev_var_617)
              sum_dn_dir = 0.0D0 + flux_dn_direct_var_607(jg_var_619, jlev_var_617)
            END DO
            flux % sw_up(jcol_var_618, jlev_var_617) = sum_up_var_615
            flux % sw_dn(jcol_var_618, jlev_var_617) = sum_dn_diff + sum_dn_dir
            flux % sw_dn_direct(jcol_var_618, jlev_var_617) = sum_dn_dir
          END DO
          DO jlev_var_617 = 1, nlev_var_587 + 1
            flux % sw_up(jcol_var_618, jlev_var_617) = total_cloud_cover_var_614 * flux % sw_up(jcol_var_618, jlev_var_617) + (1.0D0 - total_cloud_cover_var_614) * flux % sw_up_clear(jcol_var_618, jlev_var_617)
            flux % sw_dn(jcol_var_618, jlev_var_617) = total_cloud_cover_var_614 * flux % sw_dn(jcol_var_618, jlev_var_617) + (1.0D0 - total_cloud_cover_var_614) * flux % sw_dn_clear(jcol_var_618, jlev_var_617)
            flux % sw_dn_direct(jcol_var_618, jlev_var_617) = total_cloud_cover_var_614 * flux % sw_dn_direct(jcol_var_618, jlev_var_617) + (1.0D0 - total_cloud_cover_var_614) * flux % sw_dn_direct_clear(jcol_var_618, jlev_var_617)
          END DO
          DO jg_var_619 = 1, ng_var_616
            flux % sw_dn_diffuse_surf_g(jg_var_619, jcol_var_618) = flux_dn_diffuse_var_606(jg_var_619, nlev_var_587 + 1)
            flux % sw_dn_direct_surf_g(jg_var_619, jcol_var_618) = flux_dn_direct_var_607(jg_var_619, nlev_var_587 + 1)
            flux % sw_dn_diffuse_surf_g(jg_var_619, jcol_var_618) = total_cloud_cover_var_614 * flux % sw_dn_diffuse_surf_g(jg_var_619, jcol_var_618) + (1.0D0 - total_cloud_cover_var_614) * flux % sw_dn_diffuse_surf_clear_g(jg_var_619, jcol_var_618)
            flux % sw_dn_direct_surf_g(jg_var_619, jcol_var_618) = total_cloud_cover_var_614 * flux % sw_dn_direct_surf_g(jg_var_619, jcol_var_618) + (1.0D0 - total_cloud_cover_var_614) * flux % sw_dn_direct_surf_clear_g(jg_var_619, jcol_var_618)
          END DO
        ELSE
          DO jlev_var_617 = 1, nlev_var_587 + 1
            flux % sw_up(jcol_var_618, jlev_var_617) = flux % sw_up_clear(jcol_var_618, jlev_var_617)
            flux % sw_dn(jcol_var_618, jlev_var_617) = flux % sw_dn_clear(jcol_var_618, jlev_var_617)
            flux % sw_dn_direct(jcol_var_618, jlev_var_617) = flux % sw_dn_direct_clear(jcol_var_618, jlev_var_617)
          END DO
          DO jg_var_619 = 1, ng_var_616
            flux % sw_dn_diffuse_surf_g(jg_var_619, jcol_var_618) = flux % sw_dn_diffuse_surf_clear_g(jg_var_619, jcol_var_618)
            flux % sw_dn_direct_surf_g(jg_var_619, jcol_var_618) = flux % sw_dn_direct_surf_clear_g(jg_var_619, jcol_var_618)
          END DO
        END IF
      ELSE
        DO jlev_var_617 = 1, nlev_var_587 + 1
          flux % sw_up(jcol_var_618, jlev_var_617) = 0.0D0
          flux % sw_dn(jcol_var_618, jlev_var_617) = 0.0D0
          flux % sw_dn_direct(jcol_var_618, jlev_var_617) = 0.0D0
          flux % sw_up_clear(jcol_var_618, jlev_var_617) = 0.0D0
          flux % sw_dn_clear(jcol_var_618, jlev_var_617) = 0.0D0
          flux % sw_dn_direct_clear(jcol_var_618, jlev_var_617) = 0.0D0
        END DO
        DO jg_var_619 = 1, ng_var_616
          flux % sw_dn_diffuse_surf_g(jg_var_619, jcol_var_618) = 0.0D0
          flux % sw_dn_direct_surf_g(jg_var_619, jcol_var_618) = 0.0D0
          flux % sw_dn_diffuse_surf_clear_g(jg_var_619, jcol_var_618) = 0.0D0
          flux % sw_dn_direct_surf_clear_g(jg_var_619, jcol_var_618) = 0.0D0
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
  USE radiation_cloud_optics, ONLY: cloud_optics_fn_500
  USE radiation_aerosol_optics, ONLY: add_aerosol_optics
  USE radiation_mcica_lw, ONLY: solver_mcica_lw
  USE radiation_mcica_sw, ONLY: solver_mcica_sw
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE radiation(ncol_var_620, nlev_var_621, istartcol_var_622, iendcol_var_623, config, single_level, thermodynamics, gas, cloud, aerosol, flux)
    USE radiation_config, ONLY: config_type
    USE radiation_single_level, ONLY: get_albedos, single_level_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE radiation_gas, ONLY: gas_type
    USE radiation_cloud, ONLY: cloud_type, crop_cloud_fraction
    USE radiation_aerosol, ONLY: aerosol_type
    USE radiation_flux, ONLY: calc_surface_spectral, calc_toa_spectral, flux_type
    USE radiation_ifs_rrtm, ONLY: gas_optics
    USE radiation_cloud_optics, ONLY: cloud_optics_fn_500
    USE radiation_aerosol_optics, ONLY: add_aerosol_optics
    USE radiation_mcica_lw, ONLY: solver_mcica_lw
    USE radiation_mcica_sw, ONLY: solver_mcica_sw
    INTEGER, INTENT(IN) :: ncol_var_620
    INTEGER, INTENT(IN) :: nlev_var_621
    INTEGER, INTENT(IN) :: istartcol_var_622, iendcol_var_623
    TYPE(config_type), INTENT(IN) :: config
    TYPE(single_level_type), INTENT(IN) :: single_level
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics
    TYPE(gas_type), INTENT(IN) :: gas
    TYPE(cloud_type), INTENT(INOUT) :: cloud
    TYPE(aerosol_type), INTENT(IN) :: aerosol
    TYPE(flux_type), INTENT(INOUT) :: flux
    REAL(KIND = 8), DIMENSION(140, nlev_var_621, istartcol_var_622 : iendcol_var_623) :: od_lw_var_624
    REAL(KIND = 8), DIMENSION(0, nlev_var_621, istartcol_var_622 : iendcol_var_623) :: ssa_lw_var_625, g_lw_var_626
    REAL(KIND = 8), DIMENSION(16, nlev_var_621, istartcol_var_622 : iendcol_var_623) :: od_lw_cloud_var_627
    REAL(KIND = 8), DIMENSION(0, nlev_var_621, istartcol_var_622 : iendcol_var_623) :: ssa_lw_cloud_var_628, g_lw_cloud_var_629
    REAL(KIND = 8), DIMENSION(112, nlev_var_621, istartcol_var_622 : iendcol_var_623) :: od_sw_var_630, ssa_sw_var_631, g_sw_var_632
    REAL(KIND = 8), DIMENSION(14, nlev_var_621, istartcol_var_622 : iendcol_var_623) :: od_sw_cloud_var_633, ssa_sw_cloud_var_634, g_sw_cloud_var_635
    REAL(KIND = 8), DIMENSION(140, nlev_var_621 + 1, istartcol_var_622 : iendcol_var_623) :: planck_hl_var_636
    REAL(KIND = 8), DIMENSION(140, istartcol_var_622 : iendcol_var_623) :: lw_emission_var_637
    REAL(KIND = 8), DIMENSION(140, istartcol_var_622 : iendcol_var_623) :: lw_albedo_var_638
    REAL(KIND = 8), DIMENSION(112, istartcol_var_622 : iendcol_var_623) :: sw_albedo_direct_var_639
    REAL(KIND = 8), DIMENSION(112, istartcol_var_622 : iendcol_var_623) :: sw_albedo_diffuse_var_640
    REAL(KIND = 8), DIMENSION(112, istartcol_var_622 : iendcol_var_623) :: incoming_sw_var_641
    CALL global_init_fn
    CALL get_albedos(single_level, istartcol_var_622, iendcol_var_623, config, sw_albedo_direct_var_639, sw_albedo_diffuse_var_640, lw_albedo_var_638)
    CALL gas_optics(ncol_var_620, nlev_var_621, istartcol_var_622, iendcol_var_623, config, single_level, thermodynamics, gas, od_lw_var_624, od_sw_var_630, ssa_sw_var_631, lw_albedo_var_507 = lw_albedo_var_638, planck_hl_var_511 = planck_hl_var_636, lw_emission_var_512 = lw_emission_var_637, incoming_sw_var_513 = incoming_sw_var_641)
    CALL crop_cloud_fraction(cloud, istartcol_var_622, iendcol_var_623, 1D-06, 1D-09)
    CALL cloud_optics_fn_500(nlev_var_621, istartcol_var_622, iendcol_var_623, config, thermodynamics, cloud, od_lw_cloud_var_627, ssa_lw_cloud_var_628, g_lw_cloud_var_629, od_sw_cloud_var_633, ssa_sw_cloud_var_634, g_sw_cloud_var_635)
    CALL add_aerosol_optics(nlev_var_621, istartcol_var_622, iendcol_var_623, config, thermodynamics, gas, aerosol, od_lw_var_624, ssa_lw_var_625, g_lw_var_626, od_sw_var_630, ssa_sw_var_631, g_sw_var_632)
    CALL solver_mcica_lw(nlev_var_621, istartcol_var_622, iendcol_var_623, config, single_level, cloud, od_lw_var_624, ssa_lw_var_625, g_lw_var_626, od_lw_cloud_var_627, ssa_lw_cloud_var_628, g_lw_cloud_var_629, planck_hl_var_636, lw_emission_var_637, lw_albedo_var_638, flux)
    CALL solver_mcica_sw(nlev_var_621, istartcol_var_622, iendcol_var_623, config, single_level, cloud, od_sw_var_630, ssa_sw_var_631, g_sw_var_632, od_sw_cloud_var_633, ssa_sw_cloud_var_634, g_sw_cloud_var_635, sw_albedo_direct_var_639, sw_albedo_diffuse_var_640, incoming_sw_var_641, flux)
    CALL calc_surface_spectral(flux, config, istartcol_var_622, iendcol_var_623)
    CALL calc_toa_spectral(flux, config, istartcol_var_622, iendcol_var_623)
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
SUBROUTINE srtm_setcoef(kidia_var_642, kfdia_var_643, klev_var_644, pavel_var_645, ptavel_var_646, pcoldry_var_647, pwkl_var_648, klaytrop_var_649, pcolch4_var_650, pcolco2_var_651, pcolh2o_var_652, pcolmol_var_653, pcolo2_var_654, pcolo3_var_655, pforfac_var_656, pforfrac_var_657, kindfor_var_658, pselffac_var_659, pselffrac_var_660, kindself_var_661, pfac00_var_662, pfac01_var_663, pfac10_var_664, pfac11_var_665, kjp_var_666, kjt_var_667, kjt1_var_668, prmu0_var_669)
  USE yoesrtwn, ONLY: preflog_var_335, tref_var_336
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_642, kfdia_var_643
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_644
  REAL(KIND = 8), INTENT(IN) :: pavel_var_645(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(IN) :: ptavel_var_646(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(IN) :: pcoldry_var_647(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(IN) :: pwkl_var_648(kidia_var_642 : kfdia_var_643, 35, klev_var_644)
  INTEGER(KIND = 4), INTENT(INOUT) :: klaytrop_var_649(kidia_var_642 : kfdia_var_643)
  REAL(KIND = 8), INTENT(INOUT) :: pcolch4_var_650(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pcolco2_var_651(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pcolh2o_var_652(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pcolmol_var_653(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pcolo2_var_654(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pcolo3_var_655(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pforfac_var_656(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pforfrac_var_657(kidia_var_642 : kfdia_var_643, klev_var_644)
  INTEGER(KIND = 4), INTENT(INOUT) :: kindfor_var_658(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pselffac_var_659(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pselffrac_var_660(kidia_var_642 : kfdia_var_643, klev_var_644)
  INTEGER(KIND = 4), INTENT(INOUT) :: kindself_var_661(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pfac00_var_662(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pfac01_var_663(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pfac10_var_664(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(INOUT) :: pfac11_var_665(kidia_var_642 : kfdia_var_643, klev_var_644)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjp_var_666(kidia_var_642 : kfdia_var_643, klev_var_644)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjt_var_667(kidia_var_642 : kfdia_var_643, klev_var_644)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjt1_var_668(kidia_var_642 : kfdia_var_643, klev_var_644)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_669(kidia_var_642 : kfdia_var_643)
  INTEGER(KIND = 4) :: i_nlayers_var_670, jk_var_671, jl_var_672, jp1_var_673
  REAL(KIND = 8) :: z_stpfac_var_674, z_plog_var_675
  REAL(KIND = 8) :: z_fp_var_676, z_ft_var_677, z_ft1_var_678, z_water_var_679, z_scalefac_var_680
  REAL(KIND = 8) :: z_factor_var_681, z_co2reg_var_682, z_compfp_var_683
  z_stpfac_var_674 = 0.29220138203356366D0
  i_nlayers_var_670 = klev_var_644
  DO jk_var_671 = 1, klev_var_644
    DO jl_var_672 = kidia_var_642, kfdia_var_643
      pcolmol_var_653(jl_var_672, jk_var_671) = 0.0D0
    END DO
  END DO
  DO jl_var_672 = kidia_var_642, kfdia_var_643
    IF (prmu0_var_669(jl_var_672) > 0.0D0) THEN
      klaytrop_var_649(jl_var_672) = 0
    END IF
  END DO
  DO jk_var_671 = 1, i_nlayers_var_670
    DO jl_var_672 = kidia_var_642, kfdia_var_643
      IF (prmu0_var_669(jl_var_672) > 0.0D0) THEN
        z_plog_var_675 = LOG(pavel_var_645(jl_var_672, jk_var_671))
        kjp_var_666(jl_var_672, jk_var_671) = INT(36.0D0 - 5.0D0 * (z_plog_var_675 + 0.04D0))
        IF (kjp_var_666(jl_var_672, jk_var_671) < 1) THEN
          kjp_var_666(jl_var_672, jk_var_671) = 1
        ELSE IF (kjp_var_666(jl_var_672, jk_var_671) > 58) THEN
          kjp_var_666(jl_var_672, jk_var_671) = 58
        END IF
        jp1_var_673 = kjp_var_666(jl_var_672, jk_var_671) + 1
        z_fp_var_676 = 5.0 * (preflog_var_335(kjp_var_666(jl_var_672, jk_var_671)) - z_plog_var_675)
        kjt_var_667(jl_var_672, jk_var_671) = INT(3.0 + (ptavel_var_646(jl_var_672, jk_var_671) - tref_var_336(kjp_var_666(jl_var_672, jk_var_671))) / 15.0)
        IF (kjt_var_667(jl_var_672, jk_var_671) < 1) THEN
          kjt_var_667(jl_var_672, jk_var_671) = 1
        ELSE IF (kjt_var_667(jl_var_672, jk_var_671) > 4) THEN
          kjt_var_667(jl_var_672, jk_var_671) = 4
        END IF
        z_ft_var_677 = ((ptavel_var_646(jl_var_672, jk_var_671) - tref_var_336(kjp_var_666(jl_var_672, jk_var_671))) / 15.0) - REAL(kjt_var_667(jl_var_672, jk_var_671) - 3)
        kjt1_var_668(jl_var_672, jk_var_671) = INT(3.0 + (ptavel_var_646(jl_var_672, jk_var_671) - tref_var_336(jp1_var_673)) / 15.0)
        IF (kjt1_var_668(jl_var_672, jk_var_671) < 1) THEN
          kjt1_var_668(jl_var_672, jk_var_671) = 1
        ELSE IF (kjt1_var_668(jl_var_672, jk_var_671) > 4) THEN
          kjt1_var_668(jl_var_672, jk_var_671) = 4
        END IF
        z_ft1_var_678 = ((ptavel_var_646(jl_var_672, jk_var_671) - tref_var_336(jp1_var_673)) / 15.0) - REAL(kjt1_var_668(jl_var_672, jk_var_671) - 3)
        z_water_var_679 = pwkl_var_648(jl_var_672, 1, jk_var_671) / pcoldry_var_647(jl_var_672, jk_var_671)
        z_scalefac_var_680 = pavel_var_645(jl_var_672, jk_var_671) * z_stpfac_var_674 / ptavel_var_646(jl_var_672, jk_var_671)
        IF (z_plog_var_675 <= 4.56D0) GO TO 5300
        klaytrop_var_649(jl_var_672) = klaytrop_var_649(jl_var_672) + 1
        pforfac_var_656(jl_var_672, jk_var_671) = z_scalefac_var_680 / (1.0 + z_water_var_679)
        z_factor_var_681 = (332.0 - ptavel_var_646(jl_var_672, jk_var_671)) / 36.0
        kindfor_var_658(jl_var_672, jk_var_671) = MIN(2, MAX(1, INT(z_factor_var_681)))
        pforfrac_var_657(jl_var_672, jk_var_671) = z_factor_var_681 - REAL(kindfor_var_658(jl_var_672, jk_var_671))
        pselffac_var_659(jl_var_672, jk_var_671) = z_water_var_679 * pforfac_var_656(jl_var_672, jk_var_671)
        z_factor_var_681 = (ptavel_var_646(jl_var_672, jk_var_671) - 188.0) / 7.2
        kindself_var_661(jl_var_672, jk_var_671) = MIN(9, MAX(1, INT(z_factor_var_681) - 7))
        pselffrac_var_660(jl_var_672, jk_var_671) = z_factor_var_681 - REAL(kindself_var_661(jl_var_672, jk_var_671) + 7)
        pcolh2o_var_652(jl_var_672, jk_var_671) = 1E-20 * pwkl_var_648(jl_var_672, 1, jk_var_671)
        pcolco2_var_651(jl_var_672, jk_var_671) = 1E-20 * pwkl_var_648(jl_var_672, 2, jk_var_671)
        pcolo3_var_655(jl_var_672, jk_var_671) = 1E-20 * pwkl_var_648(jl_var_672, 3, jk_var_671)
        pcolch4_var_650(jl_var_672, jk_var_671) = 1E-20 * pwkl_var_648(jl_var_672, 6, jk_var_671)
        pcolo2_var_654(jl_var_672, jk_var_671) = 1E-20 * pwkl_var_648(jl_var_672, 7, jk_var_671)
        pcolmol_var_653(jl_var_672, jk_var_671) = 1E-20 * pcoldry_var_647(jl_var_672, jk_var_671) + pcolh2o_var_652(jl_var_672, jk_var_671)
        IF (pcolco2_var_651(jl_var_672, jk_var_671) == 0.0) pcolco2_var_651(jl_var_672, jk_var_671) = 1E-32 * pcoldry_var_647(jl_var_672, jk_var_671)
        IF (pcolch4_var_650(jl_var_672, jk_var_671) == 0.0) pcolch4_var_650(jl_var_672, jk_var_671) = 1E-32 * pcoldry_var_647(jl_var_672, jk_var_671)
        IF (pcolo2_var_654(jl_var_672, jk_var_671) == 0.0) pcolo2_var_654(jl_var_672, jk_var_671) = 1E-32 * pcoldry_var_647(jl_var_672, jk_var_671)
        z_co2reg_var_682 = 3.55E-24 * pcoldry_var_647(jl_var_672, jk_var_671)
        GO TO 5400
5300    CONTINUE
        pforfac_var_656(jl_var_672, jk_var_671) = z_scalefac_var_680 / (1.0 + z_water_var_679)
        z_factor_var_681 = (ptavel_var_646(jl_var_672, jk_var_671) - 188.0) / 36.0
        kindfor_var_658(jl_var_672, jk_var_671) = 3
        pforfrac_var_657(jl_var_672, jk_var_671) = z_factor_var_681 - 1.0
        pcolh2o_var_652(jl_var_672, jk_var_671) = 1E-20 * pwkl_var_648(jl_var_672, 1, jk_var_671)
        pcolco2_var_651(jl_var_672, jk_var_671) = 1E-20 * pwkl_var_648(jl_var_672, 2, jk_var_671)
        pcolo3_var_655(jl_var_672, jk_var_671) = 1E-20 * pwkl_var_648(jl_var_672, 3, jk_var_671)
        pcolch4_var_650(jl_var_672, jk_var_671) = 1E-20 * pwkl_var_648(jl_var_672, 6, jk_var_671)
        pcolo2_var_654(jl_var_672, jk_var_671) = 1E-20 * pwkl_var_648(jl_var_672, 7, jk_var_671)
        pcolmol_var_653(jl_var_672, jk_var_671) = 1E-20 * pcoldry_var_647(jl_var_672, jk_var_671) + pcolh2o_var_652(jl_var_672, jk_var_671)
        IF (pcolco2_var_651(jl_var_672, jk_var_671) == 0.0) pcolco2_var_651(jl_var_672, jk_var_671) = 1E-32 * pcoldry_var_647(jl_var_672, jk_var_671)
        IF (pcolch4_var_650(jl_var_672, jk_var_671) == 0.0) pcolch4_var_650(jl_var_672, jk_var_671) = 1E-32 * pcoldry_var_647(jl_var_672, jk_var_671)
        IF (pcolo2_var_654(jl_var_672, jk_var_671) == 0.0) pcolo2_var_654(jl_var_672, jk_var_671) = 1E-32 * pcoldry_var_647(jl_var_672, jk_var_671)
        z_co2reg_var_682 = 3.55E-24 * pcoldry_var_647(jl_var_672, jk_var_671)
        pselffac_var_659(jl_var_672, jk_var_671) = 0.0D0
        pselffrac_var_660(jl_var_672, jk_var_671) = 0.0D0
        kindself_var_661(jl_var_672, jk_var_671) = 0
5400    CONTINUE
        z_compfp_var_683 = 1.0 - z_fp_var_676
        pfac10_var_664(jl_var_672, jk_var_671) = z_compfp_var_683 * z_ft_var_677
        pfac00_var_662(jl_var_672, jk_var_671) = z_compfp_var_683 * (1.0 - z_ft_var_677)
        pfac11_var_665(jl_var_672, jk_var_671) = z_fp_var_676 * z_ft1_var_678
        pfac01_var_663(jl_var_672, jk_var_671) = z_fp_var_676 * (1.0 - z_ft1_var_678)
9000    FORMAT(1X, 2I3, 3I4, F6.1, 4F7.2, 12E9.2, 2I5)
      END IF
    END DO
  END DO
END SUBROUTINE srtm_setcoef
SUBROUTINE rrtm_gas_optical_depth(kidia_var_684, kfdia_var_685, klev_var_686, pod_var_687, pavel_var_688, pcoldry_var_689, pcolbrd, pwx_var_690, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, poneminus_var_698, pcolh2o_var_699, pcolco2_var_700, pcolo3_var_701, pcoln2o, pcolch4_var_702, pcolo2_var_703, p_co2mult_var_704, klaytrop_var_705, klayswtch, klaylow, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, kindminor, pscaleminor, pscaleminorn2, pminorfrac, prat_h2oco2_var_713, prat_h2oco2_1_var_714, prat_h2oo3_var_715, prat_h2oo3_1_var_716, prat_h2on2o_var_717, prat_h2on2o_1_var_718, prat_h2och4_var_719, prat_h2och4_1_var_720, prat_n2oco2_var_721, prat_n2oco2_1_var_722, prat_o3co2_var_723, prat_o3co2_1_var_724)
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_684
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_685
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_686
  REAL(KIND = 8), INTENT(OUT) :: pod_var_687(140, klev_var_686, kidia_var_684 : kfdia_var_685)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_688(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pcoldry_var_689(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pwx_var_690(kidia_var_684 : kfdia_var_685, 4, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: ptauaerl(kidia_var_684 : kfdia_var_685, klev_var_686, 16)
  REAL(KIND = 8), INTENT(IN) :: pfac00_var_691(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pfac01_var_692(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pfac10_var_693(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pfac11_var_694(kidia_var_684 : kfdia_var_685, klev_var_686)
  INTEGER(KIND = 4), INTENT(IN) :: kjp_var_695(kidia_var_684 : kfdia_var_685, klev_var_686)
  INTEGER(KIND = 4), INTENT(IN) :: kjt_var_696(kidia_var_684 : kfdia_var_685, klev_var_686)
  INTEGER(KIND = 4), INTENT(IN) :: kjt1_var_697(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: poneminus_var_698
  REAL(KIND = 8), INTENT(IN) :: pcolh2o_var_699(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pcolco2_var_700(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pcolo3_var_701(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pcoln2o(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pcolch4_var_702(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pcolo2_var_703(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: p_co2mult_var_704(kidia_var_684 : kfdia_var_685, klev_var_686)
  INTEGER(KIND = 4), INTENT(IN) :: klaytrop_var_705(kidia_var_684 : kfdia_var_685)
  INTEGER(KIND = 4), INTENT(IN) :: klayswtch(kidia_var_684 : kfdia_var_685)
  INTEGER(KIND = 4), INTENT(IN) :: klaylow(kidia_var_684 : kfdia_var_685)
  REAL(KIND = 8), INTENT(IN) :: pselffac_var_706(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pselffrac_var_707(kidia_var_684 : kfdia_var_685, klev_var_686)
  INTEGER(KIND = 4), INTENT(IN) :: kindself_var_708(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(OUT) :: pfrac_var_709(kidia_var_684 : kfdia_var_685, 140, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pforfac_var_710(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pforfrac_var_711(kidia_var_684 : kfdia_var_685, klev_var_686)
  INTEGER(KIND = 4), INTENT(IN) :: kindfor_var_712(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pminorfrac(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pscaleminor(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pscaleminorn2(kidia_var_684 : kfdia_var_685, klev_var_686)
  INTEGER(KIND = 4), INTENT(IN) :: kindminor(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: pcolbrd(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8), INTENT(IN) :: prat_h2oco2_var_713(kidia_var_684 : kfdia_var_685, klev_var_686), prat_h2oco2_1_var_714(kidia_var_684 : kfdia_var_685, klev_var_686), prat_h2oo3_var_715(kidia_var_684 : kfdia_var_685, klev_var_686), prat_h2oo3_1_var_716(kidia_var_684 : kfdia_var_685, klev_var_686), prat_h2on2o_var_717(kidia_var_684 : kfdia_var_685, klev_var_686), prat_h2on2o_1_var_718(kidia_var_684 : kfdia_var_685, klev_var_686), prat_h2och4_var_719(kidia_var_684 : kfdia_var_685, klev_var_686), prat_h2och4_1_var_720(kidia_var_684 : kfdia_var_685, klev_var_686), prat_n2oco2_var_721(kidia_var_684 : kfdia_var_685, klev_var_686), prat_n2oco2_1_var_722(kidia_var_684 : kfdia_var_685, klev_var_686), prat_o3co2_var_723(kidia_var_684 : kfdia_var_685, klev_var_686), prat_o3co2_1_var_724(kidia_var_684 : kfdia_var_685, klev_var_686)
  REAL(KIND = 8) :: ztau(kidia_var_684 : kfdia_var_685, 140, klev_var_686)
  INTEGER(KIND = 4) :: ji, jlev_var_725
  INTEGER(KIND = 4) :: jlon_var_726
  pfrac_var_709(:, :, :) = 0.0D0
  CALL rrtm_taumol1(kidia_var_684, kfdia_var_685, klev_var_686, ztau, pavel_var_688, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, pcolh2o_var_699, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, pminorfrac, kindminor, pscaleminorn2, pcolbrd)
  CALL rrtm_taumol2(kidia_var_684, kfdia_var_685, klev_var_686, ztau, pavel_var_688, pcoldry_var_689, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, pcolh2o_var_699, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709)
  CALL rrtm_taumol3(kidia_var_684, kfdia_var_685, klev_var_686, ztau, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, poneminus_var_698, pcolh2o_var_699, pcolco2_var_700, pcoln2o, pcoldry_var_689, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, prat_h2oco2_var_713, prat_h2oco2_1_var_714, pminorfrac, kindminor)
  CALL rrtm_taumol4(kidia_var_684, kfdia_var_685, klev_var_686, ztau, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, poneminus_var_698, pcolh2o_var_699, pcolco2_var_700, pcolo3_var_701, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, prat_h2oco2_var_713, prat_h2oco2_1_var_714, prat_o3co2_var_723, prat_o3co2_1_var_724)
  CALL rrtm_taumol5(kidia_var_684, kfdia_var_685, klev_var_686, ztau, pwx_var_690, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, poneminus_var_698, pcolh2o_var_699, pcolco2_var_700, pcolo3_var_701, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, prat_h2oco2_var_713, prat_h2oco2_1_var_714, prat_o3co2_var_723, prat_o3co2_1_var_724, pminorfrac, kindminor)
  CALL rrtm_taumol6(kidia_var_684, kfdia_var_685, klev_var_686, ztau, pwx_var_690, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, pcolh2o_var_699, pcolco2_var_700, pcoldry_var_689, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, pminorfrac, kindminor)
  CALL rrtm_taumol7(kidia_var_684, kfdia_var_685, klev_var_686, ztau, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, poneminus_var_698, pcolh2o_var_699, pcolo3_var_701, pcolco2_var_700, pcoldry_var_689, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, prat_h2oo3_var_715, prat_h2oo3_1_var_716, pminorfrac, kindminor)
  CALL rrtm_taumol8(kidia_var_684, kfdia_var_685, klev_var_686, ztau, pwx_var_690, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, pcolh2o_var_699, pcolo3_var_701, pcoln2o, pcolco2_var_700, pcoldry_var_689, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, pminorfrac, kindminor)
  CALL rrtm_taumol9(kidia_var_684, kfdia_var_685, klev_var_686, ztau, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, poneminus_var_698, pcolh2o_var_699, pcoln2o, pcolch4_var_702, pcoldry_var_689, klaytrop_var_705, klayswtch, klaylow, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, prat_h2och4_var_719, prat_h2och4_1_var_720, pminorfrac, kindminor)
  CALL rrtm_taumol10(kidia_var_684, kfdia_var_685, klev_var_686, ztau, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, pcolh2o_var_699, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709)
  CALL rrtm_taumol11(kidia_var_684, kfdia_var_685, klev_var_686, ztau, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, pcolh2o_var_699, pcolo2_var_703, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, pminorfrac, kindminor, pscaleminor)
  CALL rrtm_taumol12(kidia_var_684, kfdia_var_685, klev_var_686, ztau, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, poneminus_var_698, pcolh2o_var_699, pcolco2_var_700, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, prat_h2oco2_var_713, prat_h2oco2_1_var_714)
  CALL rrtm_taumol13(kidia_var_684, kfdia_var_685, klev_var_686, ztau, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, poneminus_var_698, pcolh2o_var_699, pcoln2o, pcolco2_var_700, pcolo3_var_701, pcoldry_var_689, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, prat_h2on2o_var_717, prat_h2on2o_1_var_718, pminorfrac, kindminor)
  CALL rrtm_taumol14(kidia_var_684, kfdia_var_685, klev_var_686, ztau, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, pcolco2_var_700, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709)
  CALL rrtm_taumol15(kidia_var_684, kfdia_var_685, klev_var_686, ztau, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, poneminus_var_698, pcolh2o_var_699, pcolco2_var_700, pcoln2o, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, prat_n2oco2_var_721, prat_n2oco2_1_var_722, pminorfrac, kindminor, pscaleminor, pcolbrd)
  CALL rrtm_taumol16(kidia_var_684, kfdia_var_685, klev_var_686, ztau, ptauaerl, pfac00_var_691, pfac01_var_692, pfac10_var_693, pfac11_var_694, pforfac_var_710, pforfrac_var_711, kindfor_var_712, kjp_var_695, kjt_var_696, kjt1_var_697, poneminus_var_698, pcolh2o_var_699, pcolch4_var_702, klaytrop_var_705, pselffac_var_706, pselffrac_var_707, kindself_var_708, pfrac_var_709, prat_h2och4_var_719, prat_h2och4_1_var_720)
  DO jlev_var_725 = 1, klev_var_686
    DO ji = 1, 140
      DO jlon_var_726 = kidia_var_684, kfdia_var_685
        pod_var_687(ji, jlev_var_725, jlon_var_726) = ztau(jlon_var_726, ji, jlev_var_725)
      END DO
    END DO
  END DO
END SUBROUTINE rrtm_gas_optical_depth
SUBROUTINE rrtm_taumol9(kidia_var_727, kfdia_var_728, klev_var_729, taug_var_730, p_tauaerl_var_731, fac00_var_732, fac01_var_733, fac10_var_734, fac11_var_735, forfac_var_754, forfrac_var_755, indfor_var_753, jp_var_736, jt_var_737, jt1_var_738, oneminus_var_739, colh2o_var_740, coln2o_var_741, colch4_var_742, coldry_var_743, laytrop_var_744, k_layswtch_var_745, k_laylow_var_746, selffac_var_747, selffrac_var_748, indself_var_749, fracs_var_750, rat_h2och4_var_751, rat_h2och4_1_var_752, minorfrac_var_756, indminor_var_757)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrtm, ONLY: ng9
  USE yoerrta9, ONLY: absa_var_228, absb_var_229, forref_var_233, fracrefa_var_226, fracrefb_var_227, ka_mn2o_var_230, kb_mn2o_var_231, selfref_var_232
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_727
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_728
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_729
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_730(kidia_var_727 : kfdia_var_728, 140, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_731(kidia_var_727 : kfdia_var_728, klev_var_729, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_732(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_733(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_734(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_735(kidia_var_727 : kfdia_var_728, klev_var_729)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_736(kidia_var_727 : kfdia_var_728, klev_var_729)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_737(kidia_var_727 : kfdia_var_728, klev_var_729)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_738(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_739
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_740(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_741(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: colch4_var_742(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_743(kidia_var_727 : kfdia_var_728, klev_var_729)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_744(kidia_var_727 : kfdia_var_728)
  INTEGER(KIND = 4), INTENT(IN) :: k_layswtch_var_745(kidia_var_727 : kfdia_var_728)
  INTEGER(KIND = 4), INTENT(IN) :: k_laylow_var_746(kidia_var_727 : kfdia_var_728)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_747(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_748(kidia_var_727 : kfdia_var_728, klev_var_729)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_749(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_750(kidia_var_727 : kfdia_var_728, 140, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_var_751(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_1_var_752(kidia_var_727 : kfdia_var_728, klev_var_729)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_753(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_754(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_755(kidia_var_727 : kfdia_var_728, klev_var_729)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_756(kidia_var_727 : kfdia_var_728, klev_var_729)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_757(kidia_var_727 : kfdia_var_728, klev_var_729)
  INTEGER(KIND = 4) :: ind0_var_758, ind1_var_759, inds_var_760, indf_var_761, indm_var_762
  INTEGER(KIND = 4) :: ig_var_763, js_var_764, lay_var_765, js1_var_766, jmn2o_var_767, jpl_var_768
  REAL(KIND = 8) :: speccomb_var_769, speccomb1_var_770, speccomb_mn2o_var_771, speccomb_planck_var_772
  REAL(KIND = 8) :: refrat_planck_a_var_773, refrat_m_a_var_774
  REAL(KIND = 8) :: fs_var_775, specmult_var_776, specparm_var_777, fs1_var_778, specmult1_var_779, specparm1_var_780, fmn2o_var_781, specmult_mn2o_var_782, specparm_mn2o_var_783, fpl_var_784, specmult_planck_var_785, specparm_planck_var_786
  REAL(KIND = 8) :: adjfac_var_787, adjcoln2o_var_788, ratn2o_var_789, chi_n2o_var_790
  REAL(KIND = 8) :: fac000_var_791, fac100_var_792, fac200_var_793, fac010_var_794, fac110_var_795, fac210_var_796, fac001_var_797, fac101_var_798, fac201_var_799, fac011_var_800, fac111_var_801, fac211_var_802
  REAL(KIND = 8) :: p_var_803, p4_var_804, fk0_var_805, fk1_var_806, fk2_var_807
  REAL(KIND = 8) :: taufor_var_808, tauself_var_809, n2om1_var_810, n2om2_var_811, absn2o_var_812, tau_major_var_813(12), tau_major1_var_814(12)
  INTEGER(KIND = 4) :: laytrop_min_var_815, laytrop_max_var_816
  INTEGER(KIND = 4) :: ixc_var_817(klev_var_729), ixlow_var_818(kfdia_var_728, klev_var_729), ixhigh_var_819(kfdia_var_728, klev_var_729)
  INTEGER(KIND = 4) :: ich_var_820, icl_var_821, ixc0_var_822, ixp_var_823, jc_var_824, jl_var_825
  laytrop_min_var_815 = MINVAL(laytrop_var_744)
  laytrop_max_var_816 = MAXVAL(laytrop_var_744)
  ixlow_var_818 = 0
  ixhigh_var_819 = 0
  ixc_var_817 = 0
  DO lay_var_765 = laytrop_min_var_815 + 1, laytrop_max_var_816
    icl_var_821 = 0
    ich_var_820 = 0
    DO jc_var_824 = kidia_var_727, kfdia_var_728
      IF (lay_var_765 <= laytrop_var_744(jc_var_824)) THEN
        icl_var_821 = icl_var_821 + 1
        ixlow_var_818(icl_var_821, lay_var_765) = jc_var_824
      ELSE
        ich_var_820 = ich_var_820 + 1
        ixhigh_var_819(ich_var_820, lay_var_765) = jc_var_824
      END IF
    END DO
    ixc_var_817(lay_var_765) = icl_var_821
  END DO
  refrat_planck_a_var_773 = chi_mls(1, 9) / chi_mls(6, 9)
  refrat_m_a_var_774 = chi_mls(1, 3) / chi_mls(6, 3)
  DO lay_var_765 = 1, laytrop_min_var_815
    DO jl_var_825 = kidia_var_727, kfdia_var_728
      speccomb_var_769 = colh2o_var_740(jl_var_825, lay_var_765) + rat_h2och4_var_751(jl_var_825, lay_var_765) * colch4_var_742(jl_var_825, lay_var_765)
      specparm_var_777 = MIN(colh2o_var_740(jl_var_825, lay_var_765) / speccomb_var_769, oneminus_var_739)
      specmult_var_776 = 8.0D0 * (specparm_var_777)
      js_var_764 = 1 + INT(specmult_var_776)
      fs_var_775 = ((specmult_var_776) - AINT((specmult_var_776)))
      speccomb1_var_770 = colh2o_var_740(jl_var_825, lay_var_765) + rat_h2och4_1_var_752(jl_var_825, lay_var_765) * colch4_var_742(jl_var_825, lay_var_765)
      specparm1_var_780 = MIN(colh2o_var_740(jl_var_825, lay_var_765) / speccomb1_var_770, oneminus_var_739)
      specmult1_var_779 = 8.0D0 * (specparm1_var_780)
      js1_var_766 = 1 + INT(specmult1_var_779)
      fs1_var_778 = ((specmult1_var_779) - AINT((specmult1_var_779)))
      speccomb_mn2o_var_771 = colh2o_var_740(jl_var_825, lay_var_765) + refrat_m_a_var_774 * colch4_var_742(jl_var_825, lay_var_765)
      specparm_mn2o_var_783 = MIN(colh2o_var_740(jl_var_825, lay_var_765) / speccomb_mn2o_var_771, oneminus_var_739)
      specmult_mn2o_var_782 = 8.0D0 * specparm_mn2o_var_783
      jmn2o_var_767 = 1 + INT(specmult_mn2o_var_782)
      fmn2o_var_781 = ((specmult_mn2o_var_782) - AINT((specmult_mn2o_var_782)))
      chi_n2o_var_790 = coln2o_var_741(jl_var_825, lay_var_765) / (coldry_var_743(jl_var_825, lay_var_765))
      ratn2o_var_789 = 1D+20 * chi_n2o_var_790 / chi_mls(4, jp_var_736(jl_var_825, lay_var_765) + 1)
      IF (ratn2o_var_789 .GT. 1.5D0) THEN
        adjfac_var_787 = 0.5D0 + (ratn2o_var_789 - 0.5D0) ** 0.65D0
        adjcoln2o_var_788 = adjfac_var_787 * chi_mls(4, jp_var_736(jl_var_825, lay_var_765) + 1) * coldry_var_743(jl_var_825, lay_var_765) * 1D-20
      ELSE
        adjcoln2o_var_788 = coln2o_var_741(jl_var_825, lay_var_765)
      END IF
      speccomb_planck_var_772 = colh2o_var_740(jl_var_825, lay_var_765) + refrat_planck_a_var_773 * colch4_var_742(jl_var_825, lay_var_765)
      specparm_planck_var_786 = MIN(colh2o_var_740(jl_var_825, lay_var_765) / speccomb_planck_var_772, oneminus_var_739)
      specmult_planck_var_785 = 8.0D0 * specparm_planck_var_786
      jpl_var_768 = 1 + INT(specmult_planck_var_785)
      fpl_var_784 = ((specmult_planck_var_785) - AINT((specmult_planck_var_785)))
      ind0_var_758 = ((jp_var_736(jl_var_825, lay_var_765) - 1) * 5 + (jt_var_737(jl_var_825, lay_var_765) - 1)) * nspa_var_236(9) + js_var_764
      ind1_var_759 = (jp_var_736(jl_var_825, lay_var_765) * 5 + (jt1_var_738(jl_var_825, lay_var_765) - 1)) * nspa_var_236(9) + js1_var_766
      inds_var_760 = indself_var_749(jl_var_825, lay_var_765)
      indf_var_761 = indfor_var_753(jl_var_825, lay_var_765)
      indm_var_762 = indminor_var_757(jl_var_825, lay_var_765)
      IF (specparm_var_777 .LT. 0.125D0) THEN
        p_var_803 = fs_var_775 - 1.0D0
        p4_var_804 = p_var_803 ** 4
        fk0_var_805 = p4_var_804
        fk1_var_806 = 1.0D0 - p_var_803 - 2.0D0 * p4_var_804
        fk2_var_807 = p_var_803 + p4_var_804
        fac000_var_791 = fk0_var_805 * fac00_var_732(jl_var_825, lay_var_765)
        fac100_var_792 = fk1_var_806 * fac00_var_732(jl_var_825, lay_var_765)
        fac200_var_793 = fk2_var_807 * fac00_var_732(jl_var_825, lay_var_765)
        fac010_var_794 = fk0_var_805 * fac10_var_734(jl_var_825, lay_var_765)
        fac110_var_795 = fk1_var_806 * fac10_var_734(jl_var_825, lay_var_765)
        fac210_var_796 = fk2_var_807 * fac10_var_734(jl_var_825, lay_var_765)
      ELSE IF (specparm_var_777 .GT. 0.875D0) THEN
        p_var_803 = - fs_var_775
        p4_var_804 = p_var_803 ** 4
        fk0_var_805 = p4_var_804
        fk1_var_806 = 1.0D0 - p_var_803 - 2.0D0 * p4_var_804
        fk2_var_807 = p_var_803 + p4_var_804
        fac000_var_791 = fk0_var_805 * fac00_var_732(jl_var_825, lay_var_765)
        fac100_var_792 = fk1_var_806 * fac00_var_732(jl_var_825, lay_var_765)
        fac200_var_793 = fk2_var_807 * fac00_var_732(jl_var_825, lay_var_765)
        fac010_var_794 = fk0_var_805 * fac10_var_734(jl_var_825, lay_var_765)
        fac110_var_795 = fk1_var_806 * fac10_var_734(jl_var_825, lay_var_765)
        fac210_var_796 = fk2_var_807 * fac10_var_734(jl_var_825, lay_var_765)
      ELSE
        fac000_var_791 = (1.0D0 - fs_var_775) * fac00_var_732(jl_var_825, lay_var_765)
        fac010_var_794 = (1.0D0 - fs_var_775) * fac10_var_734(jl_var_825, lay_var_765)
        fac100_var_792 = fs_var_775 * fac00_var_732(jl_var_825, lay_var_765)
        fac110_var_795 = fs_var_775 * fac10_var_734(jl_var_825, lay_var_765)
        fac200_var_793 = 0.0D0
        fac210_var_796 = 0.0D0
      END IF
      IF (specparm1_var_780 .LT. 0.125D0) THEN
        p_var_803 = fs1_var_778 - 1.0D0
        p4_var_804 = p_var_803 ** 4
        fk0_var_805 = p4_var_804
        fk1_var_806 = 1.0D0 - p_var_803 - 2.0D0 * p4_var_804
        fk2_var_807 = p_var_803 + p4_var_804
        fac001_var_797 = fk0_var_805 * fac01_var_733(jl_var_825, lay_var_765)
        fac101_var_798 = fk1_var_806 * fac01_var_733(jl_var_825, lay_var_765)
        fac201_var_799 = fk2_var_807 * fac01_var_733(jl_var_825, lay_var_765)
        fac011_var_800 = fk0_var_805 * fac11_var_735(jl_var_825, lay_var_765)
        fac111_var_801 = fk1_var_806 * fac11_var_735(jl_var_825, lay_var_765)
        fac211_var_802 = fk2_var_807 * fac11_var_735(jl_var_825, lay_var_765)
      ELSE IF (specparm1_var_780 .GT. 0.875D0) THEN
        p_var_803 = - fs1_var_778
        p4_var_804 = p_var_803 ** 4
        fk0_var_805 = p4_var_804
        fk1_var_806 = 1.0D0 - p_var_803 - 2.0D0 * p4_var_804
        fk2_var_807 = p_var_803 + p4_var_804
        fac001_var_797 = fk0_var_805 * fac01_var_733(jl_var_825, lay_var_765)
        fac101_var_798 = fk1_var_806 * fac01_var_733(jl_var_825, lay_var_765)
        fac201_var_799 = fk2_var_807 * fac01_var_733(jl_var_825, lay_var_765)
        fac011_var_800 = fk0_var_805 * fac11_var_735(jl_var_825, lay_var_765)
        fac111_var_801 = fk1_var_806 * fac11_var_735(jl_var_825, lay_var_765)
        fac211_var_802 = fk2_var_807 * fac11_var_735(jl_var_825, lay_var_765)
      ELSE
        fac001_var_797 = (1.0D0 - fs1_var_778) * fac01_var_733(jl_var_825, lay_var_765)
        fac011_var_800 = (1.0D0 - fs1_var_778) * fac11_var_735(jl_var_825, lay_var_765)
        fac101_var_798 = fs1_var_778 * fac01_var_733(jl_var_825, lay_var_765)
        fac111_var_801 = fs1_var_778 * fac11_var_735(jl_var_825, lay_var_765)
        fac201_var_799 = 0.0D0
        fac211_var_802 = 0.0D0
      END IF
      IF (specparm_var_777 .LT. 0.125D0) THEN
        tau_major_var_813(1 : ng9) = speccomb_var_769 * (fac000_var_791 * absa_var_228(ind0_var_758, 1 : 12) + fac100_var_792 * absa_var_228(ind0_var_758 + 1, 1 : 12) + fac200_var_793 * absa_var_228(ind0_var_758 + 2, 1 : 12) + fac010_var_794 * absa_var_228(ind0_var_758 + 9, 1 : 12) + fac110_var_795 * absa_var_228(ind0_var_758 + 10, 1 : 12) + fac210_var_796 * absa_var_228(ind0_var_758 + 11, 1 : 12))
      ELSE IF (specparm_var_777 .GT. 0.875D0) THEN
        tau_major_var_813(1 : ng9) = speccomb_var_769 * (fac200_var_793 * absa_var_228(ind0_var_758 - 1, 1 : 12) + fac100_var_792 * absa_var_228(ind0_var_758, 1 : 12) + fac000_var_791 * absa_var_228(ind0_var_758 + 1, 1 : 12) + fac210_var_796 * absa_var_228(ind0_var_758 + 8, 1 : 12) + fac110_var_795 * absa_var_228(ind0_var_758 + 9, 1 : 12) + fac010_var_794 * absa_var_228(ind0_var_758 + 10, 1 : 12))
      ELSE
        tau_major_var_813(1 : ng9) = speccomb_var_769 * (fac000_var_791 * absa_var_228(ind0_var_758, 1 : 12) + fac100_var_792 * absa_var_228(ind0_var_758 + 1, 1 : 12) + fac010_var_794 * absa_var_228(ind0_var_758 + 9, 1 : 12) + fac110_var_795 * absa_var_228(ind0_var_758 + 10, 1 : 12))
      END IF
      IF (specparm1_var_780 .LT. 0.125D0) THEN
        tau_major1_var_814(1 : ng9) = speccomb1_var_770 * (fac001_var_797 * absa_var_228(ind1_var_759, 1 : 12) + fac101_var_798 * absa_var_228(ind1_var_759 + 1, 1 : 12) + fac201_var_799 * absa_var_228(ind1_var_759 + 2, 1 : 12) + fac011_var_800 * absa_var_228(ind1_var_759 + 9, 1 : 12) + fac111_var_801 * absa_var_228(ind1_var_759 + 10, 1 : 12) + fac211_var_802 * absa_var_228(ind1_var_759 + 11, 1 : 12))
      ELSE IF (specparm1_var_780 .GT. 0.875D0) THEN
        tau_major1_var_814(1 : ng9) = speccomb1_var_770 * (fac201_var_799 * absa_var_228(ind1_var_759 - 1, 1 : 12) + fac101_var_798 * absa_var_228(ind1_var_759, 1 : 12) + fac001_var_797 * absa_var_228(ind1_var_759 + 1, 1 : 12) + fac211_var_802 * absa_var_228(ind1_var_759 + 8, 1 : 12) + fac111_var_801 * absa_var_228(ind1_var_759 + 9, 1 : 12) + fac011_var_800 * absa_var_228(ind1_var_759 + 10, 1 : 12))
      ELSE
        tau_major1_var_814(1 : ng9) = speccomb1_var_770 * (fac001_var_797 * absa_var_228(ind1_var_759, 1 : 12) + fac101_var_798 * absa_var_228(ind1_var_759 + 1, 1 : 12) + fac011_var_800 * absa_var_228(ind1_var_759 + 9, 1 : 12) + fac111_var_801 * absa_var_228(ind1_var_759 + 10, 1 : 12))
      END IF
      DO ig_var_763 = 1, 12
        tauself_var_809 = selffac_var_747(jl_var_825, lay_var_765) * (selfref_var_232(inds_var_760, ig_var_763) + selffrac_var_748(jl_var_825, lay_var_765) * (selfref_var_232(inds_var_760 + 1, ig_var_763) - selfref_var_232(inds_var_760, ig_var_763)))
        taufor_var_808 = forfac_var_754(jl_var_825, lay_var_765) * (forref_var_233(indf_var_761, ig_var_763) + forfrac_var_755(jl_var_825, lay_var_765) * (forref_var_233(indf_var_761 + 1, ig_var_763) - forref_var_233(indf_var_761, ig_var_763)))
        n2om1_var_810 = ka_mn2o_var_230(jmn2o_var_767, indm_var_762, ig_var_763) + fmn2o_var_781 * (ka_mn2o_var_230(jmn2o_var_767 + 1, indm_var_762, ig_var_763) - ka_mn2o_var_230(jmn2o_var_767, indm_var_762, ig_var_763))
        n2om2_var_811 = ka_mn2o_var_230(jmn2o_var_767, indm_var_762 + 1, ig_var_763) + fmn2o_var_781 * (ka_mn2o_var_230(jmn2o_var_767 + 1, indm_var_762 + 1, ig_var_763) - ka_mn2o_var_230(jmn2o_var_767, indm_var_762 + 1, ig_var_763))
        absn2o_var_812 = n2om1_var_810 + minorfrac_var_756(jl_var_825, lay_var_765) * (n2om2_var_811 - n2om1_var_810)
        taug_var_730(jl_var_825, 96 + ig_var_763, lay_var_765) = tau_major_var_813(ig_var_763) + tau_major1_var_814(ig_var_763) + tauself_var_809 + taufor_var_808 + adjcoln2o_var_788 * absn2o_var_812
        fracs_var_750(jl_var_825, 96 + ig_var_763, lay_var_765) = fracrefa_var_226(ig_var_763, jpl_var_768) + fpl_var_784 * (fracrefa_var_226(ig_var_763, jpl_var_768 + 1) - fracrefa_var_226(ig_var_763, jpl_var_768))
      END DO
    END DO
  END DO
  DO lay_var_765 = laytrop_max_var_816 + 1, klev_var_729
    DO jl_var_825 = kidia_var_727, kfdia_var_728
      chi_n2o_var_790 = coln2o_var_741(jl_var_825, lay_var_765) / (coldry_var_743(jl_var_825, lay_var_765))
      ratn2o_var_789 = 1D+20 * chi_n2o_var_790 / chi_mls(4, jp_var_736(jl_var_825, lay_var_765) + 1)
      IF (ratn2o_var_789 .GT. 1.5D0) THEN
        adjfac_var_787 = 0.5D0 + (ratn2o_var_789 - 0.5D0) ** 0.65D0
        adjcoln2o_var_788 = adjfac_var_787 * chi_mls(4, jp_var_736(jl_var_825, lay_var_765) + 1) * coldry_var_743(jl_var_825, lay_var_765) * 1D-20
      ELSE
        adjcoln2o_var_788 = coln2o_var_741(jl_var_825, lay_var_765)
      END IF
      ind0_var_758 = ((jp_var_736(jl_var_825, lay_var_765) - 13) * 5 + (jt_var_737(jl_var_825, lay_var_765) - 1)) * nspb_var_237(9) + 1
      ind1_var_759 = ((jp_var_736(jl_var_825, lay_var_765) - 12) * 5 + (jt1_var_738(jl_var_825, lay_var_765) - 1)) * nspb_var_237(9) + 1
      indm_var_762 = indminor_var_757(jl_var_825, lay_var_765)
      DO ig_var_763 = 1, 12
        absn2o_var_812 = kb_mn2o_var_231(indm_var_762, ig_var_763) + minorfrac_var_756(jl_var_825, lay_var_765) * (kb_mn2o_var_231(indm_var_762 + 1, ig_var_763) - kb_mn2o_var_231(indm_var_762, ig_var_763))
        taug_var_730(jl_var_825, 96 + ig_var_763, lay_var_765) = colch4_var_742(jl_var_825, lay_var_765) * (fac00_var_732(jl_var_825, lay_var_765) * absb_var_229(ind0_var_758, ig_var_763) + fac10_var_734(jl_var_825, lay_var_765) * absb_var_229(ind0_var_758 + 1, ig_var_763) + fac01_var_733(jl_var_825, lay_var_765) * absb_var_229(ind1_var_759, ig_var_763) + fac11_var_735(jl_var_825, lay_var_765) * absb_var_229(ind1_var_759 + 1, ig_var_763)) + adjcoln2o_var_788 * absn2o_var_812
        fracs_var_750(jl_var_825, 96 + ig_var_763, lay_var_765) = fracrefb_var_227(ig_var_763)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_816 /= laytrop_min_var_815) THEN
    DO lay_var_765 = laytrop_min_var_815 + 1, laytrop_max_var_816
      ixc0_var_822 = ixc_var_817(lay_var_765)
      DO ixp_var_823 = 1, ixc0_var_822
        jl_var_825 = ixlow_var_818(ixp_var_823, lay_var_765)
        speccomb_var_769 = colh2o_var_740(jl_var_825, lay_var_765) + rat_h2och4_var_751(jl_var_825, lay_var_765) * colch4_var_742(jl_var_825, lay_var_765)
        specparm_var_777 = MIN(colh2o_var_740(jl_var_825, lay_var_765) / speccomb_var_769, oneminus_var_739)
        specmult_var_776 = 8.0D0 * (specparm_var_777)
        js_var_764 = 1 + INT(specmult_var_776)
        fs_var_775 = ((specmult_var_776) - AINT((specmult_var_776)))
        speccomb1_var_770 = colh2o_var_740(jl_var_825, lay_var_765) + rat_h2och4_1_var_752(jl_var_825, lay_var_765) * colch4_var_742(jl_var_825, lay_var_765)
        specparm1_var_780 = MIN(colh2o_var_740(jl_var_825, lay_var_765) / speccomb1_var_770, oneminus_var_739)
        specmult1_var_779 = 8.0D0 * (specparm1_var_780)
        js1_var_766 = 1 + INT(specmult1_var_779)
        fs1_var_778 = ((specmult1_var_779) - AINT((specmult1_var_779)))
        speccomb_mn2o_var_771 = colh2o_var_740(jl_var_825, lay_var_765) + refrat_m_a_var_774 * colch4_var_742(jl_var_825, lay_var_765)
        specparm_mn2o_var_783 = MIN(colh2o_var_740(jl_var_825, lay_var_765) / speccomb_mn2o_var_771, oneminus_var_739)
        specmult_mn2o_var_782 = 8.0D0 * specparm_mn2o_var_783
        jmn2o_var_767 = 1 + INT(specmult_mn2o_var_782)
        fmn2o_var_781 = ((specmult_mn2o_var_782) - AINT((specmult_mn2o_var_782)))
        chi_n2o_var_790 = coln2o_var_741(jl_var_825, lay_var_765) / (coldry_var_743(jl_var_825, lay_var_765))
        ratn2o_var_789 = 1D+20 * chi_n2o_var_790 / chi_mls(4, jp_var_736(jl_var_825, lay_var_765) + 1)
        IF (ratn2o_var_789 .GT. 1.5D0) THEN
          adjfac_var_787 = 0.5D0 + (ratn2o_var_789 - 0.5D0) ** 0.65D0
          adjcoln2o_var_788 = adjfac_var_787 * chi_mls(4, jp_var_736(jl_var_825, lay_var_765) + 1) * coldry_var_743(jl_var_825, lay_var_765) * 1D-20
        ELSE
          adjcoln2o_var_788 = coln2o_var_741(jl_var_825, lay_var_765)
        END IF
        speccomb_planck_var_772 = colh2o_var_740(jl_var_825, lay_var_765) + refrat_planck_a_var_773 * colch4_var_742(jl_var_825, lay_var_765)
        specparm_planck_var_786 = MIN(colh2o_var_740(jl_var_825, lay_var_765) / speccomb_planck_var_772, oneminus_var_739)
        specmult_planck_var_785 = 8.0D0 * specparm_planck_var_786
        jpl_var_768 = 1 + INT(specmult_planck_var_785)
        fpl_var_784 = ((specmult_planck_var_785) - AINT((specmult_planck_var_785)))
        ind0_var_758 = ((jp_var_736(jl_var_825, lay_var_765) - 1) * 5 + (jt_var_737(jl_var_825, lay_var_765) - 1)) * nspa_var_236(9) + js_var_764
        ind1_var_759 = (jp_var_736(jl_var_825, lay_var_765) * 5 + (jt1_var_738(jl_var_825, lay_var_765) - 1)) * nspa_var_236(9) + js1_var_766
        inds_var_760 = indself_var_749(jl_var_825, lay_var_765)
        indf_var_761 = indfor_var_753(jl_var_825, lay_var_765)
        indm_var_762 = indminor_var_757(jl_var_825, lay_var_765)
        IF (specparm_var_777 .LT. 0.125D0) THEN
          p_var_803 = fs_var_775 - 1.0D0
          p4_var_804 = p_var_803 ** 4
          fk0_var_805 = p4_var_804
          fk1_var_806 = 1.0D0 - p_var_803 - 2.0D0 * p4_var_804
          fk2_var_807 = p_var_803 + p4_var_804
          fac000_var_791 = fk0_var_805 * fac00_var_732(jl_var_825, lay_var_765)
          fac100_var_792 = fk1_var_806 * fac00_var_732(jl_var_825, lay_var_765)
          fac200_var_793 = fk2_var_807 * fac00_var_732(jl_var_825, lay_var_765)
          fac010_var_794 = fk0_var_805 * fac10_var_734(jl_var_825, lay_var_765)
          fac110_var_795 = fk1_var_806 * fac10_var_734(jl_var_825, lay_var_765)
          fac210_var_796 = fk2_var_807 * fac10_var_734(jl_var_825, lay_var_765)
        ELSE IF (specparm_var_777 .GT. 0.875D0) THEN
          p_var_803 = - fs_var_775
          p4_var_804 = p_var_803 ** 4
          fk0_var_805 = p4_var_804
          fk1_var_806 = 1.0D0 - p_var_803 - 2.0D0 * p4_var_804
          fk2_var_807 = p_var_803 + p4_var_804
          fac000_var_791 = fk0_var_805 * fac00_var_732(jl_var_825, lay_var_765)
          fac100_var_792 = fk1_var_806 * fac00_var_732(jl_var_825, lay_var_765)
          fac200_var_793 = fk2_var_807 * fac00_var_732(jl_var_825, lay_var_765)
          fac010_var_794 = fk0_var_805 * fac10_var_734(jl_var_825, lay_var_765)
          fac110_var_795 = fk1_var_806 * fac10_var_734(jl_var_825, lay_var_765)
          fac210_var_796 = fk2_var_807 * fac10_var_734(jl_var_825, lay_var_765)
        ELSE
          fac000_var_791 = (1.0D0 - fs_var_775) * fac00_var_732(jl_var_825, lay_var_765)
          fac010_var_794 = (1.0D0 - fs_var_775) * fac10_var_734(jl_var_825, lay_var_765)
          fac100_var_792 = fs_var_775 * fac00_var_732(jl_var_825, lay_var_765)
          fac110_var_795 = fs_var_775 * fac10_var_734(jl_var_825, lay_var_765)
          fac200_var_793 = 0.0D0
          fac210_var_796 = 0.0D0
        END IF
        IF (specparm1_var_780 .LT. 0.125D0) THEN
          p_var_803 = fs1_var_778 - 1.0D0
          p4_var_804 = p_var_803 ** 4
          fk0_var_805 = p4_var_804
          fk1_var_806 = 1.0D0 - p_var_803 - 2.0D0 * p4_var_804
          fk2_var_807 = p_var_803 + p4_var_804
          fac001_var_797 = fk0_var_805 * fac01_var_733(jl_var_825, lay_var_765)
          fac101_var_798 = fk1_var_806 * fac01_var_733(jl_var_825, lay_var_765)
          fac201_var_799 = fk2_var_807 * fac01_var_733(jl_var_825, lay_var_765)
          fac011_var_800 = fk0_var_805 * fac11_var_735(jl_var_825, lay_var_765)
          fac111_var_801 = fk1_var_806 * fac11_var_735(jl_var_825, lay_var_765)
          fac211_var_802 = fk2_var_807 * fac11_var_735(jl_var_825, lay_var_765)
        ELSE IF (specparm1_var_780 .GT. 0.875D0) THEN
          p_var_803 = - fs1_var_778
          p4_var_804 = p_var_803 ** 4
          fk0_var_805 = p4_var_804
          fk1_var_806 = 1.0D0 - p_var_803 - 2.0D0 * p4_var_804
          fk2_var_807 = p_var_803 + p4_var_804
          fac001_var_797 = fk0_var_805 * fac01_var_733(jl_var_825, lay_var_765)
          fac101_var_798 = fk1_var_806 * fac01_var_733(jl_var_825, lay_var_765)
          fac201_var_799 = fk2_var_807 * fac01_var_733(jl_var_825, lay_var_765)
          fac011_var_800 = fk0_var_805 * fac11_var_735(jl_var_825, lay_var_765)
          fac111_var_801 = fk1_var_806 * fac11_var_735(jl_var_825, lay_var_765)
          fac211_var_802 = fk2_var_807 * fac11_var_735(jl_var_825, lay_var_765)
        ELSE
          fac001_var_797 = (1.0D0 - fs1_var_778) * fac01_var_733(jl_var_825, lay_var_765)
          fac011_var_800 = (1.0D0 - fs1_var_778) * fac11_var_735(jl_var_825, lay_var_765)
          fac101_var_798 = fs1_var_778 * fac01_var_733(jl_var_825, lay_var_765)
          fac111_var_801 = fs1_var_778 * fac11_var_735(jl_var_825, lay_var_765)
          fac201_var_799 = 0.0D0
          fac211_var_802 = 0.0D0
        END IF
        IF (specparm_var_777 .LT. 0.125D0) THEN
          tau_major_var_813(1 : ng9) = speccomb_var_769 * (fac000_var_791 * absa_var_228(ind0_var_758, 1 : 12) + fac100_var_792 * absa_var_228(ind0_var_758 + 1, 1 : 12) + fac200_var_793 * absa_var_228(ind0_var_758 + 2, 1 : 12) + fac010_var_794 * absa_var_228(ind0_var_758 + 9, 1 : 12) + fac110_var_795 * absa_var_228(ind0_var_758 + 10, 1 : 12) + fac210_var_796 * absa_var_228(ind0_var_758 + 11, 1 : 12))
        ELSE IF (specparm_var_777 .GT. 0.875D0) THEN
          tau_major_var_813(1 : ng9) = speccomb_var_769 * (fac200_var_793 * absa_var_228(ind0_var_758 - 1, 1 : 12) + fac100_var_792 * absa_var_228(ind0_var_758, 1 : 12) + fac000_var_791 * absa_var_228(ind0_var_758 + 1, 1 : 12) + fac210_var_796 * absa_var_228(ind0_var_758 + 8, 1 : 12) + fac110_var_795 * absa_var_228(ind0_var_758 + 9, 1 : 12) + fac010_var_794 * absa_var_228(ind0_var_758 + 10, 1 : 12))
        ELSE
          tau_major_var_813(1 : ng9) = speccomb_var_769 * (fac000_var_791 * absa_var_228(ind0_var_758, 1 : 12) + fac100_var_792 * absa_var_228(ind0_var_758 + 1, 1 : 12) + fac010_var_794 * absa_var_228(ind0_var_758 + 9, 1 : 12) + fac110_var_795 * absa_var_228(ind0_var_758 + 10, 1 : 12))
        END IF
        IF (specparm1_var_780 .LT. 0.125D0) THEN
          tau_major1_var_814(1 : ng9) = speccomb1_var_770 * (fac001_var_797 * absa_var_228(ind1_var_759, 1 : 12) + fac101_var_798 * absa_var_228(ind1_var_759 + 1, 1 : 12) + fac201_var_799 * absa_var_228(ind1_var_759 + 2, 1 : 12) + fac011_var_800 * absa_var_228(ind1_var_759 + 9, 1 : 12) + fac111_var_801 * absa_var_228(ind1_var_759 + 10, 1 : 12) + fac211_var_802 * absa_var_228(ind1_var_759 + 11, 1 : 12))
        ELSE IF (specparm1_var_780 .GT. 0.875D0) THEN
          tau_major1_var_814(1 : ng9) = speccomb1_var_770 * (fac201_var_799 * absa_var_228(ind1_var_759 - 1, 1 : 12) + fac101_var_798 * absa_var_228(ind1_var_759, 1 : 12) + fac001_var_797 * absa_var_228(ind1_var_759 + 1, 1 : 12) + fac211_var_802 * absa_var_228(ind1_var_759 + 8, 1 : 12) + fac111_var_801 * absa_var_228(ind1_var_759 + 9, 1 : 12) + fac011_var_800 * absa_var_228(ind1_var_759 + 10, 1 : 12))
        ELSE
          tau_major1_var_814(1 : ng9) = speccomb1_var_770 * (fac001_var_797 * absa_var_228(ind1_var_759, 1 : 12) + fac101_var_798 * absa_var_228(ind1_var_759 + 1, 1 : 12) + fac011_var_800 * absa_var_228(ind1_var_759 + 9, 1 : 12) + fac111_var_801 * absa_var_228(ind1_var_759 + 10, 1 : 12))
        END IF
        DO ig_var_763 = 1, 12
          tauself_var_809 = selffac_var_747(jl_var_825, lay_var_765) * (selfref_var_232(inds_var_760, ig_var_763) + selffrac_var_748(jl_var_825, lay_var_765) * (selfref_var_232(inds_var_760 + 1, ig_var_763) - selfref_var_232(inds_var_760, ig_var_763)))
          taufor_var_808 = forfac_var_754(jl_var_825, lay_var_765) * (forref_var_233(indf_var_761, ig_var_763) + forfrac_var_755(jl_var_825, lay_var_765) * (forref_var_233(indf_var_761 + 1, ig_var_763) - forref_var_233(indf_var_761, ig_var_763)))
          n2om1_var_810 = ka_mn2o_var_230(jmn2o_var_767, indm_var_762, ig_var_763) + fmn2o_var_781 * (ka_mn2o_var_230(jmn2o_var_767 + 1, indm_var_762, ig_var_763) - ka_mn2o_var_230(jmn2o_var_767, indm_var_762, ig_var_763))
          n2om2_var_811 = ka_mn2o_var_230(jmn2o_var_767, indm_var_762 + 1, ig_var_763) + fmn2o_var_781 * (ka_mn2o_var_230(jmn2o_var_767 + 1, indm_var_762 + 1, ig_var_763) - ka_mn2o_var_230(jmn2o_var_767, indm_var_762 + 1, ig_var_763))
          absn2o_var_812 = n2om1_var_810 + minorfrac_var_756(jl_var_825, lay_var_765) * (n2om2_var_811 - n2om1_var_810)
          taug_var_730(jl_var_825, 96 + ig_var_763, lay_var_765) = tau_major_var_813(ig_var_763) + tau_major1_var_814(ig_var_763) + tauself_var_809 + taufor_var_808 + adjcoln2o_var_788 * absn2o_var_812
          fracs_var_750(jl_var_825, 96 + ig_var_763, lay_var_765) = fracrefa_var_226(ig_var_763, jpl_var_768) + fpl_var_784 * (fracrefa_var_226(ig_var_763, jpl_var_768 + 1) - fracrefa_var_226(ig_var_763, jpl_var_768))
        END DO
      END DO
      ixc0_var_822 = kfdia_var_728 - kidia_var_727 + 1 - ixc0_var_822
      DO ixp_var_823 = 1, ixc0_var_822
        jl_var_825 = ixhigh_var_819(ixp_var_823, lay_var_765)
        chi_n2o_var_790 = coln2o_var_741(jl_var_825, lay_var_765) / (coldry_var_743(jl_var_825, lay_var_765))
        ratn2o_var_789 = 1D+20 * chi_n2o_var_790 / chi_mls(4, jp_var_736(jl_var_825, lay_var_765) + 1)
        IF (ratn2o_var_789 .GT. 1.5D0) THEN
          adjfac_var_787 = 0.5D0 + (ratn2o_var_789 - 0.5D0) ** 0.65D0
          adjcoln2o_var_788 = adjfac_var_787 * chi_mls(4, jp_var_736(jl_var_825, lay_var_765) + 1) * coldry_var_743(jl_var_825, lay_var_765) * 1D-20
        ELSE
          adjcoln2o_var_788 = coln2o_var_741(jl_var_825, lay_var_765)
        END IF
        ind0_var_758 = ((jp_var_736(jl_var_825, lay_var_765) - 13) * 5 + (jt_var_737(jl_var_825, lay_var_765) - 1)) * nspb_var_237(9) + 1
        ind1_var_759 = ((jp_var_736(jl_var_825, lay_var_765) - 12) * 5 + (jt1_var_738(jl_var_825, lay_var_765) - 1)) * nspb_var_237(9) + 1
        indm_var_762 = indminor_var_757(jl_var_825, lay_var_765)
        DO ig_var_763 = 1, 12
          absn2o_var_812 = kb_mn2o_var_231(indm_var_762, ig_var_763) + minorfrac_var_756(jl_var_825, lay_var_765) * (kb_mn2o_var_231(indm_var_762 + 1, ig_var_763) - kb_mn2o_var_231(indm_var_762, ig_var_763))
          taug_var_730(jl_var_825, 96 + ig_var_763, lay_var_765) = colch4_var_742(jl_var_825, lay_var_765) * (fac00_var_732(jl_var_825, lay_var_765) * absb_var_229(ind0_var_758, ig_var_763) + fac10_var_734(jl_var_825, lay_var_765) * absb_var_229(ind0_var_758 + 1, ig_var_763) + fac01_var_733(jl_var_825, lay_var_765) * absb_var_229(ind1_var_759, ig_var_763) + fac11_var_735(jl_var_825, lay_var_765) * absb_var_229(ind1_var_759 + 1, ig_var_763)) + adjcoln2o_var_788 * absn2o_var_812
          fracs_var_750(jl_var_825, 96 + ig_var_763, lay_var_765) = fracrefb_var_227(ig_var_763)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol9
SUBROUTINE rrtm_taumol11(kidia_var_826, kfdia_var_827, klev_var_828, taug_var_829, p_tauaerl_var_830, fac00_var_831, fac01_var_832, fac10_var_833, fac11_var_834, forfac_var_845, forfrac_var_846, indfor_var_844, jp_var_835, jt_var_836, jt1_var_837, colh2o_var_838, colo2, laytrop_var_839, selffac_var_840, selffrac_var_841, indself_var_842, fracs_var_843, minorfrac_var_847, indminor_var_848, scaleminor_var_849)
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrta11, ONLY: absa_var_142, absb_var_143, forref_var_145, fracrefa_var_140, fracrefb_var_141, ka_mo2, kb_mo2, selfref_var_144
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_826
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_827
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_828
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_829(kidia_var_826 : kfdia_var_827, 140, klev_var_828)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_830(kidia_var_826 : kfdia_var_827, klev_var_828, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_831(kidia_var_826 : kfdia_var_827, klev_var_828)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_832(kidia_var_826 : kfdia_var_827, klev_var_828)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_833(kidia_var_826 : kfdia_var_827, klev_var_828)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_834(kidia_var_826 : kfdia_var_827, klev_var_828)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_835(kidia_var_826 : kfdia_var_827, klev_var_828)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_836(kidia_var_826 : kfdia_var_827, klev_var_828)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_837(kidia_var_826 : kfdia_var_827, klev_var_828)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_838(kidia_var_826 : kfdia_var_827, klev_var_828)
  REAL(KIND = 8), INTENT(IN) :: colo2(kidia_var_826 : kfdia_var_827, klev_var_828)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_839(kidia_var_826 : kfdia_var_827)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_840(kidia_var_826 : kfdia_var_827, klev_var_828)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_841(kidia_var_826 : kfdia_var_827, klev_var_828)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_842(kidia_var_826 : kfdia_var_827, klev_var_828)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_843(kidia_var_826 : kfdia_var_827, 140, klev_var_828)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_844(kidia_var_826 : kfdia_var_827, klev_var_828)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_845(kidia_var_826 : kfdia_var_827, klev_var_828)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_846(kidia_var_826 : kfdia_var_827, klev_var_828)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_847(kidia_var_826 : kfdia_var_827, klev_var_828)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_848(kidia_var_826 : kfdia_var_827, klev_var_828)
  REAL(KIND = 8), INTENT(IN) :: scaleminor_var_849(kidia_var_826 : kfdia_var_827, klev_var_828)
  INTEGER(KIND = 4) :: ind0_var_850, ind1_var_851
  INTEGER(KIND = 4) :: inds_var_852, indf_var_853, indm_var_854
  INTEGER(KIND = 4) :: ig_var_855, lay_var_856
  REAL(KIND = 8) :: taufor_var_857, tauself_var_858, scaleo2, tauo2
  INTEGER(KIND = 4) :: laytrop_min_var_859, laytrop_max_var_860
  INTEGER(KIND = 4) :: ixc_var_861(klev_var_828), ixlow_var_862(kfdia_var_827, klev_var_828), ixhigh_var_863(kfdia_var_827, klev_var_828)
  INTEGER(KIND = 4) :: ich_var_864, icl_var_865, ixc0_var_866, ixp_var_867, jc_var_868, jl_var_869
  laytrop_min_var_859 = MINVAL(laytrop_var_839)
  laytrop_max_var_860 = MAXVAL(laytrop_var_839)
  ixlow_var_862 = 0
  ixhigh_var_863 = 0
  ixc_var_861 = 0
  DO lay_var_856 = laytrop_min_var_859 + 1, laytrop_max_var_860
    icl_var_865 = 0
    ich_var_864 = 0
    DO jc_var_868 = kidia_var_826, kfdia_var_827
      IF (lay_var_856 <= laytrop_var_839(jc_var_868)) THEN
        icl_var_865 = icl_var_865 + 1
        ixlow_var_862(icl_var_865, lay_var_856) = jc_var_868
      ELSE
        ich_var_864 = ich_var_864 + 1
        ixhigh_var_863(ich_var_864, lay_var_856) = jc_var_868
      END IF
    END DO
    ixc_var_861(lay_var_856) = icl_var_865
  END DO
  DO lay_var_856 = 1, laytrop_min_var_859
    DO jl_var_869 = kidia_var_826, kfdia_var_827
      ind0_var_850 = ((jp_var_835(jl_var_869, lay_var_856) - 1) * 5 + (jt_var_836(jl_var_869, lay_var_856) - 1)) * nspa_var_236(11) + 1
      ind1_var_851 = (jp_var_835(jl_var_869, lay_var_856) * 5 + (jt1_var_837(jl_var_869, lay_var_856) - 1)) * nspa_var_236(11) + 1
      inds_var_852 = indself_var_842(jl_var_869, lay_var_856)
      indf_var_853 = indfor_var_844(jl_var_869, lay_var_856)
      indm_var_854 = indminor_var_848(jl_var_869, lay_var_856)
      scaleo2 = colo2(jl_var_869, lay_var_856) * scaleminor_var_849(jl_var_869, lay_var_856)
      DO ig_var_855 = 1, 8
        tauself_var_858 = selffac_var_840(jl_var_869, lay_var_856) * (selfref_var_144(inds_var_852, ig_var_855) + selffrac_var_841(jl_var_869, lay_var_856) * (selfref_var_144(inds_var_852 + 1, ig_var_855) - selfref_var_144(inds_var_852, ig_var_855)))
        taufor_var_857 = forfac_var_845(jl_var_869, lay_var_856) * (forref_var_145(indf_var_853, ig_var_855) + forfrac_var_846(jl_var_869, lay_var_856) * (forref_var_145(indf_var_853 + 1, ig_var_855) - forref_var_145(indf_var_853, ig_var_855)))
        tauo2 = scaleo2 * (ka_mo2(indm_var_854, ig_var_855) + minorfrac_var_847(jl_var_869, lay_var_856) * (ka_mo2(indm_var_854 + 1, ig_var_855) - ka_mo2(indm_var_854, ig_var_855)))
        taug_var_829(jl_var_869, 114 + ig_var_855, lay_var_856) = colh2o_var_838(jl_var_869, lay_var_856) * (fac00_var_831(jl_var_869, lay_var_856) * absa_var_142(ind0_var_850, ig_var_855) + fac10_var_833(jl_var_869, lay_var_856) * absa_var_142(ind0_var_850 + 1, ig_var_855) + fac01_var_832(jl_var_869, lay_var_856) * absa_var_142(ind1_var_851, ig_var_855) + fac11_var_834(jl_var_869, lay_var_856) * absa_var_142(ind1_var_851 + 1, ig_var_855)) + tauself_var_858 + taufor_var_857 + tauo2
        fracs_var_843(jl_var_869, 114 + ig_var_855, lay_var_856) = fracrefa_var_140(ig_var_855)
      END DO
    END DO
  END DO
  DO lay_var_856 = laytrop_max_var_860 + 1, klev_var_828
    DO jl_var_869 = kidia_var_826, kfdia_var_827
      ind0_var_850 = ((jp_var_835(jl_var_869, lay_var_856) - 13) * 5 + (jt_var_836(jl_var_869, lay_var_856) - 1)) * nspb_var_237(11) + 1
      ind1_var_851 = ((jp_var_835(jl_var_869, lay_var_856) - 12) * 5 + (jt1_var_837(jl_var_869, lay_var_856) - 1)) * nspb_var_237(11) + 1
      indf_var_853 = indfor_var_844(jl_var_869, lay_var_856)
      indm_var_854 = indminor_var_848(jl_var_869, lay_var_856)
      scaleo2 = colo2(jl_var_869, lay_var_856) * scaleminor_var_849(jl_var_869, lay_var_856)
      DO ig_var_855 = 1, 8
        taufor_var_857 = forfac_var_845(jl_var_869, lay_var_856) * (forref_var_145(indf_var_853, ig_var_855) + forfrac_var_846(jl_var_869, lay_var_856) * (forref_var_145(indf_var_853 + 1, ig_var_855) - forref_var_145(indf_var_853, ig_var_855)))
        tauo2 = scaleo2 * (kb_mo2(indm_var_854, ig_var_855) + minorfrac_var_847(jl_var_869, lay_var_856) * (kb_mo2(indm_var_854 + 1, ig_var_855) - kb_mo2(indm_var_854, ig_var_855)))
        taug_var_829(jl_var_869, 114 + ig_var_855, lay_var_856) = colh2o_var_838(jl_var_869, lay_var_856) * (fac00_var_831(jl_var_869, lay_var_856) * absb_var_143(ind0_var_850, ig_var_855) + fac10_var_833(jl_var_869, lay_var_856) * absb_var_143(ind0_var_850 + 1, ig_var_855) + fac01_var_832(jl_var_869, lay_var_856) * absb_var_143(ind1_var_851, ig_var_855) + fac11_var_834(jl_var_869, lay_var_856) * absb_var_143(ind1_var_851 + 1, ig_var_855)) + taufor_var_857 + tauo2
        fracs_var_843(jl_var_869, 114 + ig_var_855, lay_var_856) = fracrefb_var_141(ig_var_855)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_860 /= laytrop_min_var_859) THEN
    DO lay_var_856 = laytrop_min_var_859 + 1, laytrop_max_var_860
      ixc0_var_866 = ixc_var_861(lay_var_856)
      DO ixp_var_867 = 1, ixc0_var_866
        jl_var_869 = ixlow_var_862(ixp_var_867, lay_var_856)
        ind0_var_850 = ((jp_var_835(jl_var_869, lay_var_856) - 1) * 5 + (jt_var_836(jl_var_869, lay_var_856) - 1)) * nspa_var_236(11) + 1
        ind1_var_851 = (jp_var_835(jl_var_869, lay_var_856) * 5 + (jt1_var_837(jl_var_869, lay_var_856) - 1)) * nspa_var_236(11) + 1
        inds_var_852 = indself_var_842(jl_var_869, lay_var_856)
        indf_var_853 = indfor_var_844(jl_var_869, lay_var_856)
        indm_var_854 = indminor_var_848(jl_var_869, lay_var_856)
        scaleo2 = colo2(jl_var_869, lay_var_856) * scaleminor_var_849(jl_var_869, lay_var_856)
        DO ig_var_855 = 1, 8
          tauself_var_858 = selffac_var_840(jl_var_869, lay_var_856) * (selfref_var_144(inds_var_852, ig_var_855) + selffrac_var_841(jl_var_869, lay_var_856) * (selfref_var_144(inds_var_852 + 1, ig_var_855) - selfref_var_144(inds_var_852, ig_var_855)))
          taufor_var_857 = forfac_var_845(jl_var_869, lay_var_856) * (forref_var_145(indf_var_853, ig_var_855) + forfrac_var_846(jl_var_869, lay_var_856) * (forref_var_145(indf_var_853 + 1, ig_var_855) - forref_var_145(indf_var_853, ig_var_855)))
          tauo2 = scaleo2 * (ka_mo2(indm_var_854, ig_var_855) + minorfrac_var_847(jl_var_869, lay_var_856) * (ka_mo2(indm_var_854 + 1, ig_var_855) - ka_mo2(indm_var_854, ig_var_855)))
          taug_var_829(jl_var_869, 114 + ig_var_855, lay_var_856) = colh2o_var_838(jl_var_869, lay_var_856) * (fac00_var_831(jl_var_869, lay_var_856) * absa_var_142(ind0_var_850, ig_var_855) + fac10_var_833(jl_var_869, lay_var_856) * absa_var_142(ind0_var_850 + 1, ig_var_855) + fac01_var_832(jl_var_869, lay_var_856) * absa_var_142(ind1_var_851, ig_var_855) + fac11_var_834(jl_var_869, lay_var_856) * absa_var_142(ind1_var_851 + 1, ig_var_855)) + tauself_var_858 + taufor_var_857 + tauo2
          fracs_var_843(jl_var_869, 114 + ig_var_855, lay_var_856) = fracrefa_var_140(ig_var_855)
        END DO
      END DO
      ixc0_var_866 = kfdia_var_827 - kidia_var_826 + 1 - ixc0_var_866
      DO ixp_var_867 = 1, ixc0_var_866
        jl_var_869 = ixhigh_var_863(ixp_var_867, lay_var_856)
        ind0_var_850 = ((jp_var_835(jl_var_869, lay_var_856) - 13) * 5 + (jt_var_836(jl_var_869, lay_var_856) - 1)) * nspb_var_237(11) + 1
        ind1_var_851 = ((jp_var_835(jl_var_869, lay_var_856) - 12) * 5 + (jt1_var_837(jl_var_869, lay_var_856) - 1)) * nspb_var_237(11) + 1
        indf_var_853 = indfor_var_844(jl_var_869, lay_var_856)
        indm_var_854 = indminor_var_848(jl_var_869, lay_var_856)
        scaleo2 = colo2(jl_var_869, lay_var_856) * scaleminor_var_849(jl_var_869, lay_var_856)
        DO ig_var_855 = 1, 8
          taufor_var_857 = forfac_var_845(jl_var_869, lay_var_856) * (forref_var_145(indf_var_853, ig_var_855) + forfrac_var_846(jl_var_869, lay_var_856) * (forref_var_145(indf_var_853 + 1, ig_var_855) - forref_var_145(indf_var_853, ig_var_855)))
          tauo2 = scaleo2 * (kb_mo2(indm_var_854, ig_var_855) + minorfrac_var_847(jl_var_869, lay_var_856) * (kb_mo2(indm_var_854 + 1, ig_var_855) - kb_mo2(indm_var_854, ig_var_855)))
          taug_var_829(jl_var_869, 114 + ig_var_855, lay_var_856) = colh2o_var_838(jl_var_869, lay_var_856) * (fac00_var_831(jl_var_869, lay_var_856) * absb_var_143(ind0_var_850, ig_var_855) + fac10_var_833(jl_var_869, lay_var_856) * absb_var_143(ind0_var_850 + 1, ig_var_855) + fac01_var_832(jl_var_869, lay_var_856) * absb_var_143(ind1_var_851, ig_var_855) + fac11_var_834(jl_var_869, lay_var_856) * absb_var_143(ind1_var_851 + 1, ig_var_855)) + taufor_var_857 + tauo2
          fracs_var_843(jl_var_869, 114 + ig_var_855, lay_var_856) = fracrefb_var_141(ig_var_855)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol11
SUBROUTINE rrtm_taumol10(kidia_var_870, kfdia_var_871, klev_var_872, taug_var_873, p_tauaerl_var_874, fac00_var_875, fac01_var_876, fac10_var_877, fac11_var_878, forfac_var_890, forfrac_var_889, indfor_var_888, jp_var_879, jt_var_880, jt1_var_881, colh2o_var_882, laytrop_var_883, selffac_var_885, selffrac_var_886, indself_var_887, fracs_var_884)
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrta10, ONLY: absa_var_136, absb_var_137, forref_var_139, fracrefa_var_134, fracrefb_var_135, selfref_var_138
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_870
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_871
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_872
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_873(kidia_var_870 : kfdia_var_871, 140, klev_var_872)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_874(kidia_var_870 : kfdia_var_871, klev_var_872, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_875(kidia_var_870 : kfdia_var_871, klev_var_872)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_876(kidia_var_870 : kfdia_var_871, klev_var_872)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_877(kidia_var_870 : kfdia_var_871, klev_var_872)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_878(kidia_var_870 : kfdia_var_871, klev_var_872)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_879(kidia_var_870 : kfdia_var_871, klev_var_872)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_880(kidia_var_870 : kfdia_var_871, klev_var_872)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_881(kidia_var_870 : kfdia_var_871, klev_var_872)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_882(kidia_var_870 : kfdia_var_871, klev_var_872)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_883(kidia_var_870 : kfdia_var_871)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_884(kidia_var_870 : kfdia_var_871, 140, klev_var_872)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_885(kidia_var_870 : kfdia_var_871, klev_var_872)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_886(kidia_var_870 : kfdia_var_871, klev_var_872)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_887(kidia_var_870 : kfdia_var_871, klev_var_872)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_888(kidia_var_870 : kfdia_var_871, klev_var_872)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_889(kidia_var_870 : kfdia_var_871, klev_var_872)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_890(kidia_var_870 : kfdia_var_871, klev_var_872)
  INTEGER(KIND = 4) :: ind0_var_891, ind1_var_892
  INTEGER(KIND = 4) :: inds_var_893, indf_var_894
  INTEGER(KIND = 4) :: ig_var_895, lay_var_896
  REAL(KIND = 8) :: taufor_var_897, tauself_var_898
  INTEGER(KIND = 4) :: laytrop_min_var_899, laytrop_max_var_900
  INTEGER(KIND = 4) :: ixc_var_901(klev_var_872), ixlow_var_902(kfdia_var_871, klev_var_872), ixhigh_var_903(kfdia_var_871, klev_var_872)
  INTEGER(KIND = 4) :: ich_var_904, icl_var_905, ixc0_var_906, ixp_var_907, jc_var_908, jl_var_909
  laytrop_min_var_899 = MINVAL(laytrop_var_883)
  laytrop_max_var_900 = MAXVAL(laytrop_var_883)
  ixlow_var_902 = 0
  ixhigh_var_903 = 0
  ixc_var_901 = 0
  DO lay_var_896 = laytrop_min_var_899 + 1, laytrop_max_var_900
    icl_var_905 = 0
    ich_var_904 = 0
    DO jc_var_908 = kidia_var_870, kfdia_var_871
      IF (lay_var_896 <= laytrop_var_883(jc_var_908)) THEN
        icl_var_905 = icl_var_905 + 1
        ixlow_var_902(icl_var_905, lay_var_896) = jc_var_908
      ELSE
        ich_var_904 = ich_var_904 + 1
        ixhigh_var_903(ich_var_904, lay_var_896) = jc_var_908
      END IF
    END DO
    ixc_var_901(lay_var_896) = icl_var_905
  END DO
  DO lay_var_896 = 1, laytrop_min_var_899
    DO jl_var_909 = kidia_var_870, kfdia_var_871
      ind0_var_891 = ((jp_var_879(jl_var_909, lay_var_896) - 1) * 5 + (jt_var_880(jl_var_909, lay_var_896) - 1)) * nspa_var_236(10) + 1
      ind1_var_892 = (jp_var_879(jl_var_909, lay_var_896) * 5 + (jt1_var_881(jl_var_909, lay_var_896) - 1)) * nspa_var_236(10) + 1
      inds_var_893 = indself_var_887(jl_var_909, lay_var_896)
      indf_var_894 = indfor_var_888(jl_var_909, lay_var_896)
      DO ig_var_895 = 1, 6
        tauself_var_898 = selffac_var_885(jl_var_909, lay_var_896) * (selfref_var_138(inds_var_893, ig_var_895) + selffrac_var_886(jl_var_909, lay_var_896) * (selfref_var_138(inds_var_893 + 1, ig_var_895) - selfref_var_138(inds_var_893, ig_var_895)))
        taufor_var_897 = forfac_var_890(jl_var_909, lay_var_896) * (forref_var_139(indf_var_894, ig_var_895) + forfrac_var_889(jl_var_909, lay_var_896) * (forref_var_139(indf_var_894 + 1, ig_var_895) - forref_var_139(indf_var_894, ig_var_895)))
        taug_var_873(jl_var_909, 108 + ig_var_895, lay_var_896) = colh2o_var_882(jl_var_909, lay_var_896) * (fac00_var_875(jl_var_909, lay_var_896) * absa_var_136(ind0_var_891, ig_var_895) + fac10_var_877(jl_var_909, lay_var_896) * absa_var_136(ind0_var_891 + 1, ig_var_895) + fac01_var_876(jl_var_909, lay_var_896) * absa_var_136(ind1_var_892, ig_var_895) + fac11_var_878(jl_var_909, lay_var_896) * absa_var_136(ind1_var_892 + 1, ig_var_895)) + tauself_var_898 + taufor_var_897
        fracs_var_884(jl_var_909, 108 + ig_var_895, lay_var_896) = fracrefa_var_134(ig_var_895)
      END DO
    END DO
  END DO
  DO lay_var_896 = laytrop_max_var_900 + 1, klev_var_872
    DO jl_var_909 = kidia_var_870, kfdia_var_871
      ind0_var_891 = ((jp_var_879(jl_var_909, lay_var_896) - 13) * 5 + (jt_var_880(jl_var_909, lay_var_896) - 1)) * nspb_var_237(10) + 1
      ind1_var_892 = ((jp_var_879(jl_var_909, lay_var_896) - 12) * 5 + (jt1_var_881(jl_var_909, lay_var_896) - 1)) * nspb_var_237(10) + 1
      indf_var_894 = indfor_var_888(jl_var_909, lay_var_896)
      DO ig_var_895 = 1, 6
        taufor_var_897 = forfac_var_890(jl_var_909, lay_var_896) * (forref_var_139(indf_var_894, ig_var_895) + forfrac_var_889(jl_var_909, lay_var_896) * (forref_var_139(indf_var_894 + 1, ig_var_895) - forref_var_139(indf_var_894, ig_var_895)))
        taug_var_873(jl_var_909, 108 + ig_var_895, lay_var_896) = colh2o_var_882(jl_var_909, lay_var_896) * (fac00_var_875(jl_var_909, lay_var_896) * absb_var_137(ind0_var_891, ig_var_895) + fac10_var_877(jl_var_909, lay_var_896) * absb_var_137(ind0_var_891 + 1, ig_var_895) + fac01_var_876(jl_var_909, lay_var_896) * absb_var_137(ind1_var_892, ig_var_895) + fac11_var_878(jl_var_909, lay_var_896) * absb_var_137(ind1_var_892 + 1, ig_var_895)) + taufor_var_897
        fracs_var_884(jl_var_909, 108 + ig_var_895, lay_var_896) = fracrefb_var_135(ig_var_895)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_900 /= laytrop_min_var_899) THEN
    DO lay_var_896 = laytrop_min_var_899 + 1, laytrop_max_var_900
      ixc0_var_906 = ixc_var_901(lay_var_896)
      DO ixp_var_907 = 1, ixc0_var_906
        jl_var_909 = ixlow_var_902(ixp_var_907, lay_var_896)
        ind0_var_891 = ((jp_var_879(jl_var_909, lay_var_896) - 1) * 5 + (jt_var_880(jl_var_909, lay_var_896) - 1)) * nspa_var_236(10) + 1
        ind1_var_892 = (jp_var_879(jl_var_909, lay_var_896) * 5 + (jt1_var_881(jl_var_909, lay_var_896) - 1)) * nspa_var_236(10) + 1
        inds_var_893 = indself_var_887(jl_var_909, lay_var_896)
        indf_var_894 = indfor_var_888(jl_var_909, lay_var_896)
        DO ig_var_895 = 1, 6
          tauself_var_898 = selffac_var_885(jl_var_909, lay_var_896) * (selfref_var_138(inds_var_893, ig_var_895) + selffrac_var_886(jl_var_909, lay_var_896) * (selfref_var_138(inds_var_893 + 1, ig_var_895) - selfref_var_138(inds_var_893, ig_var_895)))
          taufor_var_897 = forfac_var_890(jl_var_909, lay_var_896) * (forref_var_139(indf_var_894, ig_var_895) + forfrac_var_889(jl_var_909, lay_var_896) * (forref_var_139(indf_var_894 + 1, ig_var_895) - forref_var_139(indf_var_894, ig_var_895)))
          taug_var_873(jl_var_909, 108 + ig_var_895, lay_var_896) = colh2o_var_882(jl_var_909, lay_var_896) * (fac00_var_875(jl_var_909, lay_var_896) * absa_var_136(ind0_var_891, ig_var_895) + fac10_var_877(jl_var_909, lay_var_896) * absa_var_136(ind0_var_891 + 1, ig_var_895) + fac01_var_876(jl_var_909, lay_var_896) * absa_var_136(ind1_var_892, ig_var_895) + fac11_var_878(jl_var_909, lay_var_896) * absa_var_136(ind1_var_892 + 1, ig_var_895)) + tauself_var_898 + taufor_var_897
          fracs_var_884(jl_var_909, 108 + ig_var_895, lay_var_896) = fracrefa_var_134(ig_var_895)
        END DO
      END DO
      ixc0_var_906 = kfdia_var_871 - kidia_var_870 + 1 - ixc0_var_906
      DO ixp_var_907 = 1, ixc0_var_906
        jl_var_909 = ixhigh_var_903(ixp_var_907, lay_var_896)
        ind0_var_891 = ((jp_var_879(jl_var_909, lay_var_896) - 13) * 5 + (jt_var_880(jl_var_909, lay_var_896) - 1)) * nspb_var_237(10) + 1
        ind1_var_892 = ((jp_var_879(jl_var_909, lay_var_896) - 12) * 5 + (jt1_var_881(jl_var_909, lay_var_896) - 1)) * nspb_var_237(10) + 1
        indf_var_894 = indfor_var_888(jl_var_909, lay_var_896)
        DO ig_var_895 = 1, 6
          taufor_var_897 = forfac_var_890(jl_var_909, lay_var_896) * (forref_var_139(indf_var_894, ig_var_895) + forfrac_var_889(jl_var_909, lay_var_896) * (forref_var_139(indf_var_894 + 1, ig_var_895) - forref_var_139(indf_var_894, ig_var_895)))
          taug_var_873(jl_var_909, 108 + ig_var_895, lay_var_896) = colh2o_var_882(jl_var_909, lay_var_896) * (fac00_var_875(jl_var_909, lay_var_896) * absb_var_137(ind0_var_891, ig_var_895) + fac10_var_877(jl_var_909, lay_var_896) * absb_var_137(ind0_var_891 + 1, ig_var_895) + fac01_var_876(jl_var_909, lay_var_896) * absb_var_137(ind1_var_892, ig_var_895) + fac11_var_878(jl_var_909, lay_var_896) * absb_var_137(ind1_var_892 + 1, ig_var_895)) + taufor_var_897
          fracs_var_884(jl_var_909, 108 + ig_var_895, lay_var_896) = fracrefb_var_135(ig_var_895)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol10
SUBROUTINE rrtm_taumol8(kidia_var_910, kfdia_var_911, klev_var_912, taug_var_913, wx_var_914, p_tauaerl_var_915, fac00_var_916, fac01_var_917, fac10_var_918, fac11_var_919, forfac_var_935, forfrac_var_934, indfor_var_933, jp_var_920, jt_var_921, jt1_var_922, colh2o_var_923, colo3_var_924, coln2o_var_925, colco2_var_926, coldry_var_927, laytrop_var_928, selffac_var_929, selffrac_var_930, indself_var_931, fracs_var_932, minorfrac_var_936, indminor_var_937)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrta8, ONLY: absa_var_217, absb_var_218, cfc12_var_216, cfc22adj, forref_var_225, fracrefa_var_214, fracrefb_var_215, ka_mco2_var_219, ka_mn2o_var_220, ka_mo3_var_221, kb_mco2_var_222, kb_mn2o_var_223, selfref_var_224
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_910
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_911
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_912
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_913(kidia_var_910 : kfdia_var_911, 140, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: wx_var_914(kidia_var_910 : kfdia_var_911, 4, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_915(kidia_var_910 : kfdia_var_911, klev_var_912, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_916(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_917(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_918(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_919(kidia_var_910 : kfdia_var_911, klev_var_912)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_920(kidia_var_910 : kfdia_var_911, klev_var_912)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_921(kidia_var_910 : kfdia_var_911, klev_var_912)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_922(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_923(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_924(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_925(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_926(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_927(kidia_var_910 : kfdia_var_911, klev_var_912)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_928(kidia_var_910 : kfdia_var_911)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_929(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_930(kidia_var_910 : kfdia_var_911, klev_var_912)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_931(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_932(kidia_var_910 : kfdia_var_911, 140, klev_var_912)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_933(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_934(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_935(kidia_var_910 : kfdia_var_911, klev_var_912)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_936(kidia_var_910 : kfdia_var_911, klev_var_912)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_937(kidia_var_910 : kfdia_var_911, klev_var_912)
  INTEGER(KIND = 4) :: ind0_var_938, ind1_var_939, inds_var_940, indf_var_941, indm_var_942
  INTEGER(KIND = 4) :: ig_var_943, lay_var_944
  REAL(KIND = 8) :: chi_co2_var_945, ratco2_var_946, adjfac_var_947, adjcolco2_var_948
  REAL(KIND = 8) :: taufor_var_949, tauself_var_950, abso3_var_951, absco2_var_952, absn2o_var_953
  INTEGER(KIND = 4) :: laytrop_min_var_954, laytrop_max_var_955
  INTEGER(KIND = 4) :: ixc_var_956(klev_var_912), ixlow_var_957(kfdia_var_911, klev_var_912), ixhigh_var_958(kfdia_var_911, klev_var_912)
  INTEGER(KIND = 4) :: ich_var_959, icl_var_960, ixc0_var_961, ixp_var_962, jc_var_963, jl_var_964
  laytrop_min_var_954 = MINVAL(laytrop_var_928)
  laytrop_max_var_955 = MAXVAL(laytrop_var_928)
  ixlow_var_957 = 0
  ixhigh_var_958 = 0
  ixc_var_956 = 0
  DO lay_var_944 = laytrop_min_var_954 + 1, laytrop_max_var_955
    icl_var_960 = 0
    ich_var_959 = 0
    DO jc_var_963 = kidia_var_910, kfdia_var_911
      IF (lay_var_944 <= laytrop_var_928(jc_var_963)) THEN
        icl_var_960 = icl_var_960 + 1
        ixlow_var_957(icl_var_960, lay_var_944) = jc_var_963
      ELSE
        ich_var_959 = ich_var_959 + 1
        ixhigh_var_958(ich_var_959, lay_var_944) = jc_var_963
      END IF
    END DO
    ixc_var_956(lay_var_944) = icl_var_960
  END DO
  DO lay_var_944 = 1, laytrop_min_var_954
    DO jl_var_964 = kidia_var_910, kfdia_var_911
      chi_co2_var_945 = colco2_var_926(jl_var_964, lay_var_944) / (coldry_var_927(jl_var_964, lay_var_944))
      ratco2_var_946 = 1D+20 * chi_co2_var_945 / chi_mls(2, jp_var_920(jl_var_964, lay_var_944) + 1)
      IF (ratco2_var_946 .GT. 3.0D0) THEN
        adjfac_var_947 = 2.0D0 + (ratco2_var_946 - 2.0D0) ** 0.65D0
        adjcolco2_var_948 = adjfac_var_947 * chi_mls(2, jp_var_920(jl_var_964, lay_var_944) + 1) * coldry_var_927(jl_var_964, lay_var_944) * 1D-20
      ELSE
        adjcolco2_var_948 = colco2_var_926(jl_var_964, lay_var_944)
      END IF
      ind0_var_938 = ((jp_var_920(jl_var_964, lay_var_944) - 1) * 5 + (jt_var_921(jl_var_964, lay_var_944) - 1)) * nspa_var_236(8) + 1
      ind1_var_939 = (jp_var_920(jl_var_964, lay_var_944) * 5 + (jt1_var_922(jl_var_964, lay_var_944) - 1)) * nspa_var_236(8) + 1
      inds_var_940 = indself_var_931(jl_var_964, lay_var_944)
      indf_var_941 = indfor_var_933(jl_var_964, lay_var_944)
      indm_var_942 = indminor_var_937(jl_var_964, lay_var_944)
      DO ig_var_943 = 1, 8
        tauself_var_950 = selffac_var_929(jl_var_964, lay_var_944) * (selfref_var_224(inds_var_940, ig_var_943) + selffrac_var_930(jl_var_964, lay_var_944) * (selfref_var_224(inds_var_940 + 1, ig_var_943) - selfref_var_224(inds_var_940, ig_var_943)))
        taufor_var_949 = forfac_var_935(jl_var_964, lay_var_944) * (forref_var_225(indf_var_941, ig_var_943) + forfrac_var_934(jl_var_964, lay_var_944) * (forref_var_225(indf_var_941 + 1, ig_var_943) - forref_var_225(indf_var_941, ig_var_943)))
        absco2_var_952 = (ka_mco2_var_219(indm_var_942, ig_var_943) + minorfrac_var_936(jl_var_964, lay_var_944) * (ka_mco2_var_219(indm_var_942 + 1, ig_var_943) - ka_mco2_var_219(indm_var_942, ig_var_943)))
        abso3_var_951 = (ka_mo3_var_221(indm_var_942, ig_var_943) + minorfrac_var_936(jl_var_964, lay_var_944) * (ka_mo3_var_221(indm_var_942 + 1, ig_var_943) - ka_mo3_var_221(indm_var_942, ig_var_943)))
        absn2o_var_953 = (ka_mn2o_var_220(indm_var_942, ig_var_943) + minorfrac_var_936(jl_var_964, lay_var_944) * (ka_mn2o_var_220(indm_var_942 + 1, ig_var_943) - ka_mn2o_var_220(indm_var_942, ig_var_943)))
        taug_var_913(jl_var_964, 88 + ig_var_943, lay_var_944) = colh2o_var_923(jl_var_964, lay_var_944) * (fac00_var_916(jl_var_964, lay_var_944) * absa_var_217(ind0_var_938, ig_var_943) + fac10_var_918(jl_var_964, lay_var_944) * absa_var_217(ind0_var_938 + 1, ig_var_943) + fac01_var_917(jl_var_964, lay_var_944) * absa_var_217(ind1_var_939, ig_var_943) + fac11_var_919(jl_var_964, lay_var_944) * absa_var_217(ind1_var_939 + 1, ig_var_943)) + tauself_var_950 + taufor_var_949 + adjcolco2_var_948 * absco2_var_952 + colo3_var_924(jl_var_964, lay_var_944) * abso3_var_951 + coln2o_var_925(jl_var_964, lay_var_944) * absn2o_var_953 + wx_var_914(jl_var_964, 3, lay_var_944) * cfc12_var_216(ig_var_943) + wx_var_914(jl_var_964, 4, lay_var_944) * cfc22adj(ig_var_943)
        fracs_var_932(jl_var_964, 88 + ig_var_943, lay_var_944) = fracrefa_var_214(ig_var_943)
      END DO
    END DO
  END DO
  DO lay_var_944 = laytrop_max_var_955 + 1, klev_var_912
    DO jl_var_964 = kidia_var_910, kfdia_var_911
      chi_co2_var_945 = colco2_var_926(jl_var_964, lay_var_944) / coldry_var_927(jl_var_964, lay_var_944)
      ratco2_var_946 = 1D+20 * chi_co2_var_945 / chi_mls(2, jp_var_920(jl_var_964, lay_var_944) + 1)
      IF (ratco2_var_946 .GT. 3.0D0) THEN
        adjfac_var_947 = 2.0D0 + (ratco2_var_946 - 2.0D0) ** 0.65D0
        adjcolco2_var_948 = adjfac_var_947 * chi_mls(2, jp_var_920(jl_var_964, lay_var_944) + 1) * coldry_var_927(jl_var_964, lay_var_944) * 1D-20
      ELSE
        adjcolco2_var_948 = colco2_var_926(jl_var_964, lay_var_944)
      END IF
      ind0_var_938 = ((jp_var_920(jl_var_964, lay_var_944) - 13) * 5 + (jt_var_921(jl_var_964, lay_var_944) - 1)) * nspb_var_237(8) + 1
      ind1_var_939 = ((jp_var_920(jl_var_964, lay_var_944) - 12) * 5 + (jt1_var_922(jl_var_964, lay_var_944) - 1)) * nspb_var_237(8) + 1
      indm_var_942 = indminor_var_937(jl_var_964, lay_var_944)
      DO ig_var_943 = 1, 8
        absco2_var_952 = (kb_mco2_var_222(indm_var_942, ig_var_943) + minorfrac_var_936(jl_var_964, lay_var_944) * (kb_mco2_var_222(indm_var_942 + 1, ig_var_943) - kb_mco2_var_222(indm_var_942, ig_var_943)))
        absn2o_var_953 = (kb_mn2o_var_223(indm_var_942, ig_var_943) + minorfrac_var_936(jl_var_964, lay_var_944) * (kb_mn2o_var_223(indm_var_942 + 1, ig_var_943) - kb_mn2o_var_223(indm_var_942, ig_var_943)))
        taug_var_913(jl_var_964, 88 + ig_var_943, lay_var_944) = colo3_var_924(jl_var_964, lay_var_944) * (fac00_var_916(jl_var_964, lay_var_944) * absb_var_218(ind0_var_938, ig_var_943) + fac10_var_918(jl_var_964, lay_var_944) * absb_var_218(ind0_var_938 + 1, ig_var_943) + fac01_var_917(jl_var_964, lay_var_944) * absb_var_218(ind1_var_939, ig_var_943) + fac11_var_919(jl_var_964, lay_var_944) * absb_var_218(ind1_var_939 + 1, ig_var_943)) + adjcolco2_var_948 * absco2_var_952 + coln2o_var_925(jl_var_964, lay_var_944) * absn2o_var_953 + wx_var_914(jl_var_964, 3, lay_var_944) * cfc12_var_216(ig_var_943) + wx_var_914(jl_var_964, 4, lay_var_944) * cfc22adj(ig_var_943)
        fracs_var_932(jl_var_964, 88 + ig_var_943, lay_var_944) = fracrefb_var_215(ig_var_943)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_955 /= laytrop_min_var_954) THEN
    DO lay_var_944 = laytrop_min_var_954 + 1, laytrop_max_var_955
      ixc0_var_961 = ixc_var_956(lay_var_944)
      DO ixp_var_962 = 1, ixc0_var_961
        jl_var_964 = ixlow_var_957(ixp_var_962, lay_var_944)
        chi_co2_var_945 = colco2_var_926(jl_var_964, lay_var_944) / (coldry_var_927(jl_var_964, lay_var_944))
        ratco2_var_946 = 1D+20 * chi_co2_var_945 / chi_mls(2, jp_var_920(jl_var_964, lay_var_944) + 1)
        IF (ratco2_var_946 .GT. 3.0D0) THEN
          adjfac_var_947 = 2.0D0 + (ratco2_var_946 - 2.0D0) ** 0.65D0
          adjcolco2_var_948 = adjfac_var_947 * chi_mls(2, jp_var_920(jl_var_964, lay_var_944) + 1) * coldry_var_927(jl_var_964, lay_var_944) * 1D-20
        ELSE
          adjcolco2_var_948 = colco2_var_926(jl_var_964, lay_var_944)
        END IF
        ind0_var_938 = ((jp_var_920(jl_var_964, lay_var_944) - 1) * 5 + (jt_var_921(jl_var_964, lay_var_944) - 1)) * nspa_var_236(8) + 1
        ind1_var_939 = (jp_var_920(jl_var_964, lay_var_944) * 5 + (jt1_var_922(jl_var_964, lay_var_944) - 1)) * nspa_var_236(8) + 1
        inds_var_940 = indself_var_931(jl_var_964, lay_var_944)
        indf_var_941 = indfor_var_933(jl_var_964, lay_var_944)
        indm_var_942 = indminor_var_937(jl_var_964, lay_var_944)
        DO ig_var_943 = 1, 8
          tauself_var_950 = selffac_var_929(jl_var_964, lay_var_944) * (selfref_var_224(inds_var_940, ig_var_943) + selffrac_var_930(jl_var_964, lay_var_944) * (selfref_var_224(inds_var_940 + 1, ig_var_943) - selfref_var_224(inds_var_940, ig_var_943)))
          taufor_var_949 = forfac_var_935(jl_var_964, lay_var_944) * (forref_var_225(indf_var_941, ig_var_943) + forfrac_var_934(jl_var_964, lay_var_944) * (forref_var_225(indf_var_941 + 1, ig_var_943) - forref_var_225(indf_var_941, ig_var_943)))
          absco2_var_952 = (ka_mco2_var_219(indm_var_942, ig_var_943) + minorfrac_var_936(jl_var_964, lay_var_944) * (ka_mco2_var_219(indm_var_942 + 1, ig_var_943) - ka_mco2_var_219(indm_var_942, ig_var_943)))
          abso3_var_951 = (ka_mo3_var_221(indm_var_942, ig_var_943) + minorfrac_var_936(jl_var_964, lay_var_944) * (ka_mo3_var_221(indm_var_942 + 1, ig_var_943) - ka_mo3_var_221(indm_var_942, ig_var_943)))
          absn2o_var_953 = (ka_mn2o_var_220(indm_var_942, ig_var_943) + minorfrac_var_936(jl_var_964, lay_var_944) * (ka_mn2o_var_220(indm_var_942 + 1, ig_var_943) - ka_mn2o_var_220(indm_var_942, ig_var_943)))
          taug_var_913(jl_var_964, 88 + ig_var_943, lay_var_944) = colh2o_var_923(jl_var_964, lay_var_944) * (fac00_var_916(jl_var_964, lay_var_944) * absa_var_217(ind0_var_938, ig_var_943) + fac10_var_918(jl_var_964, lay_var_944) * absa_var_217(ind0_var_938 + 1, ig_var_943) + fac01_var_917(jl_var_964, lay_var_944) * absa_var_217(ind1_var_939, ig_var_943) + fac11_var_919(jl_var_964, lay_var_944) * absa_var_217(ind1_var_939 + 1, ig_var_943)) + tauself_var_950 + taufor_var_949 + adjcolco2_var_948 * absco2_var_952 + colo3_var_924(jl_var_964, lay_var_944) * abso3_var_951 + coln2o_var_925(jl_var_964, lay_var_944) * absn2o_var_953 + wx_var_914(jl_var_964, 3, lay_var_944) * cfc12_var_216(ig_var_943) + wx_var_914(jl_var_964, 4, lay_var_944) * cfc22adj(ig_var_943)
          fracs_var_932(jl_var_964, 88 + ig_var_943, lay_var_944) = fracrefa_var_214(ig_var_943)
        END DO
      END DO
      ixc0_var_961 = kfdia_var_911 - kidia_var_910 + 1 - ixc0_var_961
      DO ixp_var_962 = 1, ixc0_var_961
        jl_var_964 = ixhigh_var_958(ixp_var_962, lay_var_944)
        chi_co2_var_945 = colco2_var_926(jl_var_964, lay_var_944) / coldry_var_927(jl_var_964, lay_var_944)
        ratco2_var_946 = 1D+20 * chi_co2_var_945 / chi_mls(2, jp_var_920(jl_var_964, lay_var_944) + 1)
        IF (ratco2_var_946 .GT. 3.0D0) THEN
          adjfac_var_947 = 2.0D0 + (ratco2_var_946 - 2.0D0) ** 0.65D0
          adjcolco2_var_948 = adjfac_var_947 * chi_mls(2, jp_var_920(jl_var_964, lay_var_944) + 1) * coldry_var_927(jl_var_964, lay_var_944) * 1D-20
        ELSE
          adjcolco2_var_948 = colco2_var_926(jl_var_964, lay_var_944)
        END IF
        ind0_var_938 = ((jp_var_920(jl_var_964, lay_var_944) - 13) * 5 + (jt_var_921(jl_var_964, lay_var_944) - 1)) * nspb_var_237(8) + 1
        ind1_var_939 = ((jp_var_920(jl_var_964, lay_var_944) - 12) * 5 + (jt1_var_922(jl_var_964, lay_var_944) - 1)) * nspb_var_237(8) + 1
        indm_var_942 = indminor_var_937(jl_var_964, lay_var_944)
        DO ig_var_943 = 1, 8
          absco2_var_952 = (kb_mco2_var_222(indm_var_942, ig_var_943) + minorfrac_var_936(jl_var_964, lay_var_944) * (kb_mco2_var_222(indm_var_942 + 1, ig_var_943) - kb_mco2_var_222(indm_var_942, ig_var_943)))
          absn2o_var_953 = (kb_mn2o_var_223(indm_var_942, ig_var_943) + minorfrac_var_936(jl_var_964, lay_var_944) * (kb_mn2o_var_223(indm_var_942 + 1, ig_var_943) - kb_mn2o_var_223(indm_var_942, ig_var_943)))
          taug_var_913(jl_var_964, 88 + ig_var_943, lay_var_944) = colo3_var_924(jl_var_964, lay_var_944) * (fac00_var_916(jl_var_964, lay_var_944) * absb_var_218(ind0_var_938, ig_var_943) + fac10_var_918(jl_var_964, lay_var_944) * absb_var_218(ind0_var_938 + 1, ig_var_943) + fac01_var_917(jl_var_964, lay_var_944) * absb_var_218(ind1_var_939, ig_var_943) + fac11_var_919(jl_var_964, lay_var_944) * absb_var_218(ind1_var_939 + 1, ig_var_943)) + adjcolco2_var_948 * absco2_var_952 + coln2o_var_925(jl_var_964, lay_var_944) * absn2o_var_953 + wx_var_914(jl_var_964, 3, lay_var_944) * cfc12_var_216(ig_var_943) + wx_var_914(jl_var_964, 4, lay_var_944) * cfc22adj(ig_var_943)
          fracs_var_932(jl_var_964, 88 + ig_var_943, lay_var_944) = fracrefb_var_215(ig_var_943)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol8
SUBROUTINE rrtm_prepare_gases(kidia_var_966, kfdia_var_967, klon, klev_var_965, paph, pap, pth, pt, pq, pco2, pch4, pn2o, pno2, pc11, pc12, pc22, pcl4, pozn, pcoldry_var_968, pwbrodl, pwkl_var_969, pwx_var_970, pavel_var_971, ptavel_var_972, pz, ptz, kreflect)
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: klon
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_965
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_966, kfdia_var_967
  REAL(KIND = 8), INTENT(IN) :: paph(klon, klev_var_965 + 1)
  REAL(KIND = 8), INTENT(IN) :: pap(klon, klev_var_965)
  REAL(KIND = 8), INTENT(IN) :: pth(klon, klev_var_965 + 1)
  REAL(KIND = 8), INTENT(IN) :: pt(klon, klev_var_965)
  REAL(KIND = 8), INTENT(IN) :: pq(klon, klev_var_965)
  REAL(KIND = 8), INTENT(IN) :: pco2(klon, klev_var_965)
  REAL(KIND = 8), INTENT(IN) :: pch4(klon, klev_var_965)
  REAL(KIND = 8), INTENT(IN) :: pn2o(klon, klev_var_965)
  REAL(KIND = 8), INTENT(IN) :: pno2(klon, klev_var_965)
  REAL(KIND = 8), INTENT(IN) :: pc11(klon, klev_var_965)
  REAL(KIND = 8), INTENT(IN) :: pc12(klon, klev_var_965)
  REAL(KIND = 8), INTENT(IN) :: pc22(klon, klev_var_965)
  REAL(KIND = 8), INTENT(IN) :: pcl4(klon, klev_var_965)
  REAL(KIND = 8), INTENT(IN) :: pozn(klon, klev_var_965)
  REAL(KIND = 8), INTENT(OUT) :: pcoldry_var_968(kidia_var_966 : kfdia_var_967, klev_var_965)
  REAL(KIND = 8), INTENT(OUT) :: pwbrodl(kidia_var_966 : kfdia_var_967, klev_var_965)
  REAL(KIND = 8), INTENT(OUT) :: pwkl_var_969(kidia_var_966 : kfdia_var_967, 35, klev_var_965)
  REAL(KIND = 8), INTENT(OUT) :: pwx_var_970(kidia_var_966 : kfdia_var_967, 4, klev_var_965)
  REAL(KIND = 8), INTENT(OUT) :: pavel_var_971(kidia_var_966 : kfdia_var_967, klev_var_965)
  REAL(KIND = 8), INTENT(OUT) :: ptavel_var_972(kidia_var_966 : kfdia_var_967, klev_var_965)
  REAL(KIND = 8), INTENT(OUT) :: pz(kidia_var_966 : kfdia_var_967, 0 : klev_var_965)
  REAL(KIND = 8), INTENT(OUT) :: ptz(kidia_var_966 : kfdia_var_967, 0 : klev_var_965)
  INTEGER(KIND = 4), INTENT(OUT) :: kreflect(kidia_var_966 : kfdia_var_967)
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
  INTEGER(KIND = 4) :: iatm, jmol, ixmax, j1, j2, jk_var_973, jl_var_974
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
  DO jl_var_974 = kidia_var_966, kfdia_var_967
    kreflect(jl_var_974) = 0
  END DO
  DO j2 = 1, klev_var_965
    DO j1 = 1, 35
      DO jl_var_974 = kidia_var_966, kfdia_var_967
        pwkl_var_969(jl_var_974, j1, j2) = 0.0D0
      END DO
    END DO
  END DO
  iatm = 0
  ixmax = 4
  DO jl_var_974 = kidia_var_966, kfdia_var_967
    pz(jl_var_974, 0) = paph(jl_var_974, klev_var_965 + 1) / 100.0D0
    ptz(jl_var_974, 0) = pth(jl_var_974, klev_var_965 + 1)
  END DO
  DO jk_var_973 = 1, klev_var_965
    DO jl_var_974 = kidia_var_966, kfdia_var_967
      pavel_var_971(jl_var_974, jk_var_973) = pap(jl_var_974, klev_var_965 - jk_var_973 + 1) / 100.0D0
      ptavel_var_972(jl_var_974, jk_var_973) = pt(jl_var_974, klev_var_965 - jk_var_973 + 1)
      pz(jl_var_974, jk_var_973) = paph(jl_var_974, klev_var_965 - jk_var_973 + 1) / 100.0D0
      ptz(jl_var_974, jk_var_973) = pth(jl_var_974, klev_var_965 - jk_var_973 + 1)
      pwkl_var_969(jl_var_974, 1, jk_var_973) = MAX(pq(jl_var_974, klev_var_965 - jk_var_973 + 1), 1E-15) * 28.97D0 / 18.0154D0
      pwkl_var_969(jl_var_974, 2, jk_var_973) = pco2(jl_var_974, klev_var_965 - jk_var_973 + 1) * 28.97D0 / 44.011D0
      pwkl_var_969(jl_var_974, 3, jk_var_973) = pozn(jl_var_974, klev_var_965 - jk_var_973 + 1) * 28.97D0 / 47.9982D0
      pwkl_var_969(jl_var_974, 4, jk_var_973) = pn2o(jl_var_974, klev_var_965 - jk_var_973 + 1) * 28.97D0 / 44.013D0
      pwkl_var_969(jl_var_974, 6, jk_var_973) = pch4(jl_var_974, klev_var_965 - jk_var_973 + 1) * 28.97D0 / 16.043D0
      pwkl_var_969(jl_var_974, 7, jk_var_973) = 0.209488D0
      zamm = (1.0D0 - pwkl_var_969(jl_var_974, 1, jk_var_973)) * 28.97D0 + pwkl_var_969(jl_var_974, 1, jk_var_973) * 18.0154D0
      pcoldry_var_968(jl_var_974, jk_var_973) = (pz(jl_var_974, jk_var_973 - 1) - pz(jl_var_974, jk_var_973)) * 1000.0D0 * 6.02214D+23 / (980.665D0 * zamm * (1.0D0 + pwkl_var_969(jl_var_974, 1, jk_var_973)))
    END DO
  END DO
  DO j2 = 1, klev_var_965
    DO j1 = 1, 4
      DO jl_var_974 = kidia_var_966, kfdia_var_967
        pwx_var_970(jl_var_974, j1, j2) = 0.0D0
      END DO
    END DO
  END DO
  DO jk_var_973 = 1, klev_var_965
    DO jl_var_974 = kidia_var_966, kfdia_var_967
      pwx_var_970(jl_var_974, 1, jk_var_973) = pcl4(jl_var_974, klev_var_965 - jk_var_973 + 1) * 28.97D0 / 153.823D0
      pwx_var_970(jl_var_974, 2, jk_var_973) = pc11(jl_var_974, klev_var_965 - jk_var_973 + 1) * 28.97D0 / 137.3686D0
      pwx_var_970(jl_var_974, 3, jk_var_973) = pc12(jl_var_974, klev_var_965 - jk_var_973 + 1) * 28.97D0 / 120.914D0
      pwx_var_970(jl_var_974, 4, jk_var_973) = pc22(jl_var_974, klev_var_965 - jk_var_973 + 1) * 28.97D0 / 86.469D0
      pwx_var_970(jl_var_974, 1, jk_var_973) = pcoldry_var_968(jl_var_974, jk_var_973) * pwx_var_970(jl_var_974, 1, jk_var_973) * 1D-20
      pwx_var_970(jl_var_974, 2, jk_var_973) = pcoldry_var_968(jl_var_974, jk_var_973) * pwx_var_970(jl_var_974, 2, jk_var_973) * 1D-20
      pwx_var_970(jl_var_974, 3, jk_var_973) = pcoldry_var_968(jl_var_974, jk_var_973) * pwx_var_970(jl_var_974, 3, jk_var_973) * 1D-20
      pwx_var_970(jl_var_974, 4, jk_var_973) = pcoldry_var_968(jl_var_974, jk_var_973) * pwx_var_970(jl_var_974, 4, jk_var_973) * 1D-20
      zsummol = 0.0D0
      DO jmol = 2, 7
        zsummol = 0.0D0 + pwkl_var_969(jl_var_974, jmol, jk_var_973)
      END DO
      pwbrodl(jl_var_974, jk_var_973) = pcoldry_var_968(jl_var_974, jk_var_973) * (1.0D0 - zsummol)
      DO jmol = 1, 7
        pwkl_var_969(jl_var_974, jmol, jk_var_973) = pcoldry_var_968(jl_var_974, jk_var_973) * pwkl_var_969(jl_var_974, jmol, jk_var_973)
      END DO
    END DO
  END DO
END SUBROUTINE rrtm_prepare_gases
SUBROUTINE srtm_gas_optical_depth(kidia_var_975, kfdia_var_976, klev_var_977, poneminus_var_978, prmu0_var_979, klaytrop_var_980, pcolch4_var_981, pcolco2_var_982, pcolh2o_var_983, pcolmol_var_984, pcolo2_var_985, pcolo3_var_986, pforfac_var_987, pforfrac_var_988, kindfor_var_989, pselffac_var_990, pselffrac_var_991, kindself_var_992, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, pod_var_1000, pssa, pincsol)
  USE yoesrtwn, ONLY: ngc
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_975, kfdia_var_976
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_977
  REAL(KIND = 8), INTENT(IN) :: poneminus_var_978(kidia_var_975 : kfdia_var_976)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_979(kidia_var_975 : kfdia_var_976)
  INTEGER(KIND = 4), INTENT(IN) :: klaytrop_var_980(kidia_var_975 : kfdia_var_976)
  REAL(KIND = 8), INTENT(IN) :: pcolch4_var_981(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pcolco2_var_982(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pcolh2o_var_983(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pcolmol_var_984(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pcolo2_var_985(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pcolo3_var_986(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pforfac_var_987(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pforfrac_var_988(kidia_var_975 : kfdia_var_976, klev_var_977)
  INTEGER(KIND = 4), INTENT(IN) :: kindfor_var_989(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pselffac_var_990(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pselffrac_var_991(kidia_var_975 : kfdia_var_976, klev_var_977)
  INTEGER(KIND = 4), INTENT(IN) :: kindself_var_992(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pfac00_var_993(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pfac01_var_994(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pfac10_var_995(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(IN) :: pfac11_var_996(kidia_var_975 : kfdia_var_976, klev_var_977)
  INTEGER(KIND = 4), INTENT(IN) :: kjp_var_997(kidia_var_975 : kfdia_var_976, klev_var_977)
  INTEGER(KIND = 4), INTENT(IN) :: kjt_var_998(kidia_var_975 : kfdia_var_976, klev_var_977)
  INTEGER(KIND = 4), INTENT(IN) :: kjt1_var_999(kidia_var_975 : kfdia_var_976, klev_var_977)
  REAL(KIND = 8), INTENT(INOUT) :: pod_var_1000(kidia_var_975 : kfdia_var_976, klev_var_977, 112)
  REAL(KIND = 8), INTENT(INOUT) :: pssa(kidia_var_975 : kfdia_var_976, klev_var_977, 112)
  REAL(KIND = 8), INTENT(INOUT) :: pincsol(kidia_var_975 : kfdia_var_976, 112)
  INTEGER(KIND = 4) :: ib1, ib2, ibm, igt, iw(kidia_var_975 : kfdia_var_976), jb_var_1001, jg_var_1002, jk_var_1003, jl_var_1004, icount
  REAL(KIND = 8) :: ztaug(kidia_var_975 : kfdia_var_976, klev_var_977, 16)
  REAL(KIND = 8) :: ztaur(kidia_var_975 : kfdia_var_976, klev_var_977, 16)
  REAL(KIND = 8) :: zsflxzen(kidia_var_975 : kfdia_var_976, 16)
  ib1 = 16
  ib2 = 29
  icount = 0
  DO jl_var_1004 = kidia_var_975, kfdia_var_976
    IF (prmu0_var_979(jl_var_1004) > 0.0D0) THEN
      icount = icount + 1
      iw(jl_var_1004) = 0
    END IF
  END DO
  IF (icount /= 0) THEN
    DO jb_var_1001 = ib1, ib2
      ibm = jb_var_1001 - 15
      igt = ngc(ibm)
      IF (jb_var_1001 == 16) THEN
        CALL srtm_taumol16(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, poneminus_var_978, pcolh2o_var_983, pcolch4_var_981, pcolmol_var_984, klaytrop_var_980, pselffac_var_990, pselffrac_var_991, kindself_var_992, pforfac_var_987, pforfrac_var_988, kindfor_var_989, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 17) THEN
        CALL srtm_taumol17(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, poneminus_var_978, pcolh2o_var_983, pcolco2_var_982, pcolmol_var_984, klaytrop_var_980, pselffac_var_990, pselffrac_var_991, kindself_var_992, pforfac_var_987, pforfrac_var_988, kindfor_var_989, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 18) THEN
        CALL srtm_taumol18(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, poneminus_var_978, pcolh2o_var_983, pcolch4_var_981, pcolmol_var_984, klaytrop_var_980, pselffac_var_990, pselffrac_var_991, kindself_var_992, pforfac_var_987, pforfrac_var_988, kindfor_var_989, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 19) THEN
        CALL srtm_taumol19(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, poneminus_var_978, pcolh2o_var_983, pcolco2_var_982, pcolmol_var_984, klaytrop_var_980, pselffac_var_990, pselffrac_var_991, kindself_var_992, pforfac_var_987, pforfrac_var_988, kindfor_var_989, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 20) THEN
        CALL srtm_taumol20(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, pcolh2o_var_983, pcolch4_var_981, pcolmol_var_984, klaytrop_var_980, pselffac_var_990, pselffrac_var_991, kindself_var_992, pforfac_var_987, pforfrac_var_988, kindfor_var_989, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 21) THEN
        CALL srtm_taumol21(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, poneminus_var_978, pcolh2o_var_983, pcolco2_var_982, pcolmol_var_984, klaytrop_var_980, pselffac_var_990, pselffrac_var_991, kindself_var_992, pforfac_var_987, pforfrac_var_988, kindfor_var_989, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 22) THEN
        CALL srtm_taumol22(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, poneminus_var_978, pcolh2o_var_983, pcolmol_var_984, pcolo2_var_985, klaytrop_var_980, pselffac_var_990, pselffrac_var_991, kindself_var_992, pforfac_var_987, pforfrac_var_988, kindfor_var_989, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 23) THEN
        CALL srtm_taumol23(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, pcolh2o_var_983, pcolmol_var_984, klaytrop_var_980, pselffac_var_990, pselffrac_var_991, kindself_var_992, pforfac_var_987, pforfrac_var_988, kindfor_var_989, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 24) THEN
        CALL srtm_taumol24(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, poneminus_var_978, pcolh2o_var_983, pcolmol_var_984, pcolo2_var_985, pcolo3_var_986, klaytrop_var_980, pselffac_var_990, pselffrac_var_991, kindself_var_992, pforfac_var_987, pforfrac_var_988, kindfor_var_989, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 25) THEN
        CALL srtm_taumol25(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, pcolh2o_var_983, pcolmol_var_984, pcolo3_var_986, klaytrop_var_980, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 26) THEN
        CALL srtm_taumol26(kidia_var_975, kfdia_var_976, klev_var_977, pcolmol_var_984, klaytrop_var_980, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 27) THEN
        CALL srtm_taumol27(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, pcolmol_var_984, pcolo3_var_986, klaytrop_var_980, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 28) THEN
        CALL srtm_taumol28(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, poneminus_var_978, pcolmol_var_984, pcolo2_var_985, pcolo3_var_986, klaytrop_var_980, zsflxzen, ztaug, ztaur, prmu0_var_979)
      ELSE IF (jb_var_1001 == 29) THEN
        CALL srtm_taumol29(kidia_var_975, kfdia_var_976, klev_var_977, pfac00_var_993, pfac01_var_994, pfac10_var_995, pfac11_var_996, kjp_var_997, kjt_var_998, kjt1_var_999, pcolh2o_var_983, pcolco2_var_982, pcolmol_var_984, klaytrop_var_980, pselffac_var_990, pselffrac_var_991, kindself_var_992, pforfac_var_987, pforfrac_var_988, kindfor_var_989, zsflxzen, ztaug, ztaur, prmu0_var_979)
      END IF
      DO jg_var_1002 = 1, igt
        DO jl_var_1004 = kidia_var_975, kfdia_var_976
          IF (prmu0_var_979(jl_var_1004) > 0.0D0) THEN
            iw(jl_var_1004) = iw(jl_var_1004) + 1
            pincsol(jl_var_1004, iw(jl_var_1004)) = zsflxzen(jl_var_1004, jg_var_1002)
          END IF
        END DO
        DO jk_var_1003 = 1, klev_var_977
          DO jl_var_1004 = kidia_var_975, kfdia_var_976
            IF (prmu0_var_979(jl_var_1004) > 0.0D0) THEN
              pod_var_1000(jl_var_1004, jk_var_1003, iw(jl_var_1004)) = ztaur(jl_var_1004, jk_var_1003, jg_var_1002) + ztaug(jl_var_1004, jk_var_1003, jg_var_1002)
              pssa(jl_var_1004, jk_var_1003, iw(jl_var_1004)) = ztaur(jl_var_1004, jk_var_1003, jg_var_1002) / pod_var_1000(jl_var_1004, jk_var_1003, iw(jl_var_1004))
            END IF
          END DO
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE srtm_gas_optical_depth
SUBROUTINE srtm_taumol28(kidia_var_1005, kfdia_var_1006, klev_var_1007, p_fac00_var_1008, p_fac01_var_1009, p_fac10_var_1010, p_fac11_var_1011, k_jp_var_1012, k_jt_var_1013, k_jt1_var_1014, p_oneminus_var_1015, p_colmol_var_1016, p_colo2_var_1017, p_colo3_var_1018, k_laytrop_var_1019, p_sfluxzen_var_1020, p_taug_var_1021, p_taur_var_1022, prmu0_var_1023)
  USE yoesrta28, ONLY: absa_var_323, absb_var_324, layreffr_var_322, rayl_var_320, sfluxrefc_var_325, strrat_var_321
  USE yoesrtwn, ONLY: nspa_var_333, nspb_var_334
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1005, kfdia_var_1006
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1007
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1008(kidia_var_1005 : kfdia_var_1006, klev_var_1007)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1009(kidia_var_1005 : kfdia_var_1006, klev_var_1007)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1010(kidia_var_1005 : kfdia_var_1006, klev_var_1007)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1011(kidia_var_1005 : kfdia_var_1006, klev_var_1007)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1012(kidia_var_1005 : kfdia_var_1006, klev_var_1007)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1013(kidia_var_1005 : kfdia_var_1006, klev_var_1007)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1014(kidia_var_1005 : kfdia_var_1006, klev_var_1007)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1015(kidia_var_1005 : kfdia_var_1006)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1016(kidia_var_1005 : kfdia_var_1006, klev_var_1007)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_1017(kidia_var_1005 : kfdia_var_1006, klev_var_1007)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1018(kidia_var_1005 : kfdia_var_1006, klev_var_1007)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1019(kidia_var_1005 : kfdia_var_1006)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1020(kidia_var_1005 : kfdia_var_1006, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1021(kidia_var_1005 : kfdia_var_1006, klev_var_1007, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1022(kidia_var_1005 : kfdia_var_1006, klev_var_1007, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1023(kidia_var_1005 : kfdia_var_1006)
  INTEGER(KIND = 4) :: ig_var_1024, ind0_var_1025, ind1_var_1026, js_var_1027, i_lay_var_1028, i_laysolfr_var_1029(kidia_var_1005 : kfdia_var_1006), i_nlayers_var_1030, iplon_var_1031
  INTEGER(KIND = 4) :: laytrop_min_var_1032, laytrop_max_var_1033
  REAL(KIND = 8) :: z_fs_var_1034, z_speccomb_var_1035, z_specmult_var_1036, z_specparm_var_1037, z_tauray_var_1038
  laytrop_min_var_1032 = MINVAL(k_laytrop_var_1019(kidia_var_1005 : kfdia_var_1006))
  laytrop_max_var_1033 = MAXVAL(k_laytrop_var_1019(kidia_var_1005 : kfdia_var_1006))
  i_nlayers_var_1030 = klev_var_1007
  DO iplon_var_1031 = kidia_var_1005, kfdia_var_1006
    i_laysolfr_var_1029(iplon_var_1031) = i_nlayers_var_1030
  END DO
  DO i_lay_var_1028 = 1, laytrop_min_var_1032
    DO iplon_var_1031 = kidia_var_1005, kfdia_var_1006
      z_speccomb_var_1035 = p_colo3_var_1018(iplon_var_1031, i_lay_var_1028) + strrat_var_321 * p_colo2_var_1017(iplon_var_1031, i_lay_var_1028)
      z_specparm_var_1037 = p_colo3_var_1018(iplon_var_1031, i_lay_var_1028) / z_speccomb_var_1035
      z_specparm_var_1037 = MIN(p_oneminus_var_1015(iplon_var_1031), z_specparm_var_1037)
      z_specmult_var_1036 = 8.0D0 * (z_specparm_var_1037)
      js_var_1027 = 1 + INT(z_specmult_var_1036)
      z_fs_var_1034 = z_specmult_var_1036 - AINT(z_specmult_var_1036)
      ind0_var_1025 = ((k_jp_var_1012(iplon_var_1031, i_lay_var_1028) - 1) * 5 + (k_jt_var_1013(iplon_var_1031, i_lay_var_1028) - 1)) * nspa_var_333(28) + js_var_1027
      ind1_var_1026 = (k_jp_var_1012(iplon_var_1031, i_lay_var_1028) * 5 + (k_jt1_var_1014(iplon_var_1031, i_lay_var_1028) - 1)) * nspa_var_333(28) + js_var_1027
      z_tauray_var_1038 = p_colmol_var_1016(iplon_var_1031, i_lay_var_1028) * rayl_var_320
      DO ig_var_1024 = 1, 6
        p_taug_var_1021(iplon_var_1031, i_lay_var_1028, ig_var_1024) = z_speccomb_var_1035 * ((1.0D0 - z_fs_var_1034) * (absa_var_323(ind0_var_1025, ig_var_1024) * p_fac00_var_1008(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind0_var_1025 + 9, ig_var_1024) * p_fac10_var_1010(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind1_var_1026, ig_var_1024) * p_fac01_var_1009(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind1_var_1026 + 9, ig_var_1024) * p_fac11_var_1011(iplon_var_1031, i_lay_var_1028)) + z_fs_var_1034 * (absa_var_323(ind0_var_1025 + 1, ig_var_1024) * p_fac00_var_1008(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind0_var_1025 + 10, ig_var_1024) * p_fac10_var_1010(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind1_var_1026 + 1, ig_var_1024) * p_fac01_var_1009(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind1_var_1026 + 10, ig_var_1024) * p_fac11_var_1011(iplon_var_1031, i_lay_var_1028)))
        p_taur_var_1022(iplon_var_1031, i_lay_var_1028, ig_var_1024) = z_tauray_var_1038
      END DO
    END DO
  END DO
  DO i_lay_var_1028 = laytrop_min_var_1032 + 1, laytrop_max_var_1033
    DO iplon_var_1031 = kidia_var_1005, kfdia_var_1006
      IF (i_lay_var_1028 <= k_laytrop_var_1019(iplon_var_1031)) THEN
        z_speccomb_var_1035 = p_colo3_var_1018(iplon_var_1031, i_lay_var_1028) + strrat_var_321 * p_colo2_var_1017(iplon_var_1031, i_lay_var_1028)
        z_specparm_var_1037 = p_colo3_var_1018(iplon_var_1031, i_lay_var_1028) / z_speccomb_var_1035
        z_specparm_var_1037 = MIN(p_oneminus_var_1015(iplon_var_1031), z_specparm_var_1037)
        z_specmult_var_1036 = 8.0D0 * (z_specparm_var_1037)
        js_var_1027 = 1 + INT(z_specmult_var_1036)
        z_fs_var_1034 = z_specmult_var_1036 - AINT(z_specmult_var_1036)
        ind0_var_1025 = ((k_jp_var_1012(iplon_var_1031, i_lay_var_1028) - 1) * 5 + (k_jt_var_1013(iplon_var_1031, i_lay_var_1028) - 1)) * nspa_var_333(28) + js_var_1027
        ind1_var_1026 = (k_jp_var_1012(iplon_var_1031, i_lay_var_1028) * 5 + (k_jt1_var_1014(iplon_var_1031, i_lay_var_1028) - 1)) * nspa_var_333(28) + js_var_1027
        z_tauray_var_1038 = p_colmol_var_1016(iplon_var_1031, i_lay_var_1028) * rayl_var_320
        DO ig_var_1024 = 1, 6
          p_taug_var_1021(iplon_var_1031, i_lay_var_1028, ig_var_1024) = z_speccomb_var_1035 * ((1.0D0 - z_fs_var_1034) * (absa_var_323(ind0_var_1025, ig_var_1024) * p_fac00_var_1008(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind0_var_1025 + 9, ig_var_1024) * p_fac10_var_1010(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind1_var_1026, ig_var_1024) * p_fac01_var_1009(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind1_var_1026 + 9, ig_var_1024) * p_fac11_var_1011(iplon_var_1031, i_lay_var_1028)) + z_fs_var_1034 * (absa_var_323(ind0_var_1025 + 1, ig_var_1024) * p_fac00_var_1008(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind0_var_1025 + 10, ig_var_1024) * p_fac10_var_1010(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind1_var_1026 + 1, ig_var_1024) * p_fac01_var_1009(iplon_var_1031, i_lay_var_1028) + absa_var_323(ind1_var_1026 + 10, ig_var_1024) * p_fac11_var_1011(iplon_var_1031, i_lay_var_1028)))
          p_taur_var_1022(iplon_var_1031, i_lay_var_1028, ig_var_1024) = z_tauray_var_1038
        END DO
      ELSE
        IF (k_jp_var_1012(iplon_var_1031, i_lay_var_1028 - 1) < layreffr_var_322 .AND. k_jp_var_1012(iplon_var_1031, i_lay_var_1028) >= layreffr_var_322) i_laysolfr_var_1029(iplon_var_1031) = i_lay_var_1028
        z_speccomb_var_1035 = p_colo3_var_1018(iplon_var_1031, i_lay_var_1028) + strrat_var_321 * p_colo2_var_1017(iplon_var_1031, i_lay_var_1028)
        z_specparm_var_1037 = p_colo3_var_1018(iplon_var_1031, i_lay_var_1028) / z_speccomb_var_1035
        z_specparm_var_1037 = MIN(p_oneminus_var_1015(iplon_var_1031), z_specparm_var_1037)
        z_specmult_var_1036 = 4.0D0 * (z_specparm_var_1037)
        js_var_1027 = 1 + INT(z_specmult_var_1036)
        z_fs_var_1034 = z_specmult_var_1036 - AINT(z_specmult_var_1036)
        ind0_var_1025 = ((k_jp_var_1012(iplon_var_1031, i_lay_var_1028) - 13) * 5 + (k_jt_var_1013(iplon_var_1031, i_lay_var_1028) - 1)) * nspb_var_334(28) + js_var_1027
        ind1_var_1026 = ((k_jp_var_1012(iplon_var_1031, i_lay_var_1028) - 12) * 5 + (k_jt1_var_1014(iplon_var_1031, i_lay_var_1028) - 1)) * nspb_var_334(28) + js_var_1027
        z_tauray_var_1038 = p_colmol_var_1016(iplon_var_1031, i_lay_var_1028) * rayl_var_320
        DO ig_var_1024 = 1, 6
          p_taug_var_1021(iplon_var_1031, i_lay_var_1028, ig_var_1024) = z_speccomb_var_1035 * ((1.0D0 - z_fs_var_1034) * (absb_var_324(ind0_var_1025, ig_var_1024) * p_fac00_var_1008(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind0_var_1025 + 5, ig_var_1024) * p_fac10_var_1010(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind1_var_1026, ig_var_1024) * p_fac01_var_1009(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind1_var_1026 + 5, ig_var_1024) * p_fac11_var_1011(iplon_var_1031, i_lay_var_1028)) + z_fs_var_1034 * (absb_var_324(ind0_var_1025 + 1, ig_var_1024) * p_fac00_var_1008(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind0_var_1025 + 6, ig_var_1024) * p_fac10_var_1010(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind1_var_1026 + 1, ig_var_1024) * p_fac01_var_1009(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind1_var_1026 + 6, ig_var_1024) * p_fac11_var_1011(iplon_var_1031, i_lay_var_1028)))
          IF (i_lay_var_1028 == i_laysolfr_var_1029(iplon_var_1031)) p_sfluxzen_var_1020(iplon_var_1031, ig_var_1024) = sfluxrefc_var_325(ig_var_1024, js_var_1027) + z_fs_var_1034 * (sfluxrefc_var_325(ig_var_1024, js_var_1027 + 1) - sfluxrefc_var_325(ig_var_1024, js_var_1027))
          p_taur_var_1022(iplon_var_1031, i_lay_var_1028, ig_var_1024) = z_tauray_var_1038
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1028 = laytrop_max_var_1033 + 1, i_nlayers_var_1030
    DO iplon_var_1031 = kidia_var_1005, kfdia_var_1006
      IF (k_jp_var_1012(iplon_var_1031, i_lay_var_1028 - 1) < layreffr_var_322 .AND. k_jp_var_1012(iplon_var_1031, i_lay_var_1028) >= layreffr_var_322) i_laysolfr_var_1029(iplon_var_1031) = i_lay_var_1028
      z_speccomb_var_1035 = p_colo3_var_1018(iplon_var_1031, i_lay_var_1028) + strrat_var_321 * p_colo2_var_1017(iplon_var_1031, i_lay_var_1028)
      z_specparm_var_1037 = p_colo3_var_1018(iplon_var_1031, i_lay_var_1028) / z_speccomb_var_1035
      z_specparm_var_1037 = MIN(p_oneminus_var_1015(iplon_var_1031), z_specparm_var_1037)
      z_specmult_var_1036 = 4.0D0 * (z_specparm_var_1037)
      js_var_1027 = 1 + INT(z_specmult_var_1036)
      z_fs_var_1034 = z_specmult_var_1036 - AINT(z_specmult_var_1036)
      ind0_var_1025 = ((k_jp_var_1012(iplon_var_1031, i_lay_var_1028) - 13) * 5 + (k_jt_var_1013(iplon_var_1031, i_lay_var_1028) - 1)) * nspb_var_334(28) + js_var_1027
      ind1_var_1026 = ((k_jp_var_1012(iplon_var_1031, i_lay_var_1028) - 12) * 5 + (k_jt1_var_1014(iplon_var_1031, i_lay_var_1028) - 1)) * nspb_var_334(28) + js_var_1027
      z_tauray_var_1038 = p_colmol_var_1016(iplon_var_1031, i_lay_var_1028) * rayl_var_320
      DO ig_var_1024 = 1, 6
        p_taug_var_1021(iplon_var_1031, i_lay_var_1028, ig_var_1024) = z_speccomb_var_1035 * ((1.0D0 - z_fs_var_1034) * (absb_var_324(ind0_var_1025, ig_var_1024) * p_fac00_var_1008(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind0_var_1025 + 5, ig_var_1024) * p_fac10_var_1010(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind1_var_1026, ig_var_1024) * p_fac01_var_1009(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind1_var_1026 + 5, ig_var_1024) * p_fac11_var_1011(iplon_var_1031, i_lay_var_1028)) + z_fs_var_1034 * (absb_var_324(ind0_var_1025 + 1, ig_var_1024) * p_fac00_var_1008(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind0_var_1025 + 6, ig_var_1024) * p_fac10_var_1010(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind1_var_1026 + 1, ig_var_1024) * p_fac01_var_1009(iplon_var_1031, i_lay_var_1028) + absb_var_324(ind1_var_1026 + 6, ig_var_1024) * p_fac11_var_1011(iplon_var_1031, i_lay_var_1028)))
        IF (i_lay_var_1028 == i_laysolfr_var_1029(iplon_var_1031)) p_sfluxzen_var_1020(iplon_var_1031, ig_var_1024) = sfluxrefc_var_325(ig_var_1024, js_var_1027) + z_fs_var_1034 * (sfluxrefc_var_325(ig_var_1024, js_var_1027 + 1) - sfluxrefc_var_325(ig_var_1024, js_var_1027))
        p_taur_var_1022(iplon_var_1031, i_lay_var_1028, ig_var_1024) = z_tauray_var_1038
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol28
SUBROUTINE srtm_taumol29(kidia_var_1039, kfdia_var_1040, klev_var_1041, p_fac00_var_1042, p_fac01_var_1043, p_fac10_var_1044, p_fac11_var_1045, k_jp_var_1046, k_jt_var_1047, k_jt1_var_1048, p_colh2o_var_1049, p_colco2_var_1050, p_colmol_var_1051, k_laytrop_var_1052, p_selffac_var_1053, p_selffrac_var_1054, k_indself_var_1055, p_forfac_var_1056, p_forfrac_var_1057, k_indfor_var_1058, p_sfluxzen_var_1061, p_taug_var_1062, p_taur_var_1063, prmu0_var_1064)
  USE yoesrtwn, ONLY: nspa_var_333, nspb_var_334
  USE yoesrta29, ONLY: absa_var_328, absb_var_329, absco2c, absh2oc, forrefc_var_331, layreffr_var_327, rayl_var_326, selfrefc_var_330, sfluxrefc_var_332
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1039, kfdia_var_1040
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1041
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1042(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1043(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1044(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1045(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1046(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1047(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1048(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1049(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1050(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1051(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1052(kidia_var_1039 : kfdia_var_1040)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1053(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1054(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1055(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1056(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1057(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1058(kidia_var_1039 : kfdia_var_1040, klev_var_1041)
  INTEGER(KIND = 4) :: laytrop_min_var_1059, laytrop_max_var_1060
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1061(kidia_var_1039 : kfdia_var_1040, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1062(kidia_var_1039 : kfdia_var_1040, klev_var_1041, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1063(kidia_var_1039 : kfdia_var_1040, klev_var_1041, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1064(kidia_var_1039 : kfdia_var_1040)
  INTEGER(KIND = 4) :: ig_var_1065, ind0_var_1066, ind1_var_1067, inds_var_1068, indf_var_1069, i_lay_var_1070, i_laysolfr_var_1071(kidia_var_1039 : kfdia_var_1040), i_nlayers_var_1072, iplon_var_1073
  REAL(KIND = 8) :: z_tauray_var_1074
  laytrop_min_var_1059 = MINVAL(k_laytrop_var_1052(kidia_var_1039 : kfdia_var_1040))
  laytrop_max_var_1060 = MAXVAL(k_laytrop_var_1052(kidia_var_1039 : kfdia_var_1040))
  i_nlayers_var_1072 = klev_var_1041
  DO iplon_var_1073 = kidia_var_1039, kfdia_var_1040
    i_laysolfr_var_1071(iplon_var_1073) = i_nlayers_var_1072
  END DO
  DO i_lay_var_1070 = 1, laytrop_min_var_1059
    DO iplon_var_1073 = kidia_var_1039, kfdia_var_1040
      ind0_var_1066 = ((k_jp_var_1046(iplon_var_1073, i_lay_var_1070) - 1) * 5 + (k_jt_var_1047(iplon_var_1073, i_lay_var_1070) - 1)) * nspa_var_333(29) + 1
      ind1_var_1067 = (k_jp_var_1046(iplon_var_1073, i_lay_var_1070) * 5 + (k_jt1_var_1048(iplon_var_1073, i_lay_var_1070) - 1)) * nspa_var_333(29) + 1
      inds_var_1068 = k_indself_var_1055(iplon_var_1073, i_lay_var_1070)
      indf_var_1069 = k_indfor_var_1058(iplon_var_1073, i_lay_var_1070)
      z_tauray_var_1074 = p_colmol_var_1051(iplon_var_1073, i_lay_var_1070) * rayl_var_326
      DO ig_var_1065 = 1, 12
        p_taug_var_1062(iplon_var_1073, i_lay_var_1070, ig_var_1065) = p_colh2o_var_1049(iplon_var_1073, i_lay_var_1070) * ((p_fac00_var_1042(iplon_var_1073, i_lay_var_1070) * absa_var_328(ind0_var_1066, ig_var_1065) + p_fac10_var_1044(iplon_var_1073, i_lay_var_1070) * absa_var_328(ind0_var_1066 + 1, ig_var_1065) + p_fac01_var_1043(iplon_var_1073, i_lay_var_1070) * absa_var_328(ind1_var_1067, ig_var_1065) + p_fac11_var_1045(iplon_var_1073, i_lay_var_1070) * absa_var_328(ind1_var_1067 + 1, ig_var_1065)) + p_selffac_var_1053(iplon_var_1073, i_lay_var_1070) * (selfrefc_var_330(inds_var_1068, ig_var_1065) + p_selffrac_var_1054(iplon_var_1073, i_lay_var_1070) * (selfrefc_var_330(inds_var_1068 + 1, ig_var_1065) - selfrefc_var_330(inds_var_1068, ig_var_1065))) + p_forfac_var_1056(iplon_var_1073, i_lay_var_1070) * (forrefc_var_331(indf_var_1069, ig_var_1065) + p_forfrac_var_1057(iplon_var_1073, i_lay_var_1070) * (forrefc_var_331(indf_var_1069 + 1, ig_var_1065) - forrefc_var_331(indf_var_1069, ig_var_1065)))) + p_colco2_var_1050(iplon_var_1073, i_lay_var_1070) * absco2c(ig_var_1065)
        p_taur_var_1063(iplon_var_1073, i_lay_var_1070, ig_var_1065) = z_tauray_var_1074
      END DO
    END DO
  END DO
  DO i_lay_var_1070 = laytrop_min_var_1059 + 1, laytrop_max_var_1060
    DO iplon_var_1073 = kidia_var_1039, kfdia_var_1040
      IF (i_lay_var_1070 <= k_laytrop_var_1052(iplon_var_1073)) THEN
        ind0_var_1066 = ((k_jp_var_1046(iplon_var_1073, i_lay_var_1070) - 1) * 5 + (k_jt_var_1047(iplon_var_1073, i_lay_var_1070) - 1)) * nspa_var_333(29) + 1
        ind1_var_1067 = (k_jp_var_1046(iplon_var_1073, i_lay_var_1070) * 5 + (k_jt1_var_1048(iplon_var_1073, i_lay_var_1070) - 1)) * nspa_var_333(29) + 1
        inds_var_1068 = k_indself_var_1055(iplon_var_1073, i_lay_var_1070)
        indf_var_1069 = k_indfor_var_1058(iplon_var_1073, i_lay_var_1070)
        z_tauray_var_1074 = p_colmol_var_1051(iplon_var_1073, i_lay_var_1070) * rayl_var_326
        DO ig_var_1065 = 1, 12
          p_taug_var_1062(iplon_var_1073, i_lay_var_1070, ig_var_1065) = p_colh2o_var_1049(iplon_var_1073, i_lay_var_1070) * ((p_fac00_var_1042(iplon_var_1073, i_lay_var_1070) * absa_var_328(ind0_var_1066, ig_var_1065) + p_fac10_var_1044(iplon_var_1073, i_lay_var_1070) * absa_var_328(ind0_var_1066 + 1, ig_var_1065) + p_fac01_var_1043(iplon_var_1073, i_lay_var_1070) * absa_var_328(ind1_var_1067, ig_var_1065) + p_fac11_var_1045(iplon_var_1073, i_lay_var_1070) * absa_var_328(ind1_var_1067 + 1, ig_var_1065)) + p_selffac_var_1053(iplon_var_1073, i_lay_var_1070) * (selfrefc_var_330(inds_var_1068, ig_var_1065) + p_selffrac_var_1054(iplon_var_1073, i_lay_var_1070) * (selfrefc_var_330(inds_var_1068 + 1, ig_var_1065) - selfrefc_var_330(inds_var_1068, ig_var_1065))) + p_forfac_var_1056(iplon_var_1073, i_lay_var_1070) * (forrefc_var_331(indf_var_1069, ig_var_1065) + p_forfrac_var_1057(iplon_var_1073, i_lay_var_1070) * (forrefc_var_331(indf_var_1069 + 1, ig_var_1065) - forrefc_var_331(indf_var_1069, ig_var_1065)))) + p_colco2_var_1050(iplon_var_1073, i_lay_var_1070) * absco2c(ig_var_1065)
          p_taur_var_1063(iplon_var_1073, i_lay_var_1070, ig_var_1065) = z_tauray_var_1074
        END DO
      ELSE
        IF (k_jp_var_1046(iplon_var_1073, i_lay_var_1070 - 1) < layreffr_var_327 .AND. k_jp_var_1046(iplon_var_1073, i_lay_var_1070) >= layreffr_var_327) i_laysolfr_var_1071(iplon_var_1073) = i_lay_var_1070
        ind0_var_1066 = ((k_jp_var_1046(iplon_var_1073, i_lay_var_1070) - 13) * 5 + (k_jt_var_1047(iplon_var_1073, i_lay_var_1070) - 1)) * nspb_var_334(29) + 1
        ind1_var_1067 = ((k_jp_var_1046(iplon_var_1073, i_lay_var_1070) - 12) * 5 + (k_jt1_var_1048(iplon_var_1073, i_lay_var_1070) - 1)) * nspb_var_334(29) + 1
        z_tauray_var_1074 = p_colmol_var_1051(iplon_var_1073, i_lay_var_1070) * rayl_var_326
        DO ig_var_1065 = 1, 12
          p_taug_var_1062(iplon_var_1073, i_lay_var_1070, ig_var_1065) = p_colco2_var_1050(iplon_var_1073, i_lay_var_1070) * (p_fac00_var_1042(iplon_var_1073, i_lay_var_1070) * absb_var_329(ind0_var_1066, ig_var_1065) + p_fac10_var_1044(iplon_var_1073, i_lay_var_1070) * absb_var_329(ind0_var_1066 + 1, ig_var_1065) + p_fac01_var_1043(iplon_var_1073, i_lay_var_1070) * absb_var_329(ind1_var_1067, ig_var_1065) + p_fac11_var_1045(iplon_var_1073, i_lay_var_1070) * absb_var_329(ind1_var_1067 + 1, ig_var_1065)) + p_colh2o_var_1049(iplon_var_1073, i_lay_var_1070) * absh2oc(ig_var_1065)
          IF (i_lay_var_1070 == i_laysolfr_var_1071(iplon_var_1073)) p_sfluxzen_var_1061(iplon_var_1073, ig_var_1065) = sfluxrefc_var_332(ig_var_1065)
          p_taur_var_1063(iplon_var_1073, i_lay_var_1070, ig_var_1065) = z_tauray_var_1074
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1070 = laytrop_max_var_1060 + 1, i_nlayers_var_1072
    DO iplon_var_1073 = kidia_var_1039, kfdia_var_1040
      IF (k_jp_var_1046(iplon_var_1073, i_lay_var_1070 - 1) < layreffr_var_327 .AND. k_jp_var_1046(iplon_var_1073, i_lay_var_1070) >= layreffr_var_327) i_laysolfr_var_1071(iplon_var_1073) = i_lay_var_1070
      ind0_var_1066 = ((k_jp_var_1046(iplon_var_1073, i_lay_var_1070) - 13) * 5 + (k_jt_var_1047(iplon_var_1073, i_lay_var_1070) - 1)) * nspb_var_334(29) + 1
      ind1_var_1067 = ((k_jp_var_1046(iplon_var_1073, i_lay_var_1070) - 12) * 5 + (k_jt1_var_1048(iplon_var_1073, i_lay_var_1070) - 1)) * nspb_var_334(29) + 1
      z_tauray_var_1074 = p_colmol_var_1051(iplon_var_1073, i_lay_var_1070) * rayl_var_326
      DO ig_var_1065 = 1, 12
        p_taug_var_1062(iplon_var_1073, i_lay_var_1070, ig_var_1065) = p_colco2_var_1050(iplon_var_1073, i_lay_var_1070) * (p_fac00_var_1042(iplon_var_1073, i_lay_var_1070) * absb_var_329(ind0_var_1066, ig_var_1065) + p_fac10_var_1044(iplon_var_1073, i_lay_var_1070) * absb_var_329(ind0_var_1066 + 1, ig_var_1065) + p_fac01_var_1043(iplon_var_1073, i_lay_var_1070) * absb_var_329(ind1_var_1067, ig_var_1065) + p_fac11_var_1045(iplon_var_1073, i_lay_var_1070) * absb_var_329(ind1_var_1067 + 1, ig_var_1065)) + p_colh2o_var_1049(iplon_var_1073, i_lay_var_1070) * absh2oc(ig_var_1065)
        IF (i_lay_var_1070 == i_laysolfr_var_1071(iplon_var_1073)) p_sfluxzen_var_1061(iplon_var_1073, ig_var_1065) = sfluxrefc_var_332(ig_var_1065)
        p_taur_var_1063(iplon_var_1073, i_lay_var_1070, ig_var_1065) = z_tauray_var_1074
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol29
SUBROUTINE srtm_taumol17(kidia_var_1075, kfdia_var_1076, klev_var_1077, p_fac00_var_1078, p_fac01_var_1079, p_fac10_var_1080, p_fac11_var_1081, k_jp_var_1082, k_jt_var_1083, k_jt1_var_1084, p_oneminus_var_1085, p_colh2o_var_1086, p_colco2_var_1087, p_colmol_var_1088, k_laytrop_var_1089, p_selffac_var_1090, p_selffrac_var_1091, k_indself_var_1092, p_forfac_var_1093, p_forfrac_var_1094, k_indfor_var_1095, p_sfluxzen_var_1096, p_taug_var_1097, p_taur_var_1098, prmu0_var_1099)
  USE yoesrta17, ONLY: absa_var_248, absb_var_249, forrefc_var_251, layreffr_var_247, rayl_var_245, selfrefc_var_250, sfluxrefc_var_252, strrat_var_246
  USE yoesrtwn, ONLY: nspa_var_333, nspb_var_334
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1075, kfdia_var_1076
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1077
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1078(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1079(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1080(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1081(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1082(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1083(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1084(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1085(kidia_var_1075 : kfdia_var_1076)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1086(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1087(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1088(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1089(kidia_var_1075 : kfdia_var_1076)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1090(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1091(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1092(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1093(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1094(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1095(kidia_var_1075 : kfdia_var_1076, klev_var_1077)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1096(kidia_var_1075 : kfdia_var_1076, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1097(kidia_var_1075 : kfdia_var_1076, klev_var_1077, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1098(kidia_var_1075 : kfdia_var_1076, klev_var_1077, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1099(kidia_var_1075 : kfdia_var_1076)
  INTEGER(KIND = 4) :: ig_var_1100, ind0_var_1101, ind1_var_1102, inds_var_1103, indf_var_1104, js_var_1105, i_lay_var_1106, i_laysolfr_var_1107(kidia_var_1075 : kfdia_var_1076), i_nlayers_var_1108, iplon_var_1109
  INTEGER(KIND = 4) :: laytrop_min_var_1110, laytrop_max_var_1111
  REAL(KIND = 8) :: z_fs_var_1112, z_speccomb_var_1113, z_specmult_var_1114, z_specparm_var_1115, z_tauray_var_1116
  laytrop_min_var_1110 = MINVAL(k_laytrop_var_1089(kidia_var_1075 : kfdia_var_1076))
  laytrop_max_var_1111 = MAXVAL(k_laytrop_var_1089(kidia_var_1075 : kfdia_var_1076))
  i_nlayers_var_1108 = klev_var_1077
  DO iplon_var_1109 = kidia_var_1075, kfdia_var_1076
    i_laysolfr_var_1107(iplon_var_1109) = i_nlayers_var_1108
  END DO
  DO i_lay_var_1106 = 1, laytrop_min_var_1110
    DO iplon_var_1109 = kidia_var_1075, kfdia_var_1076
      z_speccomb_var_1113 = p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) + strrat_var_246 * p_colco2_var_1087(iplon_var_1109, i_lay_var_1106)
      z_specparm_var_1115 = p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) / z_speccomb_var_1113
      z_specparm_var_1115 = MIN(p_oneminus_var_1085(iplon_var_1109), z_specparm_var_1115)
      z_specmult_var_1114 = 8.0D0 * (z_specparm_var_1115)
      js_var_1105 = 1 + INT(z_specmult_var_1114)
      z_fs_var_1112 = z_specmult_var_1114 - AINT(z_specmult_var_1114)
      ind0_var_1101 = ((k_jp_var_1082(iplon_var_1109, i_lay_var_1106) - 1) * 5 + (k_jt_var_1083(iplon_var_1109, i_lay_var_1106) - 1)) * nspa_var_333(17) + js_var_1105
      ind1_var_1102 = (k_jp_var_1082(iplon_var_1109, i_lay_var_1106) * 5 + (k_jt1_var_1084(iplon_var_1109, i_lay_var_1106) - 1)) * nspa_var_333(17) + js_var_1105
      inds_var_1103 = k_indself_var_1092(iplon_var_1109, i_lay_var_1106)
      indf_var_1104 = k_indfor_var_1095(iplon_var_1109, i_lay_var_1106)
      z_tauray_var_1116 = p_colmol_var_1088(iplon_var_1109, i_lay_var_1106) * rayl_var_245
      DO ig_var_1100 = 1, 12
        p_taug_var_1097(iplon_var_1109, i_lay_var_1106, ig_var_1100) = z_speccomb_var_1113 * ((1.0D0 - z_fs_var_1112) * (absa_var_248(ind0_var_1101, ig_var_1100) * p_fac00_var_1078(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind0_var_1101 + 9, ig_var_1100) * p_fac10_var_1080(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind1_var_1102, ig_var_1100) * p_fac01_var_1079(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind1_var_1102 + 9, ig_var_1100) * p_fac11_var_1081(iplon_var_1109, i_lay_var_1106)) + z_fs_var_1112 * (absa_var_248(ind0_var_1101 + 1, ig_var_1100) * p_fac00_var_1078(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind0_var_1101 + 10, ig_var_1100) * p_fac10_var_1080(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind1_var_1102 + 1, ig_var_1100) * p_fac01_var_1079(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind1_var_1102 + 10, ig_var_1100) * p_fac11_var_1081(iplon_var_1109, i_lay_var_1106))) + p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) * (p_selffac_var_1090(iplon_var_1109, i_lay_var_1106) * (selfrefc_var_250(inds_var_1103, ig_var_1100) + p_selffrac_var_1091(iplon_var_1109, i_lay_var_1106) * (selfrefc_var_250(inds_var_1103 + 1, ig_var_1100) - selfrefc_var_250(inds_var_1103, ig_var_1100))) + p_forfac_var_1093(iplon_var_1109, i_lay_var_1106) * (forrefc_var_251(indf_var_1104, ig_var_1100) + p_forfrac_var_1094(iplon_var_1109, i_lay_var_1106) * (forrefc_var_251(indf_var_1104 + 1, ig_var_1100) - forrefc_var_251(indf_var_1104, ig_var_1100))))
        p_taur_var_1098(iplon_var_1109, i_lay_var_1106, ig_var_1100) = z_tauray_var_1116
      END DO
    END DO
  END DO
  DO i_lay_var_1106 = laytrop_min_var_1110 + 1, laytrop_max_var_1111
    DO iplon_var_1109 = kidia_var_1075, kfdia_var_1076
      IF (i_lay_var_1106 <= k_laytrop_var_1089(iplon_var_1109)) THEN
        z_speccomb_var_1113 = p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) + strrat_var_246 * p_colco2_var_1087(iplon_var_1109, i_lay_var_1106)
        z_specparm_var_1115 = p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) / z_speccomb_var_1113
        z_specparm_var_1115 = MIN(p_oneminus_var_1085(iplon_var_1109), z_specparm_var_1115)
        z_specmult_var_1114 = 8.0D0 * (z_specparm_var_1115)
        js_var_1105 = 1 + INT(z_specmult_var_1114)
        z_fs_var_1112 = z_specmult_var_1114 - AINT(z_specmult_var_1114)
        ind0_var_1101 = ((k_jp_var_1082(iplon_var_1109, i_lay_var_1106) - 1) * 5 + (k_jt_var_1083(iplon_var_1109, i_lay_var_1106) - 1)) * nspa_var_333(17) + js_var_1105
        ind1_var_1102 = (k_jp_var_1082(iplon_var_1109, i_lay_var_1106) * 5 + (k_jt1_var_1084(iplon_var_1109, i_lay_var_1106) - 1)) * nspa_var_333(17) + js_var_1105
        inds_var_1103 = k_indself_var_1092(iplon_var_1109, i_lay_var_1106)
        indf_var_1104 = k_indfor_var_1095(iplon_var_1109, i_lay_var_1106)
        z_tauray_var_1116 = p_colmol_var_1088(iplon_var_1109, i_lay_var_1106) * rayl_var_245
        DO ig_var_1100 = 1, 12
          p_taug_var_1097(iplon_var_1109, i_lay_var_1106, ig_var_1100) = z_speccomb_var_1113 * ((1.0D0 - z_fs_var_1112) * (absa_var_248(ind0_var_1101, ig_var_1100) * p_fac00_var_1078(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind0_var_1101 + 9, ig_var_1100) * p_fac10_var_1080(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind1_var_1102, ig_var_1100) * p_fac01_var_1079(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind1_var_1102 + 9, ig_var_1100) * p_fac11_var_1081(iplon_var_1109, i_lay_var_1106)) + z_fs_var_1112 * (absa_var_248(ind0_var_1101 + 1, ig_var_1100) * p_fac00_var_1078(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind0_var_1101 + 10, ig_var_1100) * p_fac10_var_1080(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind1_var_1102 + 1, ig_var_1100) * p_fac01_var_1079(iplon_var_1109, i_lay_var_1106) + absa_var_248(ind1_var_1102 + 10, ig_var_1100) * p_fac11_var_1081(iplon_var_1109, i_lay_var_1106))) + p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) * (p_selffac_var_1090(iplon_var_1109, i_lay_var_1106) * (selfrefc_var_250(inds_var_1103, ig_var_1100) + p_selffrac_var_1091(iplon_var_1109, i_lay_var_1106) * (selfrefc_var_250(inds_var_1103 + 1, ig_var_1100) - selfrefc_var_250(inds_var_1103, ig_var_1100))) + p_forfac_var_1093(iplon_var_1109, i_lay_var_1106) * (forrefc_var_251(indf_var_1104, ig_var_1100) + p_forfrac_var_1094(iplon_var_1109, i_lay_var_1106) * (forrefc_var_251(indf_var_1104 + 1, ig_var_1100) - forrefc_var_251(indf_var_1104, ig_var_1100))))
          p_taur_var_1098(iplon_var_1109, i_lay_var_1106, ig_var_1100) = z_tauray_var_1116
        END DO
      ELSE
        IF (k_jp_var_1082(iplon_var_1109, i_lay_var_1106 - 1) < layreffr_var_247 .AND. k_jp_var_1082(iplon_var_1109, i_lay_var_1106) >= layreffr_var_247) i_laysolfr_var_1107(iplon_var_1109) = i_lay_var_1106
        z_speccomb_var_1113 = p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) + strrat_var_246 * p_colco2_var_1087(iplon_var_1109, i_lay_var_1106)
        z_specparm_var_1115 = p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) / z_speccomb_var_1113
        z_specparm_var_1115 = MIN(p_oneminus_var_1085(iplon_var_1109), z_specparm_var_1115)
        z_specmult_var_1114 = 4.0D0 * (z_specparm_var_1115)
        js_var_1105 = 1 + INT(z_specmult_var_1114)
        z_fs_var_1112 = z_specmult_var_1114 - AINT(z_specmult_var_1114)
        ind0_var_1101 = ((k_jp_var_1082(iplon_var_1109, i_lay_var_1106) - 13) * 5 + (k_jt_var_1083(iplon_var_1109, i_lay_var_1106) - 1)) * nspb_var_334(17) + js_var_1105
        ind1_var_1102 = ((k_jp_var_1082(iplon_var_1109, i_lay_var_1106) - 12) * 5 + (k_jt1_var_1084(iplon_var_1109, i_lay_var_1106) - 1)) * nspb_var_334(17) + js_var_1105
        indf_var_1104 = k_indfor_var_1095(iplon_var_1109, i_lay_var_1106)
        z_tauray_var_1116 = p_colmol_var_1088(iplon_var_1109, i_lay_var_1106) * rayl_var_245
        DO ig_var_1100 = 1, 12
          p_taug_var_1097(iplon_var_1109, i_lay_var_1106, ig_var_1100) = z_speccomb_var_1113 * ((1.0D0 - z_fs_var_1112) * (absb_var_249(ind0_var_1101, ig_var_1100) * p_fac00_var_1078(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind0_var_1101 + 5, ig_var_1100) * p_fac10_var_1080(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind1_var_1102, ig_var_1100) * p_fac01_var_1079(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind1_var_1102 + 5, ig_var_1100) * p_fac11_var_1081(iplon_var_1109, i_lay_var_1106)) + z_fs_var_1112 * (absb_var_249(ind0_var_1101 + 1, ig_var_1100) * p_fac00_var_1078(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind0_var_1101 + 6, ig_var_1100) * p_fac10_var_1080(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind1_var_1102 + 1, ig_var_1100) * p_fac01_var_1079(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind1_var_1102 + 6, ig_var_1100) * p_fac11_var_1081(iplon_var_1109, i_lay_var_1106))) + p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) * p_forfac_var_1093(iplon_var_1109, i_lay_var_1106) * (forrefc_var_251(indf_var_1104, ig_var_1100) + p_forfrac_var_1094(iplon_var_1109, i_lay_var_1106) * (forrefc_var_251(indf_var_1104 + 1, ig_var_1100) - forrefc_var_251(indf_var_1104, ig_var_1100)))
          IF (i_lay_var_1106 == i_laysolfr_var_1107(iplon_var_1109)) p_sfluxzen_var_1096(iplon_var_1109, ig_var_1100) = sfluxrefc_var_252(ig_var_1100, js_var_1105) + z_fs_var_1112 * (sfluxrefc_var_252(ig_var_1100, js_var_1105 + 1) - sfluxrefc_var_252(ig_var_1100, js_var_1105))
          p_taur_var_1098(iplon_var_1109, i_lay_var_1106, ig_var_1100) = z_tauray_var_1116
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1106 = laytrop_max_var_1111 + 1, i_nlayers_var_1108
    DO iplon_var_1109 = kidia_var_1075, kfdia_var_1076
      IF (k_jp_var_1082(iplon_var_1109, i_lay_var_1106 - 1) < layreffr_var_247 .AND. k_jp_var_1082(iplon_var_1109, i_lay_var_1106) >= layreffr_var_247) i_laysolfr_var_1107(iplon_var_1109) = i_lay_var_1106
      z_speccomb_var_1113 = p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) + strrat_var_246 * p_colco2_var_1087(iplon_var_1109, i_lay_var_1106)
      z_specparm_var_1115 = p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) / z_speccomb_var_1113
      z_specparm_var_1115 = MIN(p_oneminus_var_1085(iplon_var_1109), z_specparm_var_1115)
      z_specmult_var_1114 = 4.0D0 * (z_specparm_var_1115)
      js_var_1105 = 1 + INT(z_specmult_var_1114)
      z_fs_var_1112 = z_specmult_var_1114 - AINT(z_specmult_var_1114)
      ind0_var_1101 = ((k_jp_var_1082(iplon_var_1109, i_lay_var_1106) - 13) * 5 + (k_jt_var_1083(iplon_var_1109, i_lay_var_1106) - 1)) * nspb_var_334(17) + js_var_1105
      ind1_var_1102 = ((k_jp_var_1082(iplon_var_1109, i_lay_var_1106) - 12) * 5 + (k_jt1_var_1084(iplon_var_1109, i_lay_var_1106) - 1)) * nspb_var_334(17) + js_var_1105
      indf_var_1104 = k_indfor_var_1095(iplon_var_1109, i_lay_var_1106)
      z_tauray_var_1116 = p_colmol_var_1088(iplon_var_1109, i_lay_var_1106) * rayl_var_245
      DO ig_var_1100 = 1, 12
        p_taug_var_1097(iplon_var_1109, i_lay_var_1106, ig_var_1100) = z_speccomb_var_1113 * ((1.0D0 - z_fs_var_1112) * (absb_var_249(ind0_var_1101, ig_var_1100) * p_fac00_var_1078(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind0_var_1101 + 5, ig_var_1100) * p_fac10_var_1080(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind1_var_1102, ig_var_1100) * p_fac01_var_1079(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind1_var_1102 + 5, ig_var_1100) * p_fac11_var_1081(iplon_var_1109, i_lay_var_1106)) + z_fs_var_1112 * (absb_var_249(ind0_var_1101 + 1, ig_var_1100) * p_fac00_var_1078(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind0_var_1101 + 6, ig_var_1100) * p_fac10_var_1080(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind1_var_1102 + 1, ig_var_1100) * p_fac01_var_1079(iplon_var_1109, i_lay_var_1106) + absb_var_249(ind1_var_1102 + 6, ig_var_1100) * p_fac11_var_1081(iplon_var_1109, i_lay_var_1106))) + p_colh2o_var_1086(iplon_var_1109, i_lay_var_1106) * p_forfac_var_1093(iplon_var_1109, i_lay_var_1106) * (forrefc_var_251(indf_var_1104, ig_var_1100) + p_forfrac_var_1094(iplon_var_1109, i_lay_var_1106) * (forrefc_var_251(indf_var_1104 + 1, ig_var_1100) - forrefc_var_251(indf_var_1104, ig_var_1100)))
        IF (i_lay_var_1106 == i_laysolfr_var_1107(iplon_var_1109)) p_sfluxzen_var_1096(iplon_var_1109, ig_var_1100) = sfluxrefc_var_252(ig_var_1100, js_var_1105) + z_fs_var_1112 * (sfluxrefc_var_252(ig_var_1100, js_var_1105 + 1) - sfluxrefc_var_252(ig_var_1100, js_var_1105))
        p_taur_var_1098(iplon_var_1109, i_lay_var_1106, ig_var_1100) = z_tauray_var_1116
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol17
SUBROUTINE rrtm_setcoef_140gp(kidia_var_1117, kfdia_var_1118, klev_var_1119, p_coldry, p_wbroad, p_wkl, p_fac00_var_1120, p_fac01_var_1121, p_fac10_var_1122, p_fac11_var_1123, p_forfac_var_1124, p_forfrac_var_1125, k_indfor_var_1142, k_jp_var_1126, k_jt_var_1127, k_jt1_var_1128, p_colh2o_var_1129, p_colco2_var_1130, p_colo3_var_1131, p_coln2o, p_colch4_var_1132, p_colo2_var_1133, p_co2mult_var_1134, p_colbrd, k_laytrop_var_1135, k_layswtch_var_1136, k_laylow_var_1137, pavel_var_1138, p_tavel, p_selffac_var_1139, p_selffrac_var_1140, k_indself_var_1141, k_indminor, p_scaleminor, p_scaleminorn2, p_minorfrac, prat_h2oco2_var_1143, prat_h2oco2_1_var_1144, prat_h2oo3_var_1145, prat_h2oo3_1_var_1146, prat_h2on2o_var_1147, prat_h2on2o_1_var_1148, prat_h2och4_var_1149, prat_h2och4_1_var_1150, prat_n2oco2_var_1151, prat_n2oco2_1_var_1152, prat_o3co2_var_1153, prat_o3co2_1_var_1154)
  USE yoerrtrf, ONLY: chi_mls, preflog_var_234, tref_var_235
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1117
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1118
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1119
  REAL(KIND = 8), INTENT(IN) :: p_coldry(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(IN) :: p_wbroad(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_colbrd(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(IN) :: p_wkl(kidia_var_1117 : kfdia_var_1118, 35, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_fac00_var_1120(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_fac01_var_1121(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_fac10_var_1122(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_fac11_var_1123(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_forfac_var_1124(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_forfrac_var_1125(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jp_var_1126(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jt_var_1127(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jt1_var_1128(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_colh2o_var_1129(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_colco2_var_1130(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_colo3_var_1131(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_coln2o(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_colch4_var_1132(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_colo2_var_1133(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_co2mult_var_1134(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  INTEGER(KIND = 4), INTENT(OUT) :: k_laytrop_var_1135(kidia_var_1117 : kfdia_var_1118)
  INTEGER(KIND = 4), INTENT(OUT) :: k_layswtch_var_1136(kidia_var_1117 : kfdia_var_1118)
  INTEGER(KIND = 4), INTENT(OUT) :: k_laylow_var_1137(kidia_var_1117 : kfdia_var_1118)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1138(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(IN) :: p_tavel(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_selffac_var_1139(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_selffrac_var_1140(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indself_var_1141(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indfor_var_1142(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indminor(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_scaleminor(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_scaleminorn2(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: p_minorfrac(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  REAL(KIND = 8), INTENT(OUT) :: prat_h2oco2_var_1143(kidia_var_1117 : kfdia_var_1118, klev_var_1119), prat_h2oco2_1_var_1144(kidia_var_1117 : kfdia_var_1118, klev_var_1119), prat_h2oo3_var_1145(kidia_var_1117 : kfdia_var_1118, klev_var_1119), prat_h2oo3_1_var_1146(kidia_var_1117 : kfdia_var_1118, klev_var_1119), prat_h2on2o_var_1147(kidia_var_1117 : kfdia_var_1118, klev_var_1119), prat_h2on2o_1_var_1148(kidia_var_1117 : kfdia_var_1118, klev_var_1119), prat_h2och4_var_1149(kidia_var_1117 : kfdia_var_1118, klev_var_1119), prat_h2och4_1_var_1150(kidia_var_1117 : kfdia_var_1118, klev_var_1119), prat_n2oco2_var_1151(kidia_var_1117 : kfdia_var_1118, klev_var_1119), prat_n2oco2_1_var_1152(kidia_var_1117 : kfdia_var_1118, klev_var_1119), prat_o3co2_var_1153(kidia_var_1117 : kfdia_var_1118, klev_var_1119), prat_o3co2_1_var_1154(kidia_var_1117 : kfdia_var_1118, klev_var_1119)
  INTEGER(KIND = 4) :: jp1_var_1155, jlay
  INTEGER(KIND = 4) :: jlon_var_1156
  REAL(KIND = 8) :: z_co2reg_var_1157, z_compfp_var_1158, z_factor_var_1159, z_fp_var_1160, z_ft_var_1161, z_ft1_var_1162, z_plog_var_1163, z_scalefac_var_1164, z_stpfac_var_1165, z_water_var_1166
  DO jlon_var_1156 = kidia_var_1117, kfdia_var_1118
    z_stpfac_var_1165 = 0.29220138203356366D0
    k_laytrop_var_1135(jlon_var_1156) = 0
    k_layswtch_var_1136(jlon_var_1156) = 0
    k_laylow_var_1137(jlon_var_1156) = 0
    DO jlay = 1, klev_var_1119
      z_plog_var_1163 = LOG(pavel_var_1138(jlon_var_1156, jlay))
      k_jp_var_1126(jlon_var_1156, jlay) = INT(36.0D0 - 5 * (z_plog_var_1163 + 0.04D0))
      IF (k_jp_var_1126(jlon_var_1156, jlay) < 1) THEN
        k_jp_var_1126(jlon_var_1156, jlay) = 1
      ELSE IF (k_jp_var_1126(jlon_var_1156, jlay) > 58) THEN
        k_jp_var_1126(jlon_var_1156, jlay) = 58
      END IF
      jp1_var_1155 = k_jp_var_1126(jlon_var_1156, jlay) + 1
      z_fp_var_1160 = 5.0D0 * (preflog_var_234(k_jp_var_1126(jlon_var_1156, jlay)) - z_plog_var_1163)
      z_fp_var_1160 = MAX(- 1.0D0, MIN(1.0D0, z_fp_var_1160))
      k_jt_var_1127(jlon_var_1156, jlay) = INT(3.0D0 + (p_tavel(jlon_var_1156, jlay) - tref_var_235(k_jp_var_1126(jlon_var_1156, jlay))) / 15.0D0)
      IF (k_jt_var_1127(jlon_var_1156, jlay) < 1) THEN
        k_jt_var_1127(jlon_var_1156, jlay) = 1
      ELSE IF (k_jt_var_1127(jlon_var_1156, jlay) > 4) THEN
        k_jt_var_1127(jlon_var_1156, jlay) = 4
      END IF
      z_ft_var_1161 = ((p_tavel(jlon_var_1156, jlay) - tref_var_235(k_jp_var_1126(jlon_var_1156, jlay))) / 15.0D0) - REAL(k_jt_var_1127(jlon_var_1156, jlay) - 3)
      k_jt1_var_1128(jlon_var_1156, jlay) = INT(3.0D0 + (p_tavel(jlon_var_1156, jlay) - tref_var_235(jp1_var_1155)) / 15.0D0)
      IF (k_jt1_var_1128(jlon_var_1156, jlay) < 1) THEN
        k_jt1_var_1128(jlon_var_1156, jlay) = 1
      ELSE IF (k_jt1_var_1128(jlon_var_1156, jlay) > 4) THEN
        k_jt1_var_1128(jlon_var_1156, jlay) = 4
      END IF
      z_ft1_var_1162 = ((p_tavel(jlon_var_1156, jlay) - tref_var_235(jp1_var_1155)) / 15.0D0) - REAL(k_jt1_var_1128(jlon_var_1156, jlay) - 3)
      z_water_var_1166 = p_wkl(jlon_var_1156, 1, jlay) / p_coldry(jlon_var_1156, jlay)
      z_scalefac_var_1164 = pavel_var_1138(jlon_var_1156, jlay) * 0.29220138203356366D0 / p_tavel(jlon_var_1156, jlay)
      IF (z_plog_var_1163 > 4.56D0) THEN
        k_laytrop_var_1135(jlon_var_1156) = k_laytrop_var_1135(jlon_var_1156) + 1
        p_forfac_var_1124(jlon_var_1156, jlay) = z_scalefac_var_1164 / (1.0D0 + z_water_var_1166)
        z_factor_var_1159 = (332.0D0 - p_tavel(jlon_var_1156, jlay)) / 36.0D0
        k_indfor_var_1142(jlon_var_1156, jlay) = MIN(2, MAX(1, INT(z_factor_var_1159)))
        p_forfrac_var_1125(jlon_var_1156, jlay) = z_factor_var_1159 - REAL(k_indfor_var_1142(jlon_var_1156, jlay))
        p_selffac_var_1139(jlon_var_1156, jlay) = z_water_var_1166 * p_forfac_var_1124(jlon_var_1156, jlay)
        z_factor_var_1159 = (p_tavel(jlon_var_1156, jlay) - 188.0D0) / 7.2D0
        k_indself_var_1141(jlon_var_1156, jlay) = MIN(9, MAX(1, INT(z_factor_var_1159) - 7))
        p_selffrac_var_1140(jlon_var_1156, jlay) = z_factor_var_1159 - REAL(k_indself_var_1141(jlon_var_1156, jlay) + 7)
        p_scaleminor(jlon_var_1156, jlay) = pavel_var_1138(jlon_var_1156, jlay) / p_tavel(jlon_var_1156, jlay)
        p_scaleminorn2(jlon_var_1156, jlay) = (pavel_var_1138(jlon_var_1156, jlay) / p_tavel(jlon_var_1156, jlay)) * (p_wbroad(jlon_var_1156, jlay) / (p_coldry(jlon_var_1156, jlay) + p_wkl(jlon_var_1156, 1, jlay)))
        z_factor_var_1159 = (p_tavel(jlon_var_1156, jlay) - 180.8D0) / 7.2D0
        k_indminor(jlon_var_1156, jlay) = MIN(18, MAX(1, INT(z_factor_var_1159)))
        p_minorfrac(jlon_var_1156, jlay) = z_factor_var_1159 - REAL(k_indminor(jlon_var_1156, jlay))
        prat_h2oco2_var_1143(jlon_var_1156, jlay) = chi_mls(1, k_jp_var_1126(jlon_var_1156, jlay)) / chi_mls(2, k_jp_var_1126(jlon_var_1156, jlay))
        prat_h2oco2_1_var_1144(jlon_var_1156, jlay) = chi_mls(1, k_jp_var_1126(jlon_var_1156, jlay) + 1) / chi_mls(2, k_jp_var_1126(jlon_var_1156, jlay) + 1)
        prat_h2oo3_var_1145(jlon_var_1156, jlay) = chi_mls(1, k_jp_var_1126(jlon_var_1156, jlay)) / chi_mls(3, k_jp_var_1126(jlon_var_1156, jlay))
        prat_h2oo3_1_var_1146(jlon_var_1156, jlay) = chi_mls(1, k_jp_var_1126(jlon_var_1156, jlay) + 1) / chi_mls(3, k_jp_var_1126(jlon_var_1156, jlay) + 1)
        prat_h2on2o_var_1147(jlon_var_1156, jlay) = chi_mls(1, k_jp_var_1126(jlon_var_1156, jlay)) / chi_mls(4, k_jp_var_1126(jlon_var_1156, jlay))
        prat_h2on2o_1_var_1148(jlon_var_1156, jlay) = chi_mls(1, k_jp_var_1126(jlon_var_1156, jlay) + 1) / chi_mls(4, k_jp_var_1126(jlon_var_1156, jlay) + 1)
        prat_h2och4_var_1149(jlon_var_1156, jlay) = chi_mls(1, k_jp_var_1126(jlon_var_1156, jlay)) / chi_mls(6, k_jp_var_1126(jlon_var_1156, jlay))
        prat_h2och4_1_var_1150(jlon_var_1156, jlay) = chi_mls(1, k_jp_var_1126(jlon_var_1156, jlay) + 1) / chi_mls(6, k_jp_var_1126(jlon_var_1156, jlay) + 1)
        prat_n2oco2_var_1151(jlon_var_1156, jlay) = chi_mls(4, k_jp_var_1126(jlon_var_1156, jlay)) / chi_mls(2, k_jp_var_1126(jlon_var_1156, jlay))
        prat_n2oco2_1_var_1152(jlon_var_1156, jlay) = chi_mls(4, k_jp_var_1126(jlon_var_1156, jlay) + 1) / chi_mls(2, k_jp_var_1126(jlon_var_1156, jlay) + 1)
        p_colh2o_var_1129(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 1, jlay)
        p_colco2_var_1130(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 2, jlay)
        p_colo3_var_1131(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 3, jlay)
        p_coln2o(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 4, jlay)
        p_colch4_var_1132(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 6, jlay)
        p_colo2_var_1133(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 7, jlay)
        p_colbrd(jlon_var_1156, jlay) = 1D-20 * p_wbroad(jlon_var_1156, jlay)
        IF (p_colco2_var_1130(jlon_var_1156, jlay) == 0.0D0) p_colco2_var_1130(jlon_var_1156, jlay) = 1D-32 * p_coldry(jlon_var_1156, jlay)
        IF (p_coln2o(jlon_var_1156, jlay) == 0.0D0) p_coln2o(jlon_var_1156, jlay) = 1D-32 * p_coldry(jlon_var_1156, jlay)
        IF (p_colch4_var_1132(jlon_var_1156, jlay) == 0.0D0) p_colch4_var_1132(jlon_var_1156, jlay) = 1D-32 * p_coldry(jlon_var_1156, jlay)
        z_co2reg_var_1157 = 3.55D-24 * p_coldry(jlon_var_1156, jlay)
        p_co2mult_var_1134(jlon_var_1156, jlay) = (p_colco2_var_1130(jlon_var_1156, jlay) - z_co2reg_var_1157) * 272.63D0 * EXP(- 1919.4D0 / p_tavel(jlon_var_1156, jlay)) / (0.00087604D0 * p_tavel(jlon_var_1156, jlay))
      ELSE
        p_forfac_var_1124(jlon_var_1156, jlay) = z_scalefac_var_1164 / (1.0D0 + z_water_var_1166)
        z_factor_var_1159 = (p_tavel(jlon_var_1156, jlay) - 188.0D0) / 36.0D0
        k_indfor_var_1142(jlon_var_1156, jlay) = 3
        p_forfrac_var_1125(jlon_var_1156, jlay) = z_factor_var_1159 - 1.0D0
        p_selffac_var_1139(jlon_var_1156, jlay) = z_water_var_1166 * p_forfac_var_1124(jlon_var_1156, jlay)
        p_scaleminor(jlon_var_1156, jlay) = pavel_var_1138(jlon_var_1156, jlay) / p_tavel(jlon_var_1156, jlay)
        p_scaleminorn2(jlon_var_1156, jlay) = (pavel_var_1138(jlon_var_1156, jlay) / p_tavel(jlon_var_1156, jlay)) * (p_wbroad(jlon_var_1156, jlay) / (p_coldry(jlon_var_1156, jlay) + p_wkl(jlon_var_1156, 1, jlay)))
        z_factor_var_1159 = (p_tavel(jlon_var_1156, jlay) - 180.8D0) / 7.2D0
        k_indminor(jlon_var_1156, jlay) = MIN(18, MAX(1, INT(z_factor_var_1159)))
        p_minorfrac(jlon_var_1156, jlay) = z_factor_var_1159 - REAL(k_indminor(jlon_var_1156, jlay))
        prat_h2oco2_var_1143(jlon_var_1156, jlay) = chi_mls(1, k_jp_var_1126(jlon_var_1156, jlay)) / chi_mls(2, k_jp_var_1126(jlon_var_1156, jlay))
        prat_h2oco2_1_var_1144(jlon_var_1156, jlay) = chi_mls(1, k_jp_var_1126(jlon_var_1156, jlay) + 1) / chi_mls(2, k_jp_var_1126(jlon_var_1156, jlay) + 1)
        prat_o3co2_var_1153(jlon_var_1156, jlay) = chi_mls(3, k_jp_var_1126(jlon_var_1156, jlay)) / chi_mls(2, k_jp_var_1126(jlon_var_1156, jlay))
        prat_o3co2_1_var_1154(jlon_var_1156, jlay) = chi_mls(3, k_jp_var_1126(jlon_var_1156, jlay) + 1) / chi_mls(2, k_jp_var_1126(jlon_var_1156, jlay) + 1)
        p_colh2o_var_1129(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 1, jlay)
        p_colco2_var_1130(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 2, jlay)
        p_colo3_var_1131(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 3, jlay)
        p_coln2o(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 4, jlay)
        p_colch4_var_1132(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 6, jlay)
        p_colo2_var_1133(jlon_var_1156, jlay) = 1D-20 * p_wkl(jlon_var_1156, 7, jlay)
        p_colbrd(jlon_var_1156, jlay) = 1D-20 * p_wbroad(jlon_var_1156, jlay)
        IF (p_colco2_var_1130(jlon_var_1156, jlay) == 0.0D0) p_colco2_var_1130(jlon_var_1156, jlay) = 1D-32 * p_coldry(jlon_var_1156, jlay)
        IF (p_coln2o(jlon_var_1156, jlay) == 0.0D0) p_coln2o(jlon_var_1156, jlay) = 1D-32 * p_coldry(jlon_var_1156, jlay)
        IF (p_colch4_var_1132(jlon_var_1156, jlay) == 0.0D0) p_colch4_var_1132(jlon_var_1156, jlay) = 1D-32 * p_coldry(jlon_var_1156, jlay)
        z_co2reg_var_1157 = 3.55D-24 * p_coldry(jlon_var_1156, jlay)
        p_co2mult_var_1134(jlon_var_1156, jlay) = (p_colco2_var_1130(jlon_var_1156, jlay) - z_co2reg_var_1157) * 272.63D0 * EXP(- 1919.4D0 / p_tavel(jlon_var_1156, jlay)) / (0.00087604D0 * p_tavel(jlon_var_1156, jlay))
      END IF
      z_compfp_var_1158 = 1.0D0 - z_fp_var_1160
      p_fac10_var_1122(jlon_var_1156, jlay) = z_compfp_var_1158 * z_ft_var_1161
      p_fac00_var_1120(jlon_var_1156, jlay) = z_compfp_var_1158 * (1.0D0 - z_ft_var_1161)
      p_fac11_var_1123(jlon_var_1156, jlay) = z_fp_var_1160 * z_ft1_var_1162
      p_fac01_var_1121(jlon_var_1156, jlay) = z_fp_var_1160 * (1.0D0 - z_ft1_var_1162)
      p_selffac_var_1139(jlon_var_1156, jlay) = p_colh2o_var_1129(jlon_var_1156, jlay) * p_selffac_var_1139(jlon_var_1156, jlay)
      p_forfac_var_1124(jlon_var_1156, jlay) = p_colh2o_var_1129(jlon_var_1156, jlay) * p_forfac_var_1124(jlon_var_1156, jlay)
    END DO
    IF (k_laylow_var_1137(jlon_var_1156) == 0) k_laylow_var_1137(jlon_var_1156) = 1
  END DO
END SUBROUTINE rrtm_setcoef_140gp
SUBROUTINE srtm_taumol16(kidia_var_1167, kfdia_var_1168, klev_var_1169, p_fac00_var_1170, p_fac01_var_1171, p_fac10_var_1172, p_fac11_var_1173, k_jp_var_1174, k_jt_var_1175, k_jt1_var_1176, p_oneminus_var_1177, p_colh2o_var_1178, p_colch4_var_1179, p_colmol_var_1180, k_laytrop_var_1181, p_selffac_var_1182, p_selffrac_var_1183, k_indself_var_1184, p_forfac_var_1185, p_forfrac_var_1186, k_indfor_var_1187, p_sfluxzen_var_1188, p_taug_var_1189, p_taur_var_1190, prmu0_var_1191)
  USE yoesrta16, ONLY: absa_var_240, absb_var_241, forrefc_var_243, layreffr_var_239, rayl_var_238, selfrefc_var_242, sfluxrefc_var_244, strrat1
  USE yoesrtwn, ONLY: nspa_var_333, nspb_var_334
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1167, kfdia_var_1168
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1169
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1170(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1171(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1172(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1173(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1174(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1175(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1176(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1177(kidia_var_1167 : kfdia_var_1168)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1178(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_1179(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1180(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1181(kidia_var_1167 : kfdia_var_1168)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1182(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1183(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1184(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1185(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1186(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1187(kidia_var_1167 : kfdia_var_1168, klev_var_1169)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1188(kidia_var_1167 : kfdia_var_1168, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1189(kidia_var_1167 : kfdia_var_1168, klev_var_1169, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1190(kidia_var_1167 : kfdia_var_1168, klev_var_1169, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1191(kidia_var_1167 : kfdia_var_1168)
  INTEGER(KIND = 4) :: ig_var_1192, ind0_var_1193, ind1_var_1194, inds_var_1195, indf_var_1196, js_var_1197, i_lay_var_1198, i_laysolfr_var_1199(kidia_var_1167 : kfdia_var_1168), i_nlayers_var_1200, iplon_var_1201
  INTEGER(KIND = 4) :: laytrop_min_var_1202, laytrop_max_var_1203
  REAL(KIND = 8) :: z_fs_var_1204, z_speccomb_var_1205, z_specmult_var_1206, z_specparm_var_1207, z_tauray_var_1208
  laytrop_min_var_1202 = MINVAL(k_laytrop_var_1181(kidia_var_1167 : kfdia_var_1168))
  laytrop_max_var_1203 = MAXVAL(k_laytrop_var_1181(kidia_var_1167 : kfdia_var_1168))
  i_nlayers_var_1200 = klev_var_1169
  DO iplon_var_1201 = kidia_var_1167, kfdia_var_1168
    i_laysolfr_var_1199(iplon_var_1201) = i_nlayers_var_1200
  END DO
  DO i_lay_var_1198 = 1, laytrop_min_var_1202
    DO iplon_var_1201 = kidia_var_1167, kfdia_var_1168
      z_speccomb_var_1205 = p_colh2o_var_1178(iplon_var_1201, i_lay_var_1198) + strrat1 * p_colch4_var_1179(iplon_var_1201, i_lay_var_1198)
      z_specparm_var_1207 = p_colh2o_var_1178(iplon_var_1201, i_lay_var_1198) / z_speccomb_var_1205
      z_specparm_var_1207 = MIN(p_oneminus_var_1177(iplon_var_1201), z_specparm_var_1207)
      z_specmult_var_1206 = 8.0D0 * (z_specparm_var_1207)
      js_var_1197 = 1 + INT(z_specmult_var_1206)
      z_fs_var_1204 = z_specmult_var_1206 - AINT(z_specmult_var_1206)
      ind0_var_1193 = ((k_jp_var_1174(iplon_var_1201, i_lay_var_1198) - 1) * 5 + (k_jt_var_1175(iplon_var_1201, i_lay_var_1198) - 1)) * nspa_var_333(16) + js_var_1197
      ind1_var_1194 = (k_jp_var_1174(iplon_var_1201, i_lay_var_1198) * 5 + (k_jt1_var_1176(iplon_var_1201, i_lay_var_1198) - 1)) * nspa_var_333(16) + js_var_1197
      inds_var_1195 = k_indself_var_1184(iplon_var_1201, i_lay_var_1198)
      indf_var_1196 = k_indfor_var_1187(iplon_var_1201, i_lay_var_1198)
      z_tauray_var_1208 = p_colmol_var_1180(iplon_var_1201, i_lay_var_1198) * rayl_var_238
      DO ig_var_1192 = 1, 6
        p_taug_var_1189(iplon_var_1201, i_lay_var_1198, ig_var_1192) = z_speccomb_var_1205 * ((1.0D0 - z_fs_var_1204) * (absa_var_240(ind0_var_1193, ig_var_1192) * p_fac00_var_1170(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind0_var_1193 + 9, ig_var_1192) * p_fac10_var_1172(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind1_var_1194, ig_var_1192) * p_fac01_var_1171(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind1_var_1194 + 9, ig_var_1192) * p_fac11_var_1173(iplon_var_1201, i_lay_var_1198)) + z_fs_var_1204 * (absa_var_240(ind0_var_1193 + 1, ig_var_1192) * p_fac00_var_1170(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind0_var_1193 + 10, ig_var_1192) * p_fac10_var_1172(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind1_var_1194 + 1, ig_var_1192) * p_fac01_var_1171(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind1_var_1194 + 10, ig_var_1192) * p_fac11_var_1173(iplon_var_1201, i_lay_var_1198))) + p_colh2o_var_1178(iplon_var_1201, i_lay_var_1198) * (p_selffac_var_1182(iplon_var_1201, i_lay_var_1198) * (selfrefc_var_242(inds_var_1195, ig_var_1192) + p_selffrac_var_1183(iplon_var_1201, i_lay_var_1198) * (selfrefc_var_242(inds_var_1195 + 1, ig_var_1192) - selfrefc_var_242(inds_var_1195, ig_var_1192))) + p_forfac_var_1185(iplon_var_1201, i_lay_var_1198) * (forrefc_var_243(indf_var_1196, ig_var_1192) + p_forfrac_var_1186(iplon_var_1201, i_lay_var_1198) * (forrefc_var_243(indf_var_1196 + 1, ig_var_1192) - forrefc_var_243(indf_var_1196, ig_var_1192))))
        p_taur_var_1190(iplon_var_1201, i_lay_var_1198, ig_var_1192) = z_tauray_var_1208
      END DO
    END DO
  END DO
  DO i_lay_var_1198 = laytrop_min_var_1202 + 1, laytrop_max_var_1203
    DO iplon_var_1201 = kidia_var_1167, kfdia_var_1168
      IF (i_lay_var_1198 <= k_laytrop_var_1181(iplon_var_1201)) THEN
        z_speccomb_var_1205 = p_colh2o_var_1178(iplon_var_1201, i_lay_var_1198) + strrat1 * p_colch4_var_1179(iplon_var_1201, i_lay_var_1198)
        z_specparm_var_1207 = p_colh2o_var_1178(iplon_var_1201, i_lay_var_1198) / z_speccomb_var_1205
        z_specparm_var_1207 = MIN(p_oneminus_var_1177(iplon_var_1201), z_specparm_var_1207)
        z_specmult_var_1206 = 8.0D0 * (z_specparm_var_1207)
        js_var_1197 = 1 + INT(z_specmult_var_1206)
        z_fs_var_1204 = z_specmult_var_1206 - AINT(z_specmult_var_1206)
        ind0_var_1193 = ((k_jp_var_1174(iplon_var_1201, i_lay_var_1198) - 1) * 5 + (k_jt_var_1175(iplon_var_1201, i_lay_var_1198) - 1)) * nspa_var_333(16) + js_var_1197
        ind1_var_1194 = (k_jp_var_1174(iplon_var_1201, i_lay_var_1198) * 5 + (k_jt1_var_1176(iplon_var_1201, i_lay_var_1198) - 1)) * nspa_var_333(16) + js_var_1197
        inds_var_1195 = k_indself_var_1184(iplon_var_1201, i_lay_var_1198)
        indf_var_1196 = k_indfor_var_1187(iplon_var_1201, i_lay_var_1198)
        z_tauray_var_1208 = p_colmol_var_1180(iplon_var_1201, i_lay_var_1198) * rayl_var_238
        DO ig_var_1192 = 1, 6
          p_taug_var_1189(iplon_var_1201, i_lay_var_1198, ig_var_1192) = z_speccomb_var_1205 * ((1.0D0 - z_fs_var_1204) * (absa_var_240(ind0_var_1193, ig_var_1192) * p_fac00_var_1170(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind0_var_1193 + 9, ig_var_1192) * p_fac10_var_1172(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind1_var_1194, ig_var_1192) * p_fac01_var_1171(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind1_var_1194 + 9, ig_var_1192) * p_fac11_var_1173(iplon_var_1201, i_lay_var_1198)) + z_fs_var_1204 * (absa_var_240(ind0_var_1193 + 1, ig_var_1192) * p_fac00_var_1170(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind0_var_1193 + 10, ig_var_1192) * p_fac10_var_1172(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind1_var_1194 + 1, ig_var_1192) * p_fac01_var_1171(iplon_var_1201, i_lay_var_1198) + absa_var_240(ind1_var_1194 + 10, ig_var_1192) * p_fac11_var_1173(iplon_var_1201, i_lay_var_1198))) + p_colh2o_var_1178(iplon_var_1201, i_lay_var_1198) * (p_selffac_var_1182(iplon_var_1201, i_lay_var_1198) * (selfrefc_var_242(inds_var_1195, ig_var_1192) + p_selffrac_var_1183(iplon_var_1201, i_lay_var_1198) * (selfrefc_var_242(inds_var_1195 + 1, ig_var_1192) - selfrefc_var_242(inds_var_1195, ig_var_1192))) + p_forfac_var_1185(iplon_var_1201, i_lay_var_1198) * (forrefc_var_243(indf_var_1196, ig_var_1192) + p_forfrac_var_1186(iplon_var_1201, i_lay_var_1198) * (forrefc_var_243(indf_var_1196 + 1, ig_var_1192) - forrefc_var_243(indf_var_1196, ig_var_1192))))
          p_taur_var_1190(iplon_var_1201, i_lay_var_1198, ig_var_1192) = z_tauray_var_1208
        END DO
      ELSE
        IF (k_jp_var_1174(iplon_var_1201, i_lay_var_1198 - 1) < layreffr_var_239 .AND. k_jp_var_1174(iplon_var_1201, i_lay_var_1198) >= layreffr_var_239) i_laysolfr_var_1199(iplon_var_1201) = i_lay_var_1198
        ind0_var_1193 = ((k_jp_var_1174(iplon_var_1201, i_lay_var_1198) - 13) * 5 + (k_jt_var_1175(iplon_var_1201, i_lay_var_1198) - 1)) * nspb_var_334(16) + 1
        ind1_var_1194 = ((k_jp_var_1174(iplon_var_1201, i_lay_var_1198) - 12) * 5 + (k_jt1_var_1176(iplon_var_1201, i_lay_var_1198) - 1)) * nspb_var_334(16) + 1
        z_tauray_var_1208 = p_colmol_var_1180(iplon_var_1201, i_lay_var_1198) * rayl_var_238
        DO ig_var_1192 = 1, 6
          p_taug_var_1189(iplon_var_1201, i_lay_var_1198, ig_var_1192) = p_colch4_var_1179(iplon_var_1201, i_lay_var_1198) * (p_fac00_var_1170(iplon_var_1201, i_lay_var_1198) * absb_var_241(ind0_var_1193, ig_var_1192) + p_fac10_var_1172(iplon_var_1201, i_lay_var_1198) * absb_var_241(ind0_var_1193 + 1, ig_var_1192) + p_fac01_var_1171(iplon_var_1201, i_lay_var_1198) * absb_var_241(ind1_var_1194, ig_var_1192) + p_fac11_var_1173(iplon_var_1201, i_lay_var_1198) * absb_var_241(ind1_var_1194 + 1, ig_var_1192))
          IF (i_lay_var_1198 == i_laysolfr_var_1199(iplon_var_1201)) p_sfluxzen_var_1188(iplon_var_1201, ig_var_1192) = sfluxrefc_var_244(ig_var_1192)
          p_taur_var_1190(iplon_var_1201, i_lay_var_1198, ig_var_1192) = z_tauray_var_1208
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1198 = laytrop_max_var_1203 + 1, i_nlayers_var_1200
    DO iplon_var_1201 = kidia_var_1167, kfdia_var_1168
      IF (k_jp_var_1174(iplon_var_1201, i_lay_var_1198 - 1) < layreffr_var_239 .AND. k_jp_var_1174(iplon_var_1201, i_lay_var_1198) >= layreffr_var_239) i_laysolfr_var_1199(iplon_var_1201) = i_lay_var_1198
      ind0_var_1193 = ((k_jp_var_1174(iplon_var_1201, i_lay_var_1198) - 13) * 5 + (k_jt_var_1175(iplon_var_1201, i_lay_var_1198) - 1)) * nspb_var_334(16) + 1
      ind1_var_1194 = ((k_jp_var_1174(iplon_var_1201, i_lay_var_1198) - 12) * 5 + (k_jt1_var_1176(iplon_var_1201, i_lay_var_1198) - 1)) * nspb_var_334(16) + 1
      z_tauray_var_1208 = p_colmol_var_1180(iplon_var_1201, i_lay_var_1198) * rayl_var_238
      DO ig_var_1192 = 1, 6
        p_taug_var_1189(iplon_var_1201, i_lay_var_1198, ig_var_1192) = p_colch4_var_1179(iplon_var_1201, i_lay_var_1198) * (p_fac00_var_1170(iplon_var_1201, i_lay_var_1198) * absb_var_241(ind0_var_1193, ig_var_1192) + p_fac10_var_1172(iplon_var_1201, i_lay_var_1198) * absb_var_241(ind0_var_1193 + 1, ig_var_1192) + p_fac01_var_1171(iplon_var_1201, i_lay_var_1198) * absb_var_241(ind1_var_1194, ig_var_1192) + p_fac11_var_1173(iplon_var_1201, i_lay_var_1198) * absb_var_241(ind1_var_1194 + 1, ig_var_1192))
        IF (i_lay_var_1198 == i_laysolfr_var_1199(iplon_var_1201)) p_sfluxzen_var_1188(iplon_var_1201, ig_var_1192) = sfluxrefc_var_244(ig_var_1192)
        p_taur_var_1190(iplon_var_1201, i_lay_var_1198, ig_var_1192) = z_tauray_var_1208
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol16
SUBROUTINE srtm_taumol27(kidia_var_1209, kfdia_var_1210, klev_var_1211, p_fac00_var_1212, p_fac01_var_1213, p_fac10_var_1214, p_fac11_var_1215, k_jp_var_1216, k_jt_var_1217, k_jt1_var_1218, p_colmol_var_1219, p_colo3_var_1220, k_laytrop_var_1221, p_sfluxzen_var_1222, p_taug_var_1223, p_taur_var_1224, prmu0_var_1225)
  USE yoesrtwn, ONLY: nspa_var_333, nspb_var_334
  USE yoesrta27, ONLY: absa_var_316, absb_var_317, layreffr_var_315, raylc_var_319, scalekur, sfluxrefc_var_318
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1209, kfdia_var_1210
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1211
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1212(kidia_var_1209 : kfdia_var_1210, klev_var_1211)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1213(kidia_var_1209 : kfdia_var_1210, klev_var_1211)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1214(kidia_var_1209 : kfdia_var_1210, klev_var_1211)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1215(kidia_var_1209 : kfdia_var_1210, klev_var_1211)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1216(kidia_var_1209 : kfdia_var_1210, klev_var_1211)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1217(kidia_var_1209 : kfdia_var_1210, klev_var_1211)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1218(kidia_var_1209 : kfdia_var_1210, klev_var_1211)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1219(kidia_var_1209 : kfdia_var_1210, klev_var_1211)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1220(kidia_var_1209 : kfdia_var_1210, klev_var_1211)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1221(kidia_var_1209 : kfdia_var_1210)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1222(kidia_var_1209 : kfdia_var_1210, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1223(kidia_var_1209 : kfdia_var_1210, klev_var_1211, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1224(kidia_var_1209 : kfdia_var_1210, klev_var_1211, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1225(kidia_var_1209 : kfdia_var_1210)
  INTEGER(KIND = 4) :: ig_var_1226, ind0_var_1227, ind1_var_1228, i_lay_var_1229, i_laysolfr_var_1230(kidia_var_1209 : kfdia_var_1210), i_nlayers_var_1231, iplon_var_1232
  INTEGER(KIND = 4) :: laytrop_min_var_1233, laytrop_max_var_1234
  REAL(KIND = 8) :: z_tauray_var_1235
  laytrop_min_var_1233 = MINVAL(k_laytrop_var_1221(kidia_var_1209 : kfdia_var_1210))
  laytrop_max_var_1234 = MAXVAL(k_laytrop_var_1221(kidia_var_1209 : kfdia_var_1210))
  i_nlayers_var_1231 = klev_var_1211
  DO iplon_var_1232 = kidia_var_1209, kfdia_var_1210
    i_laysolfr_var_1230(iplon_var_1232) = i_nlayers_var_1231
  END DO
  DO i_lay_var_1229 = 1, laytrop_min_var_1233
    DO iplon_var_1232 = kidia_var_1209, kfdia_var_1210
      ind0_var_1227 = ((k_jp_var_1216(iplon_var_1232, i_lay_var_1229) - 1) * 5 + (k_jt_var_1217(iplon_var_1232, i_lay_var_1229) - 1)) * nspa_var_333(27) + 1
      ind1_var_1228 = (k_jp_var_1216(iplon_var_1232, i_lay_var_1229) * 5 + (k_jt1_var_1218(iplon_var_1232, i_lay_var_1229) - 1)) * nspa_var_333(27) + 1
      DO ig_var_1226 = 1, 8
        z_tauray_var_1235 = p_colmol_var_1219(iplon_var_1232, i_lay_var_1229) * raylc_var_319(ig_var_1226)
        p_taug_var_1223(iplon_var_1232, i_lay_var_1229, ig_var_1226) = p_colo3_var_1220(iplon_var_1232, i_lay_var_1229) * (p_fac00_var_1212(iplon_var_1232, i_lay_var_1229) * absa_var_316(ind0_var_1227, ig_var_1226) + p_fac10_var_1214(iplon_var_1232, i_lay_var_1229) * absa_var_316(ind0_var_1227 + 1, ig_var_1226) + p_fac01_var_1213(iplon_var_1232, i_lay_var_1229) * absa_var_316(ind1_var_1228, ig_var_1226) + p_fac11_var_1215(iplon_var_1232, i_lay_var_1229) * absa_var_316(ind1_var_1228 + 1, ig_var_1226))
        p_taur_var_1224(iplon_var_1232, i_lay_var_1229, ig_var_1226) = z_tauray_var_1235
      END DO
    END DO
  END DO
  DO i_lay_var_1229 = laytrop_min_var_1233 + 1, laytrop_max_var_1234
    DO iplon_var_1232 = kidia_var_1209, kfdia_var_1210
      IF (i_lay_var_1229 <= k_laytrop_var_1221(iplon_var_1232)) THEN
        ind0_var_1227 = ((k_jp_var_1216(iplon_var_1232, i_lay_var_1229) - 1) * 5 + (k_jt_var_1217(iplon_var_1232, i_lay_var_1229) - 1)) * nspa_var_333(27) + 1
        ind1_var_1228 = (k_jp_var_1216(iplon_var_1232, i_lay_var_1229) * 5 + (k_jt1_var_1218(iplon_var_1232, i_lay_var_1229) - 1)) * nspa_var_333(27) + 1
        DO ig_var_1226 = 1, 8
          z_tauray_var_1235 = p_colmol_var_1219(iplon_var_1232, i_lay_var_1229) * raylc_var_319(ig_var_1226)
          p_taug_var_1223(iplon_var_1232, i_lay_var_1229, ig_var_1226) = p_colo3_var_1220(iplon_var_1232, i_lay_var_1229) * (p_fac00_var_1212(iplon_var_1232, i_lay_var_1229) * absa_var_316(ind0_var_1227, ig_var_1226) + p_fac10_var_1214(iplon_var_1232, i_lay_var_1229) * absa_var_316(ind0_var_1227 + 1, ig_var_1226) + p_fac01_var_1213(iplon_var_1232, i_lay_var_1229) * absa_var_316(ind1_var_1228, ig_var_1226) + p_fac11_var_1215(iplon_var_1232, i_lay_var_1229) * absa_var_316(ind1_var_1228 + 1, ig_var_1226))
          p_taur_var_1224(iplon_var_1232, i_lay_var_1229, ig_var_1226) = z_tauray_var_1235
        END DO
      ELSE
        IF (k_jp_var_1216(iplon_var_1232, i_lay_var_1229 - 1) < layreffr_var_315 .AND. k_jp_var_1216(iplon_var_1232, i_lay_var_1229) >= layreffr_var_315) i_laysolfr_var_1230(iplon_var_1232) = i_lay_var_1229
        ind0_var_1227 = ((k_jp_var_1216(iplon_var_1232, i_lay_var_1229) - 13) * 5 + (k_jt_var_1217(iplon_var_1232, i_lay_var_1229) - 1)) * nspb_var_334(27) + 1
        ind1_var_1228 = ((k_jp_var_1216(iplon_var_1232, i_lay_var_1229) - 12) * 5 + (k_jt1_var_1218(iplon_var_1232, i_lay_var_1229) - 1)) * nspb_var_334(27) + 1
        DO ig_var_1226 = 1, 8
          z_tauray_var_1235 = p_colmol_var_1219(iplon_var_1232, i_lay_var_1229) * raylc_var_319(ig_var_1226)
          p_taug_var_1223(iplon_var_1232, i_lay_var_1229, ig_var_1226) = p_colo3_var_1220(iplon_var_1232, i_lay_var_1229) * (p_fac00_var_1212(iplon_var_1232, i_lay_var_1229) * absb_var_317(ind0_var_1227, ig_var_1226) + p_fac10_var_1214(iplon_var_1232, i_lay_var_1229) * absb_var_317(ind0_var_1227 + 1, ig_var_1226) + p_fac01_var_1213(iplon_var_1232, i_lay_var_1229) * absb_var_317(ind1_var_1228, ig_var_1226) + p_fac11_var_1215(iplon_var_1232, i_lay_var_1229) * absb_var_317(ind1_var_1228 + 1, ig_var_1226))
          IF (i_lay_var_1229 == i_laysolfr_var_1230(iplon_var_1232)) p_sfluxzen_var_1222(iplon_var_1232, ig_var_1226) = scalekur * sfluxrefc_var_318(ig_var_1226)
          p_taur_var_1224(iplon_var_1232, i_lay_var_1229, ig_var_1226) = z_tauray_var_1235
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1229 = laytrop_max_var_1234 + 1, i_nlayers_var_1231
    DO iplon_var_1232 = kidia_var_1209, kfdia_var_1210
      IF (k_jp_var_1216(iplon_var_1232, i_lay_var_1229 - 1) < layreffr_var_315 .AND. k_jp_var_1216(iplon_var_1232, i_lay_var_1229) >= layreffr_var_315) i_laysolfr_var_1230(iplon_var_1232) = i_lay_var_1229
      ind0_var_1227 = ((k_jp_var_1216(iplon_var_1232, i_lay_var_1229) - 13) * 5 + (k_jt_var_1217(iplon_var_1232, i_lay_var_1229) - 1)) * nspb_var_334(27) + 1
      ind1_var_1228 = ((k_jp_var_1216(iplon_var_1232, i_lay_var_1229) - 12) * 5 + (k_jt1_var_1218(iplon_var_1232, i_lay_var_1229) - 1)) * nspb_var_334(27) + 1
      DO ig_var_1226 = 1, 8
        z_tauray_var_1235 = p_colmol_var_1219(iplon_var_1232, i_lay_var_1229) * raylc_var_319(ig_var_1226)
        p_taug_var_1223(iplon_var_1232, i_lay_var_1229, ig_var_1226) = p_colo3_var_1220(iplon_var_1232, i_lay_var_1229) * (p_fac00_var_1212(iplon_var_1232, i_lay_var_1229) * absb_var_317(ind0_var_1227, ig_var_1226) + p_fac10_var_1214(iplon_var_1232, i_lay_var_1229) * absb_var_317(ind0_var_1227 + 1, ig_var_1226) + p_fac01_var_1213(iplon_var_1232, i_lay_var_1229) * absb_var_317(ind1_var_1228, ig_var_1226) + p_fac11_var_1215(iplon_var_1232, i_lay_var_1229) * absb_var_317(ind1_var_1228 + 1, ig_var_1226))
        IF (i_lay_var_1229 == i_laysolfr_var_1230(iplon_var_1232)) p_sfluxzen_var_1222(iplon_var_1232, ig_var_1226) = scalekur * sfluxrefc_var_318(ig_var_1226)
        p_taur_var_1224(iplon_var_1232, i_lay_var_1229, ig_var_1226) = z_tauray_var_1235
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol27
SUBROUTINE srtm_taumol26(kidia_var_1236, kfdia_var_1237, klev_var_1238, p_colmol_var_1239, k_laytrop_var_1240, p_sfluxzen_var_1241, p_taug_var_1242, p_taur_var_1243, prmu0_var_1244)
  USE yoesrta26, ONLY: raylc_var_314, sfluxrefc_var_313
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1236, kfdia_var_1237
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1238
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1239(kidia_var_1236 : kfdia_var_1237, klev_var_1238)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1240(kidia_var_1236 : kfdia_var_1237)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1241(kidia_var_1236 : kfdia_var_1237, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1242(kidia_var_1236 : kfdia_var_1237, klev_var_1238, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1243(kidia_var_1236 : kfdia_var_1237, klev_var_1238, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1244(kidia_var_1236 : kfdia_var_1237)
  INTEGER(KIND = 4) :: ig_var_1245, i_lay_var_1246, i_laysolfr_var_1247(kidia_var_1236 : kfdia_var_1237), i_nlayers_var_1248, iplon_var_1249
  INTEGER(KIND = 4) :: laytrop_min_var_1250, laytrop_max_var_1251
  laytrop_min_var_1250 = MINVAL(k_laytrop_var_1240(kidia_var_1236 : kfdia_var_1237))
  laytrop_max_var_1251 = MAXVAL(k_laytrop_var_1240(kidia_var_1236 : kfdia_var_1237))
  i_nlayers_var_1248 = klev_var_1238
  DO iplon_var_1249 = kidia_var_1236, kfdia_var_1237
    i_laysolfr_var_1247(iplon_var_1249) = k_laytrop_var_1240(iplon_var_1249)
  END DO
  DO i_lay_var_1246 = 1, laytrop_min_var_1250
    DO iplon_var_1249 = kidia_var_1236, kfdia_var_1237
      DO ig_var_1245 = 1, 6
        IF (i_lay_var_1246 == i_laysolfr_var_1247(iplon_var_1249)) p_sfluxzen_var_1241(iplon_var_1249, ig_var_1245) = sfluxrefc_var_313(ig_var_1245)
        p_taug_var_1242(iplon_var_1249, i_lay_var_1246, ig_var_1245) = 0.0D0
        p_taur_var_1243(iplon_var_1249, i_lay_var_1246, ig_var_1245) = p_colmol_var_1239(iplon_var_1249, i_lay_var_1246) * raylc_var_314(ig_var_1245)
      END DO
    END DO
  END DO
  DO i_lay_var_1246 = laytrop_min_var_1250 + 1, laytrop_max_var_1251
    DO iplon_var_1249 = kidia_var_1236, kfdia_var_1237
      IF (i_lay_var_1246 <= k_laytrop_var_1240(iplon_var_1249)) THEN
        DO ig_var_1245 = 1, 6
          IF (i_lay_var_1246 == i_laysolfr_var_1247(iplon_var_1249)) p_sfluxzen_var_1241(iplon_var_1249, ig_var_1245) = sfluxrefc_var_313(ig_var_1245)
          p_taug_var_1242(iplon_var_1249, i_lay_var_1246, ig_var_1245) = 0.0D0
          p_taur_var_1243(iplon_var_1249, i_lay_var_1246, ig_var_1245) = p_colmol_var_1239(iplon_var_1249, i_lay_var_1246) * raylc_var_314(ig_var_1245)
        END DO
      ELSE
        DO ig_var_1245 = 1, 6
          p_taug_var_1242(iplon_var_1249, i_lay_var_1246, ig_var_1245) = 0.0D0
          p_taur_var_1243(iplon_var_1249, i_lay_var_1246, ig_var_1245) = p_colmol_var_1239(iplon_var_1249, i_lay_var_1246) * raylc_var_314(ig_var_1245)
        END DO
      END IF
    END DO
  END DO
  DO ig_var_1245 = 1, 6
    DO i_lay_var_1246 = laytrop_max_var_1251 + 1, i_nlayers_var_1248
      DO iplon_var_1249 = kidia_var_1236, kfdia_var_1237
        p_taug_var_1242(iplon_var_1249, i_lay_var_1246, ig_var_1245) = 0.0D0
        p_taur_var_1243(iplon_var_1249, i_lay_var_1246, ig_var_1245) = p_colmol_var_1239(iplon_var_1249, i_lay_var_1246) * raylc_var_314(ig_var_1245)
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol26
SUBROUTINE srtm_taumol18(kidia_var_1252, kfdia_var_1253, klev_var_1254, p_fac00_var_1255, p_fac01_var_1256, p_fac10_var_1257, p_fac11_var_1258, k_jp_var_1259, k_jt_var_1260, k_jt1_var_1261, p_oneminus_var_1262, p_colh2o_var_1263, p_colch4_var_1264, p_colmol_var_1265, k_laytrop_var_1266, p_selffac_var_1267, p_selffrac_var_1268, k_indself_var_1269, p_forfac_var_1270, p_forfrac_var_1271, k_indfor_var_1272, p_sfluxzen_var_1273, p_taug_var_1274, p_taur_var_1275, prmu0_var_1276)
  USE yoesrta18, ONLY: absa_var_256, absb_var_257, forrefc_var_259, layreffr_var_255, rayl_var_253, selfrefc_var_258, sfluxrefc_var_260, strrat_var_254
  USE yoesrtwn, ONLY: nspa_var_333, nspb_var_334
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1252, kfdia_var_1253
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1254
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1255(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1256(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1257(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1258(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1259(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1260(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1261(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1262(kidia_var_1252 : kfdia_var_1253)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1263(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_1264(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1265(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1266(kidia_var_1252 : kfdia_var_1253)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1267(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1268(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1269(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1270(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1271(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1272(kidia_var_1252 : kfdia_var_1253, klev_var_1254)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1273(kidia_var_1252 : kfdia_var_1253, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1274(kidia_var_1252 : kfdia_var_1253, klev_var_1254, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1275(kidia_var_1252 : kfdia_var_1253, klev_var_1254, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1276(kidia_var_1252 : kfdia_var_1253)
  INTEGER(KIND = 4) :: ig_var_1277, ind0_var_1278, ind1_var_1279, inds_var_1280, indf_var_1281, js_var_1282, i_lay_var_1283, i_laysolfr_var_1284(kidia_var_1252 : kfdia_var_1253), i_nlayers_var_1285, iplon_var_1286
  INTEGER(KIND = 4) :: laytrop_min_var_1287, laytrop_max_var_1288
  REAL(KIND = 8) :: z_fs_var_1289, z_speccomb_var_1290, z_specmult_var_1291, z_specparm_var_1292, z_tauray_var_1293
  laytrop_min_var_1287 = MINVAL(k_laytrop_var_1266(kidia_var_1252 : kfdia_var_1253))
  laytrop_max_var_1288 = MAXVAL(k_laytrop_var_1266(kidia_var_1252 : kfdia_var_1253))
  i_nlayers_var_1285 = klev_var_1254
  DO iplon_var_1286 = kidia_var_1252, kfdia_var_1253
    i_laysolfr_var_1284(iplon_var_1286) = k_laytrop_var_1266(iplon_var_1286)
  END DO
  DO i_lay_var_1283 = 1, laytrop_min_var_1287
    DO iplon_var_1286 = kidia_var_1252, kfdia_var_1253
      IF (k_jp_var_1259(iplon_var_1286, i_lay_var_1283) < layreffr_var_255 .AND. k_jp_var_1259(iplon_var_1286, i_lay_var_1283 + 1) >= layreffr_var_255) i_laysolfr_var_1284(iplon_var_1286) = MIN(i_lay_var_1283 + 1, k_laytrop_var_1266(iplon_var_1286))
      z_speccomb_var_1290 = p_colh2o_var_1263(iplon_var_1286, i_lay_var_1283) + strrat_var_254 * p_colch4_var_1264(iplon_var_1286, i_lay_var_1283)
      z_specparm_var_1292 = p_colh2o_var_1263(iplon_var_1286, i_lay_var_1283) / z_speccomb_var_1290
      z_specparm_var_1292 = MIN(p_oneminus_var_1262(iplon_var_1286), z_specparm_var_1292)
      z_specmult_var_1291 = 8.0D0 * (z_specparm_var_1292)
      js_var_1282 = 1 + INT(z_specmult_var_1291)
      z_fs_var_1289 = z_specmult_var_1291 - AINT(z_specmult_var_1291)
      ind0_var_1278 = ((k_jp_var_1259(iplon_var_1286, i_lay_var_1283) - 1) * 5 + (k_jt_var_1260(iplon_var_1286, i_lay_var_1283) - 1)) * nspa_var_333(18) + js_var_1282
      ind1_var_1279 = (k_jp_var_1259(iplon_var_1286, i_lay_var_1283) * 5 + (k_jt1_var_1261(iplon_var_1286, i_lay_var_1283) - 1)) * nspa_var_333(18) + js_var_1282
      inds_var_1280 = k_indself_var_1269(iplon_var_1286, i_lay_var_1283)
      indf_var_1281 = k_indfor_var_1272(iplon_var_1286, i_lay_var_1283)
      z_tauray_var_1293 = p_colmol_var_1265(iplon_var_1286, i_lay_var_1283) * rayl_var_253
      DO ig_var_1277 = 1, 8
        p_taug_var_1274(iplon_var_1286, i_lay_var_1283, ig_var_1277) = z_speccomb_var_1290 * ((1.0D0 - z_fs_var_1289) * (absa_var_256(ind0_var_1278, ig_var_1277) * p_fac00_var_1255(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind0_var_1278 + 9, ig_var_1277) * p_fac10_var_1257(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind1_var_1279, ig_var_1277) * p_fac01_var_1256(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind1_var_1279 + 9, ig_var_1277) * p_fac11_var_1258(iplon_var_1286, i_lay_var_1283)) + z_fs_var_1289 * (absa_var_256(ind0_var_1278 + 1, ig_var_1277) * p_fac00_var_1255(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind0_var_1278 + 10, ig_var_1277) * p_fac10_var_1257(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind1_var_1279 + 1, ig_var_1277) * p_fac01_var_1256(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind1_var_1279 + 10, ig_var_1277) * p_fac11_var_1258(iplon_var_1286, i_lay_var_1283))) + p_colh2o_var_1263(iplon_var_1286, i_lay_var_1283) * (p_selffac_var_1267(iplon_var_1286, i_lay_var_1283) * (selfrefc_var_258(inds_var_1280, ig_var_1277) + p_selffrac_var_1268(iplon_var_1286, i_lay_var_1283) * (selfrefc_var_258(inds_var_1280 + 1, ig_var_1277) - selfrefc_var_258(inds_var_1280, ig_var_1277))) + p_forfac_var_1270(iplon_var_1286, i_lay_var_1283) * (forrefc_var_259(indf_var_1281, ig_var_1277) + p_forfrac_var_1271(iplon_var_1286, i_lay_var_1283) * (forrefc_var_259(indf_var_1281 + 1, ig_var_1277) - forrefc_var_259(indf_var_1281, ig_var_1277))))
        IF (i_lay_var_1283 == i_laysolfr_var_1284(iplon_var_1286)) p_sfluxzen_var_1273(iplon_var_1286, ig_var_1277) = sfluxrefc_var_260(ig_var_1277, js_var_1282) + z_fs_var_1289 * (sfluxrefc_var_260(ig_var_1277, js_var_1282 + 1) - sfluxrefc_var_260(ig_var_1277, js_var_1282))
        p_taur_var_1275(iplon_var_1286, i_lay_var_1283, ig_var_1277) = z_tauray_var_1293
      END DO
    END DO
  END DO
  DO i_lay_var_1283 = laytrop_min_var_1287 + 1, laytrop_max_var_1288
    DO iplon_var_1286 = kidia_var_1252, kfdia_var_1253
      IF (i_lay_var_1283 <= k_laytrop_var_1266(iplon_var_1286)) THEN
        IF (k_jp_var_1259(iplon_var_1286, i_lay_var_1283) < layreffr_var_255 .AND. k_jp_var_1259(iplon_var_1286, i_lay_var_1283 + 1) >= layreffr_var_255) i_laysolfr_var_1284(iplon_var_1286) = MIN(i_lay_var_1283 + 1, k_laytrop_var_1266(iplon_var_1286))
        z_speccomb_var_1290 = p_colh2o_var_1263(iplon_var_1286, i_lay_var_1283) + strrat_var_254 * p_colch4_var_1264(iplon_var_1286, i_lay_var_1283)
        z_specparm_var_1292 = p_colh2o_var_1263(iplon_var_1286, i_lay_var_1283) / z_speccomb_var_1290
        z_specparm_var_1292 = MIN(p_oneminus_var_1262(iplon_var_1286), z_specparm_var_1292)
        z_specmult_var_1291 = 8.0D0 * (z_specparm_var_1292)
        js_var_1282 = 1 + INT(z_specmult_var_1291)
        z_fs_var_1289 = z_specmult_var_1291 - AINT(z_specmult_var_1291)
        ind0_var_1278 = ((k_jp_var_1259(iplon_var_1286, i_lay_var_1283) - 1) * 5 + (k_jt_var_1260(iplon_var_1286, i_lay_var_1283) - 1)) * nspa_var_333(18) + js_var_1282
        ind1_var_1279 = (k_jp_var_1259(iplon_var_1286, i_lay_var_1283) * 5 + (k_jt1_var_1261(iplon_var_1286, i_lay_var_1283) - 1)) * nspa_var_333(18) + js_var_1282
        inds_var_1280 = k_indself_var_1269(iplon_var_1286, i_lay_var_1283)
        indf_var_1281 = k_indfor_var_1272(iplon_var_1286, i_lay_var_1283)
        z_tauray_var_1293 = p_colmol_var_1265(iplon_var_1286, i_lay_var_1283) * rayl_var_253
        DO ig_var_1277 = 1, 8
          p_taug_var_1274(iplon_var_1286, i_lay_var_1283, ig_var_1277) = z_speccomb_var_1290 * ((1.0D0 - z_fs_var_1289) * (absa_var_256(ind0_var_1278, ig_var_1277) * p_fac00_var_1255(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind0_var_1278 + 9, ig_var_1277) * p_fac10_var_1257(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind1_var_1279, ig_var_1277) * p_fac01_var_1256(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind1_var_1279 + 9, ig_var_1277) * p_fac11_var_1258(iplon_var_1286, i_lay_var_1283)) + z_fs_var_1289 * (absa_var_256(ind0_var_1278 + 1, ig_var_1277) * p_fac00_var_1255(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind0_var_1278 + 10, ig_var_1277) * p_fac10_var_1257(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind1_var_1279 + 1, ig_var_1277) * p_fac01_var_1256(iplon_var_1286, i_lay_var_1283) + absa_var_256(ind1_var_1279 + 10, ig_var_1277) * p_fac11_var_1258(iplon_var_1286, i_lay_var_1283))) + p_colh2o_var_1263(iplon_var_1286, i_lay_var_1283) * (p_selffac_var_1267(iplon_var_1286, i_lay_var_1283) * (selfrefc_var_258(inds_var_1280, ig_var_1277) + p_selffrac_var_1268(iplon_var_1286, i_lay_var_1283) * (selfrefc_var_258(inds_var_1280 + 1, ig_var_1277) - selfrefc_var_258(inds_var_1280, ig_var_1277))) + p_forfac_var_1270(iplon_var_1286, i_lay_var_1283) * (forrefc_var_259(indf_var_1281, ig_var_1277) + p_forfrac_var_1271(iplon_var_1286, i_lay_var_1283) * (forrefc_var_259(indf_var_1281 + 1, ig_var_1277) - forrefc_var_259(indf_var_1281, ig_var_1277))))
          IF (i_lay_var_1283 == i_laysolfr_var_1284(iplon_var_1286)) p_sfluxzen_var_1273(iplon_var_1286, ig_var_1277) = sfluxrefc_var_260(ig_var_1277, js_var_1282) + z_fs_var_1289 * (sfluxrefc_var_260(ig_var_1277, js_var_1282 + 1) - sfluxrefc_var_260(ig_var_1277, js_var_1282))
          p_taur_var_1275(iplon_var_1286, i_lay_var_1283, ig_var_1277) = z_tauray_var_1293
        END DO
      ELSE
        ind0_var_1278 = ((k_jp_var_1259(iplon_var_1286, i_lay_var_1283) - 13) * 5 + (k_jt_var_1260(iplon_var_1286, i_lay_var_1283) - 1)) * nspb_var_334(18) + 1
        ind1_var_1279 = ((k_jp_var_1259(iplon_var_1286, i_lay_var_1283) - 12) * 5 + (k_jt1_var_1261(iplon_var_1286, i_lay_var_1283) - 1)) * nspb_var_334(18) + 1
        z_tauray_var_1293 = p_colmol_var_1265(iplon_var_1286, i_lay_var_1283) * rayl_var_253
        DO ig_var_1277 = 1, 8
          p_taug_var_1274(iplon_var_1286, i_lay_var_1283, ig_var_1277) = p_colch4_var_1264(iplon_var_1286, i_lay_var_1283) * (p_fac00_var_1255(iplon_var_1286, i_lay_var_1283) * absb_var_257(ind0_var_1278, ig_var_1277) + p_fac10_var_1257(iplon_var_1286, i_lay_var_1283) * absb_var_257(ind0_var_1278 + 1, ig_var_1277) + p_fac01_var_1256(iplon_var_1286, i_lay_var_1283) * absb_var_257(ind1_var_1279, ig_var_1277) + p_fac11_var_1258(iplon_var_1286, i_lay_var_1283) * absb_var_257(ind1_var_1279 + 1, ig_var_1277))
          p_taur_var_1275(iplon_var_1286, i_lay_var_1283, ig_var_1277) = z_tauray_var_1293
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1283 = laytrop_max_var_1288 + 1, i_nlayers_var_1285
    DO iplon_var_1286 = kidia_var_1252, kfdia_var_1253
      ind0_var_1278 = ((k_jp_var_1259(iplon_var_1286, i_lay_var_1283) - 13) * 5 + (k_jt_var_1260(iplon_var_1286, i_lay_var_1283) - 1)) * nspb_var_334(18) + 1
      ind1_var_1279 = ((k_jp_var_1259(iplon_var_1286, i_lay_var_1283) - 12) * 5 + (k_jt1_var_1261(iplon_var_1286, i_lay_var_1283) - 1)) * nspb_var_334(18) + 1
      z_tauray_var_1293 = p_colmol_var_1265(iplon_var_1286, i_lay_var_1283) * rayl_var_253
      DO ig_var_1277 = 1, 8
        p_taug_var_1274(iplon_var_1286, i_lay_var_1283, ig_var_1277) = p_colch4_var_1264(iplon_var_1286, i_lay_var_1283) * (p_fac00_var_1255(iplon_var_1286, i_lay_var_1283) * absb_var_257(ind0_var_1278, ig_var_1277) + p_fac10_var_1257(iplon_var_1286, i_lay_var_1283) * absb_var_257(ind0_var_1278 + 1, ig_var_1277) + p_fac01_var_1256(iplon_var_1286, i_lay_var_1283) * absb_var_257(ind1_var_1279, ig_var_1277) + p_fac11_var_1258(iplon_var_1286, i_lay_var_1283) * absb_var_257(ind1_var_1279 + 1, ig_var_1277))
        p_taur_var_1275(iplon_var_1286, i_lay_var_1283, ig_var_1277) = z_tauray_var_1293
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol18
SUBROUTINE srtm_taumol24(kidia_var_1294, kfdia_var_1295, klev_var_1296, p_fac00_var_1297, p_fac01_var_1298, p_fac10_var_1299, p_fac11_var_1300, k_jp_var_1301, k_jt_var_1302, k_jt1_var_1303, p_oneminus_var_1304, p_colh2o_var_1305, p_colmol_var_1306, p_colo2_var_1307, p_colo3_var_1308, k_laytrop_var_1309, p_selffac_var_1310, p_selffrac_var_1311, k_indself_var_1312, p_forfac_var_1313, p_forfrac_var_1314, k_indfor_var_1315, p_sfluxzen_var_1316, p_taug_var_1317, p_taur_var_1318, prmu0_var_1319)
  USE yoesrta24, ONLY: absa_var_300, absb_var_301, abso3ac_var_305, abso3bc_var_306, forrefc_var_303, layreffr_var_299, raylac, raylbc, selfrefc_var_302, sfluxrefc_var_304, strrat_var_298
  USE yoesrtwn, ONLY: nspa_var_333, nspb_var_334
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1294, kfdia_var_1295
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1296
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1297(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1298(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1299(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1300(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1301(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1302(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1303(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1304(kidia_var_1294 : kfdia_var_1295)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1305(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1306(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_1307(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1308(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1309(kidia_var_1294 : kfdia_var_1295)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1310(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1311(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1312(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1313(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1314(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1315(kidia_var_1294 : kfdia_var_1295, klev_var_1296)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1316(kidia_var_1294 : kfdia_var_1295, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1317(kidia_var_1294 : kfdia_var_1295, klev_var_1296, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1318(kidia_var_1294 : kfdia_var_1295, klev_var_1296, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1319(kidia_var_1294 : kfdia_var_1295)
  INTEGER(KIND = 4) :: ig_var_1320, ind0_var_1321, ind1_var_1322, inds_var_1323, indf_var_1324, js_var_1325, i_lay_var_1326, i_laysolfr_var_1327(kidia_var_1294 : kfdia_var_1295), i_nlayers_var_1328, iplon_var_1329
  INTEGER(KIND = 4) :: laytrop_min_var_1330, laytrop_max_var_1331
  REAL(KIND = 8) :: z_fs_var_1332, z_speccomb_var_1333, z_specmult_var_1334, z_specparm_var_1335, z_tauray_var_1336
  laytrop_min_var_1330 = MINVAL(k_laytrop_var_1309(kidia_var_1294 : kfdia_var_1295))
  laytrop_max_var_1331 = MAXVAL(k_laytrop_var_1309(kidia_var_1294 : kfdia_var_1295))
  i_nlayers_var_1328 = klev_var_1296
  DO iplon_var_1329 = kidia_var_1294, kfdia_var_1295
    i_laysolfr_var_1327(iplon_var_1329) = k_laytrop_var_1309(iplon_var_1329)
  END DO
  DO i_lay_var_1326 = 1, laytrop_min_var_1330
    DO iplon_var_1329 = kidia_var_1294, kfdia_var_1295
      IF (k_jp_var_1301(iplon_var_1329, i_lay_var_1326) < layreffr_var_299 .AND. k_jp_var_1301(iplon_var_1329, i_lay_var_1326 + 1) >= layreffr_var_299) i_laysolfr_var_1327(iplon_var_1329) = MIN(i_lay_var_1326 + 1, k_laytrop_var_1309(iplon_var_1329))
      z_speccomb_var_1333 = p_colh2o_var_1305(iplon_var_1329, i_lay_var_1326) + strrat_var_298 * p_colo2_var_1307(iplon_var_1329, i_lay_var_1326)
      z_specparm_var_1335 = p_colh2o_var_1305(iplon_var_1329, i_lay_var_1326) / z_speccomb_var_1333
      z_specparm_var_1335 = MIN(p_oneminus_var_1304(iplon_var_1329), z_specparm_var_1335)
      z_specmult_var_1334 = 8.0D0 * (z_specparm_var_1335)
      js_var_1325 = 1 + INT(z_specmult_var_1334)
      z_fs_var_1332 = z_specmult_var_1334 - AINT(z_specmult_var_1334)
      ind0_var_1321 = ((k_jp_var_1301(iplon_var_1329, i_lay_var_1326) - 1) * 5 + (k_jt_var_1302(iplon_var_1329, i_lay_var_1326) - 1)) * nspa_var_333(24) + js_var_1325
      ind1_var_1322 = (k_jp_var_1301(iplon_var_1329, i_lay_var_1326) * 5 + (k_jt1_var_1303(iplon_var_1329, i_lay_var_1326) - 1)) * nspa_var_333(24) + js_var_1325
      inds_var_1323 = k_indself_var_1312(iplon_var_1329, i_lay_var_1326)
      indf_var_1324 = k_indfor_var_1315(iplon_var_1329, i_lay_var_1326)
      DO ig_var_1320 = 1, 8
        z_tauray_var_1336 = p_colmol_var_1306(iplon_var_1329, i_lay_var_1326) * (raylac(ig_var_1320, js_var_1325) + z_fs_var_1332 * (raylac(ig_var_1320, js_var_1325 + 1) - raylac(ig_var_1320, js_var_1325)))
        p_taug_var_1317(iplon_var_1329, i_lay_var_1326, ig_var_1320) = z_speccomb_var_1333 * ((1.0D0 - z_fs_var_1332) * (absa_var_300(ind0_var_1321, ig_var_1320) * p_fac00_var_1297(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind0_var_1321 + 9, ig_var_1320) * p_fac10_var_1299(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind1_var_1322, ig_var_1320) * p_fac01_var_1298(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind1_var_1322 + 9, ig_var_1320) * p_fac11_var_1300(iplon_var_1329, i_lay_var_1326)) + z_fs_var_1332 * (absa_var_300(ind0_var_1321 + 1, ig_var_1320) * p_fac00_var_1297(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind0_var_1321 + 10, ig_var_1320) * p_fac10_var_1299(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind1_var_1322 + 1, ig_var_1320) * p_fac01_var_1298(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind1_var_1322 + 10, ig_var_1320) * p_fac11_var_1300(iplon_var_1329, i_lay_var_1326))) + p_colo3_var_1308(iplon_var_1329, i_lay_var_1326) * abso3ac_var_305(ig_var_1320) + p_colh2o_var_1305(iplon_var_1329, i_lay_var_1326) * (p_selffac_var_1310(iplon_var_1329, i_lay_var_1326) * (selfrefc_var_302(inds_var_1323, ig_var_1320) + p_selffrac_var_1311(iplon_var_1329, i_lay_var_1326) * (selfrefc_var_302(inds_var_1323 + 1, ig_var_1320) - selfrefc_var_302(inds_var_1323, ig_var_1320))) + p_forfac_var_1313(iplon_var_1329, i_lay_var_1326) * (forrefc_var_303(indf_var_1324, ig_var_1320) + p_forfrac_var_1314(iplon_var_1329, i_lay_var_1326) * (forrefc_var_303(indf_var_1324 + 1, ig_var_1320) - forrefc_var_303(indf_var_1324, ig_var_1320))))
        IF (i_lay_var_1326 == i_laysolfr_var_1327(iplon_var_1329)) p_sfluxzen_var_1316(iplon_var_1329, ig_var_1320) = sfluxrefc_var_304(ig_var_1320, js_var_1325) + z_fs_var_1332 * (sfluxrefc_var_304(ig_var_1320, js_var_1325 + 1) - sfluxrefc_var_304(ig_var_1320, js_var_1325))
        p_taur_var_1318(iplon_var_1329, i_lay_var_1326, ig_var_1320) = z_tauray_var_1336
      END DO
    END DO
  END DO
  DO i_lay_var_1326 = laytrop_min_var_1330 + 1, laytrop_max_var_1331
    DO iplon_var_1329 = kidia_var_1294, kfdia_var_1295
      IF (i_lay_var_1326 <= k_laytrop_var_1309(iplon_var_1329)) THEN
        IF (k_jp_var_1301(iplon_var_1329, i_lay_var_1326) < layreffr_var_299 .AND. k_jp_var_1301(iplon_var_1329, i_lay_var_1326 + 1) >= layreffr_var_299) i_laysolfr_var_1327(iplon_var_1329) = MIN(i_lay_var_1326 + 1, k_laytrop_var_1309(iplon_var_1329))
        z_speccomb_var_1333 = p_colh2o_var_1305(iplon_var_1329, i_lay_var_1326) + strrat_var_298 * p_colo2_var_1307(iplon_var_1329, i_lay_var_1326)
        z_specparm_var_1335 = p_colh2o_var_1305(iplon_var_1329, i_lay_var_1326) / z_speccomb_var_1333
        z_specparm_var_1335 = MIN(p_oneminus_var_1304(iplon_var_1329), z_specparm_var_1335)
        z_specmult_var_1334 = 8.0D0 * (z_specparm_var_1335)
        js_var_1325 = 1 + INT(z_specmult_var_1334)
        z_fs_var_1332 = z_specmult_var_1334 - AINT(z_specmult_var_1334)
        ind0_var_1321 = ((k_jp_var_1301(iplon_var_1329, i_lay_var_1326) - 1) * 5 + (k_jt_var_1302(iplon_var_1329, i_lay_var_1326) - 1)) * nspa_var_333(24) + js_var_1325
        ind1_var_1322 = (k_jp_var_1301(iplon_var_1329, i_lay_var_1326) * 5 + (k_jt1_var_1303(iplon_var_1329, i_lay_var_1326) - 1)) * nspa_var_333(24) + js_var_1325
        inds_var_1323 = k_indself_var_1312(iplon_var_1329, i_lay_var_1326)
        indf_var_1324 = k_indfor_var_1315(iplon_var_1329, i_lay_var_1326)
        DO ig_var_1320 = 1, 8
          z_tauray_var_1336 = p_colmol_var_1306(iplon_var_1329, i_lay_var_1326) * (raylac(ig_var_1320, js_var_1325) + z_fs_var_1332 * (raylac(ig_var_1320, js_var_1325 + 1) - raylac(ig_var_1320, js_var_1325)))
          p_taug_var_1317(iplon_var_1329, i_lay_var_1326, ig_var_1320) = z_speccomb_var_1333 * ((1.0D0 - z_fs_var_1332) * (absa_var_300(ind0_var_1321, ig_var_1320) * p_fac00_var_1297(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind0_var_1321 + 9, ig_var_1320) * p_fac10_var_1299(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind1_var_1322, ig_var_1320) * p_fac01_var_1298(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind1_var_1322 + 9, ig_var_1320) * p_fac11_var_1300(iplon_var_1329, i_lay_var_1326)) + z_fs_var_1332 * (absa_var_300(ind0_var_1321 + 1, ig_var_1320) * p_fac00_var_1297(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind0_var_1321 + 10, ig_var_1320) * p_fac10_var_1299(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind1_var_1322 + 1, ig_var_1320) * p_fac01_var_1298(iplon_var_1329, i_lay_var_1326) + absa_var_300(ind1_var_1322 + 10, ig_var_1320) * p_fac11_var_1300(iplon_var_1329, i_lay_var_1326))) + p_colo3_var_1308(iplon_var_1329, i_lay_var_1326) * abso3ac_var_305(ig_var_1320) + p_colh2o_var_1305(iplon_var_1329, i_lay_var_1326) * (p_selffac_var_1310(iplon_var_1329, i_lay_var_1326) * (selfrefc_var_302(inds_var_1323, ig_var_1320) + p_selffrac_var_1311(iplon_var_1329, i_lay_var_1326) * (selfrefc_var_302(inds_var_1323 + 1, ig_var_1320) - selfrefc_var_302(inds_var_1323, ig_var_1320))) + p_forfac_var_1313(iplon_var_1329, i_lay_var_1326) * (forrefc_var_303(indf_var_1324, ig_var_1320) + p_forfrac_var_1314(iplon_var_1329, i_lay_var_1326) * (forrefc_var_303(indf_var_1324 + 1, ig_var_1320) - forrefc_var_303(indf_var_1324, ig_var_1320))))
          IF (i_lay_var_1326 == i_laysolfr_var_1327(iplon_var_1329)) p_sfluxzen_var_1316(iplon_var_1329, ig_var_1320) = sfluxrefc_var_304(ig_var_1320, js_var_1325) + z_fs_var_1332 * (sfluxrefc_var_304(ig_var_1320, js_var_1325 + 1) - sfluxrefc_var_304(ig_var_1320, js_var_1325))
          p_taur_var_1318(iplon_var_1329, i_lay_var_1326, ig_var_1320) = z_tauray_var_1336
        END DO
      ELSE
        ind0_var_1321 = ((k_jp_var_1301(iplon_var_1329, i_lay_var_1326) - 13) * 5 + (k_jt_var_1302(iplon_var_1329, i_lay_var_1326) - 1)) * nspb_var_334(24) + 1
        ind1_var_1322 = ((k_jp_var_1301(iplon_var_1329, i_lay_var_1326) - 12) * 5 + (k_jt1_var_1303(iplon_var_1329, i_lay_var_1326) - 1)) * nspb_var_334(24) + 1
        DO ig_var_1320 = 1, 8
          z_tauray_var_1336 = p_colmol_var_1306(iplon_var_1329, i_lay_var_1326) * raylbc(ig_var_1320)
          p_taug_var_1317(iplon_var_1329, i_lay_var_1326, ig_var_1320) = p_colo2_var_1307(iplon_var_1329, i_lay_var_1326) * (p_fac00_var_1297(iplon_var_1329, i_lay_var_1326) * absb_var_301(ind0_var_1321, ig_var_1320) + p_fac10_var_1299(iplon_var_1329, i_lay_var_1326) * absb_var_301(ind0_var_1321 + 1, ig_var_1320) + p_fac01_var_1298(iplon_var_1329, i_lay_var_1326) * absb_var_301(ind1_var_1322, ig_var_1320) + p_fac11_var_1300(iplon_var_1329, i_lay_var_1326) * absb_var_301(ind1_var_1322 + 1, ig_var_1320)) + p_colo3_var_1308(iplon_var_1329, i_lay_var_1326) * abso3bc_var_306(ig_var_1320)
          p_taur_var_1318(iplon_var_1329, i_lay_var_1326, ig_var_1320) = z_tauray_var_1336
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1326 = laytrop_max_var_1331 + 1, i_nlayers_var_1328
    DO iplon_var_1329 = kidia_var_1294, kfdia_var_1295
      ind0_var_1321 = ((k_jp_var_1301(iplon_var_1329, i_lay_var_1326) - 13) * 5 + (k_jt_var_1302(iplon_var_1329, i_lay_var_1326) - 1)) * nspb_var_334(24) + 1
      ind1_var_1322 = ((k_jp_var_1301(iplon_var_1329, i_lay_var_1326) - 12) * 5 + (k_jt1_var_1303(iplon_var_1329, i_lay_var_1326) - 1)) * nspb_var_334(24) + 1
      DO ig_var_1320 = 1, 8
        z_tauray_var_1336 = p_colmol_var_1306(iplon_var_1329, i_lay_var_1326) * raylbc(ig_var_1320)
        p_taug_var_1317(iplon_var_1329, i_lay_var_1326, ig_var_1320) = p_colo2_var_1307(iplon_var_1329, i_lay_var_1326) * (p_fac00_var_1297(iplon_var_1329, i_lay_var_1326) * absb_var_301(ind0_var_1321, ig_var_1320) + p_fac10_var_1299(iplon_var_1329, i_lay_var_1326) * absb_var_301(ind0_var_1321 + 1, ig_var_1320) + p_fac01_var_1298(iplon_var_1329, i_lay_var_1326) * absb_var_301(ind1_var_1322, ig_var_1320) + p_fac11_var_1300(iplon_var_1329, i_lay_var_1326) * absb_var_301(ind1_var_1322 + 1, ig_var_1320)) + p_colo3_var_1308(iplon_var_1329, i_lay_var_1326) * abso3bc_var_306(ig_var_1320)
        p_taur_var_1318(iplon_var_1329, i_lay_var_1326, ig_var_1320) = z_tauray_var_1336
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol24
SUBROUTINE srtm_taumol25(kidia_var_1337, kfdia_var_1338, klev_var_1339, p_fac00_var_1340, p_fac01_var_1341, p_fac10_var_1342, p_fac11_var_1343, k_jp_var_1344, k_jt_var_1345, k_jt1_var_1346, p_colh2o_var_1347, p_colmol_var_1348, p_colo3_var_1349, k_laytrop_var_1350, p_sfluxzen_var_1351, p_taug_var_1352, p_taur_var_1353, prmu0_var_1354)
  USE yoesrta25, ONLY: absa_var_308, abso3ac_var_311, abso3bc_var_312, layreffr_var_307, raylc_var_310, sfluxrefc_var_309
  USE yoesrtwn, ONLY: nspa_var_333
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1337, kfdia_var_1338
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1339
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1340(kidia_var_1337 : kfdia_var_1338, klev_var_1339)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1341(kidia_var_1337 : kfdia_var_1338, klev_var_1339)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1342(kidia_var_1337 : kfdia_var_1338, klev_var_1339)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1343(kidia_var_1337 : kfdia_var_1338, klev_var_1339)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1344(kidia_var_1337 : kfdia_var_1338, klev_var_1339)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1345(kidia_var_1337 : kfdia_var_1338, klev_var_1339)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1346(kidia_var_1337 : kfdia_var_1338, klev_var_1339)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1347(kidia_var_1337 : kfdia_var_1338, klev_var_1339)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1348(kidia_var_1337 : kfdia_var_1338, klev_var_1339)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1349(kidia_var_1337 : kfdia_var_1338, klev_var_1339)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1350(kidia_var_1337 : kfdia_var_1338)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1351(kidia_var_1337 : kfdia_var_1338, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1352(kidia_var_1337 : kfdia_var_1338, klev_var_1339, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1353(kidia_var_1337 : kfdia_var_1338, klev_var_1339, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1354(kidia_var_1337 : kfdia_var_1338)
  INTEGER(KIND = 4) :: ig_var_1355, ind0_var_1356, ind1_var_1357, i_lay_var_1358, i_laysolfr_var_1359(kidia_var_1337 : kfdia_var_1338), i_nlayers_var_1360, iplon_var_1361
  INTEGER(KIND = 4) :: laytrop_min_var_1362, laytrop_max_var_1363
  REAL(KIND = 8) :: z_tauray_var_1364
  laytrop_min_var_1362 = MINVAL(k_laytrop_var_1350(kidia_var_1337 : kfdia_var_1338))
  laytrop_max_var_1363 = MAXVAL(k_laytrop_var_1350(kidia_var_1337 : kfdia_var_1338))
  i_nlayers_var_1360 = klev_var_1339
  DO iplon_var_1361 = kidia_var_1337, kfdia_var_1338
    i_laysolfr_var_1359(iplon_var_1361) = k_laytrop_var_1350(iplon_var_1361)
  END DO
  DO i_lay_var_1358 = 1, laytrop_min_var_1362
    DO iplon_var_1361 = kidia_var_1337, kfdia_var_1338
      IF (k_jp_var_1344(iplon_var_1361, i_lay_var_1358) < layreffr_var_307 .AND. k_jp_var_1344(iplon_var_1361, i_lay_var_1358 + 1) >= layreffr_var_307) i_laysolfr_var_1359(iplon_var_1361) = MIN(i_lay_var_1358 + 1, k_laytrop_var_1350(iplon_var_1361))
      ind0_var_1356 = ((k_jp_var_1344(iplon_var_1361, i_lay_var_1358) - 1) * 5 + (k_jt_var_1345(iplon_var_1361, i_lay_var_1358) - 1)) * nspa_var_333(25) + 1
      ind1_var_1357 = (k_jp_var_1344(iplon_var_1361, i_lay_var_1358) * 5 + (k_jt1_var_1346(iplon_var_1361, i_lay_var_1358) - 1)) * nspa_var_333(25) + 1
      DO ig_var_1355 = 1, 6
        z_tauray_var_1364 = p_colmol_var_1348(iplon_var_1361, i_lay_var_1358) * raylc_var_310(ig_var_1355)
        p_taug_var_1352(iplon_var_1361, i_lay_var_1358, ig_var_1355) = p_colh2o_var_1347(iplon_var_1361, i_lay_var_1358) * (p_fac00_var_1340(iplon_var_1361, i_lay_var_1358) * absa_var_308(ind0_var_1356, ig_var_1355) + p_fac10_var_1342(iplon_var_1361, i_lay_var_1358) * absa_var_308(ind0_var_1356 + 1, ig_var_1355) + p_fac01_var_1341(iplon_var_1361, i_lay_var_1358) * absa_var_308(ind1_var_1357, ig_var_1355) + p_fac11_var_1343(iplon_var_1361, i_lay_var_1358) * absa_var_308(ind1_var_1357 + 1, ig_var_1355)) + p_colo3_var_1349(iplon_var_1361, i_lay_var_1358) * abso3ac_var_311(ig_var_1355)
        IF (i_lay_var_1358 == i_laysolfr_var_1359(iplon_var_1361)) p_sfluxzen_var_1351(iplon_var_1361, ig_var_1355) = sfluxrefc_var_309(ig_var_1355)
        p_taur_var_1353(iplon_var_1361, i_lay_var_1358, ig_var_1355) = z_tauray_var_1364
      END DO
    END DO
  END DO
  DO i_lay_var_1358 = laytrop_min_var_1362 + 1, laytrop_max_var_1363
    DO iplon_var_1361 = kidia_var_1337, kfdia_var_1338
      IF (i_lay_var_1358 <= k_laytrop_var_1350(iplon_var_1361)) THEN
        IF (k_jp_var_1344(iplon_var_1361, i_lay_var_1358) < layreffr_var_307 .AND. k_jp_var_1344(iplon_var_1361, i_lay_var_1358 + 1) >= layreffr_var_307) i_laysolfr_var_1359(iplon_var_1361) = MIN(i_lay_var_1358 + 1, k_laytrop_var_1350(iplon_var_1361))
        ind0_var_1356 = ((k_jp_var_1344(iplon_var_1361, i_lay_var_1358) - 1) * 5 + (k_jt_var_1345(iplon_var_1361, i_lay_var_1358) - 1)) * nspa_var_333(25) + 1
        ind1_var_1357 = (k_jp_var_1344(iplon_var_1361, i_lay_var_1358) * 5 + (k_jt1_var_1346(iplon_var_1361, i_lay_var_1358) - 1)) * nspa_var_333(25) + 1
        DO ig_var_1355 = 1, 6
          z_tauray_var_1364 = p_colmol_var_1348(iplon_var_1361, i_lay_var_1358) * raylc_var_310(ig_var_1355)
          p_taug_var_1352(iplon_var_1361, i_lay_var_1358, ig_var_1355) = p_colh2o_var_1347(iplon_var_1361, i_lay_var_1358) * (p_fac00_var_1340(iplon_var_1361, i_lay_var_1358) * absa_var_308(ind0_var_1356, ig_var_1355) + p_fac10_var_1342(iplon_var_1361, i_lay_var_1358) * absa_var_308(ind0_var_1356 + 1, ig_var_1355) + p_fac01_var_1341(iplon_var_1361, i_lay_var_1358) * absa_var_308(ind1_var_1357, ig_var_1355) + p_fac11_var_1343(iplon_var_1361, i_lay_var_1358) * absa_var_308(ind1_var_1357 + 1, ig_var_1355)) + p_colo3_var_1349(iplon_var_1361, i_lay_var_1358) * abso3ac_var_311(ig_var_1355)
          IF (i_lay_var_1358 == i_laysolfr_var_1359(iplon_var_1361)) p_sfluxzen_var_1351(iplon_var_1361, ig_var_1355) = sfluxrefc_var_309(ig_var_1355)
          p_taur_var_1353(iplon_var_1361, i_lay_var_1358, ig_var_1355) = z_tauray_var_1364
        END DO
      ELSE
        DO ig_var_1355 = 1, 6
          z_tauray_var_1364 = p_colmol_var_1348(iplon_var_1361, i_lay_var_1358) * raylc_var_310(ig_var_1355)
          p_taug_var_1352(iplon_var_1361, i_lay_var_1358, ig_var_1355) = p_colo3_var_1349(iplon_var_1361, i_lay_var_1358) * abso3bc_var_312(ig_var_1355)
          p_taur_var_1353(iplon_var_1361, i_lay_var_1358, ig_var_1355) = z_tauray_var_1364
        END DO
      END IF
    END DO
  END DO
  DO ig_var_1355 = 1, 6
    DO i_lay_var_1358 = laytrop_max_var_1363 + 1, i_nlayers_var_1360
      DO iplon_var_1361 = kidia_var_1337, kfdia_var_1338
        z_tauray_var_1364 = p_colmol_var_1348(iplon_var_1361, i_lay_var_1358) * raylc_var_310(ig_var_1355)
        p_taug_var_1352(iplon_var_1361, i_lay_var_1358, ig_var_1355) = p_colo3_var_1349(iplon_var_1361, i_lay_var_1358) * abso3bc_var_312(ig_var_1355)
        p_taur_var_1353(iplon_var_1361, i_lay_var_1358, ig_var_1355) = z_tauray_var_1364
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol25
SUBROUTINE srtm_taumol19(kidia_var_1365, kfdia_var_1366, klev_var_1367, p_fac00_var_1368, p_fac01_var_1369, p_fac10_var_1370, p_fac11_var_1371, k_jp_var_1372, k_jt_var_1373, k_jt1_var_1374, p_oneminus_var_1375, p_colh2o_var_1376, p_colco2_var_1377, p_colmol_var_1378, k_laytrop_var_1379, p_selffac_var_1380, p_selffrac_var_1381, k_indself_var_1382, p_forfac_var_1383, p_forfrac_var_1384, k_indfor_var_1385, p_sfluxzen_var_1386, p_taug_var_1387, p_taur_var_1388, prmu0_var_1389)
  USE yoesrta19, ONLY: absa_var_264, absb_var_265, forrefc_var_267, layreffr_var_263, rayl_var_261, selfrefc_var_266, sfluxrefc_var_268, strrat_var_262
  USE yoesrtwn, ONLY: nspa_var_333, nspb_var_334
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1365, kfdia_var_1366
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1367
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1368(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1369(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1370(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1371(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1372(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1373(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1374(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1375(kidia_var_1365 : kfdia_var_1366)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1376(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1377(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1378(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1379(kidia_var_1365 : kfdia_var_1366)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1380(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1381(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1382(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1383(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1384(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1385(kidia_var_1365 : kfdia_var_1366, klev_var_1367)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1386(kidia_var_1365 : kfdia_var_1366, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1387(kidia_var_1365 : kfdia_var_1366, klev_var_1367, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1388(kidia_var_1365 : kfdia_var_1366, klev_var_1367, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1389(kidia_var_1365 : kfdia_var_1366)
  INTEGER(KIND = 4) :: ig_var_1390, ind0_var_1391, ind1_var_1392, inds_var_1393, indf_var_1394, js_var_1395, i_lay_var_1396, i_laysolfr_var_1397(kidia_var_1365 : kfdia_var_1366), i_nlayers_var_1398, iplon_var_1399
  INTEGER(KIND = 4) :: laytrop_min_var_1400, laytrop_max_var_1401
  REAL(KIND = 8) :: z_fs_var_1402, z_speccomb_var_1403, z_specmult_var_1404, z_specparm_var_1405, z_tauray_var_1406
  laytrop_min_var_1400 = MINVAL(k_laytrop_var_1379(kidia_var_1365 : kfdia_var_1366))
  laytrop_max_var_1401 = MAXVAL(k_laytrop_var_1379(kidia_var_1365 : kfdia_var_1366))
  i_nlayers_var_1398 = klev_var_1367
  DO iplon_var_1399 = kidia_var_1365, kfdia_var_1366
    i_laysolfr_var_1397(iplon_var_1399) = k_laytrop_var_1379(iplon_var_1399)
  END DO
  DO i_lay_var_1396 = 1, laytrop_min_var_1400
    DO iplon_var_1399 = kidia_var_1365, kfdia_var_1366
      IF (k_jp_var_1372(iplon_var_1399, i_lay_var_1396) < layreffr_var_263 .AND. k_jp_var_1372(iplon_var_1399, i_lay_var_1396 + 1) >= layreffr_var_263) i_laysolfr_var_1397(iplon_var_1399) = MIN(i_lay_var_1396 + 1, k_laytrop_var_1379(iplon_var_1399))
      z_speccomb_var_1403 = p_colh2o_var_1376(iplon_var_1399, i_lay_var_1396) + strrat_var_262 * p_colco2_var_1377(iplon_var_1399, i_lay_var_1396)
      z_specparm_var_1405 = p_colh2o_var_1376(iplon_var_1399, i_lay_var_1396) / z_speccomb_var_1403
      z_specparm_var_1405 = MIN(p_oneminus_var_1375(iplon_var_1399), z_specparm_var_1405)
      z_specmult_var_1404 = 8.0D0 * (z_specparm_var_1405)
      js_var_1395 = 1 + INT(z_specmult_var_1404)
      z_fs_var_1402 = z_specmult_var_1404 - AINT(z_specmult_var_1404)
      ind0_var_1391 = ((k_jp_var_1372(iplon_var_1399, i_lay_var_1396) - 1) * 5 + (k_jt_var_1373(iplon_var_1399, i_lay_var_1396) - 1)) * nspa_var_333(19) + js_var_1395
      ind1_var_1392 = (k_jp_var_1372(iplon_var_1399, i_lay_var_1396) * 5 + (k_jt1_var_1374(iplon_var_1399, i_lay_var_1396) - 1)) * nspa_var_333(19) + js_var_1395
      inds_var_1393 = k_indself_var_1382(iplon_var_1399, i_lay_var_1396)
      indf_var_1394 = k_indfor_var_1385(iplon_var_1399, i_lay_var_1396)
      z_tauray_var_1406 = p_colmol_var_1378(iplon_var_1399, i_lay_var_1396) * rayl_var_261
      DO ig_var_1390 = 1, 8
        p_taug_var_1387(iplon_var_1399, i_lay_var_1396, ig_var_1390) = z_speccomb_var_1403 * ((1.0D0 - z_fs_var_1402) * (absa_var_264(ind0_var_1391, ig_var_1390) * p_fac00_var_1368(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind0_var_1391 + 9, ig_var_1390) * p_fac10_var_1370(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind1_var_1392, ig_var_1390) * p_fac01_var_1369(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind1_var_1392 + 9, ig_var_1390) * p_fac11_var_1371(iplon_var_1399, i_lay_var_1396)) + z_fs_var_1402 * (absa_var_264(ind0_var_1391 + 1, ig_var_1390) * p_fac00_var_1368(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind0_var_1391 + 10, ig_var_1390) * p_fac10_var_1370(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind1_var_1392 + 1, ig_var_1390) * p_fac01_var_1369(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind1_var_1392 + 10, ig_var_1390) * p_fac11_var_1371(iplon_var_1399, i_lay_var_1396))) + p_colh2o_var_1376(iplon_var_1399, i_lay_var_1396) * (p_selffac_var_1380(iplon_var_1399, i_lay_var_1396) * (selfrefc_var_266(inds_var_1393, ig_var_1390) + p_selffrac_var_1381(iplon_var_1399, i_lay_var_1396) * (selfrefc_var_266(inds_var_1393 + 1, ig_var_1390) - selfrefc_var_266(inds_var_1393, ig_var_1390))) + p_forfac_var_1383(iplon_var_1399, i_lay_var_1396) * (forrefc_var_267(indf_var_1394, ig_var_1390) + p_forfrac_var_1384(iplon_var_1399, i_lay_var_1396) * (forrefc_var_267(indf_var_1394 + 1, ig_var_1390) - forrefc_var_267(indf_var_1394, ig_var_1390))))
        IF (i_lay_var_1396 == i_laysolfr_var_1397(iplon_var_1399)) p_sfluxzen_var_1386(iplon_var_1399, ig_var_1390) = sfluxrefc_var_268(ig_var_1390, js_var_1395) + z_fs_var_1402 * (sfluxrefc_var_268(ig_var_1390, js_var_1395 + 1) - sfluxrefc_var_268(ig_var_1390, js_var_1395))
        p_taur_var_1388(iplon_var_1399, i_lay_var_1396, ig_var_1390) = z_tauray_var_1406
      END DO
    END DO
  END DO
  DO i_lay_var_1396 = laytrop_min_var_1400 + 1, laytrop_max_var_1401
    DO iplon_var_1399 = kidia_var_1365, kfdia_var_1366
      IF (i_lay_var_1396 <= k_laytrop_var_1379(iplon_var_1399)) THEN
        IF (k_jp_var_1372(iplon_var_1399, i_lay_var_1396) < layreffr_var_263 .AND. k_jp_var_1372(iplon_var_1399, i_lay_var_1396 + 1) >= layreffr_var_263) i_laysolfr_var_1397(iplon_var_1399) = MIN(i_lay_var_1396 + 1, k_laytrop_var_1379(iplon_var_1399))
        z_speccomb_var_1403 = p_colh2o_var_1376(iplon_var_1399, i_lay_var_1396) + strrat_var_262 * p_colco2_var_1377(iplon_var_1399, i_lay_var_1396)
        z_specparm_var_1405 = p_colh2o_var_1376(iplon_var_1399, i_lay_var_1396) / z_speccomb_var_1403
        z_specparm_var_1405 = MIN(p_oneminus_var_1375(iplon_var_1399), z_specparm_var_1405)
        z_specmult_var_1404 = 8.0D0 * (z_specparm_var_1405)
        js_var_1395 = 1 + INT(z_specmult_var_1404)
        z_fs_var_1402 = z_specmult_var_1404 - AINT(z_specmult_var_1404)
        ind0_var_1391 = ((k_jp_var_1372(iplon_var_1399, i_lay_var_1396) - 1) * 5 + (k_jt_var_1373(iplon_var_1399, i_lay_var_1396) - 1)) * nspa_var_333(19) + js_var_1395
        ind1_var_1392 = (k_jp_var_1372(iplon_var_1399, i_lay_var_1396) * 5 + (k_jt1_var_1374(iplon_var_1399, i_lay_var_1396) - 1)) * nspa_var_333(19) + js_var_1395
        inds_var_1393 = k_indself_var_1382(iplon_var_1399, i_lay_var_1396)
        indf_var_1394 = k_indfor_var_1385(iplon_var_1399, i_lay_var_1396)
        z_tauray_var_1406 = p_colmol_var_1378(iplon_var_1399, i_lay_var_1396) * rayl_var_261
        DO ig_var_1390 = 1, 8
          p_taug_var_1387(iplon_var_1399, i_lay_var_1396, ig_var_1390) = z_speccomb_var_1403 * ((1.0D0 - z_fs_var_1402) * (absa_var_264(ind0_var_1391, ig_var_1390) * p_fac00_var_1368(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind0_var_1391 + 9, ig_var_1390) * p_fac10_var_1370(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind1_var_1392, ig_var_1390) * p_fac01_var_1369(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind1_var_1392 + 9, ig_var_1390) * p_fac11_var_1371(iplon_var_1399, i_lay_var_1396)) + z_fs_var_1402 * (absa_var_264(ind0_var_1391 + 1, ig_var_1390) * p_fac00_var_1368(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind0_var_1391 + 10, ig_var_1390) * p_fac10_var_1370(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind1_var_1392 + 1, ig_var_1390) * p_fac01_var_1369(iplon_var_1399, i_lay_var_1396) + absa_var_264(ind1_var_1392 + 10, ig_var_1390) * p_fac11_var_1371(iplon_var_1399, i_lay_var_1396))) + p_colh2o_var_1376(iplon_var_1399, i_lay_var_1396) * (p_selffac_var_1380(iplon_var_1399, i_lay_var_1396) * (selfrefc_var_266(inds_var_1393, ig_var_1390) + p_selffrac_var_1381(iplon_var_1399, i_lay_var_1396) * (selfrefc_var_266(inds_var_1393 + 1, ig_var_1390) - selfrefc_var_266(inds_var_1393, ig_var_1390))) + p_forfac_var_1383(iplon_var_1399, i_lay_var_1396) * (forrefc_var_267(indf_var_1394, ig_var_1390) + p_forfrac_var_1384(iplon_var_1399, i_lay_var_1396) * (forrefc_var_267(indf_var_1394 + 1, ig_var_1390) - forrefc_var_267(indf_var_1394, ig_var_1390))))
          IF (i_lay_var_1396 == i_laysolfr_var_1397(iplon_var_1399)) p_sfluxzen_var_1386(iplon_var_1399, ig_var_1390) = sfluxrefc_var_268(ig_var_1390, js_var_1395) + z_fs_var_1402 * (sfluxrefc_var_268(ig_var_1390, js_var_1395 + 1) - sfluxrefc_var_268(ig_var_1390, js_var_1395))
          p_taur_var_1388(iplon_var_1399, i_lay_var_1396, ig_var_1390) = z_tauray_var_1406
        END DO
      ELSE
        ind0_var_1391 = ((k_jp_var_1372(iplon_var_1399, i_lay_var_1396) - 13) * 5 + (k_jt_var_1373(iplon_var_1399, i_lay_var_1396) - 1)) * nspb_var_334(19) + 1
        ind1_var_1392 = ((k_jp_var_1372(iplon_var_1399, i_lay_var_1396) - 12) * 5 + (k_jt1_var_1374(iplon_var_1399, i_lay_var_1396) - 1)) * nspb_var_334(19) + 1
        z_tauray_var_1406 = p_colmol_var_1378(iplon_var_1399, i_lay_var_1396) * rayl_var_261
        DO ig_var_1390 = 1, 8
          p_taug_var_1387(iplon_var_1399, i_lay_var_1396, ig_var_1390) = p_colco2_var_1377(iplon_var_1399, i_lay_var_1396) * (p_fac00_var_1368(iplon_var_1399, i_lay_var_1396) * absb_var_265(ind0_var_1391, ig_var_1390) + p_fac10_var_1370(iplon_var_1399, i_lay_var_1396) * absb_var_265(ind0_var_1391 + 1, ig_var_1390) + p_fac01_var_1369(iplon_var_1399, i_lay_var_1396) * absb_var_265(ind1_var_1392, ig_var_1390) + p_fac11_var_1371(iplon_var_1399, i_lay_var_1396) * absb_var_265(ind1_var_1392 + 1, ig_var_1390))
          p_taur_var_1388(iplon_var_1399, i_lay_var_1396, ig_var_1390) = z_tauray_var_1406
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1396 = laytrop_max_var_1401 + 1, i_nlayers_var_1398
    DO iplon_var_1399 = kidia_var_1365, kfdia_var_1366
      ind0_var_1391 = ((k_jp_var_1372(iplon_var_1399, i_lay_var_1396) - 13) * 5 + (k_jt_var_1373(iplon_var_1399, i_lay_var_1396) - 1)) * nspb_var_334(19) + 1
      ind1_var_1392 = ((k_jp_var_1372(iplon_var_1399, i_lay_var_1396) - 12) * 5 + (k_jt1_var_1374(iplon_var_1399, i_lay_var_1396) - 1)) * nspb_var_334(19) + 1
      z_tauray_var_1406 = p_colmol_var_1378(iplon_var_1399, i_lay_var_1396) * rayl_var_261
      DO ig_var_1390 = 1, 8
        p_taug_var_1387(iplon_var_1399, i_lay_var_1396, ig_var_1390) = p_colco2_var_1377(iplon_var_1399, i_lay_var_1396) * (p_fac00_var_1368(iplon_var_1399, i_lay_var_1396) * absb_var_265(ind0_var_1391, ig_var_1390) + p_fac10_var_1370(iplon_var_1399, i_lay_var_1396) * absb_var_265(ind0_var_1391 + 1, ig_var_1390) + p_fac01_var_1369(iplon_var_1399, i_lay_var_1396) * absb_var_265(ind1_var_1392, ig_var_1390) + p_fac11_var_1371(iplon_var_1399, i_lay_var_1396) * absb_var_265(ind1_var_1392 + 1, ig_var_1390))
        p_taur_var_1388(iplon_var_1399, i_lay_var_1396, ig_var_1390) = z_tauray_var_1406
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol19
SUBROUTINE srtm_taumol21(kidia_var_1407, kfdia_var_1408, klev_var_1409, p_fac00_var_1410, p_fac01_var_1411, p_fac10_var_1412, p_fac11_var_1413, k_jp_var_1414, k_jt_var_1415, k_jt1_var_1416, p_oneminus_var_1417, p_colh2o_var_1418, p_colco2_var_1419, p_colmol_var_1420, k_laytrop_var_1421, p_selffac_var_1422, p_selffrac_var_1423, k_indself_var_1424, p_forfac_var_1425, p_forfrac_var_1426, k_indfor_var_1427, p_sfluxzen_var_1428, p_taug_var_1429, p_taur_var_1430, prmu0_var_1431)
  USE yoesrta21, ONLY: absa_var_279, absb_var_280, forrefc_var_282, layreffr_var_278, rayl_var_276, selfrefc_var_281, sfluxrefc_var_283, strrat_var_277
  USE yoesrtwn, ONLY: nspa_var_333, nspb_var_334
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1407, kfdia_var_1408
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1409
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1410(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1411(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1412(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1413(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1414(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1415(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1416(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1417(kidia_var_1407 : kfdia_var_1408)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1418(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1419(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1420(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1421(kidia_var_1407 : kfdia_var_1408)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1422(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1423(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1424(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1425(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1426(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1427(kidia_var_1407 : kfdia_var_1408, klev_var_1409)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1428(kidia_var_1407 : kfdia_var_1408, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1429(kidia_var_1407 : kfdia_var_1408, klev_var_1409, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1430(kidia_var_1407 : kfdia_var_1408, klev_var_1409, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1431(kidia_var_1407 : kfdia_var_1408)
  INTEGER(KIND = 4) :: ig_var_1432, ind0_var_1433, ind1_var_1434, inds_var_1435, indf_var_1436, js_var_1437, i_lay_var_1438, i_laysolfr_var_1439(kidia_var_1407 : kfdia_var_1408), i_nlayers_var_1440, iplon_var_1441
  INTEGER(KIND = 4) :: laytrop_min_var_1442, laytrop_max_var_1443
  REAL(KIND = 8) :: z_fs_var_1444, z_speccomb_var_1445, z_specmult_var_1446, z_specparm_var_1447, z_tauray_var_1448
  laytrop_min_var_1442 = MINVAL(k_laytrop_var_1421(kidia_var_1407 : kfdia_var_1408))
  laytrop_max_var_1443 = MAXVAL(k_laytrop_var_1421(kidia_var_1407 : kfdia_var_1408))
  i_nlayers_var_1440 = klev_var_1409
  DO iplon_var_1441 = kidia_var_1407, kfdia_var_1408
    i_laysolfr_var_1439(iplon_var_1441) = k_laytrop_var_1421(iplon_var_1441)
  END DO
  DO i_lay_var_1438 = 1, laytrop_min_var_1442
    DO iplon_var_1441 = kidia_var_1407, kfdia_var_1408
      IF (k_jp_var_1414(iplon_var_1441, i_lay_var_1438) < layreffr_var_278 .AND. k_jp_var_1414(iplon_var_1441, i_lay_var_1438 + 1) >= layreffr_var_278) i_laysolfr_var_1439(iplon_var_1441) = MIN(i_lay_var_1438 + 1, k_laytrop_var_1421(iplon_var_1441))
      z_speccomb_var_1445 = p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) + strrat_var_277 * p_colco2_var_1419(iplon_var_1441, i_lay_var_1438)
      z_specparm_var_1447 = p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) / z_speccomb_var_1445
      z_specparm_var_1447 = MIN(p_oneminus_var_1417(iplon_var_1441), z_specparm_var_1447)
      z_specmult_var_1446 = 8.0D0 * (z_specparm_var_1447)
      js_var_1437 = 1 + INT(z_specmult_var_1446)
      z_fs_var_1444 = z_specmult_var_1446 - AINT(z_specmult_var_1446)
      ind0_var_1433 = ((k_jp_var_1414(iplon_var_1441, i_lay_var_1438) - 1) * 5 + (k_jt_var_1415(iplon_var_1441, i_lay_var_1438) - 1)) * nspa_var_333(21) + js_var_1437
      ind1_var_1434 = (k_jp_var_1414(iplon_var_1441, i_lay_var_1438) * 5 + (k_jt1_var_1416(iplon_var_1441, i_lay_var_1438) - 1)) * nspa_var_333(21) + js_var_1437
      inds_var_1435 = k_indself_var_1424(iplon_var_1441, i_lay_var_1438)
      indf_var_1436 = k_indfor_var_1427(iplon_var_1441, i_lay_var_1438)
      z_tauray_var_1448 = p_colmol_var_1420(iplon_var_1441, i_lay_var_1438) * rayl_var_276
      DO ig_var_1432 = 1, 10
        p_taug_var_1429(iplon_var_1441, i_lay_var_1438, ig_var_1432) = z_speccomb_var_1445 * ((1.0D0 - z_fs_var_1444) * (absa_var_279(ind0_var_1433, ig_var_1432) * p_fac00_var_1410(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind0_var_1433 + 9, ig_var_1432) * p_fac10_var_1412(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind1_var_1434, ig_var_1432) * p_fac01_var_1411(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind1_var_1434 + 9, ig_var_1432) * p_fac11_var_1413(iplon_var_1441, i_lay_var_1438)) + z_fs_var_1444 * (absa_var_279(ind0_var_1433 + 1, ig_var_1432) * p_fac00_var_1410(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind0_var_1433 + 10, ig_var_1432) * p_fac10_var_1412(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind1_var_1434 + 1, ig_var_1432) * p_fac01_var_1411(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind1_var_1434 + 10, ig_var_1432) * p_fac11_var_1413(iplon_var_1441, i_lay_var_1438))) + p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) * (p_selffac_var_1422(iplon_var_1441, i_lay_var_1438) * (selfrefc_var_281(inds_var_1435, ig_var_1432) + p_selffrac_var_1423(iplon_var_1441, i_lay_var_1438) * (selfrefc_var_281(inds_var_1435 + 1, ig_var_1432) - selfrefc_var_281(inds_var_1435, ig_var_1432))) + p_forfac_var_1425(iplon_var_1441, i_lay_var_1438) * (forrefc_var_282(indf_var_1436, ig_var_1432) + p_forfrac_var_1426(iplon_var_1441, i_lay_var_1438) * (forrefc_var_282(indf_var_1436 + 1, ig_var_1432) - forrefc_var_282(indf_var_1436, ig_var_1432))))
        IF (i_lay_var_1438 == i_laysolfr_var_1439(iplon_var_1441)) p_sfluxzen_var_1428(iplon_var_1441, ig_var_1432) = sfluxrefc_var_283(ig_var_1432, js_var_1437) + z_fs_var_1444 * (sfluxrefc_var_283(ig_var_1432, js_var_1437 + 1) - sfluxrefc_var_283(ig_var_1432, js_var_1437))
        p_taur_var_1430(iplon_var_1441, i_lay_var_1438, ig_var_1432) = z_tauray_var_1448
      END DO
    END DO
  END DO
  DO i_lay_var_1438 = laytrop_min_var_1442 + 1, laytrop_max_var_1443
    DO iplon_var_1441 = kidia_var_1407, kfdia_var_1408
      IF (i_lay_var_1438 <= k_laytrop_var_1421(iplon_var_1441)) THEN
        IF (k_jp_var_1414(iplon_var_1441, i_lay_var_1438) < layreffr_var_278 .AND. k_jp_var_1414(iplon_var_1441, i_lay_var_1438 + 1) >= layreffr_var_278) i_laysolfr_var_1439(iplon_var_1441) = MIN(i_lay_var_1438 + 1, k_laytrop_var_1421(iplon_var_1441))
        z_speccomb_var_1445 = p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) + strrat_var_277 * p_colco2_var_1419(iplon_var_1441, i_lay_var_1438)
        z_specparm_var_1447 = p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) / z_speccomb_var_1445
        z_specparm_var_1447 = MIN(p_oneminus_var_1417(iplon_var_1441), z_specparm_var_1447)
        z_specmult_var_1446 = 8.0D0 * (z_specparm_var_1447)
        js_var_1437 = 1 + INT(z_specmult_var_1446)
        z_fs_var_1444 = z_specmult_var_1446 - AINT(z_specmult_var_1446)
        ind0_var_1433 = ((k_jp_var_1414(iplon_var_1441, i_lay_var_1438) - 1) * 5 + (k_jt_var_1415(iplon_var_1441, i_lay_var_1438) - 1)) * nspa_var_333(21) + js_var_1437
        ind1_var_1434 = (k_jp_var_1414(iplon_var_1441, i_lay_var_1438) * 5 + (k_jt1_var_1416(iplon_var_1441, i_lay_var_1438) - 1)) * nspa_var_333(21) + js_var_1437
        inds_var_1435 = k_indself_var_1424(iplon_var_1441, i_lay_var_1438)
        indf_var_1436 = k_indfor_var_1427(iplon_var_1441, i_lay_var_1438)
        z_tauray_var_1448 = p_colmol_var_1420(iplon_var_1441, i_lay_var_1438) * rayl_var_276
        DO ig_var_1432 = 1, 10
          p_taug_var_1429(iplon_var_1441, i_lay_var_1438, ig_var_1432) = z_speccomb_var_1445 * ((1.0D0 - z_fs_var_1444) * (absa_var_279(ind0_var_1433, ig_var_1432) * p_fac00_var_1410(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind0_var_1433 + 9, ig_var_1432) * p_fac10_var_1412(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind1_var_1434, ig_var_1432) * p_fac01_var_1411(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind1_var_1434 + 9, ig_var_1432) * p_fac11_var_1413(iplon_var_1441, i_lay_var_1438)) + z_fs_var_1444 * (absa_var_279(ind0_var_1433 + 1, ig_var_1432) * p_fac00_var_1410(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind0_var_1433 + 10, ig_var_1432) * p_fac10_var_1412(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind1_var_1434 + 1, ig_var_1432) * p_fac01_var_1411(iplon_var_1441, i_lay_var_1438) + absa_var_279(ind1_var_1434 + 10, ig_var_1432) * p_fac11_var_1413(iplon_var_1441, i_lay_var_1438))) + p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) * (p_selffac_var_1422(iplon_var_1441, i_lay_var_1438) * (selfrefc_var_281(inds_var_1435, ig_var_1432) + p_selffrac_var_1423(iplon_var_1441, i_lay_var_1438) * (selfrefc_var_281(inds_var_1435 + 1, ig_var_1432) - selfrefc_var_281(inds_var_1435, ig_var_1432))) + p_forfac_var_1425(iplon_var_1441, i_lay_var_1438) * (forrefc_var_282(indf_var_1436, ig_var_1432) + p_forfrac_var_1426(iplon_var_1441, i_lay_var_1438) * (forrefc_var_282(indf_var_1436 + 1, ig_var_1432) - forrefc_var_282(indf_var_1436, ig_var_1432))))
          IF (i_lay_var_1438 == i_laysolfr_var_1439(iplon_var_1441)) p_sfluxzen_var_1428(iplon_var_1441, ig_var_1432) = sfluxrefc_var_283(ig_var_1432, js_var_1437) + z_fs_var_1444 * (sfluxrefc_var_283(ig_var_1432, js_var_1437 + 1) - sfluxrefc_var_283(ig_var_1432, js_var_1437))
          p_taur_var_1430(iplon_var_1441, i_lay_var_1438, ig_var_1432) = z_tauray_var_1448
        END DO
      ELSE
        z_speccomb_var_1445 = p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) + strrat_var_277 * p_colco2_var_1419(iplon_var_1441, i_lay_var_1438)
        z_specparm_var_1447 = p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) / z_speccomb_var_1445
        z_specparm_var_1447 = MIN(p_oneminus_var_1417(iplon_var_1441), z_specparm_var_1447)
        z_specmult_var_1446 = 4.0D0 * (z_specparm_var_1447)
        js_var_1437 = 1 + INT(z_specmult_var_1446)
        z_fs_var_1444 = z_specmult_var_1446 - AINT(z_specmult_var_1446)
        ind0_var_1433 = ((k_jp_var_1414(iplon_var_1441, i_lay_var_1438) - 13) * 5 + (k_jt_var_1415(iplon_var_1441, i_lay_var_1438) - 1)) * nspb_var_334(21) + js_var_1437
        ind1_var_1434 = ((k_jp_var_1414(iplon_var_1441, i_lay_var_1438) - 12) * 5 + (k_jt1_var_1416(iplon_var_1441, i_lay_var_1438) - 1)) * nspb_var_334(21) + js_var_1437
        indf_var_1436 = k_indfor_var_1427(iplon_var_1441, i_lay_var_1438)
        z_tauray_var_1448 = p_colmol_var_1420(iplon_var_1441, i_lay_var_1438) * rayl_var_276
        DO ig_var_1432 = 1, 10
          p_taug_var_1429(iplon_var_1441, i_lay_var_1438, ig_var_1432) = z_speccomb_var_1445 * ((1.0D0 - z_fs_var_1444) * (absb_var_280(ind0_var_1433, ig_var_1432) * p_fac00_var_1410(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind0_var_1433 + 5, ig_var_1432) * p_fac10_var_1412(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind1_var_1434, ig_var_1432) * p_fac01_var_1411(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind1_var_1434 + 5, ig_var_1432) * p_fac11_var_1413(iplon_var_1441, i_lay_var_1438)) + z_fs_var_1444 * (absb_var_280(ind0_var_1433 + 1, ig_var_1432) * p_fac00_var_1410(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind0_var_1433 + 6, ig_var_1432) * p_fac10_var_1412(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind1_var_1434 + 1, ig_var_1432) * p_fac01_var_1411(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind1_var_1434 + 6, ig_var_1432) * p_fac11_var_1413(iplon_var_1441, i_lay_var_1438))) + p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) * p_forfac_var_1425(iplon_var_1441, i_lay_var_1438) * (forrefc_var_282(indf_var_1436, ig_var_1432) + p_forfrac_var_1426(iplon_var_1441, i_lay_var_1438) * (forrefc_var_282(indf_var_1436 + 1, ig_var_1432) - forrefc_var_282(indf_var_1436, ig_var_1432)))
          p_taur_var_1430(iplon_var_1441, i_lay_var_1438, ig_var_1432) = z_tauray_var_1448
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1438 = laytrop_max_var_1443 + 1, i_nlayers_var_1440
    DO iplon_var_1441 = kidia_var_1407, kfdia_var_1408
      z_speccomb_var_1445 = p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) + strrat_var_277 * p_colco2_var_1419(iplon_var_1441, i_lay_var_1438)
      z_specparm_var_1447 = p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) / z_speccomb_var_1445
      z_specparm_var_1447 = MIN(p_oneminus_var_1417(iplon_var_1441), z_specparm_var_1447)
      z_specmult_var_1446 = 4.0D0 * (z_specparm_var_1447)
      js_var_1437 = 1 + INT(z_specmult_var_1446)
      z_fs_var_1444 = z_specmult_var_1446 - AINT(z_specmult_var_1446)
      ind0_var_1433 = ((k_jp_var_1414(iplon_var_1441, i_lay_var_1438) - 13) * 5 + (k_jt_var_1415(iplon_var_1441, i_lay_var_1438) - 1)) * nspb_var_334(21) + js_var_1437
      ind1_var_1434 = ((k_jp_var_1414(iplon_var_1441, i_lay_var_1438) - 12) * 5 + (k_jt1_var_1416(iplon_var_1441, i_lay_var_1438) - 1)) * nspb_var_334(21) + js_var_1437
      indf_var_1436 = k_indfor_var_1427(iplon_var_1441, i_lay_var_1438)
      z_tauray_var_1448 = p_colmol_var_1420(iplon_var_1441, i_lay_var_1438) * rayl_var_276
      DO ig_var_1432 = 1, 10
        p_taug_var_1429(iplon_var_1441, i_lay_var_1438, ig_var_1432) = z_speccomb_var_1445 * ((1.0D0 - z_fs_var_1444) * (absb_var_280(ind0_var_1433, ig_var_1432) * p_fac00_var_1410(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind0_var_1433 + 5, ig_var_1432) * p_fac10_var_1412(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind1_var_1434, ig_var_1432) * p_fac01_var_1411(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind1_var_1434 + 5, ig_var_1432) * p_fac11_var_1413(iplon_var_1441, i_lay_var_1438)) + z_fs_var_1444 * (absb_var_280(ind0_var_1433 + 1, ig_var_1432) * p_fac00_var_1410(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind0_var_1433 + 6, ig_var_1432) * p_fac10_var_1412(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind1_var_1434 + 1, ig_var_1432) * p_fac01_var_1411(iplon_var_1441, i_lay_var_1438) + absb_var_280(ind1_var_1434 + 6, ig_var_1432) * p_fac11_var_1413(iplon_var_1441, i_lay_var_1438))) + p_colh2o_var_1418(iplon_var_1441, i_lay_var_1438) * p_forfac_var_1425(iplon_var_1441, i_lay_var_1438) * (forrefc_var_282(indf_var_1436, ig_var_1432) + p_forfrac_var_1426(iplon_var_1441, i_lay_var_1438) * (forrefc_var_282(indf_var_1436 + 1, ig_var_1432) - forrefc_var_282(indf_var_1436, ig_var_1432)))
        p_taur_var_1430(iplon_var_1441, i_lay_var_1438, ig_var_1432) = z_tauray_var_1448
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol21
SUBROUTINE srtm_taumol20(kidia_var_1449, kfdia_var_1450, klev_var_1451, p_fac00_var_1452, p_fac01_var_1453, p_fac10_var_1454, p_fac11_var_1455, k_jp_var_1456, k_jt_var_1457, k_jt1_var_1458, p_colh2o_var_1459, p_colch4_var_1460, p_colmol_var_1461, k_laytrop_var_1462, p_selffac_var_1463, p_selffrac_var_1464, k_indself_var_1465, p_forfac_var_1466, p_forfrac_var_1467, k_indfor_var_1468, p_sfluxzen_var_1469, p_taug_var_1470, p_taur_var_1471, prmu0_var_1472)
  USE yoesrta20, ONLY: absa_var_271, absb_var_272, absch4c, forrefc_var_274, layreffr_var_270, rayl_var_269, selfrefc_var_273, sfluxrefc_var_275
  USE yoesrtwn, ONLY: nspa_var_333, nspb_var_334
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1449, kfdia_var_1450
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1451
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1452(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1453(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1454(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1455(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1456(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1457(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1458(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1459(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_1460(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1461(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1462(kidia_var_1449 : kfdia_var_1450)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1463(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1464(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1465(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1466(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1467(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1468(kidia_var_1449 : kfdia_var_1450, klev_var_1451)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1469(kidia_var_1449 : kfdia_var_1450, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1470(kidia_var_1449 : kfdia_var_1450, klev_var_1451, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1471(kidia_var_1449 : kfdia_var_1450, klev_var_1451, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1472(kidia_var_1449 : kfdia_var_1450)
  INTEGER(KIND = 4) :: ig_var_1473, ind0_var_1474, ind1_var_1475, inds_var_1476, indf_var_1477, i_lay_var_1478, i_laysolfr_var_1479(kidia_var_1449 : kfdia_var_1450), i_nlayers_var_1480, iplon_var_1481
  INTEGER(KIND = 4) :: laytrop_min_var_1482, laytrop_max_var_1483
  REAL(KIND = 8) :: z_tauray_var_1484
  laytrop_min_var_1482 = MINVAL(k_laytrop_var_1462(kidia_var_1449 : kfdia_var_1450))
  laytrop_max_var_1483 = MAXVAL(k_laytrop_var_1462(kidia_var_1449 : kfdia_var_1450))
  i_nlayers_var_1480 = klev_var_1451
  DO iplon_var_1481 = kidia_var_1449, kfdia_var_1450
    i_laysolfr_var_1479(iplon_var_1481) = k_laytrop_var_1462(iplon_var_1481)
  END DO
  DO i_lay_var_1478 = 1, laytrop_min_var_1482
    DO iplon_var_1481 = kidia_var_1449, kfdia_var_1450
      IF (k_jp_var_1456(iplon_var_1481, i_lay_var_1478) < layreffr_var_270 .AND. k_jp_var_1456(iplon_var_1481, i_lay_var_1478 + 1) >= layreffr_var_270) i_laysolfr_var_1479(iplon_var_1481) = MIN(i_lay_var_1478 + 1, k_laytrop_var_1462(iplon_var_1481))
      ind0_var_1474 = ((k_jp_var_1456(iplon_var_1481, i_lay_var_1478) - 1) * 5 + (k_jt_var_1457(iplon_var_1481, i_lay_var_1478) - 1)) * nspa_var_333(20) + 1
      ind1_var_1475 = (k_jp_var_1456(iplon_var_1481, i_lay_var_1478) * 5 + (k_jt1_var_1458(iplon_var_1481, i_lay_var_1478) - 1)) * nspa_var_333(20) + 1
      inds_var_1476 = k_indself_var_1465(iplon_var_1481, i_lay_var_1478)
      indf_var_1477 = k_indfor_var_1468(iplon_var_1481, i_lay_var_1478)
      z_tauray_var_1484 = p_colmol_var_1461(iplon_var_1481, i_lay_var_1478) * rayl_var_269
      DO ig_var_1473 = 1, 10
        p_taug_var_1470(iplon_var_1481, i_lay_var_1478, ig_var_1473) = p_colh2o_var_1459(iplon_var_1481, i_lay_var_1478) * ((p_fac00_var_1452(iplon_var_1481, i_lay_var_1478) * absa_var_271(ind0_var_1474, ig_var_1473) + p_fac10_var_1454(iplon_var_1481, i_lay_var_1478) * absa_var_271(ind0_var_1474 + 1, ig_var_1473) + p_fac01_var_1453(iplon_var_1481, i_lay_var_1478) * absa_var_271(ind1_var_1475, ig_var_1473) + p_fac11_var_1455(iplon_var_1481, i_lay_var_1478) * absa_var_271(ind1_var_1475 + 1, ig_var_1473)) + p_selffac_var_1463(iplon_var_1481, i_lay_var_1478) * (selfrefc_var_273(inds_var_1476, ig_var_1473) + p_selffrac_var_1464(iplon_var_1481, i_lay_var_1478) * (selfrefc_var_273(inds_var_1476 + 1, ig_var_1473) - selfrefc_var_273(inds_var_1476, ig_var_1473))) + p_forfac_var_1466(iplon_var_1481, i_lay_var_1478) * (forrefc_var_274(indf_var_1477, ig_var_1473) + p_forfrac_var_1467(iplon_var_1481, i_lay_var_1478) * (forrefc_var_274(indf_var_1477 + 1, ig_var_1473) - forrefc_var_274(indf_var_1477, ig_var_1473)))) + p_colch4_var_1460(iplon_var_1481, i_lay_var_1478) * absch4c(ig_var_1473)
        p_taur_var_1471(iplon_var_1481, i_lay_var_1478, ig_var_1473) = z_tauray_var_1484
        IF (i_lay_var_1478 == i_laysolfr_var_1479(iplon_var_1481)) p_sfluxzen_var_1469(iplon_var_1481, ig_var_1473) = sfluxrefc_var_275(ig_var_1473)
      END DO
    END DO
  END DO
  DO i_lay_var_1478 = laytrop_min_var_1482 + 1, laytrop_max_var_1483
    DO iplon_var_1481 = kidia_var_1449, kfdia_var_1450
      IF (i_lay_var_1478 <= k_laytrop_var_1462(iplon_var_1481)) THEN
        IF (k_jp_var_1456(iplon_var_1481, i_lay_var_1478) < layreffr_var_270 .AND. k_jp_var_1456(iplon_var_1481, i_lay_var_1478 + 1) >= layreffr_var_270) i_laysolfr_var_1479(iplon_var_1481) = MIN(i_lay_var_1478 + 1, k_laytrop_var_1462(iplon_var_1481))
        ind0_var_1474 = ((k_jp_var_1456(iplon_var_1481, i_lay_var_1478) - 1) * 5 + (k_jt_var_1457(iplon_var_1481, i_lay_var_1478) - 1)) * nspa_var_333(20) + 1
        ind1_var_1475 = (k_jp_var_1456(iplon_var_1481, i_lay_var_1478) * 5 + (k_jt1_var_1458(iplon_var_1481, i_lay_var_1478) - 1)) * nspa_var_333(20) + 1
        inds_var_1476 = k_indself_var_1465(iplon_var_1481, i_lay_var_1478)
        indf_var_1477 = k_indfor_var_1468(iplon_var_1481, i_lay_var_1478)
        z_tauray_var_1484 = p_colmol_var_1461(iplon_var_1481, i_lay_var_1478) * rayl_var_269
        DO ig_var_1473 = 1, 10
          p_taug_var_1470(iplon_var_1481, i_lay_var_1478, ig_var_1473) = p_colh2o_var_1459(iplon_var_1481, i_lay_var_1478) * ((p_fac00_var_1452(iplon_var_1481, i_lay_var_1478) * absa_var_271(ind0_var_1474, ig_var_1473) + p_fac10_var_1454(iplon_var_1481, i_lay_var_1478) * absa_var_271(ind0_var_1474 + 1, ig_var_1473) + p_fac01_var_1453(iplon_var_1481, i_lay_var_1478) * absa_var_271(ind1_var_1475, ig_var_1473) + p_fac11_var_1455(iplon_var_1481, i_lay_var_1478) * absa_var_271(ind1_var_1475 + 1, ig_var_1473)) + p_selffac_var_1463(iplon_var_1481, i_lay_var_1478) * (selfrefc_var_273(inds_var_1476, ig_var_1473) + p_selffrac_var_1464(iplon_var_1481, i_lay_var_1478) * (selfrefc_var_273(inds_var_1476 + 1, ig_var_1473) - selfrefc_var_273(inds_var_1476, ig_var_1473))) + p_forfac_var_1466(iplon_var_1481, i_lay_var_1478) * (forrefc_var_274(indf_var_1477, ig_var_1473) + p_forfrac_var_1467(iplon_var_1481, i_lay_var_1478) * (forrefc_var_274(indf_var_1477 + 1, ig_var_1473) - forrefc_var_274(indf_var_1477, ig_var_1473)))) + p_colch4_var_1460(iplon_var_1481, i_lay_var_1478) * absch4c(ig_var_1473)
          p_taur_var_1471(iplon_var_1481, i_lay_var_1478, ig_var_1473) = z_tauray_var_1484
          IF (i_lay_var_1478 == i_laysolfr_var_1479(iplon_var_1481)) p_sfluxzen_var_1469(iplon_var_1481, ig_var_1473) = sfluxrefc_var_275(ig_var_1473)
        END DO
      ELSE
        ind0_var_1474 = ((k_jp_var_1456(iplon_var_1481, i_lay_var_1478) - 13) * 5 + (k_jt_var_1457(iplon_var_1481, i_lay_var_1478) - 1)) * nspb_var_334(20) + 1
        ind1_var_1475 = ((k_jp_var_1456(iplon_var_1481, i_lay_var_1478) - 12) * 5 + (k_jt1_var_1458(iplon_var_1481, i_lay_var_1478) - 1)) * nspb_var_334(20) + 1
        indf_var_1477 = k_indfor_var_1468(iplon_var_1481, i_lay_var_1478)
        z_tauray_var_1484 = p_colmol_var_1461(iplon_var_1481, i_lay_var_1478) * rayl_var_269
        DO ig_var_1473 = 1, 10
          p_taug_var_1470(iplon_var_1481, i_lay_var_1478, ig_var_1473) = p_colh2o_var_1459(iplon_var_1481, i_lay_var_1478) * (p_fac00_var_1452(iplon_var_1481, i_lay_var_1478) * absb_var_272(ind0_var_1474, ig_var_1473) + p_fac10_var_1454(iplon_var_1481, i_lay_var_1478) * absb_var_272(ind0_var_1474 + 1, ig_var_1473) + p_fac01_var_1453(iplon_var_1481, i_lay_var_1478) * absb_var_272(ind1_var_1475, ig_var_1473) + p_fac11_var_1455(iplon_var_1481, i_lay_var_1478) * absb_var_272(ind1_var_1475 + 1, ig_var_1473) + p_forfac_var_1466(iplon_var_1481, i_lay_var_1478) * (forrefc_var_274(indf_var_1477, ig_var_1473) + p_forfrac_var_1467(iplon_var_1481, i_lay_var_1478) * (forrefc_var_274(indf_var_1477 + 1, ig_var_1473) - forrefc_var_274(indf_var_1477, ig_var_1473)))) + p_colch4_var_1460(iplon_var_1481, i_lay_var_1478) * absch4c(ig_var_1473)
          p_taur_var_1471(iplon_var_1481, i_lay_var_1478, ig_var_1473) = z_tauray_var_1484
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1478 = laytrop_max_var_1483 + 1, i_nlayers_var_1480
    DO iplon_var_1481 = kidia_var_1449, kfdia_var_1450
      ind0_var_1474 = ((k_jp_var_1456(iplon_var_1481, i_lay_var_1478) - 13) * 5 + (k_jt_var_1457(iplon_var_1481, i_lay_var_1478) - 1)) * nspb_var_334(20) + 1
      ind1_var_1475 = ((k_jp_var_1456(iplon_var_1481, i_lay_var_1478) - 12) * 5 + (k_jt1_var_1458(iplon_var_1481, i_lay_var_1478) - 1)) * nspb_var_334(20) + 1
      indf_var_1477 = k_indfor_var_1468(iplon_var_1481, i_lay_var_1478)
      z_tauray_var_1484 = p_colmol_var_1461(iplon_var_1481, i_lay_var_1478) * rayl_var_269
      DO ig_var_1473 = 1, 10
        p_taug_var_1470(iplon_var_1481, i_lay_var_1478, ig_var_1473) = p_colh2o_var_1459(iplon_var_1481, i_lay_var_1478) * (p_fac00_var_1452(iplon_var_1481, i_lay_var_1478) * absb_var_272(ind0_var_1474, ig_var_1473) + p_fac10_var_1454(iplon_var_1481, i_lay_var_1478) * absb_var_272(ind0_var_1474 + 1, ig_var_1473) + p_fac01_var_1453(iplon_var_1481, i_lay_var_1478) * absb_var_272(ind1_var_1475, ig_var_1473) + p_fac11_var_1455(iplon_var_1481, i_lay_var_1478) * absb_var_272(ind1_var_1475 + 1, ig_var_1473) + p_forfac_var_1466(iplon_var_1481, i_lay_var_1478) * (forrefc_var_274(indf_var_1477, ig_var_1473) + p_forfrac_var_1467(iplon_var_1481, i_lay_var_1478) * (forrefc_var_274(indf_var_1477 + 1, ig_var_1473) - forrefc_var_274(indf_var_1477, ig_var_1473)))) + p_colch4_var_1460(iplon_var_1481, i_lay_var_1478) * absch4c(ig_var_1473)
        p_taur_var_1471(iplon_var_1481, i_lay_var_1478, ig_var_1473) = z_tauray_var_1484
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol20
SUBROUTINE srtm_taumol22(kidia_var_1485, kfdia_var_1486, klev_var_1487, p_fac00_var_1488, p_fac01_var_1489, p_fac10_var_1490, p_fac11_var_1491, k_jp_var_1492, k_jt_var_1493, k_jt1_var_1494, p_oneminus_var_1495, p_colh2o_var_1496, p_colmol_var_1497, p_colo2_var_1498, k_laytrop_var_1499, p_selffac_var_1500, p_selffrac_var_1501, k_indself_var_1502, p_forfac_var_1503, p_forfrac_var_1504, k_indfor_var_1505, p_sfluxzen_var_1506, p_taug_var_1507, p_taur_var_1508, prmu0_var_1509)
  USE yoesrta22, ONLY: absa_var_287, absb_var_288, forrefc_var_290, layreffr_var_286, rayl_var_284, selfrefc_var_289, sfluxrefc_var_291, strrat_var_285
  USE yoesrtwn, ONLY: nspa_var_333, nspb_var_334
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1485, kfdia_var_1486
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1487
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1488(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1489(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1490(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1491(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1492(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1493(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1494(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1495(kidia_var_1485 : kfdia_var_1486)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1496(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1497(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_1498(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1499(kidia_var_1485 : kfdia_var_1486)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1500(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1501(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1502(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1503(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1504(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1505(kidia_var_1485 : kfdia_var_1486, klev_var_1487)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1506(kidia_var_1485 : kfdia_var_1486, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1507(kidia_var_1485 : kfdia_var_1486, klev_var_1487, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1508(kidia_var_1485 : kfdia_var_1486, klev_var_1487, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1509(kidia_var_1485 : kfdia_var_1486)
  INTEGER(KIND = 4) :: ig_var_1510, ind0_var_1511, ind1_var_1512, inds_var_1513, indf_var_1514, js_var_1515, i_lay_var_1516, i_laysolfr_var_1517(kidia_var_1485 : kfdia_var_1486), i_nlayers_var_1518, iplon_var_1519
  INTEGER(KIND = 4) :: laytrop_min_var_1520, laytrop_max_var_1521
  REAL(KIND = 8) :: z_fs_var_1522, z_speccomb_var_1523, z_specmult_var_1524, z_specparm_var_1525, z_tauray_var_1526, z_o2adj, z_o2cont
  laytrop_min_var_1520 = MINVAL(k_laytrop_var_1499(kidia_var_1485 : kfdia_var_1486))
  laytrop_max_var_1521 = MAXVAL(k_laytrop_var_1499(kidia_var_1485 : kfdia_var_1486))
  i_nlayers_var_1518 = klev_var_1487
  z_o2adj = 1.6D0
  DO iplon_var_1519 = kidia_var_1485, kfdia_var_1486
    i_laysolfr_var_1517(iplon_var_1519) = k_laytrop_var_1499(iplon_var_1519)
  END DO
  DO i_lay_var_1516 = 1, laytrop_min_var_1520
    DO iplon_var_1519 = kidia_var_1485, kfdia_var_1486
      IF (k_jp_var_1492(iplon_var_1519, i_lay_var_1516) < layreffr_var_286 .AND. k_jp_var_1492(iplon_var_1519, i_lay_var_1516 + 1) >= layreffr_var_286) i_laysolfr_var_1517(iplon_var_1519) = MIN(i_lay_var_1516 + 1, k_laytrop_var_1499(iplon_var_1519))
      z_o2cont = 0.000435D0 * p_colo2_var_1498(iplon_var_1519, i_lay_var_1516) / (700.0D0)
      z_speccomb_var_1523 = p_colh2o_var_1496(iplon_var_1519, i_lay_var_1516) + 1.6D0 * strrat_var_285 * p_colo2_var_1498(iplon_var_1519, i_lay_var_1516)
      z_specparm_var_1525 = p_colh2o_var_1496(iplon_var_1519, i_lay_var_1516) / z_speccomb_var_1523
      z_specparm_var_1525 = MIN(p_oneminus_var_1495(iplon_var_1519), z_specparm_var_1525)
      z_specmult_var_1524 = 8.0D0 * (z_specparm_var_1525)
      js_var_1515 = 1 + INT(z_specmult_var_1524)
      z_fs_var_1522 = z_specmult_var_1524 - AINT(z_specmult_var_1524)
      ind0_var_1511 = ((k_jp_var_1492(iplon_var_1519, i_lay_var_1516) - 1) * 5 + (k_jt_var_1493(iplon_var_1519, i_lay_var_1516) - 1)) * nspa_var_333(22) + js_var_1515
      ind1_var_1512 = (k_jp_var_1492(iplon_var_1519, i_lay_var_1516) * 5 + (k_jt1_var_1494(iplon_var_1519, i_lay_var_1516) - 1)) * nspa_var_333(22) + js_var_1515
      inds_var_1513 = k_indself_var_1502(iplon_var_1519, i_lay_var_1516)
      indf_var_1514 = k_indfor_var_1505(iplon_var_1519, i_lay_var_1516)
      z_tauray_var_1526 = p_colmol_var_1497(iplon_var_1519, i_lay_var_1516) * rayl_var_284
      DO ig_var_1510 = 1, 2
        p_taug_var_1507(iplon_var_1519, i_lay_var_1516, ig_var_1510) = z_speccomb_var_1523 * ((1.0D0 - z_fs_var_1522) * (absa_var_287(ind0_var_1511, ig_var_1510) * p_fac00_var_1488(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind0_var_1511 + 9, ig_var_1510) * p_fac10_var_1490(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind1_var_1512, ig_var_1510) * p_fac01_var_1489(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind1_var_1512 + 9, ig_var_1510) * p_fac11_var_1491(iplon_var_1519, i_lay_var_1516)) + z_fs_var_1522 * (absa_var_287(ind0_var_1511 + 1, ig_var_1510) * p_fac00_var_1488(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind0_var_1511 + 10, ig_var_1510) * p_fac10_var_1490(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind1_var_1512 + 1, ig_var_1510) * p_fac01_var_1489(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind1_var_1512 + 10, ig_var_1510) * p_fac11_var_1491(iplon_var_1519, i_lay_var_1516))) + p_colh2o_var_1496(iplon_var_1519, i_lay_var_1516) * (p_selffac_var_1500(iplon_var_1519, i_lay_var_1516) * (selfrefc_var_289(inds_var_1513, ig_var_1510) + p_selffrac_var_1501(iplon_var_1519, i_lay_var_1516) * (selfrefc_var_289(inds_var_1513 + 1, ig_var_1510) - selfrefc_var_289(inds_var_1513, ig_var_1510))) + p_forfac_var_1503(iplon_var_1519, i_lay_var_1516) * (forrefc_var_290(indf_var_1514, ig_var_1510) + p_forfrac_var_1504(iplon_var_1519, i_lay_var_1516) * (forrefc_var_290(indf_var_1514 + 1, ig_var_1510) - forrefc_var_290(indf_var_1514, ig_var_1510)))) + z_o2cont
        IF (i_lay_var_1516 == i_laysolfr_var_1517(iplon_var_1519)) p_sfluxzen_var_1506(iplon_var_1519, ig_var_1510) = sfluxrefc_var_291(ig_var_1510, js_var_1515) + z_fs_var_1522 * (sfluxrefc_var_291(ig_var_1510, js_var_1515 + 1) - sfluxrefc_var_291(ig_var_1510, js_var_1515))
        p_taur_var_1508(iplon_var_1519, i_lay_var_1516, ig_var_1510) = z_tauray_var_1526
      END DO
    END DO
  END DO
  DO i_lay_var_1516 = laytrop_min_var_1520 + 1, laytrop_max_var_1521
    DO iplon_var_1519 = kidia_var_1485, kfdia_var_1486
      IF (i_lay_var_1516 <= k_laytrop_var_1499(iplon_var_1519)) THEN
        IF (k_jp_var_1492(iplon_var_1519, i_lay_var_1516) < layreffr_var_286 .AND. k_jp_var_1492(iplon_var_1519, i_lay_var_1516 + 1) >= layreffr_var_286) i_laysolfr_var_1517(iplon_var_1519) = MIN(i_lay_var_1516 + 1, k_laytrop_var_1499(iplon_var_1519))
        z_o2cont = 0.000435D0 * p_colo2_var_1498(iplon_var_1519, i_lay_var_1516) / (700.0D0)
        z_speccomb_var_1523 = p_colh2o_var_1496(iplon_var_1519, i_lay_var_1516) + z_o2adj * strrat_var_285 * p_colo2_var_1498(iplon_var_1519, i_lay_var_1516)
        z_specparm_var_1525 = p_colh2o_var_1496(iplon_var_1519, i_lay_var_1516) / z_speccomb_var_1523
        z_specparm_var_1525 = MIN(p_oneminus_var_1495(iplon_var_1519), z_specparm_var_1525)
        z_specmult_var_1524 = 8.0D0 * (z_specparm_var_1525)
        js_var_1515 = 1 + INT(z_specmult_var_1524)
        z_fs_var_1522 = z_specmult_var_1524 - AINT(z_specmult_var_1524)
        ind0_var_1511 = ((k_jp_var_1492(iplon_var_1519, i_lay_var_1516) - 1) * 5 + (k_jt_var_1493(iplon_var_1519, i_lay_var_1516) - 1)) * nspa_var_333(22) + js_var_1515
        ind1_var_1512 = (k_jp_var_1492(iplon_var_1519, i_lay_var_1516) * 5 + (k_jt1_var_1494(iplon_var_1519, i_lay_var_1516) - 1)) * nspa_var_333(22) + js_var_1515
        inds_var_1513 = k_indself_var_1502(iplon_var_1519, i_lay_var_1516)
        indf_var_1514 = k_indfor_var_1505(iplon_var_1519, i_lay_var_1516)
        z_tauray_var_1526 = p_colmol_var_1497(iplon_var_1519, i_lay_var_1516) * rayl_var_284
        DO ig_var_1510 = 1, 2
          p_taug_var_1507(iplon_var_1519, i_lay_var_1516, ig_var_1510) = z_speccomb_var_1523 * ((1.0D0 - z_fs_var_1522) * (absa_var_287(ind0_var_1511, ig_var_1510) * p_fac00_var_1488(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind0_var_1511 + 9, ig_var_1510) * p_fac10_var_1490(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind1_var_1512, ig_var_1510) * p_fac01_var_1489(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind1_var_1512 + 9, ig_var_1510) * p_fac11_var_1491(iplon_var_1519, i_lay_var_1516)) + z_fs_var_1522 * (absa_var_287(ind0_var_1511 + 1, ig_var_1510) * p_fac00_var_1488(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind0_var_1511 + 10, ig_var_1510) * p_fac10_var_1490(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind1_var_1512 + 1, ig_var_1510) * p_fac01_var_1489(iplon_var_1519, i_lay_var_1516) + absa_var_287(ind1_var_1512 + 10, ig_var_1510) * p_fac11_var_1491(iplon_var_1519, i_lay_var_1516))) + p_colh2o_var_1496(iplon_var_1519, i_lay_var_1516) * (p_selffac_var_1500(iplon_var_1519, i_lay_var_1516) * (selfrefc_var_289(inds_var_1513, ig_var_1510) + p_selffrac_var_1501(iplon_var_1519, i_lay_var_1516) * (selfrefc_var_289(inds_var_1513 + 1, ig_var_1510) - selfrefc_var_289(inds_var_1513, ig_var_1510))) + p_forfac_var_1503(iplon_var_1519, i_lay_var_1516) * (forrefc_var_290(indf_var_1514, ig_var_1510) + p_forfrac_var_1504(iplon_var_1519, i_lay_var_1516) * (forrefc_var_290(indf_var_1514 + 1, ig_var_1510) - forrefc_var_290(indf_var_1514, ig_var_1510)))) + z_o2cont
          IF (i_lay_var_1516 == i_laysolfr_var_1517(iplon_var_1519)) p_sfluxzen_var_1506(iplon_var_1519, ig_var_1510) = sfluxrefc_var_291(ig_var_1510, js_var_1515) + z_fs_var_1522 * (sfluxrefc_var_291(ig_var_1510, js_var_1515 + 1) - sfluxrefc_var_291(ig_var_1510, js_var_1515))
          p_taur_var_1508(iplon_var_1519, i_lay_var_1516, ig_var_1510) = z_tauray_var_1526
        END DO
      ELSE
        z_o2cont = 0.000435D0 * p_colo2_var_1498(iplon_var_1519, i_lay_var_1516) / (700.0D0)
        ind0_var_1511 = ((k_jp_var_1492(iplon_var_1519, i_lay_var_1516) - 13) * 5 + (k_jt_var_1493(iplon_var_1519, i_lay_var_1516) - 1)) * nspb_var_334(22) + 1
        ind1_var_1512 = ((k_jp_var_1492(iplon_var_1519, i_lay_var_1516) - 12) * 5 + (k_jt1_var_1494(iplon_var_1519, i_lay_var_1516) - 1)) * nspb_var_334(22) + 1
        z_tauray_var_1526 = p_colmol_var_1497(iplon_var_1519, i_lay_var_1516) * rayl_var_284
        DO ig_var_1510 = 1, 2
          p_taug_var_1507(iplon_var_1519, i_lay_var_1516, ig_var_1510) = p_colo2_var_1498(iplon_var_1519, i_lay_var_1516) * z_o2adj * (p_fac00_var_1488(iplon_var_1519, i_lay_var_1516) * absb_var_288(ind0_var_1511, ig_var_1510) + p_fac10_var_1490(iplon_var_1519, i_lay_var_1516) * absb_var_288(ind0_var_1511 + 1, ig_var_1510) + p_fac01_var_1489(iplon_var_1519, i_lay_var_1516) * absb_var_288(ind1_var_1512, ig_var_1510) + p_fac11_var_1491(iplon_var_1519, i_lay_var_1516) * absb_var_288(ind1_var_1512 + 1, ig_var_1510)) + z_o2cont
          p_taur_var_1508(iplon_var_1519, i_lay_var_1516, ig_var_1510) = z_tauray_var_1526
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1516 = laytrop_max_var_1521 + 1, i_nlayers_var_1518
    DO iplon_var_1519 = kidia_var_1485, kfdia_var_1486
      z_o2cont = 0.000435D0 * p_colo2_var_1498(iplon_var_1519, i_lay_var_1516) / (700.0D0)
      ind0_var_1511 = ((k_jp_var_1492(iplon_var_1519, i_lay_var_1516) - 13) * 5 + (k_jt_var_1493(iplon_var_1519, i_lay_var_1516) - 1)) * nspb_var_334(22) + 1
      ind1_var_1512 = ((k_jp_var_1492(iplon_var_1519, i_lay_var_1516) - 12) * 5 + (k_jt1_var_1494(iplon_var_1519, i_lay_var_1516) - 1)) * nspb_var_334(22) + 1
      z_tauray_var_1526 = p_colmol_var_1497(iplon_var_1519, i_lay_var_1516) * rayl_var_284
      DO ig_var_1510 = 1, 2
        p_taug_var_1507(iplon_var_1519, i_lay_var_1516, ig_var_1510) = p_colo2_var_1498(iplon_var_1519, i_lay_var_1516) * 1.6D0 * (p_fac00_var_1488(iplon_var_1519, i_lay_var_1516) * absb_var_288(ind0_var_1511, ig_var_1510) + p_fac10_var_1490(iplon_var_1519, i_lay_var_1516) * absb_var_288(ind0_var_1511 + 1, ig_var_1510) + p_fac01_var_1489(iplon_var_1519, i_lay_var_1516) * absb_var_288(ind1_var_1512, ig_var_1510) + p_fac11_var_1491(iplon_var_1519, i_lay_var_1516) * absb_var_288(ind1_var_1512 + 1, ig_var_1510)) + z_o2cont
        p_taur_var_1508(iplon_var_1519, i_lay_var_1516, ig_var_1510) = z_tauray_var_1526
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol22
SUBROUTINE srtm_taumol23(kidia_var_1527, kfdia_var_1528, klev_var_1529, p_fac00_var_1530, p_fac01_var_1531, p_fac10_var_1532, p_fac11_var_1533, k_jp_var_1534, k_jt_var_1535, k_jt1_var_1536, p_colh2o_var_1537, p_colmol_var_1538, k_laytrop_var_1539, p_selffac_var_1540, p_selffrac_var_1541, k_indself_var_1542, p_forfac_var_1543, p_forfrac_var_1544, k_indfor_var_1545, p_sfluxzen_var_1546, p_taug_var_1547, p_taur_var_1548, prmu0_var_1549)
  USE yoesrta23, ONLY: absa_var_293, forrefc_var_295, givfac, layreffr_var_292, raylc_var_297, selfrefc_var_294, sfluxrefc_var_296
  USE yoesrtwn, ONLY: nspa_var_333
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1527, kfdia_var_1528
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1529
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1530(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1531(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1532(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1533(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1534(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1535(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1536(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1537(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1538(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1539(kidia_var_1527 : kfdia_var_1528)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1540(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1541(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1542(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1543(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1544(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1545(kidia_var_1527 : kfdia_var_1528, klev_var_1529)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1546(kidia_var_1527 : kfdia_var_1528, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1547(kidia_var_1527 : kfdia_var_1528, klev_var_1529, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1548(kidia_var_1527 : kfdia_var_1528, klev_var_1529, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1549(kidia_var_1527 : kfdia_var_1528)
  INTEGER(KIND = 4) :: ig_var_1550, ind0_var_1551, ind1_var_1552, inds_var_1553, indf_var_1554, i_lay_var_1555, i_laysolfr_var_1556(kidia_var_1527 : kfdia_var_1528), i_nlayers_var_1557, iplon_var_1558
  INTEGER(KIND = 4) :: laytrop_min_var_1559, laytrop_max_var_1560
  REAL(KIND = 8) :: z_tauray_var_1561
  laytrop_min_var_1559 = MINVAL(k_laytrop_var_1539(kidia_var_1527 : kfdia_var_1528))
  laytrop_max_var_1560 = MAXVAL(k_laytrop_var_1539(kidia_var_1527 : kfdia_var_1528))
  i_nlayers_var_1557 = klev_var_1529
  DO iplon_var_1558 = kidia_var_1527, kfdia_var_1528
    i_laysolfr_var_1556(iplon_var_1558) = k_laytrop_var_1539(iplon_var_1558)
  END DO
  DO i_lay_var_1555 = 1, laytrop_min_var_1559
    DO iplon_var_1558 = kidia_var_1527, kfdia_var_1528
      IF (k_jp_var_1534(iplon_var_1558, i_lay_var_1555) < layreffr_var_292 .AND. k_jp_var_1534(iplon_var_1558, i_lay_var_1555 + 1) >= layreffr_var_292) i_laysolfr_var_1556(iplon_var_1558) = MIN(i_lay_var_1555 + 1, k_laytrop_var_1539(iplon_var_1558))
      ind0_var_1551 = ((k_jp_var_1534(iplon_var_1558, i_lay_var_1555) - 1) * 5 + (k_jt_var_1535(iplon_var_1558, i_lay_var_1555) - 1)) * nspa_var_333(23) + 1
      ind1_var_1552 = (k_jp_var_1534(iplon_var_1558, i_lay_var_1555) * 5 + (k_jt1_var_1536(iplon_var_1558, i_lay_var_1555) - 1)) * nspa_var_333(23) + 1
      inds_var_1553 = k_indself_var_1542(iplon_var_1558, i_lay_var_1555)
      indf_var_1554 = k_indfor_var_1545(iplon_var_1558, i_lay_var_1555)
      DO ig_var_1550 = 1, 10
        z_tauray_var_1561 = p_colmol_var_1538(iplon_var_1558, i_lay_var_1555) * raylc_var_297(ig_var_1550)
        p_taug_var_1547(iplon_var_1558, i_lay_var_1555, ig_var_1550) = p_colh2o_var_1537(iplon_var_1558, i_lay_var_1555) * (givfac * (p_fac00_var_1530(iplon_var_1558, i_lay_var_1555) * absa_var_293(ind0_var_1551, ig_var_1550) + p_fac10_var_1532(iplon_var_1558, i_lay_var_1555) * absa_var_293(ind0_var_1551 + 1, ig_var_1550) + p_fac01_var_1531(iplon_var_1558, i_lay_var_1555) * absa_var_293(ind1_var_1552, ig_var_1550) + p_fac11_var_1533(iplon_var_1558, i_lay_var_1555) * absa_var_293(ind1_var_1552 + 1, ig_var_1550)) + p_selffac_var_1540(iplon_var_1558, i_lay_var_1555) * (selfrefc_var_294(inds_var_1553, ig_var_1550) + p_selffrac_var_1541(iplon_var_1558, i_lay_var_1555) * (selfrefc_var_294(inds_var_1553 + 1, ig_var_1550) - selfrefc_var_294(inds_var_1553, ig_var_1550))) + p_forfac_var_1543(iplon_var_1558, i_lay_var_1555) * (forrefc_var_295(indf_var_1554, ig_var_1550) + p_forfrac_var_1544(iplon_var_1558, i_lay_var_1555) * (forrefc_var_295(indf_var_1554 + 1, ig_var_1550) - forrefc_var_295(indf_var_1554, ig_var_1550))))
        IF (i_lay_var_1555 == i_laysolfr_var_1556(iplon_var_1558)) p_sfluxzen_var_1546(iplon_var_1558, ig_var_1550) = sfluxrefc_var_296(ig_var_1550)
        p_taur_var_1548(iplon_var_1558, i_lay_var_1555, ig_var_1550) = z_tauray_var_1561
      END DO
    END DO
  END DO
  DO i_lay_var_1555 = laytrop_min_var_1559 + 1, laytrop_max_var_1560
    DO iplon_var_1558 = kidia_var_1527, kfdia_var_1528
      IF (i_lay_var_1555 <= k_laytrop_var_1539(iplon_var_1558)) THEN
        IF (k_jp_var_1534(iplon_var_1558, i_lay_var_1555) < layreffr_var_292 .AND. k_jp_var_1534(iplon_var_1558, i_lay_var_1555 + 1) >= layreffr_var_292) i_laysolfr_var_1556(iplon_var_1558) = MIN(i_lay_var_1555 + 1, k_laytrop_var_1539(iplon_var_1558))
        ind0_var_1551 = ((k_jp_var_1534(iplon_var_1558, i_lay_var_1555) - 1) * 5 + (k_jt_var_1535(iplon_var_1558, i_lay_var_1555) - 1)) * nspa_var_333(23) + 1
        ind1_var_1552 = (k_jp_var_1534(iplon_var_1558, i_lay_var_1555) * 5 + (k_jt1_var_1536(iplon_var_1558, i_lay_var_1555) - 1)) * nspa_var_333(23) + 1
        inds_var_1553 = k_indself_var_1542(iplon_var_1558, i_lay_var_1555)
        indf_var_1554 = k_indfor_var_1545(iplon_var_1558, i_lay_var_1555)
        DO ig_var_1550 = 1, 10
          z_tauray_var_1561 = p_colmol_var_1538(iplon_var_1558, i_lay_var_1555) * raylc_var_297(ig_var_1550)
          p_taug_var_1547(iplon_var_1558, i_lay_var_1555, ig_var_1550) = p_colh2o_var_1537(iplon_var_1558, i_lay_var_1555) * (givfac * (p_fac00_var_1530(iplon_var_1558, i_lay_var_1555) * absa_var_293(ind0_var_1551, ig_var_1550) + p_fac10_var_1532(iplon_var_1558, i_lay_var_1555) * absa_var_293(ind0_var_1551 + 1, ig_var_1550) + p_fac01_var_1531(iplon_var_1558, i_lay_var_1555) * absa_var_293(ind1_var_1552, ig_var_1550) + p_fac11_var_1533(iplon_var_1558, i_lay_var_1555) * absa_var_293(ind1_var_1552 + 1, ig_var_1550)) + p_selffac_var_1540(iplon_var_1558, i_lay_var_1555) * (selfrefc_var_294(inds_var_1553, ig_var_1550) + p_selffrac_var_1541(iplon_var_1558, i_lay_var_1555) * (selfrefc_var_294(inds_var_1553 + 1, ig_var_1550) - selfrefc_var_294(inds_var_1553, ig_var_1550))) + p_forfac_var_1543(iplon_var_1558, i_lay_var_1555) * (forrefc_var_295(indf_var_1554, ig_var_1550) + p_forfrac_var_1544(iplon_var_1558, i_lay_var_1555) * (forrefc_var_295(indf_var_1554 + 1, ig_var_1550) - forrefc_var_295(indf_var_1554, ig_var_1550))))
          IF (i_lay_var_1555 == i_laysolfr_var_1556(iplon_var_1558)) p_sfluxzen_var_1546(iplon_var_1558, ig_var_1550) = sfluxrefc_var_296(ig_var_1550)
          p_taur_var_1548(iplon_var_1558, i_lay_var_1555, ig_var_1550) = z_tauray_var_1561
        END DO
      ELSE
        DO ig_var_1550 = 1, 10
          p_taug_var_1547(iplon_var_1558, i_lay_var_1555, ig_var_1550) = 0.0D0
          p_taur_var_1548(iplon_var_1558, i_lay_var_1555, ig_var_1550) = p_colmol_var_1538(iplon_var_1558, i_lay_var_1555) * raylc_var_297(ig_var_1550)
        END DO
      END IF
    END DO
  END DO
  DO ig_var_1550 = 1, 10
    DO i_lay_var_1555 = laytrop_max_var_1560 + 1, i_nlayers_var_1557
      DO iplon_var_1558 = kidia_var_1527, kfdia_var_1528
        p_taug_var_1547(iplon_var_1558, i_lay_var_1555, ig_var_1550) = 0.0D0
        p_taur_var_1548(iplon_var_1558, i_lay_var_1555, ig_var_1550) = p_colmol_var_1538(iplon_var_1558, i_lay_var_1555) * raylc_var_297(ig_var_1550)
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol23
SUBROUTINE rrtm_taumol5(kidia_var_1562, kfdia_var_1563, klev_var_1564, taug_var_1565, wx_var_1566, p_tauaerl_var_1567, fac00_var_1568, fac01_var_1569, fac10_var_1570, fac11_var_1571, forfac_var_1590, forfrac_var_1589, indfor_var_1588, jp_var_1572, jt_var_1573, jt1_var_1574, oneminus_var_1575, colh2o_var_1576, colco2_var_1577, colo3_var_1578, laytrop_var_1579, selffac_var_1580, selffrac_var_1581, indself_var_1582, fracs_var_1583, rat_h2oco2_var_1584, rat_h2oco2_1_var_1585, rat_o3co2_var_1586, rat_o3co2_1_var_1587, minorfrac_var_1591, indminor_var_1592)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrtm, ONLY: ng5
  USE yoerrta5, ONLY: absa_var_195, absb_var_196, ccl4, forref_var_199, fracrefa_var_193, fracrefb_var_194, ka_mo3_var_197, selfref_var_198
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1562
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1563
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1564
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1565(kidia_var_1562 : kfdia_var_1563, 140, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: wx_var_1566(kidia_var_1562 : kfdia_var_1563, 4, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1567(kidia_var_1562 : kfdia_var_1563, klev_var_1564, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1568(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1569(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1570(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1571(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1572(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1573(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1574(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1575
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1576(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1577(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1578(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1579(kidia_var_1562 : kfdia_var_1563)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1580(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1581(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1582(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1583(kidia_var_1562 : kfdia_var_1563, 140, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1584(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1585(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_var_1586(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_1_var_1587(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1588(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1589(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1590(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1591(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1592(kidia_var_1562 : kfdia_var_1563, klev_var_1564)
  REAL(KIND = 8) :: speccomb_var_1593, speccomb1_var_1594, speccomb_mo3, speccomb_planck_var_1595
  INTEGER(KIND = 4) :: ind0_var_1596, ind1_var_1597, inds_var_1598, indf_var_1599, indm_var_1600
  INTEGER(KIND = 4) :: ig_var_1601, js_var_1602, lay_var_1603, js1_var_1604, jpl_var_1605, jmo3
  REAL(KIND = 8) :: refrat_planck_a_var_1606, refrat_planck_b_var_1607, refrat_m_a_var_1608
  REAL(KIND = 8) :: fac000_var_1609, fac100_var_1610, fac200_var_1611, fac010_var_1612, fac110_var_1613, fac210_var_1614, fac001_var_1615, fac101_var_1616, fac201_var_1617, fac011_var_1618, fac111_var_1619, fac211_var_1620
  REAL(KIND = 8) :: p_var_1621, p4_var_1622, fk0_var_1623, fk1_var_1624, fk2_var_1625
  REAL(KIND = 8) :: taufor_var_1626, tauself_var_1627, tau_major_var_1628(16), tau_major1_var_1629(16), o3m1, o3m2, abso3_var_1630
  REAL(KIND = 8) :: fs_var_1631, specmult_var_1632, specparm_var_1633, fs1_var_1634, specmult1_var_1635, specparm1_var_1636, fpl_var_1637, specmult_planck_var_1638, specparm_planck_var_1639, fmo3, specmult_mo3, specparm_mo3
  INTEGER(KIND = 4) :: laytrop_min_var_1640, laytrop_max_var_1641
  INTEGER(KIND = 4) :: ixc_var_1642(klev_var_1564), ixlow_var_1643(kfdia_var_1563, klev_var_1564), ixhigh_var_1644(kfdia_var_1563, klev_var_1564)
  INTEGER(KIND = 4) :: ich_var_1645, icl_var_1646, ixc0_var_1647, ixp_var_1648, jc_var_1649, jl_var_1650
  laytrop_min_var_1640 = MINVAL(laytrop_var_1579)
  laytrop_max_var_1641 = MAXVAL(laytrop_var_1579)
  ixlow_var_1643 = 0
  ixhigh_var_1644 = 0
  ixc_var_1642 = 0
  DO lay_var_1603 = laytrop_min_var_1640 + 1, laytrop_max_var_1641
    icl_var_1646 = 0
    ich_var_1645 = 0
    DO jc_var_1649 = kidia_var_1562, kfdia_var_1563
      IF (lay_var_1603 <= laytrop_var_1579(jc_var_1649)) THEN
        icl_var_1646 = icl_var_1646 + 1
        ixlow_var_1643(icl_var_1646, lay_var_1603) = jc_var_1649
      ELSE
        ich_var_1645 = ich_var_1645 + 1
        ixhigh_var_1644(ich_var_1645, lay_var_1603) = jc_var_1649
      END IF
    END DO
    ixc_var_1642(lay_var_1603) = icl_var_1646
  END DO
  refrat_planck_a_var_1606 = chi_mls(1, 5) / chi_mls(2, 5)
  refrat_planck_b_var_1607 = chi_mls(3, 43) / chi_mls(2, 43)
  refrat_m_a_var_1608 = chi_mls(1, 7) / chi_mls(2, 7)
  DO lay_var_1603 = 1, laytrop_min_var_1640
    DO jl_var_1650 = kidia_var_1562, kfdia_var_1563
      speccomb_var_1593 = colh2o_var_1576(jl_var_1650, lay_var_1603) + rat_h2oco2_var_1584(jl_var_1650, lay_var_1603) * colco2_var_1577(jl_var_1650, lay_var_1603)
      specparm_var_1633 = MIN(colh2o_var_1576(jl_var_1650, lay_var_1603) / speccomb_var_1593, oneminus_var_1575)
      specmult_var_1632 = 8.0D0 * (specparm_var_1633)
      js_var_1602 = 1 + INT(specmult_var_1632)
      fs_var_1631 = ((specmult_var_1632) - AINT((specmult_var_1632)))
      speccomb1_var_1594 = colh2o_var_1576(jl_var_1650, lay_var_1603) + rat_h2oco2_1_var_1585(jl_var_1650, lay_var_1603) * colco2_var_1577(jl_var_1650, lay_var_1603)
      specparm1_var_1636 = MIN(colh2o_var_1576(jl_var_1650, lay_var_1603) / speccomb1_var_1594, oneminus_var_1575)
      specmult1_var_1635 = 8.0D0 * (specparm1_var_1636)
      js1_var_1604 = 1 + INT(specmult1_var_1635)
      fs1_var_1634 = ((specmult1_var_1635) - AINT((specmult1_var_1635)))
      speccomb_mo3 = colh2o_var_1576(jl_var_1650, lay_var_1603) + refrat_m_a_var_1608 * colco2_var_1577(jl_var_1650, lay_var_1603)
      specparm_mo3 = MIN(colh2o_var_1576(jl_var_1650, lay_var_1603) / speccomb_mo3, oneminus_var_1575)
      specmult_mo3 = 8.0D0 * specparm_mo3
      jmo3 = 1 + INT(specmult_mo3)
      fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))
      speccomb_planck_var_1595 = colh2o_var_1576(jl_var_1650, lay_var_1603) + refrat_planck_a_var_1606 * colco2_var_1577(jl_var_1650, lay_var_1603)
      specparm_planck_var_1639 = MIN(colh2o_var_1576(jl_var_1650, lay_var_1603) / speccomb_planck_var_1595, oneminus_var_1575)
      specmult_planck_var_1638 = 8.0D0 * specparm_planck_var_1639
      jpl_var_1605 = 1 + INT(specmult_planck_var_1638)
      fpl_var_1637 = ((specmult_planck_var_1638) - AINT((specmult_planck_var_1638)))
      ind0_var_1596 = ((jp_var_1572(jl_var_1650, lay_var_1603) - 1) * 5 + (jt_var_1573(jl_var_1650, lay_var_1603) - 1)) * nspa_var_236(5) + js_var_1602
      ind1_var_1597 = (jp_var_1572(jl_var_1650, lay_var_1603) * 5 + (jt1_var_1574(jl_var_1650, lay_var_1603) - 1)) * nspa_var_236(5) + js1_var_1604
      inds_var_1598 = indself_var_1582(jl_var_1650, lay_var_1603)
      indf_var_1599 = indfor_var_1588(jl_var_1650, lay_var_1603)
      indm_var_1600 = indminor_var_1592(jl_var_1650, lay_var_1603)
      IF (specparm_var_1633 .LT. 0.125D0) THEN
        p_var_1621 = fs_var_1631 - 1.0D0
        p4_var_1622 = p_var_1621 ** 4
        fk0_var_1623 = p4_var_1622
        fk1_var_1624 = 1.0D0 - p_var_1621 - 2.0D0 * p4_var_1622
        fk2_var_1625 = p_var_1621 + p4_var_1622
        fac000_var_1609 = fk0_var_1623 * fac00_var_1568(jl_var_1650, lay_var_1603)
        fac100_var_1610 = fk1_var_1624 * fac00_var_1568(jl_var_1650, lay_var_1603)
        fac200_var_1611 = fk2_var_1625 * fac00_var_1568(jl_var_1650, lay_var_1603)
        fac010_var_1612 = fk0_var_1623 * fac10_var_1570(jl_var_1650, lay_var_1603)
        fac110_var_1613 = fk1_var_1624 * fac10_var_1570(jl_var_1650, lay_var_1603)
        fac210_var_1614 = fk2_var_1625 * fac10_var_1570(jl_var_1650, lay_var_1603)
      ELSE IF (specparm_var_1633 .GT. 0.875D0) THEN
        p_var_1621 = - fs_var_1631
        p4_var_1622 = p_var_1621 ** 4
        fk0_var_1623 = p4_var_1622
        fk1_var_1624 = 1.0D0 - p_var_1621 - 2.0D0 * p4_var_1622
        fk2_var_1625 = p_var_1621 + p4_var_1622
        fac000_var_1609 = fk0_var_1623 * fac00_var_1568(jl_var_1650, lay_var_1603)
        fac100_var_1610 = fk1_var_1624 * fac00_var_1568(jl_var_1650, lay_var_1603)
        fac200_var_1611 = fk2_var_1625 * fac00_var_1568(jl_var_1650, lay_var_1603)
        fac010_var_1612 = fk0_var_1623 * fac10_var_1570(jl_var_1650, lay_var_1603)
        fac110_var_1613 = fk1_var_1624 * fac10_var_1570(jl_var_1650, lay_var_1603)
        fac210_var_1614 = fk2_var_1625 * fac10_var_1570(jl_var_1650, lay_var_1603)
      ELSE
        fac000_var_1609 = (1.0D0 - fs_var_1631) * fac00_var_1568(jl_var_1650, lay_var_1603)
        fac010_var_1612 = (1.0D0 - fs_var_1631) * fac10_var_1570(jl_var_1650, lay_var_1603)
        fac100_var_1610 = fs_var_1631 * fac00_var_1568(jl_var_1650, lay_var_1603)
        fac110_var_1613 = fs_var_1631 * fac10_var_1570(jl_var_1650, lay_var_1603)
        fac200_var_1611 = 0.0D0
        fac210_var_1614 = 0.0D0
      END IF
      IF (specparm1_var_1636 .LT. 0.125D0) THEN
        p_var_1621 = fs1_var_1634 - 1.0D0
        p4_var_1622 = p_var_1621 ** 4
        fk0_var_1623 = p4_var_1622
        fk1_var_1624 = 1.0D0 - p_var_1621 - 2.0D0 * p4_var_1622
        fk2_var_1625 = p_var_1621 + p4_var_1622
        fac001_var_1615 = fk0_var_1623 * fac01_var_1569(jl_var_1650, lay_var_1603)
        fac101_var_1616 = fk1_var_1624 * fac01_var_1569(jl_var_1650, lay_var_1603)
        fac201_var_1617 = fk2_var_1625 * fac01_var_1569(jl_var_1650, lay_var_1603)
        fac011_var_1618 = fk0_var_1623 * fac11_var_1571(jl_var_1650, lay_var_1603)
        fac111_var_1619 = fk1_var_1624 * fac11_var_1571(jl_var_1650, lay_var_1603)
        fac211_var_1620 = fk2_var_1625 * fac11_var_1571(jl_var_1650, lay_var_1603)
      ELSE IF (specparm1_var_1636 .GT. 0.875D0) THEN
        p_var_1621 = - fs1_var_1634
        p4_var_1622 = p_var_1621 ** 4
        fk0_var_1623 = p4_var_1622
        fk1_var_1624 = 1.0D0 - p_var_1621 - 2.0D0 * p4_var_1622
        fk2_var_1625 = p_var_1621 + p4_var_1622
        fac001_var_1615 = fk0_var_1623 * fac01_var_1569(jl_var_1650, lay_var_1603)
        fac101_var_1616 = fk1_var_1624 * fac01_var_1569(jl_var_1650, lay_var_1603)
        fac201_var_1617 = fk2_var_1625 * fac01_var_1569(jl_var_1650, lay_var_1603)
        fac011_var_1618 = fk0_var_1623 * fac11_var_1571(jl_var_1650, lay_var_1603)
        fac111_var_1619 = fk1_var_1624 * fac11_var_1571(jl_var_1650, lay_var_1603)
        fac211_var_1620 = fk2_var_1625 * fac11_var_1571(jl_var_1650, lay_var_1603)
      ELSE
        fac001_var_1615 = (1.0D0 - fs1_var_1634) * fac01_var_1569(jl_var_1650, lay_var_1603)
        fac011_var_1618 = (1.0D0 - fs1_var_1634) * fac11_var_1571(jl_var_1650, lay_var_1603)
        fac101_var_1616 = fs1_var_1634 * fac01_var_1569(jl_var_1650, lay_var_1603)
        fac111_var_1619 = fs1_var_1634 * fac11_var_1571(jl_var_1650, lay_var_1603)
        fac201_var_1617 = 0.0D0
        fac211_var_1620 = 0.0D0
      END IF
      IF (specparm_var_1633 .LT. 0.125D0) THEN
        tau_major_var_1628(1 : ng5) = speccomb_var_1593 * (fac000_var_1609 * absa_var_195(ind0_var_1596, 1 : 16) + fac100_var_1610 * absa_var_195(ind0_var_1596 + 1, 1 : 16) + fac200_var_1611 * absa_var_195(ind0_var_1596 + 2, 1 : 16) + fac010_var_1612 * absa_var_195(ind0_var_1596 + 9, 1 : 16) + fac110_var_1613 * absa_var_195(ind0_var_1596 + 10, 1 : 16) + fac210_var_1614 * absa_var_195(ind0_var_1596 + 11, 1 : 16))
      ELSE IF (specparm_var_1633 .GT. 0.875D0) THEN
        tau_major_var_1628(1 : ng5) = speccomb_var_1593 * (fac200_var_1611 * absa_var_195(ind0_var_1596 - 1, 1 : 16) + fac100_var_1610 * absa_var_195(ind0_var_1596, 1 : 16) + fac000_var_1609 * absa_var_195(ind0_var_1596 + 1, 1 : 16) + fac210_var_1614 * absa_var_195(ind0_var_1596 + 8, 1 : 16) + fac110_var_1613 * absa_var_195(ind0_var_1596 + 9, 1 : 16) + fac010_var_1612 * absa_var_195(ind0_var_1596 + 10, 1 : 16))
      ELSE
        tau_major_var_1628(1 : ng5) = speccomb_var_1593 * (fac000_var_1609 * absa_var_195(ind0_var_1596, 1 : 16) + fac100_var_1610 * absa_var_195(ind0_var_1596 + 1, 1 : 16) + fac010_var_1612 * absa_var_195(ind0_var_1596 + 9, 1 : 16) + fac110_var_1613 * absa_var_195(ind0_var_1596 + 10, 1 : 16))
      END IF
      IF (specparm1_var_1636 .LT. 0.125D0) THEN
        tau_major1_var_1629(1 : ng5) = speccomb1_var_1594 * (fac001_var_1615 * absa_var_195(ind1_var_1597, 1 : 16) + fac101_var_1616 * absa_var_195(ind1_var_1597 + 1, 1 : 16) + fac201_var_1617 * absa_var_195(ind1_var_1597 + 2, 1 : 16) + fac011_var_1618 * absa_var_195(ind1_var_1597 + 9, 1 : 16) + fac111_var_1619 * absa_var_195(ind1_var_1597 + 10, 1 : 16) + fac211_var_1620 * absa_var_195(ind1_var_1597 + 11, 1 : 16))
      ELSE IF (specparm1_var_1636 .GT. 0.875D0) THEN
        tau_major1_var_1629(1 : ng5) = speccomb1_var_1594 * (fac201_var_1617 * absa_var_195(ind1_var_1597 - 1, 1 : 16) + fac101_var_1616 * absa_var_195(ind1_var_1597, 1 : 16) + fac001_var_1615 * absa_var_195(ind1_var_1597 + 1, 1 : 16) + fac211_var_1620 * absa_var_195(ind1_var_1597 + 8, 1 : 16) + fac111_var_1619 * absa_var_195(ind1_var_1597 + 9, 1 : 16) + fac011_var_1618 * absa_var_195(ind1_var_1597 + 10, 1 : 16))
      ELSE
        tau_major1_var_1629(1 : ng5) = speccomb1_var_1594 * (fac001_var_1615 * absa_var_195(ind1_var_1597, 1 : 16) + fac101_var_1616 * absa_var_195(ind1_var_1597 + 1, 1 : 16) + fac011_var_1618 * absa_var_195(ind1_var_1597 + 9, 1 : 16) + fac111_var_1619 * absa_var_195(ind1_var_1597 + 10, 1 : 16))
      END IF
      DO ig_var_1601 = 1, 16
        tauself_var_1627 = selffac_var_1580(jl_var_1650, lay_var_1603) * (selfref_var_198(inds_var_1598, ig_var_1601) + selffrac_var_1581(jl_var_1650, lay_var_1603) * (selfref_var_198(inds_var_1598 + 1, ig_var_1601) - selfref_var_198(inds_var_1598, ig_var_1601)))
        taufor_var_1626 = forfac_var_1590(jl_var_1650, lay_var_1603) * (forref_var_199(indf_var_1599, ig_var_1601) + forfrac_var_1589(jl_var_1650, lay_var_1603) * (forref_var_199(indf_var_1599 + 1, ig_var_1601) - forref_var_199(indf_var_1599, ig_var_1601)))
        o3m1 = ka_mo3_var_197(jmo3, indm_var_1600, ig_var_1601) + fmo3 * (ka_mo3_var_197(jmo3 + 1, indm_var_1600, ig_var_1601) - ka_mo3_var_197(jmo3, indm_var_1600, ig_var_1601))
        o3m2 = ka_mo3_var_197(jmo3, indm_var_1600 + 1, ig_var_1601) + fmo3 * (ka_mo3_var_197(jmo3 + 1, indm_var_1600 + 1, ig_var_1601) - ka_mo3_var_197(jmo3, indm_var_1600 + 1, ig_var_1601))
        abso3_var_1630 = o3m1 + minorfrac_var_1591(jl_var_1650, lay_var_1603) * (o3m2 - o3m1)
        taug_var_1565(jl_var_1650, 52 + ig_var_1601, lay_var_1603) = tau_major_var_1628(ig_var_1601) + tau_major1_var_1629(ig_var_1601) + tauself_var_1627 + taufor_var_1626 + abso3_var_1630 * colo3_var_1578(jl_var_1650, lay_var_1603) + wx_var_1566(jl_var_1650, 1, lay_var_1603) * ccl4(ig_var_1601)
        fracs_var_1583(jl_var_1650, 52 + ig_var_1601, lay_var_1603) = fracrefa_var_193(ig_var_1601, jpl_var_1605) + fpl_var_1637 * (fracrefa_var_193(ig_var_1601, jpl_var_1605 + 1) - fracrefa_var_193(ig_var_1601, jpl_var_1605))
      END DO
    END DO
  END DO
  DO lay_var_1603 = laytrop_max_var_1641 + 1, klev_var_1564
    DO jl_var_1650 = kidia_var_1562, kfdia_var_1563
      speccomb_var_1593 = colo3_var_1578(jl_var_1650, lay_var_1603) + rat_o3co2_var_1586(jl_var_1650, lay_var_1603) * colco2_var_1577(jl_var_1650, lay_var_1603)
      specparm_var_1633 = MIN(colo3_var_1578(jl_var_1650, lay_var_1603) / speccomb_var_1593, oneminus_var_1575)
      specmult_var_1632 = 4.0D0 * (specparm_var_1633)
      js_var_1602 = 1 + INT(specmult_var_1632)
      fs_var_1631 = ((specmult_var_1632) - AINT((specmult_var_1632)))
      speccomb1_var_1594 = colo3_var_1578(jl_var_1650, lay_var_1603) + rat_o3co2_1_var_1587(jl_var_1650, lay_var_1603) * colco2_var_1577(jl_var_1650, lay_var_1603)
      specparm1_var_1636 = MIN(colo3_var_1578(jl_var_1650, lay_var_1603) / speccomb1_var_1594, oneminus_var_1575)
      specmult1_var_1635 = 4.0D0 * (specparm1_var_1636)
      js1_var_1604 = 1 + INT(specmult1_var_1635)
      fs1_var_1634 = ((specmult1_var_1635) - AINT((specmult1_var_1635)))
      fac000_var_1609 = (1.0D0 - fs_var_1631) * fac00_var_1568(jl_var_1650, lay_var_1603)
      fac010_var_1612 = (1.0D0 - fs_var_1631) * fac10_var_1570(jl_var_1650, lay_var_1603)
      fac100_var_1610 = fs_var_1631 * fac00_var_1568(jl_var_1650, lay_var_1603)
      fac110_var_1613 = fs_var_1631 * fac10_var_1570(jl_var_1650, lay_var_1603)
      fac001_var_1615 = (1.0D0 - fs1_var_1634) * fac01_var_1569(jl_var_1650, lay_var_1603)
      fac011_var_1618 = (1.0D0 - fs1_var_1634) * fac11_var_1571(jl_var_1650, lay_var_1603)
      fac101_var_1616 = fs1_var_1634 * fac01_var_1569(jl_var_1650, lay_var_1603)
      fac111_var_1619 = fs1_var_1634 * fac11_var_1571(jl_var_1650, lay_var_1603)
      speccomb_planck_var_1595 = colo3_var_1578(jl_var_1650, lay_var_1603) + refrat_planck_b_var_1607 * colco2_var_1577(jl_var_1650, lay_var_1603)
      specparm_planck_var_1639 = MIN(colo3_var_1578(jl_var_1650, lay_var_1603) / speccomb_planck_var_1595, oneminus_var_1575)
      specmult_planck_var_1638 = 4.0D0 * specparm_planck_var_1639
      jpl_var_1605 = 1 + INT(specmult_planck_var_1638)
      fpl_var_1637 = ((specmult_planck_var_1638) - AINT((specmult_planck_var_1638)))
      ind0_var_1596 = ((jp_var_1572(jl_var_1650, lay_var_1603) - 13) * 5 + (jt_var_1573(jl_var_1650, lay_var_1603) - 1)) * nspb_var_237(5) + js_var_1602
      ind1_var_1597 = ((jp_var_1572(jl_var_1650, lay_var_1603) - 12) * 5 + (jt1_var_1574(jl_var_1650, lay_var_1603) - 1)) * nspb_var_237(5) + js1_var_1604
      DO ig_var_1601 = 1, 16
        taug_var_1565(jl_var_1650, 52 + ig_var_1601, lay_var_1603) = speccomb_var_1593 * (fac000_var_1609 * absb_var_196(ind0_var_1596, ig_var_1601) + fac100_var_1610 * absb_var_196(ind0_var_1596 + 1, ig_var_1601) + fac010_var_1612 * absb_var_196(ind0_var_1596 + 5, ig_var_1601) + fac110_var_1613 * absb_var_196(ind0_var_1596 + 6, ig_var_1601)) + speccomb1_var_1594 * (fac001_var_1615 * absb_var_196(ind1_var_1597, ig_var_1601) + fac101_var_1616 * absb_var_196(ind1_var_1597 + 1, ig_var_1601) + fac011_var_1618 * absb_var_196(ind1_var_1597 + 5, ig_var_1601) + fac111_var_1619 * absb_var_196(ind1_var_1597 + 6, ig_var_1601)) + wx_var_1566(jl_var_1650, 1, lay_var_1603) * ccl4(ig_var_1601)
        fracs_var_1583(jl_var_1650, 52 + ig_var_1601, lay_var_1603) = fracrefb_var_194(ig_var_1601, jpl_var_1605) + fpl_var_1637 * (fracrefb_var_194(ig_var_1601, jpl_var_1605 + 1) - fracrefb_var_194(ig_var_1601, jpl_var_1605))
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1641 /= laytrop_min_var_1640) THEN
    DO lay_var_1603 = laytrop_min_var_1640 + 1, laytrop_max_var_1641
      ixc0_var_1647 = ixc_var_1642(lay_var_1603)
      DO ixp_var_1648 = 1, ixc0_var_1647
        jl_var_1650 = ixlow_var_1643(ixp_var_1648, lay_var_1603)
        speccomb_var_1593 = colh2o_var_1576(jl_var_1650, lay_var_1603) + rat_h2oco2_var_1584(jl_var_1650, lay_var_1603) * colco2_var_1577(jl_var_1650, lay_var_1603)
        specparm_var_1633 = MIN(colh2o_var_1576(jl_var_1650, lay_var_1603) / speccomb_var_1593, oneminus_var_1575)
        specmult_var_1632 = 8.0D0 * (specparm_var_1633)
        js_var_1602 = 1 + INT(specmult_var_1632)
        fs_var_1631 = ((specmult_var_1632) - AINT((specmult_var_1632)))
        speccomb1_var_1594 = colh2o_var_1576(jl_var_1650, lay_var_1603) + rat_h2oco2_1_var_1585(jl_var_1650, lay_var_1603) * colco2_var_1577(jl_var_1650, lay_var_1603)
        specparm1_var_1636 = MIN(colh2o_var_1576(jl_var_1650, lay_var_1603) / speccomb1_var_1594, oneminus_var_1575)
        specmult1_var_1635 = 8.0D0 * (specparm1_var_1636)
        js1_var_1604 = 1 + INT(specmult1_var_1635)
        fs1_var_1634 = ((specmult1_var_1635) - AINT((specmult1_var_1635)))
        speccomb_mo3 = colh2o_var_1576(jl_var_1650, lay_var_1603) + refrat_m_a_var_1608 * colco2_var_1577(jl_var_1650, lay_var_1603)
        specparm_mo3 = MIN(colh2o_var_1576(jl_var_1650, lay_var_1603) / speccomb_mo3, oneminus_var_1575)
        specmult_mo3 = 8.0D0 * specparm_mo3
        jmo3 = 1 + INT(specmult_mo3)
        fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))
        speccomb_planck_var_1595 = colh2o_var_1576(jl_var_1650, lay_var_1603) + refrat_planck_a_var_1606 * colco2_var_1577(jl_var_1650, lay_var_1603)
        specparm_planck_var_1639 = MIN(colh2o_var_1576(jl_var_1650, lay_var_1603) / speccomb_planck_var_1595, oneminus_var_1575)
        specmult_planck_var_1638 = 8.0D0 * specparm_planck_var_1639
        jpl_var_1605 = 1 + INT(specmult_planck_var_1638)
        fpl_var_1637 = ((specmult_planck_var_1638) - AINT((specmult_planck_var_1638)))
        ind0_var_1596 = ((jp_var_1572(jl_var_1650, lay_var_1603) - 1) * 5 + (jt_var_1573(jl_var_1650, lay_var_1603) - 1)) * nspa_var_236(5) + js_var_1602
        ind1_var_1597 = (jp_var_1572(jl_var_1650, lay_var_1603) * 5 + (jt1_var_1574(jl_var_1650, lay_var_1603) - 1)) * nspa_var_236(5) + js1_var_1604
        inds_var_1598 = indself_var_1582(jl_var_1650, lay_var_1603)
        indf_var_1599 = indfor_var_1588(jl_var_1650, lay_var_1603)
        indm_var_1600 = indminor_var_1592(jl_var_1650, lay_var_1603)
        IF (specparm_var_1633 .LT. 0.125D0) THEN
          p_var_1621 = fs_var_1631 - 1.0D0
          p4_var_1622 = p_var_1621 ** 4
          fk0_var_1623 = p4_var_1622
          fk1_var_1624 = 1.0D0 - p_var_1621 - 2.0D0 * p4_var_1622
          fk2_var_1625 = p_var_1621 + p4_var_1622
          fac000_var_1609 = fk0_var_1623 * fac00_var_1568(jl_var_1650, lay_var_1603)
          fac100_var_1610 = fk1_var_1624 * fac00_var_1568(jl_var_1650, lay_var_1603)
          fac200_var_1611 = fk2_var_1625 * fac00_var_1568(jl_var_1650, lay_var_1603)
          fac010_var_1612 = fk0_var_1623 * fac10_var_1570(jl_var_1650, lay_var_1603)
          fac110_var_1613 = fk1_var_1624 * fac10_var_1570(jl_var_1650, lay_var_1603)
          fac210_var_1614 = fk2_var_1625 * fac10_var_1570(jl_var_1650, lay_var_1603)
        ELSE IF (specparm_var_1633 .GT. 0.875D0) THEN
          p_var_1621 = - fs_var_1631
          p4_var_1622 = p_var_1621 ** 4
          fk0_var_1623 = p4_var_1622
          fk1_var_1624 = 1.0D0 - p_var_1621 - 2.0D0 * p4_var_1622
          fk2_var_1625 = p_var_1621 + p4_var_1622
          fac000_var_1609 = fk0_var_1623 * fac00_var_1568(jl_var_1650, lay_var_1603)
          fac100_var_1610 = fk1_var_1624 * fac00_var_1568(jl_var_1650, lay_var_1603)
          fac200_var_1611 = fk2_var_1625 * fac00_var_1568(jl_var_1650, lay_var_1603)
          fac010_var_1612 = fk0_var_1623 * fac10_var_1570(jl_var_1650, lay_var_1603)
          fac110_var_1613 = fk1_var_1624 * fac10_var_1570(jl_var_1650, lay_var_1603)
          fac210_var_1614 = fk2_var_1625 * fac10_var_1570(jl_var_1650, lay_var_1603)
        ELSE
          fac000_var_1609 = (1.0D0 - fs_var_1631) * fac00_var_1568(jl_var_1650, lay_var_1603)
          fac010_var_1612 = (1.0D0 - fs_var_1631) * fac10_var_1570(jl_var_1650, lay_var_1603)
          fac100_var_1610 = fs_var_1631 * fac00_var_1568(jl_var_1650, lay_var_1603)
          fac110_var_1613 = fs_var_1631 * fac10_var_1570(jl_var_1650, lay_var_1603)
          fac200_var_1611 = 0.0D0
          fac210_var_1614 = 0.0D0
        END IF
        IF (specparm1_var_1636 .LT. 0.125D0) THEN
          p_var_1621 = fs1_var_1634 - 1.0D0
          p4_var_1622 = p_var_1621 ** 4
          fk0_var_1623 = p4_var_1622
          fk1_var_1624 = 1.0D0 - p_var_1621 - 2.0D0 * p4_var_1622
          fk2_var_1625 = p_var_1621 + p4_var_1622
          fac001_var_1615 = fk0_var_1623 * fac01_var_1569(jl_var_1650, lay_var_1603)
          fac101_var_1616 = fk1_var_1624 * fac01_var_1569(jl_var_1650, lay_var_1603)
          fac201_var_1617 = fk2_var_1625 * fac01_var_1569(jl_var_1650, lay_var_1603)
          fac011_var_1618 = fk0_var_1623 * fac11_var_1571(jl_var_1650, lay_var_1603)
          fac111_var_1619 = fk1_var_1624 * fac11_var_1571(jl_var_1650, lay_var_1603)
          fac211_var_1620 = fk2_var_1625 * fac11_var_1571(jl_var_1650, lay_var_1603)
        ELSE IF (specparm1_var_1636 .GT. 0.875D0) THEN
          p_var_1621 = - fs1_var_1634
          p4_var_1622 = p_var_1621 ** 4
          fk0_var_1623 = p4_var_1622
          fk1_var_1624 = 1.0D0 - p_var_1621 - 2.0D0 * p4_var_1622
          fk2_var_1625 = p_var_1621 + p4_var_1622
          fac001_var_1615 = fk0_var_1623 * fac01_var_1569(jl_var_1650, lay_var_1603)
          fac101_var_1616 = fk1_var_1624 * fac01_var_1569(jl_var_1650, lay_var_1603)
          fac201_var_1617 = fk2_var_1625 * fac01_var_1569(jl_var_1650, lay_var_1603)
          fac011_var_1618 = fk0_var_1623 * fac11_var_1571(jl_var_1650, lay_var_1603)
          fac111_var_1619 = fk1_var_1624 * fac11_var_1571(jl_var_1650, lay_var_1603)
          fac211_var_1620 = fk2_var_1625 * fac11_var_1571(jl_var_1650, lay_var_1603)
        ELSE
          fac001_var_1615 = (1.0D0 - fs1_var_1634) * fac01_var_1569(jl_var_1650, lay_var_1603)
          fac011_var_1618 = (1.0D0 - fs1_var_1634) * fac11_var_1571(jl_var_1650, lay_var_1603)
          fac101_var_1616 = fs1_var_1634 * fac01_var_1569(jl_var_1650, lay_var_1603)
          fac111_var_1619 = fs1_var_1634 * fac11_var_1571(jl_var_1650, lay_var_1603)
          fac201_var_1617 = 0.0D0
          fac211_var_1620 = 0.0D0
        END IF
        IF (specparm_var_1633 .LT. 0.125D0) THEN
          tau_major_var_1628(1 : ng5) = speccomb_var_1593 * (fac000_var_1609 * absa_var_195(ind0_var_1596, 1 : 16) + fac100_var_1610 * absa_var_195(ind0_var_1596 + 1, 1 : 16) + fac200_var_1611 * absa_var_195(ind0_var_1596 + 2, 1 : 16) + fac010_var_1612 * absa_var_195(ind0_var_1596 + 9, 1 : 16) + fac110_var_1613 * absa_var_195(ind0_var_1596 + 10, 1 : 16) + fac210_var_1614 * absa_var_195(ind0_var_1596 + 11, 1 : 16))
        ELSE IF (specparm_var_1633 .GT. 0.875D0) THEN
          tau_major_var_1628(1 : ng5) = speccomb_var_1593 * (fac200_var_1611 * absa_var_195(ind0_var_1596 - 1, 1 : 16) + fac100_var_1610 * absa_var_195(ind0_var_1596, 1 : 16) + fac000_var_1609 * absa_var_195(ind0_var_1596 + 1, 1 : 16) + fac210_var_1614 * absa_var_195(ind0_var_1596 + 8, 1 : 16) + fac110_var_1613 * absa_var_195(ind0_var_1596 + 9, 1 : 16) + fac010_var_1612 * absa_var_195(ind0_var_1596 + 10, 1 : 16))
        ELSE
          tau_major_var_1628(1 : ng5) = speccomb_var_1593 * (fac000_var_1609 * absa_var_195(ind0_var_1596, 1 : 16) + fac100_var_1610 * absa_var_195(ind0_var_1596 + 1, 1 : 16) + fac010_var_1612 * absa_var_195(ind0_var_1596 + 9, 1 : 16) + fac110_var_1613 * absa_var_195(ind0_var_1596 + 10, 1 : 16))
        END IF
        IF (specparm1_var_1636 .LT. 0.125D0) THEN
          tau_major1_var_1629(1 : ng5) = speccomb1_var_1594 * (fac001_var_1615 * absa_var_195(ind1_var_1597, 1 : 16) + fac101_var_1616 * absa_var_195(ind1_var_1597 + 1, 1 : 16) + fac201_var_1617 * absa_var_195(ind1_var_1597 + 2, 1 : 16) + fac011_var_1618 * absa_var_195(ind1_var_1597 + 9, 1 : 16) + fac111_var_1619 * absa_var_195(ind1_var_1597 + 10, 1 : 16) + fac211_var_1620 * absa_var_195(ind1_var_1597 + 11, 1 : 16))
        ELSE IF (specparm1_var_1636 .GT. 0.875D0) THEN
          tau_major1_var_1629(1 : ng5) = speccomb1_var_1594 * (fac201_var_1617 * absa_var_195(ind1_var_1597 - 1, 1 : 16) + fac101_var_1616 * absa_var_195(ind1_var_1597, 1 : 16) + fac001_var_1615 * absa_var_195(ind1_var_1597 + 1, 1 : 16) + fac211_var_1620 * absa_var_195(ind1_var_1597 + 8, 1 : 16) + fac111_var_1619 * absa_var_195(ind1_var_1597 + 9, 1 : 16) + fac011_var_1618 * absa_var_195(ind1_var_1597 + 10, 1 : 16))
        ELSE
          tau_major1_var_1629(1 : ng5) = speccomb1_var_1594 * (fac001_var_1615 * absa_var_195(ind1_var_1597, 1 : 16) + fac101_var_1616 * absa_var_195(ind1_var_1597 + 1, 1 : 16) + fac011_var_1618 * absa_var_195(ind1_var_1597 + 9, 1 : 16) + fac111_var_1619 * absa_var_195(ind1_var_1597 + 10, 1 : 16))
        END IF
        DO ig_var_1601 = 1, 16
          tauself_var_1627 = selffac_var_1580(jl_var_1650, lay_var_1603) * (selfref_var_198(inds_var_1598, ig_var_1601) + selffrac_var_1581(jl_var_1650, lay_var_1603) * (selfref_var_198(inds_var_1598 + 1, ig_var_1601) - selfref_var_198(inds_var_1598, ig_var_1601)))
          taufor_var_1626 = forfac_var_1590(jl_var_1650, lay_var_1603) * (forref_var_199(indf_var_1599, ig_var_1601) + forfrac_var_1589(jl_var_1650, lay_var_1603) * (forref_var_199(indf_var_1599 + 1, ig_var_1601) - forref_var_199(indf_var_1599, ig_var_1601)))
          o3m1 = ka_mo3_var_197(jmo3, indm_var_1600, ig_var_1601) + fmo3 * (ka_mo3_var_197(jmo3 + 1, indm_var_1600, ig_var_1601) - ka_mo3_var_197(jmo3, indm_var_1600, ig_var_1601))
          o3m2 = ka_mo3_var_197(jmo3, indm_var_1600 + 1, ig_var_1601) + fmo3 * (ka_mo3_var_197(jmo3 + 1, indm_var_1600 + 1, ig_var_1601) - ka_mo3_var_197(jmo3, indm_var_1600 + 1, ig_var_1601))
          abso3_var_1630 = o3m1 + minorfrac_var_1591(jl_var_1650, lay_var_1603) * (o3m2 - o3m1)
          taug_var_1565(jl_var_1650, 52 + ig_var_1601, lay_var_1603) = tau_major_var_1628(ig_var_1601) + tau_major1_var_1629(ig_var_1601) + tauself_var_1627 + taufor_var_1626 + abso3_var_1630 * colo3_var_1578(jl_var_1650, lay_var_1603) + wx_var_1566(jl_var_1650, 1, lay_var_1603) * ccl4(ig_var_1601)
          fracs_var_1583(jl_var_1650, 52 + ig_var_1601, lay_var_1603) = fracrefa_var_193(ig_var_1601, jpl_var_1605) + fpl_var_1637 * (fracrefa_var_193(ig_var_1601, jpl_var_1605 + 1) - fracrefa_var_193(ig_var_1601, jpl_var_1605))
        END DO
      END DO
      ixc0_var_1647 = kfdia_var_1563 - kidia_var_1562 + 1 - ixc0_var_1647
      DO ixp_var_1648 = 1, ixc0_var_1647
        jl_var_1650 = ixhigh_var_1644(ixp_var_1648, lay_var_1603)
        speccomb_var_1593 = colo3_var_1578(jl_var_1650, lay_var_1603) + rat_o3co2_var_1586(jl_var_1650, lay_var_1603) * colco2_var_1577(jl_var_1650, lay_var_1603)
        specparm_var_1633 = MIN(colo3_var_1578(jl_var_1650, lay_var_1603) / speccomb_var_1593, oneminus_var_1575)
        specmult_var_1632 = 4.0D0 * (specparm_var_1633)
        js_var_1602 = 1 + INT(specmult_var_1632)
        fs_var_1631 = ((specmult_var_1632) - AINT((specmult_var_1632)))
        speccomb1_var_1594 = colo3_var_1578(jl_var_1650, lay_var_1603) + rat_o3co2_1_var_1587(jl_var_1650, lay_var_1603) * colco2_var_1577(jl_var_1650, lay_var_1603)
        specparm1_var_1636 = MIN(colo3_var_1578(jl_var_1650, lay_var_1603) / speccomb1_var_1594, oneminus_var_1575)
        specmult1_var_1635 = 4.0D0 * (specparm1_var_1636)
        js1_var_1604 = 1 + INT(specmult1_var_1635)
        fs1_var_1634 = ((specmult1_var_1635) - AINT((specmult1_var_1635)))
        fac000_var_1609 = (1.0D0 - fs_var_1631) * fac00_var_1568(jl_var_1650, lay_var_1603)
        fac010_var_1612 = (1.0D0 - fs_var_1631) * fac10_var_1570(jl_var_1650, lay_var_1603)
        fac100_var_1610 = fs_var_1631 * fac00_var_1568(jl_var_1650, lay_var_1603)
        fac110_var_1613 = fs_var_1631 * fac10_var_1570(jl_var_1650, lay_var_1603)
        fac001_var_1615 = (1.0D0 - fs1_var_1634) * fac01_var_1569(jl_var_1650, lay_var_1603)
        fac011_var_1618 = (1.0D0 - fs1_var_1634) * fac11_var_1571(jl_var_1650, lay_var_1603)
        fac101_var_1616 = fs1_var_1634 * fac01_var_1569(jl_var_1650, lay_var_1603)
        fac111_var_1619 = fs1_var_1634 * fac11_var_1571(jl_var_1650, lay_var_1603)
        speccomb_planck_var_1595 = colo3_var_1578(jl_var_1650, lay_var_1603) + refrat_planck_b_var_1607 * colco2_var_1577(jl_var_1650, lay_var_1603)
        specparm_planck_var_1639 = MIN(colo3_var_1578(jl_var_1650, lay_var_1603) / speccomb_planck_var_1595, oneminus_var_1575)
        specmult_planck_var_1638 = 4.0D0 * specparm_planck_var_1639
        jpl_var_1605 = 1 + INT(specmult_planck_var_1638)
        fpl_var_1637 = ((specmult_planck_var_1638) - AINT((specmult_planck_var_1638)))
        ind0_var_1596 = ((jp_var_1572(jl_var_1650, lay_var_1603) - 13) * 5 + (jt_var_1573(jl_var_1650, lay_var_1603) - 1)) * nspb_var_237(5) + js_var_1602
        ind1_var_1597 = ((jp_var_1572(jl_var_1650, lay_var_1603) - 12) * 5 + (jt1_var_1574(jl_var_1650, lay_var_1603) - 1)) * nspb_var_237(5) + js1_var_1604
        DO ig_var_1601 = 1, 16
          taug_var_1565(jl_var_1650, 52 + ig_var_1601, lay_var_1603) = speccomb_var_1593 * (fac000_var_1609 * absb_var_196(ind0_var_1596, ig_var_1601) + fac100_var_1610 * absb_var_196(ind0_var_1596 + 1, ig_var_1601) + fac010_var_1612 * absb_var_196(ind0_var_1596 + 5, ig_var_1601) + fac110_var_1613 * absb_var_196(ind0_var_1596 + 6, ig_var_1601)) + speccomb1_var_1594 * (fac001_var_1615 * absb_var_196(ind1_var_1597, ig_var_1601) + fac101_var_1616 * absb_var_196(ind1_var_1597 + 1, ig_var_1601) + fac011_var_1618 * absb_var_196(ind1_var_1597 + 5, ig_var_1601) + fac111_var_1619 * absb_var_196(ind1_var_1597 + 6, ig_var_1601)) + wx_var_1566(jl_var_1650, 1, lay_var_1603) * ccl4(ig_var_1601)
          fracs_var_1583(jl_var_1650, 52 + ig_var_1601, lay_var_1603) = fracrefb_var_194(ig_var_1601, jpl_var_1605) + fpl_var_1637 * (fracrefb_var_194(ig_var_1601, jpl_var_1605 + 1) - fracrefb_var_194(ig_var_1601, jpl_var_1605))
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol5
SUBROUTINE rrtm_taumol4(kidia_var_1651, kfdia_var_1652, klev_var_1653, taug_var_1654, p_tauaerl_var_1655, fac00_var_1656, fac01_var_1657, fac10_var_1658, fac11_var_1659, forfac_var_1677, forfrac_var_1678, indfor_var_1676, jp_var_1660, jt_var_1661, jt1_var_1662, oneminus_var_1663, colh2o_var_1664, colco2_var_1665, colo3_var_1666, laytrop_var_1667, selffac_var_1668, selffrac_var_1669, indself_var_1670, fracs_var_1671, rat_h2oco2_var_1672, rat_h2oco2_1_var_1673, rat_o3co2_var_1674, rat_o3co2_1_var_1675)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrtm, ONLY: ng4
  USE yoerrta4, ONLY: absa_var_189, absb_var_190, forref_var_192, fracrefa_var_187, fracrefb_var_188, selfref_var_191
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1651
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1652
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1653
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1654(kidia_var_1651 : kfdia_var_1652, 140, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1655(kidia_var_1651 : kfdia_var_1652, klev_var_1653, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1656(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1657(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1658(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1659(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1660(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1661(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1662(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1663
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1664(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1665(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1666(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1667(kidia_var_1651 : kfdia_var_1652)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1668(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1669(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1670(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1671(kidia_var_1651 : kfdia_var_1652, 140, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1672(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1673(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_var_1674(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_1_var_1675(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1676(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1677(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1678(kidia_var_1651 : kfdia_var_1652, klev_var_1653)
  REAL(KIND = 8) :: speccomb_var_1679, speccomb1_var_1680, speccomb_planck_var_1681
  INTEGER(KIND = 4) :: ind0_var_1682, ind1_var_1683, inds_var_1684, indf_var_1685
  INTEGER(KIND = 4) :: ig_var_1686, js_var_1687, lay_var_1688, js1_var_1689, jpl_var_1690
  REAL(KIND = 8) :: refrat_planck_a_var_1691, refrat_planck_b_var_1692
  REAL(KIND = 8) :: fac000_var_1693, fac100_var_1694, fac200_var_1695, fac010_var_1696, fac110_var_1697, fac210_var_1698, fac001_var_1699, fac101_var_1700, fac201_var_1701, fac011_var_1702, fac111_var_1703, fac211_var_1704
  REAL(KIND = 8) :: p_var_1705, p4_var_1706, fk0_var_1707, fk1_var_1708, fk2_var_1709
  REAL(KIND = 8) :: taufor_var_1710, tauself_var_1711, tau_major_var_1712(14), tau_major1_var_1713(14)
  REAL(KIND = 8) :: fs_var_1714, specmult_var_1715, specparm_var_1716, fs1_var_1717, specmult1_var_1718, specparm1_var_1719, fpl_var_1720, specmult_planck_var_1721, specparm_planck_var_1722
  INTEGER(KIND = 4) :: laytrop_min_var_1723, laytrop_max_var_1724
  INTEGER(KIND = 4) :: ixc_var_1725(klev_var_1653), ixlow_var_1726(kfdia_var_1652, klev_var_1653), ixhigh_var_1727(kfdia_var_1652, klev_var_1653)
  INTEGER(KIND = 4) :: ich_var_1728, icl_var_1729, ixc0_var_1730, ixp_var_1731, jc_var_1732, jl_var_1733
  laytrop_min_var_1723 = MINVAL(laytrop_var_1667)
  laytrop_max_var_1724 = MAXVAL(laytrop_var_1667)
  ixlow_var_1726 = 0
  ixhigh_var_1727 = 0
  ixc_var_1725 = 0
  DO lay_var_1688 = laytrop_min_var_1723 + 1, laytrop_max_var_1724
    icl_var_1729 = 0
    ich_var_1728 = 0
    DO jc_var_1732 = kidia_var_1651, kfdia_var_1652
      IF (lay_var_1688 <= laytrop_var_1667(jc_var_1732)) THEN
        icl_var_1729 = icl_var_1729 + 1
        ixlow_var_1726(icl_var_1729, lay_var_1688) = jc_var_1732
      ELSE
        ich_var_1728 = ich_var_1728 + 1
        ixhigh_var_1727(ich_var_1728, lay_var_1688) = jc_var_1732
      END IF
    END DO
    ixc_var_1725(lay_var_1688) = icl_var_1729
  END DO
  refrat_planck_a_var_1691 = chi_mls(1, 11) / chi_mls(2, 11)
  refrat_planck_b_var_1692 = chi_mls(3, 13) / chi_mls(2, 13)
  DO lay_var_1688 = 1, laytrop_min_var_1723
    DO jl_var_1733 = kidia_var_1651, kfdia_var_1652
      speccomb_var_1679 = colh2o_var_1664(jl_var_1733, lay_var_1688) + rat_h2oco2_var_1672(jl_var_1733, lay_var_1688) * colco2_var_1665(jl_var_1733, lay_var_1688)
      specparm_var_1716 = MIN(colh2o_var_1664(jl_var_1733, lay_var_1688) / speccomb_var_1679, oneminus_var_1663)
      specmult_var_1715 = 8.0D0 * (specparm_var_1716)
      js_var_1687 = 1 + INT(specmult_var_1715)
      fs_var_1714 = ((specmult_var_1715) - AINT((specmult_var_1715)))
      speccomb1_var_1680 = colh2o_var_1664(jl_var_1733, lay_var_1688) + rat_h2oco2_1_var_1673(jl_var_1733, lay_var_1688) * colco2_var_1665(jl_var_1733, lay_var_1688)
      specparm1_var_1719 = MIN(colh2o_var_1664(jl_var_1733, lay_var_1688) / speccomb1_var_1680, oneminus_var_1663)
      specmult1_var_1718 = 8.0D0 * (specparm1_var_1719)
      js1_var_1689 = 1 + INT(specmult1_var_1718)
      fs1_var_1717 = ((specmult1_var_1718) - AINT((specmult1_var_1718)))
      speccomb_planck_var_1681 = colh2o_var_1664(jl_var_1733, lay_var_1688) + refrat_planck_a_var_1691 * colco2_var_1665(jl_var_1733, lay_var_1688)
      specparm_planck_var_1722 = MIN(colh2o_var_1664(jl_var_1733, lay_var_1688) / speccomb_planck_var_1681, oneminus_var_1663)
      specmult_planck_var_1721 = 8.0D0 * specparm_planck_var_1722
      jpl_var_1690 = 1 + INT(specmult_planck_var_1721)
      fpl_var_1720 = ((specmult_planck_var_1721) - AINT((specmult_planck_var_1721)))
      ind0_var_1682 = ((jp_var_1660(jl_var_1733, lay_var_1688) - 1) * 5 + (jt_var_1661(jl_var_1733, lay_var_1688) - 1)) * nspa_var_236(4) + js_var_1687
      ind1_var_1683 = (jp_var_1660(jl_var_1733, lay_var_1688) * 5 + (jt1_var_1662(jl_var_1733, lay_var_1688) - 1)) * nspa_var_236(4) + js1_var_1689
      inds_var_1684 = indself_var_1670(jl_var_1733, lay_var_1688)
      indf_var_1685 = indfor_var_1676(jl_var_1733, lay_var_1688)
      IF (specparm_var_1716 .LT. 0.125D0) THEN
        p_var_1705 = fs_var_1714 - 1.0D0
        p4_var_1706 = p_var_1705 ** 4
        fk0_var_1707 = p4_var_1706
        fk1_var_1708 = 1.0D0 - p_var_1705 - 2.0D0 * p4_var_1706
        fk2_var_1709 = p_var_1705 + p4_var_1706
        fac000_var_1693 = fk0_var_1707 * fac00_var_1656(jl_var_1733, lay_var_1688)
        fac100_var_1694 = fk1_var_1708 * fac00_var_1656(jl_var_1733, lay_var_1688)
        fac200_var_1695 = fk2_var_1709 * fac00_var_1656(jl_var_1733, lay_var_1688)
        fac010_var_1696 = fk0_var_1707 * fac10_var_1658(jl_var_1733, lay_var_1688)
        fac110_var_1697 = fk1_var_1708 * fac10_var_1658(jl_var_1733, lay_var_1688)
        fac210_var_1698 = fk2_var_1709 * fac10_var_1658(jl_var_1733, lay_var_1688)
      ELSE IF (specparm_var_1716 .GT. 0.875D0) THEN
        p_var_1705 = - fs_var_1714
        p4_var_1706 = p_var_1705 ** 4
        fk0_var_1707 = p4_var_1706
        fk1_var_1708 = 1.0D0 - p_var_1705 - 2.0D0 * p4_var_1706
        fk2_var_1709 = p_var_1705 + p4_var_1706
        fac000_var_1693 = fk0_var_1707 * fac00_var_1656(jl_var_1733, lay_var_1688)
        fac100_var_1694 = fk1_var_1708 * fac00_var_1656(jl_var_1733, lay_var_1688)
        fac200_var_1695 = fk2_var_1709 * fac00_var_1656(jl_var_1733, lay_var_1688)
        fac010_var_1696 = fk0_var_1707 * fac10_var_1658(jl_var_1733, lay_var_1688)
        fac110_var_1697 = fk1_var_1708 * fac10_var_1658(jl_var_1733, lay_var_1688)
        fac210_var_1698 = fk2_var_1709 * fac10_var_1658(jl_var_1733, lay_var_1688)
      ELSE
        fac000_var_1693 = (1.0D0 - fs_var_1714) * fac00_var_1656(jl_var_1733, lay_var_1688)
        fac010_var_1696 = (1.0D0 - fs_var_1714) * fac10_var_1658(jl_var_1733, lay_var_1688)
        fac100_var_1694 = fs_var_1714 * fac00_var_1656(jl_var_1733, lay_var_1688)
        fac110_var_1697 = fs_var_1714 * fac10_var_1658(jl_var_1733, lay_var_1688)
        fac200_var_1695 = 0.0D0
        fac210_var_1698 = 0.0D0
      END IF
      IF (specparm1_var_1719 .LT. 0.125D0) THEN
        p_var_1705 = fs1_var_1717 - 1.0D0
        p4_var_1706 = p_var_1705 ** 4
        fk0_var_1707 = p4_var_1706
        fk1_var_1708 = 1.0D0 - p_var_1705 - 2.0D0 * p4_var_1706
        fk2_var_1709 = p_var_1705 + p4_var_1706
        fac001_var_1699 = fk0_var_1707 * fac01_var_1657(jl_var_1733, lay_var_1688)
        fac101_var_1700 = fk1_var_1708 * fac01_var_1657(jl_var_1733, lay_var_1688)
        fac201_var_1701 = fk2_var_1709 * fac01_var_1657(jl_var_1733, lay_var_1688)
        fac011_var_1702 = fk0_var_1707 * fac11_var_1659(jl_var_1733, lay_var_1688)
        fac111_var_1703 = fk1_var_1708 * fac11_var_1659(jl_var_1733, lay_var_1688)
        fac211_var_1704 = fk2_var_1709 * fac11_var_1659(jl_var_1733, lay_var_1688)
      ELSE IF (specparm1_var_1719 .GT. 0.875D0) THEN
        p_var_1705 = - fs1_var_1717
        p4_var_1706 = p_var_1705 ** 4
        fk0_var_1707 = p4_var_1706
        fk1_var_1708 = 1.0D0 - p_var_1705 - 2.0D0 * p4_var_1706
        fk2_var_1709 = p_var_1705 + p4_var_1706
        fac001_var_1699 = fk0_var_1707 * fac01_var_1657(jl_var_1733, lay_var_1688)
        fac101_var_1700 = fk1_var_1708 * fac01_var_1657(jl_var_1733, lay_var_1688)
        fac201_var_1701 = fk2_var_1709 * fac01_var_1657(jl_var_1733, lay_var_1688)
        fac011_var_1702 = fk0_var_1707 * fac11_var_1659(jl_var_1733, lay_var_1688)
        fac111_var_1703 = fk1_var_1708 * fac11_var_1659(jl_var_1733, lay_var_1688)
        fac211_var_1704 = fk2_var_1709 * fac11_var_1659(jl_var_1733, lay_var_1688)
      ELSE
        fac001_var_1699 = (1.0D0 - fs1_var_1717) * fac01_var_1657(jl_var_1733, lay_var_1688)
        fac011_var_1702 = (1.0D0 - fs1_var_1717) * fac11_var_1659(jl_var_1733, lay_var_1688)
        fac101_var_1700 = fs1_var_1717 * fac01_var_1657(jl_var_1733, lay_var_1688)
        fac111_var_1703 = fs1_var_1717 * fac11_var_1659(jl_var_1733, lay_var_1688)
        fac201_var_1701 = 0.0D0
        fac211_var_1704 = 0.0D0
      END IF
      IF (specparm_var_1716 .LT. 0.125D0) THEN
        tau_major_var_1712(1 : ng4) = speccomb_var_1679 * (fac000_var_1693 * absa_var_189(ind0_var_1682, 1 : 14) + fac100_var_1694 * absa_var_189(ind0_var_1682 + 1, 1 : 14) + fac200_var_1695 * absa_var_189(ind0_var_1682 + 2, 1 : 14) + fac010_var_1696 * absa_var_189(ind0_var_1682 + 9, 1 : 14) + fac110_var_1697 * absa_var_189(ind0_var_1682 + 10, 1 : 14) + fac210_var_1698 * absa_var_189(ind0_var_1682 + 11, 1 : 14))
      ELSE IF (specparm_var_1716 .GT. 0.875D0) THEN
        tau_major_var_1712(1 : ng4) = speccomb_var_1679 * (fac200_var_1695 * absa_var_189(ind0_var_1682 - 1, 1 : 14) + fac100_var_1694 * absa_var_189(ind0_var_1682, 1 : 14) + fac000_var_1693 * absa_var_189(ind0_var_1682 + 1, 1 : 14) + fac210_var_1698 * absa_var_189(ind0_var_1682 + 8, 1 : 14) + fac110_var_1697 * absa_var_189(ind0_var_1682 + 9, 1 : 14) + fac010_var_1696 * absa_var_189(ind0_var_1682 + 10, 1 : 14))
      ELSE
        tau_major_var_1712(1 : ng4) = speccomb_var_1679 * (fac000_var_1693 * absa_var_189(ind0_var_1682, 1 : 14) + fac100_var_1694 * absa_var_189(ind0_var_1682 + 1, 1 : 14) + fac010_var_1696 * absa_var_189(ind0_var_1682 + 9, 1 : 14) + fac110_var_1697 * absa_var_189(ind0_var_1682 + 10, 1 : 14))
      END IF
      IF (specparm1_var_1719 .LT. 0.125D0) THEN
        tau_major1_var_1713(1 : ng4) = speccomb1_var_1680 * (fac001_var_1699 * absa_var_189(ind1_var_1683, 1 : 14) + fac101_var_1700 * absa_var_189(ind1_var_1683 + 1, 1 : 14) + fac201_var_1701 * absa_var_189(ind1_var_1683 + 2, 1 : 14) + fac011_var_1702 * absa_var_189(ind1_var_1683 + 9, 1 : 14) + fac111_var_1703 * absa_var_189(ind1_var_1683 + 10, 1 : 14) + fac211_var_1704 * absa_var_189(ind1_var_1683 + 11, 1 : 14))
      ELSE IF (specparm1_var_1719 .GT. 0.875D0) THEN
        tau_major1_var_1713(1 : ng4) = speccomb1_var_1680 * (fac201_var_1701 * absa_var_189(ind1_var_1683 - 1, 1 : 14) + fac101_var_1700 * absa_var_189(ind1_var_1683, 1 : 14) + fac001_var_1699 * absa_var_189(ind1_var_1683 + 1, 1 : 14) + fac211_var_1704 * absa_var_189(ind1_var_1683 + 8, 1 : 14) + fac111_var_1703 * absa_var_189(ind1_var_1683 + 9, 1 : 14) + fac011_var_1702 * absa_var_189(ind1_var_1683 + 10, 1 : 14))
      ELSE
        tau_major1_var_1713(1 : ng4) = speccomb1_var_1680 * (fac001_var_1699 * absa_var_189(ind1_var_1683, 1 : 14) + fac101_var_1700 * absa_var_189(ind1_var_1683 + 1, 1 : 14) + fac011_var_1702 * absa_var_189(ind1_var_1683 + 9, 1 : 14) + fac111_var_1703 * absa_var_189(ind1_var_1683 + 10, 1 : 14))
      END IF
      DO ig_var_1686 = 1, 14
        tauself_var_1711 = selffac_var_1668(jl_var_1733, lay_var_1688) * (selfref_var_191(inds_var_1684, ig_var_1686) + selffrac_var_1669(jl_var_1733, lay_var_1688) * (selfref_var_191(inds_var_1684 + 1, ig_var_1686) - selfref_var_191(inds_var_1684, ig_var_1686)))
        taufor_var_1710 = forfac_var_1677(jl_var_1733, lay_var_1688) * (forref_var_192(indf_var_1685, ig_var_1686) + forfrac_var_1678(jl_var_1733, lay_var_1688) * (forref_var_192(indf_var_1685 + 1, ig_var_1686) - forref_var_192(indf_var_1685, ig_var_1686)))
        taug_var_1654(jl_var_1733, 38 + ig_var_1686, lay_var_1688) = tau_major_var_1712(ig_var_1686) + tau_major1_var_1713(ig_var_1686) + tauself_var_1711 + taufor_var_1710
        fracs_var_1671(jl_var_1733, 38 + ig_var_1686, lay_var_1688) = fracrefa_var_187(ig_var_1686, jpl_var_1690) + fpl_var_1720 * (fracrefa_var_187(ig_var_1686, jpl_var_1690 + 1) - fracrefa_var_187(ig_var_1686, jpl_var_1690))
      END DO
    END DO
  END DO
  DO lay_var_1688 = laytrop_max_var_1724 + 1, klev_var_1653
    DO jl_var_1733 = kidia_var_1651, kfdia_var_1652
      speccomb_var_1679 = colo3_var_1666(jl_var_1733, lay_var_1688) + rat_o3co2_var_1674(jl_var_1733, lay_var_1688) * colco2_var_1665(jl_var_1733, lay_var_1688)
      specparm_var_1716 = MIN(colo3_var_1666(jl_var_1733, lay_var_1688) / speccomb_var_1679, oneminus_var_1663)
      specmult_var_1715 = 4.0D0 * (specparm_var_1716)
      js_var_1687 = 1 + INT(specmult_var_1715)
      fs_var_1714 = ((specmult_var_1715) - AINT((specmult_var_1715)))
      speccomb1_var_1680 = colo3_var_1666(jl_var_1733, lay_var_1688) + rat_o3co2_1_var_1675(jl_var_1733, lay_var_1688) * colco2_var_1665(jl_var_1733, lay_var_1688)
      specparm1_var_1719 = MIN(colo3_var_1666(jl_var_1733, lay_var_1688) / speccomb1_var_1680, oneminus_var_1663)
      specmult1_var_1718 = 4.0D0 * (specparm1_var_1719)
      js1_var_1689 = 1 + INT(specmult1_var_1718)
      fs1_var_1717 = ((specmult1_var_1718) - AINT((specmult1_var_1718)))
      fac000_var_1693 = (1.0D0 - fs_var_1714) * fac00_var_1656(jl_var_1733, lay_var_1688)
      fac010_var_1696 = (1.0D0 - fs_var_1714) * fac10_var_1658(jl_var_1733, lay_var_1688)
      fac100_var_1694 = fs_var_1714 * fac00_var_1656(jl_var_1733, lay_var_1688)
      fac110_var_1697 = fs_var_1714 * fac10_var_1658(jl_var_1733, lay_var_1688)
      fac001_var_1699 = (1.0D0 - fs1_var_1717) * fac01_var_1657(jl_var_1733, lay_var_1688)
      fac011_var_1702 = (1.0D0 - fs1_var_1717) * fac11_var_1659(jl_var_1733, lay_var_1688)
      fac101_var_1700 = fs1_var_1717 * fac01_var_1657(jl_var_1733, lay_var_1688)
      fac111_var_1703 = fs1_var_1717 * fac11_var_1659(jl_var_1733, lay_var_1688)
      speccomb_planck_var_1681 = colo3_var_1666(jl_var_1733, lay_var_1688) + refrat_planck_b_var_1692 * colco2_var_1665(jl_var_1733, lay_var_1688)
      specparm_planck_var_1722 = MIN(colo3_var_1666(jl_var_1733, lay_var_1688) / speccomb_planck_var_1681, oneminus_var_1663)
      specmult_planck_var_1721 = 4.0D0 * specparm_planck_var_1722
      jpl_var_1690 = 1 + INT(specmult_planck_var_1721)
      fpl_var_1720 = ((specmult_planck_var_1721) - AINT((specmult_planck_var_1721)))
      ind0_var_1682 = ((jp_var_1660(jl_var_1733, lay_var_1688) - 13) * 5 + (jt_var_1661(jl_var_1733, lay_var_1688) - 1)) * nspb_var_237(4) + js_var_1687
      ind1_var_1683 = ((jp_var_1660(jl_var_1733, lay_var_1688) - 12) * 5 + (jt1_var_1662(jl_var_1733, lay_var_1688) - 1)) * nspb_var_237(4) + js1_var_1689
      DO ig_var_1686 = 1, 14
        taug_var_1654(jl_var_1733, 38 + ig_var_1686, lay_var_1688) = speccomb_var_1679 * (fac000_var_1693 * absb_var_190(ind0_var_1682, ig_var_1686) + fac100_var_1694 * absb_var_190(ind0_var_1682 + 1, ig_var_1686) + fac010_var_1696 * absb_var_190(ind0_var_1682 + 5, ig_var_1686) + fac110_var_1697 * absb_var_190(ind0_var_1682 + 6, ig_var_1686)) + speccomb1_var_1680 * (fac001_var_1699 * absb_var_190(ind1_var_1683, ig_var_1686) + fac101_var_1700 * absb_var_190(ind1_var_1683 + 1, ig_var_1686) + fac011_var_1702 * absb_var_190(ind1_var_1683 + 5, ig_var_1686) + fac111_var_1703 * absb_var_190(ind1_var_1683 + 6, ig_var_1686))
        fracs_var_1671(jl_var_1733, 38 + ig_var_1686, lay_var_1688) = fracrefb_var_188(ig_var_1686, jpl_var_1690) + fpl_var_1720 * (fracrefb_var_188(ig_var_1686, jpl_var_1690 + 1) - fracrefb_var_188(ig_var_1686, jpl_var_1690))
      END DO
    END DO
  END DO
  DO lay_var_1688 = laytrop_max_var_1724 + 1, klev_var_1653
    DO jl_var_1733 = kidia_var_1651, kfdia_var_1652
      taug_var_1654(jl_var_1733, 46, lay_var_1688) = taug_var_1654(jl_var_1733, 46, lay_var_1688) * 0.92D0
      taug_var_1654(jl_var_1733, 47, lay_var_1688) = taug_var_1654(jl_var_1733, 47, lay_var_1688) * 0.88D0
      taug_var_1654(jl_var_1733, 48, lay_var_1688) = taug_var_1654(jl_var_1733, 48, lay_var_1688) * 1.07D0
      taug_var_1654(jl_var_1733, 49, lay_var_1688) = taug_var_1654(jl_var_1733, 49, lay_var_1688) * 1.1D0
      taug_var_1654(jl_var_1733, 50, lay_var_1688) = taug_var_1654(jl_var_1733, 50, lay_var_1688) * 0.99D0
      taug_var_1654(jl_var_1733, 51, lay_var_1688) = taug_var_1654(jl_var_1733, 51, lay_var_1688) * 0.88D0
      taug_var_1654(jl_var_1733, 52, lay_var_1688) = taug_var_1654(jl_var_1733, 52, lay_var_1688) * 0.943D0
    END DO
  END DO
  IF (laytrop_max_var_1724 /= laytrop_min_var_1723) THEN
    DO lay_var_1688 = laytrop_min_var_1723 + 1, laytrop_max_var_1724
      ixc0_var_1730 = ixc_var_1725(lay_var_1688)
      DO ixp_var_1731 = 1, ixc0_var_1730
        jl_var_1733 = ixlow_var_1726(ixp_var_1731, lay_var_1688)
        speccomb_var_1679 = colh2o_var_1664(jl_var_1733, lay_var_1688) + rat_h2oco2_var_1672(jl_var_1733, lay_var_1688) * colco2_var_1665(jl_var_1733, lay_var_1688)
        specparm_var_1716 = MIN(colh2o_var_1664(jl_var_1733, lay_var_1688) / speccomb_var_1679, oneminus_var_1663)
        specmult_var_1715 = 8.0D0 * (specparm_var_1716)
        js_var_1687 = 1 + INT(specmult_var_1715)
        fs_var_1714 = ((specmult_var_1715) - AINT((specmult_var_1715)))
        speccomb1_var_1680 = colh2o_var_1664(jl_var_1733, lay_var_1688) + rat_h2oco2_1_var_1673(jl_var_1733, lay_var_1688) * colco2_var_1665(jl_var_1733, lay_var_1688)
        specparm1_var_1719 = MIN(colh2o_var_1664(jl_var_1733, lay_var_1688) / speccomb1_var_1680, oneminus_var_1663)
        specmult1_var_1718 = 8.0D0 * (specparm1_var_1719)
        js1_var_1689 = 1 + INT(specmult1_var_1718)
        fs1_var_1717 = ((specmult1_var_1718) - AINT((specmult1_var_1718)))
        speccomb_planck_var_1681 = colh2o_var_1664(jl_var_1733, lay_var_1688) + refrat_planck_a_var_1691 * colco2_var_1665(jl_var_1733, lay_var_1688)
        specparm_planck_var_1722 = MIN(colh2o_var_1664(jl_var_1733, lay_var_1688) / speccomb_planck_var_1681, oneminus_var_1663)
        specmult_planck_var_1721 = 8.0D0 * specparm_planck_var_1722
        jpl_var_1690 = 1 + INT(specmult_planck_var_1721)
        fpl_var_1720 = ((specmult_planck_var_1721) - AINT((specmult_planck_var_1721)))
        ind0_var_1682 = ((jp_var_1660(jl_var_1733, lay_var_1688) - 1) * 5 + (jt_var_1661(jl_var_1733, lay_var_1688) - 1)) * nspa_var_236(4) + js_var_1687
        ind1_var_1683 = (jp_var_1660(jl_var_1733, lay_var_1688) * 5 + (jt1_var_1662(jl_var_1733, lay_var_1688) - 1)) * nspa_var_236(4) + js1_var_1689
        inds_var_1684 = indself_var_1670(jl_var_1733, lay_var_1688)
        indf_var_1685 = indfor_var_1676(jl_var_1733, lay_var_1688)
        IF (specparm_var_1716 .LT. 0.125D0) THEN
          p_var_1705 = fs_var_1714 - 1.0D0
          p4_var_1706 = p_var_1705 ** 4
          fk0_var_1707 = p4_var_1706
          fk1_var_1708 = 1.0D0 - p_var_1705 - 2.0D0 * p4_var_1706
          fk2_var_1709 = p_var_1705 + p4_var_1706
          fac000_var_1693 = fk0_var_1707 * fac00_var_1656(jl_var_1733, lay_var_1688)
          fac100_var_1694 = fk1_var_1708 * fac00_var_1656(jl_var_1733, lay_var_1688)
          fac200_var_1695 = fk2_var_1709 * fac00_var_1656(jl_var_1733, lay_var_1688)
          fac010_var_1696 = fk0_var_1707 * fac10_var_1658(jl_var_1733, lay_var_1688)
          fac110_var_1697 = fk1_var_1708 * fac10_var_1658(jl_var_1733, lay_var_1688)
          fac210_var_1698 = fk2_var_1709 * fac10_var_1658(jl_var_1733, lay_var_1688)
        ELSE IF (specparm_var_1716 .GT. 0.875D0) THEN
          p_var_1705 = - fs_var_1714
          p4_var_1706 = p_var_1705 ** 4
          fk0_var_1707 = p4_var_1706
          fk1_var_1708 = 1.0D0 - p_var_1705 - 2.0D0 * p4_var_1706
          fk2_var_1709 = p_var_1705 + p4_var_1706
          fac000_var_1693 = fk0_var_1707 * fac00_var_1656(jl_var_1733, lay_var_1688)
          fac100_var_1694 = fk1_var_1708 * fac00_var_1656(jl_var_1733, lay_var_1688)
          fac200_var_1695 = fk2_var_1709 * fac00_var_1656(jl_var_1733, lay_var_1688)
          fac010_var_1696 = fk0_var_1707 * fac10_var_1658(jl_var_1733, lay_var_1688)
          fac110_var_1697 = fk1_var_1708 * fac10_var_1658(jl_var_1733, lay_var_1688)
          fac210_var_1698 = fk2_var_1709 * fac10_var_1658(jl_var_1733, lay_var_1688)
        ELSE
          fac000_var_1693 = (1.0D0 - fs_var_1714) * fac00_var_1656(jl_var_1733, lay_var_1688)
          fac010_var_1696 = (1.0D0 - fs_var_1714) * fac10_var_1658(jl_var_1733, lay_var_1688)
          fac100_var_1694 = fs_var_1714 * fac00_var_1656(jl_var_1733, lay_var_1688)
          fac110_var_1697 = fs_var_1714 * fac10_var_1658(jl_var_1733, lay_var_1688)
          fac200_var_1695 = 0.0D0
          fac210_var_1698 = 0.0D0
        END IF
        IF (specparm1_var_1719 .LT. 0.125D0) THEN
          p_var_1705 = fs1_var_1717 - 1.0D0
          p4_var_1706 = p_var_1705 ** 4
          fk0_var_1707 = p4_var_1706
          fk1_var_1708 = 1.0D0 - p_var_1705 - 2.0D0 * p4_var_1706
          fk2_var_1709 = p_var_1705 + p4_var_1706
          fac001_var_1699 = fk0_var_1707 * fac01_var_1657(jl_var_1733, lay_var_1688)
          fac101_var_1700 = fk1_var_1708 * fac01_var_1657(jl_var_1733, lay_var_1688)
          fac201_var_1701 = fk2_var_1709 * fac01_var_1657(jl_var_1733, lay_var_1688)
          fac011_var_1702 = fk0_var_1707 * fac11_var_1659(jl_var_1733, lay_var_1688)
          fac111_var_1703 = fk1_var_1708 * fac11_var_1659(jl_var_1733, lay_var_1688)
          fac211_var_1704 = fk2_var_1709 * fac11_var_1659(jl_var_1733, lay_var_1688)
        ELSE IF (specparm1_var_1719 .GT. 0.875D0) THEN
          p_var_1705 = - fs1_var_1717
          p4_var_1706 = p_var_1705 ** 4
          fk0_var_1707 = p4_var_1706
          fk1_var_1708 = 1.0D0 - p_var_1705 - 2.0D0 * p4_var_1706
          fk2_var_1709 = p_var_1705 + p4_var_1706
          fac001_var_1699 = fk0_var_1707 * fac01_var_1657(jl_var_1733, lay_var_1688)
          fac101_var_1700 = fk1_var_1708 * fac01_var_1657(jl_var_1733, lay_var_1688)
          fac201_var_1701 = fk2_var_1709 * fac01_var_1657(jl_var_1733, lay_var_1688)
          fac011_var_1702 = fk0_var_1707 * fac11_var_1659(jl_var_1733, lay_var_1688)
          fac111_var_1703 = fk1_var_1708 * fac11_var_1659(jl_var_1733, lay_var_1688)
          fac211_var_1704 = fk2_var_1709 * fac11_var_1659(jl_var_1733, lay_var_1688)
        ELSE
          fac001_var_1699 = (1.0D0 - fs1_var_1717) * fac01_var_1657(jl_var_1733, lay_var_1688)
          fac011_var_1702 = (1.0D0 - fs1_var_1717) * fac11_var_1659(jl_var_1733, lay_var_1688)
          fac101_var_1700 = fs1_var_1717 * fac01_var_1657(jl_var_1733, lay_var_1688)
          fac111_var_1703 = fs1_var_1717 * fac11_var_1659(jl_var_1733, lay_var_1688)
          fac201_var_1701 = 0.0D0
          fac211_var_1704 = 0.0D0
        END IF
        IF (specparm_var_1716 .LT. 0.125D0) THEN
          tau_major_var_1712(1 : ng4) = speccomb_var_1679 * (fac000_var_1693 * absa_var_189(ind0_var_1682, 1 : 14) + fac100_var_1694 * absa_var_189(ind0_var_1682 + 1, 1 : 14) + fac200_var_1695 * absa_var_189(ind0_var_1682 + 2, 1 : 14) + fac010_var_1696 * absa_var_189(ind0_var_1682 + 9, 1 : 14) + fac110_var_1697 * absa_var_189(ind0_var_1682 + 10, 1 : 14) + fac210_var_1698 * absa_var_189(ind0_var_1682 + 11, 1 : 14))
        ELSE IF (specparm_var_1716 .GT. 0.875D0) THEN
          tau_major_var_1712(1 : ng4) = speccomb_var_1679 * (fac200_var_1695 * absa_var_189(ind0_var_1682 - 1, 1 : 14) + fac100_var_1694 * absa_var_189(ind0_var_1682, 1 : 14) + fac000_var_1693 * absa_var_189(ind0_var_1682 + 1, 1 : 14) + fac210_var_1698 * absa_var_189(ind0_var_1682 + 8, 1 : 14) + fac110_var_1697 * absa_var_189(ind0_var_1682 + 9, 1 : 14) + fac010_var_1696 * absa_var_189(ind0_var_1682 + 10, 1 : 14))
        ELSE
          tau_major_var_1712(1 : ng4) = speccomb_var_1679 * (fac000_var_1693 * absa_var_189(ind0_var_1682, 1 : 14) + fac100_var_1694 * absa_var_189(ind0_var_1682 + 1, 1 : 14) + fac010_var_1696 * absa_var_189(ind0_var_1682 + 9, 1 : 14) + fac110_var_1697 * absa_var_189(ind0_var_1682 + 10, 1 : 14))
        END IF
        IF (specparm1_var_1719 .LT. 0.125D0) THEN
          tau_major1_var_1713(1 : ng4) = speccomb1_var_1680 * (fac001_var_1699 * absa_var_189(ind1_var_1683, 1 : 14) + fac101_var_1700 * absa_var_189(ind1_var_1683 + 1, 1 : 14) + fac201_var_1701 * absa_var_189(ind1_var_1683 + 2, 1 : 14) + fac011_var_1702 * absa_var_189(ind1_var_1683 + 9, 1 : 14) + fac111_var_1703 * absa_var_189(ind1_var_1683 + 10, 1 : 14) + fac211_var_1704 * absa_var_189(ind1_var_1683 + 11, 1 : 14))
        ELSE IF (specparm1_var_1719 .GT. 0.875D0) THEN
          tau_major1_var_1713(1 : ng4) = speccomb1_var_1680 * (fac201_var_1701 * absa_var_189(ind1_var_1683 - 1, 1 : 14) + fac101_var_1700 * absa_var_189(ind1_var_1683, 1 : 14) + fac001_var_1699 * absa_var_189(ind1_var_1683 + 1, 1 : 14) + fac211_var_1704 * absa_var_189(ind1_var_1683 + 8, 1 : 14) + fac111_var_1703 * absa_var_189(ind1_var_1683 + 9, 1 : 14) + fac011_var_1702 * absa_var_189(ind1_var_1683 + 10, 1 : 14))
        ELSE
          tau_major1_var_1713(1 : ng4) = speccomb1_var_1680 * (fac001_var_1699 * absa_var_189(ind1_var_1683, 1 : 14) + fac101_var_1700 * absa_var_189(ind1_var_1683 + 1, 1 : 14) + fac011_var_1702 * absa_var_189(ind1_var_1683 + 9, 1 : 14) + fac111_var_1703 * absa_var_189(ind1_var_1683 + 10, 1 : 14))
        END IF
        DO ig_var_1686 = 1, 14
          tauself_var_1711 = selffac_var_1668(jl_var_1733, lay_var_1688) * (selfref_var_191(inds_var_1684, ig_var_1686) + selffrac_var_1669(jl_var_1733, lay_var_1688) * (selfref_var_191(inds_var_1684 + 1, ig_var_1686) - selfref_var_191(inds_var_1684, ig_var_1686)))
          taufor_var_1710 = forfac_var_1677(jl_var_1733, lay_var_1688) * (forref_var_192(indf_var_1685, ig_var_1686) + forfrac_var_1678(jl_var_1733, lay_var_1688) * (forref_var_192(indf_var_1685 + 1, ig_var_1686) - forref_var_192(indf_var_1685, ig_var_1686)))
          taug_var_1654(jl_var_1733, 38 + ig_var_1686, lay_var_1688) = tau_major_var_1712(ig_var_1686) + tau_major1_var_1713(ig_var_1686) + tauself_var_1711 + taufor_var_1710
          fracs_var_1671(jl_var_1733, 38 + ig_var_1686, lay_var_1688) = fracrefa_var_187(ig_var_1686, jpl_var_1690) + fpl_var_1720 * (fracrefa_var_187(ig_var_1686, jpl_var_1690 + 1) - fracrefa_var_187(ig_var_1686, jpl_var_1690))
        END DO
      END DO
      ixc0_var_1730 = kfdia_var_1652 - kidia_var_1651 + 1 - ixc0_var_1730
      DO ixp_var_1731 = 1, ixc0_var_1730
        jl_var_1733 = ixhigh_var_1727(ixp_var_1731, lay_var_1688)
        speccomb_var_1679 = colo3_var_1666(jl_var_1733, lay_var_1688) + rat_o3co2_var_1674(jl_var_1733, lay_var_1688) * colco2_var_1665(jl_var_1733, lay_var_1688)
        specparm_var_1716 = MIN(colo3_var_1666(jl_var_1733, lay_var_1688) / speccomb_var_1679, oneminus_var_1663)
        specmult_var_1715 = 4.0D0 * (specparm_var_1716)
        js_var_1687 = 1 + INT(specmult_var_1715)
        fs_var_1714 = ((specmult_var_1715) - AINT((specmult_var_1715)))
        speccomb1_var_1680 = colo3_var_1666(jl_var_1733, lay_var_1688) + rat_o3co2_1_var_1675(jl_var_1733, lay_var_1688) * colco2_var_1665(jl_var_1733, lay_var_1688)
        specparm1_var_1719 = MIN(colo3_var_1666(jl_var_1733, lay_var_1688) / speccomb1_var_1680, oneminus_var_1663)
        specmult1_var_1718 = 4.0D0 * (specparm1_var_1719)
        js1_var_1689 = 1 + INT(specmult1_var_1718)
        fs1_var_1717 = ((specmult1_var_1718) - AINT((specmult1_var_1718)))
        fac000_var_1693 = (1.0D0 - fs_var_1714) * fac00_var_1656(jl_var_1733, lay_var_1688)
        fac010_var_1696 = (1.0D0 - fs_var_1714) * fac10_var_1658(jl_var_1733, lay_var_1688)
        fac100_var_1694 = fs_var_1714 * fac00_var_1656(jl_var_1733, lay_var_1688)
        fac110_var_1697 = fs_var_1714 * fac10_var_1658(jl_var_1733, lay_var_1688)
        fac001_var_1699 = (1.0D0 - fs1_var_1717) * fac01_var_1657(jl_var_1733, lay_var_1688)
        fac011_var_1702 = (1.0D0 - fs1_var_1717) * fac11_var_1659(jl_var_1733, lay_var_1688)
        fac101_var_1700 = fs1_var_1717 * fac01_var_1657(jl_var_1733, lay_var_1688)
        fac111_var_1703 = fs1_var_1717 * fac11_var_1659(jl_var_1733, lay_var_1688)
        speccomb_planck_var_1681 = colo3_var_1666(jl_var_1733, lay_var_1688) + refrat_planck_b_var_1692 * colco2_var_1665(jl_var_1733, lay_var_1688)
        specparm_planck_var_1722 = MIN(colo3_var_1666(jl_var_1733, lay_var_1688) / speccomb_planck_var_1681, oneminus_var_1663)
        specmult_planck_var_1721 = 4.0D0 * specparm_planck_var_1722
        jpl_var_1690 = 1 + INT(specmult_planck_var_1721)
        fpl_var_1720 = ((specmult_planck_var_1721) - AINT((specmult_planck_var_1721)))
        ind0_var_1682 = ((jp_var_1660(jl_var_1733, lay_var_1688) - 13) * 5 + (jt_var_1661(jl_var_1733, lay_var_1688) - 1)) * nspb_var_237(4) + js_var_1687
        ind1_var_1683 = ((jp_var_1660(jl_var_1733, lay_var_1688) - 12) * 5 + (jt1_var_1662(jl_var_1733, lay_var_1688) - 1)) * nspb_var_237(4) + js1_var_1689
        DO ig_var_1686 = 1, 14
          taug_var_1654(jl_var_1733, 38 + ig_var_1686, lay_var_1688) = speccomb_var_1679 * (fac000_var_1693 * absb_var_190(ind0_var_1682, ig_var_1686) + fac100_var_1694 * absb_var_190(ind0_var_1682 + 1, ig_var_1686) + fac010_var_1696 * absb_var_190(ind0_var_1682 + 5, ig_var_1686) + fac110_var_1697 * absb_var_190(ind0_var_1682 + 6, ig_var_1686)) + speccomb1_var_1680 * (fac001_var_1699 * absb_var_190(ind1_var_1683, ig_var_1686) + fac101_var_1700 * absb_var_190(ind1_var_1683 + 1, ig_var_1686) + fac011_var_1702 * absb_var_190(ind1_var_1683 + 5, ig_var_1686) + fac111_var_1703 * absb_var_190(ind1_var_1683 + 6, ig_var_1686))
          fracs_var_1671(jl_var_1733, 38 + ig_var_1686, lay_var_1688) = fracrefb_var_188(ig_var_1686, jpl_var_1690) + fpl_var_1720 * (fracrefb_var_188(ig_var_1686, jpl_var_1690 + 1) - fracrefb_var_188(ig_var_1686, jpl_var_1690))
        END DO
      END DO
      DO ixp_var_1731 = 1, ixc0_var_1730
        jl_var_1733 = ixhigh_var_1727(ixp_var_1731, lay_var_1688)
        taug_var_1654(jl_var_1733, 46, lay_var_1688) = taug_var_1654(jl_var_1733, 46, lay_var_1688) * 0.92D0
        taug_var_1654(jl_var_1733, 47, lay_var_1688) = taug_var_1654(jl_var_1733, 47, lay_var_1688) * 0.88D0
        taug_var_1654(jl_var_1733, 48, lay_var_1688) = taug_var_1654(jl_var_1733, 48, lay_var_1688) * 1.07D0
        taug_var_1654(jl_var_1733, 49, lay_var_1688) = taug_var_1654(jl_var_1733, 49, lay_var_1688) * 1.1D0
        taug_var_1654(jl_var_1733, 50, lay_var_1688) = taug_var_1654(jl_var_1733, 50, lay_var_1688) * 0.99D0
        taug_var_1654(jl_var_1733, 51, lay_var_1688) = taug_var_1654(jl_var_1733, 51, lay_var_1688) * 0.88D0
        taug_var_1654(jl_var_1733, 52, lay_var_1688) = taug_var_1654(jl_var_1733, 52, lay_var_1688) * 0.943D0
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol4
SUBROUTINE rrtm_taumol6(kidia_var_1734, kfdia_var_1735, klev_var_1736, taug_var_1737, wx_var_1738, p_tauaerl_var_1739, fac00_var_1740, fac01_var_1741, fac10_var_1742, fac11_var_1743, forfac_var_1756, forfrac_var_1757, indfor_var_1755, jp_var_1744, jt_var_1745, jt1_var_1746, colh2o_var_1747, colco2_var_1748, coldry_var_1749, laytrop_var_1750, selffac_var_1751, selffrac_var_1752, indself_var_1753, fracs_var_1754, minorfrac_var_1758, indminor_var_1759)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_236
  USE yoerrta6, ONLY: absa_var_202, cfc11adj, cfc12_var_201, forref_var_205, fracrefa_var_200, ka_mco2_var_204, selfref_var_203
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1734
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1735
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1736
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1737(kidia_var_1734 : kfdia_var_1735, 140, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: wx_var_1738(kidia_var_1734 : kfdia_var_1735, 4, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1739(kidia_var_1734 : kfdia_var_1735, klev_var_1736, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1740(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1741(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1742(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1743(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1744(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1745(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1746(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1747(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1748(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1749(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1750(kidia_var_1734 : kfdia_var_1735)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1751(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1752(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1753(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1754(kidia_var_1734 : kfdia_var_1735, 140, klev_var_1736)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1755(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1756(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1757(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1758(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1759(kidia_var_1734 : kfdia_var_1735, klev_var_1736)
  INTEGER(KIND = 4) :: ind0_var_1760, ind1_var_1761, inds_var_1762, indf_var_1763, indm_var_1764
  INTEGER(KIND = 4) :: ig_var_1765, lay_var_1766
  REAL(KIND = 8) :: adjfac_var_1767, adjcolco2_var_1768, ratco2_var_1769, chi_co2_var_1770
  REAL(KIND = 8) :: taufor_var_1771, tauself_var_1772, absco2_var_1773
  INTEGER(KIND = 4) :: laytrop_min_var_1774, laytrop_max_var_1775
  INTEGER(KIND = 4) :: ixc_var_1776(klev_var_1736), ixlow_var_1777(kfdia_var_1735, klev_var_1736), ixhigh_var_1778(kfdia_var_1735, klev_var_1736)
  INTEGER(KIND = 4) :: ich_var_1779, icl_var_1780, ixc0_var_1781, ixp_var_1782, jc_var_1783, jl_var_1784
  laytrop_min_var_1774 = MINVAL(laytrop_var_1750)
  laytrop_max_var_1775 = MAXVAL(laytrop_var_1750)
  ixlow_var_1777 = 0
  ixhigh_var_1778 = 0
  ixc_var_1776 = 0
  DO lay_var_1766 = laytrop_min_var_1774 + 1, laytrop_max_var_1775
    icl_var_1780 = 0
    ich_var_1779 = 0
    DO jc_var_1783 = kidia_var_1734, kfdia_var_1735
      IF (lay_var_1766 <= laytrop_var_1750(jc_var_1783)) THEN
        icl_var_1780 = icl_var_1780 + 1
        ixlow_var_1777(icl_var_1780, lay_var_1766) = jc_var_1783
      ELSE
        ich_var_1779 = ich_var_1779 + 1
        ixhigh_var_1778(ich_var_1779, lay_var_1766) = jc_var_1783
      END IF
    END DO
    ixc_var_1776(lay_var_1766) = icl_var_1780
  END DO
  DO lay_var_1766 = 1, laytrop_min_var_1774
    DO jl_var_1784 = kidia_var_1734, kfdia_var_1735
      chi_co2_var_1770 = colco2_var_1748(jl_var_1784, lay_var_1766) / (coldry_var_1749(jl_var_1784, lay_var_1766))
      ratco2_var_1769 = 1D+20 * chi_co2_var_1770 / chi_mls(2, jp_var_1744(jl_var_1784, lay_var_1766) + 1)
      IF (ratco2_var_1769 .GT. 3.0D0) THEN
        adjfac_var_1767 = 2.0D0 + (ratco2_var_1769 - 2.0D0) ** 0.77D0
        adjcolco2_var_1768 = adjfac_var_1767 * chi_mls(2, jp_var_1744(jl_var_1784, lay_var_1766) + 1) * coldry_var_1749(jl_var_1784, lay_var_1766) * 1D-20
      ELSE
        adjcolco2_var_1768 = colco2_var_1748(jl_var_1784, lay_var_1766)
      END IF
      ind0_var_1760 = ((jp_var_1744(jl_var_1784, lay_var_1766) - 1) * 5 + (jt_var_1745(jl_var_1784, lay_var_1766) - 1)) * nspa_var_236(6) + 1
      ind1_var_1761 = (jp_var_1744(jl_var_1784, lay_var_1766) * 5 + (jt1_var_1746(jl_var_1784, lay_var_1766) - 1)) * nspa_var_236(6) + 1
      inds_var_1762 = indself_var_1753(jl_var_1784, lay_var_1766)
      indf_var_1763 = indfor_var_1755(jl_var_1784, lay_var_1766)
      indm_var_1764 = indminor_var_1759(jl_var_1784, lay_var_1766)
      DO ig_var_1765 = 1, 8
        tauself_var_1772 = selffac_var_1751(jl_var_1784, lay_var_1766) * (selfref_var_203(inds_var_1762, ig_var_1765) + selffrac_var_1752(jl_var_1784, lay_var_1766) * (selfref_var_203(inds_var_1762 + 1, ig_var_1765) - selfref_var_203(inds_var_1762, ig_var_1765)))
        taufor_var_1771 = forfac_var_1756(jl_var_1784, lay_var_1766) * (forref_var_205(indf_var_1763, ig_var_1765) + forfrac_var_1757(jl_var_1784, lay_var_1766) * (forref_var_205(indf_var_1763 + 1, ig_var_1765) - forref_var_205(indf_var_1763, ig_var_1765)))
        absco2_var_1773 = (ka_mco2_var_204(indm_var_1764, ig_var_1765) + minorfrac_var_1758(jl_var_1784, lay_var_1766) * (ka_mco2_var_204(indm_var_1764 + 1, ig_var_1765) - ka_mco2_var_204(indm_var_1764, ig_var_1765)))
        taug_var_1737(jl_var_1784, 68 + ig_var_1765, lay_var_1766) = colh2o_var_1747(jl_var_1784, lay_var_1766) * (fac00_var_1740(jl_var_1784, lay_var_1766) * absa_var_202(ind0_var_1760, ig_var_1765) + fac10_var_1742(jl_var_1784, lay_var_1766) * absa_var_202(ind0_var_1760 + 1, ig_var_1765) + fac01_var_1741(jl_var_1784, lay_var_1766) * absa_var_202(ind1_var_1761, ig_var_1765) + fac11_var_1743(jl_var_1784, lay_var_1766) * absa_var_202(ind1_var_1761 + 1, ig_var_1765)) + tauself_var_1772 + taufor_var_1771 + adjcolco2_var_1768 * absco2_var_1773 + wx_var_1738(jl_var_1784, 2, lay_var_1766) * cfc11adj(ig_var_1765) + wx_var_1738(jl_var_1784, 3, lay_var_1766) * cfc12_var_201(ig_var_1765)
        fracs_var_1754(jl_var_1784, 68 + ig_var_1765, lay_var_1766) = fracrefa_var_200(ig_var_1765)
      END DO
    END DO
  END DO
  DO ig_var_1765 = 1, 8
    DO lay_var_1766 = laytrop_max_var_1775 + 1, klev_var_1736
      DO jl_var_1784 = kidia_var_1734, kfdia_var_1735
        taug_var_1737(jl_var_1784, 68 + ig_var_1765, lay_var_1766) = 0.0D0 + wx_var_1738(jl_var_1784, 2, lay_var_1766) * cfc11adj(ig_var_1765) + wx_var_1738(jl_var_1784, 3, lay_var_1766) * cfc12_var_201(ig_var_1765)
        fracs_var_1754(jl_var_1784, 68 + ig_var_1765, lay_var_1766) = fracrefa_var_200(ig_var_1765)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1775 /= laytrop_min_var_1774) THEN
    DO lay_var_1766 = laytrop_min_var_1774 + 1, laytrop_max_var_1775
      ixc0_var_1781 = ixc_var_1776(lay_var_1766)
      DO ixp_var_1782 = 1, ixc0_var_1781
        jl_var_1784 = ixlow_var_1777(ixp_var_1782, lay_var_1766)
        chi_co2_var_1770 = colco2_var_1748(jl_var_1784, lay_var_1766) / (coldry_var_1749(jl_var_1784, lay_var_1766))
        ratco2_var_1769 = 1D+20 * chi_co2_var_1770 / chi_mls(2, jp_var_1744(jl_var_1784, lay_var_1766) + 1)
        IF (ratco2_var_1769 .GT. 3.0D0) THEN
          adjfac_var_1767 = 2.0D0 + (ratco2_var_1769 - 2.0D0) ** 0.77D0
          adjcolco2_var_1768 = adjfac_var_1767 * chi_mls(2, jp_var_1744(jl_var_1784, lay_var_1766) + 1) * coldry_var_1749(jl_var_1784, lay_var_1766) * 1D-20
        ELSE
          adjcolco2_var_1768 = colco2_var_1748(jl_var_1784, lay_var_1766)
        END IF
        ind0_var_1760 = ((jp_var_1744(jl_var_1784, lay_var_1766) - 1) * 5 + (jt_var_1745(jl_var_1784, lay_var_1766) - 1)) * nspa_var_236(6) + 1
        ind1_var_1761 = (jp_var_1744(jl_var_1784, lay_var_1766) * 5 + (jt1_var_1746(jl_var_1784, lay_var_1766) - 1)) * nspa_var_236(6) + 1
        inds_var_1762 = indself_var_1753(jl_var_1784, lay_var_1766)
        indf_var_1763 = indfor_var_1755(jl_var_1784, lay_var_1766)
        indm_var_1764 = indminor_var_1759(jl_var_1784, lay_var_1766)
        DO ig_var_1765 = 1, 8
          tauself_var_1772 = selffac_var_1751(jl_var_1784, lay_var_1766) * (selfref_var_203(inds_var_1762, ig_var_1765) + selffrac_var_1752(jl_var_1784, lay_var_1766) * (selfref_var_203(inds_var_1762 + 1, ig_var_1765) - selfref_var_203(inds_var_1762, ig_var_1765)))
          taufor_var_1771 = forfac_var_1756(jl_var_1784, lay_var_1766) * (forref_var_205(indf_var_1763, ig_var_1765) + forfrac_var_1757(jl_var_1784, lay_var_1766) * (forref_var_205(indf_var_1763 + 1, ig_var_1765) - forref_var_205(indf_var_1763, ig_var_1765)))
          absco2_var_1773 = (ka_mco2_var_204(indm_var_1764, ig_var_1765) + minorfrac_var_1758(jl_var_1784, lay_var_1766) * (ka_mco2_var_204(indm_var_1764 + 1, ig_var_1765) - ka_mco2_var_204(indm_var_1764, ig_var_1765)))
          taug_var_1737(jl_var_1784, 68 + ig_var_1765, lay_var_1766) = colh2o_var_1747(jl_var_1784, lay_var_1766) * (fac00_var_1740(jl_var_1784, lay_var_1766) * absa_var_202(ind0_var_1760, ig_var_1765) + fac10_var_1742(jl_var_1784, lay_var_1766) * absa_var_202(ind0_var_1760 + 1, ig_var_1765) + fac01_var_1741(jl_var_1784, lay_var_1766) * absa_var_202(ind1_var_1761, ig_var_1765) + fac11_var_1743(jl_var_1784, lay_var_1766) * absa_var_202(ind1_var_1761 + 1, ig_var_1765)) + tauself_var_1772 + taufor_var_1771 + adjcolco2_var_1768 * absco2_var_1773 + wx_var_1738(jl_var_1784, 2, lay_var_1766) * cfc11adj(ig_var_1765) + wx_var_1738(jl_var_1784, 3, lay_var_1766) * cfc12_var_201(ig_var_1765)
          fracs_var_1754(jl_var_1784, 68 + ig_var_1765, lay_var_1766) = fracrefa_var_200(ig_var_1765)
        END DO
      END DO
      ixc0_var_1781 = kfdia_var_1735 - kidia_var_1734 + 1 - ixc0_var_1781
      DO ig_var_1765 = 1, 8
        DO ixp_var_1782 = 1, ixc0_var_1781
          jl_var_1784 = ixhigh_var_1778(ixp_var_1782, lay_var_1766)
          taug_var_1737(jl_var_1784, 68 + ig_var_1765, lay_var_1766) = 0.0D0 + wx_var_1738(jl_var_1784, 2, lay_var_1766) * cfc11adj(ig_var_1765) + wx_var_1738(jl_var_1784, 3, lay_var_1766) * cfc12_var_201(ig_var_1765)
          fracs_var_1754(jl_var_1784, 68 + ig_var_1765, lay_var_1766) = fracrefa_var_200(ig_var_1765)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol6
SUBROUTINE rrtm_taumol7(kidia_var_1785, kfdia_var_1786, klev_var_1787, taug_var_1788, p_tauaerl_var_1789, fac00_var_1790, fac01_var_1791, fac10_var_1792, fac11_var_1793, forfac_var_1809, forfrac_var_1808, indfor_var_1807, jp_var_1794, jt_var_1795, jt1_var_1796, oneminus_var_1797, colh2o_var_1798, colo3_var_1799, colco2_var_1800, coldry_var_1801, laytrop_var_1802, selffac_var_1803, selffrac_var_1804, indself_var_1805, fracs_var_1806, rat_h2oo3, rat_h2oo3_1, minorfrac_var_1810, indminor_var_1811)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrtm, ONLY: ng7
  USE yoerrta7, ONLY: absa_var_208, absb_var_209, forref_var_213, fracrefa_var_206, fracrefb_var_207, ka_mco2_var_211, kb_mco2_var_212, selfref_var_210
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1785
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1786
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1787
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1788(kidia_var_1785 : kfdia_var_1786, 140, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1789(kidia_var_1785 : kfdia_var_1786, klev_var_1787, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1790(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1791(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1792(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1793(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1794(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1795(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1796(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1797
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1798(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1799(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1800(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1801(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1802(kidia_var_1785 : kfdia_var_1786)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1803(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1804(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1805(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1806(kidia_var_1785 : kfdia_var_1786, 140, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oo3(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oo3_1(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1807(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1808(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1809(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1810(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1811(kidia_var_1785 : kfdia_var_1786, klev_var_1787)
  REAL(KIND = 8) :: speccomb_var_1812, speccomb1_var_1813, speccomb_mco2_var_1814, speccomb_planck_var_1815
  INTEGER(KIND = 4) :: ind0_var_1816, ind1_var_1817, inds_var_1818, indf_var_1819, indm_var_1820
  INTEGER(KIND = 4) :: ig_var_1821, js_var_1822, lay_var_1823, js1_var_1824, jpl_var_1825, jmco2_var_1826
  REAL(KIND = 8) :: refrat_planck_a_var_1827, refrat_m_a_var_1828
  REAL(KIND = 8) :: chi_co2_var_1829, ratco2_var_1830, adjfac_var_1831, adjcolco2_var_1832
  REAL(KIND = 8) :: fac000_var_1833, fac100_var_1834, fac200_var_1835, fac010_var_1836, fac110_var_1837, fac210_var_1838, fac001_var_1839, fac101_var_1840, fac201_var_1841, fac011_var_1842, fac111_var_1843, fac211_var_1844
  REAL(KIND = 8) :: p_var_1845, p4_var_1846, fk0_var_1847, fk1_var_1848, fk2_var_1849
  REAL(KIND = 8) :: taufor_var_1850, tauself_var_1851, tau_major_var_1852(12), tau_major1_var_1853(12), co2m1_var_1854, co2m2_var_1855, absco2_var_1856
  REAL(KIND = 8) :: fs_var_1857, specmult_var_1858, specparm_var_1859, fs1_var_1860, specmult1_var_1861, specparm1_var_1862, fpl_var_1863, specmult_planck_var_1864, specparm_planck_var_1865, fmco2_var_1866, specmult_mco2_var_1867, specparm_mco2_var_1868
  INTEGER(KIND = 4) :: laytrop_min_var_1869, laytrop_max_var_1870
  INTEGER(KIND = 4) :: ixc_var_1871(klev_var_1787), ixlow_var_1872(kfdia_var_1786, klev_var_1787), ixhigh_var_1873(kfdia_var_1786, klev_var_1787)
  INTEGER(KIND = 4) :: ich_var_1874, icl_var_1875, ixc0_var_1876, ixp_var_1877, jc_var_1878, jl_var_1879
  laytrop_min_var_1869 = MINVAL(laytrop_var_1802)
  laytrop_max_var_1870 = MAXVAL(laytrop_var_1802)
  ixlow_var_1872 = 0
  ixhigh_var_1873 = 0
  ixc_var_1871 = 0
  DO lay_var_1823 = laytrop_min_var_1869 + 1, laytrop_max_var_1870
    icl_var_1875 = 0
    ich_var_1874 = 0
    DO jc_var_1878 = kidia_var_1785, kfdia_var_1786
      IF (lay_var_1823 <= laytrop_var_1802(jc_var_1878)) THEN
        icl_var_1875 = icl_var_1875 + 1
        ixlow_var_1872(icl_var_1875, lay_var_1823) = jc_var_1878
      ELSE
        ich_var_1874 = ich_var_1874 + 1
        ixhigh_var_1873(ich_var_1874, lay_var_1823) = jc_var_1878
      END IF
    END DO
    ixc_var_1871(lay_var_1823) = icl_var_1875
  END DO
  refrat_planck_a_var_1827 = chi_mls(1, 3) / chi_mls(3, 3)
  refrat_m_a_var_1828 = chi_mls(1, 3) / chi_mls(3, 3)
  DO lay_var_1823 = 1, laytrop_min_var_1869
    DO jl_var_1879 = kidia_var_1785, kfdia_var_1786
      speccomb_var_1812 = colh2o_var_1798(jl_var_1879, lay_var_1823) + rat_h2oo3(jl_var_1879, lay_var_1823) * colo3_var_1799(jl_var_1879, lay_var_1823)
      specparm_var_1859 = MIN(colh2o_var_1798(jl_var_1879, lay_var_1823) / speccomb_var_1812, oneminus_var_1797)
      specmult_var_1858 = 8.0D0 * (specparm_var_1859)
      js_var_1822 = 1 + INT(specmult_var_1858)
      fs_var_1857 = ((specmult_var_1858) - AINT((specmult_var_1858)))
      speccomb1_var_1813 = colh2o_var_1798(jl_var_1879, lay_var_1823) + rat_h2oo3_1(jl_var_1879, lay_var_1823) * colo3_var_1799(jl_var_1879, lay_var_1823)
      specparm1_var_1862 = MIN(colh2o_var_1798(jl_var_1879, lay_var_1823) / speccomb1_var_1813, oneminus_var_1797)
      specmult1_var_1861 = 8.0D0 * (specparm1_var_1862)
      js1_var_1824 = 1 + INT(specmult1_var_1861)
      fs1_var_1860 = ((specmult1_var_1861) - AINT((specmult1_var_1861)))
      speccomb_mco2_var_1814 = colh2o_var_1798(jl_var_1879, lay_var_1823) + refrat_m_a_var_1828 * colo3_var_1799(jl_var_1879, lay_var_1823)
      specparm_mco2_var_1868 = MIN(colh2o_var_1798(jl_var_1879, lay_var_1823) / speccomb_mco2_var_1814, oneminus_var_1797)
      specmult_mco2_var_1867 = 8.0D0 * specparm_mco2_var_1868
      jmco2_var_1826 = 1 + INT(specmult_mco2_var_1867)
      fmco2_var_1866 = ((specmult_mco2_var_1867) - AINT((specmult_mco2_var_1867)))
      chi_co2_var_1829 = colco2_var_1800(jl_var_1879, lay_var_1823) / (coldry_var_1801(jl_var_1879, lay_var_1823))
      ratco2_var_1830 = 1D+20 * chi_co2_var_1829 / chi_mls(2, jp_var_1794(jl_var_1879, lay_var_1823) + 1)
      IF (ratco2_var_1830 .GT. 3.0D0) THEN
        adjfac_var_1831 = 3.0D0 + (ratco2_var_1830 - 3.0D0) ** 0.79D0
        adjcolco2_var_1832 = adjfac_var_1831 * chi_mls(2, jp_var_1794(jl_var_1879, lay_var_1823) + 1) * coldry_var_1801(jl_var_1879, lay_var_1823) * 1D-20
      ELSE
        adjcolco2_var_1832 = colco2_var_1800(jl_var_1879, lay_var_1823)
      END IF
      speccomb_planck_var_1815 = colh2o_var_1798(jl_var_1879, lay_var_1823) + refrat_planck_a_var_1827 * colo3_var_1799(jl_var_1879, lay_var_1823)
      specparm_planck_var_1865 = MIN(colh2o_var_1798(jl_var_1879, lay_var_1823) / speccomb_planck_var_1815, oneminus_var_1797)
      specmult_planck_var_1864 = 8.0D0 * specparm_planck_var_1865
      jpl_var_1825 = 1 + INT(specmult_planck_var_1864)
      fpl_var_1863 = ((specmult_planck_var_1864) - AINT((specmult_planck_var_1864)))
      ind0_var_1816 = ((jp_var_1794(jl_var_1879, lay_var_1823) - 1) * 5 + (jt_var_1795(jl_var_1879, lay_var_1823) - 1)) * nspa_var_236(7) + js_var_1822
      ind1_var_1817 = (jp_var_1794(jl_var_1879, lay_var_1823) * 5 + (jt1_var_1796(jl_var_1879, lay_var_1823) - 1)) * nspa_var_236(7) + js1_var_1824
      inds_var_1818 = indself_var_1805(jl_var_1879, lay_var_1823)
      indf_var_1819 = indfor_var_1807(jl_var_1879, lay_var_1823)
      indm_var_1820 = indminor_var_1811(jl_var_1879, lay_var_1823)
      IF (specparm_var_1859 .LT. 0.125D0) THEN
        p_var_1845 = fs_var_1857 - 1.0D0
        p4_var_1846 = p_var_1845 ** 4
        fk0_var_1847 = p4_var_1846
        fk1_var_1848 = 1.0D0 - p_var_1845 - 2.0D0 * p4_var_1846
        fk2_var_1849 = p_var_1845 + p4_var_1846
        fac000_var_1833 = fk0_var_1847 * fac00_var_1790(jl_var_1879, lay_var_1823)
        fac100_var_1834 = fk1_var_1848 * fac00_var_1790(jl_var_1879, lay_var_1823)
        fac200_var_1835 = fk2_var_1849 * fac00_var_1790(jl_var_1879, lay_var_1823)
        fac010_var_1836 = fk0_var_1847 * fac10_var_1792(jl_var_1879, lay_var_1823)
        fac110_var_1837 = fk1_var_1848 * fac10_var_1792(jl_var_1879, lay_var_1823)
        fac210_var_1838 = fk2_var_1849 * fac10_var_1792(jl_var_1879, lay_var_1823)
      ELSE IF (specparm_var_1859 .GT. 0.875D0) THEN
        p_var_1845 = - fs_var_1857
        p4_var_1846 = p_var_1845 ** 4
        fk0_var_1847 = p4_var_1846
        fk1_var_1848 = 1.0D0 - p_var_1845 - 2.0D0 * p4_var_1846
        fk2_var_1849 = p_var_1845 + p4_var_1846
        fac000_var_1833 = fk0_var_1847 * fac00_var_1790(jl_var_1879, lay_var_1823)
        fac100_var_1834 = fk1_var_1848 * fac00_var_1790(jl_var_1879, lay_var_1823)
        fac200_var_1835 = fk2_var_1849 * fac00_var_1790(jl_var_1879, lay_var_1823)
        fac010_var_1836 = fk0_var_1847 * fac10_var_1792(jl_var_1879, lay_var_1823)
        fac110_var_1837 = fk1_var_1848 * fac10_var_1792(jl_var_1879, lay_var_1823)
        fac210_var_1838 = fk2_var_1849 * fac10_var_1792(jl_var_1879, lay_var_1823)
      ELSE
        fac000_var_1833 = (1.0D0 - fs_var_1857) * fac00_var_1790(jl_var_1879, lay_var_1823)
        fac010_var_1836 = (1.0D0 - fs_var_1857) * fac10_var_1792(jl_var_1879, lay_var_1823)
        fac100_var_1834 = fs_var_1857 * fac00_var_1790(jl_var_1879, lay_var_1823)
        fac110_var_1837 = fs_var_1857 * fac10_var_1792(jl_var_1879, lay_var_1823)
        fac200_var_1835 = 0.0D0
        fac210_var_1838 = 0.0D0
      END IF
      IF (specparm1_var_1862 .LT. 0.125D0) THEN
        p_var_1845 = fs1_var_1860 - 1.0D0
        p4_var_1846 = p_var_1845 ** 4
        fk0_var_1847 = p4_var_1846
        fk1_var_1848 = 1.0D0 - p_var_1845 - 2.0D0 * p4_var_1846
        fk2_var_1849 = p_var_1845 + p4_var_1846
        fac001_var_1839 = fk0_var_1847 * fac01_var_1791(jl_var_1879, lay_var_1823)
        fac101_var_1840 = fk1_var_1848 * fac01_var_1791(jl_var_1879, lay_var_1823)
        fac201_var_1841 = fk2_var_1849 * fac01_var_1791(jl_var_1879, lay_var_1823)
        fac011_var_1842 = fk0_var_1847 * fac11_var_1793(jl_var_1879, lay_var_1823)
        fac111_var_1843 = fk1_var_1848 * fac11_var_1793(jl_var_1879, lay_var_1823)
        fac211_var_1844 = fk2_var_1849 * fac11_var_1793(jl_var_1879, lay_var_1823)
      ELSE IF (specparm1_var_1862 .GT. 0.875D0) THEN
        p_var_1845 = - fs1_var_1860
        p4_var_1846 = p_var_1845 ** 4
        fk0_var_1847 = p4_var_1846
        fk1_var_1848 = 1.0D0 - p_var_1845 - 2.0D0 * p4_var_1846
        fk2_var_1849 = p_var_1845 + p4_var_1846
        fac001_var_1839 = fk0_var_1847 * fac01_var_1791(jl_var_1879, lay_var_1823)
        fac101_var_1840 = fk1_var_1848 * fac01_var_1791(jl_var_1879, lay_var_1823)
        fac201_var_1841 = fk2_var_1849 * fac01_var_1791(jl_var_1879, lay_var_1823)
        fac011_var_1842 = fk0_var_1847 * fac11_var_1793(jl_var_1879, lay_var_1823)
        fac111_var_1843 = fk1_var_1848 * fac11_var_1793(jl_var_1879, lay_var_1823)
        fac211_var_1844 = fk2_var_1849 * fac11_var_1793(jl_var_1879, lay_var_1823)
      ELSE
        fac001_var_1839 = (1.0D0 - fs1_var_1860) * fac01_var_1791(jl_var_1879, lay_var_1823)
        fac011_var_1842 = (1.0D0 - fs1_var_1860) * fac11_var_1793(jl_var_1879, lay_var_1823)
        fac101_var_1840 = fs1_var_1860 * fac01_var_1791(jl_var_1879, lay_var_1823)
        fac111_var_1843 = fs1_var_1860 * fac11_var_1793(jl_var_1879, lay_var_1823)
        fac201_var_1841 = 0.0D0
        fac211_var_1844 = 0.0D0
      END IF
      IF (specparm_var_1859 .LT. 0.125D0) THEN
        tau_major_var_1852(1 : ng7) = speccomb_var_1812 * (fac000_var_1833 * absa_var_208(ind0_var_1816, 1 : 12) + fac100_var_1834 * absa_var_208(ind0_var_1816 + 1, 1 : 12) + fac200_var_1835 * absa_var_208(ind0_var_1816 + 2, 1 : 12) + fac010_var_1836 * absa_var_208(ind0_var_1816 + 9, 1 : 12) + fac110_var_1837 * absa_var_208(ind0_var_1816 + 10, 1 : 12) + fac210_var_1838 * absa_var_208(ind0_var_1816 + 11, 1 : 12))
      ELSE IF (specparm_var_1859 .GT. 0.875D0) THEN
        tau_major_var_1852(1 : ng7) = speccomb_var_1812 * (fac200_var_1835 * absa_var_208(ind0_var_1816 - 1, 1 : 12) + fac100_var_1834 * absa_var_208(ind0_var_1816, 1 : 12) + fac000_var_1833 * absa_var_208(ind0_var_1816 + 1, 1 : 12) + fac210_var_1838 * absa_var_208(ind0_var_1816 + 8, 1 : 12) + fac110_var_1837 * absa_var_208(ind0_var_1816 + 9, 1 : 12) + fac010_var_1836 * absa_var_208(ind0_var_1816 + 10, 1 : 12))
      ELSE
        tau_major_var_1852(1 : ng7) = speccomb_var_1812 * (fac000_var_1833 * absa_var_208(ind0_var_1816, 1 : 12) + fac100_var_1834 * absa_var_208(ind0_var_1816 + 1, 1 : 12) + fac010_var_1836 * absa_var_208(ind0_var_1816 + 9, 1 : 12) + fac110_var_1837 * absa_var_208(ind0_var_1816 + 10, 1 : 12))
      END IF
      IF (specparm1_var_1862 .LT. 0.125D0) THEN
        tau_major1_var_1853(1 : ng7) = speccomb1_var_1813 * (fac001_var_1839 * absa_var_208(ind1_var_1817, 1 : 12) + fac101_var_1840 * absa_var_208(ind1_var_1817 + 1, 1 : 12) + fac201_var_1841 * absa_var_208(ind1_var_1817 + 2, 1 : 12) + fac011_var_1842 * absa_var_208(ind1_var_1817 + 9, 1 : 12) + fac111_var_1843 * absa_var_208(ind1_var_1817 + 10, 1 : 12) + fac211_var_1844 * absa_var_208(ind1_var_1817 + 11, 1 : 12))
      ELSE IF (specparm1_var_1862 .GT. 0.875D0) THEN
        tau_major1_var_1853(1 : ng7) = speccomb1_var_1813 * (fac201_var_1841 * absa_var_208(ind1_var_1817 - 1, 1 : 12) + fac101_var_1840 * absa_var_208(ind1_var_1817, 1 : 12) + fac001_var_1839 * absa_var_208(ind1_var_1817 + 1, 1 : 12) + fac211_var_1844 * absa_var_208(ind1_var_1817 + 8, 1 : 12) + fac111_var_1843 * absa_var_208(ind1_var_1817 + 9, 1 : 12) + fac011_var_1842 * absa_var_208(ind1_var_1817 + 10, 1 : 12))
      ELSE
        tau_major1_var_1853(1 : ng7) = speccomb1_var_1813 * (fac001_var_1839 * absa_var_208(ind1_var_1817, 1 : 12) + fac101_var_1840 * absa_var_208(ind1_var_1817 + 1, 1 : 12) + fac011_var_1842 * absa_var_208(ind1_var_1817 + 9, 1 : 12) + fac111_var_1843 * absa_var_208(ind1_var_1817 + 10, 1 : 12))
      END IF
      DO ig_var_1821 = 1, 12
        tauself_var_1851 = selffac_var_1803(jl_var_1879, lay_var_1823) * (selfref_var_210(inds_var_1818, ig_var_1821) + selffrac_var_1804(jl_var_1879, lay_var_1823) * (selfref_var_210(inds_var_1818 + 1, ig_var_1821) - selfref_var_210(inds_var_1818, ig_var_1821)))
        taufor_var_1850 = forfac_var_1809(jl_var_1879, lay_var_1823) * (forref_var_213(indf_var_1819, ig_var_1821) + forfrac_var_1808(jl_var_1879, lay_var_1823) * (forref_var_213(indf_var_1819 + 1, ig_var_1821) - forref_var_213(indf_var_1819, ig_var_1821)))
        co2m1_var_1854 = ka_mco2_var_211(jmco2_var_1826, indm_var_1820, ig_var_1821) + fmco2_var_1866 * (ka_mco2_var_211(jmco2_var_1826 + 1, indm_var_1820, ig_var_1821) - ka_mco2_var_211(jmco2_var_1826, indm_var_1820, ig_var_1821))
        co2m2_var_1855 = ka_mco2_var_211(jmco2_var_1826, indm_var_1820 + 1, ig_var_1821) + fmco2_var_1866 * (ka_mco2_var_211(jmco2_var_1826 + 1, indm_var_1820 + 1, ig_var_1821) - ka_mco2_var_211(jmco2_var_1826, indm_var_1820 + 1, ig_var_1821))
        absco2_var_1856 = co2m1_var_1854 + minorfrac_var_1810(jl_var_1879, lay_var_1823) * (co2m2_var_1855 - co2m1_var_1854)
        taug_var_1788(jl_var_1879, 76 + ig_var_1821, lay_var_1823) = tau_major_var_1852(ig_var_1821) + tau_major1_var_1853(ig_var_1821) + tauself_var_1851 + taufor_var_1850 + adjcolco2_var_1832 * absco2_var_1856
        fracs_var_1806(jl_var_1879, 76 + ig_var_1821, lay_var_1823) = fracrefa_var_206(ig_var_1821, jpl_var_1825) + fpl_var_1863 * (fracrefa_var_206(ig_var_1821, jpl_var_1825 + 1) - fracrefa_var_206(ig_var_1821, jpl_var_1825))
      END DO
    END DO
  END DO
  DO lay_var_1823 = laytrop_max_var_1870 + 1, klev_var_1787
    DO jl_var_1879 = kidia_var_1785, kfdia_var_1786
      chi_co2_var_1829 = colco2_var_1800(jl_var_1879, lay_var_1823) / (coldry_var_1801(jl_var_1879, lay_var_1823))
      ratco2_var_1830 = 1D+20 * chi_co2_var_1829 / chi_mls(2, jp_var_1794(jl_var_1879, lay_var_1823) + 1)
      IF (ratco2_var_1830 .GT. 3.0D0) THEN
        adjfac_var_1831 = 2.0D0 + (ratco2_var_1830 - 2.0D0) ** 0.79D0
        adjcolco2_var_1832 = adjfac_var_1831 * chi_mls(2, jp_var_1794(jl_var_1879, lay_var_1823) + 1) * coldry_var_1801(jl_var_1879, lay_var_1823) * 1D-20
      ELSE
        adjcolco2_var_1832 = colco2_var_1800(jl_var_1879, lay_var_1823)
      END IF
      ind0_var_1816 = ((jp_var_1794(jl_var_1879, lay_var_1823) - 13) * 5 + (jt_var_1795(jl_var_1879, lay_var_1823) - 1)) * nspb_var_237(7) + 1
      ind1_var_1817 = ((jp_var_1794(jl_var_1879, lay_var_1823) - 12) * 5 + (jt1_var_1796(jl_var_1879, lay_var_1823) - 1)) * nspb_var_237(7) + 1
      indm_var_1820 = indminor_var_1811(jl_var_1879, lay_var_1823)
      DO ig_var_1821 = 1, 12
        absco2_var_1856 = kb_mco2_var_212(indm_var_1820, ig_var_1821) + minorfrac_var_1810(jl_var_1879, lay_var_1823) * (kb_mco2_var_212(indm_var_1820 + 1, ig_var_1821) - kb_mco2_var_212(indm_var_1820, ig_var_1821))
        taug_var_1788(jl_var_1879, 76 + ig_var_1821, lay_var_1823) = colo3_var_1799(jl_var_1879, lay_var_1823) * (fac00_var_1790(jl_var_1879, lay_var_1823) * absb_var_209(ind0_var_1816, ig_var_1821) + fac10_var_1792(jl_var_1879, lay_var_1823) * absb_var_209(ind0_var_1816 + 1, ig_var_1821) + fac01_var_1791(jl_var_1879, lay_var_1823) * absb_var_209(ind1_var_1817, ig_var_1821) + fac11_var_1793(jl_var_1879, lay_var_1823) * absb_var_209(ind1_var_1817 + 1, ig_var_1821)) + adjcolco2_var_1832 * absco2_var_1856
        fracs_var_1806(jl_var_1879, 76 + ig_var_1821, lay_var_1823) = fracrefb_var_207(ig_var_1821)
      END DO
    END DO
  END DO
  DO lay_var_1823 = laytrop_max_var_1870 + 1, klev_var_1787
    DO jl_var_1879 = kidia_var_1785, kfdia_var_1786
      taug_var_1788(jl_var_1879, 82, lay_var_1823) = taug_var_1788(jl_var_1879, 82, lay_var_1823) * 0.92D0
      taug_var_1788(jl_var_1879, 83, lay_var_1823) = taug_var_1788(jl_var_1879, 83, lay_var_1823) * 0.88D0
      taug_var_1788(jl_var_1879, 84, lay_var_1823) = taug_var_1788(jl_var_1879, 84, lay_var_1823) * 1.07D0
      taug_var_1788(jl_var_1879, 85, lay_var_1823) = taug_var_1788(jl_var_1879, 85, lay_var_1823) * 1.1D0
      taug_var_1788(jl_var_1879, 86, lay_var_1823) = taug_var_1788(jl_var_1879, 86, lay_var_1823) * 0.99D0
      taug_var_1788(jl_var_1879, 87, lay_var_1823) = taug_var_1788(jl_var_1879, 87, lay_var_1823) * 0.855D0
    END DO
  END DO
  IF (laytrop_max_var_1870 /= laytrop_min_var_1869) THEN
    DO lay_var_1823 = laytrop_min_var_1869 + 1, laytrop_max_var_1870
      ixc0_var_1876 = ixc_var_1871(lay_var_1823)
      DO ixp_var_1877 = 1, ixc0_var_1876
        jl_var_1879 = ixlow_var_1872(ixp_var_1877, lay_var_1823)
        speccomb_var_1812 = colh2o_var_1798(jl_var_1879, lay_var_1823) + rat_h2oo3(jl_var_1879, lay_var_1823) * colo3_var_1799(jl_var_1879, lay_var_1823)
        specparm_var_1859 = MIN(colh2o_var_1798(jl_var_1879, lay_var_1823) / speccomb_var_1812, oneminus_var_1797)
        specmult_var_1858 = 8.0D0 * (specparm_var_1859)
        js_var_1822 = 1 + INT(specmult_var_1858)
        fs_var_1857 = ((specmult_var_1858) - AINT((specmult_var_1858)))
        speccomb1_var_1813 = colh2o_var_1798(jl_var_1879, lay_var_1823) + rat_h2oo3_1(jl_var_1879, lay_var_1823) * colo3_var_1799(jl_var_1879, lay_var_1823)
        specparm1_var_1862 = MIN(colh2o_var_1798(jl_var_1879, lay_var_1823) / speccomb1_var_1813, oneminus_var_1797)
        specmult1_var_1861 = 8.0D0 * (specparm1_var_1862)
        js1_var_1824 = 1 + INT(specmult1_var_1861)
        fs1_var_1860 = ((specmult1_var_1861) - AINT((specmult1_var_1861)))
        speccomb_mco2_var_1814 = colh2o_var_1798(jl_var_1879, lay_var_1823) + refrat_m_a_var_1828 * colo3_var_1799(jl_var_1879, lay_var_1823)
        specparm_mco2_var_1868 = MIN(colh2o_var_1798(jl_var_1879, lay_var_1823) / speccomb_mco2_var_1814, oneminus_var_1797)
        specmult_mco2_var_1867 = 8.0D0 * specparm_mco2_var_1868
        jmco2_var_1826 = 1 + INT(specmult_mco2_var_1867)
        fmco2_var_1866 = ((specmult_mco2_var_1867) - AINT((specmult_mco2_var_1867)))
        chi_co2_var_1829 = colco2_var_1800(jl_var_1879, lay_var_1823) / (coldry_var_1801(jl_var_1879, lay_var_1823))
        ratco2_var_1830 = 1D+20 * chi_co2_var_1829 / chi_mls(2, jp_var_1794(jl_var_1879, lay_var_1823) + 1)
        IF (ratco2_var_1830 .GT. 3.0D0) THEN
          adjfac_var_1831 = 3.0D0 + (ratco2_var_1830 - 3.0D0) ** 0.79D0
          adjcolco2_var_1832 = adjfac_var_1831 * chi_mls(2, jp_var_1794(jl_var_1879, lay_var_1823) + 1) * coldry_var_1801(jl_var_1879, lay_var_1823) * 1D-20
        ELSE
          adjcolco2_var_1832 = colco2_var_1800(jl_var_1879, lay_var_1823)
        END IF
        speccomb_planck_var_1815 = colh2o_var_1798(jl_var_1879, lay_var_1823) + refrat_planck_a_var_1827 * colo3_var_1799(jl_var_1879, lay_var_1823)
        specparm_planck_var_1865 = MIN(colh2o_var_1798(jl_var_1879, lay_var_1823) / speccomb_planck_var_1815, oneminus_var_1797)
        specmult_planck_var_1864 = 8.0D0 * specparm_planck_var_1865
        jpl_var_1825 = 1 + INT(specmult_planck_var_1864)
        fpl_var_1863 = ((specmult_planck_var_1864) - AINT((specmult_planck_var_1864)))
        ind0_var_1816 = ((jp_var_1794(jl_var_1879, lay_var_1823) - 1) * 5 + (jt_var_1795(jl_var_1879, lay_var_1823) - 1)) * nspa_var_236(7) + js_var_1822
        ind1_var_1817 = (jp_var_1794(jl_var_1879, lay_var_1823) * 5 + (jt1_var_1796(jl_var_1879, lay_var_1823) - 1)) * nspa_var_236(7) + js1_var_1824
        inds_var_1818 = indself_var_1805(jl_var_1879, lay_var_1823)
        indf_var_1819 = indfor_var_1807(jl_var_1879, lay_var_1823)
        indm_var_1820 = indminor_var_1811(jl_var_1879, lay_var_1823)
        IF (specparm_var_1859 .LT. 0.125D0) THEN
          p_var_1845 = fs_var_1857 - 1.0D0
          p4_var_1846 = p_var_1845 ** 4
          fk0_var_1847 = p4_var_1846
          fk1_var_1848 = 1.0D0 - p_var_1845 - 2.0D0 * p4_var_1846
          fk2_var_1849 = p_var_1845 + p4_var_1846
          fac000_var_1833 = fk0_var_1847 * fac00_var_1790(jl_var_1879, lay_var_1823)
          fac100_var_1834 = fk1_var_1848 * fac00_var_1790(jl_var_1879, lay_var_1823)
          fac200_var_1835 = fk2_var_1849 * fac00_var_1790(jl_var_1879, lay_var_1823)
          fac010_var_1836 = fk0_var_1847 * fac10_var_1792(jl_var_1879, lay_var_1823)
          fac110_var_1837 = fk1_var_1848 * fac10_var_1792(jl_var_1879, lay_var_1823)
          fac210_var_1838 = fk2_var_1849 * fac10_var_1792(jl_var_1879, lay_var_1823)
        ELSE IF (specparm_var_1859 .GT. 0.875D0) THEN
          p_var_1845 = - fs_var_1857
          p4_var_1846 = p_var_1845 ** 4
          fk0_var_1847 = p4_var_1846
          fk1_var_1848 = 1.0D0 - p_var_1845 - 2.0D0 * p4_var_1846
          fk2_var_1849 = p_var_1845 + p4_var_1846
          fac000_var_1833 = fk0_var_1847 * fac00_var_1790(jl_var_1879, lay_var_1823)
          fac100_var_1834 = fk1_var_1848 * fac00_var_1790(jl_var_1879, lay_var_1823)
          fac200_var_1835 = fk2_var_1849 * fac00_var_1790(jl_var_1879, lay_var_1823)
          fac010_var_1836 = fk0_var_1847 * fac10_var_1792(jl_var_1879, lay_var_1823)
          fac110_var_1837 = fk1_var_1848 * fac10_var_1792(jl_var_1879, lay_var_1823)
          fac210_var_1838 = fk2_var_1849 * fac10_var_1792(jl_var_1879, lay_var_1823)
        ELSE
          fac000_var_1833 = (1.0D0 - fs_var_1857) * fac00_var_1790(jl_var_1879, lay_var_1823)
          fac010_var_1836 = (1.0D0 - fs_var_1857) * fac10_var_1792(jl_var_1879, lay_var_1823)
          fac100_var_1834 = fs_var_1857 * fac00_var_1790(jl_var_1879, lay_var_1823)
          fac110_var_1837 = fs_var_1857 * fac10_var_1792(jl_var_1879, lay_var_1823)
          fac200_var_1835 = 0.0D0
          fac210_var_1838 = 0.0D0
        END IF
        IF (specparm1_var_1862 .LT. 0.125D0) THEN
          p_var_1845 = fs1_var_1860 - 1.0D0
          p4_var_1846 = p_var_1845 ** 4
          fk0_var_1847 = p4_var_1846
          fk1_var_1848 = 1.0D0 - p_var_1845 - 2.0D0 * p4_var_1846
          fk2_var_1849 = p_var_1845 + p4_var_1846
          fac001_var_1839 = fk0_var_1847 * fac01_var_1791(jl_var_1879, lay_var_1823)
          fac101_var_1840 = fk1_var_1848 * fac01_var_1791(jl_var_1879, lay_var_1823)
          fac201_var_1841 = fk2_var_1849 * fac01_var_1791(jl_var_1879, lay_var_1823)
          fac011_var_1842 = fk0_var_1847 * fac11_var_1793(jl_var_1879, lay_var_1823)
          fac111_var_1843 = fk1_var_1848 * fac11_var_1793(jl_var_1879, lay_var_1823)
          fac211_var_1844 = fk2_var_1849 * fac11_var_1793(jl_var_1879, lay_var_1823)
        ELSE IF (specparm1_var_1862 .GT. 0.875D0) THEN
          p_var_1845 = - fs1_var_1860
          p4_var_1846 = p_var_1845 ** 4
          fk0_var_1847 = p4_var_1846
          fk1_var_1848 = 1.0D0 - p_var_1845 - 2.0D0 * p4_var_1846
          fk2_var_1849 = p_var_1845 + p4_var_1846
          fac001_var_1839 = fk0_var_1847 * fac01_var_1791(jl_var_1879, lay_var_1823)
          fac101_var_1840 = fk1_var_1848 * fac01_var_1791(jl_var_1879, lay_var_1823)
          fac201_var_1841 = fk2_var_1849 * fac01_var_1791(jl_var_1879, lay_var_1823)
          fac011_var_1842 = fk0_var_1847 * fac11_var_1793(jl_var_1879, lay_var_1823)
          fac111_var_1843 = fk1_var_1848 * fac11_var_1793(jl_var_1879, lay_var_1823)
          fac211_var_1844 = fk2_var_1849 * fac11_var_1793(jl_var_1879, lay_var_1823)
        ELSE
          fac001_var_1839 = (1.0D0 - fs1_var_1860) * fac01_var_1791(jl_var_1879, lay_var_1823)
          fac011_var_1842 = (1.0D0 - fs1_var_1860) * fac11_var_1793(jl_var_1879, lay_var_1823)
          fac101_var_1840 = fs1_var_1860 * fac01_var_1791(jl_var_1879, lay_var_1823)
          fac111_var_1843 = fs1_var_1860 * fac11_var_1793(jl_var_1879, lay_var_1823)
          fac201_var_1841 = 0.0D0
          fac211_var_1844 = 0.0D0
        END IF
        IF (specparm_var_1859 .LT. 0.125D0) THEN
          tau_major_var_1852(1 : ng7) = speccomb_var_1812 * (fac000_var_1833 * absa_var_208(ind0_var_1816, 1 : 12) + fac100_var_1834 * absa_var_208(ind0_var_1816 + 1, 1 : 12) + fac200_var_1835 * absa_var_208(ind0_var_1816 + 2, 1 : 12) + fac010_var_1836 * absa_var_208(ind0_var_1816 + 9, 1 : 12) + fac110_var_1837 * absa_var_208(ind0_var_1816 + 10, 1 : 12) + fac210_var_1838 * absa_var_208(ind0_var_1816 + 11, 1 : 12))
        ELSE IF (specparm_var_1859 .GT. 0.875D0) THEN
          tau_major_var_1852(1 : ng7) = speccomb_var_1812 * (fac200_var_1835 * absa_var_208(ind0_var_1816 - 1, 1 : 12) + fac100_var_1834 * absa_var_208(ind0_var_1816, 1 : 12) + fac000_var_1833 * absa_var_208(ind0_var_1816 + 1, 1 : 12) + fac210_var_1838 * absa_var_208(ind0_var_1816 + 8, 1 : 12) + fac110_var_1837 * absa_var_208(ind0_var_1816 + 9, 1 : 12) + fac010_var_1836 * absa_var_208(ind0_var_1816 + 10, 1 : 12))
        ELSE
          tau_major_var_1852(1 : ng7) = speccomb_var_1812 * (fac000_var_1833 * absa_var_208(ind0_var_1816, 1 : 12) + fac100_var_1834 * absa_var_208(ind0_var_1816 + 1, 1 : 12) + fac010_var_1836 * absa_var_208(ind0_var_1816 + 9, 1 : 12) + fac110_var_1837 * absa_var_208(ind0_var_1816 + 10, 1 : 12))
        END IF
        IF (specparm1_var_1862 .LT. 0.125D0) THEN
          tau_major1_var_1853(1 : ng7) = speccomb1_var_1813 * (fac001_var_1839 * absa_var_208(ind1_var_1817, 1 : 12) + fac101_var_1840 * absa_var_208(ind1_var_1817 + 1, 1 : 12) + fac201_var_1841 * absa_var_208(ind1_var_1817 + 2, 1 : 12) + fac011_var_1842 * absa_var_208(ind1_var_1817 + 9, 1 : 12) + fac111_var_1843 * absa_var_208(ind1_var_1817 + 10, 1 : 12) + fac211_var_1844 * absa_var_208(ind1_var_1817 + 11, 1 : 12))
        ELSE IF (specparm1_var_1862 .GT. 0.875D0) THEN
          tau_major1_var_1853(1 : ng7) = speccomb1_var_1813 * (fac201_var_1841 * absa_var_208(ind1_var_1817 - 1, 1 : 12) + fac101_var_1840 * absa_var_208(ind1_var_1817, 1 : 12) + fac001_var_1839 * absa_var_208(ind1_var_1817 + 1, 1 : 12) + fac211_var_1844 * absa_var_208(ind1_var_1817 + 8, 1 : 12) + fac111_var_1843 * absa_var_208(ind1_var_1817 + 9, 1 : 12) + fac011_var_1842 * absa_var_208(ind1_var_1817 + 10, 1 : 12))
        ELSE
          tau_major1_var_1853(1 : ng7) = speccomb1_var_1813 * (fac001_var_1839 * absa_var_208(ind1_var_1817, 1 : 12) + fac101_var_1840 * absa_var_208(ind1_var_1817 + 1, 1 : 12) + fac011_var_1842 * absa_var_208(ind1_var_1817 + 9, 1 : 12) + fac111_var_1843 * absa_var_208(ind1_var_1817 + 10, 1 : 12))
        END IF
        DO ig_var_1821 = 1, 12
          tauself_var_1851 = selffac_var_1803(jl_var_1879, lay_var_1823) * (selfref_var_210(inds_var_1818, ig_var_1821) + selffrac_var_1804(jl_var_1879, lay_var_1823) * (selfref_var_210(inds_var_1818 + 1, ig_var_1821) - selfref_var_210(inds_var_1818, ig_var_1821)))
          taufor_var_1850 = forfac_var_1809(jl_var_1879, lay_var_1823) * (forref_var_213(indf_var_1819, ig_var_1821) + forfrac_var_1808(jl_var_1879, lay_var_1823) * (forref_var_213(indf_var_1819 + 1, ig_var_1821) - forref_var_213(indf_var_1819, ig_var_1821)))
          co2m1_var_1854 = ka_mco2_var_211(jmco2_var_1826, indm_var_1820, ig_var_1821) + fmco2_var_1866 * (ka_mco2_var_211(jmco2_var_1826 + 1, indm_var_1820, ig_var_1821) - ka_mco2_var_211(jmco2_var_1826, indm_var_1820, ig_var_1821))
          co2m2_var_1855 = ka_mco2_var_211(jmco2_var_1826, indm_var_1820 + 1, ig_var_1821) + fmco2_var_1866 * (ka_mco2_var_211(jmco2_var_1826 + 1, indm_var_1820 + 1, ig_var_1821) - ka_mco2_var_211(jmco2_var_1826, indm_var_1820 + 1, ig_var_1821))
          absco2_var_1856 = co2m1_var_1854 + minorfrac_var_1810(jl_var_1879, lay_var_1823) * (co2m2_var_1855 - co2m1_var_1854)
          taug_var_1788(jl_var_1879, 76 + ig_var_1821, lay_var_1823) = tau_major_var_1852(ig_var_1821) + tau_major1_var_1853(ig_var_1821) + tauself_var_1851 + taufor_var_1850 + adjcolco2_var_1832 * absco2_var_1856
          fracs_var_1806(jl_var_1879, 76 + ig_var_1821, lay_var_1823) = fracrefa_var_206(ig_var_1821, jpl_var_1825) + fpl_var_1863 * (fracrefa_var_206(ig_var_1821, jpl_var_1825 + 1) - fracrefa_var_206(ig_var_1821, jpl_var_1825))
        END DO
      END DO
      ixc0_var_1876 = kfdia_var_1786 - kidia_var_1785 + 1 - ixc0_var_1876
      DO ixp_var_1877 = 1, ixc0_var_1876
        jl_var_1879 = ixhigh_var_1873(ixp_var_1877, lay_var_1823)
        chi_co2_var_1829 = colco2_var_1800(jl_var_1879, lay_var_1823) / (coldry_var_1801(jl_var_1879, lay_var_1823))
        ratco2_var_1830 = 1D+20 * chi_co2_var_1829 / chi_mls(2, jp_var_1794(jl_var_1879, lay_var_1823) + 1)
        IF (ratco2_var_1830 .GT. 3.0D0) THEN
          adjfac_var_1831 = 2.0D0 + (ratco2_var_1830 - 2.0D0) ** 0.79D0
          adjcolco2_var_1832 = adjfac_var_1831 * chi_mls(2, jp_var_1794(jl_var_1879, lay_var_1823) + 1) * coldry_var_1801(jl_var_1879, lay_var_1823) * 1D-20
        ELSE
          adjcolco2_var_1832 = colco2_var_1800(jl_var_1879, lay_var_1823)
        END IF
        ind0_var_1816 = ((jp_var_1794(jl_var_1879, lay_var_1823) - 13) * 5 + (jt_var_1795(jl_var_1879, lay_var_1823) - 1)) * nspb_var_237(7) + 1
        ind1_var_1817 = ((jp_var_1794(jl_var_1879, lay_var_1823) - 12) * 5 + (jt1_var_1796(jl_var_1879, lay_var_1823) - 1)) * nspb_var_237(7) + 1
        indm_var_1820 = indminor_var_1811(jl_var_1879, lay_var_1823)
        DO ig_var_1821 = 1, 12
          absco2_var_1856 = kb_mco2_var_212(indm_var_1820, ig_var_1821) + minorfrac_var_1810(jl_var_1879, lay_var_1823) * (kb_mco2_var_212(indm_var_1820 + 1, ig_var_1821) - kb_mco2_var_212(indm_var_1820, ig_var_1821))
          taug_var_1788(jl_var_1879, 76 + ig_var_1821, lay_var_1823) = colo3_var_1799(jl_var_1879, lay_var_1823) * (fac00_var_1790(jl_var_1879, lay_var_1823) * absb_var_209(ind0_var_1816, ig_var_1821) + fac10_var_1792(jl_var_1879, lay_var_1823) * absb_var_209(ind0_var_1816 + 1, ig_var_1821) + fac01_var_1791(jl_var_1879, lay_var_1823) * absb_var_209(ind1_var_1817, ig_var_1821) + fac11_var_1793(jl_var_1879, lay_var_1823) * absb_var_209(ind1_var_1817 + 1, ig_var_1821)) + adjcolco2_var_1832 * absco2_var_1856
          fracs_var_1806(jl_var_1879, 76 + ig_var_1821, lay_var_1823) = fracrefb_var_207(ig_var_1821)
        END DO
      END DO
      DO ixp_var_1877 = 1, ixc0_var_1876
        jl_var_1879 = ixhigh_var_1873(ixp_var_1877, lay_var_1823)
        taug_var_1788(jl_var_1879, 82, lay_var_1823) = taug_var_1788(jl_var_1879, 82, lay_var_1823) * 0.92D0
        taug_var_1788(jl_var_1879, 83, lay_var_1823) = taug_var_1788(jl_var_1879, 83, lay_var_1823) * 0.88D0
        taug_var_1788(jl_var_1879, 84, lay_var_1823) = taug_var_1788(jl_var_1879, 84, lay_var_1823) * 1.07D0
        taug_var_1788(jl_var_1879, 85, lay_var_1823) = taug_var_1788(jl_var_1879, 85, lay_var_1823) * 1.1D0
        taug_var_1788(jl_var_1879, 86, lay_var_1823) = taug_var_1788(jl_var_1879, 86, lay_var_1823) * 0.99D0
        taug_var_1788(jl_var_1879, 87, lay_var_1823) = taug_var_1788(jl_var_1879, 87, lay_var_1823) * 0.855D0
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol7
SUBROUTINE rrtm_taumol3(kidia_var_1880, kfdia_var_1881, klev_var_1882, taug_var_1883, p_tauaerl_var_1884, fac00_var_1885, fac01_var_1886, fac10_var_1887, fac11_var_1888, forfac_var_1889, forfrac_var_1906, indfor_var_1905, jp_var_1890, jt_var_1891, jt1_var_1892, oneminus_var_1893, colh2o_var_1894, colco2_var_1895, coln2o_var_1896, coldry_var_1897, laytrop_var_1898, selffac_var_1899, selffrac_var_1900, indself_var_1901, fracs_var_1902, rat_h2oco2_var_1903, rat_h2oco2_1_var_1904, minorfrac_var_1907, indminor_var_1908)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrtm, ONLY: ng3
  USE yoerrta3, ONLY: absa_var_183, absb_var_184, forref_var_186, fracrefa_var_179, fracrefb_var_180, ka_mn2o_var_181, kb_mn2o_var_182, selfref_var_185
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1880
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1881
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1882
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1883(kidia_var_1880 : kfdia_var_1881, 140, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1884(kidia_var_1880 : kfdia_var_1881, klev_var_1882, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1885(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1886(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1887(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1888(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1889(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1890(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1891(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1892(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1893
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1894(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1895(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_1896(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1897(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1898(kidia_var_1880 : kfdia_var_1881)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1899(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1900(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1901(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1902(kidia_var_1880 : kfdia_var_1881, 140, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1903(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1904(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1905(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1906(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1907(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1908(kidia_var_1880 : kfdia_var_1881, klev_var_1882)
  REAL(KIND = 8) :: speccomb_var_1909, speccomb1_var_1910, speccomb_mn2o_var_1911, speccomb_planck_var_1912
  REAL(KIND = 8) :: refrat_planck_a_var_1913, refrat_planck_b_var_1914, refrat_m_a_var_1915, refrat_m_b
  INTEGER(KIND = 4) :: ind0_var_1916, ind1_var_1917, inds_var_1918, indf_var_1919, indm_var_1920
  INTEGER(KIND = 4) :: ig_var_1921, js_var_1922, lay_var_1923, js1_var_1924, jmn2o_var_1925, jpl_var_1926
  REAL(KIND = 8) :: fs_var_1927, specmult_var_1928, specparm_var_1929, fs1_var_1930, specmult1_var_1931, specparm1_var_1932, fmn2o_var_1933, specmult_mn2o_var_1934, specparm_mn2o_var_1935, fpl_var_1936, specmult_planck_var_1937, specparm_planck_var_1938
  REAL(KIND = 8) :: adjfac_var_1939, adjcoln2o_var_1940, ratn2o_var_1941, chi_n2o_var_1942
  REAL(KIND = 8) :: fac000_var_1943, fac100_var_1944, fac200_var_1945, fac010_var_1946, fac110_var_1947, fac210_var_1948, fac001_var_1949, fac101_var_1950, fac201_var_1951, fac011_var_1952, fac111_var_1953, fac211_var_1954
  REAL(KIND = 8) :: p_var_1955, p4_var_1956, fk0_var_1957, fk1_var_1958, fk2_var_1959
  REAL(KIND = 8) :: taufor_var_1960, tauself_var_1961, n2om1_var_1962, n2om2_var_1963, absn2o_var_1964, tau_major_var_1965(16), tau_major1_var_1966(16)
  INTEGER(KIND = 4) :: laytrop_min_var_1967, laytrop_max_var_1968
  INTEGER(KIND = 4) :: ixc_var_1969(klev_var_1882), ixlow_var_1970(kfdia_var_1881, klev_var_1882), ixhigh_var_1971(kfdia_var_1881, klev_var_1882)
  INTEGER(KIND = 4) :: ich_var_1972, icl_var_1973, ixc0_var_1974, ixp_var_1975, jc_var_1976, jl_var_1977
  laytrop_min_var_1967 = MINVAL(laytrop_var_1898)
  laytrop_max_var_1968 = MAXVAL(laytrop_var_1898)
  ixlow_var_1970 = 0
  ixhigh_var_1971 = 0
  ixc_var_1969 = 0
  DO lay_var_1923 = laytrop_min_var_1967 + 1, laytrop_max_var_1968
    icl_var_1973 = 0
    ich_var_1972 = 0
    DO jc_var_1976 = kidia_var_1880, kfdia_var_1881
      IF (lay_var_1923 <= laytrop_var_1898(jc_var_1976)) THEN
        icl_var_1973 = icl_var_1973 + 1
        ixlow_var_1970(icl_var_1973, lay_var_1923) = jc_var_1976
      ELSE
        ich_var_1972 = ich_var_1972 + 1
        ixhigh_var_1971(ich_var_1972, lay_var_1923) = jc_var_1976
      END IF
    END DO
    ixc_var_1969(lay_var_1923) = icl_var_1973
  END DO
  refrat_planck_a_var_1913 = chi_mls(1, 9) / chi_mls(2, 9)
  refrat_planck_b_var_1914 = chi_mls(1, 13) / chi_mls(2, 13)
  refrat_m_a_var_1915 = chi_mls(1, 3) / chi_mls(2, 3)
  refrat_m_b = chi_mls(1, 13) / chi_mls(2, 13)
  DO lay_var_1923 = 1, laytrop_min_var_1967
    DO jl_var_1977 = kidia_var_1880, kfdia_var_1881
      speccomb_var_1909 = colh2o_var_1894(jl_var_1977, lay_var_1923) + rat_h2oco2_var_1903(jl_var_1977, lay_var_1923) * colco2_var_1895(jl_var_1977, lay_var_1923)
      specparm_var_1929 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_var_1909, oneminus_var_1893)
      specmult_var_1928 = 8.0D0 * (specparm_var_1929)
      js_var_1922 = 1 + INT(specmult_var_1928)
      fs_var_1927 = ((specmult_var_1928) - AINT((specmult_var_1928)))
      speccomb1_var_1910 = colh2o_var_1894(jl_var_1977, lay_var_1923) + rat_h2oco2_1_var_1904(jl_var_1977, lay_var_1923) * colco2_var_1895(jl_var_1977, lay_var_1923)
      specparm1_var_1932 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb1_var_1910, oneminus_var_1893)
      specmult1_var_1931 = 8.0D0 * (specparm1_var_1932)
      js1_var_1924 = 1 + INT(specmult1_var_1931)
      fs1_var_1930 = ((specmult1_var_1931) - AINT((specmult1_var_1931)))
      speccomb_mn2o_var_1911 = colh2o_var_1894(jl_var_1977, lay_var_1923) + refrat_m_a_var_1915 * colco2_var_1895(jl_var_1977, lay_var_1923)
      specparm_mn2o_var_1935 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_mn2o_var_1911, oneminus_var_1893)
      specmult_mn2o_var_1934 = 8.0D0 * specparm_mn2o_var_1935
      jmn2o_var_1925 = 1 + INT(specmult_mn2o_var_1934)
      fmn2o_var_1933 = ((specmult_mn2o_var_1934) - AINT((specmult_mn2o_var_1934)))
      chi_n2o_var_1942 = coln2o_var_1896(jl_var_1977, lay_var_1923) / coldry_var_1897(jl_var_1977, lay_var_1923)
      ratn2o_var_1941 = 1D+20 * chi_n2o_var_1942 / chi_mls(4, jp_var_1890(jl_var_1977, lay_var_1923) + 1)
      IF (ratn2o_var_1941 .GT. 1.5D0) THEN
        adjfac_var_1939 = 0.5D0 + (ratn2o_var_1941 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1940 = adjfac_var_1939 * chi_mls(4, jp_var_1890(jl_var_1977, lay_var_1923) + 1) * coldry_var_1897(jl_var_1977, lay_var_1923) * 1D-20
      ELSE
        adjcoln2o_var_1940 = coln2o_var_1896(jl_var_1977, lay_var_1923)
      END IF
      speccomb_planck_var_1912 = colh2o_var_1894(jl_var_1977, lay_var_1923) + refrat_planck_a_var_1913 * colco2_var_1895(jl_var_1977, lay_var_1923)
      specparm_planck_var_1938 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_planck_var_1912, oneminus_var_1893)
      specmult_planck_var_1937 = 8.0D0 * specparm_planck_var_1938
      jpl_var_1926 = 1 + INT(specmult_planck_var_1937)
      fpl_var_1936 = ((specmult_planck_var_1937) - AINT((specmult_planck_var_1937)))
      ind0_var_1916 = ((jp_var_1890(jl_var_1977, lay_var_1923) - 1) * 5 + (jt_var_1891(jl_var_1977, lay_var_1923) - 1)) * nspa_var_236(3) + js_var_1922
      ind1_var_1917 = (jp_var_1890(jl_var_1977, lay_var_1923) * 5 + (jt1_var_1892(jl_var_1977, lay_var_1923) - 1)) * nspa_var_236(3) + js1_var_1924
      inds_var_1918 = indself_var_1901(jl_var_1977, lay_var_1923)
      indf_var_1919 = indfor_var_1905(jl_var_1977, lay_var_1923)
      indm_var_1920 = indminor_var_1908(jl_var_1977, lay_var_1923)
      IF (specparm_var_1929 .LT. 0.125D0) THEN
        p_var_1955 = fs_var_1927 - 1.0D0
        p4_var_1956 = p_var_1955 ** 4
        fk0_var_1957 = p4_var_1956
        fk1_var_1958 = 1.0D0 - p_var_1955 - 2.0D0 * p4_var_1956
        fk2_var_1959 = p_var_1955 + p4_var_1956
        fac000_var_1943 = fk0_var_1957 * fac00_var_1885(jl_var_1977, lay_var_1923)
        fac100_var_1944 = fk1_var_1958 * fac00_var_1885(jl_var_1977, lay_var_1923)
        fac200_var_1945 = fk2_var_1959 * fac00_var_1885(jl_var_1977, lay_var_1923)
        fac010_var_1946 = fk0_var_1957 * fac10_var_1887(jl_var_1977, lay_var_1923)
        fac110_var_1947 = fk1_var_1958 * fac10_var_1887(jl_var_1977, lay_var_1923)
        fac210_var_1948 = fk2_var_1959 * fac10_var_1887(jl_var_1977, lay_var_1923)
      ELSE IF (specparm_var_1929 .GT. 0.875D0) THEN
        p_var_1955 = - fs_var_1927
        p4_var_1956 = p_var_1955 ** 4
        fk0_var_1957 = p4_var_1956
        fk1_var_1958 = 1.0D0 - p_var_1955 - 2.0D0 * p4_var_1956
        fk2_var_1959 = p_var_1955 + p4_var_1956
        fac000_var_1943 = fk0_var_1957 * fac00_var_1885(jl_var_1977, lay_var_1923)
        fac100_var_1944 = fk1_var_1958 * fac00_var_1885(jl_var_1977, lay_var_1923)
        fac200_var_1945 = fk2_var_1959 * fac00_var_1885(jl_var_1977, lay_var_1923)
        fac010_var_1946 = fk0_var_1957 * fac10_var_1887(jl_var_1977, lay_var_1923)
        fac110_var_1947 = fk1_var_1958 * fac10_var_1887(jl_var_1977, lay_var_1923)
        fac210_var_1948 = fk2_var_1959 * fac10_var_1887(jl_var_1977, lay_var_1923)
      ELSE
        fac000_var_1943 = (1.0D0 - fs_var_1927) * fac00_var_1885(jl_var_1977, lay_var_1923)
        fac010_var_1946 = (1.0D0 - fs_var_1927) * fac10_var_1887(jl_var_1977, lay_var_1923)
        fac100_var_1944 = fs_var_1927 * fac00_var_1885(jl_var_1977, lay_var_1923)
        fac110_var_1947 = fs_var_1927 * fac10_var_1887(jl_var_1977, lay_var_1923)
        fac200_var_1945 = 0.0D0
        fac210_var_1948 = 0.0D0
      END IF
      IF (specparm1_var_1932 .LT. 0.125D0) THEN
        p_var_1955 = fs1_var_1930 - 1.0D0
        p4_var_1956 = p_var_1955 ** 4
        fk0_var_1957 = p4_var_1956
        fk1_var_1958 = 1.0D0 - p_var_1955 - 2.0D0 * p4_var_1956
        fk2_var_1959 = p_var_1955 + p4_var_1956
        fac001_var_1949 = fk0_var_1957 * fac01_var_1886(jl_var_1977, lay_var_1923)
        fac101_var_1950 = fk1_var_1958 * fac01_var_1886(jl_var_1977, lay_var_1923)
        fac201_var_1951 = fk2_var_1959 * fac01_var_1886(jl_var_1977, lay_var_1923)
        fac011_var_1952 = fk0_var_1957 * fac11_var_1888(jl_var_1977, lay_var_1923)
        fac111_var_1953 = fk1_var_1958 * fac11_var_1888(jl_var_1977, lay_var_1923)
        fac211_var_1954 = fk2_var_1959 * fac11_var_1888(jl_var_1977, lay_var_1923)
      ELSE IF (specparm1_var_1932 .GT. 0.875D0) THEN
        p_var_1955 = - fs1_var_1930
        p4_var_1956 = p_var_1955 ** 4
        fk0_var_1957 = p4_var_1956
        fk1_var_1958 = 1.0D0 - p_var_1955 - 2.0D0 * p4_var_1956
        fk2_var_1959 = p_var_1955 + p4_var_1956
        fac001_var_1949 = fk0_var_1957 * fac01_var_1886(jl_var_1977, lay_var_1923)
        fac101_var_1950 = fk1_var_1958 * fac01_var_1886(jl_var_1977, lay_var_1923)
        fac201_var_1951 = fk2_var_1959 * fac01_var_1886(jl_var_1977, lay_var_1923)
        fac011_var_1952 = fk0_var_1957 * fac11_var_1888(jl_var_1977, lay_var_1923)
        fac111_var_1953 = fk1_var_1958 * fac11_var_1888(jl_var_1977, lay_var_1923)
        fac211_var_1954 = fk2_var_1959 * fac11_var_1888(jl_var_1977, lay_var_1923)
      ELSE
        fac001_var_1949 = (1.0D0 - fs1_var_1930) * fac01_var_1886(jl_var_1977, lay_var_1923)
        fac011_var_1952 = (1.0D0 - fs1_var_1930) * fac11_var_1888(jl_var_1977, lay_var_1923)
        fac101_var_1950 = fs1_var_1930 * fac01_var_1886(jl_var_1977, lay_var_1923)
        fac111_var_1953 = fs1_var_1930 * fac11_var_1888(jl_var_1977, lay_var_1923)
        fac201_var_1951 = 0.0D0
        fac211_var_1954 = 0.0D0
      END IF
      IF (specparm_var_1929 .LT. 0.125D0) THEN
        tau_major_var_1965(1 : ng3) = speccomb_var_1909 * (fac000_var_1943 * absa_var_183(ind0_var_1916, 1 : 16) + fac100_var_1944 * absa_var_183(ind0_var_1916 + 1, 1 : 16) + fac200_var_1945 * absa_var_183(ind0_var_1916 + 2, 1 : 16) + fac010_var_1946 * absa_var_183(ind0_var_1916 + 9, 1 : 16) + fac110_var_1947 * absa_var_183(ind0_var_1916 + 10, 1 : 16) + fac210_var_1948 * absa_var_183(ind0_var_1916 + 11, 1 : 16))
      ELSE IF (specparm_var_1929 .GT. 0.875D0) THEN
        tau_major_var_1965(1 : ng3) = speccomb_var_1909 * (fac200_var_1945 * absa_var_183(ind0_var_1916 - 1, 1 : 16) + fac100_var_1944 * absa_var_183(ind0_var_1916, 1 : 16) + fac000_var_1943 * absa_var_183(ind0_var_1916 + 1, 1 : 16) + fac210_var_1948 * absa_var_183(ind0_var_1916 + 8, 1 : 16) + fac110_var_1947 * absa_var_183(ind0_var_1916 + 9, 1 : 16) + fac010_var_1946 * absa_var_183(ind0_var_1916 + 10, 1 : 16))
      ELSE
        tau_major_var_1965(1 : ng3) = speccomb_var_1909 * (fac000_var_1943 * absa_var_183(ind0_var_1916, 1 : 16) + fac100_var_1944 * absa_var_183(ind0_var_1916 + 1, 1 : 16) + fac010_var_1946 * absa_var_183(ind0_var_1916 + 9, 1 : 16) + fac110_var_1947 * absa_var_183(ind0_var_1916 + 10, 1 : 16))
      END IF
      IF (specparm1_var_1932 .LT. 0.125D0) THEN
        tau_major1_var_1966(1 : ng3) = speccomb1_var_1910 * (fac001_var_1949 * absa_var_183(ind1_var_1917, 1 : 16) + fac101_var_1950 * absa_var_183(ind1_var_1917 + 1, 1 : 16) + fac201_var_1951 * absa_var_183(ind1_var_1917 + 2, 1 : 16) + fac011_var_1952 * absa_var_183(ind1_var_1917 + 9, 1 : 16) + fac111_var_1953 * absa_var_183(ind1_var_1917 + 10, 1 : 16) + fac211_var_1954 * absa_var_183(ind1_var_1917 + 11, 1 : 16))
      ELSE IF (specparm1_var_1932 .GT. 0.875D0) THEN
        tau_major1_var_1966(1 : ng3) = speccomb1_var_1910 * (fac201_var_1951 * absa_var_183(ind1_var_1917 - 1, 1 : 16) + fac101_var_1950 * absa_var_183(ind1_var_1917, 1 : 16) + fac001_var_1949 * absa_var_183(ind1_var_1917 + 1, 1 : 16) + fac211_var_1954 * absa_var_183(ind1_var_1917 + 8, 1 : 16) + fac111_var_1953 * absa_var_183(ind1_var_1917 + 9, 1 : 16) + fac011_var_1952 * absa_var_183(ind1_var_1917 + 10, 1 : 16))
      ELSE
        tau_major1_var_1966(1 : ng3) = speccomb1_var_1910 * (fac001_var_1949 * absa_var_183(ind1_var_1917, 1 : 16) + fac101_var_1950 * absa_var_183(ind1_var_1917 + 1, 1 : 16) + fac011_var_1952 * absa_var_183(ind1_var_1917 + 9, 1 : 16) + fac111_var_1953 * absa_var_183(ind1_var_1917 + 10, 1 : 16))
      END IF
      DO ig_var_1921 = 1, 16
        tauself_var_1961 = selffac_var_1899(jl_var_1977, lay_var_1923) * (selfref_var_185(inds_var_1918, ig_var_1921) + selffrac_var_1900(jl_var_1977, lay_var_1923) * (selfref_var_185(inds_var_1918 + 1, ig_var_1921) - selfref_var_185(inds_var_1918, ig_var_1921)))
        taufor_var_1960 = forfac_var_1889(jl_var_1977, lay_var_1923) * (forref_var_186(indf_var_1919, ig_var_1921) + forfrac_var_1906(jl_var_1977, lay_var_1923) * (forref_var_186(indf_var_1919 + 1, ig_var_1921) - forref_var_186(indf_var_1919, ig_var_1921)))
        n2om1_var_1962 = ka_mn2o_var_181(jmn2o_var_1925, indm_var_1920, ig_var_1921) + fmn2o_var_1933 * (ka_mn2o_var_181(jmn2o_var_1925 + 1, indm_var_1920, ig_var_1921) - ka_mn2o_var_181(jmn2o_var_1925, indm_var_1920, ig_var_1921))
        n2om2_var_1963 = ka_mn2o_var_181(jmn2o_var_1925, indm_var_1920 + 1, ig_var_1921) + fmn2o_var_1933 * (ka_mn2o_var_181(jmn2o_var_1925 + 1, indm_var_1920 + 1, ig_var_1921) - ka_mn2o_var_181(jmn2o_var_1925, indm_var_1920 + 1, ig_var_1921))
        absn2o_var_1964 = n2om1_var_1962 + minorfrac_var_1907(jl_var_1977, lay_var_1923) * (n2om2_var_1963 - n2om1_var_1962)
        taug_var_1883(jl_var_1977, 22 + ig_var_1921, lay_var_1923) = tau_major_var_1965(ig_var_1921) + tau_major1_var_1966(ig_var_1921) + tauself_var_1961 + taufor_var_1960 + adjcoln2o_var_1940 * absn2o_var_1964
        fracs_var_1902(jl_var_1977, 22 + ig_var_1921, lay_var_1923) = fracrefa_var_179(ig_var_1921, jpl_var_1926) + fpl_var_1936 * (fracrefa_var_179(ig_var_1921, jpl_var_1926 + 1) - fracrefa_var_179(ig_var_1921, jpl_var_1926))
      END DO
    END DO
  END DO
  DO lay_var_1923 = laytrop_max_var_1968 + 1, klev_var_1882
    DO jl_var_1977 = kidia_var_1880, kfdia_var_1881
      speccomb_var_1909 = colh2o_var_1894(jl_var_1977, lay_var_1923) + rat_h2oco2_var_1903(jl_var_1977, lay_var_1923) * colco2_var_1895(jl_var_1977, lay_var_1923)
      specparm_var_1929 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_var_1909, oneminus_var_1893)
      specmult_var_1928 = 4.0D0 * (specparm_var_1929)
      js_var_1922 = 1 + INT(specmult_var_1928)
      fs_var_1927 = ((specmult_var_1928) - AINT((specmult_var_1928)))
      speccomb1_var_1910 = colh2o_var_1894(jl_var_1977, lay_var_1923) + rat_h2oco2_1_var_1904(jl_var_1977, lay_var_1923) * colco2_var_1895(jl_var_1977, lay_var_1923)
      specparm1_var_1932 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb1_var_1910, oneminus_var_1893)
      specmult1_var_1931 = 4.0D0 * (specparm1_var_1932)
      js1_var_1924 = 1 + INT(specmult1_var_1931)
      fs1_var_1930 = ((specmult1_var_1931) - AINT((specmult1_var_1931)))
      fac000_var_1943 = (1.0D0 - fs_var_1927) * fac00_var_1885(jl_var_1977, lay_var_1923)
      fac010_var_1946 = (1.0D0 - fs_var_1927) * fac10_var_1887(jl_var_1977, lay_var_1923)
      fac100_var_1944 = fs_var_1927 * fac00_var_1885(jl_var_1977, lay_var_1923)
      fac110_var_1947 = fs_var_1927 * fac10_var_1887(jl_var_1977, lay_var_1923)
      fac001_var_1949 = (1.0D0 - fs1_var_1930) * fac01_var_1886(jl_var_1977, lay_var_1923)
      fac011_var_1952 = (1.0D0 - fs1_var_1930) * fac11_var_1888(jl_var_1977, lay_var_1923)
      fac101_var_1950 = fs1_var_1930 * fac01_var_1886(jl_var_1977, lay_var_1923)
      fac111_var_1953 = fs1_var_1930 * fac11_var_1888(jl_var_1977, lay_var_1923)
      speccomb_mn2o_var_1911 = colh2o_var_1894(jl_var_1977, lay_var_1923) + refrat_m_b * colco2_var_1895(jl_var_1977, lay_var_1923)
      specparm_mn2o_var_1935 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_mn2o_var_1911, oneminus_var_1893)
      specmult_mn2o_var_1934 = 4.0D0 * specparm_mn2o_var_1935
      jmn2o_var_1925 = 1 + INT(specmult_mn2o_var_1934)
      fmn2o_var_1933 = ((specmult_mn2o_var_1934) - AINT((specmult_mn2o_var_1934)))
      chi_n2o_var_1942 = coln2o_var_1896(jl_var_1977, lay_var_1923) / coldry_var_1897(jl_var_1977, lay_var_1923)
      ratn2o_var_1941 = 1D+20 * chi_n2o_var_1942 / chi_mls(4, jp_var_1890(jl_var_1977, lay_var_1923) + 1)
      IF (ratn2o_var_1941 .GT. 1.5D0) THEN
        adjfac_var_1939 = 0.5D0 + (ratn2o_var_1941 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1940 = adjfac_var_1939 * chi_mls(4, jp_var_1890(jl_var_1977, lay_var_1923) + 1) * coldry_var_1897(jl_var_1977, lay_var_1923) * 1D-20
      ELSE
        adjcoln2o_var_1940 = coln2o_var_1896(jl_var_1977, lay_var_1923)
      END IF
      speccomb_planck_var_1912 = colh2o_var_1894(jl_var_1977, lay_var_1923) + refrat_planck_b_var_1914 * colco2_var_1895(jl_var_1977, lay_var_1923)
      specparm_planck_var_1938 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_planck_var_1912, oneminus_var_1893)
      specmult_planck_var_1937 = 4.0D0 * specparm_planck_var_1938
      jpl_var_1926 = 1 + INT(specmult_planck_var_1937)
      fpl_var_1936 = ((specmult_planck_var_1937) - AINT((specmult_planck_var_1937)))
      ind0_var_1916 = ((jp_var_1890(jl_var_1977, lay_var_1923) - 13) * 5 + (jt_var_1891(jl_var_1977, lay_var_1923) - 1)) * nspb_var_237(3) + js_var_1922
      ind1_var_1917 = ((jp_var_1890(jl_var_1977, lay_var_1923) - 12) * 5 + (jt1_var_1892(jl_var_1977, lay_var_1923) - 1)) * nspb_var_237(3) + js1_var_1924
      indf_var_1919 = indfor_var_1905(jl_var_1977, lay_var_1923)
      indm_var_1920 = indminor_var_1908(jl_var_1977, lay_var_1923)
      DO ig_var_1921 = 1, 16
        taufor_var_1960 = forfac_var_1889(jl_var_1977, lay_var_1923) * (forref_var_186(indf_var_1919, ig_var_1921) + forfrac_var_1906(jl_var_1977, lay_var_1923) * (forref_var_186(indf_var_1919 + 1, ig_var_1921) - forref_var_186(indf_var_1919, ig_var_1921)))
        n2om1_var_1962 = kb_mn2o_var_182(jmn2o_var_1925, indm_var_1920, ig_var_1921) + fmn2o_var_1933 * (kb_mn2o_var_182(jmn2o_var_1925 + 1, indm_var_1920, ig_var_1921) - kb_mn2o_var_182(jmn2o_var_1925, indm_var_1920, ig_var_1921))
        n2om2_var_1963 = kb_mn2o_var_182(jmn2o_var_1925, indm_var_1920 + 1, ig_var_1921) + fmn2o_var_1933 * (kb_mn2o_var_182(jmn2o_var_1925 + 1, indm_var_1920 + 1, ig_var_1921) - kb_mn2o_var_182(jmn2o_var_1925, indm_var_1920 + 1, ig_var_1921))
        absn2o_var_1964 = n2om1_var_1962 + minorfrac_var_1907(jl_var_1977, lay_var_1923) * (n2om2_var_1963 - n2om1_var_1962)
        taug_var_1883(jl_var_1977, 22 + ig_var_1921, lay_var_1923) = speccomb_var_1909 * (fac000_var_1943 * absb_var_184(ind0_var_1916, ig_var_1921) + fac100_var_1944 * absb_var_184(ind0_var_1916 + 1, ig_var_1921) + fac010_var_1946 * absb_var_184(ind0_var_1916 + 5, ig_var_1921) + fac110_var_1947 * absb_var_184(ind0_var_1916 + 6, ig_var_1921)) + speccomb1_var_1910 * (fac001_var_1949 * absb_var_184(ind1_var_1917, ig_var_1921) + fac101_var_1950 * absb_var_184(ind1_var_1917 + 1, ig_var_1921) + fac011_var_1952 * absb_var_184(ind1_var_1917 + 5, ig_var_1921) + fac111_var_1953 * absb_var_184(ind1_var_1917 + 6, ig_var_1921)) + taufor_var_1960 + adjcoln2o_var_1940 * absn2o_var_1964
        fracs_var_1902(jl_var_1977, 22 + ig_var_1921, lay_var_1923) = fracrefb_var_180(ig_var_1921, jpl_var_1926) + fpl_var_1936 * (fracrefb_var_180(ig_var_1921, jpl_var_1926 + 1) - fracrefb_var_180(ig_var_1921, jpl_var_1926))
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1968 /= laytrop_min_var_1967) THEN
    DO lay_var_1923 = laytrop_min_var_1967 + 1, laytrop_max_var_1968
      ixc0_var_1974 = ixc_var_1969(lay_var_1923)
      DO ixp_var_1975 = 1, ixc0_var_1974
        jl_var_1977 = ixlow_var_1970(ixp_var_1975, lay_var_1923)
        speccomb_var_1909 = colh2o_var_1894(jl_var_1977, lay_var_1923) + rat_h2oco2_var_1903(jl_var_1977, lay_var_1923) * colco2_var_1895(jl_var_1977, lay_var_1923)
        specparm_var_1929 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_var_1909, oneminus_var_1893)
        specmult_var_1928 = 8.0D0 * (specparm_var_1929)
        js_var_1922 = 1 + INT(specmult_var_1928)
        fs_var_1927 = ((specmult_var_1928) - AINT((specmult_var_1928)))
        speccomb1_var_1910 = colh2o_var_1894(jl_var_1977, lay_var_1923) + rat_h2oco2_1_var_1904(jl_var_1977, lay_var_1923) * colco2_var_1895(jl_var_1977, lay_var_1923)
        specparm1_var_1932 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb1_var_1910, oneminus_var_1893)
        specmult1_var_1931 = 8.0D0 * (specparm1_var_1932)
        js1_var_1924 = 1 + INT(specmult1_var_1931)
        fs1_var_1930 = ((specmult1_var_1931) - AINT((specmult1_var_1931)))
        speccomb_mn2o_var_1911 = colh2o_var_1894(jl_var_1977, lay_var_1923) + refrat_m_a_var_1915 * colco2_var_1895(jl_var_1977, lay_var_1923)
        specparm_mn2o_var_1935 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_mn2o_var_1911, oneminus_var_1893)
        specmult_mn2o_var_1934 = 8.0D0 * specparm_mn2o_var_1935
        jmn2o_var_1925 = 1 + INT(specmult_mn2o_var_1934)
        fmn2o_var_1933 = ((specmult_mn2o_var_1934) - AINT((specmult_mn2o_var_1934)))
        chi_n2o_var_1942 = coln2o_var_1896(jl_var_1977, lay_var_1923) / coldry_var_1897(jl_var_1977, lay_var_1923)
        ratn2o_var_1941 = 1D+20 * chi_n2o_var_1942 / chi_mls(4, jp_var_1890(jl_var_1977, lay_var_1923) + 1)
        IF (ratn2o_var_1941 .GT. 1.5D0) THEN
          adjfac_var_1939 = 0.5D0 + (ratn2o_var_1941 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1940 = adjfac_var_1939 * chi_mls(4, jp_var_1890(jl_var_1977, lay_var_1923) + 1) * coldry_var_1897(jl_var_1977, lay_var_1923) * 1D-20
        ELSE
          adjcoln2o_var_1940 = coln2o_var_1896(jl_var_1977, lay_var_1923)
        END IF
        speccomb_planck_var_1912 = colh2o_var_1894(jl_var_1977, lay_var_1923) + refrat_planck_a_var_1913 * colco2_var_1895(jl_var_1977, lay_var_1923)
        specparm_planck_var_1938 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_planck_var_1912, oneminus_var_1893)
        specmult_planck_var_1937 = 8.0D0 * specparm_planck_var_1938
        jpl_var_1926 = 1 + INT(specmult_planck_var_1937)
        fpl_var_1936 = ((specmult_planck_var_1937) - AINT((specmult_planck_var_1937)))
        ind0_var_1916 = ((jp_var_1890(jl_var_1977, lay_var_1923) - 1) * 5 + (jt_var_1891(jl_var_1977, lay_var_1923) - 1)) * nspa_var_236(3) + js_var_1922
        ind1_var_1917 = (jp_var_1890(jl_var_1977, lay_var_1923) * 5 + (jt1_var_1892(jl_var_1977, lay_var_1923) - 1)) * nspa_var_236(3) + js1_var_1924
        inds_var_1918 = indself_var_1901(jl_var_1977, lay_var_1923)
        indf_var_1919 = indfor_var_1905(jl_var_1977, lay_var_1923)
        indm_var_1920 = indminor_var_1908(jl_var_1977, lay_var_1923)
        IF (specparm_var_1929 .LT. 0.125D0) THEN
          p_var_1955 = fs_var_1927 - 1.0D0
          p4_var_1956 = p_var_1955 ** 4
          fk0_var_1957 = p4_var_1956
          fk1_var_1958 = 1.0D0 - p_var_1955 - 2.0D0 * p4_var_1956
          fk2_var_1959 = p_var_1955 + p4_var_1956
          fac000_var_1943 = fk0_var_1957 * fac00_var_1885(jl_var_1977, lay_var_1923)
          fac100_var_1944 = fk1_var_1958 * fac00_var_1885(jl_var_1977, lay_var_1923)
          fac200_var_1945 = fk2_var_1959 * fac00_var_1885(jl_var_1977, lay_var_1923)
          fac010_var_1946 = fk0_var_1957 * fac10_var_1887(jl_var_1977, lay_var_1923)
          fac110_var_1947 = fk1_var_1958 * fac10_var_1887(jl_var_1977, lay_var_1923)
          fac210_var_1948 = fk2_var_1959 * fac10_var_1887(jl_var_1977, lay_var_1923)
        ELSE IF (specparm_var_1929 .GT. 0.875D0) THEN
          p_var_1955 = - fs_var_1927
          p4_var_1956 = p_var_1955 ** 4
          fk0_var_1957 = p4_var_1956
          fk1_var_1958 = 1.0D0 - p_var_1955 - 2.0D0 * p4_var_1956
          fk2_var_1959 = p_var_1955 + p4_var_1956
          fac000_var_1943 = fk0_var_1957 * fac00_var_1885(jl_var_1977, lay_var_1923)
          fac100_var_1944 = fk1_var_1958 * fac00_var_1885(jl_var_1977, lay_var_1923)
          fac200_var_1945 = fk2_var_1959 * fac00_var_1885(jl_var_1977, lay_var_1923)
          fac010_var_1946 = fk0_var_1957 * fac10_var_1887(jl_var_1977, lay_var_1923)
          fac110_var_1947 = fk1_var_1958 * fac10_var_1887(jl_var_1977, lay_var_1923)
          fac210_var_1948 = fk2_var_1959 * fac10_var_1887(jl_var_1977, lay_var_1923)
        ELSE
          fac000_var_1943 = (1.0D0 - fs_var_1927) * fac00_var_1885(jl_var_1977, lay_var_1923)
          fac010_var_1946 = (1.0D0 - fs_var_1927) * fac10_var_1887(jl_var_1977, lay_var_1923)
          fac100_var_1944 = fs_var_1927 * fac00_var_1885(jl_var_1977, lay_var_1923)
          fac110_var_1947 = fs_var_1927 * fac10_var_1887(jl_var_1977, lay_var_1923)
          fac200_var_1945 = 0.0D0
          fac210_var_1948 = 0.0D0
        END IF
        IF (specparm1_var_1932 .LT. 0.125D0) THEN
          p_var_1955 = fs1_var_1930 - 1.0D0
          p4_var_1956 = p_var_1955 ** 4
          fk0_var_1957 = p4_var_1956
          fk1_var_1958 = 1.0D0 - p_var_1955 - 2.0D0 * p4_var_1956
          fk2_var_1959 = p_var_1955 + p4_var_1956
          fac001_var_1949 = fk0_var_1957 * fac01_var_1886(jl_var_1977, lay_var_1923)
          fac101_var_1950 = fk1_var_1958 * fac01_var_1886(jl_var_1977, lay_var_1923)
          fac201_var_1951 = fk2_var_1959 * fac01_var_1886(jl_var_1977, lay_var_1923)
          fac011_var_1952 = fk0_var_1957 * fac11_var_1888(jl_var_1977, lay_var_1923)
          fac111_var_1953 = fk1_var_1958 * fac11_var_1888(jl_var_1977, lay_var_1923)
          fac211_var_1954 = fk2_var_1959 * fac11_var_1888(jl_var_1977, lay_var_1923)
        ELSE IF (specparm1_var_1932 .GT. 0.875D0) THEN
          p_var_1955 = - fs1_var_1930
          p4_var_1956 = p_var_1955 ** 4
          fk0_var_1957 = p4_var_1956
          fk1_var_1958 = 1.0D0 - p_var_1955 - 2.0D0 * p4_var_1956
          fk2_var_1959 = p_var_1955 + p4_var_1956
          fac001_var_1949 = fk0_var_1957 * fac01_var_1886(jl_var_1977, lay_var_1923)
          fac101_var_1950 = fk1_var_1958 * fac01_var_1886(jl_var_1977, lay_var_1923)
          fac201_var_1951 = fk2_var_1959 * fac01_var_1886(jl_var_1977, lay_var_1923)
          fac011_var_1952 = fk0_var_1957 * fac11_var_1888(jl_var_1977, lay_var_1923)
          fac111_var_1953 = fk1_var_1958 * fac11_var_1888(jl_var_1977, lay_var_1923)
          fac211_var_1954 = fk2_var_1959 * fac11_var_1888(jl_var_1977, lay_var_1923)
        ELSE
          fac001_var_1949 = (1.0D0 - fs1_var_1930) * fac01_var_1886(jl_var_1977, lay_var_1923)
          fac011_var_1952 = (1.0D0 - fs1_var_1930) * fac11_var_1888(jl_var_1977, lay_var_1923)
          fac101_var_1950 = fs1_var_1930 * fac01_var_1886(jl_var_1977, lay_var_1923)
          fac111_var_1953 = fs1_var_1930 * fac11_var_1888(jl_var_1977, lay_var_1923)
          fac201_var_1951 = 0.0D0
          fac211_var_1954 = 0.0D0
        END IF
        IF (specparm_var_1929 .LT. 0.125D0) THEN
          tau_major_var_1965(1 : ng3) = speccomb_var_1909 * (fac000_var_1943 * absa_var_183(ind0_var_1916, 1 : 16) + fac100_var_1944 * absa_var_183(ind0_var_1916 + 1, 1 : 16) + fac200_var_1945 * absa_var_183(ind0_var_1916 + 2, 1 : 16) + fac010_var_1946 * absa_var_183(ind0_var_1916 + 9, 1 : 16) + fac110_var_1947 * absa_var_183(ind0_var_1916 + 10, 1 : 16) + fac210_var_1948 * absa_var_183(ind0_var_1916 + 11, 1 : 16))
        ELSE IF (specparm_var_1929 .GT. 0.875D0) THEN
          tau_major_var_1965(1 : ng3) = speccomb_var_1909 * (fac200_var_1945 * absa_var_183(ind0_var_1916 - 1, 1 : 16) + fac100_var_1944 * absa_var_183(ind0_var_1916, 1 : 16) + fac000_var_1943 * absa_var_183(ind0_var_1916 + 1, 1 : 16) + fac210_var_1948 * absa_var_183(ind0_var_1916 + 8, 1 : 16) + fac110_var_1947 * absa_var_183(ind0_var_1916 + 9, 1 : 16) + fac010_var_1946 * absa_var_183(ind0_var_1916 + 10, 1 : 16))
        ELSE
          tau_major_var_1965(1 : ng3) = speccomb_var_1909 * (fac000_var_1943 * absa_var_183(ind0_var_1916, 1 : 16) + fac100_var_1944 * absa_var_183(ind0_var_1916 + 1, 1 : 16) + fac010_var_1946 * absa_var_183(ind0_var_1916 + 9, 1 : 16) + fac110_var_1947 * absa_var_183(ind0_var_1916 + 10, 1 : 16))
        END IF
        IF (specparm1_var_1932 .LT. 0.125D0) THEN
          tau_major1_var_1966(1 : ng3) = speccomb1_var_1910 * (fac001_var_1949 * absa_var_183(ind1_var_1917, 1 : 16) + fac101_var_1950 * absa_var_183(ind1_var_1917 + 1, 1 : 16) + fac201_var_1951 * absa_var_183(ind1_var_1917 + 2, 1 : 16) + fac011_var_1952 * absa_var_183(ind1_var_1917 + 9, 1 : 16) + fac111_var_1953 * absa_var_183(ind1_var_1917 + 10, 1 : 16) + fac211_var_1954 * absa_var_183(ind1_var_1917 + 11, 1 : 16))
        ELSE IF (specparm1_var_1932 .GT. 0.875D0) THEN
          tau_major1_var_1966(1 : ng3) = speccomb1_var_1910 * (fac201_var_1951 * absa_var_183(ind1_var_1917 - 1, 1 : 16) + fac101_var_1950 * absa_var_183(ind1_var_1917, 1 : 16) + fac001_var_1949 * absa_var_183(ind1_var_1917 + 1, 1 : 16) + fac211_var_1954 * absa_var_183(ind1_var_1917 + 8, 1 : 16) + fac111_var_1953 * absa_var_183(ind1_var_1917 + 9, 1 : 16) + fac011_var_1952 * absa_var_183(ind1_var_1917 + 10, 1 : 16))
        ELSE
          tau_major1_var_1966(1 : ng3) = speccomb1_var_1910 * (fac001_var_1949 * absa_var_183(ind1_var_1917, 1 : 16) + fac101_var_1950 * absa_var_183(ind1_var_1917 + 1, 1 : 16) + fac011_var_1952 * absa_var_183(ind1_var_1917 + 9, 1 : 16) + fac111_var_1953 * absa_var_183(ind1_var_1917 + 10, 1 : 16))
        END IF
        DO ig_var_1921 = 1, 16
          tauself_var_1961 = selffac_var_1899(jl_var_1977, lay_var_1923) * (selfref_var_185(inds_var_1918, ig_var_1921) + selffrac_var_1900(jl_var_1977, lay_var_1923) * (selfref_var_185(inds_var_1918 + 1, ig_var_1921) - selfref_var_185(inds_var_1918, ig_var_1921)))
          taufor_var_1960 = forfac_var_1889(jl_var_1977, lay_var_1923) * (forref_var_186(indf_var_1919, ig_var_1921) + forfrac_var_1906(jl_var_1977, lay_var_1923) * (forref_var_186(indf_var_1919 + 1, ig_var_1921) - forref_var_186(indf_var_1919, ig_var_1921)))
          n2om1_var_1962 = ka_mn2o_var_181(jmn2o_var_1925, indm_var_1920, ig_var_1921) + fmn2o_var_1933 * (ka_mn2o_var_181(jmn2o_var_1925 + 1, indm_var_1920, ig_var_1921) - ka_mn2o_var_181(jmn2o_var_1925, indm_var_1920, ig_var_1921))
          n2om2_var_1963 = ka_mn2o_var_181(jmn2o_var_1925, indm_var_1920 + 1, ig_var_1921) + fmn2o_var_1933 * (ka_mn2o_var_181(jmn2o_var_1925 + 1, indm_var_1920 + 1, ig_var_1921) - ka_mn2o_var_181(jmn2o_var_1925, indm_var_1920 + 1, ig_var_1921))
          absn2o_var_1964 = n2om1_var_1962 + minorfrac_var_1907(jl_var_1977, lay_var_1923) * (n2om2_var_1963 - n2om1_var_1962)
          taug_var_1883(jl_var_1977, 22 + ig_var_1921, lay_var_1923) = tau_major_var_1965(ig_var_1921) + tau_major1_var_1966(ig_var_1921) + tauself_var_1961 + taufor_var_1960 + adjcoln2o_var_1940 * absn2o_var_1964
          fracs_var_1902(jl_var_1977, 22 + ig_var_1921, lay_var_1923) = fracrefa_var_179(ig_var_1921, jpl_var_1926) + fpl_var_1936 * (fracrefa_var_179(ig_var_1921, jpl_var_1926 + 1) - fracrefa_var_179(ig_var_1921, jpl_var_1926))
        END DO
      END DO
      ixc0_var_1974 = kfdia_var_1881 - kidia_var_1880 + 1 - ixc0_var_1974
      DO ixp_var_1975 = 1, ixc0_var_1974
        jl_var_1977 = ixhigh_var_1971(ixp_var_1975, lay_var_1923)
        speccomb_var_1909 = colh2o_var_1894(jl_var_1977, lay_var_1923) + rat_h2oco2_var_1903(jl_var_1977, lay_var_1923) * colco2_var_1895(jl_var_1977, lay_var_1923)
        specparm_var_1929 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_var_1909, oneminus_var_1893)
        specmult_var_1928 = 4.0D0 * (specparm_var_1929)
        js_var_1922 = 1 + INT(specmult_var_1928)
        fs_var_1927 = ((specmult_var_1928) - AINT((specmult_var_1928)))
        speccomb1_var_1910 = colh2o_var_1894(jl_var_1977, lay_var_1923) + rat_h2oco2_1_var_1904(jl_var_1977, lay_var_1923) * colco2_var_1895(jl_var_1977, lay_var_1923)
        specparm1_var_1932 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb1_var_1910, oneminus_var_1893)
        specmult1_var_1931 = 4.0D0 * (specparm1_var_1932)
        js1_var_1924 = 1 + INT(specmult1_var_1931)
        fs1_var_1930 = ((specmult1_var_1931) - AINT((specmult1_var_1931)))
        fac000_var_1943 = (1.0D0 - fs_var_1927) * fac00_var_1885(jl_var_1977, lay_var_1923)
        fac010_var_1946 = (1.0D0 - fs_var_1927) * fac10_var_1887(jl_var_1977, lay_var_1923)
        fac100_var_1944 = fs_var_1927 * fac00_var_1885(jl_var_1977, lay_var_1923)
        fac110_var_1947 = fs_var_1927 * fac10_var_1887(jl_var_1977, lay_var_1923)
        fac001_var_1949 = (1.0D0 - fs1_var_1930) * fac01_var_1886(jl_var_1977, lay_var_1923)
        fac011_var_1952 = (1.0D0 - fs1_var_1930) * fac11_var_1888(jl_var_1977, lay_var_1923)
        fac101_var_1950 = fs1_var_1930 * fac01_var_1886(jl_var_1977, lay_var_1923)
        fac111_var_1953 = fs1_var_1930 * fac11_var_1888(jl_var_1977, lay_var_1923)
        speccomb_mn2o_var_1911 = colh2o_var_1894(jl_var_1977, lay_var_1923) + refrat_m_b * colco2_var_1895(jl_var_1977, lay_var_1923)
        specparm_mn2o_var_1935 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_mn2o_var_1911, oneminus_var_1893)
        specmult_mn2o_var_1934 = 4.0D0 * specparm_mn2o_var_1935
        jmn2o_var_1925 = 1 + INT(specmult_mn2o_var_1934)
        fmn2o_var_1933 = ((specmult_mn2o_var_1934) - AINT((specmult_mn2o_var_1934)))
        chi_n2o_var_1942 = coln2o_var_1896(jl_var_1977, lay_var_1923) / coldry_var_1897(jl_var_1977, lay_var_1923)
        ratn2o_var_1941 = 1D+20 * chi_n2o_var_1942 / chi_mls(4, jp_var_1890(jl_var_1977, lay_var_1923) + 1)
        IF (ratn2o_var_1941 .GT. 1.5D0) THEN
          adjfac_var_1939 = 0.5D0 + (ratn2o_var_1941 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1940 = adjfac_var_1939 * chi_mls(4, jp_var_1890(jl_var_1977, lay_var_1923) + 1) * coldry_var_1897(jl_var_1977, lay_var_1923) * 1D-20
        ELSE
          adjcoln2o_var_1940 = coln2o_var_1896(jl_var_1977, lay_var_1923)
        END IF
        speccomb_planck_var_1912 = colh2o_var_1894(jl_var_1977, lay_var_1923) + refrat_planck_b_var_1914 * colco2_var_1895(jl_var_1977, lay_var_1923)
        specparm_planck_var_1938 = MIN(colh2o_var_1894(jl_var_1977, lay_var_1923) / speccomb_planck_var_1912, oneminus_var_1893)
        specmult_planck_var_1937 = 4.0D0 * specparm_planck_var_1938
        jpl_var_1926 = 1 + INT(specmult_planck_var_1937)
        fpl_var_1936 = ((specmult_planck_var_1937) - AINT((specmult_planck_var_1937)))
        ind0_var_1916 = ((jp_var_1890(jl_var_1977, lay_var_1923) - 13) * 5 + (jt_var_1891(jl_var_1977, lay_var_1923) - 1)) * nspb_var_237(3) + js_var_1922
        ind1_var_1917 = ((jp_var_1890(jl_var_1977, lay_var_1923) - 12) * 5 + (jt1_var_1892(jl_var_1977, lay_var_1923) - 1)) * nspb_var_237(3) + js1_var_1924
        indf_var_1919 = indfor_var_1905(jl_var_1977, lay_var_1923)
        indm_var_1920 = indminor_var_1908(jl_var_1977, lay_var_1923)
        DO ig_var_1921 = 1, 16
          taufor_var_1960 = forfac_var_1889(jl_var_1977, lay_var_1923) * (forref_var_186(indf_var_1919, ig_var_1921) + forfrac_var_1906(jl_var_1977, lay_var_1923) * (forref_var_186(indf_var_1919 + 1, ig_var_1921) - forref_var_186(indf_var_1919, ig_var_1921)))
          n2om1_var_1962 = kb_mn2o_var_182(jmn2o_var_1925, indm_var_1920, ig_var_1921) + fmn2o_var_1933 * (kb_mn2o_var_182(jmn2o_var_1925 + 1, indm_var_1920, ig_var_1921) - kb_mn2o_var_182(jmn2o_var_1925, indm_var_1920, ig_var_1921))
          n2om2_var_1963 = kb_mn2o_var_182(jmn2o_var_1925, indm_var_1920 + 1, ig_var_1921) + fmn2o_var_1933 * (kb_mn2o_var_182(jmn2o_var_1925 + 1, indm_var_1920 + 1, ig_var_1921) - kb_mn2o_var_182(jmn2o_var_1925, indm_var_1920 + 1, ig_var_1921))
          absn2o_var_1964 = n2om1_var_1962 + minorfrac_var_1907(jl_var_1977, lay_var_1923) * (n2om2_var_1963 - n2om1_var_1962)
          taug_var_1883(jl_var_1977, 22 + ig_var_1921, lay_var_1923) = speccomb_var_1909 * (fac000_var_1943 * absb_var_184(ind0_var_1916, ig_var_1921) + fac100_var_1944 * absb_var_184(ind0_var_1916 + 1, ig_var_1921) + fac010_var_1946 * absb_var_184(ind0_var_1916 + 5, ig_var_1921) + fac110_var_1947 * absb_var_184(ind0_var_1916 + 6, ig_var_1921)) + speccomb1_var_1910 * (fac001_var_1949 * absb_var_184(ind1_var_1917, ig_var_1921) + fac101_var_1950 * absb_var_184(ind1_var_1917 + 1, ig_var_1921) + fac011_var_1952 * absb_var_184(ind1_var_1917 + 5, ig_var_1921) + fac111_var_1953 * absb_var_184(ind1_var_1917 + 6, ig_var_1921)) + taufor_var_1960 + adjcoln2o_var_1940 * absn2o_var_1964
          fracs_var_1902(jl_var_1977, 22 + ig_var_1921, lay_var_1923) = fracrefb_var_180(ig_var_1921, jpl_var_1926) + fpl_var_1936 * (fracrefb_var_180(ig_var_1921, jpl_var_1926 + 1) - fracrefb_var_180(ig_var_1921, jpl_var_1926))
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol3
SUBROUTINE rrtm_taumol2(kidia_var_1978, kfdia_var_1979, klev_var_1980, taug_var_1981, pavel_var_1982, coldry_var_1983, p_tauaerl_var_1984, fac00_var_1985, fac01_var_1986, fac10_var_1987, fac11_var_1988, forfac_var_1990, forfrac_var_1989, indfor_var_1999, jp_var_1991, jt_var_1992, jt1_var_1993, colh2o_var_1994, laytrop_var_1995, selffac_var_1996, selffrac_var_1997, indself_var_1998, fracs_var_2000)
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrta2, ONLY: absa_var_175, absb_var_176, forref_var_178, fracrefa_var_173, fracrefb_var_174, selfref_var_177
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1978
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1979
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1980
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1981(kidia_var_1978 : kfdia_var_1979, 140, klev_var_1980)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1982(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1983(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1984(kidia_var_1978 : kfdia_var_1979, klev_var_1980, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1985(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1986(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1987(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1988(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1989(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1990(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1991(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1992(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1993(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1994(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1995(kidia_var_1978 : kfdia_var_1979)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1996(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1997(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1998(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1999(kidia_var_1978 : kfdia_var_1979, klev_var_1980)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2000(kidia_var_1978 : kfdia_var_1979, 140, klev_var_1980)
  INTEGER(KIND = 4) :: ind0_var_2001, ind1_var_2002, inds_var_2003, indf_var_2004
  INTEGER(KIND = 4) :: ig_var_2005, lay_var_2006
  REAL(KIND = 8) :: taufor_var_2007, tauself_var_2008, corradj_var_2009, pp_var_2010
  INTEGER(KIND = 4) :: laytrop_min_var_2011, laytrop_max_var_2012
  INTEGER(KIND = 4) :: ixc_var_2013(klev_var_1980), ixlow_var_2014(kfdia_var_1979, klev_var_1980), ixhigh_var_2015(kfdia_var_1979, klev_var_1980)
  INTEGER(KIND = 4) :: ich_var_2016, icl_var_2017, ixc0_var_2018, ixp_var_2019, jc_var_2020, jl_var_2021
  laytrop_min_var_2011 = MINVAL(laytrop_var_1995)
  laytrop_max_var_2012 = MAXVAL(laytrop_var_1995)
  ixlow_var_2014 = 0
  ixhigh_var_2015 = 0
  ixc_var_2013 = 0
  DO lay_var_2006 = laytrop_min_var_2011 + 1, laytrop_max_var_2012
    icl_var_2017 = 0
    ich_var_2016 = 0
    DO jc_var_2020 = kidia_var_1978, kfdia_var_1979
      IF (lay_var_2006 <= laytrop_var_1995(jc_var_2020)) THEN
        icl_var_2017 = icl_var_2017 + 1
        ixlow_var_2014(icl_var_2017, lay_var_2006) = jc_var_2020
      ELSE
        ich_var_2016 = ich_var_2016 + 1
        ixhigh_var_2015(ich_var_2016, lay_var_2006) = jc_var_2020
      END IF
    END DO
    ixc_var_2013(lay_var_2006) = icl_var_2017
  END DO
  DO lay_var_2006 = 1, laytrop_min_var_2011
    DO jl_var_2021 = kidia_var_1978, kfdia_var_1979
      ind0_var_2001 = ((jp_var_1991(jl_var_2021, lay_var_2006) - 1) * 5 + (jt_var_1992(jl_var_2021, lay_var_2006) - 1)) * nspa_var_236(2) + 1
      ind1_var_2002 = (jp_var_1991(jl_var_2021, lay_var_2006) * 5 + (jt1_var_1993(jl_var_2021, lay_var_2006) - 1)) * nspa_var_236(2) + 1
      inds_var_2003 = indself_var_1998(jl_var_2021, lay_var_2006)
      indf_var_2004 = indfor_var_1999(jl_var_2021, lay_var_2006)
      pp_var_2010 = pavel_var_1982(jl_var_2021, lay_var_2006)
      corradj_var_2009 = 1.0D0 - 0.05D0 * (pp_var_2010 - 100.0D0) / 900.0D0
      DO ig_var_2005 = 1, 12
        tauself_var_2008 = selffac_var_1996(jl_var_2021, lay_var_2006) * (selfref_var_177(inds_var_2003, ig_var_2005) + selffrac_var_1997(jl_var_2021, lay_var_2006) * (selfref_var_177(inds_var_2003 + 1, ig_var_2005) - selfref_var_177(inds_var_2003, ig_var_2005)))
        taufor_var_2007 = forfac_var_1990(jl_var_2021, lay_var_2006) * (forref_var_178(indf_var_2004, ig_var_2005) + forfrac_var_1989(jl_var_2021, lay_var_2006) * (forref_var_178(indf_var_2004 + 1, ig_var_2005) - forref_var_178(indf_var_2004, ig_var_2005)))
        taug_var_1981(jl_var_2021, 10 + ig_var_2005, lay_var_2006) = corradj_var_2009 * (colh2o_var_1994(jl_var_2021, lay_var_2006) * (fac00_var_1985(jl_var_2021, lay_var_2006) * absa_var_175(ind0_var_2001, ig_var_2005) + fac10_var_1987(jl_var_2021, lay_var_2006) * absa_var_175(ind0_var_2001 + 1, ig_var_2005) + fac01_var_1986(jl_var_2021, lay_var_2006) * absa_var_175(ind1_var_2002, ig_var_2005) + fac11_var_1988(jl_var_2021, lay_var_2006) * absa_var_175(ind1_var_2002 + 1, ig_var_2005)) + tauself_var_2008 + taufor_var_2007)
        fracs_var_2000(jl_var_2021, 10 + ig_var_2005, lay_var_2006) = fracrefa_var_173(ig_var_2005)
      END DO
    END DO
  END DO
  DO lay_var_2006 = laytrop_max_var_2012 + 1, klev_var_1980
    DO jl_var_2021 = kidia_var_1978, kfdia_var_1979
      ind0_var_2001 = ((jp_var_1991(jl_var_2021, lay_var_2006) - 13) * 5 + (jt_var_1992(jl_var_2021, lay_var_2006) - 1)) * nspb_var_237(2) + 1
      ind1_var_2002 = ((jp_var_1991(jl_var_2021, lay_var_2006) - 12) * 5 + (jt1_var_1993(jl_var_2021, lay_var_2006) - 1)) * nspb_var_237(2) + 1
      indf_var_2004 = indfor_var_1999(jl_var_2021, lay_var_2006)
      DO ig_var_2005 = 1, 12
        taufor_var_2007 = forfac_var_1990(jl_var_2021, lay_var_2006) * (forref_var_178(indf_var_2004, ig_var_2005) + forfrac_var_1989(jl_var_2021, lay_var_2006) * (forref_var_178(indf_var_2004 + 1, ig_var_2005) - forref_var_178(indf_var_2004, ig_var_2005)))
        taug_var_1981(jl_var_2021, 10 + ig_var_2005, lay_var_2006) = colh2o_var_1994(jl_var_2021, lay_var_2006) * (fac00_var_1985(jl_var_2021, lay_var_2006) * absb_var_176(ind0_var_2001, ig_var_2005) + fac10_var_1987(jl_var_2021, lay_var_2006) * absb_var_176(ind0_var_2001 + 1, ig_var_2005) + fac01_var_1986(jl_var_2021, lay_var_2006) * absb_var_176(ind1_var_2002, ig_var_2005) + fac11_var_1988(jl_var_2021, lay_var_2006) * absb_var_176(ind1_var_2002 + 1, ig_var_2005)) + taufor_var_2007
        fracs_var_2000(jl_var_2021, 10 + ig_var_2005, lay_var_2006) = fracrefb_var_174(ig_var_2005)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2012 /= laytrop_min_var_2011) THEN
    DO lay_var_2006 = laytrop_min_var_2011 + 1, laytrop_max_var_2012
      ixc0_var_2018 = ixc_var_2013(lay_var_2006)
      DO ixp_var_2019 = 1, ixc0_var_2018
        jl_var_2021 = ixlow_var_2014(ixp_var_2019, lay_var_2006)
        ind0_var_2001 = ((jp_var_1991(jl_var_2021, lay_var_2006) - 1) * 5 + (jt_var_1992(jl_var_2021, lay_var_2006) - 1)) * nspa_var_236(2) + 1
        ind1_var_2002 = (jp_var_1991(jl_var_2021, lay_var_2006) * 5 + (jt1_var_1993(jl_var_2021, lay_var_2006) - 1)) * nspa_var_236(2) + 1
        inds_var_2003 = indself_var_1998(jl_var_2021, lay_var_2006)
        indf_var_2004 = indfor_var_1999(jl_var_2021, lay_var_2006)
        pp_var_2010 = pavel_var_1982(jl_var_2021, lay_var_2006)
        corradj_var_2009 = 1.0D0 - 0.05D0 * (pp_var_2010 - 100.0D0) / 900.0D0
        DO ig_var_2005 = 1, 12
          tauself_var_2008 = selffac_var_1996(jl_var_2021, lay_var_2006) * (selfref_var_177(inds_var_2003, ig_var_2005) + selffrac_var_1997(jl_var_2021, lay_var_2006) * (selfref_var_177(inds_var_2003 + 1, ig_var_2005) - selfref_var_177(inds_var_2003, ig_var_2005)))
          taufor_var_2007 = forfac_var_1990(jl_var_2021, lay_var_2006) * (forref_var_178(indf_var_2004, ig_var_2005) + forfrac_var_1989(jl_var_2021, lay_var_2006) * (forref_var_178(indf_var_2004 + 1, ig_var_2005) - forref_var_178(indf_var_2004, ig_var_2005)))
          taug_var_1981(jl_var_2021, 10 + ig_var_2005, lay_var_2006) = corradj_var_2009 * (colh2o_var_1994(jl_var_2021, lay_var_2006) * (fac00_var_1985(jl_var_2021, lay_var_2006) * absa_var_175(ind0_var_2001, ig_var_2005) + fac10_var_1987(jl_var_2021, lay_var_2006) * absa_var_175(ind0_var_2001 + 1, ig_var_2005) + fac01_var_1986(jl_var_2021, lay_var_2006) * absa_var_175(ind1_var_2002, ig_var_2005) + fac11_var_1988(jl_var_2021, lay_var_2006) * absa_var_175(ind1_var_2002 + 1, ig_var_2005)) + tauself_var_2008 + taufor_var_2007)
          fracs_var_2000(jl_var_2021, 10 + ig_var_2005, lay_var_2006) = fracrefa_var_173(ig_var_2005)
        END DO
      END DO
      ixc0_var_2018 = kfdia_var_1979 - kidia_var_1978 + 1 - ixc0_var_2018
      DO ixp_var_2019 = 1, ixc0_var_2018
        jl_var_2021 = ixhigh_var_2015(ixp_var_2019, lay_var_2006)
        ind0_var_2001 = ((jp_var_1991(jl_var_2021, lay_var_2006) - 13) * 5 + (jt_var_1992(jl_var_2021, lay_var_2006) - 1)) * nspb_var_237(2) + 1
        ind1_var_2002 = ((jp_var_1991(jl_var_2021, lay_var_2006) - 12) * 5 + (jt1_var_1993(jl_var_2021, lay_var_2006) - 1)) * nspb_var_237(2) + 1
        indf_var_2004 = indfor_var_1999(jl_var_2021, lay_var_2006)
        DO ig_var_2005 = 1, 12
          taufor_var_2007 = forfac_var_1990(jl_var_2021, lay_var_2006) * (forref_var_178(indf_var_2004, ig_var_2005) + forfrac_var_1989(jl_var_2021, lay_var_2006) * (forref_var_178(indf_var_2004 + 1, ig_var_2005) - forref_var_178(indf_var_2004, ig_var_2005)))
          taug_var_1981(jl_var_2021, 10 + ig_var_2005, lay_var_2006) = colh2o_var_1994(jl_var_2021, lay_var_2006) * (fac00_var_1985(jl_var_2021, lay_var_2006) * absb_var_176(ind0_var_2001, ig_var_2005) + fac10_var_1987(jl_var_2021, lay_var_2006) * absb_var_176(ind0_var_2001 + 1, ig_var_2005) + fac01_var_1986(jl_var_2021, lay_var_2006) * absb_var_176(ind1_var_2002, ig_var_2005) + fac11_var_1988(jl_var_2021, lay_var_2006) * absb_var_176(ind1_var_2002 + 1, ig_var_2005)) + taufor_var_2007
          fracs_var_2000(jl_var_2021, 10 + ig_var_2005, lay_var_2006) = fracrefb_var_174(ig_var_2005)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol2
SUBROUTINE rrtm_taumol1(kidia_var_2022, kfdia_var_2023, klev_var_2024, taug_var_2026, pavel_var_2025, p_tauaerl_var_2027, fac00_var_2028, fac01_var_2029, fac10_var_2030, fac11_var_2031, forfac_var_2032, forfrac_var_2033, indfor_var_2044, jp_var_2034, jt_var_2035, jt1_var_2036, colh2o_var_2037, laytrop_var_2038, selffac_var_2039, selffrac_var_2040, indself_var_2042, fracs_var_2043, minorfrac_var_2041, indminor_var_2045, scaleminorn2, colbrd_var_2046)
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrta1, ONLY: absa_var_129, absb_var_130, forref_var_133, fracrefa_var_127, fracrefb_var_128, ka_mn2_var_131, kb_mn2, selfref_var_132
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2022
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2023
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2024
  REAL(KIND = 8), INTENT(IN) :: pavel_var_2025(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2026(kidia_var_2022 : kfdia_var_2023, 140, klev_var_2024)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2027(kidia_var_2022 : kfdia_var_2023, klev_var_2024, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2028(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2029(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2030(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2031(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2032(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2033(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2034(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2035(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2036(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2037(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2038(kidia_var_2022 : kfdia_var_2023)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2039(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2040(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2041(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2042(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2043(kidia_var_2022 : kfdia_var_2023, 140, klev_var_2024)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2044(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2045(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(IN) :: scaleminorn2(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  REAL(KIND = 8), INTENT(IN) :: colbrd_var_2046(kidia_var_2022 : kfdia_var_2023, klev_var_2024)
  INTEGER(KIND = 4) :: ind0_var_2047, ind1_var_2048, inds_var_2049
  INTEGER(KIND = 4) :: indf_var_2050, indm_var_2051
  INTEGER(KIND = 4) :: ig_var_2052, lay_var_2053
  REAL(KIND = 8) :: taufor_var_2054, tauself_var_2055, corradj_var_2056, pp_var_2057, scalen2_var_2058, taun2_var_2059
  INTEGER(KIND = 4) :: laytrop_min_var_2060, laytrop_max_var_2061
  INTEGER(KIND = 4) :: ixc_var_2062(klev_var_2024), ixlow_var_2063(kfdia_var_2023, klev_var_2024), ixhigh_var_2064(kfdia_var_2023, klev_var_2024)
  INTEGER(KIND = 4) :: ich_var_2065, icl_var_2066, ixc0_var_2067, ixp_var_2068, jc_var_2069, jl_var_2070
  laytrop_min_var_2060 = MINVAL(laytrop_var_2038)
  laytrop_max_var_2061 = MAXVAL(laytrop_var_2038)
  ixlow_var_2063 = 0
  ixhigh_var_2064 = 0
  ixc_var_2062 = 0
  DO lay_var_2053 = laytrop_min_var_2060 + 1, laytrop_max_var_2061
    icl_var_2066 = 0
    ich_var_2065 = 0
    DO jc_var_2069 = kidia_var_2022, kfdia_var_2023
      IF (lay_var_2053 <= laytrop_var_2038(jc_var_2069)) THEN
        icl_var_2066 = icl_var_2066 + 1
        ixlow_var_2063(icl_var_2066, lay_var_2053) = jc_var_2069
      ELSE
        ich_var_2065 = ich_var_2065 + 1
        ixhigh_var_2064(ich_var_2065, lay_var_2053) = jc_var_2069
      END IF
    END DO
    ixc_var_2062(lay_var_2053) = icl_var_2066
  END DO
  DO lay_var_2053 = 1, laytrop_min_var_2060
    DO jl_var_2070 = kidia_var_2022, kfdia_var_2023
      ind0_var_2047 = ((jp_var_2034(jl_var_2070, lay_var_2053) - 1) * 5 + (jt_var_2035(jl_var_2070, lay_var_2053) - 1)) * nspa_var_236(1) + 1
      ind1_var_2048 = (jp_var_2034(jl_var_2070, lay_var_2053) * 5 + (jt1_var_2036(jl_var_2070, lay_var_2053) - 1)) * nspa_var_236(1) + 1
      inds_var_2049 = indself_var_2042(jl_var_2070, lay_var_2053)
      indf_var_2050 = indfor_var_2044(jl_var_2070, lay_var_2053)
      indm_var_2051 = indminor_var_2045(jl_var_2070, lay_var_2053)
      pp_var_2057 = pavel_var_2025(jl_var_2070, lay_var_2053)
      corradj_var_2056 = 1.0D0
      IF (pp_var_2057 .LT. 250.0D0) THEN
        corradj_var_2056 = 1.0D0 - 0.15D0 * (250.0D0 - pp_var_2057) / 154.4D0
      END IF
      scalen2_var_2058 = colbrd_var_2046(jl_var_2070, lay_var_2053) * scaleminorn2(jl_var_2070, lay_var_2053)
      DO ig_var_2052 = 1, 10
        tauself_var_2055 = selffac_var_2039(jl_var_2070, lay_var_2053) * (selfref_var_132(inds_var_2049, ig_var_2052) + selffrac_var_2040(jl_var_2070, lay_var_2053) * (selfref_var_132(inds_var_2049 + 1, ig_var_2052) - selfref_var_132(inds_var_2049, ig_var_2052)))
        taufor_var_2054 = forfac_var_2032(jl_var_2070, lay_var_2053) * (forref_var_133(indf_var_2050, ig_var_2052) + forfrac_var_2033(jl_var_2070, lay_var_2053) * (forref_var_133(indf_var_2050 + 1, ig_var_2052) - forref_var_133(indf_var_2050, ig_var_2052)))
        taun2_var_2059 = scalen2_var_2058 * (ka_mn2_var_131(indm_var_2051, ig_var_2052) + minorfrac_var_2041(jl_var_2070, lay_var_2053) * (ka_mn2_var_131(indm_var_2051 + 1, ig_var_2052) - ka_mn2_var_131(indm_var_2051, ig_var_2052)))
        taug_var_2026(jl_var_2070, ig_var_2052, lay_var_2053) = corradj_var_2056 * (colh2o_var_2037(jl_var_2070, lay_var_2053) * (fac00_var_2028(jl_var_2070, lay_var_2053) * absa_var_129(ind0_var_2047, ig_var_2052) + fac10_var_2030(jl_var_2070, lay_var_2053) * absa_var_129(ind0_var_2047 + 1, ig_var_2052) + fac01_var_2029(jl_var_2070, lay_var_2053) * absa_var_129(ind1_var_2048, ig_var_2052) + fac11_var_2031(jl_var_2070, lay_var_2053) * absa_var_129(ind1_var_2048 + 1, ig_var_2052)) + tauself_var_2055 + taufor_var_2054 + taun2_var_2059)
        fracs_var_2043(jl_var_2070, ig_var_2052, lay_var_2053) = fracrefa_var_127(ig_var_2052)
      END DO
    END DO
  END DO
  DO lay_var_2053 = laytrop_max_var_2061 + 1, klev_var_2024
    DO jl_var_2070 = kidia_var_2022, kfdia_var_2023
      ind0_var_2047 = ((jp_var_2034(jl_var_2070, lay_var_2053) - 13) * 5 + (jt_var_2035(jl_var_2070, lay_var_2053) - 1)) * nspb_var_237(1) + 1
      ind1_var_2048 = ((jp_var_2034(jl_var_2070, lay_var_2053) - 12) * 5 + (jt1_var_2036(jl_var_2070, lay_var_2053) - 1)) * nspb_var_237(1) + 1
      indf_var_2050 = indfor_var_2044(jl_var_2070, lay_var_2053)
      indm_var_2051 = indminor_var_2045(jl_var_2070, lay_var_2053)
      pp_var_2057 = pavel_var_2025(jl_var_2070, lay_var_2053)
      corradj_var_2056 = 1.0D0 - 0.15D0 * (pp_var_2057 / 95.6D0)
      scalen2_var_2058 = colbrd_var_2046(jl_var_2070, lay_var_2053) * scaleminorn2(jl_var_2070, lay_var_2053)
      DO ig_var_2052 = 1, 10
        taufor_var_2054 = forfac_var_2032(jl_var_2070, lay_var_2053) * (forref_var_133(indf_var_2050, ig_var_2052) + forfrac_var_2033(jl_var_2070, lay_var_2053) * (forref_var_133(indf_var_2050 + 1, ig_var_2052) - forref_var_133(indf_var_2050, ig_var_2052)))
        taun2_var_2059 = scalen2_var_2058 * (kb_mn2(indm_var_2051, ig_var_2052) + minorfrac_var_2041(jl_var_2070, lay_var_2053) * (kb_mn2(indm_var_2051 + 1, ig_var_2052) - kb_mn2(indm_var_2051, ig_var_2052)))
        taug_var_2026(jl_var_2070, ig_var_2052, lay_var_2053) = corradj_var_2056 * (colh2o_var_2037(jl_var_2070, lay_var_2053) * (fac00_var_2028(jl_var_2070, lay_var_2053) * absb_var_130(ind0_var_2047, ig_var_2052) + fac10_var_2030(jl_var_2070, lay_var_2053) * absb_var_130(ind0_var_2047 + 1, ig_var_2052) + fac01_var_2029(jl_var_2070, lay_var_2053) * absb_var_130(ind1_var_2048, ig_var_2052) + fac11_var_2031(jl_var_2070, lay_var_2053) * absb_var_130(ind1_var_2048 + 1, ig_var_2052)) + taufor_var_2054 + taun2_var_2059)
        fracs_var_2043(jl_var_2070, ig_var_2052, lay_var_2053) = fracrefb_var_128(ig_var_2052)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2061 /= laytrop_min_var_2060) THEN
    DO lay_var_2053 = laytrop_min_var_2060 + 1, laytrop_max_var_2061
      ixc0_var_2067 = ixc_var_2062(lay_var_2053)
      DO ixp_var_2068 = 1, ixc0_var_2067
        jl_var_2070 = ixlow_var_2063(ixp_var_2068, lay_var_2053)
        ind0_var_2047 = ((jp_var_2034(jl_var_2070, lay_var_2053) - 1) * 5 + (jt_var_2035(jl_var_2070, lay_var_2053) - 1)) * nspa_var_236(1) + 1
        ind1_var_2048 = (jp_var_2034(jl_var_2070, lay_var_2053) * 5 + (jt1_var_2036(jl_var_2070, lay_var_2053) - 1)) * nspa_var_236(1) + 1
        inds_var_2049 = indself_var_2042(jl_var_2070, lay_var_2053)
        indf_var_2050 = indfor_var_2044(jl_var_2070, lay_var_2053)
        indm_var_2051 = indminor_var_2045(jl_var_2070, lay_var_2053)
        pp_var_2057 = pavel_var_2025(jl_var_2070, lay_var_2053)
        corradj_var_2056 = 1.0D0
        IF (pp_var_2057 .LT. 250.0D0) THEN
          corradj_var_2056 = 1.0D0 - 0.15D0 * (250.0D0 - pp_var_2057) / 154.4D0
        END IF
        scalen2_var_2058 = colbrd_var_2046(jl_var_2070, lay_var_2053) * scaleminorn2(jl_var_2070, lay_var_2053)
        DO ig_var_2052 = 1, 10
          tauself_var_2055 = selffac_var_2039(jl_var_2070, lay_var_2053) * (selfref_var_132(inds_var_2049, ig_var_2052) + selffrac_var_2040(jl_var_2070, lay_var_2053) * (selfref_var_132(inds_var_2049 + 1, ig_var_2052) - selfref_var_132(inds_var_2049, ig_var_2052)))
          taufor_var_2054 = forfac_var_2032(jl_var_2070, lay_var_2053) * (forref_var_133(indf_var_2050, ig_var_2052) + forfrac_var_2033(jl_var_2070, lay_var_2053) * (forref_var_133(indf_var_2050 + 1, ig_var_2052) - forref_var_133(indf_var_2050, ig_var_2052)))
          taun2_var_2059 = scalen2_var_2058 * (ka_mn2_var_131(indm_var_2051, ig_var_2052) + minorfrac_var_2041(jl_var_2070, lay_var_2053) * (ka_mn2_var_131(indm_var_2051 + 1, ig_var_2052) - ka_mn2_var_131(indm_var_2051, ig_var_2052)))
          taug_var_2026(jl_var_2070, ig_var_2052, lay_var_2053) = corradj_var_2056 * (colh2o_var_2037(jl_var_2070, lay_var_2053) * (fac00_var_2028(jl_var_2070, lay_var_2053) * absa_var_129(ind0_var_2047, ig_var_2052) + fac10_var_2030(jl_var_2070, lay_var_2053) * absa_var_129(ind0_var_2047 + 1, ig_var_2052) + fac01_var_2029(jl_var_2070, lay_var_2053) * absa_var_129(ind1_var_2048, ig_var_2052) + fac11_var_2031(jl_var_2070, lay_var_2053) * absa_var_129(ind1_var_2048 + 1, ig_var_2052)) + tauself_var_2055 + taufor_var_2054 + taun2_var_2059)
          fracs_var_2043(jl_var_2070, ig_var_2052, lay_var_2053) = fracrefa_var_127(ig_var_2052)
        END DO
      END DO
      ixc0_var_2067 = kfdia_var_2023 - kidia_var_2022 + 1 - ixc0_var_2067
      DO ixp_var_2068 = 1, ixc0_var_2067
        jl_var_2070 = ixhigh_var_2064(ixp_var_2068, lay_var_2053)
        ind0_var_2047 = ((jp_var_2034(jl_var_2070, lay_var_2053) - 13) * 5 + (jt_var_2035(jl_var_2070, lay_var_2053) - 1)) * nspb_var_237(1) + 1
        ind1_var_2048 = ((jp_var_2034(jl_var_2070, lay_var_2053) - 12) * 5 + (jt1_var_2036(jl_var_2070, lay_var_2053) - 1)) * nspb_var_237(1) + 1
        indf_var_2050 = indfor_var_2044(jl_var_2070, lay_var_2053)
        indm_var_2051 = indminor_var_2045(jl_var_2070, lay_var_2053)
        pp_var_2057 = pavel_var_2025(jl_var_2070, lay_var_2053)
        corradj_var_2056 = 1.0D0 - 0.15D0 * (pp_var_2057 / 95.6D0)
        scalen2_var_2058 = colbrd_var_2046(jl_var_2070, lay_var_2053) * scaleminorn2(jl_var_2070, lay_var_2053)
        DO ig_var_2052 = 1, 10
          taufor_var_2054 = forfac_var_2032(jl_var_2070, lay_var_2053) * (forref_var_133(indf_var_2050, ig_var_2052) + forfrac_var_2033(jl_var_2070, lay_var_2053) * (forref_var_133(indf_var_2050 + 1, ig_var_2052) - forref_var_133(indf_var_2050, ig_var_2052)))
          taun2_var_2059 = scalen2_var_2058 * (kb_mn2(indm_var_2051, ig_var_2052) + minorfrac_var_2041(jl_var_2070, lay_var_2053) * (kb_mn2(indm_var_2051 + 1, ig_var_2052) - kb_mn2(indm_var_2051, ig_var_2052)))
          taug_var_2026(jl_var_2070, ig_var_2052, lay_var_2053) = corradj_var_2056 * (colh2o_var_2037(jl_var_2070, lay_var_2053) * (fac00_var_2028(jl_var_2070, lay_var_2053) * absb_var_130(ind0_var_2047, ig_var_2052) + fac10_var_2030(jl_var_2070, lay_var_2053) * absb_var_130(ind0_var_2047 + 1, ig_var_2052) + fac01_var_2029(jl_var_2070, lay_var_2053) * absb_var_130(ind1_var_2048, ig_var_2052) + fac11_var_2031(jl_var_2070, lay_var_2053) * absb_var_130(ind1_var_2048 + 1, ig_var_2052)) + taufor_var_2054 + taun2_var_2059)
          fracs_var_2043(jl_var_2070, ig_var_2052, lay_var_2053) = fracrefb_var_128(ig_var_2052)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol1
SUBROUTINE rrtm_taumol14(kidia_var_2071, kfdia_var_2072, klev_var_2073, taug_var_2074, p_tauaerl_var_2075, fac00_var_2076, fac01_var_2077, fac10_var_2078, fac11_var_2079, forfac_var_2090, forfrac_var_2091, indfor_var_2089, jp_var_2080, jt_var_2081, jt1_var_2082, colco2_var_2083, laytrop_var_2084, selffac_var_2085, selffrac_var_2086, indself_var_2087, fracs_var_2088)
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrta14, ONLY: absa_var_158, absb_var_159, forref_var_161, fracrefa_var_156, fracrefb_var_157, selfref_var_160
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2071
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2072
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2073
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2074(kidia_var_2071 : kfdia_var_2072, 140, klev_var_2073)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2075(kidia_var_2071 : kfdia_var_2072, klev_var_2073, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2076(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2077(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2078(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2079(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2080(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2081(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2082(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2083(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2084(kidia_var_2071 : kfdia_var_2072)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2085(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2086(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2087(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2088(kidia_var_2071 : kfdia_var_2072, 140, klev_var_2073)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2089(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2090(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2091(kidia_var_2071 : kfdia_var_2072, klev_var_2073)
  INTEGER(KIND = 4) :: ig_var_2092, ind0_var_2093, ind1_var_2094, inds_var_2095, indf_var_2096, lay_var_2097
  REAL(KIND = 8) :: taufor_var_2098, tauself_var_2099
  INTEGER(KIND = 4) :: laytrop_min_var_2100, laytrop_max_var_2101
  INTEGER(KIND = 4) :: ixc_var_2102(klev_var_2073), ixlow_var_2103(kfdia_var_2072, klev_var_2073), ixhigh_var_2104(kfdia_var_2072, klev_var_2073)
  INTEGER(KIND = 4) :: ich_var_2105, icl_var_2106, ixc0_var_2107, ixp_var_2108, jc_var_2109, jl_var_2110
  laytrop_min_var_2100 = MINVAL(laytrop_var_2084)
  laytrop_max_var_2101 = MAXVAL(laytrop_var_2084)
  ixlow_var_2103 = 0
  ixhigh_var_2104 = 0
  ixc_var_2102 = 0
  DO lay_var_2097 = laytrop_min_var_2100 + 1, laytrop_max_var_2101
    icl_var_2106 = 0
    ich_var_2105 = 0
    DO jc_var_2109 = kidia_var_2071, kfdia_var_2072
      IF (lay_var_2097 <= laytrop_var_2084(jc_var_2109)) THEN
        icl_var_2106 = icl_var_2106 + 1
        ixlow_var_2103(icl_var_2106, lay_var_2097) = jc_var_2109
      ELSE
        ich_var_2105 = ich_var_2105 + 1
        ixhigh_var_2104(ich_var_2105, lay_var_2097) = jc_var_2109
      END IF
    END DO
    ixc_var_2102(lay_var_2097) = icl_var_2106
  END DO
  DO lay_var_2097 = 1, laytrop_min_var_2100
    DO jl_var_2110 = kidia_var_2071, kfdia_var_2072
      ind0_var_2093 = ((jp_var_2080(jl_var_2110, lay_var_2097) - 1) * 5 + (jt_var_2081(jl_var_2110, lay_var_2097) - 1)) * nspa_var_236(14) + 1
      ind1_var_2094 = (jp_var_2080(jl_var_2110, lay_var_2097) * 5 + (jt1_var_2082(jl_var_2110, lay_var_2097) - 1)) * nspa_var_236(14) + 1
      inds_var_2095 = indself_var_2087(jl_var_2110, lay_var_2097)
      indf_var_2096 = indfor_var_2089(jl_var_2110, lay_var_2097)
      DO ig_var_2092 = 1, 2
        tauself_var_2099 = selffac_var_2085(jl_var_2110, lay_var_2097) * (selfref_var_160(inds_var_2095, ig_var_2092) + selffrac_var_2086(jl_var_2110, lay_var_2097) * (selfref_var_160(inds_var_2095 + 1, ig_var_2092) - selfref_var_160(inds_var_2095, ig_var_2092)))
        taufor_var_2098 = forfac_var_2090(jl_var_2110, lay_var_2097) * (forref_var_161(indf_var_2096, ig_var_2092) + forfrac_var_2091(jl_var_2110, lay_var_2097) * (forref_var_161(indf_var_2096 + 1, ig_var_2092) - forref_var_161(indf_var_2096, ig_var_2092)))
        taug_var_2074(jl_var_2110, 134 + ig_var_2092, lay_var_2097) = colco2_var_2083(jl_var_2110, lay_var_2097) * (fac00_var_2076(jl_var_2110, lay_var_2097) * absa_var_158(ind0_var_2093, ig_var_2092) + fac10_var_2078(jl_var_2110, lay_var_2097) * absa_var_158(ind0_var_2093 + 1, ig_var_2092) + fac01_var_2077(jl_var_2110, lay_var_2097) * absa_var_158(ind1_var_2094, ig_var_2092) + fac11_var_2079(jl_var_2110, lay_var_2097) * absa_var_158(ind1_var_2094 + 1, ig_var_2092)) + tauself_var_2099 + taufor_var_2098
        fracs_var_2088(jl_var_2110, 134 + ig_var_2092, lay_var_2097) = fracrefa_var_156(ig_var_2092)
      END DO
    END DO
  END DO
  DO lay_var_2097 = laytrop_max_var_2101 + 1, klev_var_2073
    DO jl_var_2110 = kidia_var_2071, kfdia_var_2072
      ind0_var_2093 = ((jp_var_2080(jl_var_2110, lay_var_2097) - 13) * 5 + (jt_var_2081(jl_var_2110, lay_var_2097) - 1)) * nspb_var_237(14) + 1
      ind1_var_2094 = ((jp_var_2080(jl_var_2110, lay_var_2097) - 12) * 5 + (jt1_var_2082(jl_var_2110, lay_var_2097) - 1)) * nspb_var_237(14) + 1
      DO ig_var_2092 = 1, 2
        taug_var_2074(jl_var_2110, 134 + ig_var_2092, lay_var_2097) = colco2_var_2083(jl_var_2110, lay_var_2097) * (fac00_var_2076(jl_var_2110, lay_var_2097) * absb_var_159(ind0_var_2093, ig_var_2092) + fac10_var_2078(jl_var_2110, lay_var_2097) * absb_var_159(ind0_var_2093 + 1, ig_var_2092) + fac01_var_2077(jl_var_2110, lay_var_2097) * absb_var_159(ind1_var_2094, ig_var_2092) + fac11_var_2079(jl_var_2110, lay_var_2097) * absb_var_159(ind1_var_2094 + 1, ig_var_2092))
        fracs_var_2088(jl_var_2110, 134 + ig_var_2092, lay_var_2097) = fracrefb_var_157(ig_var_2092)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2101 /= laytrop_min_var_2100) THEN
    DO lay_var_2097 = laytrop_min_var_2100 + 1, laytrop_max_var_2101
      ixc0_var_2107 = ixc_var_2102(lay_var_2097)
      DO ixp_var_2108 = 1, ixc0_var_2107
        jl_var_2110 = ixlow_var_2103(ixp_var_2108, lay_var_2097)
        ind0_var_2093 = ((jp_var_2080(jl_var_2110, lay_var_2097) - 1) * 5 + (jt_var_2081(jl_var_2110, lay_var_2097) - 1)) * nspa_var_236(14) + 1
        ind1_var_2094 = (jp_var_2080(jl_var_2110, lay_var_2097) * 5 + (jt1_var_2082(jl_var_2110, lay_var_2097) - 1)) * nspa_var_236(14) + 1
        inds_var_2095 = indself_var_2087(jl_var_2110, lay_var_2097)
        indf_var_2096 = indfor_var_2089(jl_var_2110, lay_var_2097)
        DO ig_var_2092 = 1, 2
          tauself_var_2099 = selffac_var_2085(jl_var_2110, lay_var_2097) * (selfref_var_160(inds_var_2095, ig_var_2092) + selffrac_var_2086(jl_var_2110, lay_var_2097) * (selfref_var_160(inds_var_2095 + 1, ig_var_2092) - selfref_var_160(inds_var_2095, ig_var_2092)))
          taufor_var_2098 = forfac_var_2090(jl_var_2110, lay_var_2097) * (forref_var_161(indf_var_2096, ig_var_2092) + forfrac_var_2091(jl_var_2110, lay_var_2097) * (forref_var_161(indf_var_2096 + 1, ig_var_2092) - forref_var_161(indf_var_2096, ig_var_2092)))
          taug_var_2074(jl_var_2110, 134 + ig_var_2092, lay_var_2097) = colco2_var_2083(jl_var_2110, lay_var_2097) * (fac00_var_2076(jl_var_2110, lay_var_2097) * absa_var_158(ind0_var_2093, ig_var_2092) + fac10_var_2078(jl_var_2110, lay_var_2097) * absa_var_158(ind0_var_2093 + 1, ig_var_2092) + fac01_var_2077(jl_var_2110, lay_var_2097) * absa_var_158(ind1_var_2094, ig_var_2092) + fac11_var_2079(jl_var_2110, lay_var_2097) * absa_var_158(ind1_var_2094 + 1, ig_var_2092)) + tauself_var_2099 + taufor_var_2098
          fracs_var_2088(jl_var_2110, 134 + ig_var_2092, lay_var_2097) = fracrefa_var_156(ig_var_2092)
        END DO
      END DO
      ixc0_var_2107 = kfdia_var_2072 - kidia_var_2071 + 1 - ixc0_var_2107
      DO ixp_var_2108 = 1, ixc0_var_2107
        jl_var_2110 = ixhigh_var_2104(ixp_var_2108, lay_var_2097)
        ind0_var_2093 = ((jp_var_2080(jl_var_2110, lay_var_2097) - 13) * 5 + (jt_var_2081(jl_var_2110, lay_var_2097) - 1)) * nspb_var_237(14) + 1
        ind1_var_2094 = ((jp_var_2080(jl_var_2110, lay_var_2097) - 12) * 5 + (jt1_var_2082(jl_var_2110, lay_var_2097) - 1)) * nspb_var_237(14) + 1
        DO ig_var_2092 = 1, 2
          taug_var_2074(jl_var_2110, 134 + ig_var_2092, lay_var_2097) = colco2_var_2083(jl_var_2110, lay_var_2097) * (fac00_var_2076(jl_var_2110, lay_var_2097) * absb_var_159(ind0_var_2093, ig_var_2092) + fac10_var_2078(jl_var_2110, lay_var_2097) * absb_var_159(ind0_var_2093 + 1, ig_var_2092) + fac01_var_2077(jl_var_2110, lay_var_2097) * absb_var_159(ind1_var_2094, ig_var_2092) + fac11_var_2079(jl_var_2110, lay_var_2097) * absb_var_159(ind1_var_2094 + 1, ig_var_2092))
          fracs_var_2088(jl_var_2110, 134 + ig_var_2092, lay_var_2097) = fracrefb_var_157(ig_var_2092)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol14
SUBROUTINE rrtm_taumol15(kidia_var_2111, kfdia_var_2112, klev_var_2113, taug_var_2114, p_tauaerl_var_2115, fac00_var_2116, fac01_var_2117, fac10_var_2118, fac11_var_2119, forfac_var_2133, forfrac_var_2134, indfor_var_2132, jp_var_2120, jt_var_2121, jt1_var_2122, oneminus_var_2123, colh2o_var_2124, colco2_var_2125, coln2o_var_2126, laytrop_var_2127, selffac_var_2128, selffrac_var_2129, indself_var_2130, fracs_var_2131, rat_n2oco2, rat_n2oco2_1, minorfrac_var_2135, indminor_var_2136, scaleminor_var_2137, colbrd_var_2138)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_236
  USE yoerrtm, ONLY: ng15
  USE yoerrta15, ONLY: absa_var_163, forref_var_166, fracrefa_var_162, ka_mn2_var_164, selfref_var_165
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2111
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2112
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2113
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2114(kidia_var_2111 : kfdia_var_2112, 140, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2115(kidia_var_2111 : kfdia_var_2112, klev_var_2113, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2116(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2117(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2118(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2119(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2120(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2121(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2122(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2123
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2124(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2125(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_2126(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2127(kidia_var_2111 : kfdia_var_2112)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2128(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2129(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2130(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2131(kidia_var_2111 : kfdia_var_2112, 140, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: rat_n2oco2(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: rat_n2oco2_1(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2132(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2133(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2134(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2135(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2136(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: scaleminor_var_2137(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  REAL(KIND = 8), INTENT(IN) :: colbrd_var_2138(kidia_var_2111 : kfdia_var_2112, klev_var_2113)
  INTEGER(KIND = 4) :: ig_var_2139, ind0_var_2140, ind1_var_2141, inds_var_2142, indf_var_2143, indm_var_2144, js_var_2145, js1_var_2146, jpl_var_2147, jmn2, lay_var_2148
  REAL(KIND = 8) :: refrat_planck_a_var_2149, refrat_m_a_var_2150
  REAL(KIND = 8) :: taufor_var_2151, tauself_var_2152, tau_major_var_2153(2), tau_major1_var_2154(2), n2m1, n2m2, taun2_var_2155, scalen2_var_2156
  REAL(KIND = 8) :: fac000_var_2157, fac100_var_2158, fac200_var_2159, fac010_var_2160, fac110_var_2161, fac210_var_2162, fac001_var_2163, fac101_var_2164, fac201_var_2165, fac011_var_2166, fac111_var_2167, fac211_var_2168
  REAL(KIND = 8) :: p_var_2169, p4_var_2170, fk0_var_2171, fk1_var_2172, fk2_var_2173
  REAL(KIND = 8) :: fs_var_2174, specmult_var_2175, specparm_var_2176, speccomb_var_2177, fs1_var_2178, specmult1_var_2179, specparm1_var_2180, speccomb1_var_2181, fmn2, specmult_mn2, specparm_mn2, speccomb_mn2, fpl_var_2182, specmult_planck_var_2183, specparm_planck_var_2184, speccomb_planck_var_2185
  INTEGER(KIND = 4) :: laytrop_min_var_2186, laytrop_max_var_2187
  INTEGER(KIND = 4) :: ixc_var_2188(klev_var_2113), ixlow_var_2189(kfdia_var_2112, klev_var_2113), ixhigh_var_2190(kfdia_var_2112, klev_var_2113)
  INTEGER(KIND = 4) :: ich_var_2191, icl_var_2192, ixc0_var_2193, ixp_var_2194, jc_var_2195, jl_var_2196
  laytrop_min_var_2186 = MINVAL(laytrop_var_2127)
  laytrop_max_var_2187 = MAXVAL(laytrop_var_2127)
  ixlow_var_2189 = 0
  ixhigh_var_2190 = 0
  ixc_var_2188 = 0
  DO lay_var_2148 = laytrop_min_var_2186 + 1, laytrop_max_var_2187
    icl_var_2192 = 0
    ich_var_2191 = 0
    DO jc_var_2195 = kidia_var_2111, kfdia_var_2112
      IF (lay_var_2148 <= laytrop_var_2127(jc_var_2195)) THEN
        icl_var_2192 = icl_var_2192 + 1
        ixlow_var_2189(icl_var_2192, lay_var_2148) = jc_var_2195
      ELSE
        ich_var_2191 = ich_var_2191 + 1
        ixhigh_var_2190(ich_var_2191, lay_var_2148) = jc_var_2195
      END IF
    END DO
    ixc_var_2188(lay_var_2148) = icl_var_2192
  END DO
  refrat_planck_a_var_2149 = chi_mls(4, 1) / chi_mls(2, 1)
  refrat_m_a_var_2150 = chi_mls(4, 1) / chi_mls(2, 1)
  DO lay_var_2148 = 1, laytrop_min_var_2186
    DO jl_var_2196 = kidia_var_2111, kfdia_var_2112
      speccomb_var_2177 = coln2o_var_2126(jl_var_2196, lay_var_2148) + rat_n2oco2(jl_var_2196, lay_var_2148) * colco2_var_2125(jl_var_2196, lay_var_2148)
      specparm_var_2176 = MIN(coln2o_var_2126(jl_var_2196, lay_var_2148) / speccomb_var_2177, oneminus_var_2123)
      specmult_var_2175 = 8.0D0 * (specparm_var_2176)
      js_var_2145 = 1 + INT(specmult_var_2175)
      fs_var_2174 = ((specmult_var_2175) - AINT((specmult_var_2175)))
      speccomb1_var_2181 = coln2o_var_2126(jl_var_2196, lay_var_2148) + rat_n2oco2_1(jl_var_2196, lay_var_2148) * colco2_var_2125(jl_var_2196, lay_var_2148)
      specparm1_var_2180 = MIN(coln2o_var_2126(jl_var_2196, lay_var_2148) / speccomb1_var_2181, oneminus_var_2123)
      specmult1_var_2179 = 8.0D0 * (specparm1_var_2180)
      js1_var_2146 = 1 + INT(specmult1_var_2179)
      fs1_var_2178 = ((specmult1_var_2179) - AINT((specmult1_var_2179)))
      speccomb_mn2 = coln2o_var_2126(jl_var_2196, lay_var_2148) + refrat_m_a_var_2150 * colco2_var_2125(jl_var_2196, lay_var_2148)
      specparm_mn2 = MIN(coln2o_var_2126(jl_var_2196, lay_var_2148) / speccomb_mn2, oneminus_var_2123)
      specmult_mn2 = 8.0D0 * specparm_mn2
      jmn2 = 1 + INT(specmult_mn2)
      fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))
      speccomb_planck_var_2185 = coln2o_var_2126(jl_var_2196, lay_var_2148) + refrat_planck_a_var_2149 * colco2_var_2125(jl_var_2196, lay_var_2148)
      specparm_planck_var_2184 = MIN(coln2o_var_2126(jl_var_2196, lay_var_2148) / speccomb_planck_var_2185, oneminus_var_2123)
      specmult_planck_var_2183 = 8.0D0 * specparm_planck_var_2184
      jpl_var_2147 = 1 + INT(specmult_planck_var_2183)
      fpl_var_2182 = ((specmult_planck_var_2183) - AINT((specmult_planck_var_2183)))
      ind0_var_2140 = ((jp_var_2120(jl_var_2196, lay_var_2148) - 1) * 5 + (jt_var_2121(jl_var_2196, lay_var_2148) - 1)) * nspa_var_236(15) + js_var_2145
      ind1_var_2141 = (jp_var_2120(jl_var_2196, lay_var_2148) * 5 + (jt1_var_2122(jl_var_2196, lay_var_2148) - 1)) * nspa_var_236(15) + js1_var_2146
      inds_var_2142 = indself_var_2130(jl_var_2196, lay_var_2148)
      indf_var_2143 = indfor_var_2132(jl_var_2196, lay_var_2148)
      indm_var_2144 = indminor_var_2136(jl_var_2196, lay_var_2148)
      scalen2_var_2156 = colbrd_var_2138(jl_var_2196, lay_var_2148) * scaleminor_var_2137(jl_var_2196, lay_var_2148)
      IF (specparm_var_2176 .LT. 0.125D0) THEN
        p_var_2169 = fs_var_2174 - 1.0D0
        p4_var_2170 = p_var_2169 ** 4
        fk0_var_2171 = p4_var_2170
        fk1_var_2172 = 1.0D0 - p_var_2169 - 2.0D0 * p4_var_2170
        fk2_var_2173 = p_var_2169 + p4_var_2170
        fac000_var_2157 = fk0_var_2171 * fac00_var_2116(jl_var_2196, lay_var_2148)
        fac100_var_2158 = fk1_var_2172 * fac00_var_2116(jl_var_2196, lay_var_2148)
        fac200_var_2159 = fk2_var_2173 * fac00_var_2116(jl_var_2196, lay_var_2148)
        fac010_var_2160 = fk0_var_2171 * fac10_var_2118(jl_var_2196, lay_var_2148)
        fac110_var_2161 = fk1_var_2172 * fac10_var_2118(jl_var_2196, lay_var_2148)
        fac210_var_2162 = fk2_var_2173 * fac10_var_2118(jl_var_2196, lay_var_2148)
      ELSE IF (specparm_var_2176 .GT. 0.875D0) THEN
        p_var_2169 = - fs_var_2174
        p4_var_2170 = p_var_2169 ** 4
        fk0_var_2171 = p4_var_2170
        fk1_var_2172 = 1.0D0 - p_var_2169 - 2.0D0 * p4_var_2170
        fk2_var_2173 = p_var_2169 + p4_var_2170
        fac000_var_2157 = fk0_var_2171 * fac00_var_2116(jl_var_2196, lay_var_2148)
        fac100_var_2158 = fk1_var_2172 * fac00_var_2116(jl_var_2196, lay_var_2148)
        fac200_var_2159 = fk2_var_2173 * fac00_var_2116(jl_var_2196, lay_var_2148)
        fac010_var_2160 = fk0_var_2171 * fac10_var_2118(jl_var_2196, lay_var_2148)
        fac110_var_2161 = fk1_var_2172 * fac10_var_2118(jl_var_2196, lay_var_2148)
        fac210_var_2162 = fk2_var_2173 * fac10_var_2118(jl_var_2196, lay_var_2148)
      ELSE
        fac000_var_2157 = (1.0D0 - fs_var_2174) * fac00_var_2116(jl_var_2196, lay_var_2148)
        fac010_var_2160 = (1.0D0 - fs_var_2174) * fac10_var_2118(jl_var_2196, lay_var_2148)
        fac100_var_2158 = fs_var_2174 * fac00_var_2116(jl_var_2196, lay_var_2148)
        fac110_var_2161 = fs_var_2174 * fac10_var_2118(jl_var_2196, lay_var_2148)
        fac200_var_2159 = 0.0D0
        fac210_var_2162 = 0.0D0
      END IF
      IF (specparm1_var_2180 .LT. 0.125D0) THEN
        p_var_2169 = fs1_var_2178 - 1.0D0
        p4_var_2170 = p_var_2169 ** 4
        fk0_var_2171 = p4_var_2170
        fk1_var_2172 = 1.0D0 - p_var_2169 - 2.0D0 * p4_var_2170
        fk2_var_2173 = p_var_2169 + p4_var_2170
        fac001_var_2163 = fk0_var_2171 * fac01_var_2117(jl_var_2196, lay_var_2148)
        fac101_var_2164 = fk1_var_2172 * fac01_var_2117(jl_var_2196, lay_var_2148)
        fac201_var_2165 = fk2_var_2173 * fac01_var_2117(jl_var_2196, lay_var_2148)
        fac011_var_2166 = fk0_var_2171 * fac11_var_2119(jl_var_2196, lay_var_2148)
        fac111_var_2167 = fk1_var_2172 * fac11_var_2119(jl_var_2196, lay_var_2148)
        fac211_var_2168 = fk2_var_2173 * fac11_var_2119(jl_var_2196, lay_var_2148)
      ELSE IF (specparm1_var_2180 .GT. 0.875D0) THEN
        p_var_2169 = - fs1_var_2178
        p4_var_2170 = p_var_2169 ** 4
        fk0_var_2171 = p4_var_2170
        fk1_var_2172 = 1.0D0 - p_var_2169 - 2.0D0 * p4_var_2170
        fk2_var_2173 = p_var_2169 + p4_var_2170
        fac001_var_2163 = fk0_var_2171 * fac01_var_2117(jl_var_2196, lay_var_2148)
        fac101_var_2164 = fk1_var_2172 * fac01_var_2117(jl_var_2196, lay_var_2148)
        fac201_var_2165 = fk2_var_2173 * fac01_var_2117(jl_var_2196, lay_var_2148)
        fac011_var_2166 = fk0_var_2171 * fac11_var_2119(jl_var_2196, lay_var_2148)
        fac111_var_2167 = fk1_var_2172 * fac11_var_2119(jl_var_2196, lay_var_2148)
        fac211_var_2168 = fk2_var_2173 * fac11_var_2119(jl_var_2196, lay_var_2148)
      ELSE
        fac001_var_2163 = (1.0D0 - fs1_var_2178) * fac01_var_2117(jl_var_2196, lay_var_2148)
        fac011_var_2166 = (1.0D0 - fs1_var_2178) * fac11_var_2119(jl_var_2196, lay_var_2148)
        fac101_var_2164 = fs1_var_2178 * fac01_var_2117(jl_var_2196, lay_var_2148)
        fac111_var_2167 = fs1_var_2178 * fac11_var_2119(jl_var_2196, lay_var_2148)
        fac201_var_2165 = 0.0D0
        fac211_var_2168 = 0.0D0
      END IF
      IF (specparm_var_2176 .LT. 0.125D0) THEN
        tau_major_var_2153(1 : ng15) = speccomb_var_2177 * (fac000_var_2157 * absa_var_163(ind0_var_2140, 1 : 2) + fac100_var_2158 * absa_var_163(ind0_var_2140 + 1, 1 : 2) + fac200_var_2159 * absa_var_163(ind0_var_2140 + 2, 1 : 2) + fac010_var_2160 * absa_var_163(ind0_var_2140 + 9, 1 : 2) + fac110_var_2161 * absa_var_163(ind0_var_2140 + 10, 1 : 2) + fac210_var_2162 * absa_var_163(ind0_var_2140 + 11, 1 : 2))
      ELSE IF (specparm_var_2176 .GT. 0.875D0) THEN
        tau_major_var_2153(1 : ng15) = speccomb_var_2177 * (fac200_var_2159 * absa_var_163(ind0_var_2140 - 1, 1 : 2) + fac100_var_2158 * absa_var_163(ind0_var_2140, 1 : 2) + fac000_var_2157 * absa_var_163(ind0_var_2140 + 1, 1 : 2) + fac210_var_2162 * absa_var_163(ind0_var_2140 + 8, 1 : 2) + fac110_var_2161 * absa_var_163(ind0_var_2140 + 9, 1 : 2) + fac010_var_2160 * absa_var_163(ind0_var_2140 + 10, 1 : 2))
      ELSE
        tau_major_var_2153(1 : ng15) = speccomb_var_2177 * (fac000_var_2157 * absa_var_163(ind0_var_2140, 1 : 2) + fac100_var_2158 * absa_var_163(ind0_var_2140 + 1, 1 : 2) + fac010_var_2160 * absa_var_163(ind0_var_2140 + 9, 1 : 2) + fac110_var_2161 * absa_var_163(ind0_var_2140 + 10, 1 : 2))
      END IF
      IF (specparm1_var_2180 .LT. 0.125D0) THEN
        tau_major1_var_2154(1 : ng15) = speccomb1_var_2181 * (fac001_var_2163 * absa_var_163(ind1_var_2141, 1 : 2) + fac101_var_2164 * absa_var_163(ind1_var_2141 + 1, 1 : 2) + fac201_var_2165 * absa_var_163(ind1_var_2141 + 2, 1 : 2) + fac011_var_2166 * absa_var_163(ind1_var_2141 + 9, 1 : 2) + fac111_var_2167 * absa_var_163(ind1_var_2141 + 10, 1 : 2) + fac211_var_2168 * absa_var_163(ind1_var_2141 + 11, 1 : 2))
      ELSE IF (specparm1_var_2180 .GT. 0.875D0) THEN
        tau_major1_var_2154(1 : ng15) = speccomb1_var_2181 * (fac201_var_2165 * absa_var_163(ind1_var_2141 - 1, 1 : 2) + fac101_var_2164 * absa_var_163(ind1_var_2141, 1 : 2) + fac001_var_2163 * absa_var_163(ind1_var_2141 + 1, 1 : 2) + fac211_var_2168 * absa_var_163(ind1_var_2141 + 8, 1 : 2) + fac111_var_2167 * absa_var_163(ind1_var_2141 + 9, 1 : 2) + fac011_var_2166 * absa_var_163(ind1_var_2141 + 10, 1 : 2))
      ELSE
        tau_major1_var_2154(1 : ng15) = speccomb1_var_2181 * (fac001_var_2163 * absa_var_163(ind1_var_2141, 1 : 2) + fac101_var_2164 * absa_var_163(ind1_var_2141 + 1, 1 : 2) + fac011_var_2166 * absa_var_163(ind1_var_2141 + 9, 1 : 2) + fac111_var_2167 * absa_var_163(ind1_var_2141 + 10, 1 : 2))
      END IF
      DO ig_var_2139 = 1, 2
        tauself_var_2152 = selffac_var_2128(jl_var_2196, lay_var_2148) * (selfref_var_165(inds_var_2142, ig_var_2139) + selffrac_var_2129(jl_var_2196, lay_var_2148) * (selfref_var_165(inds_var_2142 + 1, ig_var_2139) - selfref_var_165(inds_var_2142, ig_var_2139)))
        taufor_var_2151 = forfac_var_2133(jl_var_2196, lay_var_2148) * (forref_var_166(indf_var_2143, ig_var_2139) + forfrac_var_2134(jl_var_2196, lay_var_2148) * (forref_var_166(indf_var_2143 + 1, ig_var_2139) - forref_var_166(indf_var_2143, ig_var_2139)))
        n2m1 = ka_mn2_var_164(jmn2, indm_var_2144, ig_var_2139) + fmn2 * (ka_mn2_var_164(jmn2 + 1, indm_var_2144, ig_var_2139) - ka_mn2_var_164(jmn2, indm_var_2144, ig_var_2139))
        n2m2 = ka_mn2_var_164(jmn2, indm_var_2144 + 1, ig_var_2139) + fmn2 * (ka_mn2_var_164(jmn2 + 1, indm_var_2144 + 1, ig_var_2139) - ka_mn2_var_164(jmn2, indm_var_2144 + 1, ig_var_2139))
        taun2_var_2155 = scalen2_var_2156 * (n2m1 + minorfrac_var_2135(jl_var_2196, lay_var_2148) * (n2m2 - n2m1))
        taug_var_2114(jl_var_2196, 136 + ig_var_2139, lay_var_2148) = tau_major_var_2153(ig_var_2139) + tau_major1_var_2154(ig_var_2139) + tauself_var_2152 + taufor_var_2151 + taun2_var_2155
        fracs_var_2131(jl_var_2196, 136 + ig_var_2139, lay_var_2148) = fracrefa_var_162(ig_var_2139, jpl_var_2147) + fpl_var_2182 * (fracrefa_var_162(ig_var_2139, jpl_var_2147 + 1) - fracrefa_var_162(ig_var_2139, jpl_var_2147))
      END DO
    END DO
  END DO
  DO ig_var_2139 = 1, 2
    DO lay_var_2148 = laytrop_max_var_2187 + 1, klev_var_2113
      DO jl_var_2196 = kidia_var_2111, kfdia_var_2112
        taug_var_2114(jl_var_2196, 136 + ig_var_2139, lay_var_2148) = 0.0D0
        fracs_var_2131(jl_var_2196, 136 + ig_var_2139, lay_var_2148) = 0.0D0
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2187 /= laytrop_min_var_2186) THEN
    DO lay_var_2148 = laytrop_min_var_2186 + 1, laytrop_max_var_2187
      ixc0_var_2193 = ixc_var_2188(lay_var_2148)
      DO ixp_var_2194 = 1, ixc0_var_2193
        jl_var_2196 = ixlow_var_2189(ixp_var_2194, lay_var_2148)
        speccomb_var_2177 = coln2o_var_2126(jl_var_2196, lay_var_2148) + rat_n2oco2(jl_var_2196, lay_var_2148) * colco2_var_2125(jl_var_2196, lay_var_2148)
        specparm_var_2176 = MIN(coln2o_var_2126(jl_var_2196, lay_var_2148) / speccomb_var_2177, oneminus_var_2123)
        specmult_var_2175 = 8.0D0 * (specparm_var_2176)
        js_var_2145 = 1 + INT(specmult_var_2175)
        fs_var_2174 = ((specmult_var_2175) - AINT((specmult_var_2175)))
        speccomb1_var_2181 = coln2o_var_2126(jl_var_2196, lay_var_2148) + rat_n2oco2_1(jl_var_2196, lay_var_2148) * colco2_var_2125(jl_var_2196, lay_var_2148)
        specparm1_var_2180 = MIN(coln2o_var_2126(jl_var_2196, lay_var_2148) / speccomb1_var_2181, oneminus_var_2123)
        specmult1_var_2179 = 8.0D0 * (specparm1_var_2180)
        js1_var_2146 = 1 + INT(specmult1_var_2179)
        fs1_var_2178 = ((specmult1_var_2179) - AINT((specmult1_var_2179)))
        speccomb_mn2 = coln2o_var_2126(jl_var_2196, lay_var_2148) + refrat_m_a_var_2150 * colco2_var_2125(jl_var_2196, lay_var_2148)
        specparm_mn2 = MIN(coln2o_var_2126(jl_var_2196, lay_var_2148) / speccomb_mn2, oneminus_var_2123)
        specmult_mn2 = 8.0D0 * specparm_mn2
        jmn2 = 1 + INT(specmult_mn2)
        fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))
        speccomb_planck_var_2185 = coln2o_var_2126(jl_var_2196, lay_var_2148) + refrat_planck_a_var_2149 * colco2_var_2125(jl_var_2196, lay_var_2148)
        specparm_planck_var_2184 = MIN(coln2o_var_2126(jl_var_2196, lay_var_2148) / speccomb_planck_var_2185, oneminus_var_2123)
        specmult_planck_var_2183 = 8.0D0 * specparm_planck_var_2184
        jpl_var_2147 = 1 + INT(specmult_planck_var_2183)
        fpl_var_2182 = ((specmult_planck_var_2183) - AINT((specmult_planck_var_2183)))
        ind0_var_2140 = ((jp_var_2120(jl_var_2196, lay_var_2148) - 1) * 5 + (jt_var_2121(jl_var_2196, lay_var_2148) - 1)) * nspa_var_236(15) + js_var_2145
        ind1_var_2141 = (jp_var_2120(jl_var_2196, lay_var_2148) * 5 + (jt1_var_2122(jl_var_2196, lay_var_2148) - 1)) * nspa_var_236(15) + js1_var_2146
        inds_var_2142 = indself_var_2130(jl_var_2196, lay_var_2148)
        indf_var_2143 = indfor_var_2132(jl_var_2196, lay_var_2148)
        indm_var_2144 = indminor_var_2136(jl_var_2196, lay_var_2148)
        scalen2_var_2156 = colbrd_var_2138(jl_var_2196, lay_var_2148) * scaleminor_var_2137(jl_var_2196, lay_var_2148)
        IF (specparm_var_2176 .LT. 0.125D0) THEN
          p_var_2169 = fs_var_2174 - 1.0D0
          p4_var_2170 = p_var_2169 ** 4
          fk0_var_2171 = p4_var_2170
          fk1_var_2172 = 1.0D0 - p_var_2169 - 2.0D0 * p4_var_2170
          fk2_var_2173 = p_var_2169 + p4_var_2170
          fac000_var_2157 = fk0_var_2171 * fac00_var_2116(jl_var_2196, lay_var_2148)
          fac100_var_2158 = fk1_var_2172 * fac00_var_2116(jl_var_2196, lay_var_2148)
          fac200_var_2159 = fk2_var_2173 * fac00_var_2116(jl_var_2196, lay_var_2148)
          fac010_var_2160 = fk0_var_2171 * fac10_var_2118(jl_var_2196, lay_var_2148)
          fac110_var_2161 = fk1_var_2172 * fac10_var_2118(jl_var_2196, lay_var_2148)
          fac210_var_2162 = fk2_var_2173 * fac10_var_2118(jl_var_2196, lay_var_2148)
        ELSE IF (specparm_var_2176 .GT. 0.875D0) THEN
          p_var_2169 = - fs_var_2174
          p4_var_2170 = p_var_2169 ** 4
          fk0_var_2171 = p4_var_2170
          fk1_var_2172 = 1.0D0 - p_var_2169 - 2.0D0 * p4_var_2170
          fk2_var_2173 = p_var_2169 + p4_var_2170
          fac000_var_2157 = fk0_var_2171 * fac00_var_2116(jl_var_2196, lay_var_2148)
          fac100_var_2158 = fk1_var_2172 * fac00_var_2116(jl_var_2196, lay_var_2148)
          fac200_var_2159 = fk2_var_2173 * fac00_var_2116(jl_var_2196, lay_var_2148)
          fac010_var_2160 = fk0_var_2171 * fac10_var_2118(jl_var_2196, lay_var_2148)
          fac110_var_2161 = fk1_var_2172 * fac10_var_2118(jl_var_2196, lay_var_2148)
          fac210_var_2162 = fk2_var_2173 * fac10_var_2118(jl_var_2196, lay_var_2148)
        ELSE
          fac000_var_2157 = (1.0D0 - fs_var_2174) * fac00_var_2116(jl_var_2196, lay_var_2148)
          fac010_var_2160 = (1.0D0 - fs_var_2174) * fac10_var_2118(jl_var_2196, lay_var_2148)
          fac100_var_2158 = fs_var_2174 * fac00_var_2116(jl_var_2196, lay_var_2148)
          fac110_var_2161 = fs_var_2174 * fac10_var_2118(jl_var_2196, lay_var_2148)
          fac200_var_2159 = 0.0D0
          fac210_var_2162 = 0.0D0
        END IF
        IF (specparm1_var_2180 .LT. 0.125D0) THEN
          p_var_2169 = fs1_var_2178 - 1.0D0
          p4_var_2170 = p_var_2169 ** 4
          fk0_var_2171 = p4_var_2170
          fk1_var_2172 = 1.0D0 - p_var_2169 - 2.0D0 * p4_var_2170
          fk2_var_2173 = p_var_2169 + p4_var_2170
          fac001_var_2163 = fk0_var_2171 * fac01_var_2117(jl_var_2196, lay_var_2148)
          fac101_var_2164 = fk1_var_2172 * fac01_var_2117(jl_var_2196, lay_var_2148)
          fac201_var_2165 = fk2_var_2173 * fac01_var_2117(jl_var_2196, lay_var_2148)
          fac011_var_2166 = fk0_var_2171 * fac11_var_2119(jl_var_2196, lay_var_2148)
          fac111_var_2167 = fk1_var_2172 * fac11_var_2119(jl_var_2196, lay_var_2148)
          fac211_var_2168 = fk2_var_2173 * fac11_var_2119(jl_var_2196, lay_var_2148)
        ELSE IF (specparm1_var_2180 .GT. 0.875D0) THEN
          p_var_2169 = - fs1_var_2178
          p4_var_2170 = p_var_2169 ** 4
          fk0_var_2171 = p4_var_2170
          fk1_var_2172 = 1.0D0 - p_var_2169 - 2.0D0 * p4_var_2170
          fk2_var_2173 = p_var_2169 + p4_var_2170
          fac001_var_2163 = fk0_var_2171 * fac01_var_2117(jl_var_2196, lay_var_2148)
          fac101_var_2164 = fk1_var_2172 * fac01_var_2117(jl_var_2196, lay_var_2148)
          fac201_var_2165 = fk2_var_2173 * fac01_var_2117(jl_var_2196, lay_var_2148)
          fac011_var_2166 = fk0_var_2171 * fac11_var_2119(jl_var_2196, lay_var_2148)
          fac111_var_2167 = fk1_var_2172 * fac11_var_2119(jl_var_2196, lay_var_2148)
          fac211_var_2168 = fk2_var_2173 * fac11_var_2119(jl_var_2196, lay_var_2148)
        ELSE
          fac001_var_2163 = (1.0D0 - fs1_var_2178) * fac01_var_2117(jl_var_2196, lay_var_2148)
          fac011_var_2166 = (1.0D0 - fs1_var_2178) * fac11_var_2119(jl_var_2196, lay_var_2148)
          fac101_var_2164 = fs1_var_2178 * fac01_var_2117(jl_var_2196, lay_var_2148)
          fac111_var_2167 = fs1_var_2178 * fac11_var_2119(jl_var_2196, lay_var_2148)
          fac201_var_2165 = 0.0D0
          fac211_var_2168 = 0.0D0
        END IF
        IF (specparm_var_2176 .LT. 0.125D0) THEN
          tau_major_var_2153(1 : ng15) = speccomb_var_2177 * (fac000_var_2157 * absa_var_163(ind0_var_2140, 1 : 2) + fac100_var_2158 * absa_var_163(ind0_var_2140 + 1, 1 : 2) + fac200_var_2159 * absa_var_163(ind0_var_2140 + 2, 1 : 2) + fac010_var_2160 * absa_var_163(ind0_var_2140 + 9, 1 : 2) + fac110_var_2161 * absa_var_163(ind0_var_2140 + 10, 1 : 2) + fac210_var_2162 * absa_var_163(ind0_var_2140 + 11, 1 : 2))
        ELSE IF (specparm_var_2176 .GT. 0.875D0) THEN
          tau_major_var_2153(1 : ng15) = speccomb_var_2177 * (fac200_var_2159 * absa_var_163(ind0_var_2140 - 1, 1 : 2) + fac100_var_2158 * absa_var_163(ind0_var_2140, 1 : 2) + fac000_var_2157 * absa_var_163(ind0_var_2140 + 1, 1 : 2) + fac210_var_2162 * absa_var_163(ind0_var_2140 + 8, 1 : 2) + fac110_var_2161 * absa_var_163(ind0_var_2140 + 9, 1 : 2) + fac010_var_2160 * absa_var_163(ind0_var_2140 + 10, 1 : 2))
        ELSE
          tau_major_var_2153(1 : ng15) = speccomb_var_2177 * (fac000_var_2157 * absa_var_163(ind0_var_2140, 1 : 2) + fac100_var_2158 * absa_var_163(ind0_var_2140 + 1, 1 : 2) + fac010_var_2160 * absa_var_163(ind0_var_2140 + 9, 1 : 2) + fac110_var_2161 * absa_var_163(ind0_var_2140 + 10, 1 : 2))
        END IF
        IF (specparm1_var_2180 .LT. 0.125D0) THEN
          tau_major1_var_2154(1 : ng15) = speccomb1_var_2181 * (fac001_var_2163 * absa_var_163(ind1_var_2141, 1 : 2) + fac101_var_2164 * absa_var_163(ind1_var_2141 + 1, 1 : 2) + fac201_var_2165 * absa_var_163(ind1_var_2141 + 2, 1 : 2) + fac011_var_2166 * absa_var_163(ind1_var_2141 + 9, 1 : 2) + fac111_var_2167 * absa_var_163(ind1_var_2141 + 10, 1 : 2) + fac211_var_2168 * absa_var_163(ind1_var_2141 + 11, 1 : 2))
        ELSE IF (specparm1_var_2180 .GT. 0.875D0) THEN
          tau_major1_var_2154(1 : ng15) = speccomb1_var_2181 * (fac201_var_2165 * absa_var_163(ind1_var_2141 - 1, 1 : 2) + fac101_var_2164 * absa_var_163(ind1_var_2141, 1 : 2) + fac001_var_2163 * absa_var_163(ind1_var_2141 + 1, 1 : 2) + fac211_var_2168 * absa_var_163(ind1_var_2141 + 8, 1 : 2) + fac111_var_2167 * absa_var_163(ind1_var_2141 + 9, 1 : 2) + fac011_var_2166 * absa_var_163(ind1_var_2141 + 10, 1 : 2))
        ELSE
          tau_major1_var_2154(1 : ng15) = speccomb1_var_2181 * (fac001_var_2163 * absa_var_163(ind1_var_2141, 1 : 2) + fac101_var_2164 * absa_var_163(ind1_var_2141 + 1, 1 : 2) + fac011_var_2166 * absa_var_163(ind1_var_2141 + 9, 1 : 2) + fac111_var_2167 * absa_var_163(ind1_var_2141 + 10, 1 : 2))
        END IF
        DO ig_var_2139 = 1, 2
          tauself_var_2152 = selffac_var_2128(jl_var_2196, lay_var_2148) * (selfref_var_165(inds_var_2142, ig_var_2139) + selffrac_var_2129(jl_var_2196, lay_var_2148) * (selfref_var_165(inds_var_2142 + 1, ig_var_2139) - selfref_var_165(inds_var_2142, ig_var_2139)))
          taufor_var_2151 = forfac_var_2133(jl_var_2196, lay_var_2148) * (forref_var_166(indf_var_2143, ig_var_2139) + forfrac_var_2134(jl_var_2196, lay_var_2148) * (forref_var_166(indf_var_2143 + 1, ig_var_2139) - forref_var_166(indf_var_2143, ig_var_2139)))
          n2m1 = ka_mn2_var_164(jmn2, indm_var_2144, ig_var_2139) + fmn2 * (ka_mn2_var_164(jmn2 + 1, indm_var_2144, ig_var_2139) - ka_mn2_var_164(jmn2, indm_var_2144, ig_var_2139))
          n2m2 = ka_mn2_var_164(jmn2, indm_var_2144 + 1, ig_var_2139) + fmn2 * (ka_mn2_var_164(jmn2 + 1, indm_var_2144 + 1, ig_var_2139) - ka_mn2_var_164(jmn2, indm_var_2144 + 1, ig_var_2139))
          taun2_var_2155 = scalen2_var_2156 * (n2m1 + minorfrac_var_2135(jl_var_2196, lay_var_2148) * (n2m2 - n2m1))
          taug_var_2114(jl_var_2196, 136 + ig_var_2139, lay_var_2148) = tau_major_var_2153(ig_var_2139) + tau_major1_var_2154(ig_var_2139) + tauself_var_2152 + taufor_var_2151 + taun2_var_2155
          fracs_var_2131(jl_var_2196, 136 + ig_var_2139, lay_var_2148) = fracrefa_var_162(ig_var_2139, jpl_var_2147) + fpl_var_2182 * (fracrefa_var_162(ig_var_2139, jpl_var_2147 + 1) - fracrefa_var_162(ig_var_2139, jpl_var_2147))
        END DO
      END DO
      ixc0_var_2193 = kfdia_var_2112 - kidia_var_2111 + 1 - ixc0_var_2193
      DO ig_var_2139 = 1, 2
        DO ixp_var_2194 = 1, ixc0_var_2193
          jl_var_2196 = ixhigh_var_2190(ixp_var_2194, lay_var_2148)
          taug_var_2114(jl_var_2196, 136 + ig_var_2139, lay_var_2148) = 0.0D0
          fracs_var_2131(jl_var_2196, 136 + ig_var_2139, lay_var_2148) = 0.0D0
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol15
SUBROUTINE rrtm_taumol16(kidia_var_2197, kfdia_var_2198, klev_var_2199, taug_var_2200, p_tauaerl_var_2201, fac00_var_2202, fac01_var_2203, fac10_var_2204, fac11_var_2205, forfac_var_2220, forfrac_var_2221, indfor_var_2219, jp_var_2206, jt_var_2207, jt1_var_2208, oneminus_var_2209, colh2o_var_2210, colch4_var_2211, laytrop_var_2212, selffac_var_2213, selffrac_var_2214, indself_var_2215, fracs_var_2216, rat_h2och4_var_2217, rat_h2och4_1_var_2218)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_236, nspb_var_237
  USE yoerrtm, ONLY: ng16
  USE yoerrta16, ONLY: absa_var_169, absb_var_170, forref_var_172, fracrefa_var_167, fracrefb_var_168, selfref_var_171
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2197
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2198
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2199
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2200(kidia_var_2197 : kfdia_var_2198, 140, klev_var_2199)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2201(kidia_var_2197 : kfdia_var_2198, klev_var_2199, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2202(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2203(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2204(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2205(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2206(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2207(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2208(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2209
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2210(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  REAL(KIND = 8), INTENT(IN) :: colch4_var_2211(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2212(kidia_var_2197 : kfdia_var_2198)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2213(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2214(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2215(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2216(kidia_var_2197 : kfdia_var_2198, 140, klev_var_2199)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_var_2217(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_1_var_2218(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2219(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2220(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2221(kidia_var_2197 : kfdia_var_2198, klev_var_2199)
  INTEGER(KIND = 4) :: ig_var_2222, ind0_var_2223, ind1_var_2224, inds_var_2225, indf_var_2226, js_var_2227, js1_var_2228, jpl_var_2229, lay_var_2230
  REAL(KIND = 8) :: fac000_var_2231, fac100_var_2232, fac200_var_2233, fac010_var_2234, fac110_var_2235, fac210_var_2236, fac001_var_2237, fac101_var_2238, fac201_var_2239, fac011_var_2240, fac111_var_2241, fac211_var_2242
  REAL(KIND = 8) :: p_var_2243, p4_var_2244, fk0_var_2245, fk1_var_2246, fk2_var_2247
  REAL(KIND = 8) :: refrat_planck_a_var_2248
  REAL(KIND = 8) :: taufor_var_2249, tauself_var_2250, tau_major_var_2251(2), tau_major1_var_2252(2)
  REAL(KIND = 8) :: fs_var_2253, specmult_var_2254, specparm_var_2255, speccomb_var_2256, fs1_var_2257, specmult1_var_2258, specparm1_var_2259, speccomb1_var_2260, fpl_var_2261, specmult_planck_var_2262, specparm_planck_var_2263, speccomb_planck_var_2264
  INTEGER(KIND = 4) :: laytrop_min_var_2265, laytrop_max_var_2266
  INTEGER(KIND = 4) :: ixc_var_2267(klev_var_2199), ixlow_var_2268(kfdia_var_2198, klev_var_2199), ixhigh_var_2269(kfdia_var_2198, klev_var_2199)
  INTEGER(KIND = 4) :: ich_var_2270, icl_var_2271, ixc0_var_2272, ixp_var_2273, jc_var_2274, jl_var_2275
  laytrop_min_var_2265 = MINVAL(laytrop_var_2212)
  laytrop_max_var_2266 = MAXVAL(laytrop_var_2212)
  ixlow_var_2268 = 0
  ixhigh_var_2269 = 0
  ixc_var_2267 = 0
  DO lay_var_2230 = laytrop_min_var_2265 + 1, laytrop_max_var_2266
    icl_var_2271 = 0
    ich_var_2270 = 0
    DO jc_var_2274 = kidia_var_2197, kfdia_var_2198
      IF (lay_var_2230 <= laytrop_var_2212(jc_var_2274)) THEN
        icl_var_2271 = icl_var_2271 + 1
        ixlow_var_2268(icl_var_2271, lay_var_2230) = jc_var_2274
      ELSE
        ich_var_2270 = ich_var_2270 + 1
        ixhigh_var_2269(ich_var_2270, lay_var_2230) = jc_var_2274
      END IF
    END DO
    ixc_var_2267(lay_var_2230) = icl_var_2271
  END DO
  refrat_planck_a_var_2248 = chi_mls(1, 6) / chi_mls(6, 6)
  DO lay_var_2230 = 1, laytrop_min_var_2265
    DO jl_var_2275 = kidia_var_2197, kfdia_var_2198
      speccomb_var_2256 = colh2o_var_2210(jl_var_2275, lay_var_2230) + rat_h2och4_var_2217(jl_var_2275, lay_var_2230) * colch4_var_2211(jl_var_2275, lay_var_2230)
      specparm_var_2255 = MIN(colh2o_var_2210(jl_var_2275, lay_var_2230) / speccomb_var_2256, oneminus_var_2209)
      specmult_var_2254 = 8.0D0 * (specparm_var_2255)
      js_var_2227 = 1 + INT(specmult_var_2254)
      fs_var_2253 = ((specmult_var_2254) - AINT((specmult_var_2254)))
      speccomb1_var_2260 = colh2o_var_2210(jl_var_2275, lay_var_2230) + rat_h2och4_1_var_2218(jl_var_2275, lay_var_2230) * colch4_var_2211(jl_var_2275, lay_var_2230)
      specparm1_var_2259 = MIN(colh2o_var_2210(jl_var_2275, lay_var_2230) / speccomb1_var_2260, oneminus_var_2209)
      specmult1_var_2258 = 8.0D0 * (specparm1_var_2259)
      js1_var_2228 = 1 + INT(specmult1_var_2258)
      fs1_var_2257 = ((specmult1_var_2258) - AINT((specmult1_var_2258)))
      speccomb_planck_var_2264 = colh2o_var_2210(jl_var_2275, lay_var_2230) + refrat_planck_a_var_2248 * colch4_var_2211(jl_var_2275, lay_var_2230)
      specparm_planck_var_2263 = MIN(colh2o_var_2210(jl_var_2275, lay_var_2230) / speccomb_planck_var_2264, oneminus_var_2209)
      specmult_planck_var_2262 = 8.0D0 * specparm_planck_var_2263
      jpl_var_2229 = 1 + INT(specmult_planck_var_2262)
      fpl_var_2261 = ((specmult_planck_var_2262) - AINT((specmult_planck_var_2262)))
      ind0_var_2223 = ((jp_var_2206(jl_var_2275, lay_var_2230) - 1) * 5 + (jt_var_2207(jl_var_2275, lay_var_2230) - 1)) * nspa_var_236(16) + js_var_2227
      ind1_var_2224 = (jp_var_2206(jl_var_2275, lay_var_2230) * 5 + (jt1_var_2208(jl_var_2275, lay_var_2230) - 1)) * nspa_var_236(16) + js1_var_2228
      inds_var_2225 = indself_var_2215(jl_var_2275, lay_var_2230)
      indf_var_2226 = indfor_var_2219(jl_var_2275, lay_var_2230)
      IF (specparm_var_2255 .LT. 0.125D0) THEN
        p_var_2243 = fs_var_2253 - 1.0D0
        p4_var_2244 = p_var_2243 ** 4
        fk0_var_2245 = p4_var_2244
        fk1_var_2246 = 1.0D0 - p_var_2243 - 2.0D0 * p4_var_2244
        fk2_var_2247 = p_var_2243 + p4_var_2244
        fac000_var_2231 = fk0_var_2245 * fac00_var_2202(jl_var_2275, lay_var_2230)
        fac100_var_2232 = fk1_var_2246 * fac00_var_2202(jl_var_2275, lay_var_2230)
        fac200_var_2233 = fk2_var_2247 * fac00_var_2202(jl_var_2275, lay_var_2230)
        fac010_var_2234 = fk0_var_2245 * fac10_var_2204(jl_var_2275, lay_var_2230)
        fac110_var_2235 = fk1_var_2246 * fac10_var_2204(jl_var_2275, lay_var_2230)
        fac210_var_2236 = fk2_var_2247 * fac10_var_2204(jl_var_2275, lay_var_2230)
      ELSE IF (specparm_var_2255 .GT. 0.875D0) THEN
        p_var_2243 = - fs_var_2253
        p4_var_2244 = p_var_2243 ** 4
        fk0_var_2245 = p4_var_2244
        fk1_var_2246 = 1.0D0 - p_var_2243 - 2.0D0 * p4_var_2244
        fk2_var_2247 = p_var_2243 + p4_var_2244
        fac000_var_2231 = fk0_var_2245 * fac00_var_2202(jl_var_2275, lay_var_2230)
        fac100_var_2232 = fk1_var_2246 * fac00_var_2202(jl_var_2275, lay_var_2230)
        fac200_var_2233 = fk2_var_2247 * fac00_var_2202(jl_var_2275, lay_var_2230)
        fac010_var_2234 = fk0_var_2245 * fac10_var_2204(jl_var_2275, lay_var_2230)
        fac110_var_2235 = fk1_var_2246 * fac10_var_2204(jl_var_2275, lay_var_2230)
        fac210_var_2236 = fk2_var_2247 * fac10_var_2204(jl_var_2275, lay_var_2230)
      ELSE
        fac000_var_2231 = (1.0D0 - fs_var_2253) * fac00_var_2202(jl_var_2275, lay_var_2230)
        fac010_var_2234 = (1.0D0 - fs_var_2253) * fac10_var_2204(jl_var_2275, lay_var_2230)
        fac100_var_2232 = fs_var_2253 * fac00_var_2202(jl_var_2275, lay_var_2230)
        fac110_var_2235 = fs_var_2253 * fac10_var_2204(jl_var_2275, lay_var_2230)
        fac200_var_2233 = 0.0D0
        fac210_var_2236 = 0.0D0
      END IF
      IF (specparm1_var_2259 .LT. 0.125D0) THEN
        p_var_2243 = fs1_var_2257 - 1.0D0
        p4_var_2244 = p_var_2243 ** 4
        fk0_var_2245 = p4_var_2244
        fk1_var_2246 = 1.0D0 - p_var_2243 - 2.0D0 * p4_var_2244
        fk2_var_2247 = p_var_2243 + p4_var_2244
        fac001_var_2237 = fk0_var_2245 * fac01_var_2203(jl_var_2275, lay_var_2230)
        fac101_var_2238 = fk1_var_2246 * fac01_var_2203(jl_var_2275, lay_var_2230)
        fac201_var_2239 = fk2_var_2247 * fac01_var_2203(jl_var_2275, lay_var_2230)
        fac011_var_2240 = fk0_var_2245 * fac11_var_2205(jl_var_2275, lay_var_2230)
        fac111_var_2241 = fk1_var_2246 * fac11_var_2205(jl_var_2275, lay_var_2230)
        fac211_var_2242 = fk2_var_2247 * fac11_var_2205(jl_var_2275, lay_var_2230)
      ELSE IF (specparm1_var_2259 .GT. 0.875D0) THEN
        p_var_2243 = - fs1_var_2257
        p4_var_2244 = p_var_2243 ** 4
        fk0_var_2245 = p4_var_2244
        fk1_var_2246 = 1.0D0 - p_var_2243 - 2.0D0 * p4_var_2244
        fk2_var_2247 = p_var_2243 + p4_var_2244
        fac001_var_2237 = fk0_var_2245 * fac01_var_2203(jl_var_2275, lay_var_2230)
        fac101_var_2238 = fk1_var_2246 * fac01_var_2203(jl_var_2275, lay_var_2230)
        fac201_var_2239 = fk2_var_2247 * fac01_var_2203(jl_var_2275, lay_var_2230)
        fac011_var_2240 = fk0_var_2245 * fac11_var_2205(jl_var_2275, lay_var_2230)
        fac111_var_2241 = fk1_var_2246 * fac11_var_2205(jl_var_2275, lay_var_2230)
        fac211_var_2242 = fk2_var_2247 * fac11_var_2205(jl_var_2275, lay_var_2230)
      ELSE
        fac001_var_2237 = (1.0D0 - fs1_var_2257) * fac01_var_2203(jl_var_2275, lay_var_2230)
        fac011_var_2240 = (1.0D0 - fs1_var_2257) * fac11_var_2205(jl_var_2275, lay_var_2230)
        fac101_var_2238 = fs1_var_2257 * fac01_var_2203(jl_var_2275, lay_var_2230)
        fac111_var_2241 = fs1_var_2257 * fac11_var_2205(jl_var_2275, lay_var_2230)
        fac201_var_2239 = 0.0D0
        fac211_var_2242 = 0.0D0
      END IF
      IF (specparm_var_2255 .LT. 0.125D0) THEN
        tau_major_var_2251(1 : ng16) = speccomb_var_2256 * (fac000_var_2231 * absa_var_169(ind0_var_2223, 1 : 2) + fac100_var_2232 * absa_var_169(ind0_var_2223 + 1, 1 : 2) + fac200_var_2233 * absa_var_169(ind0_var_2223 + 2, 1 : 2) + fac010_var_2234 * absa_var_169(ind0_var_2223 + 9, 1 : 2) + fac110_var_2235 * absa_var_169(ind0_var_2223 + 10, 1 : 2) + fac210_var_2236 * absa_var_169(ind0_var_2223 + 11, 1 : 2))
      ELSE IF (specparm_var_2255 .GT. 0.875D0) THEN
        tau_major_var_2251(1 : ng16) = speccomb_var_2256 * (fac200_var_2233 * absa_var_169(ind0_var_2223 - 1, 1 : 2) + fac100_var_2232 * absa_var_169(ind0_var_2223, 1 : 2) + fac000_var_2231 * absa_var_169(ind0_var_2223 + 1, 1 : 2) + fac210_var_2236 * absa_var_169(ind0_var_2223 + 8, 1 : 2) + fac110_var_2235 * absa_var_169(ind0_var_2223 + 9, 1 : 2) + fac010_var_2234 * absa_var_169(ind0_var_2223 + 10, 1 : 2))
      ELSE
        tau_major_var_2251(1 : ng16) = speccomb_var_2256 * (fac000_var_2231 * absa_var_169(ind0_var_2223, 1 : 2) + fac100_var_2232 * absa_var_169(ind0_var_2223 + 1, 1 : 2) + fac010_var_2234 * absa_var_169(ind0_var_2223 + 9, 1 : 2) + fac110_var_2235 * absa_var_169(ind0_var_2223 + 10, 1 : 2))
      END IF
      IF (specparm1_var_2259 .LT. 0.125D0) THEN
        tau_major1_var_2252(1 : ng16) = speccomb1_var_2260 * (fac001_var_2237 * absa_var_169(ind1_var_2224, 1 : 2) + fac101_var_2238 * absa_var_169(ind1_var_2224 + 1, 1 : 2) + fac201_var_2239 * absa_var_169(ind1_var_2224 + 2, 1 : 2) + fac011_var_2240 * absa_var_169(ind1_var_2224 + 9, 1 : 2) + fac111_var_2241 * absa_var_169(ind1_var_2224 + 10, 1 : 2) + fac211_var_2242 * absa_var_169(ind1_var_2224 + 11, 1 : 2))
      ELSE IF (specparm1_var_2259 .GT. 0.875D0) THEN
        tau_major1_var_2252(1 : ng16) = speccomb1_var_2260 * (fac201_var_2239 * absa_var_169(ind1_var_2224 - 1, 1 : 2) + fac101_var_2238 * absa_var_169(ind1_var_2224, 1 : 2) + fac001_var_2237 * absa_var_169(ind1_var_2224 + 1, 1 : 2) + fac211_var_2242 * absa_var_169(ind1_var_2224 + 8, 1 : 2) + fac111_var_2241 * absa_var_169(ind1_var_2224 + 9, 1 : 2) + fac011_var_2240 * absa_var_169(ind1_var_2224 + 10, 1 : 2))
      ELSE
        tau_major1_var_2252(1 : ng16) = speccomb1_var_2260 * (fac001_var_2237 * absa_var_169(ind1_var_2224, 1 : 2) + fac101_var_2238 * absa_var_169(ind1_var_2224 + 1, 1 : 2) + fac011_var_2240 * absa_var_169(ind1_var_2224 + 9, 1 : 2) + fac111_var_2241 * absa_var_169(ind1_var_2224 + 10, 1 : 2))
      END IF
      DO ig_var_2222 = 1, 2
        tauself_var_2250 = selffac_var_2213(jl_var_2275, lay_var_2230) * (selfref_var_171(inds_var_2225, ig_var_2222) + selffrac_var_2214(jl_var_2275, lay_var_2230) * (selfref_var_171(inds_var_2225 + 1, ig_var_2222) - selfref_var_171(inds_var_2225, ig_var_2222)))
        taufor_var_2249 = forfac_var_2220(jl_var_2275, lay_var_2230) * (forref_var_172(indf_var_2226, ig_var_2222) + forfrac_var_2221(jl_var_2275, lay_var_2230) * (forref_var_172(indf_var_2226 + 1, ig_var_2222) - forref_var_172(indf_var_2226, ig_var_2222)))
        taug_var_2200(jl_var_2275, 138 + ig_var_2222, lay_var_2230) = tau_major_var_2251(ig_var_2222) + tau_major1_var_2252(ig_var_2222) + tauself_var_2250 + taufor_var_2249
        fracs_var_2216(jl_var_2275, 138 + ig_var_2222, lay_var_2230) = fracrefa_var_167(ig_var_2222, jpl_var_2229) + fpl_var_2261 * (fracrefa_var_167(ig_var_2222, jpl_var_2229 + 1) - fracrefa_var_167(ig_var_2222, jpl_var_2229))
      END DO
    END DO
  END DO
  DO lay_var_2230 = laytrop_max_var_2266 + 1, klev_var_2199
    DO jl_var_2275 = kidia_var_2197, kfdia_var_2198
      ind0_var_2223 = ((jp_var_2206(jl_var_2275, lay_var_2230) - 13) * 5 + (jt_var_2207(jl_var_2275, lay_var_2230) - 1)) * nspb_var_237(16) + 1
      ind1_var_2224 = ((jp_var_2206(jl_var_2275, lay_var_2230) - 12) * 5 + (jt1_var_2208(jl_var_2275, lay_var_2230) - 1)) * nspb_var_237(16) + 1
      DO ig_var_2222 = 1, 2
        taug_var_2200(jl_var_2275, 138 + ig_var_2222, lay_var_2230) = colch4_var_2211(jl_var_2275, lay_var_2230) * (fac00_var_2202(jl_var_2275, lay_var_2230) * absb_var_170(ind0_var_2223, ig_var_2222) + fac10_var_2204(jl_var_2275, lay_var_2230) * absb_var_170(ind0_var_2223 + 1, ig_var_2222) + fac01_var_2203(jl_var_2275, lay_var_2230) * absb_var_170(ind1_var_2224, ig_var_2222) + fac11_var_2205(jl_var_2275, lay_var_2230) * absb_var_170(ind1_var_2224 + 1, ig_var_2222))
        fracs_var_2216(jl_var_2275, 138 + ig_var_2222, lay_var_2230) = fracrefb_var_168(ig_var_2222)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2266 /= laytrop_min_var_2265) THEN
    DO lay_var_2230 = laytrop_min_var_2265 + 1, laytrop_max_var_2266
      ixc0_var_2272 = ixc_var_2267(lay_var_2230)
      DO ixp_var_2273 = 1, ixc0_var_2272
        jl_var_2275 = ixlow_var_2268(ixp_var_2273, lay_var_2230)
        speccomb_var_2256 = colh2o_var_2210(jl_var_2275, lay_var_2230) + rat_h2och4_var_2217(jl_var_2275, lay_var_2230) * colch4_var_2211(jl_var_2275, lay_var_2230)
        specparm_var_2255 = MIN(colh2o_var_2210(jl_var_2275, lay_var_2230) / speccomb_var_2256, oneminus_var_2209)
        specmult_var_2254 = 8.0D0 * (specparm_var_2255)
        js_var_2227 = 1 + INT(specmult_var_2254)
        fs_var_2253 = ((specmult_var_2254) - AINT((specmult_var_2254)))
        speccomb1_var_2260 = colh2o_var_2210(jl_var_2275, lay_var_2230) + rat_h2och4_1_var_2218(jl_var_2275, lay_var_2230) * colch4_var_2211(jl_var_2275, lay_var_2230)
        specparm1_var_2259 = MIN(colh2o_var_2210(jl_var_2275, lay_var_2230) / speccomb1_var_2260, oneminus_var_2209)
        specmult1_var_2258 = 8.0D0 * (specparm1_var_2259)
        js1_var_2228 = 1 + INT(specmult1_var_2258)
        fs1_var_2257 = ((specmult1_var_2258) - AINT((specmult1_var_2258)))
        speccomb_planck_var_2264 = colh2o_var_2210(jl_var_2275, lay_var_2230) + refrat_planck_a_var_2248 * colch4_var_2211(jl_var_2275, lay_var_2230)
        specparm_planck_var_2263 = MIN(colh2o_var_2210(jl_var_2275, lay_var_2230) / speccomb_planck_var_2264, oneminus_var_2209)
        specmult_planck_var_2262 = 8.0D0 * specparm_planck_var_2263
        jpl_var_2229 = 1 + INT(specmult_planck_var_2262)
        fpl_var_2261 = ((specmult_planck_var_2262) - AINT((specmult_planck_var_2262)))
        ind0_var_2223 = ((jp_var_2206(jl_var_2275, lay_var_2230) - 1) * 5 + (jt_var_2207(jl_var_2275, lay_var_2230) - 1)) * nspa_var_236(16) + js_var_2227
        ind1_var_2224 = (jp_var_2206(jl_var_2275, lay_var_2230) * 5 + (jt1_var_2208(jl_var_2275, lay_var_2230) - 1)) * nspa_var_236(16) + js1_var_2228
        inds_var_2225 = indself_var_2215(jl_var_2275, lay_var_2230)
        indf_var_2226 = indfor_var_2219(jl_var_2275, lay_var_2230)
        IF (specparm_var_2255 .LT. 0.125D0) THEN
          p_var_2243 = fs_var_2253 - 1.0D0
          p4_var_2244 = p_var_2243 ** 4
          fk0_var_2245 = p4_var_2244
          fk1_var_2246 = 1.0D0 - p_var_2243 - 2.0D0 * p4_var_2244
          fk2_var_2247 = p_var_2243 + p4_var_2244
          fac000_var_2231 = fk0_var_2245 * fac00_var_2202(jl_var_2275, lay_var_2230)
          fac100_var_2232 = fk1_var_2246 * fac00_var_2202(jl_var_2275, lay_var_2230)
          fac200_var_2233 = fk2_var_2247 * fac00_var_2202(jl_var_2275, lay_var_2230)
          fac010_var_2234 = fk0_var_2245 * fac10_var_2204(jl_var_2275, lay_var_2230)
          fac110_var_2235 = fk1_var_2246 * fac10_var_2204(jl_var_2275, lay_var_2230)
          fac210_var_2236 = fk2_var_2247 * fac10_var_2204(jl_var_2275, lay_var_2230)
        ELSE IF (specparm_var_2255 .GT. 0.875D0) THEN
          p_var_2243 = - fs_var_2253
          p4_var_2244 = p_var_2243 ** 4
          fk0_var_2245 = p4_var_2244
          fk1_var_2246 = 1.0D0 - p_var_2243 - 2.0D0 * p4_var_2244
          fk2_var_2247 = p_var_2243 + p4_var_2244
          fac000_var_2231 = fk0_var_2245 * fac00_var_2202(jl_var_2275, lay_var_2230)
          fac100_var_2232 = fk1_var_2246 * fac00_var_2202(jl_var_2275, lay_var_2230)
          fac200_var_2233 = fk2_var_2247 * fac00_var_2202(jl_var_2275, lay_var_2230)
          fac010_var_2234 = fk0_var_2245 * fac10_var_2204(jl_var_2275, lay_var_2230)
          fac110_var_2235 = fk1_var_2246 * fac10_var_2204(jl_var_2275, lay_var_2230)
          fac210_var_2236 = fk2_var_2247 * fac10_var_2204(jl_var_2275, lay_var_2230)
        ELSE
          fac000_var_2231 = (1.0D0 - fs_var_2253) * fac00_var_2202(jl_var_2275, lay_var_2230)
          fac010_var_2234 = (1.0D0 - fs_var_2253) * fac10_var_2204(jl_var_2275, lay_var_2230)
          fac100_var_2232 = fs_var_2253 * fac00_var_2202(jl_var_2275, lay_var_2230)
          fac110_var_2235 = fs_var_2253 * fac10_var_2204(jl_var_2275, lay_var_2230)
          fac200_var_2233 = 0.0D0
          fac210_var_2236 = 0.0D0
        END IF
        IF (specparm1_var_2259 .LT. 0.125D0) THEN
          p_var_2243 = fs1_var_2257 - 1.0D0
          p4_var_2244 = p_var_2243 ** 4
          fk0_var_2245 = p4_var_2244
          fk1_var_2246 = 1.0D0 - p_var_2243 - 2.0D0 * p4_var_2244
          fk2_var_2247 = p_var_2243 + p4_var_2244
          fac001_var_2237 = fk0_var_2245 * fac01_var_2203(jl_var_2275, lay_var_2230)
          fac101_var_2238 = fk1_var_2246 * fac01_var_2203(jl_var_2275, lay_var_2230)
          fac201_var_2239 = fk2_var_2247 * fac01_var_2203(jl_var_2275, lay_var_2230)
          fac011_var_2240 = fk0_var_2245 * fac11_var_2205(jl_var_2275, lay_var_2230)
          fac111_var_2241 = fk1_var_2246 * fac11_var_2205(jl_var_2275, lay_var_2230)
          fac211_var_2242 = fk2_var_2247 * fac11_var_2205(jl_var_2275, lay_var_2230)
        ELSE IF (specparm1_var_2259 .GT. 0.875D0) THEN
          p_var_2243 = - fs1_var_2257
          p4_var_2244 = p_var_2243 ** 4
          fk0_var_2245 = p4_var_2244
          fk1_var_2246 = 1.0D0 - p_var_2243 - 2.0D0 * p4_var_2244
          fk2_var_2247 = p_var_2243 + p4_var_2244
          fac001_var_2237 = fk0_var_2245 * fac01_var_2203(jl_var_2275, lay_var_2230)
          fac101_var_2238 = fk1_var_2246 * fac01_var_2203(jl_var_2275, lay_var_2230)
          fac201_var_2239 = fk2_var_2247 * fac01_var_2203(jl_var_2275, lay_var_2230)
          fac011_var_2240 = fk0_var_2245 * fac11_var_2205(jl_var_2275, lay_var_2230)
          fac111_var_2241 = fk1_var_2246 * fac11_var_2205(jl_var_2275, lay_var_2230)
          fac211_var_2242 = fk2_var_2247 * fac11_var_2205(jl_var_2275, lay_var_2230)
        ELSE
          fac001_var_2237 = (1.0D0 - fs1_var_2257) * fac01_var_2203(jl_var_2275, lay_var_2230)
          fac011_var_2240 = (1.0D0 - fs1_var_2257) * fac11_var_2205(jl_var_2275, lay_var_2230)
          fac101_var_2238 = fs1_var_2257 * fac01_var_2203(jl_var_2275, lay_var_2230)
          fac111_var_2241 = fs1_var_2257 * fac11_var_2205(jl_var_2275, lay_var_2230)
          fac201_var_2239 = 0.0D0
          fac211_var_2242 = 0.0D0
        END IF
        IF (specparm_var_2255 .LT. 0.125D0) THEN
          tau_major_var_2251(1 : ng16) = speccomb_var_2256 * (fac000_var_2231 * absa_var_169(ind0_var_2223, 1 : 2) + fac100_var_2232 * absa_var_169(ind0_var_2223 + 1, 1 : 2) + fac200_var_2233 * absa_var_169(ind0_var_2223 + 2, 1 : 2) + fac010_var_2234 * absa_var_169(ind0_var_2223 + 9, 1 : 2) + fac110_var_2235 * absa_var_169(ind0_var_2223 + 10, 1 : 2) + fac210_var_2236 * absa_var_169(ind0_var_2223 + 11, 1 : 2))
        ELSE IF (specparm_var_2255 .GT. 0.875D0) THEN
          tau_major_var_2251(1 : ng16) = speccomb_var_2256 * (fac200_var_2233 * absa_var_169(ind0_var_2223 - 1, 1 : 2) + fac100_var_2232 * absa_var_169(ind0_var_2223, 1 : 2) + fac000_var_2231 * absa_var_169(ind0_var_2223 + 1, 1 : 2) + fac210_var_2236 * absa_var_169(ind0_var_2223 + 8, 1 : 2) + fac110_var_2235 * absa_var_169(ind0_var_2223 + 9, 1 : 2) + fac010_var_2234 * absa_var_169(ind0_var_2223 + 10, 1 : 2))
        ELSE
          tau_major_var_2251(1 : ng16) = speccomb_var_2256 * (fac000_var_2231 * absa_var_169(ind0_var_2223, 1 : 2) + fac100_var_2232 * absa_var_169(ind0_var_2223 + 1, 1 : 2) + fac010_var_2234 * absa_var_169(ind0_var_2223 + 9, 1 : 2) + fac110_var_2235 * absa_var_169(ind0_var_2223 + 10, 1 : 2))
        END IF
        IF (specparm1_var_2259 .LT. 0.125D0) THEN
          tau_major1_var_2252(1 : ng16) = speccomb1_var_2260 * (fac001_var_2237 * absa_var_169(ind1_var_2224, 1 : 2) + fac101_var_2238 * absa_var_169(ind1_var_2224 + 1, 1 : 2) + fac201_var_2239 * absa_var_169(ind1_var_2224 + 2, 1 : 2) + fac011_var_2240 * absa_var_169(ind1_var_2224 + 9, 1 : 2) + fac111_var_2241 * absa_var_169(ind1_var_2224 + 10, 1 : 2) + fac211_var_2242 * absa_var_169(ind1_var_2224 + 11, 1 : 2))
        ELSE IF (specparm1_var_2259 .GT. 0.875D0) THEN
          tau_major1_var_2252(1 : ng16) = speccomb1_var_2260 * (fac201_var_2239 * absa_var_169(ind1_var_2224 - 1, 1 : 2) + fac101_var_2238 * absa_var_169(ind1_var_2224, 1 : 2) + fac001_var_2237 * absa_var_169(ind1_var_2224 + 1, 1 : 2) + fac211_var_2242 * absa_var_169(ind1_var_2224 + 8, 1 : 2) + fac111_var_2241 * absa_var_169(ind1_var_2224 + 9, 1 : 2) + fac011_var_2240 * absa_var_169(ind1_var_2224 + 10, 1 : 2))
        ELSE
          tau_major1_var_2252(1 : ng16) = speccomb1_var_2260 * (fac001_var_2237 * absa_var_169(ind1_var_2224, 1 : 2) + fac101_var_2238 * absa_var_169(ind1_var_2224 + 1, 1 : 2) + fac011_var_2240 * absa_var_169(ind1_var_2224 + 9, 1 : 2) + fac111_var_2241 * absa_var_169(ind1_var_2224 + 10, 1 : 2))
        END IF
        DO ig_var_2222 = 1, 2
          tauself_var_2250 = selffac_var_2213(jl_var_2275, lay_var_2230) * (selfref_var_171(inds_var_2225, ig_var_2222) + selffrac_var_2214(jl_var_2275, lay_var_2230) * (selfref_var_171(inds_var_2225 + 1, ig_var_2222) - selfref_var_171(inds_var_2225, ig_var_2222)))
          taufor_var_2249 = forfac_var_2220(jl_var_2275, lay_var_2230) * (forref_var_172(indf_var_2226, ig_var_2222) + forfrac_var_2221(jl_var_2275, lay_var_2230) * (forref_var_172(indf_var_2226 + 1, ig_var_2222) - forref_var_172(indf_var_2226, ig_var_2222)))
          taug_var_2200(jl_var_2275, 138 + ig_var_2222, lay_var_2230) = tau_major_var_2251(ig_var_2222) + tau_major1_var_2252(ig_var_2222) + tauself_var_2250 + taufor_var_2249
          fracs_var_2216(jl_var_2275, 138 + ig_var_2222, lay_var_2230) = fracrefa_var_167(ig_var_2222, jpl_var_2229) + fpl_var_2261 * (fracrefa_var_167(ig_var_2222, jpl_var_2229 + 1) - fracrefa_var_167(ig_var_2222, jpl_var_2229))
        END DO
      END DO
      ixc0_var_2272 = kfdia_var_2198 - kidia_var_2197 + 1 - ixc0_var_2272
      DO ixp_var_2273 = 1, ixc0_var_2272
        jl_var_2275 = ixhigh_var_2269(ixp_var_2273, lay_var_2230)
        ind0_var_2223 = ((jp_var_2206(jl_var_2275, lay_var_2230) - 13) * 5 + (jt_var_2207(jl_var_2275, lay_var_2230) - 1)) * nspb_var_237(16) + 1
        ind1_var_2224 = ((jp_var_2206(jl_var_2275, lay_var_2230) - 12) * 5 + (jt1_var_2208(jl_var_2275, lay_var_2230) - 1)) * nspb_var_237(16) + 1
        DO ig_var_2222 = 1, 2
          taug_var_2200(jl_var_2275, 138 + ig_var_2222, lay_var_2230) = colch4_var_2211(jl_var_2275, lay_var_2230) * (fac00_var_2202(jl_var_2275, lay_var_2230) * absb_var_170(ind0_var_2223, ig_var_2222) + fac10_var_2204(jl_var_2275, lay_var_2230) * absb_var_170(ind0_var_2223 + 1, ig_var_2222) + fac01_var_2203(jl_var_2275, lay_var_2230) * absb_var_170(ind1_var_2224, ig_var_2222) + fac11_var_2205(jl_var_2275, lay_var_2230) * absb_var_170(ind1_var_2224 + 1, ig_var_2222))
          fracs_var_2216(jl_var_2275, 138 + ig_var_2222, lay_var_2230) = fracrefb_var_168(ig_var_2222)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol16
SUBROUTINE rrtm_taumol12(kidia_var_2276, kfdia_var_2277, klev_var_2278, taug_var_2279, p_tauaerl_var_2280, fac00_var_2281, fac01_var_2282, fac10_var_2283, fac11_var_2284, forfac_var_2300, forfrac_var_2299, indfor_var_2298, jp_var_2285, jt_var_2286, jt1_var_2287, oneminus_var_2288, colh2o_var_2289, colco2_var_2290, laytrop_var_2291, selffac_var_2292, selffrac_var_2293, indself_var_2294, fracs_var_2295, rat_h2oco2_var_2296, rat_h2oco2_1_var_2297)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_236
  USE yoerrtm, ONLY: ng12
  USE yoerrta12, ONLY: absa_var_147, forref_var_149, fracrefa_var_146, selfref_var_148
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2276
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2277
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2278
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2279(kidia_var_2276 : kfdia_var_2277, 140, klev_var_2278)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2280(kidia_var_2276 : kfdia_var_2277, klev_var_2278, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2281(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2282(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2283(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2284(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2285(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2286(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2287(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2288
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2289(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2290(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2291(kidia_var_2276 : kfdia_var_2277)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2292(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2293(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2294(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2295(kidia_var_2276 : kfdia_var_2277, 140, klev_var_2278)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_2296(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_2297(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2298(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2299(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2300(kidia_var_2276 : kfdia_var_2277, klev_var_2278)
  REAL(KIND = 8) :: speccomb_var_2301, speccomb1_var_2302, speccomb_planck_var_2303
  INTEGER(KIND = 4) :: ind0_var_2304, ind1_var_2305, inds_var_2306, indf_var_2307
  INTEGER(KIND = 4) :: ig_var_2308, js_var_2309, lay_var_2310, js1_var_2311, jpl_var_2312
  REAL(KIND = 8) :: fs_var_2313, specmult_var_2314, specparm_var_2315, fs1_var_2316, specmult1_var_2317, specparm1_var_2318, fpl_var_2319, specmult_planck_var_2320, specparm_planck_var_2321
  REAL(KIND = 8) :: fac000_var_2322, fac100_var_2323, fac200_var_2324, fac010_var_2325, fac110_var_2326, fac210_var_2327, fac001_var_2328, fac101_var_2329, fac201_var_2330, fac011_var_2331, fac111_var_2332, fac211_var_2333
  REAL(KIND = 8) :: p_var_2334, p4_var_2335, fk0_var_2336, fk1_var_2337, fk2_var_2338
  REAL(KIND = 8) :: taufor_var_2339, tauself_var_2340, tau_major_var_2341(8), tau_major1_var_2342(8)
  REAL(KIND = 8) :: refrat_planck_a_var_2343
  INTEGER(KIND = 4) :: laytrop_min_var_2344, laytrop_max_var_2345
  INTEGER(KIND = 4) :: ixc_var_2346(klev_var_2278), ixlow_var_2347(kfdia_var_2277, klev_var_2278), ixhigh_var_2348(kfdia_var_2277, klev_var_2278)
  INTEGER(KIND = 4) :: ich_var_2349, icl_var_2350, ixc0_var_2351, ixp_var_2352, jc_var_2353, jl_var_2354
  laytrop_min_var_2344 = MINVAL(laytrop_var_2291)
  laytrop_max_var_2345 = MAXVAL(laytrop_var_2291)
  ixlow_var_2347 = 0
  ixhigh_var_2348 = 0
  ixc_var_2346 = 0
  DO lay_var_2310 = laytrop_min_var_2344 + 1, laytrop_max_var_2345
    icl_var_2350 = 0
    ich_var_2349 = 0
    DO jc_var_2353 = kidia_var_2276, kfdia_var_2277
      IF (lay_var_2310 <= laytrop_var_2291(jc_var_2353)) THEN
        icl_var_2350 = icl_var_2350 + 1
        ixlow_var_2347(icl_var_2350, lay_var_2310) = jc_var_2353
      ELSE
        ich_var_2349 = ich_var_2349 + 1
        ixhigh_var_2348(ich_var_2349, lay_var_2310) = jc_var_2353
      END IF
    END DO
    ixc_var_2346(lay_var_2310) = icl_var_2350
  END DO
  refrat_planck_a_var_2343 = chi_mls(1, 10) / chi_mls(2, 10)
  DO lay_var_2310 = 1, laytrop_min_var_2344
    DO jl_var_2354 = kidia_var_2276, kfdia_var_2277
      speccomb_var_2301 = colh2o_var_2289(jl_var_2354, lay_var_2310) + rat_h2oco2_var_2296(jl_var_2354, lay_var_2310) * colco2_var_2290(jl_var_2354, lay_var_2310)
      specparm_var_2315 = MIN(colh2o_var_2289(jl_var_2354, lay_var_2310) / speccomb_var_2301, oneminus_var_2288)
      specmult_var_2314 = 8.0D0 * (specparm_var_2315)
      js_var_2309 = 1 + INT(specmult_var_2314)
      fs_var_2313 = ((specmult_var_2314) - AINT((specmult_var_2314)))
      speccomb1_var_2302 = colh2o_var_2289(jl_var_2354, lay_var_2310) + rat_h2oco2_1_var_2297(jl_var_2354, lay_var_2310) * colco2_var_2290(jl_var_2354, lay_var_2310)
      specparm1_var_2318 = MIN(colh2o_var_2289(jl_var_2354, lay_var_2310) / speccomb1_var_2302, oneminus_var_2288)
      specmult1_var_2317 = 8.0D0 * (specparm1_var_2318)
      js1_var_2311 = 1 + INT(specmult1_var_2317)
      fs1_var_2316 = ((specmult1_var_2317) - AINT((specmult1_var_2317)))
      speccomb_planck_var_2303 = colh2o_var_2289(jl_var_2354, lay_var_2310) + refrat_planck_a_var_2343 * colco2_var_2290(jl_var_2354, lay_var_2310)
      specparm_planck_var_2321 = MIN(colh2o_var_2289(jl_var_2354, lay_var_2310) / speccomb_planck_var_2303, oneminus_var_2288)
      specmult_planck_var_2320 = 8.0D0 * specparm_planck_var_2321
      jpl_var_2312 = 1 + INT(specmult_planck_var_2320)
      fpl_var_2319 = ((specmult_planck_var_2320) - AINT((specmult_planck_var_2320)))
      ind0_var_2304 = ((jp_var_2285(jl_var_2354, lay_var_2310) - 1) * 5 + (jt_var_2286(jl_var_2354, lay_var_2310) - 1)) * nspa_var_236(12) + js_var_2309
      ind1_var_2305 = (jp_var_2285(jl_var_2354, lay_var_2310) * 5 + (jt1_var_2287(jl_var_2354, lay_var_2310) - 1)) * nspa_var_236(12) + js1_var_2311
      inds_var_2306 = indself_var_2294(jl_var_2354, lay_var_2310)
      indf_var_2307 = indfor_var_2298(jl_var_2354, lay_var_2310)
      IF (specparm_var_2315 .LT. 0.125D0) THEN
        p_var_2334 = fs_var_2313 - 1.0D0
        p4_var_2335 = p_var_2334 ** 4
        fk0_var_2336 = p4_var_2335
        fk1_var_2337 = 1.0D0 - p_var_2334 - 2.0D0 * p4_var_2335
        fk2_var_2338 = p_var_2334 + p4_var_2335
        fac000_var_2322 = fk0_var_2336 * fac00_var_2281(jl_var_2354, lay_var_2310)
        fac100_var_2323 = fk1_var_2337 * fac00_var_2281(jl_var_2354, lay_var_2310)
        fac200_var_2324 = fk2_var_2338 * fac00_var_2281(jl_var_2354, lay_var_2310)
        fac010_var_2325 = fk0_var_2336 * fac10_var_2283(jl_var_2354, lay_var_2310)
        fac110_var_2326 = fk1_var_2337 * fac10_var_2283(jl_var_2354, lay_var_2310)
        fac210_var_2327 = fk2_var_2338 * fac10_var_2283(jl_var_2354, lay_var_2310)
      ELSE IF (specparm_var_2315 .GT. 0.875D0) THEN
        p_var_2334 = - fs_var_2313
        p4_var_2335 = p_var_2334 ** 4
        fk0_var_2336 = p4_var_2335
        fk1_var_2337 = 1.0D0 - p_var_2334 - 2.0D0 * p4_var_2335
        fk2_var_2338 = p_var_2334 + p4_var_2335
        fac000_var_2322 = fk0_var_2336 * fac00_var_2281(jl_var_2354, lay_var_2310)
        fac100_var_2323 = fk1_var_2337 * fac00_var_2281(jl_var_2354, lay_var_2310)
        fac200_var_2324 = fk2_var_2338 * fac00_var_2281(jl_var_2354, lay_var_2310)
        fac010_var_2325 = fk0_var_2336 * fac10_var_2283(jl_var_2354, lay_var_2310)
        fac110_var_2326 = fk1_var_2337 * fac10_var_2283(jl_var_2354, lay_var_2310)
        fac210_var_2327 = fk2_var_2338 * fac10_var_2283(jl_var_2354, lay_var_2310)
      ELSE
        fac000_var_2322 = (1.0D0 - fs_var_2313) * fac00_var_2281(jl_var_2354, lay_var_2310)
        fac010_var_2325 = (1.0D0 - fs_var_2313) * fac10_var_2283(jl_var_2354, lay_var_2310)
        fac100_var_2323 = fs_var_2313 * fac00_var_2281(jl_var_2354, lay_var_2310)
        fac110_var_2326 = fs_var_2313 * fac10_var_2283(jl_var_2354, lay_var_2310)
        fac200_var_2324 = 0.0D0
        fac210_var_2327 = 0.0D0
      END IF
      IF (specparm1_var_2318 .LT. 0.125D0) THEN
        p_var_2334 = fs1_var_2316 - 1.0D0
        p4_var_2335 = p_var_2334 ** 4
        fk0_var_2336 = p4_var_2335
        fk1_var_2337 = 1.0D0 - p_var_2334 - 2.0D0 * p4_var_2335
        fk2_var_2338 = p_var_2334 + p4_var_2335
        fac001_var_2328 = fk0_var_2336 * fac01_var_2282(jl_var_2354, lay_var_2310)
        fac101_var_2329 = fk1_var_2337 * fac01_var_2282(jl_var_2354, lay_var_2310)
        fac201_var_2330 = fk2_var_2338 * fac01_var_2282(jl_var_2354, lay_var_2310)
        fac011_var_2331 = fk0_var_2336 * fac11_var_2284(jl_var_2354, lay_var_2310)
        fac111_var_2332 = fk1_var_2337 * fac11_var_2284(jl_var_2354, lay_var_2310)
        fac211_var_2333 = fk2_var_2338 * fac11_var_2284(jl_var_2354, lay_var_2310)
      ELSE IF (specparm1_var_2318 .GT. 0.875D0) THEN
        p_var_2334 = - fs1_var_2316
        p4_var_2335 = p_var_2334 ** 4
        fk0_var_2336 = p4_var_2335
        fk1_var_2337 = 1.0D0 - p_var_2334 - 2.0D0 * p4_var_2335
        fk2_var_2338 = p_var_2334 + p4_var_2335
        fac001_var_2328 = fk0_var_2336 * fac01_var_2282(jl_var_2354, lay_var_2310)
        fac101_var_2329 = fk1_var_2337 * fac01_var_2282(jl_var_2354, lay_var_2310)
        fac201_var_2330 = fk2_var_2338 * fac01_var_2282(jl_var_2354, lay_var_2310)
        fac011_var_2331 = fk0_var_2336 * fac11_var_2284(jl_var_2354, lay_var_2310)
        fac111_var_2332 = fk1_var_2337 * fac11_var_2284(jl_var_2354, lay_var_2310)
        fac211_var_2333 = fk2_var_2338 * fac11_var_2284(jl_var_2354, lay_var_2310)
      ELSE
        fac001_var_2328 = (1.0D0 - fs1_var_2316) * fac01_var_2282(jl_var_2354, lay_var_2310)
        fac011_var_2331 = (1.0D0 - fs1_var_2316) * fac11_var_2284(jl_var_2354, lay_var_2310)
        fac101_var_2329 = fs1_var_2316 * fac01_var_2282(jl_var_2354, lay_var_2310)
        fac111_var_2332 = fs1_var_2316 * fac11_var_2284(jl_var_2354, lay_var_2310)
        fac201_var_2330 = 0.0D0
        fac211_var_2333 = 0.0D0
      END IF
      IF (specparm_var_2315 .LT. 0.125D0) THEN
        tau_major_var_2341(1 : ng12) = speccomb_var_2301 * (fac000_var_2322 * absa_var_147(ind0_var_2304, 1 : 8) + fac100_var_2323 * absa_var_147(ind0_var_2304 + 1, 1 : 8) + fac200_var_2324 * absa_var_147(ind0_var_2304 + 2, 1 : 8) + fac010_var_2325 * absa_var_147(ind0_var_2304 + 9, 1 : 8) + fac110_var_2326 * absa_var_147(ind0_var_2304 + 10, 1 : 8) + fac210_var_2327 * absa_var_147(ind0_var_2304 + 11, 1 : 8))
      ELSE IF (specparm_var_2315 .GT. 0.875D0) THEN
        tau_major_var_2341(1 : ng12) = speccomb_var_2301 * (fac200_var_2324 * absa_var_147(ind0_var_2304 - 1, 1 : 8) + fac100_var_2323 * absa_var_147(ind0_var_2304, 1 : 8) + fac000_var_2322 * absa_var_147(ind0_var_2304 + 1, 1 : 8) + fac210_var_2327 * absa_var_147(ind0_var_2304 + 8, 1 : 8) + fac110_var_2326 * absa_var_147(ind0_var_2304 + 9, 1 : 8) + fac010_var_2325 * absa_var_147(ind0_var_2304 + 10, 1 : 8))
      ELSE
        tau_major_var_2341(1 : ng12) = speccomb_var_2301 * (fac000_var_2322 * absa_var_147(ind0_var_2304, 1 : 8) + fac100_var_2323 * absa_var_147(ind0_var_2304 + 1, 1 : 8) + fac010_var_2325 * absa_var_147(ind0_var_2304 + 9, 1 : 8) + fac110_var_2326 * absa_var_147(ind0_var_2304 + 10, 1 : 8))
      END IF
      IF (specparm1_var_2318 .LT. 0.125D0) THEN
        tau_major1_var_2342(1 : ng12) = speccomb1_var_2302 * (fac001_var_2328 * absa_var_147(ind1_var_2305, 1 : 8) + fac101_var_2329 * absa_var_147(ind1_var_2305 + 1, 1 : 8) + fac201_var_2330 * absa_var_147(ind1_var_2305 + 2, 1 : 8) + fac011_var_2331 * absa_var_147(ind1_var_2305 + 9, 1 : 8) + fac111_var_2332 * absa_var_147(ind1_var_2305 + 10, 1 : 8) + fac211_var_2333 * absa_var_147(ind1_var_2305 + 11, 1 : 8))
      ELSE IF (specparm1_var_2318 .GT. 0.875D0) THEN
        tau_major1_var_2342(1 : ng12) = speccomb1_var_2302 * (fac201_var_2330 * absa_var_147(ind1_var_2305 - 1, 1 : 8) + fac101_var_2329 * absa_var_147(ind1_var_2305, 1 : 8) + fac001_var_2328 * absa_var_147(ind1_var_2305 + 1, 1 : 8) + fac211_var_2333 * absa_var_147(ind1_var_2305 + 8, 1 : 8) + fac111_var_2332 * absa_var_147(ind1_var_2305 + 9, 1 : 8) + fac011_var_2331 * absa_var_147(ind1_var_2305 + 10, 1 : 8))
      ELSE
        tau_major1_var_2342(1 : ng12) = speccomb1_var_2302 * (fac001_var_2328 * absa_var_147(ind1_var_2305, 1 : 8) + fac101_var_2329 * absa_var_147(ind1_var_2305 + 1, 1 : 8) + fac011_var_2331 * absa_var_147(ind1_var_2305 + 9, 1 : 8) + fac111_var_2332 * absa_var_147(ind1_var_2305 + 10, 1 : 8))
      END IF
      DO ig_var_2308 = 1, 8
        tauself_var_2340 = selffac_var_2292(jl_var_2354, lay_var_2310) * (selfref_var_148(inds_var_2306, ig_var_2308) + selffrac_var_2293(jl_var_2354, lay_var_2310) * (selfref_var_148(inds_var_2306 + 1, ig_var_2308) - selfref_var_148(inds_var_2306, ig_var_2308)))
        taufor_var_2339 = forfac_var_2300(jl_var_2354, lay_var_2310) * (forref_var_149(indf_var_2307, ig_var_2308) + forfrac_var_2299(jl_var_2354, lay_var_2310) * (forref_var_149(indf_var_2307 + 1, ig_var_2308) - forref_var_149(indf_var_2307, ig_var_2308)))
        taug_var_2279(jl_var_2354, 122 + ig_var_2308, lay_var_2310) = tau_major_var_2341(ig_var_2308) + tau_major1_var_2342(ig_var_2308) + tauself_var_2340 + taufor_var_2339
        fracs_var_2295(jl_var_2354, 122 + ig_var_2308, lay_var_2310) = fracrefa_var_146(ig_var_2308, jpl_var_2312) + fpl_var_2319 * (fracrefa_var_146(ig_var_2308, jpl_var_2312 + 1) - fracrefa_var_146(ig_var_2308, jpl_var_2312))
      END DO
    END DO
  END DO
  DO ig_var_2308 = 1, 8
    DO lay_var_2310 = laytrop_max_var_2345 + 1, klev_var_2278
      DO jl_var_2354 = kidia_var_2276, kfdia_var_2277
        taug_var_2279(jl_var_2354, 122 + ig_var_2308, lay_var_2310) = 0.0D0
        fracs_var_2295(jl_var_2354, 122 + ig_var_2308, lay_var_2310) = 0.0D0
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2345 /= laytrop_min_var_2344) THEN
    DO lay_var_2310 = laytrop_min_var_2344 + 1, laytrop_max_var_2345
      ixc0_var_2351 = ixc_var_2346(lay_var_2310)
      DO ixp_var_2352 = 1, ixc0_var_2351
        jl_var_2354 = ixlow_var_2347(ixp_var_2352, lay_var_2310)
        speccomb_var_2301 = colh2o_var_2289(jl_var_2354, lay_var_2310) + rat_h2oco2_var_2296(jl_var_2354, lay_var_2310) * colco2_var_2290(jl_var_2354, lay_var_2310)
        specparm_var_2315 = MIN(colh2o_var_2289(jl_var_2354, lay_var_2310) / speccomb_var_2301, oneminus_var_2288)
        specmult_var_2314 = 8.0D0 * (specparm_var_2315)
        js_var_2309 = 1 + INT(specmult_var_2314)
        fs_var_2313 = ((specmult_var_2314) - AINT((specmult_var_2314)))
        speccomb1_var_2302 = colh2o_var_2289(jl_var_2354, lay_var_2310) + rat_h2oco2_1_var_2297(jl_var_2354, lay_var_2310) * colco2_var_2290(jl_var_2354, lay_var_2310)
        specparm1_var_2318 = MIN(colh2o_var_2289(jl_var_2354, lay_var_2310) / speccomb1_var_2302, oneminus_var_2288)
        specmult1_var_2317 = 8.0D0 * (specparm1_var_2318)
        js1_var_2311 = 1 + INT(specmult1_var_2317)
        fs1_var_2316 = ((specmult1_var_2317) - AINT((specmult1_var_2317)))
        speccomb_planck_var_2303 = colh2o_var_2289(jl_var_2354, lay_var_2310) + refrat_planck_a_var_2343 * colco2_var_2290(jl_var_2354, lay_var_2310)
        specparm_planck_var_2321 = MIN(colh2o_var_2289(jl_var_2354, lay_var_2310) / speccomb_planck_var_2303, oneminus_var_2288)
        specmult_planck_var_2320 = 8.0D0 * specparm_planck_var_2321
        jpl_var_2312 = 1 + INT(specmult_planck_var_2320)
        fpl_var_2319 = ((specmult_planck_var_2320) - AINT((specmult_planck_var_2320)))
        ind0_var_2304 = ((jp_var_2285(jl_var_2354, lay_var_2310) - 1) * 5 + (jt_var_2286(jl_var_2354, lay_var_2310) - 1)) * nspa_var_236(12) + js_var_2309
        ind1_var_2305 = (jp_var_2285(jl_var_2354, lay_var_2310) * 5 + (jt1_var_2287(jl_var_2354, lay_var_2310) - 1)) * nspa_var_236(12) + js1_var_2311
        inds_var_2306 = indself_var_2294(jl_var_2354, lay_var_2310)
        indf_var_2307 = indfor_var_2298(jl_var_2354, lay_var_2310)
        IF (specparm_var_2315 .LT. 0.125D0) THEN
          p_var_2334 = fs_var_2313 - 1.0D0
          p4_var_2335 = p_var_2334 ** 4
          fk0_var_2336 = p4_var_2335
          fk1_var_2337 = 1.0D0 - p_var_2334 - 2.0D0 * p4_var_2335
          fk2_var_2338 = p_var_2334 + p4_var_2335
          fac000_var_2322 = fk0_var_2336 * fac00_var_2281(jl_var_2354, lay_var_2310)
          fac100_var_2323 = fk1_var_2337 * fac00_var_2281(jl_var_2354, lay_var_2310)
          fac200_var_2324 = fk2_var_2338 * fac00_var_2281(jl_var_2354, lay_var_2310)
          fac010_var_2325 = fk0_var_2336 * fac10_var_2283(jl_var_2354, lay_var_2310)
          fac110_var_2326 = fk1_var_2337 * fac10_var_2283(jl_var_2354, lay_var_2310)
          fac210_var_2327 = fk2_var_2338 * fac10_var_2283(jl_var_2354, lay_var_2310)
        ELSE IF (specparm_var_2315 .GT. 0.875D0) THEN
          p_var_2334 = - fs_var_2313
          p4_var_2335 = p_var_2334 ** 4
          fk0_var_2336 = p4_var_2335
          fk1_var_2337 = 1.0D0 - p_var_2334 - 2.0D0 * p4_var_2335
          fk2_var_2338 = p_var_2334 + p4_var_2335
          fac000_var_2322 = fk0_var_2336 * fac00_var_2281(jl_var_2354, lay_var_2310)
          fac100_var_2323 = fk1_var_2337 * fac00_var_2281(jl_var_2354, lay_var_2310)
          fac200_var_2324 = fk2_var_2338 * fac00_var_2281(jl_var_2354, lay_var_2310)
          fac010_var_2325 = fk0_var_2336 * fac10_var_2283(jl_var_2354, lay_var_2310)
          fac110_var_2326 = fk1_var_2337 * fac10_var_2283(jl_var_2354, lay_var_2310)
          fac210_var_2327 = fk2_var_2338 * fac10_var_2283(jl_var_2354, lay_var_2310)
        ELSE
          fac000_var_2322 = (1.0D0 - fs_var_2313) * fac00_var_2281(jl_var_2354, lay_var_2310)
          fac010_var_2325 = (1.0D0 - fs_var_2313) * fac10_var_2283(jl_var_2354, lay_var_2310)
          fac100_var_2323 = fs_var_2313 * fac00_var_2281(jl_var_2354, lay_var_2310)
          fac110_var_2326 = fs_var_2313 * fac10_var_2283(jl_var_2354, lay_var_2310)
          fac200_var_2324 = 0.0D0
          fac210_var_2327 = 0.0D0
        END IF
        IF (specparm1_var_2318 .LT. 0.125D0) THEN
          p_var_2334 = fs1_var_2316 - 1.0D0
          p4_var_2335 = p_var_2334 ** 4
          fk0_var_2336 = p4_var_2335
          fk1_var_2337 = 1.0D0 - p_var_2334 - 2.0D0 * p4_var_2335
          fk2_var_2338 = p_var_2334 + p4_var_2335
          fac001_var_2328 = fk0_var_2336 * fac01_var_2282(jl_var_2354, lay_var_2310)
          fac101_var_2329 = fk1_var_2337 * fac01_var_2282(jl_var_2354, lay_var_2310)
          fac201_var_2330 = fk2_var_2338 * fac01_var_2282(jl_var_2354, lay_var_2310)
          fac011_var_2331 = fk0_var_2336 * fac11_var_2284(jl_var_2354, lay_var_2310)
          fac111_var_2332 = fk1_var_2337 * fac11_var_2284(jl_var_2354, lay_var_2310)
          fac211_var_2333 = fk2_var_2338 * fac11_var_2284(jl_var_2354, lay_var_2310)
        ELSE IF (specparm1_var_2318 .GT. 0.875D0) THEN
          p_var_2334 = - fs1_var_2316
          p4_var_2335 = p_var_2334 ** 4
          fk0_var_2336 = p4_var_2335
          fk1_var_2337 = 1.0D0 - p_var_2334 - 2.0D0 * p4_var_2335
          fk2_var_2338 = p_var_2334 + p4_var_2335
          fac001_var_2328 = fk0_var_2336 * fac01_var_2282(jl_var_2354, lay_var_2310)
          fac101_var_2329 = fk1_var_2337 * fac01_var_2282(jl_var_2354, lay_var_2310)
          fac201_var_2330 = fk2_var_2338 * fac01_var_2282(jl_var_2354, lay_var_2310)
          fac011_var_2331 = fk0_var_2336 * fac11_var_2284(jl_var_2354, lay_var_2310)
          fac111_var_2332 = fk1_var_2337 * fac11_var_2284(jl_var_2354, lay_var_2310)
          fac211_var_2333 = fk2_var_2338 * fac11_var_2284(jl_var_2354, lay_var_2310)
        ELSE
          fac001_var_2328 = (1.0D0 - fs1_var_2316) * fac01_var_2282(jl_var_2354, lay_var_2310)
          fac011_var_2331 = (1.0D0 - fs1_var_2316) * fac11_var_2284(jl_var_2354, lay_var_2310)
          fac101_var_2329 = fs1_var_2316 * fac01_var_2282(jl_var_2354, lay_var_2310)
          fac111_var_2332 = fs1_var_2316 * fac11_var_2284(jl_var_2354, lay_var_2310)
          fac201_var_2330 = 0.0D0
          fac211_var_2333 = 0.0D0
        END IF
        IF (specparm_var_2315 .LT. 0.125D0) THEN
          tau_major_var_2341(1 : ng12) = speccomb_var_2301 * (fac000_var_2322 * absa_var_147(ind0_var_2304, 1 : 8) + fac100_var_2323 * absa_var_147(ind0_var_2304 + 1, 1 : 8) + fac200_var_2324 * absa_var_147(ind0_var_2304 + 2, 1 : 8) + fac010_var_2325 * absa_var_147(ind0_var_2304 + 9, 1 : 8) + fac110_var_2326 * absa_var_147(ind0_var_2304 + 10, 1 : 8) + fac210_var_2327 * absa_var_147(ind0_var_2304 + 11, 1 : 8))
        ELSE IF (specparm_var_2315 .GT. 0.875D0) THEN
          tau_major_var_2341(1 : ng12) = speccomb_var_2301 * (fac200_var_2324 * absa_var_147(ind0_var_2304 - 1, 1 : 8) + fac100_var_2323 * absa_var_147(ind0_var_2304, 1 : 8) + fac000_var_2322 * absa_var_147(ind0_var_2304 + 1, 1 : 8) + fac210_var_2327 * absa_var_147(ind0_var_2304 + 8, 1 : 8) + fac110_var_2326 * absa_var_147(ind0_var_2304 + 9, 1 : 8) + fac010_var_2325 * absa_var_147(ind0_var_2304 + 10, 1 : 8))
        ELSE
          tau_major_var_2341(1 : ng12) = speccomb_var_2301 * (fac000_var_2322 * absa_var_147(ind0_var_2304, 1 : 8) + fac100_var_2323 * absa_var_147(ind0_var_2304 + 1, 1 : 8) + fac010_var_2325 * absa_var_147(ind0_var_2304 + 9, 1 : 8) + fac110_var_2326 * absa_var_147(ind0_var_2304 + 10, 1 : 8))
        END IF
        IF (specparm1_var_2318 .LT. 0.125D0) THEN
          tau_major1_var_2342(1 : ng12) = speccomb1_var_2302 * (fac001_var_2328 * absa_var_147(ind1_var_2305, 1 : 8) + fac101_var_2329 * absa_var_147(ind1_var_2305 + 1, 1 : 8) + fac201_var_2330 * absa_var_147(ind1_var_2305 + 2, 1 : 8) + fac011_var_2331 * absa_var_147(ind1_var_2305 + 9, 1 : 8) + fac111_var_2332 * absa_var_147(ind1_var_2305 + 10, 1 : 8) + fac211_var_2333 * absa_var_147(ind1_var_2305 + 11, 1 : 8))
        ELSE IF (specparm1_var_2318 .GT. 0.875D0) THEN
          tau_major1_var_2342(1 : ng12) = speccomb1_var_2302 * (fac201_var_2330 * absa_var_147(ind1_var_2305 - 1, 1 : 8) + fac101_var_2329 * absa_var_147(ind1_var_2305, 1 : 8) + fac001_var_2328 * absa_var_147(ind1_var_2305 + 1, 1 : 8) + fac211_var_2333 * absa_var_147(ind1_var_2305 + 8, 1 : 8) + fac111_var_2332 * absa_var_147(ind1_var_2305 + 9, 1 : 8) + fac011_var_2331 * absa_var_147(ind1_var_2305 + 10, 1 : 8))
        ELSE
          tau_major1_var_2342(1 : ng12) = speccomb1_var_2302 * (fac001_var_2328 * absa_var_147(ind1_var_2305, 1 : 8) + fac101_var_2329 * absa_var_147(ind1_var_2305 + 1, 1 : 8) + fac011_var_2331 * absa_var_147(ind1_var_2305 + 9, 1 : 8) + fac111_var_2332 * absa_var_147(ind1_var_2305 + 10, 1 : 8))
        END IF
        DO ig_var_2308 = 1, 8
          tauself_var_2340 = selffac_var_2292(jl_var_2354, lay_var_2310) * (selfref_var_148(inds_var_2306, ig_var_2308) + selffrac_var_2293(jl_var_2354, lay_var_2310) * (selfref_var_148(inds_var_2306 + 1, ig_var_2308) - selfref_var_148(inds_var_2306, ig_var_2308)))
          taufor_var_2339 = forfac_var_2300(jl_var_2354, lay_var_2310) * (forref_var_149(indf_var_2307, ig_var_2308) + forfrac_var_2299(jl_var_2354, lay_var_2310) * (forref_var_149(indf_var_2307 + 1, ig_var_2308) - forref_var_149(indf_var_2307, ig_var_2308)))
          taug_var_2279(jl_var_2354, 122 + ig_var_2308, lay_var_2310) = tau_major_var_2341(ig_var_2308) + tau_major1_var_2342(ig_var_2308) + tauself_var_2340 + taufor_var_2339
          fracs_var_2295(jl_var_2354, 122 + ig_var_2308, lay_var_2310) = fracrefa_var_146(ig_var_2308, jpl_var_2312) + fpl_var_2319 * (fracrefa_var_146(ig_var_2308, jpl_var_2312 + 1) - fracrefa_var_146(ig_var_2308, jpl_var_2312))
        END DO
      END DO
      ixc0_var_2351 = kfdia_var_2277 - kidia_var_2276 + 1 - ixc0_var_2351
      DO ig_var_2308 = 1, 8
        DO ixp_var_2352 = 1, ixc0_var_2351
          jl_var_2354 = ixhigh_var_2348(ixp_var_2352, lay_var_2310)
          taug_var_2279(jl_var_2354, 122 + ig_var_2308, lay_var_2310) = 0.0D0
          fracs_var_2295(jl_var_2354, 122 + ig_var_2308, lay_var_2310) = 0.0D0
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol12
SUBROUTINE rrtm_taumol13(kidia_var_2355, kfdia_var_2356, klev_var_2357, taug_var_2358, p_tauaerl_var_2359, fac00_var_2360, fac01_var_2361, fac10_var_2362, fac11_var_2363, forfac_var_2379, forfrac_var_2380, indfor_var_2378, jp_var_2364, jt_var_2365, jt1_var_2366, oneminus_var_2367, colh2o_var_2368, coln2o_var_2369, colco2_var_2370, colo3_var_2371, coldry_var_2372, laytrop_var_2373, selffac_var_2374, selffrac_var_2375, indself_var_2376, fracs_var_2377, rat_h2on2o, rat_h2on2o_1, minorfrac_var_2381, indminor_var_2382)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_236
  USE yoerrtm, ONLY: ng13
  USE yoerrta13, ONLY: absa_var_152, forref_var_154, fracrefa_var_150, fracrefb_var_151, ka_mco, ka_mco2_var_155, kb_mo3, selfref_var_153
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2355
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2356
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2357
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2358(kidia_var_2355 : kfdia_var_2356, 140, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2359(kidia_var_2355 : kfdia_var_2356, klev_var_2357, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2360(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2361(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2362(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2363(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2364(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2365(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2366(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2367
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2368(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_2369(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2370(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_2371(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_2372(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2373(kidia_var_2355 : kfdia_var_2356)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2374(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2375(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2376(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2377(kidia_var_2355 : kfdia_var_2356, 140, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: rat_h2on2o(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: rat_h2on2o_1(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2378(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2379(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2380(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2381(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2382(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  REAL(KIND = 8) :: speccomb_var_2383, speccomb1_var_2384, speccomb_planck_var_2385, speccomb_mco2_var_2386, speccomb_mco
  INTEGER(KIND = 4) :: ind0_var_2387, ind1_var_2388, inds_var_2389, indf_var_2390, indm_var_2391
  INTEGER(KIND = 4) :: ig_var_2392, js_var_2393, lay_var_2394, js1_var_2395, jpl_var_2396, jmco2_var_2397, jmco
  REAL(KIND = 8) :: refrat_planck_a_var_2398, refrat_m_a_var_2399, refrat_m_a3
  REAL(KIND = 8) :: fac000_var_2400, fac100_var_2401, fac200_var_2402, fac010_var_2403, fac110_var_2404, fac210_var_2405, fac001_var_2406, fac101_var_2407, fac201_var_2408, fac011_var_2409, fac111_var_2410, fac211_var_2411
  REAL(KIND = 8) :: p_var_2412, p4_var_2413, fk0_var_2414, fk1_var_2415, fk2_var_2416
  REAL(KIND = 8) :: taufor_var_2417, tauself_var_2418, tau_major_var_2419(4), tau_major1_var_2420(4), co2m1_var_2421, co2m2_var_2422, absco2_var_2423
  REAL(KIND = 8) :: com1, com2, absco, abso3_var_2424
  REAL(KIND = 8) :: chi_co2_var_2425, ratco2_var_2426, adjfac_var_2427, adjcolco2_var_2428
  REAL(KIND = 8) :: fs_var_2429, specmult_var_2430, specparm_var_2431, fs1_var_2432, specmult1_var_2433, specparm1_var_2434, fmco2_var_2435, specmult_mco2_var_2436, specparm_mco2_var_2437, fmco, specmult_mco, specparm_mco, fpl_var_2438, specmult_planck_var_2439, specparm_planck_var_2440
  REAL(KIND = 8) :: colco(kidia_var_2355 : kfdia_var_2356, klev_var_2357)
  INTEGER(KIND = 4) :: laytrop_min_var_2441, laytrop_max_var_2442
  INTEGER(KIND = 4) :: ixc_var_2443(klev_var_2357), ixlow_var_2444(kfdia_var_2356, klev_var_2357), ixhigh_var_2445(kfdia_var_2356, klev_var_2357)
  INTEGER(KIND = 4) :: ich_var_2446, icl_var_2447, ixc0_var_2448, ixp_var_2449, jc_var_2450, jl_var_2451
  DO lay_var_2394 = 1, klev_var_2357
    DO jc_var_2450 = kidia_var_2355, kfdia_var_2356
      colco(jc_var_2450, lay_var_2394) = 0.0D0
    END DO
  END DO
  laytrop_min_var_2441 = MINVAL(laytrop_var_2373)
  laytrop_max_var_2442 = MAXVAL(laytrop_var_2373)
  ixlow_var_2444 = 0
  ixhigh_var_2445 = 0
  ixc_var_2443 = 0
  DO lay_var_2394 = laytrop_min_var_2441 + 1, laytrop_max_var_2442
    icl_var_2447 = 0
    ich_var_2446 = 0
    DO jc_var_2450 = kidia_var_2355, kfdia_var_2356
      IF (lay_var_2394 <= laytrop_var_2373(jc_var_2450)) THEN
        icl_var_2447 = icl_var_2447 + 1
        ixlow_var_2444(icl_var_2447, lay_var_2394) = jc_var_2450
      ELSE
        ich_var_2446 = ich_var_2446 + 1
        ixhigh_var_2445(ich_var_2446, lay_var_2394) = jc_var_2450
      END IF
    END DO
    ixc_var_2443(lay_var_2394) = icl_var_2447
  END DO
  refrat_planck_a_var_2398 = chi_mls(1, 5) / chi_mls(4, 5)
  refrat_m_a_var_2399 = chi_mls(1, 1) / chi_mls(4, 1)
  refrat_m_a3 = chi_mls(1, 3) / chi_mls(4, 3)
  DO lay_var_2394 = 1, laytrop_min_var_2441
    DO jl_var_2451 = kidia_var_2355, kfdia_var_2356
      speccomb_var_2383 = colh2o_var_2368(jl_var_2451, lay_var_2394) + rat_h2on2o(jl_var_2451, lay_var_2394) * coln2o_var_2369(jl_var_2451, lay_var_2394)
      specparm_var_2431 = MIN(colh2o_var_2368(jl_var_2451, lay_var_2394) / speccomb_var_2383, oneminus_var_2367)
      specmult_var_2430 = 8.0D0 * (specparm_var_2431)
      js_var_2393 = 1 + INT(specmult_var_2430)
      fs_var_2429 = ((specmult_var_2430) - AINT((specmult_var_2430)))
      speccomb1_var_2384 = colh2o_var_2368(jl_var_2451, lay_var_2394) + rat_h2on2o_1(jl_var_2451, lay_var_2394) * coln2o_var_2369(jl_var_2451, lay_var_2394)
      specparm1_var_2434 = MIN(colh2o_var_2368(jl_var_2451, lay_var_2394) / speccomb1_var_2384, oneminus_var_2367)
      specmult1_var_2433 = 8.0D0 * (specparm1_var_2434)
      js1_var_2395 = 1 + INT(specmult1_var_2433)
      fs1_var_2432 = ((specmult1_var_2433) - AINT((specmult1_var_2433)))
      speccomb_mco2_var_2386 = colh2o_var_2368(jl_var_2451, lay_var_2394) + refrat_m_a_var_2399 * coln2o_var_2369(jl_var_2451, lay_var_2394)
      specparm_mco2_var_2437 = MIN(colh2o_var_2368(jl_var_2451, lay_var_2394) / speccomb_mco2_var_2386, oneminus_var_2367)
      specmult_mco2_var_2436 = 8.0D0 * specparm_mco2_var_2437
      jmco2_var_2397 = 1 + INT(specmult_mco2_var_2436)
      fmco2_var_2435 = ((specmult_mco2_var_2436) - AINT((specmult_mco2_var_2436)))
      chi_co2_var_2425 = colco2_var_2370(jl_var_2451, lay_var_2394) / (coldry_var_2372(jl_var_2451, lay_var_2394))
      ratco2_var_2426 = 1D+20 * chi_co2_var_2425 / 0.000355D0
      IF (ratco2_var_2426 .GT. 3.0D0) THEN
        adjfac_var_2427 = 2.0D0 + (ratco2_var_2426 - 2.0D0) ** 0.68D0
        adjcolco2_var_2428 = adjfac_var_2427 * 0.000355D0 * coldry_var_2372(jl_var_2451, lay_var_2394) * 1D-20
      ELSE
        adjcolco2_var_2428 = colco2_var_2370(jl_var_2451, lay_var_2394)
      END IF
      speccomb_mco = colh2o_var_2368(jl_var_2451, lay_var_2394) + refrat_m_a3 * coln2o_var_2369(jl_var_2451, lay_var_2394)
      specparm_mco = MIN(colh2o_var_2368(jl_var_2451, lay_var_2394) / speccomb_mco, oneminus_var_2367)
      specmult_mco = 8.0D0 * specparm_mco
      jmco = 1 + INT(specmult_mco)
      fmco = ((specmult_mco) - AINT((specmult_mco)))
      speccomb_planck_var_2385 = colh2o_var_2368(jl_var_2451, lay_var_2394) + refrat_planck_a_var_2398 * coln2o_var_2369(jl_var_2451, lay_var_2394)
      specparm_planck_var_2440 = MIN(colh2o_var_2368(jl_var_2451, lay_var_2394) / speccomb_planck_var_2385, oneminus_var_2367)
      specmult_planck_var_2439 = 8.0D0 * specparm_planck_var_2440
      jpl_var_2396 = 1 + INT(specmult_planck_var_2439)
      fpl_var_2438 = ((specmult_planck_var_2439) - AINT((specmult_planck_var_2439)))
      ind0_var_2387 = ((jp_var_2364(jl_var_2451, lay_var_2394) - 1) * 5 + (jt_var_2365(jl_var_2451, lay_var_2394) - 1)) * nspa_var_236(13) + js_var_2393
      ind1_var_2388 = (jp_var_2364(jl_var_2451, lay_var_2394) * 5 + (jt1_var_2366(jl_var_2451, lay_var_2394) - 1)) * nspa_var_236(13) + js1_var_2395
      inds_var_2389 = indself_var_2376(jl_var_2451, lay_var_2394)
      indf_var_2390 = indfor_var_2378(jl_var_2451, lay_var_2394)
      indm_var_2391 = indminor_var_2382(jl_var_2451, lay_var_2394)
      IF (specparm_var_2431 .LT. 0.125D0) THEN
        p_var_2412 = fs_var_2429 - 1.0D0
        p4_var_2413 = p_var_2412 ** 4
        fk0_var_2414 = p4_var_2413
        fk1_var_2415 = 1.0D0 - p_var_2412 - 2.0D0 * p4_var_2413
        fk2_var_2416 = p_var_2412 + p4_var_2413
        fac000_var_2400 = fk0_var_2414 * fac00_var_2360(jl_var_2451, lay_var_2394)
        fac100_var_2401 = fk1_var_2415 * fac00_var_2360(jl_var_2451, lay_var_2394)
        fac200_var_2402 = fk2_var_2416 * fac00_var_2360(jl_var_2451, lay_var_2394)
        fac010_var_2403 = fk0_var_2414 * fac10_var_2362(jl_var_2451, lay_var_2394)
        fac110_var_2404 = fk1_var_2415 * fac10_var_2362(jl_var_2451, lay_var_2394)
        fac210_var_2405 = fk2_var_2416 * fac10_var_2362(jl_var_2451, lay_var_2394)
      ELSE IF (specparm_var_2431 .GT. 0.875D0) THEN
        p_var_2412 = - fs_var_2429
        p4_var_2413 = p_var_2412 ** 4
        fk0_var_2414 = p4_var_2413
        fk1_var_2415 = 1.0D0 - p_var_2412 - 2.0D0 * p4_var_2413
        fk2_var_2416 = p_var_2412 + p4_var_2413
        fac000_var_2400 = fk0_var_2414 * fac00_var_2360(jl_var_2451, lay_var_2394)
        fac100_var_2401 = fk1_var_2415 * fac00_var_2360(jl_var_2451, lay_var_2394)
        fac200_var_2402 = fk2_var_2416 * fac00_var_2360(jl_var_2451, lay_var_2394)
        fac010_var_2403 = fk0_var_2414 * fac10_var_2362(jl_var_2451, lay_var_2394)
        fac110_var_2404 = fk1_var_2415 * fac10_var_2362(jl_var_2451, lay_var_2394)
        fac210_var_2405 = fk2_var_2416 * fac10_var_2362(jl_var_2451, lay_var_2394)
      ELSE
        fac000_var_2400 = (1.0D0 - fs_var_2429) * fac00_var_2360(jl_var_2451, lay_var_2394)
        fac010_var_2403 = (1.0D0 - fs_var_2429) * fac10_var_2362(jl_var_2451, lay_var_2394)
        fac100_var_2401 = fs_var_2429 * fac00_var_2360(jl_var_2451, lay_var_2394)
        fac110_var_2404 = fs_var_2429 * fac10_var_2362(jl_var_2451, lay_var_2394)
        fac200_var_2402 = 0.0D0
        fac210_var_2405 = 0.0D0
      END IF
      IF (specparm1_var_2434 .LT. 0.125D0) THEN
        p_var_2412 = fs1_var_2432 - 1.0D0
        p4_var_2413 = p_var_2412 ** 4
        fk0_var_2414 = p4_var_2413
        fk1_var_2415 = 1.0D0 - p_var_2412 - 2.0D0 * p4_var_2413
        fk2_var_2416 = p_var_2412 + p4_var_2413
        fac001_var_2406 = fk0_var_2414 * fac01_var_2361(jl_var_2451, lay_var_2394)
        fac101_var_2407 = fk1_var_2415 * fac01_var_2361(jl_var_2451, lay_var_2394)
        fac201_var_2408 = fk2_var_2416 * fac01_var_2361(jl_var_2451, lay_var_2394)
        fac011_var_2409 = fk0_var_2414 * fac11_var_2363(jl_var_2451, lay_var_2394)
        fac111_var_2410 = fk1_var_2415 * fac11_var_2363(jl_var_2451, lay_var_2394)
        fac211_var_2411 = fk2_var_2416 * fac11_var_2363(jl_var_2451, lay_var_2394)
      ELSE IF (specparm1_var_2434 .GT. 0.875D0) THEN
        p_var_2412 = - fs1_var_2432
        p4_var_2413 = p_var_2412 ** 4
        fk0_var_2414 = p4_var_2413
        fk1_var_2415 = 1.0D0 - p_var_2412 - 2.0D0 * p4_var_2413
        fk2_var_2416 = p_var_2412 + p4_var_2413
        fac001_var_2406 = fk0_var_2414 * fac01_var_2361(jl_var_2451, lay_var_2394)
        fac101_var_2407 = fk1_var_2415 * fac01_var_2361(jl_var_2451, lay_var_2394)
        fac201_var_2408 = fk2_var_2416 * fac01_var_2361(jl_var_2451, lay_var_2394)
        fac011_var_2409 = fk0_var_2414 * fac11_var_2363(jl_var_2451, lay_var_2394)
        fac111_var_2410 = fk1_var_2415 * fac11_var_2363(jl_var_2451, lay_var_2394)
        fac211_var_2411 = fk2_var_2416 * fac11_var_2363(jl_var_2451, lay_var_2394)
      ELSE
        fac001_var_2406 = (1.0D0 - fs1_var_2432) * fac01_var_2361(jl_var_2451, lay_var_2394)
        fac011_var_2409 = (1.0D0 - fs1_var_2432) * fac11_var_2363(jl_var_2451, lay_var_2394)
        fac101_var_2407 = fs1_var_2432 * fac01_var_2361(jl_var_2451, lay_var_2394)
        fac111_var_2410 = fs1_var_2432 * fac11_var_2363(jl_var_2451, lay_var_2394)
        fac201_var_2408 = 0.0D0
        fac211_var_2411 = 0.0D0
      END IF
      IF (specparm_var_2431 .LT. 0.125D0) THEN
        tau_major_var_2419(1 : ng13) = speccomb_var_2383 * (fac000_var_2400 * absa_var_152(ind0_var_2387, 1 : 4) + fac100_var_2401 * absa_var_152(ind0_var_2387 + 1, 1 : 4) + fac200_var_2402 * absa_var_152(ind0_var_2387 + 2, 1 : 4) + fac010_var_2403 * absa_var_152(ind0_var_2387 + 9, 1 : 4) + fac110_var_2404 * absa_var_152(ind0_var_2387 + 10, 1 : 4) + fac210_var_2405 * absa_var_152(ind0_var_2387 + 11, 1 : 4))
      ELSE IF (specparm_var_2431 .GT. 0.875D0) THEN
        tau_major_var_2419(1 : ng13) = speccomb_var_2383 * (fac200_var_2402 * absa_var_152(ind0_var_2387 - 1, 1 : 4) + fac100_var_2401 * absa_var_152(ind0_var_2387, 1 : 4) + fac000_var_2400 * absa_var_152(ind0_var_2387 + 1, 1 : 4) + fac210_var_2405 * absa_var_152(ind0_var_2387 + 8, 1 : 4) + fac110_var_2404 * absa_var_152(ind0_var_2387 + 9, 1 : 4) + fac010_var_2403 * absa_var_152(ind0_var_2387 + 10, 1 : 4))
      ELSE
        tau_major_var_2419(1 : ng13) = speccomb_var_2383 * (fac000_var_2400 * absa_var_152(ind0_var_2387, 1 : 4) + fac100_var_2401 * absa_var_152(ind0_var_2387 + 1, 1 : 4) + fac010_var_2403 * absa_var_152(ind0_var_2387 + 9, 1 : 4) + fac110_var_2404 * absa_var_152(ind0_var_2387 + 10, 1 : 4))
      END IF
      IF (specparm1_var_2434 .LT. 0.125D0) THEN
        tau_major1_var_2420(1 : ng13) = speccomb1_var_2384 * (fac001_var_2406 * absa_var_152(ind1_var_2388, 1 : 4) + fac101_var_2407 * absa_var_152(ind1_var_2388 + 1, 1 : 4) + fac201_var_2408 * absa_var_152(ind1_var_2388 + 2, 1 : 4) + fac011_var_2409 * absa_var_152(ind1_var_2388 + 9, 1 : 4) + fac111_var_2410 * absa_var_152(ind1_var_2388 + 10, 1 : 4) + fac211_var_2411 * absa_var_152(ind1_var_2388 + 11, 1 : 4))
      ELSE IF (specparm1_var_2434 .GT. 0.875D0) THEN
        tau_major1_var_2420(1 : ng13) = speccomb1_var_2384 * (fac201_var_2408 * absa_var_152(ind1_var_2388 - 1, 1 : 4) + fac101_var_2407 * absa_var_152(ind1_var_2388, 1 : 4) + fac001_var_2406 * absa_var_152(ind1_var_2388 + 1, 1 : 4) + fac211_var_2411 * absa_var_152(ind1_var_2388 + 8, 1 : 4) + fac111_var_2410 * absa_var_152(ind1_var_2388 + 9, 1 : 4) + fac011_var_2409 * absa_var_152(ind1_var_2388 + 10, 1 : 4))
      ELSE
        tau_major1_var_2420(1 : ng13) = speccomb1_var_2384 * (fac001_var_2406 * absa_var_152(ind1_var_2388, 1 : 4) + fac101_var_2407 * absa_var_152(ind1_var_2388 + 1, 1 : 4) + fac011_var_2409 * absa_var_152(ind1_var_2388 + 9, 1 : 4) + fac111_var_2410 * absa_var_152(ind1_var_2388 + 10, 1 : 4))
      END IF
      DO ig_var_2392 = 1, 4
        tauself_var_2418 = selffac_var_2374(jl_var_2451, lay_var_2394) * (selfref_var_153(inds_var_2389, ig_var_2392) + selffrac_var_2375(jl_var_2451, lay_var_2394) * (selfref_var_153(inds_var_2389 + 1, ig_var_2392) - selfref_var_153(inds_var_2389, ig_var_2392)))
        taufor_var_2417 = forfac_var_2379(jl_var_2451, lay_var_2394) * (forref_var_154(indf_var_2390, ig_var_2392) + forfrac_var_2380(jl_var_2451, lay_var_2394) * (forref_var_154(indf_var_2390 + 1, ig_var_2392) - forref_var_154(indf_var_2390, ig_var_2392)))
        co2m1_var_2421 = ka_mco2_var_155(jmco2_var_2397, indm_var_2391, ig_var_2392) + fmco2_var_2435 * (ka_mco2_var_155(jmco2_var_2397 + 1, indm_var_2391, ig_var_2392) - ka_mco2_var_155(jmco2_var_2397, indm_var_2391, ig_var_2392))
        co2m2_var_2422 = ka_mco2_var_155(jmco2_var_2397, indm_var_2391 + 1, ig_var_2392) + fmco2_var_2435 * (ka_mco2_var_155(jmco2_var_2397 + 1, indm_var_2391 + 1, ig_var_2392) - ka_mco2_var_155(jmco2_var_2397, indm_var_2391 + 1, ig_var_2392))
        absco2_var_2423 = co2m1_var_2421 + minorfrac_var_2381(jl_var_2451, lay_var_2394) * (co2m2_var_2422 - co2m1_var_2421)
        com1 = ka_mco(jmco, indm_var_2391, ig_var_2392) + fmco * (ka_mco(jmco + 1, indm_var_2391, ig_var_2392) - ka_mco(jmco, indm_var_2391, ig_var_2392))
        com2 = ka_mco(jmco, indm_var_2391 + 1, ig_var_2392) + fmco * (ka_mco(jmco + 1, indm_var_2391 + 1, ig_var_2392) - ka_mco(jmco, indm_var_2391 + 1, ig_var_2392))
        absco = com1 + minorfrac_var_2381(jl_var_2451, lay_var_2394) * (com2 - com1)
        taug_var_2358(jl_var_2451, 130 + ig_var_2392, lay_var_2394) = tau_major_var_2419(ig_var_2392) + tau_major1_var_2420(ig_var_2392) + tauself_var_2418 + taufor_var_2417 + adjcolco2_var_2428 * absco2_var_2423 + colco(jl_var_2451, lay_var_2394) * absco
        fracs_var_2377(jl_var_2451, 130 + ig_var_2392, lay_var_2394) = fracrefa_var_150(ig_var_2392, jpl_var_2396) + fpl_var_2438 * (fracrefa_var_150(ig_var_2392, jpl_var_2396 + 1) - fracrefa_var_150(ig_var_2392, jpl_var_2396))
      END DO
    END DO
  END DO
  DO lay_var_2394 = laytrop_max_var_2442 + 1, klev_var_2357
    DO jl_var_2451 = kidia_var_2355, kfdia_var_2356
      indm_var_2391 = indminor_var_2382(jl_var_2451, lay_var_2394)
      DO ig_var_2392 = 1, 4
        abso3_var_2424 = kb_mo3(indm_var_2391, ig_var_2392) + minorfrac_var_2381(jl_var_2451, lay_var_2394) * (kb_mo3(indm_var_2391 + 1, ig_var_2392) - kb_mo3(indm_var_2391, ig_var_2392))
        taug_var_2358(jl_var_2451, 130 + ig_var_2392, lay_var_2394) = colo3_var_2371(jl_var_2451, lay_var_2394) * abso3_var_2424
        fracs_var_2377(jl_var_2451, 130 + ig_var_2392, lay_var_2394) = fracrefb_var_151(ig_var_2392)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2442 /= laytrop_min_var_2441) THEN
    DO lay_var_2394 = laytrop_min_var_2441 + 1, laytrop_max_var_2442
      ixc0_var_2448 = ixc_var_2443(lay_var_2394)
      DO ixp_var_2449 = 1, ixc0_var_2448
        jl_var_2451 = ixlow_var_2444(ixp_var_2449, lay_var_2394)
        speccomb_var_2383 = colh2o_var_2368(jl_var_2451, lay_var_2394) + rat_h2on2o(jl_var_2451, lay_var_2394) * coln2o_var_2369(jl_var_2451, lay_var_2394)
        specparm_var_2431 = MIN(colh2o_var_2368(jl_var_2451, lay_var_2394) / speccomb_var_2383, oneminus_var_2367)
        specmult_var_2430 = 8.0D0 * (specparm_var_2431)
        js_var_2393 = 1 + INT(specmult_var_2430)
        fs_var_2429 = ((specmult_var_2430) - AINT((specmult_var_2430)))
        speccomb1_var_2384 = colh2o_var_2368(jl_var_2451, lay_var_2394) + rat_h2on2o_1(jl_var_2451, lay_var_2394) * coln2o_var_2369(jl_var_2451, lay_var_2394)
        specparm1_var_2434 = MIN(colh2o_var_2368(jl_var_2451, lay_var_2394) / speccomb1_var_2384, oneminus_var_2367)
        specmult1_var_2433 = 8.0D0 * (specparm1_var_2434)
        js1_var_2395 = 1 + INT(specmult1_var_2433)
        fs1_var_2432 = ((specmult1_var_2433) - AINT((specmult1_var_2433)))
        speccomb_mco2_var_2386 = colh2o_var_2368(jl_var_2451, lay_var_2394) + refrat_m_a_var_2399 * coln2o_var_2369(jl_var_2451, lay_var_2394)
        specparm_mco2_var_2437 = MIN(colh2o_var_2368(jl_var_2451, lay_var_2394) / speccomb_mco2_var_2386, oneminus_var_2367)
        specmult_mco2_var_2436 = 8.0D0 * specparm_mco2_var_2437
        jmco2_var_2397 = 1 + INT(specmult_mco2_var_2436)
        fmco2_var_2435 = ((specmult_mco2_var_2436) - AINT((specmult_mco2_var_2436)))
        chi_co2_var_2425 = colco2_var_2370(jl_var_2451, lay_var_2394) / (coldry_var_2372(jl_var_2451, lay_var_2394))
        ratco2_var_2426 = 1D+20 * chi_co2_var_2425 / 0.000355D0
        IF (ratco2_var_2426 .GT. 3.0D0) THEN
          adjfac_var_2427 = 2.0D0 + (ratco2_var_2426 - 2.0D0) ** 0.68D0
          adjcolco2_var_2428 = adjfac_var_2427 * 0.000355D0 * coldry_var_2372(jl_var_2451, lay_var_2394) * 1D-20
        ELSE
          adjcolco2_var_2428 = colco2_var_2370(jl_var_2451, lay_var_2394)
        END IF
        speccomb_mco = colh2o_var_2368(jl_var_2451, lay_var_2394) + refrat_m_a3 * coln2o_var_2369(jl_var_2451, lay_var_2394)
        specparm_mco = MIN(colh2o_var_2368(jl_var_2451, lay_var_2394) / speccomb_mco, oneminus_var_2367)
        specmult_mco = 8.0D0 * specparm_mco
        jmco = 1 + INT(specmult_mco)
        fmco = ((specmult_mco) - AINT((specmult_mco)))
        speccomb_planck_var_2385 = colh2o_var_2368(jl_var_2451, lay_var_2394) + refrat_planck_a_var_2398 * coln2o_var_2369(jl_var_2451, lay_var_2394)
        specparm_planck_var_2440 = MIN(colh2o_var_2368(jl_var_2451, lay_var_2394) / speccomb_planck_var_2385, oneminus_var_2367)
        specmult_planck_var_2439 = 8.0D0 * specparm_planck_var_2440
        jpl_var_2396 = 1 + INT(specmult_planck_var_2439)
        fpl_var_2438 = ((specmult_planck_var_2439) - AINT((specmult_planck_var_2439)))
        ind0_var_2387 = ((jp_var_2364(jl_var_2451, lay_var_2394) - 1) * 5 + (jt_var_2365(jl_var_2451, lay_var_2394) - 1)) * nspa_var_236(13) + js_var_2393
        ind1_var_2388 = (jp_var_2364(jl_var_2451, lay_var_2394) * 5 + (jt1_var_2366(jl_var_2451, lay_var_2394) - 1)) * nspa_var_236(13) + js1_var_2395
        inds_var_2389 = indself_var_2376(jl_var_2451, lay_var_2394)
        indf_var_2390 = indfor_var_2378(jl_var_2451, lay_var_2394)
        indm_var_2391 = indminor_var_2382(jl_var_2451, lay_var_2394)
        IF (specparm_var_2431 .LT. 0.125D0) THEN
          p_var_2412 = fs_var_2429 - 1.0D0
          p4_var_2413 = p_var_2412 ** 4
          fk0_var_2414 = p4_var_2413
          fk1_var_2415 = 1.0D0 - p_var_2412 - 2.0D0 * p4_var_2413
          fk2_var_2416 = p_var_2412 + p4_var_2413
          fac000_var_2400 = fk0_var_2414 * fac00_var_2360(jl_var_2451, lay_var_2394)
          fac100_var_2401 = fk1_var_2415 * fac00_var_2360(jl_var_2451, lay_var_2394)
          fac200_var_2402 = fk2_var_2416 * fac00_var_2360(jl_var_2451, lay_var_2394)
          fac010_var_2403 = fk0_var_2414 * fac10_var_2362(jl_var_2451, lay_var_2394)
          fac110_var_2404 = fk1_var_2415 * fac10_var_2362(jl_var_2451, lay_var_2394)
          fac210_var_2405 = fk2_var_2416 * fac10_var_2362(jl_var_2451, lay_var_2394)
        ELSE IF (specparm_var_2431 .GT. 0.875D0) THEN
          p_var_2412 = - fs_var_2429
          p4_var_2413 = p_var_2412 ** 4
          fk0_var_2414 = p4_var_2413
          fk1_var_2415 = 1.0D0 - p_var_2412 - 2.0D0 * p4_var_2413
          fk2_var_2416 = p_var_2412 + p4_var_2413
          fac000_var_2400 = fk0_var_2414 * fac00_var_2360(jl_var_2451, lay_var_2394)
          fac100_var_2401 = fk1_var_2415 * fac00_var_2360(jl_var_2451, lay_var_2394)
          fac200_var_2402 = fk2_var_2416 * fac00_var_2360(jl_var_2451, lay_var_2394)
          fac010_var_2403 = fk0_var_2414 * fac10_var_2362(jl_var_2451, lay_var_2394)
          fac110_var_2404 = fk1_var_2415 * fac10_var_2362(jl_var_2451, lay_var_2394)
          fac210_var_2405 = fk2_var_2416 * fac10_var_2362(jl_var_2451, lay_var_2394)
        ELSE
          fac000_var_2400 = (1.0D0 - fs_var_2429) * fac00_var_2360(jl_var_2451, lay_var_2394)
          fac010_var_2403 = (1.0D0 - fs_var_2429) * fac10_var_2362(jl_var_2451, lay_var_2394)
          fac100_var_2401 = fs_var_2429 * fac00_var_2360(jl_var_2451, lay_var_2394)
          fac110_var_2404 = fs_var_2429 * fac10_var_2362(jl_var_2451, lay_var_2394)
          fac200_var_2402 = 0.0D0
          fac210_var_2405 = 0.0D0
        END IF
        IF (specparm1_var_2434 .LT. 0.125D0) THEN
          p_var_2412 = fs1_var_2432 - 1.0D0
          p4_var_2413 = p_var_2412 ** 4
          fk0_var_2414 = p4_var_2413
          fk1_var_2415 = 1.0D0 - p_var_2412 - 2.0D0 * p4_var_2413
          fk2_var_2416 = p_var_2412 + p4_var_2413
          fac001_var_2406 = fk0_var_2414 * fac01_var_2361(jl_var_2451, lay_var_2394)
          fac101_var_2407 = fk1_var_2415 * fac01_var_2361(jl_var_2451, lay_var_2394)
          fac201_var_2408 = fk2_var_2416 * fac01_var_2361(jl_var_2451, lay_var_2394)
          fac011_var_2409 = fk0_var_2414 * fac11_var_2363(jl_var_2451, lay_var_2394)
          fac111_var_2410 = fk1_var_2415 * fac11_var_2363(jl_var_2451, lay_var_2394)
          fac211_var_2411 = fk2_var_2416 * fac11_var_2363(jl_var_2451, lay_var_2394)
        ELSE IF (specparm1_var_2434 .GT. 0.875D0) THEN
          p_var_2412 = - fs1_var_2432
          p4_var_2413 = p_var_2412 ** 4
          fk0_var_2414 = p4_var_2413
          fk1_var_2415 = 1.0D0 - p_var_2412 - 2.0D0 * p4_var_2413
          fk2_var_2416 = p_var_2412 + p4_var_2413
          fac001_var_2406 = fk0_var_2414 * fac01_var_2361(jl_var_2451, lay_var_2394)
          fac101_var_2407 = fk1_var_2415 * fac01_var_2361(jl_var_2451, lay_var_2394)
          fac201_var_2408 = fk2_var_2416 * fac01_var_2361(jl_var_2451, lay_var_2394)
          fac011_var_2409 = fk0_var_2414 * fac11_var_2363(jl_var_2451, lay_var_2394)
          fac111_var_2410 = fk1_var_2415 * fac11_var_2363(jl_var_2451, lay_var_2394)
          fac211_var_2411 = fk2_var_2416 * fac11_var_2363(jl_var_2451, lay_var_2394)
        ELSE
          fac001_var_2406 = (1.0D0 - fs1_var_2432) * fac01_var_2361(jl_var_2451, lay_var_2394)
          fac011_var_2409 = (1.0D0 - fs1_var_2432) * fac11_var_2363(jl_var_2451, lay_var_2394)
          fac101_var_2407 = fs1_var_2432 * fac01_var_2361(jl_var_2451, lay_var_2394)
          fac111_var_2410 = fs1_var_2432 * fac11_var_2363(jl_var_2451, lay_var_2394)
          fac201_var_2408 = 0.0D0
          fac211_var_2411 = 0.0D0
        END IF
        IF (specparm_var_2431 .LT. 0.125D0) THEN
          tau_major_var_2419(1 : ng13) = speccomb_var_2383 * (fac000_var_2400 * absa_var_152(ind0_var_2387, 1 : 4) + fac100_var_2401 * absa_var_152(ind0_var_2387 + 1, 1 : 4) + fac200_var_2402 * absa_var_152(ind0_var_2387 + 2, 1 : 4) + fac010_var_2403 * absa_var_152(ind0_var_2387 + 9, 1 : 4) + fac110_var_2404 * absa_var_152(ind0_var_2387 + 10, 1 : 4) + fac210_var_2405 * absa_var_152(ind0_var_2387 + 11, 1 : 4))
        ELSE IF (specparm_var_2431 .GT. 0.875D0) THEN
          tau_major_var_2419(1 : ng13) = speccomb_var_2383 * (fac200_var_2402 * absa_var_152(ind0_var_2387 - 1, 1 : 4) + fac100_var_2401 * absa_var_152(ind0_var_2387, 1 : 4) + fac000_var_2400 * absa_var_152(ind0_var_2387 + 1, 1 : 4) + fac210_var_2405 * absa_var_152(ind0_var_2387 + 8, 1 : 4) + fac110_var_2404 * absa_var_152(ind0_var_2387 + 9, 1 : 4) + fac010_var_2403 * absa_var_152(ind0_var_2387 + 10, 1 : 4))
        ELSE
          tau_major_var_2419(1 : ng13) = speccomb_var_2383 * (fac000_var_2400 * absa_var_152(ind0_var_2387, 1 : 4) + fac100_var_2401 * absa_var_152(ind0_var_2387 + 1, 1 : 4) + fac010_var_2403 * absa_var_152(ind0_var_2387 + 9, 1 : 4) + fac110_var_2404 * absa_var_152(ind0_var_2387 + 10, 1 : 4))
        END IF
        IF (specparm1_var_2434 .LT. 0.125D0) THEN
          tau_major1_var_2420(1 : ng13) = speccomb1_var_2384 * (fac001_var_2406 * absa_var_152(ind1_var_2388, 1 : 4) + fac101_var_2407 * absa_var_152(ind1_var_2388 + 1, 1 : 4) + fac201_var_2408 * absa_var_152(ind1_var_2388 + 2, 1 : 4) + fac011_var_2409 * absa_var_152(ind1_var_2388 + 9, 1 : 4) + fac111_var_2410 * absa_var_152(ind1_var_2388 + 10, 1 : 4) + fac211_var_2411 * absa_var_152(ind1_var_2388 + 11, 1 : 4))
        ELSE IF (specparm1_var_2434 .GT. 0.875D0) THEN
          tau_major1_var_2420(1 : ng13) = speccomb1_var_2384 * (fac201_var_2408 * absa_var_152(ind1_var_2388 - 1, 1 : 4) + fac101_var_2407 * absa_var_152(ind1_var_2388, 1 : 4) + fac001_var_2406 * absa_var_152(ind1_var_2388 + 1, 1 : 4) + fac211_var_2411 * absa_var_152(ind1_var_2388 + 8, 1 : 4) + fac111_var_2410 * absa_var_152(ind1_var_2388 + 9, 1 : 4) + fac011_var_2409 * absa_var_152(ind1_var_2388 + 10, 1 : 4))
        ELSE
          tau_major1_var_2420(1 : ng13) = speccomb1_var_2384 * (fac001_var_2406 * absa_var_152(ind1_var_2388, 1 : 4) + fac101_var_2407 * absa_var_152(ind1_var_2388 + 1, 1 : 4) + fac011_var_2409 * absa_var_152(ind1_var_2388 + 9, 1 : 4) + fac111_var_2410 * absa_var_152(ind1_var_2388 + 10, 1 : 4))
        END IF
        DO ig_var_2392 = 1, 4
          tauself_var_2418 = selffac_var_2374(jl_var_2451, lay_var_2394) * (selfref_var_153(inds_var_2389, ig_var_2392) + selffrac_var_2375(jl_var_2451, lay_var_2394) * (selfref_var_153(inds_var_2389 + 1, ig_var_2392) - selfref_var_153(inds_var_2389, ig_var_2392)))
          taufor_var_2417 = forfac_var_2379(jl_var_2451, lay_var_2394) * (forref_var_154(indf_var_2390, ig_var_2392) + forfrac_var_2380(jl_var_2451, lay_var_2394) * (forref_var_154(indf_var_2390 + 1, ig_var_2392) - forref_var_154(indf_var_2390, ig_var_2392)))
          co2m1_var_2421 = ka_mco2_var_155(jmco2_var_2397, indm_var_2391, ig_var_2392) + fmco2_var_2435 * (ka_mco2_var_155(jmco2_var_2397 + 1, indm_var_2391, ig_var_2392) - ka_mco2_var_155(jmco2_var_2397, indm_var_2391, ig_var_2392))
          co2m2_var_2422 = ka_mco2_var_155(jmco2_var_2397, indm_var_2391 + 1, ig_var_2392) + fmco2_var_2435 * (ka_mco2_var_155(jmco2_var_2397 + 1, indm_var_2391 + 1, ig_var_2392) - ka_mco2_var_155(jmco2_var_2397, indm_var_2391 + 1, ig_var_2392))
          absco2_var_2423 = co2m1_var_2421 + minorfrac_var_2381(jl_var_2451, lay_var_2394) * (co2m2_var_2422 - co2m1_var_2421)
          com1 = ka_mco(jmco, indm_var_2391, ig_var_2392) + fmco * (ka_mco(jmco + 1, indm_var_2391, ig_var_2392) - ka_mco(jmco, indm_var_2391, ig_var_2392))
          com2 = ka_mco(jmco, indm_var_2391 + 1, ig_var_2392) + fmco * (ka_mco(jmco + 1, indm_var_2391 + 1, ig_var_2392) - ka_mco(jmco, indm_var_2391 + 1, ig_var_2392))
          absco = com1 + minorfrac_var_2381(jl_var_2451, lay_var_2394) * (com2 - com1)
          taug_var_2358(jl_var_2451, 130 + ig_var_2392, lay_var_2394) = tau_major_var_2419(ig_var_2392) + tau_major1_var_2420(ig_var_2392) + tauself_var_2418 + taufor_var_2417 + adjcolco2_var_2428 * absco2_var_2423 + colco(jl_var_2451, lay_var_2394) * absco
          fracs_var_2377(jl_var_2451, 130 + ig_var_2392, lay_var_2394) = fracrefa_var_150(ig_var_2392, jpl_var_2396) + fpl_var_2438 * (fracrefa_var_150(ig_var_2392, jpl_var_2396 + 1) - fracrefa_var_150(ig_var_2392, jpl_var_2396))
        END DO
      END DO
      ixc0_var_2448 = kfdia_var_2356 - kidia_var_2355 + 1 - ixc0_var_2448
      DO ixp_var_2449 = 1, ixc0_var_2448
        jl_var_2451 = ixhigh_var_2445(ixp_var_2449, lay_var_2394)
        indm_var_2391 = indminor_var_2382(jl_var_2451, lay_var_2394)
        DO ig_var_2392 = 1, 4
          abso3_var_2424 = kb_mo3(indm_var_2391, ig_var_2392) + minorfrac_var_2381(jl_var_2451, lay_var_2394) * (kb_mo3(indm_var_2391 + 1, ig_var_2392) - kb_mo3(indm_var_2391, ig_var_2392))
          taug_var_2358(jl_var_2451, 130 + ig_var_2392, lay_var_2394) = colo3_var_2371(jl_var_2451, lay_var_2394) * abso3_var_2424
          fracs_var_2377(jl_var_2451, 130 + ig_var_2392, lay_var_2394) = fracrefb_var_151(ig_var_2392)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol13
SUBROUTINE global_init_fn
  USE yomlun_ifsaux, ONLY: nulout
  nulout = 6
END SUBROUTINE global_init_fn