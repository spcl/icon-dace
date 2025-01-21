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
    INTEGER(KIND = 4) :: jj_var_125, jk_var_126, in_var_127, ifilled
    IF (yd_stream_var_124 % inittest /= 12345678) CALL abor1('uniform_distribution called before initialize_random_numbers')
    in_var_127 = SIZE(px)
    ifilled = 0
    DO jj_var_125 = yd_stream_var_124 % iused + 1, MIN(607, in_var_127 + yd_stream_var_124 % iused)
      px(jj_var_125 - yd_stream_var_124 % iused) = yd_stream_var_124 % ix(jj_var_125) * yd_stream_var_124 % zrm
      ifilled = ifilled + 1
    END DO
    yd_stream_var_124 % iused = yd_stream_var_124 % iused + ifilled
    IF (ifilled == in_var_127) THEN
      RETURN
    END IF
    DO WHILE (ifilled < in_var_127)
      DO jj_var_125 = 1, 273
        yd_stream_var_124 % ix(jj_var_125) = IAND(1073741823, yd_stream_var_124 % ix(jj_var_125) + yd_stream_var_124 % ix(jj_var_125 - 273 + 607))
      END DO
      DO jk_var_126 = 1, 2
        DO jj_var_125 = 274 + (jk_var_126 - 1) * 167, MIN(607, 273 + jk_var_126 * 167)
          yd_stream_var_124 % ix(jj_var_125) = IAND(1073741823, yd_stream_var_124 % ix(jj_var_125) + yd_stream_var_124 % ix(jj_var_125 - 273))
        END DO
      END DO
      yd_stream_var_124 % iused = MIN(607, in_var_127 - ifilled)
      px(ifilled + 1 : ifilled + yd_stream_var_124 % iused) = yd_stream_var_124 % ix(1 : yd_stream_var_124 % iused) * yd_stream_var_124 % zrm
      ifilled = ifilled + yd_stream_var_124 % iused
    END DO
  END SUBROUTINE uniform_distribution
END MODULE random_numbers_mix
MODULE yoerrta1
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_128(10), fracrefb_var_129(10)
  REAL(KIND = 8) :: absa_var_130(65, 10)
  REAL(KIND = 8) :: absb_var_131(235, 10)
  REAL(KIND = 8) :: ka_mn2_var_132(19, 10), kb_mn2(19, 10)
  REAL(KIND = 8) :: selfref_var_133(10, 10), forref_var_134(4, 10)
END MODULE yoerrta1
MODULE yoerrta10
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(6) :: fracrefa_var_135
  REAL(KIND = 8), DIMENSION(6) :: fracrefb_var_136
  REAL(KIND = 8) :: absa_var_137(65, 6)
  REAL(KIND = 8) :: absb_var_138(235, 6)
  REAL(KIND = 8) :: selfref_var_139(10, 6)
  REAL(KIND = 8) :: forref_var_140(4, 6)
END MODULE yoerrta10
MODULE yoerrta11
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(8) :: fracrefa_var_141
  REAL(KIND = 8), DIMENSION(8) :: fracrefb_var_142
  REAL(KIND = 8) :: absa_var_143(65, 8)
  REAL(KIND = 8) :: absb_var_144(235, 8)
  REAL(KIND = 8) :: ka_mo2(19, 8)
  REAL(KIND = 8) :: kb_mo2(19, 8)
  REAL(KIND = 8) :: selfref_var_145(10, 8)
  REAL(KIND = 8) :: forref_var_146(4, 8)
END MODULE yoerrta11
MODULE yoerrta12
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_147(8, 9)
  REAL(KIND = 8) :: absa_var_148(585, 8)
  REAL(KIND = 8) :: selfref_var_149(10, 8)
  REAL(KIND = 8) :: forref_var_150(4, 8)
END MODULE yoerrta12
MODULE yoerrta13
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_151(4, 9)
  REAL(KIND = 8), DIMENSION(4) :: fracrefb_var_152
  REAL(KIND = 8) :: absa_var_153(585, 4)
  REAL(KIND = 8) :: selfref_var_154(10, 4)
  REAL(KIND = 8) :: forref_var_155(4, 4)
  REAL(KIND = 8) :: ka_mco2_var_156(9, 19, 4)
  REAL(KIND = 8) :: ka_mco(9, 19, 4)
  REAL(KIND = 8) :: kb_mo3(19, 4)
END MODULE yoerrta13
MODULE yoerrta14
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(2) :: fracrefa_var_157
  REAL(KIND = 8), DIMENSION(2) :: fracrefb_var_158
  REAL(KIND = 8) :: absa_var_159(65, 2)
  REAL(KIND = 8) :: absb_var_160(235, 2)
  REAL(KIND = 8) :: selfref_var_161(10, 2)
  REAL(KIND = 8) :: forref_var_162(4, 2)
END MODULE yoerrta14
MODULE yoerrta15
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_163(2, 9)
  REAL(KIND = 8) :: absa_var_164(585, 2)
  REAL(KIND = 8) :: ka_mn2_var_165(9, 19, 2)
  REAL(KIND = 8) :: selfref_var_166(10, 2)
  REAL(KIND = 8) :: forref_var_167(4, 2)
END MODULE yoerrta15
MODULE yoerrta16
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_168(2, 9)
  REAL(KIND = 8), DIMENSION(2) :: fracrefb_var_169
  REAL(KIND = 8) :: absa_var_170(585, 2)
  REAL(KIND = 8) :: absb_var_171(235, 2)
  REAL(KIND = 8) :: selfref_var_172(10, 2)
  REAL(KIND = 8) :: forref_var_173(4, 2)
END MODULE yoerrta16
MODULE yoerrta2
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_174(12), fracrefb_var_175(12)
  REAL(KIND = 8) :: absa_var_176(65, 12)
  REAL(KIND = 8) :: absb_var_177(235, 12)
  REAL(KIND = 8) :: selfref_var_178(10, 12), forref_var_179(4, 12)
END MODULE yoerrta2
MODULE yoerrta3
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_180(16, 9), fracrefb_var_181(16, 5)
  REAL(KIND = 8) :: ka_mn2o_var_182(9, 19, 16), kb_mn2o_var_183(5, 19, 16)
  REAL(KIND = 8) :: absa_var_184(585, 16)
  REAL(KIND = 8) :: absb_var_185(1175, 16)
  REAL(KIND = 8) :: selfref_var_186(10, 16)
  REAL(KIND = 8) :: forref_var_187(4, 16)
END MODULE yoerrta3
MODULE yoerrta4
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_188(14, 9), fracrefb_var_189(14, 5)
  REAL(KIND = 8) :: absa_var_190(585, 14)
  REAL(KIND = 8) :: absb_var_191(1175, 14)
  REAL(KIND = 8) :: selfref_var_192(10, 14), forref_var_193(4, 14)
END MODULE yoerrta4
MODULE yoerrta5
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_194(16, 9), fracrefb_var_195(16, 5)
  REAL(KIND = 8), DIMENSION(16) :: ccl4
  REAL(KIND = 8) :: absa_var_196(585, 16)
  REAL(KIND = 8) :: absb_var_197(1175, 16)
  REAL(KIND = 8) :: ka_mo3_var_198(9, 19, 16)
  REAL(KIND = 8) :: selfref_var_199(10, 16)
  REAL(KIND = 8) :: forref_var_200(4, 16)
END MODULE yoerrta5
MODULE yoerrta6
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(8) :: fracrefa_var_201
  REAL(KIND = 8), DIMENSION(8) :: cfc11adj
  REAL(KIND = 8), DIMENSION(8) :: cfc12_var_202
  REAL(KIND = 8) :: absa_var_203(65, 8)
  REAL(KIND = 8) :: selfref_var_204(10, 8)
  REAL(KIND = 8) :: ka_mco2_var_205(19, 8)
  REAL(KIND = 8) :: forref_var_206(4, 8)
END MODULE yoerrta6
MODULE yoerrta7
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_207(12, 9)
  REAL(KIND = 8), DIMENSION(12) :: fracrefb_var_208
  REAL(KIND = 8) :: absa_var_209(585, 12)
  REAL(KIND = 8) :: absb_var_210(235, 12)
  REAL(KIND = 8) :: selfref_var_211(10, 12)
  REAL(KIND = 8) :: ka_mco2_var_212(9, 19, 12)
  REAL(KIND = 8) :: kb_mco2_var_213(19, 12)
  REAL(KIND = 8) :: forref_var_214(4, 12)
END MODULE yoerrta7
MODULE yoerrta8
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8), DIMENSION(8) :: fracrefa_var_215
  REAL(KIND = 8), DIMENSION(8) :: fracrefb_var_216
  REAL(KIND = 8), DIMENSION(8) :: cfc12_var_217
  REAL(KIND = 8), DIMENSION(8) :: cfc22adj
  REAL(KIND = 8) :: absa_var_218(65, 8)
  REAL(KIND = 8) :: absb_var_219(235, 8)
  REAL(KIND = 8) :: ka_mco2_var_220(19, 8)
  REAL(KIND = 8) :: ka_mn2o_var_221(19, 8)
  REAL(KIND = 8) :: ka_mo3_var_222(19, 8)
  REAL(KIND = 8) :: kb_mco2_var_223(19, 8)
  REAL(KIND = 8) :: kb_mn2o_var_224(19, 8)
  REAL(KIND = 8) :: selfref_var_225(10, 8)
  REAL(KIND = 8) :: forref_var_226(4, 8)
END MODULE yoerrta8
MODULE yoerrta9
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: fracrefa_var_227(12, 9)
  REAL(KIND = 8), DIMENSION(12) :: fracrefb_var_228
  REAL(KIND = 8) :: absa_var_229(585, 12)
  REAL(KIND = 8) :: absb_var_230(235, 12)
  REAL(KIND = 8) :: ka_mn2o_var_231(9, 19, 12)
  REAL(KIND = 8) :: kb_mn2o_var_232(19, 12)
  REAL(KIND = 8) :: selfref_var_233(10, 12)
  REAL(KIND = 8) :: forref_var_234(4, 12)
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
  REAL(KIND = 8), DIMENSION(59) :: preflog_var_235
  REAL(KIND = 8), DIMENSION(59) :: tref_var_236
  REAL(KIND = 8) :: chi_mls(7, 59)
END MODULE yoerrtrf
MODULE yoerrtwn
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4), DIMENSION(16) :: nspa_var_237
  INTEGER(KIND = 4), DIMENSION(16) :: nspb_var_238
  REAL(KIND = 8), DIMENSION(16) :: delwave
  REAL(KIND = 8), DIMENSION(181, 16) :: totplnk
END MODULE yoerrtwn
MODULE yoesrta16
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_239, strrat1
  INTEGER(KIND = 4) :: layreffr_var_240
  REAL(KIND = 8) :: absa_var_241(585, 16)
  REAL(KIND = 8) :: absb_var_242(235, 16)
  REAL(KIND = 8) :: selfrefc_var_243(10, 16), forrefc_var_244(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_245(16)
END MODULE yoesrta16
MODULE yoesrta17
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_246, strrat_var_247
  INTEGER(KIND = 4) :: layreffr_var_248
  REAL(KIND = 8) :: absa_var_249(585, 16)
  REAL(KIND = 8) :: absb_var_250(1175, 16)
  REAL(KIND = 8) :: selfrefc_var_251(10, 16), forrefc_var_252(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_253(16, 5)
END MODULE yoesrta17
MODULE yoesrta18
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_254, strrat_var_255
  INTEGER(KIND = 4) :: layreffr_var_256
  REAL(KIND = 8) :: absa_var_257(585, 16)
  REAL(KIND = 8) :: absb_var_258(235, 16)
  REAL(KIND = 8) :: selfrefc_var_259(10, 16), forrefc_var_260(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_261(16, 9)
END MODULE yoesrta18
MODULE yoesrta19
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_262, strrat_var_263
  INTEGER(KIND = 4) :: layreffr_var_264
  REAL(KIND = 8) :: absa_var_265(585, 16)
  REAL(KIND = 8) :: absb_var_266(235, 16)
  REAL(KIND = 8) :: selfrefc_var_267(10, 16), forrefc_var_268(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_269(16, 9)
END MODULE yoesrta19
MODULE yoesrta20
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_270
  INTEGER(KIND = 4) :: layreffr_var_271
  REAL(KIND = 8) :: absa_var_272(65, 16)
  REAL(KIND = 8) :: absb_var_273(235, 16)
  REAL(KIND = 8) :: selfrefc_var_274(10, 16), forrefc_var_275(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_276(16), absch4c(16)
END MODULE yoesrta20
MODULE yoesrta21
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_277, strrat_var_278
  INTEGER(KIND = 4) :: layreffr_var_279
  REAL(KIND = 8) :: absa_var_280(585, 16)
  REAL(KIND = 8) :: absb_var_281(1175, 16)
  REAL(KIND = 8) :: selfrefc_var_282(10, 16), forrefc_var_283(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_284(16, 9)
END MODULE yoesrta21
MODULE yoesrta22
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_285, strrat_var_286
  INTEGER(KIND = 4) :: layreffr_var_287
  REAL(KIND = 8) :: absa_var_288(585, 16)
  REAL(KIND = 8) :: absb_var_289(235, 16)
  REAL(KIND = 8) :: selfrefc_var_290(10, 16), forrefc_var_291(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_292(16, 9)
END MODULE yoesrta22
MODULE yoesrta23
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: givfac
  INTEGER(KIND = 4) :: layreffr_var_293
  REAL(KIND = 8) :: absa_var_294(65, 16)
  REAL(KIND = 8) :: selfrefc_var_295(10, 16), forrefc_var_296(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_297(16), raylc_var_298(16)
END MODULE yoesrta23
MODULE yoesrta24
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: strrat_var_299
  INTEGER(KIND = 4) :: layreffr_var_300
  REAL(KIND = 8) :: absa_var_301(585, 16)
  REAL(KIND = 8) :: absb_var_302(235, 16)
  REAL(KIND = 8) :: selfrefc_var_303(10, 16), forrefc_var_304(3, 16)
  REAL(KIND = 8) :: sfluxrefc_var_305(16, 9)
  REAL(KIND = 8) :: abso3ac_var_306(16), abso3bc_var_307(16), raylac(16, 9), raylbc(16)
END MODULE yoesrta24
MODULE yoesrta25
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4) :: layreffr_var_308
  REAL(KIND = 8) :: absa_var_309(65, 16)
  REAL(KIND = 8) :: sfluxrefc_var_310(16)
  REAL(KIND = 8) :: raylc_var_311(16), abso3ac_var_312(16), abso3bc_var_313(16)
END MODULE yoesrta25
MODULE yoesrta26
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: sfluxrefc_var_314(16), raylc_var_315(16)
END MODULE yoesrta26
MODULE yoesrta27
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: scalekur
  INTEGER(KIND = 4) :: layreffr_var_316
  REAL(KIND = 8) :: absa_var_317(65, 16)
  REAL(KIND = 8) :: absb_var_318(235, 16)
  REAL(KIND = 8) :: sfluxrefc_var_319(16), raylc_var_320(16)
END MODULE yoesrta27
MODULE yoesrta28
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_321, strrat_var_322
  INTEGER(KIND = 4) :: layreffr_var_323
  REAL(KIND = 8) :: absa_var_324(585, 16)
  REAL(KIND = 8) :: absb_var_325(1175, 16)
  REAL(KIND = 8) :: sfluxrefc_var_326(16, 5)
END MODULE yoesrta28
MODULE yoesrta29
  IMPLICIT NONE
  SAVE
  REAL(KIND = 8) :: rayl_var_327
  INTEGER(KIND = 4) :: layreffr_var_328
  REAL(KIND = 8) :: absa_var_329(65, 16)
  REAL(KIND = 8) :: absb_var_330(235, 16)
  REAL(KIND = 8) :: selfrefc_var_331(10, 16), forrefc_var_332(4, 16)
  REAL(KIND = 8) :: sfluxrefc_var_333(16), absh2oc(16), absco2c(16)
END MODULE yoesrta29
MODULE yoesrtwn
  IMPLICIT NONE
  SAVE
  INTEGER(KIND = 4), DIMENSION(16 : 29) :: nspa_var_334
  INTEGER(KIND = 4), DIMENSION(16 : 29) :: nspb_var_335
  REAL(KIND = 8), DIMENSION(59) :: preflog_var_336
  REAL(KIND = 8), DIMENSION(59) :: tref_var_337
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
  RECURSIVE SUBROUTINE assert_units_gas(this_var_339, iunits, igas, scale_factor, istatus)
    CLASS(gas_type), INTENT(IN) :: this_var_339
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
  ELEMENTAL SUBROUTINE sample_from_pdf(this_var_341, fsd_var_342, cdf_var_343, val_var_344)
    CLASS(pdf_sampler_type), INTENT(IN) :: this_var_341
    REAL(KIND = 8), INTENT(IN) :: fsd_var_342, cdf_var_343
    REAL(KIND = 8), INTENT(OUT) :: val_var_344
    INTEGER :: ifsd_var_345, icdf_var_346
    REAL(KIND = 8) :: wfsd_var_347, wcdf_var_348
    wcdf_var_348 = cdf_var_343 * (this_var_341 % ncdf - 1) + 1.0D0
    icdf_var_346 = MAX(1, MIN(INT(wcdf_var_348), this_var_341 % ncdf - 1))
    wcdf_var_348 = MAX(0.0D0, MIN(wcdf_var_348 - icdf_var_346, 1.0D0))
    wfsd_var_347 = (fsd_var_342 - this_var_341 % fsd1) * this_var_341 % inv_fsd_interval + 1.0D0
    ifsd_var_345 = MAX(1, MIN(INT(wfsd_var_347), this_var_341 % nfsd - 1))
    wfsd_var_347 = MAX(0.0D0, MIN(wfsd_var_347 - ifsd_var_345, 1.0D0))
    val_var_344 = (1.0D0 - wcdf_var_348) * (1.0D0 - wfsd_var_347) * this_var_341 % val(icdf_var_346, ifsd_var_345) + (1.0D0 - wcdf_var_348) * wfsd_var_347 * this_var_341 % val(icdf_var_346, ifsd_var_345 + 1) + wcdf_var_348 * (1.0D0 - wfsd_var_347) * this_var_341 % val(icdf_var_346 + 1, ifsd_var_345) + wcdf_var_348 * wfsd_var_347 * this_var_341 % val(icdf_var_346 + 1, ifsd_var_345 + 1)
  END SUBROUTINE sample_from_pdf
  SUBROUTINE sample_from_pdf_masked_block(this_var_349, nz, ng_var_350, fsd_var_351, cdf_var_352, val_var_353, mask_var_354)
    CLASS(pdf_sampler_type), INTENT(IN) :: this_var_349
    INTEGER, INTENT(IN) :: nz, ng_var_350
    REAL(KIND = 8), INTENT(IN) :: fsd_var_351(nz), cdf_var_352(ng_var_350, nz)
    REAL(KIND = 8), INTENT(OUT) :: val_var_353(:, :)
    LOGICAL, INTENT(IN), OPTIONAL :: mask_var_354(nz)
    INTEGER :: jz, jg_var_355
    INTEGER :: ifsd_var_356, icdf_var_357
    REAL(KIND = 8) :: wfsd_var_358, wcdf_var_359
    DO jz = 1, nz
      IF (mask_var_354(jz)) THEN
        DO jg_var_355 = 1, ng_var_350
          IF (cdf_var_352(jg_var_355, jz) > 0.0D0) THEN
            wcdf_var_359 = cdf_var_352(jg_var_355, jz) * (this_var_349 % ncdf - 1) + 1.0D0
            icdf_var_357 = MAX(1, MIN(INT(wcdf_var_359), this_var_349 % ncdf - 1))
            wcdf_var_359 = MAX(0.0D0, MIN(wcdf_var_359 - icdf_var_357, 1.0D0))
            wfsd_var_358 = (fsd_var_351(jz) - this_var_349 % fsd1) * this_var_349 % inv_fsd_interval + 1.0D0
            ifsd_var_356 = MAX(1, MIN(INT(wfsd_var_358), this_var_349 % nfsd - 1))
            wfsd_var_358 = MAX(0.0D0, MIN(wfsd_var_358 - ifsd_var_356, 1.0D0))
            val_var_353(jg_var_355, jz) = (1.0D0 - wcdf_var_359) * (1.0D0 - wfsd_var_358) * this_var_349 % val(icdf_var_357, ifsd_var_356) + (1.0D0 - wcdf_var_359) * wfsd_var_358 * this_var_349 % val(icdf_var_357, ifsd_var_356 + 1) + wcdf_var_359 * (1.0D0 - wfsd_var_358) * this_var_349 % val(icdf_var_357 + 1, ifsd_var_356) + wcdf_var_359 * wfsd_var_358 * this_var_349 % val(icdf_var_357 + 1, ifsd_var_356 + 1)
          ELSE
            val_var_353(jg_var_355, jz) = 0.0D0
          END IF
        END DO
      END IF
    END DO
  END SUBROUTINE sample_from_pdf_masked_block
END MODULE radiation_pdf_sampler
MODULE radiation_cloud_generator
  CONTAINS
  SUBROUTINE cloud_generator(ng_var_360, nlev_var_361, i_overlap_scheme, iseed_var_362, frac_threshold_var_363, frac_var_364, overlap_param_var_365, decorrelation_scaling, fractional_std_var_366, pdf_sampler_var_367, od_scaling_var_368, total_cloud_cover_var_369, use_beta_overlap, use_vectorizable_generator)
    USE radiation_pdf_sampler, ONLY: pdf_sampler_type
    USE random_numbers_mix, ONLY: initialize_random_numbers, randomnumberstream, uniform_distribution
    USE radiation_cloud_cover, ONLY: cum_cloud_cover_exp_ran
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng_var_360
    INTEGER, INTENT(IN) :: nlev_var_361
    INTEGER, INTENT(IN) :: i_overlap_scheme
    INTEGER, INTENT(IN) :: iseed_var_362
    REAL(KIND = 8), INTENT(IN) :: frac_threshold_var_363
    REAL(KIND = 8), INTENT(IN) :: frac_var_364(nlev_var_361)
    REAL(KIND = 8), INTENT(IN) :: overlap_param_var_365(nlev_var_361 - 1)
    REAL(KIND = 8), INTENT(IN) :: decorrelation_scaling
    REAL(KIND = 8), INTENT(IN) :: fractional_std_var_366(nlev_var_361)
    TYPE(pdf_sampler_type), INTENT(IN) :: pdf_sampler_var_367
    LOGICAL, INTENT(IN), OPTIONAL :: use_beta_overlap
    LOGICAL, INTENT(IN), OPTIONAL :: use_vectorizable_generator
    REAL(KIND = 8), INTENT(OUT) :: od_scaling_var_368(ng_var_360, nlev_var_361)
    REAL(KIND = 8), INTENT(OUT) :: total_cloud_cover_var_369
    REAL(KIND = 8) :: cum_cloud_cover_var_370(nlev_var_361)
    REAL(KIND = 8) :: trigger_var_371
    REAL(KIND = 8) :: rand_top(ng_var_360)
    REAL(KIND = 8) :: overlap_param_inhom_var_372(nlev_var_361 - 1)
    TYPE(randomnumberstream) :: random_stream_var_373
    INTEGER :: ibegin_var_374, iend_var_375
    INTEGER :: itrigger_var_376
    INTEGER :: jlev_var_377, jg_var_378
    REAL(KIND = 8), DIMENSION(nlev_var_361 - 1) :: pair_cloud_cover_var_379, overhang_var_380
    LOGICAL :: use_vec_gen
    CALL cum_cloud_cover_exp_ran(nlev_var_361, frac_var_364, overlap_param_var_365, cum_cloud_cover_var_370, pair_cloud_cover_var_379, .FALSE.)
    total_cloud_cover_var_369 = cum_cloud_cover_var_370(nlev_var_361)
    DO jlev_var_377 = 1, nlev_var_361 - 1
      overhang_var_380(jlev_var_377) = cum_cloud_cover_var_370(jlev_var_377 + 1) - cum_cloud_cover_var_370(jlev_var_377)
    END DO
    IF (total_cloud_cover_var_369 < 1D-06) THEN
      total_cloud_cover_var_369 = 0.0D0
    ELSE
      jlev_var_377 = 1
      DO WHILE (frac_var_364(jlev_var_377) <= 0.0D0)
        jlev_var_377 = jlev_var_377 + 1
      END DO
      ibegin_var_374 = jlev_var_377
      iend_var_375 = jlev_var_377
      DO jlev_var_377 = jlev_var_377 + 1, nlev_var_361
        IF (frac_var_364(jlev_var_377) > 0.0D0) THEN
          iend_var_375 = jlev_var_377
        END IF
      END DO
      overlap_param_inhom_var_372 = overlap_param_var_365
      DO jlev_var_377 = ibegin_var_374, iend_var_375 - 1
        IF (overlap_param_var_365(jlev_var_377) > 0.0D0) THEN
          overlap_param_inhom_var_372(jlev_var_377) = overlap_param_var_365(jlev_var_377) ** (2.0D0)
        END IF
      END DO
      od_scaling_var_368 = 0.0D0
      use_vec_gen = .FALSE.
      IF (.NOT. use_vec_gen) THEN
        CALL initialize_random_numbers(iseed_var_362, random_stream_var_373)
        CALL uniform_distribution(rand_top, random_stream_var_373)
        DO jg_var_378 = 1, ng_var_360
          trigger_var_371 = rand_top(jg_var_378) * total_cloud_cover_var_369
          jlev_var_377 = ibegin_var_374
          DO WHILE (trigger_var_371 > cum_cloud_cover_var_370(jlev_var_377) .AND. jlev_var_377 < iend_var_375)
            jlev_var_377 = jlev_var_377 + 1
          END DO
          itrigger_var_376 = jlev_var_377
          CALL generate_column_exp_ran(ng_var_360, nlev_var_361, jg_var_378, random_stream_var_373, pdf_sampler_var_367, frac_var_364, pair_cloud_cover_var_379, cum_cloud_cover_var_370, overhang_var_380, fractional_std_var_366, overlap_param_inhom_var_372, itrigger_var_376, iend_var_375, od_scaling_var_368)
        END DO
      ELSE
        CALL generate_columns_exp_ran(ng_var_360, nlev_var_361, iseed_var_362, pdf_sampler_var_367, total_cloud_cover_var_369, 1D-06, frac_var_364, pair_cloud_cover_var_379, cum_cloud_cover_var_370, overhang_var_380, fractional_std_var_366, overlap_param_inhom_var_372, ibegin_var_374, iend_var_375, od_scaling_var_368)
      END IF
    END IF
  END SUBROUTINE cloud_generator
  SUBROUTINE generate_column_exp_ran(ng_var_381, nlev_var_383, ig_var_382, random_stream_var_384, pdf_sampler_var_385, frac_var_386, pair_cloud_cover_var_389, cum_cloud_cover_var_387, overhang_var_390, fractional_std_var_388, overlap_param_inhom_var_391, itrigger_var_392, iend_var_393, od_scaling_var_394)
    USE random_numbers_mix, ONLY: randomnumberstream, uniform_distribution
    USE radiation_pdf_sampler, ONLY: pdf_sampler_type, sample_from_pdf
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng_var_381, ig_var_382
    INTEGER, INTENT(IN) :: nlev_var_383
    TYPE(randomnumberstream), INTENT(INOUT) :: random_stream_var_384
    TYPE(pdf_sampler_type), INTENT(IN) :: pdf_sampler_var_385
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_383) :: frac_var_386, cum_cloud_cover_var_387, fractional_std_var_388
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_383 - 1) :: pair_cloud_cover_var_389, overhang_var_390
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_383 - 1) :: overlap_param_inhom_var_391
    INTEGER, INTENT(IN) :: itrigger_var_392, iend_var_393
    REAL(KIND = 8), INTENT(INOUT), DIMENSION(ng_var_381, nlev_var_383) :: od_scaling_var_394
    INTEGER :: jlev_var_395, jcloud
    INTEGER :: n_layers_to_scale
    INTEGER :: iy
    LOGICAL :: do_fill_od_scaling
    REAL(KIND = 8) :: rand_cloud_var_396(nlev_var_383)
    REAL(KIND = 8) :: rand_inhom1(nlev_var_383), rand_inhom2_var_397(nlev_var_383)
    n_layers_to_scale = 1
    iy = 0
    CALL uniform_distribution(rand_cloud_var_396(1 : (iend_var_393 + 1 - itrigger_var_392)), random_stream_var_384)
    DO jlev_var_395 = itrigger_var_392 + 1, iend_var_393 + 1
      do_fill_od_scaling = .FALSE.
      IF (jlev_var_395 <= iend_var_393) THEN
        iy = iy + 1
        IF (n_layers_to_scale > 0) THEN
          IF (rand_cloud_var_396(iy) * frac_var_386(jlev_var_395 - 1) < frac_var_386(jlev_var_395) + frac_var_386(jlev_var_395 - 1) - pair_cloud_cover_var_389(jlev_var_395 - 1)) THEN
            n_layers_to_scale = n_layers_to_scale + 1
          ELSE
            do_fill_od_scaling = .TRUE.
          END IF
        ELSE
          IF (rand_cloud_var_396(iy) * (cum_cloud_cover_var_387(jlev_var_395 - 1) - frac_var_386(jlev_var_395 - 1)) < pair_cloud_cover_var_389(jlev_var_395 - 1) - overhang_var_390(jlev_var_395 - 1) - frac_var_386(jlev_var_395 - 1)) THEN
            n_layers_to_scale = 1
          END IF
        END IF
      ELSE
        do_fill_od_scaling = .TRUE.
      END IF
      IF (do_fill_od_scaling) THEN
        CALL uniform_distribution(rand_inhom1(1 : n_layers_to_scale), random_stream_var_384)
        CALL uniform_distribution(rand_inhom2_var_397(1 : n_layers_to_scale), random_stream_var_384)
        DO jcloud = 2, n_layers_to_scale
          IF (rand_inhom2_var_397(jcloud) < overlap_param_inhom_var_391(jlev_var_395 - n_layers_to_scale + jcloud - 2)) THEN
            rand_inhom1(jcloud) = rand_inhom1(jcloud - 1)
          END IF
        END DO
        CALL sample_from_pdf(pdf_sampler_var_385, fractional_std_var_388(jlev_var_395 - n_layers_to_scale : jlev_var_395 - 1), rand_inhom1(1 : n_layers_to_scale), od_scaling_var_394(ig_var_382, jlev_var_395 - n_layers_to_scale : jlev_var_395 - 1))
        n_layers_to_scale = 0
      END IF
    END DO
  END SUBROUTINE generate_column_exp_ran
  SUBROUTINE generate_columns_exp_ran(ng_var_398, nlev_var_399, iseed_var_400, pdf_sampler_var_401, total_cloud_cover_var_402, frac_threshold_var_403, frac_var_404, pair_cloud_cover_var_407, cum_cloud_cover_var_405, overhang_var_408, fractional_std_var_406, overlap_param_inhom_var_409, ibegin_var_410, iend_var_411, od_scaling_var_412)
    USE radiation_random_numbers, ONLY: initialize, rng_type, uniform_distribution_1d, uniform_distribution_2d, uniform_distribution_2d_masked
    USE radiation_pdf_sampler, ONLY: pdf_sampler_type, sample_from_pdf_masked_block
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng_var_398
    INTEGER, INTENT(IN) :: nlev_var_399
    INTEGER, INTENT(IN) :: iseed_var_400
    TYPE(rng_type) :: random_number_generator
    TYPE(pdf_sampler_type), INTENT(IN) :: pdf_sampler_var_401
    REAL(KIND = 8), INTENT(IN) :: total_cloud_cover_var_402
    REAL(KIND = 8), INTENT(IN) :: frac_threshold_var_403
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_399) :: frac_var_404, cum_cloud_cover_var_405, fractional_std_var_406
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_399 - 1) :: pair_cloud_cover_var_407, overhang_var_408
    REAL(KIND = 8), INTENT(IN), DIMENSION(nlev_var_399 - 1) :: overlap_param_inhom_var_409
    INTEGER, INTENT(INOUT) :: ibegin_var_410, iend_var_411
    REAL(KIND = 8), INTENT(INOUT), DIMENSION(ng_var_398, nlev_var_399) :: od_scaling_var_412
    INTEGER :: jlev_var_413, jg_var_414
    REAL(KIND = 8) :: rand_cloud_var_415(ng_var_398, ibegin_var_410 : iend_var_411)
    REAL(KIND = 8) :: rand_inhom(ng_var_398, ibegin_var_410 - 1 : iend_var_411), rand_inhom2_var_416(ng_var_398, ibegin_var_410 : iend_var_411)
    LOGICAL :: is_any_cloud(ibegin_var_410 : iend_var_411)
    REAL(KIND = 8) :: trigger_var_417(ng_var_398)
    LOGICAL :: is_cloud(ng_var_398)
    LOGICAL :: prev_cloud(ng_var_398)
    LOGICAL :: first_cloud(ng_var_398)
    LOGICAL :: found_cloud(ng_var_398)
    is_any_cloud = (frac_var_404(ibegin_var_410 : iend_var_411) >= frac_threshold_var_403)
    CALL initialize(random_number_generator, 1, iseed_var_90 = iseed_var_400, nmaxstreams_var_91 = ng_var_398)
    CALL uniform_distribution_1d(random_number_generator, trigger_var_417)
    CALL uniform_distribution_2d_masked(random_number_generator, rand_cloud_var_415, is_any_cloud)
    CALL uniform_distribution_2d(random_number_generator, rand_inhom)
    CALL uniform_distribution_2d_masked(random_number_generator, rand_inhom2_var_416, is_any_cloud)
    trigger_var_417 = trigger_var_417 * total_cloud_cover_var_402
    found_cloud = .FALSE.
    is_cloud = .FALSE.
    first_cloud = .FALSE.
    DO jlev_var_413 = ibegin_var_410, iend_var_411
      IF (is_any_cloud(jlev_var_413)) THEN
        DO jg_var_414 = 1, ng_var_398
          prev_cloud(jg_var_414) = is_cloud(jg_var_414)
          first_cloud(jg_var_414) = (trigger_var_417(jg_var_414) <= cum_cloud_cover_var_405(jlev_var_413) .AND. .NOT. found_cloud(jg_var_414))
          found_cloud(jg_var_414) = found_cloud(jg_var_414) .OR. first_cloud(jg_var_414)
          is_cloud(jg_var_414) = first_cloud(jg_var_414) .OR. found_cloud(jg_var_414) .AND. MERGE(rand_cloud_var_415(jg_var_414, jlev_var_413) * frac_var_404(jlev_var_413 - 1) < frac_var_404(jlev_var_413) + frac_var_404(jlev_var_413 - 1) - pair_cloud_cover_var_407(jlev_var_413 - 1), rand_cloud_var_415(jg_var_414, jlev_var_413) * (cum_cloud_cover_var_405(jlev_var_413 - 1) - frac_var_404(jlev_var_413 - 1)) < pair_cloud_cover_var_407(jlev_var_413 - 1) - overhang_var_408(jlev_var_413 - 1) - frac_var_404(jlev_var_413 - 1), prev_cloud(jg_var_414))
          rand_inhom(jg_var_414, jlev_var_413) = MERGE(MERGE(rand_inhom(jg_var_414, jlev_var_413 - 1), rand_inhom(jg_var_414, jlev_var_413), rand_inhom2_var_416(jg_var_414, jlev_var_413) < overlap_param_inhom_var_409(jlev_var_413 - 1) .AND. prev_cloud(jg_var_414)), 0.0D0, is_cloud(jg_var_414))
        END DO
      ELSE
        is_cloud = .FALSE.
      END IF
    END DO
    CALL sample_from_pdf_masked_block(pdf_sampler_var_401, iend_var_411 - ibegin_var_410 + 1, ng_var_398, fractional_std_var_406(ibegin_var_410 : iend_var_411), rand_inhom(:, ibegin_var_410 : iend_var_411), od_scaling_var_412(:, ibegin_var_410 : iend_var_411), is_any_cloud)
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
  SUBROUTINE calc_surface_spectral(this_var_424, config_var_425, istartcol_var_426, iendcol_var_427)
    USE radiation_config, ONLY: config_type
    CLASS(flux_type), INTENT(INOUT) :: this_var_424
    TYPE(config_type), INTENT(IN) :: config_var_425
    INTEGER, INTENT(IN) :: istartcol_var_426, iendcol_var_427
    INTEGER :: jcol_var_428
    DO jcol_var_428 = istartcol_var_426, iendcol_var_427
      CALL indexed_sum(this_var_424 % sw_dn_direct_surf_g(:, jcol_var_428), config_var_425 % i_band_from_reordered_g_sw, this_var_424 % sw_dn_direct_surf_band(:, jcol_var_428))
      CALL indexed_sum(this_var_424 % sw_dn_diffuse_surf_g(:, jcol_var_428), config_var_425 % i_band_from_reordered_g_sw, this_var_424 % sw_dn_surf_band(:, jcol_var_428))
      this_var_424 % sw_dn_surf_band(:, jcol_var_428) = this_var_424 % sw_dn_surf_band(:, jcol_var_428) + this_var_424 % sw_dn_direct_surf_band(:, jcol_var_428)
    END DO
    DO jcol_var_428 = istartcol_var_426, iendcol_var_427
      CALL indexed_sum(this_var_424 % sw_dn_direct_surf_clear_g(:, jcol_var_428), config_var_425 % i_band_from_reordered_g_sw, this_var_424 % sw_dn_direct_surf_clear_band(:, jcol_var_428))
      CALL indexed_sum(this_var_424 % sw_dn_diffuse_surf_clear_g(:, jcol_var_428), config_var_425 % i_band_from_reordered_g_sw, this_var_424 % sw_dn_surf_clear_band(:, jcol_var_428))
      this_var_424 % sw_dn_surf_clear_band(:, jcol_var_428) = this_var_424 % sw_dn_surf_clear_band(:, jcol_var_428) + this_var_424 % sw_dn_direct_surf_clear_band(:, jcol_var_428)
    END DO
  END SUBROUTINE calc_surface_spectral
  SUBROUTINE calc_toa_spectral(this_var_429, config_var_430, istartcol_var_431, iendcol_var_432)
    USE radiation_config, ONLY: config_type
    CLASS(flux_type), INTENT(INOUT) :: this_var_429
    TYPE(config_type), INTENT(IN) :: config_var_430
    INTEGER, INTENT(IN) :: istartcol_var_431, iendcol_var_432
  END SUBROUTINE calc_toa_spectral
  PURE SUBROUTINE indexed_sum(source_var_433, ind_var_434, dest)
    REAL(KIND = 8), INTENT(IN) :: source_var_433(:)
    INTEGER, INTENT(IN) :: ind_var_434(:)
    REAL(KIND = 8), INTENT(OUT) :: dest(:)
    INTEGER :: ig_var_435, jg_var_436, istart, iend_var_437
    dest = 0.0
    istart = LBOUND(source_var_433, 1)
    iend_var_437 = UBOUND(source_var_433, 1)
    DO jg_var_436 = istart, iend_var_437
      ig_var_435 = ind_var_434(jg_var_436)
      dest(ig_var_435) = dest(ig_var_435) + source_var_433(jg_var_436)
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
  SUBROUTINE get_albedos(this_var_441, istartcol_var_443, iendcol_var_444, config_var_442, sw_albedo_direct_var_446, sw_albedo_diffuse_var_447, lw_albedo_var_445)
    USE radiation_config, ONLY: config_type
    USE yomlun_ifsaux, ONLY: nulerr
    USE radiation_io, ONLY: radiation_abort
    CLASS(single_level_type), INTENT(IN) :: this_var_441
    TYPE(config_type), INTENT(IN) :: config_var_442
    INTEGER, INTENT(IN) :: istartcol_var_443, iendcol_var_444
    REAL(KIND = 8), INTENT(OUT), OPTIONAL :: lw_albedo_var_445(140, istartcol_var_443 : iendcol_var_444)
    REAL(KIND = 8), INTENT(OUT), DIMENSION(112, istartcol_var_443 : iendcol_var_444) :: sw_albedo_direct_var_446, sw_albedo_diffuse_var_447
    REAL(KIND = 8) :: sw_albedo_band(istartcol_var_443 : iendcol_var_444, 14)
    INTEGER :: nalbedoband
    INTEGER :: jband_var_448, jalbedoband, jcol_var_449
    nalbedoband = SIZE(config_var_442 % sw_albedo_weights, 1)
    IF (SIZE(this_var_441 % sw_albedo, 2) /= nalbedoband) THEN
      WRITE(nulerr, '(a,i0,a)') '*** error: single_level%sw_albedo does not have the expected ', nalbedoband, ' bands'
      CALL radiation_abort
    END IF
    DO jband_var_448 = 1, 14
      DO jcol_var_449 = istartcol_var_443, iendcol_var_444
        sw_albedo_band(jcol_var_449, jband_var_448) = 0.0D0
      END DO
    END DO
    DO jband_var_448 = 1, 14
      DO jalbedoband = 1, nalbedoband
        IF (config_var_442 % sw_albedo_weights(jalbedoband, jband_var_448) /= 0.0D0) THEN
          DO jcol_var_449 = istartcol_var_443, iendcol_var_444
            sw_albedo_band(jcol_var_449, jband_var_448) = sw_albedo_band(jcol_var_449, jband_var_448) + config_var_442 % sw_albedo_weights(jalbedoband, jband_var_448) * this_var_441 % sw_albedo(jcol_var_449, jalbedoband)
          END DO
        END IF
      END DO
    END DO
    sw_albedo_diffuse_var_447 = TRANSPOSE(sw_albedo_band(istartcol_var_443 : iendcol_var_444, config_var_442 % i_band_from_reordered_g_sw))
    DO jband_var_448 = 1, 14
      DO jcol_var_449 = istartcol_var_443, iendcol_var_444
        sw_albedo_band(jcol_var_449, jband_var_448) = 0.0D0
      END DO
    END DO
    DO jband_var_448 = 1, 14
      DO jalbedoband = 1, nalbedoband
        IF (config_var_442 % sw_albedo_weights(jalbedoband, jband_var_448) /= 0.0D0) THEN
          DO jcol_var_449 = istartcol_var_443, iendcol_var_444
            sw_albedo_band(jcol_var_449, jband_var_448) = sw_albedo_band(jcol_var_449, jband_var_448) + config_var_442 % sw_albedo_weights(jalbedoband, jband_var_448) * this_var_441 % sw_albedo_direct(jcol_var_449, jalbedoband)
          END DO
        END IF
      END DO
    END DO
    sw_albedo_direct_var_446 = TRANSPOSE(sw_albedo_band(istartcol_var_443 : iendcol_var_444, config_var_442 % i_band_from_reordered_g_sw))
    IF (MAXVAL(config_var_442 % i_emiss_from_band_lw) > SIZE(this_var_441 % lw_emissivity, 2)) THEN
      WRITE(nulerr, '(a,i0,a)') '*** error: single_level%lw_emissivity has fewer than required ', MAXVAL(config_var_442 % i_emiss_from_band_lw), ' bands'
      CALL radiation_abort
    END IF
    lw_albedo_var_445 = 1.0D0 - TRANSPOSE(this_var_441 % lw_emissivity(istartcol_var_443 : iendcol_var_444, config_var_442 % i_emiss_from_band_lw(config_var_442 % i_band_from_reordered_g_lw)))
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
  ELEMENTAL SUBROUTINE delta_eddington_extensive(od_var_452, scat_od_var_453, scat_od_g)
    REAL(KIND = 8), INTENT(INOUT) :: od_var_452, scat_od_var_453, scat_od_g
    REAL(KIND = 8) :: f_var_454, g_var_455
    IF (scat_od_var_453 > 0.0D0) THEN
      g_var_455 = scat_od_g / scat_od_var_453
    ELSE
      g_var_455 = 0.0
    END IF
    f_var_454 = g_var_455 * g_var_455
    od_var_452 = od_var_452 - scat_od_var_453 * f_var_454
    scat_od_var_453 = scat_od_var_453 * (1.0D0 - f_var_454)
    scat_od_g = scat_od_var_453 * g_var_455 / (1.0D0 + g_var_455)
  END SUBROUTINE delta_eddington_extensive
  SUBROUTINE add_aerosol_optics(nlev_var_456, istartcol_var_457, iendcol_var_458, config_var_459, thermodynamics_var_460, gas_var_461, aerosol_var_462, od_lw_var_463, ssa_lw_var_464, g_lw_var_465, od_sw_var_466, ssa_sw_var_467, g_sw_var_468)
    USE radiation_config, ONLY: config_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE radiation_gas, ONLY: gas_type
    USE radiation_aerosol, ONLY: aerosol_type
    INTEGER, INTENT(IN) :: nlev_var_456
    INTEGER, INTENT(IN) :: istartcol_var_457, iendcol_var_458
    TYPE(config_type), INTENT(IN), TARGET :: config_var_459
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics_var_460
    TYPE(gas_type), INTENT(IN) :: gas_var_461
    TYPE(aerosol_type), INTENT(IN) :: aerosol_var_462
    REAL(KIND = 8), DIMENSION(140, nlev_var_456, istartcol_var_457 : iendcol_var_458), INTENT(INOUT) :: od_lw_var_463
    REAL(KIND = 8), DIMENSION(0, nlev_var_456, istartcol_var_457 : iendcol_var_458), INTENT(OUT) :: ssa_lw_var_464, g_lw_var_465
    REAL(KIND = 8), DIMENSION(112, nlev_var_456, istartcol_var_457 : iendcol_var_458), INTENT(INOUT) :: od_sw_var_466, ssa_sw_var_467
    REAL(KIND = 8), DIMENSION(112, nlev_var_456, istartcol_var_457 : iendcol_var_458), INTENT(OUT) :: g_sw_var_468
    CALL add_aerosol_optics_direct(nlev_var_456, istartcol_var_457, iendcol_var_458, config_var_459, aerosol_var_462, od_lw_var_463, ssa_lw_var_464, g_lw_var_465, od_sw_var_466, ssa_sw_var_467, g_sw_var_468)
  END SUBROUTINE add_aerosol_optics
  SUBROUTINE add_aerosol_optics_direct(nlev_var_469, istartcol_var_470, iendcol_var_471, config_var_472, aerosol_var_473, od_lw_var_474, ssa_lw_var_475, g_lw_var_476, od_sw_var_477, ssa_sw_var_478, g_sw_var_479)
    USE radiation_config, ONLY: config_type
    USE radiation_aerosol, ONLY: aerosol_type
    USE yomlun_ifsaux, ONLY: nulerr
    USE radiation_io, ONLY: radiation_abort
    INTEGER, INTENT(IN) :: nlev_var_469
    INTEGER, INTENT(IN) :: istartcol_var_470, iendcol_var_471
    TYPE(config_type), INTENT(IN), TARGET :: config_var_472
    TYPE(aerosol_type), INTENT(IN) :: aerosol_var_473
    REAL(KIND = 8), DIMENSION(140, nlev_var_469, istartcol_var_470 : iendcol_var_471), INTENT(INOUT) :: od_lw_var_474
    REAL(KIND = 8), DIMENSION(0, nlev_var_469, istartcol_var_470 : iendcol_var_471), INTENT(OUT) :: ssa_lw_var_475, g_lw_var_476
    REAL(KIND = 8), DIMENSION(112, nlev_var_469, istartcol_var_470 : iendcol_var_471), INTENT(INOUT) :: od_sw_var_477, ssa_sw_var_478
    REAL(KIND = 8), DIMENSION(112, nlev_var_469, istartcol_var_470 : iendcol_var_471), INTENT(OUT) :: g_sw_var_479
    REAL(KIND = 8) :: local_od, local_scat
    REAL(KIND = 8), DIMENSION(14, nlev_var_469) :: od_sw_aerosol, scat_sw_aerosol, scat_g_sw_aerosol
    REAL(KIND = 8), DIMENSION(16, nlev_var_469) :: od_lw_aerosol
    INTEGER :: jcol_var_480, jlev_var_481, jg_var_482, jb_var_483
    INTEGER :: istartlev_var_484, iendlev_var_485
    INTEGER :: iband_var_486
    IF (UBOUND(aerosol_var_473 % od_sw, 1) /= 14) THEN
      WRITE(nulerr, '(a,i0,a,i0)') '*** error: aerosol%od_sw contains ', UBOUND(aerosol_var_473 % od_sw, 1), ' band, expected ', 14
      CALL radiation_abort
    END IF
    istartlev_var_484 = LBOUND(aerosol_var_473 % od_sw, 2)
    iendlev_var_485 = UBOUND(aerosol_var_473 % od_sw, 2)
    DO jcol_var_480 = istartcol_var_470, iendcol_var_471
      DO jlev_var_481 = 1, nlev_var_469
        DO jg_var_482 = 1, 112
          g_sw_var_479(jg_var_482, jlev_var_481, jcol_var_480) = 0.0D0
        END DO
      END DO
    END DO
    DO jcol_var_480 = istartcol_var_470, iendcol_var_471
      DO jlev_var_481 = istartlev_var_484, iendlev_var_485
        DO jb_var_483 = 1, 14
          od_sw_aerosol(jb_var_483, jlev_var_481) = aerosol_var_473 % od_sw(jb_var_483, jlev_var_481, jcol_var_480)
          scat_sw_aerosol(jb_var_483, jlev_var_481) = aerosol_var_473 % ssa_sw(jb_var_483, jlev_var_481, jcol_var_480) * od_sw_aerosol(jb_var_483, jlev_var_481)
          scat_g_sw_aerosol(jb_var_483, jlev_var_481) = aerosol_var_473 % g_sw(jb_var_483, jlev_var_481, jcol_var_480) * scat_sw_aerosol(jb_var_483, jlev_var_481)
          CALL delta_eddington_extensive(od_sw_aerosol(jb_var_483, jlev_var_481), scat_sw_aerosol(jb_var_483, jlev_var_481), scat_g_sw_aerosol(jb_var_483, jlev_var_481))
        END DO
      END DO
      DO jlev_var_481 = istartlev_var_484, iendlev_var_485
        IF (od_sw_aerosol(1, jlev_var_481) > 0.0D0) THEN
          DO jg_var_482 = 1, 112
            iband_var_486 = config_var_472 % i_band_from_reordered_g_sw(jg_var_482)
            local_od = od_sw_var_477(jg_var_482, jlev_var_481, jcol_var_480) + od_sw_aerosol(iband_var_486, jlev_var_481)
            local_scat = ssa_sw_var_478(jg_var_482, jlev_var_481, jcol_var_480) * od_sw_var_477(jg_var_482, jlev_var_481, jcol_var_480) + scat_sw_aerosol(iband_var_486, jlev_var_481)
            g_sw_var_479(jg_var_482, jlev_var_481, jcol_var_480) = scat_g_sw_aerosol(iband_var_486, jlev_var_481) / local_scat
            local_od = od_sw_var_477(jg_var_482, jlev_var_481, jcol_var_480) + od_sw_aerosol(iband_var_486, jlev_var_481)
            ssa_sw_var_478(jg_var_482, jlev_var_481, jcol_var_480) = local_scat / local_od
            od_sw_var_477(jg_var_482, jlev_var_481, jcol_var_480) = local_od
          END DO
        END IF
      END DO
    END DO
    IF (UBOUND(aerosol_var_473 % od_lw, 1) /= 16) THEN
      WRITE(nulerr, '(a,i0,a,i0)') '*** error: aerosol%od_lw contains ', UBOUND(aerosol_var_473 % od_lw, 1), ' band, expected ', 16
      CALL radiation_abort
    END IF
    istartlev_var_484 = LBOUND(aerosol_var_473 % od_lw, 2)
    iendlev_var_485 = UBOUND(aerosol_var_473 % od_lw, 2)
    DO jcol_var_480 = istartcol_var_470, iendcol_var_471
      DO jlev_var_481 = istartlev_var_484, iendlev_var_485
        DO jb_var_483 = 1, 16
          od_lw_aerosol(jb_var_483, jlev_var_481) = aerosol_var_473 % od_lw(jb_var_483, jlev_var_481, jcol_var_480) * (1.0D0 - aerosol_var_473 % ssa_lw(jb_var_483, jlev_var_481, jcol_var_480))
        END DO
      END DO
      DO jlev_var_481 = istartlev_var_484, iendlev_var_485
        DO jg_var_482 = 1, 140
          od_lw_var_474(jg_var_482, jlev_var_481, jcol_var_480) = od_lw_var_474(jg_var_482, jlev_var_481, jcol_var_480) + od_lw_aerosol(config_var_472 % i_band_from_reordered_g_lw(jg_var_482), jlev_var_481)
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
  SUBROUTINE crop_cloud_fraction(this_var_490, istartcol_var_491, iendcol_var_492, cloud_fraction_threshold, cloud_mixing_ratio_threshold)
    CLASS(cloud_type), INTENT(INOUT) :: this_var_490
    INTEGER, INTENT(IN) :: istartcol_var_491, iendcol_var_492
    INTEGER :: nlev_var_493
    INTEGER :: jcol_var_494, jlev_var_495, jh
    REAL(KIND = 8) :: cloud_fraction_threshold, cloud_mixing_ratio_threshold
    REAL(KIND = 8) :: sum_mixing_ratio(istartcol_var_491 : iendcol_var_492)
    nlev_var_493 = SIZE(this_var_490 % fraction, 2)
    DO jlev_var_495 = 1, nlev_var_493
      DO jcol_var_494 = istartcol_var_491, iendcol_var_492
        sum_mixing_ratio(jcol_var_494) = 0.0D0
      END DO
      DO jh = 1, 2
        DO jcol_var_494 = istartcol_var_491, iendcol_var_492
          sum_mixing_ratio(jcol_var_494) = sum_mixing_ratio(jcol_var_494) + this_var_490 % mixing_ratio(jcol_var_494, jlev_var_495, jh)
        END DO
      END DO
      DO jcol_var_494 = istartcol_var_491, iendcol_var_492
        IF (this_var_490 % fraction(jcol_var_494, jlev_var_495) < cloud_fraction_threshold .OR. sum_mixing_ratio(jcol_var_494) < cloud_mixing_ratio_threshold) THEN
          this_var_490 % fraction(jcol_var_494, jlev_var_495) = 0.0D0
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
  ELEMENTAL SUBROUTINE delta_eddington_scat_od(od_var_496, scat_od_var_497, g_var_498)
    REAL(KIND = 8), INTENT(INOUT) :: od_var_496, scat_od_var_497, g_var_498
    REAL(KIND = 8) :: f_var_499
    f_var_499 = g_var_498 * g_var_498
    od_var_496 = od_var_496 - scat_od_var_497 * f_var_499
    scat_od_var_497 = scat_od_var_497 * (1.0D0 - f_var_499)
    g_var_498 = g_var_498 / (1.0D0 + g_var_498)
  END SUBROUTINE delta_eddington_scat_od
  SUBROUTINE cloud_optics_fn_500(nlev_var_500, istartcol_var_501, iendcol_var_502, config_var_503, thermodynamics_var_504, cloud_var_505, od_lw_cloud_var_506, ssa_lw_cloud_var_507, g_lw_cloud_var_508, od_sw_cloud_var_509, ssa_sw_cloud_var_510, g_sw_cloud_var_511)
    USE radiation_config, ONLY: config_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE radiation_cloud, ONLY: cloud_type
    USE radiation_liquid_optics_socrates, ONLY: calc_liq_optics_socrates
    USE radiation_ice_optics_fu, ONLY: calc_ice_optics_fu_lw, calc_ice_optics_fu_sw
    USE radiation_ice_optics_yi, ONLY: calc_ice_optics_yi_lw, calc_ice_optics_yi_sw
    INTEGER, INTENT(IN) :: nlev_var_500
    INTEGER, INTENT(IN) :: istartcol_var_501, iendcol_var_502
    TYPE(config_type), INTENT(IN), TARGET :: config_var_503
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics_var_504
    TYPE(cloud_type), INTENT(IN) :: cloud_var_505
    REAL(KIND = 8), DIMENSION(16, nlev_var_500, istartcol_var_501 : iendcol_var_502), INTENT(OUT) :: od_lw_cloud_var_506
    REAL(KIND = 8), DIMENSION(0, nlev_var_500, istartcol_var_501 : iendcol_var_502), INTENT(OUT) :: ssa_lw_cloud_var_507, g_lw_cloud_var_508
    REAL(KIND = 8), DIMENSION(14, nlev_var_500, istartcol_var_501 : iendcol_var_502), INTENT(OUT) :: od_sw_cloud_var_509, ssa_sw_cloud_var_510, g_sw_cloud_var_511
    REAL(KIND = 8), DIMENSION(16) :: od_lw_liq, scat_od_lw_liq, g_lw_liq, od_lw_ice, scat_od_lw_ice, g_lw_ice
    REAL(KIND = 8), DIMENSION(14) :: od_sw_liq, scat_od_sw_liq, g_sw_liq, od_sw_ice, scat_od_sw_ice, g_sw_ice
    REAL(KIND = 8) :: lwp_in_cloud, iwp_in_cloud
    REAL(KIND = 8) :: factor_var_512
    INTEGER :: jcol_var_513, jlev_var_514, jb_var_515
    DO jcol_var_513 = istartcol_var_501, iendcol_var_502
      DO jlev_var_514 = 1, nlev_var_500
        DO jb_var_515 = 1, 14
          od_sw_cloud_var_509(jb_var_515, jlev_var_514, jcol_var_513) = 0.0D0
          ssa_sw_cloud_var_510(jb_var_515, jlev_var_514, jcol_var_513) = 0.0D0
          g_sw_cloud_var_511(jb_var_515, jlev_var_514, jcol_var_513) = 0.0D0
        END DO
        DO jb_var_515 = 1, 16
          od_lw_cloud_var_506(jb_var_515, jlev_var_514, jcol_var_513) = 0.0D0
        END DO
        DO jb_var_515 = 1, 0
          ssa_lw_cloud_var_507(jb_var_515, jlev_var_514, jcol_var_513) = 0.0D0
          g_lw_cloud_var_508(jb_var_515, jlev_var_514, jcol_var_513) = 0.0D0
        END DO
      END DO
    END DO
    DO jlev_var_514 = 1, nlev_var_500
      DO jcol_var_513 = istartcol_var_501, iendcol_var_502
        IF (cloud_var_505 % fraction(jcol_var_513, jlev_var_514) > 0.0D0) THEN
          factor_var_512 = (thermodynamics_var_504 % pressure_hl(jcol_var_513, jlev_var_514 + 1) - thermodynamics_var_504 % pressure_hl(jcol_var_513, jlev_var_514)) / (9.80665D0 * cloud_var_505 % fraction(jcol_var_513, jlev_var_514))
          lwp_in_cloud = factor_var_512 * cloud_var_505 % q_liq(jcol_var_513, jlev_var_514)
          iwp_in_cloud = factor_var_512 * cloud_var_505 % q_ice(jcol_var_513, jlev_var_514)
          IF (lwp_in_cloud > 0.0D0) THEN
            CALL calc_liq_optics_socrates(16, config_var_503 % cloud_optics % liq_coeff_lw, lwp_in_cloud, cloud_var_505 % re_liq(jcol_var_513, jlev_var_514), od_lw_liq, scat_od_lw_liq, g_lw_liq)
            CALL calc_liq_optics_socrates(14, config_var_503 % cloud_optics % liq_coeff_sw, lwp_in_cloud, cloud_var_505 % re_liq(jcol_var_513, jlev_var_514), od_sw_liq, scat_od_sw_liq, g_sw_liq)
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
            CALL calc_ice_optics_fu_lw(16, config_var_503 % cloud_optics % ice_coeff_lw, iwp_in_cloud, cloud_var_505 % re_ice(jcol_var_513, jlev_var_514), od_lw_ice, scat_od_lw_ice, g_lw_ice)
            CALL calc_ice_optics_fu_sw(14, config_var_503 % cloud_optics % ice_coeff_sw, iwp_in_cloud, cloud_var_505 % re_ice(jcol_var_513, jlev_var_514), od_sw_ice, scat_od_sw_ice, g_sw_ice)
            CALL calc_ice_optics_yi_lw(16, config_var_503 % cloud_optics % ice_coeff_lw, iwp_in_cloud, cloud_var_505 % re_ice(jcol_var_513, jlev_var_514), od_lw_ice, scat_od_lw_ice, g_lw_ice)
            CALL calc_ice_optics_yi_sw(14, config_var_503 % cloud_optics % ice_coeff_sw, iwp_in_cloud, cloud_var_505 % re_ice(jcol_var_513, jlev_var_514), od_sw_ice, scat_od_sw_ice, g_sw_ice)
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
          DO jb_var_515 = 1, 16
            od_lw_cloud_var_506(jb_var_515, jlev_var_514, jcol_var_513) = od_lw_liq(jb_var_515) - scat_od_lw_liq(jb_var_515) + od_lw_ice(jb_var_515) - scat_od_lw_ice(jb_var_515)
          END DO
          DO jb_var_515 = 1, 14
            od_sw_cloud_var_509(jb_var_515, jlev_var_514, jcol_var_513) = od_sw_liq(jb_var_515) + od_sw_ice(jb_var_515)
            g_sw_cloud_var_511(jb_var_515, jlev_var_514, jcol_var_513) = (g_sw_liq(jb_var_515) * scat_od_sw_liq(jb_var_515) + g_sw_ice(jb_var_515) * scat_od_sw_ice(jb_var_515)) / (scat_od_sw_liq(jb_var_515) + scat_od_sw_ice(jb_var_515))
            ssa_sw_cloud_var_510(jb_var_515, jlev_var_514, jcol_var_513) = (scat_od_sw_liq(jb_var_515) + scat_od_sw_ice(jb_var_515)) / (od_sw_liq(jb_var_515) + od_sw_ice(jb_var_515))
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
  SUBROUTINE gas_optics(ncol_var_516, nlev_var_517, istartcol_var_518, iendcol_var_519, config_var_520, single_level_var_521, thermodynamics_var_522, gas_var_523, od_lw_var_525, od_sw_var_526, ssa_sw_var_527, lw_albedo_var_524, planck_hl_var_528, lw_emission_var_529, incoming_sw_var_530)
    USE radiation_config, ONLY: config_type
    USE radiation_single_level, ONLY: single_level_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE radiation_gas, ONLY: assert_units_gas, gas_type
    INTEGER, INTENT(IN) :: ncol_var_516
    INTEGER, INTENT(IN) :: nlev_var_517
    INTEGER, INTENT(IN) :: istartcol_var_518, iendcol_var_519
    TYPE(config_type), INTENT(IN) :: config_var_520
    TYPE(single_level_type), INTENT(IN) :: single_level_var_521
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics_var_522
    TYPE(gas_type), INTENT(IN) :: gas_var_523
    REAL(KIND = 8), DIMENSION(140, istartcol_var_518 : iendcol_var_519), INTENT(IN), OPTIONAL :: lw_albedo_var_524
    REAL(KIND = 8), DIMENSION(140, nlev_var_517, istartcol_var_518 : iendcol_var_519), INTENT(OUT) :: od_lw_var_525
    REAL(KIND = 8), DIMENSION(112, nlev_var_517, istartcol_var_518 : iendcol_var_519), INTENT(OUT) :: od_sw_var_526, ssa_sw_var_527
    REAL(KIND = 8), DIMENSION(140, nlev_var_517 + 1, istartcol_var_518 : iendcol_var_519), INTENT(OUT), OPTIONAL :: planck_hl_var_528
    REAL(KIND = 8), DIMENSION(140, istartcol_var_518 : iendcol_var_519), INTENT(OUT), OPTIONAL :: lw_emission_var_529
    REAL(KIND = 8), DIMENSION(112, istartcol_var_518 : iendcol_var_519), INTENT(OUT), OPTIONAL :: incoming_sw_var_530
    REAL(KIND = 8) :: incoming_sw_scale(istartcol_var_518 : iendcol_var_519)
    REAL(KIND = 8) :: zod_lw(140, nlev_var_517, istartcol_var_518 : iendcol_var_519)
    REAL(KIND = 8) :: zod_sw(istartcol_var_518 : iendcol_var_519, nlev_var_517, 112)
    REAL(KIND = 8) :: zssa_sw(istartcol_var_518 : iendcol_var_519, nlev_var_517, 112)
    REAL(KIND = 8) :: zincsol(istartcol_var_518 : iendcol_var_519, 112)
    REAL(KIND = 8) :: zcolmol(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zcoldry(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zwbrodl(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zcolbrd(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zwkl(istartcol_var_518 : iendcol_var_519, 35, nlev_var_517)
    REAL(KIND = 8) :: zwx(istartcol_var_518 : iendcol_var_519, 4, nlev_var_517)
    REAL(KIND = 8) :: zfluxfac_var_531, zpi
    REAL(KIND = 8) :: ztauaerl(istartcol_var_518 : iendcol_var_519, nlev_var_517, 16)
    REAL(KIND = 8) :: zfac00(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zfac01(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zfac10(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zfac11(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zforfac(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zforfrac(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    INTEGER :: indfor_var_532(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    INTEGER :: indminor_var_533(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zscaleminor(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zscaleminorn2(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zminorfrac(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zrat_h2oco2(istartcol_var_518 : iendcol_var_519, nlev_var_517), zrat_h2oco2_1(istartcol_var_518 : iendcol_var_519, nlev_var_517), zrat_h2oo3(istartcol_var_518 : iendcol_var_519, nlev_var_517), zrat_h2oo3_1(istartcol_var_518 : iendcol_var_519, nlev_var_517), zrat_h2on2o(istartcol_var_518 : iendcol_var_519, nlev_var_517), zrat_h2on2o_1(istartcol_var_518 : iendcol_var_519, nlev_var_517), zrat_h2och4(istartcol_var_518 : iendcol_var_519, nlev_var_517), zrat_h2och4_1(istartcol_var_518 : iendcol_var_519, nlev_var_517), zrat_n2oco2(istartcol_var_518 : iendcol_var_519, nlev_var_517), zrat_n2oco2_1(istartcol_var_518 : iendcol_var_519, nlev_var_517), zrat_o3co2(istartcol_var_518 : iendcol_var_519, nlev_var_517), zrat_o3co2_1(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    INTEGER :: jp_var_534(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    INTEGER :: jt_var_535(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    INTEGER :: jt1_var_536(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zoneminus, zoneminus_array(istartcol_var_518 : iendcol_var_519)
    REAL(KIND = 8) :: zcolh2o(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zcolco2(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zcolo3(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zcoln2o(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zcolch4(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zcolo2(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zco2mult(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    INTEGER :: ilaytrop(istartcol_var_518 : iendcol_var_519)
    INTEGER :: ilayswtch(istartcol_var_518 : iendcol_var_519)
    INTEGER :: ilaylow(istartcol_var_518 : iendcol_var_519)
    REAL(KIND = 8) :: zpavel(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: ztavel(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zpz(istartcol_var_518 : iendcol_var_519, 0 : nlev_var_517)
    REAL(KIND = 8) :: ztz(istartcol_var_518 : iendcol_var_519, 0 : nlev_var_517)
    REAL(KIND = 8) :: zselffac(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zselffrac(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    INTEGER :: indself_var_537(istartcol_var_518 : iendcol_var_519, nlev_var_517)
    REAL(KIND = 8) :: zpfrac(istartcol_var_518 : iendcol_var_519, 140, nlev_var_517)
    INTEGER :: ireflect(istartcol_var_518 : iendcol_var_519)
    REAL(KIND = 8) :: pressure_fl_var_538(ncol_var_516, nlev_var_517), temperature_fl_var_539(ncol_var_516, nlev_var_517)
    INTEGER :: istartlev_var_540, iendlev_var_541
    LOGICAL :: do_sw, do_lw
    INTEGER :: jlev_var_542, jg_var_543, jcol_var_544
    do_sw = .TRUE.
    do_lw = .TRUE.
    iendlev_var_541 = UBOUND(gas_var_523 % mixing_ratio, 2)
    istartlev_var_540 = iendlev_var_541 - nlev_var_517 + 1
    zpi = 1.5707963267948966D0
    zfluxfac_var_531 = 15707.963267948966D0
    zoneminus = 0.999999D0
    DO jcol_var_544 = istartcol_var_518, iendcol_var_519
      zoneminus_array(jcol_var_544) = zoneminus
    END DO
    DO jlev_var_542 = 1, nlev_var_517
      DO jcol_var_544 = istartcol_var_518, iendcol_var_519
        pressure_fl_var_538(jcol_var_544, jlev_var_542) = thermodynamics_var_522 % pressure_fl(jcol_var_544, jlev_var_542)
        temperature_fl_var_539(jcol_var_544, jlev_var_542) = thermodynamics_var_522 % temperature_fl(jcol_var_544, jlev_var_542)
      END DO
    END DO
    CALL assert_units_gas(gas_var_523, 0)
    CALL rrtm_prepare_gases(istartcol_var_518, iendcol_var_519, ncol_var_516, nlev_var_517, thermodynamics_var_522 % pressure_hl(:, istartlev_var_540 : iendlev_var_541 + 1), pressure_fl_var_538, thermodynamics_var_522 % temperature_hl(:, istartlev_var_540 : iendlev_var_541 + 1), temperature_fl_var_539, gas_var_523 % mixing_ratio(:, istartlev_var_540 : iendlev_var_541, 1), gas_var_523 % mixing_ratio(:, istartlev_var_540 : iendlev_var_541, 2), gas_var_523 % mixing_ratio(:, istartlev_var_540 : iendlev_var_541, 6), gas_var_523 % mixing_ratio(:, istartlev_var_540 : iendlev_var_541, 4), gas_var_523 % mixing_ratio(:, istartlev_var_540 : iendlev_var_541, 12), gas_var_523 % mixing_ratio(:, istartlev_var_540 : iendlev_var_541, 8), gas_var_523 % mixing_ratio(:, istartlev_var_540 : iendlev_var_541, 9), gas_var_523 % mixing_ratio(:, istartlev_var_540 : iendlev_var_541, 10), gas_var_523 % mixing_ratio(:, istartlev_var_540 : iendlev_var_541, 11), gas_var_523 % mixing_ratio(:, istartlev_var_540 : iendlev_var_541, 3), zcoldry, zwbrodl, zwkl, zwx, zpavel, ztavel, zpz, ztz, ireflect)
    CALL rrtm_setcoef_140gp(istartcol_var_518, iendcol_var_519, nlev_var_517, zcoldry, zwbrodl, zwkl, zfac00, zfac01, zfac10, zfac11, zforfac, zforfrac, indfor_var_532, jp_var_534, jt_var_535, jt1_var_536, zcolh2o, zcolco2, zcolo3, zcoln2o, zcolch4, zcolo2, zco2mult, zcolbrd, ilaytrop, ilayswtch, ilaylow, zpavel, ztavel, zselffac, zselffrac, indself_var_537, indminor_var_533, zscaleminor, zscaleminorn2, zminorfrac, zrat_h2oco2, zrat_h2oco2_1, zrat_h2oo3, zrat_h2oo3_1, zrat_h2on2o, zrat_h2on2o_1, zrat_h2och4, zrat_h2och4_1, zrat_n2oco2, zrat_n2oco2_1, zrat_o3co2, zrat_o3co2_1)
    DO jg_var_543 = 1, 16
      DO jlev_var_542 = 1, nlev_var_517
        DO jcol_var_544 = istartcol_var_518, iendcol_var_519
          ztauaerl(jcol_var_544, jlev_var_542, jg_var_543) = 0.0D0
        END DO
      END DO
    END DO
    CALL rrtm_gas_optical_depth(istartcol_var_518, iendcol_var_519, nlev_var_517, zod_lw, zpavel, zcoldry, zcolbrd, zwx, ztauaerl, zfac00, zfac01, zfac10, zfac11, zforfac, zforfrac, indfor_var_532, jp_var_534, jt_var_535, jt1_var_536, zoneminus, zcolh2o, zcolco2, zcolo3, zcoln2o, zcolch4, zcolo2, zco2mult, ilaytrop, ilayswtch, ilaylow, zselffac, zselffrac, indself_var_537, zpfrac, indminor_var_533, zscaleminor, zscaleminorn2, zminorfrac, zrat_h2oco2, zrat_h2oco2_1, zrat_h2oo3, zrat_h2oo3_1, zrat_h2on2o, zrat_h2on2o_1, zrat_h2och4, zrat_h2och4_1, zrat_n2oco2, zrat_n2oco2_1, zrat_o3co2, zrat_o3co2_1)
    CALL planck_function_atmos(nlev_var_517, istartcol_var_518, iendcol_var_519, config_var_520, thermodynamics_var_522, zpfrac, planck_hl_var_528)
    CALL planck_function_surf(istartcol_var_518, iendcol_var_519, config_var_520, single_level_var_521 % skin_temperature, zpfrac(:, :, 1), lw_emission_var_529)
    DO jcol_var_544 = istartcol_var_518, iendcol_var_519
      DO jg_var_543 = 1, 140
        lw_emission_var_529(jg_var_543, jcol_var_544) = lw_emission_var_529(jg_var_543, jcol_var_544) * (1.0D0 - lw_albedo_var_524(jg_var_543, jcol_var_544))
      END DO
    END DO
    DO jcol_var_544 = istartcol_var_518, iendcol_var_519
      DO jlev_var_542 = 1, nlev_var_517
        DO jg_var_543 = 1, 140
          od_lw_var_525(jg_var_543, jlev_var_542, jcol_var_544) = MAX(1D-15, zod_lw(jg_var_543, nlev_var_517 + 1 - jlev_var_542, jcol_var_544))
        END DO
      END DO
    END DO
    CALL srtm_setcoef(istartcol_var_518, iendcol_var_519, nlev_var_517, zpavel, ztavel, zcoldry, zwkl, ilaytrop, zcolch4, zcolco2, zcolh2o, zcolmol, zcolo2, zcolo3, zforfac, zforfrac, indfor_var_532, zselffac, zselffrac, indself_var_537, zfac00, zfac01, zfac10, zfac11, jp_var_534, jt_var_535, jt1_var_536, single_level_var_521 % cos_sza(istartcol_var_518 : iendcol_var_519))
    DO jg_var_543 = 1, 112
      DO jlev_var_542 = 1, nlev_var_517
        DO jcol_var_544 = istartcol_var_518, iendcol_var_519
          zod_sw(jcol_var_544, jlev_var_542, jg_var_543) = 0.0D0
          zssa_sw(jcol_var_544, jlev_var_542, jg_var_543) = 0.0D0
        END DO
      END DO
    END DO
    DO jg_var_543 = 1, 112
      DO jcol_var_544 = istartcol_var_518, iendcol_var_519
        zincsol(jcol_var_544, jg_var_543) = 0.0D0
      END DO
    END DO
    CALL srtm_gas_optical_depth(istartcol_var_518, iendcol_var_519, nlev_var_517, zoneminus_array, single_level_var_521 % cos_sza(istartcol_var_518 : iendcol_var_519), ilaytrop, zcolch4, zcolco2, zcolh2o, zcolmol, zcolo2, zcolo3, zforfac, zforfrac, indfor_var_532, zselffac, zselffrac, indself_var_537, zfac00, zfac01, zfac10, zfac11, jp_var_534, jt_var_535, jt1_var_536, zod_sw, zssa_sw, zincsol)
    DO jg_var_543 = 1, 112
      DO jcol_var_544 = istartcol_var_518, iendcol_var_519
        zincsol(jcol_var_544, jg_var_543) = zincsol(jcol_var_544, jg_var_543) * single_level_var_521 % spectral_solar_scaling(config_var_520 % i_band_from_reordered_g_sw(jg_var_543))
      END DO
    END DO
    DO jcol_var_544 = istartcol_var_518, iendcol_var_519
      IF (single_level_var_521 % cos_sza(jcol_var_544) > 0.0D0) THEN
        incoming_sw_scale(jcol_var_544) = 1.0D0 / SUM(zincsol(jcol_var_544, :))
      ELSE
        incoming_sw_scale(jcol_var_544) = 1.0D0
      END IF
    END DO
    DO jcol_var_544 = istartcol_var_518, iendcol_var_519
      DO jlev_var_542 = 1, nlev_var_517
        DO jg_var_543 = 1, 112
          od_sw_var_526(jg_var_543, nlev_var_517 + 1 - jlev_var_542, jcol_var_544) = MAX(0.0D0, zod_sw(jcol_var_544, jlev_var_542, jg_var_543))
          ssa_sw_var_527(jg_var_543, nlev_var_517 + 1 - jlev_var_542, jcol_var_544) = zssa_sw(jcol_var_544, jlev_var_542, jg_var_543)
        END DO
      END DO
      DO jg_var_543 = 1, 112
        incoming_sw_var_530(jg_var_543, jcol_var_544) = incoming_sw_scale(jcol_var_544) * zincsol(jcol_var_544, jg_var_543)
      END DO
    END DO
  END SUBROUTINE gas_optics
  SUBROUTINE planck_function_atmos(nlev_var_545, istartcol_var_546, iendcol_var_547, config_var_548, thermodynamics_var_549, pfrac_var_550, planck_hl_var_551)
    USE radiation_config, ONLY: config_type
    USE radiation_thermodynamics, ONLY: thermodynamics_type
    USE yoerrtwn, ONLY: delwave, totplnk
    INTEGER, INTENT(IN) :: nlev_var_545
    INTEGER, INTENT(IN) :: istartcol_var_546, iendcol_var_547
    TYPE(config_type), INTENT(IN) :: config_var_548
    TYPE(thermodynamics_type), INTENT(IN) :: thermodynamics_var_549
    REAL(KIND = 8), INTENT(IN) :: pfrac_var_550(istartcol_var_546 : iendcol_var_547, 140, nlev_var_545)
    REAL(KIND = 8), DIMENSION(140, nlev_var_545 + 1, istartcol_var_546 : iendcol_var_547), INTENT(OUT) :: planck_hl_var_551
    REAL(KIND = 8), DIMENSION(istartcol_var_546 : iendcol_var_547, 16) :: planck_store_var_552
    REAL(KIND = 8), DIMENSION(istartcol_var_546 : iendcol_var_547) :: frac_var_553
    INTEGER, DIMENSION(istartcol_var_546 : iendcol_var_547) :: ind_var_554
    REAL(KIND = 8) :: temperature_var_555
    REAL(KIND = 8) :: factor_var_556, planck_tmp(istartcol_var_546 : iendcol_var_547, 140)
    REAL(KIND = 8) :: zfluxfac_var_557
    INTEGER :: jlev_var_558, jg_var_559, iband_var_560, jband_var_561, jcol_var_562, ilevoffset
    zfluxfac_var_557 = 15707.963267948966D0
    ilevoffset = UBOUND(thermodynamics_var_549 % temperature_hl, 2) - nlev_var_545 - 1
    DO jlev_var_558 = 1, nlev_var_545 + 1
      DO jcol_var_562 = istartcol_var_546, iendcol_var_547
        temperature_var_555 = thermodynamics_var_549 % temperature_hl(jcol_var_562, jlev_var_558 + ilevoffset)
        IF (temperature_var_555 < 339.0D0 .AND. temperature_var_555 >= 160.0D0) THEN
          ind_var_554(jcol_var_562) = INT(temperature_var_555 - 159.0D0)
          frac_var_553(jcol_var_562) = temperature_var_555 - INT(temperature_var_555)
        ELSE IF (temperature_var_555 >= 339.0D0) THEN
          ind_var_554(jcol_var_562) = 180
          frac_var_553(jcol_var_562) = temperature_var_555 - 339.0D0
        ELSE
          ind_var_554(jcol_var_562) = 1
          frac_var_553(jcol_var_562) = 0.0D0
        END IF
      END DO
      DO jband_var_561 = 1, 16
        factor_var_556 = zfluxfac_var_557 * delwave(jband_var_561)
        DO jcol_var_562 = istartcol_var_546, iendcol_var_547
          planck_store_var_552(jcol_var_562, jband_var_561) = factor_var_556 * (totplnk(ind_var_554(jcol_var_562), jband_var_561) + frac_var_553(jcol_var_562) * (totplnk(ind_var_554(jcol_var_562) + 1, jband_var_561) - totplnk(ind_var_554(jcol_var_562), jband_var_561)))
        END DO
      END DO
      IF (jlev_var_558 == 1) THEN
        DO jg_var_559 = 1, 140
          iband_var_560 = config_var_548 % i_band_from_g_lw(jg_var_559)
          DO jcol_var_562 = istartcol_var_546, iendcol_var_547
            planck_hl_var_551(jg_var_559, 1, jcol_var_562) = planck_store_var_552(jcol_var_562, iband_var_560) * pfrac_var_550(jcol_var_562, jg_var_559, nlev_var_545)
          END DO
        END DO
      ELSE
        DO jg_var_559 = 1, 140
          iband_var_560 = config_var_548 % i_band_from_g_lw(jg_var_559)
          DO jcol_var_562 = istartcol_var_546, iendcol_var_547
            planck_tmp(jcol_var_562, jg_var_559) = planck_store_var_552(jcol_var_562, iband_var_560) * pfrac_var_550(jcol_var_562, jg_var_559, nlev_var_545 + 2 - jlev_var_558)
          END DO
        END DO
        DO jcol_var_562 = istartcol_var_546, iendcol_var_547
          DO jg_var_559 = 1, 140
            planck_hl_var_551(jg_var_559, jlev_var_558, jcol_var_562) = planck_tmp(jcol_var_562, jg_var_559)
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE planck_function_atmos
  SUBROUTINE planck_function_surf(istartcol_var_563, iendcol_var_564, config_var_565, temperature_var_566, pfrac_var_567, planck_surf)
    USE radiation_config, ONLY: config_type
    USE yoerrtwn, ONLY: delwave, totplnk
    INTEGER, INTENT(IN) :: istartcol_var_563, iendcol_var_564
    TYPE(config_type), INTENT(IN) :: config_var_565
    REAL(KIND = 8), INTENT(IN) :: temperature_var_566(:)
    REAL(KIND = 8), INTENT(IN) :: pfrac_var_567(istartcol_var_563 : iendcol_var_564, 140)
    REAL(KIND = 8), DIMENSION(140, istartcol_var_563 : iendcol_var_564), INTENT(OUT) :: planck_surf
    REAL(KIND = 8), DIMENSION(istartcol_var_563 : iendcol_var_564, 16) :: planck_store_var_568
    REAL(KIND = 8), DIMENSION(istartcol_var_563 : iendcol_var_564) :: frac_var_569
    INTEGER, DIMENSION(istartcol_var_563 : iendcol_var_564) :: ind_var_570
    REAL(KIND = 8) :: tsurf
    REAL(KIND = 8) :: factor_var_571
    REAL(KIND = 8) :: zfluxfac_var_572
    INTEGER :: jg_var_573, iband_var_574, jband_var_575, jcol_var_576
    zfluxfac_var_572 = 15707.963267948966D0
    DO jcol_var_576 = istartcol_var_563, iendcol_var_564
      tsurf = temperature_var_566(jcol_var_576)
      IF (tsurf < 339.0D0 .AND. tsurf >= 160.0D0) THEN
        ind_var_570(jcol_var_576) = INT(tsurf - 159.0D0)
        frac_var_569(jcol_var_576) = tsurf - INT(tsurf)
      ELSE IF (tsurf >= 339.0D0) THEN
        ind_var_570(jcol_var_576) = 180
        frac_var_569(jcol_var_576) = tsurf - 339.0D0
      ELSE
        ind_var_570(jcol_var_576) = 1
        frac_var_569(jcol_var_576) = 0.0D0
      END IF
    END DO
    DO jband_var_575 = 1, 16
      factor_var_571 = zfluxfac_var_572 * delwave(jband_var_575)
      DO jcol_var_576 = istartcol_var_563, iendcol_var_564
        planck_store_var_568(jcol_var_576, jband_var_575) = factor_var_571 * (totplnk(ind_var_570(jcol_var_576), jband_var_575) + frac_var_569(jcol_var_576) * (totplnk(ind_var_570(jcol_var_576) + 1, jband_var_575) - totplnk(ind_var_570(jcol_var_576), jband_var_575)))
      END DO
    END DO
    DO jg_var_573 = 1, 140
      iband_var_574 = config_var_565 % i_band_from_g_lw(jg_var_573)
      DO jcol_var_576 = istartcol_var_563, iendcol_var_564
        planck_surf(jg_var_573, jcol_var_576) = planck_store_var_568(jcol_var_576, iband_var_574) * pfrac_var_567(jcol_var_576, jg_var_573)
      END DO
    END DO
  END SUBROUTINE planck_function_surf
END MODULE radiation_ifs_rrtm
MODULE radiation_mcica_lw
  CONTAINS
  SUBROUTINE solver_mcica_lw(nlev_var_577, istartcol_var_578, iendcol_var_579, config_var_580, single_level_var_581, cloud_var_582, od_var_583, ssa_var_584, g_var_585, od_cloud_var_586, ssa_cloud_var_587, g_cloud_var_588, planck_hl_var_589, emission, albedo_var_590, flux_var_591)
    USE radiation_config, ONLY: config_type
    USE radiation_single_level, ONLY: single_level_type
    USE radiation_cloud, ONLY: cloud_type
    USE radiation_flux, ONLY: flux_type
    USE radiation_two_stream, ONLY: calc_no_scattering_transmittance_lw
    USE radiation_adding_ica_lw, ONLY: calc_fluxes_no_scattering_lw
    USE radiation_cloud_generator, ONLY: cloud_generator
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlev_var_577
    INTEGER, INTENT(IN) :: istartcol_var_578, iendcol_var_579
    TYPE(config_type), INTENT(IN) :: config_var_580
    TYPE(single_level_type), INTENT(IN) :: single_level_var_581
    TYPE(cloud_type), INTENT(IN) :: cloud_var_582
    REAL(KIND = 8), INTENT(IN), DIMENSION(140, nlev_var_577, istartcol_var_578 : iendcol_var_579) :: od_var_583
    REAL(KIND = 8), INTENT(IN), DIMENSION(0, nlev_var_577, istartcol_var_578 : iendcol_var_579) :: ssa_var_584, g_var_585
    REAL(KIND = 8), INTENT(IN), DIMENSION(16, nlev_var_577, istartcol_var_578 : iendcol_var_579) :: od_cloud_var_586
    REAL(KIND = 8), INTENT(IN), DIMENSION(0, nlev_var_577, istartcol_var_578 : iendcol_var_579) :: ssa_cloud_var_587, g_cloud_var_588
    REAL(KIND = 8), INTENT(IN), DIMENSION(140, nlev_var_577 + 1, istartcol_var_578 : iendcol_var_579) :: planck_hl_var_589
    REAL(KIND = 8), INTENT(IN), DIMENSION(140, istartcol_var_578 : iendcol_var_579) :: emission, albedo_var_590
    TYPE(flux_type), INTENT(INOUT) :: flux_var_591
    REAL(KIND = 8), DIMENSION(140, nlev_var_577) :: ref_clear_var_592, trans_clear_var_593, reflectance_var_594, transmittance_var_595
    REAL(KIND = 8), DIMENSION(140, nlev_var_577) :: source_up_clear, source_dn_clear, source_up_var_596, source_dn_var_597
    REAL(KIND = 8), DIMENSION(140, nlev_var_577 + 1) :: flux_up_var_598, flux_dn_var_599
    REAL(KIND = 8), DIMENSION(140, nlev_var_577 + 1) :: flux_up_clear, flux_dn_clear
    REAL(KIND = 8), DIMENSION(140) :: od_total_var_600, ssa_total_var_601, g_total_var_602
    REAL(KIND = 8), DIMENSION(140, nlev_var_577) :: od_scaling_var_603
    REAL(KIND = 8), DIMENSION(140) :: od_cloud_new_var_604
    REAL(KIND = 8) :: total_cloud_cover_var_605
    LOGICAL :: is_clear_sky_layer(nlev_var_577)
    REAL(KIND = 8) :: sum_up_var_606, sum_dn
    INTEGER :: i_cloud_top
    INTEGER :: ng_var_607
    INTEGER :: jlev_var_608, jcol_var_609, jg_var_610
    ng_var_607 = 140
    DO jcol_var_609 = istartcol_var_578, iendcol_var_579
      CALL calc_no_scattering_transmittance_lw(ng_var_607 * nlev_var_577, od_var_583(:, :, jcol_var_609), planck_hl_var_589(:, 1 : nlev_var_577, jcol_var_609), planck_hl_var_589(:, 2 : nlev_var_577 + 1, jcol_var_609), trans_clear_var_593, source_up_clear, source_dn_clear)
      ref_clear_var_592 = 0.0D0
      CALL calc_fluxes_no_scattering_lw(ng_var_607, nlev_var_577, trans_clear_var_593, source_up_clear, source_dn_clear, emission(:, jcol_var_609), albedo_var_590(:, jcol_var_609), flux_up_clear, flux_dn_clear)
      DO jlev_var_608 = 1, nlev_var_577 + 1
        sum_up_var_606 = 0.0D0
        sum_dn = 0.0D0
        DO jg_var_610 = 1, ng_var_607
          sum_up_var_606 = sum_up_var_606 + flux_up_clear(jg_var_610, jlev_var_608)
          sum_dn = sum_dn + flux_dn_clear(jg_var_610, jlev_var_608)
        END DO
        flux_var_591 % lw_up_clear(jcol_var_609, jlev_var_608) = sum_up_var_606
        flux_var_591 % lw_dn_clear(jcol_var_609, jlev_var_608) = sum_dn
      END DO
      flux_var_591 % lw_dn_surf_clear_g(:, jcol_var_609) = flux_dn_clear(:, nlev_var_577 + 1)
      CALL cloud_generator(ng_var_607, nlev_var_577, 1, single_level_var_581 % iseed(jcol_var_609) + 997, 1D-06, cloud_var_582 % fraction(jcol_var_609, :), cloud_var_582 % overlap_param(jcol_var_609, :), 0.5D0, cloud_var_582 % fractional_std(jcol_var_609, :), config_var_580 % pdf_sampler, od_scaling_var_603, total_cloud_cover_var_605, use_beta_overlap = .FALSE., use_vectorizable_generator = .FALSE.)
      flux_var_591 % cloud_cover_lw(jcol_var_609) = total_cloud_cover_var_605
      IF (total_cloud_cover_var_605 >= 1D-06) THEN
        is_clear_sky_layer = .TRUE.
        i_cloud_top = nlev_var_577 + 1
        DO jlev_var_608 = 1, nlev_var_577
          IF (cloud_var_582 % fraction(jcol_var_609, jlev_var_608) >= 1D-06) THEN
            is_clear_sky_layer(jlev_var_608) = .FALSE.
            IF (i_cloud_top > jlev_var_608) THEN
              i_cloud_top = jlev_var_608
            END IF
            DO jg_var_610 = 1, ng_var_607
              od_cloud_new_var_604(jg_var_610) = od_scaling_var_603(jg_var_610, jlev_var_608) * od_cloud_var_586(config_var_580 % i_band_from_reordered_g_lw(jg_var_610), jlev_var_608, jcol_var_609)
              od_total_var_600(jg_var_610) = od_var_583(jg_var_610, jlev_var_608, jcol_var_609) + od_cloud_new_var_604(jg_var_610)
              ssa_total_var_601(jg_var_610) = 0.0D0
              g_total_var_602(jg_var_610) = 0.0D0
            END DO
            CALL calc_no_scattering_transmittance_lw(ng_var_607, od_total_var_600, planck_hl_var_589(:, jlev_var_608, jcol_var_609), planck_hl_var_589(:, jlev_var_608 + 1, jcol_var_609), transmittance_var_595(:, jlev_var_608), source_up_var_596(:, jlev_var_608), source_dn_var_597(:, jlev_var_608))
          ELSE
            DO jg_var_610 = 1, ng_var_607
              reflectance_var_594(jg_var_610, jlev_var_608) = ref_clear_var_592(jg_var_610, jlev_var_608)
              transmittance_var_595(jg_var_610, jlev_var_608) = trans_clear_var_593(jg_var_610, jlev_var_608)
              source_up_var_596(jg_var_610, jlev_var_608) = source_up_clear(jg_var_610, jlev_var_608)
              source_dn_var_597(jg_var_610, jlev_var_608) = source_dn_clear(jg_var_610, jlev_var_608)
            END DO
          END IF
        END DO
        CALL calc_fluxes_no_scattering_lw(ng_var_607, nlev_var_577, transmittance_var_595, source_up_var_596, source_dn_var_597, emission(:, jcol_var_609), albedo_var_590(:, jcol_var_609), flux_up_var_598, flux_dn_var_599)
        DO jlev_var_608 = 1, nlev_var_577 + 1
          sum_up_var_606 = 0.0D0
          sum_dn = 0.0D0
          DO jg_var_610 = 1, ng_var_607
            sum_up_var_606 = sum_up_var_606 + flux_up_var_598(jg_var_610, jlev_var_608)
            sum_dn = sum_dn + flux_dn_var_599(jg_var_610, jlev_var_608)
          END DO
          flux_var_591 % lw_up(jcol_var_609, jlev_var_608) = sum_up_var_606
          flux_var_591 % lw_dn(jcol_var_609, jlev_var_608) = sum_dn
        END DO
        DO jlev_var_608 = 1, nlev_var_577 + 1
          flux_var_591 % lw_up(jcol_var_609, jlev_var_608) = total_cloud_cover_var_605 * flux_var_591 % lw_up(jcol_var_609, jlev_var_608) + (1.0D0 - total_cloud_cover_var_605) * flux_var_591 % lw_up_clear(jcol_var_609, jlev_var_608)
          flux_var_591 % lw_dn(jcol_var_609, jlev_var_608) = total_cloud_cover_var_605 * flux_var_591 % lw_dn(jcol_var_609, jlev_var_608) + (1.0D0 - total_cloud_cover_var_605) * flux_var_591 % lw_dn_clear(jcol_var_609, jlev_var_608)
        END DO
        flux_var_591 % lw_dn_surf_g(:, jcol_var_609) = total_cloud_cover_var_605 * flux_dn_var_599(:, nlev_var_577 + 1) + (1.0D0 - total_cloud_cover_var_605) * flux_var_591 % lw_dn_surf_clear_g(:, jcol_var_609)
      ELSE
        DO jlev_var_608 = 1, nlev_var_577 + 1
          flux_var_591 % lw_up(jcol_var_609, jlev_var_608) = flux_var_591 % lw_up_clear(jcol_var_609, jlev_var_608)
          flux_var_591 % lw_dn(jcol_var_609, jlev_var_608) = flux_var_591 % lw_dn_clear(jcol_var_609, jlev_var_608)
        END DO
        flux_var_591 % lw_dn_surf_g(:, jcol_var_609) = flux_var_591 % lw_dn_surf_clear_g(:, jcol_var_609)
      END IF
    END DO
  END SUBROUTINE solver_mcica_lw
END MODULE radiation_mcica_lw
MODULE radiation_mcica_sw
  CONTAINS
  SUBROUTINE solver_mcica_sw(nlev_var_611, istartcol_var_612, iendcol_var_613, config_var_614, single_level_var_615, cloud_var_616, od_var_617, ssa_var_618, g_var_619, od_cloud_var_620, ssa_cloud_var_621, g_cloud_var_622, albedo_direct, albedo_diffuse, incoming_sw_var_623, flux_var_624)
    USE radiation_config, ONLY: config_type
    USE radiation_single_level, ONLY: single_level_type
    USE radiation_cloud, ONLY: cloud_type
    USE radiation_flux, ONLY: flux_type
    USE radiation_two_stream, ONLY: calc_ref_trans_sw
    USE radiation_adding_ica_sw, ONLY: adding_ica_sw
    USE radiation_cloud_generator, ONLY: cloud_generator
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlev_var_611
    INTEGER, INTENT(IN) :: istartcol_var_612, iendcol_var_613
    TYPE(config_type), INTENT(IN) :: config_var_614
    TYPE(single_level_type), INTENT(IN) :: single_level_var_615
    TYPE(cloud_type), INTENT(IN) :: cloud_var_616
    REAL(KIND = 8), INTENT(IN), DIMENSION(112, nlev_var_611, istartcol_var_612 : iendcol_var_613) :: od_var_617, ssa_var_618, g_var_619
    REAL(KIND = 8), INTENT(IN), DIMENSION(14, nlev_var_611, istartcol_var_612 : iendcol_var_613) :: od_cloud_var_620, ssa_cloud_var_621, g_cloud_var_622
    REAL(KIND = 8), INTENT(IN), DIMENSION(112, istartcol_var_612 : iendcol_var_613) :: albedo_direct, albedo_diffuse, incoming_sw_var_623
    TYPE(flux_type), INTENT(INOUT) :: flux_var_624
    REAL(KIND = 8) :: cos_sza_var_625
    REAL(KIND = 8), DIMENSION(112, nlev_var_611) :: ref_clear_var_626, trans_clear_var_627, reflectance_var_628, transmittance_var_629
    REAL(KIND = 8), DIMENSION(112, nlev_var_611) :: ref_dir_clear, trans_dir_diff_clear, ref_dir_var_630, trans_dir_diff_var_631
    REAL(KIND = 8), DIMENSION(112, nlev_var_611) :: trans_dir_dir_clear, trans_dir_dir_var_632
    REAL(KIND = 8), DIMENSION(112, nlev_var_611 + 1) :: flux_up_var_633, flux_dn_diffuse_var_634, flux_dn_direct_var_635
    REAL(KIND = 8), DIMENSION(112) :: od_total_var_636, ssa_total_var_637, g_total_var_638
    REAL(KIND = 8) :: scat_od_var_639
    REAL(KIND = 8), DIMENSION(112, nlev_var_611) :: od_scaling_var_640
    REAL(KIND = 8), DIMENSION(112) :: od_cloud_new_var_641
    REAL(KIND = 8), DIMENSION(112, nlev_var_611 + 1) :: tmp_work_albedo, tmp_work_source
    REAL(KIND = 8), DIMENSION(112, nlev_var_611) :: tmp_work_inv_denominator
    REAL(KIND = 8) :: total_cloud_cover_var_642
    REAL(KIND = 8) :: sum_up_var_643, sum_dn_diff, sum_dn_dir
    INTEGER :: ng_var_644
    INTEGER :: jlev_var_645, jcol_var_646, jg_var_647
    ng_var_644 = 112
    DO jcol_var_646 = istartcol_var_612, iendcol_var_613
      IF (single_level_var_615 % cos_sza(jcol_var_646) > 0.0D0) THEN
        cos_sza_var_625 = single_level_var_615 % cos_sza(jcol_var_646)
        CALL calc_ref_trans_sw(ng_var_644 * nlev_var_611, cos_sza_var_625, od_var_617(:, :, jcol_var_646), ssa_var_618(:, :, jcol_var_646), g_var_619(:, :, jcol_var_646), ref_clear_var_626, trans_clear_var_627, ref_dir_clear, trans_dir_diff_clear, trans_dir_dir_clear)
        CALL adding_ica_sw(ng_var_644, nlev_var_611, incoming_sw_var_623(:, jcol_var_646), albedo_diffuse(:, jcol_var_646), albedo_direct(:, jcol_var_646), cos_sza_var_625, ref_clear_var_626, trans_clear_var_627, ref_dir_clear, trans_dir_diff_clear, trans_dir_dir_clear, flux_up_var_633, flux_dn_diffuse_var_634, flux_dn_direct_var_635, albedo_var_20 = tmp_work_albedo, source_var_21 = tmp_work_source, inv_denominator = tmp_work_inv_denominator)
        DO jlev_var_645 = 1, nlev_var_611 + 1
          sum_up_var_643 = 0.0D0
          sum_dn_diff = 0.0D0
          sum_dn_dir = 0.0D0
          DO jg_var_647 = 1, ng_var_644
            sum_up_var_643 = sum_up_var_643 + flux_up_var_633(jg_var_647, jlev_var_645)
            sum_dn_diff = sum_dn_diff + flux_dn_diffuse_var_634(jg_var_647, jlev_var_645)
            sum_dn_dir = sum_dn_dir + flux_dn_direct_var_635(jg_var_647, jlev_var_645)
          END DO
          flux_var_624 % sw_up_clear(jcol_var_646, jlev_var_645) = sum_up_var_643
          flux_var_624 % sw_dn_clear(jcol_var_646, jlev_var_645) = sum_dn_diff + sum_dn_dir
          flux_var_624 % sw_dn_direct_clear(jcol_var_646, jlev_var_645) = sum_dn_dir
        END DO
        DO jg_var_647 = 1, ng_var_644
          flux_var_624 % sw_dn_diffuse_surf_clear_g(jg_var_647, jcol_var_646) = flux_dn_diffuse_var_634(jg_var_647, nlev_var_611 + 1)
          flux_var_624 % sw_dn_direct_surf_clear_g(jg_var_647, jcol_var_646) = flux_dn_direct_var_635(jg_var_647, nlev_var_611 + 1)
        END DO
        CALL cloud_generator(ng_var_644, nlev_var_611, 1, single_level_var_615 % iseed(jcol_var_646), 1D-06, cloud_var_616 % fraction(jcol_var_646, :), cloud_var_616 % overlap_param(jcol_var_646, :), 0.5D0, cloud_var_616 % fractional_std(jcol_var_646, :), config_var_614 % pdf_sampler, od_scaling_var_640, total_cloud_cover_var_642, use_beta_overlap = .FALSE., use_vectorizable_generator = .FALSE.)
        flux_var_624 % cloud_cover_sw(jcol_var_646) = total_cloud_cover_var_642
        IF (total_cloud_cover_var_642 >= 1D-06) THEN
          DO jlev_var_645 = 1, nlev_var_611
            IF (cloud_var_616 % fraction(jcol_var_646, jlev_var_645) >= 1D-06) THEN
              DO jg_var_647 = 1, ng_var_644
                od_cloud_new_var_641(jg_var_647) = od_scaling_var_640(jg_var_647, jlev_var_645) * od_cloud_var_620(config_var_614 % i_band_from_reordered_g_sw(jg_var_647), jlev_var_645, jcol_var_646)
                od_total_var_636(jg_var_647) = od_var_617(jg_var_647, jlev_var_645, jcol_var_646) + od_cloud_new_var_641(jg_var_647)
                ssa_total_var_637(jg_var_647) = 0.0D0
                g_total_var_638(jg_var_647) = 0.0D0
                IF (od_total_var_636(jg_var_647) > 0.0D0) THEN
                  scat_od_var_639 = ssa_var_618(jg_var_647, jlev_var_645, jcol_var_646) * od_var_617(jg_var_647, jlev_var_645, jcol_var_646) + ssa_cloud_var_621(config_var_614 % i_band_from_reordered_g_sw(jg_var_647), jlev_var_645, jcol_var_646) * od_cloud_new_var_641(jg_var_647)
                  ssa_total_var_637(jg_var_647) = scat_od_var_639 / od_total_var_636(jg_var_647)
                  IF (scat_od_var_639 > 0.0D0) THEN
                    g_total_var_638(jg_var_647) = (g_var_619(jg_var_647, jlev_var_645, jcol_var_646) * ssa_var_618(jg_var_647, jlev_var_645, jcol_var_646) * od_var_617(jg_var_647, jlev_var_645, jcol_var_646) + g_cloud_var_622(config_var_614 % i_band_from_reordered_g_sw(jg_var_647), jlev_var_645, jcol_var_646) * ssa_cloud_var_621(config_var_614 % i_band_from_reordered_g_sw(jg_var_647), jlev_var_645, jcol_var_646) * od_cloud_new_var_641(jg_var_647)) / scat_od_var_639
                  END IF
                END IF
              END DO
              CALL calc_ref_trans_sw(ng_var_644, cos_sza_var_625, od_total_var_636, ssa_total_var_637, g_total_var_638, reflectance_var_628(:, jlev_var_645), transmittance_var_629(:, jlev_var_645), ref_dir_var_630(:, jlev_var_645), trans_dir_diff_var_631(:, jlev_var_645), trans_dir_dir_var_632(:, jlev_var_645))
            ELSE
              DO jg_var_647 = 1, ng_var_644
                reflectance_var_628(jg_var_647, jlev_var_645) = ref_clear_var_626(jg_var_647, jlev_var_645)
                transmittance_var_629(jg_var_647, jlev_var_645) = trans_clear_var_627(jg_var_647, jlev_var_645)
                ref_dir_var_630(jg_var_647, jlev_var_645) = ref_dir_clear(jg_var_647, jlev_var_645)
                trans_dir_diff_var_631(jg_var_647, jlev_var_645) = trans_dir_diff_clear(jg_var_647, jlev_var_645)
                trans_dir_dir_var_632(jg_var_647, jlev_var_645) = trans_dir_dir_clear(jg_var_647, jlev_var_645)
              END DO
            END IF
          END DO
          CALL adding_ica_sw(ng_var_644, nlev_var_611, incoming_sw_var_623(:, jcol_var_646), albedo_diffuse(:, jcol_var_646), albedo_direct(:, jcol_var_646), cos_sza_var_625, reflectance_var_628, transmittance_var_629, ref_dir_var_630, trans_dir_diff_var_631, trans_dir_dir_var_632, flux_up_var_633, flux_dn_diffuse_var_634, flux_dn_direct_var_635, albedo_var_20 = tmp_work_albedo, source_var_21 = tmp_work_source, inv_denominator = tmp_work_inv_denominator)
          DO jlev_var_645 = 1, nlev_var_611 + 1
            sum_up_var_643 = 0.0D0
            sum_dn_diff = 0.0D0
            sum_dn_dir = 0.0D0
            DO jg_var_647 = 1, ng_var_644
              sum_up_var_643 = sum_up_var_643 + flux_up_var_633(jg_var_647, jlev_var_645)
              sum_dn_diff = sum_dn_diff + flux_dn_diffuse_var_634(jg_var_647, jlev_var_645)
              sum_dn_dir = sum_dn_dir + flux_dn_direct_var_635(jg_var_647, jlev_var_645)
            END DO
            flux_var_624 % sw_up(jcol_var_646, jlev_var_645) = sum_up_var_643
            flux_var_624 % sw_dn(jcol_var_646, jlev_var_645) = sum_dn_diff + sum_dn_dir
            flux_var_624 % sw_dn_direct(jcol_var_646, jlev_var_645) = sum_dn_dir
          END DO
          DO jlev_var_645 = 1, nlev_var_611 + 1
            flux_var_624 % sw_up(jcol_var_646, jlev_var_645) = total_cloud_cover_var_642 * flux_var_624 % sw_up(jcol_var_646, jlev_var_645) + (1.0D0 - total_cloud_cover_var_642) * flux_var_624 % sw_up_clear(jcol_var_646, jlev_var_645)
            flux_var_624 % sw_dn(jcol_var_646, jlev_var_645) = total_cloud_cover_var_642 * flux_var_624 % sw_dn(jcol_var_646, jlev_var_645) + (1.0D0 - total_cloud_cover_var_642) * flux_var_624 % sw_dn_clear(jcol_var_646, jlev_var_645)
            flux_var_624 % sw_dn_direct(jcol_var_646, jlev_var_645) = total_cloud_cover_var_642 * flux_var_624 % sw_dn_direct(jcol_var_646, jlev_var_645) + (1.0D0 - total_cloud_cover_var_642) * flux_var_624 % sw_dn_direct_clear(jcol_var_646, jlev_var_645)
          END DO
          DO jg_var_647 = 1, ng_var_644
            flux_var_624 % sw_dn_diffuse_surf_g(jg_var_647, jcol_var_646) = flux_dn_diffuse_var_634(jg_var_647, nlev_var_611 + 1)
            flux_var_624 % sw_dn_direct_surf_g(jg_var_647, jcol_var_646) = flux_dn_direct_var_635(jg_var_647, nlev_var_611 + 1)
            flux_var_624 % sw_dn_diffuse_surf_g(jg_var_647, jcol_var_646) = total_cloud_cover_var_642 * flux_var_624 % sw_dn_diffuse_surf_g(jg_var_647, jcol_var_646) + (1.0D0 - total_cloud_cover_var_642) * flux_var_624 % sw_dn_diffuse_surf_clear_g(jg_var_647, jcol_var_646)
            flux_var_624 % sw_dn_direct_surf_g(jg_var_647, jcol_var_646) = total_cloud_cover_var_642 * flux_var_624 % sw_dn_direct_surf_g(jg_var_647, jcol_var_646) + (1.0D0 - total_cloud_cover_var_642) * flux_var_624 % sw_dn_direct_surf_clear_g(jg_var_647, jcol_var_646)
          END DO
        ELSE
          DO jlev_var_645 = 1, nlev_var_611 + 1
            flux_var_624 % sw_up(jcol_var_646, jlev_var_645) = flux_var_624 % sw_up_clear(jcol_var_646, jlev_var_645)
            flux_var_624 % sw_dn(jcol_var_646, jlev_var_645) = flux_var_624 % sw_dn_clear(jcol_var_646, jlev_var_645)
            flux_var_624 % sw_dn_direct(jcol_var_646, jlev_var_645) = flux_var_624 % sw_dn_direct_clear(jcol_var_646, jlev_var_645)
          END DO
          DO jg_var_647 = 1, ng_var_644
            flux_var_624 % sw_dn_diffuse_surf_g(jg_var_647, jcol_var_646) = flux_var_624 % sw_dn_diffuse_surf_clear_g(jg_var_647, jcol_var_646)
            flux_var_624 % sw_dn_direct_surf_g(jg_var_647, jcol_var_646) = flux_var_624 % sw_dn_direct_surf_clear_g(jg_var_647, jcol_var_646)
          END DO
        END IF
      ELSE
        DO jlev_var_645 = 1, nlev_var_611 + 1
          flux_var_624 % sw_up(jcol_var_646, jlev_var_645) = 0.0D0
          flux_var_624 % sw_dn(jcol_var_646, jlev_var_645) = 0.0D0
          flux_var_624 % sw_dn_direct(jcol_var_646, jlev_var_645) = 0.0D0
          flux_var_624 % sw_up_clear(jcol_var_646, jlev_var_645) = 0.0D0
          flux_var_624 % sw_dn_clear(jcol_var_646, jlev_var_645) = 0.0D0
          flux_var_624 % sw_dn_direct_clear(jcol_var_646, jlev_var_645) = 0.0D0
        END DO
        DO jg_var_647 = 1, ng_var_644
          flux_var_624 % sw_dn_diffuse_surf_g(jg_var_647, jcol_var_646) = 0.0D0
          flux_var_624 % sw_dn_direct_surf_g(jg_var_647, jcol_var_646) = 0.0D0
          flux_var_624 % sw_dn_diffuse_surf_clear_g(jg_var_647, jcol_var_646) = 0.0D0
          flux_var_624 % sw_dn_direct_surf_clear_g(jg_var_647, jcol_var_646) = 0.0D0
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
  SUBROUTINE radiation(ncol, nlev, istartcol, iendcol, config, single_level, thermodynamics, gas, cloud, aerosol, flux)
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
    REAL(KIND = 8), DIMENSION(140, nlev, istartcol : iendcol) :: od_lw_var_648
    REAL(KIND = 8), DIMENSION(0, nlev, istartcol : iendcol) :: ssa_lw_var_649, g_lw_var_650
    REAL(KIND = 8), DIMENSION(16, nlev, istartcol : iendcol) :: od_lw_cloud_var_651
    REAL(KIND = 8), DIMENSION(0, nlev, istartcol : iendcol) :: ssa_lw_cloud_var_652, g_lw_cloud_var_653
    REAL(KIND = 8), DIMENSION(112, nlev, istartcol : iendcol) :: od_sw_var_654, ssa_sw_var_655, g_sw_var_656
    REAL(KIND = 8), DIMENSION(14, nlev, istartcol : iendcol) :: od_sw_cloud_var_657, ssa_sw_cloud_var_658, g_sw_cloud_var_659
    REAL(KIND = 8), DIMENSION(140, nlev + 1, istartcol : iendcol) :: planck_hl_var_660
    REAL(KIND = 8), DIMENSION(140, istartcol : iendcol) :: lw_emission_var_661
    REAL(KIND = 8), DIMENSION(140, istartcol : iendcol) :: lw_albedo_var_662
    REAL(KIND = 8), DIMENSION(112, istartcol : iendcol) :: sw_albedo_direct_var_663
    REAL(KIND = 8), DIMENSION(112, istartcol : iendcol) :: sw_albedo_diffuse_var_664
    REAL(KIND = 8), DIMENSION(112, istartcol : iendcol) :: incoming_sw_var_665
    CALL global_init_fn
    CALL get_albedos(single_level, istartcol, iendcol, config, sw_albedo_direct_var_663, sw_albedo_diffuse_var_664, lw_albedo_var_662)
    CALL gas_optics(ncol, nlev, istartcol, iendcol, config, single_level, thermodynamics, gas, od_lw_var_648, od_sw_var_654, ssa_sw_var_655, lw_albedo_var_524 = lw_albedo_var_662, planck_hl_var_528 = planck_hl_var_660, lw_emission_var_529 = lw_emission_var_661, incoming_sw_var_530 = incoming_sw_var_665)
    CALL crop_cloud_fraction(cloud, istartcol, iendcol, 1D-06, 1D-09)
    CALL cloud_optics_fn_500(nlev, istartcol, iendcol, config, thermodynamics, cloud, od_lw_cloud_var_651, ssa_lw_cloud_var_652, g_lw_cloud_var_653, od_sw_cloud_var_657, ssa_sw_cloud_var_658, g_sw_cloud_var_659)
    CALL add_aerosol_optics(nlev, istartcol, iendcol, config, thermodynamics, gas, aerosol, od_lw_var_648, ssa_lw_var_649, g_lw_var_650, od_sw_var_654, ssa_sw_var_655, g_sw_var_656)
    CALL solver_mcica_lw(nlev, istartcol, iendcol, config, single_level, cloud, od_lw_var_648, ssa_lw_var_649, g_lw_var_650, od_lw_cloud_var_651, ssa_lw_cloud_var_652, g_lw_cloud_var_653, planck_hl_var_660, lw_emission_var_661, lw_albedo_var_662, flux)
    CALL solver_mcica_sw(nlev, istartcol, iendcol, config, single_level, cloud, od_sw_var_654, ssa_sw_var_655, g_sw_var_656, od_sw_cloud_var_657, ssa_sw_cloud_var_658, g_sw_cloud_var_659, sw_albedo_direct_var_663, sw_albedo_diffuse_var_664, incoming_sw_var_665, flux)
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
SUBROUTINE srtm_setcoef(kidia_var_666, kfdia_var_667, klev_var_668, pavel_var_669, ptavel_var_670, pcoldry_var_671, pwkl_var_672, klaytrop_var_673, pcolch4_var_674, pcolco2_var_675, pcolh2o_var_676, pcolmol_var_677, pcolo2_var_678, pcolo3_var_679, pforfac_var_680, pforfrac_var_681, kindfor_var_682, pselffac_var_683, pselffrac_var_684, kindself_var_685, pfac00_var_686, pfac01_var_687, pfac10_var_688, pfac11_var_689, kjp_var_690, kjt_var_691, kjt1_var_692, prmu0_var_693)
  USE yoesrtwn, ONLY: preflog_var_336, tref_var_337
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_666, kfdia_var_667
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_668
  REAL(KIND = 8), INTENT(IN) :: pavel_var_669(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(IN) :: ptavel_var_670(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(IN) :: pcoldry_var_671(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(IN) :: pwkl_var_672(kidia_var_666 : kfdia_var_667, 35, klev_var_668)
  INTEGER(KIND = 4), INTENT(INOUT) :: klaytrop_var_673(kidia_var_666 : kfdia_var_667)
  REAL(KIND = 8), INTENT(INOUT) :: pcolch4_var_674(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pcolco2_var_675(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pcolh2o_var_676(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pcolmol_var_677(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pcolo2_var_678(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pcolo3_var_679(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pforfac_var_680(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pforfrac_var_681(kidia_var_666 : kfdia_var_667, klev_var_668)
  INTEGER(KIND = 4), INTENT(INOUT) :: kindfor_var_682(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pselffac_var_683(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pselffrac_var_684(kidia_var_666 : kfdia_var_667, klev_var_668)
  INTEGER(KIND = 4), INTENT(INOUT) :: kindself_var_685(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pfac00_var_686(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pfac01_var_687(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pfac10_var_688(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(INOUT) :: pfac11_var_689(kidia_var_666 : kfdia_var_667, klev_var_668)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjp_var_690(kidia_var_666 : kfdia_var_667, klev_var_668)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjt_var_691(kidia_var_666 : kfdia_var_667, klev_var_668)
  INTEGER(KIND = 4), INTENT(INOUT) :: kjt1_var_692(kidia_var_666 : kfdia_var_667, klev_var_668)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_693(kidia_var_666 : kfdia_var_667)
  INTEGER(KIND = 4) :: i_nlayers_var_694, jk_var_695, jl_var_696, jp1_var_697
  REAL(KIND = 8) :: z_stpfac_var_698, z_plog_var_699
  REAL(KIND = 8) :: z_fp_var_700, z_ft_var_701, z_ft1_var_702, z_water_var_703, z_scalefac_var_704
  REAL(KIND = 8) :: z_factor_var_705, z_co2reg_var_706, z_compfp_var_707
  z_stpfac_var_698 = 0.29220138203356366D0
  i_nlayers_var_694 = klev_var_668
  DO jk_var_695 = 1, klev_var_668
    DO jl_var_696 = kidia_var_666, kfdia_var_667
      pcolmol_var_677(jl_var_696, jk_var_695) = 0.0D0
    END DO
  END DO
  DO jl_var_696 = kidia_var_666, kfdia_var_667
    IF (prmu0_var_693(jl_var_696) > 0.0D0) THEN
      klaytrop_var_673(jl_var_696) = 0
    END IF
  END DO
  DO jk_var_695 = 1, i_nlayers_var_694
    DO jl_var_696 = kidia_var_666, kfdia_var_667
      IF (prmu0_var_693(jl_var_696) > 0.0D0) THEN
        z_plog_var_699 = LOG(pavel_var_669(jl_var_696, jk_var_695))
        kjp_var_690(jl_var_696, jk_var_695) = INT(36.0D0 - 5.0D0 * (z_plog_var_699 + 0.04D0))
        IF (kjp_var_690(jl_var_696, jk_var_695) < 1) THEN
          kjp_var_690(jl_var_696, jk_var_695) = 1
        ELSE IF (kjp_var_690(jl_var_696, jk_var_695) > 58) THEN
          kjp_var_690(jl_var_696, jk_var_695) = 58
        END IF
        jp1_var_697 = kjp_var_690(jl_var_696, jk_var_695) + 1
        z_fp_var_700 = 5.0 * (preflog_var_336(kjp_var_690(jl_var_696, jk_var_695)) - z_plog_var_699)
        kjt_var_691(jl_var_696, jk_var_695) = INT(3.0 + (ptavel_var_670(jl_var_696, jk_var_695) - tref_var_337(kjp_var_690(jl_var_696, jk_var_695))) / 15.0)
        IF (kjt_var_691(jl_var_696, jk_var_695) < 1) THEN
          kjt_var_691(jl_var_696, jk_var_695) = 1
        ELSE IF (kjt_var_691(jl_var_696, jk_var_695) > 4) THEN
          kjt_var_691(jl_var_696, jk_var_695) = 4
        END IF
        z_ft_var_701 = ((ptavel_var_670(jl_var_696, jk_var_695) - tref_var_337(kjp_var_690(jl_var_696, jk_var_695))) / 15.0) - REAL(kjt_var_691(jl_var_696, jk_var_695) - 3)
        kjt1_var_692(jl_var_696, jk_var_695) = INT(3.0 + (ptavel_var_670(jl_var_696, jk_var_695) - tref_var_337(jp1_var_697)) / 15.0)
        IF (kjt1_var_692(jl_var_696, jk_var_695) < 1) THEN
          kjt1_var_692(jl_var_696, jk_var_695) = 1
        ELSE IF (kjt1_var_692(jl_var_696, jk_var_695) > 4) THEN
          kjt1_var_692(jl_var_696, jk_var_695) = 4
        END IF
        z_ft1_var_702 = ((ptavel_var_670(jl_var_696, jk_var_695) - tref_var_337(jp1_var_697)) / 15.0) - REAL(kjt1_var_692(jl_var_696, jk_var_695) - 3)
        z_water_var_703 = pwkl_var_672(jl_var_696, 1, jk_var_695) / pcoldry_var_671(jl_var_696, jk_var_695)
        z_scalefac_var_704 = pavel_var_669(jl_var_696, jk_var_695) * z_stpfac_var_698 / ptavel_var_670(jl_var_696, jk_var_695)
        IF (z_plog_var_699 <= 4.56D0) GO TO 5300
        klaytrop_var_673(jl_var_696) = klaytrop_var_673(jl_var_696) + 1
        pforfac_var_680(jl_var_696, jk_var_695) = z_scalefac_var_704 / (1.0 + z_water_var_703)
        z_factor_var_705 = (332.0 - ptavel_var_670(jl_var_696, jk_var_695)) / 36.0
        kindfor_var_682(jl_var_696, jk_var_695) = MIN(2, MAX(1, INT(z_factor_var_705)))
        pforfrac_var_681(jl_var_696, jk_var_695) = z_factor_var_705 - REAL(kindfor_var_682(jl_var_696, jk_var_695))
        pselffac_var_683(jl_var_696, jk_var_695) = z_water_var_703 * pforfac_var_680(jl_var_696, jk_var_695)
        z_factor_var_705 = (ptavel_var_670(jl_var_696, jk_var_695) - 188.0) / 7.2
        kindself_var_685(jl_var_696, jk_var_695) = MIN(9, MAX(1, INT(z_factor_var_705) - 7))
        pselffrac_var_684(jl_var_696, jk_var_695) = z_factor_var_705 - REAL(kindself_var_685(jl_var_696, jk_var_695) + 7)
        pcolh2o_var_676(jl_var_696, jk_var_695) = 1E-20 * pwkl_var_672(jl_var_696, 1, jk_var_695)
        pcolco2_var_675(jl_var_696, jk_var_695) = 1E-20 * pwkl_var_672(jl_var_696, 2, jk_var_695)
        pcolo3_var_679(jl_var_696, jk_var_695) = 1E-20 * pwkl_var_672(jl_var_696, 3, jk_var_695)
        pcolch4_var_674(jl_var_696, jk_var_695) = 1E-20 * pwkl_var_672(jl_var_696, 6, jk_var_695)
        pcolo2_var_678(jl_var_696, jk_var_695) = 1E-20 * pwkl_var_672(jl_var_696, 7, jk_var_695)
        pcolmol_var_677(jl_var_696, jk_var_695) = 1E-20 * pcoldry_var_671(jl_var_696, jk_var_695) + pcolh2o_var_676(jl_var_696, jk_var_695)
        IF (pcolco2_var_675(jl_var_696, jk_var_695) == 0.0) pcolco2_var_675(jl_var_696, jk_var_695) = 1E-32 * pcoldry_var_671(jl_var_696, jk_var_695)
        IF (pcolch4_var_674(jl_var_696, jk_var_695) == 0.0) pcolch4_var_674(jl_var_696, jk_var_695) = 1E-32 * pcoldry_var_671(jl_var_696, jk_var_695)
        IF (pcolo2_var_678(jl_var_696, jk_var_695) == 0.0) pcolo2_var_678(jl_var_696, jk_var_695) = 1E-32 * pcoldry_var_671(jl_var_696, jk_var_695)
        z_co2reg_var_706 = 3.55E-24 * pcoldry_var_671(jl_var_696, jk_var_695)
        GO TO 5400
5300    CONTINUE
        pforfac_var_680(jl_var_696, jk_var_695) = z_scalefac_var_704 / (1.0 + z_water_var_703)
        z_factor_var_705 = (ptavel_var_670(jl_var_696, jk_var_695) - 188.0) / 36.0
        kindfor_var_682(jl_var_696, jk_var_695) = 3
        pforfrac_var_681(jl_var_696, jk_var_695) = z_factor_var_705 - 1.0
        pcolh2o_var_676(jl_var_696, jk_var_695) = 1E-20 * pwkl_var_672(jl_var_696, 1, jk_var_695)
        pcolco2_var_675(jl_var_696, jk_var_695) = 1E-20 * pwkl_var_672(jl_var_696, 2, jk_var_695)
        pcolo3_var_679(jl_var_696, jk_var_695) = 1E-20 * pwkl_var_672(jl_var_696, 3, jk_var_695)
        pcolch4_var_674(jl_var_696, jk_var_695) = 1E-20 * pwkl_var_672(jl_var_696, 6, jk_var_695)
        pcolo2_var_678(jl_var_696, jk_var_695) = 1E-20 * pwkl_var_672(jl_var_696, 7, jk_var_695)
        pcolmol_var_677(jl_var_696, jk_var_695) = 1E-20 * pcoldry_var_671(jl_var_696, jk_var_695) + pcolh2o_var_676(jl_var_696, jk_var_695)
        IF (pcolco2_var_675(jl_var_696, jk_var_695) == 0.0) pcolco2_var_675(jl_var_696, jk_var_695) = 1E-32 * pcoldry_var_671(jl_var_696, jk_var_695)
        IF (pcolch4_var_674(jl_var_696, jk_var_695) == 0.0) pcolch4_var_674(jl_var_696, jk_var_695) = 1E-32 * pcoldry_var_671(jl_var_696, jk_var_695)
        IF (pcolo2_var_678(jl_var_696, jk_var_695) == 0.0) pcolo2_var_678(jl_var_696, jk_var_695) = 1E-32 * pcoldry_var_671(jl_var_696, jk_var_695)
        z_co2reg_var_706 = 3.55E-24 * pcoldry_var_671(jl_var_696, jk_var_695)
        pselffac_var_683(jl_var_696, jk_var_695) = 0.0D0
        pselffrac_var_684(jl_var_696, jk_var_695) = 0.0D0
        kindself_var_685(jl_var_696, jk_var_695) = 0
5400    CONTINUE
        z_compfp_var_707 = 1.0 - z_fp_var_700
        pfac10_var_688(jl_var_696, jk_var_695) = z_compfp_var_707 * z_ft_var_701
        pfac00_var_686(jl_var_696, jk_var_695) = z_compfp_var_707 * (1.0 - z_ft_var_701)
        pfac11_var_689(jl_var_696, jk_var_695) = z_fp_var_700 * z_ft1_var_702
        pfac01_var_687(jl_var_696, jk_var_695) = z_fp_var_700 * (1.0 - z_ft1_var_702)
9000    FORMAT(1X, 2I3, 3I4, F6.1, 4F7.2, 12E9.2, 2I5)
      END IF
    END DO
  END DO
END SUBROUTINE srtm_setcoef
SUBROUTINE rrtm_gas_optical_depth(kidia_var_708, kfdia_var_709, klev_var_710, pod_var_711, pavel_var_712, pcoldry_var_713, pcolbrd, pwx_var_714, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, poneminus_var_722, pcolh2o_var_723, pcolco2_var_724, pcolo3_var_725, pcoln2o, pcolch4_var_726, pcolo2_var_727, p_co2mult_var_728, klaytrop_var_729, klayswtch, klaylow, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, kindminor, pscaleminor, pscaleminorn2, pminorfrac, prat_h2oco2_var_737, prat_h2oco2_1_var_738, prat_h2oo3_var_739, prat_h2oo3_1_var_740, prat_h2on2o_var_741, prat_h2on2o_1_var_742, prat_h2och4_var_743, prat_h2och4_1_var_744, prat_n2oco2_var_745, prat_n2oco2_1_var_746, prat_o3co2_var_747, prat_o3co2_1_var_748)
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_708
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_709
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_710
  REAL(KIND = 8), INTENT(OUT) :: pod_var_711(140, klev_var_710, kidia_var_708 : kfdia_var_709)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_712(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pcoldry_var_713(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pwx_var_714(kidia_var_708 : kfdia_var_709, 4, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: ptauaerl(kidia_var_708 : kfdia_var_709, klev_var_710, 16)
  REAL(KIND = 8), INTENT(IN) :: pfac00_var_715(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pfac01_var_716(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pfac10_var_717(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pfac11_var_718(kidia_var_708 : kfdia_var_709, klev_var_710)
  INTEGER(KIND = 4), INTENT(IN) :: kjp_var_719(kidia_var_708 : kfdia_var_709, klev_var_710)
  INTEGER(KIND = 4), INTENT(IN) :: kjt_var_720(kidia_var_708 : kfdia_var_709, klev_var_710)
  INTEGER(KIND = 4), INTENT(IN) :: kjt1_var_721(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: poneminus_var_722
  REAL(KIND = 8), INTENT(IN) :: pcolh2o_var_723(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pcolco2_var_724(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pcolo3_var_725(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pcoln2o(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pcolch4_var_726(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pcolo2_var_727(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: p_co2mult_var_728(kidia_var_708 : kfdia_var_709, klev_var_710)
  INTEGER(KIND = 4), INTENT(IN) :: klaytrop_var_729(kidia_var_708 : kfdia_var_709)
  INTEGER(KIND = 4), INTENT(IN) :: klayswtch(kidia_var_708 : kfdia_var_709)
  INTEGER(KIND = 4), INTENT(IN) :: klaylow(kidia_var_708 : kfdia_var_709)
  REAL(KIND = 8), INTENT(IN) :: pselffac_var_730(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pselffrac_var_731(kidia_var_708 : kfdia_var_709, klev_var_710)
  INTEGER(KIND = 4), INTENT(IN) :: kindself_var_732(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(OUT) :: pfrac_var_733(kidia_var_708 : kfdia_var_709, 140, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pforfac_var_734(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pforfrac_var_735(kidia_var_708 : kfdia_var_709, klev_var_710)
  INTEGER(KIND = 4), INTENT(IN) :: kindfor_var_736(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pminorfrac(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pscaleminor(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pscaleminorn2(kidia_var_708 : kfdia_var_709, klev_var_710)
  INTEGER(KIND = 4), INTENT(IN) :: kindminor(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: pcolbrd(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8), INTENT(IN) :: prat_h2oco2_var_737(kidia_var_708 : kfdia_var_709, klev_var_710), prat_h2oco2_1_var_738(kidia_var_708 : kfdia_var_709, klev_var_710), prat_h2oo3_var_739(kidia_var_708 : kfdia_var_709, klev_var_710), prat_h2oo3_1_var_740(kidia_var_708 : kfdia_var_709, klev_var_710), prat_h2on2o_var_741(kidia_var_708 : kfdia_var_709, klev_var_710), prat_h2on2o_1_var_742(kidia_var_708 : kfdia_var_709, klev_var_710), prat_h2och4_var_743(kidia_var_708 : kfdia_var_709, klev_var_710), prat_h2och4_1_var_744(kidia_var_708 : kfdia_var_709, klev_var_710), prat_n2oco2_var_745(kidia_var_708 : kfdia_var_709, klev_var_710), prat_n2oco2_1_var_746(kidia_var_708 : kfdia_var_709, klev_var_710), prat_o3co2_var_747(kidia_var_708 : kfdia_var_709, klev_var_710), prat_o3co2_1_var_748(kidia_var_708 : kfdia_var_709, klev_var_710)
  REAL(KIND = 8) :: ztau(kidia_var_708 : kfdia_var_709, 140, klev_var_710)
  INTEGER(KIND = 4) :: ji, jlev_var_749
  INTEGER(KIND = 4) :: jlon_var_750
  pfrac_var_733(:, :, :) = 0.0D0
  CALL rrtm_taumol1(kidia_var_708, kfdia_var_709, klev_var_710, ztau, pavel_var_712, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, pcolh2o_var_723, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, pminorfrac, kindminor, pscaleminorn2, pcolbrd)
  CALL rrtm_taumol2(kidia_var_708, kfdia_var_709, klev_var_710, ztau, pavel_var_712, pcoldry_var_713, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, pcolh2o_var_723, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733)
  CALL rrtm_taumol3(kidia_var_708, kfdia_var_709, klev_var_710, ztau, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, poneminus_var_722, pcolh2o_var_723, pcolco2_var_724, pcoln2o, pcoldry_var_713, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, prat_h2oco2_var_737, prat_h2oco2_1_var_738, pminorfrac, kindminor)
  CALL rrtm_taumol4(kidia_var_708, kfdia_var_709, klev_var_710, ztau, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, poneminus_var_722, pcolh2o_var_723, pcolco2_var_724, pcolo3_var_725, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, prat_h2oco2_var_737, prat_h2oco2_1_var_738, prat_o3co2_var_747, prat_o3co2_1_var_748)
  CALL rrtm_taumol5(kidia_var_708, kfdia_var_709, klev_var_710, ztau, pwx_var_714, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, poneminus_var_722, pcolh2o_var_723, pcolco2_var_724, pcolo3_var_725, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, prat_h2oco2_var_737, prat_h2oco2_1_var_738, prat_o3co2_var_747, prat_o3co2_1_var_748, pminorfrac, kindminor)
  CALL rrtm_taumol6(kidia_var_708, kfdia_var_709, klev_var_710, ztau, pwx_var_714, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, pcolh2o_var_723, pcolco2_var_724, pcoldry_var_713, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, pminorfrac, kindminor)
  CALL rrtm_taumol7(kidia_var_708, kfdia_var_709, klev_var_710, ztau, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, poneminus_var_722, pcolh2o_var_723, pcolo3_var_725, pcolco2_var_724, pcoldry_var_713, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, prat_h2oo3_var_739, prat_h2oo3_1_var_740, pminorfrac, kindminor)
  CALL rrtm_taumol8(kidia_var_708, kfdia_var_709, klev_var_710, ztau, pwx_var_714, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, pcolh2o_var_723, pcolo3_var_725, pcoln2o, pcolco2_var_724, pcoldry_var_713, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, pminorfrac, kindminor)
  CALL rrtm_taumol9(kidia_var_708, kfdia_var_709, klev_var_710, ztau, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, poneminus_var_722, pcolh2o_var_723, pcoln2o, pcolch4_var_726, pcoldry_var_713, klaytrop_var_729, klayswtch, klaylow, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, prat_h2och4_var_743, prat_h2och4_1_var_744, pminorfrac, kindminor)
  CALL rrtm_taumol10(kidia_var_708, kfdia_var_709, klev_var_710, ztau, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, pcolh2o_var_723, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733)
  CALL rrtm_taumol11(kidia_var_708, kfdia_var_709, klev_var_710, ztau, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, pcolh2o_var_723, pcolo2_var_727, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, pminorfrac, kindminor, pscaleminor)
  CALL rrtm_taumol12(kidia_var_708, kfdia_var_709, klev_var_710, ztau, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, poneminus_var_722, pcolh2o_var_723, pcolco2_var_724, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, prat_h2oco2_var_737, prat_h2oco2_1_var_738)
  CALL rrtm_taumol13(kidia_var_708, kfdia_var_709, klev_var_710, ztau, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, poneminus_var_722, pcolh2o_var_723, pcoln2o, pcolco2_var_724, pcolo3_var_725, pcoldry_var_713, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, prat_h2on2o_var_741, prat_h2on2o_1_var_742, pminorfrac, kindminor)
  CALL rrtm_taumol14(kidia_var_708, kfdia_var_709, klev_var_710, ztau, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, pcolco2_var_724, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733)
  CALL rrtm_taumol15(kidia_var_708, kfdia_var_709, klev_var_710, ztau, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, poneminus_var_722, pcolh2o_var_723, pcolco2_var_724, pcoln2o, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, prat_n2oco2_var_745, prat_n2oco2_1_var_746, pminorfrac, kindminor, pscaleminor, pcolbrd)
  CALL rrtm_taumol16(kidia_var_708, kfdia_var_709, klev_var_710, ztau, ptauaerl, pfac00_var_715, pfac01_var_716, pfac10_var_717, pfac11_var_718, pforfac_var_734, pforfrac_var_735, kindfor_var_736, kjp_var_719, kjt_var_720, kjt1_var_721, poneminus_var_722, pcolh2o_var_723, pcolch4_var_726, klaytrop_var_729, pselffac_var_730, pselffrac_var_731, kindself_var_732, pfrac_var_733, prat_h2och4_var_743, prat_h2och4_1_var_744)
  DO jlev_var_749 = 1, klev_var_710
    DO ji = 1, 140
      DO jlon_var_750 = kidia_var_708, kfdia_var_709
        pod_var_711(ji, jlev_var_749, jlon_var_750) = ztau(jlon_var_750, ji, jlev_var_749)
      END DO
    END DO
  END DO
END SUBROUTINE rrtm_gas_optical_depth
SUBROUTINE rrtm_taumol9(kidia_var_751, kfdia_var_752, klev_var_753, taug_var_754, p_tauaerl_var_755, fac00_var_756, fac01_var_757, fac10_var_758, fac11_var_759, forfac_var_778, forfrac_var_779, indfor_var_777, jp_var_760, jt_var_761, jt1_var_762, oneminus_var_763, colh2o_var_764, coln2o_var_765, colch4_var_766, coldry_var_767, laytrop_var_768, k_layswtch_var_769, k_laylow_var_770, selffac_var_771, selffrac_var_772, indself_var_773, fracs_var_774, rat_h2och4_var_775, rat_h2och4_1_var_776, minorfrac_var_780, indminor_var_781)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrtm, ONLY: ng9
  USE yoerrta9, ONLY: absa_var_229, absb_var_230, forref_var_234, fracrefa_var_227, fracrefb_var_228, ka_mn2o_var_231, kb_mn2o_var_232, selfref_var_233
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_751
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_752
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_753
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_754(kidia_var_751 : kfdia_var_752, 140, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_755(kidia_var_751 : kfdia_var_752, klev_var_753, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_756(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_757(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_758(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_759(kidia_var_751 : kfdia_var_752, klev_var_753)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_760(kidia_var_751 : kfdia_var_752, klev_var_753)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_761(kidia_var_751 : kfdia_var_752, klev_var_753)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_762(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_763
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_764(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_765(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: colch4_var_766(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_767(kidia_var_751 : kfdia_var_752, klev_var_753)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_768(kidia_var_751 : kfdia_var_752)
  INTEGER(KIND = 4), INTENT(IN) :: k_layswtch_var_769(kidia_var_751 : kfdia_var_752)
  INTEGER(KIND = 4), INTENT(IN) :: k_laylow_var_770(kidia_var_751 : kfdia_var_752)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_771(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_772(kidia_var_751 : kfdia_var_752, klev_var_753)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_773(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_774(kidia_var_751 : kfdia_var_752, 140, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_var_775(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_1_var_776(kidia_var_751 : kfdia_var_752, klev_var_753)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_777(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_778(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_779(kidia_var_751 : kfdia_var_752, klev_var_753)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_780(kidia_var_751 : kfdia_var_752, klev_var_753)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_781(kidia_var_751 : kfdia_var_752, klev_var_753)
  INTEGER(KIND = 4) :: ind0_var_782, ind1_var_783, inds_var_784, indf_var_785, indm_var_786
  INTEGER(KIND = 4) :: ig_var_787, js_var_788, lay_var_789, js1_var_790, jmn2o_var_791, jpl_var_792
  REAL(KIND = 8) :: speccomb_var_793, speccomb1_var_794, speccomb_mn2o_var_795, speccomb_planck_var_796
  REAL(KIND = 8) :: refrat_planck_a_var_797, refrat_m_a_var_798
  REAL(KIND = 8) :: fs_var_799, specmult_var_800, specparm_var_801, fs1_var_802, specmult1_var_803, specparm1_var_804, fmn2o_var_805, specmult_mn2o_var_806, specparm_mn2o_var_807, fpl_var_808, specmult_planck_var_809, specparm_planck_var_810
  REAL(KIND = 8) :: adjfac_var_811, adjcoln2o_var_812, ratn2o_var_813, chi_n2o_var_814
  REAL(KIND = 8) :: fac000_var_815, fac100_var_816, fac200_var_817, fac010_var_818, fac110_var_819, fac210_var_820, fac001_var_821, fac101_var_822, fac201_var_823, fac011_var_824, fac111_var_825, fac211_var_826
  REAL(KIND = 8) :: p_var_827, p4_var_828, fk0_var_829, fk1_var_830, fk2_var_831
  REAL(KIND = 8) :: taufor_var_832, tauself_var_833, n2om1_var_834, n2om2_var_835, absn2o_var_836, tau_major_var_837(12), tau_major1_var_838(12)
  INTEGER(KIND = 4) :: laytrop_min_var_839, laytrop_max_var_840
  INTEGER(KIND = 4) :: ixc_var_841(klev_var_753), ixlow_var_842(kfdia_var_752, klev_var_753), ixhigh_var_843(kfdia_var_752, klev_var_753)
  INTEGER(KIND = 4) :: ich_var_844, icl_var_845, ixc0_var_846, ixp_var_847, jc_var_848, jl_var_849
  laytrop_min_var_839 = MINVAL(laytrop_var_768)
  laytrop_max_var_840 = MAXVAL(laytrop_var_768)
  ixlow_var_842 = 0
  ixhigh_var_843 = 0
  ixc_var_841 = 0
  DO lay_var_789 = laytrop_min_var_839 + 1, laytrop_max_var_840
    icl_var_845 = 0
    ich_var_844 = 0
    DO jc_var_848 = kidia_var_751, kfdia_var_752
      IF (lay_var_789 <= laytrop_var_768(jc_var_848)) THEN
        icl_var_845 = icl_var_845 + 1
        ixlow_var_842(icl_var_845, lay_var_789) = jc_var_848
      ELSE
        ich_var_844 = ich_var_844 + 1
        ixhigh_var_843(ich_var_844, lay_var_789) = jc_var_848
      END IF
    END DO
    ixc_var_841(lay_var_789) = icl_var_845
  END DO
  refrat_planck_a_var_797 = chi_mls(1, 9) / chi_mls(6, 9)
  refrat_m_a_var_798 = chi_mls(1, 3) / chi_mls(6, 3)
  DO lay_var_789 = 1, laytrop_min_var_839
    DO jl_var_849 = kidia_var_751, kfdia_var_752
      speccomb_var_793 = colh2o_var_764(jl_var_849, lay_var_789) + rat_h2och4_var_775(jl_var_849, lay_var_789) * colch4_var_766(jl_var_849, lay_var_789)
      specparm_var_801 = MIN(colh2o_var_764(jl_var_849, lay_var_789) / speccomb_var_793, oneminus_var_763)
      specmult_var_800 = 8.0D0 * (specparm_var_801)
      js_var_788 = 1 + INT(specmult_var_800)
      fs_var_799 = ((specmult_var_800) - AINT((specmult_var_800)))
      speccomb1_var_794 = colh2o_var_764(jl_var_849, lay_var_789) + rat_h2och4_1_var_776(jl_var_849, lay_var_789) * colch4_var_766(jl_var_849, lay_var_789)
      specparm1_var_804 = MIN(colh2o_var_764(jl_var_849, lay_var_789) / speccomb1_var_794, oneminus_var_763)
      specmult1_var_803 = 8.0D0 * (specparm1_var_804)
      js1_var_790 = 1 + INT(specmult1_var_803)
      fs1_var_802 = ((specmult1_var_803) - AINT((specmult1_var_803)))
      speccomb_mn2o_var_795 = colh2o_var_764(jl_var_849, lay_var_789) + refrat_m_a_var_798 * colch4_var_766(jl_var_849, lay_var_789)
      specparm_mn2o_var_807 = MIN(colh2o_var_764(jl_var_849, lay_var_789) / speccomb_mn2o_var_795, oneminus_var_763)
      specmult_mn2o_var_806 = 8.0D0 * specparm_mn2o_var_807
      jmn2o_var_791 = 1 + INT(specmult_mn2o_var_806)
      fmn2o_var_805 = ((specmult_mn2o_var_806) - AINT((specmult_mn2o_var_806)))
      chi_n2o_var_814 = coln2o_var_765(jl_var_849, lay_var_789) / (coldry_var_767(jl_var_849, lay_var_789))
      ratn2o_var_813 = 1D+20 * chi_n2o_var_814 / chi_mls(4, jp_var_760(jl_var_849, lay_var_789) + 1)
      IF (ratn2o_var_813 .GT. 1.5D0) THEN
        adjfac_var_811 = 0.5D0 + (ratn2o_var_813 - 0.5D0) ** 0.65D0
        adjcoln2o_var_812 = adjfac_var_811 * chi_mls(4, jp_var_760(jl_var_849, lay_var_789) + 1) * coldry_var_767(jl_var_849, lay_var_789) * 1D-20
      ELSE
        adjcoln2o_var_812 = coln2o_var_765(jl_var_849, lay_var_789)
      END IF
      speccomb_planck_var_796 = colh2o_var_764(jl_var_849, lay_var_789) + refrat_planck_a_var_797 * colch4_var_766(jl_var_849, lay_var_789)
      specparm_planck_var_810 = MIN(colh2o_var_764(jl_var_849, lay_var_789) / speccomb_planck_var_796, oneminus_var_763)
      specmult_planck_var_809 = 8.0D0 * specparm_planck_var_810
      jpl_var_792 = 1 + INT(specmult_planck_var_809)
      fpl_var_808 = ((specmult_planck_var_809) - AINT((specmult_planck_var_809)))
      ind0_var_782 = ((jp_var_760(jl_var_849, lay_var_789) - 1) * 5 + (jt_var_761(jl_var_849, lay_var_789) - 1)) * nspa_var_237(9) + js_var_788
      ind1_var_783 = (jp_var_760(jl_var_849, lay_var_789) * 5 + (jt1_var_762(jl_var_849, lay_var_789) - 1)) * nspa_var_237(9) + js1_var_790
      inds_var_784 = indself_var_773(jl_var_849, lay_var_789)
      indf_var_785 = indfor_var_777(jl_var_849, lay_var_789)
      indm_var_786 = indminor_var_781(jl_var_849, lay_var_789)
      IF (specparm_var_801 .LT. 0.125D0) THEN
        p_var_827 = fs_var_799 - 1.0D0
        p4_var_828 = p_var_827 ** 4
        fk0_var_829 = p4_var_828
        fk1_var_830 = 1.0D0 - p_var_827 - 2.0D0 * p4_var_828
        fk2_var_831 = p_var_827 + p4_var_828
        fac000_var_815 = fk0_var_829 * fac00_var_756(jl_var_849, lay_var_789)
        fac100_var_816 = fk1_var_830 * fac00_var_756(jl_var_849, lay_var_789)
        fac200_var_817 = fk2_var_831 * fac00_var_756(jl_var_849, lay_var_789)
        fac010_var_818 = fk0_var_829 * fac10_var_758(jl_var_849, lay_var_789)
        fac110_var_819 = fk1_var_830 * fac10_var_758(jl_var_849, lay_var_789)
        fac210_var_820 = fk2_var_831 * fac10_var_758(jl_var_849, lay_var_789)
      ELSE IF (specparm_var_801 .GT. 0.875D0) THEN
        p_var_827 = - fs_var_799
        p4_var_828 = p_var_827 ** 4
        fk0_var_829 = p4_var_828
        fk1_var_830 = 1.0D0 - p_var_827 - 2.0D0 * p4_var_828
        fk2_var_831 = p_var_827 + p4_var_828
        fac000_var_815 = fk0_var_829 * fac00_var_756(jl_var_849, lay_var_789)
        fac100_var_816 = fk1_var_830 * fac00_var_756(jl_var_849, lay_var_789)
        fac200_var_817 = fk2_var_831 * fac00_var_756(jl_var_849, lay_var_789)
        fac010_var_818 = fk0_var_829 * fac10_var_758(jl_var_849, lay_var_789)
        fac110_var_819 = fk1_var_830 * fac10_var_758(jl_var_849, lay_var_789)
        fac210_var_820 = fk2_var_831 * fac10_var_758(jl_var_849, lay_var_789)
      ELSE
        fac000_var_815 = (1.0D0 - fs_var_799) * fac00_var_756(jl_var_849, lay_var_789)
        fac010_var_818 = (1.0D0 - fs_var_799) * fac10_var_758(jl_var_849, lay_var_789)
        fac100_var_816 = fs_var_799 * fac00_var_756(jl_var_849, lay_var_789)
        fac110_var_819 = fs_var_799 * fac10_var_758(jl_var_849, lay_var_789)
        fac200_var_817 = 0.0D0
        fac210_var_820 = 0.0D0
      END IF
      IF (specparm1_var_804 .LT. 0.125D0) THEN
        p_var_827 = fs1_var_802 - 1.0D0
        p4_var_828 = p_var_827 ** 4
        fk0_var_829 = p4_var_828
        fk1_var_830 = 1.0D0 - p_var_827 - 2.0D0 * p4_var_828
        fk2_var_831 = p_var_827 + p4_var_828
        fac001_var_821 = fk0_var_829 * fac01_var_757(jl_var_849, lay_var_789)
        fac101_var_822 = fk1_var_830 * fac01_var_757(jl_var_849, lay_var_789)
        fac201_var_823 = fk2_var_831 * fac01_var_757(jl_var_849, lay_var_789)
        fac011_var_824 = fk0_var_829 * fac11_var_759(jl_var_849, lay_var_789)
        fac111_var_825 = fk1_var_830 * fac11_var_759(jl_var_849, lay_var_789)
        fac211_var_826 = fk2_var_831 * fac11_var_759(jl_var_849, lay_var_789)
      ELSE IF (specparm1_var_804 .GT. 0.875D0) THEN
        p_var_827 = - fs1_var_802
        p4_var_828 = p_var_827 ** 4
        fk0_var_829 = p4_var_828
        fk1_var_830 = 1.0D0 - p_var_827 - 2.0D0 * p4_var_828
        fk2_var_831 = p_var_827 + p4_var_828
        fac001_var_821 = fk0_var_829 * fac01_var_757(jl_var_849, lay_var_789)
        fac101_var_822 = fk1_var_830 * fac01_var_757(jl_var_849, lay_var_789)
        fac201_var_823 = fk2_var_831 * fac01_var_757(jl_var_849, lay_var_789)
        fac011_var_824 = fk0_var_829 * fac11_var_759(jl_var_849, lay_var_789)
        fac111_var_825 = fk1_var_830 * fac11_var_759(jl_var_849, lay_var_789)
        fac211_var_826 = fk2_var_831 * fac11_var_759(jl_var_849, lay_var_789)
      ELSE
        fac001_var_821 = (1.0D0 - fs1_var_802) * fac01_var_757(jl_var_849, lay_var_789)
        fac011_var_824 = (1.0D0 - fs1_var_802) * fac11_var_759(jl_var_849, lay_var_789)
        fac101_var_822 = fs1_var_802 * fac01_var_757(jl_var_849, lay_var_789)
        fac111_var_825 = fs1_var_802 * fac11_var_759(jl_var_849, lay_var_789)
        fac201_var_823 = 0.0D0
        fac211_var_826 = 0.0D0
      END IF
      IF (specparm_var_801 .LT. 0.125D0) THEN
        tau_major_var_837(1 : ng9) = speccomb_var_793 * (fac000_var_815 * absa_var_229(ind0_var_782, 1 : 12) + fac100_var_816 * absa_var_229(ind0_var_782 + 1, 1 : 12) + fac200_var_817 * absa_var_229(ind0_var_782 + 2, 1 : 12) + fac010_var_818 * absa_var_229(ind0_var_782 + 9, 1 : 12) + fac110_var_819 * absa_var_229(ind0_var_782 + 10, 1 : 12) + fac210_var_820 * absa_var_229(ind0_var_782 + 11, 1 : 12))
      ELSE IF (specparm_var_801 .GT. 0.875D0) THEN
        tau_major_var_837(1 : ng9) = speccomb_var_793 * (fac200_var_817 * absa_var_229(ind0_var_782 - 1, 1 : 12) + fac100_var_816 * absa_var_229(ind0_var_782, 1 : 12) + fac000_var_815 * absa_var_229(ind0_var_782 + 1, 1 : 12) + fac210_var_820 * absa_var_229(ind0_var_782 + 8, 1 : 12) + fac110_var_819 * absa_var_229(ind0_var_782 + 9, 1 : 12) + fac010_var_818 * absa_var_229(ind0_var_782 + 10, 1 : 12))
      ELSE
        tau_major_var_837(1 : ng9) = speccomb_var_793 * (fac000_var_815 * absa_var_229(ind0_var_782, 1 : 12) + fac100_var_816 * absa_var_229(ind0_var_782 + 1, 1 : 12) + fac010_var_818 * absa_var_229(ind0_var_782 + 9, 1 : 12) + fac110_var_819 * absa_var_229(ind0_var_782 + 10, 1 : 12))
      END IF
      IF (specparm1_var_804 .LT. 0.125D0) THEN
        tau_major1_var_838(1 : ng9) = speccomb1_var_794 * (fac001_var_821 * absa_var_229(ind1_var_783, 1 : 12) + fac101_var_822 * absa_var_229(ind1_var_783 + 1, 1 : 12) + fac201_var_823 * absa_var_229(ind1_var_783 + 2, 1 : 12) + fac011_var_824 * absa_var_229(ind1_var_783 + 9, 1 : 12) + fac111_var_825 * absa_var_229(ind1_var_783 + 10, 1 : 12) + fac211_var_826 * absa_var_229(ind1_var_783 + 11, 1 : 12))
      ELSE IF (specparm1_var_804 .GT. 0.875D0) THEN
        tau_major1_var_838(1 : ng9) = speccomb1_var_794 * (fac201_var_823 * absa_var_229(ind1_var_783 - 1, 1 : 12) + fac101_var_822 * absa_var_229(ind1_var_783, 1 : 12) + fac001_var_821 * absa_var_229(ind1_var_783 + 1, 1 : 12) + fac211_var_826 * absa_var_229(ind1_var_783 + 8, 1 : 12) + fac111_var_825 * absa_var_229(ind1_var_783 + 9, 1 : 12) + fac011_var_824 * absa_var_229(ind1_var_783 + 10, 1 : 12))
      ELSE
        tau_major1_var_838(1 : ng9) = speccomb1_var_794 * (fac001_var_821 * absa_var_229(ind1_var_783, 1 : 12) + fac101_var_822 * absa_var_229(ind1_var_783 + 1, 1 : 12) + fac011_var_824 * absa_var_229(ind1_var_783 + 9, 1 : 12) + fac111_var_825 * absa_var_229(ind1_var_783 + 10, 1 : 12))
      END IF
      DO ig_var_787 = 1, 12
        tauself_var_833 = selffac_var_771(jl_var_849, lay_var_789) * (selfref_var_233(inds_var_784, ig_var_787) + selffrac_var_772(jl_var_849, lay_var_789) * (selfref_var_233(inds_var_784 + 1, ig_var_787) - selfref_var_233(inds_var_784, ig_var_787)))
        taufor_var_832 = forfac_var_778(jl_var_849, lay_var_789) * (forref_var_234(indf_var_785, ig_var_787) + forfrac_var_779(jl_var_849, lay_var_789) * (forref_var_234(indf_var_785 + 1, ig_var_787) - forref_var_234(indf_var_785, ig_var_787)))
        n2om1_var_834 = ka_mn2o_var_231(jmn2o_var_791, indm_var_786, ig_var_787) + fmn2o_var_805 * (ka_mn2o_var_231(jmn2o_var_791 + 1, indm_var_786, ig_var_787) - ka_mn2o_var_231(jmn2o_var_791, indm_var_786, ig_var_787))
        n2om2_var_835 = ka_mn2o_var_231(jmn2o_var_791, indm_var_786 + 1, ig_var_787) + fmn2o_var_805 * (ka_mn2o_var_231(jmn2o_var_791 + 1, indm_var_786 + 1, ig_var_787) - ka_mn2o_var_231(jmn2o_var_791, indm_var_786 + 1, ig_var_787))
        absn2o_var_836 = n2om1_var_834 + minorfrac_var_780(jl_var_849, lay_var_789) * (n2om2_var_835 - n2om1_var_834)
        taug_var_754(jl_var_849, 96 + ig_var_787, lay_var_789) = tau_major_var_837(ig_var_787) + tau_major1_var_838(ig_var_787) + tauself_var_833 + taufor_var_832 + adjcoln2o_var_812 * absn2o_var_836
        fracs_var_774(jl_var_849, 96 + ig_var_787, lay_var_789) = fracrefa_var_227(ig_var_787, jpl_var_792) + fpl_var_808 * (fracrefa_var_227(ig_var_787, jpl_var_792 + 1) - fracrefa_var_227(ig_var_787, jpl_var_792))
      END DO
    END DO
  END DO
  DO lay_var_789 = laytrop_max_var_840 + 1, klev_var_753
    DO jl_var_849 = kidia_var_751, kfdia_var_752
      chi_n2o_var_814 = coln2o_var_765(jl_var_849, lay_var_789) / (coldry_var_767(jl_var_849, lay_var_789))
      ratn2o_var_813 = 1D+20 * chi_n2o_var_814 / chi_mls(4, jp_var_760(jl_var_849, lay_var_789) + 1)
      IF (ratn2o_var_813 .GT. 1.5D0) THEN
        adjfac_var_811 = 0.5D0 + (ratn2o_var_813 - 0.5D0) ** 0.65D0
        adjcoln2o_var_812 = adjfac_var_811 * chi_mls(4, jp_var_760(jl_var_849, lay_var_789) + 1) * coldry_var_767(jl_var_849, lay_var_789) * 1D-20
      ELSE
        adjcoln2o_var_812 = coln2o_var_765(jl_var_849, lay_var_789)
      END IF
      ind0_var_782 = ((jp_var_760(jl_var_849, lay_var_789) - 13) * 5 + (jt_var_761(jl_var_849, lay_var_789) - 1)) * nspb_var_238(9) + 1
      ind1_var_783 = ((jp_var_760(jl_var_849, lay_var_789) - 12) * 5 + (jt1_var_762(jl_var_849, lay_var_789) - 1)) * nspb_var_238(9) + 1
      indm_var_786 = indminor_var_781(jl_var_849, lay_var_789)
      DO ig_var_787 = 1, 12
        absn2o_var_836 = kb_mn2o_var_232(indm_var_786, ig_var_787) + minorfrac_var_780(jl_var_849, lay_var_789) * (kb_mn2o_var_232(indm_var_786 + 1, ig_var_787) - kb_mn2o_var_232(indm_var_786, ig_var_787))
        taug_var_754(jl_var_849, 96 + ig_var_787, lay_var_789) = colch4_var_766(jl_var_849, lay_var_789) * (fac00_var_756(jl_var_849, lay_var_789) * absb_var_230(ind0_var_782, ig_var_787) + fac10_var_758(jl_var_849, lay_var_789) * absb_var_230(ind0_var_782 + 1, ig_var_787) + fac01_var_757(jl_var_849, lay_var_789) * absb_var_230(ind1_var_783, ig_var_787) + fac11_var_759(jl_var_849, lay_var_789) * absb_var_230(ind1_var_783 + 1, ig_var_787)) + adjcoln2o_var_812 * absn2o_var_836
        fracs_var_774(jl_var_849, 96 + ig_var_787, lay_var_789) = fracrefb_var_228(ig_var_787)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_840 /= laytrop_min_var_839) THEN
    DO lay_var_789 = laytrop_min_var_839 + 1, laytrop_max_var_840
      ixc0_var_846 = ixc_var_841(lay_var_789)
      DO ixp_var_847 = 1, ixc0_var_846
        jl_var_849 = ixlow_var_842(ixp_var_847, lay_var_789)
        speccomb_var_793 = colh2o_var_764(jl_var_849, lay_var_789) + rat_h2och4_var_775(jl_var_849, lay_var_789) * colch4_var_766(jl_var_849, lay_var_789)
        specparm_var_801 = MIN(colh2o_var_764(jl_var_849, lay_var_789) / speccomb_var_793, oneminus_var_763)
        specmult_var_800 = 8.0D0 * (specparm_var_801)
        js_var_788 = 1 + INT(specmult_var_800)
        fs_var_799 = ((specmult_var_800) - AINT((specmult_var_800)))
        speccomb1_var_794 = colh2o_var_764(jl_var_849, lay_var_789) + rat_h2och4_1_var_776(jl_var_849, lay_var_789) * colch4_var_766(jl_var_849, lay_var_789)
        specparm1_var_804 = MIN(colh2o_var_764(jl_var_849, lay_var_789) / speccomb1_var_794, oneminus_var_763)
        specmult1_var_803 = 8.0D0 * (specparm1_var_804)
        js1_var_790 = 1 + INT(specmult1_var_803)
        fs1_var_802 = ((specmult1_var_803) - AINT((specmult1_var_803)))
        speccomb_mn2o_var_795 = colh2o_var_764(jl_var_849, lay_var_789) + refrat_m_a_var_798 * colch4_var_766(jl_var_849, lay_var_789)
        specparm_mn2o_var_807 = MIN(colh2o_var_764(jl_var_849, lay_var_789) / speccomb_mn2o_var_795, oneminus_var_763)
        specmult_mn2o_var_806 = 8.0D0 * specparm_mn2o_var_807
        jmn2o_var_791 = 1 + INT(specmult_mn2o_var_806)
        fmn2o_var_805 = ((specmult_mn2o_var_806) - AINT((specmult_mn2o_var_806)))
        chi_n2o_var_814 = coln2o_var_765(jl_var_849, lay_var_789) / (coldry_var_767(jl_var_849, lay_var_789))
        ratn2o_var_813 = 1D+20 * chi_n2o_var_814 / chi_mls(4, jp_var_760(jl_var_849, lay_var_789) + 1)
        IF (ratn2o_var_813 .GT. 1.5D0) THEN
          adjfac_var_811 = 0.5D0 + (ratn2o_var_813 - 0.5D0) ** 0.65D0
          adjcoln2o_var_812 = adjfac_var_811 * chi_mls(4, jp_var_760(jl_var_849, lay_var_789) + 1) * coldry_var_767(jl_var_849, lay_var_789) * 1D-20
        ELSE
          adjcoln2o_var_812 = coln2o_var_765(jl_var_849, lay_var_789)
        END IF
        speccomb_planck_var_796 = colh2o_var_764(jl_var_849, lay_var_789) + refrat_planck_a_var_797 * colch4_var_766(jl_var_849, lay_var_789)
        specparm_planck_var_810 = MIN(colh2o_var_764(jl_var_849, lay_var_789) / speccomb_planck_var_796, oneminus_var_763)
        specmult_planck_var_809 = 8.0D0 * specparm_planck_var_810
        jpl_var_792 = 1 + INT(specmult_planck_var_809)
        fpl_var_808 = ((specmult_planck_var_809) - AINT((specmult_planck_var_809)))
        ind0_var_782 = ((jp_var_760(jl_var_849, lay_var_789) - 1) * 5 + (jt_var_761(jl_var_849, lay_var_789) - 1)) * nspa_var_237(9) + js_var_788
        ind1_var_783 = (jp_var_760(jl_var_849, lay_var_789) * 5 + (jt1_var_762(jl_var_849, lay_var_789) - 1)) * nspa_var_237(9) + js1_var_790
        inds_var_784 = indself_var_773(jl_var_849, lay_var_789)
        indf_var_785 = indfor_var_777(jl_var_849, lay_var_789)
        indm_var_786 = indminor_var_781(jl_var_849, lay_var_789)
        IF (specparm_var_801 .LT. 0.125D0) THEN
          p_var_827 = fs_var_799 - 1.0D0
          p4_var_828 = p_var_827 ** 4
          fk0_var_829 = p4_var_828
          fk1_var_830 = 1.0D0 - p_var_827 - 2.0D0 * p4_var_828
          fk2_var_831 = p_var_827 + p4_var_828
          fac000_var_815 = fk0_var_829 * fac00_var_756(jl_var_849, lay_var_789)
          fac100_var_816 = fk1_var_830 * fac00_var_756(jl_var_849, lay_var_789)
          fac200_var_817 = fk2_var_831 * fac00_var_756(jl_var_849, lay_var_789)
          fac010_var_818 = fk0_var_829 * fac10_var_758(jl_var_849, lay_var_789)
          fac110_var_819 = fk1_var_830 * fac10_var_758(jl_var_849, lay_var_789)
          fac210_var_820 = fk2_var_831 * fac10_var_758(jl_var_849, lay_var_789)
        ELSE IF (specparm_var_801 .GT. 0.875D0) THEN
          p_var_827 = - fs_var_799
          p4_var_828 = p_var_827 ** 4
          fk0_var_829 = p4_var_828
          fk1_var_830 = 1.0D0 - p_var_827 - 2.0D0 * p4_var_828
          fk2_var_831 = p_var_827 + p4_var_828
          fac000_var_815 = fk0_var_829 * fac00_var_756(jl_var_849, lay_var_789)
          fac100_var_816 = fk1_var_830 * fac00_var_756(jl_var_849, lay_var_789)
          fac200_var_817 = fk2_var_831 * fac00_var_756(jl_var_849, lay_var_789)
          fac010_var_818 = fk0_var_829 * fac10_var_758(jl_var_849, lay_var_789)
          fac110_var_819 = fk1_var_830 * fac10_var_758(jl_var_849, lay_var_789)
          fac210_var_820 = fk2_var_831 * fac10_var_758(jl_var_849, lay_var_789)
        ELSE
          fac000_var_815 = (1.0D0 - fs_var_799) * fac00_var_756(jl_var_849, lay_var_789)
          fac010_var_818 = (1.0D0 - fs_var_799) * fac10_var_758(jl_var_849, lay_var_789)
          fac100_var_816 = fs_var_799 * fac00_var_756(jl_var_849, lay_var_789)
          fac110_var_819 = fs_var_799 * fac10_var_758(jl_var_849, lay_var_789)
          fac200_var_817 = 0.0D0
          fac210_var_820 = 0.0D0
        END IF
        IF (specparm1_var_804 .LT. 0.125D0) THEN
          p_var_827 = fs1_var_802 - 1.0D0
          p4_var_828 = p_var_827 ** 4
          fk0_var_829 = p4_var_828
          fk1_var_830 = 1.0D0 - p_var_827 - 2.0D0 * p4_var_828
          fk2_var_831 = p_var_827 + p4_var_828
          fac001_var_821 = fk0_var_829 * fac01_var_757(jl_var_849, lay_var_789)
          fac101_var_822 = fk1_var_830 * fac01_var_757(jl_var_849, lay_var_789)
          fac201_var_823 = fk2_var_831 * fac01_var_757(jl_var_849, lay_var_789)
          fac011_var_824 = fk0_var_829 * fac11_var_759(jl_var_849, lay_var_789)
          fac111_var_825 = fk1_var_830 * fac11_var_759(jl_var_849, lay_var_789)
          fac211_var_826 = fk2_var_831 * fac11_var_759(jl_var_849, lay_var_789)
        ELSE IF (specparm1_var_804 .GT. 0.875D0) THEN
          p_var_827 = - fs1_var_802
          p4_var_828 = p_var_827 ** 4
          fk0_var_829 = p4_var_828
          fk1_var_830 = 1.0D0 - p_var_827 - 2.0D0 * p4_var_828
          fk2_var_831 = p_var_827 + p4_var_828
          fac001_var_821 = fk0_var_829 * fac01_var_757(jl_var_849, lay_var_789)
          fac101_var_822 = fk1_var_830 * fac01_var_757(jl_var_849, lay_var_789)
          fac201_var_823 = fk2_var_831 * fac01_var_757(jl_var_849, lay_var_789)
          fac011_var_824 = fk0_var_829 * fac11_var_759(jl_var_849, lay_var_789)
          fac111_var_825 = fk1_var_830 * fac11_var_759(jl_var_849, lay_var_789)
          fac211_var_826 = fk2_var_831 * fac11_var_759(jl_var_849, lay_var_789)
        ELSE
          fac001_var_821 = (1.0D0 - fs1_var_802) * fac01_var_757(jl_var_849, lay_var_789)
          fac011_var_824 = (1.0D0 - fs1_var_802) * fac11_var_759(jl_var_849, lay_var_789)
          fac101_var_822 = fs1_var_802 * fac01_var_757(jl_var_849, lay_var_789)
          fac111_var_825 = fs1_var_802 * fac11_var_759(jl_var_849, lay_var_789)
          fac201_var_823 = 0.0D0
          fac211_var_826 = 0.0D0
        END IF
        IF (specparm_var_801 .LT. 0.125D0) THEN
          tau_major_var_837(1 : ng9) = speccomb_var_793 * (fac000_var_815 * absa_var_229(ind0_var_782, 1 : 12) + fac100_var_816 * absa_var_229(ind0_var_782 + 1, 1 : 12) + fac200_var_817 * absa_var_229(ind0_var_782 + 2, 1 : 12) + fac010_var_818 * absa_var_229(ind0_var_782 + 9, 1 : 12) + fac110_var_819 * absa_var_229(ind0_var_782 + 10, 1 : 12) + fac210_var_820 * absa_var_229(ind0_var_782 + 11, 1 : 12))
        ELSE IF (specparm_var_801 .GT. 0.875D0) THEN
          tau_major_var_837(1 : ng9) = speccomb_var_793 * (fac200_var_817 * absa_var_229(ind0_var_782 - 1, 1 : 12) + fac100_var_816 * absa_var_229(ind0_var_782, 1 : 12) + fac000_var_815 * absa_var_229(ind0_var_782 + 1, 1 : 12) + fac210_var_820 * absa_var_229(ind0_var_782 + 8, 1 : 12) + fac110_var_819 * absa_var_229(ind0_var_782 + 9, 1 : 12) + fac010_var_818 * absa_var_229(ind0_var_782 + 10, 1 : 12))
        ELSE
          tau_major_var_837(1 : ng9) = speccomb_var_793 * (fac000_var_815 * absa_var_229(ind0_var_782, 1 : 12) + fac100_var_816 * absa_var_229(ind0_var_782 + 1, 1 : 12) + fac010_var_818 * absa_var_229(ind0_var_782 + 9, 1 : 12) + fac110_var_819 * absa_var_229(ind0_var_782 + 10, 1 : 12))
        END IF
        IF (specparm1_var_804 .LT. 0.125D0) THEN
          tau_major1_var_838(1 : ng9) = speccomb1_var_794 * (fac001_var_821 * absa_var_229(ind1_var_783, 1 : 12) + fac101_var_822 * absa_var_229(ind1_var_783 + 1, 1 : 12) + fac201_var_823 * absa_var_229(ind1_var_783 + 2, 1 : 12) + fac011_var_824 * absa_var_229(ind1_var_783 + 9, 1 : 12) + fac111_var_825 * absa_var_229(ind1_var_783 + 10, 1 : 12) + fac211_var_826 * absa_var_229(ind1_var_783 + 11, 1 : 12))
        ELSE IF (specparm1_var_804 .GT. 0.875D0) THEN
          tau_major1_var_838(1 : ng9) = speccomb1_var_794 * (fac201_var_823 * absa_var_229(ind1_var_783 - 1, 1 : 12) + fac101_var_822 * absa_var_229(ind1_var_783, 1 : 12) + fac001_var_821 * absa_var_229(ind1_var_783 + 1, 1 : 12) + fac211_var_826 * absa_var_229(ind1_var_783 + 8, 1 : 12) + fac111_var_825 * absa_var_229(ind1_var_783 + 9, 1 : 12) + fac011_var_824 * absa_var_229(ind1_var_783 + 10, 1 : 12))
        ELSE
          tau_major1_var_838(1 : ng9) = speccomb1_var_794 * (fac001_var_821 * absa_var_229(ind1_var_783, 1 : 12) + fac101_var_822 * absa_var_229(ind1_var_783 + 1, 1 : 12) + fac011_var_824 * absa_var_229(ind1_var_783 + 9, 1 : 12) + fac111_var_825 * absa_var_229(ind1_var_783 + 10, 1 : 12))
        END IF
        DO ig_var_787 = 1, 12
          tauself_var_833 = selffac_var_771(jl_var_849, lay_var_789) * (selfref_var_233(inds_var_784, ig_var_787) + selffrac_var_772(jl_var_849, lay_var_789) * (selfref_var_233(inds_var_784 + 1, ig_var_787) - selfref_var_233(inds_var_784, ig_var_787)))
          taufor_var_832 = forfac_var_778(jl_var_849, lay_var_789) * (forref_var_234(indf_var_785, ig_var_787) + forfrac_var_779(jl_var_849, lay_var_789) * (forref_var_234(indf_var_785 + 1, ig_var_787) - forref_var_234(indf_var_785, ig_var_787)))
          n2om1_var_834 = ka_mn2o_var_231(jmn2o_var_791, indm_var_786, ig_var_787) + fmn2o_var_805 * (ka_mn2o_var_231(jmn2o_var_791 + 1, indm_var_786, ig_var_787) - ka_mn2o_var_231(jmn2o_var_791, indm_var_786, ig_var_787))
          n2om2_var_835 = ka_mn2o_var_231(jmn2o_var_791, indm_var_786 + 1, ig_var_787) + fmn2o_var_805 * (ka_mn2o_var_231(jmn2o_var_791 + 1, indm_var_786 + 1, ig_var_787) - ka_mn2o_var_231(jmn2o_var_791, indm_var_786 + 1, ig_var_787))
          absn2o_var_836 = n2om1_var_834 + minorfrac_var_780(jl_var_849, lay_var_789) * (n2om2_var_835 - n2om1_var_834)
          taug_var_754(jl_var_849, 96 + ig_var_787, lay_var_789) = tau_major_var_837(ig_var_787) + tau_major1_var_838(ig_var_787) + tauself_var_833 + taufor_var_832 + adjcoln2o_var_812 * absn2o_var_836
          fracs_var_774(jl_var_849, 96 + ig_var_787, lay_var_789) = fracrefa_var_227(ig_var_787, jpl_var_792) + fpl_var_808 * (fracrefa_var_227(ig_var_787, jpl_var_792 + 1) - fracrefa_var_227(ig_var_787, jpl_var_792))
        END DO
      END DO
      ixc0_var_846 = kfdia_var_752 - kidia_var_751 + 1 - ixc0_var_846
      DO ixp_var_847 = 1, ixc0_var_846
        jl_var_849 = ixhigh_var_843(ixp_var_847, lay_var_789)
        chi_n2o_var_814 = coln2o_var_765(jl_var_849, lay_var_789) / (coldry_var_767(jl_var_849, lay_var_789))
        ratn2o_var_813 = 1D+20 * chi_n2o_var_814 / chi_mls(4, jp_var_760(jl_var_849, lay_var_789) + 1)
        IF (ratn2o_var_813 .GT. 1.5D0) THEN
          adjfac_var_811 = 0.5D0 + (ratn2o_var_813 - 0.5D0) ** 0.65D0
          adjcoln2o_var_812 = adjfac_var_811 * chi_mls(4, jp_var_760(jl_var_849, lay_var_789) + 1) * coldry_var_767(jl_var_849, lay_var_789) * 1D-20
        ELSE
          adjcoln2o_var_812 = coln2o_var_765(jl_var_849, lay_var_789)
        END IF
        ind0_var_782 = ((jp_var_760(jl_var_849, lay_var_789) - 13) * 5 + (jt_var_761(jl_var_849, lay_var_789) - 1)) * nspb_var_238(9) + 1
        ind1_var_783 = ((jp_var_760(jl_var_849, lay_var_789) - 12) * 5 + (jt1_var_762(jl_var_849, lay_var_789) - 1)) * nspb_var_238(9) + 1
        indm_var_786 = indminor_var_781(jl_var_849, lay_var_789)
        DO ig_var_787 = 1, 12
          absn2o_var_836 = kb_mn2o_var_232(indm_var_786, ig_var_787) + minorfrac_var_780(jl_var_849, lay_var_789) * (kb_mn2o_var_232(indm_var_786 + 1, ig_var_787) - kb_mn2o_var_232(indm_var_786, ig_var_787))
          taug_var_754(jl_var_849, 96 + ig_var_787, lay_var_789) = colch4_var_766(jl_var_849, lay_var_789) * (fac00_var_756(jl_var_849, lay_var_789) * absb_var_230(ind0_var_782, ig_var_787) + fac10_var_758(jl_var_849, lay_var_789) * absb_var_230(ind0_var_782 + 1, ig_var_787) + fac01_var_757(jl_var_849, lay_var_789) * absb_var_230(ind1_var_783, ig_var_787) + fac11_var_759(jl_var_849, lay_var_789) * absb_var_230(ind1_var_783 + 1, ig_var_787)) + adjcoln2o_var_812 * absn2o_var_836
          fracs_var_774(jl_var_849, 96 + ig_var_787, lay_var_789) = fracrefb_var_228(ig_var_787)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol9
SUBROUTINE rrtm_taumol11(kidia_var_850, kfdia_var_851, klev_var_852, taug_var_853, p_tauaerl_var_854, fac00_var_855, fac01_var_856, fac10_var_857, fac11_var_858, forfac_var_869, forfrac_var_870, indfor_var_868, jp_var_859, jt_var_860, jt1_var_861, colh2o_var_862, colo2, laytrop_var_863, selffac_var_864, selffrac_var_865, indself_var_866, fracs_var_867, minorfrac_var_871, indminor_var_872, scaleminor_var_873)
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrta11, ONLY: absa_var_143, absb_var_144, forref_var_146, fracrefa_var_141, fracrefb_var_142, ka_mo2, kb_mo2, selfref_var_145
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_850
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_851
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_852
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_853(kidia_var_850 : kfdia_var_851, 140, klev_var_852)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_854(kidia_var_850 : kfdia_var_851, klev_var_852, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_855(kidia_var_850 : kfdia_var_851, klev_var_852)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_856(kidia_var_850 : kfdia_var_851, klev_var_852)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_857(kidia_var_850 : kfdia_var_851, klev_var_852)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_858(kidia_var_850 : kfdia_var_851, klev_var_852)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_859(kidia_var_850 : kfdia_var_851, klev_var_852)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_860(kidia_var_850 : kfdia_var_851, klev_var_852)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_861(kidia_var_850 : kfdia_var_851, klev_var_852)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_862(kidia_var_850 : kfdia_var_851, klev_var_852)
  REAL(KIND = 8), INTENT(IN) :: colo2(kidia_var_850 : kfdia_var_851, klev_var_852)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_863(kidia_var_850 : kfdia_var_851)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_864(kidia_var_850 : kfdia_var_851, klev_var_852)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_865(kidia_var_850 : kfdia_var_851, klev_var_852)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_866(kidia_var_850 : kfdia_var_851, klev_var_852)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_867(kidia_var_850 : kfdia_var_851, 140, klev_var_852)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_868(kidia_var_850 : kfdia_var_851, klev_var_852)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_869(kidia_var_850 : kfdia_var_851, klev_var_852)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_870(kidia_var_850 : kfdia_var_851, klev_var_852)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_871(kidia_var_850 : kfdia_var_851, klev_var_852)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_872(kidia_var_850 : kfdia_var_851, klev_var_852)
  REAL(KIND = 8), INTENT(IN) :: scaleminor_var_873(kidia_var_850 : kfdia_var_851, klev_var_852)
  INTEGER(KIND = 4) :: ind0_var_874, ind1_var_875
  INTEGER(KIND = 4) :: inds_var_876, indf_var_877, indm_var_878
  INTEGER(KIND = 4) :: ig_var_879, lay_var_880
  REAL(KIND = 8) :: taufor_var_881, tauself_var_882, scaleo2, tauo2
  INTEGER(KIND = 4) :: laytrop_min_var_883, laytrop_max_var_884
  INTEGER(KIND = 4) :: ixc_var_885(klev_var_852), ixlow_var_886(kfdia_var_851, klev_var_852), ixhigh_var_887(kfdia_var_851, klev_var_852)
  INTEGER(KIND = 4) :: ich_var_888, icl_var_889, ixc0_var_890, ixp_var_891, jc_var_892, jl_var_893
  laytrop_min_var_883 = MINVAL(laytrop_var_863)
  laytrop_max_var_884 = MAXVAL(laytrop_var_863)
  ixlow_var_886 = 0
  ixhigh_var_887 = 0
  ixc_var_885 = 0
  DO lay_var_880 = laytrop_min_var_883 + 1, laytrop_max_var_884
    icl_var_889 = 0
    ich_var_888 = 0
    DO jc_var_892 = kidia_var_850, kfdia_var_851
      IF (lay_var_880 <= laytrop_var_863(jc_var_892)) THEN
        icl_var_889 = icl_var_889 + 1
        ixlow_var_886(icl_var_889, lay_var_880) = jc_var_892
      ELSE
        ich_var_888 = ich_var_888 + 1
        ixhigh_var_887(ich_var_888, lay_var_880) = jc_var_892
      END IF
    END DO
    ixc_var_885(lay_var_880) = icl_var_889
  END DO
  DO lay_var_880 = 1, laytrop_min_var_883
    DO jl_var_893 = kidia_var_850, kfdia_var_851
      ind0_var_874 = ((jp_var_859(jl_var_893, lay_var_880) - 1) * 5 + (jt_var_860(jl_var_893, lay_var_880) - 1)) * nspa_var_237(11) + 1
      ind1_var_875 = (jp_var_859(jl_var_893, lay_var_880) * 5 + (jt1_var_861(jl_var_893, lay_var_880) - 1)) * nspa_var_237(11) + 1
      inds_var_876 = indself_var_866(jl_var_893, lay_var_880)
      indf_var_877 = indfor_var_868(jl_var_893, lay_var_880)
      indm_var_878 = indminor_var_872(jl_var_893, lay_var_880)
      scaleo2 = colo2(jl_var_893, lay_var_880) * scaleminor_var_873(jl_var_893, lay_var_880)
      DO ig_var_879 = 1, 8
        tauself_var_882 = selffac_var_864(jl_var_893, lay_var_880) * (selfref_var_145(inds_var_876, ig_var_879) + selffrac_var_865(jl_var_893, lay_var_880) * (selfref_var_145(inds_var_876 + 1, ig_var_879) - selfref_var_145(inds_var_876, ig_var_879)))
        taufor_var_881 = forfac_var_869(jl_var_893, lay_var_880) * (forref_var_146(indf_var_877, ig_var_879) + forfrac_var_870(jl_var_893, lay_var_880) * (forref_var_146(indf_var_877 + 1, ig_var_879) - forref_var_146(indf_var_877, ig_var_879)))
        tauo2 = scaleo2 * (ka_mo2(indm_var_878, ig_var_879) + minorfrac_var_871(jl_var_893, lay_var_880) * (ka_mo2(indm_var_878 + 1, ig_var_879) - ka_mo2(indm_var_878, ig_var_879)))
        taug_var_853(jl_var_893, 114 + ig_var_879, lay_var_880) = colh2o_var_862(jl_var_893, lay_var_880) * (fac00_var_855(jl_var_893, lay_var_880) * absa_var_143(ind0_var_874, ig_var_879) + fac10_var_857(jl_var_893, lay_var_880) * absa_var_143(ind0_var_874 + 1, ig_var_879) + fac01_var_856(jl_var_893, lay_var_880) * absa_var_143(ind1_var_875, ig_var_879) + fac11_var_858(jl_var_893, lay_var_880) * absa_var_143(ind1_var_875 + 1, ig_var_879)) + tauself_var_882 + taufor_var_881 + tauo2
        fracs_var_867(jl_var_893, 114 + ig_var_879, lay_var_880) = fracrefa_var_141(ig_var_879)
      END DO
    END DO
  END DO
  DO lay_var_880 = laytrop_max_var_884 + 1, klev_var_852
    DO jl_var_893 = kidia_var_850, kfdia_var_851
      ind0_var_874 = ((jp_var_859(jl_var_893, lay_var_880) - 13) * 5 + (jt_var_860(jl_var_893, lay_var_880) - 1)) * nspb_var_238(11) + 1
      ind1_var_875 = ((jp_var_859(jl_var_893, lay_var_880) - 12) * 5 + (jt1_var_861(jl_var_893, lay_var_880) - 1)) * nspb_var_238(11) + 1
      indf_var_877 = indfor_var_868(jl_var_893, lay_var_880)
      indm_var_878 = indminor_var_872(jl_var_893, lay_var_880)
      scaleo2 = colo2(jl_var_893, lay_var_880) * scaleminor_var_873(jl_var_893, lay_var_880)
      DO ig_var_879 = 1, 8
        taufor_var_881 = forfac_var_869(jl_var_893, lay_var_880) * (forref_var_146(indf_var_877, ig_var_879) + forfrac_var_870(jl_var_893, lay_var_880) * (forref_var_146(indf_var_877 + 1, ig_var_879) - forref_var_146(indf_var_877, ig_var_879)))
        tauo2 = scaleo2 * (kb_mo2(indm_var_878, ig_var_879) + minorfrac_var_871(jl_var_893, lay_var_880) * (kb_mo2(indm_var_878 + 1, ig_var_879) - kb_mo2(indm_var_878, ig_var_879)))
        taug_var_853(jl_var_893, 114 + ig_var_879, lay_var_880) = colh2o_var_862(jl_var_893, lay_var_880) * (fac00_var_855(jl_var_893, lay_var_880) * absb_var_144(ind0_var_874, ig_var_879) + fac10_var_857(jl_var_893, lay_var_880) * absb_var_144(ind0_var_874 + 1, ig_var_879) + fac01_var_856(jl_var_893, lay_var_880) * absb_var_144(ind1_var_875, ig_var_879) + fac11_var_858(jl_var_893, lay_var_880) * absb_var_144(ind1_var_875 + 1, ig_var_879)) + taufor_var_881 + tauo2
        fracs_var_867(jl_var_893, 114 + ig_var_879, lay_var_880) = fracrefb_var_142(ig_var_879)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_884 /= laytrop_min_var_883) THEN
    DO lay_var_880 = laytrop_min_var_883 + 1, laytrop_max_var_884
      ixc0_var_890 = ixc_var_885(lay_var_880)
      DO ixp_var_891 = 1, ixc0_var_890
        jl_var_893 = ixlow_var_886(ixp_var_891, lay_var_880)
        ind0_var_874 = ((jp_var_859(jl_var_893, lay_var_880) - 1) * 5 + (jt_var_860(jl_var_893, lay_var_880) - 1)) * nspa_var_237(11) + 1
        ind1_var_875 = (jp_var_859(jl_var_893, lay_var_880) * 5 + (jt1_var_861(jl_var_893, lay_var_880) - 1)) * nspa_var_237(11) + 1
        inds_var_876 = indself_var_866(jl_var_893, lay_var_880)
        indf_var_877 = indfor_var_868(jl_var_893, lay_var_880)
        indm_var_878 = indminor_var_872(jl_var_893, lay_var_880)
        scaleo2 = colo2(jl_var_893, lay_var_880) * scaleminor_var_873(jl_var_893, lay_var_880)
        DO ig_var_879 = 1, 8
          tauself_var_882 = selffac_var_864(jl_var_893, lay_var_880) * (selfref_var_145(inds_var_876, ig_var_879) + selffrac_var_865(jl_var_893, lay_var_880) * (selfref_var_145(inds_var_876 + 1, ig_var_879) - selfref_var_145(inds_var_876, ig_var_879)))
          taufor_var_881 = forfac_var_869(jl_var_893, lay_var_880) * (forref_var_146(indf_var_877, ig_var_879) + forfrac_var_870(jl_var_893, lay_var_880) * (forref_var_146(indf_var_877 + 1, ig_var_879) - forref_var_146(indf_var_877, ig_var_879)))
          tauo2 = scaleo2 * (ka_mo2(indm_var_878, ig_var_879) + minorfrac_var_871(jl_var_893, lay_var_880) * (ka_mo2(indm_var_878 + 1, ig_var_879) - ka_mo2(indm_var_878, ig_var_879)))
          taug_var_853(jl_var_893, 114 + ig_var_879, lay_var_880) = colh2o_var_862(jl_var_893, lay_var_880) * (fac00_var_855(jl_var_893, lay_var_880) * absa_var_143(ind0_var_874, ig_var_879) + fac10_var_857(jl_var_893, lay_var_880) * absa_var_143(ind0_var_874 + 1, ig_var_879) + fac01_var_856(jl_var_893, lay_var_880) * absa_var_143(ind1_var_875, ig_var_879) + fac11_var_858(jl_var_893, lay_var_880) * absa_var_143(ind1_var_875 + 1, ig_var_879)) + tauself_var_882 + taufor_var_881 + tauo2
          fracs_var_867(jl_var_893, 114 + ig_var_879, lay_var_880) = fracrefa_var_141(ig_var_879)
        END DO
      END DO
      ixc0_var_890 = kfdia_var_851 - kidia_var_850 + 1 - ixc0_var_890
      DO ixp_var_891 = 1, ixc0_var_890
        jl_var_893 = ixhigh_var_887(ixp_var_891, lay_var_880)
        ind0_var_874 = ((jp_var_859(jl_var_893, lay_var_880) - 13) * 5 + (jt_var_860(jl_var_893, lay_var_880) - 1)) * nspb_var_238(11) + 1
        ind1_var_875 = ((jp_var_859(jl_var_893, lay_var_880) - 12) * 5 + (jt1_var_861(jl_var_893, lay_var_880) - 1)) * nspb_var_238(11) + 1
        indf_var_877 = indfor_var_868(jl_var_893, lay_var_880)
        indm_var_878 = indminor_var_872(jl_var_893, lay_var_880)
        scaleo2 = colo2(jl_var_893, lay_var_880) * scaleminor_var_873(jl_var_893, lay_var_880)
        DO ig_var_879 = 1, 8
          taufor_var_881 = forfac_var_869(jl_var_893, lay_var_880) * (forref_var_146(indf_var_877, ig_var_879) + forfrac_var_870(jl_var_893, lay_var_880) * (forref_var_146(indf_var_877 + 1, ig_var_879) - forref_var_146(indf_var_877, ig_var_879)))
          tauo2 = scaleo2 * (kb_mo2(indm_var_878, ig_var_879) + minorfrac_var_871(jl_var_893, lay_var_880) * (kb_mo2(indm_var_878 + 1, ig_var_879) - kb_mo2(indm_var_878, ig_var_879)))
          taug_var_853(jl_var_893, 114 + ig_var_879, lay_var_880) = colh2o_var_862(jl_var_893, lay_var_880) * (fac00_var_855(jl_var_893, lay_var_880) * absb_var_144(ind0_var_874, ig_var_879) + fac10_var_857(jl_var_893, lay_var_880) * absb_var_144(ind0_var_874 + 1, ig_var_879) + fac01_var_856(jl_var_893, lay_var_880) * absb_var_144(ind1_var_875, ig_var_879) + fac11_var_858(jl_var_893, lay_var_880) * absb_var_144(ind1_var_875 + 1, ig_var_879)) + taufor_var_881 + tauo2
          fracs_var_867(jl_var_893, 114 + ig_var_879, lay_var_880) = fracrefb_var_142(ig_var_879)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol11
SUBROUTINE rrtm_taumol10(kidia_var_894, kfdia_var_895, klev_var_896, taug_var_897, p_tauaerl_var_898, fac00_var_899, fac01_var_900, fac10_var_901, fac11_var_902, forfac_var_914, forfrac_var_913, indfor_var_912, jp_var_903, jt_var_904, jt1_var_905, colh2o_var_906, laytrop_var_907, selffac_var_909, selffrac_var_910, indself_var_911, fracs_var_908)
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrta10, ONLY: absa_var_137, absb_var_138, forref_var_140, fracrefa_var_135, fracrefb_var_136, selfref_var_139
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_894
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_895
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_896
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_897(kidia_var_894 : kfdia_var_895, 140, klev_var_896)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_898(kidia_var_894 : kfdia_var_895, klev_var_896, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_899(kidia_var_894 : kfdia_var_895, klev_var_896)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_900(kidia_var_894 : kfdia_var_895, klev_var_896)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_901(kidia_var_894 : kfdia_var_895, klev_var_896)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_902(kidia_var_894 : kfdia_var_895, klev_var_896)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_903(kidia_var_894 : kfdia_var_895, klev_var_896)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_904(kidia_var_894 : kfdia_var_895, klev_var_896)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_905(kidia_var_894 : kfdia_var_895, klev_var_896)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_906(kidia_var_894 : kfdia_var_895, klev_var_896)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_907(kidia_var_894 : kfdia_var_895)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_908(kidia_var_894 : kfdia_var_895, 140, klev_var_896)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_909(kidia_var_894 : kfdia_var_895, klev_var_896)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_910(kidia_var_894 : kfdia_var_895, klev_var_896)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_911(kidia_var_894 : kfdia_var_895, klev_var_896)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_912(kidia_var_894 : kfdia_var_895, klev_var_896)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_913(kidia_var_894 : kfdia_var_895, klev_var_896)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_914(kidia_var_894 : kfdia_var_895, klev_var_896)
  INTEGER(KIND = 4) :: ind0_var_915, ind1_var_916
  INTEGER(KIND = 4) :: inds_var_917, indf_var_918
  INTEGER(KIND = 4) :: ig_var_919, lay_var_920
  REAL(KIND = 8) :: taufor_var_921, tauself_var_922
  INTEGER(KIND = 4) :: laytrop_min_var_923, laytrop_max_var_924
  INTEGER(KIND = 4) :: ixc_var_925(klev_var_896), ixlow_var_926(kfdia_var_895, klev_var_896), ixhigh_var_927(kfdia_var_895, klev_var_896)
  INTEGER(KIND = 4) :: ich_var_928, icl_var_929, ixc0_var_930, ixp_var_931, jc_var_932, jl_var_933
  laytrop_min_var_923 = MINVAL(laytrop_var_907)
  laytrop_max_var_924 = MAXVAL(laytrop_var_907)
  ixlow_var_926 = 0
  ixhigh_var_927 = 0
  ixc_var_925 = 0
  DO lay_var_920 = laytrop_min_var_923 + 1, laytrop_max_var_924
    icl_var_929 = 0
    ich_var_928 = 0
    DO jc_var_932 = kidia_var_894, kfdia_var_895
      IF (lay_var_920 <= laytrop_var_907(jc_var_932)) THEN
        icl_var_929 = icl_var_929 + 1
        ixlow_var_926(icl_var_929, lay_var_920) = jc_var_932
      ELSE
        ich_var_928 = ich_var_928 + 1
        ixhigh_var_927(ich_var_928, lay_var_920) = jc_var_932
      END IF
    END DO
    ixc_var_925(lay_var_920) = icl_var_929
  END DO
  DO lay_var_920 = 1, laytrop_min_var_923
    DO jl_var_933 = kidia_var_894, kfdia_var_895
      ind0_var_915 = ((jp_var_903(jl_var_933, lay_var_920) - 1) * 5 + (jt_var_904(jl_var_933, lay_var_920) - 1)) * nspa_var_237(10) + 1
      ind1_var_916 = (jp_var_903(jl_var_933, lay_var_920) * 5 + (jt1_var_905(jl_var_933, lay_var_920) - 1)) * nspa_var_237(10) + 1
      inds_var_917 = indself_var_911(jl_var_933, lay_var_920)
      indf_var_918 = indfor_var_912(jl_var_933, lay_var_920)
      DO ig_var_919 = 1, 6
        tauself_var_922 = selffac_var_909(jl_var_933, lay_var_920) * (selfref_var_139(inds_var_917, ig_var_919) + selffrac_var_910(jl_var_933, lay_var_920) * (selfref_var_139(inds_var_917 + 1, ig_var_919) - selfref_var_139(inds_var_917, ig_var_919)))
        taufor_var_921 = forfac_var_914(jl_var_933, lay_var_920) * (forref_var_140(indf_var_918, ig_var_919) + forfrac_var_913(jl_var_933, lay_var_920) * (forref_var_140(indf_var_918 + 1, ig_var_919) - forref_var_140(indf_var_918, ig_var_919)))
        taug_var_897(jl_var_933, 108 + ig_var_919, lay_var_920) = colh2o_var_906(jl_var_933, lay_var_920) * (fac00_var_899(jl_var_933, lay_var_920) * absa_var_137(ind0_var_915, ig_var_919) + fac10_var_901(jl_var_933, lay_var_920) * absa_var_137(ind0_var_915 + 1, ig_var_919) + fac01_var_900(jl_var_933, lay_var_920) * absa_var_137(ind1_var_916, ig_var_919) + fac11_var_902(jl_var_933, lay_var_920) * absa_var_137(ind1_var_916 + 1, ig_var_919)) + tauself_var_922 + taufor_var_921
        fracs_var_908(jl_var_933, 108 + ig_var_919, lay_var_920) = fracrefa_var_135(ig_var_919)
      END DO
    END DO
  END DO
  DO lay_var_920 = laytrop_max_var_924 + 1, klev_var_896
    DO jl_var_933 = kidia_var_894, kfdia_var_895
      ind0_var_915 = ((jp_var_903(jl_var_933, lay_var_920) - 13) * 5 + (jt_var_904(jl_var_933, lay_var_920) - 1)) * nspb_var_238(10) + 1
      ind1_var_916 = ((jp_var_903(jl_var_933, lay_var_920) - 12) * 5 + (jt1_var_905(jl_var_933, lay_var_920) - 1)) * nspb_var_238(10) + 1
      indf_var_918 = indfor_var_912(jl_var_933, lay_var_920)
      DO ig_var_919 = 1, 6
        taufor_var_921 = forfac_var_914(jl_var_933, lay_var_920) * (forref_var_140(indf_var_918, ig_var_919) + forfrac_var_913(jl_var_933, lay_var_920) * (forref_var_140(indf_var_918 + 1, ig_var_919) - forref_var_140(indf_var_918, ig_var_919)))
        taug_var_897(jl_var_933, 108 + ig_var_919, lay_var_920) = colh2o_var_906(jl_var_933, lay_var_920) * (fac00_var_899(jl_var_933, lay_var_920) * absb_var_138(ind0_var_915, ig_var_919) + fac10_var_901(jl_var_933, lay_var_920) * absb_var_138(ind0_var_915 + 1, ig_var_919) + fac01_var_900(jl_var_933, lay_var_920) * absb_var_138(ind1_var_916, ig_var_919) + fac11_var_902(jl_var_933, lay_var_920) * absb_var_138(ind1_var_916 + 1, ig_var_919)) + taufor_var_921
        fracs_var_908(jl_var_933, 108 + ig_var_919, lay_var_920) = fracrefb_var_136(ig_var_919)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_924 /= laytrop_min_var_923) THEN
    DO lay_var_920 = laytrop_min_var_923 + 1, laytrop_max_var_924
      ixc0_var_930 = ixc_var_925(lay_var_920)
      DO ixp_var_931 = 1, ixc0_var_930
        jl_var_933 = ixlow_var_926(ixp_var_931, lay_var_920)
        ind0_var_915 = ((jp_var_903(jl_var_933, lay_var_920) - 1) * 5 + (jt_var_904(jl_var_933, lay_var_920) - 1)) * nspa_var_237(10) + 1
        ind1_var_916 = (jp_var_903(jl_var_933, lay_var_920) * 5 + (jt1_var_905(jl_var_933, lay_var_920) - 1)) * nspa_var_237(10) + 1
        inds_var_917 = indself_var_911(jl_var_933, lay_var_920)
        indf_var_918 = indfor_var_912(jl_var_933, lay_var_920)
        DO ig_var_919 = 1, 6
          tauself_var_922 = selffac_var_909(jl_var_933, lay_var_920) * (selfref_var_139(inds_var_917, ig_var_919) + selffrac_var_910(jl_var_933, lay_var_920) * (selfref_var_139(inds_var_917 + 1, ig_var_919) - selfref_var_139(inds_var_917, ig_var_919)))
          taufor_var_921 = forfac_var_914(jl_var_933, lay_var_920) * (forref_var_140(indf_var_918, ig_var_919) + forfrac_var_913(jl_var_933, lay_var_920) * (forref_var_140(indf_var_918 + 1, ig_var_919) - forref_var_140(indf_var_918, ig_var_919)))
          taug_var_897(jl_var_933, 108 + ig_var_919, lay_var_920) = colh2o_var_906(jl_var_933, lay_var_920) * (fac00_var_899(jl_var_933, lay_var_920) * absa_var_137(ind0_var_915, ig_var_919) + fac10_var_901(jl_var_933, lay_var_920) * absa_var_137(ind0_var_915 + 1, ig_var_919) + fac01_var_900(jl_var_933, lay_var_920) * absa_var_137(ind1_var_916, ig_var_919) + fac11_var_902(jl_var_933, lay_var_920) * absa_var_137(ind1_var_916 + 1, ig_var_919)) + tauself_var_922 + taufor_var_921
          fracs_var_908(jl_var_933, 108 + ig_var_919, lay_var_920) = fracrefa_var_135(ig_var_919)
        END DO
      END DO
      ixc0_var_930 = kfdia_var_895 - kidia_var_894 + 1 - ixc0_var_930
      DO ixp_var_931 = 1, ixc0_var_930
        jl_var_933 = ixhigh_var_927(ixp_var_931, lay_var_920)
        ind0_var_915 = ((jp_var_903(jl_var_933, lay_var_920) - 13) * 5 + (jt_var_904(jl_var_933, lay_var_920) - 1)) * nspb_var_238(10) + 1
        ind1_var_916 = ((jp_var_903(jl_var_933, lay_var_920) - 12) * 5 + (jt1_var_905(jl_var_933, lay_var_920) - 1)) * nspb_var_238(10) + 1
        indf_var_918 = indfor_var_912(jl_var_933, lay_var_920)
        DO ig_var_919 = 1, 6
          taufor_var_921 = forfac_var_914(jl_var_933, lay_var_920) * (forref_var_140(indf_var_918, ig_var_919) + forfrac_var_913(jl_var_933, lay_var_920) * (forref_var_140(indf_var_918 + 1, ig_var_919) - forref_var_140(indf_var_918, ig_var_919)))
          taug_var_897(jl_var_933, 108 + ig_var_919, lay_var_920) = colh2o_var_906(jl_var_933, lay_var_920) * (fac00_var_899(jl_var_933, lay_var_920) * absb_var_138(ind0_var_915, ig_var_919) + fac10_var_901(jl_var_933, lay_var_920) * absb_var_138(ind0_var_915 + 1, ig_var_919) + fac01_var_900(jl_var_933, lay_var_920) * absb_var_138(ind1_var_916, ig_var_919) + fac11_var_902(jl_var_933, lay_var_920) * absb_var_138(ind1_var_916 + 1, ig_var_919)) + taufor_var_921
          fracs_var_908(jl_var_933, 108 + ig_var_919, lay_var_920) = fracrefb_var_136(ig_var_919)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol10
SUBROUTINE rrtm_taumol8(kidia_var_934, kfdia_var_935, klev_var_936, taug_var_937, wx_var_938, p_tauaerl_var_939, fac00_var_940, fac01_var_941, fac10_var_942, fac11_var_943, forfac_var_959, forfrac_var_958, indfor_var_957, jp_var_944, jt_var_945, jt1_var_946, colh2o_var_947, colo3_var_948, coln2o_var_949, colco2_var_950, coldry_var_951, laytrop_var_952, selffac_var_953, selffrac_var_954, indself_var_955, fracs_var_956, minorfrac_var_960, indminor_var_961)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrta8, ONLY: absa_var_218, absb_var_219, cfc12_var_217, cfc22adj, forref_var_226, fracrefa_var_215, fracrefb_var_216, ka_mco2_var_220, ka_mn2o_var_221, ka_mo3_var_222, kb_mco2_var_223, kb_mn2o_var_224, selfref_var_225
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_934
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_935
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_936
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_937(kidia_var_934 : kfdia_var_935, 140, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: wx_var_938(kidia_var_934 : kfdia_var_935, 4, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_939(kidia_var_934 : kfdia_var_935, klev_var_936, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_940(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_941(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_942(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_943(kidia_var_934 : kfdia_var_935, klev_var_936)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_944(kidia_var_934 : kfdia_var_935, klev_var_936)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_945(kidia_var_934 : kfdia_var_935, klev_var_936)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_946(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_947(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_948(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_949(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_950(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_951(kidia_var_934 : kfdia_var_935, klev_var_936)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_952(kidia_var_934 : kfdia_var_935)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_953(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_954(kidia_var_934 : kfdia_var_935, klev_var_936)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_955(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_956(kidia_var_934 : kfdia_var_935, 140, klev_var_936)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_957(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_958(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_959(kidia_var_934 : kfdia_var_935, klev_var_936)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_960(kidia_var_934 : kfdia_var_935, klev_var_936)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_961(kidia_var_934 : kfdia_var_935, klev_var_936)
  INTEGER(KIND = 4) :: ind0_var_962, ind1_var_963, inds_var_964, indf_var_965, indm_var_966
  INTEGER(KIND = 4) :: ig_var_967, lay_var_968
  REAL(KIND = 8) :: chi_co2_var_969, ratco2_var_970, adjfac_var_971, adjcolco2_var_972
  REAL(KIND = 8) :: taufor_var_973, tauself_var_974, abso3_var_975, absco2_var_976, absn2o_var_977
  INTEGER(KIND = 4) :: laytrop_min_var_978, laytrop_max_var_979
  INTEGER(KIND = 4) :: ixc_var_980(klev_var_936), ixlow_var_981(kfdia_var_935, klev_var_936), ixhigh_var_982(kfdia_var_935, klev_var_936)
  INTEGER(KIND = 4) :: ich_var_983, icl_var_984, ixc0_var_985, ixp_var_986, jc_var_987, jl_var_988
  laytrop_min_var_978 = MINVAL(laytrop_var_952)
  laytrop_max_var_979 = MAXVAL(laytrop_var_952)
  ixlow_var_981 = 0
  ixhigh_var_982 = 0
  ixc_var_980 = 0
  DO lay_var_968 = laytrop_min_var_978 + 1, laytrop_max_var_979
    icl_var_984 = 0
    ich_var_983 = 0
    DO jc_var_987 = kidia_var_934, kfdia_var_935
      IF (lay_var_968 <= laytrop_var_952(jc_var_987)) THEN
        icl_var_984 = icl_var_984 + 1
        ixlow_var_981(icl_var_984, lay_var_968) = jc_var_987
      ELSE
        ich_var_983 = ich_var_983 + 1
        ixhigh_var_982(ich_var_983, lay_var_968) = jc_var_987
      END IF
    END DO
    ixc_var_980(lay_var_968) = icl_var_984
  END DO
  DO lay_var_968 = 1, laytrop_min_var_978
    DO jl_var_988 = kidia_var_934, kfdia_var_935
      chi_co2_var_969 = colco2_var_950(jl_var_988, lay_var_968) / (coldry_var_951(jl_var_988, lay_var_968))
      ratco2_var_970 = 1D+20 * chi_co2_var_969 / chi_mls(2, jp_var_944(jl_var_988, lay_var_968) + 1)
      IF (ratco2_var_970 .GT. 3.0D0) THEN
        adjfac_var_971 = 2.0D0 + (ratco2_var_970 - 2.0D0) ** 0.65D0
        adjcolco2_var_972 = adjfac_var_971 * chi_mls(2, jp_var_944(jl_var_988, lay_var_968) + 1) * coldry_var_951(jl_var_988, lay_var_968) * 1D-20
      ELSE
        adjcolco2_var_972 = colco2_var_950(jl_var_988, lay_var_968)
      END IF
      ind0_var_962 = ((jp_var_944(jl_var_988, lay_var_968) - 1) * 5 + (jt_var_945(jl_var_988, lay_var_968) - 1)) * nspa_var_237(8) + 1
      ind1_var_963 = (jp_var_944(jl_var_988, lay_var_968) * 5 + (jt1_var_946(jl_var_988, lay_var_968) - 1)) * nspa_var_237(8) + 1
      inds_var_964 = indself_var_955(jl_var_988, lay_var_968)
      indf_var_965 = indfor_var_957(jl_var_988, lay_var_968)
      indm_var_966 = indminor_var_961(jl_var_988, lay_var_968)
      DO ig_var_967 = 1, 8
        tauself_var_974 = selffac_var_953(jl_var_988, lay_var_968) * (selfref_var_225(inds_var_964, ig_var_967) + selffrac_var_954(jl_var_988, lay_var_968) * (selfref_var_225(inds_var_964 + 1, ig_var_967) - selfref_var_225(inds_var_964, ig_var_967)))
        taufor_var_973 = forfac_var_959(jl_var_988, lay_var_968) * (forref_var_226(indf_var_965, ig_var_967) + forfrac_var_958(jl_var_988, lay_var_968) * (forref_var_226(indf_var_965 + 1, ig_var_967) - forref_var_226(indf_var_965, ig_var_967)))
        absco2_var_976 = (ka_mco2_var_220(indm_var_966, ig_var_967) + minorfrac_var_960(jl_var_988, lay_var_968) * (ka_mco2_var_220(indm_var_966 + 1, ig_var_967) - ka_mco2_var_220(indm_var_966, ig_var_967)))
        abso3_var_975 = (ka_mo3_var_222(indm_var_966, ig_var_967) + minorfrac_var_960(jl_var_988, lay_var_968) * (ka_mo3_var_222(indm_var_966 + 1, ig_var_967) - ka_mo3_var_222(indm_var_966, ig_var_967)))
        absn2o_var_977 = (ka_mn2o_var_221(indm_var_966, ig_var_967) + minorfrac_var_960(jl_var_988, lay_var_968) * (ka_mn2o_var_221(indm_var_966 + 1, ig_var_967) - ka_mn2o_var_221(indm_var_966, ig_var_967)))
        taug_var_937(jl_var_988, 88 + ig_var_967, lay_var_968) = colh2o_var_947(jl_var_988, lay_var_968) * (fac00_var_940(jl_var_988, lay_var_968) * absa_var_218(ind0_var_962, ig_var_967) + fac10_var_942(jl_var_988, lay_var_968) * absa_var_218(ind0_var_962 + 1, ig_var_967) + fac01_var_941(jl_var_988, lay_var_968) * absa_var_218(ind1_var_963, ig_var_967) + fac11_var_943(jl_var_988, lay_var_968) * absa_var_218(ind1_var_963 + 1, ig_var_967)) + tauself_var_974 + taufor_var_973 + adjcolco2_var_972 * absco2_var_976 + colo3_var_948(jl_var_988, lay_var_968) * abso3_var_975 + coln2o_var_949(jl_var_988, lay_var_968) * absn2o_var_977 + wx_var_938(jl_var_988, 3, lay_var_968) * cfc12_var_217(ig_var_967) + wx_var_938(jl_var_988, 4, lay_var_968) * cfc22adj(ig_var_967)
        fracs_var_956(jl_var_988, 88 + ig_var_967, lay_var_968) = fracrefa_var_215(ig_var_967)
      END DO
    END DO
  END DO
  DO lay_var_968 = laytrop_max_var_979 + 1, klev_var_936
    DO jl_var_988 = kidia_var_934, kfdia_var_935
      chi_co2_var_969 = colco2_var_950(jl_var_988, lay_var_968) / coldry_var_951(jl_var_988, lay_var_968)
      ratco2_var_970 = 1D+20 * chi_co2_var_969 / chi_mls(2, jp_var_944(jl_var_988, lay_var_968) + 1)
      IF (ratco2_var_970 .GT. 3.0D0) THEN
        adjfac_var_971 = 2.0D0 + (ratco2_var_970 - 2.0D0) ** 0.65D0
        adjcolco2_var_972 = adjfac_var_971 * chi_mls(2, jp_var_944(jl_var_988, lay_var_968) + 1) * coldry_var_951(jl_var_988, lay_var_968) * 1D-20
      ELSE
        adjcolco2_var_972 = colco2_var_950(jl_var_988, lay_var_968)
      END IF
      ind0_var_962 = ((jp_var_944(jl_var_988, lay_var_968) - 13) * 5 + (jt_var_945(jl_var_988, lay_var_968) - 1)) * nspb_var_238(8) + 1
      ind1_var_963 = ((jp_var_944(jl_var_988, lay_var_968) - 12) * 5 + (jt1_var_946(jl_var_988, lay_var_968) - 1)) * nspb_var_238(8) + 1
      indm_var_966 = indminor_var_961(jl_var_988, lay_var_968)
      DO ig_var_967 = 1, 8
        absco2_var_976 = (kb_mco2_var_223(indm_var_966, ig_var_967) + minorfrac_var_960(jl_var_988, lay_var_968) * (kb_mco2_var_223(indm_var_966 + 1, ig_var_967) - kb_mco2_var_223(indm_var_966, ig_var_967)))
        absn2o_var_977 = (kb_mn2o_var_224(indm_var_966, ig_var_967) + minorfrac_var_960(jl_var_988, lay_var_968) * (kb_mn2o_var_224(indm_var_966 + 1, ig_var_967) - kb_mn2o_var_224(indm_var_966, ig_var_967)))
        taug_var_937(jl_var_988, 88 + ig_var_967, lay_var_968) = colo3_var_948(jl_var_988, lay_var_968) * (fac00_var_940(jl_var_988, lay_var_968) * absb_var_219(ind0_var_962, ig_var_967) + fac10_var_942(jl_var_988, lay_var_968) * absb_var_219(ind0_var_962 + 1, ig_var_967) + fac01_var_941(jl_var_988, lay_var_968) * absb_var_219(ind1_var_963, ig_var_967) + fac11_var_943(jl_var_988, lay_var_968) * absb_var_219(ind1_var_963 + 1, ig_var_967)) + adjcolco2_var_972 * absco2_var_976 + coln2o_var_949(jl_var_988, lay_var_968) * absn2o_var_977 + wx_var_938(jl_var_988, 3, lay_var_968) * cfc12_var_217(ig_var_967) + wx_var_938(jl_var_988, 4, lay_var_968) * cfc22adj(ig_var_967)
        fracs_var_956(jl_var_988, 88 + ig_var_967, lay_var_968) = fracrefb_var_216(ig_var_967)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_979 /= laytrop_min_var_978) THEN
    DO lay_var_968 = laytrop_min_var_978 + 1, laytrop_max_var_979
      ixc0_var_985 = ixc_var_980(lay_var_968)
      DO ixp_var_986 = 1, ixc0_var_985
        jl_var_988 = ixlow_var_981(ixp_var_986, lay_var_968)
        chi_co2_var_969 = colco2_var_950(jl_var_988, lay_var_968) / (coldry_var_951(jl_var_988, lay_var_968))
        ratco2_var_970 = 1D+20 * chi_co2_var_969 / chi_mls(2, jp_var_944(jl_var_988, lay_var_968) + 1)
        IF (ratco2_var_970 .GT. 3.0D0) THEN
          adjfac_var_971 = 2.0D0 + (ratco2_var_970 - 2.0D0) ** 0.65D0
          adjcolco2_var_972 = adjfac_var_971 * chi_mls(2, jp_var_944(jl_var_988, lay_var_968) + 1) * coldry_var_951(jl_var_988, lay_var_968) * 1D-20
        ELSE
          adjcolco2_var_972 = colco2_var_950(jl_var_988, lay_var_968)
        END IF
        ind0_var_962 = ((jp_var_944(jl_var_988, lay_var_968) - 1) * 5 + (jt_var_945(jl_var_988, lay_var_968) - 1)) * nspa_var_237(8) + 1
        ind1_var_963 = (jp_var_944(jl_var_988, lay_var_968) * 5 + (jt1_var_946(jl_var_988, lay_var_968) - 1)) * nspa_var_237(8) + 1
        inds_var_964 = indself_var_955(jl_var_988, lay_var_968)
        indf_var_965 = indfor_var_957(jl_var_988, lay_var_968)
        indm_var_966 = indminor_var_961(jl_var_988, lay_var_968)
        DO ig_var_967 = 1, 8
          tauself_var_974 = selffac_var_953(jl_var_988, lay_var_968) * (selfref_var_225(inds_var_964, ig_var_967) + selffrac_var_954(jl_var_988, lay_var_968) * (selfref_var_225(inds_var_964 + 1, ig_var_967) - selfref_var_225(inds_var_964, ig_var_967)))
          taufor_var_973 = forfac_var_959(jl_var_988, lay_var_968) * (forref_var_226(indf_var_965, ig_var_967) + forfrac_var_958(jl_var_988, lay_var_968) * (forref_var_226(indf_var_965 + 1, ig_var_967) - forref_var_226(indf_var_965, ig_var_967)))
          absco2_var_976 = (ka_mco2_var_220(indm_var_966, ig_var_967) + minorfrac_var_960(jl_var_988, lay_var_968) * (ka_mco2_var_220(indm_var_966 + 1, ig_var_967) - ka_mco2_var_220(indm_var_966, ig_var_967)))
          abso3_var_975 = (ka_mo3_var_222(indm_var_966, ig_var_967) + minorfrac_var_960(jl_var_988, lay_var_968) * (ka_mo3_var_222(indm_var_966 + 1, ig_var_967) - ka_mo3_var_222(indm_var_966, ig_var_967)))
          absn2o_var_977 = (ka_mn2o_var_221(indm_var_966, ig_var_967) + minorfrac_var_960(jl_var_988, lay_var_968) * (ka_mn2o_var_221(indm_var_966 + 1, ig_var_967) - ka_mn2o_var_221(indm_var_966, ig_var_967)))
          taug_var_937(jl_var_988, 88 + ig_var_967, lay_var_968) = colh2o_var_947(jl_var_988, lay_var_968) * (fac00_var_940(jl_var_988, lay_var_968) * absa_var_218(ind0_var_962, ig_var_967) + fac10_var_942(jl_var_988, lay_var_968) * absa_var_218(ind0_var_962 + 1, ig_var_967) + fac01_var_941(jl_var_988, lay_var_968) * absa_var_218(ind1_var_963, ig_var_967) + fac11_var_943(jl_var_988, lay_var_968) * absa_var_218(ind1_var_963 + 1, ig_var_967)) + tauself_var_974 + taufor_var_973 + adjcolco2_var_972 * absco2_var_976 + colo3_var_948(jl_var_988, lay_var_968) * abso3_var_975 + coln2o_var_949(jl_var_988, lay_var_968) * absn2o_var_977 + wx_var_938(jl_var_988, 3, lay_var_968) * cfc12_var_217(ig_var_967) + wx_var_938(jl_var_988, 4, lay_var_968) * cfc22adj(ig_var_967)
          fracs_var_956(jl_var_988, 88 + ig_var_967, lay_var_968) = fracrefa_var_215(ig_var_967)
        END DO
      END DO
      ixc0_var_985 = kfdia_var_935 - kidia_var_934 + 1 - ixc0_var_985
      DO ixp_var_986 = 1, ixc0_var_985
        jl_var_988 = ixhigh_var_982(ixp_var_986, lay_var_968)
        chi_co2_var_969 = colco2_var_950(jl_var_988, lay_var_968) / coldry_var_951(jl_var_988, lay_var_968)
        ratco2_var_970 = 1D+20 * chi_co2_var_969 / chi_mls(2, jp_var_944(jl_var_988, lay_var_968) + 1)
        IF (ratco2_var_970 .GT. 3.0D0) THEN
          adjfac_var_971 = 2.0D0 + (ratco2_var_970 - 2.0D0) ** 0.65D0
          adjcolco2_var_972 = adjfac_var_971 * chi_mls(2, jp_var_944(jl_var_988, lay_var_968) + 1) * coldry_var_951(jl_var_988, lay_var_968) * 1D-20
        ELSE
          adjcolco2_var_972 = colco2_var_950(jl_var_988, lay_var_968)
        END IF
        ind0_var_962 = ((jp_var_944(jl_var_988, lay_var_968) - 13) * 5 + (jt_var_945(jl_var_988, lay_var_968) - 1)) * nspb_var_238(8) + 1
        ind1_var_963 = ((jp_var_944(jl_var_988, lay_var_968) - 12) * 5 + (jt1_var_946(jl_var_988, lay_var_968) - 1)) * nspb_var_238(8) + 1
        indm_var_966 = indminor_var_961(jl_var_988, lay_var_968)
        DO ig_var_967 = 1, 8
          absco2_var_976 = (kb_mco2_var_223(indm_var_966, ig_var_967) + minorfrac_var_960(jl_var_988, lay_var_968) * (kb_mco2_var_223(indm_var_966 + 1, ig_var_967) - kb_mco2_var_223(indm_var_966, ig_var_967)))
          absn2o_var_977 = (kb_mn2o_var_224(indm_var_966, ig_var_967) + minorfrac_var_960(jl_var_988, lay_var_968) * (kb_mn2o_var_224(indm_var_966 + 1, ig_var_967) - kb_mn2o_var_224(indm_var_966, ig_var_967)))
          taug_var_937(jl_var_988, 88 + ig_var_967, lay_var_968) = colo3_var_948(jl_var_988, lay_var_968) * (fac00_var_940(jl_var_988, lay_var_968) * absb_var_219(ind0_var_962, ig_var_967) + fac10_var_942(jl_var_988, lay_var_968) * absb_var_219(ind0_var_962 + 1, ig_var_967) + fac01_var_941(jl_var_988, lay_var_968) * absb_var_219(ind1_var_963, ig_var_967) + fac11_var_943(jl_var_988, lay_var_968) * absb_var_219(ind1_var_963 + 1, ig_var_967)) + adjcolco2_var_972 * absco2_var_976 + coln2o_var_949(jl_var_988, lay_var_968) * absn2o_var_977 + wx_var_938(jl_var_988, 3, lay_var_968) * cfc12_var_217(ig_var_967) + wx_var_938(jl_var_988, 4, lay_var_968) * cfc22adj(ig_var_967)
          fracs_var_956(jl_var_988, 88 + ig_var_967, lay_var_968) = fracrefb_var_216(ig_var_967)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol8
SUBROUTINE rrtm_prepare_gases(kidia_var_990, kfdia_var_991, klon, klev_var_989, paph, pap, pth, pt, pq, pco2, pch4, pn2o, pno2, pc11, pc12, pc22, pcl4, pozn, pcoldry_var_992, pwbrodl, pwkl_var_993, pwx_var_994, pavel_var_995, ptavel_var_996, pz, ptz, kreflect)
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: klon
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_989
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_990, kfdia_var_991
  REAL(KIND = 8), INTENT(IN) :: paph(klon, klev_var_989 + 1)
  REAL(KIND = 8), INTENT(IN) :: pap(klon, klev_var_989)
  REAL(KIND = 8), INTENT(IN) :: pth(klon, klev_var_989 + 1)
  REAL(KIND = 8), INTENT(IN) :: pt(klon, klev_var_989)
  REAL(KIND = 8), INTENT(IN) :: pq(klon, klev_var_989)
  REAL(KIND = 8), INTENT(IN) :: pco2(klon, klev_var_989)
  REAL(KIND = 8), INTENT(IN) :: pch4(klon, klev_var_989)
  REAL(KIND = 8), INTENT(IN) :: pn2o(klon, klev_var_989)
  REAL(KIND = 8), INTENT(IN) :: pno2(klon, klev_var_989)
  REAL(KIND = 8), INTENT(IN) :: pc11(klon, klev_var_989)
  REAL(KIND = 8), INTENT(IN) :: pc12(klon, klev_var_989)
  REAL(KIND = 8), INTENT(IN) :: pc22(klon, klev_var_989)
  REAL(KIND = 8), INTENT(IN) :: pcl4(klon, klev_var_989)
  REAL(KIND = 8), INTENT(IN) :: pozn(klon, klev_var_989)
  REAL(KIND = 8), INTENT(OUT) :: pcoldry_var_992(kidia_var_990 : kfdia_var_991, klev_var_989)
  REAL(KIND = 8), INTENT(OUT) :: pwbrodl(kidia_var_990 : kfdia_var_991, klev_var_989)
  REAL(KIND = 8), INTENT(OUT) :: pwkl_var_993(kidia_var_990 : kfdia_var_991, 35, klev_var_989)
  REAL(KIND = 8), INTENT(OUT) :: pwx_var_994(kidia_var_990 : kfdia_var_991, 4, klev_var_989)
  REAL(KIND = 8), INTENT(OUT) :: pavel_var_995(kidia_var_990 : kfdia_var_991, klev_var_989)
  REAL(KIND = 8), INTENT(OUT) :: ptavel_var_996(kidia_var_990 : kfdia_var_991, klev_var_989)
  REAL(KIND = 8), INTENT(OUT) :: pz(kidia_var_990 : kfdia_var_991, 0 : klev_var_989)
  REAL(KIND = 8), INTENT(OUT) :: ptz(kidia_var_990 : kfdia_var_991, 0 : klev_var_989)
  INTEGER(KIND = 4), INTENT(OUT) :: kreflect(kidia_var_990 : kfdia_var_991)
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
  INTEGER(KIND = 4) :: iatm, jmol, ixmax, j1, j2, jk_var_997, jl_var_998
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
  DO jl_var_998 = kidia_var_990, kfdia_var_991
    kreflect(jl_var_998) = 0
  END DO
  DO j2 = 1, klev_var_989
    DO j1 = 1, 35
      DO jl_var_998 = kidia_var_990, kfdia_var_991
        pwkl_var_993(jl_var_998, j1, j2) = 0.0D0
      END DO
    END DO
  END DO
  iatm = 0
  ixmax = 4
  DO jl_var_998 = kidia_var_990, kfdia_var_991
    pz(jl_var_998, 0) = paph(jl_var_998, klev_var_989 + 1) / 100.0D0
    ptz(jl_var_998, 0) = pth(jl_var_998, klev_var_989 + 1)
  END DO
  DO jk_var_997 = 1, klev_var_989
    DO jl_var_998 = kidia_var_990, kfdia_var_991
      pavel_var_995(jl_var_998, jk_var_997) = pap(jl_var_998, klev_var_989 - jk_var_997 + 1) / 100.0D0
      ptavel_var_996(jl_var_998, jk_var_997) = pt(jl_var_998, klev_var_989 - jk_var_997 + 1)
      pz(jl_var_998, jk_var_997) = paph(jl_var_998, klev_var_989 - jk_var_997 + 1) / 100.0D0
      ptz(jl_var_998, jk_var_997) = pth(jl_var_998, klev_var_989 - jk_var_997 + 1)
      pwkl_var_993(jl_var_998, 1, jk_var_997) = MAX(pq(jl_var_998, klev_var_989 - jk_var_997 + 1), 1E-15) * zamd / zamw
      pwkl_var_993(jl_var_998, 2, jk_var_997) = pco2(jl_var_998, klev_var_989 - jk_var_997 + 1) * zamd / zamco2
      pwkl_var_993(jl_var_998, 3, jk_var_997) = pozn(jl_var_998, klev_var_989 - jk_var_997 + 1) * zamd / zamo
      pwkl_var_993(jl_var_998, 4, jk_var_997) = pn2o(jl_var_998, klev_var_989 - jk_var_997 + 1) * zamd / zamn2o
      pwkl_var_993(jl_var_998, 6, jk_var_997) = pch4(jl_var_998, klev_var_989 - jk_var_997 + 1) * zamd / zamch4
      pwkl_var_993(jl_var_998, 7, jk_var_997) = 0.209488D0
      zamm = (1.0D0 - pwkl_var_993(jl_var_998, 1, jk_var_997)) * zamd + pwkl_var_993(jl_var_998, 1, jk_var_997) * zamw
      pcoldry_var_992(jl_var_998, jk_var_997) = (pz(jl_var_998, jk_var_997 - 1) - pz(jl_var_998, jk_var_997)) * 1000.0D0 * zavgdro / (zgravit * zamm * (1.0D0 + pwkl_var_993(jl_var_998, 1, jk_var_997)))
    END DO
  END DO
  DO j2 = 1, klev_var_989
    DO j1 = 1, 4
      DO jl_var_998 = kidia_var_990, kfdia_var_991
        pwx_var_994(jl_var_998, j1, j2) = 0.0D0
      END DO
    END DO
  END DO
  DO jk_var_997 = 1, klev_var_989
    DO jl_var_998 = kidia_var_990, kfdia_var_991
      pwx_var_994(jl_var_998, 1, jk_var_997) = pcl4(jl_var_998, klev_var_989 - jk_var_997 + 1) * zamd / zamcl4
      pwx_var_994(jl_var_998, 2, jk_var_997) = pc11(jl_var_998, klev_var_989 - jk_var_997 + 1) * zamd / zamc11
      pwx_var_994(jl_var_998, 3, jk_var_997) = pc12(jl_var_998, klev_var_989 - jk_var_997 + 1) * zamd / zamc12
      pwx_var_994(jl_var_998, 4, jk_var_997) = pc22(jl_var_998, klev_var_989 - jk_var_997 + 1) * zamd / zamc22
      pwx_var_994(jl_var_998, 1, jk_var_997) = pcoldry_var_992(jl_var_998, jk_var_997) * pwx_var_994(jl_var_998, 1, jk_var_997) * 1D-20
      pwx_var_994(jl_var_998, 2, jk_var_997) = pcoldry_var_992(jl_var_998, jk_var_997) * pwx_var_994(jl_var_998, 2, jk_var_997) * 1D-20
      pwx_var_994(jl_var_998, 3, jk_var_997) = pcoldry_var_992(jl_var_998, jk_var_997) * pwx_var_994(jl_var_998, 3, jk_var_997) * 1D-20
      pwx_var_994(jl_var_998, 4, jk_var_997) = pcoldry_var_992(jl_var_998, jk_var_997) * pwx_var_994(jl_var_998, 4, jk_var_997) * 1D-20
      zsummol = 0.0D0
      DO jmol = 2, 7
        zsummol = zsummol + pwkl_var_993(jl_var_998, jmol, jk_var_997)
      END DO
      pwbrodl(jl_var_998, jk_var_997) = pcoldry_var_992(jl_var_998, jk_var_997) * (1.0D0 - zsummol)
      DO jmol = 1, 7
        pwkl_var_993(jl_var_998, jmol, jk_var_997) = pcoldry_var_992(jl_var_998, jk_var_997) * pwkl_var_993(jl_var_998, jmol, jk_var_997)
      END DO
    END DO
  END DO
END SUBROUTINE rrtm_prepare_gases
SUBROUTINE srtm_gas_optical_depth(kidia_var_999, kfdia_var_1000, klev_var_1001, poneminus_var_1002, prmu0_var_1003, klaytrop_var_1004, pcolch4_var_1005, pcolco2_var_1006, pcolh2o_var_1007, pcolmol_var_1008, pcolo2_var_1009, pcolo3_var_1010, pforfac_var_1011, pforfrac_var_1012, kindfor_var_1013, pselffac_var_1014, pselffrac_var_1015, kindself_var_1016, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, pod_var_1024, pssa, pincsol)
  USE yoesrtwn, ONLY: ngc
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_999, kfdia_var_1000
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1001
  REAL(KIND = 8), INTENT(IN) :: poneminus_var_1002(kidia_var_999 : kfdia_var_1000)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1003(kidia_var_999 : kfdia_var_1000)
  INTEGER(KIND = 4), INTENT(IN) :: klaytrop_var_1004(kidia_var_999 : kfdia_var_1000)
  REAL(KIND = 8), INTENT(IN) :: pcolch4_var_1005(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pcolco2_var_1006(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pcolh2o_var_1007(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pcolmol_var_1008(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pcolo2_var_1009(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pcolo3_var_1010(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pforfac_var_1011(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pforfrac_var_1012(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  INTEGER(KIND = 4), INTENT(IN) :: kindfor_var_1013(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pselffac_var_1014(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pselffrac_var_1015(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  INTEGER(KIND = 4), INTENT(IN) :: kindself_var_1016(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pfac00_var_1017(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pfac01_var_1018(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pfac10_var_1019(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(IN) :: pfac11_var_1020(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  INTEGER(KIND = 4), INTENT(IN) :: kjp_var_1021(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  INTEGER(KIND = 4), INTENT(IN) :: kjt_var_1022(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  INTEGER(KIND = 4), INTENT(IN) :: kjt1_var_1023(kidia_var_999 : kfdia_var_1000, klev_var_1001)
  REAL(KIND = 8), INTENT(INOUT) :: pod_var_1024(kidia_var_999 : kfdia_var_1000, klev_var_1001, 112)
  REAL(KIND = 8), INTENT(INOUT) :: pssa(kidia_var_999 : kfdia_var_1000, klev_var_1001, 112)
  REAL(KIND = 8), INTENT(INOUT) :: pincsol(kidia_var_999 : kfdia_var_1000, 112)
  INTEGER(KIND = 4) :: ib1, ib2, ibm, igt, iw(kidia_var_999 : kfdia_var_1000), jb_var_1025, jg_var_1026, jk_var_1027, jl_var_1028, icount
  REAL(KIND = 8) :: ztaug(kidia_var_999 : kfdia_var_1000, klev_var_1001, 16)
  REAL(KIND = 8) :: ztaur(kidia_var_999 : kfdia_var_1000, klev_var_1001, 16)
  REAL(KIND = 8) :: zsflxzen(kidia_var_999 : kfdia_var_1000, 16)
  ib1 = 16
  ib2 = 29
  icount = 0
  DO jl_var_1028 = kidia_var_999, kfdia_var_1000
    IF (prmu0_var_1003(jl_var_1028) > 0.0D0) THEN
      icount = icount + 1
      iw(jl_var_1028) = 0
    END IF
  END DO
  IF (icount /= 0) THEN
    DO jb_var_1025 = ib1, ib2
      ibm = jb_var_1025 - 15
      igt = ngc(ibm)
      IF (jb_var_1025 == 16) THEN
        CALL srtm_taumol16(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, poneminus_var_1002, pcolh2o_var_1007, pcolch4_var_1005, pcolmol_var_1008, klaytrop_var_1004, pselffac_var_1014, pselffrac_var_1015, kindself_var_1016, pforfac_var_1011, pforfrac_var_1012, kindfor_var_1013, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 17) THEN
        CALL srtm_taumol17(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, poneminus_var_1002, pcolh2o_var_1007, pcolco2_var_1006, pcolmol_var_1008, klaytrop_var_1004, pselffac_var_1014, pselffrac_var_1015, kindself_var_1016, pforfac_var_1011, pforfrac_var_1012, kindfor_var_1013, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 18) THEN
        CALL srtm_taumol18(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, poneminus_var_1002, pcolh2o_var_1007, pcolch4_var_1005, pcolmol_var_1008, klaytrop_var_1004, pselffac_var_1014, pselffrac_var_1015, kindself_var_1016, pforfac_var_1011, pforfrac_var_1012, kindfor_var_1013, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 19) THEN
        CALL srtm_taumol19(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, poneminus_var_1002, pcolh2o_var_1007, pcolco2_var_1006, pcolmol_var_1008, klaytrop_var_1004, pselffac_var_1014, pselffrac_var_1015, kindself_var_1016, pforfac_var_1011, pforfrac_var_1012, kindfor_var_1013, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 20) THEN
        CALL srtm_taumol20(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, pcolh2o_var_1007, pcolch4_var_1005, pcolmol_var_1008, klaytrop_var_1004, pselffac_var_1014, pselffrac_var_1015, kindself_var_1016, pforfac_var_1011, pforfrac_var_1012, kindfor_var_1013, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 21) THEN
        CALL srtm_taumol21(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, poneminus_var_1002, pcolh2o_var_1007, pcolco2_var_1006, pcolmol_var_1008, klaytrop_var_1004, pselffac_var_1014, pselffrac_var_1015, kindself_var_1016, pforfac_var_1011, pforfrac_var_1012, kindfor_var_1013, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 22) THEN
        CALL srtm_taumol22(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, poneminus_var_1002, pcolh2o_var_1007, pcolmol_var_1008, pcolo2_var_1009, klaytrop_var_1004, pselffac_var_1014, pselffrac_var_1015, kindself_var_1016, pforfac_var_1011, pforfrac_var_1012, kindfor_var_1013, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 23) THEN
        CALL srtm_taumol23(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, pcolh2o_var_1007, pcolmol_var_1008, klaytrop_var_1004, pselffac_var_1014, pselffrac_var_1015, kindself_var_1016, pforfac_var_1011, pforfrac_var_1012, kindfor_var_1013, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 24) THEN
        CALL srtm_taumol24(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, poneminus_var_1002, pcolh2o_var_1007, pcolmol_var_1008, pcolo2_var_1009, pcolo3_var_1010, klaytrop_var_1004, pselffac_var_1014, pselffrac_var_1015, kindself_var_1016, pforfac_var_1011, pforfrac_var_1012, kindfor_var_1013, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 25) THEN
        CALL srtm_taumol25(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, pcolh2o_var_1007, pcolmol_var_1008, pcolo3_var_1010, klaytrop_var_1004, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 26) THEN
        CALL srtm_taumol26(kidia_var_999, kfdia_var_1000, klev_var_1001, pcolmol_var_1008, klaytrop_var_1004, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 27) THEN
        CALL srtm_taumol27(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, pcolmol_var_1008, pcolo3_var_1010, klaytrop_var_1004, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 28) THEN
        CALL srtm_taumol28(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, poneminus_var_1002, pcolmol_var_1008, pcolo2_var_1009, pcolo3_var_1010, klaytrop_var_1004, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      ELSE IF (jb_var_1025 == 29) THEN
        CALL srtm_taumol29(kidia_var_999, kfdia_var_1000, klev_var_1001, pfac00_var_1017, pfac01_var_1018, pfac10_var_1019, pfac11_var_1020, kjp_var_1021, kjt_var_1022, kjt1_var_1023, pcolh2o_var_1007, pcolco2_var_1006, pcolmol_var_1008, klaytrop_var_1004, pselffac_var_1014, pselffrac_var_1015, kindself_var_1016, pforfac_var_1011, pforfrac_var_1012, kindfor_var_1013, zsflxzen, ztaug, ztaur, prmu0_var_1003)
      END IF
      DO jg_var_1026 = 1, igt
        DO jl_var_1028 = kidia_var_999, kfdia_var_1000
          IF (prmu0_var_1003(jl_var_1028) > 0.0D0) THEN
            iw(jl_var_1028) = iw(jl_var_1028) + 1
            pincsol(jl_var_1028, iw(jl_var_1028)) = zsflxzen(jl_var_1028, jg_var_1026)
          END IF
        END DO
        DO jk_var_1027 = 1, klev_var_1001
          DO jl_var_1028 = kidia_var_999, kfdia_var_1000
            IF (prmu0_var_1003(jl_var_1028) > 0.0D0) THEN
              pod_var_1024(jl_var_1028, jk_var_1027, iw(jl_var_1028)) = ztaur(jl_var_1028, jk_var_1027, jg_var_1026) + ztaug(jl_var_1028, jk_var_1027, jg_var_1026)
              pssa(jl_var_1028, jk_var_1027, iw(jl_var_1028)) = ztaur(jl_var_1028, jk_var_1027, jg_var_1026) / pod_var_1024(jl_var_1028, jk_var_1027, iw(jl_var_1028))
            END IF
          END DO
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE srtm_gas_optical_depth
SUBROUTINE srtm_taumol28(kidia_var_1029, kfdia_var_1030, klev_var_1031, p_fac00_var_1032, p_fac01_var_1033, p_fac10_var_1034, p_fac11_var_1035, k_jp_var_1036, k_jt_var_1037, k_jt1_var_1038, p_oneminus_var_1039, p_colmol_var_1040, p_colo2_var_1041, p_colo3_var_1042, k_laytrop_var_1043, p_sfluxzen_var_1044, p_taug_var_1045, p_taur_var_1046, prmu0_var_1047)
  USE yoesrta28, ONLY: absa_var_324, absb_var_325, layreffr_var_323, rayl_var_321, sfluxrefc_var_326, strrat_var_322
  USE yoesrtwn, ONLY: nspa_var_334, nspb_var_335
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1029, kfdia_var_1030
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1031
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1032(kidia_var_1029 : kfdia_var_1030, klev_var_1031)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1033(kidia_var_1029 : kfdia_var_1030, klev_var_1031)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1034(kidia_var_1029 : kfdia_var_1030, klev_var_1031)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1035(kidia_var_1029 : kfdia_var_1030, klev_var_1031)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1036(kidia_var_1029 : kfdia_var_1030, klev_var_1031)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1037(kidia_var_1029 : kfdia_var_1030, klev_var_1031)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1038(kidia_var_1029 : kfdia_var_1030, klev_var_1031)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1039(kidia_var_1029 : kfdia_var_1030)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1040(kidia_var_1029 : kfdia_var_1030, klev_var_1031)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_1041(kidia_var_1029 : kfdia_var_1030, klev_var_1031)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1042(kidia_var_1029 : kfdia_var_1030, klev_var_1031)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1043(kidia_var_1029 : kfdia_var_1030)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1044(kidia_var_1029 : kfdia_var_1030, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1045(kidia_var_1029 : kfdia_var_1030, klev_var_1031, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1046(kidia_var_1029 : kfdia_var_1030, klev_var_1031, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1047(kidia_var_1029 : kfdia_var_1030)
  INTEGER(KIND = 4) :: ig_var_1048, ind0_var_1049, ind1_var_1050, js_var_1051, i_lay_var_1052, i_laysolfr_var_1053(kidia_var_1029 : kfdia_var_1030), i_nlayers_var_1054, iplon_var_1055
  INTEGER(KIND = 4) :: laytrop_min_var_1056, laytrop_max_var_1057
  REAL(KIND = 8) :: z_fs_var_1058, z_speccomb_var_1059, z_specmult_var_1060, z_specparm_var_1061, z_tauray_var_1062
  laytrop_min_var_1056 = MINVAL(k_laytrop_var_1043(kidia_var_1029 : kfdia_var_1030))
  laytrop_max_var_1057 = MAXVAL(k_laytrop_var_1043(kidia_var_1029 : kfdia_var_1030))
  i_nlayers_var_1054 = klev_var_1031
  DO iplon_var_1055 = kidia_var_1029, kfdia_var_1030
    i_laysolfr_var_1053(iplon_var_1055) = i_nlayers_var_1054
  END DO
  DO i_lay_var_1052 = 1, laytrop_min_var_1056
    DO iplon_var_1055 = kidia_var_1029, kfdia_var_1030
      z_speccomb_var_1059 = p_colo3_var_1042(iplon_var_1055, i_lay_var_1052) + strrat_var_322 * p_colo2_var_1041(iplon_var_1055, i_lay_var_1052)
      z_specparm_var_1061 = p_colo3_var_1042(iplon_var_1055, i_lay_var_1052) / z_speccomb_var_1059
      z_specparm_var_1061 = MIN(p_oneminus_var_1039(iplon_var_1055), z_specparm_var_1061)
      z_specmult_var_1060 = 8.0D0 * (z_specparm_var_1061)
      js_var_1051 = 1 + INT(z_specmult_var_1060)
      z_fs_var_1058 = z_specmult_var_1060 - AINT(z_specmult_var_1060)
      ind0_var_1049 = ((k_jp_var_1036(iplon_var_1055, i_lay_var_1052) - 1) * 5 + (k_jt_var_1037(iplon_var_1055, i_lay_var_1052) - 1)) * nspa_var_334(28) + js_var_1051
      ind1_var_1050 = (k_jp_var_1036(iplon_var_1055, i_lay_var_1052) * 5 + (k_jt1_var_1038(iplon_var_1055, i_lay_var_1052) - 1)) * nspa_var_334(28) + js_var_1051
      z_tauray_var_1062 = p_colmol_var_1040(iplon_var_1055, i_lay_var_1052) * rayl_var_321
      DO ig_var_1048 = 1, 6
        p_taug_var_1045(iplon_var_1055, i_lay_var_1052, ig_var_1048) = z_speccomb_var_1059 * ((1.0D0 - z_fs_var_1058) * (absa_var_324(ind0_var_1049, ig_var_1048) * p_fac00_var_1032(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind0_var_1049 + 9, ig_var_1048) * p_fac10_var_1034(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind1_var_1050, ig_var_1048) * p_fac01_var_1033(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind1_var_1050 + 9, ig_var_1048) * p_fac11_var_1035(iplon_var_1055, i_lay_var_1052)) + z_fs_var_1058 * (absa_var_324(ind0_var_1049 + 1, ig_var_1048) * p_fac00_var_1032(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind0_var_1049 + 10, ig_var_1048) * p_fac10_var_1034(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind1_var_1050 + 1, ig_var_1048) * p_fac01_var_1033(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind1_var_1050 + 10, ig_var_1048) * p_fac11_var_1035(iplon_var_1055, i_lay_var_1052)))
        p_taur_var_1046(iplon_var_1055, i_lay_var_1052, ig_var_1048) = z_tauray_var_1062
      END DO
    END DO
  END DO
  DO i_lay_var_1052 = laytrop_min_var_1056 + 1, laytrop_max_var_1057
    DO iplon_var_1055 = kidia_var_1029, kfdia_var_1030
      IF (i_lay_var_1052 <= k_laytrop_var_1043(iplon_var_1055)) THEN
        z_speccomb_var_1059 = p_colo3_var_1042(iplon_var_1055, i_lay_var_1052) + strrat_var_322 * p_colo2_var_1041(iplon_var_1055, i_lay_var_1052)
        z_specparm_var_1061 = p_colo3_var_1042(iplon_var_1055, i_lay_var_1052) / z_speccomb_var_1059
        z_specparm_var_1061 = MIN(p_oneminus_var_1039(iplon_var_1055), z_specparm_var_1061)
        z_specmult_var_1060 = 8.0D0 * (z_specparm_var_1061)
        js_var_1051 = 1 + INT(z_specmult_var_1060)
        z_fs_var_1058 = z_specmult_var_1060 - AINT(z_specmult_var_1060)
        ind0_var_1049 = ((k_jp_var_1036(iplon_var_1055, i_lay_var_1052) - 1) * 5 + (k_jt_var_1037(iplon_var_1055, i_lay_var_1052) - 1)) * nspa_var_334(28) + js_var_1051
        ind1_var_1050 = (k_jp_var_1036(iplon_var_1055, i_lay_var_1052) * 5 + (k_jt1_var_1038(iplon_var_1055, i_lay_var_1052) - 1)) * nspa_var_334(28) + js_var_1051
        z_tauray_var_1062 = p_colmol_var_1040(iplon_var_1055, i_lay_var_1052) * rayl_var_321
        DO ig_var_1048 = 1, 6
          p_taug_var_1045(iplon_var_1055, i_lay_var_1052, ig_var_1048) = z_speccomb_var_1059 * ((1.0D0 - z_fs_var_1058) * (absa_var_324(ind0_var_1049, ig_var_1048) * p_fac00_var_1032(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind0_var_1049 + 9, ig_var_1048) * p_fac10_var_1034(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind1_var_1050, ig_var_1048) * p_fac01_var_1033(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind1_var_1050 + 9, ig_var_1048) * p_fac11_var_1035(iplon_var_1055, i_lay_var_1052)) + z_fs_var_1058 * (absa_var_324(ind0_var_1049 + 1, ig_var_1048) * p_fac00_var_1032(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind0_var_1049 + 10, ig_var_1048) * p_fac10_var_1034(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind1_var_1050 + 1, ig_var_1048) * p_fac01_var_1033(iplon_var_1055, i_lay_var_1052) + absa_var_324(ind1_var_1050 + 10, ig_var_1048) * p_fac11_var_1035(iplon_var_1055, i_lay_var_1052)))
          p_taur_var_1046(iplon_var_1055, i_lay_var_1052, ig_var_1048) = z_tauray_var_1062
        END DO
      ELSE
        IF (k_jp_var_1036(iplon_var_1055, i_lay_var_1052 - 1) < layreffr_var_323 .AND. k_jp_var_1036(iplon_var_1055, i_lay_var_1052) >= layreffr_var_323) i_laysolfr_var_1053(iplon_var_1055) = i_lay_var_1052
        z_speccomb_var_1059 = p_colo3_var_1042(iplon_var_1055, i_lay_var_1052) + strrat_var_322 * p_colo2_var_1041(iplon_var_1055, i_lay_var_1052)
        z_specparm_var_1061 = p_colo3_var_1042(iplon_var_1055, i_lay_var_1052) / z_speccomb_var_1059
        z_specparm_var_1061 = MIN(p_oneminus_var_1039(iplon_var_1055), z_specparm_var_1061)
        z_specmult_var_1060 = 4.0D0 * (z_specparm_var_1061)
        js_var_1051 = 1 + INT(z_specmult_var_1060)
        z_fs_var_1058 = z_specmult_var_1060 - AINT(z_specmult_var_1060)
        ind0_var_1049 = ((k_jp_var_1036(iplon_var_1055, i_lay_var_1052) - 13) * 5 + (k_jt_var_1037(iplon_var_1055, i_lay_var_1052) - 1)) * nspb_var_335(28) + js_var_1051
        ind1_var_1050 = ((k_jp_var_1036(iplon_var_1055, i_lay_var_1052) - 12) * 5 + (k_jt1_var_1038(iplon_var_1055, i_lay_var_1052) - 1)) * nspb_var_335(28) + js_var_1051
        z_tauray_var_1062 = p_colmol_var_1040(iplon_var_1055, i_lay_var_1052) * rayl_var_321
        DO ig_var_1048 = 1, 6
          p_taug_var_1045(iplon_var_1055, i_lay_var_1052, ig_var_1048) = z_speccomb_var_1059 * ((1.0D0 - z_fs_var_1058) * (absb_var_325(ind0_var_1049, ig_var_1048) * p_fac00_var_1032(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind0_var_1049 + 5, ig_var_1048) * p_fac10_var_1034(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind1_var_1050, ig_var_1048) * p_fac01_var_1033(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind1_var_1050 + 5, ig_var_1048) * p_fac11_var_1035(iplon_var_1055, i_lay_var_1052)) + z_fs_var_1058 * (absb_var_325(ind0_var_1049 + 1, ig_var_1048) * p_fac00_var_1032(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind0_var_1049 + 6, ig_var_1048) * p_fac10_var_1034(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind1_var_1050 + 1, ig_var_1048) * p_fac01_var_1033(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind1_var_1050 + 6, ig_var_1048) * p_fac11_var_1035(iplon_var_1055, i_lay_var_1052)))
          IF (i_lay_var_1052 == i_laysolfr_var_1053(iplon_var_1055)) p_sfluxzen_var_1044(iplon_var_1055, ig_var_1048) = sfluxrefc_var_326(ig_var_1048, js_var_1051) + z_fs_var_1058 * (sfluxrefc_var_326(ig_var_1048, js_var_1051 + 1) - sfluxrefc_var_326(ig_var_1048, js_var_1051))
          p_taur_var_1046(iplon_var_1055, i_lay_var_1052, ig_var_1048) = z_tauray_var_1062
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1052 = laytrop_max_var_1057 + 1, i_nlayers_var_1054
    DO iplon_var_1055 = kidia_var_1029, kfdia_var_1030
      IF (k_jp_var_1036(iplon_var_1055, i_lay_var_1052 - 1) < layreffr_var_323 .AND. k_jp_var_1036(iplon_var_1055, i_lay_var_1052) >= layreffr_var_323) i_laysolfr_var_1053(iplon_var_1055) = i_lay_var_1052
      z_speccomb_var_1059 = p_colo3_var_1042(iplon_var_1055, i_lay_var_1052) + strrat_var_322 * p_colo2_var_1041(iplon_var_1055, i_lay_var_1052)
      z_specparm_var_1061 = p_colo3_var_1042(iplon_var_1055, i_lay_var_1052) / z_speccomb_var_1059
      z_specparm_var_1061 = MIN(p_oneminus_var_1039(iplon_var_1055), z_specparm_var_1061)
      z_specmult_var_1060 = 4.0D0 * (z_specparm_var_1061)
      js_var_1051 = 1 + INT(z_specmult_var_1060)
      z_fs_var_1058 = z_specmult_var_1060 - AINT(z_specmult_var_1060)
      ind0_var_1049 = ((k_jp_var_1036(iplon_var_1055, i_lay_var_1052) - 13) * 5 + (k_jt_var_1037(iplon_var_1055, i_lay_var_1052) - 1)) * nspb_var_335(28) + js_var_1051
      ind1_var_1050 = ((k_jp_var_1036(iplon_var_1055, i_lay_var_1052) - 12) * 5 + (k_jt1_var_1038(iplon_var_1055, i_lay_var_1052) - 1)) * nspb_var_335(28) + js_var_1051
      z_tauray_var_1062 = p_colmol_var_1040(iplon_var_1055, i_lay_var_1052) * rayl_var_321
      DO ig_var_1048 = 1, 6
        p_taug_var_1045(iplon_var_1055, i_lay_var_1052, ig_var_1048) = z_speccomb_var_1059 * ((1.0D0 - z_fs_var_1058) * (absb_var_325(ind0_var_1049, ig_var_1048) * p_fac00_var_1032(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind0_var_1049 + 5, ig_var_1048) * p_fac10_var_1034(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind1_var_1050, ig_var_1048) * p_fac01_var_1033(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind1_var_1050 + 5, ig_var_1048) * p_fac11_var_1035(iplon_var_1055, i_lay_var_1052)) + z_fs_var_1058 * (absb_var_325(ind0_var_1049 + 1, ig_var_1048) * p_fac00_var_1032(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind0_var_1049 + 6, ig_var_1048) * p_fac10_var_1034(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind1_var_1050 + 1, ig_var_1048) * p_fac01_var_1033(iplon_var_1055, i_lay_var_1052) + absb_var_325(ind1_var_1050 + 6, ig_var_1048) * p_fac11_var_1035(iplon_var_1055, i_lay_var_1052)))
        IF (i_lay_var_1052 == i_laysolfr_var_1053(iplon_var_1055)) p_sfluxzen_var_1044(iplon_var_1055, ig_var_1048) = sfluxrefc_var_326(ig_var_1048, js_var_1051) + z_fs_var_1058 * (sfluxrefc_var_326(ig_var_1048, js_var_1051 + 1) - sfluxrefc_var_326(ig_var_1048, js_var_1051))
        p_taur_var_1046(iplon_var_1055, i_lay_var_1052, ig_var_1048) = z_tauray_var_1062
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol28
SUBROUTINE srtm_taumol29(kidia_var_1063, kfdia_var_1064, klev_var_1065, p_fac00_var_1066, p_fac01_var_1067, p_fac10_var_1068, p_fac11_var_1069, k_jp_var_1070, k_jt_var_1071, k_jt1_var_1072, p_colh2o_var_1073, p_colco2_var_1074, p_colmol_var_1075, k_laytrop_var_1076, p_selffac_var_1077, p_selffrac_var_1078, k_indself_var_1079, p_forfac_var_1080, p_forfrac_var_1081, k_indfor_var_1082, p_sfluxzen_var_1085, p_taug_var_1086, p_taur_var_1087, prmu0_var_1088)
  USE yoesrtwn, ONLY: nspa_var_334, nspb_var_335
  USE yoesrta29, ONLY: absa_var_329, absb_var_330, absco2c, absh2oc, forrefc_var_332, layreffr_var_328, rayl_var_327, selfrefc_var_331, sfluxrefc_var_333
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1063, kfdia_var_1064
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1065
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1066(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1067(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1068(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1069(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1070(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1071(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1072(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1073(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1074(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1075(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1076(kidia_var_1063 : kfdia_var_1064)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1077(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1078(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1079(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1080(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1081(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1082(kidia_var_1063 : kfdia_var_1064, klev_var_1065)
  INTEGER(KIND = 4) :: laytrop_min_var_1083, laytrop_max_var_1084
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1085(kidia_var_1063 : kfdia_var_1064, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1086(kidia_var_1063 : kfdia_var_1064, klev_var_1065, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1087(kidia_var_1063 : kfdia_var_1064, klev_var_1065, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1088(kidia_var_1063 : kfdia_var_1064)
  INTEGER(KIND = 4) :: ig_var_1089, ind0_var_1090, ind1_var_1091, inds_var_1092, indf_var_1093, i_lay_var_1094, i_laysolfr_var_1095(kidia_var_1063 : kfdia_var_1064), i_nlayers_var_1096, iplon_var_1097
  REAL(KIND = 8) :: z_tauray_var_1098
  laytrop_min_var_1083 = MINVAL(k_laytrop_var_1076(kidia_var_1063 : kfdia_var_1064))
  laytrop_max_var_1084 = MAXVAL(k_laytrop_var_1076(kidia_var_1063 : kfdia_var_1064))
  i_nlayers_var_1096 = klev_var_1065
  DO iplon_var_1097 = kidia_var_1063, kfdia_var_1064
    i_laysolfr_var_1095(iplon_var_1097) = i_nlayers_var_1096
  END DO
  DO i_lay_var_1094 = 1, laytrop_min_var_1083
    DO iplon_var_1097 = kidia_var_1063, kfdia_var_1064
      ind0_var_1090 = ((k_jp_var_1070(iplon_var_1097, i_lay_var_1094) - 1) * 5 + (k_jt_var_1071(iplon_var_1097, i_lay_var_1094) - 1)) * nspa_var_334(29) + 1
      ind1_var_1091 = (k_jp_var_1070(iplon_var_1097, i_lay_var_1094) * 5 + (k_jt1_var_1072(iplon_var_1097, i_lay_var_1094) - 1)) * nspa_var_334(29) + 1
      inds_var_1092 = k_indself_var_1079(iplon_var_1097, i_lay_var_1094)
      indf_var_1093 = k_indfor_var_1082(iplon_var_1097, i_lay_var_1094)
      z_tauray_var_1098 = p_colmol_var_1075(iplon_var_1097, i_lay_var_1094) * rayl_var_327
      DO ig_var_1089 = 1, 12
        p_taug_var_1086(iplon_var_1097, i_lay_var_1094, ig_var_1089) = p_colh2o_var_1073(iplon_var_1097, i_lay_var_1094) * ((p_fac00_var_1066(iplon_var_1097, i_lay_var_1094) * absa_var_329(ind0_var_1090, ig_var_1089) + p_fac10_var_1068(iplon_var_1097, i_lay_var_1094) * absa_var_329(ind0_var_1090 + 1, ig_var_1089) + p_fac01_var_1067(iplon_var_1097, i_lay_var_1094) * absa_var_329(ind1_var_1091, ig_var_1089) + p_fac11_var_1069(iplon_var_1097, i_lay_var_1094) * absa_var_329(ind1_var_1091 + 1, ig_var_1089)) + p_selffac_var_1077(iplon_var_1097, i_lay_var_1094) * (selfrefc_var_331(inds_var_1092, ig_var_1089) + p_selffrac_var_1078(iplon_var_1097, i_lay_var_1094) * (selfrefc_var_331(inds_var_1092 + 1, ig_var_1089) - selfrefc_var_331(inds_var_1092, ig_var_1089))) + p_forfac_var_1080(iplon_var_1097, i_lay_var_1094) * (forrefc_var_332(indf_var_1093, ig_var_1089) + p_forfrac_var_1081(iplon_var_1097, i_lay_var_1094) * (forrefc_var_332(indf_var_1093 + 1, ig_var_1089) - forrefc_var_332(indf_var_1093, ig_var_1089)))) + p_colco2_var_1074(iplon_var_1097, i_lay_var_1094) * absco2c(ig_var_1089)
        p_taur_var_1087(iplon_var_1097, i_lay_var_1094, ig_var_1089) = z_tauray_var_1098
      END DO
    END DO
  END DO
  DO i_lay_var_1094 = laytrop_min_var_1083 + 1, laytrop_max_var_1084
    DO iplon_var_1097 = kidia_var_1063, kfdia_var_1064
      IF (i_lay_var_1094 <= k_laytrop_var_1076(iplon_var_1097)) THEN
        ind0_var_1090 = ((k_jp_var_1070(iplon_var_1097, i_lay_var_1094) - 1) * 5 + (k_jt_var_1071(iplon_var_1097, i_lay_var_1094) - 1)) * nspa_var_334(29) + 1
        ind1_var_1091 = (k_jp_var_1070(iplon_var_1097, i_lay_var_1094) * 5 + (k_jt1_var_1072(iplon_var_1097, i_lay_var_1094) - 1)) * nspa_var_334(29) + 1
        inds_var_1092 = k_indself_var_1079(iplon_var_1097, i_lay_var_1094)
        indf_var_1093 = k_indfor_var_1082(iplon_var_1097, i_lay_var_1094)
        z_tauray_var_1098 = p_colmol_var_1075(iplon_var_1097, i_lay_var_1094) * rayl_var_327
        DO ig_var_1089 = 1, 12
          p_taug_var_1086(iplon_var_1097, i_lay_var_1094, ig_var_1089) = p_colh2o_var_1073(iplon_var_1097, i_lay_var_1094) * ((p_fac00_var_1066(iplon_var_1097, i_lay_var_1094) * absa_var_329(ind0_var_1090, ig_var_1089) + p_fac10_var_1068(iplon_var_1097, i_lay_var_1094) * absa_var_329(ind0_var_1090 + 1, ig_var_1089) + p_fac01_var_1067(iplon_var_1097, i_lay_var_1094) * absa_var_329(ind1_var_1091, ig_var_1089) + p_fac11_var_1069(iplon_var_1097, i_lay_var_1094) * absa_var_329(ind1_var_1091 + 1, ig_var_1089)) + p_selffac_var_1077(iplon_var_1097, i_lay_var_1094) * (selfrefc_var_331(inds_var_1092, ig_var_1089) + p_selffrac_var_1078(iplon_var_1097, i_lay_var_1094) * (selfrefc_var_331(inds_var_1092 + 1, ig_var_1089) - selfrefc_var_331(inds_var_1092, ig_var_1089))) + p_forfac_var_1080(iplon_var_1097, i_lay_var_1094) * (forrefc_var_332(indf_var_1093, ig_var_1089) + p_forfrac_var_1081(iplon_var_1097, i_lay_var_1094) * (forrefc_var_332(indf_var_1093 + 1, ig_var_1089) - forrefc_var_332(indf_var_1093, ig_var_1089)))) + p_colco2_var_1074(iplon_var_1097, i_lay_var_1094) * absco2c(ig_var_1089)
          p_taur_var_1087(iplon_var_1097, i_lay_var_1094, ig_var_1089) = z_tauray_var_1098
        END DO
      ELSE
        IF (k_jp_var_1070(iplon_var_1097, i_lay_var_1094 - 1) < layreffr_var_328 .AND. k_jp_var_1070(iplon_var_1097, i_lay_var_1094) >= layreffr_var_328) i_laysolfr_var_1095(iplon_var_1097) = i_lay_var_1094
        ind0_var_1090 = ((k_jp_var_1070(iplon_var_1097, i_lay_var_1094) - 13) * 5 + (k_jt_var_1071(iplon_var_1097, i_lay_var_1094) - 1)) * nspb_var_335(29) + 1
        ind1_var_1091 = ((k_jp_var_1070(iplon_var_1097, i_lay_var_1094) - 12) * 5 + (k_jt1_var_1072(iplon_var_1097, i_lay_var_1094) - 1)) * nspb_var_335(29) + 1
        z_tauray_var_1098 = p_colmol_var_1075(iplon_var_1097, i_lay_var_1094) * rayl_var_327
        DO ig_var_1089 = 1, 12
          p_taug_var_1086(iplon_var_1097, i_lay_var_1094, ig_var_1089) = p_colco2_var_1074(iplon_var_1097, i_lay_var_1094) * (p_fac00_var_1066(iplon_var_1097, i_lay_var_1094) * absb_var_330(ind0_var_1090, ig_var_1089) + p_fac10_var_1068(iplon_var_1097, i_lay_var_1094) * absb_var_330(ind0_var_1090 + 1, ig_var_1089) + p_fac01_var_1067(iplon_var_1097, i_lay_var_1094) * absb_var_330(ind1_var_1091, ig_var_1089) + p_fac11_var_1069(iplon_var_1097, i_lay_var_1094) * absb_var_330(ind1_var_1091 + 1, ig_var_1089)) + p_colh2o_var_1073(iplon_var_1097, i_lay_var_1094) * absh2oc(ig_var_1089)
          IF (i_lay_var_1094 == i_laysolfr_var_1095(iplon_var_1097)) p_sfluxzen_var_1085(iplon_var_1097, ig_var_1089) = sfluxrefc_var_333(ig_var_1089)
          p_taur_var_1087(iplon_var_1097, i_lay_var_1094, ig_var_1089) = z_tauray_var_1098
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1094 = laytrop_max_var_1084 + 1, i_nlayers_var_1096
    DO iplon_var_1097 = kidia_var_1063, kfdia_var_1064
      IF (k_jp_var_1070(iplon_var_1097, i_lay_var_1094 - 1) < layreffr_var_328 .AND. k_jp_var_1070(iplon_var_1097, i_lay_var_1094) >= layreffr_var_328) i_laysolfr_var_1095(iplon_var_1097) = i_lay_var_1094
      ind0_var_1090 = ((k_jp_var_1070(iplon_var_1097, i_lay_var_1094) - 13) * 5 + (k_jt_var_1071(iplon_var_1097, i_lay_var_1094) - 1)) * nspb_var_335(29) + 1
      ind1_var_1091 = ((k_jp_var_1070(iplon_var_1097, i_lay_var_1094) - 12) * 5 + (k_jt1_var_1072(iplon_var_1097, i_lay_var_1094) - 1)) * nspb_var_335(29) + 1
      z_tauray_var_1098 = p_colmol_var_1075(iplon_var_1097, i_lay_var_1094) * rayl_var_327
      DO ig_var_1089 = 1, 12
        p_taug_var_1086(iplon_var_1097, i_lay_var_1094, ig_var_1089) = p_colco2_var_1074(iplon_var_1097, i_lay_var_1094) * (p_fac00_var_1066(iplon_var_1097, i_lay_var_1094) * absb_var_330(ind0_var_1090, ig_var_1089) + p_fac10_var_1068(iplon_var_1097, i_lay_var_1094) * absb_var_330(ind0_var_1090 + 1, ig_var_1089) + p_fac01_var_1067(iplon_var_1097, i_lay_var_1094) * absb_var_330(ind1_var_1091, ig_var_1089) + p_fac11_var_1069(iplon_var_1097, i_lay_var_1094) * absb_var_330(ind1_var_1091 + 1, ig_var_1089)) + p_colh2o_var_1073(iplon_var_1097, i_lay_var_1094) * absh2oc(ig_var_1089)
        IF (i_lay_var_1094 == i_laysolfr_var_1095(iplon_var_1097)) p_sfluxzen_var_1085(iplon_var_1097, ig_var_1089) = sfluxrefc_var_333(ig_var_1089)
        p_taur_var_1087(iplon_var_1097, i_lay_var_1094, ig_var_1089) = z_tauray_var_1098
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol29
SUBROUTINE srtm_taumol17(kidia_var_1099, kfdia_var_1100, klev_var_1101, p_fac00_var_1102, p_fac01_var_1103, p_fac10_var_1104, p_fac11_var_1105, k_jp_var_1106, k_jt_var_1107, k_jt1_var_1108, p_oneminus_var_1109, p_colh2o_var_1110, p_colco2_var_1111, p_colmol_var_1112, k_laytrop_var_1113, p_selffac_var_1114, p_selffrac_var_1115, k_indself_var_1116, p_forfac_var_1117, p_forfrac_var_1118, k_indfor_var_1119, p_sfluxzen_var_1120, p_taug_var_1121, p_taur_var_1122, prmu0_var_1123)
  USE yoesrta17, ONLY: absa_var_249, absb_var_250, forrefc_var_252, layreffr_var_248, rayl_var_246, selfrefc_var_251, sfluxrefc_var_253, strrat_var_247
  USE yoesrtwn, ONLY: nspa_var_334, nspb_var_335
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1099, kfdia_var_1100
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1101
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1102(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1103(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1104(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1105(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1106(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1107(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1108(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1109(kidia_var_1099 : kfdia_var_1100)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1110(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1111(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1112(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1113(kidia_var_1099 : kfdia_var_1100)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1114(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1115(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1116(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1117(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1118(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1119(kidia_var_1099 : kfdia_var_1100, klev_var_1101)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1120(kidia_var_1099 : kfdia_var_1100, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1121(kidia_var_1099 : kfdia_var_1100, klev_var_1101, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1122(kidia_var_1099 : kfdia_var_1100, klev_var_1101, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1123(kidia_var_1099 : kfdia_var_1100)
  INTEGER(KIND = 4) :: ig_var_1124, ind0_var_1125, ind1_var_1126, inds_var_1127, indf_var_1128, js_var_1129, i_lay_var_1130, i_laysolfr_var_1131(kidia_var_1099 : kfdia_var_1100), i_nlayers_var_1132, iplon_var_1133
  INTEGER(KIND = 4) :: laytrop_min_var_1134, laytrop_max_var_1135
  REAL(KIND = 8) :: z_fs_var_1136, z_speccomb_var_1137, z_specmult_var_1138, z_specparm_var_1139, z_tauray_var_1140
  laytrop_min_var_1134 = MINVAL(k_laytrop_var_1113(kidia_var_1099 : kfdia_var_1100))
  laytrop_max_var_1135 = MAXVAL(k_laytrop_var_1113(kidia_var_1099 : kfdia_var_1100))
  i_nlayers_var_1132 = klev_var_1101
  DO iplon_var_1133 = kidia_var_1099, kfdia_var_1100
    i_laysolfr_var_1131(iplon_var_1133) = i_nlayers_var_1132
  END DO
  DO i_lay_var_1130 = 1, laytrop_min_var_1134
    DO iplon_var_1133 = kidia_var_1099, kfdia_var_1100
      z_speccomb_var_1137 = p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) + strrat_var_247 * p_colco2_var_1111(iplon_var_1133, i_lay_var_1130)
      z_specparm_var_1139 = p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) / z_speccomb_var_1137
      z_specparm_var_1139 = MIN(p_oneminus_var_1109(iplon_var_1133), z_specparm_var_1139)
      z_specmult_var_1138 = 8.0D0 * (z_specparm_var_1139)
      js_var_1129 = 1 + INT(z_specmult_var_1138)
      z_fs_var_1136 = z_specmult_var_1138 - AINT(z_specmult_var_1138)
      ind0_var_1125 = ((k_jp_var_1106(iplon_var_1133, i_lay_var_1130) - 1) * 5 + (k_jt_var_1107(iplon_var_1133, i_lay_var_1130) - 1)) * nspa_var_334(17) + js_var_1129
      ind1_var_1126 = (k_jp_var_1106(iplon_var_1133, i_lay_var_1130) * 5 + (k_jt1_var_1108(iplon_var_1133, i_lay_var_1130) - 1)) * nspa_var_334(17) + js_var_1129
      inds_var_1127 = k_indself_var_1116(iplon_var_1133, i_lay_var_1130)
      indf_var_1128 = k_indfor_var_1119(iplon_var_1133, i_lay_var_1130)
      z_tauray_var_1140 = p_colmol_var_1112(iplon_var_1133, i_lay_var_1130) * rayl_var_246
      DO ig_var_1124 = 1, 12
        p_taug_var_1121(iplon_var_1133, i_lay_var_1130, ig_var_1124) = z_speccomb_var_1137 * ((1.0D0 - z_fs_var_1136) * (absa_var_249(ind0_var_1125, ig_var_1124) * p_fac00_var_1102(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind0_var_1125 + 9, ig_var_1124) * p_fac10_var_1104(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind1_var_1126, ig_var_1124) * p_fac01_var_1103(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind1_var_1126 + 9, ig_var_1124) * p_fac11_var_1105(iplon_var_1133, i_lay_var_1130)) + z_fs_var_1136 * (absa_var_249(ind0_var_1125 + 1, ig_var_1124) * p_fac00_var_1102(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind0_var_1125 + 10, ig_var_1124) * p_fac10_var_1104(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind1_var_1126 + 1, ig_var_1124) * p_fac01_var_1103(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind1_var_1126 + 10, ig_var_1124) * p_fac11_var_1105(iplon_var_1133, i_lay_var_1130))) + p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) * (p_selffac_var_1114(iplon_var_1133, i_lay_var_1130) * (selfrefc_var_251(inds_var_1127, ig_var_1124) + p_selffrac_var_1115(iplon_var_1133, i_lay_var_1130) * (selfrefc_var_251(inds_var_1127 + 1, ig_var_1124) - selfrefc_var_251(inds_var_1127, ig_var_1124))) + p_forfac_var_1117(iplon_var_1133, i_lay_var_1130) * (forrefc_var_252(indf_var_1128, ig_var_1124) + p_forfrac_var_1118(iplon_var_1133, i_lay_var_1130) * (forrefc_var_252(indf_var_1128 + 1, ig_var_1124) - forrefc_var_252(indf_var_1128, ig_var_1124))))
        p_taur_var_1122(iplon_var_1133, i_lay_var_1130, ig_var_1124) = z_tauray_var_1140
      END DO
    END DO
  END DO
  DO i_lay_var_1130 = laytrop_min_var_1134 + 1, laytrop_max_var_1135
    DO iplon_var_1133 = kidia_var_1099, kfdia_var_1100
      IF (i_lay_var_1130 <= k_laytrop_var_1113(iplon_var_1133)) THEN
        z_speccomb_var_1137 = p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) + strrat_var_247 * p_colco2_var_1111(iplon_var_1133, i_lay_var_1130)
        z_specparm_var_1139 = p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) / z_speccomb_var_1137
        z_specparm_var_1139 = MIN(p_oneminus_var_1109(iplon_var_1133), z_specparm_var_1139)
        z_specmult_var_1138 = 8.0D0 * (z_specparm_var_1139)
        js_var_1129 = 1 + INT(z_specmult_var_1138)
        z_fs_var_1136 = z_specmult_var_1138 - AINT(z_specmult_var_1138)
        ind0_var_1125 = ((k_jp_var_1106(iplon_var_1133, i_lay_var_1130) - 1) * 5 + (k_jt_var_1107(iplon_var_1133, i_lay_var_1130) - 1)) * nspa_var_334(17) + js_var_1129
        ind1_var_1126 = (k_jp_var_1106(iplon_var_1133, i_lay_var_1130) * 5 + (k_jt1_var_1108(iplon_var_1133, i_lay_var_1130) - 1)) * nspa_var_334(17) + js_var_1129
        inds_var_1127 = k_indself_var_1116(iplon_var_1133, i_lay_var_1130)
        indf_var_1128 = k_indfor_var_1119(iplon_var_1133, i_lay_var_1130)
        z_tauray_var_1140 = p_colmol_var_1112(iplon_var_1133, i_lay_var_1130) * rayl_var_246
        DO ig_var_1124 = 1, 12
          p_taug_var_1121(iplon_var_1133, i_lay_var_1130, ig_var_1124) = z_speccomb_var_1137 * ((1.0D0 - z_fs_var_1136) * (absa_var_249(ind0_var_1125, ig_var_1124) * p_fac00_var_1102(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind0_var_1125 + 9, ig_var_1124) * p_fac10_var_1104(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind1_var_1126, ig_var_1124) * p_fac01_var_1103(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind1_var_1126 + 9, ig_var_1124) * p_fac11_var_1105(iplon_var_1133, i_lay_var_1130)) + z_fs_var_1136 * (absa_var_249(ind0_var_1125 + 1, ig_var_1124) * p_fac00_var_1102(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind0_var_1125 + 10, ig_var_1124) * p_fac10_var_1104(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind1_var_1126 + 1, ig_var_1124) * p_fac01_var_1103(iplon_var_1133, i_lay_var_1130) + absa_var_249(ind1_var_1126 + 10, ig_var_1124) * p_fac11_var_1105(iplon_var_1133, i_lay_var_1130))) + p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) * (p_selffac_var_1114(iplon_var_1133, i_lay_var_1130) * (selfrefc_var_251(inds_var_1127, ig_var_1124) + p_selffrac_var_1115(iplon_var_1133, i_lay_var_1130) * (selfrefc_var_251(inds_var_1127 + 1, ig_var_1124) - selfrefc_var_251(inds_var_1127, ig_var_1124))) + p_forfac_var_1117(iplon_var_1133, i_lay_var_1130) * (forrefc_var_252(indf_var_1128, ig_var_1124) + p_forfrac_var_1118(iplon_var_1133, i_lay_var_1130) * (forrefc_var_252(indf_var_1128 + 1, ig_var_1124) - forrefc_var_252(indf_var_1128, ig_var_1124))))
          p_taur_var_1122(iplon_var_1133, i_lay_var_1130, ig_var_1124) = z_tauray_var_1140
        END DO
      ELSE
        IF (k_jp_var_1106(iplon_var_1133, i_lay_var_1130 - 1) < layreffr_var_248 .AND. k_jp_var_1106(iplon_var_1133, i_lay_var_1130) >= layreffr_var_248) i_laysolfr_var_1131(iplon_var_1133) = i_lay_var_1130
        z_speccomb_var_1137 = p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) + strrat_var_247 * p_colco2_var_1111(iplon_var_1133, i_lay_var_1130)
        z_specparm_var_1139 = p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) / z_speccomb_var_1137
        z_specparm_var_1139 = MIN(p_oneminus_var_1109(iplon_var_1133), z_specparm_var_1139)
        z_specmult_var_1138 = 4.0D0 * (z_specparm_var_1139)
        js_var_1129 = 1 + INT(z_specmult_var_1138)
        z_fs_var_1136 = z_specmult_var_1138 - AINT(z_specmult_var_1138)
        ind0_var_1125 = ((k_jp_var_1106(iplon_var_1133, i_lay_var_1130) - 13) * 5 + (k_jt_var_1107(iplon_var_1133, i_lay_var_1130) - 1)) * nspb_var_335(17) + js_var_1129
        ind1_var_1126 = ((k_jp_var_1106(iplon_var_1133, i_lay_var_1130) - 12) * 5 + (k_jt1_var_1108(iplon_var_1133, i_lay_var_1130) - 1)) * nspb_var_335(17) + js_var_1129
        indf_var_1128 = k_indfor_var_1119(iplon_var_1133, i_lay_var_1130)
        z_tauray_var_1140 = p_colmol_var_1112(iplon_var_1133, i_lay_var_1130) * rayl_var_246
        DO ig_var_1124 = 1, 12
          p_taug_var_1121(iplon_var_1133, i_lay_var_1130, ig_var_1124) = z_speccomb_var_1137 * ((1.0D0 - z_fs_var_1136) * (absb_var_250(ind0_var_1125, ig_var_1124) * p_fac00_var_1102(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind0_var_1125 + 5, ig_var_1124) * p_fac10_var_1104(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind1_var_1126, ig_var_1124) * p_fac01_var_1103(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind1_var_1126 + 5, ig_var_1124) * p_fac11_var_1105(iplon_var_1133, i_lay_var_1130)) + z_fs_var_1136 * (absb_var_250(ind0_var_1125 + 1, ig_var_1124) * p_fac00_var_1102(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind0_var_1125 + 6, ig_var_1124) * p_fac10_var_1104(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind1_var_1126 + 1, ig_var_1124) * p_fac01_var_1103(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind1_var_1126 + 6, ig_var_1124) * p_fac11_var_1105(iplon_var_1133, i_lay_var_1130))) + p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) * p_forfac_var_1117(iplon_var_1133, i_lay_var_1130) * (forrefc_var_252(indf_var_1128, ig_var_1124) + p_forfrac_var_1118(iplon_var_1133, i_lay_var_1130) * (forrefc_var_252(indf_var_1128 + 1, ig_var_1124) - forrefc_var_252(indf_var_1128, ig_var_1124)))
          IF (i_lay_var_1130 == i_laysolfr_var_1131(iplon_var_1133)) p_sfluxzen_var_1120(iplon_var_1133, ig_var_1124) = sfluxrefc_var_253(ig_var_1124, js_var_1129) + z_fs_var_1136 * (sfluxrefc_var_253(ig_var_1124, js_var_1129 + 1) - sfluxrefc_var_253(ig_var_1124, js_var_1129))
          p_taur_var_1122(iplon_var_1133, i_lay_var_1130, ig_var_1124) = z_tauray_var_1140
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1130 = laytrop_max_var_1135 + 1, i_nlayers_var_1132
    DO iplon_var_1133 = kidia_var_1099, kfdia_var_1100
      IF (k_jp_var_1106(iplon_var_1133, i_lay_var_1130 - 1) < layreffr_var_248 .AND. k_jp_var_1106(iplon_var_1133, i_lay_var_1130) >= layreffr_var_248) i_laysolfr_var_1131(iplon_var_1133) = i_lay_var_1130
      z_speccomb_var_1137 = p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) + strrat_var_247 * p_colco2_var_1111(iplon_var_1133, i_lay_var_1130)
      z_specparm_var_1139 = p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) / z_speccomb_var_1137
      z_specparm_var_1139 = MIN(p_oneminus_var_1109(iplon_var_1133), z_specparm_var_1139)
      z_specmult_var_1138 = 4.0D0 * (z_specparm_var_1139)
      js_var_1129 = 1 + INT(z_specmult_var_1138)
      z_fs_var_1136 = z_specmult_var_1138 - AINT(z_specmult_var_1138)
      ind0_var_1125 = ((k_jp_var_1106(iplon_var_1133, i_lay_var_1130) - 13) * 5 + (k_jt_var_1107(iplon_var_1133, i_lay_var_1130) - 1)) * nspb_var_335(17) + js_var_1129
      ind1_var_1126 = ((k_jp_var_1106(iplon_var_1133, i_lay_var_1130) - 12) * 5 + (k_jt1_var_1108(iplon_var_1133, i_lay_var_1130) - 1)) * nspb_var_335(17) + js_var_1129
      indf_var_1128 = k_indfor_var_1119(iplon_var_1133, i_lay_var_1130)
      z_tauray_var_1140 = p_colmol_var_1112(iplon_var_1133, i_lay_var_1130) * rayl_var_246
      DO ig_var_1124 = 1, 12
        p_taug_var_1121(iplon_var_1133, i_lay_var_1130, ig_var_1124) = z_speccomb_var_1137 * ((1.0D0 - z_fs_var_1136) * (absb_var_250(ind0_var_1125, ig_var_1124) * p_fac00_var_1102(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind0_var_1125 + 5, ig_var_1124) * p_fac10_var_1104(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind1_var_1126, ig_var_1124) * p_fac01_var_1103(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind1_var_1126 + 5, ig_var_1124) * p_fac11_var_1105(iplon_var_1133, i_lay_var_1130)) + z_fs_var_1136 * (absb_var_250(ind0_var_1125 + 1, ig_var_1124) * p_fac00_var_1102(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind0_var_1125 + 6, ig_var_1124) * p_fac10_var_1104(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind1_var_1126 + 1, ig_var_1124) * p_fac01_var_1103(iplon_var_1133, i_lay_var_1130) + absb_var_250(ind1_var_1126 + 6, ig_var_1124) * p_fac11_var_1105(iplon_var_1133, i_lay_var_1130))) + p_colh2o_var_1110(iplon_var_1133, i_lay_var_1130) * p_forfac_var_1117(iplon_var_1133, i_lay_var_1130) * (forrefc_var_252(indf_var_1128, ig_var_1124) + p_forfrac_var_1118(iplon_var_1133, i_lay_var_1130) * (forrefc_var_252(indf_var_1128 + 1, ig_var_1124) - forrefc_var_252(indf_var_1128, ig_var_1124)))
        IF (i_lay_var_1130 == i_laysolfr_var_1131(iplon_var_1133)) p_sfluxzen_var_1120(iplon_var_1133, ig_var_1124) = sfluxrefc_var_253(ig_var_1124, js_var_1129) + z_fs_var_1136 * (sfluxrefc_var_253(ig_var_1124, js_var_1129 + 1) - sfluxrefc_var_253(ig_var_1124, js_var_1129))
        p_taur_var_1122(iplon_var_1133, i_lay_var_1130, ig_var_1124) = z_tauray_var_1140
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol17
SUBROUTINE rrtm_setcoef_140gp(kidia_var_1141, kfdia_var_1142, klev_var_1143, p_coldry, p_wbroad, p_wkl, p_fac00_var_1144, p_fac01_var_1145, p_fac10_var_1146, p_fac11_var_1147, p_forfac_var_1148, p_forfrac_var_1149, k_indfor_var_1166, k_jp_var_1150, k_jt_var_1151, k_jt1_var_1152, p_colh2o_var_1153, p_colco2_var_1154, p_colo3_var_1155, p_coln2o, p_colch4_var_1156, p_colo2_var_1157, p_co2mult_var_1158, p_colbrd, k_laytrop_var_1159, k_layswtch_var_1160, k_laylow_var_1161, pavel_var_1162, p_tavel, p_selffac_var_1163, p_selffrac_var_1164, k_indself_var_1165, k_indminor, p_scaleminor, p_scaleminorn2, p_minorfrac, prat_h2oco2_var_1167, prat_h2oco2_1_var_1168, prat_h2oo3_var_1169, prat_h2oo3_1_var_1170, prat_h2on2o_var_1171, prat_h2on2o_1_var_1172, prat_h2och4_var_1173, prat_h2och4_1_var_1174, prat_n2oco2_var_1175, prat_n2oco2_1_var_1176, prat_o3co2_var_1177, prat_o3co2_1_var_1178)
  USE yoerrtrf, ONLY: chi_mls, preflog_var_235, tref_var_236
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1141
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1142
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1143
  REAL(KIND = 8), INTENT(IN) :: p_coldry(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(IN) :: p_wbroad(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_colbrd(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(IN) :: p_wkl(kidia_var_1141 : kfdia_var_1142, 35, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_fac00_var_1144(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_fac01_var_1145(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_fac10_var_1146(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_fac11_var_1147(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_forfac_var_1148(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_forfrac_var_1149(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jp_var_1150(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jt_var_1151(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  INTEGER(KIND = 4), INTENT(OUT) :: k_jt1_var_1152(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_colh2o_var_1153(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_colco2_var_1154(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_colo3_var_1155(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_coln2o(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_colch4_var_1156(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_colo2_var_1157(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_co2mult_var_1158(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  INTEGER(KIND = 4), INTENT(OUT) :: k_laytrop_var_1159(kidia_var_1141 : kfdia_var_1142)
  INTEGER(KIND = 4), INTENT(OUT) :: k_layswtch_var_1160(kidia_var_1141 : kfdia_var_1142)
  INTEGER(KIND = 4), INTENT(OUT) :: k_laylow_var_1161(kidia_var_1141 : kfdia_var_1142)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_1162(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(IN) :: p_tavel(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_selffac_var_1163(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_selffrac_var_1164(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indself_var_1165(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indfor_var_1166(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  INTEGER(KIND = 4), INTENT(OUT) :: k_indminor(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_scaleminor(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_scaleminorn2(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: p_minorfrac(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  REAL(KIND = 8), INTENT(OUT) :: prat_h2oco2_var_1167(kidia_var_1141 : kfdia_var_1142, klev_var_1143), prat_h2oco2_1_var_1168(kidia_var_1141 : kfdia_var_1142, klev_var_1143), prat_h2oo3_var_1169(kidia_var_1141 : kfdia_var_1142, klev_var_1143), prat_h2oo3_1_var_1170(kidia_var_1141 : kfdia_var_1142, klev_var_1143), prat_h2on2o_var_1171(kidia_var_1141 : kfdia_var_1142, klev_var_1143), prat_h2on2o_1_var_1172(kidia_var_1141 : kfdia_var_1142, klev_var_1143), prat_h2och4_var_1173(kidia_var_1141 : kfdia_var_1142, klev_var_1143), prat_h2och4_1_var_1174(kidia_var_1141 : kfdia_var_1142, klev_var_1143), prat_n2oco2_var_1175(kidia_var_1141 : kfdia_var_1142, klev_var_1143), prat_n2oco2_1_var_1176(kidia_var_1141 : kfdia_var_1142, klev_var_1143), prat_o3co2_var_1177(kidia_var_1141 : kfdia_var_1142, klev_var_1143), prat_o3co2_1_var_1178(kidia_var_1141 : kfdia_var_1142, klev_var_1143)
  INTEGER(KIND = 4) :: jp1_var_1179, jlay
  INTEGER(KIND = 4) :: jlon_var_1180
  REAL(KIND = 8) :: z_co2reg_var_1181, z_compfp_var_1182, z_factor_var_1183, z_fp_var_1184, z_ft_var_1185, z_ft1_var_1186, z_plog_var_1187, z_scalefac_var_1188, z_stpfac_var_1189, z_water_var_1190
  DO jlon_var_1180 = kidia_var_1141, kfdia_var_1142
    z_stpfac_var_1189 = 0.29220138203356366D0
    k_laytrop_var_1159(jlon_var_1180) = 0
    k_layswtch_var_1160(jlon_var_1180) = 0
    k_laylow_var_1161(jlon_var_1180) = 0
    DO jlay = 1, klev_var_1143
      z_plog_var_1187 = LOG(pavel_var_1162(jlon_var_1180, jlay))
      k_jp_var_1150(jlon_var_1180, jlay) = INT(36.0D0 - 5 * (z_plog_var_1187 + 0.04D0))
      IF (k_jp_var_1150(jlon_var_1180, jlay) < 1) THEN
        k_jp_var_1150(jlon_var_1180, jlay) = 1
      ELSE IF (k_jp_var_1150(jlon_var_1180, jlay) > 58) THEN
        k_jp_var_1150(jlon_var_1180, jlay) = 58
      END IF
      jp1_var_1179 = k_jp_var_1150(jlon_var_1180, jlay) + 1
      z_fp_var_1184 = 5.0D0 * (preflog_var_235(k_jp_var_1150(jlon_var_1180, jlay)) - z_plog_var_1187)
      z_fp_var_1184 = MAX(- 1.0D0, MIN(1.0D0, z_fp_var_1184))
      k_jt_var_1151(jlon_var_1180, jlay) = INT(3.0D0 + (p_tavel(jlon_var_1180, jlay) - tref_var_236(k_jp_var_1150(jlon_var_1180, jlay))) / 15.0D0)
      IF (k_jt_var_1151(jlon_var_1180, jlay) < 1) THEN
        k_jt_var_1151(jlon_var_1180, jlay) = 1
      ELSE IF (k_jt_var_1151(jlon_var_1180, jlay) > 4) THEN
        k_jt_var_1151(jlon_var_1180, jlay) = 4
      END IF
      z_ft_var_1185 = ((p_tavel(jlon_var_1180, jlay) - tref_var_236(k_jp_var_1150(jlon_var_1180, jlay))) / 15.0D0) - REAL(k_jt_var_1151(jlon_var_1180, jlay) - 3)
      k_jt1_var_1152(jlon_var_1180, jlay) = INT(3.0D0 + (p_tavel(jlon_var_1180, jlay) - tref_var_236(jp1_var_1179)) / 15.0D0)
      IF (k_jt1_var_1152(jlon_var_1180, jlay) < 1) THEN
        k_jt1_var_1152(jlon_var_1180, jlay) = 1
      ELSE IF (k_jt1_var_1152(jlon_var_1180, jlay) > 4) THEN
        k_jt1_var_1152(jlon_var_1180, jlay) = 4
      END IF
      z_ft1_var_1186 = ((p_tavel(jlon_var_1180, jlay) - tref_var_236(jp1_var_1179)) / 15.0D0) - REAL(k_jt1_var_1152(jlon_var_1180, jlay) - 3)
      z_water_var_1190 = p_wkl(jlon_var_1180, 1, jlay) / p_coldry(jlon_var_1180, jlay)
      z_scalefac_var_1188 = pavel_var_1162(jlon_var_1180, jlay) * z_stpfac_var_1189 / p_tavel(jlon_var_1180, jlay)
      IF (z_plog_var_1187 > 4.56D0) THEN
        k_laytrop_var_1159(jlon_var_1180) = k_laytrop_var_1159(jlon_var_1180) + 1
        p_forfac_var_1148(jlon_var_1180, jlay) = z_scalefac_var_1188 / (1.0D0 + z_water_var_1190)
        z_factor_var_1183 = (332.0D0 - p_tavel(jlon_var_1180, jlay)) / 36.0D0
        k_indfor_var_1166(jlon_var_1180, jlay) = MIN(2, MAX(1, INT(z_factor_var_1183)))
        p_forfrac_var_1149(jlon_var_1180, jlay) = z_factor_var_1183 - REAL(k_indfor_var_1166(jlon_var_1180, jlay))
        p_selffac_var_1163(jlon_var_1180, jlay) = z_water_var_1190 * p_forfac_var_1148(jlon_var_1180, jlay)
        z_factor_var_1183 = (p_tavel(jlon_var_1180, jlay) - 188.0D0) / 7.2D0
        k_indself_var_1165(jlon_var_1180, jlay) = MIN(9, MAX(1, INT(z_factor_var_1183) - 7))
        p_selffrac_var_1164(jlon_var_1180, jlay) = z_factor_var_1183 - REAL(k_indself_var_1165(jlon_var_1180, jlay) + 7)
        p_scaleminor(jlon_var_1180, jlay) = pavel_var_1162(jlon_var_1180, jlay) / p_tavel(jlon_var_1180, jlay)
        p_scaleminorn2(jlon_var_1180, jlay) = (pavel_var_1162(jlon_var_1180, jlay) / p_tavel(jlon_var_1180, jlay)) * (p_wbroad(jlon_var_1180, jlay) / (p_coldry(jlon_var_1180, jlay) + p_wkl(jlon_var_1180, 1, jlay)))
        z_factor_var_1183 = (p_tavel(jlon_var_1180, jlay) - 180.8D0) / 7.2D0
        k_indminor(jlon_var_1180, jlay) = MIN(18, MAX(1, INT(z_factor_var_1183)))
        p_minorfrac(jlon_var_1180, jlay) = z_factor_var_1183 - REAL(k_indminor(jlon_var_1180, jlay))
        prat_h2oco2_var_1167(jlon_var_1180, jlay) = chi_mls(1, k_jp_var_1150(jlon_var_1180, jlay)) / chi_mls(2, k_jp_var_1150(jlon_var_1180, jlay))
        prat_h2oco2_1_var_1168(jlon_var_1180, jlay) = chi_mls(1, k_jp_var_1150(jlon_var_1180, jlay) + 1) / chi_mls(2, k_jp_var_1150(jlon_var_1180, jlay) + 1)
        prat_h2oo3_var_1169(jlon_var_1180, jlay) = chi_mls(1, k_jp_var_1150(jlon_var_1180, jlay)) / chi_mls(3, k_jp_var_1150(jlon_var_1180, jlay))
        prat_h2oo3_1_var_1170(jlon_var_1180, jlay) = chi_mls(1, k_jp_var_1150(jlon_var_1180, jlay) + 1) / chi_mls(3, k_jp_var_1150(jlon_var_1180, jlay) + 1)
        prat_h2on2o_var_1171(jlon_var_1180, jlay) = chi_mls(1, k_jp_var_1150(jlon_var_1180, jlay)) / chi_mls(4, k_jp_var_1150(jlon_var_1180, jlay))
        prat_h2on2o_1_var_1172(jlon_var_1180, jlay) = chi_mls(1, k_jp_var_1150(jlon_var_1180, jlay) + 1) / chi_mls(4, k_jp_var_1150(jlon_var_1180, jlay) + 1)
        prat_h2och4_var_1173(jlon_var_1180, jlay) = chi_mls(1, k_jp_var_1150(jlon_var_1180, jlay)) / chi_mls(6, k_jp_var_1150(jlon_var_1180, jlay))
        prat_h2och4_1_var_1174(jlon_var_1180, jlay) = chi_mls(1, k_jp_var_1150(jlon_var_1180, jlay) + 1) / chi_mls(6, k_jp_var_1150(jlon_var_1180, jlay) + 1)
        prat_n2oco2_var_1175(jlon_var_1180, jlay) = chi_mls(4, k_jp_var_1150(jlon_var_1180, jlay)) / chi_mls(2, k_jp_var_1150(jlon_var_1180, jlay))
        prat_n2oco2_1_var_1176(jlon_var_1180, jlay) = chi_mls(4, k_jp_var_1150(jlon_var_1180, jlay) + 1) / chi_mls(2, k_jp_var_1150(jlon_var_1180, jlay) + 1)
        p_colh2o_var_1153(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 1, jlay)
        p_colco2_var_1154(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 2, jlay)
        p_colo3_var_1155(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 3, jlay)
        p_coln2o(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 4, jlay)
        p_colch4_var_1156(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 6, jlay)
        p_colo2_var_1157(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 7, jlay)
        p_colbrd(jlon_var_1180, jlay) = 1D-20 * p_wbroad(jlon_var_1180, jlay)
        IF (p_colco2_var_1154(jlon_var_1180, jlay) == 0.0D0) p_colco2_var_1154(jlon_var_1180, jlay) = 1D-32 * p_coldry(jlon_var_1180, jlay)
        IF (p_coln2o(jlon_var_1180, jlay) == 0.0D0) p_coln2o(jlon_var_1180, jlay) = 1D-32 * p_coldry(jlon_var_1180, jlay)
        IF (p_colch4_var_1156(jlon_var_1180, jlay) == 0.0D0) p_colch4_var_1156(jlon_var_1180, jlay) = 1D-32 * p_coldry(jlon_var_1180, jlay)
        z_co2reg_var_1181 = 3.55D-24 * p_coldry(jlon_var_1180, jlay)
        p_co2mult_var_1158(jlon_var_1180, jlay) = (p_colco2_var_1154(jlon_var_1180, jlay) - z_co2reg_var_1181) * 272.63D0 * EXP(- 1919.4D0 / p_tavel(jlon_var_1180, jlay)) / (0.00087604D0 * p_tavel(jlon_var_1180, jlay))
      ELSE
        p_forfac_var_1148(jlon_var_1180, jlay) = z_scalefac_var_1188 / (1.0D0 + z_water_var_1190)
        z_factor_var_1183 = (p_tavel(jlon_var_1180, jlay) - 188.0D0) / 36.0D0
        k_indfor_var_1166(jlon_var_1180, jlay) = 3
        p_forfrac_var_1149(jlon_var_1180, jlay) = z_factor_var_1183 - 1.0D0
        p_selffac_var_1163(jlon_var_1180, jlay) = z_water_var_1190 * p_forfac_var_1148(jlon_var_1180, jlay)
        p_scaleminor(jlon_var_1180, jlay) = pavel_var_1162(jlon_var_1180, jlay) / p_tavel(jlon_var_1180, jlay)
        p_scaleminorn2(jlon_var_1180, jlay) = (pavel_var_1162(jlon_var_1180, jlay) / p_tavel(jlon_var_1180, jlay)) * (p_wbroad(jlon_var_1180, jlay) / (p_coldry(jlon_var_1180, jlay) + p_wkl(jlon_var_1180, 1, jlay)))
        z_factor_var_1183 = (p_tavel(jlon_var_1180, jlay) - 180.8D0) / 7.2D0
        k_indminor(jlon_var_1180, jlay) = MIN(18, MAX(1, INT(z_factor_var_1183)))
        p_minorfrac(jlon_var_1180, jlay) = z_factor_var_1183 - REAL(k_indminor(jlon_var_1180, jlay))
        prat_h2oco2_var_1167(jlon_var_1180, jlay) = chi_mls(1, k_jp_var_1150(jlon_var_1180, jlay)) / chi_mls(2, k_jp_var_1150(jlon_var_1180, jlay))
        prat_h2oco2_1_var_1168(jlon_var_1180, jlay) = chi_mls(1, k_jp_var_1150(jlon_var_1180, jlay) + 1) / chi_mls(2, k_jp_var_1150(jlon_var_1180, jlay) + 1)
        prat_o3co2_var_1177(jlon_var_1180, jlay) = chi_mls(3, k_jp_var_1150(jlon_var_1180, jlay)) / chi_mls(2, k_jp_var_1150(jlon_var_1180, jlay))
        prat_o3co2_1_var_1178(jlon_var_1180, jlay) = chi_mls(3, k_jp_var_1150(jlon_var_1180, jlay) + 1) / chi_mls(2, k_jp_var_1150(jlon_var_1180, jlay) + 1)
        p_colh2o_var_1153(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 1, jlay)
        p_colco2_var_1154(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 2, jlay)
        p_colo3_var_1155(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 3, jlay)
        p_coln2o(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 4, jlay)
        p_colch4_var_1156(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 6, jlay)
        p_colo2_var_1157(jlon_var_1180, jlay) = 1D-20 * p_wkl(jlon_var_1180, 7, jlay)
        p_colbrd(jlon_var_1180, jlay) = 1D-20 * p_wbroad(jlon_var_1180, jlay)
        IF (p_colco2_var_1154(jlon_var_1180, jlay) == 0.0D0) p_colco2_var_1154(jlon_var_1180, jlay) = 1D-32 * p_coldry(jlon_var_1180, jlay)
        IF (p_coln2o(jlon_var_1180, jlay) == 0.0D0) p_coln2o(jlon_var_1180, jlay) = 1D-32 * p_coldry(jlon_var_1180, jlay)
        IF (p_colch4_var_1156(jlon_var_1180, jlay) == 0.0D0) p_colch4_var_1156(jlon_var_1180, jlay) = 1D-32 * p_coldry(jlon_var_1180, jlay)
        z_co2reg_var_1181 = 3.55D-24 * p_coldry(jlon_var_1180, jlay)
        p_co2mult_var_1158(jlon_var_1180, jlay) = (p_colco2_var_1154(jlon_var_1180, jlay) - z_co2reg_var_1181) * 272.63D0 * EXP(- 1919.4D0 / p_tavel(jlon_var_1180, jlay)) / (0.00087604D0 * p_tavel(jlon_var_1180, jlay))
      END IF
      z_compfp_var_1182 = 1.0D0 - z_fp_var_1184
      p_fac10_var_1146(jlon_var_1180, jlay) = z_compfp_var_1182 * z_ft_var_1185
      p_fac00_var_1144(jlon_var_1180, jlay) = z_compfp_var_1182 * (1.0D0 - z_ft_var_1185)
      p_fac11_var_1147(jlon_var_1180, jlay) = z_fp_var_1184 * z_ft1_var_1186
      p_fac01_var_1145(jlon_var_1180, jlay) = z_fp_var_1184 * (1.0D0 - z_ft1_var_1186)
      p_selffac_var_1163(jlon_var_1180, jlay) = p_colh2o_var_1153(jlon_var_1180, jlay) * p_selffac_var_1163(jlon_var_1180, jlay)
      p_forfac_var_1148(jlon_var_1180, jlay) = p_colh2o_var_1153(jlon_var_1180, jlay) * p_forfac_var_1148(jlon_var_1180, jlay)
    END DO
    IF (k_laylow_var_1161(jlon_var_1180) == 0) k_laylow_var_1161(jlon_var_1180) = 1
  END DO
END SUBROUTINE rrtm_setcoef_140gp
SUBROUTINE srtm_taumol16(kidia_var_1191, kfdia_var_1192, klev_var_1193, p_fac00_var_1194, p_fac01_var_1195, p_fac10_var_1196, p_fac11_var_1197, k_jp_var_1198, k_jt_var_1199, k_jt1_var_1200, p_oneminus_var_1201, p_colh2o_var_1202, p_colch4_var_1203, p_colmol_var_1204, k_laytrop_var_1205, p_selffac_var_1206, p_selffrac_var_1207, k_indself_var_1208, p_forfac_var_1209, p_forfrac_var_1210, k_indfor_var_1211, p_sfluxzen_var_1212, p_taug_var_1213, p_taur_var_1214, prmu0_var_1215)
  USE yoesrta16, ONLY: absa_var_241, absb_var_242, forrefc_var_244, layreffr_var_240, rayl_var_239, selfrefc_var_243, sfluxrefc_var_245, strrat1
  USE yoesrtwn, ONLY: nspa_var_334, nspb_var_335
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1191, kfdia_var_1192
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1193
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1194(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1195(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1196(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1197(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1198(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1199(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1200(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1201(kidia_var_1191 : kfdia_var_1192)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1202(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_1203(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1204(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1205(kidia_var_1191 : kfdia_var_1192)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1206(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1207(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1208(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1209(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1210(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1211(kidia_var_1191 : kfdia_var_1192, klev_var_1193)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1212(kidia_var_1191 : kfdia_var_1192, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1213(kidia_var_1191 : kfdia_var_1192, klev_var_1193, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1214(kidia_var_1191 : kfdia_var_1192, klev_var_1193, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1215(kidia_var_1191 : kfdia_var_1192)
  INTEGER(KIND = 4) :: ig_var_1216, ind0_var_1217, ind1_var_1218, inds_var_1219, indf_var_1220, js_var_1221, i_lay_var_1222, i_laysolfr_var_1223(kidia_var_1191 : kfdia_var_1192), i_nlayers_var_1224, iplon_var_1225
  INTEGER(KIND = 4) :: laytrop_min_var_1226, laytrop_max_var_1227
  REAL(KIND = 8) :: z_fs_var_1228, z_speccomb_var_1229, z_specmult_var_1230, z_specparm_var_1231, z_tauray_var_1232
  laytrop_min_var_1226 = MINVAL(k_laytrop_var_1205(kidia_var_1191 : kfdia_var_1192))
  laytrop_max_var_1227 = MAXVAL(k_laytrop_var_1205(kidia_var_1191 : kfdia_var_1192))
  i_nlayers_var_1224 = klev_var_1193
  DO iplon_var_1225 = kidia_var_1191, kfdia_var_1192
    i_laysolfr_var_1223(iplon_var_1225) = i_nlayers_var_1224
  END DO
  DO i_lay_var_1222 = 1, laytrop_min_var_1226
    DO iplon_var_1225 = kidia_var_1191, kfdia_var_1192
      z_speccomb_var_1229 = p_colh2o_var_1202(iplon_var_1225, i_lay_var_1222) + strrat1 * p_colch4_var_1203(iplon_var_1225, i_lay_var_1222)
      z_specparm_var_1231 = p_colh2o_var_1202(iplon_var_1225, i_lay_var_1222) / z_speccomb_var_1229
      z_specparm_var_1231 = MIN(p_oneminus_var_1201(iplon_var_1225), z_specparm_var_1231)
      z_specmult_var_1230 = 8.0D0 * (z_specparm_var_1231)
      js_var_1221 = 1 + INT(z_specmult_var_1230)
      z_fs_var_1228 = z_specmult_var_1230 - AINT(z_specmult_var_1230)
      ind0_var_1217 = ((k_jp_var_1198(iplon_var_1225, i_lay_var_1222) - 1) * 5 + (k_jt_var_1199(iplon_var_1225, i_lay_var_1222) - 1)) * nspa_var_334(16) + js_var_1221
      ind1_var_1218 = (k_jp_var_1198(iplon_var_1225, i_lay_var_1222) * 5 + (k_jt1_var_1200(iplon_var_1225, i_lay_var_1222) - 1)) * nspa_var_334(16) + js_var_1221
      inds_var_1219 = k_indself_var_1208(iplon_var_1225, i_lay_var_1222)
      indf_var_1220 = k_indfor_var_1211(iplon_var_1225, i_lay_var_1222)
      z_tauray_var_1232 = p_colmol_var_1204(iplon_var_1225, i_lay_var_1222) * rayl_var_239
      DO ig_var_1216 = 1, 6
        p_taug_var_1213(iplon_var_1225, i_lay_var_1222, ig_var_1216) = z_speccomb_var_1229 * ((1.0D0 - z_fs_var_1228) * (absa_var_241(ind0_var_1217, ig_var_1216) * p_fac00_var_1194(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind0_var_1217 + 9, ig_var_1216) * p_fac10_var_1196(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind1_var_1218, ig_var_1216) * p_fac01_var_1195(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind1_var_1218 + 9, ig_var_1216) * p_fac11_var_1197(iplon_var_1225, i_lay_var_1222)) + z_fs_var_1228 * (absa_var_241(ind0_var_1217 + 1, ig_var_1216) * p_fac00_var_1194(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind0_var_1217 + 10, ig_var_1216) * p_fac10_var_1196(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind1_var_1218 + 1, ig_var_1216) * p_fac01_var_1195(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind1_var_1218 + 10, ig_var_1216) * p_fac11_var_1197(iplon_var_1225, i_lay_var_1222))) + p_colh2o_var_1202(iplon_var_1225, i_lay_var_1222) * (p_selffac_var_1206(iplon_var_1225, i_lay_var_1222) * (selfrefc_var_243(inds_var_1219, ig_var_1216) + p_selffrac_var_1207(iplon_var_1225, i_lay_var_1222) * (selfrefc_var_243(inds_var_1219 + 1, ig_var_1216) - selfrefc_var_243(inds_var_1219, ig_var_1216))) + p_forfac_var_1209(iplon_var_1225, i_lay_var_1222) * (forrefc_var_244(indf_var_1220, ig_var_1216) + p_forfrac_var_1210(iplon_var_1225, i_lay_var_1222) * (forrefc_var_244(indf_var_1220 + 1, ig_var_1216) - forrefc_var_244(indf_var_1220, ig_var_1216))))
        p_taur_var_1214(iplon_var_1225, i_lay_var_1222, ig_var_1216) = z_tauray_var_1232
      END DO
    END DO
  END DO
  DO i_lay_var_1222 = laytrop_min_var_1226 + 1, laytrop_max_var_1227
    DO iplon_var_1225 = kidia_var_1191, kfdia_var_1192
      IF (i_lay_var_1222 <= k_laytrop_var_1205(iplon_var_1225)) THEN
        z_speccomb_var_1229 = p_colh2o_var_1202(iplon_var_1225, i_lay_var_1222) + strrat1 * p_colch4_var_1203(iplon_var_1225, i_lay_var_1222)
        z_specparm_var_1231 = p_colh2o_var_1202(iplon_var_1225, i_lay_var_1222) / z_speccomb_var_1229
        z_specparm_var_1231 = MIN(p_oneminus_var_1201(iplon_var_1225), z_specparm_var_1231)
        z_specmult_var_1230 = 8.0D0 * (z_specparm_var_1231)
        js_var_1221 = 1 + INT(z_specmult_var_1230)
        z_fs_var_1228 = z_specmult_var_1230 - AINT(z_specmult_var_1230)
        ind0_var_1217 = ((k_jp_var_1198(iplon_var_1225, i_lay_var_1222) - 1) * 5 + (k_jt_var_1199(iplon_var_1225, i_lay_var_1222) - 1)) * nspa_var_334(16) + js_var_1221
        ind1_var_1218 = (k_jp_var_1198(iplon_var_1225, i_lay_var_1222) * 5 + (k_jt1_var_1200(iplon_var_1225, i_lay_var_1222) - 1)) * nspa_var_334(16) + js_var_1221
        inds_var_1219 = k_indself_var_1208(iplon_var_1225, i_lay_var_1222)
        indf_var_1220 = k_indfor_var_1211(iplon_var_1225, i_lay_var_1222)
        z_tauray_var_1232 = p_colmol_var_1204(iplon_var_1225, i_lay_var_1222) * rayl_var_239
        DO ig_var_1216 = 1, 6
          p_taug_var_1213(iplon_var_1225, i_lay_var_1222, ig_var_1216) = z_speccomb_var_1229 * ((1.0D0 - z_fs_var_1228) * (absa_var_241(ind0_var_1217, ig_var_1216) * p_fac00_var_1194(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind0_var_1217 + 9, ig_var_1216) * p_fac10_var_1196(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind1_var_1218, ig_var_1216) * p_fac01_var_1195(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind1_var_1218 + 9, ig_var_1216) * p_fac11_var_1197(iplon_var_1225, i_lay_var_1222)) + z_fs_var_1228 * (absa_var_241(ind0_var_1217 + 1, ig_var_1216) * p_fac00_var_1194(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind0_var_1217 + 10, ig_var_1216) * p_fac10_var_1196(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind1_var_1218 + 1, ig_var_1216) * p_fac01_var_1195(iplon_var_1225, i_lay_var_1222) + absa_var_241(ind1_var_1218 + 10, ig_var_1216) * p_fac11_var_1197(iplon_var_1225, i_lay_var_1222))) + p_colh2o_var_1202(iplon_var_1225, i_lay_var_1222) * (p_selffac_var_1206(iplon_var_1225, i_lay_var_1222) * (selfrefc_var_243(inds_var_1219, ig_var_1216) + p_selffrac_var_1207(iplon_var_1225, i_lay_var_1222) * (selfrefc_var_243(inds_var_1219 + 1, ig_var_1216) - selfrefc_var_243(inds_var_1219, ig_var_1216))) + p_forfac_var_1209(iplon_var_1225, i_lay_var_1222) * (forrefc_var_244(indf_var_1220, ig_var_1216) + p_forfrac_var_1210(iplon_var_1225, i_lay_var_1222) * (forrefc_var_244(indf_var_1220 + 1, ig_var_1216) - forrefc_var_244(indf_var_1220, ig_var_1216))))
          p_taur_var_1214(iplon_var_1225, i_lay_var_1222, ig_var_1216) = z_tauray_var_1232
        END DO
      ELSE
        IF (k_jp_var_1198(iplon_var_1225, i_lay_var_1222 - 1) < layreffr_var_240 .AND. k_jp_var_1198(iplon_var_1225, i_lay_var_1222) >= layreffr_var_240) i_laysolfr_var_1223(iplon_var_1225) = i_lay_var_1222
        ind0_var_1217 = ((k_jp_var_1198(iplon_var_1225, i_lay_var_1222) - 13) * 5 + (k_jt_var_1199(iplon_var_1225, i_lay_var_1222) - 1)) * nspb_var_335(16) + 1
        ind1_var_1218 = ((k_jp_var_1198(iplon_var_1225, i_lay_var_1222) - 12) * 5 + (k_jt1_var_1200(iplon_var_1225, i_lay_var_1222) - 1)) * nspb_var_335(16) + 1
        z_tauray_var_1232 = p_colmol_var_1204(iplon_var_1225, i_lay_var_1222) * rayl_var_239
        DO ig_var_1216 = 1, 6
          p_taug_var_1213(iplon_var_1225, i_lay_var_1222, ig_var_1216) = p_colch4_var_1203(iplon_var_1225, i_lay_var_1222) * (p_fac00_var_1194(iplon_var_1225, i_lay_var_1222) * absb_var_242(ind0_var_1217, ig_var_1216) + p_fac10_var_1196(iplon_var_1225, i_lay_var_1222) * absb_var_242(ind0_var_1217 + 1, ig_var_1216) + p_fac01_var_1195(iplon_var_1225, i_lay_var_1222) * absb_var_242(ind1_var_1218, ig_var_1216) + p_fac11_var_1197(iplon_var_1225, i_lay_var_1222) * absb_var_242(ind1_var_1218 + 1, ig_var_1216))
          IF (i_lay_var_1222 == i_laysolfr_var_1223(iplon_var_1225)) p_sfluxzen_var_1212(iplon_var_1225, ig_var_1216) = sfluxrefc_var_245(ig_var_1216)
          p_taur_var_1214(iplon_var_1225, i_lay_var_1222, ig_var_1216) = z_tauray_var_1232
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1222 = laytrop_max_var_1227 + 1, i_nlayers_var_1224
    DO iplon_var_1225 = kidia_var_1191, kfdia_var_1192
      IF (k_jp_var_1198(iplon_var_1225, i_lay_var_1222 - 1) < layreffr_var_240 .AND. k_jp_var_1198(iplon_var_1225, i_lay_var_1222) >= layreffr_var_240) i_laysolfr_var_1223(iplon_var_1225) = i_lay_var_1222
      ind0_var_1217 = ((k_jp_var_1198(iplon_var_1225, i_lay_var_1222) - 13) * 5 + (k_jt_var_1199(iplon_var_1225, i_lay_var_1222) - 1)) * nspb_var_335(16) + 1
      ind1_var_1218 = ((k_jp_var_1198(iplon_var_1225, i_lay_var_1222) - 12) * 5 + (k_jt1_var_1200(iplon_var_1225, i_lay_var_1222) - 1)) * nspb_var_335(16) + 1
      z_tauray_var_1232 = p_colmol_var_1204(iplon_var_1225, i_lay_var_1222) * rayl_var_239
      DO ig_var_1216 = 1, 6
        p_taug_var_1213(iplon_var_1225, i_lay_var_1222, ig_var_1216) = p_colch4_var_1203(iplon_var_1225, i_lay_var_1222) * (p_fac00_var_1194(iplon_var_1225, i_lay_var_1222) * absb_var_242(ind0_var_1217, ig_var_1216) + p_fac10_var_1196(iplon_var_1225, i_lay_var_1222) * absb_var_242(ind0_var_1217 + 1, ig_var_1216) + p_fac01_var_1195(iplon_var_1225, i_lay_var_1222) * absb_var_242(ind1_var_1218, ig_var_1216) + p_fac11_var_1197(iplon_var_1225, i_lay_var_1222) * absb_var_242(ind1_var_1218 + 1, ig_var_1216))
        IF (i_lay_var_1222 == i_laysolfr_var_1223(iplon_var_1225)) p_sfluxzen_var_1212(iplon_var_1225, ig_var_1216) = sfluxrefc_var_245(ig_var_1216)
        p_taur_var_1214(iplon_var_1225, i_lay_var_1222, ig_var_1216) = z_tauray_var_1232
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol16
SUBROUTINE srtm_taumol27(kidia_var_1233, kfdia_var_1234, klev_var_1235, p_fac00_var_1236, p_fac01_var_1237, p_fac10_var_1238, p_fac11_var_1239, k_jp_var_1240, k_jt_var_1241, k_jt1_var_1242, p_colmol_var_1243, p_colo3_var_1244, k_laytrop_var_1245, p_sfluxzen_var_1246, p_taug_var_1247, p_taur_var_1248, prmu0_var_1249)
  USE yoesrtwn, ONLY: nspa_var_334, nspb_var_335
  USE yoesrta27, ONLY: absa_var_317, absb_var_318, layreffr_var_316, raylc_var_320, scalekur, sfluxrefc_var_319
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1233, kfdia_var_1234
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1235
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1236(kidia_var_1233 : kfdia_var_1234, klev_var_1235)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1237(kidia_var_1233 : kfdia_var_1234, klev_var_1235)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1238(kidia_var_1233 : kfdia_var_1234, klev_var_1235)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1239(kidia_var_1233 : kfdia_var_1234, klev_var_1235)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1240(kidia_var_1233 : kfdia_var_1234, klev_var_1235)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1241(kidia_var_1233 : kfdia_var_1234, klev_var_1235)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1242(kidia_var_1233 : kfdia_var_1234, klev_var_1235)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1243(kidia_var_1233 : kfdia_var_1234, klev_var_1235)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1244(kidia_var_1233 : kfdia_var_1234, klev_var_1235)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1245(kidia_var_1233 : kfdia_var_1234)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1246(kidia_var_1233 : kfdia_var_1234, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1247(kidia_var_1233 : kfdia_var_1234, klev_var_1235, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1248(kidia_var_1233 : kfdia_var_1234, klev_var_1235, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1249(kidia_var_1233 : kfdia_var_1234)
  INTEGER(KIND = 4) :: ig_var_1250, ind0_var_1251, ind1_var_1252, i_lay_var_1253, i_laysolfr_var_1254(kidia_var_1233 : kfdia_var_1234), i_nlayers_var_1255, iplon_var_1256
  INTEGER(KIND = 4) :: laytrop_min_var_1257, laytrop_max_var_1258
  REAL(KIND = 8) :: z_tauray_var_1259
  laytrop_min_var_1257 = MINVAL(k_laytrop_var_1245(kidia_var_1233 : kfdia_var_1234))
  laytrop_max_var_1258 = MAXVAL(k_laytrop_var_1245(kidia_var_1233 : kfdia_var_1234))
  i_nlayers_var_1255 = klev_var_1235
  DO iplon_var_1256 = kidia_var_1233, kfdia_var_1234
    i_laysolfr_var_1254(iplon_var_1256) = i_nlayers_var_1255
  END DO
  DO i_lay_var_1253 = 1, laytrop_min_var_1257
    DO iplon_var_1256 = kidia_var_1233, kfdia_var_1234
      ind0_var_1251 = ((k_jp_var_1240(iplon_var_1256, i_lay_var_1253) - 1) * 5 + (k_jt_var_1241(iplon_var_1256, i_lay_var_1253) - 1)) * nspa_var_334(27) + 1
      ind1_var_1252 = (k_jp_var_1240(iplon_var_1256, i_lay_var_1253) * 5 + (k_jt1_var_1242(iplon_var_1256, i_lay_var_1253) - 1)) * nspa_var_334(27) + 1
      DO ig_var_1250 = 1, 8
        z_tauray_var_1259 = p_colmol_var_1243(iplon_var_1256, i_lay_var_1253) * raylc_var_320(ig_var_1250)
        p_taug_var_1247(iplon_var_1256, i_lay_var_1253, ig_var_1250) = p_colo3_var_1244(iplon_var_1256, i_lay_var_1253) * (p_fac00_var_1236(iplon_var_1256, i_lay_var_1253) * absa_var_317(ind0_var_1251, ig_var_1250) + p_fac10_var_1238(iplon_var_1256, i_lay_var_1253) * absa_var_317(ind0_var_1251 + 1, ig_var_1250) + p_fac01_var_1237(iplon_var_1256, i_lay_var_1253) * absa_var_317(ind1_var_1252, ig_var_1250) + p_fac11_var_1239(iplon_var_1256, i_lay_var_1253) * absa_var_317(ind1_var_1252 + 1, ig_var_1250))
        p_taur_var_1248(iplon_var_1256, i_lay_var_1253, ig_var_1250) = z_tauray_var_1259
      END DO
    END DO
  END DO
  DO i_lay_var_1253 = laytrop_min_var_1257 + 1, laytrop_max_var_1258
    DO iplon_var_1256 = kidia_var_1233, kfdia_var_1234
      IF (i_lay_var_1253 <= k_laytrop_var_1245(iplon_var_1256)) THEN
        ind0_var_1251 = ((k_jp_var_1240(iplon_var_1256, i_lay_var_1253) - 1) * 5 + (k_jt_var_1241(iplon_var_1256, i_lay_var_1253) - 1)) * nspa_var_334(27) + 1
        ind1_var_1252 = (k_jp_var_1240(iplon_var_1256, i_lay_var_1253) * 5 + (k_jt1_var_1242(iplon_var_1256, i_lay_var_1253) - 1)) * nspa_var_334(27) + 1
        DO ig_var_1250 = 1, 8
          z_tauray_var_1259 = p_colmol_var_1243(iplon_var_1256, i_lay_var_1253) * raylc_var_320(ig_var_1250)
          p_taug_var_1247(iplon_var_1256, i_lay_var_1253, ig_var_1250) = p_colo3_var_1244(iplon_var_1256, i_lay_var_1253) * (p_fac00_var_1236(iplon_var_1256, i_lay_var_1253) * absa_var_317(ind0_var_1251, ig_var_1250) + p_fac10_var_1238(iplon_var_1256, i_lay_var_1253) * absa_var_317(ind0_var_1251 + 1, ig_var_1250) + p_fac01_var_1237(iplon_var_1256, i_lay_var_1253) * absa_var_317(ind1_var_1252, ig_var_1250) + p_fac11_var_1239(iplon_var_1256, i_lay_var_1253) * absa_var_317(ind1_var_1252 + 1, ig_var_1250))
          p_taur_var_1248(iplon_var_1256, i_lay_var_1253, ig_var_1250) = z_tauray_var_1259
        END DO
      ELSE
        IF (k_jp_var_1240(iplon_var_1256, i_lay_var_1253 - 1) < layreffr_var_316 .AND. k_jp_var_1240(iplon_var_1256, i_lay_var_1253) >= layreffr_var_316) i_laysolfr_var_1254(iplon_var_1256) = i_lay_var_1253
        ind0_var_1251 = ((k_jp_var_1240(iplon_var_1256, i_lay_var_1253) - 13) * 5 + (k_jt_var_1241(iplon_var_1256, i_lay_var_1253) - 1)) * nspb_var_335(27) + 1
        ind1_var_1252 = ((k_jp_var_1240(iplon_var_1256, i_lay_var_1253) - 12) * 5 + (k_jt1_var_1242(iplon_var_1256, i_lay_var_1253) - 1)) * nspb_var_335(27) + 1
        DO ig_var_1250 = 1, 8
          z_tauray_var_1259 = p_colmol_var_1243(iplon_var_1256, i_lay_var_1253) * raylc_var_320(ig_var_1250)
          p_taug_var_1247(iplon_var_1256, i_lay_var_1253, ig_var_1250) = p_colo3_var_1244(iplon_var_1256, i_lay_var_1253) * (p_fac00_var_1236(iplon_var_1256, i_lay_var_1253) * absb_var_318(ind0_var_1251, ig_var_1250) + p_fac10_var_1238(iplon_var_1256, i_lay_var_1253) * absb_var_318(ind0_var_1251 + 1, ig_var_1250) + p_fac01_var_1237(iplon_var_1256, i_lay_var_1253) * absb_var_318(ind1_var_1252, ig_var_1250) + p_fac11_var_1239(iplon_var_1256, i_lay_var_1253) * absb_var_318(ind1_var_1252 + 1, ig_var_1250))
          IF (i_lay_var_1253 == i_laysolfr_var_1254(iplon_var_1256)) p_sfluxzen_var_1246(iplon_var_1256, ig_var_1250) = scalekur * sfluxrefc_var_319(ig_var_1250)
          p_taur_var_1248(iplon_var_1256, i_lay_var_1253, ig_var_1250) = z_tauray_var_1259
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1253 = laytrop_max_var_1258 + 1, i_nlayers_var_1255
    DO iplon_var_1256 = kidia_var_1233, kfdia_var_1234
      IF (k_jp_var_1240(iplon_var_1256, i_lay_var_1253 - 1) < layreffr_var_316 .AND. k_jp_var_1240(iplon_var_1256, i_lay_var_1253) >= layreffr_var_316) i_laysolfr_var_1254(iplon_var_1256) = i_lay_var_1253
      ind0_var_1251 = ((k_jp_var_1240(iplon_var_1256, i_lay_var_1253) - 13) * 5 + (k_jt_var_1241(iplon_var_1256, i_lay_var_1253) - 1)) * nspb_var_335(27) + 1
      ind1_var_1252 = ((k_jp_var_1240(iplon_var_1256, i_lay_var_1253) - 12) * 5 + (k_jt1_var_1242(iplon_var_1256, i_lay_var_1253) - 1)) * nspb_var_335(27) + 1
      DO ig_var_1250 = 1, 8
        z_tauray_var_1259 = p_colmol_var_1243(iplon_var_1256, i_lay_var_1253) * raylc_var_320(ig_var_1250)
        p_taug_var_1247(iplon_var_1256, i_lay_var_1253, ig_var_1250) = p_colo3_var_1244(iplon_var_1256, i_lay_var_1253) * (p_fac00_var_1236(iplon_var_1256, i_lay_var_1253) * absb_var_318(ind0_var_1251, ig_var_1250) + p_fac10_var_1238(iplon_var_1256, i_lay_var_1253) * absb_var_318(ind0_var_1251 + 1, ig_var_1250) + p_fac01_var_1237(iplon_var_1256, i_lay_var_1253) * absb_var_318(ind1_var_1252, ig_var_1250) + p_fac11_var_1239(iplon_var_1256, i_lay_var_1253) * absb_var_318(ind1_var_1252 + 1, ig_var_1250))
        IF (i_lay_var_1253 == i_laysolfr_var_1254(iplon_var_1256)) p_sfluxzen_var_1246(iplon_var_1256, ig_var_1250) = scalekur * sfluxrefc_var_319(ig_var_1250)
        p_taur_var_1248(iplon_var_1256, i_lay_var_1253, ig_var_1250) = z_tauray_var_1259
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol27
SUBROUTINE srtm_taumol26(kidia_var_1260, kfdia_var_1261, klev_var_1262, p_colmol_var_1263, k_laytrop_var_1264, p_sfluxzen_var_1265, p_taug_var_1266, p_taur_var_1267, prmu0_var_1268)
  USE yoesrta26, ONLY: raylc_var_315, sfluxrefc_var_314
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1260, kfdia_var_1261
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1262
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1263(kidia_var_1260 : kfdia_var_1261, klev_var_1262)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1264(kidia_var_1260 : kfdia_var_1261)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1265(kidia_var_1260 : kfdia_var_1261, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1266(kidia_var_1260 : kfdia_var_1261, klev_var_1262, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1267(kidia_var_1260 : kfdia_var_1261, klev_var_1262, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1268(kidia_var_1260 : kfdia_var_1261)
  INTEGER(KIND = 4) :: ig_var_1269, i_lay_var_1270, i_laysolfr_var_1271(kidia_var_1260 : kfdia_var_1261), i_nlayers_var_1272, iplon_var_1273
  INTEGER(KIND = 4) :: laytrop_min_var_1274, laytrop_max_var_1275
  laytrop_min_var_1274 = MINVAL(k_laytrop_var_1264(kidia_var_1260 : kfdia_var_1261))
  laytrop_max_var_1275 = MAXVAL(k_laytrop_var_1264(kidia_var_1260 : kfdia_var_1261))
  i_nlayers_var_1272 = klev_var_1262
  DO iplon_var_1273 = kidia_var_1260, kfdia_var_1261
    i_laysolfr_var_1271(iplon_var_1273) = k_laytrop_var_1264(iplon_var_1273)
  END DO
  DO i_lay_var_1270 = 1, laytrop_min_var_1274
    DO iplon_var_1273 = kidia_var_1260, kfdia_var_1261
      DO ig_var_1269 = 1, 6
        IF (i_lay_var_1270 == i_laysolfr_var_1271(iplon_var_1273)) p_sfluxzen_var_1265(iplon_var_1273, ig_var_1269) = sfluxrefc_var_314(ig_var_1269)
        p_taug_var_1266(iplon_var_1273, i_lay_var_1270, ig_var_1269) = 0.0D0
        p_taur_var_1267(iplon_var_1273, i_lay_var_1270, ig_var_1269) = p_colmol_var_1263(iplon_var_1273, i_lay_var_1270) * raylc_var_315(ig_var_1269)
      END DO
    END DO
  END DO
  DO i_lay_var_1270 = laytrop_min_var_1274 + 1, laytrop_max_var_1275
    DO iplon_var_1273 = kidia_var_1260, kfdia_var_1261
      IF (i_lay_var_1270 <= k_laytrop_var_1264(iplon_var_1273)) THEN
        DO ig_var_1269 = 1, 6
          IF (i_lay_var_1270 == i_laysolfr_var_1271(iplon_var_1273)) p_sfluxzen_var_1265(iplon_var_1273, ig_var_1269) = sfluxrefc_var_314(ig_var_1269)
          p_taug_var_1266(iplon_var_1273, i_lay_var_1270, ig_var_1269) = 0.0D0
          p_taur_var_1267(iplon_var_1273, i_lay_var_1270, ig_var_1269) = p_colmol_var_1263(iplon_var_1273, i_lay_var_1270) * raylc_var_315(ig_var_1269)
        END DO
      ELSE
        DO ig_var_1269 = 1, 6
          p_taug_var_1266(iplon_var_1273, i_lay_var_1270, ig_var_1269) = 0.0D0
          p_taur_var_1267(iplon_var_1273, i_lay_var_1270, ig_var_1269) = p_colmol_var_1263(iplon_var_1273, i_lay_var_1270) * raylc_var_315(ig_var_1269)
        END DO
      END IF
    END DO
  END DO
  DO ig_var_1269 = 1, 6
    DO i_lay_var_1270 = laytrop_max_var_1275 + 1, i_nlayers_var_1272
      DO iplon_var_1273 = kidia_var_1260, kfdia_var_1261
        p_taug_var_1266(iplon_var_1273, i_lay_var_1270, ig_var_1269) = 0.0D0
        p_taur_var_1267(iplon_var_1273, i_lay_var_1270, ig_var_1269) = p_colmol_var_1263(iplon_var_1273, i_lay_var_1270) * raylc_var_315(ig_var_1269)
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol26
SUBROUTINE srtm_taumol18(kidia_var_1276, kfdia_var_1277, klev_var_1278, p_fac00_var_1279, p_fac01_var_1280, p_fac10_var_1281, p_fac11_var_1282, k_jp_var_1283, k_jt_var_1284, k_jt1_var_1285, p_oneminus_var_1286, p_colh2o_var_1287, p_colch4_var_1288, p_colmol_var_1289, k_laytrop_var_1290, p_selffac_var_1291, p_selffrac_var_1292, k_indself_var_1293, p_forfac_var_1294, p_forfrac_var_1295, k_indfor_var_1296, p_sfluxzen_var_1297, p_taug_var_1298, p_taur_var_1299, prmu0_var_1300)
  USE yoesrta18, ONLY: absa_var_257, absb_var_258, forrefc_var_260, layreffr_var_256, rayl_var_254, selfrefc_var_259, sfluxrefc_var_261, strrat_var_255
  USE yoesrtwn, ONLY: nspa_var_334, nspb_var_335
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1276, kfdia_var_1277
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1278
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1279(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1280(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1281(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1282(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1283(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1284(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1285(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1286(kidia_var_1276 : kfdia_var_1277)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1287(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_1288(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1289(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1290(kidia_var_1276 : kfdia_var_1277)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1291(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1292(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1293(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1294(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1295(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1296(kidia_var_1276 : kfdia_var_1277, klev_var_1278)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1297(kidia_var_1276 : kfdia_var_1277, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1298(kidia_var_1276 : kfdia_var_1277, klev_var_1278, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1299(kidia_var_1276 : kfdia_var_1277, klev_var_1278, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1300(kidia_var_1276 : kfdia_var_1277)
  INTEGER(KIND = 4) :: ig_var_1301, ind0_var_1302, ind1_var_1303, inds_var_1304, indf_var_1305, js_var_1306, i_lay_var_1307, i_laysolfr_var_1308(kidia_var_1276 : kfdia_var_1277), i_nlayers_var_1309, iplon_var_1310
  INTEGER(KIND = 4) :: laytrop_min_var_1311, laytrop_max_var_1312
  REAL(KIND = 8) :: z_fs_var_1313, z_speccomb_var_1314, z_specmult_var_1315, z_specparm_var_1316, z_tauray_var_1317
  laytrop_min_var_1311 = MINVAL(k_laytrop_var_1290(kidia_var_1276 : kfdia_var_1277))
  laytrop_max_var_1312 = MAXVAL(k_laytrop_var_1290(kidia_var_1276 : kfdia_var_1277))
  i_nlayers_var_1309 = klev_var_1278
  DO iplon_var_1310 = kidia_var_1276, kfdia_var_1277
    i_laysolfr_var_1308(iplon_var_1310) = k_laytrop_var_1290(iplon_var_1310)
  END DO
  DO i_lay_var_1307 = 1, laytrop_min_var_1311
    DO iplon_var_1310 = kidia_var_1276, kfdia_var_1277
      IF (k_jp_var_1283(iplon_var_1310, i_lay_var_1307) < layreffr_var_256 .AND. k_jp_var_1283(iplon_var_1310, i_lay_var_1307 + 1) >= layreffr_var_256) i_laysolfr_var_1308(iplon_var_1310) = MIN(i_lay_var_1307 + 1, k_laytrop_var_1290(iplon_var_1310))
      z_speccomb_var_1314 = p_colh2o_var_1287(iplon_var_1310, i_lay_var_1307) + strrat_var_255 * p_colch4_var_1288(iplon_var_1310, i_lay_var_1307)
      z_specparm_var_1316 = p_colh2o_var_1287(iplon_var_1310, i_lay_var_1307) / z_speccomb_var_1314
      z_specparm_var_1316 = MIN(p_oneminus_var_1286(iplon_var_1310), z_specparm_var_1316)
      z_specmult_var_1315 = 8.0D0 * (z_specparm_var_1316)
      js_var_1306 = 1 + INT(z_specmult_var_1315)
      z_fs_var_1313 = z_specmult_var_1315 - AINT(z_specmult_var_1315)
      ind0_var_1302 = ((k_jp_var_1283(iplon_var_1310, i_lay_var_1307) - 1) * 5 + (k_jt_var_1284(iplon_var_1310, i_lay_var_1307) - 1)) * nspa_var_334(18) + js_var_1306
      ind1_var_1303 = (k_jp_var_1283(iplon_var_1310, i_lay_var_1307) * 5 + (k_jt1_var_1285(iplon_var_1310, i_lay_var_1307) - 1)) * nspa_var_334(18) + js_var_1306
      inds_var_1304 = k_indself_var_1293(iplon_var_1310, i_lay_var_1307)
      indf_var_1305 = k_indfor_var_1296(iplon_var_1310, i_lay_var_1307)
      z_tauray_var_1317 = p_colmol_var_1289(iplon_var_1310, i_lay_var_1307) * rayl_var_254
      DO ig_var_1301 = 1, 8
        p_taug_var_1298(iplon_var_1310, i_lay_var_1307, ig_var_1301) = z_speccomb_var_1314 * ((1.0D0 - z_fs_var_1313) * (absa_var_257(ind0_var_1302, ig_var_1301) * p_fac00_var_1279(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind0_var_1302 + 9, ig_var_1301) * p_fac10_var_1281(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind1_var_1303, ig_var_1301) * p_fac01_var_1280(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind1_var_1303 + 9, ig_var_1301) * p_fac11_var_1282(iplon_var_1310, i_lay_var_1307)) + z_fs_var_1313 * (absa_var_257(ind0_var_1302 + 1, ig_var_1301) * p_fac00_var_1279(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind0_var_1302 + 10, ig_var_1301) * p_fac10_var_1281(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind1_var_1303 + 1, ig_var_1301) * p_fac01_var_1280(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind1_var_1303 + 10, ig_var_1301) * p_fac11_var_1282(iplon_var_1310, i_lay_var_1307))) + p_colh2o_var_1287(iplon_var_1310, i_lay_var_1307) * (p_selffac_var_1291(iplon_var_1310, i_lay_var_1307) * (selfrefc_var_259(inds_var_1304, ig_var_1301) + p_selffrac_var_1292(iplon_var_1310, i_lay_var_1307) * (selfrefc_var_259(inds_var_1304 + 1, ig_var_1301) - selfrefc_var_259(inds_var_1304, ig_var_1301))) + p_forfac_var_1294(iplon_var_1310, i_lay_var_1307) * (forrefc_var_260(indf_var_1305, ig_var_1301) + p_forfrac_var_1295(iplon_var_1310, i_lay_var_1307) * (forrefc_var_260(indf_var_1305 + 1, ig_var_1301) - forrefc_var_260(indf_var_1305, ig_var_1301))))
        IF (i_lay_var_1307 == i_laysolfr_var_1308(iplon_var_1310)) p_sfluxzen_var_1297(iplon_var_1310, ig_var_1301) = sfluxrefc_var_261(ig_var_1301, js_var_1306) + z_fs_var_1313 * (sfluxrefc_var_261(ig_var_1301, js_var_1306 + 1) - sfluxrefc_var_261(ig_var_1301, js_var_1306))
        p_taur_var_1299(iplon_var_1310, i_lay_var_1307, ig_var_1301) = z_tauray_var_1317
      END DO
    END DO
  END DO
  DO i_lay_var_1307 = laytrop_min_var_1311 + 1, laytrop_max_var_1312
    DO iplon_var_1310 = kidia_var_1276, kfdia_var_1277
      IF (i_lay_var_1307 <= k_laytrop_var_1290(iplon_var_1310)) THEN
        IF (k_jp_var_1283(iplon_var_1310, i_lay_var_1307) < layreffr_var_256 .AND. k_jp_var_1283(iplon_var_1310, i_lay_var_1307 + 1) >= layreffr_var_256) i_laysolfr_var_1308(iplon_var_1310) = MIN(i_lay_var_1307 + 1, k_laytrop_var_1290(iplon_var_1310))
        z_speccomb_var_1314 = p_colh2o_var_1287(iplon_var_1310, i_lay_var_1307) + strrat_var_255 * p_colch4_var_1288(iplon_var_1310, i_lay_var_1307)
        z_specparm_var_1316 = p_colh2o_var_1287(iplon_var_1310, i_lay_var_1307) / z_speccomb_var_1314
        z_specparm_var_1316 = MIN(p_oneminus_var_1286(iplon_var_1310), z_specparm_var_1316)
        z_specmult_var_1315 = 8.0D0 * (z_specparm_var_1316)
        js_var_1306 = 1 + INT(z_specmult_var_1315)
        z_fs_var_1313 = z_specmult_var_1315 - AINT(z_specmult_var_1315)
        ind0_var_1302 = ((k_jp_var_1283(iplon_var_1310, i_lay_var_1307) - 1) * 5 + (k_jt_var_1284(iplon_var_1310, i_lay_var_1307) - 1)) * nspa_var_334(18) + js_var_1306
        ind1_var_1303 = (k_jp_var_1283(iplon_var_1310, i_lay_var_1307) * 5 + (k_jt1_var_1285(iplon_var_1310, i_lay_var_1307) - 1)) * nspa_var_334(18) + js_var_1306
        inds_var_1304 = k_indself_var_1293(iplon_var_1310, i_lay_var_1307)
        indf_var_1305 = k_indfor_var_1296(iplon_var_1310, i_lay_var_1307)
        z_tauray_var_1317 = p_colmol_var_1289(iplon_var_1310, i_lay_var_1307) * rayl_var_254
        DO ig_var_1301 = 1, 8
          p_taug_var_1298(iplon_var_1310, i_lay_var_1307, ig_var_1301) = z_speccomb_var_1314 * ((1.0D0 - z_fs_var_1313) * (absa_var_257(ind0_var_1302, ig_var_1301) * p_fac00_var_1279(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind0_var_1302 + 9, ig_var_1301) * p_fac10_var_1281(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind1_var_1303, ig_var_1301) * p_fac01_var_1280(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind1_var_1303 + 9, ig_var_1301) * p_fac11_var_1282(iplon_var_1310, i_lay_var_1307)) + z_fs_var_1313 * (absa_var_257(ind0_var_1302 + 1, ig_var_1301) * p_fac00_var_1279(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind0_var_1302 + 10, ig_var_1301) * p_fac10_var_1281(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind1_var_1303 + 1, ig_var_1301) * p_fac01_var_1280(iplon_var_1310, i_lay_var_1307) + absa_var_257(ind1_var_1303 + 10, ig_var_1301) * p_fac11_var_1282(iplon_var_1310, i_lay_var_1307))) + p_colh2o_var_1287(iplon_var_1310, i_lay_var_1307) * (p_selffac_var_1291(iplon_var_1310, i_lay_var_1307) * (selfrefc_var_259(inds_var_1304, ig_var_1301) + p_selffrac_var_1292(iplon_var_1310, i_lay_var_1307) * (selfrefc_var_259(inds_var_1304 + 1, ig_var_1301) - selfrefc_var_259(inds_var_1304, ig_var_1301))) + p_forfac_var_1294(iplon_var_1310, i_lay_var_1307) * (forrefc_var_260(indf_var_1305, ig_var_1301) + p_forfrac_var_1295(iplon_var_1310, i_lay_var_1307) * (forrefc_var_260(indf_var_1305 + 1, ig_var_1301) - forrefc_var_260(indf_var_1305, ig_var_1301))))
          IF (i_lay_var_1307 == i_laysolfr_var_1308(iplon_var_1310)) p_sfluxzen_var_1297(iplon_var_1310, ig_var_1301) = sfluxrefc_var_261(ig_var_1301, js_var_1306) + z_fs_var_1313 * (sfluxrefc_var_261(ig_var_1301, js_var_1306 + 1) - sfluxrefc_var_261(ig_var_1301, js_var_1306))
          p_taur_var_1299(iplon_var_1310, i_lay_var_1307, ig_var_1301) = z_tauray_var_1317
        END DO
      ELSE
        ind0_var_1302 = ((k_jp_var_1283(iplon_var_1310, i_lay_var_1307) - 13) * 5 + (k_jt_var_1284(iplon_var_1310, i_lay_var_1307) - 1)) * nspb_var_335(18) + 1
        ind1_var_1303 = ((k_jp_var_1283(iplon_var_1310, i_lay_var_1307) - 12) * 5 + (k_jt1_var_1285(iplon_var_1310, i_lay_var_1307) - 1)) * nspb_var_335(18) + 1
        z_tauray_var_1317 = p_colmol_var_1289(iplon_var_1310, i_lay_var_1307) * rayl_var_254
        DO ig_var_1301 = 1, 8
          p_taug_var_1298(iplon_var_1310, i_lay_var_1307, ig_var_1301) = p_colch4_var_1288(iplon_var_1310, i_lay_var_1307) * (p_fac00_var_1279(iplon_var_1310, i_lay_var_1307) * absb_var_258(ind0_var_1302, ig_var_1301) + p_fac10_var_1281(iplon_var_1310, i_lay_var_1307) * absb_var_258(ind0_var_1302 + 1, ig_var_1301) + p_fac01_var_1280(iplon_var_1310, i_lay_var_1307) * absb_var_258(ind1_var_1303, ig_var_1301) + p_fac11_var_1282(iplon_var_1310, i_lay_var_1307) * absb_var_258(ind1_var_1303 + 1, ig_var_1301))
          p_taur_var_1299(iplon_var_1310, i_lay_var_1307, ig_var_1301) = z_tauray_var_1317
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1307 = laytrop_max_var_1312 + 1, i_nlayers_var_1309
    DO iplon_var_1310 = kidia_var_1276, kfdia_var_1277
      ind0_var_1302 = ((k_jp_var_1283(iplon_var_1310, i_lay_var_1307) - 13) * 5 + (k_jt_var_1284(iplon_var_1310, i_lay_var_1307) - 1)) * nspb_var_335(18) + 1
      ind1_var_1303 = ((k_jp_var_1283(iplon_var_1310, i_lay_var_1307) - 12) * 5 + (k_jt1_var_1285(iplon_var_1310, i_lay_var_1307) - 1)) * nspb_var_335(18) + 1
      z_tauray_var_1317 = p_colmol_var_1289(iplon_var_1310, i_lay_var_1307) * rayl_var_254
      DO ig_var_1301 = 1, 8
        p_taug_var_1298(iplon_var_1310, i_lay_var_1307, ig_var_1301) = p_colch4_var_1288(iplon_var_1310, i_lay_var_1307) * (p_fac00_var_1279(iplon_var_1310, i_lay_var_1307) * absb_var_258(ind0_var_1302, ig_var_1301) + p_fac10_var_1281(iplon_var_1310, i_lay_var_1307) * absb_var_258(ind0_var_1302 + 1, ig_var_1301) + p_fac01_var_1280(iplon_var_1310, i_lay_var_1307) * absb_var_258(ind1_var_1303, ig_var_1301) + p_fac11_var_1282(iplon_var_1310, i_lay_var_1307) * absb_var_258(ind1_var_1303 + 1, ig_var_1301))
        p_taur_var_1299(iplon_var_1310, i_lay_var_1307, ig_var_1301) = z_tauray_var_1317
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol18
SUBROUTINE srtm_taumol24(kidia_var_1318, kfdia_var_1319, klev_var_1320, p_fac00_var_1321, p_fac01_var_1322, p_fac10_var_1323, p_fac11_var_1324, k_jp_var_1325, k_jt_var_1326, k_jt1_var_1327, p_oneminus_var_1328, p_colh2o_var_1329, p_colmol_var_1330, p_colo2_var_1331, p_colo3_var_1332, k_laytrop_var_1333, p_selffac_var_1334, p_selffrac_var_1335, k_indself_var_1336, p_forfac_var_1337, p_forfrac_var_1338, k_indfor_var_1339, p_sfluxzen_var_1340, p_taug_var_1341, p_taur_var_1342, prmu0_var_1343)
  USE yoesrta24, ONLY: absa_var_301, absb_var_302, abso3ac_var_306, abso3bc_var_307, forrefc_var_304, layreffr_var_300, raylac, raylbc, selfrefc_var_303, sfluxrefc_var_305, strrat_var_299
  USE yoesrtwn, ONLY: nspa_var_334, nspb_var_335
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1318, kfdia_var_1319
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1320
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1321(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1322(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1323(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1324(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1325(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1326(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1327(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1328(kidia_var_1318 : kfdia_var_1319)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1329(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1330(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_1331(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1332(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1333(kidia_var_1318 : kfdia_var_1319)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1334(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1335(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1336(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1337(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1338(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1339(kidia_var_1318 : kfdia_var_1319, klev_var_1320)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1340(kidia_var_1318 : kfdia_var_1319, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1341(kidia_var_1318 : kfdia_var_1319, klev_var_1320, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1342(kidia_var_1318 : kfdia_var_1319, klev_var_1320, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1343(kidia_var_1318 : kfdia_var_1319)
  INTEGER(KIND = 4) :: ig_var_1344, ind0_var_1345, ind1_var_1346, inds_var_1347, indf_var_1348, js_var_1349, i_lay_var_1350, i_laysolfr_var_1351(kidia_var_1318 : kfdia_var_1319), i_nlayers_var_1352, iplon_var_1353
  INTEGER(KIND = 4) :: laytrop_min_var_1354, laytrop_max_var_1355
  REAL(KIND = 8) :: z_fs_var_1356, z_speccomb_var_1357, z_specmult_var_1358, z_specparm_var_1359, z_tauray_var_1360
  laytrop_min_var_1354 = MINVAL(k_laytrop_var_1333(kidia_var_1318 : kfdia_var_1319))
  laytrop_max_var_1355 = MAXVAL(k_laytrop_var_1333(kidia_var_1318 : kfdia_var_1319))
  i_nlayers_var_1352 = klev_var_1320
  DO iplon_var_1353 = kidia_var_1318, kfdia_var_1319
    i_laysolfr_var_1351(iplon_var_1353) = k_laytrop_var_1333(iplon_var_1353)
  END DO
  DO i_lay_var_1350 = 1, laytrop_min_var_1354
    DO iplon_var_1353 = kidia_var_1318, kfdia_var_1319
      IF (k_jp_var_1325(iplon_var_1353, i_lay_var_1350) < layreffr_var_300 .AND. k_jp_var_1325(iplon_var_1353, i_lay_var_1350 + 1) >= layreffr_var_300) i_laysolfr_var_1351(iplon_var_1353) = MIN(i_lay_var_1350 + 1, k_laytrop_var_1333(iplon_var_1353))
      z_speccomb_var_1357 = p_colh2o_var_1329(iplon_var_1353, i_lay_var_1350) + strrat_var_299 * p_colo2_var_1331(iplon_var_1353, i_lay_var_1350)
      z_specparm_var_1359 = p_colh2o_var_1329(iplon_var_1353, i_lay_var_1350) / z_speccomb_var_1357
      z_specparm_var_1359 = MIN(p_oneminus_var_1328(iplon_var_1353), z_specparm_var_1359)
      z_specmult_var_1358 = 8.0D0 * (z_specparm_var_1359)
      js_var_1349 = 1 + INT(z_specmult_var_1358)
      z_fs_var_1356 = z_specmult_var_1358 - AINT(z_specmult_var_1358)
      ind0_var_1345 = ((k_jp_var_1325(iplon_var_1353, i_lay_var_1350) - 1) * 5 + (k_jt_var_1326(iplon_var_1353, i_lay_var_1350) - 1)) * nspa_var_334(24) + js_var_1349
      ind1_var_1346 = (k_jp_var_1325(iplon_var_1353, i_lay_var_1350) * 5 + (k_jt1_var_1327(iplon_var_1353, i_lay_var_1350) - 1)) * nspa_var_334(24) + js_var_1349
      inds_var_1347 = k_indself_var_1336(iplon_var_1353, i_lay_var_1350)
      indf_var_1348 = k_indfor_var_1339(iplon_var_1353, i_lay_var_1350)
      DO ig_var_1344 = 1, 8
        z_tauray_var_1360 = p_colmol_var_1330(iplon_var_1353, i_lay_var_1350) * (raylac(ig_var_1344, js_var_1349) + z_fs_var_1356 * (raylac(ig_var_1344, js_var_1349 + 1) - raylac(ig_var_1344, js_var_1349)))
        p_taug_var_1341(iplon_var_1353, i_lay_var_1350, ig_var_1344) = z_speccomb_var_1357 * ((1.0D0 - z_fs_var_1356) * (absa_var_301(ind0_var_1345, ig_var_1344) * p_fac00_var_1321(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind0_var_1345 + 9, ig_var_1344) * p_fac10_var_1323(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind1_var_1346, ig_var_1344) * p_fac01_var_1322(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind1_var_1346 + 9, ig_var_1344) * p_fac11_var_1324(iplon_var_1353, i_lay_var_1350)) + z_fs_var_1356 * (absa_var_301(ind0_var_1345 + 1, ig_var_1344) * p_fac00_var_1321(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind0_var_1345 + 10, ig_var_1344) * p_fac10_var_1323(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind1_var_1346 + 1, ig_var_1344) * p_fac01_var_1322(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind1_var_1346 + 10, ig_var_1344) * p_fac11_var_1324(iplon_var_1353, i_lay_var_1350))) + p_colo3_var_1332(iplon_var_1353, i_lay_var_1350) * abso3ac_var_306(ig_var_1344) + p_colh2o_var_1329(iplon_var_1353, i_lay_var_1350) * (p_selffac_var_1334(iplon_var_1353, i_lay_var_1350) * (selfrefc_var_303(inds_var_1347, ig_var_1344) + p_selffrac_var_1335(iplon_var_1353, i_lay_var_1350) * (selfrefc_var_303(inds_var_1347 + 1, ig_var_1344) - selfrefc_var_303(inds_var_1347, ig_var_1344))) + p_forfac_var_1337(iplon_var_1353, i_lay_var_1350) * (forrefc_var_304(indf_var_1348, ig_var_1344) + p_forfrac_var_1338(iplon_var_1353, i_lay_var_1350) * (forrefc_var_304(indf_var_1348 + 1, ig_var_1344) - forrefc_var_304(indf_var_1348, ig_var_1344))))
        IF (i_lay_var_1350 == i_laysolfr_var_1351(iplon_var_1353)) p_sfluxzen_var_1340(iplon_var_1353, ig_var_1344) = sfluxrefc_var_305(ig_var_1344, js_var_1349) + z_fs_var_1356 * (sfluxrefc_var_305(ig_var_1344, js_var_1349 + 1) - sfluxrefc_var_305(ig_var_1344, js_var_1349))
        p_taur_var_1342(iplon_var_1353, i_lay_var_1350, ig_var_1344) = z_tauray_var_1360
      END DO
    END DO
  END DO
  DO i_lay_var_1350 = laytrop_min_var_1354 + 1, laytrop_max_var_1355
    DO iplon_var_1353 = kidia_var_1318, kfdia_var_1319
      IF (i_lay_var_1350 <= k_laytrop_var_1333(iplon_var_1353)) THEN
        IF (k_jp_var_1325(iplon_var_1353, i_lay_var_1350) < layreffr_var_300 .AND. k_jp_var_1325(iplon_var_1353, i_lay_var_1350 + 1) >= layreffr_var_300) i_laysolfr_var_1351(iplon_var_1353) = MIN(i_lay_var_1350 + 1, k_laytrop_var_1333(iplon_var_1353))
        z_speccomb_var_1357 = p_colh2o_var_1329(iplon_var_1353, i_lay_var_1350) + strrat_var_299 * p_colo2_var_1331(iplon_var_1353, i_lay_var_1350)
        z_specparm_var_1359 = p_colh2o_var_1329(iplon_var_1353, i_lay_var_1350) / z_speccomb_var_1357
        z_specparm_var_1359 = MIN(p_oneminus_var_1328(iplon_var_1353), z_specparm_var_1359)
        z_specmult_var_1358 = 8.0D0 * (z_specparm_var_1359)
        js_var_1349 = 1 + INT(z_specmult_var_1358)
        z_fs_var_1356 = z_specmult_var_1358 - AINT(z_specmult_var_1358)
        ind0_var_1345 = ((k_jp_var_1325(iplon_var_1353, i_lay_var_1350) - 1) * 5 + (k_jt_var_1326(iplon_var_1353, i_lay_var_1350) - 1)) * nspa_var_334(24) + js_var_1349
        ind1_var_1346 = (k_jp_var_1325(iplon_var_1353, i_lay_var_1350) * 5 + (k_jt1_var_1327(iplon_var_1353, i_lay_var_1350) - 1)) * nspa_var_334(24) + js_var_1349
        inds_var_1347 = k_indself_var_1336(iplon_var_1353, i_lay_var_1350)
        indf_var_1348 = k_indfor_var_1339(iplon_var_1353, i_lay_var_1350)
        DO ig_var_1344 = 1, 8
          z_tauray_var_1360 = p_colmol_var_1330(iplon_var_1353, i_lay_var_1350) * (raylac(ig_var_1344, js_var_1349) + z_fs_var_1356 * (raylac(ig_var_1344, js_var_1349 + 1) - raylac(ig_var_1344, js_var_1349)))
          p_taug_var_1341(iplon_var_1353, i_lay_var_1350, ig_var_1344) = z_speccomb_var_1357 * ((1.0D0 - z_fs_var_1356) * (absa_var_301(ind0_var_1345, ig_var_1344) * p_fac00_var_1321(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind0_var_1345 + 9, ig_var_1344) * p_fac10_var_1323(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind1_var_1346, ig_var_1344) * p_fac01_var_1322(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind1_var_1346 + 9, ig_var_1344) * p_fac11_var_1324(iplon_var_1353, i_lay_var_1350)) + z_fs_var_1356 * (absa_var_301(ind0_var_1345 + 1, ig_var_1344) * p_fac00_var_1321(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind0_var_1345 + 10, ig_var_1344) * p_fac10_var_1323(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind1_var_1346 + 1, ig_var_1344) * p_fac01_var_1322(iplon_var_1353, i_lay_var_1350) + absa_var_301(ind1_var_1346 + 10, ig_var_1344) * p_fac11_var_1324(iplon_var_1353, i_lay_var_1350))) + p_colo3_var_1332(iplon_var_1353, i_lay_var_1350) * abso3ac_var_306(ig_var_1344) + p_colh2o_var_1329(iplon_var_1353, i_lay_var_1350) * (p_selffac_var_1334(iplon_var_1353, i_lay_var_1350) * (selfrefc_var_303(inds_var_1347, ig_var_1344) + p_selffrac_var_1335(iplon_var_1353, i_lay_var_1350) * (selfrefc_var_303(inds_var_1347 + 1, ig_var_1344) - selfrefc_var_303(inds_var_1347, ig_var_1344))) + p_forfac_var_1337(iplon_var_1353, i_lay_var_1350) * (forrefc_var_304(indf_var_1348, ig_var_1344) + p_forfrac_var_1338(iplon_var_1353, i_lay_var_1350) * (forrefc_var_304(indf_var_1348 + 1, ig_var_1344) - forrefc_var_304(indf_var_1348, ig_var_1344))))
          IF (i_lay_var_1350 == i_laysolfr_var_1351(iplon_var_1353)) p_sfluxzen_var_1340(iplon_var_1353, ig_var_1344) = sfluxrefc_var_305(ig_var_1344, js_var_1349) + z_fs_var_1356 * (sfluxrefc_var_305(ig_var_1344, js_var_1349 + 1) - sfluxrefc_var_305(ig_var_1344, js_var_1349))
          p_taur_var_1342(iplon_var_1353, i_lay_var_1350, ig_var_1344) = z_tauray_var_1360
        END DO
      ELSE
        ind0_var_1345 = ((k_jp_var_1325(iplon_var_1353, i_lay_var_1350) - 13) * 5 + (k_jt_var_1326(iplon_var_1353, i_lay_var_1350) - 1)) * nspb_var_335(24) + 1
        ind1_var_1346 = ((k_jp_var_1325(iplon_var_1353, i_lay_var_1350) - 12) * 5 + (k_jt1_var_1327(iplon_var_1353, i_lay_var_1350) - 1)) * nspb_var_335(24) + 1
        DO ig_var_1344 = 1, 8
          z_tauray_var_1360 = p_colmol_var_1330(iplon_var_1353, i_lay_var_1350) * raylbc(ig_var_1344)
          p_taug_var_1341(iplon_var_1353, i_lay_var_1350, ig_var_1344) = p_colo2_var_1331(iplon_var_1353, i_lay_var_1350) * (p_fac00_var_1321(iplon_var_1353, i_lay_var_1350) * absb_var_302(ind0_var_1345, ig_var_1344) + p_fac10_var_1323(iplon_var_1353, i_lay_var_1350) * absb_var_302(ind0_var_1345 + 1, ig_var_1344) + p_fac01_var_1322(iplon_var_1353, i_lay_var_1350) * absb_var_302(ind1_var_1346, ig_var_1344) + p_fac11_var_1324(iplon_var_1353, i_lay_var_1350) * absb_var_302(ind1_var_1346 + 1, ig_var_1344)) + p_colo3_var_1332(iplon_var_1353, i_lay_var_1350) * abso3bc_var_307(ig_var_1344)
          p_taur_var_1342(iplon_var_1353, i_lay_var_1350, ig_var_1344) = z_tauray_var_1360
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1350 = laytrop_max_var_1355 + 1, i_nlayers_var_1352
    DO iplon_var_1353 = kidia_var_1318, kfdia_var_1319
      ind0_var_1345 = ((k_jp_var_1325(iplon_var_1353, i_lay_var_1350) - 13) * 5 + (k_jt_var_1326(iplon_var_1353, i_lay_var_1350) - 1)) * nspb_var_335(24) + 1
      ind1_var_1346 = ((k_jp_var_1325(iplon_var_1353, i_lay_var_1350) - 12) * 5 + (k_jt1_var_1327(iplon_var_1353, i_lay_var_1350) - 1)) * nspb_var_335(24) + 1
      DO ig_var_1344 = 1, 8
        z_tauray_var_1360 = p_colmol_var_1330(iplon_var_1353, i_lay_var_1350) * raylbc(ig_var_1344)
        p_taug_var_1341(iplon_var_1353, i_lay_var_1350, ig_var_1344) = p_colo2_var_1331(iplon_var_1353, i_lay_var_1350) * (p_fac00_var_1321(iplon_var_1353, i_lay_var_1350) * absb_var_302(ind0_var_1345, ig_var_1344) + p_fac10_var_1323(iplon_var_1353, i_lay_var_1350) * absb_var_302(ind0_var_1345 + 1, ig_var_1344) + p_fac01_var_1322(iplon_var_1353, i_lay_var_1350) * absb_var_302(ind1_var_1346, ig_var_1344) + p_fac11_var_1324(iplon_var_1353, i_lay_var_1350) * absb_var_302(ind1_var_1346 + 1, ig_var_1344)) + p_colo3_var_1332(iplon_var_1353, i_lay_var_1350) * abso3bc_var_307(ig_var_1344)
        p_taur_var_1342(iplon_var_1353, i_lay_var_1350, ig_var_1344) = z_tauray_var_1360
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol24
SUBROUTINE srtm_taumol25(kidia_var_1361, kfdia_var_1362, klev_var_1363, p_fac00_var_1364, p_fac01_var_1365, p_fac10_var_1366, p_fac11_var_1367, k_jp_var_1368, k_jt_var_1369, k_jt1_var_1370, p_colh2o_var_1371, p_colmol_var_1372, p_colo3_var_1373, k_laytrop_var_1374, p_sfluxzen_var_1375, p_taug_var_1376, p_taur_var_1377, prmu0_var_1378)
  USE yoesrta25, ONLY: absa_var_309, abso3ac_var_312, abso3bc_var_313, layreffr_var_308, raylc_var_311, sfluxrefc_var_310
  USE yoesrtwn, ONLY: nspa_var_334
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1361, kfdia_var_1362
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1363
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1364(kidia_var_1361 : kfdia_var_1362, klev_var_1363)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1365(kidia_var_1361 : kfdia_var_1362, klev_var_1363)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1366(kidia_var_1361 : kfdia_var_1362, klev_var_1363)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1367(kidia_var_1361 : kfdia_var_1362, klev_var_1363)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1368(kidia_var_1361 : kfdia_var_1362, klev_var_1363)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1369(kidia_var_1361 : kfdia_var_1362, klev_var_1363)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1370(kidia_var_1361 : kfdia_var_1362, klev_var_1363)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1371(kidia_var_1361 : kfdia_var_1362, klev_var_1363)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1372(kidia_var_1361 : kfdia_var_1362, klev_var_1363)
  REAL(KIND = 8), INTENT(IN) :: p_colo3_var_1373(kidia_var_1361 : kfdia_var_1362, klev_var_1363)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1374(kidia_var_1361 : kfdia_var_1362)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1375(kidia_var_1361 : kfdia_var_1362, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1376(kidia_var_1361 : kfdia_var_1362, klev_var_1363, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1377(kidia_var_1361 : kfdia_var_1362, klev_var_1363, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1378(kidia_var_1361 : kfdia_var_1362)
  INTEGER(KIND = 4) :: ig_var_1379, ind0_var_1380, ind1_var_1381, i_lay_var_1382, i_laysolfr_var_1383(kidia_var_1361 : kfdia_var_1362), i_nlayers_var_1384, iplon_var_1385
  INTEGER(KIND = 4) :: laytrop_min_var_1386, laytrop_max_var_1387
  REAL(KIND = 8) :: z_tauray_var_1388
  laytrop_min_var_1386 = MINVAL(k_laytrop_var_1374(kidia_var_1361 : kfdia_var_1362))
  laytrop_max_var_1387 = MAXVAL(k_laytrop_var_1374(kidia_var_1361 : kfdia_var_1362))
  i_nlayers_var_1384 = klev_var_1363
  DO iplon_var_1385 = kidia_var_1361, kfdia_var_1362
    i_laysolfr_var_1383(iplon_var_1385) = k_laytrop_var_1374(iplon_var_1385)
  END DO
  DO i_lay_var_1382 = 1, laytrop_min_var_1386
    DO iplon_var_1385 = kidia_var_1361, kfdia_var_1362
      IF (k_jp_var_1368(iplon_var_1385, i_lay_var_1382) < layreffr_var_308 .AND. k_jp_var_1368(iplon_var_1385, i_lay_var_1382 + 1) >= layreffr_var_308) i_laysolfr_var_1383(iplon_var_1385) = MIN(i_lay_var_1382 + 1, k_laytrop_var_1374(iplon_var_1385))
      ind0_var_1380 = ((k_jp_var_1368(iplon_var_1385, i_lay_var_1382) - 1) * 5 + (k_jt_var_1369(iplon_var_1385, i_lay_var_1382) - 1)) * nspa_var_334(25) + 1
      ind1_var_1381 = (k_jp_var_1368(iplon_var_1385, i_lay_var_1382) * 5 + (k_jt1_var_1370(iplon_var_1385, i_lay_var_1382) - 1)) * nspa_var_334(25) + 1
      DO ig_var_1379 = 1, 6
        z_tauray_var_1388 = p_colmol_var_1372(iplon_var_1385, i_lay_var_1382) * raylc_var_311(ig_var_1379)
        p_taug_var_1376(iplon_var_1385, i_lay_var_1382, ig_var_1379) = p_colh2o_var_1371(iplon_var_1385, i_lay_var_1382) * (p_fac00_var_1364(iplon_var_1385, i_lay_var_1382) * absa_var_309(ind0_var_1380, ig_var_1379) + p_fac10_var_1366(iplon_var_1385, i_lay_var_1382) * absa_var_309(ind0_var_1380 + 1, ig_var_1379) + p_fac01_var_1365(iplon_var_1385, i_lay_var_1382) * absa_var_309(ind1_var_1381, ig_var_1379) + p_fac11_var_1367(iplon_var_1385, i_lay_var_1382) * absa_var_309(ind1_var_1381 + 1, ig_var_1379)) + p_colo3_var_1373(iplon_var_1385, i_lay_var_1382) * abso3ac_var_312(ig_var_1379)
        IF (i_lay_var_1382 == i_laysolfr_var_1383(iplon_var_1385)) p_sfluxzen_var_1375(iplon_var_1385, ig_var_1379) = sfluxrefc_var_310(ig_var_1379)
        p_taur_var_1377(iplon_var_1385, i_lay_var_1382, ig_var_1379) = z_tauray_var_1388
      END DO
    END DO
  END DO
  DO i_lay_var_1382 = laytrop_min_var_1386 + 1, laytrop_max_var_1387
    DO iplon_var_1385 = kidia_var_1361, kfdia_var_1362
      IF (i_lay_var_1382 <= k_laytrop_var_1374(iplon_var_1385)) THEN
        IF (k_jp_var_1368(iplon_var_1385, i_lay_var_1382) < layreffr_var_308 .AND. k_jp_var_1368(iplon_var_1385, i_lay_var_1382 + 1) >= layreffr_var_308) i_laysolfr_var_1383(iplon_var_1385) = MIN(i_lay_var_1382 + 1, k_laytrop_var_1374(iplon_var_1385))
        ind0_var_1380 = ((k_jp_var_1368(iplon_var_1385, i_lay_var_1382) - 1) * 5 + (k_jt_var_1369(iplon_var_1385, i_lay_var_1382) - 1)) * nspa_var_334(25) + 1
        ind1_var_1381 = (k_jp_var_1368(iplon_var_1385, i_lay_var_1382) * 5 + (k_jt1_var_1370(iplon_var_1385, i_lay_var_1382) - 1)) * nspa_var_334(25) + 1
        DO ig_var_1379 = 1, 6
          z_tauray_var_1388 = p_colmol_var_1372(iplon_var_1385, i_lay_var_1382) * raylc_var_311(ig_var_1379)
          p_taug_var_1376(iplon_var_1385, i_lay_var_1382, ig_var_1379) = p_colh2o_var_1371(iplon_var_1385, i_lay_var_1382) * (p_fac00_var_1364(iplon_var_1385, i_lay_var_1382) * absa_var_309(ind0_var_1380, ig_var_1379) + p_fac10_var_1366(iplon_var_1385, i_lay_var_1382) * absa_var_309(ind0_var_1380 + 1, ig_var_1379) + p_fac01_var_1365(iplon_var_1385, i_lay_var_1382) * absa_var_309(ind1_var_1381, ig_var_1379) + p_fac11_var_1367(iplon_var_1385, i_lay_var_1382) * absa_var_309(ind1_var_1381 + 1, ig_var_1379)) + p_colo3_var_1373(iplon_var_1385, i_lay_var_1382) * abso3ac_var_312(ig_var_1379)
          IF (i_lay_var_1382 == i_laysolfr_var_1383(iplon_var_1385)) p_sfluxzen_var_1375(iplon_var_1385, ig_var_1379) = sfluxrefc_var_310(ig_var_1379)
          p_taur_var_1377(iplon_var_1385, i_lay_var_1382, ig_var_1379) = z_tauray_var_1388
        END DO
      ELSE
        DO ig_var_1379 = 1, 6
          z_tauray_var_1388 = p_colmol_var_1372(iplon_var_1385, i_lay_var_1382) * raylc_var_311(ig_var_1379)
          p_taug_var_1376(iplon_var_1385, i_lay_var_1382, ig_var_1379) = p_colo3_var_1373(iplon_var_1385, i_lay_var_1382) * abso3bc_var_313(ig_var_1379)
          p_taur_var_1377(iplon_var_1385, i_lay_var_1382, ig_var_1379) = z_tauray_var_1388
        END DO
      END IF
    END DO
  END DO
  DO ig_var_1379 = 1, 6
    DO i_lay_var_1382 = laytrop_max_var_1387 + 1, i_nlayers_var_1384
      DO iplon_var_1385 = kidia_var_1361, kfdia_var_1362
        z_tauray_var_1388 = p_colmol_var_1372(iplon_var_1385, i_lay_var_1382) * raylc_var_311(ig_var_1379)
        p_taug_var_1376(iplon_var_1385, i_lay_var_1382, ig_var_1379) = p_colo3_var_1373(iplon_var_1385, i_lay_var_1382) * abso3bc_var_313(ig_var_1379)
        p_taur_var_1377(iplon_var_1385, i_lay_var_1382, ig_var_1379) = z_tauray_var_1388
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol25
SUBROUTINE srtm_taumol19(kidia_var_1389, kfdia_var_1390, klev_var_1391, p_fac00_var_1392, p_fac01_var_1393, p_fac10_var_1394, p_fac11_var_1395, k_jp_var_1396, k_jt_var_1397, k_jt1_var_1398, p_oneminus_var_1399, p_colh2o_var_1400, p_colco2_var_1401, p_colmol_var_1402, k_laytrop_var_1403, p_selffac_var_1404, p_selffrac_var_1405, k_indself_var_1406, p_forfac_var_1407, p_forfrac_var_1408, k_indfor_var_1409, p_sfluxzen_var_1410, p_taug_var_1411, p_taur_var_1412, prmu0_var_1413)
  USE yoesrta19, ONLY: absa_var_265, absb_var_266, forrefc_var_268, layreffr_var_264, rayl_var_262, selfrefc_var_267, sfluxrefc_var_269, strrat_var_263
  USE yoesrtwn, ONLY: nspa_var_334, nspb_var_335
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1389, kfdia_var_1390
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1391
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1392(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1393(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1394(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1395(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1396(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1397(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1398(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1399(kidia_var_1389 : kfdia_var_1390)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1400(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1401(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1402(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1403(kidia_var_1389 : kfdia_var_1390)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1404(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1405(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1406(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1407(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1408(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1409(kidia_var_1389 : kfdia_var_1390, klev_var_1391)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1410(kidia_var_1389 : kfdia_var_1390, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1411(kidia_var_1389 : kfdia_var_1390, klev_var_1391, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1412(kidia_var_1389 : kfdia_var_1390, klev_var_1391, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1413(kidia_var_1389 : kfdia_var_1390)
  INTEGER(KIND = 4) :: ig_var_1414, ind0_var_1415, ind1_var_1416, inds_var_1417, indf_var_1418, js_var_1419, i_lay_var_1420, i_laysolfr_var_1421(kidia_var_1389 : kfdia_var_1390), i_nlayers_var_1422, iplon_var_1423
  INTEGER(KIND = 4) :: laytrop_min_var_1424, laytrop_max_var_1425
  REAL(KIND = 8) :: z_fs_var_1426, z_speccomb_var_1427, z_specmult_var_1428, z_specparm_var_1429, z_tauray_var_1430
  laytrop_min_var_1424 = MINVAL(k_laytrop_var_1403(kidia_var_1389 : kfdia_var_1390))
  laytrop_max_var_1425 = MAXVAL(k_laytrop_var_1403(kidia_var_1389 : kfdia_var_1390))
  i_nlayers_var_1422 = klev_var_1391
  DO iplon_var_1423 = kidia_var_1389, kfdia_var_1390
    i_laysolfr_var_1421(iplon_var_1423) = k_laytrop_var_1403(iplon_var_1423)
  END DO
  DO i_lay_var_1420 = 1, laytrop_min_var_1424
    DO iplon_var_1423 = kidia_var_1389, kfdia_var_1390
      IF (k_jp_var_1396(iplon_var_1423, i_lay_var_1420) < layreffr_var_264 .AND. k_jp_var_1396(iplon_var_1423, i_lay_var_1420 + 1) >= layreffr_var_264) i_laysolfr_var_1421(iplon_var_1423) = MIN(i_lay_var_1420 + 1, k_laytrop_var_1403(iplon_var_1423))
      z_speccomb_var_1427 = p_colh2o_var_1400(iplon_var_1423, i_lay_var_1420) + strrat_var_263 * p_colco2_var_1401(iplon_var_1423, i_lay_var_1420)
      z_specparm_var_1429 = p_colh2o_var_1400(iplon_var_1423, i_lay_var_1420) / z_speccomb_var_1427
      z_specparm_var_1429 = MIN(p_oneminus_var_1399(iplon_var_1423), z_specparm_var_1429)
      z_specmult_var_1428 = 8.0D0 * (z_specparm_var_1429)
      js_var_1419 = 1 + INT(z_specmult_var_1428)
      z_fs_var_1426 = z_specmult_var_1428 - AINT(z_specmult_var_1428)
      ind0_var_1415 = ((k_jp_var_1396(iplon_var_1423, i_lay_var_1420) - 1) * 5 + (k_jt_var_1397(iplon_var_1423, i_lay_var_1420) - 1)) * nspa_var_334(19) + js_var_1419
      ind1_var_1416 = (k_jp_var_1396(iplon_var_1423, i_lay_var_1420) * 5 + (k_jt1_var_1398(iplon_var_1423, i_lay_var_1420) - 1)) * nspa_var_334(19) + js_var_1419
      inds_var_1417 = k_indself_var_1406(iplon_var_1423, i_lay_var_1420)
      indf_var_1418 = k_indfor_var_1409(iplon_var_1423, i_lay_var_1420)
      z_tauray_var_1430 = p_colmol_var_1402(iplon_var_1423, i_lay_var_1420) * rayl_var_262
      DO ig_var_1414 = 1, 8
        p_taug_var_1411(iplon_var_1423, i_lay_var_1420, ig_var_1414) = z_speccomb_var_1427 * ((1.0D0 - z_fs_var_1426) * (absa_var_265(ind0_var_1415, ig_var_1414) * p_fac00_var_1392(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind0_var_1415 + 9, ig_var_1414) * p_fac10_var_1394(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind1_var_1416, ig_var_1414) * p_fac01_var_1393(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind1_var_1416 + 9, ig_var_1414) * p_fac11_var_1395(iplon_var_1423, i_lay_var_1420)) + z_fs_var_1426 * (absa_var_265(ind0_var_1415 + 1, ig_var_1414) * p_fac00_var_1392(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind0_var_1415 + 10, ig_var_1414) * p_fac10_var_1394(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind1_var_1416 + 1, ig_var_1414) * p_fac01_var_1393(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind1_var_1416 + 10, ig_var_1414) * p_fac11_var_1395(iplon_var_1423, i_lay_var_1420))) + p_colh2o_var_1400(iplon_var_1423, i_lay_var_1420) * (p_selffac_var_1404(iplon_var_1423, i_lay_var_1420) * (selfrefc_var_267(inds_var_1417, ig_var_1414) + p_selffrac_var_1405(iplon_var_1423, i_lay_var_1420) * (selfrefc_var_267(inds_var_1417 + 1, ig_var_1414) - selfrefc_var_267(inds_var_1417, ig_var_1414))) + p_forfac_var_1407(iplon_var_1423, i_lay_var_1420) * (forrefc_var_268(indf_var_1418, ig_var_1414) + p_forfrac_var_1408(iplon_var_1423, i_lay_var_1420) * (forrefc_var_268(indf_var_1418 + 1, ig_var_1414) - forrefc_var_268(indf_var_1418, ig_var_1414))))
        IF (i_lay_var_1420 == i_laysolfr_var_1421(iplon_var_1423)) p_sfluxzen_var_1410(iplon_var_1423, ig_var_1414) = sfluxrefc_var_269(ig_var_1414, js_var_1419) + z_fs_var_1426 * (sfluxrefc_var_269(ig_var_1414, js_var_1419 + 1) - sfluxrefc_var_269(ig_var_1414, js_var_1419))
        p_taur_var_1412(iplon_var_1423, i_lay_var_1420, ig_var_1414) = z_tauray_var_1430
      END DO
    END DO
  END DO
  DO i_lay_var_1420 = laytrop_min_var_1424 + 1, laytrop_max_var_1425
    DO iplon_var_1423 = kidia_var_1389, kfdia_var_1390
      IF (i_lay_var_1420 <= k_laytrop_var_1403(iplon_var_1423)) THEN
        IF (k_jp_var_1396(iplon_var_1423, i_lay_var_1420) < layreffr_var_264 .AND. k_jp_var_1396(iplon_var_1423, i_lay_var_1420 + 1) >= layreffr_var_264) i_laysolfr_var_1421(iplon_var_1423) = MIN(i_lay_var_1420 + 1, k_laytrop_var_1403(iplon_var_1423))
        z_speccomb_var_1427 = p_colh2o_var_1400(iplon_var_1423, i_lay_var_1420) + strrat_var_263 * p_colco2_var_1401(iplon_var_1423, i_lay_var_1420)
        z_specparm_var_1429 = p_colh2o_var_1400(iplon_var_1423, i_lay_var_1420) / z_speccomb_var_1427
        z_specparm_var_1429 = MIN(p_oneminus_var_1399(iplon_var_1423), z_specparm_var_1429)
        z_specmult_var_1428 = 8.0D0 * (z_specparm_var_1429)
        js_var_1419 = 1 + INT(z_specmult_var_1428)
        z_fs_var_1426 = z_specmult_var_1428 - AINT(z_specmult_var_1428)
        ind0_var_1415 = ((k_jp_var_1396(iplon_var_1423, i_lay_var_1420) - 1) * 5 + (k_jt_var_1397(iplon_var_1423, i_lay_var_1420) - 1)) * nspa_var_334(19) + js_var_1419
        ind1_var_1416 = (k_jp_var_1396(iplon_var_1423, i_lay_var_1420) * 5 + (k_jt1_var_1398(iplon_var_1423, i_lay_var_1420) - 1)) * nspa_var_334(19) + js_var_1419
        inds_var_1417 = k_indself_var_1406(iplon_var_1423, i_lay_var_1420)
        indf_var_1418 = k_indfor_var_1409(iplon_var_1423, i_lay_var_1420)
        z_tauray_var_1430 = p_colmol_var_1402(iplon_var_1423, i_lay_var_1420) * rayl_var_262
        DO ig_var_1414 = 1, 8
          p_taug_var_1411(iplon_var_1423, i_lay_var_1420, ig_var_1414) = z_speccomb_var_1427 * ((1.0D0 - z_fs_var_1426) * (absa_var_265(ind0_var_1415, ig_var_1414) * p_fac00_var_1392(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind0_var_1415 + 9, ig_var_1414) * p_fac10_var_1394(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind1_var_1416, ig_var_1414) * p_fac01_var_1393(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind1_var_1416 + 9, ig_var_1414) * p_fac11_var_1395(iplon_var_1423, i_lay_var_1420)) + z_fs_var_1426 * (absa_var_265(ind0_var_1415 + 1, ig_var_1414) * p_fac00_var_1392(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind0_var_1415 + 10, ig_var_1414) * p_fac10_var_1394(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind1_var_1416 + 1, ig_var_1414) * p_fac01_var_1393(iplon_var_1423, i_lay_var_1420) + absa_var_265(ind1_var_1416 + 10, ig_var_1414) * p_fac11_var_1395(iplon_var_1423, i_lay_var_1420))) + p_colh2o_var_1400(iplon_var_1423, i_lay_var_1420) * (p_selffac_var_1404(iplon_var_1423, i_lay_var_1420) * (selfrefc_var_267(inds_var_1417, ig_var_1414) + p_selffrac_var_1405(iplon_var_1423, i_lay_var_1420) * (selfrefc_var_267(inds_var_1417 + 1, ig_var_1414) - selfrefc_var_267(inds_var_1417, ig_var_1414))) + p_forfac_var_1407(iplon_var_1423, i_lay_var_1420) * (forrefc_var_268(indf_var_1418, ig_var_1414) + p_forfrac_var_1408(iplon_var_1423, i_lay_var_1420) * (forrefc_var_268(indf_var_1418 + 1, ig_var_1414) - forrefc_var_268(indf_var_1418, ig_var_1414))))
          IF (i_lay_var_1420 == i_laysolfr_var_1421(iplon_var_1423)) p_sfluxzen_var_1410(iplon_var_1423, ig_var_1414) = sfluxrefc_var_269(ig_var_1414, js_var_1419) + z_fs_var_1426 * (sfluxrefc_var_269(ig_var_1414, js_var_1419 + 1) - sfluxrefc_var_269(ig_var_1414, js_var_1419))
          p_taur_var_1412(iplon_var_1423, i_lay_var_1420, ig_var_1414) = z_tauray_var_1430
        END DO
      ELSE
        ind0_var_1415 = ((k_jp_var_1396(iplon_var_1423, i_lay_var_1420) - 13) * 5 + (k_jt_var_1397(iplon_var_1423, i_lay_var_1420) - 1)) * nspb_var_335(19) + 1
        ind1_var_1416 = ((k_jp_var_1396(iplon_var_1423, i_lay_var_1420) - 12) * 5 + (k_jt1_var_1398(iplon_var_1423, i_lay_var_1420) - 1)) * nspb_var_335(19) + 1
        z_tauray_var_1430 = p_colmol_var_1402(iplon_var_1423, i_lay_var_1420) * rayl_var_262
        DO ig_var_1414 = 1, 8
          p_taug_var_1411(iplon_var_1423, i_lay_var_1420, ig_var_1414) = p_colco2_var_1401(iplon_var_1423, i_lay_var_1420) * (p_fac00_var_1392(iplon_var_1423, i_lay_var_1420) * absb_var_266(ind0_var_1415, ig_var_1414) + p_fac10_var_1394(iplon_var_1423, i_lay_var_1420) * absb_var_266(ind0_var_1415 + 1, ig_var_1414) + p_fac01_var_1393(iplon_var_1423, i_lay_var_1420) * absb_var_266(ind1_var_1416, ig_var_1414) + p_fac11_var_1395(iplon_var_1423, i_lay_var_1420) * absb_var_266(ind1_var_1416 + 1, ig_var_1414))
          p_taur_var_1412(iplon_var_1423, i_lay_var_1420, ig_var_1414) = z_tauray_var_1430
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1420 = laytrop_max_var_1425 + 1, i_nlayers_var_1422
    DO iplon_var_1423 = kidia_var_1389, kfdia_var_1390
      ind0_var_1415 = ((k_jp_var_1396(iplon_var_1423, i_lay_var_1420) - 13) * 5 + (k_jt_var_1397(iplon_var_1423, i_lay_var_1420) - 1)) * nspb_var_335(19) + 1
      ind1_var_1416 = ((k_jp_var_1396(iplon_var_1423, i_lay_var_1420) - 12) * 5 + (k_jt1_var_1398(iplon_var_1423, i_lay_var_1420) - 1)) * nspb_var_335(19) + 1
      z_tauray_var_1430 = p_colmol_var_1402(iplon_var_1423, i_lay_var_1420) * rayl_var_262
      DO ig_var_1414 = 1, 8
        p_taug_var_1411(iplon_var_1423, i_lay_var_1420, ig_var_1414) = p_colco2_var_1401(iplon_var_1423, i_lay_var_1420) * (p_fac00_var_1392(iplon_var_1423, i_lay_var_1420) * absb_var_266(ind0_var_1415, ig_var_1414) + p_fac10_var_1394(iplon_var_1423, i_lay_var_1420) * absb_var_266(ind0_var_1415 + 1, ig_var_1414) + p_fac01_var_1393(iplon_var_1423, i_lay_var_1420) * absb_var_266(ind1_var_1416, ig_var_1414) + p_fac11_var_1395(iplon_var_1423, i_lay_var_1420) * absb_var_266(ind1_var_1416 + 1, ig_var_1414))
        p_taur_var_1412(iplon_var_1423, i_lay_var_1420, ig_var_1414) = z_tauray_var_1430
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol19
SUBROUTINE srtm_taumol21(kidia_var_1431, kfdia_var_1432, klev_var_1433, p_fac00_var_1434, p_fac01_var_1435, p_fac10_var_1436, p_fac11_var_1437, k_jp_var_1438, k_jt_var_1439, k_jt1_var_1440, p_oneminus_var_1441, p_colh2o_var_1442, p_colco2_var_1443, p_colmol_var_1444, k_laytrop_var_1445, p_selffac_var_1446, p_selffrac_var_1447, k_indself_var_1448, p_forfac_var_1449, p_forfrac_var_1450, k_indfor_var_1451, p_sfluxzen_var_1452, p_taug_var_1453, p_taur_var_1454, prmu0_var_1455)
  USE yoesrta21, ONLY: absa_var_280, absb_var_281, forrefc_var_283, layreffr_var_279, rayl_var_277, selfrefc_var_282, sfluxrefc_var_284, strrat_var_278
  USE yoesrtwn, ONLY: nspa_var_334, nspb_var_335
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1431, kfdia_var_1432
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1433
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1434(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1435(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1436(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1437(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1438(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1439(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1440(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1441(kidia_var_1431 : kfdia_var_1432)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1442(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  REAL(KIND = 8), INTENT(IN) :: p_colco2_var_1443(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1444(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1445(kidia_var_1431 : kfdia_var_1432)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1446(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1447(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1448(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1449(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1450(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1451(kidia_var_1431 : kfdia_var_1432, klev_var_1433)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1452(kidia_var_1431 : kfdia_var_1432, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1453(kidia_var_1431 : kfdia_var_1432, klev_var_1433, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1454(kidia_var_1431 : kfdia_var_1432, klev_var_1433, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1455(kidia_var_1431 : kfdia_var_1432)
  INTEGER(KIND = 4) :: ig_var_1456, ind0_var_1457, ind1_var_1458, inds_var_1459, indf_var_1460, js_var_1461, i_lay_var_1462, i_laysolfr_var_1463(kidia_var_1431 : kfdia_var_1432), i_nlayers_var_1464, iplon_var_1465
  INTEGER(KIND = 4) :: laytrop_min_var_1466, laytrop_max_var_1467
  REAL(KIND = 8) :: z_fs_var_1468, z_speccomb_var_1469, z_specmult_var_1470, z_specparm_var_1471, z_tauray_var_1472
  laytrop_min_var_1466 = MINVAL(k_laytrop_var_1445(kidia_var_1431 : kfdia_var_1432))
  laytrop_max_var_1467 = MAXVAL(k_laytrop_var_1445(kidia_var_1431 : kfdia_var_1432))
  i_nlayers_var_1464 = klev_var_1433
  DO iplon_var_1465 = kidia_var_1431, kfdia_var_1432
    i_laysolfr_var_1463(iplon_var_1465) = k_laytrop_var_1445(iplon_var_1465)
  END DO
  DO i_lay_var_1462 = 1, laytrop_min_var_1466
    DO iplon_var_1465 = kidia_var_1431, kfdia_var_1432
      IF (k_jp_var_1438(iplon_var_1465, i_lay_var_1462) < layreffr_var_279 .AND. k_jp_var_1438(iplon_var_1465, i_lay_var_1462 + 1) >= layreffr_var_279) i_laysolfr_var_1463(iplon_var_1465) = MIN(i_lay_var_1462 + 1, k_laytrop_var_1445(iplon_var_1465))
      z_speccomb_var_1469 = p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) + strrat_var_278 * p_colco2_var_1443(iplon_var_1465, i_lay_var_1462)
      z_specparm_var_1471 = p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) / z_speccomb_var_1469
      z_specparm_var_1471 = MIN(p_oneminus_var_1441(iplon_var_1465), z_specparm_var_1471)
      z_specmult_var_1470 = 8.0D0 * (z_specparm_var_1471)
      js_var_1461 = 1 + INT(z_specmult_var_1470)
      z_fs_var_1468 = z_specmult_var_1470 - AINT(z_specmult_var_1470)
      ind0_var_1457 = ((k_jp_var_1438(iplon_var_1465, i_lay_var_1462) - 1) * 5 + (k_jt_var_1439(iplon_var_1465, i_lay_var_1462) - 1)) * nspa_var_334(21) + js_var_1461
      ind1_var_1458 = (k_jp_var_1438(iplon_var_1465, i_lay_var_1462) * 5 + (k_jt1_var_1440(iplon_var_1465, i_lay_var_1462) - 1)) * nspa_var_334(21) + js_var_1461
      inds_var_1459 = k_indself_var_1448(iplon_var_1465, i_lay_var_1462)
      indf_var_1460 = k_indfor_var_1451(iplon_var_1465, i_lay_var_1462)
      z_tauray_var_1472 = p_colmol_var_1444(iplon_var_1465, i_lay_var_1462) * rayl_var_277
      DO ig_var_1456 = 1, 10
        p_taug_var_1453(iplon_var_1465, i_lay_var_1462, ig_var_1456) = z_speccomb_var_1469 * ((1.0D0 - z_fs_var_1468) * (absa_var_280(ind0_var_1457, ig_var_1456) * p_fac00_var_1434(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind0_var_1457 + 9, ig_var_1456) * p_fac10_var_1436(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind1_var_1458, ig_var_1456) * p_fac01_var_1435(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind1_var_1458 + 9, ig_var_1456) * p_fac11_var_1437(iplon_var_1465, i_lay_var_1462)) + z_fs_var_1468 * (absa_var_280(ind0_var_1457 + 1, ig_var_1456) * p_fac00_var_1434(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind0_var_1457 + 10, ig_var_1456) * p_fac10_var_1436(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind1_var_1458 + 1, ig_var_1456) * p_fac01_var_1435(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind1_var_1458 + 10, ig_var_1456) * p_fac11_var_1437(iplon_var_1465, i_lay_var_1462))) + p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) * (p_selffac_var_1446(iplon_var_1465, i_lay_var_1462) * (selfrefc_var_282(inds_var_1459, ig_var_1456) + p_selffrac_var_1447(iplon_var_1465, i_lay_var_1462) * (selfrefc_var_282(inds_var_1459 + 1, ig_var_1456) - selfrefc_var_282(inds_var_1459, ig_var_1456))) + p_forfac_var_1449(iplon_var_1465, i_lay_var_1462) * (forrefc_var_283(indf_var_1460, ig_var_1456) + p_forfrac_var_1450(iplon_var_1465, i_lay_var_1462) * (forrefc_var_283(indf_var_1460 + 1, ig_var_1456) - forrefc_var_283(indf_var_1460, ig_var_1456))))
        IF (i_lay_var_1462 == i_laysolfr_var_1463(iplon_var_1465)) p_sfluxzen_var_1452(iplon_var_1465, ig_var_1456) = sfluxrefc_var_284(ig_var_1456, js_var_1461) + z_fs_var_1468 * (sfluxrefc_var_284(ig_var_1456, js_var_1461 + 1) - sfluxrefc_var_284(ig_var_1456, js_var_1461))
        p_taur_var_1454(iplon_var_1465, i_lay_var_1462, ig_var_1456) = z_tauray_var_1472
      END DO
    END DO
  END DO
  DO i_lay_var_1462 = laytrop_min_var_1466 + 1, laytrop_max_var_1467
    DO iplon_var_1465 = kidia_var_1431, kfdia_var_1432
      IF (i_lay_var_1462 <= k_laytrop_var_1445(iplon_var_1465)) THEN
        IF (k_jp_var_1438(iplon_var_1465, i_lay_var_1462) < layreffr_var_279 .AND. k_jp_var_1438(iplon_var_1465, i_lay_var_1462 + 1) >= layreffr_var_279) i_laysolfr_var_1463(iplon_var_1465) = MIN(i_lay_var_1462 + 1, k_laytrop_var_1445(iplon_var_1465))
        z_speccomb_var_1469 = p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) + strrat_var_278 * p_colco2_var_1443(iplon_var_1465, i_lay_var_1462)
        z_specparm_var_1471 = p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) / z_speccomb_var_1469
        z_specparm_var_1471 = MIN(p_oneminus_var_1441(iplon_var_1465), z_specparm_var_1471)
        z_specmult_var_1470 = 8.0D0 * (z_specparm_var_1471)
        js_var_1461 = 1 + INT(z_specmult_var_1470)
        z_fs_var_1468 = z_specmult_var_1470 - AINT(z_specmult_var_1470)
        ind0_var_1457 = ((k_jp_var_1438(iplon_var_1465, i_lay_var_1462) - 1) * 5 + (k_jt_var_1439(iplon_var_1465, i_lay_var_1462) - 1)) * nspa_var_334(21) + js_var_1461
        ind1_var_1458 = (k_jp_var_1438(iplon_var_1465, i_lay_var_1462) * 5 + (k_jt1_var_1440(iplon_var_1465, i_lay_var_1462) - 1)) * nspa_var_334(21) + js_var_1461
        inds_var_1459 = k_indself_var_1448(iplon_var_1465, i_lay_var_1462)
        indf_var_1460 = k_indfor_var_1451(iplon_var_1465, i_lay_var_1462)
        z_tauray_var_1472 = p_colmol_var_1444(iplon_var_1465, i_lay_var_1462) * rayl_var_277
        DO ig_var_1456 = 1, 10
          p_taug_var_1453(iplon_var_1465, i_lay_var_1462, ig_var_1456) = z_speccomb_var_1469 * ((1.0D0 - z_fs_var_1468) * (absa_var_280(ind0_var_1457, ig_var_1456) * p_fac00_var_1434(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind0_var_1457 + 9, ig_var_1456) * p_fac10_var_1436(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind1_var_1458, ig_var_1456) * p_fac01_var_1435(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind1_var_1458 + 9, ig_var_1456) * p_fac11_var_1437(iplon_var_1465, i_lay_var_1462)) + z_fs_var_1468 * (absa_var_280(ind0_var_1457 + 1, ig_var_1456) * p_fac00_var_1434(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind0_var_1457 + 10, ig_var_1456) * p_fac10_var_1436(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind1_var_1458 + 1, ig_var_1456) * p_fac01_var_1435(iplon_var_1465, i_lay_var_1462) + absa_var_280(ind1_var_1458 + 10, ig_var_1456) * p_fac11_var_1437(iplon_var_1465, i_lay_var_1462))) + p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) * (p_selffac_var_1446(iplon_var_1465, i_lay_var_1462) * (selfrefc_var_282(inds_var_1459, ig_var_1456) + p_selffrac_var_1447(iplon_var_1465, i_lay_var_1462) * (selfrefc_var_282(inds_var_1459 + 1, ig_var_1456) - selfrefc_var_282(inds_var_1459, ig_var_1456))) + p_forfac_var_1449(iplon_var_1465, i_lay_var_1462) * (forrefc_var_283(indf_var_1460, ig_var_1456) + p_forfrac_var_1450(iplon_var_1465, i_lay_var_1462) * (forrefc_var_283(indf_var_1460 + 1, ig_var_1456) - forrefc_var_283(indf_var_1460, ig_var_1456))))
          IF (i_lay_var_1462 == i_laysolfr_var_1463(iplon_var_1465)) p_sfluxzen_var_1452(iplon_var_1465, ig_var_1456) = sfluxrefc_var_284(ig_var_1456, js_var_1461) + z_fs_var_1468 * (sfluxrefc_var_284(ig_var_1456, js_var_1461 + 1) - sfluxrefc_var_284(ig_var_1456, js_var_1461))
          p_taur_var_1454(iplon_var_1465, i_lay_var_1462, ig_var_1456) = z_tauray_var_1472
        END DO
      ELSE
        z_speccomb_var_1469 = p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) + strrat_var_278 * p_colco2_var_1443(iplon_var_1465, i_lay_var_1462)
        z_specparm_var_1471 = p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) / z_speccomb_var_1469
        z_specparm_var_1471 = MIN(p_oneminus_var_1441(iplon_var_1465), z_specparm_var_1471)
        z_specmult_var_1470 = 4.0D0 * (z_specparm_var_1471)
        js_var_1461 = 1 + INT(z_specmult_var_1470)
        z_fs_var_1468 = z_specmult_var_1470 - AINT(z_specmult_var_1470)
        ind0_var_1457 = ((k_jp_var_1438(iplon_var_1465, i_lay_var_1462) - 13) * 5 + (k_jt_var_1439(iplon_var_1465, i_lay_var_1462) - 1)) * nspb_var_335(21) + js_var_1461
        ind1_var_1458 = ((k_jp_var_1438(iplon_var_1465, i_lay_var_1462) - 12) * 5 + (k_jt1_var_1440(iplon_var_1465, i_lay_var_1462) - 1)) * nspb_var_335(21) + js_var_1461
        indf_var_1460 = k_indfor_var_1451(iplon_var_1465, i_lay_var_1462)
        z_tauray_var_1472 = p_colmol_var_1444(iplon_var_1465, i_lay_var_1462) * rayl_var_277
        DO ig_var_1456 = 1, 10
          p_taug_var_1453(iplon_var_1465, i_lay_var_1462, ig_var_1456) = z_speccomb_var_1469 * ((1.0D0 - z_fs_var_1468) * (absb_var_281(ind0_var_1457, ig_var_1456) * p_fac00_var_1434(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind0_var_1457 + 5, ig_var_1456) * p_fac10_var_1436(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind1_var_1458, ig_var_1456) * p_fac01_var_1435(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind1_var_1458 + 5, ig_var_1456) * p_fac11_var_1437(iplon_var_1465, i_lay_var_1462)) + z_fs_var_1468 * (absb_var_281(ind0_var_1457 + 1, ig_var_1456) * p_fac00_var_1434(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind0_var_1457 + 6, ig_var_1456) * p_fac10_var_1436(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind1_var_1458 + 1, ig_var_1456) * p_fac01_var_1435(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind1_var_1458 + 6, ig_var_1456) * p_fac11_var_1437(iplon_var_1465, i_lay_var_1462))) + p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) * p_forfac_var_1449(iplon_var_1465, i_lay_var_1462) * (forrefc_var_283(indf_var_1460, ig_var_1456) + p_forfrac_var_1450(iplon_var_1465, i_lay_var_1462) * (forrefc_var_283(indf_var_1460 + 1, ig_var_1456) - forrefc_var_283(indf_var_1460, ig_var_1456)))
          p_taur_var_1454(iplon_var_1465, i_lay_var_1462, ig_var_1456) = z_tauray_var_1472
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1462 = laytrop_max_var_1467 + 1, i_nlayers_var_1464
    DO iplon_var_1465 = kidia_var_1431, kfdia_var_1432
      z_speccomb_var_1469 = p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) + strrat_var_278 * p_colco2_var_1443(iplon_var_1465, i_lay_var_1462)
      z_specparm_var_1471 = p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) / z_speccomb_var_1469
      z_specparm_var_1471 = MIN(p_oneminus_var_1441(iplon_var_1465), z_specparm_var_1471)
      z_specmult_var_1470 = 4.0D0 * (z_specparm_var_1471)
      js_var_1461 = 1 + INT(z_specmult_var_1470)
      z_fs_var_1468 = z_specmult_var_1470 - AINT(z_specmult_var_1470)
      ind0_var_1457 = ((k_jp_var_1438(iplon_var_1465, i_lay_var_1462) - 13) * 5 + (k_jt_var_1439(iplon_var_1465, i_lay_var_1462) - 1)) * nspb_var_335(21) + js_var_1461
      ind1_var_1458 = ((k_jp_var_1438(iplon_var_1465, i_lay_var_1462) - 12) * 5 + (k_jt1_var_1440(iplon_var_1465, i_lay_var_1462) - 1)) * nspb_var_335(21) + js_var_1461
      indf_var_1460 = k_indfor_var_1451(iplon_var_1465, i_lay_var_1462)
      z_tauray_var_1472 = p_colmol_var_1444(iplon_var_1465, i_lay_var_1462) * rayl_var_277
      DO ig_var_1456 = 1, 10
        p_taug_var_1453(iplon_var_1465, i_lay_var_1462, ig_var_1456) = z_speccomb_var_1469 * ((1.0D0 - z_fs_var_1468) * (absb_var_281(ind0_var_1457, ig_var_1456) * p_fac00_var_1434(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind0_var_1457 + 5, ig_var_1456) * p_fac10_var_1436(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind1_var_1458, ig_var_1456) * p_fac01_var_1435(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind1_var_1458 + 5, ig_var_1456) * p_fac11_var_1437(iplon_var_1465, i_lay_var_1462)) + z_fs_var_1468 * (absb_var_281(ind0_var_1457 + 1, ig_var_1456) * p_fac00_var_1434(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind0_var_1457 + 6, ig_var_1456) * p_fac10_var_1436(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind1_var_1458 + 1, ig_var_1456) * p_fac01_var_1435(iplon_var_1465, i_lay_var_1462) + absb_var_281(ind1_var_1458 + 6, ig_var_1456) * p_fac11_var_1437(iplon_var_1465, i_lay_var_1462))) + p_colh2o_var_1442(iplon_var_1465, i_lay_var_1462) * p_forfac_var_1449(iplon_var_1465, i_lay_var_1462) * (forrefc_var_283(indf_var_1460, ig_var_1456) + p_forfrac_var_1450(iplon_var_1465, i_lay_var_1462) * (forrefc_var_283(indf_var_1460 + 1, ig_var_1456) - forrefc_var_283(indf_var_1460, ig_var_1456)))
        p_taur_var_1454(iplon_var_1465, i_lay_var_1462, ig_var_1456) = z_tauray_var_1472
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol21
SUBROUTINE srtm_taumol20(kidia_var_1473, kfdia_var_1474, klev_var_1475, p_fac00_var_1476, p_fac01_var_1477, p_fac10_var_1478, p_fac11_var_1479, k_jp_var_1480, k_jt_var_1481, k_jt1_var_1482, p_colh2o_var_1483, p_colch4_var_1484, p_colmol_var_1485, k_laytrop_var_1486, p_selffac_var_1487, p_selffrac_var_1488, k_indself_var_1489, p_forfac_var_1490, p_forfrac_var_1491, k_indfor_var_1492, p_sfluxzen_var_1493, p_taug_var_1494, p_taur_var_1495, prmu0_var_1496)
  USE yoesrta20, ONLY: absa_var_272, absb_var_273, absch4c, forrefc_var_275, layreffr_var_271, rayl_var_270, selfrefc_var_274, sfluxrefc_var_276
  USE yoesrtwn, ONLY: nspa_var_334, nspb_var_335
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1473, kfdia_var_1474
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1475
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1476(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1477(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1478(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1479(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1480(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1481(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1482(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1483(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  REAL(KIND = 8), INTENT(IN) :: p_colch4_var_1484(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1485(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1486(kidia_var_1473 : kfdia_var_1474)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1487(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1488(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1489(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1490(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1491(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1492(kidia_var_1473 : kfdia_var_1474, klev_var_1475)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1493(kidia_var_1473 : kfdia_var_1474, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1494(kidia_var_1473 : kfdia_var_1474, klev_var_1475, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1495(kidia_var_1473 : kfdia_var_1474, klev_var_1475, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1496(kidia_var_1473 : kfdia_var_1474)
  INTEGER(KIND = 4) :: ig_var_1497, ind0_var_1498, ind1_var_1499, inds_var_1500, indf_var_1501, i_lay_var_1502, i_laysolfr_var_1503(kidia_var_1473 : kfdia_var_1474), i_nlayers_var_1504, iplon_var_1505
  INTEGER(KIND = 4) :: laytrop_min_var_1506, laytrop_max_var_1507
  REAL(KIND = 8) :: z_tauray_var_1508
  laytrop_min_var_1506 = MINVAL(k_laytrop_var_1486(kidia_var_1473 : kfdia_var_1474))
  laytrop_max_var_1507 = MAXVAL(k_laytrop_var_1486(kidia_var_1473 : kfdia_var_1474))
  i_nlayers_var_1504 = klev_var_1475
  DO iplon_var_1505 = kidia_var_1473, kfdia_var_1474
    i_laysolfr_var_1503(iplon_var_1505) = k_laytrop_var_1486(iplon_var_1505)
  END DO
  DO i_lay_var_1502 = 1, laytrop_min_var_1506
    DO iplon_var_1505 = kidia_var_1473, kfdia_var_1474
      IF (k_jp_var_1480(iplon_var_1505, i_lay_var_1502) < layreffr_var_271 .AND. k_jp_var_1480(iplon_var_1505, i_lay_var_1502 + 1) >= layreffr_var_271) i_laysolfr_var_1503(iplon_var_1505) = MIN(i_lay_var_1502 + 1, k_laytrop_var_1486(iplon_var_1505))
      ind0_var_1498 = ((k_jp_var_1480(iplon_var_1505, i_lay_var_1502) - 1) * 5 + (k_jt_var_1481(iplon_var_1505, i_lay_var_1502) - 1)) * nspa_var_334(20) + 1
      ind1_var_1499 = (k_jp_var_1480(iplon_var_1505, i_lay_var_1502) * 5 + (k_jt1_var_1482(iplon_var_1505, i_lay_var_1502) - 1)) * nspa_var_334(20) + 1
      inds_var_1500 = k_indself_var_1489(iplon_var_1505, i_lay_var_1502)
      indf_var_1501 = k_indfor_var_1492(iplon_var_1505, i_lay_var_1502)
      z_tauray_var_1508 = p_colmol_var_1485(iplon_var_1505, i_lay_var_1502) * rayl_var_270
      DO ig_var_1497 = 1, 10
        p_taug_var_1494(iplon_var_1505, i_lay_var_1502, ig_var_1497) = p_colh2o_var_1483(iplon_var_1505, i_lay_var_1502) * ((p_fac00_var_1476(iplon_var_1505, i_lay_var_1502) * absa_var_272(ind0_var_1498, ig_var_1497) + p_fac10_var_1478(iplon_var_1505, i_lay_var_1502) * absa_var_272(ind0_var_1498 + 1, ig_var_1497) + p_fac01_var_1477(iplon_var_1505, i_lay_var_1502) * absa_var_272(ind1_var_1499, ig_var_1497) + p_fac11_var_1479(iplon_var_1505, i_lay_var_1502) * absa_var_272(ind1_var_1499 + 1, ig_var_1497)) + p_selffac_var_1487(iplon_var_1505, i_lay_var_1502) * (selfrefc_var_274(inds_var_1500, ig_var_1497) + p_selffrac_var_1488(iplon_var_1505, i_lay_var_1502) * (selfrefc_var_274(inds_var_1500 + 1, ig_var_1497) - selfrefc_var_274(inds_var_1500, ig_var_1497))) + p_forfac_var_1490(iplon_var_1505, i_lay_var_1502) * (forrefc_var_275(indf_var_1501, ig_var_1497) + p_forfrac_var_1491(iplon_var_1505, i_lay_var_1502) * (forrefc_var_275(indf_var_1501 + 1, ig_var_1497) - forrefc_var_275(indf_var_1501, ig_var_1497)))) + p_colch4_var_1484(iplon_var_1505, i_lay_var_1502) * absch4c(ig_var_1497)
        p_taur_var_1495(iplon_var_1505, i_lay_var_1502, ig_var_1497) = z_tauray_var_1508
        IF (i_lay_var_1502 == i_laysolfr_var_1503(iplon_var_1505)) p_sfluxzen_var_1493(iplon_var_1505, ig_var_1497) = sfluxrefc_var_276(ig_var_1497)
      END DO
    END DO
  END DO
  DO i_lay_var_1502 = laytrop_min_var_1506 + 1, laytrop_max_var_1507
    DO iplon_var_1505 = kidia_var_1473, kfdia_var_1474
      IF (i_lay_var_1502 <= k_laytrop_var_1486(iplon_var_1505)) THEN
        IF (k_jp_var_1480(iplon_var_1505, i_lay_var_1502) < layreffr_var_271 .AND. k_jp_var_1480(iplon_var_1505, i_lay_var_1502 + 1) >= layreffr_var_271) i_laysolfr_var_1503(iplon_var_1505) = MIN(i_lay_var_1502 + 1, k_laytrop_var_1486(iplon_var_1505))
        ind0_var_1498 = ((k_jp_var_1480(iplon_var_1505, i_lay_var_1502) - 1) * 5 + (k_jt_var_1481(iplon_var_1505, i_lay_var_1502) - 1)) * nspa_var_334(20) + 1
        ind1_var_1499 = (k_jp_var_1480(iplon_var_1505, i_lay_var_1502) * 5 + (k_jt1_var_1482(iplon_var_1505, i_lay_var_1502) - 1)) * nspa_var_334(20) + 1
        inds_var_1500 = k_indself_var_1489(iplon_var_1505, i_lay_var_1502)
        indf_var_1501 = k_indfor_var_1492(iplon_var_1505, i_lay_var_1502)
        z_tauray_var_1508 = p_colmol_var_1485(iplon_var_1505, i_lay_var_1502) * rayl_var_270
        DO ig_var_1497 = 1, 10
          p_taug_var_1494(iplon_var_1505, i_lay_var_1502, ig_var_1497) = p_colh2o_var_1483(iplon_var_1505, i_lay_var_1502) * ((p_fac00_var_1476(iplon_var_1505, i_lay_var_1502) * absa_var_272(ind0_var_1498, ig_var_1497) + p_fac10_var_1478(iplon_var_1505, i_lay_var_1502) * absa_var_272(ind0_var_1498 + 1, ig_var_1497) + p_fac01_var_1477(iplon_var_1505, i_lay_var_1502) * absa_var_272(ind1_var_1499, ig_var_1497) + p_fac11_var_1479(iplon_var_1505, i_lay_var_1502) * absa_var_272(ind1_var_1499 + 1, ig_var_1497)) + p_selffac_var_1487(iplon_var_1505, i_lay_var_1502) * (selfrefc_var_274(inds_var_1500, ig_var_1497) + p_selffrac_var_1488(iplon_var_1505, i_lay_var_1502) * (selfrefc_var_274(inds_var_1500 + 1, ig_var_1497) - selfrefc_var_274(inds_var_1500, ig_var_1497))) + p_forfac_var_1490(iplon_var_1505, i_lay_var_1502) * (forrefc_var_275(indf_var_1501, ig_var_1497) + p_forfrac_var_1491(iplon_var_1505, i_lay_var_1502) * (forrefc_var_275(indf_var_1501 + 1, ig_var_1497) - forrefc_var_275(indf_var_1501, ig_var_1497)))) + p_colch4_var_1484(iplon_var_1505, i_lay_var_1502) * absch4c(ig_var_1497)
          p_taur_var_1495(iplon_var_1505, i_lay_var_1502, ig_var_1497) = z_tauray_var_1508
          IF (i_lay_var_1502 == i_laysolfr_var_1503(iplon_var_1505)) p_sfluxzen_var_1493(iplon_var_1505, ig_var_1497) = sfluxrefc_var_276(ig_var_1497)
        END DO
      ELSE
        ind0_var_1498 = ((k_jp_var_1480(iplon_var_1505, i_lay_var_1502) - 13) * 5 + (k_jt_var_1481(iplon_var_1505, i_lay_var_1502) - 1)) * nspb_var_335(20) + 1
        ind1_var_1499 = ((k_jp_var_1480(iplon_var_1505, i_lay_var_1502) - 12) * 5 + (k_jt1_var_1482(iplon_var_1505, i_lay_var_1502) - 1)) * nspb_var_335(20) + 1
        indf_var_1501 = k_indfor_var_1492(iplon_var_1505, i_lay_var_1502)
        z_tauray_var_1508 = p_colmol_var_1485(iplon_var_1505, i_lay_var_1502) * rayl_var_270
        DO ig_var_1497 = 1, 10
          p_taug_var_1494(iplon_var_1505, i_lay_var_1502, ig_var_1497) = p_colh2o_var_1483(iplon_var_1505, i_lay_var_1502) * (p_fac00_var_1476(iplon_var_1505, i_lay_var_1502) * absb_var_273(ind0_var_1498, ig_var_1497) + p_fac10_var_1478(iplon_var_1505, i_lay_var_1502) * absb_var_273(ind0_var_1498 + 1, ig_var_1497) + p_fac01_var_1477(iplon_var_1505, i_lay_var_1502) * absb_var_273(ind1_var_1499, ig_var_1497) + p_fac11_var_1479(iplon_var_1505, i_lay_var_1502) * absb_var_273(ind1_var_1499 + 1, ig_var_1497) + p_forfac_var_1490(iplon_var_1505, i_lay_var_1502) * (forrefc_var_275(indf_var_1501, ig_var_1497) + p_forfrac_var_1491(iplon_var_1505, i_lay_var_1502) * (forrefc_var_275(indf_var_1501 + 1, ig_var_1497) - forrefc_var_275(indf_var_1501, ig_var_1497)))) + p_colch4_var_1484(iplon_var_1505, i_lay_var_1502) * absch4c(ig_var_1497)
          p_taur_var_1495(iplon_var_1505, i_lay_var_1502, ig_var_1497) = z_tauray_var_1508
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1502 = laytrop_max_var_1507 + 1, i_nlayers_var_1504
    DO iplon_var_1505 = kidia_var_1473, kfdia_var_1474
      ind0_var_1498 = ((k_jp_var_1480(iplon_var_1505, i_lay_var_1502) - 13) * 5 + (k_jt_var_1481(iplon_var_1505, i_lay_var_1502) - 1)) * nspb_var_335(20) + 1
      ind1_var_1499 = ((k_jp_var_1480(iplon_var_1505, i_lay_var_1502) - 12) * 5 + (k_jt1_var_1482(iplon_var_1505, i_lay_var_1502) - 1)) * nspb_var_335(20) + 1
      indf_var_1501 = k_indfor_var_1492(iplon_var_1505, i_lay_var_1502)
      z_tauray_var_1508 = p_colmol_var_1485(iplon_var_1505, i_lay_var_1502) * rayl_var_270
      DO ig_var_1497 = 1, 10
        p_taug_var_1494(iplon_var_1505, i_lay_var_1502, ig_var_1497) = p_colh2o_var_1483(iplon_var_1505, i_lay_var_1502) * (p_fac00_var_1476(iplon_var_1505, i_lay_var_1502) * absb_var_273(ind0_var_1498, ig_var_1497) + p_fac10_var_1478(iplon_var_1505, i_lay_var_1502) * absb_var_273(ind0_var_1498 + 1, ig_var_1497) + p_fac01_var_1477(iplon_var_1505, i_lay_var_1502) * absb_var_273(ind1_var_1499, ig_var_1497) + p_fac11_var_1479(iplon_var_1505, i_lay_var_1502) * absb_var_273(ind1_var_1499 + 1, ig_var_1497) + p_forfac_var_1490(iplon_var_1505, i_lay_var_1502) * (forrefc_var_275(indf_var_1501, ig_var_1497) + p_forfrac_var_1491(iplon_var_1505, i_lay_var_1502) * (forrefc_var_275(indf_var_1501 + 1, ig_var_1497) - forrefc_var_275(indf_var_1501, ig_var_1497)))) + p_colch4_var_1484(iplon_var_1505, i_lay_var_1502) * absch4c(ig_var_1497)
        p_taur_var_1495(iplon_var_1505, i_lay_var_1502, ig_var_1497) = z_tauray_var_1508
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol20
SUBROUTINE srtm_taumol22(kidia_var_1509, kfdia_var_1510, klev_var_1511, p_fac00_var_1512, p_fac01_var_1513, p_fac10_var_1514, p_fac11_var_1515, k_jp_var_1516, k_jt_var_1517, k_jt1_var_1518, p_oneminus_var_1519, p_colh2o_var_1520, p_colmol_var_1521, p_colo2_var_1522, k_laytrop_var_1523, p_selffac_var_1524, p_selffrac_var_1525, k_indself_var_1526, p_forfac_var_1527, p_forfrac_var_1528, k_indfor_var_1529, p_sfluxzen_var_1530, p_taug_var_1531, p_taur_var_1532, prmu0_var_1533)
  USE yoesrta22, ONLY: absa_var_288, absb_var_289, forrefc_var_291, layreffr_var_287, rayl_var_285, selfrefc_var_290, sfluxrefc_var_292, strrat_var_286
  USE yoesrtwn, ONLY: nspa_var_334, nspb_var_335
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1509, kfdia_var_1510
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1511
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1512(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1513(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1514(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1515(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1516(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1517(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1518(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  REAL(KIND = 8), INTENT(IN) :: p_oneminus_var_1519(kidia_var_1509 : kfdia_var_1510)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1520(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1521(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  REAL(KIND = 8), INTENT(IN) :: p_colo2_var_1522(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1523(kidia_var_1509 : kfdia_var_1510)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1524(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1525(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1526(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1527(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1528(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1529(kidia_var_1509 : kfdia_var_1510, klev_var_1511)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1530(kidia_var_1509 : kfdia_var_1510, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1531(kidia_var_1509 : kfdia_var_1510, klev_var_1511, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1532(kidia_var_1509 : kfdia_var_1510, klev_var_1511, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1533(kidia_var_1509 : kfdia_var_1510)
  INTEGER(KIND = 4) :: ig_var_1534, ind0_var_1535, ind1_var_1536, inds_var_1537, indf_var_1538, js_var_1539, i_lay_var_1540, i_laysolfr_var_1541(kidia_var_1509 : kfdia_var_1510), i_nlayers_var_1542, iplon_var_1543
  INTEGER(KIND = 4) :: laytrop_min_var_1544, laytrop_max_var_1545
  REAL(KIND = 8) :: z_fs_var_1546, z_speccomb_var_1547, z_specmult_var_1548, z_specparm_var_1549, z_tauray_var_1550, z_o2adj, z_o2cont
  laytrop_min_var_1544 = MINVAL(k_laytrop_var_1523(kidia_var_1509 : kfdia_var_1510))
  laytrop_max_var_1545 = MAXVAL(k_laytrop_var_1523(kidia_var_1509 : kfdia_var_1510))
  i_nlayers_var_1542 = klev_var_1511
  z_o2adj = 1.6D0
  DO iplon_var_1543 = kidia_var_1509, kfdia_var_1510
    i_laysolfr_var_1541(iplon_var_1543) = k_laytrop_var_1523(iplon_var_1543)
  END DO
  DO i_lay_var_1540 = 1, laytrop_min_var_1544
    DO iplon_var_1543 = kidia_var_1509, kfdia_var_1510
      IF (k_jp_var_1516(iplon_var_1543, i_lay_var_1540) < layreffr_var_287 .AND. k_jp_var_1516(iplon_var_1543, i_lay_var_1540 + 1) >= layreffr_var_287) i_laysolfr_var_1541(iplon_var_1543) = MIN(i_lay_var_1540 + 1, k_laytrop_var_1523(iplon_var_1543))
      z_o2cont = 0.000435D0 * p_colo2_var_1522(iplon_var_1543, i_lay_var_1540) / (700.0D0)
      z_speccomb_var_1547 = p_colh2o_var_1520(iplon_var_1543, i_lay_var_1540) + z_o2adj * strrat_var_286 * p_colo2_var_1522(iplon_var_1543, i_lay_var_1540)
      z_specparm_var_1549 = p_colh2o_var_1520(iplon_var_1543, i_lay_var_1540) / z_speccomb_var_1547
      z_specparm_var_1549 = MIN(p_oneminus_var_1519(iplon_var_1543), z_specparm_var_1549)
      z_specmult_var_1548 = 8.0D0 * (z_specparm_var_1549)
      js_var_1539 = 1 + INT(z_specmult_var_1548)
      z_fs_var_1546 = z_specmult_var_1548 - AINT(z_specmult_var_1548)
      ind0_var_1535 = ((k_jp_var_1516(iplon_var_1543, i_lay_var_1540) - 1) * 5 + (k_jt_var_1517(iplon_var_1543, i_lay_var_1540) - 1)) * nspa_var_334(22) + js_var_1539
      ind1_var_1536 = (k_jp_var_1516(iplon_var_1543, i_lay_var_1540) * 5 + (k_jt1_var_1518(iplon_var_1543, i_lay_var_1540) - 1)) * nspa_var_334(22) + js_var_1539
      inds_var_1537 = k_indself_var_1526(iplon_var_1543, i_lay_var_1540)
      indf_var_1538 = k_indfor_var_1529(iplon_var_1543, i_lay_var_1540)
      z_tauray_var_1550 = p_colmol_var_1521(iplon_var_1543, i_lay_var_1540) * rayl_var_285
      DO ig_var_1534 = 1, 2
        p_taug_var_1531(iplon_var_1543, i_lay_var_1540, ig_var_1534) = z_speccomb_var_1547 * ((1.0D0 - z_fs_var_1546) * (absa_var_288(ind0_var_1535, ig_var_1534) * p_fac00_var_1512(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind0_var_1535 + 9, ig_var_1534) * p_fac10_var_1514(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind1_var_1536, ig_var_1534) * p_fac01_var_1513(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind1_var_1536 + 9, ig_var_1534) * p_fac11_var_1515(iplon_var_1543, i_lay_var_1540)) + z_fs_var_1546 * (absa_var_288(ind0_var_1535 + 1, ig_var_1534) * p_fac00_var_1512(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind0_var_1535 + 10, ig_var_1534) * p_fac10_var_1514(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind1_var_1536 + 1, ig_var_1534) * p_fac01_var_1513(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind1_var_1536 + 10, ig_var_1534) * p_fac11_var_1515(iplon_var_1543, i_lay_var_1540))) + p_colh2o_var_1520(iplon_var_1543, i_lay_var_1540) * (p_selffac_var_1524(iplon_var_1543, i_lay_var_1540) * (selfrefc_var_290(inds_var_1537, ig_var_1534) + p_selffrac_var_1525(iplon_var_1543, i_lay_var_1540) * (selfrefc_var_290(inds_var_1537 + 1, ig_var_1534) - selfrefc_var_290(inds_var_1537, ig_var_1534))) + p_forfac_var_1527(iplon_var_1543, i_lay_var_1540) * (forrefc_var_291(indf_var_1538, ig_var_1534) + p_forfrac_var_1528(iplon_var_1543, i_lay_var_1540) * (forrefc_var_291(indf_var_1538 + 1, ig_var_1534) - forrefc_var_291(indf_var_1538, ig_var_1534)))) + z_o2cont
        IF (i_lay_var_1540 == i_laysolfr_var_1541(iplon_var_1543)) p_sfluxzen_var_1530(iplon_var_1543, ig_var_1534) = sfluxrefc_var_292(ig_var_1534, js_var_1539) + z_fs_var_1546 * (sfluxrefc_var_292(ig_var_1534, js_var_1539 + 1) - sfluxrefc_var_292(ig_var_1534, js_var_1539))
        p_taur_var_1532(iplon_var_1543, i_lay_var_1540, ig_var_1534) = z_tauray_var_1550
      END DO
    END DO
  END DO
  DO i_lay_var_1540 = laytrop_min_var_1544 + 1, laytrop_max_var_1545
    DO iplon_var_1543 = kidia_var_1509, kfdia_var_1510
      IF (i_lay_var_1540 <= k_laytrop_var_1523(iplon_var_1543)) THEN
        IF (k_jp_var_1516(iplon_var_1543, i_lay_var_1540) < layreffr_var_287 .AND. k_jp_var_1516(iplon_var_1543, i_lay_var_1540 + 1) >= layreffr_var_287) i_laysolfr_var_1541(iplon_var_1543) = MIN(i_lay_var_1540 + 1, k_laytrop_var_1523(iplon_var_1543))
        z_o2cont = 0.000435D0 * p_colo2_var_1522(iplon_var_1543, i_lay_var_1540) / (700.0D0)
        z_speccomb_var_1547 = p_colh2o_var_1520(iplon_var_1543, i_lay_var_1540) + z_o2adj * strrat_var_286 * p_colo2_var_1522(iplon_var_1543, i_lay_var_1540)
        z_specparm_var_1549 = p_colh2o_var_1520(iplon_var_1543, i_lay_var_1540) / z_speccomb_var_1547
        z_specparm_var_1549 = MIN(p_oneminus_var_1519(iplon_var_1543), z_specparm_var_1549)
        z_specmult_var_1548 = 8.0D0 * (z_specparm_var_1549)
        js_var_1539 = 1 + INT(z_specmult_var_1548)
        z_fs_var_1546 = z_specmult_var_1548 - AINT(z_specmult_var_1548)
        ind0_var_1535 = ((k_jp_var_1516(iplon_var_1543, i_lay_var_1540) - 1) * 5 + (k_jt_var_1517(iplon_var_1543, i_lay_var_1540) - 1)) * nspa_var_334(22) + js_var_1539
        ind1_var_1536 = (k_jp_var_1516(iplon_var_1543, i_lay_var_1540) * 5 + (k_jt1_var_1518(iplon_var_1543, i_lay_var_1540) - 1)) * nspa_var_334(22) + js_var_1539
        inds_var_1537 = k_indself_var_1526(iplon_var_1543, i_lay_var_1540)
        indf_var_1538 = k_indfor_var_1529(iplon_var_1543, i_lay_var_1540)
        z_tauray_var_1550 = p_colmol_var_1521(iplon_var_1543, i_lay_var_1540) * rayl_var_285
        DO ig_var_1534 = 1, 2
          p_taug_var_1531(iplon_var_1543, i_lay_var_1540, ig_var_1534) = z_speccomb_var_1547 * ((1.0D0 - z_fs_var_1546) * (absa_var_288(ind0_var_1535, ig_var_1534) * p_fac00_var_1512(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind0_var_1535 + 9, ig_var_1534) * p_fac10_var_1514(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind1_var_1536, ig_var_1534) * p_fac01_var_1513(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind1_var_1536 + 9, ig_var_1534) * p_fac11_var_1515(iplon_var_1543, i_lay_var_1540)) + z_fs_var_1546 * (absa_var_288(ind0_var_1535 + 1, ig_var_1534) * p_fac00_var_1512(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind0_var_1535 + 10, ig_var_1534) * p_fac10_var_1514(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind1_var_1536 + 1, ig_var_1534) * p_fac01_var_1513(iplon_var_1543, i_lay_var_1540) + absa_var_288(ind1_var_1536 + 10, ig_var_1534) * p_fac11_var_1515(iplon_var_1543, i_lay_var_1540))) + p_colh2o_var_1520(iplon_var_1543, i_lay_var_1540) * (p_selffac_var_1524(iplon_var_1543, i_lay_var_1540) * (selfrefc_var_290(inds_var_1537, ig_var_1534) + p_selffrac_var_1525(iplon_var_1543, i_lay_var_1540) * (selfrefc_var_290(inds_var_1537 + 1, ig_var_1534) - selfrefc_var_290(inds_var_1537, ig_var_1534))) + p_forfac_var_1527(iplon_var_1543, i_lay_var_1540) * (forrefc_var_291(indf_var_1538, ig_var_1534) + p_forfrac_var_1528(iplon_var_1543, i_lay_var_1540) * (forrefc_var_291(indf_var_1538 + 1, ig_var_1534) - forrefc_var_291(indf_var_1538, ig_var_1534)))) + z_o2cont
          IF (i_lay_var_1540 == i_laysolfr_var_1541(iplon_var_1543)) p_sfluxzen_var_1530(iplon_var_1543, ig_var_1534) = sfluxrefc_var_292(ig_var_1534, js_var_1539) + z_fs_var_1546 * (sfluxrefc_var_292(ig_var_1534, js_var_1539 + 1) - sfluxrefc_var_292(ig_var_1534, js_var_1539))
          p_taur_var_1532(iplon_var_1543, i_lay_var_1540, ig_var_1534) = z_tauray_var_1550
        END DO
      ELSE
        z_o2cont = 0.000435D0 * p_colo2_var_1522(iplon_var_1543, i_lay_var_1540) / (700.0D0)
        ind0_var_1535 = ((k_jp_var_1516(iplon_var_1543, i_lay_var_1540) - 13) * 5 + (k_jt_var_1517(iplon_var_1543, i_lay_var_1540) - 1)) * nspb_var_335(22) + 1
        ind1_var_1536 = ((k_jp_var_1516(iplon_var_1543, i_lay_var_1540) - 12) * 5 + (k_jt1_var_1518(iplon_var_1543, i_lay_var_1540) - 1)) * nspb_var_335(22) + 1
        z_tauray_var_1550 = p_colmol_var_1521(iplon_var_1543, i_lay_var_1540) * rayl_var_285
        DO ig_var_1534 = 1, 2
          p_taug_var_1531(iplon_var_1543, i_lay_var_1540, ig_var_1534) = p_colo2_var_1522(iplon_var_1543, i_lay_var_1540) * z_o2adj * (p_fac00_var_1512(iplon_var_1543, i_lay_var_1540) * absb_var_289(ind0_var_1535, ig_var_1534) + p_fac10_var_1514(iplon_var_1543, i_lay_var_1540) * absb_var_289(ind0_var_1535 + 1, ig_var_1534) + p_fac01_var_1513(iplon_var_1543, i_lay_var_1540) * absb_var_289(ind1_var_1536, ig_var_1534) + p_fac11_var_1515(iplon_var_1543, i_lay_var_1540) * absb_var_289(ind1_var_1536 + 1, ig_var_1534)) + z_o2cont
          p_taur_var_1532(iplon_var_1543, i_lay_var_1540, ig_var_1534) = z_tauray_var_1550
        END DO
      END IF
    END DO
  END DO
  DO i_lay_var_1540 = laytrop_max_var_1545 + 1, i_nlayers_var_1542
    DO iplon_var_1543 = kidia_var_1509, kfdia_var_1510
      z_o2cont = 0.000435D0 * p_colo2_var_1522(iplon_var_1543, i_lay_var_1540) / (700.0D0)
      ind0_var_1535 = ((k_jp_var_1516(iplon_var_1543, i_lay_var_1540) - 13) * 5 + (k_jt_var_1517(iplon_var_1543, i_lay_var_1540) - 1)) * nspb_var_335(22) + 1
      ind1_var_1536 = ((k_jp_var_1516(iplon_var_1543, i_lay_var_1540) - 12) * 5 + (k_jt1_var_1518(iplon_var_1543, i_lay_var_1540) - 1)) * nspb_var_335(22) + 1
      z_tauray_var_1550 = p_colmol_var_1521(iplon_var_1543, i_lay_var_1540) * rayl_var_285
      DO ig_var_1534 = 1, 2
        p_taug_var_1531(iplon_var_1543, i_lay_var_1540, ig_var_1534) = p_colo2_var_1522(iplon_var_1543, i_lay_var_1540) * z_o2adj * (p_fac00_var_1512(iplon_var_1543, i_lay_var_1540) * absb_var_289(ind0_var_1535, ig_var_1534) + p_fac10_var_1514(iplon_var_1543, i_lay_var_1540) * absb_var_289(ind0_var_1535 + 1, ig_var_1534) + p_fac01_var_1513(iplon_var_1543, i_lay_var_1540) * absb_var_289(ind1_var_1536, ig_var_1534) + p_fac11_var_1515(iplon_var_1543, i_lay_var_1540) * absb_var_289(ind1_var_1536 + 1, ig_var_1534)) + z_o2cont
        p_taur_var_1532(iplon_var_1543, i_lay_var_1540, ig_var_1534) = z_tauray_var_1550
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol22
SUBROUTINE srtm_taumol23(kidia_var_1551, kfdia_var_1552, klev_var_1553, p_fac00_var_1554, p_fac01_var_1555, p_fac10_var_1556, p_fac11_var_1557, k_jp_var_1558, k_jt_var_1559, k_jt1_var_1560, p_colh2o_var_1561, p_colmol_var_1562, k_laytrop_var_1563, p_selffac_var_1564, p_selffrac_var_1565, k_indself_var_1566, p_forfac_var_1567, p_forfrac_var_1568, k_indfor_var_1569, p_sfluxzen_var_1570, p_taug_var_1571, p_taur_var_1572, prmu0_var_1573)
  USE yoesrta23, ONLY: absa_var_294, forrefc_var_296, givfac, layreffr_var_293, raylc_var_298, selfrefc_var_295, sfluxrefc_var_297
  USE yoesrtwn, ONLY: nspa_var_334
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1551, kfdia_var_1552
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1553
  REAL(KIND = 8), INTENT(IN) :: p_fac00_var_1554(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  REAL(KIND = 8), INTENT(IN) :: p_fac01_var_1555(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  REAL(KIND = 8), INTENT(IN) :: p_fac10_var_1556(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  REAL(KIND = 8), INTENT(IN) :: p_fac11_var_1557(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  INTEGER(KIND = 4), INTENT(IN) :: k_jp_var_1558(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt_var_1559(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  INTEGER(KIND = 4), INTENT(IN) :: k_jt1_var_1560(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  REAL(KIND = 8), INTENT(IN) :: p_colh2o_var_1561(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  REAL(KIND = 8), INTENT(IN) :: p_colmol_var_1562(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  INTEGER(KIND = 4), INTENT(IN) :: k_laytrop_var_1563(kidia_var_1551 : kfdia_var_1552)
  REAL(KIND = 8), INTENT(IN) :: p_selffac_var_1564(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  REAL(KIND = 8), INTENT(IN) :: p_selffrac_var_1565(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  INTEGER(KIND = 4), INTENT(IN) :: k_indself_var_1566(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  REAL(KIND = 8), INTENT(IN) :: p_forfac_var_1567(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  REAL(KIND = 8), INTENT(IN) :: p_forfrac_var_1568(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  INTEGER(KIND = 4), INTENT(IN) :: k_indfor_var_1569(kidia_var_1551 : kfdia_var_1552, klev_var_1553)
  REAL(KIND = 8), INTENT(INOUT) :: p_sfluxzen_var_1570(kidia_var_1551 : kfdia_var_1552, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taug_var_1571(kidia_var_1551 : kfdia_var_1552, klev_var_1553, 16)
  REAL(KIND = 8), INTENT(INOUT) :: p_taur_var_1572(kidia_var_1551 : kfdia_var_1552, klev_var_1553, 16)
  REAL(KIND = 8), INTENT(IN) :: prmu0_var_1573(kidia_var_1551 : kfdia_var_1552)
  INTEGER(KIND = 4) :: ig_var_1574, ind0_var_1575, ind1_var_1576, inds_var_1577, indf_var_1578, i_lay_var_1579, i_laysolfr_var_1580(kidia_var_1551 : kfdia_var_1552), i_nlayers_var_1581, iplon_var_1582
  INTEGER(KIND = 4) :: laytrop_min_var_1583, laytrop_max_var_1584
  REAL(KIND = 8) :: z_tauray_var_1585
  laytrop_min_var_1583 = MINVAL(k_laytrop_var_1563(kidia_var_1551 : kfdia_var_1552))
  laytrop_max_var_1584 = MAXVAL(k_laytrop_var_1563(kidia_var_1551 : kfdia_var_1552))
  i_nlayers_var_1581 = klev_var_1553
  DO iplon_var_1582 = kidia_var_1551, kfdia_var_1552
    i_laysolfr_var_1580(iplon_var_1582) = k_laytrop_var_1563(iplon_var_1582)
  END DO
  DO i_lay_var_1579 = 1, laytrop_min_var_1583
    DO iplon_var_1582 = kidia_var_1551, kfdia_var_1552
      IF (k_jp_var_1558(iplon_var_1582, i_lay_var_1579) < layreffr_var_293 .AND. k_jp_var_1558(iplon_var_1582, i_lay_var_1579 + 1) >= layreffr_var_293) i_laysolfr_var_1580(iplon_var_1582) = MIN(i_lay_var_1579 + 1, k_laytrop_var_1563(iplon_var_1582))
      ind0_var_1575 = ((k_jp_var_1558(iplon_var_1582, i_lay_var_1579) - 1) * 5 + (k_jt_var_1559(iplon_var_1582, i_lay_var_1579) - 1)) * nspa_var_334(23) + 1
      ind1_var_1576 = (k_jp_var_1558(iplon_var_1582, i_lay_var_1579) * 5 + (k_jt1_var_1560(iplon_var_1582, i_lay_var_1579) - 1)) * nspa_var_334(23) + 1
      inds_var_1577 = k_indself_var_1566(iplon_var_1582, i_lay_var_1579)
      indf_var_1578 = k_indfor_var_1569(iplon_var_1582, i_lay_var_1579)
      DO ig_var_1574 = 1, 10
        z_tauray_var_1585 = p_colmol_var_1562(iplon_var_1582, i_lay_var_1579) * raylc_var_298(ig_var_1574)
        p_taug_var_1571(iplon_var_1582, i_lay_var_1579, ig_var_1574) = p_colh2o_var_1561(iplon_var_1582, i_lay_var_1579) * (givfac * (p_fac00_var_1554(iplon_var_1582, i_lay_var_1579) * absa_var_294(ind0_var_1575, ig_var_1574) + p_fac10_var_1556(iplon_var_1582, i_lay_var_1579) * absa_var_294(ind0_var_1575 + 1, ig_var_1574) + p_fac01_var_1555(iplon_var_1582, i_lay_var_1579) * absa_var_294(ind1_var_1576, ig_var_1574) + p_fac11_var_1557(iplon_var_1582, i_lay_var_1579) * absa_var_294(ind1_var_1576 + 1, ig_var_1574)) + p_selffac_var_1564(iplon_var_1582, i_lay_var_1579) * (selfrefc_var_295(inds_var_1577, ig_var_1574) + p_selffrac_var_1565(iplon_var_1582, i_lay_var_1579) * (selfrefc_var_295(inds_var_1577 + 1, ig_var_1574) - selfrefc_var_295(inds_var_1577, ig_var_1574))) + p_forfac_var_1567(iplon_var_1582, i_lay_var_1579) * (forrefc_var_296(indf_var_1578, ig_var_1574) + p_forfrac_var_1568(iplon_var_1582, i_lay_var_1579) * (forrefc_var_296(indf_var_1578 + 1, ig_var_1574) - forrefc_var_296(indf_var_1578, ig_var_1574))))
        IF (i_lay_var_1579 == i_laysolfr_var_1580(iplon_var_1582)) p_sfluxzen_var_1570(iplon_var_1582, ig_var_1574) = sfluxrefc_var_297(ig_var_1574)
        p_taur_var_1572(iplon_var_1582, i_lay_var_1579, ig_var_1574) = z_tauray_var_1585
      END DO
    END DO
  END DO
  DO i_lay_var_1579 = laytrop_min_var_1583 + 1, laytrop_max_var_1584
    DO iplon_var_1582 = kidia_var_1551, kfdia_var_1552
      IF (i_lay_var_1579 <= k_laytrop_var_1563(iplon_var_1582)) THEN
        IF (k_jp_var_1558(iplon_var_1582, i_lay_var_1579) < layreffr_var_293 .AND. k_jp_var_1558(iplon_var_1582, i_lay_var_1579 + 1) >= layreffr_var_293) i_laysolfr_var_1580(iplon_var_1582) = MIN(i_lay_var_1579 + 1, k_laytrop_var_1563(iplon_var_1582))
        ind0_var_1575 = ((k_jp_var_1558(iplon_var_1582, i_lay_var_1579) - 1) * 5 + (k_jt_var_1559(iplon_var_1582, i_lay_var_1579) - 1)) * nspa_var_334(23) + 1
        ind1_var_1576 = (k_jp_var_1558(iplon_var_1582, i_lay_var_1579) * 5 + (k_jt1_var_1560(iplon_var_1582, i_lay_var_1579) - 1)) * nspa_var_334(23) + 1
        inds_var_1577 = k_indself_var_1566(iplon_var_1582, i_lay_var_1579)
        indf_var_1578 = k_indfor_var_1569(iplon_var_1582, i_lay_var_1579)
        DO ig_var_1574 = 1, 10
          z_tauray_var_1585 = p_colmol_var_1562(iplon_var_1582, i_lay_var_1579) * raylc_var_298(ig_var_1574)
          p_taug_var_1571(iplon_var_1582, i_lay_var_1579, ig_var_1574) = p_colh2o_var_1561(iplon_var_1582, i_lay_var_1579) * (givfac * (p_fac00_var_1554(iplon_var_1582, i_lay_var_1579) * absa_var_294(ind0_var_1575, ig_var_1574) + p_fac10_var_1556(iplon_var_1582, i_lay_var_1579) * absa_var_294(ind0_var_1575 + 1, ig_var_1574) + p_fac01_var_1555(iplon_var_1582, i_lay_var_1579) * absa_var_294(ind1_var_1576, ig_var_1574) + p_fac11_var_1557(iplon_var_1582, i_lay_var_1579) * absa_var_294(ind1_var_1576 + 1, ig_var_1574)) + p_selffac_var_1564(iplon_var_1582, i_lay_var_1579) * (selfrefc_var_295(inds_var_1577, ig_var_1574) + p_selffrac_var_1565(iplon_var_1582, i_lay_var_1579) * (selfrefc_var_295(inds_var_1577 + 1, ig_var_1574) - selfrefc_var_295(inds_var_1577, ig_var_1574))) + p_forfac_var_1567(iplon_var_1582, i_lay_var_1579) * (forrefc_var_296(indf_var_1578, ig_var_1574) + p_forfrac_var_1568(iplon_var_1582, i_lay_var_1579) * (forrefc_var_296(indf_var_1578 + 1, ig_var_1574) - forrefc_var_296(indf_var_1578, ig_var_1574))))
          IF (i_lay_var_1579 == i_laysolfr_var_1580(iplon_var_1582)) p_sfluxzen_var_1570(iplon_var_1582, ig_var_1574) = sfluxrefc_var_297(ig_var_1574)
          p_taur_var_1572(iplon_var_1582, i_lay_var_1579, ig_var_1574) = z_tauray_var_1585
        END DO
      ELSE
        DO ig_var_1574 = 1, 10
          p_taug_var_1571(iplon_var_1582, i_lay_var_1579, ig_var_1574) = 0.0D0
          p_taur_var_1572(iplon_var_1582, i_lay_var_1579, ig_var_1574) = p_colmol_var_1562(iplon_var_1582, i_lay_var_1579) * raylc_var_298(ig_var_1574)
        END DO
      END IF
    END DO
  END DO
  DO ig_var_1574 = 1, 10
    DO i_lay_var_1579 = laytrop_max_var_1584 + 1, i_nlayers_var_1581
      DO iplon_var_1582 = kidia_var_1551, kfdia_var_1552
        p_taug_var_1571(iplon_var_1582, i_lay_var_1579, ig_var_1574) = 0.0D0
        p_taur_var_1572(iplon_var_1582, i_lay_var_1579, ig_var_1574) = p_colmol_var_1562(iplon_var_1582, i_lay_var_1579) * raylc_var_298(ig_var_1574)
      END DO
    END DO
  END DO
END SUBROUTINE srtm_taumol23
SUBROUTINE rrtm_taumol5(kidia_var_1586, kfdia_var_1587, klev_var_1588, taug_var_1589, wx_var_1590, p_tauaerl_var_1591, fac00_var_1592, fac01_var_1593, fac10_var_1594, fac11_var_1595, forfac_var_1614, forfrac_var_1613, indfor_var_1612, jp_var_1596, jt_var_1597, jt1_var_1598, oneminus_var_1599, colh2o_var_1600, colco2_var_1601, colo3_var_1602, laytrop_var_1603, selffac_var_1604, selffrac_var_1605, indself_var_1606, fracs_var_1607, rat_h2oco2_var_1608, rat_h2oco2_1_var_1609, rat_o3co2_var_1610, rat_o3co2_1_var_1611, minorfrac_var_1615, indminor_var_1616)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrtm, ONLY: ng5
  USE yoerrta5, ONLY: absa_var_196, absb_var_197, ccl4, forref_var_200, fracrefa_var_194, fracrefb_var_195, ka_mo3_var_198, selfref_var_199
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1586
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1587
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1588
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1589(kidia_var_1586 : kfdia_var_1587, 140, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: wx_var_1590(kidia_var_1586 : kfdia_var_1587, 4, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1591(kidia_var_1586 : kfdia_var_1587, klev_var_1588, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1592(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1593(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1594(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1595(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1596(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1597(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1598(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1599
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1600(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1601(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1602(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1603(kidia_var_1586 : kfdia_var_1587)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1604(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1605(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1606(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1607(kidia_var_1586 : kfdia_var_1587, 140, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1608(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1609(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_var_1610(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_1_var_1611(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1612(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1613(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1614(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1615(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1616(kidia_var_1586 : kfdia_var_1587, klev_var_1588)
  REAL(KIND = 8) :: speccomb_var_1617, speccomb1_var_1618, speccomb_mo3, speccomb_planck_var_1619
  INTEGER(KIND = 4) :: ind0_var_1620, ind1_var_1621, inds_var_1622, indf_var_1623, indm_var_1624
  INTEGER(KIND = 4) :: ig_var_1625, js_var_1626, lay_var_1627, js1_var_1628, jpl_var_1629, jmo3
  REAL(KIND = 8) :: refrat_planck_a_var_1630, refrat_planck_b_var_1631, refrat_m_a_var_1632
  REAL(KIND = 8) :: fac000_var_1633, fac100_var_1634, fac200_var_1635, fac010_var_1636, fac110_var_1637, fac210_var_1638, fac001_var_1639, fac101_var_1640, fac201_var_1641, fac011_var_1642, fac111_var_1643, fac211_var_1644
  REAL(KIND = 8) :: p_var_1645, p4_var_1646, fk0_var_1647, fk1_var_1648, fk2_var_1649
  REAL(KIND = 8) :: taufor_var_1650, tauself_var_1651, tau_major_var_1652(16), tau_major1_var_1653(16), o3m1, o3m2, abso3_var_1654
  REAL(KIND = 8) :: fs_var_1655, specmult_var_1656, specparm_var_1657, fs1_var_1658, specmult1_var_1659, specparm1_var_1660, fpl_var_1661, specmult_planck_var_1662, specparm_planck_var_1663, fmo3, specmult_mo3, specparm_mo3
  INTEGER(KIND = 4) :: laytrop_min_var_1664, laytrop_max_var_1665
  INTEGER(KIND = 4) :: ixc_var_1666(klev_var_1588), ixlow_var_1667(kfdia_var_1587, klev_var_1588), ixhigh_var_1668(kfdia_var_1587, klev_var_1588)
  INTEGER(KIND = 4) :: ich_var_1669, icl_var_1670, ixc0_var_1671, ixp_var_1672, jc_var_1673, jl_var_1674
  laytrop_min_var_1664 = MINVAL(laytrop_var_1603)
  laytrop_max_var_1665 = MAXVAL(laytrop_var_1603)
  ixlow_var_1667 = 0
  ixhigh_var_1668 = 0
  ixc_var_1666 = 0
  DO lay_var_1627 = laytrop_min_var_1664 + 1, laytrop_max_var_1665
    icl_var_1670 = 0
    ich_var_1669 = 0
    DO jc_var_1673 = kidia_var_1586, kfdia_var_1587
      IF (lay_var_1627 <= laytrop_var_1603(jc_var_1673)) THEN
        icl_var_1670 = icl_var_1670 + 1
        ixlow_var_1667(icl_var_1670, lay_var_1627) = jc_var_1673
      ELSE
        ich_var_1669 = ich_var_1669 + 1
        ixhigh_var_1668(ich_var_1669, lay_var_1627) = jc_var_1673
      END IF
    END DO
    ixc_var_1666(lay_var_1627) = icl_var_1670
  END DO
  refrat_planck_a_var_1630 = chi_mls(1, 5) / chi_mls(2, 5)
  refrat_planck_b_var_1631 = chi_mls(3, 43) / chi_mls(2, 43)
  refrat_m_a_var_1632 = chi_mls(1, 7) / chi_mls(2, 7)
  DO lay_var_1627 = 1, laytrop_min_var_1664
    DO jl_var_1674 = kidia_var_1586, kfdia_var_1587
      speccomb_var_1617 = colh2o_var_1600(jl_var_1674, lay_var_1627) + rat_h2oco2_var_1608(jl_var_1674, lay_var_1627) * colco2_var_1601(jl_var_1674, lay_var_1627)
      specparm_var_1657 = MIN(colh2o_var_1600(jl_var_1674, lay_var_1627) / speccomb_var_1617, oneminus_var_1599)
      specmult_var_1656 = 8.0D0 * (specparm_var_1657)
      js_var_1626 = 1 + INT(specmult_var_1656)
      fs_var_1655 = ((specmult_var_1656) - AINT((specmult_var_1656)))
      speccomb1_var_1618 = colh2o_var_1600(jl_var_1674, lay_var_1627) + rat_h2oco2_1_var_1609(jl_var_1674, lay_var_1627) * colco2_var_1601(jl_var_1674, lay_var_1627)
      specparm1_var_1660 = MIN(colh2o_var_1600(jl_var_1674, lay_var_1627) / speccomb1_var_1618, oneminus_var_1599)
      specmult1_var_1659 = 8.0D0 * (specparm1_var_1660)
      js1_var_1628 = 1 + INT(specmult1_var_1659)
      fs1_var_1658 = ((specmult1_var_1659) - AINT((specmult1_var_1659)))
      speccomb_mo3 = colh2o_var_1600(jl_var_1674, lay_var_1627) + refrat_m_a_var_1632 * colco2_var_1601(jl_var_1674, lay_var_1627)
      specparm_mo3 = MIN(colh2o_var_1600(jl_var_1674, lay_var_1627) / speccomb_mo3, oneminus_var_1599)
      specmult_mo3 = 8.0D0 * specparm_mo3
      jmo3 = 1 + INT(specmult_mo3)
      fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))
      speccomb_planck_var_1619 = colh2o_var_1600(jl_var_1674, lay_var_1627) + refrat_planck_a_var_1630 * colco2_var_1601(jl_var_1674, lay_var_1627)
      specparm_planck_var_1663 = MIN(colh2o_var_1600(jl_var_1674, lay_var_1627) / speccomb_planck_var_1619, oneminus_var_1599)
      specmult_planck_var_1662 = 8.0D0 * specparm_planck_var_1663
      jpl_var_1629 = 1 + INT(specmult_planck_var_1662)
      fpl_var_1661 = ((specmult_planck_var_1662) - AINT((specmult_planck_var_1662)))
      ind0_var_1620 = ((jp_var_1596(jl_var_1674, lay_var_1627) - 1) * 5 + (jt_var_1597(jl_var_1674, lay_var_1627) - 1)) * nspa_var_237(5) + js_var_1626
      ind1_var_1621 = (jp_var_1596(jl_var_1674, lay_var_1627) * 5 + (jt1_var_1598(jl_var_1674, lay_var_1627) - 1)) * nspa_var_237(5) + js1_var_1628
      inds_var_1622 = indself_var_1606(jl_var_1674, lay_var_1627)
      indf_var_1623 = indfor_var_1612(jl_var_1674, lay_var_1627)
      indm_var_1624 = indminor_var_1616(jl_var_1674, lay_var_1627)
      IF (specparm_var_1657 .LT. 0.125D0) THEN
        p_var_1645 = fs_var_1655 - 1.0D0
        p4_var_1646 = p_var_1645 ** 4
        fk0_var_1647 = p4_var_1646
        fk1_var_1648 = 1.0D0 - p_var_1645 - 2.0D0 * p4_var_1646
        fk2_var_1649 = p_var_1645 + p4_var_1646
        fac000_var_1633 = fk0_var_1647 * fac00_var_1592(jl_var_1674, lay_var_1627)
        fac100_var_1634 = fk1_var_1648 * fac00_var_1592(jl_var_1674, lay_var_1627)
        fac200_var_1635 = fk2_var_1649 * fac00_var_1592(jl_var_1674, lay_var_1627)
        fac010_var_1636 = fk0_var_1647 * fac10_var_1594(jl_var_1674, lay_var_1627)
        fac110_var_1637 = fk1_var_1648 * fac10_var_1594(jl_var_1674, lay_var_1627)
        fac210_var_1638 = fk2_var_1649 * fac10_var_1594(jl_var_1674, lay_var_1627)
      ELSE IF (specparm_var_1657 .GT. 0.875D0) THEN
        p_var_1645 = - fs_var_1655
        p4_var_1646 = p_var_1645 ** 4
        fk0_var_1647 = p4_var_1646
        fk1_var_1648 = 1.0D0 - p_var_1645 - 2.0D0 * p4_var_1646
        fk2_var_1649 = p_var_1645 + p4_var_1646
        fac000_var_1633 = fk0_var_1647 * fac00_var_1592(jl_var_1674, lay_var_1627)
        fac100_var_1634 = fk1_var_1648 * fac00_var_1592(jl_var_1674, lay_var_1627)
        fac200_var_1635 = fk2_var_1649 * fac00_var_1592(jl_var_1674, lay_var_1627)
        fac010_var_1636 = fk0_var_1647 * fac10_var_1594(jl_var_1674, lay_var_1627)
        fac110_var_1637 = fk1_var_1648 * fac10_var_1594(jl_var_1674, lay_var_1627)
        fac210_var_1638 = fk2_var_1649 * fac10_var_1594(jl_var_1674, lay_var_1627)
      ELSE
        fac000_var_1633 = (1.0D0 - fs_var_1655) * fac00_var_1592(jl_var_1674, lay_var_1627)
        fac010_var_1636 = (1.0D0 - fs_var_1655) * fac10_var_1594(jl_var_1674, lay_var_1627)
        fac100_var_1634 = fs_var_1655 * fac00_var_1592(jl_var_1674, lay_var_1627)
        fac110_var_1637 = fs_var_1655 * fac10_var_1594(jl_var_1674, lay_var_1627)
        fac200_var_1635 = 0.0D0
        fac210_var_1638 = 0.0D0
      END IF
      IF (specparm1_var_1660 .LT. 0.125D0) THEN
        p_var_1645 = fs1_var_1658 - 1.0D0
        p4_var_1646 = p_var_1645 ** 4
        fk0_var_1647 = p4_var_1646
        fk1_var_1648 = 1.0D0 - p_var_1645 - 2.0D0 * p4_var_1646
        fk2_var_1649 = p_var_1645 + p4_var_1646
        fac001_var_1639 = fk0_var_1647 * fac01_var_1593(jl_var_1674, lay_var_1627)
        fac101_var_1640 = fk1_var_1648 * fac01_var_1593(jl_var_1674, lay_var_1627)
        fac201_var_1641 = fk2_var_1649 * fac01_var_1593(jl_var_1674, lay_var_1627)
        fac011_var_1642 = fk0_var_1647 * fac11_var_1595(jl_var_1674, lay_var_1627)
        fac111_var_1643 = fk1_var_1648 * fac11_var_1595(jl_var_1674, lay_var_1627)
        fac211_var_1644 = fk2_var_1649 * fac11_var_1595(jl_var_1674, lay_var_1627)
      ELSE IF (specparm1_var_1660 .GT. 0.875D0) THEN
        p_var_1645 = - fs1_var_1658
        p4_var_1646 = p_var_1645 ** 4
        fk0_var_1647 = p4_var_1646
        fk1_var_1648 = 1.0D0 - p_var_1645 - 2.0D0 * p4_var_1646
        fk2_var_1649 = p_var_1645 + p4_var_1646
        fac001_var_1639 = fk0_var_1647 * fac01_var_1593(jl_var_1674, lay_var_1627)
        fac101_var_1640 = fk1_var_1648 * fac01_var_1593(jl_var_1674, lay_var_1627)
        fac201_var_1641 = fk2_var_1649 * fac01_var_1593(jl_var_1674, lay_var_1627)
        fac011_var_1642 = fk0_var_1647 * fac11_var_1595(jl_var_1674, lay_var_1627)
        fac111_var_1643 = fk1_var_1648 * fac11_var_1595(jl_var_1674, lay_var_1627)
        fac211_var_1644 = fk2_var_1649 * fac11_var_1595(jl_var_1674, lay_var_1627)
      ELSE
        fac001_var_1639 = (1.0D0 - fs1_var_1658) * fac01_var_1593(jl_var_1674, lay_var_1627)
        fac011_var_1642 = (1.0D0 - fs1_var_1658) * fac11_var_1595(jl_var_1674, lay_var_1627)
        fac101_var_1640 = fs1_var_1658 * fac01_var_1593(jl_var_1674, lay_var_1627)
        fac111_var_1643 = fs1_var_1658 * fac11_var_1595(jl_var_1674, lay_var_1627)
        fac201_var_1641 = 0.0D0
        fac211_var_1644 = 0.0D0
      END IF
      IF (specparm_var_1657 .LT. 0.125D0) THEN
        tau_major_var_1652(1 : ng5) = speccomb_var_1617 * (fac000_var_1633 * absa_var_196(ind0_var_1620, 1 : 16) + fac100_var_1634 * absa_var_196(ind0_var_1620 + 1, 1 : 16) + fac200_var_1635 * absa_var_196(ind0_var_1620 + 2, 1 : 16) + fac010_var_1636 * absa_var_196(ind0_var_1620 + 9, 1 : 16) + fac110_var_1637 * absa_var_196(ind0_var_1620 + 10, 1 : 16) + fac210_var_1638 * absa_var_196(ind0_var_1620 + 11, 1 : 16))
      ELSE IF (specparm_var_1657 .GT. 0.875D0) THEN
        tau_major_var_1652(1 : ng5) = speccomb_var_1617 * (fac200_var_1635 * absa_var_196(ind0_var_1620 - 1, 1 : 16) + fac100_var_1634 * absa_var_196(ind0_var_1620, 1 : 16) + fac000_var_1633 * absa_var_196(ind0_var_1620 + 1, 1 : 16) + fac210_var_1638 * absa_var_196(ind0_var_1620 + 8, 1 : 16) + fac110_var_1637 * absa_var_196(ind0_var_1620 + 9, 1 : 16) + fac010_var_1636 * absa_var_196(ind0_var_1620 + 10, 1 : 16))
      ELSE
        tau_major_var_1652(1 : ng5) = speccomb_var_1617 * (fac000_var_1633 * absa_var_196(ind0_var_1620, 1 : 16) + fac100_var_1634 * absa_var_196(ind0_var_1620 + 1, 1 : 16) + fac010_var_1636 * absa_var_196(ind0_var_1620 + 9, 1 : 16) + fac110_var_1637 * absa_var_196(ind0_var_1620 + 10, 1 : 16))
      END IF
      IF (specparm1_var_1660 .LT. 0.125D0) THEN
        tau_major1_var_1653(1 : ng5) = speccomb1_var_1618 * (fac001_var_1639 * absa_var_196(ind1_var_1621, 1 : 16) + fac101_var_1640 * absa_var_196(ind1_var_1621 + 1, 1 : 16) + fac201_var_1641 * absa_var_196(ind1_var_1621 + 2, 1 : 16) + fac011_var_1642 * absa_var_196(ind1_var_1621 + 9, 1 : 16) + fac111_var_1643 * absa_var_196(ind1_var_1621 + 10, 1 : 16) + fac211_var_1644 * absa_var_196(ind1_var_1621 + 11, 1 : 16))
      ELSE IF (specparm1_var_1660 .GT. 0.875D0) THEN
        tau_major1_var_1653(1 : ng5) = speccomb1_var_1618 * (fac201_var_1641 * absa_var_196(ind1_var_1621 - 1, 1 : 16) + fac101_var_1640 * absa_var_196(ind1_var_1621, 1 : 16) + fac001_var_1639 * absa_var_196(ind1_var_1621 + 1, 1 : 16) + fac211_var_1644 * absa_var_196(ind1_var_1621 + 8, 1 : 16) + fac111_var_1643 * absa_var_196(ind1_var_1621 + 9, 1 : 16) + fac011_var_1642 * absa_var_196(ind1_var_1621 + 10, 1 : 16))
      ELSE
        tau_major1_var_1653(1 : ng5) = speccomb1_var_1618 * (fac001_var_1639 * absa_var_196(ind1_var_1621, 1 : 16) + fac101_var_1640 * absa_var_196(ind1_var_1621 + 1, 1 : 16) + fac011_var_1642 * absa_var_196(ind1_var_1621 + 9, 1 : 16) + fac111_var_1643 * absa_var_196(ind1_var_1621 + 10, 1 : 16))
      END IF
      DO ig_var_1625 = 1, 16
        tauself_var_1651 = selffac_var_1604(jl_var_1674, lay_var_1627) * (selfref_var_199(inds_var_1622, ig_var_1625) + selffrac_var_1605(jl_var_1674, lay_var_1627) * (selfref_var_199(inds_var_1622 + 1, ig_var_1625) - selfref_var_199(inds_var_1622, ig_var_1625)))
        taufor_var_1650 = forfac_var_1614(jl_var_1674, lay_var_1627) * (forref_var_200(indf_var_1623, ig_var_1625) + forfrac_var_1613(jl_var_1674, lay_var_1627) * (forref_var_200(indf_var_1623 + 1, ig_var_1625) - forref_var_200(indf_var_1623, ig_var_1625)))
        o3m1 = ka_mo3_var_198(jmo3, indm_var_1624, ig_var_1625) + fmo3 * (ka_mo3_var_198(jmo3 + 1, indm_var_1624, ig_var_1625) - ka_mo3_var_198(jmo3, indm_var_1624, ig_var_1625))
        o3m2 = ka_mo3_var_198(jmo3, indm_var_1624 + 1, ig_var_1625) + fmo3 * (ka_mo3_var_198(jmo3 + 1, indm_var_1624 + 1, ig_var_1625) - ka_mo3_var_198(jmo3, indm_var_1624 + 1, ig_var_1625))
        abso3_var_1654 = o3m1 + minorfrac_var_1615(jl_var_1674, lay_var_1627) * (o3m2 - o3m1)
        taug_var_1589(jl_var_1674, 52 + ig_var_1625, lay_var_1627) = tau_major_var_1652(ig_var_1625) + tau_major1_var_1653(ig_var_1625) + tauself_var_1651 + taufor_var_1650 + abso3_var_1654 * colo3_var_1602(jl_var_1674, lay_var_1627) + wx_var_1590(jl_var_1674, 1, lay_var_1627) * ccl4(ig_var_1625)
        fracs_var_1607(jl_var_1674, 52 + ig_var_1625, lay_var_1627) = fracrefa_var_194(ig_var_1625, jpl_var_1629) + fpl_var_1661 * (fracrefa_var_194(ig_var_1625, jpl_var_1629 + 1) - fracrefa_var_194(ig_var_1625, jpl_var_1629))
      END DO
    END DO
  END DO
  DO lay_var_1627 = laytrop_max_var_1665 + 1, klev_var_1588
    DO jl_var_1674 = kidia_var_1586, kfdia_var_1587
      speccomb_var_1617 = colo3_var_1602(jl_var_1674, lay_var_1627) + rat_o3co2_var_1610(jl_var_1674, lay_var_1627) * colco2_var_1601(jl_var_1674, lay_var_1627)
      specparm_var_1657 = MIN(colo3_var_1602(jl_var_1674, lay_var_1627) / speccomb_var_1617, oneminus_var_1599)
      specmult_var_1656 = 4.0D0 * (specparm_var_1657)
      js_var_1626 = 1 + INT(specmult_var_1656)
      fs_var_1655 = ((specmult_var_1656) - AINT((specmult_var_1656)))
      speccomb1_var_1618 = colo3_var_1602(jl_var_1674, lay_var_1627) + rat_o3co2_1_var_1611(jl_var_1674, lay_var_1627) * colco2_var_1601(jl_var_1674, lay_var_1627)
      specparm1_var_1660 = MIN(colo3_var_1602(jl_var_1674, lay_var_1627) / speccomb1_var_1618, oneminus_var_1599)
      specmult1_var_1659 = 4.0D0 * (specparm1_var_1660)
      js1_var_1628 = 1 + INT(specmult1_var_1659)
      fs1_var_1658 = ((specmult1_var_1659) - AINT((specmult1_var_1659)))
      fac000_var_1633 = (1.0D0 - fs_var_1655) * fac00_var_1592(jl_var_1674, lay_var_1627)
      fac010_var_1636 = (1.0D0 - fs_var_1655) * fac10_var_1594(jl_var_1674, lay_var_1627)
      fac100_var_1634 = fs_var_1655 * fac00_var_1592(jl_var_1674, lay_var_1627)
      fac110_var_1637 = fs_var_1655 * fac10_var_1594(jl_var_1674, lay_var_1627)
      fac001_var_1639 = (1.0D0 - fs1_var_1658) * fac01_var_1593(jl_var_1674, lay_var_1627)
      fac011_var_1642 = (1.0D0 - fs1_var_1658) * fac11_var_1595(jl_var_1674, lay_var_1627)
      fac101_var_1640 = fs1_var_1658 * fac01_var_1593(jl_var_1674, lay_var_1627)
      fac111_var_1643 = fs1_var_1658 * fac11_var_1595(jl_var_1674, lay_var_1627)
      speccomb_planck_var_1619 = colo3_var_1602(jl_var_1674, lay_var_1627) + refrat_planck_b_var_1631 * colco2_var_1601(jl_var_1674, lay_var_1627)
      specparm_planck_var_1663 = MIN(colo3_var_1602(jl_var_1674, lay_var_1627) / speccomb_planck_var_1619, oneminus_var_1599)
      specmult_planck_var_1662 = 4.0D0 * specparm_planck_var_1663
      jpl_var_1629 = 1 + INT(specmult_planck_var_1662)
      fpl_var_1661 = ((specmult_planck_var_1662) - AINT((specmult_planck_var_1662)))
      ind0_var_1620 = ((jp_var_1596(jl_var_1674, lay_var_1627) - 13) * 5 + (jt_var_1597(jl_var_1674, lay_var_1627) - 1)) * nspb_var_238(5) + js_var_1626
      ind1_var_1621 = ((jp_var_1596(jl_var_1674, lay_var_1627) - 12) * 5 + (jt1_var_1598(jl_var_1674, lay_var_1627) - 1)) * nspb_var_238(5) + js1_var_1628
      DO ig_var_1625 = 1, 16
        taug_var_1589(jl_var_1674, 52 + ig_var_1625, lay_var_1627) = speccomb_var_1617 * (fac000_var_1633 * absb_var_197(ind0_var_1620, ig_var_1625) + fac100_var_1634 * absb_var_197(ind0_var_1620 + 1, ig_var_1625) + fac010_var_1636 * absb_var_197(ind0_var_1620 + 5, ig_var_1625) + fac110_var_1637 * absb_var_197(ind0_var_1620 + 6, ig_var_1625)) + speccomb1_var_1618 * (fac001_var_1639 * absb_var_197(ind1_var_1621, ig_var_1625) + fac101_var_1640 * absb_var_197(ind1_var_1621 + 1, ig_var_1625) + fac011_var_1642 * absb_var_197(ind1_var_1621 + 5, ig_var_1625) + fac111_var_1643 * absb_var_197(ind1_var_1621 + 6, ig_var_1625)) + wx_var_1590(jl_var_1674, 1, lay_var_1627) * ccl4(ig_var_1625)
        fracs_var_1607(jl_var_1674, 52 + ig_var_1625, lay_var_1627) = fracrefb_var_195(ig_var_1625, jpl_var_1629) + fpl_var_1661 * (fracrefb_var_195(ig_var_1625, jpl_var_1629 + 1) - fracrefb_var_195(ig_var_1625, jpl_var_1629))
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1665 /= laytrop_min_var_1664) THEN
    DO lay_var_1627 = laytrop_min_var_1664 + 1, laytrop_max_var_1665
      ixc0_var_1671 = ixc_var_1666(lay_var_1627)
      DO ixp_var_1672 = 1, ixc0_var_1671
        jl_var_1674 = ixlow_var_1667(ixp_var_1672, lay_var_1627)
        speccomb_var_1617 = colh2o_var_1600(jl_var_1674, lay_var_1627) + rat_h2oco2_var_1608(jl_var_1674, lay_var_1627) * colco2_var_1601(jl_var_1674, lay_var_1627)
        specparm_var_1657 = MIN(colh2o_var_1600(jl_var_1674, lay_var_1627) / speccomb_var_1617, oneminus_var_1599)
        specmult_var_1656 = 8.0D0 * (specparm_var_1657)
        js_var_1626 = 1 + INT(specmult_var_1656)
        fs_var_1655 = ((specmult_var_1656) - AINT((specmult_var_1656)))
        speccomb1_var_1618 = colh2o_var_1600(jl_var_1674, lay_var_1627) + rat_h2oco2_1_var_1609(jl_var_1674, lay_var_1627) * colco2_var_1601(jl_var_1674, lay_var_1627)
        specparm1_var_1660 = MIN(colh2o_var_1600(jl_var_1674, lay_var_1627) / speccomb1_var_1618, oneminus_var_1599)
        specmult1_var_1659 = 8.0D0 * (specparm1_var_1660)
        js1_var_1628 = 1 + INT(specmult1_var_1659)
        fs1_var_1658 = ((specmult1_var_1659) - AINT((specmult1_var_1659)))
        speccomb_mo3 = colh2o_var_1600(jl_var_1674, lay_var_1627) + refrat_m_a_var_1632 * colco2_var_1601(jl_var_1674, lay_var_1627)
        specparm_mo3 = MIN(colh2o_var_1600(jl_var_1674, lay_var_1627) / speccomb_mo3, oneminus_var_1599)
        specmult_mo3 = 8.0D0 * specparm_mo3
        jmo3 = 1 + INT(specmult_mo3)
        fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))
        speccomb_planck_var_1619 = colh2o_var_1600(jl_var_1674, lay_var_1627) + refrat_planck_a_var_1630 * colco2_var_1601(jl_var_1674, lay_var_1627)
        specparm_planck_var_1663 = MIN(colh2o_var_1600(jl_var_1674, lay_var_1627) / speccomb_planck_var_1619, oneminus_var_1599)
        specmult_planck_var_1662 = 8.0D0 * specparm_planck_var_1663
        jpl_var_1629 = 1 + INT(specmult_planck_var_1662)
        fpl_var_1661 = ((specmult_planck_var_1662) - AINT((specmult_planck_var_1662)))
        ind0_var_1620 = ((jp_var_1596(jl_var_1674, lay_var_1627) - 1) * 5 + (jt_var_1597(jl_var_1674, lay_var_1627) - 1)) * nspa_var_237(5) + js_var_1626
        ind1_var_1621 = (jp_var_1596(jl_var_1674, lay_var_1627) * 5 + (jt1_var_1598(jl_var_1674, lay_var_1627) - 1)) * nspa_var_237(5) + js1_var_1628
        inds_var_1622 = indself_var_1606(jl_var_1674, lay_var_1627)
        indf_var_1623 = indfor_var_1612(jl_var_1674, lay_var_1627)
        indm_var_1624 = indminor_var_1616(jl_var_1674, lay_var_1627)
        IF (specparm_var_1657 .LT. 0.125D0) THEN
          p_var_1645 = fs_var_1655 - 1.0D0
          p4_var_1646 = p_var_1645 ** 4
          fk0_var_1647 = p4_var_1646
          fk1_var_1648 = 1.0D0 - p_var_1645 - 2.0D0 * p4_var_1646
          fk2_var_1649 = p_var_1645 + p4_var_1646
          fac000_var_1633 = fk0_var_1647 * fac00_var_1592(jl_var_1674, lay_var_1627)
          fac100_var_1634 = fk1_var_1648 * fac00_var_1592(jl_var_1674, lay_var_1627)
          fac200_var_1635 = fk2_var_1649 * fac00_var_1592(jl_var_1674, lay_var_1627)
          fac010_var_1636 = fk0_var_1647 * fac10_var_1594(jl_var_1674, lay_var_1627)
          fac110_var_1637 = fk1_var_1648 * fac10_var_1594(jl_var_1674, lay_var_1627)
          fac210_var_1638 = fk2_var_1649 * fac10_var_1594(jl_var_1674, lay_var_1627)
        ELSE IF (specparm_var_1657 .GT. 0.875D0) THEN
          p_var_1645 = - fs_var_1655
          p4_var_1646 = p_var_1645 ** 4
          fk0_var_1647 = p4_var_1646
          fk1_var_1648 = 1.0D0 - p_var_1645 - 2.0D0 * p4_var_1646
          fk2_var_1649 = p_var_1645 + p4_var_1646
          fac000_var_1633 = fk0_var_1647 * fac00_var_1592(jl_var_1674, lay_var_1627)
          fac100_var_1634 = fk1_var_1648 * fac00_var_1592(jl_var_1674, lay_var_1627)
          fac200_var_1635 = fk2_var_1649 * fac00_var_1592(jl_var_1674, lay_var_1627)
          fac010_var_1636 = fk0_var_1647 * fac10_var_1594(jl_var_1674, lay_var_1627)
          fac110_var_1637 = fk1_var_1648 * fac10_var_1594(jl_var_1674, lay_var_1627)
          fac210_var_1638 = fk2_var_1649 * fac10_var_1594(jl_var_1674, lay_var_1627)
        ELSE
          fac000_var_1633 = (1.0D0 - fs_var_1655) * fac00_var_1592(jl_var_1674, lay_var_1627)
          fac010_var_1636 = (1.0D0 - fs_var_1655) * fac10_var_1594(jl_var_1674, lay_var_1627)
          fac100_var_1634 = fs_var_1655 * fac00_var_1592(jl_var_1674, lay_var_1627)
          fac110_var_1637 = fs_var_1655 * fac10_var_1594(jl_var_1674, lay_var_1627)
          fac200_var_1635 = 0.0D0
          fac210_var_1638 = 0.0D0
        END IF
        IF (specparm1_var_1660 .LT. 0.125D0) THEN
          p_var_1645 = fs1_var_1658 - 1.0D0
          p4_var_1646 = p_var_1645 ** 4
          fk0_var_1647 = p4_var_1646
          fk1_var_1648 = 1.0D0 - p_var_1645 - 2.0D0 * p4_var_1646
          fk2_var_1649 = p_var_1645 + p4_var_1646
          fac001_var_1639 = fk0_var_1647 * fac01_var_1593(jl_var_1674, lay_var_1627)
          fac101_var_1640 = fk1_var_1648 * fac01_var_1593(jl_var_1674, lay_var_1627)
          fac201_var_1641 = fk2_var_1649 * fac01_var_1593(jl_var_1674, lay_var_1627)
          fac011_var_1642 = fk0_var_1647 * fac11_var_1595(jl_var_1674, lay_var_1627)
          fac111_var_1643 = fk1_var_1648 * fac11_var_1595(jl_var_1674, lay_var_1627)
          fac211_var_1644 = fk2_var_1649 * fac11_var_1595(jl_var_1674, lay_var_1627)
        ELSE IF (specparm1_var_1660 .GT. 0.875D0) THEN
          p_var_1645 = - fs1_var_1658
          p4_var_1646 = p_var_1645 ** 4
          fk0_var_1647 = p4_var_1646
          fk1_var_1648 = 1.0D0 - p_var_1645 - 2.0D0 * p4_var_1646
          fk2_var_1649 = p_var_1645 + p4_var_1646
          fac001_var_1639 = fk0_var_1647 * fac01_var_1593(jl_var_1674, lay_var_1627)
          fac101_var_1640 = fk1_var_1648 * fac01_var_1593(jl_var_1674, lay_var_1627)
          fac201_var_1641 = fk2_var_1649 * fac01_var_1593(jl_var_1674, lay_var_1627)
          fac011_var_1642 = fk0_var_1647 * fac11_var_1595(jl_var_1674, lay_var_1627)
          fac111_var_1643 = fk1_var_1648 * fac11_var_1595(jl_var_1674, lay_var_1627)
          fac211_var_1644 = fk2_var_1649 * fac11_var_1595(jl_var_1674, lay_var_1627)
        ELSE
          fac001_var_1639 = (1.0D0 - fs1_var_1658) * fac01_var_1593(jl_var_1674, lay_var_1627)
          fac011_var_1642 = (1.0D0 - fs1_var_1658) * fac11_var_1595(jl_var_1674, lay_var_1627)
          fac101_var_1640 = fs1_var_1658 * fac01_var_1593(jl_var_1674, lay_var_1627)
          fac111_var_1643 = fs1_var_1658 * fac11_var_1595(jl_var_1674, lay_var_1627)
          fac201_var_1641 = 0.0D0
          fac211_var_1644 = 0.0D0
        END IF
        IF (specparm_var_1657 .LT. 0.125D0) THEN
          tau_major_var_1652(1 : ng5) = speccomb_var_1617 * (fac000_var_1633 * absa_var_196(ind0_var_1620, 1 : 16) + fac100_var_1634 * absa_var_196(ind0_var_1620 + 1, 1 : 16) + fac200_var_1635 * absa_var_196(ind0_var_1620 + 2, 1 : 16) + fac010_var_1636 * absa_var_196(ind0_var_1620 + 9, 1 : 16) + fac110_var_1637 * absa_var_196(ind0_var_1620 + 10, 1 : 16) + fac210_var_1638 * absa_var_196(ind0_var_1620 + 11, 1 : 16))
        ELSE IF (specparm_var_1657 .GT. 0.875D0) THEN
          tau_major_var_1652(1 : ng5) = speccomb_var_1617 * (fac200_var_1635 * absa_var_196(ind0_var_1620 - 1, 1 : 16) + fac100_var_1634 * absa_var_196(ind0_var_1620, 1 : 16) + fac000_var_1633 * absa_var_196(ind0_var_1620 + 1, 1 : 16) + fac210_var_1638 * absa_var_196(ind0_var_1620 + 8, 1 : 16) + fac110_var_1637 * absa_var_196(ind0_var_1620 + 9, 1 : 16) + fac010_var_1636 * absa_var_196(ind0_var_1620 + 10, 1 : 16))
        ELSE
          tau_major_var_1652(1 : ng5) = speccomb_var_1617 * (fac000_var_1633 * absa_var_196(ind0_var_1620, 1 : 16) + fac100_var_1634 * absa_var_196(ind0_var_1620 + 1, 1 : 16) + fac010_var_1636 * absa_var_196(ind0_var_1620 + 9, 1 : 16) + fac110_var_1637 * absa_var_196(ind0_var_1620 + 10, 1 : 16))
        END IF
        IF (specparm1_var_1660 .LT. 0.125D0) THEN
          tau_major1_var_1653(1 : ng5) = speccomb1_var_1618 * (fac001_var_1639 * absa_var_196(ind1_var_1621, 1 : 16) + fac101_var_1640 * absa_var_196(ind1_var_1621 + 1, 1 : 16) + fac201_var_1641 * absa_var_196(ind1_var_1621 + 2, 1 : 16) + fac011_var_1642 * absa_var_196(ind1_var_1621 + 9, 1 : 16) + fac111_var_1643 * absa_var_196(ind1_var_1621 + 10, 1 : 16) + fac211_var_1644 * absa_var_196(ind1_var_1621 + 11, 1 : 16))
        ELSE IF (specparm1_var_1660 .GT. 0.875D0) THEN
          tau_major1_var_1653(1 : ng5) = speccomb1_var_1618 * (fac201_var_1641 * absa_var_196(ind1_var_1621 - 1, 1 : 16) + fac101_var_1640 * absa_var_196(ind1_var_1621, 1 : 16) + fac001_var_1639 * absa_var_196(ind1_var_1621 + 1, 1 : 16) + fac211_var_1644 * absa_var_196(ind1_var_1621 + 8, 1 : 16) + fac111_var_1643 * absa_var_196(ind1_var_1621 + 9, 1 : 16) + fac011_var_1642 * absa_var_196(ind1_var_1621 + 10, 1 : 16))
        ELSE
          tau_major1_var_1653(1 : ng5) = speccomb1_var_1618 * (fac001_var_1639 * absa_var_196(ind1_var_1621, 1 : 16) + fac101_var_1640 * absa_var_196(ind1_var_1621 + 1, 1 : 16) + fac011_var_1642 * absa_var_196(ind1_var_1621 + 9, 1 : 16) + fac111_var_1643 * absa_var_196(ind1_var_1621 + 10, 1 : 16))
        END IF
        DO ig_var_1625 = 1, 16
          tauself_var_1651 = selffac_var_1604(jl_var_1674, lay_var_1627) * (selfref_var_199(inds_var_1622, ig_var_1625) + selffrac_var_1605(jl_var_1674, lay_var_1627) * (selfref_var_199(inds_var_1622 + 1, ig_var_1625) - selfref_var_199(inds_var_1622, ig_var_1625)))
          taufor_var_1650 = forfac_var_1614(jl_var_1674, lay_var_1627) * (forref_var_200(indf_var_1623, ig_var_1625) + forfrac_var_1613(jl_var_1674, lay_var_1627) * (forref_var_200(indf_var_1623 + 1, ig_var_1625) - forref_var_200(indf_var_1623, ig_var_1625)))
          o3m1 = ka_mo3_var_198(jmo3, indm_var_1624, ig_var_1625) + fmo3 * (ka_mo3_var_198(jmo3 + 1, indm_var_1624, ig_var_1625) - ka_mo3_var_198(jmo3, indm_var_1624, ig_var_1625))
          o3m2 = ka_mo3_var_198(jmo3, indm_var_1624 + 1, ig_var_1625) + fmo3 * (ka_mo3_var_198(jmo3 + 1, indm_var_1624 + 1, ig_var_1625) - ka_mo3_var_198(jmo3, indm_var_1624 + 1, ig_var_1625))
          abso3_var_1654 = o3m1 + minorfrac_var_1615(jl_var_1674, lay_var_1627) * (o3m2 - o3m1)
          taug_var_1589(jl_var_1674, 52 + ig_var_1625, lay_var_1627) = tau_major_var_1652(ig_var_1625) + tau_major1_var_1653(ig_var_1625) + tauself_var_1651 + taufor_var_1650 + abso3_var_1654 * colo3_var_1602(jl_var_1674, lay_var_1627) + wx_var_1590(jl_var_1674, 1, lay_var_1627) * ccl4(ig_var_1625)
          fracs_var_1607(jl_var_1674, 52 + ig_var_1625, lay_var_1627) = fracrefa_var_194(ig_var_1625, jpl_var_1629) + fpl_var_1661 * (fracrefa_var_194(ig_var_1625, jpl_var_1629 + 1) - fracrefa_var_194(ig_var_1625, jpl_var_1629))
        END DO
      END DO
      ixc0_var_1671 = kfdia_var_1587 - kidia_var_1586 + 1 - ixc0_var_1671
      DO ixp_var_1672 = 1, ixc0_var_1671
        jl_var_1674 = ixhigh_var_1668(ixp_var_1672, lay_var_1627)
        speccomb_var_1617 = colo3_var_1602(jl_var_1674, lay_var_1627) + rat_o3co2_var_1610(jl_var_1674, lay_var_1627) * colco2_var_1601(jl_var_1674, lay_var_1627)
        specparm_var_1657 = MIN(colo3_var_1602(jl_var_1674, lay_var_1627) / speccomb_var_1617, oneminus_var_1599)
        specmult_var_1656 = 4.0D0 * (specparm_var_1657)
        js_var_1626 = 1 + INT(specmult_var_1656)
        fs_var_1655 = ((specmult_var_1656) - AINT((specmult_var_1656)))
        speccomb1_var_1618 = colo3_var_1602(jl_var_1674, lay_var_1627) + rat_o3co2_1_var_1611(jl_var_1674, lay_var_1627) * colco2_var_1601(jl_var_1674, lay_var_1627)
        specparm1_var_1660 = MIN(colo3_var_1602(jl_var_1674, lay_var_1627) / speccomb1_var_1618, oneminus_var_1599)
        specmult1_var_1659 = 4.0D0 * (specparm1_var_1660)
        js1_var_1628 = 1 + INT(specmult1_var_1659)
        fs1_var_1658 = ((specmult1_var_1659) - AINT((specmult1_var_1659)))
        fac000_var_1633 = (1.0D0 - fs_var_1655) * fac00_var_1592(jl_var_1674, lay_var_1627)
        fac010_var_1636 = (1.0D0 - fs_var_1655) * fac10_var_1594(jl_var_1674, lay_var_1627)
        fac100_var_1634 = fs_var_1655 * fac00_var_1592(jl_var_1674, lay_var_1627)
        fac110_var_1637 = fs_var_1655 * fac10_var_1594(jl_var_1674, lay_var_1627)
        fac001_var_1639 = (1.0D0 - fs1_var_1658) * fac01_var_1593(jl_var_1674, lay_var_1627)
        fac011_var_1642 = (1.0D0 - fs1_var_1658) * fac11_var_1595(jl_var_1674, lay_var_1627)
        fac101_var_1640 = fs1_var_1658 * fac01_var_1593(jl_var_1674, lay_var_1627)
        fac111_var_1643 = fs1_var_1658 * fac11_var_1595(jl_var_1674, lay_var_1627)
        speccomb_planck_var_1619 = colo3_var_1602(jl_var_1674, lay_var_1627) + refrat_planck_b_var_1631 * colco2_var_1601(jl_var_1674, lay_var_1627)
        specparm_planck_var_1663 = MIN(colo3_var_1602(jl_var_1674, lay_var_1627) / speccomb_planck_var_1619, oneminus_var_1599)
        specmult_planck_var_1662 = 4.0D0 * specparm_planck_var_1663
        jpl_var_1629 = 1 + INT(specmult_planck_var_1662)
        fpl_var_1661 = ((specmult_planck_var_1662) - AINT((specmult_planck_var_1662)))
        ind0_var_1620 = ((jp_var_1596(jl_var_1674, lay_var_1627) - 13) * 5 + (jt_var_1597(jl_var_1674, lay_var_1627) - 1)) * nspb_var_238(5) + js_var_1626
        ind1_var_1621 = ((jp_var_1596(jl_var_1674, lay_var_1627) - 12) * 5 + (jt1_var_1598(jl_var_1674, lay_var_1627) - 1)) * nspb_var_238(5) + js1_var_1628
        DO ig_var_1625 = 1, 16
          taug_var_1589(jl_var_1674, 52 + ig_var_1625, lay_var_1627) = speccomb_var_1617 * (fac000_var_1633 * absb_var_197(ind0_var_1620, ig_var_1625) + fac100_var_1634 * absb_var_197(ind0_var_1620 + 1, ig_var_1625) + fac010_var_1636 * absb_var_197(ind0_var_1620 + 5, ig_var_1625) + fac110_var_1637 * absb_var_197(ind0_var_1620 + 6, ig_var_1625)) + speccomb1_var_1618 * (fac001_var_1639 * absb_var_197(ind1_var_1621, ig_var_1625) + fac101_var_1640 * absb_var_197(ind1_var_1621 + 1, ig_var_1625) + fac011_var_1642 * absb_var_197(ind1_var_1621 + 5, ig_var_1625) + fac111_var_1643 * absb_var_197(ind1_var_1621 + 6, ig_var_1625)) + wx_var_1590(jl_var_1674, 1, lay_var_1627) * ccl4(ig_var_1625)
          fracs_var_1607(jl_var_1674, 52 + ig_var_1625, lay_var_1627) = fracrefb_var_195(ig_var_1625, jpl_var_1629) + fpl_var_1661 * (fracrefb_var_195(ig_var_1625, jpl_var_1629 + 1) - fracrefb_var_195(ig_var_1625, jpl_var_1629))
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol5
SUBROUTINE rrtm_taumol4(kidia_var_1675, kfdia_var_1676, klev_var_1677, taug_var_1678, p_tauaerl_var_1679, fac00_var_1680, fac01_var_1681, fac10_var_1682, fac11_var_1683, forfac_var_1701, forfrac_var_1702, indfor_var_1700, jp_var_1684, jt_var_1685, jt1_var_1686, oneminus_var_1687, colh2o_var_1688, colco2_var_1689, colo3_var_1690, laytrop_var_1691, selffac_var_1692, selffrac_var_1693, indself_var_1694, fracs_var_1695, rat_h2oco2_var_1696, rat_h2oco2_1_var_1697, rat_o3co2_var_1698, rat_o3co2_1_var_1699)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrtm, ONLY: ng4
  USE yoerrta4, ONLY: absa_var_190, absb_var_191, forref_var_193, fracrefa_var_188, fracrefb_var_189, selfref_var_192
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1675
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1676
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1677
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1678(kidia_var_1675 : kfdia_var_1676, 140, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1679(kidia_var_1675 : kfdia_var_1676, klev_var_1677, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1680(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1681(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1682(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1683(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1684(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1685(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1686(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1687
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1688(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1689(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1690(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1691(kidia_var_1675 : kfdia_var_1676)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1692(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1693(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1694(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1695(kidia_var_1675 : kfdia_var_1676, 140, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1696(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1697(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_var_1698(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: rat_o3co2_1_var_1699(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1700(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1701(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1702(kidia_var_1675 : kfdia_var_1676, klev_var_1677)
  REAL(KIND = 8) :: speccomb_var_1703, speccomb1_var_1704, speccomb_planck_var_1705
  INTEGER(KIND = 4) :: ind0_var_1706, ind1_var_1707, inds_var_1708, indf_var_1709
  INTEGER(KIND = 4) :: ig_var_1710, js_var_1711, lay_var_1712, js1_var_1713, jpl_var_1714
  REAL(KIND = 8) :: refrat_planck_a_var_1715, refrat_planck_b_var_1716
  REAL(KIND = 8) :: fac000_var_1717, fac100_var_1718, fac200_var_1719, fac010_var_1720, fac110_var_1721, fac210_var_1722, fac001_var_1723, fac101_var_1724, fac201_var_1725, fac011_var_1726, fac111_var_1727, fac211_var_1728
  REAL(KIND = 8) :: p_var_1729, p4_var_1730, fk0_var_1731, fk1_var_1732, fk2_var_1733
  REAL(KIND = 8) :: taufor_var_1734, tauself_var_1735, tau_major_var_1736(14), tau_major1_var_1737(14)
  REAL(KIND = 8) :: fs_var_1738, specmult_var_1739, specparm_var_1740, fs1_var_1741, specmult1_var_1742, specparm1_var_1743, fpl_var_1744, specmult_planck_var_1745, specparm_planck_var_1746
  INTEGER(KIND = 4) :: laytrop_min_var_1747, laytrop_max_var_1748
  INTEGER(KIND = 4) :: ixc_var_1749(klev_var_1677), ixlow_var_1750(kfdia_var_1676, klev_var_1677), ixhigh_var_1751(kfdia_var_1676, klev_var_1677)
  INTEGER(KIND = 4) :: ich_var_1752, icl_var_1753, ixc0_var_1754, ixp_var_1755, jc_var_1756, jl_var_1757
  laytrop_min_var_1747 = MINVAL(laytrop_var_1691)
  laytrop_max_var_1748 = MAXVAL(laytrop_var_1691)
  ixlow_var_1750 = 0
  ixhigh_var_1751 = 0
  ixc_var_1749 = 0
  DO lay_var_1712 = laytrop_min_var_1747 + 1, laytrop_max_var_1748
    icl_var_1753 = 0
    ich_var_1752 = 0
    DO jc_var_1756 = kidia_var_1675, kfdia_var_1676
      IF (lay_var_1712 <= laytrop_var_1691(jc_var_1756)) THEN
        icl_var_1753 = icl_var_1753 + 1
        ixlow_var_1750(icl_var_1753, lay_var_1712) = jc_var_1756
      ELSE
        ich_var_1752 = ich_var_1752 + 1
        ixhigh_var_1751(ich_var_1752, lay_var_1712) = jc_var_1756
      END IF
    END DO
    ixc_var_1749(lay_var_1712) = icl_var_1753
  END DO
  refrat_planck_a_var_1715 = chi_mls(1, 11) / chi_mls(2, 11)
  refrat_planck_b_var_1716 = chi_mls(3, 13) / chi_mls(2, 13)
  DO lay_var_1712 = 1, laytrop_min_var_1747
    DO jl_var_1757 = kidia_var_1675, kfdia_var_1676
      speccomb_var_1703 = colh2o_var_1688(jl_var_1757, lay_var_1712) + rat_h2oco2_var_1696(jl_var_1757, lay_var_1712) * colco2_var_1689(jl_var_1757, lay_var_1712)
      specparm_var_1740 = MIN(colh2o_var_1688(jl_var_1757, lay_var_1712) / speccomb_var_1703, oneminus_var_1687)
      specmult_var_1739 = 8.0D0 * (specparm_var_1740)
      js_var_1711 = 1 + INT(specmult_var_1739)
      fs_var_1738 = ((specmult_var_1739) - AINT((specmult_var_1739)))
      speccomb1_var_1704 = colh2o_var_1688(jl_var_1757, lay_var_1712) + rat_h2oco2_1_var_1697(jl_var_1757, lay_var_1712) * colco2_var_1689(jl_var_1757, lay_var_1712)
      specparm1_var_1743 = MIN(colh2o_var_1688(jl_var_1757, lay_var_1712) / speccomb1_var_1704, oneminus_var_1687)
      specmult1_var_1742 = 8.0D0 * (specparm1_var_1743)
      js1_var_1713 = 1 + INT(specmult1_var_1742)
      fs1_var_1741 = ((specmult1_var_1742) - AINT((specmult1_var_1742)))
      speccomb_planck_var_1705 = colh2o_var_1688(jl_var_1757, lay_var_1712) + refrat_planck_a_var_1715 * colco2_var_1689(jl_var_1757, lay_var_1712)
      specparm_planck_var_1746 = MIN(colh2o_var_1688(jl_var_1757, lay_var_1712) / speccomb_planck_var_1705, oneminus_var_1687)
      specmult_planck_var_1745 = 8.0D0 * specparm_planck_var_1746
      jpl_var_1714 = 1 + INT(specmult_planck_var_1745)
      fpl_var_1744 = ((specmult_planck_var_1745) - AINT((specmult_planck_var_1745)))
      ind0_var_1706 = ((jp_var_1684(jl_var_1757, lay_var_1712) - 1) * 5 + (jt_var_1685(jl_var_1757, lay_var_1712) - 1)) * nspa_var_237(4) + js_var_1711
      ind1_var_1707 = (jp_var_1684(jl_var_1757, lay_var_1712) * 5 + (jt1_var_1686(jl_var_1757, lay_var_1712) - 1)) * nspa_var_237(4) + js1_var_1713
      inds_var_1708 = indself_var_1694(jl_var_1757, lay_var_1712)
      indf_var_1709 = indfor_var_1700(jl_var_1757, lay_var_1712)
      IF (specparm_var_1740 .LT. 0.125D0) THEN
        p_var_1729 = fs_var_1738 - 1.0D0
        p4_var_1730 = p_var_1729 ** 4
        fk0_var_1731 = p4_var_1730
        fk1_var_1732 = 1.0D0 - p_var_1729 - 2.0D0 * p4_var_1730
        fk2_var_1733 = p_var_1729 + p4_var_1730
        fac000_var_1717 = fk0_var_1731 * fac00_var_1680(jl_var_1757, lay_var_1712)
        fac100_var_1718 = fk1_var_1732 * fac00_var_1680(jl_var_1757, lay_var_1712)
        fac200_var_1719 = fk2_var_1733 * fac00_var_1680(jl_var_1757, lay_var_1712)
        fac010_var_1720 = fk0_var_1731 * fac10_var_1682(jl_var_1757, lay_var_1712)
        fac110_var_1721 = fk1_var_1732 * fac10_var_1682(jl_var_1757, lay_var_1712)
        fac210_var_1722 = fk2_var_1733 * fac10_var_1682(jl_var_1757, lay_var_1712)
      ELSE IF (specparm_var_1740 .GT. 0.875D0) THEN
        p_var_1729 = - fs_var_1738
        p4_var_1730 = p_var_1729 ** 4
        fk0_var_1731 = p4_var_1730
        fk1_var_1732 = 1.0D0 - p_var_1729 - 2.0D0 * p4_var_1730
        fk2_var_1733 = p_var_1729 + p4_var_1730
        fac000_var_1717 = fk0_var_1731 * fac00_var_1680(jl_var_1757, lay_var_1712)
        fac100_var_1718 = fk1_var_1732 * fac00_var_1680(jl_var_1757, lay_var_1712)
        fac200_var_1719 = fk2_var_1733 * fac00_var_1680(jl_var_1757, lay_var_1712)
        fac010_var_1720 = fk0_var_1731 * fac10_var_1682(jl_var_1757, lay_var_1712)
        fac110_var_1721 = fk1_var_1732 * fac10_var_1682(jl_var_1757, lay_var_1712)
        fac210_var_1722 = fk2_var_1733 * fac10_var_1682(jl_var_1757, lay_var_1712)
      ELSE
        fac000_var_1717 = (1.0D0 - fs_var_1738) * fac00_var_1680(jl_var_1757, lay_var_1712)
        fac010_var_1720 = (1.0D0 - fs_var_1738) * fac10_var_1682(jl_var_1757, lay_var_1712)
        fac100_var_1718 = fs_var_1738 * fac00_var_1680(jl_var_1757, lay_var_1712)
        fac110_var_1721 = fs_var_1738 * fac10_var_1682(jl_var_1757, lay_var_1712)
        fac200_var_1719 = 0.0D0
        fac210_var_1722 = 0.0D0
      END IF
      IF (specparm1_var_1743 .LT. 0.125D0) THEN
        p_var_1729 = fs1_var_1741 - 1.0D0
        p4_var_1730 = p_var_1729 ** 4
        fk0_var_1731 = p4_var_1730
        fk1_var_1732 = 1.0D0 - p_var_1729 - 2.0D0 * p4_var_1730
        fk2_var_1733 = p_var_1729 + p4_var_1730
        fac001_var_1723 = fk0_var_1731 * fac01_var_1681(jl_var_1757, lay_var_1712)
        fac101_var_1724 = fk1_var_1732 * fac01_var_1681(jl_var_1757, lay_var_1712)
        fac201_var_1725 = fk2_var_1733 * fac01_var_1681(jl_var_1757, lay_var_1712)
        fac011_var_1726 = fk0_var_1731 * fac11_var_1683(jl_var_1757, lay_var_1712)
        fac111_var_1727 = fk1_var_1732 * fac11_var_1683(jl_var_1757, lay_var_1712)
        fac211_var_1728 = fk2_var_1733 * fac11_var_1683(jl_var_1757, lay_var_1712)
      ELSE IF (specparm1_var_1743 .GT. 0.875D0) THEN
        p_var_1729 = - fs1_var_1741
        p4_var_1730 = p_var_1729 ** 4
        fk0_var_1731 = p4_var_1730
        fk1_var_1732 = 1.0D0 - p_var_1729 - 2.0D0 * p4_var_1730
        fk2_var_1733 = p_var_1729 + p4_var_1730
        fac001_var_1723 = fk0_var_1731 * fac01_var_1681(jl_var_1757, lay_var_1712)
        fac101_var_1724 = fk1_var_1732 * fac01_var_1681(jl_var_1757, lay_var_1712)
        fac201_var_1725 = fk2_var_1733 * fac01_var_1681(jl_var_1757, lay_var_1712)
        fac011_var_1726 = fk0_var_1731 * fac11_var_1683(jl_var_1757, lay_var_1712)
        fac111_var_1727 = fk1_var_1732 * fac11_var_1683(jl_var_1757, lay_var_1712)
        fac211_var_1728 = fk2_var_1733 * fac11_var_1683(jl_var_1757, lay_var_1712)
      ELSE
        fac001_var_1723 = (1.0D0 - fs1_var_1741) * fac01_var_1681(jl_var_1757, lay_var_1712)
        fac011_var_1726 = (1.0D0 - fs1_var_1741) * fac11_var_1683(jl_var_1757, lay_var_1712)
        fac101_var_1724 = fs1_var_1741 * fac01_var_1681(jl_var_1757, lay_var_1712)
        fac111_var_1727 = fs1_var_1741 * fac11_var_1683(jl_var_1757, lay_var_1712)
        fac201_var_1725 = 0.0D0
        fac211_var_1728 = 0.0D0
      END IF
      IF (specparm_var_1740 .LT. 0.125D0) THEN
        tau_major_var_1736(1 : ng4) = speccomb_var_1703 * (fac000_var_1717 * absa_var_190(ind0_var_1706, 1 : 14) + fac100_var_1718 * absa_var_190(ind0_var_1706 + 1, 1 : 14) + fac200_var_1719 * absa_var_190(ind0_var_1706 + 2, 1 : 14) + fac010_var_1720 * absa_var_190(ind0_var_1706 + 9, 1 : 14) + fac110_var_1721 * absa_var_190(ind0_var_1706 + 10, 1 : 14) + fac210_var_1722 * absa_var_190(ind0_var_1706 + 11, 1 : 14))
      ELSE IF (specparm_var_1740 .GT. 0.875D0) THEN
        tau_major_var_1736(1 : ng4) = speccomb_var_1703 * (fac200_var_1719 * absa_var_190(ind0_var_1706 - 1, 1 : 14) + fac100_var_1718 * absa_var_190(ind0_var_1706, 1 : 14) + fac000_var_1717 * absa_var_190(ind0_var_1706 + 1, 1 : 14) + fac210_var_1722 * absa_var_190(ind0_var_1706 + 8, 1 : 14) + fac110_var_1721 * absa_var_190(ind0_var_1706 + 9, 1 : 14) + fac010_var_1720 * absa_var_190(ind0_var_1706 + 10, 1 : 14))
      ELSE
        tau_major_var_1736(1 : ng4) = speccomb_var_1703 * (fac000_var_1717 * absa_var_190(ind0_var_1706, 1 : 14) + fac100_var_1718 * absa_var_190(ind0_var_1706 + 1, 1 : 14) + fac010_var_1720 * absa_var_190(ind0_var_1706 + 9, 1 : 14) + fac110_var_1721 * absa_var_190(ind0_var_1706 + 10, 1 : 14))
      END IF
      IF (specparm1_var_1743 .LT. 0.125D0) THEN
        tau_major1_var_1737(1 : ng4) = speccomb1_var_1704 * (fac001_var_1723 * absa_var_190(ind1_var_1707, 1 : 14) + fac101_var_1724 * absa_var_190(ind1_var_1707 + 1, 1 : 14) + fac201_var_1725 * absa_var_190(ind1_var_1707 + 2, 1 : 14) + fac011_var_1726 * absa_var_190(ind1_var_1707 + 9, 1 : 14) + fac111_var_1727 * absa_var_190(ind1_var_1707 + 10, 1 : 14) + fac211_var_1728 * absa_var_190(ind1_var_1707 + 11, 1 : 14))
      ELSE IF (specparm1_var_1743 .GT. 0.875D0) THEN
        tau_major1_var_1737(1 : ng4) = speccomb1_var_1704 * (fac201_var_1725 * absa_var_190(ind1_var_1707 - 1, 1 : 14) + fac101_var_1724 * absa_var_190(ind1_var_1707, 1 : 14) + fac001_var_1723 * absa_var_190(ind1_var_1707 + 1, 1 : 14) + fac211_var_1728 * absa_var_190(ind1_var_1707 + 8, 1 : 14) + fac111_var_1727 * absa_var_190(ind1_var_1707 + 9, 1 : 14) + fac011_var_1726 * absa_var_190(ind1_var_1707 + 10, 1 : 14))
      ELSE
        tau_major1_var_1737(1 : ng4) = speccomb1_var_1704 * (fac001_var_1723 * absa_var_190(ind1_var_1707, 1 : 14) + fac101_var_1724 * absa_var_190(ind1_var_1707 + 1, 1 : 14) + fac011_var_1726 * absa_var_190(ind1_var_1707 + 9, 1 : 14) + fac111_var_1727 * absa_var_190(ind1_var_1707 + 10, 1 : 14))
      END IF
      DO ig_var_1710 = 1, 14
        tauself_var_1735 = selffac_var_1692(jl_var_1757, lay_var_1712) * (selfref_var_192(inds_var_1708, ig_var_1710) + selffrac_var_1693(jl_var_1757, lay_var_1712) * (selfref_var_192(inds_var_1708 + 1, ig_var_1710) - selfref_var_192(inds_var_1708, ig_var_1710)))
        taufor_var_1734 = forfac_var_1701(jl_var_1757, lay_var_1712) * (forref_var_193(indf_var_1709, ig_var_1710) + forfrac_var_1702(jl_var_1757, lay_var_1712) * (forref_var_193(indf_var_1709 + 1, ig_var_1710) - forref_var_193(indf_var_1709, ig_var_1710)))
        taug_var_1678(jl_var_1757, 38 + ig_var_1710, lay_var_1712) = tau_major_var_1736(ig_var_1710) + tau_major1_var_1737(ig_var_1710) + tauself_var_1735 + taufor_var_1734
        fracs_var_1695(jl_var_1757, 38 + ig_var_1710, lay_var_1712) = fracrefa_var_188(ig_var_1710, jpl_var_1714) + fpl_var_1744 * (fracrefa_var_188(ig_var_1710, jpl_var_1714 + 1) - fracrefa_var_188(ig_var_1710, jpl_var_1714))
      END DO
    END DO
  END DO
  DO lay_var_1712 = laytrop_max_var_1748 + 1, klev_var_1677
    DO jl_var_1757 = kidia_var_1675, kfdia_var_1676
      speccomb_var_1703 = colo3_var_1690(jl_var_1757, lay_var_1712) + rat_o3co2_var_1698(jl_var_1757, lay_var_1712) * colco2_var_1689(jl_var_1757, lay_var_1712)
      specparm_var_1740 = MIN(colo3_var_1690(jl_var_1757, lay_var_1712) / speccomb_var_1703, oneminus_var_1687)
      specmult_var_1739 = 4.0D0 * (specparm_var_1740)
      js_var_1711 = 1 + INT(specmult_var_1739)
      fs_var_1738 = ((specmult_var_1739) - AINT((specmult_var_1739)))
      speccomb1_var_1704 = colo3_var_1690(jl_var_1757, lay_var_1712) + rat_o3co2_1_var_1699(jl_var_1757, lay_var_1712) * colco2_var_1689(jl_var_1757, lay_var_1712)
      specparm1_var_1743 = MIN(colo3_var_1690(jl_var_1757, lay_var_1712) / speccomb1_var_1704, oneminus_var_1687)
      specmult1_var_1742 = 4.0D0 * (specparm1_var_1743)
      js1_var_1713 = 1 + INT(specmult1_var_1742)
      fs1_var_1741 = ((specmult1_var_1742) - AINT((specmult1_var_1742)))
      fac000_var_1717 = (1.0D0 - fs_var_1738) * fac00_var_1680(jl_var_1757, lay_var_1712)
      fac010_var_1720 = (1.0D0 - fs_var_1738) * fac10_var_1682(jl_var_1757, lay_var_1712)
      fac100_var_1718 = fs_var_1738 * fac00_var_1680(jl_var_1757, lay_var_1712)
      fac110_var_1721 = fs_var_1738 * fac10_var_1682(jl_var_1757, lay_var_1712)
      fac001_var_1723 = (1.0D0 - fs1_var_1741) * fac01_var_1681(jl_var_1757, lay_var_1712)
      fac011_var_1726 = (1.0D0 - fs1_var_1741) * fac11_var_1683(jl_var_1757, lay_var_1712)
      fac101_var_1724 = fs1_var_1741 * fac01_var_1681(jl_var_1757, lay_var_1712)
      fac111_var_1727 = fs1_var_1741 * fac11_var_1683(jl_var_1757, lay_var_1712)
      speccomb_planck_var_1705 = colo3_var_1690(jl_var_1757, lay_var_1712) + refrat_planck_b_var_1716 * colco2_var_1689(jl_var_1757, lay_var_1712)
      specparm_planck_var_1746 = MIN(colo3_var_1690(jl_var_1757, lay_var_1712) / speccomb_planck_var_1705, oneminus_var_1687)
      specmult_planck_var_1745 = 4.0D0 * specparm_planck_var_1746
      jpl_var_1714 = 1 + INT(specmult_planck_var_1745)
      fpl_var_1744 = ((specmult_planck_var_1745) - AINT((specmult_planck_var_1745)))
      ind0_var_1706 = ((jp_var_1684(jl_var_1757, lay_var_1712) - 13) * 5 + (jt_var_1685(jl_var_1757, lay_var_1712) - 1)) * nspb_var_238(4) + js_var_1711
      ind1_var_1707 = ((jp_var_1684(jl_var_1757, lay_var_1712) - 12) * 5 + (jt1_var_1686(jl_var_1757, lay_var_1712) - 1)) * nspb_var_238(4) + js1_var_1713
      DO ig_var_1710 = 1, 14
        taug_var_1678(jl_var_1757, 38 + ig_var_1710, lay_var_1712) = speccomb_var_1703 * (fac000_var_1717 * absb_var_191(ind0_var_1706, ig_var_1710) + fac100_var_1718 * absb_var_191(ind0_var_1706 + 1, ig_var_1710) + fac010_var_1720 * absb_var_191(ind0_var_1706 + 5, ig_var_1710) + fac110_var_1721 * absb_var_191(ind0_var_1706 + 6, ig_var_1710)) + speccomb1_var_1704 * (fac001_var_1723 * absb_var_191(ind1_var_1707, ig_var_1710) + fac101_var_1724 * absb_var_191(ind1_var_1707 + 1, ig_var_1710) + fac011_var_1726 * absb_var_191(ind1_var_1707 + 5, ig_var_1710) + fac111_var_1727 * absb_var_191(ind1_var_1707 + 6, ig_var_1710))
        fracs_var_1695(jl_var_1757, 38 + ig_var_1710, lay_var_1712) = fracrefb_var_189(ig_var_1710, jpl_var_1714) + fpl_var_1744 * (fracrefb_var_189(ig_var_1710, jpl_var_1714 + 1) - fracrefb_var_189(ig_var_1710, jpl_var_1714))
      END DO
    END DO
  END DO
  DO lay_var_1712 = laytrop_max_var_1748 + 1, klev_var_1677
    DO jl_var_1757 = kidia_var_1675, kfdia_var_1676
      taug_var_1678(jl_var_1757, 46, lay_var_1712) = taug_var_1678(jl_var_1757, 46, lay_var_1712) * 0.92D0
      taug_var_1678(jl_var_1757, 47, lay_var_1712) = taug_var_1678(jl_var_1757, 47, lay_var_1712) * 0.88D0
      taug_var_1678(jl_var_1757, 48, lay_var_1712) = taug_var_1678(jl_var_1757, 48, lay_var_1712) * 1.07D0
      taug_var_1678(jl_var_1757, 49, lay_var_1712) = taug_var_1678(jl_var_1757, 49, lay_var_1712) * 1.1D0
      taug_var_1678(jl_var_1757, 50, lay_var_1712) = taug_var_1678(jl_var_1757, 50, lay_var_1712) * 0.99D0
      taug_var_1678(jl_var_1757, 51, lay_var_1712) = taug_var_1678(jl_var_1757, 51, lay_var_1712) * 0.88D0
      taug_var_1678(jl_var_1757, 52, lay_var_1712) = taug_var_1678(jl_var_1757, 52, lay_var_1712) * 0.943D0
    END DO
  END DO
  IF (laytrop_max_var_1748 /= laytrop_min_var_1747) THEN
    DO lay_var_1712 = laytrop_min_var_1747 + 1, laytrop_max_var_1748
      ixc0_var_1754 = ixc_var_1749(lay_var_1712)
      DO ixp_var_1755 = 1, ixc0_var_1754
        jl_var_1757 = ixlow_var_1750(ixp_var_1755, lay_var_1712)
        speccomb_var_1703 = colh2o_var_1688(jl_var_1757, lay_var_1712) + rat_h2oco2_var_1696(jl_var_1757, lay_var_1712) * colco2_var_1689(jl_var_1757, lay_var_1712)
        specparm_var_1740 = MIN(colh2o_var_1688(jl_var_1757, lay_var_1712) / speccomb_var_1703, oneminus_var_1687)
        specmult_var_1739 = 8.0D0 * (specparm_var_1740)
        js_var_1711 = 1 + INT(specmult_var_1739)
        fs_var_1738 = ((specmult_var_1739) - AINT((specmult_var_1739)))
        speccomb1_var_1704 = colh2o_var_1688(jl_var_1757, lay_var_1712) + rat_h2oco2_1_var_1697(jl_var_1757, lay_var_1712) * colco2_var_1689(jl_var_1757, lay_var_1712)
        specparm1_var_1743 = MIN(colh2o_var_1688(jl_var_1757, lay_var_1712) / speccomb1_var_1704, oneminus_var_1687)
        specmult1_var_1742 = 8.0D0 * (specparm1_var_1743)
        js1_var_1713 = 1 + INT(specmult1_var_1742)
        fs1_var_1741 = ((specmult1_var_1742) - AINT((specmult1_var_1742)))
        speccomb_planck_var_1705 = colh2o_var_1688(jl_var_1757, lay_var_1712) + refrat_planck_a_var_1715 * colco2_var_1689(jl_var_1757, lay_var_1712)
        specparm_planck_var_1746 = MIN(colh2o_var_1688(jl_var_1757, lay_var_1712) / speccomb_planck_var_1705, oneminus_var_1687)
        specmult_planck_var_1745 = 8.0D0 * specparm_planck_var_1746
        jpl_var_1714 = 1 + INT(specmult_planck_var_1745)
        fpl_var_1744 = ((specmult_planck_var_1745) - AINT((specmult_planck_var_1745)))
        ind0_var_1706 = ((jp_var_1684(jl_var_1757, lay_var_1712) - 1) * 5 + (jt_var_1685(jl_var_1757, lay_var_1712) - 1)) * nspa_var_237(4) + js_var_1711
        ind1_var_1707 = (jp_var_1684(jl_var_1757, lay_var_1712) * 5 + (jt1_var_1686(jl_var_1757, lay_var_1712) - 1)) * nspa_var_237(4) + js1_var_1713
        inds_var_1708 = indself_var_1694(jl_var_1757, lay_var_1712)
        indf_var_1709 = indfor_var_1700(jl_var_1757, lay_var_1712)
        IF (specparm_var_1740 .LT. 0.125D0) THEN
          p_var_1729 = fs_var_1738 - 1.0D0
          p4_var_1730 = p_var_1729 ** 4
          fk0_var_1731 = p4_var_1730
          fk1_var_1732 = 1.0D0 - p_var_1729 - 2.0D0 * p4_var_1730
          fk2_var_1733 = p_var_1729 + p4_var_1730
          fac000_var_1717 = fk0_var_1731 * fac00_var_1680(jl_var_1757, lay_var_1712)
          fac100_var_1718 = fk1_var_1732 * fac00_var_1680(jl_var_1757, lay_var_1712)
          fac200_var_1719 = fk2_var_1733 * fac00_var_1680(jl_var_1757, lay_var_1712)
          fac010_var_1720 = fk0_var_1731 * fac10_var_1682(jl_var_1757, lay_var_1712)
          fac110_var_1721 = fk1_var_1732 * fac10_var_1682(jl_var_1757, lay_var_1712)
          fac210_var_1722 = fk2_var_1733 * fac10_var_1682(jl_var_1757, lay_var_1712)
        ELSE IF (specparm_var_1740 .GT. 0.875D0) THEN
          p_var_1729 = - fs_var_1738
          p4_var_1730 = p_var_1729 ** 4
          fk0_var_1731 = p4_var_1730
          fk1_var_1732 = 1.0D0 - p_var_1729 - 2.0D0 * p4_var_1730
          fk2_var_1733 = p_var_1729 + p4_var_1730
          fac000_var_1717 = fk0_var_1731 * fac00_var_1680(jl_var_1757, lay_var_1712)
          fac100_var_1718 = fk1_var_1732 * fac00_var_1680(jl_var_1757, lay_var_1712)
          fac200_var_1719 = fk2_var_1733 * fac00_var_1680(jl_var_1757, lay_var_1712)
          fac010_var_1720 = fk0_var_1731 * fac10_var_1682(jl_var_1757, lay_var_1712)
          fac110_var_1721 = fk1_var_1732 * fac10_var_1682(jl_var_1757, lay_var_1712)
          fac210_var_1722 = fk2_var_1733 * fac10_var_1682(jl_var_1757, lay_var_1712)
        ELSE
          fac000_var_1717 = (1.0D0 - fs_var_1738) * fac00_var_1680(jl_var_1757, lay_var_1712)
          fac010_var_1720 = (1.0D0 - fs_var_1738) * fac10_var_1682(jl_var_1757, lay_var_1712)
          fac100_var_1718 = fs_var_1738 * fac00_var_1680(jl_var_1757, lay_var_1712)
          fac110_var_1721 = fs_var_1738 * fac10_var_1682(jl_var_1757, lay_var_1712)
          fac200_var_1719 = 0.0D0
          fac210_var_1722 = 0.0D0
        END IF
        IF (specparm1_var_1743 .LT. 0.125D0) THEN
          p_var_1729 = fs1_var_1741 - 1.0D0
          p4_var_1730 = p_var_1729 ** 4
          fk0_var_1731 = p4_var_1730
          fk1_var_1732 = 1.0D0 - p_var_1729 - 2.0D0 * p4_var_1730
          fk2_var_1733 = p_var_1729 + p4_var_1730
          fac001_var_1723 = fk0_var_1731 * fac01_var_1681(jl_var_1757, lay_var_1712)
          fac101_var_1724 = fk1_var_1732 * fac01_var_1681(jl_var_1757, lay_var_1712)
          fac201_var_1725 = fk2_var_1733 * fac01_var_1681(jl_var_1757, lay_var_1712)
          fac011_var_1726 = fk0_var_1731 * fac11_var_1683(jl_var_1757, lay_var_1712)
          fac111_var_1727 = fk1_var_1732 * fac11_var_1683(jl_var_1757, lay_var_1712)
          fac211_var_1728 = fk2_var_1733 * fac11_var_1683(jl_var_1757, lay_var_1712)
        ELSE IF (specparm1_var_1743 .GT. 0.875D0) THEN
          p_var_1729 = - fs1_var_1741
          p4_var_1730 = p_var_1729 ** 4
          fk0_var_1731 = p4_var_1730
          fk1_var_1732 = 1.0D0 - p_var_1729 - 2.0D0 * p4_var_1730
          fk2_var_1733 = p_var_1729 + p4_var_1730
          fac001_var_1723 = fk0_var_1731 * fac01_var_1681(jl_var_1757, lay_var_1712)
          fac101_var_1724 = fk1_var_1732 * fac01_var_1681(jl_var_1757, lay_var_1712)
          fac201_var_1725 = fk2_var_1733 * fac01_var_1681(jl_var_1757, lay_var_1712)
          fac011_var_1726 = fk0_var_1731 * fac11_var_1683(jl_var_1757, lay_var_1712)
          fac111_var_1727 = fk1_var_1732 * fac11_var_1683(jl_var_1757, lay_var_1712)
          fac211_var_1728 = fk2_var_1733 * fac11_var_1683(jl_var_1757, lay_var_1712)
        ELSE
          fac001_var_1723 = (1.0D0 - fs1_var_1741) * fac01_var_1681(jl_var_1757, lay_var_1712)
          fac011_var_1726 = (1.0D0 - fs1_var_1741) * fac11_var_1683(jl_var_1757, lay_var_1712)
          fac101_var_1724 = fs1_var_1741 * fac01_var_1681(jl_var_1757, lay_var_1712)
          fac111_var_1727 = fs1_var_1741 * fac11_var_1683(jl_var_1757, lay_var_1712)
          fac201_var_1725 = 0.0D0
          fac211_var_1728 = 0.0D0
        END IF
        IF (specparm_var_1740 .LT. 0.125D0) THEN
          tau_major_var_1736(1 : ng4) = speccomb_var_1703 * (fac000_var_1717 * absa_var_190(ind0_var_1706, 1 : 14) + fac100_var_1718 * absa_var_190(ind0_var_1706 + 1, 1 : 14) + fac200_var_1719 * absa_var_190(ind0_var_1706 + 2, 1 : 14) + fac010_var_1720 * absa_var_190(ind0_var_1706 + 9, 1 : 14) + fac110_var_1721 * absa_var_190(ind0_var_1706 + 10, 1 : 14) + fac210_var_1722 * absa_var_190(ind0_var_1706 + 11, 1 : 14))
        ELSE IF (specparm_var_1740 .GT. 0.875D0) THEN
          tau_major_var_1736(1 : ng4) = speccomb_var_1703 * (fac200_var_1719 * absa_var_190(ind0_var_1706 - 1, 1 : 14) + fac100_var_1718 * absa_var_190(ind0_var_1706, 1 : 14) + fac000_var_1717 * absa_var_190(ind0_var_1706 + 1, 1 : 14) + fac210_var_1722 * absa_var_190(ind0_var_1706 + 8, 1 : 14) + fac110_var_1721 * absa_var_190(ind0_var_1706 + 9, 1 : 14) + fac010_var_1720 * absa_var_190(ind0_var_1706 + 10, 1 : 14))
        ELSE
          tau_major_var_1736(1 : ng4) = speccomb_var_1703 * (fac000_var_1717 * absa_var_190(ind0_var_1706, 1 : 14) + fac100_var_1718 * absa_var_190(ind0_var_1706 + 1, 1 : 14) + fac010_var_1720 * absa_var_190(ind0_var_1706 + 9, 1 : 14) + fac110_var_1721 * absa_var_190(ind0_var_1706 + 10, 1 : 14))
        END IF
        IF (specparm1_var_1743 .LT. 0.125D0) THEN
          tau_major1_var_1737(1 : ng4) = speccomb1_var_1704 * (fac001_var_1723 * absa_var_190(ind1_var_1707, 1 : 14) + fac101_var_1724 * absa_var_190(ind1_var_1707 + 1, 1 : 14) + fac201_var_1725 * absa_var_190(ind1_var_1707 + 2, 1 : 14) + fac011_var_1726 * absa_var_190(ind1_var_1707 + 9, 1 : 14) + fac111_var_1727 * absa_var_190(ind1_var_1707 + 10, 1 : 14) + fac211_var_1728 * absa_var_190(ind1_var_1707 + 11, 1 : 14))
        ELSE IF (specparm1_var_1743 .GT. 0.875D0) THEN
          tau_major1_var_1737(1 : ng4) = speccomb1_var_1704 * (fac201_var_1725 * absa_var_190(ind1_var_1707 - 1, 1 : 14) + fac101_var_1724 * absa_var_190(ind1_var_1707, 1 : 14) + fac001_var_1723 * absa_var_190(ind1_var_1707 + 1, 1 : 14) + fac211_var_1728 * absa_var_190(ind1_var_1707 + 8, 1 : 14) + fac111_var_1727 * absa_var_190(ind1_var_1707 + 9, 1 : 14) + fac011_var_1726 * absa_var_190(ind1_var_1707 + 10, 1 : 14))
        ELSE
          tau_major1_var_1737(1 : ng4) = speccomb1_var_1704 * (fac001_var_1723 * absa_var_190(ind1_var_1707, 1 : 14) + fac101_var_1724 * absa_var_190(ind1_var_1707 + 1, 1 : 14) + fac011_var_1726 * absa_var_190(ind1_var_1707 + 9, 1 : 14) + fac111_var_1727 * absa_var_190(ind1_var_1707 + 10, 1 : 14))
        END IF
        DO ig_var_1710 = 1, 14
          tauself_var_1735 = selffac_var_1692(jl_var_1757, lay_var_1712) * (selfref_var_192(inds_var_1708, ig_var_1710) + selffrac_var_1693(jl_var_1757, lay_var_1712) * (selfref_var_192(inds_var_1708 + 1, ig_var_1710) - selfref_var_192(inds_var_1708, ig_var_1710)))
          taufor_var_1734 = forfac_var_1701(jl_var_1757, lay_var_1712) * (forref_var_193(indf_var_1709, ig_var_1710) + forfrac_var_1702(jl_var_1757, lay_var_1712) * (forref_var_193(indf_var_1709 + 1, ig_var_1710) - forref_var_193(indf_var_1709, ig_var_1710)))
          taug_var_1678(jl_var_1757, 38 + ig_var_1710, lay_var_1712) = tau_major_var_1736(ig_var_1710) + tau_major1_var_1737(ig_var_1710) + tauself_var_1735 + taufor_var_1734
          fracs_var_1695(jl_var_1757, 38 + ig_var_1710, lay_var_1712) = fracrefa_var_188(ig_var_1710, jpl_var_1714) + fpl_var_1744 * (fracrefa_var_188(ig_var_1710, jpl_var_1714 + 1) - fracrefa_var_188(ig_var_1710, jpl_var_1714))
        END DO
      END DO
      ixc0_var_1754 = kfdia_var_1676 - kidia_var_1675 + 1 - ixc0_var_1754
      DO ixp_var_1755 = 1, ixc0_var_1754
        jl_var_1757 = ixhigh_var_1751(ixp_var_1755, lay_var_1712)
        speccomb_var_1703 = colo3_var_1690(jl_var_1757, lay_var_1712) + rat_o3co2_var_1698(jl_var_1757, lay_var_1712) * colco2_var_1689(jl_var_1757, lay_var_1712)
        specparm_var_1740 = MIN(colo3_var_1690(jl_var_1757, lay_var_1712) / speccomb_var_1703, oneminus_var_1687)
        specmult_var_1739 = 4.0D0 * (specparm_var_1740)
        js_var_1711 = 1 + INT(specmult_var_1739)
        fs_var_1738 = ((specmult_var_1739) - AINT((specmult_var_1739)))
        speccomb1_var_1704 = colo3_var_1690(jl_var_1757, lay_var_1712) + rat_o3co2_1_var_1699(jl_var_1757, lay_var_1712) * colco2_var_1689(jl_var_1757, lay_var_1712)
        specparm1_var_1743 = MIN(colo3_var_1690(jl_var_1757, lay_var_1712) / speccomb1_var_1704, oneminus_var_1687)
        specmult1_var_1742 = 4.0D0 * (specparm1_var_1743)
        js1_var_1713 = 1 + INT(specmult1_var_1742)
        fs1_var_1741 = ((specmult1_var_1742) - AINT((specmult1_var_1742)))
        fac000_var_1717 = (1.0D0 - fs_var_1738) * fac00_var_1680(jl_var_1757, lay_var_1712)
        fac010_var_1720 = (1.0D0 - fs_var_1738) * fac10_var_1682(jl_var_1757, lay_var_1712)
        fac100_var_1718 = fs_var_1738 * fac00_var_1680(jl_var_1757, lay_var_1712)
        fac110_var_1721 = fs_var_1738 * fac10_var_1682(jl_var_1757, lay_var_1712)
        fac001_var_1723 = (1.0D0 - fs1_var_1741) * fac01_var_1681(jl_var_1757, lay_var_1712)
        fac011_var_1726 = (1.0D0 - fs1_var_1741) * fac11_var_1683(jl_var_1757, lay_var_1712)
        fac101_var_1724 = fs1_var_1741 * fac01_var_1681(jl_var_1757, lay_var_1712)
        fac111_var_1727 = fs1_var_1741 * fac11_var_1683(jl_var_1757, lay_var_1712)
        speccomb_planck_var_1705 = colo3_var_1690(jl_var_1757, lay_var_1712) + refrat_planck_b_var_1716 * colco2_var_1689(jl_var_1757, lay_var_1712)
        specparm_planck_var_1746 = MIN(colo3_var_1690(jl_var_1757, lay_var_1712) / speccomb_planck_var_1705, oneminus_var_1687)
        specmult_planck_var_1745 = 4.0D0 * specparm_planck_var_1746
        jpl_var_1714 = 1 + INT(specmult_planck_var_1745)
        fpl_var_1744 = ((specmult_planck_var_1745) - AINT((specmult_planck_var_1745)))
        ind0_var_1706 = ((jp_var_1684(jl_var_1757, lay_var_1712) - 13) * 5 + (jt_var_1685(jl_var_1757, lay_var_1712) - 1)) * nspb_var_238(4) + js_var_1711
        ind1_var_1707 = ((jp_var_1684(jl_var_1757, lay_var_1712) - 12) * 5 + (jt1_var_1686(jl_var_1757, lay_var_1712) - 1)) * nspb_var_238(4) + js1_var_1713
        DO ig_var_1710 = 1, 14
          taug_var_1678(jl_var_1757, 38 + ig_var_1710, lay_var_1712) = speccomb_var_1703 * (fac000_var_1717 * absb_var_191(ind0_var_1706, ig_var_1710) + fac100_var_1718 * absb_var_191(ind0_var_1706 + 1, ig_var_1710) + fac010_var_1720 * absb_var_191(ind0_var_1706 + 5, ig_var_1710) + fac110_var_1721 * absb_var_191(ind0_var_1706 + 6, ig_var_1710)) + speccomb1_var_1704 * (fac001_var_1723 * absb_var_191(ind1_var_1707, ig_var_1710) + fac101_var_1724 * absb_var_191(ind1_var_1707 + 1, ig_var_1710) + fac011_var_1726 * absb_var_191(ind1_var_1707 + 5, ig_var_1710) + fac111_var_1727 * absb_var_191(ind1_var_1707 + 6, ig_var_1710))
          fracs_var_1695(jl_var_1757, 38 + ig_var_1710, lay_var_1712) = fracrefb_var_189(ig_var_1710, jpl_var_1714) + fpl_var_1744 * (fracrefb_var_189(ig_var_1710, jpl_var_1714 + 1) - fracrefb_var_189(ig_var_1710, jpl_var_1714))
        END DO
      END DO
      DO ixp_var_1755 = 1, ixc0_var_1754
        jl_var_1757 = ixhigh_var_1751(ixp_var_1755, lay_var_1712)
        taug_var_1678(jl_var_1757, 46, lay_var_1712) = taug_var_1678(jl_var_1757, 46, lay_var_1712) * 0.92D0
        taug_var_1678(jl_var_1757, 47, lay_var_1712) = taug_var_1678(jl_var_1757, 47, lay_var_1712) * 0.88D0
        taug_var_1678(jl_var_1757, 48, lay_var_1712) = taug_var_1678(jl_var_1757, 48, lay_var_1712) * 1.07D0
        taug_var_1678(jl_var_1757, 49, lay_var_1712) = taug_var_1678(jl_var_1757, 49, lay_var_1712) * 1.1D0
        taug_var_1678(jl_var_1757, 50, lay_var_1712) = taug_var_1678(jl_var_1757, 50, lay_var_1712) * 0.99D0
        taug_var_1678(jl_var_1757, 51, lay_var_1712) = taug_var_1678(jl_var_1757, 51, lay_var_1712) * 0.88D0
        taug_var_1678(jl_var_1757, 52, lay_var_1712) = taug_var_1678(jl_var_1757, 52, lay_var_1712) * 0.943D0
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol4
SUBROUTINE rrtm_taumol6(kidia_var_1758, kfdia_var_1759, klev_var_1760, taug_var_1761, wx_var_1762, p_tauaerl_var_1763, fac00_var_1764, fac01_var_1765, fac10_var_1766, fac11_var_1767, forfac_var_1780, forfrac_var_1781, indfor_var_1779, jp_var_1768, jt_var_1769, jt1_var_1770, colh2o_var_1771, colco2_var_1772, coldry_var_1773, laytrop_var_1774, selffac_var_1775, selffrac_var_1776, indself_var_1777, fracs_var_1778, minorfrac_var_1782, indminor_var_1783)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_237
  USE yoerrta6, ONLY: absa_var_203, cfc11adj, cfc12_var_202, forref_var_206, fracrefa_var_201, ka_mco2_var_205, selfref_var_204
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1758
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1759
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1760
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1761(kidia_var_1758 : kfdia_var_1759, 140, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: wx_var_1762(kidia_var_1758 : kfdia_var_1759, 4, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1763(kidia_var_1758 : kfdia_var_1759, klev_var_1760, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1764(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1765(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1766(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1767(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1768(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1769(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1770(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1771(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1772(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1773(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1774(kidia_var_1758 : kfdia_var_1759)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1775(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1776(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1777(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1778(kidia_var_1758 : kfdia_var_1759, 140, klev_var_1760)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1779(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1780(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1781(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1782(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1783(kidia_var_1758 : kfdia_var_1759, klev_var_1760)
  INTEGER(KIND = 4) :: ind0_var_1784, ind1_var_1785, inds_var_1786, indf_var_1787, indm_var_1788
  INTEGER(KIND = 4) :: ig_var_1789, lay_var_1790
  REAL(KIND = 8) :: adjfac_var_1791, adjcolco2_var_1792, ratco2_var_1793, chi_co2_var_1794
  REAL(KIND = 8) :: taufor_var_1795, tauself_var_1796, absco2_var_1797
  INTEGER(KIND = 4) :: laytrop_min_var_1798, laytrop_max_var_1799
  INTEGER(KIND = 4) :: ixc_var_1800(klev_var_1760), ixlow_var_1801(kfdia_var_1759, klev_var_1760), ixhigh_var_1802(kfdia_var_1759, klev_var_1760)
  INTEGER(KIND = 4) :: ich_var_1803, icl_var_1804, ixc0_var_1805, ixp_var_1806, jc_var_1807, jl_var_1808
  laytrop_min_var_1798 = MINVAL(laytrop_var_1774)
  laytrop_max_var_1799 = MAXVAL(laytrop_var_1774)
  ixlow_var_1801 = 0
  ixhigh_var_1802 = 0
  ixc_var_1800 = 0
  DO lay_var_1790 = laytrop_min_var_1798 + 1, laytrop_max_var_1799
    icl_var_1804 = 0
    ich_var_1803 = 0
    DO jc_var_1807 = kidia_var_1758, kfdia_var_1759
      IF (lay_var_1790 <= laytrop_var_1774(jc_var_1807)) THEN
        icl_var_1804 = icl_var_1804 + 1
        ixlow_var_1801(icl_var_1804, lay_var_1790) = jc_var_1807
      ELSE
        ich_var_1803 = ich_var_1803 + 1
        ixhigh_var_1802(ich_var_1803, lay_var_1790) = jc_var_1807
      END IF
    END DO
    ixc_var_1800(lay_var_1790) = icl_var_1804
  END DO
  DO lay_var_1790 = 1, laytrop_min_var_1798
    DO jl_var_1808 = kidia_var_1758, kfdia_var_1759
      chi_co2_var_1794 = colco2_var_1772(jl_var_1808, lay_var_1790) / (coldry_var_1773(jl_var_1808, lay_var_1790))
      ratco2_var_1793 = 1D+20 * chi_co2_var_1794 / chi_mls(2, jp_var_1768(jl_var_1808, lay_var_1790) + 1)
      IF (ratco2_var_1793 .GT. 3.0D0) THEN
        adjfac_var_1791 = 2.0D0 + (ratco2_var_1793 - 2.0D0) ** 0.77D0
        adjcolco2_var_1792 = adjfac_var_1791 * chi_mls(2, jp_var_1768(jl_var_1808, lay_var_1790) + 1) * coldry_var_1773(jl_var_1808, lay_var_1790) * 1D-20
      ELSE
        adjcolco2_var_1792 = colco2_var_1772(jl_var_1808, lay_var_1790)
      END IF
      ind0_var_1784 = ((jp_var_1768(jl_var_1808, lay_var_1790) - 1) * 5 + (jt_var_1769(jl_var_1808, lay_var_1790) - 1)) * nspa_var_237(6) + 1
      ind1_var_1785 = (jp_var_1768(jl_var_1808, lay_var_1790) * 5 + (jt1_var_1770(jl_var_1808, lay_var_1790) - 1)) * nspa_var_237(6) + 1
      inds_var_1786 = indself_var_1777(jl_var_1808, lay_var_1790)
      indf_var_1787 = indfor_var_1779(jl_var_1808, lay_var_1790)
      indm_var_1788 = indminor_var_1783(jl_var_1808, lay_var_1790)
      DO ig_var_1789 = 1, 8
        tauself_var_1796 = selffac_var_1775(jl_var_1808, lay_var_1790) * (selfref_var_204(inds_var_1786, ig_var_1789) + selffrac_var_1776(jl_var_1808, lay_var_1790) * (selfref_var_204(inds_var_1786 + 1, ig_var_1789) - selfref_var_204(inds_var_1786, ig_var_1789)))
        taufor_var_1795 = forfac_var_1780(jl_var_1808, lay_var_1790) * (forref_var_206(indf_var_1787, ig_var_1789) + forfrac_var_1781(jl_var_1808, lay_var_1790) * (forref_var_206(indf_var_1787 + 1, ig_var_1789) - forref_var_206(indf_var_1787, ig_var_1789)))
        absco2_var_1797 = (ka_mco2_var_205(indm_var_1788, ig_var_1789) + minorfrac_var_1782(jl_var_1808, lay_var_1790) * (ka_mco2_var_205(indm_var_1788 + 1, ig_var_1789) - ka_mco2_var_205(indm_var_1788, ig_var_1789)))
        taug_var_1761(jl_var_1808, 68 + ig_var_1789, lay_var_1790) = colh2o_var_1771(jl_var_1808, lay_var_1790) * (fac00_var_1764(jl_var_1808, lay_var_1790) * absa_var_203(ind0_var_1784, ig_var_1789) + fac10_var_1766(jl_var_1808, lay_var_1790) * absa_var_203(ind0_var_1784 + 1, ig_var_1789) + fac01_var_1765(jl_var_1808, lay_var_1790) * absa_var_203(ind1_var_1785, ig_var_1789) + fac11_var_1767(jl_var_1808, lay_var_1790) * absa_var_203(ind1_var_1785 + 1, ig_var_1789)) + tauself_var_1796 + taufor_var_1795 + adjcolco2_var_1792 * absco2_var_1797 + wx_var_1762(jl_var_1808, 2, lay_var_1790) * cfc11adj(ig_var_1789) + wx_var_1762(jl_var_1808, 3, lay_var_1790) * cfc12_var_202(ig_var_1789)
        fracs_var_1778(jl_var_1808, 68 + ig_var_1789, lay_var_1790) = fracrefa_var_201(ig_var_1789)
      END DO
    END DO
  END DO
  DO ig_var_1789 = 1, 8
    DO lay_var_1790 = laytrop_max_var_1799 + 1, klev_var_1760
      DO jl_var_1808 = kidia_var_1758, kfdia_var_1759
        taug_var_1761(jl_var_1808, 68 + ig_var_1789, lay_var_1790) = 0.0D0 + wx_var_1762(jl_var_1808, 2, lay_var_1790) * cfc11adj(ig_var_1789) + wx_var_1762(jl_var_1808, 3, lay_var_1790) * cfc12_var_202(ig_var_1789)
        fracs_var_1778(jl_var_1808, 68 + ig_var_1789, lay_var_1790) = fracrefa_var_201(ig_var_1789)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1799 /= laytrop_min_var_1798) THEN
    DO lay_var_1790 = laytrop_min_var_1798 + 1, laytrop_max_var_1799
      ixc0_var_1805 = ixc_var_1800(lay_var_1790)
      DO ixp_var_1806 = 1, ixc0_var_1805
        jl_var_1808 = ixlow_var_1801(ixp_var_1806, lay_var_1790)
        chi_co2_var_1794 = colco2_var_1772(jl_var_1808, lay_var_1790) / (coldry_var_1773(jl_var_1808, lay_var_1790))
        ratco2_var_1793 = 1D+20 * chi_co2_var_1794 / chi_mls(2, jp_var_1768(jl_var_1808, lay_var_1790) + 1)
        IF (ratco2_var_1793 .GT. 3.0D0) THEN
          adjfac_var_1791 = 2.0D0 + (ratco2_var_1793 - 2.0D0) ** 0.77D0
          adjcolco2_var_1792 = adjfac_var_1791 * chi_mls(2, jp_var_1768(jl_var_1808, lay_var_1790) + 1) * coldry_var_1773(jl_var_1808, lay_var_1790) * 1D-20
        ELSE
          adjcolco2_var_1792 = colco2_var_1772(jl_var_1808, lay_var_1790)
        END IF
        ind0_var_1784 = ((jp_var_1768(jl_var_1808, lay_var_1790) - 1) * 5 + (jt_var_1769(jl_var_1808, lay_var_1790) - 1)) * nspa_var_237(6) + 1
        ind1_var_1785 = (jp_var_1768(jl_var_1808, lay_var_1790) * 5 + (jt1_var_1770(jl_var_1808, lay_var_1790) - 1)) * nspa_var_237(6) + 1
        inds_var_1786 = indself_var_1777(jl_var_1808, lay_var_1790)
        indf_var_1787 = indfor_var_1779(jl_var_1808, lay_var_1790)
        indm_var_1788 = indminor_var_1783(jl_var_1808, lay_var_1790)
        DO ig_var_1789 = 1, 8
          tauself_var_1796 = selffac_var_1775(jl_var_1808, lay_var_1790) * (selfref_var_204(inds_var_1786, ig_var_1789) + selffrac_var_1776(jl_var_1808, lay_var_1790) * (selfref_var_204(inds_var_1786 + 1, ig_var_1789) - selfref_var_204(inds_var_1786, ig_var_1789)))
          taufor_var_1795 = forfac_var_1780(jl_var_1808, lay_var_1790) * (forref_var_206(indf_var_1787, ig_var_1789) + forfrac_var_1781(jl_var_1808, lay_var_1790) * (forref_var_206(indf_var_1787 + 1, ig_var_1789) - forref_var_206(indf_var_1787, ig_var_1789)))
          absco2_var_1797 = (ka_mco2_var_205(indm_var_1788, ig_var_1789) + minorfrac_var_1782(jl_var_1808, lay_var_1790) * (ka_mco2_var_205(indm_var_1788 + 1, ig_var_1789) - ka_mco2_var_205(indm_var_1788, ig_var_1789)))
          taug_var_1761(jl_var_1808, 68 + ig_var_1789, lay_var_1790) = colh2o_var_1771(jl_var_1808, lay_var_1790) * (fac00_var_1764(jl_var_1808, lay_var_1790) * absa_var_203(ind0_var_1784, ig_var_1789) + fac10_var_1766(jl_var_1808, lay_var_1790) * absa_var_203(ind0_var_1784 + 1, ig_var_1789) + fac01_var_1765(jl_var_1808, lay_var_1790) * absa_var_203(ind1_var_1785, ig_var_1789) + fac11_var_1767(jl_var_1808, lay_var_1790) * absa_var_203(ind1_var_1785 + 1, ig_var_1789)) + tauself_var_1796 + taufor_var_1795 + adjcolco2_var_1792 * absco2_var_1797 + wx_var_1762(jl_var_1808, 2, lay_var_1790) * cfc11adj(ig_var_1789) + wx_var_1762(jl_var_1808, 3, lay_var_1790) * cfc12_var_202(ig_var_1789)
          fracs_var_1778(jl_var_1808, 68 + ig_var_1789, lay_var_1790) = fracrefa_var_201(ig_var_1789)
        END DO
      END DO
      ixc0_var_1805 = kfdia_var_1759 - kidia_var_1758 + 1 - ixc0_var_1805
      DO ig_var_1789 = 1, 8
        DO ixp_var_1806 = 1, ixc0_var_1805
          jl_var_1808 = ixhigh_var_1802(ixp_var_1806, lay_var_1790)
          taug_var_1761(jl_var_1808, 68 + ig_var_1789, lay_var_1790) = 0.0D0 + wx_var_1762(jl_var_1808, 2, lay_var_1790) * cfc11adj(ig_var_1789) + wx_var_1762(jl_var_1808, 3, lay_var_1790) * cfc12_var_202(ig_var_1789)
          fracs_var_1778(jl_var_1808, 68 + ig_var_1789, lay_var_1790) = fracrefa_var_201(ig_var_1789)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol6
SUBROUTINE rrtm_taumol7(kidia_var_1809, kfdia_var_1810, klev_var_1811, taug_var_1812, p_tauaerl_var_1813, fac00_var_1814, fac01_var_1815, fac10_var_1816, fac11_var_1817, forfac_var_1833, forfrac_var_1832, indfor_var_1831, jp_var_1818, jt_var_1819, jt1_var_1820, oneminus_var_1821, colh2o_var_1822, colo3_var_1823, colco2_var_1824, coldry_var_1825, laytrop_var_1826, selffac_var_1827, selffrac_var_1828, indself_var_1829, fracs_var_1830, rat_h2oo3, rat_h2oo3_1, minorfrac_var_1834, indminor_var_1835)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrtm, ONLY: ng7
  USE yoerrta7, ONLY: absa_var_209, absb_var_210, forref_var_214, fracrefa_var_207, fracrefb_var_208, ka_mco2_var_212, kb_mco2_var_213, selfref_var_211
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1809
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1810
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1811
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1812(kidia_var_1809 : kfdia_var_1810, 140, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1813(kidia_var_1809 : kfdia_var_1810, klev_var_1811, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1814(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1815(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1816(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1817(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1818(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1819(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1820(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1821
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1822(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_1823(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1824(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1825(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1826(kidia_var_1809 : kfdia_var_1810)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1827(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1828(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1829(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1830(kidia_var_1809 : kfdia_var_1810, 140, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oo3(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oo3_1(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1831(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1832(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1833(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1834(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1835(kidia_var_1809 : kfdia_var_1810, klev_var_1811)
  REAL(KIND = 8) :: speccomb_var_1836, speccomb1_var_1837, speccomb_mco2_var_1838, speccomb_planck_var_1839
  INTEGER(KIND = 4) :: ind0_var_1840, ind1_var_1841, inds_var_1842, indf_var_1843, indm_var_1844
  INTEGER(KIND = 4) :: ig_var_1845, js_var_1846, lay_var_1847, js1_var_1848, jpl_var_1849, jmco2_var_1850
  REAL(KIND = 8) :: refrat_planck_a_var_1851, refrat_m_a_var_1852
  REAL(KIND = 8) :: chi_co2_var_1853, ratco2_var_1854, adjfac_var_1855, adjcolco2_var_1856
  REAL(KIND = 8) :: fac000_var_1857, fac100_var_1858, fac200_var_1859, fac010_var_1860, fac110_var_1861, fac210_var_1862, fac001_var_1863, fac101_var_1864, fac201_var_1865, fac011_var_1866, fac111_var_1867, fac211_var_1868
  REAL(KIND = 8) :: p_var_1869, p4_var_1870, fk0_var_1871, fk1_var_1872, fk2_var_1873
  REAL(KIND = 8) :: taufor_var_1874, tauself_var_1875, tau_major_var_1876(12), tau_major1_var_1877(12), co2m1_var_1878, co2m2_var_1879, absco2_var_1880
  REAL(KIND = 8) :: fs_var_1881, specmult_var_1882, specparm_var_1883, fs1_var_1884, specmult1_var_1885, specparm1_var_1886, fpl_var_1887, specmult_planck_var_1888, specparm_planck_var_1889, fmco2_var_1890, specmult_mco2_var_1891, specparm_mco2_var_1892
  INTEGER(KIND = 4) :: laytrop_min_var_1893, laytrop_max_var_1894
  INTEGER(KIND = 4) :: ixc_var_1895(klev_var_1811), ixlow_var_1896(kfdia_var_1810, klev_var_1811), ixhigh_var_1897(kfdia_var_1810, klev_var_1811)
  INTEGER(KIND = 4) :: ich_var_1898, icl_var_1899, ixc0_var_1900, ixp_var_1901, jc_var_1902, jl_var_1903
  laytrop_min_var_1893 = MINVAL(laytrop_var_1826)
  laytrop_max_var_1894 = MAXVAL(laytrop_var_1826)
  ixlow_var_1896 = 0
  ixhigh_var_1897 = 0
  ixc_var_1895 = 0
  DO lay_var_1847 = laytrop_min_var_1893 + 1, laytrop_max_var_1894
    icl_var_1899 = 0
    ich_var_1898 = 0
    DO jc_var_1902 = kidia_var_1809, kfdia_var_1810
      IF (lay_var_1847 <= laytrop_var_1826(jc_var_1902)) THEN
        icl_var_1899 = icl_var_1899 + 1
        ixlow_var_1896(icl_var_1899, lay_var_1847) = jc_var_1902
      ELSE
        ich_var_1898 = ich_var_1898 + 1
        ixhigh_var_1897(ich_var_1898, lay_var_1847) = jc_var_1902
      END IF
    END DO
    ixc_var_1895(lay_var_1847) = icl_var_1899
  END DO
  refrat_planck_a_var_1851 = chi_mls(1, 3) / chi_mls(3, 3)
  refrat_m_a_var_1852 = chi_mls(1, 3) / chi_mls(3, 3)
  DO lay_var_1847 = 1, laytrop_min_var_1893
    DO jl_var_1903 = kidia_var_1809, kfdia_var_1810
      speccomb_var_1836 = colh2o_var_1822(jl_var_1903, lay_var_1847) + rat_h2oo3(jl_var_1903, lay_var_1847) * colo3_var_1823(jl_var_1903, lay_var_1847)
      specparm_var_1883 = MIN(colh2o_var_1822(jl_var_1903, lay_var_1847) / speccomb_var_1836, oneminus_var_1821)
      specmult_var_1882 = 8.0D0 * (specparm_var_1883)
      js_var_1846 = 1 + INT(specmult_var_1882)
      fs_var_1881 = ((specmult_var_1882) - AINT((specmult_var_1882)))
      speccomb1_var_1837 = colh2o_var_1822(jl_var_1903, lay_var_1847) + rat_h2oo3_1(jl_var_1903, lay_var_1847) * colo3_var_1823(jl_var_1903, lay_var_1847)
      specparm1_var_1886 = MIN(colh2o_var_1822(jl_var_1903, lay_var_1847) / speccomb1_var_1837, oneminus_var_1821)
      specmult1_var_1885 = 8.0D0 * (specparm1_var_1886)
      js1_var_1848 = 1 + INT(specmult1_var_1885)
      fs1_var_1884 = ((specmult1_var_1885) - AINT((specmult1_var_1885)))
      speccomb_mco2_var_1838 = colh2o_var_1822(jl_var_1903, lay_var_1847) + refrat_m_a_var_1852 * colo3_var_1823(jl_var_1903, lay_var_1847)
      specparm_mco2_var_1892 = MIN(colh2o_var_1822(jl_var_1903, lay_var_1847) / speccomb_mco2_var_1838, oneminus_var_1821)
      specmult_mco2_var_1891 = 8.0D0 * specparm_mco2_var_1892
      jmco2_var_1850 = 1 + INT(specmult_mco2_var_1891)
      fmco2_var_1890 = ((specmult_mco2_var_1891) - AINT((specmult_mco2_var_1891)))
      chi_co2_var_1853 = colco2_var_1824(jl_var_1903, lay_var_1847) / (coldry_var_1825(jl_var_1903, lay_var_1847))
      ratco2_var_1854 = 1D+20 * chi_co2_var_1853 / chi_mls(2, jp_var_1818(jl_var_1903, lay_var_1847) + 1)
      IF (ratco2_var_1854 .GT. 3.0D0) THEN
        adjfac_var_1855 = 3.0D0 + (ratco2_var_1854 - 3.0D0) ** 0.79D0
        adjcolco2_var_1856 = adjfac_var_1855 * chi_mls(2, jp_var_1818(jl_var_1903, lay_var_1847) + 1) * coldry_var_1825(jl_var_1903, lay_var_1847) * 1D-20
      ELSE
        adjcolco2_var_1856 = colco2_var_1824(jl_var_1903, lay_var_1847)
      END IF
      speccomb_planck_var_1839 = colh2o_var_1822(jl_var_1903, lay_var_1847) + refrat_planck_a_var_1851 * colo3_var_1823(jl_var_1903, lay_var_1847)
      specparm_planck_var_1889 = MIN(colh2o_var_1822(jl_var_1903, lay_var_1847) / speccomb_planck_var_1839, oneminus_var_1821)
      specmult_planck_var_1888 = 8.0D0 * specparm_planck_var_1889
      jpl_var_1849 = 1 + INT(specmult_planck_var_1888)
      fpl_var_1887 = ((specmult_planck_var_1888) - AINT((specmult_planck_var_1888)))
      ind0_var_1840 = ((jp_var_1818(jl_var_1903, lay_var_1847) - 1) * 5 + (jt_var_1819(jl_var_1903, lay_var_1847) - 1)) * nspa_var_237(7) + js_var_1846
      ind1_var_1841 = (jp_var_1818(jl_var_1903, lay_var_1847) * 5 + (jt1_var_1820(jl_var_1903, lay_var_1847) - 1)) * nspa_var_237(7) + js1_var_1848
      inds_var_1842 = indself_var_1829(jl_var_1903, lay_var_1847)
      indf_var_1843 = indfor_var_1831(jl_var_1903, lay_var_1847)
      indm_var_1844 = indminor_var_1835(jl_var_1903, lay_var_1847)
      IF (specparm_var_1883 .LT. 0.125D0) THEN
        p_var_1869 = fs_var_1881 - 1.0D0
        p4_var_1870 = p_var_1869 ** 4
        fk0_var_1871 = p4_var_1870
        fk1_var_1872 = 1.0D0 - p_var_1869 - 2.0D0 * p4_var_1870
        fk2_var_1873 = p_var_1869 + p4_var_1870
        fac000_var_1857 = fk0_var_1871 * fac00_var_1814(jl_var_1903, lay_var_1847)
        fac100_var_1858 = fk1_var_1872 * fac00_var_1814(jl_var_1903, lay_var_1847)
        fac200_var_1859 = fk2_var_1873 * fac00_var_1814(jl_var_1903, lay_var_1847)
        fac010_var_1860 = fk0_var_1871 * fac10_var_1816(jl_var_1903, lay_var_1847)
        fac110_var_1861 = fk1_var_1872 * fac10_var_1816(jl_var_1903, lay_var_1847)
        fac210_var_1862 = fk2_var_1873 * fac10_var_1816(jl_var_1903, lay_var_1847)
      ELSE IF (specparm_var_1883 .GT. 0.875D0) THEN
        p_var_1869 = - fs_var_1881
        p4_var_1870 = p_var_1869 ** 4
        fk0_var_1871 = p4_var_1870
        fk1_var_1872 = 1.0D0 - p_var_1869 - 2.0D0 * p4_var_1870
        fk2_var_1873 = p_var_1869 + p4_var_1870
        fac000_var_1857 = fk0_var_1871 * fac00_var_1814(jl_var_1903, lay_var_1847)
        fac100_var_1858 = fk1_var_1872 * fac00_var_1814(jl_var_1903, lay_var_1847)
        fac200_var_1859 = fk2_var_1873 * fac00_var_1814(jl_var_1903, lay_var_1847)
        fac010_var_1860 = fk0_var_1871 * fac10_var_1816(jl_var_1903, lay_var_1847)
        fac110_var_1861 = fk1_var_1872 * fac10_var_1816(jl_var_1903, lay_var_1847)
        fac210_var_1862 = fk2_var_1873 * fac10_var_1816(jl_var_1903, lay_var_1847)
      ELSE
        fac000_var_1857 = (1.0D0 - fs_var_1881) * fac00_var_1814(jl_var_1903, lay_var_1847)
        fac010_var_1860 = (1.0D0 - fs_var_1881) * fac10_var_1816(jl_var_1903, lay_var_1847)
        fac100_var_1858 = fs_var_1881 * fac00_var_1814(jl_var_1903, lay_var_1847)
        fac110_var_1861 = fs_var_1881 * fac10_var_1816(jl_var_1903, lay_var_1847)
        fac200_var_1859 = 0.0D0
        fac210_var_1862 = 0.0D0
      END IF
      IF (specparm1_var_1886 .LT. 0.125D0) THEN
        p_var_1869 = fs1_var_1884 - 1.0D0
        p4_var_1870 = p_var_1869 ** 4
        fk0_var_1871 = p4_var_1870
        fk1_var_1872 = 1.0D0 - p_var_1869 - 2.0D0 * p4_var_1870
        fk2_var_1873 = p_var_1869 + p4_var_1870
        fac001_var_1863 = fk0_var_1871 * fac01_var_1815(jl_var_1903, lay_var_1847)
        fac101_var_1864 = fk1_var_1872 * fac01_var_1815(jl_var_1903, lay_var_1847)
        fac201_var_1865 = fk2_var_1873 * fac01_var_1815(jl_var_1903, lay_var_1847)
        fac011_var_1866 = fk0_var_1871 * fac11_var_1817(jl_var_1903, lay_var_1847)
        fac111_var_1867 = fk1_var_1872 * fac11_var_1817(jl_var_1903, lay_var_1847)
        fac211_var_1868 = fk2_var_1873 * fac11_var_1817(jl_var_1903, lay_var_1847)
      ELSE IF (specparm1_var_1886 .GT. 0.875D0) THEN
        p_var_1869 = - fs1_var_1884
        p4_var_1870 = p_var_1869 ** 4
        fk0_var_1871 = p4_var_1870
        fk1_var_1872 = 1.0D0 - p_var_1869 - 2.0D0 * p4_var_1870
        fk2_var_1873 = p_var_1869 + p4_var_1870
        fac001_var_1863 = fk0_var_1871 * fac01_var_1815(jl_var_1903, lay_var_1847)
        fac101_var_1864 = fk1_var_1872 * fac01_var_1815(jl_var_1903, lay_var_1847)
        fac201_var_1865 = fk2_var_1873 * fac01_var_1815(jl_var_1903, lay_var_1847)
        fac011_var_1866 = fk0_var_1871 * fac11_var_1817(jl_var_1903, lay_var_1847)
        fac111_var_1867 = fk1_var_1872 * fac11_var_1817(jl_var_1903, lay_var_1847)
        fac211_var_1868 = fk2_var_1873 * fac11_var_1817(jl_var_1903, lay_var_1847)
      ELSE
        fac001_var_1863 = (1.0D0 - fs1_var_1884) * fac01_var_1815(jl_var_1903, lay_var_1847)
        fac011_var_1866 = (1.0D0 - fs1_var_1884) * fac11_var_1817(jl_var_1903, lay_var_1847)
        fac101_var_1864 = fs1_var_1884 * fac01_var_1815(jl_var_1903, lay_var_1847)
        fac111_var_1867 = fs1_var_1884 * fac11_var_1817(jl_var_1903, lay_var_1847)
        fac201_var_1865 = 0.0D0
        fac211_var_1868 = 0.0D0
      END IF
      IF (specparm_var_1883 .LT. 0.125D0) THEN
        tau_major_var_1876(1 : ng7) = speccomb_var_1836 * (fac000_var_1857 * absa_var_209(ind0_var_1840, 1 : 12) + fac100_var_1858 * absa_var_209(ind0_var_1840 + 1, 1 : 12) + fac200_var_1859 * absa_var_209(ind0_var_1840 + 2, 1 : 12) + fac010_var_1860 * absa_var_209(ind0_var_1840 + 9, 1 : 12) + fac110_var_1861 * absa_var_209(ind0_var_1840 + 10, 1 : 12) + fac210_var_1862 * absa_var_209(ind0_var_1840 + 11, 1 : 12))
      ELSE IF (specparm_var_1883 .GT. 0.875D0) THEN
        tau_major_var_1876(1 : ng7) = speccomb_var_1836 * (fac200_var_1859 * absa_var_209(ind0_var_1840 - 1, 1 : 12) + fac100_var_1858 * absa_var_209(ind0_var_1840, 1 : 12) + fac000_var_1857 * absa_var_209(ind0_var_1840 + 1, 1 : 12) + fac210_var_1862 * absa_var_209(ind0_var_1840 + 8, 1 : 12) + fac110_var_1861 * absa_var_209(ind0_var_1840 + 9, 1 : 12) + fac010_var_1860 * absa_var_209(ind0_var_1840 + 10, 1 : 12))
      ELSE
        tau_major_var_1876(1 : ng7) = speccomb_var_1836 * (fac000_var_1857 * absa_var_209(ind0_var_1840, 1 : 12) + fac100_var_1858 * absa_var_209(ind0_var_1840 + 1, 1 : 12) + fac010_var_1860 * absa_var_209(ind0_var_1840 + 9, 1 : 12) + fac110_var_1861 * absa_var_209(ind0_var_1840 + 10, 1 : 12))
      END IF
      IF (specparm1_var_1886 .LT. 0.125D0) THEN
        tau_major1_var_1877(1 : ng7) = speccomb1_var_1837 * (fac001_var_1863 * absa_var_209(ind1_var_1841, 1 : 12) + fac101_var_1864 * absa_var_209(ind1_var_1841 + 1, 1 : 12) + fac201_var_1865 * absa_var_209(ind1_var_1841 + 2, 1 : 12) + fac011_var_1866 * absa_var_209(ind1_var_1841 + 9, 1 : 12) + fac111_var_1867 * absa_var_209(ind1_var_1841 + 10, 1 : 12) + fac211_var_1868 * absa_var_209(ind1_var_1841 + 11, 1 : 12))
      ELSE IF (specparm1_var_1886 .GT. 0.875D0) THEN
        tau_major1_var_1877(1 : ng7) = speccomb1_var_1837 * (fac201_var_1865 * absa_var_209(ind1_var_1841 - 1, 1 : 12) + fac101_var_1864 * absa_var_209(ind1_var_1841, 1 : 12) + fac001_var_1863 * absa_var_209(ind1_var_1841 + 1, 1 : 12) + fac211_var_1868 * absa_var_209(ind1_var_1841 + 8, 1 : 12) + fac111_var_1867 * absa_var_209(ind1_var_1841 + 9, 1 : 12) + fac011_var_1866 * absa_var_209(ind1_var_1841 + 10, 1 : 12))
      ELSE
        tau_major1_var_1877(1 : ng7) = speccomb1_var_1837 * (fac001_var_1863 * absa_var_209(ind1_var_1841, 1 : 12) + fac101_var_1864 * absa_var_209(ind1_var_1841 + 1, 1 : 12) + fac011_var_1866 * absa_var_209(ind1_var_1841 + 9, 1 : 12) + fac111_var_1867 * absa_var_209(ind1_var_1841 + 10, 1 : 12))
      END IF
      DO ig_var_1845 = 1, 12
        tauself_var_1875 = selffac_var_1827(jl_var_1903, lay_var_1847) * (selfref_var_211(inds_var_1842, ig_var_1845) + selffrac_var_1828(jl_var_1903, lay_var_1847) * (selfref_var_211(inds_var_1842 + 1, ig_var_1845) - selfref_var_211(inds_var_1842, ig_var_1845)))
        taufor_var_1874 = forfac_var_1833(jl_var_1903, lay_var_1847) * (forref_var_214(indf_var_1843, ig_var_1845) + forfrac_var_1832(jl_var_1903, lay_var_1847) * (forref_var_214(indf_var_1843 + 1, ig_var_1845) - forref_var_214(indf_var_1843, ig_var_1845)))
        co2m1_var_1878 = ka_mco2_var_212(jmco2_var_1850, indm_var_1844, ig_var_1845) + fmco2_var_1890 * (ka_mco2_var_212(jmco2_var_1850 + 1, indm_var_1844, ig_var_1845) - ka_mco2_var_212(jmco2_var_1850, indm_var_1844, ig_var_1845))
        co2m2_var_1879 = ka_mco2_var_212(jmco2_var_1850, indm_var_1844 + 1, ig_var_1845) + fmco2_var_1890 * (ka_mco2_var_212(jmco2_var_1850 + 1, indm_var_1844 + 1, ig_var_1845) - ka_mco2_var_212(jmco2_var_1850, indm_var_1844 + 1, ig_var_1845))
        absco2_var_1880 = co2m1_var_1878 + minorfrac_var_1834(jl_var_1903, lay_var_1847) * (co2m2_var_1879 - co2m1_var_1878)
        taug_var_1812(jl_var_1903, 76 + ig_var_1845, lay_var_1847) = tau_major_var_1876(ig_var_1845) + tau_major1_var_1877(ig_var_1845) + tauself_var_1875 + taufor_var_1874 + adjcolco2_var_1856 * absco2_var_1880
        fracs_var_1830(jl_var_1903, 76 + ig_var_1845, lay_var_1847) = fracrefa_var_207(ig_var_1845, jpl_var_1849) + fpl_var_1887 * (fracrefa_var_207(ig_var_1845, jpl_var_1849 + 1) - fracrefa_var_207(ig_var_1845, jpl_var_1849))
      END DO
    END DO
  END DO
  DO lay_var_1847 = laytrop_max_var_1894 + 1, klev_var_1811
    DO jl_var_1903 = kidia_var_1809, kfdia_var_1810
      chi_co2_var_1853 = colco2_var_1824(jl_var_1903, lay_var_1847) / (coldry_var_1825(jl_var_1903, lay_var_1847))
      ratco2_var_1854 = 1D+20 * chi_co2_var_1853 / chi_mls(2, jp_var_1818(jl_var_1903, lay_var_1847) + 1)
      IF (ratco2_var_1854 .GT. 3.0D0) THEN
        adjfac_var_1855 = 2.0D0 + (ratco2_var_1854 - 2.0D0) ** 0.79D0
        adjcolco2_var_1856 = adjfac_var_1855 * chi_mls(2, jp_var_1818(jl_var_1903, lay_var_1847) + 1) * coldry_var_1825(jl_var_1903, lay_var_1847) * 1D-20
      ELSE
        adjcolco2_var_1856 = colco2_var_1824(jl_var_1903, lay_var_1847)
      END IF
      ind0_var_1840 = ((jp_var_1818(jl_var_1903, lay_var_1847) - 13) * 5 + (jt_var_1819(jl_var_1903, lay_var_1847) - 1)) * nspb_var_238(7) + 1
      ind1_var_1841 = ((jp_var_1818(jl_var_1903, lay_var_1847) - 12) * 5 + (jt1_var_1820(jl_var_1903, lay_var_1847) - 1)) * nspb_var_238(7) + 1
      indm_var_1844 = indminor_var_1835(jl_var_1903, lay_var_1847)
      DO ig_var_1845 = 1, 12
        absco2_var_1880 = kb_mco2_var_213(indm_var_1844, ig_var_1845) + minorfrac_var_1834(jl_var_1903, lay_var_1847) * (kb_mco2_var_213(indm_var_1844 + 1, ig_var_1845) - kb_mco2_var_213(indm_var_1844, ig_var_1845))
        taug_var_1812(jl_var_1903, 76 + ig_var_1845, lay_var_1847) = colo3_var_1823(jl_var_1903, lay_var_1847) * (fac00_var_1814(jl_var_1903, lay_var_1847) * absb_var_210(ind0_var_1840, ig_var_1845) + fac10_var_1816(jl_var_1903, lay_var_1847) * absb_var_210(ind0_var_1840 + 1, ig_var_1845) + fac01_var_1815(jl_var_1903, lay_var_1847) * absb_var_210(ind1_var_1841, ig_var_1845) + fac11_var_1817(jl_var_1903, lay_var_1847) * absb_var_210(ind1_var_1841 + 1, ig_var_1845)) + adjcolco2_var_1856 * absco2_var_1880
        fracs_var_1830(jl_var_1903, 76 + ig_var_1845, lay_var_1847) = fracrefb_var_208(ig_var_1845)
      END DO
    END DO
  END DO
  DO lay_var_1847 = laytrop_max_var_1894 + 1, klev_var_1811
    DO jl_var_1903 = kidia_var_1809, kfdia_var_1810
      taug_var_1812(jl_var_1903, 82, lay_var_1847) = taug_var_1812(jl_var_1903, 82, lay_var_1847) * 0.92D0
      taug_var_1812(jl_var_1903, 83, lay_var_1847) = taug_var_1812(jl_var_1903, 83, lay_var_1847) * 0.88D0
      taug_var_1812(jl_var_1903, 84, lay_var_1847) = taug_var_1812(jl_var_1903, 84, lay_var_1847) * 1.07D0
      taug_var_1812(jl_var_1903, 85, lay_var_1847) = taug_var_1812(jl_var_1903, 85, lay_var_1847) * 1.1D0
      taug_var_1812(jl_var_1903, 86, lay_var_1847) = taug_var_1812(jl_var_1903, 86, lay_var_1847) * 0.99D0
      taug_var_1812(jl_var_1903, 87, lay_var_1847) = taug_var_1812(jl_var_1903, 87, lay_var_1847) * 0.855D0
    END DO
  END DO
  IF (laytrop_max_var_1894 /= laytrop_min_var_1893) THEN
    DO lay_var_1847 = laytrop_min_var_1893 + 1, laytrop_max_var_1894
      ixc0_var_1900 = ixc_var_1895(lay_var_1847)
      DO ixp_var_1901 = 1, ixc0_var_1900
        jl_var_1903 = ixlow_var_1896(ixp_var_1901, lay_var_1847)
        speccomb_var_1836 = colh2o_var_1822(jl_var_1903, lay_var_1847) + rat_h2oo3(jl_var_1903, lay_var_1847) * colo3_var_1823(jl_var_1903, lay_var_1847)
        specparm_var_1883 = MIN(colh2o_var_1822(jl_var_1903, lay_var_1847) / speccomb_var_1836, oneminus_var_1821)
        specmult_var_1882 = 8.0D0 * (specparm_var_1883)
        js_var_1846 = 1 + INT(specmult_var_1882)
        fs_var_1881 = ((specmult_var_1882) - AINT((specmult_var_1882)))
        speccomb1_var_1837 = colh2o_var_1822(jl_var_1903, lay_var_1847) + rat_h2oo3_1(jl_var_1903, lay_var_1847) * colo3_var_1823(jl_var_1903, lay_var_1847)
        specparm1_var_1886 = MIN(colh2o_var_1822(jl_var_1903, lay_var_1847) / speccomb1_var_1837, oneminus_var_1821)
        specmult1_var_1885 = 8.0D0 * (specparm1_var_1886)
        js1_var_1848 = 1 + INT(specmult1_var_1885)
        fs1_var_1884 = ((specmult1_var_1885) - AINT((specmult1_var_1885)))
        speccomb_mco2_var_1838 = colh2o_var_1822(jl_var_1903, lay_var_1847) + refrat_m_a_var_1852 * colo3_var_1823(jl_var_1903, lay_var_1847)
        specparm_mco2_var_1892 = MIN(colh2o_var_1822(jl_var_1903, lay_var_1847) / speccomb_mco2_var_1838, oneminus_var_1821)
        specmult_mco2_var_1891 = 8.0D0 * specparm_mco2_var_1892
        jmco2_var_1850 = 1 + INT(specmult_mco2_var_1891)
        fmco2_var_1890 = ((specmult_mco2_var_1891) - AINT((specmult_mco2_var_1891)))
        chi_co2_var_1853 = colco2_var_1824(jl_var_1903, lay_var_1847) / (coldry_var_1825(jl_var_1903, lay_var_1847))
        ratco2_var_1854 = 1D+20 * chi_co2_var_1853 / chi_mls(2, jp_var_1818(jl_var_1903, lay_var_1847) + 1)
        IF (ratco2_var_1854 .GT. 3.0D0) THEN
          adjfac_var_1855 = 3.0D0 + (ratco2_var_1854 - 3.0D0) ** 0.79D0
          adjcolco2_var_1856 = adjfac_var_1855 * chi_mls(2, jp_var_1818(jl_var_1903, lay_var_1847) + 1) * coldry_var_1825(jl_var_1903, lay_var_1847) * 1D-20
        ELSE
          adjcolco2_var_1856 = colco2_var_1824(jl_var_1903, lay_var_1847)
        END IF
        speccomb_planck_var_1839 = colh2o_var_1822(jl_var_1903, lay_var_1847) + refrat_planck_a_var_1851 * colo3_var_1823(jl_var_1903, lay_var_1847)
        specparm_planck_var_1889 = MIN(colh2o_var_1822(jl_var_1903, lay_var_1847) / speccomb_planck_var_1839, oneminus_var_1821)
        specmult_planck_var_1888 = 8.0D0 * specparm_planck_var_1889
        jpl_var_1849 = 1 + INT(specmult_planck_var_1888)
        fpl_var_1887 = ((specmult_planck_var_1888) - AINT((specmult_planck_var_1888)))
        ind0_var_1840 = ((jp_var_1818(jl_var_1903, lay_var_1847) - 1) * 5 + (jt_var_1819(jl_var_1903, lay_var_1847) - 1)) * nspa_var_237(7) + js_var_1846
        ind1_var_1841 = (jp_var_1818(jl_var_1903, lay_var_1847) * 5 + (jt1_var_1820(jl_var_1903, lay_var_1847) - 1)) * nspa_var_237(7) + js1_var_1848
        inds_var_1842 = indself_var_1829(jl_var_1903, lay_var_1847)
        indf_var_1843 = indfor_var_1831(jl_var_1903, lay_var_1847)
        indm_var_1844 = indminor_var_1835(jl_var_1903, lay_var_1847)
        IF (specparm_var_1883 .LT. 0.125D0) THEN
          p_var_1869 = fs_var_1881 - 1.0D0
          p4_var_1870 = p_var_1869 ** 4
          fk0_var_1871 = p4_var_1870
          fk1_var_1872 = 1.0D0 - p_var_1869 - 2.0D0 * p4_var_1870
          fk2_var_1873 = p_var_1869 + p4_var_1870
          fac000_var_1857 = fk0_var_1871 * fac00_var_1814(jl_var_1903, lay_var_1847)
          fac100_var_1858 = fk1_var_1872 * fac00_var_1814(jl_var_1903, lay_var_1847)
          fac200_var_1859 = fk2_var_1873 * fac00_var_1814(jl_var_1903, lay_var_1847)
          fac010_var_1860 = fk0_var_1871 * fac10_var_1816(jl_var_1903, lay_var_1847)
          fac110_var_1861 = fk1_var_1872 * fac10_var_1816(jl_var_1903, lay_var_1847)
          fac210_var_1862 = fk2_var_1873 * fac10_var_1816(jl_var_1903, lay_var_1847)
        ELSE IF (specparm_var_1883 .GT. 0.875D0) THEN
          p_var_1869 = - fs_var_1881
          p4_var_1870 = p_var_1869 ** 4
          fk0_var_1871 = p4_var_1870
          fk1_var_1872 = 1.0D0 - p_var_1869 - 2.0D0 * p4_var_1870
          fk2_var_1873 = p_var_1869 + p4_var_1870
          fac000_var_1857 = fk0_var_1871 * fac00_var_1814(jl_var_1903, lay_var_1847)
          fac100_var_1858 = fk1_var_1872 * fac00_var_1814(jl_var_1903, lay_var_1847)
          fac200_var_1859 = fk2_var_1873 * fac00_var_1814(jl_var_1903, lay_var_1847)
          fac010_var_1860 = fk0_var_1871 * fac10_var_1816(jl_var_1903, lay_var_1847)
          fac110_var_1861 = fk1_var_1872 * fac10_var_1816(jl_var_1903, lay_var_1847)
          fac210_var_1862 = fk2_var_1873 * fac10_var_1816(jl_var_1903, lay_var_1847)
        ELSE
          fac000_var_1857 = (1.0D0 - fs_var_1881) * fac00_var_1814(jl_var_1903, lay_var_1847)
          fac010_var_1860 = (1.0D0 - fs_var_1881) * fac10_var_1816(jl_var_1903, lay_var_1847)
          fac100_var_1858 = fs_var_1881 * fac00_var_1814(jl_var_1903, lay_var_1847)
          fac110_var_1861 = fs_var_1881 * fac10_var_1816(jl_var_1903, lay_var_1847)
          fac200_var_1859 = 0.0D0
          fac210_var_1862 = 0.0D0
        END IF
        IF (specparm1_var_1886 .LT. 0.125D0) THEN
          p_var_1869 = fs1_var_1884 - 1.0D0
          p4_var_1870 = p_var_1869 ** 4
          fk0_var_1871 = p4_var_1870
          fk1_var_1872 = 1.0D0 - p_var_1869 - 2.0D0 * p4_var_1870
          fk2_var_1873 = p_var_1869 + p4_var_1870
          fac001_var_1863 = fk0_var_1871 * fac01_var_1815(jl_var_1903, lay_var_1847)
          fac101_var_1864 = fk1_var_1872 * fac01_var_1815(jl_var_1903, lay_var_1847)
          fac201_var_1865 = fk2_var_1873 * fac01_var_1815(jl_var_1903, lay_var_1847)
          fac011_var_1866 = fk0_var_1871 * fac11_var_1817(jl_var_1903, lay_var_1847)
          fac111_var_1867 = fk1_var_1872 * fac11_var_1817(jl_var_1903, lay_var_1847)
          fac211_var_1868 = fk2_var_1873 * fac11_var_1817(jl_var_1903, lay_var_1847)
        ELSE IF (specparm1_var_1886 .GT. 0.875D0) THEN
          p_var_1869 = - fs1_var_1884
          p4_var_1870 = p_var_1869 ** 4
          fk0_var_1871 = p4_var_1870
          fk1_var_1872 = 1.0D0 - p_var_1869 - 2.0D0 * p4_var_1870
          fk2_var_1873 = p_var_1869 + p4_var_1870
          fac001_var_1863 = fk0_var_1871 * fac01_var_1815(jl_var_1903, lay_var_1847)
          fac101_var_1864 = fk1_var_1872 * fac01_var_1815(jl_var_1903, lay_var_1847)
          fac201_var_1865 = fk2_var_1873 * fac01_var_1815(jl_var_1903, lay_var_1847)
          fac011_var_1866 = fk0_var_1871 * fac11_var_1817(jl_var_1903, lay_var_1847)
          fac111_var_1867 = fk1_var_1872 * fac11_var_1817(jl_var_1903, lay_var_1847)
          fac211_var_1868 = fk2_var_1873 * fac11_var_1817(jl_var_1903, lay_var_1847)
        ELSE
          fac001_var_1863 = (1.0D0 - fs1_var_1884) * fac01_var_1815(jl_var_1903, lay_var_1847)
          fac011_var_1866 = (1.0D0 - fs1_var_1884) * fac11_var_1817(jl_var_1903, lay_var_1847)
          fac101_var_1864 = fs1_var_1884 * fac01_var_1815(jl_var_1903, lay_var_1847)
          fac111_var_1867 = fs1_var_1884 * fac11_var_1817(jl_var_1903, lay_var_1847)
          fac201_var_1865 = 0.0D0
          fac211_var_1868 = 0.0D0
        END IF
        IF (specparm_var_1883 .LT. 0.125D0) THEN
          tau_major_var_1876(1 : ng7) = speccomb_var_1836 * (fac000_var_1857 * absa_var_209(ind0_var_1840, 1 : 12) + fac100_var_1858 * absa_var_209(ind0_var_1840 + 1, 1 : 12) + fac200_var_1859 * absa_var_209(ind0_var_1840 + 2, 1 : 12) + fac010_var_1860 * absa_var_209(ind0_var_1840 + 9, 1 : 12) + fac110_var_1861 * absa_var_209(ind0_var_1840 + 10, 1 : 12) + fac210_var_1862 * absa_var_209(ind0_var_1840 + 11, 1 : 12))
        ELSE IF (specparm_var_1883 .GT. 0.875D0) THEN
          tau_major_var_1876(1 : ng7) = speccomb_var_1836 * (fac200_var_1859 * absa_var_209(ind0_var_1840 - 1, 1 : 12) + fac100_var_1858 * absa_var_209(ind0_var_1840, 1 : 12) + fac000_var_1857 * absa_var_209(ind0_var_1840 + 1, 1 : 12) + fac210_var_1862 * absa_var_209(ind0_var_1840 + 8, 1 : 12) + fac110_var_1861 * absa_var_209(ind0_var_1840 + 9, 1 : 12) + fac010_var_1860 * absa_var_209(ind0_var_1840 + 10, 1 : 12))
        ELSE
          tau_major_var_1876(1 : ng7) = speccomb_var_1836 * (fac000_var_1857 * absa_var_209(ind0_var_1840, 1 : 12) + fac100_var_1858 * absa_var_209(ind0_var_1840 + 1, 1 : 12) + fac010_var_1860 * absa_var_209(ind0_var_1840 + 9, 1 : 12) + fac110_var_1861 * absa_var_209(ind0_var_1840 + 10, 1 : 12))
        END IF
        IF (specparm1_var_1886 .LT. 0.125D0) THEN
          tau_major1_var_1877(1 : ng7) = speccomb1_var_1837 * (fac001_var_1863 * absa_var_209(ind1_var_1841, 1 : 12) + fac101_var_1864 * absa_var_209(ind1_var_1841 + 1, 1 : 12) + fac201_var_1865 * absa_var_209(ind1_var_1841 + 2, 1 : 12) + fac011_var_1866 * absa_var_209(ind1_var_1841 + 9, 1 : 12) + fac111_var_1867 * absa_var_209(ind1_var_1841 + 10, 1 : 12) + fac211_var_1868 * absa_var_209(ind1_var_1841 + 11, 1 : 12))
        ELSE IF (specparm1_var_1886 .GT. 0.875D0) THEN
          tau_major1_var_1877(1 : ng7) = speccomb1_var_1837 * (fac201_var_1865 * absa_var_209(ind1_var_1841 - 1, 1 : 12) + fac101_var_1864 * absa_var_209(ind1_var_1841, 1 : 12) + fac001_var_1863 * absa_var_209(ind1_var_1841 + 1, 1 : 12) + fac211_var_1868 * absa_var_209(ind1_var_1841 + 8, 1 : 12) + fac111_var_1867 * absa_var_209(ind1_var_1841 + 9, 1 : 12) + fac011_var_1866 * absa_var_209(ind1_var_1841 + 10, 1 : 12))
        ELSE
          tau_major1_var_1877(1 : ng7) = speccomb1_var_1837 * (fac001_var_1863 * absa_var_209(ind1_var_1841, 1 : 12) + fac101_var_1864 * absa_var_209(ind1_var_1841 + 1, 1 : 12) + fac011_var_1866 * absa_var_209(ind1_var_1841 + 9, 1 : 12) + fac111_var_1867 * absa_var_209(ind1_var_1841 + 10, 1 : 12))
        END IF
        DO ig_var_1845 = 1, 12
          tauself_var_1875 = selffac_var_1827(jl_var_1903, lay_var_1847) * (selfref_var_211(inds_var_1842, ig_var_1845) + selffrac_var_1828(jl_var_1903, lay_var_1847) * (selfref_var_211(inds_var_1842 + 1, ig_var_1845) - selfref_var_211(inds_var_1842, ig_var_1845)))
          taufor_var_1874 = forfac_var_1833(jl_var_1903, lay_var_1847) * (forref_var_214(indf_var_1843, ig_var_1845) + forfrac_var_1832(jl_var_1903, lay_var_1847) * (forref_var_214(indf_var_1843 + 1, ig_var_1845) - forref_var_214(indf_var_1843, ig_var_1845)))
          co2m1_var_1878 = ka_mco2_var_212(jmco2_var_1850, indm_var_1844, ig_var_1845) + fmco2_var_1890 * (ka_mco2_var_212(jmco2_var_1850 + 1, indm_var_1844, ig_var_1845) - ka_mco2_var_212(jmco2_var_1850, indm_var_1844, ig_var_1845))
          co2m2_var_1879 = ka_mco2_var_212(jmco2_var_1850, indm_var_1844 + 1, ig_var_1845) + fmco2_var_1890 * (ka_mco2_var_212(jmco2_var_1850 + 1, indm_var_1844 + 1, ig_var_1845) - ka_mco2_var_212(jmco2_var_1850, indm_var_1844 + 1, ig_var_1845))
          absco2_var_1880 = co2m1_var_1878 + minorfrac_var_1834(jl_var_1903, lay_var_1847) * (co2m2_var_1879 - co2m1_var_1878)
          taug_var_1812(jl_var_1903, 76 + ig_var_1845, lay_var_1847) = tau_major_var_1876(ig_var_1845) + tau_major1_var_1877(ig_var_1845) + tauself_var_1875 + taufor_var_1874 + adjcolco2_var_1856 * absco2_var_1880
          fracs_var_1830(jl_var_1903, 76 + ig_var_1845, lay_var_1847) = fracrefa_var_207(ig_var_1845, jpl_var_1849) + fpl_var_1887 * (fracrefa_var_207(ig_var_1845, jpl_var_1849 + 1) - fracrefa_var_207(ig_var_1845, jpl_var_1849))
        END DO
      END DO
      ixc0_var_1900 = kfdia_var_1810 - kidia_var_1809 + 1 - ixc0_var_1900
      DO ixp_var_1901 = 1, ixc0_var_1900
        jl_var_1903 = ixhigh_var_1897(ixp_var_1901, lay_var_1847)
        chi_co2_var_1853 = colco2_var_1824(jl_var_1903, lay_var_1847) / (coldry_var_1825(jl_var_1903, lay_var_1847))
        ratco2_var_1854 = 1D+20 * chi_co2_var_1853 / chi_mls(2, jp_var_1818(jl_var_1903, lay_var_1847) + 1)
        IF (ratco2_var_1854 .GT. 3.0D0) THEN
          adjfac_var_1855 = 2.0D0 + (ratco2_var_1854 - 2.0D0) ** 0.79D0
          adjcolco2_var_1856 = adjfac_var_1855 * chi_mls(2, jp_var_1818(jl_var_1903, lay_var_1847) + 1) * coldry_var_1825(jl_var_1903, lay_var_1847) * 1D-20
        ELSE
          adjcolco2_var_1856 = colco2_var_1824(jl_var_1903, lay_var_1847)
        END IF
        ind0_var_1840 = ((jp_var_1818(jl_var_1903, lay_var_1847) - 13) * 5 + (jt_var_1819(jl_var_1903, lay_var_1847) - 1)) * nspb_var_238(7) + 1
        ind1_var_1841 = ((jp_var_1818(jl_var_1903, lay_var_1847) - 12) * 5 + (jt1_var_1820(jl_var_1903, lay_var_1847) - 1)) * nspb_var_238(7) + 1
        indm_var_1844 = indminor_var_1835(jl_var_1903, lay_var_1847)
        DO ig_var_1845 = 1, 12
          absco2_var_1880 = kb_mco2_var_213(indm_var_1844, ig_var_1845) + minorfrac_var_1834(jl_var_1903, lay_var_1847) * (kb_mco2_var_213(indm_var_1844 + 1, ig_var_1845) - kb_mco2_var_213(indm_var_1844, ig_var_1845))
          taug_var_1812(jl_var_1903, 76 + ig_var_1845, lay_var_1847) = colo3_var_1823(jl_var_1903, lay_var_1847) * (fac00_var_1814(jl_var_1903, lay_var_1847) * absb_var_210(ind0_var_1840, ig_var_1845) + fac10_var_1816(jl_var_1903, lay_var_1847) * absb_var_210(ind0_var_1840 + 1, ig_var_1845) + fac01_var_1815(jl_var_1903, lay_var_1847) * absb_var_210(ind1_var_1841, ig_var_1845) + fac11_var_1817(jl_var_1903, lay_var_1847) * absb_var_210(ind1_var_1841 + 1, ig_var_1845)) + adjcolco2_var_1856 * absco2_var_1880
          fracs_var_1830(jl_var_1903, 76 + ig_var_1845, lay_var_1847) = fracrefb_var_208(ig_var_1845)
        END DO
      END DO
      DO ixp_var_1901 = 1, ixc0_var_1900
        jl_var_1903 = ixhigh_var_1897(ixp_var_1901, lay_var_1847)
        taug_var_1812(jl_var_1903, 82, lay_var_1847) = taug_var_1812(jl_var_1903, 82, lay_var_1847) * 0.92D0
        taug_var_1812(jl_var_1903, 83, lay_var_1847) = taug_var_1812(jl_var_1903, 83, lay_var_1847) * 0.88D0
        taug_var_1812(jl_var_1903, 84, lay_var_1847) = taug_var_1812(jl_var_1903, 84, lay_var_1847) * 1.07D0
        taug_var_1812(jl_var_1903, 85, lay_var_1847) = taug_var_1812(jl_var_1903, 85, lay_var_1847) * 1.1D0
        taug_var_1812(jl_var_1903, 86, lay_var_1847) = taug_var_1812(jl_var_1903, 86, lay_var_1847) * 0.99D0
        taug_var_1812(jl_var_1903, 87, lay_var_1847) = taug_var_1812(jl_var_1903, 87, lay_var_1847) * 0.855D0
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol7
SUBROUTINE rrtm_taumol3(kidia_var_1904, kfdia_var_1905, klev_var_1906, taug_var_1907, p_tauaerl_var_1908, fac00_var_1909, fac01_var_1910, fac10_var_1911, fac11_var_1912, forfac_var_1913, forfrac_var_1930, indfor_var_1929, jp_var_1914, jt_var_1915, jt1_var_1916, oneminus_var_1917, colh2o_var_1918, colco2_var_1919, coln2o_var_1920, coldry_var_1921, laytrop_var_1922, selffac_var_1923, selffrac_var_1924, indself_var_1925, fracs_var_1926, rat_h2oco2_var_1927, rat_h2oco2_1_var_1928, minorfrac_var_1931, indminor_var_1932)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrtm, ONLY: ng3
  USE yoerrta3, ONLY: absa_var_184, absb_var_185, forref_var_187, fracrefa_var_180, fracrefb_var_181, ka_mn2o_var_182, kb_mn2o_var_183, selfref_var_186
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_1904
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_1905
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_1906
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_1907(kidia_var_1904 : kfdia_var_1905, 140, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_1908(kidia_var_1904 : kfdia_var_1905, klev_var_1906, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_1909(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_1910(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_1911(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_1912(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_1913(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_1914(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_1915(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_1916(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_1917
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_1918(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_1919(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_1920(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_1921(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_1922(kidia_var_1904 : kfdia_var_1905)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_1923(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_1924(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_1925(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_1926(kidia_var_1904 : kfdia_var_1905, 140, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_1927(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_1928(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_1929(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_1930(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_1931(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_1932(kidia_var_1904 : kfdia_var_1905, klev_var_1906)
  REAL(KIND = 8) :: speccomb_var_1933, speccomb1_var_1934, speccomb_mn2o_var_1935, speccomb_planck_var_1936
  REAL(KIND = 8) :: refrat_planck_a_var_1937, refrat_planck_b_var_1938, refrat_m_a_var_1939, refrat_m_b
  INTEGER(KIND = 4) :: ind0_var_1940, ind1_var_1941, inds_var_1942, indf_var_1943, indm_var_1944
  INTEGER(KIND = 4) :: ig_var_1945, js_var_1946, lay_var_1947, js1_var_1948, jmn2o_var_1949, jpl_var_1950
  REAL(KIND = 8) :: fs_var_1951, specmult_var_1952, specparm_var_1953, fs1_var_1954, specmult1_var_1955, specparm1_var_1956, fmn2o_var_1957, specmult_mn2o_var_1958, specparm_mn2o_var_1959, fpl_var_1960, specmult_planck_var_1961, specparm_planck_var_1962
  REAL(KIND = 8) :: adjfac_var_1963, adjcoln2o_var_1964, ratn2o_var_1965, chi_n2o_var_1966
  REAL(KIND = 8) :: fac000_var_1967, fac100_var_1968, fac200_var_1969, fac010_var_1970, fac110_var_1971, fac210_var_1972, fac001_var_1973, fac101_var_1974, fac201_var_1975, fac011_var_1976, fac111_var_1977, fac211_var_1978
  REAL(KIND = 8) :: p_var_1979, p4_var_1980, fk0_var_1981, fk1_var_1982, fk2_var_1983
  REAL(KIND = 8) :: taufor_var_1984, tauself_var_1985, n2om1_var_1986, n2om2_var_1987, absn2o_var_1988, tau_major_var_1989(16), tau_major1_var_1990(16)
  INTEGER(KIND = 4) :: laytrop_min_var_1991, laytrop_max_var_1992
  INTEGER(KIND = 4) :: ixc_var_1993(klev_var_1906), ixlow_var_1994(kfdia_var_1905, klev_var_1906), ixhigh_var_1995(kfdia_var_1905, klev_var_1906)
  INTEGER(KIND = 4) :: ich_var_1996, icl_var_1997, ixc0_var_1998, ixp_var_1999, jc_var_2000, jl_var_2001
  laytrop_min_var_1991 = MINVAL(laytrop_var_1922)
  laytrop_max_var_1992 = MAXVAL(laytrop_var_1922)
  ixlow_var_1994 = 0
  ixhigh_var_1995 = 0
  ixc_var_1993 = 0
  DO lay_var_1947 = laytrop_min_var_1991 + 1, laytrop_max_var_1992
    icl_var_1997 = 0
    ich_var_1996 = 0
    DO jc_var_2000 = kidia_var_1904, kfdia_var_1905
      IF (lay_var_1947 <= laytrop_var_1922(jc_var_2000)) THEN
        icl_var_1997 = icl_var_1997 + 1
        ixlow_var_1994(icl_var_1997, lay_var_1947) = jc_var_2000
      ELSE
        ich_var_1996 = ich_var_1996 + 1
        ixhigh_var_1995(ich_var_1996, lay_var_1947) = jc_var_2000
      END IF
    END DO
    ixc_var_1993(lay_var_1947) = icl_var_1997
  END DO
  refrat_planck_a_var_1937 = chi_mls(1, 9) / chi_mls(2, 9)
  refrat_planck_b_var_1938 = chi_mls(1, 13) / chi_mls(2, 13)
  refrat_m_a_var_1939 = chi_mls(1, 3) / chi_mls(2, 3)
  refrat_m_b = chi_mls(1, 13) / chi_mls(2, 13)
  DO lay_var_1947 = 1, laytrop_min_var_1991
    DO jl_var_2001 = kidia_var_1904, kfdia_var_1905
      speccomb_var_1933 = colh2o_var_1918(jl_var_2001, lay_var_1947) + rat_h2oco2_var_1927(jl_var_2001, lay_var_1947) * colco2_var_1919(jl_var_2001, lay_var_1947)
      specparm_var_1953 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_var_1933, oneminus_var_1917)
      specmult_var_1952 = 8.0D0 * (specparm_var_1953)
      js_var_1946 = 1 + INT(specmult_var_1952)
      fs_var_1951 = ((specmult_var_1952) - AINT((specmult_var_1952)))
      speccomb1_var_1934 = colh2o_var_1918(jl_var_2001, lay_var_1947) + rat_h2oco2_1_var_1928(jl_var_2001, lay_var_1947) * colco2_var_1919(jl_var_2001, lay_var_1947)
      specparm1_var_1956 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb1_var_1934, oneminus_var_1917)
      specmult1_var_1955 = 8.0D0 * (specparm1_var_1956)
      js1_var_1948 = 1 + INT(specmult1_var_1955)
      fs1_var_1954 = ((specmult1_var_1955) - AINT((specmult1_var_1955)))
      speccomb_mn2o_var_1935 = colh2o_var_1918(jl_var_2001, lay_var_1947) + refrat_m_a_var_1939 * colco2_var_1919(jl_var_2001, lay_var_1947)
      specparm_mn2o_var_1959 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_mn2o_var_1935, oneminus_var_1917)
      specmult_mn2o_var_1958 = 8.0D0 * specparm_mn2o_var_1959
      jmn2o_var_1949 = 1 + INT(specmult_mn2o_var_1958)
      fmn2o_var_1957 = ((specmult_mn2o_var_1958) - AINT((specmult_mn2o_var_1958)))
      chi_n2o_var_1966 = coln2o_var_1920(jl_var_2001, lay_var_1947) / coldry_var_1921(jl_var_2001, lay_var_1947)
      ratn2o_var_1965 = 1D+20 * chi_n2o_var_1966 / chi_mls(4, jp_var_1914(jl_var_2001, lay_var_1947) + 1)
      IF (ratn2o_var_1965 .GT. 1.5D0) THEN
        adjfac_var_1963 = 0.5D0 + (ratn2o_var_1965 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1964 = adjfac_var_1963 * chi_mls(4, jp_var_1914(jl_var_2001, lay_var_1947) + 1) * coldry_var_1921(jl_var_2001, lay_var_1947) * 1D-20
      ELSE
        adjcoln2o_var_1964 = coln2o_var_1920(jl_var_2001, lay_var_1947)
      END IF
      speccomb_planck_var_1936 = colh2o_var_1918(jl_var_2001, lay_var_1947) + refrat_planck_a_var_1937 * colco2_var_1919(jl_var_2001, lay_var_1947)
      specparm_planck_var_1962 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_planck_var_1936, oneminus_var_1917)
      specmult_planck_var_1961 = 8.0D0 * specparm_planck_var_1962
      jpl_var_1950 = 1 + INT(specmult_planck_var_1961)
      fpl_var_1960 = ((specmult_planck_var_1961) - AINT((specmult_planck_var_1961)))
      ind0_var_1940 = ((jp_var_1914(jl_var_2001, lay_var_1947) - 1) * 5 + (jt_var_1915(jl_var_2001, lay_var_1947) - 1)) * nspa_var_237(3) + js_var_1946
      ind1_var_1941 = (jp_var_1914(jl_var_2001, lay_var_1947) * 5 + (jt1_var_1916(jl_var_2001, lay_var_1947) - 1)) * nspa_var_237(3) + js1_var_1948
      inds_var_1942 = indself_var_1925(jl_var_2001, lay_var_1947)
      indf_var_1943 = indfor_var_1929(jl_var_2001, lay_var_1947)
      indm_var_1944 = indminor_var_1932(jl_var_2001, lay_var_1947)
      IF (specparm_var_1953 .LT. 0.125D0) THEN
        p_var_1979 = fs_var_1951 - 1.0D0
        p4_var_1980 = p_var_1979 ** 4
        fk0_var_1981 = p4_var_1980
        fk1_var_1982 = 1.0D0 - p_var_1979 - 2.0D0 * p4_var_1980
        fk2_var_1983 = p_var_1979 + p4_var_1980
        fac000_var_1967 = fk0_var_1981 * fac00_var_1909(jl_var_2001, lay_var_1947)
        fac100_var_1968 = fk1_var_1982 * fac00_var_1909(jl_var_2001, lay_var_1947)
        fac200_var_1969 = fk2_var_1983 * fac00_var_1909(jl_var_2001, lay_var_1947)
        fac010_var_1970 = fk0_var_1981 * fac10_var_1911(jl_var_2001, lay_var_1947)
        fac110_var_1971 = fk1_var_1982 * fac10_var_1911(jl_var_2001, lay_var_1947)
        fac210_var_1972 = fk2_var_1983 * fac10_var_1911(jl_var_2001, lay_var_1947)
      ELSE IF (specparm_var_1953 .GT. 0.875D0) THEN
        p_var_1979 = - fs_var_1951
        p4_var_1980 = p_var_1979 ** 4
        fk0_var_1981 = p4_var_1980
        fk1_var_1982 = 1.0D0 - p_var_1979 - 2.0D0 * p4_var_1980
        fk2_var_1983 = p_var_1979 + p4_var_1980
        fac000_var_1967 = fk0_var_1981 * fac00_var_1909(jl_var_2001, lay_var_1947)
        fac100_var_1968 = fk1_var_1982 * fac00_var_1909(jl_var_2001, lay_var_1947)
        fac200_var_1969 = fk2_var_1983 * fac00_var_1909(jl_var_2001, lay_var_1947)
        fac010_var_1970 = fk0_var_1981 * fac10_var_1911(jl_var_2001, lay_var_1947)
        fac110_var_1971 = fk1_var_1982 * fac10_var_1911(jl_var_2001, lay_var_1947)
        fac210_var_1972 = fk2_var_1983 * fac10_var_1911(jl_var_2001, lay_var_1947)
      ELSE
        fac000_var_1967 = (1.0D0 - fs_var_1951) * fac00_var_1909(jl_var_2001, lay_var_1947)
        fac010_var_1970 = (1.0D0 - fs_var_1951) * fac10_var_1911(jl_var_2001, lay_var_1947)
        fac100_var_1968 = fs_var_1951 * fac00_var_1909(jl_var_2001, lay_var_1947)
        fac110_var_1971 = fs_var_1951 * fac10_var_1911(jl_var_2001, lay_var_1947)
        fac200_var_1969 = 0.0D0
        fac210_var_1972 = 0.0D0
      END IF
      IF (specparm1_var_1956 .LT. 0.125D0) THEN
        p_var_1979 = fs1_var_1954 - 1.0D0
        p4_var_1980 = p_var_1979 ** 4
        fk0_var_1981 = p4_var_1980
        fk1_var_1982 = 1.0D0 - p_var_1979 - 2.0D0 * p4_var_1980
        fk2_var_1983 = p_var_1979 + p4_var_1980
        fac001_var_1973 = fk0_var_1981 * fac01_var_1910(jl_var_2001, lay_var_1947)
        fac101_var_1974 = fk1_var_1982 * fac01_var_1910(jl_var_2001, lay_var_1947)
        fac201_var_1975 = fk2_var_1983 * fac01_var_1910(jl_var_2001, lay_var_1947)
        fac011_var_1976 = fk0_var_1981 * fac11_var_1912(jl_var_2001, lay_var_1947)
        fac111_var_1977 = fk1_var_1982 * fac11_var_1912(jl_var_2001, lay_var_1947)
        fac211_var_1978 = fk2_var_1983 * fac11_var_1912(jl_var_2001, lay_var_1947)
      ELSE IF (specparm1_var_1956 .GT. 0.875D0) THEN
        p_var_1979 = - fs1_var_1954
        p4_var_1980 = p_var_1979 ** 4
        fk0_var_1981 = p4_var_1980
        fk1_var_1982 = 1.0D0 - p_var_1979 - 2.0D0 * p4_var_1980
        fk2_var_1983 = p_var_1979 + p4_var_1980
        fac001_var_1973 = fk0_var_1981 * fac01_var_1910(jl_var_2001, lay_var_1947)
        fac101_var_1974 = fk1_var_1982 * fac01_var_1910(jl_var_2001, lay_var_1947)
        fac201_var_1975 = fk2_var_1983 * fac01_var_1910(jl_var_2001, lay_var_1947)
        fac011_var_1976 = fk0_var_1981 * fac11_var_1912(jl_var_2001, lay_var_1947)
        fac111_var_1977 = fk1_var_1982 * fac11_var_1912(jl_var_2001, lay_var_1947)
        fac211_var_1978 = fk2_var_1983 * fac11_var_1912(jl_var_2001, lay_var_1947)
      ELSE
        fac001_var_1973 = (1.0D0 - fs1_var_1954) * fac01_var_1910(jl_var_2001, lay_var_1947)
        fac011_var_1976 = (1.0D0 - fs1_var_1954) * fac11_var_1912(jl_var_2001, lay_var_1947)
        fac101_var_1974 = fs1_var_1954 * fac01_var_1910(jl_var_2001, lay_var_1947)
        fac111_var_1977 = fs1_var_1954 * fac11_var_1912(jl_var_2001, lay_var_1947)
        fac201_var_1975 = 0.0D0
        fac211_var_1978 = 0.0D0
      END IF
      IF (specparm_var_1953 .LT. 0.125D0) THEN
        tau_major_var_1989(1 : ng3) = speccomb_var_1933 * (fac000_var_1967 * absa_var_184(ind0_var_1940, 1 : 16) + fac100_var_1968 * absa_var_184(ind0_var_1940 + 1, 1 : 16) + fac200_var_1969 * absa_var_184(ind0_var_1940 + 2, 1 : 16) + fac010_var_1970 * absa_var_184(ind0_var_1940 + 9, 1 : 16) + fac110_var_1971 * absa_var_184(ind0_var_1940 + 10, 1 : 16) + fac210_var_1972 * absa_var_184(ind0_var_1940 + 11, 1 : 16))
      ELSE IF (specparm_var_1953 .GT. 0.875D0) THEN
        tau_major_var_1989(1 : ng3) = speccomb_var_1933 * (fac200_var_1969 * absa_var_184(ind0_var_1940 - 1, 1 : 16) + fac100_var_1968 * absa_var_184(ind0_var_1940, 1 : 16) + fac000_var_1967 * absa_var_184(ind0_var_1940 + 1, 1 : 16) + fac210_var_1972 * absa_var_184(ind0_var_1940 + 8, 1 : 16) + fac110_var_1971 * absa_var_184(ind0_var_1940 + 9, 1 : 16) + fac010_var_1970 * absa_var_184(ind0_var_1940 + 10, 1 : 16))
      ELSE
        tau_major_var_1989(1 : ng3) = speccomb_var_1933 * (fac000_var_1967 * absa_var_184(ind0_var_1940, 1 : 16) + fac100_var_1968 * absa_var_184(ind0_var_1940 + 1, 1 : 16) + fac010_var_1970 * absa_var_184(ind0_var_1940 + 9, 1 : 16) + fac110_var_1971 * absa_var_184(ind0_var_1940 + 10, 1 : 16))
      END IF
      IF (specparm1_var_1956 .LT. 0.125D0) THEN
        tau_major1_var_1990(1 : ng3) = speccomb1_var_1934 * (fac001_var_1973 * absa_var_184(ind1_var_1941, 1 : 16) + fac101_var_1974 * absa_var_184(ind1_var_1941 + 1, 1 : 16) + fac201_var_1975 * absa_var_184(ind1_var_1941 + 2, 1 : 16) + fac011_var_1976 * absa_var_184(ind1_var_1941 + 9, 1 : 16) + fac111_var_1977 * absa_var_184(ind1_var_1941 + 10, 1 : 16) + fac211_var_1978 * absa_var_184(ind1_var_1941 + 11, 1 : 16))
      ELSE IF (specparm1_var_1956 .GT. 0.875D0) THEN
        tau_major1_var_1990(1 : ng3) = speccomb1_var_1934 * (fac201_var_1975 * absa_var_184(ind1_var_1941 - 1, 1 : 16) + fac101_var_1974 * absa_var_184(ind1_var_1941, 1 : 16) + fac001_var_1973 * absa_var_184(ind1_var_1941 + 1, 1 : 16) + fac211_var_1978 * absa_var_184(ind1_var_1941 + 8, 1 : 16) + fac111_var_1977 * absa_var_184(ind1_var_1941 + 9, 1 : 16) + fac011_var_1976 * absa_var_184(ind1_var_1941 + 10, 1 : 16))
      ELSE
        tau_major1_var_1990(1 : ng3) = speccomb1_var_1934 * (fac001_var_1973 * absa_var_184(ind1_var_1941, 1 : 16) + fac101_var_1974 * absa_var_184(ind1_var_1941 + 1, 1 : 16) + fac011_var_1976 * absa_var_184(ind1_var_1941 + 9, 1 : 16) + fac111_var_1977 * absa_var_184(ind1_var_1941 + 10, 1 : 16))
      END IF
      DO ig_var_1945 = 1, 16
        tauself_var_1985 = selffac_var_1923(jl_var_2001, lay_var_1947) * (selfref_var_186(inds_var_1942, ig_var_1945) + selffrac_var_1924(jl_var_2001, lay_var_1947) * (selfref_var_186(inds_var_1942 + 1, ig_var_1945) - selfref_var_186(inds_var_1942, ig_var_1945)))
        taufor_var_1984 = forfac_var_1913(jl_var_2001, lay_var_1947) * (forref_var_187(indf_var_1943, ig_var_1945) + forfrac_var_1930(jl_var_2001, lay_var_1947) * (forref_var_187(indf_var_1943 + 1, ig_var_1945) - forref_var_187(indf_var_1943, ig_var_1945)))
        n2om1_var_1986 = ka_mn2o_var_182(jmn2o_var_1949, indm_var_1944, ig_var_1945) + fmn2o_var_1957 * (ka_mn2o_var_182(jmn2o_var_1949 + 1, indm_var_1944, ig_var_1945) - ka_mn2o_var_182(jmn2o_var_1949, indm_var_1944, ig_var_1945))
        n2om2_var_1987 = ka_mn2o_var_182(jmn2o_var_1949, indm_var_1944 + 1, ig_var_1945) + fmn2o_var_1957 * (ka_mn2o_var_182(jmn2o_var_1949 + 1, indm_var_1944 + 1, ig_var_1945) - ka_mn2o_var_182(jmn2o_var_1949, indm_var_1944 + 1, ig_var_1945))
        absn2o_var_1988 = n2om1_var_1986 + minorfrac_var_1931(jl_var_2001, lay_var_1947) * (n2om2_var_1987 - n2om1_var_1986)
        taug_var_1907(jl_var_2001, 22 + ig_var_1945, lay_var_1947) = tau_major_var_1989(ig_var_1945) + tau_major1_var_1990(ig_var_1945) + tauself_var_1985 + taufor_var_1984 + adjcoln2o_var_1964 * absn2o_var_1988
        fracs_var_1926(jl_var_2001, 22 + ig_var_1945, lay_var_1947) = fracrefa_var_180(ig_var_1945, jpl_var_1950) + fpl_var_1960 * (fracrefa_var_180(ig_var_1945, jpl_var_1950 + 1) - fracrefa_var_180(ig_var_1945, jpl_var_1950))
      END DO
    END DO
  END DO
  DO lay_var_1947 = laytrop_max_var_1992 + 1, klev_var_1906
    DO jl_var_2001 = kidia_var_1904, kfdia_var_1905
      speccomb_var_1933 = colh2o_var_1918(jl_var_2001, lay_var_1947) + rat_h2oco2_var_1927(jl_var_2001, lay_var_1947) * colco2_var_1919(jl_var_2001, lay_var_1947)
      specparm_var_1953 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_var_1933, oneminus_var_1917)
      specmult_var_1952 = 4.0D0 * (specparm_var_1953)
      js_var_1946 = 1 + INT(specmult_var_1952)
      fs_var_1951 = ((specmult_var_1952) - AINT((specmult_var_1952)))
      speccomb1_var_1934 = colh2o_var_1918(jl_var_2001, lay_var_1947) + rat_h2oco2_1_var_1928(jl_var_2001, lay_var_1947) * colco2_var_1919(jl_var_2001, lay_var_1947)
      specparm1_var_1956 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb1_var_1934, oneminus_var_1917)
      specmult1_var_1955 = 4.0D0 * (specparm1_var_1956)
      js1_var_1948 = 1 + INT(specmult1_var_1955)
      fs1_var_1954 = ((specmult1_var_1955) - AINT((specmult1_var_1955)))
      fac000_var_1967 = (1.0D0 - fs_var_1951) * fac00_var_1909(jl_var_2001, lay_var_1947)
      fac010_var_1970 = (1.0D0 - fs_var_1951) * fac10_var_1911(jl_var_2001, lay_var_1947)
      fac100_var_1968 = fs_var_1951 * fac00_var_1909(jl_var_2001, lay_var_1947)
      fac110_var_1971 = fs_var_1951 * fac10_var_1911(jl_var_2001, lay_var_1947)
      fac001_var_1973 = (1.0D0 - fs1_var_1954) * fac01_var_1910(jl_var_2001, lay_var_1947)
      fac011_var_1976 = (1.0D0 - fs1_var_1954) * fac11_var_1912(jl_var_2001, lay_var_1947)
      fac101_var_1974 = fs1_var_1954 * fac01_var_1910(jl_var_2001, lay_var_1947)
      fac111_var_1977 = fs1_var_1954 * fac11_var_1912(jl_var_2001, lay_var_1947)
      speccomb_mn2o_var_1935 = colh2o_var_1918(jl_var_2001, lay_var_1947) + refrat_m_b * colco2_var_1919(jl_var_2001, lay_var_1947)
      specparm_mn2o_var_1959 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_mn2o_var_1935, oneminus_var_1917)
      specmult_mn2o_var_1958 = 4.0D0 * specparm_mn2o_var_1959
      jmn2o_var_1949 = 1 + INT(specmult_mn2o_var_1958)
      fmn2o_var_1957 = ((specmult_mn2o_var_1958) - AINT((specmult_mn2o_var_1958)))
      chi_n2o_var_1966 = coln2o_var_1920(jl_var_2001, lay_var_1947) / coldry_var_1921(jl_var_2001, lay_var_1947)
      ratn2o_var_1965 = 1D+20 * chi_n2o_var_1966 / chi_mls(4, jp_var_1914(jl_var_2001, lay_var_1947) + 1)
      IF (ratn2o_var_1965 .GT. 1.5D0) THEN
        adjfac_var_1963 = 0.5D0 + (ratn2o_var_1965 - 0.5D0) ** 0.65D0
        adjcoln2o_var_1964 = adjfac_var_1963 * chi_mls(4, jp_var_1914(jl_var_2001, lay_var_1947) + 1) * coldry_var_1921(jl_var_2001, lay_var_1947) * 1D-20
      ELSE
        adjcoln2o_var_1964 = coln2o_var_1920(jl_var_2001, lay_var_1947)
      END IF
      speccomb_planck_var_1936 = colh2o_var_1918(jl_var_2001, lay_var_1947) + refrat_planck_b_var_1938 * colco2_var_1919(jl_var_2001, lay_var_1947)
      specparm_planck_var_1962 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_planck_var_1936, oneminus_var_1917)
      specmult_planck_var_1961 = 4.0D0 * specparm_planck_var_1962
      jpl_var_1950 = 1 + INT(specmult_planck_var_1961)
      fpl_var_1960 = ((specmult_planck_var_1961) - AINT((specmult_planck_var_1961)))
      ind0_var_1940 = ((jp_var_1914(jl_var_2001, lay_var_1947) - 13) * 5 + (jt_var_1915(jl_var_2001, lay_var_1947) - 1)) * nspb_var_238(3) + js_var_1946
      ind1_var_1941 = ((jp_var_1914(jl_var_2001, lay_var_1947) - 12) * 5 + (jt1_var_1916(jl_var_2001, lay_var_1947) - 1)) * nspb_var_238(3) + js1_var_1948
      indf_var_1943 = indfor_var_1929(jl_var_2001, lay_var_1947)
      indm_var_1944 = indminor_var_1932(jl_var_2001, lay_var_1947)
      DO ig_var_1945 = 1, 16
        taufor_var_1984 = forfac_var_1913(jl_var_2001, lay_var_1947) * (forref_var_187(indf_var_1943, ig_var_1945) + forfrac_var_1930(jl_var_2001, lay_var_1947) * (forref_var_187(indf_var_1943 + 1, ig_var_1945) - forref_var_187(indf_var_1943, ig_var_1945)))
        n2om1_var_1986 = kb_mn2o_var_183(jmn2o_var_1949, indm_var_1944, ig_var_1945) + fmn2o_var_1957 * (kb_mn2o_var_183(jmn2o_var_1949 + 1, indm_var_1944, ig_var_1945) - kb_mn2o_var_183(jmn2o_var_1949, indm_var_1944, ig_var_1945))
        n2om2_var_1987 = kb_mn2o_var_183(jmn2o_var_1949, indm_var_1944 + 1, ig_var_1945) + fmn2o_var_1957 * (kb_mn2o_var_183(jmn2o_var_1949 + 1, indm_var_1944 + 1, ig_var_1945) - kb_mn2o_var_183(jmn2o_var_1949, indm_var_1944 + 1, ig_var_1945))
        absn2o_var_1988 = n2om1_var_1986 + minorfrac_var_1931(jl_var_2001, lay_var_1947) * (n2om2_var_1987 - n2om1_var_1986)
        taug_var_1907(jl_var_2001, 22 + ig_var_1945, lay_var_1947) = speccomb_var_1933 * (fac000_var_1967 * absb_var_185(ind0_var_1940, ig_var_1945) + fac100_var_1968 * absb_var_185(ind0_var_1940 + 1, ig_var_1945) + fac010_var_1970 * absb_var_185(ind0_var_1940 + 5, ig_var_1945) + fac110_var_1971 * absb_var_185(ind0_var_1940 + 6, ig_var_1945)) + speccomb1_var_1934 * (fac001_var_1973 * absb_var_185(ind1_var_1941, ig_var_1945) + fac101_var_1974 * absb_var_185(ind1_var_1941 + 1, ig_var_1945) + fac011_var_1976 * absb_var_185(ind1_var_1941 + 5, ig_var_1945) + fac111_var_1977 * absb_var_185(ind1_var_1941 + 6, ig_var_1945)) + taufor_var_1984 + adjcoln2o_var_1964 * absn2o_var_1988
        fracs_var_1926(jl_var_2001, 22 + ig_var_1945, lay_var_1947) = fracrefb_var_181(ig_var_1945, jpl_var_1950) + fpl_var_1960 * (fracrefb_var_181(ig_var_1945, jpl_var_1950 + 1) - fracrefb_var_181(ig_var_1945, jpl_var_1950))
      END DO
    END DO
  END DO
  IF (laytrop_max_var_1992 /= laytrop_min_var_1991) THEN
    DO lay_var_1947 = laytrop_min_var_1991 + 1, laytrop_max_var_1992
      ixc0_var_1998 = ixc_var_1993(lay_var_1947)
      DO ixp_var_1999 = 1, ixc0_var_1998
        jl_var_2001 = ixlow_var_1994(ixp_var_1999, lay_var_1947)
        speccomb_var_1933 = colh2o_var_1918(jl_var_2001, lay_var_1947) + rat_h2oco2_var_1927(jl_var_2001, lay_var_1947) * colco2_var_1919(jl_var_2001, lay_var_1947)
        specparm_var_1953 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_var_1933, oneminus_var_1917)
        specmult_var_1952 = 8.0D0 * (specparm_var_1953)
        js_var_1946 = 1 + INT(specmult_var_1952)
        fs_var_1951 = ((specmult_var_1952) - AINT((specmult_var_1952)))
        speccomb1_var_1934 = colh2o_var_1918(jl_var_2001, lay_var_1947) + rat_h2oco2_1_var_1928(jl_var_2001, lay_var_1947) * colco2_var_1919(jl_var_2001, lay_var_1947)
        specparm1_var_1956 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb1_var_1934, oneminus_var_1917)
        specmult1_var_1955 = 8.0D0 * (specparm1_var_1956)
        js1_var_1948 = 1 + INT(specmult1_var_1955)
        fs1_var_1954 = ((specmult1_var_1955) - AINT((specmult1_var_1955)))
        speccomb_mn2o_var_1935 = colh2o_var_1918(jl_var_2001, lay_var_1947) + refrat_m_a_var_1939 * colco2_var_1919(jl_var_2001, lay_var_1947)
        specparm_mn2o_var_1959 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_mn2o_var_1935, oneminus_var_1917)
        specmult_mn2o_var_1958 = 8.0D0 * specparm_mn2o_var_1959
        jmn2o_var_1949 = 1 + INT(specmult_mn2o_var_1958)
        fmn2o_var_1957 = ((specmult_mn2o_var_1958) - AINT((specmult_mn2o_var_1958)))
        chi_n2o_var_1966 = coln2o_var_1920(jl_var_2001, lay_var_1947) / coldry_var_1921(jl_var_2001, lay_var_1947)
        ratn2o_var_1965 = 1D+20 * chi_n2o_var_1966 / chi_mls(4, jp_var_1914(jl_var_2001, lay_var_1947) + 1)
        IF (ratn2o_var_1965 .GT. 1.5D0) THEN
          adjfac_var_1963 = 0.5D0 + (ratn2o_var_1965 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1964 = adjfac_var_1963 * chi_mls(4, jp_var_1914(jl_var_2001, lay_var_1947) + 1) * coldry_var_1921(jl_var_2001, lay_var_1947) * 1D-20
        ELSE
          adjcoln2o_var_1964 = coln2o_var_1920(jl_var_2001, lay_var_1947)
        END IF
        speccomb_planck_var_1936 = colh2o_var_1918(jl_var_2001, lay_var_1947) + refrat_planck_a_var_1937 * colco2_var_1919(jl_var_2001, lay_var_1947)
        specparm_planck_var_1962 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_planck_var_1936, oneminus_var_1917)
        specmult_planck_var_1961 = 8.0D0 * specparm_planck_var_1962
        jpl_var_1950 = 1 + INT(specmult_planck_var_1961)
        fpl_var_1960 = ((specmult_planck_var_1961) - AINT((specmult_planck_var_1961)))
        ind0_var_1940 = ((jp_var_1914(jl_var_2001, lay_var_1947) - 1) * 5 + (jt_var_1915(jl_var_2001, lay_var_1947) - 1)) * nspa_var_237(3) + js_var_1946
        ind1_var_1941 = (jp_var_1914(jl_var_2001, lay_var_1947) * 5 + (jt1_var_1916(jl_var_2001, lay_var_1947) - 1)) * nspa_var_237(3) + js1_var_1948
        inds_var_1942 = indself_var_1925(jl_var_2001, lay_var_1947)
        indf_var_1943 = indfor_var_1929(jl_var_2001, lay_var_1947)
        indm_var_1944 = indminor_var_1932(jl_var_2001, lay_var_1947)
        IF (specparm_var_1953 .LT. 0.125D0) THEN
          p_var_1979 = fs_var_1951 - 1.0D0
          p4_var_1980 = p_var_1979 ** 4
          fk0_var_1981 = p4_var_1980
          fk1_var_1982 = 1.0D0 - p_var_1979 - 2.0D0 * p4_var_1980
          fk2_var_1983 = p_var_1979 + p4_var_1980
          fac000_var_1967 = fk0_var_1981 * fac00_var_1909(jl_var_2001, lay_var_1947)
          fac100_var_1968 = fk1_var_1982 * fac00_var_1909(jl_var_2001, lay_var_1947)
          fac200_var_1969 = fk2_var_1983 * fac00_var_1909(jl_var_2001, lay_var_1947)
          fac010_var_1970 = fk0_var_1981 * fac10_var_1911(jl_var_2001, lay_var_1947)
          fac110_var_1971 = fk1_var_1982 * fac10_var_1911(jl_var_2001, lay_var_1947)
          fac210_var_1972 = fk2_var_1983 * fac10_var_1911(jl_var_2001, lay_var_1947)
        ELSE IF (specparm_var_1953 .GT. 0.875D0) THEN
          p_var_1979 = - fs_var_1951
          p4_var_1980 = p_var_1979 ** 4
          fk0_var_1981 = p4_var_1980
          fk1_var_1982 = 1.0D0 - p_var_1979 - 2.0D0 * p4_var_1980
          fk2_var_1983 = p_var_1979 + p4_var_1980
          fac000_var_1967 = fk0_var_1981 * fac00_var_1909(jl_var_2001, lay_var_1947)
          fac100_var_1968 = fk1_var_1982 * fac00_var_1909(jl_var_2001, lay_var_1947)
          fac200_var_1969 = fk2_var_1983 * fac00_var_1909(jl_var_2001, lay_var_1947)
          fac010_var_1970 = fk0_var_1981 * fac10_var_1911(jl_var_2001, lay_var_1947)
          fac110_var_1971 = fk1_var_1982 * fac10_var_1911(jl_var_2001, lay_var_1947)
          fac210_var_1972 = fk2_var_1983 * fac10_var_1911(jl_var_2001, lay_var_1947)
        ELSE
          fac000_var_1967 = (1.0D0 - fs_var_1951) * fac00_var_1909(jl_var_2001, lay_var_1947)
          fac010_var_1970 = (1.0D0 - fs_var_1951) * fac10_var_1911(jl_var_2001, lay_var_1947)
          fac100_var_1968 = fs_var_1951 * fac00_var_1909(jl_var_2001, lay_var_1947)
          fac110_var_1971 = fs_var_1951 * fac10_var_1911(jl_var_2001, lay_var_1947)
          fac200_var_1969 = 0.0D0
          fac210_var_1972 = 0.0D0
        END IF
        IF (specparm1_var_1956 .LT. 0.125D0) THEN
          p_var_1979 = fs1_var_1954 - 1.0D0
          p4_var_1980 = p_var_1979 ** 4
          fk0_var_1981 = p4_var_1980
          fk1_var_1982 = 1.0D0 - p_var_1979 - 2.0D0 * p4_var_1980
          fk2_var_1983 = p_var_1979 + p4_var_1980
          fac001_var_1973 = fk0_var_1981 * fac01_var_1910(jl_var_2001, lay_var_1947)
          fac101_var_1974 = fk1_var_1982 * fac01_var_1910(jl_var_2001, lay_var_1947)
          fac201_var_1975 = fk2_var_1983 * fac01_var_1910(jl_var_2001, lay_var_1947)
          fac011_var_1976 = fk0_var_1981 * fac11_var_1912(jl_var_2001, lay_var_1947)
          fac111_var_1977 = fk1_var_1982 * fac11_var_1912(jl_var_2001, lay_var_1947)
          fac211_var_1978 = fk2_var_1983 * fac11_var_1912(jl_var_2001, lay_var_1947)
        ELSE IF (specparm1_var_1956 .GT. 0.875D0) THEN
          p_var_1979 = - fs1_var_1954
          p4_var_1980 = p_var_1979 ** 4
          fk0_var_1981 = p4_var_1980
          fk1_var_1982 = 1.0D0 - p_var_1979 - 2.0D0 * p4_var_1980
          fk2_var_1983 = p_var_1979 + p4_var_1980
          fac001_var_1973 = fk0_var_1981 * fac01_var_1910(jl_var_2001, lay_var_1947)
          fac101_var_1974 = fk1_var_1982 * fac01_var_1910(jl_var_2001, lay_var_1947)
          fac201_var_1975 = fk2_var_1983 * fac01_var_1910(jl_var_2001, lay_var_1947)
          fac011_var_1976 = fk0_var_1981 * fac11_var_1912(jl_var_2001, lay_var_1947)
          fac111_var_1977 = fk1_var_1982 * fac11_var_1912(jl_var_2001, lay_var_1947)
          fac211_var_1978 = fk2_var_1983 * fac11_var_1912(jl_var_2001, lay_var_1947)
        ELSE
          fac001_var_1973 = (1.0D0 - fs1_var_1954) * fac01_var_1910(jl_var_2001, lay_var_1947)
          fac011_var_1976 = (1.0D0 - fs1_var_1954) * fac11_var_1912(jl_var_2001, lay_var_1947)
          fac101_var_1974 = fs1_var_1954 * fac01_var_1910(jl_var_2001, lay_var_1947)
          fac111_var_1977 = fs1_var_1954 * fac11_var_1912(jl_var_2001, lay_var_1947)
          fac201_var_1975 = 0.0D0
          fac211_var_1978 = 0.0D0
        END IF
        IF (specparm_var_1953 .LT. 0.125D0) THEN
          tau_major_var_1989(1 : ng3) = speccomb_var_1933 * (fac000_var_1967 * absa_var_184(ind0_var_1940, 1 : 16) + fac100_var_1968 * absa_var_184(ind0_var_1940 + 1, 1 : 16) + fac200_var_1969 * absa_var_184(ind0_var_1940 + 2, 1 : 16) + fac010_var_1970 * absa_var_184(ind0_var_1940 + 9, 1 : 16) + fac110_var_1971 * absa_var_184(ind0_var_1940 + 10, 1 : 16) + fac210_var_1972 * absa_var_184(ind0_var_1940 + 11, 1 : 16))
        ELSE IF (specparm_var_1953 .GT. 0.875D0) THEN
          tau_major_var_1989(1 : ng3) = speccomb_var_1933 * (fac200_var_1969 * absa_var_184(ind0_var_1940 - 1, 1 : 16) + fac100_var_1968 * absa_var_184(ind0_var_1940, 1 : 16) + fac000_var_1967 * absa_var_184(ind0_var_1940 + 1, 1 : 16) + fac210_var_1972 * absa_var_184(ind0_var_1940 + 8, 1 : 16) + fac110_var_1971 * absa_var_184(ind0_var_1940 + 9, 1 : 16) + fac010_var_1970 * absa_var_184(ind0_var_1940 + 10, 1 : 16))
        ELSE
          tau_major_var_1989(1 : ng3) = speccomb_var_1933 * (fac000_var_1967 * absa_var_184(ind0_var_1940, 1 : 16) + fac100_var_1968 * absa_var_184(ind0_var_1940 + 1, 1 : 16) + fac010_var_1970 * absa_var_184(ind0_var_1940 + 9, 1 : 16) + fac110_var_1971 * absa_var_184(ind0_var_1940 + 10, 1 : 16))
        END IF
        IF (specparm1_var_1956 .LT. 0.125D0) THEN
          tau_major1_var_1990(1 : ng3) = speccomb1_var_1934 * (fac001_var_1973 * absa_var_184(ind1_var_1941, 1 : 16) + fac101_var_1974 * absa_var_184(ind1_var_1941 + 1, 1 : 16) + fac201_var_1975 * absa_var_184(ind1_var_1941 + 2, 1 : 16) + fac011_var_1976 * absa_var_184(ind1_var_1941 + 9, 1 : 16) + fac111_var_1977 * absa_var_184(ind1_var_1941 + 10, 1 : 16) + fac211_var_1978 * absa_var_184(ind1_var_1941 + 11, 1 : 16))
        ELSE IF (specparm1_var_1956 .GT. 0.875D0) THEN
          tau_major1_var_1990(1 : ng3) = speccomb1_var_1934 * (fac201_var_1975 * absa_var_184(ind1_var_1941 - 1, 1 : 16) + fac101_var_1974 * absa_var_184(ind1_var_1941, 1 : 16) + fac001_var_1973 * absa_var_184(ind1_var_1941 + 1, 1 : 16) + fac211_var_1978 * absa_var_184(ind1_var_1941 + 8, 1 : 16) + fac111_var_1977 * absa_var_184(ind1_var_1941 + 9, 1 : 16) + fac011_var_1976 * absa_var_184(ind1_var_1941 + 10, 1 : 16))
        ELSE
          tau_major1_var_1990(1 : ng3) = speccomb1_var_1934 * (fac001_var_1973 * absa_var_184(ind1_var_1941, 1 : 16) + fac101_var_1974 * absa_var_184(ind1_var_1941 + 1, 1 : 16) + fac011_var_1976 * absa_var_184(ind1_var_1941 + 9, 1 : 16) + fac111_var_1977 * absa_var_184(ind1_var_1941 + 10, 1 : 16))
        END IF
        DO ig_var_1945 = 1, 16
          tauself_var_1985 = selffac_var_1923(jl_var_2001, lay_var_1947) * (selfref_var_186(inds_var_1942, ig_var_1945) + selffrac_var_1924(jl_var_2001, lay_var_1947) * (selfref_var_186(inds_var_1942 + 1, ig_var_1945) - selfref_var_186(inds_var_1942, ig_var_1945)))
          taufor_var_1984 = forfac_var_1913(jl_var_2001, lay_var_1947) * (forref_var_187(indf_var_1943, ig_var_1945) + forfrac_var_1930(jl_var_2001, lay_var_1947) * (forref_var_187(indf_var_1943 + 1, ig_var_1945) - forref_var_187(indf_var_1943, ig_var_1945)))
          n2om1_var_1986 = ka_mn2o_var_182(jmn2o_var_1949, indm_var_1944, ig_var_1945) + fmn2o_var_1957 * (ka_mn2o_var_182(jmn2o_var_1949 + 1, indm_var_1944, ig_var_1945) - ka_mn2o_var_182(jmn2o_var_1949, indm_var_1944, ig_var_1945))
          n2om2_var_1987 = ka_mn2o_var_182(jmn2o_var_1949, indm_var_1944 + 1, ig_var_1945) + fmn2o_var_1957 * (ka_mn2o_var_182(jmn2o_var_1949 + 1, indm_var_1944 + 1, ig_var_1945) - ka_mn2o_var_182(jmn2o_var_1949, indm_var_1944 + 1, ig_var_1945))
          absn2o_var_1988 = n2om1_var_1986 + minorfrac_var_1931(jl_var_2001, lay_var_1947) * (n2om2_var_1987 - n2om1_var_1986)
          taug_var_1907(jl_var_2001, 22 + ig_var_1945, lay_var_1947) = tau_major_var_1989(ig_var_1945) + tau_major1_var_1990(ig_var_1945) + tauself_var_1985 + taufor_var_1984 + adjcoln2o_var_1964 * absn2o_var_1988
          fracs_var_1926(jl_var_2001, 22 + ig_var_1945, lay_var_1947) = fracrefa_var_180(ig_var_1945, jpl_var_1950) + fpl_var_1960 * (fracrefa_var_180(ig_var_1945, jpl_var_1950 + 1) - fracrefa_var_180(ig_var_1945, jpl_var_1950))
        END DO
      END DO
      ixc0_var_1998 = kfdia_var_1905 - kidia_var_1904 + 1 - ixc0_var_1998
      DO ixp_var_1999 = 1, ixc0_var_1998
        jl_var_2001 = ixhigh_var_1995(ixp_var_1999, lay_var_1947)
        speccomb_var_1933 = colh2o_var_1918(jl_var_2001, lay_var_1947) + rat_h2oco2_var_1927(jl_var_2001, lay_var_1947) * colco2_var_1919(jl_var_2001, lay_var_1947)
        specparm_var_1953 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_var_1933, oneminus_var_1917)
        specmult_var_1952 = 4.0D0 * (specparm_var_1953)
        js_var_1946 = 1 + INT(specmult_var_1952)
        fs_var_1951 = ((specmult_var_1952) - AINT((specmult_var_1952)))
        speccomb1_var_1934 = colh2o_var_1918(jl_var_2001, lay_var_1947) + rat_h2oco2_1_var_1928(jl_var_2001, lay_var_1947) * colco2_var_1919(jl_var_2001, lay_var_1947)
        specparm1_var_1956 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb1_var_1934, oneminus_var_1917)
        specmult1_var_1955 = 4.0D0 * (specparm1_var_1956)
        js1_var_1948 = 1 + INT(specmult1_var_1955)
        fs1_var_1954 = ((specmult1_var_1955) - AINT((specmult1_var_1955)))
        fac000_var_1967 = (1.0D0 - fs_var_1951) * fac00_var_1909(jl_var_2001, lay_var_1947)
        fac010_var_1970 = (1.0D0 - fs_var_1951) * fac10_var_1911(jl_var_2001, lay_var_1947)
        fac100_var_1968 = fs_var_1951 * fac00_var_1909(jl_var_2001, lay_var_1947)
        fac110_var_1971 = fs_var_1951 * fac10_var_1911(jl_var_2001, lay_var_1947)
        fac001_var_1973 = (1.0D0 - fs1_var_1954) * fac01_var_1910(jl_var_2001, lay_var_1947)
        fac011_var_1976 = (1.0D0 - fs1_var_1954) * fac11_var_1912(jl_var_2001, lay_var_1947)
        fac101_var_1974 = fs1_var_1954 * fac01_var_1910(jl_var_2001, lay_var_1947)
        fac111_var_1977 = fs1_var_1954 * fac11_var_1912(jl_var_2001, lay_var_1947)
        speccomb_mn2o_var_1935 = colh2o_var_1918(jl_var_2001, lay_var_1947) + refrat_m_b * colco2_var_1919(jl_var_2001, lay_var_1947)
        specparm_mn2o_var_1959 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_mn2o_var_1935, oneminus_var_1917)
        specmult_mn2o_var_1958 = 4.0D0 * specparm_mn2o_var_1959
        jmn2o_var_1949 = 1 + INT(specmult_mn2o_var_1958)
        fmn2o_var_1957 = ((specmult_mn2o_var_1958) - AINT((specmult_mn2o_var_1958)))
        chi_n2o_var_1966 = coln2o_var_1920(jl_var_2001, lay_var_1947) / coldry_var_1921(jl_var_2001, lay_var_1947)
        ratn2o_var_1965 = 1D+20 * chi_n2o_var_1966 / chi_mls(4, jp_var_1914(jl_var_2001, lay_var_1947) + 1)
        IF (ratn2o_var_1965 .GT. 1.5D0) THEN
          adjfac_var_1963 = 0.5D0 + (ratn2o_var_1965 - 0.5D0) ** 0.65D0
          adjcoln2o_var_1964 = adjfac_var_1963 * chi_mls(4, jp_var_1914(jl_var_2001, lay_var_1947) + 1) * coldry_var_1921(jl_var_2001, lay_var_1947) * 1D-20
        ELSE
          adjcoln2o_var_1964 = coln2o_var_1920(jl_var_2001, lay_var_1947)
        END IF
        speccomb_planck_var_1936 = colh2o_var_1918(jl_var_2001, lay_var_1947) + refrat_planck_b_var_1938 * colco2_var_1919(jl_var_2001, lay_var_1947)
        specparm_planck_var_1962 = MIN(colh2o_var_1918(jl_var_2001, lay_var_1947) / speccomb_planck_var_1936, oneminus_var_1917)
        specmult_planck_var_1961 = 4.0D0 * specparm_planck_var_1962
        jpl_var_1950 = 1 + INT(specmult_planck_var_1961)
        fpl_var_1960 = ((specmult_planck_var_1961) - AINT((specmult_planck_var_1961)))
        ind0_var_1940 = ((jp_var_1914(jl_var_2001, lay_var_1947) - 13) * 5 + (jt_var_1915(jl_var_2001, lay_var_1947) - 1)) * nspb_var_238(3) + js_var_1946
        ind1_var_1941 = ((jp_var_1914(jl_var_2001, lay_var_1947) - 12) * 5 + (jt1_var_1916(jl_var_2001, lay_var_1947) - 1)) * nspb_var_238(3) + js1_var_1948
        indf_var_1943 = indfor_var_1929(jl_var_2001, lay_var_1947)
        indm_var_1944 = indminor_var_1932(jl_var_2001, lay_var_1947)
        DO ig_var_1945 = 1, 16
          taufor_var_1984 = forfac_var_1913(jl_var_2001, lay_var_1947) * (forref_var_187(indf_var_1943, ig_var_1945) + forfrac_var_1930(jl_var_2001, lay_var_1947) * (forref_var_187(indf_var_1943 + 1, ig_var_1945) - forref_var_187(indf_var_1943, ig_var_1945)))
          n2om1_var_1986 = kb_mn2o_var_183(jmn2o_var_1949, indm_var_1944, ig_var_1945) + fmn2o_var_1957 * (kb_mn2o_var_183(jmn2o_var_1949 + 1, indm_var_1944, ig_var_1945) - kb_mn2o_var_183(jmn2o_var_1949, indm_var_1944, ig_var_1945))
          n2om2_var_1987 = kb_mn2o_var_183(jmn2o_var_1949, indm_var_1944 + 1, ig_var_1945) + fmn2o_var_1957 * (kb_mn2o_var_183(jmn2o_var_1949 + 1, indm_var_1944 + 1, ig_var_1945) - kb_mn2o_var_183(jmn2o_var_1949, indm_var_1944 + 1, ig_var_1945))
          absn2o_var_1988 = n2om1_var_1986 + minorfrac_var_1931(jl_var_2001, lay_var_1947) * (n2om2_var_1987 - n2om1_var_1986)
          taug_var_1907(jl_var_2001, 22 + ig_var_1945, lay_var_1947) = speccomb_var_1933 * (fac000_var_1967 * absb_var_185(ind0_var_1940, ig_var_1945) + fac100_var_1968 * absb_var_185(ind0_var_1940 + 1, ig_var_1945) + fac010_var_1970 * absb_var_185(ind0_var_1940 + 5, ig_var_1945) + fac110_var_1971 * absb_var_185(ind0_var_1940 + 6, ig_var_1945)) + speccomb1_var_1934 * (fac001_var_1973 * absb_var_185(ind1_var_1941, ig_var_1945) + fac101_var_1974 * absb_var_185(ind1_var_1941 + 1, ig_var_1945) + fac011_var_1976 * absb_var_185(ind1_var_1941 + 5, ig_var_1945) + fac111_var_1977 * absb_var_185(ind1_var_1941 + 6, ig_var_1945)) + taufor_var_1984 + adjcoln2o_var_1964 * absn2o_var_1988
          fracs_var_1926(jl_var_2001, 22 + ig_var_1945, lay_var_1947) = fracrefb_var_181(ig_var_1945, jpl_var_1950) + fpl_var_1960 * (fracrefb_var_181(ig_var_1945, jpl_var_1950 + 1) - fracrefb_var_181(ig_var_1945, jpl_var_1950))
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol3
SUBROUTINE rrtm_taumol2(kidia_var_2002, kfdia_var_2003, klev_var_2004, taug_var_2005, pavel_var_2006, coldry_var_2007, p_tauaerl_var_2008, fac00_var_2009, fac01_var_2010, fac10_var_2011, fac11_var_2012, forfac_var_2014, forfrac_var_2013, indfor_var_2023, jp_var_2015, jt_var_2016, jt1_var_2017, colh2o_var_2018, laytrop_var_2019, selffac_var_2020, selffrac_var_2021, indself_var_2022, fracs_var_2024)
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrta2, ONLY: absa_var_176, absb_var_177, forref_var_179, fracrefa_var_174, fracrefb_var_175, selfref_var_178
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2002
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2003
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2004
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2005(kidia_var_2002 : kfdia_var_2003, 140, klev_var_2004)
  REAL(KIND = 8), INTENT(IN) :: pavel_var_2006(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_2007(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2008(kidia_var_2002 : kfdia_var_2003, klev_var_2004, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2009(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2010(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2011(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2012(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2013(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2014(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2015(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2016(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2017(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2018(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2019(kidia_var_2002 : kfdia_var_2003)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2020(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2021(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2022(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2023(kidia_var_2002 : kfdia_var_2003, klev_var_2004)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2024(kidia_var_2002 : kfdia_var_2003, 140, klev_var_2004)
  INTEGER(KIND = 4) :: ind0_var_2025, ind1_var_2026, inds_var_2027, indf_var_2028
  INTEGER(KIND = 4) :: ig_var_2029, lay_var_2030
  REAL(KIND = 8) :: taufor_var_2031, tauself_var_2032, corradj_var_2033, pp_var_2034
  INTEGER(KIND = 4) :: laytrop_min_var_2035, laytrop_max_var_2036
  INTEGER(KIND = 4) :: ixc_var_2037(klev_var_2004), ixlow_var_2038(kfdia_var_2003, klev_var_2004), ixhigh_var_2039(kfdia_var_2003, klev_var_2004)
  INTEGER(KIND = 4) :: ich_var_2040, icl_var_2041, ixc0_var_2042, ixp_var_2043, jc_var_2044, jl_var_2045
  laytrop_min_var_2035 = MINVAL(laytrop_var_2019)
  laytrop_max_var_2036 = MAXVAL(laytrop_var_2019)
  ixlow_var_2038 = 0
  ixhigh_var_2039 = 0
  ixc_var_2037 = 0
  DO lay_var_2030 = laytrop_min_var_2035 + 1, laytrop_max_var_2036
    icl_var_2041 = 0
    ich_var_2040 = 0
    DO jc_var_2044 = kidia_var_2002, kfdia_var_2003
      IF (lay_var_2030 <= laytrop_var_2019(jc_var_2044)) THEN
        icl_var_2041 = icl_var_2041 + 1
        ixlow_var_2038(icl_var_2041, lay_var_2030) = jc_var_2044
      ELSE
        ich_var_2040 = ich_var_2040 + 1
        ixhigh_var_2039(ich_var_2040, lay_var_2030) = jc_var_2044
      END IF
    END DO
    ixc_var_2037(lay_var_2030) = icl_var_2041
  END DO
  DO lay_var_2030 = 1, laytrop_min_var_2035
    DO jl_var_2045 = kidia_var_2002, kfdia_var_2003
      ind0_var_2025 = ((jp_var_2015(jl_var_2045, lay_var_2030) - 1) * 5 + (jt_var_2016(jl_var_2045, lay_var_2030) - 1)) * nspa_var_237(2) + 1
      ind1_var_2026 = (jp_var_2015(jl_var_2045, lay_var_2030) * 5 + (jt1_var_2017(jl_var_2045, lay_var_2030) - 1)) * nspa_var_237(2) + 1
      inds_var_2027 = indself_var_2022(jl_var_2045, lay_var_2030)
      indf_var_2028 = indfor_var_2023(jl_var_2045, lay_var_2030)
      pp_var_2034 = pavel_var_2006(jl_var_2045, lay_var_2030)
      corradj_var_2033 = 1.0D0 - 0.05D0 * (pp_var_2034 - 100.0D0) / 900.0D0
      DO ig_var_2029 = 1, 12
        tauself_var_2032 = selffac_var_2020(jl_var_2045, lay_var_2030) * (selfref_var_178(inds_var_2027, ig_var_2029) + selffrac_var_2021(jl_var_2045, lay_var_2030) * (selfref_var_178(inds_var_2027 + 1, ig_var_2029) - selfref_var_178(inds_var_2027, ig_var_2029)))
        taufor_var_2031 = forfac_var_2014(jl_var_2045, lay_var_2030) * (forref_var_179(indf_var_2028, ig_var_2029) + forfrac_var_2013(jl_var_2045, lay_var_2030) * (forref_var_179(indf_var_2028 + 1, ig_var_2029) - forref_var_179(indf_var_2028, ig_var_2029)))
        taug_var_2005(jl_var_2045, 10 + ig_var_2029, lay_var_2030) = corradj_var_2033 * (colh2o_var_2018(jl_var_2045, lay_var_2030) * (fac00_var_2009(jl_var_2045, lay_var_2030) * absa_var_176(ind0_var_2025, ig_var_2029) + fac10_var_2011(jl_var_2045, lay_var_2030) * absa_var_176(ind0_var_2025 + 1, ig_var_2029) + fac01_var_2010(jl_var_2045, lay_var_2030) * absa_var_176(ind1_var_2026, ig_var_2029) + fac11_var_2012(jl_var_2045, lay_var_2030) * absa_var_176(ind1_var_2026 + 1, ig_var_2029)) + tauself_var_2032 + taufor_var_2031)
        fracs_var_2024(jl_var_2045, 10 + ig_var_2029, lay_var_2030) = fracrefa_var_174(ig_var_2029)
      END DO
    END DO
  END DO
  DO lay_var_2030 = laytrop_max_var_2036 + 1, klev_var_2004
    DO jl_var_2045 = kidia_var_2002, kfdia_var_2003
      ind0_var_2025 = ((jp_var_2015(jl_var_2045, lay_var_2030) - 13) * 5 + (jt_var_2016(jl_var_2045, lay_var_2030) - 1)) * nspb_var_238(2) + 1
      ind1_var_2026 = ((jp_var_2015(jl_var_2045, lay_var_2030) - 12) * 5 + (jt1_var_2017(jl_var_2045, lay_var_2030) - 1)) * nspb_var_238(2) + 1
      indf_var_2028 = indfor_var_2023(jl_var_2045, lay_var_2030)
      DO ig_var_2029 = 1, 12
        taufor_var_2031 = forfac_var_2014(jl_var_2045, lay_var_2030) * (forref_var_179(indf_var_2028, ig_var_2029) + forfrac_var_2013(jl_var_2045, lay_var_2030) * (forref_var_179(indf_var_2028 + 1, ig_var_2029) - forref_var_179(indf_var_2028, ig_var_2029)))
        taug_var_2005(jl_var_2045, 10 + ig_var_2029, lay_var_2030) = colh2o_var_2018(jl_var_2045, lay_var_2030) * (fac00_var_2009(jl_var_2045, lay_var_2030) * absb_var_177(ind0_var_2025, ig_var_2029) + fac10_var_2011(jl_var_2045, lay_var_2030) * absb_var_177(ind0_var_2025 + 1, ig_var_2029) + fac01_var_2010(jl_var_2045, lay_var_2030) * absb_var_177(ind1_var_2026, ig_var_2029) + fac11_var_2012(jl_var_2045, lay_var_2030) * absb_var_177(ind1_var_2026 + 1, ig_var_2029)) + taufor_var_2031
        fracs_var_2024(jl_var_2045, 10 + ig_var_2029, lay_var_2030) = fracrefb_var_175(ig_var_2029)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2036 /= laytrop_min_var_2035) THEN
    DO lay_var_2030 = laytrop_min_var_2035 + 1, laytrop_max_var_2036
      ixc0_var_2042 = ixc_var_2037(lay_var_2030)
      DO ixp_var_2043 = 1, ixc0_var_2042
        jl_var_2045 = ixlow_var_2038(ixp_var_2043, lay_var_2030)
        ind0_var_2025 = ((jp_var_2015(jl_var_2045, lay_var_2030) - 1) * 5 + (jt_var_2016(jl_var_2045, lay_var_2030) - 1)) * nspa_var_237(2) + 1
        ind1_var_2026 = (jp_var_2015(jl_var_2045, lay_var_2030) * 5 + (jt1_var_2017(jl_var_2045, lay_var_2030) - 1)) * nspa_var_237(2) + 1
        inds_var_2027 = indself_var_2022(jl_var_2045, lay_var_2030)
        indf_var_2028 = indfor_var_2023(jl_var_2045, lay_var_2030)
        pp_var_2034 = pavel_var_2006(jl_var_2045, lay_var_2030)
        corradj_var_2033 = 1.0D0 - 0.05D0 * (pp_var_2034 - 100.0D0) / 900.0D0
        DO ig_var_2029 = 1, 12
          tauself_var_2032 = selffac_var_2020(jl_var_2045, lay_var_2030) * (selfref_var_178(inds_var_2027, ig_var_2029) + selffrac_var_2021(jl_var_2045, lay_var_2030) * (selfref_var_178(inds_var_2027 + 1, ig_var_2029) - selfref_var_178(inds_var_2027, ig_var_2029)))
          taufor_var_2031 = forfac_var_2014(jl_var_2045, lay_var_2030) * (forref_var_179(indf_var_2028, ig_var_2029) + forfrac_var_2013(jl_var_2045, lay_var_2030) * (forref_var_179(indf_var_2028 + 1, ig_var_2029) - forref_var_179(indf_var_2028, ig_var_2029)))
          taug_var_2005(jl_var_2045, 10 + ig_var_2029, lay_var_2030) = corradj_var_2033 * (colh2o_var_2018(jl_var_2045, lay_var_2030) * (fac00_var_2009(jl_var_2045, lay_var_2030) * absa_var_176(ind0_var_2025, ig_var_2029) + fac10_var_2011(jl_var_2045, lay_var_2030) * absa_var_176(ind0_var_2025 + 1, ig_var_2029) + fac01_var_2010(jl_var_2045, lay_var_2030) * absa_var_176(ind1_var_2026, ig_var_2029) + fac11_var_2012(jl_var_2045, lay_var_2030) * absa_var_176(ind1_var_2026 + 1, ig_var_2029)) + tauself_var_2032 + taufor_var_2031)
          fracs_var_2024(jl_var_2045, 10 + ig_var_2029, lay_var_2030) = fracrefa_var_174(ig_var_2029)
        END DO
      END DO
      ixc0_var_2042 = kfdia_var_2003 - kidia_var_2002 + 1 - ixc0_var_2042
      DO ixp_var_2043 = 1, ixc0_var_2042
        jl_var_2045 = ixhigh_var_2039(ixp_var_2043, lay_var_2030)
        ind0_var_2025 = ((jp_var_2015(jl_var_2045, lay_var_2030) - 13) * 5 + (jt_var_2016(jl_var_2045, lay_var_2030) - 1)) * nspb_var_238(2) + 1
        ind1_var_2026 = ((jp_var_2015(jl_var_2045, lay_var_2030) - 12) * 5 + (jt1_var_2017(jl_var_2045, lay_var_2030) - 1)) * nspb_var_238(2) + 1
        indf_var_2028 = indfor_var_2023(jl_var_2045, lay_var_2030)
        DO ig_var_2029 = 1, 12
          taufor_var_2031 = forfac_var_2014(jl_var_2045, lay_var_2030) * (forref_var_179(indf_var_2028, ig_var_2029) + forfrac_var_2013(jl_var_2045, lay_var_2030) * (forref_var_179(indf_var_2028 + 1, ig_var_2029) - forref_var_179(indf_var_2028, ig_var_2029)))
          taug_var_2005(jl_var_2045, 10 + ig_var_2029, lay_var_2030) = colh2o_var_2018(jl_var_2045, lay_var_2030) * (fac00_var_2009(jl_var_2045, lay_var_2030) * absb_var_177(ind0_var_2025, ig_var_2029) + fac10_var_2011(jl_var_2045, lay_var_2030) * absb_var_177(ind0_var_2025 + 1, ig_var_2029) + fac01_var_2010(jl_var_2045, lay_var_2030) * absb_var_177(ind1_var_2026, ig_var_2029) + fac11_var_2012(jl_var_2045, lay_var_2030) * absb_var_177(ind1_var_2026 + 1, ig_var_2029)) + taufor_var_2031
          fracs_var_2024(jl_var_2045, 10 + ig_var_2029, lay_var_2030) = fracrefb_var_175(ig_var_2029)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol2
SUBROUTINE rrtm_taumol1(kidia_var_2046, kfdia_var_2047, klev_var_2048, taug_var_2050, pavel_var_2049, p_tauaerl_var_2051, fac00_var_2052, fac01_var_2053, fac10_var_2054, fac11_var_2055, forfac_var_2056, forfrac_var_2057, indfor_var_2068, jp_var_2058, jt_var_2059, jt1_var_2060, colh2o_var_2061, laytrop_var_2062, selffac_var_2063, selffrac_var_2064, indself_var_2066, fracs_var_2067, minorfrac_var_2065, indminor_var_2069, scaleminorn2, colbrd_var_2070)
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrta1, ONLY: absa_var_130, absb_var_131, forref_var_134, fracrefa_var_128, fracrefb_var_129, ka_mn2_var_132, kb_mn2, selfref_var_133
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2046
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2047
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2048
  REAL(KIND = 8), INTENT(IN) :: pavel_var_2049(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2050(kidia_var_2046 : kfdia_var_2047, 140, klev_var_2048)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2051(kidia_var_2046 : kfdia_var_2047, klev_var_2048, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2052(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2053(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2054(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2055(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2056(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2057(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2058(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2059(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2060(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2061(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2062(kidia_var_2046 : kfdia_var_2047)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2063(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2064(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2065(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2066(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2067(kidia_var_2046 : kfdia_var_2047, 140, klev_var_2048)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2068(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2069(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(IN) :: scaleminorn2(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  REAL(KIND = 8), INTENT(IN) :: colbrd_var_2070(kidia_var_2046 : kfdia_var_2047, klev_var_2048)
  INTEGER(KIND = 4) :: ind0_var_2071, ind1_var_2072, inds_var_2073
  INTEGER(KIND = 4) :: indf_var_2074, indm_var_2075
  INTEGER(KIND = 4) :: ig_var_2076, lay_var_2077
  REAL(KIND = 8) :: taufor_var_2078, tauself_var_2079, corradj_var_2080, pp_var_2081, scalen2_var_2082, taun2_var_2083
  INTEGER(KIND = 4) :: laytrop_min_var_2084, laytrop_max_var_2085
  INTEGER(KIND = 4) :: ixc_var_2086(klev_var_2048), ixlow_var_2087(kfdia_var_2047, klev_var_2048), ixhigh_var_2088(kfdia_var_2047, klev_var_2048)
  INTEGER(KIND = 4) :: ich_var_2089, icl_var_2090, ixc0_var_2091, ixp_var_2092, jc_var_2093, jl_var_2094
  laytrop_min_var_2084 = MINVAL(laytrop_var_2062)
  laytrop_max_var_2085 = MAXVAL(laytrop_var_2062)
  ixlow_var_2087 = 0
  ixhigh_var_2088 = 0
  ixc_var_2086 = 0
  DO lay_var_2077 = laytrop_min_var_2084 + 1, laytrop_max_var_2085
    icl_var_2090 = 0
    ich_var_2089 = 0
    DO jc_var_2093 = kidia_var_2046, kfdia_var_2047
      IF (lay_var_2077 <= laytrop_var_2062(jc_var_2093)) THEN
        icl_var_2090 = icl_var_2090 + 1
        ixlow_var_2087(icl_var_2090, lay_var_2077) = jc_var_2093
      ELSE
        ich_var_2089 = ich_var_2089 + 1
        ixhigh_var_2088(ich_var_2089, lay_var_2077) = jc_var_2093
      END IF
    END DO
    ixc_var_2086(lay_var_2077) = icl_var_2090
  END DO
  DO lay_var_2077 = 1, laytrop_min_var_2084
    DO jl_var_2094 = kidia_var_2046, kfdia_var_2047
      ind0_var_2071 = ((jp_var_2058(jl_var_2094, lay_var_2077) - 1) * 5 + (jt_var_2059(jl_var_2094, lay_var_2077) - 1)) * nspa_var_237(1) + 1
      ind1_var_2072 = (jp_var_2058(jl_var_2094, lay_var_2077) * 5 + (jt1_var_2060(jl_var_2094, lay_var_2077) - 1)) * nspa_var_237(1) + 1
      inds_var_2073 = indself_var_2066(jl_var_2094, lay_var_2077)
      indf_var_2074 = indfor_var_2068(jl_var_2094, lay_var_2077)
      indm_var_2075 = indminor_var_2069(jl_var_2094, lay_var_2077)
      pp_var_2081 = pavel_var_2049(jl_var_2094, lay_var_2077)
      corradj_var_2080 = 1.0D0
      IF (pp_var_2081 .LT. 250.0D0) THEN
        corradj_var_2080 = 1.0D0 - 0.15D0 * (250.0D0 - pp_var_2081) / 154.4D0
      END IF
      scalen2_var_2082 = colbrd_var_2070(jl_var_2094, lay_var_2077) * scaleminorn2(jl_var_2094, lay_var_2077)
      DO ig_var_2076 = 1, 10
        tauself_var_2079 = selffac_var_2063(jl_var_2094, lay_var_2077) * (selfref_var_133(inds_var_2073, ig_var_2076) + selffrac_var_2064(jl_var_2094, lay_var_2077) * (selfref_var_133(inds_var_2073 + 1, ig_var_2076) - selfref_var_133(inds_var_2073, ig_var_2076)))
        taufor_var_2078 = forfac_var_2056(jl_var_2094, lay_var_2077) * (forref_var_134(indf_var_2074, ig_var_2076) + forfrac_var_2057(jl_var_2094, lay_var_2077) * (forref_var_134(indf_var_2074 + 1, ig_var_2076) - forref_var_134(indf_var_2074, ig_var_2076)))
        taun2_var_2083 = scalen2_var_2082 * (ka_mn2_var_132(indm_var_2075, ig_var_2076) + minorfrac_var_2065(jl_var_2094, lay_var_2077) * (ka_mn2_var_132(indm_var_2075 + 1, ig_var_2076) - ka_mn2_var_132(indm_var_2075, ig_var_2076)))
        taug_var_2050(jl_var_2094, ig_var_2076, lay_var_2077) = corradj_var_2080 * (colh2o_var_2061(jl_var_2094, lay_var_2077) * (fac00_var_2052(jl_var_2094, lay_var_2077) * absa_var_130(ind0_var_2071, ig_var_2076) + fac10_var_2054(jl_var_2094, lay_var_2077) * absa_var_130(ind0_var_2071 + 1, ig_var_2076) + fac01_var_2053(jl_var_2094, lay_var_2077) * absa_var_130(ind1_var_2072, ig_var_2076) + fac11_var_2055(jl_var_2094, lay_var_2077) * absa_var_130(ind1_var_2072 + 1, ig_var_2076)) + tauself_var_2079 + taufor_var_2078 + taun2_var_2083)
        fracs_var_2067(jl_var_2094, ig_var_2076, lay_var_2077) = fracrefa_var_128(ig_var_2076)
      END DO
    END DO
  END DO
  DO lay_var_2077 = laytrop_max_var_2085 + 1, klev_var_2048
    DO jl_var_2094 = kidia_var_2046, kfdia_var_2047
      ind0_var_2071 = ((jp_var_2058(jl_var_2094, lay_var_2077) - 13) * 5 + (jt_var_2059(jl_var_2094, lay_var_2077) - 1)) * nspb_var_238(1) + 1
      ind1_var_2072 = ((jp_var_2058(jl_var_2094, lay_var_2077) - 12) * 5 + (jt1_var_2060(jl_var_2094, lay_var_2077) - 1)) * nspb_var_238(1) + 1
      indf_var_2074 = indfor_var_2068(jl_var_2094, lay_var_2077)
      indm_var_2075 = indminor_var_2069(jl_var_2094, lay_var_2077)
      pp_var_2081 = pavel_var_2049(jl_var_2094, lay_var_2077)
      corradj_var_2080 = 1.0D0 - 0.15D0 * (pp_var_2081 / 95.6D0)
      scalen2_var_2082 = colbrd_var_2070(jl_var_2094, lay_var_2077) * scaleminorn2(jl_var_2094, lay_var_2077)
      DO ig_var_2076 = 1, 10
        taufor_var_2078 = forfac_var_2056(jl_var_2094, lay_var_2077) * (forref_var_134(indf_var_2074, ig_var_2076) + forfrac_var_2057(jl_var_2094, lay_var_2077) * (forref_var_134(indf_var_2074 + 1, ig_var_2076) - forref_var_134(indf_var_2074, ig_var_2076)))
        taun2_var_2083 = scalen2_var_2082 * (kb_mn2(indm_var_2075, ig_var_2076) + minorfrac_var_2065(jl_var_2094, lay_var_2077) * (kb_mn2(indm_var_2075 + 1, ig_var_2076) - kb_mn2(indm_var_2075, ig_var_2076)))
        taug_var_2050(jl_var_2094, ig_var_2076, lay_var_2077) = corradj_var_2080 * (colh2o_var_2061(jl_var_2094, lay_var_2077) * (fac00_var_2052(jl_var_2094, lay_var_2077) * absb_var_131(ind0_var_2071, ig_var_2076) + fac10_var_2054(jl_var_2094, lay_var_2077) * absb_var_131(ind0_var_2071 + 1, ig_var_2076) + fac01_var_2053(jl_var_2094, lay_var_2077) * absb_var_131(ind1_var_2072, ig_var_2076) + fac11_var_2055(jl_var_2094, lay_var_2077) * absb_var_131(ind1_var_2072 + 1, ig_var_2076)) + taufor_var_2078 + taun2_var_2083)
        fracs_var_2067(jl_var_2094, ig_var_2076, lay_var_2077) = fracrefb_var_129(ig_var_2076)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2085 /= laytrop_min_var_2084) THEN
    DO lay_var_2077 = laytrop_min_var_2084 + 1, laytrop_max_var_2085
      ixc0_var_2091 = ixc_var_2086(lay_var_2077)
      DO ixp_var_2092 = 1, ixc0_var_2091
        jl_var_2094 = ixlow_var_2087(ixp_var_2092, lay_var_2077)
        ind0_var_2071 = ((jp_var_2058(jl_var_2094, lay_var_2077) - 1) * 5 + (jt_var_2059(jl_var_2094, lay_var_2077) - 1)) * nspa_var_237(1) + 1
        ind1_var_2072 = (jp_var_2058(jl_var_2094, lay_var_2077) * 5 + (jt1_var_2060(jl_var_2094, lay_var_2077) - 1)) * nspa_var_237(1) + 1
        inds_var_2073 = indself_var_2066(jl_var_2094, lay_var_2077)
        indf_var_2074 = indfor_var_2068(jl_var_2094, lay_var_2077)
        indm_var_2075 = indminor_var_2069(jl_var_2094, lay_var_2077)
        pp_var_2081 = pavel_var_2049(jl_var_2094, lay_var_2077)
        corradj_var_2080 = 1.0D0
        IF (pp_var_2081 .LT. 250.0D0) THEN
          corradj_var_2080 = 1.0D0 - 0.15D0 * (250.0D0 - pp_var_2081) / 154.4D0
        END IF
        scalen2_var_2082 = colbrd_var_2070(jl_var_2094, lay_var_2077) * scaleminorn2(jl_var_2094, lay_var_2077)
        DO ig_var_2076 = 1, 10
          tauself_var_2079 = selffac_var_2063(jl_var_2094, lay_var_2077) * (selfref_var_133(inds_var_2073, ig_var_2076) + selffrac_var_2064(jl_var_2094, lay_var_2077) * (selfref_var_133(inds_var_2073 + 1, ig_var_2076) - selfref_var_133(inds_var_2073, ig_var_2076)))
          taufor_var_2078 = forfac_var_2056(jl_var_2094, lay_var_2077) * (forref_var_134(indf_var_2074, ig_var_2076) + forfrac_var_2057(jl_var_2094, lay_var_2077) * (forref_var_134(indf_var_2074 + 1, ig_var_2076) - forref_var_134(indf_var_2074, ig_var_2076)))
          taun2_var_2083 = scalen2_var_2082 * (ka_mn2_var_132(indm_var_2075, ig_var_2076) + minorfrac_var_2065(jl_var_2094, lay_var_2077) * (ka_mn2_var_132(indm_var_2075 + 1, ig_var_2076) - ka_mn2_var_132(indm_var_2075, ig_var_2076)))
          taug_var_2050(jl_var_2094, ig_var_2076, lay_var_2077) = corradj_var_2080 * (colh2o_var_2061(jl_var_2094, lay_var_2077) * (fac00_var_2052(jl_var_2094, lay_var_2077) * absa_var_130(ind0_var_2071, ig_var_2076) + fac10_var_2054(jl_var_2094, lay_var_2077) * absa_var_130(ind0_var_2071 + 1, ig_var_2076) + fac01_var_2053(jl_var_2094, lay_var_2077) * absa_var_130(ind1_var_2072, ig_var_2076) + fac11_var_2055(jl_var_2094, lay_var_2077) * absa_var_130(ind1_var_2072 + 1, ig_var_2076)) + tauself_var_2079 + taufor_var_2078 + taun2_var_2083)
          fracs_var_2067(jl_var_2094, ig_var_2076, lay_var_2077) = fracrefa_var_128(ig_var_2076)
        END DO
      END DO
      ixc0_var_2091 = kfdia_var_2047 - kidia_var_2046 + 1 - ixc0_var_2091
      DO ixp_var_2092 = 1, ixc0_var_2091
        jl_var_2094 = ixhigh_var_2088(ixp_var_2092, lay_var_2077)
        ind0_var_2071 = ((jp_var_2058(jl_var_2094, lay_var_2077) - 13) * 5 + (jt_var_2059(jl_var_2094, lay_var_2077) - 1)) * nspb_var_238(1) + 1
        ind1_var_2072 = ((jp_var_2058(jl_var_2094, lay_var_2077) - 12) * 5 + (jt1_var_2060(jl_var_2094, lay_var_2077) - 1)) * nspb_var_238(1) + 1
        indf_var_2074 = indfor_var_2068(jl_var_2094, lay_var_2077)
        indm_var_2075 = indminor_var_2069(jl_var_2094, lay_var_2077)
        pp_var_2081 = pavel_var_2049(jl_var_2094, lay_var_2077)
        corradj_var_2080 = 1.0D0 - 0.15D0 * (pp_var_2081 / 95.6D0)
        scalen2_var_2082 = colbrd_var_2070(jl_var_2094, lay_var_2077) * scaleminorn2(jl_var_2094, lay_var_2077)
        DO ig_var_2076 = 1, 10
          taufor_var_2078 = forfac_var_2056(jl_var_2094, lay_var_2077) * (forref_var_134(indf_var_2074, ig_var_2076) + forfrac_var_2057(jl_var_2094, lay_var_2077) * (forref_var_134(indf_var_2074 + 1, ig_var_2076) - forref_var_134(indf_var_2074, ig_var_2076)))
          taun2_var_2083 = scalen2_var_2082 * (kb_mn2(indm_var_2075, ig_var_2076) + minorfrac_var_2065(jl_var_2094, lay_var_2077) * (kb_mn2(indm_var_2075 + 1, ig_var_2076) - kb_mn2(indm_var_2075, ig_var_2076)))
          taug_var_2050(jl_var_2094, ig_var_2076, lay_var_2077) = corradj_var_2080 * (colh2o_var_2061(jl_var_2094, lay_var_2077) * (fac00_var_2052(jl_var_2094, lay_var_2077) * absb_var_131(ind0_var_2071, ig_var_2076) + fac10_var_2054(jl_var_2094, lay_var_2077) * absb_var_131(ind0_var_2071 + 1, ig_var_2076) + fac01_var_2053(jl_var_2094, lay_var_2077) * absb_var_131(ind1_var_2072, ig_var_2076) + fac11_var_2055(jl_var_2094, lay_var_2077) * absb_var_131(ind1_var_2072 + 1, ig_var_2076)) + taufor_var_2078 + taun2_var_2083)
          fracs_var_2067(jl_var_2094, ig_var_2076, lay_var_2077) = fracrefb_var_129(ig_var_2076)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol1
SUBROUTINE rrtm_taumol14(kidia_var_2095, kfdia_var_2096, klev_var_2097, taug_var_2098, p_tauaerl_var_2099, fac00_var_2100, fac01_var_2101, fac10_var_2102, fac11_var_2103, forfac_var_2114, forfrac_var_2115, indfor_var_2113, jp_var_2104, jt_var_2105, jt1_var_2106, colco2_var_2107, laytrop_var_2108, selffac_var_2109, selffrac_var_2110, indself_var_2111, fracs_var_2112)
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrta14, ONLY: absa_var_159, absb_var_160, forref_var_162, fracrefa_var_157, fracrefb_var_158, selfref_var_161
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2095
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2096
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2097
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2098(kidia_var_2095 : kfdia_var_2096, 140, klev_var_2097)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2099(kidia_var_2095 : kfdia_var_2096, klev_var_2097, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2100(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2101(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2102(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2103(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2104(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2105(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2106(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2107(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2108(kidia_var_2095 : kfdia_var_2096)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2109(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2110(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2111(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2112(kidia_var_2095 : kfdia_var_2096, 140, klev_var_2097)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2113(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2114(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2115(kidia_var_2095 : kfdia_var_2096, klev_var_2097)
  INTEGER(KIND = 4) :: ig_var_2116, ind0_var_2117, ind1_var_2118, inds_var_2119, indf_var_2120, lay_var_2121
  REAL(KIND = 8) :: taufor_var_2122, tauself_var_2123
  INTEGER(KIND = 4) :: laytrop_min_var_2124, laytrop_max_var_2125
  INTEGER(KIND = 4) :: ixc_var_2126(klev_var_2097), ixlow_var_2127(kfdia_var_2096, klev_var_2097), ixhigh_var_2128(kfdia_var_2096, klev_var_2097)
  INTEGER(KIND = 4) :: ich_var_2129, icl_var_2130, ixc0_var_2131, ixp_var_2132, jc_var_2133, jl_var_2134
  laytrop_min_var_2124 = MINVAL(laytrop_var_2108)
  laytrop_max_var_2125 = MAXVAL(laytrop_var_2108)
  ixlow_var_2127 = 0
  ixhigh_var_2128 = 0
  ixc_var_2126 = 0
  DO lay_var_2121 = laytrop_min_var_2124 + 1, laytrop_max_var_2125
    icl_var_2130 = 0
    ich_var_2129 = 0
    DO jc_var_2133 = kidia_var_2095, kfdia_var_2096
      IF (lay_var_2121 <= laytrop_var_2108(jc_var_2133)) THEN
        icl_var_2130 = icl_var_2130 + 1
        ixlow_var_2127(icl_var_2130, lay_var_2121) = jc_var_2133
      ELSE
        ich_var_2129 = ich_var_2129 + 1
        ixhigh_var_2128(ich_var_2129, lay_var_2121) = jc_var_2133
      END IF
    END DO
    ixc_var_2126(lay_var_2121) = icl_var_2130
  END DO
  DO lay_var_2121 = 1, laytrop_min_var_2124
    DO jl_var_2134 = kidia_var_2095, kfdia_var_2096
      ind0_var_2117 = ((jp_var_2104(jl_var_2134, lay_var_2121) - 1) * 5 + (jt_var_2105(jl_var_2134, lay_var_2121) - 1)) * nspa_var_237(14) + 1
      ind1_var_2118 = (jp_var_2104(jl_var_2134, lay_var_2121) * 5 + (jt1_var_2106(jl_var_2134, lay_var_2121) - 1)) * nspa_var_237(14) + 1
      inds_var_2119 = indself_var_2111(jl_var_2134, lay_var_2121)
      indf_var_2120 = indfor_var_2113(jl_var_2134, lay_var_2121)
      DO ig_var_2116 = 1, 2
        tauself_var_2123 = selffac_var_2109(jl_var_2134, lay_var_2121) * (selfref_var_161(inds_var_2119, ig_var_2116) + selffrac_var_2110(jl_var_2134, lay_var_2121) * (selfref_var_161(inds_var_2119 + 1, ig_var_2116) - selfref_var_161(inds_var_2119, ig_var_2116)))
        taufor_var_2122 = forfac_var_2114(jl_var_2134, lay_var_2121) * (forref_var_162(indf_var_2120, ig_var_2116) + forfrac_var_2115(jl_var_2134, lay_var_2121) * (forref_var_162(indf_var_2120 + 1, ig_var_2116) - forref_var_162(indf_var_2120, ig_var_2116)))
        taug_var_2098(jl_var_2134, 134 + ig_var_2116, lay_var_2121) = colco2_var_2107(jl_var_2134, lay_var_2121) * (fac00_var_2100(jl_var_2134, lay_var_2121) * absa_var_159(ind0_var_2117, ig_var_2116) + fac10_var_2102(jl_var_2134, lay_var_2121) * absa_var_159(ind0_var_2117 + 1, ig_var_2116) + fac01_var_2101(jl_var_2134, lay_var_2121) * absa_var_159(ind1_var_2118, ig_var_2116) + fac11_var_2103(jl_var_2134, lay_var_2121) * absa_var_159(ind1_var_2118 + 1, ig_var_2116)) + tauself_var_2123 + taufor_var_2122
        fracs_var_2112(jl_var_2134, 134 + ig_var_2116, lay_var_2121) = fracrefa_var_157(ig_var_2116)
      END DO
    END DO
  END DO
  DO lay_var_2121 = laytrop_max_var_2125 + 1, klev_var_2097
    DO jl_var_2134 = kidia_var_2095, kfdia_var_2096
      ind0_var_2117 = ((jp_var_2104(jl_var_2134, lay_var_2121) - 13) * 5 + (jt_var_2105(jl_var_2134, lay_var_2121) - 1)) * nspb_var_238(14) + 1
      ind1_var_2118 = ((jp_var_2104(jl_var_2134, lay_var_2121) - 12) * 5 + (jt1_var_2106(jl_var_2134, lay_var_2121) - 1)) * nspb_var_238(14) + 1
      DO ig_var_2116 = 1, 2
        taug_var_2098(jl_var_2134, 134 + ig_var_2116, lay_var_2121) = colco2_var_2107(jl_var_2134, lay_var_2121) * (fac00_var_2100(jl_var_2134, lay_var_2121) * absb_var_160(ind0_var_2117, ig_var_2116) + fac10_var_2102(jl_var_2134, lay_var_2121) * absb_var_160(ind0_var_2117 + 1, ig_var_2116) + fac01_var_2101(jl_var_2134, lay_var_2121) * absb_var_160(ind1_var_2118, ig_var_2116) + fac11_var_2103(jl_var_2134, lay_var_2121) * absb_var_160(ind1_var_2118 + 1, ig_var_2116))
        fracs_var_2112(jl_var_2134, 134 + ig_var_2116, lay_var_2121) = fracrefb_var_158(ig_var_2116)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2125 /= laytrop_min_var_2124) THEN
    DO lay_var_2121 = laytrop_min_var_2124 + 1, laytrop_max_var_2125
      ixc0_var_2131 = ixc_var_2126(lay_var_2121)
      DO ixp_var_2132 = 1, ixc0_var_2131
        jl_var_2134 = ixlow_var_2127(ixp_var_2132, lay_var_2121)
        ind0_var_2117 = ((jp_var_2104(jl_var_2134, lay_var_2121) - 1) * 5 + (jt_var_2105(jl_var_2134, lay_var_2121) - 1)) * nspa_var_237(14) + 1
        ind1_var_2118 = (jp_var_2104(jl_var_2134, lay_var_2121) * 5 + (jt1_var_2106(jl_var_2134, lay_var_2121) - 1)) * nspa_var_237(14) + 1
        inds_var_2119 = indself_var_2111(jl_var_2134, lay_var_2121)
        indf_var_2120 = indfor_var_2113(jl_var_2134, lay_var_2121)
        DO ig_var_2116 = 1, 2
          tauself_var_2123 = selffac_var_2109(jl_var_2134, lay_var_2121) * (selfref_var_161(inds_var_2119, ig_var_2116) + selffrac_var_2110(jl_var_2134, lay_var_2121) * (selfref_var_161(inds_var_2119 + 1, ig_var_2116) - selfref_var_161(inds_var_2119, ig_var_2116)))
          taufor_var_2122 = forfac_var_2114(jl_var_2134, lay_var_2121) * (forref_var_162(indf_var_2120, ig_var_2116) + forfrac_var_2115(jl_var_2134, lay_var_2121) * (forref_var_162(indf_var_2120 + 1, ig_var_2116) - forref_var_162(indf_var_2120, ig_var_2116)))
          taug_var_2098(jl_var_2134, 134 + ig_var_2116, lay_var_2121) = colco2_var_2107(jl_var_2134, lay_var_2121) * (fac00_var_2100(jl_var_2134, lay_var_2121) * absa_var_159(ind0_var_2117, ig_var_2116) + fac10_var_2102(jl_var_2134, lay_var_2121) * absa_var_159(ind0_var_2117 + 1, ig_var_2116) + fac01_var_2101(jl_var_2134, lay_var_2121) * absa_var_159(ind1_var_2118, ig_var_2116) + fac11_var_2103(jl_var_2134, lay_var_2121) * absa_var_159(ind1_var_2118 + 1, ig_var_2116)) + tauself_var_2123 + taufor_var_2122
          fracs_var_2112(jl_var_2134, 134 + ig_var_2116, lay_var_2121) = fracrefa_var_157(ig_var_2116)
        END DO
      END DO
      ixc0_var_2131 = kfdia_var_2096 - kidia_var_2095 + 1 - ixc0_var_2131
      DO ixp_var_2132 = 1, ixc0_var_2131
        jl_var_2134 = ixhigh_var_2128(ixp_var_2132, lay_var_2121)
        ind0_var_2117 = ((jp_var_2104(jl_var_2134, lay_var_2121) - 13) * 5 + (jt_var_2105(jl_var_2134, lay_var_2121) - 1)) * nspb_var_238(14) + 1
        ind1_var_2118 = ((jp_var_2104(jl_var_2134, lay_var_2121) - 12) * 5 + (jt1_var_2106(jl_var_2134, lay_var_2121) - 1)) * nspb_var_238(14) + 1
        DO ig_var_2116 = 1, 2
          taug_var_2098(jl_var_2134, 134 + ig_var_2116, lay_var_2121) = colco2_var_2107(jl_var_2134, lay_var_2121) * (fac00_var_2100(jl_var_2134, lay_var_2121) * absb_var_160(ind0_var_2117, ig_var_2116) + fac10_var_2102(jl_var_2134, lay_var_2121) * absb_var_160(ind0_var_2117 + 1, ig_var_2116) + fac01_var_2101(jl_var_2134, lay_var_2121) * absb_var_160(ind1_var_2118, ig_var_2116) + fac11_var_2103(jl_var_2134, lay_var_2121) * absb_var_160(ind1_var_2118 + 1, ig_var_2116))
          fracs_var_2112(jl_var_2134, 134 + ig_var_2116, lay_var_2121) = fracrefb_var_158(ig_var_2116)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol14
SUBROUTINE rrtm_taumol15(kidia_var_2135, kfdia_var_2136, klev_var_2137, taug_var_2138, p_tauaerl_var_2139, fac00_var_2140, fac01_var_2141, fac10_var_2142, fac11_var_2143, forfac_var_2157, forfrac_var_2158, indfor_var_2156, jp_var_2144, jt_var_2145, jt1_var_2146, oneminus_var_2147, colh2o_var_2148, colco2_var_2149, coln2o_var_2150, laytrop_var_2151, selffac_var_2152, selffrac_var_2153, indself_var_2154, fracs_var_2155, rat_n2oco2, rat_n2oco2_1, minorfrac_var_2159, indminor_var_2160, scaleminor_var_2161, colbrd_var_2162)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_237
  USE yoerrtm, ONLY: ng15
  USE yoerrta15, ONLY: absa_var_164, forref_var_167, fracrefa_var_163, ka_mn2_var_165, selfref_var_166
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2135
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2136
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2137
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2138(kidia_var_2135 : kfdia_var_2136, 140, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2139(kidia_var_2135 : kfdia_var_2136, klev_var_2137, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2140(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2141(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2142(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2143(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2144(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2145(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2146(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2147
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2148(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2149(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_2150(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2151(kidia_var_2135 : kfdia_var_2136)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2152(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2153(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2154(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2155(kidia_var_2135 : kfdia_var_2136, 140, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: rat_n2oco2(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: rat_n2oco2_1(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2156(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2157(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2158(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2159(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2160(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: scaleminor_var_2161(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  REAL(KIND = 8), INTENT(IN) :: colbrd_var_2162(kidia_var_2135 : kfdia_var_2136, klev_var_2137)
  INTEGER(KIND = 4) :: ig_var_2163, ind0_var_2164, ind1_var_2165, inds_var_2166, indf_var_2167, indm_var_2168, js_var_2169, js1_var_2170, jpl_var_2171, jmn2, lay_var_2172
  REAL(KIND = 8) :: refrat_planck_a_var_2173, refrat_m_a_var_2174
  REAL(KIND = 8) :: taufor_var_2175, tauself_var_2176, tau_major_var_2177(2), tau_major1_var_2178(2), n2m1, n2m2, taun2_var_2179, scalen2_var_2180
  REAL(KIND = 8) :: fac000_var_2181, fac100_var_2182, fac200_var_2183, fac010_var_2184, fac110_var_2185, fac210_var_2186, fac001_var_2187, fac101_var_2188, fac201_var_2189, fac011_var_2190, fac111_var_2191, fac211_var_2192
  REAL(KIND = 8) :: p_var_2193, p4_var_2194, fk0_var_2195, fk1_var_2196, fk2_var_2197
  REAL(KIND = 8) :: fs_var_2198, specmult_var_2199, specparm_var_2200, speccomb_var_2201, fs1_var_2202, specmult1_var_2203, specparm1_var_2204, speccomb1_var_2205, fmn2, specmult_mn2, specparm_mn2, speccomb_mn2, fpl_var_2206, specmult_planck_var_2207, specparm_planck_var_2208, speccomb_planck_var_2209
  INTEGER(KIND = 4) :: laytrop_min_var_2210, laytrop_max_var_2211
  INTEGER(KIND = 4) :: ixc_var_2212(klev_var_2137), ixlow_var_2213(kfdia_var_2136, klev_var_2137), ixhigh_var_2214(kfdia_var_2136, klev_var_2137)
  INTEGER(KIND = 4) :: ich_var_2215, icl_var_2216, ixc0_var_2217, ixp_var_2218, jc_var_2219, jl_var_2220
  laytrop_min_var_2210 = MINVAL(laytrop_var_2151)
  laytrop_max_var_2211 = MAXVAL(laytrop_var_2151)
  ixlow_var_2213 = 0
  ixhigh_var_2214 = 0
  ixc_var_2212 = 0
  DO lay_var_2172 = laytrop_min_var_2210 + 1, laytrop_max_var_2211
    icl_var_2216 = 0
    ich_var_2215 = 0
    DO jc_var_2219 = kidia_var_2135, kfdia_var_2136
      IF (lay_var_2172 <= laytrop_var_2151(jc_var_2219)) THEN
        icl_var_2216 = icl_var_2216 + 1
        ixlow_var_2213(icl_var_2216, lay_var_2172) = jc_var_2219
      ELSE
        ich_var_2215 = ich_var_2215 + 1
        ixhigh_var_2214(ich_var_2215, lay_var_2172) = jc_var_2219
      END IF
    END DO
    ixc_var_2212(lay_var_2172) = icl_var_2216
  END DO
  refrat_planck_a_var_2173 = chi_mls(4, 1) / chi_mls(2, 1)
  refrat_m_a_var_2174 = chi_mls(4, 1) / chi_mls(2, 1)
  DO lay_var_2172 = 1, laytrop_min_var_2210
    DO jl_var_2220 = kidia_var_2135, kfdia_var_2136
      speccomb_var_2201 = coln2o_var_2150(jl_var_2220, lay_var_2172) + rat_n2oco2(jl_var_2220, lay_var_2172) * colco2_var_2149(jl_var_2220, lay_var_2172)
      specparm_var_2200 = MIN(coln2o_var_2150(jl_var_2220, lay_var_2172) / speccomb_var_2201, oneminus_var_2147)
      specmult_var_2199 = 8.0D0 * (specparm_var_2200)
      js_var_2169 = 1 + INT(specmult_var_2199)
      fs_var_2198 = ((specmult_var_2199) - AINT((specmult_var_2199)))
      speccomb1_var_2205 = coln2o_var_2150(jl_var_2220, lay_var_2172) + rat_n2oco2_1(jl_var_2220, lay_var_2172) * colco2_var_2149(jl_var_2220, lay_var_2172)
      specparm1_var_2204 = MIN(coln2o_var_2150(jl_var_2220, lay_var_2172) / speccomb1_var_2205, oneminus_var_2147)
      specmult1_var_2203 = 8.0D0 * (specparm1_var_2204)
      js1_var_2170 = 1 + INT(specmult1_var_2203)
      fs1_var_2202 = ((specmult1_var_2203) - AINT((specmult1_var_2203)))
      speccomb_mn2 = coln2o_var_2150(jl_var_2220, lay_var_2172) + refrat_m_a_var_2174 * colco2_var_2149(jl_var_2220, lay_var_2172)
      specparm_mn2 = MIN(coln2o_var_2150(jl_var_2220, lay_var_2172) / speccomb_mn2, oneminus_var_2147)
      specmult_mn2 = 8.0D0 * specparm_mn2
      jmn2 = 1 + INT(specmult_mn2)
      fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))
      speccomb_planck_var_2209 = coln2o_var_2150(jl_var_2220, lay_var_2172) + refrat_planck_a_var_2173 * colco2_var_2149(jl_var_2220, lay_var_2172)
      specparm_planck_var_2208 = MIN(coln2o_var_2150(jl_var_2220, lay_var_2172) / speccomb_planck_var_2209, oneminus_var_2147)
      specmult_planck_var_2207 = 8.0D0 * specparm_planck_var_2208
      jpl_var_2171 = 1 + INT(specmult_planck_var_2207)
      fpl_var_2206 = ((specmult_planck_var_2207) - AINT((specmult_planck_var_2207)))
      ind0_var_2164 = ((jp_var_2144(jl_var_2220, lay_var_2172) - 1) * 5 + (jt_var_2145(jl_var_2220, lay_var_2172) - 1)) * nspa_var_237(15) + js_var_2169
      ind1_var_2165 = (jp_var_2144(jl_var_2220, lay_var_2172) * 5 + (jt1_var_2146(jl_var_2220, lay_var_2172) - 1)) * nspa_var_237(15) + js1_var_2170
      inds_var_2166 = indself_var_2154(jl_var_2220, lay_var_2172)
      indf_var_2167 = indfor_var_2156(jl_var_2220, lay_var_2172)
      indm_var_2168 = indminor_var_2160(jl_var_2220, lay_var_2172)
      scalen2_var_2180 = colbrd_var_2162(jl_var_2220, lay_var_2172) * scaleminor_var_2161(jl_var_2220, lay_var_2172)
      IF (specparm_var_2200 .LT. 0.125D0) THEN
        p_var_2193 = fs_var_2198 - 1.0D0
        p4_var_2194 = p_var_2193 ** 4
        fk0_var_2195 = p4_var_2194
        fk1_var_2196 = 1.0D0 - p_var_2193 - 2.0D0 * p4_var_2194
        fk2_var_2197 = p_var_2193 + p4_var_2194
        fac000_var_2181 = fk0_var_2195 * fac00_var_2140(jl_var_2220, lay_var_2172)
        fac100_var_2182 = fk1_var_2196 * fac00_var_2140(jl_var_2220, lay_var_2172)
        fac200_var_2183 = fk2_var_2197 * fac00_var_2140(jl_var_2220, lay_var_2172)
        fac010_var_2184 = fk0_var_2195 * fac10_var_2142(jl_var_2220, lay_var_2172)
        fac110_var_2185 = fk1_var_2196 * fac10_var_2142(jl_var_2220, lay_var_2172)
        fac210_var_2186 = fk2_var_2197 * fac10_var_2142(jl_var_2220, lay_var_2172)
      ELSE IF (specparm_var_2200 .GT. 0.875D0) THEN
        p_var_2193 = - fs_var_2198
        p4_var_2194 = p_var_2193 ** 4
        fk0_var_2195 = p4_var_2194
        fk1_var_2196 = 1.0D0 - p_var_2193 - 2.0D0 * p4_var_2194
        fk2_var_2197 = p_var_2193 + p4_var_2194
        fac000_var_2181 = fk0_var_2195 * fac00_var_2140(jl_var_2220, lay_var_2172)
        fac100_var_2182 = fk1_var_2196 * fac00_var_2140(jl_var_2220, lay_var_2172)
        fac200_var_2183 = fk2_var_2197 * fac00_var_2140(jl_var_2220, lay_var_2172)
        fac010_var_2184 = fk0_var_2195 * fac10_var_2142(jl_var_2220, lay_var_2172)
        fac110_var_2185 = fk1_var_2196 * fac10_var_2142(jl_var_2220, lay_var_2172)
        fac210_var_2186 = fk2_var_2197 * fac10_var_2142(jl_var_2220, lay_var_2172)
      ELSE
        fac000_var_2181 = (1.0D0 - fs_var_2198) * fac00_var_2140(jl_var_2220, lay_var_2172)
        fac010_var_2184 = (1.0D0 - fs_var_2198) * fac10_var_2142(jl_var_2220, lay_var_2172)
        fac100_var_2182 = fs_var_2198 * fac00_var_2140(jl_var_2220, lay_var_2172)
        fac110_var_2185 = fs_var_2198 * fac10_var_2142(jl_var_2220, lay_var_2172)
        fac200_var_2183 = 0.0D0
        fac210_var_2186 = 0.0D0
      END IF
      IF (specparm1_var_2204 .LT. 0.125D0) THEN
        p_var_2193 = fs1_var_2202 - 1.0D0
        p4_var_2194 = p_var_2193 ** 4
        fk0_var_2195 = p4_var_2194
        fk1_var_2196 = 1.0D0 - p_var_2193 - 2.0D0 * p4_var_2194
        fk2_var_2197 = p_var_2193 + p4_var_2194
        fac001_var_2187 = fk0_var_2195 * fac01_var_2141(jl_var_2220, lay_var_2172)
        fac101_var_2188 = fk1_var_2196 * fac01_var_2141(jl_var_2220, lay_var_2172)
        fac201_var_2189 = fk2_var_2197 * fac01_var_2141(jl_var_2220, lay_var_2172)
        fac011_var_2190 = fk0_var_2195 * fac11_var_2143(jl_var_2220, lay_var_2172)
        fac111_var_2191 = fk1_var_2196 * fac11_var_2143(jl_var_2220, lay_var_2172)
        fac211_var_2192 = fk2_var_2197 * fac11_var_2143(jl_var_2220, lay_var_2172)
      ELSE IF (specparm1_var_2204 .GT. 0.875D0) THEN
        p_var_2193 = - fs1_var_2202
        p4_var_2194 = p_var_2193 ** 4
        fk0_var_2195 = p4_var_2194
        fk1_var_2196 = 1.0D0 - p_var_2193 - 2.0D0 * p4_var_2194
        fk2_var_2197 = p_var_2193 + p4_var_2194
        fac001_var_2187 = fk0_var_2195 * fac01_var_2141(jl_var_2220, lay_var_2172)
        fac101_var_2188 = fk1_var_2196 * fac01_var_2141(jl_var_2220, lay_var_2172)
        fac201_var_2189 = fk2_var_2197 * fac01_var_2141(jl_var_2220, lay_var_2172)
        fac011_var_2190 = fk0_var_2195 * fac11_var_2143(jl_var_2220, lay_var_2172)
        fac111_var_2191 = fk1_var_2196 * fac11_var_2143(jl_var_2220, lay_var_2172)
        fac211_var_2192 = fk2_var_2197 * fac11_var_2143(jl_var_2220, lay_var_2172)
      ELSE
        fac001_var_2187 = (1.0D0 - fs1_var_2202) * fac01_var_2141(jl_var_2220, lay_var_2172)
        fac011_var_2190 = (1.0D0 - fs1_var_2202) * fac11_var_2143(jl_var_2220, lay_var_2172)
        fac101_var_2188 = fs1_var_2202 * fac01_var_2141(jl_var_2220, lay_var_2172)
        fac111_var_2191 = fs1_var_2202 * fac11_var_2143(jl_var_2220, lay_var_2172)
        fac201_var_2189 = 0.0D0
        fac211_var_2192 = 0.0D0
      END IF
      IF (specparm_var_2200 .LT. 0.125D0) THEN
        tau_major_var_2177(1 : ng15) = speccomb_var_2201 * (fac000_var_2181 * absa_var_164(ind0_var_2164, 1 : 2) + fac100_var_2182 * absa_var_164(ind0_var_2164 + 1, 1 : 2) + fac200_var_2183 * absa_var_164(ind0_var_2164 + 2, 1 : 2) + fac010_var_2184 * absa_var_164(ind0_var_2164 + 9, 1 : 2) + fac110_var_2185 * absa_var_164(ind0_var_2164 + 10, 1 : 2) + fac210_var_2186 * absa_var_164(ind0_var_2164 + 11, 1 : 2))
      ELSE IF (specparm_var_2200 .GT. 0.875D0) THEN
        tau_major_var_2177(1 : ng15) = speccomb_var_2201 * (fac200_var_2183 * absa_var_164(ind0_var_2164 - 1, 1 : 2) + fac100_var_2182 * absa_var_164(ind0_var_2164, 1 : 2) + fac000_var_2181 * absa_var_164(ind0_var_2164 + 1, 1 : 2) + fac210_var_2186 * absa_var_164(ind0_var_2164 + 8, 1 : 2) + fac110_var_2185 * absa_var_164(ind0_var_2164 + 9, 1 : 2) + fac010_var_2184 * absa_var_164(ind0_var_2164 + 10, 1 : 2))
      ELSE
        tau_major_var_2177(1 : ng15) = speccomb_var_2201 * (fac000_var_2181 * absa_var_164(ind0_var_2164, 1 : 2) + fac100_var_2182 * absa_var_164(ind0_var_2164 + 1, 1 : 2) + fac010_var_2184 * absa_var_164(ind0_var_2164 + 9, 1 : 2) + fac110_var_2185 * absa_var_164(ind0_var_2164 + 10, 1 : 2))
      END IF
      IF (specparm1_var_2204 .LT. 0.125D0) THEN
        tau_major1_var_2178(1 : ng15) = speccomb1_var_2205 * (fac001_var_2187 * absa_var_164(ind1_var_2165, 1 : 2) + fac101_var_2188 * absa_var_164(ind1_var_2165 + 1, 1 : 2) + fac201_var_2189 * absa_var_164(ind1_var_2165 + 2, 1 : 2) + fac011_var_2190 * absa_var_164(ind1_var_2165 + 9, 1 : 2) + fac111_var_2191 * absa_var_164(ind1_var_2165 + 10, 1 : 2) + fac211_var_2192 * absa_var_164(ind1_var_2165 + 11, 1 : 2))
      ELSE IF (specparm1_var_2204 .GT. 0.875D0) THEN
        tau_major1_var_2178(1 : ng15) = speccomb1_var_2205 * (fac201_var_2189 * absa_var_164(ind1_var_2165 - 1, 1 : 2) + fac101_var_2188 * absa_var_164(ind1_var_2165, 1 : 2) + fac001_var_2187 * absa_var_164(ind1_var_2165 + 1, 1 : 2) + fac211_var_2192 * absa_var_164(ind1_var_2165 + 8, 1 : 2) + fac111_var_2191 * absa_var_164(ind1_var_2165 + 9, 1 : 2) + fac011_var_2190 * absa_var_164(ind1_var_2165 + 10, 1 : 2))
      ELSE
        tau_major1_var_2178(1 : ng15) = speccomb1_var_2205 * (fac001_var_2187 * absa_var_164(ind1_var_2165, 1 : 2) + fac101_var_2188 * absa_var_164(ind1_var_2165 + 1, 1 : 2) + fac011_var_2190 * absa_var_164(ind1_var_2165 + 9, 1 : 2) + fac111_var_2191 * absa_var_164(ind1_var_2165 + 10, 1 : 2))
      END IF
      DO ig_var_2163 = 1, 2
        tauself_var_2176 = selffac_var_2152(jl_var_2220, lay_var_2172) * (selfref_var_166(inds_var_2166, ig_var_2163) + selffrac_var_2153(jl_var_2220, lay_var_2172) * (selfref_var_166(inds_var_2166 + 1, ig_var_2163) - selfref_var_166(inds_var_2166, ig_var_2163)))
        taufor_var_2175 = forfac_var_2157(jl_var_2220, lay_var_2172) * (forref_var_167(indf_var_2167, ig_var_2163) + forfrac_var_2158(jl_var_2220, lay_var_2172) * (forref_var_167(indf_var_2167 + 1, ig_var_2163) - forref_var_167(indf_var_2167, ig_var_2163)))
        n2m1 = ka_mn2_var_165(jmn2, indm_var_2168, ig_var_2163) + fmn2 * (ka_mn2_var_165(jmn2 + 1, indm_var_2168, ig_var_2163) - ka_mn2_var_165(jmn2, indm_var_2168, ig_var_2163))
        n2m2 = ka_mn2_var_165(jmn2, indm_var_2168 + 1, ig_var_2163) + fmn2 * (ka_mn2_var_165(jmn2 + 1, indm_var_2168 + 1, ig_var_2163) - ka_mn2_var_165(jmn2, indm_var_2168 + 1, ig_var_2163))
        taun2_var_2179 = scalen2_var_2180 * (n2m1 + minorfrac_var_2159(jl_var_2220, lay_var_2172) * (n2m2 - n2m1))
        taug_var_2138(jl_var_2220, 136 + ig_var_2163, lay_var_2172) = tau_major_var_2177(ig_var_2163) + tau_major1_var_2178(ig_var_2163) + tauself_var_2176 + taufor_var_2175 + taun2_var_2179
        fracs_var_2155(jl_var_2220, 136 + ig_var_2163, lay_var_2172) = fracrefa_var_163(ig_var_2163, jpl_var_2171) + fpl_var_2206 * (fracrefa_var_163(ig_var_2163, jpl_var_2171 + 1) - fracrefa_var_163(ig_var_2163, jpl_var_2171))
      END DO
    END DO
  END DO
  DO ig_var_2163 = 1, 2
    DO lay_var_2172 = laytrop_max_var_2211 + 1, klev_var_2137
      DO jl_var_2220 = kidia_var_2135, kfdia_var_2136
        taug_var_2138(jl_var_2220, 136 + ig_var_2163, lay_var_2172) = 0.0D0
        fracs_var_2155(jl_var_2220, 136 + ig_var_2163, lay_var_2172) = 0.0D0
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2211 /= laytrop_min_var_2210) THEN
    DO lay_var_2172 = laytrop_min_var_2210 + 1, laytrop_max_var_2211
      ixc0_var_2217 = ixc_var_2212(lay_var_2172)
      DO ixp_var_2218 = 1, ixc0_var_2217
        jl_var_2220 = ixlow_var_2213(ixp_var_2218, lay_var_2172)
        speccomb_var_2201 = coln2o_var_2150(jl_var_2220, lay_var_2172) + rat_n2oco2(jl_var_2220, lay_var_2172) * colco2_var_2149(jl_var_2220, lay_var_2172)
        specparm_var_2200 = MIN(coln2o_var_2150(jl_var_2220, lay_var_2172) / speccomb_var_2201, oneminus_var_2147)
        specmult_var_2199 = 8.0D0 * (specparm_var_2200)
        js_var_2169 = 1 + INT(specmult_var_2199)
        fs_var_2198 = ((specmult_var_2199) - AINT((specmult_var_2199)))
        speccomb1_var_2205 = coln2o_var_2150(jl_var_2220, lay_var_2172) + rat_n2oco2_1(jl_var_2220, lay_var_2172) * colco2_var_2149(jl_var_2220, lay_var_2172)
        specparm1_var_2204 = MIN(coln2o_var_2150(jl_var_2220, lay_var_2172) / speccomb1_var_2205, oneminus_var_2147)
        specmult1_var_2203 = 8.0D0 * (specparm1_var_2204)
        js1_var_2170 = 1 + INT(specmult1_var_2203)
        fs1_var_2202 = ((specmult1_var_2203) - AINT((specmult1_var_2203)))
        speccomb_mn2 = coln2o_var_2150(jl_var_2220, lay_var_2172) + refrat_m_a_var_2174 * colco2_var_2149(jl_var_2220, lay_var_2172)
        specparm_mn2 = MIN(coln2o_var_2150(jl_var_2220, lay_var_2172) / speccomb_mn2, oneminus_var_2147)
        specmult_mn2 = 8.0D0 * specparm_mn2
        jmn2 = 1 + INT(specmult_mn2)
        fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))
        speccomb_planck_var_2209 = coln2o_var_2150(jl_var_2220, lay_var_2172) + refrat_planck_a_var_2173 * colco2_var_2149(jl_var_2220, lay_var_2172)
        specparm_planck_var_2208 = MIN(coln2o_var_2150(jl_var_2220, lay_var_2172) / speccomb_planck_var_2209, oneminus_var_2147)
        specmult_planck_var_2207 = 8.0D0 * specparm_planck_var_2208
        jpl_var_2171 = 1 + INT(specmult_planck_var_2207)
        fpl_var_2206 = ((specmult_planck_var_2207) - AINT((specmult_planck_var_2207)))
        ind0_var_2164 = ((jp_var_2144(jl_var_2220, lay_var_2172) - 1) * 5 + (jt_var_2145(jl_var_2220, lay_var_2172) - 1)) * nspa_var_237(15) + js_var_2169
        ind1_var_2165 = (jp_var_2144(jl_var_2220, lay_var_2172) * 5 + (jt1_var_2146(jl_var_2220, lay_var_2172) - 1)) * nspa_var_237(15) + js1_var_2170
        inds_var_2166 = indself_var_2154(jl_var_2220, lay_var_2172)
        indf_var_2167 = indfor_var_2156(jl_var_2220, lay_var_2172)
        indm_var_2168 = indminor_var_2160(jl_var_2220, lay_var_2172)
        scalen2_var_2180 = colbrd_var_2162(jl_var_2220, lay_var_2172) * scaleminor_var_2161(jl_var_2220, lay_var_2172)
        IF (specparm_var_2200 .LT. 0.125D0) THEN
          p_var_2193 = fs_var_2198 - 1.0D0
          p4_var_2194 = p_var_2193 ** 4
          fk0_var_2195 = p4_var_2194
          fk1_var_2196 = 1.0D0 - p_var_2193 - 2.0D0 * p4_var_2194
          fk2_var_2197 = p_var_2193 + p4_var_2194
          fac000_var_2181 = fk0_var_2195 * fac00_var_2140(jl_var_2220, lay_var_2172)
          fac100_var_2182 = fk1_var_2196 * fac00_var_2140(jl_var_2220, lay_var_2172)
          fac200_var_2183 = fk2_var_2197 * fac00_var_2140(jl_var_2220, lay_var_2172)
          fac010_var_2184 = fk0_var_2195 * fac10_var_2142(jl_var_2220, lay_var_2172)
          fac110_var_2185 = fk1_var_2196 * fac10_var_2142(jl_var_2220, lay_var_2172)
          fac210_var_2186 = fk2_var_2197 * fac10_var_2142(jl_var_2220, lay_var_2172)
        ELSE IF (specparm_var_2200 .GT. 0.875D0) THEN
          p_var_2193 = - fs_var_2198
          p4_var_2194 = p_var_2193 ** 4
          fk0_var_2195 = p4_var_2194
          fk1_var_2196 = 1.0D0 - p_var_2193 - 2.0D0 * p4_var_2194
          fk2_var_2197 = p_var_2193 + p4_var_2194
          fac000_var_2181 = fk0_var_2195 * fac00_var_2140(jl_var_2220, lay_var_2172)
          fac100_var_2182 = fk1_var_2196 * fac00_var_2140(jl_var_2220, lay_var_2172)
          fac200_var_2183 = fk2_var_2197 * fac00_var_2140(jl_var_2220, lay_var_2172)
          fac010_var_2184 = fk0_var_2195 * fac10_var_2142(jl_var_2220, lay_var_2172)
          fac110_var_2185 = fk1_var_2196 * fac10_var_2142(jl_var_2220, lay_var_2172)
          fac210_var_2186 = fk2_var_2197 * fac10_var_2142(jl_var_2220, lay_var_2172)
        ELSE
          fac000_var_2181 = (1.0D0 - fs_var_2198) * fac00_var_2140(jl_var_2220, lay_var_2172)
          fac010_var_2184 = (1.0D0 - fs_var_2198) * fac10_var_2142(jl_var_2220, lay_var_2172)
          fac100_var_2182 = fs_var_2198 * fac00_var_2140(jl_var_2220, lay_var_2172)
          fac110_var_2185 = fs_var_2198 * fac10_var_2142(jl_var_2220, lay_var_2172)
          fac200_var_2183 = 0.0D0
          fac210_var_2186 = 0.0D0
        END IF
        IF (specparm1_var_2204 .LT. 0.125D0) THEN
          p_var_2193 = fs1_var_2202 - 1.0D0
          p4_var_2194 = p_var_2193 ** 4
          fk0_var_2195 = p4_var_2194
          fk1_var_2196 = 1.0D0 - p_var_2193 - 2.0D0 * p4_var_2194
          fk2_var_2197 = p_var_2193 + p4_var_2194
          fac001_var_2187 = fk0_var_2195 * fac01_var_2141(jl_var_2220, lay_var_2172)
          fac101_var_2188 = fk1_var_2196 * fac01_var_2141(jl_var_2220, lay_var_2172)
          fac201_var_2189 = fk2_var_2197 * fac01_var_2141(jl_var_2220, lay_var_2172)
          fac011_var_2190 = fk0_var_2195 * fac11_var_2143(jl_var_2220, lay_var_2172)
          fac111_var_2191 = fk1_var_2196 * fac11_var_2143(jl_var_2220, lay_var_2172)
          fac211_var_2192 = fk2_var_2197 * fac11_var_2143(jl_var_2220, lay_var_2172)
        ELSE IF (specparm1_var_2204 .GT. 0.875D0) THEN
          p_var_2193 = - fs1_var_2202
          p4_var_2194 = p_var_2193 ** 4
          fk0_var_2195 = p4_var_2194
          fk1_var_2196 = 1.0D0 - p_var_2193 - 2.0D0 * p4_var_2194
          fk2_var_2197 = p_var_2193 + p4_var_2194
          fac001_var_2187 = fk0_var_2195 * fac01_var_2141(jl_var_2220, lay_var_2172)
          fac101_var_2188 = fk1_var_2196 * fac01_var_2141(jl_var_2220, lay_var_2172)
          fac201_var_2189 = fk2_var_2197 * fac01_var_2141(jl_var_2220, lay_var_2172)
          fac011_var_2190 = fk0_var_2195 * fac11_var_2143(jl_var_2220, lay_var_2172)
          fac111_var_2191 = fk1_var_2196 * fac11_var_2143(jl_var_2220, lay_var_2172)
          fac211_var_2192 = fk2_var_2197 * fac11_var_2143(jl_var_2220, lay_var_2172)
        ELSE
          fac001_var_2187 = (1.0D0 - fs1_var_2202) * fac01_var_2141(jl_var_2220, lay_var_2172)
          fac011_var_2190 = (1.0D0 - fs1_var_2202) * fac11_var_2143(jl_var_2220, lay_var_2172)
          fac101_var_2188 = fs1_var_2202 * fac01_var_2141(jl_var_2220, lay_var_2172)
          fac111_var_2191 = fs1_var_2202 * fac11_var_2143(jl_var_2220, lay_var_2172)
          fac201_var_2189 = 0.0D0
          fac211_var_2192 = 0.0D0
        END IF
        IF (specparm_var_2200 .LT. 0.125D0) THEN
          tau_major_var_2177(1 : ng15) = speccomb_var_2201 * (fac000_var_2181 * absa_var_164(ind0_var_2164, 1 : 2) + fac100_var_2182 * absa_var_164(ind0_var_2164 + 1, 1 : 2) + fac200_var_2183 * absa_var_164(ind0_var_2164 + 2, 1 : 2) + fac010_var_2184 * absa_var_164(ind0_var_2164 + 9, 1 : 2) + fac110_var_2185 * absa_var_164(ind0_var_2164 + 10, 1 : 2) + fac210_var_2186 * absa_var_164(ind0_var_2164 + 11, 1 : 2))
        ELSE IF (specparm_var_2200 .GT. 0.875D0) THEN
          tau_major_var_2177(1 : ng15) = speccomb_var_2201 * (fac200_var_2183 * absa_var_164(ind0_var_2164 - 1, 1 : 2) + fac100_var_2182 * absa_var_164(ind0_var_2164, 1 : 2) + fac000_var_2181 * absa_var_164(ind0_var_2164 + 1, 1 : 2) + fac210_var_2186 * absa_var_164(ind0_var_2164 + 8, 1 : 2) + fac110_var_2185 * absa_var_164(ind0_var_2164 + 9, 1 : 2) + fac010_var_2184 * absa_var_164(ind0_var_2164 + 10, 1 : 2))
        ELSE
          tau_major_var_2177(1 : ng15) = speccomb_var_2201 * (fac000_var_2181 * absa_var_164(ind0_var_2164, 1 : 2) + fac100_var_2182 * absa_var_164(ind0_var_2164 + 1, 1 : 2) + fac010_var_2184 * absa_var_164(ind0_var_2164 + 9, 1 : 2) + fac110_var_2185 * absa_var_164(ind0_var_2164 + 10, 1 : 2))
        END IF
        IF (specparm1_var_2204 .LT. 0.125D0) THEN
          tau_major1_var_2178(1 : ng15) = speccomb1_var_2205 * (fac001_var_2187 * absa_var_164(ind1_var_2165, 1 : 2) + fac101_var_2188 * absa_var_164(ind1_var_2165 + 1, 1 : 2) + fac201_var_2189 * absa_var_164(ind1_var_2165 + 2, 1 : 2) + fac011_var_2190 * absa_var_164(ind1_var_2165 + 9, 1 : 2) + fac111_var_2191 * absa_var_164(ind1_var_2165 + 10, 1 : 2) + fac211_var_2192 * absa_var_164(ind1_var_2165 + 11, 1 : 2))
        ELSE IF (specparm1_var_2204 .GT. 0.875D0) THEN
          tau_major1_var_2178(1 : ng15) = speccomb1_var_2205 * (fac201_var_2189 * absa_var_164(ind1_var_2165 - 1, 1 : 2) + fac101_var_2188 * absa_var_164(ind1_var_2165, 1 : 2) + fac001_var_2187 * absa_var_164(ind1_var_2165 + 1, 1 : 2) + fac211_var_2192 * absa_var_164(ind1_var_2165 + 8, 1 : 2) + fac111_var_2191 * absa_var_164(ind1_var_2165 + 9, 1 : 2) + fac011_var_2190 * absa_var_164(ind1_var_2165 + 10, 1 : 2))
        ELSE
          tau_major1_var_2178(1 : ng15) = speccomb1_var_2205 * (fac001_var_2187 * absa_var_164(ind1_var_2165, 1 : 2) + fac101_var_2188 * absa_var_164(ind1_var_2165 + 1, 1 : 2) + fac011_var_2190 * absa_var_164(ind1_var_2165 + 9, 1 : 2) + fac111_var_2191 * absa_var_164(ind1_var_2165 + 10, 1 : 2))
        END IF
        DO ig_var_2163 = 1, 2
          tauself_var_2176 = selffac_var_2152(jl_var_2220, lay_var_2172) * (selfref_var_166(inds_var_2166, ig_var_2163) + selffrac_var_2153(jl_var_2220, lay_var_2172) * (selfref_var_166(inds_var_2166 + 1, ig_var_2163) - selfref_var_166(inds_var_2166, ig_var_2163)))
          taufor_var_2175 = forfac_var_2157(jl_var_2220, lay_var_2172) * (forref_var_167(indf_var_2167, ig_var_2163) + forfrac_var_2158(jl_var_2220, lay_var_2172) * (forref_var_167(indf_var_2167 + 1, ig_var_2163) - forref_var_167(indf_var_2167, ig_var_2163)))
          n2m1 = ka_mn2_var_165(jmn2, indm_var_2168, ig_var_2163) + fmn2 * (ka_mn2_var_165(jmn2 + 1, indm_var_2168, ig_var_2163) - ka_mn2_var_165(jmn2, indm_var_2168, ig_var_2163))
          n2m2 = ka_mn2_var_165(jmn2, indm_var_2168 + 1, ig_var_2163) + fmn2 * (ka_mn2_var_165(jmn2 + 1, indm_var_2168 + 1, ig_var_2163) - ka_mn2_var_165(jmn2, indm_var_2168 + 1, ig_var_2163))
          taun2_var_2179 = scalen2_var_2180 * (n2m1 + minorfrac_var_2159(jl_var_2220, lay_var_2172) * (n2m2 - n2m1))
          taug_var_2138(jl_var_2220, 136 + ig_var_2163, lay_var_2172) = tau_major_var_2177(ig_var_2163) + tau_major1_var_2178(ig_var_2163) + tauself_var_2176 + taufor_var_2175 + taun2_var_2179
          fracs_var_2155(jl_var_2220, 136 + ig_var_2163, lay_var_2172) = fracrefa_var_163(ig_var_2163, jpl_var_2171) + fpl_var_2206 * (fracrefa_var_163(ig_var_2163, jpl_var_2171 + 1) - fracrefa_var_163(ig_var_2163, jpl_var_2171))
        END DO
      END DO
      ixc0_var_2217 = kfdia_var_2136 - kidia_var_2135 + 1 - ixc0_var_2217
      DO ig_var_2163 = 1, 2
        DO ixp_var_2218 = 1, ixc0_var_2217
          jl_var_2220 = ixhigh_var_2214(ixp_var_2218, lay_var_2172)
          taug_var_2138(jl_var_2220, 136 + ig_var_2163, lay_var_2172) = 0.0D0
          fracs_var_2155(jl_var_2220, 136 + ig_var_2163, lay_var_2172) = 0.0D0
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol15
SUBROUTINE rrtm_taumol16(kidia_var_2221, kfdia_var_2222, klev_var_2223, taug_var_2224, p_tauaerl_var_2225, fac00_var_2226, fac01_var_2227, fac10_var_2228, fac11_var_2229, forfac_var_2244, forfrac_var_2245, indfor_var_2243, jp_var_2230, jt_var_2231, jt1_var_2232, oneminus_var_2233, colh2o_var_2234, colch4_var_2235, laytrop_var_2236, selffac_var_2237, selffrac_var_2238, indself_var_2239, fracs_var_2240, rat_h2och4_var_2241, rat_h2och4_1_var_2242)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_237, nspb_var_238
  USE yoerrtm, ONLY: ng16
  USE yoerrta16, ONLY: absa_var_170, absb_var_171, forref_var_173, fracrefa_var_168, fracrefb_var_169, selfref_var_172
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2221
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2222
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2223
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2224(kidia_var_2221 : kfdia_var_2222, 140, klev_var_2223)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2225(kidia_var_2221 : kfdia_var_2222, klev_var_2223, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2226(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2227(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2228(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2229(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2230(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2231(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2232(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2233
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2234(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  REAL(KIND = 8), INTENT(IN) :: colch4_var_2235(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2236(kidia_var_2221 : kfdia_var_2222)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2237(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2238(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2239(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2240(kidia_var_2221 : kfdia_var_2222, 140, klev_var_2223)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_var_2241(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  REAL(KIND = 8), INTENT(IN) :: rat_h2och4_1_var_2242(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2243(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2244(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2245(kidia_var_2221 : kfdia_var_2222, klev_var_2223)
  INTEGER(KIND = 4) :: ig_var_2246, ind0_var_2247, ind1_var_2248, inds_var_2249, indf_var_2250, js_var_2251, js1_var_2252, jpl_var_2253, lay_var_2254
  REAL(KIND = 8) :: fac000_var_2255, fac100_var_2256, fac200_var_2257, fac010_var_2258, fac110_var_2259, fac210_var_2260, fac001_var_2261, fac101_var_2262, fac201_var_2263, fac011_var_2264, fac111_var_2265, fac211_var_2266
  REAL(KIND = 8) :: p_var_2267, p4_var_2268, fk0_var_2269, fk1_var_2270, fk2_var_2271
  REAL(KIND = 8) :: refrat_planck_a_var_2272
  REAL(KIND = 8) :: taufor_var_2273, tauself_var_2274, tau_major_var_2275(2), tau_major1_var_2276(2)
  REAL(KIND = 8) :: fs_var_2277, specmult_var_2278, specparm_var_2279, speccomb_var_2280, fs1_var_2281, specmult1_var_2282, specparm1_var_2283, speccomb1_var_2284, fpl_var_2285, specmult_planck_var_2286, specparm_planck_var_2287, speccomb_planck_var_2288
  INTEGER(KIND = 4) :: laytrop_min_var_2289, laytrop_max_var_2290
  INTEGER(KIND = 4) :: ixc_var_2291(klev_var_2223), ixlow_var_2292(kfdia_var_2222, klev_var_2223), ixhigh_var_2293(kfdia_var_2222, klev_var_2223)
  INTEGER(KIND = 4) :: ich_var_2294, icl_var_2295, ixc0_var_2296, ixp_var_2297, jc_var_2298, jl_var_2299
  laytrop_min_var_2289 = MINVAL(laytrop_var_2236)
  laytrop_max_var_2290 = MAXVAL(laytrop_var_2236)
  ixlow_var_2292 = 0
  ixhigh_var_2293 = 0
  ixc_var_2291 = 0
  DO lay_var_2254 = laytrop_min_var_2289 + 1, laytrop_max_var_2290
    icl_var_2295 = 0
    ich_var_2294 = 0
    DO jc_var_2298 = kidia_var_2221, kfdia_var_2222
      IF (lay_var_2254 <= laytrop_var_2236(jc_var_2298)) THEN
        icl_var_2295 = icl_var_2295 + 1
        ixlow_var_2292(icl_var_2295, lay_var_2254) = jc_var_2298
      ELSE
        ich_var_2294 = ich_var_2294 + 1
        ixhigh_var_2293(ich_var_2294, lay_var_2254) = jc_var_2298
      END IF
    END DO
    ixc_var_2291(lay_var_2254) = icl_var_2295
  END DO
  refrat_planck_a_var_2272 = chi_mls(1, 6) / chi_mls(6, 6)
  DO lay_var_2254 = 1, laytrop_min_var_2289
    DO jl_var_2299 = kidia_var_2221, kfdia_var_2222
      speccomb_var_2280 = colh2o_var_2234(jl_var_2299, lay_var_2254) + rat_h2och4_var_2241(jl_var_2299, lay_var_2254) * colch4_var_2235(jl_var_2299, lay_var_2254)
      specparm_var_2279 = MIN(colh2o_var_2234(jl_var_2299, lay_var_2254) / speccomb_var_2280, oneminus_var_2233)
      specmult_var_2278 = 8.0D0 * (specparm_var_2279)
      js_var_2251 = 1 + INT(specmult_var_2278)
      fs_var_2277 = ((specmult_var_2278) - AINT((specmult_var_2278)))
      speccomb1_var_2284 = colh2o_var_2234(jl_var_2299, lay_var_2254) + rat_h2och4_1_var_2242(jl_var_2299, lay_var_2254) * colch4_var_2235(jl_var_2299, lay_var_2254)
      specparm1_var_2283 = MIN(colh2o_var_2234(jl_var_2299, lay_var_2254) / speccomb1_var_2284, oneminus_var_2233)
      specmult1_var_2282 = 8.0D0 * (specparm1_var_2283)
      js1_var_2252 = 1 + INT(specmult1_var_2282)
      fs1_var_2281 = ((specmult1_var_2282) - AINT((specmult1_var_2282)))
      speccomb_planck_var_2288 = colh2o_var_2234(jl_var_2299, lay_var_2254) + refrat_planck_a_var_2272 * colch4_var_2235(jl_var_2299, lay_var_2254)
      specparm_planck_var_2287 = MIN(colh2o_var_2234(jl_var_2299, lay_var_2254) / speccomb_planck_var_2288, oneminus_var_2233)
      specmult_planck_var_2286 = 8.0D0 * specparm_planck_var_2287
      jpl_var_2253 = 1 + INT(specmult_planck_var_2286)
      fpl_var_2285 = ((specmult_planck_var_2286) - AINT((specmult_planck_var_2286)))
      ind0_var_2247 = ((jp_var_2230(jl_var_2299, lay_var_2254) - 1) * 5 + (jt_var_2231(jl_var_2299, lay_var_2254) - 1)) * nspa_var_237(16) + js_var_2251
      ind1_var_2248 = (jp_var_2230(jl_var_2299, lay_var_2254) * 5 + (jt1_var_2232(jl_var_2299, lay_var_2254) - 1)) * nspa_var_237(16) + js1_var_2252
      inds_var_2249 = indself_var_2239(jl_var_2299, lay_var_2254)
      indf_var_2250 = indfor_var_2243(jl_var_2299, lay_var_2254)
      IF (specparm_var_2279 .LT. 0.125D0) THEN
        p_var_2267 = fs_var_2277 - 1.0D0
        p4_var_2268 = p_var_2267 ** 4
        fk0_var_2269 = p4_var_2268
        fk1_var_2270 = 1.0D0 - p_var_2267 - 2.0D0 * p4_var_2268
        fk2_var_2271 = p_var_2267 + p4_var_2268
        fac000_var_2255 = fk0_var_2269 * fac00_var_2226(jl_var_2299, lay_var_2254)
        fac100_var_2256 = fk1_var_2270 * fac00_var_2226(jl_var_2299, lay_var_2254)
        fac200_var_2257 = fk2_var_2271 * fac00_var_2226(jl_var_2299, lay_var_2254)
        fac010_var_2258 = fk0_var_2269 * fac10_var_2228(jl_var_2299, lay_var_2254)
        fac110_var_2259 = fk1_var_2270 * fac10_var_2228(jl_var_2299, lay_var_2254)
        fac210_var_2260 = fk2_var_2271 * fac10_var_2228(jl_var_2299, lay_var_2254)
      ELSE IF (specparm_var_2279 .GT. 0.875D0) THEN
        p_var_2267 = - fs_var_2277
        p4_var_2268 = p_var_2267 ** 4
        fk0_var_2269 = p4_var_2268
        fk1_var_2270 = 1.0D0 - p_var_2267 - 2.0D0 * p4_var_2268
        fk2_var_2271 = p_var_2267 + p4_var_2268
        fac000_var_2255 = fk0_var_2269 * fac00_var_2226(jl_var_2299, lay_var_2254)
        fac100_var_2256 = fk1_var_2270 * fac00_var_2226(jl_var_2299, lay_var_2254)
        fac200_var_2257 = fk2_var_2271 * fac00_var_2226(jl_var_2299, lay_var_2254)
        fac010_var_2258 = fk0_var_2269 * fac10_var_2228(jl_var_2299, lay_var_2254)
        fac110_var_2259 = fk1_var_2270 * fac10_var_2228(jl_var_2299, lay_var_2254)
        fac210_var_2260 = fk2_var_2271 * fac10_var_2228(jl_var_2299, lay_var_2254)
      ELSE
        fac000_var_2255 = (1.0D0 - fs_var_2277) * fac00_var_2226(jl_var_2299, lay_var_2254)
        fac010_var_2258 = (1.0D0 - fs_var_2277) * fac10_var_2228(jl_var_2299, lay_var_2254)
        fac100_var_2256 = fs_var_2277 * fac00_var_2226(jl_var_2299, lay_var_2254)
        fac110_var_2259 = fs_var_2277 * fac10_var_2228(jl_var_2299, lay_var_2254)
        fac200_var_2257 = 0.0D0
        fac210_var_2260 = 0.0D0
      END IF
      IF (specparm1_var_2283 .LT. 0.125D0) THEN
        p_var_2267 = fs1_var_2281 - 1.0D0
        p4_var_2268 = p_var_2267 ** 4
        fk0_var_2269 = p4_var_2268
        fk1_var_2270 = 1.0D0 - p_var_2267 - 2.0D0 * p4_var_2268
        fk2_var_2271 = p_var_2267 + p4_var_2268
        fac001_var_2261 = fk0_var_2269 * fac01_var_2227(jl_var_2299, lay_var_2254)
        fac101_var_2262 = fk1_var_2270 * fac01_var_2227(jl_var_2299, lay_var_2254)
        fac201_var_2263 = fk2_var_2271 * fac01_var_2227(jl_var_2299, lay_var_2254)
        fac011_var_2264 = fk0_var_2269 * fac11_var_2229(jl_var_2299, lay_var_2254)
        fac111_var_2265 = fk1_var_2270 * fac11_var_2229(jl_var_2299, lay_var_2254)
        fac211_var_2266 = fk2_var_2271 * fac11_var_2229(jl_var_2299, lay_var_2254)
      ELSE IF (specparm1_var_2283 .GT. 0.875D0) THEN
        p_var_2267 = - fs1_var_2281
        p4_var_2268 = p_var_2267 ** 4
        fk0_var_2269 = p4_var_2268
        fk1_var_2270 = 1.0D0 - p_var_2267 - 2.0D0 * p4_var_2268
        fk2_var_2271 = p_var_2267 + p4_var_2268
        fac001_var_2261 = fk0_var_2269 * fac01_var_2227(jl_var_2299, lay_var_2254)
        fac101_var_2262 = fk1_var_2270 * fac01_var_2227(jl_var_2299, lay_var_2254)
        fac201_var_2263 = fk2_var_2271 * fac01_var_2227(jl_var_2299, lay_var_2254)
        fac011_var_2264 = fk0_var_2269 * fac11_var_2229(jl_var_2299, lay_var_2254)
        fac111_var_2265 = fk1_var_2270 * fac11_var_2229(jl_var_2299, lay_var_2254)
        fac211_var_2266 = fk2_var_2271 * fac11_var_2229(jl_var_2299, lay_var_2254)
      ELSE
        fac001_var_2261 = (1.0D0 - fs1_var_2281) * fac01_var_2227(jl_var_2299, lay_var_2254)
        fac011_var_2264 = (1.0D0 - fs1_var_2281) * fac11_var_2229(jl_var_2299, lay_var_2254)
        fac101_var_2262 = fs1_var_2281 * fac01_var_2227(jl_var_2299, lay_var_2254)
        fac111_var_2265 = fs1_var_2281 * fac11_var_2229(jl_var_2299, lay_var_2254)
        fac201_var_2263 = 0.0D0
        fac211_var_2266 = 0.0D0
      END IF
      IF (specparm_var_2279 .LT. 0.125D0) THEN
        tau_major_var_2275(1 : ng16) = speccomb_var_2280 * (fac000_var_2255 * absa_var_170(ind0_var_2247, 1 : 2) + fac100_var_2256 * absa_var_170(ind0_var_2247 + 1, 1 : 2) + fac200_var_2257 * absa_var_170(ind0_var_2247 + 2, 1 : 2) + fac010_var_2258 * absa_var_170(ind0_var_2247 + 9, 1 : 2) + fac110_var_2259 * absa_var_170(ind0_var_2247 + 10, 1 : 2) + fac210_var_2260 * absa_var_170(ind0_var_2247 + 11, 1 : 2))
      ELSE IF (specparm_var_2279 .GT. 0.875D0) THEN
        tau_major_var_2275(1 : ng16) = speccomb_var_2280 * (fac200_var_2257 * absa_var_170(ind0_var_2247 - 1, 1 : 2) + fac100_var_2256 * absa_var_170(ind0_var_2247, 1 : 2) + fac000_var_2255 * absa_var_170(ind0_var_2247 + 1, 1 : 2) + fac210_var_2260 * absa_var_170(ind0_var_2247 + 8, 1 : 2) + fac110_var_2259 * absa_var_170(ind0_var_2247 + 9, 1 : 2) + fac010_var_2258 * absa_var_170(ind0_var_2247 + 10, 1 : 2))
      ELSE
        tau_major_var_2275(1 : ng16) = speccomb_var_2280 * (fac000_var_2255 * absa_var_170(ind0_var_2247, 1 : 2) + fac100_var_2256 * absa_var_170(ind0_var_2247 + 1, 1 : 2) + fac010_var_2258 * absa_var_170(ind0_var_2247 + 9, 1 : 2) + fac110_var_2259 * absa_var_170(ind0_var_2247 + 10, 1 : 2))
      END IF
      IF (specparm1_var_2283 .LT. 0.125D0) THEN
        tau_major1_var_2276(1 : ng16) = speccomb1_var_2284 * (fac001_var_2261 * absa_var_170(ind1_var_2248, 1 : 2) + fac101_var_2262 * absa_var_170(ind1_var_2248 + 1, 1 : 2) + fac201_var_2263 * absa_var_170(ind1_var_2248 + 2, 1 : 2) + fac011_var_2264 * absa_var_170(ind1_var_2248 + 9, 1 : 2) + fac111_var_2265 * absa_var_170(ind1_var_2248 + 10, 1 : 2) + fac211_var_2266 * absa_var_170(ind1_var_2248 + 11, 1 : 2))
      ELSE IF (specparm1_var_2283 .GT. 0.875D0) THEN
        tau_major1_var_2276(1 : ng16) = speccomb1_var_2284 * (fac201_var_2263 * absa_var_170(ind1_var_2248 - 1, 1 : 2) + fac101_var_2262 * absa_var_170(ind1_var_2248, 1 : 2) + fac001_var_2261 * absa_var_170(ind1_var_2248 + 1, 1 : 2) + fac211_var_2266 * absa_var_170(ind1_var_2248 + 8, 1 : 2) + fac111_var_2265 * absa_var_170(ind1_var_2248 + 9, 1 : 2) + fac011_var_2264 * absa_var_170(ind1_var_2248 + 10, 1 : 2))
      ELSE
        tau_major1_var_2276(1 : ng16) = speccomb1_var_2284 * (fac001_var_2261 * absa_var_170(ind1_var_2248, 1 : 2) + fac101_var_2262 * absa_var_170(ind1_var_2248 + 1, 1 : 2) + fac011_var_2264 * absa_var_170(ind1_var_2248 + 9, 1 : 2) + fac111_var_2265 * absa_var_170(ind1_var_2248 + 10, 1 : 2))
      END IF
      DO ig_var_2246 = 1, 2
        tauself_var_2274 = selffac_var_2237(jl_var_2299, lay_var_2254) * (selfref_var_172(inds_var_2249, ig_var_2246) + selffrac_var_2238(jl_var_2299, lay_var_2254) * (selfref_var_172(inds_var_2249 + 1, ig_var_2246) - selfref_var_172(inds_var_2249, ig_var_2246)))
        taufor_var_2273 = forfac_var_2244(jl_var_2299, lay_var_2254) * (forref_var_173(indf_var_2250, ig_var_2246) + forfrac_var_2245(jl_var_2299, lay_var_2254) * (forref_var_173(indf_var_2250 + 1, ig_var_2246) - forref_var_173(indf_var_2250, ig_var_2246)))
        taug_var_2224(jl_var_2299, 138 + ig_var_2246, lay_var_2254) = tau_major_var_2275(ig_var_2246) + tau_major1_var_2276(ig_var_2246) + tauself_var_2274 + taufor_var_2273
        fracs_var_2240(jl_var_2299, 138 + ig_var_2246, lay_var_2254) = fracrefa_var_168(ig_var_2246, jpl_var_2253) + fpl_var_2285 * (fracrefa_var_168(ig_var_2246, jpl_var_2253 + 1) - fracrefa_var_168(ig_var_2246, jpl_var_2253))
      END DO
    END DO
  END DO
  DO lay_var_2254 = laytrop_max_var_2290 + 1, klev_var_2223
    DO jl_var_2299 = kidia_var_2221, kfdia_var_2222
      ind0_var_2247 = ((jp_var_2230(jl_var_2299, lay_var_2254) - 13) * 5 + (jt_var_2231(jl_var_2299, lay_var_2254) - 1)) * nspb_var_238(16) + 1
      ind1_var_2248 = ((jp_var_2230(jl_var_2299, lay_var_2254) - 12) * 5 + (jt1_var_2232(jl_var_2299, lay_var_2254) - 1)) * nspb_var_238(16) + 1
      DO ig_var_2246 = 1, 2
        taug_var_2224(jl_var_2299, 138 + ig_var_2246, lay_var_2254) = colch4_var_2235(jl_var_2299, lay_var_2254) * (fac00_var_2226(jl_var_2299, lay_var_2254) * absb_var_171(ind0_var_2247, ig_var_2246) + fac10_var_2228(jl_var_2299, lay_var_2254) * absb_var_171(ind0_var_2247 + 1, ig_var_2246) + fac01_var_2227(jl_var_2299, lay_var_2254) * absb_var_171(ind1_var_2248, ig_var_2246) + fac11_var_2229(jl_var_2299, lay_var_2254) * absb_var_171(ind1_var_2248 + 1, ig_var_2246))
        fracs_var_2240(jl_var_2299, 138 + ig_var_2246, lay_var_2254) = fracrefb_var_169(ig_var_2246)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2290 /= laytrop_min_var_2289) THEN
    DO lay_var_2254 = laytrop_min_var_2289 + 1, laytrop_max_var_2290
      ixc0_var_2296 = ixc_var_2291(lay_var_2254)
      DO ixp_var_2297 = 1, ixc0_var_2296
        jl_var_2299 = ixlow_var_2292(ixp_var_2297, lay_var_2254)
        speccomb_var_2280 = colh2o_var_2234(jl_var_2299, lay_var_2254) + rat_h2och4_var_2241(jl_var_2299, lay_var_2254) * colch4_var_2235(jl_var_2299, lay_var_2254)
        specparm_var_2279 = MIN(colh2o_var_2234(jl_var_2299, lay_var_2254) / speccomb_var_2280, oneminus_var_2233)
        specmult_var_2278 = 8.0D0 * (specparm_var_2279)
        js_var_2251 = 1 + INT(specmult_var_2278)
        fs_var_2277 = ((specmult_var_2278) - AINT((specmult_var_2278)))
        speccomb1_var_2284 = colh2o_var_2234(jl_var_2299, lay_var_2254) + rat_h2och4_1_var_2242(jl_var_2299, lay_var_2254) * colch4_var_2235(jl_var_2299, lay_var_2254)
        specparm1_var_2283 = MIN(colh2o_var_2234(jl_var_2299, lay_var_2254) / speccomb1_var_2284, oneminus_var_2233)
        specmult1_var_2282 = 8.0D0 * (specparm1_var_2283)
        js1_var_2252 = 1 + INT(specmult1_var_2282)
        fs1_var_2281 = ((specmult1_var_2282) - AINT((specmult1_var_2282)))
        speccomb_planck_var_2288 = colh2o_var_2234(jl_var_2299, lay_var_2254) + refrat_planck_a_var_2272 * colch4_var_2235(jl_var_2299, lay_var_2254)
        specparm_planck_var_2287 = MIN(colh2o_var_2234(jl_var_2299, lay_var_2254) / speccomb_planck_var_2288, oneminus_var_2233)
        specmult_planck_var_2286 = 8.0D0 * specparm_planck_var_2287
        jpl_var_2253 = 1 + INT(specmult_planck_var_2286)
        fpl_var_2285 = ((specmult_planck_var_2286) - AINT((specmult_planck_var_2286)))
        ind0_var_2247 = ((jp_var_2230(jl_var_2299, lay_var_2254) - 1) * 5 + (jt_var_2231(jl_var_2299, lay_var_2254) - 1)) * nspa_var_237(16) + js_var_2251
        ind1_var_2248 = (jp_var_2230(jl_var_2299, lay_var_2254) * 5 + (jt1_var_2232(jl_var_2299, lay_var_2254) - 1)) * nspa_var_237(16) + js1_var_2252
        inds_var_2249 = indself_var_2239(jl_var_2299, lay_var_2254)
        indf_var_2250 = indfor_var_2243(jl_var_2299, lay_var_2254)
        IF (specparm_var_2279 .LT. 0.125D0) THEN
          p_var_2267 = fs_var_2277 - 1.0D0
          p4_var_2268 = p_var_2267 ** 4
          fk0_var_2269 = p4_var_2268
          fk1_var_2270 = 1.0D0 - p_var_2267 - 2.0D0 * p4_var_2268
          fk2_var_2271 = p_var_2267 + p4_var_2268
          fac000_var_2255 = fk0_var_2269 * fac00_var_2226(jl_var_2299, lay_var_2254)
          fac100_var_2256 = fk1_var_2270 * fac00_var_2226(jl_var_2299, lay_var_2254)
          fac200_var_2257 = fk2_var_2271 * fac00_var_2226(jl_var_2299, lay_var_2254)
          fac010_var_2258 = fk0_var_2269 * fac10_var_2228(jl_var_2299, lay_var_2254)
          fac110_var_2259 = fk1_var_2270 * fac10_var_2228(jl_var_2299, lay_var_2254)
          fac210_var_2260 = fk2_var_2271 * fac10_var_2228(jl_var_2299, lay_var_2254)
        ELSE IF (specparm_var_2279 .GT. 0.875D0) THEN
          p_var_2267 = - fs_var_2277
          p4_var_2268 = p_var_2267 ** 4
          fk0_var_2269 = p4_var_2268
          fk1_var_2270 = 1.0D0 - p_var_2267 - 2.0D0 * p4_var_2268
          fk2_var_2271 = p_var_2267 + p4_var_2268
          fac000_var_2255 = fk0_var_2269 * fac00_var_2226(jl_var_2299, lay_var_2254)
          fac100_var_2256 = fk1_var_2270 * fac00_var_2226(jl_var_2299, lay_var_2254)
          fac200_var_2257 = fk2_var_2271 * fac00_var_2226(jl_var_2299, lay_var_2254)
          fac010_var_2258 = fk0_var_2269 * fac10_var_2228(jl_var_2299, lay_var_2254)
          fac110_var_2259 = fk1_var_2270 * fac10_var_2228(jl_var_2299, lay_var_2254)
          fac210_var_2260 = fk2_var_2271 * fac10_var_2228(jl_var_2299, lay_var_2254)
        ELSE
          fac000_var_2255 = (1.0D0 - fs_var_2277) * fac00_var_2226(jl_var_2299, lay_var_2254)
          fac010_var_2258 = (1.0D0 - fs_var_2277) * fac10_var_2228(jl_var_2299, lay_var_2254)
          fac100_var_2256 = fs_var_2277 * fac00_var_2226(jl_var_2299, lay_var_2254)
          fac110_var_2259 = fs_var_2277 * fac10_var_2228(jl_var_2299, lay_var_2254)
          fac200_var_2257 = 0.0D0
          fac210_var_2260 = 0.0D0
        END IF
        IF (specparm1_var_2283 .LT. 0.125D0) THEN
          p_var_2267 = fs1_var_2281 - 1.0D0
          p4_var_2268 = p_var_2267 ** 4
          fk0_var_2269 = p4_var_2268
          fk1_var_2270 = 1.0D0 - p_var_2267 - 2.0D0 * p4_var_2268
          fk2_var_2271 = p_var_2267 + p4_var_2268
          fac001_var_2261 = fk0_var_2269 * fac01_var_2227(jl_var_2299, lay_var_2254)
          fac101_var_2262 = fk1_var_2270 * fac01_var_2227(jl_var_2299, lay_var_2254)
          fac201_var_2263 = fk2_var_2271 * fac01_var_2227(jl_var_2299, lay_var_2254)
          fac011_var_2264 = fk0_var_2269 * fac11_var_2229(jl_var_2299, lay_var_2254)
          fac111_var_2265 = fk1_var_2270 * fac11_var_2229(jl_var_2299, lay_var_2254)
          fac211_var_2266 = fk2_var_2271 * fac11_var_2229(jl_var_2299, lay_var_2254)
        ELSE IF (specparm1_var_2283 .GT. 0.875D0) THEN
          p_var_2267 = - fs1_var_2281
          p4_var_2268 = p_var_2267 ** 4
          fk0_var_2269 = p4_var_2268
          fk1_var_2270 = 1.0D0 - p_var_2267 - 2.0D0 * p4_var_2268
          fk2_var_2271 = p_var_2267 + p4_var_2268
          fac001_var_2261 = fk0_var_2269 * fac01_var_2227(jl_var_2299, lay_var_2254)
          fac101_var_2262 = fk1_var_2270 * fac01_var_2227(jl_var_2299, lay_var_2254)
          fac201_var_2263 = fk2_var_2271 * fac01_var_2227(jl_var_2299, lay_var_2254)
          fac011_var_2264 = fk0_var_2269 * fac11_var_2229(jl_var_2299, lay_var_2254)
          fac111_var_2265 = fk1_var_2270 * fac11_var_2229(jl_var_2299, lay_var_2254)
          fac211_var_2266 = fk2_var_2271 * fac11_var_2229(jl_var_2299, lay_var_2254)
        ELSE
          fac001_var_2261 = (1.0D0 - fs1_var_2281) * fac01_var_2227(jl_var_2299, lay_var_2254)
          fac011_var_2264 = (1.0D0 - fs1_var_2281) * fac11_var_2229(jl_var_2299, lay_var_2254)
          fac101_var_2262 = fs1_var_2281 * fac01_var_2227(jl_var_2299, lay_var_2254)
          fac111_var_2265 = fs1_var_2281 * fac11_var_2229(jl_var_2299, lay_var_2254)
          fac201_var_2263 = 0.0D0
          fac211_var_2266 = 0.0D0
        END IF
        IF (specparm_var_2279 .LT. 0.125D0) THEN
          tau_major_var_2275(1 : ng16) = speccomb_var_2280 * (fac000_var_2255 * absa_var_170(ind0_var_2247, 1 : 2) + fac100_var_2256 * absa_var_170(ind0_var_2247 + 1, 1 : 2) + fac200_var_2257 * absa_var_170(ind0_var_2247 + 2, 1 : 2) + fac010_var_2258 * absa_var_170(ind0_var_2247 + 9, 1 : 2) + fac110_var_2259 * absa_var_170(ind0_var_2247 + 10, 1 : 2) + fac210_var_2260 * absa_var_170(ind0_var_2247 + 11, 1 : 2))
        ELSE IF (specparm_var_2279 .GT. 0.875D0) THEN
          tau_major_var_2275(1 : ng16) = speccomb_var_2280 * (fac200_var_2257 * absa_var_170(ind0_var_2247 - 1, 1 : 2) + fac100_var_2256 * absa_var_170(ind0_var_2247, 1 : 2) + fac000_var_2255 * absa_var_170(ind0_var_2247 + 1, 1 : 2) + fac210_var_2260 * absa_var_170(ind0_var_2247 + 8, 1 : 2) + fac110_var_2259 * absa_var_170(ind0_var_2247 + 9, 1 : 2) + fac010_var_2258 * absa_var_170(ind0_var_2247 + 10, 1 : 2))
        ELSE
          tau_major_var_2275(1 : ng16) = speccomb_var_2280 * (fac000_var_2255 * absa_var_170(ind0_var_2247, 1 : 2) + fac100_var_2256 * absa_var_170(ind0_var_2247 + 1, 1 : 2) + fac010_var_2258 * absa_var_170(ind0_var_2247 + 9, 1 : 2) + fac110_var_2259 * absa_var_170(ind0_var_2247 + 10, 1 : 2))
        END IF
        IF (specparm1_var_2283 .LT. 0.125D0) THEN
          tau_major1_var_2276(1 : ng16) = speccomb1_var_2284 * (fac001_var_2261 * absa_var_170(ind1_var_2248, 1 : 2) + fac101_var_2262 * absa_var_170(ind1_var_2248 + 1, 1 : 2) + fac201_var_2263 * absa_var_170(ind1_var_2248 + 2, 1 : 2) + fac011_var_2264 * absa_var_170(ind1_var_2248 + 9, 1 : 2) + fac111_var_2265 * absa_var_170(ind1_var_2248 + 10, 1 : 2) + fac211_var_2266 * absa_var_170(ind1_var_2248 + 11, 1 : 2))
        ELSE IF (specparm1_var_2283 .GT. 0.875D0) THEN
          tau_major1_var_2276(1 : ng16) = speccomb1_var_2284 * (fac201_var_2263 * absa_var_170(ind1_var_2248 - 1, 1 : 2) + fac101_var_2262 * absa_var_170(ind1_var_2248, 1 : 2) + fac001_var_2261 * absa_var_170(ind1_var_2248 + 1, 1 : 2) + fac211_var_2266 * absa_var_170(ind1_var_2248 + 8, 1 : 2) + fac111_var_2265 * absa_var_170(ind1_var_2248 + 9, 1 : 2) + fac011_var_2264 * absa_var_170(ind1_var_2248 + 10, 1 : 2))
        ELSE
          tau_major1_var_2276(1 : ng16) = speccomb1_var_2284 * (fac001_var_2261 * absa_var_170(ind1_var_2248, 1 : 2) + fac101_var_2262 * absa_var_170(ind1_var_2248 + 1, 1 : 2) + fac011_var_2264 * absa_var_170(ind1_var_2248 + 9, 1 : 2) + fac111_var_2265 * absa_var_170(ind1_var_2248 + 10, 1 : 2))
        END IF
        DO ig_var_2246 = 1, 2
          tauself_var_2274 = selffac_var_2237(jl_var_2299, lay_var_2254) * (selfref_var_172(inds_var_2249, ig_var_2246) + selffrac_var_2238(jl_var_2299, lay_var_2254) * (selfref_var_172(inds_var_2249 + 1, ig_var_2246) - selfref_var_172(inds_var_2249, ig_var_2246)))
          taufor_var_2273 = forfac_var_2244(jl_var_2299, lay_var_2254) * (forref_var_173(indf_var_2250, ig_var_2246) + forfrac_var_2245(jl_var_2299, lay_var_2254) * (forref_var_173(indf_var_2250 + 1, ig_var_2246) - forref_var_173(indf_var_2250, ig_var_2246)))
          taug_var_2224(jl_var_2299, 138 + ig_var_2246, lay_var_2254) = tau_major_var_2275(ig_var_2246) + tau_major1_var_2276(ig_var_2246) + tauself_var_2274 + taufor_var_2273
          fracs_var_2240(jl_var_2299, 138 + ig_var_2246, lay_var_2254) = fracrefa_var_168(ig_var_2246, jpl_var_2253) + fpl_var_2285 * (fracrefa_var_168(ig_var_2246, jpl_var_2253 + 1) - fracrefa_var_168(ig_var_2246, jpl_var_2253))
        END DO
      END DO
      ixc0_var_2296 = kfdia_var_2222 - kidia_var_2221 + 1 - ixc0_var_2296
      DO ixp_var_2297 = 1, ixc0_var_2296
        jl_var_2299 = ixhigh_var_2293(ixp_var_2297, lay_var_2254)
        ind0_var_2247 = ((jp_var_2230(jl_var_2299, lay_var_2254) - 13) * 5 + (jt_var_2231(jl_var_2299, lay_var_2254) - 1)) * nspb_var_238(16) + 1
        ind1_var_2248 = ((jp_var_2230(jl_var_2299, lay_var_2254) - 12) * 5 + (jt1_var_2232(jl_var_2299, lay_var_2254) - 1)) * nspb_var_238(16) + 1
        DO ig_var_2246 = 1, 2
          taug_var_2224(jl_var_2299, 138 + ig_var_2246, lay_var_2254) = colch4_var_2235(jl_var_2299, lay_var_2254) * (fac00_var_2226(jl_var_2299, lay_var_2254) * absb_var_171(ind0_var_2247, ig_var_2246) + fac10_var_2228(jl_var_2299, lay_var_2254) * absb_var_171(ind0_var_2247 + 1, ig_var_2246) + fac01_var_2227(jl_var_2299, lay_var_2254) * absb_var_171(ind1_var_2248, ig_var_2246) + fac11_var_2229(jl_var_2299, lay_var_2254) * absb_var_171(ind1_var_2248 + 1, ig_var_2246))
          fracs_var_2240(jl_var_2299, 138 + ig_var_2246, lay_var_2254) = fracrefb_var_169(ig_var_2246)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol16
SUBROUTINE rrtm_taumol12(kidia_var_2300, kfdia_var_2301, klev_var_2302, taug_var_2303, p_tauaerl_var_2304, fac00_var_2305, fac01_var_2306, fac10_var_2307, fac11_var_2308, forfac_var_2324, forfrac_var_2323, indfor_var_2322, jp_var_2309, jt_var_2310, jt1_var_2311, oneminus_var_2312, colh2o_var_2313, colco2_var_2314, laytrop_var_2315, selffac_var_2316, selffrac_var_2317, indself_var_2318, fracs_var_2319, rat_h2oco2_var_2320, rat_h2oco2_1_var_2321)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_237
  USE yoerrtm, ONLY: ng12
  USE yoerrta12, ONLY: absa_var_148, forref_var_150, fracrefa_var_147, selfref_var_149
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2300
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2301
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2302
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2303(kidia_var_2300 : kfdia_var_2301, 140, klev_var_2302)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2304(kidia_var_2300 : kfdia_var_2301, klev_var_2302, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2305(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2306(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2307(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2308(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2309(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2310(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2311(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2312
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2313(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2314(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2315(kidia_var_2300 : kfdia_var_2301)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2316(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2317(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2318(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2319(kidia_var_2300 : kfdia_var_2301, 140, klev_var_2302)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_var_2320(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  REAL(KIND = 8), INTENT(IN) :: rat_h2oco2_1_var_2321(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2322(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2323(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2324(kidia_var_2300 : kfdia_var_2301, klev_var_2302)
  REAL(KIND = 8) :: speccomb_var_2325, speccomb1_var_2326, speccomb_planck_var_2327
  INTEGER(KIND = 4) :: ind0_var_2328, ind1_var_2329, inds_var_2330, indf_var_2331
  INTEGER(KIND = 4) :: ig_var_2332, js_var_2333, lay_var_2334, js1_var_2335, jpl_var_2336
  REAL(KIND = 8) :: fs_var_2337, specmult_var_2338, specparm_var_2339, fs1_var_2340, specmult1_var_2341, specparm1_var_2342, fpl_var_2343, specmult_planck_var_2344, specparm_planck_var_2345
  REAL(KIND = 8) :: fac000_var_2346, fac100_var_2347, fac200_var_2348, fac010_var_2349, fac110_var_2350, fac210_var_2351, fac001_var_2352, fac101_var_2353, fac201_var_2354, fac011_var_2355, fac111_var_2356, fac211_var_2357
  REAL(KIND = 8) :: p_var_2358, p4_var_2359, fk0_var_2360, fk1_var_2361, fk2_var_2362
  REAL(KIND = 8) :: taufor_var_2363, tauself_var_2364, tau_major_var_2365(8), tau_major1_var_2366(8)
  REAL(KIND = 8) :: refrat_planck_a_var_2367
  INTEGER(KIND = 4) :: laytrop_min_var_2368, laytrop_max_var_2369
  INTEGER(KIND = 4) :: ixc_var_2370(klev_var_2302), ixlow_var_2371(kfdia_var_2301, klev_var_2302), ixhigh_var_2372(kfdia_var_2301, klev_var_2302)
  INTEGER(KIND = 4) :: ich_var_2373, icl_var_2374, ixc0_var_2375, ixp_var_2376, jc_var_2377, jl_var_2378
  laytrop_min_var_2368 = MINVAL(laytrop_var_2315)
  laytrop_max_var_2369 = MAXVAL(laytrop_var_2315)
  ixlow_var_2371 = 0
  ixhigh_var_2372 = 0
  ixc_var_2370 = 0
  DO lay_var_2334 = laytrop_min_var_2368 + 1, laytrop_max_var_2369
    icl_var_2374 = 0
    ich_var_2373 = 0
    DO jc_var_2377 = kidia_var_2300, kfdia_var_2301
      IF (lay_var_2334 <= laytrop_var_2315(jc_var_2377)) THEN
        icl_var_2374 = icl_var_2374 + 1
        ixlow_var_2371(icl_var_2374, lay_var_2334) = jc_var_2377
      ELSE
        ich_var_2373 = ich_var_2373 + 1
        ixhigh_var_2372(ich_var_2373, lay_var_2334) = jc_var_2377
      END IF
    END DO
    ixc_var_2370(lay_var_2334) = icl_var_2374
  END DO
  refrat_planck_a_var_2367 = chi_mls(1, 10) / chi_mls(2, 10)
  DO lay_var_2334 = 1, laytrop_min_var_2368
    DO jl_var_2378 = kidia_var_2300, kfdia_var_2301
      speccomb_var_2325 = colh2o_var_2313(jl_var_2378, lay_var_2334) + rat_h2oco2_var_2320(jl_var_2378, lay_var_2334) * colco2_var_2314(jl_var_2378, lay_var_2334)
      specparm_var_2339 = MIN(colh2o_var_2313(jl_var_2378, lay_var_2334) / speccomb_var_2325, oneminus_var_2312)
      specmult_var_2338 = 8.0D0 * (specparm_var_2339)
      js_var_2333 = 1 + INT(specmult_var_2338)
      fs_var_2337 = ((specmult_var_2338) - AINT((specmult_var_2338)))
      speccomb1_var_2326 = colh2o_var_2313(jl_var_2378, lay_var_2334) + rat_h2oco2_1_var_2321(jl_var_2378, lay_var_2334) * colco2_var_2314(jl_var_2378, lay_var_2334)
      specparm1_var_2342 = MIN(colh2o_var_2313(jl_var_2378, lay_var_2334) / speccomb1_var_2326, oneminus_var_2312)
      specmult1_var_2341 = 8.0D0 * (specparm1_var_2342)
      js1_var_2335 = 1 + INT(specmult1_var_2341)
      fs1_var_2340 = ((specmult1_var_2341) - AINT((specmult1_var_2341)))
      speccomb_planck_var_2327 = colh2o_var_2313(jl_var_2378, lay_var_2334) + refrat_planck_a_var_2367 * colco2_var_2314(jl_var_2378, lay_var_2334)
      specparm_planck_var_2345 = MIN(colh2o_var_2313(jl_var_2378, lay_var_2334) / speccomb_planck_var_2327, oneminus_var_2312)
      specmult_planck_var_2344 = 8.0D0 * specparm_planck_var_2345
      jpl_var_2336 = 1 + INT(specmult_planck_var_2344)
      fpl_var_2343 = ((specmult_planck_var_2344) - AINT((specmult_planck_var_2344)))
      ind0_var_2328 = ((jp_var_2309(jl_var_2378, lay_var_2334) - 1) * 5 + (jt_var_2310(jl_var_2378, lay_var_2334) - 1)) * nspa_var_237(12) + js_var_2333
      ind1_var_2329 = (jp_var_2309(jl_var_2378, lay_var_2334) * 5 + (jt1_var_2311(jl_var_2378, lay_var_2334) - 1)) * nspa_var_237(12) + js1_var_2335
      inds_var_2330 = indself_var_2318(jl_var_2378, lay_var_2334)
      indf_var_2331 = indfor_var_2322(jl_var_2378, lay_var_2334)
      IF (specparm_var_2339 .LT. 0.125D0) THEN
        p_var_2358 = fs_var_2337 - 1.0D0
        p4_var_2359 = p_var_2358 ** 4
        fk0_var_2360 = p4_var_2359
        fk1_var_2361 = 1.0D0 - p_var_2358 - 2.0D0 * p4_var_2359
        fk2_var_2362 = p_var_2358 + p4_var_2359
        fac000_var_2346 = fk0_var_2360 * fac00_var_2305(jl_var_2378, lay_var_2334)
        fac100_var_2347 = fk1_var_2361 * fac00_var_2305(jl_var_2378, lay_var_2334)
        fac200_var_2348 = fk2_var_2362 * fac00_var_2305(jl_var_2378, lay_var_2334)
        fac010_var_2349 = fk0_var_2360 * fac10_var_2307(jl_var_2378, lay_var_2334)
        fac110_var_2350 = fk1_var_2361 * fac10_var_2307(jl_var_2378, lay_var_2334)
        fac210_var_2351 = fk2_var_2362 * fac10_var_2307(jl_var_2378, lay_var_2334)
      ELSE IF (specparm_var_2339 .GT. 0.875D0) THEN
        p_var_2358 = - fs_var_2337
        p4_var_2359 = p_var_2358 ** 4
        fk0_var_2360 = p4_var_2359
        fk1_var_2361 = 1.0D0 - p_var_2358 - 2.0D0 * p4_var_2359
        fk2_var_2362 = p_var_2358 + p4_var_2359
        fac000_var_2346 = fk0_var_2360 * fac00_var_2305(jl_var_2378, lay_var_2334)
        fac100_var_2347 = fk1_var_2361 * fac00_var_2305(jl_var_2378, lay_var_2334)
        fac200_var_2348 = fk2_var_2362 * fac00_var_2305(jl_var_2378, lay_var_2334)
        fac010_var_2349 = fk0_var_2360 * fac10_var_2307(jl_var_2378, lay_var_2334)
        fac110_var_2350 = fk1_var_2361 * fac10_var_2307(jl_var_2378, lay_var_2334)
        fac210_var_2351 = fk2_var_2362 * fac10_var_2307(jl_var_2378, lay_var_2334)
      ELSE
        fac000_var_2346 = (1.0D0 - fs_var_2337) * fac00_var_2305(jl_var_2378, lay_var_2334)
        fac010_var_2349 = (1.0D0 - fs_var_2337) * fac10_var_2307(jl_var_2378, lay_var_2334)
        fac100_var_2347 = fs_var_2337 * fac00_var_2305(jl_var_2378, lay_var_2334)
        fac110_var_2350 = fs_var_2337 * fac10_var_2307(jl_var_2378, lay_var_2334)
        fac200_var_2348 = 0.0D0
        fac210_var_2351 = 0.0D0
      END IF
      IF (specparm1_var_2342 .LT. 0.125D0) THEN
        p_var_2358 = fs1_var_2340 - 1.0D0
        p4_var_2359 = p_var_2358 ** 4
        fk0_var_2360 = p4_var_2359
        fk1_var_2361 = 1.0D0 - p_var_2358 - 2.0D0 * p4_var_2359
        fk2_var_2362 = p_var_2358 + p4_var_2359
        fac001_var_2352 = fk0_var_2360 * fac01_var_2306(jl_var_2378, lay_var_2334)
        fac101_var_2353 = fk1_var_2361 * fac01_var_2306(jl_var_2378, lay_var_2334)
        fac201_var_2354 = fk2_var_2362 * fac01_var_2306(jl_var_2378, lay_var_2334)
        fac011_var_2355 = fk0_var_2360 * fac11_var_2308(jl_var_2378, lay_var_2334)
        fac111_var_2356 = fk1_var_2361 * fac11_var_2308(jl_var_2378, lay_var_2334)
        fac211_var_2357 = fk2_var_2362 * fac11_var_2308(jl_var_2378, lay_var_2334)
      ELSE IF (specparm1_var_2342 .GT. 0.875D0) THEN
        p_var_2358 = - fs1_var_2340
        p4_var_2359 = p_var_2358 ** 4
        fk0_var_2360 = p4_var_2359
        fk1_var_2361 = 1.0D0 - p_var_2358 - 2.0D0 * p4_var_2359
        fk2_var_2362 = p_var_2358 + p4_var_2359
        fac001_var_2352 = fk0_var_2360 * fac01_var_2306(jl_var_2378, lay_var_2334)
        fac101_var_2353 = fk1_var_2361 * fac01_var_2306(jl_var_2378, lay_var_2334)
        fac201_var_2354 = fk2_var_2362 * fac01_var_2306(jl_var_2378, lay_var_2334)
        fac011_var_2355 = fk0_var_2360 * fac11_var_2308(jl_var_2378, lay_var_2334)
        fac111_var_2356 = fk1_var_2361 * fac11_var_2308(jl_var_2378, lay_var_2334)
        fac211_var_2357 = fk2_var_2362 * fac11_var_2308(jl_var_2378, lay_var_2334)
      ELSE
        fac001_var_2352 = (1.0D0 - fs1_var_2340) * fac01_var_2306(jl_var_2378, lay_var_2334)
        fac011_var_2355 = (1.0D0 - fs1_var_2340) * fac11_var_2308(jl_var_2378, lay_var_2334)
        fac101_var_2353 = fs1_var_2340 * fac01_var_2306(jl_var_2378, lay_var_2334)
        fac111_var_2356 = fs1_var_2340 * fac11_var_2308(jl_var_2378, lay_var_2334)
        fac201_var_2354 = 0.0D0
        fac211_var_2357 = 0.0D0
      END IF
      IF (specparm_var_2339 .LT. 0.125D0) THEN
        tau_major_var_2365(1 : ng12) = speccomb_var_2325 * (fac000_var_2346 * absa_var_148(ind0_var_2328, 1 : 8) + fac100_var_2347 * absa_var_148(ind0_var_2328 + 1, 1 : 8) + fac200_var_2348 * absa_var_148(ind0_var_2328 + 2, 1 : 8) + fac010_var_2349 * absa_var_148(ind0_var_2328 + 9, 1 : 8) + fac110_var_2350 * absa_var_148(ind0_var_2328 + 10, 1 : 8) + fac210_var_2351 * absa_var_148(ind0_var_2328 + 11, 1 : 8))
      ELSE IF (specparm_var_2339 .GT. 0.875D0) THEN
        tau_major_var_2365(1 : ng12) = speccomb_var_2325 * (fac200_var_2348 * absa_var_148(ind0_var_2328 - 1, 1 : 8) + fac100_var_2347 * absa_var_148(ind0_var_2328, 1 : 8) + fac000_var_2346 * absa_var_148(ind0_var_2328 + 1, 1 : 8) + fac210_var_2351 * absa_var_148(ind0_var_2328 + 8, 1 : 8) + fac110_var_2350 * absa_var_148(ind0_var_2328 + 9, 1 : 8) + fac010_var_2349 * absa_var_148(ind0_var_2328 + 10, 1 : 8))
      ELSE
        tau_major_var_2365(1 : ng12) = speccomb_var_2325 * (fac000_var_2346 * absa_var_148(ind0_var_2328, 1 : 8) + fac100_var_2347 * absa_var_148(ind0_var_2328 + 1, 1 : 8) + fac010_var_2349 * absa_var_148(ind0_var_2328 + 9, 1 : 8) + fac110_var_2350 * absa_var_148(ind0_var_2328 + 10, 1 : 8))
      END IF
      IF (specparm1_var_2342 .LT. 0.125D0) THEN
        tau_major1_var_2366(1 : ng12) = speccomb1_var_2326 * (fac001_var_2352 * absa_var_148(ind1_var_2329, 1 : 8) + fac101_var_2353 * absa_var_148(ind1_var_2329 + 1, 1 : 8) + fac201_var_2354 * absa_var_148(ind1_var_2329 + 2, 1 : 8) + fac011_var_2355 * absa_var_148(ind1_var_2329 + 9, 1 : 8) + fac111_var_2356 * absa_var_148(ind1_var_2329 + 10, 1 : 8) + fac211_var_2357 * absa_var_148(ind1_var_2329 + 11, 1 : 8))
      ELSE IF (specparm1_var_2342 .GT. 0.875D0) THEN
        tau_major1_var_2366(1 : ng12) = speccomb1_var_2326 * (fac201_var_2354 * absa_var_148(ind1_var_2329 - 1, 1 : 8) + fac101_var_2353 * absa_var_148(ind1_var_2329, 1 : 8) + fac001_var_2352 * absa_var_148(ind1_var_2329 + 1, 1 : 8) + fac211_var_2357 * absa_var_148(ind1_var_2329 + 8, 1 : 8) + fac111_var_2356 * absa_var_148(ind1_var_2329 + 9, 1 : 8) + fac011_var_2355 * absa_var_148(ind1_var_2329 + 10, 1 : 8))
      ELSE
        tau_major1_var_2366(1 : ng12) = speccomb1_var_2326 * (fac001_var_2352 * absa_var_148(ind1_var_2329, 1 : 8) + fac101_var_2353 * absa_var_148(ind1_var_2329 + 1, 1 : 8) + fac011_var_2355 * absa_var_148(ind1_var_2329 + 9, 1 : 8) + fac111_var_2356 * absa_var_148(ind1_var_2329 + 10, 1 : 8))
      END IF
      DO ig_var_2332 = 1, 8
        tauself_var_2364 = selffac_var_2316(jl_var_2378, lay_var_2334) * (selfref_var_149(inds_var_2330, ig_var_2332) + selffrac_var_2317(jl_var_2378, lay_var_2334) * (selfref_var_149(inds_var_2330 + 1, ig_var_2332) - selfref_var_149(inds_var_2330, ig_var_2332)))
        taufor_var_2363 = forfac_var_2324(jl_var_2378, lay_var_2334) * (forref_var_150(indf_var_2331, ig_var_2332) + forfrac_var_2323(jl_var_2378, lay_var_2334) * (forref_var_150(indf_var_2331 + 1, ig_var_2332) - forref_var_150(indf_var_2331, ig_var_2332)))
        taug_var_2303(jl_var_2378, 122 + ig_var_2332, lay_var_2334) = tau_major_var_2365(ig_var_2332) + tau_major1_var_2366(ig_var_2332) + tauself_var_2364 + taufor_var_2363
        fracs_var_2319(jl_var_2378, 122 + ig_var_2332, lay_var_2334) = fracrefa_var_147(ig_var_2332, jpl_var_2336) + fpl_var_2343 * (fracrefa_var_147(ig_var_2332, jpl_var_2336 + 1) - fracrefa_var_147(ig_var_2332, jpl_var_2336))
      END DO
    END DO
  END DO
  DO ig_var_2332 = 1, 8
    DO lay_var_2334 = laytrop_max_var_2369 + 1, klev_var_2302
      DO jl_var_2378 = kidia_var_2300, kfdia_var_2301
        taug_var_2303(jl_var_2378, 122 + ig_var_2332, lay_var_2334) = 0.0D0
        fracs_var_2319(jl_var_2378, 122 + ig_var_2332, lay_var_2334) = 0.0D0
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2369 /= laytrop_min_var_2368) THEN
    DO lay_var_2334 = laytrop_min_var_2368 + 1, laytrop_max_var_2369
      ixc0_var_2375 = ixc_var_2370(lay_var_2334)
      DO ixp_var_2376 = 1, ixc0_var_2375
        jl_var_2378 = ixlow_var_2371(ixp_var_2376, lay_var_2334)
        speccomb_var_2325 = colh2o_var_2313(jl_var_2378, lay_var_2334) + rat_h2oco2_var_2320(jl_var_2378, lay_var_2334) * colco2_var_2314(jl_var_2378, lay_var_2334)
        specparm_var_2339 = MIN(colh2o_var_2313(jl_var_2378, lay_var_2334) / speccomb_var_2325, oneminus_var_2312)
        specmult_var_2338 = 8.0D0 * (specparm_var_2339)
        js_var_2333 = 1 + INT(specmult_var_2338)
        fs_var_2337 = ((specmult_var_2338) - AINT((specmult_var_2338)))
        speccomb1_var_2326 = colh2o_var_2313(jl_var_2378, lay_var_2334) + rat_h2oco2_1_var_2321(jl_var_2378, lay_var_2334) * colco2_var_2314(jl_var_2378, lay_var_2334)
        specparm1_var_2342 = MIN(colh2o_var_2313(jl_var_2378, lay_var_2334) / speccomb1_var_2326, oneminus_var_2312)
        specmult1_var_2341 = 8.0D0 * (specparm1_var_2342)
        js1_var_2335 = 1 + INT(specmult1_var_2341)
        fs1_var_2340 = ((specmult1_var_2341) - AINT((specmult1_var_2341)))
        speccomb_planck_var_2327 = colh2o_var_2313(jl_var_2378, lay_var_2334) + refrat_planck_a_var_2367 * colco2_var_2314(jl_var_2378, lay_var_2334)
        specparm_planck_var_2345 = MIN(colh2o_var_2313(jl_var_2378, lay_var_2334) / speccomb_planck_var_2327, oneminus_var_2312)
        specmult_planck_var_2344 = 8.0D0 * specparm_planck_var_2345
        jpl_var_2336 = 1 + INT(specmult_planck_var_2344)
        fpl_var_2343 = ((specmult_planck_var_2344) - AINT((specmult_planck_var_2344)))
        ind0_var_2328 = ((jp_var_2309(jl_var_2378, lay_var_2334) - 1) * 5 + (jt_var_2310(jl_var_2378, lay_var_2334) - 1)) * nspa_var_237(12) + js_var_2333
        ind1_var_2329 = (jp_var_2309(jl_var_2378, lay_var_2334) * 5 + (jt1_var_2311(jl_var_2378, lay_var_2334) - 1)) * nspa_var_237(12) + js1_var_2335
        inds_var_2330 = indself_var_2318(jl_var_2378, lay_var_2334)
        indf_var_2331 = indfor_var_2322(jl_var_2378, lay_var_2334)
        IF (specparm_var_2339 .LT. 0.125D0) THEN
          p_var_2358 = fs_var_2337 - 1.0D0
          p4_var_2359 = p_var_2358 ** 4
          fk0_var_2360 = p4_var_2359
          fk1_var_2361 = 1.0D0 - p_var_2358 - 2.0D0 * p4_var_2359
          fk2_var_2362 = p_var_2358 + p4_var_2359
          fac000_var_2346 = fk0_var_2360 * fac00_var_2305(jl_var_2378, lay_var_2334)
          fac100_var_2347 = fk1_var_2361 * fac00_var_2305(jl_var_2378, lay_var_2334)
          fac200_var_2348 = fk2_var_2362 * fac00_var_2305(jl_var_2378, lay_var_2334)
          fac010_var_2349 = fk0_var_2360 * fac10_var_2307(jl_var_2378, lay_var_2334)
          fac110_var_2350 = fk1_var_2361 * fac10_var_2307(jl_var_2378, lay_var_2334)
          fac210_var_2351 = fk2_var_2362 * fac10_var_2307(jl_var_2378, lay_var_2334)
        ELSE IF (specparm_var_2339 .GT. 0.875D0) THEN
          p_var_2358 = - fs_var_2337
          p4_var_2359 = p_var_2358 ** 4
          fk0_var_2360 = p4_var_2359
          fk1_var_2361 = 1.0D0 - p_var_2358 - 2.0D0 * p4_var_2359
          fk2_var_2362 = p_var_2358 + p4_var_2359
          fac000_var_2346 = fk0_var_2360 * fac00_var_2305(jl_var_2378, lay_var_2334)
          fac100_var_2347 = fk1_var_2361 * fac00_var_2305(jl_var_2378, lay_var_2334)
          fac200_var_2348 = fk2_var_2362 * fac00_var_2305(jl_var_2378, lay_var_2334)
          fac010_var_2349 = fk0_var_2360 * fac10_var_2307(jl_var_2378, lay_var_2334)
          fac110_var_2350 = fk1_var_2361 * fac10_var_2307(jl_var_2378, lay_var_2334)
          fac210_var_2351 = fk2_var_2362 * fac10_var_2307(jl_var_2378, lay_var_2334)
        ELSE
          fac000_var_2346 = (1.0D0 - fs_var_2337) * fac00_var_2305(jl_var_2378, lay_var_2334)
          fac010_var_2349 = (1.0D0 - fs_var_2337) * fac10_var_2307(jl_var_2378, lay_var_2334)
          fac100_var_2347 = fs_var_2337 * fac00_var_2305(jl_var_2378, lay_var_2334)
          fac110_var_2350 = fs_var_2337 * fac10_var_2307(jl_var_2378, lay_var_2334)
          fac200_var_2348 = 0.0D0
          fac210_var_2351 = 0.0D0
        END IF
        IF (specparm1_var_2342 .LT. 0.125D0) THEN
          p_var_2358 = fs1_var_2340 - 1.0D0
          p4_var_2359 = p_var_2358 ** 4
          fk0_var_2360 = p4_var_2359
          fk1_var_2361 = 1.0D0 - p_var_2358 - 2.0D0 * p4_var_2359
          fk2_var_2362 = p_var_2358 + p4_var_2359
          fac001_var_2352 = fk0_var_2360 * fac01_var_2306(jl_var_2378, lay_var_2334)
          fac101_var_2353 = fk1_var_2361 * fac01_var_2306(jl_var_2378, lay_var_2334)
          fac201_var_2354 = fk2_var_2362 * fac01_var_2306(jl_var_2378, lay_var_2334)
          fac011_var_2355 = fk0_var_2360 * fac11_var_2308(jl_var_2378, lay_var_2334)
          fac111_var_2356 = fk1_var_2361 * fac11_var_2308(jl_var_2378, lay_var_2334)
          fac211_var_2357 = fk2_var_2362 * fac11_var_2308(jl_var_2378, lay_var_2334)
        ELSE IF (specparm1_var_2342 .GT. 0.875D0) THEN
          p_var_2358 = - fs1_var_2340
          p4_var_2359 = p_var_2358 ** 4
          fk0_var_2360 = p4_var_2359
          fk1_var_2361 = 1.0D0 - p_var_2358 - 2.0D0 * p4_var_2359
          fk2_var_2362 = p_var_2358 + p4_var_2359
          fac001_var_2352 = fk0_var_2360 * fac01_var_2306(jl_var_2378, lay_var_2334)
          fac101_var_2353 = fk1_var_2361 * fac01_var_2306(jl_var_2378, lay_var_2334)
          fac201_var_2354 = fk2_var_2362 * fac01_var_2306(jl_var_2378, lay_var_2334)
          fac011_var_2355 = fk0_var_2360 * fac11_var_2308(jl_var_2378, lay_var_2334)
          fac111_var_2356 = fk1_var_2361 * fac11_var_2308(jl_var_2378, lay_var_2334)
          fac211_var_2357 = fk2_var_2362 * fac11_var_2308(jl_var_2378, lay_var_2334)
        ELSE
          fac001_var_2352 = (1.0D0 - fs1_var_2340) * fac01_var_2306(jl_var_2378, lay_var_2334)
          fac011_var_2355 = (1.0D0 - fs1_var_2340) * fac11_var_2308(jl_var_2378, lay_var_2334)
          fac101_var_2353 = fs1_var_2340 * fac01_var_2306(jl_var_2378, lay_var_2334)
          fac111_var_2356 = fs1_var_2340 * fac11_var_2308(jl_var_2378, lay_var_2334)
          fac201_var_2354 = 0.0D0
          fac211_var_2357 = 0.0D0
        END IF
        IF (specparm_var_2339 .LT. 0.125D0) THEN
          tau_major_var_2365(1 : ng12) = speccomb_var_2325 * (fac000_var_2346 * absa_var_148(ind0_var_2328, 1 : 8) + fac100_var_2347 * absa_var_148(ind0_var_2328 + 1, 1 : 8) + fac200_var_2348 * absa_var_148(ind0_var_2328 + 2, 1 : 8) + fac010_var_2349 * absa_var_148(ind0_var_2328 + 9, 1 : 8) + fac110_var_2350 * absa_var_148(ind0_var_2328 + 10, 1 : 8) + fac210_var_2351 * absa_var_148(ind0_var_2328 + 11, 1 : 8))
        ELSE IF (specparm_var_2339 .GT. 0.875D0) THEN
          tau_major_var_2365(1 : ng12) = speccomb_var_2325 * (fac200_var_2348 * absa_var_148(ind0_var_2328 - 1, 1 : 8) + fac100_var_2347 * absa_var_148(ind0_var_2328, 1 : 8) + fac000_var_2346 * absa_var_148(ind0_var_2328 + 1, 1 : 8) + fac210_var_2351 * absa_var_148(ind0_var_2328 + 8, 1 : 8) + fac110_var_2350 * absa_var_148(ind0_var_2328 + 9, 1 : 8) + fac010_var_2349 * absa_var_148(ind0_var_2328 + 10, 1 : 8))
        ELSE
          tau_major_var_2365(1 : ng12) = speccomb_var_2325 * (fac000_var_2346 * absa_var_148(ind0_var_2328, 1 : 8) + fac100_var_2347 * absa_var_148(ind0_var_2328 + 1, 1 : 8) + fac010_var_2349 * absa_var_148(ind0_var_2328 + 9, 1 : 8) + fac110_var_2350 * absa_var_148(ind0_var_2328 + 10, 1 : 8))
        END IF
        IF (specparm1_var_2342 .LT. 0.125D0) THEN
          tau_major1_var_2366(1 : ng12) = speccomb1_var_2326 * (fac001_var_2352 * absa_var_148(ind1_var_2329, 1 : 8) + fac101_var_2353 * absa_var_148(ind1_var_2329 + 1, 1 : 8) + fac201_var_2354 * absa_var_148(ind1_var_2329 + 2, 1 : 8) + fac011_var_2355 * absa_var_148(ind1_var_2329 + 9, 1 : 8) + fac111_var_2356 * absa_var_148(ind1_var_2329 + 10, 1 : 8) + fac211_var_2357 * absa_var_148(ind1_var_2329 + 11, 1 : 8))
        ELSE IF (specparm1_var_2342 .GT. 0.875D0) THEN
          tau_major1_var_2366(1 : ng12) = speccomb1_var_2326 * (fac201_var_2354 * absa_var_148(ind1_var_2329 - 1, 1 : 8) + fac101_var_2353 * absa_var_148(ind1_var_2329, 1 : 8) + fac001_var_2352 * absa_var_148(ind1_var_2329 + 1, 1 : 8) + fac211_var_2357 * absa_var_148(ind1_var_2329 + 8, 1 : 8) + fac111_var_2356 * absa_var_148(ind1_var_2329 + 9, 1 : 8) + fac011_var_2355 * absa_var_148(ind1_var_2329 + 10, 1 : 8))
        ELSE
          tau_major1_var_2366(1 : ng12) = speccomb1_var_2326 * (fac001_var_2352 * absa_var_148(ind1_var_2329, 1 : 8) + fac101_var_2353 * absa_var_148(ind1_var_2329 + 1, 1 : 8) + fac011_var_2355 * absa_var_148(ind1_var_2329 + 9, 1 : 8) + fac111_var_2356 * absa_var_148(ind1_var_2329 + 10, 1 : 8))
        END IF
        DO ig_var_2332 = 1, 8
          tauself_var_2364 = selffac_var_2316(jl_var_2378, lay_var_2334) * (selfref_var_149(inds_var_2330, ig_var_2332) + selffrac_var_2317(jl_var_2378, lay_var_2334) * (selfref_var_149(inds_var_2330 + 1, ig_var_2332) - selfref_var_149(inds_var_2330, ig_var_2332)))
          taufor_var_2363 = forfac_var_2324(jl_var_2378, lay_var_2334) * (forref_var_150(indf_var_2331, ig_var_2332) + forfrac_var_2323(jl_var_2378, lay_var_2334) * (forref_var_150(indf_var_2331 + 1, ig_var_2332) - forref_var_150(indf_var_2331, ig_var_2332)))
          taug_var_2303(jl_var_2378, 122 + ig_var_2332, lay_var_2334) = tau_major_var_2365(ig_var_2332) + tau_major1_var_2366(ig_var_2332) + tauself_var_2364 + taufor_var_2363
          fracs_var_2319(jl_var_2378, 122 + ig_var_2332, lay_var_2334) = fracrefa_var_147(ig_var_2332, jpl_var_2336) + fpl_var_2343 * (fracrefa_var_147(ig_var_2332, jpl_var_2336 + 1) - fracrefa_var_147(ig_var_2332, jpl_var_2336))
        END DO
      END DO
      ixc0_var_2375 = kfdia_var_2301 - kidia_var_2300 + 1 - ixc0_var_2375
      DO ig_var_2332 = 1, 8
        DO ixp_var_2376 = 1, ixc0_var_2375
          jl_var_2378 = ixhigh_var_2372(ixp_var_2376, lay_var_2334)
          taug_var_2303(jl_var_2378, 122 + ig_var_2332, lay_var_2334) = 0.0D0
          fracs_var_2319(jl_var_2378, 122 + ig_var_2332, lay_var_2334) = 0.0D0
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE rrtm_taumol12
SUBROUTINE rrtm_taumol13(kidia_var_2379, kfdia_var_2380, klev_var_2381, taug_var_2382, p_tauaerl_var_2383, fac00_var_2384, fac01_var_2385, fac10_var_2386, fac11_var_2387, forfac_var_2403, forfrac_var_2404, indfor_var_2402, jp_var_2388, jt_var_2389, jt1_var_2390, oneminus_var_2391, colh2o_var_2392, coln2o_var_2393, colco2_var_2394, colo3_var_2395, coldry_var_2396, laytrop_var_2397, selffac_var_2398, selffrac_var_2399, indself_var_2400, fracs_var_2401, rat_h2on2o, rat_h2on2o_1, minorfrac_var_2405, indminor_var_2406)
  USE yoerrtrf, ONLY: chi_mls
  USE yoerrtwn, ONLY: nspa_var_237
  USE yoerrtm, ONLY: ng13
  USE yoerrta13, ONLY: absa_var_153, forref_var_155, fracrefa_var_151, fracrefb_var_152, ka_mco, ka_mco2_var_156, kb_mo3, selfref_var_154
  IMPLICIT NONE
  INTEGER(KIND = 4), INTENT(IN) :: kidia_var_2379
  INTEGER(KIND = 4), INTENT(IN) :: kfdia_var_2380
  INTEGER(KIND = 4), INTENT(IN) :: klev_var_2381
  REAL(KIND = 8), INTENT(INOUT) :: taug_var_2382(kidia_var_2379 : kfdia_var_2380, 140, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: p_tauaerl_var_2383(kidia_var_2379 : kfdia_var_2380, klev_var_2381, 16)
  REAL(KIND = 8), INTENT(IN) :: fac00_var_2384(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: fac01_var_2385(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: fac10_var_2386(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: fac11_var_2387(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  INTEGER(KIND = 4), INTENT(IN) :: jp_var_2388(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  INTEGER(KIND = 4), INTENT(IN) :: jt_var_2389(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  INTEGER(KIND = 4), INTENT(IN) :: jt1_var_2390(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: oneminus_var_2391
  REAL(KIND = 8), INTENT(IN) :: colh2o_var_2392(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: coln2o_var_2393(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: colco2_var_2394(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: colo3_var_2395(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: coldry_var_2396(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  INTEGER(KIND = 4), INTENT(IN) :: laytrop_var_2397(kidia_var_2379 : kfdia_var_2380)
  REAL(KIND = 8), INTENT(IN) :: selffac_var_2398(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: selffrac_var_2399(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  INTEGER(KIND = 4), INTENT(IN) :: indself_var_2400(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(INOUT) :: fracs_var_2401(kidia_var_2379 : kfdia_var_2380, 140, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: rat_h2on2o(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: rat_h2on2o_1(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  INTEGER(KIND = 4), INTENT(IN) :: indfor_var_2402(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: forfac_var_2403(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: forfrac_var_2404(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8), INTENT(IN) :: minorfrac_var_2405(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  INTEGER(KIND = 4), INTENT(IN) :: indminor_var_2406(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  REAL(KIND = 8) :: speccomb_var_2407, speccomb1_var_2408, speccomb_planck_var_2409, speccomb_mco2_var_2410, speccomb_mco
  INTEGER(KIND = 4) :: ind0_var_2411, ind1_var_2412, inds_var_2413, indf_var_2414, indm_var_2415
  INTEGER(KIND = 4) :: ig_var_2416, js_var_2417, lay_var_2418, js1_var_2419, jpl_var_2420, jmco2_var_2421, jmco
  REAL(KIND = 8) :: refrat_planck_a_var_2422, refrat_m_a_var_2423, refrat_m_a3
  REAL(KIND = 8) :: fac000_var_2424, fac100_var_2425, fac200_var_2426, fac010_var_2427, fac110_var_2428, fac210_var_2429, fac001_var_2430, fac101_var_2431, fac201_var_2432, fac011_var_2433, fac111_var_2434, fac211_var_2435
  REAL(KIND = 8) :: p_var_2436, p4_var_2437, fk0_var_2438, fk1_var_2439, fk2_var_2440
  REAL(KIND = 8) :: taufor_var_2441, tauself_var_2442, tau_major_var_2443(4), tau_major1_var_2444(4), co2m1_var_2445, co2m2_var_2446, absco2_var_2447
  REAL(KIND = 8) :: com1, com2, absco, abso3_var_2448
  REAL(KIND = 8) :: chi_co2_var_2449, ratco2_var_2450, adjfac_var_2451, adjcolco2_var_2452
  REAL(KIND = 8) :: fs_var_2453, specmult_var_2454, specparm_var_2455, fs1_var_2456, specmult1_var_2457, specparm1_var_2458, fmco2_var_2459, specmult_mco2_var_2460, specparm_mco2_var_2461, fmco, specmult_mco, specparm_mco, fpl_var_2462, specmult_planck_var_2463, specparm_planck_var_2464
  REAL(KIND = 8) :: colco(kidia_var_2379 : kfdia_var_2380, klev_var_2381)
  INTEGER(KIND = 4) :: laytrop_min_var_2465, laytrop_max_var_2466
  INTEGER(KIND = 4) :: ixc_var_2467(klev_var_2381), ixlow_var_2468(kfdia_var_2380, klev_var_2381), ixhigh_var_2469(kfdia_var_2380, klev_var_2381)
  INTEGER(KIND = 4) :: ich_var_2470, icl_var_2471, ixc0_var_2472, ixp_var_2473, jc_var_2474, jl_var_2475
  DO lay_var_2418 = 1, klev_var_2381
    DO jc_var_2474 = kidia_var_2379, kfdia_var_2380
      colco(jc_var_2474, lay_var_2418) = 0.0D0
    END DO
  END DO
  laytrop_min_var_2465 = MINVAL(laytrop_var_2397)
  laytrop_max_var_2466 = MAXVAL(laytrop_var_2397)
  ixlow_var_2468 = 0
  ixhigh_var_2469 = 0
  ixc_var_2467 = 0
  DO lay_var_2418 = laytrop_min_var_2465 + 1, laytrop_max_var_2466
    icl_var_2471 = 0
    ich_var_2470 = 0
    DO jc_var_2474 = kidia_var_2379, kfdia_var_2380
      IF (lay_var_2418 <= laytrop_var_2397(jc_var_2474)) THEN
        icl_var_2471 = icl_var_2471 + 1
        ixlow_var_2468(icl_var_2471, lay_var_2418) = jc_var_2474
      ELSE
        ich_var_2470 = ich_var_2470 + 1
        ixhigh_var_2469(ich_var_2470, lay_var_2418) = jc_var_2474
      END IF
    END DO
    ixc_var_2467(lay_var_2418) = icl_var_2471
  END DO
  refrat_planck_a_var_2422 = chi_mls(1, 5) / chi_mls(4, 5)
  refrat_m_a_var_2423 = chi_mls(1, 1) / chi_mls(4, 1)
  refrat_m_a3 = chi_mls(1, 3) / chi_mls(4, 3)
  DO lay_var_2418 = 1, laytrop_min_var_2465
    DO jl_var_2475 = kidia_var_2379, kfdia_var_2380
      speccomb_var_2407 = colh2o_var_2392(jl_var_2475, lay_var_2418) + rat_h2on2o(jl_var_2475, lay_var_2418) * coln2o_var_2393(jl_var_2475, lay_var_2418)
      specparm_var_2455 = MIN(colh2o_var_2392(jl_var_2475, lay_var_2418) / speccomb_var_2407, oneminus_var_2391)
      specmult_var_2454 = 8.0D0 * (specparm_var_2455)
      js_var_2417 = 1 + INT(specmult_var_2454)
      fs_var_2453 = ((specmult_var_2454) - AINT((specmult_var_2454)))
      speccomb1_var_2408 = colh2o_var_2392(jl_var_2475, lay_var_2418) + rat_h2on2o_1(jl_var_2475, lay_var_2418) * coln2o_var_2393(jl_var_2475, lay_var_2418)
      specparm1_var_2458 = MIN(colh2o_var_2392(jl_var_2475, lay_var_2418) / speccomb1_var_2408, oneminus_var_2391)
      specmult1_var_2457 = 8.0D0 * (specparm1_var_2458)
      js1_var_2419 = 1 + INT(specmult1_var_2457)
      fs1_var_2456 = ((specmult1_var_2457) - AINT((specmult1_var_2457)))
      speccomb_mco2_var_2410 = colh2o_var_2392(jl_var_2475, lay_var_2418) + refrat_m_a_var_2423 * coln2o_var_2393(jl_var_2475, lay_var_2418)
      specparm_mco2_var_2461 = MIN(colh2o_var_2392(jl_var_2475, lay_var_2418) / speccomb_mco2_var_2410, oneminus_var_2391)
      specmult_mco2_var_2460 = 8.0D0 * specparm_mco2_var_2461
      jmco2_var_2421 = 1 + INT(specmult_mco2_var_2460)
      fmco2_var_2459 = ((specmult_mco2_var_2460) - AINT((specmult_mco2_var_2460)))
      chi_co2_var_2449 = colco2_var_2394(jl_var_2475, lay_var_2418) / (coldry_var_2396(jl_var_2475, lay_var_2418))
      ratco2_var_2450 = 1D+20 * chi_co2_var_2449 / 0.000355D0
      IF (ratco2_var_2450 .GT. 3.0D0) THEN
        adjfac_var_2451 = 2.0D0 + (ratco2_var_2450 - 2.0D0) ** 0.68D0
        adjcolco2_var_2452 = adjfac_var_2451 * 0.000355D0 * coldry_var_2396(jl_var_2475, lay_var_2418) * 1D-20
      ELSE
        adjcolco2_var_2452 = colco2_var_2394(jl_var_2475, lay_var_2418)
      END IF
      speccomb_mco = colh2o_var_2392(jl_var_2475, lay_var_2418) + refrat_m_a3 * coln2o_var_2393(jl_var_2475, lay_var_2418)
      specparm_mco = MIN(colh2o_var_2392(jl_var_2475, lay_var_2418) / speccomb_mco, oneminus_var_2391)
      specmult_mco = 8.0D0 * specparm_mco
      jmco = 1 + INT(specmult_mco)
      fmco = ((specmult_mco) - AINT((specmult_mco)))
      speccomb_planck_var_2409 = colh2o_var_2392(jl_var_2475, lay_var_2418) + refrat_planck_a_var_2422 * coln2o_var_2393(jl_var_2475, lay_var_2418)
      specparm_planck_var_2464 = MIN(colh2o_var_2392(jl_var_2475, lay_var_2418) / speccomb_planck_var_2409, oneminus_var_2391)
      specmult_planck_var_2463 = 8.0D0 * specparm_planck_var_2464
      jpl_var_2420 = 1 + INT(specmult_planck_var_2463)
      fpl_var_2462 = ((specmult_planck_var_2463) - AINT((specmult_planck_var_2463)))
      ind0_var_2411 = ((jp_var_2388(jl_var_2475, lay_var_2418) - 1) * 5 + (jt_var_2389(jl_var_2475, lay_var_2418) - 1)) * nspa_var_237(13) + js_var_2417
      ind1_var_2412 = (jp_var_2388(jl_var_2475, lay_var_2418) * 5 + (jt1_var_2390(jl_var_2475, lay_var_2418) - 1)) * nspa_var_237(13) + js1_var_2419
      inds_var_2413 = indself_var_2400(jl_var_2475, lay_var_2418)
      indf_var_2414 = indfor_var_2402(jl_var_2475, lay_var_2418)
      indm_var_2415 = indminor_var_2406(jl_var_2475, lay_var_2418)
      IF (specparm_var_2455 .LT. 0.125D0) THEN
        p_var_2436 = fs_var_2453 - 1.0D0
        p4_var_2437 = p_var_2436 ** 4
        fk0_var_2438 = p4_var_2437
        fk1_var_2439 = 1.0D0 - p_var_2436 - 2.0D0 * p4_var_2437
        fk2_var_2440 = p_var_2436 + p4_var_2437
        fac000_var_2424 = fk0_var_2438 * fac00_var_2384(jl_var_2475, lay_var_2418)
        fac100_var_2425 = fk1_var_2439 * fac00_var_2384(jl_var_2475, lay_var_2418)
        fac200_var_2426 = fk2_var_2440 * fac00_var_2384(jl_var_2475, lay_var_2418)
        fac010_var_2427 = fk0_var_2438 * fac10_var_2386(jl_var_2475, lay_var_2418)
        fac110_var_2428 = fk1_var_2439 * fac10_var_2386(jl_var_2475, lay_var_2418)
        fac210_var_2429 = fk2_var_2440 * fac10_var_2386(jl_var_2475, lay_var_2418)
      ELSE IF (specparm_var_2455 .GT. 0.875D0) THEN
        p_var_2436 = - fs_var_2453
        p4_var_2437 = p_var_2436 ** 4
        fk0_var_2438 = p4_var_2437
        fk1_var_2439 = 1.0D0 - p_var_2436 - 2.0D0 * p4_var_2437
        fk2_var_2440 = p_var_2436 + p4_var_2437
        fac000_var_2424 = fk0_var_2438 * fac00_var_2384(jl_var_2475, lay_var_2418)
        fac100_var_2425 = fk1_var_2439 * fac00_var_2384(jl_var_2475, lay_var_2418)
        fac200_var_2426 = fk2_var_2440 * fac00_var_2384(jl_var_2475, lay_var_2418)
        fac010_var_2427 = fk0_var_2438 * fac10_var_2386(jl_var_2475, lay_var_2418)
        fac110_var_2428 = fk1_var_2439 * fac10_var_2386(jl_var_2475, lay_var_2418)
        fac210_var_2429 = fk2_var_2440 * fac10_var_2386(jl_var_2475, lay_var_2418)
      ELSE
        fac000_var_2424 = (1.0D0 - fs_var_2453) * fac00_var_2384(jl_var_2475, lay_var_2418)
        fac010_var_2427 = (1.0D0 - fs_var_2453) * fac10_var_2386(jl_var_2475, lay_var_2418)
        fac100_var_2425 = fs_var_2453 * fac00_var_2384(jl_var_2475, lay_var_2418)
        fac110_var_2428 = fs_var_2453 * fac10_var_2386(jl_var_2475, lay_var_2418)
        fac200_var_2426 = 0.0D0
        fac210_var_2429 = 0.0D0
      END IF
      IF (specparm1_var_2458 .LT. 0.125D0) THEN
        p_var_2436 = fs1_var_2456 - 1.0D0
        p4_var_2437 = p_var_2436 ** 4
        fk0_var_2438 = p4_var_2437
        fk1_var_2439 = 1.0D0 - p_var_2436 - 2.0D0 * p4_var_2437
        fk2_var_2440 = p_var_2436 + p4_var_2437
        fac001_var_2430 = fk0_var_2438 * fac01_var_2385(jl_var_2475, lay_var_2418)
        fac101_var_2431 = fk1_var_2439 * fac01_var_2385(jl_var_2475, lay_var_2418)
        fac201_var_2432 = fk2_var_2440 * fac01_var_2385(jl_var_2475, lay_var_2418)
        fac011_var_2433 = fk0_var_2438 * fac11_var_2387(jl_var_2475, lay_var_2418)
        fac111_var_2434 = fk1_var_2439 * fac11_var_2387(jl_var_2475, lay_var_2418)
        fac211_var_2435 = fk2_var_2440 * fac11_var_2387(jl_var_2475, lay_var_2418)
      ELSE IF (specparm1_var_2458 .GT. 0.875D0) THEN
        p_var_2436 = - fs1_var_2456
        p4_var_2437 = p_var_2436 ** 4
        fk0_var_2438 = p4_var_2437
        fk1_var_2439 = 1.0D0 - p_var_2436 - 2.0D0 * p4_var_2437
        fk2_var_2440 = p_var_2436 + p4_var_2437
        fac001_var_2430 = fk0_var_2438 * fac01_var_2385(jl_var_2475, lay_var_2418)
        fac101_var_2431 = fk1_var_2439 * fac01_var_2385(jl_var_2475, lay_var_2418)
        fac201_var_2432 = fk2_var_2440 * fac01_var_2385(jl_var_2475, lay_var_2418)
        fac011_var_2433 = fk0_var_2438 * fac11_var_2387(jl_var_2475, lay_var_2418)
        fac111_var_2434 = fk1_var_2439 * fac11_var_2387(jl_var_2475, lay_var_2418)
        fac211_var_2435 = fk2_var_2440 * fac11_var_2387(jl_var_2475, lay_var_2418)
      ELSE
        fac001_var_2430 = (1.0D0 - fs1_var_2456) * fac01_var_2385(jl_var_2475, lay_var_2418)
        fac011_var_2433 = (1.0D0 - fs1_var_2456) * fac11_var_2387(jl_var_2475, lay_var_2418)
        fac101_var_2431 = fs1_var_2456 * fac01_var_2385(jl_var_2475, lay_var_2418)
        fac111_var_2434 = fs1_var_2456 * fac11_var_2387(jl_var_2475, lay_var_2418)
        fac201_var_2432 = 0.0D0
        fac211_var_2435 = 0.0D0
      END IF
      IF (specparm_var_2455 .LT. 0.125D0) THEN
        tau_major_var_2443(1 : ng13) = speccomb_var_2407 * (fac000_var_2424 * absa_var_153(ind0_var_2411, 1 : 4) + fac100_var_2425 * absa_var_153(ind0_var_2411 + 1, 1 : 4) + fac200_var_2426 * absa_var_153(ind0_var_2411 + 2, 1 : 4) + fac010_var_2427 * absa_var_153(ind0_var_2411 + 9, 1 : 4) + fac110_var_2428 * absa_var_153(ind0_var_2411 + 10, 1 : 4) + fac210_var_2429 * absa_var_153(ind0_var_2411 + 11, 1 : 4))
      ELSE IF (specparm_var_2455 .GT. 0.875D0) THEN
        tau_major_var_2443(1 : ng13) = speccomb_var_2407 * (fac200_var_2426 * absa_var_153(ind0_var_2411 - 1, 1 : 4) + fac100_var_2425 * absa_var_153(ind0_var_2411, 1 : 4) + fac000_var_2424 * absa_var_153(ind0_var_2411 + 1, 1 : 4) + fac210_var_2429 * absa_var_153(ind0_var_2411 + 8, 1 : 4) + fac110_var_2428 * absa_var_153(ind0_var_2411 + 9, 1 : 4) + fac010_var_2427 * absa_var_153(ind0_var_2411 + 10, 1 : 4))
      ELSE
        tau_major_var_2443(1 : ng13) = speccomb_var_2407 * (fac000_var_2424 * absa_var_153(ind0_var_2411, 1 : 4) + fac100_var_2425 * absa_var_153(ind0_var_2411 + 1, 1 : 4) + fac010_var_2427 * absa_var_153(ind0_var_2411 + 9, 1 : 4) + fac110_var_2428 * absa_var_153(ind0_var_2411 + 10, 1 : 4))
      END IF
      IF (specparm1_var_2458 .LT. 0.125D0) THEN
        tau_major1_var_2444(1 : ng13) = speccomb1_var_2408 * (fac001_var_2430 * absa_var_153(ind1_var_2412, 1 : 4) + fac101_var_2431 * absa_var_153(ind1_var_2412 + 1, 1 : 4) + fac201_var_2432 * absa_var_153(ind1_var_2412 + 2, 1 : 4) + fac011_var_2433 * absa_var_153(ind1_var_2412 + 9, 1 : 4) + fac111_var_2434 * absa_var_153(ind1_var_2412 + 10, 1 : 4) + fac211_var_2435 * absa_var_153(ind1_var_2412 + 11, 1 : 4))
      ELSE IF (specparm1_var_2458 .GT. 0.875D0) THEN
        tau_major1_var_2444(1 : ng13) = speccomb1_var_2408 * (fac201_var_2432 * absa_var_153(ind1_var_2412 - 1, 1 : 4) + fac101_var_2431 * absa_var_153(ind1_var_2412, 1 : 4) + fac001_var_2430 * absa_var_153(ind1_var_2412 + 1, 1 : 4) + fac211_var_2435 * absa_var_153(ind1_var_2412 + 8, 1 : 4) + fac111_var_2434 * absa_var_153(ind1_var_2412 + 9, 1 : 4) + fac011_var_2433 * absa_var_153(ind1_var_2412 + 10, 1 : 4))
      ELSE
        tau_major1_var_2444(1 : ng13) = speccomb1_var_2408 * (fac001_var_2430 * absa_var_153(ind1_var_2412, 1 : 4) + fac101_var_2431 * absa_var_153(ind1_var_2412 + 1, 1 : 4) + fac011_var_2433 * absa_var_153(ind1_var_2412 + 9, 1 : 4) + fac111_var_2434 * absa_var_153(ind1_var_2412 + 10, 1 : 4))
      END IF
      DO ig_var_2416 = 1, 4
        tauself_var_2442 = selffac_var_2398(jl_var_2475, lay_var_2418) * (selfref_var_154(inds_var_2413, ig_var_2416) + selffrac_var_2399(jl_var_2475, lay_var_2418) * (selfref_var_154(inds_var_2413 + 1, ig_var_2416) - selfref_var_154(inds_var_2413, ig_var_2416)))
        taufor_var_2441 = forfac_var_2403(jl_var_2475, lay_var_2418) * (forref_var_155(indf_var_2414, ig_var_2416) + forfrac_var_2404(jl_var_2475, lay_var_2418) * (forref_var_155(indf_var_2414 + 1, ig_var_2416) - forref_var_155(indf_var_2414, ig_var_2416)))
        co2m1_var_2445 = ka_mco2_var_156(jmco2_var_2421, indm_var_2415, ig_var_2416) + fmco2_var_2459 * (ka_mco2_var_156(jmco2_var_2421 + 1, indm_var_2415, ig_var_2416) - ka_mco2_var_156(jmco2_var_2421, indm_var_2415, ig_var_2416))
        co2m2_var_2446 = ka_mco2_var_156(jmco2_var_2421, indm_var_2415 + 1, ig_var_2416) + fmco2_var_2459 * (ka_mco2_var_156(jmco2_var_2421 + 1, indm_var_2415 + 1, ig_var_2416) - ka_mco2_var_156(jmco2_var_2421, indm_var_2415 + 1, ig_var_2416))
        absco2_var_2447 = co2m1_var_2445 + minorfrac_var_2405(jl_var_2475, lay_var_2418) * (co2m2_var_2446 - co2m1_var_2445)
        com1 = ka_mco(jmco, indm_var_2415, ig_var_2416) + fmco * (ka_mco(jmco + 1, indm_var_2415, ig_var_2416) - ka_mco(jmco, indm_var_2415, ig_var_2416))
        com2 = ka_mco(jmco, indm_var_2415 + 1, ig_var_2416) + fmco * (ka_mco(jmco + 1, indm_var_2415 + 1, ig_var_2416) - ka_mco(jmco, indm_var_2415 + 1, ig_var_2416))
        absco = com1 + minorfrac_var_2405(jl_var_2475, lay_var_2418) * (com2 - com1)
        taug_var_2382(jl_var_2475, 130 + ig_var_2416, lay_var_2418) = tau_major_var_2443(ig_var_2416) + tau_major1_var_2444(ig_var_2416) + tauself_var_2442 + taufor_var_2441 + adjcolco2_var_2452 * absco2_var_2447 + colco(jl_var_2475, lay_var_2418) * absco
        fracs_var_2401(jl_var_2475, 130 + ig_var_2416, lay_var_2418) = fracrefa_var_151(ig_var_2416, jpl_var_2420) + fpl_var_2462 * (fracrefa_var_151(ig_var_2416, jpl_var_2420 + 1) - fracrefa_var_151(ig_var_2416, jpl_var_2420))
      END DO
    END DO
  END DO
  DO lay_var_2418 = laytrop_max_var_2466 + 1, klev_var_2381
    DO jl_var_2475 = kidia_var_2379, kfdia_var_2380
      indm_var_2415 = indminor_var_2406(jl_var_2475, lay_var_2418)
      DO ig_var_2416 = 1, 4
        abso3_var_2448 = kb_mo3(indm_var_2415, ig_var_2416) + minorfrac_var_2405(jl_var_2475, lay_var_2418) * (kb_mo3(indm_var_2415 + 1, ig_var_2416) - kb_mo3(indm_var_2415, ig_var_2416))
        taug_var_2382(jl_var_2475, 130 + ig_var_2416, lay_var_2418) = colo3_var_2395(jl_var_2475, lay_var_2418) * abso3_var_2448
        fracs_var_2401(jl_var_2475, 130 + ig_var_2416, lay_var_2418) = fracrefb_var_152(ig_var_2416)
      END DO
    END DO
  END DO
  IF (laytrop_max_var_2466 /= laytrop_min_var_2465) THEN
    DO lay_var_2418 = laytrop_min_var_2465 + 1, laytrop_max_var_2466
      ixc0_var_2472 = ixc_var_2467(lay_var_2418)
      DO ixp_var_2473 = 1, ixc0_var_2472
        jl_var_2475 = ixlow_var_2468(ixp_var_2473, lay_var_2418)
        speccomb_var_2407 = colh2o_var_2392(jl_var_2475, lay_var_2418) + rat_h2on2o(jl_var_2475, lay_var_2418) * coln2o_var_2393(jl_var_2475, lay_var_2418)
        specparm_var_2455 = MIN(colh2o_var_2392(jl_var_2475, lay_var_2418) / speccomb_var_2407, oneminus_var_2391)
        specmult_var_2454 = 8.0D0 * (specparm_var_2455)
        js_var_2417 = 1 + INT(specmult_var_2454)
        fs_var_2453 = ((specmult_var_2454) - AINT((specmult_var_2454)))
        speccomb1_var_2408 = colh2o_var_2392(jl_var_2475, lay_var_2418) + rat_h2on2o_1(jl_var_2475, lay_var_2418) * coln2o_var_2393(jl_var_2475, lay_var_2418)
        specparm1_var_2458 = MIN(colh2o_var_2392(jl_var_2475, lay_var_2418) / speccomb1_var_2408, oneminus_var_2391)
        specmult1_var_2457 = 8.0D0 * (specparm1_var_2458)
        js1_var_2419 = 1 + INT(specmult1_var_2457)
        fs1_var_2456 = ((specmult1_var_2457) - AINT((specmult1_var_2457)))
        speccomb_mco2_var_2410 = colh2o_var_2392(jl_var_2475, lay_var_2418) + refrat_m_a_var_2423 * coln2o_var_2393(jl_var_2475, lay_var_2418)
        specparm_mco2_var_2461 = MIN(colh2o_var_2392(jl_var_2475, lay_var_2418) / speccomb_mco2_var_2410, oneminus_var_2391)
        specmult_mco2_var_2460 = 8.0D0 * specparm_mco2_var_2461
        jmco2_var_2421 = 1 + INT(specmult_mco2_var_2460)
        fmco2_var_2459 = ((specmult_mco2_var_2460) - AINT((specmult_mco2_var_2460)))
        chi_co2_var_2449 = colco2_var_2394(jl_var_2475, lay_var_2418) / (coldry_var_2396(jl_var_2475, lay_var_2418))
        ratco2_var_2450 = 1D+20 * chi_co2_var_2449 / 0.000355D0
        IF (ratco2_var_2450 .GT. 3.0D0) THEN
          adjfac_var_2451 = 2.0D0 + (ratco2_var_2450 - 2.0D0) ** 0.68D0
          adjcolco2_var_2452 = adjfac_var_2451 * 0.000355D0 * coldry_var_2396(jl_var_2475, lay_var_2418) * 1D-20
        ELSE
          adjcolco2_var_2452 = colco2_var_2394(jl_var_2475, lay_var_2418)
        END IF
        speccomb_mco = colh2o_var_2392(jl_var_2475, lay_var_2418) + refrat_m_a3 * coln2o_var_2393(jl_var_2475, lay_var_2418)
        specparm_mco = MIN(colh2o_var_2392(jl_var_2475, lay_var_2418) / speccomb_mco, oneminus_var_2391)
        specmult_mco = 8.0D0 * specparm_mco
        jmco = 1 + INT(specmult_mco)
        fmco = ((specmult_mco) - AINT((specmult_mco)))
        speccomb_planck_var_2409 = colh2o_var_2392(jl_var_2475, lay_var_2418) + refrat_planck_a_var_2422 * coln2o_var_2393(jl_var_2475, lay_var_2418)
        specparm_planck_var_2464 = MIN(colh2o_var_2392(jl_var_2475, lay_var_2418) / speccomb_planck_var_2409, oneminus_var_2391)
        specmult_planck_var_2463 = 8.0D0 * specparm_planck_var_2464
        jpl_var_2420 = 1 + INT(specmult_planck_var_2463)
        fpl_var_2462 = ((specmult_planck_var_2463) - AINT((specmult_planck_var_2463)))
        ind0_var_2411 = ((jp_var_2388(jl_var_2475, lay_var_2418) - 1) * 5 + (jt_var_2389(jl_var_2475, lay_var_2418) - 1)) * nspa_var_237(13) + js_var_2417
        ind1_var_2412 = (jp_var_2388(jl_var_2475, lay_var_2418) * 5 + (jt1_var_2390(jl_var_2475, lay_var_2418) - 1)) * nspa_var_237(13) + js1_var_2419
        inds_var_2413 = indself_var_2400(jl_var_2475, lay_var_2418)
        indf_var_2414 = indfor_var_2402(jl_var_2475, lay_var_2418)
        indm_var_2415 = indminor_var_2406(jl_var_2475, lay_var_2418)
        IF (specparm_var_2455 .LT. 0.125D0) THEN
          p_var_2436 = fs_var_2453 - 1.0D0
          p4_var_2437 = p_var_2436 ** 4
          fk0_var_2438 = p4_var_2437
          fk1_var_2439 = 1.0D0 - p_var_2436 - 2.0D0 * p4_var_2437
          fk2_var_2440 = p_var_2436 + p4_var_2437
          fac000_var_2424 = fk0_var_2438 * fac00_var_2384(jl_var_2475, lay_var_2418)
          fac100_var_2425 = fk1_var_2439 * fac00_var_2384(jl_var_2475, lay_var_2418)
          fac200_var_2426 = fk2_var_2440 * fac00_var_2384(jl_var_2475, lay_var_2418)
          fac010_var_2427 = fk0_var_2438 * fac10_var_2386(jl_var_2475, lay_var_2418)
          fac110_var_2428 = fk1_var_2439 * fac10_var_2386(jl_var_2475, lay_var_2418)
          fac210_var_2429 = fk2_var_2440 * fac10_var_2386(jl_var_2475, lay_var_2418)
        ELSE IF (specparm_var_2455 .GT. 0.875D0) THEN
          p_var_2436 = - fs_var_2453
          p4_var_2437 = p_var_2436 ** 4
          fk0_var_2438 = p4_var_2437
          fk1_var_2439 = 1.0D0 - p_var_2436 - 2.0D0 * p4_var_2437
          fk2_var_2440 = p_var_2436 + p4_var_2437
          fac000_var_2424 = fk0_var_2438 * fac00_var_2384(jl_var_2475, lay_var_2418)
          fac100_var_2425 = fk1_var_2439 * fac00_var_2384(jl_var_2475, lay_var_2418)
          fac200_var_2426 = fk2_var_2440 * fac00_var_2384(jl_var_2475, lay_var_2418)
          fac010_var_2427 = fk0_var_2438 * fac10_var_2386(jl_var_2475, lay_var_2418)
          fac110_var_2428 = fk1_var_2439 * fac10_var_2386(jl_var_2475, lay_var_2418)
          fac210_var_2429 = fk2_var_2440 * fac10_var_2386(jl_var_2475, lay_var_2418)
        ELSE
          fac000_var_2424 = (1.0D0 - fs_var_2453) * fac00_var_2384(jl_var_2475, lay_var_2418)
          fac010_var_2427 = (1.0D0 - fs_var_2453) * fac10_var_2386(jl_var_2475, lay_var_2418)
          fac100_var_2425 = fs_var_2453 * fac00_var_2384(jl_var_2475, lay_var_2418)
          fac110_var_2428 = fs_var_2453 * fac10_var_2386(jl_var_2475, lay_var_2418)
          fac200_var_2426 = 0.0D0
          fac210_var_2429 = 0.0D0
        END IF
        IF (specparm1_var_2458 .LT. 0.125D0) THEN
          p_var_2436 = fs1_var_2456 - 1.0D0
          p4_var_2437 = p_var_2436 ** 4
          fk0_var_2438 = p4_var_2437
          fk1_var_2439 = 1.0D0 - p_var_2436 - 2.0D0 * p4_var_2437
          fk2_var_2440 = p_var_2436 + p4_var_2437
          fac001_var_2430 = fk0_var_2438 * fac01_var_2385(jl_var_2475, lay_var_2418)
          fac101_var_2431 = fk1_var_2439 * fac01_var_2385(jl_var_2475, lay_var_2418)
          fac201_var_2432 = fk2_var_2440 * fac01_var_2385(jl_var_2475, lay_var_2418)
          fac011_var_2433 = fk0_var_2438 * fac11_var_2387(jl_var_2475, lay_var_2418)
          fac111_var_2434 = fk1_var_2439 * fac11_var_2387(jl_var_2475, lay_var_2418)
          fac211_var_2435 = fk2_var_2440 * fac11_var_2387(jl_var_2475, lay_var_2418)
        ELSE IF (specparm1_var_2458 .GT. 0.875D0) THEN
          p_var_2436 = - fs1_var_2456
          p4_var_2437 = p_var_2436 ** 4
          fk0_var_2438 = p4_var_2437
          fk1_var_2439 = 1.0D0 - p_var_2436 - 2.0D0 * p4_var_2437
          fk2_var_2440 = p_var_2436 + p4_var_2437
          fac001_var_2430 = fk0_var_2438 * fac01_var_2385(jl_var_2475, lay_var_2418)
          fac101_var_2431 = fk1_var_2439 * fac01_var_2385(jl_var_2475, lay_var_2418)
          fac201_var_2432 = fk2_var_2440 * fac01_var_2385(jl_var_2475, lay_var_2418)
          fac011_var_2433 = fk0_var_2438 * fac11_var_2387(jl_var_2475, lay_var_2418)
          fac111_var_2434 = fk1_var_2439 * fac11_var_2387(jl_var_2475, lay_var_2418)
          fac211_var_2435 = fk2_var_2440 * fac11_var_2387(jl_var_2475, lay_var_2418)
        ELSE
          fac001_var_2430 = (1.0D0 - fs1_var_2456) * fac01_var_2385(jl_var_2475, lay_var_2418)
          fac011_var_2433 = (1.0D0 - fs1_var_2456) * fac11_var_2387(jl_var_2475, lay_var_2418)
          fac101_var_2431 = fs1_var_2456 * fac01_var_2385(jl_var_2475, lay_var_2418)
          fac111_var_2434 = fs1_var_2456 * fac11_var_2387(jl_var_2475, lay_var_2418)
          fac201_var_2432 = 0.0D0
          fac211_var_2435 = 0.0D0
        END IF
        IF (specparm_var_2455 .LT. 0.125D0) THEN
          tau_major_var_2443(1 : ng13) = speccomb_var_2407 * (fac000_var_2424 * absa_var_153(ind0_var_2411, 1 : 4) + fac100_var_2425 * absa_var_153(ind0_var_2411 + 1, 1 : 4) + fac200_var_2426 * absa_var_153(ind0_var_2411 + 2, 1 : 4) + fac010_var_2427 * absa_var_153(ind0_var_2411 + 9, 1 : 4) + fac110_var_2428 * absa_var_153(ind0_var_2411 + 10, 1 : 4) + fac210_var_2429 * absa_var_153(ind0_var_2411 + 11, 1 : 4))
        ELSE IF (specparm_var_2455 .GT. 0.875D0) THEN
          tau_major_var_2443(1 : ng13) = speccomb_var_2407 * (fac200_var_2426 * absa_var_153(ind0_var_2411 - 1, 1 : 4) + fac100_var_2425 * absa_var_153(ind0_var_2411, 1 : 4) + fac000_var_2424 * absa_var_153(ind0_var_2411 + 1, 1 : 4) + fac210_var_2429 * absa_var_153(ind0_var_2411 + 8, 1 : 4) + fac110_var_2428 * absa_var_153(ind0_var_2411 + 9, 1 : 4) + fac010_var_2427 * absa_var_153(ind0_var_2411 + 10, 1 : 4))
        ELSE
          tau_major_var_2443(1 : ng13) = speccomb_var_2407 * (fac000_var_2424 * absa_var_153(ind0_var_2411, 1 : 4) + fac100_var_2425 * absa_var_153(ind0_var_2411 + 1, 1 : 4) + fac010_var_2427 * absa_var_153(ind0_var_2411 + 9, 1 : 4) + fac110_var_2428 * absa_var_153(ind0_var_2411 + 10, 1 : 4))
        END IF
        IF (specparm1_var_2458 .LT. 0.125D0) THEN
          tau_major1_var_2444(1 : ng13) = speccomb1_var_2408 * (fac001_var_2430 * absa_var_153(ind1_var_2412, 1 : 4) + fac101_var_2431 * absa_var_153(ind1_var_2412 + 1, 1 : 4) + fac201_var_2432 * absa_var_153(ind1_var_2412 + 2, 1 : 4) + fac011_var_2433 * absa_var_153(ind1_var_2412 + 9, 1 : 4) + fac111_var_2434 * absa_var_153(ind1_var_2412 + 10, 1 : 4) + fac211_var_2435 * absa_var_153(ind1_var_2412 + 11, 1 : 4))
        ELSE IF (specparm1_var_2458 .GT. 0.875D0) THEN
          tau_major1_var_2444(1 : ng13) = speccomb1_var_2408 * (fac201_var_2432 * absa_var_153(ind1_var_2412 - 1, 1 : 4) + fac101_var_2431 * absa_var_153(ind1_var_2412, 1 : 4) + fac001_var_2430 * absa_var_153(ind1_var_2412 + 1, 1 : 4) + fac211_var_2435 * absa_var_153(ind1_var_2412 + 8, 1 : 4) + fac111_var_2434 * absa_var_153(ind1_var_2412 + 9, 1 : 4) + fac011_var_2433 * absa_var_153(ind1_var_2412 + 10, 1 : 4))
        ELSE
          tau_major1_var_2444(1 : ng13) = speccomb1_var_2408 * (fac001_var_2430 * absa_var_153(ind1_var_2412, 1 : 4) + fac101_var_2431 * absa_var_153(ind1_var_2412 + 1, 1 : 4) + fac011_var_2433 * absa_var_153(ind1_var_2412 + 9, 1 : 4) + fac111_var_2434 * absa_var_153(ind1_var_2412 + 10, 1 : 4))
        END IF
        DO ig_var_2416 = 1, 4
          tauself_var_2442 = selffac_var_2398(jl_var_2475, lay_var_2418) * (selfref_var_154(inds_var_2413, ig_var_2416) + selffrac_var_2399(jl_var_2475, lay_var_2418) * (selfref_var_154(inds_var_2413 + 1, ig_var_2416) - selfref_var_154(inds_var_2413, ig_var_2416)))
          taufor_var_2441 = forfac_var_2403(jl_var_2475, lay_var_2418) * (forref_var_155(indf_var_2414, ig_var_2416) + forfrac_var_2404(jl_var_2475, lay_var_2418) * (forref_var_155(indf_var_2414 + 1, ig_var_2416) - forref_var_155(indf_var_2414, ig_var_2416)))
          co2m1_var_2445 = ka_mco2_var_156(jmco2_var_2421, indm_var_2415, ig_var_2416) + fmco2_var_2459 * (ka_mco2_var_156(jmco2_var_2421 + 1, indm_var_2415, ig_var_2416) - ka_mco2_var_156(jmco2_var_2421, indm_var_2415, ig_var_2416))
          co2m2_var_2446 = ka_mco2_var_156(jmco2_var_2421, indm_var_2415 + 1, ig_var_2416) + fmco2_var_2459 * (ka_mco2_var_156(jmco2_var_2421 + 1, indm_var_2415 + 1, ig_var_2416) - ka_mco2_var_156(jmco2_var_2421, indm_var_2415 + 1, ig_var_2416))
          absco2_var_2447 = co2m1_var_2445 + minorfrac_var_2405(jl_var_2475, lay_var_2418) * (co2m2_var_2446 - co2m1_var_2445)
          com1 = ka_mco(jmco, indm_var_2415, ig_var_2416) + fmco * (ka_mco(jmco + 1, indm_var_2415, ig_var_2416) - ka_mco(jmco, indm_var_2415, ig_var_2416))
          com2 = ka_mco(jmco, indm_var_2415 + 1, ig_var_2416) + fmco * (ka_mco(jmco + 1, indm_var_2415 + 1, ig_var_2416) - ka_mco(jmco, indm_var_2415 + 1, ig_var_2416))
          absco = com1 + minorfrac_var_2405(jl_var_2475, lay_var_2418) * (com2 - com1)
          taug_var_2382(jl_var_2475, 130 + ig_var_2416, lay_var_2418) = tau_major_var_2443(ig_var_2416) + tau_major1_var_2444(ig_var_2416) + tauself_var_2442 + taufor_var_2441 + adjcolco2_var_2452 * absco2_var_2447 + colco(jl_var_2475, lay_var_2418) * absco
          fracs_var_2401(jl_var_2475, 130 + ig_var_2416, lay_var_2418) = fracrefa_var_151(ig_var_2416, jpl_var_2420) + fpl_var_2462 * (fracrefa_var_151(ig_var_2416, jpl_var_2420 + 1) - fracrefa_var_151(ig_var_2416, jpl_var_2420))
        END DO
      END DO
      ixc0_var_2472 = kfdia_var_2380 - kidia_var_2379 + 1 - ixc0_var_2472
      DO ixp_var_2473 = 1, ixc0_var_2472
        jl_var_2475 = ixhigh_var_2469(ixp_var_2473, lay_var_2418)
        indm_var_2415 = indminor_var_2406(jl_var_2475, lay_var_2418)
        DO ig_var_2416 = 1, 4
          abso3_var_2448 = kb_mo3(indm_var_2415, ig_var_2416) + minorfrac_var_2405(jl_var_2475, lay_var_2418) * (kb_mo3(indm_var_2415 + 1, ig_var_2416) - kb_mo3(indm_var_2415, ig_var_2416))
          taug_var_2382(jl_var_2475, 130 + ig_var_2416, lay_var_2418) = colo3_var_2395(jl_var_2475, lay_var_2418) * abso3_var_2448
          fracs_var_2401(jl_var_2475, 130 + ig_var_2416, lay_var_2418) = fracrefb_var_152(ig_var_2416)
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