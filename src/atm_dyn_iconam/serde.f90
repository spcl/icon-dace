MODULE f90_glue_vt_serde
  USE, INTRINSIC :: iso_c_binding
    USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
  USE mo_intp_data_strc, ONLY: t_int_state
  USE mo_model_domain, ONLY: t_grid_cells
  USE mo_model_domain, ONLY: t_grid_edges
  USE mo_model_domain, ONLY: t_grid_vertices
  USE mo_model_domain, ONLY: t_patch
  USE mo_nonhydro_types, ONLY: t_nh_diag
  USE mo_nonhydro_types, ONLY: t_nh_metrics
  USE mo_nonhydro_types, ONLY: t_nh_prog
  IMPLICIT NONE
  TYPE, BIND(C) :: glue_global_data_type
    REAL(KIND = c_double) :: m_divdamp_fac
    REAL(KIND = c_double) :: m_divdamp_fac_o2
    REAL(KIND = c_double) :: m_iau_wgt_dyn
    INTEGER(KIND = c_int) :: m_divdamp_order
    INTEGER(KIND = c_int) :: m_divdamp_type
    INTEGER(KIND = c_int) :: m_grf_intmethod_e
    INTEGER(KIND = c_int) :: m_i_am_accel_node
    INTEGER(KIND = c_int) :: m_iadv_rhotheta
    INTEGER(KIND = c_int) :: m_igradp_method
    INTEGER(KIND = c_int) :: m_is_iau_active
    INTEGER(KIND = c_int) :: m_itime_scheme
    INTEGER(KIND = c_int) :: m_l_limited_area
    INTEGER(KIND = c_int) :: m_ldeepatmo
    INTEGER(KIND = c_int) :: m_lextra_diffu
    INTEGER(KIND = c_int) :: m_lvert_nest
    INTEGER(KIND = c_int) :: m_nproma
    INTEGER(KIND = c_int) :: m_rayleigh_type
    INTEGER(KIND = c_int) :: m_timer_intp
    INTEGER(KIND = c_int) :: m_timer_solve_nh_cellcomp
    INTEGER(KIND = c_int) :: m_timer_solve_nh_edgecomp
    INTEGER(KIND = c_int) :: m_timer_solve_nh_veltend
    INTEGER(KIND = c_int) :: m_timer_solve_nh_vimpl
    INTEGER(KIND = c_int) :: m_timer_solve_nh_vnupd
    INTEGER(KIND = c_int) :: m_timers_level
    TYPE(c_ptr) :: m_kstart_dd3d
    TYPE(c_ptr) :: m_kstart_moist
    TYPE(c_ptr) :: m_ndyn_substeps_var
    TYPE(c_ptr) :: m_nflat_gradp
    TYPE(c_ptr) :: m_nflatlev
    TYPE(c_ptr) :: m_nrdmax
  END TYPE glue_global_data_type
  TYPE, BIND(C) :: glue_t_patch
    TYPE(c_ptr) :: m_cells
    TYPE(c_ptr) :: m_edges
    INTEGER(KIND = c_int) :: m_id
    INTEGER(KIND = c_int) :: m_n_childdom
    INTEGER(KIND = c_int) :: m_nblks_c
    INTEGER(KIND = c_int) :: m_nblks_e
    INTEGER(KIND = c_int) :: m_nblks_v
    INTEGER(KIND = c_int) :: m_nlev
    INTEGER(KIND = c_int) :: m_nlevp1
    INTEGER(KIND = c_int) :: m_nshift
    TYPE(c_ptr) :: m_verts
  END TYPE glue_t_patch
  TYPE, BIND(C) :: glue_t_int_state
    INTEGER(KIND = c_int) :: m___f2dace_SA_c_lin_e_d_0_s_25
    INTEGER(KIND = c_int) :: m___f2dace_SA_c_lin_e_d_1_s_26
    INTEGER(KIND = c_int) :: m___f2dace_SA_c_lin_e_d_2_s_27
    INTEGER(KIND = c_int) :: m___f2dace_SA_cells_aw_verts_d_0_s_31
    INTEGER(KIND = c_int) :: m___f2dace_SA_cells_aw_verts_d_1_s_32
    INTEGER(KIND = c_int) :: m___f2dace_SA_cells_aw_verts_d_2_s_33
    INTEGER(KIND = c_int) :: m___f2dace_SA_e_bln_c_s_d_0_s_28
    INTEGER(KIND = c_int) :: m___f2dace_SA_e_bln_c_s_d_1_s_29
    INTEGER(KIND = c_int) :: m___f2dace_SA_e_bln_c_s_d_2_s_30
    INTEGER(KIND = c_int) :: m___f2dace_SA_geofac_grdiv_d_0_s_37
    INTEGER(KIND = c_int) :: m___f2dace_SA_geofac_grdiv_d_1_s_38
    INTEGER(KIND = c_int) :: m___f2dace_SA_geofac_grdiv_d_2_s_39
    INTEGER(KIND = c_int) :: m___f2dace_SA_geofac_n2s_d_0_s_43
    INTEGER(KIND = c_int) :: m___f2dace_SA_geofac_n2s_d_1_s_44
    INTEGER(KIND = c_int) :: m___f2dace_SA_geofac_n2s_d_2_s_45
    INTEGER(KIND = c_int) :: m___f2dace_SA_geofac_rot_d_0_s_40
    INTEGER(KIND = c_int) :: m___f2dace_SA_geofac_rot_d_1_s_41
    INTEGER(KIND = c_int) :: m___f2dace_SA_geofac_rot_d_2_s_42
    INTEGER(KIND = c_int) :: m___f2dace_SA_rbf_vec_coeff_e_d_0_s_34
    INTEGER(KIND = c_int) :: m___f2dace_SA_rbf_vec_coeff_e_d_1_s_35
    INTEGER(KIND = c_int) :: m___f2dace_SA_rbf_vec_coeff_e_d_2_s_36
    INTEGER(KIND = c_int) :: m___f2dace_SOA_c_lin_e_d_0_s_25
    INTEGER(KIND = c_int) :: m___f2dace_SOA_c_lin_e_d_1_s_26
    INTEGER(KIND = c_int) :: m___f2dace_SOA_c_lin_e_d_2_s_27
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cells_aw_verts_d_0_s_31
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cells_aw_verts_d_1_s_32
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cells_aw_verts_d_2_s_33
    INTEGER(KIND = c_int) :: m___f2dace_SOA_e_bln_c_s_d_0_s_28
    INTEGER(KIND = c_int) :: m___f2dace_SOA_e_bln_c_s_d_1_s_29
    INTEGER(KIND = c_int) :: m___f2dace_SOA_e_bln_c_s_d_2_s_30
    INTEGER(KIND = c_int) :: m___f2dace_SOA_geofac_grdiv_d_0_s_37
    INTEGER(KIND = c_int) :: m___f2dace_SOA_geofac_grdiv_d_1_s_38
    INTEGER(KIND = c_int) :: m___f2dace_SOA_geofac_grdiv_d_2_s_39
    INTEGER(KIND = c_int) :: m___f2dace_SOA_geofac_n2s_d_0_s_43
    INTEGER(KIND = c_int) :: m___f2dace_SOA_geofac_n2s_d_1_s_44
    INTEGER(KIND = c_int) :: m___f2dace_SOA_geofac_n2s_d_2_s_45
    INTEGER(KIND = c_int) :: m___f2dace_SOA_geofac_rot_d_0_s_40
    INTEGER(KIND = c_int) :: m___f2dace_SOA_geofac_rot_d_1_s_41
    INTEGER(KIND = c_int) :: m___f2dace_SOA_geofac_rot_d_2_s_42
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rbf_vec_coeff_e_d_0_s_34
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rbf_vec_coeff_e_d_1_s_35
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rbf_vec_coeff_e_d_2_s_36
    TYPE(c_ptr) :: m_c_lin_e
    TYPE(c_ptr) :: m_cells_aw_verts
    TYPE(c_ptr) :: m_e_bln_c_s
    TYPE(c_ptr) :: m_e_flx_avg
    TYPE(c_ptr) :: m_geofac_div
    TYPE(c_ptr) :: m_geofac_grdiv
    TYPE(c_ptr) :: m_geofac_grg
    TYPE(c_ptr) :: m_geofac_n2s
    TYPE(c_ptr) :: m_geofac_rot
    TYPE(c_ptr) :: m_nudgecoeff_e
    TYPE(c_ptr) :: m_pos_on_tplane_e
    TYPE(c_ptr) :: m_rbf_vec_coeff_e
  END TYPE glue_t_int_state
  TYPE, BIND(C) :: glue_t_nh_prog
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_d_0_s_288
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_d_1_s_289
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_d_2_s_290
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_d_0_s_285
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_d_1_s_286
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_d_2_s_287
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_d_0_s_288
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_d_1_s_289
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_d_2_s_290
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_d_0_s_285
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_d_1_s_286
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_d_2_s_287
    TYPE(c_ptr) :: m_exner
    TYPE(c_ptr) :: m_rho
    TYPE(c_ptr) :: m_theta_v
    TYPE(c_ptr) :: m_vn
    TYPE(c_ptr) :: m_w
  END TYPE glue_t_nh_prog
  TYPE, BIND(C) :: glue_t_nh_metrics
    INTEGER(KIND = c_int) :: m___f2dace_SA_bdy_mflx_e_blk_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_bdy_mflx_e_idx_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff1_dwdz_d_0_s_332
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff1_dwdz_d_1_s_333
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff1_dwdz_d_2_s_334
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff2_dwdz_d_0_s_335
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff2_dwdz_d_1_s_336
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff2_dwdz_d_2_s_337
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff_gradekin_d_0_s_329
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff_gradekin_d_1_s_330
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff_gradekin_d_2_s_331
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff_gradp_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff_gradp_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff_gradp_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_coeff_gradp_d_3_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_d2dexdz2_fac1_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_d2dexdz2_fac1_mc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_d2dexdz2_fac1_mc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_d2dexdz2_fac2_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_d2dexdz2_fac2_mc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_d2dexdz2_fac2_mc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_d_exner_dz_ref_ic_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_d_exner_dz_ref_ic_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_d_exner_dz_ref_ic_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddqz_z_full_e_d_0_s_314
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddqz_z_full_e_d_1_s_315
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddqz_z_full_e_d_2_s_316
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddqz_z_half_d_0_s_317
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddqz_z_half_d_1_s_318
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddqz_z_half_d_2_s_319
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddxn_z_full_d_0_s_308
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddxn_z_full_d_1_s_309
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddxn_z_full_d_2_s_310
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddxt_z_full_d_0_s_311
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddxt_z_full_d_1_s_312
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddxt_z_full_d_2_s_313
    INTEGER(KIND = c_int) :: m___f2dace_SA_deepatmo_divh_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_deepatmo_divzl_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_deepatmo_divzu_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_deepatmo_gradh_ifc_d_0_s_340
    INTEGER(KIND = c_int) :: m___f2dace_SA_deepatmo_gradh_mc_d_0_s_338
    INTEGER(KIND = c_int) :: m___f2dace_SA_deepatmo_invr_ifc_d_0_s_341
    INTEGER(KIND = c_int) :: m___f2dace_SA_deepatmo_invr_mc_d_0_s_339
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_exfac_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_exfac_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_exfac_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_ref_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_ref_mc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_ref_mc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_hmask_dd3d_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_hmask_dd3d_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_inv_ddqz_z_full_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_inv_ddqz_z_full_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_inv_ddqz_z_full_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_pg_edgeblk_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_pg_edgeidx_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_pg_exdist_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_pg_vertidx_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rayleigh_vn_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rayleigh_w_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ref_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ref_mc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ref_mc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ref_me_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ref_me_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ref_me_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_scalfac_dd3d_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_ref_ic_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_ref_ic_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_ref_ic_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_ref_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_ref_mc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_ref_mc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_ref_me_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_ref_me_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_ref_me_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vertidx_gradp_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vertidx_gradp_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vertidx_gradp_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vertidx_gradp_d_3_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vwind_expl_wgt_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vwind_expl_wgt_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vwind_impl_wgt_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vwind_impl_wgt_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfac_c_d_0_s_320
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfac_c_d_1_s_321
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfac_c_d_2_s_322
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfac_e_d_0_s_323
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfac_e_d_1_s_324
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfac_e_d_2_s_325
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfacq1_c_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfacq1_c_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfacq1_c_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfacq_c_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfacq_c_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfacq_c_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfacq_e_d_0_s_326
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfacq_e_d_1_s_327
    INTEGER(KIND = c_int) :: m___f2dace_SA_wgtfacq_e_d_2_s_328
    INTEGER(KIND = c_int) :: m___f2dace_SA_zdiff_gradp_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_zdiff_gradp_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_zdiff_gradp_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_zdiff_gradp_d_3_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_bdy_mflx_e_blk_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_bdy_mflx_e_idx_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff1_dwdz_d_0_s_332
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff1_dwdz_d_1_s_333
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff1_dwdz_d_2_s_334
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff2_dwdz_d_0_s_335
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff2_dwdz_d_1_s_336
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff2_dwdz_d_2_s_337
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff_gradekin_d_0_s_329
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff_gradekin_d_1_s_330
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff_gradekin_d_2_s_331
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff_gradp_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff_gradp_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff_gradp_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_coeff_gradp_d_3_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_d2dexdz2_fac1_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_d2dexdz2_fac1_mc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_d2dexdz2_fac1_mc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_d2dexdz2_fac2_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_d2dexdz2_fac2_mc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_d2dexdz2_fac2_mc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_d_exner_dz_ref_ic_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_d_exner_dz_ref_ic_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_d_exner_dz_ref_ic_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddqz_z_full_e_d_0_s_314
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddqz_z_full_e_d_1_s_315
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddqz_z_full_e_d_2_s_316
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddqz_z_half_d_0_s_317
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddqz_z_half_d_1_s_318
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddqz_z_half_d_2_s_319
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddxn_z_full_d_0_s_308
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddxn_z_full_d_1_s_309
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddxn_z_full_d_2_s_310
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddxt_z_full_d_0_s_311
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddxt_z_full_d_1_s_312
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddxt_z_full_d_2_s_313
    INTEGER(KIND = c_int) :: m___f2dace_SOA_deepatmo_divh_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_deepatmo_divzl_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_deepatmo_divzu_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_deepatmo_gradh_ifc_d_0_s_340
    INTEGER(KIND = c_int) :: m___f2dace_SOA_deepatmo_gradh_mc_d_0_s_338
    INTEGER(KIND = c_int) :: m___f2dace_SOA_deepatmo_invr_ifc_d_0_s_341
    INTEGER(KIND = c_int) :: m___f2dace_SOA_deepatmo_invr_mc_d_0_s_339
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_exfac_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_exfac_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_exfac_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_ref_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_ref_mc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_ref_mc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_hmask_dd3d_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_hmask_dd3d_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_inv_ddqz_z_full_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_inv_ddqz_z_full_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_inv_ddqz_z_full_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_pg_edgeblk_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_pg_edgeidx_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_pg_exdist_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_pg_vertidx_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rayleigh_vn_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rayleigh_w_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ref_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ref_mc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ref_mc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ref_me_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ref_me_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ref_me_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_scalfac_dd3d_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_ref_ic_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_ref_ic_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_ref_ic_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_ref_mc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_ref_mc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_ref_mc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_ref_me_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_ref_me_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_ref_me_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vertidx_gradp_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vertidx_gradp_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vertidx_gradp_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vertidx_gradp_d_3_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vwind_expl_wgt_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vwind_expl_wgt_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vwind_impl_wgt_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vwind_impl_wgt_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfac_c_d_0_s_320
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfac_c_d_1_s_321
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfac_c_d_2_s_322
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfac_e_d_0_s_323
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfac_e_d_1_s_324
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfac_e_d_2_s_325
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfacq1_c_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfacq1_c_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfacq1_c_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfacq_c_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfacq_c_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfacq_c_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfacq_e_d_0_s_326
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfacq_e_d_1_s_327
    INTEGER(KIND = c_int) :: m___f2dace_SOA_wgtfacq_e_d_2_s_328
    INTEGER(KIND = c_int) :: m___f2dace_SOA_zdiff_gradp_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_zdiff_gradp_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_zdiff_gradp_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_zdiff_gradp_d_3_s
    TYPE(c_ptr) :: m_bdy_mflx_e_blk
    INTEGER(KIND = c_int) :: m_bdy_mflx_e_dim
    TYPE(c_ptr) :: m_bdy_mflx_e_idx
    TYPE(c_ptr) :: m_coeff1_dwdz
    TYPE(c_ptr) :: m_coeff2_dwdz
    TYPE(c_ptr) :: m_coeff_gradekin
    TYPE(c_ptr) :: m_coeff_gradp
    TYPE(c_ptr) :: m_d2dexdz2_fac1_mc
    TYPE(c_ptr) :: m_d2dexdz2_fac1_mc_d_1_s
    TYPE(c_ptr) :: m_d2dexdz2_fac1_mc_d_2_s
    TYPE(c_ptr) :: m_d2dexdz2_fac2_mc
    TYPE(c_ptr) :: m_d2dexdz2_fac2_mc_d_1_s
    TYPE(c_ptr) :: m_d2dexdz2_fac2_mc_d_2_s
    TYPE(c_ptr) :: m_d_exner_dz_ref_ic
    TYPE(c_ptr) :: m_ddqz_z_full_e
    TYPE(c_ptr) :: m_ddqz_z_half
    TYPE(c_ptr) :: m_ddxn_z_full
    TYPE(c_ptr) :: m_ddxt_z_full
    TYPE(c_ptr) :: m_deepatmo_divh_mc
    TYPE(c_ptr) :: m_deepatmo_divzl_mc
    TYPE(c_ptr) :: m_deepatmo_divzu_mc
    TYPE(c_ptr) :: m_deepatmo_gradh_ifc
    TYPE(c_ptr) :: m_deepatmo_gradh_mc
    TYPE(c_ptr) :: m_deepatmo_invr_ifc
    TYPE(c_ptr) :: m_deepatmo_invr_mc
    TYPE(c_ptr) :: m_exner_exfac
    TYPE(c_ptr) :: m_exner_ref_mc
    TYPE(c_ptr) :: m_hmask_dd3d
    TYPE(c_ptr) :: m_inv_ddqz_z_full
    TYPE(c_ptr) :: m_pg_edgeblk
    TYPE(c_ptr) :: m_pg_edgeidx
    TYPE(c_ptr) :: m_pg_exdist
    INTEGER(KIND = c_int) :: m_pg_listdim
    TYPE(c_ptr) :: m_pg_vertidx
    TYPE(c_ptr) :: m_rayleigh_vn
    TYPE(c_ptr) :: m_rayleigh_w
    TYPE(c_ptr) :: m_rho_ref_mc
    TYPE(c_ptr) :: m_rho_ref_me
    TYPE(c_ptr) :: m_scalfac_dd3d
    TYPE(c_ptr) :: m_theta_ref_ic
    TYPE(c_ptr) :: m_theta_ref_mc
    TYPE(c_ptr) :: m_theta_ref_me
    TYPE(c_ptr) :: m_vertidx_gradp
    TYPE(c_ptr) :: m_vwind_expl_wgt
    TYPE(c_ptr) :: m_vwind_impl_wgt
    TYPE(c_ptr) :: m_wgtfac_c
    TYPE(c_ptr) :: m_wgtfac_e
    TYPE(c_ptr) :: m_wgtfacq1_c
    TYPE(c_ptr) :: m_wgtfacq_c
    TYPE(c_ptr) :: m_wgtfacq_e
    TYPE(c_ptr) :: m_zdiff_gradp
  END TYPE glue_t_nh_metrics
  TYPE, BIND(C) :: glue_t_nh_diag
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_exner_phy_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_exner_phy_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_exner_phy_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_adv_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_adv_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_adv_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_apc_pc_d_0_s_300
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_apc_pc_d_1_s_301
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_apc_pc_d_2_s_302
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_apc_pc_d_3_s_303
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_cor_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_cor_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_cor_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_cor_pc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_cor_pc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_cor_pc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_cor_pc_d_3_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_dmp_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_dmp_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_dmp_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_dyn_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_dyn_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_dyn_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_grf_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_grf_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_grf_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_iau_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_iau_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_iau_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_pgr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_pgr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_pgr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_phd_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_phd_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_phd_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_phy_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_phy_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_phy_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_ray_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_ray_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_vn_ray_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_w_adv_pc_d_0_s_304
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_w_adv_pc_d_1_s_305
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_w_adv_pc_d_2_s_306
    INTEGER(KIND = c_int) :: m___f2dace_SA_ddt_w_adv_pc_d_3_s_307
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_dyn_incr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_dyn_incr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_dyn_incr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_incr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_incr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_incr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_pr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_pr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_exner_pr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_bdy_mflx_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_bdy_mflx_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_bdy_mflx_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_mflx_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_mflx_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_mflx_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_rho_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_rho_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_rho_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_thv_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_thv_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_thv_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_vn_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_vn_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_vn_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_w_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_w_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_grf_tend_w_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mass_fl_e_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mass_fl_e_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mass_fl_e_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mass_fl_e_sv_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mass_fl_e_sv_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mass_fl_e_sv_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mflx_ic_int_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mflx_ic_int_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mflx_ic_int_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mflx_ic_ubc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mflx_ic_ubc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_mflx_ic_ubc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ic_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ic_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ic_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ic_int_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ic_int_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ic_int_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ic_ubc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ic_ubc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_ic_ubc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_incr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_incr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_rho_incr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_ic_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_ic_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_ic_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_ic_int_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_ic_int_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_ic_int_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_ic_ubc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_ic_ubc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_theta_v_ic_ubc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_ie_d_0_s_294
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_ie_d_1_s_295
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_ie_d_2_s_296
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_ie_int_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_ie_int_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_ie_int_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_ie_ubc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_ie_ubc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_ie_ubc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_incr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_incr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vn_incr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_vt_d_0_s_291
    INTEGER(KIND = c_int) :: m___f2dace_SA_vt_d_1_s_292
    INTEGER(KIND = c_int) :: m___f2dace_SA_vt_d_2_s_293
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_concorr_c_d_0_s_297
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_concorr_c_d_1_s_298
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_concorr_c_d_2_s_299
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_int_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_int_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_int_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_ubc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_ubc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_w_ubc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_exner_phy_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_exner_phy_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_exner_phy_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_adv_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_adv_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_adv_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_apc_pc_d_0_s_300
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_apc_pc_d_1_s_301
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_apc_pc_d_2_s_302
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_apc_pc_d_3_s_303
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_cor_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_cor_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_cor_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_cor_pc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_cor_pc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_cor_pc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_cor_pc_d_3_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_dmp_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_dmp_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_dmp_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_dyn_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_dyn_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_dyn_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_grf_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_grf_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_grf_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_iau_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_iau_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_iau_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_pgr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_pgr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_pgr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_phd_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_phd_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_phd_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_phy_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_phy_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_phy_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_ray_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_ray_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_vn_ray_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_w_adv_pc_d_0_s_304
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_w_adv_pc_d_1_s_305
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_w_adv_pc_d_2_s_306
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ddt_w_adv_pc_d_3_s_307
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_dyn_incr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_dyn_incr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_dyn_incr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_incr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_incr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_incr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_pr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_pr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_exner_pr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_bdy_mflx_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_bdy_mflx_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_bdy_mflx_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_mflx_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_mflx_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_mflx_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_rho_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_rho_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_rho_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_thv_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_thv_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_thv_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_vn_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_vn_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_vn_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_w_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_w_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_grf_tend_w_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mass_fl_e_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mass_fl_e_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mass_fl_e_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mass_fl_e_sv_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mass_fl_e_sv_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mass_fl_e_sv_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mflx_ic_int_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mflx_ic_int_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mflx_ic_int_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mflx_ic_ubc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mflx_ic_ubc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_mflx_ic_ubc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ic_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ic_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ic_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ic_int_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ic_int_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ic_int_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ic_ubc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ic_ubc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_ic_ubc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_incr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_incr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_rho_incr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_ic_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_ic_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_ic_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_ic_int_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_ic_int_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_ic_int_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_ic_ubc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_ic_ubc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_theta_v_ic_ubc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_ie_d_0_s_294
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_ie_d_1_s_295
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_ie_d_2_s_296
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_ie_int_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_ie_int_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_ie_int_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_ie_ubc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_ie_ubc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_ie_ubc_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_incr_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_incr_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vn_incr_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vt_d_0_s_291
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vt_d_1_s_292
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vt_d_2_s_293
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_concorr_c_d_0_s_297
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_concorr_c_d_1_s_298
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_concorr_c_d_2_s_299
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_int_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_int_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_int_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_ubc_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_ubc_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_w_ubc_d_2_s
    TYPE(c_ptr) :: m_ddt_exner_phy
    TYPE(c_ptr) :: m_ddt_vn_adv
    INTEGER(KIND = c_int) :: m_ddt_vn_adv_is_associated
    TYPE(c_ptr) :: m_ddt_vn_apc_pc
    TYPE(c_ptr) :: m_ddt_vn_cor
    INTEGER(KIND = c_int) :: m_ddt_vn_cor_is_associated
    TYPE(c_ptr) :: m_ddt_vn_cor_pc
    TYPE(c_ptr) :: m_ddt_vn_dmp
    INTEGER(KIND = c_int) :: m_ddt_vn_dmp_is_associated
    TYPE(c_ptr) :: m_ddt_vn_dyn
    INTEGER(KIND = c_int) :: m_ddt_vn_dyn_is_associated
    TYPE(c_ptr) :: m_ddt_vn_grf
    INTEGER(KIND = c_int) :: m_ddt_vn_grf_is_associated
    TYPE(c_ptr) :: m_ddt_vn_iau
    INTEGER(KIND = c_int) :: m_ddt_vn_iau_is_associated
    TYPE(c_ptr) :: m_ddt_vn_pgr
    INTEGER(KIND = c_int) :: m_ddt_vn_pgr_is_associated
    TYPE(c_ptr) :: m_ddt_vn_phd
    INTEGER(KIND = c_int) :: m_ddt_vn_phd_is_associated
    TYPE(c_ptr) :: m_ddt_vn_phy
    TYPE(c_ptr) :: m_ddt_vn_ray
    INTEGER(KIND = c_int) :: m_ddt_vn_ray_is_associated
    TYPE(c_ptr) :: m_ddt_w_adv_pc
    TYPE(c_ptr) :: m_exner_dyn_incr
    TYPE(c_ptr) :: m_exner_incr
    TYPE(c_ptr) :: m_exner_pr
    TYPE(c_ptr) :: m_grf_bdy_mflx
    TYPE(c_ptr) :: m_grf_tend_mflx
    TYPE(c_ptr) :: m_grf_tend_rho
    TYPE(c_ptr) :: m_grf_tend_thv
    TYPE(c_ptr) :: m_grf_tend_vn
    TYPE(c_ptr) :: m_grf_tend_w
    TYPE(c_ptr) :: m_mass_fl_e
    TYPE(c_ptr) :: m_mass_fl_e_sv
    REAL(KIND = c_double) :: m_max_vcfl_dyn
    TYPE(c_ptr) :: m_mflx_ic_int
    TYPE(c_ptr) :: m_mflx_ic_ubc
    TYPE(c_ptr) :: m_rho_ic
    TYPE(c_ptr) :: m_rho_ic_int
    TYPE(c_ptr) :: m_rho_ic_ubc
    TYPE(c_ptr) :: m_rho_incr
    TYPE(c_ptr) :: m_theta_v_ic
    TYPE(c_ptr) :: m_theta_v_ic_int
    TYPE(c_ptr) :: m_theta_v_ic_ubc
    TYPE(c_ptr) :: m_vn_ie
    TYPE(c_ptr) :: m_vn_ie_int
    TYPE(c_ptr) :: m_vn_ie_ubc
    TYPE(c_ptr) :: m_vn_incr
    TYPE(c_ptr) :: m_vt
    TYPE(c_ptr) :: m_w_concorr_c
    TYPE(c_ptr) :: m_w_int
    TYPE(c_ptr) :: m_w_ubc
  END TYPE glue_t_nh_diag
  TYPE, BIND(C) :: glue_t_grid_edges
    INTEGER(KIND = c_int) :: m___f2dace_SA_area_edge_d_0_s_188
    INTEGER(KIND = c_int) :: m___f2dace_SA_area_edge_d_1_s_189
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_blk_d_0_s_167
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_blk_d_1_s_168
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_blk_d_2_s_169
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_idx_d_0_s_164
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_idx_d_1_s_165
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_idx_d_2_s_166
    INTEGER(KIND = c_int) :: m___f2dace_SA_dual_normal_cell_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_dual_normal_cell_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_dual_normal_cell_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_end_block_d_0_s_199
    INTEGER(KIND = c_int) :: m___f2dace_SA_end_index_d_0_s_197
    INTEGER(KIND = c_int) :: m___f2dace_SA_f_e_d_0_s_190
    INTEGER(KIND = c_int) :: m___f2dace_SA_f_e_d_1_s_191
    INTEGER(KIND = c_int) :: m___f2dace_SA_fn_e_d_0_s_192
    INTEGER(KIND = c_int) :: m___f2dace_SA_fn_e_d_1_s_193
    INTEGER(KIND = c_int) :: m___f2dace_SA_ft_e_d_0_s_194
    INTEGER(KIND = c_int) :: m___f2dace_SA_ft_e_d_1_s_195
    INTEGER(KIND = c_int) :: m___f2dace_SA_inv_dual_edge_length_d_0_s_186
    INTEGER(KIND = c_int) :: m___f2dace_SA_inv_dual_edge_length_d_1_s_187
    INTEGER(KIND = c_int) :: m___f2dace_SA_inv_primal_edge_length_d_0_s_184
    INTEGER(KIND = c_int) :: m___f2dace_SA_inv_primal_edge_length_d_1_s_185
    INTEGER(KIND = c_int) :: m___f2dace_SA_primal_normal_cell_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_primal_normal_cell_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_primal_normal_cell_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_quad_blk_d_0_s_181
    INTEGER(KIND = c_int) :: m___f2dace_SA_quad_blk_d_1_s_182
    INTEGER(KIND = c_int) :: m___f2dace_SA_quad_blk_d_2_s_183
    INTEGER(KIND = c_int) :: m___f2dace_SA_quad_idx_d_0_s_178
    INTEGER(KIND = c_int) :: m___f2dace_SA_quad_idx_d_1_s_179
    INTEGER(KIND = c_int) :: m___f2dace_SA_quad_idx_d_2_s_180
    INTEGER(KIND = c_int) :: m___f2dace_SA_refin_ctrl_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_refin_ctrl_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_start_block_d_0_s_198
    INTEGER(KIND = c_int) :: m___f2dace_SA_start_index_d_0_s_196
    INTEGER(KIND = c_int) :: m___f2dace_SA_tangent_orientation_d_0_s_176
    INTEGER(KIND = c_int) :: m___f2dace_SA_tangent_orientation_d_1_s_177
    INTEGER(KIND = c_int) :: m___f2dace_SA_vertex_blk_d_0_s_173
    INTEGER(KIND = c_int) :: m___f2dace_SA_vertex_blk_d_1_s_174
    INTEGER(KIND = c_int) :: m___f2dace_SA_vertex_blk_d_2_s_175
    INTEGER(KIND = c_int) :: m___f2dace_SA_vertex_idx_d_0_s_170
    INTEGER(KIND = c_int) :: m___f2dace_SA_vertex_idx_d_1_s_171
    INTEGER(KIND = c_int) :: m___f2dace_SA_vertex_idx_d_2_s_172
    INTEGER(KIND = c_int) :: m___f2dace_SOA_area_edge_d_0_s_188
    INTEGER(KIND = c_int) :: m___f2dace_SOA_area_edge_d_1_s_189
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_blk_d_0_s_167
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_blk_d_1_s_168
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_blk_d_2_s_169
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_idx_d_0_s_164
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_idx_d_1_s_165
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_idx_d_2_s_166
    INTEGER(KIND = c_int) :: m___f2dace_SOA_dual_normal_cell_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_dual_normal_cell_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_dual_normal_cell_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_end_block_d_0_s_199
    INTEGER(KIND = c_int) :: m___f2dace_SOA_end_index_d_0_s_197
    INTEGER(KIND = c_int) :: m___f2dace_SOA_f_e_d_0_s_190
    INTEGER(KIND = c_int) :: m___f2dace_SOA_f_e_d_1_s_191
    INTEGER(KIND = c_int) :: m___f2dace_SOA_fn_e_d_0_s_192
    INTEGER(KIND = c_int) :: m___f2dace_SOA_fn_e_d_1_s_193
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ft_e_d_0_s_194
    INTEGER(KIND = c_int) :: m___f2dace_SOA_ft_e_d_1_s_195
    INTEGER(KIND = c_int) :: m___f2dace_SOA_inv_dual_edge_length_d_0_s_186
    INTEGER(KIND = c_int) :: m___f2dace_SOA_inv_dual_edge_length_d_1_s_187
    INTEGER(KIND = c_int) :: m___f2dace_SOA_inv_primal_edge_length_d_0_s_184
    INTEGER(KIND = c_int) :: m___f2dace_SOA_inv_primal_edge_length_d_1_s_185
    INTEGER(KIND = c_int) :: m___f2dace_SOA_primal_normal_cell_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_primal_normal_cell_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_primal_normal_cell_d_2_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_quad_blk_d_0_s_181
    INTEGER(KIND = c_int) :: m___f2dace_SOA_quad_blk_d_1_s_182
    INTEGER(KIND = c_int) :: m___f2dace_SOA_quad_blk_d_2_s_183
    INTEGER(KIND = c_int) :: m___f2dace_SOA_quad_idx_d_0_s_178
    INTEGER(KIND = c_int) :: m___f2dace_SOA_quad_idx_d_1_s_179
    INTEGER(KIND = c_int) :: m___f2dace_SOA_quad_idx_d_2_s_180
    INTEGER(KIND = c_int) :: m___f2dace_SOA_refin_ctrl_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_refin_ctrl_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_start_block_d_0_s_198
    INTEGER(KIND = c_int) :: m___f2dace_SOA_start_index_d_0_s_196
    INTEGER(KIND = c_int) :: m___f2dace_SOA_tangent_orientation_d_0_s_176
    INTEGER(KIND = c_int) :: m___f2dace_SOA_tangent_orientation_d_1_s_177
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vertex_blk_d_0_s_173
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vertex_blk_d_1_s_174
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vertex_blk_d_2_s_175
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vertex_idx_d_0_s_170
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vertex_idx_d_1_s_171
    INTEGER(KIND = c_int) :: m___f2dace_SOA_vertex_idx_d_2_s_172
    TYPE(c_ptr) :: m_area_edge
    TYPE(c_ptr) :: m_cell_blk
    TYPE(c_ptr) :: m_cell_idx
    TYPE(c_ptr) :: m_dual_normal_cell
    TYPE(c_ptr) :: m_end_block
    TYPE(c_ptr) :: m_end_index
    TYPE(c_ptr) :: m_f_e
    TYPE(c_ptr) :: m_fn_e
    TYPE(c_ptr) :: m_ft_e
    TYPE(c_ptr) :: m_inv_dual_edge_length
    TYPE(c_ptr) :: m_inv_primal_edge_length
    TYPE(c_ptr) :: m_primal_normal_cell
    TYPE(c_ptr) :: m_quad_blk
    TYPE(c_ptr) :: m_quad_idx
    TYPE(c_ptr) :: m_refin_ctrl
    TYPE(c_ptr) :: m_start_block
    TYPE(c_ptr) :: m_start_index
    TYPE(c_ptr) :: m_tangent_orientation
    TYPE(c_ptr) :: m_vertex_blk
    TYPE(c_ptr) :: m_vertex_idx
  END TYPE glue_t_grid_edges
  TYPE, BIND(C) :: glue_t_grid_cells
    INTEGER(KIND = c_int) :: m___f2dace_SA_area_d_0_s_158
    INTEGER(KIND = c_int) :: m___f2dace_SA_area_d_1_s_159
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_blk_d_0_s_155
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_blk_d_1_s_156
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_blk_d_2_s_157
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_idx_d_0_s_152
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_idx_d_1_s_153
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_idx_d_2_s_154
    INTEGER(KIND = c_int) :: m___f2dace_SA_end_blk_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_end_blk_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_end_block_d_0_s_163
    INTEGER(KIND = c_int) :: m___f2dace_SA_end_index_d_0_s_161
    INTEGER(KIND = c_int) :: m___f2dace_SA_neighbor_blk_d_0_s_149
    INTEGER(KIND = c_int) :: m___f2dace_SA_neighbor_blk_d_1_s_150
    INTEGER(KIND = c_int) :: m___f2dace_SA_neighbor_blk_d_2_s_151
    INTEGER(KIND = c_int) :: m___f2dace_SA_neighbor_idx_d_0_s_146
    INTEGER(KIND = c_int) :: m___f2dace_SA_neighbor_idx_d_1_s_147
    INTEGER(KIND = c_int) :: m___f2dace_SA_neighbor_idx_d_2_s_148
    INTEGER(KIND = c_int) :: m___f2dace_SA_start_blk_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_start_blk_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SA_start_block_d_0_s_162
    INTEGER(KIND = c_int) :: m___f2dace_SA_start_index_d_0_s_160
    INTEGER(KIND = c_int) :: m___f2dace_SOA_area_d_0_s_158
    INTEGER(KIND = c_int) :: m___f2dace_SOA_area_d_1_s_159
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_blk_d_0_s_155
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_blk_d_1_s_156
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_blk_d_2_s_157
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_idx_d_0_s_152
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_idx_d_1_s_153
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_idx_d_2_s_154
    INTEGER(KIND = c_int) :: m___f2dace_SOA_end_blk_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_end_blk_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_end_block_d_0_s_163
    INTEGER(KIND = c_int) :: m___f2dace_SOA_end_index_d_0_s_161
    INTEGER(KIND = c_int) :: m___f2dace_SOA_neighbor_blk_d_0_s_149
    INTEGER(KIND = c_int) :: m___f2dace_SOA_neighbor_blk_d_1_s_150
    INTEGER(KIND = c_int) :: m___f2dace_SOA_neighbor_blk_d_2_s_151
    INTEGER(KIND = c_int) :: m___f2dace_SOA_neighbor_idx_d_0_s_146
    INTEGER(KIND = c_int) :: m___f2dace_SOA_neighbor_idx_d_1_s_147
    INTEGER(KIND = c_int) :: m___f2dace_SOA_neighbor_idx_d_2_s_148
    INTEGER(KIND = c_int) :: m___f2dace_SOA_start_blk_d_0_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_start_blk_d_1_s
    INTEGER(KIND = c_int) :: m___f2dace_SOA_start_block_d_0_s_162
    INTEGER(KIND = c_int) :: m___f2dace_SOA_start_index_d_0_s_160
    TYPE(c_ptr) :: m_area
    TYPE(c_ptr) :: m_decomp_info
    TYPE(c_ptr) :: m_edge_blk
    TYPE(c_ptr) :: m_edge_idx
    TYPE(c_ptr) :: m_end_blk
    TYPE(c_ptr) :: m_end_block
    TYPE(c_ptr) :: m_end_index
    TYPE(c_ptr) :: m_neighbor_blk
    TYPE(c_ptr) :: m_neighbor_idx
    TYPE(c_ptr) :: m_start_blk
    TYPE(c_ptr) :: m_start_block
    TYPE(c_ptr) :: m_start_index
  END TYPE glue_t_grid_cells
  TYPE, BIND(C) :: glue_t_grid_domain_decomp_info
    INTEGER(KIND = c_int) :: m___f2dace_SA_owner_mask_d_0_s_2
    INTEGER(KIND = c_int) :: m___f2dace_SA_owner_mask_d_1_s_3
    INTEGER(KIND = c_int) :: m___f2dace_SOA_owner_mask_d_0_s_2
    INTEGER(KIND = c_int) :: m___f2dace_SOA_owner_mask_d_1_s_3
    TYPE(c_ptr) :: m_owner_mask
  END TYPE glue_t_grid_domain_decomp_info
  TYPE, BIND(C) :: glue_t_grid_vertices
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_blk_d_0_s_167
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_blk_d_1_s_168
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_blk_d_2_s_169
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_idx_d_0_s_164
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_idx_d_1_s_165
    INTEGER(KIND = c_int) :: m___f2dace_SA_cell_idx_d_2_s_166
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_blk_d_0_s_155
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_blk_d_1_s_156
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_blk_d_2_s_157
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_idx_d_0_s_152
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_idx_d_1_s_153
    INTEGER(KIND = c_int) :: m___f2dace_SA_edge_idx_d_2_s_154
    INTEGER(KIND = c_int) :: m___f2dace_SA_end_block_d_0_s_163
    INTEGER(KIND = c_int) :: m___f2dace_SA_end_index_d_0_s_161
    INTEGER(KIND = c_int) :: m___f2dace_SA_start_block_d_0_s_162
    INTEGER(KIND = c_int) :: m___f2dace_SA_start_index_d_0_s_160
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_blk_d_0_s_167
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_blk_d_1_s_168
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_blk_d_2_s_169
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_idx_d_0_s_164
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_idx_d_1_s_165
    INTEGER(KIND = c_int) :: m___f2dace_SOA_cell_idx_d_2_s_166
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_blk_d_0_s_155
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_blk_d_1_s_156
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_blk_d_2_s_157
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_idx_d_0_s_152
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_idx_d_1_s_153
    INTEGER(KIND = c_int) :: m___f2dace_SOA_edge_idx_d_2_s_154
    INTEGER(KIND = c_int) :: m___f2dace_SOA_end_block_d_0_s_163
    INTEGER(KIND = c_int) :: m___f2dace_SOA_end_index_d_0_s_161
    INTEGER(KIND = c_int) :: m___f2dace_SOA_start_block_d_0_s_162
    INTEGER(KIND = c_int) :: m___f2dace_SOA_start_index_d_0_s_160
    TYPE(c_ptr) :: m_cell_blk
    TYPE(c_ptr) :: m_cell_idx
    TYPE(c_ptr) :: m_edge_blk
    TYPE(c_ptr) :: m_edge_idx
    TYPE(c_ptr) :: m_end_block
    TYPE(c_ptr) :: m_end_index
    TYPE(c_ptr) :: m_start_block
    TYPE(c_ptr) :: m_start_index
  END TYPE glue_t_grid_vertices
  INTERFACE ctor
    MODULE PROCEDURE :: ctor_t_patch
    MODULE PROCEDURE :: ctor_t_int_state
    MODULE PROCEDURE :: ctor_t_nh_prog
    MODULE PROCEDURE :: ctor_t_nh_metrics
    MODULE PROCEDURE :: ctor_t_nh_diag
    MODULE PROCEDURE :: ctor_t_grid_edges
    MODULE PROCEDURE :: ctor_t_grid_cells
    MODULE PROCEDURE :: ctor_t_grid_domain_decomp_info
    MODULE PROCEDURE :: ctor_t_grid_vertices
  END INTERFACE ctor
  CONTAINS
  SUBROUTINE ctor_t_patch(inp, out, initalloc)
    TYPE(t_patch), INTENT(IN) :: inp
    TYPE(glue_t_patch), INTENT(INOUT) :: out
    LOGICAL, INTENT(IN) :: initalloc
    TYPE(glue_t_grid_cells), ALLOCATABLE, TARGET, SAVE :: a_cells
    TYPE(glue_t_grid_edges), ALLOCATABLE, TARGET, SAVE :: a_edges
    TYPE(glue_t_grid_vertices), ALLOCATABLE, TARGET, SAVE :: a_verts
    out % m_id = 0
    out % m_n_childdom = 0
    out % m_nlev = 0
    out % m_nlevp1 = 0
    out % m_nshift = 0
    IF (initalloc .AND. .NOT. ALLOCATED(a_cells)) ALLOCATE(a_cells)
    CALL ctor(inp % cells, a_cells, initalloc)
    out % m_cells = c_loc(a_cells)
    IF (initalloc .AND. .NOT. ALLOCATED(a_edges)) ALLOCATE(a_edges)
    CALL ctor(inp % edges, a_edges, initalloc)
    out % m_edges = c_loc(a_edges)
    out % m_nblks_c = inp % nblks_c
    out % m_nblks_e = inp % nblks_e
    out % m_nblks_v = inp % nblks_v
    IF (initalloc .AND. .NOT. ALLOCATED(a_verts)) ALLOCATE(a_verts)
    CALL ctor(inp % verts, a_verts, initalloc)
    out % m_verts = c_loc(a_verts)
  END SUBROUTINE ctor_t_patch
  SUBROUTINE ctor_t_int_state(inp, out, initalloc)
    TYPE(t_int_state), INTENT(IN) :: inp
    TYPE(glue_t_int_state), INTENT(INOUT) :: out
    LOGICAL, INTENT(IN) :: initalloc
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_c_lin_e(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_cells_aw_verts(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_e_bln_c_s(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_geofac_grdiv(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_geofac_n2s(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_geofac_rot(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_rbf_vec_coeff_e(:, :, :)
    out % m_e_flx_avg = c_null_ptr
    out % m_geofac_div = c_null_ptr
    out % m_geofac_grg = c_null_ptr
    out % m_nudgecoeff_e = c_null_ptr
    out % m_pos_on_tplane_e = c_null_ptr
    IF (initalloc .AND. .NOT. ALLOCATED(a_c_lin_e)) ALLOCATE(a_c_lin_e(SIZE(inp % c_lin_e, 1), SIZE(inp % c_lin_e, 2), SIZE(inp % c_lin_e, 3)))
    a_c_lin_e = inp % c_lin_e
    out % m_c_lin_e = c_loc(a_c_lin_e)
    out % m___f2dace_SA_c_lin_e_d_0_s_25 = SIZE(inp % c_lin_e, 1)
    out % m___f2dace_SA_c_lin_e_d_1_s_26 = SIZE(inp % c_lin_e, 2)
    out % m___f2dace_SA_c_lin_e_d_2_s_27 = SIZE(inp % c_lin_e, 3)
    out % m___f2dace_SOA_c_lin_e_d_0_s_25 = LBOUND(inp % c_lin_e, 1)
    out % m___f2dace_SOA_c_lin_e_d_1_s_26 = LBOUND(inp % c_lin_e, 2)
    out % m___f2dace_SOA_c_lin_e_d_2_s_27 = LBOUND(inp % c_lin_e, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_cells_aw_verts)) ALLOCATE(a_cells_aw_verts(SIZE(inp % cells_aw_verts, 1), SIZE(inp % cells_aw_verts, 2), SIZE(inp % cells_aw_verts, 3)))
    a_cells_aw_verts = inp % cells_aw_verts
    out % m_cells_aw_verts = c_loc(a_cells_aw_verts)
    out % m___f2dace_SA_cells_aw_verts_d_0_s_31 = SIZE(inp % cells_aw_verts, 1)
    out % m___f2dace_SA_cells_aw_verts_d_1_s_32 = SIZE(inp % cells_aw_verts, 2)
    out % m___f2dace_SA_cells_aw_verts_d_2_s_33 = SIZE(inp % cells_aw_verts, 3)
    out % m___f2dace_SOA_cells_aw_verts_d_0_s_31 = LBOUND(inp % cells_aw_verts, 1)
    out % m___f2dace_SOA_cells_aw_verts_d_1_s_32 = LBOUND(inp % cells_aw_verts, 2)
    out % m___f2dace_SOA_cells_aw_verts_d_2_s_33 = LBOUND(inp % cells_aw_verts, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_e_bln_c_s)) ALLOCATE(a_e_bln_c_s(SIZE(inp % e_bln_c_s, 1), SIZE(inp % e_bln_c_s, 2), SIZE(inp % e_bln_c_s, 3)))
    a_e_bln_c_s = inp % e_bln_c_s
    out % m_e_bln_c_s = c_loc(a_e_bln_c_s)
    out % m___f2dace_SA_e_bln_c_s_d_0_s_28 = SIZE(inp % e_bln_c_s, 1)
    out % m___f2dace_SA_e_bln_c_s_d_1_s_29 = SIZE(inp % e_bln_c_s, 2)
    out % m___f2dace_SA_e_bln_c_s_d_2_s_30 = SIZE(inp % e_bln_c_s, 3)
    out % m___f2dace_SOA_e_bln_c_s_d_0_s_28 = LBOUND(inp % e_bln_c_s, 1)
    out % m___f2dace_SOA_e_bln_c_s_d_1_s_29 = LBOUND(inp % e_bln_c_s, 2)
    out % m___f2dace_SOA_e_bln_c_s_d_2_s_30 = LBOUND(inp % e_bln_c_s, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_geofac_grdiv)) ALLOCATE(a_geofac_grdiv(SIZE(inp % geofac_grdiv, 1), SIZE(inp % geofac_grdiv, 2), SIZE(inp % geofac_grdiv, 3)))
    a_geofac_grdiv = inp % geofac_grdiv
    out % m_geofac_grdiv = c_loc(a_geofac_grdiv)
    out % m___f2dace_SA_geofac_grdiv_d_0_s_37 = SIZE(inp % geofac_grdiv, 1)
    out % m___f2dace_SA_geofac_grdiv_d_1_s_38 = SIZE(inp % geofac_grdiv, 2)
    out % m___f2dace_SA_geofac_grdiv_d_2_s_39 = SIZE(inp % geofac_grdiv, 3)
    out % m___f2dace_SOA_geofac_grdiv_d_0_s_37 = LBOUND(inp % geofac_grdiv, 1)
    out % m___f2dace_SOA_geofac_grdiv_d_1_s_38 = LBOUND(inp % geofac_grdiv, 2)
    out % m___f2dace_SOA_geofac_grdiv_d_2_s_39 = LBOUND(inp % geofac_grdiv, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_geofac_n2s)) ALLOCATE(a_geofac_n2s(SIZE(inp % geofac_n2s, 1), SIZE(inp % geofac_n2s, 2), SIZE(inp % geofac_n2s, 3)))
    a_geofac_n2s = inp % geofac_n2s
    out % m_geofac_n2s = c_loc(a_geofac_n2s)
    out % m___f2dace_SA_geofac_n2s_d_0_s_43 = SIZE(inp % geofac_n2s, 1)
    out % m___f2dace_SA_geofac_n2s_d_1_s_44 = SIZE(inp % geofac_n2s, 2)
    out % m___f2dace_SA_geofac_n2s_d_2_s_45 = SIZE(inp % geofac_n2s, 3)
    out % m___f2dace_SOA_geofac_n2s_d_0_s_43 = LBOUND(inp % geofac_n2s, 1)
    out % m___f2dace_SOA_geofac_n2s_d_1_s_44 = LBOUND(inp % geofac_n2s, 2)
    out % m___f2dace_SOA_geofac_n2s_d_2_s_45 = LBOUND(inp % geofac_n2s, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_geofac_rot)) ALLOCATE(a_geofac_rot(SIZE(inp % geofac_rot, 1), SIZE(inp % geofac_rot, 2), SIZE(inp % geofac_rot, 3)))
    a_geofac_rot = inp % geofac_rot
    out % m_geofac_rot = c_loc(a_geofac_rot)
    out % m___f2dace_SA_geofac_rot_d_0_s_40 = SIZE(inp % geofac_rot, 1)
    out % m___f2dace_SA_geofac_rot_d_1_s_41 = SIZE(inp % geofac_rot, 2)
    out % m___f2dace_SA_geofac_rot_d_2_s_42 = SIZE(inp % geofac_rot, 3)
    out % m___f2dace_SOA_geofac_rot_d_0_s_40 = LBOUND(inp % geofac_rot, 1)
    out % m___f2dace_SOA_geofac_rot_d_1_s_41 = LBOUND(inp % geofac_rot, 2)
    out % m___f2dace_SOA_geofac_rot_d_2_s_42 = LBOUND(inp % geofac_rot, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_rbf_vec_coeff_e)) ALLOCATE(a_rbf_vec_coeff_e(SIZE(inp % rbf_vec_coeff_e, 1), SIZE(inp % rbf_vec_coeff_e, 2), SIZE(inp % rbf_vec_coeff_e, 3)))
    a_rbf_vec_coeff_e = inp % rbf_vec_coeff_e
    out % m_rbf_vec_coeff_e = c_loc(a_rbf_vec_coeff_e)
    out % m___f2dace_SA_rbf_vec_coeff_e_d_0_s_34 = SIZE(inp % rbf_vec_coeff_e, 1)
    out % m___f2dace_SA_rbf_vec_coeff_e_d_1_s_35 = SIZE(inp % rbf_vec_coeff_e, 2)
    out % m___f2dace_SA_rbf_vec_coeff_e_d_2_s_36 = SIZE(inp % rbf_vec_coeff_e, 3)
    out % m___f2dace_SOA_rbf_vec_coeff_e_d_0_s_34 = LBOUND(inp % rbf_vec_coeff_e, 1)
    out % m___f2dace_SOA_rbf_vec_coeff_e_d_1_s_35 = LBOUND(inp % rbf_vec_coeff_e, 2)
    out % m___f2dace_SOA_rbf_vec_coeff_e_d_2_s_36 = LBOUND(inp % rbf_vec_coeff_e, 3)
  END SUBROUTINE ctor_t_int_state
  SUBROUTINE ctor_t_nh_prog(inp, out, initalloc)
    TYPE(t_nh_prog), INTENT(IN) :: inp
    TYPE(glue_t_nh_prog), INTENT(INOUT) :: out
    LOGICAL, INTENT(IN) :: initalloc
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_vn(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_w(:, :, :)
    out % m___f2dace_SA_exner_d_0_s = 0
    out % m___f2dace_SA_exner_d_1_s = 0
    out % m___f2dace_SA_exner_d_2_s = 0
    out % m___f2dace_SA_rho_d_0_s = 0
    out % m___f2dace_SA_rho_d_1_s = 0
    out % m___f2dace_SA_rho_d_2_s = 0
    out % m___f2dace_SA_theta_v_d_0_s = 0
    out % m___f2dace_SA_theta_v_d_1_s = 0
    out % m___f2dace_SA_theta_v_d_2_s = 0
    out % m___f2dace_SOA_exner_d_0_s = 0
    out % m___f2dace_SOA_exner_d_1_s = 0
    out % m___f2dace_SOA_exner_d_2_s = 0
    out % m___f2dace_SOA_rho_d_0_s = 0
    out % m___f2dace_SOA_rho_d_1_s = 0
    out % m___f2dace_SOA_rho_d_2_s = 0
    out % m___f2dace_SOA_theta_v_d_0_s = 0
    out % m___f2dace_SOA_theta_v_d_1_s = 0
    out % m___f2dace_SOA_theta_v_d_2_s = 0
    out % m_exner = c_null_ptr
    out % m_rho = c_null_ptr
    out % m_theta_v = c_null_ptr
    IF (initalloc .AND. .NOT. ALLOCATED(a_vn)) ALLOCATE(a_vn(SIZE(inp % vn, 1), SIZE(inp % vn, 2), SIZE(inp % vn, 3)))
    a_vn = inp % vn
    out % m_vn = c_loc(a_vn)
    out % m___f2dace_SA_vn_d_0_s_288 = SIZE(inp % vn, 1)
    out % m___f2dace_SA_vn_d_1_s_289 = SIZE(inp % vn, 2)
    out % m___f2dace_SA_vn_d_2_s_290 = SIZE(inp % vn, 3)
    out % m___f2dace_SOA_vn_d_0_s_288 = LBOUND(inp % vn, 1)
    out % m___f2dace_SOA_vn_d_1_s_289 = LBOUND(inp % vn, 2)
    out % m___f2dace_SOA_vn_d_2_s_290 = LBOUND(inp % vn, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_w)) ALLOCATE(a_w(SIZE(inp % w, 1), SIZE(inp % w, 2), SIZE(inp % w, 3)))
    a_w = inp % w
    out % m_w = c_loc(a_w)
    out % m___f2dace_SA_w_d_0_s_285 = SIZE(inp % w, 1)
    out % m___f2dace_SA_w_d_1_s_286 = SIZE(inp % w, 2)
    out % m___f2dace_SA_w_d_2_s_287 = SIZE(inp % w, 3)
    out % m___f2dace_SOA_w_d_0_s_285 = LBOUND(inp % w, 1)
    out % m___f2dace_SOA_w_d_1_s_286 = LBOUND(inp % w, 2)
    out % m___f2dace_SOA_w_d_2_s_287 = LBOUND(inp % w, 3)
  END SUBROUTINE ctor_t_nh_prog
  SUBROUTINE ctor_t_nh_metrics(inp, out, initalloc)
    TYPE(t_nh_metrics), INTENT(IN) :: inp
    TYPE(glue_t_nh_metrics), INTENT(INOUT) :: out
    LOGICAL, INTENT(IN) :: initalloc
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_coeff1_dwdz(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_coeff2_dwdz(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_coeff_gradekin(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_ddqz_z_full_e(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_ddqz_z_half(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_ddxn_z_full(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_ddxt_z_full(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_deepatmo_gradh_ifc(:)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_deepatmo_gradh_mc(:)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_deepatmo_invr_ifc(:)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_deepatmo_invr_mc(:)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_wgtfac_c(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_wgtfac_e(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_wgtfacq_e(:, :, :)
    out % m___f2dace_SA_bdy_mflx_e_blk_d_0_s = 0
    out % m___f2dace_SA_bdy_mflx_e_idx_d_0_s = 0
    out % m___f2dace_SA_coeff_gradp_d_0_s = 0
    out % m___f2dace_SA_coeff_gradp_d_1_s = 0
    out % m___f2dace_SA_coeff_gradp_d_2_s = 0
    out % m___f2dace_SA_coeff_gradp_d_3_s = 0
    out % m___f2dace_SA_d2dexdz2_fac1_mc_d_0_s = 0
    out % m___f2dace_SA_d2dexdz2_fac1_mc_d_1_s = 0
    out % m___f2dace_SA_d2dexdz2_fac1_mc_d_2_s = 0
    out % m___f2dace_SA_d2dexdz2_fac2_mc_d_0_s = 0
    out % m___f2dace_SA_d2dexdz2_fac2_mc_d_1_s = 0
    out % m___f2dace_SA_d2dexdz2_fac2_mc_d_2_s = 0
    out % m___f2dace_SA_d_exner_dz_ref_ic_d_0_s = 0
    out % m___f2dace_SA_d_exner_dz_ref_ic_d_1_s = 0
    out % m___f2dace_SA_d_exner_dz_ref_ic_d_2_s = 0
    out % m___f2dace_SA_deepatmo_divh_mc_d_0_s = 0
    out % m___f2dace_SA_deepatmo_divzl_mc_d_0_s = 0
    out % m___f2dace_SA_deepatmo_divzu_mc_d_0_s = 0
    out % m___f2dace_SA_exner_exfac_d_0_s = 0
    out % m___f2dace_SA_exner_exfac_d_1_s = 0
    out % m___f2dace_SA_exner_exfac_d_2_s = 0
    out % m___f2dace_SA_exner_ref_mc_d_0_s = 0
    out % m___f2dace_SA_exner_ref_mc_d_1_s = 0
    out % m___f2dace_SA_exner_ref_mc_d_2_s = 0
    out % m___f2dace_SA_hmask_dd3d_d_0_s = 0
    out % m___f2dace_SA_hmask_dd3d_d_1_s = 0
    out % m___f2dace_SA_inv_ddqz_z_full_d_0_s = 0
    out % m___f2dace_SA_inv_ddqz_z_full_d_1_s = 0
    out % m___f2dace_SA_inv_ddqz_z_full_d_2_s = 0
    out % m___f2dace_SA_pg_edgeblk_d_0_s = 0
    out % m___f2dace_SA_pg_edgeidx_d_0_s = 0
    out % m___f2dace_SA_pg_exdist_d_0_s = 0
    out % m___f2dace_SA_pg_vertidx_d_0_s = 0
    out % m___f2dace_SA_rayleigh_vn_d_0_s = 0
    out % m___f2dace_SA_rayleigh_w_d_0_s = 0
    out % m___f2dace_SA_rho_ref_mc_d_0_s = 0
    out % m___f2dace_SA_rho_ref_mc_d_1_s = 0
    out % m___f2dace_SA_rho_ref_mc_d_2_s = 0
    out % m___f2dace_SA_rho_ref_me_d_0_s = 0
    out % m___f2dace_SA_rho_ref_me_d_1_s = 0
    out % m___f2dace_SA_rho_ref_me_d_2_s = 0
    out % m___f2dace_SA_scalfac_dd3d_d_0_s = 0
    out % m___f2dace_SA_theta_ref_ic_d_0_s = 0
    out % m___f2dace_SA_theta_ref_ic_d_1_s = 0
    out % m___f2dace_SA_theta_ref_ic_d_2_s = 0
    out % m___f2dace_SA_theta_ref_mc_d_0_s = 0
    out % m___f2dace_SA_theta_ref_mc_d_1_s = 0
    out % m___f2dace_SA_theta_ref_mc_d_2_s = 0
    out % m___f2dace_SA_theta_ref_me_d_0_s = 0
    out % m___f2dace_SA_theta_ref_me_d_1_s = 0
    out % m___f2dace_SA_theta_ref_me_d_2_s = 0
    out % m___f2dace_SA_vertidx_gradp_d_0_s = 0
    out % m___f2dace_SA_vertidx_gradp_d_1_s = 0
    out % m___f2dace_SA_vertidx_gradp_d_2_s = 0
    out % m___f2dace_SA_vertidx_gradp_d_3_s = 0
    out % m___f2dace_SA_vwind_expl_wgt_d_0_s = 0
    out % m___f2dace_SA_vwind_expl_wgt_d_1_s = 0
    out % m___f2dace_SA_vwind_impl_wgt_d_0_s = 0
    out % m___f2dace_SA_vwind_impl_wgt_d_1_s = 0
    out % m___f2dace_SA_wgtfacq1_c_d_0_s = 0
    out % m___f2dace_SA_wgtfacq1_c_d_1_s = 0
    out % m___f2dace_SA_wgtfacq1_c_d_2_s = 0
    out % m___f2dace_SA_wgtfacq_c_d_0_s = 0
    out % m___f2dace_SA_wgtfacq_c_d_1_s = 0
    out % m___f2dace_SA_wgtfacq_c_d_2_s = 0
    out % m___f2dace_SA_zdiff_gradp_d_0_s = 0
    out % m___f2dace_SA_zdiff_gradp_d_1_s = 0
    out % m___f2dace_SA_zdiff_gradp_d_2_s = 0
    out % m___f2dace_SA_zdiff_gradp_d_3_s = 0
    out % m___f2dace_SOA_bdy_mflx_e_blk_d_0_s = 0
    out % m___f2dace_SOA_bdy_mflx_e_idx_d_0_s = 0
    out % m___f2dace_SOA_coeff_gradp_d_0_s = 0
    out % m___f2dace_SOA_coeff_gradp_d_1_s = 0
    out % m___f2dace_SOA_coeff_gradp_d_2_s = 0
    out % m___f2dace_SOA_coeff_gradp_d_3_s = 0
    out % m___f2dace_SOA_d2dexdz2_fac1_mc_d_0_s = 0
    out % m___f2dace_SOA_d2dexdz2_fac1_mc_d_1_s = 0
    out % m___f2dace_SOA_d2dexdz2_fac1_mc_d_2_s = 0
    out % m___f2dace_SOA_d2dexdz2_fac2_mc_d_0_s = 0
    out % m___f2dace_SOA_d2dexdz2_fac2_mc_d_1_s = 0
    out % m___f2dace_SOA_d2dexdz2_fac2_mc_d_2_s = 0
    out % m___f2dace_SOA_d_exner_dz_ref_ic_d_0_s = 0
    out % m___f2dace_SOA_d_exner_dz_ref_ic_d_1_s = 0
    out % m___f2dace_SOA_d_exner_dz_ref_ic_d_2_s = 0
    out % m___f2dace_SOA_deepatmo_divh_mc_d_0_s = 0
    out % m___f2dace_SOA_deepatmo_divzl_mc_d_0_s = 0
    out % m___f2dace_SOA_deepatmo_divzu_mc_d_0_s = 0
    out % m___f2dace_SOA_exner_exfac_d_0_s = 0
    out % m___f2dace_SOA_exner_exfac_d_1_s = 0
    out % m___f2dace_SOA_exner_exfac_d_2_s = 0
    out % m___f2dace_SOA_exner_ref_mc_d_0_s = 0
    out % m___f2dace_SOA_exner_ref_mc_d_1_s = 0
    out % m___f2dace_SOA_exner_ref_mc_d_2_s = 0
    out % m___f2dace_SOA_hmask_dd3d_d_0_s = 0
    out % m___f2dace_SOA_hmask_dd3d_d_1_s = 0
    out % m___f2dace_SOA_inv_ddqz_z_full_d_0_s = 0
    out % m___f2dace_SOA_inv_ddqz_z_full_d_1_s = 0
    out % m___f2dace_SOA_inv_ddqz_z_full_d_2_s = 0
    out % m___f2dace_SOA_pg_edgeblk_d_0_s = 0
    out % m___f2dace_SOA_pg_edgeidx_d_0_s = 0
    out % m___f2dace_SOA_pg_exdist_d_0_s = 0
    out % m___f2dace_SOA_pg_vertidx_d_0_s = 0
    out % m___f2dace_SOA_rayleigh_vn_d_0_s = 0
    out % m___f2dace_SOA_rayleigh_w_d_0_s = 0
    out % m___f2dace_SOA_rho_ref_mc_d_0_s = 0
    out % m___f2dace_SOA_rho_ref_mc_d_1_s = 0
    out % m___f2dace_SOA_rho_ref_mc_d_2_s = 0
    out % m___f2dace_SOA_rho_ref_me_d_0_s = 0
    out % m___f2dace_SOA_rho_ref_me_d_1_s = 0
    out % m___f2dace_SOA_rho_ref_me_d_2_s = 0
    out % m___f2dace_SOA_scalfac_dd3d_d_0_s = 0
    out % m___f2dace_SOA_theta_ref_ic_d_0_s = 0
    out % m___f2dace_SOA_theta_ref_ic_d_1_s = 0
    out % m___f2dace_SOA_theta_ref_ic_d_2_s = 0
    out % m___f2dace_SOA_theta_ref_mc_d_0_s = 0
    out % m___f2dace_SOA_theta_ref_mc_d_1_s = 0
    out % m___f2dace_SOA_theta_ref_mc_d_2_s = 0
    out % m___f2dace_SOA_theta_ref_me_d_0_s = 0
    out % m___f2dace_SOA_theta_ref_me_d_1_s = 0
    out % m___f2dace_SOA_theta_ref_me_d_2_s = 0
    out % m___f2dace_SOA_vertidx_gradp_d_0_s = 0
    out % m___f2dace_SOA_vertidx_gradp_d_1_s = 0
    out % m___f2dace_SOA_vertidx_gradp_d_2_s = 0
    out % m___f2dace_SOA_vertidx_gradp_d_3_s = 0
    out % m___f2dace_SOA_vwind_expl_wgt_d_0_s = 0
    out % m___f2dace_SOA_vwind_expl_wgt_d_1_s = 0
    out % m___f2dace_SOA_vwind_impl_wgt_d_0_s = 0
    out % m___f2dace_SOA_vwind_impl_wgt_d_1_s = 0
    out % m___f2dace_SOA_wgtfacq1_c_d_0_s = 0
    out % m___f2dace_SOA_wgtfacq1_c_d_1_s = 0
    out % m___f2dace_SOA_wgtfacq1_c_d_2_s = 0
    out % m___f2dace_SOA_wgtfacq_c_d_0_s = 0
    out % m___f2dace_SOA_wgtfacq_c_d_1_s = 0
    out % m___f2dace_SOA_wgtfacq_c_d_2_s = 0
    out % m___f2dace_SOA_zdiff_gradp_d_0_s = 0
    out % m___f2dace_SOA_zdiff_gradp_d_1_s = 0
    out % m___f2dace_SOA_zdiff_gradp_d_2_s = 0
    out % m___f2dace_SOA_zdiff_gradp_d_3_s = 0
    out % m_bdy_mflx_e_blk = c_null_ptr
    out % m_bdy_mflx_e_dim = 0
    out % m_bdy_mflx_e_idx = c_null_ptr
    out % m_coeff_gradp = c_null_ptr
    out % m_d2dexdz2_fac1_mc = c_null_ptr
    out % m_d2dexdz2_fac1_mc_d_1_s = c_null_ptr
    out % m_d2dexdz2_fac1_mc_d_2_s = c_null_ptr
    out % m_d2dexdz2_fac2_mc = c_null_ptr
    out % m_d2dexdz2_fac2_mc_d_1_s = c_null_ptr
    out % m_d2dexdz2_fac2_mc_d_2_s = c_null_ptr
    out % m_d_exner_dz_ref_ic = c_null_ptr
    out % m_deepatmo_divh_mc = c_null_ptr
    out % m_deepatmo_divzl_mc = c_null_ptr
    out % m_deepatmo_divzu_mc = c_null_ptr
    out % m_exner_exfac = c_null_ptr
    out % m_exner_ref_mc = c_null_ptr
    out % m_hmask_dd3d = c_null_ptr
    out % m_inv_ddqz_z_full = c_null_ptr
    out % m_pg_edgeblk = c_null_ptr
    out % m_pg_edgeidx = c_null_ptr
    out % m_pg_exdist = c_null_ptr
    out % m_pg_listdim = 0
    out % m_pg_vertidx = c_null_ptr
    out % m_rayleigh_vn = c_null_ptr
    out % m_rayleigh_w = c_null_ptr
    out % m_rho_ref_mc = c_null_ptr
    out % m_rho_ref_me = c_null_ptr
    out % m_scalfac_dd3d = c_null_ptr
    out % m_theta_ref_ic = c_null_ptr
    out % m_theta_ref_mc = c_null_ptr
    out % m_theta_ref_me = c_null_ptr
    out % m_vertidx_gradp = c_null_ptr
    out % m_vwind_expl_wgt = c_null_ptr
    out % m_vwind_impl_wgt = c_null_ptr
    out % m_wgtfacq1_c = c_null_ptr
    out % m_wgtfacq_c = c_null_ptr
    out % m_zdiff_gradp = c_null_ptr
    IF (initalloc .AND. .NOT. ALLOCATED(a_coeff1_dwdz)) ALLOCATE(a_coeff1_dwdz(SIZE(inp % coeff1_dwdz, 1), SIZE(inp % coeff1_dwdz, 2), SIZE(inp % coeff1_dwdz, 3)))
    a_coeff1_dwdz = inp % coeff1_dwdz
    out % m_coeff1_dwdz = c_loc(a_coeff1_dwdz)
    out % m___f2dace_SA_coeff1_dwdz_d_0_s_332 = SIZE(inp % coeff1_dwdz, 1)
    out % m___f2dace_SA_coeff1_dwdz_d_1_s_333 = SIZE(inp % coeff1_dwdz, 2)
    out % m___f2dace_SA_coeff1_dwdz_d_2_s_334 = SIZE(inp % coeff1_dwdz, 3)
    out % m___f2dace_SOA_coeff1_dwdz_d_0_s_332 = LBOUND(inp % coeff1_dwdz, 1)
    out % m___f2dace_SOA_coeff1_dwdz_d_1_s_333 = LBOUND(inp % coeff1_dwdz, 2)
    out % m___f2dace_SOA_coeff1_dwdz_d_2_s_334 = LBOUND(inp % coeff1_dwdz, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_coeff2_dwdz)) ALLOCATE(a_coeff2_dwdz(SIZE(inp % coeff2_dwdz, 1), SIZE(inp % coeff2_dwdz, 2), SIZE(inp % coeff2_dwdz, 3)))
    a_coeff2_dwdz = inp % coeff2_dwdz
    out % m_coeff2_dwdz = c_loc(a_coeff2_dwdz)
    out % m___f2dace_SA_coeff2_dwdz_d_0_s_335 = SIZE(inp % coeff2_dwdz, 1)
    out % m___f2dace_SA_coeff2_dwdz_d_1_s_336 = SIZE(inp % coeff2_dwdz, 2)
    out % m___f2dace_SA_coeff2_dwdz_d_2_s_337 = SIZE(inp % coeff2_dwdz, 3)
    out % m___f2dace_SOA_coeff2_dwdz_d_0_s_335 = LBOUND(inp % coeff2_dwdz, 1)
    out % m___f2dace_SOA_coeff2_dwdz_d_1_s_336 = LBOUND(inp % coeff2_dwdz, 2)
    out % m___f2dace_SOA_coeff2_dwdz_d_2_s_337 = LBOUND(inp % coeff2_dwdz, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_coeff_gradekin)) ALLOCATE(a_coeff_gradekin(SIZE(inp % coeff_gradekin, 1), SIZE(inp % coeff_gradekin, 2), SIZE(inp % coeff_gradekin, 3)))
    a_coeff_gradekin = inp % coeff_gradekin
    out % m_coeff_gradekin = c_loc(a_coeff_gradekin)
    out % m___f2dace_SA_coeff_gradekin_d_0_s_329 = SIZE(inp % coeff_gradekin, 1)
    out % m___f2dace_SA_coeff_gradekin_d_1_s_330 = SIZE(inp % coeff_gradekin, 2)
    out % m___f2dace_SA_coeff_gradekin_d_2_s_331 = SIZE(inp % coeff_gradekin, 3)
    out % m___f2dace_SOA_coeff_gradekin_d_0_s_329 = LBOUND(inp % coeff_gradekin, 1)
    out % m___f2dace_SOA_coeff_gradekin_d_1_s_330 = LBOUND(inp % coeff_gradekin, 2)
    out % m___f2dace_SOA_coeff_gradekin_d_2_s_331 = LBOUND(inp % coeff_gradekin, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_ddqz_z_full_e)) ALLOCATE(a_ddqz_z_full_e(SIZE(inp % ddqz_z_full_e, 1), SIZE(inp % ddqz_z_full_e, 2), SIZE(inp % ddqz_z_full_e, 3)))
    a_ddqz_z_full_e = inp % ddqz_z_full_e
    out % m_ddqz_z_full_e = c_loc(a_ddqz_z_full_e)
    out % m___f2dace_SA_ddqz_z_full_e_d_0_s_314 = SIZE(inp % ddqz_z_full_e, 1)
    out % m___f2dace_SA_ddqz_z_full_e_d_1_s_315 = SIZE(inp % ddqz_z_full_e, 2)
    out % m___f2dace_SA_ddqz_z_full_e_d_2_s_316 = SIZE(inp % ddqz_z_full_e, 3)
    out % m___f2dace_SOA_ddqz_z_full_e_d_0_s_314 = LBOUND(inp % ddqz_z_full_e, 1)
    out % m___f2dace_SOA_ddqz_z_full_e_d_1_s_315 = LBOUND(inp % ddqz_z_full_e, 2)
    out % m___f2dace_SOA_ddqz_z_full_e_d_2_s_316 = LBOUND(inp % ddqz_z_full_e, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_ddqz_z_half)) ALLOCATE(a_ddqz_z_half(SIZE(inp % ddqz_z_half, 1), SIZE(inp % ddqz_z_half, 2), SIZE(inp % ddqz_z_half, 3)))
    a_ddqz_z_half = inp % ddqz_z_half
    out % m_ddqz_z_half = c_loc(a_ddqz_z_half)
    out % m___f2dace_SA_ddqz_z_half_d_0_s_317 = SIZE(inp % ddqz_z_half, 1)
    out % m___f2dace_SA_ddqz_z_half_d_1_s_318 = SIZE(inp % ddqz_z_half, 2)
    out % m___f2dace_SA_ddqz_z_half_d_2_s_319 = SIZE(inp % ddqz_z_half, 3)
    out % m___f2dace_SOA_ddqz_z_half_d_0_s_317 = LBOUND(inp % ddqz_z_half, 1)
    out % m___f2dace_SOA_ddqz_z_half_d_1_s_318 = LBOUND(inp % ddqz_z_half, 2)
    out % m___f2dace_SOA_ddqz_z_half_d_2_s_319 = LBOUND(inp % ddqz_z_half, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_ddxn_z_full)) ALLOCATE(a_ddxn_z_full(SIZE(inp % ddxn_z_full, 1), SIZE(inp % ddxn_z_full, 2), SIZE(inp % ddxn_z_full, 3)))
    a_ddxn_z_full = inp % ddxn_z_full
    out % m_ddxn_z_full = c_loc(a_ddxn_z_full)
    out % m___f2dace_SA_ddxn_z_full_d_0_s_308 = SIZE(inp % ddxn_z_full, 1)
    out % m___f2dace_SA_ddxn_z_full_d_1_s_309 = SIZE(inp % ddxn_z_full, 2)
    out % m___f2dace_SA_ddxn_z_full_d_2_s_310 = SIZE(inp % ddxn_z_full, 3)
    out % m___f2dace_SOA_ddxn_z_full_d_0_s_308 = LBOUND(inp % ddxn_z_full, 1)
    out % m___f2dace_SOA_ddxn_z_full_d_1_s_309 = LBOUND(inp % ddxn_z_full, 2)
    out % m___f2dace_SOA_ddxn_z_full_d_2_s_310 = LBOUND(inp % ddxn_z_full, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_ddxt_z_full)) ALLOCATE(a_ddxt_z_full(SIZE(inp % ddxt_z_full, 1), SIZE(inp % ddxt_z_full, 2), SIZE(inp % ddxt_z_full, 3)))
    a_ddxt_z_full = inp % ddxt_z_full
    out % m_ddxt_z_full = c_loc(a_ddxt_z_full)
    out % m___f2dace_SA_ddxt_z_full_d_0_s_311 = SIZE(inp % ddxt_z_full, 1)
    out % m___f2dace_SA_ddxt_z_full_d_1_s_312 = SIZE(inp % ddxt_z_full, 2)
    out % m___f2dace_SA_ddxt_z_full_d_2_s_313 = SIZE(inp % ddxt_z_full, 3)
    out % m___f2dace_SOA_ddxt_z_full_d_0_s_311 = LBOUND(inp % ddxt_z_full, 1)
    out % m___f2dace_SOA_ddxt_z_full_d_1_s_312 = LBOUND(inp % ddxt_z_full, 2)
    out % m___f2dace_SOA_ddxt_z_full_d_2_s_313 = LBOUND(inp % ddxt_z_full, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_deepatmo_gradh_ifc)) ALLOCATE(a_deepatmo_gradh_ifc(SIZE(inp % deepatmo_gradh_ifc, 1)))
    a_deepatmo_gradh_ifc = inp % deepatmo_gradh_ifc
    out % m_deepatmo_gradh_ifc = c_loc(a_deepatmo_gradh_ifc)
    out % m___f2dace_SA_deepatmo_gradh_ifc_d_0_s_340 = SIZE(inp % deepatmo_gradh_ifc, 1)
    out % m___f2dace_SOA_deepatmo_gradh_ifc_d_0_s_340 = LBOUND(inp % deepatmo_gradh_ifc, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_deepatmo_gradh_mc)) ALLOCATE(a_deepatmo_gradh_mc(SIZE(inp % deepatmo_gradh_mc, 1)))
    a_deepatmo_gradh_mc = inp % deepatmo_gradh_mc
    out % m_deepatmo_gradh_mc = c_loc(a_deepatmo_gradh_mc)
    out % m___f2dace_SA_deepatmo_gradh_mc_d_0_s_338 = SIZE(inp % deepatmo_gradh_mc, 1)
    out % m___f2dace_SOA_deepatmo_gradh_mc_d_0_s_338 = LBOUND(inp % deepatmo_gradh_mc, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_deepatmo_invr_ifc)) ALLOCATE(a_deepatmo_invr_ifc(SIZE(inp % deepatmo_invr_ifc, 1)))
    a_deepatmo_invr_ifc = inp % deepatmo_invr_ifc
    out % m_deepatmo_invr_ifc = c_loc(a_deepatmo_invr_ifc)
    out % m___f2dace_SA_deepatmo_invr_ifc_d_0_s_341 = SIZE(inp % deepatmo_invr_ifc, 1)
    out % m___f2dace_SOA_deepatmo_invr_ifc_d_0_s_341 = LBOUND(inp % deepatmo_invr_ifc, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_deepatmo_invr_mc)) ALLOCATE(a_deepatmo_invr_mc(SIZE(inp % deepatmo_invr_mc, 1)))
    a_deepatmo_invr_mc = inp % deepatmo_invr_mc
    out % m_deepatmo_invr_mc = c_loc(a_deepatmo_invr_mc)
    out % m___f2dace_SA_deepatmo_invr_mc_d_0_s_339 = SIZE(inp % deepatmo_invr_mc, 1)
    out % m___f2dace_SOA_deepatmo_invr_mc_d_0_s_339 = LBOUND(inp % deepatmo_invr_mc, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_wgtfac_c)) ALLOCATE(a_wgtfac_c(SIZE(inp % wgtfac_c, 1), SIZE(inp % wgtfac_c, 2), SIZE(inp % wgtfac_c, 3)))
    a_wgtfac_c = inp % wgtfac_c
    out % m_wgtfac_c = c_loc(a_wgtfac_c)
    out % m___f2dace_SA_wgtfac_c_d_0_s_320 = SIZE(inp % wgtfac_c, 1)
    out % m___f2dace_SA_wgtfac_c_d_1_s_321 = SIZE(inp % wgtfac_c, 2)
    out % m___f2dace_SA_wgtfac_c_d_2_s_322 = SIZE(inp % wgtfac_c, 3)
    out % m___f2dace_SOA_wgtfac_c_d_0_s_320 = LBOUND(inp % wgtfac_c, 1)
    out % m___f2dace_SOA_wgtfac_c_d_1_s_321 = LBOUND(inp % wgtfac_c, 2)
    out % m___f2dace_SOA_wgtfac_c_d_2_s_322 = LBOUND(inp % wgtfac_c, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_wgtfac_e)) ALLOCATE(a_wgtfac_e(SIZE(inp % wgtfac_e, 1), SIZE(inp % wgtfac_e, 2), SIZE(inp % wgtfac_e, 3)))
    a_wgtfac_e = inp % wgtfac_e
    out % m_wgtfac_e = c_loc(a_wgtfac_e)
    out % m___f2dace_SA_wgtfac_e_d_0_s_323 = SIZE(inp % wgtfac_e, 1)
    out % m___f2dace_SA_wgtfac_e_d_1_s_324 = SIZE(inp % wgtfac_e, 2)
    out % m___f2dace_SA_wgtfac_e_d_2_s_325 = SIZE(inp % wgtfac_e, 3)
    out % m___f2dace_SOA_wgtfac_e_d_0_s_323 = LBOUND(inp % wgtfac_e, 1)
    out % m___f2dace_SOA_wgtfac_e_d_1_s_324 = LBOUND(inp % wgtfac_e, 2)
    out % m___f2dace_SOA_wgtfac_e_d_2_s_325 = LBOUND(inp % wgtfac_e, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_wgtfacq_e)) ALLOCATE(a_wgtfacq_e(SIZE(inp % wgtfacq_e, 1), SIZE(inp % wgtfacq_e, 2), SIZE(inp % wgtfacq_e, 3)))
    a_wgtfacq_e = inp % wgtfacq_e
    out % m_wgtfacq_e = c_loc(a_wgtfacq_e)
    out % m___f2dace_SA_wgtfacq_e_d_0_s_326 = SIZE(inp % wgtfacq_e, 1)
    out % m___f2dace_SA_wgtfacq_e_d_1_s_327 = SIZE(inp % wgtfacq_e, 2)
    out % m___f2dace_SA_wgtfacq_e_d_2_s_328 = SIZE(inp % wgtfacq_e, 3)
    out % m___f2dace_SOA_wgtfacq_e_d_0_s_326 = LBOUND(inp % wgtfacq_e, 1)
    out % m___f2dace_SOA_wgtfacq_e_d_1_s_327 = LBOUND(inp % wgtfacq_e, 2)
    out % m___f2dace_SOA_wgtfacq_e_d_2_s_328 = LBOUND(inp % wgtfacq_e, 3)
  END SUBROUTINE ctor_t_nh_metrics
  SUBROUTINE ctor_t_nh_diag(inp, out, initalloc)
    TYPE(t_nh_diag), INTENT(IN) :: inp
    TYPE(glue_t_nh_diag), INTENT(INOUT) :: out
    LOGICAL, INTENT(IN) :: initalloc
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_ddt_vn_apc_pc(:, :, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_ddt_w_adv_pc(:, :, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_vn_ie(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_vt(:, :, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_w_concorr_c(:, :, :)
    out % m___f2dace_SA_ddt_exner_phy_d_0_s = 0
    out % m___f2dace_SA_ddt_exner_phy_d_1_s = 0
    out % m___f2dace_SA_ddt_exner_phy_d_2_s = 0
    out % m___f2dace_SA_ddt_vn_adv_d_0_s = 0
    out % m___f2dace_SA_ddt_vn_adv_d_1_s = 0
    out % m___f2dace_SA_ddt_vn_adv_d_2_s = 0
    out % m___f2dace_SA_ddt_vn_cor_d_0_s = 0
    out % m___f2dace_SA_ddt_vn_cor_d_1_s = 0
    out % m___f2dace_SA_ddt_vn_cor_d_2_s = 0
    out % m___f2dace_SA_ddt_vn_cor_pc_d_0_s = 0
    out % m___f2dace_SA_ddt_vn_cor_pc_d_1_s = 0
    out % m___f2dace_SA_ddt_vn_cor_pc_d_2_s = 0
    out % m___f2dace_SA_ddt_vn_cor_pc_d_3_s = 0
    out % m___f2dace_SA_ddt_vn_dmp_d_0_s = 0
    out % m___f2dace_SA_ddt_vn_dmp_d_1_s = 0
    out % m___f2dace_SA_ddt_vn_dmp_d_2_s = 0
    out % m___f2dace_SA_ddt_vn_dyn_d_0_s = 0
    out % m___f2dace_SA_ddt_vn_dyn_d_1_s = 0
    out % m___f2dace_SA_ddt_vn_dyn_d_2_s = 0
    out % m___f2dace_SA_ddt_vn_grf_d_0_s = 0
    out % m___f2dace_SA_ddt_vn_grf_d_1_s = 0
    out % m___f2dace_SA_ddt_vn_grf_d_2_s = 0
    out % m___f2dace_SA_ddt_vn_iau_d_0_s = 0
    out % m___f2dace_SA_ddt_vn_iau_d_1_s = 0
    out % m___f2dace_SA_ddt_vn_iau_d_2_s = 0
    out % m___f2dace_SA_ddt_vn_pgr_d_0_s = 0
    out % m___f2dace_SA_ddt_vn_pgr_d_1_s = 0
    out % m___f2dace_SA_ddt_vn_pgr_d_2_s = 0
    out % m___f2dace_SA_ddt_vn_phd_d_0_s = 0
    out % m___f2dace_SA_ddt_vn_phd_d_1_s = 0
    out % m___f2dace_SA_ddt_vn_phd_d_2_s = 0
    out % m___f2dace_SA_ddt_vn_phy_d_0_s = 0
    out % m___f2dace_SA_ddt_vn_phy_d_1_s = 0
    out % m___f2dace_SA_ddt_vn_phy_d_2_s = 0
    out % m___f2dace_SA_ddt_vn_ray_d_0_s = 0
    out % m___f2dace_SA_ddt_vn_ray_d_1_s = 0
    out % m___f2dace_SA_ddt_vn_ray_d_2_s = 0
    out % m___f2dace_SA_exner_dyn_incr_d_0_s = 0
    out % m___f2dace_SA_exner_dyn_incr_d_1_s = 0
    out % m___f2dace_SA_exner_dyn_incr_d_2_s = 0
    out % m___f2dace_SA_exner_incr_d_0_s = 0
    out % m___f2dace_SA_exner_incr_d_1_s = 0
    out % m___f2dace_SA_exner_incr_d_2_s = 0
    out % m___f2dace_SA_exner_pr_d_0_s = 0
    out % m___f2dace_SA_exner_pr_d_1_s = 0
    out % m___f2dace_SA_exner_pr_d_2_s = 0
    out % m___f2dace_SA_grf_bdy_mflx_d_0_s = 0
    out % m___f2dace_SA_grf_bdy_mflx_d_1_s = 0
    out % m___f2dace_SA_grf_bdy_mflx_d_2_s = 0
    out % m___f2dace_SA_grf_tend_mflx_d_0_s = 0
    out % m___f2dace_SA_grf_tend_mflx_d_1_s = 0
    out % m___f2dace_SA_grf_tend_mflx_d_2_s = 0
    out % m___f2dace_SA_grf_tend_rho_d_0_s = 0
    out % m___f2dace_SA_grf_tend_rho_d_1_s = 0
    out % m___f2dace_SA_grf_tend_rho_d_2_s = 0
    out % m___f2dace_SA_grf_tend_thv_d_0_s = 0
    out % m___f2dace_SA_grf_tend_thv_d_1_s = 0
    out % m___f2dace_SA_grf_tend_thv_d_2_s = 0
    out % m___f2dace_SA_grf_tend_vn_d_0_s = 0
    out % m___f2dace_SA_grf_tend_vn_d_1_s = 0
    out % m___f2dace_SA_grf_tend_vn_d_2_s = 0
    out % m___f2dace_SA_grf_tend_w_d_0_s = 0
    out % m___f2dace_SA_grf_tend_w_d_1_s = 0
    out % m___f2dace_SA_grf_tend_w_d_2_s = 0
    out % m___f2dace_SA_mass_fl_e_d_0_s = 0
    out % m___f2dace_SA_mass_fl_e_d_1_s = 0
    out % m___f2dace_SA_mass_fl_e_d_2_s = 0
    out % m___f2dace_SA_mass_fl_e_sv_d_0_s = 0
    out % m___f2dace_SA_mass_fl_e_sv_d_1_s = 0
    out % m___f2dace_SA_mass_fl_e_sv_d_2_s = 0
    out % m___f2dace_SA_mflx_ic_int_d_0_s = 0
    out % m___f2dace_SA_mflx_ic_int_d_1_s = 0
    out % m___f2dace_SA_mflx_ic_int_d_2_s = 0
    out % m___f2dace_SA_mflx_ic_ubc_d_0_s = 0
    out % m___f2dace_SA_mflx_ic_ubc_d_1_s = 0
    out % m___f2dace_SA_mflx_ic_ubc_d_2_s = 0
    out % m___f2dace_SA_rho_ic_d_0_s = 0
    out % m___f2dace_SA_rho_ic_d_1_s = 0
    out % m___f2dace_SA_rho_ic_d_2_s = 0
    out % m___f2dace_SA_rho_ic_int_d_0_s = 0
    out % m___f2dace_SA_rho_ic_int_d_1_s = 0
    out % m___f2dace_SA_rho_ic_int_d_2_s = 0
    out % m___f2dace_SA_rho_ic_ubc_d_0_s = 0
    out % m___f2dace_SA_rho_ic_ubc_d_1_s = 0
    out % m___f2dace_SA_rho_ic_ubc_d_2_s = 0
    out % m___f2dace_SA_rho_incr_d_0_s = 0
    out % m___f2dace_SA_rho_incr_d_1_s = 0
    out % m___f2dace_SA_rho_incr_d_2_s = 0
    out % m___f2dace_SA_theta_v_ic_d_0_s = 0
    out % m___f2dace_SA_theta_v_ic_d_1_s = 0
    out % m___f2dace_SA_theta_v_ic_d_2_s = 0
    out % m___f2dace_SA_theta_v_ic_int_d_0_s = 0
    out % m___f2dace_SA_theta_v_ic_int_d_1_s = 0
    out % m___f2dace_SA_theta_v_ic_int_d_2_s = 0
    out % m___f2dace_SA_theta_v_ic_ubc_d_0_s = 0
    out % m___f2dace_SA_theta_v_ic_ubc_d_1_s = 0
    out % m___f2dace_SA_theta_v_ic_ubc_d_2_s = 0
    out % m___f2dace_SA_vn_ie_int_d_0_s = 0
    out % m___f2dace_SA_vn_ie_int_d_1_s = 0
    out % m___f2dace_SA_vn_ie_int_d_2_s = 0
    out % m___f2dace_SA_vn_ie_ubc_d_0_s = 0
    out % m___f2dace_SA_vn_ie_ubc_d_1_s = 0
    out % m___f2dace_SA_vn_ie_ubc_d_2_s = 0
    out % m___f2dace_SA_vn_incr_d_0_s = 0
    out % m___f2dace_SA_vn_incr_d_1_s = 0
    out % m___f2dace_SA_vn_incr_d_2_s = 0
    out % m___f2dace_SA_w_int_d_0_s = 0
    out % m___f2dace_SA_w_int_d_1_s = 0
    out % m___f2dace_SA_w_int_d_2_s = 0
    out % m___f2dace_SA_w_ubc_d_0_s = 0
    out % m___f2dace_SA_w_ubc_d_1_s = 0
    out % m___f2dace_SA_w_ubc_d_2_s = 0
    out % m___f2dace_SOA_ddt_exner_phy_d_0_s = 0
    out % m___f2dace_SOA_ddt_exner_phy_d_1_s = 0
    out % m___f2dace_SOA_ddt_exner_phy_d_2_s = 0
    out % m___f2dace_SOA_ddt_vn_adv_d_0_s = 0
    out % m___f2dace_SOA_ddt_vn_adv_d_1_s = 0
    out % m___f2dace_SOA_ddt_vn_adv_d_2_s = 0
    out % m___f2dace_SOA_ddt_vn_cor_d_0_s = 0
    out % m___f2dace_SOA_ddt_vn_cor_d_1_s = 0
    out % m___f2dace_SOA_ddt_vn_cor_d_2_s = 0
    out % m___f2dace_SOA_ddt_vn_cor_pc_d_0_s = 0
    out % m___f2dace_SOA_ddt_vn_cor_pc_d_1_s = 0
    out % m___f2dace_SOA_ddt_vn_cor_pc_d_2_s = 0
    out % m___f2dace_SOA_ddt_vn_cor_pc_d_3_s = 0
    out % m___f2dace_SOA_ddt_vn_dmp_d_0_s = 0
    out % m___f2dace_SOA_ddt_vn_dmp_d_1_s = 0
    out % m___f2dace_SOA_ddt_vn_dmp_d_2_s = 0
    out % m___f2dace_SOA_ddt_vn_dyn_d_0_s = 0
    out % m___f2dace_SOA_ddt_vn_dyn_d_1_s = 0
    out % m___f2dace_SOA_ddt_vn_dyn_d_2_s = 0
    out % m___f2dace_SOA_ddt_vn_grf_d_0_s = 0
    out % m___f2dace_SOA_ddt_vn_grf_d_1_s = 0
    out % m___f2dace_SOA_ddt_vn_grf_d_2_s = 0
    out % m___f2dace_SOA_ddt_vn_iau_d_0_s = 0
    out % m___f2dace_SOA_ddt_vn_iau_d_1_s = 0
    out % m___f2dace_SOA_ddt_vn_iau_d_2_s = 0
    out % m___f2dace_SOA_ddt_vn_pgr_d_0_s = 0
    out % m___f2dace_SOA_ddt_vn_pgr_d_1_s = 0
    out % m___f2dace_SOA_ddt_vn_pgr_d_2_s = 0
    out % m___f2dace_SOA_ddt_vn_phd_d_0_s = 0
    out % m___f2dace_SOA_ddt_vn_phd_d_1_s = 0
    out % m___f2dace_SOA_ddt_vn_phd_d_2_s = 0
    out % m___f2dace_SOA_ddt_vn_phy_d_0_s = 0
    out % m___f2dace_SOA_ddt_vn_phy_d_1_s = 0
    out % m___f2dace_SOA_ddt_vn_phy_d_2_s = 0
    out % m___f2dace_SOA_ddt_vn_ray_d_0_s = 0
    out % m___f2dace_SOA_ddt_vn_ray_d_1_s = 0
    out % m___f2dace_SOA_ddt_vn_ray_d_2_s = 0
    out % m___f2dace_SOA_exner_dyn_incr_d_0_s = 0
    out % m___f2dace_SOA_exner_dyn_incr_d_1_s = 0
    out % m___f2dace_SOA_exner_dyn_incr_d_2_s = 0
    out % m___f2dace_SOA_exner_incr_d_0_s = 0
    out % m___f2dace_SOA_exner_incr_d_1_s = 0
    out % m___f2dace_SOA_exner_incr_d_2_s = 0
    out % m___f2dace_SOA_exner_pr_d_0_s = 0
    out % m___f2dace_SOA_exner_pr_d_1_s = 0
    out % m___f2dace_SOA_exner_pr_d_2_s = 0
    out % m___f2dace_SOA_grf_bdy_mflx_d_0_s = 0
    out % m___f2dace_SOA_grf_bdy_mflx_d_1_s = 0
    out % m___f2dace_SOA_grf_bdy_mflx_d_2_s = 0
    out % m___f2dace_SOA_grf_tend_mflx_d_0_s = 0
    out % m___f2dace_SOA_grf_tend_mflx_d_1_s = 0
    out % m___f2dace_SOA_grf_tend_mflx_d_2_s = 0
    out % m___f2dace_SOA_grf_tend_rho_d_0_s = 0
    out % m___f2dace_SOA_grf_tend_rho_d_1_s = 0
    out % m___f2dace_SOA_grf_tend_rho_d_2_s = 0
    out % m___f2dace_SOA_grf_tend_thv_d_0_s = 0
    out % m___f2dace_SOA_grf_tend_thv_d_1_s = 0
    out % m___f2dace_SOA_grf_tend_thv_d_2_s = 0
    out % m___f2dace_SOA_grf_tend_vn_d_0_s = 0
    out % m___f2dace_SOA_grf_tend_vn_d_1_s = 0
    out % m___f2dace_SOA_grf_tend_vn_d_2_s = 0
    out % m___f2dace_SOA_grf_tend_w_d_0_s = 0
    out % m___f2dace_SOA_grf_tend_w_d_1_s = 0
    out % m___f2dace_SOA_grf_tend_w_d_2_s = 0
    out % m___f2dace_SOA_mass_fl_e_d_0_s = 0
    out % m___f2dace_SOA_mass_fl_e_d_1_s = 0
    out % m___f2dace_SOA_mass_fl_e_d_2_s = 0
    out % m___f2dace_SOA_mass_fl_e_sv_d_0_s = 0
    out % m___f2dace_SOA_mass_fl_e_sv_d_1_s = 0
    out % m___f2dace_SOA_mass_fl_e_sv_d_2_s = 0
    out % m___f2dace_SOA_mflx_ic_int_d_0_s = 0
    out % m___f2dace_SOA_mflx_ic_int_d_1_s = 0
    out % m___f2dace_SOA_mflx_ic_int_d_2_s = 0
    out % m___f2dace_SOA_mflx_ic_ubc_d_0_s = 0
    out % m___f2dace_SOA_mflx_ic_ubc_d_1_s = 0
    out % m___f2dace_SOA_mflx_ic_ubc_d_2_s = 0
    out % m___f2dace_SOA_rho_ic_d_0_s = 0
    out % m___f2dace_SOA_rho_ic_d_1_s = 0
    out % m___f2dace_SOA_rho_ic_d_2_s = 0
    out % m___f2dace_SOA_rho_ic_int_d_0_s = 0
    out % m___f2dace_SOA_rho_ic_int_d_1_s = 0
    out % m___f2dace_SOA_rho_ic_int_d_2_s = 0
    out % m___f2dace_SOA_rho_ic_ubc_d_0_s = 0
    out % m___f2dace_SOA_rho_ic_ubc_d_1_s = 0
    out % m___f2dace_SOA_rho_ic_ubc_d_2_s = 0
    out % m___f2dace_SOA_rho_incr_d_0_s = 0
    out % m___f2dace_SOA_rho_incr_d_1_s = 0
    out % m___f2dace_SOA_rho_incr_d_2_s = 0
    out % m___f2dace_SOA_theta_v_ic_d_0_s = 0
    out % m___f2dace_SOA_theta_v_ic_d_1_s = 0
    out % m___f2dace_SOA_theta_v_ic_d_2_s = 0
    out % m___f2dace_SOA_theta_v_ic_int_d_0_s = 0
    out % m___f2dace_SOA_theta_v_ic_int_d_1_s = 0
    out % m___f2dace_SOA_theta_v_ic_int_d_2_s = 0
    out % m___f2dace_SOA_theta_v_ic_ubc_d_0_s = 0
    out % m___f2dace_SOA_theta_v_ic_ubc_d_1_s = 0
    out % m___f2dace_SOA_theta_v_ic_ubc_d_2_s = 0
    out % m___f2dace_SOA_vn_ie_int_d_0_s = 0
    out % m___f2dace_SOA_vn_ie_int_d_1_s = 0
    out % m___f2dace_SOA_vn_ie_int_d_2_s = 0
    out % m___f2dace_SOA_vn_ie_ubc_d_0_s = 0
    out % m___f2dace_SOA_vn_ie_ubc_d_1_s = 0
    out % m___f2dace_SOA_vn_ie_ubc_d_2_s = 0
    out % m___f2dace_SOA_vn_incr_d_0_s = 0
    out % m___f2dace_SOA_vn_incr_d_1_s = 0
    out % m___f2dace_SOA_vn_incr_d_2_s = 0
    out % m___f2dace_SOA_w_int_d_0_s = 0
    out % m___f2dace_SOA_w_int_d_1_s = 0
    out % m___f2dace_SOA_w_int_d_2_s = 0
    out % m___f2dace_SOA_w_ubc_d_0_s = 0
    out % m___f2dace_SOA_w_ubc_d_1_s = 0
    out % m___f2dace_SOA_w_ubc_d_2_s = 0
    out % m_ddt_exner_phy = c_null_ptr
    out % m_ddt_vn_adv = c_null_ptr
    out % m_ddt_vn_adv_is_associated = 0
    out % m_ddt_vn_cor = c_null_ptr
    out % m_ddt_vn_cor_is_associated = 0
    out % m_ddt_vn_cor_pc = c_null_ptr
    out % m_ddt_vn_dmp = c_null_ptr
    out % m_ddt_vn_dmp_is_associated = 0
    out % m_ddt_vn_dyn = c_null_ptr
    out % m_ddt_vn_dyn_is_associated = 0
    out % m_ddt_vn_grf = c_null_ptr
    out % m_ddt_vn_grf_is_associated = 0
    out % m_ddt_vn_iau = c_null_ptr
    out % m_ddt_vn_iau_is_associated = 0
    out % m_ddt_vn_pgr = c_null_ptr
    out % m_ddt_vn_pgr_is_associated = 0
    out % m_ddt_vn_phd = c_null_ptr
    out % m_ddt_vn_phd_is_associated = 0
    out % m_ddt_vn_phy = c_null_ptr
    out % m_ddt_vn_ray = c_null_ptr
    out % m_ddt_vn_ray_is_associated = 0
    out % m_exner_dyn_incr = c_null_ptr
    out % m_exner_incr = c_null_ptr
    out % m_exner_pr = c_null_ptr
    out % m_grf_bdy_mflx = c_null_ptr
    out % m_grf_tend_mflx = c_null_ptr
    out % m_grf_tend_rho = c_null_ptr
    out % m_grf_tend_thv = c_null_ptr
    out % m_grf_tend_vn = c_null_ptr
    out % m_grf_tend_w = c_null_ptr
    out % m_mass_fl_e = c_null_ptr
    out % m_mass_fl_e_sv = c_null_ptr
    out % m_mflx_ic_int = c_null_ptr
    out % m_mflx_ic_ubc = c_null_ptr
    out % m_rho_ic = c_null_ptr
    out % m_rho_ic_int = c_null_ptr
    out % m_rho_ic_ubc = c_null_ptr
    out % m_rho_incr = c_null_ptr
    out % m_theta_v_ic = c_null_ptr
    out % m_theta_v_ic_int = c_null_ptr
    out % m_theta_v_ic_ubc = c_null_ptr
    out % m_vn_ie_int = c_null_ptr
    out % m_vn_ie_ubc = c_null_ptr
    out % m_vn_incr = c_null_ptr
    out % m_w_int = c_null_ptr
    out % m_w_ubc = c_null_ptr
    IF (initalloc .AND. .NOT. ALLOCATED(a_ddt_vn_apc_pc)) ALLOCATE(a_ddt_vn_apc_pc(SIZE(inp % ddt_vn_apc_pc, 1), SIZE(inp % ddt_vn_apc_pc, 2), SIZE(inp % ddt_vn_apc_pc, 3), SIZE(inp % ddt_vn_apc_pc, 4)))
    a_ddt_vn_apc_pc = inp % ddt_vn_apc_pc
    out % m_ddt_vn_apc_pc = c_loc(a_ddt_vn_apc_pc)
    out % m___f2dace_SA_ddt_vn_apc_pc_d_0_s_300 = SIZE(inp % ddt_vn_apc_pc, 1)
    out % m___f2dace_SA_ddt_vn_apc_pc_d_1_s_301 = SIZE(inp % ddt_vn_apc_pc, 2)
    out % m___f2dace_SA_ddt_vn_apc_pc_d_2_s_302 = SIZE(inp % ddt_vn_apc_pc, 3)
    out % m___f2dace_SA_ddt_vn_apc_pc_d_3_s_303 = SIZE(inp % ddt_vn_apc_pc, 4)
    out % m___f2dace_SOA_ddt_vn_apc_pc_d_0_s_300 = LBOUND(inp % ddt_vn_apc_pc, 1)
    out % m___f2dace_SOA_ddt_vn_apc_pc_d_1_s_301 = LBOUND(inp % ddt_vn_apc_pc, 2)
    out % m___f2dace_SOA_ddt_vn_apc_pc_d_2_s_302 = LBOUND(inp % ddt_vn_apc_pc, 3)
    out % m___f2dace_SOA_ddt_vn_apc_pc_d_3_s_303 = LBOUND(inp % ddt_vn_apc_pc, 4)
    IF (initalloc .AND. .NOT. ALLOCATED(a_ddt_w_adv_pc)) ALLOCATE(a_ddt_w_adv_pc(SIZE(inp % ddt_w_adv_pc, 1), SIZE(inp % ddt_w_adv_pc, 2), SIZE(inp % ddt_w_adv_pc, 3), SIZE(inp % ddt_w_adv_pc, 4)))
    a_ddt_w_adv_pc = inp % ddt_w_adv_pc
    out % m_ddt_w_adv_pc = c_loc(a_ddt_w_adv_pc)
    out % m___f2dace_SA_ddt_w_adv_pc_d_0_s_304 = SIZE(inp % ddt_w_adv_pc, 1)
    out % m___f2dace_SA_ddt_w_adv_pc_d_1_s_305 = SIZE(inp % ddt_w_adv_pc, 2)
    out % m___f2dace_SA_ddt_w_adv_pc_d_2_s_306 = SIZE(inp % ddt_w_adv_pc, 3)
    out % m___f2dace_SA_ddt_w_adv_pc_d_3_s_307 = SIZE(inp % ddt_w_adv_pc, 4)
    out % m___f2dace_SOA_ddt_w_adv_pc_d_0_s_304 = LBOUND(inp % ddt_w_adv_pc, 1)
    out % m___f2dace_SOA_ddt_w_adv_pc_d_1_s_305 = LBOUND(inp % ddt_w_adv_pc, 2)
    out % m___f2dace_SOA_ddt_w_adv_pc_d_2_s_306 = LBOUND(inp % ddt_w_adv_pc, 3)
    out % m___f2dace_SOA_ddt_w_adv_pc_d_3_s_307 = LBOUND(inp % ddt_w_adv_pc, 4)
    out % m_max_vcfl_dyn = inp % max_vcfl_dyn
    IF (initalloc .AND. .NOT. ALLOCATED(a_vn_ie)) ALLOCATE(a_vn_ie(SIZE(inp % vn_ie, 1), SIZE(inp % vn_ie, 2), SIZE(inp % vn_ie, 3)))
    a_vn_ie = inp % vn_ie
    out % m_vn_ie = c_loc(a_vn_ie)
    out % m___f2dace_SA_vn_ie_d_0_s_294 = SIZE(inp % vn_ie, 1)
    out % m___f2dace_SA_vn_ie_d_1_s_295 = SIZE(inp % vn_ie, 2)
    out % m___f2dace_SA_vn_ie_d_2_s_296 = SIZE(inp % vn_ie, 3)
    out % m___f2dace_SOA_vn_ie_d_0_s_294 = LBOUND(inp % vn_ie, 1)
    out % m___f2dace_SOA_vn_ie_d_1_s_295 = LBOUND(inp % vn_ie, 2)
    out % m___f2dace_SOA_vn_ie_d_2_s_296 = LBOUND(inp % vn_ie, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_vt)) ALLOCATE(a_vt(SIZE(inp % vt, 1), SIZE(inp % vt, 2), SIZE(inp % vt, 3)))
    a_vt = inp % vt
    out % m_vt = c_loc(a_vt)
    out % m___f2dace_SA_vt_d_0_s_291 = SIZE(inp % vt, 1)
    out % m___f2dace_SA_vt_d_1_s_292 = SIZE(inp % vt, 2)
    out % m___f2dace_SA_vt_d_2_s_293 = SIZE(inp % vt, 3)
    out % m___f2dace_SOA_vt_d_0_s_291 = LBOUND(inp % vt, 1)
    out % m___f2dace_SOA_vt_d_1_s_292 = LBOUND(inp % vt, 2)
    out % m___f2dace_SOA_vt_d_2_s_293 = LBOUND(inp % vt, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_w_concorr_c)) ALLOCATE(a_w_concorr_c(SIZE(inp % w_concorr_c, 1), SIZE(inp % w_concorr_c, 2), SIZE(inp % w_concorr_c, 3)))
    a_w_concorr_c = inp % w_concorr_c
    out % m_w_concorr_c = c_loc(a_w_concorr_c)
    out % m___f2dace_SA_w_concorr_c_d_0_s_297 = SIZE(inp % w_concorr_c, 1)
    out % m___f2dace_SA_w_concorr_c_d_1_s_298 = SIZE(inp % w_concorr_c, 2)
    out % m___f2dace_SA_w_concorr_c_d_2_s_299 = SIZE(inp % w_concorr_c, 3)
    out % m___f2dace_SOA_w_concorr_c_d_0_s_297 = LBOUND(inp % w_concorr_c, 1)
    out % m___f2dace_SOA_w_concorr_c_d_1_s_298 = LBOUND(inp % w_concorr_c, 2)
    out % m___f2dace_SOA_w_concorr_c_d_2_s_299 = LBOUND(inp % w_concorr_c, 3)
  END SUBROUTINE ctor_t_nh_diag
  SUBROUTINE ctor_t_grid_edges(inp, out, initalloc)
    TYPE(t_grid_edges), INTENT(IN) :: inp
    TYPE(glue_t_grid_edges), INTENT(INOUT) :: out
    LOGICAL, INTENT(IN) :: initalloc
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_area_edge(:, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_cell_blk(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_cell_idx(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_end_block(:)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_end_index(:)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_f_e(:, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_fn_e(:, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_ft_e(:, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_inv_dual_edge_length(:, :)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_inv_primal_edge_length(:, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_quad_blk(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_quad_idx(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_start_block(:)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_start_index(:)
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_tangent_orientation(:, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_vertex_blk(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_vertex_idx(:, :, :)
    out % m___f2dace_SA_dual_normal_cell_d_0_s = 0
    out % m___f2dace_SA_dual_normal_cell_d_1_s = 0
    out % m___f2dace_SA_dual_normal_cell_d_2_s = 0
    out % m___f2dace_SA_primal_normal_cell_d_0_s = 0
    out % m___f2dace_SA_primal_normal_cell_d_1_s = 0
    out % m___f2dace_SA_primal_normal_cell_d_2_s = 0
    out % m___f2dace_SA_refin_ctrl_d_0_s = 0
    out % m___f2dace_SA_refin_ctrl_d_1_s = 0
    out % m___f2dace_SOA_dual_normal_cell_d_0_s = 0
    out % m___f2dace_SOA_dual_normal_cell_d_1_s = 0
    out % m___f2dace_SOA_dual_normal_cell_d_2_s = 0
    out % m___f2dace_SOA_primal_normal_cell_d_0_s = 0
    out % m___f2dace_SOA_primal_normal_cell_d_1_s = 0
    out % m___f2dace_SOA_primal_normal_cell_d_2_s = 0
    out % m___f2dace_SOA_refin_ctrl_d_0_s = 0
    out % m___f2dace_SOA_refin_ctrl_d_1_s = 0
    out % m_dual_normal_cell = c_null_ptr
    out % m_primal_normal_cell = c_null_ptr
    out % m_refin_ctrl = c_null_ptr
    IF (initalloc .AND. .NOT. ALLOCATED(a_area_edge)) ALLOCATE(a_area_edge(SIZE(inp % area_edge, 1), SIZE(inp % area_edge, 2)))
    a_area_edge = inp % area_edge
    out % m_area_edge = c_loc(a_area_edge)
    out % m___f2dace_SA_area_edge_d_0_s_188 = SIZE(inp % area_edge, 1)
    out % m___f2dace_SA_area_edge_d_1_s_189 = SIZE(inp % area_edge, 2)
    out % m___f2dace_SOA_area_edge_d_0_s_188 = LBOUND(inp % area_edge, 1)
    out % m___f2dace_SOA_area_edge_d_1_s_189 = LBOUND(inp % area_edge, 2)
    IF (initalloc .AND. .NOT. ALLOCATED(a_cell_blk)) ALLOCATE(a_cell_blk(SIZE(inp % cell_blk, 1), SIZE(inp % cell_blk, 2), SIZE(inp % cell_blk, 3)))
    a_cell_blk = inp % cell_blk
    out % m_cell_blk = c_loc(a_cell_blk)
    out % m___f2dace_SA_cell_blk_d_0_s_167 = SIZE(inp % cell_blk, 1)
    out % m___f2dace_SA_cell_blk_d_1_s_168 = SIZE(inp % cell_blk, 2)
    out % m___f2dace_SA_cell_blk_d_2_s_169 = SIZE(inp % cell_blk, 3)
    out % m___f2dace_SOA_cell_blk_d_0_s_167 = LBOUND(inp % cell_blk, 1)
    out % m___f2dace_SOA_cell_blk_d_1_s_168 = LBOUND(inp % cell_blk, 2)
    out % m___f2dace_SOA_cell_blk_d_2_s_169 = LBOUND(inp % cell_blk, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_cell_idx)) ALLOCATE(a_cell_idx(SIZE(inp % cell_idx, 1), SIZE(inp % cell_idx, 2), SIZE(inp % cell_idx, 3)))
    a_cell_idx = inp % cell_idx
    out % m_cell_idx = c_loc(a_cell_idx)
    out % m___f2dace_SA_cell_idx_d_0_s_164 = SIZE(inp % cell_idx, 1)
    out % m___f2dace_SA_cell_idx_d_1_s_165 = SIZE(inp % cell_idx, 2)
    out % m___f2dace_SA_cell_idx_d_2_s_166 = SIZE(inp % cell_idx, 3)
    out % m___f2dace_SOA_cell_idx_d_0_s_164 = LBOUND(inp % cell_idx, 1)
    out % m___f2dace_SOA_cell_idx_d_1_s_165 = LBOUND(inp % cell_idx, 2)
    out % m___f2dace_SOA_cell_idx_d_2_s_166 = LBOUND(inp % cell_idx, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_end_block)) ALLOCATE(a_end_block(SIZE(inp % end_block, 1)))
    a_end_block = inp % end_block
    out % m_end_block = c_loc(a_end_block)
    out % m___f2dace_SA_end_block_d_0_s_199 = SIZE(inp % end_block, 1)
    out % m___f2dace_SOA_end_block_d_0_s_199 = LBOUND(inp % end_block, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_end_index)) ALLOCATE(a_end_index(SIZE(inp % end_index, 1)))
    a_end_index = inp % end_index
    out % m_end_index = c_loc(a_end_index)
    out % m___f2dace_SA_end_index_d_0_s_197 = SIZE(inp % end_index, 1)
    out % m___f2dace_SOA_end_index_d_0_s_197 = LBOUND(inp % end_index, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_f_e)) ALLOCATE(a_f_e(SIZE(inp % f_e, 1), SIZE(inp % f_e, 2)))
    a_f_e = inp % f_e
    out % m_f_e = c_loc(a_f_e)
    out % m___f2dace_SA_f_e_d_0_s_190 = SIZE(inp % f_e, 1)
    out % m___f2dace_SA_f_e_d_1_s_191 = SIZE(inp % f_e, 2)
    out % m___f2dace_SOA_f_e_d_0_s_190 = LBOUND(inp % f_e, 1)
    out % m___f2dace_SOA_f_e_d_1_s_191 = LBOUND(inp % f_e, 2)
    IF (initalloc .AND. .NOT. ALLOCATED(a_fn_e)) ALLOCATE(a_fn_e(SIZE(inp % fn_e, 1), SIZE(inp % fn_e, 2)))
    a_fn_e = inp % fn_e
    out % m_fn_e = c_loc(a_fn_e)
    out % m___f2dace_SA_fn_e_d_0_s_192 = SIZE(inp % fn_e, 1)
    out % m___f2dace_SA_fn_e_d_1_s_193 = SIZE(inp % fn_e, 2)
    out % m___f2dace_SOA_fn_e_d_0_s_192 = LBOUND(inp % fn_e, 1)
    out % m___f2dace_SOA_fn_e_d_1_s_193 = LBOUND(inp % fn_e, 2)
    IF (initalloc .AND. .NOT. ALLOCATED(a_ft_e)) ALLOCATE(a_ft_e(SIZE(inp % ft_e, 1), SIZE(inp % ft_e, 2)))
    a_ft_e = inp % ft_e
    out % m_ft_e = c_loc(a_ft_e)
    out % m___f2dace_SA_ft_e_d_0_s_194 = SIZE(inp % ft_e, 1)
    out % m___f2dace_SA_ft_e_d_1_s_195 = SIZE(inp % ft_e, 2)
    out % m___f2dace_SOA_ft_e_d_0_s_194 = LBOUND(inp % ft_e, 1)
    out % m___f2dace_SOA_ft_e_d_1_s_195 = LBOUND(inp % ft_e, 2)
    IF (initalloc .AND. .NOT. ALLOCATED(a_inv_dual_edge_length)) ALLOCATE(a_inv_dual_edge_length(SIZE(inp % inv_dual_edge_length, 1), SIZE(inp % inv_dual_edge_length, 2)))
    a_inv_dual_edge_length = inp % inv_dual_edge_length
    out % m_inv_dual_edge_length = c_loc(a_inv_dual_edge_length)
    out % m___f2dace_SA_inv_dual_edge_length_d_0_s_186 = SIZE(inp % inv_dual_edge_length, 1)
    out % m___f2dace_SA_inv_dual_edge_length_d_1_s_187 = SIZE(inp % inv_dual_edge_length, 2)
    out % m___f2dace_SOA_inv_dual_edge_length_d_0_s_186 = LBOUND(inp % inv_dual_edge_length, 1)
    out % m___f2dace_SOA_inv_dual_edge_length_d_1_s_187 = LBOUND(inp % inv_dual_edge_length, 2)
    IF (initalloc .AND. .NOT. ALLOCATED(a_inv_primal_edge_length)) ALLOCATE(a_inv_primal_edge_length(SIZE(inp % inv_primal_edge_length, 1), SIZE(inp % inv_primal_edge_length, 2)))
    a_inv_primal_edge_length = inp % inv_primal_edge_length
    out % m_inv_primal_edge_length = c_loc(a_inv_primal_edge_length)
    out % m___f2dace_SA_inv_primal_edge_length_d_0_s_184 = SIZE(inp % inv_primal_edge_length, 1)
    out % m___f2dace_SA_inv_primal_edge_length_d_1_s_185 = SIZE(inp % inv_primal_edge_length, 2)
    out % m___f2dace_SOA_inv_primal_edge_length_d_0_s_184 = LBOUND(inp % inv_primal_edge_length, 1)
    out % m___f2dace_SOA_inv_primal_edge_length_d_1_s_185 = LBOUND(inp % inv_primal_edge_length, 2)
    IF (initalloc .AND. .NOT. ALLOCATED(a_quad_blk)) ALLOCATE(a_quad_blk(SIZE(inp % quad_blk, 1), SIZE(inp % quad_blk, 2), SIZE(inp % quad_blk, 3)))
    a_quad_blk = inp % quad_blk
    out % m_quad_blk = c_loc(a_quad_blk)
    out % m___f2dace_SA_quad_blk_d_0_s_181 = SIZE(inp % quad_blk, 1)
    out % m___f2dace_SA_quad_blk_d_1_s_182 = SIZE(inp % quad_blk, 2)
    out % m___f2dace_SA_quad_blk_d_2_s_183 = SIZE(inp % quad_blk, 3)
    out % m___f2dace_SOA_quad_blk_d_0_s_181 = LBOUND(inp % quad_blk, 1)
    out % m___f2dace_SOA_quad_blk_d_1_s_182 = LBOUND(inp % quad_blk, 2)
    out % m___f2dace_SOA_quad_blk_d_2_s_183 = LBOUND(inp % quad_blk, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_quad_idx)) ALLOCATE(a_quad_idx(SIZE(inp % quad_idx, 1), SIZE(inp % quad_idx, 2), SIZE(inp % quad_idx, 3)))
    a_quad_idx = inp % quad_idx
    out % m_quad_idx = c_loc(a_quad_idx)
    out % m___f2dace_SA_quad_idx_d_0_s_178 = SIZE(inp % quad_idx, 1)
    out % m___f2dace_SA_quad_idx_d_1_s_179 = SIZE(inp % quad_idx, 2)
    out % m___f2dace_SA_quad_idx_d_2_s_180 = SIZE(inp % quad_idx, 3)
    out % m___f2dace_SOA_quad_idx_d_0_s_178 = LBOUND(inp % quad_idx, 1)
    out % m___f2dace_SOA_quad_idx_d_1_s_179 = LBOUND(inp % quad_idx, 2)
    out % m___f2dace_SOA_quad_idx_d_2_s_180 = LBOUND(inp % quad_idx, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_start_block)) ALLOCATE(a_start_block(SIZE(inp % start_block, 1)))
    a_start_block = inp % start_block
    out % m_start_block = c_loc(a_start_block)
    out % m___f2dace_SA_start_block_d_0_s_198 = SIZE(inp % start_block, 1)
    out % m___f2dace_SOA_start_block_d_0_s_198 = LBOUND(inp % start_block, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_start_index)) ALLOCATE(a_start_index(SIZE(inp % start_index, 1)))
    a_start_index = inp % start_index
    out % m_start_index = c_loc(a_start_index)
    out % m___f2dace_SA_start_index_d_0_s_196 = SIZE(inp % start_index, 1)
    out % m___f2dace_SOA_start_index_d_0_s_196 = LBOUND(inp % start_index, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_tangent_orientation)) ALLOCATE(a_tangent_orientation(SIZE(inp % tangent_orientation, 1), SIZE(inp % tangent_orientation, 2)))
    a_tangent_orientation = inp % tangent_orientation
    out % m_tangent_orientation = c_loc(a_tangent_orientation)
    out % m___f2dace_SA_tangent_orientation_d_0_s_176 = SIZE(inp % tangent_orientation, 1)
    out % m___f2dace_SA_tangent_orientation_d_1_s_177 = SIZE(inp % tangent_orientation, 2)
    out % m___f2dace_SOA_tangent_orientation_d_0_s_176 = LBOUND(inp % tangent_orientation, 1)
    out % m___f2dace_SOA_tangent_orientation_d_1_s_177 = LBOUND(inp % tangent_orientation, 2)
    IF (initalloc .AND. .NOT. ALLOCATED(a_vertex_blk)) ALLOCATE(a_vertex_blk(SIZE(inp % vertex_blk, 1), SIZE(inp % vertex_blk, 2), SIZE(inp % vertex_blk, 3)))
    a_vertex_blk = inp % vertex_blk
    out % m_vertex_blk = c_loc(a_vertex_blk)
    out % m___f2dace_SA_vertex_blk_d_0_s_173 = SIZE(inp % vertex_blk, 1)
    out % m___f2dace_SA_vertex_blk_d_1_s_174 = SIZE(inp % vertex_blk, 2)
    out % m___f2dace_SA_vertex_blk_d_2_s_175 = SIZE(inp % vertex_blk, 3)
    out % m___f2dace_SOA_vertex_blk_d_0_s_173 = LBOUND(inp % vertex_blk, 1)
    out % m___f2dace_SOA_vertex_blk_d_1_s_174 = LBOUND(inp % vertex_blk, 2)
    out % m___f2dace_SOA_vertex_blk_d_2_s_175 = LBOUND(inp % vertex_blk, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_vertex_idx)) ALLOCATE(a_vertex_idx(SIZE(inp % vertex_idx, 1), SIZE(inp % vertex_idx, 2), SIZE(inp % vertex_idx, 3)))
    a_vertex_idx = inp % vertex_idx
    out % m_vertex_idx = c_loc(a_vertex_idx)
    out % m___f2dace_SA_vertex_idx_d_0_s_170 = SIZE(inp % vertex_idx, 1)
    out % m___f2dace_SA_vertex_idx_d_1_s_171 = SIZE(inp % vertex_idx, 2)
    out % m___f2dace_SA_vertex_idx_d_2_s_172 = SIZE(inp % vertex_idx, 3)
    out % m___f2dace_SOA_vertex_idx_d_0_s_170 = LBOUND(inp % vertex_idx, 1)
    out % m___f2dace_SOA_vertex_idx_d_1_s_171 = LBOUND(inp % vertex_idx, 2)
    out % m___f2dace_SOA_vertex_idx_d_2_s_172 = LBOUND(inp % vertex_idx, 3)
  END SUBROUTINE ctor_t_grid_edges
  SUBROUTINE ctor_t_grid_cells(inp, out, initalloc)
    TYPE(t_grid_cells), INTENT(IN) :: inp
    TYPE(glue_t_grid_cells), INTENT(INOUT) :: out
    LOGICAL, INTENT(IN) :: initalloc
    REAL(KIND = c_double), ALLOCATABLE, TARGET, SAVE :: a_area(:, :)
    TYPE(glue_t_grid_domain_decomp_info), ALLOCATABLE, TARGET, SAVE :: a_decomp_info
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_edge_blk(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_edge_idx(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_end_block(:)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_end_index(:)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_neighbor_blk(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_neighbor_idx(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_start_block(:)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_start_index(:)
    out % m___f2dace_SA_end_blk_d_0_s = 0
    out % m___f2dace_SA_end_blk_d_1_s = 0
    out % m___f2dace_SA_start_blk_d_0_s = 0
    out % m___f2dace_SA_start_blk_d_1_s = 0
    out % m___f2dace_SOA_end_blk_d_0_s = 0
    out % m___f2dace_SOA_end_blk_d_1_s = 0
    out % m___f2dace_SOA_start_blk_d_0_s = 0
    out % m___f2dace_SOA_start_blk_d_1_s = 0
    out % m_end_blk = c_null_ptr
    out % m_start_blk = c_null_ptr
    IF (initalloc .AND. .NOT. ALLOCATED(a_area)) ALLOCATE(a_area(SIZE(inp % area, 1), SIZE(inp % area, 2)))
    a_area = inp % area
    out % m_area = c_loc(a_area)
    out % m___f2dace_SA_area_d_0_s_158 = SIZE(inp % area, 1)
    out % m___f2dace_SA_area_d_1_s_159 = SIZE(inp % area, 2)
    out % m___f2dace_SOA_area_d_0_s_158 = LBOUND(inp % area, 1)
    out % m___f2dace_SOA_area_d_1_s_159 = LBOUND(inp % area, 2)
    IF (initalloc .AND. .NOT. ALLOCATED(a_decomp_info)) ALLOCATE(a_decomp_info)
    CALL ctor(inp % decomp_info, a_decomp_info, initalloc)
    out % m_decomp_info = c_loc(a_decomp_info)
    IF (initalloc .AND. .NOT. ALLOCATED(a_edge_blk)) ALLOCATE(a_edge_blk(SIZE(inp % edge_blk, 1), SIZE(inp % edge_blk, 2), SIZE(inp % edge_blk, 3)))
    a_edge_blk = inp % edge_blk
    out % m_edge_blk = c_loc(a_edge_blk)
    out % m___f2dace_SA_edge_blk_d_0_s_155 = SIZE(inp % edge_blk, 1)
    out % m___f2dace_SA_edge_blk_d_1_s_156 = SIZE(inp % edge_blk, 2)
    out % m___f2dace_SA_edge_blk_d_2_s_157 = SIZE(inp % edge_blk, 3)
    out % m___f2dace_SOA_edge_blk_d_0_s_155 = LBOUND(inp % edge_blk, 1)
    out % m___f2dace_SOA_edge_blk_d_1_s_156 = LBOUND(inp % edge_blk, 2)
    out % m___f2dace_SOA_edge_blk_d_2_s_157 = LBOUND(inp % edge_blk, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_edge_idx)) ALLOCATE(a_edge_idx(SIZE(inp % edge_idx, 1), SIZE(inp % edge_idx, 2), SIZE(inp % edge_idx, 3)))
    a_edge_idx = inp % edge_idx
    out % m_edge_idx = c_loc(a_edge_idx)
    out % m___f2dace_SA_edge_idx_d_0_s_152 = SIZE(inp % edge_idx, 1)
    out % m___f2dace_SA_edge_idx_d_1_s_153 = SIZE(inp % edge_idx, 2)
    out % m___f2dace_SA_edge_idx_d_2_s_154 = SIZE(inp % edge_idx, 3)
    out % m___f2dace_SOA_edge_idx_d_0_s_152 = LBOUND(inp % edge_idx, 1)
    out % m___f2dace_SOA_edge_idx_d_1_s_153 = LBOUND(inp % edge_idx, 2)
    out % m___f2dace_SOA_edge_idx_d_2_s_154 = LBOUND(inp % edge_idx, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_end_block)) ALLOCATE(a_end_block(SIZE(inp % end_block, 1)))
    a_end_block = inp % end_block
    out % m_end_block = c_loc(a_end_block)
    out % m___f2dace_SA_end_block_d_0_s_163 = SIZE(inp % end_block, 1)
    out % m___f2dace_SOA_end_block_d_0_s_163 = LBOUND(inp % end_block, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_end_index)) ALLOCATE(a_end_index(SIZE(inp % end_index, 1)))
    a_end_index = inp % end_index
    out % m_end_index = c_loc(a_end_index)
    out % m___f2dace_SA_end_index_d_0_s_161 = SIZE(inp % end_index, 1)
    out % m___f2dace_SOA_end_index_d_0_s_161 = LBOUND(inp % end_index, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_neighbor_blk)) ALLOCATE(a_neighbor_blk(SIZE(inp % neighbor_blk, 1), SIZE(inp % neighbor_blk, 2), SIZE(inp % neighbor_blk, 3)))
    a_neighbor_blk = inp % neighbor_blk
    out % m_neighbor_blk = c_loc(a_neighbor_blk)
    out % m___f2dace_SA_neighbor_blk_d_0_s_149 = SIZE(inp % neighbor_blk, 1)
    out % m___f2dace_SA_neighbor_blk_d_1_s_150 = SIZE(inp % neighbor_blk, 2)
    out % m___f2dace_SA_neighbor_blk_d_2_s_151 = SIZE(inp % neighbor_blk, 3)
    out % m___f2dace_SOA_neighbor_blk_d_0_s_149 = LBOUND(inp % neighbor_blk, 1)
    out % m___f2dace_SOA_neighbor_blk_d_1_s_150 = LBOUND(inp % neighbor_blk, 2)
    out % m___f2dace_SOA_neighbor_blk_d_2_s_151 = LBOUND(inp % neighbor_blk, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_neighbor_idx)) ALLOCATE(a_neighbor_idx(SIZE(inp % neighbor_idx, 1), SIZE(inp % neighbor_idx, 2), SIZE(inp % neighbor_idx, 3)))
    a_neighbor_idx = inp % neighbor_idx
    out % m_neighbor_idx = c_loc(a_neighbor_idx)
    out % m___f2dace_SA_neighbor_idx_d_0_s_146 = SIZE(inp % neighbor_idx, 1)
    out % m___f2dace_SA_neighbor_idx_d_1_s_147 = SIZE(inp % neighbor_idx, 2)
    out % m___f2dace_SA_neighbor_idx_d_2_s_148 = SIZE(inp % neighbor_idx, 3)
    out % m___f2dace_SOA_neighbor_idx_d_0_s_146 = LBOUND(inp % neighbor_idx, 1)
    out % m___f2dace_SOA_neighbor_idx_d_1_s_147 = LBOUND(inp % neighbor_idx, 2)
    out % m___f2dace_SOA_neighbor_idx_d_2_s_148 = LBOUND(inp % neighbor_idx, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_start_block)) ALLOCATE(a_start_block(SIZE(inp % start_block, 1)))
    a_start_block = inp % start_block
    out % m_start_block = c_loc(a_start_block)
    out % m___f2dace_SA_start_block_d_0_s_162 = SIZE(inp % start_block, 1)
    out % m___f2dace_SOA_start_block_d_0_s_162 = LBOUND(inp % start_block, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_start_index)) ALLOCATE(a_start_index(SIZE(inp % start_index, 1)))
    a_start_index = inp % start_index
    out % m_start_index = c_loc(a_start_index)
    out % m___f2dace_SA_start_index_d_0_s_160 = SIZE(inp % start_index, 1)
    out % m___f2dace_SOA_start_index_d_0_s_160 = LBOUND(inp % start_index, 1)
  END SUBROUTINE ctor_t_grid_cells
  SUBROUTINE ctor_t_grid_domain_decomp_info(inp, out, initalloc)
    TYPE(t_grid_domain_decomp_info), INTENT(IN) :: inp
    TYPE(glue_t_grid_domain_decomp_info), INTENT(INOUT) :: out
    LOGICAL, INTENT(IN) :: initalloc
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_owner_mask(:, :)
    IF (initalloc .AND. .NOT. ALLOCATED(a_owner_mask)) ALLOCATE(a_owner_mask(SIZE(inp % owner_mask, 1), SIZE(inp % owner_mask, 2)))
    a_owner_mask = inp % owner_mask
    out % m_owner_mask = c_loc(a_owner_mask)
    out % m___f2dace_SA_owner_mask_d_0_s_2 = SIZE(inp % owner_mask, 1)
    out % m___f2dace_SA_owner_mask_d_1_s_3 = SIZE(inp % owner_mask, 2)
    out % m___f2dace_SOA_owner_mask_d_0_s_2 = LBOUND(inp % owner_mask, 1)
    out % m___f2dace_SOA_owner_mask_d_1_s_3 = LBOUND(inp % owner_mask, 2)
  END SUBROUTINE ctor_t_grid_domain_decomp_info
  SUBROUTINE ctor_t_grid_vertices(inp, out, initalloc)
    TYPE(t_grid_vertices), INTENT(IN) :: inp
    TYPE(glue_t_grid_vertices), INTENT(INOUT) :: out
    LOGICAL, INTENT(IN) :: initalloc
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_cell_blk(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_cell_idx(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_edge_blk(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_edge_idx(:, :, :)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_end_block(:)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_end_index(:)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_start_block(:)
    INTEGER(KIND = c_int), ALLOCATABLE, TARGET, SAVE :: a_start_index(:)
    IF (initalloc .AND. .NOT. ALLOCATED(a_cell_blk)) ALLOCATE(a_cell_blk(SIZE(inp % cell_blk, 1), SIZE(inp % cell_blk, 2), SIZE(inp % cell_blk, 3)))
    a_cell_blk = inp % cell_blk
    out % m_cell_blk = c_loc(a_cell_blk)
    out % m___f2dace_SA_cell_blk_d_0_s_167 = SIZE(inp % cell_blk, 1)
    out % m___f2dace_SA_cell_blk_d_1_s_168 = SIZE(inp % cell_blk, 2)
    out % m___f2dace_SA_cell_blk_d_2_s_169 = SIZE(inp % cell_blk, 3)
    out % m___f2dace_SOA_cell_blk_d_0_s_167 = LBOUND(inp % cell_blk, 1)
    out % m___f2dace_SOA_cell_blk_d_1_s_168 = LBOUND(inp % cell_blk, 2)
    out % m___f2dace_SOA_cell_blk_d_2_s_169 = LBOUND(inp % cell_blk, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_cell_idx)) ALLOCATE(a_cell_idx(SIZE(inp % cell_idx, 1), SIZE(inp % cell_idx, 2), SIZE(inp % cell_idx, 3)))
    a_cell_idx = inp % cell_idx
    out % m_cell_idx = c_loc(a_cell_idx)
    out % m___f2dace_SA_cell_idx_d_0_s_164 = SIZE(inp % cell_idx, 1)
    out % m___f2dace_SA_cell_idx_d_1_s_165 = SIZE(inp % cell_idx, 2)
    out % m___f2dace_SA_cell_idx_d_2_s_166 = SIZE(inp % cell_idx, 3)
    out % m___f2dace_SOA_cell_idx_d_0_s_164 = LBOUND(inp % cell_idx, 1)
    out % m___f2dace_SOA_cell_idx_d_1_s_165 = LBOUND(inp % cell_idx, 2)
    out % m___f2dace_SOA_cell_idx_d_2_s_166 = LBOUND(inp % cell_idx, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_edge_blk)) ALLOCATE(a_edge_blk(SIZE(inp % edge_blk, 1), SIZE(inp % edge_blk, 2), SIZE(inp % edge_blk, 3)))
    a_edge_blk = inp % edge_blk
    out % m_edge_blk = c_loc(a_edge_blk)
    out % m___f2dace_SA_edge_blk_d_0_s_155 = SIZE(inp % edge_blk, 1)
    out % m___f2dace_SA_edge_blk_d_1_s_156 = SIZE(inp % edge_blk, 2)
    out % m___f2dace_SA_edge_blk_d_2_s_157 = SIZE(inp % edge_blk, 3)
    out % m___f2dace_SOA_edge_blk_d_0_s_155 = LBOUND(inp % edge_blk, 1)
    out % m___f2dace_SOA_edge_blk_d_1_s_156 = LBOUND(inp % edge_blk, 2)
    out % m___f2dace_SOA_edge_blk_d_2_s_157 = LBOUND(inp % edge_blk, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_edge_idx)) ALLOCATE(a_edge_idx(SIZE(inp % edge_idx, 1), SIZE(inp % edge_idx, 2), SIZE(inp % edge_idx, 3)))
    a_edge_idx = inp % edge_idx
    out % m_edge_idx = c_loc(a_edge_idx)
    out % m___f2dace_SA_edge_idx_d_0_s_152 = SIZE(inp % edge_idx, 1)
    out % m___f2dace_SA_edge_idx_d_1_s_153 = SIZE(inp % edge_idx, 2)
    out % m___f2dace_SA_edge_idx_d_2_s_154 = SIZE(inp % edge_idx, 3)
    out % m___f2dace_SOA_edge_idx_d_0_s_152 = LBOUND(inp % edge_idx, 1)
    out % m___f2dace_SOA_edge_idx_d_1_s_153 = LBOUND(inp % edge_idx, 2)
    out % m___f2dace_SOA_edge_idx_d_2_s_154 = LBOUND(inp % edge_idx, 3)
    IF (initalloc .AND. .NOT. ALLOCATED(a_end_block)) ALLOCATE(a_end_block(SIZE(inp % end_block, 1)))
    a_end_block = inp % end_block
    out % m_end_block = c_loc(a_end_block)
    out % m___f2dace_SA_end_block_d_0_s_163 = SIZE(inp % end_block, 1)
    out % m___f2dace_SOA_end_block_d_0_s_163 = LBOUND(inp % end_block, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_end_index)) ALLOCATE(a_end_index(SIZE(inp % end_index, 1)))
    a_end_index = inp % end_index
    out % m_end_index = c_loc(a_end_index)
    out % m___f2dace_SA_end_index_d_0_s_161 = SIZE(inp % end_index, 1)
    out % m___f2dace_SOA_end_index_d_0_s_161 = LBOUND(inp % end_index, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_start_block)) ALLOCATE(a_start_block(SIZE(inp % start_block, 1)))
    a_start_block = inp % start_block
    out % m_start_block = c_loc(a_start_block)
    out % m___f2dace_SA_start_block_d_0_s_162 = SIZE(inp % start_block, 1)
    out % m___f2dace_SOA_start_block_d_0_s_162 = LBOUND(inp % start_block, 1)
    IF (initalloc .AND. .NOT. ALLOCATED(a_start_index)) ALLOCATE(a_start_index(SIZE(inp % start_index, 1)))
    a_start_index = inp % start_index
    out % m_start_index = c_loc(a_start_index)
    out % m___f2dace_SA_start_index_d_0_s_160 = SIZE(inp % start_index, 1)
    out % m___f2dace_SOA_start_index_d_0_s_160 = LBOUND(inp % start_index, 1)
  END SUBROUTINE ctor_t_grid_vertices
END MODULE f90_glue_vt_serde
MODULE vt_serde
  IMPLICIT NONE
  INTERFACE serialize
    MODULE PROCEDURE W_integer1, W_integer2, W_integer4, W_integer8, W_integer__4_R_1, W_integer__4_R_3, W_logical, W_logical_R_2, W_real4, W_real8, W_real__8_R_1, W_real__8_R_2, W_real__8_R_3, W_real__8_R_4, W_string, W_t_grid_cells, W_t_grid_domain_decomp_info, W_t_grid_edges, W_t_grid_vertices, W_t_int_state, W_t_nh_diag, W_t_nh_metrics, W_t_nh_prog, W_t_patch
  END INTERFACE serialize
  INTEGER :: generation = 0
  INTEGER :: vt_generation = 0
  INTEGER :: dycore_generation = 0
  INTEGER :: physics_generation = 0
  INTEGER :: dyn_substeps = 0
  LOGICAL :: do_serialize = .false.
  ! Runtime-tunable knobs (overridable via env vars at first physics_tic)
  INTEGER :: ndyn_substeps_override = 0       ! 0 = don't override; >0 = force this value
  INTEGER :: serde_gen_start        = 0
  INTEGER :: serde_gen_end          = 51
  INTEGER :: serde_gen_stride       = 1
  INTEGER, PARAMETER :: serde_gen_maxset = 64
  INTEGER :: serde_gen_nset         = 0   ! >0 => use explicit gen list, ignore stride window
  INTEGER :: serde_gen_set(64)      = 0
  LOGICAL :: use_vt_gpu             = .true.  ! dispatch solve_nh's velocity_tendencies to VT's libvelocity.so
  LOGICAL, PRIVATE :: vt_serde_inited = .false.
  CONTAINS
  ! SUBROUTINE tic
  !   generation = generation + 1
  ! END SUBROUTINE tic
  SUBROUTINE vt_serde_init
    CHARACTER(LEN = 64) :: v
    INTEGER :: ios
    IF (vt_serde_inited) RETURN
    vt_serde_inited = .true.
    CALL GET_ENVIRONMENT_VARIABLE('NDYN_SUBSTEPS_OVERRIDE', v)
    IF (LEN_TRIM(v) > 0) READ (v, *, IOSTAT = ios) ndyn_substeps_override
    CALL GET_ENVIRONMENT_VARIABLE('SERDE_GEN_START',        v)
    IF (LEN_TRIM(v) > 0) READ (v, *, IOSTAT = ios) serde_gen_start
    CALL GET_ENVIRONMENT_VARIABLE('SERDE_GEN_END',          v)
    IF (LEN_TRIM(v) > 0) READ (v, *, IOSTAT = ios) serde_gen_end
    CALL GET_ENVIRONMENT_VARIABLE('SERDE_GEN_STRIDE',       v)
    IF (LEN_TRIM(v) > 0) READ (v, *, IOSTAT = ios) serde_gen_stride
    CALL parse_serde_gen_list
    CALL GET_ENVIRONMENT_VARIABLE('USE_VT_GPU',             v)
    IF (LEN_TRIM(v) > 0) THEN
      SELECT CASE (TRIM(ADJUSTL(v)))
        CASE ('0', 'F', 'f', 'false', 'FALSE', '.false.', '.FALSE.')
          use_vt_gpu = .false.
        CASE DEFAULT
          use_vt_gpu = .true.
      END SELECT
    END IF
  END SUBROUTINE vt_serde_init
  SUBROUTINE parse_serde_gen_list
    ! Read SERDE_GEN_LIST (comma/space separated ints) into serde_gen_set.
    ! When serde_gen_nset>0, the gate uses set-membership (logspace dumps).
    CHARACTER(LEN = 2048) :: vl
    INTEGER :: i, n, ios2
    LOGICAL :: prev_blank, cur_blank
    serde_gen_nset = 0
    CALL GET_ENVIRONMENT_VARIABLE('SERDE_GEN_LIST', vl)
    IF (LEN_TRIM(vl) == 0) RETURN
    DO i = 1, LEN_TRIM(vl)
      IF (vl(i:i) == ',') vl(i:i) = ' '
    END DO
    n = 0
    prev_blank = .true.
    DO i = 1, LEN_TRIM(vl)
      cur_blank = (vl(i:i) == ' ')
      IF (prev_blank .AND. .NOT. cur_blank) n = n + 1
      prev_blank = cur_blank
    END DO
    IF (n > serde_gen_maxset) n = serde_gen_maxset
    IF (n > 0) READ (vl, *, IOSTAT = ios2) (serde_gen_set(i), i = 1, n)
    IF (n > 0 .AND. ios2 == 0) serde_gen_nset = n
  END SUBROUTINE parse_serde_gen_list
  SUBROUTINE vt_tic
    vt_generation = vt_generation + 1
  END SUBROUTINE vt_tic
  SUBROUTINE dycore_tic
    dycore_generation = dycore_generation + 1
    vt_generation = 0
  END SUBROUTINE dycore_tic
  SUBROUTINE physics_tic
    CALL vt_serde_init
    physics_generation = physics_generation + 1
    dycore_generation = 0
  END SUBROUTINE physics_tic
  FUNCTION cat(prefix, asis) RESULT(path)
    CHARACTER(LEN = *), INTENT(IN) :: prefix
    CHARACTER(LEN = :), ALLOCATABLE :: path
    CHARACTER(LEN = 50) :: gen
    LOGICAL, INTENT(IN) :: asis
    IF (asis) THEN
      path = prefix
    ELSE
      WRITE(gen, '(A,I0,A,I0,A,I0,A,I0)') 'p', physics_generation, '.d', dycore_generation, '.vt', vt_generation, '.ss', dyn_substeps
      path = prefix // '.' // TRIM(gen) // ".data"
    END IF
  END FUNCTION cat
  FUNCTION at(prefix, asis) RESULT(io)
    CHARACTER(LEN = *), INTENT(IN) :: prefix
    INTEGER :: io
    LOGICAL, OPTIONAL, INTENT(IN) :: asis
    LOGICAL :: asis_local
    asis_local = .FALSE.
    IF (PRESENT(asis)) asis_local = asis
    OPEN(NEWUNIT = io, FILE = cat(prefix, asis_local), STATUS = "replace", ACTION = "write")
  END FUNCTION at
  SUBROUTINE W_string(io, x, cleanup, nline)
    INTEGER :: io
    CHARACTER(LEN = *), INTENT(IN) :: x
    INTEGER :: i, xend
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    xend = LEN(x)
    DO i = 1, LEN(x)
      IF (x(i : i) == CHAR(0)) THEN
        xend = i - 1
        EXIT
      END IF
    END DO
    WRITE(io, '(A)', ADVANCE = 'no') TRIM(x(1 : xend))
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_string
  SUBROUTINE W_logical(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    LOGICAL, INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    y = MERGE(1, 0, x)
    WRITE(io, '(g0)', ADVANCE = 'no') y
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_logical
  SUBROUTINE W_integer1(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    INTEGER(KIND = 1), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(io, '(g0)', ADVANCE = 'no') x
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer1
  SUBROUTINE W_integer2(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    INTEGER(KIND = 2), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(io, '(g0)', ADVANCE = 'no') x
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer2
  SUBROUTINE W_integer4(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    INTEGER(KIND = 4), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(io, '(g0)', ADVANCE = 'no') x
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer4
  SUBROUTINE W_integer8(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    INTEGER(KIND = 8), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(io, '(g0)', ADVANCE = 'no') x
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer8
  SUBROUTINE W_real4(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    REAL(KIND = 4), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(buf, '(e28.20)') x
    WRITE(io, '(A)', ADVANCE = 'no') TRIM(ADJUSTL(buf))
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real4
  SUBROUTINE W_real8(io, x, cleanup, nline)
    CHARACTER(LEN = 50) :: buf
    INTEGER :: io
    REAL(KIND = 8), INTENT(IN) :: x
    INTEGER :: y
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    WRITE(buf, '(e28.20)') x
    WRITE(io, '(A)', ADVANCE = 'no') TRIM(ADJUSTL(buf))
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real8
  SUBROUTINE W_t_grid_domain_decomp_info(io, x, cleanup, nline)
    USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
    INTEGER :: io
    TYPE(t_grid_domain_decomp_info), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# owner_mask', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 2, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, SIZE(x % owner_mask, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, LBOUND(x % owner_mask, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % owner_mask, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_grid_domain_decomp_info
  SUBROUTINE W_t_int_state(io, x, cleanup, nline)
    USE mo_intp_data_strc, ONLY: t_int_state
    INTEGER :: io
    TYPE(t_int_state), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# c_lin_e', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % c_lin_e, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % c_lin_e, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % c_lin_e, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# e_bln_c_s', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % e_bln_c_s, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % e_bln_c_s, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % e_bln_c_s, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# cells_aw_verts', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % cells_aw_verts, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % cells_aw_verts, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % cells_aw_verts, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# rbf_vec_coeff_e', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % rbf_vec_coeff_e, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % rbf_vec_coeff_e, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % rbf_vec_coeff_e, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# geofac_grdiv', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % geofac_grdiv, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % geofac_grdiv, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % geofac_grdiv, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# geofac_rot', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % geofac_rot, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % geofac_rot, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % geofac_rot, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# geofac_n2s', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % geofac_n2s, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % geofac_n2s, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % geofac_n2s, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_int_state
  SUBROUTINE W_t_grid_cells(io, x, cleanup, nline)
    USE mo_model_domain, ONLY: t_grid_cells
    INTEGER :: io
    TYPE(t_grid_cells), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# neighbor_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % neighbor_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % neighbor_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % neighbor_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# neighbor_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % neighbor_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % neighbor_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % neighbor_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# edge_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % edge_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % edge_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % edge_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# edge_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % edge_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % edge_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % edge_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# area', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % area), cleanup = .FALSE.)
    IF (ASSOCIATED(x % area)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % area, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# start_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % start_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % start_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % start_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# end_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % end_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % end_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % end_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# start_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % start_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % start_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % start_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# end_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % end_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % end_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % end_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# decomp_info', cleanup = .FALSE.)
    CALL serialize(io, x % decomp_info, cleanup = .FALSE.)
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_grid_cells
  SUBROUTINE W_t_grid_edges(io, x, cleanup, nline)
    USE mo_model_domain, ONLY: t_grid_edges
    INTEGER :: io
    TYPE(t_grid_edges), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# cell_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % cell_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % cell_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % cell_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# cell_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % cell_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % cell_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % cell_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# vertex_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % vertex_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % vertex_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % vertex_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# vertex_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % vertex_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % vertex_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % vertex_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# tangent_orientation', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 2, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, SIZE(x % tangent_orientation, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, LBOUND(x % tangent_orientation, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % tangent_orientation, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# quad_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % quad_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % quad_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % quad_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# quad_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % quad_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % quad_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % quad_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# inv_primal_edge_length', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 2, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, SIZE(x % inv_primal_edge_length, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, LBOUND(x % inv_primal_edge_length, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % inv_primal_edge_length, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# inv_dual_edge_length', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 2, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, SIZE(x % inv_dual_edge_length, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, LBOUND(x % inv_dual_edge_length, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % inv_dual_edge_length, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# area_edge', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 2, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, SIZE(x % area_edge, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, LBOUND(x % area_edge, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % area_edge, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# f_e', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 2, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, SIZE(x % f_e, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, LBOUND(x % f_e, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % f_e, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# fn_e', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 2, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, SIZE(x % fn_e, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, LBOUND(x % fn_e, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % fn_e, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# ft_e', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 2, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, SIZE(x % ft_e, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 2
      CALL serialize(io, LBOUND(x % ft_e, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % ft_e, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# start_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % start_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % start_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % start_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# end_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % end_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % end_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % end_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# start_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % start_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % start_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % start_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# end_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % end_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % end_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % end_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_grid_edges
  SUBROUTINE W_t_grid_vertices(io, x, cleanup, nline)
    USE mo_model_domain, ONLY: t_grid_vertices
    INTEGER :: io
    TYPE(t_grid_vertices), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# cell_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % cell_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % cell_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % cell_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# cell_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % cell_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % cell_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % cell_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# edge_idx', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % edge_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % edge_idx, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % edge_idx, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# edge_blk', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 3, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, SIZE(x % edge_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 3
      CALL serialize(io, LBOUND(x % edge_blk, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % edge_blk, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# start_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % start_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % start_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % start_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# end_index', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % end_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % end_index, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % end_index, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# start_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % start_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % start_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % start_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    CALL serialize(io, '# end_block', cleanup = .FALSE.)
    CALL serialize(io, '# alloc', cleanup = .FALSE.)
    CALL serialize(io, .TRUE., cleanup = .FALSE.)
    CALL serialize(io, "# rank", cleanup = .FALSE.)
    CALL serialize(io, 1, cleanup = .FALSE.)
    CALL serialize(io, "# size", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, SIZE(x % end_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, "# lbound", cleanup = .FALSE.)
    DO kmeta = 1, 1
      CALL serialize(io, LBOUND(x % end_block, kmeta), cleanup = .FALSE.)
    END DO
    CALL serialize(io, x % end_block, cleanup = .FALSE., nline = .TRUE., meta = .FALSE.)
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_grid_vertices
  SUBROUTINE W_t_patch(io, x, cleanup, nline)
    USE mo_model_domain, ONLY: t_patch
    INTEGER :: io
    TYPE(t_patch), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# nblks_c', cleanup = .FALSE.)
    CALL serialize(io, x % nblks_c, cleanup = .FALSE.)
    CALL serialize(io, '# nblks_e', cleanup = .FALSE.)
    CALL serialize(io, x % nblks_e, cleanup = .FALSE.)
    CALL serialize(io, '# nblks_v', cleanup = .FALSE.)
    CALL serialize(io, x % nblks_v, cleanup = .FALSE.)
    CALL serialize(io, '# cells', cleanup = .FALSE.)
    CALL serialize(io, x % cells, cleanup = .FALSE.)
    CALL serialize(io, '# edges', cleanup = .FALSE.)
    CALL serialize(io, x % edges, cleanup = .FALSE.)
    CALL serialize(io, '# verts', cleanup = .FALSE.)
    CALL serialize(io, x % verts, cleanup = .FALSE.)
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_patch
  SUBROUTINE W_t_nh_prog(io, x, cleanup, nline)
    USE mo_nonhydro_types, ONLY: t_nh_prog
    INTEGER :: io
    TYPE(t_nh_prog), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# w', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % w), cleanup = .FALSE.)
    IF (ASSOCIATED(x % w)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % w, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# vn', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % vn), cleanup = .FALSE.)
    IF (ASSOCIATED(x % vn)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % vn, cleanup = .FALSE.)
    END IF
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_nh_prog
  SUBROUTINE W_t_nh_diag(io, x, cleanup, nline)
    USE mo_nonhydro_types, ONLY: t_nh_diag
    INTEGER :: io
    TYPE(t_nh_diag), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# vt', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % vt), cleanup = .FALSE.)
    IF (ASSOCIATED(x % vt)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % vt, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# vn_ie', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % vn_ie), cleanup = .FALSE.)
    IF (ASSOCIATED(x % vn_ie)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % vn_ie, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# w_concorr_c', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % w_concorr_c), cleanup = .FALSE.)
    IF (ASSOCIATED(x % w_concorr_c)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % w_concorr_c, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# ddt_vn_apc_pc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddt_vn_apc_pc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddt_vn_apc_pc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddt_vn_apc_pc, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# ddt_w_adv_pc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddt_w_adv_pc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddt_w_adv_pc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddt_w_adv_pc, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# max_vcfl_dyn', cleanup = .FALSE.)
    CALL serialize(io, x % max_vcfl_dyn, cleanup = .FALSE.)
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_nh_diag
  SUBROUTINE W_t_nh_metrics(io, x, cleanup, nline)
    USE mo_nonhydro_types, ONLY: t_nh_metrics
    INTEGER :: io
    TYPE(t_nh_metrics), TARGET, INTENT(IN) :: x
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline
    INTEGER :: kmeta, kmeta_0, kmeta_1, kmeta_2, kmeta_3, kmeta_4, kmeta_5, kmeta_6, kmeta_7, kmeta_8, kmeta_9
    LOGICAL :: cleanup_local, nline_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    CALL serialize(io, '# ddxn_z_full', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddxn_z_full), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddxn_z_full)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddxn_z_full, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# ddxt_z_full', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddxt_z_full), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddxt_z_full)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddxt_z_full, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# ddqz_z_full_e', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddqz_z_full_e), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddqz_z_full_e)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddqz_z_full_e, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# ddqz_z_half', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % ddqz_z_half), cleanup = .FALSE.)
    IF (ASSOCIATED(x % ddqz_z_half)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % ddqz_z_half, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# wgtfac_c', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % wgtfac_c), cleanup = .FALSE.)
    IF (ASSOCIATED(x % wgtfac_c)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % wgtfac_c, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# wgtfac_e', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % wgtfac_e), cleanup = .FALSE.)
    IF (ASSOCIATED(x % wgtfac_e)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % wgtfac_e, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# wgtfacq_e', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % wgtfacq_e), cleanup = .FALSE.)
    IF (ASSOCIATED(x % wgtfacq_e)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % wgtfacq_e, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# coeff_gradekin', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % coeff_gradekin), cleanup = .FALSE.)
    IF (ASSOCIATED(x % coeff_gradekin)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % coeff_gradekin, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# coeff1_dwdz', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % coeff1_dwdz), cleanup = .FALSE.)
    IF (ASSOCIATED(x % coeff1_dwdz)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % coeff1_dwdz, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# coeff2_dwdz', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % coeff2_dwdz), cleanup = .FALSE.)
    IF (ASSOCIATED(x % coeff2_dwdz)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % coeff2_dwdz, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# deepatmo_gradh_mc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % deepatmo_gradh_mc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % deepatmo_gradh_mc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % deepatmo_gradh_mc, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# deepatmo_invr_mc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % deepatmo_invr_mc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % deepatmo_invr_mc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % deepatmo_invr_mc, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# deepatmo_gradh_ifc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % deepatmo_gradh_ifc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % deepatmo_gradh_ifc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % deepatmo_gradh_ifc, cleanup = .FALSE.)
    END IF
    CALL serialize(io, '# deepatmo_invr_ifc', cleanup = .FALSE.)
    CALL serialize(io, '# assoc', cleanup = .FALSE.)
    CALL serialize(io, ASSOCIATED(x % deepatmo_invr_ifc), cleanup = .FALSE.)
    IF (ASSOCIATED(x % deepatmo_invr_ifc)) THEN
      kmeta = 0
      CALL serialize(io, "# missing", cleanup = .FALSE.)
      CALL serialize(io, (kmeta == 0), cleanup = .FALSE.)
      CALL serialize(io, x % deepatmo_invr_ifc, cleanup = .FALSE.)
    END IF
    IF (nline_local) WRITE(io, '(g0)', ADVANCE = 'no') NEW_LINE('A')
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_t_nh_metrics
  SUBROUTINE W_logical_R_2(io, x, cleanup, nline, meta)
    INTEGER :: io
    LOGICAL, INTENT(IN) :: x(:, :)
    INTEGER :: k, kmeta, k1, k2
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k2 = LBOUND(x, 2), UBOUND(x, 2)
      DO k1 = LBOUND(x, 1), UBOUND(x, 1)
        CALL serialize(io, x(k1, k2), cleanup = .FALSE.)
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_logical_R_2
  SUBROUTINE write_real8_bulk(io, arr, n)
    INTEGER, INTENT(IN) :: io, n
    REAL(KIND = 8), INTENT(IN) :: arr(n)
    INTEGER, PARAMETER :: CHUNK = 4096
    INTEGER, PARAMETER :: LINELEN = 29  ! "e28.20" + newline
    CHARACTER(LEN = CHUNK * LINELEN) :: buf
    CHARACTER(LEN = 50) :: tmp
    INTEGER :: i, pos, blk_end
    i = 1
    DO WHILE (i <= n)
      blk_end = MIN(i + CHUNK - 1, n)
      pos = 1
      DO WHILE (i <= blk_end)
        WRITE(tmp, '(e28.20)') arr(i)
        buf(pos:pos+27) = ADJUSTL(tmp(1:28))
        ! Find actual length after ADJUSTL and place newline
        pos = pos + LEN_TRIM(buf(pos:pos+27))
        buf(pos:pos) = NEW_LINE('A')
        pos = pos + 1
        i = i + 1
      END DO
      WRITE(io, '(A)', ADVANCE = 'no') buf(1:pos-1)
    END DO
  END SUBROUTINE write_real8_bulk

  SUBROUTINE W_real__8_R_3(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 8), INTENT(IN) :: x(:, :, :)
    INTEGER :: k, kmeta, k1, k2, k3
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    REAL(KIND = 8), ALLOCATABLE :: flat(:)
    INTEGER :: ntotal, idx
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    ntotal = SIZE(x)
    ALLOCATE(flat(ntotal))
    idx = 0
    DO k3 = LBOUND(x, 3), UBOUND(x, 3)
      DO k2 = LBOUND(x, 2), UBOUND(x, 2)
        DO k1 = LBOUND(x, 1), UBOUND(x, 1)
          idx = idx + 1
          flat(idx) = x(k1, k2, k3)
        END DO
      END DO
    END DO
    CALL write_real8_bulk(io, flat, ntotal)
    DEALLOCATE(flat)
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__8_R_3
  SUBROUTINE W_integer__4_R_3(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 4), INTENT(IN) :: x(:, :, :)
    INTEGER :: k, kmeta, k1, k2, k3
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 3, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 3
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k3 = LBOUND(x, 3), UBOUND(x, 3)
      DO k2 = LBOUND(x, 2), UBOUND(x, 2)
        DO k1 = LBOUND(x, 1), UBOUND(x, 1)
          CALL serialize(io, x(k1, k2, k3), cleanup = .FALSE.)
        END DO
      END DO
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__4_R_3
  SUBROUTINE W_real__8_R_2(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 8), INTENT(IN) :: x(:, :)
    INTEGER :: k, kmeta, k1, k2
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    REAL(KIND = 8), ALLOCATABLE :: flat(:)
    INTEGER :: ntotal, idx
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 2, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 2
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    ntotal = SIZE(x)
    ALLOCATE(flat(ntotal))
    idx = 0
    DO k2 = LBOUND(x, 2), UBOUND(x, 2)
      DO k1 = LBOUND(x, 1), UBOUND(x, 1)
        idx = idx + 1
        flat(idx) = x(k1, k2)
      END DO
    END DO
    CALL write_real8_bulk(io, flat, ntotal)
    DEALLOCATE(flat)
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__8_R_2
  SUBROUTINE W_integer__4_R_1(io, x, cleanup, nline, meta)
    INTEGER :: io
    INTEGER(KIND = 4), INTENT(IN) :: x(:)
    INTEGER :: k, kmeta, k1
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    DO k1 = LBOUND(x, 1), UBOUND(x, 1)
      CALL serialize(io, x(k1), cleanup = .FALSE.)
    END DO
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_integer__4_R_1
  SUBROUTINE W_real__8_R_4(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 8), INTENT(IN) :: x(:, :, :, :)
    INTEGER :: k, kmeta, k1, k2, k3, k4
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    REAL(KIND = 8), ALLOCATABLE :: flat(:)
    INTEGER :: ntotal, idx
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 4, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 4
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    ntotal = SIZE(x)
    ALLOCATE(flat(ntotal))
    idx = 0
    DO k4 = LBOUND(x, 4), UBOUND(x, 4)
      DO k3 = LBOUND(x, 3), UBOUND(x, 3)
        DO k2 = LBOUND(x, 2), UBOUND(x, 2)
          DO k1 = LBOUND(x, 1), UBOUND(x, 1)
            idx = idx + 1
            flat(idx) = x(k1, k2, k3, k4)
          END DO
        END DO
      END DO
    END DO
    CALL write_real8_bulk(io, flat, ntotal)
    DEALLOCATE(flat)
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__8_R_4
  SUBROUTINE W_real__8_R_1(io, x, cleanup, nline, meta)
    INTEGER :: io
    REAL(KIND = 8), INTENT(IN) :: x(:)
    INTEGER :: k, kmeta, k1
    LOGICAL, OPTIONAL, INTENT(IN) :: cleanup, nline, meta
    LOGICAL :: cleanup_local, nline_local, meta_local
    cleanup_local = .TRUE.
    nline_local = .TRUE.
    meta_local = .TRUE.
    IF (PRESENT(cleanup)) cleanup_local = cleanup
    IF (PRESENT(nline)) nline_local = nline
    IF (PRESENT(meta)) meta_local = meta
    IF (meta_local) THEN
      CALL serialize(io, "# rank", cleanup = .FALSE.)
      CALL serialize(io, 1, cleanup = .FALSE.)
      CALL serialize(io, "# size", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, SIZE(x, kmeta), cleanup = .FALSE.)
      END DO
      CALL serialize(io, "# lbound", cleanup = .FALSE.)
      DO kmeta = 1, 1
        CALL serialize(io, LBOUND(x, kmeta), cleanup = .FALSE.)
      END DO
    END IF
    CALL serialize(io, "# entries", cleanup = .FALSE.)
    CALL write_real8_bulk(io, x, SIZE(x))
    IF (cleanup_local) CLOSE(UNIT = io)
  END SUBROUTINE W_real__8_R_1
  SUBROUTINE serialize_global_data(io)
    USE mo_init_vgrid, ONLY: nflatlev => nflatlev
    USE mo_mpi, ONLY: i_am_accel_node => i_am_accel_node
    USE mo_nonhydrostatic_config, ONLY: lextra_diffu => lextra_diffu
    USE mo_parallel_config, ONLY: nproma => nproma
    USE mo_run_config, ONLY: timers_level => timers_level
    USE mo_timer, ONLY: timer_intp => timer_intp
    USE mo_timer, ONLY: timer_solve_nh_veltend => timer_solve_nh_veltend
    USE mo_vertical_grid, ONLY: nrdmax => nrdmax
    INTEGER :: io
    CALL serialize(io, '# nflatlev', cleanup = .FALSE.)
    CALL serialize(io, nflatlev, cleanup = .FALSE.)
    CALL serialize(io, '# i_am_accel_node', cleanup = .FALSE.)
    CALL serialize(io, i_am_accel_node, cleanup = .FALSE.)
    CALL serialize(io, '# lextra_diffu', cleanup = .FALSE.)
    CALL serialize(io, lextra_diffu, cleanup = .FALSE.)
    CALL serialize(io, '# nproma', cleanup = .FALSE.)
    CALL serialize(io, nproma, cleanup = .FALSE.)
    CALL serialize(io, '# timers_level', cleanup = .FALSE.)
    CALL serialize(io, timers_level, cleanup = .FALSE.)
    CALL serialize(io, '# timer_solve_nh_veltend', cleanup = .FALSE.)
    CALL serialize(io, timer_solve_nh_veltend, cleanup = .FALSE.)
    CALL serialize(io, '# timer_intp', cleanup = .FALSE.)
    CALL serialize(io, timer_intp, cleanup = .FALSE.)
    CALL serialize(io, '# nrdmax', cleanup = .FALSE.)
    CALL serialize(io, nrdmax, cleanup = .FALSE.)
    CLOSE(UNIT = io)
  END SUBROUTINE serialize_global_data
END MODULE vt_serde