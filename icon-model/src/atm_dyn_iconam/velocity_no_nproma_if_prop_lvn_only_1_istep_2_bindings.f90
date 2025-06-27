! Auto-generated file by "/data/ben/spcl/icon-dace/upstream-repo/./icon-model/sdfgs/utils/dace-genfi.py"
module mo_velocity_no_nproma_if_prop_lvn_only_1_istep_2_bindings

  use iso_c_binding

  use mo_nonhydro_types, only: &
    t_nh_diag, &
    t_nh_metrics, &
    t_nh_prog
  use mo_intp_data_strc, only: &
    t_int_state
  use mo_model_domain, only: &
    t_patch, &
    t_grid_cells, &
    t_grid_edges, &
    t_grid_vertices
  use mo_decomposition_tools, only: &
    t_grid_domain_decomp_info
  use mo_init_vgrid, only: &
    nflatlev
  use mo_mpi, only: &
    i_am_accel_node
  use mo_nonhydrostatic_config, only: &
    lextra_diffu
  use mo_parallel_config, only: &
    nproma
  use mo_run_config, only: &
    timers_level
  use mo_timer, only: &
    timer_solve_nh_veltend, &
    timer_intp
  use mo_vertical_grid, only: &
    nrdmax

  implicit none

  private
  public :: run_velocity_no_nproma_if_prop_lvn_only_1_istep_2
  public :: run_velocity_no_nproma_if_prop_lvn_only_1_istep_2_verification
  public :: verify_velocity_no_nproma_if_prop_lvn_only_1_istep_2
  public :: dace_init_velocity_no_nproma_if_prop_lvn_only_1_istep_2
  public :: dace_exit_velocity_no_nproma_if_prop_lvn_only_1_istep_2
  public :: dace_program_velocity_no_nproma_if_prop_lvn_only_1_istep_2


  type, bind(c) :: dace_global_data_type
    integer(kind=c_int) :: i_am_accel_node
    integer(kind=c_int) :: lextra_diffu
    type(c_ptr) :: nflatlev
    integer(kind=c_int) :: nproma
    type(c_ptr) :: nrdmax
    integer(kind=c_int) :: timer_intp
    integer(kind=c_int) :: timer_solve_nh_veltend
    integer(kind=c_int) :: timers_level

  end type dace_global_data_type

  type, bind(c) :: dace_t_nh_diag
    integer(kind=c_int) :: f2dace_SA_ddt_vn_apc_pc_d_0_s_300
    integer(kind=c_int) :: f2dace_SA_ddt_vn_apc_pc_d_1_s_301
    integer(kind=c_int) :: f2dace_SA_ddt_vn_apc_pc_d_2_s_302
    integer(kind=c_int) :: f2dace_SA_ddt_vn_apc_pc_d_3_s_303
    integer(kind=c_int) :: f2dace_SA_ddt_w_adv_pc_d_0_s_304
    integer(kind=c_int) :: f2dace_SA_ddt_w_adv_pc_d_1_s_305
    integer(kind=c_int) :: f2dace_SA_ddt_w_adv_pc_d_2_s_306
    integer(kind=c_int) :: f2dace_SA_ddt_w_adv_pc_d_3_s_307
    integer(kind=c_int) :: f2dace_SA_vn_ie_d_0_s_294
    integer(kind=c_int) :: f2dace_SA_vn_ie_d_1_s_295
    integer(kind=c_int) :: f2dace_SA_vn_ie_d_2_s_296
    integer(kind=c_int) :: f2dace_SA_vt_d_0_s_291
    integer(kind=c_int) :: f2dace_SA_vt_d_1_s_292
    integer(kind=c_int) :: f2dace_SA_vt_d_2_s_293
    integer(kind=c_int) :: f2dace_SA_w_concorr_c_d_0_s_297
    integer(kind=c_int) :: f2dace_SA_w_concorr_c_d_1_s_298
    integer(kind=c_int) :: f2dace_SA_w_concorr_c_d_2_s_299
    integer(kind=c_int) :: f2dace_SOA_ddt_vn_apc_pc_d_0_s_300
    integer(kind=c_int) :: f2dace_SOA_ddt_vn_apc_pc_d_1_s_301
    integer(kind=c_int) :: f2dace_SOA_ddt_vn_apc_pc_d_2_s_302
    integer(kind=c_int) :: f2dace_SOA_ddt_vn_apc_pc_d_3_s_303
    integer(kind=c_int) :: f2dace_SOA_ddt_w_adv_pc_d_0_s_304
    integer(kind=c_int) :: f2dace_SOA_ddt_w_adv_pc_d_1_s_305
    integer(kind=c_int) :: f2dace_SOA_ddt_w_adv_pc_d_2_s_306
    integer(kind=c_int) :: f2dace_SOA_ddt_w_adv_pc_d_3_s_307
    integer(kind=c_int) :: f2dace_SOA_vn_ie_d_0_s_294
    integer(kind=c_int) :: f2dace_SOA_vn_ie_d_1_s_295
    integer(kind=c_int) :: f2dace_SOA_vn_ie_d_2_s_296
    integer(kind=c_int) :: f2dace_SOA_vt_d_0_s_291
    integer(kind=c_int) :: f2dace_SOA_vt_d_1_s_292
    integer(kind=c_int) :: f2dace_SOA_vt_d_2_s_293
    integer(kind=c_int) :: f2dace_SOA_w_concorr_c_d_0_s_297
    integer(kind=c_int) :: f2dace_SOA_w_concorr_c_d_1_s_298
    integer(kind=c_int) :: f2dace_SOA_w_concorr_c_d_2_s_299
    type(c_ptr) :: ddt_vn_apc_pc
    type(c_ptr) :: ddt_w_adv_pc
    real(kind=c_double) :: max_vcfl_dyn
    type(c_ptr) :: vn_ie
    type(c_ptr) :: vt
    type(c_ptr) :: w_concorr_c

  end type dace_t_nh_diag

  type, bind(c) :: dace_t_int_state
    integer(kind=c_int) :: f2dace_SA_c_lin_e_d_0_s_25
    integer(kind=c_int) :: f2dace_SA_c_lin_e_d_1_s_26
    integer(kind=c_int) :: f2dace_SA_c_lin_e_d_2_s_27
    integer(kind=c_int) :: f2dace_SA_cells_aw_verts_d_0_s_31
    integer(kind=c_int) :: f2dace_SA_cells_aw_verts_d_1_s_32
    integer(kind=c_int) :: f2dace_SA_cells_aw_verts_d_2_s_33
    integer(kind=c_int) :: f2dace_SA_e_bln_c_s_d_0_s_28
    integer(kind=c_int) :: f2dace_SA_e_bln_c_s_d_1_s_29
    integer(kind=c_int) :: f2dace_SA_e_bln_c_s_d_2_s_30
    integer(kind=c_int) :: f2dace_SA_geofac_grdiv_d_0_s_37
    integer(kind=c_int) :: f2dace_SA_geofac_grdiv_d_1_s_38
    integer(kind=c_int) :: f2dace_SA_geofac_grdiv_d_2_s_39
    integer(kind=c_int) :: f2dace_SA_geofac_n2s_d_0_s_43
    integer(kind=c_int) :: f2dace_SA_geofac_n2s_d_1_s_44
    integer(kind=c_int) :: f2dace_SA_geofac_n2s_d_2_s_45
    integer(kind=c_int) :: f2dace_SA_geofac_rot_d_0_s_40
    integer(kind=c_int) :: f2dace_SA_geofac_rot_d_1_s_41
    integer(kind=c_int) :: f2dace_SA_geofac_rot_d_2_s_42
    integer(kind=c_int) :: f2dace_SA_rbf_vec_coeff_e_d_0_s_34
    integer(kind=c_int) :: f2dace_SA_rbf_vec_coeff_e_d_1_s_35
    integer(kind=c_int) :: f2dace_SA_rbf_vec_coeff_e_d_2_s_36
    integer(kind=c_int) :: f2dace_SOA_c_lin_e_d_0_s_25
    integer(kind=c_int) :: f2dace_SOA_c_lin_e_d_1_s_26
    integer(kind=c_int) :: f2dace_SOA_c_lin_e_d_2_s_27
    integer(kind=c_int) :: f2dace_SOA_cells_aw_verts_d_0_s_31
    integer(kind=c_int) :: f2dace_SOA_cells_aw_verts_d_1_s_32
    integer(kind=c_int) :: f2dace_SOA_cells_aw_verts_d_2_s_33
    integer(kind=c_int) :: f2dace_SOA_e_bln_c_s_d_0_s_28
    integer(kind=c_int) :: f2dace_SOA_e_bln_c_s_d_1_s_29
    integer(kind=c_int) :: f2dace_SOA_e_bln_c_s_d_2_s_30
    integer(kind=c_int) :: f2dace_SOA_geofac_grdiv_d_0_s_37
    integer(kind=c_int) :: f2dace_SOA_geofac_grdiv_d_1_s_38
    integer(kind=c_int) :: f2dace_SOA_geofac_grdiv_d_2_s_39
    integer(kind=c_int) :: f2dace_SOA_geofac_n2s_d_0_s_43
    integer(kind=c_int) :: f2dace_SOA_geofac_n2s_d_1_s_44
    integer(kind=c_int) :: f2dace_SOA_geofac_n2s_d_2_s_45
    integer(kind=c_int) :: f2dace_SOA_geofac_rot_d_0_s_40
    integer(kind=c_int) :: f2dace_SOA_geofac_rot_d_1_s_41
    integer(kind=c_int) :: f2dace_SOA_geofac_rot_d_2_s_42
    integer(kind=c_int) :: f2dace_SOA_rbf_vec_coeff_e_d_0_s_34
    integer(kind=c_int) :: f2dace_SOA_rbf_vec_coeff_e_d_1_s_35
    integer(kind=c_int) :: f2dace_SOA_rbf_vec_coeff_e_d_2_s_36
    type(c_ptr) :: c_lin_e
    type(c_ptr) :: cells_aw_verts
    type(c_ptr) :: e_bln_c_s
    type(c_ptr) :: geofac_grdiv
    type(c_ptr) :: geofac_n2s
    type(c_ptr) :: geofac_rot
    type(c_ptr) :: rbf_vec_coeff_e

  end type dace_t_int_state

  type, bind(c) :: dace_t_nh_metrics
    integer(kind=c_int) :: f2dace_SA_coeff1_dwdz_d_0_s_332
    integer(kind=c_int) :: f2dace_SA_coeff1_dwdz_d_1_s_333
    integer(kind=c_int) :: f2dace_SA_coeff1_dwdz_d_2_s_334
    integer(kind=c_int) :: f2dace_SA_coeff2_dwdz_d_0_s_335
    integer(kind=c_int) :: f2dace_SA_coeff2_dwdz_d_1_s_336
    integer(kind=c_int) :: f2dace_SA_coeff2_dwdz_d_2_s_337
    integer(kind=c_int) :: f2dace_SA_coeff_gradekin_d_0_s_329
    integer(kind=c_int) :: f2dace_SA_coeff_gradekin_d_1_s_330
    integer(kind=c_int) :: f2dace_SA_coeff_gradekin_d_2_s_331
    integer(kind=c_int) :: f2dace_SA_ddqz_z_full_e_d_0_s_314
    integer(kind=c_int) :: f2dace_SA_ddqz_z_full_e_d_1_s_315
    integer(kind=c_int) :: f2dace_SA_ddqz_z_full_e_d_2_s_316
    integer(kind=c_int) :: f2dace_SA_ddqz_z_half_d_0_s_317
    integer(kind=c_int) :: f2dace_SA_ddqz_z_half_d_1_s_318
    integer(kind=c_int) :: f2dace_SA_ddqz_z_half_d_2_s_319
    integer(kind=c_int) :: f2dace_SA_ddxn_z_full_d_0_s_308
    integer(kind=c_int) :: f2dace_SA_ddxn_z_full_d_1_s_309
    integer(kind=c_int) :: f2dace_SA_ddxn_z_full_d_2_s_310
    integer(kind=c_int) :: f2dace_SA_ddxt_z_full_d_0_s_311
    integer(kind=c_int) :: f2dace_SA_ddxt_z_full_d_1_s_312
    integer(kind=c_int) :: f2dace_SA_ddxt_z_full_d_2_s_313
    integer(kind=c_int) :: f2dace_SA_deepatmo_gradh_ifc_d_0_s_340
    integer(kind=c_int) :: f2dace_SA_deepatmo_gradh_mc_d_0_s_338
    integer(kind=c_int) :: f2dace_SA_deepatmo_invr_ifc_d_0_s_341
    integer(kind=c_int) :: f2dace_SA_deepatmo_invr_mc_d_0_s_339
    integer(kind=c_int) :: f2dace_SA_wgtfac_c_d_0_s_320
    integer(kind=c_int) :: f2dace_SA_wgtfac_c_d_1_s_321
    integer(kind=c_int) :: f2dace_SA_wgtfac_c_d_2_s_322
    integer(kind=c_int) :: f2dace_SA_wgtfac_e_d_0_s_323
    integer(kind=c_int) :: f2dace_SA_wgtfac_e_d_1_s_324
    integer(kind=c_int) :: f2dace_SA_wgtfac_e_d_2_s_325
    integer(kind=c_int) :: f2dace_SA_wgtfacq_e_d_0_s_326
    integer(kind=c_int) :: f2dace_SA_wgtfacq_e_d_1_s_327
    integer(kind=c_int) :: f2dace_SA_wgtfacq_e_d_2_s_328
    integer(kind=c_int) :: f2dace_SOA_coeff1_dwdz_d_0_s_332
    integer(kind=c_int) :: f2dace_SOA_coeff1_dwdz_d_1_s_333
    integer(kind=c_int) :: f2dace_SOA_coeff1_dwdz_d_2_s_334
    integer(kind=c_int) :: f2dace_SOA_coeff2_dwdz_d_0_s_335
    integer(kind=c_int) :: f2dace_SOA_coeff2_dwdz_d_1_s_336
    integer(kind=c_int) :: f2dace_SOA_coeff2_dwdz_d_2_s_337
    integer(kind=c_int) :: f2dace_SOA_coeff_gradekin_d_0_s_329
    integer(kind=c_int) :: f2dace_SOA_coeff_gradekin_d_1_s_330
    integer(kind=c_int) :: f2dace_SOA_coeff_gradekin_d_2_s_331
    integer(kind=c_int) :: f2dace_SOA_ddqz_z_full_e_d_0_s_314
    integer(kind=c_int) :: f2dace_SOA_ddqz_z_full_e_d_1_s_315
    integer(kind=c_int) :: f2dace_SOA_ddqz_z_full_e_d_2_s_316
    integer(kind=c_int) :: f2dace_SOA_ddqz_z_half_d_0_s_317
    integer(kind=c_int) :: f2dace_SOA_ddqz_z_half_d_1_s_318
    integer(kind=c_int) :: f2dace_SOA_ddqz_z_half_d_2_s_319
    integer(kind=c_int) :: f2dace_SOA_ddxn_z_full_d_0_s_308
    integer(kind=c_int) :: f2dace_SOA_ddxn_z_full_d_1_s_309
    integer(kind=c_int) :: f2dace_SOA_ddxn_z_full_d_2_s_310
    integer(kind=c_int) :: f2dace_SOA_ddxt_z_full_d_0_s_311
    integer(kind=c_int) :: f2dace_SOA_ddxt_z_full_d_1_s_312
    integer(kind=c_int) :: f2dace_SOA_ddxt_z_full_d_2_s_313
    integer(kind=c_int) :: f2dace_SOA_deepatmo_gradh_ifc_d_0_s_340
    integer(kind=c_int) :: f2dace_SOA_deepatmo_gradh_mc_d_0_s_338
    integer(kind=c_int) :: f2dace_SOA_deepatmo_invr_ifc_d_0_s_341
    integer(kind=c_int) :: f2dace_SOA_deepatmo_invr_mc_d_0_s_339
    integer(kind=c_int) :: f2dace_SOA_wgtfac_c_d_0_s_320
    integer(kind=c_int) :: f2dace_SOA_wgtfac_c_d_1_s_321
    integer(kind=c_int) :: f2dace_SOA_wgtfac_c_d_2_s_322
    integer(kind=c_int) :: f2dace_SOA_wgtfac_e_d_0_s_323
    integer(kind=c_int) :: f2dace_SOA_wgtfac_e_d_1_s_324
    integer(kind=c_int) :: f2dace_SOA_wgtfac_e_d_2_s_325
    integer(kind=c_int) :: f2dace_SOA_wgtfacq_e_d_0_s_326
    integer(kind=c_int) :: f2dace_SOA_wgtfacq_e_d_1_s_327
    integer(kind=c_int) :: f2dace_SOA_wgtfacq_e_d_2_s_328
    type(c_ptr) :: coeff1_dwdz
    type(c_ptr) :: coeff2_dwdz
    type(c_ptr) :: coeff_gradekin
    type(c_ptr) :: ddqz_z_full_e
    type(c_ptr) :: ddqz_z_half
    type(c_ptr) :: ddxn_z_full
    type(c_ptr) :: ddxt_z_full
    type(c_ptr) :: deepatmo_gradh_ifc
    type(c_ptr) :: deepatmo_gradh_mc
    type(c_ptr) :: deepatmo_invr_ifc
    type(c_ptr) :: deepatmo_invr_mc
    type(c_ptr) :: wgtfac_c
    type(c_ptr) :: wgtfac_e
    type(c_ptr) :: wgtfacq_e

  end type dace_t_nh_metrics

  type, bind(c) :: dace_t_patch
    type(c_ptr) :: cells
    type(c_ptr) :: edges
    integer(kind=c_int) :: nblks_c
    integer(kind=c_int) :: nblks_e
    integer(kind=c_int) :: nblks_v
    type(c_ptr) :: verts

  end type dace_t_patch

  type, bind(c) :: dace_t_grid_cells
    integer(kind=c_int) :: f2dace_SA_area_d_0_s_158
    integer(kind=c_int) :: f2dace_SA_area_d_1_s_159
    integer(kind=c_int) :: f2dace_SA_edge_blk_d_0_s_155
    integer(kind=c_int) :: f2dace_SA_edge_blk_d_1_s_156
    integer(kind=c_int) :: f2dace_SA_edge_blk_d_2_s_157
    integer(kind=c_int) :: f2dace_SA_edge_idx_d_0_s_152
    integer(kind=c_int) :: f2dace_SA_edge_idx_d_1_s_153
    integer(kind=c_int) :: f2dace_SA_edge_idx_d_2_s_154
    integer(kind=c_int) :: f2dace_SA_end_block_d_0_s_163
    integer(kind=c_int) :: f2dace_SA_end_index_d_0_s_161
    integer(kind=c_int) :: f2dace_SA_neighbor_blk_d_0_s_149
    integer(kind=c_int) :: f2dace_SA_neighbor_blk_d_1_s_150
    integer(kind=c_int) :: f2dace_SA_neighbor_blk_d_2_s_151
    integer(kind=c_int) :: f2dace_SA_neighbor_idx_d_0_s_146
    integer(kind=c_int) :: f2dace_SA_neighbor_idx_d_1_s_147
    integer(kind=c_int) :: f2dace_SA_neighbor_idx_d_2_s_148
    integer(kind=c_int) :: f2dace_SA_start_block_d_0_s_162
    integer(kind=c_int) :: f2dace_SA_start_index_d_0_s_160
    integer(kind=c_int) :: f2dace_SOA_area_d_0_s_158
    integer(kind=c_int) :: f2dace_SOA_area_d_1_s_159
    integer(kind=c_int) :: f2dace_SOA_edge_blk_d_0_s_155
    integer(kind=c_int) :: f2dace_SOA_edge_blk_d_1_s_156
    integer(kind=c_int) :: f2dace_SOA_edge_blk_d_2_s_157
    integer(kind=c_int) :: f2dace_SOA_edge_idx_d_0_s_152
    integer(kind=c_int) :: f2dace_SOA_edge_idx_d_1_s_153
    integer(kind=c_int) :: f2dace_SOA_edge_idx_d_2_s_154
    integer(kind=c_int) :: f2dace_SOA_end_block_d_0_s_163
    integer(kind=c_int) :: f2dace_SOA_end_index_d_0_s_161
    integer(kind=c_int) :: f2dace_SOA_neighbor_blk_d_0_s_149
    integer(kind=c_int) :: f2dace_SOA_neighbor_blk_d_1_s_150
    integer(kind=c_int) :: f2dace_SOA_neighbor_blk_d_2_s_151
    integer(kind=c_int) :: f2dace_SOA_neighbor_idx_d_0_s_146
    integer(kind=c_int) :: f2dace_SOA_neighbor_idx_d_1_s_147
    integer(kind=c_int) :: f2dace_SOA_neighbor_idx_d_2_s_148
    integer(kind=c_int) :: f2dace_SOA_start_block_d_0_s_162
    integer(kind=c_int) :: f2dace_SOA_start_index_d_0_s_160
    type(c_ptr) :: area
    type(c_ptr) :: decomp_info
    type(c_ptr) :: edge_blk
    type(c_ptr) :: edge_idx
    type(c_ptr) :: end_block
    type(c_ptr) :: end_index
    type(c_ptr) :: neighbor_blk
    type(c_ptr) :: neighbor_idx
    type(c_ptr) :: start_block
    type(c_ptr) :: start_index

  end type dace_t_grid_cells

  type, bind(c) :: dace_t_grid_domain_decomp_info
    integer(kind=c_int) :: f2dace_SA_owner_mask_d_0_s_2
    integer(kind=c_int) :: f2dace_SA_owner_mask_d_1_s_3
    integer(kind=c_int) :: f2dace_SOA_owner_mask_d_0_s_2
    integer(kind=c_int) :: f2dace_SOA_owner_mask_d_1_s_3
    type(c_ptr) :: owner_mask

  end type dace_t_grid_domain_decomp_info

  type, bind(c) :: dace_t_grid_edges
    integer(kind=c_int) :: f2dace_SA_area_edge_d_0_s_188
    integer(kind=c_int) :: f2dace_SA_area_edge_d_1_s_189
    integer(kind=c_int) :: f2dace_SA_cell_blk_d_0_s_167
    integer(kind=c_int) :: f2dace_SA_cell_blk_d_1_s_168
    integer(kind=c_int) :: f2dace_SA_cell_blk_d_2_s_169
    integer(kind=c_int) :: f2dace_SA_cell_idx_d_0_s_164
    integer(kind=c_int) :: f2dace_SA_cell_idx_d_1_s_165
    integer(kind=c_int) :: f2dace_SA_cell_idx_d_2_s_166
    integer(kind=c_int) :: f2dace_SA_end_block_d_0_s_199
    integer(kind=c_int) :: f2dace_SA_end_index_d_0_s_197
    integer(kind=c_int) :: f2dace_SA_f_e_d_0_s_190
    integer(kind=c_int) :: f2dace_SA_f_e_d_1_s_191
    integer(kind=c_int) :: f2dace_SA_fn_e_d_0_s_192
    integer(kind=c_int) :: f2dace_SA_fn_e_d_1_s_193
    integer(kind=c_int) :: f2dace_SA_ft_e_d_0_s_194
    integer(kind=c_int) :: f2dace_SA_ft_e_d_1_s_195
    integer(kind=c_int) :: f2dace_SA_inv_dual_edge_length_d_0_s_186
    integer(kind=c_int) :: f2dace_SA_inv_dual_edge_length_d_1_s_187
    integer(kind=c_int) :: f2dace_SA_inv_primal_edge_length_d_0_s_184
    integer(kind=c_int) :: f2dace_SA_inv_primal_edge_length_d_1_s_185
    integer(kind=c_int) :: f2dace_SA_quad_blk_d_0_s_181
    integer(kind=c_int) :: f2dace_SA_quad_blk_d_1_s_182
    integer(kind=c_int) :: f2dace_SA_quad_blk_d_2_s_183
    integer(kind=c_int) :: f2dace_SA_quad_idx_d_0_s_178
    integer(kind=c_int) :: f2dace_SA_quad_idx_d_1_s_179
    integer(kind=c_int) :: f2dace_SA_quad_idx_d_2_s_180
    integer(kind=c_int) :: f2dace_SA_start_block_d_0_s_198
    integer(kind=c_int) :: f2dace_SA_start_index_d_0_s_196
    integer(kind=c_int) :: f2dace_SA_tangent_orientation_d_0_s_176
    integer(kind=c_int) :: f2dace_SA_tangent_orientation_d_1_s_177
    integer(kind=c_int) :: f2dace_SA_vertex_blk_d_0_s_173
    integer(kind=c_int) :: f2dace_SA_vertex_blk_d_1_s_174
    integer(kind=c_int) :: f2dace_SA_vertex_blk_d_2_s_175
    integer(kind=c_int) :: f2dace_SA_vertex_idx_d_0_s_170
    integer(kind=c_int) :: f2dace_SA_vertex_idx_d_1_s_171
    integer(kind=c_int) :: f2dace_SA_vertex_idx_d_2_s_172
    integer(kind=c_int) :: f2dace_SOA_area_edge_d_0_s_188
    integer(kind=c_int) :: f2dace_SOA_area_edge_d_1_s_189
    integer(kind=c_int) :: f2dace_SOA_cell_blk_d_0_s_167
    integer(kind=c_int) :: f2dace_SOA_cell_blk_d_1_s_168
    integer(kind=c_int) :: f2dace_SOA_cell_blk_d_2_s_169
    integer(kind=c_int) :: f2dace_SOA_cell_idx_d_0_s_164
    integer(kind=c_int) :: f2dace_SOA_cell_idx_d_1_s_165
    integer(kind=c_int) :: f2dace_SOA_cell_idx_d_2_s_166
    integer(kind=c_int) :: f2dace_SOA_end_block_d_0_s_199
    integer(kind=c_int) :: f2dace_SOA_end_index_d_0_s_197
    integer(kind=c_int) :: f2dace_SOA_f_e_d_0_s_190
    integer(kind=c_int) :: f2dace_SOA_f_e_d_1_s_191
    integer(kind=c_int) :: f2dace_SOA_fn_e_d_0_s_192
    integer(kind=c_int) :: f2dace_SOA_fn_e_d_1_s_193
    integer(kind=c_int) :: f2dace_SOA_ft_e_d_0_s_194
    integer(kind=c_int) :: f2dace_SOA_ft_e_d_1_s_195
    integer(kind=c_int) :: f2dace_SOA_inv_dual_edge_length_d_0_s_186
    integer(kind=c_int) :: f2dace_SOA_inv_dual_edge_length_d_1_s_187
    integer(kind=c_int) :: f2dace_SOA_inv_primal_edge_length_d_0_s_184
    integer(kind=c_int) :: f2dace_SOA_inv_primal_edge_length_d_1_s_185
    integer(kind=c_int) :: f2dace_SOA_quad_blk_d_0_s_181
    integer(kind=c_int) :: f2dace_SOA_quad_blk_d_1_s_182
    integer(kind=c_int) :: f2dace_SOA_quad_blk_d_2_s_183
    integer(kind=c_int) :: f2dace_SOA_quad_idx_d_0_s_178
    integer(kind=c_int) :: f2dace_SOA_quad_idx_d_1_s_179
    integer(kind=c_int) :: f2dace_SOA_quad_idx_d_2_s_180
    integer(kind=c_int) :: f2dace_SOA_start_block_d_0_s_198
    integer(kind=c_int) :: f2dace_SOA_start_index_d_0_s_196
    integer(kind=c_int) :: f2dace_SOA_tangent_orientation_d_0_s_176
    integer(kind=c_int) :: f2dace_SOA_tangent_orientation_d_1_s_177
    integer(kind=c_int) :: f2dace_SOA_vertex_blk_d_0_s_173
    integer(kind=c_int) :: f2dace_SOA_vertex_blk_d_1_s_174
    integer(kind=c_int) :: f2dace_SOA_vertex_blk_d_2_s_175
    integer(kind=c_int) :: f2dace_SOA_vertex_idx_d_0_s_170
    integer(kind=c_int) :: f2dace_SOA_vertex_idx_d_1_s_171
    integer(kind=c_int) :: f2dace_SOA_vertex_idx_d_2_s_172
    type(c_ptr) :: area_edge
    type(c_ptr) :: cell_blk
    type(c_ptr) :: cell_idx
    type(c_ptr) :: end_block
    type(c_ptr) :: end_index
    type(c_ptr) :: f_e
    type(c_ptr) :: fn_e
    type(c_ptr) :: ft_e
    type(c_ptr) :: inv_dual_edge_length
    type(c_ptr) :: inv_primal_edge_length
    type(c_ptr) :: quad_blk
    type(c_ptr) :: quad_idx
    type(c_ptr) :: start_block
    type(c_ptr) :: start_index
    type(c_ptr) :: tangent_orientation
    type(c_ptr) :: vertex_blk
    type(c_ptr) :: vertex_idx

  end type dace_t_grid_edges

  type, bind(c) :: dace_t_grid_vertices
    integer(kind=c_int) :: f2dace_SA_cell_blk_d_0_s_203
    integer(kind=c_int) :: f2dace_SA_cell_blk_d_1_s_204
    integer(kind=c_int) :: f2dace_SA_cell_blk_d_2_s_205
    integer(kind=c_int) :: f2dace_SA_cell_idx_d_0_s_200
    integer(kind=c_int) :: f2dace_SA_cell_idx_d_1_s_201
    integer(kind=c_int) :: f2dace_SA_cell_idx_d_2_s_202
    integer(kind=c_int) :: f2dace_SA_edge_blk_d_0_s_209
    integer(kind=c_int) :: f2dace_SA_edge_blk_d_1_s_210
    integer(kind=c_int) :: f2dace_SA_edge_blk_d_2_s_211
    integer(kind=c_int) :: f2dace_SA_edge_idx_d_0_s_206
    integer(kind=c_int) :: f2dace_SA_edge_idx_d_1_s_207
    integer(kind=c_int) :: f2dace_SA_edge_idx_d_2_s_208
    integer(kind=c_int) :: f2dace_SA_end_block_d_0_s_215
    integer(kind=c_int) :: f2dace_SA_end_index_d_0_s_213
    integer(kind=c_int) :: f2dace_SA_start_block_d_0_s_214
    integer(kind=c_int) :: f2dace_SA_start_index_d_0_s_212
    integer(kind=c_int) :: f2dace_SOA_cell_blk_d_0_s_203
    integer(kind=c_int) :: f2dace_SOA_cell_blk_d_1_s_204
    integer(kind=c_int) :: f2dace_SOA_cell_blk_d_2_s_205
    integer(kind=c_int) :: f2dace_SOA_cell_idx_d_0_s_200
    integer(kind=c_int) :: f2dace_SOA_cell_idx_d_1_s_201
    integer(kind=c_int) :: f2dace_SOA_cell_idx_d_2_s_202
    integer(kind=c_int) :: f2dace_SOA_edge_blk_d_0_s_209
    integer(kind=c_int) :: f2dace_SOA_edge_blk_d_1_s_210
    integer(kind=c_int) :: f2dace_SOA_edge_blk_d_2_s_211
    integer(kind=c_int) :: f2dace_SOA_edge_idx_d_0_s_206
    integer(kind=c_int) :: f2dace_SOA_edge_idx_d_1_s_207
    integer(kind=c_int) :: f2dace_SOA_edge_idx_d_2_s_208
    integer(kind=c_int) :: f2dace_SOA_end_block_d_0_s_215
    integer(kind=c_int) :: f2dace_SOA_end_index_d_0_s_213
    integer(kind=c_int) :: f2dace_SOA_start_block_d_0_s_214
    integer(kind=c_int) :: f2dace_SOA_start_index_d_0_s_212
    type(c_ptr) :: cell_blk
    type(c_ptr) :: cell_idx
    type(c_ptr) :: edge_blk
    type(c_ptr) :: edge_idx
    type(c_ptr) :: end_block
    type(c_ptr) :: end_index
    type(c_ptr) :: start_block
    type(c_ptr) :: start_index

  end type dace_t_grid_vertices

  type, bind(c) :: dace_t_nh_prog
    integer(kind=c_int) :: f2dace_SA_vn_d_0_s_288
    integer(kind=c_int) :: f2dace_SA_vn_d_1_s_289
    integer(kind=c_int) :: f2dace_SA_vn_d_2_s_290
    integer(kind=c_int) :: f2dace_SA_w_d_0_s_285
    integer(kind=c_int) :: f2dace_SA_w_d_1_s_286
    integer(kind=c_int) :: f2dace_SA_w_d_2_s_287
    integer(kind=c_int) :: f2dace_SOA_vn_d_0_s_288
    integer(kind=c_int) :: f2dace_SOA_vn_d_1_s_289
    integer(kind=c_int) :: f2dace_SOA_vn_d_2_s_290
    integer(kind=c_int) :: f2dace_SOA_w_d_0_s_285
    integer(kind=c_int) :: f2dace_SOA_w_d_1_s_286
    integer(kind=c_int) :: f2dace_SOA_w_d_2_s_287
    type(c_ptr) :: vn
    type(c_ptr) :: w

  end type dace_t_nh_prog


  logical :: is_initialized = .false.
  type(c_ptr) :: dace_state = C_NULL_PTR

  type(c_ptr) :: cached_shallow_copy_global_data = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_p_diag = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_p_int = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_p_metrics = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_p_patch = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_p_prog = C_NULL_PTR

  type(c_ptr) :: copy_or_ptr_global_data = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_p_diag = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_p_int = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_p_metrics = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_p_patch = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_p_prog = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_kin_hor_e = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_vt_ie = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_w_concorr_me = C_NULL_PTR

interface

  type(c_ptr) function malloc(size) &
    bind(c, name="malloc")
    use iso_c_binding
    integer(kind=c_size_t), value :: size
  end function malloc

  subroutine free(ptr) &
    bind(c, name="free")
    use iso_c_binding
    type(c_ptr), value :: ptr
  end subroutine free

#ifdef _OPENACC

  type(c_ptr) function c_acc_malloc(size) &
    bind(c, name="acc_malloc")
    use iso_c_binding
    integer(kind=c_size_t), value :: size
  end function c_acc_malloc

  subroutine c_acc_free(device_ptr) &
    bind(c, name="acc_free")
    use iso_c_binding
    type(c_ptr), value :: device_ptr
  end subroutine c_acc_free

  type(c_ptr) function c_acc_deviceptr(host_ptr) &
    bind(c, name="acc_deviceptr")
    use iso_c_binding
    type(c_ptr), value :: host_ptr
  end function c_acc_deviceptr

  subroutine c_acc_memcpy_device(dst, src, size) &
    bind(c, name="acc_memcpy_device")
    use iso_c_binding
    type(c_ptr), value :: dst
    type(c_ptr), value :: src
    integer(kind=c_size_t), value :: size
  end subroutine c_acc_memcpy_device

#endif


  type(c_ptr) function dace_init_velocity_no_nproma_if_prop_lvn_only_1_istep_2( &
    global_data, &
    p_diag, &
    p_int, &
    p_metrics, &
    p_patch, &
    p_prog, &
    z_kin_hor_e, &
    z_vt_ie, &
    z_w_concorr_me, &
    f2dace_A_z_kin_hor_e_d_0_s_363, &
    f2dace_A_z_kin_hor_e_d_1_s_364, &
    f2dace_A_z_kin_hor_e_d_2_s_365, &
    f2dace_A_z_vt_ie_d_0_s_366, &
    f2dace_A_z_vt_ie_d_1_s_367, &
    f2dace_A_z_vt_ie_d_2_s_368, &
    f2dace_A_z_w_concorr_me_d_0_s_360, &
    f2dace_A_z_w_concorr_me_d_1_s_361, &
    f2dace_A_z_w_concorr_me_d_2_s_362, &
    f2dace_OA_z_kin_hor_e_d_0_s_363, &
    f2dace_OA_z_kin_hor_e_d_1_s_364, &
    f2dace_OA_z_kin_hor_e_d_2_s_365, &
    dt_linintp_ubc, &
    dtime, &
    istep, &
    ldeepatmo, &
    lvn_only, &
    ntnd &
  ) &
    bind(c, name="__dace_init_velocity_no_nproma_if_prop_lvn_only_1_istep_2")
    use iso_c_binding

    type(c_ptr), value :: global_data
    type(c_ptr), value :: p_diag
    type(c_ptr), value :: p_int
    type(c_ptr), value :: p_metrics
    type(c_ptr), value :: p_patch
    type(c_ptr), value :: p_prog
    type(c_ptr), value :: z_kin_hor_e
    type(c_ptr), value :: z_vt_ie
    type(c_ptr), value :: z_w_concorr_me
    integer(kind=c_int), value :: f2dace_A_z_kin_hor_e_d_0_s_363
    integer(kind=c_int), value :: f2dace_A_z_kin_hor_e_d_1_s_364
    integer(kind=c_int), value :: f2dace_A_z_kin_hor_e_d_2_s_365
    integer(kind=c_int), value :: f2dace_A_z_vt_ie_d_0_s_366
    integer(kind=c_int), value :: f2dace_A_z_vt_ie_d_1_s_367
    integer(kind=c_int), value :: f2dace_A_z_vt_ie_d_2_s_368
    integer(kind=c_int), value :: f2dace_A_z_w_concorr_me_d_0_s_360
    integer(kind=c_int), value :: f2dace_A_z_w_concorr_me_d_1_s_361
    integer(kind=c_int), value :: f2dace_A_z_w_concorr_me_d_2_s_362
    integer(kind=c_int), value :: f2dace_OA_z_kin_hor_e_d_0_s_363
    integer(kind=c_int), value :: f2dace_OA_z_kin_hor_e_d_1_s_364
    integer(kind=c_int), value :: f2dace_OA_z_kin_hor_e_d_2_s_365
    real(kind=c_double), value :: dt_linintp_ubc
    real(kind=c_double), value :: dtime
    integer(kind=c_int), value :: istep
    integer(kind=c_int), value :: ldeepatmo
    integer(kind=c_int), value :: lvn_only
    integer(kind=c_int), value :: ntnd
  end function dace_init_velocity_no_nproma_if_prop_lvn_only_1_istep_2

  integer(c_int) function dace_exit_velocity_no_nproma_if_prop_lvn_only_1_istep_2(state) &
    bind(c, name="__dace_exit_velocity_no_nproma_if_prop_lvn_only_1_istep_2")
    use iso_c_binding

    type(c_ptr), value :: state
  end function dace_exit_velocity_no_nproma_if_prop_lvn_only_1_istep_2

  subroutine dace_program_velocity_no_nproma_if_prop_lvn_only_1_istep_2( &
    state, &
    global_data, &
    p_diag, &
    p_int, &
    p_metrics, &
    p_patch, &
    p_prog, &
    z_kin_hor_e, &
    z_vt_ie, &
    z_w_concorr_me, &
    f2dace_A_z_kin_hor_e_d_0_s_363, &
    f2dace_A_z_kin_hor_e_d_1_s_364, &
    f2dace_A_z_kin_hor_e_d_2_s_365, &
    f2dace_A_z_vt_ie_d_0_s_366, &
    f2dace_A_z_vt_ie_d_1_s_367, &
    f2dace_A_z_vt_ie_d_2_s_368, &
    f2dace_A_z_w_concorr_me_d_0_s_360, &
    f2dace_A_z_w_concorr_me_d_1_s_361, &
    f2dace_A_z_w_concorr_me_d_2_s_362, &
    f2dace_OA_z_kin_hor_e_d_0_s_363, &
    f2dace_OA_z_kin_hor_e_d_1_s_364, &
    f2dace_OA_z_kin_hor_e_d_2_s_365, &
    dt_linintp_ubc, &
    dtime, &
    istep, &
    ldeepatmo, &
    lvn_only, &
    ntnd &
  ) &
    bind(c, name="__program_velocity_no_nproma_if_prop_lvn_only_1_istep_2")
    use iso_c_binding

    type(c_ptr), value :: state
    type(c_ptr), value :: global_data
    type(c_ptr), value :: p_diag
    type(c_ptr), value :: p_int
    type(c_ptr), value :: p_metrics
    type(c_ptr), value :: p_patch
    type(c_ptr), value :: p_prog
    type(c_ptr), value :: z_kin_hor_e
    type(c_ptr), value :: z_vt_ie
    type(c_ptr), value :: z_w_concorr_me
    integer(kind=c_int), value :: f2dace_A_z_kin_hor_e_d_0_s_363
    integer(kind=c_int), value :: f2dace_A_z_kin_hor_e_d_1_s_364
    integer(kind=c_int), value :: f2dace_A_z_kin_hor_e_d_2_s_365
    integer(kind=c_int), value :: f2dace_A_z_vt_ie_d_0_s_366
    integer(kind=c_int), value :: f2dace_A_z_vt_ie_d_1_s_367
    integer(kind=c_int), value :: f2dace_A_z_vt_ie_d_2_s_368
    integer(kind=c_int), value :: f2dace_A_z_w_concorr_me_d_0_s_360
    integer(kind=c_int), value :: f2dace_A_z_w_concorr_me_d_1_s_361
    integer(kind=c_int), value :: f2dace_A_z_w_concorr_me_d_2_s_362
    integer(kind=c_int), value :: f2dace_OA_z_kin_hor_e_d_0_s_363
    integer(kind=c_int), value :: f2dace_OA_z_kin_hor_e_d_1_s_364
    integer(kind=c_int), value :: f2dace_OA_z_kin_hor_e_d_2_s_365
    real(kind=c_double), value :: dt_linintp_ubc
    real(kind=c_double), value :: dtime
    integer(kind=c_int), value :: istep
    integer(kind=c_int), value :: ldeepatmo
    integer(kind=c_int), value :: lvn_only
    integer(kind=c_int), value :: ntnd
  end subroutine dace_program_velocity_no_nproma_if_prop_lvn_only_1_istep_2

end interface

interface logical_fix_1d
  module procedure logical_to_int_1d
  module procedure int_to_int_1d
end interface logical_fix_1d

interface logical_fix_2d
  module procedure logical_to_int_2d
  module procedure int_to_int_2d
end interface logical_fix_2d

interface logical_fix_3d
  module procedure logical_to_int_3d
  module procedure int_to_int_3d
end interface logical_fix_3d


  real(8), parameter :: float32_default_rel_threshold = 1.0e-8
  real(8), parameter :: float32_default_abs_threshold = 0.0

  real(8), parameter :: float64_default_rel_threshold = 1.0e-12
  real(8), parameter :: float64_default_abs_threshold = 0.0

  real(8), parameter :: int32_default_rel_threshold = 1.0e-8
  real(8), parameter :: int32_default_abs_threshold = 0.0

  real(8), parameter :: int64_default_rel_threshold = 1.0e-12
  real(8), parameter :: int64_default_abs_threshold = 0.0


contains

  function copy_in_global_data_type(steal_arrays, minimal_structs) result(dace_obj_ptr)
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_global_data_type), pointer :: dace_rich_obj
    type(dace_global_data_type) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%nflatlev = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(nflatlev), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
#if defined(DACE_SUBST_VERIFY)
    if (10 /= size(nflatlev, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 'global_data.nflatlev'"//char(10), &
        "    - actual = (", &
        size(nflatlev, dim=1), &
        "), config propagated = (10)"
    end if
#endif
    dace_rich_obj%i_am_accel_node = i_am_accel_node
    dace_rich_obj%lextra_diffu = lextra_diffu
    dace_rich_obj%nproma = nproma
    dace_rich_obj%timers_level = timers_level
    dace_rich_obj%timer_solve_nh_veltend = timer_solve_nh_veltend
    dace_rich_obj%timer_intp = timer_intp
    dace_rich_obj%nrdmax = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(nrdmax), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
#if defined(DACE_SUBST_VERIFY)
    if (10 /= size(nrdmax, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 'global_data.nrdmax'"//char(10), &
        "    - actual = (", &
        size(nrdmax, dim=1), &
        "), config propagated = (10)"
    end if
#endif

  end function copy_in_global_data_type

  function copy_in_t_nh_diag(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_nh_diag), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_nh_diag), pointer :: dace_rich_obj
    type(dace_t_nh_diag) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_vt_d_0_s_291 = size(fortran_obj%vt, dim=1)
    dace_rich_obj%f2dace_SOA_vt_d_0_s_291 = lbound(fortran_obj%vt, dim=1)
    dace_rich_obj%f2dace_SA_vt_d_1_s_292 = size(fortran_obj%vt, dim=2)
    dace_rich_obj%f2dace_SOA_vt_d_1_s_292 = lbound(fortran_obj%vt, dim=2)
    dace_rich_obj%f2dace_SA_vt_d_2_s_293 = size(fortran_obj%vt, dim=3)
    dace_rich_obj%f2dace_SOA_vt_d_2_s_293 = lbound(fortran_obj%vt, dim=3)
    dace_rich_obj%f2dace_SA_vn_ie_d_0_s_294 = size(fortran_obj%vn_ie, dim=1)
    dace_rich_obj%f2dace_SOA_vn_ie_d_0_s_294 = lbound(fortran_obj%vn_ie, dim=1)
    dace_rich_obj%f2dace_SA_vn_ie_d_1_s_295 = size(fortran_obj%vn_ie, dim=2)
    dace_rich_obj%f2dace_SOA_vn_ie_d_1_s_295 = lbound(fortran_obj%vn_ie, dim=2)
    dace_rich_obj%f2dace_SA_vn_ie_d_2_s_296 = size(fortran_obj%vn_ie, dim=3)
    dace_rich_obj%f2dace_SOA_vn_ie_d_2_s_296 = lbound(fortran_obj%vn_ie, dim=3)
    dace_rich_obj%f2dace_SA_w_concorr_c_d_0_s_297 = size(fortran_obj%w_concorr_c, dim=1)
    dace_rich_obj%f2dace_SOA_w_concorr_c_d_0_s_297 = lbound(fortran_obj%w_concorr_c, dim=1)
    dace_rich_obj%f2dace_SA_w_concorr_c_d_1_s_298 = size(fortran_obj%w_concorr_c, dim=2)
    dace_rich_obj%f2dace_SOA_w_concorr_c_d_1_s_298 = lbound(fortran_obj%w_concorr_c, dim=2)
    dace_rich_obj%f2dace_SA_w_concorr_c_d_2_s_299 = size(fortran_obj%w_concorr_c, dim=3)
    dace_rich_obj%f2dace_SOA_w_concorr_c_d_2_s_299 = lbound(fortran_obj%w_concorr_c, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_apc_pc_d_0_s_300 = size(fortran_obj%ddt_vn_apc_pc, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_apc_pc_d_0_s_300 = lbound(fortran_obj%ddt_vn_apc_pc, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_apc_pc_d_1_s_301 = size(fortran_obj%ddt_vn_apc_pc, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_apc_pc_d_1_s_301 = lbound(fortran_obj%ddt_vn_apc_pc, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_apc_pc_d_2_s_302 = size(fortran_obj%ddt_vn_apc_pc, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_apc_pc_d_2_s_302 = lbound(fortran_obj%ddt_vn_apc_pc, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_apc_pc_d_3_s_303 = size(fortran_obj%ddt_vn_apc_pc, dim=4)
    dace_rich_obj%f2dace_SOA_ddt_vn_apc_pc_d_3_s_303 = lbound(fortran_obj%ddt_vn_apc_pc, dim=4)
    dace_rich_obj%f2dace_SA_ddt_w_adv_pc_d_0_s_304 = size(fortran_obj%ddt_w_adv_pc, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_w_adv_pc_d_0_s_304 = lbound(fortran_obj%ddt_w_adv_pc, dim=1)
    dace_rich_obj%f2dace_SA_ddt_w_adv_pc_d_1_s_305 = size(fortran_obj%ddt_w_adv_pc, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_w_adv_pc_d_1_s_305 = lbound(fortran_obj%ddt_w_adv_pc, dim=2)
    dace_rich_obj%f2dace_SA_ddt_w_adv_pc_d_2_s_306 = size(fortran_obj%ddt_w_adv_pc, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_w_adv_pc_d_2_s_306 = lbound(fortran_obj%ddt_w_adv_pc, dim=3)
    dace_rich_obj%f2dace_SA_ddt_w_adv_pc_d_3_s_307 = size(fortran_obj%ddt_w_adv_pc, dim=4)
    dace_rich_obj%f2dace_SOA_ddt_w_adv_pc_d_3_s_307 = lbound(fortran_obj%ddt_w_adv_pc, dim=4)
    dace_rich_obj%max_vcfl_dyn = fortran_obj%max_vcfl_dyn
#ifndef _OPENACC
    dace_rich_obj%vt = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vt, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vt = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vt, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (90 /= size(fortran_obj%vt, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_diag.vt'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%vt, dim=1), ",", &
        size(fortran_obj%vt, dim=2), ",", &
        size(fortran_obj%vt, dim=3), &
        "), config propagated = (__f2dace_SA_vt_d_0_s_291_p_diag_9, 90, __f2dace_SA_vt_d_2_s_293_p_diag_9)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%vn_ie = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_ie, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vn_ie = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_ie, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (91 /= size(fortran_obj%vn_ie, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_diag.vn_ie'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%vn_ie, dim=1), ",", &
        size(fortran_obj%vn_ie, dim=2), ",", &
        size(fortran_obj%vn_ie, dim=3), &
        "), config propagated = (__f2dace_SA_vn_ie_d_0_s_294_p_diag_9, 91, __f2dace_SA_vn_ie_d_2_s_296_p_diag_9)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%w_concorr_c = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%w_concorr_c, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%w_concorr_c = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%w_concorr_c, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (91 /= size(fortran_obj%w_concorr_c, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_diag.w_concorr_c'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%w_concorr_c, dim=1), ",", &
        size(fortran_obj%w_concorr_c, dim=2), ",", &
        size(fortran_obj%w_concorr_c, dim=3), &
        "), config propagated = (__f2dace_SA_w_concorr_c_d_0_s_297_p_diag_9, 91, __f2dace_SA_w_concorr_c_d_2_s_299_p_diag_9)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_apc_pc = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%ddt_vn_apc_pc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_apc_pc = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%ddt_vn_apc_pc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (90 /= size(fortran_obj%ddt_vn_apc_pc, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_diag.ddt_vn_apc_pc'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%ddt_vn_apc_pc, dim=1), ",", &
        size(fortran_obj%ddt_vn_apc_pc, dim=2), ",", &
        size(fortran_obj%ddt_vn_apc_pc, dim=3), ",", &
        size(fortran_obj%ddt_vn_apc_pc, dim=4), &
        "), config propagated = (__f2dace_SA_ddt_vn_apc_pc_d_0_s_300_p_diag_9, 90, __f2dace_SA_ddt_vn_apc_pc_d_2_s_302_p_diag_9, __f2dace_SA_ddt_vn_apc_pc_d_3_s_303_p_diag_9)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%ddt_w_adv_pc = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%ddt_w_adv_pc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_w_adv_pc = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%ddt_w_adv_pc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (91 /= size(fortran_obj%ddt_w_adv_pc, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_diag.ddt_w_adv_pc'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%ddt_w_adv_pc, dim=1), ",", &
        size(fortran_obj%ddt_w_adv_pc, dim=2), ",", &
        size(fortran_obj%ddt_w_adv_pc, dim=3), ",", &
        size(fortran_obj%ddt_w_adv_pc, dim=4), &
        "), config propagated = (__f2dace_SA_ddt_w_adv_pc_d_0_s_304_p_diag_9, 91, __f2dace_SA_ddt_w_adv_pc_d_2_s_306_p_diag_9, __f2dace_SA_ddt_w_adv_pc_d_3_s_307_p_diag_9)"
    end if
#endif

  end function copy_in_t_nh_diag

  function copy_in_t_int_state(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_int_state), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_int_state), pointer :: dace_rich_obj
    type(dace_t_int_state) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_c_lin_e_d_0_s_25 = size(fortran_obj%c_lin_e, dim=1)
    dace_rich_obj%f2dace_SOA_c_lin_e_d_0_s_25 = lbound(fortran_obj%c_lin_e, dim=1)
    dace_rich_obj%f2dace_SA_c_lin_e_d_1_s_26 = size(fortran_obj%c_lin_e, dim=2)
    dace_rich_obj%f2dace_SOA_c_lin_e_d_1_s_26 = lbound(fortran_obj%c_lin_e, dim=2)
    dace_rich_obj%f2dace_SA_c_lin_e_d_2_s_27 = size(fortran_obj%c_lin_e, dim=3)
    dace_rich_obj%f2dace_SOA_c_lin_e_d_2_s_27 = lbound(fortran_obj%c_lin_e, dim=3)
    dace_rich_obj%f2dace_SA_e_bln_c_s_d_0_s_28 = size(fortran_obj%e_bln_c_s, dim=1)
    dace_rich_obj%f2dace_SOA_e_bln_c_s_d_0_s_28 = lbound(fortran_obj%e_bln_c_s, dim=1)
    dace_rich_obj%f2dace_SA_e_bln_c_s_d_1_s_29 = size(fortran_obj%e_bln_c_s, dim=2)
    dace_rich_obj%f2dace_SOA_e_bln_c_s_d_1_s_29 = lbound(fortran_obj%e_bln_c_s, dim=2)
    dace_rich_obj%f2dace_SA_e_bln_c_s_d_2_s_30 = size(fortran_obj%e_bln_c_s, dim=3)
    dace_rich_obj%f2dace_SOA_e_bln_c_s_d_2_s_30 = lbound(fortran_obj%e_bln_c_s, dim=3)
    dace_rich_obj%f2dace_SA_cells_aw_verts_d_0_s_31 = size(fortran_obj%cells_aw_verts, dim=1)
    dace_rich_obj%f2dace_SOA_cells_aw_verts_d_0_s_31 = lbound(fortran_obj%cells_aw_verts, dim=1)
    dace_rich_obj%f2dace_SA_cells_aw_verts_d_1_s_32 = size(fortran_obj%cells_aw_verts, dim=2)
    dace_rich_obj%f2dace_SOA_cells_aw_verts_d_1_s_32 = lbound(fortran_obj%cells_aw_verts, dim=2)
    dace_rich_obj%f2dace_SA_cells_aw_verts_d_2_s_33 = size(fortran_obj%cells_aw_verts, dim=3)
    dace_rich_obj%f2dace_SOA_cells_aw_verts_d_2_s_33 = lbound(fortran_obj%cells_aw_verts, dim=3)
    dace_rich_obj%f2dace_SA_rbf_vec_coeff_e_d_0_s_34 = size(fortran_obj%rbf_vec_coeff_e, dim=1)
    dace_rich_obj%f2dace_SOA_rbf_vec_coeff_e_d_0_s_34 = lbound(fortran_obj%rbf_vec_coeff_e, dim=1)
    dace_rich_obj%f2dace_SA_rbf_vec_coeff_e_d_1_s_35 = size(fortran_obj%rbf_vec_coeff_e, dim=2)
    dace_rich_obj%f2dace_SOA_rbf_vec_coeff_e_d_1_s_35 = lbound(fortran_obj%rbf_vec_coeff_e, dim=2)
    dace_rich_obj%f2dace_SA_rbf_vec_coeff_e_d_2_s_36 = size(fortran_obj%rbf_vec_coeff_e, dim=3)
    dace_rich_obj%f2dace_SOA_rbf_vec_coeff_e_d_2_s_36 = lbound(fortran_obj%rbf_vec_coeff_e, dim=3)
    dace_rich_obj%f2dace_SA_geofac_grdiv_d_0_s_37 = size(fortran_obj%geofac_grdiv, dim=1)
    dace_rich_obj%f2dace_SOA_geofac_grdiv_d_0_s_37 = lbound(fortran_obj%geofac_grdiv, dim=1)
    dace_rich_obj%f2dace_SA_geofac_grdiv_d_1_s_38 = size(fortran_obj%geofac_grdiv, dim=2)
    dace_rich_obj%f2dace_SOA_geofac_grdiv_d_1_s_38 = lbound(fortran_obj%geofac_grdiv, dim=2)
    dace_rich_obj%f2dace_SA_geofac_grdiv_d_2_s_39 = size(fortran_obj%geofac_grdiv, dim=3)
    dace_rich_obj%f2dace_SOA_geofac_grdiv_d_2_s_39 = lbound(fortran_obj%geofac_grdiv, dim=3)
    dace_rich_obj%f2dace_SA_geofac_rot_d_0_s_40 = size(fortran_obj%geofac_rot, dim=1)
    dace_rich_obj%f2dace_SOA_geofac_rot_d_0_s_40 = lbound(fortran_obj%geofac_rot, dim=1)
    dace_rich_obj%f2dace_SA_geofac_rot_d_1_s_41 = size(fortran_obj%geofac_rot, dim=2)
    dace_rich_obj%f2dace_SOA_geofac_rot_d_1_s_41 = lbound(fortran_obj%geofac_rot, dim=2)
    dace_rich_obj%f2dace_SA_geofac_rot_d_2_s_42 = size(fortran_obj%geofac_rot, dim=3)
    dace_rich_obj%f2dace_SOA_geofac_rot_d_2_s_42 = lbound(fortran_obj%geofac_rot, dim=3)
    dace_rich_obj%f2dace_SA_geofac_n2s_d_0_s_43 = size(fortran_obj%geofac_n2s, dim=1)
    dace_rich_obj%f2dace_SOA_geofac_n2s_d_0_s_43 = lbound(fortran_obj%geofac_n2s, dim=1)
    dace_rich_obj%f2dace_SA_geofac_n2s_d_1_s_44 = size(fortran_obj%geofac_n2s, dim=2)
    dace_rich_obj%f2dace_SOA_geofac_n2s_d_1_s_44 = lbound(fortran_obj%geofac_n2s, dim=2)
    dace_rich_obj%f2dace_SA_geofac_n2s_d_2_s_45 = size(fortran_obj%geofac_n2s, dim=3)
    dace_rich_obj%f2dace_SOA_geofac_n2s_d_2_s_45 = lbound(fortran_obj%geofac_n2s, dim=3)
#ifndef _OPENACC
    dace_rich_obj%c_lin_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%c_lin_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%c_lin_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%c_lin_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%e_bln_c_s = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%e_bln_c_s, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%e_bln_c_s = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%e_bln_c_s, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%cells_aw_verts = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%cells_aw_verts, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%cells_aw_verts = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%cells_aw_verts, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%rbf_vec_coeff_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rbf_vec_coeff_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%rbf_vec_coeff_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rbf_vec_coeff_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%geofac_grdiv = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%geofac_grdiv, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%geofac_grdiv = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%geofac_grdiv, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%geofac_rot = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%geofac_rot, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%geofac_rot = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%geofac_rot, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%geofac_n2s = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%geofac_n2s, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%geofac_n2s = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%geofac_n2s, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif

  end function copy_in_t_int_state

  function copy_in_t_nh_metrics(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_nh_metrics), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_nh_metrics), pointer :: dace_rich_obj
    type(dace_t_nh_metrics) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_ddxn_z_full_d_0_s_308 = size(fortran_obj%ddxn_z_full, dim=1)
    dace_rich_obj%f2dace_SOA_ddxn_z_full_d_0_s_308 = lbound(fortran_obj%ddxn_z_full, dim=1)
    dace_rich_obj%f2dace_SA_ddxn_z_full_d_1_s_309 = size(fortran_obj%ddxn_z_full, dim=2)
    dace_rich_obj%f2dace_SOA_ddxn_z_full_d_1_s_309 = lbound(fortran_obj%ddxn_z_full, dim=2)
    dace_rich_obj%f2dace_SA_ddxn_z_full_d_2_s_310 = size(fortran_obj%ddxn_z_full, dim=3)
    dace_rich_obj%f2dace_SOA_ddxn_z_full_d_2_s_310 = lbound(fortran_obj%ddxn_z_full, dim=3)
    dace_rich_obj%f2dace_SA_ddxt_z_full_d_0_s_311 = size(fortran_obj%ddxt_z_full, dim=1)
    dace_rich_obj%f2dace_SOA_ddxt_z_full_d_0_s_311 = lbound(fortran_obj%ddxt_z_full, dim=1)
    dace_rich_obj%f2dace_SA_ddxt_z_full_d_1_s_312 = size(fortran_obj%ddxt_z_full, dim=2)
    dace_rich_obj%f2dace_SOA_ddxt_z_full_d_1_s_312 = lbound(fortran_obj%ddxt_z_full, dim=2)
    dace_rich_obj%f2dace_SA_ddxt_z_full_d_2_s_313 = size(fortran_obj%ddxt_z_full, dim=3)
    dace_rich_obj%f2dace_SOA_ddxt_z_full_d_2_s_313 = lbound(fortran_obj%ddxt_z_full, dim=3)
    dace_rich_obj%f2dace_SA_ddqz_z_full_e_d_0_s_314 = size(fortran_obj%ddqz_z_full_e, dim=1)
    dace_rich_obj%f2dace_SOA_ddqz_z_full_e_d_0_s_314 = lbound(fortran_obj%ddqz_z_full_e, dim=1)
    dace_rich_obj%f2dace_SA_ddqz_z_full_e_d_1_s_315 = size(fortran_obj%ddqz_z_full_e, dim=2)
    dace_rich_obj%f2dace_SOA_ddqz_z_full_e_d_1_s_315 = lbound(fortran_obj%ddqz_z_full_e, dim=2)
    dace_rich_obj%f2dace_SA_ddqz_z_full_e_d_2_s_316 = size(fortran_obj%ddqz_z_full_e, dim=3)
    dace_rich_obj%f2dace_SOA_ddqz_z_full_e_d_2_s_316 = lbound(fortran_obj%ddqz_z_full_e, dim=3)
    dace_rich_obj%f2dace_SA_ddqz_z_half_d_0_s_317 = size(fortran_obj%ddqz_z_half, dim=1)
    dace_rich_obj%f2dace_SOA_ddqz_z_half_d_0_s_317 = lbound(fortran_obj%ddqz_z_half, dim=1)
    dace_rich_obj%f2dace_SA_ddqz_z_half_d_1_s_318 = size(fortran_obj%ddqz_z_half, dim=2)
    dace_rich_obj%f2dace_SOA_ddqz_z_half_d_1_s_318 = lbound(fortran_obj%ddqz_z_half, dim=2)
    dace_rich_obj%f2dace_SA_ddqz_z_half_d_2_s_319 = size(fortran_obj%ddqz_z_half, dim=3)
    dace_rich_obj%f2dace_SOA_ddqz_z_half_d_2_s_319 = lbound(fortran_obj%ddqz_z_half, dim=3)
    dace_rich_obj%f2dace_SA_wgtfac_c_d_0_s_320 = size(fortran_obj%wgtfac_c, dim=1)
    dace_rich_obj%f2dace_SOA_wgtfac_c_d_0_s_320 = lbound(fortran_obj%wgtfac_c, dim=1)
    dace_rich_obj%f2dace_SA_wgtfac_c_d_1_s_321 = size(fortran_obj%wgtfac_c, dim=2)
    dace_rich_obj%f2dace_SOA_wgtfac_c_d_1_s_321 = lbound(fortran_obj%wgtfac_c, dim=2)
    dace_rich_obj%f2dace_SA_wgtfac_c_d_2_s_322 = size(fortran_obj%wgtfac_c, dim=3)
    dace_rich_obj%f2dace_SOA_wgtfac_c_d_2_s_322 = lbound(fortran_obj%wgtfac_c, dim=3)
    dace_rich_obj%f2dace_SA_wgtfac_e_d_0_s_323 = size(fortran_obj%wgtfac_e, dim=1)
    dace_rich_obj%f2dace_SOA_wgtfac_e_d_0_s_323 = lbound(fortran_obj%wgtfac_e, dim=1)
    dace_rich_obj%f2dace_SA_wgtfac_e_d_1_s_324 = size(fortran_obj%wgtfac_e, dim=2)
    dace_rich_obj%f2dace_SOA_wgtfac_e_d_1_s_324 = lbound(fortran_obj%wgtfac_e, dim=2)
    dace_rich_obj%f2dace_SA_wgtfac_e_d_2_s_325 = size(fortran_obj%wgtfac_e, dim=3)
    dace_rich_obj%f2dace_SOA_wgtfac_e_d_2_s_325 = lbound(fortran_obj%wgtfac_e, dim=3)
    dace_rich_obj%f2dace_SA_wgtfacq_e_d_0_s_326 = size(fortran_obj%wgtfacq_e, dim=1)
    dace_rich_obj%f2dace_SOA_wgtfacq_e_d_0_s_326 = lbound(fortran_obj%wgtfacq_e, dim=1)
    dace_rich_obj%f2dace_SA_wgtfacq_e_d_1_s_327 = size(fortran_obj%wgtfacq_e, dim=2)
    dace_rich_obj%f2dace_SOA_wgtfacq_e_d_1_s_327 = lbound(fortran_obj%wgtfacq_e, dim=2)
    dace_rich_obj%f2dace_SA_wgtfacq_e_d_2_s_328 = size(fortran_obj%wgtfacq_e, dim=3)
    dace_rich_obj%f2dace_SOA_wgtfacq_e_d_2_s_328 = lbound(fortran_obj%wgtfacq_e, dim=3)
    dace_rich_obj%f2dace_SA_coeff_gradekin_d_0_s_329 = size(fortran_obj%coeff_gradekin, dim=1)
    dace_rich_obj%f2dace_SOA_coeff_gradekin_d_0_s_329 = lbound(fortran_obj%coeff_gradekin, dim=1)
    dace_rich_obj%f2dace_SA_coeff_gradekin_d_1_s_330 = size(fortran_obj%coeff_gradekin, dim=2)
    dace_rich_obj%f2dace_SOA_coeff_gradekin_d_1_s_330 = lbound(fortran_obj%coeff_gradekin, dim=2)
    dace_rich_obj%f2dace_SA_coeff_gradekin_d_2_s_331 = size(fortran_obj%coeff_gradekin, dim=3)
    dace_rich_obj%f2dace_SOA_coeff_gradekin_d_2_s_331 = lbound(fortran_obj%coeff_gradekin, dim=3)
    dace_rich_obj%f2dace_SA_coeff1_dwdz_d_0_s_332 = size(fortran_obj%coeff1_dwdz, dim=1)
    dace_rich_obj%f2dace_SOA_coeff1_dwdz_d_0_s_332 = lbound(fortran_obj%coeff1_dwdz, dim=1)
    dace_rich_obj%f2dace_SA_coeff1_dwdz_d_1_s_333 = size(fortran_obj%coeff1_dwdz, dim=2)
    dace_rich_obj%f2dace_SOA_coeff1_dwdz_d_1_s_333 = lbound(fortran_obj%coeff1_dwdz, dim=2)
    dace_rich_obj%f2dace_SA_coeff1_dwdz_d_2_s_334 = size(fortran_obj%coeff1_dwdz, dim=3)
    dace_rich_obj%f2dace_SOA_coeff1_dwdz_d_2_s_334 = lbound(fortran_obj%coeff1_dwdz, dim=3)
    dace_rich_obj%f2dace_SA_coeff2_dwdz_d_0_s_335 = size(fortran_obj%coeff2_dwdz, dim=1)
    dace_rich_obj%f2dace_SOA_coeff2_dwdz_d_0_s_335 = lbound(fortran_obj%coeff2_dwdz, dim=1)
    dace_rich_obj%f2dace_SA_coeff2_dwdz_d_1_s_336 = size(fortran_obj%coeff2_dwdz, dim=2)
    dace_rich_obj%f2dace_SOA_coeff2_dwdz_d_1_s_336 = lbound(fortran_obj%coeff2_dwdz, dim=2)
    dace_rich_obj%f2dace_SA_coeff2_dwdz_d_2_s_337 = size(fortran_obj%coeff2_dwdz, dim=3)
    dace_rich_obj%f2dace_SOA_coeff2_dwdz_d_2_s_337 = lbound(fortran_obj%coeff2_dwdz, dim=3)
    dace_rich_obj%f2dace_SA_deepatmo_gradh_mc_d_0_s_338 = size(fortran_obj%deepatmo_gradh_mc, dim=1)
    dace_rich_obj%f2dace_SOA_deepatmo_gradh_mc_d_0_s_338 = lbound(fortran_obj%deepatmo_gradh_mc, dim=1)
    dace_rich_obj%deepatmo_gradh_mc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_gradh_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
#if defined(DACE_SUBST_VERIFY)
    if (90 /= size(fortran_obj%deepatmo_gradh_mc, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.deepatmo_gradh_mc'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%deepatmo_gradh_mc, dim=1), &
        "), config propagated = (90)"
    end if
#endif
    dace_rich_obj%f2dace_SA_deepatmo_invr_mc_d_0_s_339 = size(fortran_obj%deepatmo_invr_mc, dim=1)
    dace_rich_obj%f2dace_SOA_deepatmo_invr_mc_d_0_s_339 = lbound(fortran_obj%deepatmo_invr_mc, dim=1)
    dace_rich_obj%deepatmo_invr_mc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_invr_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
#if defined(DACE_SUBST_VERIFY)
    if (90 /= size(fortran_obj%deepatmo_invr_mc, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.deepatmo_invr_mc'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%deepatmo_invr_mc, dim=1), &
        "), config propagated = (90)"
    end if
#endif
    dace_rich_obj%f2dace_SA_deepatmo_gradh_ifc_d_0_s_340 = size(fortran_obj%deepatmo_gradh_ifc, dim=1)
    dace_rich_obj%f2dace_SOA_deepatmo_gradh_ifc_d_0_s_340 = lbound(fortran_obj%deepatmo_gradh_ifc, dim=1)
    dace_rich_obj%deepatmo_gradh_ifc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_gradh_ifc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
#if defined(DACE_SUBST_VERIFY)
    if (91 /= size(fortran_obj%deepatmo_gradh_ifc, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.deepatmo_gradh_ifc'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%deepatmo_gradh_ifc, dim=1), &
        "), config propagated = (91)"
    end if
#endif
    dace_rich_obj%f2dace_SA_deepatmo_invr_ifc_d_0_s_341 = size(fortran_obj%deepatmo_invr_ifc, dim=1)
    dace_rich_obj%f2dace_SOA_deepatmo_invr_ifc_d_0_s_341 = lbound(fortran_obj%deepatmo_invr_ifc, dim=1)
    dace_rich_obj%deepatmo_invr_ifc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_invr_ifc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
#if defined(DACE_SUBST_VERIFY)
    if (91 /= size(fortran_obj%deepatmo_invr_ifc, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.deepatmo_invr_ifc'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%deepatmo_invr_ifc, dim=1), &
        "), config propagated = (91)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%ddxn_z_full = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddxn_z_full, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddxn_z_full = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddxn_z_full, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (90 /= size(fortran_obj%ddxn_z_full, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.ddxn_z_full'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%ddxn_z_full, dim=1), ",", &
        size(fortran_obj%ddxn_z_full, dim=2), ",", &
        size(fortran_obj%ddxn_z_full, dim=3), &
        "), config propagated = (__f2dace_SA_ddxn_z_full_d_0_s_308_p_metrics_8, 90, __f2dace_SA_ddxn_z_full_d_2_s_310_p_metrics_8)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%ddxt_z_full = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddxt_z_full, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddxt_z_full = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddxt_z_full, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (90 /= size(fortran_obj%ddxt_z_full, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.ddxt_z_full'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%ddxt_z_full, dim=1), ",", &
        size(fortran_obj%ddxt_z_full, dim=2), ",", &
        size(fortran_obj%ddxt_z_full, dim=3), &
        "), config propagated = (__f2dace_SA_ddxt_z_full_d_0_s_311_p_metrics_8, 90, __f2dace_SA_ddxt_z_full_d_2_s_313_p_metrics_8)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%ddqz_z_full_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddqz_z_full_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddqz_z_full_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddqz_z_full_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (90 /= size(fortran_obj%ddqz_z_full_e, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.ddqz_z_full_e'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%ddqz_z_full_e, dim=1), ",", &
        size(fortran_obj%ddqz_z_full_e, dim=2), ",", &
        size(fortran_obj%ddqz_z_full_e, dim=3), &
        "), config propagated = (__f2dace_SA_ddqz_z_full_e_d_0_s_314_p_metrics_8, 90, __f2dace_SA_ddqz_z_full_e_d_2_s_316_p_metrics_8)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%ddqz_z_half = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddqz_z_half, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddqz_z_half = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddqz_z_half, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (91 /= size(fortran_obj%ddqz_z_half, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.ddqz_z_half'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%ddqz_z_half, dim=1), ",", &
        size(fortran_obj%ddqz_z_half, dim=2), ",", &
        size(fortran_obj%ddqz_z_half, dim=3), &
        "), config propagated = (__f2dace_SA_ddqz_z_half_d_0_s_317_p_metrics_8, 91, __f2dace_SA_ddqz_z_half_d_2_s_319_p_metrics_8)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%wgtfac_c = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%wgtfac_c, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%wgtfac_c = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%wgtfac_c, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (91 /= size(fortran_obj%wgtfac_c, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.wgtfac_c'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%wgtfac_c, dim=1), ",", &
        size(fortran_obj%wgtfac_c, dim=2), ",", &
        size(fortran_obj%wgtfac_c, dim=3), &
        "), config propagated = (__f2dace_SA_wgtfac_c_d_0_s_320_p_metrics_8, 91, __f2dace_SA_wgtfac_c_d_2_s_322_p_metrics_8)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%wgtfac_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%wgtfac_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%wgtfac_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%wgtfac_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (91 /= size(fortran_obj%wgtfac_e, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.wgtfac_e'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%wgtfac_e, dim=1), ",", &
        size(fortran_obj%wgtfac_e, dim=2), ",", &
        size(fortran_obj%wgtfac_e, dim=3), &
        "), config propagated = (__f2dace_SA_wgtfac_e_d_0_s_323_p_metrics_8, 91, __f2dace_SA_wgtfac_e_d_2_s_325_p_metrics_8)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%wgtfacq_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%wgtfacq_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%wgtfacq_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%wgtfacq_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%coeff_gradekin = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%coeff_gradekin, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%coeff_gradekin = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%coeff_gradekin, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%coeff1_dwdz = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%coeff1_dwdz, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%coeff1_dwdz = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%coeff1_dwdz, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (90 /= size(fortran_obj%coeff1_dwdz, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.coeff1_dwdz'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%coeff1_dwdz, dim=1), ",", &
        size(fortran_obj%coeff1_dwdz, dim=2), ",", &
        size(fortran_obj%coeff1_dwdz, dim=3), &
        "), config propagated = (__f2dace_SA_coeff1_dwdz_d_0_s_332_p_metrics_8, 90, __f2dace_SA_coeff1_dwdz_d_2_s_334_p_metrics_8)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%coeff2_dwdz = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%coeff2_dwdz, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%coeff2_dwdz = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%coeff2_dwdz, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (90 /= size(fortran_obj%coeff2_dwdz, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_metrics.coeff2_dwdz'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%coeff2_dwdz, dim=1), ",", &
        size(fortran_obj%coeff2_dwdz, dim=2), ",", &
        size(fortran_obj%coeff2_dwdz, dim=3), &
        "), config propagated = (__f2dace_SA_coeff2_dwdz_d_0_s_335_p_metrics_8, 90, __f2dace_SA_coeff2_dwdz_d_2_s_337_p_metrics_8)"
    end if
#endif

  end function copy_in_t_nh_metrics

  function copy_in_t_patch(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_patch), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_patch), pointer :: dace_rich_obj
    type(dace_t_patch) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%nblks_c = fortran_obj%nblks_c
    dace_rich_obj%nblks_e = fortran_obj%nblks_e
    dace_rich_obj%nblks_v = fortran_obj%nblks_v
    dace_rich_obj%cells = copy_in_t_grid_cells( &
    fortran_obj=fortran_obj%cells, &
    steal_arrays=steal_arrays, &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%edges = copy_in_t_grid_edges( &
    fortran_obj=fortran_obj%edges, &
    steal_arrays=steal_arrays, &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%verts = copy_in_t_grid_vertices( &
    fortran_obj=fortran_obj%verts, &
    steal_arrays=steal_arrays, &
    minimal_structs=minimal_structs &
  )

  end function copy_in_t_patch

  function copy_in_t_grid_cells(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_grid_cells), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_grid_cells), pointer :: dace_rich_obj
    type(dace_t_grid_cells) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_neighbor_idx_d_0_s_146 = size(fortran_obj%neighbor_idx, dim=1)
    dace_rich_obj%f2dace_SOA_neighbor_idx_d_0_s_146 = lbound(fortran_obj%neighbor_idx, dim=1)
    dace_rich_obj%f2dace_SA_neighbor_idx_d_1_s_147 = size(fortran_obj%neighbor_idx, dim=2)
    dace_rich_obj%f2dace_SOA_neighbor_idx_d_1_s_147 = lbound(fortran_obj%neighbor_idx, dim=2)
    dace_rich_obj%f2dace_SA_neighbor_idx_d_2_s_148 = size(fortran_obj%neighbor_idx, dim=3)
    dace_rich_obj%f2dace_SOA_neighbor_idx_d_2_s_148 = lbound(fortran_obj%neighbor_idx, dim=3)
    dace_rich_obj%f2dace_SA_neighbor_blk_d_0_s_149 = size(fortran_obj%neighbor_blk, dim=1)
    dace_rich_obj%f2dace_SOA_neighbor_blk_d_0_s_149 = lbound(fortran_obj%neighbor_blk, dim=1)
    dace_rich_obj%f2dace_SA_neighbor_blk_d_1_s_150 = size(fortran_obj%neighbor_blk, dim=2)
    dace_rich_obj%f2dace_SOA_neighbor_blk_d_1_s_150 = lbound(fortran_obj%neighbor_blk, dim=2)
    dace_rich_obj%f2dace_SA_neighbor_blk_d_2_s_151 = size(fortran_obj%neighbor_blk, dim=3)
    dace_rich_obj%f2dace_SOA_neighbor_blk_d_2_s_151 = lbound(fortran_obj%neighbor_blk, dim=3)
    dace_rich_obj%f2dace_SA_edge_idx_d_0_s_152 = size(fortran_obj%edge_idx, dim=1)
    dace_rich_obj%f2dace_SOA_edge_idx_d_0_s_152 = lbound(fortran_obj%edge_idx, dim=1)
    dace_rich_obj%f2dace_SA_edge_idx_d_1_s_153 = size(fortran_obj%edge_idx, dim=2)
    dace_rich_obj%f2dace_SOA_edge_idx_d_1_s_153 = lbound(fortran_obj%edge_idx, dim=2)
    dace_rich_obj%f2dace_SA_edge_idx_d_2_s_154 = size(fortran_obj%edge_idx, dim=3)
    dace_rich_obj%f2dace_SOA_edge_idx_d_2_s_154 = lbound(fortran_obj%edge_idx, dim=3)
    dace_rich_obj%f2dace_SA_edge_blk_d_0_s_155 = size(fortran_obj%edge_blk, dim=1)
    dace_rich_obj%f2dace_SOA_edge_blk_d_0_s_155 = lbound(fortran_obj%edge_blk, dim=1)
    dace_rich_obj%f2dace_SA_edge_blk_d_1_s_156 = size(fortran_obj%edge_blk, dim=2)
    dace_rich_obj%f2dace_SOA_edge_blk_d_1_s_156 = lbound(fortran_obj%edge_blk, dim=2)
    dace_rich_obj%f2dace_SA_edge_blk_d_2_s_157 = size(fortran_obj%edge_blk, dim=3)
    dace_rich_obj%f2dace_SOA_edge_blk_d_2_s_157 = lbound(fortran_obj%edge_blk, dim=3)
    dace_rich_obj%f2dace_SA_area_d_0_s_158 = size(fortran_obj%area, dim=1)
    dace_rich_obj%f2dace_SOA_area_d_0_s_158 = lbound(fortran_obj%area, dim=1)
    dace_rich_obj%f2dace_SA_area_d_1_s_159 = size(fortran_obj%area, dim=2)
    dace_rich_obj%f2dace_SOA_area_d_1_s_159 = lbound(fortran_obj%area, dim=2)
    dace_rich_obj%f2dace_SA_start_index_d_0_s_160 = size(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SOA_start_index_d_0_s_160 = lbound(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SA_end_index_d_0_s_161 = size(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SOA_end_index_d_0_s_161 = lbound(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SA_start_block_d_0_s_162 = size(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SOA_start_block_d_0_s_162 = lbound(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SA_end_block_d_0_s_163 = size(fortran_obj%end_block, dim=1)
    dace_rich_obj%f2dace_SOA_end_block_d_0_s_163 = lbound(fortran_obj%end_block, dim=1)
    dace_rich_obj%decomp_info = copy_in_t_grid_domain_decomp_info( &
    fortran_obj=fortran_obj%decomp_info, &
    steal_arrays=steal_arrays, &
    minimal_structs=minimal_structs &
  )
#ifndef _OPENACC
    dace_rich_obj%neighbor_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%neighbor_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%neighbor_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%neighbor_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%neighbor_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%neighbor_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%neighbor_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%neighbor_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%edge_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%edge_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%edge_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%edge_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%edge_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%edge_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%edge_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%edge_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%area = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%area, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%area = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%area, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
    dace_rich_obj%start_index = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%start_index), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%end_index = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%end_index), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%start_block = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%start_block), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%end_block = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%end_block), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

  end function copy_in_t_grid_cells

  function copy_in_t_grid_domain_decomp_info(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_grid_domain_decomp_info), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_grid_domain_decomp_info), pointer :: dace_rich_obj
    type(dace_t_grid_domain_decomp_info) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_owner_mask_d_0_s_2 = size(fortran_obj%owner_mask, dim=1)
    dace_rich_obj%f2dace_SOA_owner_mask_d_0_s_2 = lbound(fortran_obj%owner_mask, dim=1)
    dace_rich_obj%f2dace_SA_owner_mask_d_1_s_3 = size(fortran_obj%owner_mask, dim=2)
    dace_rich_obj%f2dace_SOA_owner_mask_d_1_s_3 = lbound(fortran_obj%owner_mask, dim=2)
#ifndef _OPENACC
    dace_rich_obj%owner_mask = copy_in_int32_2d_array( &
    fortran_array=logical_fix_2d(fortran_obj%owner_mask), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%owner_mask = copy_in_int32_2d_array( &
    fortran_array=logical_fix_2d(fortran_obj%owner_mask), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif

  end function copy_in_t_grid_domain_decomp_info

  function copy_in_t_grid_edges(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_grid_edges), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_grid_edges), pointer :: dace_rich_obj
    type(dace_t_grid_edges) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_cell_idx_d_0_s_164 = size(fortran_obj%cell_idx, dim=1)
    dace_rich_obj%f2dace_SOA_cell_idx_d_0_s_164 = lbound(fortran_obj%cell_idx, dim=1)
    dace_rich_obj%f2dace_SA_cell_idx_d_1_s_165 = size(fortran_obj%cell_idx, dim=2)
    dace_rich_obj%f2dace_SOA_cell_idx_d_1_s_165 = lbound(fortran_obj%cell_idx, dim=2)
    dace_rich_obj%f2dace_SA_cell_idx_d_2_s_166 = size(fortran_obj%cell_idx, dim=3)
    dace_rich_obj%f2dace_SOA_cell_idx_d_2_s_166 = lbound(fortran_obj%cell_idx, dim=3)
    dace_rich_obj%f2dace_SA_cell_blk_d_0_s_167 = size(fortran_obj%cell_blk, dim=1)
    dace_rich_obj%f2dace_SOA_cell_blk_d_0_s_167 = lbound(fortran_obj%cell_blk, dim=1)
    dace_rich_obj%f2dace_SA_cell_blk_d_1_s_168 = size(fortran_obj%cell_blk, dim=2)
    dace_rich_obj%f2dace_SOA_cell_blk_d_1_s_168 = lbound(fortran_obj%cell_blk, dim=2)
    dace_rich_obj%f2dace_SA_cell_blk_d_2_s_169 = size(fortran_obj%cell_blk, dim=3)
    dace_rich_obj%f2dace_SOA_cell_blk_d_2_s_169 = lbound(fortran_obj%cell_blk, dim=3)
    dace_rich_obj%f2dace_SA_vertex_idx_d_0_s_170 = size(fortran_obj%vertex_idx, dim=1)
    dace_rich_obj%f2dace_SOA_vertex_idx_d_0_s_170 = lbound(fortran_obj%vertex_idx, dim=1)
    dace_rich_obj%f2dace_SA_vertex_idx_d_1_s_171 = size(fortran_obj%vertex_idx, dim=2)
    dace_rich_obj%f2dace_SOA_vertex_idx_d_1_s_171 = lbound(fortran_obj%vertex_idx, dim=2)
    dace_rich_obj%f2dace_SA_vertex_idx_d_2_s_172 = size(fortran_obj%vertex_idx, dim=3)
    dace_rich_obj%f2dace_SOA_vertex_idx_d_2_s_172 = lbound(fortran_obj%vertex_idx, dim=3)
    dace_rich_obj%f2dace_SA_vertex_blk_d_0_s_173 = size(fortran_obj%vertex_blk, dim=1)
    dace_rich_obj%f2dace_SOA_vertex_blk_d_0_s_173 = lbound(fortran_obj%vertex_blk, dim=1)
    dace_rich_obj%f2dace_SA_vertex_blk_d_1_s_174 = size(fortran_obj%vertex_blk, dim=2)
    dace_rich_obj%f2dace_SOA_vertex_blk_d_1_s_174 = lbound(fortran_obj%vertex_blk, dim=2)
    dace_rich_obj%f2dace_SA_vertex_blk_d_2_s_175 = size(fortran_obj%vertex_blk, dim=3)
    dace_rich_obj%f2dace_SOA_vertex_blk_d_2_s_175 = lbound(fortran_obj%vertex_blk, dim=3)
    dace_rich_obj%f2dace_SA_tangent_orientation_d_0_s_176 = size(fortran_obj%tangent_orientation, dim=1)
    dace_rich_obj%f2dace_SOA_tangent_orientation_d_0_s_176 = lbound(fortran_obj%tangent_orientation, dim=1)
    dace_rich_obj%f2dace_SA_tangent_orientation_d_1_s_177 = size(fortran_obj%tangent_orientation, dim=2)
    dace_rich_obj%f2dace_SOA_tangent_orientation_d_1_s_177 = lbound(fortran_obj%tangent_orientation, dim=2)
    dace_rich_obj%f2dace_SA_quad_idx_d_0_s_178 = size(fortran_obj%quad_idx, dim=1)
    dace_rich_obj%f2dace_SOA_quad_idx_d_0_s_178 = lbound(fortran_obj%quad_idx, dim=1)
    dace_rich_obj%f2dace_SA_quad_idx_d_1_s_179 = size(fortran_obj%quad_idx, dim=2)
    dace_rich_obj%f2dace_SOA_quad_idx_d_1_s_179 = lbound(fortran_obj%quad_idx, dim=2)
    dace_rich_obj%f2dace_SA_quad_idx_d_2_s_180 = size(fortran_obj%quad_idx, dim=3)
    dace_rich_obj%f2dace_SOA_quad_idx_d_2_s_180 = lbound(fortran_obj%quad_idx, dim=3)
    dace_rich_obj%f2dace_SA_quad_blk_d_0_s_181 = size(fortran_obj%quad_blk, dim=1)
    dace_rich_obj%f2dace_SOA_quad_blk_d_0_s_181 = lbound(fortran_obj%quad_blk, dim=1)
    dace_rich_obj%f2dace_SA_quad_blk_d_1_s_182 = size(fortran_obj%quad_blk, dim=2)
    dace_rich_obj%f2dace_SOA_quad_blk_d_1_s_182 = lbound(fortran_obj%quad_blk, dim=2)
    dace_rich_obj%f2dace_SA_quad_blk_d_2_s_183 = size(fortran_obj%quad_blk, dim=3)
    dace_rich_obj%f2dace_SOA_quad_blk_d_2_s_183 = lbound(fortran_obj%quad_blk, dim=3)
    dace_rich_obj%f2dace_SA_inv_primal_edge_length_d_0_s_184 = size(fortran_obj%inv_primal_edge_length, dim=1)
    dace_rich_obj%f2dace_SOA_inv_primal_edge_length_d_0_s_184 = lbound(fortran_obj%inv_primal_edge_length, dim=1)
    dace_rich_obj%f2dace_SA_inv_primal_edge_length_d_1_s_185 = size(fortran_obj%inv_primal_edge_length, dim=2)
    dace_rich_obj%f2dace_SOA_inv_primal_edge_length_d_1_s_185 = lbound(fortran_obj%inv_primal_edge_length, dim=2)
    dace_rich_obj%f2dace_SA_inv_dual_edge_length_d_0_s_186 = size(fortran_obj%inv_dual_edge_length, dim=1)
    dace_rich_obj%f2dace_SOA_inv_dual_edge_length_d_0_s_186 = lbound(fortran_obj%inv_dual_edge_length, dim=1)
    dace_rich_obj%f2dace_SA_inv_dual_edge_length_d_1_s_187 = size(fortran_obj%inv_dual_edge_length, dim=2)
    dace_rich_obj%f2dace_SOA_inv_dual_edge_length_d_1_s_187 = lbound(fortran_obj%inv_dual_edge_length, dim=2)
    dace_rich_obj%f2dace_SA_area_edge_d_0_s_188 = size(fortran_obj%area_edge, dim=1)
    dace_rich_obj%f2dace_SOA_area_edge_d_0_s_188 = lbound(fortran_obj%area_edge, dim=1)
    dace_rich_obj%f2dace_SA_area_edge_d_1_s_189 = size(fortran_obj%area_edge, dim=2)
    dace_rich_obj%f2dace_SOA_area_edge_d_1_s_189 = lbound(fortran_obj%area_edge, dim=2)
    dace_rich_obj%f2dace_SA_f_e_d_0_s_190 = size(fortran_obj%f_e, dim=1)
    dace_rich_obj%f2dace_SOA_f_e_d_0_s_190 = lbound(fortran_obj%f_e, dim=1)
    dace_rich_obj%f2dace_SA_f_e_d_1_s_191 = size(fortran_obj%f_e, dim=2)
    dace_rich_obj%f2dace_SOA_f_e_d_1_s_191 = lbound(fortran_obj%f_e, dim=2)
    dace_rich_obj%f2dace_SA_fn_e_d_0_s_192 = size(fortran_obj%fn_e, dim=1)
    dace_rich_obj%f2dace_SOA_fn_e_d_0_s_192 = lbound(fortran_obj%fn_e, dim=1)
    dace_rich_obj%f2dace_SA_fn_e_d_1_s_193 = size(fortran_obj%fn_e, dim=2)
    dace_rich_obj%f2dace_SOA_fn_e_d_1_s_193 = lbound(fortran_obj%fn_e, dim=2)
    dace_rich_obj%f2dace_SA_ft_e_d_0_s_194 = size(fortran_obj%ft_e, dim=1)
    dace_rich_obj%f2dace_SOA_ft_e_d_0_s_194 = lbound(fortran_obj%ft_e, dim=1)
    dace_rich_obj%f2dace_SA_ft_e_d_1_s_195 = size(fortran_obj%ft_e, dim=2)
    dace_rich_obj%f2dace_SOA_ft_e_d_1_s_195 = lbound(fortran_obj%ft_e, dim=2)
    dace_rich_obj%f2dace_SA_start_index_d_0_s_196 = size(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SOA_start_index_d_0_s_196 = lbound(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SA_end_index_d_0_s_197 = size(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SOA_end_index_d_0_s_197 = lbound(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SA_start_block_d_0_s_198 = size(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SOA_start_block_d_0_s_198 = lbound(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SA_end_block_d_0_s_199 = size(fortran_obj%end_block, dim=1)
    dace_rich_obj%f2dace_SOA_end_block_d_0_s_199 = lbound(fortran_obj%end_block, dim=1)
#ifndef _OPENACC
    dace_rich_obj%cell_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%cell_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%cell_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%cell_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%cell_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%cell_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%cell_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%cell_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%vertex_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%vertex_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vertex_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%vertex_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%vertex_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%vertex_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vertex_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%vertex_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%tangent_orientation = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%tangent_orientation, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%tangent_orientation = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%tangent_orientation, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%quad_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%quad_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%quad_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%quad_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%quad_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%quad_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%quad_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%quad_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%inv_primal_edge_length = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%inv_primal_edge_length, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%inv_primal_edge_length = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%inv_primal_edge_length, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%inv_dual_edge_length = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%inv_dual_edge_length, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%inv_dual_edge_length = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%inv_dual_edge_length, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%area_edge = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%area_edge, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%area_edge = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%area_edge, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%f_e = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%f_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%f_e = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%f_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
    dace_rich_obj%fn_e = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%fn_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%ft_e = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%ft_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%start_index = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%start_index), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%end_index = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%end_index), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%start_block = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%start_block), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%end_block = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%end_block), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

  end function copy_in_t_grid_edges

  function copy_in_t_grid_vertices(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_grid_vertices), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_grid_vertices), pointer :: dace_rich_obj
    type(dace_t_grid_vertices) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_cell_idx_d_0_s_200 = size(fortran_obj%cell_idx, dim=1)
    dace_rich_obj%f2dace_SOA_cell_idx_d_0_s_200 = lbound(fortran_obj%cell_idx, dim=1)
    dace_rich_obj%f2dace_SA_cell_idx_d_1_s_201 = size(fortran_obj%cell_idx, dim=2)
    dace_rich_obj%f2dace_SOA_cell_idx_d_1_s_201 = lbound(fortran_obj%cell_idx, dim=2)
    dace_rich_obj%f2dace_SA_cell_idx_d_2_s_202 = size(fortran_obj%cell_idx, dim=3)
    dace_rich_obj%f2dace_SOA_cell_idx_d_2_s_202 = lbound(fortran_obj%cell_idx, dim=3)
    dace_rich_obj%f2dace_SA_cell_blk_d_0_s_203 = size(fortran_obj%cell_blk, dim=1)
    dace_rich_obj%f2dace_SOA_cell_blk_d_0_s_203 = lbound(fortran_obj%cell_blk, dim=1)
    dace_rich_obj%f2dace_SA_cell_blk_d_1_s_204 = size(fortran_obj%cell_blk, dim=2)
    dace_rich_obj%f2dace_SOA_cell_blk_d_1_s_204 = lbound(fortran_obj%cell_blk, dim=2)
    dace_rich_obj%f2dace_SA_cell_blk_d_2_s_205 = size(fortran_obj%cell_blk, dim=3)
    dace_rich_obj%f2dace_SOA_cell_blk_d_2_s_205 = lbound(fortran_obj%cell_blk, dim=3)
    dace_rich_obj%f2dace_SA_edge_idx_d_0_s_206 = size(fortran_obj%edge_idx, dim=1)
    dace_rich_obj%f2dace_SOA_edge_idx_d_0_s_206 = lbound(fortran_obj%edge_idx, dim=1)
    dace_rich_obj%f2dace_SA_edge_idx_d_1_s_207 = size(fortran_obj%edge_idx, dim=2)
    dace_rich_obj%f2dace_SOA_edge_idx_d_1_s_207 = lbound(fortran_obj%edge_idx, dim=2)
    dace_rich_obj%f2dace_SA_edge_idx_d_2_s_208 = size(fortran_obj%edge_idx, dim=3)
    dace_rich_obj%f2dace_SOA_edge_idx_d_2_s_208 = lbound(fortran_obj%edge_idx, dim=3)
    dace_rich_obj%f2dace_SA_edge_blk_d_0_s_209 = size(fortran_obj%edge_blk, dim=1)
    dace_rich_obj%f2dace_SOA_edge_blk_d_0_s_209 = lbound(fortran_obj%edge_blk, dim=1)
    dace_rich_obj%f2dace_SA_edge_blk_d_1_s_210 = size(fortran_obj%edge_blk, dim=2)
    dace_rich_obj%f2dace_SOA_edge_blk_d_1_s_210 = lbound(fortran_obj%edge_blk, dim=2)
    dace_rich_obj%f2dace_SA_edge_blk_d_2_s_211 = size(fortran_obj%edge_blk, dim=3)
    dace_rich_obj%f2dace_SOA_edge_blk_d_2_s_211 = lbound(fortran_obj%edge_blk, dim=3)
    dace_rich_obj%f2dace_SA_start_index_d_0_s_212 = size(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SOA_start_index_d_0_s_212 = lbound(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SA_end_index_d_0_s_213 = size(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SOA_end_index_d_0_s_213 = lbound(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SA_start_block_d_0_s_214 = size(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SOA_start_block_d_0_s_214 = lbound(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SA_end_block_d_0_s_215 = size(fortran_obj%end_block, dim=1)
    dace_rich_obj%f2dace_SOA_end_block_d_0_s_215 = lbound(fortran_obj%end_block, dim=1)
#ifndef _OPENACC
    dace_rich_obj%cell_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%cell_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%cell_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%cell_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%cell_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%cell_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%cell_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%cell_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%edge_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%edge_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%edge_idx = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%edge_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%edge_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%edge_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%edge_blk = copy_in_int32_3d_array( &
    fortran_array=logical_fix_3d(fortran_obj%edge_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
    dace_rich_obj%start_index = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%start_index), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%end_index = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%end_index), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%start_block = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%start_block), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%end_block = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%end_block), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

  end function copy_in_t_grid_vertices

  function copy_in_t_nh_prog(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_nh_prog), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_nh_prog), pointer :: dace_rich_obj
    type(dace_t_nh_prog) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_w_d_0_s_285 = size(fortran_obj%w, dim=1)
    dace_rich_obj%f2dace_SOA_w_d_0_s_285 = lbound(fortran_obj%w, dim=1)
    dace_rich_obj%f2dace_SA_w_d_1_s_286 = size(fortran_obj%w, dim=2)
    dace_rich_obj%f2dace_SOA_w_d_1_s_286 = lbound(fortran_obj%w, dim=2)
    dace_rich_obj%f2dace_SA_w_d_2_s_287 = size(fortran_obj%w, dim=3)
    dace_rich_obj%f2dace_SOA_w_d_2_s_287 = lbound(fortran_obj%w, dim=3)
    dace_rich_obj%f2dace_SA_vn_d_0_s_288 = size(fortran_obj%vn, dim=1)
    dace_rich_obj%f2dace_SOA_vn_d_0_s_288 = lbound(fortran_obj%vn, dim=1)
    dace_rich_obj%f2dace_SA_vn_d_1_s_289 = size(fortran_obj%vn, dim=2)
    dace_rich_obj%f2dace_SOA_vn_d_1_s_289 = lbound(fortran_obj%vn, dim=2)
    dace_rich_obj%f2dace_SA_vn_d_2_s_290 = size(fortran_obj%vn, dim=3)
    dace_rich_obj%f2dace_SOA_vn_d_2_s_290 = lbound(fortran_obj%vn, dim=3)
#ifndef _OPENACC
    dace_rich_obj%w = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%w, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%w = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%w, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (91 /= size(fortran_obj%w, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_prog.w'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%w, dim=1), ",", &
        size(fortran_obj%w, dim=2), ",", &
        size(fortran_obj%w, dim=3), &
        "), config propagated = (__f2dace_SA_w_d_0_s_285_p_prog_7, 91, __f2dace_SA_w_d_2_s_287_p_prog_7)"
    end if
#endif
#ifndef _OPENACC
    dace_rich_obj%vn = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vn = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (90 /= size(fortran_obj%vn, dim=2)) then
      print *, &
        "Array size conflicts with config propagation for array 't_nh_prog.vn'"//char(10), &
        "    - actual = (", &
        size(fortran_obj%vn, dim=1), ",", &
        size(fortran_obj%vn, dim=2), ",", &
        size(fortran_obj%vn, dim=3), &
        "), config propagated = (__f2dace_SA_vn_d_0_s_288_p_prog_7, 90, __f2dace_SA_vn_d_2_s_290_p_prog_7)"
    end if
#endif

  end function copy_in_t_nh_prog

  function logical_to_int_1d(inp) result(out)
    logical(4), dimension(:), target :: inp
    integer(kind=c_int), dimension(:), pointer :: out

    call c_f_pointer(c_loc(inp), out, shape=shape(inp))
  end function logical_to_int_1d

  function int_to_int_1d(inp) result(out)
    integer(kind=c_int), dimension(:), target :: inp
    integer(kind=c_int), dimension(:), pointer :: out

    out => inp
  end function int_to_int_1d

  function copy_in_int32_1d_array( &
    fortran_array, &
    steal_arrays, &
    use_openacc, &
    minimal_structs &
   ) &
   result(dace_array_ptr)
    integer(kind=c_int), dimension(:), target :: fortran_array
    logical :: steal_arrays, use_openacc, minimal_structs
    type(c_ptr) :: dace_array_ptr
    integer(kind=c_int), dimension(:), pointer :: dace_rich_array

    integer :: i0
#ifdef _OPENACC
    integer(kind=c_size_t) :: size_bytes
#endif

    if (use_openacc) then
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#endif
    end if

    if (.not. c_associated(c_loc(fortran_array))) then
      dace_array_ptr = c_null_ptr
      return
    end if

    if (steal_arrays .eqv. .true.) then
      if (use_openacc .eqv. .false.) then
          dace_array_ptr = c_loc(fortran_array)
      else
#ifdef _OPENACC
          dace_array_ptr = c_acc_deviceptr(c_loc(fortran_array))
#endif
      end if
      return
    end if

    if (use_openacc .eqv. .false.) then
        dace_array_ptr = malloc(c_sizeof(dace_array_ptr) * size(fortran_array))
        call c_f_pointer(dace_array_ptr, dace_rich_array, shape=shape(fortran_array))

    do i0 = 1, size(fortran_array, dim=1)
         dace_rich_array(i0) = fortran_array(i0)
    end do

    else
#ifdef _OPENACC
    size_bytes = size(fortran_array) * c_sizeof(fortran_array(1))
    dace_array_ptr = c_acc_malloc(size_bytes)
    call c_acc_memcpy_device(dace_array_ptr, c_acc_deviceptr(c_loc(fortran_array)), size_bytes)
#endif
    end if
  end function copy_in_int32_1d_array

  function copy_in_float64_3d_array( &
    fortran_array, &
    steal_arrays, &
    use_openacc, &
    minimal_structs &
   ) &
   result(dace_array_ptr)
    real(kind=c_double), dimension(:,:,:), target :: fortran_array
    logical :: steal_arrays, use_openacc, minimal_structs
    type(c_ptr) :: dace_array_ptr
    real(kind=c_double), dimension(:,:,:), pointer :: dace_rich_array

    integer :: i0, i1, i2
#ifdef _OPENACC
    integer(kind=c_size_t) :: size_bytes
#endif

    if (use_openacc) then
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#endif
    end if

    if (.not. c_associated(c_loc(fortran_array))) then
      dace_array_ptr = c_null_ptr
      return
    end if

    if (steal_arrays .eqv. .true.) then
      if (use_openacc .eqv. .false.) then
          dace_array_ptr = c_loc(fortran_array)
      else
#ifdef _OPENACC
          dace_array_ptr = c_acc_deviceptr(c_loc(fortran_array))
#endif
      end if
      return
    end if

    if (use_openacc .eqv. .false.) then
        dace_array_ptr = malloc(c_sizeof(dace_array_ptr) * size(fortran_array))
        call c_f_pointer(dace_array_ptr, dace_rich_array, shape=shape(fortran_array))

    do i0 = 1, size(fortran_array, dim=1)
      do i1 = 1, size(fortran_array, dim=2)
        do i2 = 1, size(fortran_array, dim=3)
                   dace_rich_array(i0, i1, i2) = fortran_array(i0, i1, i2)
        end do
      end do
    end do

    else
#ifdef _OPENACC
    size_bytes = size(fortran_array) * c_sizeof(fortran_array(1, 1, 1))
    dace_array_ptr = c_acc_malloc(size_bytes)
    call c_acc_memcpy_device(dace_array_ptr, c_acc_deviceptr(c_loc(fortran_array)), size_bytes)
#endif
    end if
  end function copy_in_float64_3d_array

  function copy_in_float64_4d_array( &
    fortran_array, &
    steal_arrays, &
    use_openacc, &
    minimal_structs &
   ) &
   result(dace_array_ptr)
    real(kind=c_double), dimension(:,:,:,:), target :: fortran_array
    logical :: steal_arrays, use_openacc, minimal_structs
    type(c_ptr) :: dace_array_ptr
    real(kind=c_double), dimension(:,:,:,:), pointer :: dace_rich_array

    integer :: i0, i1, i2, i3
#ifdef _OPENACC
    integer(kind=c_size_t) :: size_bytes
#endif

    if (use_openacc) then
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#endif
    end if

    if (.not. c_associated(c_loc(fortran_array))) then
      dace_array_ptr = c_null_ptr
      return
    end if

    if (steal_arrays .eqv. .true.) then
      if (use_openacc .eqv. .false.) then
          dace_array_ptr = c_loc(fortran_array)
      else
#ifdef _OPENACC
          dace_array_ptr = c_acc_deviceptr(c_loc(fortran_array))
#endif
      end if
      return
    end if

    if (use_openacc .eqv. .false.) then
        dace_array_ptr = malloc(c_sizeof(dace_array_ptr) * size(fortran_array))
        call c_f_pointer(dace_array_ptr, dace_rich_array, shape=shape(fortran_array))

    do i0 = 1, size(fortran_array, dim=1)
      do i1 = 1, size(fortran_array, dim=2)
        do i2 = 1, size(fortran_array, dim=3)
          do i3 = 1, size(fortran_array, dim=4)
                        dace_rich_array(i0, i1, i2, i3) = fortran_array(i0, i1, i2, i3)
          end do
        end do
      end do
    end do

    else
#ifdef _OPENACC
    size_bytes = size(fortran_array) * c_sizeof(fortran_array(1, 1, 1, 1))
    dace_array_ptr = c_acc_malloc(size_bytes)
    call c_acc_memcpy_device(dace_array_ptr, c_acc_deviceptr(c_loc(fortran_array)), size_bytes)
#endif
    end if
  end function copy_in_float64_4d_array

  function copy_in_float64_1d_array( &
    fortran_array, &
    steal_arrays, &
    use_openacc, &
    minimal_structs &
   ) &
   result(dace_array_ptr)
    real(kind=c_double), dimension(:), target :: fortran_array
    logical :: steal_arrays, use_openacc, minimal_structs
    type(c_ptr) :: dace_array_ptr
    real(kind=c_double), dimension(:), pointer :: dace_rich_array

    integer :: i0
#ifdef _OPENACC
    integer(kind=c_size_t) :: size_bytes
#endif

    if (use_openacc) then
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#endif
    end if

    if (.not. c_associated(c_loc(fortran_array))) then
      dace_array_ptr = c_null_ptr
      return
    end if

    if (steal_arrays .eqv. .true.) then
      if (use_openacc .eqv. .false.) then
          dace_array_ptr = c_loc(fortran_array)
      else
#ifdef _OPENACC
          dace_array_ptr = c_acc_deviceptr(c_loc(fortran_array))
#endif
      end if
      return
    end if

    if (use_openacc .eqv. .false.) then
        dace_array_ptr = malloc(c_sizeof(dace_array_ptr) * size(fortran_array))
        call c_f_pointer(dace_array_ptr, dace_rich_array, shape=shape(fortran_array))

    do i0 = 1, size(fortran_array, dim=1)
         dace_rich_array(i0) = fortran_array(i0)
    end do

    else
#ifdef _OPENACC
    size_bytes = size(fortran_array) * c_sizeof(fortran_array(1))
    dace_array_ptr = c_acc_malloc(size_bytes)
    call c_acc_memcpy_device(dace_array_ptr, c_acc_deviceptr(c_loc(fortran_array)), size_bytes)
#endif
    end if
  end function copy_in_float64_1d_array

  function logical_to_int_2d(inp) result(out)
    logical(4), dimension(:,:), target :: inp
    integer(kind=c_int), dimension(:,:), pointer :: out

    call c_f_pointer(c_loc(inp), out, shape=shape(inp))
  end function logical_to_int_2d

  function int_to_int_2d(inp) result(out)
    integer(kind=c_int), dimension(:,:), target :: inp
    integer(kind=c_int), dimension(:,:), pointer :: out

    out => inp
  end function int_to_int_2d

  function copy_in_int32_2d_array( &
    fortran_array, &
    steal_arrays, &
    use_openacc, &
    minimal_structs &
   ) &
   result(dace_array_ptr)
    integer(kind=c_int), dimension(:,:), target :: fortran_array
    logical :: steal_arrays, use_openacc, minimal_structs
    type(c_ptr) :: dace_array_ptr
    integer(kind=c_int), dimension(:,:), pointer :: dace_rich_array

    integer :: i0, i1
#ifdef _OPENACC
    integer(kind=c_size_t) :: size_bytes
#endif

    if (use_openacc) then
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#endif
    end if

    if (.not. c_associated(c_loc(fortran_array))) then
      dace_array_ptr = c_null_ptr
      return
    end if

    if (steal_arrays .eqv. .true.) then
      if (use_openacc .eqv. .false.) then
          dace_array_ptr = c_loc(fortran_array)
      else
#ifdef _OPENACC
          dace_array_ptr = c_acc_deviceptr(c_loc(fortran_array))
#endif
      end if
      return
    end if

    if (use_openacc .eqv. .false.) then
        dace_array_ptr = malloc(c_sizeof(dace_array_ptr) * size(fortran_array))
        call c_f_pointer(dace_array_ptr, dace_rich_array, shape=shape(fortran_array))

    do i0 = 1, size(fortran_array, dim=1)
      do i1 = 1, size(fortran_array, dim=2)
              dace_rich_array(i0, i1) = fortran_array(i0, i1)
      end do
    end do

    else
#ifdef _OPENACC
    size_bytes = size(fortran_array) * c_sizeof(fortran_array(1, 1))
    dace_array_ptr = c_acc_malloc(size_bytes)
    call c_acc_memcpy_device(dace_array_ptr, c_acc_deviceptr(c_loc(fortran_array)), size_bytes)
#endif
    end if
  end function copy_in_int32_2d_array

  function logical_to_int_3d(inp) result(out)
    logical(4), dimension(:,:,:), target :: inp
    integer(kind=c_int), dimension(:,:,:), pointer :: out

    call c_f_pointer(c_loc(inp), out, shape=shape(inp))
  end function logical_to_int_3d

  function int_to_int_3d(inp) result(out)
    integer(kind=c_int), dimension(:,:,:), target :: inp
    integer(kind=c_int), dimension(:,:,:), pointer :: out

    out => inp
  end function int_to_int_3d

  function copy_in_int32_3d_array( &
    fortran_array, &
    steal_arrays, &
    use_openacc, &
    minimal_structs &
   ) &
   result(dace_array_ptr)
    integer(kind=c_int), dimension(:,:,:), target :: fortran_array
    logical :: steal_arrays, use_openacc, minimal_structs
    type(c_ptr) :: dace_array_ptr
    integer(kind=c_int), dimension(:,:,:), pointer :: dace_rich_array

    integer :: i0, i1, i2
#ifdef _OPENACC
    integer(kind=c_size_t) :: size_bytes
#endif

    if (use_openacc) then
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#endif
    end if

    if (.not. c_associated(c_loc(fortran_array))) then
      dace_array_ptr = c_null_ptr
      return
    end if

    if (steal_arrays .eqv. .true.) then
      if (use_openacc .eqv. .false.) then
          dace_array_ptr = c_loc(fortran_array)
      else
#ifdef _OPENACC
          dace_array_ptr = c_acc_deviceptr(c_loc(fortran_array))
#endif
      end if
      return
    end if

    if (use_openacc .eqv. .false.) then
        dace_array_ptr = malloc(c_sizeof(dace_array_ptr) * size(fortran_array))
        call c_f_pointer(dace_array_ptr, dace_rich_array, shape=shape(fortran_array))

    do i0 = 1, size(fortran_array, dim=1)
      do i1 = 1, size(fortran_array, dim=2)
        do i2 = 1, size(fortran_array, dim=3)
                   dace_rich_array(i0, i1, i2) = fortran_array(i0, i1, i2)
        end do
      end do
    end do

    else
#ifdef _OPENACC
    size_bytes = size(fortran_array) * c_sizeof(fortran_array(1, 1, 1))
    dace_array_ptr = c_acc_malloc(size_bytes)
    call c_acc_memcpy_device(dace_array_ptr, c_acc_deviceptr(c_loc(fortran_array)), size_bytes)
#endif
    end if
  end function copy_in_int32_3d_array

  function copy_in_float64_2d_array( &
    fortran_array, &
    steal_arrays, &
    use_openacc, &
    minimal_structs &
   ) &
   result(dace_array_ptr)
    real(kind=c_double), dimension(:,:), target :: fortran_array
    logical :: steal_arrays, use_openacc, minimal_structs
    type(c_ptr) :: dace_array_ptr
    real(kind=c_double), dimension(:,:), pointer :: dace_rich_array

    integer :: i0, i1
#ifdef _OPENACC
    integer(kind=c_size_t) :: size_bytes
#endif

    if (use_openacc) then
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#endif
    end if

    if (.not. c_associated(c_loc(fortran_array))) then
      dace_array_ptr = c_null_ptr
      return
    end if

    if (steal_arrays .eqv. .true.) then
      if (use_openacc .eqv. .false.) then
          dace_array_ptr = c_loc(fortran_array)
      else
#ifdef _OPENACC
          dace_array_ptr = c_acc_deviceptr(c_loc(fortran_array))
#endif
      end if
      return
    end if

    if (use_openacc .eqv. .false.) then
        dace_array_ptr = malloc(c_sizeof(dace_array_ptr) * size(fortran_array))
        call c_f_pointer(dace_array_ptr, dace_rich_array, shape=shape(fortran_array))

    do i0 = 1, size(fortran_array, dim=1)
      do i1 = 1, size(fortran_array, dim=2)
              dace_rich_array(i0, i1) = fortran_array(i0, i1)
      end do
    end do

    else
#ifdef _OPENACC
    size_bytes = size(fortran_array) * c_sizeof(fortran_array(1, 1))
    dace_array_ptr = c_acc_malloc(size_bytes)
    call c_acc_memcpy_device(dace_array_ptr, c_acc_deviceptr(c_loc(fortran_array)), size_bytes)
#endif
    end if
  end function copy_in_float64_2d_array


  subroutine copy_back_global_data_type(dace_obj_ptr)
    type(c_ptr) :: dace_obj_ptr

    type(dace_global_data_type), pointer :: dace_rich_obj

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    i_am_accel_node = transfer(dace_rich_obj%i_am_accel_node, mold=i_am_accel_node)
    lextra_diffu = transfer(dace_rich_obj%lextra_diffu, mold=lextra_diffu)
    nproma = transfer(dace_rich_obj%nproma, mold=nproma)
    timers_level = transfer(dace_rich_obj%timers_level, mold=timers_level)
    timer_solve_nh_veltend = transfer(dace_rich_obj%timer_solve_nh_veltend, mold=timer_solve_nh_veltend)
    timer_intp = transfer(dace_rich_obj%timer_intp, mold=timer_intp)


    call free(dace_obj_ptr)

  end subroutine copy_back_global_data_type
  subroutine copy_back_t_nh_diag(fortran_obj, dace_obj_ptr)
    type(t_nh_diag), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_nh_diag), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_nh_diag: Invalid allocation of t_nh_diag by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    fortran_obj%max_vcfl_dyn = dace_rich_obj%max_vcfl_dyn


    call free(dace_obj_ptr)

  end subroutine copy_back_t_nh_diag
  subroutine copy_back_t_int_state(fortran_obj, dace_obj_ptr)
    type(t_int_state), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_int_state), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_int_state: Invalid allocation of t_int_state by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



    call free(dace_obj_ptr)

  end subroutine copy_back_t_int_state
  subroutine copy_back_t_nh_metrics(fortran_obj, dace_obj_ptr)
    type(t_nh_metrics), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_nh_metrics), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_nh_metrics: Invalid allocation of t_nh_metrics by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



    call free(dace_obj_ptr)

  end subroutine copy_back_t_nh_metrics
  subroutine copy_back_t_patch(fortran_obj, dace_obj_ptr)
    type(t_patch), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_patch), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_patch: Invalid allocation of t_patch by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    fortran_obj%nblks_c = transfer(dace_rich_obj%nblks_c, mold=fortran_obj%nblks_c)
    fortran_obj%nblks_e = transfer(dace_rich_obj%nblks_e, mold=fortran_obj%nblks_e)
    fortran_obj%nblks_v = transfer(dace_rich_obj%nblks_v, mold=fortran_obj%nblks_v)
    call copy_back_t_grid_cells(fortran_obj%cells, dace_rich_obj%cells)
    call copy_back_t_grid_edges(fortran_obj%edges, dace_rich_obj%edges)
    call copy_back_t_grid_vertices(fortran_obj%verts, dace_rich_obj%verts)


    call free(dace_obj_ptr)

  end subroutine copy_back_t_patch
  subroutine copy_back_t_grid_cells(fortran_obj, dace_obj_ptr)
    type(t_grid_cells), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_grid_cells), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_grid_cells: Invalid allocation of t_grid_cells by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    call copy_back_t_grid_domain_decomp_info(fortran_obj%decomp_info, dace_rich_obj%decomp_info)


    call free(dace_obj_ptr)

  end subroutine copy_back_t_grid_cells
  subroutine copy_back_t_grid_domain_decomp_info(fortran_obj, dace_obj_ptr)
    type(t_grid_domain_decomp_info), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_grid_domain_decomp_info), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_grid_domain_decomp_info: Invalid allocation of t_grid_domain_decomp_info by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



    call free(dace_obj_ptr)

  end subroutine copy_back_t_grid_domain_decomp_info
  subroutine copy_back_t_grid_edges(fortran_obj, dace_obj_ptr)
    type(t_grid_edges), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_grid_edges), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_grid_edges: Invalid allocation of t_grid_edges by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



    call free(dace_obj_ptr)

  end subroutine copy_back_t_grid_edges
  subroutine copy_back_t_grid_vertices(fortran_obj, dace_obj_ptr)
    type(t_grid_vertices), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_grid_vertices), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_grid_vertices: Invalid allocation of t_grid_vertices by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



    call free(dace_obj_ptr)

  end subroutine copy_back_t_grid_vertices
  subroutine copy_back_t_nh_prog(fortran_obj, dace_obj_ptr)
    type(t_nh_prog), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_nh_prog), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_nh_prog: Invalid allocation of t_nh_prog by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



    call free(dace_obj_ptr)

  end subroutine copy_back_t_nh_prog

  subroutine check_initializations()
    CHARACTER(len=5000) :: message_text = ''

  end subroutine check_initializations



  subroutine compare_float32_scalar( &
    actual, &
    ref, &
    result, &
    rel_threshold, &
    abs_threshold, &
    scalar_expr &
  )
    real(kind=c_float), intent(in) :: actual, ref
    logical, intent(out) :: result
    real(8), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in), optional :: scalar_expr

    real(8) :: rel_error, abs_error, threshold_ratio
    real(8) :: actual_rel_threshold, actual_abs_threshold
    CHARACTER(len=5000) :: message_text = ''

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = float32_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = float32_default_abs_threshold
    end if

    result = abs(ref - actual) <= max(actual_rel_threshold * abs(ref), actual_abs_threshold)

    if (present(scalar_expr) .and. .not. result) then

      threshold_ratio = real(abs(ref - actual), kind=8) / max(actual_rel_threshold * abs(ref), actual_abs_threshold)
      rel_error = abs(real(ref - actual, kind=8)/ref)
      abs_error = abs(ref - actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,e28.20,a,e28.20)') &
        "Verification failed for scalar '", &
          trim(scalar_expr), &
        "'"//char(10)//"    - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
          ", threshold_ratio = ", &
          threshold_ratio, &
        char(10)//"    - ref = ", &
          ref, &
          ", actual = ", &
          actual
        print *, "compare_float32_scalar"
        print *, trim(message_text)

    end if

  end subroutine compare_float32_scalar


  subroutine compare_float64_scalar( &
    actual, &
    ref, &
    result, &
    rel_threshold, &
    abs_threshold, &
    scalar_expr &
  )
    real(kind=c_double), intent(in) :: actual, ref
    logical, intent(out) :: result
    real(8), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in), optional :: scalar_expr

    real(8) :: rel_error, abs_error, threshold_ratio
    real(8) :: actual_rel_threshold, actual_abs_threshold
    CHARACTER(len=5000) :: message_text = ''

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = float64_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = float64_default_abs_threshold
    end if

    result = abs(ref - actual) <= max(actual_rel_threshold * abs(ref), actual_abs_threshold)

    if (present(scalar_expr) .and. .not. result) then

      threshold_ratio = real(abs(ref - actual), kind=8) / max(actual_rel_threshold * abs(ref), actual_abs_threshold)
      rel_error = abs(real(ref - actual, kind=8)/ref)
      abs_error = abs(ref - actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,e28.20,a,e28.20)') &
        "Verification failed for scalar '", &
          trim(scalar_expr), &
        "'"//char(10)//"    - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
          ", threshold_ratio = ", &
          threshold_ratio, &
        char(10)//"    - ref = ", &
          ref, &
          ", actual = ", &
          actual
        print *, "compare_float64_scalar"
        print *, trim(message_text)

    end if

  end subroutine compare_float64_scalar


  subroutine compare_int32_scalar( &
    actual, &
    ref, &
    result, &
    rel_threshold, &
    abs_threshold, &
    scalar_expr &
  )
    integer(kind=c_int), intent(in) :: actual, ref
    logical, intent(out) :: result
    real(8), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in), optional :: scalar_expr

    real(8) :: rel_error, abs_error, threshold_ratio
    real(8) :: actual_rel_threshold, actual_abs_threshold
    CHARACTER(len=5000) :: message_text = ''

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = int32_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = int32_default_abs_threshold
    end if

    result = abs(ref - actual) <= max(actual_rel_threshold * abs(ref), actual_abs_threshold)

    if (present(scalar_expr) .and. .not. result) then

      threshold_ratio = real(abs(ref - actual), kind=8) / max(actual_rel_threshold * abs(ref), actual_abs_threshold)
      rel_error = abs(real(ref - actual, kind=8)/ref)
      abs_error = abs(ref - actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,e28.20,a,e28.20)') &
        "Verification failed for scalar '", &
          trim(scalar_expr), &
        "'"//char(10)//"    - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
          ", threshold_ratio = ", &
          threshold_ratio, &
        char(10)//"    - ref = ", &
          ref, &
          ", actual = ", &
          actual
        print *, "compare_int32_scalar"
        print *, trim(message_text)

    end if

  end subroutine compare_int32_scalar


  subroutine compare_int64_scalar( &
    actual, &
    ref, &
    result, &
    rel_threshold, &
    abs_threshold, &
    scalar_expr &
  )
    integer(kind=c_long), intent(in) :: actual, ref
    logical, intent(out) :: result
    real(8), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in), optional :: scalar_expr

    real(8) :: rel_error, abs_error, threshold_ratio
    real(8) :: actual_rel_threshold, actual_abs_threshold
    CHARACTER(len=5000) :: message_text = ''

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = int64_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = int64_default_abs_threshold
    end if

    result = abs(ref - actual) <= max(actual_rel_threshold * abs(ref), actual_abs_threshold)

    if (present(scalar_expr) .and. .not. result) then

      threshold_ratio = real(abs(ref - actual), kind=8) / max(actual_rel_threshold * abs(ref), actual_abs_threshold)
      rel_error = abs(real(ref - actual, kind=8)/ref)
      abs_error = abs(ref - actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,e28.20,a,e28.20)') &
        "Verification failed for scalar '", &
          trim(scalar_expr), &
        "'"//char(10)//"    - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
          ", threshold_ratio = ", &
          threshold_ratio, &
        char(10)//"    - ref = ", &
          ref, &
          ", actual = ", &
          actual
        print *, "compare_int64_scalar"
        print *, trim(message_text)

    end if

  end subroutine compare_int64_scalar

  subroutine compare_global_data_type_struct( &
    actual, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_global_data_type), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%nflatlev"

    call compare_int32_1d_array( &
        actual=actual_rich%nflatlev, &
        ref=logical_fix_1d(nflatlev), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%i_am_accel_node"

    call compare_int32_scalar( &
      actual=actual_rich%i_am_accel_node, &
      ref=transfer(i_am_accel_node, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%lextra_diffu"

    call compare_int32_scalar( &
      actual=actual_rich%lextra_diffu, &
      ref=transfer(lextra_diffu, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%nproma"

    call compare_int32_scalar( &
      actual=actual_rich%nproma, &
      ref=transfer(nproma, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%timers_level"

    call compare_int32_scalar( &
      actual=actual_rich%timers_level, &
      ref=transfer(timers_level, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%timer_solve_nh_veltend"

    call compare_int32_scalar( &
      actual=actual_rich%timer_solve_nh_veltend, &
      ref=transfer(timer_solve_nh_veltend, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%timer_intp"

    call compare_int32_scalar( &
      actual=actual_rich%timer_intp, &
      ref=transfer(timer_intp, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%nrdmax"

    call compare_int32_1d_array( &
        actual=actual_rich%nrdmax, &
        ref=logical_fix_1d(nrdmax), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result


    call free(actual)

  end subroutine compare_global_data_type_struct

  subroutine compare_t_nh_diag_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_nh_diag), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_nh_diag), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%max_vcfl_dyn"

    call compare_float64_scalar( &
      actual=actual_rich%max_vcfl_dyn, &
      ref=ref%max_vcfl_dyn, &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vt"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%vt, &
        ref=ref%vt, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%vt, &
        ref=ref%vt, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vn_ie"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%vn_ie, &
        ref=ref%vn_ie, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%vn_ie, &
        ref=ref%vn_ie, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%w_concorr_c"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%w_concorr_c, &
        ref=ref%w_concorr_c, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%w_concorr_c, &
        ref=ref%w_concorr_c, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_apc_pc"
#ifndef _OPENACC

    call compare_float64_4d_array( &
        actual=actual_rich%ddt_vn_apc_pc, &
        ref=ref%ddt_vn_apc_pc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_4d_array( &
        actual=actual_rich%ddt_vn_apc_pc, &
        ref=ref%ddt_vn_apc_pc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_w_adv_pc"
#ifndef _OPENACC

    call compare_float64_4d_array( &
        actual=actual_rich%ddt_w_adv_pc, &
        ref=ref%ddt_w_adv_pc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_4d_array( &
        actual=actual_rich%ddt_w_adv_pc, &
        ref=ref%ddt_w_adv_pc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_nh_diag_struct

  subroutine compare_t_int_state_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_int_state), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_int_state), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%c_lin_e"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%c_lin_e, &
        ref=ref%c_lin_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%c_lin_e, &
        ref=ref%c_lin_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%e_bln_c_s"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%e_bln_c_s, &
        ref=ref%e_bln_c_s, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%e_bln_c_s, &
        ref=ref%e_bln_c_s, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%cells_aw_verts"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%cells_aw_verts, &
        ref=ref%cells_aw_verts, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%cells_aw_verts, &
        ref=ref%cells_aw_verts, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%rbf_vec_coeff_e"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%rbf_vec_coeff_e, &
        ref=ref%rbf_vec_coeff_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%rbf_vec_coeff_e, &
        ref=ref%rbf_vec_coeff_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%geofac_grdiv"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%geofac_grdiv, &
        ref=ref%geofac_grdiv, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%geofac_grdiv, &
        ref=ref%geofac_grdiv, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%geofac_rot"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%geofac_rot, &
        ref=ref%geofac_rot, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%geofac_rot, &
        ref=ref%geofac_rot, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%geofac_n2s"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%geofac_n2s, &
        ref=ref%geofac_n2s, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%geofac_n2s, &
        ref=ref%geofac_n2s, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_int_state_struct

  subroutine compare_t_nh_metrics_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_nh_metrics), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_nh_metrics), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%deepatmo_gradh_mc"

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_gradh_mc, &
        ref=ref%deepatmo_gradh_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%deepatmo_invr_mc"

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_invr_mc, &
        ref=ref%deepatmo_invr_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%deepatmo_gradh_ifc"

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_gradh_ifc, &
        ref=ref%deepatmo_gradh_ifc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%deepatmo_invr_ifc"

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_invr_ifc, &
        ref=ref%deepatmo_invr_ifc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddxn_z_full"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddxn_z_full, &
        ref=ref%ddxn_z_full, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddxn_z_full, &
        ref=ref%ddxn_z_full, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddxt_z_full"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddxt_z_full, &
        ref=ref%ddxt_z_full, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddxt_z_full, &
        ref=ref%ddxt_z_full, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddqz_z_full_e"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddqz_z_full_e, &
        ref=ref%ddqz_z_full_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddqz_z_full_e, &
        ref=ref%ddqz_z_full_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddqz_z_half"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddqz_z_half, &
        ref=ref%ddqz_z_half, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddqz_z_half, &
        ref=ref%ddqz_z_half, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%wgtfac_c"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%wgtfac_c, &
        ref=ref%wgtfac_c, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%wgtfac_c, &
        ref=ref%wgtfac_c, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%wgtfac_e"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%wgtfac_e, &
        ref=ref%wgtfac_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%wgtfac_e, &
        ref=ref%wgtfac_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%wgtfacq_e"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%wgtfacq_e, &
        ref=ref%wgtfacq_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%wgtfacq_e, &
        ref=ref%wgtfacq_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%coeff_gradekin"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%coeff_gradekin, &
        ref=ref%coeff_gradekin, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%coeff_gradekin, &
        ref=ref%coeff_gradekin, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%coeff1_dwdz"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%coeff1_dwdz, &
        ref=ref%coeff1_dwdz, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%coeff1_dwdz, &
        ref=ref%coeff1_dwdz, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%coeff2_dwdz"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%coeff2_dwdz, &
        ref=ref%coeff2_dwdz, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%coeff2_dwdz, &
        ref=ref%coeff2_dwdz, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_nh_metrics_struct

  subroutine compare_t_patch_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_patch), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_patch), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%nblks_c"

    call compare_int32_scalar( &
      actual=actual_rich%nblks_c, &
      ref=transfer(ref%nblks_c, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%nblks_e"

    call compare_int32_scalar( &
      actual=actual_rich%nblks_e, &
      ref=transfer(ref%nblks_e, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%nblks_v"

    call compare_int32_scalar( &
      actual=actual_rich%nblks_v, &
      ref=transfer(ref%nblks_v, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%cells"

    call compare_t_grid_cells_struct( &
        actual=actual_rich%cells, &
        ref=ref%cells, &
        result=local_result, &
        struct_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%edges"

    call compare_t_grid_edges_struct( &
        actual=actual_rich%edges, &
        ref=ref%edges, &
        result=local_result, &
        struct_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%verts"

    call compare_t_grid_vertices_struct( &
        actual=actual_rich%verts, &
        ref=ref%verts, &
        result=local_result, &
        struct_expr=member_expr &
    )

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_patch_struct

  subroutine compare_t_grid_cells_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_grid_cells), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_grid_cells), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%decomp_info"

    call compare_t_grid_domain_decomp_info_struct( &
        actual=actual_rich%decomp_info, &
        ref=ref%decomp_info, &
        result=local_result, &
        struct_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%neighbor_idx"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%neighbor_idx, &
        ref=logical_fix_3d(ref%neighbor_idx), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%neighbor_idx, &
        ref=logical_fix_3d(ref%neighbor_idx), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%neighbor_blk"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%neighbor_blk, &
        ref=logical_fix_3d(ref%neighbor_blk), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%neighbor_blk, &
        ref=logical_fix_3d(ref%neighbor_blk), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%edge_idx"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%edge_idx, &
        ref=logical_fix_3d(ref%edge_idx), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%edge_idx, &
        ref=logical_fix_3d(ref%edge_idx), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%edge_blk"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%edge_blk, &
        ref=logical_fix_3d(ref%edge_blk), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%edge_blk, &
        ref=logical_fix_3d(ref%edge_blk), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%area"
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=actual_rich%area, &
        ref=ref%area, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_2d_array( &
        actual=actual_rich%area, &
        ref=ref%area, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%start_index"

    call compare_int32_1d_array( &
        actual=actual_rich%start_index, &
        ref=logical_fix_1d(ref%start_index), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%end_index"

    call compare_int32_1d_array( &
        actual=actual_rich%end_index, &
        ref=logical_fix_1d(ref%end_index), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%start_block"

    call compare_int32_1d_array( &
        actual=actual_rich%start_block, &
        ref=logical_fix_1d(ref%start_block), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%end_block"

    call compare_int32_1d_array( &
        actual=actual_rich%end_block, &
        ref=logical_fix_1d(ref%end_block), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_grid_cells_struct

  subroutine compare_t_grid_domain_decomp_info_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_grid_domain_decomp_info), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_grid_domain_decomp_info), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%owner_mask"
#ifndef _OPENACC

    call compare_int32_2d_array( &
        actual=actual_rich%owner_mask, &
        ref=logical_fix_2d(ref%owner_mask), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_2d_array( &
        actual=actual_rich%owner_mask, &
        ref=logical_fix_2d(ref%owner_mask), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_grid_domain_decomp_info_struct

  subroutine compare_t_grid_edges_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_grid_edges), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_grid_edges), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%cell_idx"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%cell_idx, &
        ref=logical_fix_3d(ref%cell_idx), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%cell_idx, &
        ref=logical_fix_3d(ref%cell_idx), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%cell_blk"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%cell_blk, &
        ref=logical_fix_3d(ref%cell_blk), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%cell_blk, &
        ref=logical_fix_3d(ref%cell_blk), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vertex_idx"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%vertex_idx, &
        ref=logical_fix_3d(ref%vertex_idx), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%vertex_idx, &
        ref=logical_fix_3d(ref%vertex_idx), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vertex_blk"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%vertex_blk, &
        ref=logical_fix_3d(ref%vertex_blk), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%vertex_blk, &
        ref=logical_fix_3d(ref%vertex_blk), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%tangent_orientation"
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=actual_rich%tangent_orientation, &
        ref=ref%tangent_orientation, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_2d_array( &
        actual=actual_rich%tangent_orientation, &
        ref=ref%tangent_orientation, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%quad_idx"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%quad_idx, &
        ref=logical_fix_3d(ref%quad_idx), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%quad_idx, &
        ref=logical_fix_3d(ref%quad_idx), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%quad_blk"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%quad_blk, &
        ref=logical_fix_3d(ref%quad_blk), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%quad_blk, &
        ref=logical_fix_3d(ref%quad_blk), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%inv_primal_edge_length"
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=actual_rich%inv_primal_edge_length, &
        ref=ref%inv_primal_edge_length, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_2d_array( &
        actual=actual_rich%inv_primal_edge_length, &
        ref=ref%inv_primal_edge_length, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%inv_dual_edge_length"
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=actual_rich%inv_dual_edge_length, &
        ref=ref%inv_dual_edge_length, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_2d_array( &
        actual=actual_rich%inv_dual_edge_length, &
        ref=ref%inv_dual_edge_length, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%area_edge"
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=actual_rich%area_edge, &
        ref=ref%area_edge, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_2d_array( &
        actual=actual_rich%area_edge, &
        ref=ref%area_edge, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%f_e"
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=actual_rich%f_e, &
        ref=ref%f_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_2d_array( &
        actual=actual_rich%f_e, &
        ref=ref%f_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%fn_e"

    call compare_float64_2d_array( &
        actual=actual_rich%fn_e, &
        ref=ref%fn_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ft_e"

    call compare_float64_2d_array( &
        actual=actual_rich%ft_e, &
        ref=ref%ft_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%start_index"

    call compare_int32_1d_array( &
        actual=actual_rich%start_index, &
        ref=logical_fix_1d(ref%start_index), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%end_index"

    call compare_int32_1d_array( &
        actual=actual_rich%end_index, &
        ref=logical_fix_1d(ref%end_index), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%start_block"

    call compare_int32_1d_array( &
        actual=actual_rich%start_block, &
        ref=logical_fix_1d(ref%start_block), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%end_block"

    call compare_int32_1d_array( &
        actual=actual_rich%end_block, &
        ref=logical_fix_1d(ref%end_block), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_grid_edges_struct

  subroutine compare_t_grid_vertices_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_grid_vertices), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_grid_vertices), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%cell_idx"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%cell_idx, &
        ref=logical_fix_3d(ref%cell_idx), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%cell_idx, &
        ref=logical_fix_3d(ref%cell_idx), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%cell_blk"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%cell_blk, &
        ref=logical_fix_3d(ref%cell_blk), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%cell_blk, &
        ref=logical_fix_3d(ref%cell_blk), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%edge_idx"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%edge_idx, &
        ref=logical_fix_3d(ref%edge_idx), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%edge_idx, &
        ref=logical_fix_3d(ref%edge_idx), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%edge_blk"
#ifndef _OPENACC

    call compare_int32_3d_array( &
        actual=actual_rich%edge_blk, &
        ref=logical_fix_3d(ref%edge_blk), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_3d_array( &
        actual=actual_rich%edge_blk, &
        ref=logical_fix_3d(ref%edge_blk), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%start_index"

    call compare_int32_1d_array( &
        actual=actual_rich%start_index, &
        ref=logical_fix_1d(ref%start_index), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%end_index"

    call compare_int32_1d_array( &
        actual=actual_rich%end_index, &
        ref=logical_fix_1d(ref%end_index), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%start_block"

    call compare_int32_1d_array( &
        actual=actual_rich%start_block, &
        ref=logical_fix_1d(ref%start_block), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%end_block"

    call compare_int32_1d_array( &
        actual=actual_rich%end_block, &
        ref=logical_fix_1d(ref%end_block), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_grid_vertices_struct

  subroutine compare_t_nh_prog_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_nh_prog), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_nh_prog), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%w"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%w, &
        ref=ref%w, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%w, &
        ref=ref%w, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vn"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%vn, &
        ref=ref%vn, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%vn, &
        ref=ref%vn, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_nh_prog_struct

  subroutine compare_int32_1d_array( &
    actual, &
    ref, &
    use_openacc, &
    result, &
    array_expr, &
    rel_threshold, &
    abs_threshold &
  )
    type(c_ptr) :: actual
    integer(kind=c_int), dimension(:), target, intent(in) :: ref
    logical, intent(in) :: use_openacc
    logical, intent(out) :: result
    real(kind=c_double), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in) :: array_expr

    real(kind=c_double) :: actual_rel_threshold, actual_abs_threshold
    logical :: local_result
    integer :: i0
    integer(kind=c_int), dimension(:), pointer :: actual_rich
    CHARACTER(len=5000) :: message_text = ''

    integer(kind=c_int) :: error_ref, error_actual
    integer, dimension(0:0) :: max_threshold_ratio_loc
    real(8) :: rel_error, abs_error, threshold_ratio

    integer :: max_threshold_ratio_i0
    integer :: first_fail_i0
    integer :: last_fail_i0
    integer :: dim_i0
    integer :: total_fails
    integer :: total_indices
    integer :: first_fail

    first_fail_i0 = -1
    last_fail_i0 = -1
    dim_i0 = size(ref, dim=1)
    total_fails = 0
    total_indices = 0
    first_fail = -1

    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not!"
        print *, "compare_int32_1d_array"
        print *, trim(message_text)
      end if

      return
    end if

    result = .true.

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = int32_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = int32_default_abs_threshold
    end if

    call c_f_pointer(actual, actual_rich, shape=shape(ref))

    if (use_openacc .eqv. .false.) then

    do i0 = 1, size(ref, dim=1)

    call compare_int32_scalar( &
      actual=actual_rich(i0), &
      ref=transfer(ref(i0), mold=int(1, kind=4)), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

      result = result .and. local_result
    end do



    do i0 = 1, size(ref, dim=1)
       
    call compare_int32_scalar( &
      actual=actual_rich(i0), &
      ref=transfer(ref(i0), mold=int(1, kind=4)), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

       if (.not. local_result) then
           if (first_fail == -1) then
               first_fail_i0 = i0
               first_fail = 1
           endif
           last_fail_i0 = i0
           total_fails = total_fails + 1
       endif
    end do


      total_indices =  size(ref, dim=1)

      if (.not. result) then
        max_threshold_ratio_loc = maxloc(abs(ref - actual_rich) / max(actual_rel_threshold * abs(ref), actual_abs_threshold))
        max_threshold_ratio_i0 = max_threshold_ratio_loc(0)

        error_ref = ref(max_threshold_ratio_i0)
        error_actual = actual_rich(max_threshold_ratio_i0)

        threshold_ratio = real(abs(error_ref - error_actual), kind=8) / max(actual_rel_threshold * abs(error_ref), actual_abs_threshold)
        rel_error = abs(real(error_ref - error_actual, kind=8)/error_ref)
        abs_error = abs(error_ref - error_actual)

        write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,"(",i0,")",a,e28.20,a,e28.20,a,"(",i0,")",a,"(",i0,")",a,i0,a,i0,a,i0,a,"(",i0,")")') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - max_threshold_ratio = ", &
            threshold_ratio, &
            ", rel_error = ", &
            rel_error, &
            ", abs_error = ", &
            abs_error, &
          char(10)//"    - at (", &
            max_threshold_ratio_i0, &
            "), ref = ", &
            error_ref, &
            ", actual = ", &
            error_actual, &
            char(10)//"    - first_fail_index: ", first_fail_i0, &
            " last_fail_index: ", last_fail_i0, &
            " total_fails: ", total_fails, &
            " total_indices: ", total_indices, &
            " call_to_size: ", size(ref), &
            " shape: ", dim_i0
        print *, "compare_int32_1d_array"
        print *, trim(message_text)

      end if

      call free(actual)

    else
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#else

    rel_error = 0
    abs_error = 0

    !$ACC PARALLEL &
    !$ACC   DEFAULT(PRESENT) &
    !$ACC   DEVICEPTR(actual_rich) &
    !$ACC   REDUCTION(.AND.:result) &
    !$ACC   REDUCTION(MAX:rel_error) &
    !$ACC   REDUCTION(MAX:abs_error)
    !$ACC LOOP GANG VECTOR COLLAPSE(1)

    do i0 = 1, size(ref, dim=1)
      result = abs(ref(i0) - actual_rich(i0)) <= max(actual_rel_threshold * abs(ref(i0)), actual_abs_threshold)

      if (.not. result) then
        rel_error = abs(real(ref(i0) - actual_rich(i0), kind=8)/ref(i0))
        abs_error = abs(ref(i0) - actual_rich(i0))
      else
        rel_error = 0
        abs_error = 0
      end if
    end do

    !$ACC END PARALLEL

    if (.not. result) then
      write (message_text, '(a,a,a,e28.20,a,e28.20,a,i0,a,"(",i0,")")') &
        "Verification failed for array '", &
          trim(array_expr), &
        "':"//char(10)//"   - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
        char(10)//"   - call_to_size: ", size(ref), &
          " shape: ", dim_i0
      print *, "compare_float64_3d_array"
      print *, trim(message_text)
    end if

    call c_acc_free(actual)
#endif
    end if

  end subroutine compare_int32_1d_array

  subroutine compare_float64_3d_array( &
    actual, &
    ref, &
    use_openacc, &
    result, &
    array_expr, &
    rel_threshold, &
    abs_threshold &
  )
    type(c_ptr) :: actual
    real(kind=c_double), dimension(:,:,:), target, intent(in) :: ref
    logical, intent(in) :: use_openacc
    logical, intent(out) :: result
    real(kind=c_double), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in) :: array_expr

    real(kind=c_double) :: actual_rel_threshold, actual_abs_threshold
    logical :: local_result
    integer :: i0, i1, i2
    real(kind=c_double), dimension(:,:,:), pointer :: actual_rich
    CHARACTER(len=5000) :: message_text = ''

    real(kind=c_double) :: error_ref, error_actual
    integer, dimension(0:2) :: max_threshold_ratio_loc
    real(8) :: rel_error, abs_error, threshold_ratio

    integer :: max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2
    integer :: first_fail_i0, first_fail_i1, first_fail_i2
    integer :: last_fail_i0, last_fail_i1, last_fail_i2
    integer :: dim_i0, dim_i1, dim_i2
    integer :: total_fails
    integer :: total_indices
    integer :: first_fail

    first_fail_i0 = -1
    first_fail_i1 = -1
    first_fail_i2 = -1
    last_fail_i0 = -1
    last_fail_i1 = -1
    last_fail_i2 = -1
    dim_i0 = size(ref, dim=1)
    dim_i1 = size(ref, dim=2)
    dim_i2 = size(ref, dim=3)
    total_fails = 0
    total_indices = 0
    first_fail = -1

    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not!"
        print *, "compare_float64_3d_array"
        print *, trim(message_text)
      end if

      return
    end if

    result = .true.

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = float64_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = float64_default_abs_threshold
    end if

    call c_f_pointer(actual, actual_rich, shape=shape(ref))

    if (use_openacc .eqv. .false.) then

    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
        do i2 = 1, size(ref, dim=3)

    call compare_float64_scalar( &
      actual=actual_rich(i0, i1, i2), &
      ref=ref(i0, i1, i2), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

      result = result .and. local_result
        end do
      end do
    end do



    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
        do i2 = 1, size(ref, dim=3)
       
    call compare_float64_scalar( &
      actual=actual_rich(i0, i1, i2), &
      ref=ref(i0, i1, i2), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

       if (.not. local_result) then
           if (first_fail == -1) then
               first_fail_i0 = i0
            first_fail_i1 = i1
            first_fail_i2 = i2
               first_fail = 1
           endif
           last_fail_i0 = i0
        last_fail_i1 = i1
        last_fail_i2 = i2
           total_fails = total_fails + 1
       endif
        end do
      end do
    end do


      total_indices =  size(ref, dim=1) * size(ref, dim=2) * size(ref, dim=3)

      if (.not. result) then
        max_threshold_ratio_loc = maxloc(abs(ref - actual_rich) / max(actual_rel_threshold * abs(ref), actual_abs_threshold))
        max_threshold_ratio_i0 = max_threshold_ratio_loc(0); max_threshold_ratio_i1 = max_threshold_ratio_loc(1); max_threshold_ratio_i2 = max_threshold_ratio_loc(2)

        error_ref = ref(max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2)
        error_actual = actual_rich(max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2)

        threshold_ratio = real(abs(error_ref - error_actual), kind=8) / max(actual_rel_threshold * abs(error_ref), actual_abs_threshold)
        rel_error = abs(real(error_ref - error_actual, kind=8)/error_ref)
        abs_error = abs(error_ref - error_actual)

        write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,"(",i0,", ",i0,", ",i0,")",a,e28.20,a,e28.20,a,"(",i0,", ",i0,", ",i0,")",a,"(",i0,", ",i0,", ",i0,")",a,i0,a,i0,a,i0,a,"(",i0,", ",i0,", ",i0,")")') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - max_threshold_ratio = ", &
            threshold_ratio, &
            ", rel_error = ", &
            rel_error, &
            ", abs_error = ", &
            abs_error, &
          char(10)//"    - at (", &
            max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2, &
            "), ref = ", &
            error_ref, &
            ", actual = ", &
            error_actual, &
            char(10)//"    - first_fail_index: ", first_fail_i0, first_fail_i1, first_fail_i2, &
            " last_fail_index: ", last_fail_i0, last_fail_i1, last_fail_i2, &
            " total_fails: ", total_fails, &
            " total_indices: ", total_indices, &
            " call_to_size: ", size(ref), &
            " shape: ", dim_i0, dim_i1, dim_i2
        print *, "compare_float64_3d_array"
        print *, trim(message_text)

      end if

      call free(actual)

    else
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#else

    rel_error = 0
    abs_error = 0

    !$ACC PARALLEL &
    !$ACC   DEFAULT(PRESENT) &
    !$ACC   DEVICEPTR(actual_rich) &
    !$ACC   REDUCTION(.AND.:result) &
    !$ACC   REDUCTION(MAX:rel_error) &
    !$ACC   REDUCTION(MAX:abs_error)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)

    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
        do i2 = 1, size(ref, dim=3)
      result = abs(ref(i0, i1, i2) - actual_rich(i0, i1, i2)) <= max(actual_rel_threshold * abs(ref(i0, i1, i2)), actual_abs_threshold)

      if (.not. result) then
        rel_error = abs(real(ref(i0, i1, i2) - actual_rich(i0, i1, i2), kind=8)/ref(i0, i1, i2))
        abs_error = abs(ref(i0, i1, i2) - actual_rich(i0, i1, i2))
      else
        rel_error = 0
        abs_error = 0
      end if
        end do
      end do
    end do

    !$ACC END PARALLEL

    if (.not. result) then
      write (message_text, '(a,a,a,e28.20,a,e28.20,a,i0,a,"(",i0,", ",i0,", ",i0,")")') &
        "Verification failed for array '", &
          trim(array_expr), &
        "':"//char(10)//"   - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
        char(10)//"   - call_to_size: ", size(ref), &
          " shape: ", dim_i0, dim_i1, dim_i2
      print *, "compare_float64_3d_array"
      print *, trim(message_text)
    end if

    call c_acc_free(actual)
#endif
    end if

  end subroutine compare_float64_3d_array

  subroutine compare_float64_4d_array( &
    actual, &
    ref, &
    use_openacc, &
    result, &
    array_expr, &
    rel_threshold, &
    abs_threshold &
  )
    type(c_ptr) :: actual
    real(kind=c_double), dimension(:,:,:,:), target, intent(in) :: ref
    logical, intent(in) :: use_openacc
    logical, intent(out) :: result
    real(kind=c_double), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in) :: array_expr

    real(kind=c_double) :: actual_rel_threshold, actual_abs_threshold
    logical :: local_result
    integer :: i0, i1, i2, i3
    real(kind=c_double), dimension(:,:,:,:), pointer :: actual_rich
    CHARACTER(len=5000) :: message_text = ''

    real(kind=c_double) :: error_ref, error_actual
    integer, dimension(0:3) :: max_threshold_ratio_loc
    real(8) :: rel_error, abs_error, threshold_ratio

    integer :: max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2, max_threshold_ratio_i3
    integer :: first_fail_i0, first_fail_i1, first_fail_i2, first_fail_i3
    integer :: last_fail_i0, last_fail_i1, last_fail_i2, last_fail_i3
    integer :: dim_i0, dim_i1, dim_i2, dim_i3
    integer :: total_fails
    integer :: total_indices
    integer :: first_fail

    first_fail_i0 = -1
    first_fail_i1 = -1
    first_fail_i2 = -1
    first_fail_i3 = -1
    last_fail_i0 = -1
    last_fail_i1 = -1
    last_fail_i2 = -1
    last_fail_i3 = -1
    dim_i0 = size(ref, dim=1)
    dim_i1 = size(ref, dim=2)
    dim_i2 = size(ref, dim=3)
    dim_i3 = size(ref, dim=4)
    total_fails = 0
    total_indices = 0
    first_fail = -1

    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not!"
        print *, "compare_float64_4d_array"
        print *, trim(message_text)
      end if

      return
    end if

    result = .true.

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = float64_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = float64_default_abs_threshold
    end if

    call c_f_pointer(actual, actual_rich, shape=shape(ref))

    if (use_openacc .eqv. .false.) then

    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
        do i2 = 1, size(ref, dim=3)
          do i3 = 1, size(ref, dim=4)

    call compare_float64_scalar( &
      actual=actual_rich(i0, i1, i2, i3), &
      ref=ref(i0, i1, i2, i3), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

      result = result .and. local_result
          end do
        end do
      end do
    end do



    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
        do i2 = 1, size(ref, dim=3)
          do i3 = 1, size(ref, dim=4)
       
    call compare_float64_scalar( &
      actual=actual_rich(i0, i1, i2, i3), &
      ref=ref(i0, i1, i2, i3), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

       if (.not. local_result) then
           if (first_fail == -1) then
               first_fail_i0 = i0
            first_fail_i1 = i1
            first_fail_i2 = i2
            first_fail_i3 = i3
               first_fail = 1
           endif
           last_fail_i0 = i0
        last_fail_i1 = i1
        last_fail_i2 = i2
        last_fail_i3 = i3
           total_fails = total_fails + 1
       endif
          end do
        end do
      end do
    end do


      total_indices =  size(ref, dim=1) * size(ref, dim=2) * size(ref, dim=3) * size(ref, dim=4)

      if (.not. result) then
        max_threshold_ratio_loc = maxloc(abs(ref - actual_rich) / max(actual_rel_threshold * abs(ref), actual_abs_threshold))
        max_threshold_ratio_i0 = max_threshold_ratio_loc(0); max_threshold_ratio_i1 = max_threshold_ratio_loc(1); max_threshold_ratio_i2 = max_threshold_ratio_loc(2); max_threshold_ratio_i3 = max_threshold_ratio_loc(3)

        error_ref = ref(max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2, max_threshold_ratio_i3)
        error_actual = actual_rich(max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2, max_threshold_ratio_i3)

        threshold_ratio = real(abs(error_ref - error_actual), kind=8) / max(actual_rel_threshold * abs(error_ref), actual_abs_threshold)
        rel_error = abs(real(error_ref - error_actual, kind=8)/error_ref)
        abs_error = abs(error_ref - error_actual)

        write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,"(",i0,", ",i0,", ",i0,", ",i0,")",a,e28.20,a,e28.20,a,"(",i0,", ",i0,", ",i0,", ",i0,")",a,"(",i0,", ",i0,", ",i0,", ",i0,")",a,i0,a,i0,a,i0,a,"(",i0,", ",i0,", ",i0,", ",i0,")")') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - max_threshold_ratio = ", &
            threshold_ratio, &
            ", rel_error = ", &
            rel_error, &
            ", abs_error = ", &
            abs_error, &
          char(10)//"    - at (", &
            max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2, max_threshold_ratio_i3, &
            "), ref = ", &
            error_ref, &
            ", actual = ", &
            error_actual, &
            char(10)//"    - first_fail_index: ", first_fail_i0, first_fail_i1, first_fail_i2, first_fail_i3, &
            " last_fail_index: ", last_fail_i0, last_fail_i1, last_fail_i2, last_fail_i3, &
            " total_fails: ", total_fails, &
            " total_indices: ", total_indices, &
            " call_to_size: ", size(ref), &
            " shape: ", dim_i0, dim_i1, dim_i2, dim_i3
        print *, "compare_float64_4d_array"
        print *, trim(message_text)

      end if

      call free(actual)

    else
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#else

    rel_error = 0
    abs_error = 0

    !$ACC PARALLEL &
    !$ACC   DEFAULT(PRESENT) &
    !$ACC   DEVICEPTR(actual_rich) &
    !$ACC   REDUCTION(.AND.:result) &
    !$ACC   REDUCTION(MAX:rel_error) &
    !$ACC   REDUCTION(MAX:abs_error)
    !$ACC LOOP GANG VECTOR COLLAPSE(4)

    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
        do i2 = 1, size(ref, dim=3)
          do i3 = 1, size(ref, dim=4)
      result = abs(ref(i0, i1, i2, i3) - actual_rich(i0, i1, i2, i3)) <= max(actual_rel_threshold * abs(ref(i0, i1, i2, i3)), actual_abs_threshold)

      if (.not. result) then
        rel_error = abs(real(ref(i0, i1, i2, i3) - actual_rich(i0, i1, i2, i3), kind=8)/ref(i0, i1, i2, i3))
        abs_error = abs(ref(i0, i1, i2, i3) - actual_rich(i0, i1, i2, i3))
      else
        rel_error = 0
        abs_error = 0
      end if
          end do
        end do
      end do
    end do

    !$ACC END PARALLEL

    if (.not. result) then
      write (message_text, '(a,a,a,e28.20,a,e28.20,a,i0,a,"(",i0,", ",i0,", ",i0,", ",i0,")")') &
        "Verification failed for array '", &
          trim(array_expr), &
        "':"//char(10)//"   - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
        char(10)//"   - call_to_size: ", size(ref), &
          " shape: ", dim_i0, dim_i1, dim_i2, dim_i3
      print *, "compare_float64_3d_array"
      print *, trim(message_text)
    end if

    call c_acc_free(actual)
#endif
    end if

  end subroutine compare_float64_4d_array

  subroutine compare_float64_1d_array( &
    actual, &
    ref, &
    use_openacc, &
    result, &
    array_expr, &
    rel_threshold, &
    abs_threshold &
  )
    type(c_ptr) :: actual
    real(kind=c_double), dimension(:), target, intent(in) :: ref
    logical, intent(in) :: use_openacc
    logical, intent(out) :: result
    real(kind=c_double), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in) :: array_expr

    real(kind=c_double) :: actual_rel_threshold, actual_abs_threshold
    logical :: local_result
    integer :: i0
    real(kind=c_double), dimension(:), pointer :: actual_rich
    CHARACTER(len=5000) :: message_text = ''

    real(kind=c_double) :: error_ref, error_actual
    integer, dimension(0:0) :: max_threshold_ratio_loc
    real(8) :: rel_error, abs_error, threshold_ratio

    integer :: max_threshold_ratio_i0
    integer :: first_fail_i0
    integer :: last_fail_i0
    integer :: dim_i0
    integer :: total_fails
    integer :: total_indices
    integer :: first_fail

    first_fail_i0 = -1
    last_fail_i0 = -1
    dim_i0 = size(ref, dim=1)
    total_fails = 0
    total_indices = 0
    first_fail = -1

    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not!"
        print *, "compare_float64_1d_array"
        print *, trim(message_text)
      end if

      return
    end if

    result = .true.

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = float64_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = float64_default_abs_threshold
    end if

    call c_f_pointer(actual, actual_rich, shape=shape(ref))

    if (use_openacc .eqv. .false.) then

    do i0 = 1, size(ref, dim=1)

    call compare_float64_scalar( &
      actual=actual_rich(i0), &
      ref=ref(i0), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

      result = result .and. local_result
    end do



    do i0 = 1, size(ref, dim=1)
       
    call compare_float64_scalar( &
      actual=actual_rich(i0), &
      ref=ref(i0), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

       if (.not. local_result) then
           if (first_fail == -1) then
               first_fail_i0 = i0
               first_fail = 1
           endif
           last_fail_i0 = i0
           total_fails = total_fails + 1
       endif
    end do


      total_indices =  size(ref, dim=1)

      if (.not. result) then
        max_threshold_ratio_loc = maxloc(abs(ref - actual_rich) / max(actual_rel_threshold * abs(ref), actual_abs_threshold))
        max_threshold_ratio_i0 = max_threshold_ratio_loc(0)

        error_ref = ref(max_threshold_ratio_i0)
        error_actual = actual_rich(max_threshold_ratio_i0)

        threshold_ratio = real(abs(error_ref - error_actual), kind=8) / max(actual_rel_threshold * abs(error_ref), actual_abs_threshold)
        rel_error = abs(real(error_ref - error_actual, kind=8)/error_ref)
        abs_error = abs(error_ref - error_actual)

        write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,"(",i0,")",a,e28.20,a,e28.20,a,"(",i0,")",a,"(",i0,")",a,i0,a,i0,a,i0,a,"(",i0,")")') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - max_threshold_ratio = ", &
            threshold_ratio, &
            ", rel_error = ", &
            rel_error, &
            ", abs_error = ", &
            abs_error, &
          char(10)//"    - at (", &
            max_threshold_ratio_i0, &
            "), ref = ", &
            error_ref, &
            ", actual = ", &
            error_actual, &
            char(10)//"    - first_fail_index: ", first_fail_i0, &
            " last_fail_index: ", last_fail_i0, &
            " total_fails: ", total_fails, &
            " total_indices: ", total_indices, &
            " call_to_size: ", size(ref), &
            " shape: ", dim_i0
        print *, "compare_float64_1d_array"
        print *, trim(message_text)

      end if

      call free(actual)

    else
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#else

    rel_error = 0
    abs_error = 0

    !$ACC PARALLEL &
    !$ACC   DEFAULT(PRESENT) &
    !$ACC   DEVICEPTR(actual_rich) &
    !$ACC   REDUCTION(.AND.:result) &
    !$ACC   REDUCTION(MAX:rel_error) &
    !$ACC   REDUCTION(MAX:abs_error)
    !$ACC LOOP GANG VECTOR COLLAPSE(1)

    do i0 = 1, size(ref, dim=1)
      result = abs(ref(i0) - actual_rich(i0)) <= max(actual_rel_threshold * abs(ref(i0)), actual_abs_threshold)

      if (.not. result) then
        rel_error = abs(real(ref(i0) - actual_rich(i0), kind=8)/ref(i0))
        abs_error = abs(ref(i0) - actual_rich(i0))
      else
        rel_error = 0
        abs_error = 0
      end if
    end do

    !$ACC END PARALLEL

    if (.not. result) then
      write (message_text, '(a,a,a,e28.20,a,e28.20,a,i0,a,"(",i0,")")') &
        "Verification failed for array '", &
          trim(array_expr), &
        "':"//char(10)//"   - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
        char(10)//"   - call_to_size: ", size(ref), &
          " shape: ", dim_i0
      print *, "compare_float64_3d_array"
      print *, trim(message_text)
    end if

    call c_acc_free(actual)
#endif
    end if

  end subroutine compare_float64_1d_array

  subroutine compare_int32_2d_array( &
    actual, &
    ref, &
    use_openacc, &
    result, &
    array_expr, &
    rel_threshold, &
    abs_threshold &
  )
    type(c_ptr) :: actual
    integer(kind=c_int), dimension(:,:), target, intent(in) :: ref
    logical, intent(in) :: use_openacc
    logical, intent(out) :: result
    real(kind=c_double), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in) :: array_expr

    real(kind=c_double) :: actual_rel_threshold, actual_abs_threshold
    logical :: local_result
    integer :: i0, i1
    integer(kind=c_int), dimension(:,:), pointer :: actual_rich
    CHARACTER(len=5000) :: message_text = ''

    integer(kind=c_int) :: error_ref, error_actual
    integer, dimension(0:1) :: max_threshold_ratio_loc
    real(8) :: rel_error, abs_error, threshold_ratio

    integer :: max_threshold_ratio_i0, max_threshold_ratio_i1
    integer :: first_fail_i0, first_fail_i1
    integer :: last_fail_i0, last_fail_i1
    integer :: dim_i0, dim_i1
    integer :: total_fails
    integer :: total_indices
    integer :: first_fail

    first_fail_i0 = -1
    first_fail_i1 = -1
    last_fail_i0 = -1
    last_fail_i1 = -1
    dim_i0 = size(ref, dim=1)
    dim_i1 = size(ref, dim=2)
    total_fails = 0
    total_indices = 0
    first_fail = -1

    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not!"
        print *, "compare_int32_2d_array"
        print *, trim(message_text)
      end if

      return
    end if

    result = .true.

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = int32_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = int32_default_abs_threshold
    end if

    call c_f_pointer(actual, actual_rich, shape=shape(ref))

    if (use_openacc .eqv. .false.) then

    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)

    call compare_int32_scalar( &
      actual=actual_rich(i0, i1), &
      ref=transfer(ref(i0, i1), mold=int(1, kind=4)), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

      result = result .and. local_result
      end do
    end do



    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
       
    call compare_int32_scalar( &
      actual=actual_rich(i0, i1), &
      ref=transfer(ref(i0, i1), mold=int(1, kind=4)), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

       if (.not. local_result) then
           if (first_fail == -1) then
               first_fail_i0 = i0
            first_fail_i1 = i1
               first_fail = 1
           endif
           last_fail_i0 = i0
        last_fail_i1 = i1
           total_fails = total_fails + 1
       endif
      end do
    end do


      total_indices =  size(ref, dim=1) * size(ref, dim=2)

      if (.not. result) then
        max_threshold_ratio_loc = maxloc(abs(ref - actual_rich) / max(actual_rel_threshold * abs(ref), actual_abs_threshold))
        max_threshold_ratio_i0 = max_threshold_ratio_loc(0); max_threshold_ratio_i1 = max_threshold_ratio_loc(1)

        error_ref = ref(max_threshold_ratio_i0, max_threshold_ratio_i1)
        error_actual = actual_rich(max_threshold_ratio_i0, max_threshold_ratio_i1)

        threshold_ratio = real(abs(error_ref - error_actual), kind=8) / max(actual_rel_threshold * abs(error_ref), actual_abs_threshold)
        rel_error = abs(real(error_ref - error_actual, kind=8)/error_ref)
        abs_error = abs(error_ref - error_actual)

        write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,"(",i0,", ",i0,")",a,e28.20,a,e28.20,a,"(",i0,", ",i0,")",a,"(",i0,", ",i0,")",a,i0,a,i0,a,i0,a,"(",i0,", ",i0,")")') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - max_threshold_ratio = ", &
            threshold_ratio, &
            ", rel_error = ", &
            rel_error, &
            ", abs_error = ", &
            abs_error, &
          char(10)//"    - at (", &
            max_threshold_ratio_i0, max_threshold_ratio_i1, &
            "), ref = ", &
            error_ref, &
            ", actual = ", &
            error_actual, &
            char(10)//"    - first_fail_index: ", first_fail_i0, first_fail_i1, &
            " last_fail_index: ", last_fail_i0, last_fail_i1, &
            " total_fails: ", total_fails, &
            " total_indices: ", total_indices, &
            " call_to_size: ", size(ref), &
            " shape: ", dim_i0, dim_i1
        print *, "compare_int32_2d_array"
        print *, trim(message_text)

      end if

      call free(actual)

    else
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#else

    rel_error = 0
    abs_error = 0

    !$ACC PARALLEL &
    !$ACC   DEFAULT(PRESENT) &
    !$ACC   DEVICEPTR(actual_rich) &
    !$ACC   REDUCTION(.AND.:result) &
    !$ACC   REDUCTION(MAX:rel_error) &
    !$ACC   REDUCTION(MAX:abs_error)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)

    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
      result = abs(ref(i0, i1) - actual_rich(i0, i1)) <= max(actual_rel_threshold * abs(ref(i0, i1)), actual_abs_threshold)

      if (.not. result) then
        rel_error = abs(real(ref(i0, i1) - actual_rich(i0, i1), kind=8)/ref(i0, i1))
        abs_error = abs(ref(i0, i1) - actual_rich(i0, i1))
      else
        rel_error = 0
        abs_error = 0
      end if
      end do
    end do

    !$ACC END PARALLEL

    if (.not. result) then
      write (message_text, '(a,a,a,e28.20,a,e28.20,a,i0,a,"(",i0,", ",i0,")")') &
        "Verification failed for array '", &
          trim(array_expr), &
        "':"//char(10)//"   - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
        char(10)//"   - call_to_size: ", size(ref), &
          " shape: ", dim_i0, dim_i1
      print *, "compare_float64_3d_array"
      print *, trim(message_text)
    end if

    call c_acc_free(actual)
#endif
    end if

  end subroutine compare_int32_2d_array

  subroutine compare_int32_3d_array( &
    actual, &
    ref, &
    use_openacc, &
    result, &
    array_expr, &
    rel_threshold, &
    abs_threshold &
  )
    type(c_ptr) :: actual
    integer(kind=c_int), dimension(:,:,:), target, intent(in) :: ref
    logical, intent(in) :: use_openacc
    logical, intent(out) :: result
    real(kind=c_double), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in) :: array_expr

    real(kind=c_double) :: actual_rel_threshold, actual_abs_threshold
    logical :: local_result
    integer :: i0, i1, i2
    integer(kind=c_int), dimension(:,:,:), pointer :: actual_rich
    CHARACTER(len=5000) :: message_text = ''

    integer(kind=c_int) :: error_ref, error_actual
    integer, dimension(0:2) :: max_threshold_ratio_loc
    real(8) :: rel_error, abs_error, threshold_ratio

    integer :: max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2
    integer :: first_fail_i0, first_fail_i1, first_fail_i2
    integer :: last_fail_i0, last_fail_i1, last_fail_i2
    integer :: dim_i0, dim_i1, dim_i2
    integer :: total_fails
    integer :: total_indices
    integer :: first_fail

    first_fail_i0 = -1
    first_fail_i1 = -1
    first_fail_i2 = -1
    last_fail_i0 = -1
    last_fail_i1 = -1
    last_fail_i2 = -1
    dim_i0 = size(ref, dim=1)
    dim_i1 = size(ref, dim=2)
    dim_i2 = size(ref, dim=3)
    total_fails = 0
    total_indices = 0
    first_fail = -1

    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not!"
        print *, "compare_int32_3d_array"
        print *, trim(message_text)
      end if

      return
    end if

    result = .true.

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = int32_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = int32_default_abs_threshold
    end if

    call c_f_pointer(actual, actual_rich, shape=shape(ref))

    if (use_openacc .eqv. .false.) then

    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
        do i2 = 1, size(ref, dim=3)

    call compare_int32_scalar( &
      actual=actual_rich(i0, i1, i2), &
      ref=transfer(ref(i0, i1, i2), mold=int(1, kind=4)), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

      result = result .and. local_result
        end do
      end do
    end do



    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
        do i2 = 1, size(ref, dim=3)
       
    call compare_int32_scalar( &
      actual=actual_rich(i0, i1, i2), &
      ref=transfer(ref(i0, i1, i2), mold=int(1, kind=4)), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

       if (.not. local_result) then
           if (first_fail == -1) then
               first_fail_i0 = i0
            first_fail_i1 = i1
            first_fail_i2 = i2
               first_fail = 1
           endif
           last_fail_i0 = i0
        last_fail_i1 = i1
        last_fail_i2 = i2
           total_fails = total_fails + 1
       endif
        end do
      end do
    end do


      total_indices =  size(ref, dim=1) * size(ref, dim=2) * size(ref, dim=3)

      if (.not. result) then
        max_threshold_ratio_loc = maxloc(abs(ref - actual_rich) / max(actual_rel_threshold * abs(ref), actual_abs_threshold))
        max_threshold_ratio_i0 = max_threshold_ratio_loc(0); max_threshold_ratio_i1 = max_threshold_ratio_loc(1); max_threshold_ratio_i2 = max_threshold_ratio_loc(2)

        error_ref = ref(max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2)
        error_actual = actual_rich(max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2)

        threshold_ratio = real(abs(error_ref - error_actual), kind=8) / max(actual_rel_threshold * abs(error_ref), actual_abs_threshold)
        rel_error = abs(real(error_ref - error_actual, kind=8)/error_ref)
        abs_error = abs(error_ref - error_actual)

        write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,"(",i0,", ",i0,", ",i0,")",a,e28.20,a,e28.20,a,"(",i0,", ",i0,", ",i0,")",a,"(",i0,", ",i0,", ",i0,")",a,i0,a,i0,a,i0,a,"(",i0,", ",i0,", ",i0,")")') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - max_threshold_ratio = ", &
            threshold_ratio, &
            ", rel_error = ", &
            rel_error, &
            ", abs_error = ", &
            abs_error, &
          char(10)//"    - at (", &
            max_threshold_ratio_i0, max_threshold_ratio_i1, max_threshold_ratio_i2, &
            "), ref = ", &
            error_ref, &
            ", actual = ", &
            error_actual, &
            char(10)//"    - first_fail_index: ", first_fail_i0, first_fail_i1, first_fail_i2, &
            " last_fail_index: ", last_fail_i0, last_fail_i1, last_fail_i2, &
            " total_fails: ", total_fails, &
            " total_indices: ", total_indices, &
            " call_to_size: ", size(ref), &
            " shape: ", dim_i0, dim_i1, dim_i2
        print *, "compare_int32_3d_array"
        print *, trim(message_text)

      end if

      call free(actual)

    else
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#else

    rel_error = 0
    abs_error = 0

    !$ACC PARALLEL &
    !$ACC   DEFAULT(PRESENT) &
    !$ACC   DEVICEPTR(actual_rich) &
    !$ACC   REDUCTION(.AND.:result) &
    !$ACC   REDUCTION(MAX:rel_error) &
    !$ACC   REDUCTION(MAX:abs_error)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)

    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
        do i2 = 1, size(ref, dim=3)
      result = abs(ref(i0, i1, i2) - actual_rich(i0, i1, i2)) <= max(actual_rel_threshold * abs(ref(i0, i1, i2)), actual_abs_threshold)

      if (.not. result) then
        rel_error = abs(real(ref(i0, i1, i2) - actual_rich(i0, i1, i2), kind=8)/ref(i0, i1, i2))
        abs_error = abs(ref(i0, i1, i2) - actual_rich(i0, i1, i2))
      else
        rel_error = 0
        abs_error = 0
      end if
        end do
      end do
    end do

    !$ACC END PARALLEL

    if (.not. result) then
      write (message_text, '(a,a,a,e28.20,a,e28.20,a,i0,a,"(",i0,", ",i0,", ",i0,")")') &
        "Verification failed for array '", &
          trim(array_expr), &
        "':"//char(10)//"   - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
        char(10)//"   - call_to_size: ", size(ref), &
          " shape: ", dim_i0, dim_i1, dim_i2
      print *, "compare_float64_3d_array"
      print *, trim(message_text)
    end if

    call c_acc_free(actual)
#endif
    end if

  end subroutine compare_int32_3d_array

  subroutine compare_float64_2d_array( &
    actual, &
    ref, &
    use_openacc, &
    result, &
    array_expr, &
    rel_threshold, &
    abs_threshold &
  )
    type(c_ptr) :: actual
    real(kind=c_double), dimension(:,:), target, intent(in) :: ref
    logical, intent(in) :: use_openacc
    logical, intent(out) :: result
    real(kind=c_double), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in) :: array_expr

    real(kind=c_double) :: actual_rel_threshold, actual_abs_threshold
    logical :: local_result
    integer :: i0, i1
    real(kind=c_double), dimension(:,:), pointer :: actual_rich
    CHARACTER(len=5000) :: message_text = ''

    real(kind=c_double) :: error_ref, error_actual
    integer, dimension(0:1) :: max_threshold_ratio_loc
    real(8) :: rel_error, abs_error, threshold_ratio

    integer :: max_threshold_ratio_i0, max_threshold_ratio_i1
    integer :: first_fail_i0, first_fail_i1
    integer :: last_fail_i0, last_fail_i1
    integer :: dim_i0, dim_i1
    integer :: total_fails
    integer :: total_indices
    integer :: first_fail

    first_fail_i0 = -1
    first_fail_i1 = -1
    last_fail_i0 = -1
    last_fail_i1 = -1
    dim_i0 = size(ref, dim=1)
    dim_i1 = size(ref, dim=2)
    total_fails = 0
    total_indices = 0
    first_fail = -1

    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not!"
        print *, "compare_float64_2d_array"
        print *, trim(message_text)
      end if

      return
    end if

    result = .true.

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = float64_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = float64_default_abs_threshold
    end if

    call c_f_pointer(actual, actual_rich, shape=shape(ref))

    if (use_openacc .eqv. .false.) then

    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)

    call compare_float64_scalar( &
      actual=actual_rich(i0, i1), &
      ref=ref(i0, i1), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

      result = result .and. local_result
      end do
    end do



    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
       
    call compare_float64_scalar( &
      actual=actual_rich(i0, i1), &
      ref=ref(i0, i1), &
      result=local_result, &
      rel_threshold=actual_rel_threshold, &
      abs_threshold=actual_abs_threshold &
    )

       if (.not. local_result) then
           if (first_fail == -1) then
               first_fail_i0 = i0
            first_fail_i1 = i1
               first_fail = 1
           endif
           last_fail_i0 = i0
        last_fail_i1 = i1
           total_fails = total_fails + 1
       endif
      end do
    end do


      total_indices =  size(ref, dim=1) * size(ref, dim=2)

      if (.not. result) then
        max_threshold_ratio_loc = maxloc(abs(ref - actual_rich) / max(actual_rel_threshold * abs(ref), actual_abs_threshold))
        max_threshold_ratio_i0 = max_threshold_ratio_loc(0); max_threshold_ratio_i1 = max_threshold_ratio_loc(1)

        error_ref = ref(max_threshold_ratio_i0, max_threshold_ratio_i1)
        error_actual = actual_rich(max_threshold_ratio_i0, max_threshold_ratio_i1)

        threshold_ratio = real(abs(error_ref - error_actual), kind=8) / max(actual_rel_threshold * abs(error_ref), actual_abs_threshold)
        rel_error = abs(real(error_ref - error_actual, kind=8)/error_ref)
        abs_error = abs(error_ref - error_actual)

        write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,"(",i0,", ",i0,")",a,e28.20,a,e28.20,a,"(",i0,", ",i0,")",a,"(",i0,", ",i0,")",a,i0,a,i0,a,i0,a,"(",i0,", ",i0,")")') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - max_threshold_ratio = ", &
            threshold_ratio, &
            ", rel_error = ", &
            rel_error, &
            ", abs_error = ", &
            abs_error, &
          char(10)//"    - at (", &
            max_threshold_ratio_i0, max_threshold_ratio_i1, &
            "), ref = ", &
            error_ref, &
            ", actual = ", &
            error_actual, &
            char(10)//"    - first_fail_index: ", first_fail_i0, first_fail_i1, &
            " last_fail_index: ", last_fail_i0, last_fail_i1, &
            " total_fails: ", total_fails, &
            " total_indices: ", total_indices, &
            " call_to_size: ", size(ref), &
            " shape: ", dim_i0, dim_i1
        print *, "compare_float64_2d_array"
        print *, trim(message_text)

      end if

      call free(actual)

    else
#ifndef _OPENACC
      print *, "!!!ERROR!!! Requested OpenACC, but built without OpenACC (SDFG bindings file)"
      return
#else

    rel_error = 0
    abs_error = 0

    !$ACC PARALLEL &
    !$ACC   DEFAULT(PRESENT) &
    !$ACC   DEVICEPTR(actual_rich) &
    !$ACC   REDUCTION(.AND.:result) &
    !$ACC   REDUCTION(MAX:rel_error) &
    !$ACC   REDUCTION(MAX:abs_error)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)

    do i0 = 1, size(ref, dim=1)
      do i1 = 1, size(ref, dim=2)
      result = abs(ref(i0, i1) - actual_rich(i0, i1)) <= max(actual_rel_threshold * abs(ref(i0, i1)), actual_abs_threshold)

      if (.not. result) then
        rel_error = abs(real(ref(i0, i1) - actual_rich(i0, i1), kind=8)/ref(i0, i1))
        abs_error = abs(ref(i0, i1) - actual_rich(i0, i1))
      else
        rel_error = 0
        abs_error = 0
      end if
      end do
    end do

    !$ACC END PARALLEL

    if (.not. result) then
      write (message_text, '(a,a,a,e28.20,a,e28.20,a,i0,a,"(",i0,", ",i0,")")') &
        "Verification failed for array '", &
          trim(array_expr), &
        "':"//char(10)//"   - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
        char(10)//"   - call_to_size: ", size(ref), &
          " shape: ", dim_i0, dim_i1
      print *, "compare_float64_3d_array"
      print *, trim(message_text)
    end if

    call c_acc_free(actual)
#endif
    end if

  end subroutine compare_float64_2d_array



  subroutine run_velocity_no_nproma_if_prop_lvn_only_1_istep_2( &
    p_diag, &
    p_int, &
    p_metrics, &
    p_patch, &
    p_prog, &
    z_kin_hor_e, &
    z_vt_ie, &
    z_w_concorr_me, &
    dt_linintp_ubc, &
    dtime, &
    istep, &
    ldeepatmo, &
    lvn_only, &
    ntnd &
  )
    type(t_nh_diag), target :: p_diag
    type(t_int_state), target :: p_int
    type(t_nh_metrics), target :: p_metrics
    type(t_patch), target :: p_patch
    type(t_nh_prog), target :: p_prog
    real(kind=c_double), dimension(:,:,:), target :: z_kin_hor_e
    real(kind=c_double), dimension(:,:,:), target :: z_vt_ie
    real(kind=c_double), dimension(:,:,:), target :: z_w_concorr_me
    real(kind=c_double) :: dt_linintp_ubc
    real(kind=c_double) :: dtime
    integer(kind=c_int) :: istep
    integer(kind=c_int) :: ldeepatmo
    integer(kind=c_int) :: lvn_only
    integer(kind=c_int) :: ntnd





    copy_or_ptr_global_data = copy_in_global_data_type( &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_diag = copy_in_t_nh_diag( &
    fortran_obj=p_diag, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_int = copy_in_t_int_state( &
    fortran_obj=p_int, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_metrics = copy_in_t_nh_metrics( &
    fortran_obj=p_metrics, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_patch = copy_in_t_patch( &
    fortran_obj=p_patch, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_prog = copy_in_t_nh_prog( &
    fortran_obj=p_prog, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
#ifndef _OPENACC
    copy_or_ptr_z_kin_hor_e = copy_in_float64_3d_array( &
    fortran_array=z_kin_hor_e, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_kin_hor_e = copy_in_float64_3d_array( &
    fortran_array=z_kin_hor_e, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_vt_ie = copy_in_float64_3d_array( &
    fortran_array=z_vt_ie, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_vt_ie = copy_in_float64_3d_array( &
    fortran_array=z_vt_ie, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_w_concorr_me = copy_in_float64_3d_array( &
    fortran_array=z_w_concorr_me, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_w_concorr_me = copy_in_float64_3d_array( &
    fortran_array=z_w_concorr_me, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif


    if (is_initialized .eqv. .false.) then
      is_initialized = .true.
      dace_state = dace_init_velocity_no_nproma_if_prop_lvn_only_1_istep_2( &
      global_data = copy_or_ptr_global_data, &
      p_diag = copy_or_ptr_p_diag, &
      p_int = copy_or_ptr_p_int, &
      p_metrics = copy_or_ptr_p_metrics, &
      p_patch = copy_or_ptr_p_patch, &
      p_prog = copy_or_ptr_p_prog, &
      z_kin_hor_e = copy_or_ptr_z_kin_hor_e, &
      z_vt_ie = copy_or_ptr_z_vt_ie, &
      z_w_concorr_me = copy_or_ptr_z_w_concorr_me, &
      dt_linintp_ubc = dt_linintp_ubc, &
      dtime = dtime, &
      istep = istep, &
      ldeepatmo = ldeepatmo, &
      lvn_only = lvn_only, &
      ntnd = ntnd, &
      f2dace_A_z_kin_hor_e_d_0_s_363 = size(z_kin_hor_e, dim=1), &
      f2dace_A_z_kin_hor_e_d_1_s_364 = size(z_kin_hor_e, dim=2), &
      f2dace_A_z_kin_hor_e_d_2_s_365 = size(z_kin_hor_e, dim=3), &
      f2dace_A_z_vt_ie_d_0_s_366 = size(z_vt_ie, dim=1), &
      f2dace_A_z_vt_ie_d_1_s_367 = size(z_vt_ie, dim=2), &
      f2dace_A_z_vt_ie_d_2_s_368 = size(z_vt_ie, dim=3), &
      f2dace_A_z_w_concorr_me_d_0_s_360 = size(z_w_concorr_me, dim=1), &
      f2dace_A_z_w_concorr_me_d_1_s_361 = size(z_w_concorr_me, dim=2), &
      f2dace_A_z_w_concorr_me_d_2_s_362 = size(z_w_concorr_me, dim=3), &
      f2dace_OA_z_kin_hor_e_d_0_s_363 = lbound(z_kin_hor_e, dim=1), &
      f2dace_OA_z_kin_hor_e_d_1_s_364 = lbound(z_kin_hor_e, dim=2), &
      f2dace_OA_z_kin_hor_e_d_2_s_365 = lbound(z_kin_hor_e, dim=3) &
    )
    end if

    call dace_program_velocity_no_nproma_if_prop_lvn_only_1_istep_2( &
      state = dace_state, &
      global_data = copy_or_ptr_global_data, &
      p_diag = copy_or_ptr_p_diag, &
      p_int = copy_or_ptr_p_int, &
      p_metrics = copy_or_ptr_p_metrics, &
      p_patch = copy_or_ptr_p_patch, &
      p_prog = copy_or_ptr_p_prog, &
      z_kin_hor_e = copy_or_ptr_z_kin_hor_e, &
      z_vt_ie = copy_or_ptr_z_vt_ie, &
      z_w_concorr_me = copy_or_ptr_z_w_concorr_me, &
      dt_linintp_ubc = dt_linintp_ubc, &
      dtime = dtime, &
      istep = istep, &
      ldeepatmo = ldeepatmo, &
      lvn_only = lvn_only, &
      ntnd = ntnd, &
      f2dace_A_z_kin_hor_e_d_0_s_363 = size(z_kin_hor_e, dim=1), &
      f2dace_A_z_kin_hor_e_d_1_s_364 = size(z_kin_hor_e, dim=2), &
      f2dace_A_z_kin_hor_e_d_2_s_365 = size(z_kin_hor_e, dim=3), &
      f2dace_A_z_vt_ie_d_0_s_366 = size(z_vt_ie, dim=1), &
      f2dace_A_z_vt_ie_d_1_s_367 = size(z_vt_ie, dim=2), &
      f2dace_A_z_vt_ie_d_2_s_368 = size(z_vt_ie, dim=3), &
      f2dace_A_z_w_concorr_me_d_0_s_360 = size(z_w_concorr_me, dim=1), &
      f2dace_A_z_w_concorr_me_d_1_s_361 = size(z_w_concorr_me, dim=2), &
      f2dace_A_z_w_concorr_me_d_2_s_362 = size(z_w_concorr_me, dim=3), &
      f2dace_OA_z_kin_hor_e_d_0_s_363 = lbound(z_kin_hor_e, dim=1), &
      f2dace_OA_z_kin_hor_e_d_1_s_364 = lbound(z_kin_hor_e, dim=2), &
      f2dace_OA_z_kin_hor_e_d_2_s_365 = lbound(z_kin_hor_e, dim=3) &
    )

    call copy_back_global_data_type(copy_or_ptr_global_data)
    call copy_back_t_nh_diag(p_diag, copy_or_ptr_p_diag)
    call copy_back_t_int_state(p_int, copy_or_ptr_p_int)
    call copy_back_t_nh_metrics(p_metrics, copy_or_ptr_p_metrics)
    call copy_back_t_patch(p_patch, copy_or_ptr_p_patch)
    call copy_back_t_nh_prog(p_prog, copy_or_ptr_p_prog)


  end subroutine run_velocity_no_nproma_if_prop_lvn_only_1_istep_2

  subroutine run_velocity_no_nproma_if_prop_lvn_only_1_istep_2_verification( &
    p_diag, &
    p_int, &
    p_metrics, &
    p_patch, &
    p_prog, &
    z_kin_hor_e, &
    z_vt_ie, &
    z_w_concorr_me, &
    dt_linintp_ubc, &
    dtime, &
    istep, &
    ldeepatmo, &
    lvn_only, &
    ntnd &
  )
    type(t_nh_diag), target :: p_diag
    type(t_int_state), target :: p_int
    type(t_nh_metrics), target :: p_metrics
    type(t_patch), target :: p_patch
    type(t_nh_prog), target :: p_prog
    real(kind=c_double), dimension(:,:,:), target :: z_kin_hor_e
    real(kind=c_double), dimension(:,:,:), target :: z_vt_ie
    real(kind=c_double), dimension(:,:,:), target :: z_w_concorr_me
    real(kind=c_double) :: dt_linintp_ubc
    real(kind=c_double) :: dtime
    integer(kind=c_int) :: istep
    integer(kind=c_int) :: ldeepatmo
    integer(kind=c_int) :: lvn_only
    integer(kind=c_int) :: ntnd



    !$ACC WAIT



    call check_initializations()

    copy_or_ptr_global_data = copy_in_global_data_type( &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_diag = copy_in_t_nh_diag( &
    fortran_obj=p_diag, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_int = copy_in_t_int_state( &
    fortran_obj=p_int, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_metrics = copy_in_t_nh_metrics( &
    fortran_obj=p_metrics, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_patch = copy_in_t_patch( &
    fortran_obj=p_patch, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_prog = copy_in_t_nh_prog( &
    fortran_obj=p_prog, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
#ifndef _OPENACC
    copy_or_ptr_z_kin_hor_e = copy_in_float64_3d_array( &
    fortran_array=z_kin_hor_e, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_kin_hor_e = copy_in_float64_3d_array( &
    fortran_array=z_kin_hor_e, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_vt_ie = copy_in_float64_3d_array( &
    fortran_array=z_vt_ie, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_vt_ie = copy_in_float64_3d_array( &
    fortran_array=z_vt_ie, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_w_concorr_me = copy_in_float64_3d_array( &
    fortran_array=z_w_concorr_me, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_w_concorr_me = copy_in_float64_3d_array( &
    fortran_array=z_w_concorr_me, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif


    if (is_initialized .eqv. .false.) then
      is_initialized = .true.
      dace_state = dace_init_velocity_no_nproma_if_prop_lvn_only_1_istep_2( &
      global_data = copy_or_ptr_global_data, &
      p_diag = copy_or_ptr_p_diag, &
      p_int = copy_or_ptr_p_int, &
      p_metrics = copy_or_ptr_p_metrics, &
      p_patch = copy_or_ptr_p_patch, &
      p_prog = copy_or_ptr_p_prog, &
      z_kin_hor_e = copy_or_ptr_z_kin_hor_e, &
      z_vt_ie = copy_or_ptr_z_vt_ie, &
      z_w_concorr_me = copy_or_ptr_z_w_concorr_me, &
      dt_linintp_ubc = dt_linintp_ubc, &
      dtime = dtime, &
      istep = istep, &
      ldeepatmo = ldeepatmo, &
      lvn_only = lvn_only, &
      ntnd = ntnd, &
      f2dace_A_z_kin_hor_e_d_0_s_363 = size(z_kin_hor_e, dim=1), &
      f2dace_A_z_kin_hor_e_d_1_s_364 = size(z_kin_hor_e, dim=2), &
      f2dace_A_z_kin_hor_e_d_2_s_365 = size(z_kin_hor_e, dim=3), &
      f2dace_A_z_vt_ie_d_0_s_366 = size(z_vt_ie, dim=1), &
      f2dace_A_z_vt_ie_d_1_s_367 = size(z_vt_ie, dim=2), &
      f2dace_A_z_vt_ie_d_2_s_368 = size(z_vt_ie, dim=3), &
      f2dace_A_z_w_concorr_me_d_0_s_360 = size(z_w_concorr_me, dim=1), &
      f2dace_A_z_w_concorr_me_d_1_s_361 = size(z_w_concorr_me, dim=2), &
      f2dace_A_z_w_concorr_me_d_2_s_362 = size(z_w_concorr_me, dim=3), &
      f2dace_OA_z_kin_hor_e_d_0_s_363 = lbound(z_kin_hor_e, dim=1), &
      f2dace_OA_z_kin_hor_e_d_1_s_364 = lbound(z_kin_hor_e, dim=2), &
      f2dace_OA_z_kin_hor_e_d_2_s_365 = lbound(z_kin_hor_e, dim=3) &
    )
    end if

    call dace_program_velocity_no_nproma_if_prop_lvn_only_1_istep_2( &
      state = dace_state, &
      global_data = copy_or_ptr_global_data, &
      p_diag = copy_or_ptr_p_diag, &
      p_int = copy_or_ptr_p_int, &
      p_metrics = copy_or_ptr_p_metrics, &
      p_patch = copy_or_ptr_p_patch, &
      p_prog = copy_or_ptr_p_prog, &
      z_kin_hor_e = copy_or_ptr_z_kin_hor_e, &
      z_vt_ie = copy_or_ptr_z_vt_ie, &
      z_w_concorr_me = copy_or_ptr_z_w_concorr_me, &
      dt_linintp_ubc = dt_linintp_ubc, &
      dtime = dtime, &
      istep = istep, &
      ldeepatmo = ldeepatmo, &
      lvn_only = lvn_only, &
      ntnd = ntnd, &
      f2dace_A_z_kin_hor_e_d_0_s_363 = size(z_kin_hor_e, dim=1), &
      f2dace_A_z_kin_hor_e_d_1_s_364 = size(z_kin_hor_e, dim=2), &
      f2dace_A_z_kin_hor_e_d_2_s_365 = size(z_kin_hor_e, dim=3), &
      f2dace_A_z_vt_ie_d_0_s_366 = size(z_vt_ie, dim=1), &
      f2dace_A_z_vt_ie_d_1_s_367 = size(z_vt_ie, dim=2), &
      f2dace_A_z_vt_ie_d_2_s_368 = size(z_vt_ie, dim=3), &
      f2dace_A_z_w_concorr_me_d_0_s_360 = size(z_w_concorr_me, dim=1), &
      f2dace_A_z_w_concorr_me_d_1_s_361 = size(z_w_concorr_me, dim=2), &
      f2dace_A_z_w_concorr_me_d_2_s_362 = size(z_w_concorr_me, dim=3), &
      f2dace_OA_z_kin_hor_e_d_0_s_363 = lbound(z_kin_hor_e, dim=1), &
      f2dace_OA_z_kin_hor_e_d_1_s_364 = lbound(z_kin_hor_e, dim=2), &
      f2dace_OA_z_kin_hor_e_d_2_s_365 = lbound(z_kin_hor_e, dim=3) &
    )
  end subroutine run_velocity_no_nproma_if_prop_lvn_only_1_istep_2_verification

  subroutine verify_velocity_no_nproma_if_prop_lvn_only_1_istep_2( &
    p_diag, &
    p_int, &
    p_metrics, &
    p_patch, &
    p_prog, &
    z_kin_hor_e, &
    z_vt_ie, &
    z_w_concorr_me, &
    dt_linintp_ubc, &
    dtime, &
    istep, &
    ldeepatmo, &
    lvn_only, &
    ntnd &
  )
    type(t_nh_diag), target :: p_diag
    type(t_int_state), target :: p_int
    type(t_nh_metrics), target :: p_metrics
    type(t_patch), target :: p_patch
    type(t_nh_prog), target :: p_prog
    real(kind=c_double), dimension(:,:,:), target :: z_kin_hor_e
    real(kind=c_double), dimension(:,:,:), target :: z_vt_ie
    real(kind=c_double), dimension(:,:,:), target :: z_w_concorr_me
    real(kind=c_double) :: dt_linintp_ubc
    real(kind=c_double) :: dtime
    integer(kind=c_int) :: istep
    integer(kind=c_int) :: ldeepatmo
    integer(kind=c_int) :: lvn_only
    integer(kind=c_int) :: ntnd


    logical :: local_result, result

    !$ACC WAIT

    result = .true.
    local_result = .true.


    if (is_initialized .eqv. .false.) then
      print *, "verify_velocity_no_nproma_if_prop_lvn_only_1_istep_2: dace state is not initialized"
    end if

    call check_initializations()


    call compare_global_data_type_struct( &
        actual=copy_or_ptr_global_data, &
        result=local_result, &
        struct_expr="global_data" &
    )

    result = result .and. local_result

    call compare_t_nh_diag_struct( &
        actual=copy_or_ptr_p_diag, &
        ref=p_diag, &
        result=local_result, &
        struct_expr="p_diag" &
    )

    result = result .and. local_result

    call compare_t_int_state_struct( &
        actual=copy_or_ptr_p_int, &
        ref=p_int, &
        result=local_result, &
        struct_expr="p_int" &
    )

    result = result .and. local_result

    call compare_t_nh_metrics_struct( &
        actual=copy_or_ptr_p_metrics, &
        ref=p_metrics, &
        result=local_result, &
        struct_expr="p_metrics" &
    )

    result = result .and. local_result

    call compare_t_patch_struct( &
        actual=copy_or_ptr_p_patch, &
        ref=p_patch, &
        result=local_result, &
        struct_expr="p_patch" &
    )

    result = result .and. local_result

    call compare_t_nh_prog_struct( &
        actual=copy_or_ptr_p_prog, &
        ref=p_prog, &
        result=local_result, &
        struct_expr="p_prog" &
    )

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_kin_hor_e, &
        ref=z_kin_hor_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_kin_hor_e" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_kin_hor_e, &
        ref=z_kin_hor_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_kin_hor_e" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_vt_ie, &
        ref=z_vt_ie, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_vt_ie" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_vt_ie, &
        ref=z_vt_ie, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_vt_ie" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_w_concorr_me, &
        ref=z_w_concorr_me, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_w_concorr_me" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_w_concorr_me, &
        ref=z_w_concorr_me, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_w_concorr_me" &
    )

#endif

    result = result .and. local_result

    if (.not. result) then
      print *, "verify_velocity_no_nproma_if_prop_lvn_only_1_istep_2: Failed verification"
    else
      print *, "verify_velocity_no_nproma_if_prop_lvn_only_1_istep_2: Verification successful :)"
    end if

  end subroutine verify_velocity_no_nproma_if_prop_lvn_only_1_istep_2

end module mo_velocity_no_nproma_if_prop_lvn_only_1_istep_2_bindings
