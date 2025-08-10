! Auto-generated file by "/data/ben/spcl/icon-dace/upstream-repo/./icon-model/sdfgs/utils/dace-genfi.py"
module mo_solve_nh_predictor_pre_bindings

  use iso_c_binding

  ! use external struct definitions
  use mo_solve_nh_dace_structs

  use mo_intp_data_strc, only: &
    t_int_state
  use mo_nonhydro_types, only: &
    t_nh_state, &
    t_nh_diag, &
    t_nh_ref, &
    t_nh_metrics, &
    t_nh_prog
  use mo_model_domain, only: &
    t_patch, &
    t_grid_cells, &
    t_grid_edges, &
    t_tangent_vectors, &
    t_grid_vertices
  use mo_decomposition_tools, only: &
    t_grid_domain_decomp_info
  use mo_prepadv_types, only: &
    t_prepare_adv
  use mo_dynamics_config, only: &
    ldeepatmo
  use mo_grid_config, only: &
    l_limited_area
  use mo_gridref_config, only: &
    grf_intmethod_e
  use mo_init_vgrid, only: &
    nflatlev
  use mo_initicon_config, only: &
    is_iau_active, &
    iau_wgt_dyn
  use mo_mpi, only: &
    i_am_accel_node
  use mo_nonhydrostatic_config, only: &
    itime_scheme, &
    lextra_diffu, &
    rayleigh_type, &
    iadv_rhotheta, &
    igradp_method, &
    kstart_dd3d
  use mo_parallel_config, only: &
    nproma
  use mo_run_config, only: &
    lvert_nest, &
    timers_level
  use mo_timer, only: &
    timer_solve_nh_veltend, &
    timer_solve_nh_cellcomp, &
    timer_solve_nh_vnupd, &
    timer_intp
  use mo_vertical_grid, only: &
    nrdmax, &
    nflat_gradp

  implicit none

  private
  public :: run_solve_nh_predictor_pre
  public :: run_solve_nh_predictor_pre_verification
  public :: verify_solve_nh_predictor_pre
  public :: dace_init_solve_nh_predictor_pre
  public :: dace_exit_solve_nh_predictor_pre
  public :: dace_program_solve_nh_predictor_pre


  type(c_ptr) :: properly_cached_t_patch = c_null_ptr


  logical :: is_initialized = .false.
  type(c_ptr) :: dace_state = C_NULL_PTR

  type(c_ptr) :: cached_shallow_copy_global_data = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_p_int = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_p_nh = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_p_nh_prog_nnew = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_p_nh_prog_nnow = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_p_patch = C_NULL_PTR
  type(c_ptr) :: cached_shallow_copy_prep_adv = C_NULL_PTR

  type(c_ptr) :: copy_or_ptr_bdy_divdamp = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_enh_divdamp_fac = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_global_data = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_p_int = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_p_nh = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_p_nh_prog_nnew = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_p_nh_prog_nnow = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_p_patch = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_prep_adv = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_scal_divdamp = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_alpha = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_beta = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_contr_w_fl_l = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_dexner_dz_c = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_dwdz_dd = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_exner_ex_pr = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_exner_expl = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_exner_ic = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_flxdiv_mass = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_flxdiv_theta = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_grad_rth = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_graddiv2_vn = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_graddiv_vn = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_gradh_exner = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_hydro_corr = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_kin_hor_e = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_mflx_top = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_q = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_raylfac = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_rho_e = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_rho_expl = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_rho_v = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_rth_pr = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_th_ddz_exner_c = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_theta_v_e = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_theta_v_fl_e = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_theta_v_pr_ic = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_theta_v_v = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_vn_avg = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_vt_ie = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_w_concorr_mc = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_w_concorr_me = C_NULL_PTR
  type(c_ptr) :: copy_or_ptr_z_w_expl = C_NULL_PTR

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


  type(c_ptr) function dace_init_solve_nh_predictor_pre( &
    bdy_divdamp, &
    enh_divdamp_fac, &
    global_data, &
    p_int, &
    p_nh, &
    p_nh_prog_nnew, &
    p_nh_prog_nnow, &
    p_patch, &
    prep_adv, &
    scal_divdamp, &
    z_alpha, &
    z_beta, &
    z_contr_w_fl_l, &
    z_dexner_dz_c, &
    z_dwdz_dd, &
    z_exner_ex_pr, &
    z_exner_expl, &
    z_exner_ic, &
    z_flxdiv_mass, &
    z_flxdiv_theta, &
    z_grad_rth, &
    z_graddiv2_vn, &
    z_graddiv_vn, &
    z_gradh_exner, &
    z_hydro_corr, &
    z_kin_hor_e, &
    z_mflx_top, &
    z_q, &
    z_raylfac, &
    z_rho_e, &
    z_rho_expl, &
    z_rho_v, &
    z_rth_pr, &
    z_th_ddz_exner_c, &
    z_theta_v_e, &
    z_theta_v_fl_e, &
    z_theta_v_pr_ic, &
    z_theta_v_v, &
    z_vn_avg, &
    z_vt_ie, &
    z_w_concorr_mc, &
    z_w_concorr_me, &
    z_w_expl, &
    f2dace_OPTIONAL_lacc, &
    alin, &
    aqdr, &
    bqdr, &
    df32, &
    df42, &
    distv_bary_1, &
    distv_bary_2, &
    dt_linintp_ubc, &
    dt_linintp_ubc_nnew, &
    dt_linintp_ubc_nnow, &
    dt_shift, &
    dthalf, &
    dtime, &
    dz32, &
    dz42, &
    dzlin, &
    dzqdr, &
    i_endblk, &
    i_endidx, &
    i_startblk, &
    i_startidx, &
    idyn_timestep, &
    ishift, &
    istep, &
    jb, &
    jc, &
    je, &
    jg, &
    jk, &
    jk_start, &
    jks, &
    jstep, &
    l_child_vertnest, &
    l_init, &
    l_recompute, &
    l_vert_nested, &
    lacc, &
    lclean_mflx, &
    lprep_adv, &
    lsave_mflx, &
    lvn_only, &
    lvn_pos, &
    nblks_gradp, &
    nlen_gradp, &
    nlev, &
    nlevp1, &
    nnew, &
    nnow, &
    nproma_gradp, &
    npromz_gradp, &
    nshift, &
    nshift_total, &
    ntl1, &
    ntl2, &
    nvar, &
    r_dtimensubsteps, &
    r_nsubsteps, &
    rl_end, &
    rl_start, &
    scal_divdamp_o2, &
    wgt_nnew_rth, &
    wgt_nnew_vel, &
    wgt_nnow_rth, &
    wgt_nnow_vel, &
    z_a, &
    z_b, &
    z_c, &
    z_d_vn_dmp, &
    z_d_vn_iau, &
    z_ddt_vn_apc, &
    z_ddt_vn_cor, &
    z_ddt_vn_dyn, &
    z_ddt_vn_pgr, &
    z_ddt_vn_ray, &
    z_g, &
    z_gamma, &
    z_ntdistv_bary_1, &
    z_ntdistv_bary_2, &
    z_rho_tavg, &
    z_rho_tavg_m1, &
    z_theta1, &
    z_theta2, &
    z_theta_tavg, &
    z_theta_tavg_m1, &
    z_theta_v_pr_mc, &
    z_theta_v_pr_mc_m1, &
    z_w_backtraj, &
    zf &
  ) &
    bind(c, name="__dace_init_solve_nh_predictor_pre")
    use iso_c_binding

    type(c_ptr), value :: bdy_divdamp
    type(c_ptr), value :: enh_divdamp_fac
    type(c_ptr), value :: global_data
    type(c_ptr), value :: p_int
    type(c_ptr), value :: p_nh
    type(c_ptr), value :: p_nh_prog_nnew
    type(c_ptr), value :: p_nh_prog_nnow
    type(c_ptr), value :: p_patch
    type(c_ptr), value :: prep_adv
    type(c_ptr), value :: scal_divdamp
    type(c_ptr), value :: z_alpha
    type(c_ptr), value :: z_beta
    type(c_ptr), value :: z_contr_w_fl_l
    type(c_ptr), value :: z_dexner_dz_c
    type(c_ptr), value :: z_dwdz_dd
    type(c_ptr), value :: z_exner_ex_pr
    type(c_ptr), value :: z_exner_expl
    type(c_ptr), value :: z_exner_ic
    type(c_ptr), value :: z_flxdiv_mass
    type(c_ptr), value :: z_flxdiv_theta
    type(c_ptr), value :: z_grad_rth
    type(c_ptr), value :: z_graddiv2_vn
    type(c_ptr), value :: z_graddiv_vn
    type(c_ptr), value :: z_gradh_exner
    type(c_ptr), value :: z_hydro_corr
    type(c_ptr), value :: z_kin_hor_e
    type(c_ptr), value :: z_mflx_top
    type(c_ptr), value :: z_q
    type(c_ptr), value :: z_raylfac
    type(c_ptr), value :: z_rho_e
    type(c_ptr), value :: z_rho_expl
    type(c_ptr), value :: z_rho_v
    type(c_ptr), value :: z_rth_pr
    type(c_ptr), value :: z_th_ddz_exner_c
    type(c_ptr), value :: z_theta_v_e
    type(c_ptr), value :: z_theta_v_fl_e
    type(c_ptr), value :: z_theta_v_pr_ic
    type(c_ptr), value :: z_theta_v_v
    type(c_ptr), value :: z_vn_avg
    type(c_ptr), value :: z_vt_ie
    type(c_ptr), value :: z_w_concorr_mc
    type(c_ptr), value :: z_w_concorr_me
    type(c_ptr), value :: z_w_expl
    integer(kind=c_int), value :: f2dace_OPTIONAL_lacc
    real(kind=c_double), value :: alin
    real(kind=c_double), value :: aqdr
    real(kind=c_double), value :: bqdr
    real(kind=c_double), value :: df32
    real(kind=c_double), value :: df42
    real(kind=c_double), value :: distv_bary_1
    real(kind=c_double), value :: distv_bary_2
    real(kind=c_double), value :: dt_linintp_ubc
    real(kind=c_double), value :: dt_linintp_ubc_nnew
    real(kind=c_double), value :: dt_linintp_ubc_nnow
    real(kind=c_double), value :: dt_shift
    real(kind=c_double), value :: dthalf
    real(kind=c_double), value :: dtime
    real(kind=c_double), value :: dz32
    real(kind=c_double), value :: dz42
    real(kind=c_double), value :: dzlin
    real(kind=c_double), value :: dzqdr
    integer(kind=c_int), value :: i_endblk
    integer(kind=c_int), value :: i_endidx
    integer(kind=c_int), value :: i_startblk
    integer(kind=c_int), value :: i_startidx
    integer(kind=c_int), value :: idyn_timestep
    integer(kind=c_int), value :: ishift
    integer(kind=c_int), value :: istep
    integer(kind=c_int), value :: jb
    integer(kind=c_int), value :: jc
    integer(kind=c_int), value :: je
    integer(kind=c_int), value :: jg
    integer(kind=c_int), value :: jk
    integer(kind=c_int), value :: jk_start
    integer(kind=c_int), value :: jks
    integer(kind=c_int), value :: jstep
    integer(kind=c_int), value :: l_child_vertnest
    integer(kind=c_int), value :: l_init
    integer(kind=c_int), value :: l_recompute
    integer(kind=c_int), value :: l_vert_nested
    integer(kind=c_int), value :: lacc
    integer(kind=c_int), value :: lclean_mflx
    integer(kind=c_int), value :: lprep_adv
    integer(kind=c_int), value :: lsave_mflx
    integer(kind=c_int), value :: lvn_only
    integer(kind=c_int), value :: lvn_pos
    integer(kind=c_int), value :: nblks_gradp
    integer(kind=c_int), value :: nlen_gradp
    integer(kind=c_int), value :: nlev
    integer(kind=c_int), value :: nlevp1
    integer(kind=c_int), value :: nnew
    integer(kind=c_int), value :: nnow
    integer(kind=c_int), value :: nproma_gradp
    integer(kind=c_int), value :: npromz_gradp
    integer(kind=c_int), value :: nshift
    integer(kind=c_int), value :: nshift_total
    integer(kind=c_int), value :: ntl1
    integer(kind=c_int), value :: ntl2
    integer(kind=c_int), value :: nvar
    real(kind=c_double), value :: r_dtimensubsteps
    real(kind=c_double), value :: r_nsubsteps
    integer(kind=c_int), value :: rl_end
    integer(kind=c_int), value :: rl_start
    real(kind=c_double), value :: scal_divdamp_o2
    real(kind=c_double), value :: wgt_nnew_rth
    real(kind=c_double), value :: wgt_nnew_vel
    real(kind=c_double), value :: wgt_nnow_rth
    real(kind=c_double), value :: wgt_nnow_vel
    real(kind=c_double), value :: z_a
    real(kind=c_double), value :: z_b
    real(kind=c_double), value :: z_c
    real(kind=c_double), value :: z_d_vn_dmp
    real(kind=c_double), value :: z_d_vn_iau
    real(kind=c_double), value :: z_ddt_vn_apc
    real(kind=c_double), value :: z_ddt_vn_cor
    real(kind=c_double), value :: z_ddt_vn_dyn
    real(kind=c_double), value :: z_ddt_vn_pgr
    real(kind=c_double), value :: z_ddt_vn_ray
    real(kind=c_double), value :: z_g
    real(kind=c_double), value :: z_gamma
    real(kind=c_double), value :: z_ntdistv_bary_1
    real(kind=c_double), value :: z_ntdistv_bary_2
    real(kind=c_double), value :: z_rho_tavg
    real(kind=c_double), value :: z_rho_tavg_m1
    real(kind=c_double), value :: z_theta1
    real(kind=c_double), value :: z_theta2
    real(kind=c_double), value :: z_theta_tavg
    real(kind=c_double), value :: z_theta_tavg_m1
    real(kind=c_double), value :: z_theta_v_pr_mc
    real(kind=c_double), value :: z_theta_v_pr_mc_m1
    real(kind=c_double), value :: z_w_backtraj
    real(kind=c_double), value :: zf
  end function dace_init_solve_nh_predictor_pre

  integer(c_int) function dace_exit_solve_nh_predictor_pre(state) &
    bind(c, name="__dace_exit_solve_nh_predictor_pre")
    use iso_c_binding

    type(c_ptr), value :: state
  end function dace_exit_solve_nh_predictor_pre

  subroutine dace_program_solve_nh_predictor_pre( &
    state, &
    bdy_divdamp, &
    enh_divdamp_fac, &
    global_data, &
    p_int, &
    p_nh, &
    p_nh_prog_nnew, &
    p_nh_prog_nnow, &
    p_patch, &
    prep_adv, &
    scal_divdamp, &
    z_alpha, &
    z_beta, &
    z_contr_w_fl_l, &
    z_dexner_dz_c, &
    z_dwdz_dd, &
    z_exner_ex_pr, &
    z_exner_expl, &
    z_exner_ic, &
    z_flxdiv_mass, &
    z_flxdiv_theta, &
    z_grad_rth, &
    z_graddiv2_vn, &
    z_graddiv_vn, &
    z_gradh_exner, &
    z_hydro_corr, &
    z_kin_hor_e, &
    z_mflx_top, &
    z_q, &
    z_raylfac, &
    z_rho_e, &
    z_rho_expl, &
    z_rho_v, &
    z_rth_pr, &
    z_th_ddz_exner_c, &
    z_theta_v_e, &
    z_theta_v_fl_e, &
    z_theta_v_pr_ic, &
    z_theta_v_v, &
    z_vn_avg, &
    z_vt_ie, &
    z_w_concorr_mc, &
    z_w_concorr_me, &
    z_w_expl, &
    f2dace_OPTIONAL_lacc, &
    alin, &
    aqdr, &
    bqdr, &
    df32, &
    df42, &
    distv_bary_1, &
    distv_bary_2, &
    dt_linintp_ubc, &
    dt_linintp_ubc_nnew, &
    dt_linintp_ubc_nnow, &
    dt_shift, &
    dthalf, &
    dtime, &
    dz32, &
    dz42, &
    dzlin, &
    dzqdr, &
    i_endblk, &
    i_endidx, &
    i_startblk, &
    i_startidx, &
    idyn_timestep, &
    ishift, &
    istep, &
    jb, &
    jc, &
    je, &
    jg, &
    jk, &
    jk_start, &
    jks, &
    jstep, &
    l_child_vertnest, &
    l_init, &
    l_recompute, &
    l_vert_nested, &
    lacc, &
    lclean_mflx, &
    lprep_adv, &
    lsave_mflx, &
    lvn_only, &
    lvn_pos, &
    nblks_gradp, &
    nlen_gradp, &
    nlev, &
    nlevp1, &
    nnew, &
    nnow, &
    nproma_gradp, &
    npromz_gradp, &
    nshift, &
    nshift_total, &
    ntl1, &
    ntl2, &
    nvar, &
    r_dtimensubsteps, &
    r_nsubsteps, &
    rl_end, &
    rl_start, &
    scal_divdamp_o2, &
    wgt_nnew_rth, &
    wgt_nnew_vel, &
    wgt_nnow_rth, &
    wgt_nnow_vel, &
    z_a, &
    z_b, &
    z_c, &
    z_d_vn_dmp, &
    z_d_vn_iau, &
    z_ddt_vn_apc, &
    z_ddt_vn_cor, &
    z_ddt_vn_dyn, &
    z_ddt_vn_pgr, &
    z_ddt_vn_ray, &
    z_g, &
    z_gamma, &
    z_ntdistv_bary_1, &
    z_ntdistv_bary_2, &
    z_rho_tavg, &
    z_rho_tavg_m1, &
    z_theta1, &
    z_theta2, &
    z_theta_tavg, &
    z_theta_tavg_m1, &
    z_theta_v_pr_mc, &
    z_theta_v_pr_mc_m1, &
    z_w_backtraj, &
    zf &
  ) &
    bind(c, name="__program_solve_nh_predictor_pre")
    use iso_c_binding

    type(c_ptr), value :: state
    type(c_ptr), value :: bdy_divdamp
    type(c_ptr), value :: enh_divdamp_fac
    type(c_ptr), value :: global_data
    type(c_ptr), value :: p_int
    type(c_ptr), value :: p_nh
    type(c_ptr), value :: p_nh_prog_nnew
    type(c_ptr), value :: p_nh_prog_nnow
    type(c_ptr), value :: p_patch
    type(c_ptr), value :: prep_adv
    type(c_ptr), value :: scal_divdamp
    type(c_ptr), value :: z_alpha
    type(c_ptr), value :: z_beta
    type(c_ptr), value :: z_contr_w_fl_l
    type(c_ptr), value :: z_dexner_dz_c
    type(c_ptr), value :: z_dwdz_dd
    type(c_ptr), value :: z_exner_ex_pr
    type(c_ptr), value :: z_exner_expl
    type(c_ptr), value :: z_exner_ic
    type(c_ptr), value :: z_flxdiv_mass
    type(c_ptr), value :: z_flxdiv_theta
    type(c_ptr), value :: z_grad_rth
    type(c_ptr), value :: z_graddiv2_vn
    type(c_ptr), value :: z_graddiv_vn
    type(c_ptr), value :: z_gradh_exner
    type(c_ptr), value :: z_hydro_corr
    type(c_ptr), value :: z_kin_hor_e
    type(c_ptr), value :: z_mflx_top
    type(c_ptr), value :: z_q
    type(c_ptr), value :: z_raylfac
    type(c_ptr), value :: z_rho_e
    type(c_ptr), value :: z_rho_expl
    type(c_ptr), value :: z_rho_v
    type(c_ptr), value :: z_rth_pr
    type(c_ptr), value :: z_th_ddz_exner_c
    type(c_ptr), value :: z_theta_v_e
    type(c_ptr), value :: z_theta_v_fl_e
    type(c_ptr), value :: z_theta_v_pr_ic
    type(c_ptr), value :: z_theta_v_v
    type(c_ptr), value :: z_vn_avg
    type(c_ptr), value :: z_vt_ie
    type(c_ptr), value :: z_w_concorr_mc
    type(c_ptr), value :: z_w_concorr_me
    type(c_ptr), value :: z_w_expl
    integer(kind=c_int), value :: f2dace_OPTIONAL_lacc
    real(kind=c_double), value :: alin
    real(kind=c_double), value :: aqdr
    real(kind=c_double), value :: bqdr
    real(kind=c_double), value :: df32
    real(kind=c_double), value :: df42
    real(kind=c_double), value :: distv_bary_1
    real(kind=c_double), value :: distv_bary_2
    real(kind=c_double), value :: dt_linintp_ubc
    real(kind=c_double), value :: dt_linintp_ubc_nnew
    real(kind=c_double), value :: dt_linintp_ubc_nnow
    real(kind=c_double), value :: dt_shift
    real(kind=c_double), value :: dthalf
    real(kind=c_double), value :: dtime
    real(kind=c_double), value :: dz32
    real(kind=c_double), value :: dz42
    real(kind=c_double), value :: dzlin
    real(kind=c_double), value :: dzqdr
    integer(kind=c_int), value :: i_endblk
    integer(kind=c_int), value :: i_endidx
    integer(kind=c_int), value :: i_startblk
    integer(kind=c_int), value :: i_startidx
    integer(kind=c_int), value :: idyn_timestep
    integer(kind=c_int), value :: ishift
    integer(kind=c_int), value :: istep
    integer(kind=c_int), value :: jb
    integer(kind=c_int), value :: jc
    integer(kind=c_int), value :: je
    integer(kind=c_int), value :: jg
    integer(kind=c_int), value :: jk
    integer(kind=c_int), value :: jk_start
    integer(kind=c_int), value :: jks
    integer(kind=c_int), value :: jstep
    integer(kind=c_int), value :: l_child_vertnest
    integer(kind=c_int), value :: l_init
    integer(kind=c_int), value :: l_recompute
    integer(kind=c_int), value :: l_vert_nested
    integer(kind=c_int), value :: lacc
    integer(kind=c_int), value :: lclean_mflx
    integer(kind=c_int), value :: lprep_adv
    integer(kind=c_int), value :: lsave_mflx
    integer(kind=c_int), value :: lvn_only
    integer(kind=c_int), value :: lvn_pos
    integer(kind=c_int), value :: nblks_gradp
    integer(kind=c_int), value :: nlen_gradp
    integer(kind=c_int), value :: nlev
    integer(kind=c_int), value :: nlevp1
    integer(kind=c_int), value :: nnew
    integer(kind=c_int), value :: nnow
    integer(kind=c_int), value :: nproma_gradp
    integer(kind=c_int), value :: npromz_gradp
    integer(kind=c_int), value :: nshift
    integer(kind=c_int), value :: nshift_total
    integer(kind=c_int), value :: ntl1
    integer(kind=c_int), value :: ntl2
    integer(kind=c_int), value :: nvar
    real(kind=c_double), value :: r_dtimensubsteps
    real(kind=c_double), value :: r_nsubsteps
    integer(kind=c_int), value :: rl_end
    integer(kind=c_int), value :: rl_start
    real(kind=c_double), value :: scal_divdamp_o2
    real(kind=c_double), value :: wgt_nnew_rth
    real(kind=c_double), value :: wgt_nnew_vel
    real(kind=c_double), value :: wgt_nnow_rth
    real(kind=c_double), value :: wgt_nnow_vel
    real(kind=c_double), value :: z_a
    real(kind=c_double), value :: z_b
    real(kind=c_double), value :: z_c
    real(kind=c_double), value :: z_d_vn_dmp
    real(kind=c_double), value :: z_d_vn_iau
    real(kind=c_double), value :: z_ddt_vn_apc
    real(kind=c_double), value :: z_ddt_vn_cor
    real(kind=c_double), value :: z_ddt_vn_dyn
    real(kind=c_double), value :: z_ddt_vn_pgr
    real(kind=c_double), value :: z_ddt_vn_ray
    real(kind=c_double), value :: z_g
    real(kind=c_double), value :: z_gamma
    real(kind=c_double), value :: z_ntdistv_bary_1
    real(kind=c_double), value :: z_ntdistv_bary_2
    real(kind=c_double), value :: z_rho_tavg
    real(kind=c_double), value :: z_rho_tavg_m1
    real(kind=c_double), value :: z_theta1
    real(kind=c_double), value :: z_theta2
    real(kind=c_double), value :: z_theta_tavg
    real(kind=c_double), value :: z_theta_tavg_m1
    real(kind=c_double), value :: z_theta_v_pr_mc
    real(kind=c_double), value :: z_theta_v_pr_mc_m1
    real(kind=c_double), value :: z_w_backtraj
    real(kind=c_double), value :: zf
  end subroutine dace_program_solve_nh_predictor_pre

end interface

interface logical_fix_1d
  module procedure logical_to_int_1d
  module procedure int_to_int_1d
end interface logical_fix_1d

interface logical_fix_4d
  module procedure logical_to_int_4d
  module procedure int_to_int_4d
end interface logical_fix_4d

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

    dace_rich_obj%ldeepatmo = ldeepatmo
    dace_rich_obj%l_limited_area = l_limited_area
    dace_rich_obj%grf_intmethod_e = grf_intmethod_e
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
    dace_rich_obj%is_iau_active = is_iau_active
    dace_rich_obj%iau_wgt_dyn = iau_wgt_dyn
    dace_rich_obj%i_am_accel_node = i_am_accel_node
    dace_rich_obj%itime_scheme = itime_scheme
    dace_rich_obj%lextra_diffu = lextra_diffu
    dace_rich_obj%rayleigh_type = rayleigh_type
    dace_rich_obj%iadv_rhotheta = iadv_rhotheta
    dace_rich_obj%igradp_method = igradp_method
    dace_rich_obj%kstart_dd3d = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(kstart_dd3d), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
#if defined(DACE_SUBST_VERIFY)
    if (10 /= size(kstart_dd3d, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 'global_data.kstart_dd3d'"//char(10), &
        "    - actual = (", &
        size(kstart_dd3d, dim=1), &
        "), config propagated = (10)"
    end if
#endif
    dace_rich_obj%nproma = nproma
    dace_rich_obj%lvert_nest = lvert_nest
    dace_rich_obj%timers_level = timers_level
    dace_rich_obj%timer_solve_nh_veltend = timer_solve_nh_veltend
    dace_rich_obj%timer_solve_nh_cellcomp = timer_solve_nh_cellcomp
    dace_rich_obj%timer_solve_nh_vnupd = timer_solve_nh_vnupd
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
    dace_rich_obj%nflat_gradp = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(nflat_gradp), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
#if defined(DACE_SUBST_VERIFY)
    if (10 /= size(nflat_gradp, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 'global_data.nflat_gradp'"//char(10), &
        "    - actual = (", &
        size(nflat_gradp, dim=1), &
        "), config propagated = (10)"
    end if
#endif

  end function copy_in_global_data_type

  function copy_in_t_int_state(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_int_state), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_int_state), pointer :: dace_rich_obj
    type(dace_t_int_state) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_c_lin_e_d_0_s = size(fortran_obj%c_lin_e, dim=1)
    dace_rich_obj%f2dace_SOA_c_lin_e_d_0_s = lbound(fortran_obj%c_lin_e, dim=1)
    dace_rich_obj%f2dace_SA_c_lin_e_d_1_s = size(fortran_obj%c_lin_e, dim=2)
    dace_rich_obj%f2dace_SOA_c_lin_e_d_1_s = lbound(fortran_obj%c_lin_e, dim=2)
    dace_rich_obj%f2dace_SA_c_lin_e_d_2_s = size(fortran_obj%c_lin_e, dim=3)
    dace_rich_obj%f2dace_SOA_c_lin_e_d_2_s = lbound(fortran_obj%c_lin_e, dim=3)
    dace_rich_obj%f2dace_SA_e_bln_c_s_d_0_s = size(fortran_obj%e_bln_c_s, dim=1)
    dace_rich_obj%f2dace_SOA_e_bln_c_s_d_0_s = lbound(fortran_obj%e_bln_c_s, dim=1)
    dace_rich_obj%f2dace_SA_e_bln_c_s_d_1_s = size(fortran_obj%e_bln_c_s, dim=2)
    dace_rich_obj%f2dace_SOA_e_bln_c_s_d_1_s = lbound(fortran_obj%e_bln_c_s, dim=2)
    dace_rich_obj%f2dace_SA_e_bln_c_s_d_2_s = size(fortran_obj%e_bln_c_s, dim=3)
    dace_rich_obj%f2dace_SOA_e_bln_c_s_d_2_s = lbound(fortran_obj%e_bln_c_s, dim=3)
    dace_rich_obj%f2dace_SA_e_flx_avg_d_0_s = size(fortran_obj%e_flx_avg, dim=1)
    dace_rich_obj%f2dace_SOA_e_flx_avg_d_0_s = lbound(fortran_obj%e_flx_avg, dim=1)
    dace_rich_obj%f2dace_SA_e_flx_avg_d_1_s = size(fortran_obj%e_flx_avg, dim=2)
    dace_rich_obj%f2dace_SOA_e_flx_avg_d_1_s = lbound(fortran_obj%e_flx_avg, dim=2)
    dace_rich_obj%f2dace_SA_e_flx_avg_d_2_s = size(fortran_obj%e_flx_avg, dim=3)
    dace_rich_obj%f2dace_SOA_e_flx_avg_d_2_s = lbound(fortran_obj%e_flx_avg, dim=3)
    dace_rich_obj%f2dace_SA_cells_aw_verts_d_0_s = size(fortran_obj%cells_aw_verts, dim=1)
    dace_rich_obj%f2dace_SOA_cells_aw_verts_d_0_s = lbound(fortran_obj%cells_aw_verts, dim=1)
    dace_rich_obj%f2dace_SA_cells_aw_verts_d_1_s = size(fortran_obj%cells_aw_verts, dim=2)
    dace_rich_obj%f2dace_SOA_cells_aw_verts_d_1_s = lbound(fortran_obj%cells_aw_verts, dim=2)
    dace_rich_obj%f2dace_SA_cells_aw_verts_d_2_s = size(fortran_obj%cells_aw_verts, dim=3)
    dace_rich_obj%f2dace_SOA_cells_aw_verts_d_2_s = lbound(fortran_obj%cells_aw_verts, dim=3)
    dace_rich_obj%f2dace_SA_rbf_vec_coeff_e_d_0_s = size(fortran_obj%rbf_vec_coeff_e, dim=1)
    dace_rich_obj%f2dace_SOA_rbf_vec_coeff_e_d_0_s = lbound(fortran_obj%rbf_vec_coeff_e, dim=1)
    dace_rich_obj%f2dace_SA_rbf_vec_coeff_e_d_1_s = size(fortran_obj%rbf_vec_coeff_e, dim=2)
    dace_rich_obj%f2dace_SOA_rbf_vec_coeff_e_d_1_s = lbound(fortran_obj%rbf_vec_coeff_e, dim=2)
    dace_rich_obj%f2dace_SA_rbf_vec_coeff_e_d_2_s = size(fortran_obj%rbf_vec_coeff_e, dim=3)
    dace_rich_obj%f2dace_SOA_rbf_vec_coeff_e_d_2_s = lbound(fortran_obj%rbf_vec_coeff_e, dim=3)
    dace_rich_obj%f2dace_SA_geofac_div_d_0_s = size(fortran_obj%geofac_div, dim=1)
    dace_rich_obj%f2dace_SOA_geofac_div_d_0_s = lbound(fortran_obj%geofac_div, dim=1)
    dace_rich_obj%f2dace_SA_geofac_div_d_1_s = size(fortran_obj%geofac_div, dim=2)
    dace_rich_obj%f2dace_SOA_geofac_div_d_1_s = lbound(fortran_obj%geofac_div, dim=2)
    dace_rich_obj%f2dace_SA_geofac_div_d_2_s = size(fortran_obj%geofac_div, dim=3)
    dace_rich_obj%f2dace_SOA_geofac_div_d_2_s = lbound(fortran_obj%geofac_div, dim=3)
    dace_rich_obj%f2dace_SA_geofac_grdiv_d_0_s = size(fortran_obj%geofac_grdiv, dim=1)
    dace_rich_obj%f2dace_SOA_geofac_grdiv_d_0_s = lbound(fortran_obj%geofac_grdiv, dim=1)
    dace_rich_obj%f2dace_SA_geofac_grdiv_d_1_s = size(fortran_obj%geofac_grdiv, dim=2)
    dace_rich_obj%f2dace_SOA_geofac_grdiv_d_1_s = lbound(fortran_obj%geofac_grdiv, dim=2)
    dace_rich_obj%f2dace_SA_geofac_grdiv_d_2_s = size(fortran_obj%geofac_grdiv, dim=3)
    dace_rich_obj%f2dace_SOA_geofac_grdiv_d_2_s = lbound(fortran_obj%geofac_grdiv, dim=3)
    dace_rich_obj%f2dace_SA_geofac_rot_d_0_s = size(fortran_obj%geofac_rot, dim=1)
    dace_rich_obj%f2dace_SOA_geofac_rot_d_0_s = lbound(fortran_obj%geofac_rot, dim=1)
    dace_rich_obj%f2dace_SA_geofac_rot_d_1_s = size(fortran_obj%geofac_rot, dim=2)
    dace_rich_obj%f2dace_SOA_geofac_rot_d_1_s = lbound(fortran_obj%geofac_rot, dim=2)
    dace_rich_obj%f2dace_SA_geofac_rot_d_2_s = size(fortran_obj%geofac_rot, dim=3)
    dace_rich_obj%f2dace_SOA_geofac_rot_d_2_s = lbound(fortran_obj%geofac_rot, dim=3)
    dace_rich_obj%f2dace_SA_geofac_n2s_d_0_s = size(fortran_obj%geofac_n2s, dim=1)
    dace_rich_obj%f2dace_SOA_geofac_n2s_d_0_s = lbound(fortran_obj%geofac_n2s, dim=1)
    dace_rich_obj%f2dace_SA_geofac_n2s_d_1_s = size(fortran_obj%geofac_n2s, dim=2)
    dace_rich_obj%f2dace_SOA_geofac_n2s_d_1_s = lbound(fortran_obj%geofac_n2s, dim=2)
    dace_rich_obj%f2dace_SA_geofac_n2s_d_2_s = size(fortran_obj%geofac_n2s, dim=3)
    dace_rich_obj%f2dace_SOA_geofac_n2s_d_2_s = lbound(fortran_obj%geofac_n2s, dim=3)
    dace_rich_obj%f2dace_SA_geofac_grg_d_0_s = size(fortran_obj%geofac_grg, dim=1)
    dace_rich_obj%f2dace_SOA_geofac_grg_d_0_s = lbound(fortran_obj%geofac_grg, dim=1)
    dace_rich_obj%f2dace_SA_geofac_grg_d_1_s = size(fortran_obj%geofac_grg, dim=2)
    dace_rich_obj%f2dace_SOA_geofac_grg_d_1_s = lbound(fortran_obj%geofac_grg, dim=2)
    dace_rich_obj%f2dace_SA_geofac_grg_d_2_s = size(fortran_obj%geofac_grg, dim=3)
    dace_rich_obj%f2dace_SOA_geofac_grg_d_2_s = lbound(fortran_obj%geofac_grg, dim=3)
    dace_rich_obj%f2dace_SA_geofac_grg_d_3_s = size(fortran_obj%geofac_grg, dim=4)
    dace_rich_obj%f2dace_SOA_geofac_grg_d_3_s = lbound(fortran_obj%geofac_grg, dim=4)
    dace_rich_obj%f2dace_SA_pos_on_tplane_e_d_0_s = size(fortran_obj%pos_on_tplane_e, dim=1)
    dace_rich_obj%f2dace_SOA_pos_on_tplane_e_d_0_s = lbound(fortran_obj%pos_on_tplane_e, dim=1)
    dace_rich_obj%f2dace_SA_pos_on_tplane_e_d_1_s = size(fortran_obj%pos_on_tplane_e, dim=2)
    dace_rich_obj%f2dace_SOA_pos_on_tplane_e_d_1_s = lbound(fortran_obj%pos_on_tplane_e, dim=2)
    dace_rich_obj%f2dace_SA_pos_on_tplane_e_d_2_s = size(fortran_obj%pos_on_tplane_e, dim=3)
    dace_rich_obj%f2dace_SOA_pos_on_tplane_e_d_2_s = lbound(fortran_obj%pos_on_tplane_e, dim=3)
    dace_rich_obj%f2dace_SA_pos_on_tplane_e_d_3_s = size(fortran_obj%pos_on_tplane_e, dim=4)
    dace_rich_obj%f2dace_SOA_pos_on_tplane_e_d_3_s = lbound(fortran_obj%pos_on_tplane_e, dim=4)
    dace_rich_obj%f2dace_SA_nudgecoeff_e_d_0_s = size(fortran_obj%nudgecoeff_e, dim=1)
    dace_rich_obj%f2dace_SOA_nudgecoeff_e_d_0_s = lbound(fortran_obj%nudgecoeff_e, dim=1)
    dace_rich_obj%f2dace_SA_nudgecoeff_e_d_1_s = size(fortran_obj%nudgecoeff_e, dim=2)
    dace_rich_obj%f2dace_SOA_nudgecoeff_e_d_1_s = lbound(fortran_obj%nudgecoeff_e, dim=2)
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
    dace_rich_obj%e_flx_avg = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%e_flx_avg, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%e_flx_avg = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%e_flx_avg, &
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
    dace_rich_obj%geofac_div = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%geofac_div, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%geofac_div = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%geofac_div, &
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
#ifndef _OPENACC
    dace_rich_obj%geofac_grg = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%geofac_grg, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%geofac_grg = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%geofac_grg, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%pos_on_tplane_e = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%pos_on_tplane_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%pos_on_tplane_e = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%pos_on_tplane_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%nudgecoeff_e = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%nudgecoeff_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%nudgecoeff_e = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%nudgecoeff_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif

  end function copy_in_t_int_state

  function copy_in_t_nh_state(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_nh_state), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_nh_state), pointer :: dace_rich_obj
    type(dace_t_nh_state) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%diag = copy_in_t_nh_diag( &
    fortran_obj=fortran_obj%diag, &
    steal_arrays=steal_arrays, &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%ref = copy_in_t_nh_ref( &
    fortran_obj=fortran_obj%ref, &
    steal_arrays=steal_arrays, &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%metrics = copy_in_t_nh_metrics( &
    fortran_obj=fortran_obj%metrics, &
    steal_arrays=steal_arrays, &
    minimal_structs=minimal_structs &
  )

  end function copy_in_t_nh_state

  function copy_in_t_nh_diag(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_nh_diag), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_nh_diag), pointer :: dace_rich_obj
    type(dace_t_nh_diag) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_exner_pr_d_0_s = size(fortran_obj%exner_pr, dim=1)
    dace_rich_obj%f2dace_SOA_exner_pr_d_0_s = lbound(fortran_obj%exner_pr, dim=1)
    dace_rich_obj%f2dace_SA_exner_pr_d_1_s = size(fortran_obj%exner_pr, dim=2)
    dace_rich_obj%f2dace_SOA_exner_pr_d_1_s = lbound(fortran_obj%exner_pr, dim=2)
    dace_rich_obj%f2dace_SA_exner_pr_d_2_s = size(fortran_obj%exner_pr, dim=3)
    dace_rich_obj%f2dace_SOA_exner_pr_d_2_s = lbound(fortran_obj%exner_pr, dim=3)
    dace_rich_obj%f2dace_SA_mass_fl_e_d_0_s = size(fortran_obj%mass_fl_e, dim=1)
    dace_rich_obj%f2dace_SOA_mass_fl_e_d_0_s = lbound(fortran_obj%mass_fl_e, dim=1)
    dace_rich_obj%f2dace_SA_mass_fl_e_d_1_s = size(fortran_obj%mass_fl_e, dim=2)
    dace_rich_obj%f2dace_SOA_mass_fl_e_d_1_s = lbound(fortran_obj%mass_fl_e, dim=2)
    dace_rich_obj%f2dace_SA_mass_fl_e_d_2_s = size(fortran_obj%mass_fl_e, dim=3)
    dace_rich_obj%f2dace_SOA_mass_fl_e_d_2_s = lbound(fortran_obj%mass_fl_e, dim=3)
    dace_rich_obj%f2dace_SA_rho_ic_d_0_s = size(fortran_obj%rho_ic, dim=1)
    dace_rich_obj%f2dace_SOA_rho_ic_d_0_s = lbound(fortran_obj%rho_ic, dim=1)
    dace_rich_obj%f2dace_SA_rho_ic_d_1_s = size(fortran_obj%rho_ic, dim=2)
    dace_rich_obj%f2dace_SOA_rho_ic_d_1_s = lbound(fortran_obj%rho_ic, dim=2)
    dace_rich_obj%f2dace_SA_rho_ic_d_2_s = size(fortran_obj%rho_ic, dim=3)
    dace_rich_obj%f2dace_SOA_rho_ic_d_2_s = lbound(fortran_obj%rho_ic, dim=3)
    dace_rich_obj%f2dace_SA_theta_v_ic_d_0_s = size(fortran_obj%theta_v_ic, dim=1)
    dace_rich_obj%f2dace_SOA_theta_v_ic_d_0_s = lbound(fortran_obj%theta_v_ic, dim=1)
    dace_rich_obj%f2dace_SA_theta_v_ic_d_1_s = size(fortran_obj%theta_v_ic, dim=2)
    dace_rich_obj%f2dace_SOA_theta_v_ic_d_1_s = lbound(fortran_obj%theta_v_ic, dim=2)
    dace_rich_obj%f2dace_SA_theta_v_ic_d_2_s = size(fortran_obj%theta_v_ic, dim=3)
    dace_rich_obj%f2dace_SOA_theta_v_ic_d_2_s = lbound(fortran_obj%theta_v_ic, dim=3)
    dace_rich_obj%f2dace_SA_grf_tend_vn_d_0_s = size(fortran_obj%grf_tend_vn, dim=1)
    dace_rich_obj%f2dace_SOA_grf_tend_vn_d_0_s = lbound(fortran_obj%grf_tend_vn, dim=1)
    dace_rich_obj%f2dace_SA_grf_tend_vn_d_1_s = size(fortran_obj%grf_tend_vn, dim=2)
    dace_rich_obj%f2dace_SOA_grf_tend_vn_d_1_s = lbound(fortran_obj%grf_tend_vn, dim=2)
    dace_rich_obj%f2dace_SA_grf_tend_vn_d_2_s = size(fortran_obj%grf_tend_vn, dim=3)
    dace_rich_obj%f2dace_SOA_grf_tend_vn_d_2_s = lbound(fortran_obj%grf_tend_vn, dim=3)
    dace_rich_obj%f2dace_SA_grf_tend_w_d_0_s = size(fortran_obj%grf_tend_w, dim=1)
    dace_rich_obj%f2dace_SOA_grf_tend_w_d_0_s = lbound(fortran_obj%grf_tend_w, dim=1)
    dace_rich_obj%f2dace_SA_grf_tend_w_d_1_s = size(fortran_obj%grf_tend_w, dim=2)
    dace_rich_obj%f2dace_SOA_grf_tend_w_d_1_s = lbound(fortran_obj%grf_tend_w, dim=2)
    dace_rich_obj%f2dace_SA_grf_tend_w_d_2_s = size(fortran_obj%grf_tend_w, dim=3)
    dace_rich_obj%f2dace_SOA_grf_tend_w_d_2_s = lbound(fortran_obj%grf_tend_w, dim=3)
    dace_rich_obj%f2dace_SA_grf_tend_rho_d_0_s = size(fortran_obj%grf_tend_rho, dim=1)
    dace_rich_obj%f2dace_SOA_grf_tend_rho_d_0_s = lbound(fortran_obj%grf_tend_rho, dim=1)
    dace_rich_obj%f2dace_SA_grf_tend_rho_d_1_s = size(fortran_obj%grf_tend_rho, dim=2)
    dace_rich_obj%f2dace_SOA_grf_tend_rho_d_1_s = lbound(fortran_obj%grf_tend_rho, dim=2)
    dace_rich_obj%f2dace_SA_grf_tend_rho_d_2_s = size(fortran_obj%grf_tend_rho, dim=3)
    dace_rich_obj%f2dace_SOA_grf_tend_rho_d_2_s = lbound(fortran_obj%grf_tend_rho, dim=3)
    dace_rich_obj%f2dace_SA_grf_tend_mflx_d_0_s = size(fortran_obj%grf_tend_mflx, dim=1)
    dace_rich_obj%f2dace_SOA_grf_tend_mflx_d_0_s = lbound(fortran_obj%grf_tend_mflx, dim=1)
    dace_rich_obj%f2dace_SA_grf_tend_mflx_d_1_s = size(fortran_obj%grf_tend_mflx, dim=2)
    dace_rich_obj%f2dace_SOA_grf_tend_mflx_d_1_s = lbound(fortran_obj%grf_tend_mflx, dim=2)
    dace_rich_obj%f2dace_SA_grf_tend_mflx_d_2_s = size(fortran_obj%grf_tend_mflx, dim=3)
    dace_rich_obj%f2dace_SOA_grf_tend_mflx_d_2_s = lbound(fortran_obj%grf_tend_mflx, dim=3)
    dace_rich_obj%f2dace_SA_grf_bdy_mflx_d_0_s = size(fortran_obj%grf_bdy_mflx, dim=1)
    dace_rich_obj%f2dace_SOA_grf_bdy_mflx_d_0_s = lbound(fortran_obj%grf_bdy_mflx, dim=1)
    dace_rich_obj%f2dace_SA_grf_bdy_mflx_d_1_s = size(fortran_obj%grf_bdy_mflx, dim=2)
    dace_rich_obj%f2dace_SOA_grf_bdy_mflx_d_1_s = lbound(fortran_obj%grf_bdy_mflx, dim=2)
    dace_rich_obj%f2dace_SA_grf_bdy_mflx_d_2_s = size(fortran_obj%grf_bdy_mflx, dim=3)
    dace_rich_obj%f2dace_SOA_grf_bdy_mflx_d_2_s = lbound(fortran_obj%grf_bdy_mflx, dim=3)
    dace_rich_obj%f2dace_SA_grf_tend_thv_d_0_s = size(fortran_obj%grf_tend_thv, dim=1)
    dace_rich_obj%f2dace_SOA_grf_tend_thv_d_0_s = lbound(fortran_obj%grf_tend_thv, dim=1)
    dace_rich_obj%f2dace_SA_grf_tend_thv_d_1_s = size(fortran_obj%grf_tend_thv, dim=2)
    dace_rich_obj%f2dace_SOA_grf_tend_thv_d_1_s = lbound(fortran_obj%grf_tend_thv, dim=2)
    dace_rich_obj%f2dace_SA_grf_tend_thv_d_2_s = size(fortran_obj%grf_tend_thv, dim=3)
    dace_rich_obj%f2dace_SOA_grf_tend_thv_d_2_s = lbound(fortran_obj%grf_tend_thv, dim=3)
    dace_rich_obj%f2dace_SA_vn_ie_int_d_0_s = size(fortran_obj%vn_ie_int, dim=1)
    dace_rich_obj%f2dace_SOA_vn_ie_int_d_0_s = lbound(fortran_obj%vn_ie_int, dim=1)
    dace_rich_obj%f2dace_SA_vn_ie_int_d_1_s = size(fortran_obj%vn_ie_int, dim=2)
    dace_rich_obj%f2dace_SOA_vn_ie_int_d_1_s = lbound(fortran_obj%vn_ie_int, dim=2)
    dace_rich_obj%f2dace_SA_vn_ie_int_d_2_s = size(fortran_obj%vn_ie_int, dim=3)
    dace_rich_obj%f2dace_SOA_vn_ie_int_d_2_s = lbound(fortran_obj%vn_ie_int, dim=3)
    dace_rich_obj%f2dace_SA_vn_ie_ubc_d_0_s = size(fortran_obj%vn_ie_ubc, dim=1)
    dace_rich_obj%f2dace_SOA_vn_ie_ubc_d_0_s = lbound(fortran_obj%vn_ie_ubc, dim=1)
    dace_rich_obj%f2dace_SA_vn_ie_ubc_d_1_s = size(fortran_obj%vn_ie_ubc, dim=2)
    dace_rich_obj%f2dace_SOA_vn_ie_ubc_d_1_s = lbound(fortran_obj%vn_ie_ubc, dim=2)
    dace_rich_obj%f2dace_SA_vn_ie_ubc_d_2_s = size(fortran_obj%vn_ie_ubc, dim=3)
    dace_rich_obj%f2dace_SOA_vn_ie_ubc_d_2_s = lbound(fortran_obj%vn_ie_ubc, dim=3)
    dace_rich_obj%f2dace_SA_w_int_d_0_s = size(fortran_obj%w_int, dim=1)
    dace_rich_obj%f2dace_SOA_w_int_d_0_s = lbound(fortran_obj%w_int, dim=1)
    dace_rich_obj%f2dace_SA_w_int_d_1_s = size(fortran_obj%w_int, dim=2)
    dace_rich_obj%f2dace_SOA_w_int_d_1_s = lbound(fortran_obj%w_int, dim=2)
    dace_rich_obj%f2dace_SA_w_int_d_2_s = size(fortran_obj%w_int, dim=3)
    dace_rich_obj%f2dace_SOA_w_int_d_2_s = lbound(fortran_obj%w_int, dim=3)
    dace_rich_obj%f2dace_SA_w_ubc_d_0_s = size(fortran_obj%w_ubc, dim=1)
    dace_rich_obj%f2dace_SOA_w_ubc_d_0_s = lbound(fortran_obj%w_ubc, dim=1)
    dace_rich_obj%f2dace_SA_w_ubc_d_1_s = size(fortran_obj%w_ubc, dim=2)
    dace_rich_obj%f2dace_SOA_w_ubc_d_1_s = lbound(fortran_obj%w_ubc, dim=2)
    dace_rich_obj%f2dace_SA_w_ubc_d_2_s = size(fortran_obj%w_ubc, dim=3)
    dace_rich_obj%f2dace_SOA_w_ubc_d_2_s = lbound(fortran_obj%w_ubc, dim=3)
    dace_rich_obj%f2dace_SA_theta_v_ic_int_d_0_s = size(fortran_obj%theta_v_ic_int, dim=1)
    dace_rich_obj%f2dace_SOA_theta_v_ic_int_d_0_s = lbound(fortran_obj%theta_v_ic_int, dim=1)
    dace_rich_obj%f2dace_SA_theta_v_ic_int_d_1_s = size(fortran_obj%theta_v_ic_int, dim=2)
    dace_rich_obj%f2dace_SOA_theta_v_ic_int_d_1_s = lbound(fortran_obj%theta_v_ic_int, dim=2)
    dace_rich_obj%f2dace_SA_theta_v_ic_int_d_2_s = size(fortran_obj%theta_v_ic_int, dim=3)
    dace_rich_obj%f2dace_SOA_theta_v_ic_int_d_2_s = lbound(fortran_obj%theta_v_ic_int, dim=3)
    dace_rich_obj%f2dace_SA_theta_v_ic_ubc_d_0_s = size(fortran_obj%theta_v_ic_ubc, dim=1)
    dace_rich_obj%f2dace_SOA_theta_v_ic_ubc_d_0_s = lbound(fortran_obj%theta_v_ic_ubc, dim=1)
    dace_rich_obj%f2dace_SA_theta_v_ic_ubc_d_1_s = size(fortran_obj%theta_v_ic_ubc, dim=2)
    dace_rich_obj%f2dace_SOA_theta_v_ic_ubc_d_1_s = lbound(fortran_obj%theta_v_ic_ubc, dim=2)
    dace_rich_obj%f2dace_SA_theta_v_ic_ubc_d_2_s = size(fortran_obj%theta_v_ic_ubc, dim=3)
    dace_rich_obj%f2dace_SOA_theta_v_ic_ubc_d_2_s = lbound(fortran_obj%theta_v_ic_ubc, dim=3)
    dace_rich_obj%f2dace_SA_rho_ic_int_d_0_s = size(fortran_obj%rho_ic_int, dim=1)
    dace_rich_obj%f2dace_SOA_rho_ic_int_d_0_s = lbound(fortran_obj%rho_ic_int, dim=1)
    dace_rich_obj%f2dace_SA_rho_ic_int_d_1_s = size(fortran_obj%rho_ic_int, dim=2)
    dace_rich_obj%f2dace_SOA_rho_ic_int_d_1_s = lbound(fortran_obj%rho_ic_int, dim=2)
    dace_rich_obj%f2dace_SA_rho_ic_int_d_2_s = size(fortran_obj%rho_ic_int, dim=3)
    dace_rich_obj%f2dace_SOA_rho_ic_int_d_2_s = lbound(fortran_obj%rho_ic_int, dim=3)
    dace_rich_obj%f2dace_SA_rho_ic_ubc_d_0_s = size(fortran_obj%rho_ic_ubc, dim=1)
    dace_rich_obj%f2dace_SOA_rho_ic_ubc_d_0_s = lbound(fortran_obj%rho_ic_ubc, dim=1)
    dace_rich_obj%f2dace_SA_rho_ic_ubc_d_1_s = size(fortran_obj%rho_ic_ubc, dim=2)
    dace_rich_obj%f2dace_SOA_rho_ic_ubc_d_1_s = lbound(fortran_obj%rho_ic_ubc, dim=2)
    dace_rich_obj%f2dace_SA_rho_ic_ubc_d_2_s = size(fortran_obj%rho_ic_ubc, dim=3)
    dace_rich_obj%f2dace_SOA_rho_ic_ubc_d_2_s = lbound(fortran_obj%rho_ic_ubc, dim=3)
    dace_rich_obj%f2dace_SA_mflx_ic_int_d_0_s = size(fortran_obj%mflx_ic_int, dim=1)
    dace_rich_obj%f2dace_SOA_mflx_ic_int_d_0_s = lbound(fortran_obj%mflx_ic_int, dim=1)
    dace_rich_obj%f2dace_SA_mflx_ic_int_d_1_s = size(fortran_obj%mflx_ic_int, dim=2)
    dace_rich_obj%f2dace_SOA_mflx_ic_int_d_1_s = lbound(fortran_obj%mflx_ic_int, dim=2)
    dace_rich_obj%f2dace_SA_mflx_ic_int_d_2_s = size(fortran_obj%mflx_ic_int, dim=3)
    dace_rich_obj%f2dace_SOA_mflx_ic_int_d_2_s = lbound(fortran_obj%mflx_ic_int, dim=3)
    dace_rich_obj%f2dace_SA_mflx_ic_ubc_d_0_s = size(fortran_obj%mflx_ic_ubc, dim=1)
    dace_rich_obj%f2dace_SOA_mflx_ic_ubc_d_0_s = lbound(fortran_obj%mflx_ic_ubc, dim=1)
    dace_rich_obj%f2dace_SA_mflx_ic_ubc_d_1_s = size(fortran_obj%mflx_ic_ubc, dim=2)
    dace_rich_obj%f2dace_SOA_mflx_ic_ubc_d_1_s = lbound(fortran_obj%mflx_ic_ubc, dim=2)
    dace_rich_obj%f2dace_SA_mflx_ic_ubc_d_2_s = size(fortran_obj%mflx_ic_ubc, dim=3)
    dace_rich_obj%f2dace_SOA_mflx_ic_ubc_d_2_s = lbound(fortran_obj%mflx_ic_ubc, dim=3)
    dace_rich_obj%f2dace_SA_vn_incr_d_0_s = size(fortran_obj%vn_incr, dim=1)
    dace_rich_obj%f2dace_SOA_vn_incr_d_0_s = lbound(fortran_obj%vn_incr, dim=1)
    dace_rich_obj%f2dace_SA_vn_incr_d_1_s = size(fortran_obj%vn_incr, dim=2)
    dace_rich_obj%f2dace_SOA_vn_incr_d_1_s = lbound(fortran_obj%vn_incr, dim=2)
    dace_rich_obj%f2dace_SA_vn_incr_d_2_s = size(fortran_obj%vn_incr, dim=3)
    dace_rich_obj%f2dace_SOA_vn_incr_d_2_s = lbound(fortran_obj%vn_incr, dim=3)
    dace_rich_obj%f2dace_SA_exner_incr_d_0_s = size(fortran_obj%exner_incr, dim=1)
    dace_rich_obj%f2dace_SOA_exner_incr_d_0_s = lbound(fortran_obj%exner_incr, dim=1)
    dace_rich_obj%f2dace_SA_exner_incr_d_1_s = size(fortran_obj%exner_incr, dim=2)
    dace_rich_obj%f2dace_SOA_exner_incr_d_1_s = lbound(fortran_obj%exner_incr, dim=2)
    dace_rich_obj%f2dace_SA_exner_incr_d_2_s = size(fortran_obj%exner_incr, dim=3)
    dace_rich_obj%f2dace_SOA_exner_incr_d_2_s = lbound(fortran_obj%exner_incr, dim=3)
    dace_rich_obj%f2dace_SA_rho_incr_d_0_s = size(fortran_obj%rho_incr, dim=1)
    dace_rich_obj%f2dace_SOA_rho_incr_d_0_s = lbound(fortran_obj%rho_incr, dim=1)
    dace_rich_obj%f2dace_SA_rho_incr_d_1_s = size(fortran_obj%rho_incr, dim=2)
    dace_rich_obj%f2dace_SOA_rho_incr_d_1_s = lbound(fortran_obj%rho_incr, dim=2)
    dace_rich_obj%f2dace_SA_rho_incr_d_2_s = size(fortran_obj%rho_incr, dim=3)
    dace_rich_obj%f2dace_SOA_rho_incr_d_2_s = lbound(fortran_obj%rho_incr, dim=3)
    dace_rich_obj%f2dace_SA_vt_d_0_s = size(fortran_obj%vt, dim=1)
    dace_rich_obj%f2dace_SOA_vt_d_0_s = lbound(fortran_obj%vt, dim=1)
    dace_rich_obj%f2dace_SA_vt_d_1_s = size(fortran_obj%vt, dim=2)
    dace_rich_obj%f2dace_SOA_vt_d_1_s = lbound(fortran_obj%vt, dim=2)
    dace_rich_obj%f2dace_SA_vt_d_2_s = size(fortran_obj%vt, dim=3)
    dace_rich_obj%f2dace_SOA_vt_d_2_s = lbound(fortran_obj%vt, dim=3)
    dace_rich_obj%f2dace_SA_ddt_exner_phy_d_0_s = size(fortran_obj%ddt_exner_phy, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_exner_phy_d_0_s = lbound(fortran_obj%ddt_exner_phy, dim=1)
    dace_rich_obj%f2dace_SA_ddt_exner_phy_d_1_s = size(fortran_obj%ddt_exner_phy, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_exner_phy_d_1_s = lbound(fortran_obj%ddt_exner_phy, dim=2)
    dace_rich_obj%f2dace_SA_ddt_exner_phy_d_2_s = size(fortran_obj%ddt_exner_phy, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_exner_phy_d_2_s = lbound(fortran_obj%ddt_exner_phy, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_phy_d_0_s = size(fortran_obj%ddt_vn_phy, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_phy_d_0_s = lbound(fortran_obj%ddt_vn_phy, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_phy_d_1_s = size(fortran_obj%ddt_vn_phy, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_phy_d_1_s = lbound(fortran_obj%ddt_vn_phy, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_phy_d_2_s = size(fortran_obj%ddt_vn_phy, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_phy_d_2_s = lbound(fortran_obj%ddt_vn_phy, dim=3)
    dace_rich_obj%f2dace_SA_exner_dyn_incr_d_0_s = size(fortran_obj%exner_dyn_incr, dim=1)
    dace_rich_obj%f2dace_SOA_exner_dyn_incr_d_0_s = lbound(fortran_obj%exner_dyn_incr, dim=1)
    dace_rich_obj%f2dace_SA_exner_dyn_incr_d_1_s = size(fortran_obj%exner_dyn_incr, dim=2)
    dace_rich_obj%f2dace_SOA_exner_dyn_incr_d_1_s = lbound(fortran_obj%exner_dyn_incr, dim=2)
    dace_rich_obj%f2dace_SA_exner_dyn_incr_d_2_s = size(fortran_obj%exner_dyn_incr, dim=3)
    dace_rich_obj%f2dace_SOA_exner_dyn_incr_d_2_s = lbound(fortran_obj%exner_dyn_incr, dim=3)
    dace_rich_obj%f2dace_SA_vn_ie_d_0_s = size(fortran_obj%vn_ie, dim=1)
    dace_rich_obj%f2dace_SOA_vn_ie_d_0_s = lbound(fortran_obj%vn_ie, dim=1)
    dace_rich_obj%f2dace_SA_vn_ie_d_1_s = size(fortran_obj%vn_ie, dim=2)
    dace_rich_obj%f2dace_SOA_vn_ie_d_1_s = lbound(fortran_obj%vn_ie, dim=2)
    dace_rich_obj%f2dace_SA_vn_ie_d_2_s = size(fortran_obj%vn_ie, dim=3)
    dace_rich_obj%f2dace_SOA_vn_ie_d_2_s = lbound(fortran_obj%vn_ie, dim=3)
    dace_rich_obj%f2dace_SA_w_concorr_c_d_0_s = size(fortran_obj%w_concorr_c, dim=1)
    dace_rich_obj%f2dace_SOA_w_concorr_c_d_0_s = lbound(fortran_obj%w_concorr_c, dim=1)
    dace_rich_obj%f2dace_SA_w_concorr_c_d_1_s = size(fortran_obj%w_concorr_c, dim=2)
    dace_rich_obj%f2dace_SOA_w_concorr_c_d_1_s = lbound(fortran_obj%w_concorr_c, dim=2)
    dace_rich_obj%f2dace_SA_w_concorr_c_d_2_s = size(fortran_obj%w_concorr_c, dim=3)
    dace_rich_obj%f2dace_SOA_w_concorr_c_d_2_s = lbound(fortran_obj%w_concorr_c, dim=3)
    dace_rich_obj%f2dace_SA_mass_fl_e_sv_d_0_s = size(fortran_obj%mass_fl_e_sv, dim=1)
    dace_rich_obj%f2dace_SOA_mass_fl_e_sv_d_0_s = lbound(fortran_obj%mass_fl_e_sv, dim=1)
    dace_rich_obj%f2dace_SA_mass_fl_e_sv_d_1_s = size(fortran_obj%mass_fl_e_sv, dim=2)
    dace_rich_obj%f2dace_SOA_mass_fl_e_sv_d_1_s = lbound(fortran_obj%mass_fl_e_sv, dim=2)
    dace_rich_obj%f2dace_SA_mass_fl_e_sv_d_2_s = size(fortran_obj%mass_fl_e_sv, dim=3)
    dace_rich_obj%f2dace_SOA_mass_fl_e_sv_d_2_s = lbound(fortran_obj%mass_fl_e_sv, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_apc_pc_d_0_s = size(fortran_obj%ddt_vn_apc_pc, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_apc_pc_d_0_s = lbound(fortran_obj%ddt_vn_apc_pc, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_apc_pc_d_1_s = size(fortran_obj%ddt_vn_apc_pc, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_apc_pc_d_1_s = lbound(fortran_obj%ddt_vn_apc_pc, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_apc_pc_d_2_s = size(fortran_obj%ddt_vn_apc_pc, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_apc_pc_d_2_s = lbound(fortran_obj%ddt_vn_apc_pc, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_apc_pc_d_3_s = size(fortran_obj%ddt_vn_apc_pc, dim=4)
    dace_rich_obj%f2dace_SOA_ddt_vn_apc_pc_d_3_s = lbound(fortran_obj%ddt_vn_apc_pc, dim=4)
    dace_rich_obj%f2dace_SA_ddt_vn_cor_pc_d_0_s = size(fortran_obj%ddt_vn_cor_pc, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_cor_pc_d_0_s = lbound(fortran_obj%ddt_vn_cor_pc, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_cor_pc_d_1_s = size(fortran_obj%ddt_vn_cor_pc, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_cor_pc_d_1_s = lbound(fortran_obj%ddt_vn_cor_pc, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_cor_pc_d_2_s = size(fortran_obj%ddt_vn_cor_pc, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_cor_pc_d_2_s = lbound(fortran_obj%ddt_vn_cor_pc, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_cor_pc_d_3_s = size(fortran_obj%ddt_vn_cor_pc, dim=4)
    dace_rich_obj%f2dace_SOA_ddt_vn_cor_pc_d_3_s = lbound(fortran_obj%ddt_vn_cor_pc, dim=4)
    dace_rich_obj%f2dace_SA_ddt_w_adv_pc_d_0_s = size(fortran_obj%ddt_w_adv_pc, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_w_adv_pc_d_0_s = lbound(fortran_obj%ddt_w_adv_pc, dim=1)
    dace_rich_obj%f2dace_SA_ddt_w_adv_pc_d_1_s = size(fortran_obj%ddt_w_adv_pc, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_w_adv_pc_d_1_s = lbound(fortran_obj%ddt_w_adv_pc, dim=2)
    dace_rich_obj%f2dace_SA_ddt_w_adv_pc_d_2_s = size(fortran_obj%ddt_w_adv_pc, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_w_adv_pc_d_2_s = lbound(fortran_obj%ddt_w_adv_pc, dim=3)
    dace_rich_obj%f2dace_SA_ddt_w_adv_pc_d_3_s = size(fortran_obj%ddt_w_adv_pc, dim=4)
    dace_rich_obj%f2dace_SOA_ddt_w_adv_pc_d_3_s = lbound(fortran_obj%ddt_w_adv_pc, dim=4)
    dace_rich_obj%f2dace_SA_ddt_vn_dyn_d_0_s = size(fortran_obj%ddt_vn_dyn, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_dyn_d_0_s = lbound(fortran_obj%ddt_vn_dyn, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_dyn_d_1_s = size(fortran_obj%ddt_vn_dyn, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_dyn_d_1_s = lbound(fortran_obj%ddt_vn_dyn, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_dyn_d_2_s = size(fortran_obj%ddt_vn_dyn, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_dyn_d_2_s = lbound(fortran_obj%ddt_vn_dyn, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_dmp_d_0_s = size(fortran_obj%ddt_vn_dmp, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_dmp_d_0_s = lbound(fortran_obj%ddt_vn_dmp, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_dmp_d_1_s = size(fortran_obj%ddt_vn_dmp, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_dmp_d_1_s = lbound(fortran_obj%ddt_vn_dmp, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_dmp_d_2_s = size(fortran_obj%ddt_vn_dmp, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_dmp_d_2_s = lbound(fortran_obj%ddt_vn_dmp, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_adv_d_0_s = size(fortran_obj%ddt_vn_adv, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_adv_d_0_s = lbound(fortran_obj%ddt_vn_adv, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_adv_d_1_s = size(fortran_obj%ddt_vn_adv, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_adv_d_1_s = lbound(fortran_obj%ddt_vn_adv, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_adv_d_2_s = size(fortran_obj%ddt_vn_adv, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_adv_d_2_s = lbound(fortran_obj%ddt_vn_adv, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_cor_d_0_s = size(fortran_obj%ddt_vn_cor, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_cor_d_0_s = lbound(fortran_obj%ddt_vn_cor, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_cor_d_1_s = size(fortran_obj%ddt_vn_cor, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_cor_d_1_s = lbound(fortran_obj%ddt_vn_cor, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_cor_d_2_s = size(fortran_obj%ddt_vn_cor, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_cor_d_2_s = lbound(fortran_obj%ddt_vn_cor, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_pgr_d_0_s = size(fortran_obj%ddt_vn_pgr, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_pgr_d_0_s = lbound(fortran_obj%ddt_vn_pgr, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_pgr_d_1_s = size(fortran_obj%ddt_vn_pgr, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_pgr_d_1_s = lbound(fortran_obj%ddt_vn_pgr, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_pgr_d_2_s = size(fortran_obj%ddt_vn_pgr, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_pgr_d_2_s = lbound(fortran_obj%ddt_vn_pgr, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_phd_d_0_s = size(fortran_obj%ddt_vn_phd, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_phd_d_0_s = lbound(fortran_obj%ddt_vn_phd, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_phd_d_1_s = size(fortran_obj%ddt_vn_phd, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_phd_d_1_s = lbound(fortran_obj%ddt_vn_phd, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_phd_d_2_s = size(fortran_obj%ddt_vn_phd, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_phd_d_2_s = lbound(fortran_obj%ddt_vn_phd, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_iau_d_0_s = size(fortran_obj%ddt_vn_iau, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_iau_d_0_s = lbound(fortran_obj%ddt_vn_iau, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_iau_d_1_s = size(fortran_obj%ddt_vn_iau, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_iau_d_1_s = lbound(fortran_obj%ddt_vn_iau, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_iau_d_2_s = size(fortran_obj%ddt_vn_iau, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_iau_d_2_s = lbound(fortran_obj%ddt_vn_iau, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_ray_d_0_s = size(fortran_obj%ddt_vn_ray, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_ray_d_0_s = lbound(fortran_obj%ddt_vn_ray, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_ray_d_1_s = size(fortran_obj%ddt_vn_ray, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_ray_d_1_s = lbound(fortran_obj%ddt_vn_ray, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_ray_d_2_s = size(fortran_obj%ddt_vn_ray, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_ray_d_2_s = lbound(fortran_obj%ddt_vn_ray, dim=3)
    dace_rich_obj%f2dace_SA_ddt_vn_grf_d_0_s = size(fortran_obj%ddt_vn_grf, dim=1)
    dace_rich_obj%f2dace_SOA_ddt_vn_grf_d_0_s = lbound(fortran_obj%ddt_vn_grf, dim=1)
    dace_rich_obj%f2dace_SA_ddt_vn_grf_d_1_s = size(fortran_obj%ddt_vn_grf, dim=2)
    dace_rich_obj%f2dace_SOA_ddt_vn_grf_d_1_s = lbound(fortran_obj%ddt_vn_grf, dim=2)
    dace_rich_obj%f2dace_SA_ddt_vn_grf_d_2_s = size(fortran_obj%ddt_vn_grf, dim=3)
    dace_rich_obj%f2dace_SOA_ddt_vn_grf_d_2_s = lbound(fortran_obj%ddt_vn_grf, dim=3)
    dace_rich_obj%ddt_vn_dyn_is_associated = fortran_obj%ddt_vn_dyn_is_associated
    dace_rich_obj%ddt_vn_dmp_is_associated = fortran_obj%ddt_vn_dmp_is_associated
    dace_rich_obj%ddt_vn_adv_is_associated = fortran_obj%ddt_vn_adv_is_associated
    dace_rich_obj%ddt_vn_cor_is_associated = fortran_obj%ddt_vn_cor_is_associated
    dace_rich_obj%ddt_vn_pgr_is_associated = fortran_obj%ddt_vn_pgr_is_associated
    dace_rich_obj%ddt_vn_phd_is_associated = fortran_obj%ddt_vn_phd_is_associated
    dace_rich_obj%ddt_vn_iau_is_associated = fortran_obj%ddt_vn_iau_is_associated
    dace_rich_obj%ddt_vn_ray_is_associated = fortran_obj%ddt_vn_ray_is_associated
    dace_rich_obj%ddt_vn_grf_is_associated = fortran_obj%ddt_vn_grf_is_associated
    dace_rich_obj%max_vcfl_dyn = fortran_obj%max_vcfl_dyn
#ifndef _OPENACC
    dace_rich_obj%exner_pr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner_pr, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%exner_pr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner_pr, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%mass_fl_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mass_fl_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%mass_fl_e = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mass_fl_e, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%rho_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%rho_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%theta_v_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_v_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%theta_v_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_v_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%grf_tend_vn = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_tend_vn, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%grf_tend_vn = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_tend_vn, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%grf_tend_w = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_tend_w, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%grf_tend_w = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_tend_w, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%grf_tend_rho = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_tend_rho, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%grf_tend_rho = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_tend_rho, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%grf_tend_mflx = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_tend_mflx, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%grf_tend_mflx = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_tend_mflx, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%grf_bdy_mflx = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_bdy_mflx, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%grf_bdy_mflx = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_bdy_mflx, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%grf_tend_thv = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_tend_thv, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%grf_tend_thv = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%grf_tend_thv, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%vn_ie_int = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_ie_int, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vn_ie_int = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_ie_int, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%vn_ie_ubc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_ie_ubc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vn_ie_ubc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_ie_ubc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%w_int = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%w_int, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%w_int = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%w_int, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%w_ubc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%w_ubc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%w_ubc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%w_ubc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%theta_v_ic_int = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_v_ic_int, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%theta_v_ic_int = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_v_ic_int, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%theta_v_ic_ubc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_v_ic_ubc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%theta_v_ic_ubc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_v_ic_ubc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%rho_ic_int = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_ic_int, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%rho_ic_int = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_ic_int, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%rho_ic_ubc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_ic_ubc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%rho_ic_ubc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_ic_ubc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%mflx_ic_int = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mflx_ic_int, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%mflx_ic_int = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mflx_ic_int, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%mflx_ic_ubc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mflx_ic_ubc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%mflx_ic_ubc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mflx_ic_ubc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%vn_incr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_incr, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vn_incr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_incr, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%exner_incr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner_incr, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%exner_incr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner_incr, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%rho_incr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_incr, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%rho_incr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_incr, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
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
#ifndef _OPENACC
    dace_rich_obj%ddt_exner_phy = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_exner_phy, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_exner_phy = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_exner_phy, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_phy = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_phy, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_phy = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_phy, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%exner_dyn_incr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner_dyn_incr, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%exner_dyn_incr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner_dyn_incr, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
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
#ifndef _OPENACC
    dace_rich_obj%mass_fl_e_sv = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mass_fl_e_sv, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%mass_fl_e_sv = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mass_fl_e_sv, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
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
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_cor_pc = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%ddt_vn_cor_pc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_cor_pc = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%ddt_vn_cor_pc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
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
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_dyn = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_dyn, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_dyn = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_dyn, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_dmp = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_dmp, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_dmp = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_dmp, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_adv = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_adv, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_adv = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_adv, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_cor = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_cor, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_cor = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_cor, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_pgr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_pgr, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_pgr = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_pgr, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_phd = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_phd, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_phd = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_phd, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_iau = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_iau, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_iau = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_iau, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_ray = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_ray, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_ray = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_ray, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%ddt_vn_grf = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_grf, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%ddt_vn_grf = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%ddt_vn_grf, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif

  end function copy_in_t_nh_diag

  function copy_in_t_nh_ref(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_nh_ref), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_nh_ref), pointer :: dace_rich_obj
    type(dace_t_nh_ref) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_vn_ref_d_0_s = size(fortran_obj%vn_ref, dim=1)
    dace_rich_obj%f2dace_SOA_vn_ref_d_0_s = lbound(fortran_obj%vn_ref, dim=1)
    dace_rich_obj%f2dace_SA_vn_ref_d_1_s = size(fortran_obj%vn_ref, dim=2)
    dace_rich_obj%f2dace_SOA_vn_ref_d_1_s = lbound(fortran_obj%vn_ref, dim=2)
    dace_rich_obj%f2dace_SA_vn_ref_d_2_s = size(fortran_obj%vn_ref, dim=3)
    dace_rich_obj%f2dace_SOA_vn_ref_d_2_s = lbound(fortran_obj%vn_ref, dim=3)
    dace_rich_obj%f2dace_SA_w_ref_d_0_s = size(fortran_obj%w_ref, dim=1)
    dace_rich_obj%f2dace_SOA_w_ref_d_0_s = lbound(fortran_obj%w_ref, dim=1)
    dace_rich_obj%f2dace_SA_w_ref_d_1_s = size(fortran_obj%w_ref, dim=2)
    dace_rich_obj%f2dace_SOA_w_ref_d_1_s = lbound(fortran_obj%w_ref, dim=2)
    dace_rich_obj%f2dace_SA_w_ref_d_2_s = size(fortran_obj%w_ref, dim=3)
    dace_rich_obj%f2dace_SOA_w_ref_d_2_s = lbound(fortran_obj%w_ref, dim=3)
#ifndef _OPENACC
    dace_rich_obj%vn_ref = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_ref, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vn_ref = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_ref, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%w_ref = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%w_ref, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%w_ref = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%w_ref, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif

  end function copy_in_t_nh_ref

  function copy_in_t_nh_metrics(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_nh_metrics), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_nh_metrics), pointer :: dace_rich_obj
    type(dace_t_nh_metrics) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_rayleigh_w_d_0_s = size(fortran_obj%rayleigh_w, dim=1)
    dace_rich_obj%f2dace_SOA_rayleigh_w_d_0_s = lbound(fortran_obj%rayleigh_w, dim=1)
    dace_rich_obj%f2dace_SA_rayleigh_vn_d_0_s = size(fortran_obj%rayleigh_vn, dim=1)
    dace_rich_obj%f2dace_SOA_rayleigh_vn_d_0_s = lbound(fortran_obj%rayleigh_vn, dim=1)
    dace_rich_obj%f2dace_SA_scalfac_dd3d_d_0_s = size(fortran_obj%scalfac_dd3d, dim=1)
    dace_rich_obj%f2dace_SOA_scalfac_dd3d_d_0_s = lbound(fortran_obj%scalfac_dd3d, dim=1)
    dace_rich_obj%f2dace_SA_hmask_dd3d_d_0_s = size(fortran_obj%hmask_dd3d, dim=1)
    dace_rich_obj%f2dace_SOA_hmask_dd3d_d_0_s = lbound(fortran_obj%hmask_dd3d, dim=1)
    dace_rich_obj%f2dace_SA_hmask_dd3d_d_1_s = size(fortran_obj%hmask_dd3d, dim=2)
    dace_rich_obj%f2dace_SOA_hmask_dd3d_d_1_s = lbound(fortran_obj%hmask_dd3d, dim=2)
    dace_rich_obj%f2dace_SA_vwind_expl_wgt_d_0_s = size(fortran_obj%vwind_expl_wgt, dim=1)
    dace_rich_obj%f2dace_SOA_vwind_expl_wgt_d_0_s = lbound(fortran_obj%vwind_expl_wgt, dim=1)
    dace_rich_obj%f2dace_SA_vwind_expl_wgt_d_1_s = size(fortran_obj%vwind_expl_wgt, dim=2)
    dace_rich_obj%f2dace_SOA_vwind_expl_wgt_d_1_s = lbound(fortran_obj%vwind_expl_wgt, dim=2)
    dace_rich_obj%f2dace_SA_vwind_impl_wgt_d_0_s = size(fortran_obj%vwind_impl_wgt, dim=1)
    dace_rich_obj%f2dace_SOA_vwind_impl_wgt_d_0_s = lbound(fortran_obj%vwind_impl_wgt, dim=1)
    dace_rich_obj%f2dace_SA_vwind_impl_wgt_d_1_s = size(fortran_obj%vwind_impl_wgt, dim=2)
    dace_rich_obj%f2dace_SOA_vwind_impl_wgt_d_1_s = lbound(fortran_obj%vwind_impl_wgt, dim=2)
    dace_rich_obj%f2dace_SA_ddxn_z_full_d_0_s = size(fortran_obj%ddxn_z_full, dim=1)
    dace_rich_obj%f2dace_SOA_ddxn_z_full_d_0_s = lbound(fortran_obj%ddxn_z_full, dim=1)
    dace_rich_obj%f2dace_SA_ddxn_z_full_d_1_s = size(fortran_obj%ddxn_z_full, dim=2)
    dace_rich_obj%f2dace_SOA_ddxn_z_full_d_1_s = lbound(fortran_obj%ddxn_z_full, dim=2)
    dace_rich_obj%f2dace_SA_ddxn_z_full_d_2_s = size(fortran_obj%ddxn_z_full, dim=3)
    dace_rich_obj%f2dace_SOA_ddxn_z_full_d_2_s = lbound(fortran_obj%ddxn_z_full, dim=3)
    dace_rich_obj%f2dace_SA_ddxt_z_full_d_0_s = size(fortran_obj%ddxt_z_full, dim=1)
    dace_rich_obj%f2dace_SOA_ddxt_z_full_d_0_s = lbound(fortran_obj%ddxt_z_full, dim=1)
    dace_rich_obj%f2dace_SA_ddxt_z_full_d_1_s = size(fortran_obj%ddxt_z_full, dim=2)
    dace_rich_obj%f2dace_SOA_ddxt_z_full_d_1_s = lbound(fortran_obj%ddxt_z_full, dim=2)
    dace_rich_obj%f2dace_SA_ddxt_z_full_d_2_s = size(fortran_obj%ddxt_z_full, dim=3)
    dace_rich_obj%f2dace_SOA_ddxt_z_full_d_2_s = lbound(fortran_obj%ddxt_z_full, dim=3)
    dace_rich_obj%f2dace_SA_ddqz_z_full_e_d_0_s = size(fortran_obj%ddqz_z_full_e, dim=1)
    dace_rich_obj%f2dace_SOA_ddqz_z_full_e_d_0_s = lbound(fortran_obj%ddqz_z_full_e, dim=1)
    dace_rich_obj%f2dace_SA_ddqz_z_full_e_d_1_s = size(fortran_obj%ddqz_z_full_e, dim=2)
    dace_rich_obj%f2dace_SOA_ddqz_z_full_e_d_1_s = lbound(fortran_obj%ddqz_z_full_e, dim=2)
    dace_rich_obj%f2dace_SA_ddqz_z_full_e_d_2_s = size(fortran_obj%ddqz_z_full_e, dim=3)
    dace_rich_obj%f2dace_SOA_ddqz_z_full_e_d_2_s = lbound(fortran_obj%ddqz_z_full_e, dim=3)
    dace_rich_obj%f2dace_SA_ddqz_z_half_d_0_s = size(fortran_obj%ddqz_z_half, dim=1)
    dace_rich_obj%f2dace_SOA_ddqz_z_half_d_0_s = lbound(fortran_obj%ddqz_z_half, dim=1)
    dace_rich_obj%f2dace_SA_ddqz_z_half_d_1_s = size(fortran_obj%ddqz_z_half, dim=2)
    dace_rich_obj%f2dace_SOA_ddqz_z_half_d_1_s = lbound(fortran_obj%ddqz_z_half, dim=2)
    dace_rich_obj%f2dace_SA_ddqz_z_half_d_2_s = size(fortran_obj%ddqz_z_half, dim=3)
    dace_rich_obj%f2dace_SOA_ddqz_z_half_d_2_s = lbound(fortran_obj%ddqz_z_half, dim=3)
    dace_rich_obj%f2dace_SA_inv_ddqz_z_full_d_0_s = size(fortran_obj%inv_ddqz_z_full, dim=1)
    dace_rich_obj%f2dace_SOA_inv_ddqz_z_full_d_0_s = lbound(fortran_obj%inv_ddqz_z_full, dim=1)
    dace_rich_obj%f2dace_SA_inv_ddqz_z_full_d_1_s = size(fortran_obj%inv_ddqz_z_full, dim=2)
    dace_rich_obj%f2dace_SOA_inv_ddqz_z_full_d_1_s = lbound(fortran_obj%inv_ddqz_z_full, dim=2)
    dace_rich_obj%f2dace_SA_inv_ddqz_z_full_d_2_s = size(fortran_obj%inv_ddqz_z_full, dim=3)
    dace_rich_obj%f2dace_SOA_inv_ddqz_z_full_d_2_s = lbound(fortran_obj%inv_ddqz_z_full, dim=3)
    dace_rich_obj%f2dace_SA_wgtfac_c_d_0_s = size(fortran_obj%wgtfac_c, dim=1)
    dace_rich_obj%f2dace_SOA_wgtfac_c_d_0_s = lbound(fortran_obj%wgtfac_c, dim=1)
    dace_rich_obj%f2dace_SA_wgtfac_c_d_1_s = size(fortran_obj%wgtfac_c, dim=2)
    dace_rich_obj%f2dace_SOA_wgtfac_c_d_1_s = lbound(fortran_obj%wgtfac_c, dim=2)
    dace_rich_obj%f2dace_SA_wgtfac_c_d_2_s = size(fortran_obj%wgtfac_c, dim=3)
    dace_rich_obj%f2dace_SOA_wgtfac_c_d_2_s = lbound(fortran_obj%wgtfac_c, dim=3)
    dace_rich_obj%f2dace_SA_wgtfac_e_d_0_s = size(fortran_obj%wgtfac_e, dim=1)
    dace_rich_obj%f2dace_SOA_wgtfac_e_d_0_s = lbound(fortran_obj%wgtfac_e, dim=1)
    dace_rich_obj%f2dace_SA_wgtfac_e_d_1_s = size(fortran_obj%wgtfac_e, dim=2)
    dace_rich_obj%f2dace_SOA_wgtfac_e_d_1_s = lbound(fortran_obj%wgtfac_e, dim=2)
    dace_rich_obj%f2dace_SA_wgtfac_e_d_2_s = size(fortran_obj%wgtfac_e, dim=3)
    dace_rich_obj%f2dace_SOA_wgtfac_e_d_2_s = lbound(fortran_obj%wgtfac_e, dim=3)
    dace_rich_obj%f2dace_SA_wgtfacq_c_d_0_s = size(fortran_obj%wgtfacq_c, dim=1)
    dace_rich_obj%f2dace_SOA_wgtfacq_c_d_0_s = lbound(fortran_obj%wgtfacq_c, dim=1)
    dace_rich_obj%f2dace_SA_wgtfacq_c_d_1_s = size(fortran_obj%wgtfacq_c, dim=2)
    dace_rich_obj%f2dace_SOA_wgtfacq_c_d_1_s = lbound(fortran_obj%wgtfacq_c, dim=2)
    dace_rich_obj%f2dace_SA_wgtfacq_c_d_2_s = size(fortran_obj%wgtfacq_c, dim=3)
    dace_rich_obj%f2dace_SOA_wgtfacq_c_d_2_s = lbound(fortran_obj%wgtfacq_c, dim=3)
    dace_rich_obj%f2dace_SA_wgtfacq_e_d_0_s = size(fortran_obj%wgtfacq_e, dim=1)
    dace_rich_obj%f2dace_SOA_wgtfacq_e_d_0_s = lbound(fortran_obj%wgtfacq_e, dim=1)
    dace_rich_obj%f2dace_SA_wgtfacq_e_d_1_s = size(fortran_obj%wgtfacq_e, dim=2)
    dace_rich_obj%f2dace_SOA_wgtfacq_e_d_1_s = lbound(fortran_obj%wgtfacq_e, dim=2)
    dace_rich_obj%f2dace_SA_wgtfacq_e_d_2_s = size(fortran_obj%wgtfacq_e, dim=3)
    dace_rich_obj%f2dace_SOA_wgtfacq_e_d_2_s = lbound(fortran_obj%wgtfacq_e, dim=3)
    dace_rich_obj%f2dace_SA_wgtfacq1_c_d_0_s = size(fortran_obj%wgtfacq1_c, dim=1)
    dace_rich_obj%f2dace_SOA_wgtfacq1_c_d_0_s = lbound(fortran_obj%wgtfacq1_c, dim=1)
    dace_rich_obj%f2dace_SA_wgtfacq1_c_d_1_s = size(fortran_obj%wgtfacq1_c, dim=2)
    dace_rich_obj%f2dace_SOA_wgtfacq1_c_d_1_s = lbound(fortran_obj%wgtfacq1_c, dim=2)
    dace_rich_obj%f2dace_SA_wgtfacq1_c_d_2_s = size(fortran_obj%wgtfacq1_c, dim=3)
    dace_rich_obj%f2dace_SOA_wgtfacq1_c_d_2_s = lbound(fortran_obj%wgtfacq1_c, dim=3)
    dace_rich_obj%f2dace_SA_coeff_gradekin_d_0_s = size(fortran_obj%coeff_gradekin, dim=1)
    dace_rich_obj%f2dace_SOA_coeff_gradekin_d_0_s = lbound(fortran_obj%coeff_gradekin, dim=1)
    dace_rich_obj%f2dace_SA_coeff_gradekin_d_1_s = size(fortran_obj%coeff_gradekin, dim=2)
    dace_rich_obj%f2dace_SOA_coeff_gradekin_d_1_s = lbound(fortran_obj%coeff_gradekin, dim=2)
    dace_rich_obj%f2dace_SA_coeff_gradekin_d_2_s = size(fortran_obj%coeff_gradekin, dim=3)
    dace_rich_obj%f2dace_SOA_coeff_gradekin_d_2_s = lbound(fortran_obj%coeff_gradekin, dim=3)
    dace_rich_obj%f2dace_SA_coeff1_dwdz_d_0_s = size(fortran_obj%coeff1_dwdz, dim=1)
    dace_rich_obj%f2dace_SOA_coeff1_dwdz_d_0_s = lbound(fortran_obj%coeff1_dwdz, dim=1)
    dace_rich_obj%f2dace_SA_coeff1_dwdz_d_1_s = size(fortran_obj%coeff1_dwdz, dim=2)
    dace_rich_obj%f2dace_SOA_coeff1_dwdz_d_1_s = lbound(fortran_obj%coeff1_dwdz, dim=2)
    dace_rich_obj%f2dace_SA_coeff1_dwdz_d_2_s = size(fortran_obj%coeff1_dwdz, dim=3)
    dace_rich_obj%f2dace_SOA_coeff1_dwdz_d_2_s = lbound(fortran_obj%coeff1_dwdz, dim=3)
    dace_rich_obj%f2dace_SA_coeff2_dwdz_d_0_s = size(fortran_obj%coeff2_dwdz, dim=1)
    dace_rich_obj%f2dace_SOA_coeff2_dwdz_d_0_s = lbound(fortran_obj%coeff2_dwdz, dim=1)
    dace_rich_obj%f2dace_SA_coeff2_dwdz_d_1_s = size(fortran_obj%coeff2_dwdz, dim=2)
    dace_rich_obj%f2dace_SOA_coeff2_dwdz_d_1_s = lbound(fortran_obj%coeff2_dwdz, dim=2)
    dace_rich_obj%f2dace_SA_coeff2_dwdz_d_2_s = size(fortran_obj%coeff2_dwdz, dim=3)
    dace_rich_obj%f2dace_SOA_coeff2_dwdz_d_2_s = lbound(fortran_obj%coeff2_dwdz, dim=3)
    dace_rich_obj%f2dace_SA_zdiff_gradp_d_0_s = size(fortran_obj%zdiff_gradp, dim=1)
    dace_rich_obj%f2dace_SOA_zdiff_gradp_d_0_s = lbound(fortran_obj%zdiff_gradp, dim=1)
    dace_rich_obj%f2dace_SA_zdiff_gradp_d_1_s = size(fortran_obj%zdiff_gradp, dim=2)
    dace_rich_obj%f2dace_SOA_zdiff_gradp_d_1_s = lbound(fortran_obj%zdiff_gradp, dim=2)
    dace_rich_obj%f2dace_SA_zdiff_gradp_d_2_s = size(fortran_obj%zdiff_gradp, dim=3)
    dace_rich_obj%f2dace_SOA_zdiff_gradp_d_2_s = lbound(fortran_obj%zdiff_gradp, dim=3)
    dace_rich_obj%f2dace_SA_zdiff_gradp_d_3_s = size(fortran_obj%zdiff_gradp, dim=4)
    dace_rich_obj%f2dace_SOA_zdiff_gradp_d_3_s = lbound(fortran_obj%zdiff_gradp, dim=4)
    dace_rich_obj%f2dace_SA_coeff_gradp_d_0_s = size(fortran_obj%coeff_gradp, dim=1)
    dace_rich_obj%f2dace_SOA_coeff_gradp_d_0_s = lbound(fortran_obj%coeff_gradp, dim=1)
    dace_rich_obj%f2dace_SA_coeff_gradp_d_1_s = size(fortran_obj%coeff_gradp, dim=2)
    dace_rich_obj%f2dace_SOA_coeff_gradp_d_1_s = lbound(fortran_obj%coeff_gradp, dim=2)
    dace_rich_obj%f2dace_SA_coeff_gradp_d_2_s = size(fortran_obj%coeff_gradp, dim=3)
    dace_rich_obj%f2dace_SOA_coeff_gradp_d_2_s = lbound(fortran_obj%coeff_gradp, dim=3)
    dace_rich_obj%f2dace_SA_coeff_gradp_d_3_s = size(fortran_obj%coeff_gradp, dim=4)
    dace_rich_obj%f2dace_SOA_coeff_gradp_d_3_s = lbound(fortran_obj%coeff_gradp, dim=4)
    dace_rich_obj%f2dace_SA_exner_exfac_d_0_s = size(fortran_obj%exner_exfac, dim=1)
    dace_rich_obj%f2dace_SOA_exner_exfac_d_0_s = lbound(fortran_obj%exner_exfac, dim=1)
    dace_rich_obj%f2dace_SA_exner_exfac_d_1_s = size(fortran_obj%exner_exfac, dim=2)
    dace_rich_obj%f2dace_SOA_exner_exfac_d_1_s = lbound(fortran_obj%exner_exfac, dim=2)
    dace_rich_obj%f2dace_SA_exner_exfac_d_2_s = size(fortran_obj%exner_exfac, dim=3)
    dace_rich_obj%f2dace_SOA_exner_exfac_d_2_s = lbound(fortran_obj%exner_exfac, dim=3)
    dace_rich_obj%f2dace_SA_theta_ref_mc_d_0_s = size(fortran_obj%theta_ref_mc, dim=1)
    dace_rich_obj%f2dace_SOA_theta_ref_mc_d_0_s = lbound(fortran_obj%theta_ref_mc, dim=1)
    dace_rich_obj%f2dace_SA_theta_ref_mc_d_1_s = size(fortran_obj%theta_ref_mc, dim=2)
    dace_rich_obj%f2dace_SOA_theta_ref_mc_d_1_s = lbound(fortran_obj%theta_ref_mc, dim=2)
    dace_rich_obj%f2dace_SA_theta_ref_mc_d_2_s = size(fortran_obj%theta_ref_mc, dim=3)
    dace_rich_obj%f2dace_SOA_theta_ref_mc_d_2_s = lbound(fortran_obj%theta_ref_mc, dim=3)
    dace_rich_obj%f2dace_SA_theta_ref_me_d_0_s = size(fortran_obj%theta_ref_me, dim=1)
    dace_rich_obj%f2dace_SOA_theta_ref_me_d_0_s = lbound(fortran_obj%theta_ref_me, dim=1)
    dace_rich_obj%f2dace_SA_theta_ref_me_d_1_s = size(fortran_obj%theta_ref_me, dim=2)
    dace_rich_obj%f2dace_SOA_theta_ref_me_d_1_s = lbound(fortran_obj%theta_ref_me, dim=2)
    dace_rich_obj%f2dace_SA_theta_ref_me_d_2_s = size(fortran_obj%theta_ref_me, dim=3)
    dace_rich_obj%f2dace_SOA_theta_ref_me_d_2_s = lbound(fortran_obj%theta_ref_me, dim=3)
    dace_rich_obj%f2dace_SA_theta_ref_ic_d_0_s = size(fortran_obj%theta_ref_ic, dim=1)
    dace_rich_obj%f2dace_SOA_theta_ref_ic_d_0_s = lbound(fortran_obj%theta_ref_ic, dim=1)
    dace_rich_obj%f2dace_SA_theta_ref_ic_d_1_s = size(fortran_obj%theta_ref_ic, dim=2)
    dace_rich_obj%f2dace_SOA_theta_ref_ic_d_1_s = lbound(fortran_obj%theta_ref_ic, dim=2)
    dace_rich_obj%f2dace_SA_theta_ref_ic_d_2_s = size(fortran_obj%theta_ref_ic, dim=3)
    dace_rich_obj%f2dace_SOA_theta_ref_ic_d_2_s = lbound(fortran_obj%theta_ref_ic, dim=3)
    dace_rich_obj%f2dace_SA_exner_ref_mc_d_0_s = size(fortran_obj%exner_ref_mc, dim=1)
    dace_rich_obj%f2dace_SOA_exner_ref_mc_d_0_s = lbound(fortran_obj%exner_ref_mc, dim=1)
    dace_rich_obj%f2dace_SA_exner_ref_mc_d_1_s = size(fortran_obj%exner_ref_mc, dim=2)
    dace_rich_obj%f2dace_SOA_exner_ref_mc_d_1_s = lbound(fortran_obj%exner_ref_mc, dim=2)
    dace_rich_obj%f2dace_SA_exner_ref_mc_d_2_s = size(fortran_obj%exner_ref_mc, dim=3)
    dace_rich_obj%f2dace_SOA_exner_ref_mc_d_2_s = lbound(fortran_obj%exner_ref_mc, dim=3)
    dace_rich_obj%f2dace_SA_rho_ref_mc_d_0_s = size(fortran_obj%rho_ref_mc, dim=1)
    dace_rich_obj%f2dace_SOA_rho_ref_mc_d_0_s = lbound(fortran_obj%rho_ref_mc, dim=1)
    dace_rich_obj%f2dace_SA_rho_ref_mc_d_1_s = size(fortran_obj%rho_ref_mc, dim=2)
    dace_rich_obj%f2dace_SOA_rho_ref_mc_d_1_s = lbound(fortran_obj%rho_ref_mc, dim=2)
    dace_rich_obj%f2dace_SA_rho_ref_mc_d_2_s = size(fortran_obj%rho_ref_mc, dim=3)
    dace_rich_obj%f2dace_SOA_rho_ref_mc_d_2_s = lbound(fortran_obj%rho_ref_mc, dim=3)
    dace_rich_obj%f2dace_SA_rho_ref_me_d_0_s = size(fortran_obj%rho_ref_me, dim=1)
    dace_rich_obj%f2dace_SOA_rho_ref_me_d_0_s = lbound(fortran_obj%rho_ref_me, dim=1)
    dace_rich_obj%f2dace_SA_rho_ref_me_d_1_s = size(fortran_obj%rho_ref_me, dim=2)
    dace_rich_obj%f2dace_SOA_rho_ref_me_d_1_s = lbound(fortran_obj%rho_ref_me, dim=2)
    dace_rich_obj%f2dace_SA_rho_ref_me_d_2_s = size(fortran_obj%rho_ref_me, dim=3)
    dace_rich_obj%f2dace_SOA_rho_ref_me_d_2_s = lbound(fortran_obj%rho_ref_me, dim=3)
    dace_rich_obj%f2dace_SA_d_exner_dz_ref_ic_d_0_s = size(fortran_obj%d_exner_dz_ref_ic, dim=1)
    dace_rich_obj%f2dace_SOA_d_exner_dz_ref_ic_d_0_s = lbound(fortran_obj%d_exner_dz_ref_ic, dim=1)
    dace_rich_obj%f2dace_SA_d_exner_dz_ref_ic_d_1_s = size(fortran_obj%d_exner_dz_ref_ic, dim=2)
    dace_rich_obj%f2dace_SOA_d_exner_dz_ref_ic_d_1_s = lbound(fortran_obj%d_exner_dz_ref_ic, dim=2)
    dace_rich_obj%f2dace_SA_d_exner_dz_ref_ic_d_2_s = size(fortran_obj%d_exner_dz_ref_ic, dim=3)
    dace_rich_obj%f2dace_SOA_d_exner_dz_ref_ic_d_2_s = lbound(fortran_obj%d_exner_dz_ref_ic, dim=3)
    dace_rich_obj%f2dace_SA_d2dexdz2_fac1_mc_d_0_s = size(fortran_obj%d2dexdz2_fac1_mc, dim=1)
    dace_rich_obj%f2dace_SOA_d2dexdz2_fac1_mc_d_0_s = lbound(fortran_obj%d2dexdz2_fac1_mc, dim=1)
    dace_rich_obj%f2dace_SA_d2dexdz2_fac1_mc_d_1_s = size(fortran_obj%d2dexdz2_fac1_mc, dim=2)
    dace_rich_obj%f2dace_SOA_d2dexdz2_fac1_mc_d_1_s = lbound(fortran_obj%d2dexdz2_fac1_mc, dim=2)
    dace_rich_obj%f2dace_SA_d2dexdz2_fac1_mc_d_2_s = size(fortran_obj%d2dexdz2_fac1_mc, dim=3)
    dace_rich_obj%f2dace_SOA_d2dexdz2_fac1_mc_d_2_s = lbound(fortran_obj%d2dexdz2_fac1_mc, dim=3)
    dace_rich_obj%f2dace_SA_d2dexdz2_fac2_mc_d_0_s = size(fortran_obj%d2dexdz2_fac2_mc, dim=1)
    dace_rich_obj%f2dace_SOA_d2dexdz2_fac2_mc_d_0_s = lbound(fortran_obj%d2dexdz2_fac2_mc, dim=1)
    dace_rich_obj%f2dace_SA_d2dexdz2_fac2_mc_d_1_s = size(fortran_obj%d2dexdz2_fac2_mc, dim=2)
    dace_rich_obj%f2dace_SOA_d2dexdz2_fac2_mc_d_1_s = lbound(fortran_obj%d2dexdz2_fac2_mc, dim=2)
    dace_rich_obj%f2dace_SA_d2dexdz2_fac2_mc_d_2_s = size(fortran_obj%d2dexdz2_fac2_mc, dim=3)
    dace_rich_obj%f2dace_SOA_d2dexdz2_fac2_mc_d_2_s = lbound(fortran_obj%d2dexdz2_fac2_mc, dim=3)
    dace_rich_obj%f2dace_SA_pg_exdist_d_0_s = size(fortran_obj%pg_exdist, dim=1)
    dace_rich_obj%f2dace_SOA_pg_exdist_d_0_s = lbound(fortran_obj%pg_exdist, dim=1)
    dace_rich_obj%f2dace_SA_vertidx_gradp_d_0_s = size(fortran_obj%vertidx_gradp, dim=1)
    dace_rich_obj%f2dace_SOA_vertidx_gradp_d_0_s = lbound(fortran_obj%vertidx_gradp, dim=1)
    dace_rich_obj%f2dace_SA_vertidx_gradp_d_1_s = size(fortran_obj%vertidx_gradp, dim=2)
    dace_rich_obj%f2dace_SOA_vertidx_gradp_d_1_s = lbound(fortran_obj%vertidx_gradp, dim=2)
    dace_rich_obj%f2dace_SA_vertidx_gradp_d_2_s = size(fortran_obj%vertidx_gradp, dim=3)
    dace_rich_obj%f2dace_SOA_vertidx_gradp_d_2_s = lbound(fortran_obj%vertidx_gradp, dim=3)
    dace_rich_obj%f2dace_SA_vertidx_gradp_d_3_s = size(fortran_obj%vertidx_gradp, dim=4)
    dace_rich_obj%f2dace_SOA_vertidx_gradp_d_3_s = lbound(fortran_obj%vertidx_gradp, dim=4)
    dace_rich_obj%f2dace_SA_pg_edgeidx_d_0_s = size(fortran_obj%pg_edgeidx, dim=1)
    dace_rich_obj%f2dace_SOA_pg_edgeidx_d_0_s = lbound(fortran_obj%pg_edgeidx, dim=1)
    dace_rich_obj%f2dace_SA_pg_edgeblk_d_0_s = size(fortran_obj%pg_edgeblk, dim=1)
    dace_rich_obj%f2dace_SOA_pg_edgeblk_d_0_s = lbound(fortran_obj%pg_edgeblk, dim=1)
    dace_rich_obj%f2dace_SA_pg_vertidx_d_0_s = size(fortran_obj%pg_vertidx, dim=1)
    dace_rich_obj%f2dace_SOA_pg_vertidx_d_0_s = lbound(fortran_obj%pg_vertidx, dim=1)
    dace_rich_obj%f2dace_SA_bdy_mflx_e_idx_d_0_s = size(fortran_obj%bdy_mflx_e_idx, dim=1)
    dace_rich_obj%f2dace_SOA_bdy_mflx_e_idx_d_0_s = lbound(fortran_obj%bdy_mflx_e_idx, dim=1)
    dace_rich_obj%f2dace_SA_bdy_mflx_e_blk_d_0_s = size(fortran_obj%bdy_mflx_e_blk, dim=1)
    dace_rich_obj%f2dace_SOA_bdy_mflx_e_blk_d_0_s = lbound(fortran_obj%bdy_mflx_e_blk, dim=1)
    dace_rich_obj%f2dace_SA_deepatmo_gradh_mc_d_0_s = size(fortran_obj%deepatmo_gradh_mc, dim=1)
    dace_rich_obj%f2dace_SOA_deepatmo_gradh_mc_d_0_s = lbound(fortran_obj%deepatmo_gradh_mc, dim=1)
    dace_rich_obj%f2dace_SA_deepatmo_divh_mc_d_0_s = size(fortran_obj%deepatmo_divh_mc, dim=1)
    dace_rich_obj%f2dace_SOA_deepatmo_divh_mc_d_0_s = lbound(fortran_obj%deepatmo_divh_mc, dim=1)
    dace_rich_obj%f2dace_SA_deepatmo_invr_mc_d_0_s = size(fortran_obj%deepatmo_invr_mc, dim=1)
    dace_rich_obj%f2dace_SOA_deepatmo_invr_mc_d_0_s = lbound(fortran_obj%deepatmo_invr_mc, dim=1)
    dace_rich_obj%f2dace_SA_deepatmo_divzu_mc_d_0_s = size(fortran_obj%deepatmo_divzu_mc, dim=1)
    dace_rich_obj%f2dace_SOA_deepatmo_divzu_mc_d_0_s = lbound(fortran_obj%deepatmo_divzu_mc, dim=1)
    dace_rich_obj%f2dace_SA_deepatmo_divzl_mc_d_0_s = size(fortran_obj%deepatmo_divzl_mc, dim=1)
    dace_rich_obj%f2dace_SOA_deepatmo_divzl_mc_d_0_s = lbound(fortran_obj%deepatmo_divzl_mc, dim=1)
    dace_rich_obj%f2dace_SA_deepatmo_gradh_ifc_d_0_s = size(fortran_obj%deepatmo_gradh_ifc, dim=1)
    dace_rich_obj%f2dace_SOA_deepatmo_gradh_ifc_d_0_s = lbound(fortran_obj%deepatmo_gradh_ifc, dim=1)
    dace_rich_obj%f2dace_SA_deepatmo_invr_ifc_d_0_s = size(fortran_obj%deepatmo_invr_ifc, dim=1)
    dace_rich_obj%f2dace_SOA_deepatmo_invr_ifc_d_0_s = lbound(fortran_obj%deepatmo_invr_ifc, dim=1)
    dace_rich_obj%pg_listdim = fortran_obj%pg_listdim
    dace_rich_obj%bdy_mflx_e_dim = fortran_obj%bdy_mflx_e_dim
#ifndef _OPENACC
    dace_rich_obj%rayleigh_w = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%rayleigh_w, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%rayleigh_w = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%rayleigh_w, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%rayleigh_vn = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%rayleigh_vn, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%rayleigh_vn = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%rayleigh_vn, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%scalfac_dd3d = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%scalfac_dd3d, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%scalfac_dd3d = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%scalfac_dd3d, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%hmask_dd3d = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%hmask_dd3d, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%hmask_dd3d = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%hmask_dd3d, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%vwind_expl_wgt = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%vwind_expl_wgt, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vwind_expl_wgt = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%vwind_expl_wgt, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%vwind_impl_wgt = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%vwind_impl_wgt, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vwind_impl_wgt = copy_in_float64_2d_array( &
    fortran_array=fortran_obj%vwind_impl_wgt, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
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
#ifndef _OPENACC
    dace_rich_obj%inv_ddqz_z_full = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%inv_ddqz_z_full, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%inv_ddqz_z_full = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%inv_ddqz_z_full, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
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
#ifndef _OPENACC
    dace_rich_obj%wgtfacq_c = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%wgtfacq_c, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%wgtfacq_c = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%wgtfacq_c, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
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
    dace_rich_obj%wgtfacq1_c = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%wgtfacq1_c, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%wgtfacq1_c = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%wgtfacq1_c, &
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
#ifndef _OPENACC
    dace_rich_obj%zdiff_gradp = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%zdiff_gradp, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%zdiff_gradp = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%zdiff_gradp, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%coeff_gradp = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%coeff_gradp, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%coeff_gradp = copy_in_float64_4d_array( &
    fortran_array=fortran_obj%coeff_gradp, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%exner_exfac = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner_exfac, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%exner_exfac = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner_exfac, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%theta_ref_mc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_ref_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%theta_ref_mc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_ref_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%theta_ref_me = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_ref_me, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%theta_ref_me = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_ref_me, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%theta_ref_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_ref_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%theta_ref_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_ref_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%exner_ref_mc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner_ref_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%exner_ref_mc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner_ref_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%rho_ref_mc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_ref_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%rho_ref_mc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_ref_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%rho_ref_me = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_ref_me, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%rho_ref_me = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho_ref_me, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%d_exner_dz_ref_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%d_exner_dz_ref_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%d_exner_dz_ref_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%d_exner_dz_ref_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%d2dexdz2_fac1_mc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%d2dexdz2_fac1_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%d2dexdz2_fac1_mc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%d2dexdz2_fac1_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%d2dexdz2_fac2_mc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%d2dexdz2_fac2_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%d2dexdz2_fac2_mc = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%d2dexdz2_fac2_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
    dace_rich_obj%pg_exdist = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%pg_exdist, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
#ifndef _OPENACC
    dace_rich_obj%vertidx_gradp = copy_in_int32_4d_array( &
    fortran_array=logical_fix_4d(fortran_obj%vertidx_gradp), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vertidx_gradp = copy_in_int32_4d_array( &
    fortran_array=logical_fix_4d(fortran_obj%vertidx_gradp), &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
    dace_rich_obj%pg_edgeidx = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%pg_edgeidx), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%pg_edgeblk = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%pg_edgeblk), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%pg_vertidx = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%pg_vertidx), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%bdy_mflx_e_idx = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%bdy_mflx_e_idx), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%bdy_mflx_e_blk = copy_in_int32_1d_array( &
    fortran_array=logical_fix_1d(fortran_obj%bdy_mflx_e_blk), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
#ifndef _OPENACC
    dace_rich_obj%deepatmo_gradh_mc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_gradh_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%deepatmo_gradh_mc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_gradh_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%deepatmo_divh_mc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_divh_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%deepatmo_divh_mc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_divh_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
    dace_rich_obj%deepatmo_invr_mc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_invr_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
#ifndef _OPENACC
    dace_rich_obj%deepatmo_divzu_mc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_divzu_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%deepatmo_divzu_mc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_divzu_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%deepatmo_divzl_mc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_divzl_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%deepatmo_divzl_mc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_divzl_mc, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
    dace_rich_obj%deepatmo_gradh_ifc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_gradh_ifc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%deepatmo_invr_ifc = copy_in_float64_1d_array( &
    fortran_array=fortran_obj%deepatmo_invr_ifc, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

  end function copy_in_t_nh_metrics

  function copy_in_t_nh_prog(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_nh_prog), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_nh_prog), pointer :: dace_rich_obj
    type(dace_t_nh_prog) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_w_d_0_s = size(fortran_obj%w, dim=1)
    dace_rich_obj%f2dace_SOA_w_d_0_s = lbound(fortran_obj%w, dim=1)
    dace_rich_obj%f2dace_SA_w_d_1_s = size(fortran_obj%w, dim=2)
    dace_rich_obj%f2dace_SOA_w_d_1_s = lbound(fortran_obj%w, dim=2)
    dace_rich_obj%f2dace_SA_w_d_2_s = size(fortran_obj%w, dim=3)
    dace_rich_obj%f2dace_SOA_w_d_2_s = lbound(fortran_obj%w, dim=3)
    dace_rich_obj%f2dace_SA_vn_d_0_s = size(fortran_obj%vn, dim=1)
    dace_rich_obj%f2dace_SOA_vn_d_0_s = lbound(fortran_obj%vn, dim=1)
    dace_rich_obj%f2dace_SA_vn_d_1_s = size(fortran_obj%vn, dim=2)
    dace_rich_obj%f2dace_SOA_vn_d_1_s = lbound(fortran_obj%vn, dim=2)
    dace_rich_obj%f2dace_SA_vn_d_2_s = size(fortran_obj%vn, dim=3)
    dace_rich_obj%f2dace_SOA_vn_d_2_s = lbound(fortran_obj%vn, dim=3)
    dace_rich_obj%f2dace_SA_rho_d_0_s = size(fortran_obj%rho, dim=1)
    dace_rich_obj%f2dace_SOA_rho_d_0_s = lbound(fortran_obj%rho, dim=1)
    dace_rich_obj%f2dace_SA_rho_d_1_s = size(fortran_obj%rho, dim=2)
    dace_rich_obj%f2dace_SOA_rho_d_1_s = lbound(fortran_obj%rho, dim=2)
    dace_rich_obj%f2dace_SA_rho_d_2_s = size(fortran_obj%rho, dim=3)
    dace_rich_obj%f2dace_SOA_rho_d_2_s = lbound(fortran_obj%rho, dim=3)
    dace_rich_obj%f2dace_SA_exner_d_0_s = size(fortran_obj%exner, dim=1)
    dace_rich_obj%f2dace_SOA_exner_d_0_s = lbound(fortran_obj%exner, dim=1)
    dace_rich_obj%f2dace_SA_exner_d_1_s = size(fortran_obj%exner, dim=2)
    dace_rich_obj%f2dace_SOA_exner_d_1_s = lbound(fortran_obj%exner, dim=2)
    dace_rich_obj%f2dace_SA_exner_d_2_s = size(fortran_obj%exner, dim=3)
    dace_rich_obj%f2dace_SOA_exner_d_2_s = lbound(fortran_obj%exner, dim=3)
    dace_rich_obj%f2dace_SA_theta_v_d_0_s = size(fortran_obj%theta_v, dim=1)
    dace_rich_obj%f2dace_SOA_theta_v_d_0_s = lbound(fortran_obj%theta_v, dim=1)
    dace_rich_obj%f2dace_SA_theta_v_d_1_s = size(fortran_obj%theta_v, dim=2)
    dace_rich_obj%f2dace_SOA_theta_v_d_1_s = lbound(fortran_obj%theta_v, dim=2)
    dace_rich_obj%f2dace_SA_theta_v_d_2_s = size(fortran_obj%theta_v, dim=3)
    dace_rich_obj%f2dace_SOA_theta_v_d_2_s = lbound(fortran_obj%theta_v, dim=3)
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
#ifndef _OPENACC
    dace_rich_obj%rho = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%rho = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%rho, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%exner = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%exner = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%exner, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%theta_v = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_v, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%theta_v = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%theta_v, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif

  end function copy_in_t_nh_prog

  function copy_in_t_patch(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_patch), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_patch), pointer :: dace_rich_obj
    type(dace_t_patch) :: dace_c_obj

    if (properly_cached_t_patch /= c_null_ptr) then
      return properly_cached_t_patch
    end if

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    properly_cached_t_patch = dace_obj_ptr
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%id = fortran_obj%id
    dace_rich_obj%n_childdom = fortran_obj%n_childdom
    dace_rich_obj%nblks_c = fortran_obj%nblks_c
    dace_rich_obj%nblks_e = fortran_obj%nblks_e
    dace_rich_obj%nblks_v = fortran_obj%nblks_v
    dace_rich_obj%nlev = fortran_obj%nlev
    dace_rich_obj%nlevp1 = fortran_obj%nlevp1
    dace_rich_obj%nshift = fortran_obj%nshift
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

    dace_rich_obj%f2dace_SA_neighbor_idx_d_0_s = size(fortran_obj%neighbor_idx, dim=1)
    dace_rich_obj%f2dace_SOA_neighbor_idx_d_0_s = lbound(fortran_obj%neighbor_idx, dim=1)
    dace_rich_obj%f2dace_SA_neighbor_idx_d_1_s = size(fortran_obj%neighbor_idx, dim=2)
    dace_rich_obj%f2dace_SOA_neighbor_idx_d_1_s = lbound(fortran_obj%neighbor_idx, dim=2)
    dace_rich_obj%f2dace_SA_neighbor_idx_d_2_s = size(fortran_obj%neighbor_idx, dim=3)
    dace_rich_obj%f2dace_SOA_neighbor_idx_d_2_s = lbound(fortran_obj%neighbor_idx, dim=3)
    dace_rich_obj%f2dace_SA_neighbor_blk_d_0_s = size(fortran_obj%neighbor_blk, dim=1)
    dace_rich_obj%f2dace_SOA_neighbor_blk_d_0_s = lbound(fortran_obj%neighbor_blk, dim=1)
    dace_rich_obj%f2dace_SA_neighbor_blk_d_1_s = size(fortran_obj%neighbor_blk, dim=2)
    dace_rich_obj%f2dace_SOA_neighbor_blk_d_1_s = lbound(fortran_obj%neighbor_blk, dim=2)
    dace_rich_obj%f2dace_SA_neighbor_blk_d_2_s = size(fortran_obj%neighbor_blk, dim=3)
    dace_rich_obj%f2dace_SOA_neighbor_blk_d_2_s = lbound(fortran_obj%neighbor_blk, dim=3)
    dace_rich_obj%f2dace_SA_edge_idx_d_0_s = size(fortran_obj%edge_idx, dim=1)
    dace_rich_obj%f2dace_SOA_edge_idx_d_0_s = lbound(fortran_obj%edge_idx, dim=1)
    dace_rich_obj%f2dace_SA_edge_idx_d_1_s = size(fortran_obj%edge_idx, dim=2)
    dace_rich_obj%f2dace_SOA_edge_idx_d_1_s = lbound(fortran_obj%edge_idx, dim=2)
    dace_rich_obj%f2dace_SA_edge_idx_d_2_s = size(fortran_obj%edge_idx, dim=3)
    dace_rich_obj%f2dace_SOA_edge_idx_d_2_s = lbound(fortran_obj%edge_idx, dim=3)
    dace_rich_obj%f2dace_SA_edge_blk_d_0_s = size(fortran_obj%edge_blk, dim=1)
    dace_rich_obj%f2dace_SOA_edge_blk_d_0_s = lbound(fortran_obj%edge_blk, dim=1)
    dace_rich_obj%f2dace_SA_edge_blk_d_1_s = size(fortran_obj%edge_blk, dim=2)
    dace_rich_obj%f2dace_SOA_edge_blk_d_1_s = lbound(fortran_obj%edge_blk, dim=2)
    dace_rich_obj%f2dace_SA_edge_blk_d_2_s = size(fortran_obj%edge_blk, dim=3)
    dace_rich_obj%f2dace_SOA_edge_blk_d_2_s = lbound(fortran_obj%edge_blk, dim=3)
    dace_rich_obj%f2dace_SA_area_d_0_s = size(fortran_obj%area, dim=1)
    dace_rich_obj%f2dace_SOA_area_d_0_s = lbound(fortran_obj%area, dim=1)
    dace_rich_obj%f2dace_SA_area_d_1_s = size(fortran_obj%area, dim=2)
    dace_rich_obj%f2dace_SOA_area_d_1_s = lbound(fortran_obj%area, dim=2)
    dace_rich_obj%f2dace_SA_start_index_d_0_s = size(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SOA_start_index_d_0_s = lbound(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SA_end_index_d_0_s = size(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SOA_end_index_d_0_s = lbound(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SA_start_blk_d_0_s = size(fortran_obj%start_blk, dim=1)
    dace_rich_obj%f2dace_SOA_start_blk_d_0_s = lbound(fortran_obj%start_blk, dim=1)
    dace_rich_obj%f2dace_SA_start_blk_d_1_s = size(fortran_obj%start_blk, dim=2)
    dace_rich_obj%f2dace_SOA_start_blk_d_1_s = lbound(fortran_obj%start_blk, dim=2)
    dace_rich_obj%f2dace_SA_start_block_d_0_s = size(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SOA_start_block_d_0_s = lbound(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SA_end_blk_d_0_s = size(fortran_obj%end_blk, dim=1)
    dace_rich_obj%f2dace_SOA_end_blk_d_0_s = lbound(fortran_obj%end_blk, dim=1)
    dace_rich_obj%f2dace_SA_end_blk_d_1_s = size(fortran_obj%end_blk, dim=2)
    dace_rich_obj%f2dace_SOA_end_blk_d_1_s = lbound(fortran_obj%end_blk, dim=2)
    dace_rich_obj%f2dace_SA_end_block_d_0_s = size(fortran_obj%end_block, dim=1)
    dace_rich_obj%f2dace_SOA_end_block_d_0_s = lbound(fortran_obj%end_block, dim=1)
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
    dace_rich_obj%start_blk = copy_in_int32_2d_array( &
    fortran_array=logical_fix_2d(fortran_obj%start_blk), &
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
    dace_rich_obj%end_blk = copy_in_int32_2d_array( &
    fortran_array=logical_fix_2d(fortran_obj%end_blk), &
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

    dace_rich_obj%f2dace_SA_owner_mask_d_0_s = size(fortran_obj%owner_mask, dim=1)
    dace_rich_obj%f2dace_SOA_owner_mask_d_0_s = lbound(fortran_obj%owner_mask, dim=1)
    dace_rich_obj%f2dace_SA_owner_mask_d_1_s = size(fortran_obj%owner_mask, dim=2)
    dace_rich_obj%f2dace_SOA_owner_mask_d_1_s = lbound(fortran_obj%owner_mask, dim=2)
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

    dace_rich_obj%f2dace_SA_cell_idx_d_0_s = size(fortran_obj%cell_idx, dim=1)
    dace_rich_obj%f2dace_SOA_cell_idx_d_0_s = lbound(fortran_obj%cell_idx, dim=1)
    dace_rich_obj%f2dace_SA_cell_idx_d_1_s = size(fortran_obj%cell_idx, dim=2)
    dace_rich_obj%f2dace_SOA_cell_idx_d_1_s = lbound(fortran_obj%cell_idx, dim=2)
    dace_rich_obj%f2dace_SA_cell_idx_d_2_s = size(fortran_obj%cell_idx, dim=3)
    dace_rich_obj%f2dace_SOA_cell_idx_d_2_s = lbound(fortran_obj%cell_idx, dim=3)
    dace_rich_obj%f2dace_SA_cell_blk_d_0_s = size(fortran_obj%cell_blk, dim=1)
    dace_rich_obj%f2dace_SOA_cell_blk_d_0_s = lbound(fortran_obj%cell_blk, dim=1)
    dace_rich_obj%f2dace_SA_cell_blk_d_1_s = size(fortran_obj%cell_blk, dim=2)
    dace_rich_obj%f2dace_SOA_cell_blk_d_1_s = lbound(fortran_obj%cell_blk, dim=2)
    dace_rich_obj%f2dace_SA_cell_blk_d_2_s = size(fortran_obj%cell_blk, dim=3)
    dace_rich_obj%f2dace_SOA_cell_blk_d_2_s = lbound(fortran_obj%cell_blk, dim=3)
    dace_rich_obj%f2dace_SA_vertex_idx_d_0_s = size(fortran_obj%vertex_idx, dim=1)
    dace_rich_obj%f2dace_SOA_vertex_idx_d_0_s = lbound(fortran_obj%vertex_idx, dim=1)
    dace_rich_obj%f2dace_SA_vertex_idx_d_1_s = size(fortran_obj%vertex_idx, dim=2)
    dace_rich_obj%f2dace_SOA_vertex_idx_d_1_s = lbound(fortran_obj%vertex_idx, dim=2)
    dace_rich_obj%f2dace_SA_vertex_idx_d_2_s = size(fortran_obj%vertex_idx, dim=3)
    dace_rich_obj%f2dace_SOA_vertex_idx_d_2_s = lbound(fortran_obj%vertex_idx, dim=3)
    dace_rich_obj%f2dace_SA_vertex_blk_d_0_s = size(fortran_obj%vertex_blk, dim=1)
    dace_rich_obj%f2dace_SOA_vertex_blk_d_0_s = lbound(fortran_obj%vertex_blk, dim=1)
    dace_rich_obj%f2dace_SA_vertex_blk_d_1_s = size(fortran_obj%vertex_blk, dim=2)
    dace_rich_obj%f2dace_SOA_vertex_blk_d_1_s = lbound(fortran_obj%vertex_blk, dim=2)
    dace_rich_obj%f2dace_SA_vertex_blk_d_2_s = size(fortran_obj%vertex_blk, dim=3)
    dace_rich_obj%f2dace_SOA_vertex_blk_d_2_s = lbound(fortran_obj%vertex_blk, dim=3)
    dace_rich_obj%f2dace_SA_tangent_orientation_d_0_s = size(fortran_obj%tangent_orientation, dim=1)
    dace_rich_obj%f2dace_SOA_tangent_orientation_d_0_s = lbound(fortran_obj%tangent_orientation, dim=1)
    dace_rich_obj%f2dace_SA_tangent_orientation_d_1_s = size(fortran_obj%tangent_orientation, dim=2)
    dace_rich_obj%f2dace_SOA_tangent_orientation_d_1_s = lbound(fortran_obj%tangent_orientation, dim=2)
    dace_rich_obj%f2dace_SA_quad_idx_d_0_s = size(fortran_obj%quad_idx, dim=1)
    dace_rich_obj%f2dace_SOA_quad_idx_d_0_s = lbound(fortran_obj%quad_idx, dim=1)
    dace_rich_obj%f2dace_SA_quad_idx_d_1_s = size(fortran_obj%quad_idx, dim=2)
    dace_rich_obj%f2dace_SOA_quad_idx_d_1_s = lbound(fortran_obj%quad_idx, dim=2)
    dace_rich_obj%f2dace_SA_quad_idx_d_2_s = size(fortran_obj%quad_idx, dim=3)
    dace_rich_obj%f2dace_SOA_quad_idx_d_2_s = lbound(fortran_obj%quad_idx, dim=3)
    dace_rich_obj%f2dace_SA_quad_blk_d_0_s = size(fortran_obj%quad_blk, dim=1)
    dace_rich_obj%f2dace_SOA_quad_blk_d_0_s = lbound(fortran_obj%quad_blk, dim=1)
    dace_rich_obj%f2dace_SA_quad_blk_d_1_s = size(fortran_obj%quad_blk, dim=2)
    dace_rich_obj%f2dace_SOA_quad_blk_d_1_s = lbound(fortran_obj%quad_blk, dim=2)
    dace_rich_obj%f2dace_SA_quad_blk_d_2_s = size(fortran_obj%quad_blk, dim=3)
    dace_rich_obj%f2dace_SOA_quad_blk_d_2_s = lbound(fortran_obj%quad_blk, dim=3)
    dace_rich_obj%f2dace_SA_primal_normal_cell_d_0_s = size(fortran_obj%primal_normal_cell, dim=1)
    dace_rich_obj%f2dace_SOA_primal_normal_cell_d_0_s = lbound(fortran_obj%primal_normal_cell, dim=1)
    dace_rich_obj%f2dace_SA_primal_normal_cell_d_1_s = size(fortran_obj%primal_normal_cell, dim=2)
    dace_rich_obj%f2dace_SOA_primal_normal_cell_d_1_s = lbound(fortran_obj%primal_normal_cell, dim=2)
    dace_rich_obj%f2dace_SA_primal_normal_cell_d_2_s = size(fortran_obj%primal_normal_cell, dim=3)
    dace_rich_obj%f2dace_SOA_primal_normal_cell_d_2_s = lbound(fortran_obj%primal_normal_cell, dim=3)
    dace_rich_obj%f2dace_SA_dual_normal_cell_d_0_s = size(fortran_obj%dual_normal_cell, dim=1)
    dace_rich_obj%f2dace_SOA_dual_normal_cell_d_0_s = lbound(fortran_obj%dual_normal_cell, dim=1)
    dace_rich_obj%f2dace_SA_dual_normal_cell_d_1_s = size(fortran_obj%dual_normal_cell, dim=2)
    dace_rich_obj%f2dace_SOA_dual_normal_cell_d_1_s = lbound(fortran_obj%dual_normal_cell, dim=2)
    dace_rich_obj%f2dace_SA_dual_normal_cell_d_2_s = size(fortran_obj%dual_normal_cell, dim=3)
    dace_rich_obj%f2dace_SOA_dual_normal_cell_d_2_s = lbound(fortran_obj%dual_normal_cell, dim=3)
    dace_rich_obj%f2dace_SA_inv_primal_edge_length_d_0_s = size(fortran_obj%inv_primal_edge_length, dim=1)
    dace_rich_obj%f2dace_SOA_inv_primal_edge_length_d_0_s = lbound(fortran_obj%inv_primal_edge_length, dim=1)
    dace_rich_obj%f2dace_SA_inv_primal_edge_length_d_1_s = size(fortran_obj%inv_primal_edge_length, dim=2)
    dace_rich_obj%f2dace_SOA_inv_primal_edge_length_d_1_s = lbound(fortran_obj%inv_primal_edge_length, dim=2)
    dace_rich_obj%f2dace_SA_inv_dual_edge_length_d_0_s = size(fortran_obj%inv_dual_edge_length, dim=1)
    dace_rich_obj%f2dace_SOA_inv_dual_edge_length_d_0_s = lbound(fortran_obj%inv_dual_edge_length, dim=1)
    dace_rich_obj%f2dace_SA_inv_dual_edge_length_d_1_s = size(fortran_obj%inv_dual_edge_length, dim=2)
    dace_rich_obj%f2dace_SOA_inv_dual_edge_length_d_1_s = lbound(fortran_obj%inv_dual_edge_length, dim=2)
    dace_rich_obj%f2dace_SA_area_edge_d_0_s = size(fortran_obj%area_edge, dim=1)
    dace_rich_obj%f2dace_SOA_area_edge_d_0_s = lbound(fortran_obj%area_edge, dim=1)
    dace_rich_obj%f2dace_SA_area_edge_d_1_s = size(fortran_obj%area_edge, dim=2)
    dace_rich_obj%f2dace_SOA_area_edge_d_1_s = lbound(fortran_obj%area_edge, dim=2)
    dace_rich_obj%f2dace_SA_f_e_d_0_s = size(fortran_obj%f_e, dim=1)
    dace_rich_obj%f2dace_SOA_f_e_d_0_s = lbound(fortran_obj%f_e, dim=1)
    dace_rich_obj%f2dace_SA_f_e_d_1_s = size(fortran_obj%f_e, dim=2)
    dace_rich_obj%f2dace_SOA_f_e_d_1_s = lbound(fortran_obj%f_e, dim=2)
    dace_rich_obj%f2dace_SA_fn_e_d_0_s = size(fortran_obj%fn_e, dim=1)
    dace_rich_obj%f2dace_SOA_fn_e_d_0_s = lbound(fortran_obj%fn_e, dim=1)
    dace_rich_obj%f2dace_SA_fn_e_d_1_s = size(fortran_obj%fn_e, dim=2)
    dace_rich_obj%f2dace_SOA_fn_e_d_1_s = lbound(fortran_obj%fn_e, dim=2)
    dace_rich_obj%f2dace_SA_ft_e_d_0_s = size(fortran_obj%ft_e, dim=1)
    dace_rich_obj%f2dace_SOA_ft_e_d_0_s = lbound(fortran_obj%ft_e, dim=1)
    dace_rich_obj%f2dace_SA_ft_e_d_1_s = size(fortran_obj%ft_e, dim=2)
    dace_rich_obj%f2dace_SOA_ft_e_d_1_s = lbound(fortran_obj%ft_e, dim=2)
    dace_rich_obj%f2dace_SA_refin_ctrl_d_0_s = size(fortran_obj%refin_ctrl, dim=1)
    dace_rich_obj%f2dace_SOA_refin_ctrl_d_0_s = lbound(fortran_obj%refin_ctrl, dim=1)
    dace_rich_obj%f2dace_SA_refin_ctrl_d_1_s = size(fortran_obj%refin_ctrl, dim=2)
    dace_rich_obj%f2dace_SOA_refin_ctrl_d_1_s = lbound(fortran_obj%refin_ctrl, dim=2)
    dace_rich_obj%f2dace_SA_start_index_d_0_s = size(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SOA_start_index_d_0_s = lbound(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SA_end_index_d_0_s = size(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SOA_end_index_d_0_s = lbound(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SA_start_block_d_0_s = size(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SOA_start_block_d_0_s = lbound(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SA_end_block_d_0_s = size(fortran_obj%end_block, dim=1)
    dace_rich_obj%f2dace_SOA_end_block_d_0_s = lbound(fortran_obj%end_block, dim=1)
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
    dace_rich_obj%primal_normal_cell = copy_in_t_tangent_vectors_3d_array( &
    fortran_array=fortran_obj%primal_normal_cell, &
    steal_arrays=steal_arrays, &
    minimal_structs=minimal_structs &
  )
    dace_rich_obj%dual_normal_cell = copy_in_t_tangent_vectors_3d_array( &
    fortran_array=fortran_obj%dual_normal_cell, &
    steal_arrays=steal_arrays, &
    minimal_structs=minimal_structs &
  )
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
#ifndef _OPENACC
    dace_rich_obj%refin_ctrl = copy_in_int32_2d_array( &
    fortran_array=logical_fix_2d(fortran_obj%refin_ctrl), &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%refin_ctrl = copy_in_int32_2d_array( &
    fortran_array=logical_fix_2d(fortran_obj%refin_ctrl), &
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

  end function copy_in_t_grid_edges

  function copy_in_t_tangent_vectors(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_tangent_vectors), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_tangent_vectors), pointer :: dace_rich_obj
    type(dace_t_tangent_vectors) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%v1 = fortran_obj%v1
    dace_rich_obj%v2 = fortran_obj%v2

  end function copy_in_t_tangent_vectors

  function copy_in_t_grid_vertices(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_grid_vertices), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_grid_vertices), pointer :: dace_rich_obj
    type(dace_t_grid_vertices) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_cell_idx_d_0_s = size(fortran_obj%cell_idx, dim=1)
    dace_rich_obj%f2dace_SOA_cell_idx_d_0_s = lbound(fortran_obj%cell_idx, dim=1)
    dace_rich_obj%f2dace_SA_cell_idx_d_1_s = size(fortran_obj%cell_idx, dim=2)
    dace_rich_obj%f2dace_SOA_cell_idx_d_1_s = lbound(fortran_obj%cell_idx, dim=2)
    dace_rich_obj%f2dace_SA_cell_idx_d_2_s = size(fortran_obj%cell_idx, dim=3)
    dace_rich_obj%f2dace_SOA_cell_idx_d_2_s = lbound(fortran_obj%cell_idx, dim=3)
    dace_rich_obj%f2dace_SA_cell_blk_d_0_s = size(fortran_obj%cell_blk, dim=1)
    dace_rich_obj%f2dace_SOA_cell_blk_d_0_s = lbound(fortran_obj%cell_blk, dim=1)
    dace_rich_obj%f2dace_SA_cell_blk_d_1_s = size(fortran_obj%cell_blk, dim=2)
    dace_rich_obj%f2dace_SOA_cell_blk_d_1_s = lbound(fortran_obj%cell_blk, dim=2)
    dace_rich_obj%f2dace_SA_cell_blk_d_2_s = size(fortran_obj%cell_blk, dim=3)
    dace_rich_obj%f2dace_SOA_cell_blk_d_2_s = lbound(fortran_obj%cell_blk, dim=3)
    dace_rich_obj%f2dace_SA_edge_idx_d_0_s = size(fortran_obj%edge_idx, dim=1)
    dace_rich_obj%f2dace_SOA_edge_idx_d_0_s = lbound(fortran_obj%edge_idx, dim=1)
    dace_rich_obj%f2dace_SA_edge_idx_d_1_s = size(fortran_obj%edge_idx, dim=2)
    dace_rich_obj%f2dace_SOA_edge_idx_d_1_s = lbound(fortran_obj%edge_idx, dim=2)
    dace_rich_obj%f2dace_SA_edge_idx_d_2_s = size(fortran_obj%edge_idx, dim=3)
    dace_rich_obj%f2dace_SOA_edge_idx_d_2_s = lbound(fortran_obj%edge_idx, dim=3)
    dace_rich_obj%f2dace_SA_edge_blk_d_0_s = size(fortran_obj%edge_blk, dim=1)
    dace_rich_obj%f2dace_SOA_edge_blk_d_0_s = lbound(fortran_obj%edge_blk, dim=1)
    dace_rich_obj%f2dace_SA_edge_blk_d_1_s = size(fortran_obj%edge_blk, dim=2)
    dace_rich_obj%f2dace_SOA_edge_blk_d_1_s = lbound(fortran_obj%edge_blk, dim=2)
    dace_rich_obj%f2dace_SA_edge_blk_d_2_s = size(fortran_obj%edge_blk, dim=3)
    dace_rich_obj%f2dace_SOA_edge_blk_d_2_s = lbound(fortran_obj%edge_blk, dim=3)
    dace_rich_obj%f2dace_SA_start_index_d_0_s = size(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SOA_start_index_d_0_s = lbound(fortran_obj%start_index, dim=1)
    dace_rich_obj%f2dace_SA_end_index_d_0_s = size(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SOA_end_index_d_0_s = lbound(fortran_obj%end_index, dim=1)
    dace_rich_obj%f2dace_SA_start_block_d_0_s = size(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SOA_start_block_d_0_s = lbound(fortran_obj%start_block, dim=1)
    dace_rich_obj%f2dace_SA_end_block_d_0_s = size(fortran_obj%end_block, dim=1)
    dace_rich_obj%f2dace_SOA_end_block_d_0_s = lbound(fortran_obj%end_block, dim=1)
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

  function copy_in_t_prepare_adv(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    type(t_prepare_adv), target :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type(dace_t_prepare_adv), pointer :: dace_rich_obj
    type(dace_t_prepare_adv) :: dace_c_obj

    dace_obj_ptr = malloc(c_sizeof(dace_c_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    dace_rich_obj%f2dace_SA_mass_flx_me_d_0_s = size(fortran_obj%mass_flx_me, dim=1)
    dace_rich_obj%f2dace_SOA_mass_flx_me_d_0_s = lbound(fortran_obj%mass_flx_me, dim=1)
    dace_rich_obj%f2dace_SA_mass_flx_me_d_1_s = size(fortran_obj%mass_flx_me, dim=2)
    dace_rich_obj%f2dace_SOA_mass_flx_me_d_1_s = lbound(fortran_obj%mass_flx_me, dim=2)
    dace_rich_obj%f2dace_SA_mass_flx_me_d_2_s = size(fortran_obj%mass_flx_me, dim=3)
    dace_rich_obj%f2dace_SOA_mass_flx_me_d_2_s = lbound(fortran_obj%mass_flx_me, dim=3)
    dace_rich_obj%f2dace_SA_mass_flx_ic_d_0_s = size(fortran_obj%mass_flx_ic, dim=1)
    dace_rich_obj%f2dace_SOA_mass_flx_ic_d_0_s = lbound(fortran_obj%mass_flx_ic, dim=1)
    dace_rich_obj%f2dace_SA_mass_flx_ic_d_1_s = size(fortran_obj%mass_flx_ic, dim=2)
    dace_rich_obj%f2dace_SOA_mass_flx_ic_d_1_s = lbound(fortran_obj%mass_flx_ic, dim=2)
    dace_rich_obj%f2dace_SA_mass_flx_ic_d_2_s = size(fortran_obj%mass_flx_ic, dim=3)
    dace_rich_obj%f2dace_SOA_mass_flx_ic_d_2_s = lbound(fortran_obj%mass_flx_ic, dim=3)
    dace_rich_obj%f2dace_SA_vol_flx_ic_d_0_s = size(fortran_obj%vol_flx_ic, dim=1)
    dace_rich_obj%f2dace_SOA_vol_flx_ic_d_0_s = lbound(fortran_obj%vol_flx_ic, dim=1)
    dace_rich_obj%f2dace_SA_vol_flx_ic_d_1_s = size(fortran_obj%vol_flx_ic, dim=2)
    dace_rich_obj%f2dace_SOA_vol_flx_ic_d_1_s = lbound(fortran_obj%vol_flx_ic, dim=2)
    dace_rich_obj%f2dace_SA_vol_flx_ic_d_2_s = size(fortran_obj%vol_flx_ic, dim=3)
    dace_rich_obj%f2dace_SOA_vol_flx_ic_d_2_s = lbound(fortran_obj%vol_flx_ic, dim=3)
    dace_rich_obj%f2dace_SA_vn_traj_d_0_s = size(fortran_obj%vn_traj, dim=1)
    dace_rich_obj%f2dace_SOA_vn_traj_d_0_s = lbound(fortran_obj%vn_traj, dim=1)
    dace_rich_obj%f2dace_SA_vn_traj_d_1_s = size(fortran_obj%vn_traj, dim=2)
    dace_rich_obj%f2dace_SOA_vn_traj_d_1_s = lbound(fortran_obj%vn_traj, dim=2)
    dace_rich_obj%f2dace_SA_vn_traj_d_2_s = size(fortran_obj%vn_traj, dim=3)
    dace_rich_obj%f2dace_SOA_vn_traj_d_2_s = lbound(fortran_obj%vn_traj, dim=3)
#ifndef _OPENACC
    dace_rich_obj%mass_flx_me = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mass_flx_me, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%mass_flx_me = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mass_flx_me, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%mass_flx_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mass_flx_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%mass_flx_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%mass_flx_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%vol_flx_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vol_flx_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vol_flx_ic = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vol_flx_ic, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif
#ifndef _OPENACC
    dace_rich_obj%vn_traj = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_traj, &
    steal_arrays=steal_arrays, &
    use_openacc=.false., &
    minimal_structs=minimal_structs &
  )

#else
    dace_rich_obj%vn_traj = copy_in_float64_3d_array( &
    fortran_array=fortran_obj%vn_traj, &
    steal_arrays=steal_arrays, &
    use_openacc=.true., &
    minimal_structs=minimal_structs &
  )
#endif

  end function copy_in_t_prepare_adv

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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      .not. c_associated(c_loc(fortran_array)) &
      .or. &
      transfer(c_loc(fortran_array), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      .not. c_associated(c_loc(fortran_array)) &
      .or. &
      transfer(c_loc(fortran_array), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      .not. c_associated(c_loc(fortran_array)) &
      .or. &
      transfer(c_loc(fortran_array), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      .not. c_associated(c_loc(fortran_array)) &
      .or. &
      transfer(c_loc(fortran_array), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      .not. c_associated(c_loc(fortran_array)) &
      .or. &
      transfer(c_loc(fortran_array), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
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

  function logical_to_int_4d(inp) result(out)
    logical(4), dimension(:,:,:,:), target :: inp
    integer(kind=c_int), dimension(:,:,:,:), pointer :: out

    call c_f_pointer(c_loc(inp), out, shape=shape(inp))
  end function logical_to_int_4d

  function int_to_int_4d(inp) result(out)
    integer(kind=c_int), dimension(:,:,:,:), target :: inp
    integer(kind=c_int), dimension(:,:,:,:), pointer :: out

    out => inp
  end function int_to_int_4d

  function copy_in_int32_4d_array( &
    fortran_array, &
    steal_arrays, &
    use_openacc, &
    minimal_structs &
   ) &
   result(dace_array_ptr)
    integer(kind=c_int), dimension(:,:,:,:), target :: fortran_array
    logical :: steal_arrays, use_openacc, minimal_structs
    type(c_ptr) :: dace_array_ptr
    integer(kind=c_int), dimension(:,:,:,:), pointer :: dace_rich_array

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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      .not. c_associated(c_loc(fortran_array)) &
      .or. &
      transfer(c_loc(fortran_array), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
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
  end function copy_in_int32_4d_array

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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      .not. c_associated(c_loc(fortran_array)) &
      .or. &
      transfer(c_loc(fortran_array), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      .not. c_associated(c_loc(fortran_array)) &
      .or. &
      transfer(c_loc(fortran_array), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
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

  ! requires special handling
  function copy_in_t_tangent_vectors_3d_array( &
    fortran_array, &
    steal_arrays, &
    minimal_structs &
  ) &
  result(dace_array_ptr)
    type(t_tangent_vectors), dimension(:,:,:), target :: fortran_array
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_array_ptr
    real(kind=c_double), dimension(:,:,:,:), pointer :: dace_rich_array
    integer(kind=c_size_t) :: size_bytes

    integer :: i0, i1, i2

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      .not. c_associated(c_loc(fortran_array)) &
      .or. &
      transfer(c_loc(fortran_array), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
      dace_array_ptr = c_null_ptr
      return
    end if

    size_bytes = 2 * size(fortran_array) * c_sizeof(dace_rich_array(1, 1, 1, 1))
#ifndef _OPENACC
    dace_array_ptr = malloc(size_bytes)
#else
    dace_array_ptr = c_acc_malloc(size_bytes)
#endif

    call c_f_pointer(dace_array_ptr, dace_rich_array, &
      shape=[ &
        size(fortran_array, dim=1), &
        size(fortran_array, dim=2), &
        size(fortran_array, dim=3), &
        2 &
      ] &
    )

    !$ACC PARALLEL &
    !$ACC   DEFAULT(PRESENT) &
    !$ACC   DEVICEPTR(dace_rich_array)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    do i0 = 1, size(fortran_array, dim=1)
      do i1 = 1, size(fortran_array, dim=2)
        do i2 = 1, size(fortran_array, dim=3)
          dace_rich_array(i0, i1, i2, 1) = fortran_array(i0, i1, i2)%v1
          dace_rich_array(i0, i1, i2, 2) = fortran_array(i0, i1, i2)%v2
        end do
      end do
    end do
    !$ACC END PARALLEL

  end function copy_in_t_tangent_vectors_3d_array

  subroutine copy_back_global_data_type(dace_obj_ptr)
    type(c_ptr) :: dace_obj_ptr

    type(dace_global_data_type), pointer :: dace_rich_obj

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    ldeepatmo = transfer(dace_rich_obj%ldeepatmo, mold=ldeepatmo)
    l_limited_area = transfer(dace_rich_obj%l_limited_area, mold=l_limited_area)
    grf_intmethod_e = transfer(dace_rich_obj%grf_intmethod_e, mold=grf_intmethod_e)
    is_iau_active = transfer(dace_rich_obj%is_iau_active, mold=is_iau_active)
    iau_wgt_dyn = dace_rich_obj%iau_wgt_dyn
    i_am_accel_node = transfer(dace_rich_obj%i_am_accel_node, mold=i_am_accel_node)
    itime_scheme = transfer(dace_rich_obj%itime_scheme, mold=itime_scheme)
    lextra_diffu = transfer(dace_rich_obj%lextra_diffu, mold=lextra_diffu)
    rayleigh_type = transfer(dace_rich_obj%rayleigh_type, mold=rayleigh_type)
    iadv_rhotheta = transfer(dace_rich_obj%iadv_rhotheta, mold=iadv_rhotheta)
    igradp_method = transfer(dace_rich_obj%igradp_method, mold=igradp_method)
    nproma = transfer(dace_rich_obj%nproma, mold=nproma)
    lvert_nest = transfer(dace_rich_obj%lvert_nest, mold=lvert_nest)
    timers_level = transfer(dace_rich_obj%timers_level, mold=timers_level)
    timer_solve_nh_veltend = transfer(dace_rich_obj%timer_solve_nh_veltend, mold=timer_solve_nh_veltend)
    timer_solve_nh_cellcomp = transfer(dace_rich_obj%timer_solve_nh_cellcomp, mold=timer_solve_nh_cellcomp)
    timer_solve_nh_vnupd = transfer(dace_rich_obj%timer_solve_nh_vnupd, mold=timer_solve_nh_vnupd)
    timer_intp = transfer(dace_rich_obj%timer_intp, mold=timer_intp)


    call free(dace_obj_ptr)

  end subroutine copy_back_global_data_type
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
  subroutine copy_back_t_nh_state(fortran_obj, dace_obj_ptr)
    type(t_nh_state), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_nh_state), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_nh_state: Invalid allocation of t_nh_state by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    call copy_back_t_nh_diag(fortran_obj%diag, dace_rich_obj%diag)
    call copy_back_t_nh_ref(fortran_obj%ref, dace_rich_obj%ref)
    call copy_back_t_nh_metrics(fortran_obj%metrics, dace_rich_obj%metrics)


    call free(dace_obj_ptr)

  end subroutine copy_back_t_nh_state
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

    fortran_obj%ddt_vn_dyn_is_associated = transfer(dace_rich_obj%ddt_vn_dyn_is_associated, mold=fortran_obj%ddt_vn_dyn_is_associated)
    fortran_obj%ddt_vn_dmp_is_associated = transfer(dace_rich_obj%ddt_vn_dmp_is_associated, mold=fortran_obj%ddt_vn_dmp_is_associated)
    fortran_obj%ddt_vn_adv_is_associated = transfer(dace_rich_obj%ddt_vn_adv_is_associated, mold=fortran_obj%ddt_vn_adv_is_associated)
    fortran_obj%ddt_vn_cor_is_associated = transfer(dace_rich_obj%ddt_vn_cor_is_associated, mold=fortran_obj%ddt_vn_cor_is_associated)
    fortran_obj%ddt_vn_pgr_is_associated = transfer(dace_rich_obj%ddt_vn_pgr_is_associated, mold=fortran_obj%ddt_vn_pgr_is_associated)
    fortran_obj%ddt_vn_phd_is_associated = transfer(dace_rich_obj%ddt_vn_phd_is_associated, mold=fortran_obj%ddt_vn_phd_is_associated)
    fortran_obj%ddt_vn_iau_is_associated = transfer(dace_rich_obj%ddt_vn_iau_is_associated, mold=fortran_obj%ddt_vn_iau_is_associated)
    fortran_obj%ddt_vn_ray_is_associated = transfer(dace_rich_obj%ddt_vn_ray_is_associated, mold=fortran_obj%ddt_vn_ray_is_associated)
    fortran_obj%ddt_vn_grf_is_associated = transfer(dace_rich_obj%ddt_vn_grf_is_associated, mold=fortran_obj%ddt_vn_grf_is_associated)
    fortran_obj%max_vcfl_dyn = dace_rich_obj%max_vcfl_dyn


    call free(dace_obj_ptr)

  end subroutine copy_back_t_nh_diag
  subroutine copy_back_t_nh_ref(fortran_obj, dace_obj_ptr)
    type(t_nh_ref), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_nh_ref), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_nh_ref: Invalid allocation of t_nh_ref by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



    call free(dace_obj_ptr)

  end subroutine copy_back_t_nh_ref
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

    fortran_obj%pg_listdim = transfer(dace_rich_obj%pg_listdim, mold=fortran_obj%pg_listdim)
    fortran_obj%bdy_mflx_e_dim = transfer(dace_rich_obj%bdy_mflx_e_dim, mold=fortran_obj%bdy_mflx_e_dim)


    call free(dace_obj_ptr)

  end subroutine copy_back_t_nh_metrics
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
  subroutine copy_back_t_patch(fortran_obj, dace_obj_ptr)
    type(t_patch), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_patch), pointer :: dace_rich_obj

    ! No writes, so no copy back
    return

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_patch: Invalid allocation of t_patch by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    fortran_obj%id = transfer(dace_rich_obj%id, mold=fortran_obj%id)
    fortran_obj%n_childdom = transfer(dace_rich_obj%n_childdom, mold=fortran_obj%n_childdom)
    fortran_obj%nblks_c = transfer(dace_rich_obj%nblks_c, mold=fortran_obj%nblks_c)
    fortran_obj%nblks_e = transfer(dace_rich_obj%nblks_e, mold=fortran_obj%nblks_e)
    fortran_obj%nblks_v = transfer(dace_rich_obj%nblks_v, mold=fortran_obj%nblks_v)
    fortran_obj%nlev = transfer(dace_rich_obj%nlev, mold=fortran_obj%nlev)
    fortran_obj%nlevp1 = transfer(dace_rich_obj%nlevp1, mold=fortran_obj%nlevp1)
    fortran_obj%nshift = transfer(dace_rich_obj%nshift, mold=fortran_obj%nshift)
    call copy_back_t_grid_cells(fortran_obj%cells, dace_rich_obj%cells)
    call copy_back_t_grid_edges(fortran_obj%edges, dace_rich_obj%edges)
    call copy_back_t_grid_vertices(fortran_obj%verts, dace_rich_obj%verts)


    ! Keep for proper caching
    !call free(dace_obj_ptr)

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

    call copy_back_t_tangent_vectors_3d_array(fortran_obj%primal_normal_cell, dace_rich_obj%primal_normal_cell)
    call copy_back_t_tangent_vectors_3d_array(fortran_obj%dual_normal_cell, dace_rich_obj%dual_normal_cell)


    call free(dace_obj_ptr)

  end subroutine copy_back_t_grid_edges
  subroutine copy_back_t_tangent_vectors(fortran_obj, dace_obj_ptr)
    type(t_tangent_vectors), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_tangent_vectors), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_tangent_vectors: Invalid allocation of t_tangent_vectors by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

    fortran_obj%v1 = dace_rich_obj%v1
    fortran_obj%v2 = dace_rich_obj%v2


    call free(dace_obj_ptr)

  end subroutine copy_back_t_tangent_vectors
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
  subroutine copy_back_t_prepare_adv(fortran_obj, dace_obj_ptr)
    type(t_prepare_adv), target :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_t_prepare_adv), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_t_prepare_adv: Invalid allocation of t_prepare_adv by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)



    call free(dace_obj_ptr)

  end subroutine copy_back_t_prepare_adv
  ! requires special handling
  subroutine copy_back_t_tangent_vectors_3d_array(fortran_struct_array, dace_struct_array_ptr)
    type(t_tangent_vectors), dimension(:,:,:), target, intent(in) :: fortran_struct_array
    type(c_ptr), intent(in) :: dace_struct_array_ptr

    integer :: i0, i1, i2
    type(c_ptr), dimension(:,:,:), pointer :: dace_struct_array_rich

    if (.not. c_associated(c_loc(fortran_struct_array))) then
      if (c_associated(dace_struct_array_ptr)) then
        print *, "copy_back_t_tangent_vectors_3d_array: Invalid allocation of t_tangent_vectors array by DaCe!"
      end if
      return
    end if

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      transfer(c_loc(fortran_struct_array), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
      if (dace_struct_array_ptr /= c_null_ptr) then
        print *, "copy_back_{base_name}_{rank}d_array: Invalid allocation of {base_name} array by DaCe (ff..fff8)!"
      end if
      return
    end if

    ! skip copy back of `t_tangent_vectors`

#ifndef _OPENACC
    call free(dace_struct_array_ptr)
#else
    call c_acc_free(dace_struct_array_ptr)
#endif

  end subroutine copy_back_t_tangent_vectors_3d_array

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
      "%ldeepatmo"

    call compare_int32_scalar( &
      actual=actual_rich%ldeepatmo, &
      ref=transfer(ldeepatmo, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%l_limited_area"

    call compare_int32_scalar( &
      actual=actual_rich%l_limited_area, &
      ref=transfer(l_limited_area, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%grf_intmethod_e"

    call compare_int32_scalar( &
      actual=actual_rich%grf_intmethod_e, &
      ref=transfer(grf_intmethod_e, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

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
      "%is_iau_active"

    call compare_int32_scalar( &
      actual=actual_rich%is_iau_active, &
      ref=transfer(is_iau_active, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%iau_wgt_dyn"

    call compare_float64_scalar( &
      actual=actual_rich%iau_wgt_dyn, &
      ref=iau_wgt_dyn, &
      result=local_result, &
      scalar_expr=member_expr &
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
      "%itime_scheme"

    call compare_int32_scalar( &
      actual=actual_rich%itime_scheme, &
      ref=transfer(itime_scheme, mold=int(1, kind=4)), &
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
      "%rayleigh_type"

    call compare_int32_scalar( &
      actual=actual_rich%rayleigh_type, &
      ref=transfer(rayleigh_type, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%iadv_rhotheta"

    call compare_int32_scalar( &
      actual=actual_rich%iadv_rhotheta, &
      ref=transfer(iadv_rhotheta, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%igradp_method"

    call compare_int32_scalar( &
      actual=actual_rich%igradp_method, &
      ref=transfer(igradp_method, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%kstart_dd3d"

    call compare_int32_1d_array( &
        actual=actual_rich%kstart_dd3d, &
        ref=logical_fix_1d(kstart_dd3d), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
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
      "%lvert_nest"

    call compare_int32_scalar( &
      actual=actual_rich%lvert_nest, &
      ref=transfer(lvert_nest, mold=int(1, kind=4)), &
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
      "%timer_solve_nh_cellcomp"

    call compare_int32_scalar( &
      actual=actual_rich%timer_solve_nh_cellcomp, &
      ref=transfer(timer_solve_nh_cellcomp, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%timer_solve_nh_vnupd"

    call compare_int32_scalar( &
      actual=actual_rich%timer_solve_nh_vnupd, &
      ref=transfer(timer_solve_nh_vnupd, mold=int(1, kind=4)), &
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

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%nflat_gradp"

    call compare_int32_1d_array( &
        actual=actual_rich%nflat_gradp, &
        ref=logical_fix_1d(nflat_gradp), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result


    call free(actual)

  end subroutine compare_global_data_type_struct

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
      "%e_flx_avg"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%e_flx_avg, &
        ref=ref%e_flx_avg, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%e_flx_avg, &
        ref=ref%e_flx_avg, &
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
      "%geofac_div"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%geofac_div, &
        ref=ref%geofac_div, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%geofac_div, &
        ref=ref%geofac_div, &
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

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%geofac_grg"
#ifndef _OPENACC

    call compare_float64_4d_array( &
        actual=actual_rich%geofac_grg, &
        ref=ref%geofac_grg, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_4d_array( &
        actual=actual_rich%geofac_grg, &
        ref=ref%geofac_grg, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%pos_on_tplane_e"
#ifndef _OPENACC

    call compare_float64_4d_array( &
        actual=actual_rich%pos_on_tplane_e, &
        ref=ref%pos_on_tplane_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_4d_array( &
        actual=actual_rich%pos_on_tplane_e, &
        ref=ref%pos_on_tplane_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%nudgecoeff_e"
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=actual_rich%nudgecoeff_e, &
        ref=ref%nudgecoeff_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_2d_array( &
        actual=actual_rich%nudgecoeff_e, &
        ref=ref%nudgecoeff_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_int_state_struct

  subroutine compare_t_nh_state_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_nh_state), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_nh_state), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%diag"

    call compare_t_nh_diag_struct( &
        actual=actual_rich%diag, &
        ref=ref%diag, &
        result=local_result, &
        struct_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ref"

    call compare_t_nh_ref_struct( &
        actual=actual_rich%ref, &
        ref=ref%ref, &
        result=local_result, &
        struct_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%metrics"

    call compare_t_nh_metrics_struct( &
        actual=actual_rich%metrics, &
        ref=ref%metrics, &
        result=local_result, &
        struct_expr=member_expr &
    )

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_nh_state_struct

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
      "%ddt_vn_dyn_is_associated"

    call compare_int32_scalar( &
      actual=actual_rich%ddt_vn_dyn_is_associated, &
      ref=transfer(ref%ddt_vn_dyn_is_associated, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_dmp_is_associated"

    call compare_int32_scalar( &
      actual=actual_rich%ddt_vn_dmp_is_associated, &
      ref=transfer(ref%ddt_vn_dmp_is_associated, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_adv_is_associated"

    call compare_int32_scalar( &
      actual=actual_rich%ddt_vn_adv_is_associated, &
      ref=transfer(ref%ddt_vn_adv_is_associated, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_cor_is_associated"

    call compare_int32_scalar( &
      actual=actual_rich%ddt_vn_cor_is_associated, &
      ref=transfer(ref%ddt_vn_cor_is_associated, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_pgr_is_associated"

    call compare_int32_scalar( &
      actual=actual_rich%ddt_vn_pgr_is_associated, &
      ref=transfer(ref%ddt_vn_pgr_is_associated, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_phd_is_associated"

    call compare_int32_scalar( &
      actual=actual_rich%ddt_vn_phd_is_associated, &
      ref=transfer(ref%ddt_vn_phd_is_associated, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_iau_is_associated"

    call compare_int32_scalar( &
      actual=actual_rich%ddt_vn_iau_is_associated, &
      ref=transfer(ref%ddt_vn_iau_is_associated, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_ray_is_associated"

    call compare_int32_scalar( &
      actual=actual_rich%ddt_vn_ray_is_associated, &
      ref=transfer(ref%ddt_vn_ray_is_associated, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_grf_is_associated"

    call compare_int32_scalar( &
      actual=actual_rich%ddt_vn_grf_is_associated, &
      ref=transfer(ref%ddt_vn_grf_is_associated, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

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
      "%exner_pr"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%exner_pr, &
        ref=ref%exner_pr, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%exner_pr, &
        ref=ref%exner_pr, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%mass_fl_e"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%mass_fl_e, &
        ref=ref%mass_fl_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%mass_fl_e, &
        ref=ref%mass_fl_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%rho_ic"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%rho_ic, &
        ref=ref%rho_ic, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%rho_ic, &
        ref=ref%rho_ic, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%theta_v_ic"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%theta_v_ic, &
        ref=ref%theta_v_ic, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%theta_v_ic, &
        ref=ref%theta_v_ic, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%grf_tend_vn"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%grf_tend_vn, &
        ref=ref%grf_tend_vn, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%grf_tend_vn, &
        ref=ref%grf_tend_vn, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%grf_tend_w"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%grf_tend_w, &
        ref=ref%grf_tend_w, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%grf_tend_w, &
        ref=ref%grf_tend_w, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%grf_tend_rho"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%grf_tend_rho, &
        ref=ref%grf_tend_rho, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%grf_tend_rho, &
        ref=ref%grf_tend_rho, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%grf_tend_mflx"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%grf_tend_mflx, &
        ref=ref%grf_tend_mflx, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%grf_tend_mflx, &
        ref=ref%grf_tend_mflx, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%grf_bdy_mflx"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%grf_bdy_mflx, &
        ref=ref%grf_bdy_mflx, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%grf_bdy_mflx, &
        ref=ref%grf_bdy_mflx, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%grf_tend_thv"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%grf_tend_thv, &
        ref=ref%grf_tend_thv, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%grf_tend_thv, &
        ref=ref%grf_tend_thv, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vn_ie_int"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%vn_ie_int, &
        ref=ref%vn_ie_int, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%vn_ie_int, &
        ref=ref%vn_ie_int, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vn_ie_ubc"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%vn_ie_ubc, &
        ref=ref%vn_ie_ubc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%vn_ie_ubc, &
        ref=ref%vn_ie_ubc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%w_int"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%w_int, &
        ref=ref%w_int, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%w_int, &
        ref=ref%w_int, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%w_ubc"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%w_ubc, &
        ref=ref%w_ubc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%w_ubc, &
        ref=ref%w_ubc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%theta_v_ic_int"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%theta_v_ic_int, &
        ref=ref%theta_v_ic_int, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%theta_v_ic_int, &
        ref=ref%theta_v_ic_int, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%theta_v_ic_ubc"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%theta_v_ic_ubc, &
        ref=ref%theta_v_ic_ubc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%theta_v_ic_ubc, &
        ref=ref%theta_v_ic_ubc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%rho_ic_int"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%rho_ic_int, &
        ref=ref%rho_ic_int, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%rho_ic_int, &
        ref=ref%rho_ic_int, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%rho_ic_ubc"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%rho_ic_ubc, &
        ref=ref%rho_ic_ubc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%rho_ic_ubc, &
        ref=ref%rho_ic_ubc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%mflx_ic_int"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%mflx_ic_int, &
        ref=ref%mflx_ic_int, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%mflx_ic_int, &
        ref=ref%mflx_ic_int, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%mflx_ic_ubc"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%mflx_ic_ubc, &
        ref=ref%mflx_ic_ubc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%mflx_ic_ubc, &
        ref=ref%mflx_ic_ubc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vn_incr"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%vn_incr, &
        ref=ref%vn_incr, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%vn_incr, &
        ref=ref%vn_incr, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%exner_incr"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%exner_incr, &
        ref=ref%exner_incr, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%exner_incr, &
        ref=ref%exner_incr, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%rho_incr"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%rho_incr, &
        ref=ref%rho_incr, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%rho_incr, &
        ref=ref%rho_incr, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

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
      "%ddt_exner_phy"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_exner_phy, &
        ref=ref%ddt_exner_phy, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_exner_phy, &
        ref=ref%ddt_exner_phy, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_phy"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_phy, &
        ref=ref%ddt_vn_phy, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_phy, &
        ref=ref%ddt_vn_phy, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%exner_dyn_incr"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%exner_dyn_incr, &
        ref=ref%exner_dyn_incr, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%exner_dyn_incr, &
        ref=ref%exner_dyn_incr, &
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
      "%mass_fl_e_sv"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%mass_fl_e_sv, &
        ref=ref%mass_fl_e_sv, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%mass_fl_e_sv, &
        ref=ref%mass_fl_e_sv, &
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
      "%ddt_vn_cor_pc"
#ifndef _OPENACC

    call compare_float64_4d_array( &
        actual=actual_rich%ddt_vn_cor_pc, &
        ref=ref%ddt_vn_cor_pc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_4d_array( &
        actual=actual_rich%ddt_vn_cor_pc, &
        ref=ref%ddt_vn_cor_pc, &
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

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_dyn"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_dyn, &
        ref=ref%ddt_vn_dyn, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_dyn, &
        ref=ref%ddt_vn_dyn, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_dmp"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_dmp, &
        ref=ref%ddt_vn_dmp, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_dmp, &
        ref=ref%ddt_vn_dmp, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_adv"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_adv, &
        ref=ref%ddt_vn_adv, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_adv, &
        ref=ref%ddt_vn_adv, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_cor"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_cor, &
        ref=ref%ddt_vn_cor, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_cor, &
        ref=ref%ddt_vn_cor, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_pgr"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_pgr, &
        ref=ref%ddt_vn_pgr, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_pgr, &
        ref=ref%ddt_vn_pgr, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_phd"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_phd, &
        ref=ref%ddt_vn_phd, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_phd, &
        ref=ref%ddt_vn_phd, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_iau"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_iau, &
        ref=ref%ddt_vn_iau, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_iau, &
        ref=ref%ddt_vn_iau, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_ray"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_ray, &
        ref=ref%ddt_vn_ray, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_ray, &
        ref=ref%ddt_vn_ray, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%ddt_vn_grf"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_grf, &
        ref=ref%ddt_vn_grf, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%ddt_vn_grf, &
        ref=ref%ddt_vn_grf, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_nh_diag_struct

  subroutine compare_t_nh_ref_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_nh_ref), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_nh_ref), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vn_ref"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%vn_ref, &
        ref=ref%vn_ref, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%vn_ref, &
        ref=ref%vn_ref, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%w_ref"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%w_ref, &
        ref=ref%w_ref, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%w_ref, &
        ref=ref%w_ref, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_nh_ref_struct

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
      "%pg_listdim"

    call compare_int32_scalar( &
      actual=actual_rich%pg_listdim, &
      ref=transfer(ref%pg_listdim, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%bdy_mflx_e_dim"

    call compare_int32_scalar( &
      actual=actual_rich%bdy_mflx_e_dim, &
      ref=transfer(ref%bdy_mflx_e_dim, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%rayleigh_w"
#ifndef _OPENACC

    call compare_float64_1d_array( &
        actual=actual_rich%rayleigh_w, &
        ref=ref%rayleigh_w, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_1d_array( &
        actual=actual_rich%rayleigh_w, &
        ref=ref%rayleigh_w, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%rayleigh_vn"
#ifndef _OPENACC

    call compare_float64_1d_array( &
        actual=actual_rich%rayleigh_vn, &
        ref=ref%rayleigh_vn, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_1d_array( &
        actual=actual_rich%rayleigh_vn, &
        ref=ref%rayleigh_vn, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%scalfac_dd3d"
#ifndef _OPENACC

    call compare_float64_1d_array( &
        actual=actual_rich%scalfac_dd3d, &
        ref=ref%scalfac_dd3d, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_1d_array( &
        actual=actual_rich%scalfac_dd3d, &
        ref=ref%scalfac_dd3d, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%hmask_dd3d"
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=actual_rich%hmask_dd3d, &
        ref=ref%hmask_dd3d, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_2d_array( &
        actual=actual_rich%hmask_dd3d, &
        ref=ref%hmask_dd3d, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vwind_expl_wgt"
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=actual_rich%vwind_expl_wgt, &
        ref=ref%vwind_expl_wgt, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_2d_array( &
        actual=actual_rich%vwind_expl_wgt, &
        ref=ref%vwind_expl_wgt, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vwind_impl_wgt"
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=actual_rich%vwind_impl_wgt, &
        ref=ref%vwind_impl_wgt, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_2d_array( &
        actual=actual_rich%vwind_impl_wgt, &
        ref=ref%vwind_impl_wgt, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

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
      "%inv_ddqz_z_full"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%inv_ddqz_z_full, &
        ref=ref%inv_ddqz_z_full, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%inv_ddqz_z_full, &
        ref=ref%inv_ddqz_z_full, &
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
      "%wgtfacq_c"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%wgtfacq_c, &
        ref=ref%wgtfacq_c, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%wgtfacq_c, &
        ref=ref%wgtfacq_c, &
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
      "%wgtfacq1_c"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%wgtfacq1_c, &
        ref=ref%wgtfacq1_c, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%wgtfacq1_c, &
        ref=ref%wgtfacq1_c, &
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

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%zdiff_gradp"
#ifndef _OPENACC

    call compare_float64_4d_array( &
        actual=actual_rich%zdiff_gradp, &
        ref=ref%zdiff_gradp, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_4d_array( &
        actual=actual_rich%zdiff_gradp, &
        ref=ref%zdiff_gradp, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%coeff_gradp"
#ifndef _OPENACC

    call compare_float64_4d_array( &
        actual=actual_rich%coeff_gradp, &
        ref=ref%coeff_gradp, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_4d_array( &
        actual=actual_rich%coeff_gradp, &
        ref=ref%coeff_gradp, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%exner_exfac"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%exner_exfac, &
        ref=ref%exner_exfac, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%exner_exfac, &
        ref=ref%exner_exfac, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%theta_ref_mc"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%theta_ref_mc, &
        ref=ref%theta_ref_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%theta_ref_mc, &
        ref=ref%theta_ref_mc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%theta_ref_me"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%theta_ref_me, &
        ref=ref%theta_ref_me, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%theta_ref_me, &
        ref=ref%theta_ref_me, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%theta_ref_ic"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%theta_ref_ic, &
        ref=ref%theta_ref_ic, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%theta_ref_ic, &
        ref=ref%theta_ref_ic, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%exner_ref_mc"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%exner_ref_mc, &
        ref=ref%exner_ref_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%exner_ref_mc, &
        ref=ref%exner_ref_mc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%rho_ref_mc"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%rho_ref_mc, &
        ref=ref%rho_ref_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%rho_ref_mc, &
        ref=ref%rho_ref_mc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%rho_ref_me"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%rho_ref_me, &
        ref=ref%rho_ref_me, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%rho_ref_me, &
        ref=ref%rho_ref_me, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%d_exner_dz_ref_ic"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%d_exner_dz_ref_ic, &
        ref=ref%d_exner_dz_ref_ic, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%d_exner_dz_ref_ic, &
        ref=ref%d_exner_dz_ref_ic, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%d2dexdz2_fac1_mc"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%d2dexdz2_fac1_mc, &
        ref=ref%d2dexdz2_fac1_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%d2dexdz2_fac1_mc, &
        ref=ref%d2dexdz2_fac1_mc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%d2dexdz2_fac2_mc"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%d2dexdz2_fac2_mc, &
        ref=ref%d2dexdz2_fac2_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%d2dexdz2_fac2_mc, &
        ref=ref%d2dexdz2_fac2_mc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%pg_exdist"

    call compare_float64_1d_array( &
        actual=actual_rich%pg_exdist, &
        ref=ref%pg_exdist, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vertidx_gradp"
#ifndef _OPENACC

    call compare_int32_4d_array( &
        actual=actual_rich%vertidx_gradp, &
        ref=logical_fix_4d(ref%vertidx_gradp), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_4d_array( &
        actual=actual_rich%vertidx_gradp, &
        ref=logical_fix_4d(ref%vertidx_gradp), &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%pg_edgeidx"

    call compare_int32_1d_array( &
        actual=actual_rich%pg_edgeidx, &
        ref=logical_fix_1d(ref%pg_edgeidx), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%pg_edgeblk"

    call compare_int32_1d_array( &
        actual=actual_rich%pg_edgeblk, &
        ref=logical_fix_1d(ref%pg_edgeblk), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%pg_vertidx"

    call compare_int32_1d_array( &
        actual=actual_rich%pg_vertidx, &
        ref=logical_fix_1d(ref%pg_vertidx), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%bdy_mflx_e_idx"

    call compare_int32_1d_array( &
        actual=actual_rich%bdy_mflx_e_idx, &
        ref=logical_fix_1d(ref%bdy_mflx_e_idx), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%bdy_mflx_e_blk"

    call compare_int32_1d_array( &
        actual=actual_rich%bdy_mflx_e_blk, &
        ref=logical_fix_1d(ref%bdy_mflx_e_blk), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%deepatmo_gradh_mc"
#ifndef _OPENACC

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_gradh_mc, &
        ref=ref%deepatmo_gradh_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_gradh_mc, &
        ref=ref%deepatmo_gradh_mc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%deepatmo_divh_mc"
#ifndef _OPENACC

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_divh_mc, &
        ref=ref%deepatmo_divh_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_divh_mc, &
        ref=ref%deepatmo_divh_mc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

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
      "%deepatmo_divzu_mc"
#ifndef _OPENACC

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_divzu_mc, &
        ref=ref%deepatmo_divzu_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_divzu_mc, &
        ref=ref%deepatmo_divzu_mc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%deepatmo_divzl_mc"
#ifndef _OPENACC

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_divzl_mc, &
        ref=ref%deepatmo_divzl_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_1d_array( &
        actual=actual_rich%deepatmo_divzl_mc, &
        ref=ref%deepatmo_divzl_mc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

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


    call free(actual)

  end subroutine compare_t_nh_metrics_struct

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

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%rho"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%rho, &
        ref=ref%rho, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%rho, &
        ref=ref%rho, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%exner"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%exner, &
        ref=ref%exner, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%exner, &
        ref=ref%exner, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%theta_v"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%theta_v, &
        ref=ref%theta_v, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%theta_v, &
        ref=ref%theta_v, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_nh_prog_struct

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
      "%id"

    call compare_int32_scalar( &
      actual=actual_rich%id, &
      ref=transfer(ref%id, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%n_childdom"

    call compare_int32_scalar( &
      actual=actual_rich%n_childdom, &
      ref=transfer(ref%n_childdom, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

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
      "%nlev"

    call compare_int32_scalar( &
      actual=actual_rich%nlev, &
      ref=transfer(ref%nlev, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%nlevp1"

    call compare_int32_scalar( &
      actual=actual_rich%nlevp1, &
      ref=transfer(ref%nlevp1, mold=int(1, kind=4)), &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%nshift"

    call compare_int32_scalar( &
      actual=actual_rich%nshift, &
      ref=transfer(ref%nshift, mold=int(1, kind=4)), &
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


    ! Keep for proper caching
    !call free(actual)

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
      "%start_blk"

    call compare_int32_2d_array( &
        actual=actual_rich%start_blk, &
        ref=logical_fix_2d(ref%start_blk), &
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
      "%end_blk"

    call compare_int32_2d_array( &
        actual=actual_rich%end_blk, &
        ref=logical_fix_2d(ref%end_blk), &
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
      "%primal_normal_cell"

    call compare_t_tangent_vectors_3d_array( &
        actual=actual_rich%primal_normal_cell, &
        ref=ref%primal_normal_cell, &
        result=local_result, &
        struct_array_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%dual_normal_cell"

    call compare_t_tangent_vectors_3d_array( &
        actual=actual_rich%dual_normal_cell, &
        ref=ref%dual_normal_cell, &
        result=local_result, &
        struct_array_expr=member_expr &
    )

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
      "%refin_ctrl"
#ifndef _OPENACC

    call compare_int32_2d_array( &
        actual=actual_rich%refin_ctrl, &
        ref=logical_fix_2d(ref%refin_ctrl), &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_int32_2d_array( &
        actual=actual_rich%refin_ctrl, &
        ref=logical_fix_2d(ref%refin_ctrl), &
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

  end subroutine compare_t_grid_edges_struct

  subroutine compare_t_tangent_vectors_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_tangent_vectors), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_tangent_vectors), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%v1"

    call compare_float64_scalar( &
      actual=actual_rich%v1, &
      ref=ref%v1, &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%v2"

    call compare_float64_scalar( &
      actual=actual_rich%v2, &
      ref=ref%v2, &
      result=local_result, &
      scalar_expr=member_expr &
    )

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_tangent_vectors_struct

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

  subroutine compare_t_prepare_adv_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_prepare_adv), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_t_prepare_adv), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%mass_flx_me"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%mass_flx_me, &
        ref=ref%mass_flx_me, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%mass_flx_me, &
        ref=ref%mass_flx_me, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%mass_flx_ic"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%mass_flx_ic, &
        ref=ref%mass_flx_ic, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%mass_flx_ic, &
        ref=ref%mass_flx_ic, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vol_flx_ic"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%vol_flx_ic, &
        ref=ref%vol_flx_ic, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%vol_flx_ic, &
        ref=ref%vol_flx_ic, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result

    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%vn_traj"
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=actual_rich%vn_traj, &
        ref=ref%vn_traj, &
        result=local_result, &
        use_openacc=.false., &
        array_expr=member_expr &
    )

#else

    call compare_float64_3d_array( &
        actual=actual_rich%vn_traj, &
        ref=ref%vn_traj, &
        result=local_result, &
        use_openacc=.true., &
        array_expr=member_expr &
    )

#endif

    result = result .and. local_result


    call free(actual)

  end subroutine compare_t_prepare_adv_struct

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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      transfer(c_loc(ref), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
      result = actual == c_null_ptr

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not (ff..fff8)!"
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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      transfer(c_loc(ref), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
      result = actual == c_null_ptr

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not (ff..fff8)!"
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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      transfer(c_loc(ref), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
      result = actual == c_null_ptr

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not (ff..fff8)!"
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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      transfer(c_loc(ref), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
      result = actual == c_null_ptr

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not (ff..fff8)!"
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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      transfer(c_loc(ref), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
      result = actual == c_null_ptr

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not (ff..fff8)!"
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

  subroutine compare_int32_4d_array( &
    actual, &
    ref, &
    use_openacc, &
    result, &
    array_expr, &
    rel_threshold, &
    abs_threshold &
  )
    type(c_ptr) :: actual
    integer(kind=c_int), dimension(:,:,:,:), target, intent(in) :: ref
    logical, intent(in) :: use_openacc
    logical, intent(out) :: result
    real(kind=c_double), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in) :: array_expr

    real(kind=c_double) :: actual_rel_threshold, actual_abs_threshold
    logical :: local_result
    integer :: i0, i1, i2, i3
    integer(kind=c_int), dimension(:,:,:,:), pointer :: actual_rich
    CHARACTER(len=5000) :: message_text = ''

    integer(kind=c_int) :: error_ref, error_actual
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
        print *, "compare_int32_4d_array"
        print *, trim(message_text)
      end if

      return
    end if

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      transfer(c_loc(ref), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
      result = actual == c_null_ptr

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not (ff..fff8)!"
        print *, "compare_int32_4d_array"
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
          do i3 = 1, size(ref, dim=4)

    call compare_int32_scalar( &
      actual=actual_rich(i0, i1, i2, i3), &
      ref=transfer(ref(i0, i1, i2, i3), mold=int(1, kind=4)), &
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
       
    call compare_int32_scalar( &
      actual=actual_rich(i0, i1, i2, i3), &
      ref=transfer(ref(i0, i1, i2, i3), mold=int(1, kind=4)), &
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
        print *, "compare_int32_4d_array"
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

  end subroutine compare_int32_4d_array

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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      transfer(c_loc(ref), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
      result = actual == c_null_ptr

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not (ff..fff8)!"
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

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      transfer(c_loc(ref), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
      result = actual == c_null_ptr

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not (ff..fff8)!"
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
  ! requires special handling
  subroutine compare_t_tangent_vectors_3d_array( &
    actual, &
    ref, &
    result, &
    struct_array_expr &
  )
    type(c_ptr), intent(in) :: actual
    type(t_tangent_vectors), dimension(:,:,:), target, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_array_expr

    CHARACTER(len=5000) :: member_expr = ''
    CHARACTER(len=5000) :: message_text = ''
    logical :: local_result
    integer :: i0, i1, i2
    type(c_ptr), dimension(:,:,:), pointer :: actual_rich

    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(struct_array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not!"
        print *, "compare_t_tangent_vectors_3d_array"
        print *, trim(message_text)
      end if

      return
    end if

    ! ff..fff8 seems to be some kind of magic value for nvfortran (possibly together with OpenACC)
    if ( &
      transfer(c_loc(ref), mold=int(1, kind=c_intptr_t)) == int(Z'fffffffffffffff8', kind=c_intptr_t) &
    ) then
      result = actual == c_null_ptr

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(struct_array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not (ff..fff8)!"
        print *, "compare_t_tangent_vectors_3d_array"
        print *, trim(message_text)
      end if

      return
    end if

    ! skip comparison of `t_tangent_vectors`

#ifndef _OPENACC
    call free(actual)
#else
    call c_acc_free(actual)
#endif

  end subroutine compare_t_tangent_vectors_3d_array



  subroutine run_solve_nh_predictor_pre( &
    bdy_divdamp, &
    enh_divdamp_fac, &
    p_int, &
    p_nh, &
    p_nh_prog_nnew, &
    p_nh_prog_nnow, &
    p_patch, &
    prep_adv, &
    scal_divdamp, &
    z_alpha, &
    z_beta, &
    z_contr_w_fl_l, &
    z_dexner_dz_c, &
    z_dwdz_dd, &
    z_exner_ex_pr, &
    z_exner_expl, &
    z_exner_ic, &
    z_flxdiv_mass, &
    z_flxdiv_theta, &
    z_grad_rth, &
    z_graddiv2_vn, &
    z_graddiv_vn, &
    z_gradh_exner, &
    z_hydro_corr, &
    z_kin_hor_e, &
    z_mflx_top, &
    z_q, &
    z_raylfac, &
    z_rho_e, &
    z_rho_expl, &
    z_rho_v, &
    z_rth_pr, &
    z_th_ddz_exner_c, &
    z_theta_v_e, &
    z_theta_v_fl_e, &
    z_theta_v_pr_ic, &
    z_theta_v_v, &
    z_vn_avg, &
    z_vt_ie, &
    z_w_concorr_mc, &
    z_w_concorr_me, &
    z_w_expl, &
    alin, &
    aqdr, &
    bqdr, &
    df32, &
    df42, &
    distv_bary_1, &
    distv_bary_2, &
    dt_linintp_ubc, &
    dt_linintp_ubc_nnew, &
    dt_linintp_ubc_nnow, &
    dt_shift, &
    dthalf, &
    dtime, &
    dz32, &
    dz42, &
    dzlin, &
    dzqdr, &
    i_endblk, &
    i_endidx, &
    i_startblk, &
    i_startidx, &
    idyn_timestep, &
    ishift, &
    istep, &
    jb, &
    jc, &
    je, &
    jg, &
    jk, &
    jk_start, &
    jks, &
    jstep, &
    l_child_vertnest, &
    l_init, &
    l_recompute, &
    l_vert_nested, &
    lacc, &
    lclean_mflx, &
    lprep_adv, &
    lsave_mflx, &
    lvn_only, &
    lvn_pos, &
    nblks_gradp, &
    nlen_gradp, &
    nlev, &
    nlevp1, &
    nnew, &
    nnow, &
    nproma_gradp, &
    npromz_gradp, &
    nshift, &
    nshift_total, &
    ntl1, &
    ntl2, &
    nvar, &
    r_dtimensubsteps, &
    r_nsubsteps, &
    rl_end, &
    rl_start, &
    scal_divdamp_o2, &
    wgt_nnew_rth, &
    wgt_nnew_vel, &
    wgt_nnow_rth, &
    wgt_nnow_vel, &
    z_a, &
    z_b, &
    z_c, &
    z_d_vn_dmp, &
    z_d_vn_iau, &
    z_ddt_vn_apc, &
    z_ddt_vn_cor, &
    z_ddt_vn_dyn, &
    z_ddt_vn_pgr, &
    z_ddt_vn_ray, &
    z_g, &
    z_gamma, &
    z_ntdistv_bary_1, &
    z_ntdistv_bary_2, &
    z_rho_tavg, &
    z_rho_tavg_m1, &
    z_theta1, &
    z_theta2, &
    z_theta_tavg, &
    z_theta_tavg_m1, &
    z_theta_v_pr_mc, &
    z_theta_v_pr_mc_m1, &
    z_w_backtraj, &
    zf &
  )
! SOLVE_NH PART TIMERS : PRATYAI
real :: t0, t1
    real(kind=c_double), dimension(:), target :: bdy_divdamp
    real(kind=c_double), dimension(:), target :: enh_divdamp_fac
    type(t_int_state), target :: p_int
    type(t_nh_state), target :: p_nh
    type(t_nh_prog), target :: p_nh_prog_nnew
    type(t_nh_prog), target :: p_nh_prog_nnow
    type(t_patch), target :: p_patch
    type(t_prepare_adv), target :: prep_adv
    real(kind=c_double), dimension(:), target :: scal_divdamp
    real(kind=c_double), dimension(:,:), target :: z_alpha
    real(kind=c_double), dimension(:,:), target :: z_beta
    real(kind=c_double), dimension(:,:), target :: z_contr_w_fl_l
    real(kind=c_double), dimension(:,:,:,:), target :: z_dexner_dz_c
    real(kind=c_double), dimension(:,:,:), target :: z_dwdz_dd
    real(kind=c_double), dimension(:,:,:), target :: z_exner_ex_pr
    real(kind=c_double), dimension(:,:), target :: z_exner_expl
    real(kind=c_double), dimension(:,:), target :: z_exner_ic
    real(kind=c_double), dimension(:,:), target :: z_flxdiv_mass
    real(kind=c_double), dimension(:,:), target :: z_flxdiv_theta
    real(kind=c_double), dimension(:,:,:,:), target :: z_grad_rth
    real(kind=c_double), dimension(:,:), target :: z_graddiv2_vn
    real(kind=c_double), dimension(:,:,:), target :: z_graddiv_vn
    real(kind=c_double), dimension(:,:,:), target :: z_gradh_exner
    real(kind=c_double), dimension(:,:), target :: z_hydro_corr
    real(kind=c_double), dimension(:,:,:), target :: z_kin_hor_e
    real(kind=c_double), dimension(:,:), target :: z_mflx_top
    real(kind=c_double), dimension(:,:), target :: z_q
    real(kind=c_double), dimension(:), target :: z_raylfac
    real(kind=c_double), dimension(:,:,:), target :: z_rho_e
    real(kind=c_double), dimension(:,:), target :: z_rho_expl
    real(kind=c_double), dimension(:,:,:), target :: z_rho_v
    real(kind=c_double), dimension(:,:,:,:), target :: z_rth_pr
    real(kind=c_double), dimension(:,:,:), target :: z_th_ddz_exner_c
    real(kind=c_double), dimension(:,:,:), target :: z_theta_v_e
    real(kind=c_double), dimension(:,:,:), target :: z_theta_v_fl_e
    real(kind=c_double), dimension(:,:), target :: z_theta_v_pr_ic
    real(kind=c_double), dimension(:,:,:), target :: z_theta_v_v
    real(kind=c_double), dimension(:,:), target :: z_vn_avg
    real(kind=c_double), dimension(:,:,:), target :: z_vt_ie
    real(kind=c_double), dimension(:,:), target :: z_w_concorr_mc
    real(kind=c_double), dimension(:,:,:), target :: z_w_concorr_me
    real(kind=c_double), dimension(:,:), target :: z_w_expl
    real(kind=c_double) :: alin
    real(kind=c_double) :: aqdr
    real(kind=c_double) :: bqdr
    real(kind=c_double) :: df32
    real(kind=c_double) :: df42
    real(kind=c_double) :: distv_bary_1
    real(kind=c_double) :: distv_bary_2
    real(kind=c_double) :: dt_linintp_ubc
    real(kind=c_double) :: dt_linintp_ubc_nnew
    real(kind=c_double) :: dt_linintp_ubc_nnow
    real(kind=c_double) :: dt_shift
    real(kind=c_double) :: dthalf
    real(kind=c_double) :: dtime
    real(kind=c_double) :: dz32
    real(kind=c_double) :: dz42
    real(kind=c_double) :: dzlin
    real(kind=c_double) :: dzqdr
    integer(kind=c_int) :: i_endblk
    integer(kind=c_int) :: i_endidx
    integer(kind=c_int) :: i_startblk
    integer(kind=c_int) :: i_startidx
    integer(kind=c_int) :: idyn_timestep
    integer(kind=c_int) :: ishift
    integer(kind=c_int) :: istep
    integer(kind=c_int) :: jb
    integer(kind=c_int) :: jc
    integer(kind=c_int) :: je
    integer(kind=c_int) :: jg
    integer(kind=c_int) :: jk
    integer(kind=c_int) :: jk_start
    integer(kind=c_int) :: jks
    integer(kind=c_int) :: jstep
    integer(kind=c_int) :: l_child_vertnest
    integer(kind=c_int) :: l_init
    integer(kind=c_int) :: l_recompute
    integer(kind=c_int) :: l_vert_nested
    integer(kind=c_int), optional :: lacc
    integer(kind=c_int) :: lclean_mflx
    integer(kind=c_int) :: lprep_adv
    integer(kind=c_int) :: lsave_mflx
    integer(kind=c_int) :: lvn_only
    integer(kind=c_int) :: lvn_pos
    integer(kind=c_int) :: nblks_gradp
    integer(kind=c_int) :: nlen_gradp
    integer(kind=c_int) :: nlev
    integer(kind=c_int) :: nlevp1
    integer(kind=c_int) :: nnew
    integer(kind=c_int) :: nnow
    integer(kind=c_int) :: nproma_gradp
    integer(kind=c_int) :: npromz_gradp
    integer(kind=c_int) :: nshift
    integer(kind=c_int) :: nshift_total
    integer(kind=c_int) :: ntl1
    integer(kind=c_int) :: ntl2
    integer(kind=c_int) :: nvar
    real(kind=c_double) :: r_dtimensubsteps
    real(kind=c_double) :: r_nsubsteps
    integer(kind=c_int) :: rl_end
    integer(kind=c_int) :: rl_start
    real(kind=c_double) :: scal_divdamp_o2
    real(kind=c_double) :: wgt_nnew_rth
    real(kind=c_double) :: wgt_nnew_vel
    real(kind=c_double) :: wgt_nnow_rth
    real(kind=c_double) :: wgt_nnow_vel
    real(kind=c_double) :: z_a
    real(kind=c_double) :: z_b
    real(kind=c_double) :: z_c
    real(kind=c_double) :: z_d_vn_dmp
    real(kind=c_double) :: z_d_vn_iau
    real(kind=c_double) :: z_ddt_vn_apc
    real(kind=c_double) :: z_ddt_vn_cor
    real(kind=c_double) :: z_ddt_vn_dyn
    real(kind=c_double) :: z_ddt_vn_pgr
    real(kind=c_double) :: z_ddt_vn_ray
    real(kind=c_double) :: z_g
    real(kind=c_double) :: z_gamma
    real(kind=c_double) :: z_ntdistv_bary_1
    real(kind=c_double) :: z_ntdistv_bary_2
    real(kind=c_double) :: z_rho_tavg
    real(kind=c_double) :: z_rho_tavg_m1
    real(kind=c_double) :: z_theta1
    real(kind=c_double) :: z_theta2
    real(kind=c_double) :: z_theta_tavg
    real(kind=c_double) :: z_theta_tavg_m1
    real(kind=c_double) :: z_theta_v_pr_mc
    real(kind=c_double) :: z_theta_v_pr_mc_m1
    real(kind=c_double) :: z_w_backtraj
    real(kind=c_double) :: zf
    !Optional helper parameter
    !WARNING: HACKFIX, POSSIBLE EXPLOSION
    integer(kind=c_int) :: f2dace_OPTIONAL_lacc

    integer(kind=c_int) :: DACE_OPT_PROXY_lacc



    if (present(lacc)) then
      f2dace_OPTIONAL_lacc = 1
      DACE_OPT_PROXY_lacc = lacc
      
    else
      f2dace_OPTIONAL_lacc = 0
    end if


#ifndef _OPENACC
    copy_or_ptr_bdy_divdamp = copy_in_float64_1d_array( &
    fortran_array=bdy_divdamp, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_bdy_divdamp = copy_in_float64_1d_array( &
    fortran_array=bdy_divdamp, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
    copy_or_ptr_enh_divdamp_fac = copy_in_float64_1d_array( &
    fortran_array=enh_divdamp_fac, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_global_data = copy_in_global_data_type( &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_int = copy_in_t_int_state( &
    fortran_obj=p_int, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_nh = copy_in_t_nh_state( &
    fortran_obj=p_nh, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_nh_prog_nnew = copy_in_t_nh_prog( &
    fortran_obj=p_nh_prog_nnew, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_nh_prog_nnow = copy_in_t_nh_prog( &
    fortran_obj=p_nh_prog_nnow, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_patch = copy_in_t_patch( &
    fortran_obj=p_patch, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
    copy_or_ptr_prep_adv = copy_in_t_prepare_adv( &
    fortran_obj=prep_adv, &
    steal_arrays=.true., &
    minimal_structs=.false. &
  )
#ifndef _OPENACC
    copy_or_ptr_scal_divdamp = copy_in_float64_1d_array( &
    fortran_array=scal_divdamp, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_scal_divdamp = copy_in_float64_1d_array( &
    fortran_array=scal_divdamp, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_alpha = copy_in_float64_2d_array( &
    fortran_array=z_alpha, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_alpha = copy_in_float64_2d_array( &
    fortran_array=z_alpha, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_beta = copy_in_float64_2d_array( &
    fortran_array=z_beta, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_beta = copy_in_float64_2d_array( &
    fortran_array=z_beta, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_contr_w_fl_l = copy_in_float64_2d_array( &
    fortran_array=z_contr_w_fl_l, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_contr_w_fl_l = copy_in_float64_2d_array( &
    fortran_array=z_contr_w_fl_l, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_dexner_dz_c = copy_in_float64_4d_array( &
    fortran_array=z_dexner_dz_c, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_dexner_dz_c = copy_in_float64_4d_array( &
    fortran_array=z_dexner_dz_c, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (2 /= size(z_dexner_dz_c, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 'z_dexner_dz_c'"//char(10), &
        "    - actual = (", &
        size(z_dexner_dz_c, dim=1), ",", &
        size(z_dexner_dz_c, dim=2), ",", &
        size(z_dexner_dz_c, dim=3), ",", &
        size(z_dexner_dz_c, dim=4), &
        "), config propagated = (2, tmp_struct_symbol_18, tmp_struct_symbol_19, tmp_struct_symbol_20)"
    end if
#endif
#ifndef _OPENACC
    copy_or_ptr_z_dwdz_dd = copy_in_float64_3d_array( &
    fortran_array=z_dwdz_dd, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_dwdz_dd = copy_in_float64_3d_array( &
    fortran_array=z_dwdz_dd, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_exner_ex_pr = copy_in_float64_3d_array( &
    fortran_array=z_exner_ex_pr, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_exner_ex_pr = copy_in_float64_3d_array( &
    fortran_array=z_exner_ex_pr, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_exner_expl = copy_in_float64_2d_array( &
    fortran_array=z_exner_expl, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_exner_expl = copy_in_float64_2d_array( &
    fortran_array=z_exner_expl, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_exner_ic = copy_in_float64_2d_array( &
    fortran_array=z_exner_ic, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_exner_ic = copy_in_float64_2d_array( &
    fortran_array=z_exner_ic, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_flxdiv_mass = copy_in_float64_2d_array( &
    fortran_array=z_flxdiv_mass, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_flxdiv_mass = copy_in_float64_2d_array( &
    fortran_array=z_flxdiv_mass, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_flxdiv_theta = copy_in_float64_2d_array( &
    fortran_array=z_flxdiv_theta, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_flxdiv_theta = copy_in_float64_2d_array( &
    fortran_array=z_flxdiv_theta, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_grad_rth = copy_in_float64_4d_array( &
    fortran_array=z_grad_rth, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_grad_rth = copy_in_float64_4d_array( &
    fortran_array=z_grad_rth, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (4 /= size(z_grad_rth, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 'z_grad_rth'"//char(10), &
        "    - actual = (", &
        size(z_grad_rth, dim=1), ",", &
        size(z_grad_rth, dim=2), ",", &
        size(z_grad_rth, dim=3), ",", &
        size(z_grad_rth, dim=4), &
        "), config propagated = (4, tmp_struct_symbol_36, tmp_struct_symbol_37, tmp_struct_symbol_38)"
    end if
#endif
#ifndef _OPENACC
    copy_or_ptr_z_graddiv2_vn = copy_in_float64_2d_array( &
    fortran_array=z_graddiv2_vn, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_graddiv2_vn = copy_in_float64_2d_array( &
    fortran_array=z_graddiv2_vn, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_graddiv_vn = copy_in_float64_3d_array( &
    fortran_array=z_graddiv_vn, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_graddiv_vn = copy_in_float64_3d_array( &
    fortran_array=z_graddiv_vn, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_gradh_exner = copy_in_float64_3d_array( &
    fortran_array=z_gradh_exner, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_gradh_exner = copy_in_float64_3d_array( &
    fortran_array=z_gradh_exner, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_hydro_corr = copy_in_float64_2d_array( &
    fortran_array=z_hydro_corr, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_hydro_corr = copy_in_float64_2d_array( &
    fortran_array=z_hydro_corr, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
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
    copy_or_ptr_z_mflx_top = copy_in_float64_2d_array( &
    fortran_array=z_mflx_top, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_mflx_top = copy_in_float64_2d_array( &
    fortran_array=z_mflx_top, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_q = copy_in_float64_2d_array( &
    fortran_array=z_q, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_q = copy_in_float64_2d_array( &
    fortran_array=z_q, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_raylfac = copy_in_float64_1d_array( &
    fortran_array=z_raylfac, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_raylfac = copy_in_float64_1d_array( &
    fortran_array=z_raylfac, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_rho_e = copy_in_float64_3d_array( &
    fortran_array=z_rho_e, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_rho_e = copy_in_float64_3d_array( &
    fortran_array=z_rho_e, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_rho_expl = copy_in_float64_2d_array( &
    fortran_array=z_rho_expl, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_rho_expl = copy_in_float64_2d_array( &
    fortran_array=z_rho_expl, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_rho_v = copy_in_float64_3d_array( &
    fortran_array=z_rho_v, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_rho_v = copy_in_float64_3d_array( &
    fortran_array=z_rho_v, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_rth_pr = copy_in_float64_4d_array( &
    fortran_array=z_rth_pr, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_rth_pr = copy_in_float64_4d_array( &
    fortran_array=z_rth_pr, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (2 /= size(z_rth_pr, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 'z_rth_pr'"//char(10), &
        "    - actual = (", &
        size(z_rth_pr, dim=1), ",", &
        size(z_rth_pr, dim=2), ",", &
        size(z_rth_pr, dim=3), ",", &
        size(z_rth_pr, dim=4), &
        "), config propagated = (2, tmp_struct_symbol_33, tmp_struct_symbol_34, tmp_struct_symbol_35)"
    end if
#endif
#ifndef _OPENACC
    copy_or_ptr_z_th_ddz_exner_c = copy_in_float64_3d_array( &
    fortran_array=z_th_ddz_exner_c, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_th_ddz_exner_c = copy_in_float64_3d_array( &
    fortran_array=z_th_ddz_exner_c, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_theta_v_e = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_e, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_theta_v_e = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_e, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_theta_v_fl_e = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_fl_e, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_theta_v_fl_e = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_fl_e, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_theta_v_pr_ic = copy_in_float64_2d_array( &
    fortran_array=z_theta_v_pr_ic, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_theta_v_pr_ic = copy_in_float64_2d_array( &
    fortran_array=z_theta_v_pr_ic, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_theta_v_v = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_v, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_theta_v_v = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_v, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_vn_avg = copy_in_float64_2d_array( &
    fortran_array=z_vn_avg, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_vn_avg = copy_in_float64_2d_array( &
    fortran_array=z_vn_avg, &
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
    copy_or_ptr_z_w_concorr_mc = copy_in_float64_2d_array( &
    fortran_array=z_w_concorr_mc, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_w_concorr_mc = copy_in_float64_2d_array( &
    fortran_array=z_w_concorr_mc, &
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
#ifndef _OPENACC
    copy_or_ptr_z_w_expl = copy_in_float64_2d_array( &
    fortran_array=z_w_expl, &
    steal_arrays=.true., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_w_expl = copy_in_float64_2d_array( &
    fortran_array=z_w_expl, &
    steal_arrays=.true., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif


    if (is_initialized .eqv. .false.) then
      is_initialized = .true.
      dace_state = dace_init_solve_nh_predictor_pre( &
      bdy_divdamp = copy_or_ptr_bdy_divdamp, &
      enh_divdamp_fac = copy_or_ptr_enh_divdamp_fac, &
      global_data = copy_or_ptr_global_data, &
      p_int = copy_or_ptr_p_int, &
      p_nh = copy_or_ptr_p_nh, &
      p_nh_prog_nnew = copy_or_ptr_p_nh_prog_nnew, &
      p_nh_prog_nnow = copy_or_ptr_p_nh_prog_nnow, &
      p_patch = copy_or_ptr_p_patch, &
      prep_adv = copy_or_ptr_prep_adv, &
      scal_divdamp = copy_or_ptr_scal_divdamp, &
      z_alpha = copy_or_ptr_z_alpha, &
      z_beta = copy_or_ptr_z_beta, &
      z_contr_w_fl_l = copy_or_ptr_z_contr_w_fl_l, &
      z_dexner_dz_c = copy_or_ptr_z_dexner_dz_c, &
      z_dwdz_dd = copy_or_ptr_z_dwdz_dd, &
      z_exner_ex_pr = copy_or_ptr_z_exner_ex_pr, &
      z_exner_expl = copy_or_ptr_z_exner_expl, &
      z_exner_ic = copy_or_ptr_z_exner_ic, &
      z_flxdiv_mass = copy_or_ptr_z_flxdiv_mass, &
      z_flxdiv_theta = copy_or_ptr_z_flxdiv_theta, &
      z_grad_rth = copy_or_ptr_z_grad_rth, &
      z_graddiv2_vn = copy_or_ptr_z_graddiv2_vn, &
      z_graddiv_vn = copy_or_ptr_z_graddiv_vn, &
      z_gradh_exner = copy_or_ptr_z_gradh_exner, &
      z_hydro_corr = copy_or_ptr_z_hydro_corr, &
      z_kin_hor_e = copy_or_ptr_z_kin_hor_e, &
      z_mflx_top = copy_or_ptr_z_mflx_top, &
      z_q = copy_or_ptr_z_q, &
      z_raylfac = copy_or_ptr_z_raylfac, &
      z_rho_e = copy_or_ptr_z_rho_e, &
      z_rho_expl = copy_or_ptr_z_rho_expl, &
      z_rho_v = copy_or_ptr_z_rho_v, &
      z_rth_pr = copy_or_ptr_z_rth_pr, &
      z_th_ddz_exner_c = copy_or_ptr_z_th_ddz_exner_c, &
      z_theta_v_e = copy_or_ptr_z_theta_v_e, &
      z_theta_v_fl_e = copy_or_ptr_z_theta_v_fl_e, &
      z_theta_v_pr_ic = copy_or_ptr_z_theta_v_pr_ic, &
      z_theta_v_v = copy_or_ptr_z_theta_v_v, &
      z_vn_avg = copy_or_ptr_z_vn_avg, &
      z_vt_ie = copy_or_ptr_z_vt_ie, &
      z_w_concorr_mc = copy_or_ptr_z_w_concorr_mc, &
      z_w_concorr_me = copy_or_ptr_z_w_concorr_me, &
      z_w_expl = copy_or_ptr_z_w_expl, &
      f2dace_OPTIONAL_lacc = f2dace_OPTIONAL_lacc, &
      alin = alin, &
      aqdr = aqdr, &
      bqdr = bqdr, &
      df32 = df32, &
      df42 = df42, &
      distv_bary_1 = distv_bary_1, &
      distv_bary_2 = distv_bary_2, &
      dt_linintp_ubc = dt_linintp_ubc, &
      dt_linintp_ubc_nnew = dt_linintp_ubc_nnew, &
      dt_linintp_ubc_nnow = dt_linintp_ubc_nnow, &
      dt_shift = dt_shift, &
      dthalf = dthalf, &
      dtime = dtime, &
      dz32 = dz32, &
      dz42 = dz42, &
      dzlin = dzlin, &
      dzqdr = dzqdr, &
      i_endblk = i_endblk, &
      i_endidx = i_endidx, &
      i_startblk = i_startblk, &
      i_startidx = i_startidx, &
      idyn_timestep = idyn_timestep, &
      ishift = ishift, &
      istep = istep, &
      jb = jb, &
      jc = jc, &
      je = je, &
      jg = jg, &
      jk = jk, &
      jk_start = jk_start, &
      jks = jks, &
      jstep = jstep, &
      l_child_vertnest = l_child_vertnest, &
      l_init = l_init, &
      l_recompute = l_recompute, &
      l_vert_nested = l_vert_nested, &
      lacc = DACE_OPT_PROXY_lacc, &
      lclean_mflx = lclean_mflx, &
      lprep_adv = lprep_adv, &
      lsave_mflx = lsave_mflx, &
      lvn_only = lvn_only, &
      lvn_pos = lvn_pos, &
      nblks_gradp = nblks_gradp, &
      nlen_gradp = nlen_gradp, &
      nlev = nlev, &
      nlevp1 = nlevp1, &
      nnew = nnew, &
      nnow = nnow, &
      nproma_gradp = nproma_gradp, &
      npromz_gradp = npromz_gradp, &
      nshift = nshift, &
      nshift_total = nshift_total, &
      ntl1 = ntl1, &
      ntl2 = ntl2, &
      nvar = nvar, &
      r_dtimensubsteps = r_dtimensubsteps, &
      r_nsubsteps = r_nsubsteps, &
      rl_end = rl_end, &
      rl_start = rl_start, &
      scal_divdamp_o2 = scal_divdamp_o2, &
      wgt_nnew_rth = wgt_nnew_rth, &
      wgt_nnew_vel = wgt_nnew_vel, &
      wgt_nnow_rth = wgt_nnow_rth, &
      wgt_nnow_vel = wgt_nnow_vel, &
      z_a = z_a, &
      z_b = z_b, &
      z_c = z_c, &
      z_d_vn_dmp = z_d_vn_dmp, &
      z_d_vn_iau = z_d_vn_iau, &
      z_ddt_vn_apc = z_ddt_vn_apc, &
      z_ddt_vn_cor = z_ddt_vn_cor, &
      z_ddt_vn_dyn = z_ddt_vn_dyn, &
      z_ddt_vn_pgr = z_ddt_vn_pgr, &
      z_ddt_vn_ray = z_ddt_vn_ray, &
      z_g = z_g, &
      z_gamma = z_gamma, &
      z_ntdistv_bary_1 = z_ntdistv_bary_1, &
      z_ntdistv_bary_2 = z_ntdistv_bary_2, &
      z_rho_tavg = z_rho_tavg, &
      z_rho_tavg_m1 = z_rho_tavg_m1, &
      z_theta1 = z_theta1, &
      z_theta2 = z_theta2, &
      z_theta_tavg = z_theta_tavg, &
      z_theta_tavg_m1 = z_theta_tavg_m1, &
      z_theta_v_pr_mc = z_theta_v_pr_mc, &
      z_theta_v_pr_mc_m1 = z_theta_v_pr_mc_m1, &
      z_w_backtraj = z_w_backtraj, &
      zf = zf &
    )
    end if

! SOLVE_NH PART TIMERS : PRATYAI
call cpu_time(t0)

    call dace_program_solve_nh_predictor_pre( &
      state = dace_state, &
      bdy_divdamp = copy_or_ptr_bdy_divdamp, &
      enh_divdamp_fac = copy_or_ptr_enh_divdamp_fac, &
      global_data = copy_or_ptr_global_data, &
      p_int = copy_or_ptr_p_int, &
      p_nh = copy_or_ptr_p_nh, &
      p_nh_prog_nnew = copy_or_ptr_p_nh_prog_nnew, &
      p_nh_prog_nnow = copy_or_ptr_p_nh_prog_nnow, &
      p_patch = copy_or_ptr_p_patch, &
      prep_adv = copy_or_ptr_prep_adv, &
      scal_divdamp = copy_or_ptr_scal_divdamp, &
      z_alpha = copy_or_ptr_z_alpha, &
      z_beta = copy_or_ptr_z_beta, &
      z_contr_w_fl_l = copy_or_ptr_z_contr_w_fl_l, &
      z_dexner_dz_c = copy_or_ptr_z_dexner_dz_c, &
      z_dwdz_dd = copy_or_ptr_z_dwdz_dd, &
      z_exner_ex_pr = copy_or_ptr_z_exner_ex_pr, &
      z_exner_expl = copy_or_ptr_z_exner_expl, &
      z_exner_ic = copy_or_ptr_z_exner_ic, &
      z_flxdiv_mass = copy_or_ptr_z_flxdiv_mass, &
      z_flxdiv_theta = copy_or_ptr_z_flxdiv_theta, &
      z_grad_rth = copy_or_ptr_z_grad_rth, &
      z_graddiv2_vn = copy_or_ptr_z_graddiv2_vn, &
      z_graddiv_vn = copy_or_ptr_z_graddiv_vn, &
      z_gradh_exner = copy_or_ptr_z_gradh_exner, &
      z_hydro_corr = copy_or_ptr_z_hydro_corr, &
      z_kin_hor_e = copy_or_ptr_z_kin_hor_e, &
      z_mflx_top = copy_or_ptr_z_mflx_top, &
      z_q = copy_or_ptr_z_q, &
      z_raylfac = copy_or_ptr_z_raylfac, &
      z_rho_e = copy_or_ptr_z_rho_e, &
      z_rho_expl = copy_or_ptr_z_rho_expl, &
      z_rho_v = copy_or_ptr_z_rho_v, &
      z_rth_pr = copy_or_ptr_z_rth_pr, &
      z_th_ddz_exner_c = copy_or_ptr_z_th_ddz_exner_c, &
      z_theta_v_e = copy_or_ptr_z_theta_v_e, &
      z_theta_v_fl_e = copy_or_ptr_z_theta_v_fl_e, &
      z_theta_v_pr_ic = copy_or_ptr_z_theta_v_pr_ic, &
      z_theta_v_v = copy_or_ptr_z_theta_v_v, &
      z_vn_avg = copy_or_ptr_z_vn_avg, &
      z_vt_ie = copy_or_ptr_z_vt_ie, &
      z_w_concorr_mc = copy_or_ptr_z_w_concorr_mc, &
      z_w_concorr_me = copy_or_ptr_z_w_concorr_me, &
      z_w_expl = copy_or_ptr_z_w_expl, &
      f2dace_OPTIONAL_lacc = f2dace_OPTIONAL_lacc, &
      alin = alin, &
      aqdr = aqdr, &
      bqdr = bqdr, &
      df32 = df32, &
      df42 = df42, &
      distv_bary_1 = distv_bary_1, &
      distv_bary_2 = distv_bary_2, &
      dt_linintp_ubc = dt_linintp_ubc, &
      dt_linintp_ubc_nnew = dt_linintp_ubc_nnew, &
      dt_linintp_ubc_nnow = dt_linintp_ubc_nnow, &
      dt_shift = dt_shift, &
      dthalf = dthalf, &
      dtime = dtime, &
      dz32 = dz32, &
      dz42 = dz42, &
      dzlin = dzlin, &
      dzqdr = dzqdr, &
      i_endblk = i_endblk, &
      i_endidx = i_endidx, &
      i_startblk = i_startblk, &
      i_startidx = i_startidx, &
      idyn_timestep = idyn_timestep, &
      ishift = ishift, &
      istep = istep, &
      jb = jb, &
      jc = jc, &
      je = je, &
      jg = jg, &
      jk = jk, &
      jk_start = jk_start, &
      jks = jks, &
      jstep = jstep, &
      l_child_vertnest = l_child_vertnest, &
      l_init = l_init, &
      l_recompute = l_recompute, &
      l_vert_nested = l_vert_nested, &
      lacc = DACE_OPT_PROXY_lacc, &
      lclean_mflx = lclean_mflx, &
      lprep_adv = lprep_adv, &
      lsave_mflx = lsave_mflx, &
      lvn_only = lvn_only, &
      lvn_pos = lvn_pos, &
      nblks_gradp = nblks_gradp, &
      nlen_gradp = nlen_gradp, &
      nlev = nlev, &
      nlevp1 = nlevp1, &
      nnew = nnew, &
      nnow = nnow, &
      nproma_gradp = nproma_gradp, &
      npromz_gradp = npromz_gradp, &
      nshift = nshift, &
      nshift_total = nshift_total, &
      ntl1 = ntl1, &
      ntl2 = ntl2, &
      nvar = nvar, &
      r_dtimensubsteps = r_dtimensubsteps, &
      r_nsubsteps = r_nsubsteps, &
      rl_end = rl_end, &
      rl_start = rl_start, &
      scal_divdamp_o2 = scal_divdamp_o2, &
      wgt_nnew_rth = wgt_nnew_rth, &
      wgt_nnew_vel = wgt_nnew_vel, &
      wgt_nnow_rth = wgt_nnow_rth, &
      wgt_nnow_vel = wgt_nnow_vel, &
      z_a = z_a, &
      z_b = z_b, &
      z_c = z_c, &
      z_d_vn_dmp = z_d_vn_dmp, &
      z_d_vn_iau = z_d_vn_iau, &
      z_ddt_vn_apc = z_ddt_vn_apc, &
      z_ddt_vn_cor = z_ddt_vn_cor, &
      z_ddt_vn_dyn = z_ddt_vn_dyn, &
      z_ddt_vn_pgr = z_ddt_vn_pgr, &
      z_ddt_vn_ray = z_ddt_vn_ray, &
      z_g = z_g, &
      z_gamma = z_gamma, &
      z_ntdistv_bary_1 = z_ntdistv_bary_1, &
      z_ntdistv_bary_2 = z_ntdistv_bary_2, &
      z_rho_tavg = z_rho_tavg, &
      z_rho_tavg_m1 = z_rho_tavg_m1, &
      z_theta1 = z_theta1, &
      z_theta2 = z_theta2, &
      z_theta_tavg = z_theta_tavg, &
      z_theta_tavg_m1 = z_theta_tavg_m1, &
      z_theta_v_pr_mc = z_theta_v_pr_mc, &
      z_theta_v_pr_mc_m1 = z_theta_v_pr_mc_m1, &
      z_w_backtraj = z_w_backtraj, &
      zf = zf &
    )

! SOLVE_NH PART TIMERS : PRATYAI
call cpu_time(t1)
print *, 'PREPRE INTERNAL (s): ', (t1-t0)

    call copy_back_global_data_type(copy_or_ptr_global_data)
    call copy_back_t_int_state(p_int, copy_or_ptr_p_int)
    call copy_back_t_nh_state(p_nh, copy_or_ptr_p_nh)
    call copy_back_t_nh_prog(p_nh_prog_nnew, copy_or_ptr_p_nh_prog_nnew)
    call copy_back_t_nh_prog(p_nh_prog_nnow, copy_or_ptr_p_nh_prog_nnow)
    call copy_back_t_patch(p_patch, copy_or_ptr_p_patch)
    call copy_back_t_prepare_adv(prep_adv, copy_or_ptr_prep_adv)


  end subroutine run_solve_nh_predictor_pre

  subroutine run_solve_nh_predictor_pre_verification( &
    bdy_divdamp, &
    enh_divdamp_fac, &
    p_int, &
    p_nh, &
    p_nh_prog_nnew, &
    p_nh_prog_nnow, &
    p_patch, &
    prep_adv, &
    scal_divdamp, &
    z_alpha, &
    z_beta, &
    z_contr_w_fl_l, &
    z_dexner_dz_c, &
    z_dwdz_dd, &
    z_exner_ex_pr, &
    z_exner_expl, &
    z_exner_ic, &
    z_flxdiv_mass, &
    z_flxdiv_theta, &
    z_grad_rth, &
    z_graddiv2_vn, &
    z_graddiv_vn, &
    z_gradh_exner, &
    z_hydro_corr, &
    z_kin_hor_e, &
    z_mflx_top, &
    z_q, &
    z_raylfac, &
    z_rho_e, &
    z_rho_expl, &
    z_rho_v, &
    z_rth_pr, &
    z_th_ddz_exner_c, &
    z_theta_v_e, &
    z_theta_v_fl_e, &
    z_theta_v_pr_ic, &
    z_theta_v_v, &
    z_vn_avg, &
    z_vt_ie, &
    z_w_concorr_mc, &
    z_w_concorr_me, &
    z_w_expl, &
    alin, &
    aqdr, &
    bqdr, &
    df32, &
    df42, &
    distv_bary_1, &
    distv_bary_2, &
    dt_linintp_ubc, &
    dt_linintp_ubc_nnew, &
    dt_linintp_ubc_nnow, &
    dt_shift, &
    dthalf, &
    dtime, &
    dz32, &
    dz42, &
    dzlin, &
    dzqdr, &
    i_endblk, &
    i_endidx, &
    i_startblk, &
    i_startidx, &
    idyn_timestep, &
    ishift, &
    istep, &
    jb, &
    jc, &
    je, &
    jg, &
    jk, &
    jk_start, &
    jks, &
    jstep, &
    l_child_vertnest, &
    l_init, &
    l_recompute, &
    l_vert_nested, &
    lacc, &
    lclean_mflx, &
    lprep_adv, &
    lsave_mflx, &
    lvn_only, &
    lvn_pos, &
    nblks_gradp, &
    nlen_gradp, &
    nlev, &
    nlevp1, &
    nnew, &
    nnow, &
    nproma_gradp, &
    npromz_gradp, &
    nshift, &
    nshift_total, &
    ntl1, &
    ntl2, &
    nvar, &
    r_dtimensubsteps, &
    r_nsubsteps, &
    rl_end, &
    rl_start, &
    scal_divdamp_o2, &
    wgt_nnew_rth, &
    wgt_nnew_vel, &
    wgt_nnow_rth, &
    wgt_nnow_vel, &
    z_a, &
    z_b, &
    z_c, &
    z_d_vn_dmp, &
    z_d_vn_iau, &
    z_ddt_vn_apc, &
    z_ddt_vn_cor, &
    z_ddt_vn_dyn, &
    z_ddt_vn_pgr, &
    z_ddt_vn_ray, &
    z_g, &
    z_gamma, &
    z_ntdistv_bary_1, &
    z_ntdistv_bary_2, &
    z_rho_tavg, &
    z_rho_tavg_m1, &
    z_theta1, &
    z_theta2, &
    z_theta_tavg, &
    z_theta_tavg_m1, &
    z_theta_v_pr_mc, &
    z_theta_v_pr_mc_m1, &
    z_w_backtraj, &
    zf &
  )
    real(kind=c_double), dimension(:), target :: bdy_divdamp
    real(kind=c_double), dimension(:), target :: enh_divdamp_fac
    type(t_int_state), target :: p_int
    type(t_nh_state), target :: p_nh
    type(t_nh_prog), target :: p_nh_prog_nnew
    type(t_nh_prog), target :: p_nh_prog_nnow
    type(t_patch), target :: p_patch
    type(t_prepare_adv), target :: prep_adv
    real(kind=c_double), dimension(:), target :: scal_divdamp
    real(kind=c_double), dimension(:,:), target :: z_alpha
    real(kind=c_double), dimension(:,:), target :: z_beta
    real(kind=c_double), dimension(:,:), target :: z_contr_w_fl_l
    real(kind=c_double), dimension(:,:,:,:), target :: z_dexner_dz_c
    real(kind=c_double), dimension(:,:,:), target :: z_dwdz_dd
    real(kind=c_double), dimension(:,:,:), target :: z_exner_ex_pr
    real(kind=c_double), dimension(:,:), target :: z_exner_expl
    real(kind=c_double), dimension(:,:), target :: z_exner_ic
    real(kind=c_double), dimension(:,:), target :: z_flxdiv_mass
    real(kind=c_double), dimension(:,:), target :: z_flxdiv_theta
    real(kind=c_double), dimension(:,:,:,:), target :: z_grad_rth
    real(kind=c_double), dimension(:,:), target :: z_graddiv2_vn
    real(kind=c_double), dimension(:,:,:), target :: z_graddiv_vn
    real(kind=c_double), dimension(:,:,:), target :: z_gradh_exner
    real(kind=c_double), dimension(:,:), target :: z_hydro_corr
    real(kind=c_double), dimension(:,:,:), target :: z_kin_hor_e
    real(kind=c_double), dimension(:,:), target :: z_mflx_top
    real(kind=c_double), dimension(:,:), target :: z_q
    real(kind=c_double), dimension(:), target :: z_raylfac
    real(kind=c_double), dimension(:,:,:), target :: z_rho_e
    real(kind=c_double), dimension(:,:), target :: z_rho_expl
    real(kind=c_double), dimension(:,:,:), target :: z_rho_v
    real(kind=c_double), dimension(:,:,:,:), target :: z_rth_pr
    real(kind=c_double), dimension(:,:,:), target :: z_th_ddz_exner_c
    real(kind=c_double), dimension(:,:,:), target :: z_theta_v_e
    real(kind=c_double), dimension(:,:,:), target :: z_theta_v_fl_e
    real(kind=c_double), dimension(:,:), target :: z_theta_v_pr_ic
    real(kind=c_double), dimension(:,:,:), target :: z_theta_v_v
    real(kind=c_double), dimension(:,:), target :: z_vn_avg
    real(kind=c_double), dimension(:,:,:), target :: z_vt_ie
    real(kind=c_double), dimension(:,:), target :: z_w_concorr_mc
    real(kind=c_double), dimension(:,:,:), target :: z_w_concorr_me
    real(kind=c_double), dimension(:,:), target :: z_w_expl
    real(kind=c_double) :: alin
    real(kind=c_double) :: aqdr
    real(kind=c_double) :: bqdr
    real(kind=c_double) :: df32
    real(kind=c_double) :: df42
    real(kind=c_double) :: distv_bary_1
    real(kind=c_double) :: distv_bary_2
    real(kind=c_double) :: dt_linintp_ubc
    real(kind=c_double) :: dt_linintp_ubc_nnew
    real(kind=c_double) :: dt_linintp_ubc_nnow
    real(kind=c_double) :: dt_shift
    real(kind=c_double) :: dthalf
    real(kind=c_double) :: dtime
    real(kind=c_double) :: dz32
    real(kind=c_double) :: dz42
    real(kind=c_double) :: dzlin
    real(kind=c_double) :: dzqdr
    integer(kind=c_int) :: i_endblk
    integer(kind=c_int) :: i_endidx
    integer(kind=c_int) :: i_startblk
    integer(kind=c_int) :: i_startidx
    integer(kind=c_int) :: idyn_timestep
    integer(kind=c_int) :: ishift
    integer(kind=c_int) :: istep
    integer(kind=c_int) :: jb
    integer(kind=c_int) :: jc
    integer(kind=c_int) :: je
    integer(kind=c_int) :: jg
    integer(kind=c_int) :: jk
    integer(kind=c_int) :: jk_start
    integer(kind=c_int) :: jks
    integer(kind=c_int) :: jstep
    integer(kind=c_int) :: l_child_vertnest
    integer(kind=c_int) :: l_init
    integer(kind=c_int) :: l_recompute
    integer(kind=c_int) :: l_vert_nested
    integer(kind=c_int), optional :: lacc
    integer(kind=c_int) :: lclean_mflx
    integer(kind=c_int) :: lprep_adv
    integer(kind=c_int) :: lsave_mflx
    integer(kind=c_int) :: lvn_only
    integer(kind=c_int) :: lvn_pos
    integer(kind=c_int) :: nblks_gradp
    integer(kind=c_int) :: nlen_gradp
    integer(kind=c_int) :: nlev
    integer(kind=c_int) :: nlevp1
    integer(kind=c_int) :: nnew
    integer(kind=c_int) :: nnow
    integer(kind=c_int) :: nproma_gradp
    integer(kind=c_int) :: npromz_gradp
    integer(kind=c_int) :: nshift
    integer(kind=c_int) :: nshift_total
    integer(kind=c_int) :: ntl1
    integer(kind=c_int) :: ntl2
    integer(kind=c_int) :: nvar
    real(kind=c_double) :: r_dtimensubsteps
    real(kind=c_double) :: r_nsubsteps
    integer(kind=c_int) :: rl_end
    integer(kind=c_int) :: rl_start
    real(kind=c_double) :: scal_divdamp_o2
    real(kind=c_double) :: wgt_nnew_rth
    real(kind=c_double) :: wgt_nnew_vel
    real(kind=c_double) :: wgt_nnow_rth
    real(kind=c_double) :: wgt_nnow_vel
    real(kind=c_double) :: z_a
    real(kind=c_double) :: z_b
    real(kind=c_double) :: z_c
    real(kind=c_double) :: z_d_vn_dmp
    real(kind=c_double) :: z_d_vn_iau
    real(kind=c_double) :: z_ddt_vn_apc
    real(kind=c_double) :: z_ddt_vn_cor
    real(kind=c_double) :: z_ddt_vn_dyn
    real(kind=c_double) :: z_ddt_vn_pgr
    real(kind=c_double) :: z_ddt_vn_ray
    real(kind=c_double) :: z_g
    real(kind=c_double) :: z_gamma
    real(kind=c_double) :: z_ntdistv_bary_1
    real(kind=c_double) :: z_ntdistv_bary_2
    real(kind=c_double) :: z_rho_tavg
    real(kind=c_double) :: z_rho_tavg_m1
    real(kind=c_double) :: z_theta1
    real(kind=c_double) :: z_theta2
    real(kind=c_double) :: z_theta_tavg
    real(kind=c_double) :: z_theta_tavg_m1
    real(kind=c_double) :: z_theta_v_pr_mc
    real(kind=c_double) :: z_theta_v_pr_mc_m1
    real(kind=c_double) :: z_w_backtraj
    real(kind=c_double) :: zf
    !Optional helper parameter
    !WARNING: HACKFIX, POSSIBLE EXPLOSION
    integer(kind=c_int) :: f2dace_OPTIONAL_lacc

    integer(kind=c_int) :: DACE_OPT_PROXY_lacc


    !$ACC WAIT


    if (present(lacc)) then
      f2dace_OPTIONAL_lacc = 1
      DACE_OPT_PROXY_lacc = lacc
      
    else
      f2dace_OPTIONAL_lacc = 0
    end if


    call check_initializations()

#ifndef _OPENACC
    copy_or_ptr_bdy_divdamp = copy_in_float64_1d_array( &
    fortran_array=bdy_divdamp, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_bdy_divdamp = copy_in_float64_1d_array( &
    fortran_array=bdy_divdamp, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
    copy_or_ptr_enh_divdamp_fac = copy_in_float64_1d_array( &
    fortran_array=enh_divdamp_fac, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_global_data = copy_in_global_data_type( &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_int = copy_in_t_int_state( &
    fortran_obj=p_int, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_nh = copy_in_t_nh_state( &
    fortran_obj=p_nh, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_nh_prog_nnew = copy_in_t_nh_prog( &
    fortran_obj=p_nh_prog_nnew, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_nh_prog_nnow = copy_in_t_nh_prog( &
    fortran_obj=p_nh_prog_nnow, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_p_patch = copy_in_t_patch( &
    fortran_obj=p_patch, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
    copy_or_ptr_prep_adv = copy_in_t_prepare_adv( &
    fortran_obj=prep_adv, &
    steal_arrays=.false., &
    minimal_structs=.false. &
  )
#ifndef _OPENACC
    copy_or_ptr_scal_divdamp = copy_in_float64_1d_array( &
    fortran_array=scal_divdamp, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_scal_divdamp = copy_in_float64_1d_array( &
    fortran_array=scal_divdamp, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_alpha = copy_in_float64_2d_array( &
    fortran_array=z_alpha, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_alpha = copy_in_float64_2d_array( &
    fortran_array=z_alpha, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_beta = copy_in_float64_2d_array( &
    fortran_array=z_beta, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_beta = copy_in_float64_2d_array( &
    fortran_array=z_beta, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_contr_w_fl_l = copy_in_float64_2d_array( &
    fortran_array=z_contr_w_fl_l, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_contr_w_fl_l = copy_in_float64_2d_array( &
    fortran_array=z_contr_w_fl_l, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_dexner_dz_c = copy_in_float64_4d_array( &
    fortran_array=z_dexner_dz_c, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_dexner_dz_c = copy_in_float64_4d_array( &
    fortran_array=z_dexner_dz_c, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (2 /= size(z_dexner_dz_c, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 'z_dexner_dz_c'"//char(10), &
        "    - actual = (", &
        size(z_dexner_dz_c, dim=1), ",", &
        size(z_dexner_dz_c, dim=2), ",", &
        size(z_dexner_dz_c, dim=3), ",", &
        size(z_dexner_dz_c, dim=4), &
        "), config propagated = (2, tmp_struct_symbol_18, tmp_struct_symbol_19, tmp_struct_symbol_20)"
    end if
#endif
#ifndef _OPENACC
    copy_or_ptr_z_dwdz_dd = copy_in_float64_3d_array( &
    fortran_array=z_dwdz_dd, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_dwdz_dd = copy_in_float64_3d_array( &
    fortran_array=z_dwdz_dd, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_exner_ex_pr = copy_in_float64_3d_array( &
    fortran_array=z_exner_ex_pr, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_exner_ex_pr = copy_in_float64_3d_array( &
    fortran_array=z_exner_ex_pr, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_exner_expl = copy_in_float64_2d_array( &
    fortran_array=z_exner_expl, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_exner_expl = copy_in_float64_2d_array( &
    fortran_array=z_exner_expl, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_exner_ic = copy_in_float64_2d_array( &
    fortran_array=z_exner_ic, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_exner_ic = copy_in_float64_2d_array( &
    fortran_array=z_exner_ic, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_flxdiv_mass = copy_in_float64_2d_array( &
    fortran_array=z_flxdiv_mass, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_flxdiv_mass = copy_in_float64_2d_array( &
    fortran_array=z_flxdiv_mass, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_flxdiv_theta = copy_in_float64_2d_array( &
    fortran_array=z_flxdiv_theta, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_flxdiv_theta = copy_in_float64_2d_array( &
    fortran_array=z_flxdiv_theta, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_grad_rth = copy_in_float64_4d_array( &
    fortran_array=z_grad_rth, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_grad_rth = copy_in_float64_4d_array( &
    fortran_array=z_grad_rth, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (4 /= size(z_grad_rth, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 'z_grad_rth'"//char(10), &
        "    - actual = (", &
        size(z_grad_rth, dim=1), ",", &
        size(z_grad_rth, dim=2), ",", &
        size(z_grad_rth, dim=3), ",", &
        size(z_grad_rth, dim=4), &
        "), config propagated = (4, tmp_struct_symbol_36, tmp_struct_symbol_37, tmp_struct_symbol_38)"
    end if
#endif
#ifndef _OPENACC
    copy_or_ptr_z_graddiv2_vn = copy_in_float64_2d_array( &
    fortran_array=z_graddiv2_vn, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_graddiv2_vn = copy_in_float64_2d_array( &
    fortran_array=z_graddiv2_vn, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_graddiv_vn = copy_in_float64_3d_array( &
    fortran_array=z_graddiv_vn, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_graddiv_vn = copy_in_float64_3d_array( &
    fortran_array=z_graddiv_vn, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_gradh_exner = copy_in_float64_3d_array( &
    fortran_array=z_gradh_exner, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_gradh_exner = copy_in_float64_3d_array( &
    fortran_array=z_gradh_exner, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_hydro_corr = copy_in_float64_2d_array( &
    fortran_array=z_hydro_corr, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_hydro_corr = copy_in_float64_2d_array( &
    fortran_array=z_hydro_corr, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
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
    copy_or_ptr_z_mflx_top = copy_in_float64_2d_array( &
    fortran_array=z_mflx_top, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_mflx_top = copy_in_float64_2d_array( &
    fortran_array=z_mflx_top, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_q = copy_in_float64_2d_array( &
    fortran_array=z_q, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_q = copy_in_float64_2d_array( &
    fortran_array=z_q, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_raylfac = copy_in_float64_1d_array( &
    fortran_array=z_raylfac, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_raylfac = copy_in_float64_1d_array( &
    fortran_array=z_raylfac, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_rho_e = copy_in_float64_3d_array( &
    fortran_array=z_rho_e, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_rho_e = copy_in_float64_3d_array( &
    fortran_array=z_rho_e, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_rho_expl = copy_in_float64_2d_array( &
    fortran_array=z_rho_expl, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_rho_expl = copy_in_float64_2d_array( &
    fortran_array=z_rho_expl, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_rho_v = copy_in_float64_3d_array( &
    fortran_array=z_rho_v, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_rho_v = copy_in_float64_3d_array( &
    fortran_array=z_rho_v, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_rth_pr = copy_in_float64_4d_array( &
    fortran_array=z_rth_pr, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_rth_pr = copy_in_float64_4d_array( &
    fortran_array=z_rth_pr, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#if defined(DACE_SUBST_VERIFY)
    if (2 /= size(z_rth_pr, dim=1)) then
      print *, &
        "Array size conflicts with config propagation for array 'z_rth_pr'"//char(10), &
        "    - actual = (", &
        size(z_rth_pr, dim=1), ",", &
        size(z_rth_pr, dim=2), ",", &
        size(z_rth_pr, dim=3), ",", &
        size(z_rth_pr, dim=4), &
        "), config propagated = (2, tmp_struct_symbol_33, tmp_struct_symbol_34, tmp_struct_symbol_35)"
    end if
#endif
#ifndef _OPENACC
    copy_or_ptr_z_th_ddz_exner_c = copy_in_float64_3d_array( &
    fortran_array=z_th_ddz_exner_c, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_th_ddz_exner_c = copy_in_float64_3d_array( &
    fortran_array=z_th_ddz_exner_c, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_theta_v_e = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_e, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_theta_v_e = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_e, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_theta_v_fl_e = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_fl_e, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_theta_v_fl_e = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_fl_e, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_theta_v_pr_ic = copy_in_float64_2d_array( &
    fortran_array=z_theta_v_pr_ic, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_theta_v_pr_ic = copy_in_float64_2d_array( &
    fortran_array=z_theta_v_pr_ic, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_theta_v_v = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_v, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_theta_v_v = copy_in_float64_3d_array( &
    fortran_array=z_theta_v_v, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif
#ifndef _OPENACC
    copy_or_ptr_z_vn_avg = copy_in_float64_2d_array( &
    fortran_array=z_vn_avg, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_vn_avg = copy_in_float64_2d_array( &
    fortran_array=z_vn_avg, &
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
    copy_or_ptr_z_w_concorr_mc = copy_in_float64_2d_array( &
    fortran_array=z_w_concorr_mc, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_w_concorr_mc = copy_in_float64_2d_array( &
    fortran_array=z_w_concorr_mc, &
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
#ifndef _OPENACC
    copy_or_ptr_z_w_expl = copy_in_float64_2d_array( &
    fortran_array=z_w_expl, &
    steal_arrays=.false., &
    use_openacc=.false., &
    minimal_structs=.false. &
  )

#else
    copy_or_ptr_z_w_expl = copy_in_float64_2d_array( &
    fortran_array=z_w_expl, &
    steal_arrays=.false., &
    use_openacc=.true., &
    minimal_structs=.false. &
  )
#endif


    if (is_initialized .eqv. .false.) then
      is_initialized = .true.
      dace_state = dace_init_solve_nh_predictor_pre( &
      bdy_divdamp = copy_or_ptr_bdy_divdamp, &
      enh_divdamp_fac = copy_or_ptr_enh_divdamp_fac, &
      global_data = copy_or_ptr_global_data, &
      p_int = copy_or_ptr_p_int, &
      p_nh = copy_or_ptr_p_nh, &
      p_nh_prog_nnew = copy_or_ptr_p_nh_prog_nnew, &
      p_nh_prog_nnow = copy_or_ptr_p_nh_prog_nnow, &
      p_patch = copy_or_ptr_p_patch, &
      prep_adv = copy_or_ptr_prep_adv, &
      scal_divdamp = copy_or_ptr_scal_divdamp, &
      z_alpha = copy_or_ptr_z_alpha, &
      z_beta = copy_or_ptr_z_beta, &
      z_contr_w_fl_l = copy_or_ptr_z_contr_w_fl_l, &
      z_dexner_dz_c = copy_or_ptr_z_dexner_dz_c, &
      z_dwdz_dd = copy_or_ptr_z_dwdz_dd, &
      z_exner_ex_pr = copy_or_ptr_z_exner_ex_pr, &
      z_exner_expl = copy_or_ptr_z_exner_expl, &
      z_exner_ic = copy_or_ptr_z_exner_ic, &
      z_flxdiv_mass = copy_or_ptr_z_flxdiv_mass, &
      z_flxdiv_theta = copy_or_ptr_z_flxdiv_theta, &
      z_grad_rth = copy_or_ptr_z_grad_rth, &
      z_graddiv2_vn = copy_or_ptr_z_graddiv2_vn, &
      z_graddiv_vn = copy_or_ptr_z_graddiv_vn, &
      z_gradh_exner = copy_or_ptr_z_gradh_exner, &
      z_hydro_corr = copy_or_ptr_z_hydro_corr, &
      z_kin_hor_e = copy_or_ptr_z_kin_hor_e, &
      z_mflx_top = copy_or_ptr_z_mflx_top, &
      z_q = copy_or_ptr_z_q, &
      z_raylfac = copy_or_ptr_z_raylfac, &
      z_rho_e = copy_or_ptr_z_rho_e, &
      z_rho_expl = copy_or_ptr_z_rho_expl, &
      z_rho_v = copy_or_ptr_z_rho_v, &
      z_rth_pr = copy_or_ptr_z_rth_pr, &
      z_th_ddz_exner_c = copy_or_ptr_z_th_ddz_exner_c, &
      z_theta_v_e = copy_or_ptr_z_theta_v_e, &
      z_theta_v_fl_e = copy_or_ptr_z_theta_v_fl_e, &
      z_theta_v_pr_ic = copy_or_ptr_z_theta_v_pr_ic, &
      z_theta_v_v = copy_or_ptr_z_theta_v_v, &
      z_vn_avg = copy_or_ptr_z_vn_avg, &
      z_vt_ie = copy_or_ptr_z_vt_ie, &
      z_w_concorr_mc = copy_or_ptr_z_w_concorr_mc, &
      z_w_concorr_me = copy_or_ptr_z_w_concorr_me, &
      z_w_expl = copy_or_ptr_z_w_expl, &
      f2dace_OPTIONAL_lacc = f2dace_OPTIONAL_lacc, &
      alin = alin, &
      aqdr = aqdr, &
      bqdr = bqdr, &
      df32 = df32, &
      df42 = df42, &
      distv_bary_1 = distv_bary_1, &
      distv_bary_2 = distv_bary_2, &
      dt_linintp_ubc = dt_linintp_ubc, &
      dt_linintp_ubc_nnew = dt_linintp_ubc_nnew, &
      dt_linintp_ubc_nnow = dt_linintp_ubc_nnow, &
      dt_shift = dt_shift, &
      dthalf = dthalf, &
      dtime = dtime, &
      dz32 = dz32, &
      dz42 = dz42, &
      dzlin = dzlin, &
      dzqdr = dzqdr, &
      i_endblk = i_endblk, &
      i_endidx = i_endidx, &
      i_startblk = i_startblk, &
      i_startidx = i_startidx, &
      idyn_timestep = idyn_timestep, &
      ishift = ishift, &
      istep = istep, &
      jb = jb, &
      jc = jc, &
      je = je, &
      jg = jg, &
      jk = jk, &
      jk_start = jk_start, &
      jks = jks, &
      jstep = jstep, &
      l_child_vertnest = l_child_vertnest, &
      l_init = l_init, &
      l_recompute = l_recompute, &
      l_vert_nested = l_vert_nested, &
      lacc = DACE_OPT_PROXY_lacc, &
      lclean_mflx = lclean_mflx, &
      lprep_adv = lprep_adv, &
      lsave_mflx = lsave_mflx, &
      lvn_only = lvn_only, &
      lvn_pos = lvn_pos, &
      nblks_gradp = nblks_gradp, &
      nlen_gradp = nlen_gradp, &
      nlev = nlev, &
      nlevp1 = nlevp1, &
      nnew = nnew, &
      nnow = nnow, &
      nproma_gradp = nproma_gradp, &
      npromz_gradp = npromz_gradp, &
      nshift = nshift, &
      nshift_total = nshift_total, &
      ntl1 = ntl1, &
      ntl2 = ntl2, &
      nvar = nvar, &
      r_dtimensubsteps = r_dtimensubsteps, &
      r_nsubsteps = r_nsubsteps, &
      rl_end = rl_end, &
      rl_start = rl_start, &
      scal_divdamp_o2 = scal_divdamp_o2, &
      wgt_nnew_rth = wgt_nnew_rth, &
      wgt_nnew_vel = wgt_nnew_vel, &
      wgt_nnow_rth = wgt_nnow_rth, &
      wgt_nnow_vel = wgt_nnow_vel, &
      z_a = z_a, &
      z_b = z_b, &
      z_c = z_c, &
      z_d_vn_dmp = z_d_vn_dmp, &
      z_d_vn_iau = z_d_vn_iau, &
      z_ddt_vn_apc = z_ddt_vn_apc, &
      z_ddt_vn_cor = z_ddt_vn_cor, &
      z_ddt_vn_dyn = z_ddt_vn_dyn, &
      z_ddt_vn_pgr = z_ddt_vn_pgr, &
      z_ddt_vn_ray = z_ddt_vn_ray, &
      z_g = z_g, &
      z_gamma = z_gamma, &
      z_ntdistv_bary_1 = z_ntdistv_bary_1, &
      z_ntdistv_bary_2 = z_ntdistv_bary_2, &
      z_rho_tavg = z_rho_tavg, &
      z_rho_tavg_m1 = z_rho_tavg_m1, &
      z_theta1 = z_theta1, &
      z_theta2 = z_theta2, &
      z_theta_tavg = z_theta_tavg, &
      z_theta_tavg_m1 = z_theta_tavg_m1, &
      z_theta_v_pr_mc = z_theta_v_pr_mc, &
      z_theta_v_pr_mc_m1 = z_theta_v_pr_mc_m1, &
      z_w_backtraj = z_w_backtraj, &
      zf = zf &
    )
    end if

    call dace_program_solve_nh_predictor_pre( &
      state = dace_state, &
      bdy_divdamp = copy_or_ptr_bdy_divdamp, &
      enh_divdamp_fac = copy_or_ptr_enh_divdamp_fac, &
      global_data = copy_or_ptr_global_data, &
      p_int = copy_or_ptr_p_int, &
      p_nh = copy_or_ptr_p_nh, &
      p_nh_prog_nnew = copy_or_ptr_p_nh_prog_nnew, &
      p_nh_prog_nnow = copy_or_ptr_p_nh_prog_nnow, &
      p_patch = copy_or_ptr_p_patch, &
      prep_adv = copy_or_ptr_prep_adv, &
      scal_divdamp = copy_or_ptr_scal_divdamp, &
      z_alpha = copy_or_ptr_z_alpha, &
      z_beta = copy_or_ptr_z_beta, &
      z_contr_w_fl_l = copy_or_ptr_z_contr_w_fl_l, &
      z_dexner_dz_c = copy_or_ptr_z_dexner_dz_c, &
      z_dwdz_dd = copy_or_ptr_z_dwdz_dd, &
      z_exner_ex_pr = copy_or_ptr_z_exner_ex_pr, &
      z_exner_expl = copy_or_ptr_z_exner_expl, &
      z_exner_ic = copy_or_ptr_z_exner_ic, &
      z_flxdiv_mass = copy_or_ptr_z_flxdiv_mass, &
      z_flxdiv_theta = copy_or_ptr_z_flxdiv_theta, &
      z_grad_rth = copy_or_ptr_z_grad_rth, &
      z_graddiv2_vn = copy_or_ptr_z_graddiv2_vn, &
      z_graddiv_vn = copy_or_ptr_z_graddiv_vn, &
      z_gradh_exner = copy_or_ptr_z_gradh_exner, &
      z_hydro_corr = copy_or_ptr_z_hydro_corr, &
      z_kin_hor_e = copy_or_ptr_z_kin_hor_e, &
      z_mflx_top = copy_or_ptr_z_mflx_top, &
      z_q = copy_or_ptr_z_q, &
      z_raylfac = copy_or_ptr_z_raylfac, &
      z_rho_e = copy_or_ptr_z_rho_e, &
      z_rho_expl = copy_or_ptr_z_rho_expl, &
      z_rho_v = copy_or_ptr_z_rho_v, &
      z_rth_pr = copy_or_ptr_z_rth_pr, &
      z_th_ddz_exner_c = copy_or_ptr_z_th_ddz_exner_c, &
      z_theta_v_e = copy_or_ptr_z_theta_v_e, &
      z_theta_v_fl_e = copy_or_ptr_z_theta_v_fl_e, &
      z_theta_v_pr_ic = copy_or_ptr_z_theta_v_pr_ic, &
      z_theta_v_v = copy_or_ptr_z_theta_v_v, &
      z_vn_avg = copy_or_ptr_z_vn_avg, &
      z_vt_ie = copy_or_ptr_z_vt_ie, &
      z_w_concorr_mc = copy_or_ptr_z_w_concorr_mc, &
      z_w_concorr_me = copy_or_ptr_z_w_concorr_me, &
      z_w_expl = copy_or_ptr_z_w_expl, &
      f2dace_OPTIONAL_lacc = f2dace_OPTIONAL_lacc, &
      alin = alin, &
      aqdr = aqdr, &
      bqdr = bqdr, &
      df32 = df32, &
      df42 = df42, &
      distv_bary_1 = distv_bary_1, &
      distv_bary_2 = distv_bary_2, &
      dt_linintp_ubc = dt_linintp_ubc, &
      dt_linintp_ubc_nnew = dt_linintp_ubc_nnew, &
      dt_linintp_ubc_nnow = dt_linintp_ubc_nnow, &
      dt_shift = dt_shift, &
      dthalf = dthalf, &
      dtime = dtime, &
      dz32 = dz32, &
      dz42 = dz42, &
      dzlin = dzlin, &
      dzqdr = dzqdr, &
      i_endblk = i_endblk, &
      i_endidx = i_endidx, &
      i_startblk = i_startblk, &
      i_startidx = i_startidx, &
      idyn_timestep = idyn_timestep, &
      ishift = ishift, &
      istep = istep, &
      jb = jb, &
      jc = jc, &
      je = je, &
      jg = jg, &
      jk = jk, &
      jk_start = jk_start, &
      jks = jks, &
      jstep = jstep, &
      l_child_vertnest = l_child_vertnest, &
      l_init = l_init, &
      l_recompute = l_recompute, &
      l_vert_nested = l_vert_nested, &
      lacc = DACE_OPT_PROXY_lacc, &
      lclean_mflx = lclean_mflx, &
      lprep_adv = lprep_adv, &
      lsave_mflx = lsave_mflx, &
      lvn_only = lvn_only, &
      lvn_pos = lvn_pos, &
      nblks_gradp = nblks_gradp, &
      nlen_gradp = nlen_gradp, &
      nlev = nlev, &
      nlevp1 = nlevp1, &
      nnew = nnew, &
      nnow = nnow, &
      nproma_gradp = nproma_gradp, &
      npromz_gradp = npromz_gradp, &
      nshift = nshift, &
      nshift_total = nshift_total, &
      ntl1 = ntl1, &
      ntl2 = ntl2, &
      nvar = nvar, &
      r_dtimensubsteps = r_dtimensubsteps, &
      r_nsubsteps = r_nsubsteps, &
      rl_end = rl_end, &
      rl_start = rl_start, &
      scal_divdamp_o2 = scal_divdamp_o2, &
      wgt_nnew_rth = wgt_nnew_rth, &
      wgt_nnew_vel = wgt_nnew_vel, &
      wgt_nnow_rth = wgt_nnow_rth, &
      wgt_nnow_vel = wgt_nnow_vel, &
      z_a = z_a, &
      z_b = z_b, &
      z_c = z_c, &
      z_d_vn_dmp = z_d_vn_dmp, &
      z_d_vn_iau = z_d_vn_iau, &
      z_ddt_vn_apc = z_ddt_vn_apc, &
      z_ddt_vn_cor = z_ddt_vn_cor, &
      z_ddt_vn_dyn = z_ddt_vn_dyn, &
      z_ddt_vn_pgr = z_ddt_vn_pgr, &
      z_ddt_vn_ray = z_ddt_vn_ray, &
      z_g = z_g, &
      z_gamma = z_gamma, &
      z_ntdistv_bary_1 = z_ntdistv_bary_1, &
      z_ntdistv_bary_2 = z_ntdistv_bary_2, &
      z_rho_tavg = z_rho_tavg, &
      z_rho_tavg_m1 = z_rho_tavg_m1, &
      z_theta1 = z_theta1, &
      z_theta2 = z_theta2, &
      z_theta_tavg = z_theta_tavg, &
      z_theta_tavg_m1 = z_theta_tavg_m1, &
      z_theta_v_pr_mc = z_theta_v_pr_mc, &
      z_theta_v_pr_mc_m1 = z_theta_v_pr_mc_m1, &
      z_w_backtraj = z_w_backtraj, &
      zf = zf &
    )
  end subroutine run_solve_nh_predictor_pre_verification

  subroutine verify_solve_nh_predictor_pre( &
    bdy_divdamp, &
    enh_divdamp_fac, &
    p_int, &
    p_nh, &
    p_nh_prog_nnew, &
    p_nh_prog_nnow, &
    p_patch, &
    prep_adv, &
    scal_divdamp, &
    z_alpha, &
    z_beta, &
    z_contr_w_fl_l, &
    z_dexner_dz_c, &
    z_dwdz_dd, &
    z_exner_ex_pr, &
    z_exner_expl, &
    z_exner_ic, &
    z_flxdiv_mass, &
    z_flxdiv_theta, &
    z_grad_rth, &
    z_graddiv2_vn, &
    z_graddiv_vn, &
    z_gradh_exner, &
    z_hydro_corr, &
    z_kin_hor_e, &
    z_mflx_top, &
    z_q, &
    z_raylfac, &
    z_rho_e, &
    z_rho_expl, &
    z_rho_v, &
    z_rth_pr, &
    z_th_ddz_exner_c, &
    z_theta_v_e, &
    z_theta_v_fl_e, &
    z_theta_v_pr_ic, &
    z_theta_v_v, &
    z_vn_avg, &
    z_vt_ie, &
    z_w_concorr_mc, &
    z_w_concorr_me, &
    z_w_expl, &
    alin, &
    aqdr, &
    bqdr, &
    df32, &
    df42, &
    distv_bary_1, &
    distv_bary_2, &
    dt_linintp_ubc, &
    dt_linintp_ubc_nnew, &
    dt_linintp_ubc_nnow, &
    dt_shift, &
    dthalf, &
    dtime, &
    dz32, &
    dz42, &
    dzlin, &
    dzqdr, &
    i_endblk, &
    i_endidx, &
    i_startblk, &
    i_startidx, &
    idyn_timestep, &
    ishift, &
    istep, &
    jb, &
    jc, &
    je, &
    jg, &
    jk, &
    jk_start, &
    jks, &
    jstep, &
    l_child_vertnest, &
    l_init, &
    l_recompute, &
    l_vert_nested, &
    lacc, &
    lclean_mflx, &
    lprep_adv, &
    lsave_mflx, &
    lvn_only, &
    lvn_pos, &
    nblks_gradp, &
    nlen_gradp, &
    nlev, &
    nlevp1, &
    nnew, &
    nnow, &
    nproma_gradp, &
    npromz_gradp, &
    nshift, &
    nshift_total, &
    ntl1, &
    ntl2, &
    nvar, &
    r_dtimensubsteps, &
    r_nsubsteps, &
    rl_end, &
    rl_start, &
    scal_divdamp_o2, &
    wgt_nnew_rth, &
    wgt_nnew_vel, &
    wgt_nnow_rth, &
    wgt_nnow_vel, &
    z_a, &
    z_b, &
    z_c, &
    z_d_vn_dmp, &
    z_d_vn_iau, &
    z_ddt_vn_apc, &
    z_ddt_vn_cor, &
    z_ddt_vn_dyn, &
    z_ddt_vn_pgr, &
    z_ddt_vn_ray, &
    z_g, &
    z_gamma, &
    z_ntdistv_bary_1, &
    z_ntdistv_bary_2, &
    z_rho_tavg, &
    z_rho_tavg_m1, &
    z_theta1, &
    z_theta2, &
    z_theta_tavg, &
    z_theta_tavg_m1, &
    z_theta_v_pr_mc, &
    z_theta_v_pr_mc_m1, &
    z_w_backtraj, &
    zf &
  )
    real(kind=c_double), dimension(:), target :: bdy_divdamp
    real(kind=c_double), dimension(:), target :: enh_divdamp_fac
    type(t_int_state), target :: p_int
    type(t_nh_state), target :: p_nh
    type(t_nh_prog), target :: p_nh_prog_nnew
    type(t_nh_prog), target :: p_nh_prog_nnow
    type(t_patch), target :: p_patch
    type(t_prepare_adv), target :: prep_adv
    real(kind=c_double), dimension(:), target :: scal_divdamp
    real(kind=c_double), dimension(:,:), target :: z_alpha
    real(kind=c_double), dimension(:,:), target :: z_beta
    real(kind=c_double), dimension(:,:), target :: z_contr_w_fl_l
    real(kind=c_double), dimension(:,:,:,:), target :: z_dexner_dz_c
    real(kind=c_double), dimension(:,:,:), target :: z_dwdz_dd
    real(kind=c_double), dimension(:,:,:), target :: z_exner_ex_pr
    real(kind=c_double), dimension(:,:), target :: z_exner_expl
    real(kind=c_double), dimension(:,:), target :: z_exner_ic
    real(kind=c_double), dimension(:,:), target :: z_flxdiv_mass
    real(kind=c_double), dimension(:,:), target :: z_flxdiv_theta
    real(kind=c_double), dimension(:,:,:,:), target :: z_grad_rth
    real(kind=c_double), dimension(:,:), target :: z_graddiv2_vn
    real(kind=c_double), dimension(:,:,:), target :: z_graddiv_vn
    real(kind=c_double), dimension(:,:,:), target :: z_gradh_exner
    real(kind=c_double), dimension(:,:), target :: z_hydro_corr
    real(kind=c_double), dimension(:,:,:), target :: z_kin_hor_e
    real(kind=c_double), dimension(:,:), target :: z_mflx_top
    real(kind=c_double), dimension(:,:), target :: z_q
    real(kind=c_double), dimension(:), target :: z_raylfac
    real(kind=c_double), dimension(:,:,:), target :: z_rho_e
    real(kind=c_double), dimension(:,:), target :: z_rho_expl
    real(kind=c_double), dimension(:,:,:), target :: z_rho_v
    real(kind=c_double), dimension(:,:,:,:), target :: z_rth_pr
    real(kind=c_double), dimension(:,:,:), target :: z_th_ddz_exner_c
    real(kind=c_double), dimension(:,:,:), target :: z_theta_v_e
    real(kind=c_double), dimension(:,:,:), target :: z_theta_v_fl_e
    real(kind=c_double), dimension(:,:), target :: z_theta_v_pr_ic
    real(kind=c_double), dimension(:,:,:), target :: z_theta_v_v
    real(kind=c_double), dimension(:,:), target :: z_vn_avg
    real(kind=c_double), dimension(:,:,:), target :: z_vt_ie
    real(kind=c_double), dimension(:,:), target :: z_w_concorr_mc
    real(kind=c_double), dimension(:,:,:), target :: z_w_concorr_me
    real(kind=c_double), dimension(:,:), target :: z_w_expl
    real(kind=c_double) :: alin
    real(kind=c_double) :: aqdr
    real(kind=c_double) :: bqdr
    real(kind=c_double) :: df32
    real(kind=c_double) :: df42
    real(kind=c_double) :: distv_bary_1
    real(kind=c_double) :: distv_bary_2
    real(kind=c_double) :: dt_linintp_ubc
    real(kind=c_double) :: dt_linintp_ubc_nnew
    real(kind=c_double) :: dt_linintp_ubc_nnow
    real(kind=c_double) :: dt_shift
    real(kind=c_double) :: dthalf
    real(kind=c_double) :: dtime
    real(kind=c_double) :: dz32
    real(kind=c_double) :: dz42
    real(kind=c_double) :: dzlin
    real(kind=c_double) :: dzqdr
    integer(kind=c_int) :: i_endblk
    integer(kind=c_int) :: i_endidx
    integer(kind=c_int) :: i_startblk
    integer(kind=c_int) :: i_startidx
    integer(kind=c_int) :: idyn_timestep
    integer(kind=c_int) :: ishift
    integer(kind=c_int) :: istep
    integer(kind=c_int) :: jb
    integer(kind=c_int) :: jc
    integer(kind=c_int) :: je
    integer(kind=c_int) :: jg
    integer(kind=c_int) :: jk
    integer(kind=c_int) :: jk_start
    integer(kind=c_int) :: jks
    integer(kind=c_int) :: jstep
    integer(kind=c_int) :: l_child_vertnest
    integer(kind=c_int) :: l_init
    integer(kind=c_int) :: l_recompute
    integer(kind=c_int) :: l_vert_nested
    integer(kind=c_int), optional :: lacc
    integer(kind=c_int) :: lclean_mflx
    integer(kind=c_int) :: lprep_adv
    integer(kind=c_int) :: lsave_mflx
    integer(kind=c_int) :: lvn_only
    integer(kind=c_int) :: lvn_pos
    integer(kind=c_int) :: nblks_gradp
    integer(kind=c_int) :: nlen_gradp
    integer(kind=c_int) :: nlev
    integer(kind=c_int) :: nlevp1
    integer(kind=c_int) :: nnew
    integer(kind=c_int) :: nnow
    integer(kind=c_int) :: nproma_gradp
    integer(kind=c_int) :: npromz_gradp
    integer(kind=c_int) :: nshift
    integer(kind=c_int) :: nshift_total
    integer(kind=c_int) :: ntl1
    integer(kind=c_int) :: ntl2
    integer(kind=c_int) :: nvar
    real(kind=c_double) :: r_dtimensubsteps
    real(kind=c_double) :: r_nsubsteps
    integer(kind=c_int) :: rl_end
    integer(kind=c_int) :: rl_start
    real(kind=c_double) :: scal_divdamp_o2
    real(kind=c_double) :: wgt_nnew_rth
    real(kind=c_double) :: wgt_nnew_vel
    real(kind=c_double) :: wgt_nnow_rth
    real(kind=c_double) :: wgt_nnow_vel
    real(kind=c_double) :: z_a
    real(kind=c_double) :: z_b
    real(kind=c_double) :: z_c
    real(kind=c_double) :: z_d_vn_dmp
    real(kind=c_double) :: z_d_vn_iau
    real(kind=c_double) :: z_ddt_vn_apc
    real(kind=c_double) :: z_ddt_vn_cor
    real(kind=c_double) :: z_ddt_vn_dyn
    real(kind=c_double) :: z_ddt_vn_pgr
    real(kind=c_double) :: z_ddt_vn_ray
    real(kind=c_double) :: z_g
    real(kind=c_double) :: z_gamma
    real(kind=c_double) :: z_ntdistv_bary_1
    real(kind=c_double) :: z_ntdistv_bary_2
    real(kind=c_double) :: z_rho_tavg
    real(kind=c_double) :: z_rho_tavg_m1
    real(kind=c_double) :: z_theta1
    real(kind=c_double) :: z_theta2
    real(kind=c_double) :: z_theta_tavg
    real(kind=c_double) :: z_theta_tavg_m1
    real(kind=c_double) :: z_theta_v_pr_mc
    real(kind=c_double) :: z_theta_v_pr_mc_m1
    real(kind=c_double) :: z_w_backtraj
    real(kind=c_double) :: zf
    !Optional helper parameter
    !WARNING: HACKFIX, POSSIBLE EXPLOSION
    integer(kind=c_int) :: f2dace_OPTIONAL_lacc

    integer(kind=c_int) :: DACE_OPT_PROXY_lacc

    logical :: local_result, result

    !$ACC WAIT

    result = .true.
    local_result = .true.

    if (present(lacc)) then
      f2dace_OPTIONAL_lacc = 1
      DACE_OPT_PROXY_lacc = lacc
      
    else
      f2dace_OPTIONAL_lacc = 0
    end if


    if (is_initialized .eqv. .false.) then
      print *, "verify_solve_nh_predictor_pre: dace state is not initialized"
    end if

    call check_initializations()

#ifndef _OPENACC

    call compare_float64_1d_array( &
        actual=copy_or_ptr_bdy_divdamp, &
        ref=bdy_divdamp, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="bdy_divdamp" &
    )

#else

    call compare_float64_1d_array( &
        actual=copy_or_ptr_bdy_divdamp, &
        ref=bdy_divdamp, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="bdy_divdamp" &
    )

#endif

    result = result .and. local_result

    call compare_float64_1d_array( &
        actual=copy_or_ptr_enh_divdamp_fac, &
        ref=enh_divdamp_fac, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="enh_divdamp_fac" &
    )

    result = result .and. local_result

    call compare_global_data_type_struct( &
        actual=copy_or_ptr_global_data, &
        result=local_result, &
        struct_expr="global_data" &
    )

    result = result .and. local_result

    call compare_t_int_state_struct( &
        actual=copy_or_ptr_p_int, &
        ref=p_int, &
        result=local_result, &
        struct_expr="p_int" &
    )

    result = result .and. local_result

    call compare_t_nh_state_struct( &
        actual=copy_or_ptr_p_nh, &
        ref=p_nh, &
        result=local_result, &
        struct_expr="p_nh" &
    )

    result = result .and. local_result

    call compare_t_nh_prog_struct( &
        actual=copy_or_ptr_p_nh_prog_nnew, &
        ref=p_nh_prog_nnew, &
        result=local_result, &
        struct_expr="p_nh_prog_nnew" &
    )

    result = result .and. local_result

    call compare_t_nh_prog_struct( &
        actual=copy_or_ptr_p_nh_prog_nnow, &
        ref=p_nh_prog_nnow, &
        result=local_result, &
        struct_expr="p_nh_prog_nnow" &
    )

    result = result .and. local_result

    call compare_t_patch_struct( &
        actual=copy_or_ptr_p_patch, &
        ref=p_patch, &
        result=local_result, &
        struct_expr="p_patch" &
    )

    result = result .and. local_result

    call compare_t_prepare_adv_struct( &
        actual=copy_or_ptr_prep_adv, &
        ref=prep_adv, &
        result=local_result, &
        struct_expr="prep_adv" &
    )

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_1d_array( &
        actual=copy_or_ptr_scal_divdamp, &
        ref=scal_divdamp, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="scal_divdamp" &
    )

#else

    call compare_float64_1d_array( &
        actual=copy_or_ptr_scal_divdamp, &
        ref=scal_divdamp, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="scal_divdamp" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_alpha, &
        ref=z_alpha, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_alpha" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_alpha, &
        ref=z_alpha, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_alpha" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_beta, &
        ref=z_beta, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_beta" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_beta, &
        ref=z_beta, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_beta" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_contr_w_fl_l, &
        ref=z_contr_w_fl_l, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_contr_w_fl_l" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_contr_w_fl_l, &
        ref=z_contr_w_fl_l, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_contr_w_fl_l" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_4d_array( &
        actual=copy_or_ptr_z_dexner_dz_c, &
        ref=z_dexner_dz_c, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_dexner_dz_c" &
    )

#else

    call compare_float64_4d_array( &
        actual=copy_or_ptr_z_dexner_dz_c, &
        ref=z_dexner_dz_c, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_dexner_dz_c" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_dwdz_dd, &
        ref=z_dwdz_dd, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_dwdz_dd" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_dwdz_dd, &
        ref=z_dwdz_dd, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_dwdz_dd" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_exner_ex_pr, &
        ref=z_exner_ex_pr, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_exner_ex_pr" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_exner_ex_pr, &
        ref=z_exner_ex_pr, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_exner_ex_pr" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_exner_expl, &
        ref=z_exner_expl, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_exner_expl" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_exner_expl, &
        ref=z_exner_expl, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_exner_expl" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_exner_ic, &
        ref=z_exner_ic, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_exner_ic" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_exner_ic, &
        ref=z_exner_ic, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_exner_ic" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_flxdiv_mass, &
        ref=z_flxdiv_mass, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_flxdiv_mass" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_flxdiv_mass, &
        ref=z_flxdiv_mass, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_flxdiv_mass" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_flxdiv_theta, &
        ref=z_flxdiv_theta, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_flxdiv_theta" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_flxdiv_theta, &
        ref=z_flxdiv_theta, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_flxdiv_theta" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_4d_array( &
        actual=copy_or_ptr_z_grad_rth, &
        ref=z_grad_rth, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_grad_rth" &
    )

#else

    call compare_float64_4d_array( &
        actual=copy_or_ptr_z_grad_rth, &
        ref=z_grad_rth, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_grad_rth" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_graddiv2_vn, &
        ref=z_graddiv2_vn, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_graddiv2_vn" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_graddiv2_vn, &
        ref=z_graddiv2_vn, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_graddiv2_vn" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_graddiv_vn, &
        ref=z_graddiv_vn, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_graddiv_vn" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_graddiv_vn, &
        ref=z_graddiv_vn, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_graddiv_vn" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_gradh_exner, &
        ref=z_gradh_exner, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_gradh_exner" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_gradh_exner, &
        ref=z_gradh_exner, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_gradh_exner" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_hydro_corr, &
        ref=z_hydro_corr, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_hydro_corr" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_hydro_corr, &
        ref=z_hydro_corr, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_hydro_corr" &
    )

#endif

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

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_mflx_top, &
        ref=z_mflx_top, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_mflx_top" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_mflx_top, &
        ref=z_mflx_top, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_mflx_top" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_q, &
        ref=z_q, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_q" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_q, &
        ref=z_q, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_q" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_1d_array( &
        actual=copy_or_ptr_z_raylfac, &
        ref=z_raylfac, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_raylfac" &
    )

#else

    call compare_float64_1d_array( &
        actual=copy_or_ptr_z_raylfac, &
        ref=z_raylfac, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_raylfac" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_rho_e, &
        ref=z_rho_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_rho_e" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_rho_e, &
        ref=z_rho_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_rho_e" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_rho_expl, &
        ref=z_rho_expl, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_rho_expl" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_rho_expl, &
        ref=z_rho_expl, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_rho_expl" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_rho_v, &
        ref=z_rho_v, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_rho_v" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_rho_v, &
        ref=z_rho_v, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_rho_v" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_4d_array( &
        actual=copy_or_ptr_z_rth_pr, &
        ref=z_rth_pr, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_rth_pr" &
    )

#else

    call compare_float64_4d_array( &
        actual=copy_or_ptr_z_rth_pr, &
        ref=z_rth_pr, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_rth_pr" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_th_ddz_exner_c, &
        ref=z_th_ddz_exner_c, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_th_ddz_exner_c" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_th_ddz_exner_c, &
        ref=z_th_ddz_exner_c, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_th_ddz_exner_c" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_theta_v_e, &
        ref=z_theta_v_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_theta_v_e" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_theta_v_e, &
        ref=z_theta_v_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_theta_v_e" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_theta_v_fl_e, &
        ref=z_theta_v_fl_e, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_theta_v_fl_e" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_theta_v_fl_e, &
        ref=z_theta_v_fl_e, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_theta_v_fl_e" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_theta_v_pr_ic, &
        ref=z_theta_v_pr_ic, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_theta_v_pr_ic" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_theta_v_pr_ic, &
        ref=z_theta_v_pr_ic, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_theta_v_pr_ic" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_theta_v_v, &
        ref=z_theta_v_v, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_theta_v_v" &
    )

#else

    call compare_float64_3d_array( &
        actual=copy_or_ptr_z_theta_v_v, &
        ref=z_theta_v_v, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_theta_v_v" &
    )

#endif

    result = result .and. local_result
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_vn_avg, &
        ref=z_vn_avg, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_vn_avg" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_vn_avg, &
        ref=z_vn_avg, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_vn_avg" &
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

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_w_concorr_mc, &
        ref=z_w_concorr_mc, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_w_concorr_mc" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_w_concorr_mc, &
        ref=z_w_concorr_mc, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_w_concorr_mc" &
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
#ifndef _OPENACC

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_w_expl, &
        ref=z_w_expl, &
        result=local_result, &
        use_openacc=.false., &
        array_expr="z_w_expl" &
    )

#else

    call compare_float64_2d_array( &
        actual=copy_or_ptr_z_w_expl, &
        ref=z_w_expl, &
        result=local_result, &
        use_openacc=.true., &
        array_expr="z_w_expl" &
    )

#endif

    result = result .and. local_result

    if (.not. result) then
      print *, "verify_solve_nh_predictor_pre: Failed verification"
    else
      print *, "verify_solve_nh_predictor_pre: Verification successful :)"
    end if

  end subroutine verify_solve_nh_predictor_pre

end module mo_solve_nh_predictor_pre_bindings
