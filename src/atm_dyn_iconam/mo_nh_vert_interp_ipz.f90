!
! This module contains routines for the vertical interpolation of
! atmospheric data provided by external analyses to the ICON grid
!
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

!----------------------------
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

!----------------------------

MODULE mo_nh_vert_interp_ipz

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_opt_diagnostics,     ONLY: t_vcoeff, vcoeff_allocate
  USE mo_intp_data_strc,      ONLY: t_int_state, p_int_state
  USE mo_intp,                ONLY: cell_avg, cells2edges_scalar
  USE mo_parallel_config,     ONLY: nproma
  USE mo_physical_constants,  ONLY: grav, rd, dtdz_standardatm
  USE mo_run_config,          ONLY: iforcing, num_lev
  USE mo_io_config,           ONLY: itype_pres_msl
  USE mo_grid_config,         ONLY: l_limited_area, grid_sphere_radius
  USE mo_impl_constants,      ONLY: inwp, iaes, PRES_MSL_METHOD_GME, PRES_MSL_METHOD_IFS,   &
    &                               PRES_MSL_METHOD_DWD, PRES_MSL_METHOD_IFS_CORR,          &
    &                               SUCCESS
  USE mo_exception,           ONLY: finish
  USE mo_fortran_tools,       ONLY: copy, init, assert_acc_device_only, assert_lacc_equals_i_am_accel_node, &
    &                               set_acc_host_or_device, assert_acc_host_only, minval_1d
  USE mo_initicon_config,     ONLY: zpbl1, zpbl2
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_sync,                ONLY: SYNC_E, sync_patch_array_mult
  USE mo_nh_vert_interp,      ONLY: prepare_lin_intp, prepare_extrap, prepare_cubic_intp, &
    &                               temperature_intp, prepare_extrap_ifspp, pressure_intp
  USE mo_dynamics_config,     ONLY: ldeepatmo
  USE mo_deepatmo,            ONLY: deepatmo_htrafo

  IMPLICIT NONE
  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nh_vert_interp_ipz'


  ! Threshold for switching between analytical formulas for constant temperature and
  ! constant vertical gradient of temperature, respectively
  REAL(wp), PARAMETER :: dtvdz_thresh = 1.e-4_wp ! 0.1 K/km

  ! Artificial limits on temperature profile used for extrapolation below the ground
  REAL(wp), PARAMETER :: t_low  = 255.0_wp
  REAL(wp), PARAMETER :: t_high = 290.5_wp

  PUBLIC :: prepare_vert_interp_z
  PUBLIC :: prepare_vert_interp_p
  PUBLIC :: prepare_vert_interp_i
  PUBLIC :: z_at_plevels


CONTAINS

  !-------------
  !>
  !! Compute coefficients for vertical interpolation of model-level
  !! fields to constant-height levels.
  !!
  !! @note Computation of coefficients already performs interpolation
  !! of temperature, QV (prognostic+diagnostic) and pressure onto z-levels.
  !!

  SUBROUTINE prepare_vert_interp_z(p_patch, p_diag, p_metrics, intp_hrz, nzlev,  &
    &                              temp_z_out, pres_z_out, p_z3d_out, vcoeff_z, lacc)

    TYPE(t_patch),        TARGET,      INTENT(IN)    :: p_patch
    TYPE(t_nh_diag),      POINTER                    :: p_diag
    TYPE(t_nh_metrics),                INTENT(IN)    :: p_metrics
    TYPE(t_int_state),                 INTENT(IN)    :: intp_hrz          !< pointer to data structure for interpolation
    INTEGER,                           INTENT(IN)    :: nzlev             !< number of output levels (height)
    REAL(wp),                          INTENT(INOUT) :: temp_z_out(:,:,:),          &
      &                                                 pres_z_out(:,:,:)
    REAL(wp) ,            POINTER                    :: p_z3d_out(:,:,:)  !< height field
    TYPE(t_vcoeff),                    INTENT(INOUT) :: vcoeff_z
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    ! LOCAL VARIABLES
    CHARACTER(*), PARAMETER :: routine = modname//":prepare_vert_interp_z"

    INTEGER :: nlev, nlevp1, jg, i_endblk, nblks_c,nblks_e, npromz_c,npromz_e !< blocking parameters
    ! Auxiliary field for output data
    REAL(wp), DIMENSION(nproma,nzlev,p_patch%nblks_c)        :: z_auxz
    REAL(wp), DIMENSION(nproma,nzlev,p_patch%nblks_e)        :: p_z3d_edge
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_me
    ! Pointer to virtual temperature / temperature, depending on whether the run is moist or dry
    REAL(wp), POINTER, DIMENSION(:,:,:) :: ptr_tempv

    !-------------------------------------------------------------------------
    CALL assert_acc_device_only(routine, lacc)
    CALL assert_lacc_equals_i_am_accel_node(routine, lacc)

    vcoeff_z%l_initialized = .TRUE.
    IF (p_patch%n_patch_cells==0) RETURN

    nlev     = p_patch%nlev
    nlevp1   = p_patch%nlevp1
    nblks_c  = p_patch%nblks_c
    npromz_c = p_patch%npromz_c
    nblks_e  = p_patch%nblks_e
    npromz_e = p_patch%npromz_e
    jg       = p_patch%id

    IF (  iforcing == inwp .OR. iforcing == iaes  ) THEN
      ptr_tempv => p_diag%tempv(:,:,:)
    ELSE
      ptr_tempv => p_diag%temp(:,:,:)
    END IF
    !$ACC DATA PRESENT(p_patch, p_diag, p_metrics, temp_z_out, pres_z_out, p_z3d_out) &
    !$ACC   CREATE(z_auxz, p_z3d_edge, z_me)

    !--- Coefficients: Interpolation to z-level fields

    ! allocate coefficient table:
    CALL vcoeff_allocate(nblks_c, nblks_e, nzlev, vcoeff_z, lacc=.TRUE.)

    ! Prepare interpolation coefficients
    CALL prepare_lin_intp(p_metrics%z_mc, p_z3d_out,                                         & !in
      &                   nblks_c, npromz_c, nlev, nzlev,                                    & !in
      &                   vcoeff_z%lin_cell%wfac_lin,                                        & !out
      &                   vcoeff_z%lin_cell%idx0_lin,                                        & !out
      &                   vcoeff_z%lin_cell%bot_idx_lin,                                     & !out
      &                   lacc=.TRUE.)
    IF ( ANY((/PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR, PRES_MSL_METHOD_DWD/)== itype_pres_msl) ) THEN
      CALL prepare_extrap_ifspp(p_metrics%z_ifc, p_metrics%z_mc, nblks_c, npromz_c, nlev,    & !in
        &                       vcoeff_z%lin_cell%kpbl1, vcoeff_z%lin_cell%zextrap,          & !out
        &                       vcoeff_z%lin_cell%wfacpbl1,                                  & !out
        &                       lacc=.TRUE.)
      !$OMP PARALLEL
      CALL init(vcoeff_z%lin_cell%kpbl2(:,:), opt_acc_async=.TRUE.)
      CALL init(vcoeff_z%lin_cell%wfacpbl2(:,:), opt_acc_async=.TRUE.)
      !$OMP END PARALLEL
    ELSE
      CALL prepare_extrap(p_metrics%z_mc, nblks_c, npromz_c, nlev,                           & !in
        &                 vcoeff_z%lin_cell%kpbl1, vcoeff_z%lin_cell%wfacpbl1,               & !out
        &                 vcoeff_z%lin_cell%kpbl2, vcoeff_z%lin_cell%wfacpbl2,               & !out
        &                  lacc=.TRUE.)
    ENDIF
    CALL prepare_cubic_intp(p_metrics%z_mc, p_z3d_out, nblks_c, npromz_c, nlev, nzlev,       & !in
      &                     vcoeff_z%cub_cell%coef1, vcoeff_z%cub_cell%coef2,                & !out
      &                     vcoeff_z%cub_cell%coef3,                                         & !out
      &                     vcoeff_z%cub_cell%idx0_cub, vcoeff_z%cub_cell%bot_idx_cub,       & !out
      &                     lacc=.TRUE.)

    ! Perform vertical interpolation

    CALL temperature_intp(p_diag%temp, z_auxz, p_metrics%z_mc, p_z3d_out,                   & !in,out
      &                   nblks_c, npromz_c, nlev, nzlev,                                   & !in
      &                   vcoeff_z%cub_cell%coef1, vcoeff_z%cub_cell%coef2,                 & !in
      &                   vcoeff_z%cub_cell%coef3, vcoeff_z%lin_cell%wfac_lin,              & !in
      &                   vcoeff_z%cub_cell%idx0_cub, vcoeff_z%lin_cell%idx0_lin,           & !in
      &                   vcoeff_z%cub_cell%bot_idx_cub, vcoeff_z%lin_cell%bot_idx_lin,     & !in
      &                   vcoeff_z%lin_cell%wfacpbl1, vcoeff_z%lin_cell%kpbl1,              & !in
      &                   vcoeff_z%lin_cell%wfacpbl2, vcoeff_z%lin_cell%kpbl2,              & !in
      &                   l_restore_sfcinv=.FALSE., l_hires_corr=.FALSE.,                   & !in
      &                   extrapol_dist=0._wp, l_pz_mode=.TRUE.,                            & !in
      &                   zextrap=vcoeff_z%lin_cell%zextrap, lacc=.TRUE.)


    IF (jg > 1 .OR. l_limited_area) THEN ! copy outermost nest boundary row in order to avoid missing values
      i_endblk = p_patch%cells%end_blk(1,1)
      !$OMP PARALLEL
      CALL copy(z_auxz(:,:,1:i_endblk), temp_z_out(:,:,1:i_endblk), opt_acc_async=.TRUE.)
      !$OMP END PARALLEL
    ENDIF
   
    CALL cell_avg(z_auxz, p_patch, p_int_state(jg)%c_bln_avg, temp_z_out, lacc=.TRUE.)

    ! Interpolate pressure on z-levels
    CALL pressure_intp(p_diag%pres, ptr_tempv, p_metrics%z_mc,                              & !in
      &                z_auxz, p_z3d_out,                                                   & !out,in
      &                nblks_c, npromz_c, nlev, nzlev,                                      & !in
      &                vcoeff_z%lin_cell%wfac_lin,    vcoeff_z%lin_cell%idx0_lin,           & !in
      &                vcoeff_z%lin_cell%bot_idx_lin, vcoeff_z%lin_cell%wfacpbl1,           & !in
      &                vcoeff_z%lin_cell%kpbl1, vcoeff_z%lin_cell%wfacpbl2,                 & !in
      &                vcoeff_z%lin_cell%kpbl2,                                             & !in
      &                zextrap=vcoeff_z%lin_cell%zextrap,                                   & !optin
      &                lacc=.TRUE.                                                          ) !optin

    IF (jg > 1 .OR. l_limited_area) THEN ! copy outermost nest boundary row in order to avoid missing values
      i_endblk = p_patch%cells%end_blk(1,1)
      !$OMP PARALLEL
      CALL copy(z_auxz(:,:,1:i_endblk), pres_z_out(:,:,1:i_endblk), opt_acc_async=.TRUE.)
      !$OMP END PARALLEL
    ENDIF

    CALL cell_avg(z_auxz, p_patch, p_int_state(jg)%c_bln_avg, pres_z_out, lacc=.TRUE.)

    IF ( ANY((/PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR, PRES_MSL_METHOD_DWD/)== itype_pres_msl) ) THEN
      CALL prepare_extrap(p_metrics%z_mc, nblks_c, npromz_c, nlev,                           & !in
        &                 vcoeff_z%lin_cell%kpbl1, vcoeff_z%lin_cell%wfacpbl1,               & !out
        &                 vcoeff_z%lin_cell%kpbl2, vcoeff_z%lin_cell%wfacpbl2,               & !out
        &                 lacc=.TRUE.)
    ENDIF

    !--- Interpolation data for the vertical interface of cells, "nlevp1"

    CALL prepare_lin_intp(p_metrics%z_ifc, p_z3d_out, nblks_c, npromz_c, nlevp1, nzlev,     & !in
      &                   vcoeff_z%lin_cell_nlevp1%wfac_lin,                                & !out
      &                   vcoeff_z%lin_cell_nlevp1%idx0_lin,                                & !out
      &                   vcoeff_z%lin_cell_nlevp1%bot_idx_lin,                             & !out
      &                   lacc=.TRUE. )
    CALL prepare_extrap(p_metrics%z_ifc, nblks_c, npromz_c, nlevp1,                         & !in
      &                 vcoeff_z%lin_cell_nlevp1%kpbl1, vcoeff_z%lin_cell_nlevp1%wfacpbl1,  & !out
      &                 vcoeff_z%lin_cell_nlevp1%kpbl2, vcoeff_z%lin_cell_nlevp1%wfacpbl2,  & !out
      &                 lacc=.TRUE.)

    !--- Compute data for interpolation of edge-based fields

    ! Compute geometric height at edge points
    CALL cells2edges_scalar(p_metrics%z_mc, p_patch, intp_hrz%c_lin_e, z_me, opt_fill_latbc=.TRUE., lacc=.TRUE.)
    CALL cells2edges_scalar(p_z3d_out, p_patch, intp_hrz%c_lin_e, p_z3d_edge, opt_fill_latbc=.TRUE., lacc=.TRUE.)
    CALL sync_patch_array_mult(SYNC_E,p_patch,2,z_me,p_z3d_edge)

    CALL prepare_lin_intp(z_me, p_z3d_edge, nblks_e, npromz_e, nlev, nzlev,                 & !in
      &                   vcoeff_z%lin_edge%wfac_lin, vcoeff_z%lin_edge%idx0_lin,           & !out
      &                   vcoeff_z%lin_edge%bot_idx_lin,                                    & !out
      &                   lacc=.TRUE. )
    CALL prepare_extrap(z_me, nblks_e, npromz_e, nlev,                                      & !in
      &                 vcoeff_z%lin_edge%kpbl1, vcoeff_z%lin_edge%wfacpbl1,                & !out
      &                 vcoeff_z%lin_edge%kpbl2, vcoeff_z%lin_edge%wfacpbl2,                & !out
      &                 lacc=.TRUE.)
    CALL prepare_cubic_intp(z_me, p_z3d_edge, nblks_e, npromz_e, nlev, nzlev,               & !in
      &                     vcoeff_z%cub_edge%coef1, vcoeff_z%cub_edge%coef2,               & !out
      &                     vcoeff_z%cub_edge%coef3,                                        & !out
      &                     vcoeff_z%cub_edge%idx0_cub, vcoeff_z%cub_edge%bot_idx_cub,      & !out
      &                     lacc=.TRUE.)

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE prepare_vert_interp_z


  !-------------
  !>
  !! Compute coefficients for vertical interpolation of model-level
  !! fields to constant-pressure levels.
  !!
  SUBROUTINE prepare_vert_interp_p(p_patch, p_diag, p_metrics, intp_hrz, nplev,     &
    &                              gh_p_out, temp_p_out, p_p3d_out, vcoeff_p, lacc)

    TYPE(t_patch),        TARGET,      INTENT(IN)    :: p_patch
    TYPE(t_nh_diag),      POINTER                    :: p_diag
    TYPE(t_nh_metrics),                INTENT(IN)    :: p_metrics
    TYPE(t_int_state),                 INTENT(IN)    :: intp_hrz          !< pointer to data structure for interpolation
    INTEGER,                           INTENT(IN)    :: nplev             !< number of output levels (pres)
    REAL(wp),                          INTENT(INOUT) :: gh_p_out(:,:,:), &
      &                                                 temp_p_out(:,:,:)
    REAL(wp),             POINTER                    :: p_p3d_out(:,:,:)  !< pressure field
    TYPE(t_vcoeff),                    INTENT(INOUT) :: vcoeff_p          !< out
    LOGICAL,              OPTIONAL,    INTENT(IN)    :: lacc              !< If true, use openacc

    ! LOCAL VARIABLES
    CHARACTER(*), PARAMETER :: routine = modname//":prepare_vert_interp_p"

    INTEGER :: nlev, nlevp1, jg, i_endblk, nblks_c,nblks_e, npromz_c,npromz_e !< blocking parameters
    ! Auxiliary field for smoothing temperature and geopotential on pressure levels
    REAL(wp), DIMENSION(nproma,nplev,p_patch%nblks_c)        :: z_auxp
    REAL(wp), DIMENSION(nproma,nplev,p_patch%nblks_e)        :: gh_p_edge
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_me
    ! Pointer to virtual temperature / temperature, depending on whether the run is moist or dry
    REAL(wp), POINTER, DIMENSION(:,:,:) :: ptr_tempv

    CALL assert_acc_device_only(routine, lacc)
    CALL assert_lacc_equals_i_am_accel_node(routine, lacc)

    vcoeff_p%l_initialized = .TRUE.
    IF (p_patch%n_patch_cells==0) RETURN

    nlev     = p_patch%nlev
    nlevp1   = p_patch%nlevp1
    nblks_c  = p_patch%nblks_c
    npromz_c = p_patch%npromz_c
    nblks_e  = p_patch%nblks_e
    npromz_e = p_patch%npromz_e
    jg       = p_patch%id

    IF (  iforcing == inwp .OR. iforcing == iaes  ) THEN
      ptr_tempv => p_diag%tempv
    ELSE
      ptr_tempv => p_diag%temp
    ENDIF

    !$ACC DATA CREATE(z_auxp, gh_p_edge, z_me)

    ! allocate coefficient table:
    CALL vcoeff_allocate(nblks_c, nblks_e, nplev, vcoeff_p, lacc=.TRUE.)

    !--- Coefficients: Interpolation to pressure-level fields
    IF ( ANY((/PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR, PRES_MSL_METHOD_DWD/)== itype_pres_msl) ) THEN
      CALL prepare_extrap_ifspp(p_metrics%z_ifc, p_metrics%z_mc, nblks_c, npromz_c, nlev, & !in
        &                       vcoeff_p%lin_cell%kpbl1, vcoeff_p%lin_cell%zextrap,       & !out
        &                       vcoeff_p%lin_cell%wfacpbl1,                               & !out
        &                       lacc=.TRUE.)
      vcoeff_p%lin_cell%kpbl2(:,:) = 0
      vcoeff_p%lin_cell%wfacpbl2(:,:) = 0._wp
    ELSE
      CALL prepare_extrap(p_metrics%z_mc, nblks_c, npromz_c, nlev,                          & !in
        &                 vcoeff_p%lin_cell%kpbl1, vcoeff_p%lin_cell%wfacpbl1,              & !out
        &                 vcoeff_p%lin_cell%kpbl2, vcoeff_p%lin_cell%wfacpbl2,              & !out
        &                 lacc=.TRUE.)
    ENDIF

    ! Compute height at pressure levels (i.e. geopot/g); this height
    ! field is afterwards also used as target coordinate for vertical
    ! interpolation
    CALL z_at_plevels(p_diag%pres, ptr_tempv, p_metrics%z_mc,                               & !in
      &               p_p3d_out, z_auxp, nblks_c, npromz_c, nlev, nplev,                    & !in,out,in
      &               vcoeff_p%lin_cell%kpbl1, vcoeff_p%lin_cell%wfacpbl1,                  & !in
      &               vcoeff_p%lin_cell%kpbl2, vcoeff_p%lin_cell%wfacpbl2,                  & !in
      &               zextrap=vcoeff_p%lin_cell%zextrap,                                    & !optin
      &               lacc=.TRUE. )

    IF (jg > 1 .OR. l_limited_area) THEN ! copy outermost nest boundary row in order to avoid missing values
      i_endblk = p_patch%cells%end_blk(1,1)
      !$OMP PARALLEL
      CALL copy(z_auxp(:,:,1:i_endblk), gh_p_out(:,:,1:i_endblk), opt_acc_async=.TRUE.)
      !$OMP END PARALLEL
    ENDIF

    CALL cell_avg(z_auxp, p_patch, p_int_state(jg)%c_bln_avg, gh_p_out, lacc=.TRUE.)

    ! Prepare again interpolation coefficients (now for pressure levels)
    ! Note: the coefficients are partly identical to the z-level
    !       setup, but they are computed independently to avoid some
    !       complicated code.
    CALL prepare_lin_intp(p_metrics%z_mc, gh_p_out, nblks_c, npromz_c, nlev, nplev,         & !in
      &                   vcoeff_p%lin_cell%wfac_lin, vcoeff_p%lin_cell%idx0_lin,           & !out
      &                   vcoeff_p%lin_cell%bot_idx_lin,                                    & !out
      &                   lacc=.TRUE. )
    CALL prepare_cubic_intp(p_metrics%z_mc, gh_p_out, nblks_c, npromz_c, nlev, nplev,       & !in
      &                     vcoeff_p%cub_cell%coef1, vcoeff_p%cub_cell%coef2,               & !out
      &                     vcoeff_p%cub_cell%coef3,                                        & !out
      &                     vcoeff_p%cub_cell%idx0_cub, vcoeff_p%cub_cell%bot_idx_cub,      & !out
      &                     lacc=.TRUE.)

    ! Perform vertical interpolation
    CALL temperature_intp(p_diag%temp, z_auxp, p_metrics%z_mc, gh_p_out,                    & !in,out
      &                   nblks_c, npromz_c, nlev, nplev,                                   & !in
      &                   vcoeff_p%cub_cell%coef1, vcoeff_p%cub_cell%coef2,                 & !in
      &                   vcoeff_p%cub_cell%coef3,                                          & !in
      &                   vcoeff_p%lin_cell%wfac_lin, vcoeff_p%cub_cell%idx0_cub,           & !in
      &                   vcoeff_p%lin_cell%idx0_lin, vcoeff_p%cub_cell%bot_idx_cub,        & !in
      &                   vcoeff_p%lin_cell%bot_idx_lin, vcoeff_p%lin_cell%wfacpbl1,        & !in
      &                   vcoeff_p%lin_cell%kpbl1, vcoeff_p%lin_cell%wfacpbl2,              & !in
      &                   vcoeff_p%lin_cell%kpbl2, l_restore_sfcinv=.FALSE.,                & !in
      &                   l_hires_corr=.FALSE., extrapol_dist=0._wp, l_pz_mode=.TRUE.,      & !in
      &                   zextrap=vcoeff_p%lin_cell%zextrap, lacc=.TRUE.)

    IF (jg > 1 .OR. l_limited_area) THEN ! copy outermost nest boundary row in order to avoid missing values
      i_endblk = p_patch%cells%end_blk(1,1)
      !$OMP PARALLEL
      CALL copy(z_auxp(:,:,1:i_endblk), temp_p_out(:,:,1:i_endblk), opt_acc_async=.TRUE.)
      !$OMP END PARALLEL
    ENDIF

    CALL cell_avg(z_auxp, p_patch, p_int_state(jg)%c_bln_avg, temp_p_out, lacc=.TRUE.)

    IF ( ANY((/PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR, PRES_MSL_METHOD_DWD/)== itype_pres_msl) ) THEN
      CALL prepare_extrap(p_metrics%z_mc, nblks_c, npromz_c, nlev,                          & !in
        &                 vcoeff_p%lin_cell%kpbl1, vcoeff_p%lin_cell%wfacpbl1,              & !out
        &                 vcoeff_p%lin_cell%kpbl2, vcoeff_p%lin_cell%wfacpbl2,              & !out
        &                 lacc=.TRUE.)
    ENDIF

    !--- Interpolation data for the vertical interface of cells, "nlevp1"

    CALL prepare_lin_intp(p_metrics%z_ifc, gh_p_out, nblks_c, npromz_c, nlevp1, nplev,      & !in
      &                   vcoeff_p%lin_cell_nlevp1%wfac_lin,                                & !out
      &                   vcoeff_p%lin_cell_nlevp1%idx0_lin,                                & !out
      &                   vcoeff_p%lin_cell_nlevp1%bot_idx_lin,                             & !out
      &                   lacc=.TRUE. )
    CALL prepare_extrap(p_metrics%z_ifc, nblks_c, npromz_c, nlevp1,                         & !in
      &                 vcoeff_p%lin_cell_nlevp1%kpbl1, vcoeff_p%lin_cell_nlevp1%wfacpbl1,  & !out
      &                 vcoeff_p%lin_cell_nlevp1%kpbl2, vcoeff_p%lin_cell_nlevp1%wfacpbl2,  & !out
      &                 lacc=.TRUE.)

    !--- Compute data for interpolation of edge-based fields

    CALL cells2edges_scalar(p_metrics%z_mc, p_patch, intp_hrz%c_lin_e, z_me, opt_fill_latbc=.TRUE., lacc=.TRUE.)
    CALL cells2edges_scalar(gh_p_out, p_patch, intp_hrz%c_lin_e, gh_p_edge, opt_fill_latbc=.TRUE., lacc=.TRUE.)
    CALL sync_patch_array_mult(SYNC_E,p_patch,2,z_me,gh_p_edge)

    CALL prepare_lin_intp(z_me, gh_p_edge, nblks_e, npromz_e, nlev, nplev,                   & !in
      &                   vcoeff_p%lin_edge%wfac_lin, vcoeff_p%lin_edge%idx0_lin,            & !out
      &                   vcoeff_p%lin_edge%bot_idx_lin,                                     & !out
      &                   lacc=.TRUE. )
    CALL prepare_extrap(z_me, nblks_e, npromz_e, nlev,                                       & !in
      &                 vcoeff_p%lin_edge%kpbl1, vcoeff_p%lin_edge%wfacpbl1,                 & !out
      &                 vcoeff_p%lin_edge%kpbl2, vcoeff_p%lin_edge%wfacpbl2,                 & !out
      &                 lacc=.TRUE.)
    CALL prepare_cubic_intp(z_me, gh_p_edge, nblks_e, npromz_e, nlev, nplev,                 & !in
      &                     vcoeff_p%cub_edge%coef1, vcoeff_p%cub_edge%coef2,                & !out
      &                     vcoeff_p%cub_edge%coef3,                                         & !out
      &                     vcoeff_p%cub_edge%idx0_cub, vcoeff_p%cub_edge%bot_idx_cub,       & !out
      &                     lacc=.TRUE.)
    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE prepare_vert_interp_p


  !-------------
  !>
  !! Compute coefficients for vertical interpolation of model-level
  !! fields to isentropes.
  !!
  SUBROUTINE prepare_vert_interp_i(p_patch, p_prog, p_diag, p_metrics, intp_hrz, nilev,     &
    &                              gh_i_out, temp_i_out, p_i3d_out, vcoeff_i, lacc)

    TYPE(t_patch),        TARGET,      INTENT(IN)    :: p_patch
    TYPE(t_nh_prog),      POINTER                    :: p_prog
    TYPE(t_nh_diag),      POINTER                    :: p_diag
    TYPE(t_nh_metrics),                INTENT(IN)    :: p_metrics
    TYPE(t_int_state),                 INTENT(IN)    :: intp_hrz          !< pointer to data structure for interpolation
    INTEGER,                           INTENT(IN)    :: nilev             !< number of output levels (isentropes)
    REAL(wp),                          INTENT(INOUT) :: gh_i_out(:,:,:), &
      &                                                 temp_i_out(:,:,:)
    REAL(wp),             POINTER                    :: p_i3d_out(:,:,:)  !< theta field
    TYPE(t_vcoeff),                    INTENT(INOUT) :: vcoeff_i          !< out
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    ! LOCAL VARIABLES
    CHARACTER(*), PARAMETER :: routine = modname//":prepare_vert_interp_i"
    INTEGER :: nlev, nlevp1, nblks_c, nblks_e, npromz_c,npromz_e !< blocking parameters
    REAL(wp), DIMENSION(nproma,nilev,p_patch%nblks_e)        :: gh_i_edge
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_me

    vcoeff_i%l_initialized = .TRUE.

    nlev     = p_patch%nlev
    nlevp1   = p_patch%nlevp1
    nblks_c  = p_patch%nblks_c
    npromz_c = p_patch%npromz_c
    nblks_e  = p_patch%nblks_e
    npromz_e = p_patch%npromz_e

    !$ACC DATA CREATE(gh_i_edge, z_me)

    CALL assert_acc_device_only(routine, lacc)

    ! allocate coefficient table:
    CALL vcoeff_allocate(nblks_c, nblks_e, nilev, vcoeff_i, lacc=.TRUE.)

    !--- Coefficients: Interpolation to isentropes

    ! Compute height at isentropic levels (i.e. geopot/g); this height
    ! field is afterwards also used as target coordinate for vertical
    ! interpolation

    ! Note: the following coefficients are identical to the z-level
    !       setup, but they are computed independently to avoid some
    !       complicated code:
    CALL prepare_extrap(p_metrics%z_mc, nblks_c, npromz_c, nlev,                            & !in
      &                 vcoeff_i%lin_cell%kpbl1, vcoeff_i%lin_cell%wfacpbl1,                & !out
      &                 vcoeff_i%lin_cell%kpbl2, vcoeff_i%lin_cell%wfacpbl2,                & !out
      &                 lacc=.TRUE.)

    CALL z_at_theta_levels(p_prog%theta_v, p_metrics%z_mc,                                  & !in
      &                    p_i3d_out, gh_i_out,                                             & !out
      &                    vcoeff_i%lin_cell%wfacpbl1, vcoeff_i%lin_cell%kpbl1,             & !in
      &                    vcoeff_i%lin_cell%wfacpbl2, vcoeff_i%lin_cell%kpbl2, nblks_c,    & !in
      &                    npromz_c, nlev, nilev, lacc=.TRUE.)

    ! Prepare again interpolation coefficients (now for isentropic levels)
    CALL prepare_lin_intp(p_metrics%z_mc, gh_i_out, nblks_c, npromz_c, nlev, nilev,         & !in
      &                   vcoeff_i%lin_cell%wfac_lin, vcoeff_i%lin_cell%idx0_lin,           & !out
      &                   vcoeff_i%lin_cell%bot_idx_lin,                                    & !out
      &                   lacc=.TRUE.)

    CALL prepare_cubic_intp(p_metrics%z_mc, gh_i_out, nblks_c, npromz_c, nlev, nilev,       & !in
      &                     vcoeff_i%cub_cell%coef1, vcoeff_i%cub_cell%coef2,               & !out
      &                     vcoeff_i%cub_cell%coef3,                                        & !out
      &                     vcoeff_i%cub_cell%idx0_cub, vcoeff_i%cub_cell%bot_idx_cub,      & !out
      &                     lacc=.TRUE.)

    ! Perform vertical interpolation
    CALL temperature_intp(p_diag%temp, temp_i_out, p_metrics%z_mc, gh_i_out,                & !in,out
      &                   nblks_c, npromz_c, nlev, nilev,                                   & !in
      &                   vcoeff_i%cub_cell%coef1, vcoeff_i%cub_cell%coef2,                 & !in
      &                   vcoeff_i%cub_cell%coef3,                                          & !in
      &                   vcoeff_i%lin_cell%wfac_lin, vcoeff_i%cub_cell%idx0_cub,           & !in
      &                   vcoeff_i%lin_cell%idx0_lin, vcoeff_i%cub_cell%bot_idx_cub,        & !in
      &                   vcoeff_i%lin_cell%bot_idx_lin, vcoeff_i%lin_cell%wfacpbl1,        & !in
      &                   vcoeff_i%lin_cell%kpbl1, vcoeff_i%lin_cell%wfacpbl2,              & !in
      &                   vcoeff_i%lin_cell%kpbl2, l_restore_sfcinv=.TRUE.,                 & !in
      &                   l_hires_corr=.FALSE., extrapol_dist=0._wp, l_pz_mode=.TRUE.,      & !in
      &                   lacc=.TRUE.)

    !--- Interpolation data for the vertical interface of cells, "nlevp1"

    CALL prepare_lin_intp(p_metrics%z_ifc, gh_i_out,                                        & !in
      &                   nblks_c, npromz_c, nlevp1, nilev,                                 & !in
      &                   vcoeff_i%lin_cell_nlevp1%wfac_lin,                                & !out
      &                   vcoeff_i%lin_cell_nlevp1%idx0_lin,                                & !out
      &                   vcoeff_i%lin_cell_nlevp1%bot_idx_lin,                             & !out
      &                   lacc=.TRUE. )
    CALL prepare_extrap(p_metrics%z_ifc, nblks_c, npromz_c, nlevp1,                         & !in
      &                 vcoeff_i%lin_cell_nlevp1%kpbl1, vcoeff_i%lin_cell_nlevp1%wfacpbl1,  & !out
      &                 vcoeff_i%lin_cell_nlevp1%kpbl2, vcoeff_i%lin_cell_nlevp1%wfacpbl2,  & !out
      &                 lacc=.TRUE.)

    !--- Compute data for interpolation of edge-based fields

    CALL cells2edges_scalar(p_metrics%z_mc, p_patch, intp_hrz%c_lin_e, z_me, opt_fill_latbc=.TRUE., lacc=.TRUE.)
    CALL cells2edges_scalar(gh_i_out, p_patch, intp_hrz%c_lin_e, gh_i_edge, opt_fill_latbc=.TRUE., lacc=.TRUE.)
    CALL sync_patch_array_mult(SYNC_E,p_patch,2,z_me,gh_i_edge)

    CALL prepare_lin_intp(z_me, gh_i_edge, nblks_e, npromz_e, nlev, nilev,                   & !in
      &                   vcoeff_i%lin_edge%wfac_lin, vcoeff_i%lin_edge%idx0_lin,            & !out
      &                   vcoeff_i%lin_edge%bot_idx_lin,                                     & !out
      &                   lacc=.TRUE. )
    CALL prepare_extrap(z_me, nblks_e, npromz_e, nlev,                                       & !in
      &                 vcoeff_i%lin_edge%kpbl1, vcoeff_i%lin_edge%wfacpbl1,                 & !out
      &                 vcoeff_i%lin_edge%kpbl2, vcoeff_i%lin_edge%wfacpbl2,                 & !out
      &                 lacc=.TRUE.)
    CALL prepare_cubic_intp(z_me, gh_i_edge, nblks_e, npromz_e, nlev, nilev,                 & !in
      &                     vcoeff_i%cub_edge%coef1, vcoeff_i%cub_edge%coef2,                & !out
      &                     vcoeff_i%cub_edge%coef3,                                         & !out
      &                     vcoeff_i%cub_edge%idx0_cub, vcoeff_i%cub_edge%bot_idx_cub,       & !out
      &                     lacc=.TRUE. )
    !$ACC WAIT
    !$ACC END DATA
  END SUBROUTINE prepare_vert_interp_i


  !-------------
  !>
  !! SUBROUTINE z_at_ple
  !! Computes height for a given set of pressure levels based on the same extrapolation
  !! assumptions as used in ECMWF's IFS.
  !! The purpose of this is to prepare diagnostic output on pressure levels.
  !!
  !! Required input fields: pressure, virtual temperature and 3D height coordinate
  !! on model levels, pressure of output levels
  !! Output: 3D height field of pressure levels
  !!
  !! Method: piecewise analytical integration of the hydrostatic equation
  !!
  SUBROUTINE z_at_plevels(pres_ml, tempv_ml, z3d_ml, pres_pl, z3d_pl, &
                          nblks, npromz, nlevs_ml, nlevs_pl,          &
                          kpbl1, wfacpbl1, kpbl2, wfacpbl2,           &
                          zextrap, lacc                               )

    ! Input fields
    REAL(wp),         INTENT(IN)  :: pres_ml  (:,:,:) ! pressure field on model levels
    REAL(wp),         INTENT(IN)  :: tempv_ml (:,:,:) ! virtual temperature on model levels
    REAL(wp), TARGET, INTENT(IN)  :: z3d_ml   (:,:,:) ! 3D height coordinate field on model levels
    REAL(wp),         INTENT(IN)  :: pres_pl  (:,:,:) ! pressure of output levels

    ! Comment: for z3d_zl and pres_pl, 1D fields would actually be sufficient, but having
    ! everything as 3D fields simplifies programming

    ! Output
    REAL(wp), INTENT(OUT) :: z3d_pl  (:,:,:) ! 3D height coordinate field of pressure levels

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_ml   ! Number of model levels
    INTEGER , INTENT(IN) :: nlevs_pl   ! Number of pressure levels

    ! Coefficients needed for removing the surface inversion before extrapolating downward
    INTEGER , INTENT(IN) :: kpbl1(:,:)     ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)  ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)     ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)  ! corresponding interpolation coefficient

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! if true use OpenACC

    REAL(wp), OPTIONAL, INTENT(IN) :: zextrap(:,:)   ! AGL height from which downward extrapolation starts (in postprocesing mode)

    ! LOCAL VARIABLES

    CHARACTER(*), PARAMETER :: routine = modname//":z_at_plevels"
    INTEGER  :: jb, jkm, jkp, jc, jkm_start, jkp_start
    INTEGER  :: nlen
    INTEGER , DIMENSION(nproma)          :: bot_idx_ml
    INTEGER , DIMENSION(nproma,nlevs_pl) :: idx0_ml

    REAL(wp) :: z_up, z_down
    REAL(wp), DIMENSION(nproma,nlevs_pl) :: wfac_ml, dtvdz
    ! temporary to store pre-computed inverse of differences
    REAL(wp), DIMENSION(nproma,nlevs_ml) :: z3d_ml_di
    REAL(wp), DIMENSION(nproma)          :: tmsl, tsfc_mod, tempv1, tempv2, vtgrad_up, sfc_inv
    REAL(wp), ALLOCATABLE, TARGET        :: zgpot_ml(:,:,:)
    REAL(wp),              POINTER       :: z_ml(:,:,:)

    LOGICAL :: l_found(nproma),lfound_all,lzextrap, lerror
    INTEGER :: istat
    LOGICAL :: lzacc ! non-optional version of lacc

!-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (PRESENT(zextrap)) THEN
      lzextrap = .TRUE.
    ELSE
      lzextrap = .FALSE.
    ENDIF

    IF (.NOT. ldeepatmo) THEN
      z_ml => z3d_ml
    ELSE
      ALLOCATE(zgpot_ml(nproma, nlevs_ml, nblks), STAT=istat)
      IF (istat /= SUCCESS) CALL finish(routine, 'Allocation of zgpot failed')
      !$ACC ENTER DATA &
      !$ACC   CREATE(zgpot_ml) &
      !$ACC   IF(lzacc)
      CALL deepatmo_htrafo(z_in=z3d_ml, z_out=zgpot_ml, nblks_nproma_npromz=[nblks,nproma,npromz], &
        & start_end_levels=[1,nlevs_ml], radius=grid_sphere_radius, trafo_type='z2zgpot', ierror=istat, lacc=lzacc)
      IF (istat /= SUCCESS) CALL finish(routine, 'deepatmo_htrafo failed')
      z_ml => zgpot_ml
      ! Note: the heights above ground level, 'zpbl1', 'zpbl2', 'zextrap' and heights derived from them 
      ! have relatively low values (~ 1 km), so no deep-atmosphere modification is applied to them. 
      ! (Put another way: we regard 'zpbl1', 'zpbl2', and 'zextrap' to represent geopotential heights.)
    ENDIF

    !$ACC DATA IF(lzacc) &
    !$ACC   PRESENT(pres_ml, tempv_ml, z3d_ml, pres_pl, z3d_pl, kpbl1, wfacpbl1, kpbl2, wfacpbl2, zextrap, vct_a) &
    !$ACC   PRESENT(num_lev, z_ml) &
    !$ACC   CREATE(bot_idx_ml, idx0_ml, wfac_ml, dtvdz, z3d_ml_di, tmsl, tsfc_mod, tempv1, tempv2, vtgrad_up) &
    !$ACC   CREATE(sfc_inv, l_found)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jkm,jkp,jc,nlen,jkm_start,jkp_start,bot_idx_ml,idx0_ml,z_up,z_down,&
!$OMP wfac_ml,tmsl,tsfc_mod,dtvdz,tempv1,tempv2,vtgrad_up,sfc_inv,l_found,&
!$OMP lfound_all,z3d_ml_di,lerror) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jkp = 1, nlevs_pl
          DO jc = nlen+1, nproma
            z3d_pl(jc,jkp,jb) = 0.0_wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF
      lerror = .FALSE.

      ! First compute coefficients for those target levels/points for which interpolation
      ! between model-level data is possible
      jkm_start = 1
      DO jkp = 1, nlevs_pl
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          l_found(jc) = .FALSE.
        ENDDO
        lfound_all = .FALSE.
        !$ACC LOOP SEQ
        DO jkm = jkm_start,nlevs_ml-1
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen
            IF( .NOT. l_found(jc)) THEN
              IF (pres_pl(jc,jkp,jb) >= pres_ml(jc,jkm,jb) .AND. &
                  pres_pl(jc,jkp,jb) <  pres_ml(jc,jkm+1,jb)) THEN
                idx0_ml(jc,jkp) = jkm
                wfac_ml(jc,jkp) = (pres_pl(jc,jkp,jb)-pres_ml(jc,jkm+1,jb))/&
                                  (pres_ml(jc,jkm,jb)-pres_ml(jc,jkm+1,jb))
                bot_idx_ml(jc) = jkp
                l_found(jc) = .TRUE.
              ELSE IF (pres_pl(jc,jkp,jb) >= pres_ml(jc,nlevs_ml,jb)) THEN
                l_found(jc) = .TRUE.
                idx0_ml(jc,jkp) = nlevs_ml
                IF (jkp == 1) bot_idx_ml(jc) = 0
              ELSE IF (pres_pl(jc,jkp,jb) < pres_ml(jc,1,jb)) THEN
                idx0_ml(jc,jkp) = 1
                wfac_ml(jc,jkp) = (pres_pl(jc,jkp,jb)-pres_ml(jc,2,jb))/&
                                  (pres_ml(jc,1  ,jb)-pres_ml(jc,2,jb))
                bot_idx_ml(jc) = jkp
                l_found(jc) = .TRUE.
              ENDIF
            ENDIF
          ENDDO
          ! this optimization is not supported by OpenACC but also less important
          IF (ALL(l_found(1:nlen))) THEN
            lfound_all = .TRUE.
            EXIT
          ENDIF
        ENDDO
        !$ACC END PARALLEL








        IF (lfound_all) THEN
          jkm_start = MIN(minval_1d(idx0_ml(1:nlen, jkp), lacc=lzacc), nlevs_ml-1)
        ELSE
          lerror = .TRUE.
          EXIT
        ENDIF
      ENDDO

      IF (lerror) CALL finish(routine, "Error in computing interpolation coefficients")

      IF (itype_pres_msl == PRES_MSL_METHOD_GME) THEN

        ! Preparations for extrapolation below the ground: calculate temperature
        ! profile with artificial limits like in IFS. This temperature profile
        ! is NOT consistent with the temperature on pressure levels computed
        ! afterwards, implying that temperature and geopotential on pressure
        ! levels are not in hydrostatic balance. We are not happy with this
        ! solution, but it is used for the time being in order to be consistent
        ! with IFS and GME.

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = 1, nlen
          tsfc_mod(jc) = tempv_ml(jc,nlevs_ml,jb)
          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z_ml(jc,nlevs_ml,jb)
          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF
        ENDDO
        !$ACC END PARALLEL

      ELSE IF ( lzextrap .AND. itype_pres_msl >= 3 ) THEN

        ! Similar method to option 1, but we use the temperature at 150 m AGL
        ! for extrapolation (kpbl1 contains the required index in this case)

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = 1, nlen

          tsfc_mod(jc) = wfacpbl1(jc,jb) *tempv_ml(jc,kpbl1(jc,jb),jb  ) + &
                  (1._wp-wfacpbl1(jc,jb))*tempv_ml(jc,kpbl1(jc,jb)+1,jb) - &
                  dtdz_standardatm*(zextrap(jc,jb)-0.5_wp*vct_a(num_lev(1)))

          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z_ml(jc,nlevs_ml,jb)
          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF
        ENDDO
        !$ACC END PARALLEL

      ELSE
        ! Use a temperature profile that is (roughly) consistent with
        ! that used for temperature extrapolation

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = 1, nlen
          ! Virtual temperature at height zpbl1
          tempv1(jc) = wfacpbl1(jc,jb) *tempv_ml(jc,kpbl1(jc,jb),jb  ) + &
                (1._wp-wfacpbl1(jc,jb))*tempv_ml(jc,kpbl1(jc,jb)+1,jb)

          ! Virtual temperature at height zpbl2
          tempv2(jc) = wfacpbl2(jc,jb) *tempv_ml(jc,kpbl2(jc,jb),jb  ) + &
                (1._wp-wfacpbl2(jc,jb))*tempv_ml(jc,kpbl2(jc,jb)+1,jb)

          ! Vertical gradient between zpbl1 and zpbl2
          vtgrad_up(jc) = (tempv2(jc) - tempv1(jc))/(zpbl2 - zpbl1)

          ! Modified (extrapolated) surface temperature without inversion
          tsfc_mod(jc) = tempv1(jc) - zpbl1*vtgrad_up(jc)

          ! "surface inversion", defined by the difference between the extrapolated
          ! extrapolated temperature from above and the original input temperature
          sfc_inv(jc) = tsfc_mod(jc) - tempv_ml(jc,nlevs_ml,jb)

          ! Reduction of the surface inversion depending on the extrapolation
          ! distance. The surface inversion is fully restored for extrapolation distances
          ! up to zpbl1 and disregarded for distances larger than 3*zpbl1
          IF (z_ml(jc,nlevs_ml,jb) > 3._wp*zpbl1) THEN
            sfc_inv(jc) = 0._wp
          ELSE IF (z_ml(jc,nlevs_ml,jb) > zpbl1) THEN
            sfc_inv(jc) = sfc_inv(jc)*(1._wp - (z_ml(jc,nlevs_ml,jb)-zpbl1)/(2._wp*zpbl1))
          ENDIF

          tsfc_mod(jc) = tsfc_mod(jc) - sfc_inv(jc)

          ! Limitation of very cold temperatures according to GME method
          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))

          ! Estimated temperature at mean sea level
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z_ml(jc,nlevs_ml,jb)

          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF

        ENDDO
        !$ACC END PARALLEL

      ENDIF

      jkp_start = minval_1d(bot_idx_ml(1: nlen), lacc=lzacc) + 1


      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP SEQ
      DO jkp = jkp_start, nlevs_pl
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          IF (jkp >= bot_idx_ml(jc)+1) THEN
            ! Store extrapolation distance on wfac_ml if target point is below the
            ! surface of the input data (note: this is a positive quantity)
            idx0_ml(jc,jkp) = nlevs_ml
            wfac_ml(jc,jkp) = pres_pl(jc,jkp,jb)-pres_ml(jc,nlevs_ml,jb)
          ENDIF
        ENDDO
      ENDDO

      !$ACC LOOP SEQ
      DO jkm = 1, nlevs_ml - 1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          z3d_ml_di(jc, jkm) = 1.0_wp &
               / (z_ml(jc,jkm,jb) - z_ml(jc,jkm+1,jb))
        END DO
      END DO

      ! Now, compute vertical gradients of virtual potential temperature
      !$ACC LOOP SEQ
      DO jkp = 1, nlevs_pl
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jkm)
        DO jc = 1, nlen
          IF (jkp <= bot_idx_ml(jc)) THEN
            jkm = idx0_ml(jc,jkp)
            dtvdz(jc,jkp) = (tempv_ml(jc,jkm,jb)-tempv_ml(jc,jkm+1,jb)) &
                 * z3d_ml_di(jc, jkm)
          ELSE IF (z_ml(jc,nlevs_ml,jb) > 1._wp) THEN ! extrapolation below lowest model level required
            dtvdz(jc,jkp) = (tsfc_mod(jc)-tmsl(jc))/z_ml(jc,nlevs_ml,jb)

          ELSE ! avoid pathological results at grid points below sea level
            dtvdz(jc,jkp) = dtdz_standardatm
          ENDIF
        ENDDO
      ENDDO

      ! Finally, compute height of pressure levels
      !$ACC LOOP SEQ
      DO jkp = 1, nlevs_pl
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jkm, z_up, z_down)
        DO jc = 1, nlen
          IF (jkp <= bot_idx_ml(jc)) THEN

            ! interpolation based on piecewise analytical integration of the hydrostatic equation
            jkm = idx0_ml(jc,jkp)
            IF (ABS(dtvdz(jc,jkp)) > dtvdz_thresh) THEN
              z_up   = z_ml(jc,jkm,jb) + tempv_ml(jc,jkm,jb)*((pres_pl(jc,jkp,jb) /     &
                       pres_ml(jc,jkm,jb))**(-rd*dtvdz(jc,jkp)/grav)-1._wp)/dtvdz(jc,jkp)
              z_down = z_ml(jc,jkm+1,jb) + tempv_ml(jc,jkm+1,jb)*((pres_pl(jc,jkp,jb) /   &
                       pres_ml(jc,jkm+1,jb))**(-rd*dtvdz(jc,jkp)/grav)-1._wp)/dtvdz(jc,jkp)

            ELSE
              z_up   = z_ml(jc,jkm,jb) + (rd*0.5_wp*(tempv_ml(jc,jkm,jb) +                   &
                       tempv_ml(jc,jkm+1,jb))/grav)*LOG(pres_ml(jc,jkm,jb)/pres_pl(jc,jkp,jb))
              z_down = z_ml(jc,jkm+1,jb) + (rd*0.5_wp*(tempv_ml(jc,jkm,jb) +                   &
                       tempv_ml(jc,jkm+1,jb))/grav)*LOG(pres_ml(jc,jkm+1,jb)/pres_pl(jc,jkp,jb))
            ENDIF

            ! Finally, apply inverse-distance weighting between top-down and bottom-up integrated value
            ! except in the case of extrapolation above the top of the input data
            IF (jkm == 1 .AND. wfac_ml(jc,jkp) > 1._wp) THEN
              z3d_pl(jc,jkp,jb) = z_up
            ELSE
              z3d_pl(jc,jkp,jb) = wfac_ml(jc,jkp)*z_up + (1._wp-wfac_ml(jc,jkp))*z_down
            ENDIF

          ELSE !  extrapolation below lowest height level

            jkm = nlevs_ml
            IF (ABS(dtvdz(jc,jkp)) > dtvdz_thresh) THEN
              z_up   = z_ml(jc,jkm,jb) + tsfc_mod(jc)*((pres_pl(jc,jkp,jb) /            &
                       pres_ml(jc,jkm,jb))**(-rd*dtvdz(jc,jkp)/grav)-1._wp)/dtvdz(jc,jkp)
            ELSE
              z_up   = z_ml(jc,jkm,jb) + (rd*tsfc_mod(jc)/grav) * &
                       LOG(pres_ml(jc,jkm,jb)/pres_pl(jc,jkp,jb))
            ENDIF

            z3d_pl(jc,jkp,jb) = z_up

          ENDIF

        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT
    !$ACC END DATA
    NULLIFY(z_ml)
    IF (ldeepatmo) THEN
      !$ACC EXIT DATA &
      !$ACC   DELETE(zgpot_ml) &
      !$ACC   IF(lzacc)
      DEALLOCATE(zgpot_ml, STAT=istat)
      IF (istat /= SUCCESS) CALL finish(routine, 'Deallocation of zgpot failed')
      ! 'z3d_pl' contains geopotential heights: transform it to geometric heights
      CALL deepatmo_htrafo(z_inout=z3d_pl, nblks_nproma_npromz=[nblks,nproma,npromz], &
        & start_end_levels=[1,nlevs_pl], radius=grid_sphere_radius,trafo_type='zgpot2z', ierror=istat, lacc=lzacc)
      IF (istat /= SUCCESS) CALL finish(routine, "deepatmo_htrafo failed.")
    ENDIF

  END SUBROUTINE z_at_plevels

  !-------------
  !>
  !! SUBROUTINE z_at_theta_levels
  !! Computes height for a given set of potential temperature levels
  !! The purpose of this is to prepare diagnostic output on isentropic levels.
  !!
  !! Required input fields: potential temperature and 3D height coordinate
  !! on model levels, coefficients for computing the boundary layer temperature gradients
  !! Output: 3D height field of potential temperature levels
  !!
  !! Method: linear vertical interpolation / extrapolation
  !!
  !! Note: this routine should not be used for large extrapolation distances below
  !! the model orography. The target values of the interpolation (theta_thl)
  !! must be in descending order. It the input theta field does not increase monotonically
  !! with height, the first occurrence of the target value - scanning from top to bottom - is selected
  !!
  SUBROUTINE z_at_theta_levels(theta_ml, z3d_ml, theta_thl, z3d_thl, wfacpbl1, kpbl1, &
                               wfacpbl2, kpbl2, nblks, npromz, nlevs_ml, nlevs_thl, lacc)


    ! Input fields
    REAL(wp), INTENT(IN)  :: theta_ml (:,:,:) ! potential temperature on model levels
    REAL(wp), INTENT(IN)  :: z3d_ml   (:,:,:) ! 3D height coordinate field on model levels
    REAL(wp), INTENT(IN)  :: theta_thl(:,:,:) ! Potential temperature on theta levels
                                              ! (target field of interpolation)

    ! Comment: for theta_thl, a 1D field would actually be sufficient, but having
    ! everything as 3D fields simplifies programming

    ! Input coefficients
    INTEGER , INTENT(IN) :: kpbl1(:,:)      ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)   ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)      ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)   ! corresponding interpolation coefficient

    ! Output
    REAL(wp), INTENT(OUT)  :: z3d_thl  (:,:,:) ! 3D height coordinate field on theta levels

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_ml   ! Number of model levels
    INTEGER , INTENT(IN) :: nlevs_thl  ! Number of theta levels
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc
    ! LOCAL VARIABLES


    INTEGER  :: jb, jkm, jkt, jc, jkm_start
    INTEGER  :: nlen
    INTEGER , DIMENSION(nproma)          :: bot_idx_ml
    INTEGER , DIMENSION(nproma,nlevs_thl) :: idx0_ml

    REAL(wp), DIMENSION(nproma,nlevs_thl) :: wfac_ml

    REAL(wp) :: theta1, theta2, vthgrad

    LOGICAL :: l_found(nproma),lfound_all, lerror

!-------------------------------------------------------------------------
    CALL assert_acc_device_only("z_at_theta_levels", lacc)

    !$ACC DATA PRESENT(theta_ml, z3d_ml, theta_thl, kpbl1, wfacpbl1, kpbl2, wfacpbl2, z3d_thl) &
    !$ACC   CREATE(bot_idx_ml, idx0_ml, wfac_ml, l_found)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jkm,jkt,jc,nlen,jkm_start,bot_idx_ml,idx0_ml,wfac_ml,theta1,theta2,vthgrad, &
!$OMP            l_found,lfound_all, lerror) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1)
        DO jkt = 1, nlevs_thl
          DO jc = nlen+1, nproma
            z3d_thl(jc,jkt,jb) = 0.0_wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF
      lerror = .FALSE.

      ! Compute coefficients for those target levels/points for which interpolation
      ! between model-level data is possible
      jkm_start = 1
      !$ACC PARALLEL LOOP ASYNC(1) DEFAULT(PRESENT)
      DO jc = 1, nlen
        bot_idx_ml(jc) = 0
      ENDDO
      !$ACC END PARALLEL LOOP
      DO jkt = 1, nlevs_thl
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          l_found(jc) = .FALSE.
        ENDDO
        lfound_all = .FALSE.
        !$ACC LOOP SEQ
        DO jkm = jkm_start,nlevs_ml-1
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen
            IF (.NOT. l_found(jc)) THEN
              IF (theta_thl(jc,jkt,jb) <= theta_ml(jc,jkm,jb) .AND. &
                  theta_thl(jc,jkt,jb) >  theta_ml(jc,jkm+1,jb)) THEN
                idx0_ml(jc,jkt) = jkm
                wfac_ml(jc,jkt) = (theta_thl(jc,jkt,jb)-theta_ml(jc,jkm+1,jb))/&
                                  (theta_ml(jc,jkm,jb) -theta_ml(jc,jkm+1,jb))
                bot_idx_ml(jc) = jkt
                l_found(jc) = .TRUE.
              ! The second criterion here is needed to prevent running into the extrapolation branch
              ! if a pair of levels fulfiling the interpolation condition has already been found
              ! (relevant if theta does not increase monotonically with height)
              ELSE IF (theta_thl(jc,jkt,jb) <= theta_ml(jc,nlevs_ml,jb) .AND. bot_idx_ml(jc) < jkt ) THEN
                l_found(jc) = .TRUE.
                idx0_ml(jc,jkt) = nlevs_ml
                ! Store "extrapolation distance" on wfac_ml if target point is below the
                ! surface of the input data
                wfac_ml(jc,jkt) = theta_thl(jc,jkt,jb)-theta_ml(jc,nlevs_ml,jb)
              ELSE IF (theta_thl(jc,jkt,jb) > theta_ml(jc,1,jb)) THEN
                idx0_ml(jc,jkt) = 1
                wfac_ml(jc,jkt) = (theta_thl(jc,jkt,jb)-theta_ml(jc,2,jb))/&
                                  (theta_ml (jc,1  ,jb)-theta_ml(jc,2,jb))
                bot_idx_ml(jc) = jkt
                l_found(jc) = .TRUE.
              ENDIF
            ENDIF
          ENDDO

          ! this optimization is not supported by OpenACC but also less important
          IF (ALL(l_found(1:nlen))) THEN
            lfound_all = .TRUE.
            EXIT
          ENDIF

        ENDDO
        !$ACC END PARALLEL
        IF (lfound_all) THEN
          jkm_start = MIN(minval_1d(idx0_ml(1:nlen, jkt), lacc=.TRUE.), nlevs_ml-1)
        ELSE
          lerror = .TRUE.
          EXIT
        ENDIF
      ENDDO

      IF (lerror) CALL finish("z_at_theta_levels:",&
        "Calculation of interpolation coefficients failed")

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP SEQ
      DO jkt = 1, nlevs_thl
        !$ACC LOOP GANG VECTOR PRIVATE(theta1, theta2, vthgrad)
        DO jc = 1, nlen
          IF (jkt <= bot_idx_ml(jc)) THEN

            ! linear interpolation
            z3d_thl(jc,jkt,jb) = wfac_ml(jc,jkt)*z3d_ml(jc,idx0_ml(jc,jkt),jb) + &
              (1._wp-wfac_ml(jc,jkt))*z3d_ml(jc,idx0_ml(jc,jkt)+1,jb)

          ELSE

            ! linear extrapolation, using the gradient between (by default) 1000 m and 500 m AGL

            ! Potential temperature at height zpbl1
            theta1 = wfacpbl1(jc,jb) *theta_ml(jc,kpbl1(jc,jb),jb  ) + &
              (1._wp-wfacpbl1(jc,jb))*theta_ml(jc,kpbl1(jc,jb)+1,jb)

            ! Potential temperature at height zpbl2
            theta2 = wfacpbl2(jc,jb) *theta_ml(jc,kpbl2(jc,jb),jb  ) + &
              (1._wp-wfacpbl2(jc,jb))*theta_ml(jc,kpbl2(jc,jb)+1,jb)

            ! Vertical gradient between zpbl1 and zpbl2
            vthgrad = (theta2 - theta1)/(zpbl2 - zpbl1)

            ! Set reasonable limits
            vthgrad = MAX(vthgrad,3.0e-3_wp)
            vthgrad = MIN(vthgrad,8.5e-3_wp)

            ! wfac carries the theta difference in this case
            z3d_thl(jc,jkt,jb) = z3d_ml(jc,nlevs_ml,jb) + wfac_ml(jc,jkt)/vthgrad

          ENDIF

        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE z_at_theta_levels

END MODULE mo_nh_vert_interp_ipz
