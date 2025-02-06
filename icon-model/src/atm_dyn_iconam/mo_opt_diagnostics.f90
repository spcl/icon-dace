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

! Definition of optional (diagnostic) model variables
!
! In the metadata of each model variable ("add_var") it is specified
! *how* certain post-processing tasks, e.g. vertical interpolation
! onto p/z-levels, are treated. In the namelist, users can then
! specify *if* computations for these variables are performed.
!
! If so, the resulting model fields are appended to the list of
! internal post-processing tasks (each field forms its own task). As
! we do not know in advance the contents of this list, we call them
! "optional diagnostics".

MODULE mo_opt_diagnostics

  USE mo_kind,                 ONLY: wp
  USE mo_parallel_config,      ONLY: nproma
  USE mo_master_control,       ONLY: get_my_process_name
  USE mo_model_domain,         ONLY: t_patch, t_subset_range
  USE mo_nonhydro_types,       ONLY: t_nh_diag,t_nh_prog,      &
                                     t_nh_state_lists
  USE mo_impl_constants,       ONLY: success,     &
    &                                VINTP_METHOD_QV,                    &
    &                                VINTP_METHOD_PRES,                  &
    &                                VINTP_METHOD_LIN,                   &
    &                                VINTP_METHOD_LIN_NLEVP1,            &
    &                                HINTP_TYPE_NONE
  USE mo_physical_constants,   ONLY: earth_radius
  USE mo_exception,            ONLY: finish!!$, message, message_text
  USE mo_fortran_tools,        ONLY: init, assert_acc_device_only
  USE mo_grid_config,          ONLY: n_dom
  USE mo_run_config,           ONLY: ntracer,iqv,iqc,iqi
  USE mo_advection_config,     ONLY: t_advection_config, advection_config
  USE mo_zaxis_type,           ONLY: ZA_REFERENCE, ZA_REFERENCE_HALF, ZA_SURFACE, &
    &                                ZA_MEANSEA
  USE mo_cdi,                  ONLY: DATATYPE_FLT32, DATATYPE_PACK16,                  &
    &                                DATATYPE_PACK24,                                  &
    &                                DATATYPE_FLT64, GRID_UNSTRUCTURED,                &
    &                                TSTEP_CONSTANT
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_CELL,                           &
    &                                GRID_CELL, GRID_REGULAR_LONLAT, LONLAT_PREFIX
  USE mo_var_list,             ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var_list_register,    ONLY: vlr_add, vlr_del
  USE mo_var, ONLY: level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_name_list_output_config, ONLY: is_variable_in_output
  USE mo_io_config,            ONLY: lnetcdf_flt64_output
  USE mo_gribout_config,       ONLY: gribout_config
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var, grib2_var
  USE mo_var_groups,           ONLY: groups
  USE mo_var_metadata,         ONLY: create_vert_interp_metadata,                      &
    &                                create_hor_interp_metadata,                       &
    &                                vintp_types
  USE mo_var,                  ONLY: t_var
  USE mo_tracer_metadata,      ONLY: create_tracer_metadata
  USE mo_statistics,           ONLY: add_fields
  USE mo_util_dbg_prnt,        ONLY: dbg_print
  USE mo_lonlat_grid,          ONLY: t_lon_lat_grid, latlon_compute_area_weights
  USE mo_intp_lonlat_types,    ONLY: t_lon_lat_intp, t_lon_lat_list

  IMPLICIT NONE

  PRIVATE


  ! data types
  PUBLIC :: t_nh_opt_diag         ! optional diagnostic variables (data type)
  PUBLIC :: t_nh_acc
  PUBLIC :: p_nh_opt_diag         ! state vector of optional diagnostic variables
                                  ! e.g. variables on p- and/or z-levels
  PUBLIC :: t_nh_diag_pz
  PUBLIC :: t_vcoeff, t_vcoeff_lin, t_vcoeff_cub
  ! subroutines
  PUBLIC :: vcoeff_allocate, vcoeff_deallocate
  PUBLIC :: construct_opt_diag
  PUBLIC :: destruct_opt_diag
  PUBLIC :: update_opt_acc, reset_opt_acc, calc_mean_opt_acc
  PUBLIC :: compute_lonlat_area_weights

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_opt_diagnostics'

  ! Sub-type of "t_vcoeff" containing linear interpolation
  ! coefficients
  TYPE t_vcoeff_lin
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: wfac_lin            ! (nproma, nlev, nblks)
    INTEGER,  ALLOCATABLE, DIMENSION(:,:,:) :: idx0_lin            ! (nproma, nlev, nblks)
    INTEGER,  ALLOCATABLE, DIMENSION(:,:)   :: bot_idx_lin         ! (nproma, nblks)
    INTEGER,  ALLOCATABLE, DIMENSION(:,:)   :: kpbl1, kpbl2        ! (nproma, nblks)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: wfacpbl1, wfacpbl2  ! (nproma, nblks)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: zextrap             ! (nproma, nblks)
  END TYPE t_vcoeff_lin

  ! Sub-type of "t_vcoeff" containing cubic interpolation
  ! coefficients
  TYPE t_vcoeff_cub
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: coef1, coef2, coef3 ! (nproma, nlev, nblks)
    INTEGER,  ALLOCATABLE, DIMENSION(:,:,:) :: idx0_cub            ! (nproma, nlev, nblks)
    INTEGER,  ALLOCATABLE, DIMENSION(:,:)   :: bot_idx_cub         ! (nproma, nblks)
  END TYPE t_vcoeff_cub

  ! Derived type containing coefficient tables for vertical
  ! interpolation.
  TYPE t_vcoeff
    LOGICAL :: &
      & l_initialized = .FALSE.,  &
      & l_allocated   = .FALSE.

    ! LINEAR interpolation data
    TYPE (t_vcoeff_lin) ::    &
      &    lin_cell,          &  !< cell centers: interpolation data for the model levels
      &    lin_cell_nlevp1,   &  !< cell centers: interpolation data for the vertical interface of cells, "nlevp1"
      &    lin_edge              !< edge midpts:  interpolation data for the model levels

    ! CUBIC interpolation (model levels)
    TYPE (t_vcoeff_cub) ::    &
      &    cub_cell,          &  !< cell centers: interpolation data for the model levels
      &    cub_edge

  END TYPE t_vcoeff

  TYPE t_pointer_3d_wp
    REAL(wp),POINTER :: p(:,:,:)  ! pointer to 3D array
  END TYPE t_pointer_3d_wp
  ! variable to be accumulated manually
  TYPE t_nh_acc
    REAL(wp), POINTER, CONTIGUOUS :: &
    !
    ! dynamics
    &  rho(:,:,:),      &
    &  qv(:,:,:),       &
    &  qc(:,:,:),       &
    &  qi(:,:,:),       &
    &  temp(:,:,:),     &
    &  pres_sfc(:,:),   &
    &  pres_msl(:,:),   &
    &  pres(:,:,:),     &
    &  pres_ifc(:,:,:), &
    &  u(:,:,:),        &
    &  v(:,:,:),        &
    &  w(:,:,:),        &
    &  omega(:,:,:),    &
    !
    ! tracers container
    &  tracer(:,:,:,:)

    TYPE(t_pointer_3d_wp),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer

    ! Internal counter for accumulation operations
    INTEGER :: numberOfAccumulations

    ! logicals for presence of time mean output variables in the output name lists
    !
    !  inidcate if any time averaged variable is requested for the output
    LOGICAL :: l_any_m
    !
    !  dynamics
    LOGICAL :: l_ua_m
    LOGICAL :: l_va_m
    LOGICAL :: l_wa_m
    LOGICAL :: l_rho_m
    LOGICAL :: l_ta_m
    LOGICAL :: l_ps_m
    LOGICAL :: l_psl_m
    LOGICAL :: l_pfull_m
    LOGICAL :: l_phalf_m
    LOGICAL :: l_wap_m
    !
    !  tracers
    LOGICAL :: l_tracer_m

  END TYPE t_nh_acc


  ! State vector for diagnostic variables on p-, z- and/or i-levels
  !
  ! @note The pointers which are collected in this derived type
  !       constitute only the minimum set of fields that are required
  !       for i/p/z-level interpolation. All other variables are
  !       stored inside the "opt_diag_list_p", "opt_diag_list_z"
  !       variable lists.
  TYPE t_nh_diag_pz

    REAL(wp), POINTER ::    &
      !--- cells (nproma,nlev,nblks)
      ! fields that are essential for z-level interpolation:
      &  z_temp(:,:,:),        & ! temperature                  [K]
      &  z_pres(:,:,:),        & ! pressure                     [Pa]
      ! fields that are essential for p-level interpolation only:
      &  p_gh    (:,:,:),      & ! geopotential height          [m]
      &  p_temp(:,:,:),        & ! temperature                  [K]
      ! fields that are essential for interpolation on isentropes only:
      &  i_gh    (:,:,:),      & ! geopotential height          [m]
      &  i_temp(:,:,:)           ! temperature                  [K]

    ! coefficient tables for vertical interpolation. There exist
    ! different kinds of coefficients: For p-, z-,and for
    ! i-level-interpolation; for cells and for edges.
    TYPE(t_vcoeff) :: vcoeff_z, vcoeff_p, vcoeff_i

  END TYPE t_nh_diag_pz


  ! List of optional diagnostics + necessary meta data
  TYPE t_nh_opt_diag

    ! diag_pz: data structure containing coefficient tables and
    ! pointers to a few number of fields which are required for
    ! interpolation of model variables to p/z-levels
    TYPE(t_nh_diag_pz) :: diag_pz

    TYPE(t_nh_acc)     :: acc

    ! opt_diag_list: List of optional diagnostics variables.
    !
    ! The "opt_diag_list_*" lists contain all variables that have been
    ! interpolated onto p/z-levels
    TYPE(t_var_list_ptr)   :: opt_diag_list,   opt_diag_list_p, &
      &                   opt_diag_list_z, opt_diag_list_i, &
      &                   opt_acc_list

  END TYPE t_nh_opt_diag


  ! Actual instantiation of optional diagnostics type "t_nh_opt_diag"
  TYPE(t_nh_opt_diag), TARGET, ALLOCATABLE :: p_nh_opt_diag(:)


CONTAINS

  ! setup of accumulation variables
  SUBROUTINE construct_opt_acc(p_patch,list,p_acc)
    TYPE(t_patch),        INTENT(IN) :: p_patch
    TYPE(t_var_list_ptr)                 :: list
    TYPE(t_nh_acc)                   :: p_acc

    ! LOCAL ===================================================================
    INTEGER :: nblks_c       !< number of cell blocks to allocate
!!$    INTEGER :: nblks_e       !< number of edge blocks to allocate
!!$    INTEGER :: nblks_v       !< number of vertex blocks to allocate

    INTEGER :: nlev
    INTEGER :: nlevp1

    INTEGER :: jt

    INTEGER :: shape2d  (2)
    INTEGER :: shape2d_c(2), shape3d_c(3), shape3d_chalf(3), shape4d_c(4)
!!$    INTEGER :: shape2d_e(2), shape3d_e(3)
!!$    INTEGER ::               shape3d_v(3)

    INTEGER :: ibits,iextbits     !< "entropy" of horizontal slice
    INTEGER :: DATATYPE_PACK_VAR  !< variable "entropy" for some thermodynamic fields
    INTEGER :: datatype_flt       !< floating point accuracy in NetCDF output

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc
    TYPE(t_advection_config), POINTER :: advconf
    ! =========================================================================

    !determine size of arrays
    nblks_c = p_patch%nblks_c
!!$    nblks_e = p_patch%nblks_e
!!$    nblks_v = p_patch%nblks_v

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice
    iextbits = DATATYPE_PACK24

    IF (gribout_config(p_patch%id)%lgribout_24bit) THEN  ! analysis
      ! higher accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK24
    ELSE
      ! standard accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK16
    ENDIF

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ! pointer to advection_config(jg) to save some paperwork
    advconf => advection_config(p_patch%id)

    ! predefined array shapes
    shape2d_c     = (/nproma,          nblks_c    /)
    shape2d       = shape2d_c
    shape3d_c     = (/nproma, nlev   , nblks_c    /)
    shape3d_chalf = (/nproma, nlevp1 , nblks_c    /)
    shape4d_c     = (/nproma, nlev   , nblks_c, ntracer     /)
!!$    shape2d_e     = (/nproma,          nblks_e    /)
!!$    shape3d_e     = (/nproma, nlev   , nblks_e    /)
!!$    shape3d_v     = (/nproma, nlev   , nblks_v    /)

    p_acc%l_any_m = .FALSE.

    ! PROGS {{{
    p_acc%l_ua_m  = is_variable_in_output(var_name="ua_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_ua_m
    IF (p_acc%l_ua_m) THEN
       cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'Zonal wind (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( list, 'ua_m', p_acc%u,                                        &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                   & ldims=shape3d_c,                                              &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &   vert_intp_type=vintp_types("P","Z","I") ),                  &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_va_m  = is_variable_in_output(var_name="va_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_va_m
    IF (p_acc%l_va_m) THEN
       cf_desc    = t_cf_var('northward_wind', 'm s-1', 'Meridional wind (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( list, 'va_m', p_acc%v,                                        &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                   & ldims=shape3d_c,                                              &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &   vert_intp_type=vintp_types("P","Z","I") ),                  &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_wa_m  = is_variable_in_output(var_name="wa_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_wa_m
    IF (p_acc%l_wa_m) THEN
       cf_desc    = t_cf_var('upward_air_velocity', 'm s-1', 'Vertical velocity (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 2, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( list, 'wa_m', p_acc%w,                                        &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                   & ldims=shape3d_chalf,                                          &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &   vert_intp_type=vintp_types("P","Z","I"),                    &
                   &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),                 &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_rho_m = is_variable_in_output(var_name="rho_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_rho_m
    IF (p_acc%l_rho_m) THEN
       cf_desc    = t_cf_var('air_density', 'kg m-3', 'density (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 3, 10, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( list, 'rho_m', p_acc%rho,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                   & ldims=shape3d_c,                                              &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &   vert_intp_type=vintp_types("P","Z","I"),                    &
                   &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_ta_m  = is_variable_in_output(var_name="ta_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_ta_m
    IF (p_acc%l_ta_m) THEN
       cf_desc    = t_cf_var('air temperature', 'K', 'Temperature', datatype_flt)
       grib2_desc = grib2_var(0, 0, 0, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( list, 'ta_m', p_acc%temp,                                     &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                   & ldims=shape3d_c,                                              &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &             vert_intp_type=vintp_types("P","Z","I"),          &
                   &             vert_intp_method=VINTP_METHOD_LIN ),              &
                   & in_group=groups("prog_timemean","atmo_timemean"))
    END IF

    p_acc%l_ps_m  = is_variable_in_output(var_name="ps_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_ps_m
    IF (p_acc%l_ps_m) THEN
       cf_desc    = t_cf_var('surface_air_pressure', 'Pa', 'surface pressure (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( list, 'ps_m', p_acc%pres_sfc,                                 &
                   & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                   & ldims=shape2d_c,                                              &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_psl_m = is_variable_in_output(var_name="psl_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_psl_m
    IF (p_acc%l_psl_m) THEN
       cf_desc    = t_cf_var('mean sea level pressure', 'Pa',                      &
         &                   'mean sea level pressure (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 3, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( list, 'psl_m', p_acc%pres_msl,                                &
                   & GRID_UNSTRUCTURED_CELL, ZA_MEANSEA, cf_desc, grib2_desc,      &
                   & ldims=shape2d_c,                                              &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_pfull_m = is_variable_in_output(var_name="pfull_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_pfull_m
    IF (p_acc%l_pfull_m) THEN
       cf_desc    = t_cf_var('air_pressure', 'Pa', 'pressure at full level (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 3, 0, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( list, 'pfull_m', p_acc%pres,                                  &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                   & ldims=shape3d_c, lrestart=.FALSE. ,                           &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &             vert_intp_type=vintp_types("P","Z","I"),          &
                   &             vert_intp_method=VINTP_METHOD_PRES ),             &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_phalf_m = is_variable_in_output(var_name="phalf_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_phalf_m
    IF (p_acc%l_phalf_m) THEN
       cf_desc    = t_cf_var('air_pressure', 'Pa', 'pressure at half level (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( list, 'phalf_m', p_acc%pres_ifc,                              &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                   & ldims=shape3d_chalf, lrestart=.FALSE.,                        &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &             vert_intp_type=vintp_types("P","Z","I"),          &
                   &             vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),       &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF

    p_acc%l_wap_m = is_variable_in_output(var_name="wap_m")
    p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_wap_m
    IF (p_acc%l_wap_m) THEN
       cf_desc    = t_cf_var('omega', 'Pa/s', 'vertical velocity (time mean)', datatype_flt)
       grib2_desc = grib2_var(0, 2, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( list,"wap_m", p_acc%omega,                                    &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                            &
                   & cf_desc, grib2_desc,                                          &
                   & ldims=shape3d_c,                                              &
                   & vert_interp=create_vert_interp_metadata(                      &
                   &             vert_intp_type=vintp_types("P","Z","I"),          &
                   &             vert_intp_method=VINTP_METHOD_LIN,                &
                   &             l_loglin=.FALSE., l_extrapol=.FALSE.),            &
                   & in_group=groups("prog_timemean","atmo_timemean") )
    END IF
    ! }}}

    ! TRACERS {{{
    ! support qv,qc,qi because they are always there
    IF (ntracer > 0) THEN
       p_acc%l_tracer_m = is_variable_in_output(var_name="hus_m") .OR. &
                        & is_variable_in_output(var_name="clw_m") .OR. &
                        & is_variable_in_output(var_name="cli_m")
       p_acc%l_any_m = p_acc%l_any_m .OR. p_acc%l_tracer_m
       IF (p_acc%l_tracer_m) THEN
          cf_desc    = t_cf_var('tracer', 'kg kg-1', 'air tracer (time mean)', datatype_flt)
          grib2_desc = grib2_var(0,20,2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( list, 'tracer_m', p_acc%tracer,                         &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                      & ldims=shape4d_c ,                                       &
                      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

          ALLOCATE(p_acc%tracer_ptr(ntracer))
          DO jt=1,ntracer
             IF (jt == iqv ) CALL add_ref(                                          &
                  &  list, 'tracer_m', 'hus_m', p_acc%tracer_ptr(jt)%p,             &
                  &  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                  &  t_cf_var('specific_humidity', 'kg kg-1',                       &
                  &           'specific humidity (time mean)', datatype_flt),       &
                  &  grib2_var( 0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),         &
                  &  ref_idx=jt,                                                    &
                  &  ldims=shape3d_c,                                               &
                  &  tlev_source=1,                                                 &
                  &  tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                  &              ihadv_tracer=advconf%ihadv_tracer(iqv),            &
                  &              ivadv_tracer=advconf%ivadv_tracer(iqv)),           &
                  &  vert_interp=create_vert_interp_metadata(                       &
                  &              vert_intp_type=vintp_types("P","Z","I"),           &
                  &              vert_intp_method=VINTP_METHOD_QV,                  &
                  &              l_satlimit=.FALSE.,                                &
                  &              lower_limit=2.5e-7_wp, l_restore_pbldev=.FALSE. ), &
                  &  in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                  &                  "tracer_timemean","atmo_timemean"))

             IF ( jt == iqc )  CALL add_ref(                                        &
                  &  list, 'tracer_m', 'clw_m', p_acc%tracer_ptr(jt)%p,             &
                  &  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                  &  t_cf_var('specific_cloud_water_content', 'kg kg-1',            &
                  &           'specific cloud water content (time mean)',datatype_flt), &
                  &  grib2_var(0, 1, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL),         &
                  &  ref_idx=jt,                                                    &
                  &  ldims=shape3d_c,                                               &
                  &  tlev_source=1,                                                 &
                  &  tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                  &              ihadv_tracer=advconf%ihadv_tracer(iqc),            &
                  &              ivadv_tracer=advconf%ivadv_tracer(iqc)),           &
                  &  vert_interp=create_vert_interp_metadata(                       &
                  &              vert_intp_type=vintp_types("P","Z","I"),           &
                  &              vert_intp_method=VINTP_METHOD_LIN,                 &
                  &              l_loglin=.FALSE.,                                  &
                  &              l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                  &              lower_limit=0._wp  ),                              &
                  &  in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                  &                  "tracer_timemean","atmo_timemean"))

             IF ( jt == iqi ) CALL add_ref(                                         &
                  &  list, 'tracer_m', 'cli_m', p_acc%tracer_ptr(jt)%p,             &
                  &  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                  &  t_cf_var('specific_cloud_ice_content', 'kg kg-1',              &
                  &           'specific cloud ice content (time mean)', datatype_flt),  &
                  &  grib2_var(0, 1, 82, ibits, GRID_UNSTRUCTURED, GRID_CELL),         &
                  &  ref_idx=jt,                                                    &
                  &  ldims=shape3d_c,                                               &
                  &  tlev_source=1,                                                 &
                  &  tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                  &              ihadv_tracer=advconf%ihadv_tracer(iqi),            &
                  &              ivadv_tracer=advconf%ivadv_tracer(iqi)),           &
                  &  vert_interp=create_vert_interp_metadata(                       &
                  &              vert_intp_type=vintp_types("P","Z","I"),           &
                  &              vert_intp_method=VINTP_METHOD_LIN,                 &
                  &              l_loglin=.FALSE.,                                  &
                  &              l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                  &              lower_limit=0._wp  ),                              &
                  &  in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                  &                  "tracer_timemean","atmo_timemean"))
          END DO
       END IF
    END IF
    ! }}}



    p_acc%numberOfAccumulations = 0
  END SUBROUTINE construct_opt_acc


  SUBROUTINE update_opt_acc(acc, nh_prog, rho, nh_diag, subset, levels)
    TYPE(t_nh_acc),  INTENT(INOUT)   :: acc
    TYPE(t_nh_prog), INTENT(IN)      :: nh_prog     ! for jg=1
    REAL(wp), INTENT(IN)             :: rho(:,:,:)
    TYPE(t_nh_diag), INTENT(IN)      :: nh_diag     ! for jg=1
    TYPE(t_subset_range), INTENT(IN) :: subset
    INTEGER , INTENT(IN)             :: levels

    INTEGER :: jt

    !WRITE(message_text,'(a,i2)') '(pre ): numberOfAccumulations:',acc%numberOfAccumulations
    !CALL message('update_opt_nh_acc', message_text)
    IF (acc%l_ua_m)    CALL add_fields(acc%u       , nh_diag%u       , subset, levels=levels)
    IF (acc%l_va_m)    CALL add_fields(acc%v       , nh_diag%v       , subset, levels=levels)
    IF (acc%l_wa_m)    CALL add_fields(acc%w       , nh_prog%w       , subset, levels=levels+1)
    IF (acc%l_rho_m) THEN
      CALL add_fields(acc%rho     , rho             , subset, levels=levels)
      CALL dbg_print('RHO Update FROM',nh_prog%rho  ,'opt_diag',5, in_subset=subset)
      CALL dbg_print('RHO Update TO  ',acc%rho      ,'opt_diag',5, in_subset=subset)
    ENDIF
    IF (acc%l_ta_m)    CALL add_fields(acc%temp    , nh_diag%temp    , subset, levels=levels)
    IF (acc%l_ps_m)    CALL add_fields(acc%pres_sfc, nh_diag%pres_sfc, subset)
    IF (acc%l_psl_m)   CALL add_fields(acc%pres_msl, nh_diag%pres_msl, subset)
    IF (acc%l_pfull_m) CALL add_fields(acc%pres    , nh_diag%pres    , subset, levels=levels)
    IF (acc%l_phalf_m) CALL add_fields(acc%pres_ifc, nh_diag%pres_ifc, subset, levels=levels+1)
    IF (acc%l_wap_m)   CALL add_fields(acc%omega   , nh_diag%omega   , subset, levels=levels)

    IF (ntracer > 0) THEN
       IF (acc%l_tracer_m) THEN
          DO jt=1,ntracer
             CALL add_fields(acc%tracer(:,:,:,jt) ,nh_prog%tracer(:,:,:,jt),subset,levels=levels)
             CALL dbg_print('Tracer Update FROM',nh_prog%tracer(:,:,:,jt),'opt_diag',5, in_subset=subset)
             CALL dbg_print('Tracer Update TO  ',acc%tracer(:,:,:,jt),    'opt_diag',5, in_subset=subset)
          END DO
       END IF
    END IF

    acc%numberOfAccumulations = acc%numberOfAccumulations + 1
    !WRITE(message_text,'(a,i2)') '(post): numberOfAccumulations:',acc%numberOfAccumulations
    !CALL message('update_opt_nh_acc', message_text)

  END SUBROUTINE update_opt_acc


  SUBROUTINE reset_opt_acc(acc)
    TYPE(t_nh_acc) :: acc!(n_dom)

    INTEGER :: jt

    IF (acc%l_ua_m)    acc%u        = 0.0_wp
    IF (acc%l_va_m)    acc%v        = 0.0_wp
    IF (acc%l_wa_m)    acc%w        = 0.0_wp
    IF (acc%l_rho_m)   acc%rho      = 0.0_wp
    IF (acc%l_ta_m)    acc%temp     = 0.0_wp
    IF (acc%l_ps_m)    acc%pres_sfc = 0.0_wp
    IF (acc%l_psl_m)   acc%pres_msl = 0.0_wp
    IF (acc%l_pfull_m) acc%pres     = 0.0_wp
    IF (acc%l_phalf_m) acc%pres_ifc = 0.0_wp
    IF (acc%l_wap_m)   acc%omega    = 0.0_wp

    IF (ntracer > 0) THEN
       IF (acc%l_tracer_m) THEN
          DO jt=1,ntracer
             acc%tracer(:,:,:,jt) = 0.0_wp
          END DO
       END IF
    END IF

    acc%numberOfAccumulations = 0

  END SUBROUTINE reset_opt_acc


  SUBROUTINE calc_mean_opt_acc(acc)
    TYPE(t_nh_acc) :: acc!(n_dom)

    INTEGER  :: jt
    REAL(wp) :: xfactor

    xfactor = 1._wp/REAL(acc%numberOfAccumulations,wp)

    IF (acc%l_ua_m)    acc%u        = acc%u        *xfactor
    IF (acc%l_va_m)    acc%v        = acc%v        *xfactor
    IF (acc%l_wa_m)    acc%w        = acc%w        *xfactor
    IF (acc%l_rho_m)   acc%rho      = acc%rho      *xfactor
    IF (acc%l_ta_m)    acc%temp     = acc%temp     *xfactor
    IF (acc%l_ps_m)    acc%pres_sfc = acc%pres_sfc *xfactor
    IF (acc%l_psl_m)   acc%pres_msl = acc%pres_msl *xfactor
    IF (acc%l_pfull_m) acc%pres     = acc%pres     *xfactor
    IF (acc%l_phalf_m) acc%pres_ifc = acc%pres_ifc *xfactor
    IF (acc%l_wap_m)   acc%omega    = acc%omega    *xfactor

    IF (ntracer > 0) THEN
       IF (acc%l_tracer_m) THEN
          DO jt=1,ntracer
             acc%tracer(:,:,:,jt) = acc%tracer(:,:,:,jt) *xfactor
          END DO
       END IF
    END IF

  END SUBROUTINE calc_mean_opt_acc

  !-------------
  !
  !> Add optional diagnostic variable lists (might remain empty)
  !
  SUBROUTINE construct_opt_diag(p_patch, l_init_pz)
    TYPE(t_patch),        INTENT(IN)   :: p_patch(n_dom)
    LOGICAL,              INTENT(IN)   :: l_init_pz

    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":construct_opt_diag"
    INTEGER                   :: jg, ist
    CHARACTER(LEN=2)          :: dom_str
    CHARACTER(:), ALLOCATABLE :: model_type

    ! initialize data structure for optional diagnostics
    ALLOCATE(p_nh_opt_diag(n_dom), STAT=ist)
    IF (ist /= SUCCESS) &
      CALL finish (routine, 'Allocation of optional diagnostics failed')

    model_type = get_my_process_name()

    DO jg = 1, n_dom
      WRITE(dom_str, "(i2.2)") jg

      CALL vlr_add(p_nh_opt_diag(jg)%opt_diag_list, &
        & 'nh_state_opt_diag_of_domain_'//dom_str, &
        & patch_id=p_patch(jg)%id, vlevel_type=level_type_ml, lrestart=.FALSE., &
        & model_type=model_type)

      IF (.NOT. l_init_pz) CYCLE

      CALL vlr_add(p_nh_opt_diag(jg)%opt_diag_list_z, &
        & 'nh_state_opt_diag_z_of_domain_'//dom_str, &
        & patch_id=p_patch(jg)%id, vlevel_type=level_type_hl, lrestart=.FALSE., &
        & model_type=model_type)

      CALL vlr_add(p_nh_opt_diag(jg)%opt_diag_list_p, &
        & 'nh_state_opt_diag_p_of_domain_'//dom_str, &
        & patch_id=p_patch(jg)%id, vlevel_type=level_type_pl, lrestart=.FALSE., &
        & model_type=model_type)

      CALL vlr_add(p_nh_opt_diag(jg)%opt_diag_list_i, &
        & 'nh_state_opt_diag_i_of_domain_'//dom_str, &
        & patch_id=p_patch(jg)%id, vlevel_type=level_type_il, lrestart=.FALSE., &
        & model_type=model_type)

      CALL vlr_add(p_nh_opt_diag(jg)%opt_acc_list, &
        & 'nh_accumulation_for_ProgAndDiag_of_domain_'//dom_str, &
        & patch_id=p_patch(jg)%id, vlevel_type=level_type_ml,    &
        & lrestart=.FALSE.,loutput=.TRUE.,                       &
        & model_type=model_type)
    ENDDO ! jg

    ! provisional construction of memory for a hardwired set of variables on domain 1
    CALL construct_opt_acc( p_patch(1),                    &
                          & p_nh_opt_diag(1)%opt_acc_list, &
                          & p_nh_opt_diag(1)%acc           )

  END SUBROUTINE construct_opt_diag


  !-------------
  !
  !> Clear optional diagnostic variable lists
  !
  SUBROUTINE destruct_opt_diag()
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":destruct_opt_diag"
    INTEGER :: jg, ist

    DO jg = 1, n_dom
      CALL vlr_del(p_nh_opt_diag(jg)%opt_diag_list_z)
      CALL vlr_del(p_nh_opt_diag(jg)%opt_diag_list_p)
      CALL vlr_del(p_nh_opt_diag(jg)%opt_diag_list_i)
      CALL vlr_del(p_nh_opt_diag(jg)%opt_diag_list  )
      CALL vlr_del(p_nh_opt_diag(jg)%opt_acc_list   )
    ENDDO ! jg

    ! Delete optional diagnostics
    DEALLOCATE(p_nh_opt_diag, STAT=ist)
    IF (ist /= SUCCESS) &
      CALL finish(routine,'Deallocation for optional diagnostics failed.')

  END SUBROUTINE destruct_opt_diag


  !-------------
  !>
  ! Initialize a variable containing coefficient tables for vertical
  ! interpolation. There exist to different kinds of coefficients: For
  ! p- and for z-level-interpolation.
  SUBROUTINE vcoeff_lin_allocate(nblks, nlev, vcoeff_lin)
    INTEGER,                   INTENT(IN)    :: nblks
    INTEGER,                   INTENT(IN)    :: nlev
    TYPE(t_vcoeff_lin),        INTENT(INOUT) :: vcoeff_lin

    CHARACTER(*), PARAMETER :: routine = modname//":vcoeff_lin_allocate"
    INTEGER :: ierrstat

    ! real(wp)
    ALLOCATE(vcoeff_lin%wfac_lin(nproma,nlev,nblks), STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! integer
    ALLOCATE( vcoeff_lin%idx0_lin(nproma,nlev,nblks), vcoeff_lin%bot_idx_lin(nproma,nblks),   &
      &       STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ALLOCATE( vcoeff_lin%wfacpbl1(nproma,nblks), vcoeff_lin%wfacpbl2(nproma,nblks),           &
      &       STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ALLOCATE( vcoeff_lin%kpbl1(nproma,nblks), vcoeff_lin%kpbl2(nproma,nblks),                 &
      &       vcoeff_lin%zextrap(nproma,nblks), STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    !$ACC ENTER DATA ASYNC(1) &
    !$ACC   CREATE(vcoeff_lin%wfac_lin, vcoeff_lin%idx0_lin, vcoeff_lin%bot_idx_lin) &
    !$ACC   CREATE(vcoeff_lin%wfacpbl1, vcoeff_lin%wfacpbl2, vcoeff_lin%kpbl1) &
    !$ACC   CREATE(vcoeff_lin%kpbl2, vcoeff_lin%zextrap)

    ! Initialization
    !$OMP PARALLEL
    CALL init(vcoeff_lin%wfac_lin, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(vcoeff_lin%idx0_lin, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(vcoeff_lin%bot_idx_lin, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(vcoeff_lin%wfacpbl1, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(vcoeff_lin%wfacpbl2, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(vcoeff_lin%kpbl1, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(vcoeff_lin%kpbl2, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(vcoeff_lin%zextrap, lacc=.TRUE., opt_acc_async=.TRUE.)
    !$OMP END PARALLEL
  END SUBROUTINE vcoeff_lin_allocate


  !-------------
  !>
  ! Initialize a variable containing coefficient tables for cubic
  ! vertical interpolation.
  !
  SUBROUTINE vcoeff_cub_allocate(nblks, nlev, vcoeff_cub)
    INTEGER,                   INTENT(IN)    :: nblks
    INTEGER,                   INTENT(IN)    :: nlev
    TYPE(t_vcoeff_cub),        INTENT(INOUT) :: vcoeff_cub

    CHARACTER(*), PARAMETER :: routine = modname//":vcoeff_cub_allocate"
    INTEGER :: ierrstat

    ! real(wp)
    ALLOCATE( vcoeff_cub%coef1(nproma,nlev,nblks), vcoeff_cub%coef2(nproma,nlev,nblks),     &
      &       vcoeff_cub%coef3(nproma,nlev,nblks), STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! integer
    ALLOCATE( vcoeff_cub%idx0_cub(nproma,nlev,nblks), vcoeff_cub%bot_idx_cub(nproma,nblks), &
      &       STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    !$ACC ENTER DATA ASYNC(1) &
    !$ACC   CREATE(vcoeff_cub%coef1, vcoeff_cub%coef2, vcoeff_cub%coef3) &
    !$ACC   CREATE(vcoeff_cub%idx0_cub, vcoeff_cub%bot_idx_cub)

    ! Initialization
    !$OMP PARALLEL
    CALL init(vcoeff_cub%coef1, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(vcoeff_cub%coef2, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(vcoeff_cub%coef3, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(vcoeff_cub%idx0_cub, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(vcoeff_cub%bot_idx_cub, lacc=.TRUE., opt_acc_async=.TRUE.)
    !$OMP END PARALLEL
  END SUBROUTINE vcoeff_cub_allocate


  !-------------
  !>
  ! Initialize a variable containing coefficient tables for vertical
  ! interpolation. There exist to different kinds of coefficients: For
  ! p- and for z-level-interpolation.
  SUBROUTINE vcoeff_allocate(nblks_c, nblks_e, nlev, vcoeff, lacc)
    INTEGER,                           INTENT(IN)    :: nblks_c, nblks_e
    INTEGER,                           INTENT(IN)    :: nlev
    TYPE(t_vcoeff),                    INTENT(INOUT) :: vcoeff
    LOGICAL, OPTIONAL,                 INTENT(IN)    :: lacc ! If true, use openacc
    CHARACTER(*), PARAMETER :: routine = modname//":vcoeff_allocate"

    CALL assert_acc_device_only(routine, lacc)
    IF (.NOT. vcoeff%l_allocated) THEN

      CALL vcoeff_lin_allocate(nblks_c, nlev, vcoeff%lin_cell)
      CALL vcoeff_lin_allocate(nblks_c, nlev, vcoeff%lin_cell_nlevp1)
      CALL vcoeff_lin_allocate(nblks_e, nlev, vcoeff%lin_edge)

      ! CUBIC interpolation coefficients:
      CALL vcoeff_cub_allocate(nblks_c, nlev, vcoeff%cub_cell)
      CALL vcoeff_cub_allocate(nblks_e, nlev, vcoeff%cub_edge)

      ! MJ: At the moment, we don't have to attach the allocatable components
      !     vcoeff%*%{coef1, coef2, coef3, idx0_cub, bot_idx_cub} (!$ACC ATTACH)
      !$ACC ENTER DATA CREATE(vcoeff) ASYNC(1)

      vcoeff%l_allocated = .TRUE.
    END IF
  END SUBROUTINE vcoeff_allocate


  !-------------
  !>
  ! Clear a coefficient tables for linear vertical interpolation.
  SUBROUTINE vcoeff_lin_deallocate(vcoeff_lin)
    TYPE(t_vcoeff_lin), INTENT(INOUT) :: vcoeff_lin

    CHARACTER(*), PARAMETER :: routine = modname//":vcoeff_lin_deallocate"
    INTEGER :: ierrstat

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(vcoeff_lin%wfac_lin, vcoeff_lin%idx0_lin, vcoeff_lin%bot_idx_lin) &
    !$ACC   DELETE(vcoeff_lin%wfacpbl1, vcoeff_lin%wfacpbl2, vcoeff_lin%kpbl1) &
    !$ACC   DELETE(vcoeff_lin%kpbl2, vcoeff_lin%zextrap)

    ! real(wp)
    DEALLOCATE( vcoeff_lin%wfac_lin, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    ! integer
    DEALLOCATE( vcoeff_lin%idx0_lin, vcoeff_lin%bot_idx_lin, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    DEALLOCATE( vcoeff_lin%wfacpbl1, vcoeff_lin%wfacpbl2, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    DEALLOCATE( vcoeff_lin%kpbl1, vcoeff_lin%kpbl2, vcoeff_lin%zextrap, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE vcoeff_lin_deallocate


  !-------------
  !>
  ! Clear a coefficient tables for cubic vertical interpolation.
  SUBROUTINE vcoeff_cub_deallocate(vcoeff_cub)
    TYPE(t_vcoeff_cub), INTENT(INOUT) :: vcoeff_cub

    CHARACTER(*), PARAMETER :: routine = modname//":vcoeff_cub_deallocate"
    INTEGER :: ierrstat

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(vcoeff_cub%coef1, vcoeff_cub%coef2, vcoeff_cub%coef3) &
    !$ACC   DELETE(vcoeff_cub%idx0_cub, vcoeff_cub%bot_idx_cub)

    ! CUBIC interpolation coefficients:
    ! real(wp)
    DEALLOCATE( vcoeff_cub%coef1, vcoeff_cub%coef2, vcoeff_cub%coef3, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    ! integer
    DEALLOCATE( vcoeff_cub%idx0_cub, vcoeff_cub%bot_idx_cub, STAT=ierrstat )
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE vcoeff_cub_deallocate


  !-------------
  !>
  ! Clear a variable containing coefficient tables for vertical
  ! interpolation. There exist to different kinds of coefficients: For
  ! p-, z- and for i-level-interpolation.
  SUBROUTINE vcoeff_deallocate(vcoeff, lacc)
    TYPE(t_vcoeff), INTENT(INOUT) :: vcoeff
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    CHARACTER(*), PARAMETER :: routine = modname//":vcoeff_deallocate"

    CALL assert_acc_device_only(routine, lacc)

    ! deallocate coefficient tables:
    IF (vcoeff%l_allocated) THEN
      !$ACC WAIT ! subsequent code is synchronous.
      ! MJ: If the sub components of vcoeff were ATTACHed in vcoeff_allocate, they
      !     would have to be DETACHed here. (!$ACC DETACH)
      CALL vcoeff_lin_deallocate(vcoeff%lin_cell)
      CALL vcoeff_lin_deallocate(vcoeff%lin_cell_nlevp1)
      CALL vcoeff_lin_deallocate(vcoeff%lin_edge)

      call vcoeff_cub_deallocate(vcoeff%cub_cell)
      call vcoeff_cub_deallocate(vcoeff%cub_edge)

      vcoeff%l_allocated = .FALSE.
#if defined(_CRAYFTN) && _RELEASE_MAJOR <= 16
      !ACCWA: Cray compiler (16.0.1) is too eager on the optimization
      !$ACC WAIT
#endif
      !$ACC EXIT DATA DELETE(vcoeff)
    END IF

    vcoeff%l_initialized = .FALSE.
  END SUBROUTINE vcoeff_deallocate


  !---------------------------------------------------------------
  !> Adds a special metrics variable containing the area weights of
  !  the regular lon-lat grid.
  !
  !  @note This new variable is time-constant!
  !
  SUBROUTINE compute_lonlat_area_weights(lonlat_data, p_nh_state_lists)
    TYPE (t_lon_lat_list), TARGET, INTENT(IN) :: lonlat_data
    TYPE (t_nh_state_lists), INTENT(INOUT) :: p_nh_state_lists(:)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::compute_lonlat_area_weights"
    TYPE(t_cf_var)       :: cf_desc
    TYPE(t_grib2_var)    :: grib2_desc
    INTEGER              :: var_shape(3), nblks_lonlat, i, jg, ierrstat, &
      &                     i_lat, idx_glb, jc, jb, j, datatype_flt
    TYPE (t_lon_lat_grid), POINTER :: grid
    TYPE (t_lon_lat_intp), POINTER :: ptr_int_lonlat
!DR      REAL(wp),              POINTER :: area_weights(:), p_dummy(:,:,:)
!DR !!! Using the POINTER attribute for area_weights(:) mysteriously leads
!DR !!! to a bus error on NEC SX9 (tested with compiler revision 450). However,
!DR !!! using the ALLOCATABLE attribute, instead, works.
    REAL(wp), ALLOCATABLE :: area_weights(:)
    REAL(wp), POINTER :: p_dummy(:,:,:)
    TYPE(t_var), POINTER :: new_element

    ! define NetCDF output precision
    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
    ! Add area weights
    DO i=1, lonlat_data%ngrids
      DO jg=1,n_dom
        IF (lonlat_data%list(i)%l_dom(jg)) THEN
          grid           => lonlat_data%list(i)%grid
          ptr_int_lonlat => lonlat_data%list(i)%intp(jg)
          nblks_lonlat   =  ptr_int_lonlat%nblks_lonlat(nproma)
          var_shape = (/ nproma, 1, nblks_lonlat /)
          cf_desc    = t_cf_var('aw', '1', 'area weights for regular lat-lon grid', datatype_flt)
          grib2_desc = grib2_var(0, 191, 193, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)

          ALLOCATE(area_weights(grid%lat_dim), STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

          CALL add_var( p_nh_state_lists(jg)%diag_list,                       &
            &           LONLAT_PREFIX//"aw", p_dummy,                         &
            &           GRID_REGULAR_LONLAT, ZA_SURFACE, cf_desc, grib2_desc, &
            &           ldims=var_shape, lrestart=.FALSE.,                    &
            &           loutput=.TRUE., new_element=new_element,              &
            &           hor_interp=create_hor_interp_metadata(                &
            &             hor_intp_type=HINTP_TYPE_NONE, lonlat_id=i ),       &
            &           isteptype=TSTEP_CONSTANT )
          ! compute area weights:
!CDIR NOIEXPAND
          CALL latlon_compute_area_weights(grid, earth_radius, area_weights)
          ! for each local lon-lat point on this PE:
          DO j=1, ptr_int_lonlat%nthis_local_pts
            ! determine block, index
            jb = (j-1)/nproma + 1
            jc = j - (jb-1)*nproma
            ! determine latitude index:
            idx_glb = ptr_int_lonlat%global_idx(j)
            i_lat   = (idx_glb-1)/grid%lon_dim + 1
            ! set area weight:
            p_dummy(jc,1,jb) = area_weights(i_lat)
          END DO
          DEALLOCATE(area_weights, STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
        END IF
      END DO
    END DO
  END SUBROUTINE compute_lonlat_area_weights

END MODULE mo_opt_diagnostics
