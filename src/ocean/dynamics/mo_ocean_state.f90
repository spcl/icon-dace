!  Contains the data structures for the hydrostatic ocean model.
!
!  Contains the data structures to store the hydrostatic & boussinesq ocean model state.
!  Implementation is based on ICON-Shallow-Water model
!  to store the shallow water model state and other auxiliary variables.
!  Constructors and destructors for these data structures are also defined here.
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

MODULE mo_ocean_state
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



  !-------------------------------------------------------------------------
  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants,      ONLY: success, TLEV_NNEW
  USE mo_ocean_nml,           ONLY: n_zlev, dzlev_m, no_tracer, use_tracer_x_height, cfl_write,&
    &                               Cartesian_Mixing , &
    &                               k_tracer_dianeutral_parameter,                                &
    &                               k_tracer_isoneutral_parameter, k_tracer_GM_kappa_parameter,   &
    &                               GMRedi_configuration,GMRedi_combined,                         &
    &                               GM_only,Redi_only, type_3dimrelax_salt, type_3dimrelax_temp,  &
    &                               GMREDI_COMBINED_DIAGNOSTIC,GM_INDIVIDUAL_DIAGNOSTIC,          &
    &                               REDI_INDIVIDUAL_DIAGNOSTIC, eddydiag,                         &
    &                               diagnose_for_tendencies, diagnose_for_heat_content, lhamocc,  &
    &                               use_tides_SAL, vert_cor_type
  USE mo_run_config,          ONLY: test_mode
  USE mo_ocean_types,         ONLY: t_hydro_ocean_base ,t_hydro_ocean_state ,t_hydro_ocean_prog ,t_hydro_ocean_diag, &
    &                               t_hydro_ocean_aux , t_oce_config
  USE mo_ocean_nudging_types, ONLY: t_ocean_nudge
  USE mo_ocean_nudging,       ONLY: ocean_nudge
  USE mo_mpi,                 ONLY: get_my_global_mpi_id, global_mpi_barrier,my_process_is_mpi_test
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants,      ONLY: land, land_boundary, boundary, sea_boundary, sea,     &
    &                               success, MIN_DOLIC,                                   &
    &                               full_coriolis, beta_plane_coriolis,                   &
    &                               f_plane_coriolis, zero_coriolis, halo_levels_ceiling
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge, grid_unstructured_cell,         &
    &                               grid_unstructured_edge, grid_unstructured_vert,       &
    &                               grid_vertex
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_model_domain,        ONLY: t_patch,t_patch_3d, t_grid_cells, t_grid_edges
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, grid_sphere_radius, grid_angular_velocity, &
    & use_dummy_cell_closure
  USE mo_dynamics_config,     ONLY: nnew, nold, nnow
  USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_var_list_register,   ONLY: vlr_add, vlr_del
  USE mo_var_list,            ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var_metadata,        ONLY: get_timelevel_string
  USE mo_var_groups,          ONLY: groups, MAX_GROUPS
  USE mo_cf_convention
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_grib2,               ONLY: grib2_var, t_grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_FLT32 => CDI_DATATYPE_FLT32, &
    &                               DATATYPE_FLT64 => CDI_DATATYPE_FLT64, &
    &                               DATATYPE_INT8 => CDI_DATATYPE_INT8, &
    &                               DATATYPE_PACK16 => CDI_DATATYPE_PACK16, &
    &                               tstep_constant, GRID_LONLAT, GRID_UNSTRUCTURED
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge, grid_unstructured_cell, grid_unstructured_edge, &
    &                               grid_unstructured_vert, grid_vertex, GRID_ZONAL
  USE mo_zaxis_type,          ONLY: za_depth_below_sea, za_depth_below_sea_half, za_surface
  !  USE mo_ocean_config,        ONLY: ignore_land_points
  USE mo_io_config,           ONLY: lnetcdf_flt64_output

  USE mo_ocean_tracer_transport_types
  USE mo_var_list_gpu,        ONLY: gpu_update_var_list
  USE mo_ocean_types,              ONLY: t_operator_coeff

  IMPLICIT NONE
  PRIVATE


  !public interface
  !
  ! subroutines
  PUBLIC :: construct_hydro_ocean_base
  PUBLIC :: destruct_hydro_ocean_base
  PUBLIC :: construct_hydro_ocean_state
  PUBLIC :: destruct_hydro_ocean_state
  PUBLIC :: construct_patch_3d, destruct_patch_3d
  PUBLIC :: construct_ocean_var_lists
  PUBLIC :: transfer_ocean_state
  !

  PUBLIC :: ocean_restart_list
  PUBLIC :: ocean_default_list
  PUBLIC :: v_base
  PUBLIC :: oce_config

  !constructors
  PRIVATE :: construct_hydro_ocean_diag
  PRIVATE :: construct_hydro_ocean_prog
  PRIVATE :: construct_hydro_ocean_aux
  PUBLIC  :: construct_ocean_nudge
  !destructors
  PRIVATE :: destruct_hydro_ocean_diag
  PRIVATE :: destruct_hydro_ocean_aux
  PRIVATE :: destruct_ocean_nudge

  !----------------------------------------------------------------------------
  PUBLIC  :: ocean_state
  !----------------------------------------------------------------------------

  ! variables
  TYPE(t_var_list_ptr)                              :: ocean_restart_list
  TYPE(t_var_list_ptr)                              :: ocean_default_list
  TYPE(t_hydro_ocean_base) ,TARGET :: v_base
  TYPE(t_oce_config)                            :: oce_config
  INTEGER, PARAMETER :: max_oce_tracer = 50
  !-------------------------------------------------------------------------

  TYPE(t_hydro_ocean_state), ALLOCATABLE, TARGET  :: ocean_state(:)

CONTAINS

  !-------------------------------------------------------------------------
  !
  !
!<Optimize:inUse>
  SUBROUTINE construct_ocean_var_lists(patch_2d)
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2d
    CHARACTER(:), ALLOCATABLE :: model_name

    model_name = TRIM(get_my_process_name())

    ! IMO the number of variable lists should be as small as possible
    !
    ! Restart list: everything belonging to that list will be written to the
    ! restart file and is ready for output
    CALL vlr_add(ocean_restart_list, 'ocean_restart_list', &
      & patch_id=patch_2d%id, lrestart=.TRUE., loutput=.TRUE.,  &
      & model_type=model_name)

    ! default list: elements can be written to disk, but not to the restart file
    CALL vlr_add(ocean_default_list, 'ocean_default_list', &
      & patch_id=patch_2d%id, lrestart=.FALSE., loutput=.TRUE., &
      & model_type=model_name)
  END SUBROUTINE construct_ocean_var_lists
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Constructor for hydrostatic ocean state + diagnostic and auxiliary  states.
  !!
  !! Constructor for hydrostatic ocean state
  !! It calls  constructors to single time level
  !! auxiliary and diagnostic states. Then it constructs state array,
  !! whose components (representing multiple time levels).
  !! Initialization of all components with zero.
  !!
  !
  !
!<Optimize:inUse>
  SUBROUTINE construct_hydro_ocean_state( patch_3d, ocean_state)

    TYPE(t_patch_3D), TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state(n_dom)

    ! local variables
    TYPE(t_patch), POINTER :: patch_2d

    INTEGER :: i_status, timelevel, prlength ! local prognostic array length
    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_state:construct_hydro_ocean_state'

    patch_2d => patch_3d%p_patch_2d(1)
    CALL message(routine, 'start to construct hydro_ocean state' )
    ocean_state(1)%patch_3d => patch_3d
    ! Using Adams-Bashforth semi-implicit timestepping with 3 prognostic time levels:
    prlength = 3

    CALL setup_tracer_info(oce_config)

    !$ACC ENTER DATA COPYIN(ocean_state)

    !create state array for each domain
    ALLOCATE(ocean_state(1)%p_prog(1:prlength), stat=i_status)
    IF (i_status/=success) THEN
      CALL finish(routine, 'allocation of progn. state array failed')
    END IF

    !$ACC ENTER DATA COPYIN(ocean_state(1)%p_prog)

    DO timelevel = 1,prlength
        ! timelevel nnow is not used by the ocean - therefore we dont allocate
        ! variables for this timelevel
        IF (timelevel .ne. nnow(1)) THEN
          CALL construct_hydro_ocean_prog(patch_3d, &
            &                             ocean_state(1)%p_prog(timelevel), &
            &                             TRIM(get_timelevel_string(timelevel)))

      END IF
    END DO

    !$ACC ENTER DATA COPYIN(ocean_state(1)%p_diag, ocean_state(1)%p_aux, ocean_state(1)%transport_state)
    CALL construct_hydro_ocean_diag(patch_2d, ocean_state(1)%p_diag)
    CALL construct_hydro_ocean_aux(patch_2d,  ocean_state(1)%p_aux, ocean_state(1)%transport_state)

    CALL message(routine,'construction of hydrostatic ocean state finished')

  END SUBROUTINE construct_hydro_ocean_state
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  !!               Destructor for hydrostatic ocean state.
  !
  !
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_hydro_ocean_state(ocean_state)
    TYPE(t_hydro_ocean_state), TARGET,INTENT(inout)   :: ocean_state(n_dom)

    ! local variables

    INTEGER :: jg, prlength, ist

    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_state:destruct_hydro_ocean_state'

    !-------------------------------------------------------------------------
    CALL message(routine, 'start to destruct hydro ocean state ')

    prlength = SIZE(ocean_state(1)%p_prog)

    IF (prlength==0) THEN
      CALL finish(routine, 'prog array has length zero')
    END IF

    CALL vlr_del(ocean_restart_list)
    CALL vlr_del(ocean_default_list)

    DO jg = 1, n_dom
      CALL destruct_hydro_ocean_diag(ocean_state(jg)%p_diag)
      CALL destruct_hydro_ocean_aux (ocean_state(jg)%p_aux)

      ! destruct state array
      ist = 1
      DEALLOCATE(ocean_state(jg)%p_prog, stat=ist)
      IF (ist/=success) THEN
        CALL finish(routine, 'deallocation of state array failed')
      END IF
    END DO

    CALL message(routine, 'destruction of hydrostatic ocean state finished')

  END SUBROUTINE destruct_hydro_ocean_state

  !-------------------------------------------------------------------------
  !>
  !! Allocation of basic 3-dimensional structure components of hydrostatic ocean state.
  !! Initialization of components with zero.
  !
  !
  !!

!<Optimize:inUse>
  SUBROUTINE construct_hydro_ocean_base(patch_2d, v_base)

    TYPE(t_patch), TARGET, INTENT(in)          :: patch_2d
    TYPE(t_hydro_ocean_base), INTENT(inout)    :: v_base

    ! local variables

    INTEGER :: ist
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v, n_zlvp
    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_state:construct_hydro_ocean_base'

    !-------------------------------------------------------------------------

    !CALL message(routine, 'start to construct basic hydro ocean state')

    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v
    n_zlvp = n_zlev + 1

    ! allocate and set vertical level thickness from the namelist
    ALLOCATE(v_base%del_zlev_m(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating del_zlev_m failed')
    ENDIF
    ALLOCATE(v_base%zlev_m(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating zlev_m failed')
    ENDIF
    ALLOCATE(v_base%zlev_i(n_zlvp),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating zlev_i failed')
    ENDIF
    ALLOCATE(v_base%del_zlev_i(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating del_zlev_i failed')
    ENDIF

    !
    !! 3-dim land-sea-mask at cells, edges and vertices
    !
    ! cells
    ALLOCATE(v_base%lsm_c(nproma,n_zlev,alloc_cell_blocks),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating lsm_c failed')
    ENDIF
    ! edges
    ALLOCATE(v_base%lsm_e(nproma,n_zlev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating lsm_e failed')
    ENDIF
    ! deepest ocean layer in column
    ALLOCATE(v_base%dolic_c(nproma,alloc_cell_blocks),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating dolic_c failed')
    ENDIF
    ALLOCATE(v_base%dolic_e(nproma,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating dolic_e failed')
    ENDIF
    ! 2-dim basins and areas
    ALLOCATE(v_base%basin_c(nproma,alloc_cell_blocks),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating basin_c failed')
    ENDIF
    ALLOCATE(v_base%regio_c(nproma,alloc_cell_blocks),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating regio_c failed')
    ENDIF
    ! 3-dim real land-sea-mask
    ! cells
    ALLOCATE(v_base%wet_c(nproma,n_zlev,alloc_cell_blocks),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating wet_c failed')
    ENDIF
    ! edges
    ALLOCATE(v_base%wet_e(nproma,n_zlev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating wet_e failed')
    ENDIF

    v_base%del_zlev_m = 0._wp
    v_base%del_zlev_i = 0._wp
    v_base%zlev_m     = 0._wp
    v_base%zlev_i     = 0._wp

    v_base%lsm_c = 0
    v_base%lsm_e = 0
    v_base%dolic_c = 0
    v_base%dolic_e = 0
    v_base%basin_c = 0
    v_base%regio_c = 0

    v_base%wet_c = 0.0_wp
    v_base%wet_e = 0.0_wp

  END SUBROUTINE construct_hydro_ocean_base

  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of diagnostic hydrostatic ocean state.
  !
  !!
  SUBROUTINE destruct_hydro_ocean_base(v_base)

    TYPE(t_hydro_ocean_base), INTENT(inout) :: v_base

    ! local variables

    INTEGER :: ist

    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_state:destruct_hydro_ocean_base'

    CALL message(routine,' start to destruct hydrostatic ocean basic state')

    DEALLOCATE(v_base%zlev_m,stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'deallocating zlev_m failed')
    ENDIF

  END SUBROUTINE destruct_hydro_ocean_base

  !-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean prognostic state.
  !!               Initialization of components with zero.
  !
  !
  !!
!<Optimize:inUse>
  SUBROUTINE construct_hydro_ocean_prog(patch_3d, ocean_state_prog, var_suffix)

    TYPE(t_patch_3D), TARGET, INTENT(in)      :: patch_3d
    TYPE(t_hydro_ocean_prog), INTENT(inout)   :: ocean_state_prog
    CHARACTER(LEN=*)                          :: var_suffix

    INTEGER :: alloc_cell_blocks, nblks_e !, nblks_v
    INTEGER :: jtrc
    INTEGER :: datatype_flt
    TYPE(t_ocean_tracer), POINTER :: tracer
    TYPE(t_patch), POINTER         :: patch_2d
    LOGICAL :: oce_tr_groups(MAX_GROUPS)

    oce_tr_groups = groups("oce_default", "oce_essentials","oce_prog")
    patch_2d => patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e

    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)

      IF (vert_cor_type == 0) THEN
        ! height
        CALL add_var(ocean_restart_list, 'zos'//TRIM(var_suffix), ocean_state_prog%h , &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
          & t_cf_var('zos'//TRIM(var_suffix), 'm', 'surface elevation at cell center', DATATYPE_FLT,'zos'),&
          & grib2_var(10, 3, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
          & ldims=(/nproma,alloc_cell_blocks/), tlev_source=TLEV_NNEW,&
          & in_group=groups("oce_default", "oce_essentials","oce_prog"), lopenacc = .TRUE.)
        
        ! zstar stretching
        CALL add_var(ocean_default_list, 'stretch_c'//TRIM(var_suffix), ocean_state_prog%stretch_c , &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
          & t_cf_var('stretch_c'//TRIM(var_suffix), 'm', 'zstar surface stretch at cell center', &
          & DATATYPE_FLT,'stretch_c'),&
          & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
          & ldims=(/nproma,alloc_cell_blocks/), tlev_source=TLEV_NNEW, lopenacc = .TRUE., initval=0.0_wp)
        
      END IF

      IF (vert_cor_type == 1) THEN
        ! zstar height
        CALL add_var(ocean_restart_list, 'zos'//TRIM(var_suffix), ocean_state_prog%eta_c , &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
          & t_cf_var('zos'//TRIM(var_suffix), 'm', 'zstar sfc elevation at cell center', DATATYPE_FLT,'zos'),&
          & grib2_var(10, 3, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
          & ldims=(/nproma,alloc_cell_blocks/), tlev_source=TLEV_NNEW,&
          & in_group=groups("oce_default", "oce_essentials","oce_prog"))

        ! zstar stretching
        CALL add_var(ocean_restart_list, 'stretch_c'//TRIM(var_suffix), ocean_state_prog%stretch_c , &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
          & t_cf_var('stretch_c'//TRIM(var_suffix), '1', 'zstar surface stretch at cell center', &
          & DATATYPE_FLT,'stretch_c'),&
          & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
          & ldims=(/nproma,alloc_cell_blocks/), tlev_source=TLEV_NNEW,&
          & in_group=groups("oce_default", "oce_essentials","oce_prog"), lrestart_cont=.TRUE.)
        
        ! height
        CALL add_var(ocean_default_list, 'z_ht'//TRIM(var_suffix), ocean_state_prog%h , &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
          & t_cf_var('z_ht'//TRIM(var_suffix), 'm', 'dummy elevation at cell center', DATATYPE_FLT,'z_ht'),&
          & grib2_var(255, 255, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
          & ldims=(/nproma,alloc_cell_blocks/), tlev_source=TLEV_NNEW,&
          & in_group=groups("oce_default", "oce_essentials","oce_prog"))

      ENDIF

      !! normal velocity component
      CALL add_var(ocean_restart_list,'normal_velocity'//var_suffix,ocean_state_prog%vn,grid_unstructured_edge, &
        & za_depth_below_sea, &
        & t_cf_var('vn'//var_suffix, 'm s-1', 'normal velocity on edge', DATATYPE_FLT),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/), tlev_source=TLEV_NNEW, lopenacc = .TRUE.)
      

      !! Tracers
      IF ( no_tracer > 0 ) THEN
        CALL add_var(ocean_restart_list, 'tracers'//var_suffix, ocean_state_prog%tracer , &
          & grid_unstructured_cell, za_depth_below_sea, &
          & t_cf_var('tracers'//var_suffix, '', '1:temperature 2:salinity', &
          & DATATYPE_FLT),&
          & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
          & ldims=(/nproma,n_zlev,alloc_cell_blocks,no_tracer/), &
          & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc = .TRUE.)
        

        ! Reference to individual tracer, for I/O
        ALLOCATE(ocean_state_prog%tracer_ptr(no_tracer))

        DO jtrc = 1,no_tracer
        CALL add_ref( ocean_restart_list, 'tracers'//var_suffix,   &
            & TRIM(oce_config%tracer_shortnames(jtrc))//var_suffix,        &
            & ocean_state_prog%tracer_ptr(jtrc)%p,                         &
            & grid_unstructured_cell, za_depth_below_sea,                  &
            & t_cf_var(TRIM(oce_config%tracer_stdnames(jtrc)), &
            & TRIM(oce_config%tracer_units(jtrc)), &
            & TRIM(oce_config%tracer_longnames(jtrc)), &
            & DATATYPE_FLT, &
            & TRIM(oce_config%tracer_shortnames(jtrc))), &
            ! GRIB codes for tracer are the three lower bytes of oce_config%tracer_codes:
            ! discipline, parameterCategory, parameterNumber
            & grib2_var(ISHFT(oce_config%tracer_codes(jtrc), -16), &
              & IAND(ISHFT(oce_config%tracer_codes(jtrc), -8), INT(255)), &  ! 255 = hex FF to mask out right-most byte
              & IAND(oce_config%tracer_codes(jtrc), INT(255)), DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell), &
            & ref_idx=jtrc, &
            & ldims=(/nproma,n_zlev,alloc_cell_blocks/), tlev_source=TLEV_NNEW, &
            & in_group=oce_tr_groups)
        END DO

        !--------------------------------------------------------------------------
        ! use of the ocean_tracers structure
        ALLOCATE(ocean_state_prog%tracer_collection%tracer(no_tracer))

        !$ACC ENTER DATA COPYIN(ocean_state_prog%tracer_collection, ocean_state_prog%tracer_collection%tracer)

        ocean_state_prog%tracer_collection%typeOfTracers = "ocean"
        ocean_state_prog%tracer_collection%no_of_tracers = no_tracer
        ocean_state_prog%tracer_collection%patch_3d => patch_3d

        DO jtrc = 1,no_tracer
          tracer => ocean_state_prog%tracer_collection%tracer(jtrc)
          ! point the concentration to the 4D tracer
          ! this is a tmeporary solution until the whole code is cleaned
          tracer%concentration => ocean_state_prog%tracer(:,:,:,jtrc)
          NULLIFY(tracer%top_bc)
          NULLIFY(tracer%bottom_bc)
          IF (jtrc <= no_tracer) THEN
            tracer%is_advected = .true.
          ELSE
            tracer%is_advected = .false.
          ENDIF

          ! allocate a
  !         IF (use_tracer_x_height) THEN
  !           !
  !           CALL add_var(ocean_restart_list, 'ocean_tracers'//var_suffix, &
  !             & ocean_state_prog%ocean_tracers(jtrc)%concentration_x_height , &
  !             & grid_unstructured_cell, za_depth_below_sea, &
  !             & t_cf_var(TRIM(oce_tracer_names(jtrc))//"_x_height",  &
  !             & oce_tracer_units(jtrc),                            &
  !             & oce_tracer_longnames(jtrc), DATATYPE_FLT64),       &
  !             & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
  !             & lopenacc=.TRUE., ldims=(/nproma,n_zlev,alloc_cell_blocks/))
  !
  !         ENDIF ! use_tracer_x_height
            !$ACC ENTER DATA COPYIN(tracer%top_bc, tracer%bottom_bc, tracer%is_advected)
        ENDDO
        !---------------------------------------------------------------------------
      ELSE

        ocean_state_prog%tracer_collection%no_of_tracers = 0

        !$ACC ENTER DATA COPYIN(ocean_state_prog%tracer_collection, ocean_state_prog%tracer_collection%no_of_tracers)

      ENDIF ! no_tracer > 0

  END SUBROUTINE construct_hydro_ocean_prog
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean diagnostic state.
  !!               Initialization of components with zero.
  !
  !!

!<Optimize:inUse>
  SUBROUTINE construct_hydro_ocean_diag(patch_2d,ocean_state_diag)

    TYPE(t_patch), TARGET, INTENT(in)          :: patch_2d
    TYPE(t_hydro_ocean_diag), INTENT(inout)    :: ocean_state_diag

    ! local variables

    INTEGER :: ist
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v
    !INTEGER :: jb, jc, jk, je
    !INTEGER :: i_startidx_c, i_endidx_c, i_startidx_e, i_endidx_e
    !INTEGER, PARAMETER :: max_oce_tracer = 2
    INTEGER :: oce_tracer_codes(max_oce_tracer)
    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_state:construct_hydro_ocean_diag'
    INTEGER :: datatype_flt, jc, blockNo, start_cell_index, end_cell_index
    REAL(wp), PARAMETER :: equator = 0.00001_wp
    TYPE(t_subset_range), POINTER :: owned_cells
    LOGICAL, DIMENSION(MAX_GROUPS) :: groups_oce_monitor, groups_oce_flows, &
         groups_oce_eddy, groups_oce_diag, groups_oce_default, groups_oce_moc, &
         groups_oce_dde
    TYPE(t_grib2_var) :: dflt_g2_decl_lonlat, dflt_g2_decl_cell, dflt_g2_decl_edge

    dflt_g2_decl_cell = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell)
    dflt_g2_decl_edge = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge)
    dflt_g2_decl_lonlat = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat)

    groups_oce_monitor = groups("ocean_monitor")
    groups_oce_flows = groups("ocean_flows")
    groups_oce_diag = groups("oce_diag")
    groups_oce_eddy = groups("oce_eddy")
    groups_oce_moc = groups("ocean_moc")
    groups_oce_default = groups("oce_default")
    groups_oce_dde = groups("oce_diag","oce_default","oce_essentials")

    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)

    !-------------------------------------------------------------------------
    CALL message(routine, 'construct diagnostic hydro ocean state...')

    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v

    ! add monitoring {{{

   CALL add_var(ocean_default_list, 'total_salt', ocean_state_diag%monitor%total_salt , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('total_salt', 'kg', 'total_salt', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

   CALL add_var(ocean_default_list, 'total_saltinseaice', ocean_state_diag%monitor%total_saltinseaice , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('total_saltinseaice', 'kg', 'total_saltinseaice', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'total_saltinliquidwater', ocean_state_diag%monitor%total_saltinliquidwater , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('total_saltinliquidwater', 'kg', 'total_saltinliquidwater', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))
    CALL add_var(ocean_default_list, 'amoc26n', ocean_state_diag%monitor%amoc26n , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('amoc26n', 'kg s-1', 'amoc26n', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'kin_energy_global', ocean_state_diag%monitor%kin_energy , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('kin_energy', 'J', 'kin_energy', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'pot_energy_global', ocean_state_diag%monitor%pot_energy , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('pot_energy', 'J', 'pot_energy', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'total_energy_global', ocean_state_diag%monitor%total_energy , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('total_energy', 'J', 'total_energy', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'ssh_global', ocean_state_diag%monitor%ssh_global , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ssh_global', 'm', 'ssh_global', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'sst_global', ocean_state_diag%monitor%sst_global , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('sea_surface_temperature', 'C', 'global mean sea surface temperature', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'sss_global', ocean_state_diag%monitor%sss_global , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('sea_surface_salinity', 'kg m-3', 'global mean sea surface salinity', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'potential_enstrophy_global', ocean_state_diag%monitor%potential_enstrophy , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('potential_enstrophy', 'm2 s-2', 'potential_enstrophy', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'HeatFlux_Total_global', ocean_state_diag%monitor%HeatFlux_Total , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_Total', 'W m-2', 'HeatFlux_Total', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_Precipitation_global', ocean_state_diag%monitor%FrshFlux_Precipitation , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_Precipitation', 'm s-1', 'FrshFlux_Precipitation', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_SnowFall_global', ocean_state_diag%monitor%FrshFlux_SnowFall , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_SnowFall', 'm s-1', 'FrshFlux_SnowFall_global', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_Evaporation_global', ocean_state_diag%monitor%FrshFlux_Evaporation , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_Evaporation', 'm s-1', 'FrshFlux_Evaporation_global', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_Runoff_global', ocean_state_diag%monitor%FrshFlux_Runoff , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_Runoff', 'm s-1', 'FrshFlux_Runoff', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_VolumeIce_global', ocean_state_diag%monitor%FrshFlux_VolumeIce , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_VolumeIce', 'm s-1', 'FrshFlux_VolumeIce', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_TotalOcean_global', ocean_state_diag%monitor%FrshFlux_TotalOcean , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_TotalOcean', 'm s-1', 'FrshFlux_TotalOcean', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_TotalIce_global', ocean_state_diag%monitor%FrshFlux_TotalIce , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_TotalIce', 'm s-1', 'FrshFlux_TotalIce', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_VolumeTotal_global', ocean_state_diag%monitor%FrshFlux_VolumeTotal , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_VolumeTotal', 'm s-1', 'FrshFlux_VolumeTotal', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'totalsnowfall_global', ocean_state_diag%monitor%totalsnowfall , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('totalsnowfall', 'm s-1', 'totalsnowfall', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_volume_nh', ocean_state_diag%monitor%ice_volume_nh , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_volume_nh', 'km3', 'ice_volume_nh', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_volume_sh', ocean_state_diag%monitor%ice_volume_sh , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_volume_sh', 'km3', 'ice_volume_sh', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_extent_nh', ocean_state_diag%monitor%ice_extent_nh , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_extent_nh', 'km2', 'ice_extent_nh', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'ice_extent_sh', ocean_state_diag%monitor%ice_extent_sh , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_extent_sh', 'km2', 'ice_extent_sh', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))
    ! }}}

! { not yet in monitoring
    CALL add_var(ocean_default_list, 'FrshFlux_TotalSalt_global', ocean_state_diag%monitor%FrshFlux_TotalSalt , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_TotalSalt', 'm s-1', 'FrshFlux_TotalSalt', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & ldims=(/1/))
!     & in_group=groups_oce_monitor,ldims=(/1/))
! }

    ! ocean through-flows {{{
    CALL add_var(ocean_default_list, 'ice_framStrait', ocean_state_diag%monitor%ice_framStrait , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('ice_framStrait', 'm3 s-1', 'ice_framStrait', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'gibraltar', ocean_state_diag%monitor%gibraltar , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('gibraltar', 'Sv', 'gibraltar', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'denmark_strait', ocean_state_diag%monitor%denmark_strait , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('denmark_strait', 'Sv', 'denmark_strait', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'drake_passage', ocean_state_diag%monitor%drake_passage , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('drake_passage', 'Sv', 'drake_passage', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'indonesian_throughflow', ocean_state_diag%monitor%indonesian_throughflow , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('indonesian_throughflow', 'Sv', 'indonesian_throughflow', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'scotland_iceland', ocean_state_diag%monitor%scotland_iceland , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('scotland_iceland', 'Sv', 'scotland_iceland', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'mozambique', ocean_state_diag%monitor%mozambique , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('mozambique', 'Sv', 'mozambique', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'framStrait', ocean_state_diag%monitor%framStrait , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('framStrait', 'Sv', 'framStrait', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'beringStrait', ocean_state_diag%monitor%beringStrait , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('beringStrait', 'Sv', 'beringStrait', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'barentsOpening', ocean_state_diag%monitor%barentsOpening , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('barentsOpening', 'Sv', 'barentsOpening', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'agulhas', ocean_state_diag%monitor%agulhas , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('agulhas', 'Sv', 'agulhas', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'agulhas_long', ocean_state_diag%monitor%agulhas_long , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('agulhas_long', 'Sv', 'agulhas_long', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'agulhas_longer', ocean_state_diag%monitor%agulhas_longer , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('agulhas_longer', 'Sv', 'agulhas_longer', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))

    CALL add_var(ocean_default_list, 'florida_strait', ocean_state_diag%monitor%florida_strait , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('florida_strait', 'Sv', 'florida_strait', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_flows,ldims=(/1/))
    ! }}}

    ! currently NOT maintained monitoring quantities {{{

    CALL add_var(ocean_default_list, 'HeatFlux_Relax_global', ocean_state_diag%monitor%HeatFlux_Relax , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_Relax', 'W m-2', 'HeatFlux_Relax', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & ldims=(/1/))

    CALL add_var(ocean_default_list, 'FrshFlux_Relax_global', ocean_state_diag%monitor%FrshFlux_Relax , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('FrshFlux_Relax', 'm s-1', 'FrshFlux_Relax', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & ldims=(/1/))

    CALL add_var(ocean_default_list, 'TempFlux_Relax_global', ocean_state_diag%monitor%TempFlux_Relax , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('TempFlux_Relax', 'K s-1', 'TempFlux_Relax', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & ldims=(/1/))

    CALL add_var(ocean_default_list, 'SaltFlux_Relax_global', ocean_state_diag%monitor%SaltFlux_Relax , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('SaltFlux_Relax', 'kg m-3s-1', 'SaltFlux_Relax', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & ldims=(/1/))

    CALL add_var(ocean_default_list, 'volume_global', ocean_state_diag%monitor%volume , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('volume', 'm3', 'volume', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & ldims=(/1/))

    CALL add_var(ocean_default_list, 'vorticity_global', ocean_state_diag%monitor%vorticity , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('vorticity', 's-1', 'vorticity', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & ldims=(/1/))

    CALL add_var(ocean_default_list, 'HeatFlux_ShortWave_global', ocean_state_diag%monitor%HeatFlux_ShortWave , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_ShortWave', 'W m-2', 'HeatFlux_ShortWave', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & ldims=(/1/))

    CALL add_var(ocean_default_list, 'HeatFlux_LongWave_global', ocean_state_diag%monitor%HeatFlux_LongWave , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_LongWave', 'W m-2', 'HeatFlux_LongWave_global', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & ldims=(/1/))

    CALL add_var(ocean_default_list, 'HeatFlux_Sensible_global', ocean_state_diag%monitor%HeatFlux_Sensible , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_Sensible', 'W m-2', 'HeatFlux_Sensible', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & ldims=(/1/))

    CALL add_var(ocean_default_list, 'HeatFlux_Latent_global', ocean_state_diag%monitor%HeatFlux_Latent , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('HeatFlux_Latent', 'W m-2', 'HeatFlux_Latent', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & ldims=(/1/))

    ! add 2d/3d diagnostic variables
    CALL add_var(ocean_default_list, 'rho', ocean_state_diag%rho , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('rho', 'kg m-3', 'sea water in-situ density', datatype_flt),&
      & grib2_var(10, 4, 16, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/), in_group=groups_oce_dde, lopenacc=.TRUE.)
    

    CALL add_var(ocean_default_list, 'rhopot', ocean_state_diag%rhopot , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('rhopot', 'kg m-3', 'sea water potential density', datatype_flt),&
      & grib2_var(10, 4, 19, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_dde, lopenacc=.TRUE.)
    

    ! why is this a diagnostic and in the restart list ???
    CALL add_var(ocean_default_list, 'grad_rho_PP_vert', ocean_state_diag%grad_rho_PP_vert, grid_unstructured_cell, &
      & za_depth_below_sea_half, &
      & t_cf_var('grad_rho_PP_vert','kg m-4','vertical density gradient at cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_diag)
    !is this usefull ? 
    CALL add_var(ocean_default_list, 'zgrad_rho', ocean_state_diag%zgrad_rho , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('zgrad_rho', 'kg m-4', 'vertical density gradiant', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag, lopenacc=.TRUE.)
    

!   CALL add_var(ocean_default_list, 'vt', ocean_state_diag%vt, grid_unstructured_edge, &
!     & za_depth_below_sea, &
!     & t_cf_var('vt','m/s','tangential velocity at edges', datatype_flt),&
!     & dflt_g2_decl_edge,&
!     & lopenacc=.TRUE., ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag)

    CALL add_var(ocean_default_list, 'h_e', ocean_state_diag%h_e, grid_unstructured_edge,&
      & za_surface, &
      & t_cf_var('h_e','m','surface height at edges', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,nblks_e/), in_group=groups_oce_diag, lopenacc=.TRUE.)
    

    ! thicknesses
    CALL add_var(ocean_default_list, 'thick_c', ocean_state_diag%thick_c,  &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('thick_c','m','fluid column thickness at cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/), in_group=groups_oce_diag, lopenacc=.TRUE.)
    

    CALL add_var(ocean_default_list, 'thick_e', ocean_state_diag%thick_e, &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('thick_e','m','fluid column thickness at edges', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,nblks_e/),in_group=groups_oce_diag, lopenacc=.TRUE.)
    

    CALL add_var(ocean_default_list, 'div_mass_flx_c', ocean_state_diag%div_mass_flx_c,&
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('div mass flux','kg s-1','divergence mass flux at cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag, lopenacc=.TRUE.)
    

    !check the unit volume or mass ?
    CALL add_var(ocean_default_list, 'mass_flux', ocean_state_diag%mass_flx_e, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('mass_flux','kg s-1',' mass flux at edges', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_default, &
      & lopenacc=.TRUE., initval=0.0_wp)
    
      
      CALL add_var(ocean_default_list, 'verticallyTotal_mass_flux_e', &
      & ocean_state_diag%verticallyTotal_mass_flux_e, grid_unstructured_edge,&
      & za_surface, &
      & t_cf_var('verticallyTotal_mass_flux_e','kg s-1','vertically integrated mass flux at edges', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,nblks_e/),in_group=groups_oce_diag)

    ! velocities
    CALL add_var(ocean_restart_list, 'w', ocean_state_diag%w, grid_unstructured_cell, &
      & za_depth_below_sea_half, &
      & t_cf_var('upward_sea_water_velocity','m s-1','vertical velocity at cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_diag","oce_default"), &
      & lrestart_cont=.TRUE., lopenacc=.TRUE.)
    

    !fixme is this usefull ?
    CALL add_var(ocean_default_list, 'w_old', ocean_state_diag%w_old, grid_unstructured_cell, &
      & za_depth_below_sea_half,&
      & t_cf_var('w_old','m s-1','vertical velocity at cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_diag,lopenacc=.TRUE.)
    

    IF ( ( vert_cor_type == 1 ) ) THEN
      CALL add_var(ocean_default_list, 'w_deriv', ocean_state_diag%w_deriv, grid_unstructured_cell, &
       & za_depth_below_sea_half, &
       & t_cf_var('w_der','m s-1','physical vertical velocity at cells for zstar', datatype_flt),&
       & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_diag","oce_default"))
    ENDIF

    CALL add_var(ocean_default_list, 'ssh', ocean_state_diag%ssh , &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,    &
         & t_cf_var('ssh', 'm', 'surface elevation at cell center', DATATYPE_FLT),&
         & grib2_var(10, 3, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
         & ldims=(/nproma,alloc_cell_blocks/), &
         & in_group=groups("oce_default"))

    ! tendencies per timestep used to calculate heat and salt transport in the moc diagnostic
    ! tendency of snow on ice
      CALL add_var(ocean_default_list, 'delta_snow', ocean_state_diag%delta_snow , &
       &         grid_unstructured_cell, za_surface,&
       & t_cf_var('tendency_of_snow_mass_expressed_as_heat_content', &
       & 'W m-2', 'tendency_of_snow_mass_expressed_as_heat_content', datatype_flt),&
       &         dflt_g2_decl_cell,&
       &         ldims=(/nproma,alloc_cell_blocks/), in_group=groups_oce_diag, lopenacc=.TRUE.)
      

    ! tendency of ice
      CALL add_var(ocean_default_list, 'delta_ice', ocean_state_diag%delta_ice, &
       &         grid_unstructured_cell, za_surface,&
       & t_cf_var('tendency_of_sea_ice_mass_expressed_as_heat_content', &
       & 'W m-2', 'tendency_of_sea_ice_mass_expressed_as_heat_content', datatype_flt),&
       &         dflt_g2_decl_cell,&
       &         ldims=(/nproma,alloc_cell_blocks/), in_group=groups_oce_diag, lopenacc=.TRUE.)
      

      CALL add_var(ocean_default_list, 'delta_thetao', ocean_state_diag%delta_thetao,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('tendency_of_sea_water_potential_temperature_expressed_as_heat_content',&
       & 'W m-2','tendency_of_sea_water_potential_temperature_expressed_as_heat_content', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/), in_group=groups_oce_diag, lopenacc=.TRUE.)
      

      CALL add_var(ocean_default_list, 'delta_so', ocean_state_diag%delta_so,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('tendency_of_sea_water_salinity_expressed_as_salt_content', &
       & 'kg m-2 s-1','tendency_of_sea_water_salinity_expressed_as_salt_content', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/), in_group=groups_oce_diag, lopenacc=.TRUE.)
      


    IF (diagnose_for_tendencies .or. GMRedi_configuration/=Cartesian_mixing) THEN

      CALL add_var(ocean_default_list, 'opottemptend', ocean_state_diag%opottemptend,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('tendency_of_sea_water_potential_temperature','K s-1','complete temperature tendency at cells', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

      CALL add_var(ocean_default_list, 'osalttend', ocean_state_diag%osalttend,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('tendency_of_sea_water_salinity','kg m-3 s-1','complete salinity tendency at cells', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

      CALL add_var(ocean_default_list, 'odensitytend', ocean_state_diag%odensitytend,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('odensitytend','kg m-3 s-1','complete density tendency at cells', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

    ENDIF ! diagnose_for_tendencies

      IF (eddydiag) THEN

      CALL add_var(ocean_default_list, 'w_prismcenter ', ocean_state_diag%w_prismcenter,&
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
       & t_cf_var('w prism center','m s-1','vertical velocity at prism center', DATATYPE_FLT32),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)
       CALL add_var(ocean_default_list, 'sigma0', ocean_state_diag%sigma0,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('sigma0','kg m-3','density anomaly', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'hflR', ocean_state_diag%hflR,&
       & grid_unstructured_cell, &
       & za_surface, &
       & t_cf_var('hflR','kg2ms-3','product of net heatflux and density', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'fwR', ocean_state_diag%fwR,&
       & grid_unstructured_cell, &
       & za_surface, &
       & t_cf_var('fwR','kg m-2s-1','product of fw flux and density', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'tauxu', ocean_state_diag%tauxu,&
       & grid_unstructured_cell, &
       & za_surface, &
       & t_cf_var('tauxu','Pa m s-1','product of x windstress and u-velocity', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'tauyv', ocean_state_diag%tauyv,&
       & grid_unstructured_cell, &
       & za_surface, &
       & t_cf_var('tauyv','Pa m s-1','product of y windstress and v-velocity', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'uT', ocean_state_diag%uT,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('uT','ms-1K','product of zonal velocity and temperature', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'uS', ocean_state_diag%uS,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('uS','m s-1 kg m-3','product of zonal velocity and salinity', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'uR', ocean_state_diag%uR,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('uR','ms-1 kg m-3','product of zonal velocity and density', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'uu', ocean_state_diag%uu,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('uu','m2s-2','square of zonal velocity', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'vT', ocean_state_diag%vT,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('vT','ms-1K','product of meridional velocity and temperature', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'vS', ocean_state_diag%vS,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('vS','m s-1 kg m-3','product of meridional velocity and salinity', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'vR', ocean_state_diag%vR,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('vR','ms-1 kg m-3','product of meridional velocity and density', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'vv', ocean_state_diag%vv,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('vv','m2s-2','square of meridional velocity', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)


       CALL add_var(ocean_default_list, 'wT', ocean_state_diag%wT,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('wT','ms-1K','product of vertical velocity and temperature', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'wS', ocean_state_diag%wS,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('wS','m s-1 kg m-3','product of vertical velocity and salinity', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'wR', ocean_state_diag%wR,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('wR','ms-1 kg m-3','product of vertical velocity and density', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'ww', ocean_state_diag%ww,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('ww','m2s-2','square of vertical velocity', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'RR', ocean_state_diag%rr,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('RR','kg2m-6','square of density', datatype_flt),&
       & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_eddy"))

       CALL add_var(ocean_default_list, 'SS', ocean_state_diag%ss,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('SS','kg2 m-6','square of salinity', datatype_flt),&
       & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_eddy"))

       CALL add_var(ocean_default_list, 'TT', ocean_state_diag%tt,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('TT','K2','square of temperature', datatype_flt),&
       & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_eddy"))

       CALL add_var(ocean_default_list, 'uv', ocean_state_diag%uv,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('uv','m2s-2','product of zonal velocity and meridional velocity', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'uw', ocean_state_diag%uw,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('uw','m2 s-2','product of zonal velocity and vertical velocity', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)

       CALL add_var(ocean_default_list, 'vw', ocean_state_diag%vw,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('vw','m2 s-2','product of meridional velocity and vertical velocity', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_eddy)
    ENDIF ! eddydiag

!   CALL add_var(ocean_restart_list, 'w_e', ocean_state_diag%w_e, grid_unstructured_cell, &
!     & za_depth_below_sea_half, &
!     & t_cf_var('w_e','m/s','vertical velocity at edges', datatype_flt),&
!     & dflt_g2_decl_edge,&
!     & lopenacc=.TRUE., ldims=(/nproma,n_zlev+1,nblks_e/),lrestart_cont=.TRUE.)
!   CALL add_var(ocean_default_list, 'w_prev', ocean_state_diag%w_prev, &
!     & grid_unstructured_edge, za_depth_below_sea_half, &
!     & t_cf_var('w_prev','m/s','vertical velocity at edges', datatype_flt),&
!     & dflt_g2_decl_edge,&
!     & lopenacc=.TRUE., ldims=(/nproma,n_zlev+1,nblks_e/))
    ! reconstructed u velocity component

    CALL add_var(ocean_default_list, 'u', ocean_state_diag%u, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('sea_water_x_velocity','m s-1','u zonal velocity component', datatype_flt),&
      & grib2_var(10, 4, 23, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_dde, lopenacc=.TRUE.)
    
    ! reconstructed v velocity component
    CALL add_var(ocean_default_list, 'v', ocean_state_diag%v, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('sea_water_y_velocity','m s-1','v meridional velocity component', datatype_flt),&
      & grib2_var(10, 4, 24, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_dde, lopenacc=.TRUE.)
    
    ! reconstrcuted velocity in cartesian coordinates
    !   CALL add_var(ocean_restart_list, 'p_vn', ocean_state_diag%p_vn, GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_SEA, &
    !   &            t_cf_var('p_vn','m/s','normal velocity in cartesian coordinates', datatype_flt),&
    !   &            dflt_g2_decl_cell,&
    !   &            lopenacc=.TRUE., ldims=(/nproma,n_zlev,alloc_cell_blocks/))
    ! integrated barotropic stream function

! check the unit
    CALL add_var(ocean_default_list, 'u_vint', ocean_state_diag%u_vint, grid_unstructured_cell, &
      & za_surface, &
      & t_cf_var('u_vint','m*m/s','barotropic zonal velocity', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/), &
      & in_group=groups_oce_dde, &
      & lopenacc=.TRUE.)
    
    CALL add_var(ocean_default_list, 'v_vint', ocean_state_diag%v_vint, grid_unstructured_cell, &
      & za_surface, &
      & t_cf_var('v_vint','m2 s-1','barotropic meridional velocity', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),&
      & in_group=groups_oce_dde)

    CALL add_var(ocean_default_list, 'ptp_vn', ocean_state_diag%ptp_vn, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('ptp_vn','m s-1','ptp_vn', &
      & datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),loutput=.TRUE., lopenacc=.TRUE., initval=0.0_wp)
    
    ! predicted vn normal velocity component

!     CALL add_var(ocean_restart_list, 'zlim', ocean_state_diag%zlim, &
!       & grid_unstructured_cell, za_depth_below_sea, &
!       & t_cf_var('zlim','1','zalesak limiter factor', &
!       & datatype_flt),&
!       & dflt_g2_decl_edge,&
!       & lopenacc=.TRUE., ldims=(/nproma,n_zlev,nblks_e/),loutput=.true., lrestart_cont=.false.)

    CALL add_var(ocean_default_list, 'vn_pred', ocean_state_diag%vn_pred, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('vn_pred','m s-1','predicted vn normal velocity component', &
      & datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag,lopenacc=.TRUE.)
    

    CALL add_var(ocean_default_list, 'vn_pred_ptp', ocean_state_diag%vn_pred_ptp, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('vn_pred_ptp','m s-1','transformed predicted vn normal velocity component', &
      & datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag)


    ! predicted vn normal velocity component
!    CALL add_var(ocean_restart_list, 'vn_impl_vert_diff', ocean_state_diag%vn_impl_vert_diff,&
!      & grid_unstructured_edge, za_depth_below_sea, &
!      & t_cf_var('vn_impl_vert_diff','m/s','predicted vn normal velocity component', &
!      & datatype_flt),&
!      & dflt_g2_decl_edge,&
!      & lopenacc=.TRUE., ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag,lrestart_cont=.TRUE.)

      CALL add_var(ocean_default_list, 'vn_time_weighted', ocean_state_diag%vn_time_weighted, &
      & grid_unstructured_edge, za_depth_below_sea, &
      & t_cf_var('vn_pred','m s-1','average vn normal velocity component', &
      & datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag,lopenacc=.TRUE.)
    

!     CALL add_var(ocean_default_list, 'w_time_weighted', ocean_state_diag%w_time_weighted, &
!       & grid_unstructured_cell, za_depth_below_sea_half, &
!       & t_cf_var('w_prev','m/s','vertical velocity at cells', datatype_flt),&
!       & dflt_g2_decl_cell,&
!       & lopenacc=.TRUE., ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_diag)

    ! vorticity
    CALL add_var(ocean_default_list, 'vort', ocean_state_diag%vort, &
      & grid_unstructured_vert, za_depth_below_sea, &
      & t_cf_var('vort','s-1','vorticity', datatype_flt),&
      & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_vertex),&
      & ldims=(/nproma,n_zlev,nblks_v/),in_group=groups_oce_diag,lopenacc=.TRUE.)
    

!    CALL add_var(ocean_restart_list, 'potential_vort_e', ocean_state_diag%potential_vort_e, &
!      & grid_unstructured_edge, za_depth_below_sea, &
!      & t_cf_var('vort_e','1/s','potential vorticity at edges', datatype_flt),&
!      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_unstructured, grid_edge),&
!      & lopenacc=.TRUE., ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_essentials"),lrestart_cont=.TRUE.)

!    CALL add_var(ocean_restart_list, 'potential vort_c', ocean_state_diag%potential_vort_c, &
!      & grid_unstructured_cell, za_depth_below_sea, &
!      & t_cf_var('vort_e','1/s','potential vorticity at cells', datatype_flt),&
!      & dflt_g2_decl_cell,&
!      & lopenacc=.TRUE., ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_essentials"),lrestart_cont=.TRUE.)

    ! kinetic energy component
    CALL add_var(ocean_default_list, 'kin', ocean_state_diag%kin, grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('kin','J','kinetic energy', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag, lopenacc=.TRUE.)
    

    ! gradient term
    CALL add_var(ocean_default_list, 'grad', ocean_state_diag%grad, grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('grad','','gradient', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag, lopenacc=.TRUE.)
    

    ! divergence component
!     CALL add_var(ocean_default_list, 'div', ocean_state_diag%div, grid_unstructured_cell, &
!       & za_depth_below_sea, &
!       & t_cf_var('div','','divergence', datatype_flt),&
!       & dflt_g2_decl_cell,&
!       & lopenacc=.TRUE., ldims=(/nproma,n_zlev,alloc_cell_blocks/),lrestart_cont=.FALSE.)

    ! pressures
    CALL add_var(ocean_default_list, 'press_hyd', ocean_state_diag%press_hyd, grid_unstructured_cell,&
      & za_depth_below_sea, t_cf_var('press_hyd','Pa','hydrostatic pressure', &
      & datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag,lrestart_cont=.TRUE., lopenacc=.TRUE.)
    

    CALL add_var(ocean_default_list, 'press_grad', ocean_state_diag%press_grad, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('press_grad','Pa m-1',' pressure gradient', &
      & datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag,lopenacc=.TRUE.)
    

    ! horizontal velocity advection
    CALL add_var(ocean_default_list, 'veloc_adv_horz', ocean_state_diag%veloc_adv_horz, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('veloc_adv_horz','fixme','horizontal velocity advection', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag,lopenacc=.TRUE.)
    

    ! vertical velocity advection
    CALL add_var(ocean_default_list, 'veloc_adv_vert', ocean_state_diag%veloc_adv_vert, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('veloc_adv_vert','fixme','vertical velocity advection', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag,lopenacc=.TRUE.)
    

    ! horizontal diffusion
    CALL add_var(ocean_default_list, 'laplacian_horz', ocean_state_diag%laplacian_horz, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('laplacian_horz','fixme','horizontal diffusion', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag,lrestart_cont=.FALSE., lopenacc=.TRUE.)
    

    ! vertical diffusion
    CALL add_var(ocean_default_list, 'laplacian_vert', ocean_state_diag%laplacian_vert, &
      & grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('laplacian_vert','fixme','vertical diffusion', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),lrestart_cont=.FALSE.)
    ! mixed layer depths
    CALL add_var(ocean_default_list, 'mld', ocean_state_diag%mld , grid_unstructured_cell,za_surface, &
      &          t_cf_var('ocean_mixed_layer_thickness', 'm', 'mixed layer depth', datatype_flt),&
      &          dflt_g2_decl_cell,&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_dde)

    ! CMIP6
    CALL add_var(ocean_default_list, 'mlotst', ocean_state_diag%mlotst , grid_unstructured_cell,za_surface, &
      &          t_cf_var('ocean_mixed_layer_thickness_defined_by_sigma_t', 'm', 'ocean_mixed_layer_thickness_defined_by_sigma_t', datatype_flt),&
      &          dflt_g2_decl_cell,&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_dde)
    ! CMIP6
    CALL add_var(ocean_default_list, 'mlotstsq', ocean_state_diag%mlotstsq , grid_unstructured_cell,za_surface, &
      &          t_cf_var('square_of_ocean_mixed_layer_thickness_defined_by_sigma_t', 'm2', 'square_of_ocean_mixed_layer_thickness_defined_by_sigma_t', datatype_flt),&
      &          dflt_g2_decl_cell,&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_dde)

    ! EERIE
    CALL add_var(ocean_default_list, 'mlotst10', ocean_state_diag%mlotst10 , grid_unstructured_cell,za_surface, &
      &          t_cf_var('mlotst10', 'm', 'ocean_mixed_layer_thickness_defined_by_sigma_t_10m', datatype_flt),&
      &          dflt_g2_decl_cell,&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_dde)
    ! EERIE
    CALL add_var(ocean_default_list, 'mlotst10sq', ocean_state_diag%mlotst10sq , grid_unstructured_cell,za_surface, &
      &          t_cf_var('mlotst10sq','m','square_of_ocean_mixed_layer_thickness_defined_by_sigma_t_10m', datatype_flt),&
      &          dflt_g2_decl_cell,&
      &          ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_dde)

    IF (diagnose_for_heat_content) THEN

    CALL add_var(ocean_default_list, 'global_heat_content', ocean_state_diag%monitor%global_heat_content , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_heat_content', 'J', 'global_heat_content', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    CALL add_var(ocean_default_list, 'global_heat_content_solid', ocean_state_diag%monitor%global_heat_content_solid , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_heat_content_solid', 'J', 'global_heat_content_solid', datatype_flt),&
      & dflt_g2_decl_lonlat,&
      & in_group=groups_oce_monitor,ldims=(/1/))

    ! heat content of liquid water
      CALL add_var(ocean_default_list, 'heat_content_liquid_water', ocean_state_diag%heat_content_liquid_water,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('heat_content_liquid_water','J m-2','heat_content_liquid_water', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_default)

    ! heat content of snow
    CALL add_var(ocean_default_list, 'heat_content_snow', ocean_state_diag%heat_content_snow , &
      &         grid_unstructured_cell, za_surface,&
      &         t_cf_var('heat_content_snow', 'J m-2', 'heat_conten_snow', datatype_flt),&
      &         dflt_g2_decl_cell,&
      &         ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_default)

   ! heat content of seaice
    CALL add_var(ocean_default_list, 'heat_content_seaice', ocean_state_diag%heat_content_seaice , &
      &         grid_unstructured_cell, za_surface,&
      &         t_cf_var('heat_content_seaice', 'J m-2', 'heat_content_seaice', datatype_flt),&
      &         dflt_g2_decl_cell,&
      &         ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_default)

   ! total heat content per column
    CALL add_var(ocean_default_list, 'heat_content_total', ocean_state_diag%heat_content_total , &
      &         grid_unstructured_cell, za_surface,&
      &         t_cf_var('heat_content_total', 'J m-2', 'heat_content_total', datatype_flt),&
      &         dflt_g2_decl_cell,&
      &         ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_default)
    ENDIF ! diagnose_for_heat_content

   ! swr fraction absorbed in the surface layer
    CALL add_var(ocean_default_list, 'swsum', ocean_state_diag%swsum , &
      &         grid_unstructured_cell, za_surface,&
      &         t_cf_var('swsum', '1', 'swsum', datatype_flt),&
      &         dflt_g2_decl_cell,&
      &         ldims=(/nproma,alloc_cell_blocks/), in_group=groups_oce_default, lopenacc=.TRUE.)
    

   ! (total subsurface ) heating due to absorption [J m-2]
    CALL add_var(ocean_default_list, 'heatabs', ocean_state_diag%heatabs , &
      &         grid_unstructured_cell, za_surface,&
      &         t_cf_var('heatabs', 'J m-2', 'heatabs', datatype_flt),&
      &         dflt_g2_decl_cell,&
      &         ldims=(/nproma,alloc_cell_blocks/), in_group=groups_oce_default, lopenacc=.TRUE.)
    

   ! heating due to fw flux [W m-2]
    CALL add_var(ocean_default_list, 'heatflux_rainevaprunoff', ocean_state_diag%heatflux_rainevaprunoff, &
      &         grid_unstructured_cell, za_surface,&
      &         t_cf_var('heatflux_rainevaprunof', 'W m-2', 'heatflux_rainevaprunof', datatype_flt),&
      &         dflt_g2_decl_cell,&
      &         ldims=(/nproma,alloc_cell_blocks/), in_group=groups_oce_default, lopenacc=.TRUE.)
    



   ! dummy2d
!    CALL add_var(ocean_default_list, 'dummy2d', ocean_state_diag%dummy2d , &
!      &         grid_unstructured_cell, za_surface,&
!      &         t_cf_var('dummy2d', '', 'dummy2d', datatype_flt),&
!      &         dflt_g2_decl_cell,&
!      &         lopenacc=.TRUE., ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_default)


   ! relative swr absorption factor
      CALL add_var(ocean_default_list, 'swrab', ocean_state_diag%swrab,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('swrab','1','swrab', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_default, lopenacc=.TRUE.)
      

   ! relative swr absorption factor from hamocc (LFB_BGC_OCE)
    IF ( lhamocc ) THEN
      CALL add_var(ocean_restart_list, 'swr_frac', ocean_state_diag%swr_frac,&
       & grid_unstructured_cell, &
       & za_depth_below_sea, &
       & t_cf_var('swr_frac','1','swr_frac', datatype_flt),&
       & dflt_g2_decl_cell,&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_default, lrestart_cont=.FALSE.,&
       & lopenacc=.TRUE.,initval=0.1_wp)
      
    ENDIF

   ! CMIP6 Net Rate of Absorption of Shortwave Energy in Ocean Layer
      CALL add_var(ocean_default_list, 'rsdoabsorb', ocean_state_diag%rsdoabsorb,&
        & grid_unstructured_cell, &
        & za_depth_below_sea, &
        & t_cf_var('rsdoabsorb','W m-2','Net Rate of Absorption of Shortwave Energy in Ocean Layer', &
        & datatype_flt), dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/), in_group=groups_oce_default, lopenacc=.TRUE.)
      

   ! dummy3d
!      CALL add_var(ocean_default_list, 'dummy3d', ocean_state_diag%dummy3d,&
!       & grid_unstructured_cell, &
!       & za_depth_below_sea, &
!       & t_cf_var('dummy3d','','dummy3d', &
!       & datatype_flt), dflt_g2_decl_cell,&
!       & lopenacc=.TRUE., ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_default)


    ! number of deepest convection layer
    CALL add_var(ocean_default_list, 'condep', ocean_state_diag%condep , grid_unstructured_cell, za_surface,&
      &         t_cf_var('condep', '1', 'convection depth index', datatype_flt),&
      &         dflt_g2_decl_cell,&
      &         ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_dde)
    IF (cfl_write) THEN
    CALL add_var(ocean_default_list, 'cfl_vert', ocean_state_diag%cfl_vert , &
      &          grid_unstructured_cell, za_depth_below_sea_half,&
      &          t_cf_var('cdf_vert', '1', 'vertical cfl relation', datatype_flt),&
      &          dflt_g2_decl_cell,&
      &          ldims=(/nproma,n_zlev+1,alloc_cell_blocks/), &
      &          in_group=groups_oce_diag)
    CALL add_var(ocean_default_list, 'cfl_horz', ocean_state_diag%cfl_horz, &
      &          grid_unstructured_edge, za_depth_below_sea,&
      &          t_cf_var('cfl_horz', '1', 'horizontal cfl relation', datatype_flt),&
      &          dflt_g2_decl_cell,&
      &          ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag)
    ENDIF
    
      !reconstrcuted velocity in cartesian coordinates
    ALLOCATE(ocean_state_diag%p_vn(nproma,n_zlev,alloc_cell_blocks), stat=ist)
    IF (ist/=success) THEN
      CALL finish(routine, 'allocation for p_vn at cells failed')
    END IF

    ALLOCATE(ocean_state_diag%p_vn_dual(nproma,n_zlev,nblks_v), stat=ist)
    IF (ist/=success) THEN
      CALL finish(routine, 'allocation for p_vn at verts failed')
    END IF

    !reconstrcuted velocity at edges in cartesian coordinates
    ALLOCATE(ocean_state_diag%p_mass_flux_sfc_cc(nproma,alloc_cell_blocks), stat=ist)
    IF (ist/=success) THEN
      CALL finish(routine, 'allocation for p_mass_flux_sfc_cc at cells failed')
    END IF
    ! set all values - incl. last block - of cartesian coordinates to zero (NAG compiler)
    ocean_state_diag%p_vn     (:,:,:)%x(1)=0.0_wp
    ocean_state_diag%p_vn     (:,:,:)%x(2)=0.0_wp
    ocean_state_diag%p_vn     (:,:,:)%x(3)=0.0_wp
    ocean_state_diag%p_vn_dual(:,:,:)%x(1)=0.0_wp
    ocean_state_diag%p_vn_dual(:,:,:)%x(2)=0.0_wp
    ocean_state_diag%p_vn_dual(:,:,:)%x(3)=0.0_wp
    ocean_state_diag%p_mass_flux_sfc_cc(:,:)%x(1)=0.0_wp
    ocean_state_diag%p_mass_flux_sfc_cc(:,:)%x(2)=0.0_wp
    ocean_state_diag%p_mass_flux_sfc_cc(:,:)%x(3)=0.0_wp

    !$ACC ENTER DATA COPYIN(ocean_state_diag%p_vn, ocean_state_diag%p_vn_dual)

    !remapped velocity at cell edges
!     ALLOCATE(ocean_state_diag%ptp_vn(nproma,n_zlev,nblks_e), stat=ist)
!     IF (ist/=success) THEN
!       CALL finish(routine, 'allocation for ptp_vn at edges failed')
!     END IF

!     CALL add_var(ocean_default_list,'temp_insitu',ocean_state_diag%temp_insitu,grid_unstructured_cell,&
!       & za_depth_below_sea, &
!       & t_cf_var('temp_insitu', 'K', 'in situ temperature', datatype_flt),&
!       & dflt_g2_decl_cell,&
!       & lopenacc=.TRUE., ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

    CALL add_var(ocean_default_list,'zos_square',ocean_state_diag%zos_square,grid_unstructured_cell,&
      & za_surface, &
      & t_cf_var('zos_square', 'm2', 'square of sea surface hight', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_default)

     CALL add_var(ocean_default_list,'Rossby_Radius',ocean_state_diag%Rossby_Radius,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('Rossby_Radius', 'm', 'Rossby Radius', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_diag)

    CALL add_var(ocean_default_list,'Richardson_Number',ocean_state_diag%Richardson_Number,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('richardson_number_in_sea_water', '1', 'Richardson Number', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag, lopenacc=.TRUE.)
    


    IF (.false.) THEN
    CALL add_var(ocean_default_list,'Buoyancy_Freq',ocean_state_diag%Buoyancy_Freq,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('Buoyancy_Freq', '1 s-1', 'Buoyancy Frequency', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

   CALL add_var(ocean_default_list,'Wavespeed_baroclinic',ocean_state_diag%Wavespeed_baroclinic,grid_unstructured_cell,&
      & za_surface, &
      & t_cf_var('Wavespeed_baroclinic', 'm s-1', 'Baroclinic wavespeed', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_diag)
  ENDIF

   ! masks for northern and southern part of the earth {{{
   CALL add_var(ocean_default_list,'northernHemisphere',ocean_state_diag%northernHemisphere,&
      & grid_unstructured_cell,&
      & za_surface, &
      & t_cf_var('northern_hemisphere', '1', 'northern hemisphere ', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/))

   CALL add_var(ocean_default_list,'southernHemisphere',ocean_state_diag%southernHemisphere, &
      & grid_unstructured_cell,&
      & za_surface, &
      & t_cf_var('southern_hemisphere', '1', 'southern hemisphere ', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/))

   owned_cells    => patch_2d%cells%owned
   DO blockNo = owned_cells%start_block, owned_cells%end_block
     CALL get_index_range(owned_cells, blockNo, start_cell_index, end_cell_index)
     DO jc =  start_cell_index, end_cell_index
       IF (patch_2d%cells%center(jc,blockNo)%lat > equator) THEN
         ocean_state_diag%northernHemisphere(jc,blockNo) = 1.0_wp
         ocean_state_diag%southernHemisphere(jc,blockNo) = 0.0_wp
       ELSE
         ocean_state_diag%northernHemisphere(jc,blockNo) = 0.0_wp
         ocean_state_diag%southernHemisphere(jc,blockNo) = 1.0_wp
       END IF
     END DO
   END DO

   

!     CALL add_var(ocean_restart_list,'temp_horDiffused',ocean_state_diag%temp_horizontally_diffused, grid_unstructured_cell,&
!       & za_depth_below_sea, &
!       & t_cf_var('temp_insitu', 'K', 'horizonatlly diffused temperature', datatype_flt),&
!       & dflt_g2_decl_cell,&
!       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag,lrestart_cont=.TRUE.)

    !--------------------------------------------------------------------------
    ! MOCs on a zonal 1deg grid
    CALL add_var(ocean_default_list, 'global_moc',ocean_state_diag%global_moc,    &
      & GRID_ZONAL, za_depth_below_sea,&
      & t_cf_var('global_moc','kg s-1','global meridional overturning', datatype_flt), &
      & grib2_var(255, 255, 147, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/n_zlev,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)
    CALL add_var(ocean_default_list, 'atlantic_moc',ocean_state_diag%atlantic_moc,    &
      & GRID_ZONAL, za_depth_below_sea,&
      & t_cf_var('atlantic_moc','kg s-1','atlantic meridional overturning', datatype_flt), &
      & grib2_var(255, 255, 148, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/n_zlev,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)
    CALL add_var(ocean_default_list, 'pacific_moc',ocean_state_diag%pacific_moc,    &
      & GRID_ZONAL, za_depth_below_sea,&
      & t_cf_var('pacific_moc','kg s-1','indopacific meridional overturning', datatype_flt), &
      & grib2_var(255, 255, 149, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/n_zlev,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)

    ! Implied ocean heat transport
    CALL add_var(ocean_default_list, 'global_hfl',ocean_state_diag%global_hfl,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('global_hfl','W','global implied heat transport', datatype_flt), &
      & grib2_var(255, 255, 147, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)
    CALL add_var(ocean_default_list, 'atlantic_hfl',ocean_state_diag%atlantic_hfl,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('atlantic_hfl','W','atlantic implied heat transport', datatype_flt), &
      & grib2_var(255, 255, 148, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)

    CALL add_var(ocean_default_list, 'pacific_hfl',ocean_state_diag%pacific_hfl,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('pacific_hfl','W','indopacific implied heat transport', datatype_flt), &
      & grib2_var(255, 255, 149, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)

    ! Implied ocean fw transport
    CALL add_var(ocean_default_list, 'global_wfl',ocean_state_diag%global_wfl,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('global_wfl','m3s-1','global implied freshwater transport', datatype_flt), &
      & grib2_var(255, 255, 147, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)
    CALL add_var(ocean_default_list, 'atlantic_wfl',ocean_state_diag%atlantic_wfl,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('atlantic_wfl','m3s-1','atlantic implied freshwater transport', datatype_flt), &
      & grib2_var(255, 255, 148, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)
    CALL add_var(ocean_default_list, 'pacific_wfl',ocean_state_diag%pacific_wfl,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('pacific_wfl','m3 s-1','indopacific implied freshwater transport', datatype_flt), &
      & grib2_var(255, 255, 149, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)

    ! hfbasin ocean heat transport
    CALL add_var(ocean_default_list, 'global_hfbasin',ocean_state_diag%global_hfbasin,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('global_hfbasin','W','global northward ocean heat transport', datatype_flt), &
      & grib2_var(255, 255, 147, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)
    CALL add_var(ocean_default_list, 'atlantic_hfbasin',ocean_state_diag%atlantic_hfbasin,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('atlantic_hfbasin','W','atlantic northward ocean heat transport', datatype_flt), &
      & grib2_var(255, 255, 148, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)
    CALL add_var(ocean_default_list, 'pacific_hfbasin',ocean_state_diag%pacific_hfbasin,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('pacific_hfbasin','W','indopacific northward ocean heat transport', datatype_flt), &
      & grib2_var(255, 255, 149, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)

    ! hfbasin ocean salt transport
    CALL add_var(ocean_default_list, 'global_sltbasin',ocean_state_diag%global_sltbasin,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('global_sltbasin','kg s-1','global northward ocean salt transport', datatype_flt), &
      & grib2_var(255, 255, 147, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)
    CALL add_var(ocean_default_list, 'atlantic_sltbasin',ocean_state_diag%atlantic_sltbasin,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('atlantic_sltbasin','kg s-1','atlantic northward ocean salt transport', datatype_flt), &
      & grib2_var(255, 255, 148, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)
    CALL add_var(ocean_default_list, 'pacific_sltbasin',ocean_state_diag%pacific_sltbasin,    &
      & GRID_ZONAL, za_surface,&
      & t_cf_var('pacific_sltbasin','kg s-1','indopacific northward ocean salt transport', datatype_flt), &
      & grib2_var(255, 255, 149, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/1,180/),in_group=groups_oce_moc,&
      & loutput=.TRUE.)



!       CALL add_var(ocean_restart_list, 'tracers'//TRIM(var_suffix), ocean_state_prog%tracer , &
!         & grid_unstructured_cell, za_depth_below_sea, &
!         & t_cf_var('tracers'//TRIM(var_suffix), '', '1:temperature 2:salinity', &
!         & DATATYPE_FLT64),&
!         & dflt_g2_decl_cell,&
!         & lopenacc=.TRUE., ldims=(/nproma,n_zlev,alloc_cell_blocks,no_tracer/), &
!         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
!
!       ! Reference to individual tracer, for I/O
!       ALLOCATE(ocean_state_prog%tracer_ptr(no_tracer))
!       DO jtrc = 1,no_tracer
!         CALL add_ref( ocean_restart_list, 'tracers'//TRIM(var_suffix),              &
!           & oce_tracer_names(jtrc),                 &
!           & ocean_state_prog%tracer_ptr(jtrc)%p,                             &
!           & grid_unstructured_cell, za_depth_below_sea,               &
!           & t_cf_var(TRIM(oce_tracer_names(jtrc))//TRIM(var_suffix), &
!           & oce_tracer_units(jtrc), &
!           & oce_tracer_longnames(jtrc), DATATYPE_FLT64), &
!           & dflt_g2_decl_cell,&
!           & ref_idx=jtrc, &
!           & ldims=(/nproma,n_zlev,alloc_cell_blocks/))
!       END DO
!
!       ! use of the ocean_tracers structure
!       ALLOCATE(ocean_state_prog%ocean_tracers(no_tracer))
!       DO jtrc = 1,no_tracer
!         ocean_state_prog%ocean_tracers(jtrc)%concentration =>  ocean_state_prog%tracer(:,:,:,jtrc)
!       ENDDO
!

   IF(GMRedi_configuration/=Cartesian_Mixing)THEN

!     CALL add_var(ocean_restart_list, 'vn_bolus', ocean_state_diag%vn_bolus, &
!       & grid_unstructured_edge, za_depth_below_sea, &
!       & t_cf_var('vn_bolus','m/s',' bolus vn normal velocity component', &
!       & datatype_flt),&
!       & dflt_g2_decl_edge,&
!       & lopenacc=.TRUE., ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_diag,lrestart_cont=.TRUE.)

!       CALL add_var(ocean_restart_list, 'w_bolus', ocean_state_diag%w_bolus, grid_unstructured_cell, &
!        & za_depth_below_sea, &
!        & t_cf_var('w bolus','m/s','vertical bolus velocity at prism center', DATATYPE_FLT32),&
!        & dflt_g2_decl_cell,&
!        & lopenacc=.TRUE., ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag,lrestart_cont=.TRUE.)

     CALL add_var(ocean_default_list, 'rhopo_GM', ocean_state_diag%rho_GM , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('rhopot', 'kg m-3', 'potential density', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

!_by_Oliver
     CALL add_var(ocean_default_list, 'kappa_GM', ocean_state_diag%kappa_GM , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('kappa_GM', 'm2 s-1', 'kappa_GM', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

! what is this ?
     CALL add_var(ocean_default_list, 'Nratio', ocean_state_diag%Nratio , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('Nratio', '', 'Nratio', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

     CALL add_var(ocean_default_list, 'N2_ref', ocean_state_diag%N2_ref , grid_unstructured_cell,&
      & za_surface, &
      & t_cf_var('N2_ref', '1/s', 'N2_ref', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_diag)

!are these needed ?
    CALL add_var(ocean_default_list,'osaltGMRedi',ocean_state_diag%osaltGMRedi,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('osaltGMRedi', 'kg m-3 s-1', 'osaltGMRedi', datatype_flt),&
      & dflt_g2_decl_cell,&
     & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

   CALL add_var(ocean_default_list,'opottempGMRedi',ocean_state_diag%opottempGMRedi,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('opottempGMRedi', 'K s-1', 'opottempGMRedi', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

  CALL add_var(ocean_default_list,'div_of_GMRedi_flux',ocean_state_diag%div_of_GMRedi_flux,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('temp_insitu', 'm', 'div_of_GMRedi_flux', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

   CALL add_var(ocean_default_list,'div_of_GMRedi_flux_horizontal',&
      &ocean_state_diag%div_of_GMRedi_flux_horizontal,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('temp_insitu', 'm', 'div_of_GMRedi_flux_horizontal', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

   CALL add_var(ocean_default_list,'div_of_GMRedi_flux_vertical',&
      &ocean_state_diag%div_of_GMRedi_flux_vertical,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('temp_insitu', 'm', 'div_of_GMRedi_flux_vertical', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

   CALL add_var(ocean_default_list,'vertical_mixing_coeff_GMRedi_implicit',&
   &ocean_state_diag%vertical_mixing_coeff_GMRedi_implicit,grid_unstructured_cell,&
      & za_depth_below_sea_half, &
      & t_cf_var('temp_insitu', 'm', 'vertical_mixing_coeff_GMRedi_implicit', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_diag)


   CALL add_var(ocean_default_list,'div_of_GM_flux',ocean_state_diag%div_of_GM_flux,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('temp_insitu', 'm', 'div_of_GM_flux', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

   CALL add_var(ocean_default_list,'div_of_Redi_flux',ocean_state_diag%div_of_Redi_flux,grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('temp_insitu', 'm', 'div_of_Redi_flux', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_diag)

   CALL add_var(ocean_default_list, 'GMRedi_flux_horz',ocean_state_diag%GMRedi_flux_horz, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('GMRedi_flux_horz', '', '1:temperature 2:salinity', &
        & DATATYPE_FLT),&
        & dflt_g2_decl_edge,&
        & ldims=(/nproma,n_zlev,nblks_e,no_tracer/), &
        & lcontainer=.TRUE., initval=0.0_wp)

   CALL add_var(ocean_default_list, 'GMRedi_flux_vert',ocean_state_diag%GMRedi_flux_vert, &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('GMRedi_flux_vert', '', '1:temperature 2:salinity', &
        & DATATYPE_FLT),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev+1,alloc_cell_blocks,no_tracer/), &
        & lcontainer=.TRUE., initval=0.0_wp)
    ENDIF

    SELECT CASE (test_mode)

      CASE (103:113) ! testbed_div

        CALL add_var(ocean_default_list,'div_model',ocean_state_diag%div_model,grid_unstructured_cell,&
          & za_depth_below_sea, &
          & t_cf_var('div_model', 'm', 'div_model', datatype_flt),&
          & dflt_g2_decl_cell,&
          & ldims=(/nproma,n_zlev,alloc_cell_blocks/),lrestart_cont=.FALSE., loutput=.TRUE.)

       CALL add_var(ocean_default_list,'div_diff',ocean_state_diag%div_diff,grid_unstructured_cell,&
          & za_depth_below_sea, &
          & t_cf_var('div_diff', 'm', 'div_diff', datatype_flt),&
          & dflt_g2_decl_cell,&
          & ldims=(/nproma,n_zlev,alloc_cell_blocks/),lrestart_cont=.FALSE., loutput=.TRUE.)

       CALL add_var(ocean_default_list,'divPtP',ocean_state_diag%divPtP,grid_unstructured_cell,&
          & za_depth_below_sea, &
          & t_cf_var('divPtP', 'm', 'divPtP', datatype_flt),&
          & dflt_g2_decl_cell,&
          & ldims=(/nproma,n_zlev,alloc_cell_blocks/),lrestart_cont=.FALSE., loutput=.TRUE.)

       CALL add_var(ocean_default_list,'divPtP_diff',ocean_state_diag%divPtP_diff,grid_unstructured_cell,&
          & za_depth_below_sea, &
          & t_cf_var('divPtP_diff', 'm', 'divPtP_diff', datatype_flt),&
          & dflt_g2_decl_cell,&
          & ldims=(/nproma,n_zlev,alloc_cell_blocks/),lrestart_cont=.FALSE., loutput=.TRUE.)

    END SELECT

  END SUBROUTINE construct_hydro_ocean_diag


  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of diagnostic hydrostatic ocean state.
  !
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_hydro_ocean_diag(ocean_state_diag)

    TYPE(t_hydro_ocean_diag), INTENT(inout) :: ocean_state_diag

    ! local variables

    INTEGER :: ist

    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_state:destruct_hydro_ocean_diag'

    DEALLOCATE(ocean_state_diag%p_vn, stat=ist)
    IF (ist/=success) THEN
      CALL finish(routine, 'deallocation for p_vn failed')
    END IF
    DEALLOCATE(ocean_state_diag%p_vn_dual, stat=ist)
    IF (ist/=success) THEN
      CALL finish(routine, 'deallocation for p_vn_dual failed')
    END IF
!     DEALLOCATE(ocean_state_diag%ptp_vn, stat=ist)
!     IF (ist/=success) THEN
!       CALL finish(routine, 'deallocation for ptp_vn failed')
!     END IF

  END SUBROUTINE destruct_hydro_ocean_diag
 !-------------------------------------------------------------------------
  !>
  !!               Allocation of components for 3dim ocean nudging.
  !!               Initialization of components with zero.
  !
  !!
!<Optimize:inUse>
  SUBROUTINE construct_ocean_nudge(patch_2d, ocean_nudge)

    TYPE(t_patch),TARGET, INTENT(in)                :: patch_2d
    TYPE(t_ocean_nudge), TARGET,INTENT(inout)   :: ocean_nudge

    ! local variables

    INTEGER ::  ist  !, jtrc
    INTEGER ::  alloc_cell_blocks, nblks_e, nblks_v

    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_state:construct_ocean_nudge'
    INTEGER :: datatype_flt
    LOGICAL :: groups_oce_nudge(MAX_GROUPS)
    TYPE(t_grib2_var) :: dflt_g2_decl_cell

    dflt_g2_decl_cell = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell)

    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)

    groups_oce_nudge = groups("oce_nudge")
    !-------------------------------------------------------------------------
    CALL message(routine, 'construct  ocean nudge...')

    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v
    

    !$ACC ENTER DATA COPYIN(ocean_nudge)

    ! allocation of 3-dim tracer relaxation:
    IF (no_tracer>=1 .AND. type_3dimrelax_temp >0) THEN
      CALL add_var(ocean_default_list,'data_3dimRelax_Temp',ocean_nudge%data_3dimRelax_Temp,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('data_3dimRelax_Temp','C','relaxation temperatrure', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_nudge,lopenacc = .TRUE.)
      
      CALL add_var(ocean_default_list,'forc_3dimRelax_Temp',ocean_nudge%forc_3dimRelax_Temp,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('forc_3dimRelax_Temp','C','forc_3dimRelax_Temp', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_nudge,lopenacc = .TRUE.)
      
      CALL add_var(ocean_default_list,'relax_3dim_coefficient',ocean_nudge%relax_3dim_coefficient,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('relax_3dim_coefficient','m s-1','relax_3dim_coefficient', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_nudge)
!       ocean_state_aux%relax_3dim_coefficient(:,:,:) = 1.0_wp 
    END IF
    IF (no_tracer==2 .AND. type_3dimrelax_salt >0) THEN
      CALL add_var(ocean_default_list,'data_3dimRelax_Salt',ocean_nudge%data_3dimRelax_Salt,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('data_3dimRelax_Salt','psu','data_3dimRelax_Salt', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_nudge, lopenacc = .TRUE.)
      
      CALL add_var(ocean_default_list,'forc_3dimRelax_Salt',ocean_nudge%forc_3dimRelax_Salt,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('forc_3dimRelax_Salt','psu','forc_3dimRelax_Salt', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_nudge, lopenacc = .TRUE.)
      
    END IF

  END SUBROUTINE construct_ocean_nudge

  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of auxilliary hydrostatic ocean state.
  !
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_ocean_nudge(ocean_nudge)

    TYPE(t_ocean_nudge), INTENT(inout)      :: ocean_nudge

    ! local variables

    INTEGER :: ist

    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_state:destruct_hydro_ocean_aux'

    DEALLOCATE(ocean_nudge%data_3dimRelax_Temp, stat=ist)
    IF (ist/=success) THEN
      CALL finish(routine,'deallocation of data_3dimRelax_Temp failed')
    END IF

  END SUBROUTINE destruct_ocean_nudge
  !-------------------------------------------------------------------------




!-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean auxiliary state.
  !!               Initialization of components with zero.
  !
  !!
!<Optimize:inUse>
  SUBROUTINE construct_hydro_ocean_aux(patch_2d, ocean_state_aux, ocean_transport_state)

    TYPE(t_patch),TARGET, INTENT(in)                :: patch_2d
    TYPE(t_hydro_ocean_aux), TARGET,INTENT(inout)   :: ocean_state_aux
    TYPE(t_ocean_transport_state), TARGET,INTENT(inout)  :: ocean_transport_state

    ! local variables
    INTEGER ::  ist, alloc_cell_blocks, nblks_e, nblks_v
    CHARACTER(*), PARAMETER :: routine = 'mo_ocean_state:construct_hydro_ocean_aux'
    INTEGER :: datatype_flt
    LOGICAL, DIMENSION(MAX_GROUPS) :: groups_oce_diag, groups_oce_aux
    TYPE(t_grib2_var) :: dflt_g2_decl_cell, dflt_g2_decl_edge

    dflt_g2_decl_cell = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell)
    dflt_g2_decl_edge = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge)

    groups_oce_diag = groups("oce_diag")
    groups_oce_aux = groups("oce_aux")

    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)

    !-------------------------------------------------------------------------
    CALL message(routine, 'construct hydro ocean auxiliary state...')

    ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v

    ! allocation for Adam-Bashford time stepping
    CALL add_var(ocean_restart_list,'g_n',ocean_state_aux%g_n, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('g_n','fixme','g_n', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_aux,loutput=.TRUE.,lrestart_cont=.TRUE., lopenacc=.TRUE.)
    

    CALL add_var(ocean_restart_list,'g_nm1',ocean_state_aux%g_nm1, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('g_nm1','fixme','g_nm1', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_aux,loutput=.TRUE.,lrestart_cont=.TRUE., lopenacc=.TRUE.)
    

    CALL add_var(ocean_default_list,'g_nimd',ocean_state_aux%g_nimd, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('g_nimd','fixme','g_nimd', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_aux,loutput=.TRUE.,lrestart_cont=.TRUE., lopenacc=.TRUE.)
    

    CALL add_var(ocean_default_list,'p_rhs_sfc_eq',ocean_state_aux%p_rhs_sfc_eq, grid_unstructured_cell,&
      & za_surface, t_cf_var('p_rhs_sfc_eq','fixme','p_rhs_sfc_eq', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_aux, lopenacc=.TRUE.)
    

    ! allocation for boundary conditions
    CALL add_var(ocean_default_list,'bc_top_u',ocean_state_aux%bc_top_u, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_u','fixme','bc_top_u', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_aux, lopenacc=.TRUE.)
    

    CALL add_var(ocean_default_list,'bc_top_v',ocean_state_aux%bc_top_v, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_v','fixme','bc_top_v', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_aux,lrestart_cont=.TRUE., lopenacc=.TRUE.)
    

    CALL add_var(ocean_default_list,'bc_top_vn',ocean_state_aux%bc_top_vn, grid_unstructured_edge,&
      & za_surface, t_cf_var('bc_top_vn','fixme','bc_top_vn', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,nblks_e/),in_group=groups_oce_aux, lopenacc=.TRUE.)
    

    CALL add_var(ocean_default_list,'bc_top_WindStress',ocean_state_aux%bc_top_WindStress, grid_unstructured_edge,&
      & za_surface, t_cf_var('bc_top_WindStress','fixme','bc_top_WindStress', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,nblks_e/),in_group=groups_oce_aux)

    CALL add_var(ocean_default_list,'bc_bot_vn',ocean_state_aux%bc_bot_vn, grid_unstructured_edge,&
      & za_surface, t_cf_var('bc_bot_vn','fixme','bc_bot_vn', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,nblks_e/),in_group=groups_oce_aux, lopenacc=.TRUE.)
    

    CALL add_var(ocean_default_list,'bc_bot_w',ocean_state_aux%bc_bot_w, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_bot_w','fixme','bc_bot_w', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_aux)

    CALL add_var(ocean_default_list,'bc_top_w',ocean_state_aux%bc_top_w, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_w','fixme','bc_top_w', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_aux)

    CALL add_var(ocean_default_list,'bc_tides_potential',ocean_state_aux%bc_tides_potential, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_tides_potential','fixme','bc_tides_potential', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_aux, initval=0.0_wp)

!     this is used by default in the calculation of the total potential
!    IF (use_tides_SAL) THEN
!      CALL add_var(ocean_default_list,'bc_tides_load',ocean_state_aux%bc_tides_load, grid_unstructured_cell,&
!        & za_surface, t_cf_var('bc_tides_load','','', datatype_flt),&
!        & dflt_g2_decl_cell,&
!        & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_aux,lrestart_cont=.FALSE., lrestart=.FALSE.)
!      ocean_state_aux%bc_tides_load = 0.0_wp
!    ENDIF

    CALL add_var(ocean_default_list,'bc_total_top_potential',ocean_state_aux%bc_total_top_potential, grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_total_top_potential','fixme','bc_total_top_potential', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_aux, lopenacc=.TRUE., initval=0.0_wp)
    

    CALL add_var(ocean_default_list,'bc_SAL_potential',ocean_state_aux%bc_SAL_potential, &
      & grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_SAL_potential','fixme','bc_SAL_potential', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_aux, initval=0.0_wp)

    CALL add_var(ocean_default_list,'bc_bot_tracer',ocean_state_aux%bc_bot_tracer,grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_bot_tracer','fixme','bc_bot_tracer', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks,no_tracer/),in_group=groups_oce_aux)
    CALL add_var(ocean_default_list,'bc_top_tracer',ocean_state_aux%bc_top_tracer,grid_unstructured_cell,&
      & za_surface, t_cf_var('bc_top_tracer','fixme','bc_top_tracer', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks,no_tracer/),in_group=groups_oce_aux)

    ALLOCATE(ocean_state_aux%bc_top_veloc_cc(nproma,alloc_cell_blocks), stat=ist)
    IF (ist/=success) THEN
      CALL finish(routine,'allocation of top boundary cond cc failed')
    END IF

    ocean_state_aux%bc_top_veloc_cc(:,:)%x(1) = 0.0_wp
    ocean_state_aux%bc_top_veloc_cc(:,:)%x(2) = 0.0_wp
    ocean_state_aux%bc_top_veloc_cc(:,:)%x(3) = 0.0_wp

    !$ACC ENTER DATA COPYIN(ocean_state_aux%bc_top_veloc_cc)


   IF(GMRedi_configuration/=Cartesian_mixing)THEN

     ALLOCATE(ocean_state_aux%slopes(nproma,n_zlev,alloc_cell_blocks), stat=ist)
     IF (ist/=success) THEN
      CALL finish(routine, 'allocation for slopes at cells failed')
     END IF
     ocean_state_aux%slopes    (:,:,:)%x(1)=0.0_wp
     ocean_state_aux%slopes    (:,:,:)%x(2)=0.0_wp
     ocean_state_aux%slopes    (:,:,:)%x(3)=0.0_wp
     ALLOCATE(ocean_state_aux%PgradTemperature_horz_center(nproma,n_zlev,alloc_cell_blocks), stat=ist)
     IF (ist/=success) THEN
      CALL finish(routine, 'allocation for PgradTemperature_horz_center at cells failed')
     END IF
     ocean_state_aux%PgradTemperature_horz_center    (:,:,:)%x(1)=0.0_wp
     ocean_state_aux%PgradTemperature_horz_center    (:,:,:)%x(2)=0.0_wp
     ocean_state_aux%PgradTemperature_horz_center    (:,:,:)%x(3)=0.0_wp
     ALLOCATE(ocean_state_aux%PgradTracer_horz_center(nproma,n_zlev,alloc_cell_blocks), stat=ist)
     IF (ist/=success) THEN
      CALL finish(routine, 'allocation for PgradTracer_horz_center at cells failed')
     END IF
     ocean_state_aux%PgradTracer_horz_center    (:,:,:)%x(1)=0.0_wp
     ocean_state_aux%PgradTracer_horz_center    (:,:,:)%x(2)=0.0_wp
     ocean_state_aux%PgradTracer_horz_center    (:,:,:)%x(3)=0.0_wp

     ALLOCATE(ocean_state_aux%PgradSalinity_horz_center(nproma,n_zlev,alloc_cell_blocks), stat=ist)
     IF (ist/=success) THEN
      CALL finish(routine, 'allocation for PgradSalinity_horz_center at cells failed')
     END IF
     ocean_state_aux%PgradSalinity_horz_center       (:,:,:)%x(1)=0.0_wp
     ocean_state_aux%PgradSalinity_horz_center       (:,:,:)%x(2)=0.0_wp
     ocean_state_aux%PgradSalinity_horz_center       (:,:,:)%x(3)=0.0_wp

     ALLOCATE(ocean_state_aux%PgradDensity_horz_center(nproma,n_zlev,alloc_cell_blocks), stat=ist)
     IF (ist/=success) THEN
      CALL finish(routine, 'allocation for PgradSalinity_horz_center at cells failed')
     END IF
     ocean_state_aux%PgradDensity_horz_center       (:,:,:)%x(1)=0.0_wp
     ocean_state_aux%PgradDensity_horz_center       (:,:,:)%x(2)=0.0_wp
     ocean_state_aux%PgradDensity_horz_center       (:,:,:)%x(3)=0.0_wp

     ALLOCATE(ocean_state_aux%diagnose_Redi_flux_temp(nproma,n_zlev,alloc_cell_blocks), stat=ist)
     IF (ist/=success) THEN
      CALL finish(routine, 'allocation for diagnose_Redi_flux_temp at cells failed')
     END IF
     ocean_state_aux%diagnose_Redi_flux_temp       (:,:,:)%x(1)=0.0_wp
     ocean_state_aux%diagnose_Redi_flux_temp       (:,:,:)%x(2)=0.0_wp
     ocean_state_aux%diagnose_Redi_flux_temp       (:,:,:)%x(3)=0.0_wp

    ALLOCATE(ocean_state_aux%diagnose_Redi_flux_sal(nproma,n_zlev,alloc_cell_blocks), stat=ist)
     IF (ist/=success) THEN
      CALL finish(routine, 'allocation for diagnose_Redi_flux_sal at cells failed')
     END IF
     ocean_state_aux%diagnose_Redi_flux_sal       (:,:,:)%x(1)=0.0_wp
     ocean_state_aux%diagnose_Redi_flux_sal       (:,:,:)%x(2)=0.0_wp
     ocean_state_aux%diagnose_Redi_flux_sal       (:,:,:)%x(3)=0.0_wp

     CALL add_var(ocean_default_list,'DerivTemperature_vert',ocean_state_aux%DerivTemperature_vert_center,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('DerivTemperature_vert','fixme','DerivTemperature_vert', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_aux)
        ! & lopenacc=.TRUE., ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_aux,loutput=.FALSE.)

     CALL add_var(ocean_default_list,'DerivSalinity_vert',ocean_state_aux%DerivSalinity_vert_center,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('DerivSalinity_vert','fixme','DerivSalinity_vert', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_aux)

     CALL add_var(ocean_default_list,'DerivDensity_vert',ocean_state_aux%DerivDensity_vert_center,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('DerivDensity_vert','fixme','DerivDensity_vert', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_aux)

     CALL add_var(ocean_default_list,'DerivTracer_vert_center',ocean_state_aux%DerivTracer_vert_center,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('DerivTracer_vert:center','fixme','DerivTracer_vert_center', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_aux)
        ! & lopenacc=.TRUE., ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_aux,loutput=.FALSE.)

    CALL add_var(ocean_default_list, 'tracer_deriv_vert', ocean_state_aux%tracer_deriv_vert, &
      & grid_unstructured_cell, za_depth_below_sea_half, &
      & t_cf_var('tracer_deriv_vert','fixme','tracer_deriv_vert', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_diag)

    CALL add_var(ocean_default_list, 'temperature_deriv_vert', ocean_state_aux%temperature_deriv_vert, &
      & grid_unstructured_cell, za_depth_below_sea_half, &
      & t_cf_var('temperature_deriv_vert','fixme','temperature_deriv_vert', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_diag)

    CALL add_var(ocean_default_list, 'salinity_deriv_vert', ocean_state_aux%salinity_deriv_vert, &
      & grid_unstructured_cell, za_depth_below_sea_half, &
      & t_cf_var('temperature_deriv_vert','fixme','temperature_deriv_vert', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_diag)

!     CALL add_var(ocean_default_list,'DerivTracer_vert',ocean_state_aux%DerivTracer_vert_center,&
!        & grid_unstructured_cell,&
!        & za_depth_below_sea, t_cf_var('DerivTracer_vert','','', datatype_flt),&
!        & dflt_g2_decl_cell,&
!        & lopenacc=.TRUE., ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_aux,loutput=.true.)

    CALL add_var(ocean_default_list,'tracergrad_horz',ocean_state_aux%tracer_grad_horz, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('tracergrad_horz','fixme','tracergrad_horz', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_aux)

    CALL add_var(ocean_default_list,'temperaturegrad_horz',ocean_state_aux%temperature_grad_horz, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('temperaturegrad_horz','fixme','temperaturegrad_horz', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_aux)

    CALL add_var(ocean_default_list,'salinitygrad_horz',ocean_state_aux%salinity_grad_horz, grid_unstructured_edge,&
      & za_depth_below_sea, t_cf_var('salinitygrad_horz','fixme','salinitygrad_horz', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_aux)


     CALL add_var(ocean_default_list,'slopes_squared',ocean_state_aux%slopes_squared,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('slopes_squared','fixme','slopes_squared', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_aux)

     CALL add_var(ocean_default_list,'slopes_drdz',ocean_state_aux%slopes_drdz,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('slopes_drdz','fixme','slopes_drdz', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_aux)

     CALL add_var(ocean_default_list,'slopes_drdx',ocean_state_aux%slopes_drdx,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('slopes_drdx','fixme','slopes_drdx', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_aux)

     CALL add_var(ocean_default_list,'taper_function_1',ocean_state_aux%taper_function_1,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('taper_function_1','1','taper_function_1', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_aux)

     CALL add_var(ocean_default_list,'taper_function_2',ocean_state_aux%taper_function_2,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('taper_function_2','1','taper_function_2', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_aux)

     CALL add_var(ocean_default_list,'diagnose_Redi_flux_vert',ocean_state_aux%diagnose_Redi_flux_vert,&
        & grid_unstructured_cell,&
        & za_depth_below_sea, t_cf_var('diagnose_Redi_flux_vert','fixme','diagnose_Redi_flux_vert', datatype_flt),&
        & dflt_g2_decl_cell,&
        & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_aux)

    ENDIF !(GMRedi)

     ! set all values - incl. last block - of cartesian coordinates to zero (NAG compiler)
   CALL add_var(ocean_default_list,'transport_h_old',ocean_transport_state%h_old, grid_unstructured_cell,&
      & za_surface, t_cf_var('transport_h_old','fixme','transport_h_old', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_aux, lopenacc=.TRUE.)
   

   CALL add_var(ocean_default_list,'transport_h_new',ocean_transport_state%h_new, grid_unstructured_cell,&
      & za_surface, t_cf_var('transport_h_new','fixme','transport_h_new', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_aux, lopenacc=.TRUE.)
   

!      ocean_state_aux%slopes_squared=0.0_wp
  END SUBROUTINE construct_hydro_ocean_aux

  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of auxilliary hydrostatic ocean state.
  !
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_hydro_ocean_aux(ocean_state_aux)

    TYPE(t_hydro_ocean_aux), INTENT(inout)      :: ocean_state_aux

    ! local variables

    INTEGER :: ist

    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_state:destruct_hydro_ocean_aux'

    !$ACC EXIT DATA DELETE(ocean_state_aux%bc_top_veloc_cc)

    DEALLOCATE(ocean_state_aux%bc_top_veloc_cc, stat=ist)
    IF (ist/=success) THEN
      CALL finish(routine,'deallocation of top boundary cond cc failed')
    END IF

  END SUBROUTINE destruct_hydro_ocean_aux
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
!<Optimize:inUse>
  SUBROUTINE destruct_patch_3d(patch_3d)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)    :: patch_3d

    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%del_zlev_m)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%inv_del_zlev_m)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%zlev_m)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%zlev_i)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%del_zlev_i)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%ocean_area)
    DEALLOCATE(patch_3d%p_patch_1d(n_dom)%ocean_volume)
    DEALLOCATE(patch_3d%p_patch_1d)

  END SUBROUTINE destruct_patch_3d

  !-------------------------------------------------------------------------
  !>
  !! Allocation of basic 3-dimensional patch structure. This sbr assumes that
  !! the 2D horizontal patch components is already initialized.
  !
  !
  !!
!<Optimize:inUse>
  SUBROUTINE construct_patch_3d(patch_3d)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)    :: patch_3d

    ! local variables
    INTEGER :: ist
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v, n_zlvp, n_zlvm!, ie
    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_state:construct_patch_3D'
    INTEGER :: datatype_flt
    LOGICAL, DIMENSION(MAX_GROUPS) :: groups_oce_geometry
    TYPE(t_grib2_var) :: dflt_g2_decl_lonlat, dflt_g2_decl_cell, &
         dflt_g2_decl_edge, dflt_g2_decl_vrtx

    dflt_g2_decl_cell = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell)
    dflt_g2_decl_edge = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge)
    dflt_g2_decl_vrtx = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_vertex)
    dflt_g2_decl_lonlat = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_lonlat)

    groups_oce_geometry = groups("oce_geometry")

    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)

    !-------------------------------------------------------------------------

    !CALL message(routine, 'start to construct basic hydro ocean state')

    ! determine size of arrays
    alloc_cell_blocks = patch_3d%p_patch_2d(n_dom)%alloc_cell_blocks
    nblks_e = patch_3d%p_patch_2d(n_dom)%nblks_e
    nblks_v = patch_3d%p_patch_2d(n_dom)%nblks_v
    n_zlvp = n_zlev + 1
    n_zlvm = n_zlev - 1

    ALLOCATE(patch_3d%p_patch_1d(n_dom_start:n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating p_patch_1D failed')
    ENDIF
    !$ACC ENTER DATA COPYIN(patch_3d, patch_3d%p_patch_1d)

    ! allocate and set vertical level thickness from the namelist
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%del_zlev_m(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating del_zlev_m failed')
    ENDIF
    ! allocate the inverse of the above
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%inv_del_zlev_m(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating inv_del_zlev_m failed')
    ENDIF
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%zlev_m(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating zlev_m failed')
    ENDIF
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%zlev_i(n_zlvp),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating zlev_i failed')
    ENDIF
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%del_zlev_i(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating del_zlev_i failed')
    ENDIF
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%ocean_area(n_zlev),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating ocean_area failed')
    ENDIF
    ALLOCATE(patch_3d%p_patch_1d(n_dom)%ocean_volume(n_zlvp),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating ocean_volume failed')
    ENDIF
    
    !
    !! 3-dim land-sea-mask at cells, edges and vertices
    !
    ! cells
    CALL add_var(ocean_default_list, 'lsm_c', patch_3d%lsm_c, &
      & grid_unstructured_cell, za_depth_below_sea, &
      & t_cf_var('lsm_c','1','3d lsm on cells', DATATYPE_INT8),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_geometry,&
      & isteptype=tstep_constant)
    ! edges
    CALL add_var(ocean_default_list, 'lsm_e', patch_3d%lsm_e, &
      & grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('lsm_e','1','3d lsm on edges', DATATYPE_INT8),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_geometry,&
      & isteptype=tstep_constant)
    ! surface cells
    CALL add_var(ocean_default_list, 'surface_cell_sea_land_mask', patch_3d%surface_cell_sea_land_mask , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('surface_cell_sea_land_mask', '1', 'surface_cell_sea_land_mask', DATATYPE_INT8),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)
    ! surface vertices
    CALL add_var(ocean_default_list, 'surface_edge_sea_land_mask', patch_3d%surface_edge_sea_land_mask , &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('surface_edge_sea_land_mask', '1', 'surface_edge_sea_land_mask', DATATYPE_INT8),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,nblks_e/),in_group=groups_oce_geometry,isteptype=tstep_constant)
    ! surface vertices
    CALL add_var(ocean_default_list, 'surface_vertex_sea_land_mask', patch_3d%surface_vertex_sea_land_mask , &
      & grid_unstructured_vert, za_surface, &
      & t_cf_var('surface_vertex_sea_land_mask', '1', 'surface_vertex_sea_land_mask', DATATYPE_INT8),&
      & dflt_g2_decl_vrtx,&
      & ldims=(/nproma,nblks_v/),in_group=groups_oce_geometry,isteptype=tstep_constant)
    ! deepest ocean layer in column
    CALL add_var(ocean_default_list, 'dolic_c', patch_3d%p_patch_1d(n_dom)%dolic_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('dolic_c', '1', 'dolic_c', DATATYPE_INT8),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant, lopenacc=.TRUE.)
    
    CALL add_var(ocean_default_list, 'dolic_e', patch_3d%p_patch_1d(n_dom)%dolic_e , &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('dolic_e', '1', 'dolic_e', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,nblks_e/),in_group=groups_oce_geometry,isteptype=tstep_constant, lopenacc=.TRUE.)
    
    CALL add_var(ocean_default_list, 'vertex_bottomLevel', patch_3d%p_patch_1d(n_dom)%vertex_bottomLevel , &
      & grid_unstructured_vert, za_surface, &
      & t_cf_var('vertex_bottomLevel', '1', 'vertex_bottomLevel', datatype_flt),&
      & dflt_g2_decl_vrtx,&
      & ldims=(/nproma,nblks_v/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    ! 
    ! 2-dim basins and areas
    CALL add_var(ocean_default_list, 'basin_c', patch_3d%basin_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('basin_c', '1', 'basin_c', DATATYPE_INT8),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_geometry, &
      & isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'regio_c', patch_3d%regio_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('regio_c', '1', 'regio_c', DATATYPE_INT8),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_geometry, &
      & isteptype=tstep_constant)
    ! 2-dim bottom and column thickness
    CALL add_var(ocean_default_list, 'bottom_thick_c', patch_3d%bottom_thick_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('bottom_thick_c', 'm', 'bottom_thick_c', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'bottom_thick_e', patch_3d%bottom_thick_e , &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('bottom_thick_e', 'm', 'bottom_thick_e', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,nblks_e/),in_group=groups_oce_geometry,isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'column_thick_c', patch_3d%column_thick_c , &
      & grid_unstructured_cell, za_surface, &
      & t_cf_var('column_thick_c', 'm', 'column_thick_c', datatype_flt),&
      & grib2_var(10, 4, 7, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'column_thick_e', patch_3d%column_thick_e , &
      & grid_unstructured_edge, za_surface, &
      & t_cf_var('column_thick_e', 'm', 'column_thick_e', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,nblks_e/),in_group=groups_oce_geometry,isteptype=tstep_constant)
    ! 3-dim real land-sea-mask
    ! cells
    CALL add_var(ocean_default_list, 'wet_c', patch_3d%wet_c , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('wet_c', '1', '3d lsm on cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("oce_geometry","oce_default","oce_essentials"), &
      & isteptype=tstep_constant)
    ! edges
    CALL add_var(ocean_default_list, 'wet_e', patch_3d%wet_e , grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('wet_e', '1', '3d lsm on edges', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_geometry","oce_default"), &
      & isteptype=tstep_constant)
    ! 3-dim real land-sea-mask with zero on halos
    ! cells
    CALL add_var(ocean_default_list, 'wet_halo_zero_c', patch_3d%wet_halo_zero_c , grid_unstructured_cell,&
      & za_depth_below_sea, &
      & t_cf_var('wet_c_halo_zero', '1', '3d lsm with halo zero on cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)
    ! edges
    CALL add_var(ocean_default_list, 'wet_halo_zero_e', patch_3d%wet_halo_zero_e , grid_unstructured_edge,&
      & za_depth_below_sea, &
      & t_cf_var('wet_e_halo_zero', '1', '3d lsm with halo zero on edges', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_geometry,isteptype=tstep_constant)

    patch_3d%p_patch_1d(n_dom)%del_zlev_m = 0._wp
    patch_3d%p_patch_1d(n_dom)%del_zlev_i = 0._wp
    patch_3d%p_patch_1d(n_dom)%zlev_m     = 0._wp
    patch_3d%p_patch_1d(n_dom)%zlev_i     = 0._wp

    patch_3d%p_patch_1d(n_dom)%ocean_area(:)   = 0._wp
    patch_3d%p_patch_1d(n_dom)%ocean_volume(:) = 0._wp

    CALL add_var(ocean_default_list, 'prism_thick_c', patch_3d%p_patch_1d(1)%prism_thick_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('cons thick','m','prism thickness at cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'invConstantPrismThickness', patch_3d%p_patch_1d(1)%invConstantPrismThickness, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('inv cons thick','m','inverse prism thickness at cells', DATATYPE_FLT32),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)
    CALL add_var(ocean_default_list, 'prism_volume', patch_3d%p_patch_1d(1)%prism_volume, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('cons thick','m','prism volume (cells)', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant, lopenacc=.TRUE.)
    
    CALL add_var(ocean_default_list, 'prism_thick_e', patch_3d%p_patch_1d(n_dom)%prism_thick_e, &
      & grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('cons thick','m','prism thickness at edges', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'prism_thick_flat_sfc_c', patch_3d%p_patch_1d(n_dom)%prism_thick_flat_sfc_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('prism_thick_flat_sfc_c','m','time independent depth at cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'prism_thick_flat_sfc_e', patch_3d%p_patch_1d(n_dom)%prism_thick_flat_sfc_e, &
      & grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('prism_thick_flat_sfc_e','m','time independent depth at edges', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'inverse prism_thick_c', patch_3d%p_patch_1d(n_dom)%inv_prism_thick_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea, &
      & t_cf_var('inverse prism_thick_c','m','time dependent depth at cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'prism_center_dist_c', patch_3d%p_patch_1d(n_dom)%prism_center_dist_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea_half, &
      & t_cf_var('prism_center_dist_c','m','time dependent distance between prism centers', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'constantPrismCenters_Zdistance', &
      & patch_3d%p_patch_1d(n_dom)%constantPrismCenters_Zdistance, &
      & grid_unstructured_cell, &
      & za_depth_below_sea_half, &
      & t_cf_var('constantPrismCenters_Zdistance','m','constant distance between prism centers', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'constantPrismCenters_invZdistance', &
      & patch_3d%p_patch_1d(n_dom)%constantPrismCenters_invZdistance, &
      & grid_unstructured_cell, &
      & za_depth_below_sea_half, &
      & t_cf_var('constantPrismCenters_invZdistance','m','inverse constant distance between prism centers', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'inv_prism_thick_e', patch_3d%p_patch_1d(n_dom)%inv_prism_thick_e, &
      & grid_unstructured_edge, &
      & za_depth_below_sea, &
      & t_cf_var('inv_prism_thick_e','m','time dependent inverse thickeness at edges', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'inv_prism_center_dist_c', &
      & patch_3d%p_patch_1d(n_dom)%inv_prism_center_dist_c, &
      & grid_unstructured_cell, &
      & za_depth_below_sea_half, &
      & t_cf_var('inv_prism_center_dist_c','1/m','inverse of dist between prism centers at cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'inv_prism_center_dist_e', &
      & patch_3d%p_patch_1d(n_dom)%inv_prism_center_dist_e, &
      & grid_unstructured_edge, &
      & za_depth_below_sea_half, &
      & t_cf_var('inv_prism_center_dist_e','1/m','inverse of dist between prism centers at edges', datatype_flt),&
      & dflt_g2_decl_edge,&
      & ldims=(/nproma,n_zlev+1,nblks_e/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'depth_CellMiddle',   &
      & patch_3d%p_patch_1d(n_dom)%depth_CellMiddle,       &
      & grid_unstructured_cell,                             &
      & za_depth_below_sea,                                 &
      & t_cf_var('depth_CellMiddle','m','depth at the middle of the cells', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !
    CALL add_var(ocean_default_list, 'depth_CellInterface',   &
      & patch_3d%p_patch_1d(n_dom)%depth_CellInterface,       &
      & grid_unstructured_cell,                                &
      & za_depth_below_sea_half,                               &
      & t_cf_var('depth_CellInterface','m','depth at cell interfaces', datatype_flt),&
      & dflt_g2_decl_cell,&
      & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups_oce_geometry,isteptype=tstep_constant)!, lopenacc=.TRUE.)
    !

  END SUBROUTINE construct_patch_3d

  !------------------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE setup_tracer_info(oce_config)
      TYPE(t_oce_config) :: oce_config
    oce_config%tracer_shortnames(1) = 'to'
    oce_config%tracer_stdnames(1)   = 'sea_water_potential_temperature'
    oce_config%tracer_longnames(1)  = 'sea water potential temperature'
    oce_config%tracer_units(1)      = 'C'
    ! discipline=10, parameterCategory=4, parameterNumber=18 encoded in one integer
    oce_config%tracer_codes(1)      = ISHFT(10,16)+ISHFT(4,8)+18

    oce_config%tracer_shortnames(2) = 'so'
    oce_config%tracer_stdnames(2)   = 'sea_water_salinity'
    oce_config%tracer_longnames(2)  = 'sea water salinity'
    oce_config%tracer_units(2)      = 'psu'
    ! discipline=10, parameterCategory=4, parameterNumber=21 encoded in one integer
    oce_config%tracer_codes(2)      = ISHFT(10,16)+ISHFT(4,8)+21
  END SUBROUTINE setup_tracer_info

  SUBROUTINE transfer_ocean_state( patch_3d, operators_coefficients )

    TYPE(t_patch_3d), POINTER, INTENT(in) :: patch_3d
    TYPE(t_operator_coeff),   INTENT(inout), TARGET :: operators_coefficients

    INTEGER :: blockNo

    !$ACC ENTER DATA COPYIN(patch_3d%p_patch_2d) &
    !$ACC   COPYIN(patch_3d%lsm_e, patch_3d%surface_cell_sea_land_mask) &
    !$ACC   COPYIN(patch_3d%basin_c, patch_3d%wet_c, patch_3d%wet_halo_zero_c, patch_3d%lsm_c) &
    !$ACC   COPYIN(patch_3d%column_thick_e, patch_3d%column_thick_c) &
    !$ACC   COPYIN(patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%zlev_i) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%zlev_m) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%del_zlev_i) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%del_zlev_m) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%vertex_bottomlevel) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%constantprismcenters_zdistance) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%constantprismcenters_invzdistance) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%depth_cellmiddle) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%depth_cellinterface) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%prism_center_dist_c) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%prism_thick_c) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%inv_prism_thick_c) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%inv_prism_center_dist_c) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%inv_prism_center_dist_e) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%prism_thick_e) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%inv_prism_thick_e) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%alloc_cell_blocks) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%nblks_v, patch_3d%p_patch_2D(1)%nblks_e) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%cells, patch_3d%p_patch_2d(1)%edges) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%edges%cell_idx, patch_3d%p_patch_2d(1)%edges%cell_blk) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%edges%in_domain) &
    !$ACC   COPYIN(patch_3D%p_patch_2d(1)%edges%area_edge) &
    !$ACC   COPYIN(patch_3D%p_patch_2D(1)%edges%vertex_idx) &
    !$ACC   COPYIN(patch_3D%p_patch_2D(1)%edges%vertex_blk) &
    !$ACC   COPYIN(patch_3D%p_patch_2D(1)%edges%primal_cart_normal) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%edges%inv_dual_edge_length) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%edges%primal_edge_length) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%edges%tangent_orientation) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%edges%inv_primal_edge_length) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%cells%edge_idx, patch_3d%p_patch_2d(1)%cells%edge_blk) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%cells%neighbor_idx, patch_3d%p_patch_2d(1)%cells%neighbor_blk) &
    !$ACC   COPYIN(patch_3D%p_patch_2D(1)%cells%center) &
    !$ACC   COPYIN(patch_3D%p_patch_2D(1)%cells%owned, patch_3D%p_patch_2D(1)%cells%owned%vertical_levels) &
    !$ACC   COPYIN(patch_3D%p_patch_2D(1)%cells%all) &
    !$ACC   COPYIN(patch_3D%p_patch_2D(1)%cells%all%vertical_levels, patch_3D%p_patch_2D(1)%cells%area) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%cells%in_domain) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%cells%in_domain%vertical_levels) &
    !$ACC   COPYIN(patch_3D%p_patch_2D(1)%verts, patch_3D%p_patch_2D(1)%verts%cell_idx) &
    !$ACC   COPYIN(patch_3D%p_patch_2D(1)%verts%cell_blk, patch_3D%p_patch_2D(1)%verts%num_edges) &
    !$ACC   COPYIN(patch_3D%p_patch_2D(1)%verts%edge_idx, patch_3D%p_patch_2D(1)%verts%edge_blk) &
    !$ACC   COPYIN(patch_3d%p_patch_2d(1)%verts%f_v) &
    !$ACC   COPYIN(patch_3D%p_patch_2D(1)%verts%all) &
    !$ACC   COPYIN(nold, nnew)
    
    !$ACC ENTER DATA COPYIN(operators_coefficients, operators_coefficients%verticaladvectionppmcoeffs)
    DO blockNo = patch_3D%p_patch_2D(1)%cells%ALL%start_block, patch_3D%p_patch_2D(1)%cells%ALL%end_block
      !$ACC ENTER DATA COPYIN(operators_coefficients%verticaladvectionppmcoeffs(blockNo)%cellheightratio_this_tobelow) &
      !$ACC   COPYIN(operators_coefficients%verticaladvectionppmcoeffs(blockNo)%cellheightratio_this_tothisbelow) &
      !$ACC   COPYIN(operators_coefficients%verticaladvectionppmcoeffs(blockNo)%cellheight_2xbelow_x_ratiothis_tothisbelow) &
      !$ACC   COPYIN(operators_coefficients%verticaladvectionppmcoeffs(blockNo)%cellheightratio_this_tothisabovebelow) &
      !$ACC   COPYIN(operators_coefficients%verticaladvectionppmcoeffs(blockNo)%cellheightratio_2xaboveplusthis_tothisbelow) &
      !$ACC   COPYIN(operators_coefficients%verticaladvectionppmcoeffs(blockNo)%cellheightratio_2xbelowplusthis_tothisabove) &
      !$ACC   COPYIN(operators_coefficients%verticaladvectionppmcoeffs(blockNo)%cellheightratio_thisabove_to2xthisplusbelow) &
      !$ACC   COPYIN(operators_coefficients%verticaladvectionppmcoeffs(blockNo)%cellheightratio_thisbelow_to2xthisplusabove) &
      !$ACC   COPYIN(operators_coefficients%verticaladvectionppmcoeffs(blockNo)%cellheight_inv_thisabovebelow2below)
    END DO

    !$ACC ENTER DATA &
    !$ACC   COPYIN(operators_coefficients, operators_coefficients%vertex_bnd_edge_idx, operators_coefficients%vertex_bnd_edge_blk) &
    !$ACC   COPYIN(operators_coefficients%edge2cell_coeff_cc_t) &
    !$ACC   COPYIN(operators_coefficients%edge2vert_coeff_cc) &
    !$ACC   COPYIN(operators_coefficients%edge2cell_coeff_cc) &
    !$ACC   COPYIN(operators_coefficients%edge2vert_coeff_cc_t) &
    !$ACC   COPYIN(operators_coefficients%div_coeff) &
    !$ACC   COPYIN(operators_coefficients%edge2edge_viacell_coeff_all) &
    !$ACC   COPYIN(operators_coefficients%edge2edge_viacell_coeff_integrated) &
    !$ACC   COPYIN(operators_coefficients%edge2edge_viacell_coeff_top) &
    !$ACC   COPYIN(operators_coefficients%lhs_all) &
    !$ACC   COPYIN(operators_coefficients%grad_coeff) &
    !$ACC   COPYIN(operators_coefficients%lhs_celltocell_block) &
    !$ACC   COPYIN(operators_coefficients%lhs_celltocell_index) &
    !$ACC   COPYIN(operators_coefficients%fixed_vol_norm) &
    !$ACC   COPYIN(operators_coefficients%edge2edge_viacell_coeff) &
    !$ACC   COPYIN(operators_coefficients%rot_coeff) &
    !$ACC   COPYIN(operators_coefficients%edge2edge_viavert_coeff) &
    !$ACC   COPYIN(operators_coefficients%bnd_edges_per_vertex, operators_coefficients%boundaryedge_coefficient_index) &
    !$ACC   COPYIN(operators_coefficients%edges_seaboundarylevel)

    CALL gpu_update_var_list('ocean_restart_list' , .true., lacc=.true.)
    CALL gpu_update_var_list('ocean_default_list' , .true., lacc=.true.)
    CALL gpu_update_var_list('ocean_params_list' , .true., lacc=.true.)
    CALL gpu_update_var_list('ext_data_oce_D01' , .true., lacc=.true.)
    CALL gpu_update_var_list('ext_data_bgc_D01' , .true., lacc=.true.)

  END SUBROUTINE transfer_ocean_state

END MODULE mo_ocean_state
