! Contains basic diagnostics for ICON ocean model.
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
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_diagnostics
  USE mo_master_control,     ONLY: get_my_process_name
  USE mo_kind,               ONLY: wp, dp, i8
#ifdef _OPENMP
  USE omp_lib
#endif
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range, t_subset_indexed
  USE mo_grid_tools,         ONLY: get_oriented_edges_from_global_vertices, check_global_indexes
  USE mo_mpi,                ONLY: my_process_is_stdio, p_field_sum, get_my_mpi_work_id, &
    & p_comm_work_test, p_comm_work, p_io, p_bcast, my_process_is_mpi_workroot, p_sum, &
    & my_process_is_mpi_parallel
  USE mo_sync,               ONLY: global_sum_array, disable_sync_checks, enable_sync_checks, &
    &                              sync_c, sync_e, sync_patch_array
  USE mo_math_types,         ONLY: t_cartesian_coordinates
  USE mo_math_utilities,     ONLY: cvec2gvec
  USE mo_advection_utils,    ONLY: laxfr_upflux
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_dbg_nml,            ONLY: idbg_val
  USE mo_math_constants,     ONLY: rad2deg, dbl_eps
  USE mo_impl_constants,     ONLY: sea_boundary,sea, &
    & min_rlcell, min_rledge, min_rlcell, &
    & min_dolic, nlat_moc
  USE mo_timer,              ONLY: timer_calc_moc, timer_start, timer_stop
  USE mo_cdi_constants,      ONLY: GRID_EDGE, GRID_CELL, GRID_UNSTRUCTURED_EDGE, &
    & GRID_UNSTRUCTURED_CELL
  USE mo_ocean_nml,          ONLY: n_zlev, no_tracer, &
    & gibraltar, &
    & denmark_strait, &
    & drake_passage, &
    & florida_strait, &
    & indonesian_throughflow,&
    & scotland_iceland, &
    & mozambique, &
    & framStrait, &
    & beringStrait, &
    & barentsOpening, &
    & agulhas, &
    & agulhas_long, &
    & agulhas_longer, &
    & ab_const, ab_beta, ab_gam, iswm_oce, discretization_scheme, &
    & iforc_oce, No_Forcing, i_sea_ice, diagnostics_level, &
    & diagnose_for_horizontalVelocity, OceanReferenceDensity, &
    & eddydiag, &
    & vert_cor_type, check_total_volume
  USE mo_sea_ice_nml,        ONLY: kice, sice
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_parallel_config,    ONLY: nproma, p_test_run
  USE mo_run_config,         ONLY: dtime, nsteps
  USE mo_physical_constants, ONLY: grav, rhos, rhoi, rho_ref, clw, alf
  USE mo_model_domain,       ONLY: t_patch, t_patch_3d,t_patch_vert, t_grid_edges
  USE mo_ocean_types,        ONLY: t_hydro_ocean_state, t_hydro_ocean_diag
  USE mo_ocean_diagnostics_types,  ONLY: t_ocean_regions, t_ocean_region_volumes, &
    &  t_ocean_region_areas, t_ocean_monitor
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_exception,          ONLY: message, finish, message_text, warning
  USE mo_sea_ice_types,      ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_ocean_surface_types,ONLY: t_ocean_surface
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_scalar_product,     ONLY: map_edges2cell_3d
  USE mo_io_units,           ONLY: find_next_free_unit
  USE mo_util_file,          ONLY: util_symlink, util_rename, util_islink, util_unlink
  USE mo_statistics,         ONLY: subset_sum, levels_horizontal_mean, total_mean, gather_sums, &
    & verticallyIntegrated_field
  USE mo_var_list,           ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var_list_register,  ONLY: vlr_add, vlr_del
  USE mo_var_groups,         ONLY: groups
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, GRID_UNSTRUCTURED
  USE mo_cdi_constants,      ONLY: GRID_EDGE, GRID_CELL, GRID_UNSTRUCTURED_EDGE, &
    &                              GRID_UNSTRUCTURED_CELL
  USE mo_zaxis_type,         ONLY: ZA_DEPTH_BELOW_SEA
  USE mo_io_config,          ONLY: lnetcdf_flt64_output
  USE mo_name_list_output_init, ONLY: isRegistered

  USE mtime,                 ONLY: datetime, MAX_DATETIME_STR_LEN, datetimeToPosixString
  USE mo_ocean_check_total_content , ONLY : calc_total_salt_content, &
    & calc_total_salt_content_zstar
  USE mo_fortran_tools,      ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=12)           :: str_module    = 'oceDiag     '  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

  INTEGER :: moc_unit  = -1 ! file handle for the global timeseries output

  !
  ! PUBLIC INTERFACE
  !
  PUBLIC :: calc_fast_oce_diagnostics
  PUBLIC :: construct_oce_diagnostics
  PUBLIC :: destruct_oce_diagnostics
  PUBLIC :: calc_moc
  PUBLIC :: calc_psi
  PUBLIC :: diag_heat_salt_tendency

  INTERFACE calc_moc
! 2023-11 dzo-DKRZ: Temporary until the difference in atlantic_moc and global_moc between both functions is solved
#if defined(__LVECTOR__) && !defined(__LVEC_BITID__)
    MODULE PROCEDURE calc_moc_hfl_internal_vector
#else
    MODULE PROCEDURE calc_moc_hfl_internal_scalar
#endif
  END INTERFACE

  TYPE t_oce_section
    TYPE(t_subset_indexed) :: subset
    REAL(wp), POINTER :: orientation(:)
  END TYPE t_oce_section

  INTEGER, PARAMETER  :: oce_section_count = 13
  PRIVATE             :: oce_section_count
  TYPE(t_oce_section) :: oce_sections(oce_section_count)
  PRIVATE             :: oce_sections

  TYPE(t_ocean_region_volumes),SAVE :: ocean_region_volumes
  PRIVATE                           :: ocean_region_volumes
  TYPE(t_ocean_region_areas),SAVE   :: ocean_region_areas
  PRIVATE                           :: ocean_region_areas


  TYPE(t_var_list_ptr) :: horizontal_velocity_diagnostics
  ! addtional diagnostics
  REAL(wp), POINTER :: veloc_adv_horz_u(:,:,:),  veloc_adv_horz_v(:,:,:), &
    & laplacian_horz_u(:,:,:), laplacian_horz_v(:,:,:), vn_u(:,:,:), vn_v(:,:,:), &
    & mass_flx_e_u(:,:,:), mass_flx_e_v(:,:,:), pressure_grad_u(:,:,:), pressure_grad_v(:,:,:), &
    & potential_vort_e(:,:,:), potential_vort_c(:,:,:)

  CHARACTER(*), PARAMETER :: module_name="mo_ocean_statistics"

CONTAINS

  !-------------------------------------------------------------------------
  !  The constructor of the types related to ocean diagnostics
  !>
  !!
!<Optimize:inUse>
  SUBROUTINE construct_oce_diagnostics( patch_3D, ocean_state )
    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET      :: ocean_state

    !local variable
    INTEGER :: i,ist
    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_diagnostics:construct_oce_diagnostics'
    !-----------------------------------------------------------------------
    INTEGER  :: nblks_e,blockNo,jc,jk, region_index,start_index,end_index
    REAL(wp) :: surface_area, surface_height, prism_vol, prism_area, column_volume

    TYPE(t_patch), POINTER        :: patch_2d
    TYPE(t_subset_range), POINTER :: owned_cells
    INTEGER, POINTER              :: regions(:,:)
    TYPE(t_ocean_regions)         :: ocean_regions
    INTEGER                       :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    CALL message (routine, 'start')
    !-----------------------------------------------------------------------
    patch_2d => patch_3D%p_patch_2d(1)
    regions => patch_3D%regio_c
    !-----------------------------------------------------------------------
    owned_cells => patch_2d%cells%owned
    nblks_e = patch_2d%nblks_e
    !-----------------------------------------------------------------------
    CALL vlr_add(horizontal_velocity_diagnostics, 'horizontal_velocity_diagnostics', &
      & patch_id=patch_2d%id, lrestart=.FALSE., loutput=.TRUE.,                           &
      & model_type=TRIM(get_my_process_name()))
    !-----------------------------------------------------------------------
    IF (diagnose_for_horizontalVelocity) THEN
      CALL add_var(horizontal_velocity_diagnostics, 'veloc_adv_horz_u', veloc_adv_horz_u, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('veloc_adv_horz_u','m/s','velocity advection zonal', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'veloc_adv_horz_v', veloc_adv_horz_v, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('veloc_adv_horz_v','m/s','velocity advection meridional', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'laplacian_horz_u', laplacian_horz_u, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('laplacian_horz_u','m/s','velocity laplacian zonal', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'laplacian_horz_v', laplacian_horz_v, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('laplacian_horz_v','m/s','velocity laplacian meridional', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'vn_u', vn_u, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('vn_u','m/s','edge velocity zonal', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'vn_v', vn_v, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('vn_v','m/s','edge velocity meridional', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'mass_flx_e_u', mass_flx_e_u, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('mass_flx_e_u','m*m/s','mass flux zonal', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'mass_flx_e_v', mass_flx_e_v, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('mass_flx_e_v','m*m/s','mass flux meridional', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'pressure_grad_u', pressure_grad_u, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('pressure_grad_u','N','pressure gradient zonal', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'pressure_grad_v', pressure_grad_v, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('pressure_grad_v','N','pressure gradient meridional', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'potential_vort_e', potential_vort_e, &
        & grid_unstructured_edge, za_depth_below_sea, &
        & t_cf_var('vn_v','1/s','potential vorticity at edges', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_edge),&
        & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

      CALL add_var(horizontal_velocity_diagnostics, 'potential_vort_c', potential_vort_c, &
        & grid_unstructured_cell, za_depth_below_sea, &
        & t_cf_var('vn_v','1/s','potential vorticity at cells', &
        & datatype_flt),&
        & grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, grid_cell),&
        & ldims=(/nproma,n_zlev,patch_2d%alloc_cell_blocks/),in_group=groups("oce_diag"),lrestart_cont=.FALSE.)

    ENDIF

    ! compute subsets for given sections path allong edges
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(1)%subset,      &
      & orientation = oce_sections(1)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array = gibraltar,            &
      & subset_name = 'gibraltar')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(2)%subset,      &
      & orientation = oce_sections(2)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =denmark_strait,        &
      & subset_name = 'denmark_strait')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(3)%subset,      &
      & orientation = oce_sections(3)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =drake_passage,         &
      & subset_name = 'drake_passage')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(4)%subset,      &
      & orientation = oce_sections(4)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =indonesian_throughflow,&
      & subset_name = 'indonesian_throughflow')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(5)%subset,      &
      & orientation = oce_sections(5)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =scotland_iceland,      &
      & subset_name = 'scotland_iceland')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(6)%subset,      &
      & orientation = oce_sections(6)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =mozambique,      &
      & subset_name = 'mozambique')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(7)%subset,      &
      & orientation = oce_sections(7)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =framStrait,      &
      & subset_name = 'framStrait')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(8)%subset,      &
      & orientation = oce_sections(8)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =beringStrait,      &
      & subset_name = 'beringStrait')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(9)%subset,      &
      & orientation = oce_sections(9)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =barentsOpening,      &
      & subset_name = 'barentsOpening')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(10)%subset,      &
      & orientation = oce_sections(10)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =agulhas,      &
      & subset_name = 'agulhas')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(11)%subset,      &
      & orientation = oce_sections(11)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =agulhas_long,      &
      & subset_name = 'agulhas_long')
    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(12)%subset,      &
      & orientation = oce_sections(12)%orientation, &
      & patch_3d = patch_3D,                     &
      & global_vertex_array =agulhas_longer,      &
      & subset_name = 'agulhas_longer')

    CALL get_oriented_edges_from_global_vertices(    &
      & edge_subset = oce_sections(13)%subset,       &
      & orientation = oce_sections(13)%orientation, &
      & patch_3d = patch_3D,                         &
      & global_vertex_array =florida_strait,         &
      & subset_name = 'florida_strait')
!     CALL finish("","")

    surface_area   = 0.0_wp
    surface_height = 0.0_wp
    prism_vol      = 0.0_wp
    prism_area     = 0.0_wp
    ocean_region_areas%total = 0.0_wp
    ! compute regional ocean volumes
    DO blockNo = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index

        ! area
        prism_area               = patch_2d%cells%area(jc,blockNo)
        ocean_region_areas%total = ocean_region_areas%total + prism_area

        ! volume
        CALL compute_vertical_volume(blockNo,jc, &
          & prism_area, &
          & ocean_state%p_prog(nnew(1))%h(jc,blockNo), &
          & patch_3D%p_patch_1d(1)%prism_thick_c(jc,:,blockNo), &
          & patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo), &
          & column_volume)
        ocean_region_volumes%total = ocean_region_volumes%total + column_volume

        region_index = regions(jc,blockNo)
        IF (ocean_regions%greenland_iceland_norwegian_sea == region_index) THEN
          ocean_region_volumes%greenland_iceland_norwegian_sea = &
            & ocean_region_volumes%greenland_iceland_norwegian_sea + column_volume
          ocean_region_areas%greenland_iceland_norwegian_sea   = &
            & ocean_region_areas%greenland_iceland_norwegian_sea + prism_area

        ELSEIF (ocean_regions%arctic_ocean == region_index) THEN
          ocean_region_volumes%arctic_ocean                    = ocean_region_volumes%arctic_ocean + column_volume
          ocean_region_areas%arctic_ocean                      = ocean_region_areas%arctic_ocean + prism_area

        ELSEIF (ocean_regions%labrador_sea == region_index) THEN
          ocean_region_volumes%labrador_sea                    = ocean_region_volumes%labrador_sea + column_volume
          ocean_region_areas%labrador_sea                      = ocean_region_areas%labrador_sea + prism_area

        ELSEIF (ocean_regions%north_atlantic == region_index) THEN
          ocean_region_volumes%north_atlantic                  = ocean_region_volumes%north_atlantic + column_volume
          ocean_region_areas%north_atlantic                    = ocean_region_areas%north_atlantic + prism_area

        ELSEIF (ocean_regions%tropical_atlantic == region_index) THEN
          ocean_region_volumes%tropical_atlantic               = ocean_region_volumes%tropical_atlantic + column_volume
          ocean_region_areas%tropical_atlantic                 = ocean_region_areas%tropical_atlantic + prism_area

        ELSEIF (ocean_regions%southern_ocean == region_index) THEN
          ocean_region_volumes%southern_ocean                  = ocean_region_volumes%southern_ocean + column_volume
          ocean_region_areas%southern_ocean                    = ocean_region_areas%southern_ocean + prism_area

        ELSEIF (ocean_regions%indian_ocean == region_index) THEN
          ocean_region_volumes%indian_ocean                    = ocean_region_volumes%indian_ocean + column_volume
          ocean_region_areas%indian_ocean                      = ocean_region_areas%indian_ocean + prism_area

        ELSEIF (ocean_regions%tropical_pacific == region_index) THEN
          ocean_region_volumes%tropical_pacific                = ocean_region_volumes%tropical_pacific + column_volume
          ocean_region_areas%tropical_pacific                  = ocean_region_areas%tropical_pacific + prism_area

        ELSEIF (ocean_regions%north_pacific == region_index) THEN
          ocean_region_volumes%north_pacific                   = ocean_region_volumes%north_pacific + column_volume
          ocean_region_areas%north_pacific                     = ocean_region_areas%north_pacific + prism_area

        ELSEIF (ocean_regions%caribbean == region_index) THEN
          ocean_region_volumes%caribbean                       = ocean_region_volumes%caribbean + column_volume
          ocean_region_areas%caribbean                         = ocean_region_areas%caribbean + prism_area
        END IF

      END DO
    END DO
    ! compute global values

    CALL disable_sync_checks()
    ocean_region_volumes%total = global_sum_array(ocean_region_volumes%total)
    ocean_region_areas%total   = global_sum_array(ocean_region_areas%total)
    CALL enable_sync_checks()

    CALL message (routine, 'end')
  END SUBROUTINE construct_oce_diagnostics
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !  !The destructor of the types related to ocean diagnostics
  !>
  !!
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_oce_diagnostics()
    !
    !
    !local variables
    INTEGER :: i,iret

    CHARACTER(LEN=*), PARAMETER :: &
      & routine = 'mo_ocean_diagnostics:destruct_oce_diagnostics'
    !-----------------------------------------------------------------------
    CALL vlr_del(horizontal_velocity_diagnostics)

    IF (diagnostics_level <= 0) RETURN

    CALL message (routine, 'end')
  END SUBROUTINE destruct_oce_diagnostics
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE compute_vertical_volume(blockNo,jc,prism_area,surface_height,thicknesses,max_vertical_level,volume)
    INTEGER,  INTENT(in)  :: blockNo,jc,max_vertical_level
    REAL(wp), INTENT(in)  :: prism_area, surface_height, thicknesses(:)
    REAL(wp), INTENT(inout) :: volume

    INTEGER :: jk
    REAL(wp) :: surface_height_,prism_vol_

    volume  = 0.0_wp
    DO jk = 1,max_vertical_level
      !local volume
      surface_height_ = MERGE(surface_height,0.0_wp, 1 == jk)
      prism_vol_      = prism_area * (thicknesses(jk) + surface_height_)
      !Fluid volume wrt lsm
      volume          = volume + prism_vol_
    END DO
  END SUBROUTINE compute_vertical_volume
  !-------------------------------------------------------------------------


!<Optimize:inUse>
  SUBROUTINE calc_fast_oce_diagnostics(patch_2d, patch_3d, ocean_state, dolic, prism_thickness, depths, &
          &  p_diag, sea_surface_height, normal_veloc, tracers, p_atm_f, p_oce_sfc, ice, lacc)
    TYPE(t_patch ),TARGET :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(inout)     :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout)    :: ocean_state
    INTEGER,  POINTER                           :: dolic(:,:)
    REAL(wp), POINTER                           :: prism_thickness(:,:,:)
    REAL(wp), INTENT(in)                        :: depths(:)
    TYPE(t_hydro_ocean_diag), TARGET            :: p_diag
    REAL(wp), POINTER                           :: sea_surface_height(:,:)
    REAL(wp), POINTER                           :: normal_veloc(:,:,:)
    REAL(wp), POINTER                           :: tracers(:,:,:,:)
     TYPE(t_atmos_fluxes ),    INTENT(IN)        :: p_atm_f
    TYPE(t_ocean_surface), INTENT(IN)           :: p_oce_sfc
    TYPE(t_sea_ice),          INTENT(inout)     :: ice
    LOGICAL, INTENT(IN), OPTIONAL                     :: lacc

    REAL(wp) :: w(nproma, n_zlev + 1, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL :: lzacc

    !Local variables
    INTEGER :: start_cell_index, end_cell_index,i
    INTEGER :: jk,jc,blockNo!,je
    REAL(wp):: ssh_global_mean,sst_global,sss_global,total_runoff_flux,total_heat_flux, &
      &        total_precipitation_flux,total_evaporation_flux, atmos_snowfall_flux, &
      &        ice_volume_nh, ice_volume_sh, ice_extent_nh, ice_extent_sh, &
      &        global_mean_potEnergy, global_mean_kinEnergy, global_mean_totalEnergy, &
      &        global_mean_potEnstrophy,global_heat_content, global_heat_content_solid, &
      &        VolumeIce_flux, TotalOcean_flux, TotalIce_flux, VolumeTotal_flux, totalsnowfall_flux
!   REAL(wp) :: sflux
      REAL(wp) ::total_salt, total_saltinseaice, total_saltinliquidwater
    CHARACTER(LEN=*), PARAMETER :: method_name='mo_ocean_diagnostics:calc_fast_oce_diagnostics'

    TYPE(t_subset_range), POINTER :: owned_cells, owned_edges
    TYPE(t_ocean_monitor),  POINTER :: monitor
    !-----------------------------------------------------------------------
    owned_cells    => patch_2d%cells%owned
    owned_edges    => patch_2d%edges%owned
    monitor        => p_diag%monitor

    CALL set_acc_host_or_device(lzacc, lacc)

    !cell loop to calculate cell based monitored fields volume, kinetic energy and tracer content
    SELECT CASE (iswm_oce)
    CASE (1) ! shallow water mode

    CASE default !3D model

      ! {{{ compute global mean values of:
      ! total_salt
      total_salt = 0.0_wp
      total_saltinseaice = 0.0_wp
      total_saltinliquidwater = 0.0_wp

      IF (isRegistered('total_salt')) THEN
        IF (vert_cor_type .EQ. 0) THEN
          call calc_total_salt_content(tracers(:,:,:,2), patch_2d, &
            & sea_surface_height(:,:),&
            & patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
            & ice, 0 , total_salt, total_saltinseaice, &
            & total_saltinliquidwater, lacc=lzacc )
        ELSE IF (vert_cor_type .EQ. 1) THEN
          call calc_total_salt_content_zstar(tracers(:,:,:,2), patch_2d, &
            & ocean_state%p_prog(nnew(1))%stretch_c(:, :), &
            & patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:),&
            & ice, total_salt, total_saltinseaice, &
            & total_saltinliquidwater, lacc=lzacc )
        END IF

        monitor%total_salt = total_salt
        monitor%total_saltinseaice = total_saltinseaice
        monitor%total_saltinliquidwater = total_saltinliquidwater

        !write(0,*)'total_salt_content =' , monitor%total_salt

      END IF

      ! {{{ compute global mean values of:
      ! sea surface height
      ssh_global_mean = 0.0_wp
      IF (isRegistered('ssh_global')) THEN
        CALL levels_horizontal_mean( sea_surface_height, &
            & patch_2d%cells%area(:,:), &
            & owned_cells, &
            & ssh_global_mean, lopenacc=lzacc)
      END IF
      monitor%ssh_global = ssh_global_mean
      IF (my_process_is_stdio() .and. check_total_volume) THEN
        WRITE(0,*) ' -- monitor%ssh_global:', monitor%ssh_global
      ENDIF


      ! sea surface temperature
      sst_global = 0.0_wp
      IF (isRegistered('sst_global')) THEN
!       CALL levels_horizontal_mean( p_oce_sfc%sst, &
        CALL levels_horizontal_mean( tracers(:,1,:,1), &
            & patch_2d%cells%area(:,:), &
            & owned_cells, &
            & sst_global, lopenacc=lzacc)
      END IF
      monitor%sst_global = sst_global

      ! sea surface salinity
      sss_global = 0.0_wp
      IF (isRegistered('sss_global')) THEN
!       CALL levels_horizontal_mean( p_oce_sfc%sss, &
        CALL levels_horizontal_mean( tracers(:,1,:,2), &
            & patch_2d%cells%area(:,:), &
            & owned_cells, &
            & sss_global, lopenacc=lzacc)
      END IF
      monitor%sss_global = sss_global

      ! total heat flux
      total_heat_flux = 0.0_wp
      IF (isRegistered('HeatFlux_Total_global')) THEN
        CALL levels_horizontal_mean( p_oce_sfc%HeatFlux_Total, &
            & patch_2d%cells%area(:,:), &
            & owned_cells, &
            & total_heat_flux, lopenacc=lzacc)
      END IF
      monitor%HeatFlux_Total = total_heat_flux

      ! total precipitation flux
      total_precipitation_flux = 0.0_wp
      IF (isRegistered('FrshFlux_Precipitation_Global')) THEN
      call levels_horizontal_mean( p_oce_sfc%FrshFlux_Precipitation, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
            & total_precipitation_flux, lopenacc=lzacc)
      END IF
      monitor%FrshFlux_Precipitation = total_precipitation_flux

      ! total evaporation
      total_evaporation_flux = 0.0_wp
      IF (isRegistered('FrshFlux_Evaporation_Global')) THEN
      call levels_horizontal_mean( p_oce_sfc%FrshFlux_Evaporation, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & total_evaporation_flux, lopenacc=lzacc)
      END IF
      monitor%FrshFlux_Evaporation = total_evaporation_flux

      ! total runoff
      total_runoff_flux = 0.0_wp
      IF (isRegistered('FrshFlux_Runoff_Global')) THEN
      call levels_horizontal_mean( p_oce_sfc%FrshFlux_Runoff, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & total_runoff_flux, lopenacc=lzacc)
      END IF
      monitor%FrshFlux_Runoff = total_runoff_flux

      ! total (atmospheric) snowfall
      atmos_snowfall_flux = 0.0_wp
      IF (isRegistered('FrshFlux_SnowFall_Global')) THEN
      call levels_horizontal_mean( p_oce_sfc%FrshFlux_Snowfall, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & atmos_snowfall_flux, lopenacc=lzacc)
      END IF
      monitor%FrshFlux_SnowFall = atmos_snowfall_flux

      ! VolumeIce
      VolumeIce_flux = 0.0_wp
      IF (isRegistered('FrshFlux_VolumeIce_Global')) THEN
      call levels_horizontal_mean( p_oce_sfc%FrshFlux_VolumeIce, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & VolumeIce_flux, lopenacc=lzacc)
      END IF
      monitor%FrshFlux_VolumeIce = VolumeIce_flux

      ! TotalOcean
      TotalOcean_flux = 0.0_wp
      IF (isRegistered('FrshFlux_TotalOcean_Global')) THEN
      call levels_horizontal_mean( p_oce_sfc%FrshFlux_TotalOcean, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & TotalOcean_flux, lopenacc=lzacc)
      END IF
      monitor%FrshFlux_TotalOcean = TotalOcean_flux

      ! TotalIce
      TotalIce_flux = 0.0_wp
      IF (isRegistered('FrshFlux_TotalIce_Global')) THEN
      call levels_horizontal_mean( p_oce_sfc%FrshFlux_TotalIce, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & TotalIce_flux, lopenacc=lzacc)
      END IF
      monitor%FrshFlux_TotalIce = TotalIce_flux

      ! VolumeTotal
      VolumeTotal_flux = 0.0_wp
      IF (isRegistered('FrshFlux_VolumeTotal_Global')) THEN
      call levels_horizontal_mean( p_oce_sfc%FrshFlux_VolumeTotal, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & VolumeTotal_flux, lopenacc=lzacc)
      END IF
      monitor%FrshFlux_VolumeTotal = VolumeTotal_flux

      ! totalsnowfall
      totalsnowfall_flux = 0.0_wp
      IF (isRegistered('totalsnowfall_Global')) THEN
      call levels_horizontal_mean( ice%totalsnowfall, &
          & patch_2d%cells%area(:,:), &
          & owned_cells, &
          & totalsnowfall_flux, lopenacc=lzacc)
      END IF
      monitor%totalsnowfall = totalsnowfall_flux

      ! ice volume and extend
      ice_volume_nh = 0.0_wp
      IF (isRegistered('ice_volume_nh')) THEN
      ice_volume_nh = subset_sum( ice%vol(:,1,:)*p_diag%northernHemisphere(:,:), &
          & owned_cells, lopenacc=lzacc)
      END IF
      monitor%ice_volume_nh = ice_volume_nh/1.0e9_wp !scaling to km^3

      ice_volume_sh = 0.0_wp
      IF (isRegistered('ice_volume_sh')) THEN
      ice_volume_sh = subset_sum( ice%vol(:,1,:)*p_diag%southernHemisphere(:,:), &
          & owned_cells, lopenacc=lzacc)
      END IF
      monitor%ice_volume_sh = ice_volume_sh/1.0e9_wp !scaling to km^3

      ice_extent_nh = 0.0_wp
      IF (isRegistered('ice_extent_nh')) THEN
      ice_extent_nh = subset_sum( ice%concsum*p_diag%northernHemisphere*patch_2d%cells%area, &
          & owned_cells, lopenacc=lzacc)
      END IF
      monitor%ice_extent_nh = ice_extent_nh/1.0e6_wp !scaling to km^2

      ice_extent_sh = 0.0_wp
      IF (isRegistered('ice_extent_sh')) THEN
      ice_extent_sh = subset_sum( ice%concsum*p_diag%southernHemisphere*patch_2d%cells%area, &
          & owned_cells, lopenacc=lzacc)
      END IF
      monitor%ice_extent_sh = ice_extent_sh/1.0e6_wp !scaling to km^2

      w = p_diag%w
      IF ( ( vert_cor_type == 1 ) ) THEN
        w = p_diag%w_deriv
      ENDIF

      ! energy/enstrophy
      global_mean_potEnergy = 0.0_wp
      IF (isRegistered('pot_energy_global')) THEN
        IF (vert_cor_type .EQ. 0) THEN
          global_mean_potEnergy = potential_energy(&
              & w, &
  !TODO       & p_prog(nold(1))%h,&
              & sea_surface_height , & ! this is h_new, the old implementation used h_old
              & p_diag%rho, &
              & patch_3D%p_patch_1d(1)%del_zlev_i, &
              & patch_3D%p_patch_1d(1)%prism_volume, &
              & owned_cells, lacc=lzacc)
        ELSEIF (vert_cor_type .EQ. 1) THEN
          global_mean_potEnergy = potential_energy_zstar(&
              & w, &
              & sea_surface_height , & ! this is h_new, the old implementation used h_old
              & p_diag%rho, &
              & patch_3D%p_patch_1d(1)%del_zlev_i, &
              & ocean_state%p_prog(nnew(1))%stretch_c(:, :), &
              & patch_3D%p_patch_1d(1)%prism_volume, &
              & owned_cells, lacc=lzacc)
         END IF

      END IF
      monitor%pot_energy = global_mean_potEnergy

      global_mean_kinEnergy = 0.0_wp
      IF (isRegistered('kin_energy_global')) THEN
         global_mean_kinEnergy = total_mean( p_diag%kin, &
          & patch_3d%p_patch_1d(1)%prism_volume, &
          & owned_cells, lopenacc=lzacc )
      END IF
      monitor%kin_energy = global_mean_kinEnergy

      global_mean_totalEnergy = 0.0_wp
      IF (isRegistered('total_energy_global')) THEN
        ! use precomputed variables
        IF (isRegistered('kin_energy_global') .AND. isRegistered('pot_energy_global')) THEN
          global_mean_totalEnergy = global_mean_kinEnergy + global_mean_potEnergy
        END IF
      END IF
      monitor%total_energy = global_mean_totalEnergy

      global_mean_potEnstrophy = 0.0_wp
      IF (isRegistered('potential_enstrophy_global')) THEN
      END IF
      monitor%potential_enstrophy = global_mean_potEnstrophy
      !}}}

      IF ( isRegistered('delta_ice') .OR. isRegistered('delta_snow') .OR. &
           isRegistered('delta_thetao') .OR. &
           isRegistered('delta_so') .OR. &
           isRegistered('global_sltbasin') .OR. isRegistered('atlant_sltbasin') .OR. &
           isRegistered('pacind_sltbasin') .OR. &
           isRegistered('global_hfbasin') .OR. isRegistered('atlant_hfbasin') .OR. &
           isRegistered('pacind_hfbasin') ) THEN

      	IF (vert_cor_type .EQ. 0) THEN

          CALL diag_heat_salt_tendency(patch_3d, 2, ice, tracers(:,:,:,1), tracers(:,:,:,2), &
             p_diag%delta_ice,                                      &
             p_diag%delta_snow,                                     &
             p_diag%delta_thetao,                                   &
             p_diag%delta_so, lacc=lzacc)

        ELSEIF (vert_cor_type .EQ. 1) THEN

          CALL diag_heat_salt_tendency(patch_3d, 2, ice, tracers(:,:,:,1), tracers(:,:,:,2), &
             p_diag%delta_ice,                                      &
             p_diag%delta_snow,                                     &
             p_diag%delta_thetao,                                   &
             p_diag%delta_so, ocean_state%p_prog(nnew(1))%stretch_c(:, :), lacc=lzacc)

        ENDIF

      ENDIF

      ! calc moc each timestep from non-accumulated vertical veloc
      IF ( isRegistered('global_moc') .OR. isRegistered('atlant_moc') .OR. isRegistered('pacind_moc') .OR. &
           isRegistered('amoc26n') .OR. &
           isRegistered('global_hfl') .OR. isRegistered('atlant_hfl') .OR. isRegistered('pacind_hfl') .OR. &
           isRegistered('global_wfl') .OR. isRegistered('atlant_wfl') .OR. isRegistered('pacind_wfl') .OR. &
           isRegistered('global_sltbasin') .OR. isRegistered('atlant_sltbasin') .OR. isRegistered('pacind_sltbasin') .OR. &
           isRegistered('global_hfbasin') .OR. isRegistered('atlant_hfbasin') .OR. isRegistered('pacind_hfbasin') ) THEN
        CALL timer_start(timer_calc_moc)
        CALL calc_moc(patch_2d, patch_3d, &
             & w, &
             & p_oce_sfc%heatflux_total, &
             & p_oce_sfc%frshflux_volumetotal, &
             & p_diag%delta_thetao, &
             & p_diag%delta_so, &
             & p_diag%delta_snow, &
             & p_diag%delta_ice, &
             & p_diag%global_moc, &
             & p_diag%atlantic_moc, &
             & p_diag%pacific_moc,&
             & p_diag%global_hfl, &
             & p_diag%atlantic_hfl, &
             & p_diag%pacific_hfl, &
             & p_diag%global_wfl, &
             & p_diag%atlantic_wfl, &
             & p_diag%pacific_wfl, &
             & p_diag%global_hfbasin, &
             & p_diag%atlantic_hfbasin, &
             & p_diag%pacific_hfbasin, &
             & p_diag%global_sltbasin, &
             & p_diag%atlantic_sltbasin, &
             & p_diag%pacific_sltbasin, &
             & monitor%amoc26n, lacc=lzacc)


        CALL timer_stop(timer_calc_moc)
      ENDIF

      IF ( isRegistered('heat_content_liquid_water') .OR. isRegistered('heat_content_seaice') &
           .OR. isRegistered('heat_content_snow')   .OR. isRegistered('heat_content_total') &
           .OR. isRegistered('global_heat_content') .OR. isRegistered('global_heat_content_solid') ) THEN

      	IF (vert_cor_type .EQ. 0) THEN

          CALL calc_heat_content(patch_3d, prism_thickness, ice, tracers, &
             p_diag%heat_content_liquid_water, &
             p_diag%heat_content_seaice, &
             p_diag%heat_content_snow,&
             p_diag%heat_content_total, lacc=lzacc )

        ELSEIF (vert_cor_type .EQ. 1) THEN

          CALL calc_heat_content(patch_3d, prism_thickness, ice, tracers, &
             p_diag%heat_content_liquid_water, &
             p_diag%heat_content_seaice, &
             p_diag%heat_content_snow,&
             p_diag%heat_content_total, &
             ocean_state%p_prog(nnew(1))%stretch_c(:, :), lacc=lzacc )

        ENDIF

        ! global_heat_content for monitoring
        IF (isRegistered('global_heat_content')) THEN
          global_heat_content = 0.0_wp
          global_heat_content = global_sum_array(patch_2d%cells%area(:,:) * p_diag%heat_content_total(:,:) )
          monitor%global_heat_content = global_heat_content
        END IF

        ! global_heat_content_solid (snow and ice heat content) for monitoring
        IF (isRegistered('global_heat_content_solid')) THEN
          global_heat_content_solid = 0.0_wp
          global_heat_content_solid = global_sum_array( patch_2d%cells%area(:,:)* &
            &                      (p_diag%heat_content_seaice(:,:) + p_diag%heat_content_snow(:,:)) )
          monitor%global_heat_content_solid = global_heat_content_solid
        END IF


      ENDIF



      IF ( eddydiag .AND. &
         ( isRegistered('uT') .OR. isRegistered('uS') .OR. isRegistered('uR') .OR. &
           isRegistered('vT') .OR. isRegistered('vS') .OR. isRegistered('vR') .OR. &
           isRegistered('wT') .OR. isRegistered('wS') .OR. isRegistered('wR') .OR. &
           isRegistered('uu') .OR. isRegistered('uv') .OR. isRegistered('uw') .OR. &
           isRegistered('RR') .OR. isRegistered('SS') .OR. isRegistered('TT') .OR. &
           isRegistered('vv') .OR. isRegistered('ww') .OR. isRegistered('vw') .OR. &
           isRegistered('sigma0') .OR. isRegistered('hflR') .OR. isRegistered('fwR') .OR. &
           isRegistered('tauxU') .OR. isRegistered('tauyV') ) &
          ) THEN

        CALL calc_eddydiag(patch_3d, p_diag%u, p_diag%v, p_diag%w, p_diag%w_prismcenter  &
               ,tracers(:,:,:,1), tracers(:,:,:,2),p_diag%rho &  ! use in-situ density in the eddy products
               ,p_diag%uT, p_diag%uS, p_diag%uR, p_diag%uu    &
               ,p_diag%vT, p_diag%vS, p_diag%vR, p_diag%vv    &
               ,p_diag%wT, p_diag%wS, p_diag%wR, p_diag%ww    &
               ,p_diag%uv, p_diag%uw, p_diag%vw               &
               ,p_diag%RR, p_diag%SS, p_diag%TT, p_diag%sigma0 &
               ,p_diag%hflR, p_diag%fwR, p_diag%tauxU, p_diag%tauyV &
               ,p_oce_sfc%topbc_windstress_u, p_oce_sfc%topbc_windstress_v &
               ,p_oce_sfc%heatflux_total, p_oce_sfc%frshflux_volumetotal, lacc=lzacc )

      ENDIF

      IF (isRegistered('mld')) THEN

        CALL calc_mld(patch_3d, ocean_state%p_diag%mld, &
             ocean_state%p_diag%zgrad_rho,1,0.125_wp, lacc=lzacc)

        CALL dbg_print('Diag: mld',ocean_state%p_diag%mld, &
             str_module,4,in_subset=owned_cells)
      ENDIF

      IF (isRegistered('mlotst') .OR. isRegistered('mlotstsq') ) THEN

        CALL calc_mld(patch_3d, ocean_state%p_diag%mlotst, &
             ocean_state%p_diag%zgrad_rho,1,0.03_wp, lacc=lzacc)

        CALL dbg_print('Diag: mlotst',ocean_state%p_diag%mlotst, &
             str_module,4,in_subset=owned_cells)

        IF (isRegistered('mlotstsq')) THEN

          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          ocean_state%p_diag%mlotstsq= &
               ocean_state%p_diag%mlotst*ocean_state%p_diag%mlotst
          !$ACC END KERNELS
          !$ACC WAIT(1)

          CALL dbg_print('Diag: mlotstsq',ocean_state%p_diag%mlotstsq, &
               str_module,4,in_subset=owned_cells)
        ENDIF

      ENDIF

      IF (isRegistered('mlotst10') .OR. isRegistered('mlotst10sq') ) THEN

        CALL calc_mld(patch_3d, ocean_state%p_diag%mlotst10, &
             ocean_state%p_diag%zgrad_rho,get_level_index_by_depth(patch_3d, 10.0_wp),0.03_wp, lacc=lzacc)

        CALL dbg_print('Diag: mlotst10',ocean_state%p_diag%mlotst10, &
             str_module,4,in_subset=owned_cells)

        IF (isRegistered('mlotst10sq')) THEN

          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          ocean_state%p_diag%mlotst10sq= &
               ocean_state%p_diag%mlotst10*ocean_state%p_diag%mlotst10
          !$ACC END KERNELS
          !$ACC WAIT(1)

          CALL dbg_print('Diag: mlotst10sq',ocean_state%p_diag%mlotst10sq, &
               str_module,4,in_subset=owned_cells)
        ENDIF

      ENDIF

      IF (isRegistered('condep')) THEN

        CALL calc_condep(patch_3d, ocean_state%p_diag%condep, &
             ocean_state%p_diag%zgrad_rho, lacc=lzacc)

        CALL dbg_print('Diag: condep',ocean_state%p_diag%condep,str_module,4, &
             in_subset=owned_cells)
      ENDIF

      IF (isRegistered('ssh')) THEN
        IF (vert_cor_type .EQ. 1) THEN
          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          p_diag%ssh = sea_surface_height + ice%draftave
          !$ACC END KERNELS
          !$ACC WAIT(1)
        ELSE
          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          p_diag%ssh = sea_surface_height
          !$ACC END KERNELS
          !$ACC WAIT(1)
        ENDIF
      ENDIF

      CALL dbg_print('Diag: mld',p_diag%mld,str_module,4,in_subset=owned_cells)

      ! square of ssh
      p_diag%zos_square = merge(sea_surface_height*sea_surface_height,0.0_wp,isRegistered('zos_square'))

      monitor%gibraltar = merge( section_flux(oce_sections(1),normal_veloc)*OceanReferenceDensity, &
          &                      0.0_wp, &
          &                      isRegistered('gibraltar'))
      monitor%denmark_strait = merge( section_flux(oce_sections(2),normal_veloc)*OceanReferenceDensity, &
          &                      0.0_wp, &
          &                      isRegistered('denmark_strait'))
      monitor%drake_passage = merge( section_flux(oce_sections(3),normal_veloc)*OceanReferenceDensity, &
          &                      0.0_wp, &
          &                      isRegistered('drake_passage'))
      monitor%indonesian_throughflow = merge( section_flux(oce_sections(4),normal_veloc)*OceanReferenceDensity, &
          &                      0.0_wp, &
          &                      isRegistered('indonesian_throughflow'))
      monitor%scotland_iceland = merge( section_flux(oce_sections(5),normal_veloc)*OceanReferenceDensity, &
          &                      0.0_wp, &
          &                      isRegistered('scotland_iceland'))
      monitor%mozambique = merge( section_flux(oce_sections(6),normal_veloc)*OceanReferenceDensity, &
          &                      0.0_wp, &
          &                      isRegistered('mozambique'))
      monitor%framStrait = merge( section_flux(oce_sections(7),normal_veloc)*OceanReferenceDensity, &
          &                      0.0_wp, &
          &                      isRegistered('framStrait'))
      monitor%beringStrait = merge( section_flux(oce_sections(8),normal_veloc)*OceanReferenceDensity, &
          &                      0.0_wp, &
          &                      isRegistered('beringStrait'))
      monitor%barentsOpening = merge( section_flux(oce_sections(9),normal_veloc)*OceanReferenceDensity, &
          &                      0.0_wp, &
          &                      isRegistered('barentsOpening'))
      monitor%ice_framStrait = merge(section_ice_flux(oce_sections(7), ice%hi*ice%conc, ice%vn_e), &
          &                      0.0_wp, &
          &                      isRegistered('ice_framStrait'))

      IF (isRegistered('verticallyTotal_mass_flux_e')) THEN
        IF(lzacc) CALL finish("calc_fast_oce_diagnostic", "verticallyIntegrated_field is ported to GPU &
                              but not checked whether it gives correct results")
        CALL verticallyIntegrated_field(ocean_state%p_diag%verticallyTotal_mass_flux_e, &
          & ocean_state%p_diag%mass_flx_e, owned_edges, lacc=lzacc)
        CALL dbg_print('Total_mass_flux_e ', ocean_state%p_diag%verticallyTotal_mass_flux_e, &
           str_module, 1, in_subset=owned_edges)

      ENDIF

!TODO       CASE (10)
!TODO         monitor%agulhas                = sflux*OceanReferenceDensity
!TODO       CASE (11)
!TODO         monitor%agulhas_long           = sflux*OceanReferenceDensity
!TODO       CASE (12)
!TODO         monitor%agulhas_longer         = sflux*OceanReferenceDensity
!TODO       CASE (13)
!TODO         monitor%florida_strait         = sflux*OceanReferenceDensity
    END SELECT
  END SUBROUTINE calc_fast_oce_diagnostics
  !-------------------------------------------------------------------------

  !TODO potential_energy(&
  !TODO     & p_diag%w,p_prog(nold(1))%h,&
  !TODO     & p_diag%rho, &
  !TODO     & patch_3D%p_patch_1d(1))%del_zlev_i, &
  !TODO     & patch_3D%p_patch_1d(1)%prism_volume, &
  !TODO     & owned_cells)
  !-------------------------------------------------------------------------
  FUNCTION potential_energy(w,h,rho,del_zlev_i,weights,in_subset,lacc)
    REAL(wp), INTENT(IN) :: w(:,:,:)
    REAL(wp), INTENT(IN) :: h(:,:)
    REAL(wp), INTENT(IN) :: rho(:,:,:)
    REAL(wp), INTENT(IN) :: del_zlev_i(:)
    REAL(wp), INTENT(IN) :: weights(:,:,:)
    TYPE(t_subset_range), INTENT(IN) :: in_subset
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    REAL(wp) :: potential_energy
#define VerticalDim_Position 2

    REAL(wp), ALLOCATABLE :: sum_value(:,:), sum_weight(:,:), total_weight(:), total_sum(:)
    INTEGER :: block, level, start_index, end_index, idx, start_vertical, end_vertical
    INTEGER :: allocated_levels, no_of_threads, myThreadNo
    REAL(wp) :: z_w, totalSum, totalWeight

    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':potential_energy'
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (in_subset%no_of_holes > 0) CALL warning(method_name, "there are holes in the subset")

    no_of_threads = 1
    myThreadNo = 0
#ifdef _OPENMP
    no_of_threads = omp_get_max_threads()
#endif

    allocated_levels = SIZE(w,VerticalDim_Position)
    ALLOCATE( sum_value(allocated_levels, 0:no_of_threads-1), &
      & sum_weight(allocated_levels, 0:no_of_threads-1), &
      & total_weight(allocated_levels), total_sum(allocated_levels) )

    start_vertical = 1
    end_vertical = SIZE(w, VerticalDim_Position)

    IF (start_vertical > end_vertical) &
      & CALL finish(method_name, "start_vertical > end_vertical")
    IF ( allocated_levels < end_vertical) &
      & CALL finish(method_name, "allocated_levels < end_vertical")

    !$ACC DATA PRESENT(w, h, rho, del_zlev_i, weights, in_subset%vertical_levels) &
    !$ACC   CREATE(sum_value, sum_weight, total_sum, total_weight) IF(lzacc)

!ICON_OMP_PARALLEL PRIVATE(myThreadNo)
#ifdef _OPENMP
    myThreadNo = omp_get_thread_num()
#endif
!ICON_OMP_SINGLE
#ifdef _OPENMP
    no_of_threads = OMP_GET_NUM_THREADS()
#endif
!ICON_OMP_END_SINGLE NOWAIT
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    sum_value(:,  myThreadNo) = 0.0_wp
    sum_weight(:,  myThreadNo) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)
    IF (ASSOCIATED(in_subset%vertical_levels)) THEN
!ICON_OMP_DO PRIVATE(block, start_index, end_index, idx)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO idx = start_index, end_index
          !$ACC LOOP SEQ
          DO level = start_vertical, MIN(end_vertical, in_subset%vertical_levels(idx,block)) - 1
            z_w = MERGE( &
              & (w(idx,level,block)*h(idx,block) &
              &  + w(idx,level+1,block)*0.5_wp*del_zlev_i(level)) &
              & /(0.5_wp*del_zlev_i(level)+h(idx,block)) &
              & , &
              & (w(idx,level,block)*del_zlev_i(level) &
              &  + w(idx,level+1,block)*del_zlev_i(level+1)) &
              & /(del_zlev_i(level)+del_zlev_i(level+1)) &
              & , &
              & 1 .EQ. level)

            sum_value(level, myThreadNo)  = sum_value(level, myThreadNo) + &
              & grav*z_w*rho(idx, level, block) * weights(idx, level, block)

            sum_weight(level, myThreadNo)  = sum_weight(level, myThreadNo) + weights(idx, level, block)

          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
      !$ACC WAIT(1)
!ICON_OMP_END_DO

    ELSE ! no in_subset%vertical_levels

!ICON_OMP_DO PRIVATE(block, start_index, end_index)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO idx = start_index, end_index
          ! since we have the same numbder of vertical layers, the weight is the same
          ! for all levels. Compute it only for the first level, and then copy it
          !$ACC LOOP SEQ
          DO level = start_vertical, end_vertical - 1
            z_w = MERGE( &
              & (w(idx,level,block)*h(idx,block) &
              &  + w(idx,level+1,block)*0.5_wp*del_zlev_i(level)) &
              & /(0.5_wp*del_zlev_i(level)+h(idx,block)) &
              & , &
              & (w(idx,level,block)*del_zlev_i(level) &
              &  + w(idx,level+1,block)*del_zlev_i(level+1)) &
              & /(del_zlev_i(level)+del_zlev_i(level+1)) &
              & , &
              & 1 .EQ. level)

            sum_value(level, myThreadNo)  = sum_value(level, myThreadNo) + &
              & grav*z_w*rho(idx, level, block) * weights(idx,level, block)
            sum_weight(level, myThreadNo)  = sum_weight(start_vertical, myThreadNo) + weights(idx, level, block)
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
      !$ACC WAIT(1)
!ICON_OMP_END_DO

    ENDIF
!ICON_OMP_END_PARALLEL

    ! gather the total level sum of this process in total_sum(level)
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    total_sum(:)     = 0.0_wp
    total_weight(:) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)
    DO myThreadNo=0, no_of_threads-1
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO level = start_vertical, end_vertical - 1
        ! write(0,*) myThreadNo, level, " sum=", sum_value(level, myThreadNo), sum_weight(level, myThreadNo)
        total_sum(level)    = total_sum(level)    + sum_value(level, myThreadNo)
        total_weight(level) = total_weight(level) + sum_weight(level, myThreadNo)
      ENDDO
      !$ACC END PARALLEL LOOP
    ENDDO
    !$ACC WAIT(1)

    ! Collect the value and weight sums (at all procs)
    CALL gather_sums(total_sum, total_weight, lopenacc=lzacc)


    totalSum = 0.0_wp
    totalWeight = 0.0_wp
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO level = start_vertical, end_vertical
      totalSum    = totalSum    + total_sum(level)
      totalWeight = totalWeight + total_weight(level)
    ENDDO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    !$ACC END DATA

    DEALLOCATE(sum_value, sum_weight)
    DEALLOCATE(total_weight)
    DEALLOCATE(total_sum)

    potential_energy = totalSum / totalWeight
  END FUNCTION potential_energy


  FUNCTION potential_energy_zstar(w,h,rho,del_zlev_i, stretch, weights,in_subset,lacc)
    REAL(wp), INTENT(IN) :: w(:,:,:)
    REAL(wp), INTENT(IN) :: h(:,:)
    REAL(wp), INTENT(IN) :: rho(:,:,:)
    REAL(wp), INTENT(IN) :: del_zlev_i(:)
    REAL(wp), INTENT(IN) :: stretch(:, :)
    REAL(wp), INTENT(IN) :: weights(:,:,:)
    TYPE(t_subset_range), INTENT(IN) :: in_subset
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    REAL(wp) :: potential_energy_zstar
#define VerticalDim_Position 2

    REAL(wp), ALLOCATABLE :: sum_value(:,:), sum_weight(:,:), total_weight(:), total_sum(:)
    INTEGER :: block, level, start_index, end_index, idx, start_vertical, end_vertical
    INTEGER :: allocated_levels, no_of_threads, myThreadNo
    REAL(wp) :: z_w, totalSum, totalWeight

    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':potential_energy'
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (in_subset%no_of_holes > 0) CALL warning(method_name, "there are holes in the subset")

    no_of_threads = 1
    myThreadNo = 0
#ifdef _OPENMP
    no_of_threads = omp_get_max_threads()
#endif

    allocated_levels = SIZE(w,VerticalDim_Position)
    ALLOCATE( sum_value(allocated_levels, 0:no_of_threads-1), &
      & sum_weight(allocated_levels, 0:no_of_threads-1), &
      & total_weight(allocated_levels), total_sum(allocated_levels) )

    start_vertical = 1
    end_vertical = SIZE(w, VerticalDim_Position)

    IF (start_vertical > end_vertical) &
      & CALL finish(method_name, "start_vertical > end_vertical")
    IF ( allocated_levels < end_vertical) &
      & CALL finish(method_name, "allocated_levels < end_vertical")

    !$ACC DATA PRESENT(w, h, rho, del_zlev_i, stretch, weights, in_subset%vertical_levels) &
    !$ACC   CREATE(sum_value, sum_weight, total_sum, total_weight) IF(lzacc)

    !ICON_OMP_PARALLEL PRIVATE(myThreadNo)
#ifdef _OPENMP
    myThreadNo = omp_get_thread_num()
#endif
    !ICON_OMP_SINGLE
#ifdef _OPENMP
    no_of_threads = OMP_GET_NUM_THREADS()
#endif
    !ICON_OMP_END_SINGLE NOWAIT
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    sum_value(:,  myThreadNo) = 0.0_wp
    sum_weight(:,  myThreadNo) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)
    IF (ASSOCIATED(in_subset%vertical_levels)) THEN
    !ICON_OMP_DO PRIVATE(block, start_index, end_index, idx)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO idx = start_index, end_index
          !$ACC LOOP SEQ
          DO level = start_vertical, MIN(end_vertical, in_subset%vertical_levels(idx,block)) - 1
            z_w = (w(idx,level,block)*del_zlev_i(level) &
              &  + w(idx,level+1,block)*del_zlev_i(level+1)) &
              & /(del_zlev_i(level)+del_zlev_i(level+1))

            sum_value(level, myThreadNo)  = sum_value(level, myThreadNo) + &
              & grav*z_w*rho(idx, level, block) * weights(idx, level, block)*stretch(idx, block)

            sum_weight(level, myThreadNo)  = sum_weight(level, myThreadNo) + weights(idx, level, block)

          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
      !$ACC WAIT(1)
  !ICON_OMP_END_DO

    ELSE ! no in_subset%vertical_levels

  !ICON_OMP_DO PRIVATE(block, start_index, end_index)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO idx = start_index, end_index
          ! since we have the same numbder of vertical layers, the weight is the same
          ! for all levels. Compute it only for the first level, and then copy it
          !$ACC LOOP SEQ
          DO level = start_vertical, end_vertical - 1
            z_w = (w(idx,level,block)*del_zlev_i(level) &
              &  + w(idx,level+1,block)*del_zlev_i(level+1)) &
              & /(del_zlev_i(level)+del_zlev_i(level+1))

            sum_value(level, myThreadNo)  = sum_value(level, myThreadNo) + &
              & grav*z_w*rho(idx, level, block) * weights(idx,level,block)*stretch(idx, block)
            sum_weight(level, myThreadNo)  = sum_weight(start_vertical, myThreadNo) + weights(idx, level, block)
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
      !$ACC WAIT(1)
  !ICON_OMP_END_DO

    ENDIF
  !ICON_OMP_END_PARALLEL

    ! gather the total level sum of this process in total_sum(level)
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    total_sum(:)     = 0.0_wp
    total_weight(:) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)
    DO myThreadNo=0, no_of_threads-1
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO level = start_vertical, end_vertical - 1
        ! write(0,*) myThreadNo, level, " sum=", sum_value(level, myThreadNo), sum_weight(level, myThreadNo)
        total_sum(level)    = total_sum(level)    + sum_value(level, myThreadNo)
        total_weight(level) = total_weight(level) + sum_weight(level, myThreadNo)
      ENDDO
      !$ACC END PARALLEL LOOP
    ENDDO
    !$ACC WAIT(1)

    ! Collect the value and weight sums (at all procs)
    CALL gather_sums(total_sum, total_weight, lopenacc=lzacc)


    totalSum = 0.0_wp
    totalWeight = 0.0_wp
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO level = start_vertical, end_vertical
      totalSum    = totalSum    + total_sum(level)
      totalWeight = totalWeight + total_weight(level)
    ENDDO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    !$ACC END DATA

    DEALLOCATE(sum_value, sum_weight)
    DEALLOCATE(total_weight)
    DEALLOCATE(total_sum)

    potential_energy_zstar = totalSum / totalWeight
  END FUNCTION potential_energy_zstar


  !-------------------------------------------------------------------------
  REAL(wp) FUNCTION section_flux(in_oce_section, velocity_values)
    TYPE(t_oce_section) :: in_oce_section
    REAL(wp), POINTER :: velocity_values(:,:,:)

    INTEGER :: i, k, edge_idx, edge_block
    REAL(wp) :: oriented_length
    REAL(wp), ALLOCATABLE :: flux_weights(:,:)
    TYPE(t_grid_edges), POINTER ::  edges
    TYPE(t_patch_vert),POINTER :: patch_vertical

    CHARACTER(LEN=*), PARAMETER :: method_name='mo_ocean_diagnostics:section_flux'

    edges          => in_oce_section%subset%patch%edges
    patch_vertical => in_oce_section%subset%patch_3d%p_patch_1d(1)

    ! calculate weights
    ! flux_weights can also be preallocated
    ALLOCATE(flux_weights(n_zlev, MAX(in_oce_section%subset%SIZE, 1)))
    flux_weights(:,:) = 0.0_wp
    DO i=1, in_oce_section%subset%SIZE

      edge_idx   = in_oce_section%subset%idx(i)
      edge_block = in_oce_section%subset%BLOCK(i)
      oriented_length = edges%primal_edge_length(edge_idx, edge_block) * &
        & in_oce_section%orientation(i) ! this can also be pre-calculated and stored in in_oce_section%orientation

      !write(0,*) "oriented_length:",  oriented_length

      DO k=1, n_zlev
        flux_weights(k, i) = patch_vertical%prism_thick_e(edge_idx, k, edge_block) * oriented_length ! maybe also use slm
        !write(0,*) i, k, in_oce_section%subset%name, " flux_weights:",  flux_weights(k, i), &
        !  & patch_vertical%prism_thick_e(edge_idx, k, edge_block)
        !write(0,*) i, k, in_oce_section%subset%name, " velocity_value:", velocity_values(edge_idx, k, edge_block)
      ENDDO

    ENDDO


    section_flux = subset_sum(                           &
      & values                 = velocity_values,        &
      & indexed_subset         = in_oce_section%subset,  &
      & subset_indexed_weights = flux_weights)

    DEALLOCATE(flux_weights)

    !write(0,*) get_my_mpi_work_id(), ": section_flux on subset ", in_oce_section%subset%name, ":", &
    !  & section_flux, in_oce_section%subset%size

  END FUNCTION section_flux
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  REAL(wp) FUNCTION section_ice_flux(in_oce_section, ice_hmean, ice_vel)
    TYPE(t_oce_section)  :: in_oce_section
    REAL(wp), INTENT(IN) :: ice_hmean(:,:,:)
    REAL(wp), POINTER    :: ice_vel(:,:)

    TYPE(t_grid_edges), POINTER ::  edges
    INTEGER :: i, k, edge_idx, edge_block, cell_idx(2), cell_block(2)
    REAL(wp) :: vel_vals, sum_value, oriented_length
    REAL(wp), ALLOCATABLE :: fluxes(:,:)
    INTEGER :: communicator

    CHARACTER(LEN=*), PARAMETER :: method_name='mo_ocean_diagnostics:section_ice_flux'

    edges          => in_oce_section%subset%patch%edges

    ! calculate weights
    ! fluxes can also be preallocated
    ALLOCATE(fluxes(kice, MAX(in_oce_section%subset%SIZE, 1)))
    fluxes(:,:) = 0.0_wp
    sum_value  = 0.0_wp

    DO i=1, in_oce_section%subset%SIZE

      edge_idx   = in_oce_section%subset%idx(i)
      edge_block = in_oce_section%subset%BLOCK(i)

      cell_idx   = edges%cell_idx(edge_idx, edge_block,:)
      cell_block = edges%cell_blk(edge_idx, edge_block,:)

      oriented_length = edges%primal_edge_length(edge_idx, edge_block) * &
        & in_oce_section%orientation(i) ! this can also be pre-calculated and stored in in_oce_section%orientation

      vel_vals = ice_vel(edge_idx, edge_block)

      ! compute the first order upwind flux using cell-centered ice_hmean vals and edge-centered vel_vals
      ! Same as in upwind_hflux_ice which is used to calculate ice advection
      DO k=1,kice
          fluxes(k, i) = laxfr_upflux( vel_vals, ice_hmean(cell_idx(1),k,cell_block(1)), &
          &                    ice_hmean(cell_idx(2),k,cell_block(2)) ) * oriented_length
          !write(0,*) i, k, in_oce_section%subset%name, " fluxes:",  fluxes(k, i), &
          !  & patch_vertical%prism_thick_e(edge_idx, k, edge_block)
          !write(0,*) i, k, in_oce_section%subset%name, " velocity_value:", velocity_values(edge_idx, k, edge_block)
      ENDDO

    ENDDO

    sum_value = sum(fluxes(:,:))

    DEALLOCATE(fluxes)

    ! the global min, max is avaliable only to stdio process
    IF (my_process_is_mpi_parallel()) THEN

      communicator = in_oce_section%subset%patch%work_communicator
      ! these are avaliable to all processes
      section_ice_flux = p_sum( sum_value,  comm=communicator)

    ELSE

      section_ice_flux = sum_value

    ENDIF

!    write(0,*) get_my_mpi_work_id(), ": section_ice_flux on subset ", in_oce_section%subset%name, ":", &
!      & section_ice_flux, in_oce_section%subset%size

  END FUNCTION section_ice_flux
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !
  !!  Calculation of meridional overturning circulation (MOC)
  !
  !   Calculation of meridional overturning circulation for different basins
  !   (Atlantic, Pacific, Indian, global)
  !>
  !!
  !!  based on code from MPIOM
  !
  ! TODO: implement variable output dimension (1 deg resolution) and smoothing extent
  ! TODO: calculate the 1 deg resolution meridional distance
  !!

  SUBROUTINE calc_moc_hfl_internal_vector (patch_2d, patch_3d, w, heatflux_total, &
             frshflux_volumetotal, delta_thetao, delta_so, delta_snow, delta_ice, &
             global_moc, atlant_moc, pacind_moc, global_hfl, atlant_hfl, pacind_hfl, &
             global_wfl, atlant_wfl, pacind_wfl, &
             global_hfbasin, atlant_hfbasin, pacind_hfbasin, &
             global_sltbasin, atlant_sltbasin, pacind_sltbasin, amoc26n, lacc)

    TYPE(t_patch),    TARGET, INTENT(in)  :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d

    REAL(wp), INTENT(in)  :: w(:,:,:)   ! vertical velocity (nproma,nlev+1,alloc_cell_blocks)
    REAL(wp), INTENT(in)  :: heatflux_total(:,:)   ! heatflux_total (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(in)  :: frshflux_volumetotal(:,:)   ! fw flux (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(in)  :: delta_snow(:,:)       ! tendency of snow thickness (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(in)  :: delta_ice(:,:)        ! tendendy of ice  thickness (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(in)  :: delta_thetao(:,:,:)   ! temperature tendency (nproma,nlev+1,alloc_cell_blocks)
    REAL(wp), INTENT(in)  :: delta_so(:,:,:)   ! salinity tendency (nproma,nlev+1,alloc_cell_blocks)
    REAL(wp), INTENT(inout)  :: amoc26n(:)

    REAL(wp), INTENT(out) :: global_moc(:,:), atlant_moc(:,:), pacind_moc(:,:) ! (n_zlev,nlat_moc)

    ! implied ocean heat transport calculated from surface fluxes
    REAL(wp), INTENT(out) :: global_hfl(:,:), atlant_hfl(:,:), pacind_hfl(:,:) ! (1,nlat_moc)

    ! implied ocean fw transport calculated from surface fluxes
    REAL(wp), INTENT(out) :: global_wfl(:,:), atlant_wfl(:,:), pacind_wfl(:,:) ! (1,nlat_moc)

    ! northward ocean heat transport calculated from tendencies
    REAL(wp), INTENT(out) :: global_hfbasin(:,:), atlant_hfbasin(:,:), pacind_hfbasin(:,:) ! (1,nlat_moc)

    ! northward ocean salt transport calculated from tendencies
    REAL(wp), INTENT(out) :: global_sltbasin(:,:), atlant_sltbasin(:,:), pacind_sltbasin(:,:) ! (1,nlat_moc)
    LOGICAL, INTENT(in), OPTIONAL :: lacc

    ! local variables
    INTEGER, PARAMETER ::  latSmooth = 3   !  latitudinal smoothing area is 2*jbrei-1 rows of 1 deg
    REAL(wp), PARAMETER :: smoothWeight = 1.0_wp / REAL(2*latSmooth + 1, wp)
    INTEGER :: jb, level, start_index, end_index, jc, ilat, l, ic
    INTEGER :: mpi_comm

    REAL(wp) :: lat, deltaMoc, deltahfl, deltawfl, deltahfbasin, deltasltbasin
    REAL(wp) :: mocs_lat(nproma,patch_2d%nblks_c,9)
    REAL(wp) :: mocs_2d(nlat_moc,12)
    REAL(wp) :: mocs_3d(n_zlev,nlat_moc,3)

    INTEGER :: dist, weight

    REAL(wp) :: factor_to_sv

    TYPE(t_subset_range), POINTER :: cells
    INTEGER, POINTER, CONTIGUOUS :: vertical_levels(:,:)
    INTEGER, POINTER, CONTIGUOUS :: basin_c(:,:)
    REAL(wp), POINTER, CONTIGUOUS :: area(:,:)
    REAL(wp), POINTER, CONTIGUOUS :: wet_c(:,:,:)

    TYPE t_idxlist_weighted
      INTEGER, ALLOCATABLE :: idx(:)
      INTEGER, ALLOCATABLE :: blk(:)
      REAL(wp), ALLOCATABLE :: wgt(:)
      INTEGER :: len
    END TYPE

    INTEGER, ALLOCATABLE :: itmp(:)
    INTEGER, ALLOCATABLE :: btmp(:)
    REAL(wp), ALLOCATABLE :: wtmp(:)

    TYPE(t_idxlist_weighted), ALLOCATABLE, SAVE :: idxlist(:)
    TYPE(t_patch_3d), POINTER, SAVE :: idxlist_patch

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_ocean_diagnostics:calc_moc'

    LOGICAL  :: lzacc
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish(routine, 'OpenACC version currently not implemented')
#endif

    mpi_comm = MERGE(p_comm_work_test, p_comm_work, p_test_run)

    mocs_2d(:,:) = 0._wp
    mocs_3d(:,:,:) = 0._wp

    global_moc(:,:) = 0.0_wp
    pacind_moc(:,:) = 0.0_wp
    atlant_moc(:,:) = 0.0_wp

    global_hfl(:,:) = 0.0_wp
    pacind_hfl(:,:) = 0.0_wp
    atlant_hfl(:,:) = 0.0_wp

    global_wfl(:,:) = 0.0_wp
    pacind_wfl(:,:) = 0.0_wp
    atlant_wfl(:,:) = 0.0_wp

    global_hfbasin(:,:) = 0.0_wp
    pacind_hfbasin(:,:) = 0.0_wp
    atlant_hfbasin(:,:) = 0.0_wp

    global_sltbasin(:,:) = 0.0_wp
    pacind_sltbasin(:,:) = 0.0_wp
    atlant_sltbasin(:,:) = 0.0_wp


    ! limit cells to in-domain because of summation
    cells   => patch_2d%cells%in_domain

    vertical_levels => cells%vertical_levels
    area => patch_2d%cells%area
    basin_c => patch_3d%basin_c
    wet_c => patch_3d%wet_c

    ! lat: corresponding latitude row of 1 deg extension
    !            1 south pole
    ! nlat_moc=180 north pole

    ! distribute MOC over (2*jbrei)+1 latitude rows
    !  - no weighting with latitudes done
    !  - lat: index of 180 X 1 deg meridional resolution
    IF (.NOT. ALLOCATED(idxlist)) THEN
      ALLOCATE(idxlist(nlat_moc))
      DO ilat = 1, nlat_moc

        ALLOCATE( &
            & idxlist(ilat)%idx(nproma * patch_2d%nblks_c), &
            & idxlist(ilat)%blk(nproma * patch_2d%nblks_c), &
            & idxlist(ilat)%wgt(nproma * patch_2d%nblks_c) &
          )

        ic = 0

        DO jb = cells%start_block, cells%end_block
          CALL get_index_range(cells, jb, start_index, end_index)

          DO jc = start_index, end_index
            IF (vertical_levels(jc,jb) < 1) CYCLE

            lat = patch_2d%cells%center(jc,jb)%lat * rad2deg

            dist = ABS(ilat - NINT(REAL(nlat_moc, wp)*0.5_wp + lat))

            IF (dist <= latSmooth) THEN
              ! Smeared contributions from points close to the poles are added multiple times in the
              ! scalar version. We replicate this here, although a change of the smoothWeight might
              ! be more appropriate, as it does not bias the smoothing towards the poles.
              IF (ilat == 1 .OR. ilat == nlat_moc) THEN
                weight = ABS(dist - latSmooth) + 1
              ELSE
                weight = 1
              END IF
              ic = ic + 1
              idxlist(ilat)%idx(ic) = jc
              idxlist(ilat)%blk(ic) = jb
              idxlist(ilat)%wgt(ic) = weight * smoothWeight
            END IF
          END DO
        END DO

        ALLOCATE(itmp(ic), btmp(ic), wtmp(ic))
        itmp(:) = idxlist(ilat)%idx(1:ic)
        btmp(:) = idxlist(ilat)%blk(1:ic)
        wtmp(:) = idxlist(ilat)%wgt(1:ic)

        CALL MOVE_ALLOC(itmp, idxlist(ilat)%idx)
        CALL MOVE_ALLOC(btmp, idxlist(ilat)%blk)
        CALL MOVE_ALLOC(wtmp, idxlist(ilat)%wgt)
        idxlist(ilat)%len = ic

      END DO

      idxlist_patch => patch_3d
    END IF

    IF (.NOT. ASSOCIATED(idxlist_patch, target=patch_3d)) THEN
      CALL finish(routine, 'Can only be called for a single patch due to index-list caching.')
    END IF


    DO jb = cells%start_block, cells%end_block
      CALL get_index_range(cells, jb, start_index, end_index)

      ! Surface variables
      !NEC$ nolstval
      DO jc = start_index, end_index
        IF (1 <= vertical_levels(jc,jb)) THEN
          deltahfl = area(jc,jb) * heatflux_total(jc,jb) * wet_c(jc,1,jb)

          deltawfl = area(jc,jb) * frshflux_volumetotal(jc,jb) * wet_c(jc,1,jb)

          ! Global HFL
          mocs_lat(jc,jb,1) = - deltahfl
          ! Global WFL
          mocs_lat(jc,jb,4) = - deltawfl

          IF (basin_c(jc,jb) == 1) THEN
            ! Atlantic HFL
            mocs_lat(jc,jb,2) = - deltahfl
            ! Atlantic WFL
            mocs_lat(jc,jb,5) = - deltawfl
          ELSE
            mocs_lat(jc,jb,2) = 0._wp
            mocs_lat(jc,jb,5) = 0._wp
          END IF

          IF (basin_c(jc,jb) >= 2) THEN
            ! Pacind HFL
            mocs_lat(jc,jb,3) = - deltahfl
            ! Pacind WFL
            mocs_lat(jc,jb,6) = - deltawfl
          ELSE
            mocs_lat(jc,jb,3) = 0._wp
            mocs_lat(jc,jb,6) = 0._wp
          ENDIF
        ELSE
          mocs_lat(jc,jb,1:6) = 0._wp
        END IF
      END DO ! jc
    END DO ! jb

    ! Collapse the contribution of each cell.
    DO ilat = 1, nlat_moc
      DO ic = 1, idxlist(ilat)%len
        !NEC$ assoc
        DO l = 1, 6
          mocs_2d(ilat,l) = mocs_2d(ilat,l) &
              & + mocs_lat(idxlist(ilat)%idx(ic), idxlist(ilat)%blk(ic), l) * idxlist(ilat)%wgt(ic)
        END DO
      END DO
    END DO

    DO level = 1, n_zlev
      DO jb = cells%start_block, cells%end_block
        CALL get_index_range(cells, jb, start_index, end_index)

        !NEC$ nolstval
        DO jc = start_index, end_index
          IF (level <= vertical_levels(jc,jb)) THEN
            deltaMoc = area(jc,jb) * OceanReferenceDensity * w(jc,level,jb)

            deltahfbasin = area(jc,jb) * delta_thetao(jc,level,jb)

            deltasltbasin = area(jc,jb) * delta_so(jc,level,jb)

            IF (level .EQ. 1) THEN
              deltahfbasin = deltahfbasin                                &
                   + area(jc,jb) * ( delta_ice(jc,jb) + delta_snow(jc,jb) )
            ENDIF

            ! Global MOC
            mocs_lat(jc,jb,1) = - deltaMoc
            ! Global HFBasin
            mocs_lat(jc,jb,4) = - deltahfbasin
            ! Global SLTBasin
            mocs_lat(jc,jb,7) = - deltasltbasin

            IF (basin_c(jc,jb) == 1) THEN
              ! Atlantic MOC
              mocs_lat(jc,jb,2) = - deltaMoc
              ! Atlantic HFBasin
              mocs_lat(jc,jb,5) = - deltahfbasin
              ! Atlantic SLTBasin
              mocs_lat(jc,jb,8) = - deltasltbasin
            ELSE
              mocs_lat(jc,jb,2) = 0._wp
              mocs_lat(jc,jb,5) = 0._wp
              mocs_lat(jc,jb,8) = 0._wp
            END IF

            IF (basin_c(jc,jb) >= 2) THEN
              ! Pacind MOC
              mocs_lat(jc,jb,3) = - deltaMoc
              ! Pacind HFBasin
              mocs_lat(jc,jb,6) = - deltahfbasin
              ! Pacind SLTBasin
              mocs_lat(jc,jb,9) = - deltasltbasin
            ELSE
              mocs_lat(jc,jb,3) = 0._wp
              mocs_lat(jc,jb,6) = 0._wp
              mocs_lat(jc,jb,9) = 0._wp
            END IF
          ELSE
            !NEC$ unroll_complete
            mocs_lat(jc,jb,1:9) = 0._wp
          END IF
        END DO ! jc
      END DO ! jb

      ! Collapse the contribution of each cell.
      DO ilat = 1, nlat_moc
        DO ic = 1, idxlist(ilat)%len
          ! 3D fields
          !NEC$ assoc
          DO l = 1, 3
            mocs_3d(level,ilat,l) = mocs_3d(level,ilat,l) &
                & + mocs_lat(idxlist(ilat)%idx(ic), idxlist(ilat)%blk(ic), l) * idxlist(ilat)%wgt(ic)
          END DO

          ! The hfbasin and sltbasin vars are summed over levels.
          ! Results are located in mocs_2d(1:nlat_moc,7:12)
          !NEC$ assoc
          DO l = 4, 9
            mocs_2d(ilat,l+3) = mocs_2d(ilat,l+3) &
                & + mocs_lat(idxlist(ilat)%idx(ic), idxlist(ilat)%blk(ic), l) * idxlist(ilat)%wgt(ic)
          END DO
        END DO
      END DO

    END DO ! level

    ! compute point-wise sum over all mpi ranks and store results
    mocs_2d = p_sum(mocs_2d,mpi_comm)
    mocs_3d = p_sum(mocs_3d,mpi_comm)

    global_moc(1:n_zlev,:) = mocs_3d(:,:,1)
    atlant_moc(1:n_zlev,:) = mocs_3d(:,:,2)
    pacind_moc(1:n_zlev,:) = mocs_3d(:,:,3)

    global_hfl(1,:) = mocs_2d(:,1)
    atlant_hfl(1,:) = mocs_2d(:,2)
    pacind_hfl(1,:) = mocs_2d(:,3)
    global_wfl(1,:) = mocs_2d(:,4)
    atlant_wfl(1,:) = mocs_2d(:,5)
    pacind_wfl(1,:) = mocs_2d(:,6)
    global_hfbasin(1,:) = mocs_2d(:,7)
    atlant_hfbasin(1,:) = mocs_2d(:,8)
    pacind_hfbasin(1,:) = mocs_2d(:,9)
    global_sltbasin(1,:) = mocs_2d(:,10)
    atlant_sltbasin(1,:) = mocs_2d(:,11)
    pacind_sltbasin(1,:) = mocs_2d(:,12)

    ! compute partial sums along meridian
    DO l=nlat_moc-1,1,-1   ! fixed to 1 deg meridional resolution
      global_moc(:,l)=global_moc(:,l+1)+global_moc(:,l)
      atlant_moc(:,l)=atlant_moc(:,l+1)+atlant_moc(:,l)
      pacind_moc(:,l)=pacind_moc(:,l+1)+pacind_moc(:,l)
      global_hfl(:,l)=global_hfl(:,l+1)+global_hfl(:,l)
      atlant_hfl(:,l)=atlant_hfl(:,l+1)+atlant_hfl(:,l)
      pacind_hfl(:,l)=pacind_hfl(:,l+1)+pacind_hfl(:,l)
      global_hfbasin(:,l)=global_hfbasin(:,l+1)+global_hfbasin(:,l)
      atlant_hfbasin(:,l)=atlant_hfbasin(:,l+1)+atlant_hfbasin(:,l)
      pacind_hfbasin(:,l)=pacind_hfbasin(:,l+1)+pacind_hfbasin(:,l)
      global_wfl(:,l)=global_wfl(:,l+1)+global_wfl(:,l)
      atlant_wfl(:,l)=atlant_wfl(:,l+1)+atlant_wfl(:,l)
      pacind_wfl(:,l)=pacind_wfl(:,l+1)+pacind_wfl(:,l)
      global_sltbasin(:,l)=global_sltbasin(:,l+1)+global_sltbasin(:,l)
      atlant_sltbasin(:,l)=atlant_sltbasin(:,l+1)+atlant_sltbasin(:,l)
      pacind_sltbasin(:,l)=pacind_sltbasin(:,l+1)+pacind_sltbasin(:,l)
    END DO

    !find atlantic moc at 26n , depth=1000m
    factor_to_sv=1.0_wp/OceanReferenceDensity*1e-6_wp
    amoc26n(1)=atlant_moc(get_level_index_by_depth(patch_3d, 1000.0_wp),116)*factor_to_sv


    ! calculate ocean heat transport as residual from the tendency in heat content (dH/dt)
    ! minus the integral of surface heat flux

    global_hfbasin(:,:)=global_hfl(:,:)-global_hfbasin(:,:)
    atlant_hfbasin(:,:)=atlant_hfl(:,:)-atlant_hfbasin(:,:)
    pacind_hfbasin(:,:)=pacind_hfl(:,:)-pacind_hfbasin(:,:)

  END SUBROUTINE calc_moc_hfl_internal_vector
  !-------------------------------------------------------------------------


  SUBROUTINE calc_moc_hfl_internal_scalar (patch_2d, patch_3d, w, heatflux_total, &
             frshflux_volumetotal, delta_thetao, delta_so, delta_snow, delta_ice, &
             global_moc, atlant_moc, pacind_moc, global_hfl, atlant_hfl, pacind_hfl, &
             global_wfl, atlant_wfl, pacind_wfl, &
             global_hfbasin, atlant_hfbasin, pacind_hfbasin, &
             global_sltbasin, atlant_sltbasin, pacind_sltbasin, amoc26n, lacc)

    TYPE(t_patch),    TARGET, INTENT(in)  :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3d

    REAL(wp), INTENT(in)  :: w(:,:,:)   ! vertical velocity (nproma,nlev+1,alloc_cell_blocks)
    REAL(wp), INTENT(in)  :: heatflux_total(:,:)   ! heatflux_total (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(in)  :: frshflux_volumetotal(:,:)   ! fw flux (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(inout)  :: delta_snow(:,:)       ! tendency of snow thickness (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(inout)  :: delta_ice(:,:)        ! tendendy of ice  thickness (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(inout)  :: delta_thetao(:,:,:)   ! temperature tendency (nproma,nlev+1,alloc_cell_blocks)
    REAL(wp), INTENT(inout)  :: delta_so(:,:,:)   ! salinity tendency (nproma,nlev+1,alloc_cell_blocks)
    REAL(wp), INTENT(inout)  :: amoc26n(:)

    REAL(wp), INTENT(inout) :: global_moc(:,:), atlant_moc(:,:), pacind_moc(:,:) ! (n_zlev,nlat_moc)

    ! implied ocean heat transport calculated from surface fluxes
    REAL(wp), INTENT(inout) :: global_hfl(:,:), atlant_hfl(:,:), pacind_hfl(:,:) ! (1,nlat_moc)

    ! implied ocean fw transport calculated from surface fluxes
    REAL(wp), INTENT(inout) :: global_wfl(:,:), atlant_wfl(:,:), pacind_wfl(:,:) ! (1,nlat_moc)

    ! northward ocean heat transport calculated from tendencies
    REAL(wp), INTENT(inout) :: global_hfbasin(:,:), atlant_hfbasin(:,:), pacind_hfbasin(:,:) ! (1,nlat_moc)

    ! northward ocean salt transport calculated from tendencies
    REAL(wp), INTENT(inout) :: global_sltbasin(:,:), atlant_sltbasin(:,:), pacind_sltbasin(:,:) ! (1,nlat_moc)
    LOGICAL, INTENT(in), OPTIONAL :: lacc

    ! local variables
    INTEGER, PARAMETER ::  latSmooth = 3   !  latitudinal smoothing area is 2*jbrei-1 rows of 1 deg
    INTEGER :: BLOCK, level, start_index, end_index, idx, ilat, l, n
    INTEGER :: mpi_comm

    REAL(wp) :: lat, deltaMoc, deltahfl, deltawfl, deltahfbasin, deltasltbasin, smoothWeight
    REAL(wp), ALLOCATABLE :: allmocs(:,:,:)

    REAL(wp) :: factor_to_sv

    TYPE(t_subset_range), POINTER :: cells

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_ocean_diagnostics:calc_moc'
    LOGICAL  :: lzacc

    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    mpi_comm = MERGE(p_comm_work_test, p_comm_work, p_test_run)

    n=MAX(12,n_zlev) !needs at leat 12 levels to store the wfl/hfl/hfbasin variables
    ALLOCATE(allmocs(4,n,nlat_moc))

    !$ACC DATA CREATE(allmocs) IF(lzacc)

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    allmocs(:,:,:)  = 0.0_wp

    global_moc(:,:) = 0.0_wp
    pacind_moc(:,:) = 0.0_wp
    atlant_moc(:,:) = 0.0_wp

    global_hfl(:,:) = 0.0_wp
    pacind_hfl(:,:) = 0.0_wp
    atlant_hfl(:,:) = 0.0_wp

    global_wfl(:,:) = 0.0_wp
    pacind_wfl(:,:) = 0.0_wp
    atlant_wfl(:,:) = 0.0_wp

    global_hfbasin(:,:) = 0.0_wp
    pacind_hfbasin(:,:) = 0.0_wp
    atlant_hfbasin(:,:) = 0.0_wp

    global_sltbasin(:,:) = 0.0_wp
    pacind_sltbasin(:,:) = 0.0_wp
    atlant_sltbasin(:,:) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)

    ! limit cells to in-domain because of summation
    cells   => patch_2d%cells%in_domain

    smoothWeight = 1.0_wp / REAL(2*latSmooth + 1, wp)

    DO BLOCK = cells%start_block, cells%end_block
      CALL get_index_range(cells, BLOCK, start_index, end_index)
      ! 2023-08 psam-DKRZ: An alternate implementation with ACC LOOP GANG needs the use of atomic update,
      ! but it does not give bit-identical results compared to the CPU results
      !$ACC PARALLEL DEFAULT(PRESENT) &
      !$ACC   PRIVATE(deltaMoc, deltahfbasin, deltasltbasin, deltahfl, deltawfl, ilat) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO idx = start_index, end_index
        lat = patch_2d%cells%center(idx,BLOCK)%lat*rad2deg
        !$ACC LOOP SEQ
        DO level = 1, cells%vertical_levels(idx,BLOCK)

          deltaMoc = patch_2d%cells%area(idx,BLOCK) * OceanReferenceDensity * w(idx,level,BLOCK)

          deltahfbasin = patch_2d%cells%area(idx,BLOCK) * delta_thetao(idx,level,BLOCK)

          deltasltbasin = patch_2d%cells%area(idx,BLOCK) * delta_so(idx,level,BLOCK)

          IF (level .EQ. 1) THEN
            deltahfl = patch_2d%cells%area(idx,BLOCK) * heatflux_total(idx,BLOCK) &
                    * patch_3D%wet_c(idx,1,BLOCK)

            deltawfl = patch_2d%cells%area(idx,BLOCK) * frshflux_volumetotal(idx,BLOCK) &
                     * patch_3D%wet_c(idx,1,BLOCK)

            deltahfbasin = deltahfbasin                                &
                 + patch_2d%cells%area(idx,BLOCK) * ( delta_ice(idx,BLOCK) + delta_snow(idx,BLOCK) )

          ENDIF

          ! lat: corresponding latitude row of 1 deg extension
          !            1 south pole
          ! nlat_moc=180 north pole
          ilat     = NINT(REAL(nlat_moc, wp)*0.5_wp + lat)
          ilat     = MAX(1,MIN(ilat,nlat_moc))

          ! distribute MOC over (2*jbrei)+1 latitude rows
          !  - no weighting with latitudes done
          !  - lat: index of 180 X 1 deg meridional resolution
          !$ACC LOOP SEQ
          DO l = -latSmooth, latSmooth
            ilat = NINT(REAL(nlat_moc, wp)*0.5_wp + lat + REAL(l, wp))
            ilat = MAX(1,MIN(ilat,nlat_moc))

            global_moc(level,ilat) =       global_moc(level,ilat) - deltaMoc*smoothWeight
            atlant_moc(level,ilat) = atlant_moc(level,ilat) - MERGE(deltaMoc*smoothWeight, &
                 0.0_wp, patch_3D%basin_c(idx,BLOCK) == 1)
            pacind_moc(level,ilat) = pacind_moc(level,ilat) - MERGE(deltaMoc*smoothWeight, &
                 0.0_wp, patch_3D%basin_c(idx,BLOCK) >= 2)

            global_hfbasin(1,ilat) =       global_hfbasin(1,ilat) - deltahfbasin*smoothWeight
            atlant_hfbasin(1,ilat) = atlant_hfbasin(1,ilat) - MERGE(deltahfbasin*smoothWeight, &
               0.0_wp, patch_3D%basin_c(idx,BLOCK) == 1)
            pacind_hfbasin(1,ilat) = pacind_hfbasin(1,ilat) - MERGE(deltahfbasin*smoothWeight, &
               0.0_wp, patch_3D%basin_c(idx,BLOCK) >= 2)

            global_sltbasin(1,ilat) =       global_sltbasin(1,ilat) - deltasltbasin*smoothWeight
            atlant_sltbasin(1,ilat) = atlant_sltbasin(1,ilat) - MERGE(deltasltbasin*smoothWeight, &
               0.0_wp, patch_3D%basin_c(idx,BLOCK) == 1)
            pacind_sltbasin(1,ilat) = pacind_sltbasin(1,ilat) - MERGE(deltasltbasin*smoothWeight, &
               0.0_wp, patch_3D%basin_c(idx,BLOCK) >= 2)

            IF (level .EQ. 1) THEN
              global_hfl(level,ilat) =       global_hfl(level,ilat) - deltahfl*smoothWeight
              atlant_hfl(level,ilat) = atlant_hfl(level,ilat) - MERGE(deltahfl*smoothWeight, &
                   0.0_wp, patch_3D%basin_c(idx,BLOCK) == 1)
              pacind_hfl(level,ilat) = pacind_hfl(level,ilat) - MERGE(deltahfl*smoothWeight, &
                   0.0_wp, patch_3D%basin_c(idx,BLOCK) >= 2)

              global_wfl(level,ilat) =       global_wfl(level,ilat) - deltawfl*smoothWeight
              atlant_wfl(level,ilat) = atlant_wfl(level,ilat) - MERGE(deltawfl*smoothWeight, &
                   0.0_wp, patch_3D%basin_c(idx,BLOCK) == 1)
              pacind_wfl(level,ilat) = pacind_wfl(level,ilat) - MERGE(deltawfl*smoothWeight, &
                   0.0_wp, patch_3D%basin_c(idx,BLOCK) >= 2)
            END IF

          END DO

        END DO
      END DO
      !$ACC END PARALLEL
    END DO
    !$ACC WAIT(1)

    ! compute point-wise sum over all mpi ranks and store results
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    allmocs(1,1:n_zlev,:) = global_moc(1:n_zlev,:)
    allmocs(2,1:n_zlev,:) = atlant_moc(1:n_zlev,:)
    allmocs(3,1:n_zlev,:) = pacind_moc(1:n_zlev,:)

    allmocs(4,1,:) = global_hfl(1,:)
    allmocs(4,2,:) = atlant_hfl(1,:)
    allmocs(4,3,:) = pacind_hfl(1,:)

    allmocs(4,4,:) = global_hfbasin(1,:)
    allmocs(4,5,:) = atlant_hfbasin(1,:)
    allmocs(4,6,:) = pacind_hfbasin(1,:)

    allmocs(4,7,:) = global_wfl(1,:)
    allmocs(4,8,:) = atlant_wfl(1,:)
    allmocs(4,9,:) = pacind_wfl(1,:)

    allmocs(4,10,:) = global_sltbasin(1,:)
    allmocs(4,11,:) = atlant_sltbasin(1,:)
    allmocs(4,12,:) = pacind_sltbasin(1,:)
    !$ACC END KERNELS
    !$ACC WAIT(1)

    !$ACC UPDATE HOST(allmocs) ASYNC(1) IF(lzacc)
    !$ACC WAIT(1) IF(lzacc)
    allmocs = p_sum(allmocs,mpi_comm)
    !$ACC UPDATE DEVICE(allmocs) ASYNC(1) IF(lzacc)
    !$ACC WAIT(1) IF(lzacc) ! can be removed when all ACC compute regions are ASYNC(1)


    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    global_moc(1:n_zlev,:) = allmocs(1,1:n_zlev,:)
    atlant_moc(1:n_zlev,:) = allmocs(2,1:n_zlev,:)
    pacind_moc(1:n_zlev,:) = allmocs(3,1:n_zlev,:)
    global_hfl(1,:) = allmocs(4,1,:)
    atlant_hfl(1,:) = allmocs(4,2,:)
    pacind_hfl(1,:) = allmocs(4,3,:)
    global_hfbasin(1,:) = allmocs(4,4,:)
    atlant_hfbasin(1,:) = allmocs(4,5,:)
    pacind_hfbasin(1,:) = allmocs(4,6,:)
    global_wfl(1,:) = allmocs(4,7,:)
    atlant_wfl(1,:) = allmocs(4,8,:)
    pacind_wfl(1,:) = allmocs(4,9,:)
    global_sltbasin(1,:) = allmocs(4,10,:)
    atlant_sltbasin(1,:) = allmocs(4,11,:)
    pacind_sltbasin(1,:) = allmocs(4,12,:)
    !$ACC END KERNELS
    !$ACC WAIT(1)

    ! compute partial sums along meridian
    DO l=nlat_moc-1,1,-1   ! fixed to 1 deg meridional resolution
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      global_moc(:,l)=global_moc(:,l+1)+global_moc(:,l)
      atlant_moc(:,l)=atlant_moc(:,l+1)+atlant_moc(:,l)
      pacind_moc(:,l)=pacind_moc(:,l+1)+pacind_moc(:,l)
      global_hfl(:,l)=global_hfl(:,l+1)+global_hfl(:,l)
      atlant_hfl(:,l)=atlant_hfl(:,l+1)+atlant_hfl(:,l)
      pacind_hfl(:,l)=pacind_hfl(:,l+1)+pacind_hfl(:,l)
      global_hfbasin(:,l)=global_hfbasin(:,l+1)+global_hfbasin(:,l)
      atlant_hfbasin(:,l)=atlant_hfbasin(:,l+1)+atlant_hfbasin(:,l)
      pacind_hfbasin(:,l)=pacind_hfbasin(:,l+1)+pacind_hfbasin(:,l)
      global_wfl(:,l)=global_wfl(:,l+1)+global_wfl(:,l)
      atlant_wfl(:,l)=atlant_wfl(:,l+1)+atlant_wfl(:,l)
      pacind_wfl(:,l)=pacind_wfl(:,l+1)+pacind_wfl(:,l)
      global_sltbasin(:,l)=global_sltbasin(:,l+1)+global_sltbasin(:,l)
      atlant_sltbasin(:,l)=atlant_sltbasin(:,l+1)+atlant_sltbasin(:,l)
      pacind_sltbasin(:,l)=pacind_sltbasin(:,l+1)+pacind_sltbasin(:,l)
      !$ACC END KERNELS
      !$ACC WAIT(1)
    END DO

    !$ACC UPDATE HOST(atlant_moc) ASYNC(1) IF(lzacc)
    !$ACC WAIT(1) IF(lzacc)
    !find atlantic moc at 26n , depth=1000m
    factor_to_sv=1.0_wp/OceanReferenceDensity*1e-6_wp
    amoc26n(1)=atlant_moc(get_level_index_by_depth(patch_3d, 1000.0_wp),116)*factor_to_sv
    !$ACC UPDATE DEVICE(amoc26n) ASYNC(1) IF(lzacc)
    !$ACC WAIT(1) IF(lzacc) ! can be removed when all ACC compute regions are ASYNC(1)

    ! calculate ocean heat transport as residual from the tendency in heat content (dH/dt)
    ! minus the integral of surface heat flux

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    global_hfbasin(:,:)=global_hfl(:,:)-global_hfbasin(:,:)
    atlant_hfbasin(:,:)=atlant_hfl(:,:)-atlant_hfbasin(:,:)
    pacind_hfbasin(:,:)=pacind_hfl(:,:)-pacind_hfbasin(:,:)
    !$ACC END KERNELS
    !$ACC WAIT(1)

    !$ACC END DATA
    DEALLOCATE (allmocs)

  END SUBROUTINE calc_moc_hfl_internal_scalar
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !
  !
  !!  Calculation of horizontal stream function
  !
  !>
  !!
  !!  based on code from MPIOM
  !
  ! TODO: implement variable output dimension (1 deg resolution) and smoothing extent
  !!
!<Optimize:inUse>
  SUBROUTINE calc_psi (patch_3D, u, prism_thickness, u_vint, this_datetime, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3D
    REAL(wp), INTENT(in)               :: u(:,:,:)     ! zonal velocity at cell centers
    REAL(wp), INTENT(in)               :: prism_thickness(:,:,:)       ! elevation on cell centers
    ! dims: (nproma,nlev,alloc_cell_blocks)
    REAL(wp), INTENT(inout)            :: u_vint(:,:)  ! barotropic zonal velocity on icon grid
    TYPE(datetime), POINTER            :: this_datetime
    !
    ! local variables
    ! INTEGER :: i

    ! switch for writing stream function (not yet in namelist); 1: icon-grid; 2: regular grid output
    INTEGER, PARAMETER ::  idiag_psi = 1

    INTEGER, PARAMETER ::  nlat = 180                    ! meridional dimension of regular grid
    INTEGER, PARAMETER ::  nlon = 360                    ! zonal dimension of regular grid

    ! smoothing area is 2*jsmth-1 lat/lon areas of 1 deg
    INTEGER, PARAMETER ::  jsmth = 3
    INTEGER :: blockNo, jc, jk, start_index, end_index
    INTEGER :: jlat, jlon, jlt, jln, jltx, jlnx, jsmth2
    INTEGER(i8)        :: idate, iextra(4)

    REAL(wp) :: z_lat_deg, z_lon_deg, z_lat_dist, delta_z, rsmth
    REAL(wp), ALLOCATABLE :: z_uint_reg(:,:)                     ! vertical integral on regular grid
    REAL(wp), ALLOCATABLE :: psi_reg(:,:)                        ! horizontal stream function

    TYPE(t_patch), POINTER  :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells, dom_cells

    LOGICAL, INTENT(in), OPTIONAL :: lacc
    LOGICAL :: lzacc

    !CHARACTER(len=*), PARAMETER :: routine = 'mo_ocean_diagnostics:calc_psi'

    !-----------------------------------------------------------------------

    jsmth2          = 2*jsmth + 1
    rsmth           = REAL(jsmth2*jsmth2, wp)

    ! with all cells no sync is necessary
    patch_2d  => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    dom_cells => patch_2d%cells%in_domain

    CALL set_acc_host_or_device(lzacc, lacc)

    ! (1) barotropic system:
    !     vertical integration of zonal velocity times vertical layer thickness [m/s*m]

!ICON_OMP_PARALLEL_DO PRIVATE(jc, jk, start_index, end_index) SCHEDULE(dynamic)
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      u_vint(:,blockNo)     = 0.0_wp
      !$ACC END KERNELS
      !$ACC WAIT(1)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index

        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)

          u_vint(jc,blockNo) = u_vint(jc,blockNo) - u(jc,jk,blockNo) * prism_thickness(jc,jk,blockNo)

        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('calc_psi:    u_vint    ', u_vint        , str_module, idt_src, in_subset=patch_2d%cells%owned)
    !---------------------------------------------------------------------

    IF (idiag_psi == 1) RETURN

    ALLOCATE(z_uint_reg(nlon,nlat), psi_reg(nlon,nlat))
    !$OMP PARALLEL WORKSHARE
    z_uint_reg(:,:) = 0.0_wp
    !$OMP END PARALLEL WORKSHARE
    ! (2) distribute integrated zonal velocity (u*dz) on 1x1 deg grid
    !     this code is not mature yet

    ! in domain: count all cells only once
    !$OMP PARALLEL DO PRIVATE(start_index,end_index,jc,z_lat_deg,z_lon_deg, &
    !$OMP & jlat,jlon,jltx,jlt,jlnx,jln)
    DO blockNo = dom_cells%start_block, dom_cells%end_block
      CALL get_index_range(dom_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        z_lat_deg = patch_2d%cells%center(jc,blockNo)%lat * rad2deg
        z_lon_deg = patch_2d%cells%center(jc,blockNo)%lon * rad2deg

        !  ! 0 <= lon <= 360 deg
        !  z_lon_deg = z_lon_deg + 180.0_wp

        ! jlat/jlon: corresponding latitude/longitude coordinates of 1 deg extension
        ! jlat: 1 = south of 89.0S; 89 = 1S-Eq.; 90 = Eq-1N;  180 = north of 89N
        ! jlon: 1 = 180W-179W; 180 = 1-0 deg west; 360 = 179E-180E

        jlat = NINT(91.0_wp + z_lat_deg)
        jlon = NINT(z_lon_deg + 180.5_wp)

        ! distribute stream function over rsmth=(2*jsmth+1)**2 lat/lon regular grid points
        !  - no weighting with latitudes done
        !  - no correction with regular lsm done
        DO jltx = jlat-jsmth, jlat+jsmth

          jlt = jltx
          IF (jlt <    1) jlt =      1-jlt  ! apply equatorwards
          IF (jlt > nlat) jlt = 2*nlat-jlt  ! apply equatorwards
          DO jlnx = jlon-jsmth, jlon+jsmth

            jln = jlnx
            IF (jln <    1) jln = jln+nlon  ! circular boundary
            IF (jln > nlon) jln = jln-nlon  ! circular boundary
            !$OMP ATOMIC
            z_uint_reg(jln,jlt) = z_uint_reg(jln,jlt) + u_vint(jc,blockNo) / rsmth


          END DO
        END DO

      END DO
    END DO
    !$OMP END PARALLEL DO

    ! (3) calculate meridional integral on regular grid starting from south pole:

    DO jlt = nlat-1, 1, -1
      z_uint_reg(:,jlt) = z_uint_reg(:,jlt) + z_uint_reg(:,jlt+1)
    END DO

    ! (4) calculate stream function: scale with length of 1 deg*rho [m/s*m*m*kg/m3=kg/s]

    ! meridional distance of 1 deg
    ! ATTENTION - fixed 1 deg resolution should be related to icon-resolution
    z_lat_dist = 111111.0_wp  ! * 1.3_wp ??
    !$OMP PARALLEL WORKSHARE
    psi_reg(:,:) = z_uint_reg(:,:) * z_lat_dist * OceanReferenceDensity
    !$OMP END PARALLEL WORKSHARE
    ! stream function on icon grid without calculation of meridional integral
    !  - tbd after interpolation to regular grid externally
    !  psi    (:,:) = u_vint    (:,:)              * OceanReferenceDensity


    ! write out in extra format - integer*8
    idate = INT(this_datetime%date%month*1000000+this_datetime%date%day*10000 &
         &     +this_datetime%time%hour*100+this_datetime%time%minute,i8)
    WRITE(0,*) 'write global PSI at iyear, idate:',this_datetime%date%year, idate

    iextra(1) = INT(idate,i8)
    iextra(2) = INT(780,i8)
    iextra(3) = INT(0,i8)
    iextra(4) = INT(nlon*nlat,i8)

    WRITE(80) (iextra(blockNo),blockNo=1,4)
    WRITE(80) ((psi_reg(jln,jlt),jln=1,nlon),jlt=1,nlat)

    DO jlat=1,nlat
      WRITE(82,*) 'jlat=',jlat
      WRITE(82,'(1p10e12.3)') (psi_reg(jlon,jlat),jlon=1,nlon)
    ENDDO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
 !  idt_src=3  ! output print level (1-5, fix)
 !  CALL dbg_print('calc_psi:    psi_reg   ', psi_reg       , str_module, idt_src, in_subset=patch_2d%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE calc_psi
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !
  !!  Calculation of horizontal stream function using normal velocity on edges
  !
  !>
  !>  Calculation of horizontal stream function using normal velocity on edges
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE calc_psi_vn (patch_3D, vn, prism_thickness_e, op_coeff, u_vint, v_vint, this_datetime)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3D
    REAL(wp), INTENT(in)               :: vn(:,:,:)                 ! normal velocity at cell edges
    REAL(wp), INTENT(in)               :: prism_thickness_e(:,:,:)  ! elevation on cell edges
    ! dims: (nproma,nlev,alloc_edge_blocks)
    REAL(wp), INTENT(inout)            :: u_vint(:,:)               ! barotropic zonal velocity on cell centers
    REAL(wp), INTENT(inout)            :: v_vint(:,:)               ! barotropic meridional velocity on cell centers
    TYPE(t_operator_coeff),INTENT(in)  :: op_coeff
    TYPE(datetime), POINTER            :: this_datetime
    !
    INTEGER  :: blockNo, jc, je, jk, start_index, end_index
    ! vertical integral vn on edges and in cartesian coordinates - no 2-dim mapping is available
    REAL(wp) :: vn_vint(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: u_2d(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)   !  scratch arrays for test
    REAL(wp) :: v_2d(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)   !  scratch arrays for test
    TYPE(t_cartesian_coordinates) :: vint_cc(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    TYPE(t_patch), POINTER  :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges, all_cells

    !CHARACTER(len=*), PARAMETER :: routine = 'mo_ocean_diagnostics:calc_psi_vn'

    !-----------------------------------------------------------------------

    patch_2d  => patch_3d%p_patch_2d(1)
    all_edges => patch_2d%edges%ALL
    all_cells => patch_2d%cells%ALL
    !$OMP PARALLEL WORKSHARE
    vn_vint  (:,:,:)    = 0.0_wp
    vint_cc(:,:,:)%x(1) = 0.0_wp
    vint_cc(:,:,:)%x(2) = 0.0_wp
    vint_cc(:,:,:)%x(3) = 0.0_wp
    u_2d       (:,:)    = 0.0_wp
    v_2d       (:,:)    = 0.0_wp
    !$OMP END PARALLEL WORKSHARE
    ! (1) barotropic system:
    !     vertical integration of normal velocity times vertical layer thickness [m/s*m]
!ICON_OMP_PARALLEL_DO PRIVATE(je, jk, start_index, end_index) SCHEDULE(dynamic)
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_index, end_index)
      vn_vint(:,1,blockNo) = 0.0_wp
      DO je = start_index, end_index

        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

          vn_vint(je,1,blockNo) = vn_vint(je,1,blockNo) - vn(je,jk,blockNo) * prism_thickness_e(je,jk,blockNo)

        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

    ! (2) remapping normal velocity to zonal and meridional velocity at cell centers
    CALL sync_patch_array(sync_e, patch_2d, vn_vint)

    CALL map_edges2cell_3d(patch_3D, vn_vint, op_coeff, vint_cc)

    CALL sync_patch_array(sync_c, patch_2d, vint_cc(:,:,:)%x(1))
    CALL sync_patch_array(sync_c, patch_2d, vint_cc(:,:,:)%x(2))
    CALL sync_patch_array(sync_c, patch_2d, vint_cc(:,:,:)%x(3))

    ! calculate zonal and meridional velocity:
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        CALL cvec2gvec(vint_cc(jc,1,blockNo)%x(1), vint_cc(jc,1,blockNo)%x(2), vint_cc(jc,1,blockNo)%x(3), &
          &            patch_2d%cells%center(jc,blockNo)%lon, patch_2d%cells%center(jc,blockNo)%lat,  &
          &            u_2d(jc,blockNo), v_2d(jc,blockNo))
!         &            u_vint(jc,blockNo), v_vint(jc,blockNo))
      END DO
    END DO
    !CALL sync_patch_array(sync_c, patch_2d, u_vint)
    !CALL sync_patch_array(sync_c, patch_2d, v_vint)
    CALL sync_patch_array(sync_c, patch_2d, u_2d)
    CALL sync_patch_array(sync_c, patch_2d, v_2d)

    ! hack for test: calc_psy for u_vint, calc_psi_vn for v_vint - accumulated and written out
    v_vint(:,:) = u_2d(:,:)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('calc_psi_vn: u_2d        ', u_2d        , str_module, idt_src, in_subset=patch_2d%cells%owned)
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('calc_psi_vn: v_2d        ', v_2d        , str_module, idt_src, in_subset=patch_2d%cells%owned)
    CALL dbg_print('calc_psi_vn: vint_cc%x(1)' ,vint_cc%x(1), str_module, idt_src, in_subset=patch_2d%cells%owned)
    !---------------------------------------------------------------------


  END SUBROUTINE calc_psi_vn
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!taken from MPIOM
!<Optimize:inUse>
  FUNCTION calc_max_condep(vertical_density_gradient,max_lev) RESULT(condep)
    !$ACC ROUTINE SEQ
    REAL(wp),INTENT(in)  :: vertical_density_gradient(n_zlev)
    INTEGER, INTENT(in)  :: max_lev
    INTEGER :: condep

    INTEGER :: jk
    INTEGER :: maxcondep     !< maximum convective penetration level
    REAL(wp) :: masked_vertical_density_gradient(n_zlev)

    condep = 1

    ! remove dbl_eps, which  is added in the vertical gradient computation
    masked_vertical_density_gradient = MAX(vertical_density_gradient - dbl_eps,0.0_wp)

    !! diagnose maximum convection level
    !! condep = maximum model level penetrated by vertically continous
    !! convection from the surface downward
    !! calculated over integration period ; it should be written out
    !! as snapshot at the end of the run
   maxcondep=1
    DO jk=2,max_lev
      IF (masked_vertical_density_gradient(jk) .ne. 0.0_wp) THEN
        maxcondep = jk
        EXIT
      ENDIF
    ENDDO

    condep = max0(maxcondep,condep)
  END FUNCTION calc_max_condep
!<Optimize:inUse>
  FUNCTION calc_mixed_layer_depth(vertical_density_gradient,critical_value,min_lev,max_lev,thickness, depth_of_first_layer) &
    & result(mixed_layer_depth)
    !$ACC ROUTINE SEQ
#if defined(__NVCOMPILER_MAJOR__) && __NVCOMPILER_MAJOR__ <= 21
    REAL(wp), INTENT(in)  :: vertical_density_gradient(n_zlev)
#else
    REAL(wp), INTENT(in)  :: vertical_density_gradient(:)
#endif
    REAL(wp), INTENT(in)  :: critical_value
    INTEGER,  INTENT(in)  :: min_lev
    INTEGER,  INTENT(in)  :: max_lev
#if defined(__NVCOMPILER_MAJOR__) && __NVCOMPILER_MAJOR__ <= 21
    REAL(wp), INTENT(in)  :: thickness(n_zlev)
#else
    REAL(wp), INTENT(in)  :: thickness(:)
#endif
    REAL(wp), INTENT(in)  :: depth_of_first_layer

    REAL(wp) :: sigh        ,zzz
    REAL(wp) :: mixed_layer_depth
    REAL(wp) :: masked_vertical_density_gradient
    INTEGER :: jk

    sigh              = critical_value
    mixed_layer_depth = depth_of_first_layer

    ! This diagnostic calculates the mixed layer depth.
    ! It uses the incremental density increase between two
    ! levels and substracts it from the initial density criterion (sigcrit)
    ! and adds the level thickness (zzz) to zmld. This is done till
    ! the accumulated density increase between the surface and
    ! layer k is sigcrit or sigh = O, respectively.

    ! stabio(k) = insitu density gradient
    ! sigh = remaining density difference

    DO jk = min_lev+1, max_lev
      IF (sigh .GT. 1.e-6_wp) THEN
        masked_vertical_density_gradient = MAX(vertical_density_gradient(jk),0.0_wp)
        zzz               = MIN(sigh / (masked_vertical_density_gradient + 1.0E-19_wp), thickness(jk))
        sigh              = sigh - zzz*masked_vertical_density_gradient
        mixed_layer_depth = mixed_layer_depth + zzz
      ELSE
        sigh = 0._wp
      ENDIF
    ENDDO

  END FUNCTION calc_mixed_layer_depth


  SUBROUTINE diag_heat_salt_tendency(patch_3d, n, ice, thetao, so, delta_ice, delta_snow, &
       delta_thetao, delta_so, stretch_c, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(in)     :: patch_3D

    REAL(wp), INTENT(in)                     :: thetao(:,:,:)   ! temperature
    REAL(wp), INTENT(in)                     :: so(:,:,:)   ! salinity
    REAL(wp), INTENT(in),OPTIONAL            :: stretch_c(:, :)
    TYPE(t_sea_ice), INTENT(in)              :: ice
    TYPE(t_subset_range), POINTER            :: subset

    REAL(wp), INTENT(inout)  :: delta_ice(:,:)
    REAL(wp), INTENT(inout)  :: delta_snow(:,:)
    REAL(wp), INTENT(inout)  :: delta_thetao(:,:,:)
    REAL(wp), INTENT(inout)  :: delta_so(:,:,:)

    LOGICAL, INTENT(in), OPTIONAL :: lacc

    INTEGER  :: n, blk, cell, cellStart,cellEnd, level, dz
    REAL(wp) :: dti, rhoicwa, rhosnic, rhosnwa, tfreeze, tmelt,           &
                tref, entmel,  sithk, snthk
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    rhoicwa = rhoi / rho_ref
    rhosnwa = rhos / rho_ref
    rhosnic = rhos / rhoi

    tfreeze = -1.9
    tmelt = 273.15
    tref = 273.15
    entmel = rhoi * alf

    subset => patch_3d%p_patch_2d(1)%cells%owned

    IF ( n .EQ. 1) THEN

      DO blk = subset%start_block, subset%end_block
        CALL get_index_range(subset, blk, cellStart, cellEnd)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO cell = cellStart, cellEnd

          delta_ice(cell,blk) = SUM(ice%hi(cell,:,blk)*ice%conc(cell,:,blk))
          delta_snow(cell,blk) = SUM(ice%hs(cell,:,blk)*ice%conc(cell,:,blk))

          !$ACC LOOP SEQ
          DO level = 1,subset%vertical_levels(cell,blk)
            delta_thetao(cell,level,blk) = thetao(cell,level,blk)
            delta_so(cell,level,blk) = so(cell,level,blk)
          END DO ! level

        END DO ! cell
        !$ACC END PARALLEL LOOP
      END DO ! blk
      !$ACC WAIT(1)

    ENDIF

    IF ( n .EQ. 2) THEN

      dti = 1.0_wp / dtime

      DO blk = subset%start_block, subset%end_block
        CALL get_index_range(subset, blk, cellStart, cellEnd)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(sithk, snthk, dz) ASYNC(1) IF(lzacc)
        DO cell = cellStart, cellEnd

          ! tendency of equivalent thickness of sea ice
          sithk = SUM(ice%hi(cell,:,blk)*ice%conc(cell,:,blk)) - delta_ice(cell,blk)

          ! converted to heat content
          delta_ice(cell,blk) = (( rhoicwa * clw * OceanReferenceDensity * sithk           &
               * ( tfreeze + tmelt - tref )  )                                            &
               - ( sithk * entmel )) * dti

          ! tendency of equivalent thickness of snow
          snthk = SUM(ice%hs(cell,:,blk)*ice%conc(cell,:,blk)) - delta_snow(cell,blk)

          ! converted to heat content
          delta_snow(cell,blk) = (( rhosnwa * clw * OceanReferenceDensity * snthk          &
               * ( tmelt - tref )  )                                                      &
               - ( rhosnic * snthk * entmel )) * dti

          !$ACC LOOP SEQ
          DO level = 1,subset%vertical_levels(cell,blk)

            IF (vert_cor_type .EQ. 1) THEN
              dz=MERGE(ice%zunderice(cell,blk)                     &
                   ,patch_3D%p_patch_1d(1)%prism_thick_c(cell,level,blk)*stretch_c(cell,blk),level.EQ.1)
            ELSE
              dz=MERGE(ice%zunderice(cell,blk), &
                   patch_3D%p_patch_1d(1)%prism_thick_c(cell,level,blk),level.EQ.1)
            ENDIF

            delta_thetao(cell,level,blk) = ( thetao(cell,level,blk) - delta_thetao(cell,level,blk) ) &
                 * clw * OceanReferenceDensity * dz * dti
            delta_so(cell,level,blk) = ( so(cell,level,blk) - delta_so(cell,level,blk) ) &
                 * dz * dti

          END DO ! level

        END DO ! cell
        !$ACC END PARALLEL LOOP
      END DO ! blk
      !$ACC WAIT(1)

    ENDIF
    !$ACC UPDATE HOST(delta_ice, delta_snow, delta_so, delta_thetao) &
    !$ACC   ASYNC(1) IF(lzacc)
    !$ACC WAIT(1) IF(lzacc)

  END SUBROUTINE diag_heat_salt_tendency

  SUBROUTINE calc_heat_content(patch_3d, thickness, ice, tracers, &
       heat_content_liquid_water, heat_content_seaice,            &
       heat_content_snow, heat_content_total, stretch_c, lacc)

    TYPE(t_patch_3d), TARGET, INTENT(in)  :: patch_3d

    REAL(wp), INTENT(IN)   :: thickness(:,:,:)
    REAL(wp), INTENT(IN)   :: tracers(:,:,:,:)
    REAL(wp), INTENT(INOUT)  :: heat_content_liquid_water(:,:,:)
    REAL(wp), INTENT(INOUT)  :: heat_content_seaice(:,:)
    REAL(wp), INTENT(INOUT)  :: heat_content_snow(:,:)
    REAL(wp), INTENT(INOUT)  :: heat_content_total(:,:)
    REAL(wp), INTENT(IN), OPTIONAL  :: stretch_c(:,:)

    TYPE(t_sea_ice), INTENT(IN)              :: ice
    TYPE(t_subset_range), POINTER            :: subset

    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER  :: blk, cell, cellStart,cellEnd, level
    REAL(wp) :: rhoicwa, rhosnic, rhosnwa, tfreeze, tmelt, &
         tref, entmel, rocp, sithk, snthk, dz
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    rhoicwa = rhoi / rho_ref
    rhosnwa = rhos / rho_ref
    rhosnic = rhos / rhoi
    rocp = rho_ref * clw
    tfreeze = -1.9
    tmelt = 273.15
    tref = 273.15
    entmel = rhoi * alf

    subset => patch_3d%p_patch_2d(1)%cells%owned
    DO blk = subset%start_block, subset%end_block
      CALL get_index_range(subset, blk, cellStart, cellEnd)
      !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) PRIVATE(sithk, snthk, dz) ASYNC(1) IF(lzacc) ! 2023-07 psam-DKRZ: Use of GANG VECTOR introduces error here
      DO cell = cellStart, cellEnd
        ! surface:
        ! heat of ice : heat of water equivalent at tfreeze - latent heat of fusion

        sithk = SUM(ice%hi(cell,:,blk)*ice%conc(cell,:,blk)) ! equivalent thickness of sea ice equally distributed over the cell area

        heat_content_seaice(cell,blk) = ( rhoicwa * rocp * sithk  &
             * ( tfreeze + tmelt - tref )  )                      &
             - ( sithk * entmel )

        ! heat of snow : heat of water equivalent at tmelt - latent heat of fusion

        snthk = SUM(ice%hs(cell,:,blk)*ice%conc(cell,:,blk)) ! equivalent thickness of snow on sea ice equally distributed over the cell area

        heat_content_snow(cell,blk) = ( rhosnwa * rocp * snthk  &
             * ( tmelt - tref )  )                              &
             - ( rhosnic * snthk * entmel )

        ! liquid water heat
        ! surface : tho * rho * cp * draft

        heat_content_liquid_water(cell,1,blk) = (tmelt - tref &
             + tracers(cell,1,blk,1) ) * rocp                 &
             * ice%zUnderIce(cell,blk)


        !$ACC LOOP SEQ
        DO level=2,subset%vertical_levels(cell,blk)

          ! 2023-07 dzo-DKRZ: The following MERGE command does not work as intended with NVIDIA compiler
          !dz=MERGE(thickness(cell,level,blk)*stretch_c(cell,blk) &
          !  ,thickness(cell,level,blk),vert_cor_type .EQ. 1)  !check for vert_cor_type
          IF (vert_cor_type .EQ. 1) THEN
             dz=thickness(cell,level,blk)*stretch_c(cell,blk)
          ELSE
             dz=thickness(cell,level,blk)
          END IF

          heat_content_liquid_water(cell,level,blk) = (tmelt - tref &
               + tracers(cell,level,blk,1) ) * rocp                 &
               * dz
        END DO

        ! total heat per column
        heat_content_total(cell,blk) = heat_content_snow(cell,blk) &
             + heat_content_seaice(cell, blk)                      &
             + SUM(heat_content_liquid_water(cell,1:subset%vertical_levels(cell,blk),blk))

        ! rest of the underwater world
      END DO ! cell
        !$ACC END PARALLEL LOOP
    END DO !block
    !$ACC WAIT(1)
    ! 2023-07 psam-DKRZ: The following UPDATE SELF directive is necessary as the updated arrays are required elsewhere
    ! for CPU-operations. This should not be necessary, I guess, when all subroutines are ported to GPU
    !$ACC UPDATE HOST(heat_content_liquid_water, heat_content_seaice, heat_content_snow, heat_content_total) &
    !$ACC   ASYNC(1) IF(lzacc)
    !$ACC WAIT(1) IF(lzacc)

  END SUBROUTINE calc_heat_content


  SUBROUTINE calc_eddydiag(patch_3d,u,v,w,w_prismcenter,T,S,R &
               ,uT, uS, uR, uu    &
               ,vT, vS, vR, vv    &
               ,wT, wS, wR, ww, uv, uw, vw, rr, ss, tt, sigma0   &
               ,hflr, fwr, tauxu, tauyv &
               , topbc_windstress_u, topbc_windstress_v &
               ,heatflux_total, frshflux_volumetotal, lacc)

    TYPE(t_patch_3d), TARGET, INTENT(in)  :: patch_3d

    REAL(wp), INTENT(IN)   :: u(:,:,:) !< zonal velocity at cell center
    REAL(wp), INTENT(IN)   :: v(:,:,:) !< meridional velocity at cell center
    REAL(wp), INTENT(IN)   :: w(:,:,:) !< vertical velocity at cell interface
    REAL(wp), INTENT(INOUT):: w_prismcenter(:,:,:) !< vertical velocity at cell center
    REAL(wp), INTENT(IN)   :: T(:,:,:) !< temerature
    REAL(wp), INTENT(IN)   :: S(:,:,:) !< salinity
    REAL(wp), INTENT(IN)   :: R(:,:,:) !< density
    REAL(wp), INTENT(in)  :: heatflux_total(:,:)   !< net heatflux
    REAL(wp), INTENT(in)  :: frshflux_volumetotal(:,:)   !< net fresh water flux
    REAL(wp), INTENT(in)  :: topbc_windstress_u(:,:)  !< windstress x
    REAL(wp), INTENT(in)  :: topbc_windstress_v(:,:)  !< windstress y


    REAL(wp), INTENT(INOUT)  :: sigma0(:,:,:) !< density - 1000

    REAL(wp), INTENT(INOUT)  :: hflR(:,:) !< product of netheatflux and density
    REAL(wp), INTENT(INOUT)  :: fwR(:,:) !< product of fw flux and density
    REAL(wp), INTENT(INOUT)  :: tauxU(:,:) !< product of x-windstress and u-velocity
    REAL(wp), INTENT(INOUT)  :: tauyV(:,:) !< product of y-windstress and v-velocity
    REAL(wp), INTENT(INOUT)  :: uT(:,:,:) !< product of temperature and u-velocity
    REAL(wp), INTENT(INOUT)  :: uS(:,:,:) !< product of salinity and u-velocity
    REAL(wp), INTENT(INOUT)  :: uR(:,:,:) !< product of density and u-velocity
    REAL(wp), INTENT(INOUT)  :: uu(:,:,:) !< square of u-velocity

    REAL(wp), INTENT(INOUT)  :: RR(:,:,:) !< square of density
    REAL(wp), INTENT(INOUT)  :: SS(:,:,:) !< square of salinity
    REAL(wp), INTENT(INOUT)  :: TT(:,:,:) !< square of temperature
    REAL(wp), INTENT(INOUT)  :: vT(:,:,:) !< product of temperature and v-velocity
    REAL(wp), INTENT(INOUT)  :: vS(:,:,:) !< product of salinity and v-velocity
    REAL(wp), INTENT(INOUT)  :: vR(:,:,:) !< product of density and v-velocity
    REAL(wp), INTENT(INOUT)  :: vv(:,:,:) !< square of  v-velocity

    REAL(wp), INTENT(INOUT)  :: wT(:,:,:) !< product of temperature and w-velocity
    REAL(wp), INTENT(INOUT)  :: wS(:,:,:) !< product of salinity and w-velocity
    REAL(wp), INTENT(INOUT)  :: wR(:,:,:) !< product of density and w-velocity
    REAL(wp), INTENT(INOUT)  :: ww(:,:,:) !< square of w-velocity

    REAL(wp), INTENT(INOUT)  :: uv(:,:,:) !< product of u-velocity and w-velocity
    REAL(wp), INTENT(INOUT)  :: uw(:,:,:) !< product of v-velocity and w-velocity
    REAL(wp), INTENT(INOUT)  :: vw(:,:,:) !< product of u-velocity and v-velocity

    LOGICAL, INTENT(IN), OPTIONAL :: lacc


    TYPE(t_subset_range), POINTER            :: subset

    INTEGER  :: blk, cell, cellStart,cellEnd, level
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    subset => patch_3d%p_patch_2d(1)%cells%owned
    DO blk = subset%start_block, subset%end_block
      CALL get_index_range(subset, blk, cellStart, cellEnd)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO cell = cellStart, cellEnd

          hflR(cell,blk) = heatflux_total(cell,blk) * ( R(cell,1,blk) -1000.0_wp )
          fwR(cell,blk) = frshflux_volumetotal(cell,blk) * ( R(cell,1,blk) -1000.0_wp )
          tauxu(cell,blk) = topbc_windstress_u(cell,blk) * U(cell,1,blk)
          tauyv(cell,blk) = topbc_windstress_v(cell,blk) * V(cell,1,blk)


        !$ACC LOOP SEQ
        DO level=1,subset%vertical_levels(cell,blk)


          w_prismcenter(cell,level,blk) = 0.5*(w(cell,level,blk) + w(cell,level+1,blk))
          sigma0(cell,level,blk) = R(cell,level,blk) -1000.0_wp
          uT(cell,level,blk) = T(cell,level,blk) * u(cell,level,blk)
          uS(cell,level,blk) = S(cell,level,blk) * u(cell,level,blk)
          uR(cell,level,blk) = sigma0(cell,level,blk) * u(cell,level,blk)
          uu(cell,level,blk) = u(cell,level,blk) * u(cell,level,blk)

          vT(cell,level,blk) = T(cell,level,blk) * v(cell,level,blk)
          vS(cell,level,blk) = S(cell,level,blk) * v(cell,level,blk)
          vR(cell,level,blk) = sigma0(cell,level,blk) * v(cell,level,blk)
          vv(cell,level,blk) = v(cell,level,blk) * v(cell,level,blk)

          wT(cell,level,blk) = T(cell,level,blk) * w_prismcenter(cell,level,blk)
          wS(cell,level,blk) = S(cell,level,blk) * w_prismcenter(cell,level,blk)
          wR(cell,level,blk) = sigma0(cell,level,blk) * w_prismcenter(cell,level,blk)
          ww(cell,level,blk) = w_prismcenter(cell,level,blk) * w_prismcenter(cell,level,blk)

          uv(cell,level,blk) = u(cell,level,blk) * v(cell,level,blk)
          uw(cell,level,blk) = u(cell,level,blk) * w_prismcenter(cell,level,blk)
          vw(cell,level,blk) = v(cell,level,blk) * w_prismcenter(cell,level,blk)

          RR(cell,level,blk) = sigma0(cell,level,blk) * sigma0(cell,level,blk)
          SS(cell,level,blk) = S(cell,level,blk) * S(cell,level,blk)
          TT(cell,level,blk) = T(cell,level,blk) * T(cell,level,blk)


        END DO ! level
      END DO ! cell
      !$ACC END PARALLEL LOOP
    END DO !block
    !$ACC WAIT(1)

  END SUBROUTINE calc_eddydiag


  SUBROUTINE reset_ocean_monitor(monitor)
    TYPE(t_ocean_monitor) :: monitor
    monitor%volume(:)                     = 0.0_wp
!   monitor%kin_energy(:)                 = 0.0_wp
!   monitor%pot_energy(:)                 = 0.0_wp
!   monitor%total_energy(:)               = 0.0_wp
    monitor%total_salt(:)                 = 0.0_wp
!   monitor%vorticity(:)                  = 0.0_wp
!   monitor%enstrophy(:)                  = 0.0_wp
!   monitor%potential_enstrophy(:)        = 0.0_wp
    monitor%HeatFlux_ShortWave(:)         = 0.0_wp
    monitor%HeatFlux_LongWave(:)          = 0.0_wp
    monitor%HeatFlux_Sensible(:)          = 0.0_wp
    monitor%HeatFlux_Latent(:)            = 0.0_wp
    monitor%FrshFlux_SnowFall(:)          = 0.0_wp
    monitor%FrshFlux_TotalSalt(:)         = 0.0_wp
    monitor%FrshFlux_TotalOcean(:)        = 0.0_wp
    monitor%FrshFlux_TotalIce(:)          = 0.0_wp
    monitor%FrshFlux_VolumeIce(:)         = 0.0_wp
    monitor%FrshFlux_VolumeTotal(:)       = 0.0_wp
    monitor%HeatFlux_Relax(:)             = 0.0_wp
    monitor%FrshFlux_Relax(:)             = 0.0_wp
    monitor%TempFlux_Relax(:)             = 0.0_wp
    monitor%SaltFlux_Relax(:)             = 0.0_wp
    monitor%ice_framStrait(:)             = 0.0_wp
    monitor%florida_strait(:)             = 0.0_wp
    monitor%gibraltar(:)                  = 0.0_wp
    monitor%denmark_strait(:)             = 0.0_wp
    monitor%drake_passage(:)              = 0.0_wp
    monitor%indonesian_throughflow(:)     = 0.0_wp
    monitor%scotland_iceland(:)           = 0.0_wp
    monitor%mozambique(:)                 = 0.0_wp
    monitor%framStrait(:)                 = 0.0_wp
    monitor%beringStrait(:)               = 0.0_wp
    monitor%barentsOpening(:)             = 0.0_wp
    monitor%agulhas(:)                    = 0.0_wp
    monitor%agulhas_long(:)               = 0.0_wp
    monitor%agulhas_longer(:)             = 0.0_wp
  END SUBROUTINE reset_ocean_monitor

  SUBROUTINE calc_condep(patch_3d, condep, zgrad_rho, lacc)

    TYPE(t_patch_3d), TARGET, INTENT(in)  :: patch_3d
!    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    TYPE(t_patch), POINTER        :: patch_2d
    TYPE(t_subset_range), POINTER :: owned_cells

    REAL(wp), INTENT(inout) ::   condep(:,:)
    REAL(wp), INTENT(in) ::      zgrad_rho(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc


    INTEGER  :: blockNo, jc, start_index, end_index
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)


    patch_2d => patch_3D%p_patch_2d(1)
    owned_cells => patch_2d%cells%owned

    !ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index) SCHEDULE(dynamic)
    DO blockNo = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, blockNo, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc =  start_index, end_index
                condep(jc,blockNo) = &
             REAL(calc_max_condep(zgrad_rho(jc,:,blockNo), &
             patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)),KIND=wp)

      ENDDO
      !$ACC END PARALLEL LOOP
    ENDDO
    !$ACC WAIT(1)
    !ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE calc_condep

  SUBROUTINE calc_mld(patch_3d, mld, zgrad_rho, min_lev,sigcrit, lacc)

    TYPE(t_patch_3d), TARGET, INTENT(in)     :: patch_3d
    REAL(wp), INTENT(inout)                  :: mld(:,:)
    REAL(wp), INTENT(in)                     :: zgrad_rho(:,:,:)
    INTEGER,  INTENT(in)                     :: min_lev
    REAL(wp)                                 :: sigcrit
    LOGICAL, INTENT(IN), OPTIONAL :: lacc


    TYPE(t_patch), POINTER                   :: patch_2d
    TYPE(t_subset_range), POINTER            :: owned_cells

#ifdef __LVECTOR__
    REAL(wp) :: sigh(nproma)
    REAL(wp) :: masked_vertical_density_gradient
    REAL(wp) :: delta_h
    INTEGER :: jk
    INTEGER :: max_lev
#endif

    INTEGER  :: blockNo, jc, start_index, end_index
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)


    patch_2d => patch_3D%p_patch_2d(1)
    owned_cells => patch_2d%cells%owned

#ifndef __LVECTOR__
    ! Non-vector variant

    !ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index) SCHEDULE(dynamic)
    DO blockNo = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, blockNo, start_index, end_index)
      ! 2023-08 psam-DKRZ: use of GANG VECTOR here gives runtime error
      !$ACC PARALLEL LOOP VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc =  start_index, end_index

         mld(jc,blockNo) = calc_mixed_layer_depth(zgrad_rho(jc,:,blockNo),&
             sigcrit, &
             min_lev, &
             patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), &
             patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,:,blockNo), &
             patch_3d%p_patch_1d(1)%zlev_m(min_lev))

      ENDDO
      !$ACC END PARALLEL LOOP
    ENDDO
    !$ACC WAIT(1)
    !ICON_OMP_END_PARALLEL_DO

#else
    ! Vector variant

    !$ACC DATA CREATE(sigh) IF(lzacc)

    !ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, jc, jk, max_lev, &
    !ICON_OMP   & masked_vertical_density_gradient, delta_h, sigh) SCHEDULE(dynamic)
    DO blockNo = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, blockNo, start_index, end_index)

      ! This diagnostic calculates the mixed layer depth.
      ! It uses the incremental density increase between two
      ! levels and substracts it from the initial density criterion (sigcrit)
      ! and adds the level thickness (zzz) to zmld. This is done till
      ! the accumulated density increase between the surface and
      ! layer k is sigcrit or sigh = O, respectively.

      ! stabio(k) = insitu density gradient
      ! sigh = remaining density difference

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      mld(start_index:end_index,blockNo) = patch_3d%p_patch_1d(1)%zlev_m(min_lev)
      sigh(start_index:end_index) = sigcrit

      max_lev = MAXVAL(patch_3d%p_patch_1d(1)%dolic_c(start_index:end_index,blockNo))
      !$ACC END KERNELS
      !$ACC WAIT(1)

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jk = min_lev+1, max_lev
        DO jc = start_index, end_index
          IF (jk <= patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)) THEN
            IF (sigh(jc) > 1.e-6_wp) THEN
              masked_vertical_density_gradient = MAX(zgrad_rho(jc,jk,blockNo), 0.0_wp)
              delta_h = MIN( &
                  & sigh(jc) / (masked_vertical_density_gradient + 1.0E-19_wp), &
                  & patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,blockNo) &
                )
              sigh(jc) = sigh(jc) - delta_h * masked_vertical_density_gradient
              mld(jc,blockNo) = mld(jc,blockNo) + delta_h
            ELSE
              sigh(jc) = 0._wp
            END IF
          END IF
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)
    !ICON_OMP_END_PARALLEL_DO

    !$ACC END DATA
#endif
  END SUBROUTINE calc_mld

  !>
  !! Find level index of layer containing given depth.
  !!
  !! If depth is exactly at layer boundary, take layer above given depth.
  !! Return top layer for zero depth. Use bottom layer for too large depths.
  !!
  FUNCTION get_level_index_by_depth(patch_3d, depth) RESULT(level_index)


    TYPE(t_patch_3d ),TARGET, INTENT(in)     :: patch_3D
    REAL(dp), INTENT(in) :: depth
    INTEGER :: level_index

    level_index = 1
    DO WHILE(patch_3d%p_patch_1d(1)%zlev_i(level_index+1) < depth .AND. level_index < n_zlev)
      level_index = level_index + 1
    END DO

  END FUNCTION get_level_index_by_depth

END MODULE mo_ocean_diagnostics
