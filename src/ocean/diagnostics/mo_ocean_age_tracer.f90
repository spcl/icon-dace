! Contains code for age tacer dynamics
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
MODULE mo_ocean_age_tracer
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: &
    & no_tracer,                                              &
    & n_zlev, bottom_drag_coeff,                              &
    & OceanReferenceDensity,                                  &
    & i_sea_ice,                                              &
    & ReferencePressureIndbars,                               &
    & age_tracer_inv_relax_time,                              &
    & l_relaxage_ice,                                         &
    & use_lbound_dirichlet!,        &

  USE mo_ocean_physics_types, ONLY: t_ho_params, v_params
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_dynamics_config,     ONLY: nold, nnew
  USE mo_impl_constants,      ONLY: success, max_char_length, min_dolic, sea
  USE mo_var_list,            ONLY: add_var
  USE mo_grib2,               ONLY: grib2_var, t_grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_FLT32 => CDI_DATATYPE_FLT32, &
    &                               DATATYPE_FLT64 => CDI_DATATYPE_FLT64, &
    &                               DATATYPE_INT8 => CDI_DATATYPE_INT8, &
    &                               DATATYPE_PACK16 => CDI_DATATYPE_PACK16, &
    &                               tstep_constant, GRID_LONLAT, GRID_UNSTRUCTURED
  USE mo_cdi_constants,       ONLY: grid_cell, grid_edge, grid_unstructured_cell, grid_unstructured_edge, &
    &                               grid_unstructured_vert, grid_vertex, GRID_ZONAL
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_var_groups,          ONLY: groups, max_groups
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_util_dbg_prnt,       ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state, t_onEdges_Pointer_3d_wp, t_onCells_HalfLevels_Pointer_wp, t_operator_coeff, t_hydro_ocean_diag
  USE mo_ocean_state,         ONLY: oce_config, ocean_default_list
  USE mo_physical_constants,  ONLY: grav
  USE mo_cf_convention
  USE mo_zaxis_type,          ONLY: &
    & za_depth_below_sea, za_depth_below_sea_half, za_surface
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sea_ice_types,       ONLY: t_sea_ice
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, global_max, sync_patch_array_mult
  USE mo_ocean_thermodyn,     ONLY: calculate_density_onColumn
  USE mo_ocean_math_operators,ONLY: div_oce_3d
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, &
    & timer_extra10, timer_extra11
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

!  PUBLIC :: init_age_tracer
  PUBLIC :: calc_age_tracer

CONTAINS

!  SUBROUTINE init_age_tracer(patch_3d, ocean_state_diag)
!    TYPE(t_patch_3D), TARGET, INTENT(in)      :: patch_3d
!    TYPE(t_hydro_ocean_diag), INTENT(inout)   :: ocean_state_diag
!
!    INTEGER :: jd
!    CHARACTER(LEN = *), PARAMETER :: routine = 'mo_ocean_age_tracer:init_age_tracer'
!    INTEGER :: status,i, alloc_cell_blocks, nblks_e, nblks_v
!    INTEGER :: datatype_flt
!
!    ! --- some information to print out as log
!    CALL message (TRIM(routine), "Age tracer setup:")
!    WRITE(message_text,'(a,i4)') "mode_age_tracer = ", mode_age_tracer
!    CALL message (TRIM(routine), message_text)
!    WRITE(message_text,'(a,i4)') "n_dlev = ", n_dlev
!    CALL message (TRIM(routine), message_text)
!    DO jd = 1, n_dlev+1
!      WRITE(message_text,'(a,i4,a,f10.3)') "rho_lev(", jd, ") = ", rho_lev(jd)
!      CALL message (TRIM(routine), message_text)
!    END DO
!    DO jd = 1, n_dlev
!      WRITE(message_text,'(a,i4,a,f10.3)') "rho_lev_cent(", jd, ") = ", rho_lev_cent(jd)
!      CALL message (TRIM(routine), message_text)
!    END DO
!
!  END SUBROUTINE init_age_tracer

  SUBROUTINE calc_age_tracer(patch_3d, ocean_state, jstep, p_ice, lacc)
    TYPE(t_patch_3d ), TARGET, INTENT(in)     :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(in) :: ocean_state
    INTEGER                                   :: jstep
    TYPE (t_sea_ice),              INTENT(IN) :: p_ice
    LOGICAL, INTENT(in), OPTIONAL             :: lacc
    REAL(wp), DIMENSION(:,:,:), POINTER       :: age_tracer
    REAL(wp), DIMENSION(:,:,:), POINTER       :: age_tracer_squared
    REAL(wp)                                  :: relax_strength

    INTEGER                                   :: jc, jb, je, jk, jd, blockNo, tracer_index, kk, elev
    INTEGER                                   :: i_startidx_c, i_endidx_c
    INTEGER                                   :: start_index, end_index
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: p_patch
    LOGICAL                       :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !write(*,*) "Before age tracer."

    ! First we increment the age tracers
    ! --- age_treacer
    age_tracer => ocean_state%p_prog(nnew(1))%tracer(:,:,:,3)
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    age_tracer(:,:,:) = age_tracer(:,:,:) + patch_3d%wet_c(:,:,:) * dtime    
    !$ACC END KERNELS
    ! --- age_tracer_squared
    IF (no_tracer>=4) THEN
      age_tracer_squared => ocean_state%p_prog(nnew(1))%tracer(:,:,:,4)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      age_tracer_squared(:,:,:) = age_tracer_squared(:,:,:) + age_tracer(:,:,:)*dtime*2._wp
      !$ACC END KERNELS
    ENDIF
    

    ! Now we relax the tracers at the surface to zero
    p_patch   => patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = i_startidx_c, i_endidx_c

        ! --- Calculate the relaxation strength
        IF (l_relaxage_ice .AND. i_sea_ice >=1) THEN
          relax_strength = (1.0_wp - p_ice%concsum(jc,jb)) * age_tracer_inv_relax_time  * dtime
        ELSE
          relax_strength = age_tracer_inv_relax_time  * dtime
        ENDIF
        
        ! --- relax the age
        age_tracer(jc,1,jb) = (1._wp - relax_strength) * age_tracer(jc,1,jb)
        
        ! --- relax the age squared
        ! relaxation strength is doubled for the squared tracer
        IF (no_tracer >= 4) THEN
          age_tracer_squared(jc,1,jb) = (1._wp - 2._wp * relax_strength) * age_tracer_squared(jc,1,jb)
        ENDIF
      
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)
    
    ! --- Check everything is still positive
    DO jb = all_cells%start_block, all_cells%end_block
      elev = p_patch%nlev
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jk = 1, elev
        DO jc = i_startidx_c, i_endidx_c
          ! --- age tracer
          IF (age_tracer(jc,jk,jb) < 0) THEN
            age_tracer(jc,jk,jb) = 0
          ENDIF
          
          ! --- age tracer squared
          IF (no_tracer >= 4) THEN
            IF (age_tracer_squared(jc,jk,jb) < 0)  THEN
              age_tracer_squared(jc,jk,jb) = 0
            ENDIF
          ENDIF
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)
    
    !write(*,*) "After age tracer."

  END SUBROUTINE calc_age_tracer

END MODULE mo_ocean_age_tracer
