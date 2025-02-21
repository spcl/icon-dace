! Implementation of tides by computation of the Sun's and Moon's full tidal potential
! This will be used in the pressure gradient calculation
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


!--------------------------------------------------
! timers definition
!needs:
!   USE mo_timer, ONLY: timer_start, timer_stop, timers_level, <timers_names>...
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

! DSL definitions 










!---------------------
! block definitions

!---------------------
! mappings






!---------------------
! connectivity









!---------------------
! generic types




!---------------------
! shortcuts







































!---------------------
! Upper-lower case
!----------------------------
MODULE mo_ocean_pressure_bc_conditions
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp, dp
   USE mtime,                    ONLY: datetime
  USE mo_exception,              ONLY: finish
  USE mo_model_domain,           ONLY: t_patch_3d, t_patch
  USE mo_ocean_nml,              ONLY: use_tides,   &
    & use_tides_SAL, atm_pressure_included_in_ocedyn,       &
    & OceanReferenceDensity_inv, vert_cor_type
  USE mo_physical_constants,     ONLY: grav
 ! USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
 ! USE mo_parallel_config,        ONLY: nproma
  USE mo_impl_constants,         ONLY: sea_boundary
  USE mo_util_dbg_prnt,          ONLY: dbg_print, debug_printValue
  USE mo_ocean_tides,            ONLY: calculate_tides_potential
  USE mo_ocean_surface_types,    ONLY: t_atmos_for_ocean
  USE mo_ocean_types,            ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,          ONLY: t_sea_ice
  USE mo_dynamics_config,        ONLY: nold
  USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
  USE mo_fortran_tools,          ONLY: set_acc_host_or_device


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: create_pressure_bc_conditions

  CHARACTER(LEN=12)  :: str_module = 'mo_ocean_pressure_bc_conditions'  ! Output of module for 1 line debug
  !-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE create_pressure_bc_conditions(patch_3d,ocean_state, p_as,sea_ice, current_time, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    TYPE(t_atmos_for_ocean), TARGET, INTENT(in)      :: p_as
    TYPE (t_sea_ice), TARGET, INTENT(in)             :: sea_ice
    TYPE(datetime), POINTER, INTENT(in)              :: current_time
    LOGICAL, INTENT(IN), OPTIONAL                    :: lacc

    REAL(wp) :: switch_atm_pressure, switch_vert_cor_type, switch_tides
    INTEGER :: jc, jb, start_index, end_index
    LOGICAL  :: lzacc
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    CHARACTER(len=*), PARAMETER :: routine = 'create_pressure_bc_conditions'

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    patch_2d  => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%all
    !------------------------------------------------------------------------
    IF (use_tides .OR. use_tides_SAL) THEN
      ! compute tidal potential
      CALL calculate_tides_potential(patch_3d,current_time,ocean_state%p_diag%rho, ocean_state%p_prog(nold(1))%h, &
           ocean_state%p_aux%bc_tides_potential, ocean_state%p_aux%bc_SAL_potential)
    ENDIF

    !------------------------------------------------------------------------
    ! total top potential (individiual terms)

    IF (atm_pressure_included_in_ocedyn) THEN
      switch_atm_pressure=1.0_wp
    ELSE
      switch_atm_pressure=0.0_wp
    ENDIF

    IF (vert_cor_type .EQ. 1 ) THEN
      switch_vert_cor_type=1.0_wp
    ELSE
      switch_vert_cor_type=0.0_wp
    ENDIF

    IF (use_tides .OR. use_tides_SAL) THEN
      switch_tides=1.0_wp
    ELSE
      switch_tides=0.0_wp
    ENDIF

    !------------------------------------------------------------------------
    ! total top potential

    IF ( (switch_atm_pressure + switch_vert_cor_type + switch_tides) > 0.0_wp ) THEN
    
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = start_index, end_index
          ocean_state%p_aux%bc_total_top_potential(jc,jb) =  &
            & ocean_state%p_aux%bc_tides_potential(jc,jb) &
            & + ocean_state%p_aux%bc_SAL_potential(jc,jb) &
            & + p_as%pao(jc,jb) * OceanReferenceDensity_inv * switch_atm_pressure & ! add acceleration by air pressure
            & + grav * sea_ice%draftave(jc,jb) * switch_vert_cor_type ! only zstar: pressure of sea ice on top of the first layer (divided by rhoref to create an acceleration)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)
    ENDIF

    IF (use_tides .OR. use_tides_SAL) THEN
      CALL dbg_print('tides_potential',  ocean_state%p_aux%bc_tides_potential, &
           str_module, 3, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
      CALL dbg_print('tides_SAL',      ocean_state%p_aux%bc_SAL_potential, &
           str_module, 3, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    ENDIF

  END SUBROUTINE create_pressure_bc_conditions
  !-------------------------------------------------------------------------


END MODULE mo_ocean_pressure_bc_conditions
!=============================================================================
