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

MODULE mo_ocean_coupling_frame

  USE mo_exception,           ONLY: message, finish
  USE mo_run_config,          ONLY: ltimer
  USE mo_timer,               ONLY: timer_start, timer_stop, &
       &                            timer_coupling_init
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mtime,                  ONLY: timedeltaToString, MAX_TIMEDELTA_STR_LEN

  !-------------------------------------------------------------
  ! For the coupling
  !
  USE mo_coupling_utils,      ONLY: cpl_def_main, cpl_enddef
  USE mo_coupling_config,     ONLY: is_coupled_run, is_coupled_to_atmo, &
    &                               is_coupled_to_output
  USE mo_output_coupling,     ONLY: construct_output_coupling, &
    &                               construct_output_coupling_finalize
  USE mo_ocean_atmo_coupling, ONLY: construct_ocean_atmo_coupling
  USE mo_time_config,         ONLY: time_config

  !-------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_ocean_coupling_frame' ! Output of module for debug

  PUBLIC :: construct_ocean_coupling, destruct_ocean_coupling
  PUBLIC :: nbr_inner_cells

  INTEGER, SAVE :: nbr_inner_cells

CONTAINS

  !--------------------------------------------------------------------------
  ! Prepare the coupling
  !
  ! For the time being this could all go into a subroutine which is
  ! common to atmo and ocean. Does this make sense if the setup deviates
  ! too much in future.
  !------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE construct_ocean_coupling(patch_3d)
    TYPE(t_patch_3d ), TARGET, INTENT(in)    :: patch_3d

    INTEGER                :: patch_no
    TYPE(t_patch), POINTER :: patch_horz

    INTEGER :: comp_id, output_comp_id
    INTEGER :: grid_id
    INTEGER :: cell_point_id, vertex_point_id

    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: timestepstring

    CHARACTER(LEN=*), PARAMETER :: routine = str_module // ':construct_ocean_coupling'

    IF (.NOT. is_coupled_run()) RETURN

    IF (ltimer) CALL timer_start(timer_coupling_init)

    CALL message(str_module, 'Constructing the ocean coupling frame.')

    patch_no = 1
    patch_horz => patch_3d%p_patch_2d(patch_no)

    ! Do basic initialisation of the component
    IF( is_coupled_to_output() ) THEN
      CALL cpl_def_main(routine,           & !in
                        patch_horz,        & !in
                        "icon_ocean_grid", & !in
                        comp_id,           & !out
                        output_comp_id,    & !out
                        grid_id,           & !out
                        cell_point_id,     & !out
                        vertex_point_id,   & !out
                        nbr_inner_cells)     !out
    ELSE
      CALL cpl_def_main(routine,           & !in
                        patch_horz,        & !in
                        "icon_ocean_grid", & !in
                        comp_id,           & !out
                        grid_id,           & !out
                        cell_point_id,     & !out
                        nbr_inner_cells)     !out
    ENDIF

    ! get model timestep
    CALL timedeltaToString(time_config%tc_dt_model, timestepstring)

    IF( is_coupled_to_output() ) THEN

      CALL message(str_module, 'Constructing the coupling frame ocean-output.')

      CALL construct_output_coupling ( &
        patch_3d%p_patch_2d(1:), output_comp_id, cell_point_id, &
        vertex_point_id, timestepstring)

    END IF

    IF ( is_coupled_to_atmo() ) THEN

      ! Construct coupling frame for ocean-atmosphere
      CALL message(str_module, 'Constructing the coupling frame ocean-atmosphere.')

      CALL construct_ocean_atmo_coupling( &
        patch_3d, comp_id, grid_id, cell_point_id, timestepstring, &
        nbr_inner_cells)

    END IF

    CALL cpl_enddef(routine)

    ! finalizes the output coupling
    IF( is_coupled_to_output() ) CALL construct_output_coupling_finalize()

    IF (ltimer) CALL timer_stop(timer_coupling_init)

  END SUBROUTINE construct_ocean_coupling

  !--------------------------------------------------------------------------

!<Optimize:inUse>
  SUBROUTINE destruct_ocean_coupling()

    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':destruct_ocean_coupling'

    IF (is_coupled_run()) THEN

      CALL message(str_module, 'Destructing the coean coupling frame.')

    END IF

  END SUBROUTINE destruct_ocean_coupling

  !--------------------------------------------------------------------------

END MODULE mo_ocean_coupling_frame

