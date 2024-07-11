! Initialisation of wave-atmosphere coupling
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

MODULE mo_wave_coupling_frame

  USE mo_exception,       ONLY: finish, message
  USE mo_model_domain,    ONLY: t_patch
  USE mo_run_config,      ONLY: ltimer
  USE mo_time_config,     ONLY: time_config
  USE mtime,              ONLY: timedeltaToString, MAX_TIMEDELTA_STR_LEN
  USE mo_coupling_config, ONLY: is_coupled_run, is_coupled_to_atmo
  USE mo_wave_atmo_coupling, ONLY: construct_wave_atmo_coupling
  USE mo_coupling_utils,  ONLY: cpl_def_main, cpl_enddef
  USE mo_timer,           ONLY: timer_start, timer_stop, timer_coupling_init

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_wave_coupling, destruct_wave_coupling
  PUBLIC :: nbr_inner_cells

  ! Output of module for debug
  CHARACTER(len=*), PARAMETER :: str_module = 'mo_wave_coupling_frame'

  INTEGER, SAVE         :: nbr_inner_cells

CONTAINS

  !>
  !! SUBROUTINE construct_wave_coupling -- the initialisation for the coupling
  !! of wave model and the atmosphere, through a coupler
  !!
  !! Note that the corresponding routine for the ATMO model construct_atmo_wave_coupling
  !! stored in src/atm_coupling/mo_atmo_waves_coupling_frame.f90
  !!
  SUBROUTINE construct_wave_coupling (p_patch)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)

    TYPE(t_patch), POINTER :: patch_horz

    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------

    INTEGER :: comp_id              ! component identifier
    INTEGER :: grid_id              ! grid identifier
    INTEGER :: cell_point_id

    INTEGER :: jg

    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN):: timestepstring

    CHARACTER(len=*), PARAMETER :: routine = str_module//':construct_wave_coupling'

    IF ( .NOT. is_coupled_run() ) RETURN

    IF (ltimer) CALL timer_start (timer_coupling_init)

    CALL message(str_module, 'Constructing the wave coupling frame.')

    jg = 1
    patch_horz => p_patch(jg)

    ! do basic initialisation of the component
    CALL cpl_def_main(routine,           & !in
                      patch_horz,        & !in
                      "icon_waves_grid", & !in
                      comp_id,           & !out
                      grid_id,           & !out
                      cell_point_id,     & !out
                      nbr_inner_cells)     !out

    ! get model timestep
    CALL timedeltaToString(time_config%tc_dt_model, timestepstring)

    IF ( is_coupled_to_atmo() ) THEN

      CALL message(str_module, 'Constructing the coupling frame wave-atmosphere.')

      CALL construct_wave_atmo_coupling( &
        comp_id, cell_point_id, timestepstring)

    END IF

    ! End definition of coupling fields and search
    CALL cpl_enddef(routine)

    IF (ltimer) CALL timer_stop(timer_coupling_init)

  END SUBROUTINE construct_wave_coupling

  !>
  !! SUBROUTINE destruct_wave_coupling -- the finalization for the coupling
  !! of wave model and the atmosphere, through a coupler
  !!
  SUBROUTINE destruct_wave_coupling ()

    CHARACTER(len=*), PARAMETER :: routine = str_module//':destruct_wave_coupling'

    IF ( is_coupled_run() ) THEN

      CALL message(str_module, 'Destructing the wave coupling frame.')

    END IF

  END SUBROUTINE destruct_wave_coupling

END MODULE mo_wave_coupling_frame
