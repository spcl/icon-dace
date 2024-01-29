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

MODULE mo_coupling_config

  IMPLICIT NONE
  PRIVATE

  !>
  !! Namelist input to steer the coupling modes
  !! note that default is potentially overwritten in corresponding Namelist routine(s)
  !!
  LOGICAL :: config_coupled_to_ocean     = .FALSE.
  LOGICAL :: config_coupled_to_waves     = .FALSE.
  LOGICAL :: config_coupled_to_atmo      = .FALSE.
  LOGICAL :: config_coupled_to_hydrodisc = .FALSE.
  LOGICAL :: config_coupled_to_output    = .FALSE.


  ! variables
  PUBLIC :: config_coupled_to_ocean
  PUBLIC :: config_coupled_to_waves
  PUBLIC :: config_coupled_to_atmo
  PUBLIC :: config_coupled_to_hydrodisc
  PUBLIC :: config_coupled_to_output

  ! functions
  PUBLIC :: is_coupled_run
  PUBLIC :: is_coupled_to_ocean
  PUBLIC :: is_coupled_to_waves
  PUBLIC :: is_coupled_to_atmo
  PUBLIC :: is_coupled_to_hydrodisc
  PUBLIC :: is_coupled_to_output

CONTAINS

  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_run()

    is_coupled_run = config_coupled_to_ocean .OR.     &
      &              config_coupled_to_waves .OR.     &
      &              config_coupled_to_atmo  .OR.     &
      &              config_coupled_to_hydrodisc .OR. &
      &              config_coupled_to_output

  END FUNCTION is_coupled_run

  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_to_ocean()

    is_coupled_to_ocean = config_coupled_to_ocean

  END FUNCTION is_coupled_to_ocean

  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_to_waves()

    is_coupled_to_waves = config_coupled_to_waves

  END FUNCTION is_coupled_to_waves

  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_to_atmo()

    is_coupled_to_atmo = config_coupled_to_atmo

  END FUNCTION is_coupled_to_atmo

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_to_hydrodisc()

    is_coupled_to_hydrodisc = config_coupled_to_hydrodisc

  END FUNCTION is_coupled_to_hydrodisc

  LOGICAL FUNCTION is_coupled_to_output()

    is_coupled_to_output = config_coupled_to_output

  END FUNCTION is_coupled_to_output

  !------------------------------------------------------------------------

END MODULE mo_coupling_config
