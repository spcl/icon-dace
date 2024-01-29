!
! Classes and functions for the turbulent mixing package (tmx)
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

MODULE mo_tmx_time_integration_class

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish
  USE mo_surrogate_class,   ONLY: t_surrogate
  ! USE mo_variable_list, ONLY: t_variable_list

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: explicit_euler, implicit
  PUBLIC :: t_time_scheme

  ! Timestepping scheme
  ENUM, BIND(C)
    ENUMERATOR :: explicit_euler=1, implicit
  END ENUM

  TYPE, ABSTRACT :: t_time_scheme
    ! TYPE(t_variable_list) :: config_list
    INTEGER :: type
  CONTAINS
    PROCEDURE(step_iface), NOPASS, DEFERRED :: Step_forward
  END TYPE t_time_scheme

  ABSTRACT INTERFACE
    SUBROUTINE step_iface(process, dt)
      IMPORT :: t_surrogate, wp
      CLASS(t_surrogate), INTENT(inout) :: process
      REAL(wp), INTENT(in) :: dt
    END SUBROUTINE
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_time_integration_class'
  
! CONTAINS

END MODULE mo_tmx_time_integration_class
