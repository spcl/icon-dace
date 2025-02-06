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

MODULE mo_comin_config

  USE mo_exception,            ONLY: message, message_text
  USE mo_impl_constants,       ONLY: max_dom, vname_len
#ifndef __NO_ICON_COMIN__
  USE comin_host_interface,    ONLY: t_comin_plugin_description
#endif

  IMPLICIT NONE
  PUBLIC

  TYPE :: t_comin_tracer_info
    CHARACTER(LEN=vname_len)           :: name
    INTEGER                            :: idx_tracer = -1 !< Index in ICON's tracer array
    INTEGER                            :: idx_turb   = -1 !< Index in ICON's ddt_tracer_turb array
    INTEGER                            :: idx_conv   = -1 !< Index in ICON's ddt_tracer_conv array
    TYPE(t_comin_tracer_info), POINTER :: next => NULL()
  END TYPE t_comin_tracer_info

  TYPE :: t_comin_icon_domain_config
    INTEGER :: nturb_tracer = 0
    INTEGER :: nconv_tracer = 0
    TYPE(t_comin_tracer_info), POINTER :: tracer_info_head
  END TYPE t_comin_icon_domain_config

  TYPE :: t_comin_config
    INTEGER :: nplugins = 0
#ifndef __NO_ICON_COMIN__
    TYPE(t_comin_plugin_description) :: plugin_list(16) !< list of dynamic libs (max: 16)
#endif
    TYPE(t_comin_icon_domain_config) :: comin_icon_domain_config(max_dom)
  END TYPE t_comin_config

  TYPE(t_comin_config), TARGET :: comin_config


CONTAINS

  SUBROUTINE configure_comin()
    CHARACTER(*), PARAMETER :: routine = "mo_comin_config::configure_comin"
#ifndef __NO_ICON_COMIN__
    INTEGER :: i

    WRITE (message_text,'(a)') "ICON Community Interface (ComIn)."
    CALL message(TRIM(routine),message_text)

    WRITE (message_text,'(i0,a)') comin_config%nplugins, " plugin(s) enabled."
    CALL message(TRIM(routine),message_text)

    DO i=1,comin_config%nplugins
       WRITE (message_text,'(i0,2a,2a)') i, ": ", TRIM(comin_config%plugin_list(i)%plugin_library), &
            TRIM(comin_config%plugin_list(i)%primary_constructor)
      CALL message(TRIM(routine),message_text)
    END DO
#endif
  END SUBROUTINE configure_comin

END MODULE mo_comin_config
