!
! mo_art_unit_conversion
! This module provides unit conversion routines
! Based on: Philipp Gasch (2016):
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

MODULE mo_art_unit_conversion
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  USE mo_var_metadata_types,            ONLY: t_var_metadata, POST_OP_RHO  
  USE mo_math_constants,                ONLY: rad2deg, deg2rad
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: art_massmix2density
  PUBLIC :: geocentric_latitude

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_massmix2density(prog_list, tracer_now, tracer_new, rho)
!<
! SUBROUTINE art_massmix2density
! Converts mass mixing ratio to concentrations
! Based on: -
! Part of Module: mo_art_unit_conversion
! Author: Daniel Rieger, KIT
! Initial Release: 2013-12-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  TYPE(t_var_list_ptr),TARGET,INTENT(in) :: &
    &  prog_list                          !< prognostic state list
  REAL(wp),INTENT(in)                :: &
    &  tracer_now(:,:,:,:),             & !< tracer concentrations at timelevel nnow
    &  rho(:,:,:)                         !< Density
  REAL(wp),INTENT(inout)             :: &
    &  tracer_new(:,:,:,:)                !< tracer concentrations at timelevel nnew
! Local variables
  TYPE(t_var_metadata), POINTER :: &  !< returns reference to tracer
    &  info                           !< metadata of current element
  INTEGER                       :: &
    &  iv                             !< loop index
  
  DO iv = 1, prog_list%p%nvars
    
    info=>prog_list%p%vl(iv)%p%info
    IF (info%tlev_source == 3 .AND. info%post_op%ipost_op_type == POST_OP_RHO) THEN
      tracer_new(:,:,:,info%ncontained) = tracer_now(:,:,:,info%ncontained) * rho(:,:,:)
    END IF
  END DO
  
END SUBROUTINE art_massmix2density
!!
!!-------------------------------------------------------------------------
!!
FUNCTION geocentric_latitude( geodetic_latitude, opt_squared_eccentricity )

  REAL(wp), INTENT(in)            ::  geodetic_latitude         !< in degrees
  REAL(wp), INTENT(in), OPTIONAL  ::  opt_squared_eccentricity

! Result
  REAL(wp)                        ::  geocentric_latitude       !< in degrees
! Local variables
  REAL(wp)                        ::  squared_eccentricity

  IF ( PRESENT(opt_squared_eccentricity) ) THEN
    squared_eccentricity = opt_squared_eccentricity
  ELSE
    squared_eccentricity = 0.00669438_wp
  END IF

  geocentric_latitude =  &
    &  rad2deg * ATAN((1.0_wp-squared_eccentricity)*TAN(deg2rad*geodetic_latitude))

END FUNCTION geocentric_latitude
!!
!!-------------------------------------------------------------------------
!!

END MODULE mo_art_unit_conversion
