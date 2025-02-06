!
! mo_art_integration
! This module performs numerical time integrations. So far, only an explicit
! solver is included
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

MODULE mo_art_integration
! ICON Routines
  USE mo_kind,                          ONLY: wp
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_integrate_explicit
  
INTERFACE art_integrate_explicit
  MODULE PROCEDURE art_integrate_explicit_1d
  MODULE PROCEDURE art_integrate_explicit_2d
END INTERFACE art_integrate_explicit
  
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_integrate_explicit_1d(field, update_rate, dtime, istart,iend, &
              &                   opt_rho, opt_fac)
!<
! SUBROUTINE art_integrate_explicit_1d
! This subroutine performs a numerical integration in time
! using an explicit solution for a 1-dimensional field.
! Part of Module: mo_art_integration
! Based on: -
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-11
! Modifications:
! YYYY-MM-DD: <name>, <instituttion>
! - ...
!>
  REAL(wp),INTENT(inout) :: &
    &  field(:)            !< Field to be updated [UNIT kg-1 s-1], UNIT might be mug, kg
                           !    or just a number
  REAL(wp),INTENT(in) :: &
    &  update_rate(:),   & !< Update rate [UNIT kg-1 s-1] or [UNIT m-3 s-1], if opt_rho is passed!
                           !    (UNIT has to be the same as above)
    &  dtime               !< Integration time step [s]
  INTEGER,INTENT(in)  :: &
    &  istart,iend         !< Loop indizes (start and end)
  REAL(wp),OPTIONAL   :: &
    &  opt_rho(:)          !< Air density, needs to be passed if a volume specific update_rate is 
                           !    passed (see above) [kg m-3]
  REAL(wp),OPTIONAL   :: &
    &  opt_fac             !< Optional scalar factor, e.g. to transform mass emissions into number
                           !    emissions or vice versa
  ! Local variables
  INTEGER             :: &
    &  i                   !< Loop index
  REAL(wp)            :: &
    &  factor              !< Local instance of opt_fac if present
  
  IF(PRESENT(opt_fac)) THEN
    factor = opt_fac
  ELSE
    factor = 1.0_wp
  ENDIF
  
  IF(PRESENT(opt_rho)) THEN
    DO i = istart, iend
      field(i) = field(i) + update_rate(i) * dtime * factor * (1._wp/opt_rho(i))
    ENDDO
  ELSE
    DO i = istart, iend
      field(i) = field(i) + update_rate(i) * dtime * factor
    ENDDO
  ENDIF
  
END SUBROUTINE art_integrate_explicit_1d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_integrate_explicit_2d(field, update_rate, dtime, istart,iend, nlev, &
              &                      opt_rho, opt_fac, opt_kstart)
!<
! SUBROUTINE art_integrate_explicit_2d
! This subroutine performs a numerical integration in time
! using an explicit solution for a 2-dimensional field.
! Part of Module: mo_art_integration
! Based on: -
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-11
! Modifications:
! YYYY-MM-DD: <name>, <instituttion>
! - ...
!>
  REAL(wp),INTENT(inout) :: &
    &  field(:,:)          !< Field to be updated [UNIT kg-1 s-1], UNIT might be mug, kg 
                           !    or just a number
  REAL(wp),INTENT(in) :: &
    &  update_rate(:,:), & !< Update rate [UNIT kg-1 s-1] or [UNIT m-3 s-1], if opt_rho is passed!
                           !    (UNIT has to be the same as above)
    &  dtime               !< Integration time step [s]
  INTEGER,INTENT(in)  :: &
    &  istart,iend,      & !< Loop indizes (start and end)
    &  nlev
  REAL(wp),INTENT(in),OPTIONAL   :: &
    &  opt_rho(:,:)        !< Air density, needs to be passed if a volume specific update_rate is
                           !    passed (see above) [kg m-3]
  REAL(wp),INTENT(in),OPTIONAL   :: &
    &  opt_fac             !< Optional scalar factor, e.g. to transform mass emissions into number
                           !    emissions or vice versa
  INTEGER,INTENT(in),OPTIONAL :: &
    &  opt_kstart

  ! Local variables
  INTEGER             :: &
    &  i, k,             & !< Loop indices
    &  kstart
  REAL(wp)            :: &
    &  factor              !< Local instance of opt_fac if present
  
  IF(PRESENT(opt_fac)) THEN
    factor = opt_fac
  ELSE
    factor = 1.0_wp
  ENDIF
  
  IF(PRESENT(opt_kstart)) THEN
    kstart = opt_kstart
  ELSE
    kstart = 1
  ENDIF
  
  IF(PRESENT(opt_rho)) THEN
    DO k = kstart, nlev
      DO i = istart, iend
        field(i,k) = field(i,k) + update_rate(i,k) * dtime * factor * (1._wp/opt_rho(i,k))
      ENDDO
    ENDDO
  ELSE
    DO k = kstart,nlev
      DO i = istart, iend
        field(i,k) = field(i,k) + update_rate(i,k) * dtime * factor
      ENDDO
    ENDDO
  ENDIF
  
END SUBROUTINE art_integrate_explicit_2d
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_integration
