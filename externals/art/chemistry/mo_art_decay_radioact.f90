!
! mo_art_decay_radioact
! This module provides the routine for the decay of radioactive particles.
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

MODULE mo_art_decay_radioact
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: ln2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_decay_radioact

CONTAINS

SUBROUTINE art_decay_radioact(dtime,istart,iend,nlev,halflife_tracer,tracer)
!<
! SUBROUTINE art_decay_radioact
! Calculates decay of radioactive particles
! Based on: -
! Part of Module: mo_art_decay_radioact
! Author: Daniel Rieger, KIT
! Initial Release: 2013-02-28
! Modifications:
! 2016-11-07: Daniel Rieger, KIT
! - Adaption to ART physics coupling concept, cleanup
!>
  REAL(wp), INTENT(IN)              :: &
    &  dtime,                          & !< time step
    &  halflife_tracer                   !< halflife time
  INTEGER, INTENT(IN)               :: &
    &  istart, iend,                   & !< start and end index of jc loop
    &  nlev                              !< number of vertical levels
  REAL(wp),INTENT(INOUT)            :: &
    &  tracer(:,:)                       !< concentration of the radioactive species
! Local variables
  REAL(wp)                          :: &
    &  decay_rate                      !< rate of radioactive decay
  INTEGER                           :: &
    &  jc, jk                            !< loop indizes

  IF(halflife_tracer > 0.0_wp) THEN

    decay_rate = ln2 / halflife_tracer

    DO jk=1,nlev
      DO jc=istart,iend
        tracer(jc,jk) = tracer(jc,jk) - tracer(jc,jk)*decay_rate*dtime
      ENDDO! jc
    ENDDO! jk

  ENDIF ! halflife_tracer
END SUBROUTINE art_decay_radioact

END MODULE mo_art_decay_radioact
