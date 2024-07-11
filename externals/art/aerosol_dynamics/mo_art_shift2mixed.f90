!
! mo_art_shift2mixed
! In case of a soluble mass content >5% shift content of insoluble mode to
! according mixed mode.
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

MODULE mo_art_shift2mixed
! ICON
  USE mo_kind,                          ONLY: wp
! ART

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_shift2mixed'

  PUBLIC :: art_shift2mixed

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_shift2mixed(solmass, totmass, tracer,                           &
  &                        istart, iend, kstart, kend, njsp3, itr0, itr3)
!<
! SUBROUTINE art_shift2mixed
! In case of a soluble mass content >5% shift content of insoluble mode to
! according mixed mode. 
! Based on: -
! Part of Module: mo_art_shift2mixed
! Author: Daniel Rieger, KIT
! Initial Release: 2017-05-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  INTEGER, INTENT(in)            :: &
    &  istart,iend,                 & !< Start and end index of nproma loop
    &  kstart, kend,                & !< Start and end index of vertical loop
    &  njsp3,                       & !< Number of species for mass mixing ratios
    &  itr0(:), itr3(:,:)             !< Index pairs to (dim1: shift from/to, dim2: njsp3)
  REAL(wp), INTENT(in)           :: &
    &  solmass(istart:iend,kstart:kend),                & !< Total soluble mass (mug kg-1)
    &  totmass(istart:iend,kstart:kend)                   !< Total insoluble mass (mug kg-1)
  REAL(wp), INTENT(inout)        :: &
    &  tracer(:,:,:)                  !< Tracer container (nproma,nlev,jsp)
!Local variables
  REAL(wp)                            :: &
    &  fr_sol                              !< Fraction of soluble mass
  INTEGER                             :: &
    &  jc, jk , jsp                        !< Loop indices

  DO jk= kstart, kend
!NEC$ ivdep
    DO jc= istart, iend
      IF (totmass(jc,jk) > 0._wp) THEN
        fr_sol  = solmass(jc,jk) / totmass(jc,jk)
        IF (fr_sol > 0.05_wp) THEN !perform shift
          tracer(jc,jk,itr0(2)) = tracer(jc,jk,itr0(2)) + tracer(jc,jk,itr0(1))
          tracer(jc,jk,itr0(1)) = 0._wp
          DO jsp = 1, njsp3
            tracer(jc,jk,itr3(2,jsp)) = tracer(jc,jk,itr3(2,jsp)) + tracer(jc,jk,itr3(1,jsp))
            tracer(jc,jk,itr3(1,jsp)) = 0._wp
          ENDDO
        ENDIF
      ENDIF
    ENDDO !jc
  ENDDO !jk

END SUBROUTINE art_shift2mixed
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_shift2mixed
