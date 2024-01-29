!
! mo_art_prepare_aerosol
! Prepares mineral dust aerosol concentration for the KL06 scheme
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

MODULE mo_art_prepare_aerosol
  USE mo_kind,                          ONLY: wp
  USE mo_impl_constants,                ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,                     ONLY: finish
! ART ROUTINES
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_emiss_types,               ONLY: t_art_emiss2tracer
  USE mo_art_modes_linked_list,          ONLY: t_mode
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_prepare_dust_kl06

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_prepare_dust_kl06(jg, jb, istart, iend, kstart, kend, rho, tracer)
!<
! SUBROUTINE art_prepare_dust_kl06
! Prepares mineral dust aerosol concentration for the KL06 scheme
! Based on: -
! Part of Module: mo_art_prepare_aerosol
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-10
! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!>
  INTEGER,  INTENT(IN)            :: &
    &  jg, jb,                       & !< Domain and block index
    &  istart, iend,                 & !< Start/end index jc loop
    &  kstart, kend                    !< Start/end index jk loop
  REAL(wp), INTENT(in), TARGET    :: &
    &  rho(:,:)                        !< Density
  REAL(wp), INTENT(INOUT), TARGET :: &
    &  tracer(:,:,:)                   !< Tracer fields
! Local Variables
  INTEGER                         :: &
    &  jc, jk,                       & !< Loop indices
    &  imod, itr                       !< Loop indices
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_prepare_aerosol:art_prepare_dust_kl06"
  TYPE(t_mode), POINTER  :: &
    & current

  IF (p_art_data(jg)%tracer2aeroemiss%lisinit) THEN
    current=>p_art_data(jg)%tracer2aeroemiss%e2t_list%p%first_mode
    DO WHILE(ASSOCIATED(current))
      SELECT TYPE(this=>current%fields)
        TYPE IS(t_art_emiss2tracer)
          SELECT CASE(this%name)
            CASE('dust')
              IF(this%lcalcemiss) THEN
                p_art_data(jg)%diag%ndust_tot(:,:,jb) = 0._wp
                DO imod = 1, this%nmodes
                  DO jk = kstart, kend
                    DO jc = istart, iend
                      p_art_data(jg)%diag%ndust_tot(jc,jk,jb) =                                     &
                        &                                   p_art_data(jg)%diag%ndust_tot(jc,jk,jb) &
                        &                                 + tracer(jc,jk,this%itr0(imod))*rho(jc,jk)
                    ENDDO !jc
                  ENDDO !jk
                ENDDO ! imod
              ELSE
                CALL finish (thisroutine, 'no dust emissions found.')
              ENDIF ! this%lcalcemiss
            CASE DEFAULT
              !nothing to do
          END SELECT
      END SELECT
      current=>current%next_mode
    END DO
  END IF
  
END SUBROUTINE art_prepare_dust_kl06
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_prepare_aerosol
