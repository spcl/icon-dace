!
! mo_art_feedback_icon
! This module provides the subroutines for chemical tracer processes
! that include feedback to ICON variables, especially ozone and water.
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

MODULE mo_art_feedback_icon
  ! ICON
  USE mo_kind,                        ONLY: wp
  USE mo_exception,                   ONLY: message
  USE mo_impl_constants,              ONLY: SUCCESS
  USE mo_var_list,                    ONLY: t_var_list_ptr
  USE mo_var,                         ONLY: t_var
  USE mo_art_config,                  ONLY: art_config
  USE mo_key_value_store,             ONLY: t_key_value_store
  USE mo_tracer_metadata_types,       ONLY: t_chem_meta
  USE mo_radiation_config,            ONLY: irad_o3
  ! ART
  USE mo_art_data,                    ONLY: p_art_data
  USE mo_art_atmo_data,               ONLY: t_art_atmo
  USE mo_art_wrapper_routines,        ONLY: art_get_indices_c

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_feedback_o3
  PUBLIC :: art_set_artconfig_o3_feedback
  PUBLIC :: art_dry_freezing_H2O

CONTAINS

SUBROUTINE art_feedback_o3(jg,o3_mmr_icon,o3_mmr_art)
!<
! SUBROUTINE art_feedback_o3
! This subroutine overwrites the ozone mass mixing ratio of ICON by the ART
! tracer
! Part of Module: mo_art_feedback_icon
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-15
!>
  INTEGER, INTENT(in)       ::  &
    &  jg                       !< patch on which computation is performed
  REAL(wp), INTENT(inout) ::  &
    &  o3_mmr_icon(:,:,:)       !< mass mixing ratio of ozone in ICON (kg / kg)
  REAL(wp), INTENT(in) :: &
    &  o3_mmr_art(:,:,:)       !< mass mixing ratio of ozone in ART (kg / kg)
  ! local variables
  INTEGER :: &
    &  jb,jk,jc                 !< loop indices etc.
  INTEGER :: &
    &  i_startidx, i_endidx
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  
  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
    DO jk = 1,art_atmo%nlev
!NEC$ ivdep
      DO jc = i_startidx,i_endidx
        o3_mmr_icon(jc,jk,jb) = o3_mmr_art(jc,jk,jb)
      END DO
    END DO
  END DO
END SUBROUTINE art_feedback_o3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_set_artconfig_o3_feedback(p_prog_list,jg,name_O3, dict_tracer)
!<
! SUBROUTINE art_feedback_o3
! This subroutine reads the XML element "feedback" from the key_value_store of O3 and
! sets it into art_config%O3_feedback
! Part of Module: mo_art_feedback_icon
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-15
!>
  TYPE(t_var_list_ptr), INTENT(in) :: &
    &  p_prog_list
  INTEGER, INTENT(in) :: &
    &  jg
  CHARACTER(LEN=*), INTENT(in) :: &
    &  name_O3
  TYPE(t_key_value_store), INTENT(in) :: &
    &  dict_tracer
  ! local variables
  INTEGER :: iTRO3, O3_feed
  INTEGER :: ierror, iv
  INTEGER, POINTER :: jsp

  CALL dict_tracer%get(name_O3,iTRO3,ierror)

  art_config(jg)%O3_feedback = 0

  IF (ierror == SUCCESS) THEN
     
    DO iv = 1, p_prog_list%p%nvars
  
      jsp=>p_prog_list%p%vl(iv)%p%info%ncontained
      IF (jsp == iTRO3) THEN
        SELECT TYPE(meta => p_prog_list%p%vl(iv)%p%info_dyn%tracer)
          CLASS IS (t_chem_meta)
            CALL meta%opt_meta%get('feedback',O3_feed,ierror)
            IF (ierror == SUCCESS) THEN
              IF (irad_o3 == 10) THEN
                art_config(jg)%O3_feedback = O3_feed
                IF (O3_feed == 1) THEN
                  CALL message('mo_art_feedback_icon:art_set_artconfig_o3_feedback',   &
                           &   'using ART Ozone for radiation')
                END IF
              END IF
            END IF

            EXIT
        END SELECT
      END IF

    END DO
  END IF
  
END SUBROUTINE art_set_artconfig_o3_feedback
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_dry_freezing_H2O(jg, p_tracer_H2O, vmr2Nconc)
!<
! SUBROUTINE art_dry_freezing_H2O
! This subroutine reads calculates the dry freezing of water vapour
! Part of Module: mo_art_feedback_icon
! Based on:
! - KASIMA
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-04
!>
  INTEGER, INTENT(in) :: &
    &  jg                  !< patch id on which computation is performed
  REAL(wp), INTENT(inout) :: &
    &  p_tracer_H2O(:,:,:) !< H2O number concentration (#/cm3)
  REAL(wp), INTENT(in) :: &
    &  vmr2Nconc(:,:,:)    !< conversion from volume mixing ratio (mol / mol) to number
                           !  concentration (#/cm3)
  ! local variables
  INTEGER :: &
    &  jc,jk,jb, i_startidx, i_endidx
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  ! ----------------------------------
  ! --- Dry freezing of water vapour
  ! ----------------------------------
  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
    DO jk = 1,art_atmo%nlev
!NEC$ ivdep
      DO jc = i_startidx,i_endidx
        p_tracer_H2O(jc,jk,jb) = p_tracer_H2O(jc,jk,jb) / vmr2Nconc(jc,jk,jb)

        p_tracer_H2O(jc,jk,jb) = min(                                          &
                 &                 p_tracer_H2O(jc,jk,jb),                     &
                 &                 exp( 24.306_wp                              &
                 &                    - 6144.9_wp / (art_atmo%temp(jc,jk,jb))  &
                 &                    ) / (art_atmo%pres(jc,jk,jb)/100._wp )   &
                 &                 )

        p_tracer_H2O(jc,jk,jb) = p_tracer_H2O(jc,jk,jb) * vmr2Nconc(jc,jk,jb)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE art_dry_freezing_H2O

END MODULE mo_art_feedback_icon
