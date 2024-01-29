!
! mo_art_drydepo_radioact
! This module provides the routine for the dry deposition of radioactive particles.
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

MODULE mo_art_drydepo_radioact
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_art_config,                    ONLY: t_art_config
  USE mo_fortran_tools,                 ONLY: t_ptr_tracer
! ART
  USE mo_art_data,                      ONLY: t_art_data

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_drydepo_radioact'

  PUBLIC :: art_drydepo_radioact

CONTAINS


SUBROUTINE art_drydepo_radioact(turb_tracer, art_config, itracer, vdep_const, &
  &                             istart, iend, jb, p_art_data)
!<
! SUBROUTINE art_drydepo_radioact
! This subroutine sets constant deposition velocities for
! radioactive tracers
! Based on: -
! Part of Module: mo_art_drydepo_radioact
! Author: Daniel Rieger, KIT
! Initial Release: 2013-02-28
! Modifications:
! 2016-11-04: Daniel Rieger, KIT
! - Adapted to unified ART-physics structure
!>
  TYPE(t_ptr_tracer),INTENT(in) :: &
    &  turb_tracer(:)                !< List of turbulent tracers
  TYPE(t_art_config)            :: &
    &  art_config                    !< ART configuration state
  INTEGER,INTENT(in)            :: &
    &  itracer,                    & !< Index of tracer in p_tracer_new(:,:,:,itracer)
    &  istart, iend,               & !< Start and end indices of nproma loop
    &  jb                            !< block index
  REAL(wp),INTENT(in)           :: &
    &  vdep_const                    !< constant deposition velocity for this tracer
  TYPE(t_art_data),INTENT(inout):: &
    &  p_art_data
! Local variables
  INTEGER                       :: &
    &  iturb, jc,                  & !< Loop indexes
    &  jsp                           !< index of current turbulent tracer field in turb_tracer structure

  DO iturb = 1, art_config%nturb_tracer
    IF (turb_tracer(iturb)%idx_tracer == itracer) THEN
      ! Get the index of the turbulent tracer field
      jsp = iturb
    ENDIF
  ENDDO
  DO jc = istart, iend
    p_art_data%turb_fields%vdep(jc,jb,jsp) = vdep_const
  ENDDO

END SUBROUTINE art_drydepo_radioact
END MODULE mo_art_drydepo_radioact
