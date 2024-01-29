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

MODULE mo_ser_debug

  USE mo_kind,               ONLY: vp, wp
  USE mo_exception,          ONLY: warning
  !$ser verbatim USE mo_ser_common,         ONLY: init
  USE mtime,                 ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_ser_nml,            ONLY: ser_debug

  IMPLICIT NONE

  PUBLIC :: serialize_debug_output, ser_debug_on
  
  LOGICAL :: ser_debug_on = .FALSE.

  CONTAINS

  SUBROUTINE serialize_debug_output(nproma, nlev, jb, id, lcpu_only, r1d1, r1d2, r1d3, r1d4, r1d5, &
                                                                 r2d1, r2d2, r2d3, r2d4, r2d5, r2d6, r2d7, &
                                                                 r3d1, r3d2, r3d3, r3d4, r3d5, r3d6)

    INTEGER, INTENT(IN) :: nproma, nlev, jb, id
    LOGICAL :: lcpu_only
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r1d1(:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r1d2(:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r1d3(:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r1d4(:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r1d5(:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r2d1(:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r2d2(:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r2d3(:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r2d4(:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r2d5(:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r2d6(:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r2d7(:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r3d1(:,:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r3d2(:,:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r3d3(:,:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r3d4(:,:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r3d5(:,:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: r3d6(:,:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date

    !$ser verbatim IF(ser_debug .and. ser_debug_on) THEN
    !$ser verbatim   CALL warning('SER:mo_ser_debug:debug-output','Serialization is active!')
#if defined(_OPENACC)
    !$ser verbatim IF(.NOT. lcpu_only) THEN
    !$ser verbatim   CALL warning('GPU:mo_ser_debug:debug-output','GPU HOST synchronization forced by serialization!')
    !$ser verbatim   !$ACC UPDATE HOST(r1d1) ASYNC(1) IF(PRESENT(r1d1))
    !$ser verbatim   !$ACC UPDATE HOST(r1d2) ASYNC(1) IF(PRESENT(r1d2))
    !$ser verbatim   !$ACC UPDATE HOST(r1d3) ASYNC(1) IF(PRESENT(r1d3))
    !$ser verbatim   !$ACC UPDATE HOST(r1d4) ASYNC(1) IF(PRESENT(r1d4))
    !$ser verbatim   !$ACC UPDATE HOST(r1d5) ASYNC(1) IF(PRESENT(r1d5))
    !$ser verbatim   !$ACC UPDATE HOST(r2d1) ASYNC(1) IF(PRESENT(r2d1))
    !$ser verbatim   !$ACC UPDATE HOST(r2d2) ASYNC(1) IF(PRESENT(r2d2))
    !$ser verbatim   !$ACC UPDATE HOST(r2d3) ASYNC(1) IF(PRESENT(r2d3))
    !$ser verbatim   !$ACC UPDATE HOST(r2d4) ASYNC(1) IF(PRESENT(r2d4))
    !$ser verbatim   !$ACC UPDATE HOST(r2d5) ASYNC(1) IF(PRESENT(r2d5))
    !$ser verbatim   !$ACC UPDATE HOST(r2d6) ASYNC(1) IF(PRESENT(r2d6))
    !$ser verbatim   !$ACC UPDATE HOST(r2d7) ASYNC(1) IF(PRESENT(r2d7))
    !$ser verbatim   !$ACC UPDATE HOST(r3d1) ASYNC(1) IF(PRESENT(r3d1))
    !$ser verbatim   !$ACC UPDATE HOST(r3d2) ASYNC(1) IF(PRESENT(r3d2))
    !$ser verbatim   !$ACC UPDATE HOST(r3d3) ASYNC(1) IF(PRESENT(r3d3))
    !$ser verbatim   !$ACC UPDATE HOST(r3d4) ASYNC(1) IF(PRESENT(r3d4))
    !$ser verbatim   !$ACC UPDATE HOST(r3d5) ASYNC(1) IF(PRESENT(r3d5))
    !$ser verbatim   !$ACC UPDATE HOST(r3d6) ASYNC(1) IF(PRESENT(r3d6))
    !$ser verbatim   !$ACC WAIT(1)
    !$ser verbatim ENDIF
#endif

    !$ser verbatim   CALL datetimeToString(time_config%tc_current_date, date)
    !$ser verbatim   CALL init('icon')
    !$ser savepoint nwp_debug-output nproma=nproma nlev=nlev jb=jb date=TRIM(date) id=id
    !$ser mode write
    !$ser data nwp_debug_r1d1=r1d1(:) IF (PRESENT(r1d1))
    !$ser data nwp_debug_r1d2=r1d2(:) IF (PRESENT(r1d2))
    !$ser data nwp_debug_r1d3=r1d3(:) IF (PRESENT(r1d3))
    !$ser data nwp_debug_r1d4=r1d4(:) IF (PRESENT(r1d4))
    !$ser data nwp_debug_r1d5=r1d5(:) IF (PRESENT(r1d5))
    !$ser data nwp_debug_r2d1=r2d1(:,:) IF (PRESENT(r2d1))
    !$ser data nwp_debug_r2d2=r2d2(:,:) IF (PRESENT(r2d2))
    !$ser data nwp_debug_r2d3=r2d3(:,:) IF (PRESENT(r2d3))
    !$ser data nwp_debug_r2d4=r2d4(:,:) IF (PRESENT(r2d4))
    !$ser data nwp_debug_r2d5=r2d5(:,:) IF (PRESENT(r2d5))
    !$ser data nwp_debug_r2d6=r2d6(:,:) IF (PRESENT(r2d6))
    !$ser data nwp_debug_r2d7=r2d7(:,:) IF (PRESENT(r2d7))
    !$ser data nwp_debug_r3d1=r3d1(:,:,:) IF (PRESENT(r3d1))
    !$ser data nwp_debug_r3d2=r3d2(:,:,:) IF (PRESENT(r3d2))
    !$ser data nwp_debug_r3d3=r3d3(:,:,:) IF (PRESENT(r3d3))
    !$ser data nwp_debug_r3d4=r3d4(:,:,:) IF (PRESENT(r3d4))
    !$ser data nwp_debug_r3d5=r3d5(:,:,:) IF (PRESENT(r3d5))
    !$ser data nwp_debug_r3d6=r3d6(:,:,:) IF (PRESENT(r3d6))

    !$ser verbatim END IF

  END SUBROUTINE serialize_debug_output

END MODULE mo_ser_debug
