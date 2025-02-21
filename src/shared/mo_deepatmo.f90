! This module contains routines related to the deep atmosphere.
!
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

!----------------------------
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

!----------------------------

MODULE mo_deepatmo

  USE mo_kind,            ONLY: wp
  USE mo_impl_constants,  ONLY: SUCCESS
  USE mo_fortran_tools,   ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: deepatmo_htrafo

  INTERFACE deepatmo_htrafo
    MODULE PROCEDURE deepatmo_htrafo_a
    MODULE PROCEDURE deepatmo_htrafo_b
  END INTERFACE deepatmo_htrafo

CONTAINS !..................................................

  !>
  !! Transformation of the height coordinate
  !!
  !! Variant a: in- is also out-field
  !!
  SUBROUTINE deepatmo_htrafo_a( z_inout,             & !in/out
    &                           nblks_nproma_npromz, & !in
    &                           start_end_levels,    & !in
    &                           radius,              & !in
    &                           trafo_type,          & !in
    &                           ierror,              & !optout
    &                           lacc)                  !optin

    ! In/out variables
    REAL(wp),         INTENT(INOUT)         :: z_inout(:,:,:)         ! Height coordinate field
    INTEGER,          INTENT(IN)            :: nblks_nproma_npromz(3) ! Number of blocks, block length
                                                                      ! and length of last block
    INTEGER,          INTENT(IN)            :: start_end_levels(2)    ! Start (upermost) and end (lowermost) levels
    REAL(wp),         INTENT(IN)            :: radius                 ! Radius of Earth
    CHARACTER(LEN=7), INTENT(IN)            :: trafo_type             ! Type of transformation
    INTEGER,          INTENT(OUT), OPTIONAL :: ierror                 ! Optional error flag
    LOGICAL,          INTENT(IN),  OPTIONAL :: lacc                   ! Optional flag for use of OpenACC

    ! Local variables
    REAL(wp) :: trafo_factor
    INTEGER  :: z_inout_shape(3)
    INTEGER  :: jb, jk, jc  ! (jc is habitual placeholder for jc, je, jv)
    INTEGER  :: nlen
    LOGICAL  :: is_ierror_present
    LOGICAL  :: lzacc

    !-----------------------------------------------------------------------

    ! Unfortunately, there are subroutines, 
    ! which are conditionally calling this subroutine,
    ! that allow for execution on GPUs. Therefore, we have no other choice 
    ! but to implement OpenACC directives here, too.
    CALL set_acc_host_or_device(lzacc, lacc)

    is_ierror_present = PRESENT(ierror)

    ! Initialize error flag with error:
     IF (is_ierror_present) ierror = 1

    ! Some consistency checks
    z_inout_shape = SHAPE(z_inout)
    IF (.NOT. (radius > 0._wp)) THEN
      ! We need a positive-definite radius
      RETURN
    ELSEIF(ANY(start_end_levels(:) < 1)) THEN
      ! Invalid level indices
      RETURN
    ELSEIF(start_end_levels(1) > start_end_levels(2)) THEN
      ! Start level is greater than end level: nothing needs to be done
      RETURN
    ELSEIF (ANY(nblks_nproma_npromz(:) < [1,1,0])) THEN
      ! Invalid boundaries
      RETURN
    ELSEIF (nblks_nproma_npromz(3) > nblks_nproma_npromz(2)) THEN
      ! We expect npromz <(=) nproma
      RETURN
    ELSEIF (z_inout_shape(1) /= nblks_nproma_npromz(2)) THEN
      ! We expect z_inout_shape(1) = nproma
      RETURN
    ELSEIF (z_inout_shape(3) /= nblks_nproma_npromz(1)) THEN
      ! We expect z_inout_shape(3) = nblks
      RETURN
    ELSEIF (start_end_levels(2) > z_inout_shape(2)) THEN
      ! End level index is out of bounds
      RETURN
    ENDIF

    SELECT CASE(trafo_type) 
    CASE('zgpot2z')  
      ! Transform geopotential height z_gpot into 
      ! geometric height z by means of 
      ! z = z_gpot / ( 1 - z_gpot / a ), 
      ! where a is radius of Earth
      trafo_factor = -1._wp / radius
    CASE('z2zgpot')
      ! Transform geometric height z into 
      ! geopotential height z_gpot by means of 
      ! z_gpot = z / ( 1 + z / a ), 
      ! where a is radius of Earth
      trafo_factor = 1._wp / radius
    CASE DEFAULT
      ! Return in case of an invalid trafo type specifier
      RETURN
    END SELECT

    !$ACC DATA &
    !$ACC   PRESENT(z_inout, nblks_nproma_npromz, start_end_levels) &
    !$ACC   IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_nproma_npromz(1)

      IF (jb /= nblks_nproma_npromz(1)) THEN
        nlen = nblks_nproma_npromz(2)
      ELSE
        nlen = nblks_nproma_npromz(3)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = start_end_levels(1), start_end_levels(2)
          DO jc = nlen+1, nblks_nproma_npromz(2)
            z_inout(jc,jk,jb) = 0.0_wp
          ENDDO  !jc
        ENDDO  !jk
        !$ACC END PARALLEL
      ENDIF
      
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = start_end_levels(1), start_end_levels(2)
        DO jc = 1, nlen
          z_inout(jc,jk,jb) = z_inout(jc,jk,jb) &
            &               / ( 1._wp + trafo_factor * z_inout(jc,jk,jb) )
        ENDDO  !jc
      ENDDO  !jk
    !$ACC END PARALLEL
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

    IF (is_ierror_present) ierror = SUCCESS

  END SUBROUTINE deepatmo_htrafo_a

  !----------------------------------------------------------------------------

  !>
  !! Transformation of the height coordinate
  !!
  !! Variant b: in- and out-fields differ
  !!
  SUBROUTINE deepatmo_htrafo_b( z_in,                & !in
    &                           z_out,               & !out
    &                           nblks_nproma_npromz, & !in
    &                           start_end_levels,    & !in
    &                           radius,              & !in
    &                           trafo_type,          & !in
    &                           ierror,              & !optout
    &                           lacc)                  !optin

    ! In/out variables
    REAL(wp),         INTENT(IN)            :: z_in(:,:,:)            ! Input height coordinate field
    REAL(wp),         INTENT(OUT)           :: z_out(:,:,:)           ! Output height coord. field
    INTEGER,          INTENT(IN)            :: nblks_nproma_npromz(3) ! Number of blocks, block length
                                                                      ! and length of last block
    INTEGER,          INTENT(IN)            :: start_end_levels(2)    ! Start (upermost) and end (lowermost) levels
    REAL(wp),         INTENT(IN)            :: radius                 ! Radius of Earth
    CHARACTER(LEN=7), INTENT(IN)            :: trafo_type             ! Type of transformation
    INTEGER,          INTENT(OUT), OPTIONAL :: ierror                 ! Optional error flag
    LOGICAL,          INTENT(IN),  OPTIONAL :: lacc                   ! Optional flag for use of OpenACC

    ! Local variables
    REAL(wp) :: trafo_factor
    INTEGER  :: z_in_shape(3), z_out_shape(3)
    INTEGER  :: jb, jk, jc  ! (jc is habitual placeholder for jc, je, jv)
    INTEGER  :: nlen
    LOGICAL  :: is_ierror_present
    LOGICAL  :: lzacc

    !-----------------------------------------------------------------------

    ! Unfortunately, there are subroutines, 
    ! which are conditionally calling this subroutine,
    ! that allow for execution on GPUs. Therefore, we have no other choice 
    ! but to implement OpenACC directives here, too.
    CALL set_acc_host_or_device(lzacc, lacc)

    is_ierror_present = PRESENT(ierror)

    ! Initialize error flag with error:
    IF (is_ierror_present) ierror = 1

    ! Some consistency checks
    z_in_shape  = SHAPE(z_in)
    z_out_shape = SHAPE(z_out)
    IF (.NOT. (radius > 0._wp)) THEN
      ! We need a positive-definite radius
      RETURN
    ELSEIF(ANY(start_end_levels(:) < 1)) THEN
      ! Invalid level indices
      RETURN
    ELSEIF(start_end_levels(1) > start_end_levels(2)) THEN
      ! Start level is greater than end level: nothing needs to be done
      RETURN
    ELSEIF (ANY(nblks_nproma_npromz(:) < [1,1,0])) THEN
      ! Invalid boundaries
      RETURN
    ELSEIF (nblks_nproma_npromz(3) > nblks_nproma_npromz(2)) THEN
      ! We expect npromz <(=) nproma
      RETURN
    ELSEIF (ANY(z_in_shape(:) /= z_out_shape(:))) THEN
      ! We expect z_in and z_out to be of the same shape
      RETURN
    ELSEIF (z_in_shape(1) /= nblks_nproma_npromz(2)) THEN
      ! We expect z_in(out)_shape(1) = nproma
      RETURN
    ELSEIF (z_in_shape(3) /= nblks_nproma_npromz(1)) THEN
      ! We expect z_in(out)_shape(3) = nblks
      RETURN
    ELSEIF (start_end_levels(2) > z_in_shape(2)) THEN
      ! End level index is out of bounds
      RETURN
    ENDIF

    SELECT CASE(trafo_type) 
    CASE('zgpot2z')  
      ! Transform geopotential height z_gpot into 
      ! geometric height z by means of 
      ! z = z_gpot / ( 1 - z_gpot / a ), 
      ! where a is radius of Earth
      trafo_factor = -1._wp / radius
    CASE('z2zgpot')
      ! Transform geometric height z into 
      ! geopotential height z_gpot by means of 
      ! z_gpot = z / ( 1 + z / a ), 
      ! where a is radius of Earth
      trafo_factor = 1._wp / radius
    CASE DEFAULT
      ! Return in case of an invalid trafo type specifier
      RETURN
    END SELECT

    !$ACC DATA &
    !$ACC   PRESENT(z_in, z_out, nblks_nproma_npromz, start_end_levels) &
    !$ACC   IF(lzacc)
      
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_nproma_npromz(1)

      IF (jb /= nblks_nproma_npromz(1)) THEN
        nlen = nblks_nproma_npromz(2)
      ELSE
        nlen = nblks_nproma_npromz(3)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = start_end_levels(1), start_end_levels(2)
          DO jc = nlen+1, nblks_nproma_npromz(2)
            ! z_out has INTENT(OUT)
            z_out(jc,jk,jb) = 0.0_wp
          ENDDO  !jc
        ENDDO  !jk
        !$ACC END PARALLEL
      ENDIF

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = start_end_levels(1), start_end_levels(2)
        DO jc = 1, nlen
          z_out(jc,jk,jb) = z_in(jc,jk,jb) &
            &             / ( 1._wp + trafo_factor * z_in(jc,jk,jb) )
        ENDDO  !jc
      ENDDO  !jk
    !$ACC END PARALLEL
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

    IF (is_ierror_present) ierror = SUCCESS

  END SUBROUTINE deepatmo_htrafo_b

END MODULE mo_deepatmo
