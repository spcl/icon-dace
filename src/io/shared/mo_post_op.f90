! This module defines small arithmetic operations ("post-ops") as
! post-processing tasks.
!
! These post-processing tasks are restricted to point-wise
! operations (no halo synchronization) of a single field, like
! value scaling.
!
! @note The "post-ops" are performed at output time and DO NOT
!       MODIFY THE FIELD ITSELF.
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

MODULE mo_post_op

  USE mo_kind,                ONLY: dp, sp
  USE mo_var_metadata_types,  ONLY: t_post_op_meta, POST_OP_NONE,    &
    &                               POST_OP_SCALE, POST_OP_LUC, POST_OP_LIN2DBZ
#ifndef __NO_ICON_ATMO__
  USE mo_lnd_nwp_config,      ONLY: convert_luc_ICON2GRIB
#endif
  USE mo_exception,           ONLY: finish
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device
#ifdef _OPENACC
  USE openacc,                  ONLY: acc_is_present
#endif
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: perform_post_op

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_post_op')

  INTERFACE perform_post_op
    MODULE PROCEDURE perform_post_op_r2D
    MODULE PROCEDURE perform_post_op_s2D
    MODULE PROCEDURE perform_post_op_i2D
    MODULE PROCEDURE perform_post_op_r3D
    MODULE PROCEDURE perform_post_op_s3D
    MODULE PROCEDURE perform_post_op_i3D
  END INTERFACE perform_post_op

CONTAINS

  !> Performs small arithmetic operations ("post-ops") on 2D REAL field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_r2D(post_op, field2D, opt_inverse, lacc)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    REAL(dp),              INTENT(INOUT) :: field2D(:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    LOGICAL , OPTIONAL   , INTENT(IN)    :: lacc
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_r2D"
    REAL(dp) :: scalfac, lowlim   ! scale factor, lower limit
    INTEGER  :: idim(2), l1,l2
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    idim = SHAPE(field2D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          scalfac = 1._dp/post_op%arg1%rval
        ELSE
          scalfac = post_op%arg1%rval
        ENDIF
      ELSE
        scalfac = post_op%arg1%rval
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2), SCHEDULE(runtime)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) COPYIN(idim) IF(acc_is_present(field2D) .AND. lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO l2=1,idim(2)
        DO l1=1,idim(1)
          field2D(l1,l2) = field2D(l1,l2) * scalfac
        END DO
      END DO
      !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL
      !
    CASE (POST_OP_LIN2DBZ)

      lowlim  = post_op%arg1%rval
      scalfac = 10.0_dp/LOG(10._dp)

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2), SCHEDULE(runtime)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) COPYIN(idim) IF(acc_is_present(field2D) .AND. lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO l2=1,idim(2)
        DO l1=1,idim(1)
          field2D(l1,l2) = scalfac * LOG( MAX( field2D(l1,l2), lowlim) )
        END DO
      END DO
      !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL

    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    !$ACC WAIT(1)
  END SUBROUTINE perform_post_op_r2D


  !> Performs small arithmetic operations ("post-ops") on 2D REAL field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_s2D(post_op, field2D, opt_inverse, lacc)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    REAL(sp),              INTENT(INOUT) :: field2D(:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    LOGICAL , OPTIONAL   , INTENT(IN)    :: lacc

    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_s2D"
    REAL(sp) :: scalfac, lowlim   ! scale factor, lower limit
    INTEGER  :: idim(2), l1,l2
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    idim = SHAPE(field2D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          scalfac = 1._sp/post_op%arg1%sval
        ELSE
          scalfac = post_op%arg1%sval
        ENDIF
      ELSE
        scalfac = post_op%arg1%sval
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2), SCHEDULE(runtime)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) COPYIN(idim) IF(acc_is_present(field2D) .AND. lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO l2=1,idim(2)
        DO l1=1,idim(1)
          field2D(l1,l2) = field2D(l1,l2) * scalfac
        END DO
      END DO
      !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL
      !
    CASE (POST_OP_LIN2DBZ)

      lowlim  = post_op%arg1%sval
      scalfac = 10.0_sp/LOG(10._sp)

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2), SCHEDULE(runtime)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) COPYIN(idim) IF(acc_is_present(field2D) .AND. lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO l2=1,idim(2)
        DO l1=1,idim(1)
          field2D(l1,l2) = scalfac * LOG( MAX( field2D(l1,l2), lowlim) )
        END DO
      END DO
      !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL 

    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    !$ACC WAIT(1)
  END SUBROUTINE perform_post_op_s2D


  !> Performs small arithmetic operations ("post-ops") on 2D INTEGER field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_i2D(post_op, field2D, opt_inverse, lacc)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    INTEGER              , INTENT(INOUT) :: field2D(:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    LOGICAL , OPTIONAL   , INTENT(IN)    :: lacc
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_i2D"
    INTEGER  :: scalfac           ! scale factor
    INTEGER  :: idim(2), l1,l2
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    idim = SHAPE(field2D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          CALL finish(routine, "Inverse POST_OP_SCALE not applicable to INTEGER field!")
        ELSE
          scalfac = post_op%arg1%ival
        ENDIF
      ELSE
        scalfac = post_op%arg1%ival
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2), SCHEDULE(runtime)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) COPYIN(idim) IF(acc_is_present(field2D) .AND. lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO l2=1,idim(2)
        DO l1=1,idim(1)
          field2D(l1,l2) = field2D(l1,l2) * scalfac
        END DO
      END DO
      !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL

#ifndef __NO_ICON_ATMO__
    CASE(POST_OP_LUC)
#ifdef _OPENACC
      CALL finish(routine,'convert_luc_ICON2GRIB not supported on GPU')
#endif
!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2), SCHEDULE(runtime)
      DO l2=1,idim(2)
        DO l1=1,idim(1)
          field2D(l1,l2) = convert_luc_ICON2GRIB(lc_datbase  = post_op%arg1%ival, &
            &                                    iluc_in     = field2D(l1,l2),    &
            &                                    opt_linverse= opt_inverse)
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
#endif
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

  END SUBROUTINE perform_post_op_i2D


  !> Performs small arithmetic operations ("post-ops") on 3D REAL field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_r3D(post_op, field3D, opt_inverse, lacc)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    REAL(dp),              INTENT(INOUT) :: field3D(:,:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    LOGICAL , OPTIONAL   , INTENT(IN)    :: lacc
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_r3D"
    REAL(dp) :: scalfac, lowlim   ! scale factor, lower limit
    INTEGER  :: idim(3), l1,l2,l3
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    idim = SHAPE(field3D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          scalfac = 1._dp/post_op%arg1%rval
        ELSE
          scalfac = post_op%arg1%rval
        ENDIF
      ELSE
        scalfac = post_op%arg1%rval
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2,l3), SCHEDULE(runtime)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) COPYIN(idim) IF(acc_is_present(field3D) .AND. lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO l3=1,idim(3)
        DO l2=1,idim(2)
          DO l1=1,idim(1)
            field3D(l1,l2,l3) = field3D(l1,l2,l3) * scalfac
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL

    CASE (POST_OP_LIN2DBZ)

      lowlim  = post_op%arg1%rval
      scalfac = 10.0_dp/LOG(10._dp)

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2,l3), SCHEDULE(runtime)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) COPYIN(idim) IF(acc_is_present(field3D) .AND. lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO l3=1,idim(3)
        DO l2=1,idim(2)
          DO l1=1,idim(1)
            field3D(l1,l2,l3) = scalfac * LOG( MAX( field3D(l1,l2,l3), lowlim) )
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL

    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    !$ACC WAIT(1)

  END SUBROUTINE perform_post_op_r3D


  !> Performs small arithmetic operations ("post-ops") on 3D REAL field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_s3D(post_op, field3D, opt_inverse, lacc)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    REAL(sp),              INTENT(INOUT) :: field3D(:,:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    LOGICAL , OPTIONAL   , INTENT(IN)    :: lacc
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_s3D"
    REAL(sp) :: scalfac, lowlim   ! scale factor, lower limit
    INTEGER  :: idim(3), l1,l2,l3
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    idim = SHAPE(field3D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          scalfac = 1._sp/post_op%arg1%sval
        ELSE
          scalfac = post_op%arg1%sval
        ENDIF
      ELSE
        scalfac = post_op%arg1%sval
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2,l3), SCHEDULE(runtime)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) COPYIN(idim) IF(acc_is_present(field3D) .AND. lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO l3=1,idim(3)
        DO l2=1,idim(2)
          DO l1=1,idim(1)
            field3D(l1,l2,l3) = field3D(l1,l2,l3) * scalfac
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)
!$OMP END DO
!$OMP END PARALLEL
      !
    CASE (POST_OP_LIN2DBZ)

      lowlim  = post_op%arg1%sval
      scalfac = 10.0_sp/LOG(10._sp)

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2,l3), SCHEDULE(runtime)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) COPYIN(idim) IF(acc_is_present(field3D) .AND. lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO l3=1,idim(3)
        DO l2=1,idim(2)
          DO l1=1,idim(1)
            field3D(l1,l2,l3) = scalfac * LOG( MAX( field3D(l1,l2,l3), lowlim) )
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)
!$OMP END DO
!$OMP END PARALLEL

    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

  END SUBROUTINE perform_post_op_s3D


  !> Performs small arithmetic operations ("post-ops") on 3D INTEGER field 
  !  as post-processing tasks.
  SUBROUTINE perform_post_op_i3D(post_op, field3D, opt_inverse, lacc)
    TYPE (t_post_op_meta), INTENT(IN)    :: post_op
    INTEGER ,              INTENT(INOUT) :: field3D(:,:,:)
    LOGICAL , OPTIONAL   , INTENT(IN)    :: opt_inverse   ! .TRUE.: inverse operation
    LOGICAL , OPTIONAL   , INTENT(IN)    :: lacc
    !
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM(modname)//":perform_post_op_i3D"
    INTEGER  :: scalfac           ! scale factor
    INTEGER  :: idim(3), l1,l2,l3
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    idim = SHAPE(field3D)

    SELECT CASE(post_op%ipost_op_type)
    CASE(POST_OP_NONE)
      ! do nothing
    CASE(POST_OP_SCALE)
      IF (PRESENT(opt_inverse)) THEN
        IF (opt_inverse) THEN
          CALL finish(routine, "Inverse POST_OP_SCALE not applicable to INTEGER field!")
        ELSE
          scalfac = post_op%arg1%ival
        ENDIF
      ELSE
        scalfac = post_op%arg1%ival
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2,l3), SCHEDULE(runtime)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) COPYIN(idim) IF(acc_is_present(field3D) .AND. lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO l3=1,idim(3)
        DO l2=1,idim(2)
          DO l1=1,idim(1)
            field3D(l1,l2,l3) = field3D(l1,l2,l3) * scalfac
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL
      !
#ifndef __NO_ICON_ATMO__
    CASE(POST_OP_LUC)
#ifdef _OPENACC
      CALL finish(routine,'convert_luc_ICON2GRIB not supported on GPU')
#endif
!$OMP PARALLEL
!$OMP DO PRIVATE(l1,l2,l3), SCHEDULE(runtime)
      DO l3=1,idim(3)
        DO l2=1,idim(2)
          DO l1=1,idim(1)
            field3D(l1,l2,l3) = convert_luc_ICON2GRIB(lc_datbase  = post_op%arg1%ival, &
              &                                       iluc_in     = field3D(l1,l2,l3), &
              &                                       opt_linverse= opt_inverse)
          END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
#endif
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    !$ACC WAIT(1)

  END SUBROUTINE perform_post_op_i3D


END MODULE mo_post_op
