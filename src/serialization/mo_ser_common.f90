! Common serialization routines using Serialbox2
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

MODULE mo_ser_common

#ifdef SERIALIZE
  USE m_serialize,  ONLY: fs_write_field, &
                          fs_read_field
  USE utils_ppser,  ONLY: ppser_savepoint, &
                          ppser_serializer, &
                          ppser_serializer_ref, &
                          ppser_zrperturb
  USE mo_exception, ONLY: finish
  USE mo_kind,      ONLY: wp, sp, dp
  USE mo_mpi,       ONLY: get_my_mpi_work_id
  USE mo_ser_nml,   ONLY: ser_nreport, ser_nfail

  IMPLICIT NONE

  PUBLIC :: init, ser_component, t_ser_options, open_compare_file, close_compare_file

  PRIVATE

  LOGICAL :: linitialize = .TRUE.

  INTEGER :: unit_sum = 0 ! file unit for summary txt
  INTEGER :: unit_long = 0 ! file unit for detailed report

  REAL(wp) :: eps_r = EPSILON(1._wp)
  REAL(sp) :: eps_s = EPSILON(1._sp)
  !$ACC DECLARE COPYIN(eps_r, eps_s)

  INTERFACE compare
    MODULE PROCEDURE compare_r_0d
    MODULE PROCEDURE compare_r_1d
    MODULE PROCEDURE compare_r_2d
    MODULE PROCEDURE compare_r_3d
    MODULE PROCEDURE compare_r_4d
    MODULE PROCEDURE compare_s_0d
    MODULE PROCEDURE compare_s_1d
    MODULE PROCEDURE compare_s_2d
    MODULE PROCEDURE compare_s_3d
    MODULE PROCEDURE compare_s_4d
    MODULE PROCEDURE compare_i_0d
    MODULE PROCEDURE compare_i_1d
    MODULE PROCEDURE compare_i_2d
    MODULE PROCEDURE compare_i_3d
    MODULE PROCEDURE compare_i_4d
    MODULE PROCEDURE compare_l_0d
    MODULE PROCEDURE compare_l_1d
    MODULE PROCEDURE compare_l_2d
    MODULE PROCEDURE compare_l_3d
    MODULE PROCEDURE compare_l_4d
  END INTERFACE compare

  INTERFACE is_close
    MODULE PROCEDURE is_close_r
    MODULE PROCEDURE is_close_s
    MODULE PROCEDURE is_close_i
  END INTERFACE is_close

  INTERFACE ser_component
    MODULE PROCEDURE ser_scalar_d
    MODULE PROCEDURE ser_scalar_s
    MODULE PROCEDURE ser_scalar_i
    MODULE PROCEDURE ser_scalar_l
    MODULE PROCEDURE ser_array_d_1d
    MODULE PROCEDURE ser_array_d_2d
    MODULE PROCEDURE ser_array_d_3d
    MODULE PROCEDURE ser_array_d_4d
    MODULE PROCEDURE ser_array_s_1d
    MODULE PROCEDURE ser_array_s_2d
    MODULE PROCEDURE ser_array_s_3d
    MODULE PROCEDURE ser_array_s_4d
    MODULE PROCEDURE ser_array_i_1d
    MODULE PROCEDURE ser_array_i_2d
    MODULE PROCEDURE ser_array_i_3d
    MODULE PROCEDURE ser_array_i_4d
    MODULE PROCEDURE ser_array_l_1d
    MODULE PROCEDURE ser_array_l_2d
    MODULE PROCEDURE ser_array_l_3d
    MODULE PROCEDURE ser_array_l_4d
  END INTERFACE ser_component

  TYPE :: t_ser_options
    INTEGER :: abs_threshold = 12 ! absolute threshold
    INTEGER :: rel_threshold = 12 ! relative threshold
    LOGICAL :: lupdate_cpu = .FALSE. ! Update data on CPU from GPU if .TRUE.,
                                     ! used in write and compare modes.
    LOGICAL :: lopenacc = .FALSE. ! Update data on GPU from CPU if .TRUE., 
                                  ! used in read and perturb modes. Usually,
                                  ! lopenacc has to be set .TRUE. if 
                                  ! avariable is present on GPU.
    INTEGER :: ser_mode = -1 ! write(0), read(1), read perturbed(2), or compare(3)
    INTEGER :: domain = -1 ! domain index to identify field
  END TYPE t_ser_options

  CONTAINS

  SUBROUTINE init(suffix)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: suffix
    REAL(KIND=wp) :: rprecision
    rprecision = 10.0**(-PRECISION(1.0))

    !$ser verbatim    IF (linitialize) THEN

#if defined( SERIALIZE_CREATE_REFERENCE )
    !$ser init directory='./ser_data' &
    !$ser&     prefix='reference_'//TRIM(suffix) &
    !$ser&     mpi_rank=get_my_mpi_work_id() &
    !$ser&     rprecision=rprecision &
    !$ser&     rperturb=1.0e-14_wp
#else
    !$ser init directory='./ser_data' &
    !$ser&     prefix='current_'//TRIM(suffix) &
    !$ser&     prefix_ref='reference_'//TRIM(suffix) &
    !$ser&     mpi_rank=get_my_mpi_work_id() &
    !$ser&     rprecision=rprecision &
    !$ser&     rperturb=1.0e-14_wp
#endif

    !$ser verbatim     linitialize = .FALSE.
    !$ser verbatim     END IF

  END SUBROUTINE init

  SUBROUTINE open_compare_file(compare_file_name)
    CHARACTER(len=*) :: compare_file_name
    OPEN( newunit=unit_sum, file=TRIM(compare_file_name)//"_sum.txt", action="WRITE", status="REPLACE")
    OPEN( newunit=unit_long, file=TRIM(compare_file_name)//".txt", action="WRITE", status="REPLACE")
    WRITE(unit_sum, "(A, T40,A7, T50,A7, T60,A7, T70,A7, T80,A7)") "field", "rel", "abs", "%", "nfail", "ntot"
  END SUBROUTINE open_compare_file

  SUBROUTINE close_compare_file
    CLOSE( unit=unit_sum )
    CLOSE( unit=unit_long )
  END SUBROUTINE close_compare_file

  SUBROUTINE is_close_r(ref, cur, abs_threshold, rel_threshold, rel_diff, abs_diff, out)
    REAL(wp), INTENT(IN) :: ref, cur
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold
    REAL(wp), INTENT(OUT) :: rel_diff, abs_diff
    LOGICAL, INTENT(OUT) :: out

    REAL(wp) :: maxval, at, rt
    !$ACC ROUTINE SEQ

    abs_diff = ABS(cur - ref)
    maxval = MAX(ABS(cur), ABS(ref))

    ! compute relative difference for report
    IF(maxval <  eps_r) THEN
      rel_diff = 0
    ELSE
      rel_diff = abs_diff / maxval
    END IF

    ! thresholds given as negative exponents
    rt = 10._wp**(-rel_threshold)
    at = 10._wp**(-abs_threshold)

    ! threshold computed as in python's math.isclose
    out = abs_diff <= MAX(rt*maxval, at)

  END SUBROUTINE is_close_r

  SUBROUTINE is_close_s(ref, cur, abs_threshold, rel_threshold, rel_diff, abs_diff, out)
    REAL(sp), INTENT(IN) :: ref, cur
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold
    REAL(sp), INTENT(OUT) :: rel_diff, abs_diff
    LOGICAL, INTENT(OUT) :: out

    REAL(sp) :: maxval, at, rt
    !$ACC ROUTINE SEQ

    abs_diff = ABS(cur - ref)
    maxval = MAX(ABS(cur), ABS(ref))

    ! compute relative difference for report
    IF(maxval < eps_s) THEN
      rel_diff = 0
    ELSE
      rel_diff = abs_diff / maxval
    END IF

    ! thresholds given as negative exponents
    rt = 10._sp**(-rel_threshold)
    at = 10._sp**(-abs_threshold)

    ! threshold computed as in python's math.isclose
    out = abs_diff <= MAX(rt*maxval, at)

  END SUBROUTINE is_close_s

  SUBROUTINE is_close_i(ref, cur, abs_threshold, rel_threshold, rel_diff, abs_diff, out)
    INTEGER, INTENT(IN) :: ref, cur
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold
    INTEGER, INTENT(OUT) :: abs_diff
    REAL(sp), INTENT(OUT) :: rel_diff
    LOGICAL, INTENT(OUT) :: out

    INTEGER  :: maxval
    REAL(sp) :: at, rt
    !$ACC ROUTINE SEQ

    abs_diff = ABS(cur - ref)
    maxval = MAX(ABS(cur), ABS(ref))

    ! compute relative difference for report
    IF(maxval ==  0) THEN
      rel_diff = 0
    ELSE
      rel_diff = REAL(abs_diff, sp) / REAL(maxval, sp)
    END IF

    ! thresholds given as negative exponents
    rt = 10._sp**(-rel_threshold)
    at = 10._sp**(-abs_threshold)

    ! threshold computed as in python's math.isclose
    out = abs_diff <= MAX(rt*maxval, at)

  END SUBROUTINE is_close_i

  SUBROUTINE report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, n_tot)
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(wp), DIMENSION(:), INTENT(IN) :: report_rel_diff, report_abs_diff, report_cur, report_ref
    CHARACTER(len=*), DIMENSION(:), INTENT(IN) :: report_idx
    INTEGER, INTENT(IN) :: n_fail, n_tot
    REAL(wp) :: q
    INTEGER :: z

    IF(name(1:8) == '__TEST__') THEN
      IF(name(9:15) == 'error__') THEN
        ! internal test, fail expected (n_fail==n_tot)
        IF(n_fail /= n_tot) THEN
          WRITE(unit_sum, "(A, T40,A,I5,I5)") TRIM(name), " Internal test failed."
        ENDIF
      ELSEIF(name(9:17) == 'noerror__') THEN
        ! internal test, no error expected (n_fail==0)
        IF(n_fail /= 0) THEN
          WRITE(unit_sum, "(A, T40,A)") TRIM(name), " Internal test failed. (E.g. ACC update failed.)"
        ENDIF
      ELSE
        WRITE(unit_sum, "(A, T40,A)") TRIM(name), " Unexpected internal test."
      ENDIF
      RETURN
    ENDIF

    q = REAL(n_fail, wp) / REAL(n_tot, wp) * 100

    IF(q > ser_nfail) THEN
      WRITE(unit_sum, "(A, T40,E7.1E2, T50,E7.1E2, T60,F7.3, T70,I7, T80,I7)") TRIM(name), report_rel_diff(1), report_abs_diff(1), q, n_fail, n_tot
      WRITE(unit_long, "(A, A, I7, A, I7, A, F7.3, A)") TRIM(name), ": ", n_fail, " out of ", n_tot, " elements are off (", q, " %)"
      WRITE(unit_long, "(T10,A, T30,A, T50,A, T70,A, T90,A)") "rel diff", "abs diff", "current", "reference", "index"
      DO z=1, MIN(n_fail, SIZE(report_rel_diff, 1))
        WRITE(unit_long, "(T10,E14.8E2, T30,E14.8E2, T50,E14.8E2, T70,E14.8E2, T90,A)") report_rel_diff(z), report_abs_diff(z), report_cur(z), report_ref(z), TRIM(report_idx(z))
      END DO
    ENDIF
  END SUBROUTINE report

  SUBROUTINE compare_r_0d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN), TARGET :: cur
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(wp), TARGET :: ref, cur_cpy, rel_diff, abs_diff
    REAL(wp) :: report_rel_diff(1), report_abs_diff(1), report_cur(1), report_ref(1)
    CHARACTER(len=60) :: report_idx(1)
    LOGICAL :: out
    INTEGER :: n_fail

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    call is_close(ref, cur, abs_threshold, rel_threshold, rel_diff, abs_diff, out)
    IF (.NOT. out) THEN
      n_fail = 1
    ENDIF
    !$ACC END PARALLEL

    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      report_idx(1) = "(REAL(wp) scalar)"
      report_abs_diff(1) = abs_diff
      report_rel_diff(1) = rel_diff
      report_cur(1) = cur_cpy
      report_ref(1) = ref
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, 1)

  END SUBROUTINE compare_r_0d

  SUBROUTINE compare_r_1d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN) :: cur(:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(wp), DIMENSION(size(cur, 1)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(1)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(1) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      call is_close(ref(i), cur(i), abs_threshold, rel_threshold, rel_diff(i), abs_diff(i), out)
      IF (.NOT. out) THEN
        n_fail = n_fail + 1
      ENDIF
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        WRITE(report_idx(z), "(A,I6,A)") "(", i, ")"
        report_abs_diff(z) = abs_diff(i)
        report_rel_diff(z) = rel_diff(i)
        report_cur(z) = cur_cpy(i)
        report_ref(z) = ref(i)
        mask(i) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_r_1d

  SUBROUTINE compare_r_2d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN) :: cur(:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(wp), DIMENSION(size(cur, 1), size(cur, 2)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(2)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        call is_close(ref(i,j), cur(i,j), abs_threshold, rel_threshold, rel_diff(i,j), abs_diff(i,j), out)
        IF (.NOT. out) THEN
          n_fail = n_fail + 1
        ENDIF
      END DO
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        WRITE(report_idx(z), "(A,I6,A,I6,A)") "(", i, ",", j, ")"
        report_abs_diff(z) = abs_diff(i,j)
        report_rel_diff(z) = rel_diff(i,j)
        report_cur(z) = cur_cpy(i,j)
        report_ref(z) = ref(i,j)
        mask(i,j) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_r_2d

  SUBROUTINE compare_r_3d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN) :: cur(:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(wp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(3)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(3) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          call is_close(ref(i,j,k), cur(i,j,k), abs_threshold, rel_threshold, rel_diff(i,j,k), abs_diff(i,j,k), out)
          IF (.NOT. out) THEN
            n_fail = n_fail + 1
          ENDIF
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ")"
        report_abs_diff(z) = abs_diff(i,j,k)
        report_rel_diff(z) = rel_diff(i,j,k)
        report_cur(z) = cur_cpy(i,j,k)
        report_ref(z) = ref(i,j,k)
        mask(i,j,k) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_r_3d

  SUBROUTINE compare_r_4d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN) :: cur(:,:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(wp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(4)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, l, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(4) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          DO l=1,size(cur, 4)
            call is_close(ref(i,j,k,l), cur(i,j,k,l), abs_threshold, rel_threshold, rel_diff(i,j,k,l), abs_diff(i,j,k,l), out)
            IF (.NOT. out) THEN
              n_fail = n_fail + 1
            ENDIF
          END DO
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        l = idx(4)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ",", l, ")"
        report_abs_diff(z) = abs_diff(i,j,k,l)
        report_rel_diff(z) = rel_diff(i,j,k,l)
        report_cur(z) = cur_cpy(i,j,k,l)
        report_ref(z) = ref(i,j,k,l)
        mask(i,j,k,l) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_r_4d

  SUBROUTINE compare_s_0d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name

    REAL(sp), INTENT(IN), TARGET :: cur
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(sp), TARGET :: ref, cur_cpy, rel_diff, abs_diff
    REAL(wp) :: report_rel_diff(1), report_abs_diff(1), report_cur(1), report_ref(1)
    CHARACTER(len=60) :: report_idx(1)
    LOGICAL :: out
    INTEGER :: n_fail

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    call is_close(ref, cur, abs_threshold, rel_threshold, rel_diff, abs_diff, out)
    IF (.NOT. out) THEN
      n_fail = 1
    ENDIF
    !$ACC END PARALLEL

    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      report_idx(1) = "(REAL(sp) scalar)"
      report_abs_diff(1) = REAL(abs_diff, wp)
      report_rel_diff(1) = REAL(rel_diff, wp)
      report_cur(1) = cur_cpy
      report_ref(1) = ref
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, 1)

  END SUBROUTINE compare_s_0d

  SUBROUTINE compare_s_1d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(sp), INTENT(IN) :: cur(:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(sp), DIMENSION(size(cur, 1)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(1)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(1) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      call is_close(ref(i), cur(i), abs_threshold, rel_threshold, rel_diff(i), abs_diff(i), out)
      IF (.NOT. out) THEN
        n_fail = n_fail + 1
      ENDIF
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        WRITE(report_idx(z), "(A,I6,A)") "(", i, ")"
        report_abs_diff(z) = abs_diff(i)
        report_rel_diff(z) = rel_diff(i)
        report_cur(z) = cur_cpy(i)
        report_ref(z) = ref(i)
        mask(i) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_s_1d

  SUBROUTINE compare_s_2d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(sp), INTENT(IN) :: cur(:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(2)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        call is_close(ref(i,j), cur(i,j), abs_threshold, rel_threshold, rel_diff(i,j), abs_diff(i,j), out)
        IF (.NOT. out) THEN
          n_fail = n_fail + 1
        ENDIF
      END DO
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        WRITE(report_idx(z), "(A,I6,A,I6,A)") "(", i, ",", j, ")"
        report_abs_diff(z) = abs_diff(i,j)
        report_rel_diff(z) = rel_diff(i,j)
        report_cur(z) = cur_cpy(i,j)
        report_ref(z) = ref(i,j)
        mask(i,j) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_s_2d

  SUBROUTINE compare_s_3d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(sp), INTENT(IN) :: cur(:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(3)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(3) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          call is_close(ref(i,j,k), cur(i,j,k), abs_threshold, rel_threshold, rel_diff(i,j,k), abs_diff(i,j,k), out)
          IF (.NOT. out) THEN
            n_fail = n_fail + 1
          ENDIF
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ")"
        report_abs_diff(z) = abs_diff(i,j,k)
        report_rel_diff(z) = rel_diff(i,j,k)
        report_cur(z) = cur_cpy(i,j,k)
        report_ref(z) = ref(i,j,k)
        mask(i,j,k) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_s_3d

  SUBROUTINE compare_s_4d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(sp), INTENT(IN) :: cur(:,:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4)) :: ref, cur_cpy, rel_diff, abs_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4))
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(4)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, l, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(4) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          DO l=1,size(cur, 4)
            call is_close(ref(i,j,k,l), cur(i,j,k,l), abs_threshold, rel_threshold, rel_diff(i,j,k,l), abs_diff(i,j,k,l), out)
            IF (.NOT. out) THEN
              n_fail = n_fail + 1
            ENDIF
          END DO
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        l = idx(4)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ",", l, ")"
        report_abs_diff(z) = abs_diff(i,j,k,l)
        report_rel_diff(z) = rel_diff(i,j,k,l)
        report_cur(z) = cur_cpy(i,j,k,l)
        report_ref(z) = ref(i,j,k,l)
        mask(i,j,k,l) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))

  END SUBROUTINE compare_s_4d

  SUBROUTINE compare_i_0d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN), TARGET :: cur
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    INTEGER, TARGET :: ref, cur_cpy, abs_diff
    REAL(sp), TARGET :: rel_diff
    REAL(wp) :: report_rel_diff(1)
    INTEGER :: report_abs_diff(1), report_cur(1), report_ref(1)
    CHARACTER(len=60) :: report_idx(1)
    LOGICAL :: out
    INTEGER :: n_fail

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    call is_close(ref, cur, abs_threshold, rel_threshold, rel_diff, abs_diff, out)
    IF (.NOT. out) THEN
      n_fail = 1
    ENDIF
    !$ACC END PARALLEL

    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      report_idx(1) = "(INTEGER scalar)"
      report_abs_diff(1) = abs_diff
      report_rel_diff(1) = rel_diff
      report_cur(1) = cur_cpy
      report_ref(1) = ref
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, REAL(report_abs_diff, wp), REAL(report_cur, wp), REAL(report_ref, wp), report_idx, n_fail, 1)

  END SUBROUTINE compare_i_0d

  SUBROUTINE compare_i_1d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: cur(:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    INTEGER, DIMENSION(size(cur, 1)) :: ref, cur_cpy, abs_diff
    REAL(sp), DIMENSION(size(cur, 1)) :: rel_diff
    LOGICAL :: mask(size(cur, 1))
    INTEGER :: report_cur(ser_nreport), report_ref(ser_nreport), report_abs_diff(ser_nreport)
    REAL(wp) :: report_rel_diff(ser_nreport)
    INTEGER :: idx(1)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(1) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      call is_close(ref(i), cur(i), abs_threshold, rel_threshold, rel_diff(i), abs_diff(i), out)
      IF (.NOT. out) THEN
        n_fail = n_fail + 1
      ENDIF
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        WRITE(report_idx(z), "(A,I6,A)") "(", i, ")"
        report_abs_diff(z) = abs_diff(i)
        report_rel_diff(z) = rel_diff(i)
        report_cur(z) = cur_cpy(i)
        report_ref(z) = ref(i)
        mask(i) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, REAL(report_abs_diff, wp), REAL(report_cur, wp), REAL(report_ref, wp), report_idx, n_fail, size(cur))

  END SUBROUTINE compare_i_1d

  SUBROUTINE compare_i_2d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: cur(:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    INTEGER, DIMENSION(size(cur, 1), size(cur, 2)) :: ref, cur_cpy, abs_diff
    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2)) :: rel_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2))
    INTEGER :: report_cur(ser_nreport), report_ref(ser_nreport), report_abs_diff(ser_nreport)
    REAL(wp) :: report_rel_diff(ser_nreport)
    INTEGER :: idx(2)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        call is_close(ref(i,j), cur(i,j), abs_threshold, rel_threshold, rel_diff(i,j), abs_diff(i,j), out)
        IF (.NOT. out) THEN
          n_fail = n_fail + 1
        ENDIF
      END DO
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        WRITE(report_idx(z), "(A,I6,A,I6,A)") "(", i, ",", j, ")"
        report_abs_diff(z) = abs_diff(i,j)
        report_rel_diff(z) = rel_diff(i,j)
        report_cur(z) = cur_cpy(i,j)
        report_ref(z) = ref(i,j)
        mask(i,j) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, REAL(report_abs_diff, wp), REAL(report_cur, wp), REAL(report_ref, wp), report_idx, n_fail, size(cur))

  END SUBROUTINE compare_i_2d

  SUBROUTINE compare_i_3d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: cur(:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    INTEGER, DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3)) :: ref, cur_cpy, abs_diff
    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3)) :: rel_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3))
    INTEGER :: report_cur(ser_nreport), report_ref(ser_nreport), report_abs_diff(ser_nreport)
    REAL(wp) :: report_rel_diff(ser_nreport)
    INTEGER :: idx(3)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(3) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          call is_close(ref(i,j,k), cur(i,j,k), abs_threshold, rel_threshold, rel_diff(i,j,k), abs_diff(i,j,k), out)
          IF (.NOT. out) THEN
            n_fail = n_fail + 1
          ENDIF
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ")"
        report_abs_diff(z) = abs_diff(i,j,k)
        report_rel_diff(z) = rel_diff(i,j,k)
        report_cur(z) = cur_cpy(i,j,k)
        report_ref(z) = ref(i,j,k)
        mask(i,j,k) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, REAL(report_abs_diff, wp), REAL(report_cur, wp), REAL(report_ref, wp), report_idx, n_fail, size(cur))

  END SUBROUTINE compare_i_3d

  SUBROUTINE compare_i_4d(name, cur, lopenacc, abs_threshold, rel_threshold)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: cur(:,:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc
    INTEGER, INTENT(IN) :: abs_threshold, rel_threshold

    INTEGER, DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4)) :: ref, cur_cpy, abs_diff
    REAL(sp), DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4)) :: rel_diff
    LOGICAL :: mask(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4))
    INTEGER :: report_cur(ser_nreport), report_ref(ser_nreport), report_abs_diff(ser_nreport)
    REAL(wp) :: report_rel_diff(ser_nreport)
    INTEGER :: idx(4)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, l, z

    n_fail = 0
    !$ACC DATA CREATE(ref, rel_diff, abs_diff) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(4) REDUCTION(+: n_fail)
    DO i=1,size(cur, 1)
      DO j=1,size(cur, 2)
        DO k=1,size(cur, 3)
          DO l=1,size(cur, 4)
            call is_close(ref(i,j,k,l), cur(i,j,k,l), abs_threshold, rel_threshold, rel_diff(i,j,k,l), abs_diff(i,j,k,l), out)
            IF (.NOT. out) THEN
              n_fail = n_fail + 1
            ENDIF
          END DO
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    ! compute additional info on CPU
    IF (n_fail > 0) THEN
      !$ACC UPDATE HOST(rel_diff, abs_diff) ASYNC(1) IF(lopenacc)
      mask = .TRUE.
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      !$ACC WAIT(1) IF(lopenacc)
      DO z=1,ser_nreport
        idx(:) = MAXLOC(rel_diff, mask)
        i = idx(1)
        j = idx(2)
        k = idx(3)
        l = idx(4)
        WRITE(report_idx(z), "(A,I6,A,I6,A,I6,A,I6,A)") "(", i, ",", j, ",", k, ",", l, ")"
        report_abs_diff(z) = abs_diff(i,j,k,l)
        report_rel_diff(z) = rel_diff(i,j,k,l)
        report_cur(z) = cur_cpy(i,j,k,l)
        report_ref(z) = ref(i,j,k,l)
        mask(i,j,k,l) = .FALSE.
      END DO
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, REAL(report_abs_diff, wp), REAL(report_cur, wp), REAL(report_ref, wp), report_idx, n_fail, size(cur))

  END SUBROUTINE compare_i_4d

  SUBROUTINE compare_l_0d(name, cur, lopenacc)
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(IN), TARGET :: cur
    LOGICAL, INTENT(IN) :: lopenacc

    LOGICAL, TARGET :: ref, cur_cpy
    REAL(wp) :: report_rel_diff(1), report_abs_diff(1), report_cur(1), report_ref(1)
    CHARACTER(len=60) :: report_idx(1)
    LOGICAL :: out
    INTEGER :: n_fail

    n_fail = 0
    !$ACC DATA CREATE(ref) PRESENT(cur) IF(lopenacc)

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC UPDATE DEVICE(ref) ASYNC(1) IF(lopenacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lopenacc)
    IF (ref .NEQV. cur) THEN
      n_fail = 1
    ENDIF
    !$ACC END PARALLEL

    IF (n_fail > 0) THEN
      !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
      cur_cpy = cur
      !$ACC END KERNELS
      report_idx(1) = "(LOGICAL scalar)"
      report_abs_diff(1) = 1._wp
      report_rel_diff(1) = 1._wp
      report_cur(1) = MERGE(1._wp, 0._wp, cur_cpy)
      report_ref(1) = MERGE(1._wp, 0._wp, ref)
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, 1)

  END SUBROUTINE compare_l_0d

  SUBROUTINE compare_l_1d(name, cur, lopenacc)
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(IN) :: cur(:)
    LOGICAL, INTENT(IN) :: lopenacc

    LOGICAL, DIMENSION(size(cur, 1)), TARGET :: ref, cur_cpy
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(1)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i

    n_fail = 0

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
    cur_cpy = cur
    !$ACC END KERNELS

    DO i=1,size(cur, 1)
      IF (ref(i) .NEQV. cur_cpy(i)) THEN
        n_fail = n_fail + 1
        IF (n_fail <= ser_nreport) THEN
          WRITE(report_idx(n_fail), "(A,I6,A)") "(", i, ")"
          report_abs_diff(n_fail) = 1._wp
          report_rel_diff(n_fail) = 1._wp
          report_cur(n_fail) = MERGE(1._wp, 0._wp, cur_cpy(i))
          report_ref(n_fail) = MERGE(1._wp, 0._wp, ref(i))
        ENDIF
      ENDIF
    END DO

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))
    !$ACC WAIT(1)
  END SUBROUTINE compare_l_1d

  SUBROUTINE compare_l_2d(name, cur, lopenacc)
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(IN) :: cur(:,:)
    LOGICAL, INTENT(IN) :: lopenacc

    LOGICAL, DIMENSION(size(cur, 1), size(cur, 2)), TARGET :: ref, cur_cpy
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(2)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j

    n_fail = 0

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
    cur_cpy = cur
    !$ACC END KERNELS

    DO j=1,size(cur, 2)
      DO i=1,size(cur, 1)
        IF (ref(i,j) .NEQV. cur_cpy(i,j)) THEN
          n_fail = n_fail + 1
          IF (n_fail <= ser_nreport) THEN
            WRITE(report_idx(n_fail), "(A,I6,',',I6,A)") "(", i, j, ")"
            report_abs_diff(n_fail) = 1._wp
            report_rel_diff(n_fail) = 1._wp
            report_cur(n_fail) = MERGE(1._wp, 0._wp, cur_cpy(i,j))
            report_ref(n_fail) = MERGE(1._wp, 0._wp, ref(i,j))
          ENDIF
        ENDIF
      END DO
    END DO

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))
    !$ACC WAIT(1)
  END SUBROUTINE compare_l_2d

  SUBROUTINE compare_l_3d(name, cur, lopenacc)
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(IN) :: cur(:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc

    LOGICAL, DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3)), TARGET :: ref, cur_cpy
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(3)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k

    n_fail = 0

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
    cur_cpy = cur
    !$ACC END KERNELS

    DO k=1,size(cur, 3)
      DO j=1,size(cur, 2)
        DO i=1,size(cur, 1)
          IF (ref(i,j,k) .NEQV. cur_cpy(i,j,k)) THEN
            n_fail = n_fail + 1
            IF (n_fail <= ser_nreport) THEN
              WRITE(report_idx(n_fail), "(A,I6,',',I6,',',I6,A)") "(", i, j, k, ")"
              report_abs_diff(n_fail) = 1._wp
              report_rel_diff(n_fail) = 1._wp
              report_cur(n_fail) = MERGE(1._wp, 0._wp, cur_cpy(i,j,k))
              report_ref(n_fail) = MERGE(1._wp, 0._wp, ref(i,j,k))
            ENDIF
          ENDIF
        END DO
      END DO
    END DO

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))
    !$ACC WAIT(1)
  END SUBROUTINE compare_l_3d

  SUBROUTINE compare_l_4d(name, cur, lopenacc)
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(IN) :: cur(:,:,:,:)
    LOGICAL, INTENT(IN) :: lopenacc

    LOGICAL, DIMENSION(size(cur, 1), size(cur, 2), size(cur, 3), size(cur, 4)), TARGET :: ref, cur_cpy
    REAL(wp) :: report_rel_diff(ser_nreport), report_abs_diff(ser_nreport), report_cur(ser_nreport), report_ref(ser_nreport)
    INTEGER :: idx(4)
    CHARACTER(len=60) :: report_idx(ser_nreport)
    LOGICAL :: out
    INTEGER :: n_fail
    INTEGER :: i, j, k, l

    n_fail = 0

    call fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name), ref)
    !$ACC KERNELS COPYOUT(cur_cpy) ASYNC(1) IF(lopenacc)
    cur_cpy = cur
    !$ACC END KERNELS

    DO l=1,size(cur, 4)
      DO k=1,size(cur, 3)
        DO j=1,size(cur, 2)
          DO i=1,size(cur, 1)
            IF (ref(i,j,k,l) .NEQV. cur_cpy(i,j,k,l)) THEN
              n_fail = n_fail + 1
              IF (n_fail <= ser_nreport) THEN
                WRITE(report_idx(n_fail), "(A,I6,',',I6,',',I6,',',I6,A)") "(", i, j, k, l, ")"
                report_abs_diff(n_fail) = 1._wp
                report_rel_diff(n_fail) = 1._wp
                report_cur(n_fail) = MERGE(1._wp, 0._wp, cur_cpy(i,j,k,l))
                report_ref(n_fail) = MERGE(1._wp, 0._wp, ref(i,j,k,l))
              ENDIF
            ENDIF
          END DO
        END DO
      END DO
    END DO

    call report(name, report_rel_diff, report_abs_diff, report_cur, report_ref, report_idx, n_fail, size(cur))
    !$ACC WAIT(1)
  END SUBROUTINE compare_l_4d

! For ser_component INTERFACE
  SUBROUTINE ser_scalar_d(o, name, p)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(dp), TARGET, INTENT(INOUT) :: p
    REAL(dp), POINTER :: pp

    pp => p ! we have to pass a pointer to ppser_*

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(p) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, name, pp)
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, name,  pp)
        !$ACC UPDATE DEVICE(p) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, name,  pp,  ppser_zrperturb)
        !$ACC UPDATE DEVICE(p) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(p) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(name,  p, o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_scalar_d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_scalar_d

  SUBROUTINE ser_scalar_s(o, name, p)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(sp), TARGET, INTENT(INOUT) :: p
    REAL(sp), POINTER :: pp

    pp => p ! we have to pass a pointer to ppser_*

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(p) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, name, pp)
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, name,  pp)
        !$ACC UPDATE DEVICE(p) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, name,  pp,  ppser_zrperturb)
        !$ACC UPDATE DEVICE(p) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(p) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(name,  p, o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_scalar_s', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_scalar_s

  SUBROUTINE ser_scalar_i(o, name, p)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, TARGET, INTENT(INOUT) :: p
    INTEGER, POINTER :: pp

    pp => p ! we have to pass a pointer to ppser_*

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(p) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, name, pp)
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, name,  pp)
        !$ACC UPDATE DEVICE(p) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, name,  pp,  ppser_zrperturb)
        !$ACC UPDATE DEVICE(p) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(p) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(name,  p, o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_scalar_i', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_scalar_i

  SUBROUTINE ser_scalar_l(o, name, p)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    LOGICAL, TARGET, INTENT(INOUT) :: p
    LOGICAL, POINTER :: pp

    pp => p ! we have to pass a pointer to ppser_*

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(p) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, name, pp)
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, name,  pp)
        !$ACC UPDATE DEVICE(p) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, name,  pp,  ppser_zrperturb)
        !$ACC UPDATE DEVICE(p) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(p) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(name,  p, o%lopenacc)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_scalar_l', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_scalar_l

  SUBROUTINE ser_array_d_1d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(dp), TARGET, INTENT(INOUT) :: ptr(:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_d_1d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_d_1d

  SUBROUTINE ser_array_d_2d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(dp), TARGET, INTENT(INOUT) :: ptr(:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_d_2d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_d_2d

  SUBROUTINE ser_array_d_3d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(dp), TARGET, INTENT(INOUT) :: ptr(:,:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:,:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_d_3d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_d_3d

  SUBROUTINE ser_array_d_4d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(dp), TARGET, INTENT(INOUT) :: ptr(:,:,:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:,:,:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_d_4d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_d_4d

  SUBROUTINE ser_array_s_1d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(sp), TARGET, INTENT(INOUT) :: ptr(:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_s_1d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_s_1d

  SUBROUTINE ser_array_s_2d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(sp), TARGET, INTENT(INOUT) :: ptr(:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_s_2d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_s_2d

  SUBROUTINE ser_array_s_3d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(sp), TARGET, INTENT(INOUT) :: ptr(:,:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:,:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_s_3d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_s_3d

  SUBROUTINE ser_array_s_4d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(sp), TARGET, INTENT(INOUT) :: ptr(:,:,:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:,:,:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_s_4d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_s_4d

  SUBROUTINE ser_array_i_1d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, TARGET, INTENT(INOUT) :: ptr(:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_i_1d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_i_1d

  SUBROUTINE ser_array_i_2d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, TARGET, INTENT(INOUT) :: ptr(:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_i_2d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_i_2d

  SUBROUTINE ser_array_i_3d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, TARGET, INTENT(INOUT) :: ptr(:,:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:,:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_i_3d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_i_3d

  SUBROUTINE ser_array_i_4d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, TARGET, INTENT(INOUT) :: ptr(:,:,:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:,:,:), o%lopenacc, o%abs_threshold, o%rel_threshold)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_i_4d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_i_4d

  SUBROUTINE ser_array_l_1d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    LOGICAL, TARGET, INTENT(INOUT) :: ptr(:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:), o%lopenacc)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_l_1d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_l_1d

  SUBROUTINE ser_array_l_2d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    LOGICAL, TARGET, INTENT(INOUT) :: ptr(:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:), o%lopenacc)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_l_2d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_l_2d

  SUBROUTINE ser_array_l_3d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    LOGICAL, TARGET, INTENT(INOUT) :: ptr(:,:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:,:), o%lopenacc)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_l_3d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_l_3d

  SUBROUTINE ser_array_l_4d(o, name, ptr)
    TYPE(t_ser_options), INTENT(IN) :: o
    CHARACTER(len=*), INTENT(IN) :: name
    LOGICAL, TARGET, INTENT(INOUT) :: ptr(:,:,:,:)
    CHARACTER(len=3) :: c="   "

    IF(o%domain >= 0) WRITE(c, '("_",i2.2)') o%domain

    SELECT CASE ( o%ser_mode )
      CASE(0) ! write
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL fs_write_field(ppser_serializer, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:))
      CASE(1) ! read
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:))
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(2) ! read perturb
        CALL fs_read_field(ppser_serializer_ref, ppser_savepoint, TRIM(name//c), ptr(:,:,:,:), ppser_zrperturb)
        !$ACC UPDATE DEVICE(ptr) ASYNC(1) IF(o%lopenacc)
      CASE(3) ! compare
        !$ACC UPDATE HOST(ptr) ASYNC(1) IF(o%lupdate_cpu)
        !$ACC WAIT(1) IF(o%lupdate_cpu)
        CALL compare(TRIM(name//c), ptr(:,:,:,:), o%lopenacc)
      CASE DEFAULT
        CALL finish('mo_ser_manually:ser_array_l_4d', 'Internal error. Unknown ser_mode.')
    END SELECT

  END SUBROUTINE ser_array_l_4d
#endif
END MODULE mo_ser_common
